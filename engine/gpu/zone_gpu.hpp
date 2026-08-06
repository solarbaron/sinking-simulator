// SPDX-License-Identifier: MIT
//
// Vulkan compute back-end for the explicit **solid-shell** zone solver -- the
// Tier-2 element kernel of `engine/sim/zone.{hpp,cpp}` on the GPU.
//
// `fem_gpu.{hpp,cpp}` is the pattern, and the two are deliberately the same shape:
// per-element force slots, a companion node kernel that gathers them over a CSR
// adjacency in a fixed order, and every substep recorded into one command buffer so
// the CPU is not in the loop. What differs is the element -- eight nodes and 24
// degrees of freedom against a tet's four and twelve, seven enhanced assumed strains
// condensed at element level, and eight integration points of plastic history that
// live on the device.
//
// **The one part of the tet's pattern that does not carry over is one invocation per
// element**, and that took a measurement to establish rather than an argument. A
// linear tet's whole state fits in registers -- the driver reports zero spill for
// `tet_forces.comp`. This element's does not: one thread per element spills 1 936
// bytes, 484 floats, and Pascal serves that out of global memory. Both mappings are
// built and `Mapping` selects between them; see below and
// `07-fem-spike-findings.md` §8.
//
// --- Why this exists in the shape it does, which was measured first --------------
//
// The roadmap names the tet back-end as "the pattern to follow" for a solid-shell
// GPU element solver, and before writing any Vulkan the question was whether the
// pattern applies at all -- a GPU element kernel is worth nothing if element
// evaluation is not where the time goes. Profiled on the ferry's own side patch
// (`tools/zone_probe`), driving a punch into 192 elements for 6 608 steps:
//
//     element evaluation   98.5%      one worker
//     CSR nodal gather      0.3%
//     integration           0.2%
//     energy accounting     1.0%
//
// So the pattern applies. But **inside** that 98.5% the profile found something
// the roadmap did not: half of it was `computeForms`, rebuilding each element's
// strain-displacement matrices from its *rest* configuration, which an explicit
// solve never moves. Hoisting that out (`solidshell::RestForms`) is a 2.0x
// end-to-end speedup on the CPU for no change in arithmetic -- the two answers are
// bit-identical -- and it is what makes this kernel tractable, because the ANS/EAS
// geometry pipeline never has to exist in GLSL. The tet has always had this: it
// uploads `restInv` and `restVolume` and reads them from a buffer.
//
// **The GPU baseline is therefore the fixed CPU, not the one that was wasting half
// its time.** Quoting a speedup against the unfixed path would have doubled the
// number for free.
//
// --- What lives on the device ----------------------------------------------------
//
// Per element: 1 505 floats of `RestForms` (B, G, the Gauss weights and the rest
// Jacobian inverse), 24 force slots, and 129 floats of state -- eight integration
// points of plastic strain, back stress, accumulated plastic strain and damage,
// plus the **seven enhanced parameters**. The enhanced parameters are per-element
// internal variables and they are condensed out inside the element every step, by
// a Newton on `int G^T sigma dV = 0` solved with a 7x7 Cholesky in the shader --
// the same algorithm the CPU runs, with a float tolerance rather than a double
// one. They are read at the start of a step as a warm start and written back at the
// end, so they are device-resident state and never cross the bus.
//
// Plasticity is **in scope**, and it is the reason the kernel is the size it is.
// It is also the only version worth building: the elastic element is 0.27 µs
// against the elastoplastic 7.3, so a GPU path for the elastic one would be
// accelerating a twenty-seventh of the bill.
//
// --- Precision -------------------------------------------------------------------
//
// Float, because the target is a GTX 1070 Ti where fp64 runs at 1/32 rate. Whether
// that is enough is a measurement and not an assumption -- `tools/zone_gpu_probe`
// runs this against the CPU double reference on the same patch and reports the
// divergence in the quantities that mean anything.
//
// **It is not enough, and the remap did not change that by so much as a digit.**
// Where it fails is the torn set: at 768 and 3 072 elements this kernel tears 41 and
// 248 elements where the double reference tears 32 and 162, while the negative
// control -- the same double solver on a mesh jittered by the size of float's own
// representation error -- tears exactly 32 and 162. That is the whole argument, and
// it is why the CPU is still the path a zone's answer comes from.
//
// **Keeping `alpha` in double was the last thing left to try, and it does not work
// either.** `EasPrecision` below compiles the enhanced block in fp64 at three depths,
// all from `solidshell_forces_wg.comp` so that nothing but the precision differs. The
// five resulting kernels tear between 40 and 49 elements at 768 and between 205 and 268
// at 3 072 -- a spread the size of the gap they were meant to close, and not monotone in
// precision. Holding the stopping rule fixed, the whole block in fp64 moves alpha by
// 1/810 of the amount float is already wrong by, because Kaa's inputs are a float
// tangent and a float stress and widening what consumes them recovers nothing they never
// had. It costs 5-10x on the kernel. `docs/07-fem-spike-findings.md` §8 has the tables.
//
// It inherits §2's reproducibility bound: an explicit scheme at the CFL limit
// amplifies float rounding in the modes at its stability boundary, so this is not
// bit-reproducible against the CPU and will not be across drivers either.
#pragma once

#include "../sim/plasticity.hpp"
#include "../sim/zone.hpp"

#include <string>
#include <vector>

namespace gpu {

// What the device produced, in the CPU's own types so `zone::Solver::adopt` can
// take it and the validated energy, tearing and panel-reporting paths run on it
// unchanged.
struct ZoneGpuState {
    std::vector<double> position;
    std::vector<double> velocity;
    std::vector<sim::solidshell::ElementPlasticState> plastic;
    double work = 0;         // J the punch put in, accumulated on the device
    double dissipation = 0;  // J plastic, accumulated on the device
    int    steps = 0;
    double time = 0;
    double penetration = 0;
    int    tornElements = 0;
};

// How the element kernel is mapped onto the device. **Both compute the same thing**
// -- the same element, the same arithmetic, the same accumulation order, the same 24
// force slots -- and they exist together so the mapping can be A/B'd on one run
// rather than argued about.
//
// `Invocation` is the tet back-end's pattern: one thread per element. It needs about
// five hundred floats of thread-private, dynamically indexed state, which is twice
// Pascal's register file per thread, so all of it spills to local memory.
// `Workgroup` is 32 threads per element with that state in shared memory and the
// eight Gauss points as threads rather than a loop. `07-fem-spike-findings.md` §8
// has what each measures.
enum class Mapping {
    Invocation,  // solidshell_forces.comp
    Workgroup,   // solidshell_forces_wg.comp
};

// **How much of the 7x7 enhanced-strain block runs in fp64.** `07-fem-spike-findings.md`
// §8's last open item was "keep alpha in double -- the only part that has been shown to
// need the digits". That claim bundles three separable things together, so these are the
// three, each a superset of the one before, and all of them are the **workgroup**
// mapping only: they are variants of `solidshell_forces_wg.comp` compiled from one
// source, so any difference between them is the precision and nothing else.
//
// Pascal runs fp64 at 1/32 rate, so what each costs is a measurement rather than an
// argument; §8 has the table.
enum class EasPrecision {
    Float,      // as shipped -- the whole element, the enhanced block included, in float
    // **A control, not a step on the ladder.** Float arithmetic with the CPU's own
    // stopping rule (|delta.r| <= 1e-16 sigma_y V, 40 iterations) instead of the
    // float kernel's 1e-9 and 12. `Newton` below changes the tolerance and the
    // arithmetic together; without this there is no way to say which of the two moved
    // the answer.
    FloatTight,
    Solve,      // + the equilibration, the 7x7 Cholesky and its substitutions in double
    Condense,   // + Kaa and the residual accumulated in double, over a double G
    Newton,     // + alpha, its correction, its persistent storage and the CPU's rule
};

// **What the driver's compiler actually did with each mapping**, via
// `VK_KHR_pipeline_executable_properties`: registers per thread, and how many bytes
// of local memory the shader spills to. `07-fem-spike-findings.md` §8 attributed the
// one-invocation kernel's throughput curve to spilling, which was a diagnosis from
// the *shape* of the curve rather than from the compiler; this reads it out.
//
// Builds its own throwaway device and pipelines so it cannot perturb the path being
// timed -- capturing statistics needs a pipeline-creation flag the spec allows to
// change codegen, so it must not be set on a pipeline anyone measures. Returns a
// human-readable report, or a line saying why it could not (no device, or a driver
// without the extension). **Never fails; it is an instrument, not a dependency.**
std::string describeElementPipelines(const std::string& shaderDirectory);

class ZoneGpuSolver {
public:
    ~ZoneGpuSolver();
    ZoneGpuSolver() = default;
    ZoneGpuSolver(const ZoneGpuSolver&) = delete;
    ZoneGpuSolver& operator=(const ZoneGpuSolver&) = delete;

    // Uploads the patch and builds both pipelines. Returns false with a reason on
    // any failure, **including "no Vulkan device" -- this must degrade, not
    // abort**, which is why every caller checks the return rather than the device.
    //
    // `cpu` supplies the lumped mass, the CSR adjacency, the driven-node list, the
    // pinned degrees of freedom and the rest forms. They are taken from it rather
    // than rebuilt so that the two paths cannot silently stop being the same
    // solver; `cpu` must therefore have been constructed with
    // `SolveParams::cacheRestForms` on.
    // An `eas` other than `Float` needs the workgroup mapping and a device with
    // `shaderFloat64`; both are reported through `error` and both **skip rather than
    // fail**, like a missing device.
    bool initialise(const sim::zone::Patch& patch, const sim::zone::Solver& cpu,
                    const sim::plasticity::Material& material,
                    const sim::zone::SolveParams& params, const std::string& shaderDirectory,
                    std::string& error, Mapping mapping = Mapping::Workgroup,
                    EasPrecision eas = EasPrecision::Float);

    // Records and submits `substeps` steps as one command buffer and returns the
    // **GPU-measured** milliseconds. Batching is the whole point: an explicit step
    // on a few hundred elements is far shorter than a submission round trip, so the
    // CPU must not be in the loop. The punch's ramp varies per step and is pushed
    // between dispatches inside the same buffer.
    double run(int substeps);

    // Everything the device holds, back in double. Costs a full round trip, so it
    // is called once at the end of a run and not per step.
    ZoneGpuState readback() const;

    const std::string& deviceName() const { return deviceName_; }
    // Which element shader this solver actually loaded. Exposed for one reason:
    // `tests/test_zone_gpu.cpp` compares the two mappings against each other, and
    // that comparison is **vacuous if both of them loaded the same kernel** -- two
    // runs of one shader agree perfectly and every tolerance passes. Mutation
    // testing proved it: a mutant collapsing the ternary below to a single filename
    // survived the whole suite. The guard is now an equality on this string.
    const std::string& elementShader() const { return elementShader_; }
    bool valid() const { return device_ != nullptr; }
    // Steps recorded so far, and the simulated time and penetration they represent.
    int steps() const { return steps_; }

private:
    struct Impl;
    Impl* impl_ = nullptr;
    void* device_ = nullptr;  // mirrors impl_->device, for valid()
    std::string deviceName_;
    std::string elementShader_;
    int steps_ = 0;
};

}  // namespace gpu
