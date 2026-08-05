// SPDX-License-Identifier: MIT
//
// Vulkan compute back-end for the explicit **solid-shell** zone solver -- the
// Tier-2 element kernel of `engine/sim/zone.{hpp,cpp}` on the GPU.
//
// `fem_gpu.{hpp,cpp}` is the pattern, and the two are deliberately the same shape:
// one invocation per element writing to a per-element force slot, a companion node
// kernel that gathers them over a CSR adjacency in a fixed order, and every substep
// recorded into one command buffer so the CPU is not in the loop. What differs is
// the element -- eight nodes and 24 degrees of freedom against a tet's four and
// twelve, seven enhanced assumed strains condensed at element level, and eight
// integration points of plastic history that live on the device.
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
// divergence in the quantities that mean anything. Two places were expected to
// hurt and one of them does; the numbers are in `docs/07-fem-spike-findings.md` §8.
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
    bool initialise(const sim::zone::Patch& patch, const sim::zone::Solver& cpu,
                    const sim::plasticity::Material& material,
                    const sim::zone::SolveParams& params, const std::string& shaderDirectory,
                    std::string& error);

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
    bool valid() const { return device_ != nullptr; }
    // Steps recorded so far, and the simulated time and penetration they represent.
    int steps() const { return steps_; }

private:
    struct Impl;
    Impl* impl_ = nullptr;
    void* device_ = nullptr;  // mirrors impl_->device, for valid()
    std::string deviceName_;
    int steps_ = 0;
};

}  // namespace gpu
