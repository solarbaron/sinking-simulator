// SPDX-License-Identifier: MIT
//
// The GPU solid-shell back-end. **Skips rather than fails without a Vulkan
// device**, per CLAUDE.md, so this file's first job is to say nothing at all on a
// machine with no GPU.
//
// What it can assert cheaply is the part that is not physics: that the buffer
// layout the host writes is the layout the shader reads. That is a pure
// bookkeeping fact and it is the single most likely thing to be silently wrong --
// a stride off by one in either file gives an answer that is plausible everywhere
// and correct nowhere. The physics comparison against the CPU double reference
// costs core-minutes and lives in `tools/zone_gpu_probe`, for the same reason
// `zone_probe` is a tool.
#include "engine/gpu/zone_gpu.hpp"
#include "engine/sim/zone.hpp"
#include "harness.hpp"

#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

using testing::expectTrue;
using namespace sim;

namespace {

// The shader declares its own copies of the strides. They are read back out of the
// GLSL source and checked against the C++ ones, because two constants that must
// agree and live in two files are two constants that will not.
bool constantFromShader(const std::string& source, const std::string& name, long long& out) {
    const std::string key = "const uint " + name + " = ";
    const std::size_t at = source.find(key);
    if (at == std::string::npos) return false;
    out = std::atoll(source.c_str() + at + key.size());
    return true;
}

// Two more spellings the reader above does not know, and both had a C++ twin that
// nothing compared: `const float NAME = ` for the fixed-point work scale, and
// `layout(local_size_x = N)` for the dispatch geometry.
bool floatConstantFromShader(const std::string& source, const std::string& name, double& out) {
    const std::string key = "const float " + name + " = ";
    const std::size_t at = source.find(key);
    if (at == std::string::npos) return false;
    out = std::atof(source.c_str() + at + key.size());
    return true;
}

bool localSizeXFromShader(const std::string& source, long long& out) {
    const std::string key = "local_size_x = ";
    const std::size_t at = source.find(key);
    if (at == std::string::npos) return false;
    out = std::atoll(source.c_str() + at + key.size());
    return true;
}

// **The dispatch geometry and the work scale, which the layout check above skips.**
// A wrong `local_size_x` dispatches the wrong number of groups -- `groupsFor` divides
// by the host's copy -- and a wrong `kWorkScale` mis-scales every `atomicAdd` into the
// work accumulator. The GPU/CPU work comparison further down would catch the second,
// but it *skips without a Vulkan device* and this does not.
void testTheShadersAndTheHostAgreeOnTheirDispatch() {
    struct Pair { const char* shader; long long group; };
    const Pair pairs[] = {{"solidshell_forces.comp", gpu::kElementGroup},
                          {"solidshell_integrate.comp", gpu::kNodeGroup}};
    for (const Pair& p : pairs) {
        const std::string path = std::string(SHIPSIM_SHADER_SOURCE_DIR) + "/" + p.shader;
        std::ifstream file(path);
        expectTrue(std::string(p.shader) + " is readable", file.good());
        if (!file.good()) continue;
        const std::string source((std::istreambuf_iterator<char>(file)),
                                 std::istreambuf_iterator<char>());
        long long got = -1;
        expectTrue(std::string(p.shader) + " declares a workgroup size",
                   localSizeXFromShader(source, got));
        testing::expectEqual(std::string(p.shader) + "'s workgroup matches the host's dispatch", got,
                    p.group);
    }

    const std::string path =
        std::string(SHIPSIM_SHADER_SOURCE_DIR) + "/solidshell_integrate.comp";
    std::ifstream file(path);
    if (file.good()) {
        const std::string source((std::istreambuf_iterator<char>(file)),
                                 std::istreambuf_iterator<char>());
        double scale = -1;
        expectTrue("the integrate shader declares kWorkScale",
                   floatConstantFromShader(source, "kWorkScale", scale));
        // Exact: the host divides by its copy what the shader multiplied by its own,
        // so anything but equality is a scaled answer nobody would question.
        testing::expectNear("and it is the scale the host divides by", scale, gpu::kWorkScale, 0.0);
    }
}

// **Both** force shaders, because there are now two mappings of the same element and
// they read the same buffers. A stride that drifts in one file and not the other is
// exactly the kind of difference that would show up as "the new mapping disagrees
// with the CPU" and be blamed on the mapping.
void testTheShadersAndTheHostAgreeOnTheLayout() {
    // Derived here from the element's own constants, so a change to kGauss, kDof or
    // kEas moves this side and both shaders have to follow.
    const long long b = solidshell::kGauss * 6 * solidshell::kDof;
    const long long g = b + solidshell::kGauss * 6 * solidshell::kEas;
    const long long w = g + solidshell::kGauss;
    const long long stride = w + 9;
    const long long enhanced = solidshell::kGauss * 15;
    const long long failure = enhanced + solidshell::kEas;

    struct Check { const char* name; long long want; };
    const Check checks[] = {
        {"kFormG", b}, {"kFormW", g}, {"kFormJ", w}, {"kFormStride", stride},
        {"kPointStride", 15}, {"kStateEnhanced", enhanced}, {"kStateFailure", failure},
        {"kStateTorn", failure + 1}, {"kStateStride", failure + 2},
    };
    const char* shaders[] = {"solidshell_forces.comp", "solidshell_forces_wg.comp"};
    for (const char* shader : shaders) {
        const std::string path = std::string(SHIPSIM_SHADER_SOURCE_DIR) + "/" + shader;
        std::ifstream file(path);
        expectTrue(std::string(shader) + " is readable", file.good());
        if (!file.good()) continue;
        const std::string source((std::istreambuf_iterator<char>(file)),
                                 std::istreambuf_iterator<char>());
        for (const Check& c : checks) {
            long long got = -1;
            const bool found = constantFromShader(source, c.name, got);
            expectTrue(std::string(shader) + " declares " + c.name, found);
            expectTrue(std::string("and its ") + c.name + " matches the element's",
                       found && got == c.want);
        }
    }
    // The fp64 enhanced forms are the workgroup kernel's alone -- a second, double copy
    // of G and the Gauss weights, with B deliberately left out because the enhanced
    // block never reads it. Same argument as above: two strides in two files.
    {
        const long long g = solidshell::kGauss * 6 * solidshell::kEas;   // 336
        const Check easChecks[] = {
            {"kEasFormG", 0}, {"kEasFormW", g}, {"kEasFormStride", g + solidshell::kGauss},
        };
        const std::string path =
            std::string(SHIPSIM_SHADER_SOURCE_DIR) + "/solidshell_forces_wg.comp";
        std::ifstream file(path);
        if (file.good()) {
            const std::string source((std::istreambuf_iterator<char>(file)),
                                     std::istreambuf_iterator<char>());
            for (const Check& c : easChecks) {
                long long got = -1;
                const bool found = constantFromShader(source, c.name, got);
                expectTrue(std::string("the workgroup shader declares ") + c.name, found);
                expectTrue(std::string("and its ") + c.name + " matches the element's",
                           found && got == c.want);
            }
        }
    }
    // Guard against a vacuous pass: if the parser silently returned zero for
    // everything the loop above would compare zeros to zeros on a mis-sized element.
    expectTrue("and the strides being compared are not degenerate",
               stride > 1000 && enhanced == 120);
}

// The back-end must decline, not crash, on the two configurations it cannot serve
// -- and it must say which. Both are reachable without a device, because they are
// checked before Vulkan is touched.
void testItDeclinesWhatItCannotDo() {
    solidshell::HexMesh flat = solidshell::makePlateMesh(0.4, 0.2, 0.01, 4, 2, 1);
    zone::Patch patch;
    patch.mesh = flat;
    patch.axis = Vec3{0, 0, 1};
    patch.elementArea.assign(flat.elementCount(), 0.02);
    patch.panelOf.assign(flat.elementCount(), 0);
    patch.material = ah36Steel();
    patch.thickness = 0.01;
    patch.criticalTimestep =
        solidshell::criticalTimestep(flat, patch.material, solidshell::Formulation::SolidShell);

    {
        zone::SolveParams params;
        params.cacheRestForms = false;
        zone::Solver cpu(patch, plasticity::shipSteel(), params);
        gpu::ZoneGpuSolver device;
        std::string error;
        const bool ok = device.initialise(patch, cpu, plasticity::shipSteel(), params,
                                          SHIPSIM_SHADER_DIR, error);
        expectTrue("with no cached rest forms the GPU back-end declines",
                   !ok && error.find("cacheRestForms") != std::string::npos);
    }
    {
        zone::SolveParams params;
        params.plastic = false;
        zone::Solver cpu(patch, plasticity::shipSteel(), params);
        gpu::ZoneGpuSolver device;
        std::string error;
        const bool ok = device.initialise(patch, cpu, plasticity::shipSteel(), params,
                                          SHIPSIM_SHADER_DIR, error);
        expectTrue("and it declines the elastic path, saying why",
                   !ok && error.find("elastoplastic") != std::string::npos);
    }
    {
        // The fp64 enhanced-strain variants are compiled for the workgroup mapping
        // only. Asking for one with the other mapping has to be refused rather than
        // quietly served the float kernel: that would be a precision comparison in
        // which both sides were float, which is this repo's most common shape of
        // vacuous test and the exact mutant that survived §8's last pass.
        zone::SolveParams params;
        zone::Solver cpu(patch, plasticity::shipSteel(), params);
        gpu::ZoneGpuSolver device;
        std::string error;
        const bool ok = device.initialise(patch, cpu, plasticity::shipSteel(), params,
                                          SHIPSIM_SHADER_DIR, error, gpu::Mapping::Invocation,
                                          gpu::EasPrecision::Newton);
        expectTrue("and it declines an fp64 enhanced block on the one-invocation mapping",
                   !ok && error.find("workgroup mapping") != std::string::npos);
    }
}

// End to end on a small patch, if there is a device. Not a tolerance on the
// physics -- that is `zone_gpu_probe`'s job and it costs minutes -- but the
// statement that the device path runs, advances the state and comes back with
// something finite and in the right direction.
void testItRunsIfThereIsADevice(gpu::Mapping mapping) {
    solidshell::HexMesh flat = solidshell::makePlateMesh(0.4, 0.2, 0.01, 8, 4, 1);
    for (std::size_t n = 0; n < flat.nodeCount(); ++n) {
        const double x = flat.position[n * 3], y = flat.position[n * 3 + 1];
        if (x <= 1e-9 || x >= 0.4 - 1e-9 || y <= 1e-9 || y >= 0.2 - 1e-9)
            for (int k = 0; k < 3; ++k) flat.pin(n, k);
    }
    zone::Patch patch;
    patch.mesh = flat;
    patch.axis = Vec3{0, 0, 1};
    patch.right = Vec3{1, 0, 0};
    patch.up = Vec3{0, 1, 0};
    // `makePlateMesh` puts the mid-surface at z = 0, so the outer face is z > 0 --
    // not z > t/2, which is satisfied by no node at all and leaves the punch
    // gripping nothing. That was the first version of this fixture and it is the
    // reason `drivenNodes()` is asserted non-empty below.
    patch.centre = Vec3{0.2, 0.1, 0.0};
    patch.outerFace.assign(flat.nodeCount(), 0u);
    for (std::size_t n = 0; n < flat.nodeCount(); ++n)
        patch.outerFace[n] = flat.position[n * 3 + 2] > 0.0 ? 1u : 0u;
    patch.elementArea.assign(flat.elementCount(), 0.4 * 0.2 / 32.0);
    patch.panelOf.assign(flat.elementCount(), 0);
    patch.panels.push_back(0);
    patch.panelArea.push_back(0.08);
    patch.material = ah36Steel();
    patch.thickness = 0.01;
    patch.criticalTimestep =
        solidshell::criticalTimestep(flat, patch.material, solidshell::Formulation::SolidShell);

    zone::SolveParams params;
    params.indenter.halfLength = 0.05;
    params.indenter.halfWidth = 0.05;
    params.indenter.speed = 20.0;
    params.indenter.rampTime = 1.0e-4;
    params.indenter.stopAt = 1.0;

    zone::Solver seed(patch, plasticity::shipSteel(), params);
    // Guard: everything below compares against a punch that is supposed to be
    // pushing. A footprint that caught no free node makes every assertion pass on a
    // solver that did nothing.
    expectTrue("the punch grips some free nodes", !seed.drivenNodes().empty());
    gpu::ZoneGpuSolver device;
    std::string error;
    if (!device.initialise(patch, seed, plasticity::shipSteel(), params, SHIPSIM_SHADER_DIR,
                           error, mapping)) {
        std::printf("     skipped: %s\n", error.c_str());
        return;
    }
    std::printf("     device: %s\n", device.deviceName().c_str());

    const int steps = 3000;
    device.run(steps);
    const gpu::ZoneGpuState state = device.readback();

    bool finite = true;
    for (std::size_t i = 0; i < state.position.size(); ++i)
        finite = finite && std::isfinite(state.position[i]) && std::isfinite(state.velocity[i]);
    double peak = 0;
    for (const solidshell::ElementPlasticState& s : state.plastic)
        for (int gp = 0; gp < solidshell::kGauss; ++gp)
            peak = std::max(peak, s.point[gp].equivalentPlasticStrain);

    // **A gripped node is exactly where the punch is.** The kinematics are
    // prescribed, so this is a closed form and not an eyeballed number: every
    // driven node must have moved `state.penetration` along -axis and nothing
    // across it. It checks the punch condition, the velocity update and the
    // position update at once, and it is insensitive to whether the plate around it
    // has torn -- which the peak plastic strain below says it has.
    // Only the component **along the axis** is prescribed. The in-plane velocity of
    // a driven node is left alone, exactly as `zone::Solver::step` leaves it, so the
    // node is free to be dragged sideways by the plate it is part of -- and it is,
    // by 0.2 m once the bay round it has torn. Asserting it had not moved across
    // would be asserting a promise the indenter does not make.
    double worstAlong = 0, worstAcross = 0;
    for (std::uint32_t n : seed.drivenNodes()) {
        double along = 0, across = 0;
        for (int k = 0; k < 3; ++k) {
            const double d = state.position[n * 3 + static_cast<std::size_t>(k)] -
                             patch.mesh.position[n * 3 + static_cast<std::size_t>(k)];
            along += d * patch.axis[k];
            across += d * d;
        }
        worstAlong = std::max(worstAlong, std::abs(along + state.penetration));
        worstAcross = std::max(worstAcross, std::sqrt(std::max(0.0, across - along * along)));
    }
    std::printf("     %d steps, punch at %.4f m: %zu driven node(s) off by %.2e m along"
                " (%.2e across, unconstrained); peak eps_p %.4f, %d torn, work %.1f J\n",
                state.steps, state.penetration, seed.drivenNodes().size(), worstAlong,
                worstAcross, peak, state.tornElements, state.work);
    expectTrue("the device run stayed finite", finite);
    // Float carries ~7 digits, so 0.09 m of travel accumulated over 3 000 steps
    // lands within a few times 1e-7 m of the punch. Tight enough to catch a
    // dropped ramp, a doubled step or a lost constraint; loose enough not to be a
    // statement about float rounding.
    expectTrue("every gripped node is exactly where the punch is, to float",
               worstAlong < 5e-6);
    expectTrue("and the plate yielded, so the return map ran", peak > 1e-4);
    expectTrue("and it tore, so element deletion ran on the device too",
               state.tornElements > 0);
    expectTrue("and the punch did positive work on it", state.work > 0.0);
    expectTrue("and the punch went somewhere, so those are not four zeros",
               state.penetration > 1e-3);
}

// **The device against the CPU double reference, at unit scale.**
//
// The kinematic test above says the punch, the gather and the integrator work. It
// says nothing at all about the constitutive path, and mutation testing proved it:
// five single-edit mutants of the shader -- a wrong sqrt(2/3), a lost
// engineering-shear factor of two in the yield norm, the enhanced modes kept on a
// degraded element, a transposed rest Jacobian, and a host-side state offset --
// **all survived a suite that only asked whether the plate moved and yielded**.
//
// So this compares the physics, over a run short enough that float still tracks
// double. That qualification is not a hedge: `07-fem-spike-findings.md` §8 measures
// the two parting company over a full 5 500-step run, and the tolerances here are
// set from the agreement at 600 steps rather than from what would be nice.
void testItAgreesWithTheCpuOverAShortRun(gpu::Mapping mapping) {
    solidshell::HexMesh flat = solidshell::makePlateMesh(0.4, 0.2, 0.01, 8, 4, 1);
    // **Sheared in plane, and that is load-bearing.** On an axis-aligned
    // rectangular mesh the rest Jacobian is diagonal, so a kernel that read it
    // transposed would be indistinguishable from one that did not -- and mutation
    // testing confirmed that mutant surviving on the rectangular version of this
    // fixture. The shear is a function of y alone, so both faces of a node pair move
    // together and the element stays prismatic through its thickness, which is the
    // geometry `07-fem-spike-findings.md` §6 limit 1 requires.
    for (std::size_t n = 0; n < flat.nodeCount(); ++n) {
        const double x = flat.position[n * 3], y = flat.position[n * 3 + 1];
        if (x <= 1e-9 || x >= 0.4 - 1e-9 || y <= 1e-9 || y >= 0.2 - 1e-9)
            for (int k = 0; k < 3; ++k) flat.pin(n, k);
    }
    // Sheared *after* the perimeter is identified, so the boundary is still the
    // boundary. The shear depends on y alone, so both faces of a mid-surface node
    // move together and the element stays prismatic through its thickness.
    for (std::size_t n = 0; n < flat.nodeCount(); ++n)
        flat.position[n * 3] += 0.35 * flat.position[n * 3 + 1];
    zone::Patch patch;
    patch.mesh = flat;
    patch.axis = Vec3{0, 0, 1};
    patch.right = Vec3{1, 0, 0};
    patch.up = Vec3{0, 1, 0};
    patch.centre = Vec3{0.2 + 0.35 * 0.1, 0.1, 0.0};  // the sheared plate's centre
    patch.outerFace.assign(flat.nodeCount(), 0u);
    for (std::size_t n = 0; n < flat.nodeCount(); ++n)
        patch.outerFace[n] = flat.position[n * 3 + 2] > 0.0 ? 1u : 0u;
    patch.elementArea.assign(flat.elementCount(), 0.4 * 0.2 / 32.0);
    patch.panelOf.assign(flat.elementCount(), 0);
    patch.panels.push_back(0);
    patch.panelArea.push_back(0.08);
    patch.material = ah36Steel();
    patch.thickness = 0.01;
    patch.criticalTimestep =
        solidshell::criticalTimestep(flat, patch.material, solidshell::Formulation::SolidShell);

    // Slow and short on purpose. The plate has to yield -- otherwise there is no
    // constitutive path under comparison -- but it must stay far from failure,
    // because localisation is precisely where a float and a double solver stop
    // agreeing and the comparison would then be measuring §8's divergence rather
    // than the kernel.
    const int steps = 400;
    const double dt = patch.criticalTimestep * 1.0;  // timestepSafety / 0.9 at the default
    zone::SolveParams params;
    params.indenter.halfLength = 0.05;
    params.indenter.halfWidth = 0.05;
    params.indenter.speed = 4.0;
    params.indenter.rampTime = 1.0e-4;
    params.indenter.stopAt = 0.0;          // run to the step count
    params.duration = steps * dt;
    params.maxSteps = steps;

    zone::Solver seed(patch, plasticity::shipSteel(), params);
    gpu::ZoneGpuSolver device;
    std::string error;
    if (!device.initialise(patch, seed, plasticity::shipSteel(), params, SHIPSIM_SHADER_DIR,
                           error, mapping)) {
        std::printf("     skipped: %s\n", error.c_str());
        return;
    }
    device.run(steps);
    const gpu::ZoneGpuState state = device.readback();

    zone::Solver cpu(patch, plasticity::shipSteel(), params);
    const zone::SolveResult& reference = cpu.run();
    expectTrue("both paths ran the same number of steps", reference.steps == state.steps);

    // Integral quantities, per `07-fem-spike-findings.md` §2: never a max norm on a
    // nodal velocity, which reads over 200% between two correct explicit solvers.
    const double work = std::abs(state.work - reference.work) /
                        std::max(std::abs(reference.work), 1e-30);
    const double dissipation = std::abs(state.dissipation - reference.dissipation) /
                               std::max(std::abs(reference.dissipation), 1e-30);

    double sum = 0, travel = 0, peak = 0, worstStrain = 0;
    for (std::size_t i = 0; i < state.position.size(); ++i) {
        const double d = state.position[i] - cpu.position()[i];
        sum += d * d;
        travel = std::max(travel, std::abs(cpu.position()[i] - patch.mesh.position[i]));
    }
    const double rms = std::sqrt(sum / static_cast<double>(state.position.size()));
    for (std::size_t e = 0; e < patch.elementCount(); ++e)
        for (int gp = 0; gp < solidshell::kGauss; ++gp) {
            const double a = cpu.elementState()[e].point[gp].equivalentPlasticStrain;
            peak = std::max(peak, a);
            worstStrain =
                std::max(worstStrain, std::abs(a - state.plastic[e].point[gp].equivalentPlasticStrain));
        }

    std::printf("     %d steps: work %.3e, dissipation %.3e, node RMS %.2e m of %.4f m"
                " travelled, eps_p %.2e of a peak %.4f\n",
                steps, work, dissipation, rms, travel, worstStrain, peak);

    // **These tolerances are the float path's own noise floor, not a specification,
    // and the difference matters.** Measured here: work 9.3e-3, dissipation 5.8e-2,
    // displacement 1.2e-2 of the travel -- after four hundred steps of a centred,
    // barely-yielded plate, where the CPU reproduces itself to 1e-6 under the same
    // size of perturbation. They are set three to four times above that so the test
    // is not flaky, and what they can catch is bounded by it: an error in the yield
    // surface of a *percent* is invisible underneath 5.8e-2 of dissipation noise,
    // and mutation testing confirms one such mutant survives. That is a statement
    // about the float kernel and is recorded as one in
    // `07-fem-spike-findings.md` §8 rather than papered over here.
    expectTrue("the device's work agrees with the CPU's to 3e-2", work < 3.0e-2);
    expectTrue("and its plastic dissipation to 2e-1", dissipation < 2.0e-1);
    expectTrue("and the displacement field to 5e-2 of the travel", rms < 5.0e-2 * travel);
    // The plastic strain is deliberately *not* compared relatively. Just past first
    // yield it is a threshold quantity -- whether a point yielded at all turns on
    // crossing the surface -- so its relative difference is large for reasons that
    // have nothing to do with either solver being wrong.
    expectTrue("and the plastic strain absolutely, to 1e-2", worstStrain < 1.0e-2);
    // Guards. Two solvers that both did nothing agree perfectly, and a run that
    // never yielded compares no constitutive path at all.
    expectTrue("the run travelled somewhere", travel > 1e-4);
    expectTrue("and yielded, so the return map was compared", peak > 1e-3);
    expectTrue("and dissipated something to compare", reference.dissipation > 1.0);
}

// A degraded element -- one with any integration point failed -- must come back
// with **all seven** enhanced parameters at exactly zero, not with the values it
// held before it started tearing. It is an exact statement about the committed
// state rather than a tolerance on a trajectory, so it survives the divergence a
// long float run develops. Mutation testing put it here: keeping the enhanced
// modes on a degraded element passed everything.
void testADegradedElementReportsZeroedEnhancedParameters(gpu::Mapping mapping) {
    solidshell::HexMesh flat = solidshell::makePlateMesh(0.4, 0.2, 0.01, 8, 4, 1);
    for (std::size_t n = 0; n < flat.nodeCount(); ++n) {
        const double x = flat.position[n * 3], y = flat.position[n * 3 + 1];
        if (x <= 1e-9 || x >= 0.4 - 1e-9 || y <= 1e-9 || y >= 0.2 - 1e-9)
            for (int k = 0; k < 3; ++k) flat.pin(n, k);
    }
    zone::Patch patch;
    patch.mesh = flat;
    patch.axis = Vec3{0, 0, 1};
    patch.right = Vec3{1, 0, 0};
    patch.up = Vec3{0, 1, 0};
    patch.centre = Vec3{0.2, 0.1, 0.0};
    patch.outerFace.assign(flat.nodeCount(), 0u);
    for (std::size_t n = 0; n < flat.nodeCount(); ++n)
        patch.outerFace[n] = flat.position[n * 3 + 2] > 0.0 ? 1u : 0u;
    patch.elementArea.assign(flat.elementCount(), 0.4 * 0.2 / 32.0);
    patch.panelOf.assign(flat.elementCount(), 0);
    patch.panels.push_back(0);
    patch.panelArea.push_back(0.08);
    patch.material = ah36Steel();
    patch.thickness = 0.01;
    patch.criticalTimestep =
        solidshell::criticalTimestep(flat, patch.material, solidshell::Formulation::SolidShell);

    zone::SolveParams params;
    params.indenter.halfLength = 0.05;
    params.indenter.halfWidth = 0.05;
    params.indenter.speed = 20.0;
    params.indenter.rampTime = 1.0e-4;
    params.indenter.stopAt = 1.0;

    zone::Solver seed(patch, plasticity::shipSteel(), params);
    gpu::ZoneGpuSolver device;
    std::string error;
    if (!device.initialise(patch, seed, plasticity::shipSteel(), params, SHIPSIM_SHADER_DIR,
                           error, mapping))
        return;
    device.run(3000);
    const gpu::ZoneGpuState state = device.readback();

    int degraded = 0, live = 0, offenders = 0;
    double liveAlpha = 0;
    for (const solidshell::ElementPlasticState& s : state.plastic) {
        bool any = false;
        for (int gp = 0; gp < solidshell::kGauss; ++gp) any = any || s.point[gp].failed;
        if (any) {
            ++degraded;
            for (int k = 0; k < solidshell::kEas; ++k)
                if (s.enhanced[k] != 0.0) ++offenders;
        } else {
            ++live;
            for (int k = 0; k < solidshell::kEas; ++k)
                liveAlpha = std::max(liveAlpha, std::abs(s.enhanced[k]));
        }
    }
    std::printf("     %d degraded element(s), %d live; worst live |alpha| %.3e\n", degraded, live,
                liveAlpha);
    expectTrue("some elements degraded, so the rule had something to apply to", degraded > 0);
    expectTrue("every degraded element reports all seven enhanced parameters as zero",
               offenders == 0);
    // Guard: if the enhanced parameters were zero *everywhere* the check above is
    // satisfied by a kernel that never computes them at all.
    expectTrue("and the live elements carry non-zero ones, so zero means something",
               live > 0 && liveAlpha > 0.0);
}

// **The two mappings against each other**, which is a far sharper instrument than
// either against the CPU.
//
// The float-versus-double comparison above can only assert to the float path's own
// noise floor -- 5.8e-2 on the plastic dissipation after four hundred steps -- and
// `07-fem-spike-findings.md` §8 records that a 0.4% error in `kRoot23` hides
// underneath it. Two *float* kernels have no such floor between them: they run the
// same arithmetic on the same numbers in the same order, so they agree to a few ULP
// or one of them is wrong. That is what makes this the test that can see a mutant
// the CPU comparison cannot.
//
// It is not asserted as bit-identity, and the reason is measured rather than
// assumed: the two are different GLSL texts, so the driver's compiler is free to
// contract a multiply-add in one and not the other. Measured, the disagreement after
// 400 steps is **2.8e-5 on the plastic dissipation**, where the test above measures
// either kernel against the double CPU at 6.8e-2 on the same quantity and the same
// fixture. Three orders of magnitude between the two comparisons is the signature of
// the same arithmetic compiled twice rather than of a different computation, and it
// is the whole reason this test can assert what the other one cannot.
//
// **The material is given kinematic hardening on purpose.** `plasticity::shipSteel`
// leaves `kinematicModulus` at its default of zero, so the back stress is
// identically zero everywhere and every comparison of it is vacuous -- which is what
// the first version of this test discovered, by way of a guard that failed. Turning
// it on is the only condition under which the shader's back-stress arithmetic is
// exercised on the device at all.
void testTheTwoMappingsComputeTheSameElement() {
    solidshell::HexMesh flat = solidshell::makePlateMesh(0.4, 0.2, 0.01, 8, 4, 1);
    for (std::size_t n = 0; n < flat.nodeCount(); ++n) {
        const double x = flat.position[n * 3], y = flat.position[n * 3 + 1];
        if (x <= 1e-9 || x >= 0.4 - 1e-9 || y <= 1e-9 || y >= 0.2 - 1e-9)
            for (int k = 0; k < 3; ++k) flat.pin(n, k);
    }
    // Sheared, for the same reason the CPU comparison is: on an axis-aligned mesh the
    // rest Jacobian is diagonal and a transpose is not a change.
    for (std::size_t n = 0; n < flat.nodeCount(); ++n)
        flat.position[n * 3] += 0.35 * flat.position[n * 3 + 1];
    zone::Patch patch;
    patch.mesh = flat;
    patch.axis = Vec3{0, 0, 1};
    patch.right = Vec3{1, 0, 0};
    patch.up = Vec3{0, 1, 0};
    patch.centre = Vec3{0.2 + 0.35 * 0.1, 0.1, 0.0};
    patch.outerFace.assign(flat.nodeCount(), 0u);
    for (std::size_t n = 0; n < flat.nodeCount(); ++n)
        patch.outerFace[n] = flat.position[n * 3 + 2] > 0.0 ? 1u : 0u;
    patch.elementArea.assign(flat.elementCount(), 0.4 * 0.2 / 32.0);
    patch.panelOf.assign(flat.elementCount(), 0);
    patch.panels.push_back(0);
    patch.panelArea.push_back(0.08);
    patch.material = ah36Steel();
    patch.thickness = 0.01;
    patch.criticalTimestep =
        solidshell::criticalTimestep(flat, patch.material, solidshell::Formulation::SolidShell);

    // Short of localisation on purpose. Once a point crosses its damage limit the
    // element's response is a **threshold** event, so two runs that differ by one ULP
    // can tear on different steps and diverge macroscopically -- which is a true
    // statement about the physics and a useless one about the kernel. `zone_gpu_probe`
    // measures what happens past that point; this asserts the arithmetic.
    const int steps = 400;
    const double dt = patch.criticalTimestep * 1.0;
    zone::SolveParams params;
    params.indenter.halfLength = 0.05;
    params.indenter.halfWidth = 0.05;
    params.indenter.speed = 4.0;
    params.indenter.rampTime = 1.0e-4;
    params.indenter.stopAt = 0.0;
    params.duration = steps * dt;
    params.maxSteps = steps;

    // Kinematic hardening on, so the back stress is a live quantity. 2 GPa is about
    // 1% of this steel's Young's modulus, which is the order the Armstrong-Frederick
    // literature puts it at and, more to the point here, large enough that the back
    // stress reaches a few tens of MPa rather than staying at the zero the default
    // pins it to.
    plasticity::Material material = plasticity::shipSteel();
    material.flow.kinematicModulus = 2.0e9;

    // Three runs, not two: the third repeats the workgroup mapping. **That repeat is
    // the only test in this suite that can see a missing `barrier()`**, and it is the
    // one failure mode the remapped kernel has that the one-invocation kernel cannot.
    // A race between the Gauss-point threads and the Kaa threads would leave the
    // answer dependent on warp scheduling, so it would differ run to run while still
    // being perfectly plausible -- and every comparison against the CPU, and against
    // the other mapping, would pass on whichever of the two answers it happened to
    // get. Bit-identity across two independent submissions is what rules it out, and
    // it is asserted as bit-identity rather than a tolerance because a correctly
    // synchronised kernel has no reason to differ by even one bit.
    gpu::ZoneGpuState got[3];
    std::string shader[3];
    const gpu::Mapping mappings[3] = {gpu::Mapping::Invocation, gpu::Mapping::Workgroup,
                                      gpu::Mapping::Workgroup};
    for (int which = 0; which < 3; ++which) {
        zone::Solver seed(patch, material, params);
        gpu::ZoneGpuSolver device;
        std::string error;
        if (!device.initialise(patch, seed, material, params, SHIPSIM_SHADER_DIR,
                               error, mappings[which])) {
            std::printf("     skipped: %s\n", error.c_str());
            return;
        }
        device.run(steps);
        got[which] = device.readback();
        shader[which] = device.elementShader();
    }

    // **The guard this whole test rests on, and it was missing.** Everything below
    // compares two mappings; if `Mapping` ever stopped selecting between two shaders
    // and both loaded the same one, the two runs would be the same run and every
    // tolerance would pass by construction. Mutation testing found exactly that --
    // collapsing the host's ternary to a single filename survived the entire suite --
    // which is `CLAUDE.md`'s "two solvers that both did nothing agree perfectly" in
    // its most literal form. Asserted first, so a failure here explains the rest.
    expectTrue("the two mappings really did load different element shaders",
               shader[0] != shader[1] && !shader[0].empty() && !shader[1].empty());
    expectTrue("and the repeated run loaded the same one it is being compared against",
               shader[1] == shader[2]);

    bool repeatable = got[1].position == got[2].position &&
                      got[1].velocity == got[2].velocity &&
                      got[1].work == got[2].work && got[1].dissipation == got[2].dissipation &&
                      got[1].tornElements == got[2].tornElements;
    for (std::size_t e = 0; e < got[1].plastic.size() && repeatable; ++e) {
        for (int gp = 0; gp < solidshell::kGauss; ++gp) {
            const plasticity::State& a = got[1].plastic[e].point[gp];
            const plasticity::State& b = got[2].plastic[e].point[gp];
            repeatable = repeatable && a.equivalentPlasticStrain == b.equivalentPlasticStrain &&
                         a.damage == b.damage && a.failed == b.failed;
            for (int i = 0; i < 6; ++i)
                repeatable = repeatable && a.plasticStrain[i] == b.plasticStrain[i] &&
                             a.backStress[i] == b.backStress[i];
        }
        for (int k = 0; k < solidshell::kEas; ++k)
            repeatable = repeatable && got[1].plastic[e].enhanced[k] == got[2].plastic[e].enhanced[k];
    }
    expectTrue("the workgroup kernel is bit-identical run to run, so nothing races a barrier",
               repeatable);

    // Every quantity the device owns, not just the ones the energy account reads: a
    // mutant that got the back stress wrong while leaving the dissipation right is
    // exactly the shape this repo keeps finding.
    double worstPosition = 0, worstVelocity = 0, worstStrain = 0, worstBack = 0;
    double worstAlpha = 0, worstDamage = 0;
    double scalePosition = 0, scaleBack = 0, scaleAlpha = 0, scaleVelocity = 0;
    double velocitySum = 0, velocityScaleSum = 0;
    for (std::size_t i = 0; i < got[0].position.size(); ++i) {
        worstPosition = std::max(worstPosition, std::abs(got[0].position[i] - got[1].position[i]));
        worstVelocity = std::max(worstVelocity, std::abs(got[0].velocity[i] - got[1].velocity[i]));
        scaleVelocity = std::max(scaleVelocity, std::abs(got[0].velocity[i]));
        const double dv = got[0].velocity[i] - got[1].velocity[i];
        velocitySum += dv * dv;
        velocityScaleSum += got[0].velocity[i] * got[0].velocity[i];
        scalePosition = std::max(scalePosition,
                                 std::abs(got[0].position[i] - patch.mesh.position[i]));
    }
    const double velocityRms = std::sqrt(velocitySum / std::max<double>(got[0].velocity.size(), 1));
    const double velocityScaleRms =
        std::sqrt(velocityScaleSum / std::max<double>(got[0].velocity.size(), 1));
    for (std::size_t e = 0; e < got[0].plastic.size(); ++e) {
        for (int gp = 0; gp < solidshell::kGauss; ++gp) {
            const plasticity::State& a = got[0].plastic[e].point[gp];
            const plasticity::State& b = got[1].plastic[e].point[gp];
            worstStrain = std::max(worstStrain, std::abs(a.equivalentPlasticStrain -
                                                         b.equivalentPlasticStrain));
            worstDamage = std::max(worstDamage, std::abs(a.damage - b.damage));
            for (int i = 0; i < 6; ++i) {
                worstBack = std::max(worstBack, std::abs(a.backStress[i] - b.backStress[i]));
                scaleBack = std::max(scaleBack, std::abs(a.backStress[i]));
            }
        }
        for (int k = 0; k < solidshell::kEas; ++k) {
            worstAlpha = std::max(worstAlpha, std::abs(got[0].plastic[e].enhanced[k] -
                                                      got[1].plastic[e].enhanced[k]));
            scaleAlpha = std::max(scaleAlpha, std::abs(got[0].plastic[e].enhanced[k]));
        }
    }
    const double work = std::abs(got[0].work - got[1].work) / std::max(std::abs(got[0].work), 1e-30);
    const double dissipation = std::abs(got[0].dissipation - got[1].dissipation) /
                               std::max(std::abs(got[0].dissipation), 1e-30);
    std::printf("     invocation vs workgroup after %d steps: work %.2e, dissipation %.2e,"
                " position %.2e m of %.4f, alpha %.2e of %.2e, back stress %.2e of %.2e,"
                " eps_p %.2e, damage %.2e, velocity RMS %.2e of %.3f m/s (max %.2e of %.3f)\n",
                steps, work, dissipation, worstPosition, scalePosition, worstAlpha, scaleAlpha,
                worstBack, scaleBack, worstStrain, worstDamage, velocityRms, velocityScaleRms,
                worstVelocity, scaleVelocity);

    // **Tolerances set from what was measured, not from what would be comfortable**,
    // per CLAUDE.md's "a loose assertion is nearly a vacuous one". Measured on a
    // 1070 Ti after these 400 steps: work 4.0e-4, dissipation 2.6e-5, position 3.4e-5
    // of the travel, alpha 1.6e-4 of the peak, eps_p 3.0e-6, damage 2.7e-5. Each
    // bound below is that figure rounded up by about a factor of three -- enough that
    // a driver update which re-contracts one shader's multiply-adds does not turn the
    // suite red, and no more.
    //
    // What that buys, which is the point of the test: the same fixture compared
    // against the **double CPU** can only assert dissipation to 2e-1, because that is
    // the float path's own noise floor. Here it is 1e-4 -- two thousand times
    // sharper, on the identical quantity -- and that is the margin in which the five
    // shader mutants §8 records as surviving the device suite now die.
    expectTrue("the two mappings agree on the work to 2e-3", work < 2.0e-3);
    expectTrue("and on the plastic dissipation to 1e-4", dissipation < 1.0e-4);
    expectTrue("and on every node position to 2e-4 of the travel",
               worstPosition < 2.0e-4 * std::max(scalePosition, 1e-30));
    // **The nodal velocity is compared in RMS and deliberately not in max norm.**
    // Measured here, the two kernels' worst single nodal velocity differs by 7.9e-3
    // of the peak while their worst node *position* differs by 4.1e-5 of the travel
    // -- two hundred times better. That is not one of them being wrong: it is
    // `07-fem-spike-findings.md` §2's finding that an explicit scheme carries no
    // accuracy at all in the modes at its stability boundary, and velocity is where
    // those modes live. §2 records two *correct* explicit solvers reading over 200%
    // apart on a peak nodal velocity. Asserting a max norm on it would be asserting
    // a property neither kernel has, so what is asserted is the RMS.
    //
    // Even in RMS it is the loosest bound in this test -- measured at 7.3e-3 of the
    // RMS speed, against 4.1e-5 for the position field. **That ratio is itself the
    // evidence**: position is the 400-step integral of velocity, so a difference this
    // much smaller in the integral than in the integrand is a difference living in
    // modes that oscillate at close to 2*dt and cancel when summed. A difference in
    // anything the element actually computes would show up in both. It is asserted at
    // 2e-2 and named here as the weakest assertion in the file, because the
    // load-bearing ones are the six above it.
    expectTrue("and on the RMS nodal velocity to 2e-2",
               velocityRms < 2.0e-2 * std::max(velocityScaleRms, 1e-30));
    expectTrue("and on every enhanced parameter to 1e-3 of the peak",
               worstAlpha < 1.0e-3 * std::max(scaleAlpha, 1e-30));
    expectTrue("and on every back stress to 1e-3 of the peak",
               worstBack < 1.0e-3 * std::max(scaleBack, 1e-30));
    expectTrue("and on every accumulated plastic strain to 1e-5", worstStrain < 1.0e-5);
    expectTrue("and on every damage to 1e-4", worstDamage < 1.0e-4);
    expectTrue("and they report the same torn count",
               got[0].tornElements == got[1].tornElements);

    // **Vacuity guards, because two kernels that both did nothing agree perfectly.**
    // Each names a quantity that would be identically zero on a do-nothing kernel, so
    // every field compared above has to have had something in it to compare.
    expectTrue("the run moved the patch, so the positions compared are not the rest ones",
               scalePosition > 1e-4);
    expectTrue("and the nodes were moving, so the velocities compared are not zeros",
               scaleVelocity > 1e-3);
    expectTrue("and it yielded, so the return map and the back stress were compared",
               scaleBack > 1.0);
    expectTrue("and the enhanced modes were live, so the condensation was compared",
               scaleAlpha > 0.0);
    expectTrue("and it did work, so the integrator was compared", got[0].work > 0.0);
}

// **The enhanced block in fp64, and what actually moves when it is.**
//
// `07-fem-spike-findings.md` §8's last open item was "keep alpha in double -- the only
// part that has been shown to need the digits". Five kernels are compared here, all
// compiled from one GLSL source so that nothing but the precision differs:
//
//   Float        as shipped
//   FloatTight   float arithmetic, the CPU's stopping rule (1e-16 sigma_y V, 40 its)
//   Solve        the 7x7 equilibration, Cholesky and substitutions in fp64
//   Condense     + Kaa and the residual accumulated in fp64 over a fp64 G
//   Newton       + alpha, its correction and its persistent state, and the CPU's rule
//
// The pairing is what makes this a measurement rather than a demonstration. `Newton`
// changes the tolerance and the arithmetic together, so `FloatTight` is the control
// that holds the tolerance fixed and asks what the *arithmetic* was worth; `Solve` and
// `Condense` hold the tolerance at the shipped one and ask the same question from the
// other side.
//
// Two regimes, because they say different things.
void testTheEnhancedBlockInDouble() {
    solidshell::HexMesh flat = solidshell::makePlateMesh(0.4, 0.2, 0.01, 8, 4, 1);
    for (std::size_t n = 0; n < flat.nodeCount(); ++n) {
        const double x = flat.position[n * 3], y = flat.position[n * 3 + 1];
        if (x <= 1e-9 || x >= 0.4 - 1e-9 || y <= 1e-9 || y >= 0.2 - 1e-9)
            for (int k = 0; k < 3; ++k) flat.pin(n, k);
    }
    // Sheared, for the same reason the other comparisons are: on an axis-aligned mesh
    // the rest Jacobian is diagonal and a transpose is not a change.
    for (std::size_t n = 0; n < flat.nodeCount(); ++n)
        flat.position[n * 3] += 0.35 * flat.position[n * 3 + 1];
    zone::Patch patch;
    patch.mesh = flat;
    patch.axis = Vec3{0, 0, 1};
    patch.right = Vec3{1, 0, 0};
    patch.up = Vec3{0, 1, 0};
    patch.centre = Vec3{0.2 + 0.35 * 0.1, 0.1, 0.0};
    patch.outerFace.assign(flat.nodeCount(), 0u);
    for (std::size_t n = 0; n < flat.nodeCount(); ++n)
        patch.outerFace[n] = flat.position[n * 3 + 2] > 0.0 ? 1u : 0u;
    patch.elementArea.assign(flat.elementCount(), 0.4 * 0.2 / 32.0);
    patch.panelOf.assign(flat.elementCount(), 0);
    patch.panels.push_back(0);
    patch.panelArea.push_back(0.08);
    patch.material = ah36Steel();
    patch.thickness = 0.01;
    patch.criticalTimestep =
        solidshell::criticalTimestep(flat, patch.material, solidshell::Formulation::SolidShell);

    const int steps = 400;
    const double dt = patch.criticalTimestep * 1.0;
    const gpu::EasPrecision ladder[5] = {
        gpu::EasPrecision::Float, gpu::EasPrecision::FloatTight, gpu::EasPrecision::Solve,
        gpu::EasPrecision::Condense, gpu::EasPrecision::Newton};
    const char* names[5] = {"float", "tight", "solve", "condense", "fp64"};

    // What one run of the ladder measures: the peak enhanced parameter each kernel
    // committed, and the worst it is off the double CPU by.
    struct Ladder {
        double peak[5] = {0, 0, 0, 0, 0};
        double againstCpu[5] = {0, 0, 0, 0, 0};
        double cpuPeak = 0;
        double travel = 0;
        double cpuDissipation = 0;
        double dissipation[5] = {0, 0, 0, 0, 0};
        double alphaSpread[5] = {0, 0, 0, 0, 0};   // against kernel 0's own answer
        double positionSpread[5] = {0, 0, 0, 0, 0};
        bool ran = false;
        std::string shader[5];
    };
    const auto sweep = [&](double speed, int steps) {
        Ladder out;
        zone::SolveParams params;
        params.indenter.halfLength = 0.05;
        params.indenter.halfWidth = 0.05;
        params.indenter.speed = speed;
        params.indenter.rampTime = 1.0e-4;
        params.indenter.stopAt = 0.0;
        params.duration = steps * dt;
        params.maxSteps = steps;

        gpu::ZoneGpuState got[5];
        for (int which = 0; which < 5; ++which) {
            zone::Solver seed(patch, plasticity::shipSteel(), params);
            gpu::ZoneGpuSolver device;
            std::string error;
            if (!device.initialise(patch, seed, plasticity::shipSteel(), params,
                                   SHIPSIM_SHADER_DIR, error, gpu::Mapping::Workgroup,
                                   ladder[which])) {
                std::printf("     skipped: %s\n", error.c_str());
                return out;
            }
            device.run(steps);
            got[which] = device.readback();
            out.shader[which] = device.elementShader();
            out.dissipation[which] = got[which].dissipation;
        }
        zone::Solver cpu(patch, plasticity::shipSteel(), params);
        out.cpuDissipation = cpu.run().dissipation;
        for (std::size_t e = 0; e < patch.elementCount(); ++e)
            for (int k = 0; k < solidshell::kEas; ++k) {
                const double want = cpu.elementState()[e].enhanced[k];
                out.cpuPeak = std::max(out.cpuPeak, std::abs(want));
                for (int which = 0; which < 5; ++which) {
                    const double a = got[which].plastic[e].enhanced[k];
                    out.peak[which] = std::max(out.peak[which], std::abs(a));
                    out.againstCpu[which] = std::max(out.againstCpu[which], std::abs(a - want));
                    out.alphaSpread[which] =
                        std::max(out.alphaSpread[which],
                                 std::abs(a - got[which == 4 ? 1 : 0].plastic[e].enhanced[k]));
                }
            }
        for (std::size_t i = 0; i < got[0].position.size(); ++i) {
            out.travel = std::max(out.travel, std::abs(got[0].position[i] - patch.mesh.position[i]));
            for (int which = 0; which < 5; ++which)
                out.positionSpread[which] =
                    std::max(out.positionSpread[which],
                             std::abs(got[which].position[i] - got[which == 4 ? 1 : 0].position[i]));
        }
        out.ran = true;
        return out;
    };

    // --- Regime one: lightly loaded, and the enhanced modes are simply off ----------
    //
    // At 0.005 m/s over 120 steps the plate stays elastic and `|delta . r| <= 1e-9
    // sigma_y V` is met on the **first** iteration of every element. The shader stops
    // *before applying* that first correction, so alpha never leaves its warm start of
    // zero and the enhanced modes are, on the device, switched off.
    //
    // **That is the mechanism behind §8's "the GPU's enhanced parameters were exactly
    // zero", and it is not a precision failure.** `Solve` and `Condense` compute the
    // same correction in fp64 and then discard it at the same gate, so they return the
    // identical bit-zero. Only the stopping rule turns them back on, and once it does
    // float alone gets within 20% of the CPU's alpha on the worst parameter and 5.5% on
    // the peak -- which is the ceiling the float tangent and stress impose, and not one
    // that fp64 lifts: the fp64 kernel lands on the same figure.
    const Ladder light = sweep(0.005, 120);
    if (!light.ran) return;
    std::printf("     lightly loaded, 120 steps: cpu |alpha| %.3e; device", light.cpuPeak);
    for (int which = 0; which < 5; ++which)
        std::printf(" %s %.3e", names[which], light.peak[which]);
    std::printf(" (worst off the cpu: float %.2e, tight %.2e, fp64 %.2e)\n",
                light.againstCpu[0], light.againstCpu[1], light.againstCpu[4]);

    expectTrue("the five enhanced-block precisions loaded five different kernels",
               light.shader[0] != light.shader[1] && light.shader[1] != light.shader[2] &&
                   light.shader[2] != light.shader[3] && light.shader[3] != light.shader[4] &&
                   light.shader[0] != light.shader[4] && !light.shader[0].empty());
    // Exact, not a tolerance: bit-zero on all seven parameters of all thirty-two
    // elements is not something rounding produces by accident.
    expectTrue("under the shipped work gate the device's enhanced modes are exactly off",
               light.peak[0] == 0.0);
    expectTrue("and fp64 in the 7x7 solve does not turn them on", light.peak[2] == 0.0);
    expectTrue("nor fp64 in the condensation", light.peak[3] == 0.0);
    expectTrue("but the CPU's own stopping rule does, in float", light.peak[1] > 0.0);
    expectTrue("and in fp64", light.peak[4] > 0.0);
    // The guard that stops those three zeros from being a run that did nothing: the
    // double reference solved the same element and got enhanced parameters out of it.
    expectTrue("and the CPU's are non-zero, so that zero is the gate and not the fixture",
               light.cpuPeak > 0.0);
    expectTrue("and the light run still deformed the patch", light.travel > 1e-8);
    // **Quantified, so "off" and "on" are a difference and not two adjectives.** With
    // the modes off the worst enhanced parameter is wrong by the whole of the CPU's --
    // exactly, since it is zero. With the rule tightened, float alone brings that to
    // 1.10e-9 of a 5.51e-9 peak, i.e. 20.0%, and the peak itself lands within 5.5% of
    // the CPU's. Asserted at 30%: the measurement with room for a driver that contracts
    // a multiply-add differently, and still a fifth of the 100% above it.
    expectTrue("with the modes off the device's alpha is wrong by the whole of the CPU's",
               light.againstCpu[0] == light.cpuPeak);
    expectTrue("and with the rule tightened, float alone recovers it to 30%",
               light.againstCpu[1] < 0.3 * light.cpuPeak);
    expectTrue("and fp64 recovers it no better than float does",
               light.againstCpu[4] >= 0.5 * light.againstCpu[1]);

    // --- Regime two: loaded and yielding, where the modes are live on every kernel ---
    //
    // Here the gate never fires early and all five kernels carry live enhanced modes.
    // What they show is the ceiling: every one of them is off the double CPU by
    // 3.8e-6 of an alpha of 3.6e-5, and the *whole* fp64 enhanced block moves alpha by
    // 4.6e-9 against its float twin -- eight hundred times less than the error it was
    // supposed to fix. That ratio is the finding, and it is asserted as a ratio so
    // that it fails if fp64 ever does start mattering here.
    const Ladder loaded = sweep(4.0, steps);
    if (!loaded.ran) return;
    std::printf("     loaded, %d steps: cpu |alpha| %.3e; device", steps, loaded.cpuPeak);
    for (int which = 0; which < 5; ++which)
        std::printf(" %s %.3e (off cpu %.2e, off its float twin %.2e)", names[which],
                    loaded.peak[which], loaded.againstCpu[which], loaded.alphaSpread[which]);
    std::printf("; worst node against the float twin");
    for (int which = 0; which < 5; ++which)
        std::printf(" %.2e", loaded.positionSpread[which]);
    std::printf(" of %.4f m\n", loaded.travel);

    expectTrue("with the plate loaded every kernel's enhanced modes are live",
               loaded.peak[0] > 0.0 && loaded.peak[4] > 0.0 && loaded.cpuPeak > 0.0);
    // Measured: 4.03e-9, 4.24e-9 and 4.64e-9 against an error of 3.77e-6, so ratios of
    // 940, 890 and 810. Asserted at 100, which is a factor of eight of headroom.
    expectTrue("fp64 in the 7x7 solve moves alpha by under 1% of how wrong it already is",
               loaded.alphaSpread[2] * 100.0 < loaded.againstCpu[0]);
    expectTrue("and fp64 in the condensation likewise",
               loaded.alphaSpread[3] * 100.0 < loaded.againstCpu[0]);
    expectTrue("and the whole enhanced block in fp64 likewise, against its float twin",
               loaded.alphaSpread[4] * 100.0 < loaded.againstCpu[1]);
    // And the same statement on the field the energies are integrals of: 1.04e-7 m of
    // 2.2e-3 m travelled, asserted at ten times that.
    expectTrue("and it moves no node position by more than 5e-4 of the travel",
               loaded.positionSpread[4] < 5.0e-4 * std::max(loaded.travel, 1e-30));
    // Vacuity. Every ratio above is a difference divided by a difference, and both
    // halves have to be real: the kernels have to disagree with the CPU at all, and the
    // fp64 one has to have differed from its float twin rather than being the same
    // pipeline twice.
    expectTrue("the kernels really do differ from the CPU, so those ratios have a "
               "denominator", loaded.againstCpu[0] > 0.0 && loaded.againstCpu[1] > 0.0);
    expectTrue("and fp64 really did change the answer, so they have a numerator",
               loaded.alphaSpread[4] > 0.0 && loaded.alphaSpread[2] > 0.0);
    expectTrue("and the loaded run yielded, so the return map was exercised",
               loaded.cpuDissipation > 1.0 && loaded.dissipation[4] > 1.0);
}

}  // namespace

void runZoneGpuTests() {
    std::printf("\n--- the solid-shell zone on the GPU ---\n");
    testTheShadersAndTheHostAgreeOnTheLayout();
    testTheShadersAndTheHostAgreeOnTheirDispatch();
    testItDeclinesWhatItCannotDo();
    // **Every device test runs on both mappings.** They are two kernels, and a suite
    // that exercised only the default would leave the other one shipping untested --
    // which is `CLAUDE.md`'s "two functions on the caller's own path shipped with no
    // test at all, in a commit whose headline feature was well tested".
    const gpu::Mapping mappings[2] = {gpu::Mapping::Invocation, gpu::Mapping::Workgroup};
    const char* names[2] = {"one invocation per element", "one workgroup per element"};
    for (int which = 0; which < 2; ++which) {
        std::printf("   [%s]\n", names[which]);
        testItRunsIfThereIsADevice(mappings[which]);
        testItAgreesWithTheCpuOverAShortRun(mappings[which]);
        testADegradedElementReportsZeroedEnhancedParameters(mappings[which]);
    }
    testTheTwoMappingsComputeTheSameElement();
    testTheEnhancedBlockInDouble();
}
