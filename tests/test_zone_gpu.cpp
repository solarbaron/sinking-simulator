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

void testTheShaderAndTheHostAgreeOnTheLayout() {
    const std::string path =
        std::string(SHIPSIM_SHADER_SOURCE_DIR) + "/solidshell_forces.comp";
    std::ifstream file(path);
    expectTrue("the solid-shell force shader source is readable", file.good());
    if (!file.good()) return;
    const std::string source((std::istreambuf_iterator<char>(file)),
                             std::istreambuf_iterator<char>());

    // Derived here from the element's own constants, so a change to kGauss, kDof or
    // kEas moves this side and the shader has to follow.
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
    for (const Check& c : checks) {
        long long got = -1;
        const bool found = constantFromShader(source, c.name, got);
        expectTrue(std::string("the shader declares ") + c.name, found);
        expectTrue(std::string("and the shader's ") + c.name + " matches the element's",
                   found && got == c.want);
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
}

// End to end on a small patch, if there is a device. Not a tolerance on the
// physics -- that is `zone_gpu_probe`'s job and it costs minutes -- but the
// statement that the device path runs, advances the state and comes back with
// something finite and in the right direction.
void testItRunsIfThereIsADevice() {
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
                           error)) {
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
void testItAgreesWithTheCpuOverAShortRun() {
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
                           error)) {
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
void testADegradedElementReportsZeroedEnhancedParameters() {
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
                           error))
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

}  // namespace

void runZoneGpuTests() {
    std::printf("\n--- the solid-shell zone on the GPU ---\n");
    testTheShaderAndTheHostAgreeOnTheLayout();
    testItDeclinesWhatItCannotDo();
    testItRunsIfThereIsADevice();
    testItAgreesWithTheCpuOverAShortRun();
    testADegradedElementReportsZeroedEnhancedParameters();
}
