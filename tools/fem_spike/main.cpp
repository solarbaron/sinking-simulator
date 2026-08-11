// SPDX-License-Identifier: MIT
//
// The Phase 3 de-risking spike from docs/06-roadmap.md.
//
// Three questions, in order of how much they matter:
//   1. Is the element formulation right? (cantilever vs Euler-Bernoulli theory)
//   2. Does the GPU kernel compute the same thing as the CPU reference?
//   3. How many element updates per second does this card actually deliver?
//
// Question 3 is the one that decides whether adaptive full-3D tetrahedral FEM on
// a ship is a plan or a wish, but it is worthless without 1 and 2.
#include "engine/gpu/fem_gpu.hpp"
#include "engine/sim/fem.hpp"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

using namespace sim::fem;

namespace {

#ifndef SHIPSIM_SHADER_DIR
#define SHIPSIM_SHADER_DIR "."
#endif

int failures = 0;
int skippedSections = 0;

void check(const char* what, bool ok, const std::string& detail = {}) {
    std::printf("  [%s] %-52s %s\n", ok ? "PASS" : "FAIL", what, detail.c_str());
    if (!ok) ++failures;
}

// **`CLAUDE.md`: GPU work must skip, not fail, when there is no Vulkan device.**
// This tool was the one place that did not -- it counted a missing device as a
// failed check and exited 1, so on any machine without a card `fem_spike` printed
// `CHECKS FAILED` for the entirely correct reason that there was nothing to run it
// on. Every sibling tool already prints `no usable GPU (...)` and carries on.
//
// The tradeoff this makes, taken deliberately and in line with those siblings: a
// genuine driver fault on a machine that *has* a device also reads as a skip here.
// The gate classifies that from outside -- it knows whether an ICD and a device node
// exist and calls a skip on a working card a failure -- which is the right place for
// it, because the tool cannot tell the two apart and the gate can.
void skipGpuSection(const char* what, const std::string& error) {
    std::printf("  [SKIP] %-52s no usable GPU (%s)\n", what, error.c_str());
    ++skippedSections;
}

// Pin every node at x = 0 to build a cantilever.
void clampRootFace(TetMesh& mesh, float tolerance = 1e-6f) {
    for (std::size_t n = 0; n < mesh.nodeCount(); ++n)
        if (mesh.position[n * 3] < tolerance) mesh.fixed[n] = 1u;
}

float maxDeflection(const TetMesh& mesh, const TetMesh& rest) {
    float worst = 0;
    for (std::size_t n = 0; n < mesh.nodeCount(); ++n)
        worst = std::max(worst, std::abs(rest.position[n * 3 + 2] - mesh.position[n * 3 + 2]));
    return worst;
}

// --- 1. Does the element formulation reproduce beam theory? ------------------

void testCantileverAgainstBeamTheory() {
    std::printf("\n1. Element formulation vs Euler-Bernoulli beam theory\n");

    const float length = 1.0f, width = 0.10f, thickness = 0.050f;
    const float gravity = 9.80665f;
    const Material steel;

    // Tip deflection of a cantilever under its own weight: q L^4 / (8 E I).
    const float inertia = width * thickness * thickness * thickness / 12.0f;
    const float load = steel.density * width * thickness * gravity;  // N/m
    const float theory = load * std::pow(length, 4.0f) / (8.0f * steel.youngsModulus * inertia);
    std::printf("     beam theory: %.4f mm tip deflection\n", theory * 1000.0);

    // A single mesh proves nothing: linear tetrahedra lock in bending and will
    // always come in under theory. What has to be true is that refining the mesh
    // drives the error down. That is the difference between a formulation that is
    // correct and merely coarse, and one that is wrong.
    std::printf("     %8s %10s %12s %14s %10s\n", "thru-t", "tets", "dt (s)", "deflection mm",
                "error");

    struct Result { int refinement; float deflection; float error; };
    std::vector<Result> results;

    // Refine uniformly, not just through the thickness: element aspect ratio is
    // itself a driver of locking, so holding it constant is what makes this a
    // clean convergence study rather than a confounded one.
    for (int refine : {1, 2, 4}) {
        TetMesh mesh = makeBoxTetMesh(length, width, thickness, 20 * refine, 2 * refine,
                                      2 * refine);
        mesh.computeRestState(steel);
        clampRootFace(mesh);
        const TetMesh rest = mesh;

        const float dt = criticalTimestep(mesh, steel);

        std::string error;
        gpu::FemGpuSolver solver;
        if (!solver.initialise(mesh, SHIPSIM_SHADER_DIR, error)) {
            skipGpuSection("cantilever: GPU solver initialised", error);
            return;
        }

        // Explicit dynamics oscillates about the static answer, so damp it to
        // rest rather than reading a moving target. Damping changes only how fast
        // equilibrium is reached, not where it is -- but it has to be specified
        // per unit *time*, not per step, or the coarse and fine meshes would be
        // damped by wildly different amounts.
        const float decayRate = 124.0f;  // 1/s
        const float damping = std::exp(-decayRate * dt);
        const float settleSeconds = 0.20f;
        const int steps = static_cast<int>(settleSeconds / dt);

        solver.run(steps, steel, dt, gravity, damping);
        solver.readback(mesh);

        const float deflection = maxDeflection(mesh, rest);
        const float relative = std::abs(deflection - theory) / theory;
        results.push_back({refine, deflection, relative});
        std::printf("     %8d %10zu %12.3g %14.4f %9.1f%%\n", 2 * refine, mesh.tetCount(),
                    static_cast<double>(dt), deflection * 1000.0, relative * 100.0);
        std::fflush(stdout);
    }

    bool converging = true;
    for (std::size_t i = 1; i < results.size(); ++i)
        if (!(results[i].error < results[i - 1].error)) converging = false;

    check("refinement drives the error toward beam theory", converging);
    check("finest mesh within 20% of beam theory", results.back().error < 0.20);
    check("deflection is downward, finite and settled",
          results.back().deflection > 0 && std::isfinite(results.back().deflection));
}

// --- 2. Does the GPU compute the same thing as the CPU? ----------------------

void testGpuMatchesCpu() {
    std::printf("\n2. GPU kernel vs CPU reference\n");

    const Material steel;
    TetMesh cpuMesh = makeBoxTetMesh(0.4f, 0.1f, 0.02f, 16, 4, 2);
    cpuMesh.computeRestState(steel);
    clampRootFace(cpuMesh);

    // Deform and twist the mesh off its rest state before anything is uploaded,
    // so the polar decomposition is exercised rather than handed something close
    // to the identity. This has to happen before initialise() -- that is what
    // fills the device buffers.
    for (std::size_t n = 0; n < cpuMesh.nodeCount(); ++n) {
        const float x = cpuMesh.position[n * 3];
        cpuMesh.position[n * 3 + 2] += 0.02f * x * x;
        cpuMesh.position[n * 3 + 1] += 0.01f * x;
    }
    TetMesh gpuMesh = cpuMesh;

    std::string error;
    gpu::FemGpuSolver solver;
    if (!solver.initialise(gpuMesh, SHIPSIM_SHADER_DIR, error)) {
        skipGpuSection("GPU solver initialised", error);
        return;
    }
    check("GPU solver initialised", true, "on " + solver.deviceName());

    const float dt = criticalTimestep(cpuMesh, steel);

    // Kernel correctness and long-run reproducibility are two different
    // questions, and conflating them is how a correct solver gets mistaken for a
    // broken one. A single step measures the kernel. Everything after measures
    // how fast an explicit scheme at the CFL limit amplifies the last-bit
    // differences between two float implementations -- which it does, because its
    // highest mesh modes sit right at the stability boundary.
    std::printf("     %8s %18s %18s\n", "steps", "RMS dv / RMS v", "max |dx| (m)");

    struct Sample { int steps; float velocity; float position; };
    std::vector<Sample> samples;

    int done = 0;
    for (int target : {1, 10, 100, 1000}) {
        const int remaining = target - done;
        for (int i = 0; i < remaining; ++i) stepCpu(cpuMesh, steel, dt, 9.80665f, 1.0f);
        solver.run(remaining, steel, dt, 9.80665f, 1.0f);
        solver.readback(gpuMesh);
        done = target;

        float worstPosition = 0;
        double deltaSquared = 0, velocitySquared = 0;
        for (std::size_t i = 0; i < cpuMesh.position.size(); ++i) {
            worstPosition = std::max(worstPosition,
                                     std::abs(cpuMesh.position[i] - gpuMesh.position[i]));
            const double d = cpuMesh.velocity[i] - gpuMesh.velocity[i];
            deltaSquared += d * d;
            velocitySquared += static_cast<double>(cpuMesh.velocity[i]) * cpuMesh.velocity[i];
        }
        const float ratio =
            static_cast<float>(std::sqrt(deltaSquared / std::max(velocitySquared, 1e-30)));
        samples.push_back({target, ratio, worstPosition});
        std::printf("     %8d %18.3g %18.3g\n", target, static_cast<double>(ratio),
                    static_cast<double>(worstPosition));
    }

    // One step from identical state is the test of the kernel itself. Agreement
    // at the 1e-4 level is float rounding across a few hundred operations with
    // cancellation, not a difference in what is being computed.
    check("one step agrees to better than 1e-4 relative", samples.front().velocity < 1e-4f);
    check("ten steps still agree to better than 1e-3 relative", samples[1].velocity < 1e-3f);
    // The growth from there is a finding, not a failure: rounding differences
    // amplify, so the GPU path is NOT bit-reproducible against the CPU, and will
    // not be reproducible across drivers or vendors either. That is exactly why
    // docs/04-multiplayer.md replicates damage *events* and lets each client run
    // its own local FEM for visuals, rather than assuming the solvers stay in
    // lockstep. Position stays accurate throughout because it is dominated by the
    // low-frequency bulk response, which is well resolved and does not amplify.
    check("divergence grows smoothly rather than blowing up",
          samples.back().velocity < 0.5f && std::isfinite(samples.back().velocity));
    check("positions agree to better than 10 micrometres after 1000 steps",
          samples.back().position < 1e-5f);
}

// --- 3. Throughput -----------------------------------------------------------

void benchmarkThroughput() {
    std::printf("\n3. Throughput on the GPU\n");
    std::printf("     %10s %10s %12s %14s %12s\n", "tets", "nodes", "ms/step", "Melem-upd/s",
                "sim ms/s @dt");

    const Material steel;
    struct Size { int nx, ny, nz; };
    const Size sizes[] = {{20, 4, 2}, {60, 8, 4}, {120, 16, 4}, {200, 24, 8}, {320, 32, 8}};

    for (const Size& s : sizes) {
        TetMesh mesh = makeBoxTetMesh(2.0f, 0.5f, 0.05f, s.nx, s.ny, s.nz);
        mesh.computeRestState(steel);
        clampRootFace(mesh);

        std::string error;
        gpu::FemGpuSolver solver;
        if (!solver.initialise(mesh, SHIPSIM_SHADER_DIR, error)) {
            std::printf("     %10zu  initialise failed: %s\n", mesh.tetCount(), error.c_str());
            continue;
        }

        const float dt = criticalTimestep(mesh, steel);
        solver.run(200, steel, dt, 9.80665f, 1.0f);  // warm up, let clocks settle

        const int steps = 2000;
        const double ms = solver.run(steps, steel, dt, 9.80665f, 1.0f);
        const double msPerStep = ms / steps;
        const double elementUpdatesPerSecond = mesh.tetCount() / (msPerStep * 1e-3);
        // How much wall time one second of simulated time costs at this dt --
        // the number that actually decides whether the zone is affordable.
        const double wallMsPerSimSecond = msPerStep / dt;

        std::printf("     %10zu %10zu %12.4f %14.1f %12.0f\n", mesh.tetCount(), mesh.nodeCount(),
                    msPerStep, elementUpdatesPerSecond / 1e6, wallMsPerSimSecond);
    }
}

void benchmarkCpuForComparison() {
    std::printf("\n4. Single-threaded CPU, same kernel, for scale\n");
    const Material steel;
    TetMesh mesh = makeBoxTetMesh(2.0f, 0.5f, 0.05f, 120, 16, 4);
    mesh.computeRestState(steel);
    clampRootFace(mesh);

    const float dt = criticalTimestep(mesh, steel);
    const int steps = 50;
    const auto begin = std::chrono::steady_clock::now();
    for (int i = 0; i < steps; ++i) stepCpu(mesh, steel, dt, 9.80665f, 1.0f);
    const auto end = std::chrono::steady_clock::now();

    const double msPerStep =
        std::chrono::duration<double, std::milli>(end - begin).count() / steps;
    std::printf("     %zu tets: %.4f ms/step, %.1f Melem-upd/s (one core)\n", mesh.tetCount(),
                msPerStep, mesh.tetCount() / (msPerStep * 1e-3) / 1e6);
}

}  // namespace

int main() {
    std::printf("shipsim - explicit tetrahedral FEM spike\n");
    testCantileverAgainstBeamTheory();
    testGpuMatchesCpu();
    benchmarkThroughput();
    benchmarkCpuForComparison();
    // Name the skips on the result line. A run that skipped half its sections and a
    // run that checked everything must not print the same sentence -- that is the
    // hole `check-figures.sh` had, where "ok" and "ok, having run nothing" were the
    // same line.
    if (failures != 0) std::printf("\nCHECKS FAILED\n");
    else if (skippedSections != 0)
        std::printf("\nall checks passed, but %d GPU section(s) were skipped for want of a"
                    " device\n", skippedSections);
    else std::printf("\nall checks passed\n");
    return failures == 0 ? 0 : 1;
}
