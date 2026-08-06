// SPDX-License-Identifier: MIT
//
// **The Tier-2 zone on the GPU, against the CPU double reference on the same
// patch.** Two questions, and neither can be answered by the kernel alone:
//
//   1. **Is it faster end to end?** Not "how many element updates per second" --
//      `07-fem-spike-findings.md` §3 already answers that for a tet and it is the
//      wrong question. What a promoted zone costs is wall time for a whole run, so
//      that is what is compared: the same patch, the same number of steps, the
//      same punch, one on 23 threads of double and one on a 1070 Ti of float.
//   2. **Does float suffice for a co-rotational solid-shell?** The CPU is double
//      and Pascal runs fp64 at 1/32 rate, so a double kernel is pointless -- but
//      "float is fine" is a claim and this measures it. The two places to expect
//      trouble are the polar decomposition, which is a fixed-point iteration whose
//      CPU tolerance (1e-16) float cannot represent, and the near-incompressible
//      bending response, where the enhanced thickness modes are scaled by 1/t^2
//      and their Newton is conditioned on the square of that.
//
// The comparison is on **integral quantities**: work in, plastic dissipation,
// penetration, the displacement field and the accumulated plastic strain.
// `07-fem-spike-findings.md` §2 records why -- comparing two explicit solvers on a
// peak nodal velocity reads over 200% difference and means nothing, because an
// explicit scheme has no accuracy at all in the modes near its stability limit.
//
//   ./zone_gpu_probe [--radius=M] [--sub=N] [--depth=M] [--speed=M_PER_S]
//                    [--aim=M] [--height=M] [--threads=N] [--steps=N]
//
// **Skips rather than fails with no Vulkan device.**
#include "engine/core/jobs.hpp"
#include "engine/gpu/zone_gpu.hpp"
#include "engine/sim/promotion.hpp"
#include "engine/sim/scantlings.hpp"
#include "engine/sim/zone.hpp"
#include "game/prototype/ferry.hpp"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <string>
#include <vector>

namespace {

struct Options {
    double radius = 2.5;
    double depth = 0.06;
    double speed = 6.0;
    double aim = 0.0;
    double height = 8.0;
    int subdivision = 4;
    int threads = 0;
    int steps = 0;  // 0 derives the count from `depth`
    // The negative control. Perturbs the patch geometry by this many metres and
    // runs the **double** CPU solver on it, so the divergence a float kernel
    // produces can be compared against the divergence the same size of perturbation
    // produces in double. Without it there is no way to tell a wrong kernel from a
    // problem that amplifies any perturbation at all.
    double jitter = 0.0;
    // Which device mapping to measure. `07-fem-spike-findings.md` §8's table was
    // taken on `invocation` -- one thread per element -- and attributed its shape to
    // register spilling rather than to the element, so the two have to be comparable
    // on one command line or the re-measurement is not a re-measurement.
    gpu::Mapping mapping = gpu::Mapping::Workgroup;
    bool stats = false;   // report the compiler's register and spill counts, then stop
};

bool parse(int argc, char** argv, Options& o) {
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        const auto value = [&](const char* key) -> const char* {
            const std::string prefix = std::string("--") + key + "=";
            return a.rfind(prefix, 0) == 0 ? a.c_str() + prefix.size() : nullptr;
        };
        if (const char* v = value("radius")) o.radius = std::atof(v);
        else if (const char* v = value("depth")) o.depth = std::atof(v);
        else if (const char* v = value("speed")) o.speed = std::atof(v);
        else if (const char* v = value("aim")) o.aim = std::atof(v);
        else if (const char* v = value("height")) o.height = std::atof(v);
        else if (const char* v = value("sub")) o.subdivision = std::atoi(v);
        else if (const char* v = value("threads")) o.threads = std::atoi(v);
        else if (const char* v = value("steps")) o.steps = std::atoi(v);
        else if (a == "--stats") o.stats = true;
        else if (const char* v = value("jitter")) o.jitter = std::atof(v);
        else if (const char* v = value("mapping")) {
            const std::string m = v;
            if (m == "workgroup") o.mapping = gpu::Mapping::Workgroup;
            else if (m == "invocation") o.mapping = gpu::Mapping::Invocation;
            else { std::printf("--mapping must be workgroup or invocation\n"); return false; }
        }
        else {
            std::printf("unknown option %s\n", a.c_str());
            return false;
        }
    }
    return true;
}

double relative(double got, double want) {
    const double scale = std::max(std::abs(want), 1e-30);
    return std::abs(got - want) / scale;
}

}  // namespace

int main(int argc, char** argv) {
    Options options;
    if (!parse(argc, argv, options)) return 2;

    // What the driver's compiler did with each mapping, before anything is timed.
    // §8 diagnosed register spilling from the *shape* of the throughput curve; this
    // asks the compiler directly, so the diagnosis is checked rather than inherited.
    if (options.stats) {
        std::printf("pipeline statistics, %s:\n", SHIPSIM_SHADER_DIR);
        std::printf("%s", gpu::describeElementPipelines(SHIPSIM_SHADER_DIR).c_str());
        return 0;
    }

    sim::Ship ferry = game::buildFerry();
    ferry.initialise(0.0);
    const sim::Scantlings scantlings = sim::ferryScantlings();
    const sim::StructuralMesh structure = sim::makeStructuralMesh(ferry.hull, scantlings);

    sim::zone::MeshParams mesh;
    mesh.radius = options.radius;
    mesh.subdivision = options.subdivision;
    const sim::Vec3 impact{options.aim, -9.9, options.height};
    const sim::zone::Patch patch = sim::zone::buildPatch(structure, impact, mesh);
    if (patch.empty()) {
        std::printf("no zone at that impact point\n");
        return 1;
    }
    std::printf("zone   : %zu panels -> %zu elements, %zu nodes, %.1f m2 of %.0f mm plating\n",
                patch.panels.size(), patch.elementCount(), patch.nodeCount(), patch.area,
                patch.thickness * 1000.0);

    sim::zone::SolveParams solve;
    solve.plastic = true;
    solve.indenter.halfLength = 1.0;
    solve.indenter.halfWidth = 1.0;
    solve.indenter.speed = options.speed;
    solve.indenter.rampTime = 4.0e-3;
    solve.indenter.stopAt = options.depth;

    const double timestep = patch.criticalTimestep * (solve.timestepSafety / 0.9);
    const int steps = options.steps > 0
                          ? options.steps
                          : static_cast<int>(options.depth / (options.speed * timestep));
    std::printf("run    : dt %.4f us, %d steps, punch %.1f m/s to %.3f m\n", timestep * 1e6, steps,
                options.speed, options.depth);

    // --- The GPU, first, so a missing device skips before anything expensive ------
    //
    // The device solver is handed the CPU solver's own lumped mass, CSR adjacency,
    // driven-node list, pinned degrees of freedom and rest forms. Nothing is
    // rebuilt: a second mass lumping would be a difference no force could show.
    sim::zone::Solver seed(patch, sim::plasticity::shipSteel(), solve);
    gpu::ZoneGpuSolver device;
    std::string error;
    if (!device.initialise(patch, seed, sim::plasticity::shipSteel(), solve, SHIPSIM_SHADER_DIR,
                           error, options.mapping)) {
        std::printf("skipped: %s\n", error.c_str());
        return 0;
    }
    std::printf("device : %s, %s mapping\n", device.deviceName().c_str(),
                options.mapping == gpu::Mapping::Workgroup ? "one workgroup per element"
                                                           : "one invocation per element");

    const auto gpuBegin = std::chrono::steady_clock::now();
    const double kernelMs = device.run(steps);
    const gpu::ZoneGpuState state = device.readback();
    const auto gpuEnd = std::chrono::steady_clock::now();
    const double gpuWall = std::chrono::duration<double>(gpuEnd - gpuBegin).count();

    // Continue through the CPU's own energy account, tearing rules and panel
    // reporting, so the two answers are read by the same code.
    sim::zone::Solver adopted(patch, sim::plasticity::shipSteel(), solve);
    adopted.adopt(state.position, state.velocity, state.plastic, state.steps, state.time,
                  state.penetration, state.work, state.dissipation);
    const sim::zone::SolveResult& gpuResult = adopted.result();

    // --- The CPU reference --------------------------------------------------------
    core::JobSystem jobs(options.threads > 0 ? static_cast<unsigned>(options.threads)
                                             : core::JobSystem::defaultWorkerCount());
    sim::zone::SolveParams cpuSolve = solve;
    cpuSolve.jobs = &jobs;
    cpuSolve.maxSteps = steps;
    cpuSolve.indenter.stopAt = 0.0;      // run to the step count, not to a depth
    cpuSolve.duration = steps * timestep;
    sim::zone::Solver cpu(patch, sim::plasticity::shipSteel(), cpuSolve);
    const sim::zone::SolveResult& cpuResult = cpu.run();

    const unsigned workers = options.threads > 0 ? static_cast<unsigned>(options.threads)
                                                 : core::JobSystem::defaultWorkerCount();
    std::printf("\n%-28s %16s %16s %12s\n", "", "CPU (double)", "GPU (float)", "relative");
    std::printf("%-28s %16d %16d\n", "steps", cpuResult.steps, gpuResult.steps);
    std::printf("%-28s %16.4f %16.4f %12s\n", "wall seconds", cpuResult.wallSeconds, gpuWall, "");
    std::printf("%-28s %16.4f %16.4f %12s\n", "  of which kernel (ms)",
                cpuResult.profile.element * 1e3, kernelMs, "");
    std::printf("%-28s %16.6f %16.6f %12.2e\n", "penetration (m)", cpuResult.penetration,
                gpuResult.penetration, relative(gpuResult.penetration, cpuResult.penetration));
    std::printf("%-28s %16.4f %16.4f %12.2e\n", "work in (MJ)", cpuResult.work / 1e6,
                gpuResult.work / 1e6, relative(gpuResult.work, cpuResult.work));
    std::printf("%-28s %16.4f %16.4f %12.2e\n", "plastic dissipation (MJ)",
                cpuResult.dissipation / 1e6, gpuResult.dissipation / 1e6,
                relative(gpuResult.dissipation, cpuResult.dissipation));
    std::printf("%-28s %16.4f %16.4f %12.2e\n", "strain energy (MJ)", cpuResult.strainEnergy / 1e6,
                gpuResult.strainEnergy / 1e6,
                relative(gpuResult.strainEnergy, cpuResult.strainEnergy));
    std::printf("%-28s %16d %16d\n", "torn elements", cpuResult.tornElements,
                gpuResult.tornElements);

    // The displacement field, which is what the energies are integrals of. RMS and
    // max over every node, against the travel the patch actually made -- a relative
    // measure, because an absolute one says nothing until the patch has moved.
    const std::vector<double>& a = cpu.position();
    const std::vector<double>& b = adopted.position();
    double sum = 0, worst = 0, travel = 0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        const double d = b[i] - a[i];
        sum += d * d;
        worst = std::max(worst, std::abs(d));
        travel = std::max(travel, std::abs(a[i] - patch.mesh.position[i]));
    }
    const double rms = a.empty() ? 0.0 : std::sqrt(sum / static_cast<double>(a.size()));
    std::printf("%-28s %16.3e %16.3e %12.2e\n", "node position RMS / max (m)", rms, worst,
                travel > 0 ? worst / travel : 0.0);
    std::printf("%-28s %16.4f\n", "  patch travel (m)", travel);

    // The plastic strain, which is the quantity the tearing criterion reads and the
    // one a float return map is most likely to drift on.
    double worstStrain = 0, peakStrain = 0;
    for (std::size_t e = 0; e < patch.elementCount(); ++e)
        for (int gp = 0; gp < sim::solidshell::kGauss; ++gp) {
            const double p = cpu.elementState()[e].point[gp].equivalentPlasticStrain;
            const double q = state.plastic[e].point[gp].equivalentPlasticStrain;
            peakStrain = std::max(peakStrain, p);
            worstStrain = std::max(worstStrain, std::abs(p - q));
        }
    std::printf("%-28s %16.4f %16.3e %12.2e\n", "peak eps_p / worst diff", peakStrain, worstStrain,
                peakStrain > 0 ? worstStrain / peakStrain : 0.0);

    // **How close either path came to tearing at all**, which is what decides
    // whether a torn-element count is a measurement or a coincidence. Damage is the
    // integral the tearing criterion reads: a point fails at 1.0 and an element is
    // torn when all eight of its points have. Reporting the peak on both sides turns
    // "0 torn against 0 torn" from an agreement that might be two zeros for
    // unrelated reasons into a margin that can be quoted.
    double peakDamage = 0, peakDamageGpu = 0;
    int failedPoints = 0, failedPointsGpu = 0;
    for (std::size_t e = 0; e < patch.elementCount(); ++e)
        for (int gp = 0; gp < sim::solidshell::kGauss; ++gp) {
            const auto& p = cpu.elementState()[e].point[gp];
            const auto& q = state.plastic[e].point[gp];
            peakDamage = std::max(peakDamage, p.damage);
            peakDamageGpu = std::max(peakDamageGpu, q.damage);
            if (p.failed) ++failedPoints;
            if (q.failed) ++failedPointsGpu;
        }
    std::printf("%-28s %16.4f %16.4f %12.2e\n", "peak damage (1.0 fails)", peakDamage,
                peakDamageGpu, relative(peakDamageGpu, peakDamage));
    std::printf("%-28s %16d %16d\n", "failed Gauss points", failedPoints, failedPointsGpu);

    // The seven enhanced parameters, which are the per-element internal variables
    // the GPU condenses out every step. They are the part of this element most
    // likely to suffer in float: the enhanced thickness modes carry E_zeta,zeta, so
    // their columns of G are scaled by the Voigt transform's 1/t^2 -- of order 7e3
    // for 12 mm plate -- and Kaa, which the 7x7 Cholesky factors, inherits the
    // square of that. Double has sixteen digits to spend on it and float has seven.
    double worstAlpha = 0, peakAlpha = 0, peakAlphaGpu = 0;
    for (std::size_t e = 0; e < patch.elementCount(); ++e)
        for (int k = 0; k < sim::solidshell::kEas; ++k) {
            const double p = cpu.elementState()[e].enhanced[k];
            const double q = state.plastic[e].enhanced[k];
            peakAlpha = std::max(peakAlpha, std::abs(p));
            peakAlphaGpu = std::max(peakAlphaGpu, std::abs(q));
            worstAlpha = std::max(worstAlpha, std::abs(p - q));
        }
    std::printf("%-28s %16.3e %16.3e %12.2e\n", "peak |alpha|", peakAlpha, peakAlphaGpu,
                peakAlpha > 0 ? worstAlpha / peakAlpha : 0.0);

    // --- The negative control ----------------------------------------------------
    //
    // **Is the float kernel wrong, or is the problem this sensitive?** The only way
    // to tell is to perturb the double solver by the same amount float rounding
    // perturbs it and see how far *it* moves. A geometric jitter of the size of the
    // float mesh resolution is the right perturbation: that is exactly what
    // narrowing the node positions to float does to the element geometry.
    if (options.jitter > 0) {
        sim::zone::Patch shaken = patch;
        // Deterministic, and not a random-number generator: an integer hash of the
        // coordinate index, so the same run always produces the same perturbation
        // and the comparison is repeatable.
        for (std::size_t i = 0; i < shaken.mesh.position.size(); ++i) {
            std::uint32_t h = static_cast<std::uint32_t>(i) * 2654435761u + 1013904223u;
            h ^= h >> 15;
            const double unit = static_cast<double>(h % 2000001u) / 1000000.0 - 1.0;
            shaken.mesh.position[i] += options.jitter * unit;
        }
        sim::zone::SolveParams shakenSolve = cpuSolve;
        sim::zone::Solver shakenCpu(shaken, sim::plasticity::shipSteel(), shakenSolve);
        const sim::zone::SolveResult& shakenResult = shakenCpu.run();

        double shakenSum = 0, shakenWorst = 0, shakenStrain = 0;
        for (std::size_t i = 0; i < a.size(); ++i) {
            const double d = shakenCpu.position()[i] - a[i];
            shakenSum += d * d;
            shakenWorst = std::max(shakenWorst, std::abs(d));
        }
        for (std::size_t e = 0; e < patch.elementCount(); ++e)
            for (int gp = 0; gp < sim::solidshell::kGauss; ++gp)
                shakenStrain = std::max(
                    shakenStrain,
                    std::abs(cpu.elementState()[e].point[gp].equivalentPlasticStrain -
                             shakenCpu.elementState()[e].point[gp].equivalentPlasticStrain));
        const double shakenRms =
            a.empty() ? 0.0 : std::sqrt(shakenSum / static_cast<double>(a.size()));
        std::printf("\ncontrol: the same double solver on a mesh jittered by %.1e m\n",
                    options.jitter);
        std::printf("%-28s %16.4f %16.4f %12.2e\n", "  plastic dissipation (MJ)",
                    cpuResult.dissipation / 1e6, shakenResult.dissipation / 1e6,
                    relative(shakenResult.dissipation, cpuResult.dissipation));
        std::printf("%-28s %16d %16d\n", "  torn elements", cpuResult.tornElements,
                    shakenResult.tornElements);
        std::printf("%-28s %16.3e %16.3e %12.2e\n", "  node RMS / max (m)", shakenRms,
                    shakenWorst, travel > 0 ? shakenWorst / travel : 0.0);
        std::printf("%-28s %16.4f %16.3e %12.2e\n", "  peak eps_p / worst diff", peakStrain,
                    shakenStrain, peakStrain > 0 ? shakenStrain / peakStrain : 0.0);
    }

    std::printf("\nspeedup: %.2fx end to end against %u CPU workers"
                " (%.2fx against the element kernel alone)\n",
                gpuWall > 0 ? cpuResult.wallSeconds / gpuWall : 0.0, workers,
                kernelMs > 0 ? cpuResult.profile.element * 1e3 / kernelMs : 0.0);

    // A guard against a vacuous comparison: two solvers that both did nothing agree
    // perfectly. The run has to have deformed the patch and yielded it.
    if (!(travel > 1e-4) || !(peakStrain > 1e-4)) {
        std::printf("       ! the run neither deformed nor yielded the patch, so the"
                    " comparison compared nothing\n");
        return 1;
    }
    std::printf("\nok\n");
    return 0;
}
