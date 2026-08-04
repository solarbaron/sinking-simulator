// SPDX-License-Identifier: MIT
//
// **Drive a real bow into the ferry's real plating, with elements.**
//
// `ram_view` is the Phase 3 milestone end to end, and the deformation in it comes
// from `indentation.hpp` -- rigid-plastic membrane stretching of one bay, closed
// form, microseconds. This is the same question answered by
// `engine/sim/zone.{hpp,cpp}`: solid-shell elements over the struck plating, J2
// plasticity with ductile failure, explicit dynamics, and the panels that tore
// handed to `breachesFromFailedPanels` unchanged.
//
// It lives here rather than in `tests/test_zone.cpp` because it costs core-minutes
// and the unit gate has to stay honest. The unit tests check the pieces at unit
// scale -- the mesher against closed forms, the solver against conservation, the
// answer against the membrane model on a bay whose idealisation both models can
// hold. This is the only thing that runs the whole chain at ship scale.
//
// **What it is for, beyond exercising the chain.** The two models disagree about
// something structural and the disagreement is worth more than either answer:
// `impactDamage()` takes the membrane span as the **frame** spacing, 2.4 m on this
// ship, and the plating between two longitudinals spans the **longitudinal**
// spacing, 0.70 m. Which one is right changes the energy a bay absorbs by an order
// of magnitude. The FEM has no span in it at all -- it has plating and boundaries
// -- so it is the instrument that settles it, and the number it produces is
// printed against both readings.
//
//   ./zone_probe [--speed=M_PER_S] [--depth=METRES] [--radius=METRES]
//                [--sub=N] [--aim=X_METRES] [--height=Z_METRES] [--threads=N]
//                [--elastic]
#include "engine/core/jobs.hpp"
#include "engine/sim/breach.hpp"
#include "engine/sim/indentation.hpp"
#include "engine/sim/scantlings.hpp"
#include "engine/sim/zone.hpp"
#include "game/prototype/ferry.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

namespace {

struct Options {
    double speed = 6.0;    // m/s, the punch's approach
    double depth = 0.45;   // m of penetration to drive to
    double radius = 4.0;   // m of plating to promote
    double aim = 0.0;      // m along the ship
    double height = 8.0;   // m above the baseline
    int subdivision = 4;
    int threads = 0;       // 0 takes the job system's default
    bool elastic = false;
};

bool parse(int argc, char** argv, Options& o) {
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        const auto value = [&](const char* key) -> const char* {
            const std::string prefix = std::string("--") + key + "=";
            return a.rfind(prefix, 0) == 0 ? a.c_str() + prefix.size() : nullptr;
        };
        if (const char* v = value("speed")) o.speed = std::atof(v);
        else if (const char* v = value("depth")) o.depth = std::atof(v);
        else if (const char* v = value("radius")) o.radius = std::atof(v);
        else if (const char* v = value("aim")) o.aim = std::atof(v);
        else if (const char* v = value("height")) o.height = std::atof(v);
        else if (const char* v = value("sub")) o.subdivision = std::atoi(v);
        else if (const char* v = value("threads")) o.threads = std::atoi(v);
        else if (a == "--elastic") o.elastic = true;
        else {
            std::printf("unknown option %s\n", a.c_str());
            return false;
        }
    }
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    Options options;
    if (!parse(argc, argv, options)) return 2;

    const sim::Ship ferry = game::buildFerry();
    const sim::Scantlings scantlings = sim::ferryScantlings();
    const sim::StructuralMesh structure = sim::makeStructuralMesh(ferry.hull, scantlings);
    std::printf("ferry  : %zu plating panels, %zu members, frames at %.2f m\n",
                structure.panels.size(), structure.members.size(), structure.frameSpacing);

    // --- 1. Promote a zone -----------------------------------------------------
    const sim::Vec3 impact{options.aim, -9.9, options.height};
    sim::zone::MeshParams mesh;
    mesh.radius = options.radius;
    mesh.subdivision = options.subdivision;
    const sim::zone::Patch patch = sim::zone::buildPatch(structure, impact, mesh);
    if (patch.empty()) {
        std::printf("no zone: nothing to promote at that impact point\n");
        for (const std::string& problem : patch.problems) std::printf("       ! %s\n",
                                                                     problem.c_str());
        return 1;
    }
    std::printf("zone   : %zu panels -> %zu elements, %zu nodes, %.1f m2 of %.0f mm plating,"
                " %.1f t\n",
                patch.panels.size(), patch.elementCount(), patch.nodeCount(), patch.area,
                patch.thickness * 1000.0, patch.mass / 1000.0);
    std::printf("         %d node(s) held by a stiffener line; %.0f%% of the zone's degrees of"
                " freedom are free\n", patch.stiffenerNodes, 100.0 * patch.freeFraction);
    std::printf("         faces %.4f rad from parallel on %d element(s) -> +%.0f%% bending"
                " stiffness at worst; aspect %.2f\n",
                patch.worstNormalSpread, patch.distortedElements,
                100.0 * patch.spuriousStiffness, patch.worstAspect);
    for (const std::string& problem : patch.problems)
        std::printf("       ! %s\n", problem.c_str());

    const double cost = sim::zone::estimatedCost(patch, !options.elastic);
    std::printf("cost   : dt = %.3f us (t/c_p, in-plane size irrelevant), %.0f steps per"
                " simulated second\n", patch.criticalTimestep * 1e6, 1.0 / patch.criticalTimestep);
    std::printf("         predicted %.0f core-seconds per simulated second; this run is"
                " %.3f s of it\n", cost, options.depth / options.speed);

    // --- 2. Drive a punch into it ----------------------------------------------
    //
    // A rigid flat punch two metres across, which is a bow shoulder rather than a
    // bulbous stem. The zone is clamped at its perimeter, so the answer stiffens
    // as the radius shrinks -- run it at two radii to see by how much.
    core::JobSystem jobs(options.threads > 0 ? static_cast<unsigned>(options.threads)
                                             : core::JobSystem::defaultWorkerCount());
    sim::zone::SolveParams solve;
    solve.plastic = !options.elastic;
    solve.jobs = &jobs;
    solve.indenter.halfLength = 1.0;   // m along the ship
    solve.indenter.halfWidth = 1.0;    // m up her side
    solve.indenter.speed = options.speed;
    solve.indenter.rampTime = 4.0e-3;
    solve.indenter.stopAt = options.depth;
    solve.historyStride = static_cast<int>(0.005 / patch.criticalTimestep);

    sim::zone::Solver solver(patch, sim::plasticity::shipSteel(), solve);
    const sim::zone::SolveResult& result = solver.run();

    // The membrane model on the same plating, at both readings of its span, sample
    // by sample -- because the interesting comparison is *before* anything tears
    // and the final numbers are all after.
    const auto membrane = [&](double span, double penetration, bool energy) {
        sim::IndentedPanel model;
        model.span = span;
        model.contactWidth = 2.0 * solve.indenter.halfLength;
        model.thickness = patch.thickness;
        model.yieldStrength = patch.material.yieldStrength;
        model.failureStrain = sim::plasticity::regularisedFailureStrain(
            sim::plasticity::shipSteel().failure, span, model.thickness);
        // The punch spans several bays of the shorter pitch, and each resists.
        const double bays = std::max(1.0, 2.0 * solve.indenter.halfWidth / span);
        return bays * (energy ? sim::indentationEnergy(model, penetration)
                              : sim::indentationForce(model, penetration));
    };
    const double shortSpan = scantlings.longitudinalSpacing;
    const double longSpan = structure.frameSpacing;

    std::printf("\n%9s %9s %12s %12s %12s %8s %12s %12s\n", "t (ms)", "depth", "force (MN)",
                "work (MJ)", "plastic (MJ)", "torn", "F: L=0.70", "F: L=2.40");
    for (const sim::zone::Sample& s : result.history)
        std::printf("%9.2f %9.3f %12.2f %12.3f %12.3f %8d %12.2f %12.2f\n", s.time * 1e3,
                    s.penetration, s.force / 1e6, s.work / 1e6, s.dissipation / 1e6,
                    s.tornElements, membrane(shortSpan, s.penetration, false) / 1e6,
                    membrane(longSpan, s.penetration, false) / 1e6);

    std::printf("\nsolve  : %d steps in %.2f s wall on %u workers, %.2f us/element/step\n",
                result.steps, result.wallSeconds,
                options.threads > 0 ? static_cast<unsigned>(options.threads)
                                    : core::JobSystem::defaultWorkerCount(),
                result.microsecondsPerElementStep);
    std::printf("energy : in %.3f MJ = strain %.4f + plastic %.3f + kinetic %.4f"
                "  (residual %+.2f%%)\n",
                result.work / 1e6, result.strainEnergy / 1e6, result.dissipation / 1e6,
                result.kinetic / 1e6,
                result.work > 0 ? 100.0 * result.energyResidual() / result.work : 0.0);
    for (const std::string& problem : result.problems)
        std::printf("       ! %s\n", problem.c_str());

    // --- 3. What tore ----------------------------------------------------------
    std::printf("damage : %d of %zu elements deleted (%.2f m2), %zu of %zu panels torn"
                " (%.2f m2)\n",
                result.tornElements, patch.elementCount(), result.tornArea,
                result.tornPanels.size(), patch.panels.size(), result.tornPanelArea);

    const sim::BreachSet breaches =
        sim::breachesFromFailedPanels(ferry, structure, result.tornPanels);
    double open = 0;
    for (const sim::Breach& b : breaches.breaches) open += b.opening.area;
    std::printf("breach : %zu opening(s), %.2f m2 reaching a compartment\n",
                breaches.breaches.size(), open);
    for (const std::string& problem : breaches.problems)
        std::printf("       ! %s\n", problem.c_str());

    // --- 4. Against the membrane model -----------------------------------------
    //
    // Both readings of the span, because they differ by an order of magnitude and
    // the FEM is the thing that can say which is right. Taken at the last sample
    // before anything tore, since the membrane model has nothing to say past that.
    double before = 0, forceBefore = 0, workBefore = 0;
    for (const sim::zone::Sample& s : result.history)
        if (s.tornElements == 0) {
            before = s.penetration;
            forceBefore = s.force;
            workBefore = s.work;
        }
    std::printf("\nat %.3f m, the last sample with nothing torn:\n", before);
    std::printf("%-38s %10s %12s %12s\n", "", "span (m)", "force (MN)", "energy (MJ)");
    std::printf("%-38s %10.2f %12.2f %12.3f\n", "membrane, span = longitudinal spacing",
                shortSpan, membrane(shortSpan, before, false) / 1e6,
                membrane(shortSpan, before, true) / 1e6);
    std::printf("%-38s %10.2f %12.2f %12.3f\n", "membrane, span = frame spacing", longSpan,
                membrane(longSpan, before, false) / 1e6, membrane(longSpan, before, true) / 1e6);
    std::printf("%-38s %10s %12.2f %12.3f\n", "solid-shell FEM", "-", forceBefore / 1e6,
                workBefore / 1e6);

    std::printf("\nok\n");
    return 0;
}
