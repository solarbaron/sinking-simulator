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
// **And it drives the same collision twice.** `zone.hpp` §6: a collision delivers
// joules, and until `Drive::Inertial` existed the zone consumed a prescribed
// travel, so the depth of the hole was an assumption. The tool now runs the punch
// to `--depth`, takes the energy that cost, and hands it straight back as a
// striking body carrying exactly those joules -- and prints where *that* stops,
// beside what the membrane model says the same joules would do. Three answers to
// one question in the units a collision poses it in.
//
//   ./zone_probe [--speed=M_PER_S] [--depth=METRES] [--radius=METRES]
//                [--sub=N] [--aim=X_METRES] [--height=Z_METRES] [--threads=N]
//                [--elastic] [--no-preload] [--no-energy] [--time-budget=X]
#include "engine/core/jobs.hpp"
#include "engine/sim/breach.hpp"
#include "engine/sim/indentation.hpp"
#include "engine/sim/promotion.hpp"
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
    // Solve the zone as if the hull girder were putting nothing through it, which
    // is what it did before `promotion.hpp` existed. Here so the pre-load's effect
    // at ship scale is a measurement rather than a claim.
    bool noPreload = false;
    // Keep the step-invariant element forms, or rebuild them every step. Here as a
    // switch because it is a 2x cost decision whose two answers must agree to the
    // last bit, and the only way to say that is to run both.
    std::string formsCache = "auto";
    // Mesh and solve at the requested point even when the promotion criterion
    // declines. The criterion holds an element budget, so it is the thing that
    // stops a cost sweep from reaching the sizes the budget exists to refuse --
    // which are exactly the sizes a cost sweep wants.
    bool force = false;
    // Skip the second, energy-driven solve. It roughly doubles the run, and a cost
    // sweep or a mesh-convergence study wants the prescribed drive alone -- two
    // meshes compared at one penetration are comparable, where two meshes compared
    // at one energy stop at two different penetrations.
    bool noEnergy = false;
    // How much longer than the prescribed run the energy-driven one may take. It is
    // a **cost** bound and not a travel one, and it is needed because the two are
    // not the same: a striker that has nearly stopped crawls, so the last
    // centimetres of an arresting run are the expensive ones. Measured on this very
    // patch, an inertial punch that perforated and then coasted at 0.5 m/s took
    // 216 000 steps against the prescribed run's 21 000.
    double timeBudget = 3.0;
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
        else if (const char* v = value("forms-cache")) o.formsCache = v;
        else if (const char* v = value("time-budget")) o.timeBudget = std::atof(v);
        else if (a == "--no-energy") o.noEnergy = true;
        else if (a == "--force") o.force = true;
        else if (a == "--elastic") o.elastic = true;
        else if (a == "--no-preload") o.noPreload = true;
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

    sim::Ship ferry = game::buildFerry();
    ferry.initialise(0.0);
    const sim::Scantlings scantlings = sim::ferryScantlings();
    // **The mesher's own account of what it could not build, which this tool used to
    // throw away.** `makeStructuralMesh` always returns *something* -- a description
    // whose frames will not lay out yields an empty mesh and the reason goes only
    // into `problems` -- so passing no pointer meshes a zone out of a structure
    // nobody checked. Everything below is a function of that structure: the patch is
    // cut from its panels, the pre-load is read off its section, and the reaction at
    // the end re-sections it damaged.
    std::vector<std::string> meshProblems;
    const sim::StructuralMesh structure =
        sim::makeStructuralMesh(ferry.hull, scantlings, &meshProblems);
    std::printf("ferry  : %zu plating panels, %zu members, frames at %.2f m\n",
                structure.panels.size(), structure.members.size(), structure.frameSpacing);
    for (const std::string& problem : meshProblems)
        std::printf("       ! %s\n", problem.c_str());

    // --- 0. Tier 0, and what it decides ----------------------------------------
    //
    // The whole point of `promotion.hpp`: the zone below is not promoted because a
    // command line asked for it, it is promoted because a load arrived that a beam
    // cannot represent. Both halves are exercised -- the ship on a wave, which
    // promotes nothing, and the same ship with a bow against her side, which does.
    const double length = ferry.hullHi.x - ferry.hullLo.x;
    const double omega = std::sqrt(sim::kGravity * 2.0 * sim::kPi / length);
    const sim::WaveField wave = sim::WaveField::regular(3.0, omega, 0.0);
    sim::Sea crest;
    crest.waves = &wave;
    {
        const double period = 2.0 * sim::kPi / omega;
        double best = -1e30;
        for (int i = 0; i < 720; ++i) {
            const double t = period * i / 720.0;
            const double eta = wave.elevation(ferry.state.position.x, 0.0, t);
            if (eta > best) { best = eta; crest.time = t; }
        }
    }
    const sim::promotion::TierZero tier =
        sim::promotion::tierZero(ferry, crest, structure, scantlings);
    std::printf("tier-0 : on a 3 m crest -- yield %.3f, buckling %.3f, collapse %.3f;"
                " %.0f ms to know\n", tier.yieldUtilisation, tier.buckleUtilisation,
                tier.collapseUtilisation, tier.seconds * 1e3);
    // **The three numbers above are the ones the promotion criterion fires on, and
    // `tier.problems` is the only place that says whether they mean anything.** Two
    // things arrive here and both are unsafe to drop. `validateGirder` puts the
    // shear and moment closure at the perpendiculars in it, and it only speaks above
    // 5% -- so a few percent of imbalance moves the peak bending moment with nothing
    // else lit anywhere, and the headline reads exactly as it does on a ship that
    // balanced. And `tierZero` records how many girder stations produced no strength
    // curve at all: those stations are absent from `collapseUtilisation` rather than
    // zero in it, so a station that went missing is one the collapse trigger cannot
    // fire on, and the error is in the unsafe direction.
    for (const std::string& problem : tier.problems)
        std::printf("       ! %s\n", problem.c_str());

    const sim::Vec3 impact{options.aim, -9.9, options.height};
    sim::promotion::Criterion criterion;
    criterion.mesh.radius = options.radius;
    criterion.mesh.subdivision = options.subdivision;

    sim::promotion::ContactPatch bow;
    bow.centre = impact;
    bow.radius = 1.5;
    // A 5175 t bow at the run's approach speed, stopped in a tenth of a second:
    // the order `ram_view` reports, and enough to be a load rather than a touch.
    bow.force = 5.175e6 * options.speed / 0.1;

    sim::promotion::Promoter promoter(criterion);
    for (int review = 0; review < 4; ++review) {
        const sim::promotion::Review r =
            promoter.review(structure, tier, review == 0 ? std::vector<sim::promotion::ContactPatch>{}
                                                         : std::vector{bow});
        std::printf("promote: review %d -- %zu candidate(s), %zu promoted, %d element(s) active,"
                    " %.0f core-s/s, decided in %.0f us\n", review, r.considered.size(),
                    r.promoted.size(), r.elementsActive, r.costActive, r.microseconds);
        for (const sim::promotion::Candidate& c : r.considered)
            std::printf("         panel %d, %s, score %.2f: %s\n", c.panel,
                        sim::promotion::name(c.trigger), c.score, c.why.c_str());
        // **A budget refusal and a criterion that did not fire are the same output
        // until this line prints.** `Promoter::review` drops a candidate the element
        // budget had no room for into `problems` and nowhere else: it is absent from
        // `promoted`, nothing already running was evicted for it, and the candidate
        // line above still carries the `why` that says it qualified. Without the
        // channel the review reads as hysteresis doing its job -- one number smaller
        // on the headline -- and the verdict a few lines down then states the wrong
        // cause outright, "the criterion did not fire" about a criterion that fired
        // and was overruled. It is printed inside the loop so it stands *before* that
        // verdict rather than after it. Same shape as the refusals `gasCandidates`
        // reports through its own `problems`, found there first.
        for (const std::string& problem : r.problems)
            std::printf("       ! %s\n", problem.c_str());
    }
    if (promoter.active().empty() && !options.force) {
        std::printf("nothing promoted: the criterion did not fire\n");
        return 1;
    }

    // --- 1. Mesh the zone the criterion asked for ------------------------------
    sim::zone::MeshParams mesh = criterion.mesh;
    sim::Vec3 meshAt = impact;
    if (!promoter.active().empty()) {
        const sim::promotion::Active& zone = promoter.active().front();
        mesh.role = zone.role;
        meshAt = zone.impact;
    } else {
        std::printf("forced : the criterion declined and --force meshed at the impact"
                    " point anyway\n");
    }
    const sim::zone::Patch patch = sim::zone::buildPatch(structure, meshAt, mesh);
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

    // **Every input that moves a printed number, on one line, so a table can carry
    // its own provenance.** `zone_gpu_probe` has always printed a `run :` line and
    // this tool never did, and the asymmetry cost two published tables their
    // parameters: §8's profile named `--radius=2.5` and omitted `--depth`, whose
    // default of 0.45 m is 42 396 steps and half the patch torn rather than the
    // 6 608 steps the table was taken at -- a different experiment wearing the same
    // invocation. The `RestForms` 17 800-element row lost its mesh entirely and had
    // to be found by search, and its published time still does not reproduce.
    //
    // The rule this encodes: a tool whose answer depends on a default must print
    // the default it used, because the person transcribing the number is copying
    // the command line and the command line is exactly where the default is not.
    std::printf("run    : speed %.3f m/s, depth %.3f m, radius %.3f m, aim %.3f m,"
                " height %.3f m, sub %d, forms-cache %s%s%s%s\n",
                options.speed, options.depth, options.radius, options.aim,
                options.height, options.subdivision, options.formsCache.c_str(),
                options.elastic ? ", elastic" : "", options.force ? ", force" : "",
                options.noPreload ? ", no-preload" : "");

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
    // The girder's own stress through this patch, which is what the plating is
    // already carrying before the bow arrives.
    const sim::promotion::PreloadCheck preload =
        sim::promotion::preloadFor(tier.girder, structure, patch);
    // **`applied` is the check's own verdict, not the solve's.** `preloadFor` says
    // whether the girder *has* a pre-stress worth applying; whether the solve below
    // receives it is `--no-preload`, decided ten lines down. The two used to be
    // printed as one field, so the control run announced `applied 1` on a solve that
    // got `Preload{}` -- a diagnostic naming a different quantity from the one
    // computed, which is a defect shape this repo has already been caught by. Both
    // are printed now, because both are true and they are not the same thing.
    std::printf("preload: M %.3e N m, neutral axis %.2f m -> %.1f MPa through the patch"
                " (%.0f%% of yield spent before contact); obliquity %.4f rad,"
                " available %d, given to the solve %d\n",
                preload.moment, preload.neutralAxis, preload.surfaceStress / 1e6,
                100.0 * std::abs(preload.surfaceStress) / patch.material.yieldStrength,
                preload.obliquity, static_cast<int>(preload.applied),
                static_cast<int>(preload.applied && !options.noPreload));
    for (const std::string& problem : preload.problems) std::printf("       ! %s\n",
                                                                    problem.c_str());

    sim::zone::SolveParams solve;
    solve.plastic = !options.elastic;
    solve.cacheRestForms = options.formsCache != "never";
    solve.preload = options.noPreload ? sim::zone::Preload{} : preload.preload;
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

    std::printf("\nsolve  : %d steps in %.2f s wall on %u workers, %.2f us/element/step,"
                " rest forms %s\n",
                result.steps, result.wallSeconds,
                options.threads > 0 ? static_cast<unsigned>(options.threads)
                                    : core::JobSystem::defaultWorkerCount(),
                result.microsecondsPerElementStep,
                result.cachedRestForms ? "cached" : "rebuilt each step");
    // A bit-exact digest of the answer. Two solver configurations that claim to be
    // the same arithmetic have to produce the same bits, and a printed force to two
    // decimals cannot say so -- the hex is what makes the claim checkable from a
    // shell. `%a` is exact for a double.
    std::printf("digest : force %a peak %a work %a strain %a plastic %a kinetic %a"
                " torn %d steps %d\n",
                result.force, result.peakForce, result.work, result.strainEnergy,
                result.dissipation, result.kinetic, result.tornElements, result.steps);

    // Where the time went. **This is the number that decides whether an
    // accelerator is worth building**, because Amdahl's law bounds what moving one
    // phase can buy however fast the replacement is. The element kernel's share is
    // the ceiling on a GPU element solver that leaves the rest on the CPU.
    {
        const sim::zone::SolveResult::Profile& p = result.profile;
        const double total = p.total();
        const auto share = [&](double s) { return total > 0 ? 100.0 * s / total : 0.0; };
        std::printf("profile: %.2f s accounted of %.2f s wall\n", total, result.wallSeconds);
        std::printf("         element   %8.3f s  %5.1f%%   (per-element force / return map)\n",
                    p.element, share(p.element));
        std::printf("         gather    %8.3f s  %5.1f%%   (CSR nodal gather + per-element"
                    " reduction)\n", p.gather, share(p.gather));
        std::printf("         integrate %8.3f s  %5.1f%%   (velocity, punch, position)\n",
                    p.integrate, share(p.integrate));
        std::printf("         energy    %8.3f s  %5.1f%%   (strain + kinetic energy"
                    " accumulation)\n", p.energy, share(p.energy));
        std::printf("         other     %8.3f s  %5.1f%%\n", p.other, share(p.other));
    }

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

    // --- 3b. And back to Tier 0 -------------------------------------------------
    //
    // The other half of the coupling: what the zone lost, as a section the beam
    // already knows how to read. Nothing in `girder.hpp`, `buckling.hpp` or
    // `collapse.hpp` is reimplemented here -- the damaged ship is a thinner ship.
    const sim::promotion::SectionReduction reduction =
        sim::promotion::reactionOf(structure, patch, solver);
    std::printf("react  : %zu panel(s) reduced, worst effectiveness %.3f, worst dent %.3f m,"
                " %.2f m2 and %.0f kg of steel no longer carrying\n", reduction.panels.size(),
                reduction.worstEffectiveness, reduction.worstOutOfPlane, reduction.lostPlateArea,
                reduction.lostSteelMass);
    for (const std::string& problem : reduction.problems) std::printf("       ! %s\n",
                                                                      problem.c_str());
    if (!reduction.empty()) {
        const sim::StructuralMesh damaged = sim::promotion::reduce(structure, reduction);
        const double station = patch.centre.x;
        const sim::HullGirderSection before = sim::hullGirderSection(structure, station);
        const sim::HullGirderSection after = sim::hullGirderSection(damaged, station);
        // **Both curves swept over the same range, at the same step count.**
        // `collapseCurve` sizes its own sweep from its own first-yield curvature,
        // and a damaged section has a different first yield -- so asking it twice
        // and comparing the peaks compares two maxima found at *different*
        // resolutions. The peak is picked off the samples, so that difference is a
        // sampling floor of order a tenth of a percent.
        //
        // It matters here because this damage barely moves the section: the panels
        // a contact tears out sit near the neutral axis, so `I` falls 0.08% and the
        // ultimate moment with it. The ordering was being decided below the
        // sampling floor, and it read -0.02% -- correct by a hundredth of a
        // percent, which is not the same as correct.
        //
        // Sharing the reach makes the sampling error common-mode, and 600 steps
        // rather than 150 puts the floor well under the signal.
        const double sign = tier.girder.hogging() ? 1.0 : -1.0;
        const std::vector<sim::CollapseElement> intactElements =
            sim::collapseElementsAt(structure, scantlings, station);
        const double reach = 6.0 * sim::firstYieldCurvature(intactElements);
        const auto ultimate = [&](const sim::StructuralMesh& m) {
            return sim::progressiveCollapse(sim::collapseElementsAt(m, scantlings, station),
                                            sign * reach, 600)
                .ultimateMoment;
        };
        const double intact = ultimate(structure), hurt = ultimate(damaged);
        std::printf("         at x = %.2f: area %.4f -> %.4f m2, I %.3f -> %.3f m4,"
                    " modulus %.4f -> %.4f m3\n", station, before.area, after.area,
                    before.secondMoment, after.secondMoment, before.modulusDeck, after.modulusDeck);
        std::printf("         ultimate moment %.4e -> %.4e N m (%+.2f%%), against an applied"
                    " %.4e; margin %.2f -> %.2f\n", intact, hurt,
                    100.0 * (hurt / intact - 1.0), tier.girder.maxMoment,
                    std::abs(intact / tier.girder.maxMoment),
                    std::abs(hurt / tier.girder.maxMoment));
        // **A strict ordering is not a property this method guarantees**, and
        // asserting one here was reading a sign off a difference of hundredths of a
        // percent. The panels a contact tears out sit near the neutral axis: `I`
        // falls 0.08% and the section modulus 0.26%, so the ultimate moment is
        // essentially unmoved, and which way it moves is decided by how the
        // instantaneous neutral axis rebalances once shedding compression elements
        // are removed -- which a Smith sweep can genuinely resolve either way.
        //
        // Measured across punch depths 0.16 to 0.25 m, +0.05% to +0.08%, rising
        // monotonically with the damage. Under the previous scantling error -- the
        // column check using the transverse frame's profile, so stiffeners could
        // not buckle at all -- it read a flat -0.02% at every depth, which is to say
        // the section barely responded to damage in either direction. That flatness
        // was the artefact; this sensitivity is real and its sign near the neutral
        // axis is a known limit of the method rather than a defect in it.
        //
        // So what is checked is that damage does not make her *materially* stronger.
        // Half a per cent is six times the largest reading here and far below
        // anything that would change a decision, and a section whose modulus has
        // genuinely collapsed cannot hide under it.
        if (hurt > intact * 1.005) {
            std::printf("       ! the damaged section is materially stronger than the intact one\n");
            return 1;
        }
    }

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

    // --- 5. The same collision, delivered as energy ------------------------------
    //
    // Everything above took a travel and reported the joules. A collision arrives
    // the other way round, so the joules go back in: a striking body whose kinetic
    // energy is exactly what the run above absorbed, released at the same approach
    // speed, and stopped by the plating rather than by a number in a struct. See
    // `zone.hpp` §6 for why the entry point is a mass and a velocity and not the
    // energy itself.
    //
    // The membrane model is asked the same question in the same units --
    // `penetrationForEnergy` is the closed-form inverse of its energy integral --
    // so the three answers are comparable without either model being run twice.
    if (!options.noEnergy && result.work > 0) {
        const double energy = result.work;
        const double mass = 2.0 * energy / (options.speed * options.speed);
        sim::zone::SolveParams spending = solve;
        spending.indenter.drive = sim::zone::Drive::Inertial;
        spending.indenter.mass = mass;
        spending.indenter.rampTime = 0.0;   // a body with a mass arrives travelling
        // The travel cap stays where the prescribed run stopped, so the two answers
        // are the same question; the *time* budget is the one that matters for cost,
        // because a striker that has nearly stopped crawls and the last centimetres
        // are the expensive ones.
        spending.duration = options.timeBudget * options.depth / options.speed;
        spending.maxSteps = 8000000;

        sim::zone::Solver striker(patch, sim::plasticity::shipSteel(), spending);
        const sim::zone::SolveResult& spent = striker.run();

        std::printf("\nenergy : the %.3f MJ that travel cost, handed back as a %.0f t striker at"
                    " %.1f m/s\n", energy / 1e6, mass / 1000.0, options.speed);
        std::printf("         %s at %.4f m against the prescribed %.4f m; %.4f MJ of %.4f left"
                    " unspent\n",
                    spent.indenterArrested ? "stopped" : "still moving",
                    spent.penetration, result.penetration, spent.indenterKinetic / 1e6,
                    spent.indenterEnergy / 1e6);
        std::printf("         %d element(s) torn against %d, %zu panel(s) against %zu;"
                    " %d steps against %d, %.2f s wall\n",
                    spent.tornElements, result.tornElements, spent.tornPanels.size(),
                    result.tornPanels.size(), spent.steps, result.steps, spent.wallSeconds);
        // The striker's own ledger. It closes far tighter than the patch's does --
        // it involves only the punch's kinematics and the force it was handed, where
        // the patch's runs through the whole constitutive path -- so the two are
        // printed together rather than one standing for the other.
        std::printf("         striker ledger %+.4g J on %.4g (%.4f%%); patch account"
                    " %+.2f%%\n", spent.indenterResidual(), spent.work,
                    spent.work > 0 ? 100.0 * spent.indenterResidual() / spent.work : 0.0,
                    spent.work > 0 ? 100.0 * spent.energyResidual() / spent.work : 0.0);
        for (const std::string& problem : spent.problems)
            std::printf("       ! %s\n", problem.c_str());

        // And the membrane model on the same joules, at both readings of the span.
        // `membrane(span, d, true)` is the energy at a depth; this is its inverse,
        // taken through `indentation.hpp`'s own closed form rather than by search.
        //
        // **Whether it tore is printed beside the depth, because otherwise the
        // number is not a penetration at all.** `penetrationForEnergy` returns the
        // *tearing* penetration once the energy exceeds what a bay can absorb
        // intact, so a torn row is a bay that let go and not a bay that stopped the
        // strike -- the two read identically in a column of metres.
        std::printf("\n%.3f MJ into her side, three ways:\n", energy / 1e6);
        std::printf("%-38s %10s %14s %10s\n", "", "span (m)", "depth (m)", "tore");
        const char* label[2] = {"membrane, span = longitudinal spacing",
                                "membrane, span = frame spacing"};
        const double spans[2] = {shortSpan, longSpan};
        for (int which = 0; which < 2; ++which) {
            const double span = spans[which];
            sim::IndentedPanel model;
            model.span = span;
            model.contactWidth = 2.0 * solve.indenter.halfLength;
            model.thickness = patch.thickness;
            model.yieldStrength = patch.material.yieldStrength;
            model.failureStrain = sim::plasticity::regularisedFailureStrain(
                sim::plasticity::shipSteel().failure, span, model.thickness);
            const double bays = std::max(1.0, 2.0 * solve.indenter.halfWidth / span);
            const double perBay = energy / bays;
            const bool tore = perBay >= sim::energyToTear(model);
            std::printf("%-38s %10.2f %14.4f %10s\n", label[which], span,
                        sim::penetrationForEnergy(model, perBay), tore ? "yes, capped" : "no");
        }
        std::printf("%-38s %10s %14.4f %10s%s\n", "solid-shell FEM, energy-driven", "-",
                    spent.penetration,
                    spent.tornElements > 0 ? "yes" : "no",
                    spent.indenterArrested ? "" : "   (depth is a cap, not the answer)");
        // The bias, stated at the point of use rather than left in a header: this
        // striker is rigid, so all of it went into her plating. A real bow crushes
        // and takes a share, which makes every depth above an upper bound.
        std::printf("the striking bow is rigid here, so every joule tore her plating: an upper"
                    " bound on all three\n");
    }

    std::printf("\nok\n");
    return 0;
}
