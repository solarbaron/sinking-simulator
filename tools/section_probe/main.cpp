// SPDX-License-Identifier: MIT
//
// **Cut the reference ferry into a hold and reduce it.**
//
// `tests/test_section.cpp` checks the mesher against closed forms at unit scale and
// against `hullGirderSection` on two frame bays, and it has to stay inside the unit
// gate. This is the thing that runs at ship scale: eleven bays of the ferry, the
// resolution sweep that says what a Tier-1 section actually needs, and the
// Craig-Bampton reduction of the result.
//
// **Every *ship-scale* figure `docs/02-simulation.md` §3 publishes about the section
// mesher comes out of this program**, so re-running it is how a claim there is
// checked rather than quoted -- which is the failure mode this repo has hit three
// times, most recently with a factor of ten that three documents repeated from each
// other while the code beneath them had been right for weeks. The box-girder figures
// there -- the cost of welding a corner, the Bredt ratios, what the thickness taper
// costs -- come from `tests/test_section.cpp`, where they are assertions rather than
// output.
//
//   ./section_probe [--from=X] [--to=X] [--sub=N] [--sweep=N] [--reduce] [--modes=N]
#include "engine/sim/girder.hpp"
#include "engine/sim/reduction.hpp"
#include "engine/sim/scantlings.hpp"
#include "engine/sim/section.hpp"
#include "game/prototype/ferry.hpp"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <numbers>
#include <string>
#include <vector>

namespace {

double now() {
    using namespace std::chrono;
    return duration<double>(steady_clock::now().time_since_epoch()).count();
}

struct Options {
    // Eleven frame bays about midship. Frame stations are multiples of 2.4 m; the
    // watertight bulkheads are at -44, -38, -8, 20, 44 and none of them is, so
    // cutting "at a bulkhead" would put the plane through 188 panels. See
    // `section.hpp` and the `--from=-8 --to=20` run, which reports exactly that.
    double from = -7.2, to = 19.2;
    int subdivision = 1;
    int sweep = 3;       // refine to this subdivision; 0 skips the sweep
    bool reduce = true;  // build the substructure and reduce it
    int modes = 0;       // -1 takes ReduceParams' frequency cutoff
};

bool parse(int argc, char** argv, Options& o) {
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        const auto value = [&](const char* key) -> const char* {
            const std::string prefix = std::string("--") + key + "=";
            return a.rfind(prefix, 0) == 0 ? a.c_str() + prefix.size() : nullptr;
        };
        if (const char* v = value("from")) o.from = std::atof(v);
        else if (const char* v = value("to")) o.to = std::atof(v);
        else if (const char* v = value("sub")) o.subdivision = std::atoi(v);
        else if (const char* v = value("sweep")) o.sweep = std::atoi(v);
        else if (const char* v = value("modes")) o.modes = std::atoi(v);
        else if (a == "--no-reduce") o.reduce = false;
        else if (a == "--reduce") o.reduce = true;
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
    std::vector<std::string> problems;
    const sim::StructuralMesh structure =
        sim::makeStructuralMesh(ferry.hull, sim::ferryScantlings(), &problems);
    const sim::StructuralMaterial material = sim::ah36Steel();
    const double youngs = material.youngsModulus;
    const double middle = 0.5 * (options.from + options.to);
    const sim::HullGirderSection girder = sim::hullGirderSection(structure, middle);

    std::printf("ferry structure: %zu panels, %zu members, frame spacing %.3f m\n",
                structure.panels.size(), structure.members.size(), structure.frameSpacing);
    std::printf("Tier 0 at x = %.2f m: A = %.5f m^2, neutral axis %.5f m, I = %.5f m^4\n",
                middle, girder.area, girder.neutralAxis, girder.secondMoment);
    std::printf("             so EA = %.5e N, EI = %.5e N m^2\n\n", youngs * girder.area,
                youngs * girder.secondMoment);

    // --- What the mesher built -----------------------------------------------------

    sim::section::SectionParams params;
    params.xFrom = options.from;
    params.xTo = options.to;
    params.subdivision = options.subdivision;

    const double built = now();
    const sim::section::Section hold = sim::section::buildSection(structure, params);
    const double meshSeconds = now() - built;
    std::printf("=== section [%.2f, %.2f] m, subdivision %d, meshed in %.3f s ===\n", options.from,
                options.to, options.subdivision, meshSeconds);
    if (hold.empty()) {
        for (const std::string& problem : hold.problems) std::printf("  ! %s\n", problem.c_str());
        return 1;
    }
    std::printf("  %zu elements, %zu nodes, DOF half-bandwidth %zu\n", hold.elementCount(),
                hold.nodeCount(), hold.halfBandwidth);
    std::printf("  %zu interface nodes (%zu boundary DOF), %d surfaces, %d components,"
                " %d floating, %d spanning\n",
                hold.interfaceNodes.size(), 3 * hold.interfaceNodes.size(), hold.surfaces,
                hold.components, hold.floatingComponents, hold.spanningComponents);
    std::printf("  %.1f m^2 of mid-surface, %.0f kg of plate + %.0f kg of member = %.0f kg\n",
                hold.area, hold.plateMass, hold.memberMass, hold.mass());
    std::printf("  worst Jacobian %.3e, worst aspect %.2f, normal spread %.4f rad on %d elements"
                " (%.0f%% excess plate bending)\n",
                hold.worstJacobian, hold.worstAspect, hold.worstNormalSpread,
                hold.distortedElements, 100.0 * hold.spuriousStiffness);
    std::printf("  %d tapered elements, worst dt/t %.4f (%.0f%% excess plate bending)\n",
                hold.taperedElements, hold.worstTaper, 100.0 * hold.taperStiffness);
    std::printf("  free edge %.1f m, of which %.1f m sits on plating it is not welded to"
                " (worst gap %.4f m)\n",
                hold.freeEdgeLength, hold.junctionEdges, hold.worstJunctionGap);
    std::printf("  junction ties: %d nodes tied joining %.1f m of that edge; %d left on a cut"
                " plane, %d refused as a chain, %d outside a master face\n",
                hold.junctionTies, hold.tiedEdges, hold.junctionsOnInterface,
                hold.junctionsChained, hold.junctionsOutsideFace);
    std::printf("  worst tie: overshoot %.4f of a face, through-thickness weight %.4f\n",
                hold.worstJunctionOvershoot, hold.worstJunctionWeight);
    std::printf("  panels straddling a cut plane: %d\n", hold.straddlingPanels);
    std::printf("  members: %d attached, %d refused, %d missed; effective area attached %.5f"
                " + missed %.5f = %.5f m^2\n",
                hold.membersAttached, hold.membersRefused, hold.membersMissed,
                hold.attachedMemberArea, hold.missedMemberArea,
                hold.attachedMemberArea + hold.missedMemberArea);
    for (const std::string& problem : hold.problems) std::printf("  ! %s\n", problem.c_str());

    // --- Against Tier 0 --------------------------------------------------------------

    const auto report = [&](const sim::section::Section& section, const char* label) {
        sim::section::BeamLoad axial;
        axial.strain = 1e-6;
        axial.reference = girder.neutralAxis;
        const sim::section::BeamResponse stretched =
            sim::section::applyBeamLoad(section, material, axial);
        sim::section::BeamLoad bending;
        bending.curvature = 1e-6;
        bending.reference = girder.neutralAxis;
        const sim::section::BeamResponse bent =
            sim::section::applyBeamLoad(section, material, bending);
        if (!stretched.ok || !bent.ok) {
            std::printf("  %-16s refused: %s %s\n", label, stretched.problem.c_str(),
                        bent.problem.c_str());
            return;
        }
        const double area = stretched.axialStiffness / youngs;
        const double neutralAxis =
            stretched.bendingMoment / stretched.axialForce + girder.neutralAxis;
        const double second = bent.bendingStiffness / youngs;
        const double predicted = girder.area - section.missedMemberArea;
        std::printf("  %-16s A_eff %.5f m^2 (%+.3f%% of Tier 0, %+.3f%% of Tier 0 less the"
                    " members it could not attach)\n",
                    label, area, 100.0 * (area / girder.area - 1.0),
                    100.0 * (area / predicted - 1.0));
        std::printf("  %-16s z_na  %.5f m (Tier 0 %.5f)   I_eff %.5f m^4 (%+.3f%% of Tier 0)\n", "",
                    neutralAxis, girder.neutralAxis, second,
                    100.0 * (second / girder.secondMoment - 1.0));
        std::printf("  %-16s residual %.2e N, rigid-body restraint reaction %.2e N\n", "",
                    stretched.residual, stretched.restraintReaction);
    };

    std::printf("\n=== against Tier 0 ===\n");
    report(hold, "with members");
    sim::section::SectionParams bare = params;
    bare.members = false;
    const sim::section::Section plating = sim::section::buildSection(structure, bare);
    report(plating, "bare plating");
    double stiffenerArea = 0, plateArea = 0;
    for (const sim::SectionElement& element : sim::sectionElements(structure, middle))
        (element.stiffener ? stiffenerArea : plateArea) += element.area;
    std::printf("  Tier 0 says the stiffeners are %.5f m^2 of %.5f, %.1f%% -- which is what"
                " bare plating must be short by\n",
                stiffenerArea, girder.area, 100.0 * stiffenerArea / girder.area);

    // --- Resolution ------------------------------------------------------------------

    // --- What the junction tie costs and buys, at one element per panel -------------
    //
    // The two quantities a prescribed plane-sections field cannot see. `EA` and `EI`
    // are above and they move 0.19% and 0.12%; these move by factors.
    {
        std::printf("\n=== the junction tie: cut against tied ===\n");
        std::printf("  %-6s %7s %6s %9s %9s %8s %12s %10s %10s\n", "", "band", "comps", "tied m",
                    "A_eff", "solve s", "GJ", "first Hz", "reduce s");
        for (int tie = 0; tie < 2; ++tie) {
            sim::section::SectionParams part = params;
            part.junctions = tie != 0;
            const sim::section::Section piece = sim::section::buildSection(structure, part);
            const double solving = now();
            sim::section::BeamLoad axial;
            axial.strain = 1e-6;
            const sim::section::BeamResponse stretched =
                sim::section::applyBeamLoad(piece, material, axial);
            const sim::section::TorsionResponse twisted =
                sim::section::applyTwist(piece, material, 1e-6, girder.neutralAxis);
            const double solveSeconds = now() - solving;
            const double reducing = now();
            const sim::reduction::Substructure substructure(piece.mesh, piece.material,
                                                            piece.interfaceNodes, piece.attachment);
            double hz = 0;
            if (substructure.ready()) {
                const sim::reduction::Eigenpairs modes = substructure.fixedInterfaceModes(1);
                if (!modes.value.empty())
                    hz = std::sqrt(std::max(0.0, modes.value[0])) / (2.0 * std::numbers::pi);
            }
            std::printf("  %-6s %7zu %6d %9.1f %9.5f %8.2f %12.4e %10.4f %10.2f\n",
                        tie ? "tied" : "cut", piece.halfBandwidth, piece.components,
                        piece.tiedEdges, stretched.axialStiffness / youngs, solveSeconds,
                        twisted.torsionalStiffness, hz, now() - reducing);
            std::fflush(stdout);
        }
    }

    // The resolution sweep runs **untied**. It is a study of the mesher's own
    // convergence, and the tie's cost is a band: 146 against 1 520 at subdivision 1,
    // which is 0.16 s of banded factorisation against 5.34, and 278 against 3 188 at
    // subdivision 2, which is 1.1 s against 149. Sweeping tied would take an hour and
    // measure the same four columns to within 0.2%.
    if (options.sweep >= 1) {
        std::printf("\n=== resolution: refining the reduced answer (junctions untied) ===\n");
        std::printf("  %4s %9s %7s %11s %10s %12s %13s %9s\n", "sub", "elements", "band", "A_eff",
                    "z_na", "I_eff", "GJ", "solve s");
        for (int subdivision = 1; subdivision <= options.sweep; ++subdivision) {
            sim::section::SectionParams refined = params;
            refined.subdivision = subdivision;
            refined.junctions = false;
            const sim::section::Section section = sim::section::buildSection(structure, refined);
            const double start = now();
            sim::section::BeamLoad axial;
            axial.strain = 1e-6;
            axial.reference = girder.neutralAxis;
            const sim::section::BeamResponse stretched =
                sim::section::applyBeamLoad(section, material, axial);
            sim::section::BeamLoad bending;
            bending.curvature = 1e-6;
            bending.reference = girder.neutralAxis;
            const sim::section::BeamResponse bent =
                sim::section::applyBeamLoad(section, material, bending);
            const sim::section::TorsionResponse twisted =
                sim::section::applyTwist(section, material, 1e-6, girder.neutralAxis);
            std::printf("  %4d %9zu %7zu %11.5f %10.5f %12.5f %13.4e %9.2f\n", subdivision,
                        section.elementCount(), section.halfBandwidth,
                        stretched.axialStiffness / youngs,
                        stretched.bendingMoment / stretched.axialForce + girder.neutralAxis,
                        bent.bendingStiffness / youngs, twisted.torsionalStiffness, now() - start);
            std::fflush(stdout);
        }
    }

    // --- The reduction ----------------------------------------------------------------

    if (options.reduce) {
        std::printf("\n=== Craig-Bampton ===\n");
        // Which piece of the section owns the softest fixed-interface mode. Untied,
        // the answer is the decks and it is the *same number* with and without the
        // shell they should be welded to -- the junctions' cost stated as a
        // frequency. Tied, the whole section is stiffer than either piece, which is
        // what a joined structure does. See `section.hpp` §2 and §5.
        struct Case {
            const char* label;
            bool shell, deck, bulkhead, junctions;
        };
        const Case cases[] = {{"shell only", true, false, false, false},
                              {"decks only", false, true, false, false},
                              {"bulkheads only", false, false, true, false},
                              {"whole, untied", true, true, true, false},
                              {"whole, tied", true, true, true, true}};
        for (const Case& one : cases) {
            sim::section::SectionParams part = params;
            part.shell = one.shell;
            part.deck = one.deck;
            part.bulkhead = one.bulkhead;
            part.junctions = one.junctions;
            const sim::section::Section piece = sim::section::buildSection(structure, part);
            if (piece.empty()) continue;
            const double assembled = now();
            const sim::reduction::Substructure substructure(piece.mesh, piece.material,
                                                            piece.interfaceNodes, piece.attachment);
            const double assembleSeconds = now() - assembled;
            if (!substructure.ready()) {
                std::printf("  %-15s substructure refused\n", one.label);
                for (const std::string& problem : substructure.problems())
                    std::printf("      ! %s\n", problem.c_str());
                continue;
            }
            sim::reduction::ReduceParams reduceParams;
            reduceParams.modes = one.shell && one.deck && one.bulkhead ? options.modes : 0;
            reduceParams.verifyModes = reduceParams.modes != 0;
            const double reducing = now();
            const sim::reduction::Reduction reduced =
                sim::reduction::craigBampton(substructure, reduceParams);
            const double reduceSeconds = now() - reducing;
            std::printf("  %-15s %zu boundary + %zu interior DOF, band %zu, mass %.0f kg"
                        " (attached %.0f), assemble %.2f s\n",
                        one.label, substructure.boundaryCount(), substructure.interiorCount(),
                        substructure.halfBandwidth(), substructure.totalMass(),
                        substructure.attachedMass(), assembleSeconds);
            std::printf("  %-15s first fixed-interface mode %.4f rad/s = %.4f Hz;"
                        " %d modes kept in %.2f s\n",
                        "", reduced.firstFixedFrequency,
                        reduced.firstFixedFrequency / (2.0 * std::numbers::pi), reduced.modes,
                        reduceSeconds);
            // **The subspace iteration says it did not converge, so the frequency
            // above is not evidence on its own.** `eigenvaluesBelow` counts rather
            // than converges -- the inertia of an LDL^T factorisation, Sylvester's
            // law -- so bracketing the reported value with it is a different
            // instrument answering the same question, which is the only kind of
            // agreement worth having here.
            const double omega = reduced.firstFixedFrequency;
            bool exactBelow = false, exactAbove = false;
            const int below = substructure.eigenvaluesBelow(0.99 * 0.99 * omega * omega, &exactBelow);
            const int above = substructure.eigenvaluesBelow(1.01 * 1.01 * omega * omega, &exactAbove);
            std::printf("  %-15s inertia count brackets it: %d modes below %.4f rad/s, %d below"
                        " %.4f (exact %d/%d)\n",
                        "", below, 0.99 * omega, above, 1.01 * omega, exactBelow, exactAbove);
            for (const std::string& problem : reduced.problems)
                std::printf("      ! %s\n", problem.c_str());
            std::fflush(stdout);
        }
    }
    return 0;
}
