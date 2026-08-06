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
//                   [--chain=N] [--match=M] [--scan=BAYS] [--invariance=BAYS]
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
    // Cut [from, to] into this many sections, reduce each and assemble them, against
    // the same length in one piece. 0 skips it. It is off by default because the
    // assembled model is dense at (N+1) x 1 170 boundary DOF on this ship and the
    // Cholesky at the end of it is the most expensive thing this program does.
    int chain = 0;
    // How far apart two boundary DOF on a shared cut plane may be and still be the
    // same unknown. `matchBoundaries`' 1e-9 default used to be right amidships and to
    // leave a third of the interface unmatched at the ends -- see `section.hpp` §6
    // note 1 and the `--from=36 --to=48` run, which is what measured it. Since the
    // halo of §8 the two sections put that plane's nodes at the same double, so it is
    // not a tolerance any more and the default matches everything. It stays an option
    // because `--invariance` runs the mesher with the halo off, and that is the
    // configuration it was ever needed for.
    double match = 1e-9;
    // **The reach.** Slide a window this many frame bays long along the whole ship
    // and report, station by station, whether the mesher delivers something that
    // reduces and solves. 0 skips it. This is the measurement that says how much of
    // the ship a Tier-1 model can be built on, and before the collapsed-element
    // work (`section.hpp` §7) the answer was 21 windows of 49 -- 62.4 m of 120 in
    // five islands, whose longest unbroken run was 26.4 m.
    int scan = 0;
    // **The invariance.** At every frame station of the ship, mesh the window aft of
    // it, the window forward of it and the window spanning it, and ask whether the
    // three agree about the station they share. They did not: a node is the
    // mid-surface offset by t/2 along the nodal normal, both of which `buildSection`
    // averaged over the sub-quads *inside* the window, so the same node came out in a
    // different place from either side -- up to 1.0 mm at the bow -- and the
    // stiffener steel of a window depended on where it was cut, by 25.2% at the bow
    // shoulder. `SectionParams::halo` is the fix; this is the sweep that says so
    // along the whole hull rather than at the two stations the unit suite can afford.
    // 0 skips it. Meshing only: no solve, no reduction.
    int invariance = 0;
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
        else if (const char* v = value("chain")) o.chain = std::atoi(v);
        else if (const char* v = value("match")) o.match = std::atof(v);
        else if (const char* v = value("scan")) o.scan = std::atoi(v);
        else if (const char* v = value("invariance")) o.invariance = std::atoi(v);
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

    // --- The invariance: a section is worth what it is worth wherever it was cut ----
    //
    // Slide along the hull a station at a time and, at each, mesh the window aft of
    // it, the window forward of it and the window spanning both. Two questions, and
    // they fail in different places:
    //
    //   * **Do the two halves put the shared plane's nodes in the same place?** They
    //     did not, wherever the hull turns: the nodal normal was averaged over the
    //     sub-quads inside the window, so it leaned aft for one section and forward
    //     for the other. Reported as the furthest any node on the plane is from its
    //     partner across it -- the same question `matchBoundaries` asks, asked
    //     without a tolerance.
    //   * **Does the spanning window carry the steel its two halves carry?** It did
    //     not, wherever the plating steps: the nodal thickness was averaged the same
    //     way, so a run of stiffener broke in the spanning window and not in either
    //     half. Reported as the relative disagreement in `memberMass`.
    //
    // `SectionParams::halo` averages both over one bay beyond each plane and meshes
    // only what is inside. The node gap then reads **exactly zero**, and "exactly" is
    // the claim rather than a rounding: the same node is reached by the same panels in
    // the same order whichever window contains it, so its normal and its thickness are
    // formed from the same terms and come out at the same double. The steel is the
    // same fibres summed in a different order and agrees to 3.9e-14 of itself, which
    // is why that one is counted against 1e-9 and the node against nothing at all.
    if (options.invariance > 0) {
        double lo = 1e300, hi = -1e300;
        for (const sim::PlatePanel& p : structure.panels)
            for (int c = 0; c < 4; ++c) {
                lo = std::min(lo, p.corner[c].x);
                hi = std::max(hi, p.corner[c].x);
            }
        const double bay = structure.frameSpacing;
        const double span = options.invariance * bay;
        std::printf("=== invariance: %g m windows either side of every station, halo on and off"
                    " ===\n", span);
        std::printf("  %8s %6s %7s %13s %13s %11s %11s\n", "station", "halo", "nodes", "worst gap m",
                    "steel halves", "one piece", "relative");
        const auto cut = [&](double from, double to, bool halo) {
            sim::section::SectionParams p;
            p.xFrom = from;
            p.xTo = to;
            p.subdivision = options.subdivision;
            p.junctions = false;  // the tie moves no node and this is about nodes
            p.halo = halo;
            return sim::section::buildSection(structure, p);
        };
        int stations = 0, moved[2] = {0, 0}, lost[2] = {0, 0}, withHalo = 0;
        double worstGap[2] = {0, 0}, worstSteel[2] = {0, 0};
        for (double x = std::floor(lo / bay + 0.5) * bay + span; x + span <= hi + 1e-9; x += bay) {
            ++stations;
            for (int halo = 1; halo >= 0; --halo) {
                const sim::section::Section aft = cut(x - span, x, halo != 0);
                const sim::section::Section forward = cut(x, x + span, halo != 0);
                const sim::section::Section whole = cut(x - span, x + span, halo != 0);
                if (aft.empty() || forward.empty() || whole.empty()) continue;
                if (halo == 1 && whole.haloPanels > 0) ++withHalo;
                // Every node of the aft section's forward plane, against the nearest
                // node of the forward section's aft plane. Nearest and not "the same
                // index": the two meshes are numbered for bandwidth and share no
                // numbering at all.
                const auto place = [](const sim::section::Section& s, std::uint32_t node) {
                    return sim::Vec3{s.mesh.position[node * 3], s.mesh.position[node * 3 + 1],
                                     s.mesh.position[node * 3 + 2]};
                };
                double gap = 0;
                for (std::uint32_t a : aft.forwardNodes) {
                    double nearest = 1e300;
                    for (std::uint32_t b : forward.aftNodes)
                        nearest = std::min(nearest, sim::length(place(aft, a) - place(forward, b)));
                    gap = std::max(gap, nearest);
                }
                const double halves = aft.memberMass + forward.memberMass;
                const double steel = halves > 0 ? std::abs(whole.memberMass / halves - 1.0) : 0.0;
                const auto slot = static_cast<std::size_t>(halo);
                worstGap[slot] = std::max(worstGap[slot], gap);
                worstSteel[slot] = std::max(worstSteel[slot], steel);
                if (gap > 1e-9) ++moved[slot];
                if (steel > 1e-9) ++lost[slot];
                if (halo == 0 && (gap > 1e-9 || steel > 1e-9))
                    std::printf("  %8.1f %6d %7zu %13.3e %13.1f %11.1f %+11.3e\n", x, halo,
                                aft.forwardNodes.size(), gap, halves, whole.memberMass,
                                halves > 0 ? whole.memberMass / halves - 1.0 : 0.0);
            }
            std::fflush(stdout);
        }
        std::printf("\n  %d stations swept at subdivision %d\n", stations, options.subdivision);
        for (int halo = 1; halo >= 0; --halo) {
            const auto slot = static_cast<std::size_t>(halo);
            std::printf("  halo %s: %3d stations move a node (worst %.3e m), %3d lose stiffener"
                        " steel (worst %.3e)\n",
                        halo ? "on " : "off", moved[slot], worstGap[slot], lost[slot],
                        worstSteel[slot]);
        }

        // --- The success contract -------------------------------------------------
        //
        // Zero, and two guards against zero being free. The first is that the halo
        // was made of something; the second is that the *same sweep without it* is
        // not zero, because a hull with no shoulder and no strake step would report a
        // clean sweep from a mesher that had never been fixed.
        if (moved[1] != 0 || lost[1] != 0) {
            std::printf("       ! %d stations put a node in a different place depending on the"
                        " window (worst %.3e m) and %d carry different stiffener steel (worst"
                        " %.3e): a section is not a property of the ship\n",
                        moved[1], worstGap[1], lost[1], worstSteel[1]);
            return 1;
        }
        if (withHalo != stations) {
            std::printf("       ! only %d of %d stations had any plating beyond the cut to"
                        " average, so the sweep is a statement about nothing\n",
                        withHalo, stations);
            return 1;
        }
        if (moved[0] == 0 || lost[0] == 0) {
            std::printf("       ! without the halo the same sweep moves %d nodes and loses steel"
                        " at %d stations. Both have to be non-zero or this ship cannot tell a"
                        " fixed mesher from an unfixed one\n",
                        moved[0], lost[0]);
            return 1;
        }
        if (stations < 8) {
            std::printf("       ! %d stations is not a ship\n", stations);
            return 1;
        }
        std::printf("\nok\n");
        return 0;
    }

    // --- The reach: how much of the ship can be meshed at all ----------------------
    //
    // Everything below this measures one section. This measures how many sections
    // there are: a window of `--scan` frame bays slid along the whole hull, asked
    // for a mesh, a reduction and two solves at every station. A window counts as
    // reached only when all four succeed *and* the mesh is one connected piece --
    // `EA` alone would score a mesh that joined nothing (see `section.hpp` §2), so
    // it is reported but is not the criterion.
    if (options.scan > 0) {
        double lo = 1e300, hi = -1e300;
        for (const sim::PlatePanel& p : structure.panels)
            for (int c = 0; c < 4; ++c) {
                lo = std::min(lo, p.corner[c].x);
                hi = std::max(hi, p.corner[c].x);
            }
        const double bay = structure.frameSpacing;
        const double span = options.scan * bay;
        std::printf("=== reach: a %g m window (%d bays) along a hull of %.1f m, subdivision %d ===\n",
                    span, options.scan, hi - lo, options.subdivision);
        std::printf("  %8s %8s %8s %7s %7s %6s %11s %9s %9s %10s %9s %8s\n", "xFrom", "xTo",
                    "elems", "collps", "invert", "comps", "minGaussJ", "A_eff m2", "Tier0 m2",
                    "GJ", "first Hz", "verdict");
        // Which bays are inside at least one working window. Consecutive windows
        // overlap by `span - bay`, so summing their lengths would count the middle
        // body several times over; a union is the only honest reach.
        const int bays = static_cast<int>(std::lround((hi - lo) / bay));
        std::vector<bool> covered(static_cast<std::size_t>(std::max(bays, 0)), false);
        int windows = 0, good = 0, usable = 0, collapsedSeen = 0;
        double firstGood = 1e300, lastGood = -1e300, bestTorsion = 0;
        for (double x = std::floor(lo / bay + 0.5) * bay; x + span <= hi + 1e-9; x += bay) {
            sim::section::SectionParams window;
            window.xFrom = x;
            window.xTo = x + span;
            window.subdivision = options.subdivision;
            const sim::section::Section piece = sim::section::buildSection(structure, window);
            ++windows;
            // **Tier 0 at the centre of a bay, not at the window's midpoint.** An
            // even number of bays puts that midpoint on a frame station, and
            // `hullGirderSection` sampled exactly on a panel seam loses every panel
            // either side of it: 0.42932 m^2 against 1.80133 amidships, a 76%
            // shortfall that is a knife-edge on the half-open `straddles` test in
            // `sectionElements` rather than anything about the ship. It is not this
            // file's defect to fix -- and it is why the column below reports both
            // areas rather than their ratio.
            const sim::HullGirderSection tier0 =
                sim::hullGirderSection(structure, x + 0.5 * (span - bay));
            double effective = 0, gj = 0, hz = 0;
            bool solved = false, reduces = false;
            if (!piece.empty()) {
                sim::section::BeamLoad axial;
                axial.strain = 1e-6;
                axial.reference = tier0.neutralAxis;
                const sim::section::BeamResponse stretched =
                    sim::section::applyBeamLoad(piece, material, axial);
                const sim::section::TorsionResponse twisted =
                    sim::section::applyTwist(piece, material, 1e-6, tier0.neutralAxis);
                solved = stretched.ok && twisted.ok;
                if (stretched.ok) effective = stretched.axialStiffness / youngs;
                if (twisted.ok) gj = twisted.torsionalStiffness;
                const sim::reduction::Substructure substructure(
                    piece.mesh, piece.material, piece.interfaceNodes, piece.attachment);
                reduces = substructure.ready();
                if (reduces) {
                    const sim::reduction::Eigenpairs modes = substructure.fixedInterfaceModes(1);
                    if (!modes.value.empty())
                        hz = std::sqrt(std::max(0.0, modes.value[0])) / (2.0 * std::numbers::pi);
                }
            }
            // Two verdicts, because they answer different questions. **Usable** is
            // what a Tier-1 model needs: it meshes, it reduces, it solves, and no
            // piece of it floats free of the interface (a floating component is a
            // mechanism `reduction::Substructure` does not catch -- `section.hpp`).
            // **One piece** is stricter and is what a chain wants; a section in two
            // spanning pieces still reduces and still carries load, it just has a
            // junction the mesher declined to close.
            const bool works = !piece.empty() && piece.invertedElements == 0 && solved && reduces &&
                               piece.floatingComponents == 0;
            const bool ok = works && piece.components == 1;
            if (works) ++usable;
            if (piece.collapsedElements > 0) ++collapsedSeen;
            bestTorsion = std::max(bestTorsion, gj);
            if (ok) {
                ++good;
                firstGood = std::min(firstGood, window.xFrom);
                lastGood = std::max(lastGood, window.xTo);
                for (int b = 0; b < options.scan; ++b) {
                    const auto index = static_cast<std::size_t>(std::lround((x - lo) / bay) + b);
                    if (index < covered.size()) covered[index] = true;
                }
            }
            std::printf("  %8.2f %8.2f %8zu %7d %7d %6d %11.3e %9.5f %9.5f %10.3e %9.4f %8s\n",
                        window.xFrom, window.xTo, piece.elementCount(), piece.collapsedElements,
                        piece.invertedElements, piece.components, piece.worstJacobian, effective,
                        tier0.area, gj, hz,
                        ok ? "ok" : (piece.empty() ? "empty" : (works ? "npieces" : "REFUSED")));
            std::fflush(stdout);
        }
        double reached = 0;
        for (bool b : covered)
            if (b) reached += bay;
        std::printf("\n  %d of %d windows mesh, reduce and solve with nothing floating (%.1f%%)\n",
                    usable, windows, 100.0 * usable / windows);
        std::printf("  %d of %d are also a single connected piece (%.1f%%)\n", good, windows,
                    100.0 * good / windows);
        std::printf("  outermost working cut planes x = %.1f .. %.1f m\n",
                    good > 0 ? firstGood : 0.0, good > 0 ? lastGood : 0.0);
        std::printf("  reach: %.1f m of a %.1f m ship (%.1f%%)\n", reached, hi - lo,
                    100.0 * reached / (hi - lo));

        // --- The success contract -------------------------------------------------
        //
        // The reach is the claim this mode exists to make, so it is the thing the
        // gate checks. **With two guards against it being vacuous**, because "every
        // window worked" is also what an empty ship, or a ship with no degenerate
        // panels on it, would report -- and the second of those is exactly the state
        // this code was in before the collapsed element was understood.
        if (usable != windows) {
            std::printf("       ! %d of %d windows do not mesh, reduce and solve\n",
                        windows - usable, windows);
            return 1;
        }
        if (reached < hi - lo - 0.5 * bay) {
            std::printf("       ! the mesher reaches %.1f m of a %.1f m ship\n", reached, hi - lo);
            return 1;
        }
        if (collapsedSeen == 0) {
            std::printf("       ! not one window contained a collapsed element, so the reach"
                        " measured nothing: the ship this ran on has no degenerate plating and"
                        " the whole comparison is vacuous\n");
            return 1;
        }
        if (windows < 8 || !(bestTorsion > 1e11)) {
            std::printf("       ! %d windows and a best GJ of %.3e: the sections carry no real"
                        " stiffness, so 'every window worked' is a statement about nothing\n",
                        windows, bestTorsion);
            return 1;
        }
        std::printf("\nok\n");
        return 0;
    }

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
    std::printf("  panels straddling a cut plane: %d; halo panels averaged into the nodal"
                " normals and thicknesses and then dropped: %d\n",
                hold.straddlingPanels, hold.haloPanels);
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
    // What the junction tie is worth, kept for the success contract at the end of
    // main: `EA` cannot see whether the plating is joined (§2 of `section.hpp`) and
    // torsion and the lowest fixed-interface frequency can, so those two are what
    // the gate is allowed to assert on.
    double torsionCut = 0, torsionTied = 0, hertzCut = 0, hertzTied = 0;
    int componentsCut = 0, componentsTied = 0;
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
            (tie ? torsionTied : torsionCut) = twisted.torsionalStiffness;
            (tie ? hertzTied : hertzCut) = hz;
            (tie ? componentsTied : componentsCut) = piece.components;
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

    // --- A ship: the same length as a chain of sections ------------------------------
    //
    // The end-to-end claim of the tier, at ship scale: cut a length into N pieces,
    // reduce each once, assemble, and get the model the same length in one piece
    // gives. The reference here owes nothing to the assembly -- it is
    // `applyBeamLoad` and `applyTwist` on the monolithic section, through
    // `solidshell::solveStatic`.
    //
    // Everything is untied unless asked otherwise, because tied at ship scale is a
    // band of 1 520 and a 5.3 s banded factorisation per section, and the two
    // questions -- does the assembly reproduce the monolith, and what do the cut
    // planes cost the ties -- are answered separately below.
    if (options.chain > 0) {
        std::printf("\n=== a chain of %d sections against the same length in one piece ===\n",
                    options.chain);
        for (int tie = 0; tie < 2; ++tie) {
            sim::section::ChainParams chainParams;
            chainParams.section = params;
            chainParams.section.junctions = tie != 0;
            chainParams.reduce.modes = 0;
            chainParams.reduce.cutoffFrequency = 0;
            chainParams.matchTolerance = options.match;
            for (int i = 0; i <= options.chain; ++i)
                chainParams.station.push_back(options.from + (options.to - options.from) * i /
                                                                 options.chain);

            const double start = now();
            const sim::section::Chain chain = sim::section::buildChain(structure, chainParams);
            std::printf("  %-5s built in %.2f s (mesh %.2f, reduce %.2f, assemble %.2f):"
                        " %d pieces, %d assembled boundary DOF, size %d\n",
                        tie ? "tied" : "cut", now() - start, chain.meshSeconds,
                        chain.reduceSeconds, chain.assembleSeconds, chain.components,
                        chain.assembly.boundary, chain.assembly.size());
            std::size_t unmatched = 0;
            for (std::size_t u : chain.unmatched) unmatched += u;
            std::printf("  %-5s interior planes: %zu shared DOF each, %zu unmatched, worst gap"
                        " %.3e m; ties %.1f m of %.1f m of junction edge\n",
                        "", chain.shared.empty() ? 0 : chain.shared.front(), unmatched,
                        chain.worstGap, chain.tiedEdges, chain.junctionEdges);
            for (const std::string& problem : chain.problems)
                std::printf("      ! %s\n", problem.c_str());
            std::fflush(stdout);
            if (!chain.ready()) continue;

            sim::section::SectionParams monoParams = params;
            monoParams.junctions = tie != 0;
            const sim::section::Section mono = sim::section::buildSection(structure, monoParams);
            std::printf("  %-5s one piece: %d components, ties %.1f m of %.1f m\n", "",
                        mono.components, mono.tiedEdges, mono.junctionEdges);

            sim::section::BeamLoad axial;
            axial.strain = 1e-6;
            axial.reference = girder.neutralAxis;
            sim::section::BeamLoad bending;
            bending.curvature = 1e-6;
            bending.reference = girder.neutralAxis;

            double at = now();
            const sim::section::BeamResponse chainAxial =
                sim::section::applyBeamLoad(chain, axial);
            const double chainSolve = now() - at;
            at = now();
            const sim::section::BeamResponse monoAxial =
                sim::section::applyBeamLoad(mono, material, axial);
            const double monoSolve = now() - at;
            const sim::section::BeamResponse chainBend = sim::section::applyBeamLoad(chain, bending);
            const sim::section::BeamResponse monoBend =
                sim::section::applyBeamLoad(mono, material, bending);
            const sim::section::TorsionResponse chainTwist =
                sim::section::applyTwist(chain, 1e-6, girder.neutralAxis);
            const sim::section::TorsionResponse monoTwist =
                sim::section::applyTwist(mono, material, 1e-6, girder.neutralAxis);
            if (!chainAxial.ok || !monoAxial.ok) {
                std::printf("      ! %s / %s\n", chainAxial.problem.c_str(),
                            monoAxial.problem.c_str());
                continue;
            }
            std::printf("  %-5s %-8s %14s %14s %10s\n", "", "", "chain", "one piece", "relative");
            const auto line = [&](const char* what, double a, double b) {
                std::printf("  %-5s %-8s %14.6e %14.6e %+10.3e\n", "", what, a, b,
                            b != 0 ? a / b - 1.0 : 0.0);
            };
            line("EA", chainAxial.axialStiffness, monoAxial.axialStiffness);
            line("z_na", chainAxial.bendingMoment / chainAxial.axialForce + girder.neutralAxis,
                 monoAxial.bendingMoment / monoAxial.axialForce + girder.neutralAxis);
            line("EI", chainBend.bendingStiffness, monoBend.bendingStiffness);
            line("GJ", chainTwist.torsionalStiffness, monoTwist.torsionalStiffness);
            std::printf("  %-5s dense assembled solve %.2f s against a banded monolithic %.2f s;"
                        " chain restraint reaction %.2e N, residual %.2e N\n",
                        "", chainSolve, monoSolve, chainAxial.restraintReaction,
                        chainAxial.residual);
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

    // --- The success contract ------------------------------------------------------
    //
    // **This program had none until it was put in the gate**, and it published a
    // ship-scale frequency that nothing ever re-ran. What it asserts is chosen the
    // same way `section.hpp` §2 chooses what a test may assert on: `EA` and `EI` are
    // exact on a section whose plating is joined to nothing, so they are checked
    // only loosely and against Tier 0, while everything about *joining* is checked
    // on torsion, on the lowest fixed-interface frequency and on the component
    // count -- the three quantities a prescribed plane-sections field cannot see.
    if (hold.invertedElements != 0 || !(hold.worstJacobian > 0)) {
        std::printf("       ! the section has %d inverted elements and a worst Gauss determinant"
                    " of %.3e\n", hold.invertedElements, hold.worstJacobian);
        return 1;
    }
    if (hold.elementCount() < 1000 || !(hold.area > 100.0)) {
        std::printf("       ! %zu elements over %.1f m2 is not a ship section, so nothing below"
                    " means anything\n", hold.elementCount(), hold.area);
        return 1;
    }
    if (!(torsionCut > 0) || !(torsionTied > 1.3 * torsionCut)) {
        std::printf("       ! tying the junctions moved GJ from %.4e to %.4e. A tie that does not"
                    " close the cell has not joined the plating, and EA would not have said so\n",
                    torsionCut, torsionTied);
        return 1;
    }
    if (!(hertzCut > 0) || !(hertzTied > 1.3 * hertzCut)) {
        std::printf("       ! the first fixed-interface mode went %.4f -> %.4f Hz. Untied it is"
                    " the decks' own frequency; a tie that does not raise it has joined nothing\n",
                    hertzCut, hertzTied);
        return 1;
    }
    if (componentsTied != 1 || componentsCut <= 1) {
        std::printf("       ! components went %d untied -> %d tied, and the pair has to be"
                    " many -> one or the tie is being credited with a mesh that was already"
                    " joined\n", componentsCut, componentsTied);
        return 1;
    }
    std::printf("\nok\n");
    return 0;
}
