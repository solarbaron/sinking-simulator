// SPDX-License-Identifier: MIT
//
// Validation of the Tier-1 section mesher.
//
// Two reference structures, and they answer different questions on purpose.
//
//   * A **rectangular box girder**, authored here, four flat plates that really do
//     share their corners. Its area, neutral axis, second moment and Bredt torsion
//     constant are closed forms in `B`, `H`, `t` and `L`, so almost nothing in the
//     first half of this file is a tolerance on an eyeballed number. It is also the
//     only conforming multi-plate input that exists -- `makeStructuralMesh` shares
//     no corner between two panel roles at all -- so it is the only way to ask what
//     welding a junction would do if welding one were possible.
//   * The **reference ferry**, where the point is that the answer survives real
//     geometry and where the independent reference is `hullGirderSection`, which
//     reaches the same numbers by summing `A`, `A z` and `A z^2` over a transverse
//     cut and shares no code with anything here.
//
// **The vacuity guards carry as much of the weight as the assertions**, because the
// central finding of `section.hpp` §2 is that the obvious test proves the wrong
// thing: a section's `EA` and `EI` come out *exact* on a mesh whose plates are not
// joined at all, so a mesher that welded nothing would pass a section-properties
// test with the best score in the file. Every junction claim here is therefore
// checked against torsion or against a fixed-interface frequency as well, and every
// "these two agree" is checked against a third case where they visibly do not.
#include "engine/sim/section.hpp"
#include "engine/sim/constraint.hpp"
#include "engine/sim/girder.hpp"
#include "engine/sim/reduction.hpp"
#include "engine/sim/scantlings.hpp"
#include "game/prototype/ferry.hpp"
#include "harness.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <numbers>
#include <string>
#include <vector>

using namespace sim;
using testing::expectEqual;
using testing::expectNear;
using testing::expectTrue;

namespace {

void expectEqualCount(const std::string& what, std::size_t got, std::size_t want) {
    expectEqual(what, static_cast<long long>(got), static_cast<long long>(want));
}

// --- The box girder --------------------------------------------------------------
//
// Mid-surface dimensions, so every closed form below is exact for the idealisation
// the mesher is built on: corners lie on the plate's mid-surface (`scantlings.hpp`).
// The flanges may carry a different thickness over their outboard half, which is how
// the thickness-seam tests get a seam.

constexpr double kB = 2.0, kH = 1.0, kL = 8.0, kT = 0.010;
constexpr int kNx = 16, kNy = 4, kNz = 2;

// `alternateWinding` reverses every other bay's corner order **within each of the
// four surfaces**, which is what `makeStructuralMesh` does when it mirrors the
// starboard side: the panels of one surface do not all wind the same way round, and
// the mesher has to orient them against one another before averaging their normals.
// Reversing whole plates instead would prove nothing -- the surface would simply
// face the other way and every normal would still agree.
StructuralMesh makeBox(double thicknessInner, double thicknessOuter,
                       bool alternateWinding = false) {
    StructuralMesh mesh;
    mesh.materials = {ah36Steel()};
    mesh.frameSpacing = kL / kNx;
    const auto quad = [&](const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d, double t,
                          bool reversed) {
        PlatePanel p;
        p.corner[0] = a;
        p.corner[1] = reversed ? d : b;
        p.corner[2] = c;
        p.corner[3] = reversed ? b : d;
        p.thickness = t;
        p.material = 0;
        p.role = PanelRole::Shell;
        p.source = 0;
        mesh.panels.push_back(p);
    };
    for (int i = 0; i < kNx; ++i) {
        const bool reversed = alternateWinding && (i % 2 == 1);
        const double x0 = kL * i / kNx, x1 = kL * (i + 1) / kNx;
        for (int j = 0; j < kNy; ++j) {
            const double y0 = -kB / 2 + kB * j / kNy, y1 = -kB / 2 + kB * (j + 1) / kNy;
            const double t = j < kNy / 2 ? thicknessInner : thicknessOuter;
            quad({x0, y0, 0}, {x1, y0, 0}, {x1, y1, 0}, {x0, y1, 0}, t, reversed);
            quad({x0, y0, kH}, {x1, y0, kH}, {x1, y1, kH}, {x0, y1, kH}, t, reversed);
        }
        for (int k = 0; k < kNz; ++k) {
            const double z0 = kH * k / kNz, z1 = kH * (k + 1) / kNz;
            quad({x0, -kB / 2, z0}, {x1, -kB / 2, z0}, {x1, -kB / 2, z1}, {x0, -kB / 2, z1},
                 thicknessInner, reversed);
            quad({x0, kB / 2, z0}, {x1, kB / 2, z0}, {x1, kB / 2, z1}, {x0, kB / 2, z1},
                 thicknessInner, reversed);
        }
    }
    for (int i = 0; i <= kNx; ++i) mesh.frameStations.push_back(kL * i / kNx);
    return mesh;
}

// **The same box with every panel cut down its own diagonal**, which turns each
// quad into two `PlatePanel`s whose fourth corner sits on their first. That is
// exactly the shape `makeStructuralMesh` emits where a bulkhead runs into the keel
// or a deck strake runs out at the stem -- 166 of the reference ferry's 8 900
// panels -- and extruding one gives a **collapsed hexahedron**: a wedge, whose
// Jacobian is zero on the closed edge and sound everywhere it is integrated. See
// `section.hpp` §7.
//
// The split conforms at any subdivision. Panel `(c0, c1, c2, c0)` divides its
// `c0->c1`, `c1->c2` and `c0->c2` edges into `n`, and `(c0, c2, c3, c0)` divides
// `c0->c2` the same way round, so the diagonal the two share is cut identically and
// the mesh needs no seam matching -- the same argument `SectionParams::subdivision`
// makes for whole panels.
StructuralMesh triangulated(const StructuralMesh& source) {
    StructuralMesh out = source;
    out.panels.clear();
    for (const PlatePanel& p : source.panels) {
        PlatePanel a = p, b = p;
        a.corner[0] = p.corner[0];
        a.corner[1] = p.corner[1];
        a.corner[2] = p.corner[2];
        a.corner[3] = p.corner[0];
        b.corner[0] = p.corner[0];
        b.corner[1] = p.corner[2];
        b.corner[2] = p.corner[3];
        b.corner[3] = p.corner[0];
        out.panels.push_back(a);
        out.panels.push_back(b);
    }
    return out;
}

section::SectionParams boxParams(int subdivision = 1) {
    section::SectionParams p;
    p.xFrom = 0.0;
    p.xTo = kL;
    p.subdivision = subdivision;
    p.members = false;
    return p;
}

// The same box with the junction tie switched off -- what the mesher delivered
// before §5 existed, and the negative control every junction claim below is made
// against. It is also the *right* configuration for anything that is not about the
// junctions: the tie joins the four corners, so a test of the welder or of the
// closed-form `EI` would otherwise be measuring two things at once.
section::SectionParams openBoxParams(int subdivision = 1) {
    section::SectionParams p = boxParams(subdivision);
    p.junctions = false;
    return p;
}

// And with the corners welded instead: one node pair where two plates meet, which
// is what `section.hpp` §1 refuses to build and §5 prices the tie against.
section::SectionParams weldedBoxParams(int subdivision = 1) {
    section::SectionParams p = openBoxParams(subdivision);
    p.foldLimit = 2.0;
    return p;
}

// Uniform-thickness closed forms.
constexpr double kBoxArea = 2 * (kB + kH) * kT;
// The plates' own second moment about their own centroids is in this and is *not*
// negligible at the tolerance the mesh reaches: it is 2.857e-5 of the total, which
// is 58 times the error the mesh actually has. The mid-surface form without it is
// what a hand calculation would produce and it is used below as the vacuity guard.
constexpr double kBoxSecondMomentMidSurface = 2 * (kB * kT * kH * kH / 4) + 2 * (kT * kH * kH * kH / 12);
constexpr double kBoxSecondMoment = kBoxSecondMomentMidSurface + 2 * (kB * kT * kT * kT / 12);
// Bredt: GJ = G * 4 A_enclosed^2 / integral(ds/t).
constexpr double kBoxTorsionConstant = 2 * kB * kB * kH * kH * kT / (kB + kH);

// The tolerances below are what was **measured**, with a factor of ten or so of
// margin for a different optimisation level -- the gate builds the engine at -O3 and
// the sanitizers at -O1. A loose tolerance is nearly a vacuous one: `EA` at 1e-10
// would pass on a model that had lost a whole plate and was merely well converged.
//
//   EA on the cut box            6.1e-13   asserted at 1e-11
//   EI against the exact form    4.9e-07   asserted at 2e-06 (guard: 2.9e-05)
//   taper against split, EA      4.5e-07   asserted at 2e-06
//   Guyan energy identity        1.4e-10   asserted at 1e-08
//   reduced rigid mass           1.8e-12   asserted at 1e-10
//   a fibre at its own offset    1.8e-13 m asserted at 1e-11 m
constexpr double kAxialTolerance = 1e-11;
constexpr double kBendingTolerance = 2e-6;

// --- 1. What the mesher built ----------------------------------------------------

void testBoxMesh() {
    const StructuralMesh structure = makeBox(kT, kT);
    const section::Section cut = section::buildSection(structure, openBoxParams());

    expectEqualCount("box elements at subdivision 1", cut.elementCount(),
                     static_cast<std::size_t>(kNx) * (2 * kNy + 2 * kNz));
    expectNear("box mid-surface area", cut.area, 2 * (kB + kH) * kL, 1e-12);
    expectNear("box plate mass", cut.plateMass, 2 * (kB + kH) * kL * kT * ah36Steel().density, 1e-6);
    expectTrue("no element is inverted", cut.worstJacobian > 0);
    expectEqual("no element tapers on a uniform box", cut.taperedElements, 0);

    // Four plates meeting at right angles: a fold of pi/2 is past `foldLimit`, so
    // the mesher leaves them as four surfaces rather than welding a node pair that
    // would point 45 degrees out of both plates. See `section.hpp` §1.
    expectEqual("box surfaces at the default fold limit", cut.surfaces, 4);
    expectEqual("box components at the default fold limit", cut.components, 4);
    expectEqual("nothing floats free of the interface", cut.floatingComponents, 0);
    expectEqual("every piece reaches both cut planes", cut.spanningComponents, 4);
    // Four corner lines, each `kL` long, each carrying one free edge from each of the
    // two plates that meet there.
    expectNear("free edge length is the four corner lines, twice", cut.freeEdgeLength, 8 * kL, 1e-9);
    expectNear("and all of it is sitting on plating it is not joined to", cut.junctionEdges,
               8 * kL, 1e-9);
    expectTrue("the junction gap is a weld tolerance, not a real gap", cut.worstJunctionGap < 1e-6);

    // The interface. Both nodes of every through-thickness pair, chosen on the
    // mid-surface, which is exact where choosing on the nodes is not.
    expectTrue("the interface is not empty", !cut.interfaceNodes.empty());
    expectEqualCount("the two cut planes carry the same number of nodes", cut.aftNodes.size(),
                     cut.forwardNodes.size());
    expectEqualCount("the interface is the two planes together", cut.interfaceNodes.size(),
                     cut.aftNodes.size() + cut.forwardNodes.size());
    double worstOffPlane = 0;
    for (std::uint32_t node : cut.interfaceNodes) {
        const double x = cut.mesh.position[static_cast<std::size_t>(node) * 3];
        worstOffPlane = std::max(worstOffPlane, std::min(std::abs(x - 0.0), std::abs(x - kL)));
    }
    expectTrue("every interface node lies on a cut plane", worstOffPlane < 1e-12);

    // Welding the corners: one surface, no free edge, and the price is in §1.
    const section::Section joined = section::buildSection(structure, weldedBoxParams());
    expectEqual("raising the fold limit welds the box into one surface", joined.surfaces, 1);
    expectEqual("and one component", joined.components, 1);
    expectNear("with no free edge left", joined.freeEdgeLength, 0.0, 1e-12);
    expectEqualCount("the two meshes have the same elements", joined.elementCount(),
                     cut.elementCount());
    expectNear("and the same mid-surface area", joined.area, cut.area, 1e-12);
    // The corner node's normal is the mean of two plates 90 degrees apart, so it is
    // 45 degrees out of each; an element with two such corners and two ordinary ones
    // has a mean normal 22.5 degrees from both, and the chord of 22.5 degrees is
    // `2 sin(pi/16)`. A closed form for the geometry §1 refuses to build.
    expectNear("the welded corner turns the thickness direction 22.5 degrees off the element",
               joined.worstNormalSpread, 2.0 * std::sin(std::numbers::pi / 16), 1e-9);

    // --- And the third way, which is the one this file is now about --------------
    //
    // Tying moves no geometry at all: the same nodes, in the same places, still four
    // surfaces, still every corner edge free in the element sense -- and one
    // connected component, because the constraint joins what the weld would have
    // merged. That combination is the whole claim: `surfaces` counts what shares
    // nodes, `components` counts what is joined, and a tie separates the two
    // questions for the first time.
    const section::Section tied = section::buildSection(structure, boxParams());
    expectEqual("tying leaves the box four surfaces", tied.surfaces, 4);
    expectEqual("and one component", tied.components, 1);
    expectEqual("which spans the section", tied.spanningComponents, 1);
    expectEqualCount("with exactly the nodes the untied mesh had", tied.nodeCount(),
                     cut.nodeCount());
    {
        // The same points, not the same indices: tying changes the numbering,
        // because the ordering is chosen knowing which nodes the ties couple. So the
        // two position arrays are compared as sorted sets.
        std::vector<std::array<double, 3>> a, b;
        for (std::size_t n = 0; n < tied.nodeCount(); ++n) {
            a.push_back({tied.mesh.position[n * 3], tied.mesh.position[n * 3 + 1],
                         tied.mesh.position[n * 3 + 2]});
            b.push_back({cut.mesh.position[n * 3], cut.mesh.position[n * 3 + 1],
                         cut.mesh.position[n * 3 + 2]});
        }
        std::sort(a.begin(), a.end());
        std::sort(b.begin(), b.end());
        double worst = 0;
        for (std::size_t n = 0; n < a.size(); ++n)
            for (int k = 0; k < 3; ++k)
                worst = std::max(worst, std::abs(a[n][static_cast<std::size_t>(k)] -
                                                 b[n][static_cast<std::size_t>(k)]));
        expectNear("in exactly the same places", worst, 0.0, 0.0);
    }
    expectNear("and the same mid-surface area", tied.area, cut.area, 1e-12);
    expectNear("and the same steel", tied.plateMass, cut.plateMass, 1e-9);
    // The welded mesh is the contrast: same elements, but its nodal normals have
    // turned, which is exactly the thinning §1 measures.
    expectNear("where welding turns the thickness direction and tying does not",
               tied.worstNormalSpread, 0.0, 1e-12);

    // Four corner lines, 17 stations each, of which the two on the cut planes get no
    // *face* tie: half a face's nodes are interior to the section, and a boundary
    // degree of freedom written as a function of an interior one is not one a
    // reduction keeps exactly.
    expectEqual("every corner station off the cut planes is tied", tied.junctionTies,
                4 * (kNx - 1));
    // They get a line tie instead -- §9 -- and one per corner per plane, because the
    // two coincident nodes at a corner are one junction and the second is the first's
    // master. Nothing is left over: `junctionsOnInterface` counts what no line reached.
    expectEqual("the stations on the cut planes are tied to the line in the plane",
                tied.planeTieNodes, 4 * 2);
    expectEqual("and nothing on a cut plane is left open", tied.junctionsOnInterface, 0);
    // The negative control, and the figure this file used to publish: with the line
    // tie off, every junction node on a cut plane is one the interface leaves open.
    {
        section::SectionParams open = boxParams();
        open.interfaceTies = false;
        const section::Section noLine = section::buildSection(structure, open);
        expectEqual("with the line tie off they are counted, not silently dropped",
                    noLine.junctionsOnInterface, 4 * 2 * 2);
        expectEqual("and none of them is tied", noLine.planeTieNodes, 0);
        expectEqualCount("so the section carries no constraint for them",
                         noLine.planeTies.size(), 0u);
        // And the switch changes nothing else: a line tie is data a section carries,
        // not something it applies, so the mesh and its own ties are bit-identical.
        expectEqual("the face ties are untouched", noLine.junctionTies, tied.junctionTies);
        expectNear("and so is what they join", noLine.tiedEdges, tied.tiedEdges, 0.0);
    }
    expectEqualCount("six constraints per tied node, one per axis per extruded node",
                     tied.planeTies.size(), static_cast<std::size_t>(6 * 4 * 2));
    expectEqual("nothing was refused as a chain", tied.junctionsChained, 0);
    expectEqual("nor as a line-tie chain", tied.planeTiesChained, 0);
    expectEqual("nor for landing off the end of a line", tied.planeTiesOutsideLine, 0);
    expectEqual("nor for having no line to land on", tied.planeTiesUnreached, 0);
    expectEqual("nor for landing outside a face", tied.junctionsOutsideFace, 0);
    // A butt corner puts the flange's node half a plate thickness past the end of
    // the side plate's mid-surface: `t / 2` over a half-metre element is 0.02 of the
    // face's natural coordinate, and it is a closed form rather than a reading.
    expectNear("the overshoot is half a plate thickness over half an element",
               tied.worstJunctionOvershoot, (kT / 2) / (0.5 * (kH / kNz)), 1e-12);
    // Both plates' mid-surfaces meet at the corner line, so there is no
    // through-thickness offset to carry and the weight is the mid-surface.
    expectNear("and the through-thickness weight is the mid-surface", tied.worstJunctionWeight, 0.5,
               1e-12);
    // The free edge is unchanged -- no element grew a neighbour -- and all but the
    // twelve half-metre segments touching a cut plane is now joined.
    expectNear("the free edge is what it was", tied.freeEdgeLength, cut.freeEdgeLength, 1e-12);
    expectNear("and so is the junction census", tied.junctionEdges, cut.junctionEdges, 1e-12);
    expectNear("of which all but the segments against a cut plane is tied", tied.tiedEdges,
               8 * kL - 12 * (kL / kNx), 1e-9);
    expectNear("where the untied mesh joins none of it", cut.tiedEdges, 0.0, 0.0);
    // And the twelve segments the section itself leaves open are exactly the ones a
    // line tie on its two planes would close. A chain applies them at its *interior*
    // planes only, which is why they are reported per plane rather than summed.
    //
    // **Four against the aft plane and eight against the forward one, and the
    // asymmetry is real rather than an accounting slip.** Of the two segments meeting
    // at a corner on a plane, whether the side plate's is already joined depends on
    // which of its two sub-quads §5 chose as the master face for the slave one station
    // in -- and that is always the lower-x one, so at the aft plane the side plate's
    // segment is already a master's and at the forward plane it is not.
    expectNear("what a line tie on the aft plane would add", tied.planeTiedEdgesAft,
               4 * (kL / kNx), 1e-9);
    expectNear("and on the forward plane", tied.planeTiedEdgesForward, 8 * (kL / kNx), 1e-9);
    expectNear("with nothing needing both", tied.planeTiedEdgesBoth, 0.0, 0.0);
    expectNear("which together is the whole junction census", tied.junctionEdges,
               tied.tiedEdges + tied.planeTiedEdgesAft + tied.planeTiedEdgesForward, 1e-9);
}

// **What an in-plane line tie is, checked on the constraint itself.**
//
// The two identities a tie has to have, on the four masters of a line rather than the
// eight of a face: the weights are a partition of unity, so a rigid translation is
// reproduced exactly, and they interpolate the slave's own point, so a rigid rotation
// is too. Checked here rather than through a solve because a tie wrong in either is
// wrong in a way an energy comparison reports as a small stiffness change.
//
// **And the property that makes it a line tie at all**: every degree of freedom in it,
// slave and master alike, lies on the cut plane it belongs to. That is the whole of
// why it survives the plane being an interface, and it is one loop.
void testPlaneTieLiesInThePlane() {
    std::printf("\n--- section: the in-plane line tie is in the plane ---\n");
    const StructuralMesh structure = makeBox(kT, kT);
    const section::Section tied = section::buildSection(structure, boxParams());
    expectTrue("the box built some line ties", !tied.planeTies.empty());

    // Which mesh nodes are on which cut plane, from the mesher's own answer rather
    // than from a position test -- `aftNodes` and `forwardNodes` are chosen on the
    // mid-surface, which is the distinction §6 note 1 cost a day.
    std::vector<int> plane(tied.nodeCount(), 0);
    for (std::uint32_t n : tied.aftNodes) plane[n] = 1;
    for (std::uint32_t n : tied.forwardNodes) plane[n] = 2;

    double worstUnity = 0, worstPoint = 0;
    int offPlane = 0;
    for (const section::PlaneTie& tie : tied.planeTies) {
        double sum = 0;
        Vec3 interpolated{0, 0, 0};
        const std::size_t slave = tie.mpc.slave / 3;
        if (plane[slave] != tie.plane) ++offPlane;
        for (std::size_t a = 0; a < tie.mpc.master.size(); ++a) {
            const std::size_t node = tie.mpc.master[a] / 3;
            if (plane[node] != tie.plane) ++offPlane;
            // The masters have to carry the same axis as the slave: a tie that
            // coupled x to y would solve and would be a different structure.
            if (tie.mpc.master[a] % 3 != tie.mpc.slave % 3) ++offPlane;
            sum += tie.mpc.weight[a];
            interpolated += Vec3{tied.mesh.position[node * 3], tied.mesh.position[node * 3 + 1],
                                 tied.mesh.position[node * 3 + 2]} *
                            tie.mpc.weight[a];
        }
        worstUnity = std::max(worstUnity, std::abs(sum - 1.0));
        const Vec3 at{tied.mesh.position[slave * 3], tied.mesh.position[slave * 3 + 1],
                      tied.mesh.position[slave * 3 + 2]};
        worstPoint = std::max(worstPoint, length(interpolated - at));
    }
    std::printf("     line tie: |sum w - 1| %.2e, |sum w X - X_slave| %.2e m, %d DOF off plane\n",
                worstUnity, worstPoint, offPlane);
    // Every degree of freedom on the plane is not a tolerance, it is the construction.
    expectEqual("every degree of freedom of a line tie is on its own cut plane", offPlane, 0);
    // A partition of unity to rounding: the shape functions are `0.5 (1 -/+ s)` times
    // `(1 - w)` and `w`, and both pairs sum to one exactly in exact arithmetic.
    expectNear("the weights are a partition of unity", worstUnity, 0.0, 1e-15);
    // On this box the masters *interpolate* the slave rather than extrapolating to it:
    // the corner nodes coincide, so the tie is exact to rounding. On a real hull the
    // line drops whatever the slave's normal leans along the ship, which is
    // `worstPlaneTieGap` and is measured on the ferry rather than here.
    expectNear("and they interpolate the slave's own point", worstPoint, 0.0, 1e-15);
    expectNear("which the mesher agrees it did", tied.worstPlaneTieSlip, 0.0, 1e-15);
    // The vacuity guard: a tie of one master to itself would pass both identities
    // above. There are four distinct masters, and they are four distinct nodes.
    std::size_t distinct = 0;
    for (const section::PlaneTie& tie : tied.planeTies) {
        std::vector<std::uint32_t> nodes;
        for (std::uint32_t d : tie.mpc.master) nodes.push_back(d / 3);
        std::sort(nodes.begin(), nodes.end());
        nodes.erase(std::unique(nodes.begin(), nodes.end()), nodes.end());
        distinct = std::max(distinct, nodes.size());
        expectTrue("no line tie names its own slave as a master",
                   std::find(nodes.begin(), nodes.end(), tie.mpc.slave / 3) == nodes.end());
    }
    expectEqualCount("a line tie has four masters: two mid-surface nodes, both extruded",
                     distinct, 4u);
}

// --- 1c. What the tie is, before what it does ------------------------------------
//
// Two identities, checked on the constraint set itself rather than through a solve,
// because a tie that is wrong in either is wrong in a way an energy comparison
// would report as a small stiffness change.

// --- 1b. Collapsed panels: the whole ship, not just its parallel middle body ------
//
// **This is the defect that kept the mesher off three quarters of the reference
// ferry's two-bay windows.** A
// degenerate plate panel -- a quad with two coincident corners -- extrudes to a
// triangular prism, and `solidshell::smallestJacobian` samples the *corners*, where
// such an element's determinant is exactly zero. `buildSection` reported "an element
// came out inverted", `applyBeamLoad` and `reduction::Substructure` refused it, and
// nothing forward of x = 19.2 m or aft of x = -26.4 m could be modelled. Nothing was
// ever inverted: over all 49 two-bay windows of that ship there is **not one negative
// Jacobian**, and the worst determinant over the *sound* elements is 2.7e-5 at the
// stem as well as amidships.
//
// The box is where every quantity has a closed form, so the claim is made here first
// and on the ferry below. Triangulating every panel doubles the element count and
// changes nothing that has an exact answer.
void testCollapsedPanelsMeshAsWedges() {
    std::printf("\n--- section: a panel with two coincident corners is a wedge, not a fold ---\n");
    const StructuralMesh quads = makeBox(kT, kT);
    const StructuralMesh triangles = triangulated(quads);
    const StructuralMaterial material = ah36Steel();
    const double youngs = material.youngsModulus;
    const double shear = youngs / (2 * (1 + material.poissonRatio));
    const double bredt = shear * kBoxTorsionConstant;

    const section::Section plain = section::buildSection(quads, openBoxParams());
    const section::Section wedges = section::buildSection(triangles, openBoxParams());

    // The mesher builds it, says what it built, and does not call it inverted.
    expectEqual("no quad panel is collapsed", plain.collapsedElements, 0);
    expectEqualCount("triangulating doubles the elements", wedges.elementCount(),
                     2 * plain.elementCount());
    expectEqual("and every one of them is a collapsed hexahedron", wedges.collapsedElements,
                static_cast<int>(wedges.elementCount()));
    expectEqual("none of which is inverted", wedges.invertedElements, 0);
    expectTrue("so the section's worst determinant is positive", wedges.worstJacobian > 0);
    // The guard that says the old rule really would have refused this. Without it
    // the test passes on a mesher that never had the problem.
    {
        double worstNodal = 1e300;
        for (std::size_t e = 0; e < wedges.elementCount(); ++e) {
            double nodes[solidshell::kDof];
            wedges.mesh.gather(e, wedges.mesh.position, nodes);
            worstNodal = std::min(worstNodal, solidshell::smallestJacobian(nodes));
        }
        std::printf("     worst determinant: %.3e at the corners, %.3e at the Gauss points\n",
                    worstNodal, wedges.worstJacobian);
        expectNear("the corner-sampling rule sees exactly zero on this mesh", worstNodal, 0.0, 0.0);
    }

    // --- The steel is all there, which "drop the element instead" would lose -------
    expectNear("the wedges cover the same mid-surface", wedges.area, plain.area, 1e-12);
    expectNear("and carry the same plating", wedges.plateMass, plain.plateMass, 1e-9);
    expectNear("which is the closed form", wedges.area, 2 * (kB + kH) * kL, 1e-12);

    // --- And they integrate ---------------------------------------------------------
    section::BeamLoad axial;
    axial.strain = 1e-6;
    section::BeamLoad bending;
    bending.curvature = 1e-6;
    bending.reference = kH / 2;
    const section::BeamResponse stretched = section::applyBeamLoad(wedges, material, axial);
    const section::BeamResponse bent = section::applyBeamLoad(wedges, material, bending);
    expectTrue("a section of wedges takes an axial load: " + stretched.problem, stretched.ok);
    expectTrue("and a curvature: " + bent.problem, bent.ok);
    // `EA` is a closed form and a wedge that passes the patch test reproduces it
    // exactly -- a mis-integrated one does not, because `EA` is the *integral* of
    // sigma over the cut and a wrong Jacobian is a wrong area.
    expectNear("EA is the closed form, to the same exactness the quad mesh reaches",
               stretched.axialStiffness, youngs * kBoxArea, kAxialTolerance * youngs * kBoxArea);
    // **`EI` is where a wedge costs something, and the cost is a discretisation
    // error rather than the formulation.** Two collapsed hexes are a different and
    // stiffer approximation than the one quad they replace -- the same reason a
    // constant-strain triangle is stiffer than a bilinear quad -- and at one element
    // per panel, with *every* panel triangulated, that is 13.7% of `EI`. What makes
    // it a mesh error and not a broken element is that it converges away, so that is
    // what is asserted rather than the number.
    const double second = bent.bendingStiffness / (youngs * kBoxSecondMoment);
    std::printf("     wedge mesh: EA/closed form %.10f, EI/closed form %.6f\n",
                stretched.axialStiffness / (youngs * kBoxArea), second);
    {
        double previous = std::abs(second - 1.0);
        std::printf("     EI error under refinement: %.3e", previous);
        bool falling = true;
        double finest = previous;
        for (int subdivision = 2; subdivision <= 3; ++subdivision) {
            const section::Section refined =
                section::buildSection(triangles, openBoxParams(subdivision));
            const section::BeamResponse fine =
                section::applyBeamLoad(refined, material, bending);
            expectTrue("the refined wedge mesh solves", fine.ok);
            const double error =
                std::abs(fine.bendingStiffness / (youngs * kBoxSecondMoment) - 1.0);
            std::printf(" %.3e", error);
            if (!(error < previous)) falling = false;
            previous = error;
            finest = error;
        }
        std::printf("\n");
        // 0.137 -> 0.0068 -> 0.0022: it falls at every step and by twenty fold on the
        // first, which is what a discretisation error does and a wrong element does
        // not. Asserted as *convergence*, because a fixed tolerance on the coarse
        // value would pass on an element that was simply wrong by a constant.
        expectTrue("the wedge mesh's EI error falls at every refinement", falling);
        expectTrue("and is under a per cent one refinement in", finest < 0.01);
        expectTrue("having started above ten per cent, so there was something to converge",
                   std::abs(second - 1.0) > 0.10);
    }

    // --- The one that says the wedges carry shear ------------------------------------
    //
    // `EA` and `EI` are prescribed by the cut planes and would come out right on a
    // mesh of wedges that touched nothing (§2). Torsion is not: Bredt needs the cell
    // closed, and closing it means shear flowing *through* the collapsed elements.
    const section::Section tiedWedges = section::buildSection(triangles, boxParams());
    const section::Section tiedQuads = section::buildSection(quads, boxParams());
    const section::TorsionResponse wedgeTwist =
        section::applyTwist(tiedWedges, material, 1e-6, kH / 2);
    const section::TorsionResponse quadTwist =
        section::applyTwist(tiedQuads, material, 1e-6, kH / 2);
    const section::TorsionResponse openTwist = section::applyTwist(wedges, material, 1e-6, kH / 2);
    expectTrue("all three twists ran",
               wedgeTwist.ok && quadTwist.ok && openTwist.ok);
    std::printf("     GJ/Bredt: wedges tied %.4f, quads tied %.4f, wedges untied %.4f\n",
                wedgeTwist.torsionalStiffness / bredt, quadTwist.torsionalStiffness / bredt,
                openTwist.torsionalStiffness / bredt);
    expectEqual("the tied wedge mesh is one piece", tiedWedges.components, 1);
    // Bredt plus the open-section term, as `testJunctionTieClosesTheCell` explains,
    // and the wedge mesh is over-stiff at one element per panel by the same 10-14%
    // it is in `EI`. The band is round the sum and is wide enough for that; what it
    // is not wide enough for is an *open* cell, which is a factor of ten below.
    expectTrue("and carries Bredt's shear flow",
               wedgeTwist.torsionalStiffness > 0.9 * bredt &&
                   wedgeTwist.torsionalStiffness < 1.3 * bredt);
    // Vacuity: an open cell is an order of magnitude softer, so this is a
    // measurement of something. The same guard `testJunctionTieClosesTheCell` uses.
    expectTrue("where the untied wedge mesh carries a fraction of it",
               openTwist.torsionalStiffness < 0.2 * bredt);
    // And the two discretisations converge on each other: they differ by 9.4% at one
    // element per panel and 0.03% at two, which is the same statement `EI` makes.
    const section::Section fineWedges = section::buildSection(triangles, boxParams(2));
    const section::Section fineQuads = section::buildSection(quads, boxParams(2));
    const section::TorsionResponse fineWedgeTwist =
        section::applyTwist(fineWedges, material, 1e-6, kH / 2);
    const section::TorsionResponse fineQuadTwist =
        section::applyTwist(fineQuads, material, 1e-6, kH / 2);
    expectTrue("the refined twists ran", fineWedgeTwist.ok && fineQuadTwist.ok);
    const double coarseGap =
        std::abs(wedgeTwist.torsionalStiffness / quadTwist.torsionalStiffness - 1.0);
    const double fineGap =
        std::abs(fineWedgeTwist.torsionalStiffness / fineQuadTwist.torsionalStiffness - 1.0);
    std::printf("     wedges against quads in GJ: %.4f at one element per panel, %.4f at two\n",
                coarseGap, fineGap);
    expectTrue("the two discretisations converge on the same torsion", fineGap < 0.1 * coarseGap);
    expectTrue("having started apart, so the comparison had something in it", coarseGap > 0.05);

    // --- Nothing floats, and the reduction takes it -----------------------------------
    const reduction::Substructure substructure(tiedWedges.mesh, tiedWedges.material,
                                               tiedWedges.interfaceNodes, tiedWedges.attachment);
    expectTrue("a substructure of collapsed elements is ready", substructure.ready());
    expectEqual("with nothing floating", tiedWedges.floatingComponents, 0);
    // The condensed mass has to stay positive at the apex of every wedge, where two
    // of the element's own nodes land on one mesh node. Nothing else here would
    // notice a lumping that dropped one of them.
    expectNear("and the same total mass as the quad mesh", substructure.totalMass(),
               tiedQuads.mass(), 1e-9 * tiedQuads.mass());
}

// **A wedge's apex is named twice by its own sub-quad, and the nodal averages must
// count it once.** `nodeThickness` and `nodeNormal` are area-weighted over the
// sub-quads that reach a node; a collapsed sub-quad lists its apex in two of its four
// corners, so counting per corner gives that one panel double weight at exactly the
// node where several panels meet. Mutation testing found this: removing the
// deduplication survived everything else.
//
// Two coplanar panels sharing an edge, one of them a triangle, with *different*
// thicknesses -- so the apex node's mean thickness is a closed form that the two
// countings disagree about.
void testWedgeApexIsCountedOnce() {
    std::printf("\n--- section: the nodal average counts a wedge's apex once ---\n");
    const double square = 0.010, triangle = 0.020;
    StructuralMesh mesh;
    mesh.materials = {ah36Steel()};
    mesh.frameSpacing = 1.0;
    const auto panel = [&](const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d, double t) {
        PlatePanel p;
        p.corner[0] = a;
        p.corner[1] = b;
        p.corner[2] = c;
        p.corner[3] = d;
        p.thickness = t;
        p.material = 0;
        p.role = PanelRole::Shell;
        p.source = 0;
        mesh.panels.push_back(p);
    };
    // A unit square, and to starboard of it a triangle whose collapsed corner sits on
    // the square's own corner at (1, 0, 0).
    panel({0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, square);
    panel({1, 0, 0}, {2, 0, 0}, {1, 1, 0}, {1, 0, 0}, triangle);
    mesh.frameStations = {0.0, 1.0, 2.0};

    section::SectionParams params;
    params.xFrom = 0.0;
    params.xTo = 2.0;
    params.subdivision = 1;
    params.members = false;
    const section::Section flat = section::buildSection(mesh, params);
    expectEqualCount("two panels, two elements", flat.elementCount(), 2u);
    expectEqual("one of which is a wedge", flat.collapsedElements, 1);
    expectEqual("and neither is inverted", flat.invertedElements, 0);

    // Areas: the square is 1 m^2, the triangle is 1/2. The shared node at (1, 0, 0)
    // is reached by both, so its area-weighted thickness is
    //     (1 * 10 mm + 0.5 * 20 mm) / 1.5 = 13.333 mm,
    // where counting the triangle twice would give (1 * 10 + 1.0 * 20) / 2 = 15 mm.
    const double correct = (1.0 * square + 0.5 * triangle) / 1.5;
    const double doubled = (1.0 * square + 1.0 * triangle) / 2.0;
    double apex = -1;
    for (std::size_t n = 0; n < flat.nodeCount(); ++n)
        if (std::abs(flat.mesh.position[n * 3] - 1.0) < 1e-9 &&
            std::abs(flat.mesh.position[n * 3 + 1]) < 1e-9)
            apex = flat.nodeThickness[n];
    std::printf("     apex thickness %.6f m; counted once %.6f, counted twice %.6f\n", apex,
                correct, doubled);
    expectTrue("the apex node was found", apex > 0);
    expectNear("and carries the mean of the two panels it belongs to", apex, correct, 1e-12);
    // The guard: the two answers have to differ, or the assertion is about nothing.
    expectTrue("the two countings really do disagree", std::abs(correct - doubled) > 1e-4);

    // And the mass follows, because the plate mass is the same area-weighted
    // thickness integrated: 1 m^2 at a mean of the two nodes on the seam plus the
    // triangle's own half. Asserted through the total, which is what a caller sees.
    expectNear("the section covers both panels", flat.area, 1.5, 1e-12);

    // **And the other half of the classification, which nothing else here exercises.**
    // Every section in this file has `invertedElements == 0`, so a counter that never
    // fired would pass all of them -- mutation testing said exactly that.
    //
    // The case the check exists for: **plating thick against the panel it is meshed
    // on, over a fold.** Two panels 0.10 m long meeting at 40 degrees -- inside
    // `foldLimit`, so they weld into one surface -- carrying 1.00 m of thickness. The
    // node on the seam is extruded half a thickness along the bisector, which is five
    // times the panel's own length, so the inner face turns through itself. Nothing
    // here is collapsed; the element is simply inside out.
    const auto foldedPlate = [&](double panelLength, double thickness) {
        const double angle = 40.0 * 3.14159265358979323846 / 180.0;
        StructuralMesh bent;
        bent.materials = {ah36Steel()};
        bent.frameSpacing = panelLength;
        const Vec3 far{panelLength + panelLength * std::cos(angle), 0,
                       panelLength * std::sin(angle)};
        const auto add = [&](const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d) {
            PlatePanel p;
            p.corner[0] = a;
            p.corner[1] = b;
            p.corner[2] = c;
            p.corner[3] = d;
            p.thickness = thickness;
            p.material = 0;
            p.role = PanelRole::Shell;
            p.source = 0;
            bent.panels.push_back(p);
        };
        add({0, 0, 0}, {panelLength, 0, 0}, {panelLength, 1, 0}, {0, 1, 0});
        add({panelLength, 0, 0}, far, {far.x, 1, far.z}, {panelLength, 1, 0});
        section::SectionParams p;
        p.xFrom = 0.0;
        p.xTo = far.x;
        p.subdivision = 1;
        p.members = false;
        return section::buildSection(bent, p);
    };

    const section::Section inverted = foldedPlate(0.10, 1.00);
    std::printf("     1.00 m of plating over a 0.10 m fold: %d inverted, %d collapsed,"
                " %d surfaces, worst Gauss determinant %.3e\n",
                inverted.invertedElements, inverted.collapsedElements, inverted.surfaces,
                inverted.worstJacobian);
    expectEqual("the two panels weld into one surface", inverted.surfaces, 1);
    expectEqual("both elements come out inverted", inverted.invertedElements, 2);
    expectEqual("and neither is a wedge", inverted.collapsedElements, 0);
    expectTrue("so the section's worst determinant is not positive",
               !(inverted.worstJacobian > 0));
    expectTrue("and the mesher says inverted rather than collapsed", [&] {
        for (const std::string& problem : inverted.problems)
            if (problem.find("came out inverted") != std::string::npos) return true;
        return false;
    }());
    const section::BeamResponse refused =
        section::applyBeamLoad(inverted, ah36Steel(), section::BeamLoad{1e-6, 0, 0});
    expectTrue("and loading it is refused rather than solved", !refused.ok);

    // **The corner rule earns its place here too.** At 0.15 m the same fold leaves
    // every *quadrature* point positive and still turns two corners through
    // themselves, so a check that looked only at the Gauss points would accept it.
    const section::Section corners = foldedPlate(0.15, 1.00);
    std::printf("     the same over a 0.15 m fold: %d inverted, worst Gauss determinant %.3e\n",
                corners.invertedElements, corners.worstJacobian);
    expectTrue("the quadrature alone would accept the shallower fold",
               corners.worstJacobian > 0);
    expectEqual("but its corners are inside out", corners.invertedElements, 2);
    expectTrue("so it is refused as well",
               !section::applyBeamLoad(corners, ah36Steel(), section::BeamLoad{1e-6, 0, 0}).ok);

    // The vacuity guard: the same fold on ordinary plating is fine, so this is about
    // the thickness against the panel and not about folds.
    const section::Section sound = foldedPlate(0.10, 0.10);
    expectEqual("the same fold on plating a tenth as thick is sound",
                sound.invertedElements, 0);
    expectTrue("and solves", section::applyBeamLoad(sound, ah36Steel(),
                                                    section::BeamLoad{1e-6, 0, 0}).ok);
}

void testJunctionTieIsAPartitionOfUnity() {
    const section::Section tied = section::buildSection(makeBox(kT, kT), boxParams());
    expectTrue("the box built some ties", !tied.attachment.constrained.empty());
    expectEqualCount("three constraints per extruded node, eight masters each",
                     tied.attachment.constrained.size(),
                     static_cast<std::size_t>(tied.junctionTies) * 2 * 3);

    double worstSum = 0, worstPoint = 0, worstAxis = 0;
    for (const solidshell::Mpc& mpc : tied.attachment.constrained) {
        expectEqualCount("a junction tie has eight masters", mpc.master.size(), 8u);
        double sum = 0;
        Vec3 interpolated{0, 0, 0};
        for (std::size_t a = 0; a < mpc.master.size(); ++a) {
            sum += mpc.weight[a];
            // Every master must be the same axis as the slave: the constraint is
            // between like components, and a tie that mixed them would rotate the
            // plating it joined.
            worstAxis = std::max<double>(worstAxis,
                                         std::abs(static_cast<int>(mpc.master[a] % 3) -
                                                  static_cast<int>(mpc.slave % 3)));
            const std::size_t node = mpc.master[a] / 3;
            for (int k = 0; k < 3; ++k)
                interpolated[k] += mpc.weight[a] * tied.mesh.position[node * 3 + static_cast<std::size_t>(k)];
        }
        worstSum = std::max(worstSum, std::abs(sum - 1.0));
        const std::size_t slaveNode = mpc.slave / 3;
        worstPoint = std::max(
            worstPoint,
            length(interpolated - Vec3{tied.mesh.position[slaveNode * 3],
                                       tied.mesh.position[slaveNode * 3 + 1],
                                       tied.mesh.position[slaveNode * 3 + 2]}));
    }
    std::printf("     junction tie: |sum w - 1| %.2e, |sum w X - X_slave| %.2e m\n", worstSum,
                worstPoint);
    expectNear("a master and its slave are the same axis", worstAxis, 0.0, 0.0);
    // Sum to one: the tie reproduces a rigid **translation** exactly. 1.1e-16
    // measured, so the tolerance is rounding and not a band drawn round zero.
    expectTrue("the weights are a partition of unity", worstSum < 4e-16);
    // And the weighted master positions are the slave's own position, so the tie
    // reproduces a rigid **rotation** exactly too: `sum w R X = R sum w X = R X`.
    // 8.9e-16 m measured on a box two metres across.
    expectTrue("and they interpolate the slave's own position", worstPoint < 1e-14);
}

// --- 1b. The weld is a distance, the winding is not a promise --------------------
//
// Both of these exist because mutation testing found nothing else that could see
// them. A weld radius that is not a radius, and a bucket probe that does not reach
// as far as the tolerance does, both leave a crack down the middle of a mesh that
// is invisible in every aggregate -- which is the mutant `zone.cpp`'s own header
// records surviving everything it had.

void testWeldIsADistance() {
    const StructuralMesh structure = makeBox(kT, kT);
    // Untied throughout: the junction tie joins the four corners into one component
    // whatever the welder does, so a component count taken with it on would be
    // measuring the tie and reporting it as the weld.
    const section::Section reference = section::buildSection(structure, openBoxParams());

    // A displacement of 0.3 um is inside the 1 um weld tolerance, so the panel must
    // stay joined -- **and it is chosen to straddle a bucket boundary**. The welder
    // buckets on a 2 um grid and the bay seam sits at x = 0.5, an exact multiple of
    // it, so moving one panel back by 0.3 um puts the two copies of that corner in
    // *adjacent* cells. A probe that only looks in its own cell welds every mesh
    // whose duplicates are bit-identical and cracks this one.
    for (double delta : {-3e-7, 3e-7}) {
        StructuralMesh nudged = structure;
        for (int c = 0; c < 4; ++c)
            if (std::abs(nudged.panels[0].corner[c].x - 0.5) < 1e-12)
                nudged.panels[0].corner[c].x = 0.5 + delta;
        const section::Section joined = section::buildSection(nudged, openBoxParams());
        expectEqual("a sub-tolerance nudge across a bucket boundary still welds",
                    joined.components, reference.components);
        expectEqualCount("and adds no nodes", joined.nodeCount(), reference.nodeCount());
    }

    // Past the tolerance the panel must come away. Three microns is the one that
    // matters: it is **inside the bucket probe's reach** and outside the tolerance,
    // so only the distance test can refuse it. A weld that compared the squared
    // distance against the tolerance rather than against its square would have a
    // radius of a millimetre and would still pass at ten microns, because the probe
    // would not have looked that far -- one bug masking another.
    for (double delta : {3e-6, 1e-5}) {
        StructuralMesh torn = structure;
        for (int c = 0; c < 4; ++c)
            if (std::abs(torn.panels[0].corner[c].x - 0.5) < 1e-12)
                torn.panels[0].corner[c].x = 0.5 + delta;
        const section::Section split = section::buildSection(torn, openBoxParams());
        expectTrue("a nudge past the tolerance does not weld",
                   split.components > reference.components);
        expectTrue("and shows up as extra nodes", split.nodeCount() > reference.nodeCount());
    }
}

void testReversedWindingIsOriented() {
    // `makeStructuralMesh` mirrors the starboard side, so panels on one surface do
    // not all wind the same way round. The surface walk has to orient them against
    // one another; without it the nodal normals cancel, the extrusion collapses and
    // the elements are degenerate -- which is a mesh that computes numbers.
    const section::Section section = section::buildSection(makeBox(kT, kT, true), openBoxParams());
    const section::Section reference = section::buildSection(makeBox(kT, kT), openBoxParams());

    expectEqual("alternating winding is still four surfaces", section.surfaces, 4);
    expectTrue("no element is inverted", section.worstJacobian > 0);
    // The tell-tale: an unoriented mesh averages `+n` against `-n` at every seam, so
    // the nodal normals scatter instead of agreeing. On a flat plate they agree
    // exactly, which is a stronger statement than the Jacobian's.
    expectNear("and the nodal normals agree exactly, as they do on flat plating",
               section.worstNormalSpread, 0.0, 1e-12);
    expectNear("the mid-surface area is unchanged", section.area, reference.area, 1e-12);
    expectEqualCount("and so is the node count", section.nodeCount(), reference.nodeCount());
    const StructuralMaterial material = ah36Steel();
    section::BeamLoad axial;
    axial.strain = 1e-6;
    expectNear("and so is EA",
               section::applyBeamLoad(section, material, axial).axialStiffness /
                   (material.youngsModulus * kBoxArea),
               1.0, kAxialTolerance);
}

// --- 2. Section properties against closed forms ----------------------------------

void testBoxSectionProperties() {
    const StructuralMesh structure = makeBox(kT, kT);
    // The closed forms below are beam theory, which is exact for this problem: a
    // linear `sigma_xx` field is self-equilibrated and compatible, so it is the
    // three-dimensional answer too. They are asserted on the **untied** box because
    // that mesh reproduces it to rounding; what the tie costs the same numbers is
    // `testJunctionTieClosesTheCell` and it is a consistency error that converges.
    const section::Section box = section::buildSection(structure, openBoxParams());
    const StructuralMaterial material = ah36Steel();
    const double youngs = material.youngsModulus;

    section::BeamLoad axial;
    axial.strain = 1e-6;
    axial.reference = kH / 2;
    const section::BeamResponse stretched = section::applyBeamLoad(box, material, axial);
    expectTrue("the axial solve ran: " + stretched.problem, stretched.ok);
    // 6e-13 measured. The field is linear in x and the element reproduces a linear
    // field exactly, so this is an identity and not a convergence: asserting it at
    // anything looser would pass on a mesh that had lost a plate and was merely
    // close.
    expectNear("EA against 2(B+H)t E", stretched.axialStiffness / (youngs * kBoxArea), 1.0,
               kAxialTolerance);
    // A doubly symmetric box has its neutral axis at mid-height, so a pure axial
    // strain puts no moment about it.
    expectTrue("a pure axial strain carries no moment about mid-height",
               std::abs(stretched.bendingMoment) < 1e-7);
    // The restraints that remove the rigid translations and the roll are statically
    // determinate, so they carry exactly nothing. This is the only assertion in the
    // file that would notice `applyBeamLoad` holding a degree of freedom the exact
    // solution wants to move.
    // Against an axial force of 1.2e4 N, so 1e-7 N is 1e-11 of the load and not a
    // band drawn round zero.
    expectTrue("the rigid-body restraints carry no force", stretched.restraintReaction < 1e-7);
    expectTrue("the free degrees of freedom are in equilibrium", stretched.residual < 1e-7);

    section::BeamLoad bending;
    bending.curvature = 1e-6;
    bending.reference = kH / 2;
    const section::BeamResponse bent = section::applyBeamLoad(box, material, bending);
    expectTrue("the bending solve ran: " + bent.problem, bent.ok);
    expectNear("EI against the closed form", bent.bendingStiffness / (youngs * kBoxSecondMoment),
               1.0, kBendingTolerance);
    // The guard that makes the line above mean something: the mid-surface second
    // moment -- the one a hand calculation gives, without the plates' own `b t^3/12`
    // -- differs from the exact one by 2.9e-5, which is six times the tolerance
    // asserted. So the test can tell the two apart, and it picked the right one.
    expectTrue("the tolerance is finer than the plates' own thickness term",
               std::abs(kBoxSecondMoment / kBoxSecondMomentMidSurface - 1.0) > kBendingTolerance);
    expectTrue("pure bending about the true neutral axis carries no axial force",
               std::abs(bent.axialForce) < 1e-7);

    // ... and the guard for *that*: move the reference 100 mm and the axial force is
    // not small at all, so "no axial force" is a statement about the neutral axis
    // rather than about the loading being too weak to see.
    section::BeamLoad offset = bending;
    offset.reference = kH / 2 + 0.1;
    const section::BeamResponse leaning = section::applyBeamLoad(box, material, offset);
    expectTrue("a curvature about the wrong axis does carry one",
               std::abs(leaning.axialForce) > 1e3);
}

// --- 3. The junctions: invisible to bending, loud in torsion ---------------------
//
// **The three meshes this file can build over the same panels, and there is exactly
// one question that separates them.** `EA` and `EI` are prescribed by the two cut
// planes -- every longitudinally continuous strip carries `sigma = E eps` whatever
// it is attached to -- so a mesh that joins nothing scores *best* on them. Torsion
// is the one that needs the cell closed, and it is the whole of §5's case.

void testJunctionTieClosesTheCell() {
    const StructuralMesh structure = makeBox(kT, kT);
    const StructuralMaterial material = ah36Steel();
    const double youngs = material.youngsModulus;
    const double shear = youngs / (2 * (1 + material.poissonRatio));
    const double bredt = shear * kBoxTorsionConstant;

    const section::Section cut = section::buildSection(structure, openBoxParams());
    const section::Section welded = section::buildSection(structure, weldedBoxParams());
    const section::Section tied = section::buildSection(structure, boxParams());

    section::BeamLoad axial;
    axial.strain = 1e-6;
    section::BeamLoad bending;
    bending.curvature = 1e-6;
    bending.reference = kH / 2;

    struct Row {
        const char* label;
        const section::Section& section;
        double area, second, torsion;
    };
    Row rows[] = {{"corners cut", cut, 0, 0, 0},
                  {"corners welded", welded, 0, 0, 0},
                  {"corners tied", tied, 0, 0, 0}};
    for (Row& row : rows) {
        const section::BeamResponse stretched =
            section::applyBeamLoad(row.section, material, axial);
        const section::BeamResponse bent = section::applyBeamLoad(row.section, material, bending);
        const section::TorsionResponse twisted =
            section::applyTwist(row.section, material, 1e-6, kH / 2);
        expectTrue(std::string(row.label) + " solves: " + stretched.problem + bent.problem +
                       twisted.problem,
                   stretched.ok && bent.ok && twisted.ok);
        row.area = stretched.axialStiffness / (youngs * kBoxArea);
        row.second = bent.bendingStiffness / (youngs * kBoxSecondMoment);
        row.torsion = twisted.torsionalStiffness / bredt;
        std::printf("     %-15s %d component(s), band %4zu: EA %.8f  EI %.8f  GJ/Bredt %.4f\n",
                    row.label, row.section.components, row.section.halfBandwidth, row.area,
                    row.second, row.torsion);
    }
    const Row& open = rows[0];
    const Row& joined = rows[1];
    const Row& tie = rows[2];

    // --- What each one costs the quantities the two cut planes prescribe ---
    //
    // Welding a right-angled corner extrudes the shared node along the mean normal
    // and thins the plating towards the corner by `cos 45`, so the section simply
    // loses steel. A tie moves no node, so it loses none: `EA` comes out at the
    // untied mesh's own exactness, which is machine precision.
    expectNear("EA is exact on the section whose corners are cut", open.area, 1.0,
               kAxialTolerance);
    expectNear("and on the one whose corners are tied, because a tie moves no steel", tie.area,
               1.0, kAxialTolerance);
    expectTrue("where welding them loses at least 5% of EA", joined.area < 0.95);
    expectTrue("and at least 5% of EI", joined.second < 0.95);

    // What the tie *does* cost `EI` is a consistency error at the butt corner, where
    // the tied node sits half a plate thickness past the end of the other plate's
    // mid-surface and the bilinear extrapolation of a quadratic transverse
    // contraction is out by `O(h^2)`. 2.16e-3 measured at one element per panel;
    // `testResolutionConvergence` is the evidence that it is a discretisation error
    // and not the formulation, because it falls by twenty-five fold over a
    // three-fold refinement.
    expectNear("the tie costs EI 0.22% at one element per panel", tie.second, 1.0021578, 1e-6);
    expectTrue("which is forty times less than welding costs it",
               (1.0 - joined.second) / (tie.second - 1.0) > 40.0);

    // --- And the one quantity that can tell them apart ---
    //
    // Bredt's shear flow needs the cell closed. It is not closed by proximity: the
    // untied mesh has its four corner lines coincident to a weld tolerance and
    // carries a twelfth of Bredt.
    expectTrue("the cut box carries less than a fifth of Bredt", open.torsion < 0.2);
    expectNear("the welded box carries Bredt's torsion", joined.torsion, 1.0, 0.05);
    // The tied box carries Bredt **plus** the open-section term, which is the
    // `sum s t^3 / 3` the cut box is made of: a closed cell does not stop having
    // plate torsion. The two add to 1.083 and the mesh gives 1.099 at one element
    // per panel, converging downwards on refinement. Asserted as a band round the
    // sum rather than round 1, because 1 is the number a section that had *lost*
    // its open-section term would give.
    expectNear("and the tied box carries Bredt plus the plating's own torsion",
               tie.torsion, 1.0 + open.torsion, 0.02);
    expectTrue("which is more than the welded box manages, because tying loses no steel",
               tie.torsion > joined.torsion);
    expectTrue("and thirteenfold what the cut one does", tie.torsion / open.torsion > 13.0);

    // **The tied section has to be in equilibrium, and only the tie's transpose puts
    // it there.** A junction eliminates its slave, so `K u` lands on a degree of
    // freedom the solved system does not have; until it is moved to the masters by
    // the transpose of the constraint, the masters' rows are short by whatever the
    // junction next to them carries. Nothing else here would see it -- `EA` comes
    // from the *interface* reaction, which is a long way from any tie.
    {
        const section::BeamResponse stretched = section::applyBeamLoad(tied, material, axial);
        const section::BeamResponse bent = section::applyBeamLoad(tied, material, bending);
        std::printf("     tied box equilibrium: residual %.2e N and %.2e N, restraint %.2e N, "
                    "against an axial force of %.3e N\n",
                    stretched.residual, bent.residual, stretched.restraintReaction,
                    stretched.axialForce);
        // 3.7e-9 and 5.6e-8 measured against an axial force of 1.24e4 N, so 1e-6 N is
        // 1e-10 of the load and not a band drawn round zero.
        expectTrue("the free degrees of freedom of the tied section are in equilibrium",
                   stretched.residual < 1e-6 && bent.residual < 1e-6);
        expectTrue("and its rigid-body restraints still carry nothing",
                   stretched.restraintReaction < 1e-5);
        expectTrue("against a load worth measuring against", stretched.axialForce > 1e4);
    }

    // **The overshoot limit is a real gate, not a formality.** A butt corner puts the
    // tied node 0.02 of a face outside it; a limit either side of that number turns
    // every tie in the box on or off, so the parameter is exercised rather than
    // merely present.
    for (const auto& [limit, want] : std::vector<std::pair<double, int>>{{0.019, 0}, {0.021, 60}}) {
        section::SectionParams bounded = boxParams();
        bounded.junctionOvershoot = limit;
        const section::Section gated = section::buildSection(structure, bounded);
        expectEqual("an overshoot limit of " + std::to_string(limit) + " leaves " +
                        std::to_string(want) + " ties",
                    gated.junctionTies, want);
        expectEqual("and the refusals are counted rather than dropped",
                    gated.junctionsOutsideFace, want == 0 ? 120 : 0);
        expectEqual("so the components follow", gated.components, want == 0 ? 4 : 1);
    }

    // The vacuity guard for the whole comparison: these are three meshes of the
    // *same* elements over the *same* panels. If they differed in size the torsion
    // ratio would be a statement about resolution.
    expectEqualCount("all three meshes have the same elements", tie.section.elementCount(),
                     open.section.elementCount());
    expectEqualCount("and the welded one too", joined.section.elementCount(),
                     open.section.elementCount());
    expectEqualCount("the tie adds no nodes either", tie.section.nodeCount(),
                     open.section.nodeCount());
    expectTrue("where welding removes some", joined.section.nodeCount() < open.section.nodeCount());
}

// --- 4. Thickness seams -----------------------------------------------------------

void testThicknessSeam() {
    const double inner = 0.010, outer = 0.016;
    const StructuralMesh structure = makeBox(inner, outer);
    const StructuralMaterial material = ah36Steel();
    const double youngs = material.youngsModulus;
    const double area = 2 * ((kB / 2) * inner + (kB / 2) * outer) + 2 * kH * inner;

    section::SectionParams tapered = boxParams();
    section::SectionParams split = boxParams();
    split.thicknessSeam = section::ThicknessSeam::Split;
    const section::Section withTaper = section::buildSection(structure, tapered);
    const section::Section withSplit = section::buildSection(structure, split);

    expectTrue("the seam produced tapered elements", withTaper.taperedElements > 0);
    expectEqual("splitting the seam leaves every element prismatic", withSplit.taperedElements, 0);
    // The node on the seam takes the area-weighted mean of the two strakes, 13 mm,
    // so the element between it and the 10 mm side has `dt/t = 3/11.5`. A closed
    // form, not a reading.
    expectNear("the worst taper is the half-step over the mean", withTaper.worstTaper,
               (0.5 * (inner + outer) - inner) / (0.5 * (inner + 0.5 * (inner + outer))), 1e-12);
    // And the geometry the spike's rule calls badly wrong really is present: 90 dt^2
    // is more than five, so this is not a test of a taper too small to matter.
    expectTrue("the taper is well past what the element is exact on",
               withTaper.taperStiffness > 5.0);
    // The price of splitting instead: the section comes apart along the strake seam
    // as well as at the corners.
    expectTrue("splitting costs connectivity", withSplit.components > withTaper.components);

    section::BeamLoad axial;
    axial.strain = 1e-6;
    const double taperArea =
        section::applyBeamLoad(withTaper, material, axial).axialStiffness / youngs;
    const double splitArea =
        section::applyBeamLoad(withSplit, material, axial).axialStiffness / youngs;
    expectNear("the tapered section still has the right area", taperArea / area, 1.0,
               kBendingTolerance);
    expectNear("the split one too", splitArea / area, 1.0, kAxialTolerance);
    // 4.5e-7 measured. Six hundred per cent of excess *plate bending* stiffness is
    // worth 4e-7 of a membrane quantity, which is the whole argument of §3.
    expectNear("and the taper costs the membrane answer nothing", taperArea / splitArea, 1.0,
               kBendingTolerance);
}

// --- 4b. A member run stops where the plate steps --------------------------------

void testMemberRunsStopAtAThicknessStep() {
    // One transverse frame across the bottom flange of the two-thickness box, so it
    // crosses the 10/16 mm seam with two stations of each thickness either side and
    // the seam's own averaged node in the middle.
    //
    // `constraint::addStiffener` takes one `plateThickness` for a whole run and turns
    // it into one tie weight per fibre; a weight is only right for the pair
    // separation it was computed against, so a run has to stop where the separation
    // changes. What that costs is one segment per seam -- the seam station is left in
    // a run of its own -- and it is the honest direction to fail in, because missing
    // steel shows up in a mass and a misplaced eccentricity does not.
    // It stops short of the corners on purpose. A member that reaches one picks up
    // the side plate's node pair as well, whose thickness direction is 90 degrees
    // from the flange's, and the mesher then refuses the whole member for having no
    // single eccentricity -- which is right, and is not what this test is about.
    StructuralMesh structure = makeBox(0.010, 0.016);
    {
        StructuralMember frame;
        frame.a = {kL / 2, -0.9 * kB / 2, 0.0};
        frame.b = {kL / 2, 0.9 * kB / 2, 0.0};
        frame.rise = {0, 0, 1};
        frame.profile = flatBar(0.200, 0.010);
        frame.attachedPlateThickness = 0.010;
        frame.role = MemberRole::Frame;
        structure.members.push_back(frame);
    }
    // Subdivision 2, so the frame crosses three stations of 10 mm, the seam's own
    // averaged 13 mm station, and three of 16 mm.
    section::SectionParams params = boxParams(2);
    params.members = true;
    const section::Section section = section::buildSection(structure, params);

    expectEqual("the frame was attached", section.membersAttached, 1);
    expectTrue("and its run was broken at the step", section.memberRunsSplitByThickness > 0);
    // Two segments survive either side of the seam and a flat bar is two Gauss
    // fibres, so eight. The seam station is left in a run of one and adds nothing:
    // that is the two segments the rule costs.
    expectEqualCount("two runs of three stations survive the step",
                     section.stiffening.fiberCount(), 8);

    // And the invariant that is the whole reason for the rule, on a mesh where the
    // step is 10 mm against 16 mm rather than a fraction of a millimetre: every
    // fibre sits at the offset it records, from its own pair's mid-surface.
    double worst = 0, largest = 0;
    for (const constraint::Fiber& fiber : section.stiffening.fiber)
        for (int end = 0; end < 2; ++end) {
            const constraint::Tie& tie = fiber.end[end];
            const auto at = [&](std::uint32_t node) {
                return Vec3{section.mesh.position[static_cast<std::size_t>(node) * 3],
                            section.mesh.position[static_cast<std::size_t>(node) * 3 + 1],
                            section.mesh.position[static_cast<std::size_t>(node) * 3 + 2]};
            };
            const Vec3 tied = constraint::tiedPoint(tie, section.mesh.position);
            const Vec3 middle = (at(tie.bottom) + at(tie.top)) * 0.5;
            worst = std::max(worst, std::abs(length(tied - middle) - std::abs(fiber.offset)));
            largest = std::max(largest, std::abs(fiber.offset));
        }
    expectTrue("every fibre sits at its own offset", worst < 1e-11);
    expectTrue("on offsets of tens of millimetres, so that is a statement", largest > 0.02);
}

// --- 5. Resolution ----------------------------------------------------------------

void testResolutionConvergence() {
    const StructuralMesh structure = makeBox(kT, kT);
    const StructuralMaterial material = ah36Steel();
    const double youngs = material.youngsModulus;

    section::BeamLoad axial;
    axial.strain = 1e-6;
    section::BeamLoad bending;
    bending.curvature = 1e-6;
    bending.reference = kH / 2;

    std::vector<double> torsion, tiedBendingError, tiedTorsion;
    std::size_t coarseElements = 0, fineElements = 0;
    for (int subdivision = 1; subdivision <= 3; ++subdivision) {
        const section::Section box =
            section::buildSection(structure, openBoxParams(subdivision));
        if (subdivision == 1) coarseElements = box.elementCount();
        if (subdivision == 3) fineElements = box.elementCount();
        expectNear("EA does not move with resolution, subdivision " + std::to_string(subdivision),
                   section::applyBeamLoad(box, material, axial).axialStiffness /
                       (youngs * kBoxArea),
                   1.0, kAxialTolerance);
        expectNear("nor does EI, subdivision " + std::to_string(subdivision),
                   section::applyBeamLoad(box, material, bending).bendingStiffness /
                       (youngs * kBoxSecondMoment),
                   1.0, kBendingTolerance);
        torsion.push_back(section::applyTwist(box, material, 1e-6, kH / 2).torsionalStiffness);

        // The same sweep with the junctions tied. `EA` stays exact at every
        // resolution -- the tie moves no steel -- while `EI` carries the butt
        // corner's extrapolation error, and **that** is what has to converge.
        const section::Section joined = section::buildSection(structure, boxParams(subdivision));
        expectNear("EA is exact whatever the resolution, tied, subdivision " +
                       std::to_string(subdivision),
                   section::applyBeamLoad(joined, material, axial).axialStiffness /
                       (youngs * kBoxArea),
                   1.0, kAxialTolerance);
        tiedBendingError.push_back(
            std::abs(section::applyBeamLoad(joined, material, bending).bendingStiffness /
                         (youngs * kBoxSecondMoment) -
                     1.0));
        tiedTorsion.push_back(section::applyTwist(joined, material, 1e-6, kH / 2).torsionalStiffness);
    }
    // The vacuity guard for the two assertions above: the mesh really did refine
    // ninefold, so "did not move" is a property of the quantity and not of the sweep.
    expectEqualCount("the sweep refined the mesh ninefold", fineElements, coarseElements * 9);
    // And the guard that the sweep can see a change at all: the one quantity here
    // with plate bending in it does move, monotonically downwards, because a
    // reduction can only stiffen and refining removes stiffness.
    expectTrue("torsion falls monotonically under refinement",
               torsion[1] < torsion[0] && torsion[2] < torsion[1]);
    expectTrue("by more than the tolerance the others are asserted at",
               (torsion[0] - torsion[2]) / torsion[0] > 1e-3);
    expectTrue("but by less than three per cent, so one element per panel is converged",
               (torsion[0] - torsion[2]) / torsion[0] < 0.03);

    std::printf("     tied box under refinement: EI error %.2e %.2e %.2e, GJ %.4e %.4e %.4e\n",
                tiedBendingError[0], tiedBendingError[1], tiedBendingError[2], tiedTorsion[0],
                tiedTorsion[1], tiedTorsion[2]);
    // **This is the assertion that says the tie is consistent rather than merely
    // close.** A formulation error would sit still under refinement; a
    // discretisation error goes away. It falls monotonically and by twenty-five
    // fold over a threefold refinement, which is faster than first order.
    expectTrue("the tie's EI error falls monotonically under refinement",
               tiedBendingError[1] < tiedBendingError[0] &&
                   tiedBendingError[2] < tiedBendingError[1]);
    expectTrue("and by more than tenfold over a threefold refinement",
               tiedBendingError[0] / tiedBendingError[2] > 10.0);
    // The guard that keeps that from being vacuous: the coarse error has to be
    // something worth converging from, or "it fell" is a statement about rounding.
    expectTrue("from a coarse error worth converging from", tiedBendingError[0] > 1e-3);
    // And the tied torsion converges too -- downwards, like the open one, and to a
    // value that is still an order of magnitude above the open section's.
    expectTrue("the tied torsion also falls monotonically",
               tiedTorsion[1] < tiedTorsion[0] && tiedTorsion[2] < tiedTorsion[1]);
    expectTrue("and stays more than ten times the open section's",
               tiedTorsion[2] / torsion[2] > 10.0);
}

// --- 6. What the mesher refuses ---------------------------------------------------

void testRefusals() {
    const StructuralMesh structure = makeBox(kT, kT);

    // **`empty()` is `elementCount() == 0` and never reads `problems`**, so the
    // three refusals below are three ways of reaching one assertion. Deleting all
    // three `report(...)` calls while leaving the returns -- the header's promise of
    // "a full account in `problems`" removed entirely, the refusals kept -- left the
    // suite at 200 420 checks and zero failures. So each refusal now says which one
    // it is. A section that came back empty for an unrelated reason, or a guard that
    // fired with another guard's text, is a different bug from the one under test
    // and is what these strings separate.
    section::SectionParams backwards = boxParams();
    backwards.xTo = backwards.xFrom;
    const section::Section outOfOrder = section::buildSection(structure, backwards);
    expectTrue("cut planes out of order are refused", outOfOrder.empty());
    expectTrue("and say it is the planes, not the plating",
               !outOfOrder.problems.empty() &&
                   outOfOrder.problems.front().find("not in order") != std::string::npos);

    section::SectionParams noSubdivision = boxParams();
    noSubdivision.subdivision = 0;
    const section::Section unsubdivided = section::buildSection(structure, noSubdivision);
    expectTrue("a subdivision below one is refused", unsubdivided.empty());
    expectTrue("and say it is the subdivision",
               !unsubdivided.problems.empty() &&
                   unsubdivided.problems.front().find("subdivision") != std::string::npos);

    section::SectionParams noRoles = boxParams();
    noRoles.shell = noRoles.deck = noRoles.bulkhead = false;
    const section::Section roleless = section::buildSection(structure, noRoles);
    expectTrue("a section of no roles at all is refused", roleless.empty());
    // This one is the reason the strings matter. "No plating of the requested roles"
    // and "the cut planes are empty of ship" are the same empty section and
    // different bugs, and the box fixture has plating across the whole of its
    // length -- so a reader who saw only `empty()` could not tell which happened.
    expectTrue("and say it is the roles that emptied it, not the station",
               !roleless.problems.empty() &&
                   roleless.problems.front().find("requested roles") != std::string::npos);

    // A cut plane that is not a panel seam. Panels the plane passes through have no
    // corner set that belongs to either side, so they are dropped -- and *counted*,
    // because a section quietly missing a bay is indistinguishable from a shorter
    // ship.
    section::SectionParams offSeam = boxParams();
    offSeam.xFrom = 0.25;  // a quarter of the way into the first bay
    const section::Section ragged = section::buildSection(structure, offSeam);
    expectTrue("a cut that is not on a seam reports the panels it went through",
               ragged.straddlingPanels > 0);
    expectTrue("and a cut that is on one reports none",
               section::buildSection(structure, boxParams()).straddlingPanels == 0);

    // A plate floating inside the section, touching neither cut plane. It is a
    // mechanism in `K_ii` that `reduction::Substructure` accepts -- its precondition
    // check is about the interface, and a free component leaves a tiny *positive*
    // pivot rather than a zero one -- so the mesher has to be the thing that
    // notices.
    StructuralMesh floating = structure;
    {
        PlatePanel p;
        p.corner[0] = {3.0, -0.5, 0.4};
        p.corner[1] = {3.5, -0.5, 0.4};
        p.corner[2] = {3.5, 0.5, 0.4};
        p.corner[3] = {3.0, 0.5, 0.4};
        p.thickness = kT;
        p.role = PanelRole::Deck;
        floating.panels.push_back(p);
    }
    const section::Section adrift = section::buildSection(floating, boxParams());
    expectEqual("a plate touching neither cut plane is reported as floating",
                adrift.floatingComponents, 1);
    // And the tie did not reach it: it is half a metre inboard of the side plates,
    // which is twenty times `junctionTolerance`, so nothing joined it and the
    // mechanism survives to be reported. Without this the test would pass on a tie
    // that joined everything to everything.
    expectEqual("a plate in fresh air is tied to nothing", adrift.junctionTies,
                section::buildSection(structure, boxParams()).junctionTies);
    const section::BeamResponse refused =
        section::applyBeamLoad(adrift, ah36Steel(), section::BeamLoad{1e-6, 0, 0});
    expectTrue("and loading a section with one in it is refused rather than solved", !refused.ok);

    // And the case in between: a plate reaching one cut plane and stopping. It is
    // restrained, so nothing refuses it, and it carries none of the section's axial
    // or bending stiffness -- a cantilever hanging off the interface. `floating`
    // cannot see it and `spanning` has to.
    StructuralMesh halfLength = structure;
    {
        PlatePanel p;
        p.corner[0] = {0.0, -0.5, 0.4};
        p.corner[1] = {kL / 2, -0.5, 0.4};
        p.corner[2] = {kL / 2, 0.5, 0.4};
        p.corner[3] = {0.0, 0.5, 0.4};
        p.thickness = kT;
        p.role = PanelRole::Deck;
        halfLength.panels.push_back(p);
    }
    const section::Section stub = section::buildSection(halfLength, boxParams());
    expectEqual("a plate reaching one plane does not float", stub.floatingComponents, 0);
    // Two components, not five: the box's own four are tied into one, and the stub
    // is the second. Taken against `openBoxParams` below, which is the count this
    // assertion had before the tie existed.
    expectEqual("but it is a component", stub.components, 2);
    expectEqual("and it does not span the section", stub.spanningComponents, 1);
    const section::Section untiedStub = section::buildSection(halfLength, openBoxParams());
    expectEqual("where untied the same input is five pieces", untiedStub.components, 5);
    expectEqual("four of which span", untiedStub.spanningComponents, 4);
}

void testBandwidthReducingOrder() {
    // `reduction::bandwidthReducingOrder` was internal to the reduction until this
    // file needed it, and a function newly on a caller's path with no test of its own
    // is a shape this repo has been bitten by. Two things are its contract: it
    // returns a **permutation** -- every node once, or the caller's renumbering
    // silently drops part of the mesh -- and it narrows the band it is given.
    //
    // A path graph numbered by the bit-reversal of its index: adjacent nodes are far
    // apart in the numbering it arrives in, and one apart in any traversal order.
    constexpr int kBits = 6, kNodes = 1 << kBits;
    const auto reversed = [&](int i) {
        int out = 0;
        for (int b = 0; b < kBits; ++b) out |= ((i >> b) & 1) << (kBits - 1 - b);
        return out;
    };
    std::vector<std::vector<std::uint32_t>> adjacency(kNodes);
    for (int i = 0; i + 1 < kNodes; ++i) {
        adjacency[static_cast<std::size_t>(reversed(i))].push_back(
            static_cast<std::uint32_t>(reversed(i + 1)));
        adjacency[static_cast<std::size_t>(reversed(i + 1))].push_back(
            static_cast<std::uint32_t>(reversed(i)));
    }
    const std::vector<std::uint32_t> order = reduction::bandwidthReducingOrder(adjacency);
    expectEqualCount("the ordering covers every node", order.size(), kNodes);
    std::vector<std::uint8_t> seen(kNodes, 0);
    for (std::uint32_t node : order) seen[node] = 1;
    expectEqual("exactly once each",
                static_cast<long long>(std::count(seen.begin(), seen.end(), 1)), kNodes);

    std::vector<std::size_t> place(kNodes, 0);
    for (std::size_t i = 0; i < order.size(); ++i) place[order[i]] = i;
    std::size_t before = 0, after = 0;
    for (std::size_t i = 0; i < adjacency.size(); ++i)
        for (std::uint32_t j : adjacency[i]) {
            before = std::max(before, i > j ? i - j : j - i);
            after = std::max(after, place[i] > place[j] ? place[i] - place[j]
                                                       : place[j] - place[i]);
        }
    // A path has a bandwidth of one under any traversal of it, and the bit-reversed
    // numbering it arrives in has a bandwidth of half the graph. Asserting the
    // *value* rather than "it got better" is what makes this a test of the
    // algorithm rather than of the direction of an inequality.
    expectEqualCount("a path graph reorders to a bandwidth of one", after, 1);
    expectTrue("from a numbering that was far worse", before > kNodes / 4);
}

void testFreeEdgesInFreshAirAreNotJunctions() {
    // A single flat plate. Its two long edges are free and there is nothing behind
    // them, so `junctionEdges` must be zero -- otherwise the junction census is
    // measuring "this edge is free" rather than "this edge is sitting on plating it
    // is not welded to", and the ferry's 370 m would mean nothing.
    StructuralMesh flat;
    flat.materials = {ah36Steel()};
    for (int i = 0; i < kNx; ++i) {
        const double x0 = kL * i / kNx, x1 = kL * (i + 1) / kNx;
        for (int j = 0; j < kNy; ++j) {
            const double y0 = -1 + 2.0 * j / kNy, y1 = -1 + 2.0 * (j + 1) / kNy;
            PlatePanel p;
            p.corner[0] = {x0, y0, 0};
            p.corner[1] = {x1, y0, 0};
            p.corner[2] = {x1, y1, 0};
            p.corner[3] = {x0, y1, 0};
            p.thickness = kT;
            p.role = PanelRole::Shell;
            flat.panels.push_back(p);
        }
    }
    const section::Section plate = section::buildSection(flat, boxParams());
    expectEqual("a flat plate is one surface", plate.surfaces, 1);
    expectNear("with its two long edges free", plate.freeEdgeLength, 2 * kL, 1e-9);
    expectNear("and not one metre of junction", plate.junctionEdges, 0.0, 1e-12);
    // The box is the control: there, every free edge *is* a junction.
    expectNear("where the box's free edge is all junction",
               section::buildSection(makeBox(kT, kT), boxParams()).junctionEdges, 8 * kL, 1e-9);
    // The node numbering is chosen by comparing candidate orderings, and the wrong
    // choice is a hundredfold in the factorisation rather than a wrong answer. 62 is
    // what the box delivers; the lexicographic ordering that loses here gives more.
    expectTrue("the mesher picks a numbering with a narrow band",
               section::buildSection(makeBox(kT, kT), openBoxParams()).halfBandwidth <= 80);
    // **And the ordering has to be chosen knowing about the ties**, because a tie
    // couples nodes no element edge joins and an ordering scored on the sub-quads
    // alone is scored on the wrong graph. Scoring it on the element graph gave the
    // box 62 and the ferry hold 146, and both were fictions: the assembled band was
    // 1 337 and 10 769. Joining the box genuinely costs a factor of four here; the
    // assertion is that it costs four and not twenty.
    const section::Section tied = section::buildSection(makeBox(kT, kT), boxParams());
    expectTrue("and one that knows the junction ties are part of the graph",
               tied.halfBandwidth <= 250);
    expectTrue("which is a real cost, not a free lunch", tied.halfBandwidth > 150);
    // The band the mesher reports has to be the band the assembly delivers, or it is
    // a decoration. Rebuild it here from the mesh and the constraints, the way
    // `solveStatic` does.
    std::vector<std::vector<std::uint32_t>> stands(tied.mesh.nodeCount() * 3);
    for (std::size_t d = 0; d < stands.size(); ++d)
        stands[d].push_back(static_cast<std::uint32_t>(d));
    for (const solidshell::Mpc& mpc : tied.attachment.constrained) stands[mpc.slave] = mpc.master;
    std::size_t assembled = 0;
    for (std::size_t e = 0; e < tied.mesh.elementCount(); ++e) {
        std::uint32_t lo = 0xffffffffu, hi = 0;
        for (int a = 0; a < 8; ++a) {
            const std::uint32_t n = tied.mesh.index[e * 8 + static_cast<std::size_t>(a)];
            for (int k = 0; k < 3; ++k)
                for (std::uint32_t t : stands[n * 3 + static_cast<std::size_t>(k)]) {
                    lo = std::min(lo, t);
                    hi = std::max(hi, t);
                }
        }
        assembled = std::max(assembled, static_cast<std::size_t>(hi - lo));
    }
    std::printf("     box band: reported %zu, assembled %zu\n", tied.halfBandwidth, assembled);
    expectTrue("the reported band covers what the assembly needs",
               assembled <= tied.halfBandwidth + 3);
}

// --- 7. The reduction consumes it, and Guyan is the same solve -------------------

void testSubstructureConsumesTheSection() {
    const StructuralMesh structure = makeBox(kT, kT);
    const section::Section box = section::buildSection(structure, boxParams());
    const StructuralMaterial material = ah36Steel();

    reduction::Substructure substructure(box.mesh, material, box.interfaceNodes, box.attachment);
    expectTrue("the section is a substructure", substructure.ready());
    expectEqualCount("every interface node became boundary degrees of freedom",
                     substructure.boundaryCount(), box.interfaceNodes.size() * 3);
    // The element mass lumping and `area * thickness * density` are two different
    // routes to the same steel; on a flat box they agree to rounding.
    expectNear("the substructure weighs what the section does", substructure.totalMass(),
               box.mass(), 1e-6 * box.mass());

    const reduction::ReduceParams params{.modes = 0, .verifyModes = false};
    const reduction::Reduction reduced = reduction::craigBampton(substructure, params);
    expectTrue("it reduces", !reduced.empty());
    expectEqual("with no modes, so this is Guyan", reduced.modes, 0);

    // The identity `section.hpp` claims for `applyBeamLoad` and `applyTwist`: they
    // prescribe the whole interface and relax the interior, which *is* static
    // condensation, so the strain energy they report is `0.5 u_b^T K_r u_b` for this
    // reduction. Checked on the twist, because a twist prescribes every interface
    // degree of freedom and so `u_b` can be written down from the geometry alone --
    // there is nothing to take from the solve being checked.
    const double twist = 1e-6;
    const section::TorsionResponse twisted = section::applyTwist(box, material, twist, kH / 2);
    expectTrue("the twist solved: " + twisted.problem, twisted.ok);

    std::vector<double> boundary(static_cast<std::size_t>(reduced.size()), 0.0);
    for (int b = 0; b < reduced.boundary; ++b) {
        const std::uint32_t dof = substructure.boundaryDof()[static_cast<std::size_t>(b)];
        const std::uint32_t node = dof / 3, axis = dof % 3;
        const double x = box.mesh.position[static_cast<std::size_t>(node) * 3];
        if (x < 0.5 * kL) continue;  // the aft plane is held at zero
        const double y = box.mesh.position[static_cast<std::size_t>(node) * 3 + 1];
        const double z = box.mesh.position[static_cast<std::size_t>(node) * 3 + 2] - kH / 2;
        if (axis == 1) boundary[static_cast<std::size_t>(b)] = -twist * z;
        if (axis == 2) boundary[static_cast<std::size_t>(b)] = twist * y;
    }
    double energy = 0;
    for (int i = 0; i < reduced.size(); ++i)
        for (int j = 0; j < reduced.size(); ++j)
            energy += 0.5 * boundary[static_cast<std::size_t>(i)] *
                      reduced.stiffness[static_cast<std::size_t>(i) *
                                            static_cast<std::size_t>(reduced.size()) +
                                        static_cast<std::size_t>(j)] *
                      boundary[static_cast<std::size_t>(j)];
    expectTrue("the twist stored something to compare", twisted.strainEnergy > 1e-6);
    expectNear("the Guyan reduced energy is the section's own", energy / twisted.strainEnergy, 1.0,
               1e-8);

    // Rigid body: a translation of the whole interface carries no strain energy and
    // the reduced mass of it is the substructure's own mass. Both are exact and both
    // go through the same `T`, so they check the mass reduction rather than the
    // stiffness one.
    std::vector<double> rigid(static_cast<std::size_t>(reduced.size()), 0.0);
    for (int b = 0; b < reduced.boundary; ++b)
        if (substructure.boundaryDof()[static_cast<std::size_t>(b)] % 3 == 0)
            rigid[static_cast<std::size_t>(b)] = 1.0;
    double rigidMass = 0;
    for (int i = 0; i < reduced.size(); ++i)
        for (int j = 0; j < reduced.size(); ++j)
            rigidMass += rigid[static_cast<std::size_t>(i)] *
                         reduced.mass[static_cast<std::size_t>(i) *
                                          static_cast<std::size_t>(reduced.size()) +
                                      static_cast<std::size_t>(j)] *
                         rigid[static_cast<std::size_t>(j)];
    expectNear("a rigid translation of the reduced model weighs the substructure", rigidMass,
               substructure.totalMass(), 1e-10 * substructure.totalMass());
}

// --- 8. The ferry, against `hullGirderSection` -----------------------------------

StructuralMesh ferryStructure() {
    Ship ferry = game::buildFerry();
    return makeStructuralMesh(ferry.hull, ferryScantlings());
}

section::SectionParams ferryParams() {
    section::SectionParams p;
    // Two frame bays about midship. Frame stations are multiples of 2.4 m, so this
    // is a cut on a panel seam; the ship's watertight bulkheads are at -44, -38, -8,
    // 20 and 44 and none of them is, which `testFerryBulkheadCutIsNotASeam` is about.
    p.xFrom = -2.4;
    p.xTo = 2.4;
    p.subdivision = 1;
    return p;
}

void testFerrySection() {
    const StructuralMesh structure = ferryStructure();
    const section::Section hold = section::buildSection(structure, ferryParams());
    const HullGirderSection girder = hullGirderSection(structure, 0.0);
    const StructuralMaterial material = ah36Steel();
    const double youngs = material.youngsModulus;

    expectTrue("the ferry section meshed something", !hold.empty());
    expectEqual("a frame-station cut takes no panel with it", hold.straddlingPanels, 0);
    expectEqual("nothing floats free of the interface", hold.floatingComponents, 0);
    expectEqual("and every piece reaches both cut planes", hold.spanningComponents, hold.components);

    // The finding this whole file was built around: `makeStructuralMesh` shares no
    // corner between two panel roles, so a section of a real ship arrives in pieces
    // however it is meshed, and hundreds of metres of its free edge are lying on
    // plating they are not joined to. That is still true of the *mesh* -- and the
    // junction tie is what stops it being true of the model.
    const section::Section untied =
        section::buildSection(structure, [] {
            section::SectionParams p = ferryParams();
            p.junctions = false;
            return p;
        }());
    expectTrue("the ferry's mesh is in several disconnected pieces", untied.components > 1);
    expectTrue("with tens of metres of unwelded junction", untied.junctionEdges > 10.0);
    expectTrue("and the junction really is a junction rather than a gap",
               untied.worstJunctionGap < 0.02);
    expectEqual("the tie makes it one", hold.components, 1);
    expectEqual("which reaches both cut planes", hold.spanningComponents, 1);
    expectNear("over the same plating", hold.junctionEdges, untied.junctionEdges, 1e-9);
    expectEqualCount("with the same elements", hold.elementCount(), untied.elementCount());
    expectEqualCount("and the same nodes", hold.nodeCount(), untied.nodeCount());
    expectTrue("some of that junction is now tied", hold.tiedEdges > 0.0);
    expectNear("and none of it was before", untied.tiedEdges, 0.0, 0.0);
    // Two bays is three stations, of which two are cut planes, so most of this
    // section's junction nodes are on the interface and cannot be tied -- which is
    // why the frequency measurement below runs on eight bays instead. It is asserted
    // here so that the small `tiedEdges` is a stated property rather than a
    // surprise.
    // Two bays is three stations, of which two are cut planes, so most of this
    // section's junction nodes are on the interface -- twice as many as are not. They
    // are tied to the line in their own plane rather than left open (§9), which is
    // what the chain applies and a lone section does not; `tiedEdges` is small here
    // for that reason and it is asserted so that it is a stated property.
    expectTrue("most of a two-bay section's junction nodes are on a cut plane",
               hold.planeTieNodes > hold.junctionTies);
    expectEqual("and none of them is left open", hold.junctionsOnInterface, 0);
    expectNear("what the section itself joins plus what its two planes would is all of it",
               hold.tiedEdges + hold.planeTiedEdgesAft + hold.planeTiedEdgesForward +
                   hold.planeTiedEdgesBoth,
               hold.junctionEdges, 1e-9);
    // A line drops whatever the slave's own normal leans along the ship -- the one
    // thing a face tie has no analogue of -- bounded above by half a plate thickness,
    // 7.8 mm on this hull's bilge strake. Measured at **1.9e-5 m**, four hundred times
    // below that bound, because a transverse cut of a ship is very nearly square to
    // the plating it passes through. Asserted at 1e-4, which is what was measured with
    // a factor of five; the half-thickness bound would be nearly a vacuous assertion.
    expectTrue("and the line drops a fiftieth of a millimetre, not half a thickness",
               hold.worstPlaneTieSlip < 1e-4);
    expectTrue("which is not zero, so there was something to drop", hold.worstPlaneTieSlip > 0);
    // A real ship's junction has a real gap -- 9 mm at worst on this hull -- so the
    // through-thickness weight extrapolates rather than interpolating, and the tie
    // has to survive a weight outside [0, 1] without producing a negative nodal
    // mass. `reduction::Substructure` is the thing that would refuse it.
    expectTrue("a real junction extrapolates through the master's thickness",
               hold.worstJunctionWeight > 1.0);
    const reduction::Substructure substructure(hold.mesh, hold.material, hold.interfaceNodes,
                                               hold.attachment);
    expectTrue("and the substructure still finds every nodal mass positive",
               substructure.ready());
    // Not against `hold.mass()`, which is `area * thickness * density` and differs
    // from the element lumping by 0.07% on a curved hull -- against the *same*
    // substructure without the ties, which is the property being claimed: an
    // eliminated degree of freedom's steel goes to its masters and none of it is
    // lost. Exact, because the weights are a partition of unity.
    const reduction::Substructure loose(untied.mesh, untied.material, untied.interfaceNodes,
                                        untied.attachment);
    expectTrue("the untied section is a substructure too", loose.ready());
    expectNear("with the section's steel, none of it lost to the elimination",
               substructure.totalMass(), loose.totalMass(), 1e-9 * loose.totalMass());

    // The member accounting. Every longitudinally effective member is either
    // attached or reported missing, and the two together are exactly what
    // `sectionElements` says the stiffeners are worth -- so the section coming out
    // short against the girder below is an accounting rather than a discrepancy.
    double stiffenerArea = 0, girderArea = 0;
    for (const SectionElement& element : sectionElements(structure, 0.0))
        if (element.stiffener) stiffenerArea += element.area;
    for (const StructuralMember& member : structure.members) {
        if (member.role != MemberRole::Girder) continue;
        const double lo = std::min(member.a.x, member.b.x), hi = std::max(member.a.x, member.b.x);
        if (!(lo <= 0.0 && 0.0 < hi)) continue;
        girderArea += profileSection(member.profile).area;
    }
    expectNear("attached plus missed is the stiffener area of the cut",
               hold.attachedMemberArea + hold.missedMemberArea, stiffenerArea, 1e-9);
    expectTrue("and the stiffeners are a quarter of the section, so that is not a small claim",
               stiffenerArea / girder.area > 0.2);
    // What is missing is the girders, and nothing else: they sit off the
    // longitudinal spacing, so a mesh whose nodes are panel seams has no node on
    // them.
    expectTrue("the girders exist to be missed", girderArea > 0.05);
    expectNear("and they are exactly what was missed", hold.missedMemberArea, girderArea, 1e-9);

    section::BeamLoad axial;
    axial.strain = 1e-6;
    axial.reference = girder.neutralAxis;
    const section::BeamResponse stretched = section::applyBeamLoad(hold, material, axial);
    expectTrue("the ferry section takes an axial load: " + stretched.problem, stretched.ok);
    expectTrue("its rigid-body restraints carry nothing", stretched.restraintReaction < 1e-3);

    // The comparison the whole exercise is for: `hullGirderSection` reaches `A`,
    // `z_na` and `I` by summing over a transverse cut and shares no line of code
    // with the mesher, the element or the solver.
    const double effectiveArea = stretched.axialStiffness / youngs;
    const double predictedArea = girder.area - hold.missedMemberArea;
    expectNear("EA against the girder's own area, less the members it could not attach",
               effectiveArea / predictedArea, 1.0, 0.01);
    // **0.356% measured on this window, not the 0.44% this comment carried**, and it is
    // not noise: the frames and deck beams restrain the section's Poisson contraction,
    // which a beam model has no way to represent. They are worth 0.439% of the area
    // here, measured by omitting them in `testFramesRestrainThePoissonContraction`
    // below -- where both figures come out of a run rather than out of a comment.
    //
    // **The old pair, 0.44% and 0.52%, were the *hold*'s figures written into a
    // two-bay test, and they were the hold's figures before §8's halo.** The hold now
    // gives +0.411% and +0.508% (`section_probe --from=-7.2 --to=19.2`, with and
    // without `--no-frames`), against +0.420% before the halo moved
    // `missedMemberArea` -- so one of the two numbers was near enough to be believed
    // and neither belonged to the section this test builds. The effect is a function
    // of section length: 0.31% at one bay, 0.44% at two, 0.48% at four, tending to the
    // hold's 0.51%, because a cut plane is free in y and z and a longer section has
    // proportionally less of itself next to one. Nothing tests a comment.
    expectTrue("the agreement is not vacuous: the girders alone are four per cent",
               hold.missedMemberArea / girder.area > 0.03);

    const double neutralAxis = stretched.bendingMoment / stretched.axialForce + girder.neutralAxis;
    expectNear("the neutral axis the section finds for itself", neutralAxis, 6.8534, 0.02);
    expectTrue("which is not simply the reference it was given",
               std::abs(neutralAxis - girder.neutralAxis) > 0.1);

    section::BeamLoad bending;
    bending.curvature = 1e-6;
    bending.reference = girder.neutralAxis;
    const section::BeamResponse bent = section::applyBeamLoad(hold, material, bending);
    expectTrue("the ferry section takes a curvature: " + bent.problem, bent.ok);
    // The girders carry 2.459 m^4 of the girder's 46.205, by the same decomposition
    // that produced `girderArea` above -- taken about the section's own neutral axis,
    // so it neglects the shift of that axis when they are removed. Carrying that
    // through too gives 43.70 rather than 43.75, and the measured value sits inside
    // one per cent of either, which is why the tolerance is not tighter.
    expectNear("EI against the girder's own second moment, less the girders",
               bent.bendingStiffness / youngs / (girder.secondMoment - 2.459), 1.0, 0.01);
}

// --- 8b. The measurement the junctions were costing --------------------------------
//
// **`EA` and `EI` cannot see a junction and this can.** Prescribing plane sections
// at two cuts makes every longitudinally continuous strip carry `sigma = E eps`
// whatever it is attached to, so an unjoined section reports the right hull girder.
// A fixed-interface frequency has no such crutch: it is the softest thing the
// section can do with its interface held, and on this ship that is a deck spanning
// twenty-six metres on two edges instead of four.
//
// Eight bays rather than the eleven `tools/section_probe` uses, because the effect
// grows with length and eight is where it is still unmistakable inside a unit gate:
// 1.84x here against 2.96x there.
void testFerryJunctionMovesTheLowestFrequency() {
    const StructuralMesh structure = ferryStructure();
    const auto lowest = [&](bool shell, bool deck, bool bulkhead, bool junctions,
                            const char* label) {
        section::SectionParams p;
        p.xFrom = -7.2;
        p.xTo = -7.2 + 2.4 * 8;
        p.subdivision = 1;
        p.shell = shell;
        p.deck = deck;
        p.bulkhead = bulkhead;
        p.junctions = junctions;
        const section::Section piece = section::buildSection(structure, p);
        const reduction::Substructure substructure(piece.mesh, piece.material,
                                                   piece.interfaceNodes, piece.attachment);
        expectTrue(std::string(label) + " is a substructure", substructure.ready());
        const reduction::Eigenpairs modes = substructure.fixedInterfaceModes(1);
        expectTrue(std::string(label) + " has a lowest mode", !modes.value.empty());
        const double omega = modes.value.empty() ? 0.0 : std::sqrt(std::max(0.0, modes.value[0]));
        // Subspace iteration converges to *a* mode; `eigenvaluesBelow` counts by the
        // inertia of an LDL^T factorisation, which is a different instrument
        // answering the same question, and it is what says nothing was skipped.
        expectEqual(std::string(label) + " has nothing below what was found",
                    substructure.eigenvaluesBelow(0.99 * 0.99 * omega * omega), 0);
        expectTrue(std::string(label) + " has something at it",
                   substructure.eigenvaluesBelow(1.01 * 1.01 * omega * omega) >= 1);
        std::printf("     %-24s %d component(s), band %3zu (mesh %4zu): %.4f Hz\n", label,
                    piece.components, substructure.halfBandwidth(), piece.halfBandwidth,
                    omega / (2 * std::numbers::pi));
        if (junctions) {
            // **The ordering has to be chosen on the tied graph, and on real geometry
            // that is worth more than it is on the box.** A tie couples nodes no
            // element edge joins; leaving those edges out of the Cuthill-McKee
            // adjacency gives 1 910 here against 740, and the box cannot see the
            // difference at all because its four corner lines are already adjacent in
            // every candidate ordering. Measured 740; asserted at 1 200 so a change of
            // ordering heuristic is allowed and a loss of the tie edges is not.
            expectTrue("the ordering knows about the ties on real geometry",
                       piece.halfBandwidth < 1200);
            expectTrue("which is a real cost and not a free lunch", piece.halfBandwidth > 400);
            // What the tie's geometry costs, on plating that is genuinely warped: the
            // bilinear surface through a master face's four corners does not pass
            // exactly through the node being tied, and the gap is the modelling error
            // the tie carries. 1.94e-05 m on this ship, against a 2.4 m element.
            double worstOffset = 0, worstSum = 0;
            for (const solidshell::Mpc& mpc : piece.attachment.constrained) {
                double sum = 0;
                Vec3 interpolated{0, 0, 0};
                for (std::size_t a = 0; a < mpc.master.size(); ++a) {
                    sum += mpc.weight[a];
                    const std::size_t node = mpc.master[a] / 3;
                    for (int k = 0; k < 3; ++k)
                        interpolated[k] += mpc.weight[a] *
                                           piece.mesh.position[node * 3 + static_cast<std::size_t>(k)];
                }
                worstSum = std::max(worstSum, std::abs(sum - 1.0));
                const std::size_t slave = mpc.slave / 3;
                worstOffset = std::max(
                    worstOffset,
                    length(interpolated - Vec3{piece.mesh.position[slave * 3],
                                               piece.mesh.position[slave * 3 + 1],
                                               piece.mesh.position[slave * 3 + 2]}));
            }
            std::printf("     %-24s tie geometry: |sum w - 1| %.2e, |interp - slave| %.2e m\n", "",
                        worstSum, worstOffset);
            expectTrue("the weights are a partition of unity on real geometry too",
                       worstSum < 4e-16);
            expectTrue("and the tied point is where the node is, to a fiftieth of a millimetre",
                       worstOffset < 5e-5);
        }
        return omega / (2 * std::numbers::pi);
    };

    const double decks = lowest(false, true, false, false, "decks alone");
    const double shell = lowest(true, false, false, false, "shell alone");
    const double open = lowest(true, true, true, false, "whole section, untied");
    const double tied = lowest(true, true, true, true, "whole section, tied");

    // **The defect, as a number.** The untied section's softest mode is the decks'
    // own, to four figures: adding the shell they are welded to in the ship changes
    // it by nothing, because in the model they are not welded to it.
    expectNear("untied, the whole section is exactly as soft as its decks alone", open / decks,
               1.0, 1e-4);
    // The vacuity guard: the shell is a genuinely different structure, so the
    // agreement above is a statement about the junction and not about the two pieces
    // happening to be alike.
    expectTrue("where the shell alone is a different answer entirely", shell / decks > 1.5);

    // And what tying it buys. 1.8415 measured; asserted as a factor rather than a
    // frequency so it does not have to be re-measured when the ferry's scantlings
    // move.
    expectTrue("tying the junctions stiffens the softest mode by more than half again",
               tied / open > 1.5);
    expectTrue("and it is no longer the decks' own mode", std::abs(tied / decks - 1.0) > 0.5);
    // A tie can only add constraint, so it can only raise a fixed-interface
    // frequency -- there is no configuration in which joining the section softens
    // it, and a tie that had gone in with the wrong sign somewhere could.
    expectTrue("a constraint cannot soften the section", tied > open);
}

void testFerryMembersAreWorthWhatTheSectionSays() {
    const StructuralMesh structure = ferryStructure();
    const HullGirderSection girder = hullGirderSection(structure, 0.0);
    const StructuralMaterial material = ah36Steel();
    const double youngs = material.youngsModulus;
    section::BeamLoad axial;
    axial.strain = 1e-6;

    // Bare plating. `reduction.hpp` §8 measured a single flat bar on one patch at
    // 7.9% of its displacement field and warned that longitudinals carry a large
    // share of a hull girder; here is the share, on a real section.
    section::SectionParams bare = ferryParams();
    bare.members = false;
    const section::Section plating = section::buildSection(structure, bare);
    const double platingArea =
        section::applyBeamLoad(plating, material, axial).axialStiffness / youngs;

    double stiffenerArea = 0, plateArea = 0;
    for (const SectionElement& element : sectionElements(structure, 0.0))
        (element.stiffener ? stiffenerArea : plateArea) += element.area;
    expectNear("plating alone is the girder's plating", platingArea / plateArea, 1.0, 0.005);
    expectNear("which is short of the whole section by exactly the stiffeners",
               platingArea / girder.area, 1.0 - stiffenerArea / girder.area, 0.005);
    expectTrue("and that shortfall is 20% or more, so it is not a rounding claim",
               stiffenerArea / girder.area > 0.2);

    // Now put back only the members that run along the ship. The difference has to
    // be exactly the area the `Attachment` says it added -- an athwartships member
    // contributes to neither, so leaving frames and deck beams out isolates the
    // identity from the Poisson restraint they add.
    StructuralMesh longitudinalsOnly = structure;
    {
        std::vector<StructuralMember> kept;
        for (const StructuralMember& member : longitudinalsOnly.members)
            if (member.role == MemberRole::Longitudinal ||
                member.role == MemberRole::DeckLongitudinal ||
                member.role == MemberRole::BulkheadStiffener)
                kept.push_back(member);
        longitudinalsOnly.members = std::move(kept);
    }
    const section::Section stiffened =
        section::buildSection(longitudinalsOnly, ferryParams());
    const double stiffenedArea =
        section::applyBeamLoad(stiffened, material, axial).axialStiffness / youngs;
    expectTrue("the fibres added something", stiffened.attachedMemberArea > 0.3);
    expectNear("and it is exactly the cross-section they stand for",
               (stiffenedArea - platingArea) / stiffened.attachedMemberArea, 1.0, 1e-3);
}

void testFerryFibresAndTheirMass() {
    const StructuralMesh structure = ferryStructure();
    const section::Section hold = section::buildSection(structure, ferryParams());

    expectTrue("the ferry section carries fibres", hold.stiffening.fiberCount() > 100);
    // Some of this ship's members have a web that is not along the plating's
    // thickness direction -- the tie has no single eccentricity for those and the
    // mesher declines them rather than projecting one, which would be a plausible
    // wrong second moment.
    expectTrue("and declines the ones it cannot tie", hold.membersRefused > 0);
    expectTrue("but only a small minority of them",
               hold.membersRefused < hold.membersAttached / 10);
    // A run stops where the plating under it changes thickness, because one tie
    // weight is only right for the pair separation it was computed against.
    expectTrue("and breaks a run where the plate steps", hold.memberRunsSplitByThickness > 0);

    // **The invariant that says the tie weights are right**: every fibre end sits at
    // the offset it records, measured from its own node pair's mid-surface. Feeding
    // `addStiffener` any thickness other than the pair's own separation puts the
    // fibre at `e * t_pair / t_given` -- 47 mm out on this ship, a quarter of a
    // 700 mm frame's Steiner term, and invisible in `EA`, in `EI` and in the mass.
    double worst = 0;
    for (const constraint::Fiber& fiber : hold.stiffening.fiber)
        for (int end = 0; end < 2; ++end) {
            const constraint::Tie& tie = fiber.end[end];
            const Vec3 tied = constraint::tiedPoint(tie, hold.mesh.position);
            const auto at = [&](std::uint32_t node) {
                return Vec3{hold.mesh.position[static_cast<std::size_t>(node) * 3],
                            hold.mesh.position[static_cast<std::size_t>(node) * 3 + 1],
                            hold.mesh.position[static_cast<std::size_t>(node) * 3 + 2]};
            };
            const Vec3 middle = (at(tie.bottom) + at(tie.top)) * 0.5;
            worst = std::max(worst, std::abs(length(tied - middle) - std::abs(fiber.offset)));
        }
    expectTrue("every fibre sits at the offset it records", worst < 1e-11);
    // Not vacuous: the offsets are of order a hundred millimetres, so 1e-9 is nine
    // orders below the thing being measured rather than a loose band round zero.
    double largest = 0;
    for (const constraint::Fiber& fiber : hold.stiffening.fiber)
        largest = std::max(largest, std::abs(fiber.offset));
    expectTrue("and the offsets are large enough for that to be a statement", largest > 0.05);

    // The `Attachment` the reduction consumes has to carry the steel as well as the
    // stiffness: a stiffened model that is stiffer and no heavier puts every
    // frequency high, and frequencies are most of what Tier 1 is for
    // (`reduction.hpp` §8).
    reduction::Substructure substructure(hold.mesh, hold.material, hold.interfaceNodes,
                                         hold.attachment);
    expectTrue("the ferry section is a substructure", substructure.ready());
    expectNear("and the fibres' steel arrived with them", substructure.attachedMass(),
               hold.stiffening.mass, 1e-9);
    expectTrue("which is a fifth of the section, so that is not a small claim",
               hold.stiffening.mass > 0.2 * hold.mass());
    expectNear("the substructure weighs what the section does", substructure.totalMass(),
               hold.mass(), 0.002 * hold.mass());
}

void testFerryBulkheadCutIsNotASeam() {
    // The obvious place to cut a ship is at a bulkhead, and on this ship it is the
    // wrong place: the bulkheads are at -44, -38, -8, 20 and 44 and the frames are
    // multiples of 2.4, so a bulkhead plane passes through the middle of the shell
    // and deck panels either side of it.
    const StructuralMesh structure = ferryStructure();
    section::SectionParams atBulkheads;
    atBulkheads.xFrom = -8.0;
    atBulkheads.xTo = 20.0;
    atBulkheads.members = false;
    const section::Section hold = section::buildSection(structure, atBulkheads);
    expectTrue("a bulkhead cut goes through panels", hold.straddlingPanels > 100);

    section::SectionParams atFrames = atBulkheads;
    atFrames.xFrom = -7.2;
    atFrames.xTo = 19.2;
    const section::Section clean = section::buildSection(structure, atFrames);
    expectEqual("moving it to the nearest frames takes none", clean.straddlingPanels, 0);
    // The bulkhead cut is the *longer* piece of ship and still ends up with less
    // plating per metre of it, because the panels the planes went through are simply
    // gone. That is the guard: a section short of 2.5% of its own plating would
    // otherwise look like a slightly smaller ship.
    expectTrue("the bulkhead cut is the longer piece", hold.length() > clean.length());
    expectTrue("and yet it carries less plating per metre",
               hold.area / hold.length() < 0.99 * clean.area / clean.length());
}

void testInterfaceIsChosenOnTheMidSurface() {
    // `reduction::nodesNearPlanes` at anything like its default tolerance keeps only
    // the node of each through-thickness pair that happened to land on the plane, and
    // the plating's normal leans out of the transverse direction wherever the hull is
    // not parallel-sided. Half an interface is not a cut, it is a hinge, so the
    // mesher chooses the interface on the mid-surface and takes both nodes -- and
    // this is the measurement that says the two are different.
    const StructuralMesh structure = ferryStructure();
    section::SectionParams params = ferryParams();
    params.xFrom = 36.0;   // forward, where the hull is narrowing fast
    params.xTo = 40.8;
    params.members = false;
    const section::Section bow = section::buildSection(structure, params);
    expectTrue("the bow section meshed something", !bow.empty());

    const std::vector<reduction::Plane> planes{{{params.xFrom, 0, 0}, {1, 0, 0}},
                                               {{params.xTo, 0, 0}, {1, 0, 0}}};
    const std::vector<std::uint32_t> byNode = reduction::nodesNearPlanes(bow.mesh, planes, 1e-9);
    expectTrue("choosing the interface by node position loses some of it",
               byNode.size() < bow.interfaceNodes.size());
    // Every node the node-based rule found is one the mid-surface rule found too, so
    // the difference is nodes it *missed* rather than a different interface.
    bool subset = true;
    for (std::uint32_t node : byNode)
        subset = subset && std::binary_search(bow.interfaceNodes.begin(),
                                              bow.interfaceNodes.end(), node);
    expectTrue("and what it does find is a subset of what the mesher found", subset);
}

// --- 6. A ship: a chain of sections ----------------------------------------------
//
// The end-to-end claim of the tier. **The reference is the same length meshed as one
// section**, which owes nothing to the assembly -- it goes through
// `solidshell::solveStatic` on the monolithic mesh -- and the box girder is used
// because every quantity in it also has a closed form.
//
// The trap this file already knows about is the one that decides how these tests are
// written. `EA` comes out exact on a chain of eight sections that tie nothing to
// anything: prescribing plane sections at the ends makes every longitudinal strip
// carry `sigma = E eps` whatever it is joined to (§2), and a chain is no different.
// So every claim below that is about **joining** is made on torsion or on a frequency,
// and `EA` appears only as the thing that is exact while `GJ` is 19% out.

section::ChainParams boxChain(int sections, bool junctions) {
    section::ChainParams p;
    p.section.subdivision = 1;
    p.section.members = false;
    p.section.junctions = junctions;
    p.reduce.modes = 0;
    p.reduce.cutoffFrequency = 0;
    for (int i = 0; i <= sections; ++i) p.station.push_back(kL * i / sections);
    return p;
}

constexpr double kBoxTorsion = kBoxTorsionConstant;  // GJ / G, m^4

double shearModulus() {
    const StructuralMaterial steel = ah36Steel();
    return steel.youngsModulus / (2.0 * (1.0 + steel.poissonRatio));
}

// **Like for like: N sections and one section, with the junctions off on both.**
//
// Off on both, because the chain and the monolith are then the *same structure* and
// any difference between them is the assembly. Tied they are not the same structure
// -- an interior cut plane is an interface, an interface DOF is prescribed rather
// than derived, and the chain therefore ties one station fewer per cut -- and mixing
// the two questions is how a 19% error would come to look like round-off.
void testChainOfSectionsIsTheWholeSection() {
    std::printf("\n--- section: a chain of sections against the same length in one piece ---\n");
    const StructuralMesh structure = makeBox(kT, kT);
    const StructuralMaterial steel = ah36Steel();

    const section::Chain chain = section::buildChain(structure, boxChain(4, false));
    const section::Section whole = section::buildSection(structure, openBoxParams());
    expectTrue("the chain built", chain.ready());
    expectEqual("four sections went in", chain.assembly.parts, 4);
    expectEqualCount("the sections together have the elements the whole box has",
                     chain.section[0].elementCount() + chain.section[1].elementCount() +
                         chain.section[2].elementCount() + chain.section[3].elementCount(),
                     whole.elementCount());
    double plate = 0;
    for (const section::Section& s : chain.section) plate += s.plateMass;
    expectNear("and the steel the whole box has", plate, whole.plateMass, 1e-9);

    // Five cut planes, three of them interior and shared. A shared plane must be
    // shared *entirely*: a chain assembles and solves with part of a cut unmatched.
    expectEqualCount("three interior cut planes", chain.shared.size(), 3u);
    for (std::size_t i = 0; i < chain.shared.size(); ++i) {
        expectEqualCount("each shares a whole plane", chain.shared[i],
                         3 * chain.section[i].forwardNodes.size());
        expectEqualCount("with nothing left over", chain.unmatched[i], 0u);
    }
    expectNear("the box's plates are prismatic, so the planes coincide exactly",
               chain.worstGap, 0.0, 0.0);
    // The union of five planes counted once, not the sum of eight section planes.
    expectEqualCount("the assembled boundary is the five planes", chain.boundaryDof(),
                     5 * 3 * chain.section[0].forwardNodes.size());

    // Untied, the box is four plates that touch nothing -- in the chain and in the
    // monolith alike. `components` has to see both: that each section is four pieces,
    // *and* that a piece of one section is the same piece of ship as the piece of the
    // next it shares a cut plane with.
    expectEqual("the untied chain is the four plates the untied box is", chain.components,
                whole.components);
    expectEqual("which is four and not sixteen", chain.components, 4);
    // And the reduced model alone cannot see that, which is why there are two counts.
    expectEqual("the reduced model is one piece even so", chain.reducedComponents, 1);

    section::BeamLoad axial;
    axial.strain = 1e-4;
    section::BeamLoad bending;
    bending.curvature = 1e-5;
    bending.reference = kH / 2;
    const section::BeamResponse chainAxial = section::applyBeamLoad(chain, axial);
    const section::BeamResponse wholeAxial = section::applyBeamLoad(whole, steel, axial);
    const section::BeamResponse chainBend = section::applyBeamLoad(chain, bending);
    const section::BeamResponse wholeBend = section::applyBeamLoad(whole, steel, bending);
    const section::TorsionResponse chainTwist = section::applyTwist(chain, 1e-4, kH / 2);
    const section::TorsionResponse wholeTwist = section::applyTwist(whole, steel, 1e-4, kH / 2);
    expectTrue("every solve ran", chainAxial.ok && wholeAxial.ok && chainBend.ok &&
                                      wholeBend.ok && chainTwist.ok && wholeTwist.ok);

    // Exact, not converged. The assembled boundary is the union of the cut planes,
    // Guyan condensation is exact there for a load-free interior, and the scatter-add
    // approximates nothing -- so the chain is the monolith to the conditioning of two
    // independent solves. Measured at 9.9e-13, 1.4e-10 and 1.3e-10; asserted a decade
    // or two above, because a tolerance far looser than the measurement would pass on
    // a model that had lost the property and was merely well converged.
    std::printf("     untied, 4 sections against 1: EA %+.2e, EI %+.2e, GJ %+.2e\n",
                chainAxial.axialStiffness / wholeAxial.axialStiffness - 1,
                chainBend.bendingStiffness / wholeBend.bendingStiffness - 1,
                chainTwist.torsionalStiffness / wholeTwist.torsionalStiffness - 1);
    expectNear("the chain's EA is the whole section's", chainAxial.axialStiffness,
               wholeAxial.axialStiffness, 1e-11 * wholeAxial.axialStiffness);
    expectNear("and its EI", chainBend.bendingStiffness, wholeBend.bendingStiffness,
               1e-8 * wholeBend.bendingStiffness);
    expectNear("and its GJ", chainTwist.torsionalStiffness, wholeTwist.torsionalStiffness,
               1e-8 * wholeTwist.torsionalStiffness);
    // Against the closed forms, so the pair are not merely agreeing with each other.
    expectNear("and EA is the closed form", chainAxial.axialStiffness,
               steel.youngsModulus * kBoxArea, kAxialTolerance * steel.youngsModulus * kBoxArea);
    // The rigid body restraints pick one of a family of zero-energy motions rather
    // than resisting anything, so their reaction is a defect and not a tolerance --
    // and there are three of them per piece here, twelve in all.
    expectTrue("the restraints carry nothing",
               chainAxial.restraintReaction < 1e-9 * std::fabs(chainAxial.axialForce));

    // **The vacuity guard.** All three agreeing to 1e-10 would also be what a chain
    // that had silently collapsed onto one section produced, so the comparison has to
    // be shown to have something in it: the same box tied is a *different* answer in
    // torsion by an order of magnitude, and the chain is not the monolith's own solve.
    const section::Section tiedWhole = section::buildSection(structure, boxParams());
    const section::TorsionResponse tiedTwist = section::applyTwist(tiedWhole, steel, 1e-4, kH / 2);
    expectTrue("a closed cell carries an order of magnitude more torque than an open one",
               tiedTwist.torsionalStiffness > 10.0 * wholeTwist.torsionalStiffness);
    expectTrue("and the chain solved a model of its own size, not the section's",
               chain.assembly.size() > static_cast<int>(3 * chain.section[0].forwardNodes.size()));
}

// **What a cut plane costs, which `EA` cannot see and torsion can -- and what the
// in-plane line tie of §9 gives back.**
//
// A junction node on a cut plane cannot be tied to a master *face*: half a face's
// nodes are interior to the section, and a boundary degree of freedom written as a
// function of an interior one is not one a reduction keeps exactly. So cutting a
// length into N pieces used to turn N-1 interior stations into open rings of
// junctions, and `GJ` fell by 19% at N = 8 while `EA` stayed exact at every N.
//
// A tie to the *line* the other surface draws on that same plane has every master on
// the plane, so it is a relation between assembled unknowns that both sections either
// side of the cut have. It is applied to the assembled model rather than inside a
// section, and it puts every figure below back on the monolith's.
//
// **The negative control is the whole of the evidence.** `interfaceTies = false` is
// what this file did before and it still reproduces the table it published, so the
// comparison is between two measurements and not between a measurement and a memory.
void testChainCutPlanesUntieTheJunctions() {
    std::printf("\n--- section: what an interior cut plane costs the junctions ---\n");
    const StructuralMesh structure = makeBox(kT, kT);
    const StructuralMaterial steel = ah36Steel();
    const section::Section whole = section::buildSection(structure, boxParams());
    const section::TorsionResponse wholeTwist = section::applyTwist(whole, steel, 1e-4, kH / 2);
    const double bredt = shearModulus() * kBoxTorsion;
    expectTrue("the monolithic box is closed in torsion",
               wholeTwist.torsionalStiffness > 0.9 * bredt);

    double previousOpenTorsion = 2.0 * wholeTwist.torsionalStiffness;
    double previousOpenTied = 2.0 * whole.tiedEdges;
    for (int sections : {1, 2, 4, 8}) {
        const section::Chain chain = section::buildChain(structure, boxChain(sections, true));
        section::ChainParams openParams = boxChain(sections, true);
        openParams.section.interfaceTies = false;
        const section::Chain open = section::buildChain(structure, openParams);
        expectTrue("both chains built", chain.ready() && open.ready());
        expectEqual("and each is one piece of ship", chain.components, 1);
        expectEqual("with or without the line tie", open.components, 1);
        section::BeamLoad axial;
        axial.strain = 1e-4;
        const section::BeamResponse stretched = section::applyBeamLoad(chain, axial);
        const section::TorsionResponse twisted = section::applyTwist(chain, 1e-4, kH / 2);
        const section::BeamResponse openStretched = section::applyBeamLoad(open, axial);
        const section::TorsionResponse openTwisted = section::applyTwist(open, 1e-4, kH / 2);
        expectTrue("every solve ran",
                   stretched.ok && twisted.ok && openStretched.ok && openTwisted.ok);
        std::printf("     %d sections: tied %5.1f m of %.1f, EA %+.2e, GJ %+.2e"
                    "   |  open: %5.1f m, EA %+.2e, GJ %+.3f%%\n",
                    sections, chain.tiedEdges, chain.junctionEdges,
                    stretched.axialStiffness / (steel.youngsModulus * kBoxArea) - 1,
                    twisted.torsionalStiffness / wholeTwist.torsionalStiffness - 1,
                    open.tiedEdges,
                    openStretched.axialStiffness / (steel.youngsModulus * kBoxArea) - 1,
                    100.0 * (openTwisted.torsionalStiffness / wholeTwist.torsionalStiffness - 1));

        // **`EA` says nothing, in either column, at any N.** It is here only as the
        // thing that is exact while `GJ` is 19% out, which is the whole reason a
        // validation of this that stopped at `EA` would have proved nothing.
        expectNear("EA is the closed form at every N", stretched.axialStiffness,
                   steel.youngsModulus * kBoxArea, kAxialTolerance * steel.youngsModulus * kBoxArea);
        expectNear("and it is the closed form with the ties off too", openStretched.axialStiffness,
                   steel.youngsModulus * kBoxArea, kAxialTolerance * steel.youngsModulus * kBoxArea);

        // The line ties: four corners on each of the N-1 interior planes, and the two
        // sections either side of every one of them derived the same constraint.
        expectEqual("four corners tied on each interior plane", chain.planeTieNodes,
                    4 * (sections - 1));
        expectEqual("and the two sides of every cut agreed about all of it",
                    chain.planeTiesDisagreeing, 0);
        expectNear("to the last bit", chain.worstPlaneTieDisagreement, 0.0, 0.0);
        expectEqual("with the line tie off there are none at all", open.planeTieNodes, 0);

        // **The claim.** A chain ties what the same length in one piece ties -- 58 m
        // of the box's 64, the missing 6 being the two outermost planes, which a
        // monolith leaves open as well because they are what a load is prescribed on.
        expectNear("a chain joins what the monolith joins, at every N", chain.tiedEdges,
                   whole.tiedEdges, 1e-9);
        // And torsion, which is the only quantity here that can tell. Measured at
        // 4e-13, 2.2e-13 and 1.7e-13 relative for N = 2, 4, 8 -- the conditioning of
        // two independent solves and not a truncation. Asserted at 1e-9, which is four
        // orders above the measurement and ten below the 3.3% the first cut costs.
        expectNear("and it carries the torsion the monolith carries", twisted.torsionalStiffness,
                   wholeTwist.torsionalStiffness, 1e-9 * wholeTwist.torsionalStiffness);

        // --- The negative control, which is what makes the line above mean anything --
        expectTrue("with the line tie off, every extra cut plane unties another station",
                   open.tiedEdges < previousOpenTied);
        expectTrue("and costs torsional stiffness",
                   openTwisted.torsionalStiffness < previousOpenTorsion);
        previousOpenTied = open.tiedEdges;
        previousOpenTorsion = openTwisted.torsionalStiffness;
        if (sections == 1) {
            // One section is the monolith reduced and reassembled: there is no interior
            // plane, so the two columns are the same model and both must reproduce it.
            expectNear("a chain of one is the section it was cut from",
                       twisted.torsionalStiffness, wholeTwist.torsionalStiffness,
                       1e-9 * wholeTwist.torsionalStiffness);
            expectNear("with the line tie off as well", openTwisted.torsionalStiffness,
                       wholeTwist.torsionalStiffness, 1e-9 * wholeTwist.torsionalStiffness);
        }
        if (sections == 8) {
            // Measured at -19.20% with the ties off. Asserted as a band rather than a
            // point because it is a property of this box's bay count, but tight enough
            // that a control which had stopped losing ties would fail it -- which is
            // the failure mode that would make the line above vacuous.
            expectTrue("eight sections lose about a fifth of the torsional stiffness",
                       openTwisted.torsionalStiffness < 0.85 * wholeTwist.torsionalStiffness &&
                           openTwisted.torsionalStiffness > 0.75 * wholeTwist.torsionalStiffness);
            expectNear("and the same eight sections tied lose 16 m of junction", open.tiedEdges,
                       whole.tiedEdges - 7 * 6.0, 1e-9);
        }
    }
}

// **A junction that is *off* the master's mid-surface, which is where the
// through-thickness half of a line tie lives -- and the box's own corners are not.**
//
// Every corner of a rectangular box is a butt joint: the two mid-surfaces meet on the
// corner line, so the offset along the master's normal is zero and the weight is 0.5
// whatever the plating is. Mutation testing is what said so — replacing the split with
// a flat 0.5, and taking the master's thickness from one end of the segment instead of
// interpolating it, both survived the whole suite. A real hull is not like that: a deck
// edge sits on the shell's *outer* face plus whatever gap two plates of different
// thickness leave, and the reference ferry's junctions run at `w = 1.69`.
//
// So: a box with a deck laid **inboard of the side plating**, stopping `kDeckGap` short
// of its mid-surface, at a height where it lands in the *middle* of a side element
// rather than on one of its nodes — and side plating that steps in thickness at
// mid-height, so the two ends of the master segment carry different thicknesses and
// interpolating them is not the same as taking either. Every number below is then a
// closed form in the gap and the two thicknesses.
constexpr double kDeckGap = 0.012;    // m inboard of the side plating's mid-surface
constexpr double kSideLow = 0.010;    // m, the side plating below mid-height
constexpr double kSideHigh = 0.020;   // m, and above it
constexpr double kDeckZ = kH / 4;     // mid-way along the lower side element

StructuralMesh makeBoxWithInboardDeck(double deckZ = kDeckZ) {
    StructuralMesh mesh;
    mesh.materials = {ah36Steel()};
    mesh.frameSpacing = kL / kNx;
    const auto quad = [&](const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d, double t) {
        PlatePanel p;
        p.corner[0] = a;
        p.corner[1] = b;
        p.corner[2] = c;
        p.corner[3] = d;
        p.thickness = t;
        p.material = 0;
        p.role = PanelRole::Shell;
        p.source = 0;
        mesh.panels.push_back(p);
    };
    for (int i = 0; i < kNx; ++i) {
        const double x0 = kL * i / kNx, x1 = kL * (i + 1) / kNx;
        for (int j = 0; j < kNy; ++j) {
            const double y0 = -kB / 2 + kB * j / kNy, y1 = -kB / 2 + kB * (j + 1) / kNy;
            quad({x0, y0, 0}, {x1, y0, 0}, {x1, y1, 0}, {x0, y1, 0}, kSideLow);
            quad({x0, y0, kH}, {x1, y0, kH}, {x1, y1, kH}, {x0, y1, kH}, kSideLow);
        }
        // The sides, thin below mid-height and thick above it. The node they share at
        // z = kH/2 therefore carries the area-weighted mean of the two, which is what
        // makes the master segment's thickness vary along its own length.
        for (int k = 0; k < kNz; ++k) {
            const double z0 = kH * k / kNz, z1 = kH * (k + 1) / kNz;
            const double t = k == 0 ? kSideLow : kSideHigh;
            quad({x0, -kB / 2, z0}, {x1, -kB / 2, z0}, {x1, -kB / 2, z1}, {x0, -kB / 2, z1}, t);
            quad({x0, kB / 2, z0}, {x1, kB / 2, z0}, {x1, kB / 2, z1}, {x0, kB / 2, z1}, t);
        }
        // The deck, stopping short of the sides on both hands.
        const double edge = kB / 2 - kDeckGap;
        quad({x0, -edge, deckZ}, {x1, -edge, deckZ}, {x1, 0, deckZ}, {x0, 0, deckZ}, kSideLow);
        quad({x0, 0, deckZ}, {x1, 0, deckZ}, {x1, edge, deckZ}, {x0, edge, deckZ}, kSideLow);
    }
    for (int i = 0; i <= kNx; ++i) mesh.frameStations.push_back(kL * i / kNx);
    return mesh;
}

void testPlaneTieSplitsTheMastersThickness() {
    std::printf("\n--- section: a line tie through a master of two thicknesses ---\n");
    const section::Section tied =
        section::buildSection(makeBoxWithInboardDeck(), boxParams());
    expectTrue("the deck box meshed", !tied.empty());
    expectTrue("and tied its cut planes", tied.planeTieNodes > 0);
    expectEqual("with nothing left open on them", tied.junctionsOnInterface, 0);

    // **The closed form, and getting it wrong first is the useful part.** The obvious
    // answer -- the master's thickness at the deck's own height, half way along the
    // side's lower element -- is 1.46 and the mesher says 1.463855. The difference is
    // that a tie is made at the slave's two *extruded* nodes and not at its
    // mid-surface: they sit half the deck's own thickness above and below it, so they
    // land at `s = -0.02` and `+0.02` along the segment rather than at its middle, and
    // the master's interpolated thickness is a shade different at each.
    //
    // So: the side's node at z = 0 is reached only by the thin band and carries
    // `kSideLow`; the node at z = kH/2 is reached by both and carries their mean; the
    // thickness anywhere between is the linear interpolation of those two; and
    // `constraint::tieWeight` puts the deck `kDeckGap` off that mid-surface along the
    // side's normal. `worstPlaneTieWeight` is `max(|w|, |1 - w|)` over both faces,
    // which comes out the same number whichever way round the side's normal points and
    // is set by the *thinner* of the two.
    const double element = kH / kNz;
    const double atZero = kSideLow, atMid = 0.5 * (kSideLow + kSideHigh);
    const auto masterThickness = [&](double z) {
        const double s = 2 * (z / element) - 1;
        return 0.5 * (1 - s) * atZero + 0.5 * (1 + s) * atMid;
    };
    const double lower = masterThickness(kDeckZ - 0.5 * kSideLow);
    const double upper = masterThickness(kDeckZ + 0.5 * kSideLow);
    const double expected = kDeckGap / std::min(lower, upper) + 0.5;
    std::printf("     master thickness %.6f / %.6f m at the two extruded nodes, over a"
                " %.3f/%.3f m step: worst weight %.6f against a closed form of %.6f\n",
                lower, upper, kSideLow, kSideHigh, tied.worstPlaneTieWeight, expected);
    // Exact to rounding: every term is a mean of doubles that are exactly representable
    // multiples of a millimetre, so there is nothing here to converge.
    expectNear("the line tie splits the master's own interpolated thickness",
               tied.worstPlaneTieWeight, expected, 1e-12);
    // **Three vacuity guards, one per way this could have passed while proving
    // nothing.** The measurement has to be away from the mid-surface at all; it has to
    // differ from what either end of the master segment alone would give, because
    // taking one end instead of interpolating is the plausible wrong implementation;
    // and the box's own butt corners -- which are what every other junction test here
    // measures -- have to be the 0.5 that says they could not have produced it.
    expectTrue("which is not the mid-surface", std::abs(expected - 0.5) > 0.5);
    expectTrue("and not what either end of the segment alone would give",
               std::abs(expected - (kDeckGap / atZero + 0.5)) > 0.1 &&
                   std::abs(expected - (kDeckGap / atMid + 0.5)) > 0.1);
    // And the two faces really did land at different points along the segment, which is
    // the whole reason the figure is 1.4639 and not 1.46.
    expectTrue("the slave's two extruded nodes land at different thicknesses",
               std::abs(upper - lower) > 1e-5);
    const section::Section butt = section::buildSection(makeBox(kT, kT), boxParams());
    expectNear("where a butt corner has no offset to split and gives exactly the mid-surface",
               butt.worstPlaneTieWeight, 0.5, 1e-15);
    // And the deck really is a junction the section could not weld, rather than plating
    // that happens to touch: it is a separate surface, and the tie is what joins it.
    expectTrue("the deck is a surface of its own", tied.surfaces > 4);
    expectNear("and the line drops nothing on prismatic plating", tied.worstPlaneTieSlip, 0.0,
               1e-15);
}

// **A chain whose line ties actually split a master's thickness, solved.**
//
// Every junction on the plain box is a butt joint at `w = 1/2`, where the two masters
// of a through-thickness pair carry the same weight and swapping them is a no-op --
// mutation testing said so, by surviving that swap on the whole suite. The inboard
// deck is at `w = 1.464`, so `-0.464` on the master's other face, and now the two
// halves are worth different things.
//
// The reference is the same length in one piece, where the same junction is closed by
// a *face* tie instead. That the two agree exactly is a claim in its own right: a
// slave whose own normal is square to the cut lands on the master face's edge, where
// the bilinear shape functions collapse onto the two nodes the line tie uses, so the
// line tie is the face tie restricted rather than an approximation to it.
void testDeckChainSplitsTheSameThicknessTheMonolithDoes() {
    std::printf("\n--- section: a chain of the deck box, where the split is not a half ---\n");
    const StructuralMesh structure = makeBoxWithInboardDeck();
    const StructuralMaterial steel = ah36Steel();
    const section::Section whole = section::buildSection(structure, boxParams());
    const section::TorsionResponse wholeTwist = section::applyTwist(whole, steel, 1e-4, kH / 2);
    expectTrue("the monolith solved", wholeTwist.ok);
    expectTrue("and its junctions really do extrapolate through the master",
               whole.worstPlaneTieWeight > 1.4 && whole.worstJunctionWeight > 1.4);

    section::ChainParams params = boxChain(4, true);
    section::ChainParams open = params;
    open.section.interfaceTies = false;
    const section::Chain chain = section::buildChain(structure, params);
    const section::Chain cut = section::buildChain(structure, open);
    expectTrue("both chains built", chain.ready() && cut.ready());
    const section::TorsionResponse twisted = section::applyTwist(chain, 1e-4, kH / 2);
    const section::TorsionResponse cutTwist = section::applyTwist(cut, 1e-4, kH / 2);
    expectTrue("both solved", twisted.ok && cutTwist.ok);
    std::printf("     deck box, 4 sections: tied %.1f m of %.1f (one piece %.1f), GJ %+.2e"
                "  |  cut planes open: %.1f m, GJ %+.3f%%\n",
                chain.tiedEdges, chain.junctionEdges, whole.tiedEdges,
                twisted.torsionalStiffness / wholeTwist.torsionalStiffness - 1, cut.tiedEdges,
                100.0 * (cutTwist.torsionalStiffness / wholeTwist.torsionalStiffness - 1));
    expectNear("the chain joins what the monolith joins", chain.tiedEdges, whole.tiedEdges, 1e-9);
    // The face tie and the line tie are the same constraint here, so this is exact to
    // the conditioning of two independent solves rather than converged.
    expectNear("and carries the torsion the monolith carries", twisted.torsionalStiffness,
               wholeTwist.torsionalStiffness, 1e-9 * wholeTwist.torsionalStiffness);
    // The guard, without which the line above would pass on a chain that had never lost
    // anything: with the cut planes left open the same model is visibly softer.
    expectTrue("where leaving the cut planes open costs it torsional stiffness",
               cutTwist.torsionalStiffness < 0.97 * wholeTwist.torsionalStiffness);
    expectTrue("and leaves junction edge unjoined", cut.tiedEdges < whole.tiedEdges - 1.0);
}

// **The three ways a line tie is refused, each on an input whose verdict it changes.**
//
// Mutation testing is what said these needed writing: deleting the overshoot bound,
// deleting the weight bound and deleting the one-segment rule all survived the whole
// suite, because every junction on the box and on the ferry is comfortably inside all
// three. A predicate is only tested by an input it says no to.
void testPlaneTieRefusals() {
    std::printf("\n--- section: what a line tie on a cut plane refuses ---\n");

    // **1. A slave that straddles the node between two segments.** Put the deck at
    // exactly mid-height, where the side plating's two thickness bands meet: the deck's
    // pair straddles that node, so its lower half lands on the segment below and its
    // upper half on the one above. Two segments would split the slave against two
    // plates -- and against two thicknesses here -- so it is refused whole.
    const section::Section straddled =
        section::buildSection(makeBoxWithInboardDeck(kH / 2), boxParams());
    const section::Section clear = section::buildSection(makeBoxWithInboardDeck(), boxParams());
    std::printf("     deck on the node between two segments: %d tied, %d unreached;"
                " mid-element: %d tied, %d unreached\n",
                straddled.planeTieNodes, straddled.planeTiesUnreached, clear.planeTieNodes,
                clear.planeTiesUnreached);
    expectTrue("the deck on a segment node is refused", straddled.planeTiesUnreached > 0);
    // The vacuity guard: the same deck a quarter of the way up lands inside one segment
    // and is tied, so the refusal is about where it landed and not about the fixture.
    expectEqual("where the same deck inside a segment is not", clear.planeTiesUnreached, 0);
    expectTrue("and is tied instead", clear.planeTieNodes > straddled.planeTieNodes);
    // A refused junction is still counted as open rather than forgotten.
    expectTrue("and a refused one is counted open", straddled.junctionsOnInterface > 0);

    // **2. A face overshoot past `junctionOvershoot`.** The box's butt corners put the
    // flange's node half a plate thickness past the end of the side's mid-surface,
    // which is `t / 2` over half an element -- 0.02 of the segment's own coordinate,
    // a closed form and not a reading. Bound it below that and every line tie goes.
    const double overshoot = (kT / 2) / (0.5 * (kH / kNz));
    const section::Section tied = section::buildSection(makeBox(kT, kT), boxParams());
    expectNear("the box's line ties overshoot by half a thickness over half an element",
               tied.worstPlaneTieOvershoot, overshoot, 1e-12);
    section::SectionParams tight = boxParams();
    tight.junctionOvershoot = 0.5 * overshoot;
    const section::Section bounded = section::buildSection(makeBox(kT, kT), tight);
    expectEqual("bounding the overshoot below that refuses every one of them",
                bounded.planeTieNodes, 0);
    // Sixteen and not eight: with no tie made, the node that would have been the
    // *master* half of each corner is no longer absorbed by one, so it comes round in
    // its own turn and is refused on the same ground. Four corners, two planes, two
    // plates at each.
    expectEqual("and counts every one of them", bounded.planeTiesOutsideLine, 4 * 2 * 2);
    expectEqual("which is what the section then reports as still open",
                bounded.junctionsOnInterface, 4 * 2 * 2);
    expectTrue("where the default admits them all", tied.planeTieNodes > 0);
    expectEqual("and refuses none", tied.planeTiesOutsideLine, 0);

    // **3. A through-thickness split past `junctionWeightLimit`.** The inboard deck
    // asks for `w = 1.464`, so `-0.464` on the master's other face -- inside the
    // default one full share and outside a fifth of one.
    section::SectionParams mean = boxParams();
    mean.junctionWeightLimit = 0.2;
    const section::Section refused = section::buildSection(makeBoxWithInboardDeck(), mean);
    std::printf("     deck weight %.4f: %d tied at a limit of 1.0, %d at 0.2 (%d over weight)\n",
                clear.worstPlaneTieWeight, clear.planeTieNodes, refused.planeTieNodes,
                refused.planeTiesThroughThickness);
    expectTrue("the deck's split is past a fifth of a share",
               clear.worstPlaneTieWeight - 1.0 > 0.2);
    expectTrue("so a limit of a fifth refuses some of them",
               refused.planeTiesThroughThickness > 0);
    expectTrue("and ties fewer", refused.planeTieNodes < clear.planeTieNodes);
    expectEqual("where the default refuses none", clear.planeTiesThroughThickness, 0);
    // And the guard that the limit is about the *split* and not about everything: the
    // box's own butt corners are at 0.5 and survive the same limit untouched.
    const section::Section buttAtLimit = section::buildSection(makeBox(kT, kT), mean);
    expectEqual("a butt corner is unaffected by the same limit", buttAtLimit.planeTieNodes,
                tied.planeTieNodes);
}

// **The same station named by two sub-quads that wind opposite ways.**
//
// `makeStructuralMesh` reverses every other bay's corner order when it mirrors the
// starboard side, so the sub-quad aft of a cut plane and the one forward of it can name
// the segment they share in opposite directions. Both sections would then measure `s`
// from opposite ends of the same line and interpolate its masters backwards for one of
// them — which is why the segment is ordered by *position* before anything is measured
// along it. Mutation testing is what said this needed a test: deleting that ordering
// survived every fixture wound one way.
void testPlaneTieSurvivesReversedWinding() {
    std::printf("\n--- section: a cut plane whose two sides wind opposite ways ---\n");
    const StructuralMesh reversed = makeBox(kT, kT, /*alternateWinding=*/true);
    const StructuralMaterial steel = ah36Steel();
    const section::Section whole = section::buildSection(reversed, boxParams());
    const section::TorsionResponse wholeTwist = section::applyTwist(whole, steel, 1e-4, kH / 2);
    expectTrue("the monolith solved", wholeTwist.ok);

    // Four sections, so three interior planes, each of which has a bay of one winding
    // on one side and a bay of the other on the other.
    const section::Chain chain = section::buildChain(reversed, boxChain(4, true));
    expectTrue("the chain built", chain.ready());
    const section::TorsionResponse twisted = section::applyTwist(chain, 1e-4, kH / 2);
    expectTrue("and solved", twisted.ok);
    std::printf("     alternating winding, 4 sections: %d line ties, %d planes disagreed,"
                " tied %.1f m, GJ %+.2e of the monolith's\n",
                chain.planeTieNodes, chain.planeTiesDisagreeing, chain.tiedEdges,
                twisted.torsionalStiffness / wholeTwist.torsionalStiffness - 1);
    expectEqual("every interior plane is tied", chain.planeTieNodes, 4 * 3);
    expectEqual("and the two sides of every one of them agree", chain.planeTiesDisagreeing, 0);
    expectNear("to the last bit", chain.worstPlaneTieDisagreement, 0.0, 0.0);
    expectNear("so the chain joins what the monolith joins", chain.tiedEdges, whole.tiedEdges,
               1e-9);
    expectNear("and carries its torsion", twisted.torsionalStiffness,
               wholeTwist.torsionalStiffness, 1e-9 * wholeTwist.torsionalStiffness);
    // The guard that the winding was actually alternating: the same box wound one way
    // is a different input, and a fixture that had quietly stopped reversing would make
    // this test a duplicate of the one above it.
    expectTrue("and the fixture really does alternate", [] {
        const StructuralMesh a = makeBox(kT, kT, true), b = makeBox(kT, kT, false);
        for (std::size_t p = 0; p < a.panels.size() && p < b.panels.size(); ++p)
            for (int c = 0; c < 4; ++c)
                if (length(a.panels[p].corner[c] - b.panels[p].corner[c]) > 0) return true;
        return false;
    }());
}

// **The constraint on the assembled model, checked as arithmetic rather than through
// a solve.**
//
// `GJ` coming back to the monolith's is the finding, and it is an aggregate: a fold
// that dropped one of the four terms of `TᵀKT` would move it by a little and look like
// conditioning. The recurring shape of a defect in this repository is exactly that --
// an error that cancels when asked globally -- so the fold is asked about alone.
//
// A **rigid translation of the whole chain stores no energy**, and it has to survive
// the constraint: the tie's weights are a partition of unity, so the slave moves with
// its masters, and `TᵀKT` applied to the translation with the eliminated entries at
// zero must give back zero force. That dies on a weight that does not sum to one, on a
// master named wrongly, and on the single-pass fold that drops the `w_a w_b K[s][s]`
// term -- none of which any energy comparison separates.
void testChainPlaneTiesAreExactOnTheAssembly() {
    std::printf("\n--- section: a rigid translation of a tied chain stores nothing ---\n");
    const StructuralMesh structure = makeBox(kT, kT);
    const section::Chain chain = section::buildChain(structure, boxChain(4, true));
    expectTrue("the chain built", chain.ready());
    expectTrue("and it applied some line ties", !chain.planeTies.empty());
    expectEqualCount("six constraints per tied node", chain.planeTies.size(),
                     static_cast<std::size_t>(6 * chain.planeTieNodes));
    expectEqualCount("and one eliminated row each", chain.planeTieDof.size(),
                     chain.planeTies.size());

    const auto n = static_cast<std::size_t>(chain.assembly.size());
    expectEqualCount("the assembly is square", chain.assembly.stiffness.size(), n * n);
    // The fold's postcondition: an eliminated row and column carry nothing but a unit
    // diagonal, so holding the row costs the model nothing. A fold that left the slave
    // coupled would double-count it, and the energy would still look plausible.
    double worstLeak = 0, worstDiagonal = 0;
    for (std::uint32_t d : chain.planeTieDof)
        for (std::size_t j = 0; j < n; ++j) {
            if (j == d) {
                worstDiagonal = std::max(worstDiagonal,
                                         std::abs(chain.assembly.stiffness[d * n + j] - 1.0));
                continue;
            }
            worstLeak = std::max(worstLeak, std::abs(chain.assembly.stiffness[d * n + j]));
            worstLeak = std::max(worstLeak, std::abs(chain.assembly.stiffness[j * n + d]));
        }
    expectNear("an eliminated row and column are exactly empty", worstLeak, 0.0, 0.0);
    expectNear("with a unit diagonal", worstDiagonal, 0.0, 0.0);

    // A rigid translation along each axis, written on the boundary rows the constraint
    // did *not* eliminate. The slaves stay at zero, which is what `TᵀKT` expects: the
    // masters carry them, and the weights being a partition of unity is what makes the
    // slave arrive at the same 1 m as everything else.
    std::vector<std::uint8_t> eliminated(n, 0u);
    for (std::uint32_t d : chain.planeTieDof) eliminated[d] = 1u;
    const std::vector<reduction::BoundaryDof>& point = chain.assembly.boundaryPoint;
    expectEqualCount("the assembly carries its boundary identity", point.size(),
                     static_cast<std::size_t>(chain.assembly.boundary));
    double scale = 0;
    for (std::size_t i = 0; i < n; ++i)
        scale = std::max(scale, std::abs(chain.assembly.stiffness[i * n + i]));
    expectTrue("the assembled stiffness is not empty", scale > 0);

    double worstForce = 0, worstEnergy = 0;
    for (int axis = 0; axis < 3; ++axis) {
        std::vector<double> u(n, 0.0);
        for (int b = 0; b < chain.assembly.boundary; ++b) {
            const auto d = static_cast<std::size_t>(b);
            if (point[d].axis == static_cast<std::uint32_t>(axis) && !eliminated[d]) u[d] = 1.0;
        }
        double energy = 0;
        for (std::size_t i = 0; i < n; ++i) {
            double force = 0;
            for (std::size_t j = 0; j < n; ++j) force += chain.assembly.stiffness[i * n + j] * u[j];
            energy += 0.5 * u[i] * force;
            if (!eliminated[i]) worstForce = std::max(worstForce, std::abs(force));
        }
        worstEnergy = std::max(worstEnergy, std::abs(energy));
    }
    std::printf("     rigid translation of the tied chain: worst force %.3e N of a %.3e N/m"
                " diagonal, energy %.3e J\n",
                worstForce, scale, worstEnergy);
    // Measured at 1.2e-1 N against a 2.27e+12 N/m diagonal -- **5.3e-14 relative**,
    // which is the assembly's own conditioning on a dense Guyan model this size.
    // Asserted at 1e-12 of the diagonal: nineteen times the measurement, which is the
    // margin a different optimisation level wants, and many orders below what dropping
    // one of the four terms of the fold would cost.
    expectTrue("a rigid translation of the tied chain carries no force",
               worstForce < 1e-12 * scale);
    expectTrue("and stores no energy", worstEnergy < 1e-12 * scale);

    // **The vacuity guard.** A zero matrix would pass everything above, and so would a
    // model whose boundary rows had all been eliminated. A field that is *not* rigid
    // has to cost something, on the same matrix, through the same product.
    std::vector<double> stretched(n, 0.0);
    for (int b = 0; b < chain.assembly.boundary; ++b) {
        const auto d = static_cast<std::size_t>(b);
        if (point[d].axis == 0 && !eliminated[d]) stretched[d] = point[d].position.x;
    }
    double stretchEnergy = 0;
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j)
            stretchEnergy += 0.5 * stretched[i] * chain.assembly.stiffness[i * n + j] * stretched[j];
    std::printf("     against %.3e J for stretching it, on the same matrix\n", stretchEnergy);
    // 5.2e+10 J. The guard exists because a zero matrix, or one whose boundary rows had
    // all been eliminated, would pass every assertion above.
    expectTrue("where stretching the same chain costs a great deal", stretchEnergy > 1e6);
}

// **The other quantity that sees a junction: the lowest fixed-interface frequency.**
//
// `EA` and `EI` are exact on a mesh whose plates are not joined at all (§2), so the
// two instruments that can tell a joined chain from an open one are torsion and this.
// The reference is the same length in one piece with the same plating tied, through
// none of the assembly -- so the chain owes it nothing.
void testChainFrequencyReachesTheMonolith() {
    std::printf("\n--- section: a tied chain's first frequency against the monolith's ---\n");
    const StructuralMesh structure = makeBox(kT, kT);
    const section::Section whole = section::buildSection(structure, boxParams());
    reduction::Substructure sWhole(whole.mesh, whole.material, whole.interfaceNodes,
                                   whole.attachment);
    expectTrue("the monolithic substructure is ready", sWhole.ready());
    const reduction::Eigenpairs exact = sWhole.fixedInterfaceModes(1);
    expectTrue("and has a first mode", !exact.value.empty());
    const double reference = std::sqrt(std::max(0.0, exact.value[0]));
    expectTrue("which is not zero", reference > 0);

    for (int sections : {2, 4}) {
        section::ChainParams tiedParams = boxChain(sections, true);
        tiedParams.reduce.modes = 6;
        section::ChainParams openParams = tiedParams;
        openParams.section.interfaceTies = false;
        const section::Chain tied = section::buildChain(structure, tiedParams);
        const section::Chain open = section::buildChain(structure, openParams);
        expectTrue("both chains built", tied.ready() && open.ready());
        const std::vector<double> tiedOmega = section::chainFrequencies(tied);
        const std::vector<double> openOmega = section::chainFrequencies(open);
        expectTrue("both have a spectrum", !tiedOmega.empty() && !openOmega.empty());
        std::printf("     %d sections at 6 modes: %.6f Hz tied, %.6f Hz with the cut planes open,"
                    " monolith %.6f Hz\n",
                    sections, tiedOmega[0] / (2 * std::numbers::pi),
                    openOmega[0] / (2 * std::numbers::pi), reference / (2 * std::numbers::pi));

        // **The two sit on opposite sides of the monolith, and that is the whole
        // discrimination -- no tolerance can fudge a sign.** A tied chain is the same
        // structure as the monolith with each piece reduced, so its Rayleigh quotient
        // is an upper bound and it approaches from above. A chain whose cut planes are
        // open is a *different, softer* structure -- rings of junctions carrying no
        // shear -- so it falls through the monolith rather than converging on it, and
        // no number of modes brings it back.
        expectTrue("an assembled frequency is an upper bound on the true one",
                   tiedOmega[0] > reference * (1.0 - 1e-9));
        expectTrue("where an open cut plane makes it a softer structure, and it falls below",
                   openOmega[0] < reference * (1.0 - 1e-3));
        // Measured 1.9e-3 above at N = 2 and 7e-6 at N = 4 -- the modal truncation,
        // which four sections carry twice as many modes against. The open chain is
        // 2.3e-3 and 6.2e-3 *below*, so the two are separated by their sign and not by
        // the width of this tolerance.
        expectNear("a tied chain reaches the monolith's own first frequency", tiedOmega[0],
                   reference, 5e-3 * reference);
        expectTrue("and it is stiffer than the open one at every N", tiedOmega[0] > openOmega[0]);
    }
}

// **The lowest fixed-interface frequency, which is the instrument that sees whether
// the sections are joined.**
//
// Both end planes held is exactly the boundary condition
// `Substructure::fixedInterfaceModes` applies to the same length in one piece, so the
// monolithic answer is an independent reference. Three things are checked against it,
// and the two negative controls are what make the first mean anything:
//
//   * untied on both sides, the chain converges on the monolith's own frequency from
//     above as modes are added -- which is the reduction's upper-bound property
//     surviving assembly;
//   * a chain whose plating is not tied reads a fourteenth of the tied frequency,
//     which is the failure §2 of this file exists to catch;
//   * a chain whose *interfaces* did not match leaves its interior sections held at
//     nothing, and opens with six rigid body modes per loose section instead of an
//     elastic one. Neither shows up in `EA`.
void testChainFrequencySeesTheJoins() {
    std::printf("\n--- section: the lowest frequency of a chain, and what it sees ---\n");
    const StructuralMesh structure = makeBox(kT, kT);

    // The reference: the same length in one piece, its own fixed-interface modes,
    // through none of the assembly.
    const section::Section whole = section::buildSection(structure, openBoxParams());
    reduction::Substructure sWhole(whole.mesh, whole.material, whole.interfaceNodes,
                                   whole.attachment);
    expectTrue("the monolithic substructure is ready", sWhole.ready());
    const reduction::Eigenpairs exact = sWhole.fixedInterfaceModes(1);
    expectTrue("and has a first mode", !exact.value.empty());
    const double reference = std::sqrt(std::max(0.0, exact.value[0]));
    expectTrue("which is not zero", reference > 0);

    double previous = 0, atZeroModes = 0, converged = 0;
    for (int modes : {0, 2, 6}) {
        section::ChainParams params = boxChain(4, false);
        params.reduce.modes = modes;
        const section::Chain chain = section::buildChain(structure, params);
        expectTrue("the chain built", chain.ready());
        const std::vector<double> omega = section::chainFrequencies(chain);
        expectTrue("it has a spectrum", !omega.empty());
        std::printf("     %d modes a section: %.9f Hz against the monolith's %.9f\n", modes,
                    omega[0] / (2 * std::numbers::pi), reference / (2 * std::numbers::pi));
        expectTrue("an assembled frequency is an upper bound on the true one",
                   omega[0] > reference * (1.0 - 1e-9));
        if (modes == 0) atZeroModes = omega[0];
        // Non-increasing to 1e-7. Enlarging each component's Ritz space can only
        // lower an assembled Rayleigh quotient, so this is a theorem rather than a
        // hope -- but by six modes the sequence has stopped moving and what is left
        // is the dense generalised eigensolve's own accuracy, measured at 1.7e-9
        // relative between the 2- and 6-mode runs. The slack is a hundred times that
        // and forty times below the truncation the sequence is still carrying, so a
        // chain that had stopped converging would still fail this.
        if (modes > 0)
            expectTrue("and adding modes brings it down", omega[0] <= previous * (1.0 + 1e-7));
        previous = omega[0];
        converged = omega[0];
    }
    // Measured: 0.844705 Hz at zero modes falling to 0.843943, which is the monolith's
    // own 0.843941 to six figures. Asserted at 1e-5 relative because that is what was
    // measured; 1e-2 would pass on a chain that had lost a whole section.
    expectNear("with modes the chain reaches the monolith's own frequency", converged, reference,
               1e-5 * reference);
    // And the guard that there was anything to converge: Guyan alone is visibly above
    // it, so this is a convergence rather than an identity that would hold on any
    // assembly at all.
    expectTrue("zero modes was visibly above the monolith's own frequency",
               atZeroModes > reference * (1.0 + 1e-4));

    // Negative control 1: the plating not tied. A chain of untied sections is four
    // loose plates and reads a plate's frequency, not a box's.
    const section::Chain tied = section::buildChain(structure, boxChain(4, true));
    const section::Chain loosePlating = section::buildChain(structure, boxChain(4, false));
    const std::vector<double> tiedOmega = section::chainFrequencies(tied);
    const std::vector<double> looseOmega = section::chainFrequencies(loosePlating);
    expectTrue("both chains have a spectrum", !tiedOmega.empty() && !looseOmega.empty());
    std::printf("     tied %.4f Hz against untied %.4f Hz (x%.1f)\n",
                tiedOmega[0] / (2 * std::numbers::pi), looseOmega[0] / (2 * std::numbers::pi),
                tiedOmega[0] / looseOmega[0]);
    expectTrue("tying the plating moves the chain's first frequency by an order of magnitude",
               tiedOmega[0] > 10.0 * looseOmega[0]);

    // Negative control 2: the sections not joined to *each other*. The same four
    // components, assembled with no joints at all -- which is four perfectly good
    // reduced models sitting in one matrix and produces a full set of plausible
    // frequencies. What gives it away is that the two interior sections are then held
    // at nothing and float free, six rigid body modes each.
    {
        std::vector<reduction::Component> parts(tied.section.size());
        for (std::size_t i = 0; i < parts.size(); ++i)
            parts[i] = {&tied.substructure[i], &tied.reduced[i]};
        const reduction::Assembly apart = reduction::assemble(parts, {});
        expectEqual("unjoined, the assembled model is four pieces",
                    reduction::assembledComponents(apart), 4);
        expectEqualCount("and its boundary is the sum rather than the union",
                         static_cast<std::size_t>(apart.boundary),
                         static_cast<std::size_t>(4 * tied.reduced[0].boundary));

        // Held exactly as `chainFrequencies` holds them: the two outermost planes,
        // named off the mesher's own plane lists rather than by position.
        std::vector<std::uint32_t> held;
        const auto holdPlane = [&](std::size_t i, const std::vector<std::uint32_t>& nodes) {
            std::vector<std::uint8_t> onPlane(tied.section[i].mesh.nodeCount(), 0u);
            for (std::uint32_t n : nodes) onPlane[n] = 1u;
            const std::vector<std::uint32_t>& dof = tied.substructure[i].boundaryDof();
            for (std::size_t b = 0; b < dof.size(); ++b)
                if (onPlane[dof[b] / 3u])
                    held.push_back(static_cast<std::uint32_t>(apart.from[i][b]));
        };
        holdPlane(0, tied.section.front().aftNodes);
        holdPlane(parts.size() - 1, tied.section.back().forwardNodes);
        expectEqualCount("two whole planes are held", held.size(),
                         3 * (tied.section.front().aftNodes.size() +
                              tied.section.back().forwardNodes.size()));

        const std::vector<double> loose = reduction::assembledFrequencies(apart, held);
        expectTrue("the unjoined model has a spectrum", loose.size() > 13);
        // The cutoff is scaled off an independently derived frequency rather than
        // fixed: rigid modes are near zero but not at machine precision, and a
        // constant landing between the translations and the rotations is a mistake
        // this repository has now made twice.
        // **Count the rigid modes by the gap, not by a threshold.** A rigid mode is
        // near zero and not at machine precision -- these run up to 7e-3 rad/s -- and
        // a fixed cutoff landing between the translations and the rotations is a
        // mistake this repository has recorded twice. So the split is the largest
        // ratio between consecutive frequencies in the bottom of the spectrum, which
        // is derived from the answer rather than assumed about it.
        std::size_t split = 1;
        double widest = 0;
        for (std::size_t i = 1; i < 20 && i < loose.size(); ++i)
            if (loose[i - 1] > 0 && loose[i] / loose[i - 1] > widest) {
                widest = loose[i] / loose[i - 1];
                split = i;
            }
        std::printf("     unjoined: %zu modes below a gap of x%.0f, first elastic %.4f Hz\n",
                    split, widest, loose[split] / (2 * std::numbers::pi));
        expectTrue("the rigid and elastic parts of the spectrum are orders of magnitude apart",
                   widest > 100.0);
        expectEqualCount("two of the four sections float free and bring six rigid modes each",
                         split, 12u);
        // And the elastic part of an unjoined chain is *higher*, not lower, which is
        // the other half of why this control is worth having: what is left after the
        // rigid modes is the two end sections held at one plane each, a quarter of the
        // span, and a quarter-length cantilever is stiff. Measured 6.31 Hz against the
        // joined chain's 0.84 -- so a chain that failed to join would show a spectrum
        // that is wrong in both directions at once and plausible in neither.
        expectTrue("what is left is two quarter-length pieces, and they are stiffer",
                   loose[split] > 3.0 * looseOmega[0]);
    }
}

// --- The ferry: the same identities on real geometry ------------------------------
//
// Mesh level only. A chain at ship scale is a dense reduced model of
// `(N+1) x 1 170` DOF and the Cholesky at the end of it is core-seconds, so the
// solved comparison lives in `tools/section_probe --chain=N`. What is here is what a
// solve would never notice: whether cutting a length in two changes how much steel
// there is.
void testFerryChainConservesTheStructure() {
    std::printf("\n--- section: two ferry sections carry the steel one section carries ---\n");
    const StructuralMesh structure = ferryStructure();

    section::SectionParams params;
    params.subdivision = 1;
    params.junctions = false;
    const auto cut = [&](double from, double to) {
        section::SectionParams p = params;
        p.xFrom = from;
        p.xTo = to;
        return section::buildSection(structure, p);
    };
    const section::Section aft = cut(-7.2, -2.4);
    const section::Section forward = cut(-2.4, 2.4);
    const section::Section whole = cut(-7.2, 2.4);
    expectTrue("all three meshed", !aft.empty() && !forward.empty() && !whole.empty());

    expectEqualCount("the two sections have the elements the one section has",
                     aft.elementCount() + forward.elementCount(), whole.elementCount());
    expectNear("and its plating, to the last bit", aft.plateMass + forward.plateMass,
               whole.plateMass, 1e-9);

    // **The one that was wrong.** A member with no extent along x -- a frame, a deck
    // beam -- lies *on* a station, and a station is a cut plane for the sections on
    // both sides of it. Taken by both, a chain carries that ring of frames twice: 608
    // fibres and 10.1% of two sections' stiffener mass, landing on shared interface
    // DOF where nothing but a total would see it. The half-open rule the transverse
    // *plating* already used is what fixes it, and this is the identity that says so.
    expectNear("and its stiffener steel", aft.memberMass + forward.memberMass, whole.memberMass,
               1e-9 * whole.memberMass);
    expectEqualCount("fibre for fibre", aft.stiffening.fiber.size() + forward.stiffening.fiber.size(),
                     whole.stiffening.fiber.size());
    expectEqual("member for member", aft.membersAttached + forward.membersAttached,
                whole.membersAttached);
    // The guard: the rule has to have bitten, or the identity above is the identity
    // of two things neither of which has a frame on the shared plane at all.
    expectTrue("the shared cut plane really does carry members",
               aft.membersOnForwardPlane > 0 && aft.membersOnForwardPlane == whole.membersOnForwardPlane);
    expectTrue("and the stiffener is a real fraction of the section",
               whole.memberMass > 0.3 * whole.plateMass);

    // And the interfaces this chain would be assembled through: the aft section's
    // forward plane is the forward section's aft plane, node for node.
    expectEqualCount("the shared plane has the same nodes on both sides",
                     aft.forwardNodes.size(), forward.aftNodes.size());
    reduction::Substructure sAft(aft.mesh, aft.material, aft.interfaceNodes, aft.attachment);
    reduction::Substructure sForward(forward.mesh, forward.material, forward.interfaceNodes,
                                     forward.attachment);
    expectTrue("both reduce", sAft.ready() && sForward.ready());
    const reduction::InterfaceMap map = reduction::matchBoundaries(sAft, sForward);
    expectEqualCount("and amidships every one of them matches at the default tolerance",
                     map.shared, 3 * aft.forwardNodes.size());
    expectNear("exactly", map.worstGap, 0.0, 0.0);
}

// --- The reach: the ends of the ship, which is what §7 is about ---------------------
//
// Before the collapsed element was understood, 21 of this ship's 49 two-bay windows
// meshed, reduced and solved -- five islands covering 62.4 m of 120 m, whose longest
// unbroken run was the 26.4 m hold every other measurement in this file is made on.
// That made every downstream result -- the chain, the N-way assembly, the hull-girder
// agreement, the torsion -- a statement about the parallel middle body only.
//
// The stations below are the ones that could not be meshed: the extreme stern, the
// extreme bow, a bulkhead amidships, and the shoulder. What is asserted is not
// "it meshed" -- an empty section meshes -- but that it reduces, solves, holds
// together, and carries the steel an independent count of the same panels says it
// should.
void testMesherReachesTheEndsOfTheShip() {
    std::printf("\n--- section: the bow, the stern and the bulkheads, which used to refuse ---\n");
    const StructuralMesh structure = ferryStructure();
    const StructuralMaterial material = ah36Steel();

    // The plating an independent count of the same panels finds in a window. This
    // is `buildSection`'s ownership rule -- a panel with extent along x belongs when
    // it lies wholly between the planes, a transverse plate is taken at the aft
    // plane only -- applied to `structure.panels` directly, so it shares no code
    // with the mesher. A wedge that had been *dropped* rather than meshed would show
    // up here and nowhere else: `EA` cannot see 0.3% of the plating and neither can
    // a component count.
    const auto expectedPlating = [&](double from, double to, double& area, double& mass) {
        area = 0;
        mass = 0;
        for (const PlatePanel& p : structure.panels) {
            double lo = 1e300, hi = -1e300;
            for (int c = 0; c < 4; ++c) {
                lo = std::min(lo, p.corner[c].x);
                hi = std::max(hi, p.corner[c].x);
            }
            const bool transverse = hi - lo <= 1e-6;
            const bool owned = transverse ? (lo >= from - 1e-6 && lo < to - 1e-6)
                                          : (lo >= from - 1e-6 && hi <= to + 1e-6);
            if (!owned) continue;
            // A quad's area as two triangles, which is right for a triangle too --
            // one of the two is degenerate and contributes nothing.
            const double quad =
                0.5 * length(cross(p.corner[1] - p.corner[0], p.corner[2] - p.corner[0])) +
                0.5 * length(cross(p.corner[2] - p.corner[0], p.corner[3] - p.corner[0]));
            area += quad;
            mass += quad * p.thickness * material.density;
        }
    };

    struct Station {
        const char* label;
        double from, to;
    };
    const Station stations[] = {{"extreme stern", -60.0, -55.2},
                                {"stern shoulder", -45.6, -40.8},
                                {"bulkhead at -8", -12.0, -7.2},
                                {"bulkhead at 20", 16.8, 21.6},
                                {"bow shoulder", 40.8, 45.6},
                                {"extreme bow", 55.2, 60.0}};
    std::printf("     %-16s %8s %7s %7s %6s %11s %11s %10s\n", "", "elements", "collps", "invert",
                "comps", "minGaussJ", "plate kg", "GJ");
    int reached = 0, withWedges = 0;
    for (const Station& station : stations) {
        section::SectionParams params;
        params.xFrom = station.from;
        params.xTo = station.to;
        params.subdivision = 1;
        const section::Section piece = section::buildSection(structure, params);

        const section::BeamResponse stretched =
            section::applyBeamLoad(piece, material, section::BeamLoad{1e-6, 0, 0});
        const section::TorsionResponse twisted = section::applyTwist(piece, material, 1e-6, 6.86);
        const reduction::Substructure substructure(piece.mesh, piece.material,
                                                   piece.interfaceNodes, piece.attachment);
        std::printf("     %-16s %8zu %7d %7d %6d %11.3e %11.0f %10.3e\n", station.label,
                    piece.elementCount(), piece.collapsedElements, piece.invertedElements,
                    piece.components, piece.worstJacobian, piece.plateMass,
                    twisted.torsionalStiffness);

        expectTrue(std::string(station.label) + " meshes something", !piece.empty());
        expectEqual(std::string(station.label) + " has no inverted element",
                    piece.invertedElements, 0);
        expectTrue(std::string(station.label) + "'s quadrature is positive",
                   piece.worstJacobian > 0);
        expectEqual(std::string(station.label) + " cuts no panel in half",
                    piece.straddlingPanels, 0);
        expectEqual(std::string(station.label) + " has nothing floating free of the interface",
                    piece.floatingComponents, 0);
        expectTrue(std::string(station.label) + " takes an axial load: " + stretched.problem,
                   stretched.ok);
        expectTrue(std::string(station.label) + " twists: " + twisted.problem, twisted.ok);
        expectTrue(std::string(station.label) + " reduces", substructure.ready());
        expectTrue(std::string(station.label) + " carries a real torque",
                   twisted.torsionalStiffness > 1e11);

        // The independent count. Exact, not approximate: the mesher's mid-surface is
        // the panel and the two are the same integral.
        double area = 0, mass = 0;
        expectedPlating(station.from, station.to, area, mass);
        expectNear(std::string(station.label) + " meshes the plating those panels have",
                   piece.area, area, 1e-9 * area);
        expectTrue(std::string(station.label) + " is a real piece of ship", area > 100.0);

        // **The steel is exact against that count only without the halo, and saying
        // why is the point.** `A_e * mean(t at its corners)` telescopes back to
        // `sum A_p t_p` precisely when the nodal thickness is averaged over the
        // sub-quads the section meshes. The halo of §8 averages over one bay beyond
        // each plane as well, so a cut plane where the strake steps carries the
        // *mean* of the two -- §3's taper, now reaching the plane -- and this window
        // hands its neighbour whatever it takes. What replaces the identity is
        // conservation across the cut, which is exact and is asserted below.
        section::SectionParams unaveraged = params;
        unaveraged.halo = false;
        const section::Section owned = section::buildSection(structure, unaveraged);
        expectNear(std::string(station.label) + " carries the steel those panels have",
                   owned.plateMass, mass, 1e-9 * mass);
        // And with the halo it is off by the seam it moved across the plane, which is
        // parts per million: three orders of magnitude below the 0.3% of plating a
        // dropped wedge would be, which is what this check exists to catch.
        expectNear(std::string(station.label) + " carries it to a seam's redistribution",
                   piece.plateMass, mass, 3e-6 * mass);

        if (piece.collapsedElements > 0) ++withWedges;
        ++reached;
    }
    expectEqual("every station reached", reached, 6);
    // **The vacuity guard, and it is the load-bearing one.** If none of these windows
    // contained a collapsed element then none of them was ever refused, and the whole
    // test is a statement about stations that always worked.
    expectEqual("and every one of them contains the shape that used to refuse", withWedges, 6);

    // The middle body is *not* one of them, which is what says the new path cannot
    // have moved it. Amidships there is no degenerate panel at all, so
    // `testFerryChainConservesTheStructure` and the chain tests below are measuring
    // exactly the mesh they measured before.
    section::SectionParams amidships;
    amidships.xFrom = -2.4;
    amidships.xTo = 2.4;
    const section::Section middle = section::buildSection(structure, amidships);
    expectEqual("amidships there is no collapsed element to be affected by any of this",
                middle.collapsedElements, 0);
    expectEqual("nor an inverted one", middle.invertedElements, 0);

    // --- What the reach exposed, and what the halo does about it ---------------------
    //
    // **A section is worth what it is worth wherever it was cut.** Cutting a window in
    // two used to conserve the plating exactly and the stiffener steel only where the
    // plating does not change thickness: `nodeThickness` was an area-weighted mean over
    // the sub-quads *inside the section*, so a station where the strake steps carried
    // one thickness to the section aft of it, another to the one forward of it, and the
    // mean of the two to a section that spanned it. A member run stops at a thickness
    // change (§3), so the spanning section stopped runs the two halves never saw and
    // dropped the seam node's run of one. **The spanning section is the one that
    // lost**, which is the direction the header states the mechanism in and the
    // opposite of the direction its summary states the figures in: one piece carried
    // 1.9% less than its two halves at x = -24 .. -19.2 -- inside the range that always
    // worked -- and **25.2%** less at the bow shoulder.
    //
    // `SectionParams::halo` averages both nodal fields over one bay beyond each plane
    // as well (§8), so the seam is in the same place whatever window contains it, and
    // the identity below is exact everywhere. The `halo = false` column is the negative
    // control: it is the old mesher, it still shows the whole defect, and it is what
    // keeps this from being a test of two numbers that were never going to differ.
    struct Split {
        const char* label;
        double from, middle, to;
        bool uniform;  // is the plating the same thickness across the shared station
    };
    for (const Split& split : {Split{"amidships", -7.2, -2.4, 2.4, true},
                               Split{"old reach", -24.0, -21.6, -19.2, false},
                               Split{"bow shoulder", 40.8, 43.2, 45.6, false}}) {
        const auto cut = [&](double from, double to, bool halo) {
            section::SectionParams p;
            p.xFrom = from;
            p.xTo = to;
            p.junctions = false;
            p.halo = halo;
            return section::buildSection(structure, p);
        };
        double loose = 0;
        for (int halo = 1; halo >= 0; --halo) {
            const section::Section aft = cut(split.from, split.middle, halo != 0);
            const section::Section forward = cut(split.middle, split.to, halo != 0);
            const section::Section whole = cut(split.from, split.to, halo != 0);
            // The plating, the elements and the wedges conserve whatever the thickness
            // does, with the halo or without it. Only the members were ever at issue.
            expectEqualCount(std::string(split.label) + ": two sections have one section's elements",
                             aft.elementCount() + forward.elementCount(), whole.elementCount());
            expectNear(std::string(split.label) + ": and its plating",
                       aft.plateMass + forward.plateMass, whole.plateMass, 1e-9);
            expectEqual(std::string(split.label) + ": and its collapsed elements",
                        aft.collapsedElements + forward.collapsedElements, whole.collapsedElements);

            const double halves = aft.memberMass + forward.memberMass;
            const double shortfall = 1.0 - whole.memberMass / halves;
            const int splitHalves =
                aft.memberRunsSplitByThickness + forward.memberRunsSplitByThickness;
            std::printf("     %-13s halo %d: stiffener steel halves %8.1f kg, one piece %8.1f"
                        " (%+.1f%%); runs stopped at a seam %d against %d, halo panels %d\n",
                        split.label, halo, halves, whole.memberMass, -100.0 * shortfall,
                        whole.memberRunsSplitByThickness, splitHalves, whole.haloPanels);
            if (halo == 0) {
                loose = shortfall;
                continue;
            }
            // **The property, asserted directly**: the same steel, and the same runs
            // stopped, however the length was cut. Exact, because the halo makes the
            // nodal thickness a function of the ship rather than of the window and
            // every run therefore breaks in the same places.
            expectNear(std::string(split.label) + ": the stiffener steel conserves exactly",
                       whole.memberMass, halves, 1e-9 * halves);
            expectEqual(std::string(split.label) + ": and no extra run is stopped",
                        whole.memberRunsSplitByThickness, splitHalves);
            // The halo has to have been made of something, or the two builds are the
            // same build and the control below proves nothing.
            expectTrue(std::string(split.label) + ": and the halo had plating to average",
                       whole.haloPanels > 0 && aft.haloPanels > 0);
        }
        if (split.uniform) {
            // The guard against the whole comparison being vacuous: amidships the old
            // mesher conserved too, so a station where it did not is what says the
            // identity above is about the fix and not about the ship.
            expectTrue(std::string(split.label) + ": and the old mesher conserved here as well",
                       std::abs(loose) < 1e-9);
        } else {
            expectTrue(std::string(split.label) + ": which the old mesher did not",
                       loose > 0.01);
        }
    }
}

// **A tie that would take more of a slave's mass than the slave has.** The junction
// search accepts a master whose *mid-surface* is within `junctionTolerance` -- 25 mm,
// an absolute figure -- and the through-thickness split is `(e + t/2) / t`, which is
// relative. On 10 mm plating those disagree: a slave 21.5 mm off the master's
// mid-surface is inside the tolerance and asks for a weight of 2.65, so `-1.65` of
// the slave's steel comes off one master face and `reduction::Substructure` refuses
// the whole section for a non-positive nodal mass.
//
// `SectionParams::junctionWeightLimit` refuses the **junction** instead, which is the
// choice `section.hpp` §1 makes everywhere else: left open and reported rather than
// closed wrongly and not.
void testJunctionWeightLimitRefusesTheNodeAndNotTheSection() {
    std::printf("\n--- section: a tie that would make the condensed mass negative ---\n");
    const StructuralMesh structure = ferryStructure();
    section::SectionParams params;
    params.xFrom = 45.6;
    params.xTo = 50.4;
    params.subdivision = 1;

    const section::Section limited = section::buildSection(structure, params);
    section::SectionParams unlimited = params;
    unlimited.junctionWeightLimit = 1e9;
    const section::Section loose = section::buildSection(structure, unlimited);

    const reduction::Substructure withLimit(limited.mesh, limited.material,
                                            limited.interfaceNodes, limited.attachment);
    const reduction::Substructure without(loose.mesh, loose.material, loose.interfaceNodes,
                                          loose.attachment);
    std::printf("     with the limit: %d ties, %d refused, worst weight %.4f, reduces %d\n",
                limited.junctionTies, limited.junctionsThroughThickness,
                limited.worstJunctionWeight, static_cast<int>(withLimit.ready()));
    std::printf("     without it:     %d ties, %d refused, worst weight %.4f, reduces %d\n",
                loose.junctionTies, loose.junctionsThroughThickness, loose.worstJunctionWeight,
                static_cast<int>(without.ready()));

    // **The negative control comes first**, because without it "the section reduces"
    // is what a mesher with no ties at all would report.
    expectEqual("without the limit nothing is refused", loose.junctionsThroughThickness, 0);
    expectTrue("and a tie takes more than a full share of the slave's mass off a face",
               loose.worstJunctionWeight > 2.0);
    expectTrue("so the substructure refuses the whole section", !without.ready());
    expectTrue("for the reason the limit is about",
               !without.problems().empty() &&
                   without.problems().front().find("lumped mass") != std::string::npos);

    // And with it, the two junctions are refused and the section is usable.
    expectEqual("the limit refuses exactly the two junctions that do it",
                limited.junctionsThroughThickness, 2);
    expectTrue("and says so", [&] {
        for (const std::string& problem : limited.problems)
            if (problem.find("junctionWeightLimit") != std::string::npos) return true;
        return false;
    }());
    expectTrue("every tie that is left is inside the limit", limited.worstJunctionWeight <= 2.0);
    expectTrue("and the section reduces", withLimit.ready());
    expectEqual("in one piece", limited.components, 1);
    // The ties it did keep are still doing something -- the section is not merely
    // usable because it stopped tying anything.
    expectTrue("having kept most of them", limited.junctionTies > 0 &&
                                               limited.junctionTies >= loose.junctionTies - 2);

    // Amidships nothing is refused, so the limit is not quietly costing the middle
    // body its junctions. The ferry's ordinary deck-to-shell junction runs at 1.69,
    // which is a *negative* 0.69 on the inner face and is where the steel is.
    section::SectionParams midship = params;
    midship.xFrom = -2.4;
    midship.xTo = 2.4;
    const section::Section middle = section::buildSection(structure, midship);
    expectEqual("amidships the limit refuses nothing", middle.junctionsThroughThickness, 0);
    expectTrue("though its junctions do extrapolate through the master's thickness",
               middle.worstJunctionWeight > 1.5 && middle.worstJunctionWeight <= 2.0);
}

// **The interface is coincident because the halo makes it so, and the parallel middle
// body is what used to hide that it was not.** See `section.hpp` §6 note 1 and §8. A
// node is the mid-surface offset by `t/2` along the nodal normal; averaged over the
// sub-quads *inside* the section, both that normal and that thickness are statements
// about the window, so two sections cut on the same station disagreed about where the
// node was wherever the hull is not prismatic -- 3.4 µm at x = -21.6, three orders of
// magnitude above `matchBoundaries`' default, and 336 of the plane's 1 170 boundary
// DOF finding no partner at all.
//
// Averaged over one bay beyond the plane as well they are statements about the ship,
// and the two sections agree **bit for bit**: not to a tolerance, because they sum the
// same contributions from the same panels in the same order. `halo = false` below is
// the old mesher and is what keeps this from being a comparison of two zeros.
void testHaloMakesTheCutPlanesCoincide() {
    std::printf("\n--- section: two sections agreeing about the cut plane ---\n");
    const StructuralMesh structure = ferryStructure();
    section::ChainParams params;
    params.section.subdivision = 1;
    params.section.members = false;
    params.section.junctions = false;
    params.reduce.modes = 0;
    params.reduce.cutoffFrequency = 0;
    params.station = {-26.4, -21.6, -16.8};

    const section::Chain joined = section::buildChain(structure, params);
    expectEqualCount("the chain meshed two sections", joined.section.size(), 2u);
    expectEqualCount("one interior plane", joined.unmatched.size(), 1u);
    const std::size_t planeDof = 3 * joined.section[0].forwardNodes.size();
    expectEqualCount("of 1170 boundary DOF", planeDof, 1170u);
    std::printf("     halo on:  %zu of %zu unmatched at 1e-9, worst matched gap %.4e m\n",
                joined.unmatched[0], planeDof, joined.worstGap);
    // **The sharpest check available, and it needs no new instrument.**
    expectEqualCount("at matchBoundaries' own default nothing is left over", joined.unmatched[0],
                     0u);
    expectEqualCount("and the whole plane is shared", joined.shared[0], planeDof);
    expectTrue("so the chain is ready without a tolerance being chosen for it", joined.ready());
    // Exactly zero, not nearly: the two sections form the same sums from the same
    // panels in the same order, so the node is at the same double.
    expectNear("the two sections put the node at the same double", joined.worstGap, 0.0, 0.0);
    // The assembly's own account of it, which is a different measurement of the same
    // thing: `Chain::worstGap` comes from the interface map, `worstMergedGap` from
    // comparing the identities of the DOF that actually landed on one assembled row.
    expectNear("and the assembly reports the same gap on the rows it merged",
               joined.assembly.worstMergedGap, joined.worstGap, 1e-15);
    expectEqual("with no axis crossed", joined.assembly.axisDisagreements, 0);
    expectTrue("and the halo had plating beyond the cut to average",
               joined.section[0].haloPanels > 0 && joined.section[1].haloPanels > 0);

    // --- The negative control: the same chain built the way it used to be -----------
    section::ChainParams unaveraged = params;
    unaveraged.section.halo = false;
    const section::Chain strict = section::buildChain(structure, unaveraged);
    std::printf("     halo off: %zu of %zu unmatched at 1e-9\n", strict.unmatched[0], planeDof);
    expectEqualCount("without the halo 336 of the plane's 1170 DOF find no partner",
                     strict.unmatched[0], 336u);
    expectTrue("and the chain is refused rather than solved", !strict.ready());
    section::ChainParams loose = unaveraged;
    loose.matchTolerance = 1e-5;
    const section::Chain tolerated = section::buildChain(structure, loose);
    expectTrue("a tolerance the size of the disagreement joins it instead", tolerated.ready());
    // Measured at 3.4021e-06 m. Bracketed rather than asserted at a point, because it
    // is a property of this hull's shoulder -- but bracketed tightly enough that a
    // mesher which had started averaging its normals across the cut would fail it,
    // which is exactly what the halo does and what the assertion above says it did.
    expectTrue("and the disagreement it is covering is microns, not the plate thickness",
               tolerated.worstGap > 1e-6 && tolerated.worstGap < 1e-5);

    // --- And what the halo buys the in-plane line ties of §9 -------------------------
    //
    // The two sections either side of an interior cut plane derive that plane's line
    // ties **independently**, and the two derivations have to come out identical or the
    // ship is tied to itself twice in two different ways. What makes them identical is
    // the halo: the plane's nodal normals and nodal thicknesses are then properties of
    // the ship rather than of the window, so both sides start from the same doubles and
    // run the same arithmetic over them.
    //
    // The masters are the same either way -- they are chosen on the mid-surface, which
    // comes from panel corners and never depended on the cut. What moves is the
    // **weights**, so this is a difference no structural comparison would find, and it
    // is why `worstPlaneTieDisagreement` is compared against zero rather than against a
    // tolerance.
    section::ChainParams withTies = params;
    withTies.section.junctions = true;
    const section::Chain tiedJoined = section::buildChain(structure, withTies);
    expectTrue("the tied chain built", tiedJoined.ready());
    expectTrue("and tied the shoulder's interior plane", tiedJoined.planeTieNodes > 0);
    expectEqual("with both sides agreeing about all of it", tiedJoined.planeTiesDisagreeing, 0);
    expectNear("to the last bit", tiedJoined.worstPlaneTieDisagreement, 0.0, 0.0);

    section::ChainParams tiedNoHalo = withTies;
    tiedNoHalo.section.halo = false;
    tiedNoHalo.matchTolerance = 1e-5;  // or the plane does not match at all, as above
    const section::Chain tiedCut = section::buildChain(structure, tiedNoHalo);
    std::printf("     line ties on the shoulder plane: %d with the halo (%d planes disagreed),"
                " %d without (%d disagreed, worst %.3e)\n",
                tiedJoined.planeTieNodes, tiedJoined.planeTiesDisagreeing, tiedCut.planeTieNodes,
                tiedCut.planeTiesDisagreeing, tiedCut.worstPlaneTieDisagreement);
    expectTrue("without the halo the two sides derive a different constraint",
               tiedCut.planeTiesDisagreeing > 0);
    expectTrue("and it is a difference in the weights, not a tie one side simply lacks",
               tiedCut.worstPlaneTieDisagreement > 0 && tiedCut.worstPlaneTieDisagreement < 1.0);
    expectEqual("so the chain applies nothing on that plane rather than picking a side",
                tiedCut.planeTieNodes, 0);
    expectTrue("and says which plane", !tiedCut.problems.empty());

    // The guard against the whole test being about nothing: amidships the same
    // measurement is exactly zero, so this is the hull's shape and not the mesher's
    // arithmetic.
    section::ChainParams amidships = params;
    amidships.station = {-7.2, -2.4, 2.4};
    const section::Chain prismatic = section::buildChain(structure, amidships);
    expectTrue("amidships the same chain needs no tolerance at all", prismatic.ready());
    expectNear("and the planes coincide exactly", prismatic.worstGap, 0.0, 0.0);
    expectNear("by both measurements", prismatic.assembly.worstMergedGap, 0.0, 0.0);

    // --- And the one solved comparison at ship scale ------------------------------
    //
    // The box says the assembly is exact; this says it on real geometry, with 1 170
    // boundary DOF a plane and a section that comes apart into seven pieces. Torsion
    // and not `EA`, for the reason this file exists to make: `EA` is exact on a chain
    // that joined nothing. The rest of the ship-scale matrix -- `EA`, the neutral
    // axis, `EI`, and what all of it costs tied -- is in `tools/section_probe
    // --chain=N`, because a dense Cholesky at (N+1) x 1 170 is core-seconds and this
    // one already is.
    const section::Section mono = section::buildSection(structure, [&] {
        section::SectionParams p = amidships.section;
        p.xFrom = -7.2;
        p.xTo = 2.4;
        return p;
    }());
    expectEqual("the chain and the monolith are the same seven pieces", prismatic.components,
                mono.components);
    const section::TorsionResponse chainTwist = section::applyTwist(prismatic, 1e-6, 6.86);
    const section::TorsionResponse monoTwist =
        section::applyTwist(mono, mono.material, 1e-6, 6.86);
    expectTrue("both twists ran", chainTwist.ok && monoTwist.ok);
    std::printf("     ship scale: GJ chain %.6e against one piece %.6e (%+.2e)\n",
                chainTwist.torsionalStiffness, monoTwist.torsionalStiffness,
                chainTwist.torsionalStiffness / monoTwist.torsionalStiffness - 1);
    // Measured at 1.4e-11 relative, with the halo and without it -- the comment this
    // replaces said 3.9e-12 and no run of this file has produced that figure, which is
    // what a tolerance three decades wide lets happen to a number nothing asserts.
    // Asserted at 1e-9, two decades above the measurement and six below anything a
    // lost interface could hide in.
    expectNear("two ferry sections carry the torque one ferry section carries",
               chainTwist.torsionalStiffness, monoTwist.torsionalStiffness,
               1e-9 * monoTwist.torsionalStiffness);
    // Vacuity: an open cell and a closed one differ by orders of magnitude, so this
    // is a comparison of something rather than of two zeros.
    expectTrue("and it is a real torque", monoTwist.torsionalStiffness > 1e11);

    // --- And the same, on the shoulder, which is what the halo buys ----------------
    //
    // The comparison above is the middle body, where the interface always matched.
    // This is the one that could not be made at the default tolerance at all: the
    // chain across x = -21.6 against the same length in one piece. Torsion again,
    // because `EA` would be exact on a chain torn along 29% of its cut.
    const section::Section shoulder = section::buildSection(structure, [&] {
        section::SectionParams p = params.section;
        p.xFrom = -26.4;
        p.xTo = -16.8;
        return p;
    }());
    expectEqual("the shoulder chain and its monolith are the same pieces", joined.components,
                shoulder.components);
    const section::TorsionResponse shoulderChain = section::applyTwist(joined, 1e-6, 6.86);
    const section::TorsionResponse shoulderMono =
        section::applyTwist(shoulder, shoulder.material, 1e-6, 6.86);
    expectTrue("both shoulder twists ran", shoulderChain.ok && shoulderMono.ok);
    std::printf("     shoulder:   GJ chain %.6e against one piece %.6e (%+.2e)\n",
                shoulderChain.torsionalStiffness, shoulderMono.torsionalStiffness,
                shoulderChain.torsionalStiffness / shoulderMono.torsionalStiffness - 1);
    // Measured at 1.4e-9 relative, against 1.4e-11 amidships, and the hundredfold is
    // the solve rather than the assembly: the same two loads report a residual of
    // 3.7e-5 N here against 1.0e-7 N amidships, because a shoulder section is a worse
    // conditioned mesh than a prismatic one. Asserted at 1e-8 -- seven times the
    // measurement and five decades below the 29% of a cut that used to be torn.
    expectNear("two shoulder sections carry the torque one shoulder section carries",
               shoulderChain.torsionalStiffness, shoulderMono.torsionalStiffness,
               1e-8 * shoulderMono.torsionalStiffness);
    expectTrue("and it too is a real torque", shoulderMono.torsionalStiffness > 1e11);
}

// **The halo reaches exactly as far as the weld does, and nothing on the ferry could
// say so.** Whether a panel joins the halo is decided by welding its corners against
// the section's own, and the reference ship's panel corners are *bit* identical where
// they meet -- they come from the same evaluation of the same station -- so the
// distance in that test is never exercised and any tolerance at all would pass. That
// is the shape of the mutant that survived everything in `zone.cpp`: a probe that
// reaches less far than the tolerance leaves a crack down the middle of any mesh whose
// duplicates sit between 0.75 and 1 tolerance apart.
//
// So: a strip of plating that **kinks** at its middle station -- so a node there has
// two different normals to average and its position depends on whether the halo
// reached -- with the forward bay's aft corners nudged sideways by a fraction of
// `weldTolerance`. The weld joins them; the halo has to agree, and at 0.75 it is the
// only thing being asked.
StructuralMesh makeKinkedStrip(double nudge, double slope = 0.5) {
    StructuralMesh mesh;
    // A second material for the forward bay, differing in the two things the section
    // census reads. It is the halo's material and never the section's, which is the
    // only way to ask whether a halo panel can make a section report spanning a
    // material it does not contain.
    StructuralMaterial soft = ah36Steel();
    soft.name = "halo steel";
    soft.youngsModulus *= 0.5;
    soft.density *= 0.5;
    mesh.materials = {ah36Steel(), soft};
    mesh.frameSpacing = 1.0;
    // Flat aft of x = 1 and rising forward of it. At `slope = 0.5` that is 26.6
    // degrees -- inside `foldLimit`, so the two bays are one surface and their nodes
    // weld. A node at the kink is then 0.44 rad of normal away from either answer,
    // which is what makes its position depend on whether the halo reached it. At
    // `slope = 4` it is 76 degrees and the two bays are two surfaces.
    const auto height = [&](double x) { return x <= 1.0 ? 0.0 : (x - 1.0) * slope; };
    for (int i = 0; i < 2; ++i) {
        // Both panels of the forward bay carry the same nudge, so they still share
        // their own corners exactly and only the *seam between the bays* is opened.
        const double x0 = i, x1 = i + 1, shift = i == 1 ? nudge : 0.0;
        for (int j = 0; j < 3; ++j) {
            const double y0 = 0.5 * j, y1 = 0.5 * (j + 1);
            PlatePanel p;
            p.corner[0] = {x0, y0 + shift, height(x0)};
            p.corner[1] = {x1, y0, height(x1)};
            p.corner[2] = {x1, y1, height(x1)};
            p.corner[3] = {x0, y1 + shift, height(x0)};
            p.thickness = 0.010;
            p.material = i;
            p.role = PanelRole::Shell;
            p.source = 0;
            mesh.panels.push_back(p);
        }
    }
    mesh.frameStations = {0.0, 1.0, 2.0};
    return mesh;
}

void testHaloReachesAsFarAsTheWeld() {
    std::printf("\n--- section: the halo reaches as far as the weld, not as far as a cell ---\n");
    const auto planeGap = [](double nudge, int& haloPanels) {
        const StructuralMesh strip = makeKinkedStrip(nudge);
        const auto cut = [&](double from, double to) {
            section::SectionParams p;
            p.xFrom = from;
            p.xTo = to;
            p.subdivision = 1;
            p.members = false;
            p.junctions = false;
            return section::buildSection(strip, p);
        };
        const section::Section aft = cut(0.0, 1.0), forward = cut(1.0, 2.0);
        haloPanels = aft.haloPanels;
        double worst = 0;
        for (std::uint32_t a : aft.forwardNodes) {
            const Vec3 at{aft.mesh.position[a * 3], aft.mesh.position[a * 3 + 1],
                          aft.mesh.position[a * 3 + 2]};
            double nearest = 1e300;
            for (std::uint32_t b : forward.aftNodes)
                nearest = std::min(
                    nearest, length(at - Vec3{forward.mesh.position[b * 3],
                                              forward.mesh.position[b * 3 + 1],
                                              forward.mesh.position[b * 3 + 2]}));
            worst = std::max(worst, nearest);
        }
        return worst;
    };

    // The default `weldTolerance` is 1e-6 m.
    int panels = 0;
    const double joined = planeGap(0.0, panels);
    std::printf("     corners coincident:      halo %d panels, plane gap %.3e m\n", panels, joined);
    expectTrue("coincident corners put the halo in reach", panels > 0);
    expectNear("and the kink's nodes land at the same double from either side", joined, 0.0, 0.0);

    // **0.75 of the tolerance**, which is the gap a probe reaching one cell instead of
    // one tolerance would miss. Nothing else in this file puts a corner there.
    const double nudged = planeGap(0.75e-6, panels);
    std::printf("     corners 0.75e-6 m apart: halo %d panels, plane gap %.3e m\n", panels, nudged);
    expectTrue("a corner three quarters of a tolerance away is still in reach", panels > 0);
    expectNear("and the kink's nodes still land at the same double", nudged, 0.0, 0.0);

    // The control, and it is what says the two above measure a distance rather than
    // passing by construction: past the tolerance the corners are genuinely two
    // corners, the weld does not join them either, and the plane comes apart.
    const double apart = planeGap(4e-6, panels);
    std::printf("     corners 4e-6 m apart:    halo %d panels, plane gap %.3e m\n", panels, apart);
    expectEqual("past the tolerance the panel is not in the halo at all", panels, 0);
    expectTrue("and the two sections stop agreeing about the plane", apart > 1e-4);
    // A millimetre and not a micron: the kink is what makes the disagreement a
    // *normal* rather than a corner offset, so this is the defect the halo removes
    // rather than the nudge itself showing through.
    expectTrue("by the plating's own half-thickness, not by the nudge", apart > 1e-3);

    // --- What the halo is looked at for, and nothing else ---------------------------
    //
    // The forward bay carries a second material, half the modulus and half the
    // density. It is the halo's and never the section's, so the aft section must
    // report one material and mesh at the first -- a section that took its material
    // census over the halo would report spanning two and hand `reduction::Substructure`
    // a stiffness that is not its own. Nothing on the ferry could say so: its two
    // materials differ in yield alone, which the census does not read.
    const StructuralMesh twoMaterials = makeKinkedStrip(0.0);
    section::SectionParams aftHalf;
    aftHalf.xFrom = 0.0;
    aftHalf.xTo = 1.0;
    aftHalf.subdivision = 1;
    aftHalf.members = false;
    aftHalf.junctions = false;
    const section::Section aft = section::buildSection(twoMaterials, aftHalf);
    expectTrue("the halo is there to be miscounted", aft.haloPanels > 0);
    expectNear("the section is meshed at its own plating's modulus",
               aft.material.youngsModulus, ah36Steel().youngsModulus, 0.0);
    bool spanned = false;
    for (const std::string& problem : aft.problems)
        spanned = spanned || problem.find("spans materials") != std::string::npos;
    expectTrue("and does not report spanning the halo's", !spanned);
    // The control: the section that really does contain both reports it, so the
    // assertion above is about the halo and not about a check that never fires.
    section::SectionParams both = aftHalf;
    both.xTo = 2.0;
    const section::Section whole = section::buildSection(twoMaterials, both);
    bool reported = false;
    for (const std::string& problem : whole.problems)
        reported = reported || problem.find("spans materials") != std::string::npos;
    expectTrue("a section that does contain both says so", reported);

    // --- And a fold the halo is on the far side of ----------------------------------
    //
    // At 76 degrees the two bays are two surfaces, so the halo cannot weld to the
    // section's nodes and cannot -- and must not -- make them agree: a shared node
    // pair has one thickness direction and two plates at that angle need two. What it
    // must still do is *say* so. The fold count is a fact about the section's own
    // boundary, so an edge with the section on one side of it and the halo on the
    // other belongs in it; on the ferry no fold on a cut plane exists to check that.
    const StructuralMesh folded = makeKinkedStrip(0.0, 4.0);
    const section::Section overFold = section::buildSection(folded, aftHalf);
    bool saidFold = false;
    for (const std::string& problem : overFold.problems)
        saidFold = saidFold || problem.find("fold further than") != std::string::npos;
    std::printf("     across a 76 degree fold: halo %d panels, %d surfaces, fold reported %d\n",
                overFold.haloPanels, overFold.surfaces, static_cast<int>(saidFold));
    expectTrue("a fold between the section and its halo is the section's fold", saidFold);
    // And the halo on the far side of it reaches nothing, which is what a weld class
    // keyed on the surface is for.
    expectEqual("while the plating past it joins no node of this section",
                overFold.haloPanels, 0);
    expectEqual("so the section is the one surface it is made of", overFold.surfaces, 1);
}

// **The property the halo exists for, asserted directly**: a section's nodes are where
// they are because of the ship, not because of where it was cut. `section.hpp` §8.
//
// The two symptoms this repository had recorded -- an interface 3.4 µm out at
// x = -21.6, and 25.2% of the stiffener steel appearing or disappearing at the bow
// shoulder depending on the window -- are both consequences of its absence, so this
// asserts the cause and they follow. It is asserted at **zero**: the same node is
// reached by the same panels whichever window contains it, the halo makes `candidate`
// hold all of them in ascending panel order whichever window it is, and the nodal sums
// are therefore formed from the same terms in the same order. Anything less than a bit
// would mean one of those three had stopped being true.
void testNodePositionsDoNotDependOnWhereTheSectionWasCut() {
    std::printf("\n--- section: a node is where the ship puts it, not where the cut does ---\n");
    const StructuralMesh structure = ferryStructure();

    // Every mesh node position of a section, sorted, so two sections built with
    // different numberings can be compared term by term. Positions and not indices:
    // the numbering is chosen for bandwidth and owes nothing to the ship.
    const auto positions = [](const section::Section& s, const std::vector<std::uint32_t>& nodes) {
        std::vector<std::array<double, 3>> out;
        out.reserve(nodes.size());
        for (std::uint32_t n : nodes)
            out.push_back({s.mesh.position[static_cast<std::size_t>(n) * 3],
                           s.mesh.position[static_cast<std::size_t>(n) * 3 + 1],
                           s.mesh.position[static_cast<std::size_t>(n) * 3 + 2]});
        std::sort(out.begin(), out.end());
        return out;
    };

    struct Where {
        const char* label;
        double station;   // the plane three windows are cut on
        double bay;       // how far either side each window reaches
        int subdivision;
    };
    // Two stations, and they fail in different ways without the halo: the hull turns
    // at x = -21.6 and the nodal *normal* is what disagrees, while the plating steps
    // at x = 43.2 and the nodal *thickness* is. A halo that averaged only one of the
    // two would pass one of these and fail the other.
    //
    // The third is the second station refined, and it asks a question one element per
    // panel cannot. The halo is "every panel sharing a welded corner with one
    // inside", which is exactly the set that can reach a node **when every node is a
    // panel corner**. At `subdivision = 2` a node can also sit part way along a panel
    // edge, where the panels that reach it are the ones sharing that edge -- two
    // corners rather than one, so still inside the rule, but nothing tested it.
    for (const Where& where : {Where{"stern shoulder, the normal turns", -21.6, 4.8, 1},
                               Where{"bow shoulder, the plating steps", 43.2, 2.4, 1},
                               Where{"bow shoulder, refined", 43.2, 2.4, 2}}) {
        double loose = 0;
        for (int halo = 1; halo >= 0; --halo) {
            const auto cut = [&](double from, double to) {
                section::SectionParams p;
                p.xFrom = from;
                p.xTo = to;
                p.subdivision = where.subdivision;
                p.junctions = false;
                p.halo = halo != 0;
                return section::buildSection(structure, p);
            };
            const section::Section aft = cut(where.station - where.bay, where.station);
            const section::Section forward = cut(where.station, where.station + where.bay);
            const section::Section whole =
                cut(where.station - where.bay, where.station + where.bay);
            expectTrue(std::string(where.label) + ": all three meshed",
                       !aft.empty() && !forward.empty() && !whole.empty());

            // The aft section's forward plane, the forward section's aft plane, and
            // the same station seen from inside the section that spans it. Three
            // different windows, one physical ring of nodes.
            const std::vector<std::array<double, 3>> fromAft = positions(aft, aft.forwardNodes);
            const std::vector<std::array<double, 3>> fromForward =
                positions(forward, forward.aftNodes);
            expectEqualCount(std::string(where.label) + ": the plane has the same nodes both ways",
                             fromAft.size(), fromForward.size());
            // **Nearest neighbour, not term by term.** Sorting by coordinate is only a
            // usable pairing once the two agree; a micron of disagreement in x
            // reorders the sort and the term-by-term difference becomes the width of
            // the ship rather than the error. This is the same question
            // `matchBoundaries` asks, asked without a tolerance.
            double worst = 0;
            for (const std::array<double, 3>& a : fromAft) {
                double nearest = 1e300;
                for (const std::array<double, 3>& b : fromForward) {
                    const double dx = a[0] - b[0], dy = a[1] - b[1], dz = a[2] - b[2];
                    nearest = std::min(nearest, std::sqrt(dx * dx + dy * dy + dz * dz));
                }
                worst = std::max(worst, nearest);
            }
            std::printf("     %-34s halo %d: %3zu nodes on the plane, furthest from its"
                        " partner %.3e m, halo panels %d\n",
                        where.label, halo, fromAft.size(), worst, whole.haloPanels);
            if (halo == 0) {
                loose = worst;
                continue;
            }
            expectNear(std::string(where.label) + ": and puts them at the same doubles", worst,
                       0.0, 0.0);
            expectTrue(std::string(where.label) + ": on a plane that carries a ship's worth"
                                                  " of nodes",
                       fromAft.size() > 100);
            expectTrue(std::string(where.label) + ": with a halo to average", whole.haloPanels > 0);

            // And the stronger statement, which the plane alone does not make: the
            // section that *spans* the station agrees with both halves about every
            // node either of them has, not only the ones on the shared cut.
            std::vector<std::uint32_t> all(whole.mesh.nodeCount());
            for (std::uint32_t i = 0; i < whole.mesh.nodeCount(); ++i) all[i] = i;
            std::vector<std::array<double, 3>> spanning = positions(whole, all);
            std::vector<std::uint32_t> everyAft(aft.mesh.nodeCount()),
                everyForward(forward.mesh.nodeCount());
            for (std::uint32_t i = 0; i < aft.mesh.nodeCount(); ++i) everyAft[i] = i;
            for (std::uint32_t i = 0; i < forward.mesh.nodeCount(); ++i) everyForward[i] = i;
            std::vector<std::array<double, 3>> halves = positions(aft, everyAft);
            const std::vector<std::array<double, 3>> other = positions(forward, everyForward);
            halves.insert(halves.end(), other.begin(), other.end());
            std::sort(halves.begin(), halves.end());
            // The shared plane is in both halves and once in the monolith, so the
            // halves carry it twice. Deduplicating by value would hide a genuine
            // duplicate; removing exactly the plane's nodes once does not.
            bool inBoth = true;
            for (const std::array<double, 3>& p : fromForward) {
                const auto found = std::lower_bound(halves.begin(), halves.end(), p);
                if (found == halves.end() || *found != p) {
                    inBoth = false;
                    break;
                }
                halves.erase(found);
            }
            expectTrue(std::string(where.label) + ": the shared plane is in both halves", inBoth);
            expectEqualCount(std::string(where.label) + ": two halves have one section's nodes",
                             halves.size(), spanning.size());
            bool identical = halves.size() == spanning.size();
            for (std::size_t i = 0; i < halves.size() && i < spanning.size(); ++i)
                identical = identical && halves[i] == spanning[i];
            expectTrue(std::string(where.label) + ": every one of them at the same double",
                       identical);
        }
        // The guard, and it is the load-bearing one: without the halo this station
        // disagrees by microns to millimetres, so "they agree" is a statement about
        // the fix rather than about a hull that was prismatic all along.
        std::printf("     %-34s without the halo the same plane is %.3e m out\n", where.label,
                    loose);
        expectTrue(std::string(where.label) + ": and the old mesher did not agree",
                   loose > 1e-9);
    }

    // Amidships, where it was always true. The two builds have to be the *same* there
    // or the halo has moved a mesh that had nothing wrong with it -- which is the one
    // way this change could regress the part of the ship every published figure in
    // `section.hpp` was measured on.
    const auto midship = [&](bool halo) {
        section::SectionParams p;
        p.xFrom = -7.2;
        p.xTo = 2.4;
        p.subdivision = 1;
        p.junctions = false;
        p.halo = halo;
        return section::buildSection(structure, p);
    };
    const section::Section withHalo = midship(true), without = midship(false);
    expectEqualCount("amidships the halo changes no node count", withHalo.mesh.nodeCount(),
                     without.mesh.nodeCount());
    double moved = 0;
    for (std::size_t d = 0; d < withHalo.mesh.position.size() &&
                            d < without.mesh.position.size();
         ++d)
        moved = std::max(moved, std::abs(withHalo.mesh.position[d] - without.mesh.position[d]));
    std::printf("     amidships the halo moves the mesh by %.3e m\n", moved);
    // **Rounding, and the bound is derived rather than chosen.** The middle body is
    // prismatic, so the halo's sub-quads carry the same normal and the same thickness
    // as the ones inside: the average is the same vector formed from twice as many
    // terms, which differs from it by a few ulps before `normalize`. A node is that
    // unit vector times `t/2`, so the move is bounded by a handful of ulps of half the
    // plating -- 8 eps * 0.008 m = 1.4e-17. Measured at 1.7e-18, which is one.
    //
    // It is fourteen orders of magnitude below the millimetre the halo removes at the
    // ends and eight below anything `matchBoundaries` would notice, so "the middle
    // body did not move" is a statement and not a tolerance.
    expectTrue("and moves no node amidships beyond the rounding of a longer sum",
               moved < 2e-17);
    expectEqual("nor changes what is collapsed there", withHalo.collapsedElements, 0);
    expectTrue("though it did average a halo there too", withHalo.haloPanels > 0);
}

// --- 10. The stress the section carries, against the only field it can be -----------
//
// `section.hpp` §2 is that `EA` and `EI` are the wrong instrument for the junctions
// because a prescribed plane-sections field makes them come out right on a mesh that
// joined nothing. **The stress path has the mirror-image hazard**: on a *prismatic*
// section the exact three-dimensional answer to that same prescribed field **is**
// plane sections, so a stress reader that returned `E eps` from a lookup table would
// pass everything below. What stops that being vacuous is that the same reader is
// asked for a field it must get wrong -- the ferry, where the residual is not zero --
// and that the box is asked at Gauss points whose `z` the fit has to discover.

void testAxialStressIsThePatchTest() {
    const StructuralMesh structure = makeBox(kT, kT);
    const section::Section box = section::buildSection(structure, openBoxParams());
    const StructuralMaterial material = ah36Steel();
    expectTrue("the box meshed", !box.empty());

    section::BeamLoad axial;
    axial.strain = 1e-4;
    std::vector<double> field;
    const section::BeamResponse stretched =
        section::applyBeamLoad(box, material, axial, &field);
    expectTrue("the box takes an axial load: " + stretched.problem, stretched.ok);
    expectEqualCount("and hands back a field over its own mesh", field.size(),
                     box.nodeCount() * 3);

    const std::vector<section::StressSample> samples =
        section::axialStress(box, material, field, 0.5 * kL);
    // Twelve sub-quads cross a bay of this box -- four in each flange, two in each
    // side -- so ninety-six Gauss points, and the bound is set off that rather than off
    // a round number.
    expectTrue("a transverse plane crosses elements", samples.size() >= 96);

    // The patch test. A prismatic section stretched by `eps` carries `sigma = E eps`
    // at every point, because every longitudinal strip is free to contract and the
    // ends alone say what its strain is. Nothing here is a discretisation: it is what
    // the element is supposed to reproduce exactly.
    const double exact = material.youngsModulus * axial.strain;
    const section::BeamFit fit = section::fitBeam(samples, 0.5 * kH);
    const double fitVolume = fit.volume;
    double worst = 0, volume = 0;
    for (const section::StressSample& s : samples) {
        worst = std::max(worst, std::abs(s.sigmaXX / exact - 1.0));
        volume += s.volume;
        expectTrue("every Gauss point carries a volume", s.volume > 0);
    }
    std::printf("     box under pure strain: worst |sigma/E eps - 1| over %zu Gauss points"
                " = %.3e\n", samples.size(), worst);
    expectTrue("sigma_xx is E eps at every Gauss point on the cut, to 5e-11 (measured"
               " 1.16e-11)", worst < 5e-11);
    // And it is not vacuous: `exact` is 20.6 MPa, so a reader returning zero would
    // fail.
    expectTrue("against a stress worth reading", exact > 1e7);
    // **The volume is a closed form and it is asserted as one**, because it is the
    // weight in every average `fitBeam` and `fibreStress` take and a mutation run said
    // nothing here could tell it from a constant. One bay of this box is its perimeter
    // times the plate thickness times the bay length -- the Gauss weights of the
    // elements the plane crosses sum to exactly the steel in that slab, which is what
    // makes a volume-weighted mean a mean over steel rather than over mesh.
    expectNear("the Gauss weights sum to the slab's own steel", volume,
               2 * (kB + kH) * kT * (kL / kNx), 1e-12 * (2 * (kB + kH) * kT * (kL / kNx)));
    expectNear("which the fit agrees about", fitVolume, volume, 1e-15 * volume);

    expectTrue("the fit ran", fit.ok);
    expectNear("a uniform field fits a constant", fit.axial, exact, 1e-11 * exact);
    expectTrue("with no gradient in z", std::abs(fit.gradient) < 1e-6 * exact / kH);
    expectTrue("and nothing left over", fit.residualRms < 1e-11 * exact);
    // The guard that says the fit had a `z` range to find a gradient in. Without it a
    // field sampled at one height would report "no gradient" and mean nothing.
    expectNear("the samples reach the keel", fit.zLo, 0.0, 0.1 * kH);
    expectNear("and the deck", fit.zHi, kH, 0.1 * kH);

    // **And the same field fitted about an axis that is not the section's own.** Every
    // fit above turns about the neutral axis, where the axial term is zero and the
    // cross term `s1 t0` in the normal equations vanishes with it -- so a sign error
    // there is invisible, which a mutation run said in as many words. About the keel
    // both are non-zero and the gradient still has to come out at nothing: a uniform
    // field has no gradient whatever it is measured from.
    const section::BeamFit low = section::fitBeam(samples, 0.0);
    expectTrue("the fit about the keel ran", low.ok);
    expectNear("a uniform field is uniform from any axis", low.axial, exact, 1e-11 * exact);
    expectTrue("and still has no gradient in z, to 1e-9 of the extreme fibre (measured 2e-11)",
               std::abs(low.gradient) * kH < 1e-9 * exact);
    // Vacuity: the cross term has to be non-zero, or "the sign of it does not matter"
    // is what is being asserted rather than tested. `sum w z` over a box about its keel
    // is half its depth times its steel.
    expectTrue("about the keel the first moment of the samples is not zero",
               std::abs(fit.volume * 0.5 * kH) > 1e-4);

    // **And the same field in compression, because a ship sags as well as hogs.**
    // `fitBeam` is arithmetic over a sample list, so the sign case is free: negate every
    // stress and every reported quantity has to follow. A `peak` taken as the largest
    // *signed* value instead of the largest magnitude reports zero here, and nothing
    // that only ever looks at a hogging deck would notice.
    std::vector<section::StressSample> compressed = samples;
    for (section::StressSample& c : compressed) c.sigmaXX = -c.sigmaXX;
    const section::BeamFit sagging = section::fitBeam(compressed, 0.5 * kH);
    expectTrue("the mirrored fit ran", sagging.ok);
    expectNear("a uniform compression fits its own constant", sagging.axial, -exact,
               1e-11 * exact);
    expectNear("and the peak carries the sign", sagging.peak, -exact, 1e-11 * exact);
    expectNear("as does the unmirrored one", fit.peak, exact, 1e-11 * exact);
}

void testAxialStressUnderCurvatureIsTheBeamsOwnLine() {
    const StructuralMesh structure = makeBox(kT, kT);
    const section::Section box = section::buildSection(structure, openBoxParams());
    const StructuralMaterial material = ah36Steel();

    // About the box's own neutral axis, which is mid-height by symmetry. Turning about
    // anything else would superpose an axial strain and the two would have to be
    // separated before anything could be asserted.
    section::BeamLoad bending;
    bending.curvature = 1e-4;
    bending.reference = 0.5 * kH;
    std::vector<double> field;
    const section::BeamResponse bent = section::applyBeamLoad(box, material, bending, &field);
    expectTrue("the box takes a curvature: " + bent.problem, bent.ok);

    const std::vector<section::StressSample> samples =
        section::axialStress(box, material, field, 0.5 * kL);
    const double gradient = material.youngsModulus * bending.curvature;
    double worst = 0;
    for (const section::StressSample& s : samples)
        worst = std::max(worst, std::abs(s.sigmaXX - gradient * (s.at.z - 0.5 * kH)));
    std::printf("     box under pure curvature: worst |sigma - E kappa (z - z_na)| = %.3e Pa"
                " against a %.3e Pa extreme fibre\n", worst, gradient * 0.5 * kH);
    expectTrue("sigma_xx is E kappa (z - z_na) at every Gauss point, to 100 Pa of a 10.3 MPa"
               " extreme fibre (measured 28 Pa)", worst < 100.0);

    const section::BeamFit fit = section::fitBeam(samples, 0.5 * kH);
    expectTrue("the fit ran", fit.ok);
    std::printf("     fitted gradient %.6e against E kappa %.6e (%.3e relative), fitted"
                " neutral axis %.9f\n", fit.gradient, gradient, fit.gradient / gradient - 1.0,
                fit.neutralAxis);
    expectNear("and recovers the curvature it was given, to 2e-6 (measured 4.9e-7)",
               fit.gradient, gradient, 2e-6 * gradient);
    expectNear("and the neutral axis, which it was not told", fit.neutralAxis, 0.5 * kH, 1e-6);
    // The rms is an average and can hide one bad point; the worst residual cannot, and
    // on a prismatic box under pure bending there is no bad point to hide.
    expectTrue("with no single Gauss point off the line either", fit.residualWorst < 100.0);
    expectTrue("and the worst is at least the rms, which is what a worst means",
               fit.residualWorst >= fit.residualRms);
    // **From an axis that is not the neutral one**, which is the case that has any
    // content: fitted about its own axis the axial term is zero and `about - axial /
    // gradient` is `about` whichever sign it carries. About the keel the axial term is
    // `-E kappa H / 2` and the answer has to come back at mid-height anyway.
    const section::BeamFit low = section::fitBeam(samples, 0.0);
    expectTrue("the fit about the keel ran", low.ok);
    expectNear("the neutral axis is found from the keel too", low.neutralAxis, 0.5 * kH, 1e-6);
    expectNear("with the same gradient", low.gradient, gradient, 2e-6 * gradient);
    expectTrue("and a non-zero axial term, so the sign of the shift is being tested",
               std::abs(low.axial) > 0.1 * gradient * kH);
    expectTrue("with nothing left over", fit.residualRms < 1e-6 * gradient * kH);
    // Vacuity: the field this is agreeing with has to be a real gradient, or "linear
    // in z" is a statement about a constant. The extreme fibre is 10.3 MPa.
    expectTrue("over a gradient worth fitting", gradient * 0.5 * kH > 1e7);

    // The two extreme fibres, which is what a section modulus is taken at. On a box
    // the deck and the keel are flat plates at one height apiece, so the worst and the
    // mean are the *same number* -- which is the negative control for shear lag: a
    // prismatic box has none, and a ferry does.
    const section::FibreStress deck = section::fibreStress(samples, true);
    const section::FibreStress keel = section::fibreStress(samples, false);
    expectTrue("there is a deck band and a keel band", deck.ok && keel.ok);
    expectNear("the deck carries E kappa H/2", deck.mean, gradient * 0.5 * kH,
               1e-3 * gradient * kH);
    expectNear("and the keel the same in compression", keel.mean, -gradient * 0.5 * kH,
               1e-3 * gradient * kH);
    // **Not exactly one, and the reason is the band and not the ship.** A fibre is a
    // line and a band has a depth, so under a gradient the worst point in it is higher
    // than the mean point in it by the band's own 2% of the section. That is 0.6% here
    // and it is the floor this instrument can resolve -- against 1.4 on the ferry, which
    // is what makes the ferry's number shear lag rather than arithmetic.
    std::printf("     box deck band: worst/mean %.6f over %zu samples\n",
                deck.worst / deck.mean, deck.samples);
    expectNear("a flat deck has no shear lag beyond its own band depth: measured 1.0058,"
               " asserted at 1.01", deck.worst / deck.mean, 1.0, 0.01);
    expectNear("nor a flat keel", keel.worst / keel.mean, 1.0, 0.01);
    expectTrue("and the two bands are on opposite sides of the axis", deck.z > keel.z);
    // **The sign, which a magnitude cannot carry.** Hogging tensions the deck and
    // compresses the keel; deck plating in compression buckles well below yield, so
    // reporting the worst fibre as a magnitude would lose the one thing that tells the
    // two failures apart. `girder.hpp` makes the same point about `M / Z`.
    expectTrue("the worst deck fibre is in tension", deck.worst > 0);
    expectTrue("and the worst keel fibre is in compression", keel.worst < 0);

    // **A cut with no depth has no fibres and is refused rather than answered.** Every
    // sample at one height makes the band the whole field, so a "deck stress" would come
    // back as the section's mean and read as a perfectly ordinary number. Nothing this
    // mesher builds produces one -- which is exactly why it has to be constructed here,
    // and why a mutation run found the guard untested.
    std::vector<section::StressSample> flat;
    for (int i = 0; i < 8; ++i)
        flat.push_back({Vec3{0.0, 0.1 * i, 1.0}, 1e-3, 1e6 * (i + 1)});
    expectTrue("a cut with no depth has no deck", !section::fibreStress(flat, true).ok);
    expectTrue("and no keel", !section::fibreStress(flat, false).ok);
    // Nor a beam: there is no gradient to find and the normal equations are singular.
    expectTrue("and no beam to fit through it", !section::fitBeam(flat, 0.0).ok);
    // The guard that says the list was otherwise usable, so "refused" is about the depth
    // and not about the samples.
    std::vector<section::StressSample> spread = flat;
    for (std::size_t i = 0; i < spread.size(); ++i) spread[i].at.z = 0.1 * i;
    expectTrue("the same samples spread over a depth do have a deck",
               section::fibreStress(spread, true).ok);
    expectTrue("and a beam", section::fitBeam(spread, 0.0).ok);
}

// **A wide deck does, and this is the fixture that has one.** The box's flanges are
// 2 m across on a 1 m depth and are fed by side plating at both edges, so a plane
// section is very nearly right; the reference ferry's strength deck is 20 m across on
// a 15 m depth and it is not. The claim is not a number -- it is that `fibreStress`
// separates a mean from a worst on a real hull and reports the same number on a
// prismatic one, which is the only way the ratio means anything.
void testFerryDeckStressIsNotOneNumber() {
    const StructuralMesh structure = ferryStructure();
    const section::Section hold = section::buildSection(structure, ferryParams());
    const StructuralMaterial material = ah36Steel();
    const HullGirderSection girder = hullGirderSection(structure, 0.0);
    expectTrue("the ferry section meshed", !hold.empty());

    section::BeamLoad bending;
    bending.curvature = 1e-6;
    bending.reference = girder.neutralAxis;
    std::vector<double> field;
    const section::BeamResponse bent = section::applyBeamLoad(hold, material, bending, &field);
    expectTrue("the ferry section takes a curvature: " + bent.problem, bent.ok);

    const std::vector<section::StressSample> samples =
        section::axialStress(hold, material, field, 0.0);
    expectTrue("a midship cut crosses hundreds of elements", samples.size() > 500);
    const section::BeamFit fit = section::fitBeam(samples, girder.neutralAxis);
    expectTrue("the fit ran", fit.ok);

    // The beam's own answer, from `hullGirderSection`, which shares no code with any
    // of this: `sigma = M (z - z_na) / I`, and `M / I` is `E kappa` if the two tiers
    // agree about `EI`. They agree to a fraction of a per cent amidships.
    const double beam = material.youngsModulus * bending.curvature;
    expectNear("the fitted gradient is E kappa", fit.gradient, beam, 0.02 * beam);
    expectNear("and the fitted neutral axis is the girder's", fit.neutralAxis,
               girder.neutralAxis, 0.15);
    // **And the residual is not zero, which is the whole point.** A beam has exactly
    // none of it. Asserted as a band rather than a value: below 1e-3 of the peak the
    // ferry would be behaving like a prismatic box and this test would be measuring
    // the box twice; above a tenth the fit would not be a beam at all and the
    // comparison against Tier 0 would be meaningless.
    expectTrue("a real hull does not carry a beam's stress field", fit.residualRms >
                                                                      1e-3 * std::abs(fit.peak));
    expectTrue("but it is still recognisably one", fit.residualRms < 0.2 * std::abs(fit.peak));
    // And the worst single point is worse than the rms by a real factor, which is what
    // says the departure is *local* -- shear lag crowding stress into the deck edge --
    // rather than a uniform offset the fit could have absorbed.
    expectTrue("the worst point is well above the rms", fit.residualWorst > 1.5 * fit.residualRms);

    const section::FibreStress deck = section::fibreStress(samples, true);
    expectTrue("the deck band has samples", deck.ok && deck.samples > 4);
    expectTrue("the worst deck stress is above its own mean", std::abs(deck.worst) >
                                                                  std::abs(deck.mean));
    expectTrue("and it is the *same sign*, so this is a spread and not two structures",
               deck.worst * deck.mean > 0);
}

// The half-open rule, which is the one `sectionElements` records a defect for.
//
// A plane on a frame station is served by the bay forward of it and not by both.
// Counting both would double the steel on the cut, and every average taken over it --
// the fit's gradient, the deck's mean -- is a *ratio*, so it would look perfectly
// correct while the volume behind it was twice what it should be.
void testAxialStressCountsASeamOnce() {
    const StructuralMesh structure = makeBox(kT, kT);
    const section::Section box = section::buildSection(structure, openBoxParams());
    const StructuralMaterial material = ah36Steel();
    section::BeamLoad axial;
    axial.strain = 1e-4;
    std::vector<double> field;
    const section::BeamResponse stretched = section::applyBeamLoad(box, material, axial, &field);
    expectTrue("the box solved", stretched.ok);

    const double bay = kL / kNx;
    const std::vector<section::StressSample> inside =
        section::axialStress(box, material, field, 4.0 * bay + 0.5 * bay);
    const std::vector<section::StressSample> onSeam =
        section::axialStress(box, material, field, 4.0 * bay);
    expectEqualCount("a seam is served by one bay, not two", onSeam.size(), inside.size());
    double insideVolume = 0, seamVolume = 0;
    for (const section::StressSample& s : inside) insideVolume += s.volume;
    for (const section::StressSample& s : onSeam) seamVolume += s.volume;
    expectNear("and carries the same steel", seamVolume, insideVolume, 1e-12 * insideVolume);
    // Vacuity: the plane really is on a seam, and a bay really is a distinct set of
    // elements -- otherwise "the same count" is the same elements twice.
    expectTrue("the box has bays to land between", kNx > 4);
    const std::vector<section::StressSample> next =
        section::axialStress(box, material, field, 5.0 * bay + 0.5 * bay);
    expectEqualCount("each bay contributes the same population", next.size(), inside.size());
}

// --- 11. A whole-ship model can answer a stress question ---------------------------
//
// `sectionDisplacement` is the one route from an assembled chain back to a mesh, and
// without it a Tier-1 ship reports resultants -- which is what Tier 0 already had.
// **What makes it exact is `reduction.hpp` property 1**, and it is worth checking
// rather than quoting: with zero modes the boundary response of a chain is the same
// solve as eliminating every interior, so a section's interface displacement is not an
// approximation of the monolith's, it is the monolith's.
void testChainRecoversTheStressItsSectionsCarry() {
    const StructuralMesh structure = makeBox(kT, kT);
    const StructuralMaterial material = ah36Steel();
    section::ChainParams params;
    // Tied, not open: a chain has to be one piece before `applyBeamLoad` will solve it,
    // and an open box is four plates. The corner ties cost `EA` nothing -- `section.hpp`
    // §5's table is 1.000000 either way -- so the patch test below is unaffected by the
    // choice and the component count is not.
    params.section = boxParams();
    params.reduce.modes = 0;
    params.reduce.cutoffFrequency = 0;
    for (int i = 0; i <= 4; ++i) params.station.push_back(kL * i / 4);
    const section::Chain chain = section::buildChain(structure, params);
    expectTrue("the chain of four is ready", chain.ready());
    expectEqual("in one piece", chain.components, 1);

    section::BeamLoad axial;
    axial.strain = 1e-4;
    std::vector<double> state;
    const section::BeamResponse stretched = section::applyBeamLoad(chain, axial, &state);
    expectTrue("the chain takes an axial load: " + stretched.problem, stretched.ok);
    expectEqualCount("and hands back its assembled state", state.size(),
                     static_cast<std::size_t>(chain.assembly.size()));

    // Every section of the chain, expanded onto its own mesh and asked for its stress.
    // The answer is the patch test again -- `E eps` everywhere -- and the point is that
    // it now comes through `componentState`, `recover` and the constraint expansion
    // rather than out of a single solve.
    const double exact = material.youngsModulus * axial.strain;
    double worst = 0;
    std::size_t total = 0;
    for (std::size_t i = 0; i < chain.section.size(); ++i) {
        const std::vector<double> field = section::sectionDisplacement(chain, i, state);
        expectEqualCount("section " + std::to_string(i) + " expands onto its own mesh",
                         field.size(), chain.section[i].nodeCount() * 3);
        const double middle = 0.5 * (chain.section[i].xFrom + chain.section[i].xTo);
        const std::vector<section::StressSample> samples =
            section::axialStress(chain.section[i], material, field, middle);
        total += samples.size();
        for (const section::StressSample& s : samples)
            worst = std::max(worst, std::abs(s.sigmaXX / exact - 1.0));
    }
    std::printf("     chain of four: worst |sigma/E eps - 1| over %zu Gauss points = %.3e\n",
                total, worst);
    expectTrue("every section of the assembled ship carries E eps, to 5e-9 (measured 5.8e-10 --"
               " three orders above the single section's 1.2e-11, which is the assembled"
               " Cholesky and not the recovery)", worst < 5e-9);
    expectTrue("over the whole length, not one section of it", total >= 4 * 96);

    // **The negative control for the recovery.** A state that is not this assembly's
    // comes back empty rather than as a plausible field: every index would otherwise be
    // in range and the modal content would quietly be missing.
    expectTrue("a foreign state is refused",
               section::sectionDisplacement(chain, 0, std::vector<double>(3, 0.0)).empty());
    expectTrue("and so is a section index the chain does not have",
               section::sectionDisplacement(chain, chain.section.size(), state).empty());
}

// --- 12. The member accounting reaches the second moment too -------------------------
//
// `missedMemberArea` was the whole of the accounting and it is not enough to correct
// an `I` comparison with. The ferry's girders are 4.4% of her area and 5.3% of her
// second moment -- they sit low in a double bottom, and a second moment is a lever arm
// squared -- so subtracting the area alone leaves the two tiers looking 5% apart
// amidships where they agree to 0.4%. That is the difference between "a mesher that
// cannot reach three girders" and "a mesher that is wrong".
void testMissedMembersAccountForTheirMomentsToo() {
    const StructuralMesh structure = ferryStructure();
    const section::Section hold = section::buildSection(structure, ferryParams());
    expectTrue("the ferry section meshed", !hold.empty());

    // `sectionElements` at the section's own mid-station, which is where a member's
    // share is evaluated. It shares no code with the mesher.
    double area = 0, first = 0, second = 0;
    for (const SectionElement& element : sectionElements(structure, 0.0)) {
        if (!element.stiffener) continue;
        area += element.area;
        first += element.area * element.height;
        second += element.ownSecondMoment + element.area * element.height * element.height;
    }
    expectNear("attached plus missed is the stiffener area of the cut",
               hold.attachedMemberArea + hold.missedMemberArea, area, 1e-9);
    expectNear("and its first moment about the baseline",
               hold.attachedMemberFirstMoment + hold.missedMemberFirstMoment, first, 1e-9 * first);
    expectNear("and its second moment about the baseline",
               hold.attachedMemberSecondMoment + hold.missedMemberSecondMoment, second,
               1e-9 * second);

    // What is missed is worth more of the second moment than of the area, and that is
    // the reason this exists rather than a curiosity: correcting one and not the other
    // is correcting the wrong one.
    const HullGirderSection girder = hullGirderSection(structure, 0.0);
    const double missedAboutAxis =
        hold.missedMemberSecondMoment - 2.0 * girder.neutralAxis * hold.missedMemberFirstMoment +
        girder.neutralAxis * girder.neutralAxis * hold.missedMemberArea;
    expectTrue("the girders are a few per cent of the ferry's area",
               hold.missedMemberArea / girder.area > 0.03);
    expectTrue("and more of her second moment",
               missedAboutAxis / girder.secondMoment > hold.missedMemberArea / girder.area);
    expectNear("2.46 m^4 of her 46.2, which is the figure `section.hpp` quotes",
               missedAboutAxis, 2.459, 0.05);

    // And with both corrections the two tiers agree, where with one of them they do
    // not. Asserting the *pair* is the point: the area agreement alone would pass on a
    // section that had lost the girders' whole Steiner term.
    const StructuralMaterial material = ah36Steel();
    section::BeamLoad bending;
    bending.curvature = 1e-6;
    bending.reference = girder.neutralAxis;
    const section::BeamResponse bent = section::applyBeamLoad(hold, material, bending);
    expectTrue("the section takes a curvature", bent.ok);
    const double effective = bent.bendingStiffness / material.youngsModulus;
    expectNear("EI against the girder's second moment less the girders' own",
               effective / (girder.secondMoment - missedAboutAxis), 1.0, 0.01);
    expectTrue("where subtracting only the area would be five per cent out",
               std::abs(effective / girder.secondMoment - 1.0) > 0.04);
}

// --- 13. The frames a beam cannot see ------------------------------------------------
//
// `hullGirderSection` drops every member with no extent along x, and it is right to:
// an athwartships member carries no longitudinal stress, so a frame, a deck beam or a
// bulkhead stiffener is worth nothing to a beam idealisation. **Tier 0 therefore scores
// a structure with no frames in it identically, and Tier 1 does not.**
//
// What is left is the transverse restraint a ring of frames puts on the section's own
// Poisson contraction: a strip that cannot contract in y and z carries more than
// `E eps` for the same strain. It is real, it is worth a few tenths of a per cent, and
// **the FEM is right about it and the beam is wrong** -- which is the one place in this
// comparison where a difference is not the finer model's error.
void testFramesRestrainThePoissonContraction() {
    const StructuralMesh structure = ferryStructure();
    StructuralMesh bare = structure;
    bare.members.clear();
    for (const StructuralMember& m : structure.members)
        if (std::abs(m.b.x - m.a.x) > 1e-9) bare.members.push_back(m);
    expectTrue("the ferry has athwartships members to remove",
               bare.members.size() < structure.members.size());
    expectTrue("and longitudinal ones to keep", !bare.members.empty());

    // Tier 0 scores the two identically, and that is checked rather than assumed --
    // it is the whole reason the difference below belongs to Tier 1.
    const HullGirderSection full = hullGirderSection(structure, 0.0);
    const HullGirderSection stripped = hullGirderSection(bare, 0.0);
    expectNear("Tier 0 cannot tell the two structures apart", stripped.area, full.area, 1e-12);
    expectNear("in area or in second moment", stripped.secondMoment, full.secondMoment, 1e-12);

    const StructuralMaterial material = ah36Steel();
    const auto effective = [&](const StructuralMesh& mesh) {
        const section::Section piece = section::buildSection(mesh, ferryParams());
        section::BeamLoad axial;
        axial.strain = 1e-6;
        axial.reference = full.neutralAxis;
        const section::BeamResponse stretchedOne = section::applyBeamLoad(piece, material, axial);
        expectTrue("the section takes an axial load: " + stretchedOne.problem, stretchedOne.ok);
        return stretchedOne.axialStiffness / material.youngsModulus /
               (full.area - piece.missedMemberArea);
    };
    const double with = effective(structure), without = effective(bare);
    std::printf("     EA against Tier 0 less the members it could not attach: %+.4f%% with the"
                " frames, %+.4f%% without, so they are worth %+.4f%%\n", 100.0 * (with - 1.0),
                100.0 * (without - 1.0), 100.0 * (with - without));

    // The direction is the claim: restraint can only stiffen. The size is reported.
    expectTrue("the frames stiffen the section rather than soften it", with > without);
    // Measured 0.44% on this two-bay window, and it rises with length towards 0.5% as
    // the free cut planes become a smaller share of the section. Asserted as a band,
    // because a value would be asserting the ferry's frame spacing.
    expectTrue("by a few tenths of a per cent, not by a factor",
               with - without > 1e-3 && with - without < 1e-2);
    // **Without them the two tiers agree to a tenth of a per cent, and that is what
    // says the rest of the accounting is right.** With them Tier 1 reads high, and
    // reading high is the direction a restraint has to move it.
    expectTrue("and without them the beam and the mesh agree to a tenth of a per cent",
               std::abs(without - 1.0) < 1e-3);
    expectTrue("while with them the mesh reads above the beam", with > 1.0);
}

}  // namespace

void runSectionTests() {
    std::printf("\n=== Tier-1 section mesher ===\n");
    testBoxMesh();
    testCollapsedPanelsMeshAsWedges();
    testWedgeApexIsCountedOnce();
    testJunctionTieIsAPartitionOfUnity();
    testPlaneTieLiesInThePlane();
    testWeldIsADistance();
    testReversedWindingIsOriented();
    testBoxSectionProperties();
    testJunctionTieClosesTheCell();
    testThicknessSeam();
    testMemberRunsStopAtAThicknessStep();
    testResolutionConvergence();
    testRefusals();
    testBandwidthReducingOrder();
    testFreeEdgesInFreshAirAreNotJunctions();
    testSubstructureConsumesTheSection();
    testFerrySection();
    testFerryJunctionMovesTheLowestFrequency();
    testFerryMembersAreWorthWhatTheSectionSays();
    testFerryFibresAndTheirMass();
    testFerryBulkheadCutIsNotASeam();
    testInterfaceIsChosenOnTheMidSurface();
    testChainOfSectionsIsTheWholeSection();
    testChainCutPlanesUntieTheJunctions();
    testPlaneTieSplitsTheMastersThickness();
    testDeckChainSplitsTheSameThicknessTheMonolithDoes();
    testPlaneTieRefusals();
    testPlaneTieSurvivesReversedWinding();
    testChainPlaneTiesAreExactOnTheAssembly();
    testChainFrequencyReachesTheMonolith();
    testChainFrequencySeesTheJoins();
    testMesherReachesTheEndsOfTheShip();
    testJunctionWeightLimitRefusesTheNodeAndNotTheSection();
    testFerryChainConservesTheStructure();
    testHaloReachesAsFarAsTheWeld();
    testNodePositionsDoNotDependOnWhereTheSectionWasCut();
    testHaloMakesTheCutPlanesCoincide();
    testAxialStressIsThePatchTest();
    testAxialStressUnderCurvatureIsTheBeamsOwnLine();
    testAxialStressCountsASeamOnce();
    testFerryDeckStressIsNotOneNumber();
    testChainRecoversTheStressItsSectionsCarry();
    testMissedMembersAccountForTheirMomentsToo();
    testFramesRestrainThePoissonContraction();
}
