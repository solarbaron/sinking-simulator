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

    // Four corner lines, 17 stations each, of which the two on the cut planes cannot
    // be tied: an interface degree of freedom is prescribed, and prescribing it and
    // deriving it are two different claims on one unknown.
    expectEqual("every corner station off the cut planes is tied", tied.junctionTies,
                4 * (kNx - 1));
    expectEqual("and the ones on them are counted, not silently dropped",
                tied.junctionsOnInterface, 4 * 2 * 2);
    expectEqual("nothing was refused as a chain", tied.junctionsChained, 0);
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
}

// --- 1c. What the tie is, before what it does ------------------------------------
//
// Two identities, checked on the constraint set itself rather than through a solve,
// because a tie that is wrong in either is wrong in a way an energy comparison
// would report as a small stiffness change.

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

    section::SectionParams backwards = boxParams();
    backwards.xTo = backwards.xFrom;
    expectTrue("cut planes out of order are refused",
               section::buildSection(structure, backwards).empty());

    section::SectionParams noSubdivision = boxParams();
    noSubdivision.subdivision = 0;
    expectTrue("a subdivision below one is refused",
               section::buildSection(structure, noSubdivision).empty());

    section::SectionParams noRoles = boxParams();
    noRoles.shell = noRoles.deck = noRoles.bulkhead = false;
    expectTrue("a section of no roles at all is refused",
               section::buildSection(structure, noRoles).empty());

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
    expectTrue("most of a two-bay section's junction nodes are on a cut plane",
               hold.junctionsOnInterface > hold.junctionTies);
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
    // 0.44% measured, and it is not noise: the frames and deck beams restrain the
    // section's Poisson contraction, which a beam model has no way to represent. It
    // is worth 0.52% of the area, measured by omitting them.
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

}  // namespace

void runSectionTests() {
    std::printf("\n=== Tier-1 section mesher ===\n");
    testBoxMesh();
    testJunctionTieIsAPartitionOfUnity();
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
}
