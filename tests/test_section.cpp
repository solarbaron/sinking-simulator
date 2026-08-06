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

// **What a cut plane costs, which `EA` cannot see and torsion can.**
//
// A junction node on a cut plane is an interface degree of freedom, and one of those
// is prescribed rather than derived -- so it cannot also be tied. Cutting a length
// into N pieces turns N-1 interior stations into interfaces and unties them. `EA` is
// exact at every N regardless; `GJ` falls monotonically, and by 19% at N = 8.
void testChainCutPlanesUntieTheJunctions() {
    std::printf("\n--- section: what an interior cut plane costs the junctions ---\n");
    const StructuralMesh structure = makeBox(kT, kT);
    const StructuralMaterial steel = ah36Steel();
    const section::Section whole = section::buildSection(structure, boxParams());
    const section::TorsionResponse wholeTwist = section::applyTwist(whole, steel, 1e-4, kH / 2);
    const double bredt = shearModulus() * kBoxTorsion;
    expectTrue("the monolithic box is closed in torsion",
               wholeTwist.torsionalStiffness > 0.9 * bredt);

    double previousTorsion = 2.0 * wholeTwist.torsionalStiffness;
    double previousTied = 2.0 * whole.tiedEdges;
    for (int sections : {1, 2, 4, 8}) {
        const section::Chain chain = section::buildChain(structure, boxChain(sections, true));
        expectTrue("the tied chain built", chain.ready());
        expectEqual("and is one piece of ship", chain.components, 1);
        section::BeamLoad axial;
        axial.strain = 1e-4;
        const section::BeamResponse stretched = section::applyBeamLoad(chain, axial);
        const section::TorsionResponse twisted = section::applyTwist(chain, 1e-4, kH / 2);
        expectTrue("both solves ran", stretched.ok && twisted.ok);
        std::printf("     %d sections: tied %5.1f m of %.1f, EA %+.2e of the whole, GJ %+.3f%%,"
                    " %.4f of Bredt\n",
                    sections, chain.tiedEdges, chain.junctionEdges,
                    stretched.axialStiffness / (steel.youngsModulus * kBoxArea) - 1,
                    100.0 * (twisted.torsionalStiffness / wholeTwist.torsionalStiffness - 1),
                    twisted.torsionalStiffness / bredt);

        // The whole point: axial says nothing.
        expectNear("EA is the closed form at every N", stretched.axialStiffness,
                   steel.youngsModulus * kBoxArea, kAxialTolerance * steel.youngsModulus * kBoxArea);
        // And torsion says it every time.
        expectTrue("every extra cut plane unties another station", chain.tiedEdges < previousTied);
        expectTrue("and costs torsional stiffness", twisted.torsionalStiffness < previousTorsion);
        previousTied = chain.tiedEdges;
        previousTorsion = twisted.torsionalStiffness;
        if (sections == 1)
            // One section is the monolith reduced and reassembled, so it is the
            // control for everything above: it must reproduce it exactly.
            expectNear("a chain of one is the section it was cut from",
                       twisted.torsionalStiffness, wholeTwist.torsionalStiffness,
                       1e-9 * wholeTwist.torsionalStiffness);
        if (sections == 8)
            // Measured at -19.20%. Asserted as a band rather than a point because it
            // is a property of this box's bay count, but a band tight enough that a
            // chain which had stopped losing ties would fail it.
            expectTrue("eight sections lose about a fifth of the torsional stiffness",
                       twisted.torsionalStiffness < 0.85 * wholeTwist.torsionalStiffness &&
                           twisted.torsionalStiffness > 0.75 * wholeTwist.torsionalStiffness);
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
        expectNear(std::string(station.label) + " carries the steel those panels have",
                   piece.plateMass, mass, 1e-9 * mass);
        expectTrue(std::string(station.label) + " is a real piece of ship", area > 100.0);

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

    // --- What the reach exposes, which is not the reach ------------------------------
    //
    // Cutting a window in two conserves the plating exactly and the **stiffener steel
    // only where the plating does not change thickness**, and reaching the ends is
    // what makes that visible. `nodeThickness` is an area-weighted mean over the
    // sub-quads *inside the section*, so a station where the strake steps carries one
    // thickness to the section aft of it, another to the section forward of it, and
    // the mean of the two to a section that spans it. A member run stops at a
    // thickness change (§3), so the spanning section stops runs the two halves never
    // see, and the seam node's run of one is dropped.
    //
    // **It is not new and it is not the collapsed element.** It is the nodal-thickness
    // twin of `section.hpp` §6 note 1's nodal normal, and it is measurable inside the
    // range that always worked -- 1.9% at x = -24 .. -19.2. The ends are where it
    // stops being a rounding: at the bow shoulder the spanning section carries **25%
    // less** stiffener steel than the two halves it is made of.
    const auto cut = [&](double from, double to) {
        section::SectionParams p;
        p.xFrom = from;
        p.xTo = to;
        p.junctions = false;
        return section::buildSection(structure, p);
    };
    struct Split {
        const char* label;
        double from, middle, to;
        bool uniform;  // is the plating the same thickness across the shared station
    };
    for (const Split& split : {Split{"amidships", -7.2, -2.4, 2.4, true},
                               Split{"old reach", -24.0, -21.6, -19.2, false},
                               Split{"bow shoulder", 40.8, 43.2, 45.6, false}}) {
        const section::Section aft = cut(split.from, split.middle);
        const section::Section forward = cut(split.middle, split.to);
        const section::Section whole = cut(split.from, split.to);
        // The plating, the elements and the wedges conserve whatever the thickness
        // does. Only the members are at issue.
        expectEqualCount(std::string(split.label) + ": two sections have one section's elements",
                         aft.elementCount() + forward.elementCount(), whole.elementCount());
        expectNear(std::string(split.label) + ": and its plating", aft.plateMass + forward.plateMass,
                   whole.plateMass, 1e-9);
        expectEqual(std::string(split.label) + ": and its collapsed elements",
                    aft.collapsedElements + forward.collapsedElements, whole.collapsedElements);

        const double halves = aft.memberMass + forward.memberMass;
        const double shortfall = 1.0 - whole.memberMass / halves;
        const int splitHalves =
            aft.memberRunsSplitByThickness + forward.memberRunsSplitByThickness;
        std::printf("     %-13s stiffener steel: halves %8.1f kg, one piece %8.1f (%+.1f%%);"
                    " runs stopped at a seam %d against %d\n",
                    split.label, halves, whole.memberMass, -100.0 * shortfall,
                    whole.memberRunsSplitByThickness, splitHalves);
        if (split.uniform) {
            expectNear(std::string(split.label) + ": the stiffener steel conserves exactly",
                       whole.memberMass, halves, 1e-9 * halves);
            expectEqual(std::string(split.label) + ": and no extra run is stopped",
                        whole.memberRunsSplitByThickness, splitHalves);
        } else {
            // The direction is the whole point: the *spanning* section is the one
            // that loses, because it is the one that sees the seam. A test that
            // asserted only "they differ" would pass on the opposite sign, which
            // would be double-counted steel rather than dropped steel.
            expectTrue(std::string(split.label) + ": the spanning section stops more runs",
                       whole.memberRunsSplitByThickness > splitHalves);
            expectTrue(std::string(split.label) + ": and is the one short of steel",
                       whole.memberMass < halves);
            expectTrue(std::string(split.label) + ": by a real fraction of it", shortfall > 0.01);
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

// **The interface is not coincident by construction, and the parallel middle body is
// what hides it.** See `section.hpp` §6 note 1. A node is the mid-surface offset by
// `t/2` along a normal this file averages over the sub-quads *inside* the section, so
// two sections cut on the same station disagree about where that node is wherever the
// hull is not prismatic. It is 3.4 µm at the outermost plane a chain can still be
// built on -- and 3.4 µm is three orders of magnitude above the default tolerance.
void testFerryChainInterfaceIsNotCoincidentByConstruction() {
    std::printf("\n--- section: where two sections stop agreeing about the cut plane ---\n");
    const StructuralMesh structure = ferryStructure();
    section::ChainParams params;
    params.section.subdivision = 1;
    params.section.members = false;
    params.section.junctions = false;
    params.reduce.modes = 0;
    params.reduce.cutoffFrequency = 0;
    params.station = {-26.4, -21.6, -16.8};

    const section::Chain strict = section::buildChain(structure, params);
    expectEqualCount("the chain meshed two sections", strict.section.size(), 2u);
    expectEqualCount("one interior plane", strict.unmatched.size(), 1u);
    const std::size_t planeDof = 3 * strict.section[0].forwardNodes.size();
    expectEqualCount("of 1170 boundary DOF", planeDof, 1170u);
    std::printf("     at 1e-9: %zu of %zu unmatched, worst matched gap %.4e m\n",
                strict.unmatched[0], planeDof, strict.worstGap);
    expectTrue("at the default tolerance a large part of the plane finds no partner",
               strict.unmatched[0] > planeDof / 4);
    expectTrue("so the chain is refused rather than solved", !strict.ready());

    section::ChainParams loose = params;
    loose.matchTolerance = 1e-5;
    const section::Chain joined = section::buildChain(structure, loose);
    expectTrue("a tolerance the size of the disagreement joins the whole plane",
               joined.ready());
    expectEqualCount("with nothing left over", joined.unmatched[0], 0u);
    expectEqualCount("and the whole plane shared", joined.shared[0], planeDof);
    std::printf("     at 1e-5: 0 unmatched, worst matched gap %.4e m\n", joined.worstGap);
    // The assembly's own account of it, which is a different measurement of the same
    // thing: `Chain::worstGap` comes from the interface map, `worstMergedGap` from
    // comparing the identities of the DOF that actually landed on one assembled row.
    // Nothing else in the suite ever merges two DOF that are not exactly coincident,
    // so without this the two are indistinguishable from zero and from each other.
    expectNear("and the assembly reports the same gap on the rows it merged",
               joined.assembly.worstMergedGap, joined.worstGap, 1e-15);
    expectEqual("with no axis crossed", joined.assembly.axisDisagreements, 0);

    // Measured at 3.4021e-06 m. Bracketed rather than asserted at a point, because it
    // is a property of this hull's shoulder -- but bracketed tightly enough that a
    // mesher which had started averaging its normals across the cut would fail it,
    // which is the fix §6 names.
    expectTrue("and the disagreement is microns, not the plate thickness",
               joined.worstGap > 1e-6 && joined.worstGap < 1e-5);
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
    // Measured at 3.9e-12 relative. Asserted at 1e-9 -- three decades above the
    // measurement and six below anything a lost interface could hide in.
    expectNear("two ferry sections carry the torque one ferry section carries",
               chainTwist.torsionalStiffness, monoTwist.torsionalStiffness,
               1e-9 * monoTwist.torsionalStiffness);
    // Vacuity: an open cell and a closed one differ by orders of magnitude, so this
    // is a comparison of something rather than of two zeros.
    expectTrue("and it is a real torque", monoTwist.torsionalStiffness > 1e11);
}

}  // namespace

void runSectionTests() {
    std::printf("\n=== Tier-1 section mesher ===\n");
    testBoxMesh();
    testCollapsedPanelsMeshAsWedges();
    testWedgeApexIsCountedOnce();
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
    testChainOfSectionsIsTheWholeSection();
    testChainCutPlanesUntieTheJunctions();
    testChainFrequencySeesTheJoins();
    testMesherReachesTheEndsOfTheShip();
    testJunctionWeightLimitRefusesTheNodeAndNotTheSection();
    testFerryChainConservesTheStructure();
    testFerryChainInterfaceIsNotCoincidentByConstruction();
}
