// SPDX-License-Identifier: MIT
//
// Validation of the Tier-2 zone: `StructuralMesh` -> solid-shell elements -> an
// explicit solve -> torn panels.
//
// The element and the return map are already validated against closed forms in
// `test_solid_shell.cpp` and `test_plasticity.cpp`, so nothing here re-checks
// them. What is new, and therefore what is aimed at, is:
//
//   * **the mesher**, which carries the element's binding limit (`07-fem-spike-
//     findings.md` §6 limit 1 -- the element must be prismatic through its
//     thickness). Every geometric claim here is an *identity* on flat plating,
//     because on flat plating the extrusion is exact and "nearly prismatic" would
//     hide a defect that only shows on a curve;
//   * **the explicit solver**, checked by conservation and by invariance rather
//     than against a reference run: nothing at rest starts moving, a rigid
//     translation carries no force, and the energy in equals the energy out;
//   * **the answer**, against `indentation.hpp` on the same bay. Two models built
//     from different physics agreeing to a factor is worth more than either
//     agreeing with itself -- and the factor is *predicted* here, from the
//     hardening curve and the plane-strain constraint the membrane model does not
//     have, rather than accepted.
//
// The expensive runs are in `tools/zone_probe`, which `verify.sh full` runs. What
// is left here costs a few seconds, and every one of those seconds is spent on
// a check that would fail against a plausible wrong implementation.
#include "engine/sim/breach.hpp"
#include "engine/sim/indentation.hpp"
#include "engine/sim/plasticity.hpp"
#include "engine/sim/scantlings.hpp"
#include "engine/sim/solid_shell.hpp"
#include "engine/sim/zone.hpp"
#include "engine/core/jobs.hpp"
#include "game/prototype/ferry.hpp"
#include "harness.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <string>
#include <utility>
#include <vector>

using namespace sim;
using testing::expectNear;
using testing::expectEqual;
using testing::expectTrue;

namespace {

// `expectEqual` takes long long; sizes are size_t. One place to cast rather than
// twenty.
void expectEqualCount(const std::string& what, std::size_t got, std::size_t want) {
    testing::expectEqual(what, static_cast<long long>(got), static_cast<long long>(want));
}

// --- Structures to mesh --------------------------------------------------------
//
// Built as `StructuralMesh` directly rather than through `makeStructuralMesh`,
// because the point of most of these is a geometry whose answer is known exactly:
// a flat rectangle's area is `lengthX * spanY` and nothing else.

// A flat strip: `lengthX` along x, `spanY` along y, divided into nx by ny panels.
// The tent forms across y, so `spanY` is the membrane model's span and `lengthX`
// its struck width.
StructuralMesh flatStrip(double lengthX, double spanY, double thickness, int nx, int ny,
                         bool stiffened = false) {
    StructuralMesh mesh;
    mesh.materials.push_back(ah36Steel());
    mesh.frameSpacing = lengthX;
    for (int i = 0; i < nx; ++i)
        for (int j = 0; j < ny; ++j) {
            PlatePanel p;
            const double x0 = -0.5 * lengthX + lengthX * i / nx;
            const double x1 = -0.5 * lengthX + lengthX * (i + 1) / nx;
            const double y0 = -0.5 * spanY + spanY * j / ny;
            const double y1 = -0.5 * spanY + spanY * (j + 1) / ny;
            p.corner[0] = {x0, y0, 0};
            p.corner[1] = {x1, y0, 0};
            p.corner[2] = {x1, y1, 0};
            p.corner[3] = {x0, y1, 0};
            p.thickness = thickness;
            p.material = 0;
            p.role = PanelRole::Shell;
            mesh.panels.push_back(p);
        }
    // Two longitudinals at the quarter points, on panel seams, exactly where
    // `makeStructuralMesh` puts them -- and clear of the middle, so a punch on the
    // centreline is not landing on a support.
    if (stiffened)
        for (double side : {-0.25, 0.25}) {
            StructuralMember member;
            member.a = {-0.5 * lengthX, side * spanY, 0.0};
            member.b = {0.5 * lengthX, side * spanY, 0.0};
            member.rise = {0, 0, 1};
            member.profile = flatBar(0.200, 0.010);
            member.attachedPlateThickness = thickness;
            member.role = MemberRole::Longitudinal;
            mesh.members.push_back(member);
        }
    return mesh;
}

// A cylindrical shell of radius R, spanning `sweep` radians in `bands` flat
// chords, `nx` panels along its axis. The panels chord the surface exactly as
// `makeStructuralMesh`'s do, so the facet angle is `sweep / bands` and the
// mesher's normal spread has a closed form to be checked against.
StructuralMesh cylinder(double radius, double sweep, double lengthX, double thickness, int nx,
                        int bands) {
    StructuralMesh mesh;
    mesh.materials.push_back(ah36Steel());
    mesh.frameSpacing = lengthX / nx;
    for (int i = 0; i < nx; ++i)
        for (int j = 0; j < bands; ++j) {
            const double x0 = -0.5 * lengthX + lengthX * i / nx;
            const double x1 = -0.5 * lengthX + lengthX * (i + 1) / nx;
            const double a0 = -0.5 * sweep + sweep * j / bands;
            const double a1 = -0.5 * sweep + sweep * (j + 1) / bands;
            const Vec3 p0{x0, radius * std::sin(a0), radius * std::cos(a0)};
            const Vec3 p1{x1, radius * std::sin(a0), radius * std::cos(a0)};
            const Vec3 p2{x1, radius * std::sin(a1), radius * std::cos(a1)};
            const Vec3 p3{x0, radius * std::sin(a1), radius * std::cos(a1)};
            PlatePanel p;
            p.corner[0] = p0;
            p.corner[1] = p1;
            p.corner[2] = p2;
            p.corner[3] = p3;
            p.thickness = thickness;
            p.material = 0;
            p.role = PanelRole::Shell;
            mesh.panels.push_back(p);
        }
    return mesh;
}

// A single panel that is a *trapezoid*: one end `spanY` wide, the other `taper`
// times that. Every rectangle in this file is a parallelogram, and a parallelogram
// is the one quad on which halving a diagonal and doubling it gives the right area
// -- the same trap `test_scantlings.cpp` records finding in `PlatePanel::area()`.
StructuralMesh trapezoid(double lengthX, double spanY, double taper, double thickness) {
    StructuralMesh mesh;
    mesh.materials.push_back(ah36Steel());
    PlatePanel p;
    p.corner[0] = {-0.5 * lengthX, -0.5 * spanY, 0};
    p.corner[1] = {0.5 * lengthX, -0.5 * taper * spanY, 0};
    p.corner[2] = {0.5 * lengthX, 0.5 * taper * spanY, 0};
    p.corner[3] = {-0.5 * lengthX, 0.5 * spanY, 0};
    p.thickness = thickness;
    p.material = 0;
    p.role = PanelRole::Shell;
    mesh.panels.push_back(p);
    return mesh;
}

// Two flat panels meeting along y = 0 at an angle, with the second `ratio` times
// the first's area. The nodal normal on the seam is the *area-weighted* mean of
// the two face normals, and on a graded pair the weighted and unweighted answers
// differ -- which is the only condition under which the weighting is testable.
StructuralMesh gradedFold(double lengthX, double wide, double narrow, double angle,
                          double thickness) {
    StructuralMesh mesh;
    mesh.materials.push_back(ah36Steel());
    const Vec3 tilt{0.0, std::cos(angle), std::sin(angle)};
    const Vec3 corners[2][4] = {
        {{-0.5 * lengthX, -wide, 0}, {0.5 * lengthX, -wide, 0},
         {0.5 * lengthX, 0, 0}, {-0.5 * lengthX, 0, 0}},
        {{-0.5 * lengthX, 0, 0}, {0.5 * lengthX, 0, 0},
         {0.5 * lengthX, narrow * tilt.y, narrow * tilt.z},
         {-0.5 * lengthX, narrow * tilt.y, narrow * tilt.z}},
    };
    for (const auto& quad : corners) {
        PlatePanel p;
        for (int c = 0; c < 4; ++c) p.corner[c] = quad[c];
        p.thickness = thickness;
        p.material = 0;
        p.role = PanelRole::Shell;
        mesh.panels.push_back(p);
    }
    return mesh;
}

zone::MeshParams flatParams(int subdivision) {
    zone::MeshParams params;
    params.radius = 1e3;
    params.subdivision = subdivision;
    params.outward = {0, 0, 1};
    return params;
}

// A material that cannot tear, so the flow rule is exercised on its own -- the
// same device `test_plasticity.cpp` uses with `kNeverFails`.
plasticity::Material unbreakableSteel() {
    plasticity::Material material = plasticity::shipSteel();
    material.failure.uniformStrain = 1e9;
    material.failure.fractureStrain = 1e9;
    return material;
}

// And one that cannot yield either, so the *plastic* code path has to reproduce
// the elastic one exactly.
//
// 20 GPa rather than something astronomical, and the reason is a trap worth
// recording: `elementPlasticUpdate` scales its enhanced-strain convergence
// tolerance by `sigma_y * V`, the element's yield energy, because that is the only
// scale in the problem that does not move as the element unloads. Handing it an
// absurd yield strength therefore loosens the tolerance in proportion, and at
// 1e14 Pa the Newton stopped two orders early and the two paths disagreed by 1.9%
// -- which reads exactly like a defect in the strain-energy expression.
plasticity::Material rigidPlasticSteel() {
    plasticity::Material material = unbreakableSteel();
    material.flow = plasticity::linearHardening(2.0e10, 0.0);
    return material;
}

// --- 1. The mesher, on flat plating, where every answer is exact ---------------

void testTheMesherIsExactOnFlatPlating() {
    std::printf("\n   meshing flat plating: the answers are identities, not tolerances\n");
    const double lengthX = 1.6, spanY = 0.8, thickness = 0.020;
    const StructuralMesh strip = flatStrip(lengthX, spanY, thickness, 8, 4);

    for (int subdivision : {1, 2, 3}) {
        const zone::Patch patch = zone::buildPatch(strip, {0, 0, 0}, flatParams(subdivision));
        const std::size_t panels = strip.panels.size();
        const auto want = static_cast<std::size_t>(subdivision) *
                          static_cast<std::size_t>(subdivision) * panels;
        expectEqualCount("every panel becomes subdivision^2 elements", patch.elementCount(), want);

        // Conforming: the nodes are *welded*, so the count is the grid's and not
        // the sum of the panels' own corners. Without the weld it would be
        // 4 * subdivision^2 * panels * 2, which is the negative control below.
        const auto across = static_cast<std::size_t>(8 * subdivision + 1);
        const auto down = static_cast<std::size_t>(4 * subdivision + 1);
        expectEqualCount("and the nodes are shared between them", patch.nodeCount(),
                        across * down * 2);

        expectNear("the meshed area is the plating's area, exactly", patch.area, lengthX * spanY,
                   1e-12 * lengthX * spanY);
        expectNear("and the mass is that area times thickness times density", patch.mass,
                   lengthX * spanY * thickness * ah36Steel().density, 1e-9 * patch.mass);
        expectNear("the patch carries the plating's thickness", patch.thickness, thickness, 1e-15);

        // The identity that says the elements are *exactly* prismatic rather than
        // nearly so. On flat plating every nodal normal is the same vector, so the
        // face offset is identically zero -- and 90 * 0^2 is zero spurious
        // stiffness. Anything else means the extrusion is not along one direction.
        expectTrue("on flat plating the faces are exactly parallel",
                   patch.worstNormalSpread == 0.0);
        expectEqual("so no element is distorted", patch.distortedElements, 0LL);
        expectTrue("and no element is inverted", patch.worstJacobian > 0);

        // Each element measures its own thickness as the plating's, to rounding.
        // `elementSize` takes it as volume over mid-surface area, so this catches
        // an extrusion that is the right length in the wrong direction.
        double worstThickness = 0;
        for (std::size_t e = 0; e < patch.elementCount(); ++e) {
            double nodes[solidshell::kDof];
            patch.mesh.gather(e, patch.mesh.position, nodes);
            double inPlane = 0, measured = 0;
            solidshell::elementSize(nodes, &inPlane, &measured);
            worstThickness = std::max(worstThickness, std::abs(measured - thickness));
        }
        expectTrue("every element measures the plating's own thickness",
                   worstThickness < 1e-14 * thickness);
    }

    // The negative control for the weld: two panels that do not touch cannot share
    // a node, so the count has to be the unwelded one. Without this the "nodes are
    // shared" check above would pass on a mesher that welded everything within a
    // metre of everything else.
    StructuralMesh apart = flatStrip(0.4, 0.4, thickness, 1, 1);
    PlatePanel far = apart.panels[0];
    for (Vec3& corner : far.corner) corner.x += 10.0;
    apart.panels.push_back(far);
    const zone::Patch split = zone::buildPatch(apart, {0, 0, 0}, flatParams(1));
    expectEqualCount("a panel ten metres away is not welded on, and is not reached either",
                    split.elementCount(), std::size_t{1});

    // The struck panel decides the patch, so a strike at the far panel finds that
    // one instead. Without this the previous check could pass on a mesher that
    // always returned the first panel.
    const zone::Patch other = zone::buildPatch(apart, {10.0, 0, 0}, flatParams(1));
    expectEqualCount("a strike at the far panel meshes the far panel", other.elementCount(),
                    std::size_t{1});
    expectTrue("and it is a different panel", other.panels != split.panels);
}

// --- 2. Curvature: what it costs, and where the cure has to be applied ---------

void testTheMesherReportsWhatCurvatureCosts() {
    std::printf("\n   curvature: the face offset is the facet angle, and it is the"
                " Scantlings' to fix\n");
    const double radius = 4.0, sweep = 1.2, lengthX = 2.0, thickness = 0.020;

    std::printf("     %6s %6s %14s %14s %12s\n", "bands", "sub", "facet angle", "face offset",
                "excess stiff");
    double previousSpread = 0, previousStiffness = 0;
    for (int bands : {4, 8, 16}) {
        const StructuralMesh shell = cylinder(radius, sweep, lengthX, thickness, 4, bands);
        zone::MeshParams params = flatParams(2);
        params.outward = {0, 0, 1};
        params.normalSpreadWarning = 1e-3;
        const zone::Patch patch = zone::buildPatch(shell, {0, 0, radius}, params);
        const double facet = sweep / bands;
        std::printf("     %6d %6d %14.4f %14.4f %11.1f%%\n", bands, 2, facet,
                    patch.worstNormalSpread, 100.0 * patch.spuriousStiffness);

        // A node on a seam takes the mean of the two facet normals, so an element
        // touching one seam has two corners at the facet normal and two half way
        // to the neighbour's: the spread is a quarter of the facet angle. That is
        // a closed form, and it is what makes this a test of the extrusion rather
        // than a record of whatever it produced.
        expectNear("the face offset is a quarter of the facet angle", patch.worstNormalSpread,
                   0.25 * facet, 0.06 * facet);
        if (previousSpread > 0) {
            expectNear("halving the girth band halves it", patch.worstNormalSpread,
                       0.5 * previousSpread, 0.03 * previousSpread);
            // And quarters the penalty, because the parasitic strain is linear in the
            // distortion and the energy is its square. A penalty that moved linearly
            // would be reporting the offset twice rather than its consequence.
            expectNear("and quarters the spurious stiffness", patch.spuriousStiffness,
                       0.25 * previousStiffness, 0.06 * previousStiffness);
        }
        previousSpread = patch.worstNormalSpread;
        previousStiffness = patch.spuriousStiffness;
        expectTrue("and the patch says so", !patch.problems.empty());
        expectTrue("naming how many elements are affected", patch.distortedElements > 0);
    }

    // **Subdividing does not help**, and that is the whole point: a panel is a flat
    // facet, so all the turning is at the seam whatever the subdivision, and the
    // cure is a finer girth layout in `Scantlings` rather than a finer zone. A
    // mesher that reported a falling offset under subdivision would be reporting
    // a cure it had not applied.
    const StructuralMesh shell = cylinder(radius, sweep, lengthX, thickness, 4, 8);
    double coarse = 0;
    for (int subdivision : {1, 2, 4}) {
        zone::MeshParams params = flatParams(subdivision);
        params.normalSpreadWarning = 1e-3;
        const zone::Patch patch = zone::buildPatch(shell, {0, 0, radius}, params);
        if (subdivision == 1) coarse = patch.worstNormalSpread;
        // Subdivision 1 puts the whole element between two seams, so its spread is
        // half the facet angle; 2 and beyond leave it at a quarter and stay there.
        if (subdivision > 1)
            expectNear("subdividing the zone does not flatten the facet",
                       patch.worstNormalSpread, 0.25 * sweep / 8, 0.06 * sweep / 8);
    }
    expectTrue("and a single element per panel is worse still, not better",
               coarse > 0.3 * sweep / 8);

    // A flat patch is the zero of the same measurement, which is the guard against
    // this whole test passing on a mesher that reports a constant.
    const zone::Patch flat =
        zone::buildPatch(flatStrip(2.0, 2.0, thickness, 4, 4), {0, 0, 0}, flatParams(2));
    expectTrue("flat plating reports no offset at all", flat.worstNormalSpread == 0.0);
    expectTrue("and raises nothing", flat.problems.empty());
}

// --- 3. Where a zone stops ------------------------------------------------------

void testTheZoneStopsAtASeamAndAtAFold() {
    std::printf("\n   where a zone stops: thickness seams, folds, and the radius\n");

    // Two strakes butted together, 12 mm and 20 mm.
    StructuralMesh strakes = flatStrip(1.6, 1.6, 0.012, 4, 4);
    for (PlatePanel& p : strakes.panels)
        if (p.centroid().y > 0) p.thickness = 0.020;

    const zone::Patch stopped = zone::buildPatch(strakes, {0, -0.6, 0}, flatParams(1));
    expectEqualCount("the zone contains only the struck strake", stopped.elementCount(),
                    std::size_t{8});
    expectNear("at the struck strake's thickness", stopped.thickness, 0.012, 1e-15);
    expectTrue("and says it stopped at a seam", !stopped.problems.empty());
    for (int index : stopped.panels)
        expectTrue("no panel of the other thickness is in it",
                   strakes.panels[static_cast<std::size_t>(index)].thickness == 0.012);

    // Told to cross it, it meshes the whole thing at the *struck* thickness --
    // because one patch has one extrusion distance -- and says that too. Reporting
    // the truncation and reporting the thickness error are different messages and
    // both have to exist, or a caller cannot tell which trade they got.
    zone::MeshParams crossing = flatParams(1);
    crossing.singleThickness = false;
    const zone::Patch crossed = zone::buildPatch(strakes, {0, -0.6, 0}, crossing);
    expectEqualCount("crossing the seam reaches the whole plate", crossed.elementCount(),
                    std::size_t{16});
    expectNear("at one thickness throughout", crossed.thickness, 0.012, 1e-15);
    expectTrue("and says which thickness it used", !crossed.problems.empty());

    // A fold: two flat plates meeting at 90 degrees, like a chine or a deck edge.
    StructuralMesh folded = flatStrip(1.6, 0.8, 0.012, 4, 2);
    for (PlatePanel& p : folded.panels)
        if (p.centroid().y > 0)
            for (Vec3& corner : p.corner)
                if (corner.y > 1e-9) {
                    corner.z = corner.y;
                    corner.y = 0.0;
                }
    const zone::Patch unfolded = zone::buildPatch(folded, {0, -0.2, 0}, flatParams(1));
    expectEqualCount("the zone does not turn a right-angle corner", unfolded.elementCount(),
                    std::size_t{4});
    expectTrue("and says so", !unfolded.problems.empty());

    // The radius is a bound, and it is the one thing that decides how big the zone
    // is. Nothing here should depend on the panel ordering, so a strike at the far
    // end must reach the same number of panels as one at the near end.
    const StructuralMesh wide = flatStrip(4.0, 0.8, 0.012, 10, 2);
    const zone::Patch near = zone::buildPatch(wide, {-1.8, 0, 0}, [] {
        zone::MeshParams p = flatParams(1);
        p.radius = 0.9;
        return p;
    }());
    const zone::Patch far = zone::buildPatch(wide, {1.8, 0, 0}, [] {
        zone::MeshParams p = flatParams(1);
        p.radius = 0.9;
        return p;
    }());
    expectEqualCount("a strike at either end reaches the same amount of plating",
                    near.elementCount(), far.elementCount());
    expectTrue("and the radius really bounded it", near.elementCount() < wide.panels.size());
    expectTrue("and they are different panels", near.panels != far.panels);
}

// --- 4. The meshed patch, against a closed form --------------------------------
//
// `solveStatic` and the element are already validated; what is under test is the
// mesh that `buildPatch` produced. A clamped square plate under uniform pressure
// has Timoshenko's tabulated central deflection 0.00126 q a^4 / D, which is the
// same number `test_solid_shell.cpp` checks its own generator against -- so this
// says the zone mesher builds the same plate that generator does, from panels.

void testAMeshedPatchReproducesAClampedPlate() {
    std::printf("\n   a meshed patch under pressure, against Timoshenko's clamped plate\n");
    const double side = 2.0, thickness = 0.02, pressure = 5.0e4;
    const StructuralMaterial steel = ah36Steel();
    const double rigidity = steel.youngsModulus * thickness * thickness * thickness /
                            (12.0 * (1.0 - steel.poissonRatio * steel.poissonRatio));
    const double want = 0.00126 * pressure * std::pow(side, 4) / rigidity;

    std::printf("     Timoshenko: %.6g m at the centre\n", want);
    double previous = 0;
    for (int subdivision : {2, 4, 8}) {
        const StructuralMesh plate = flatStrip(side, side, thickness, 2, 2);
        const zone::Patch patch = zone::buildPatch(plate, {0, 0, 0}, flatParams(subdivision));
        std::vector<double> displacement;
        const std::vector<double> load = solidshell::uniformPressureLoad(patch.mesh, pressure);
        std::string problem;
        const bool ok = solidshell::solveStatic(patch.mesh, patch.material,
                                                solidshell::Formulation::SolidShell, load,
                                                displacement, &problem);
        expectTrue("the meshed patch is solvable: " + problem, ok);
        if (!ok) return;

        // The centre node, found by position rather than by index, because the
        // mesher's numbering is its own business.
        double centre = 0;
        for (std::size_t node = 0; node < patch.nodeCount(); ++node) {
            const double x = patch.mesh.position[node * 3];
            const double y = patch.mesh.position[node * 3 + 1];
            if (std::abs(x) < 1e-9 && std::abs(y) < 1e-9)
                centre = std::min(centre, displacement[node * 3 + 2]);
        }
        std::printf("     %2d x %2d elements: %.6g m  (%+.2f%%)\n", 2 * subdivision,
                    2 * subdivision, -centre, 100.0 * (-centre / want - 1.0));
        if (previous != 0)
            expectTrue("refining the zone moves the answer towards the closed form",
                       std::abs(-centre - want) < std::abs(previous - want));
        previous = -centre;
    }
    expectNear("and reaches it", previous, want, 0.03 * want);

    // Guard against the whole thing being satisfied by a plate so thin the answer
    // is dominated by the pressure load rather than by the mesh: the deflection
    // has to be a real number, not machine noise, and it has to be a *plate*
    // deflection rather than a membrane one.
    expectTrue("the deflection is a real, small-deflection plate deflection",
               want > 1e-4 && want < 0.2 * thickness * 20.0);
}

// --- 5. Nothing happens to a patch nothing is happening to ---------------------

void testAPatchAtRestStaysAtRest() {
    std::printf("\n   a patch under no load\n");
    const StructuralMesh strip = flatStrip(1.6, 0.8, 0.02, 4, 2);

    for (zone::Edge edge : {zone::Edge::Clamped, zone::Edge::Free}) {
        zone::MeshParams params = flatParams(2);
        params.edge = edge;
        const zone::Patch patch = zone::buildPatch(strip, {0, 0, 0}, params);

        zone::SolveParams solve;
        solve.indenter.halfLength = 0.0;   // no punch at all
        solve.duration = 400.0 * patch.criticalTimestep;
        zone::Solver solver(patch, plasticity::shipSteel(), solve);
        const zone::SolveResult& result = solver.run();

        // Not `== 0`, and the reason is worth stating: the polar decomposition and
        // the strain-displacement product leave denormal dust at 1e-35 m and
        // 1e-55 J, twenty orders below anything physical. A bound at 1e-20 is an
        // identity in every sense that matters and is not hostage to the last bit
        // of a matrix inverse.
        expectTrue("nothing moves", solver.largestDisplacement() < 1e-20);
        expectTrue("no strain energy is stored", result.strainEnergy < 1e-20);
        expectTrue("nothing is dissipated", result.dissipation == 0.0);
        expectTrue("nothing is moving", result.kinetic < 1e-20);
        expectTrue("no work was done", result.work == 0.0);
        expectTrue("and it really took the steps", result.steps > 100);
        expectTrue("and nothing tore", result.tornPanels.empty());
    }

    // Non-vacuous: the same patch, the same number of steps, with the punch
    // switched on, does move. Without this the checks above pass on a solver that
    // does nothing at all.
    const zone::Patch patch = zone::buildPatch(strip, {0, 0, 0}, flatParams(2));
    zone::SolveParams solve;
    solve.indenter.halfLength = 0.06;
    solve.indenter.speed = 4.0;
    solve.indenter.stopAt = 1e9;
    solve.duration = 400.0 * patch.criticalTimestep;
    zone::Solver driven(patch, plasticity::shipSteel(), solve);
    driven.run();
    expectTrue("but a driven patch moves", driven.largestDisplacement() > 1e-6);
    expectTrue("and stores energy", driven.result().strainEnergy > 0);
}

// --- 6. Rigid body motion ------------------------------------------------------

void testRigidTranslationCarriesNoForce() {
    std::printf("\n   a rigid translation of the whole patch\n");
    const StructuralMesh strip = flatStrip(1.6, 0.8, 0.02, 4, 2);
    zone::MeshParams params = flatParams(2);
    params.edge = zone::Edge::Free;   // clamped, it would not be free to translate
    const zone::Patch patch = zone::buildPatch(strip, {0, 0, 0}, params);

    for (bool plastic : {false, true}) {
        zone::SolveParams solve;
        solve.plastic = plastic;
        solve.indenter.halfLength = 0.0;   // no punch: a held node is not free to travel
        solve.duration = 2000.0 * patch.criticalTimestep;
        zone::Solver solver(patch, plasticity::shipSteel(), solve);
        solver.translate({3.0, -2.0, 7.0});
        const double kinetic = solver.result().kinetic;
        solver.run();

        const double travel = solver.largestDisplacement();
        expectTrue(std::string(plastic ? "plastic" : "elastic") +
                       ": the patch actually moved a long way",
                   travel > 0.01);
        expectTrue(std::string(plastic ? "plastic" : "elastic") +
                       ": and stored no strain energy doing it",
                   solver.result().strainEnergy < 1e-9 * kinetic);
        expectTrue(std::string(plastic ? "plastic" : "elastic") + ": and dissipated nothing",
                   solver.result().dissipation == 0.0);
        expectNear(std::string(plastic ? "plastic" : "elastic") +
                       ": and kept every joule of its kinetic energy",
                   solver.result().kinetic, kinetic, 1e-9 * kinetic);
        // The energy that would have been stored had the patch been strained by
        // the same distance, as the scale the zero above is a zero against.
        expectTrue("the kinetic energy is a real number", kinetic > 1.0);
    }
}

// --- 7. The energy account -----------------------------------------------------

void testTheEnergyAccountCloses() {
    std::printf("\n   energy in = strain energy + plastic dissipation + kinetic + damping\n");
    const StructuralMesh strip = flatStrip(1.6, 0.8, 0.02, 8, 4);
    const zone::Patch patch = zone::buildPatch(strip, {0, 0, 0}, flatParams(1));

    const auto parameters = [](bool plastic, double depth, double damping) {
        zone::SolveParams solve;
        solve.plastic = plastic;
        solve.damping = damping;
        solve.indenter.halfLength = 0.06;
        solve.indenter.halfWidth = 1e3;
        solve.indenter.speed = 20.0;
        solve.indenter.rampTime = 1.0e-3;
        solve.indenter.stopAt = depth;
        return solve;
    };

    // The elastic path first, where the only approximation is the integrator.
    zone::Solver elasticSolver(patch, plasticity::shipSteel(),
                               parameters(false, 0.004, 1.0));
    const zone::SolveResult elastic = elasticSolver.run();
    std::printf("     elastic, 4 mm:            W %10.4g J, residual %+.3f%%\n", elastic.work,
                100.0 * elastic.energyResidual() / elastic.work);
    expectTrue("the elastic account closes", std::abs(elastic.energyResidual()) < 0.01 * elastic.work);
    expectTrue("and it had a real amount of energy to account for", elastic.work > 1.0e3);
    expectTrue("all of which is recoverable strain energy and motion",
               elastic.dissipation == 0.0 && elastic.strainEnergy > 0.5 * elastic.work);

    // The plastic path with a material that cannot yield has to reproduce it. This
    // is the tie between two whole solver paths -- different force routine,
    // different strain-energy expression, different state -- and it is the check
    // that would catch a Gauss weight or a compliance inverted the wrong way,
    // which no conservation test can see.
    zone::Solver rigidSolver(patch, rigidPlasticSteel(), parameters(true, 0.004, 1.0));
    const zone::SolveResult rigid = rigidSolver.run();

    // The tie is on the **trajectory**, not on an energy. A 0.2 ms run at 20 m/s
    // rings, so the instantaneous strain energy is an oscillating quantity and two
    // runs a hair out of phase disagree by 2% while being the same simulation --
    // which is what an energy comparison here reports, and it reads exactly like a
    // defect. The node positions do not oscillate against each other.
    double worst = 0, travelled = 0;
    for (std::size_t i = 0; i < elasticSolver.position().size(); ++i) {
        worst = std::max(worst, std::abs(elasticSolver.position()[i] - rigidSolver.position()[i]));
        travelled = std::max(travelled,
                             std::abs(elasticSolver.position()[i] - patch.mesh.position[i]));
    }
    std::printf("     the same through the plastic path: W %10.4g J, worst node %.2e m of"
                " %.2e m travelled\n", rigid.work, worst, travelled);
    expectTrue("a material that cannot yield follows the elastic path node for node",
               worst < 1e-9 * travelled);
    expectTrue("and it really went somewhere", travelled > 1e-3);
    expectNear("with the same work done on it", rigid.work, elastic.work, 1e-9 * elastic.work);
    expectTrue("through the plastic path, with nothing dissipated", rigid.dissipation == 0.0);
    expectTrue("and a real amount of energy stored", rigid.strainEnergy > 1.0e3);

    // And with real steel, yielding.
    zone::Solver yieldingSolver(patch, unbreakableSteel(), parameters(true, 0.010, 1.0));
    const zone::SolveResult yielding = yieldingSolver.run();
    std::printf("     plastic, 10 mm:           W %10.4g J, D/W %.3f, residual %+.3f%%\n",
                yielding.work, yielding.dissipation / yielding.work,
                100.0 * yielding.energyResidual() / yielding.work);
    expectTrue("it really is dissipating", yielding.dissipation > 0.3 * yielding.work);
    expectTrue("the plastic account closes",
               std::abs(yielding.energyResidual()) < 0.02 * yielding.work);

    // Damping is a sink like any other, and it has to be in the account or every
    // balance above would be satisfied by a solver that quietly removed energy.
    zone::Solver dampedSolver(patch, unbreakableSteel(), parameters(true, 0.010, 0.9995));
    const zone::SolveResult damped = dampedSolver.run();
    std::printf("     the same, damped:         W %10.4g J, damping took %.3f of it,"
                " residual %+.3f%%\n",
                damped.work, damped.dampingLoss / damped.work,
                100.0 * damped.energyResidual() / damped.work);
    expectTrue("damping removed a measurable share", damped.dampingLoss > 0.02 * damped.work);
    expectTrue("and the account still closes",
               std::abs(damped.energyResidual()) < 0.02 * damped.work);

    // **The limit, measured rather than assumed away.** The residual is not an
    // integrator error -- it does not move when the step is quartered -- it is the
    // co-rotational formulation's small-strain measure, and it grows with the
    // rotation the elements have to carry. Recorded here so that a future change
    // that made it worse would be visible.
    std::printf("     %10s %14s %14s\n", "depth", "element turn", "residual");
    double previousResidual = 0;
    for (double depth : {0.004, 0.02, 0.06}) {
        zone::Solver turning(patch, rigidPlasticSteel(), parameters(true, depth, 1.0));
        const zone::SolveResult run = turning.run();
        const double turn = std::atan2(depth, 0.4);
        std::printf("     %10.3f %13.3f  %+13.2f%%\n", depth, turn,
                    100.0 * run.energyResidual() / run.work);
        if (depth > 0.05)
            expectTrue("the account degrades with the rotation, and by this much",
                       std::abs(run.energyResidual()) > 3.0 * std::abs(previousResidual));
        previousResidual = run.energyResidual() / run.work * run.work;
        if (depth < 0.005) previousResidual = run.energyResidual();
    }

    // And it is not the timestep: quartering it leaves the residual where it was.
    zone::SolveParams fine;
    fine.plastic = true;
    fine.timestepSafety = 0.225;
    fine.indenter.halfLength = 0.06;
    fine.indenter.halfWidth = 1e3;
    fine.indenter.speed = 20.0;
    fine.indenter.rampTime = 1.0e-3;
    fine.indenter.stopAt = 0.010;
    zone::Solver fineSolver(patch, unbreakableSteel(), fine);
    const zone::SolveResult refined = fineSolver.run();
    expectNear("quartering the step does not move the residual, so it is not the integrator",
               refined.energyResidual() / refined.work, yielding.energyResidual() / yielding.work,
               0.004);
    expectTrue("and the refined run really took four times the steps",
               refined.steps > 3 * yielding.steps);
}

// --- 8. Against the membrane model ---------------------------------------------
//
// The load-bearing test. `indentation.hpp` is rigid-plastic membrane stretching
// of one bay: no bending, no hardening, no transverse constraint. Set up so that
// its own idealisation holds -- a long strip, a rigid line punch across it,
// boundaries held -- the two must agree to a factor, and the factor must be the
// one the two models' differences predict.

void testAgainstTheMembraneModel() {
    std::printf("\n   against indentation.hpp on the same bay\n");
    const double lengthX = 1.6, spanY = 0.8, thickness = 0.020, depth = 0.24;
    const StructuralMesh strip = flatStrip(lengthX, spanY, thickness, 8, 4);
    const zone::Patch patch = zone::buildPatch(strip, {0, 0, 0}, flatParams(1));

    zone::SolveParams solve;
    solve.indenter.halfLength = 0.06;  // a line across the strip, along `right` = -y
    solve.indenter.halfWidth = 1e3;    // the full length, along `up` = +x
    solve.indenter.speed = 20.0;
    solve.indenter.rampTime = 1.0e-3;
    solve.indenter.stopAt = depth;
    zone::Solver solver(patch, unbreakableSteel(), solve);
    const zone::SolveResult& fem = solver.run();

    IndentedPanel membrane;
    membrane.span = spanY;
    membrane.thickness = thickness;
    membrane.contactWidth = lengthX;
    membrane.yieldStrength = ah36Steel().yieldStrength;
    membrane.failureStrain = 1.0;
    expectTrue("the membrane model is being used inside what it claims",
               validateIndentation(membrane).empty());

    const double membraneForce = indentationForce(membrane, depth);
    const double membraneEnergy = indentationEnergy(membrane, depth);
    const double strain = membraneStrain(spanY, depth);

    std::printf("     at %.3f m: membrane strain %.4f\n", depth, strain);
    std::printf("     %-22s %14s %14s\n", "", "force (N)", "energy (J)");
    std::printf("     %-22s %14.4g %14.4g\n", "membrane", membraneForce, membraneEnergy);
    std::printf("     %-22s %14.4g %14.4g\n", "solid-shell FEM", fem.force, fem.work);
    std::printf("     %-22s %14.2f %14.2f\n", "ratio", fem.force / membraneForce,
                fem.work / membraneEnergy);

    expectTrue("the two models agree to a factor on force",
               fem.force > membraneForce && fem.force < 3.0 * membraneForce);
    expectTrue("and on energy", fem.work > membraneEnergy && fem.work < 3.0 * membraneEnergy);

    // **Why they differ, as an assertion rather than as a remark.** The membrane
    // model carries the tension at sigma_y flat and in one direction. The FEM
    // carries it at the *hardening* flow stress averaged over the path, and under
    // the plane-strain constraint the clamped side edges impose, which raises the
    // yield stress to 2/sqrt(3) of it. Both corrections are closed forms.
    const plasticity::FlowCurve& flow = unbreakableSteel().flow;
    double meanFlow = 0;
    const int samples = 2000;
    for (int i = 0; i < samples; ++i)
        meanFlow += plasticity::flowStress(flow, strain * (i + 0.5) / samples) / samples;
    const double predicted = membraneEnergy * (meanFlow / membrane.yieldStrength) * 2.0 /
                             std::sqrt(3.0);
    std::printf("     mean flow stress over the path %.4g Pa against sigma_y %.4g;"
                " plane strain 2/sqrt(3)\n", meanFlow, membrane.yieldStrength);
    std::printf("     corrected membrane energy %14.4g J, FEM / corrected %.3f\n", predicted,
                fem.work / predicted);
    // 25%, and what is left is bending and the strip's own ends: the membrane model
    // is a section of an infinite strip and this one is only twice as long as it
    // is wide, so material is drawn in from the ends as well as from the sides.
    expectNear("the disagreement is hardening and plane strain to within a quarter",
               fem.work / predicted, 1.0, 0.25);
    // Guard: the correction has to be doing real work, or the check above is the
    // uncorrected one wearing a hat.
    expectTrue("and the correction is a substantial one", predicted > 1.5 * membraneEnergy);

    // Quasi-static: whatever the punch speed, the answer is the same, so the
    // comparison against a rate-independent model is fair. Asserted through the
    // kinetic energy rather than by running it twice, which would double the cost
    // of the most expensive test in the file.
    std::printf("     kinetic energy is %.4f of the work, so the run is quasi-static\n",
                fem.kinetic / fem.work);
    expectTrue("the run is quasi-static", fem.kinetic < 0.01 * fem.work);
    expectTrue("and the plate is genuinely plastic, not ringing elastically",
               fem.dissipation > 0.9 * fem.strainEnergy * 5.0);

    std::printf("     cost: %d steps, %.3f s, %.2f us/element/step over %zu elements;"
                " %.0f core-s per simulated second\n",
                fem.steps, fem.wallSeconds, fem.microsecondsPerElementStep, patch.elementCount(),
                zone::estimatedCost(patch));
    // Two orders loose, per the rule test_plasticity.cpp wrote down: a timing
    // assertion tight enough to be interesting is a flaky test on a shared machine.
    expectTrue("the elastoplastic element costs microseconds, not nanoseconds or milliseconds",
               fem.microsecondsPerElementStep > 0.1 && fem.microsecondsPerElementStep < 500.0);
}

// --- 9. Torn panels, as indices into the structure `breach.hpp` takes -----------

void testTearingReportsPanelIndicesForBreach() {
    std::printf("\n   tearing: which panels, and in whose numbering\n");
    // Three small panels, four elements each, so a panel can be *partly* torn --
    // which is the only condition under which `tearFraction` means anything, and
    // the condition a one-element-per-panel mesh silently never reaches.
    const StructuralMesh strip = flatStrip(0.6, 0.3, 0.006, 3, 1);
    // Struck off-centre, so the punch straddles a panel seam and one panel is torn
    // through while its neighbour is torn part way. A symmetric strike tears whole
    // panels at a time and `tearFraction` never has anything to decide.
    const zone::Patch patch = zone::buildPatch(strip, {0.05, 0, 0}, flatParams(2));

    zone::SolveParams solve;
    solve.indenter.halfLength = 0.03;   // a line across the span
    solve.indenter.halfWidth = 0.08;    // over one seam, so the ends survive
    solve.indenter.speed = 20.0;
    solve.indenter.rampTime = 5.0e-4;
    solve.indenter.stopAt = 0.07;
    solve.historyStride = 100;
    zone::Solver solver(patch, plasticity::shipSteel(), solve);
    const zone::SolveResult& result = solver.run();

    std::printf("     %d of %zu elements deleted, %.4f m2; %zu panels torn of %zu meshed\n",
                result.tornElements, patch.elementCount(), result.tornArea,
                result.tornPanels.size(), patch.panels.size());
    expectTrue("elements tore", result.tornElements > 0);
    expectTrue("but not all of them, so the threshold has something to decide",
               result.tornElements < static_cast<int>(patch.elementCount()));
    expectTrue("panels are reported torn", !result.tornPanels.empty());
    // Where the first element let go, against the membrane model's tearing depth
    // on the same span with the same regularised failure strain. The two failure
    // criteria are the same one -- `plasticity::regularisedFailureStrain` -- fed a
    // strain the two models compute completely differently, so agreeing to a
    // factor is a statement about the *kinematics*.
    double firstTear = 0;
    for (const zone::Sample& sample : result.history)
        if (sample.tornElements > 0 && firstTear == 0) firstTear = sample.penetration;
    const double failureStrain =
        plasticity::regularisedFailureStrain(plasticity::shipSteel().failure, 0.15, 0.006);
    const double membraneDepth = penetrationForStrain(0.3, failureStrain);
    std::printf("     first element let go at %.4f m; the membrane model tears at %.4f m\n",
                firstTear, membraneDepth);
    expectTrue("something tore during the recorded history", firstTear > 0);
    // It tears *earlier*, and the direction is predictable rather than lucky: the
    // membrane model spreads one strain over the whole leg, where the FEM has the
    // punch edge and the clamped support concentrating it.
    expectTrue("and it tore earlier than the membrane model, as bending says it must",
               firstTear < membraneDepth);
    expectTrue("but within a factor of it", firstTear > 0.4 * membraneDepth);

    // The indices are into the *original* `StructuralMesh`, ascending and unique,
    // which is exactly the contract `breachesFromFailedPanels` is written against.
    for (std::size_t i = 0; i < result.tornPanels.size(); ++i) {
        expectTrue("every torn panel is a real panel of the structure",
                   result.tornPanels[i] >= 0 &&
                       result.tornPanels[i] < static_cast<int>(strip.panels.size()));
        if (i > 0)
            expectTrue("and they come back ascending and distinct",
                       result.tornPanels[i] > result.tornPanels[i - 1]);
    }
    for (int index : result.tornPanels) {
        bool meshed = false;
        for (int covered : patch.panels) meshed = meshed || covered == index;
        expectTrue("and only panels the zone actually meshed can tear", meshed);
    }

    // The threshold is a decision, so it has to change the answer in the direction
    // it claims. Demanding the whole panel reports fewer than demanding a quarter.
    const std::size_t strict = solver.tornPanelsAt(1.0).size();
    const std::size_t half = solver.tornPanelsAt(0.5).size();
    const std::size_t loose = solver.tornPanelsAt(0.2).size();
    std::printf("     tearFraction 0.2 / 0.5 / 1.0 -> %zu / %zu / %zu panels\n", loose, half,
                strict);
    expectEqualCount("the default threshold is 0.5 and that is what was reported",
                     result.tornPanels.size(), half);
    expectTrue("a stricter threshold reports no more panels", strict <= half);
    expectTrue("a looser one reports no fewer", loose >= half);
    expectTrue("and the threshold is not decorative", loose > strict);

    // The other end of the chain: `breach.hpp` takes them unchanged. It is run
    // against the reference ferry, whose structure is the one a real zone would be
    // cut from, so the indices have to survive the round trip through a mesh with
    // nine thousand panels in it.
    const Ship ferry = game::buildFerry();
    const StructuralMesh structure = makeStructuralMesh(ferry.hull, ferryScantlings());
    zone::MeshParams params;
    params.radius = 2.0;
    params.subdivision = 1;
    const zone::Patch side = zone::buildPatch(structure, {0.0, -8.5, 8.0}, params);
    expectTrue("a zone on the ferry's side meshes something", !side.empty());
    const BreachSet breaches = breachesFromFailedPanels(ferry, structure, side.panels);
    std::printf("     the ferry's own panels, fed to breach.hpp: %zu openings, %.2f m2\n",
                breaches.breaches.size(), breaches.totalArea());
    expectTrue("and breach.hpp opens holes from them", !breaches.breaches.empty());
    expectTrue("whose area is the plating's", breaches.totalArea() > 0.5 * side.area);
}

// --- 10. The parallel element loop is the serial one ---------------------------

void testTheParallelSolveIsBitIdentical() {
    std::printf("\n   the threaded element loop, against the serial one\n");
    const StructuralMesh strip = flatStrip(2.0, 1.0, 0.02, 8, 4);
    const zone::Patch patch = zone::buildPatch(strip, {0, 0, 0}, flatParams(2));
    expectTrue("the patch is big enough for the solver to bother dispatching",
               patch.elementCount() >= 64);

    const auto drive = [&](core::JobSystem* jobs) {
        zone::SolveParams solve;
        solve.jobs = jobs;
        solve.indenter.halfLength = 0.1;
        solve.indenter.halfWidth = 1e3;
        solve.indenter.speed = 20.0;
        solve.indenter.rampTime = 5.0e-4;
        solve.indenter.stopAt = 0.008;
        zone::Solver solver(patch, unbreakableSteel(), solve);
        solver.run();
        return solver.position();
    };

    const std::vector<double> serial = drive(nullptr);
    core::JobSystem jobs(4);
    const std::vector<double> threaded = drive(&jobs);

    expectEqualCount("the two runs produced the same mesh", threaded.size(), serial.size());
    bool identical = true;
    double moved = 0;
    for (std::size_t i = 0; i < serial.size(); ++i) {
        identical = identical && serial[i] == threaded[i];
        moved = std::max(moved, std::abs(serial[i] - patch.mesh.position[i]));
    }
    // Bit-identical, not near: the nodal gather is by CSR in a fixed order, which
    // is the whole reason `fem.cpp` gathers rather than scatters. "Near" would
    // pass on a solver that had quietly become worker-count dependent, which is
    // the property multiplayer and replays rest on.
    expectTrue("four workers give bit-identical positions", identical);
    expectTrue("and the run actually deformed the patch", moved > 1e-4);
}

// --- 11. The reference ferry ---------------------------------------------------

void testAZoneOnTheFerry() {
    std::printf("\n   a zone cut from the reference ferry\n");
    const Ship ferry = game::buildFerry();
    const StructuralMesh structure = makeStructuralMesh(ferry.hull, ferryScantlings());

    std::printf("     %-26s %7s %8s %8s %10s %9s %12s\n", "impact", "panels", "elems", "area",
                "offset/t", "excess", "core-s/s");
    struct Case {
        const char* name;
        Vec3 point;
        double radius;
    };
    for (const Case& probe : {Case{"side, flat of side", {0.0, -8.5, 8.0}, 4.0},
                              Case{"side, over the shoulder", {0.0, -8.5, 5.0}, 4.0},
                              Case{"bottom", {0.0, -2.0, 0.05}, 4.0}}) {
        zone::MeshParams params;
        params.radius = probe.radius;
        params.subdivision = 2;
        params.stiffeners = zone::Stiffeners::Ignored;
        const zone::Patch patch = zone::buildPatch(structure, probe.point, params);
        std::printf("     %-26s %7zu %8zu %8.1f %10.4f %8.0f%% %12.0f\n", probe.name,
                    patch.panels.size(), patch.elementCount(), patch.area,
                    patch.worstNormalSpread, 100.0 * patch.spuriousStiffness,
                    zone::estimatedCost(patch));

        expectTrue(std::string(probe.name) + ": the zone meshed real plating", !patch.empty());
        expectTrue(std::string(probe.name) + ": no element is inverted", patch.worstJacobian > 0);
        expectNear(std::string(probe.name) + ": the meshed area is the panels' own area",
                   patch.area, [&] {
                       double total = 0;
                       for (int index : patch.panels)
                           total += structure.panels[static_cast<std::size_t>(index)].area();
                       return total;
                   }(), 1e-9 * patch.area);
        expectTrue(std::string(probe.name) + ": every panel carries the patch's thickness", [&] {
            for (int index : patch.panels)
                if (std::abs(structure.panels[static_cast<std::size_t>(index)].thickness -
                             patch.thickness) > 1e-12)
                    return false;
            return true;
        }());
        // The stable step is thickness governed, so it has to be the plating's own
        // dilatational transit time however big the elements are in plane.
        const double waveSpeed = std::sqrt(patch.material.youngsModulus *
                                           (1.0 - patch.material.poissonRatio) /
                                           (patch.material.density *
                                            (1.0 + patch.material.poissonRatio) *
                                            (1.0 - 2.0 * patch.material.poissonRatio)));
        expectNear(std::string(probe.name) + ": the stable step is t / c_p",
                   patch.criticalTimestep, 0.9 * patch.thickness / waveSpeed,
                   0.1 * patch.criticalTimestep);
    }

    // The zone is *connected*, so a strike on the starboard side cannot reach the
    // port side however narrow the ship is there. A radius test alone would.
    zone::MeshParams params;
    params.radius = 12.0;
    params.subdivision = 1;
    params.stiffeners = zone::Stiffeners::Ignored;
    const zone::Patch wide = zone::buildPatch(structure, {0.0, -9.9, 10.0}, params);
    int wrongSide = 0;
    for (int index : wide.panels)
        if (structure.panels[static_cast<std::size_t>(index)].centroid().y > 0) ++wrongSide;
    expectTrue("a twelve-metre zone on the starboard side reaches a long way",
               wide.area > 40.0);
    expectEqual("and none of it is on the port side", wrongSide, 0LL);

    // The flat of side is flat; the shoulder is not, and the difference is the
    // whole content of `07-fem-spike-findings.md` §6 limit 1 on a real ship.
    zone::MeshParams tight;
    tight.radius = 2.0;
    tight.subdivision = 2;
    tight.stiffeners = zone::Stiffeners::Ignored;
    const zone::Patch flat = zone::buildPatch(structure, {0.0, -10.0, 9.0}, tight);
    const zone::Patch shoulder = zone::buildPatch(structure, {0.0, -9.8, 4.2}, tight);
    std::printf("     flat of side offset/t %.4f (+%.0f%%); over the shoulder %.4f (+%.0f%%)\n",
                flat.worstNormalSpread, 100.0 * flat.spuriousStiffness,
                shoulder.worstNormalSpread, 100.0 * shoulder.spuriousStiffness);
    expectTrue("a zone on the flat of side is prismatic to within a degree",
               flat.worstNormalSpread < 0.02);
    expectTrue("and raises nothing", flat.problems.empty());
    expectTrue("a zone over the shoulder is not, by a factor of five or more",
               shoulder.worstNormalSpread > 5.0 * std::max(flat.worstNormalSpread, 1e-4));
    expectTrue("and says so", !shoulder.problems.empty());

    // The stiffener lines are the panel seams, so a subdivision that leaves nothing
    // between them leaves nothing free to bend -- and the patch has to say so
    // rather than solve a zone that cannot move.
    for (int subdivision : {1, 2, 4}) {
        zone::MeshParams params;
        params.radius = 3.0;
        params.subdivision = subdivision;
        const zone::Patch held = zone::buildPatch(structure, {0.0, -9.9, 8.0}, params);
        std::printf("     subdivision %d with rigid stiffeners: %d node(s) held,"
                    " %.0f%% of the DOF free\n", subdivision, held.stiffenerNodes,
                    100.0 * held.freeFraction);
        expectTrue("the ferry's longitudinals and frames are found and held",
                   held.stiffenerNodes > 0);
        if (subdivision < 3)
            expectTrue("and a coarse zone is reported as having nothing free to deform",
                       held.freeFraction < 0.25 && !held.problems.empty());
        else
            expectTrue("while a fine one has plating between the supports",
                       held.freeFraction > 0.4);
    }
    zone::MeshParams loose;
    loose.radius = 3.0;
    loose.subdivision = 4;
    loose.stiffeners = zone::Stiffeners::Ignored;
    const zone::Patch ignored = zone::buildPatch(structure, {0.0, -9.9, 8.0}, loose);
    expectEqual("ignoring the stiffeners holds no node at all", ignored.stiffenerNodes, 0LL);
}

// --- 12. Stiffeners: the bracket, and the span it sets --------------------------
//
// The zone meshes plating and not stiffeners (`zone.hpp` §3), so the two things it
// can honestly do with one are nothing and hold it rigidly. Those bracket the
// answer, and the bracket has to be *wide* -- if it were narrow the omission would
// not matter and the header's argument would be wrong.

void testStiffenersBracketTheAnswer() {
    std::printf("\n   stiffeners: the two bounds, and the span they set\n");
    const double lengthX = 1.6, spanY = 0.8, thickness = 0.020, depth = 0.05;
    const StructuralMesh strip = flatStrip(lengthX, spanY, thickness, 4, 4, true);
    expectEqualCount("the strip carries two longitudinals", strip.members.size(),
                     std::size_t{2});

    const auto run = [&](zone::Stiffeners stiffeners) {
        zone::MeshParams params = flatParams(2);
        params.stiffeners = stiffeners;
        const zone::Patch patch = zone::buildPatch(strip, {0, 0, 0}, params);
        zone::SolveParams solve;
        solve.indenter.halfLength = 0.06;
        solve.indenter.halfWidth = 1e3;
        solve.indenter.speed = 20.0;
        solve.indenter.rampTime = 1.0e-3;
        solve.indenter.stopAt = depth;
        zone::Solver solver(patch, unbreakableSteel(), solve);
        const zone::SolveResult result = solver.run();
        return std::pair<int, zone::SolveResult>{patch.stiffenerNodes, result};
    };

    const std::pair<int, zone::SolveResult> loose = run(zone::Stiffeners::Ignored);
    const std::pair<int, zone::SolveResult> held = run(zone::Stiffeners::RigidSupport);
    expectEqual("ignoring the stiffener holds no node",
                static_cast<long long>(loose.first), 0LL);
    expectTrue("holding it holds the nodes along its line", held.first > 0);

    // With the stiffener held the plating spans half as far, so at the same depth
    // the membrane strain is about four times as large and the resistance rises
    // with it. Each reading is compared against `indentation.hpp` at *its own*
    // span, which is the point: the model's span is the distance between supports,
    // and a zone with no supports in it has only its own perimeter.
    const auto membrane = [&](double span) {
        IndentedPanel model;
        model.span = span;
        model.thickness = thickness;
        model.contactWidth = lengthX;
        model.yieldStrength = ah36Steel().yieldStrength;
        model.failureStrain = 1.0;
        return indentationForce(model, depth);
    };
    std::printf("     %-28s %14s %14s %8s\n", "", "force (N)", "membrane (N)", "ratio");
    std::printf("     %-28s %14.4g %14.4g %8.2f\n", "stiffener ignored", loose.second.force,
                membrane(spanY), loose.second.force / membrane(spanY));
    std::printf("     %-28s %14.4g %14.4g %8.2f\n", "stiffener held rigid", held.second.force,
                membrane(0.5 * spanY), held.second.force / membrane(0.5 * spanY));

    const double bracket = held.second.force / loose.second.force;
    const double predicted = membrane(0.5 * spanY) / membrane(spanY);
    std::printf("     the bracket is %.2f wide on force and %.2f on energy; halving the span"
                " predicts %.2f\n", bracket, held.second.work / loose.second.work, predicted);
    expectTrue("holding the stiffener stiffens the plating substantially", bracket > 1.4);
    // And by the amount halving the span says it should, which is what makes this
    // a statement about the support rather than about the solver.
    expectNear("by the factor halving the membrane span predicts", bracket, predicted,
               0.3 * predicted);
    // Each end of the bracket agrees with the membrane model on *its own* span, to
    // the hardening-and-plane-strain factor the previous test pinned at about two.
    expectTrue("ignored agrees with a membrane spanning the whole zone",
               loose.second.force > 0.7 * membrane(spanY) &&
                   loose.second.force < 3.0 * membrane(spanY));
    expectTrue("held agrees with a membrane spanning between supports",
               held.second.force > 0.7 * membrane(0.5 * spanY) &&
                   held.second.force < 3.0 * membrane(0.5 * spanY));
    // And the two readings really are different, or the checks above are satisfied
    // by a model that does not care what its span is.
    expectTrue("the two spans are not the same question",
               membrane(0.5 * spanY) > 1.5 * membrane(spanY));
}

// --- 13. What the weld tolerance means, and what the bucket grid is for ---------
//
// Mutation testing found this whole class untested: squaring the tolerance,
// shrinking the bucket, and dropping the neighbour probe all passed, because every
// mesh in the file has its duplicate points either bit-identical or centimetres
// apart. The tolerance is a *distance*, and the grid has to find a pair that
// straddles a cell boundary.

void testTheWeldIsADistanceAndNotABucket() {
    std::printf("\n   welding: the tolerance is a distance, and the bucket grid overlaps\n");
    const double thickness = 0.012;

    // Two rows of panels with a gap between them. Well inside the tolerance they
    // are one patch; well outside it they are two, and the zone reaches only the
    // struck one. Both directions, because a tolerance that is too loose and one
    // that is too tight both pass a one-sided check.
    const auto splitAt = [&](double gap) {
        StructuralMesh mesh = flatStrip(0.8, 0.4, thickness, 2, 1);
        for (PlatePanel& p : mesh.panels)
            if (p.centroid().x > 0)
                for (Vec3& corner : p.corner) corner.x += gap;
        return zone::buildPatch(mesh, {-0.2, 0, 0}, flatParams(1)).elementCount();
    };
    // The tolerance is 1 um. Half of it welds; one and a half times it does not,
    // and that is the band the check has to straddle -- a gap of 0.1 mm is refused
    // by any comparison anyone might write, including one that has forgotten to
    // square the tolerance and is really working at a millimetre.
    expectEqualCount("a seam half a tolerance wide is one weld, so both panels are in the zone",
                     splitAt(5.0e-7), std::size_t{2});
    expectEqualCount("a seam one and a half tolerances wide is a crack, and the zone stops",
                     splitAt(1.5e-6), std::size_t{1});
    expectEqualCount("and so is one at a hundred tolerances", splitAt(1.0e-4), std::size_t{1});

    // A shared corner reached by two different arithmetic routes, landing one ulp
    // apart *at a coordinate where the weld's bucket grid splits*. Only a probe
    // that looks in the neighbouring cell finds the pair, and without one the mesh
    // comes apart down that seam -- silently, because everything else about it is
    // right. The mirrored starboard panels are the real case: they reach their
    // shared corners by a different expression from the port ones.
    const double cell = 2.0 * 1e-6;   // the weld's cell is twice its tolerance
    double above = 0.8, below = std::nextafter(above, 0.0);
    for (int step = 0; step < 64 && std::floor(above / cell) == std::floor(below / cell); ++step) {
        above = below;
        below = std::nextafter(above, 0.0);
    }
    expectTrue("the two coordinates are one ulp apart", above - below < 1e-15 && above > below);
    expectTrue("and a bucket boundary runs between them",
               std::floor(above / cell) != std::floor(below / cell));

    // Built directly, so the seam really sits on those two coordinates: the after
    // panel ends at `above` and the forward one begins at `below`.
    StructuralMesh straddle;
    straddle.materials.push_back(ah36Steel());
    const double ends[2][2] = {{above - 0.8, above}, {below, below + 0.8}};
    for (const auto& span : ends) {
        PlatePanel p;
        p.corner[0] = {span[0], -0.2, 0};
        p.corner[1] = {span[1], -0.2, 0};
        p.corner[2] = {span[1], 0.2, 0};
        p.corner[3] = {span[0], 0.2, 0};
        p.thickness = thickness;
        p.material = 0;
        p.role = PanelRole::Shell;
        straddle.panels.push_back(p);
    }
    const zone::Patch welded = zone::buildPatch(straddle, {above - 0.4, 0, 0}, flatParams(1));
    expectEqualCount("a corner one ulp across a bucket boundary still welds",
                     welded.elementCount(), std::size_t{2});
    expectEqualCount("into one grid of nodes", welded.nodeCount(), std::size_t{12});
}

// --- 14. Which way is out ------------------------------------------------------

void testTheAxisPointsOutward() {
    std::printf("\n   orientation: which way is out, and which face the punch touches\n");
    const double thickness = 0.02;
    const StructuralMesh forward = flatStrip(1.6, 0.8, thickness, 4, 2);
    StructuralMesh reversed = forward;
    for (PlatePanel& p : reversed.panels) std::swap(p.corner[1], p.corner[3]);
    expectTrue("the two meshes really are wound opposite ways",
               dot(forward.panels[0].normal(), reversed.panels[0].normal()) < -0.99);

    const zone::Patch a = zone::buildPatch(forward, {0, 0, 0}, flatParams(2));
    const zone::Patch b = zone::buildPatch(reversed, {0, 0, 0}, flatParams(2));
    expectTrue("the requested outward direction is the patch's axis, whichever way the"
               " panels are wound",
               dot(a.axis, Vec3{0, 0, 1}) > 0.99 && dot(b.axis, Vec3{0, 0, 1}) > 0.99);
    // The node grid is numbered by position, so the two patches must come out with
    // *identical* coordinates -- the strongest statement available that rewinding a
    // panel changed nothing but the corner order inside an element.
    expectEqualCount("and the same number of nodes", b.nodeCount(), a.nodeCount());
    double worst = 0;
    for (std::size_t i = 0; i < a.mesh.position.size() && i < b.mesh.position.size(); ++i)
        worst = std::max(worst, std::abs(a.mesh.position[i] - b.mesh.position[i]));
    // Not bit-identical: the bilinear grid reaches the same point by a different
    // route when the corners are relabelled, so the two agree to rounding on a
    // 1.6 m strip and not to the last bit. Twelve orders below the element size is
    // an identity in every sense that matters here.
    expectTrue("in identical positions", worst < 1e-12);
    expectTrue("with the same face called the outer one", a.outerFace == b.outerFace);
    expectTrue("and no element inverted either way", a.worstJacobian > 0 && b.worstJacobian > 0);
    // The outer face is the one on the +axis side, which is what the punch touches.
    for (std::size_t node = 0; node < a.nodeCount(); ++node)
        if (a.outerFace[node])
            expectTrue("every outer-face node is on the +axis side",
                       a.mesh.position[node * 3 + 2] > 0);

    // With no `outward` given the direction comes from the structure's own centroid,
    // which is the path a real caller takes. On the ferry it has to point out of her
    // on both sides and downward on the bottom.
    const StructuralMesh structure = makeStructuralMesh(game::buildFerry().hull,
                                                        ferryScantlings());
    zone::MeshParams derived;
    derived.radius = 2.0;
    derived.subdivision = 1;
    derived.stiffeners = zone::Stiffeners::Ignored;
    const zone::Patch starboard = zone::buildPatch(structure, {0.0, -9.9, 8.0}, derived);
    const zone::Patch port = zone::buildPatch(structure, {0.0, 9.9, 8.0}, derived);
    const zone::Patch bottom = zone::buildPatch(structure, {0.0, -2.0, 0.05}, derived);
    std::printf("     ferry axes: starboard %+.2f, port %+.2f (y);  bottom %+.2f (z)\n",
                starboard.axis.y, port.axis.y, bottom.axis.z);
    expectTrue("a strike on her starboard side looks outboard to starboard",
               starboard.axis.y < -0.9);
    expectTrue("and one on her port side outboard to port", port.axis.y > 0.9);
    expectTrue("and one on her bottom looks down", bottom.axis.z < -0.9);
}

// --- 15. Areas, on the one quad where the shortcut is wrong ---------------------

void testAreasAreRightOnATrapezoid() {
    std::printf("\n   area: a trapezoid, where halving one diagonal is not enough\n");
    const double lengthX = 2.0, spanY = 1.0, taper = 0.4, thickness = 0.012;
    const StructuralMesh mesh = trapezoid(lengthX, spanY, taper, thickness);
    // Analytic: a trapezoid of height lengthX with parallel sides spanY and
    // taper*spanY.
    const double want = 0.5 * lengthX * spanY * (1.0 + taper);

    for (int subdivision : {1, 2, 4}) {
        const zone::Patch patch = zone::buildPatch(mesh, {0, 0, 0}, flatParams(subdivision));
        expectNear("the meshed area is the trapezoid's, at every subdivision", patch.area, want,
                   1e-12 * want);
        double total = 0;
        for (double area : patch.elementArea) total += area;
        expectNear("and the elements' own areas sum to it", total, want, 1e-12 * want);
    }
    // The guard: on a rectangle the wrong construction would give the same answer,
    // so the taper has to be doing work.
    expectTrue("the trapezoid is genuinely not a parallelogram", taper < 0.9);
    std::printf("     %.6f m2 against %.6f m2 analytic\n",
                zone::buildPatch(mesh, {0, 0, 0}, flatParams(2)).area, want);
}

// --- 16. Nodal normals are area weighted ---------------------------------------

void testNodalNormalsAreAreaWeighted() {
    std::printf("\n   nodal normals: area weighted, which only a graded mesh can see\n");
    const double angle = 0.30, wide = 0.8, narrow = 0.2, thickness = 0.02;
    const StructuralMesh mesh = gradedFold(1.0, wide, narrow, angle, thickness);
    zone::MeshParams params = flatParams(1);
    params.normalSpreadWarning = 1e-3;
    const zone::Patch patch = zone::buildPatch(mesh, {0, -0.4, 0}, params);
    expectEqualCount("the fold meshed two elements", patch.elementCount(), std::size_t{2});

    // The seam node's normal is the normalised sum of the two faces' area vectors.
    // Weighted, it leans towards the wide panel by wide:narrow; unweighted it sits
    // half way. Each element then has two corners on its own face normal and two on
    // the seam, and the patch reports the worse of the two -- which is the *narrow*
    // panel's, because the seam normal has been pulled away from it.
    const auto spreadFor = [&](double weightWide, double weightNarrow) {
        const Vec3 flat{0, 0, 1};
        const Vec3 tilted{0, -std::sin(angle), std::cos(angle)};
        const Vec3 seam = normalize(flat * weightWide + tilted * weightNarrow);
        const Vec3 wideMean = normalize(flat * 2.0 + seam * 2.0);
        const Vec3 narrowMean = normalize(seam * 2.0 + tilted * 2.0);
        return std::max(length(flat - wideMean), length(seam - narrowMean));
    };
    const double weighted = spreadFor(wide, narrow);
    const double even = spreadFor(1.0, 1.0);
    std::printf("     area weighted predicts %.5f, an even average %.5f; measured %.5f\n",
                weighted, even, patch.worstNormalSpread);
    expectTrue("the two predictions differ, so the weighting is testable at all",
               std::abs(weighted - even) > 0.2 * even);
    expectNear("the measured offset is the area-weighted one", patch.worstNormalSpread, weighted,
               0.02 * weighted);
}

// --- 17. What `estimatedCost` is a function of ---------------------------------

void testTheCostEstimateIsAFunctionOfTheRightThings() {
    std::printf("\n   the cost estimate: elements, the timestep, and the constitutive path\n");
    const StructuralMesh thin = flatStrip(1.6, 0.8, 0.010, 4, 2);
    const StructuralMesh thick = flatStrip(1.6, 0.8, 0.020, 4, 2);
    const zone::Patch coarse = zone::buildPatch(thin, {0, 0, 0}, flatParams(1));
    const zone::Patch fine = zone::buildPatch(thin, {0, 0, 0}, flatParams(2));
    const zone::Patch stout = zone::buildPatch(thick, {0, 0, 0}, flatParams(1));

    std::printf("     %2zu elements at 10 mm: %8.1f core-s/s;  %2zu at 10 mm: %8.1f;"
                "  %2zu at 20 mm: %8.1f\n",
                coarse.elementCount(), zone::estimatedCost(coarse), fine.elementCount(),
                zone::estimatedCost(fine), stout.elementCount(), zone::estimatedCost(stout));

    // Four times the elements, four times the cost: the step does not care about the
    // in-plane size, which is the whole reason the zone is bounded by element count
    // rather than by area.
    expectNear("four times the elements is four times the cost", zone::estimatedCost(fine),
               4.0 * zone::estimatedCost(coarse), 0.02 * zone::estimatedCost(fine));
    // Twice the plate is twice the step, so half the cost. Without the step count in
    // it the estimate would not move at all.
    expectNear("twice the plate thickness is half the cost", zone::estimatedCost(stout),
               0.5 * zone::estimatedCost(coarse), 0.02 * zone::estimatedCost(coarse));
    expectTrue("and the two thicknesses really do differ in step",
               stout.criticalTimestep > 1.8 * coarse.criticalTimestep);
    // And the elastic path is cheaper by the ratio the measurements say. The
    // plastic figure moved from 7.3 µs to 3.1 when the step-invariant element forms
    // stopped being rebuilt every step, so this ratio moved from 27 to 11 -- the
    // elastic path barely uses those forms, so only one of the two numbers changed.
    expectNear("the elastic path is 11 times cheaper, as measured",
               zone::estimatedCost(coarse, true) / zone::estimatedCost(coarse, false),
               3.1 / 0.273, 0.1);
    expectTrue("an empty patch costs nothing", zone::estimatedCost(zone::Patch{}) == 0.0);
}

// --- 18. A stiffener holds the plating it reaches, and no more ------------------

void testStiffenersHoldOnlyWhatTheyReach() {
    std::printf("\n   stiffeners: a member holds the nodes on it, not the ones beyond it\n");
    const double lengthX = 1.6, spanY = 0.8, thickness = 0.02;
    StructuralMesh mesh = flatStrip(lengthX, spanY, thickness, 4, 4, true);
    // Cut both longitudinals back to the after half of the strip. The nodes on the
    // forward half are on the member's *line* but past its end, and a projection
    // that is not clamped to the segment holds them anyway.
    for (StructuralMember& member : mesh.members) member.b.x = 0.0;

    zone::MeshParams params = flatParams(2);
    params.radius = 2.0;   // a real radius, so the members' own proximity filter has work to do
    const zone::Patch half = zone::buildPatch(mesh, {0, 0, 0}, params);
    StructuralMesh whole = flatStrip(lengthX, spanY, thickness, 4, 4, true);
    const zone::Patch full = zone::buildPatch(whole, {0, 0, 0}, params);
    std::printf("     members over the whole length hold %d nodes; over half, %d\n",
                full.stiffenerNodes, half.stiffenerNodes);
    expectTrue("a member reaching half way holds about half as many nodes",
               half.stiffenerNodes > 0 && half.stiffenerNodes < full.stiffenerNodes);
    // Nine node columns across the strip at subdivision 2, five of them on the after
    // half including the one on the cut, times two seams.
    expectEqual("exactly the nodes on the segment", static_cast<long long>(half.stiffenerNodes),
                10LL);
    expectEqual("and all of them when it runs the whole way",
                static_cast<long long>(full.stiffenerNodes), 18LL);

    // A member far outside the zone changes nothing, which is what the proximity
    // filter is for -- and it has to be measured from the *segment*, because a
    // girder running the length of a ship has its midpoint nowhere near the zone it
    // passes through.
    StructuralMesh distant = whole;
    StructuralMember runner = distant.members[0];
    // Its *midpoint* is two hundred metres away and it runs straight through the
    // zone, which is what a centre girder does to a zone amidships. A proximity
    // filter measured from the midpoint drops it.
    runner.a = {-0.9, 0.0, 0.0};
    runner.b = {400.0, 0.0, 0.0};
    distant.members.push_back(runner);
    const zone::Patch reached = zone::buildPatch(distant, {0, 0, 0}, params);
    std::printf("     a girder whose midpoint is 200 m away but which runs through the zone"
                " takes the count from %d to %d\n", full.stiffenerNodes, reached.stiffenerNodes);
    expectEqual("a member whose midpoint is far away but which runs through the zone is"
                " still found",
                static_cast<long long>(reached.stiffenerNodes),
                static_cast<long long>(full.stiffenerNodes) + 9LL);
}


// The rest forms are cached for the life of a solve, which is a 2x cost decision
// and must be no other kind of decision at all. Everything the solver reports has
// to come back bit for bit the same with the cache off.
//
// Two guards against a vacuous pass. The patch has to have **deformed and torn**,
// because an undisturbed patch agrees with anything; and the run has to be long
// enough that the state has moved many times, because a cache that is stale rather
// than wrong agrees on the first step and diverges later.
void testTheRestFormsCacheChangesNothing() {
    std::printf("\n   the cached rest forms, against rebuilding them every step\n");
    // The same small, thin strip `testTearingReportsPanelIndicesForBreach` uses,
    // driven past first tear -- so the comparison covers the element-deletion path,
    // the dropped enhanced modes and the in-step retry, which are the parts of the
    // element that a merely-deforming patch never reaches.
    const StructuralMesh strip = flatStrip(0.6, 0.3, 0.006, 3, 1);
    const zone::Patch patch = zone::buildPatch(strip, {0.05, 0, 0}, flatParams(2));

    struct Answer {
        std::vector<double> position;
        zone::SolveResult result;
    };
    const auto drive = [&](bool cache, bool preload = false) {
        zone::SolveParams solve;
        solve.cacheRestForms = cache;
        if (preload) {
            // A real pre-strain: 120 MPa of the ship's own hogging stress through
            // the patch, which is a third of yield and moves `rest_` away from
            // `position_` by ~6e-4 of every coordinate.
            solve.preload.stress = 120.0e6;
            solve.preload.gradient = 8.0e6;
            solve.preload.reference = 0.0;
        }
        solve.indenter.halfLength = 0.03;
        solve.indenter.halfWidth = 0.08;
        solve.indenter.speed = 20.0;
        solve.indenter.rampTime = 5.0e-4;
        solve.indenter.stopAt = 0.07;
        zone::Solver solver(patch, plasticity::shipSteel(), solve);
        Answer answer;
        answer.result = solver.run();
        answer.position = solver.position();
        return answer;
    };

    const Answer cached = drive(true);
    const Answer rebuilt = drive(false);
    // And the same again **with a pre-load**, which is the only condition under
    // which the rest configuration and the initial position differ. Without one
    // they are the same array, and building the cached forms from the wrong one of
    // the two is then invisible -- mutation testing found exactly that edit
    // surviving the whole suite, and this is what closes it.
    const Answer preloadedCache = drive(true, true);
    const Answer preloadedRebuild = drive(false, true);
    expectTrue("the cache reports itself on", cached.result.cachedRestForms);
    expectTrue("and off", !rebuilt.result.cachedRestForms);
    expectEqualCount("both runs took the same number of steps",
                     static_cast<std::size_t>(cached.result.steps),
                     static_cast<std::size_t>(rebuilt.result.steps));

    bool identical = cached.result.force == rebuilt.result.force &&
                     cached.result.peakForce == rebuilt.result.peakForce &&
                     cached.result.work == rebuilt.result.work &&
                     cached.result.strainEnergy == rebuilt.result.strainEnergy &&
                     cached.result.dissipation == rebuilt.result.dissipation &&
                     cached.result.kinetic == rebuilt.result.kinetic &&
                     cached.result.tornElements == rebuilt.result.tornElements;
    double moved = 0;
    for (std::size_t i = 0; i < cached.position.size(); ++i) {
        identical = identical && cached.position[i] == rebuilt.position[i];
        moved = std::max(moved, std::abs(cached.position[i] - patch.mesh.position[i]));
    }
    expectTrue("caching the rest forms is bit-identical to rebuilding them", identical);

    bool preloadedIdentical =
        preloadedCache.result.work == preloadedRebuild.result.work &&
        preloadedCache.result.dissipation == preloadedRebuild.result.dissipation &&
        preloadedCache.result.strainEnergy == preloadedRebuild.result.strainEnergy &&
        preloadedCache.result.tornElements == preloadedRebuild.result.tornElements;
    for (std::size_t i = 0; i < preloadedCache.position.size(); ++i)
        preloadedIdentical =
            preloadedIdentical && preloadedCache.position[i] == preloadedRebuild.position[i];
    expectTrue("and bit-identical under a pre-load, where rest and position differ",
               preloadedIdentical);
    // Guard: the pre-load has to have actually moved the rest configuration, or the
    // case above is the case below it wearing a hat.
    expectTrue("and the pre-load stored energy, so rest really did move",
               preloadedCache.result.initialStrainEnergy > 0.0);
    expectTrue("and moved the answer, so it was not a no-op",
               preloadedCache.result.work != cached.result.work);
    std::printf("     %d steps, %.4f m of travel, %d element(s) torn, peak %.2f MN\n",
                cached.result.steps, moved, cached.result.tornElements,
                cached.result.peakForce / 1e6);
    expectTrue("and the run deformed the patch, so there was something to compare",
               moved > 1e-3);
    expectTrue("and tore it, so the element-deletion path was compared too",
               cached.result.tornElements > 0);
}

// `Solver::adopt` exists so a state computed elsewhere -- a GPU run -- is read by
// the energy account, the tearing rules and the panel reporting that are already
// validated, rather than by a second copy of them. It shipped with no test of its
// own, on the caller's path, which is the failure `CLAUDE.md` records as "asking
// what in the diff was not exercised".
//
// The statement it has to satisfy is an identity, not a tolerance: adopting the
// state a run *did* produce must reproduce that run's own report.
void testAdoptingAStateReproducesTheRunThatMadeIt() {
    std::printf("\n   adopting a state computed elsewhere\n");
    const StructuralMesh strip = flatStrip(0.6, 0.3, 0.006, 3, 1);
    const zone::Patch patch = zone::buildPatch(strip, {0.05, 0, 0}, flatParams(2));

    zone::SolveParams solve;
    solve.indenter.halfLength = 0.03;
    solve.indenter.halfWidth = 0.08;
    solve.indenter.speed = 20.0;
    solve.indenter.rampTime = 5.0e-4;
    solve.indenter.stopAt = 0.07;
    zone::Solver ran(patch, plasticity::shipSteel(), solve);
    const zone::SolveResult& original = ran.run();

    zone::Solver fresh(patch, plasticity::shipSteel(), solve);
    fresh.adopt(ran.position(), ran.velocity(), ran.elementState(), original.steps, original.time,
                original.penetration, original.work, original.dissipation);
    const zone::SolveResult& adopted = fresh.result();

    std::printf("     %d steps, %d torn, strain %.4f MJ; adopted reports %d torn, %.4f MJ\n",
                original.steps, original.tornElements, original.strainEnergy / 1e6,
                adopted.tornElements, adopted.strainEnergy / 1e6);
    expectEqual("adopting reports the same step count",
                static_cast<long long>(adopted.steps), static_cast<long long>(original.steps));
    expectTrue("and the same torn element count", adopted.tornElements == original.tornElements);
    expectEqualCount("and the same torn panels", adopted.tornPanels.size(),
                     original.tornPanels.size());
    // The strain energy agrees closely but **not bit for bit**, and the reason is a
    // property of `step()` rather than of `adopt()`: a step runs
    // computeForces -> integrate -> accumulateEnergy, so the energy a run finishes
    // with is the stress from *before* its last integration, one step stale against
    // the positions it also reports. `adopt` recomputes it from the positions it was
    // handed, so it is the un-stale version of the same quantity. A tolerance of one
    // step's worth is the honest statement; equality would be asserting a defect.
    expectTrue("and a strain energy within one step of it",
               std::abs(adopted.strainEnergy - original.strainEnergy) <
                   0.05 * std::abs(original.strainEnergy));
    expectTrue("and the same torn area", adopted.tornArea == original.tornArea);
    bool samePosition = true;
    for (std::size_t i = 0; i < ran.position().size(); ++i)
        samePosition = samePosition && fresh.position()[i] == ran.position()[i];
    expectTrue("and holds the positions it was given", samePosition);
    // Guards. A run that neither moved nor tore makes every line above vacuous, and
    // a `fresh` that had happened to reach the same state on its own would too --
    // it has taken no steps at all, so its strain energy before adopting is the
    // pre-load's, which here is zero.
    expectTrue("the run tore something, so the tearing path was reproduced",
               original.tornElements > 0);
    expectTrue("and stored strain energy, so the quadrature had something to do",
               original.strainEnergy > 0.0);
    // And adopting must not *advance* the state it was given: the plastic history
    // comes back exactly as handed over, not one increment further on. Recovering
    // the stress runs the element update, which commits, so this is a real hazard
    // and not a hypothetical -- a point just under its damage limit would be tipped
    // over it by the act of being read.
    bool historyUntouched = true;
    for (std::size_t e = 0; e < patch.elementCount(); ++e)
        for (int gp = 0; gp < solidshell::kGauss; ++gp)
            historyUntouched =
                historyUntouched &&
                fresh.elementState()[e].point[gp].equivalentPlasticStrain ==
                    ran.elementState()[e].point[gp].equivalentPlasticStrain &&
                fresh.elementState()[e].point[gp].damage ==
                    ran.elementState()[e].point[gp].damage;
    expectTrue("and adopting does not advance the history it was handed",
               historyUntouched);
}

}  // namespace

void runZoneTests() {
    std::printf("\n--- damage zone: structure to torn panels ---\n");
    testTheMesherIsExactOnFlatPlating();
    testTheMesherReportsWhatCurvatureCosts();
    testTheZoneStopsAtASeamAndAtAFold();
    testAMeshedPatchReproducesAClampedPlate();
    testAPatchAtRestStaysAtRest();
    testRigidTranslationCarriesNoForce();
    testTheEnergyAccountCloses();
    testAgainstTheMembraneModel();
    testTearingReportsPanelIndicesForBreach();
    testTheParallelSolveIsBitIdentical();
    testAZoneOnTheFerry();
    testStiffenersBracketTheAnswer();
    testTheWeldIsADistanceAndNotABucket();
    testTheAxisPointsOutward();
    testAreasAreRightOnATrapezoid();
    testNodalNormalsAreAreaWeighted();
    testTheCostEstimateIsAFunctionOfTheRightThings();
    testStiffenersHoldOnlyWhatTheyReach();
    testTheRestFormsCacheChangesNothing();
    testAdoptingAStateReproducesTheRunThatMadeIt();
}
