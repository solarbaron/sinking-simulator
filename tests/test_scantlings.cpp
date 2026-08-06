// SPDX-License-Identifier: MIT
//
// Validation of scantlings -> structural mesh.
//
// Naval architecture hands out exact answers here, so nothing below is eyeballed.
// A plate-and-bar section is an integral anyone can do by hand; a box hull's shell
// area is L*B + 2*L*D and its steel weight follows from it; the neutral axis of a
// section made of three rectangles has one value and it is not a matter of taste.
//
// Three instruments do most of the work:
//
//   * **A box hull**, where every quantity is closed form and every panel is
//     planar, so tiling and weight are *exact* rather than approximate. This is
//     what caught the section routine counting both bays either side of a frame
//     station, which doubled area and second moment while leaving the neutral
//     axis -- a ratio -- looking perfectly correct.
//   * **Quadrature over the section's own width profile**, which reaches the same
//     second moment without ever forming a parallel-axis term, so agreement means
//     the algebra is right rather than merely self-consistent. The same device
//     `tests/test_hullform.cpp` uses on the area curve.
//   * **Refinement**, because a small constant error and a discretisation error
//     look identical at one resolution.
#include "engine/sim/scantlings.hpp"
#include "engine/core/geometry.hpp"
#include "game/prototype/ferry.hpp"
#include "harness.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

using namespace sim;
using testing::expectEqual;
using testing::expectNear;
using testing::expectTrue;

namespace {

// --- Shared fixtures ----------------------------------------------------------

// A rectangular box hull: L x B x D with flat sides, a flat bottom and a flat
// deck. Every area, weight and section property of a structure built on it is a
// closed form with no tessellation error at all, which is the only way to tell a
// tiling *approximation* from a tiling *mistake*.
constexpr double kBoxLength = 100.0;
constexpr double kBoxBeam = 20.0;
constexpr double kBoxDepth = 10.0;

TriMesh boxHull(int stations = 5, int waterlines = 3) {
    std::vector<double> levels;
    for (int k = 0; k < waterlines; ++k)
        levels.push_back(kBoxDepth * k / (waterlines - 1));
    std::vector<Station> out;
    for (int i = 0; i < stations; ++i) {
        Station s;
        s.x = -0.5 * kBoxLength + kBoxLength * i / (stations - 1);
        s.halfBeam.assign(levels.size(), 0.5 * kBoxBeam);
        out.push_back(s);
    }
    return makeHullFromStations(out, levels);
}

// Shell area of a hull mesh, measured straight off its triangles: everything
// except the flat caps at the extreme stations and the weather deck. The keel
// strip is *not* a cap -- it is the flat of bottom, and it is plating.
//
// Deliberately independent of anything in scantlings.cpp, so it is a second
// opinion rather than an echo.
double hullShellArea(const TriMesh& mesh) {
    double xLo = 1e300, xHi = -1e300, zHi = -1e300;
    for (const Vec3& v : mesh.verts) {
        xLo = std::min(xLo, v.x);
        xHi = std::max(xHi, v.x);
        zHi = std::max(zHi, v.z);
    }
    double total = 0;
    for (const Tri& t : mesh.tris) {
        const Vec3 a = mesh.verts[t.a], b = mesh.verts[t.b], c = mesh.verts[t.c];
        const auto flatAt = [&](double Vec3::*axis, double value) {
            return std::abs(a.*axis - value) < 1e-9 && std::abs(b.*axis - value) < 1e-9 &&
                   std::abs(c.*axis - value) < 1e-9;
        };
        if (flatAt(&Vec3::x, xLo) || flatAt(&Vec3::x, xHi) || flatAt(&Vec3::z, zHi)) continue;
        total += 0.5 * length(cross(b - a, c - a));
    }
    return total;
}

// Uniform plating over the whole girth, with everything optional switched off, so
// a test can add back exactly the one thing it means to measure.
Scantlings bareBoxScantlings(double thickness = 0.020) {
    Scantlings s;
    s.frameSpacing = 5.0;
    s.longitudinalSpacing = 1.0;
    s.framed = false;
    ShellRegion r;
    r.name = "all";
    r.thickness = thickness;
    r.stiffened = false;
    s.shell = {r};
    return s;
}

double shellPanelArea(const StructuralMesh& mesh) {
    double total = 0;
    for (const PlatePanel& p : mesh.panels)
        if (p.role == PanelRole::Shell) total += p.area();
    return total;
}

// --- Section properties, against hand integration -----------------------------

// Second moments of a stack of rectangles, by midpoint quadrature over the
// section's own width profile. The strips are aligned to the material
// boundaries, so the only error is the quadrature's own, and -- crucially -- the
// route never forms a parallel-axis term. It computes the raw moments about zero
// and subtracts A*NA^2 once at the end, which is a different piece of algebra
// from the composite sum the implementation does.
struct Strip {
    double from = 0, to = 0, width = 0;
};

void integrateStrips(const std::vector<Strip>& strips, double& area, double& neutralAxis,
                     double& secondMoment) {
    constexpr int kSteps = 200000;
    double a = 0, q = 0, m = 0;
    for (const Strip& s : strips) {
        const double h = (s.to - s.from) / kSteps;
        for (int i = 0; i < kSteps; ++i) {
            const double z = s.from + h * (i + 0.5);
            const double dA = s.width * h;
            a += dA;
            q += dA * z;
            m += dA * z * z;
        }
    }
    area = a;
    neutralAxis = a > 0 ? q / a : 0.0;
    secondMoment = m - a * neutralAxis * neutralAxis;
}

// A tee bar on its attached plating. Every number below is an integral a naval
// architect does on paper, so the assertion is the analytic value and not a
// recorded output.
void testStiffenedSectionAgainstHandIntegration() {
    const double webHeight = 0.300, webThickness = 0.012;
    const double flangeWidth = 0.100, flangeThickness = 0.016;
    const double plateThickness = 0.014, plateWidth = 0.700;
    const StiffenerProfile profile = tee(webHeight, webThickness, flangeWidth, flangeThickness);
    const StiffenedSection s = stiffenedSection(profile, plateThickness, plateWidth);

    // Route 1: quadrature over the width profile, measured from the plate
    // mid-surface, which is where the header says the datum is.
    const double half = 0.5 * plateThickness;
    double qArea = 0, qNeutral = 0, qSecond = 0;
    integrateStrips({{-half, half, plateWidth},
                     {half, half + webHeight, webThickness},
                     {half + webHeight, half + webHeight + flangeThickness, flangeWidth}},
                    qArea, qNeutral, qSecond);
    expectNear("stiffened section area matches quadrature", s.area, qArea, 1e-9 * qArea);
    expectNear("stiffened section neutral axis matches quadrature", s.neutralAxis, qNeutral,
               1e-9 * qNeutral);
    expectNear("stiffened section second moment matches quadrature", s.secondMoment, qSecond,
               1e-8 * qSecond);

    // Route 2: the composite sum, written out longhand. Two independent routes to
    // one number, as tests/test_hullform.cpp does with the area curve.
    const double plateArea = plateWidth * plateThickness;
    const double webArea = webHeight * webThickness;
    const double flangeArea = flangeWidth * flangeThickness;
    const double webCentre = half + 0.5 * webHeight;
    const double flangeCentre = half + webHeight + 0.5 * flangeThickness;
    const double totalArea = plateArea + webArea + flangeArea;
    const double neutral = (webArea * webCentre + flangeArea * flangeCentre) / totalArea;
    const double second =
        plateWidth * plateThickness * plateThickness * plateThickness / 12.0 +
        plateArea * neutral * neutral +
        webThickness * webHeight * webHeight * webHeight / 12.0 +
        webArea * (webCentre - neutral) * (webCentre - neutral) +
        flangeWidth * flangeThickness * flangeThickness * flangeThickness / 12.0 +
        flangeArea * (flangeCentre - neutral) * (flangeCentre - neutral);
    expectNear("area equals the sum of three rectangles", s.area, totalArea, 1e-12);
    expectNear("neutral axis equals the hand composite", s.neutralAxis, neutral, 1e-12);
    expectNear("second moment equals the hand composite", s.secondMoment, second, 1e-15);

    // And the section moduli, which are what a rule check actually reads.
    const double top = half + webHeight + flangeThickness;
    expectNear("section modulus to the flange", s.modulusStiffener, second / (top - neutral),
               1e-12);
    expectNear("section modulus to the plate", s.modulusPlate, second / (neutral + half), 1e-12);

    // Guards. The plate is the larger area but the stiffener is what moves the
    // neutral axis, so a section whose neutral axis sat in the plate -- which is
    // what a smeared model produces -- would pass a sloppier version of this test.
    expectTrue("the neutral axis is lifted clear of the plating", s.neutralAxis > 3.0 * half);
    expectTrue("the two section moduli genuinely differ", s.modulusPlate > 2.0 * s.modulusStiffener);
    expectTrue("the plate is the larger area", plateArea > webArea + flangeArea);
}

// A flat bar has no flange, so its section is one rectangle and every property is
// a schoolbook formula. If this is wrong nothing above it can be right.
void testFlatBarProfileIsTheSchoolbookRectangle() {
    const double h = 0.200, t = 0.010;
    const ProfileSection s = profileSection(flatBar(h, t));
    expectNear("flat bar area", s.area, h * t, 1e-15);
    expectNear("flat bar centroid is at half its height", s.centroid, 0.5 * h, 1e-15);
    expectNear("flat bar height", s.height, h, 1e-15);
    expectNear("flat bar strong axis is t h^3 / 12", s.secondMoment, t * h * h * h / 12.0, 1e-18);
    expectNear("flat bar weak axis is h t^3 / 12", s.secondMomentWeak, h * t * t * t / 12.0, 1e-18);
    expectNear("flat bar torsion constant is h t^3 / 3", s.torsionConstant, h * t * t * t / 3.0,
               1e-18);
    expectTrue("the strong axis really is the strong one",
               s.secondMoment > 100.0 * s.secondMomentWeak);
}

// An angle and a tee of the same dimensions bend identically in the plane of the
// web -- the flange sits at the same height either way -- and differently about
// the axis at right angles to it. A generator that treated them as the same
// profile would pass every hull girder test in this file, so the difference is
// asserted where it lives.
//
// The direction of that difference is the opposite of the intuition that an
// unsymmetrical section must be weaker: about the *geometric* axis parallel to
// the web, an angle's flange hangs off to one side and so sits further from the
// section centroid than a tee's straddling one, and its second moment is
// therefore larger. It is the principal minimum that is smaller, and that needs
// the product of inertia, which ProfileSection deliberately does not carry.
void testAngleAndTeeDifferAboutTheAxisAcrossTheWeb() {
    const double h = 0.200, tw = 0.010, bf = 0.090, tf = 0.012;
    const ProfileSection t = profileSection(tee(h, tw, bf, tf));
    const ProfileSection a = profileSection(angle(h, tw, bf, tf));
    expectNear("same area", a.area, t.area, 1e-15);
    expectNear("same centroid", a.centroid, t.centroid, 1e-15);
    expectNear("same strong axis", a.secondMoment, t.secondMoment, 1e-15 + 1e-12 * t.secondMoment);

    // Both are hand-computable. The tee's parts are centred on the web, so no
    // parallel-axis term arises at all.
    expectNear("a tee's cross-web moment is the two rectangles about the web",
               t.secondMomentWeak, h * tw * tw * tw / 12.0 + tf * bf * bf * bf / 12.0, 1e-15);
    // The angle's flange runs from the web out to bf, so the section centroid is
    // offset and both parts pick one up.
    const double webArea = h * tw, flangeArea = bf * tf;
    const double centre = flangeArea * 0.5 * bf / (webArea + flangeArea);
    const double expected = h * tw * tw * tw / 12.0 + webArea * centre * centre +
                            tf * bf * bf * bf / 12.0 +
                            flangeArea * (0.5 * bf - centre) * (0.5 * bf - centre);
    expectNear("an angle's is the same two rectangles about their own centroid",
               a.secondMomentWeak, expected, 1e-15);
    expectTrue("which is the larger of the two", a.secondMomentWeak > 2.0 * t.secondMomentWeak);
    expectTrue("and both are non-trivial", t.secondMomentWeak > 1e-9);
}

// The reason stiffeners are discrete elements and not smeared into the plate.
//
// A smeared panel has the same area and the same axial stiffness by construction;
// what it loses is the stiffener's lever arm, and with it the panel's bending
// stiffness. The header claims a factor of 130 for a representative panel. That
// is a claim about the physics, so it is measured rather than asserted in prose.
void testSmearingLosesTheStiffenersLeverArm() {
    const StiffenerProfile profile = flatBar(0.200, 0.010);
    const double plate = 0.012, spacing = 0.700;

    const StiffenedSection discrete = stiffenedSection(profile, plate, spacing);
    const double equivalent = smearedThickness(plate, profile, spacing);
    // The homogenised panel is a bare plate of the equivalent thickness.
    const StiffenedSection smeared = stiffenedSection(flatBar(0.0, 0.0), equivalent, spacing);

    expectNear("smearing preserves the cross-sectional area", smeared.area, discrete.area,
               1e-12 * discrete.area);
    expectNear("smearing throws away the bending stiffness, by a factor of 130",
               discrete.secondMoment / smeared.secondMoment, 130.25, 0.05);
    expectNear("and puts the neutral axis back in the plate", smeared.neutralAxis, 0.0, 1e-12);
    expectTrue("where the real one is not", discrete.neutralAxis > 0.5 * plate);

    // Guard against the comparison being vacuous: the equivalent thickness has to
    // be a real thickening, not a rounding of the plate.
    expectNear("the equivalent thickness is t + A/s", equivalent,
               plate + profileSection(profile).area / spacing, 1e-15);
    expectTrue("and it is materially thicker than the bare plate", equivalent > 1.2 * plate);
}

// A stiffener's eccentricity has a sign. Welding the same bar to the other face
// of the same plate moves the neutral axis the other way, and any model that
// carried only the profile's area would report both as identical.
void testStiffenerEccentricityHasASign() {
    const TriMesh hull = boxHull();
    Scantlings up = bareBoxScantlings();
    up.shell[0].stiffened = true;
    up.shell[0].longitudinal = flatBar(0.250, 0.012);

    const StructuralMesh mesh = makeStructuralMesh(hull, up);
    const HullGirderSection section = hullGirderSection(mesh, 2.5);

    // Same structure with every web pointing the other way -- outboard, through
    // the shell, which is not a ship but is the control this needs.
    StructuralMesh flipped = mesh;
    for (StructuralMember& m : flipped.members) m.rise = -m.rise;
    const HullGirderSection other = hullGirderSection(flipped, 2.5);

    expectNear("flipping the webs leaves the material area alone", other.area, section.area,
               1e-12 * section.area);
    expectTrue("but moves the neutral axis", std::abs(other.neutralAxis - section.neutralAxis) >
                                                 1e-3);
    expectTrue("and changes the second moment",
               std::abs(other.secondMoment - section.secondMoment) > 1e-3 * section.secondMoment);
}

// A panel is a quad and a quad on a curved shell is a trapezoid at best and
// warped at worst, so its area is the sum of two triangles rather than anything
// simpler. Doubling one triangle is exactly right for a rectangle and wrong for
// everything else, which is why it survives a suite that only ever builds
// rectangles.
void testPanelAreaOfATrapezoidAndAWarpedQuad() {
    PlatePanel trapezoid;
    trapezoid.corner[0] = {0, 0, 0};
    trapezoid.corner[1] = {10, 0, 0};
    trapezoid.corner[2] = {10, 4, 0};
    trapezoid.corner[3] = {0, 2, 0};
    expectNear("a trapezoid's area is its mean width times its length", trapezoid.area(),
               10.0 * 0.5 * (2.0 + 4.0), 1e-12);
    // Guard: the doubled-single-triangle answer must differ, or the assertion
    // above is satisfied by the wrong formula too.
    expectTrue("and that is not twice one triangle",
               std::abs(trapezoid.area() - 2.0 * 0.5 * 10.0 * 4.0) > 1.0);

    PlatePanel warped;
    warped.corner[0] = {0, 0, 0};
    warped.corner[1] = {1, 0, 0};
    warped.corner[2] = {1, 1, 1};
    warped.corner[3] = {0, 1, 0};
    expectNear("a warped quad is two triangles", warped.area(), std::sqrt(2.0), 1e-12);
    expectNear("its centroid is the area-weighted mean of the two", warped.centroid().z, 1.0 / 3.0,
               1e-12);

    // The normal comes off the diagonals, which is the only construction that
    // stays sensible on a warped quad -- one triangle's edges give whichever of
    // the two facets happens to be first.
    PlatePanel flat;
    flat.corner[0] = {0, 0, 0};
    flat.corner[1] = {2, 0, 0};
    flat.corner[2] = {2, 3, 0};
    flat.corner[3] = {0, 3, 0};
    const Vec3 n = flat.normal();
    expectNear("a flat panel's normal is a unit vector", length(n), 1.0, 1e-12);
    expectNear("and points along the axis it should", std::abs(n.z), 1.0, 1e-12);
    const Vec3 warpedNormal = warped.normal();
    expectNear("a warped panel's normal is still a unit vector", length(warpedNormal), 1.0, 1e-12);
    expectNear("and is perpendicular to both diagonals",
               dot(warpedNormal, warped.corner[2] - warped.corner[0]), 0.0, 1e-12);
    expectNear("both of them", dot(warpedNormal, warped.corner[3] - warped.corner[1]), 0.0, 1e-12);
}

// --- Tiling -------------------------------------------------------------------

// A box hull's shell is L*B + 2*L*D exactly: flat faces, no tessellation error,
// no argument. Panels that tile it must sum to that and no more, and -- the part
// that distinguishes tiling from a lucky total -- every panel must be the same
// bay by band rectangle, and there must be exactly as many of them as the
// spacings imply.
void testBoxPanelsTileExactly() {
    const TriMesh hull = boxHull();
    expectTrue("the box hull is a closed manifold", isClosedManifold(hull));

    const double analytic = kBoxLength * kBoxBeam + 2.0 * kBoxLength * kBoxDepth;
    expectNear("the mesh's own shell area is the closed form", hullShellArea(hull), analytic,
               1e-9 * analytic);

    const Scantlings s = bareBoxScantlings();
    std::vector<std::string> problems;
    const StructuralMesh mesh = makeStructuralMesh(hull, s, &problems);
    expectTrue("an ordinary box raises no problems", problems.empty());

    // Girth per side is B/2 across the bottom plus D up the side; the layout puts
    // one band per metre of it.
    const int bays = static_cast<int>(std::lround(kBoxLength / s.frameSpacing));
    const int bands =
        static_cast<int>(std::lround((0.5 * kBoxBeam + kBoxDepth) / s.longitudinalSpacing));
    expectEqual("one panel per bay, band and side", static_cast<long long>(mesh.panels.size()),
                2LL * bays * bands);
    expectEqual("frame stations bound the bays",
                static_cast<long long>(mesh.frameStations.size()), bays + 1LL);

    const double expectedPanel = (kBoxLength / bays) * ((0.5 * kBoxBeam + kBoxDepth) / bands);
    double worst = 0;
    for (const PlatePanel& p : mesh.panels) worst = std::max(worst, std::abs(p.area() - expectedPanel));
    expectTrue("every panel is the same bay by band rectangle", worst < 1e-5 * expectedPanel);

    expectNear("the panels sum to the hull's shell area", mesh.plateArea(), analytic,
               1e-6 * analytic);

    // Every corner is on the box: either on the bottom or on a side.
    bool onSurface = true;
    for (const PlatePanel& p : mesh.panels)
        for (const Vec3& c : p.corner) {
            const bool bottom = std::abs(c.z) < 1e-4 && std::abs(c.y) <= 0.5 * kBoxBeam + 1e-6;
            const bool side = std::abs(std::abs(c.y) - 0.5 * kBoxBeam) < 1e-4 && c.z >= -1e-6 &&
                              c.z <= kBoxDepth + 1e-6;
            if (!bottom && !side) onSurface = false;
        }
    expectTrue("every panel corner lies on the hull surface", onSurface);
}

// Steel weight, against plate area times thickness times density and stiffener
// length times section area times density, both summed from the closed form
// rather than from the mesh. Everything about the box is exact, so this is an
// equality and not an estimate.
void testBoxSteelWeightIsTheClosedForm() {
    const TriMesh hull = boxHull();
    Scantlings s = bareBoxScantlings(0.020);
    s.framed = true;
    s.frameProfile = flatBar(0.300, 0.012);
    s.shell[0].stiffened = true;
    s.shell[0].longitudinal = flatBar(0.200, 0.010);

    const StructuralMesh mesh = makeStructuralMesh(hull, s);
    const double density = 7850.0;

    const double shell = kBoxLength * kBoxBeam + 2.0 * kBoxLength * kBoxDepth;
    expectNear("plating weight is area x thickness x density", mesh.plateMass(),
               shell * 0.020 * density, 1e-5 * shell * 0.020 * density);

    // Longitudinals: one on every interior girth seam, both sides, full length.
    const double girth = 0.5 * kBoxBeam + kBoxDepth;
    const int bands = static_cast<int>(std::lround(girth / s.longitudinalSpacing));
    const double longitudinalLength = 2.0 * (bands - 1) * kBoxLength;
    // Frames: the whole girth of both sides, at every frame station.
    const int bays = static_cast<int>(std::lround(kBoxLength / s.frameSpacing));
    const double frameLength = (bays + 1) * 2.0 * girth;

    expectEqual("one longitudinal segment per bay, seam and side",
                static_cast<long long>(mesh.memberCount(MemberRole::Longitudinal)),
                2LL * bays * (bands - 1));
    expectEqual("one frame segment per station, band and side",
                static_cast<long long>(mesh.memberCount(MemberRole::Frame)),
                2LL * (bays + 1) * bands);
    expectNear("total member length is the closed form", mesh.memberLength(),
               longitudinalLength + frameLength, 1e-5 * (longitudinalLength + frameLength));

    const double expected = longitudinalLength * profileSection(s.shell[0].longitudinal).area +
                            frameLength * profileSection(s.frameProfile).area;
    expectNear("stiffener weight is length x section area x density", mesh.memberMass(),
               expected * density, 1e-5 * expected * density);

    // Guard: the stiffening must be a real fraction of the whole, or the plating
    // term alone would carry the test.
    expectTrue("stiffening is a meaningful share of the steel",
               mesh.memberMass() > 0.10 * mesh.plateMass());
    expectNear("the two parts sum to the whole", mesh.steelMass(),
               mesh.plateMass() + mesh.memberMass(), 1e-6 * mesh.steelMass());
}

// On a curved hull a panel is a flat chord across a curved shell, so the panels
// sum to slightly *less* than the surface they cover. That is a discretisation
// error rather than a mistake only if it shrinks under refinement -- and the
// direction it shrinks from matters too, because a chord can only lose area.
void testFerryPanelsApproachTheHullAreaFromBelow() {
    const TriMesh hull = game::buildFerry().hull;
    const double shell = hullShellArea(hull);
    expectTrue("the ferry hull has a shell to cover", shell > 4000.0);

    double previousError = 1e30;
    int improved = 0;
    double finestError = 1.0;
    for (double refine : {1.0, 0.5, 0.25, 0.125}) {
        Scantlings s = ferryScantlings();
        s.frameSpacing *= refine;
        s.longitudinalSpacing *= refine;
        for (ShellRegion& r : s.shell) r.longitudinalSpacing *= refine;
        const StructuralMesh mesh = makeStructuralMesh(hull, s);

        const double area = shellPanelArea(mesh);
        expectTrue("chorded panels never exceed the surface they cover", area < shell + 1e-6);
        const double error = (shell - area) / shell;
        // A rate, not just a direction: a sequence that merely improved could
        // still be converging to the wrong answer.
        if (error < 0.8 * previousError) ++improved;
        previousError = error;
        finestError = error;
    }
    expectEqual("every halving cuts the missing area by at least a fifth", improved, 4);
    expectTrue("and the finest is within a tenth of a percent", finestError < 1e-3);

    // Guard: the coarsest must be materially worse, or "converges" would pass on
    // four identical numbers.
    Scantlings coarse = ferryScantlings();
    coarse.longitudinalSpacing = 2.8;
    for (ShellRegion& r : coarse.shell) r.longitudinalSpacing = 2.8;
    const double coarseArea = shellPanelArea(makeStructuralMesh(hull, coarse));
    expectTrue("the coarse mesh really does lose area", (shell - coarseArea) / shell > 4e-3);
}

// Distance from a point to the hull's triangulated surface: the smallest
// perpendicular distance to any triangle whose footprint contains the point's
// projection. Slower than a ray cast and completely unlike one, which is the
// point -- it checks the sampler's output without repeating its method.
double distanceToMesh(const TriMesh& mesh, const Vec3& p) {
    double best = 1e300;
    for (const Tri& t : mesh.tris) {
        const Vec3 a = mesh.verts[t.a], b = mesh.verts[t.b], c = mesh.verts[t.c];
        const Vec3 n = cross(b - a, c - a);
        const double nn = dot(n, n);
        if (nn < 1e-24) continue;
        const double offset = dot(n, p - a) / std::sqrt(nn);
        if (std::abs(offset) >= best) continue;
        const Vec3 projected = p - n * (dot(n, p - a) / nn);
        const double wc = dot(cross(b - a, projected - a), n) / nn;
        const double wb = dot(cross(projected - a, c - a), n) / nn;
        if (wb < -1e-9 || wc < -1e-9 || wb + wc > 1.0 + 1e-9) continue;
        best = std::abs(offset);
    }
    return best;
}

// Every shell panel corner has to sit on the hull, not near it.
//
// It does not sit on it exactly, and the reason is worth stating: a corner at a
// given girth fraction is interpolated along a polyline of sampled heights, and
// where two consecutive samples straddle a waterline or a triangulation diagonal
// the chord between them leaves the surface. So the right assertion is not a
// fixed tolerance but a *convergence*: the deviation must fall as the section is
// sampled more finely, which a mis-aimed ray or a wrong interpolation would not
// do. Measured here at 41, 28, 12, 3.8 and 1.0 mm.
void testPanelCornersConvergeOntoTheHull() {
    const TriMesh hull = game::buildFerry().hull;

    double previous = 1e300;
    int improved = 0;
    double finest = 1e300;
    double outermost = 0;
    for (int samples : {48, 96, 192, 384, 768}) {
        Scantlings s = ferryScantlings();
        s.girthSamples = samples;
        const StructuralMesh mesh = makeStructuralMesh(hull, s);

        double worst = 0;
        // Every seventh panel: the check is quadratic in the hull mesh and the
        // sample is far larger than it needs to be to find a systematic error.
        for (std::size_t i = 0; i < mesh.panels.size(); i += 7) {
            if (mesh.panels[i].role != PanelRole::Shell) continue;
            for (const Vec3& c : mesh.panels[i].corner)
                worst = std::max(worst, distanceToMesh(hull, c));
            for (const Vec3& c : mesh.panels[i].corner) outermost = std::max(outermost, std::abs(c.y));
        }
        if (worst < previous) ++improved;
        previous = worst;
        finest = worst;
    }
    expectEqual("finer sampling puts the corners closer to the hull every time", improved, 5);
    expectTrue("and the finest is inside a millimetre and a half", finest < 1.5e-3);
    expectNear("the widest corner reaches the moulded half beam exactly", outermost, 10.0, 1e-9);
}

// Frame spacing is nominal: the hull is divided into a whole number of bays so
// that the ends land on frames instead of on a ragged short one. The rounding
// has to go to the *nearest* whole number and the difference has to be reported,
// because a delivered spacing four percent off the requested one is a different
// ship's structure and nothing else would say so.
void testFrameSpacingIsRoundedToAWholeNumberOfBays() {
    const TriMesh hull = boxHull();
    Scantlings s = bareBoxScantlings();
    s.frameSpacing = 8.0;   // 100 / 8 = 12.5 bays, which rounds up rather than down

    std::vector<std::string> problems;
    const StructuralMesh mesh = makeStructuralMesh(hull, s, &problems);
    expectEqual("12.5 bays becomes 13, not 12", static_cast<long long>(mesh.frameStations.size()),
                14);
    expectNear("and the spacing delivered says so", mesh.frameSpacing, kBoxLength / 13.0, 1e-12);
    expectTrue("a spacing that had to move is reported", !problems.empty());
    expectNear("the stations span the whole hull", mesh.frameStations.front(),
               -0.5 * kBoxLength, 1e-12);
    expectNear("both ends", mesh.frameStations.back(), 0.5 * kBoxLength, 1e-12);

    // A spacing that divides exactly is delivered exactly, and silently.
    std::vector<std::string> quiet;
    Scantlings exact = bareBoxScantlings();
    const StructuralMesh clean = makeStructuralMesh(hull, exact, &quiet);
    expectNear("a spacing that divides is delivered as asked", clean.frameSpacing,
               exact.frameSpacing, 1e-12);
    expectTrue("without a word", quiet.empty());
}

// Later regions override earlier ones, so a thickened patch can be declared on
// top of the strake it sits in rather than having to be cut out of it. The
// reference ferry relies on this for its thinner end plating, and getting it
// backwards changes the steel weight by only a percent or two -- far too little
// for a weight check to notice.
void testALaterShellRegionOverridesAnEarlierOne() {
    const TriMesh hull = boxHull();
    Scantlings s = bareBoxScantlings(0.010);
    ShellRegion patch;
    patch.name = "midship_thickening";
    patch.thickness = 0.030;
    patch.xFrom = -25.0;
    patch.xTo = 25.0;
    patch.stiffened = false;
    s.shell.push_back(patch);

    const StructuralMesh mesh = makeStructuralMesh(hull, s);
    // Girth of both sides is B + 2D per metre of length; the patch covers the
    // middle 50 m exactly, because the frame stations fall on -25 and 25.
    const double perMetre = kBoxBeam + 2.0 * kBoxDepth;
    const double expected = (50.0 * perMetre * 0.010 + 50.0 * perMetre * 0.030) * 7850.0;
    expectNear("the patch's thickness applies where it was declared", mesh.plateMass(), expected,
               1e-5 * expected);

    // Guard: without the override the answer is a long way off, so this is not a
    // test that passes on any thickness at all.
    const double unpatched = 100.0 * perMetre * 0.010 * 7850.0;
    expectTrue("and the un-patched weight is materially different",
               std::abs(expected - unpatched) > 0.5 * unpatched);
}

// A deck is bounded by the hull, not by its own stiffener spacing. The last band
// outboard is whatever is left, so the deck's area is the hull's breadth however
// the spacing happens to divide it.
void testDeckPlatingIsClippedToTheHull() {
    const TriMesh hull = boxHull();
    Scantlings s = bareBoxScantlings();
    s.materials = {ah36Steel(), mildSteel()};
    Deck deck;
    deck.name = "weather";
    deck.z = kBoxDepth;
    deck.xFrom = -0.5 * kBoxLength;
    deck.xTo = 0.5 * kBoxLength;
    deck.thickness = 0.010;
    deck.longitudinalSpacing = 0.7;   // 10 / 0.7 is not a whole number of bands
    deck.material = 1;
    deck.stiffened = false;
    deck.beam = flatBar(0.300, 0.010);
    deck.beamed = true;
    s.decks = {deck};

    const StructuralMesh mesh = makeStructuralMesh(hull, s);
    double area = 0;
    for (const PlatePanel& p : mesh.panels)
        if (p.role == PanelRole::Deck) area += p.area();
    expectNear("the deck covers the hull's breadth exactly", area, kBoxBeam * kBoxLength,
               1e-6 * kBoxBeam * kBoxLength);

    // Guard: the bands must genuinely overhang, or the clip is never exercised.
    const int bands = static_cast<int>(std::ceil(0.5 * kBoxBeam / deck.longitudinalSpacing));
    expectTrue("the outermost band really does overhang the side",
               bands * deck.longitudinalSpacing > 0.5 * kBoxBeam + 0.1);
    // And no panel corner is outside the hull.
    double widest = 0;
    for (const PlatePanel& p : mesh.panels)
        if (p.role == PanelRole::Deck)
            for (const Vec3& c : p.corner) widest = std::max(widest, std::abs(c.y));
    expectNear("no deck panel hangs over the side", widest, 0.5 * kBoxBeam, 1e-9);

    // Deck beams run right across, both sides of the centreline, at every station
    // the deck reaches. Two beams a station, each half the beam long.
    double beamLength = 0;
    for (const StructuralMember& m : mesh.members)
        if (m.role == MemberRole::DeckBeam) beamLength += m.length();
    const long long stations = static_cast<long long>(mesh.frameStations.size());
    expectEqual("two deck beams at every station",
                static_cast<long long>(mesh.memberCount(MemberRole::DeckBeam)), 2 * stations);
    expectNear("and together they span the full breadth at each",
               beamLength, static_cast<double>(stations) * kBoxBeam, 1e-6 * stations * kBoxBeam);

    // The material a deck names is the material it is built from. Both steels
    // here have the same density, so weight cannot tell them apart and nothing
    // but the index itself will.
    bool deckMaterial = true, shellMaterial = true;
    for (const PlatePanel& p : mesh.panels) {
        if (p.role == PanelRole::Deck && p.material != 1) deckMaterial = false;
        if (p.role == PanelRole::Shell && p.material != 0) shellMaterial = false;
    }
    expectTrue("deck panels carry the deck's material", deckMaterial);
    expectTrue("and shell panels the shell's", shellMaterial);
    expectNear("which weight could not have distinguished", mesh.materials[1].density,
               mesh.materials[0].density, 1e-12);
    expectTrue("though they are different steels",
               mesh.materials[1].yieldStrength != mesh.materials[0].yieldStrength);
}

// A bulkhead is bounded by the hull in both directions: clipped in height to the
// hull's own depth, and in breadth to the shell at each height. A structure that
// built the height it was given would hang panels in the air above the ship,
// which weighs something and looks like nothing.
void testBulkheadsAreClippedToTheHull() {
    const TriMesh hull = boxHull();
    Scantlings s = bareBoxScantlings();
    Bulkhead tall;
    tall.name = "over_tall";
    tall.transverse = true;
    tall.position = 0.0;
    tall.zFrom = -5.0;          // below the keel
    tall.zTo = 100.0;           // and far above the deck
    tall.thickness = 0.010;
    tall.stiffenerSpacing = 0.7;   // does not divide the half beam
    tall.stiffened = false;
    s.bulkheads = {tall};

    const StructuralMesh mesh = makeStructuralMesh(hull, s);
    double area = 0, highest = -1e300, lowest = 1e300, widest = 0;
    for (const PlatePanel& p : mesh.panels) {
        if (p.role != PanelRole::Bulkhead) continue;
        area += p.area();
        for (const Vec3& c : p.corner) {
            highest = std::max(highest, c.z);
            lowest = std::min(lowest, c.z);
            widest = std::max(widest, std::abs(c.y));
        }
    }
    expectNear("a bulkhead is the hull's section, however tall it was declared", area,
               kBoxBeam * kBoxDepth, 1e-6 * kBoxBeam * kBoxDepth);
    expectNear("clipped to the deck", highest, kBoxDepth, 1e-9);
    expectNear("and to the keel", lowest, 0.0, 1e-9);
    expectNear("and to the shell", widest, 0.5 * kBoxBeam, 1e-9);

    // A longitudinal bulkhead outside the shell builds nothing, and says so.
    Scantlings outside = bareBoxScantlings();
    Bulkhead wing;
    wing.name = "outside_the_shell";
    wing.transverse = false;
    wing.position = 0.5 * kBoxBeam + 2.0;
    wing.zFrom = 0.0;
    wing.zTo = kBoxDepth;
    wing.xFrom = -20.0;
    wing.xTo = 20.0;
    wing.thickness = 0.010;
    wing.stiffened = false;
    outside.bulkheads = {wing};
    std::vector<std::string> problems;
    const StructuralMesh none = makeStructuralMesh(hull, outside, &problems);
    expectEqual("a bulkhead outside the shell builds no panels",
                static_cast<long long>(none.panelCount(PanelRole::Bulkhead)), 0);
    expectTrue("and is reported", !problems.empty());

    // So does a girder, which is positioned by hand and so is the easiest thing
    // in the description to put in the wrong place.
    Scantlings stray = bareBoxScantlings();
    Girder girder;
    girder.name = "stray";
    girder.y = 0.5 * kBoxBeam + 2.0;
    girder.z = 0.0;
    girder.xFrom = -20.0;
    girder.xTo = 20.0;
    girder.profile = flatBar(0.500, 0.012);
    stray.girders = {girder};
    std::vector<std::string> strayProblems;
    const StructuralMesh nothing = makeStructuralMesh(hull, stray, &strayProblems);
    expectEqual("a girder outside the shell builds nothing",
                static_cast<long long>(nothing.memberCount(MemberRole::Girder)), 0);
    expectTrue("and is reported too", !strayProblems.empty());
}

// Two regions can claim the same strake over different stretches of the length,
// which is how a thinner end plating is described. The seams have to be
// continuous end to end, so one spacing has to win -- and it must be the
// *tightest*, because a stiffener spacing is a maximum and exceeding it anywhere
// is a scantling failure.
void testTheTightestSpacingWinsAStrake() {
    const TriMesh hull = boxHull();
    Scantlings s = bareBoxScantlings();
    s.shell[0].longitudinalSpacing = 1.0;
    s.shell[0].xTo = 0.0;
    ShellRegion forward = s.shell[0];
    forward.name = "forward";
    forward.xFrom = 0.0;
    forward.xTo = 1e9;
    forward.longitudinalSpacing = 0.5;
    s.shell.push_back(forward);

    const StructuralMesh mesh = makeStructuralMesh(hull, s);
    const double girth = 0.5 * kBoxBeam + kBoxDepth;
    const int tight = static_cast<int>(std::lround(girth / 0.5));
    const int loose = static_cast<int>(std::lround(girth / 1.0));
    const int bays = static_cast<int>(std::lround(kBoxLength / s.frameSpacing));
    expectTrue("the two spacings really do give different band counts", tight == 2 * loose);
    expectEqual("the tighter spacing sets the seams for the whole strake",
                static_cast<long long>(mesh.panels.size()), 2LL * bays * tight);
    // And the plating still tiles the hull exactly, whichever spacing won.
    expectNear("with the shell still exactly covered", mesh.plateArea(),
               kBoxLength * kBoxBeam + 2.0 * kBoxLength * kBoxDepth,
               1e-6 * kBoxLength * kBoxBeam);
}

// Each strake carries the longitudinal on its upper seam. Shifting them all one
// band down leaves the count and the total length untouched -- so weight says
// nothing -- but puts one stiffener on the centreline keel, where the centre
// girder goes, and leaves the seam below the deck edge bare.
void testLongitudinalsSitOnStrakeSeamsNotOnTheCentreline() {
    const TriMesh hull = boxHull();
    Scantlings s = bareBoxScantlings();
    s.shell[0].stiffened = true;
    s.shell[0].longitudinal = flatBar(0.200, 0.010);
    const StructuralMesh mesh = makeStructuralMesh(hull, s);

    double nearestToKeel = 1e300, highest = -1e300;
    int count = 0;
    for (const StructuralMember& m : mesh.members) {
        if (m.role != MemberRole::Longitudinal) continue;
        ++count;
        nearestToKeel = std::min(nearestToKeel, std::abs(m.a.y) + m.a.z);
        highest = std::max(highest, m.a.z);
    }
    expectTrue("there are longitudinals to look at", count > 20);
    // The bands are 1 m of girth, so the lowest seam is a metre outboard of the
    // centreline and the highest is a metre below the deck edge.
    expectTrue("none sits on the centreline keel", nearestToKeel > 0.5);
    expectTrue("and the seam below the deck edge carries one",
               highest > kBoxDepth - 1.5 && highest < kBoxDepth - 0.5);
}

// A stiffener's web grows *into* the hull. Getting the direction wrong on one
// side is invisible to weight and almost invisible to the hull girder -- a web
// lying horizontally contributes the same second moment whichever way it points
// -- but it is the eccentricity the FEM will attach its beam elements with, so
// it is checked directly.
void testStiffenerWebsPointIntoTheHull() {
    const StructuralMesh mesh = makeStructuralMesh(game::buildFerry().hull, ferryScantlings());

    double worstOutboard = 0, worstDownward = 0;
    int checked = 0;
    for (const StructuralMember& m : mesh.members) {
        if (m.role != MemberRole::Longitudinal && m.role != MemberRole::Frame) continue;
        ++checked;
        // A web pointing outboard has its rise on the same side as the plating it
        // is welded to.
        worstOutboard = std::max(worstOutboard, m.a.y * m.rise.y);
        worstDownward = std::min(worstDownward, m.rise.z);
    }
    expectTrue("there are shell stiffeners to check", checked > 5000);
    expectTrue("no web points out through the shell", worstOutboard < 1e-9);
    expectTrue("and none points down through the bottom", worstDownward > -1e-9);

    // Guard: the check would be vacuous if every rise were vertical, which it is
    // on the flat of bottom. The bilge and the side must contribute webs that
    // really do lean inboard.
    int leaning = 0;
    for (const StructuralMember& m : mesh.members)
        if (m.role == MemberRole::Longitudinal && std::abs(m.rise.y) > 0.5) ++leaning;
    expectTrue("and most webs are not simply vertical", leaning > 1000);
}

// --- Hull girder ---------------------------------------------------------------

// A box hull with a deck: bottom plate, two side plates and a deck plate, all
// rectangles at known heights. The neutral axis and second moment are one
// evening's arithmetic and there is exactly one right answer.
void testHullGirderSectionAgainstHandCalculation() {
    const TriMesh hull = boxHull();
    const double shellThickness = 0.020, deckThickness = 0.010;

    Scantlings s = bareBoxScantlings(shellThickness);
    Deck deck;
    deck.name = "weather";
    deck.z = kBoxDepth;
    deck.xFrom = -0.5 * kBoxLength;
    deck.xTo = 0.5 * kBoxLength;
    deck.thickness = deckThickness;
    deck.stiffened = false;
    deck.beamed = false;
    s.decks = {deck};

    const StructuralMesh mesh = makeStructuralMesh(hull, s);
    // 2.5 m is deliberately *between* frame stations, which are every 5 m.
    const HullGirderSection got = hullGirderSection(mesh, 2.5);

    // By hand, about the baseline. Panel corners are on the plate mid-surface, so
    // the bottom's material centre is at z = 0 and the deck's at z = D.
    const double bottomArea = kBoxBeam * shellThickness;
    const double sideArea = kBoxDepth * shellThickness;   // one side
    const double deckArea = kBoxBeam * deckThickness;
    const double area = bottomArea + 2.0 * sideArea + deckArea;
    const double neutral =
        (bottomArea * 0.0 + 2.0 * sideArea * 0.5 * kBoxDepth + deckArea * kBoxDepth) / area;
    const double second =
        kBoxBeam * std::pow(shellThickness, 3) / 12.0 + bottomArea * neutral * neutral +
        2.0 * (shellThickness * std::pow(kBoxDepth, 3) / 12.0 +
               sideArea * std::pow(0.5 * kBoxDepth - neutral, 2)) +
        kBoxBeam * std::pow(deckThickness, 3) / 12.0 +
        deckArea * std::pow(kBoxDepth - neutral, 2);

    expectNear("hull girder area is the sum of four plates", got.area, area, 1e-6 * area);
    expectNear("the neutral axis is where hand calculation puts it", got.neutralAxis, neutral,
               1e-6 * neutral);
    expectNear("and so is the second moment", got.secondMoment, second, 1e-5 * second);
    expectNear("extreme fibre at the deck", got.zDeck, kBoxDepth + 0.5 * deckThickness, 1e-4);
    expectNear("extreme fibre at the keel", got.zKeel, -0.5 * shellThickness, 1e-4);
    expectNear("section modulus to the deck", got.modulusDeck, second / (got.zDeck - neutral),
               1e-5 * got.modulusDeck);
    expectNear("section modulus to the keel", got.modulusKeel, second / (neutral - got.zKeel),
               1e-5 * got.modulusKeel);

    // Guards. The bottom is 20 mm over the full beam and the deck 10 mm, so the
    // neutral axis sits well below mid-depth -- a section that came out symmetric,
    // with the two moduli equal, is exactly the degenerate case a careless version
    // of this test would accept.
    expectTrue("the section is not symmetric about mid-depth", neutral < 0.45 * kBoxDepth);
    expectTrue("so the two section moduli differ", got.modulusKeel > 1.3 * got.modulusDeck);
    expectTrue("and the deck is really carrying some of it", deckArea > 0.15 * area);
}

// A stiffener contributes its *strong* second moment to the hull girder when its
// web stands vertically and its *weak* one when the web lies flat, and which of
// the two applies is decided by the direction the web rises. Two deep girders --
// one on the keel with its web upright, one on the side with its web horizontal
// -- put a measurable amount of the section's second moment into that choice, so
// swapping the two axes is caught rather than lost in the Steiner terms.
void testHullGirderTakesTheAxisTheWebActuallyPresents() {
    const TriMesh hull = boxHull();
    const double plate = 0.006;
    Scantlings s = bareBoxScantlings(plate);

    Girder centre;
    centre.name = "centre";
    centre.y = 0.0;
    centre.z = 0.0;
    centre.rise = {0, 0, 1};   // web upright: strong axis about the horizontal
    centre.profile = flatBar(2.000, 0.040);
    centre.attachedPlateThickness = plate;
    Girder side;
    side.name = "side_stringer";
    side.y = 0.5 * kBoxBeam;
    side.z = 0.5 * kBoxDepth;
    side.bothSides = true;
    side.rise = {0, -1, 0};    // web flat: weak axis about the horizontal
    side.profile = flatBar(1.500, 0.030);
    side.attachedPlateThickness = plate;
    s.girders = {centre, side};

    const StructuralMesh mesh = makeStructuralMesh(hull, s);
    const HullGirderSection got = hullGirderSection(mesh, 2.5);

    const ProfileSection centreSection = profileSection(centre.profile);
    const ProfileSection sideSection = profileSection(side.profile);

    const double bottomArea = kBoxBeam * plate;
    const double sideArea = kBoxDepth * plate;
    const double centreHeight = 0.5 * plate + centreSection.centroid;
    const double area = bottomArea + 2.0 * sideArea + centreSection.area + 2.0 * sideSection.area;
    const double neutral = (2.0 * sideArea * 0.5 * kBoxDepth + centreSection.area * centreHeight +
                            2.0 * sideSection.area * 0.5 * kBoxDepth) /
                           area;
    const double second =
        kBoxBeam * std::pow(plate, 3) / 12.0 + bottomArea * neutral * neutral +
        2.0 * (plate * std::pow(kBoxDepth, 3) / 12.0 +
               sideArea * std::pow(0.5 * kBoxDepth - neutral, 2)) +
        // Web upright: the strong axis.
        centreSection.secondMoment +
        centreSection.area * std::pow(centreHeight - neutral, 2) +
        // Webs flat: the weak one.
        2.0 * (sideSection.secondMomentWeak +
               sideSection.area * std::pow(0.5 * kBoxDepth - neutral, 2));

    expectNear("area with girders", got.area, area, 1e-6 * area);
    expectNear("neutral axis with girders", got.neutralAxis, neutral, 1e-6 * neutral);
    expectNear("second moment takes each web's own axis", got.secondMoment, second, 1e-5 * second);

    // Guards. The two own-axis terms have to be a real share of the section, or
    // swapping them would be undetectable here as well.
    expectTrue("the upright web's own second moment matters",
               centreSection.secondMoment > 0.005 * second);
    expectTrue("and the flat webs' does not, which is the point",
               2.0 * sideSection.secondMomentWeak < 1e-4 * second);
    expectTrue("the two axes differ by orders",
               sideSection.secondMoment > 1000.0 * sideSection.secondMomentWeak);
}

// A cut that lands on a frame station sits on the seam between two bays. Counting
// both doubles the section -- and leaves the neutral axis, which is a ratio,
// looking exactly right. The box hand calculation is what caught this; the check
// below is what keeps it caught, without needing the hand arithmetic.
void testSectionAtAFrameStationIsNotDoubleCounted() {
    const TriMesh hull = boxHull();
    Scantlings s = bareBoxScantlings();
    s.shell[0].stiffened = true;
    s.shell[0].longitudinal = flatBar(0.200, 0.010);
    const StructuralMesh mesh = makeStructuralMesh(hull, s);

    // 0 is a frame station; 0.01 either side of it is not.
    expectTrue("the cut really is on a station",
               std::any_of(mesh.frameStations.begin(), mesh.frameStations.end(),
                           [](double x) { return std::abs(x) < 1e-9; }));

    const HullGirderSection on = hullGirderSection(mesh, 0.0);
    const HullGirderSection before = hullGirderSection(mesh, -0.01);
    const HullGirderSection after = hullGirderSection(mesh, 0.01);

    expectTrue("the section is non-trivial", on.area > 0.1);
    expectNear("a cut on a frame station matches the bay just aft of it", on.area, before.area,
               1e-9 * on.area);
    expectNear("and the bay just forward of it", on.area, after.area, 1e-9 * on.area);
    expectNear("second moment likewise", on.secondMoment, after.secondMoment,
               1e-9 * on.secondMoment);

    // The forward-most station has no bay ahead of it and must still report a
    // section rather than nothing at all.
    const HullGirderSection stem = hullGirderSection(mesh, mesh.frameStations.back());
    expectNear("the forward end still has a section", stem.area, on.area, 1e-6 * on.area);
}

// One cut, with the two element populations kept apart. A total on its own cannot
// tell "the plating vanished" from "the taper moved": the defect these tests were
// written for took the ferry's section at x = 19.2 down to 23.8% of its neighbours,
// which is not a wrong number so much as exactly the stiffeners' share of the
// section with all the plating gone.
struct SectionCut {
    double area = 0, neutralAxis = 0, secondMoment = 0;
    double plateArea = 0;
    int plate = 0, stiffener = 0;
};

SectionCut cutSection(const StructuralMesh& mesh, double x) {
    SectionCut c;
    const HullGirderSection g = hullGirderSection(mesh, x);
    c.area = g.area;
    c.neutralAxis = g.neutralAxis;
    c.secondMoment = g.secondMoment;
    for (const SectionElement& e : sectionElements(mesh, x)) {
        if (e.stiffener) {
            ++c.stiffener;
        } else {
            ++c.plate;
            c.plateArea += e.area;
        }
    }
    return c;
}

// The same frame station, named four ways, must give the same section.
//
// The ferry is 120 m over 50 bays, and `-60 + 120*33/50` comes out at
// 19.200000000000003 -- one unit in the last place above the 19.2 a drawing says.
// A cut asked for at 19.2 therefore sits 3.6e-15 m *aft* of the bay that owns the
// seam. The membership test in `sectionElements` tolerates that by 1e-9 and admits
// the bay; the crossing search that follows it was exact, found no panel edge
// changing sign across the plane, and dropped all 188 plate panels. All 181
// stiffeners survived, because the member branch clamps its interpolation
// parameter into the member and so was never exposed to the same question. What
// came back was 23.8% of the area with a neutral axis -- a ratio -- still looking
// perfectly correct.
//
// Eleven of the ferry's 51 stations were affected and forty were not, and which is
// which is decided by nothing more than the direction the division rounds. So the
// station is asked for every way a caller might write it.
void testTheSectionDoesNotDependOnHowAStationIsSpelled() {
    const StructuralMesh mesh = makeStructuralMesh(game::buildFerry().hull, ferryScantlings());

    int spelledDifferently = 0;
    for (double station : mesh.frameStations) {
        // A micrometre is finer than any station is ever placed, so this is the
        // number a drawing would carry.
        const double drawing = std::round(station * 1e6) / 1e6;
        if (drawing != station) ++spelledDifferently;

        const SectionCut want = cutSection(mesh, station);
        expectTrue("every station has both plating and stiffening",
                   want.plate > 60 && want.stiffener > 60);

        for (double x : {drawing, std::nextafter(station, -1e300),
                         std::nextafter(station, 1e300)}) {
            const SectionCut got = cutSection(mesh, x);
            expectEqual("the same station names the same plate panels", got.plate, want.plate);
            expectEqual("and the same stiffeners", got.stiffener, want.stiffener);
            // The spellings differ by at most 7e-15 m and the section's own taper is
            // under 0.02 m2 per metre along this hull, so the answers can only
            // differ by round-off. Measured, the worst disagreement over all 51
            // stations is 1.3e-15 relative; most are bit-identical, because a cut
            // within the seam tolerance is taken at the seam itself.
            expectNear("and the same area", got.area, want.area, 1e-13 * want.area);
            expectNear("the same second moment", got.secondMoment, want.secondMoment,
                       1e-13 * want.secondMoment);
            expectNear("and the same neutral axis", got.neutralAxis, want.neutralAxis,
                       1e-13 * want.neutralAxis);
        }
    }

    // Otherwise vacuous: if every station's stored value already were the decimal a
    // drawing gives, the loop above would have asked the same question twice. 23 of
    // the ferry's 51 do not land on their decimal.
    expectTrue("the stored stations and a drawing's decimals really do differ",
               spelledDifferently >= 5);
}

// A sweep, not a point.
//
// A single cut cannot say whether a station is wrong, because there is nothing to
// compare it with; the ferry's midship area is 1.80 m2 and so was the value at
// every station either side of the broken one. Walking x finely along the whole
// length and requiring the section to move no faster than the hull's own taper is
// what turns "this number looks plausible" into a test.
void testTheSectionSweepIsContinuousAlongTheLength() {
    const Scantlings description = ferryScantlings();
    const StructuralMesh mesh = makeStructuralMesh(game::buildFerry().hull, description);

    // Every 100 mm, plus every frame station both as stored and as a drawing writes
    // it, because that pair is exactly where the seam arithmetic bites.
    std::vector<double> xs;
    const double xLo = mesh.frameStations.front(), xHi = mesh.frameStations.back();
    const int steps = 1200;
    for (int i = 0; i <= steps; ++i) xs.push_back(xLo + (xHi - xLo) * i / steps);
    for (double station : mesh.frameStations) {
        xs.push_back(station);
        xs.push_back(std::round(station * 1e6) / 1e6);
    }
    std::sort(xs.begin(), xs.end());
    xs.erase(std::unique(xs.begin(), xs.end()), xs.end());

    // Lipschitz bounds on the section, measured over this hull's own bays with the
    // element population held fixed: 0.0198 m2/m, 1.00 m4/m and 0.125 m/m. They are
    // properties of the ferry rather than converging quantities, so they carry
    // enough room for a re-tessellation and no more. The failure they exist to catch
    // moved the area by 1.37 m2 across a 30 mm step -- a slope of 45 m2/m.
    constexpr double kAreaSlope = 0.03;    // m2 per m
    constexpr double kSecondSlope = 2.0;   // m4 per m
    constexpr double kAxisSlope = 0.2;     // m per m

    SectionCut previous{};
    double previousX = 0;
    double worstAreaSlope = 0;
    for (std::size_t i = 0; i < xs.size(); ++i) {
        const SectionCut c = cutSection(mesh, xs[i]);

        // Count, do not just total. A section with stiffeners and no plating at all
        // is what the defect produced, and it is not a ship's section: plating is
        // three quarters of the area of every station on this hull, the thinnest of
        // them at the stem being 71.3%.
        expectTrue("every cut has plating", c.plate > 0);
        expectTrue("every cut has stiffening", c.stiffener > 0);
        expectTrue("and the plating carries most of the section",
                   c.area > 0 && c.plateArea > 0.65 * c.area);

        if (i > 0 && previous.plate == c.plate && previous.stiffener == c.stiffener) {
            const double h = xs[i] - previousX;
            if (h > 0) {
                worstAreaSlope = std::max(worstAreaSlope, std::abs(c.area - previous.area) / h);
                expectTrue("area moves no faster than the taper",
                           std::abs(c.area - previous.area) <= kAreaSlope * h);
                expectTrue("nor does the second moment",
                           std::abs(c.secondMoment - previous.secondMoment) <= kSecondSlope * h);
                expectTrue("nor the neutral axis",
                           std::abs(c.neutralAxis - previous.neutralAxis) <= kAxisSlope * h);
            }
        }
        previous = c;
        previousX = xs[i];
    }

    // The bound has to be doing work: if the section were constant along the whole
    // length the loop above would pass on any implementation that returned the same
    // answer everywhere, including a wrong one.
    expectTrue("the hull really does taper", worstAreaSlope > 1e-3);
    expectTrue("the sweep is fine enough to straddle every station",
               xs.size() > 4 * mesh.frameStations.size());

    // --- The parallel middle body, where the population may not change at all ------
    //
    // Nothing in the description starts or stops inside |x| < 34 -- the innermost
    // boundary is the weather deck's forward end and the wing bulkhead's -- so over
    // a range well inside that, every station must present the same elements as
    // midship. This is the flat form of the assertion above and it is what the
    // defect tripped over: at x = 19.2 the plate count went from 188 to 0.
    constexpr double kParallel = 26.4;   // m either side of midship: 23 frame stations
    for (const Deck& d : description.decks)
        expectTrue("no deck begins or ends inside the parallel body",
                   std::abs(d.xFrom) > kParallel && std::abs(d.xTo) > kParallel);
    for (const Girder& g : description.girders)
        expectTrue("no girder does either",
                   std::abs(g.xFrom) > kParallel && std::abs(g.xTo) > kParallel);
    for (const ShellRegion& r : description.shell)
        expectTrue("nor a plating region",
                   std::abs(r.xFrom) > kParallel && std::abs(r.xTo) > kParallel);
    // Bulkheads are not in that list on purpose: the ferry has transverse ones at
    // -8 m and 20 m, and a transverse plate has no extent along x, so the hull
    // girder never sees it however close to midship it stands.

    const SectionCut midship = cutSection(mesh, 0.0);
    expectTrue("the midship section is a full one", midship.plate > 150 && midship.stiffener > 150);

    int stationsInside = 0, sampled = 0;
    for (double station : mesh.frameStations)
        if (std::abs(station) <= kParallel) ++stationsInside;
    for (double x : xs) {
        if (std::abs(x) > kParallel) continue;
        ++sampled;
        const SectionCut c = cutSection(mesh, x);
        expectEqual("the parallel body presents the same plate panels everywhere", c.plate,
                    midship.plate);
        expectEqual("and the same stiffeners", c.stiffener, midship.stiffener);
    }
    expectTrue("the parallel body really does contain a run of stations", stationsInside >= 20);
    expectTrue("and the sweep samples it densely", sampled >= 400);
}

// Exactly once, across a chain: the property the half-open rule exists for.
//
// `sectionElements` gives a cut landing on a frame seam to the bay *forward* of it,
// so that the two bays either side do not both build it. Asserted directly here by
// making the bays tell themselves apart -- one plating thickness per bay, reported
// back on every element -- so a cut that took two bays, or none, is visible rather
// than inferred from a total. Both mistakes have happened: counting both doubled
// the section, and counting neither dropped every plate panel.
void testEachBayServesExactlyOneFrameStation() {
    const TriMesh hull = boxHull();

    // The stations have to be known before the description can name them, so the
    // mesh is built twice: once to lay out the frames, once with a plating band per
    // bay. The bands are widened by a micron so the girth stays covered at the
    // seams; a band is chosen by the bay's midpoint, which is unambiguous.
    const std::vector<double> stations = makeStructuralMesh(hull, bareBoxScantlings()).frameStations;
    const std::size_t bays = stations.size() - 1;

    Scantlings banded = bareBoxScantlings();
    banded.shell.clear();
    for (std::size_t k = 0; k < bays; ++k) {
        ShellRegion r;
        r.name = "bay_" + std::to_string(k);
        r.xFrom = stations[k] - 1e-6;
        r.xTo = stations[k + 1] + 1e-6;
        r.thickness = 0.010 + 0.001 * static_cast<double>(k);
        r.stiffened = false;
        banded.shell.push_back(r);
    }

    std::vector<std::string> problems;
    const StructuralMesh mesh = makeStructuralMesh(hull, banded, &problems);
    expectTrue("a plating band per bay builds without complaint", problems.empty());
    expectTrue("there are enough bays for the chain to mean anything", bays >= 10);

    // The box's girth is a closed form: the flat of bottom plus both sides.
    const double girth = kBoxBeam + 2.0 * kBoxDepth;

    // Which bay served the cut, read off the elements rather than assumed.
    const auto served = [&](double x) {
        std::vector<double> thickness;
        double area = 0;
        for (const SectionElement& e : sectionElements(mesh, x)) {
            area += e.area;
            if (std::find_if(thickness.begin(), thickness.end(), [&](double t) {
                    return std::abs(t - e.thickness) < 1e-12;
                }) == thickness.end())
                thickness.push_back(e.thickness);
        }
        return std::pair<std::vector<double>, double>{thickness, area};
    };

    for (std::size_t k = 0; k < stations.size(); ++k) {
        // The forward-most station has no bay ahead of it, so the cut interval
        // closes there and the last bay serves it.
        const std::size_t owner = std::min(k, bays - 1);
        const double want = banded.shell[owner].thickness;

        // Every way of landing on the seam: on it, a hair aft of it but inside the
        // seam tolerance, and clear of it on the forward side.
        std::vector<double> spellings{stations[k], stations[k] - 1e-12, stations[k] + 1e-6};
        if (k + 1 == stations.size()) spellings.pop_back();   // past the stem
        for (double x : spellings) {
            const auto [thickness, area] = served(x);
            expectEqual("a cut on a seam is served by exactly one bay",
                        static_cast<long long>(thickness.size()), 1);
            if (thickness.size() != 1) continue;
            expectNear("and it is the bay forward of the seam", thickness[0], want, 1e-12);
            // The area is the closed form, so neither a doubled bay nor a dropped
            // one can hide inside it. The 1e-7 is the ray nudge `sampleSection`
            // uses to keep off the transom and the deck cap.
            expectNear("with the whole girth plated exactly once", area, girth * want,
                       3e-7 * girth * want);
        }

        // And a hair aft of the seam, outside the tolerance, must be the bay astern
        // -- otherwise "exactly one" would be satisfied by one bay owning everything.
        if (k > 0) {
            const double astern = banded.shell[k - 1].thickness;
            const auto [thickness, area] = served(stations[k] - 1e-6);
            expectEqual("just aft of a seam it is the other bay",
                        static_cast<long long>(thickness.size()), 1);
            if (thickness.size() == 1)
                expectNear("namely the one astern", thickness[0], astern, 1e-12);
            expectNear("which plates the whole girth in its turn", area, girth * astern,
                       3e-7 * girth * astern);
        }
    }

    // Vacuous unless the bays are actually distinguishable.
    expectTrue("the bays really do differ",
               banded.shell.front().thickness < 0.5 * banded.shell.back().thickness);
}

// A strip of plating narrow enough to look like a rounding error is still plating.
//
// The cut through a panel is rejected below a length, because a chord of zero
// across a girth band the hull has closed to nothing is not a section element. That
// guard sits between "numerically degenerate" and "thin", and where it is put is a
// decision: raised to a centimetre it silently eats a gunwale strake and moves
// nothing else enough to notice. `SectionElement::width` is the field that says so,
// and it is checked here because nothing else in the engine reads it -- a written
// and never-read field is exactly where a wrong value survives.
void testANarrowStrakeStillReachesTheSection() {
    Scantlings s = bareBoxScantlings();
    ShellRegion main = s.shell[0];
    main.girthTo = 0.99975;              // 19.995 m of the box's 20 m half-girth
    ShellRegion gunwale = s.shell[0];
    gunwale.name = "gunwale";
    gunwale.girthFrom = 0.99975;         // and 5 mm at the deck edge
    gunwale.thickness = 0.030;
    s.shell = {main, gunwale};

    std::vector<std::string> problems;
    const StructuralMesh mesh = makeStructuralMesh(boxHull(), s, &problems);
    expectTrue("a hull with a gunwale strake builds without complaint", problems.empty());

    double plateWidth = 0, plateArea = 0, widthTimesThickness = 0;
    double narrowest = 1e300, widest = 0;
    int narrow = 0, plate = 0, stiffenerWithWidth = 0;
    for (const SectionElement& e : sectionElements(mesh, 2.5)) {
        if (e.stiffener) {
            if (e.width != 0.0) ++stiffenerWithWidth;
            continue;
        }
        ++plate;
        plateWidth += e.width;
        plateArea += e.area;
        widthTimesThickness += e.width * e.thickness;
        narrowest = std::min(narrowest, e.width);
        widest = std::max(widest, e.width);
        if (e.width < 0.01) ++narrow;
    }

    // The contract of the field: a plating strip's area *is* its girth times its
    // thickness, and a stiffener spans no girth at all. Exact, not approximate --
    // both come from the same two numbers.
    expectEqual("a stiffener spans no girth", stiffenerWithWidth, 0);
    expectNear("a strip's area is its girth times its thickness", widthTimesThickness, plateArea,
               1e-12 * plateArea);

    // The gunwale strake is 5 mm of a 40 m girth: one strip each side, and no other
    // strip anywhere near that narrow.
    expectEqual("the gunwale strake reaches the section, both sides", narrow, 2);
    expectNear("at the width it was asked for", narrowest, 0.005, 1e-9);
    expectTrue("and it really is the odd one out", widest > 100.0 * narrowest);
    expectTrue("the rest of the girth is plated too", plate > 20);

    // The box's girth is B + 2D, and the strips must tile it: never long, because
    // that would mean a strip counted twice, and short only by what one strip per
    // side loses chording across the right angle at the bilge. A band of `w` across
    // a square corner is at worst w/sqrt(2) of it, so the whole girth cannot fall
    // below B + 2D - 2*w*(1 - 1/sqrt(2)); it is measured at 5 mm short of 40 m,
    // because the knuckle happens to land near a band seam rather than mid-band.
    const double girth = kBoxBeam + 2.0 * kBoxDepth;
    const double chordLoss = 2.0 * widest * (1.0 - 1.0 / std::sqrt(2.0));
    expectTrue("no girth is plated twice", plateWidth <= girth + 1e-9);
    expectTrue("and none is left bare beyond the chord across the bilge",
               plateWidth >= girth - chordLoss);
}

// Transverse members carry no longitudinal stress, so nothing athwartships may
// appear in the hull girder. The routine decides that from geometry rather than
// from a label, which is checked here by adding a great deal of transverse steel
// and requiring the section not to notice.
void testTransverseStructureStaysOutOfTheHullGirder() {
    const TriMesh hull = boxHull();
    Scantlings bare = bareBoxScantlings();
    const HullGirderSection without = hullGirderSection(makeStructuralMesh(hull, bare), 2.5);

    Scantlings heavy = bare;
    heavy.framed = true;
    heavy.frameProfile = flatBar(1.200, 0.030);   // absurdly deep, on purpose
    Bulkhead bulkhead;
    bulkhead.name = "transverse";
    bulkhead.transverse = true;
    bulkhead.position = 0.0;
    bulkhead.zFrom = 0.0;
    bulkhead.zTo = kBoxDepth;
    bulkhead.thickness = 0.040;
    bulkhead.stiffened = true;
    bulkhead.stiffener = flatBar(0.300, 0.020);
    heavy.bulkheads = {bulkhead};

    const StructuralMesh mesh = makeStructuralMesh(hull, heavy);
    const HullGirderSection with = hullGirderSection(mesh, 2.5);

    expectTrue("the transverse steel was actually built",
               mesh.memberCount(MemberRole::Frame) > 100 &&
                   mesh.panelCount(PanelRole::Bulkhead) > 10);
    expectTrue("and it weighs a great deal", mesh.steelMass() > 1.2 * makeStructuralMesh(hull, bare).steelMass());
    expectNear("but none of it works in the hull girder", with.area, without.area,
               1e-9 * without.area);
    expectNear("nor in its second moment", with.secondMoment, without.secondMoment,
               1e-9 * without.secondMoment);

    // The cut taken *through* the bulkhead must not pick it up either.
    const HullGirderSection atBulkhead = hullGirderSection(mesh, 0.0);
    expectNear("a cut in the plane of a bulkhead still ignores it", atBulkhead.area, without.area,
               1e-9 * without.area);

    // Nor at the forward end, which is the one station where the cut interval has
    // to close rather than stay half-open -- and so the one place a zero-extent
    // member could slip in.
    Scantlings stemmed = heavy;
    Bulkhead collision = bulkhead;
    collision.name = "fore_peak";
    collision.position = 0.5 * kBoxLength;
    stemmed.bulkheads.push_back(collision);
    const StructuralMesh withStem = makeStructuralMesh(hull, stemmed);
    const HullGirderSection foreEnd = hullGirderSection(withStem, 0.5 * kBoxLength);
    expectTrue("the fore end still reports a section", foreEnd.area > 0.1);
    expectNear("and the stem frame and fore peak bulkhead stay out of it", foreEnd.area,
               without.area, 1e-9 * without.area);

    // And the control: longitudinal steel must change the section, or the test
    // above would pass on a routine that ignored everything.
    Scantlings stiffened = bare;
    stiffened.shell[0].stiffened = true;
    stiffened.shell[0].longitudinal = flatBar(0.200, 0.010);
    const HullGirderSection longitudinally =
        hullGirderSection(makeStructuralMesh(hull, stiffened), 2.5);
    expectTrue("longitudinal steel does change the section",
               longitudinally.area > 1.05 * without.area);
    expectTrue("and its second moment", longitudinally.secondMoment > 1.05 * without.secondMoment);
}

// Section properties must *converge* under refinement, not merely be plausible.
// The cut is taken in the tapering forward body, where a coarse bay chords across
// real change in the section; amidships the hull is parallel and the answer is
// exact at any spacing, which would make this test say nothing.
void testSectionPropertiesConvergeUnderFrameRefinement() {
    const TriMesh hull = game::buildFerry().hull;
    const double cut = 36.5;   // forward of the parallel middle body, off any station

    std::vector<double> second, neutral;
    for (double spacing : {9.6, 4.8, 2.4, 1.2, 0.6}) {
        Scantlings s = ferryScantlings();
        s.frameSpacing = spacing;
        const HullGirderSection g = hullGirderSection(makeStructuralMesh(hull, s), cut);
        expectTrue("every refinement produces a section", g.area > 0.5);
        second.push_back(g.secondMoment);
        neutral.push_back(g.neutralAxis);
    }

    // Successive differences must shrink, and the last pair must agree far more
    // closely than the first.
    const std::size_t last = second.size() - 1;
    const double coarse = std::abs(second[0] - second[last]);
    const double fine = std::abs(second[last - 1] - second[last]);
    expectTrue("the coarse mesh is measurably wrong", coarse / second[last] > 5e-3);
    expectTrue("refinement converges rather than wandering", fine < 0.02 * coarse);
    expectTrue("and the neutral axis converges with it",
               std::abs(neutral[last - 1] - neutral[last]) < 1e-4);

    int improved = 0;
    for (std::size_t i = 1; i < second.size(); ++i)
        if (std::abs(second[i] - second[last]) <= std::abs(second[i - 1] - second[last]) + 1e-12)
            ++improved;
    expectEqual("every refinement gets closer", improved, static_cast<int>(second.size()) - 1);
}

// --- The reference ferry --------------------------------------------------------

// Whether a set of scantlings is a ship or a guess has one published answer: the
// IACS minimum midship section modulus. A structure that cannot reach it would
// not be classed, whatever it weighed.
void testFerryMidshipSectionMeetsTheRuleMinimum() {
    const TriMesh hull = game::buildFerry().hull;
    const Scantlings s = ferryScantlings();
    std::vector<std::string> problems;
    const StructuralMesh mesh = makeStructuralMesh(hull, s, &problems);
    for (const std::string& p : problems) std::printf("     ferry scantlings: %s\n", p.c_str());
    expectTrue("the reference scantlings build without complaint", problems.empty());

    const HullGirderSection midship = hullGirderSection(mesh, 0.0);
    const double minimum = ruleMinimumSectionModulus(120.0, 20.0, 0.66);
    std::printf("     ferry midship: A = %.3f m2, NA = %.3f m, I = %.2f m4, "
                "Z deck %.2f m3, Z keel %.2f m3 (rule minimum %.2f m3)\n",
                midship.area, midship.neutralAxis, midship.secondMoment, midship.modulusDeck,
                midship.modulusKeel, minimum);

    expectTrue("the rule minimum is a real number", minimum > 3.0 && minimum < 3.5);
    expectTrue("section modulus to the deck clears the rule minimum",
               midship.modulusDeck > minimum);
    expectTrue("and to the keel", midship.modulusKeel > minimum);
    // Two to three times the minimum would mean the arrangement is not a ship's.
    expectTrue("without being absurdly over-built", midship.modulusDeck < 3.0 * minimum);

    // The bottom carries more steel than the deck, so the neutral axis sits below
    // mid-depth. A section that came out symmetric would mean the girders, the
    // double bottom and the strake thicknesses had all failed to land.
    expectTrue("the neutral axis sits below mid-depth", midship.neutralAxis < 7.2);
    expectTrue("but not on the bottom", midship.neutralAxis > 5.0);
    expectTrue("the section area is a 120 m ship's", midship.area > 1.2 && midship.area < 3.0);
    expectTrue("and the extreme fibres are the hull's own",
               midship.zDeck > 14.9 && midship.zKeel < 0.01);
}

// Steel weight, and whether it is a credible hull for a 120 m ro-pax. It is
// checked against the ship's own displacement rather than a remembered figure,
// and the bound is deliberately loose: this is plating and primary stiffening
// only, so it *should* come in under a real hull steel weight.
void testFerrySteelWeightIsCredible() {
    const Ship ferry = game::buildFerry();
    const StructuralMesh mesh = makeStructuralMesh(ferry.hull, ferryScantlings());

    const double displacement = integrateBelowPlane(ferry.hull, {0, 0, 1}, 5.5).volume * 1.025;
    const double steel = mesh.steelMass() / 1000.0;
    std::printf("     ferry structure: %zu panels, %zu members, plating %.0f t, "
                "stiffening %.0f t, total %.0f t against %.0f t displacement at the design draft\n",
                mesh.panels.size(), mesh.members.size(), mesh.plateMass() / 1000.0,
                mesh.memberMass() / 1000.0, steel, displacement);

    // Independent of the mesh accessors: re-derive the mass from the panels and
    // members directly, which is what the accessors are supposed to be doing.
    double byHand = 0;
    for (const PlatePanel& p : mesh.panels)
        byHand += p.area() * p.thickness * mesh.materials[static_cast<std::size_t>(p.material)].density;
    for (const StructuralMember& m : mesh.members)
        byHand += m.length() * profileSection(m.profile).area *
                  mesh.materials[static_cast<std::size_t>(m.material)].density;
    expectNear("steelMass() is the sum of its parts", mesh.steelMass(), byHand, 1e-9 * byHand);

    // Both sides in tonnes. Plating and primary stiffening only, so this *should*
    // land under a real hull steel weight, which also carries brackets, floors,
    // pillars, foundations, the stem and stern frames and every piece of minor
    // structure -- another 15 to 30 percent.
    expectTrue("the structure weighs a credible fraction of the ship",
               steel > 0.12 * displacement && steel < 0.35 * displacement);
    // Real hulls run 0.4 to 0.7 tonnes of stiffening per tonne of plating. Outside
    // that band the arrangement is either unstiffened plate or a forest of bars.
    const double ratio = mesh.memberMass() / mesh.plateMass();
    expectTrue("stiffening to plating is in the shipbuilding band", ratio > 0.35 && ratio < 0.75);
    // Mean plating thickness, as a sanity check on the strake thicknesses.
    const double meanThickness = mesh.plateMass() / (mesh.plateArea() * 7850.0);
    expectTrue("mean plating is between 8 and 16 mm",
               meanThickness > 0.008 && meanThickness < 0.016);
}

// The structural mesh is derived from the hull *surface*, not from the hull
// mesh's vertices. Re-tessellating the hull -- which changes every triangle in
// it -- must therefore leave the structure essentially unchanged, and the frame
// stations must be the ship's, not the tessellation's.
void testStructureIsIndependentOfTheHullTessellation() {
    const TriMesh coarse = boxHull(5, 3);
    const TriMesh fine = boxHull(37, 19);
    expectTrue("the two tessellations really do differ",
               fine.tris.size() > 8 * coarse.tris.size());

    Scantlings s = bareBoxScantlings();
    s.framed = true;
    s.frameProfile = flatBar(0.300, 0.012);
    s.shell[0].stiffened = true;
    s.shell[0].longitudinal = flatBar(0.200, 0.010);

    const StructuralMesh a = makeStructuralMesh(coarse, s);
    const StructuralMesh b = makeStructuralMesh(fine, s);
    expectEqual("the same number of panels either way",
                static_cast<long long>(a.panels.size()), static_cast<long long>(b.panels.size()));
    expectNear("the same steel weight", b.steelMass(), a.steelMass(), 1e-6 * a.steelMass());
    expectNear("the same frame spacing", b.frameSpacing, a.frameSpacing, 1e-12);

    // And the frame stations are the ones the scantlings asked for, not the ones
    // the hull mesh happened to be drawn at.
    expectNear("frame spacing is the ship's, not the mesh's", a.frameSpacing, s.frameSpacing,
               1e-9);
    expectEqual("frame count follows the spacing", static_cast<long long>(a.frameStations.size()),
                static_cast<long long>(std::lround(kBoxLength / s.frameSpacing)) + 1);

    // A tapered hull is the harder case: there the two tessellations describe
    // slightly different surfaces, so the structure may differ a little -- but
    // only by the tessellation's own error, not by an order.
    Scantlings ferry = ferryScantlings();
    const StructuralMesh c = makeStructuralMesh(game::buildFerry().hull, ferry);
    expectTrue("the ferry structure is built", c.panels.size() > 5000);
    expectEqual("with a frame at every station the spacing implies",
                static_cast<long long>(c.frameStations.size()),
                static_cast<long long>(std::lround(120.0 / ferry.frameSpacing)) + 1);
}

// --- Refusals --------------------------------------------------------------------

void testUnreasonableScantlingsAreReported() {
    expectTrue("the reference arrangement raises nothing",
               validateScantlings(ferryScantlings()).empty());

    Scantlings gap = ferryScantlings();
    // Punch a hole in the girth coverage: the side strake now starts above where
    // the bilge ends. Nothing about the thicknesses looks wrong, and the gap is
    // *not* everywhere -- the thinner end regions still cover that girth forward
    // of 44 m and aft of -44 m, so a check that ignored x would see full coverage.
    gap.shell[2].girthFrom = 0.40;
    const std::vector<std::string> reported = validateScantlings(gap);
    expectTrue("a gap in the girth coverage is reported", !reported.empty());
    expectTrue("even though the ends still cover that girth",
               gap.shell[4].girthFrom < 0.35 && gap.shell[5].girthFrom < 0.35);

    // And the generator leaves that band unplated rather than guessing at it.
    std::vector<std::string> problems;
    const StructuralMesh mesh = makeStructuralMesh(game::buildFerry().hull, gap, &problems);
    expectTrue("the generator says so too", !problems.empty());
    expectTrue("and the shell really is short of plating",
               shellPanelArea(mesh) < 0.95 * shellPanelArea(makeStructuralMesh(
                                          game::buildFerry().hull, ferryScantlings())));

    Scantlings bad;
    bad.frameSpacing = 0.0;
    bad.longitudinalSpacing = -1.0;
    bad.materials.clear();
    expectTrue("a description with no spacing and no material is refused",
               validateScantlings(bad).size() >= 3);

    Scantlings dangling = ferryScantlings();
    dangling.decks[0].material = 7;
    expectTrue("a deck naming a material that does not exist is refused",
               !validateScantlings(dangling).empty());

    // An empty hull must not crash or invent structure.
    std::vector<std::string> empty;
    const StructuralMesh nothing = makeStructuralMesh(TriMesh{}, ferryScantlings(), &empty);
    expectTrue("an empty hull is reported", !empty.empty());
    expectEqual("and produces no panels", static_cast<long long>(nothing.panels.size()), 0);
    expectNear("and no steel", nothing.steelMass(), 0.0, 1e-12);
}

// The rule formula has a published shape, and its own edges are worth pinning:
// it must rise with length and beam, and it must not reward a fine hull with a
// lower requirement than the floor allows.
void testRuleMinimumSectionModulus() {
    const double a = ruleMinimumSectionModulus(120.0, 20.0, 0.66);
    // C = 10.75 - ((300 - 120)/100)^1.5 = 10.75 - 1.8^1.5.
    const double c = 10.75 - std::pow(1.8, 1.5);
    expectNear("the rule formula is C L^2 B (Cb + 0.7) x 1e-6", a,
               c * 120.0 * 120.0 * 20.0 * (0.66 + 0.7) * 1e-6, 1e-12);
    expectTrue("a longer ship needs more", ruleMinimumSectionModulus(200.0, 20.0, 0.66) > 3.0 * a);
    expectTrue("a beamier ship needs more", ruleMinimumSectionModulus(120.0, 30.0, 0.66) > a);
    expectNear("the block coefficient has a floor of 0.60",
               ruleMinimumSectionModulus(120.0, 20.0, 0.45),
               ruleMinimumSectionModulus(120.0, 20.0, 0.60), 1e-12);
    expectTrue("a fuller ship needs more than the floor",
               ruleMinimumSectionModulus(120.0, 20.0, 0.80) >
                   ruleMinimumSectionModulus(120.0, 20.0, 0.60));
}

}  // namespace

void runScantlingTests() {
    std::printf("\n--- scantlings and structural mesh ---\n");
    testFlatBarProfileIsTheSchoolbookRectangle();
    testAngleAndTeeDifferAboutTheAxisAcrossTheWeb();
    testStiffenedSectionAgainstHandIntegration();
    testSmearingLosesTheStiffenersLeverArm();
    testStiffenerEccentricityHasASign();
    testPanelAreaOfATrapezoidAndAWarpedQuad();
    testBoxPanelsTileExactly();
    testBoxSteelWeightIsTheClosedForm();
    testFerryPanelsApproachTheHullAreaFromBelow();
    testPanelCornersConvergeOntoTheHull();
    testFrameSpacingIsRoundedToAWholeNumberOfBays();
    testALaterShellRegionOverridesAnEarlierOne();
    testDeckPlatingIsClippedToTheHull();
    testBulkheadsAreClippedToTheHull();
    testTheTightestSpacingWinsAStrake();
    testLongitudinalsSitOnStrakeSeamsNotOnTheCentreline();
    testStiffenerWebsPointIntoTheHull();
    testHullGirderSectionAgainstHandCalculation();
    testHullGirderTakesTheAxisTheWebActuallyPresents();
    testSectionAtAFrameStationIsNotDoubleCounted();
    testTheSectionDoesNotDependOnHowAStationIsSpelled();
    testTheSectionSweepIsContinuousAlongTheLength();
    testEachBayServesExactlyOneFrameStation();
    testANarrowStrakeStillReachesTheSection();
    testTransverseStructureStaysOutOfTheHullGirder();
    testSectionPropertiesConvergeUnderFrameRefinement();
    testFerryMidshipSectionMeetsTheRuleMinimum();
    testFerrySteelWeightIsCredible();
    testStructureIsIndependentOfTheHullTessellation();
    testUnreasonableScantlingsAreReported();
    testRuleMinimumSectionModulus();
}
