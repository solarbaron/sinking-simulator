// SPDX-License-Identifier: MIT
//
// Validation of the geometric and hydrostatic core against closed-form answers.
// Everything downstream -- flooding rates, stability, capsize -- is only as good
// as these integrals, so they get checked against algebra rather than eyeballed.
#include "engine/core/geometry.hpp"
#include "engine/sim/ship.hpp"
#include "harness.hpp"

#include <cstdio>
#include <string>

using namespace sim;
using testing::expectNear;
using testing::expectTrue;

namespace {

// --- Closed-mesh integration -----------------------------------------------

void testBoxIntegrals() {
    const TriMesh box = makeBox({-2, -3, -4}, {6, 5, 4});  // 8 x 8 x 8
    const VolumeIntegral v = integrate(box);
    expectNear("box volume", v.volume, 512.0, 1e-9);
    expectNear("box centroid x", v.centroid.x, 2.0, 1e-9);
    expectNear("box centroid y", v.centroid.y, 1.0, 1e-9);
    expectNear("box centroid z", v.centroid.z, 0.0, 1e-9);
}

void testAxisAlignedClip() {
    const TriMesh box = makeBox({0, 0, 0}, {10, 4, 6});  // 240 m^3
    // Cut horizontally at z = 1.5: a 10 x 4 x 1.5 slab.
    const VolumeIntegral lower = integrateBelowPlane(box, {0, 0, 1}, 1.5);
    expectNear("axis-aligned clip volume", lower.volume, 60.0, 1e-9);
    expectNear("axis-aligned clip centroid z", lower.centroid.z, 0.75, 1e-9);
    expectNear("axis-aligned clip centroid x", lower.centroid.x, 5.0, 1e-9);

    // Degenerate ends must be exact, not merely close.
    expectNear("clip below everything", integrateBelowPlane(box, {0, 0, 1}, -1).volume, 0.0, 1e-12);
    expectNear("clip above everything", integrateBelowPlane(box, {0, 0, 1}, 99).volume, 240.0, 1e-9);
}

void testTiltedClipAgainstAlgebra() {
    // A unit cube cut by a plane through the origin with normal (1,1,1)/sqrt(3)
    // keeps exactly one corner tetrahedron of the cube: volume 1/6.
    const TriMesh cube = makeBox({0, 0, 0}, {1, 1, 1});
    const Vec3 n = normalize(Vec3{1, 1, 1});
    const VolumeIntegral cut = integrateBelowPlane(cube, n, dot(n, Vec3{1, 0, 0}));
    expectNear("tilted clip: corner tetrahedron volume", cut.volume, 1.0 / 6.0, 1e-12);
    // Centroid of that tetrahedron is the average of its four vertices.
    expectNear("tilted clip: centroid x", cut.centroid.x, 0.25, 1e-12);
}

void testVolumeSolveRoundTrip() {
    const TriMesh box = makeBox({-30, -8, 0}, {30, 8, 10});
    const Vec3 n = normalize(Vec3{0.12, -0.25, 1.0});  // an arbitrary heel and trim
    for (double target : {50.0, 900.0, 4800.0, 9000.0}) {
        const double off = solvePlaneOffsetForVolume(box, n, target, -1e30, 1e30);
        const double got = integrateBelowPlane(box, n, off).volume;
        expectNear("volume solve round trip @ " + std::to_string(target), got, target,
                   1e-6 * target);
    }
}

// --- Constructive solid geometry --------------------------------------------

// A hull shaped like the real thing: fine at the ends, full amidships, with a
// turn of the bilge. Non-convex enough to be a fair test of the clipper.
TriMesh testHull() {
    const std::vector<double> waterlines{0.0, 1.0, 2.0, 3.2, 5.0, 7.0, 10.0};
    std::vector<Station> stations;
    for (int i = 0; i <= 20; ++i) {
        Station s;
        s.x = -50.0 + 100.0 * i / 20.0;
        const double u = std::abs(s.x) / 50.0;
        const double fx = u < 0.6 ? 1.0 : 1.0 - 0.8 * std::pow((u - 0.6) / 0.4, 2.0);
        for (double z : waterlines)
            s.halfBeam.push_back(9.0 * fx * (z >= 3.2 ? 1.0 : 0.5 + 0.5 * z / 3.2));
        stations.push_back(s);
    }
    return makeHullFromStations(stations, waterlines);
}

// Watertightness is the precondition for every integral in the engine, and it is
// silently violated rather than loudly: a mesh with one strip wound backwards
// still produces plausible-looking numbers. Check the generators directly.
void testMeshesAreWatertight() {
    expectTrue("makeBox produces a closed manifold", isClosedManifold(makeBox({0, 0, 0}, {2, 3, 4})));
    expectTrue("makeHullFromStations produces a closed manifold", isClosedManifold(testHull()));

    const TriMesh hull = testHull();
    expectTrue("clipByPlane output is a closed manifold",
               isClosedManifold(clipByPlane(hull, {0, 0, 1}, 4.0)));
    expectTrue("clipToBox output is a closed manifold",
               isClosedManifold(clipToBox(hull, {-10, -3, 1}, {15, 3, 6})));

    // The integral is reference-point independent only on a closed mesh, so
    // disagreement between the two routines is the symptom to watch for.
    const double whole = integrate(hull).volume;
    for (double off : {40.0, 80.0, 200.0})
        expectNear("integrateBelowPlane above the mesh equals the whole volume",
                   integrateBelowPlane(hull, {0, 0, 1}, off).volume, whole, 1e-9 * whole);
}

// The capped solid produced by clipByPlane() must enclose exactly the volume the
// cap-free integrator reports for the same half-space. If the cap leaks, or is
// wound backwards, or a loop is dropped, these disagree immediately.
void testClipMatchesIntegral() {
    const TriMesh hull = testHull();
    const Vec3 normals[] = {{1, 0, 0}, {0, 0, 1}, {0, 1, 0},
                            normalize(Vec3{0.3, -0.8, 0.5})};
    const double offsets[] = {-20.0, 0.0, 3.7, 12.0};

    for (const Vec3& n : normals)
        for (double off : offsets) {
            const double want = integrateBelowPlane(hull, n, off).volume;
            const TriMesh cut = clipByPlane(hull, n, off);
            const VolumeIntegral got = integrate(cut);
            expectNear("clip volume matches integral", got.volume, want,
                       std::max(1e-6, 1e-6 * std::abs(want)));
            if (want > 1.0) {
                const Vec3 wantC = integrateBelowPlane(hull, n, off).centroid;
                expectNear("clip centroid matches integral (z)", got.centroid.z, wantC.z, 1e-4);
            }
        }
}

// Compartments carved out of the hull must tile it: cut the same hull into a grid
// of boxes and the volumes have to add back up to the whole, with nothing lost in
// the cracks and nothing counted twice.
void testSubdivisionTiles() {
    const TriMesh hull = testHull();
    const double whole = integrate(hull).volume;

    const double xs[] = {-60, -30, -5, 20, 60};
    const double ys[] = {-20, -4, 0, 4, 20};
    const double zs[] = {-1, 2.5, 6.0, 20};

    double sum = 0;
    int nonEmpty = 0;
    for (int i = 0; i + 1 < 5; ++i)
        for (int j = 0; j + 1 < 5; ++j)
            for (int k = 0; k + 1 < 4; ++k) {
                const TriMesh cell =
                    clipToBox(hull, {xs[i], ys[j], zs[k]}, {xs[i + 1], ys[j + 1], zs[k + 1]});
                if (cell.tris.empty()) continue;
                const double v = integrate(cell).volume;
                expectTrue("subdivision cell has non-negative volume", v > -1e-6);
                sum += v;
                ++nonEmpty;
            }

    expectTrue("subdivision produced cells", nonEmpty > 20);
    expectNear("subdivision volumes sum to the hull volume", sum, whole, 1e-4 * whole);
}

// A compartment carved from the hull must actually be inside the hull -- the
// failure the hand-authored boxes had.
void testClippedCompartmentStaysInsideHull() {
    const TriMesh hull = testHull();
    // A box far wider than the hull at the bow, where the hull is finest.
    const TriMesh naive = makeBox({35, -9, 0}, {50, 9, 7});
    const TriMesh carved = clipToBox(hull, {35, -9, 0}, {50, 9, 7});

    const double naiveV = integrate(naive).volume;
    const double carvedV = integrate(carved).volume;
    expectTrue("naive box overstates a fine-ended compartment", naiveV > carvedV * 1.3);

    // Every vertex of the carved compartment must lie within the hull's own
    // half-breadth at its station.
    bool inside = true;
    for (const Vec3& p : carved.verts) {
        const double u = std::abs(p.x) / 50.0;
        const double fx = u < 0.6 ? 1.0 : 1.0 - 0.8 * std::pow((u - 0.6) / 0.4, 2.0);
        const double halfBeam = 9.0 * fx * (p.z >= 3.2 ? 1.0 : 0.5 + 0.5 * std::max(p.z, 0.0) / 3.2);
        if (std::abs(p.y) > halfBeam + 0.05) inside = false;
    }
    expectTrue("carved compartment stays within the hull envelope", inside);
}

// --- Hydrostatics -----------------------------------------------------------

// A homogeneous box floats at a draft of (rho_body / rho_water) * height.
// This is Archimedes' principle and the engine has no excuse to miss it.
void testArchimedes() {
    Ship s;
    s.hull = makeBox({-25, -6, 0}, {25, 6, 8});
    s.deckEdgeZ = 8.0;
    const double expectedDraft = 3.0;
    const double volume = 50.0 * 12.0 * expectedDraft;
    s.lightshipMass = volume * kRhoSeawater;
    s.lightshipCog = {0, 0, 4.0};
    s.gyradii = {4.0, 14.0, 14.0};
    s.initialise(0.0);

    // Let any residual transient settle.
    for (int i = 0; i < 20000; ++i) s.step(0.005, 0.0);

    const Diagnostics d = s.diagnostics(0.0);
    expectNear("box barge floats at Archimedean draft", d.draftMidship, expectedDraft, 0.01);
    expectNear("box barge stays upright (heel)", d.heelDeg, 0.0, 0.05);
    expectNear("box barge stays upright (trim)", d.trimDeg, 0.0, 0.05);

    // Metacentric radius of a rectangular waterplane: BM = I/V = L*B^3/12 / V.
    const double bm = (50.0 * 12.0 * 12.0 * 12.0 / 12.0) / volume;
    const double kb = expectedDraft / 2.0;
    expectNear("box barge GM matches KB + BM - KG", d.gmTransverse,
               kb + bm - s.lightshipCog.z, 0.02);
}

// Free surface effect: the same mass of water, loose in a wide tank instead of
// bolted down as solid ballast, must reduce GM by rho*i/displacement.
void testFreeSurfaceEffect() {
    auto makeBarge = [](bool liquid) {
        Ship s;
        s.hull = makeBox({-25, -6, 0}, {25, 6, 8});
        s.deckEdgeZ = 8.0;
        s.compartments = {[&] {
            Compartment c;
            c.name = "tank";
            c.mesh = makeBox({-20, -5, 0.0}, {20, 5, 4.0});
            c.permeability = 1.0;
            c.ventedToAtmosphere = true;
            c.waterVolume = liquid ? 200.0 : 0.0;  // 200 m^3, a shallow layer
            return c;
        }()};
        // Keep total displacement identical between the two cases by moving the
        // equivalent mass into the lightship when the tank is dry.
        const double waterMass = 200.0 * kRhoSeawater;
        s.lightshipMass = 50.0 * 12.0 * 3.0 * kRhoSeawater - (liquid ? waterMass : 0.0);
        // Solid ballast sits at the same height the loose water would.
        const double solidZ = 0.5;
        s.lightshipCog = {0, 0, liquid ? 4.0
                                       : (4.0 * (50.0 * 12.0 * 3.0 * kRhoSeawater - waterMass)
                                          + solidZ * waterMass)
                                             / (50.0 * 12.0 * 3.0 * kRhoSeawater)};
        if (!liquid) s.lightshipMass = 50.0 * 12.0 * 3.0 * kRhoSeawater;
        s.gyradii = {4.0, 14.0, 14.0};
        s.initialise(0.0);
        return s;
    };

    const Ship loose = makeBarge(true);
    const Ship solid = makeBarge(false);
    const double gmLoose = loose.diagnostics(0.0).gmTransverse;
    const double gmSolid = solid.diagnostics(0.0).gmTransverse;

    // i = l*b^3/12 for the 40 x 10 m tank surface; the correction is rho*i/Delta.
    const double i = 40.0 * 10.0 * 10.0 * 10.0 / 12.0;
    const double displacementVolume = 50.0 * 12.0 * 3.0;
    const double expectedLoss = i / displacementVolume;

    expectTrue("free surface reduces GM", gmLoose < gmSolid);
    expectNear("free surface loss matches rho*i/displacement", gmSolid - gmLoose,
               expectedLoss, 0.15 * expectedLoss);
}

// Air trapped in a sealed compartment must arrest flooding once its pressure
// balances the outside head -- the reason capsized hulls stay up for hours.
void testTrappedAirArrestsFlooding() {
    Ship s;
    s.hull = makeBox({-25, -6, 0}, {25, 6, 8});
    s.deckEdgeZ = 8.0;

    Compartment sealed;
    sealed.name = "sealed_void";
    sealed.mesh = makeBox({-10, -4, 0.0}, {10, 4, 3.0});
    sealed.permeability = 1.0;
    s.compartments = {sealed};

    Opening hole;
    hole.name = "hole";
    hole.a = kSea;
    hole.b = 0;
    hole.pos = {0, -4, 0.2};
    hole.area = 0.5;
    s.openings = {hole};  // no vent anywhere: the air has nowhere to go

    s.lightshipMass = 50.0 * 12.0 * 3.0 * kRhoSeawater;
    s.lightshipCog = {0, 0, 4.0};
    s.gyradii = {4.0, 14.0, 14.0};
    s.initialise(0.0);

    for (int i = 0; i < 120000; ++i) s.step(0.005, 0.0);

    const Compartment& c = s.compartments[0];
    expectTrue("sealed compartment took some water", c.fillFraction() > 0.05);
    expectTrue("trapped air stopped it filling", c.fillFraction() < 0.95);
    expectTrue("trapped air is above atmospheric", c.airPressure > kPatm * 1.02);

    // Boyle's law check: the air was compressed from the full compartment volume
    // to the remaining void, isothermally.
    const double p0V0 = kPatm * c.grossVolume;
    expectNear("isothermal compression conserves pV", c.airPressure * c.airVolume(),
               p0V0, 0.02 * p0V0);
}

// Water must not appear or vanish: the sum over compartments has to equal what
// crossed the hull boundary.
void testMassConservation() {
    Ship s;
    s.hull = makeBox({-25, -6, 0}, {25, 6, 8});
    s.deckEdgeZ = 8.0;

    auto vented = [](const char* name, Vec3 lo, Vec3 hi) {
        Compartment c;
        c.name = name;
        c.mesh = makeBox(lo, hi);
        c.permeability = 1.0;
        c.ventedToAtmosphere = true;
        return c;
    };
    s.compartments = {vented("a", {-20, -5, 0}, {0, 5, 4}),
                      vented("b", {0, -5, 0}, {20, 5, 4})};

    Opening breach;
    breach.name = "breach";
    breach.a = kSea;
    breach.b = 0;
    breach.pos = {-10, -5, 0.5};
    breach.area = 0.3;

    Opening door;
    door.name = "door";
    door.a = 0;
    door.b = 1;
    door.pos = {0, 0, 0.1};
    door.area = 1.0;
    s.openings = {breach, door};

    s.lightshipMass = 50.0 * 12.0 * 3.0 * kRhoSeawater;
    s.lightshipCog = {0, 0, 4.0};
    s.gyradii = {4.0, 14.0, 14.0};
    s.initialise(0.0);

    double crossedBoundary = 0.0;
    const double dt = 0.005;
    for (int i = 0; i < 40000; ++i) {
        s.step(dt, 0.0);
        for (const Opening& o : s.openings)
            if (o.name == "breach" && o.lastFlowWasWater) crossedBoundary += o.lastFlow * dt;
    }

    double held = 0.0;
    for (const Compartment& c : s.compartments) held += c.waterVolume;
    expectNear("water in compartments equals water through the breach", held, crossedBoundary,
               1e-3 * std::max(held, 1.0));

    // While the breach is admitting water the two spaces do *not* sit level: an
    // orifice only passes flow when there is a head difference across it, so the
    // damaged side runs a few centimetres high. Shut the breach and the offset
    // must decay to nothing.
    const double headWhileFlooding =
        std::abs(s.compartments[0].surfaceWorldZ - s.compartments[1].surfaceWorldZ);
    expectTrue("a head difference drives flow through the door", headWhileFlooding > 0.01);

    for (Opening& o : s.openings)
        if (o.name == "breach") o.open = false;
    for (int i = 0; i < 40000; ++i) s.step(dt, 0.0);

    expectTrue("compartments level off once inflow stops",
               std::abs(s.compartments[0].surfaceWorldZ - s.compartments[1].surfaceWorldZ) < 0.01);

    double heldAfter = 0.0;
    for (const Compartment& c : s.compartments) heldAfter += c.waterVolume;
    expectNear("no water created or destroyed while settling", heldAfter, held, 1e-6 * held);
}

}  // namespace

void runCoreTests() {
    std::printf("\n--- geometry, hydrostatics and flooding ---\n");
    testBoxIntegrals();
    testAxisAlignedClip();
    testTiltedClipAgainstAlgebra();
    testVolumeSolveRoundTrip();
    testMeshesAreWatertight();
    testClipMatchesIntegral();
    testSubdivisionTiles();
    testClippedCompartmentStaysInsideHull();
    testArchimedes();
    testFreeSurfaceEffect();
    testTrappedAirArrestsFlooding();
    testMassConservation();
}
