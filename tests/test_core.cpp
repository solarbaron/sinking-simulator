// SPDX-License-Identifier: MIT
//
// Validation of the geometric and hydrostatic core against closed-form answers.
// Everything downstream -- flooding rates, stability, capsize -- is only as good
// as these integrals, so they get checked against algebra rather than eyeballed.
#include "engine/core/geometry.hpp"
#include "engine/sim/ship.hpp"
#include "engine/sim/waves.hpp"
#include "game/prototype/ferry.hpp"
#include "harness.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <string>
#include <vector>

using namespace sim;
using testing::expectEqual;
using testing::expectNear;
using testing::expectTrue;

namespace {

// The adiabatic limit: gas that gives its heat to nothing at all. A compartment
// with a fire in it is nearer this than it is to the isothermal default over the
// seconds that decide where the smoke goes.
constexpr double kInfinity = std::numeric_limits<double>::infinity();

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
// The buoyancy force and the reported diagnostics use *different* implementations
// of one integral, and nothing compared them.
//
// `integrateBelowPlane` and `PlaneSweep::below` carry the same contract, word for
// word — "volume and centroid of the region satisfying dot(n, x) <= offset" — and
// `ship.cpp` splits between them: the force path (`integrateRigidBody`, the
// floodwater centroid, the cached waterplane area, the roll/pitch stiffness) goes
// through `PlaneSweep`, while `diagnostics` and `rightingArmAtHeel` call the free
// function. No test constructed a `PlaneSweep` at all.
//
// They agree to the bit on a unit normal, which is what makes the assertion below
// tight rather than tolerant. **They do not agree otherwise**, and that is the
// point of the second half: `PlaneSweep` normalises in its constructor and
// `integrateBelowPlane` does not, so `offset` means "distance along the unit
// normal" to one and "dot(n,x) with the caller's own n" to the other. Neither
// header states a precondition. At |n| = 0.5 they return answers 1.7x apart; at
// |n| = 10 the free function reports eighteen times the hull's entire volume,
// with no error and no assert.
//
// No live caller is broken — every one passes a unit normal or a rotation of one.
// This pins the agreement where it is relied on and states the divergence where it
// is not, so the next caller to pass an unnormalised vector finds out here.
void testTheTwoPlaneIntegralsAgree() {
    const TriMesh hull = testHull();
    const Vec3 normals[] = {{0, 0, 1}, {1, 0, 0},
                            normalize(Vec3{0.3, -0.8, 0.5}),
                            normalize(Vec3{0, std::sin(0.6), std::cos(0.6)})};
    const double offsets[] = {-20.0, 0.0, 3.7, 5.5, 12.0, 40.0};

    double worstVolume = 0, worstCentroid = 0, largest = 0;
    for (const Vec3& n : normals) {
        const PlaneSweep sweep(hull, n);
        for (double off : offsets) {
            const VolumeIntegral loose = integrateBelowPlane(hull, n, off);
            const VolumeIntegral cached = sweep.below(off);
            worstVolume = std::max(worstVolume, std::abs(loose.volume - cached.volume));
            largest = std::max(largest, std::abs(loose.volume));
            if (loose.volume > 1e-9)
                worstCentroid = std::max(worstCentroid,
                                         length(loose.centroid - cached.centroid));
        }
    }
    std::printf("     the two plane integrals: worst volume %.3e m3, worst centroid %.3e m\n",
                worstVolume, worstCentroid);

    // Bit-identical on a unit normal but for the ULP `normalize` introduces, so
    // this is asserted at what was measured rather than at a comfortable band.
    expectTrue("the cached sweep and the free integral agree on volume",
               worstVolume <= 1e-14 * std::max(largest, 1.0));
    expectTrue("and on centroid", worstCentroid <= 1e-12);
    // The guard: a hull that clipped to nothing would satisfy both trivially.
    expectTrue("and there was a real volume to compare", largest > 1e3);

    // The divergence, stated rather than merely known. `PlaneSweep` normalises and
    // `integrateBelowPlane` does not, so on a non-unit normal the same arguments
    // describe different planes — and the free function's reference point
    // `o = n * offset` stops lying on `dot(n,x) == offset`, which is the whole
    // reason its cap tetrahedra are supposed to vanish.
    const Vec3 half{0, 0, 0.5};
    const VolumeIntegral looseHalf = integrateBelowPlane(hull, half, 2.5);
    const VolumeIntegral sweepHalf = PlaneSweep(hull, half).below(2.5);
    expectTrue("a non-unit normal makes the two disagree, which is why callers"
               " must pass a unit one",
               std::abs(looseHalf.volume - sweepHalf.volume) > 0.1 * sweepHalf.volume);
}

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

// --- Buoyancy under a wavy free surface --------------------------------------

// Flat water is the degenerate case of a height field, so the general routine
// must reproduce the plane routine exactly. If it does not, one of them is wrong
// and every wave result would inherit the error invisibly.
void testWavySurfaceReducesToThePlaneCase() {
    const TriMesh hull = testHull();
    for (double level : {-1.0, 0.0, 2.5, 4.0, 6.5, 9.0}) {
        const VolumeIntegral plane = integrateBelowPlane(hull, {0, 0, 1}, level);
        const VolumeIntegral surface =
            integrateBelowSurface(hull, [level](double, double) { return level; });
        expectNear("flat height field matches the plane volume at " + std::to_string(level),
                   surface.volume, plane.volume, std::max(1e-9, 1e-9 * plane.volume));
        if (plane.volume > 1.0) {
            expectNear("flat height field matches the plane centroid z", surface.centroid.z,
                       plane.centroid.z, 1e-7);
            expectNear("flat height field matches the plane centroid x", surface.centroid.x,
                       plane.centroid.x, 1e-6);
        }
    }
}

// A box under a sinusoidal surface has an analytic submerged volume, which is
// what makes this a real check rather than a comparison of two approximations.
//
// For a box spanning x in [x0, x1] with its bottom below the trough and its top
// above the crest, the submerged volume is W * integral of (h(x) - zBottom) dx,
// and for h = A cos(k x) that integral is closed form.
//
// **cos, not sin.** The first version of this test used a sine over a symmetric
// range, where the wave's contribution integrates to exactly zero -- so the
// analytic answer did not depend on the wave at all, the routine returned it to
// machine precision at every resolution, and the convergence assertion was
// comparing zero against zero. A cosine over the same range does not cancel.
void testSinusoidalSurfaceAgainstAnalyticVolume() {
    const double x0 = -10.0, x1 = 10.0, y1 = 3.0;
    const double bottom = -5.0, top = 5.0;
    const double amplitude = 2.0;
    const double wavelength = 40.0;
    const double k = 2.0 * kPi / wavelength;

    auto wave = [&](double x, double) { return amplitude * std::cos(k * x); };

    const double width = 2.0 * y1;
    const double integral = (amplitude / k * std::sin(k * x1) - bottom * x1) -
                            (amplitude / k * std::sin(k * x0) - bottom * x0);
    const double analytic = width * integral;

    // Guard against the vacuity that caught the first version: the wave must
    // actually move the answer away from the still-water volume.
    const double stillWater = width * (-bottom) * (x1 - x0);
    expectTrue("the analytic answer genuinely depends on the wave",
               std::abs(analytic - stillWater) > 0.1 * stillWater);

    double previousError = 1e30;
    double ratio = 0;
    for (int divisions : {8, 16, 32, 64}) {
        std::vector<Station> stations;
        const std::vector<double> waterlines{bottom, top};
        for (int i = 0; i <= divisions; ++i) {
            Station station;
            station.x = x0 + (x1 - x0) * i / divisions;
            station.halfBeam = {y1, y1};
            stations.push_back(station);
        }
        const TriMesh box = makeHullFromStations(stations, waterlines);

        const VolumeIntegral got = integrateBelowSurface(box, wave);
        const double error = std::abs(got.volume - analytic) / analytic;
        expectTrue("sinusoidal volume converges, divisions=" + std::to_string(divisions),
                   error < previousError);
        if (previousError < 1e29) ratio = previousError / error;
        previousError = error;
    }
    expectTrue("the finest sinusoidal mesh is within 1e-8 of the analytic volume",
               previousError < 1e-8);
    // Measured fourth order: the edge-midpoint rule is exact for quadratics and
    // its leading error term is O(h^4) for a smooth surface. Asserting the *rate*
    // catches a scheme that converges to the right answer for the wrong reason.
    expectTrue("convergence is fourth order, as the quadrature rule implies",
               ratio > 10.0 && ratio < 24.0);
}

// A wave must actually change the answer, or the routine is silently ignoring
// its height field and every check above would still pass.
void testWavySurfaceIsNotSecretlyFlat() {
    const TriMesh hull = testHull();
    const VolumeIntegral flat =
        integrateBelowSurface(hull, [](double, double) { return 5.0; });
    const VolumeIntegral crest =
        integrateBelowSurface(hull, [](double x, double) { return 5.0 + 2.0 * std::cos(x * 0.1); });

    expectTrue("a wave changes the displaced volume", std::abs(crest.volume - flat.volume) > 1.0);
    // A wave crest amidships and troughs at the ends shifts buoyancy toward the
    // middle, but the total is bounded by the flat-water volumes at crest and
    // trough level.
    const double atCrest = integrateBelowPlane(hull, {0, 0, 1}, 7.0).volume;
    const double atTrough = integrateBelowPlane(hull, {0, 0, 1}, 3.0).volume;
    expectTrue("wave displacement lies between the trough and crest still-water values",
               crest.volume > atTrough && crest.volume < atCrest);

    // And a wave whose crest is entirely above the hull must submerge it wholly.
    const VolumeIntegral drowned =
        integrateBelowSurface(hull, [](double, double) { return 100.0; });
    expectNear("a surface above the whole hull submerges all of it", drowned.volume,
               integrate(hull).volume, 1e-6 * integrate(hull).volume);
}

// The height field is evaluated once per vertex, not once per triangle corner.
// That is a performance property with no effect on any answer, so nothing else
// in this file would notice if a refactor undid it -- and it is worth about a
// factor of two on a wavy tick, because with a 128-component spectrum the
// surface query *is* the tick.
//
// Two configurations pin it exactly rather than approximately: with the surface
// below the whole hull every triangle takes the dry early-out, so the corner
// form asks 3*tris and the vertex form asks exactly verts; with the surface
// above it, both add three quadrature points per triangle on top.
void testSurfaceIsQueriedOncePerVertex() {
    const TriMesh hull = testHull();
    const auto verts = static_cast<long long>(hull.verts.size());
    const auto tris = static_cast<long long>(hull.tris.size());
    expectTrue("the mesh has more corners than vertices, or there is nothing to save",
               3 * tris > verts);

    long long dryQueries = 0;
    const VolumeIntegral dry =
        integrateBelowSurface(hull, [&](double, double) { ++dryQueries; return -100.0; });
    expectEqual("a hull entirely above the surface is queried once per vertex", dryQueries, verts);
    expectEqual("and displaces nothing", static_cast<long long>(dry.volume), 0LL);

    // Submerged, every triangle is wet -- but not every triangle is quadratured.
    // The integrand F = (0, 0, z - h) points straight up, so a panel whose normal
    // is horizontal carries no flux and is dropped before the quadrature runs.
    // On a ship that is not a corner case: the side plating is vertical, and here
    // it is 45% of the hull.
    long long upright = 0;
    for (const Tri& t : hull.tris)
        if (0.5 * cross(hull.verts[t.b] - hull.verts[t.a],
                        hull.verts[t.c] - hull.verts[t.a]).z != 0.0)
            ++upright;
    expectTrue("the hull really does have vertical plating to skip", upright < tris * 3 / 4);

    long long wetQueries = 0;
    const VolumeIntegral wet =
        integrateBelowSurface(hull, [&](double, double) { ++wetQueries; return 100.0; });
    expectEqual("a submerged hull adds three quadrature points per sloping triangle, no more",
                wetQueries, verts + 3 * upright);
    expectNear("and displaces its whole volume", wet.volume, integrate(hull).volume,
               1e-6 * integrate(hull).volume);

    // Evaluating eagerly must not let vertices no triangle uses into the answer.
    TriMesh padded = hull;
    for (int i = 0; i < 32; ++i) padded.verts.push_back({1e4 + i, -1e4, -1e4});
    const auto wave = [](double x, double) { return 5.0 + 2.0 * std::cos(x * 0.1); };
    expectNear("unreferenced vertices do not perturb the integral",
               integrateBelowSurface(padded, wave).volume,
               integrateBelowSurface(hull, wave).volume, 1e-9);
}

// A station short of a half-breadth column used to be padded with zeros, which
// turns the missing levels into a knife edge. The mesh still closes, still passes
// the manifold check, and still integrates -- to a displacement that is simply
// wrong. That is the exact shape of every defect in CLAUDE.md's table: a wrong
// answer that survives every check aimed at something else.
void testShortStationsAreRepairedAndReported() {
    const std::vector<double> waterlines{0.0, 2.0, 4.0, 6.0};
    const auto build = [&](std::size_t columnsOnOne, std::vector<std::string>* problems) {
        std::vector<Station> stations;
        for (int i = 0; i <= 8; ++i) {
            Station s;
            s.x = -20.0 + 40.0 * i / 8.0;
            s.halfBeam = {5.0, 5.0, 5.0, 5.0};
            if (i == 4) s.halfBeam.resize(columnsOnOne);
            stations.push_back(s);
        }
        return makeHullFromStations(stations, waterlines, problems);
    };

    std::vector<std::string> problems;
    const TriMesh whole = build(4, &problems);
    expectTrue("a complete station list raises nothing", problems.empty());

    problems.clear();
    const TriMesh short2 = build(2, &problems);
    expectTrue("a short station is reported", problems.size() == 1);

    // The dangerous part, asserted directly: the damaged mesh is still closed, so
    // isClosedManifold() -- which caught a 40% displacement error once -- cannot
    // see this one at all.
    expectTrue("and the damaged hull is still a closed manifold", isClosedManifold(short2));

    // With the repair carrying the last breadth upward, a station missing its top
    // columns is wall-sided there, so the volume is unchanged for this prismatic
    // hull. Zero padding would have cut a notch out of it.
    expectNear("carrying the last breadth upward preserves a prismatic hull's volume",
               integrate(short2).volume, integrate(whole).volume, 1e-9 * integrate(whole).volume);

    // Guard the guard: the hull has to be big enough that a notch would have
    // shown, or "unchanged" says nothing.
    expectTrue("the hull has a volume worth comparing", integrate(whole).volume > 1000.0);

    // An empty station is still zero-width, and still reported rather than
    // producing a silently smaller ship.
    problems.clear();
    const TriMesh empty = build(0, &problems);
    expectTrue("a station with no half-breadths at all is reported", problems.size() == 1);
    expectTrue("and it does change the hull", integrate(empty).volume < integrate(whole).volume);
}

// --- Ships in waves ----------------------------------------------------------

namespace {

// A plain box barge, so every hydrostatic quantity has a closed form -- but
// tessellated along its length rather than built from makeBox().
//
// This is not cosmetic. makeBox() gives two triangles per face, and a single
// 60 m panel under a wave shorter than itself picks up a *phase-dependent*
// buoyancy error of around 6% of displacement, sign and all. Under-tessellation
// does not merely blur the answer, it invents or destroys displacement, and the
// error would ride silently through every seakeeping result. See
// testHullMustResolveTheWavelength below.
Ship makeBarge(int lengthDivisions = 24) {
    Ship s;
    std::vector<Station> stations;
    for (int i = 0; i <= lengthDivisions; ++i) {
        Station station;
        station.x = -30.0 + 60.0 * i / lengthDivisions;
        station.halfBeam = {8.0, 8.0};
        stations.push_back(station);
    }
    s.hull = makeHullFromStations(stations, {0.0, 12.0});
    s.deckEdgeZ = 12.0;
    s.lightshipMass = 60.0 * 16.0 * 4.0 * kRhoSeawater;  // floats at 4 m
    s.lightshipCog = {0, 0, 5.0};
    s.gyradii = {5.0, 16.0, 16.0};
    return s;
}

// A single-component sea of known amplitude and wavelength. Phase is seeded, so
// the caller finds a crest by scanning time rather than by assuming one.
SeaState oneWave(double waveHeight, double wavelength) {
    SeaState state;
    state.significantHeight = waveHeight;
    // Deep-water dispersion: lambda = g T^2 / (2 pi).
    state.peakPeriod = std::sqrt(2.0 * kPi * wavelength / kGravity);
    state.frequencyCount = 1;
    state.directionCount = 1;
    state.spreadingExponent = 1e6;  // long-crested
    state.meanDirection = 0.0;      // travelling along +x
    state.seed = 20240607;
    return state;
}

// Time at which the surface at the origin is highest / lowest, found by scan.
void findCrestAndTrough(const WaveField& field, double period, double& crest,
                        double& trough) {
    double best = -1e30, worst = 1e30;
    for (int i = 0; i < 720; ++i) {
        const double t = period * i / 720.0;
        const double e = field.elevation(0, 0, t);
        if (e > best) { best = e; crest = t; }
        if (e < worst) { worst = e; trough = t; }
    }
}

void settle(Ship& ship, const Sea& sea, int steps) {
    for (int i = 0; i < steps; ++i) ship.step(0.02, sea);
}

}  // namespace

// A wave field of zero height must reproduce still water exactly. This is the
// load-bearing check: it drives the wavy code path -- world-space hull,
// integrateBelowSurface, per-quadrature height evaluation -- against an answer
// the flat path already computes a completely different way.
void testZeroAmplitudeWaveFieldMatchesFlatWater() {
    const SeaState calm = oneWave(0.0, 200.0);
    const WaveField field(calm);

    Ship flat = makeBarge();
    flat.initialise(0.0);
    Ship wavy = makeBarge();
    wavy.initialise(0.0);

    Sea sea;
    sea.waves = &field;
    sea.level = 0.0;

    settle(flat, 0.0, 1500);
    settle(wavy, sea, 1500);

    expectNear("a zero-amplitude wave field floats the ship at the still-water draft",
               wavy.diagnostics(sea).draftMidship, flat.diagnostics(0.0).draftMidship, 1e-4);
    expectNear("and gives the same displaced volume",
               wavy.diagnostics(sea).buoyantVolume, flat.diagnostics(0.0).buoyantVolume,
               1e-3 * flat.diagnostics(0.0).buoyantVolume);
    expectTrue("and leaves the ship upright",
               std::abs(wavy.diagnostics(sea).heelDeg) < 0.05);
}

// Held at a fixed position, a hull under a crest must displace more water than
// the same hull under a trough. If it does not, the height field is being
// ignored somewhere between the Sea and the integral.
void testCrestDisplacesMoreThanTrough() {
    const SeaState state = oneWave(4.0, 300.0);
    const WaveField field(state);
    double crest = 0, trough = 0;
    findCrestAndTrough(field, state.peakPeriod, crest, trough);

    Ship ship = makeBarge();
    ship.initialise(0.0);

    Sea atCrest;
    atCrest.waves = &field;
    atCrest.time = crest;
    Sea atTrough = atCrest;
    atTrough.time = trough;

    const double crestVolume = ship.diagnostics(atCrest).buoyantVolume;
    const double troughVolume = ship.diagnostics(atTrough).buoyantVolume;
    const double still = ship.diagnostics(0.0).buoyantVolume;

    expectTrue("a crest displaces more than still water", crestVolume > still);
    expectTrue("a trough displaces less than still water", troughVolume < still);

    // Closed form: over a 300 m wave the 60 m barge is nearly under a uniform
    // elevation, so the extra displacement is close to waterplane area times the
    // local surface rise. Within 20%, which is the wave's curvature over the hull.
    const double area = 60.0 * 16.0;
    const double rise = field.elevation(0, 0, crest);
    expectNear("the extra displacement is the waterplane area times the surface rise",
               crestVolume - still, area * rise, 0.20 * area * rise);
}

// A wave far longer than the ship is locally a slowly tilting flat surface, so
// the ship contours it: at equilibrium it sits a full wave amplitude higher on a
// crest than in still water. A wave far shorter than the ship averages out along
// the hull and barely moves it. Getting these two limits right is most of what
// makes a seakeeping model believable.
void testLongWavesLiftTheShipAndShortWavesDoNot() {
    Ship still = makeBarge();
    still.initialise(0.0);
    settle(still, 0.0, 1500);
    const double stillZ = still.state.position.z;

    // Long wave: 600 m against a 60 m hull.
    {
        const SeaState state = oneWave(3.0, 600.0);
        const WaveField field(state);
        double crest = 0, trough = 0;
        findCrestAndTrough(field, state.peakPeriod, crest, trough);
        const double amplitude = field.elevation(0, 0, crest);

        Ship ship = makeBarge();
        ship.initialise(0.0);
        Sea sea;
        sea.waves = &field;
        sea.time = crest;  // frozen at the crest, so equilibrium is well defined
        settle(ship, sea, 3000);

        const double lift = ship.state.position.z - stillZ;
        expectTrue("the amplitude is worth measuring against", amplitude > 0.5);
        expectNear("a wave ten times the ship length lifts it by the full amplitude",
                   lift, amplitude, 0.25 * amplitude);
    }

    // Short wave: 12 m against the same 60 m hull, five wavelengths along it.
    {
        const SeaState state = oneWave(3.0, 12.0);
        const WaveField field(state);
        double crest = 0, trough = 0;
        findCrestAndTrough(field, state.peakPeriod, crest, trough);
        const double amplitude = field.elevation(0, 0, crest);

        Ship ship = makeBarge();
        ship.initialise(0.0);
        Sea sea;
        sea.waves = &field;
        sea.time = crest;
        settle(ship, sea, 3000);

        const double lift = std::abs(ship.state.position.z - stillZ);
        expectTrue("a wave a fifth of the ship length barely lifts it",
                   lift < 0.30 * amplitude);
    }
}

// The hull mesh must resolve the wave, and the failure mode when it does not is
// the dangerous kind: not noise, but a systematic phase-dependent gain or loss of
// buoyancy. A single panel spanning several wavelengths samples the surface at
// three points and reports whatever those three points happen to say.
void testHullMustResolveTheWavelength() {
    const double amplitude = 1.06, wavelength = 12.0;
    const double k = 2.0 * kPi / wavelength;
    // Phase offset chosen so the error does not cancel by symmetry -- with a
    // symmetric phase a coarse mesh looks deceptively fine.
    auto wave = [&](double x, double) {
        return amplitude * std::cos(k * (x + 0.37 * wavelength));
    };

    // Five whole wavelengths along a 60 m hull: the exact extra displacement is
    // zero, because the crests and troughs cancel over the waterplane.
    double coarse = 0, fine = 0;
    for (int divisions : {1, 2, 8, 64}) {
        Ship barge = makeBarge(divisions);
        TriMesh hull = barge.hull;
        for (Vec3& v : hull.verts) v.z -= 4.0;  // float it at 4 m draft
        const double flat = integrateBelowSurface(hull, [](double, double) { return 0.0; }).volume;
        const double wavy = integrateBelowSurface(hull, wave).volume;
        if (divisions == 1) coarse = wavy - flat;
        if (divisions == 64) fine = wavy - flat;
    }

    expectTrue("a single-panel hull invents a large spurious displacement",
               std::abs(coarse) > 100.0);
    expectNear("a resolved hull correctly cancels a whole number of wavelengths", fine, 0.0,
               1e-6);
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

// --- Where a free surface stops being linear, and what GM is sampled at -------
//
// A metacentric height is the slope of the righting arm *at the origin*, so a
// finite difference only delivers one while the arm is linear across the angle it
// is taken at. A shallow layer in a wide compartment stops being linear almost
// immediately: past `atan(2h/b)` the water has pulled off the high side, the
// surface is narrower than the compartment, and the free-surface moment collapses.
//
// The six tests below establish, in order: that the collapse follows a closed
// form; that the angle it starts at scales as depth over half-breadth; that
// `Diagnostics::gmTransverse` is sampled inside the linear region and is the
// actual slope there; that the ferry's 50 t case reads -3.77 m and not the
// +0.59 m it used to; that a linear region too narrow to resolve is reported as
// one; and that a ship with no free surface is untouched, bit for bit.

// A box barge with a shallow layer in one tank. Everything about it is
// arithmetic: the free-surface moment is `rho mu b^3 l / 12`, the pocketing angle
// is `atan(2h/b)`, and both are independent of the mesh integration that computes
// the answer.
//
// **Permeability is deliberately not 1.** It enters the two closed forms
// differently -- the water occupies the geometric region `V/mu`, so it sets the
// depth and hence the pocketing angle, while the moment is the *water's* mass
// times that region's centroid shift and so carries `mu` linearly. A model that
// dropped it from either would be caught by only one of the assertions.
constexpr double kFsL = 60.0, kFsB = 14.0, kFsD = 9.0;
constexpr double kFsTankL = 48.0, kFsTankB = 11.0;
constexpr double kFsMu = 0.85;
constexpr double kFsTankI = kFsTankB * kFsTankB * kFsTankB * kFsTankL / 12.0;

Ship freeSurfaceBarge() {
    Ship s;
    s.hull = makeBox({-kFsL / 2, -kFsB / 2, 0}, {kFsL / 2, kFsB / 2, kFsD});
    s.deckEdgeZ = kFsD;
    Compartment c;
    c.name = "tank";
    c.permeability = kFsMu;
    c.ventedToAtmosphere = true;
    c.mesh = makeBox({-kFsTankL / 2, -kFsTankB / 2, 3.0}, {kFsTankL / 2, kFsTankB / 2, 6.0});
    s.compartments.push_back(c);
    s.lightshipMass = kFsL * kFsB * 3.5 * kRhoSeawater;
    s.lightshipCog = {0, 0, 5.0};
    s.gyradii = {5, 16, 16};
    s.initialise(0.0);
    return s;
}

// The same barge twice: once with the water loose, once with the identical mass
// bolted down at the identical centre of gravity. Their GZ curves differ by the
// free surface and by nothing else, which is what makes the difference assertable
// against a closed form.
struct FreeSurfacePair {
    Ship liquid, solid;
    double layerDepth = 0;    // m, over the *geometric* region
    double pocketingRad = 0;  // atan(2h/b)
    double correction = 0;    // m, rho mu I / displacement
};

FreeSurfacePair freeSurfacePair(double tonnes) {
    FreeSurfacePair p{freeSurfaceBarge(), freeSurfaceBarge(), 0, 0, 0};
    p.liquid.compartments[0].waterVolume = tonnes * 1000.0 / kRhoSeawater;
    p.liquid.step(1e-9, 0.0);

    const double m = tonnes * 1000.0, lm = p.solid.lightshipMass;
    p.solid.lightshipCog =
        (p.solid.lightshipCog * lm + p.liquid.compartments[0].waterCentroid * m) / (lm + m);
    p.solid.lightshipMass = lm + m;
    p.solid.step(1e-9, 0.0);

    p.layerDepth = p.liquid.compartments[0].waterVolume / kFsMu / (kFsTankL * kFsTankB);
    p.pocketingRad = std::atan(2.0 * p.layerDepth / kFsTankB);
    p.correction =
        kRhoSeawater * kFsMu * kFsTankI / p.liquid.diagnostics(0.0).displacementMass;
    return p;
}

double slopeAtHeel(const Ship& s, const Sea& sea, double eps) {
    return (s.rightingArmAtHeel(eps, sea) - s.rightingArmAtHeel(-eps, sea)) / (2 * eps);
}

// The free-surface loss at `heel`, as a fraction of the loss the linear theory
// predicts. Both ships are sampled at the *same* angle, so the hull's own O(eps^2)
// curvature cancels out of the difference and what is left is the free surface.
double freeSurfaceFraction(const FreeSurfacePair& p, double heel) {
    return (slopeAtHeel(p.solid, 0.0, heel) - slopeAtHeel(p.liquid, 0.0, heel)) / p.correction;
}

// The wedge closed form. Once `tan(heel) > 2h/b` the water is a triangular prism
// of the same area, whose base is `b/sqrt(tau)` wide; its centroid sits
// `b/2 - b/(3 sqrt(tau))` off the upright one, and dividing that by the linear
// theory's `b^2 tan(heel) / 12h` gives this. It is 1 at tau = 1 -- the two regimes
// meet without a step -- and falls away as tau^-1 thereafter.
double pocketedMomentFraction(double tau) {
    return tau <= 1.0 ? 1.0 : 3.0 / tau - 2.0 / (tau * std::sqrt(tau));
}

void testShallowLayerPocketingFollowsTheWedgeClosedForm() {
    const FreeSurfacePair p = freeSurfacePair(30.0);

    // Vacuity: the layer has to be shallow enough that it pockets well inside the
    // angles anyone samples GM at, or every assertion below is about the linear
    // regime and proves nothing.
    expectTrue("the layer is a few centimetres deep",
               p.layerDepth > 0.06 && p.layerDepth < 0.07);
    expectTrue("so it spans the tank only out to 0.012 rad, well inside 0.03",
               p.pocketingRad > 0.011 && p.pocketingRad < 0.0125);

    // Inside the linear region the loss is the textbook one to five decimals. The
    // deviation grows as eps^2 -- 5.0e-7 at a quarter of the pocketing angle,
    // 3.8e-5 at nine tenths -- so both tolerances are what was measured.
    expectTrue("well inside the pocketing angle the loss is rho mu I / displacement",
               std::abs(freeSurfaceFraction(p, 0.25 * p.pocketingRad) - 1.0) < 1e-5);
    expectTrue("and still at nine tenths of it",
               std::abs(freeSurfaceFraction(p, 0.9 * p.pocketingRad) - 1.0) < 1e-4);

    // **The permeability is in the closed form and it is not optional.** Dropping
    // it under-states the loss by 15% here, which is four orders of magnitude
    // outside the tolerance above and would still look entirely plausible.
    expectTrue("with the compartment's permeability genuinely in the answer",
               std::abs(freeSurfaceFraction(p, 0.25 * p.pocketingRad) * kFsMu - 1.0) > 0.14);

    // Past it, against the wedge. This is the part no free-surface *correction*
    // knows about, and it is the whole reason a ship with no GM lolls rather than
    // capsizing.
    double worst = 0, previous = 1.0;
    for (double tau : {1.1, 1.5, 2.0, 3.0, 5.0, 10.0}) {
        // tau is defined on the tangent, so invert it there rather than scaling
        // the angle -- at ten times the pocketing angle the two differ by 8%.
        const double heel = std::atan(tau * 2.0 * p.layerDepth / kFsTankB);
        const double measured = freeSurfaceFraction(p, heel);
        worst = std::max(worst, std::abs(measured / pocketedMomentFraction(tau) - 1.0));
        expectTrue("and the moment falls away monotonically", measured < previous);
        previous = measured;
    }
    expectTrue("the pocketed free-surface moment is 3/tau - 2/tau^1.5", worst < 1e-3);

    // Vacuity again: a closed form that stayed near 1 would make the loop above
    // pass on a model with no pocketing in it at all.
    expectTrue("and that closed form has to actually collapse",
               pocketedMomentFraction(10.0) < 0.25);
}

// The pocketing angle is geometric -- the surface reaches the high side when the
// layer's depth equals half the breadth times the tangent -- so it must scale as
// depth over half-breadth and nothing else.
void testTheLinearRegionScalesAsDepthOverHalfBreadth() {
    // Where the closed form itself crosses 99%, by bisection on tau. Asserting the
    // measured departure against *this* rather than against a remembered constant
    // is what makes the claim "the mesh follows the wedge" instead of "the mesh
    // gives the number it gave last time".
    double lo = 1.0, hi = 4.0;
    for (int k = 0; k < 60; ++k) {
        const double mid = 0.5 * (lo + hi);
        if (pocketedMomentFraction(mid) > 0.99) lo = mid; else hi = mid;
    }
    const double tauAt99 = 0.5 * (lo + hi);
    expectNear("the wedge form is 99% of linear at tau = 1.1291", tauAt99, 1.12912, 1e-4);

    double firstRatio = 0, spread = 0, firstPocket = 0;
    for (double tonnes : {7.5, 15.0, 30.0, 60.0}) {
        const FreeSurfacePair p = freeSurfacePair(tonnes);
        expectNear("the pocketing angle is the depth over the half breadth",
                   p.pocketingRad, std::atan(p.layerDepth / (kFsTankB / 2)), 1e-12);

        // The measured departure: the angle at which the loss has fallen 1% below
        // the linear theory, bisected geometrically.
        double a = 1e-6, b = 0.5;
        for (int k = 0; k < 45; ++k) {
            const double mid = std::sqrt(a * b);
            if (freeSurfaceFraction(p, mid) > 0.99) a = mid; else b = mid;
        }
        const double departure = std::sqrt(a * b);
        expectNear("and the measured 1% departure sits at the closed form's own root",
                   departure / p.pocketingRad, tauAt99, 3e-3);

        if (firstRatio == 0) { firstRatio = departure / p.layerDepth; firstPocket = p.pocketingRad; }
        else expectNear("so the departure angle is proportional to the layer depth",
                        departure / p.layerDepth, firstRatio, 0.02 * firstRatio);
        spread = p.pocketingRad / firstPocket;
    }
    // Vacuity: the sweep has to have actually swept. Eight-fold in depth, and the
    // proportionality above is only worth asserting across that.
    expectTrue("over an eight-fold range of layer depth", spread > 7.9 && spread < 8.1);
}

// What `Diagnostics` publishes has to be the slope of the righting arm near zero,
// and it has to be sampled somewhere that slope means something. Both halves are
// checked against a least-squares fit to the arm itself, which shares no code with
// the central difference in `Ship::diagnostics()`.
void testGmIsSampledInsideTheLinearRegionAndIsTheSlopeThere() {
    // dGZ/dheel by least squares through ten points either side of upright,
    // deliberately over half the angle the ship reports so that the fit is inside
    // whatever window it chose rather than straddling its edge.
    auto fittedSlope = [](const Ship& s, const Sea& sea, double window) {
        double sxx = 0, sxy = 0;
        for (int k = -5; k <= 5; ++k) {
            if (k == 0) continue;
            const double heel = k * window / 5.0;
            sxx += heel * heel;
            sxy += heel * s.rightingArmAtHeel(heel, sea);
        }
        return sxy / sxx;
    };

    const FreeSurfacePair p = freeSurfacePair(30.0);
    const Diagnostics barge = p.liquid.diagnostics(0.0);
    expectTrue("the barge samples GM inside its own pocketing angle",
               barge.gmSampledAtRad > 0 && barge.gmSampledAtRad < p.pocketingRad);
    expectTrue("and says so", barge.gmSlopeConverged);
    // Agreement is asserted in metres, because that is the scale the refinement
    // works to: it stops when halving the angle moves GM by less than 1e-6 m.
    // Measured 3.5e-7 m on the barge and 4.7e-7 m on the ferry.
    expectTrue("and the GM it publishes is the fitted slope of its own righting arm",
               std::abs(barge.gmTransverse - fittedSlope(p.liquid, 0.0, barge.gmSampledAtRad / 2))
                   < 2e-6);

    Sea sea{0.0};
    Ship ferry = game::buildFerry();
    ferry.initialise(sea);
    const int deck = ferry.findCompartment("vehicle_deck");
    ferry.compartments[static_cast<std::size_t>(deck)].waterVolume = 50.0e3 / kRhoSeawater;
    ferry.step(1e-9, sea);
    const Diagnostics d = ferry.diagnostics(sea);
    expectTrue("the ferry too", std::abs(d.gmTransverse -
                                         fittedSlope(ferry, sea, d.gmSampledAtRad / 2)) < 2e-6);

    // The other half of the claim, and the reason any of this was worth doing: the
    // slope at 0.03 rad is a different number, and it is a number of the wrong
    // sign. A ship that went back to sampling there would fail here.
    const double atThirtyMilliradians = slopeAtHeel(ferry, sea, 0.03);
    expectTrue("while the slope at 0.03 rad is positive and four metres away",
               atThirtyMilliradians > 0.4 && atThirtyMilliradians - d.gmTransverse > 4.0);
}

// The 50 t case, which is the one that started this: a 2.9 cm layer on an
// undivided 100 x 19 m vehicle deck. The published GM read **+0.59 m where the
// initial GM is -3.77 m** -- positive, and positive is the unsafe direction on the
// one number every stability judgement in this project keys off.
void testTheFerrysShallowLayerGmIsTheNegativeOne() {
    Sea sea{0.0};
    Ship s = game::buildFerry();
    s.initialise(sea);
    const int deck = s.findCompartment("vehicle_deck");
    s.compartments[static_cast<std::size_t>(deck)].waterVolume = 50.0e3 / kRhoSeawater;
    s.step(1e-9, sea);

    const Diagnostics d = s.diagnostics(sea);
    expectNear("with 50 t on the vehicle deck the ferry's GM is -3.77 m", d.gmTransverse,
               -3.7696, 0.01);
    expectTrue("and the sampling converged", d.gmSlopeConverged);

    // The mechanism, so this is a statement about pocketing rather than about two
    // numbers. The deck's plan area comes off its own mesh: a 5 cm slab on the
    // floor, one-sided upward so it is not measured through a slab half outside
    // the space.
    const Compartment& c = s.compartments[static_cast<std::size_t>(deck)];
    const double planArea =
        integrate(clipToBox(c.mesh, c.bboxLo, {c.bboxHi.x, c.bboxHi.y, c.bboxLo.z + 0.05}))
            .volume / 0.05;
    const double depth = c.waterVolume / c.permeability / planArea;
    const double pocketing = std::atan(2.0 * depth / (planArea / (c.bboxHi.x - c.bboxLo.x)));
    expectTrue("because the layer is under three centimetres deep",
               depth > 0.028 && depth < 0.030);
    expectTrue("and spans the deck only out to a fifth of a degree",
               pocketing > 0.0030 && pocketing < 0.0032);
    expectTrue("so GM is sampled inside that and not at 0.03 rad",
               d.gmSampledAtRad < pocketing && d.gmSampledAtRad > 0.5 * pocketing);

    // **This is the assertion that fails if the sampling angle is put back.** The
    // old fixed 0.03 rad is ten times the angle this layer spans the deck to and
    // delivers a quarter of the free-surface effect; it reported +0.59 m.
    expectTrue("where the answer would be positive and wrong",
               slopeAtHeel(s, sea, 0.03) > 0.5 && slopeAtHeel(s, sea, 0.03) < 0.7);

    // Vacuity: the two GMs have to be measuring the same ship, and the ship has to
    // have been stable before the water arrived.
    Ship dry = game::buildFerry();
    dry.initialise(sea);
    const Diagnostics intact = dry.diagnostics(sea);
    expectNear("intact she has 2.00 m of GM", intact.gmTransverse, 2.0, 0.01);
    expectTrue("sampled at the unrefined angle, there being no free surface to pocket",
               intact.gmSampledAtRad == 0.03);
}

// The free-surface moment `rho mu I` does not depend on how *much* water there is,
// only on the shape of its surface -- while the angle that surface survives to
// does. So the initial GM is **discontinuous** in the amount of loose water, and
// the ferry's vehicle deck straddles that discontinuity inside a factor of two.
//
// Half a litre spread over 1868 m2 is a third of a micron deep and pockets at
// 3e-8 rad. There is no angle anything can measure at which that water has a free
// surface, and the ship correctly reports her own +2.00 m over a window six
// decades wider than the layer. One litre is the same statement one step further
// in: no window is wide enough to establish anything, and `gmSlopeConverged` says
// so rather than publishing a slope with a domain of two nanoradians.
//
// **The number behind an unset flag is deliberately not asserted here.** It is
// whatever the refinement was holding when it ran out of room, and asserting it
// would be asserting that an unreliable number stays unreliable in the same way.
void testAnUnresolvableLinearRegionIsReportedAsOne() {
    Sea sea{0.0};
    const int deck = game::buildFerry().findCompartment("vehicle_deck");

    auto ferryHolding = [&](double cubicMetres) {
        Ship s = game::buildFerry();
        s.initialise(sea);
        s.compartments[static_cast<std::size_t>(deck)].waterVolume = cubicMetres;
        s.step(1e-9, sea);
        return s.diagnostics(sea);
    };

    const Diagnostics litre = ferryHolding(0.001);
    expectTrue("one litre on the vehicle deck has no resolvable linear region",
               !litre.gmSlopeConverged);
    expectTrue("so the refinement bottomed out at 1.8 nanoradians",
               litre.gmSampledAtRad < 2e-9);

    // Half of it, and the ship is back to a perfectly ordinary answer: the layer
    // is too thin to have a free surface at any measurable angle, so what she
    // reports is her own intact GM over a window wider than she will ever roll in.
    const Diagnostics half = ferryHolding(0.0005);
    expectTrue("half a litre converges", half.gmSlopeConverged);
    expectNear("to her own intact GM", half.gmTransverse, 2.0, 0.01);
    expectTrue("over a window six decades wider than that layer pockets at",
               half.gmSampledAtRad > 1e-3);

    // And a tonne, which is still only a 0.6 mm layer, resolves the other way --
    // so the flag distinguishes three regimes rather than firing on any free
    // surface at all.
    const Diagnostics tonne = ferryHolding(1.0);
    expectTrue("a tonne of it converges too", tonne.gmSlopeConverged);
    expectNear("to a GM of -3.78 m", tonne.gmTransverse, -3.7819, 0.01);
    expectTrue("over a window between the two", tonne.gmSampledAtRad > 5e-5 &&
                                                    tonne.gmSampledAtRad < 7e-5);

    // Vacuity: the discontinuity has to actually be a discontinuity. Half a litre
    // and a tonne differ by a factor of 2000 in water and by 5.8 m of GM, in
    // opposite signs, and that is the whole reason a bare GM is not enough.
    expectTrue("and the two answers straddle zero by nearly six metres",
               half.gmTransverse - tonne.gmTransverse > 5.7);
}

// --- What a verdict may be drawn from -----------------------------------------

// The ferry's vehicle deck holding `cubicMetres` of loose water, with the
// geometry of the layer it makes alongside: everything the wedge closed form
// needs, measured off the compartment's own mesh rather than assumed. The deck's
// plan area comes from a 5 cm slab on its floor, one-sided upward so it is not
// measured through a slab half outside the space.
struct DeckLayer {
    Ship   ship;
    double depth = 0;         // m, over the geometric region V/mu
    double breadth = 0;       // m, plan area / deck length
    double pocketingRad = 0;  // atan(2h/b): where the surface leaves the high side
};

DeckLayer ferryDeckLayer(double cubicMetres) {
    DeckLayer layer{game::buildFerry(), 0, 0, 0};
    layer.ship.initialise(0.0);
    const auto deck = static_cast<std::size_t>(layer.ship.findCompartment("vehicle_deck"));
    layer.ship.compartments[deck].waterVolume = cubicMetres;
    layer.ship.step(1e-9, 0.0);

    const Compartment& c = layer.ship.compartments[deck];
    const double planArea =
        integrate(clipToBox(c.mesh, c.bboxLo, {c.bboxHi.x, c.bboxHi.y, c.bboxLo.z + 0.05}))
            .volume / 0.05;
    layer.depth = c.waterVolume / c.permeability / planArea;
    layer.breadth = planArea / (c.bboxHi.x - c.bboxLo.x);
    layer.pocketingRad = std::atan(2.0 * layer.depth / layer.breadth);
    return layer;
}

// The refinement's own ladder of sampling angles, walked from the righting arm
// instead of read off the flag: 0.03 rad, halved, until two successive slopes
// agree to `max(1e-6, 1e-4 |GM|)`. Nothing here touches `Diagnostics`, so an
// assertion that the ladder cannot agree is a statement about the *arm* rather
// than about the code that reports on it.
struct Ladder {
    int    agreedAtRung = -1;    // -1 when none of the twenty-four pairs agree
    double agreedAt = 0;         // rad, the coarser angle of the agreeing pair
    double tightestMargin = 0;   // min over the rungs walked of |change| / tolerance
    std::vector<double> eps;     // the finer angle of each pair, in order
    std::vector<double> change;  // |slope(eps) - slope(2 eps)| at each
};

Ladder walkTheLadder(const Ship& s, const Sea& sea) {
    Ladder l;
    l.tightestMargin = kInfinity;
    double eps = 0.03;
    double gm = slopeAtHeel(s, sea, eps);
    for (int i = 0; i < 24; ++i) {
        const double halved = 0.5 * eps;
        const double gmHalved = slopeAtHeel(s, sea, halved);
        const double moved = std::abs(gmHalved - gm);
        const double tolerance = std::max(1e-6, 1e-4 * std::abs(gm));
        l.eps.push_back(halved);
        l.change.push_back(moved);
        l.tightestMargin = std::min(l.tightestMargin, moved / tolerance);
        if (moved <= tolerance) {
            l.agreedAtRung = i;
            l.agreedAt = eps;
            return l;
        }
        eps = halved;
        gm = gmHalved;
    }
    return l;
}

// **Why one litre cannot be resolved, as algebra rather than as an observation.**
//
// The flag is a claim that no two sampling angles agree, and "we looked and they
// did not" is the weakest possible support for it -- it would hold just as well on
// a refinement that was merely noisy. The wedge closed form says something much
// stronger. With the layer pocketed the slope carries `C f(tau)` of free-surface
// loss, `tau = tan(eps) b / 2h`, so *halving the angle moves the slope by a known
// amount* -- `C [f(tau/2) - f(tau)]`, which is `3C/tau` in the deep-wedge limit and
// therefore **doubles with every halving**. The iteration is not converging slowly
// there; it is walking away from agreement at a rate the algebra predicts, and it
// can only stop when it reaches the linear region below `atan(2h/b)`. A litre puts
// that region at 6.4e-8 rad, five rungs above the refinement's floor, and the
// movement is still 1.5e-3 m -- four times the tolerance -- when the floor arrives.
//
// `C` is not fitted here: it is the full free-surface loss, and it is taken from
// the *neighbouring case that does resolve*. `rho mu I / displacement` is a
// property of the deck rather than of what is standing on it, so the case that has
// an answer is what predicts the case that has none.
void testTheLitreCaseProvablyCannotConverge() {
    const Sea sea{0.0};
    const DeckLayer litre = ferryDeckLayer(0.001);
    const DeckLayer tonne = ferryDeckLayer(1.0);
    const Diagnostics intact = game::buildFerry().diagnostics(sea);

    // The full loss, from the resolved neighbour: she reads +2.00 m dry and
    // -3.78 m once a tonne has spread out, and the difference is the deck's own
    // free-surface moment over her displacement.
    const Diagnostics resolved = tonne.ship.diagnostics(sea);
    const double fullLoss = intact.gmTransverse - resolved.gmTransverse;
    expectNear("the resolved neighbour puts the deck's free-surface loss at 5.78 m",
               fullLoss, 5.7823, 0.02);

    const Ladder l = walkTheLadder(litre.ship, sea);
    expectEqual("with a litre on the deck no pair of angles on the ladder agrees",
                l.agreedAtRung, -1);
    expectEqual("all twenty-four rungs were walked", static_cast<int>(l.eps.size()), 24);
    // Measured 1.385, at the third rung, which is where the free surface first
    // outruns the hull's own curvature. The margin is asserted because "they did
    // not agree" is worth much less if they nearly did.
    expectTrue("and the closest of them misses the tolerance by 38%",
               l.tightestMargin > 1.3 && l.tightestMargin < 1.45);

    // The closed form, rung by rung, over the fourteen where the free surface is
    // what moves the slope. Above them the movement is the hull's own O(eps^2)
    // curvature -- 1.2e-3 m of GM between 0.03 and 0.0075 rad, which swamps a wedge
    // that has barely started -- and below them the layer is a fifth of a micron
    // deep and the volume solve that re-levels it is not exact at that scale either.
    const auto tauOf = [&](double e) {
        return std::tan(e) * litre.breadth / (2.0 * litre.depth);
    };
    int compared = 0;
    double worst = 0;
    for (std::size_t i = 0; i < l.eps.size(); ++i) {
        const double eps = l.eps[i];
        if (eps > 1e-3 || eps < 1.1e-7) continue;
        const double predicted = fullLoss * (pocketedMomentFraction(tauOf(eps)) -
                                             pocketedMomentFraction(tauOf(2 * eps)));
        worst = std::max(worst, std::abs(l.change[i] / predicted - 1.0));
        expectTrue("every rung moves the slope by more than the tolerance it is judged by",
                   l.change[i] > 1e-4 * 3.8);
        ++compared;
    }
    expectEqual("fourteen rungs are governed by the wedge", compared, 14);
    expectTrue("and on every one the movement is the wedge's own, to 5%", worst < 0.05);

    // The doubling itself, which is why refining cannot help: it is a property of
    // `3/tau` and needs no constant at all. Only in the deep wedge, though --
    // `f = 3/tau - 2/tau^1.5` gives a ratio of `2(1 - 1.72/sqrt tau)/(1 - 1.22/sqrt
    // tau)`, which is 1.91 at tau = 157 and falls under 1 as the wedge reaches the
    // full breadth. The window below runs tau from 157 to 1.5e4, where the closed
    // form's own ratio is 1.91 to 1.99 and the measurement is 1.95 to 1.99.
    double slowest = kInfinity, fastest = 0;
    int ratios = 0;
    for (std::size_t i = 1; i < l.eps.size(); ++i) {
        if (l.eps[i] > 1e-3 || l.eps[i] < 1e-5 || l.eps[i - 1] > 1e-3) continue;
        slowest = std::min(slowest, l.change[i] / l.change[i - 1]);
        fastest = std::max(fastest, l.change[i] / l.change[i - 1]);
        ++ratios;
    }
    expectEqual("over six halvings of the deep wedge", ratios, 6);
    expectTrue("each moves the slope about twice as far as the one before",
               slowest > 1.94 && fastest < 2.0);

    // And where it would have to stop, which the floor does not reach: the ladder's
    // finest angle is 0.03 / 2^24 = 1.79e-9 rad.
    const double floorAngle = 0.03 / std::pow(2.0, 24);
    expectTrue("the linear region begins at 6.4e-8 rad",
               litre.pocketingRad > 6.3e-8 && litre.pocketingRad < 6.5e-8);
    expectTrue("which is five rungs above the floor the refinement bottoms out at",
               litre.pocketingRad > 32.0 * floorAngle && litre.pocketingRad < 64.0 * floorAngle);

    // So the ship says so, and that is the only line here that reads the flag.
    expectTrue("so the ship reports her GM as unresolved",
               !litre.ship.diagnostics(sea).gmSlopeConverged);

    // --- The neighbours, which resolve, and for two different reasons -----------
    //
    // A tonne agrees *inside* the linear region: both angles of the agreeing pair
    // are below its own pocketing angle, where the surface spans the deck and the
    // slope is constant. That is a metacentric height in the textbook sense.
    const Ladder resolvedLadder = walkTheLadder(tonne.ship, sea);
    expectEqual("a tonne agrees at the tenth rung", resolvedLadder.agreedAtRung, 9);
    expectTrue("and it agrees inside its own linear region, which is where GM lives",
               resolvedLadder.agreedAt < tonne.pocketingRad &&
                   resolvedLadder.agreedAt > 0.5 * tonne.pocketingRad);
    expectTrue("the ship agrees that it resolved", tonne.ship.diagnostics(sea).gmSlopeConverged);

    // Half a litre agrees for the opposite reason, and five decades the other side
    // of its own pocketing angle: the wedge it makes at the angle it agreed at
    // carries less loss than the tolerance, so there is nothing there to resolve.
    // The closed form is what says so -- 1.47e-4 m of free surface against a
    // 2.0e-4 m tolerance -- rather than the fact that the loop stopped. The pair
    // actually moves 1.83e-4 m, the balance being the hull's own O(eps^2)
    // curvature, which is on the dry ship too and is not a free surface at all.
    const DeckLayer half = ferryDeckLayer(0.0005);
    const Ladder halfLadder = walkTheLadder(half.ship, sea);
    expectEqual("half a litre agrees at the fourth rung", halfLadder.agreedAtRung, 3);
    expectTrue("five decades above its own pocketing angle",
               halfLadder.agreedAt > 1e5 * half.pocketingRad);
    const auto halfTau = [&](double e) {
        return std::tan(e) * half.breadth / (2.0 * half.depth);
    };
    const double wedgeThere =
        fullLoss * (pocketedMomentFraction(halfTau(0.5 * halfLadder.agreedAt)) -
                    pocketedMomentFraction(halfTau(halfLadder.agreedAt)));
    expectTrue("because the wedge there is worth 1.47e-4 m of GM",
               wedgeThere > 1.4e-4 && wedgeThere < 1.55e-4);
    expectTrue("which is under the tolerance the pair is judged by", wedgeThere < 1e-4 * 2.0);
    expectTrue("so that one resolves", half.ship.diagnostics(sea).gmSlopeConverged);

    // Vacuity: the three cases have to be three cases. Half a litre and a tonne
    // straddle zero by 5.8 m with a factor of 2000 of water between them, and the
    // litre sits between them resolving to nothing.
    expectTrue("and the two that resolve disagree about her by 5.8 m",
               half.ship.diagnostics(sea).gmTransverse -
                       tonne.ship.diagnostics(sea).gmTransverse > 5.7);
}

// **The consumer.** A verdict drawn off `gmTransverse < 0` scores an unresolved GM
// exactly as confidently as a resolved one, which is how `--scenario=full`
// published SURVIVED off a 6 mm puddle. Both tools now go through
// `sim::judgeStability`, and an unresolved GM buys a refusal rather than a sign.
void testAnUnresolvedGmIsRefusedAVerdictRatherThanGivenOne() {
    const Sea sea{0.0};
    const Diagnostics half  = ferryDeckLayer(0.0005).ship.diagnostics(sea);
    const Diagnostics litre = ferryDeckLayer(0.001).ship.diagnostics(sea);
    const Diagnostics tonne = ferryDeckLayer(1.0).ship.diagnostics(sea);

    expectTrue("half a litre supports a positive verdict",
               judgeStability(half) == StabilityJudgement::Positive);
    expectTrue("one litre supports none at all",
               judgeStability(litre) == StabilityJudgement::Unresolved);
    expectTrue("a tonne supports a negative one",
               judgeStability(tonne) == StabilityJudgement::Negative);

    // The strings, exactly: `scripts/verify.sh` compares the compiled ferry's
    // outcome line with `ships/ferry.ship`'s as exact strings, so they are asserted
    // here as exact strings too.
    expectTrue("and the outcome line says so",
               std::string(game::floodingOutcome(half)) == "SURVIVED - positive GM, deck edge dry");
    expectTrue("refuses to say so",
               std::string(game::floodingOutcome(litre)) ==
                   "UNDETERMINED - the righting arm is not linear at any angle it can be "
                   "sampled at, so there is no GM to judge her by");
    expectTrue("and says so the other way",
               std::string(game::floodingOutcome(tonne)) == "LOST - negative GM, loll imminent");
    expectTrue("the one-word form too", std::string(game::stabilityWord(litre)) == "UNDETERMINED");
    expectTrue("SURVIVED", std::string(game::stabilityWord(half)) == "SURVIVED");
    expectTrue("LOST", std::string(game::stabilityWord(tonne)) == "LOST");

    // **The vacuity guard, and the whole argument in one line.** The refusal is not
    // merely different from a sign; it replaces a sign that a factor of two in
    // water flips. The old rule read the sign of a number the ship had already
    // flagged as meaningless, and it would have reached a definite verdict from it
    // -- one of the two, and which one is not the point.
    expectTrue("the two resolvable neighbours reach opposite verdicts",
               std::string(game::stabilityWord(half)) != std::string(game::stabilityWord(tonne)));
    expectTrue("while the sign the old rule keyed off is perfectly definite",
               std::isfinite(litre.gmTransverse) && litre.gmTransverse != 0.0);
    expectTrue("and the new verdict is not that sign under another name",
               std::string(game::stabilityWord(litre)) !=
                   std::string(litre.gmTransverse < 0 ? "LOST" : "SURVIVED"));

    // --- What the refusal does not swallow -------------------------------------
    //
    // `afloat` and the heel are observations. They hold whether or not a
    // metacentric height could be measured, they are asked first, and an
    // unresolved GM must not turn a hull with no reserve buoyancy left into an
    // open question.
    Diagnostics d;
    d.afloat = false;
    d.gmSlopeConverged = false;
    d.gmTransverse = 1.0;
    expectTrue("a foundered ship needs no GM to be lost",
               std::string(game::floodingOutcome(d)) == "FOUNDERED");
    expectTrue("in either vocabulary", std::string(game::stabilityWord(d)) == "LOST");
    d.gmSlopeConverged = true;
    expectTrue("and that holds with a perfectly good positive GM",
               std::string(game::stabilityWord(d)) == "LOST");

    // The boundary is pinned on purpose. A GM of exactly zero is *neutral* -- the
    // last state before she loses her initial stability rather than the first
    // state after -- so it takes the positive branch. Nothing in this suite
    // reaches it by accident, which is exactly why it is asserted here rather
    // than left to whichever way the comparison happened to be written.
    d = Diagnostics{};
    d.gmTransverse = 0.0;
    expectTrue("a GM of exactly zero is neutral, and neutral is not negative",
               judgeStability(d) == StabilityJudgement::Positive);
    d.gmTransverse = -std::numeric_limits<double>::denorm_min();
    expectTrue("while the smallest representable step below it is",
               judgeStability(d) == StabilityJudgement::Negative);

    // The rest of the branch table, which is the published vocabulary and is
    // otherwise reachable only by running a 900 s scenario to its end.
    d = Diagnostics{};
    d.gmTransverse = -0.5;
    d.heelDeg = 21.0;
    expectTrue("negative GM past 20 degrees is a loll that has already happened",
               std::string(game::floodingOutcome(d)) ==
                   "LOST - lolled over with negative GM, flooding continuing");
    d.heelDeg = 19.0;
    expectTrue("and short of it, one that has not",
               std::string(game::floodingOutcome(d)) == "LOST - negative GM, loll imminent");
    d.heelDeg = -21.0;
    expectTrue("to port as readily as to starboard",
               std::string(game::floodingOutcome(d)) ==
                   "LOST - lolled over with negative GM, flooding continuing");
    d.gmTransverse = 0.5;
    d.heelDeg = 0.0;
    d.freeboardMin = -0.1;
    expectTrue("a positive GM with the deck edge under is survival with no margin",
               std::string(game::floodingOutcome(d)) ==
                   "SURVIVED but the deck edge is under; no margin left");
    d.freeboardMin = 0.1;
    expectTrue("and with it dry, survival",
               std::string(game::floodingOutcome(d)) == "SURVIVED - positive GM, deck edge dry");
    // The deck edge is geometry rather than stability: it must not become a verdict
    // of its own when there is no GM for it to qualify.
    d.gmSlopeConverged = false;
    d.freeboardMin = -0.1;
    expectTrue("an unresolved GM outranks the deck edge, which is reported beneath it",
               std::string(game::floodingOutcome(d)) ==
                   "UNDETERMINED - the righting arm is not linear at any angle it can be "
                   "sampled at, so there is no GM to judge her by");
}

// **What `--gm-detail` publishes, against the box it can be computed by hand on.**
//
// `sim::largestFreeSurface` reports the geometry the pocketing angle is made of,
// and the ferry's own numbers for it are on the front page. The barge is where it
// can be checked exactly: the tank is a box, so its floor area is `48 x 11 m` to
// the last bit, the layer is `V/mu` over that area, and the angle at which the
// surface leaves the high side is `atan(2h/b)` and nothing else.
//
// Every assertion here is against that arithmetic rather than against a
// remembered number, which is what makes the *ferry's* published angle a
// measurement of the same quantity rather than a coincidence.
void testTheReportedLayerGeometryIsTheBoxsOwnArithmetic() {
    const FreeSurfacePair p = freeSurfacePair(30.0);
    const FreeSurfaceLayer layer = largestFreeSurface(p.liquid);

    expectEqual("the wet tank is the one reported", layer.compartment, 0);
    // The mesh integration has to return the box's plan area, not merely something
    // close to it: this is a clip and an integral over a prism, both exact.
    expectNear("its floor area is the box's own", layer.planArea, kFsTankL * kFsTankB, 1e-9);
    expectNear("its mean breadth is that area over its length", layer.breadth, kFsTankB, 1e-9);
    expectNear("the layer is the water over the geometric region",
               layer.depth, p.layerDepth, 1e-12);
    expectNear("and it pockets at atan(2h/b)", layer.pocketingRad, p.pocketingRad, 1e-12);

    // The fraction is the wedge closed form and is asserted against the test's own
    // copy of it, which shares no code with the engine's.
    int compared = 0;
    for (double tau : {0.5, 1.0, 1.5, 3.0, 10.0, 44.0}) {
        const double heel = std::atan(tau * 2.0 * layer.depth / layer.breadth);
        expectNear("the reported moment fraction is 3/tau - 2/tau^1.5",
                   layer.momentFractionAt(heel), pocketedMomentFraction(tau), 1e-12);
        ++compared;
    }
    expectEqual("over six values of tau either side of the pocketing angle", compared, 6);
    expectTrue("which is a collapse and not a constant: 2.8% left at a hundred times over",
               layer.momentFractionAt(std::atan(100.0 * 2.0 * layer.depth / layer.breadth)) < 0.03);

    // **And the fraction means what it claims to mean**, which no amount of
    // self-consistent algebra would establish: at 0.03 rad the barge's *measured*
    // free-surface loss, taken as the difference between the loose-water ship and
    // the same mass bolted down, has to be that fraction of the linear theory's
    // `rho mu I / displacement`. This is the assertion that makes the front page's
    // "a fixed sample sees 6% of the free surface" a statement about the ship.
    // Measured 0.6888 against the closed form's 0.6886, on a 6.5 cm layer that
    // pockets at 0.0118 rad -- so a fixed 0.03 rad sample on *this* ship sees 69%
    // of its free surface, where the ferry's 6 mm layer leaves 6%.
    expectNear("and the loss measured at 0.03 rad is that fraction of the linear theory",
               freeSurfaceFraction(p, 0.03), layer.momentFractionAt(0.03), 3e-4);
    expectTrue("with the fraction a real departure from linear theory, or this proves nothing",
               layer.momentFractionAt(0.03) > 0.65 && layer.momentFractionAt(0.03) < 0.72);

    // A dry ship has no free surface to report, on the same test `diagnostics()`
    // uses to decide whether to refine the sampling angle at all.
    const FreeSurfaceLayer none = largestFreeSurface(p.solid);
    expectEqual("a ship with the same mass bolted down reports no free surface",
                none.compartment, -1);
    expectTrue("and no geometry with it", none.planArea == 0 && none.pocketingRad == 0);

    // The selection is by floor area, and the ferry is where that matters: a wing
    // tank and the vehicle deck are both wet, and it is the deck -- eighteen times
    // the area -- that sets the angle GM can be measured at.
    Ship ferry = game::buildFerry();
    ferry.initialise(0.0);
    const auto deck = static_cast<std::size_t>(ferry.findCompartment("vehicle_deck"));
    const auto wing = static_cast<std::size_t>(ferry.findCompartment("wing_tank_fwd_p"));
    ferry.compartments[wing].waterVolume = 50.0;
    ferry.compartments[deck].waterVolume = 1.0;
    ferry.step(1e-9, 0.0);
    const FreeSurfaceLayer chosen = largestFreeSurface(ferry);
    expectEqual("the deck is chosen over a wing tank holding fifty times the water",
                chosen.compartment, static_cast<int>(deck));

    // **And the breadth is the deck's *mean*, not the widest part of it.** The
    // bounding box is 20.00 m across -- the deck at its widest, amidships -- while
    // the surface that actually pockets averages 18.68 m over a hull with a fine
    // bow and a fuller stern. The pocketing angle is inversely proportional to
    // this, so the two differ by 7% in the one number the front page publishes
    // about where a metacentric height stops meaning anything, and a breadth taken
    // off the bounding box would be plausible, convenient and wrong.
    const Compartment& deckC = ferry.compartments[deck];
    expectNear("the deck's mean breadth is its area over its length", chosen.breadth,
               chosen.planArea / (deckC.bboxHi.x - deckC.bboxLo.x), 1e-12);
    expectTrue("which is 18.68 m", chosen.breadth > 18.6 && chosen.breadth < 18.8);
    expectTrue("and is not the 20.00 m the bounding box would have given",
               deckC.bboxHi.y - deckC.bboxLo.y - chosen.breadth > 1.3);
    // 1868 m2 of deck against 193 m2 of tank -- and the tank holds fifty times the
    // water. Area is what decides which surface can be resolved, not volume.
    expectTrue("because it is nearly ten times the surface, holding a fiftieth of the water",
               chosen.planArea > 9.0 * largestFreeSurface([&] {
                   Ship only = game::buildFerry();
                   only.initialise(0.0);
                   only.compartments[wing].waterVolume = 50.0;
                   only.step(1e-9, 0.0);
                   return only;
               }()).planArea);
}

// The halving count is the sampling angle said the other way round, and the two
// must not be able to drift apart: `eps = 0.03 * 2^-n` exactly, because halving a
// double is exact in binary.
void testTheHalvingCountAndTheAngleAreTheSameStatement() {
    const Sea sea{0.0};
    int distinct = 0;
    double previous = -1;
    for (double cubicMetres : {0.0, 0.0005, 0.001, 1.0, 50.0e3 / kRhoSeawater}) {
        const Diagnostics d = cubicMetres == 0.0
                                  ? game::buildFerry().diagnostics(sea)
                                  : ferryDeckLayer(cubicMetres).ship.diagnostics(sea);
        expectTrue("the sampling angle is 0.03 halved as many times as it says",
                   d.gmSampledAtRad == 0.03 / std::pow(2.0, d.gmHalvings));
        expectTrue("and the count never exceeds the bound the search is given",
                   d.gmHalvings >= 0 && d.gmHalvings <= 24);
        if (d.gmHalvings != previous) ++distinct;
        previous = d.gmHalvings;
    }
    // Vacuity: an implementation that always reported zero halvings would satisfy
    // the identity above on every ship, so the sweep has to have moved it.
    expectTrue("over a sweep that actually moves the count", distinct >= 4);

    // The two ends of it, named. A dry ship is not refined at all; the litre runs
    // the search to its floor and reports that it did.
    expectEqual("a dry ship is sampled at 0.03 rad with no halvings at all",
                game::buildFerry().diagnostics(sea).gmHalvings, 0);
    const Diagnostics litre = ferryDeckLayer(0.001).ship.diagnostics(sea);
    expectEqual("a litre on the deck runs the search to the bound", litre.gmHalvings, 24);
    expectTrue("which is the case that does not converge", !litre.gmSlopeConverged);
}

// **The exact control.** Refining the sampling angle is a repair to a case that
// only arises with a free surface, so a ship without one has to come out not
// merely close but bit-identical -- which is what lets every published figure
// taken on an intact ship stand un-re-derived. `scripts/check-figures.sh` is where
// the count of those lives; putting it here too would be a second copy to rot.
//
// Asserted as an identity against the expression `diagnostics()` used before this
// existed, rather than against a remembered constant, so it holds under any
// compiler and says what it means.
void testAnIntactShipsGmIsTheUnrefinedDifferenceExactly() {
    Sea sea{0.0};
    Ship ferry = game::buildFerry();
    ferry.initialise(sea);

    Ship barge = freeSurfaceBarge();  // same hull, tank left dry

    int checked = 0;
    for (Ship* s : {&ferry, &barge}) {
        for (double heelDeg : {0.0, 8.0}) {
            s->state.orientation = Quat::fromAxisAngle(Vec3{1, 0, 0}, heelDeg * kDegToRad);
            const Diagnostics d = s->diagnostics(sea);
            const double unrefined =
                (s->rightingArmAtHeel(0.03, sea) - s->rightingArmAtHeel(-0.03, sea)) / (2 * 0.03);
            expectTrue("with no free surface GM is the 0.03 rad difference, bit for bit",
                       d.gmTransverse == unrefined);
            expectTrue("taken at 0.03 rad exactly", d.gmSampledAtRad == 0.03);
            ++checked;
        }
    }
    expectEqual("over four dry configurations", checked, 4);

    // Vacuity: put water in the same barge and the identity must break, or the
    // assertions above would pass on a build that never refines anything.
    Ship wet = freeSurfaceBarge();
    wet.compartments[0].waterVolume = 30.0e3 / kRhoSeawater;
    wet.step(1e-9, 0.0);
    const Diagnostics d = wet.diagnostics(0.0);
    expectTrue("while with a shallow layer in it, it does not",
               d.gmTransverse != slopeAtHeel(wet, 0.0, 0.03) && d.gmSampledAtRad < 0.03);
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

// --- Gas temperature --------------------------------------------------------

// A barge with one sealed compartment and one hole low in its side, so the
// rising water compresses the air it traps. The shape the trapped-air tests all
// need.
Ship makeTrappedAirBarge(double thermalTime) {
    Ship s;
    s.hull = makeBox({-25, -6, 0}, {25, 6, 8});
    s.deckEdgeZ = 8.0;

    Compartment sealed;
    sealed.name = "sealed_void";
    sealed.mesh = makeBox({-10, -4, 0.0}, {10, 4, 3.0});
    sealed.permeability = 1.0;
    sealed.gasThermalTime = thermalTime;
    s.compartments = {sealed};

    Opening hole;
    hole.name = "hole";
    hole.a = kSea;
    hole.b = 0;
    hole.pos = {0, -4, 0.2};
    hole.area = 0.5;
    s.openings = {hole};

    s.lightshipMass = 50.0 * 12.0 * 3.0 * kRhoSeawater;
    s.lightshipCog = {0, 0, 4.0};
    s.gyradii = {4.0, 14.0, 14.0};
    s.initialise(0.0);
    return s;
}

// The closed form. Air trapped in a compartment and compressed by the water
// rising into it, with nowhere to send the heat, follows the isentrope:
//
//     T / T0 = (p / p0) ^ ((gamma - 1) / gamma)
//
// This is the sharp test that the energy is being *carried* rather than reset
// every tick, because an implementation that recomputed the temperature from
// anything else would land on the isotherm instead, and the two differ by 30% by
// the time a compartment is half full.
//
// The same barge with the default thermal time is the control, and it must land
// on Boyle instead -- which is not a second-best answer but the model the
// flooding scenarios are validated against, so both are asserted here.
void testAdiabaticCompressionFollowsTheIsentrope() {
    Ship s = makeTrappedAirBarge(kInfinity);
    const Compartment& c = s.compartments[0];

    const double p0 = c.airPressure;
    const double t0 = c.gasTemperature;
    const double v0 = c.airVolume();
    for (int i = 0; i < 120000; ++i) s.step(0.005, 0.0);

    // Non-vacuous: the compartment has to have actually been compressed, or the
    // isentrope and the isotherm agree trivially at a ratio of one. It settles at
    // p/p0 = 1.264 on a void 3 m deep with the hole 0.2 m off its floor, which is
    // all the head there is to compress with.
    expectTrue("the trapped air was really compressed", c.airPressure > 1.2 * p0);
    expectTrue("and the gas space really shrank", c.airVolume() < 0.9 * v0);
    expectTrue("adiabatic compression heated the trapped air", c.gasTemperature > t0 + 15.0);

    const double ratio = c.airPressure / p0;
    const double expected = t0 * std::pow(ratio, (kGammaAir - 1.0) / kGammaAir);
    // 1e-12 and not something comfortable: this is an exact relation between two
    // stored numbers, not a converged one, and it comes back at 1e-14. A
    // tolerance loose enough to be safe would also be loose enough to pass on a
    // model that had lost the property and was merely in the right region.
    expectNear("trapped air follows T/T0 = (p/p0)^((gamma-1)/gamma)", c.gasTemperature,
               expected, 1e-12 * expected);

    // And the mass never moved, so pV^gamma is the same statement seen from the
    // other side. Asserting it as well is what distinguishes carrying the energy
    // from having fitted the temperature to the pressure after the fact -- the
    // first form would still hold if the temperature were being *derived* from
    // the pressure each tick, and this one would not.
    const double pvg0 = p0 * std::pow(v0, kGammaAir);
    expectNear("and equivalently conserves p V^gamma",
               c.airPressure * std::pow(c.airVolume(), kGammaAir), pvg0, 1e-12 * pvg0);

    // The control: the default gas is isothermal, so the *same* barge must sit on
    // Boyle's law and nowhere near the isentrope. The two answers are far apart,
    // which is the point -- this is the difference the state variable makes.
    Ship cold = makeTrappedAirBarge(0.0);
    for (int i = 0; i < 120000; ++i) cold.step(0.005, 0.0);
    const Compartment& cc = cold.compartments[0];
    expectNear("the default gas stays exactly at ambient", cc.gasTemperature, kTAmbient, 0.0);
    expectNear("and compresses isothermally", cc.airPressure * cc.airVolume(),
               kPatm * cc.grossVolume, 0.02 * kPatm * cc.grossVolume);
    // The two models are being asked the same question and must give visibly
    // different answers, or neither assertion above proves anything. They do, and
    // by more than a temperature: air that is allowed to heat as it is compressed
    // pushes back harder, so the adiabatic void arrests the flooding at 15.4% full
    // where the isothermal one goes to 20.2%. A **24% difference in how much water
    // gets in** is the measurement behind keeping the isothermal model the default
    // -- it is far too large to have absorbed into the published scenarios
    // silently, and far too large for the choice to be a matter of taste.
    expectTrue("the isentrope and the isotherm are far apart here",
               c.gasTemperature > 1.05 * cc.gasTemperature);
    expectTrue("and adiabatic air keeps materially more water out",
               c.fillFraction() < 0.85 * cc.fillFraction());
}

// Two compartments at different temperatures, joined by a high opening and a low
// one, must exchange gas: hot out the top, cold in the bottom, with no net mass
// going anywhere once it settles. This is the whole mechanism of smoke spread,
// and it is what a network at a single temperature structurally cannot express.
//
// **The direction is the assertion, not the magnitude.** A transport term with
// the wrong sign, or a buoyancy head taken about the wrong datum, produces
// exchange of entirely plausible size flowing the wrong way through both holes.
void testBuoyancyDrivesGasThroughAVerticalOpeningPair() {
    auto makePair = [](double tHot) {
        Ship s;
        s.hull = makeBox({-25, -6, 0}, {25, 6, 12});
        s.deckEdgeZ = 12.0;
        auto room = [&](const char* name, double x0, double x1) {
            Compartment c;
            c.name = name;
            c.mesh = makeBox({x0, -4, 2.0}, {x1, 4, 10.0});
            c.permeability = 1.0;
            c.gasThermalTime = kInfinity;  // adiabatic: nothing to cool the smoke
            return c;
        };
        s.compartments = {room("hot", -20, 0), room("cold", 0, 20)};
        s.compartments[0].gasTemperature = tHot;

        // A doorway is not one orifice. `Opening` has a position and an area and
        // no height at all, so a single one of them samples the pressure
        // difference at exactly one elevation and moves gas one way through it --
        // at its own neutral plane it transports nothing whatever. The two-way
        // formulation the docs specify is expressed here as two openings, one
        // each side of where the neutral plane falls, which is also how a
        // two-layer fire model discretises the same doorway.
        auto door = [&](const char* name, double z) {
            Opening o;
            o.name = name;
            o.a = 0;
            o.b = 1;
            o.pos = {0, 0, z};
            o.area = 0.9;
            o.kind = OpeningKind::Door;
            return o;
        };
        s.openings = {door("door_top", 8.5), door("door_bottom", 3.5)};

        s.lightshipMass = 50.0 * 12.0 * 3.0 * kRhoSeawater;
        s.lightshipCog = {0, 0, 4.0};
        s.gyradii = {4.0, 14.0, 14.0};
        s.initialise(0.0);
        return s;
    };

    Ship s = makePair(900.0);
    // Hold the two temperatures where they were put, the way a sustained fire in
    // one space and open sea air in the other would. Left to mix they equalise
    // and the exchange correctly stops, which would make the steady state this
    // test is about unreachable.
    for (int i = 0; i < 4000; ++i) {
        s.compartments[0].gasTemperature = 900.0;
        s.compartments[1].gasTemperature = kTAmbient;
        s.step(0.002, 0.0);
    }

    const double top = s.openings[0].lastGasMassFlow;     // + is hot -> cold
    const double bottom = s.openings[1].lastGasMassFlow;

    expectTrue("hot gas leaves the burning space through the high opening", top > 0);
    expectTrue("cool air is drawn back in through the low opening", bottom < 0);
    expectTrue("the exchange is not a rounding artefact", top > 1e-3);

    // Zero *net* mass flow: the compartment is neither filling nor emptying, it is
    // circulating. Scaled against the size of the exchange itself, because an
    // absolute tolerance would pass on an exchange of zero.
    const double exchange = 0.5 * (std::abs(top) + std::abs(bottom));
    expectNear("no net mass crosses the bulkhead at equilibrium", (top + bottom) / exchange,
               0.0, 1e-3);

    // Guard against vacuity from the other end: the temperatures really did
    // differ, and with them equal the same geometry must exchange *exactly*
    // nothing. Exactly, not nearly -- the buoyancy head is built to vanish by IEEE
    // arithmetic at ambient rather than merely to be small, which is what keeps a
    // cold ship bit-identical to the model that had no temperature in it.
    expectTrue("the two spaces were really at different temperatures",
               s.compartments[0].gasTemperature > 3.0 * s.compartments[1].gasTemperature);

    Ship flat = makePair(kTAmbient);
    for (int i = 0; i < 4000; ++i) flat.step(0.002, 0.0);
    expectNear("at one temperature the high opening moves exactly nothing",
               flat.openings[0].lastGasMassFlow, 0.0, 0.0);
    expectNear("and so does the low one", flat.openings[1].lastGasMassFlow, 0.0, 0.0);
}

// Energy in the gas must balance what crossed the openings, the way the flooding
// solve already balances mass -- and on the ferry rather than on a toy, because a
// toy closes at roundoff while hiding whatever the real network's stiffness does.
// The ferry's air pipes are 0.02 m2 against compartments of hundreds of cubic
// metres, which is the stiffest corner of this solve.
//
// The account is taken over the *sealed* compartments only. A vented space is
// held at atmospheric by exchanging mass with an atmosphere this ledger cannot
// see, so it belongs outside the boundary along with the sea; every opening with
// one endpoint inside and one outside is then a boundary flux.
void testGasEnergyBalancesWhatCrossedTheOpenings() {
    Ship s = game::buildFerry();
    for (Compartment& c : s.compartments) c.gasThermalTime = kInfinity;
    // The hull breaches are shut, so no water rises against the trapped air and
    // the only energy in the ledger is what came through a hole. That is
    // deliberate and it is not a retreat to a toy: this is the ferry's whole
    // network -- eighteen spaces, thirty-odd openings, and air pipes of 0.02 m2
    // serving compartments of hundreds of cubic metres, which is the stiffest
    // corner of this solve and the one a two-box fixture does not have. The pdV
    // work the flooding case adds is a *closed-system* term with an exact answer
    // of its own, and it is asserted against that answer in the isentrope test
    // above rather than folded into this sum where it could only be estimated.
    for (Opening& o : s.openings)
        if (o.kind == OpeningKind::Breach) o.open = false;
    s.initialise(0.0);

    // Set one machinery space alight, in the sense that matters here: hot gas,
    // which is less dense, higher pressure and free to go looking for a way out.
    const int seat = s.findCompartment("engine_room_s");
    expectTrue("the ferry has the compartment this test heats", seat != kSea);
    s.compartments[static_cast<std::size_t>(seat)].gasTemperature = 1100.0;

    const auto& comps = s.compartments;
    auto inside = [&](int i) { return i != kSea && !comps[static_cast<std::size_t>(i)].ventedToAtmosphere; };
    auto energy = [&] {
        double e = 0;
        for (const Compartment& c : comps)
            if (!c.ventedToAtmosphere) e += c.airMass * c.gasTemperature;
        return e;
    };

    const double e0 = energy();
    double crossed = 0.0;   // enthalpy in, in the same kg*K units as m*T
    const double dt = 0.005;
    double worst = 0.0;
    for (int i = 0; i < 40000; ++i) {
        s.step(dt, 0.0);
        for (const Opening& o : s.openings) {
            if (o.lastGasMassFlow == 0.0) continue;
            const double dm = o.lastGasMassFlow * dt;      // + is a -> b
            const double h = kGammaAir * o.lastGasDonorTemperature * dm;
            if (inside(o.b) && !inside(o.a)) crossed += h;
            if (inside(o.a) && !inside(o.b)) crossed -= h;
        }
        worst = std::max(worst, std::abs(energy() - e0 - crossed));
    }

    // Non-vacuous: real energy has to have crossed the boundary, or a solve that
    // moved nothing at all would balance perfectly.
    expectTrue("gas really did cross the ferry's boundary", std::abs(crossed) > 0.02 * e0);
    expectTrue("and the compartment cooled as it vented",
               comps[static_cast<std::size_t>(seat)].gasTemperature < 1000.0);

    // Tight, because this is an exactly-conservative update rather than a
    // converged one: every opening adds the same gamma*T*dm to one side that it
    // takes from the other, so the only error is floating-point summation over
    // 40 000 steps and 30-odd openings.
    expectNear("gas energy balances the enthalpy that crossed the openings", worst / e0, 0.0,
               1e-9);
}

// A barge shell the gas fixtures below hang compartments inside.
Ship makeGasShell() {
    Ship s;
    s.hull = makeBox({-25, -6, 0}, {25, 6, 12});
    s.deckEdgeZ = 12.0;
    s.lightshipMass = 50.0 * 12.0 * 3.0 * kRhoSeawater;
    s.lightshipCog = {0, 0, 4.0};
    s.gyradii = {4.0, 14.0, 14.0};
    return s;
}

// Heat leaves the gas on the time constant it was given, and the closed form for
// that is an exponential. Asserted against exp() rather than against a shape,
// because the whole reason this is an exact relaxation rather than an explicit
// step is that the fire's time constant is seconds while the flooding solve steps
// at milliseconds.
void testGasRelaxesToTheStructureExponentially() {
    Ship s = makeGasShell();
    Compartment c;
    c.name = "hot";
    c.mesh = makeBox({-10, -4, 2.0}, {10, 4, 8.0});
    c.permeability = 1.0;
    c.gasThermalTime = 5.0;
    c.gasTemperature = 800.0;
    s.compartments = {c};   // sealed, no openings: nothing but the relaxation acts
    s.initialise(0.0);

    const double t0 = s.compartments[0].gasTemperature;
    const double dt = 0.01;
    // **An odd number of steps, deliberately.** Flip the sign of the increment in
    // the relaxation and the map becomes T <- Tamb - (T - Tamb) a, whose *square*
    // is exactly the correct two-step map T <- Tamb + (T - Tamb) a^2. Over an even
    // number of steps that mutation therefore lands on the right answer to the
    // last bit and the test proves nothing. This is the error that cancels when it
    // is asked globally, and an odd count is what refuses to let it.
    const int n = 1501;
    for (int i = 0; i < n; ++i) s.step(dt, 0.0);

    const double want = kTAmbient + (t0 - kTAmbient) * std::exp(-n * dt / 5.0);
    expectTrue("the gas started well away from ambient", t0 > 2.0 * kTAmbient);
    expectTrue("and has not finished relaxing", want > kTAmbient + 20.0);
    // 3e-15 measured over 1 500 steps, so the tolerance is set at what a chain of
    // that many exact exponentials is entitled to and no looser.
    expectNear("gas cools to the structure as exp(-t/tau)", s.compartments[0].gasTemperature,
               want, 1e-12 * want);
}

// Charging a rigid space from the atmosphere *heats* it, even though the
// atmosphere is at ambient and the space started at ambient. That is not a quirk:
// the arriving gas is pushed in by the reservoir behind it, so it brings cp T
// while only cv T is needed to sit there, and the difference is the flow work.
//
// It is also the sharpest available check that gas leaving the *sea* leaves at
// ambient and that the transport is enthalpy: the closed form
// T1 = T0 (m0 + gamma dm) / m1 contains gamma explicitly, and an implementation
// that moved internal energy instead would land on T1 = T0 exactly.
void testChargingFromTheSeaHeatsTheCompartment() {
    Ship s = makeGasShell();
    Compartment c;
    c.name = "void";
    c.mesh = makeBox({-10, -4, 2.0}, {10, 4, 8.0});
    c.permeability = 1.0;
    c.gasThermalTime = kInfinity;
    s.compartments = {c};
    Opening o;
    o.name = "pipe";
    o.a = kSea;
    o.b = 0;
    o.pos = {0, 0, 9.0};
    o.area = 0.05;
    o.kind = OpeningKind::Vent;
    s.openings = {o};
    s.initialise(0.0);

    // Pump half the air out and let the sea push it back.
    s.compartments[0].airMass *= 0.5;
    const double m0 = s.compartments[0].airMass;
    const double t0 = s.compartments[0].gasTemperature;
    for (int i = 0; i < 20000; ++i) s.step(0.005, 0.0);

    const Compartment& g = s.compartments[0];
    expectTrue("air really was drawn in from outside", g.airMass > 1.5 * m0);

    // It equalises with the atmosphere *at the vent*, which is not the same
    // number as `airPressure`: that is quoted at the middle of the gas space, and
    // this vent is four metres above it through a gas at 336 K. The gap is 6.9 Pa
    // and it is the datum convention rather than an error, so it is asserted as
    // the convention rather than tolerated as slop.
    const double head = -(kPatm / (kRAir * g.gasTemperature) - kPatm / (kRAir * kTAmbient)) *
                        kGravity * ((o.pos.z + s.state.position.z) - g.gasCentroidWorldZ);
    expectTrue("the datum correction is real and not noise", std::abs(head) > 1.0);
    expectNear("the gas equalises with the atmosphere at the vent", g.airPressure + head, kPatm,
               0.05);

    const double want = t0 * (m0 + kGammaAir * (g.airMass - m0)) / g.airMass;
    expectNear("charging from the sea gives T1 = T0 (m0 + gamma dm) / m1", g.gasTemperature,
               want, 1e-9 * want);
    // Non-vacuous, and the discriminating half: internal-energy transport would
    // have left this at exactly ambient.
    expectTrue("and that is materially hotter than ambient", g.gasTemperature > t0 + 30.0);
}

// The clamp, under the condition it exists for: a big door onto a space whose gas
// volume is almost nothing, so an explicit Torricelli step wants to move far more
// mass than equilibrium can absorb. Master clamped to Boyle's equalising mass;
// with a temperature the equalising mass depends on the temperature the arriving
// gas produces, and the question this test settles is whether that circularity
// needs an implicit solve. It does not -- pressure is linear in the energy
// delivered -- and the evidence is that the pressures converge instead of ringing.
void testStiffGasTransferStaysStableWithATemperature() {
    auto makeStiff = [](double tt) {
        Ship s = makeGasShell();
        auto room = [&](const char* n, double x0, double x1, double fill) {
            Compartment c;
            c.name = n;
            c.mesh = makeBox({x0, -4, 0.0}, {x1, 4, 8.0});
            c.permeability = 1.0;
            c.gasThermalTime = tt;
            c.waterVolume = fill * (x1 - x0) * 8.0 * 8.0;
            return c;
        };
        s.compartments = {room("nearly_full", -20, 0, 0.97), room("pressurised", 0, 20, 0.0)};
        // **Two** doors into the same sliver of gas, not one. With a single
        // opening the Jacobi accumulators are still zero when that opening's
        // clamp is evaluated, so whether the clamp reads them at all -- or reads
        // them with the wrong sign -- makes no difference to anything. A second
        // opening into the same compartment is what gives them a value to be
        // wrong about.
        auto door = [&](const char* n, double z) {
            Opening o;
            o.name = n;
            o.a = 1;
            o.b = 0;
            o.pos = {0, 0, z};      // in the thin gas layer, not under its water
            o.area = 3.0;
            o.dischargeCoeff = 0.9;
            o.kind = OpeningKind::Door;
            return o;
        };
        s.openings = {door("upper", 7.95), door("lower", 7.80)};
        s.initialise(0.0);
        s.compartments[1].airMass *= 3.0;   // three atmospheres against a sliver of gas
        return s;
    };

    for (double tt : {0.0, kInfinity}) {
        Ship s = makeStiff(tt);
        const char* label = tt == 0.0 ? "isothermal" : "adiabatic";
        double peak = 0;
        bool sane = true;
        for (int i = 0; i < 20000; ++i) {
            s.step(0.005, 0.0);
            for (const Compartment& c : s.compartments) {
                sane = sane && std::isfinite(c.airPressure) && c.airPressure > 0;
                peak = std::max(peak, c.airPressure);
            }
        }
        expectTrue(std::string("stiff ") + label + " gas pressure stays finite and positive",
                   sane);
        const double pa = s.compartments[0].airPressure;
        const double pb = s.compartments[1].airPressure;
        // It equalises, and it never went anywhere near the excursion an unclamped
        // explicit step produces on a gas volume this small.
        expectTrue(std::string("stiff ") + label + " pressures equalise",
                   std::abs(pa - pb) < 1e-3 * pa);
        expectTrue(std::string("stiff ") + label + " never rang past the starting pressure",
                   peak < 1.02 * 3.0 * kPatm);
    }

    // And the two regimes are distinguishable in exactly the way the enthalpy
    // argument says they must be: with somewhere to put the heat, nothing moves
    // off ambient; with nowhere, the space being charged ends up hotter than it
    // started and the one blowing down ends up cooler. Both, in one fixture --
    // a global check of "temperature changed" would be satisfied by either sign.
    Ship cold = makeStiff(0.0);
    Ship hot = makeStiff(kInfinity);
    for (int i = 0; i < 20000; ++i) { cold.step(0.005, 0.0); hot.step(0.005, 0.0); }
    expectNear("with a heat sink the receiver stays exactly at ambient",
               cold.compartments[0].gasTemperature, kTAmbient, 0.0);
    expectNear("and so does the donor", cold.compartments[1].gasTemperature, kTAmbient, 0.0);
    expectTrue("without one the charged space heats", hot.compartments[0].gasTemperature > kTAmbient + 40.0);
    expectTrue("and the space that blew down cools", hot.compartments[1].gasTemperature < kTAmbient - 1.0);
}

// The gas head has to reach the *water* too. A lighter column of gas presses down
// less, so a hot space sits at a lower pressure on its own water surface than a
// cold one at the same well-mixed pressure -- and water therefore runs towards
// the fire. Aimed at the branch of sideStateAt that a dry fixture never enters.
void testHotGasPullsWaterThroughASubmergedOpening() {
    Ship s = makeGasShell();
    auto room = [&](const char* n, double x0, double x1, double t) {
        Compartment c;
        c.name = n;
        c.mesh = makeBox({x0, -4, 0.0}, {x1, 4, 9.0});
        c.permeability = 1.0;
        c.gasThermalTime = kInfinity;
        c.waterVolume = 0.5 * (x1 - x0) * 8.0 * 9.0;
        c.gasTemperature = t;
        return c;
    };
    s.compartments = {room("hot", -20, 0, 1200.0), room("cold", 0, 20, kTAmbient)};
    Opening o;
    o.name = "low";
    o.a = 0;
    o.b = 1;
    o.pos = {0, 0, 0.5};           // well below both internal free surfaces
    o.area = 0.4;
    o.kind = OpeningKind::Door;
    s.openings = {o};
    s.initialise(0.0);

    const double w0 = s.compartments[0].waterVolume;
    for (int i = 0; i < 6000; ++i) {
        s.compartments[0].gasTemperature = 1200.0;
        s.compartments[1].gasTemperature = kTAmbient;
        s.step(0.002, 0.0);
    }

    expectTrue("the submerged opening carried water, not gas", s.openings[0].lastFlowWasWater);
    expectTrue("water runs towards the lighter gas column",
               s.compartments[0].waterVolume > w0 + 1e-3);
    expectNear("a water opening reports exactly no gas mass flow",
               s.openings[0].lastGasMassFlow, 0.0, 0.0);
    expectNear("and no donor temperature", s.openings[0].lastGasDonorTemperature, 0.0, 0.0);
}

// An opening that carries gas and is then drowned must stop reporting gas. The
// two assertions above are satisfied by a field that was simply never written, so
// they say nothing about the per-tick reset; this drives the same opening through
// both phases and asks the question in the order that can fail. Without the
// reset, a downstream energy or species account reads last tick's mass against
// this tick's phase and books the same joules twice.
void testADrownedOpeningStopsReportingGasFlow() {
    Ship s = makeGasShell();
    auto room = [&](const char* n, double x0, double x1, double t) {
        Compartment c;
        c.name = n;
        c.mesh = makeBox({x0, -4, 0.0}, {x1, 4, 9.0});
        c.permeability = 1.0;
        c.gasThermalTime = kInfinity;
        c.gasTemperature = t;
        return c;
    };
    s.compartments = {room("hot", -20, 0, 900.0), room("cold", 0, 20, kTAmbient)};
    Opening o;
    o.name = "door";
    o.a = 0;
    o.b = 1;
    o.pos = {0, 0, 2.0};
    o.area = 0.8;
    o.kind = OpeningKind::Door;
    s.openings = {o};
    s.initialise(0.0);

    // Phase one: dry, so the door moves gas. Taken as the largest rate over the
    // phase rather than the last one, because a single opening equalises the two
    // pressures and then correctly falls quiet -- reading only the final tick
    // would confuse "has finished" with "never started".
    double movedGas = 0;
    for (int i = 0; i < 500; ++i) {
        s.compartments[0].gasTemperature = 900.0;
        s.compartments[1].gasTemperature = kTAmbient;
        s.step(0.002, 0.0);
        movedGas = std::max(movedGas, std::abs(s.openings[0].lastGasMassFlow));
    }
    expectTrue("the door was moving gas before it was drowned", movedGas > 1e-6);

    // Phase two: drown it, and it must go quiet on the gas channel immediately.
    for (Compartment& c : s.compartments) c.waterVolume = 0.7 * c.floodableVolume();
    for (int i = 0; i < 200; ++i) {
        s.compartments[0].gasTemperature = 900.0;
        s.compartments[1].gasTemperature = kTAmbient;
        s.step(0.002, 0.0);
    }
    expectTrue("the drowned door now carries water", s.openings[0].lastFlowWasWater);
    expectNear("and reports exactly no gas mass flow once drowned",
               s.openings[0].lastGasMassFlow, 0.0, 0.0);
    expectNear("nor a stale donor temperature", s.openings[0].lastGasDonorTemperature, 0.0, 0.0);
}

// The datum the gas pressure is quoted at, against geometry derived here rather
// than read back out of the ship.
//
// **This is the assertion the rest of the buoyancy tests structurally cannot
// make.** A constant error in the datum is invisible to any pair of openings
// between the same two compartments: it shifts both openings' heads by the same
// amount, the well-mixed pressure absorbs the shift, and the circulation comes
// out identical. What it does change is every comparison against a space with a
// *different* datum -- the sea, or a compartment on another deck. So the datum
// has to be pinned to an independently computed number, and checking it with the
// ship's own `gasCentroidWorldZ` on both sides of the equation would prove
// nothing at all.
//
// Checked heeled as well as upright, because the gas space is bounded along the
// ship's `up`, and while she is upright the transverse terms of that are
// multiplied by a component of exactly zero -- so an error in either of them is
// perfectly hidden until she lists.
void testTheGasDatumSitsAtTheMiddleOfTheGasSpace() {
    // Both signs of both angles. The gas space is bounded by picking, per axis,
    // whichever end of the box lies further along `up` -- so at any one attitude
    // half of those choices are the branch not taken, and an error in it is
    // invisible. Listing to port and to starboard, by the bow and by the stern,
    // is what makes every branch load-bearing in at least one case.
    const double attitudes[5][2] = {{0.0, 0.0}, {0.28, 0.11}, {-0.28, 0.11},
                                    {0.28, -0.11}, {-0.28, -0.11}};
    for (const auto& att : attitudes) {
        const double heel = att[0], trim = att[1];
        Ship s = makeGasShell();
        Compartment c;
        c.name = "space";
        // Deliberately off-centre in y, so a transverse term that is wrong shows
        // up as soon as there is any heel at all to expose it.
        c.mesh = makeBox({-10, 1.0, 2.0}, {10, 4.5, 8.0});
        c.permeability = 1.0;
        c.gasThermalTime = kInfinity;
        c.gasTemperature = 700.0;
        s.compartments = {c};
        // Heel *and* trim. Heeling alone leaves `up` with an x component of
        // exactly zero, and the longitudinal term of the gas-space bound is then
        // multiplied by zero however wrong it is.
        s.state.orientation = Quat::fromAxisAngle(Vec3{1, 0, 0}, heel) *
                              Quat::fromAxisAngle(Vec3{0, 1, 0}, trim);
        s.initialise(0.0);

        // Independently derived: the gas fills the whole space (it is dry), so its
        // middle is halfway between the lowest and highest corners of the box
        // measured along the ship's own up, plus wherever the origin has floated to.
        const Mat3 R = s.state.orientation.toMat3();
        const Vec3 up = R.transposed() * Vec3{0, 0, 1};
        const Vec3 lo{-10, 1.0, 2.0}, hi{10, 4.5, 8.0};
        double bottom = 1e30, top = -1e30;
        for (int k = 0; k < 8; ++k) {
            const Vec3 corner{(k & 1) ? hi.x : lo.x, (k & 2) ? hi.y : lo.y, (k & 4) ? hi.z : lo.z};
            const double d = dot(up, corner);
            bottom = std::min(bottom, d);
            top = std::max(top, d);
        }
        const double want = 0.5 * (bottom + top) + s.state.position.z;
        expectNear(heel == 0.0 ? "the gas datum is the mid-height of the gas space"
                               : "and still is at every heel and trim, either way",
                   s.compartments[0].gasCentroidWorldZ, want, 1e-12);
    }

    // And the heel really did move it, so the pair above is two measurements
    // rather than one repeated.
    Ship a = makeGasShell(), b = makeGasShell();
    Compartment c;
    c.name = "space";
    c.mesh = makeBox({-10, 1.0, 2.0}, {10, 4.5, 8.0});
    c.permeability = 1.0;
    a.compartments = {c};
    b.compartments = {c};
    b.state.orientation = Quat::fromAxisAngle(Vec3{1, 0, 0}, 0.28) *
                          Quat::fromAxisAngle(Vec3{0, 1, 0}, 0.112);
    a.initialise(0.0);
    b.initialise(0.0);
    expectTrue("heeling actually moves the gas datum",
               std::abs(a.compartments[0].gasCentroidWorldZ -
                        b.compartments[0].gasCentroidWorldZ) > 0.05);
}

// A vented space is held at atmospheric by swapping mass with the outside air,
// and that air arrives at ambient. So a hot vented compartment whose gas space is
// growing -- water draining out of it -- must be *diluted* by what comes in and
// cool towards ambient, without ever overshooting past it. Left out, a vented
// space could be filled from the atmosphere while keeping all its heat, which is
// energy from nowhere.
void testAVentedSpaceIsCooledByTheAirItDrawsIn() {
    Ship s = makeGasShell();
    Compartment c;
    c.name = "hold";
    c.mesh = makeBox({-15, -4, 0.0}, {15, 4, 9.0});
    c.permeability = 1.0;
    c.ventedToAtmosphere = true;
    c.gasThermalTime = kInfinity;   // nothing but the incoming air can cool it
    c.gasTemperature = 900.0;
    c.waterVolume = 0.97 * 30.0 * 8.0 * 9.0;   // nearly full: room for a 33x dilution
    s.compartments = {c};
    Pump p;
    p.name = "bilge";
    p.compartment = 0;
    p.capacity = 6.0;
    p.maxHead = 40.0;
    p.on = true;
    s.pumps = {p};
    s.initialise(0.0);

    const double m0 = s.compartments[0].airMass;
    const double t0 = s.compartments[0].gasTemperature;
    double previous = t0;
    bool monotone = true, aboveAmbient = true;
    for (int i = 0; i < 120000; ++i) {
        s.step(0.005, 0.0);
        const double t = s.compartments[0].gasTemperature;
        monotone = monotone && t <= previous + 1e-9;
        aboveAmbient = aboveAmbient && t >= kTAmbient - 1e-9;
        previous = t;
    }

    const Compartment& g = s.compartments[0];
    expectTrue("the pump really did open the gas space up", g.airMass > 20.0 * m0);
    expectTrue("the vented space was diluted and cooled", g.gasTemperature < t0 - 100.0);
    // Diluted thirty-three fold, so it must have come down to *ambient* and not
    // merely downwards. This is what pins the temperature the incoming air is
    // carrying: mixing towards anything else settles visibly short of here.
    expectTrue("dilution drives it towards ambient, not just downwards",
               g.gasTemperature < kTAmbient + 25.0);
    expectTrue("it cooled monotonically", monotone);
    // The sharp half: mixing with a reservoir cannot take it past the reservoir.
    expectTrue("and never overshot below the air it was mixing with", aboveAmbient);
    expectNear("a vented space stays at atmospheric throughout", g.airPressure, kPatm, 1e-9);
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
    testTheTwoPlaneIntegralsAgree();
    testClipMatchesIntegral();
    testSubdivisionTiles();
    testClippedCompartmentStaysInsideHull();
    testWavySurfaceReducesToThePlaneCase();
    testSinusoidalSurfaceAgainstAnalyticVolume();
    testWavySurfaceIsNotSecretlyFlat();
    testSurfaceIsQueriedOncePerVertex();
    testShortStationsAreRepairedAndReported();
    testZeroAmplitudeWaveFieldMatchesFlatWater();
    testCrestDisplacesMoreThanTrough();
    testLongWavesLiftTheShipAndShortWavesDoNot();
    testHullMustResolveTheWavelength();
    testArchimedes();
    testFreeSurfaceEffect();
    testShallowLayerPocketingFollowsTheWedgeClosedForm();
    testTheLinearRegionScalesAsDepthOverHalfBreadth();
    testGmIsSampledInsideTheLinearRegionAndIsTheSlopeThere();
    testTheFerrysShallowLayerGmIsTheNegativeOne();
    testAnUnresolvableLinearRegionIsReportedAsOne();
    testTheLitreCaseProvablyCannotConverge();
    testAnUnresolvedGmIsRefusedAVerdictRatherThanGivenOne();
    testTheReportedLayerGeometryIsTheBoxsOwnArithmetic();
    testTheHalvingCountAndTheAngleAreTheSameStatement();
    testAnIntactShipsGmIsTheUnrefinedDifferenceExactly();
    testTrappedAirArrestsFlooding();
    testAdiabaticCompressionFollowsTheIsentrope();
    testBuoyancyDrivesGasThroughAVerticalOpeningPair();
    testGasEnergyBalancesWhatCrossedTheOpenings();
    testGasRelaxesToTheStructureExponentially();
    testChargingFromTheSeaHeatsTheCompartment();
    testStiffGasTransferStaysStableWithATemperature();
    testHotGasPullsWaterThroughASubmergedOpening();
    testADrownedOpeningStopsReportingGasFlow();
    testTheGasDatumSitsAtTheMiddleOfTheGasSpace();
    testAVentedSpaceIsCooledByTheAirItDrawsIn();
    testMassConservation();
}
