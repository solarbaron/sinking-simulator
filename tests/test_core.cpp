// SPDX-License-Identifier: MIT
//
// Validation of the geometric and hydrostatic core against closed-form answers.
// Everything downstream -- flooding rates, stability, capsize -- is only as good
// as these integrals, so they get checked against algebra rather than eyeballed.
#include "engine/core/geometry.hpp"
#include "engine/sim/ship.hpp"
#include "engine/sim/waves.hpp"
#include "harness.hpp"

#include <cstdio>
#include <string>

using namespace sim;
using testing::expectEqual;
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
    testTrappedAirArrestsFlooding();
    testMassConservation();
}
