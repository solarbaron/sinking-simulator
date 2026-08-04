// SPDX-License-Identifier: MIT
//
// Validation of hull-to-hull contact.
//
// A collision is an unusually testable thing, because almost everything about it
// has a closed form:
//
//   * the interpenetration of two boxes is a box, and its volume, centroid,
//     principal extents, patch areas, patch centroids and patch normals are all
//     arithmetic. Two boxes at 45 degrees intersect in a regular octagon whose
//     area is 2(sqrt 2 - 1) times the square's -- the same machinery asked a
//     question no axis alignment can answer by accident;
//   * the lens between two spheres has an analytic volume, so a tessellated pair
//     gives a *convergence* test rather than a tolerance;
//   * linear momentum is conserved exactly in any collision, and angular momentum
//     about any origin is conserved exactly when both bodies are loaded at one
//     shared point -- which is the thing that is easy to get wrong and that still
//     looks like a collision;
//   * an elastic collision conserves kinetic energy exactly and an inelastic one
//     loses exactly (1 - e^2)/2 times the reduced mass times the approach speed
//     squared, at both ends of the restitution range rather than in the middle;
//   * a penalty contact between two flat faces *is* a linear oscillator, so its
//     duration is pi sqrt(m/K), its peak force u sqrt(K m) and its maximum
//     penetration u sqrt(m/K), all three exactly.
//
// The other half of the suite is negative. A detector that never fires passes
// every test that only asks it not to fire, so every "these do not touch" check
// here has a companion that moves the same two hulls together and requires that
// they do.
#include "engine/sim/collision.hpp"
#include "engine/sim/hullform.hpp"
#include "harness.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <string>
#include <utility>
#include <vector>

using namespace sim;
using testing::expectEqual;
using testing::expectNear;
using testing::expectTrue;

namespace {

constexpr double kPiLocal = 3.14159265358979323846;

// --- Fixtures ----------------------------------------------------------------

TriMesh rotatedBox(const Vec3& lo, const Vec3& hi, const Vec3& axis, double angle,
                   const Vec3& translation = {}) {
    return transformed(makeBox(lo, hi), Quat::fromAxisAngle(axis, angle).toMat3(), translation);
}

// A UV sphere. Closed and outward-wound, and it *inscribes* the true sphere, so
// every derived volume approaches the analytic one from below -- which makes the
// convergence one-sided and therefore readable.
TriMesh makeSphere(double radius, const Vec3& centre, int latitudes, int longitudes) {
    TriMesh m;
    m.verts.push_back(centre + Vec3{0, 0, radius});
    for (int i = 1; i < latitudes; ++i) {
        const double theta = kPiLocal * i / latitudes;
        for (int j = 0; j < longitudes; ++j) {
            const double phi = 2.0 * kPiLocal * j / longitudes;
            m.verts.push_back(centre + Vec3{radius * std::sin(theta) * std::cos(phi),
                                            radius * std::sin(theta) * std::sin(phi),
                                            radius * std::cos(theta)});
        }
    }
    m.verts.push_back(centre + Vec3{0, 0, -radius});
    const auto ring = [&](int i, int j) {
        return static_cast<std::uint32_t>(1 + (i - 1) * longitudes + (j % longitudes));
    };
    const auto south = static_cast<std::uint32_t>(m.verts.size() - 1);
    for (int j = 0; j < longitudes; ++j) m.tris.push_back({0, ring(1, j), ring(1, j + 1)});
    for (int i = 1; i + 1 < latitudes; ++i)
        for (int j = 0; j < longitudes; ++j) {
            m.tris.push_back({ring(i, j), ring(i + 1, j), ring(i + 1, j + 1)});
            m.tris.push_back({ring(i, j), ring(i + 1, j + 1), ring(i, j + 1)});
        }
    for (int j = 0; j < longitudes; ++j)
        m.tris.push_back({south, ring(latitudes - 1, j + 1), ring(latitudes - 1, j)});
    return m;
}

// An L-shaped prism: the cheapest closed mesh that is genuinely **not convex**,
// and whose notch can be made to swallow the point the tetrahedral decomposition
// is taken about. Cross-section (0,0) (3,0) (3,1) (1,1) (1,3) (0,3), area 5.
TriMesh makeLPrism(double zLo, double zHi) {
    const double xy[6][2] = {{0, 0}, {3, 0}, {3, 1}, {1, 1}, {1, 3}, {0, 3}};
    TriMesh m;
    for (const auto& p : xy) m.verts.push_back({p[0], p[1], zLo});
    for (const auto& p : xy) m.verts.push_back({p[0], p[1], zHi});
    // The reflex corner is at (1,1) and every other vertex is visible from (0,0),
    // so a fan from vertex 0 triangulates the L without leaving the polygon.
    for (std::uint32_t i = 1; i + 1 < 6; ++i) {
        m.tris.push_back({0, static_cast<std::uint32_t>(i + 1), i});             // -z
        m.tris.push_back({6, static_cast<std::uint32_t>(6 + i),
                          static_cast<std::uint32_t>(7 + i)});                    // +z
    }
    for (std::uint32_t i = 0; i < 6; ++i) {
        const std::uint32_t j = (i + 1) % 6;
        m.tris.push_back({i, j, static_cast<std::uint32_t>(j + 6)});
        m.tris.push_back({i, static_cast<std::uint32_t>(j + 6),
                          static_cast<std::uint32_t>(i + 6)});
    }
    return m;
}

HullParticulars coarseParticulars() {
    HullParticulars p = s175Particulars();
    // Coarse on purpose: the contact sweep is quadratic in the triangles that
    // survive the overlap box, and nothing asserted here is a function of
    // tessellation. The hull-form suite is where fineness is checked.
    p.stationCount = 15;
    p.waterlineCount = 9;
    return p;
}

// A ship with nothing attached but its own hydrostatics: no radiation, no
// propulsion, no Ikeda. Contact is the thing under test and every other coupling
// is a confounder.
Ship plainShip(const HullParticulars& p) {
    Ship ship;
    ship.hull = makeHullFromParticulars(p);
    ship.deckEdgeZ = p.depth;
    ship.lightshipMass = p.blockCoefficient * p.lengthPp * p.beam * p.draft * kRhoSeawater;
    ship.lightshipCog = {0.0, 0.0, 0.55 * p.depth};
    ship.gyradii = {0.35 * p.beam, 0.25 * p.lengthPp, 0.25 * p.lengthPp};
    return ship;
}

// --- 1. Moments of a closed mesh ---------------------------------------------

void testSolidMomentsOfABox() {
    const Vec3 lo{-1, 2, -3}, hi{1, 6, 3};   // 2 x 4 x 6
    const TriMesh box = makeBox(lo, hi);
    const Vec3 centre = (lo + hi) * 0.5;

    const SolidMoments m = solidMoments(box);
    expectNear("box volume", m.volume, 48.0, 1e-12);
    expectNear("box centroid x", m.centroid().x, centre.x, 1e-12);
    expectNear("box centroid y", m.centroid().y, centre.y, 1e-12);
    expectNear("box centroid z", m.centroid().z, centre.z, 1e-12);

    // Covariance of a uniform box is diag(d^2 / 12), about the centroid, so it
    // must not depend on where the moments were taken from. Taking them about a
    // point 500 m away is the check that the shift is there at all: an
    // implementation that ignored `reference` passes with a reference of zero.
    const Mat3 c0 = m.covariance();
    const SolidMoments far = solidMoments(box, Vec3{500, -300, 250});
    const Mat3 c1 = far.covariance();
    expectNear("covariance xx", c0(0, 0), 4.0 / 12.0, 1e-12);
    expectNear("covariance yy", c0(1, 1), 16.0 / 12.0, 1e-12);
    expectNear("covariance zz", c0(2, 2), 36.0 / 12.0, 1e-12);
    expectNear("covariance xy is zero", c0(0, 1), 0.0, 1e-12);
    double worst = 0;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) worst = std::max(worst, std::abs(c0(i, j) - c1(i, j)));
    expectNear("covariance is independent of the reference point", worst, 0.0, 1e-8);
    expectNear("volume is independent of the reference point", far.volume, 48.0, 1e-9);
    expectNear("centroid is independent of the reference point",
               length(far.centroid() - centre), 0.0, 1e-10);
}

void testSymmetricEigenSolvesAKnownMatrix() {
    // A diagonal matrix conjugated by a rotation: the eigenvalues are known
    // exactly and the eigenvectors are the rotated axes.
    const Mat3 r = Quat::fromAxisAngle(Vec3{1, 2, 3}, 0.7).toMat3();
    Mat3 d = Mat3::zero();
    d(0, 0) = 5.0;
    d(1, 1) = 2.0;
    d(2, 2) = 11.0;
    const Mat3 a = r * d * r.transposed();

    std::array<double, 3> values{};
    std::array<Vec3, 3> vectors{};
    symmetricEigen(a, values, vectors);
    expectNear("largest eigenvalue", values[0], 11.0, 1e-10);
    expectNear("middle eigenvalue", values[1], 5.0, 1e-10);
    expectNear("smallest eigenvalue", values[2], 2.0, 1e-10);
    for (int i = 0; i < 3; ++i) {
        const Vec3 v = vectors[static_cast<std::size_t>(i)];
        const double lambda = values[static_cast<std::size_t>(i)];
        expectNear("eigenvector " + std::to_string(i) + " residual",
                   length(a * v - v * lambda), 0.0, 1e-9);
        expectNear("eigenvector " + std::to_string(i) + " is unit", length(v), 1.0, 1e-12);
    }
}

void testPrincipalBoxRecoversABox() {
    const TriMesh box = makeBox({-1, 2, -3}, {1, 6, 3});
    const PrincipalBox p = principalBox(solidMoments(box));
    expectNear("principal extent 0", p.extent[0], 6.0, 1e-9);
    expectNear("principal extent 1", p.extent[1], 4.0, 1e-9);
    expectNear("principal extent 2", p.extent[2], 2.0, 1e-9);
    expectNear("principal centre y", p.centre.y, 4.0, 1e-12);
    expectTrue("the thinnest axis is the box's own short side", std::abs(p.axis[2].x) > 0.999);
    expectTrue("the axes are right-handed",
               dot(cross(p.axis[0], p.axis[1]), p.axis[2]) > 0.999);
}

// --- 2. Detection against closed forms ---------------------------------------

void testTwoBoxesOverlapExactly() {
    // A = [0,10]^3, B = [9,20] x [2,8] x [3,9].  Overlap = [9,10] x [2,8] x [3,9].
    const TriMesh a = makeBox({0, 0, 0}, {10, 10, 10});
    const TriMesh b = makeBox({9, 2, 3}, {20, 8, 9});
    const HullContact c = hullContact(a, b);

    expectTrue("boxes that overlap are reported as touching", c.touching);
    expectNear("overlap volume", c.volume, 36.0, 1e-9);
    expectNear("overlap centroid x", c.centroid.x, 9.5, 1e-9);
    expectNear("overlap centroid y", c.centroid.y, 5.0, 1e-9);
    expectNear("overlap centroid z", c.centroid.z, 6.0, 1e-9);

    // The overlap is a 6 x 6 x 1 slab, and the equivalent box inverts that exactly.
    expectNear("overlap extent 0", c.extent.x, 6.0, 1e-7);
    expectNear("overlap extent 1", c.extent.y, 6.0, 1e-7);
    expectNear("penetration depth", c.depth, 1.0, 1e-7);
    expectNear("projected patch area", c.patchArea, 36.0, 1e-6);

    // A's +x face inside B is 6 x 6, and nothing else of A is inside B at all.
    expectNear("patch area on A", c.onA.area, 36.0, 1e-9);
    expectNear("patch centroid on A, x", c.onA.centroid.x, 10.0, 1e-9);
    expectNear("patch centroid on A, y", c.onA.centroid.y, 5.0, 1e-9);
    expectNear("patch normal on A is +x", c.onA.normal.x, 1.0, 1e-9);

    // B's share is its -x face (36 m2) plus four 1 x 6 strips of its side faces
    // that the corner of A has swallowed (24 m2). Their normals cancel in pairs,
    // so the mean normal is still exactly -x -- which is the property the contact
    // normal is built on, and it is not the property a single face would show.
    expectNear("patch area on B", c.onB.area, 60.0, 1e-9);
    expectNear("patch centroid on B, x", c.onB.centroid.x, 9.2, 1e-9);
    expectNear("patch centroid on B, y", c.onB.centroid.y, 5.0, 1e-9);
    expectNear("patch normal on B is -x", c.onB.normal.x, -1.0, 1e-9);
    expectNear("the two patch normals oppose exactly", c.normalAgreement, 1.0, 1e-9);
    expectNear("patch separation", c.patchSeparation, 0.8, 1e-9);

    // The contact normal points out of B and into A: A must be pushed towards -x.
    expectNear("contact normal x", c.normal.x, -1.0, 1e-9);
    expectNear("contact normal y", c.normal.y, 0.0, 1e-9);
    expectNear("contact normal z", c.normal.z, 0.0, 1e-9);
    expectTrue("no complaints about the geometry", c.problems.empty());
}

void testOverlapIsSymmetricInItsArguments() {
    // Two hull chunks at a general angle. Swapping the arguments recomputes the
    // whole thing with the roles of "surface" and "solid" exchanged, so agreement
    // is an independent evaluation rather than a restatement.
    const TriMesh a = rotatedBox({-6, -3, -2}, {6, 3, 2}, {1, 2, 0.5}, 0.4);
    const TriMesh b = rotatedBox({-5, -4, -3}, {5, 4, 3}, {0.3, -1, 2}, -0.9, {7, 2, 1});

    const HullContact ab = hullContact(a, b);
    const HullContact ba = hullContact(b, a);
    expectTrue("the two hulls do overlap", ab.touching && ab.volume > 1.0);
    expectNear("volume is the same either way round", ab.volume - ba.volume, 0.0,
               1e-9 * ab.volume);
    expectNear("centroid is the same either way round", length(ab.centroid - ba.centroid), 0.0,
               1e-8);
    expectNear("depth is the same either way round", ab.depth - ba.depth, 0.0, 1e-8);
    expectNear("the normal reverses when the roles do", length(ab.normal + ba.normal), 0.0, 1e-8);
    expectNear("A's patch is A's patch either way round", ab.onA.area - ba.onB.area, 0.0, 1e-8);
}

void testARotatedBoxGivesTheOctagon() {
    // A square and the same square turned 45 degrees intersect in a regular
    // octagon of area 2(sqrt 2 - 1) s^2. Nothing about it is axis aligned, so a
    // detector that quietly assumed alignment cannot reach this answer.
    const double octagon = 2.0 * (std::sqrt(2.0) - 1.0);

    const TriMesh tall = makeBox({-0.5, -0.5, -1}, {0.5, 0.5, 1});
    const TriMesh turned = rotatedBox({-0.5, -0.5, -1}, {0.5, 0.5, 1}, {0, 0, 1}, kPiLocal / 4);
    const HullContact c = hullContact(tall, turned);
    expectNear("octagonal prism volume", c.volume, 2.0 * octagon, 1e-9);
    expectNear("octagonal prism centroid", length(c.centroid), 0.0, 1e-9);

    // Thin one of them down and the smallest principal extent becomes the slab
    // thickness exactly, which makes the projected patch area the octagon itself.
    const TriMesh slab = makeBox({-0.5, -0.5, -0.05}, {0.5, 0.5, 0.05});
    const HullContact s = hullContact(slab, turned);
    expectNear("slab overlap volume", s.volume, 0.1 * octagon, 1e-10);
    expectNear("slab penetration depth", s.depth, 0.1, 1e-9);
    expectNear("slab projected patch area", s.patchArea, octagon, 1e-8);
    expectTrue("the slab's normal is vertical", std::abs(s.normal.z) > 0.999);
}

void testContainmentIsExact() {
    // A box entirely inside a ship hull. The overlap is the box, its patch on the
    // box is the box's whole surface, and the hull contributes no patch at all --
    // three answers a partial-overlap test would never distinguish.
    const HullParticulars p = coarseParticulars();
    const TriMesh hull = makeHullFromParticulars(p);
    expectTrue("the test hull is a closed manifold", isClosedManifold(hull));

    const Vec3 lo{-2, -1, 4}, hi{2, 1, 6};
    const TriMesh box = makeBox(lo, hi);
    const double boxVolume = 4.0 * 2.0 * 2.0;
    const double boxArea = 2 * (4 * 2 + 4 * 2 + 2 * 2);

    const HullContact c = hullContact(box, hull);
    expectNear("a contained box overlaps by its own volume", c.volume, boxVolume, 1e-9);
    expectNear("its whole surface is the patch", c.onA.area, boxArea, 1e-9);
    expectNear("the hull contributes no patch at all", c.onB.area, 0.0, 1e-9);
    expectNear("the overlap centroid is the box centroid",
               length(c.centroid - (lo + hi) * 0.5), 0.0, 1e-9);
    expectNear("the overlap's longest extent is the box's", c.extent.x, 4.0, 1e-7);
    expectNear("and its depth is the box's shortest side", c.depth, 2.0, 1e-7);

    // The normal is meaningless when one solid swallows the other, and the model
    // says so rather than inventing one: a closed surface's mean normal is zero.
    expectTrue("a swallowed solid reports that it has no usable normal", !c.problems.empty());
}

void testANonConvexSolidNeedsNegativeTetrahedra() {
    // Every fixture above is convex, and in a convex overlap the point the
    // tetrahedral decomposition is taken about always lands *inside* the solid --
    // so every tetrahedron is positively oriented and the signs never do any work.
    // Replacing the sign with a constant +1 passes all of them. This one it cannot
    // pass: the L's bounding box centre sits in the notch, outside the solid, so
    // the decomposition is a difference of overlapping cones and the answer is
    // wrong without the signs.
    const TriMesh ell = makeLPrism(0.0, 2.0);
    expectTrue("the L-prism is a closed manifold", isClosedManifold(ell));
    expectNear("and encloses its 5 m2 section over 2 m", integrate(ell).volume, 10.0, 1e-12);

    // The notch really is empty, so the apex really is outside: a small box placed
    // where the bounding box centre falls overlaps the L not at all.
    const HullContact notch = hullContact(makeBox({1.45, 1.45, 0.9}, {1.55, 1.55, 1.1}), ell);
    expectNear("the bounding box centre lies in the notch, outside the solid", notch.volume, 0.0,
               0.0);

    // A box across the notch meets the L in two rectangles: 2 x 0.5 of the
    // horizontal arm and 0.5 x 1.5 of the vertical one, over 2 m of height.
    const TriMesh box = makeBox({0.5, 0.5, -1.0}, {2.5, 2.5, 3.0});
    const double expected = 2.0 * (2.0 * 0.5 + 0.5 * 1.5);
    const double centre = 2.0 * (1.0 * 1.5 + 0.75 * 0.75) / 3.5;
    const HullContact c = hullContact(box, ell);
    expectNear("a non-convex overlap volume", c.volume, expected, 1e-9);
    expectNear("its centroid, x", c.centroid.x, centre, 1e-9);
    expectNear("its centroid, y", c.centroid.y, centre, 1e-9);
    expectNear("its centroid, z", c.centroid.z, 1.0, 1e-9);
    expectNear("and the same answer with the arguments swapped",
               hullContact(ell, box).volume - c.volume, 0.0, 1e-9);
}

void testHullsTouchingExactlyDoNotOverlap() {
    // The other half of the coincident-face convention. Two boxes sharing the
    // plane x = 10 with *opposing* outward normals are touching from outside and
    // bound no overlap at all, so the shared face must contribute nothing. Give it
    // half weight -- the treatment the aligned case needs -- and a pair of hulls
    // laid alongside each other reports a large phantom volume out of a contact
    // patch of exactly zero thickness.
    const TriMesh a = makeBox({0, 0, 0}, {10, 10, 10});
    const TriMesh flush = makeBox({10, 2, 3}, {20, 8, 9});
    const HullContact c = hullContact(a, flush);
    // Not *identically* zero, and it should not be asserted as such. Two hulls a
    // beam apart give exactly 0 because no triangle of either is inside the other
    // and every polygon clip comes back empty; here the triangles lie exactly on
    // each other, the clips produce real polygons of zero thickness, and their
    // flux cancels down to the last bits rather than to nothing -- measured at
    // 1e-30 m3, which is twenty-four orders below the volume at which this reports
    // contact at all.
    expectTrue("hulls touching exactly enclose no volume worth the name", c.volume < 1e-20);
    expectTrue("and are not reported as in contact", !c.touching);
    // The volume alone does not pin the convention down: give the shared face half
    // weight on both sides and the two fluxes, being equal and opposite, still
    // cancel to zero. The patches do not cancel, and 36 m2 of reported contact
    // patch between two hulls that are merely alongside each other is exactly the
    // load case a structural model must never be handed.
    expectNear("nor is any of it reported as contact patch on A", c.onA.area, 0.0, 0.0);
    expectNear("nor on B", c.onB.area, 0.0, 0.0);

    // The companion, a millimetre in, so this is not a test of a detector that
    // has stopped firing.
    const TriMesh inside = makeBox({9.999, 2, 3}, {20, 8, 9});
    const HullContact d = hullContact(a, inside);
    expectTrue("a millimetre of interference does register", d.touching);
    expectNear("with the volume that millimetre implies", d.volume, 0.001 * 36.0, 1e-9);
    expectTrue("and a millimetre of real patch on each hull",
               d.onA.area > 35.0 && d.onB.area > 35.0);

    // The two boxes above do not actually reach the opposed-normal branch: their
    // shared plane is also where the overlap box's centre lands, so the
    // tetrahedra over that face have their apex in their own base plane and are
    // dropped as degenerate before the convention is ever consulted. A fixture
    // that *does* reach it needs the flush faces somewhere other than the middle
    // of the overlap box, which needs a non-convex solid: here a box laid into the
    // L-prism's notch, its -x face flush against the notch's inner wall at x = 1,
    // with the overlap box centred at (2, 2, 1) two metres away from it.
    const TriMesh ell = makeLPrism(0.0, 2.0);
    const HullContact e = hullContact(ell, makeBox({1.0, 1.5, 0.5}, {4.0, 2.5, 1.5}));
    expectTrue("a box laid flush into the notch encloses nothing", e.volume < 1e-20);
    expectTrue("and is not in contact", !e.touching);
    expectNear("with no patch on the notch wall", e.onA.area, 0.0, 0.0);
    expectNear("and none on the box", e.onB.area, 0.0, 0.0);

    // A box that is flush against the notch wall *and* overlaps the other arm.
    // The refusal of the flush face is then not the whole story -- part of that
    // face is claimed by a tetrahedron it is not coplanar with, which no
    // per-face rule can undo -- so the answer there is a convention and the
    // routine has to say so rather than present it as a measurement.
    const HullContact g = hullContact(ell, makeBox({1.0, 0.2, 0.5}, {4.0, 2.5, 1.5}));
    expectTrue("a real overlap alongside a flush face is still contact", g.touching);
    expectTrue("and the coincident opposed faces are reported", !g.problems.empty());
    // The L's own share of the patch *is* exact here, and it is what pins the
    // refusal down: 0.8 m2 of its outer face at x = 3, 2.0 m2 of the notch floor
    // at y = 1, and nothing at all from the 1.5 m2 of notch wall lying flush
    // against the box. Count that wall and this reads 3.55.
    expectNear("the refused wall contributes no patch on the L", g.onA.area, 2.8, 1e-12);
    expectTrue("a hull pair with no coincident faces reports nothing",
               hullContact(ell, makeBox({1.2, 0.2, 0.5}, {4.0, 2.5, 1.5})).problems.empty());

    // Its companion, a centimetre into the wall: 1 m2 of wall, one square metre
    // by a centimetre of box, and the four slivers of the box's own sides.
    const HullContact f = hullContact(ell, makeBox({0.99, 1.5, 0.5}, {4.0, 2.5, 1.5}));
    expectTrue("a centimetre into the wall does register", f.touching);
    expectNear("with the volume that centimetre implies", f.volume, 0.01, 1e-12);
    expectNear("one square metre of notch wall", f.onA.area, 1.0, 1e-12);
    expectNear("and the box's face plus its four slivers", f.onB.area, 1.04, 1e-12);
}

void testAnInsideOutHullIsRefusedRatherThanMisintegrated() {
    // Every volume integral in this engine returns nonsense on a mesh wound the
    // wrong way, and it returns it quietly -- the hull that shipped with a 40%
    // displacement error was found by a manifold check added while doing something
    // else. The overlap of two solids can never be bigger than either of them, so
    // that is checked rather than assumed, and an inverted hull trips it.
    const TriMesh good = makeBox({0, 0, 0}, {10, 10, 10});
    TriMesh inverted = makeBox({9, 2, 3}, {20, 8, 9});
    for (Tri& t : inverted.tris) std::swap(t.b, t.c);
    expectTrue("the fixture really is inside out", integrate(inverted).volume < 0);

    for (int order = 0; order < 2; ++order) {
        const HullContact c =
            order == 0 ? hullContact(good, inverted) : hullContact(inverted, good);
        expectTrue("an inside-out hull is reported, not silently integrated, order " +
                       std::to_string(order),
                   c.problems.size() >= 2);
        expectTrue("and nothing is claimed to be in contact, order " + std::to_string(order),
                   !c.touching);
    }
    // Guard: the same two boxes wound correctly say nothing at all.
    expectTrue("the correctly wound pair has no complaints",
               hullContact(good, makeBox({9, 2, 3}, {20, 8, 9})).problems.empty());
}

void testTheAnswerDoesNotDependOnTheOverlapBoxRounds() {
    // Tightening the overlap box is a cost optimisation and must not be anything
    // else. Each round moves the point the moments are taken about, so agreement
    // between one round and four is a genuinely separate evaluation of the same
    // integral rather than a restatement of it.
    const HullParticulars p = coarseParticulars();
    const TriMesh hull = makeHullFromParticulars(p);
    const TriMesh other = transformed(hull, Quat::fromAxisAngle({0, 0, 1}, -kPiLocal / 2).toMat3(),
                                      {30.0, 0.5 * p.beam + 0.5 * p.lengthPp - 6.0, 0.0});

    ContactParams one, many;
    one.boxRounds = 1;
    many.boxRounds = 4;
    const HullContact a = hullContact(hull, other, one);
    const HullContact b = hullContact(hull, other, many);
    expectTrue("the fixture is a real overlap", a.touching && a.volume > 1.0);
    expectNear("volume does not depend on how hard the box was tightened",
               (a.volume - b.volume) / a.volume, 0.0, 1e-9);
    expectNear("nor does the centroid", length(a.centroid - b.centroid), 0.0, 1e-6);
    expectNear("nor the patch area", (a.patchArea - b.patchArea) / a.patchArea, 0.0, 1e-6);
}

void testSphereLensConverges() {
    // Two equal spheres of radius R with centres d apart overlap in a lens of
    // volume 2 pi h^2 (3R - h) / 3 with h = R - d/2. A tessellated sphere
    // inscribes the real one, so the error is one-sided and must fall like the
    // square of the facet size.
    const double radius = 5.0, separation = 8.0;
    const double h = radius - 0.5 * separation;
    const double exact = 2.0 * kPiLocal * h * h * (3.0 * radius - h) / 3.0;
    expectTrue("the lens is not a trivially small quantity", exact > 20.0);

    double previous = 0;
    for (int level = 0; level < 3; ++level) {
        const int lat = 8 << level, lon = 16 << level;
        const TriMesh a = makeSphere(radius, {0, 0, 0}, lat, lon);
        const TriMesh b = makeSphere(radius, {separation, 0, 0}, lat, lon);
        const HullContact c = hullContact(a, b);
        const double error = exact - c.volume;
        expectTrue("a tessellated lens never exceeds the analytic one, level " +
                       std::to_string(level),
                   error > -1e-9);
        // Second order, asserted as an order rather than as a tolerance: halving
        // the facet size must quarter the error. Measured 3.73 and 3.92.
        if (level > 0)
            expectNear("the lens error falls fourfold per refinement, level " +
                           std::to_string(level),
                       previous / error, 4.0, 0.6);
        previous = error;
    }
    expectNear("the finest lens matches the closed form to 2%", previous / exact, 0.0, 0.02);
}

// --- 3. Detection must not fire when it should not ----------------------------

void testHullsApartDoNotTouchAndHullsTogetherDo() {
    const HullParticulars p = coarseParticulars();
    const TriMesh hull = makeHullFromParticulars(p);

    // Two hulls abeam, swept past one another over their whole length, at gaps
    // from a full beam of clear water down to a quarter of a metre. Nothing may
    // fire anywhere, and the reported volume must be *exactly* zero rather than
    // small -- with no triangle of either hull inside the other, every polygon
    // clip is empty and there is nothing for a rounding error to come from.
    for (double gap : {p.beam, 5.0, 1.0, 0.25}) {
        int fired = 0;
        double worst = 0;
        for (int i = -20; i <= 20; ++i) {
            const double surge = i * (p.lengthPp / 20.0);
            const TriMesh other =
                transformed(hull, Mat3::identity(), {surge, p.beam + gap, 0.0});
            const HullContact c = hullContact(hull, other);
            worst = std::max(worst, c.volume);
            if (c.touching) ++fired;
        }
        expectEqual("clear water is never contact, gap " + std::to_string(gap), fired, 0);
        expectNear("and the overlap is exactly zero, gap " + std::to_string(gap), worst, 0.0, 0.0);
    }

    // The companion. The same sweep, moved together until the shells interfere,
    // must fire -- otherwise the checks above are a test of a detector that never
    // fires at all.
    int touched = 0;
    double deepest = 0;
    for (int i = -12; i <= 12; ++i) {
        const double surge = i * (p.lengthPp / 40.0);
        const TriMesh other = transformed(hull, Mat3::identity(), {surge, p.beam - 3.0, 0.0});
        const HullContact c = hullContact(hull, other);
        if (c.touching) {
            ++touched;
            deepest = std::max(deepest, c.volume);
        }
    }
    expectTrue("overlapping the same hulls by three metres does fire", touched >= 20);
    expectTrue("and the overlap is a real volume, not a rounding artefact", deepest > 500.0);

    // Crossing overhead: the bounding boxes overlap in plan, so the broad phase
    // cannot reject and the answer has to come from the geometry.
    const TriMesh above = transformed(hull, Quat::fromAxisAngle({0, 0, 1}, kPiLocal / 2).toMat3(),
                                      {0, 0, p.depth + 2.0});
    const HullContact crossing = hullContact(hull, above);
    expectTrue("a hull passing clean overhead is not contact", !crossing.touching);
    expectNear("and its overlap is exactly zero", crossing.volume, 0.0, 0.0);

    // ... while the same crossing hull lowered into the deck does fire.
    const TriMesh into = transformed(hull, Quat::fromAxisAngle({0, 0, 1}, kPiLocal / 2).toMat3(),
                                     {0, 0, p.depth - 3.0});
    expectTrue("and lowering it into the deck does", hullContact(hull, into).touching);
}

// --- 4. The impulse solver: conservation ---------------------------------------

ContactBody heavyBody(double mass, const Vec3& gyradii, const Vec3& cog, const Vec3& velocity,
                      const Vec3& omega) {
    ContactBody body;
    body.mass = mass;
    body.inertia = Mat3::zero();
    body.inertia(0, 0) = mass * gyradii.x * gyradii.x;
    body.inertia(1, 1) = mass * gyradii.y * gyradii.y;
    body.inertia(2, 2) = mass * gyradii.z * gyradii.z;
    body.cog = cog;
    body.velocity = velocity;
    body.angularVelocity = omega;
    return body;
}

struct Pair {
    ContactBody a, b;
    Vec3 point;
    Vec3 normal;
};

// Deliberately unequal, deliberately spinning, and deliberately not aligned with
// anything: a symmetric fixture would let a wrong lever arm cancel.
Pair genericPair() {
    Pair p;
    p.a = heavyBody(2.4e7, {9.0, 44.0, 44.0}, {-40, 90, -2}, {0.3, -6.0, 0.05},
                    {0.01, -0.02, 0.004});
    p.b = heavyBody(1.7e7, {8.0, 38.0, 38.0}, {12, 0, -1}, {2.5, 0.2, -0.03},
                    {-0.005, 0.01, 0.02});
    p.point = {30.0, 12.4, 1.5};
    p.normal = normalize(Vec3{0.12, 1.0, -0.05});
    return p;
}

Vec3 totalMomentum(const Pair& p) { return p.a.momentum() + p.b.momentum(); }
Vec3 totalAngular(const Pair& p, const Vec3& about) {
    return p.a.angularMomentum(about) + p.b.angularMomentum(about);
}
double totalKinetic(const Pair& p) { return p.a.kineticEnergy() + p.b.kineticEnergy(); }
double momentumScale(const Pair& p) {
    return p.a.mass * length(p.a.velocity) + p.b.mass * length(p.b.velocity);
}
double angularScale(const Pair& p, const Vec3& about) {
    return length(p.a.inertia * p.a.angularVelocity) + length(p.b.inertia * p.b.angularVelocity) +
           p.a.mass * length(p.a.cog - about) * length(p.a.velocity) +
           p.b.mass * length(p.b.cog - about) * length(p.b.velocity);
}

void testImpulseConservesMomentumAndAngularMomentum() {
    for (double e : {0.0, 0.35, 1.0}) {
        const Pair before = genericPair();
        Pair p = before;
        const Vec3 systemCog = (p.a.cog * p.a.mass + p.b.cog * p.b.mass) / (p.a.mass + p.b.mass);
        const Vec3 elsewhere{1234.0, -567.0, 89.0};

        const Vec3 p0 = totalMomentum(p);
        const Vec3 l0 = totalAngular(p, systemCog);
        const Vec3 lElse0 = totalAngular(p, elsewhere);
        const double scaleP = momentumScale(p);
        const double scaleL = angularScale(p, systemCog);
        const double scaleElse = angularScale(p, elsewhere);
        const std::string at = ", e = " + std::to_string(e);

        const ImpulseSolution s = normalImpulse(p.a, p.b, p.point, p.normal, e);
        expectTrue("the pair is closing" + at, s.approachSpeed > 1.0);
        expectTrue("the impulse is not trivial" + at, length(s.impulse) > 1e6);
        applyImpulse(p.a, p.b, p.point, s.impulse);

        expectNear("linear momentum is conserved to machine precision" + at,
                   length(totalMomentum(p) - p0) / scaleP, 0.0, 1e-14);
        expectNear("angular momentum about the system centre of mass" + at,
                   length(totalAngular(p, systemCog) - l0) / scaleL, 0.0, 1e-14);
        // And about a point a kilometre away, which is the same statement only if
        // the two loads really did act at one shared point.
        expectNear("angular momentum about an arbitrary origin" + at,
                   length(totalAngular(p, elsewhere) - lElse0) / scaleElse, 0.0, 1e-14);

        // Guard against a vacuous conservation check: both bodies must actually
        // have changed, in translation *and* in rotation. Conservation of nothing
        // is exact.
        expectTrue("body A's velocity changed" + at,
                   length(p.a.velocity - before.a.velocity) > 0.01);
        expectTrue("body B's velocity changed" + at,
                   length(p.b.velocity - before.b.velocity) > 0.01);
        expectTrue("body A's spin changed" + at,
                   length(p.a.angularVelocity - before.a.angularVelocity) > 1e-6);
        expectTrue("body B's spin changed" + at,
                   length(p.b.angularVelocity - before.b.angularVelocity) > 1e-6);
    }
}

void testAnImpulseOnTheWrongLeverBreaksAngularMomentum() {
    // The instrument, demonstrated. The failure that still looks like a collision
    // is a load applied to the two bodies at *different* points -- momentum stays
    // perfect, the ships still separate, and only angular momentum notices. If
    // this negative control ever stops failing, the angular check above has
    // stopped being a check.
    Pair p = genericPair();
    const Vec3 systemCog = (p.a.cog * p.a.mass + p.b.cog * p.b.mass) / (p.a.mass + p.b.mass);
    const Vec3 p0 = totalMomentum(p);
    const Vec3 l0 = totalAngular(p, systemCog);

    const ImpulseSolution s = normalImpulse(p.a, p.b, p.point, p.normal, 0.4);
    const Vec3 wrong = p.point + Vec3{5.0, 0, 0};   // five metres along the hull
    if (p.a.mass > 0) p.a.velocity += s.impulse / p.a.mass;
    if (p.b.mass > 0) p.b.velocity -= s.impulse / p.b.mass;
    p.a.angularVelocity += inverse(p.a.inertia) * cross(p.point - p.a.cog, s.impulse);
    p.b.angularVelocity -= inverse(p.b.inertia) * cross(wrong - p.b.cog, s.impulse);

    expectNear("a mismatched contact point still conserves linear momentum exactly",
               length(totalMomentum(p) - p0) / momentumScale(p), 0.0, 1e-14);
    expectTrue("but angular momentum sees it",
               length(totalAngular(p, systemCog) - l0) / angularScale(p, systemCog) > 1e-4);
}

void testRestitutionSetsTheEnergyLossExactly() {
    // Both ends of the range, not only somewhere convenient in between.
    for (double e : {0.0, 0.25, 0.5, 0.75, 1.0}) {
        Pair p = genericPair();
        const std::string at = " at e = " + std::to_string(e);
        const double k0 = totalKinetic(p);
        const ImpulseSolution s = normalImpulse(p.a, p.b, p.point, p.normal, e);
        applyImpulse(p.a, p.b, p.point, s.impulse);
        const double k1 = totalKinetic(p);

        const double predicted =
            0.5 * (1.0 - e * e) * s.effectiveMass * s.approachSpeed * s.approachSpeed;
        expectNear("energy lost matches (1 - e^2)/2 m u^2" + at, (k0 - k1 - predicted) / k0, 0.0,
                   1e-13);
        expectNear("and the solver's own figure agrees" + at,
                   (s.energyLost - predicted) / std::max(k0, 1.0), 0.0, 1e-15);

        // Separation speed is exactly e times the approach speed, which is the
        // definition the energy figure follows from.
        const Vec3 relative = p.a.velocityAt(p.point) - p.b.velocityAt(p.point);
        expectNear("separation speed is e times approach" + at, dot(relative, p.normal),
                   e * s.approachSpeed, 1e-9);
        if (e == 0.0)
            expectNear("a fully plastic contact leaves no normal closing at all",
                       dot(relative, p.normal), 0.0, 1e-12);
        if (e == 1.0)
            expectNear("a perfectly elastic contact conserves kinetic energy exactly",
                       (k1 - k0) / k0, 0.0, 1e-14);
    }
}

void testEffectiveMassIsTheReducedMass() {
    // Contact on the line joining the two centres of gravity, normal along it:
    // the rotational terms vanish identically and the effective mass must be the
    // textbook reduced mass. A spurious r x n contribution shows up here.
    ContactBody a = heavyBody(3.0e7, {10, 40, 40}, {-50, 0, 0}, {5, 0, 0}, {0, 0, 0});
    ContactBody b = heavyBody(1.0e7, {10, 40, 40}, {50, 0, 0}, {-1, 0, 0}, {0, 0, 0});

    // The normal points out of B into A. A sits at -x, so a normal of -x is the
    // one that points out of B towards A, and the pair is closing at 6 m/s.
    const ImpulseSolution s = normalImpulse(a, b, {0, 0, 0}, {-1, 0, 0}, 0.0);
    expectNear("approach speed", s.approachSpeed, 6.0, 1e-12);
    expectNear("reduced mass", s.effectiveMass, 3.0e7 * 1.0e7 / 4.0e7, 1e-3);

    // Reversing the sign of the normal reverses the sense of "approach", which is
    // what makes the solver refuse to act on a pair that is separating.
    const ImpulseSolution back = normalImpulse(a, b, {0, 0, 0}, {1, 0, 0}, 0.0);
    expectNear("the reversed normal reports separation", back.approachSpeed, -6.0, 1e-12);
    expectNear("and delivers no impulse at all", length(back.impulse), 0.0, 0.0);

    // An immovable body -- zero mass, zero inertia -- is the infinite-mass limit,
    // so the effective mass collapses to the moving body's own.
    ContactBody wall;
    ContactBody hit = heavyBody(1.0e7, {10, 40, 40}, {-50, 0, 0}, {5, 0, 0}, {0, 0, 0});
    const ImpulseSolution w = normalImpulse(hit, wall, {0, 0, 0}, {-1, 0, 0}, 0.0);
    expectNear("against an immovable body the effective mass is the ship's", w.effectiveMass,
               1.0e7, 1e-3);
}

void testHeadOnIsSymmetricAndYawFree() {
    // Two identical ships, mirror images about x = 0, closing along their common
    // axis with the contact exactly on it.
    ContactBody a = heavyBody(2.0e7, {9, 42, 42}, {-100, 0, 0}, {6, 0, 0}, {0, 0, 0});
    ContactBody b = heavyBody(2.0e7, {9, 42, 42}, {100, 0, 0}, {-6, 0, 0}, {0, 0, 0});
    const ImpulseSolution s = normalImpulse(a, b, {0, 0, 0}, {-1, 0, 0}, 0.4);
    applyImpulse(a, b, {0, 0, 0}, s.impulse);

    expectTrue("the head-on collision happened", length(s.impulse) > 1e7);
    expectNear("the two ships end mirror images in surge", a.velocity.x + b.velocity.x, 0.0, 1e-9);
    expectNear("outgoing speed is e times incoming", a.velocity.x, -0.4 * 6.0, 1e-9);
    expectNear("no sway appears from nowhere", a.velocity.y, 0.0, 0.0);
    expectNear("no heave appears from nowhere", a.velocity.z, 0.0, 0.0);
    expectNear("and no yaw at all", a.angularVelocity.z, 0.0, 0.0);
    expectNear("nor on the struck ship", b.angularVelocity.z, 0.0, 0.0);
}

void testGlancingBlowYawsTheRightWay() {
    // A striker hits the port side of a ship, forward of her centre of gravity.
    // The blow must push her bow to starboard: a negative yaw rate about +z, with
    // a magnitude the lever and the impulse fix exactly.
    const double lever = 45.0;
    const Vec3 normal{0, 1, 0};   // out of the struck ship's port side, into the striker

    const auto blow = [&](double alongship) {
        ContactBody striker = heavyBody(2.0e7, {9, 42, 42}, {0, 90, 0}, {0, -6, 0}, {0, 0, 0});
        ContactBody struck = heavyBody(2.0e7, {9, 42, 42}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0});
        const Vec3 point{alongship, 12.7, 0};
        const ImpulseSolution s = normalImpulse(striker, struck, point, normal, 0.2);
        applyImpulse(striker, struck, point, s.impulse);
        return std::make_pair(struck, s);
    };

    const auto forward = blow(lever);
    expectTrue("the glancing blow landed", length(forward.second.impulse) > 1e7);
    expectTrue("the struck ship is pushed to starboard", forward.first.velocity.y < -0.1);
    expectTrue("and yaws bow-to-starboard", forward.first.angularVelocity.z < 0);
    expectNear("with the yaw rate the lever demands", forward.first.angularVelocity.z,
               -lever * length(forward.second.impulse) / forward.first.inertia(2, 2), 1e-9);

    // The mirror image of the same blow, aft of the cog, yaws the other way. If it
    // did not, the sign above would be a property of the solver, not the geometry.
    const auto aft = blow(-lever);
    expectTrue("a blow aft of the cog yaws the other way", aft.first.angularVelocity.z > 0);
    expectNear("by exactly as much", aft.first.angularVelocity.z + forward.first.angularVelocity.z,
               0.0, 1e-12);

    // And a blow on the centre of gravity yaws not at all.
    const auto centred = blow(0.0);
    expectNear("a blow through the cog produces no yaw", centred.first.angularVelocity.z, 0.0,
               1e-15);
    expectTrue("though it still pushes her sideways", centred.first.velocity.y < -0.1);
}

// --- 5. The penalty contact ----------------------------------------------------

// A free rigid body carrying its own shape, integrated semi-implicitly. Small
// enough to be obviously right, which matters when it is the instrument.
struct FreeBody {
    ContactBody body;
    Vec3 position{};
    Quat orientation{};
    TriMesh shape;

    TriMesh worldShape() const { return transformed(shape, orientation.toMat3(), position); }
    void integrate(const Vec3& force, const Vec3& moment, double dt) {
        if (body.mass > 0) body.velocity += force * (dt / body.mass);
        body.angularVelocity += inverse(body.inertia) * moment * dt;
        position += body.velocity * dt;
        orientation = orientation.integrated(body.angularVelocity, dt);
        body.cog = position;
    }
};

void setBlock(FreeBody& f, double mass) {
    f.shape = makeBox({-10, -5, -5}, {10, 5, 5});
    f.body.mass = mass;
    f.body.inertia = Mat3::zero();
    f.body.inertia(0, 0) = mass * 16.0;
    f.body.inertia(1, 1) = mass * 36.0;
    f.body.inertia(2, 2) = mass * 36.0;
}

struct PenaltyRun {
    double duration = 0;
    double peakForce = 0;
    double maxPenetration = 0;
    double restitution = 0;
    double worstMomentumError = 0;
    double worstSpin = 0;
    double workDone = 0;
    double kineticLost = 0;
};

// Two flat-faced blocks, closing head on. The overlap of two faces is area times
// penetration exactly, so the contact is a linear oscillator with a known period,
// a known peak force and a known maximum penetration.
PenaltyRun runFlatCollision(double dt, const ContactMaterial& material, double closingSpeed) {
    const double mass = 2.0e7;
    FreeBody a, b;
    setBlock(a, mass);
    setBlock(b, mass);
    a.position = {-10.02, 0, 0};
    b.position = {10.02, 0, 0};
    a.body.cog = a.position;
    b.body.cog = b.position;
    a.body.velocity = {0.5 * closingSpeed, 0, 0};
    b.body.velocity = {-0.5 * closingSpeed, 0, 0};

    PenaltyRun run;
    const double kineticStart = a.body.kineticEnergy() + b.body.kineticEnergy();

    bool wasTouching = false;
    for (int step = 0; step * dt < 5.0; ++step) {
        const Vec3 momentum0 = a.body.momentum() + b.body.momentum();
        const double scaleP = a.body.mass * length(a.body.velocity) +
                              b.body.mass * length(b.body.velocity) + 1.0;

        const HullContact contact = hullContact(a.worldShape(), b.worldShape());
        ContactLoad load;
        if (contact.touching) {
            load = contactLoad(contact, a.body, b.body, material);
            run.duration += dt;
            run.peakForce = std::max(run.peakForce, load.normalForce);
            run.maxPenetration = std::max(run.maxPenetration, contact.depth);
            const Vec3 relative = a.body.velocityAt(load.point) - b.body.velocityAt(load.point);
            run.workDone += -dot(load.force, relative) * dt;
            wasTouching = true;
        } else if (wasTouching) {
            break;
        }

        a.integrate(load.force, cross(load.point - a.body.cog, load.force), dt);
        b.integrate(-load.force, cross(load.point - b.body.cog, -load.force), dt);

        run.worstMomentumError = std::max(
            run.worstMomentumError,
            length(a.body.momentum() + b.body.momentum() - momentum0) / scaleP);
        run.worstSpin = std::max(run.worstSpin, length(a.body.angularVelocity) +
                                                    length(b.body.angularVelocity));
    }
    run.restitution = (b.body.velocity.x - a.body.velocity.x) / closingSpeed;
    run.kineticLost = kineticStart - (a.body.kineticEnergy() + b.body.kineticEnergy());
    return run;
}

void testPenaltyContactIsTheOscillatorItClaimsToBe() {
    ContactMaterial material;
    material.dissipation = 0.0;
    material.friction = 0.0;
    material.stiffness = 2.0e7;

    const double closing = 6.0;
    const double area = 10.0 * 10.0;                 // the blocks' 10 x 10 face
    const double springRate = material.stiffness * area;
    const double reduced = 2.0e7 * 2.0e7 / 4.0e7;
    const double omega = std::sqrt(springRate / reduced);

    const PenaltyRun run = runFlatCollision(1e-3, material, closing);
    expectTrue("the blocks actually collided", run.peakForce > 1e8);
    expectNear("contact duration is pi sqrt(m/K)", run.duration, kPiLocal / omega, 5e-3);
    expectNear("maximum penetration is u sqrt(m/K)", run.maxPenetration, closing / omega, 5e-3);
    expectNear("peak force is u sqrt(K m)",
               run.peakForce / (closing * std::sqrt(springRate * reduced)), 1.0, 5e-3);
    expectNear("an undamped penalty contact is perfectly elastic", run.restitution, 1.0, 5e-3);
    expectNear("and returns all the energy it took",
               run.kineticLost / (0.5 * reduced * closing * closing), 0.0, 1e-2);
    expectNear("linear momentum is conserved every step", run.worstMomentumError, 0.0, 1e-14);
    // A contact resultant applied anywhere but on the line of centres would spin
    // these blocks. It is on it to within the flux sum's own rounding, which puts
    // the overlap centroid about 1e-17 m off axis and the resulting spin fifteen
    // orders of magnitude below anything that matters.
    expectNear("a centred contact produces no spin", run.worstSpin, 0.0, 1e-14);
    expectTrue("the contact was resolved over many steps, not one", run.duration / 1e-3 > 100);
}

void testPenaltyRestitutionConvergesWithTheTimestep() {
    // The restitution error is not monotone step by step -- where the last step
    // lands relative to separation moves it around -- so the convergence is
    // asserted over a fourfold refinement, where second order has to show through
    // the noise. Measured 1.04e-4 at dt = 4 ms and 6.8e-6 at 1 ms, a factor of 15
    // against the 16 a second-order method owes.
    ContactMaterial material;
    material.dissipation = 0.0;
    material.friction = 0.0;
    const double coarse = std::abs(runFlatCollision(4e-3, material, 6.0).restitution - 1.0);
    const double fine = std::abs(runFlatCollision(1e-3, material, 6.0).restitution - 1.0);
    expectTrue("the coarse step is already elastic to a part in three thousand", coarse < 3e-4);
    expectTrue("the fine step is elastic to a part in fifty thousand", fine < 2e-5);
    expectTrue("and refining the step fourfold cuts the error by much more than fourfold",
               fine * 8.0 < coarse);
}

void testDissipationSetsRestitution() {
    // Hunt-Crossley's small-dissipation law is e = 1 - (2/3) * dissipation * u.
    // It is asymptotic, so the honest assertion is on the *slope* as the
    // dissipation goes to zero rather than on the value at some convenient point:
    // measured 3.92 against the predicted 4.00 at 0.005 s/m, falling away as the
    // higher-order terms take over. Asserting the value at 0.04 s/m instead would
    // be asserting the size of the term the law neglects.
    const double closing = 6.0;
    const double reduced = 2.0e7 * 2.0e7 / 4.0e7;
    const double closingEnergy = 0.5 * reduced * closing * closing;

    double last = 2.0, lastSlope = 1e30;
    for (double dissipation : {0.0, 0.005, 0.01, 0.02, 0.04}) {
        ContactMaterial material;
        material.friction = 0.0;
        material.dissipation = dissipation;
        const PenaltyRun run = runFlatCollision(1e-3, material, closing);
        const std::string at = " at dissipation " + std::to_string(dissipation);
        expectTrue("more dissipation is always less restitution" + at, run.restitution < last);
        if (dissipation > 0) {
            const double slope = (1.0 - run.restitution) / dissipation;
            if (dissipation <= 0.005)
                expectNear("the small-dissipation slope is (2/3) u", slope, 4.0, 0.15);
            expectTrue("and the law softens rather than steepens" + at, slope < lastSlope);
            lastSlope = slope;
        }
        last = run.restitution;
    }

    // The energy bookkeeping. `work` is a rectangle rule on the contact power at
    // the start of each step, so it is first order: the discrepancy against the
    // kinetic energy actually lost must halve when the step does, and does --
    // 6.19%, 3.07%, 1.53% of the closing energy at 4, 2 and 1 ms.
    ContactMaterial defaults;
    defaults.friction = 0.0;
    double previousGap = 0;
    for (int level = 0; level < 3; ++level) {
        const double dt = 4e-3 / (1 << level);
        const PenaltyRun run = runFlatCollision(dt, defaults, closing);
        const double gap = std::abs(run.workDone - run.kineticLost) / closingEnergy;
        if (level > 0)
            expectNear("the work integral's error halves with the step, level " +
                           std::to_string(level),
                       previousGap / gap, 2.0, 0.2);
        previousGap = gap;
    }

    // The default material at the speed a ram happens at. Not the "nearly
    // plastic" the first-order law would suggest, because that law is not valid
    // this far out -- but inelastic enough that the great majority of the closing
    // energy is absorbed, which is what a ship collision does.
    const PenaltyRun run = runFlatCollision(1e-3, defaults, closing);
    expectTrue("the default material is strongly inelastic at 6 m/s",
               run.restitution > 0.15 && run.restitution < 0.27);
    expectTrue("and absorbs the great majority of the closing energy",
               run.kineticLost > 0.93 * closingEnergy);

    // And it is speed dependent in the right direction: a gentle touch is more
    // elastic than a ram, because less of it goes into damage.
    const PenaltyRun gentle = runFlatCollision(1e-3, defaults, 1.0);
    expectTrue("a gentle touch is markedly more elastic than a ram",
               gentle.restitution > run.restitution + 0.3);
}

void testFrictionRemovesEnergyWithoutBreakingConservation() {
    // The same blocks, but sliding across each other while they press together.
    const double mass = 2.0e7;
    ContactMaterial material;
    material.dissipation = 0.0;
    material.friction = 0.4;

    FreeBody a, b;
    setBlock(a, mass);
    setBlock(b, mass);
    // Started *apart*, not already interfering. Beginning inside the spring hands
    // the pair stored energy they did not arrive with, and their kinetic energy
    // then rises across the contact however much friction takes out of it -- which
    // is a fixture that quietly inverts the thing being measured.
    a.position = {-10.05, 0, 0};
    b.position = {10.05, 0, 0};
    a.body.cog = a.position;
    b.body.cog = b.position;
    a.body.velocity = {3.0, 2.0, 0};   // closing and sliding
    b.body.velocity = {-3.0, 0.0, 0};

    // An origin off every axis of the problem, so the orbital part of the angular
    // momentum is genuinely non-zero and the check is not 0 == 0.
    const Vec3 origin{5, -7, 3};
    const double kinetic0 = a.body.kineticEnergy() + b.body.kineticEnergy();
    const double angularStart =
        length(a.body.angularMomentum(origin) + b.body.angularMomentum(origin));
    expectTrue("the fixture carries real angular momentum", angularStart > 1e8);

    double worstMomentum = 0, worstAngular = 0, slipWork = 0;
    const double dt = 1e-3;
    for (int step = 0; step < 3000; ++step) {
        const Vec3 momentum0 = a.body.momentum() + b.body.momentum();
        const Vec3 angular0 = a.body.angularMomentum(origin) + b.body.angularMomentum(origin);

        const HullContact contact = hullContact(a.worldShape(), b.worldShape());
        ContactLoad load;
        if (contact.touching) {
            load = contactLoad(contact, a.body, b.body, material);
            slipWork += length(load.force - contact.normal * load.normalForce) *
                        length(load.slip) * dt;
        }
        a.integrate(load.force, cross(load.point - a.body.cog, load.force), dt);
        b.integrate(-load.force, cross(load.point - b.body.cog, -load.force), dt);

        const Vec3 momentum1 = a.body.momentum() + b.body.momentum();
        const Vec3 angular1 = a.body.angularMomentum(origin) + b.body.angularMomentum(origin);
        const double scaleP = mass * (length(a.body.velocity) + length(b.body.velocity)) + 1.0;
        const double scaleL = length(angular0) + length(angular1) + 1.0;
        worstMomentum = std::max(worstMomentum, length(momentum1 - momentum0) / scaleP);
        worstAngular = std::max(worstAngular, length(angular1 - angular0) / scaleL);
        if (!contact.touching && step * dt > 0.05) break;
    }
    const double kinetic1 = a.body.kineticEnergy() + b.body.kineticEnergy();
    expectTrue("friction did real work", slipWork > 1e6);
    expectTrue("friction removed energy", kinetic1 < kinetic0);
    expectTrue("the sliding was slowed", a.body.velocity.y < 1.99);
    expectNear("friction still conserves linear momentum exactly", worstMomentum, 0.0, 1e-14);
    expectNear("friction still conserves angular momentum exactly", worstAngular, 0.0, 1e-13);
}

void testAContactNeverPulls() {
    // The dissipation term is proportional to the approach rate, so a pair that is
    // separating faster than 1 / dissipation drives the bracket negative. Left
    // unclamped, the contact then *sucks the hulls together* -- an attractive force
    // between two ships that are already coming apart, which looks like nothing at
    // all in a summary and quietly adds energy.
    const TriMesh a = makeBox({0, 0, 0}, {10, 10, 10});
    const TriMesh b = makeBox({9.5, 2, 3}, {20, 8, 9});
    const HullContact c = hullContact(a, b);
    expectTrue("the fixture is overlapping", c.touching);

    ContactMaterial material;
    material.dissipation = 0.8;
    material.friction = 0.0;
    // The contact normal is -x, so A moving further -x is A separating.
    ContactBody moving = heavyBody(1e7, {5, 20, 20}, {5, 5, 5}, {-3.0, 0, 0}, {0, 0, 0});
    ContactBody still = heavyBody(1e7, {5, 20, 20}, {15, 5, 6}, {0, 0, 0}, {0, 0, 0});
    const ContactLoad load = contactLoad(c, moving, still, material);
    expectTrue("the pair really is separating fast enough to invert the bracket",
               load.approachRate < -1.0 / material.dissipation);
    expectNear("a separating contact pushes with exactly nothing", load.normalForce, 0.0, 0.0);
    expectNear("and applies no force at all", length(load.force), 0.0, 0.0);

    // Still compressive while closing, so the clamp has not simply switched the
    // contact off.
    ContactBody closing = heavyBody(1e7, {5, 20, 20}, {5, 5, 5}, {3.0, 0, 0}, {0, 0, 0});
    expectTrue("but a closing contact still pushes",
               contactLoad(c, closing, still, material).normalForce > 1e6);
}

void testTheLoadActsAtTheCentreOfPressure() {
    // Contact pressure goes as the local penetration, so the resultant of that
    // pressure field acts at the *depth-weighted* centroid of the patch -- which is
    // the overlap solid's own centroid, not the patch's area centroid. On this
    // fixture the two are 0.3 m apart, which on a ship-length lever is a moment
    // error of the same order as the one the contact is trying to produce.
    const TriMesh a = makeBox({0, 0, 0}, {10, 10, 10});
    const TriMesh b = makeBox({9, 2, 3}, {20, 8, 9});
    const HullContact c = hullContact(a, b);
    const ContactBody moving = heavyBody(1e7, {5, 20, 20}, {5, 5, 5}, {2.0, 0, 0}, {0, 0, 0});
    const ContactBody still = heavyBody(1e7, {5, 20, 20}, {15, 5, 6}, {0, 0, 0}, {0, 0, 0});
    const ContactLoad load = contactLoad(c, moving, still, ContactMaterial{});

    expectNear("the load acts at the overlap's own centroid",
               length(load.point - c.centroid), 0.0, 0.0);
    // The guard: the three candidate points must actually be distinguishable here,
    // or this asserts nothing.
    expectTrue("and that is not the same point as either patch centroid",
               length(c.centroid - c.onA.centroid) > 0.4 &&
                   length(c.centroid - c.onB.centroid) > 0.25);
}

// --- 6. Two ships ---------------------------------------------------------------

void testContactBodyIsReferredToTheCentreOfGravity() {
    // Two conversions live in contactBodyOf() and both are invisible on a ship that
    // happens to be upright and pointing along +x -- which is what every scenario
    // in this file would otherwise have.
    HullParticulars p = coarseParticulars();
    Ship ship = plainShip(p);
    ship.gyradii = {9.0, 40.0, 44.0};   // three distinct radii, so no axis is degenerate
    ship.initialise(Sea(0.0));
    ship.state.orientation = Quat::fromAxisAngle(Vec3{0.3, -0.8, 0.5}, 0.9);
    ship.state.velocity = {4.0, -1.0, 0.3};
    ship.state.angularVelocity = {0.05, -0.03, 0.08};

    const Mat3 r = ship.state.orientation.toMat3();
    expectTrue("the ship is at a genuinely non-trivial attitude", std::abs(r(0, 0)) < 0.9);

    const ContactBody body = contactBodyOf(ship);

    // `state.velocity` is the *body origin's*. Whatever velocity the solver
    // carries, evaluating it back at the body origin has to reproduce that
    // exactly -- and it does not if the cog's velocity is taken to be the
    // origin's, which is the same class of error as reading state.velocity.x as
    // the ship's speed.
    expectTrue("the ship is spinning, so the two points differ",
               length(body.velocity - ship.state.velocity) > 0.05);
    expectNear("the solver's body reproduces the origin's velocity",
               length(body.velocityAt(ship.state.position) - ship.state.velocity), 0.0, 1e-12);

    // The inertia must be in world axes. The ship's own x axis is a principal axis
    // of her body-frame inertia, so spinning about its world image must give an
    // angular momentum along that same image. Left in the body frame it does not.
    const Vec3 bow = r * Vec3{1, 0, 0};
    const Vec3 momentum = body.inertia * bow;
    expectNear("angular momentum about the ship's own bow axis is along it",
               length(normalize(momentum) - bow), 0.0, 1e-9);
    const Vec3 unrotated = ship.massProperties().inertiaAboutCog * bow;
    expectTrue("and the body-frame tensor would not have said so",
               length(normalize(unrotated) - bow) > 1e-3);
}

void testTheExternalLoadIsConsumedExactlyOnce() {
    // The accumulator is a per-step channel, not a setting. A version that never
    // cleared it would push a ship for the rest of the run off one call, which
    // reads as a ship that will not stop rather than as a bug in contact.
    const HullParticulars p = coarseParticulars();
    const Sea sea(0.0);
    const double dt = 0.05;
    const Vec3 push{0, -3.0e7, 0};

    const auto run = [&](int applications) {
        Ship ship = plainShip(p);
        ship.initialise(sea);
        for (int step = 0; step < 20; ++step) {
            if (step < applications) {
                ship.externalForce += push;
                ship.externalMoment += Vec3{0, 0, 4.0e8};
            }
            ship.step(dt, sea);
        }
        return ship;
    };

    const Ship once = run(1);
    expectNear("the force accumulator is empty after a step", length(once.externalForce), 0.0, 0.0);
    expectNear("and so is the moment accumulator", length(once.externalMoment), 0.0, 0.0);

    const Ship never = run(0);
    const Ship twice = run(2);
    const double fromOne = std::abs(once.state.velocity.y - never.state.velocity.y);
    const double fromTwo = std::abs(twice.state.velocity.y - never.state.velocity.y);
    expectTrue("one application does move the ship", fromOne > 1e-3);
    expectTrue("two applications move her about twice as far, not twenty times",
               fromTwo > 1.8 * fromOne && fromTwo < 2.2 * fromOne);
    const double spinOne = std::abs(once.state.angularVelocity.z - never.state.angularVelocity.z);
    const double spinTwo = std::abs(twice.state.angularVelocity.z - never.state.angularVelocity.z);
    expectTrue("and the moment behaves the same way",
               spinOne > 1e-6 && spinTwo > 1.8 * spinOne && spinTwo < 2.2 * spinOne);
}

void testThePatchIsReportedInEachShipsOwnFrame() {
    // The body-frame patch is the number the structural model consumes, and the
    // transform into it is the identity on a ship that is upright and heading
    // along +x. Both ships here are at real attitudes, so a rotation used in place
    // of its transpose lands the load on the wrong side of the wrong ship.
    const HullParticulars p = coarseParticulars();
    const Sea sea(0.0);
    Ship a = plainShip(p);
    Ship b = plainShip(p);
    a.initialise(sea);
    b.initialise(sea);
    b.state.orientation = Quat::fromAxisAngle(Vec3{0, 0, 1}, 0.6);
    a.state.orientation = Quat::fromAxisAngle(Vec3{0, 0, 1}, 0.6 - kPiLocal / 2);
    const Mat3 rb = b.state.orientation.toMat3();
    a.state.position = rb * Vec3{30.0, 0.5 * p.beam + 0.5 * p.lengthPp - 6.0, 0.0};

    const HullContact c = shipContact(a, b);
    expectTrue("the two ships are in contact", c.touching);
    expectTrue("and both are at non-trivial headings", std::abs(rb(0, 0)) < 0.9);

    const Mat3 ra = a.state.orientation.toMat3();
    expectNear("the struck ship's patch maps back to the world point it came from",
               length(rb * c.onB.centroidBody + b.state.position - c.onB.centroid), 0.0, 1e-9);
    expectNear("and so does its normal",
               length(rb * c.onB.normalBody - c.onB.normal), 0.0, 1e-12);
    expectNear("the striker's patch too",
               length(ra * c.onA.centroidBody + a.state.position - c.onA.centroid), 0.0, 1e-9);
    expectNear("and its normal", length(ra * c.onA.normalBody - c.onA.normal), 0.0, 1e-12);
    // The guard: the body-frame answer must differ from the world one, or the
    // round trip above is the identity twice over.
    expectTrue("the two frames genuinely differ",
               length(c.onB.centroidBody - c.onB.centroid) > 1.0);
}


void testTwoShipsPassingCloseDoNotInteract() {
    const HullParticulars p = coarseParticulars();
    const Sea sea(0.0);

    // Two ships abeam with one beam of clear water between them, one overtaking.
    // Run them with the contact solver in the loop and again without, and require
    // the two tracks to be identical: contact must be *inert* when nothing
    // touches, not merely small.
    const auto run = [&](bool withContact) {
        Ship a = plainShip(p);
        Ship b = plainShip(p);
        a.initialise(sea);
        b.initialise(sea);
        b.state.position.y = 2.0 * p.beam;
        a.state.position.x = -60;
        a.state.velocity = {8.0, 0, 0};
        b.state.velocity = {2.0, 0, 0};
        const ContactMaterial material;
        struct Result { Vec3 track; double firedVolume = 0; double closestAbeam = 1e30; };
        Result result;
        for (int step = 0; step < 300; ++step) {
            if (withContact)
                result.firedVolume =
                    std::max(result.firedVolume, applyContact(a, b, material, 0.05).volume);
            a.step(0.05, sea);
            b.step(0.05, sea);
            result.closestAbeam =
                std::min(result.closestAbeam, std::abs(a.state.position.x - b.state.position.x));
        }
        result.track = a.state.position;
        return result;
    };

    const auto with = run(true);
    const auto without = run(false);
    expectNear("a beam of clear water never registers as contact", with.firedVolume, 0.0, 0.0);
    expectNear("and the ship's track is bit-for-bit what it was without the solver",
               length(with.track - without.track), 0.0, 0.0);
    // The guard against a vacuous negative: the overtaking ship has to actually
    // come abreast, or this is a test of two hulls that were never near each other.
    expectTrue("the overtaking ship did come abreast", with.closestAbeam < 5.0);
}

void testRammingAShip() {
    const HullParticulars p = coarseParticulars();
    const Sea sea(0.0);

    Ship struck = plainShip(p);
    Ship striker = plainShip(p);
    struck.initialise(sea);
    striker.initialise(sea);

    // The struck ship lies along +x at rest. The striker comes in on her port beam
    // at 6 m/s, aimed 30 m forward of midship: bow first, since a rotation of -90
    // degrees about z carries the striker's own +x onto world -y.
    const double impactX = 30.0;
    const double closing = 6.0;
    striker.state.orientation = Quat::fromAxisAngle({0, 0, 1}, -kPiLocal / 2);
    striker.state.position.x = impactX;
    striker.state.position.y = 0.5 * p.beam + 0.5 * p.lengthPp + 1.0;
    striker.state.velocity = {0, -closing, 0};

    const ContactMaterial material;
    ContactHistory history;
    const double dt = 0.01;

    std::string csv =
        "t,volume,depth,patchArea,normalForce,pressure,approach,workJ,bodyX,bodyY,bodyZ\n";
    double time = 0, firstTouch = -1;
    int touchingSteps = 0;
    for (int step = 0; step < 500; ++step) {
        const HullContact c = applyContact(striker, struck, material, dt, &history);
        if (c.touching) {
            ++touchingSteps;
            if (firstTouch < 0) firstTouch = time;
            const ContactLoad load =
                contactLoad(c, contactBodyOf(striker), contactBodyOf(struck), material);
            char line[512];
            std::snprintf(line, sizeof(line),
                          "%.3f,%.4f,%.4f,%.3f,%.6g,%.6g,%.4f,%.6g,%.3f,%.3f,%.3f\n", time,
                          c.volume, c.depth, c.patchArea, load.normalForce, load.pressure,
                          load.approachRate, history.work, c.onB.centroidBody.x,
                          c.onB.centroidBody.y, c.onB.centroidBody.z);
            csv += line;
        }
        striker.step(dt, sea);
        struck.step(dt, sea);
        time += dt;
        if (touchingSteps > 0 && !c.touching && time > firstTouch + 0.5) break;
    }

    const std::string path = testing::scratchDir() + "collision_ram.csv";
    if (std::FILE* f = std::fopen(path.c_str(), "w")) {
        std::fwrite(csv.data(), 1, csv.size(), f);
        std::fclose(f);
    }

    std::printf("     ram: closing %.1f m/s, %d steps in contact (%.2f s), peak %.1f MN,\n"
                "          %.3f m penetration over %.1f m2 at %.2f MPa, %.1f MJ absorbed,\n"
                "          patch at x = %+.1f m, y = %+.1f m, z = %+.1f m in the struck ship\n",
                closing, touchingSteps, history.duration, history.peakForce * 1e-6,
                history.atPeak.depth, history.atPeak.patchArea, history.loadAtPeak.pressure * 1e-6,
                history.work * 1e-6, history.atPeak.onB.centroidBody.x,
                history.atPeak.onB.centroidBody.y, history.atPeak.onB.centroidBody.z);

    expectTrue("the ram made contact", touchingSteps > 10);
    expectTrue("over many steps, so the force history is resolved rather than impulsive",
               history.duration > 0.15);
    expectTrue("the peak force is in the tens to hundreds of meganewtons",
               history.peakForce > 2e7 && history.peakForce < 5e9);
    expectTrue("energy was absorbed", history.work > 1e6);

    // The pair cannot give up more energy than it brought. **What it brought is
    // not the rigid closing energy**: each ship must also accelerate the water it
    // shoves aside, and the contact removes that entrained momentum too. The
    // striker decelerates along its own bow, so it carries surge added mass; the
    // struck ship is driven along her beam and carries sway added mass, which is
    // nearly her own displacement again. At floating equilibrium rho times the
    // submerged volume *is* the ship's mass, so the coefficients multiply it
    // directly.
    //
    // Bounding against the rigid masses instead gives 223 MJ against a measured
    // 233 MJ, and the model reads as though it were creating energy when it is
    // only accounting for the water. That version passed, by half a percent.
    const double massA = striker.massProperties().mass;
    const double massB = struck.massProperties().mass;
    const double effectiveA = massA * (1.0 + striker.addedMassSurge);
    const double effectiveB = massB * (1.0 + struck.addedMassSway);
    const double reduced = effectiveA * effectiveB / (effectiveA + effectiveB);
    expectTrue("and not more than the pair brought with them, entrained water included",
               history.work < 0.5 * reduced * closing * closing);
    expectTrue("but more than the rigid masses alone could have supplied, which is"
               " the entrained water showing up",
               history.work > 0.5 * (massA * massB / (massA + massB)) * closing * closing);

    // Where the load landed, in the struck ship's own frame: on her port side,
    // near the aimed station, and inside her depth. This is the number the
    // structural model consumes, so it is the number worth asserting.
    expectNear("the patch is at the station that was aimed at",
               history.atPeak.onB.centroidBody.x, impactX, 12.0);
    expectTrue("on the port side", history.atPeak.onB.centroidBody.y > 0.3 * p.beam);
    expectTrue("and inside the hull's depth", history.atPeak.onB.centroidBody.z > 0 &&
                                                  history.atPeak.onB.centroidBody.z < p.depth);
    expectTrue("the patch normal points out of the struck ship's port side",
               history.atPeak.onB.normalBody.y > 0.5);
    expectTrue("the two surfaces faced each other", history.atPeak.normalAgreement > 0.5);
    expectTrue("and it was one contact region, not two",
               history.atPeak.patchSeparation < 5.0 * history.atPeak.depth + 2.0);

    // What the collision did to the ships. The struck ship is driven to starboard
    // and yaws bow-to-starboard, because the blow landed forward of her cog.
    expectTrue("the struck ship is driven to starboard", struck.state.velocity.y < -0.05);
    expectTrue("and yaws bow-to-starboard from a blow forward of her cog",
               struck.state.angularVelocity.z < 0);
    expectTrue("the striker is slowed", striker.state.velocity.y > -closing + 0.05);
    expectTrue("neither ship went numerically mad",
               std::isfinite(struck.state.position.z) && std::isfinite(striker.state.position.z) &&
                   std::abs(struck.state.position.z) < 50.0);
}

}  // namespace

void runCollisionTests() {
    std::printf("\n--- collision ---\n");
    testSolidMomentsOfABox();
    testSymmetricEigenSolvesAKnownMatrix();
    testPrincipalBoxRecoversABox();
    testTwoBoxesOverlapExactly();
    testOverlapIsSymmetricInItsArguments();
    testARotatedBoxGivesTheOctagon();
    testContainmentIsExact();
    testANonConvexSolidNeedsNegativeTetrahedra();
    testHullsTouchingExactlyDoNotOverlap();
    testAnInsideOutHullIsRefusedRatherThanMisintegrated();
    testTheAnswerDoesNotDependOnTheOverlapBoxRounds();
    testSphereLensConverges();
    testHullsApartDoNotTouchAndHullsTogetherDo();
    testImpulseConservesMomentumAndAngularMomentum();
    testAnImpulseOnTheWrongLeverBreaksAngularMomentum();
    testRestitutionSetsTheEnergyLossExactly();
    testEffectiveMassIsTheReducedMass();
    testHeadOnIsSymmetricAndYawFree();
    testGlancingBlowYawsTheRightWay();
    testPenaltyContactIsTheOscillatorItClaimsToBe();
    testPenaltyRestitutionConvergesWithTheTimestep();
    testDissipationSetsRestitution();
    testFrictionRemovesEnergyWithoutBreakingConservation();
    testAContactNeverPulls();
    testTheLoadActsAtTheCentreOfPressure();
    testContactBodyIsReferredToTheCentreOfGravity();
    testTheExternalLoadIsConsumedExactlyOnce();
    testThePatchIsReportedInEachShipsOwnFrame();
    testTwoShipsPassingCloseDoNotInteract();
    testRammingAShip();
}
