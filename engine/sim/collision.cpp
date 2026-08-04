// SPDX-License-Identifier: MIT
#include "collision.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace sim {
namespace {

constexpr double kInf = std::numeric_limits<double>::infinity();

Mat3 outer(const Vec3& a, const Vec3& b) {
    Mat3 m;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) m(i, j) = a[i] * b[j];
    return m;
}

// --- Axis-aligned bounds ------------------------------------------------------

struct Aabb {
    Vec3 lo{kInf, kInf, kInf};
    Vec3 hi{-kInf, -kInf, -kInf};

    void add(const Vec3& p) {
        lo = {std::min(lo.x, p.x), std::min(lo.y, p.y), std::min(lo.z, p.z)};
        hi = {std::max(hi.x, p.x), std::max(hi.y, p.y), std::max(hi.z, p.z)};
    }
    bool empty() const { return !(lo.x <= hi.x && lo.y <= hi.y && lo.z <= hi.z); }
    Vec3 centre() const { return (lo + hi) * 0.5; }
    double diagonal() const { return empty() ? 0.0 : length(hi - lo); }
    Aabb intersected(const Aabb& o) const {
        Aabb r;
        r.lo = {std::max(lo.x, o.lo.x), std::max(lo.y, o.lo.y), std::max(lo.z, o.lo.z)};
        r.hi = {std::min(hi.x, o.hi.x), std::min(hi.y, o.hi.y), std::min(hi.z, o.hi.z)};
        return r;
    }
    Aabb padded(double d) const {
        Aabb r;
        r.lo = lo - Vec3{d, d, d};
        r.hi = hi + Vec3{d, d, d};
        return r;
    }
};

Aabb boundsOf(const TriMesh& m) {
    Aabb b;
    for (const Vec3& v : m.verts) b.add(v);
    return b;
}

// --- Convex polygon clipping --------------------------------------------------

// Sutherland-Hodgman against dot(n, x) <= offset. A convex k-gon clipped by a
// plane has at most k + 1 vertices, so a triangle through four planes never
// exceeds seven and the eight-slot buffers below are never overrun.
int clipPolygon(const Vec3* in, int n, Vec3* out, const Vec3& normal, double offset) {
    int m = 0;
    for (int i = 0; i < n; ++i) {
        const Vec3& p = in[i];
        const Vec3& q = in[(i + 1) % n];
        const double dp = dot(normal, p) - offset;
        const double dq = dot(normal, q) - offset;
        if (dp <= 0) out[m++] = p;
        if ((dp < 0 && dq > 0) || (dp > 0 && dq < 0))
            out[m++] = p + (q - p) * (dp / (dp - dq));
    }
    return m;
}

// --- Accumulators -------------------------------------------------------------
//
// Both take polygons already expressed *relative to the reference point*, which
// keeps the cancellation in the flux sum local to the contact rather than
// proportional to the ships' world positions.

// The divergence-theorem flux of a piece of the overlap solid's boundary. The
// boundary of A ∩ B is (∂A inside B) plus (∂B inside A) and nothing else, so
// summing this over both gives the overlap's own volume and moments exactly.
struct FluxAccumulator {
    double volume = 0;
    Vec3   first{};
    // **Not `Mat3 second{}`.** Mat3 carries a default member initialiser of the
    // identity, so a value-initialised one is I and not 0 -- an accumulator
    // started that way reports every solid's second moment one too large on the
    // diagonal, which for a 2 x 4 x 6 box is a 6% error in the shortest extent and
    // reads as a slightly fat box rather than as a broken integral.
    Mat3   second = Mat3::zero();

    void addPolygon(const Vec3* p, int n, double weight) {
        for (int i = 1; i + 1 < n; ++i) {
            const Vec3& a = p[0];
            const Vec3& b = p[i];
            const Vec3& c = p[i + 1];
            // Signed volume of the tetrahedron (reference, a, b, c).
            const double v = dot(a, cross(b, c)) / 6.0;
            if (v == 0.0) continue;
            const double w = weight * v;
            volume += w;
            const Vec3 s = a + b + c;
            first += s * (0.25 * w);
            // Second moment of a simplex: with one vertex at the origin,
            // integral of x (x) x = V/20 * (sum p_k (x) p_k + S (x) S).
            const Mat3 m = outer(a, a) + outer(b, b) + outer(c, c) + outer(s, s);
            second = second + m * (w / 20.0);
        }
    }
};

// One hull's surface inside the other: the contact patch.
struct PatchAccumulator {
    double area = 0;
    Vec3   moment{};      // area-weighted centroid numerator, reference-relative
    Vec3   normalMoment{};  // area-weighted outward normal

    void addPolygon(const Vec3* p, int n, const Vec3& unitNormal, double weight) {
        for (int i = 1; i + 1 < n; ++i) {
            const Vec3& a = p[0];
            const Vec3& b = p[i];
            const Vec3& c = p[i + 1];
            // Signed against the parent triangle's own normal, so a fragment that
            // the clip left wound the other way cannot inflate the area.
            const double da = 0.5 * dot(cross(b - a, c - a), unitNormal) * weight;
            if (da == 0.0) continue;
            area += da;
            moment += (a + b + c) * (da / 3.0);
            normalMoment += unitNormal * da;
        }
    }
};

// --- The overlap sweep --------------------------------------------------------

struct SurfaceTri {
    Vec3 v[3];        // reference-relative
    Vec3 unitNormal;  // outward
    Vec3 lo, hi;      // reference-relative bounds
};

void buildSurface(const TriMesh& mesh, const Vec3& reference, std::vector<SurfaceTri>& out) {
    out.clear();
    out.reserve(mesh.tris.size());
    for (const Tri& t : mesh.tris) {
        SurfaceTri s;
        s.v[0] = mesh.verts[t.a] - reference;
        s.v[1] = mesh.verts[t.b] - reference;
        s.v[2] = mesh.verts[t.c] - reference;
        const Vec3 areaVector = cross(s.v[1] - s.v[0], s.v[2] - s.v[0]);
        const double twiceArea = length(areaVector);
        if (twiceArea <= 0) continue;  // degenerate triangle carries no flux
        s.unitNormal = areaVector / twiceArea;
        s.lo = {std::min({s.v[0].x, s.v[1].x, s.v[2].x}),
                std::min({s.v[0].y, s.v[1].y, s.v[2].y}),
                std::min({s.v[0].z, s.v[1].z, s.v[2].z})};
        s.hi = {std::max({s.v[0].x, s.v[1].x, s.v[2].x}),
                std::max({s.v[0].y, s.v[1].y, s.v[2].y}),
                std::max({s.v[0].z, s.v[1].z, s.v[2].z})};
        out.push_back(s);
    }
}

// Accumulate "the part of `surface` that lies inside `solid`".
//
// `solid` is decomposed into signed tetrahedra sharing the origin -- which is the
// reference point, since everything here is reference-relative. That decomposition
// reproduces the solid's indicator function pointwise, so clipping each surface
// triangle against each tetrahedron and summing with the tetrahedron's sign gives
// the exact intersection of the triangle with the solid. Four half-spaces per
// tetrahedron, three of which pass through the origin.
//
// **Coincident faces get half weight.** The boundary of A ∩ B is (∂A inside B)
// plus (∂B inside A), and those two sets are disjoint *except* where a face of one
// hull lies exactly in a face of the other -- which both sweeps then claim, and the
// flux of that face is counted twice. It is not an exotic case: two boxes meeting
// face to face share four side planes, and the overlap volume comes out 5/3 too
// large with no other symptom, because everything else about the answer stays
// plausible. Halving the weight in *both* sweeps counts the shared face exactly
// once and keeps the routine symmetric in its two arguments. Faces that coincide
// with *opposing* normals are two hulls touching from outside, which bounds no
// overlap at all, and get weight zero.
void sweepSurfaceInSolid(const std::vector<SurfaceTri>& surface, const TriMesh& solid,
                         const Vec3& reference, double coplanarTolerance, FluxAccumulator& flux,
                         PatchAccumulator& patch, double& opposedCoincidentArea) {
    for (const Tri& t : solid.tris) {
        const Vec3 a = solid.verts[t.a] - reference;
        const Vec3 b = solid.verts[t.b] - reference;
        const Vec3 c = solid.verts[t.c] - reference;

        const double six = dot(a, cross(b, c));
        // A tetrahedron whose apex is coplanar with its base has zero volume, so
        // its indicator is zero almost everywhere and dropping it is exact. It has
        // to be dropped rather than clipped against, because its face normals are
        // then degenerate and a degenerate plane clips nothing.
        if (six * six <= 1e-24 * length2(a) * length2(b) * length2(c)) continue;
        const double sign = six > 0 ? 1.0 : -1.0;

        Vec3 planeNormal[4];
        double planeOffset[4];
        planeNormal[0] = cross(a, b);
        planeOffset[0] = 0;
        if (dot(planeNormal[0], c) > 0) planeNormal[0] = -planeNormal[0];
        planeNormal[1] = cross(b, c);
        planeOffset[1] = 0;
        if (dot(planeNormal[1], a) > 0) planeNormal[1] = -planeNormal[1];
        planeNormal[2] = cross(c, a);
        planeOffset[2] = 0;
        if (dot(planeNormal[2], b) > 0) planeNormal[2] = -planeNormal[2];
        // The solid's own outward normal at this triangle, kept separately: the
        // base half-space below is oriented away from the apex, which is the same
        // direction only when the apex lies inside the solid.
        const Vec3 baseArea = cross(b - a, c - a);
        const double baseAreaLength = length(baseArea);
        planeNormal[3] = baseArea;
        if (dot(planeNormal[3], a) < 0) planeNormal[3] = -planeNormal[3];
        planeOffset[3] = dot(planeNormal[3], a);

        const Vec3 tetLo{std::min({0.0, a.x, b.x, c.x}), std::min({0.0, a.y, b.y, c.y}),
                         std::min({0.0, a.z, b.z, c.z})};
        const Vec3 tetHi{std::max({0.0, a.x, b.x, c.x}), std::max({0.0, a.y, b.y, c.y}),
                         std::max({0.0, a.z, b.z, c.z})};

        for (const SurfaceTri& s : surface) {
            if (s.hi.x < tetLo.x || s.lo.x > tetHi.x) continue;
            if (s.hi.y < tetLo.y || s.lo.y > tetHi.y) continue;
            if (s.hi.z < tetLo.z || s.lo.z > tetHi.z) continue;

            Vec3 buffer[2][8];
            buffer[0][0] = s.v[0];
            buffer[0][1] = s.v[1];
            buffer[0][2] = s.v[2];
            int n = 3;
            int cur = 0;
            for (int k = 0; k < 4; ++k) {
                n = clipPolygon(buffer[cur], n, buffer[cur ^ 1], planeNormal[k], planeOffset[k]);
                cur ^= 1;
                if (n < 3) break;
            }
            if (n < 3) continue;

            double weight = sign;
            if (baseAreaLength > 0) {
                double furthest = 0;
                for (const Vec3& v : s.v)
                    furthest = std::max(furthest, std::abs(dot(baseArea, v - a)) / baseAreaLength);
                if (furthest <= coplanarTolerance) {
                    const bool aligned = dot(s.unitNormal, baseArea) > 0;
                    if (!aligned) {
                        // Measure what was refused, so hullContact() can say that
                        // this case arose at all. See the note on its limits there.
                        Vec3 area{};
                        for (int i = 1; i + 1 < n; ++i)
                            area += cross(buffer[cur][i] - buffer[cur][0],
                                          buffer[cur][i + 1] - buffer[cur][0]);
                        opposedCoincidentArea += 0.5 * length(area);
                    }
                    weight *= aligned ? 0.5 : 0.0;
                }
            }
            if (weight == 0.0) continue;

            flux.addPolygon(buffer[cur], n, weight);
            patch.addPolygon(buffer[cur], n, s.unitNormal, weight);
        }
    }
}

// A box clip that can be trusted: non-empty, positively oriented, and still a
// closed manifold. `clipByPlane` caps its cuts by ear clipping, which assumes the
// intersection loops are simple and non-nested; when that assumption fails the
// result integrates to a plausible wrong volume rather than to an obvious one.
bool usableClip(const TriMesh& mesh, bool verify, double& volume) {
    volume = 0;
    if (mesh.tris.empty()) return false;
    volume = integrate(mesh).volume;
    if (!(volume > 0)) return false;
    return !verify || isClosedManifold(mesh);
}

}  // namespace

// --- Moments ------------------------------------------------------------------

Vec3 SolidMoments::centroid() const {
    return volume != 0 ? reference + first / volume : reference;
}

Mat3 SolidMoments::covariance() const {
    if (!(volume > 0)) return Mat3::zero();
    const Vec3 c = first / volume;
    return second * (1.0 / volume) + outer(c, c) * -1.0;
}

SolidMoments solidMoments(const TriMesh& mesh, const Vec3& reference) {
    FluxAccumulator flux;
    Vec3 poly[3];
    for (const Tri& t : mesh.tris) {
        poly[0] = mesh.verts[t.a] - reference;
        poly[1] = mesh.verts[t.b] - reference;
        poly[2] = mesh.verts[t.c] - reference;
        flux.addPolygon(poly, 3, 1.0);
    }
    SolidMoments m;
    m.volume = flux.volume;
    m.first = flux.first;
    m.second = flux.second;
    m.reference = reference;
    return m;
}

void symmetricEigen(const Mat3& in, std::array<double, 3>& values,
                    std::array<Vec3, 3>& vectors) {
    double a[3][3];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) a[i][j] = 0.5 * (in(i, j) + in(j, i));

    double v[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    const double scale = std::max({std::abs(a[0][0]), std::abs(a[1][1]), std::abs(a[2][2]), 1e-300});

    // Cyclic Jacobi. Twenty sweeps is far past convergence for 3x3; the loop
    // exists only so that a pathological input terminates.
    for (int sweep = 0; sweep < 20; ++sweep) {
        const double off = a[0][1] * a[0][1] + a[0][2] * a[0][2] + a[1][2] * a[1][2];
        if (off <= 1e-32 * scale * scale) break;
        for (int p = 0; p < 2; ++p)
            for (int q = p + 1; q < 3; ++q) {
                if (std::abs(a[p][q]) <= 1e-18 * scale) continue;
                const double theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
                const double sign = theta >= 0 ? 1.0 : -1.0;
                const double tangent = sign / (std::abs(theta) + std::sqrt(theta * theta + 1.0));
                const double cosine = 1.0 / std::sqrt(tangent * tangent + 1.0);
                const double sine = tangent * cosine;
                for (int k = 0; k < 3; ++k) {
                    const double akp = a[k][p], akq = a[k][q];
                    a[k][p] = cosine * akp - sine * akq;
                    a[k][q] = sine * akp + cosine * akq;
                }
                for (int k = 0; k < 3; ++k) {
                    const double apk = a[p][k], aqk = a[q][k];
                    a[p][k] = cosine * apk - sine * aqk;
                    a[q][k] = sine * apk + cosine * aqk;
                }
                for (int k = 0; k < 3; ++k) {
                    const double vkp = v[k][p], vkq = v[k][q];
                    v[k][p] = cosine * vkp - sine * vkq;
                    v[k][q] = sine * vkp + cosine * vkq;
                }
            }
    }

    int order[3] = {0, 1, 2};
    std::sort(order, order + 3, [&](int i, int j) { return a[i][i] > a[j][j]; });
    for (int i = 0; i < 3; ++i) {
        values[static_cast<std::size_t>(i)] = a[order[i]][order[i]];
        vectors[static_cast<std::size_t>(i)] =
            normalize(Vec3{v[0][order[i]], v[1][order[i]], v[2][order[i]]});
    }
}

PrincipalBox principalBox(const SolidMoments& moments) {
    PrincipalBox box;
    box.axis = {Vec3{1, 0, 0}, Vec3{0, 1, 0}, Vec3{0, 0, 1}};
    if (!(moments.volume > 0)) return box;
    box.centre = moments.centroid();

    std::array<double, 3> lambda{};
    symmetricEigen(moments.covariance(), lambda, box.axis);
    // A uniform box of side d has covariance d^2 / 12 about that axis, so this
    // inverts exactly for a box and reads as an equivalent thickness otherwise.
    for (int i = 0; i < 3; ++i)
        box.extent[static_cast<std::size_t>(i)] =
            std::sqrt(std::max(0.0, 12.0 * lambda[static_cast<std::size_t>(i)]));
    // Right-handed, so the axes can be used as a frame without a surprise.
    if (dot(cross(box.axis[0], box.axis[1]), box.axis[2]) < 0) box.axis[2] = -box.axis[2];
    return box;
}

// --- Detection ----------------------------------------------------------------

HullContact hullContact(const TriMesh& worldA, const TriMesh& worldB,
                        const ContactParams& params) {
    HullContact out;

    // Broad phase. Two ships passing a beam apart fail here and cost nothing.
    Aabb box = boundsOf(worldA).intersected(boundsOf(worldB));
    if (box.empty()) return out;

    const auto pad = [&](const Aabb& b) {
        return b.padded(std::max(params.boxPadFloor, params.boxPadFraction * b.diagonal()));
    };

    // Cut both hulls down to the overlap region. Legitimate because A ∩ B lies
    // inside every box used here, so (A ∩ box) ∩ (B ∩ box) is still A ∩ B; and
    // because the box is *padded*, A ∩ B is strictly interior to it and the cut
    // faces the clip introduces contribute nothing to either the volume or the
    // contact patches.
    double volumeA = 0, volumeB = 0;
    TriMesh a = clipToBox(worldA, pad(box).lo, pad(box).hi);
    TriMesh b = clipToBox(worldB, pad(box).lo, pad(box).hi);
    if (!usableClip(a, params.verifyClips, volumeA) ||
        !usableClip(b, params.verifyClips, volumeB)) {
        // Fall back to the whole hull: slower by the ratio of triangle counts, and
        // right, which is the correct way round for a fallback.
        if (a.tris.empty() || b.tris.empty()) return out;
        out.problems.push_back("overlap box clip is not a closed solid; using whole hulls");
        a = worldA;
        b = worldB;
        volumeA = integrate(a).volume;
        volumeB = integrate(b).volume;
    } else {
        for (int round = 1; round < params.boxRounds; ++round) {
            const Aabb next = boundsOf(a).intersected(boundsOf(b));
            if (next.empty()) return out;
            if (next.diagonal() >= 0.95 * box.diagonal()) break;
            box = next;
            TriMesh a2 = clipToBox(a, pad(box).lo, pad(box).hi);
            TriMesh b2 = clipToBox(b, pad(box).lo, pad(box).hi);
            double v2a = 0, v2b = 0;
            if (!usableClip(a2, params.verifyClips, v2a) ||
                !usableClip(b2, params.verifyClips, v2b))
                break;
            a = std::move(a2);
            b = std::move(b2);
            volumeA = v2a;
            volumeB = v2b;
        }
    }

    const Vec3 reference = box.centre();
    // Coincident faces are recognised relative to the size of the region being
    // integrated, not absolutely: the meshes arrive through a rigid transform, so
    // two planes that are meant to be the same plane agree to a part in 1e15.
    const double coplanarTolerance = 1e-9 * std::max(box.diagonal(), 1e-6);

    std::vector<SurfaceTri> surface;
    FluxAccumulator flux;
    PatchAccumulator patchA, patchB;

    double opposedCoincidentArea = 0;
    buildSurface(a, reference, surface);
    sweepSurfaceInSolid(surface, b, reference, coplanarTolerance, flux, patchA,
                        opposedCoincidentArea);
    buildSurface(b, reference, surface);
    sweepSurfaceInSolid(surface, a, reference, coplanarTolerance, flux, patchB,
                        opposedCoincidentArea);

    out.volume = flux.volume;
    // Two hulls with faces coincident in *opposing* directions are touching from
    // outside, and the convention refuses those faces -- which is exactly right
    // when the coincident tetrahedron is the only one that claims them, as it is
    // for two convex hulls laid flush. It is not right in general: a face lying
    // on a *non-convex* solid's own skin can also be claimed by a tetrahedron it
    // is not coplanar with, and that claim survives the refusal. Deciding
    // membership per point on a coincident face is the coincident-face problem
    // every mesh boolean has, and it is not solved here -- so the case is
    // reported instead of guessed at. It needs exact coplanarity, to a part in
    // 1e9 of the contact region, which two independently placed hulls do not
    // produce.
    if (opposedCoincidentArea > 0)
        out.problems.push_back("hulls have faces coincident in opposing directions; the overlap "
                               "there is a convention, not a measurement");
    // An intersection cannot be larger than either of the things intersected. This
    // costs one integral of each clipped chunk -- already paid for above -- and it
    // is the check that would have caught the coincident-face double count on its
    // first run rather than on a closed-form test that happened to be pointed at
    // it. Reported rather than clamped: a wrong answer that says so is a different
    // thing from a wrong answer that does not.
    if (flux.volume > 1.000001 * std::min(volumeA, volumeB) + 1e-9)
        out.problems.push_back("overlap exceeds the volume of the hulls it lies in");
    if (!(flux.volume > params.minVolume)) return out;
    out.touching = true;

    SolidMoments moments;
    moments.volume = flux.volume;
    moments.first = flux.first;
    moments.second = flux.second;
    moments.reference = reference;
    out.centroid = moments.centroid();

    const PrincipalBox shape = principalBox(moments);
    out.extent = {shape.extent[0], shape.extent[1], shape.extent[2]};
    out.depth = shape.extent[2];
    out.patchArea = out.depth > 0 ? out.volume / out.depth : 0.0;

    out.onA.area = patchA.area;
    out.onB.area = patchB.area;
    if (patchA.area > 0) out.onA.centroid = reference + patchA.moment / patchA.area;
    if (patchB.area > 0) out.onB.centroid = reference + patchB.moment / patchB.area;
    out.onA.normal = normalize(patchA.normalMoment);
    out.onB.normal = normalize(patchB.normalMoment);
    out.normalAgreement = -dot(out.onA.normal, out.onB.normal);
    out.patchSeparation = length(out.onA.centroid - out.onB.centroid);

    // The contact normal, out of B and into A. The two patches' area-weighted
    // normals are the surfaces' own and are nearly antiparallel for any contact
    // this model claims to represent, so their difference is both directions'
    // evidence at once. Where they cancel -- a grazing touch, or an overlap deep
    // enough to see both sides of the struck hull -- the overlap solid's own
    // thinnest principal axis is the fallback, signed so that pushing A along it
    // separates the two hull centres.
    const Vec3 opposed = patchB.normalMoment - patchA.normalMoment;
    const double patchScale = std::max(patchA.area, patchB.area);
    if (length(opposed) > 1e-6 * std::max(patchScale, 1e-12)) {
        out.normal = normalize(opposed);
    } else {
        const Vec3 separation = boundsOf(worldA).centre() - boundsOf(worldB).centre();
        out.normal = shape.axis[2];
        if (dot(out.normal, separation) < 0) out.normal = -out.normal;
        out.problems.push_back("contact patches gave no usable normal; used the overlap's own axis");
    }
    return out;
}

HullContact shipContact(const Ship& a, const Ship& b, const ContactParams& params) {
    const Mat3 ra = a.state.orientation.toMat3();
    const Mat3 rb = b.state.orientation.toMat3();
    const TriMesh worldA = transformed(a.hull, ra, a.state.position);
    const TriMesh worldB = transformed(b.hull, rb, b.state.position);

    HullContact contact = hullContact(worldA, worldB, params);
    if (!contact.touching) return contact;

    const Mat3 raT = ra.transposed();
    const Mat3 rbT = rb.transposed();
    contact.onA.centroidBody = raT * (contact.onA.centroid - a.state.position);
    contact.onA.normalBody = raT * contact.onA.normal;
    contact.onB.centroidBody = rbT * (contact.onB.centroid - b.state.position);
    contact.onB.normalBody = rbT * contact.onB.normal;
    return contact;
}

// --- Response -----------------------------------------------------------------

Vec3 ContactBody::velocityAt(const Vec3& worldPoint) const {
    return velocity + cross(angularVelocity, worldPoint - cog);
}

double ContactBody::kineticEnergy() const {
    return 0.5 * mass * length2(velocity) + 0.5 * dot(angularVelocity, inertia * angularVelocity);
}

Vec3 ContactBody::momentum() const { return velocity * mass; }

Vec3 ContactBody::angularMomentum(const Vec3& about) const {
    return inertia * angularVelocity + cross(cog - about, velocity * mass);
}

ImpulseSolution normalImpulse(const ContactBody& a, const ContactBody& b, const Vec3& point,
                              const Vec3& normal, double restitution) {
    ImpulseSolution solution;
    const Vec3 n = normalize(normal);
    if (length2(n) <= 0) return solution;

    const Vec3 ra = point - a.cog;
    const Vec3 rb = point - b.cog;
    const Vec3 relative = a.velocityAt(point) - b.velocityAt(point);
    // The normal points out of B into A, so A closing on B moves along -n.
    solution.approachSpeed = -dot(relative, n);

    // A body with non-positive mass and a singular inertia is immovable: both
    // reciprocals below vanish, which is exactly the infinite-mass limit.
    double inverseMass = 0;
    if (a.mass > 0) inverseMass += 1.0 / a.mass;
    if (b.mass > 0) inverseMass += 1.0 / b.mass;
    inverseMass += dot(n, cross(inverse(a.inertia) * cross(ra, n), ra));
    inverseMass += dot(n, cross(inverse(b.inertia) * cross(rb, n), rb));
    if (!(inverseMass > 0)) return solution;
    solution.effectiveMass = 1.0 / inverseMass;

    if (solution.approachSpeed <= 0) return solution;  // already separating
    const double e = std::clamp(restitution, 0.0, 1.0);
    solution.impulse = n * ((1.0 + e) * solution.approachSpeed * solution.effectiveMass);
    solution.energyLost = 0.5 * (1.0 - e * e) * solution.effectiveMass * solution.approachSpeed *
                          solution.approachSpeed;
    return solution;
}

void applyImpulse(ContactBody& a, ContactBody& b, const Vec3& point, const Vec3& impulse) {
    if (a.mass > 0) a.velocity += impulse / a.mass;
    if (b.mass > 0) b.velocity -= impulse / b.mass;
    a.angularVelocity += inverse(a.inertia) * cross(point - a.cog, impulse);
    b.angularVelocity -= inverse(b.inertia) * cross(point - b.cog, impulse);
}

ContactLoad contactLoad(const HullContact& contact, const ContactBody& a, const ContactBody& b,
                        const ContactMaterial& material) {
    ContactLoad load;
    if (!contact.touching || !(contact.volume > 0)) return load;

    // The overlap solid's centroid is the depth-weighted centroid of the contact
    // patch, which for a pressure proportional to local penetration is the centre
    // of pressure. Applying the resultant anywhere else would put the right force
    // on the wrong lever.
    load.point = contact.centroid;
    const Vec3 n = contact.normal;
    const Vec3 relative = a.velocityAt(load.point) - b.velocityAt(load.point);
    load.approachRate = -dot(relative, n);

    double normalForce = material.stiffness * contact.volume *
                         (1.0 + material.dissipation * load.approachRate);
    if (normalForce < 0) normalForce = 0;  // never pull the hulls together
    load.normalForce = normalForce;
    load.pressure = contact.patchArea > 0 ? normalForce / contact.patchArea : 0.0;
    load.slip = relative - n * dot(relative, n);

    Vec3 force = n * normalForce;
    if (material.friction > 0 && normalForce > 0) {
        const double slipSpeed = length(load.slip);
        if (slipSpeed > 0)
            force -= load.slip * (material.friction * normalForce /
                                  std::max(slipSpeed, material.frictionSpeed));
    }
    load.force = force;
    return load;
}

void ContactHistory::accumulate(const HullContact& contact, const ContactLoad& load,
                                const ContactBody& a, const ContactBody& b, double dt) {
    if (!contact.touching || !(dt > 0)) return;
    duration += dt;
    impulse += load.force * dt;
    // The pair feels +force on A and -force on B at one shared point, so the rate
    // at which their combined kinetic energy changes is +force . relativeVelocity.
    // Work is reported as energy *removed*, which at separation -- when the spring
    // holds nothing -- is the energy the structure had to absorb.
    const Vec3 relative = a.velocityAt(load.point) - b.velocityAt(load.point);
    work += -dot(load.force, relative) * dt;
    if (load.normalForce > peakForce) {
        peakForce = load.normalForce;
        atPeak = contact;
        loadAtPeak = load;
    }
}

// --- Ships --------------------------------------------------------------------

ContactBody contactBodyOf(const Ship& ship) {
    const Ship::MassProperties mp = ship.massProperties();
    const Mat3 r = ship.state.orientation.toMat3();
    ContactBody body;
    body.mass = mp.mass;
    body.inertia = r * mp.inertiaAboutCog * r.transposed();
    body.cog = r * mp.cog + ship.state.position;
    body.angularVelocity = ship.state.angularVelocity;
    // state.velocity is the *body origin's*. The contact solver takes moments
    // about the centre of gravity and needs that point's velocity.
    body.velocity = ship.state.velocity + cross(ship.state.angularVelocity, r * mp.cog);
    return body;
}

HullContact applyContact(Ship& a, Ship& b, const ContactMaterial& material, double dt,
                         ContactHistory* history, const ContactParams& params) {
    const HullContact contact = shipContact(a, b, params);
    if (!contact.touching) return contact;

    const ContactBody bodyA = contactBodyOf(a);
    const ContactBody bodyB = contactBodyOf(b);
    const ContactLoad load = contactLoad(contact, bodyA, bodyB, material);

    a.externalForce += load.force;
    a.externalMoment += cross(load.point - bodyA.cog, load.force);
    b.externalForce -= load.force;
    b.externalMoment -= cross(load.point - bodyB.cog, load.force);

    if (history) history->accumulate(contact, load, bodyA, bodyB, dt);
    return contact;
}

}  // namespace sim
