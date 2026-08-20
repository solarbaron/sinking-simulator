// SPDX-License-Identifier: MIT
// Closed-triangle-mesh volume integration with half-space clipping.
//
// Everything hydrostatic in this engine reduces to one question: given a closed
// mesh and a plane, what is the volume and centroid of the part below the plane?
// Buoyancy asks it of the hull against the sea surface. Floodwater asks it of a
// compartment against the internal free surface. Both use integrateBelowPlane().
#pragma once

#include "math.hpp"

#include <cstdint>
#include <string>
#include <vector>

namespace sim {

struct Tri {
    std::uint32_t a, b, c;
};

struct TriMesh {
    std::vector<Vec3> verts;
    std::vector<Tri>  tris;

    void append(const TriMesh& other);
};

// Result of integrating a (sub)region of a closed mesh.
struct VolumeIntegral {
    double volume = 0;    // m^3
    Vec3   centroid{};    // m, valid only when volume > 0
};

// Full closed-mesh volume and centroid. Assumes outward-facing winding.
VolumeIntegral integrate(const TriMesh& mesh);

// Volume and centroid of the region satisfying dot(n, x) <= planeOffset.
//
// No cap geometry is constructed: the divergence-theorem integral is taken about
// a reference point lying *on* the cutting plane, so every cap tetrahedron is
// degenerate and contributes exactly zero. This makes the routine correct for
// meshes of any topology, including ones the plane cuts into several loops.
VolumeIntegral integrateBelowPlane(const TriMesh& mesh, const Vec3& n, double planeOffset);

// A mesh bound to one plane *direction*, with every vertex's projection onto that
// direction computed once.
//
// The free-surface solve clips the same mesh against a dozen parallel planes to
// find one water level, and then does it again next tick, and the plane normal is
// the same for every compartment in the ship because they all share the same
// gravity. Under those conditions the per-triangle dot products dominate --
// a vertex shared by six triangles gets projected twelve times per clip. Hoisting
// them out is the single largest win available in the hydrostatic inner loop.
class PlaneSweep {
public:
    PlaneSweep(const TriMesh& mesh, const Vec3& n);
    // **A temporary is a compile error, not a dangling read.** `mesh_` below is a
    // pointer captured in the constructor, so `PlaneSweep(makeBox(...), n).below(z)`
    // would bind it to something destroyed at the end of the full expression and
    // every later call would read freed heap. All eight construction sites bind
    // something longer-lived today; deleting this overload is what keeps that true,
    // where a comment would only record it.
    //
    // It cannot be covered by a test instead: a test that constructed one would *be*
    // the undefined behaviour rather than a check on it. Prevention and detection are
    // not interchangeable here, and only one of them is available.
    PlaneSweep(TriMesh&& mesh, const Vec3& n) = delete;

    // Volume and centroid of the region satisfying dot(n, x) <= offset.
    VolumeIntegral below(double offset) const;

    // Offset at which below() returns targetVolume. `guess` warm-starts a tight
    // bracket; pass the mesh's total volume to avoid recomputing it.
    double solveOffsetForVolume(double targetVolume, double fullVolume, double guess,
                                double bracket = 0.35, int iterations = 40) const;

    double loOffset() const { return lo_; }
    double hiOffset() const { return hi_; }

private:
    double illinois(double targetVolume, double lo, double hi, double flo, double fhi,
                    double tol, int iterations) const;

    const TriMesh*      mesh_;
    Vec3                n_;
    std::vector<double> projected_;  // dot(n, verts[i])
    double              lo_ = 0, hi_ = 0;
};

// Second moment of area of the plane-mesh intersection (the "waterplane"), about
// axes through the intersection's own centroid. i.x is the transverse moment that
// drives free-surface effect; i.y is the longitudinal moment. Returns area in .z.
// Approximated by finite difference of the clipped volume centroid, which is exact
// for prismatic sections and adequate elsewhere.
Vec3 waterplaneMoments(const TriMesh& mesh, const Vec3& n, double planeOffset);

// Volume and centroid of the closed mesh below a general free surface z = h(x, y).
//
// integrateBelowPlane() gets to pick a reference point *on* the cutting plane, so
// every cap tetrahedron is degenerate and the open boundary needs no geometry. A
// wavy surface has no such point, and building the cap explicitly would mean
// stitching an intersection curve -- fragile, and wrong wherever the curve leaves
// the mesh's own tessellation.
//
// So this uses a different device: integrate a vector field that vanishes *on the
// surface*. With F = (0, 0, z - h(x, y)) the divergence is 1, so by the divergence
// theorem the enclosed volume equals the flux of F through the boundary -- and the
// free-surface part of that boundary contributes exactly nothing, because F is
// zero there. Only the wetted hull is integrated. The moments come from the same
// trick with fields whose divergences are x, y and z and which also vanish at
// z = h.
//
// The one approximation left is where a triangle crosses the surface: the crossing
// is found by linear interpolation of (z - h) along each edge, so the waterline is
// piecewise-linear over the mesh's own tessellation. That is exact for a flat
// surface, and **measured fourth-order** against a sinusoid -- the error falls by
// a factor of sixteen per halving of triangle size, which is the edge-midpoint
// quadrature's leading term. A box under a cosine surface reaches 1e-11 relative
// error at 1028 triangles.
//
// **The mesh must resolve the surface.** This is the one way to use the routine
// badly, and it fails in the dangerous direction: a panel spanning several
// wavelengths samples the surface at three points and reports whatever those
// three points say, so the error is a systematic, phase-dependent gain or loss
// of volume rather than noise. Measured on a 60 m barge under a 12 m wave, a
// single-panel bottom face invents plus or minus 6% of the ship's displacement
// depending only on where the crests happen to fall; two panels along the length
// already give the exact answer. Several panels per wavelength, always.
//
// `height` is any callable with signature double(double x, double y).
template <typename HeightField>
VolumeIntegral integrateBelowSurface(const TriMesh& mesh, HeightField&& height) {
    double volume = 0;
    double magnitude = 0;
    Vec3 moment{};

    // Evaluate the surface once per *vertex*, not once per triangle corner. On a
    // closed mesh Euler's formula gives tris ~ 2*verts, so the corner form asks
    // the same question about six times over -- measured at 3588 queries against
    // 600 distinct points on the ferry hull. That matters because the surface
    // query is not incidental: with a 128-component wave spectrum it accounts for
    // essentially the entire tick, so removing the redundancy is worth more than
    // vectorising what remains.
    std::vector<double> vertexHeight(mesh.verts.size());
    for (std::size_t i = 0; i < mesh.verts.size(); ++i)
        vertexHeight[i] = height(mesh.verts[i].x, mesh.verts[i].y);

    Vec3 poly[8];
    for (const Tri& t : mesh.tris) {
        const Vec3 tri[3] = {mesh.verts[t.a], mesh.verts[t.b], mesh.verts[t.c]};
        const double surfaceAt[3] = {vertexHeight[t.a], vertexHeight[t.b], vertexHeight[t.c]};
        double depth[3];
        for (int i = 0; i < 3; ++i) depth[i] = tri[i].z - surfaceAt[i];
        if (depth[0] > 0 && depth[1] > 0 && depth[2] > 0) continue;  // wholly dry

        int count = 0;
        for (int i = 0; i < 3; ++i) {
            const int j = (i + 1) % 3;
            if (depth[i] <= 0) poly[count++] = tri[i];
            if ((depth[i] < 0 && depth[j] > 0) || (depth[i] > 0 && depth[j] < 0))
                poly[count++] = tri[i] + (tri[j] - tri[i]) * (depth[i] / (depth[i] - depth[j]));
        }
        if (count < 3) continue;

        for (int i = 1; i + 1 < count; ++i) {
            const Vec3& a = poly[0];
            const Vec3& b = poly[i];
            const Vec3& c = poly[i + 1];
            // Outward normal times area. Only the z component matters: F points
            // along z, so the sides of the hull contribute through their vertical
            // projection alone.
            const double areaZ = 0.5 * cross(b - a, c - a).z;
            if (areaZ == 0.0) continue;

            // Three-point edge-midpoint quadrature, exact for quadratic
            // integrands over a triangle -- which covers the z-moment field
            // exactly and the rest to the same order as the surface sampling.
            const Vec3 quadrature[3] = {(a + b) * 0.5, (b + c) * 0.5, (c + a) * 0.5};
            double meanDepth = 0;
            Vec3 meanMoment{};
            for (const Vec3& q : quadrature) {
                const double surface = height(q.x, q.y);
                const double below = q.z - surface;
                meanDepth += below;
                meanMoment.x += q.x * below;
                meanMoment.y += q.y * below;
                meanMoment.z += 0.5 * (q.z * q.z - surface * surface);
            }
            const double contribution = areaZ * meanDepth / 3.0;
            volume += contribution;
            magnitude += std::abs(contribution);
            moment += meanMoment * (areaZ / 3.0);
        }
    }

    VolumeIntegral result;
    result.volume = volume;
    // Same guard, same derivation as integrateBelowPlane(): the centroid is
    // meaningful only above the round-off floor of its own accumulation, which
    // scales as L^3 and so cannot be a constant in m^3. `magnitude` is the same
    // per-panel contributions summed without their signs -- also m^3, so the
    // ratio is nondimensional -- and it bounds that floor, each `volume +=` being
    // in error by at most eps * |running sum|. Panels above and below the surface
    // enter with opposite signs here, so the cancellation this protects against
    // is if anything larger than in the plane case.
    if (volume > 1e-12 * magnitude) result.centroid = moment / volume;
    return result;
}

// Solve for the plane offset that makes integrateBelowPlane() return targetVolume --
// i.e. "where does the water surface sit inside this compartment?". Bracketed
// regula falsi; volume is monotone in the offset, so it always converges.
double solvePlaneOffsetForVolume(const TriMesh& mesh, const Vec3& n, double targetVolume,
                                 double loOffset, double hiOffset, int iterations = 40);

// Warm-started form, for the per-tick case where the answer moves only slightly
// between calls. `guess` is the previous solution and `fullVolume` the mesh's
// cached total, which saves a whole-mesh pass on every call. Brackets tightly
// around the guess and widens only if that bracket does not contain the root, so
// the common case costs a handful of mesh clips instead of a dozen or more.
double solvePlaneOffsetForVolumeWarm(const TriMesh& mesh, const Vec3& n, double targetVolume,
                                     double fullVolume, double guess, double bracket = 0.35,
                                     int iterations = 40);

// --- Constructive solid geometry --------------------------------------------
//
// integrateBelowPlane() answers "how much is below this plane" without ever
// building the cut surface, which is what makes it fast enough to run per tick.
// Authoring a compartment is the opposite problem: it happens once, at load, and
// it needs the actual capped solid so the result can be clipped again later.

// True when every directed edge occurs exactly once -- that is, the mesh is a
// closed, consistently wound 2-manifold. Every volume integral in this engine
// silently returns nonsense on a mesh that fails this, so anything built by a
// generator, an importer or a boolean should be checked once at load.
bool isClosedManifold(const TriMesh& mesh, double weldEpsilon = 1e-6);

// A retained half-space: the region satisfying dot(n, x) <= offset.
struct HalfSpace {
    Vec3   n;
    double offset = 0;
};

// Clip a closed mesh by a half-space, returning a closed mesh.
//
// The cut is capped: boundary edges left by the clip are welded, chained into
// loops and triangulated by ear clipping. Assumes the intersection loops are
// simple and non-nested, which holds for ship sections cut by bulkhead and deck
// planes but would not for, say, a plane slicing through a hollow mast.
TriMesh clipByPlane(const TriMesh& mesh, const Vec3& n, double offset,
                    double weldEpsilon = 1e-6);

TriMesh clipByHalfSpaces(const TriMesh& mesh, const std::vector<HalfSpace>& halfSpaces);

// The common authoring case: carve a compartment out of the hull interior with
// six axis-aligned bulkhead and deck planes. Unlike makeBox(), the result follows
// the hull form, so a wing tank tapers into the turn of the bilge and a forepeak
// narrows into the stem instead of poking out through the shell.
TriMesh clipToBox(const TriMesh& mesh, const Vec3& lo, const Vec3& hi);

// --- Primitive builders -----------------------------------------------------

// Axis-aligned box from min to max corner, outward winding.
TriMesh makeBox(const Vec3& lo, const Vec3& hi);

// A ship hull generated from transverse station offsets, the way real hull forms
// are tabulated. Each station is a half-breadth curve sampled at waterlines; the
// mesh is mirrored about the centreplane and capped at bow, stern and keel.
struct Station {
    double x = 0;                  // longitudinal position, m
    std::vector<double> halfBeam;  // half-breadth at each waterline, m
};
//
// **A station must carry one half-breadth per waterline.** A short one used to be
// padded silently with zeros, which turns the missing levels into a knife edge --
// and the mesh still closes, still passes the manifold check, and still
// integrates, to a displacement that is simply wrong. The padding now carries the
// last supplied value upward instead, which is at worst a wall-sided extension,
// and `problems` reports any station that needed it. Passing nullptr accepts the
// repair without being told, which is why the loaders check the count themselves.
TriMesh makeHullFromStations(const std::vector<Station>& stations,
                             const std::vector<double>& waterlines,
                             std::vector<std::string>* problems = nullptr);

// Transform every vertex of a mesh by R*v + t. Used to bring compartment meshes
// from their authoring frame into the ship body frame.
TriMesh transformed(const TriMesh& mesh, const Mat3& R, const Vec3& t);

}  // namespace sim
