// SPDX-License-Identifier: MIT
#include "geometry.hpp"

#include <limits>

namespace sim {
namespace {

// Signed volume and first moment of the tetrahedron (o, a, b, c).
inline void accumulateTet(const Vec3& o, const Vec3& a, const Vec3& b, const Vec3& c,
                          double& volume, Vec3& moment) {
    const Vec3 pa = a - o, pb = b - o, pc = c - o;
    const double v = dot(pa, cross(pb, pc)) / 6.0;
    volume += v;
    moment += (o + a + b + c) * (0.25 * v);
}

void flipWinding(TriMesh& mesh) {
    for (Tri& t : mesh.tris) std::swap(t.b, t.c);
}

// Ensure outward winding so integrate() yields a positive volume.
void orientOutward(TriMesh& mesh) {
    if (integrate(mesh).volume < 0) flipWinding(mesh);
}

void addQuad(TriMesh& m, std::uint32_t a, std::uint32_t b, std::uint32_t c, std::uint32_t d) {
    m.tris.push_back({a, b, c});
    m.tris.push_back({a, c, d});
}

// Extent of the mesh projected onto n, used to bracket the bisection.
void projectedExtent(const TriMesh& mesh, const Vec3& n, double& lo, double& hi) {
    lo = std::numeric_limits<double>::infinity();
    hi = -std::numeric_limits<double>::infinity();
    for (const Vec3& v : mesh.verts) {
        const double d = dot(n, v);
        lo = std::min(lo, d);
        hi = std::max(hi, d);
    }
}

}  // namespace

void TriMesh::append(const TriMesh& other) {
    const auto base = static_cast<std::uint32_t>(verts.size());
    verts.insert(verts.end(), other.verts.begin(), other.verts.end());
    for (const Tri& t : other.tris)
        tris.push_back({t.a + base, t.b + base, t.c + base});
}

VolumeIntegral integrate(const TriMesh& mesh) {
    double volume = 0;
    Vec3 moment{};
    const Vec3 origin{0, 0, 0};
    for (const Tri& t : mesh.tris)
        accumulateTet(origin, mesh.verts[t.a], mesh.verts[t.b], mesh.verts[t.c], volume, moment);

    VolumeIntegral r;
    r.volume = volume;
    if (std::abs(volume) > 1e-12) r.centroid = moment / volume;
    return r;
}

VolumeIntegral integrateBelowPlane(const TriMesh& mesh, const Vec3& n, double planeOffset) {
    // Reference point on the cutting plane: cap tetrahedra collapse to zero volume,
    // so the open boundary needs no explicit geometry.
    const Vec3 o = n * planeOffset;

    double volume = 0;
    Vec3 moment{};

    Vec3 poly[8];
    for (const Tri& t : mesh.tris) {
        const Vec3 tri[3] = {mesh.verts[t.a], mesh.verts[t.b], mesh.verts[t.c]};

        // Sutherland-Hodgman against the half-space dot(n, x) <= planeOffset.
        int count = 0;
        for (int i = 0; i < 3; ++i) {
            const Vec3& p = tri[i];
            const Vec3& q = tri[(i + 1) % 3];
            const double dp = dot(n, p) - planeOffset;
            const double dq = dot(n, q) - planeOffset;
            if (dp <= 0) poly[count++] = p;
            if ((dp < 0 && dq > 0) || (dp > 0 && dq < 0))
                poly[count++] = p + (q - p) * (dp / (dp - dq));
        }
        for (int i = 1; i + 1 < count; ++i)
            accumulateTet(o, poly[0], poly[i], poly[i + 1], volume, moment);
    }

    VolumeIntegral r;
    r.volume = volume;
    if (volume > 1e-12) r.centroid = moment / volume;
    return r;
}

double solvePlaneOffsetForVolume(const TriMesh& mesh, const Vec3& n, double targetVolume,
                                 double loOffset, double hiOffset, int iterations) {
    double meshLo, meshHi;
    projectedExtent(mesh, n, meshLo, meshHi);
    // A hair of slack so the fully-empty and fully-full ends stay bracketed.
    const double pad = 1e-6 * std::max(1.0, meshHi - meshLo);
    double lo = std::max(loOffset, meshLo - pad);
    double hi = std::min(hiOffset, meshHi + pad);
    if (!(lo < hi)) { lo = meshLo - pad; hi = meshHi + pad; }

    if (targetVolume <= 0) return lo;
    const double full = integrate(mesh).volume;
    if (targetVolume >= full) return hi;

    // Illinois-modified regula falsi. Volume is monotone and smooth in the offset,
    // so this converges superlinearly while keeping the root bracketed -- roughly a
    // third of the mesh clips plain bisection needs for the same precision, which
    // matters because this runs per compartment per tick.
    double flo = -targetVolume;         // volume(lo) is 0 by construction
    double fhi = full - targetVolume;   // volume(hi) is the full mesh
    const double tol = 1e-9 * std::max(1.0, full);

    double x = 0.5 * (lo + hi);
    for (int i = 0; i < iterations; ++i) {
        x = lo - flo * (hi - lo) / (fhi - flo);
        const double fx = integrateBelowPlane(mesh, n, x).volume - targetVolume;
        if (std::abs(fx) < tol || hi - lo < 1e-12) break;
        if (fx < 0) { lo = x; flo = fx; fhi *= 0.5; }
        else        { hi = x; fhi = fx; flo *= 0.5; }
    }
    return x;
}

Vec3 waterplaneMoments(const TriMesh& mesh, const Vec3& n, double planeOffset) {
    const VolumeIntegral base = integrateBelowPlane(mesh, n, planeOffset);
    if (base.volume <= 1e-9) return {0, 0, 0};

    // Area from the derivative of displaced volume with respect to plane offset.
    const double h = 1e-3;
    const double area =
        (integrateBelowPlane(mesh, n, planeOffset + h).volume -
         integrateBelowPlane(mesh, n, planeOffset - h).volume) / (2 * h);

    // Metacentric radii from the centroid shift under a small constant-volume tilt.
    // BM = I / V is the definition; measuring BM directly avoids ever having to
    // extract and integrate the intersection polygon.
    const double delta = 1e-4;
    auto bmAbout = [&](const Vec3& axis, const Vec3& shiftDir) {
        const Mat3 R = Quat::fromAxisAngle(axis, delta).toMat3();
        const Vec3 n2 = R * n;
        const double off2 = solvePlaneOffsetForVolume(
            mesh, n2, base.volume, -std::numeric_limits<double>::infinity(),
            std::numeric_limits<double>::infinity());
        const VolumeIntegral tilted = integrateBelowPlane(mesh, n2, off2);
        return std::abs(dot(tilted.centroid - base.centroid, shiftDir)) / delta;
    };

    const Vec3 ex{1, 0, 0}, ey{0, 1, 0};
    const double bmTransverse   = bmAbout(ex, ey);
    const double bmLongitudinal = bmAbout(ey, ex);
    return {bmTransverse * base.volume, bmLongitudinal * base.volume, area};
}

TriMesh makeBox(const Vec3& lo, const Vec3& hi) {
    TriMesh m;
    m.verts = {
        {lo.x, lo.y, lo.z}, {hi.x, lo.y, lo.z}, {hi.x, hi.y, lo.z}, {lo.x, hi.y, lo.z},
        {lo.x, lo.y, hi.z}, {hi.x, lo.y, hi.z}, {hi.x, hi.y, hi.z}, {lo.x, hi.y, hi.z},
    };
    addQuad(m, 0, 3, 2, 1);  // bottom (-z)
    addQuad(m, 4, 5, 6, 7);  // top (+z)
    addQuad(m, 0, 1, 5, 4);  // -y
    addQuad(m, 2, 3, 7, 6);  // +y
    addQuad(m, 1, 2, 6, 5);  // +x
    addQuad(m, 3, 0, 4, 7);  // -x
    orientOutward(m);
    return m;
}

TriMesh makeHullFromStations(const std::vector<Station>& stations,
                             const std::vector<double>& waterlines) {
    TriMesh m;
    const auto ns = stations.size();
    const auto nw = waterlines.size();
    if (ns < 2 || nw < 2) return m;

    // Two vertex sheets: port (+y) and starboard (-y).
    const auto port = [&](std::size_t i, std::size_t k) {
        return static_cast<std::uint32_t>(i * nw + k);
    };
    const auto stbd = [&](std::size_t i, std::size_t k) {
        return static_cast<std::uint32_t>(ns * nw + i * nw + k);
    };

    m.verts.resize(2 * ns * nw);
    for (std::size_t i = 0; i < ns; ++i)
        for (std::size_t k = 0; k < nw; ++k) {
            const double hb = k < stations[i].halfBeam.size() ? stations[i].halfBeam[k] : 0.0;
            m.verts[port(i, k)] = {stations[i].x, hb, waterlines[k]};
            m.verts[stbd(i, k)] = {stations[i].x, -hb, waterlines[k]};
        }

    // Shell plating, both sides.
    for (std::size_t i = 0; i + 1 < ns; ++i)
        for (std::size_t k = 0; k + 1 < nw; ++k) {
            addQuad(m, port(i, k), port(i + 1, k), port(i + 1, k + 1), port(i, k + 1));
            addQuad(m, stbd(i, k), stbd(i, k + 1), stbd(i + 1, k + 1), stbd(i + 1, k));
        }

    // Keel strip and weather deck close the sheets to each other.
    for (std::size_t i = 0; i + 1 < ns; ++i) {
        addQuad(m, port(i, 0), port(i + 1, 0), stbd(i + 1, 0), stbd(i, 0));
        addQuad(m, port(i, nw - 1), stbd(i, nw - 1), stbd(i + 1, nw - 1), port(i + 1, nw - 1));
    }

    // Transom and stem caps.
    for (std::size_t k = 0; k + 1 < nw; ++k) {
        addQuad(m, port(0, k), port(0, k + 1), stbd(0, k + 1), stbd(0, k));
        addQuad(m, port(ns - 1, k), stbd(ns - 1, k), stbd(ns - 1, k + 1), port(ns - 1, k + 1));
    }

    orientOutward(m);
    return m;
}

TriMesh transformed(const TriMesh& mesh, const Mat3& R, const Vec3& t) {
    TriMesh out = mesh;
    for (Vec3& v : out.verts) v = R * v + t;
    return out;
}

}  // namespace sim
