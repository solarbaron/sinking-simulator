// SPDX-License-Identifier: MIT
#include "geometry.hpp"

#include <string>

#include <limits>
#include <unordered_map>
#include <utility>

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

// Shared tail of both offset solvers: Illinois-modified regula falsi on an
// already-established bracket, where flo and fhi are the residuals at lo and hi.
double illinois(const TriMesh& mesh, const Vec3& n, double targetVolume, double lo, double hi,
                double flo, double fhi, double tol, int iterations) {
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

    // Volume is monotone and smooth in the offset, so Illinois-modified regula
    // falsi converges superlinearly while keeping the root bracketed.
    return illinois(mesh, n, targetVolume, lo, hi, -targetVolume, full - targetVolume,
                    1e-9 * std::max(1.0, full), iterations);
}

// --- PlaneSweep -------------------------------------------------------------

PlaneSweep::PlaneSweep(const TriMesh& mesh, const Vec3& n) : mesh_(&mesh), n_(normalize(n)) {
    projected_.resize(mesh.verts.size());
    lo_ = std::numeric_limits<double>::infinity();
    hi_ = -std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < mesh.verts.size(); ++i) {
        const double d = dot(n_, mesh.verts[i]);
        projected_[i] = d;
        lo_ = std::min(lo_, d);
        hi_ = std::max(hi_, d);
    }
    const double pad = 1e-6 * std::max(1.0, hi_ - lo_);
    lo_ -= pad;
    hi_ += pad;
}

VolumeIntegral PlaneSweep::below(double offset) const {
    const Vec3 o = n_ * offset;
    double volume = 0;
    Vec3 moment{};

    Vec3 poly[8];
    for (const Tri& t : mesh_->tris) {
        const std::uint32_t vi[3] = {t.a, t.b, t.c};
        // The only arithmetic per vertex is a subtraction; the projections were
        // done once in the constructor.
        const double d[3] = {projected_[t.a] - offset, projected_[t.b] - offset,
                             projected_[t.c] - offset};
        if (d[0] > 0 && d[1] > 0 && d[2] > 0) continue;  // wholly above the plane

        int count = 0;
        for (int i = 0; i < 3; ++i) {
            const int j = (i + 1) % 3;
            const Vec3& p = mesh_->verts[vi[i]];
            if (d[i] <= 0) poly[count++] = p;
            if ((d[i] < 0 && d[j] > 0) || (d[i] > 0 && d[j] < 0))
                poly[count++] = p + (mesh_->verts[vi[j]] - p) * (d[i] / (d[i] - d[j]));
        }
        for (int i = 1; i + 1 < count; ++i)
            accumulateTet(o, poly[0], poly[i], poly[i + 1], volume, moment);
    }

    VolumeIntegral r;
    r.volume = volume;
    if (volume > 1e-12) r.centroid = moment / volume;
    return r;
}

double PlaneSweep::illinois(double targetVolume, double lo, double hi, double flo, double fhi,
                            double tol, int iterations) const {
    double x = 0.5 * (lo + hi);
    for (int i = 0; i < iterations; ++i) {
        x = lo - flo * (hi - lo) / (fhi - flo);
        const double fx = below(x).volume - targetVolume;
        if (std::abs(fx) < tol || hi - lo < 1e-12) break;
        if (fx < 0) { lo = x; flo = fx; fhi *= 0.5; }
        else        { hi = x; fhi = fx; flo *= 0.5; }
    }
    return x;
}

double PlaneSweep::solveOffsetForVolume(double targetVolume, double fullVolume, double guess,
                                        double bracket, int iterations) const {
    if (targetVolume <= 0) return lo_;
    if (targetVolume >= fullVolume) return hi_;
    const double tol = 1e-9 * std::max(1.0, fullVolume);

    // Two probes either side of last tick's answer. Whenever the level moved less
    // than `bracket`, the root is boxed in for the cost of two clips.
    const double lo = std::max(lo_, guess - bracket);
    const double hi = std::min(hi_, guess + bracket);
    if (lo < hi) {
        const double flo = below(lo).volume - targetVolume;
        if (flo <= 0) {
            const double fhi = below(hi).volume - targetVolume;
            if (fhi >= 0) return illinois(targetVolume, lo, hi, flo, fhi, tol, iterations);
        }
    }
    return illinois(targetVolume, lo_, hi_, -targetVolume, fullVolume - targetVolume, tol,
                    iterations);
}

double solvePlaneOffsetForVolumeWarm(const TriMesh& mesh, const Vec3& n, double targetVolume,
                                     double fullVolume, double guess, double bracket,
                                     int iterations) {
    double meshLo, meshHi;
    projectedExtent(mesh, n, meshLo, meshHi);
    const double pad = 1e-6 * std::max(1.0, meshHi - meshLo);
    const double extentLo = meshLo - pad, extentHi = meshHi + pad;

    if (targetVolume <= 0) return extentLo;
    if (targetVolume >= fullVolume) return extentHi;

    const double tol = 1e-9 * std::max(1.0, fullVolume);

    // Two probes either side of last tick's answer. If they straddle the target
    // -- which they do whenever the water level moved less than `bracket` -- the
    // root is boxed in for the cost of two clips.
    const double lo = std::max(extentLo, guess - bracket);
    const double hi = std::min(extentHi, guess + bracket);
    if (lo < hi) {
        const double flo = integrateBelowPlane(mesh, n, lo).volume - targetVolume;
        if (flo <= 0) {
            const double fhi = integrateBelowPlane(mesh, n, hi).volume - targetVolume;
            if (fhi >= 0) return illinois(mesh, n, targetVolume, lo, hi, flo, fhi, tol, iterations);
        }
    }
    return illinois(mesh, n, targetVolume, extentLo, extentHi, -targetVolume,
                    fullVolume - targetVolume, tol, iterations);
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

// --- Constructive solid geometry --------------------------------------------

namespace {

// Merges coincident positions so that intersection points produced independently
// by neighbouring triangles become one vertex, which is what lets the cut edges
// chain into loops at all.
class VertexWelder {
public:
    explicit VertexWelder(double epsilon) : eps_(epsilon), inv_(1.0 / epsilon) {}

    std::uint32_t add(const Vec3& p) {
        const long long kx = std::llround(p.x * inv_);
        const long long ky = std::llround(p.y * inv_);
        const long long kz = std::llround(p.z * inv_);
        // A point can land either side of a cell boundary, so search the 27-cell
        // neighbourhood rather than just its own bucket.
        for (int dx = -1; dx <= 1; ++dx)
            for (int dy = -1; dy <= 1; ++dy)
                for (int dz = -1; dz <= 1; ++dz) {
                    auto it = grid_.find(key(kx + dx, ky + dy, kz + dz));
                    if (it == grid_.end()) continue;
                    for (std::uint32_t i : it->second)
                        if (length2(verts[i] - p) <= eps_ * eps_) return i;
                }
        const auto index = static_cast<std::uint32_t>(verts.size());
        verts.push_back(p);
        grid_[key(kx, ky, kz)].push_back(index);
        return index;
    }

    std::vector<Vec3> verts;

private:
    static std::uint64_t key(long long x, long long y, long long z) {
        // Cheap spatial hash; collisions are handled by the bucket scan above.
        return static_cast<std::uint64_t>(x) * 0x9E3779B97F4A7C15ull
             ^ static_cast<std::uint64_t>(y) * 0xC2B2AE3D27D4EB4Full
             ^ static_cast<std::uint64_t>(z) * 0x165667B19E3779F9ull;
    }

    double eps_, inv_;
    std::unordered_map<std::uint64_t, std::vector<std::uint32_t>> grid_;
};

struct Vec2 {
    double x = 0, y = 0;
};

double signedArea(const std::vector<Vec2>& p) {
    double a = 0;
    for (std::size_t i = 0, n = p.size(); i < n; ++i) {
        const Vec2& u = p[i];
        const Vec2& v = p[(i + 1) % n];
        a += u.x * v.y - v.x * u.y;
    }
    return 0.5 * a;
}

double cross2(const Vec2& a, const Vec2& b, const Vec2& c) {
    return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

bool pointInTriangle(const Vec2& p, const Vec2& a, const Vec2& b, const Vec2& c) {
    const double d1 = cross2(a, b, p), d2 = cross2(b, c, p), d3 = cross2(c, a, p);
    return (d1 >= 0 && d2 >= 0 && d3 >= 0) || (d1 <= 0 && d2 <= 0 && d3 <= 0);
}

// Ear clipping. O(n^2), which is irrelevant: these loops are ship sections with a
// few dozen vertices and this runs once at load.
void earClip(std::vector<Vec2> poly, std::vector<std::uint32_t> idx, TriMesh& out) {
    if (poly.size() < 3) return;
    if (signedArea(poly) < 0) {  // force counter-clockwise so the cap faces +n
        std::reverse(poly.begin(), poly.end());
        std::reverse(idx.begin(), idx.end());
    }

    while (poly.size() > 3) {
        bool clipped = false;
        const std::size_t n = poly.size();
        for (std::size_t i = 0; i < n; ++i) {
            const std::size_t prev = (i + n - 1) % n;
            const std::size_t next = (i + 1) % n;
            if (cross2(poly[prev], poly[i], poly[next]) <= 0) continue;  // reflex

            bool contains = false;
            for (std::size_t j = 0; j < n && !contains; ++j) {
                if (j == prev || j == i || j == next) continue;
                contains = pointInTriangle(poly[j], poly[prev], poly[i], poly[next]);
            }
            if (contains) continue;

            out.tris.push_back({idx[prev], idx[i], idx[next]});
            poly.erase(poly.begin() + static_cast<long>(i));
            idx.erase(idx.begin() + static_cast<long>(i));
            clipped = true;
            break;
        }
        // Degenerate or self-touching loop: fall back to a fan rather than
        // spinning forever or dropping the cap and leaking the solid.
        if (!clipped) break;
    }
    for (std::size_t i = 1; i + 1 < idx.size(); ++i)
        out.tris.push_back({idx[0], idx[i], idx[i + 1]});
}

}  // namespace

bool isClosedManifold(const TriMesh& mesh, double weldEpsilon) {
    if (mesh.tris.empty()) return false;

    // Weld first: two triangles can reference distinct vertex records holding the
    // same position, and topologically they still share that edge.
    VertexWelder welder(weldEpsilon);
    std::vector<std::uint32_t> remap(mesh.verts.size());
    for (std::size_t i = 0; i < mesh.verts.size(); ++i)
        remap[i] = welder.add(mesh.verts[i]);

    std::unordered_map<std::uint64_t, int> directed;
    auto key = [](std::uint32_t a, std::uint32_t b) {
        return (static_cast<std::uint64_t>(a) << 32) | b;
    };

    for (const Tri& t : mesh.tris) {
        const std::uint32_t v[3] = {remap[t.a], remap[t.b], remap[t.c]};
        if (v[0] == v[1] || v[1] == v[2] || v[2] == v[0]) continue;  // degenerate
        for (int i = 0; i < 3; ++i) ++directed[key(v[i], v[(i + 1) % 3])];
    }

    for (const auto& [k, count] : directed) {
        if (count != 1) return false;  // an edge traversed twice the same way
        const auto a = static_cast<std::uint32_t>(k >> 32);
        const auto b = static_cast<std::uint32_t>(k & 0xFFFFFFFFull);
        auto opposite = directed.find(key(b, a));
        if (opposite == directed.end() || opposite->second != 1) return false;
    }
    return true;
}

TriMesh clipByPlane(const TriMesh& mesh, const Vec3& nIn, double offset, double weldEpsilon) {
    const Vec3 n = normalize(nIn);
    double meshLo, meshHi;
    projectedExtent(mesh, n, meshLo, meshHi);
    if (meshLo >= offset) return {};        // nothing survives
    if (meshHi <= offset) return mesh;      // the plane misses the mesh entirely

    const double planeTol = 1e-9 * std::max(1.0, meshHi - meshLo);

    VertexWelder welder(weldEpsilon);
    TriMesh out;
    std::vector<std::pair<std::uint32_t, std::uint32_t>> cutEdges;

    Vec3 poly[8];
    bool onPlane[8];
    for (const Tri& t : mesh.tris) {
        const Vec3 tri[3] = {mesh.verts[t.a], mesh.verts[t.b], mesh.verts[t.c]};

        int count = 0;
        for (int i = 0; i < 3; ++i) {
            const Vec3& p = tri[i];
            const Vec3& q = tri[(i + 1) % 3];
            const double dp = dot(n, p) - offset;
            const double dq = dot(n, q) - offset;
            if (dp <= 0) {
                poly[count] = p;
                onPlane[count] = std::abs(dp) <= planeTol;
                ++count;
            }
            if ((dp < 0 && dq > 0) || (dp > 0 && dq < 0)) {
                poly[count] = p + (q - p) * (dp / (dp - dq));
                onPlane[count] = true;
                ++count;
            }
        }
        if (count < 3) continue;

        std::uint32_t idx[8];
        for (int i = 0; i < count; ++i) idx[i] = welder.add(poly[i]);

        // A boundary edge of the clipped surface is one whose endpoints both lie
        // in the cutting plane. Collect them now; they become the cap.
        for (int i = 0; i < count; ++i) {
            const int j = (i + 1) % count;
            if (onPlane[i] && onPlane[j] && idx[i] != idx[j])
                cutEdges.emplace_back(idx[i], idx[j]);
        }

        for (int i = 1; i + 1 < count; ++i)
            out.tris.push_back({idx[0], idx[i], idx[i + 1]});
    }
    out.verts = welder.verts;

    // A face lying exactly in the cutting plane contributes an edge in each
    // direction. Those pairs are interior, not boundary, and must cancel.
    std::vector<bool> dropped(cutEdges.size(), false);
    for (std::size_t i = 0; i < cutEdges.size(); ++i) {
        if (dropped[i]) continue;
        for (std::size_t j = i + 1; j < cutEdges.size(); ++j) {
            if (dropped[j]) continue;
            if (cutEdges[j].first == cutEdges[i].second &&
                cutEdges[j].second == cutEdges[i].first) {
                dropped[i] = dropped[j] = true;
                break;
            }
        }
    }

    // The cap must traverse the boundary opposite to the clipped surface for the
    // combined solid to stay consistently wound, so chain the reversed edges.
    std::unordered_map<std::uint32_t, std::vector<std::uint32_t>> successors;
    std::size_t edgeCount = 0;
    for (std::size_t i = 0; i < cutEdges.size(); ++i) {
        if (dropped[i]) continue;
        successors[cutEdges[i].second].push_back(cutEdges[i].first);
        ++edgeCount;
    }

    // Basis in the cutting plane, chosen so that u x v == n and a counter-clockwise
    // loop in (u, v) has its normal along +n -- outward for the retained side.
    const Vec3 seed = std::abs(n.x) < 0.9 ? Vec3{1, 0, 0} : Vec3{0, 1, 0};
    const Vec3 u = normalize(cross(n, seed));
    const Vec3 v = cross(n, u);

    std::size_t consumed = 0;
    while (consumed < edgeCount) {
        // Find any vertex that still has an outgoing edge.
        std::uint32_t start = 0;
        bool found = false;
        for (auto& [from, tos] : successors)
            if (!tos.empty()) { start = from; found = true; break; }
        if (!found) break;

        std::vector<std::uint32_t> loop;
        std::uint32_t at = start;
        while (true) {
            auto it = successors.find(at);
            if (it == successors.end() || it->second.empty()) break;
            const std::uint32_t next = it->second.back();
            it->second.pop_back();
            ++consumed;
            loop.push_back(at);
            at = next;
            if (at == start) break;
            if (loop.size() > edgeCount + 1) break;  // malformed; bail out
        }
        if (loop.size() < 3) continue;

        std::vector<Vec2> flat;
        flat.reserve(loop.size());
        for (std::uint32_t i : loop)
            flat.push_back({dot(out.verts[i], u), dot(out.verts[i], v)});
        earClip(std::move(flat), loop, out);
    }

    return out;
}

TriMesh clipByHalfSpaces(const TriMesh& mesh, const std::vector<HalfSpace>& halfSpaces) {
    TriMesh result = mesh;
    for (const HalfSpace& h : halfSpaces) {
        if (result.tris.empty()) break;
        result = clipByPlane(result, h.n, h.offset);
    }
    return result;
}

TriMesh clipToBox(const TriMesh& mesh, const Vec3& lo, const Vec3& hi) {
    return clipByHalfSpaces(mesh, {{{ 1, 0, 0},  hi.x}, {{-1, 0, 0}, -lo.x},
                                   {{ 0, 1, 0},  hi.y}, {{ 0,-1, 0}, -lo.y},
                                   {{ 0, 0, 1},  hi.z}, {{ 0, 0,-1}, -lo.z}});
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
                             const std::vector<double>& waterlines,
                             std::vector<std::string>* problems) {
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
    if (problems)
        for (std::size_t i = 0; i < ns; ++i) {
            if (stations[i].halfBeam.size() == nw) continue;
            problems->push_back("station at x = " + std::to_string(stations[i].x) + " has " +
                                std::to_string(stations[i].halfBeam.size()) +
                                " half-breadths for " + std::to_string(nw) + " waterlines");
        }
    for (std::size_t i = 0; i < ns; ++i)
        for (std::size_t k = 0; k < nw; ++k) {
            // Carry the last supplied half-breadth upward rather than dropping to
            // zero. A missing column is an authoring error either way, but a
            // wall-sided extension is a visible one and a knife edge is not.
            double hb = 0.0;
            if (!stations[i].halfBeam.empty())
                hb = stations[i].halfBeam[std::min(k, stations[i].halfBeam.size() - 1)];
            m.verts[port(i, k)] = {stations[i].x, hb, waterlines[k]};
            m.verts[stbd(i, k)] = {stations[i].x, -hb, waterlines[k]};
        }

    // Shell plating, both sides.
    for (std::size_t i = 0; i + 1 < ns; ++i)
        for (std::size_t k = 0; k + 1 < nw; ++k) {
            addQuad(m, port(i, k), port(i + 1, k), port(i + 1, k + 1), port(i, k + 1));
            addQuad(m, stbd(i, k), stbd(i, k + 1), stbd(i + 1, k + 1), stbd(i + 1, k));
        }

    // Keel strip and weather deck close the sheets to each other. Both are wound
    // to match the plating sheets above -- orientOutward() flips the whole mesh at
    // the end, and it can only do that correctly if every face already agrees.
    for (std::size_t i = 0; i + 1 < ns; ++i) {
        addQuad(m, stbd(i, 0), stbd(i + 1, 0), port(i + 1, 0), port(i, 0));
        addQuad(m, port(i + 1, nw - 1), stbd(i + 1, nw - 1), stbd(i, nw - 1), port(i, nw - 1));
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
