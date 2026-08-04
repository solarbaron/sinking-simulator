// SPDX-License-Identifier: MIT
#include "damage.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <map>

namespace gpu {
namespace {

using sim::Vec3;

constexpr double kInfiniteDistance = 1e300;

// A right-handed in-plane frame about `normal`. Only used where the caller has no
// frame of its own -- the zone brings `Patch::right` and `Patch::up` with it. The
// seed is chosen off the smaller component so the cross product is never near
// degenerate: at |n.z| = 0.9 the seed is still 26 degrees off the normal.
void planeFrame(const Vec3& normal, Vec3& right, Vec3& up) {
    const Vec3 n = normalize(normal);
    const Vec3 seed = std::abs(n.z) < 0.9 ? Vec3{0, 0, 1} : Vec3{1, 0, 0};
    right = normalize(cross(seed, n));
    up = cross(n, right);
}

// Barycentric weights of `p` in the triangle (a, b, c), all in 2D. Returns false
// when the point is outside or the triangle is degenerate.
//
// The epsilon admits a point exactly on an edge, and **which side answers cannot
// matter**: the interpolant is linear along a shared edge and both triangles carry
// the same two nodal values there, so they agree. Mutation testing bears that out
// rather more strongly than intended -- loosening it from 1e-9 to 1e-2 is
// invisible to every assertion in the suite, because a mesh node is a vertex of
// every triangle that could claim it and the answer at one is that node's own
// value whichever triangle is asked. It is a guard against *missing* a point, not
// against getting it wrong.
bool barycentric(double px, double py, double ax, double ay, double bx, double by, double cx,
                 double cy, double w[3]) {
    const double area = (bx - ax) * (cy - ay) - (by - ay) * (cx - ax);
    if (std::abs(area) < 1e-18) return false;
    const double inverse = 1.0 / area;
    w[1] = ((px - ax) * (cy - ay) - (py - ay) * (cx - ax)) * inverse;
    w[2] = ((bx - ax) * (py - ay) - (by - ay) * (px - ax)) * inverse;
    w[0] = 1.0 - w[1] - w[2];
    const double tolerance = 1e-9;
    return w[0] >= -tolerance && w[1] >= -tolerance && w[2] >= -tolerance;
}

// Distance from a point to a segment, in 2D.
double distanceToSegment(double px, double py, double ax, double ay, double bx, double by) {
    const double dx = bx - ax, dy = by - ay;
    const double len2 = dx * dx + dy * dy;
    double t = 0.0;
    if (len2 > 0) t = std::clamp(((px - ax) * dx + (py - ay) * dy) / len2, 0.0, 1.0);
    const double qx = ax + dx * t, qy = ay + dy * t;
    return std::sqrt((px - qx) * (px - qx) + (py - qy) * (py - qy));
}

std::array<long long, 3> weldKey(const Vec3& p, double epsilon) {
    const double inverse = 1.0 / epsilon;
    return {std::llround(p.x * inverse), std::llround(p.y * inverse),
            std::llround(p.z * inverse)};
}

double triangleArea(const Vec3& a, const Vec3& b, const Vec3& c) {
    return 0.5 * length(cross(b - a, c - a));
}

}  // namespace

// --- The deformation field ------------------------------------------------------

bool HullDamage::addZone(const sim::zone::Patch& patch, const std::vector<double>& deformed,
                         std::string& error) {
    if (patch.empty()) {
        error = "the patch has no elements";
        return false;
    }
    if (deformed.size() != patch.mesh.position.size()) {
        error = "the deformed node array does not match the patch (" +
                std::to_string(deformed.size()) + " values against " +
                std::to_string(patch.mesh.position.size()) + ")";
        return false;
    }
    if (patch.outerFace.size() != patch.nodeCount()) {
        error = "the patch does not carry one outer-face flag per node";
        return false;
    }

    Zone zone;
    zone.centre = patch.centre;
    zone.right = patch.right;
    zone.up = patch.up;
    zone.axis = patch.axis;

    // A mid-surface point is the pair (2m, 2m+1): one node on each face of the
    // plate. Its motion is the average of the two, which is the mid-surface
    // displacement with the thickness stretch taken out -- and the mid-surface is
    // where `scantlings.hpp` puts the plating, so it is also where the hull mesh
    // is.
    const std::size_t mids = patch.nodeCount() / 2;
    zone.rest.resize(mids);
    zone.displacement.resize(mids);
    zone.u.resize(mids);
    zone.v.resize(mids);
    for (std::size_t m = 0; m < mids; ++m) {
        Vec3 rest{}, moved{};
        for (int k = 0; k < 3; ++k) {
            const std::size_t lo = (2 * m) * 3 + static_cast<std::size_t>(k);
            const std::size_t hi = (2 * m + 1) * 3 + static_cast<std::size_t>(k);
            rest[k] = 0.5 * (patch.mesh.position[lo] + patch.mesh.position[hi]);
            moved[k] = 0.5 * (deformed[lo] + deformed[hi]);
        }
        zone.rest[m] = rest;
        zone.displacement[m] = moved - rest;
        const Vec3 offset = rest - zone.centre;
        zone.u[m] = dot(offset, zone.right);
        zone.v[m] = dot(offset, zone.up);
    }

    // Each element quad, split into two triangles. Splitting rather than inverting
    // the bilinear map is deliberate: linear interpolation over the two halves is
    // exactly C0 (both halves are linear in the same parameter along the shared
    // diagonal, and along an element edge both neighbours interpolate the same two
    // nodal values), and it is interpolating at the nodes, which is the identity
    // the tests hold the picture to.
    const std::size_t elements = patch.elementCount();
    zone.tri.reserve(elements * 6);
    for (std::size_t e = 0; e < elements; ++e) {
        std::uint32_t corner[4];
        for (int k = 0; k < 4; ++k) corner[k] = patch.mesh.index[e * 8 + static_cast<std::size_t>(k)] / 2u;
        zone.tri.push_back(corner[0]);
        zone.tri.push_back(corner[1]);
        zone.tri.push_back(corner[2]);
        zone.tri.push_back(corner[0]);
        zone.tri.push_back(corner[2]);
        zone.tri.push_back(corner[3]);
    }

    // The boundary loop: an edge of the triangulation seen once. On a clamped
    // patch every node on it is pinned, so the largest displacement over the loop
    // is zero and the field reaches zero before its support ends -- which is the
    // whole reason there is no seam. Measured rather than assumed, because a patch
    // solved with a free edge would not have the property and should say so.
    {
        std::map<std::uint64_t, int> seen;
        for (std::size_t t = 0; t + 2 < zone.tri.size(); t += 3)
            for (int e = 0; e < 3; ++e) {
                const std::uint32_t a = zone.tri[t + static_cast<std::size_t>(e)];
                const std::uint32_t b = zone.tri[t + static_cast<std::size_t>((e + 1) % 3)];
                if (a == b) continue;
                const std::uint64_t low = std::min(a, b), high = std::max(a, b);
                ++seen[(low << 32) | high];
            }
        for (const auto& [key, count] : seen) {
            if (count != 1) continue;
            const auto low = static_cast<std::uint32_t>(key >> 32);
            const auto high = static_cast<std::uint32_t>(key & 0xffffffffu);
            zone.boundary = std::max(zone.boundary, length(zone.displacement[low]));
            zone.boundary = std::max(zone.boundary, length(zone.displacement[high]));
        }
    }

    double uLo = kInfiniteDistance, uHi = -kInfiniteDistance;
    double vLo = kInfiniteDistance, vHi = -kInfiniteDistance;
    for (std::size_t m = 0; m < mids; ++m) {
        uLo = std::min(uLo, zone.u[m]);
        uHi = std::max(uHi, zone.u[m]);
        vLo = std::min(vLo, zone.v[m]);
        vHi = std::max(vHi, zone.v[m]);
        zone.largest = std::max(zone.largest, length(zone.displacement[m]));
        zone.radius = std::max(zone.radius, length(zone.rest[m] - zone.centre));
    }

    // Bucket grid. A triangle is registered in **every** bucket its bounding box
    // overlaps, so a lookup is correct for any cell size and only its occupancy
    // depends on the choice -- the largest triangle extent keeps a bucket to a
    // handful of entries. Mutation testing established that: fixing the cell at a
    // metre is invisible to every assertion in the suite, which is right, and the
    // comment that used to be here claimed it was a correctness condition.
    double widest = 0;
    for (std::size_t t = 0; t + 2 < zone.tri.size(); t += 3) {
        const std::uint32_t a = zone.tri[t], b = zone.tri[t + 1], c = zone.tri[t + 2];
        widest = std::max(widest,
                          std::max({zone.u[a], zone.u[b], zone.u[c]}) -
                              std::min({zone.u[a], zone.u[b], zone.u[c]}));
        widest = std::max(widest,
                          std::max({zone.v[a], zone.v[b], zone.v[c]}) -
                              std::min({zone.v[a], zone.v[b], zone.v[c]}));
    }
    zone.cell = std::max(widest, 1e-6);
    zone.uLo = uLo;
    zone.vLo = vLo;
    zone.nu = std::max(1, static_cast<int>((uHi - uLo) / zone.cell) + 1);
    zone.nv = std::max(1, static_cast<int>((vHi - vLo) / zone.cell) + 1);

    const std::size_t cells = static_cast<std::size_t>(zone.nu) * static_cast<std::size_t>(zone.nv);
    std::vector<std::uint32_t> count(cells + 1, 0);
    const auto cellRange = [&](std::size_t t, int& u0, int& u1, int& v0, int& v1) {
        const std::uint32_t a = zone.tri[t], b = zone.tri[t + 1], c = zone.tri[t + 2];
        const double lo0 = std::min({zone.u[a], zone.u[b], zone.u[c]});
        const double hi0 = std::max({zone.u[a], zone.u[b], zone.u[c]});
        const double lo1 = std::min({zone.v[a], zone.v[b], zone.v[c]});
        const double hi1 = std::max({zone.v[a], zone.v[b], zone.v[c]});
        u0 = std::clamp(static_cast<int>((lo0 - zone.uLo) / zone.cell), 0, zone.nu - 1);
        u1 = std::clamp(static_cast<int>((hi0 - zone.uLo) / zone.cell), 0, zone.nu - 1);
        v0 = std::clamp(static_cast<int>((lo1 - zone.vLo) / zone.cell), 0, zone.nv - 1);
        v1 = std::clamp(static_cast<int>((hi1 - zone.vLo) / zone.cell), 0, zone.nv - 1);
    };
    for (std::size_t t = 0; t + 2 < zone.tri.size(); t += 3) {
        int u0, u1, v0, v1;
        cellRange(t, u0, u1, v0, v1);
        for (int iv = v0; iv <= v1; ++iv)
            for (int iu = u0; iu <= u1; ++iu)
                ++count[static_cast<std::size_t>(iv) * static_cast<std::size_t>(zone.nu) +
                        static_cast<std::size_t>(iu) + 1];
    }
    zone.bucketStart.assign(cells + 1, 0);
    for (std::size_t i = 0; i < cells; ++i) zone.bucketStart[i + 1] = zone.bucketStart[i] + count[i + 1];
    zone.bucketItem.assign(zone.bucketStart[cells], 0);
    std::vector<std::uint32_t> cursor(zone.bucketStart.begin(), zone.bucketStart.end() - 1);
    for (std::size_t t = 0; t + 2 < zone.tri.size(); t += 3) {
        int u0, u1, v0, v1;
        cellRange(t, u0, u1, v0, v1);
        for (int iv = v0; iv <= v1; ++iv)
            for (int iu = u0; iu <= u1; ++iu) {
                const std::size_t cell = static_cast<std::size_t>(iv) *
                                             static_cast<std::size_t>(zone.nu) +
                                         static_cast<std::size_t>(iu);
                zone.bucketItem[cursor[cell]++] = static_cast<std::uint32_t>(t / 3);
            }
    }

    zones_.push_back(std::move(zone));
    return true;
}

void HullDamage::addTent(const Vec3& centre, const Vec3& outward, double halfSpan, double depth) {
    if (!(halfSpan > 0)) return;
    Tent tent;
    tent.centre = centre;
    tent.outward = normalize(outward);
    if (length2(tent.outward) < 0.5) return;
    tent.halfSpan = halfSpan;
    tent.depth = depth;
    tents_.push_back(tent);
}

void HullDamage::addTornPanels(const sim::StructuralMesh& structure,
                               const std::vector<int>& panels) {
    for (int index : panels) {
        if (index < 0 || static_cast<std::size_t>(index) >= structure.panels.size()) continue;
        const sim::PlatePanel& panel = structure.panels[static_cast<std::size_t>(index)];
        Quad quad;
        for (int k = 0; k < 4; ++k) quad.corner[k] = panel.corner[k];
        quad.normal = panel.normal();
        quad.centroid = panel.centroid();
        if (length2(quad.normal) < 0.5) continue;
        planeFrame(quad.normal, quad.right, quad.up);
        for (int k = 0; k < 4; ++k) {
            const Vec3 offset = quad.corner[k] - quad.centroid;
            quad.u[k] = dot(offset, quad.right);
            quad.v[k] = dot(offset, quad.up);
            quad.radius = std::max(quad.radius, length(offset));
        }
        tornArea_ += panel.area();
        torn_.push_back(quad);
    }
}

Vec3 HullDamage::displacementAt(const Vec3& rest) const {
    Vec3 total{0, 0, 0};

    for (const Zone& zone : zones_) {
        const Vec3 offset = rest - zone.centre;
        // Through-thickness guard. A hull surface wraps, so a point on the far
        // side of the ship can project into the same in-plane cell; without this
        // the port side would dent when the starboard side was struck.
        if (std::abs(dot(offset, zone.axis)) > params.panelReach) continue;
        const double u = dot(offset, zone.right), v = dot(offset, zone.up);
        const int iu = static_cast<int>(std::floor((u - zone.uLo) / zone.cell));
        const int iv = static_cast<int>(std::floor((v - zone.vLo) / zone.cell));
        if (iu < 0 || iu >= zone.nu || iv < 0 || iv >= zone.nv) continue;
        const std::size_t cell = static_cast<std::size_t>(iv) * static_cast<std::size_t>(zone.nu) +
                                 static_cast<std::size_t>(iu);
        for (std::uint32_t i = zone.bucketStart[cell]; i < zone.bucketStart[cell + 1]; ++i) {
            const std::size_t t = static_cast<std::size_t>(zone.bucketItem[i]) * 3;
            const std::uint32_t a = zone.tri[t], b = zone.tri[t + 1], c = zone.tri[t + 2];
            double w[3];
            if (!barycentric(u, v, zone.u[a], zone.v[a], zone.u[b], zone.v[b], zone.u[c],
                             zone.v[c], w))
                continue;
            total += zone.displacement[a] * w[0] + zone.displacement[b] * w[1] +
                     zone.displacement[c] * w[2];
            break;
        }
    }

    for (const Tent& tent : tents_) {
        // The falloff is on the **full 3D distance**, so the field reaches zero at
        // the edge of its own support and needs no guard to stop it -- and a tent
        // whose radius is smaller than a ship's half breadth therefore cannot
        // reach her far side, which is what the zone's axial guard is for. An
        // earlier version ran the falloff on the in-plane radius and cut the
        // support with a slab: continuous over the flat of a side, and a step at
        // the slab's face, which is the seam this file exists to prevent hiding
        // inside the mechanism meant to keep the other shell out of it.
        const double radius = length(rest - tent.centre);
        if (radius >= tent.halfSpan) continue;
        total += tent.outward * (-tent.depth * (1.0 - radius / tent.halfSpan));
    }
    return total;
}

double HullDamage::largestNodeDisplacement() const {
    double largest = 0;
    for (const Zone& zone : zones_) largest = std::max(largest, zone.largest);
    return largest;
}

double HullDamage::boundaryDisplacement() const {
    double largest = 0;
    for (const Zone& zone : zones_) largest = std::max(largest, zone.boundary);
    return largest;
}

double HullDamage::targetEdgeSize(const Vec3& at) const {
    double distance = kInfiniteDistance;
    for (const Zone& zone : zones_)
        distance = std::min(distance, std::max(0.0, length(at - zone.centre) - zone.radius));
    for (const Tent& tent : tents_)
        distance = std::min(distance, std::max(0.0, length(at - tent.centre) - tent.halfSpan));
    for (const Quad& quad : torn_)
        distance = std::min(distance, std::max(0.0, length(at - quad.centroid) - quad.radius));
    if (distance >= kInfiniteDistance) return params.coarseSize;
    return std::clamp(params.fineSize + params.grading * distance, params.fineSize,
                      params.coarseSize);
}

bool HullDamage::insideTear(const Vec3& at) const {
    for (const Quad& quad : torn_) {
        const Vec3 offset = at - quad.centroid;
        if (std::abs(dot(offset, quad.normal)) > params.panelReach) continue;
        const double u = dot(offset, quad.right), v = dot(offset, quad.up);
        double w[3];
        if (barycentric(u, v, quad.u[0], quad.v[0], quad.u[1], quad.v[1], quad.u[2], quad.v[2],
                        w) ||
            barycentric(u, v, quad.u[0], quad.v[0], quad.u[2], quad.v[2], quad.u[3], quad.v[3],
                        w))
            return true;
    }
    return false;
}

double HullDamage::distanceToTearEdge(const Vec3& at) const {
    double best = kInfiniteDistance;
    for (const Quad& quad : torn_) {
        const Vec3 offset = at - quad.centroid;
        if (std::abs(dot(offset, quad.normal)) > params.panelReach) continue;
        const double u = dot(offset, quad.right), v = dot(offset, quad.up);
        for (int k = 0; k < 4; ++k) {
            const int next = (k + 1) % 4;
            best = std::min(best, distanceToSegment(u, v, quad.u[k], quad.v[k], quad.u[next],
                                                    quad.v[next]));
        }
    }
    return best;
}

// --- Refinement -----------------------------------------------------------------

namespace {

// The recursive subdivision. Everything that keeps it crack-free is here:
// `needsSplit` reads only the two welded endpoint positions, the midpoint is
// interned against the welded endpoint pair, and every child is examined at
// depth + 1 whatever template produced it.
struct Refiner {
    const HullDamage* damage = nullptr;
    HullDamageParams params;

    std::vector<Vec3>* verts = nullptr;                // output positions
    std::vector<sim::Tri>* tris = nullptr;             // output triangles
    std::vector<Vec3> repPos;                          // welded representative positions
    std::vector<std::uint32_t> outRep;                 // output vertex -> representative
    std::map<std::uint64_t, std::uint32_t> midpoint_;  // welded edge -> output vertex

    int deepest = 0;
    bool depthLimited = false;

    bool needsSplit(std::uint32_t ra, std::uint32_t rb) const {
        // Sorted, so the two sides of an edge cannot even in principle evaluate a
        // different expression. `length(b - a)` and `(a + b) * 0.5` are already
        // symmetric bit for bit; this removes the need to argue it.
        if (ra > rb) std::swap(ra, rb);
        const Vec3& a = repPos[ra];
        const Vec3& b = repPos[rb];
        return length(b - a) > damage->targetEdgeSize((a + b) * 0.5);
    }

    std::uint32_t midpointOf(std::uint32_t ra, std::uint32_t rb) {
        if (ra > rb) std::swap(ra, rb);
        const std::uint64_t key = (static_cast<std::uint64_t>(ra) << 32) | rb;
        const auto found = midpoint_.find(key);
        if (found != midpoint_.end()) return found->second;
        const Vec3 mid = (repPos[ra] + repPos[rb]) * 0.5;
        const auto vertex = static_cast<std::uint32_t>(verts->size());
        verts->push_back(mid);
        outRep.push_back(static_cast<std::uint32_t>(repPos.size()));
        repPos.push_back(mid);
        midpoint_.emplace(key, vertex);
        return vertex;
    }

    void emit(std::uint32_t a, std::uint32_t b, std::uint32_t c, int depth) {
        deepest = std::max(deepest, depth);
        const std::uint32_t ra = outRep[a], rb = outRep[b], rc = outRep[c];
        const bool ab = needsSplit(ra, rb);
        const bool bc = needsSplit(rb, rc);
        const bool ca = needsSplit(rc, ra);
        const int count = (ab ? 1 : 0) + (bc ? 1 : 0) + (ca ? 1 : 0);

        if (count == 0) {
            tris->push_back({a, b, c});
            return;
        }
        if (depth >= params.maxDepth) {
            // Applied at the same depth on both sides of a shared edge -- the
            // recursion is uniform -- so stopping here cannot open a crack. It can
            // leave the mesh coarser than asked for, which is why it is reported.
            depthLimited = true;
            tris->push_back({a, b, c});
            return;
        }
        // The negative control: without the green cases a triangle whose neighbour
        // split an edge interpolates straight past the new vertex, which is exactly
        // a T-junction. Nothing in the engine takes this path.
        if (!params.stitch && count < 3) {
            tris->push_back({a, b, c});
            return;
        }

        if (count == 3) {
            const std::uint32_t mab = midpointOf(ra, rb);
            const std::uint32_t mbc = midpointOf(rb, rc);
            const std::uint32_t mca = midpointOf(rc, ra);
            emit(a, mab, mca, depth + 1);
            emit(mab, b, mbc, depth + 1);
            emit(mca, mbc, c, depth + 1);
            emit(mab, mbc, mca, depth + 1);
            return;
        }
        if (count == 1) {
            // Rotate so the split edge is ab, then two children sharing the new
            // midpoint. Winding is preserved in every rotation.
            if (bc) { emit2(b, c, a, depth); return; }
            if (ca) { emit2(c, a, b, depth); return; }
            emit2(a, b, c, depth);
            return;
        }
        // Two split edges: rotate so they are ab and bc.
        if (ab && bc) { emit3(a, b, c, depth); return; }
        if (bc && ca) { emit3(b, c, a, depth); return; }
        emit3(c, a, b, depth);
    }

    // ab split, ca and bc whole.
    void emit2(std::uint32_t a, std::uint32_t b, std::uint32_t c, int depth) {
        const std::uint32_t mab = midpointOf(outRep[a], outRep[b]);
        emit(a, mab, c, depth + 1);
        emit(mab, b, c, depth + 1);
    }

    // ab and bc split, ca whole.
    void emit3(std::uint32_t a, std::uint32_t b, std::uint32_t c, int depth) {
        const std::uint32_t mab = midpointOf(outRep[a], outRep[b]);
        const std::uint32_t mbc = midpointOf(outRep[b], outRep[c]);
        emit(mab, b, mbc, depth + 1);
        emit(a, mab, mbc, depth + 1);
        emit(a, mbc, c, depth + 1);
    }
};

}  // namespace

DamagedHull buildDamagedHull(const sim::TriMesh& hull, const HullDamage& damage) {
    const auto started = std::chrono::steady_clock::now();
    DamagedHull out;
    out.sourceTriangles = hull.tris.size();

    // The output vertex array starts as the hull's own and is never rewritten, so
    // an undamaged hull comes back bit-identical rather than nearly so.
    out.rest.verts = hull.verts;

    Refiner refiner;
    refiner.damage = &damage;
    refiner.params = damage.params;
    refiner.verts = &out.rest.verts;
    refiner.tris = &out.rest.tris;

    // Weld for the *decisions* only. Two triangles can reference distinct vertex
    // records holding the same position -- a mirrored or clipped hull carries them
    // all along its seam -- and if the two sides of that seam split their shared
    // edge differently, or address two different midpoints, the seam becomes a
    // crack. Welding the representatives removes both possibilities without
    // touching the emitted geometry.
    {
        std::map<std::array<long long, 3>, std::uint32_t> lookup;
        refiner.outRep.resize(hull.verts.size());
        for (std::size_t i = 0; i < hull.verts.size(); ++i) {
            const auto inserted =
                lookup.emplace(weldKey(hull.verts[i], 1e-6),
                               static_cast<std::uint32_t>(refiner.repPos.size()));
            refiner.outRep[i] = inserted.first->second;
            if (inserted.second) refiner.repPos.push_back(hull.verts[i]);
        }
    }

    out.rest.tris.reserve(hull.tris.size());
    for (const sim::Tri& tri : hull.tris) refiner.emit(tri.a, tri.b, tri.c, 0);
    out.deepestSplit = refiner.deepest;
    out.depthLimited = refiner.depthLimited;

    // Displace. `-0.0 + 0.0` is `+0.0`, so a zero displacement must not go through
    // the addition at all or "bit-identical when undamaged" quietly stops being
    // true for any vertex that carries a negative zero.
    out.deformed.verts.resize(out.rest.verts.size());
    for (std::size_t i = 0; i < out.rest.verts.size(); ++i) {
        const Vec3& rest = out.rest.verts[i];
        const Vec3 offset = damage.displacementAt(rest);
        if (offset.x == 0.0 && offset.y == 0.0 && offset.z == 0.0) {
            out.deformed.verts[i] = rest;
        } else {
            out.deformed.verts[i] = rest + offset;
            out.largestDisplacement = std::max(out.largestDisplacement, length(offset));
        }
    }

    // Cut. A torn panel is not recoloured: its triangles are removed, so what is
    // behind them is what the frame shows. The band around the hole takes the
    // exposed-metal material, decided on the *undeformed* centroid for the same
    // reason paint bands are.
    out.deformed.tris.reserve(out.rest.tris.size());
    out.exposed.reserve(out.rest.tris.size());
    std::vector<sim::Tri> kept;
    kept.reserve(out.rest.tris.size());
    for (const sim::Tri& tri : out.rest.tris) {
        const Vec3& a = out.rest.verts[tri.a];
        const Vec3& b = out.rest.verts[tri.b];
        const Vec3& c = out.rest.verts[tri.c];
        const Vec3 centroid = (a + b + c) / 3.0;
        if (damage.insideTear(centroid)) {
            ++out.droppedTriangles;
            out.holeArea += triangleArea(a, b, c);
            out.removed.push_back(tri);
            continue;
        }
        kept.push_back(tri);
        out.exposed.push_back(
            damage.distanceToTearEdge(centroid) <= damage.params.exposedWidth ? 1u : 0u);
    }
    out.rest.tris.swap(kept);
    out.deformed.tris = out.rest.tris;

    out.buildSeconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
    return out;
}

}  // namespace gpu
