// SPDX-License-Identifier: MIT
#include "section.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <functional>
#include <cmath>
#include <cstdio>
#include <map>
#include <set>
#include <string>
#include <utility>

namespace sim::section {
namespace {

using solidshell::kDof;
using solidshell::kGauss;
using solidshell::kNodes;

// --- Welding -------------------------------------------------------------------
//
// The same bucket-grid weld `zone.cpp` uses, and for the same reason: adjacent
// panels share their corners exactly, but the bilinear points along a shared edge
// agree bit for bit only when both panels number that edge the same way round, and
// the mirrored starboard panels do not. The probe reaches as far as the tolerance
// does whatever the cell size is -- hard-coding +/-1 ties the two together and
// leaves a crack down the middle of any mesh whose duplicates sit between 0.75 and
// 1 tolerance apart, which is the mutant that survived everything in `zone.cpp`.
class Welder {
public:
    explicit Welder(double tolerance)
        : tolerance_(tolerance),
          cell_(2.0 * tolerance),
          // The probe reaches as far as the tolerance does, whatever the cell size
          // is. Hard-coding +/-1 ties the two together silently.
          reach_(std::max<long long>(1, static_cast<long long>(std::ceil(tolerance / cell_)))) {}

    std::uint32_t weld(const Vec3& p, std::vector<Vec3>& points) {
        const long long bx = cellOf(p.x), by = cellOf(p.y), bz = cellOf(p.z);
        for (long long dx = -reach_; dx <= reach_; ++dx)
            for (long long dy = -reach_; dy <= reach_; ++dy)
                for (long long dz = -reach_; dz <= reach_; ++dz) {
                    auto found = buckets_.find(Key{bx + dx, by + dy, bz + dz});
                    if (found == buckets_.end()) continue;
                    for (std::uint32_t candidate : found->second)
                        if (length2(points[candidate] - p) <= tolerance_ * tolerance_)
                            return candidate;
                }
        const auto index = static_cast<std::uint32_t>(points.size());
        points.push_back(p);
        buckets_[Key{bx, by, bz}].push_back(index);
        return index;
    }

private:
    using Key = std::array<long long, 3>;
    long long cellOf(double v) const { return static_cast<long long>(std::floor(v / cell_)); }
    double tolerance_, cell_;
    long long reach_;
    std::map<Key, std::vector<std::uint32_t>> buckets_;
};

Vec3 quadNormal(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d) {
    return cross(c - a, d - b);  // twice the area vector; direction is what matters
}

double quadArea(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d) {
    return 0.5 * length(cross(b - a, c - a)) + 0.5 * length(cross(c - a, d - a));
}

// Squared distance from a point to a triangle. Used by the junction census: a free
// edge that is sitting *on* another plate rather than in fresh air is a joint the
// mesher could not make, and the only way to tell the two apart is to measure the
// distance to the other plate's surface rather than to its nodes -- a deck edge
// clipped to the hull lands in the middle of a shell face, so the nearest node can
// be half a bay away while the surface is at zero.
double distanceSquaredToTriangle(const Vec3& p, const Vec3& a, const Vec3& b, const Vec3& c) {
    const Vec3 ab = b - a, ac = c - a, ap = p - a;
    const double d1 = dot(ab, ap), d2 = dot(ac, ap);
    if (d1 <= 0 && d2 <= 0) return length2(ap);
    const Vec3 bp = p - b;
    const double d3 = dot(ab, bp), d4 = dot(ac, bp);
    if (d3 >= 0 && d4 <= d3) return length2(bp);
    const double vc = d1 * d4 - d3 * d2;
    if (vc <= 0 && d1 >= 0 && d3 <= 0) {
        const double v = d1 / (d1 - d3);
        return length2(ap - ab * v);
    }
    const Vec3 cp = p - c;
    const double d5 = dot(ab, cp), d6 = dot(ac, cp);
    if (d6 >= 0 && d5 <= d6) return length2(cp);
    const double vb = d5 * d2 - d1 * d6;
    if (vb <= 0 && d2 >= 0 && d6 <= 0) {
        const double w = d2 / (d2 - d6);
        return length2(ap - ac * w);
    }
    const double va = d3 * d6 - d5 * d4;
    if (va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0) {
        const double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return length2(bp - (c - b) * w);
    }
    const double denominator = 1.0 / (va + vb + vc);
    const double v = vb * denominator, w = vc * denominator;
    return length2(ap - ab * v - ac * w);
}

// The four bilinear shape functions of a quad wound
// (-1,-1) (+1,-1) (+1,+1) (-1,+1), which is the order a sub-quad's corners are
// emitted in. They sum to one everywhere, inside the face and outside it, which is
// what makes a tie built on them reproduce a rigid translation exactly.
void bilinear(double xi, double eta, double out[4]) {
    out[0] = 0.25 * (1 - xi) * (1 - eta);
    out[1] = 0.25 * (1 + xi) * (1 - eta);
    out[2] = 0.25 * (1 + xi) * (1 + eta);
    out[3] = 0.25 * (1 - xi) * (1 + eta);
}

// Where a point sits in a sub-quad's own coordinates: the (xi, eta) minimising the
// distance to the bilinear surface through its four corners, by Gauss-Newton on
// `dX/dxi . (X - p) = 0`.
//
// **Not clamped to the face**, and that is deliberate. A deck edge lands inside a
// shell face and comes back with |xi|, |eta| <= 1; two plates butting at a corner
// put a node half a plate thickness past the end of the other's mid-surface, and
// the honest answer there is a small extrapolation rather than a point moved onto
// the boundary. `SectionParams::junctionOvershoot` is what bounds it. Bilinear
// interpolation reproduces a linear field exactly at any (xi, eta), inside or out,
// so the tie stays exact for the plane-sections field either way -- what the
// overshoot costs is the sign of a weight, not the accuracy of the tie.
struct FaceCoordinate {
    double xi = 0, eta = 0;
    double distance = 0;   // m from the point to the surface it landed on
    double overshoot = 0;  // how far outside [-1, 1] the worse coordinate is
};

FaceCoordinate projectOntoQuad(const Vec3& p, const Vec3& a, const Vec3& b, const Vec3& c,
                               const Vec3& d) {
    FaceCoordinate out;
    for (int iteration = 0; iteration < 24; ++iteration) {
        double n[4];
        bilinear(out.xi, out.eta, n);
        const Vec3 x = a * n[0] + b * n[1] + c * n[2] + d * n[3];
        const Vec3 dxi = ((b - a) * (1 - out.eta) + (c - d) * (1 + out.eta)) * 0.25;
        const Vec3 deta = ((d - a) * (1 - out.xi) + (c - b) * (1 + out.xi)) * 0.25;
        const Vec3 r = x - p;
        const double g0 = dot(dxi, r), g1 = dot(deta, r);
        const double h00 = dot(dxi, dxi), h01 = dot(dxi, deta), h11 = dot(deta, deta);
        const double det = h00 * h11 - h01 * h01;
        if (!(std::abs(det) > 0)) break;  // a degenerate quad has no coordinates
        const double stepXi = -(h11 * g0 - h01 * g1) / det;
        const double stepEta = -(h00 * g1 - h01 * g0) / det;
        out.xi += stepXi;
        out.eta += stepEta;
        if (std::abs(stepXi) + std::abs(stepEta) < 1e-15) break;
    }
    double n[4];
    bilinear(out.xi, out.eta, n);
    out.distance = length(a * n[0] + b * n[1] + c * n[2] + d * n[3] - p);
    out.overshoot = std::max(0.0, std::max(std::abs(out.xi), std::abs(out.eta)) - 1.0);
    return out;
}

// A union-find over sub-quads, for the connected component count. It is not a
// diagnostic in the ordinary sense: a component touching neither cut plane is a
// mechanism in `K_ii`, and `reduction::Substructure` will *not* catch it -- its
// precondition check is geometric and about the interface, while a floating
// component leaves a tiny positive pivot that `BandedSpd::factor` accepts.
struct DisjointSet {
    std::vector<int> parent;
    explicit DisjointSet(std::size_t n) : parent(n) {
        for (std::size_t i = 0; i < n; ++i) parent[i] = static_cast<int>(i);
    }
    int find(int a) {
        while (parent[static_cast<std::size_t>(a)] != a) {
            parent[static_cast<std::size_t>(a)] =
                parent[static_cast<std::size_t>(parent[static_cast<std::size_t>(a)])];
            a = parent[static_cast<std::size_t>(a)];
        }
        return a;
    }
    void join(int a, int b) {
        a = find(a);
        b = find(b);
        if (a != b) parent[static_cast<std::size_t>(a)] = b;
    }
};

struct SubQuad {
    std::uint32_t node[4];
    int panel;
    int surface;
    // Which way round the surface walk decided this panel's winding should be read.
    // Carried per sub-quad rather than looked up per sub-quad, because the lookup is
    // over the candidate list and both are in the thousands.
    double orient;
};

}  // namespace

// --- Meshing -------------------------------------------------------------------

Section buildSection(const StructuralMesh& structure, const SectionParams& params) {
    Section section;
    const auto report = [&](std::string message) { section.problems.push_back(std::move(message)); };
    section.xFrom = params.xFrom;
    section.xTo = params.xTo;

    if (!(params.xTo > params.xFrom)) {
        report("the two cut planes are not in order: xFrom must be below xTo");
        return section;
    }
    if (params.subdivision < 1) {
        report("subdivision below one element per panel");
        return section;
    }

    // --- Which panels the section owns ------------------------------------------
    //
    // A panel with extent along x belongs when it lies wholly between the planes; a
    // panel a plane passes through is dropped and counted, because half a panel has
    // no corner set. A transverse plate -- a bulkhead, whose x extent is zero -- is
    // taken at the aft plane and not at the forward one, the same half-open rule
    // `sectionElements` uses so that two adjacent sections do not both own it.
    const double eps = params.weldTolerance;
    std::vector<int> candidate;
    for (std::size_t i = 0; i < structure.panels.size(); ++i) {
        const PlatePanel& p = structure.panels[i];
        if (p.role == PanelRole::Shell && !params.shell) continue;
        if (p.role == PanelRole::Deck && !params.deck) continue;
        if (p.role == PanelRole::Bulkhead && !params.bulkhead) continue;
        double lo = 1e300, hi = -1e300;
        for (int c = 0; c < 4; ++c) {
            lo = std::min(lo, p.corner[c].x);
            hi = std::max(hi, p.corner[c].x);
        }
        if (hi - lo <= eps) {  // a transverse plate: half-open at the aft plane
            if (lo >= params.xFrom - eps && lo < params.xTo - eps)
                candidate.push_back(static_cast<int>(i));
            continue;
        }
        if (lo >= params.xFrom - eps && hi <= params.xTo + eps) {
            candidate.push_back(static_cast<int>(i));
        } else if (hi > params.xFrom + eps && lo < params.xTo - eps) {
            ++section.straddlingPanels;
        }
    }
    if (section.straddlingPanels > 0)
        report(std::to_string(section.straddlingPanels) +
               " panels straddle a cut plane and were dropped, leaving a hole in the section:"
               " cut on a frame station, where the panel seams are");
    if (candidate.empty()) {
        report("no plating of the requested roles lies between the cut planes");
        return section;
    }

    // Material. One `StructuralMaterial` goes to `reduction::Substructure`, so a
    // section spanning two of them is only honest when they differ in nothing the
    // stiffness or the mass reads -- which on this ship is true, the weather deck's
    // mild steel differing from AH36 in yield alone.
    {
        std::set<int> used;
        for (int index : candidate)
            used.insert(structure.panels[static_cast<std::size_t>(index)].material);
        const int first = *used.begin();
        section.material = first >= 0 && first < static_cast<int>(structure.materials.size())
                               ? structure.materials[static_cast<std::size_t>(first)]
                               : ah36Steel();
        for (int index : used) {
            if (index == first) continue;
            const StructuralMaterial& other =
                index >= 0 && index < static_cast<int>(structure.materials.size())
                    ? structure.materials[static_cast<std::size_t>(index)]
                    : ah36Steel();
            if (other.youngsModulus != section.material.youngsModulus ||
                other.density != section.material.density ||
                other.poissonRatio != section.material.poissonRatio)
                report("the section spans materials '" + section.material.name + "' and '" +
                       other.name +
                       "' whose stiffness or density differ, and it is meshed entirely as the"
                       " first: a reduction takes one material");
        }
    }

    // --- Panel adjacency, by welded corner ---------------------------------------

    std::vector<Vec3> cornerPoint;
    Welder cornerWeld(params.weldTolerance);
    std::vector<std::array<std::uint32_t, 4>> panelCorner(candidate.size());
    for (std::size_t s = 0; s < candidate.size(); ++s) {
        const PlatePanel& p = structure.panels[static_cast<std::size_t>(candidate[s])];
        for (int c = 0; c < 4; ++c)
            panelCorner[s][static_cast<std::size_t>(c)] = cornerWeld.weld(p.corner[c], cornerPoint);
    }
    std::map<std::pair<std::uint32_t, std::uint32_t>, std::vector<std::size_t>> panelEdge;
    for (std::size_t s = 0; s < candidate.size(); ++s)
        for (int e = 0; e < 4; ++e) {
            const std::uint32_t a = panelCorner[s][static_cast<std::size_t>(e)];
            const std::uint32_t b = panelCorner[s][static_cast<std::size_t>((e + 1) % 4)];
            if (a == b) continue;  // a girth band the hull has closed to nothing
            panelEdge[{std::min(a, b), std::max(a, b)}].push_back(s);
        }

    // --- Surfaces -----------------------------------------------------------------
    //
    // A surface is a maximal set of panels reachable from one another across shared
    // edges without folding more than `foldLimit`. It is the unit the weld works
    // in: two panels in different surfaces never share a node however close their
    // corners are, because a shared node pair has **one** thickness direction and
    // two plates at an angle have two. See the header §1.
    std::vector<int> surfaceOf(candidate.size(), -1);
    std::vector<double> orientation(candidate.size(), 1.0);
    int surfaces = 0;
    int foldStops = 0;
    for (std::size_t seed = 0; seed < candidate.size(); ++seed) {
        if (surfaceOf[seed] >= 0) continue;
        surfaceOf[seed] = surfaces;
        orientation[seed] = 1.0;
        std::vector<std::size_t> order{seed};
        for (std::size_t head = 0; head < order.size(); ++head) {
            const std::size_t s = order[head];
            const Vec3 ns =
                normalize(structure.panels[static_cast<std::size_t>(candidate[s])].normal()) *
                orientation[s];
            for (int e = 0; e < 4; ++e) {
                const std::uint32_t a = panelCorner[s][static_cast<std::size_t>(e)];
                const std::uint32_t b = panelCorner[s][static_cast<std::size_t>((e + 1) % 4)];
                if (a == b) continue;
                auto users = panelEdge.find({std::min(a, b), std::max(a, b)});
                if (users == panelEdge.end()) continue;
                for (std::size_t other : users->second) {
                    if (other == s || surfaceOf[other] >= 0) continue;
                    const Vec3 nq =
                        normalize(structure.panels[static_cast<std::size_t>(candidate[other])].normal());
                    const double aligned = dot(nq, ns);
                    if (std::abs(aligned) < std::cos(params.foldLimit)) {
                        ++foldStops;
                        continue;
                    }
                    surfaceOf[other] = surfaces;
                    orientation[other] = aligned > 0 ? 1.0 : -1.0;
                    order.push_back(other);
                }
            }
        }
        ++surfaces;
    }
    section.surfaces = surfaces;
    if (foldStops > 0)
        report(std::to_string(foldStops) + " panel edges fold further than " +
               std::to_string(params.foldLimit) +
               " rad and were left unwelded: a solid-shell node pair carries one thickness"
               " direction, so a shared node at a corner would point between the two plates"
               " and be wrong for both");

    // --- Mid-surface grid ---------------------------------------------------------
    //
    // One `Welder` per weld class, all writing into one point array: a Welder only
    // ever finds points it put there itself, so keying the class on the surface --
    // and, under `ThicknessSeam::Split`, on the thickness as well -- is the whole of
    // how the mesher refuses to join what it should not join.
    const int n = params.subdivision;
    std::vector<Vec3> mid;
    std::map<std::pair<int, long long>, Welder> welders;
    std::vector<SubQuad> quads;
    quads.reserve(candidate.size() * static_cast<std::size_t>(n) * static_cast<std::size_t>(n));
    std::vector<std::uint32_t> grid(static_cast<std::size_t>(n + 1) * static_cast<std::size_t>(n + 1));

    for (std::size_t s = 0; s < candidate.size(); ++s) {
        const int panelIndex = candidate[s];
        const PlatePanel& p = structure.panels[static_cast<std::size_t>(panelIndex)];
        const long long thicknessKey = params.thicknessSeam == ThicknessSeam::Split
                                           ? std::llround(p.thickness * 1e9)
                                           : 0;
        auto found = welders.find({surfaceOf[s], thicknessKey});
        if (found == welders.end())
            found = welders.emplace(std::pair<int, long long>{surfaceOf[s], thicknessKey},
                                    Welder(params.weldTolerance))
                        .first;
        Welder& welder = found->second;
        for (int i = 0; i <= n; ++i)
            for (int j = 0; j <= n; ++j) {
                const double u = static_cast<double>(i) / n, v = static_cast<double>(j) / n;
                const Vec3 point = p.corner[0] * ((1.0 - u) * (1.0 - v)) +
                                   p.corner[1] * (u * (1.0 - v)) + p.corner[2] * (u * v) +
                                   p.corner[3] * ((1.0 - u) * v);
                grid[static_cast<std::size_t>(i) * static_cast<std::size_t>(n + 1) +
                     static_cast<std::size_t>(j)] = welder.weld(point, mid);
            }
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) {
                const auto at = [&](int a, int b) {
                    return grid[static_cast<std::size_t>(a) * static_cast<std::size_t>(n + 1) +
                                static_cast<std::size_t>(b)];
                };
                quads.push_back({{at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1)},
                                 panelIndex, surfaceOf[s], orientation[s]});
            }
    }
    if (quads.empty()) {
        report("the section meshed no elements");
        return section;
    }

    // --- Nodal normals and nodal thickness ----------------------------------------
    //
    // Area weighted, over the sub-quads that reached this node -- which, by the weld
    // classes above, are all on one surface and so all have a normal worth
    // averaging. The thickness is averaged the same way, and **that is the thickness
    // seam answer**: a node between a 12 mm and a 15.5 mm strake sits at 13.75 mm and
    // the taper lands inside the one element either side of it, measured below. See
    // the header §3 for why stopping instead -- which is what `zone.cpp` does -- is
    // not available to a section.
    std::vector<Vec3> nodeNormal(mid.size(), Vec3{0, 0, 0});
    std::vector<double> nodeThickness(mid.size(), 0.0), nodeWeight(mid.size(), 0.0);
    // Which surface a node belongs to. Single-valued by construction -- the weld
    // classes above are keyed on the surface, so a node is only ever reached by
    // sub-quads of one -- and the junction tie needs it to refuse to tie a free
    // edge to its own plating.
    std::vector<int> nodeSurface(mid.size(), -1);
    for (const SubQuad& q : quads) {
        const Vec3 raw = quadNormal(mid[q.node[0]], mid[q.node[1]], mid[q.node[2]], mid[q.node[3]]);
        const double area = quadArea(mid[q.node[0]], mid[q.node[1]], mid[q.node[2]], mid[q.node[3]]);
        const double thickness = structure.panels[static_cast<std::size_t>(q.panel)].thickness;
        // The panel's own winding decides the sub-quad's normal; the surface walk
        // above decided which way round that winding should be read.
        const double sign =
            dot(raw, structure.panels[static_cast<std::size_t>(q.panel)].normal()) >= 0 ? 1.0 : -1.0;
        const Vec3 oriented = raw * (sign * q.orient);
        for (int a = 0; a < 4; ++a) {
            const std::uint32_t node = q.node[static_cast<std::size_t>(a)];
            // **Once per sub-quad, not once per corner of it.** A collapsed
            // sub-quad -- see §7 -- names its apex twice, and counting it twice
            // would give that panel double weight in the area-weighted normal and
            // thickness at exactly the node where several panels meet. Where no
            // sub-quad is collapsed this loop is what it always was.
            bool repeat = false;
            for (int b = 0; b < a; ++b)
                if (q.node[static_cast<std::size_t>(b)] == node) repeat = true;
            if (repeat) continue;
            nodeNormal[node] += oriented;
            nodeThickness[node] += thickness * area;
            nodeWeight[node] += area;
            nodeSurface[node] = q.surface;
        }
    }
    for (std::size_t i = 0; i < mid.size(); ++i) {
        nodeNormal[i] = normalize(nodeNormal[i]);
        if (nodeWeight[i] > 0) nodeThickness[i] /= nodeWeight[i];
    }

    // --- Element edges, and which of them nothing is on the other side of ----------
    //
    // Computed here rather than with the rest of the topology because the junction
    // tie below needs it, and the tie has to be *paired* before the nodes are
    // numbered: a tie couples two nodes no element edge joins, and an ordering
    // chosen without knowing that is an ordering chosen for the wrong graph. It
    // cost the ferry hold a DOF half-bandwidth of 10 769 against 149 to find out.
    std::map<std::pair<std::uint32_t, std::uint32_t>, int> edgeUse;
    for (const SubQuad& q : quads)
        for (int e = 0; e < 4; ++e) {
            const std::uint32_t a = q.node[static_cast<std::size_t>(e)];
            const std::uint32_t b = q.node[static_cast<std::size_t>((e + 1) % 4)];
            if (a == b) continue;
            ++edgeUse[{std::min(a, b), std::max(a, b)}];
        }

    // Which mid-surface nodes lie on a cut plane. The interface is chosen on the
    // **mid-surface** and both extruded nodes of a pair taken, which is exact.
    // Choosing it on the nodes instead -- which is what `reduction::nodesNearPlanes`
    // at its default tolerance does -- keeps only the node of each pair that
    // happened to land on the plane wherever the plating's normal leans out of the
    // transverse direction, and half an interface is not a cut, it is a hinge.
    std::vector<std::uint8_t> onPlane(mid.size(), 0u);
    for (std::uint32_t i = 0; i < mid.size(); ++i) {
        const bool aft = std::abs(mid[i].x - params.xFrom) <= params.planeTolerance;
        const bool forward = std::abs(mid[i].x - params.xTo) <= params.planeTolerance;
        if (aft || forward) onPlane[i] = aft ? 1u : 2u;
    }

    // A grid over the sub-quads, so that "which other surface is this node lying
    // on" is a local search. Shared by the junction tie below and by the free-edge
    // census much further down, because the two have to agree about what counts as
    // lying on another surface or the census would report an edge as unjoined that
    // the tie had just joined.
    const double cell = 1.0;
    std::map<std::array<long long, 3>, std::vector<std::size_t>> bucket;
    const auto cellOf = [&](double v) { return static_cast<long long>(std::floor(v / cell)); };
    for (std::size_t k = 0; k < quads.size(); ++k) {
        Vec3 lo = mid[quads[k].node[0]], hi = lo;
        for (int c = 1; c < 4; ++c) {
            const Vec3& p = mid[quads[k].node[c]];
            lo = Vec3{std::min(lo.x, p.x), std::min(lo.y, p.y), std::min(lo.z, p.z)};
            hi = Vec3{std::max(hi.x, p.x), std::max(hi.y, p.y), std::max(hi.z, p.z)};
        }
        for (long long ix = cellOf(lo.x - params.junctionTolerance);
             ix <= cellOf(hi.x + params.junctionTolerance); ++ix)
            for (long long iy = cellOf(lo.y - params.junctionTolerance);
                 iy <= cellOf(hi.y + params.junctionTolerance); ++iy)
                for (long long iz = cellOf(lo.z - params.junctionTolerance);
                     iz <= cellOf(hi.z + params.junctionTolerance); ++iz)
                    bucket[{ix, iy, iz}].push_back(k);
    }

    // --- The junction tie, section 5, paired ---------------------------------------
    //
    // A free-edge node lying on another surface is tied to the point of that surface
    // it lands in: bilinear in the master face's two in-plane coordinates, and
    // through the master's thickness by the weight `constraint.hpp` §1 derives.
    // Eight masters, three constraints -- one per axis -- and the slave keeps no
    // unknown of its own.
    //
    // **This is not a weld and the difference is the whole point.** A weld merges
    // the two node pairs into one, which has one thickness direction where two
    // plates at an angle need two, and pays for it in steel: 9.4% of `EA` on the
    // box. A tie leaves both pairs where they are and moves no geometry at all; it
    // constrains the deck's *displacement* to the shell's, which is what a weld does
    // physically.
    //
    // The pairing happens here and the constraints are written after the numbering,
    // because the numbering has to see them.
    struct PairedTie {
        std::uint32_t slave;    // mid-surface node
        std::size_t quad[2];    // the master face each of its two extruded nodes landed in
        double xi[2], eta[2];   // where in that face
        // The through-thickness split of the master, `constraint::tieWeight`.
        // Computed once, when the pairing decides whether the tie is admissible at
        // all, and carried to the write-out -- because a tie refused on one number
        // and written with another is exactly the kind of disagreement that never
        // shows up in a total.
        double weight[2];
    };
    // Where a slave's extruded node lands on a master face, and how the master's
    // thickness splits it. Shared by the pairing and the write-out for the same
    // reason.
    const auto tieOnFace = [&](std::uint32_t slave, int face, std::size_t quad, double xi,
                               double eta, double shape[4]) {
        const SubQuad& q = quads[quad];
        bilinear(xi, eta, shape);
        Vec3 surface{0, 0, 0}, normal{0, 0, 0};
        double thickness = 0;
        for (int a = 0; a < 4; ++a) {
            surface += mid[q.node[a]] * shape[a];
            normal += nodeNormal[q.node[a]] * shape[a];
            thickness += nodeThickness[q.node[a]] * shape[a];
        }
        normal = normalize(normal);
        const Vec3 at =
            mid[slave] + nodeNormal[slave] * ((face == 0 ? -0.5 : 0.5) * nodeThickness[slave]);
        return constraint::tieWeight(dot(at - surface, normal), thickness);
    };
    std::vector<PairedTie> paired;
    std::vector<std::uint8_t> tieRole(mid.size(), 0u);  // 1 slave, 2 master
    if (params.junctions) {
        std::vector<std::uint8_t> onFreeEdge(mid.size(), 0u);
        for (const auto& [edge, count] : edgeUse) {
            if (count != 1) continue;
            onFreeEdge[edge.first] = 1u;
            onFreeEdge[edge.second] = 1u;
        }
        // Ascending mid-surface order, so which side of a junction becomes the slave
        // is a property of the mesh rather than of a map's iteration.
        for (std::uint32_t i = 0; i < mid.size(); ++i) {
            if (!onFreeEdge[i] || nodeWeight[i] <= 0) continue;
            // The two mesh nodes of the pair, and the two faces they land in --
            // separately, because they are a plate thickness apart along the
            // *master's* in-plane direction and their difference is exactly what
            // carries this plate's rotation about the junction line into the other's
            // own surface slope. Tying both to one point would clamp that rotation.
            const Vec3 through[2] = {mid[i] - nodeNormal[i] * (0.5 * nodeThickness[i]),
                                     mid[i] + nodeNormal[i] * (0.5 * nodeThickness[i])};
            PairedTie tie{i, {quads.size(), quads.size()}, {0, 0}, {0, 0}, {0.5, 0.5}};
            FaceCoordinate landed[2];
            for (int face = 0; face < 2; ++face) {
                auto found = bucket.find(
                    {cellOf(through[face].x), cellOf(through[face].y), cellOf(through[face].z)});
                if (found == bucket.end()) continue;
                for (std::size_t k : found->second) {
                    const SubQuad& q = quads[k];
                    if (q.surface == nodeSurface[i]) continue;  // its own plating
                    const FaceCoordinate at =
                        projectOntoQuad(through[face], mid[q.node[0]], mid[q.node[1]],
                                        mid[q.node[2]], mid[q.node[3]]);
                    if (at.distance > params.junctionTolerance) continue;
                    // The face the point is most *inside* wins, not the nearest one.
                    // Along a butt corner several faces are the same distance away
                    // and only one of them contains the point; picking by distance
                    // alone would extrapolate off the end of a neighbour.
                    if (tie.quad[face] == quads.size() ||
                        at.overshoot < landed[face].overshoot - 1e-12 ||
                        (at.overshoot < landed[face].overshoot + 1e-12 &&
                         at.distance < landed[face].distance)) {
                        tie.quad[face] = k;
                        landed[face] = at;
                    }
                }
            }
            if (tie.quad[0] == quads.size() || tie.quad[1] == quads.size()) continue;
            if (onPlane[i] != 0u) {
                ++section.junctionsOnInterface;
                continue;
            }
            if (landed[0].overshoot > params.junctionOvershoot ||
                landed[1].overshoot > params.junctionOvershoot) {
                ++section.junctionsOutsideFace;
                continue;
            }
            // And the through-thickness half of the same guarantee. A weight
            // outside [0, 1] is normal -- a deck edge sits on the far face of the
            // shell -- but one that takes more than `junctionWeightLimit` shares of
            // the slave's mass off a master face leaves `T^T M T` with a
            // non-positive diagonal, and `reduction::Substructure` then refuses the
            // whole section rather than this one junction. Better to leave the
            // junction open and say which.
            bool overWeight = false;
            for (int face = 0; face < 2; ++face) {
                double shape[4];
                tie.weight[face] = tieOnFace(i, face, tie.quad[face], landed[face].xi,
                                             landed[face].eta, shape);
                const double share = std::min(tie.weight[face], 1.0 - tie.weight[face]);
                if (share < -params.junctionWeightLimit) overWeight = true;
            }
            if (overWeight) {
                ++section.junctionsThroughThickness;
                continue;
            }
            // Already the master side of a junction made from the other plate. The
            // joint is closed -- one constraint joins two edges -- and adding a
            // second tie the other way round would be the same statement twice, in a
            // form `solidshell::DofExpansion` would have to refuse as a chain.
            if (tieRole[i] == 2u) continue;
            // A master of mine is already somebody's slave. That **is** a chain, and
            // it is refused rather than composed: the mesher would rather leave a
            // junction open and say so than resolve one silently.
            bool chained = false;
            for (int face = 0; face < 2 && !chained; ++face)
                for (std::uint32_t node : quads[tie.quad[face]].node)
                    chained = chained || tieRole[node] == 1u;
            if (chained) {
                ++section.junctionsChained;
                continue;
            }
            for (int face = 0; face < 2; ++face) {
                tie.xi[face] = landed[face].xi;
                tie.eta[face] = landed[face].eta;
                section.worstJunctionOvershoot =
                    std::max(section.worstJunctionOvershoot, landed[face].overshoot);
                for (std::uint32_t node : quads[tie.quad[face]].node) tieRole[node] = 2u;
            }
            tieRole[i] = 1u;
            paired.push_back(tie);
            ++section.junctionTies;
        }
        if (section.junctionsChained > 0)
            report(std::to_string(section.junctionsChained) +
                   " junction nodes could not be tied because a master of theirs is already a"
                   " slave: chained constraints are refused rather than composed, so those"
                   " junctions are still open");
        if (section.junctionsOnInterface > 0)
            report(std::to_string(section.junctionsOnInterface) +
                   " junction nodes lie on a cut plane and were left untied: the interface is"
                   " prescribed, and a prescribed degree of freedom cannot also be a function of"
                   " others");
        if (section.junctionsOutsideFace > 0)
            report(std::to_string(section.junctionsOutsideFace) +
                   " junction nodes land further outside a master face than"
                   " SectionParams::junctionOvershoot allows and were left untied");
        if (section.junctionsThroughThickness > 0)
            report(std::to_string(section.junctionsThroughThickness) +
                   " junction nodes sit so far off their master's mid-surface that the tie would"
                   " take more than SectionParams::junctionWeightLimit of the slave's mass off one"
                   " master face; they were left untied rather than left to make the section's"
                   " condensed mass non-positive");
    }

    // --- Node numbering -----------------------------------------------------------
    //
    // `solidshell::solveStatic` numbers its free degrees of freedom in the **mesh's
    // own order** and takes whatever band that delivers -- it has no renumbering
    // pass, unlike `reduction::Substructure`, which reverse-Cuthill-McKees its
    // interior and does not care what arrives.
    //
    // So the ordering matters here, and *which* ordering is right is a property of
    // the section rather than a rule: a hold-length piece of the reference ferry has
    // twelve stations along x and hundreds round the girth, where numbering x
    // fastest is what makes a slab of the band twelve nodes deep; a long slender box
    // is the other way round, and getting it backwards cost a 6 240 DOF solve **9.53
    // seconds instead of 0.06**, while leaving the ferry hold's half-bandwidth at
    // 1 382 instead of 146 cost it 5.3 s instead of 0.14. Three orderings are built
    // -- the two lexicographic ones and `reduction::bandwidthReducingOrder` -- and
    // the narrowest is kept.
    // Comparing is free, a bandwidth being one pass over the sub-quads, and it makes
    // the ordering incapable of being a pessimisation, which is the same argument
    // `reduction.cpp` makes for keeping the better of RCM and natural.
    std::vector<std::uint32_t> rank(mid.size());
    {
        const auto ranking = [&](bool alongFastest) {
            std::vector<std::uint32_t> byPosition(mid.size());
            for (std::uint32_t i = 0; i < mid.size(); ++i) byPosition[i] = i;
            const auto key = [&](std::uint32_t i) {
                const long long x = std::llround(mid[i].x * 1e6);
                const long long y = std::llround(mid[i].y * 1e6);
                const long long z = std::llround(mid[i].z * 1e6);
                return alongFastest ? std::array<long long, 3>{z, y, x}
                                    : std::array<long long, 3>{x, z, y};
            };
            std::sort(byPosition.begin(), byPosition.end(),
                      [&](std::uint32_t a, std::uint32_t b) { return key(a) < key(b); });
            std::vector<std::uint32_t> out(mid.size());
            for (std::uint32_t r = 0; r < byPosition.size(); ++r) out[byPosition[r]] = r;
            return out;
        };
        // **The spread is taken over the tied graph, not over the element graph.** An
        // element with a tied corner reaches that corner's eight masters, which are on
        // another surface entirely, so an ordering scored on the sub-quads alone is
        // scored on a graph the solver does not have. Measured on the ferry hold: the
        // tie takes the assembled DOF half-bandwidth from 149 to 10 769 -- essentially
        // dense on 14 040 degrees of freedom -- unless the numbering is chosen knowing
        // about it.
        const auto reachOf = [&](const std::vector<std::uint32_t>& order, std::uint32_t node,
                                 std::uint32_t& lo, std::uint32_t& hi) {
            lo = std::min(lo, order[node]);
            hi = std::max(hi, order[node]);
        };
        std::vector<std::vector<std::uint32_t>> tiedTo(mid.size());
        for (const PairedTie& tie : paired)
            for (int face = 0; face < 2; ++face)
                for (std::uint32_t node : quads[tie.quad[face]].node)
                    tiedTo[tie.slave].push_back(node);
        const auto spreadOf = [&](const std::vector<std::uint32_t>& order) {
            std::size_t worst = 0;
            for (const SubQuad& q : quads) {
                std::uint32_t lo = order[q.node[0]], hi = lo;
                for (int c = 0; c < 4; ++c) {
                    reachOf(order, q.node[c], lo, hi);
                    for (std::uint32_t master : tiedTo[q.node[c]]) reachOf(order, master, lo, hi);
                }
                worst = std::max(worst, static_cast<std::size_t>(hi - lo));
            }
            return worst;
        };
        std::vector<std::vector<std::uint32_t>> adjacency(mid.size());
        {
            std::set<std::pair<std::uint32_t, std::uint32_t>> seen;
            const auto join = [&](std::uint32_t a, std::uint32_t b) {
                if (a == b || !seen.insert({std::min(a, b), std::max(a, b)}).second) return;
                adjacency[a].push_back(b);
                adjacency[b].push_back(a);
            };
            for (const SubQuad& q : quads)
                for (int e = 0; e < 4; ++e)
                    for (int f = e + 1; f < 4; ++f)
                        join(q.node[static_cast<std::size_t>(e)], q.node[static_cast<std::size_t>(f)]);
            // The tie's own edges, so the Cuthill-McKee sweep crosses a junction the
            // way it crosses an element and numbers the two surfaces together instead
            // of one after the other.
            for (const PairedTie& tie : paired)
                for (int face = 0; face < 2; ++face)
                    for (std::uint32_t node : quads[tie.quad[face]].node) join(tie.slave, node);
        }
        std::vector<std::uint32_t> cuthill(mid.size());
        {
            const std::vector<std::uint32_t> order = reduction::bandwidthReducingOrder(adjacency);
            if (order.size() == mid.size())
                for (std::uint32_t r = 0; r < order.size(); ++r) cuthill[order[r]] = r;
            else
                for (std::uint32_t r = 0; r < mid.size(); ++r) cuthill[r] = r;
        }
        const std::vector<std::uint32_t> candidates[3] = {ranking(true), ranking(false), cuthill};
        std::size_t best = spreadOf(candidates[0]);
        int chosen = 0;
        for (int k = 1; k < 3; ++k) {
            const std::size_t spread = spreadOf(candidates[static_cast<std::size_t>(k)]);
            if (spread < best) {
                best = spread;
                chosen = k;
            }
        }
        rank = candidates[static_cast<std::size_t>(chosen)];
        // A node pair is two mesh nodes, so a mid-surface spread of `s` is a nodal
        // spread of `2s + 1` and a DOF half-bandwidth of `3 * (2s + 1) - 1`.
        section.halfBandwidth = 3 * (2 * best + 1) - 1;
    }

    section.mesh.position.assign(mid.size() * 2 * 3, 0.0);
    section.outerFace.assign(mid.size() * 2, 0u);
    section.nodeThickness.assign(mid.size() * 2, 0.0);
    for (std::uint32_t i = 0; i < mid.size(); ++i) {
        const double half = 0.5 * nodeThickness[i];
        const Vec3 lower = mid[i] - nodeNormal[i] * half;
        const Vec3 upper = mid[i] + nodeNormal[i] * half;
        const std::size_t low = static_cast<std::size_t>(rank[i]) * 2;
        for (int k = 0; k < 3; ++k) {
            section.mesh.position[low * 3 + static_cast<std::size_t>(k)] = lower[k];
            section.mesh.position[(low + 1) * 3 + static_cast<std::size_t>(k)] = upper[k];
        }
        section.outerFace[low + 1] = 1u;
        section.nodeThickness[low] = nodeThickness[i];
        section.nodeThickness[low + 1] = nodeThickness[i];
    }

    // --- Elements ------------------------------------------------------------------

    section.worstJacobian = 1e300;
    double worstSpread = 0;
    DisjointSet connected(mid.size());
    for (const SubQuad& q : quads) {
        const Vec3& a = mid[q.node[0]];
        const Vec3& b = mid[q.node[1]];
        const Vec3& c = mid[q.node[2]];
        const Vec3& d = mid[q.node[3]];
        const double area = quadArea(a, b, c, d);
        if (!(area > 0)) continue;

        std::uint32_t corner[4] = {q.node[0], q.node[1], q.node[2], q.node[3]};
        // Node ordering is a contract: 0-3 on the -zeta face, counter-clockwise seen
        // from +zeta, with 4 above 0. zeta is the thickness direction, so a quad
        // wound the other way has to be reversed or the assumed strains would cure
        // locking *across* the plate instead of through it.
        Vec3 meanNormal{0, 0, 0};
        for (std::uint32_t node : corner) meanNormal += nodeNormal[node];
        meanNormal = normalize(meanNormal);
        if (dot(quadNormal(a, b, c, d), meanNormal) < 0) std::swap(corner[1], corner[3]);

        double spread = 0;
        for (std::uint32_t node : corner)
            spread = std::max(spread, length(nodeNormal[node] - meanNormal));
        worstSpread = std::max(worstSpread, spread);
        if (spread > params.normalSpreadWarning) ++section.distortedElements;

        double thin = 1e300, thick = 0;
        for (std::uint32_t node : corner) {
            thin = std::min(thin, nodeThickness[node]);
            thick = std::max(thick, nodeThickness[node]);
        }
        if (thick - thin > 1e-12) {
            ++section.taperedElements;
            const double mean = 0.5 * (thin + thick);
            if (mean > 0) section.worstTaper = std::max(section.worstTaper, (thick - thin) / mean);
        }

        double edge[4];
        for (int e = 0; e < 4; ++e)
            edge[e] = length(mid[corner[static_cast<std::size_t>((e + 1) % 4)]] -
                             mid[corner[static_cast<std::size_t>(e)]]);
        const double shortest = *std::min_element(edge, edge + 4);
        const double longest = *std::max_element(edge, edge + 4);
        if (shortest > 0) section.worstAspect = std::max(section.worstAspect, longest / shortest);

        for (int k = 0; k < 4; ++k)
            section.mesh.index.push_back(static_cast<std::uint32_t>(rank[corner[k]]) * 2);
        for (int k = 0; k < 4; ++k)
            section.mesh.index.push_back(static_cast<std::uint32_t>(rank[corner[k]]) * 2 + 1);

        double nodes[kDof];
        section.mesh.gather(section.mesh.elementCount() - 1, section.mesh.position, nodes);
        // The Gauss-point determinant, not the nodal one, because that is the one a
        // solve requires and the one whose sign is fatal. A collapsed element's
        // nodal determinant is zero by construction and says nothing about it; it
        // is counted instead. See `solidshell::ElementShape` and §7.
        const solidshell::ElementShape shape = solidshell::elementShape(nodes);
        section.worstJacobian = std::min(section.worstJacobian, shape.smallestGauss);
        if (shape.collapsedNodes > 0) ++section.collapsedElements;
        if (!shape.integrable) ++section.invertedElements;

        for (int k = 1; k < 4; ++k)
            connected.join(static_cast<int>(corner[0]), static_cast<int>(corner[k]));

        section.panelOf.push_back(q.panel);
        section.elementArea.push_back(area);
        section.area += area;
        section.plateMass += area * 0.25 *
                             (nodeThickness[corner[0]] + nodeThickness[corner[1]] +
                              nodeThickness[corner[2]] + nodeThickness[corner[3]]) *
                             section.material.density;
    }
    section.worstNormalSpread = worstSpread;
    section.spuriousStiffness = 90.0 * worstSpread * worstSpread;
    section.taperStiffness = 90.0 * section.worstTaper * section.worstTaper;
    if (worstSpread > params.normalSpreadWarning)
        report("the plating turns " + std::to_string(worstSpread) + " rad across " +
               std::to_string(section.distortedElements) + " of " +
               std::to_string(section.elementCount()) +
               " elements, so their faces are that far from parallel and the worst is about " +
               std::to_string(100.0 * section.spuriousStiffness) + "% too stiff in bending");
    if (section.worstTaper > params.taperWarning)
        report(std::to_string(section.taperedElements) + " of " +
               std::to_string(section.elementCount()) +
               " elements taper across a plate thickness seam, the worst by dt/t = " +
               std::to_string(section.worstTaper) + ", which is about " +
               std::to_string(100.0 * section.taperStiffness) +
               "% excess bending stiffness in those elements alone. It is the plating's own"
               " bending about its own mid-surface, which is 0.03% of a hull girder's second"
               " moment; ThicknessSeam::Split is the control that measures it");
    if (section.invertedElements > 0)
        report(std::to_string(section.invertedElements) + " of " +
               std::to_string(section.elementCount()) +
               " elements came out inverted; nothing computed on this section means anything");
    else if (section.collapsedElements > 0)
        report(std::to_string(section.collapsedElements) + " of " +
               std::to_string(section.elementCount()) +
               " elements are collapsed hexahedra -- the wedge a plate panel with two coincident"
               " corners extrudes to, which is what a bulkhead running into the keel or a deck"
               " strake running out at the stem is. They integrate; the worst Gauss determinant"
               " over the whole section is " +
               std::to_string(section.worstJacobian));

    {
        std::set<int> covered;
        for (int index : section.panelOf) covered.insert(index);
        section.panels.assign(covered.begin(), covered.end());
    }

    // --- The interface: the two cut planes -----------------------------------------
    //
    // `onPlane` was decided on the mid-surface further up, before the numbering; this
    // is the same set written out in the numbering's own terms, both extruded nodes
    // of every pair.
    section.mesh.fixed.assign(section.mesh.nodeCount() * 3, 0u);
    section.mesh.prescribed.assign(section.mesh.nodeCount() * 3, 0.0);
    for (std::uint32_t i = 0; i < mid.size(); ++i) {
        if (onPlane[i] == 0u) continue;
        const auto low = static_cast<std::uint32_t>(rank[i]) * 2;
        std::vector<std::uint32_t>& into =
            onPlane[i] == 1u ? section.aftNodes : section.forwardNodes;
        into.push_back(low);
        into.push_back(low + 1);
        section.interfaceNodes.push_back(low);
        section.interfaceNodes.push_back(low + 1);
    }
    std::sort(section.aftNodes.begin(), section.aftNodes.end());
    std::sort(section.forwardNodes.begin(), section.forwardNodes.end());
    std::sort(section.interfaceNodes.begin(), section.interfaceNodes.end());
    if (section.aftNodes.empty() || section.forwardNodes.empty())
        report("one of the cut planes carries no nodes, so the section is attached to the ship"
               " at one end only and its reduction has no interface to be exact at");

    // --- The junction tie, written out ---------------------------------------------
    //
    // The pairing is above; this turns each pair into three constraints per extruded
    // node, one per axis, over the eight mesh nodes of the master face. The weights
    // are `N_a(xi, eta)` times the through-thickness split `constraint::tieWeight`
    // gives, and they sum to one -- so a rigid translation of the whole section is
    // reproduced exactly, and so is a rigid rotation, because bilinear interpolation
    // of a linear field is that field.
    for (const PairedTie& tie : paired) {
        const std::size_t low = static_cast<std::size_t>(rank[tie.slave]) * 2;
        for (int face = 0; face < 2; ++face) {
            const SubQuad& q = quads[tie.quad[face]];
            double shape[4];
            bilinear(tie.xi[face], tie.eta[face], shape);
            // The weight the pairing admitted this tie on, not a second computation
            // of it.
            const double weight = tie.weight[face];
            section.worstJunctionWeight = std::max(
                section.worstJunctionWeight, std::max(std::abs(weight), std::abs(1.0 - weight)));
            for (int k = 0; k < 3; ++k) {
                solidshell::Mpc mpc;
                mpc.slave = static_cast<std::uint32_t>(
                    (low + static_cast<std::size_t>(face)) * 3 + static_cast<std::size_t>(k));
                mpc.master.reserve(8);
                mpc.weight.reserve(8);
                for (int a = 0; a < 4; ++a) {
                    const auto pair = static_cast<std::uint32_t>(rank[q.node[a]]) * 2;
                    mpc.master.push_back(pair * 3 + static_cast<std::uint32_t>(k));
                    mpc.weight.push_back(shape[a] * (1.0 - weight));
                    mpc.master.push_back((pair + 1) * 3 + static_cast<std::uint32_t>(k));
                    mpc.weight.push_back(shape[a] * weight);
                }
                section.attachment.constrained.push_back(std::move(mpc));
            }
        }
        for (int face = 0; face < 2; ++face)
            for (std::uint32_t node : quads[tie.quad[face]].node)
                connected.join(static_cast<int>(tie.slave), static_cast<int>(node));
    }

    // --- Topology: components, free edges, and junctions that are still open --------

    // Which surface each free edge belongs to, so the census does not report an edge
    // as sitting on its own plating.
    std::map<std::pair<std::uint32_t, std::uint32_t>, int> edgeSurface;
    for (const SubQuad& q : quads)
        for (int e = 0; e < 4; ++e) {
            const std::uint32_t a = q.node[static_cast<std::size_t>(e)];
            const std::uint32_t b = q.node[static_cast<std::size_t>((e + 1) % 4)];
            if (a == b) continue;
            edgeSurface[{std::min(a, b), std::max(a, b)}] = q.surface;
        }

    // Components, compacted to 0..components-1 and recorded per node. A caller has
    // to know: `applyBeamLoad` prescribes only `u_x` at the interface and every
    // *separate* piece then keeps three rigid body motions of its own, so the
    // restraints that remove them are per component and not per section. Getting
    // that wrong is not a small error -- the whole solve goes singular, which is at
    // least loud.
    section.componentOf.assign(section.mesh.nodeCount(), -1);
    {
        std::map<int, int> compact;
        std::map<int, int> reaches;  // 0 neither plane, 1 one, 2 both
        for (std::uint32_t i = 0; i < mid.size(); ++i) {
            if (nodeWeight[i] <= 0) continue;
            const int root = connected.find(static_cast<int>(i));
            auto found = compact.find(root);
            if (found == compact.end())
                found = compact.emplace(root, static_cast<int>(compact.size())).first;
            const auto low = static_cast<std::size_t>(rank[i]) * 2;
            section.componentOf[low] = found->second;
            section.componentOf[low + 1] = found->second;
            if (onPlane[i] == 1u) reaches[root] |= 1;
            if (onPlane[i] == 2u) reaches[root] |= 2;
        }
        section.components = static_cast<int>(compact.size());
        for (const auto& [root, index] : compact) {
            const int touched = reaches.count(root) ? reaches[root] : 0;
            if (touched == 0) ++section.floatingComponents;
            if (touched == 3) ++section.spanningComponents;
        }
    }
    if (section.floatingComponents > 0)
        report(std::to_string(section.floatingComponents) + " of " +
               std::to_string(section.components) +
               " connected components touch neither cut plane. Each is a free rigid body inside"
               " K_ii, which factors to a tiny positive pivot rather than a zero one, so"
               " reduction::Substructure will accept it and return nonsense -- drop the role"
               " or move the planes");
    if (section.components > section.spanningComponents)
        report(std::to_string(section.components - section.spanningComponents) + " of " +
               std::to_string(section.components) +
               " connected components do not reach both cut planes, so they carry none of the"
               " section's axial or bending stiffness and are held at one end or not at all");
    if (section.components > 1)
        report("the section is in " + std::to_string(section.components) +
               " disconnected pieces: either the input shares no corner across the junction --"
               " which is what makeStructuralMesh delivers, see section.hpp section 1 -- or the"
               " fold there is past SectionParams::foldLimit and welding it would have put one"
               " thickness direction between two plates");

    // Free edges, and how many of them are lying on another surface they are not
    // joined to. A grid over the sub-quads keeps the search local; a free edge is
    // only interesting when it is *not* on a cut plane, because those are the cut.
    for (const auto& [edge, count] : edgeUse) {
        if (count != 1) continue;
        if (onPlane[edge.first] && onPlane[edge.second] &&
            onPlane[edge.first] == onPlane[edge.second])
            continue;  // this edge is the cut
        const double span = length(mid[edge.second] - mid[edge.first]);
        section.freeEdgeLength += span;
        const Vec3 centre = (mid[edge.first] + mid[edge.second]) * 0.5;
        const int own = edgeSurface[edge];
        double nearest = 1e300;
        auto found = bucket.find({cellOf(centre.x), cellOf(centre.y), cellOf(centre.z)});
        if (found != bucket.end())
            for (std::size_t k : found->second) {
                if (quads[k].surface == own) continue;
                const Vec3& a = mid[quads[k].node[0]];
                const Vec3& b = mid[quads[k].node[1]];
                const Vec3& c = mid[quads[k].node[2]];
                const Vec3& d = mid[quads[k].node[3]];
                nearest = std::min(nearest, distanceSquaredToTriangle(centre, a, b, c));
                nearest = std::min(nearest, distanceSquaredToTriangle(centre, a, c, d));
            }
        const double gap = nearest < 1e300 ? std::sqrt(nearest) : 1e300;
        if (gap <= params.junctionTolerance) {
            section.junctionEdges += span;
            section.worstJunctionGap = std::max(section.worstJunctionGap, gap);
            // An edge is joined when **both** of its nodes take part in a tie, on
            // either side of it. Both, because a tie is a statement about a node and
            // one end tied leaves the edge free to open at the other -- counting a
            // half-tied edge would make `tiedEdges` reach `junctionEdges` on a
            // section stitched at every other station. Either side, because one
            // constraint joins *two* edges: at a butt corner the flange's edge nodes
            // are the slaves and the side plate's are the masters, and reporting only
            // the slaves would call four closed corners half-open.
            if (tieRole[edge.first] != 0u && tieRole[edge.second] != 0u)
                section.tiedEdges += span;
        }
    }
    if (section.junctionEdges - section.tiedEdges > 1e-9)
        report(std::to_string(section.junctionEdges - section.tiedEdges) + " m of " +
               std::to_string(section.junctionEdges) +
               " m of free element edge lies within " + std::to_string(params.junctionTolerance) +
               " m of another surface's plating without being joined to it: those junctions"
               " carry no shear");

    // --- The members ---------------------------------------------------------------
    //
    // The same construction `zone.cpp` uses under `Stiffeners::Modelled`: the member
    // runs along a panel seam, so the mid-surface grid points on that seam lie on its
    // own segment exactly, and the fibres are tied to the plating by `constraint.hpp`
    // rather than meshed.
    if (params.members) {
        section.stiffening.material = section.material;
        for (const StructuralMember& member : structure.members) {
            const Vec3 edge = member.b - member.a;
            const double lengthSquared = length2(edge);
            if (!(lengthSquared > 0)) continue;
            const double lo = std::min(member.a.x, member.b.x), hi = std::max(member.a.x, member.b.x);
            // A member with no extent along x -- a frame, a deck beam -- lies *on* a
            // station, and a station is a cut plane for the section on either side of
            // it. Taken at both, the chain of sections a ship is built from carries
            // that frame's steel and its stiffness **twice**, once from each
            // neighbour: measured on the reference ferry, two four-bay sections come
            // to 10.1% more stiffener mass than the same length in one piece, and the
            // doubled tie lands on shared interface DOF where it is invisible in
            // anything but a total. So it is half-open at the aft plane, which is the
            // rule the transverse *plating* above already uses and for the same
            // reason. The forward-most section of a ship loses the frame on its
            // forward plane, exactly as it loses a bulkhead there.
            if (hi - lo <= eps) {
                if (!(lo >= params.xFrom - eps && lo < params.xTo - eps)) {
                    if (lo >= params.xTo - eps && lo <= params.xTo + eps) ++section.membersOnForwardPlane;
                    continue;
                }
            } else if (hi < params.xFrom - eps || lo > params.xTo + eps) {
                continue;
            }

            // A member with no extent along x -- a frame, a deck beam, a bulkhead
            // stiffener -- carries no longitudinal stress, so leaving it out costs a
            // hull girder nothing; one with extent costs it its own area over the
            // fraction of the section it spans, which is what an average transverse
            // cut of the section would see.
            const double spanX = std::min(hi, params.xTo) - std::max(lo, params.xFrom);
            const double share = spanX > 0 ? profileSection(member.profile).area * spanX /
                                                 (params.xTo - params.xFrom)
                                           : 0.0;
            const auto lost = [&]() { section.missedMemberArea += share; };

            std::vector<std::pair<double, std::uint32_t>> onMember;
            for (std::size_t i = 0; i < mid.size(); ++i) {
                const double t = std::clamp(dot(mid[i] - member.a, edge) / lengthSquared, 0.0, 1.0);
                if (length(mid[i] - (member.a + edge * t)) > params.weldTolerance) continue;
                onMember.push_back({t, static_cast<std::uint32_t>(i)});
            }
            if (onMember.size() < 2) {
                ++section.membersMissed;
                lost();
                continue;
            }
            std::sort(onMember.begin(), onMember.end());

            // The tie measures eccentricity along the plating's own thickness
            // direction, so a web that is not perpendicular to the plate has no single
            // offset and is refused rather than projected -- a projected one would be
            // a plausible wrong second moment.
            const Vec3 rise = normalize(member.rise);
            double alignment = 0;
            for (const auto& [t, node] : onMember) alignment += dot(rise, nodeNormal[node]);
            alignment /= static_cast<double>(onMember.size());
            if (std::abs(alignment) < 0.9) {
                ++section.membersRefused;
                lost();
                continue;
            }

            // Runs of stations actually joined by an element edge **and carrying one
            // plate thickness**. Both breaks matter and the second is not obvious:
            // `constraint::addStiffener` takes a single `plateThickness` for the whole
            // run and turns it into one tie weight per fibre, `(e + t/2) / t`, while
            // the tie is applied to a node pair whose real separation is the *local*
            // nodal thickness. Feed it the run's mean and every fibre lands at
            // `e * t_local / t_mean` instead of at `e` -- measured on the reference
            // ferry at **47 mm** on members that cross a strake seam, which is a
            // quarter of a 700 mm frame's Steiner term put in the wrong place. So a
            // run stops at a thickness change, exactly as `zone.hpp` §2 stops a patch
            // at one, and each run carries the thickness its own nodes have.
            constraint::SeamRun run;
            run.sign = alignment > 0 ? 1.0 : -1.0;
            std::uint32_t previous = 0;
            double runThickness = 0;
            bool have = false, contributed = false;
            const auto flush = [&]() {
                if (run.bottom.size() >= 2 &&
                    constraint::addStiffener(run, member.profile, runThickness,
                                             section.mesh.position, section.stiffening) > 0)
                    contributed = true;
                run.bottom.clear();
                run.top.clear();
            };
            for (const auto& [t, node] : onMember) {
                const bool broken =
                    have && edgeUse.find({std::min(previous, node), std::max(previous, node)}) ==
                                edgeUse.end();
                const bool stepped = have && std::abs(nodeThickness[node] - runThickness) > 1e-12;
                if (broken || stepped) {
                    flush();
                    if (stepped && !broken) ++section.memberRunsSplitByThickness;
                }
                if (run.bottom.empty()) runThickness = nodeThickness[node];
                run.bottom.push_back(static_cast<std::uint32_t>(rank[node]) * 2);
                run.top.push_back(static_cast<std::uint32_t>(rank[node]) * 2 + 1);
                previous = node;
                have = true;
            }
            flush();
            if (contributed) {
                ++section.membersAttached;
                ++section.stiffening.members;
                section.attachedMemberArea += share;
            } else {
                ++section.membersMissed;
                lost();
            }
        }
        if (section.membersRefused > 0)
            report(std::to_string(section.membersRefused) +
                   " members have a web that is not along the plating's thickness direction and"
                   " were left out: they have no single eccentricity to tie at");
        if (section.membersMissed > 0)
            report(std::to_string(section.membersMissed) +
                   " members lie on no run of two or more mesh nodes and were left out. A member"
                   " is only picked up where it runs along a panel seam, so a girder positioned"
                   " off the longitudinal spacing is invisible to a mesh built from panels");
        if (section.missedMemberArea > 0)
            report(std::to_string(section.missedMemberArea) +
                   " m^2 of longitudinally effective member area was left out, which is what the"
                   " section will be short by against hullGirderSection");

        if (!section.stiffening.empty()) {
            const constraint::RestFibers forms =
                constraint::restFibers(section.stiffening, section.mesh.position);
            if (!forms.ok) report("a stiffener fibre came out with no length");
            // Both halves, from `attachedForms` rather than `stiffnessBlocks`
            // alone: a block says what a member adds to the stiffness and cannot
            // be asked what the member is *carrying*, so a section built without
            // the stress forms is one whose longitudinals are invisible to
            // `reduction::checkValidity` -- and the plating is the softer half, so
            // the promotion trigger then reads low. See `reduction.hpp` §9.
            constraint::AttachedForms built =
                constraint::attachedForms(section.stiffening, section.mesh.position, forms,
                                          section.material.youngsModulus);
            section.attachment.stiffness = std::move(built.stiffness);
            section.attachment.stress = std::move(built.stress);
            section.attachment.mass.assign(section.mesh.nodeCount(), 0.0);
            constraint::lumpFiberMass(section.stiffening, forms, section.material.density,
                                      section.attachment.mass);
            section.memberMass = section.stiffening.mass;
        }
    }

    return section;
}

// --- Loading a section like a beam ----------------------------------------------

namespace {

// K u, and the strain energy that goes with it, over the plating and whatever the
// `Attachment` adds. Formed from the same `elementStiffness` and the same
// `DofBlock`s `solveStatic` assembles, so the reaction reported here is the
// reaction of the system that was solved and not of a second one.
struct Applied {
    std::vector<double> force;  // K u, per global DOF
    double energy = 0;          // 0.5 u^T K u
};

Applied stiffnessTimes(const Section& section, const StructuralMaterial& material,
                       const std::vector<double>& u) {
    Applied out;
    out.force.assign(section.mesh.nodeCount() * 3, 0.0);
    double stiffness[kDof * kDof];
    double local[kDof], result[kDof];
    for (std::size_t e = 0; e < section.mesh.elementCount(); ++e) {
        double nodes[kDof];
        section.mesh.gather(e, section.mesh.position, nodes);
        solidshell::elementStiffness(nodes, material, solidshell::Formulation::SolidShell, stiffness);
        section.mesh.gather(e, u, local);
        for (int a = 0; a < kDof; ++a) {
            double sum = 0;
            for (int b = 0; b < kDof; ++b)
                sum += stiffness[static_cast<std::size_t>(a) * kDof + static_cast<std::size_t>(b)] *
                       local[b];
            result[a] = sum;
        }
        for (int a = 0; a < kNodes; ++a) {
            const std::uint32_t node = section.mesh.index[e * kNodes + static_cast<std::size_t>(a)];
            for (int k = 0; k < 3; ++k)
                out.force[static_cast<std::size_t>(node) * 3 + static_cast<std::size_t>(k)] +=
                    result[a * 3 + k];
        }
        for (int a = 0; a < kDof; ++a) out.energy += 0.5 * local[a] * result[a];
    }
    for (const solidshell::DofBlock& block : section.attachment.stiffness) {
        const std::size_t m = block.dof.size();
        for (std::size_t i = 0; i < m; ++i) {
            double sum = 0;
            for (std::size_t j = 0; j < m; ++j) sum += block.stiffness[i * m + j] * u[block.dof[j]];
            out.force[block.dof[i]] += sum;
            out.energy += 0.5 * u[block.dof[i]] * sum;
        }
    }
    // A junction tie eliminates its slave, so the force `K u` lands on a degree of
    // freedom the solved system does not have. It belongs to the masters, by the
    // transpose of the constraint -- exactly as `constraint::applyTiedForce`
    // distributes a fibre's force -- and until it is moved there the reaction at the
    // forward plane is short by whatever the junction next to it carries. **The
    // energy needs no such correction**: `u` already satisfies the constraint, so
    // `0.5 u^T K u` is `0.5 u_r^T T^T K T u_r` term for term.
    // No master is ever a slave -- `solidshell::DofExpansion` refuses a chain -- so
    // one pass distributes everything and nothing cascades.
    for (const solidshell::Mpc& mpc : section.attachment.constrained) {
        const double carried = out.force[mpc.slave];
        for (std::size_t a = 0; a < mpc.master.size(); ++a)
            out.force[mpc.master[a]] += mpc.weight[a] * carried;
        out.force[mpc.slave] = 0.0;
    }
    return out;
}

}  // namespace

BeamResponse applyBeamLoad(const Section& section, const StructuralMaterial& material,
                           const BeamLoad& load) {
    BeamResponse out;
    if (section.empty() || section.aftNodes.empty() || section.forwardNodes.empty()) {
        out.problem = "the section has no elements or no interface";
        return out;
    }

    solidshell::HexMesh mesh = section.mesh;
    const double mid = 0.5 * (section.xFrom + section.xTo);
    const auto at = [&](std::uint32_t node, int axis) {
        return mesh.position[static_cast<std::size_t>(node) * 3 + static_cast<std::size_t>(axis)];
    };

    // Only `u_x` is prescribed on the two planes. Leaving `u_y` and `u_z` free is
    // what makes this a section property rather than a section property with the
    // Poisson contraction of both ends held out of it.
    for (const std::vector<std::uint32_t>* plane : {&section.aftNodes, &section.forwardNodes})
        for (std::uint32_t node : *plane) {
            const double x = at(node, 0) - mid, z = at(node, 2);
            mesh.pin(node, 0, load.strain * x + load.curvature * x * (z - load.reference));
        }

    // Three more **per connected component**, to take out the rigid translations in
    // y and z and the rotation about x that prescribing `u_x` alone leaves. Per
    // component, not per section: a section of a real ship comes apart into pieces
    // (`section.hpp` §1) and each piece keeps its own three. Restraining only one of
    // them leaves the rest as mechanisms and the factorisation fails -- which is at
    // least loud, and was how this was found.
    //
    // They are a *statically determinate* restraint, so their reaction is exactly
    // zero on the exact solution, and `restraintReaction` says so rather than
    // assuming it.
    if (section.floatingComponents > 0) {
        out.problem = "the section has " + std::to_string(section.floatingComponents) +
                      " components touching neither cut plane: nothing prescribes their motion"
                      " and the answer would be a mechanism solved as though it were not";
        return out;
    }
    std::vector<std::uint32_t> restrained;
    {
        // A tied node has no unknown of its own, so pinning one says nothing and
        // `solveStatic` refuses it outright rather than letting the pin and the tie
        // both claim the same degree of freedom. The rigid-body restraints therefore
        // go on nodes the solve still owns.
        std::vector<std::uint8_t> tied(mesh.nodeCount(), 0u);
        for (const solidshell::Mpc& mpc : section.attachment.constrained)
            tied[mpc.slave / 3] = 1u;

        const int components = std::max(section.components, 1);
        std::vector<int> anchorOf(static_cast<std::size_t>(components), -1);
        std::vector<int> leverOf(static_cast<std::size_t>(components), -1);
        std::vector<double> armOf(static_cast<std::size_t>(components), 0.0);
        for (std::uint32_t node = 0; node < mesh.nodeCount(); ++node) {
            const int component = section.componentOf[node];
            if (component < 0 || tied[node]) continue;
            const auto slot = static_cast<std::size_t>(component);
            if (anchorOf[slot] < 0) {
                anchorOf[slot] = static_cast<int>(node);
                continue;
            }
            const auto anchor = static_cast<std::uint32_t>(anchorOf[slot]);
            const double dy = at(node, 1) - at(anchor, 1), dz = at(node, 2) - at(anchor, 2);
            const double arm = std::sqrt(dy * dy + dz * dz);
            if (arm > armOf[slot]) {
                armOf[slot] = arm;
                leverOf[slot] = static_cast<int>(node);
            }
        }
        for (int component = 0; component < components; ++component) {
            const auto slot = static_cast<std::size_t>(component);
            if (anchorOf[slot] < 0 || leverOf[slot] < 0 || !(armOf[slot] > 0)) {
                out.problem = "a connected component has no lever arm in the transverse plane, so"
                              " its rotation about x cannot be held";
                return out;
            }
            const auto anchor = static_cast<std::uint32_t>(anchorOf[slot]);
            const auto lever = static_cast<std::uint32_t>(leverOf[slot]);
            mesh.pin(anchor, 1, 0.0);
            mesh.pin(anchor, 2, 0.0);
            restrained.push_back(anchor * 3 + 1);
            restrained.push_back(anchor * 3 + 2);
            // Hold whichever of the lever's two transverse degrees of freedom the arm
            // is *perpendicular* to. A plate lying in the x-z plane has every node at
            // the same y to within its own thickness, so pinning `u_z` there would
            // hold the rotation about x through a 12 mm arm and the restraint would be
            // a conditioning problem rather than a choice.
            const double dy = at(lever, 1) - at(anchor, 1), dz = at(lever, 2) - at(anchor, 2);
            const int axis = std::abs(dy) >= std::abs(dz) ? 2 : 1;
            mesh.pin(lever, axis, 0.0);
            restrained.push_back(lever * 3 + static_cast<std::uint32_t>(axis));
        }
    }

    std::vector<double> displacement;
    const std::vector<double> zeroLoad(mesh.nodeCount() * 3, 0.0);
    std::string problem;
    if (!solidshell::solveStatic(mesh, material, solidshell::Formulation::SolidShell,
                                 section.attachment.stiffness, section.attachment.constrained,
                                 zeroLoad, displacement, &problem)) {
        out.problem = problem;
        return out;
    }

    const Applied applied = stiffnessTimes(section, material, displacement);
    for (std::uint32_t node : section.forwardNodes) {
        const double reaction = applied.force[static_cast<std::size_t>(node) * 3];
        out.axialForce += reaction;
        out.bendingMoment += reaction * (at(node, 2) - load.reference);
    }
    out.strainEnergy = applied.energy;
    for (std::size_t d = 0; d < displacement.size(); ++d) {
        if (mesh.fixed[d]) continue;
        out.residual = std::max(out.residual, std::abs(applied.force[d]));
    }
    // The rigid-body restraints carry **exactly zero** on the exact solution, because
    // they pick one of a family of zero-energy motions rather than holding anything.
    // A non-zero reading here is a mesher or a solver defect and not a modelling
    // choice, which is what makes it worth reporting separately from `residual`.
    for (std::uint32_t d : restrained)
        out.restraintReaction = std::max(out.restraintReaction, std::abs(applied.force[d]));
    for (std::size_t node = 0; node < mesh.nodeCount(); ++node)
        out.peakDisplacement = std::max(
            out.peakDisplacement,
            length(Vec3{displacement[node * 3], displacement[node * 3 + 1], displacement[node * 3 + 2]}));
    if (load.strain != 0) out.axialStiffness = out.axialForce / load.strain;
    if (load.curvature != 0) out.bendingStiffness = out.bendingMoment / load.curvature;
    out.ok = true;
    return out;
}

TorsionResponse applyTwist(const Section& section, const StructuralMaterial& material, double twist,
                           double reference) {
    TorsionResponse out;
    if (section.empty() || section.aftNodes.empty() || section.forwardNodes.empty()) {
        out.problem = "the section has no elements or no interface";
        return out;
    }
    solidshell::HexMesh mesh = section.mesh;
    const auto at = [&](std::uint32_t node, int axis) {
        return mesh.position[static_cast<std::size_t>(node) * 3 + static_cast<std::size_t>(axis)];
    };
    for (std::uint32_t node : section.aftNodes)
        for (int k = 0; k < 3; ++k) mesh.pin(node, k, 0.0);
    for (std::uint32_t node : section.forwardNodes) {
        const double y = at(node, 1), z = at(node, 2) - reference;
        mesh.pin(node, 0, 0.0);
        mesh.pin(node, 1, -twist * z);
        mesh.pin(node, 2, twist * y);
    }

    std::vector<double> displacement;
    const std::vector<double> zeroLoad(mesh.nodeCount() * 3, 0.0);
    std::string problem;
    if (!solidshell::solveStatic(mesh, material, solidshell::Formulation::SolidShell,
                                 section.attachment.stiffness, section.attachment.constrained,
                                 zeroLoad, displacement, &problem)) {
        out.problem = problem;
        return out;
    }
    const Applied applied = stiffnessTimes(section, material, displacement);
    for (std::uint32_t node : section.forwardNodes) {
        const double y = at(node, 1), z = at(node, 2) - reference;
        const double fy = applied.force[static_cast<std::size_t>(node) * 3 + 1];
        const double fz = applied.force[static_cast<std::size_t>(node) * 3 + 2];
        out.torque += y * fz - z * fy;
    }
    out.strainEnergy = applied.energy;
    if (twist != 0) out.torsionalStiffness = out.torque * section.length() / twist;
    out.ok = true;
    return out;
}

// --- 6. A ship: a chain of sections ----------------------------------------------

namespace {

double nowSeconds() {
    using namespace std::chrono;
    return duration<double>(steady_clock::now().time_since_epoch()).count();
}

// The boundary indices of a substructure whose node is one of `nodes`, ascending.
std::vector<std::size_t> boundaryOnPlane(const reduction::Substructure& substructure,
                                         const std::vector<std::uint32_t>& nodes) {
    std::vector<std::size_t> out;
    if (nodes.empty()) return out;
    std::uint32_t highest = 0;
    for (std::uint32_t n : nodes) highest = std::max(highest, n);
    std::vector<std::uint8_t> onPlane(static_cast<std::size_t>(highest) + 1, 0u);
    for (std::uint32_t n : nodes) onPlane[n] = 1u;
    const std::vector<std::uint32_t>& dof = substructure.boundaryDof();
    for (std::size_t b = 0; b < dof.size(); ++b) {
        const std::size_t node = dof[b] / 3u;
        if (node < onPlane.size() && onPlane[node]) out.push_back(b);
    }
    return out;
}

}  // namespace

bool Chain::ready() const {
    // **Not `components == 1`.** A chain of sections whose plating is not tied is
    // several pieces and is still a model worth solving -- it is the negative control
    // every junction claim is made against, and `applyBeamLoad` restrains each piece.
    // What makes a chain unusable is an interface that did not match, because that is
    // two unknowns where the ship has one.
    if (section.empty() || assembly.empty()) return false;
    if (reduced.size() != section.size() || substructure.size() != section.size()) return false;
    for (const reduction::Substructure& s : substructure)
        if (!s.ready()) return false;
    for (const reduction::Reduction& r : reduced)
        if (r.empty()) return false;
    for (std::size_t u : unmatched)
        if (u != 0) return false;
    return !aftDof.empty() && !forwardDof.empty();
}

Chain buildChain(const StructuralMesh& structure, const ChainParams& params) {
    Chain chain;
    const auto report = [&](std::string message) { chain.problems.push_back(std::move(message)); };
    if (params.station.size() < 2) {
        report("a chain needs at least two cut planes");
        return chain;
    }
    for (std::size_t i = 0; i + 1 < params.station.size(); ++i)
        if (!(params.station[i + 1] > params.station[i])) {
            report("the cut planes are not ascending");
            return chain;
        }
    chain.xFrom = params.station.front();
    chain.xTo = params.station.back();

    // Every section first, and the vector sized once: a `reduction::Substructure`
    // holds a reference to the mesh it was built from, so a reallocation here would
    // leave every substructure pointing at a freed mesh.
    const std::size_t pieces = params.station.size() - 1;
    chain.section.reserve(pieces);
    double t = nowSeconds();
    for (std::size_t i = 0; i < pieces; ++i) {
        SectionParams p = params.section;
        p.xFrom = params.station[i];
        p.xTo = params.station[i + 1];
        chain.section.push_back(buildSection(structure, p));
        const Section& s = chain.section.back();
        chain.tiedEdges += s.tiedEdges;
        chain.junctionEdges += s.junctionEdges;
        if (s.empty()) {
            report("section " + std::to_string(i) + " meshed no elements");
            return chain;
        }
        if (s.floatingComponents > 0)
            report("section " + std::to_string(i) + " has " +
                   std::to_string(s.floatingComponents) +
                   " components touching neither cut plane, which nothing in the assembled model"
                   " restrains");
    }
    chain.meshSeconds = nowSeconds() - t;

    t = nowSeconds();
    chain.substructure.reserve(pieces);
    chain.reduced.reserve(pieces);
    for (std::size_t i = 0; i < pieces; ++i) {
        const Section& s = chain.section[i];
        chain.substructure.emplace_back(s.mesh, s.material, s.interfaceNodes, s.attachment);
        if (!chain.substructure.back().ready()) {
            report("section " + std::to_string(i) + " did not reduce: " +
                   (chain.substructure.back().problems().empty()
                        ? std::string("no reason given")
                        : chain.substructure.back().problems().front()));
            return chain;
        }
    }
    for (std::size_t i = 0; i < pieces; ++i) {
        chain.reduced.push_back(reduction::craigBampton(chain.substructure[i], params.reduce));
        if (chain.reduced.back().empty()) {
            report("section " + std::to_string(i) + " reduced to nothing");
            return chain;
        }
    }
    chain.reduceSeconds = nowSeconds() - t;

    // Neighbours only. Two sections two bays apart share no cut plane, and asking
    // `matchComponents` for every pair would cost k^2 boundary matches to find that
    // out -- 1 170 x 1 170 distance tests each on this ship.
    t = nowSeconds();
    std::vector<reduction::Component> parts(pieces);
    for (std::size_t i = 0; i < pieces; ++i)
        parts[i] = {&chain.substructure[i], &chain.reduced[i]};
    const std::vector<reduction::Joint> joints =
        reduction::matchNeighbours(parts, params.matchTolerance);
    if (pieces > 1 && joints.size() != pieces - 1)
        report("only " + std::to_string(joints.size()) + " of " + std::to_string(pieces - 1) +
               " interior cut planes joined anything at all: consecutive sections that share no"
               " boundary DOF assemble into independent models standing side by side");

    // What each interior plane actually matched. The count that matters is the one
    // that did *not*: a chain assembles and solves with unmatched interface DOF, and
    // what it has built is a ship torn along part of a cut.
    chain.shared.assign(pieces > 0 ? pieces - 1 : 0, 0);
    chain.unmatched.assign(pieces > 0 ? pieces - 1 : 0, 0);
    for (const reduction::Joint& joint : joints) {
        if (joint.b != joint.a + 1) continue;
        const auto plane = static_cast<std::size_t>(joint.a);
        if (plane >= chain.shared.size()) continue;
        chain.shared[plane] = joint.map.shared;
        chain.worstGap = std::max(chain.worstGap, joint.map.worstGap);
    }
    // A boundary DOF of section i that lies on the shared plane and found no partner
    // in section i+1.
    //
    // **Which DOF are on the plane comes from the mesher, not from a position test.**
    // A mesh node is the mid-surface offset by `t/2` along the plating normal, so a
    // node of a panel whose normal leans out of the transverse plane is up to `t/2`
    // off it -- 8 mm on the reference ferry, where `SectionParams::planeTolerance` is
    // 1e-6. `Section::aftNodes` and `forwardNodes` are chosen on the mid-surface and
    // take both of each pair, which is exact; asking `x == plane` instead silently
    // drops every strip of plating that is not square to the cut, and a chain whose
    // end planes were selected that way came out 33% short in `EA` on a shoulder
    // while reporting nothing wrong at all.
    for (std::size_t i = 0; i + 1 < pieces; ++i) {
        std::vector<int> matched(chain.substructure[i].boundaryDof().size(), -1);
        for (const reduction::Joint& joint : joints)
            if (joint.a == static_cast<int>(i) && joint.b == static_cast<int>(i) + 1)
                for (std::size_t b = 0; b < joint.map.aToB.size() && b < matched.size(); ++b)
                    matched[b] = joint.map.aToB[b];
        for (std::size_t b : boundaryOnPlane(chain.substructure[i], chain.section[i].forwardNodes))
            if (matched[b] < 0) ++chain.unmatched[i];
    }
    chain.assembly = reduction::assemble(parts, joints);
    chain.assembleSeconds = nowSeconds() - t;
    for (const std::string& p : chain.assembly.problems) report(p);
    if (chain.assembly.empty()) {
        report("the chain assembled into nothing");
        return chain;
    }
    chain.reducedComponents = reduction::assembledComponents(chain.assembly);
    if (chain.reducedComponents != 1)
        report("the assembled reduced model is in " + std::to_string(chain.reducedComponents) +
               " pieces: consecutive sections do not share the cut plane between them");

    // Pieces of the *structure*, which is a different count and the one a ship is
    // judged by. `reduction::assembledComponents` cannot see inside a section -- a
    // component's reduced pair is dense over its own DOF whether the mesh behind it
    // is one piece or seven -- and `Section::components` cannot see across a cut
    // plane. This joins the two: a mesh component of one section and a mesh component
    // of the next are the same piece of ship when they share an assembled row.
    {
        std::vector<std::size_t> base(pieces + 1, 0);
        for (std::size_t i = 0; i < pieces; ++i)
            base[i + 1] = base[i] + static_cast<std::size_t>(std::max(chain.section[i].components, 1));
        chain.pieceOf.assign(base[pieces], -1);
        std::vector<int> parent(base[pieces]);
        for (std::size_t i = 0; i < parent.size(); ++i) parent[i] = static_cast<int>(i);
        const std::function<int(int)> find = [&](int a) {
            while (parent[static_cast<std::size_t>(a)] != a)
                a = parent[static_cast<std::size_t>(a)] =
                    parent[static_cast<std::size_t>(parent[static_cast<std::size_t>(a)])];
            return a;
        };
        std::vector<int> claim(static_cast<std::size_t>(chain.assembly.boundary), -1);
        for (std::size_t i = 0; i < pieces; ++i) {
            const std::vector<std::uint32_t>& dof = chain.substructure[i].boundaryDof();
            const std::vector<int>& from = chain.assembly.from[i];
            for (std::size_t b = 0; b < dof.size() && b < from.size(); ++b) {
                const std::size_t node = dof[b] / 3u;
                if (node >= chain.section[i].componentOf.size()) continue;
                const int component = chain.section[i].componentOf[node];
                if (component < 0 || from[b] < 0) continue;
                const int key = static_cast<int>(base[i]) + component;
                const auto row = static_cast<std::size_t>(from[b]);
                if (claim[row] < 0) {
                    claim[row] = key;
                    continue;
                }
                const int x = find(claim[row]), y = find(key);
                if (x != y) parent[static_cast<std::size_t>(x)] = y;
            }
        }
        std::map<int, int> label;
        for (std::size_t d = 0; d < chain.pieceOf.size(); ++d) {
            const int root = find(static_cast<int>(d));
            auto found = label.find(root);
            if (found == label.end()) found = label.emplace(root, static_cast<int>(label.size())).first;
            chain.pieceOf[d] = found->second;
        }
        chain.componentBase = base;
        chain.components = static_cast<int>(label.size());
    }
    if (chain.components != 1)
        report("the chain is in " + std::to_string(chain.components) +
               " structurally disconnected pieces, counting the mesh inside each section and the"
               " cut planes between them");
    for (std::size_t i = 0; i < chain.unmatched.size(); ++i)
        if (chain.unmatched[i] > 0)
            report(std::to_string(chain.unmatched[i]) +
                   " boundary DOF on the cut plane at x = " + std::to_string(params.station[i + 1]) +
                   " found no partner in the next section and are two unknowns where the ship has"
                   " one: raise ChainParams::matchTolerance, and see section.hpp section 6 note 1"
                   " for why they are not coincident");

    // The two end planes, in assembled numbering: the aft-most section's own aft plane
    // and the forward-most section's own forward plane, mapped through `from`. See
    // `boundaryOnPlane` for why this is not a position test.
    if (chain.assembly.boundaryPoint.size() != static_cast<std::size_t>(chain.assembly.boundary))
        report("the assembly carries no boundary DOF identity, so nothing can be prescribed on it"
               " by position");
    {
        const auto collect = [&](std::size_t i, const std::vector<std::uint32_t>& nodes,
                                 std::vector<std::uint32_t>& into) {
            const std::vector<int>& from = chain.assembly.from[i];
            for (std::size_t b : boundaryOnPlane(chain.substructure[i], nodes))
                if (b < from.size() && from[b] >= 0)
                    into.push_back(static_cast<std::uint32_t>(from[b]));
        };
        collect(0, chain.section.front().aftNodes, chain.aftDof);
        collect(pieces - 1, chain.section.back().forwardNodes, chain.forwardDof);
    }
    if (chain.aftDof.empty() || chain.forwardDof.empty())
        report("one of the chain's end planes carries no assembled boundary DOF");
    return chain;
}

namespace {

// The reaction at every assembled DOF: `K x` over the dense assembled stiffness.
// This is the one thing that would notice a stiffness assembled one way and a
// reaction taken another, which is why it is formed from `assembly.stiffness` rather
// than from the sections' element matrices.
std::vector<double> assembledReaction(const reduction::Assembly& assembly,
                                      const std::vector<double>& state) {
    const std::size_t n = static_cast<std::size_t>(assembly.size());
    std::vector<double> out(n, 0.0);
    if (state.size() != n || assembly.stiffness.size() != n * n) return out;
    for (std::size_t i = 0; i < n; ++i) {
        const double* row = &assembly.stiffness[i * n];
        double sum = 0;
        for (std::size_t j = 0; j < n; ++j) sum += row[j] * state[j];
        out[i] = sum;
    }
    return out;
}

// Which assembled boundary rows belong to each structural piece of the chain.
std::vector<std::vector<std::uint32_t>> rowsPerPiece(const Chain& chain) {
    std::vector<std::vector<std::uint32_t>> out(static_cast<std::size_t>(std::max(chain.components, 0)));
    if (out.empty() || chain.componentBase.size() != chain.section.size() + 1) return {};
    std::vector<std::uint8_t> taken(static_cast<std::size_t>(chain.assembly.boundary), 0u);
    for (std::size_t i = 0; i < chain.section.size(); ++i) {
        const std::vector<std::uint32_t>& dof = chain.substructure[i].boundaryDof();
        const std::vector<int>& from = chain.assembly.from[i];
        for (std::size_t b = 0; b < dof.size() && b < from.size(); ++b) {
            const std::size_t node = dof[b] / 3u;
            if (node >= chain.section[i].componentOf.size() || from[b] < 0) continue;
            const int component = chain.section[i].componentOf[node];
            if (component < 0) continue;
            const auto key = chain.componentBase[i] + static_cast<std::size_t>(component);
            if (key >= chain.pieceOf.size()) continue;
            const int piece = chain.pieceOf[key];
            const auto row = static_cast<std::size_t>(from[b]);
            if (piece < 0 || static_cast<std::size_t>(piece) >= out.size() || taken[row]) continue;
            taken[row] = 1u;
            out[static_cast<std::size_t>(piece)].push_back(static_cast<std::uint32_t>(row));
        }
    }
    return out;
}

// Three restraints **per structural piece**, which is what takes out the rigid body
// motions prescribing `u_x` alone leaves: translation in y and z, and rotation about
// x. Per piece and not per chain, for exactly the reason `applyBeamLoad` gives for
// one section -- a section of a real ship comes apart into pieces, each keeps its
// own three, and restraining only one leaves the rest as mechanisms and the
// factorisation fails. An untied chain of the reference ferry is seven pieces.
//
// They are a *statically determinate* restraint, so their reaction is exactly zero on
// the exact solution, and `restraintReaction` says so rather than assuming it.
//
// Chosen off the assembly's own boundary identity rather than off a mesh: an anchor
// with both transverse DOF in the model, and whichever row of the piece has the
// longest arm from it in the transverse plane, held on the axis the arm is
// **perpendicular** to -- a deck node's arm is almost all in y, so holding its `u_y`
// would carry the rotation about x through a lever of nothing.
bool chainRestraints(const Chain& chain, std::vector<std::uint32_t>& held) {
    held.clear();
    const std::vector<reduction::BoundaryDof>& point = chain.assembly.boundaryPoint;
    if (point.size() != static_cast<std::size_t>(chain.assembly.boundary)) return false;
    const std::vector<std::vector<std::uint32_t>> rows = rowsPerPiece(chain);
    if (rows.empty()) return false;

    for (const std::vector<std::uint32_t>& piece : rows) {
        int anchorY = -1, anchorZ = -1;
        Vec3 anchor{};
        for (std::uint32_t b : piece) {
            if (point[b].axis != 1) continue;
            for (std::uint32_t c : piece)
                if (point[c].axis == 2 && length(point[c].position - point[b].position) < 1e-12) {
                    anchorY = static_cast<int>(b);
                    anchorZ = static_cast<int>(c);
                    anchor = point[b].position;
                    break;
                }
            if (anchorY >= 0) break;
        }
        if (anchorY < 0 || anchorZ < 0) return false;

        int lever = -1;
        double best = 0;
        for (std::uint32_t b : piece) {
            const reduction::BoundaryDof& d = point[b];
            const double dy = d.position.y - anchor.y, dz = d.position.z - anchor.z;
            const double arm = std::sqrt(dy * dy + dz * dz);
            const std::uint32_t axis = std::fabs(dy) >= std::fabs(dz) ? 2u : 1u;
            if (d.axis != axis || !(arm > best)) continue;
            best = arm;
            lever = static_cast<int>(b);
        }
        if (lever < 0 || !(best > 0)) return false;
        held.push_back(static_cast<std::uint32_t>(anchorY));
        held.push_back(static_cast<std::uint32_t>(anchorZ));
        held.push_back(static_cast<std::uint32_t>(lever));
    }
    return true;
}

}  // namespace

BeamResponse applyBeamLoad(const Chain& chain, const BeamLoad& load) {
    BeamResponse out;
    if (!chain.ready()) {
        out.problem = "the chain is not ready: " +
                      (chain.problems.empty() ? std::string("no sections") : chain.problems.front());
        return out;
    }
    const std::vector<reduction::BoundaryDof>& point = chain.assembly.boundaryPoint;
    const double mid = 0.5 * (chain.xFrom + chain.xTo);

    // Only `u_x` on the two end planes, exactly as on one section: pinning `u_y` and
    // `u_z` there would suppress the Poisson contraction and report a ship stiffer
    // than it is by roughly `1 / (1 - nu^2)`.
    std::vector<reduction::Prescribed> prescribed;
    for (const std::vector<std::uint32_t>* plane : {&chain.aftDof, &chain.forwardDof})
        for (std::uint32_t b : *plane) {
            const reduction::BoundaryDof& d = point[b];
            if (d.axis != 0) continue;
            const double x = d.position.x - mid, z = d.position.z;
            prescribed.push_back({b, load.strain * x + load.curvature * x * (z - load.reference)});
        }
    std::vector<std::uint32_t> held;
    if (!chainRestraints(chain, held)) {
        out.problem = "the chain's end planes carry no lever arm in the transverse plane, so the"
                      " rotation about x cannot be held";
        return out;
    }

    std::vector<double> state;
    const std::vector<double> zeroLoad(static_cast<std::size_t>(chain.assembly.size()), 0.0);
    std::string problem;
    if (!reduction::assembledStaticSolve(chain.assembly, zeroLoad, held, prescribed, state,
                                         &problem)) {
        out.problem = problem;
        return out;
    }

    const std::vector<double> reaction = assembledReaction(chain.assembly, state);
    for (std::uint32_t b : chain.forwardDof) {
        const reduction::BoundaryDof& d = point[b];
        if (d.axis != 0) continue;
        out.axialForce += reaction[b];
        out.bendingMoment += reaction[b] * (d.position.z - load.reference);
    }
    for (std::size_t i = 0; i < state.size(); ++i) out.strainEnergy += 0.5 * state[i] * reaction[i];

    std::vector<std::uint8_t> fixed(state.size(), 0u);
    for (const reduction::Prescribed& p : prescribed) fixed[p.dof] = 1u;
    for (std::uint32_t d : held) fixed[d] = 1u;
    for (std::size_t i = 0; i < state.size(); ++i)
        if (!fixed[i]) out.residual = std::max(out.residual, std::fabs(reaction[i]));
    for (std::uint32_t d : held)
        out.restraintReaction = std::max(out.restraintReaction, std::fabs(reaction[d]));
    // Over the boundary rows only. A modal coordinate is not a displacement -- it is
    // an amplitude against a mass-normalised mode shape -- so taking the largest of
    // the whole reduced state would report metres that are not metres.
    for (int b = 0; b < chain.assembly.boundary; ++b)
        out.peakDisplacement =
            std::max(out.peakDisplacement, std::fabs(state[static_cast<std::size_t>(b)]));
    if (load.strain != 0) out.axialStiffness = out.axialForce / load.strain;
    if (load.curvature != 0) out.bendingStiffness = out.bendingMoment / load.curvature;
    out.ok = true;
    return out;
}

TorsionResponse applyTwist(const Chain& chain, double twist, double reference) {
    TorsionResponse out;
    if (!chain.ready()) {
        out.problem = "the chain is not ready: " +
                      (chain.problems.empty() ? std::string("no sections") : chain.problems.front());
        return out;
    }
    const std::vector<reduction::BoundaryDof>& point = chain.assembly.boundaryPoint;

    // A rigid disc at each end: every degree of freedom prescribed, the aft plane
    // held and the forward one turned. A twist is a statement about the section's
    // shape and leaving anything free there would let the ship shear instead of turn.
    std::vector<reduction::Prescribed> prescribed;
    for (std::uint32_t b : chain.aftDof) prescribed.push_back({b, 0.0});
    for (std::uint32_t b : chain.forwardDof) {
        const reduction::BoundaryDof& d = point[b];
        const double y = d.position.y, z = d.position.z - reference;
        prescribed.push_back({b, d.axis == 0 ? 0.0 : (d.axis == 1 ? -twist * z : twist * y)});
    }

    std::vector<double> state;
    const std::vector<double> zeroLoad(static_cast<std::size_t>(chain.assembly.size()), 0.0);
    std::string problem;
    if (!reduction::assembledStaticSolve(chain.assembly, zeroLoad, {}, prescribed, state,
                                         &problem)) {
        out.problem = problem;
        return out;
    }
    const std::vector<double> reaction = assembledReaction(chain.assembly, state);
    // The torque is the moment of the forward plane's reaction about x, and it needs
    // both transverse components of the same node -- so it is summed per DOF with the
    // arm each one has, which is the same sum written the other way round.
    for (std::uint32_t b : chain.forwardDof) {
        const reduction::BoundaryDof& d = point[b];
        if (d.axis == 1) out.torque += -(d.position.z - reference) * reaction[b];
        else if (d.axis == 2) out.torque += d.position.y * reaction[b];
    }
    for (std::size_t i = 0; i < state.size(); ++i) out.strainEnergy += 0.5 * state[i] * reaction[i];
    if (twist != 0) out.torsionalStiffness = out.torque * chain.length() / twist;
    out.ok = true;
    return out;
}

std::vector<double> chainFrequencies(const Chain& chain) {
    if (chain.assembly.empty()) return {};
    std::vector<std::uint32_t> held = chain.aftDof;
    held.insert(held.end(), chain.forwardDof.begin(), chain.forwardDof.end());
    return reduction::assembledFrequencies(chain.assembly, held);
}

}  // namespace sim::section
