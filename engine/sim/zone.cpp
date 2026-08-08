// SPDX-License-Identifier: MIT
#include "zone.hpp"

#include "../core/jobs.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <map>
#include <string>
#include <utility>

namespace sim::zone {
namespace {

using solidshell::kDof;
using solidshell::kGauss;
using solidshell::kNodes;

// Per-element cost, from the measurements in `02-simulation.md` §3. Used only to
// predict; `SolveResult::microsecondsPerElementStep` is what actually happened.
// 3.1 rather than the 7.3 `02-simulation.md` §3 first published: the step-invariant
// element forms are now built once at promotion rather than every step, which is
// the whole of the difference. Measured both ways on the same element, 7.18 µs
// against 3.06.
constexpr double kPlasticMicroseconds = 3.1;
constexpr double kElasticMicroseconds = 0.273;

// --- Welding -------------------------------------------------------------------
//
// Adjacent panels from `makeStructuralMesh` share their corners *exactly* -- both
// come from one `Section::at(fraction)` call -- and the bilinear points along a
// shared edge agree to the last bit whenever the two panels number that edge the
// same way round. They do not when one of them is the mirrored starboard panel, so
// the merge is by distance rather than by equality. A bucket grid, probed as far as
// the tolerance reaches, because quantising alone splits two points that straddle a
// cell boundary -- which is the classic way a weld leaves a crack down the middle of
// a mesh, and the crack is invisible in everything except the connectivity.
class Welder {
public:
    explicit Welder(double tolerance)
        : tolerance_(tolerance),
          cell_(2.0 * tolerance),
          // The probe reaches as far as the tolerance does, whatever the cell size
          // happens to be. Hard-coding +/-1 ties the two together silently: shrink
          // the cell and the weld quietly stops being a distance, which is the
          // failure mutation testing found here -- a mutant that halved the cell
          // passed everything while leaving a crack down the middle of any mesh
          // whose duplicates sat between 0.75 and 1 tolerance apart.
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
    long long cellOf(double v) const {
        return static_cast<long long>(std::floor(v / cell_));
    }
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

// The two in-plane directions of a patch whose outward normal is `axis`. `right`
// follows the ship's length for a vertical side shell, which is the direction the
// frame spacing -- the membrane model's span -- is measured along.
void inPlaneFrame(const Vec3& axis, Vec3& right, Vec3& up) {
    Vec3 seed{0, 0, 1};
    if (std::abs(dot(seed, axis)) > 0.9) seed = Vec3{1, 0, 0};
    right = normalize(cross(seed, axis));
    up = cross(axis, right);
}

}  // namespace

// --- Meshing -------------------------------------------------------------------

Patch buildPatch(const StructuralMesh& structure, const Vec3& impact, const MeshParams& params) {
    Patch patch;
    const auto report = [&](std::string message) { patch.problems.push_back(std::move(message)); };

    if (!(params.radius > 0)) {
        report("zone radius is not positive");
        return patch;
    }
    if (params.subdivision < 1) {
        report("subdivision below one element per panel");
        return patch;
    }

    // Candidates: plating of the requested role whose centroid the zone reaches.
    // The flood fill below never leaves this set, so adjacency is only ever built
    // over it -- which is the difference between a few hundred panels and the
    // ferry's nine thousand.
    std::vector<int> candidate;
    for (std::size_t i = 0; i < structure.panels.size(); ++i) {
        const PlatePanel& p = structure.panels[i];
        if (p.role != params.role) continue;
        if (length(p.centroid() - impact) > params.radius) continue;
        candidate.push_back(static_cast<int>(i));
    }
    if (candidate.empty()) {
        report("no plating of that role within " + std::to_string(params.radius) +
               " m of the impact");
        return patch;
    }

    int struck = candidate.front();
    double nearest = length(structure.panels[static_cast<std::size_t>(struck)].centroid() - impact);
    for (int index : candidate) {
        const double d = length(structure.panels[static_cast<std::size_t>(index)].centroid() - impact);
        if (d < nearest) {
            nearest = d;
            struck = index;
        }
    }

    const PlatePanel& struckPanel = structure.panels[static_cast<std::size_t>(struck)];
    patch.struckPanel = struck;
    patch.thickness = struckPanel.thickness;
    if (!(patch.thickness > 0)) {
        report("the struck panel has no thickness");
        return patch;
    }

    // Which way is out. Away from the structure's own centroid is right for shell
    // plating and arbitrary for a bulkhead, so a caller who knows says so.
    Vec3 axis = struckPanel.normal();
    if (length2(params.outward) > 0) {
        if (dot(axis, params.outward) < 0) axis = -axis;
    } else {
        Vec3 centroid{0, 0, 0};
        double total = 0;
        for (const PlatePanel& p : structure.panels) {
            const double a = p.area();
            centroid += p.centroid() * a;
            total += a;
        }
        if (total > 0) centroid = centroid / total;
        if (dot(axis, struckPanel.centroid() - centroid) < 0) axis = -axis;
    }
    patch.axis = axis;
    inPlaneFrame(axis, patch.right, patch.up);
    patch.centre = impact - axis * dot(impact - struckPanel.centroid(), axis);

    const int materialIndex = struckPanel.material;
    patch.material = materialIndex >= 0 &&
                             materialIndex < static_cast<int>(structure.materials.size())
                         ? structure.materials[static_cast<std::size_t>(materialIndex)]
                         : ah36Steel();

    // --- Adjacency over the candidate set, by welded corner ---------------------

    std::vector<Vec3> cornerPoint;
    Welder cornerWeld(params.weldTolerance);
    std::vector<std::array<std::uint32_t, 4>> panelCorner(candidate.size());
    std::map<int, std::size_t> slotOf;
    for (std::size_t s = 0; s < candidate.size(); ++s) {
        const PlatePanel& p = structure.panels[static_cast<std::size_t>(candidate[s])];
        for (int c = 0; c < 4; ++c)
            panelCorner[s][static_cast<std::size_t>(c)] = cornerWeld.weld(p.corner[c], cornerPoint);
        slotOf[candidate[s]] = s;
    }
    std::map<std::pair<std::uint32_t, std::uint32_t>, std::vector<std::size_t>> edgeUsers;
    for (std::size_t s = 0; s < candidate.size(); ++s)
        for (int e = 0; e < 4; ++e) {
            const std::uint32_t a = panelCorner[s][static_cast<std::size_t>(e)];
            const std::uint32_t b = panelCorner[s][static_cast<std::size_t>((e + 1) % 4)];
            if (a == b) continue;  // a band the hull has closed to nothing
            edgeUsers[{std::min(a, b), std::max(a, b)}].push_back(s);
        }

    // --- Flood fill --------------------------------------------------------------
    //
    // Breadth first from the struck panel, so the patch is *connected*: a radius
    // test alone would pick up the mirrored panel on the other side of the ship
    // wherever the hull is narrower than the zone.
    std::vector<std::uint8_t> taken(candidate.size(), 0u);
    std::vector<Vec3> oriented(candidate.size());
    std::vector<std::size_t> order;
    const std::size_t seed = slotOf[struck];
    taken[seed] = 1u;
    oriented[seed] = axis;
    order.push_back(seed);
    int foldStops = 0, thicknessStops = 0, materialStops = 0, thicknessCrossings = 0;

    for (std::size_t head = 0; head < order.size(); ++head) {
        const std::size_t s = order[head];
        for (int e = 0; e < 4; ++e) {
            const std::uint32_t a = panelCorner[s][static_cast<std::size_t>(e)];
            const std::uint32_t b = panelCorner[s][static_cast<std::size_t>((e + 1) % 4)];
            if (a == b) continue;
            auto users = edgeUsers.find({std::min(a, b), std::max(a, b)});
            if (users == edgeUsers.end()) continue;
            for (std::size_t other : users->second) {
                if (other == s || taken[other]) continue;
                const PlatePanel& q = structure.panels[static_cast<std::size_t>(candidate[other])];
                if (std::abs(q.thickness - patch.thickness) > 1e-9) {
                    if (params.singleThickness) {
                        ++thicknessStops;
                        continue;
                    }
                    ++thicknessCrossings;
                }
                if (q.material != materialIndex) {
                    ++materialStops;
                    continue;
                }
                Vec3 n = q.normal();
                if (dot(n, oriented[s]) < 0) n = -n;
                if (dot(n, oriented[s]) < std::cos(params.foldLimit)) {
                    ++foldStops;
                    continue;
                }
                taken[other] = 1u;
                oriented[other] = n;
                order.push_back(other);
            }
        }
    }
    if (thicknessStops > 0)
        report("the zone stopped at a plate thickness seam on " +
               std::to_string(thicknessStops) +
               " edges: the elements either side of a seam cannot share a node without"
               " putting a taper inside one of them");
    if (thicknessCrossings > 0)
        report(std::to_string(thicknessCrossings) +
               " edges crossed a plate thickness seam and were meshed at the struck panel's " +
               std::to_string(patch.thickness) + " m anyway, because one patch has one extrusion"
               " distance");
    if (foldStops > 0)
        report(std::to_string(foldStops) + " edges folded further than the zone follows");
    if (materialStops > 0)
        report(std::to_string(materialStops) + " edges changed material and the zone stopped");

    std::sort(order.begin(), order.end());

    // --- Mid-surface grid --------------------------------------------------------

    const int n = params.subdivision;
    std::vector<Vec3> mid;
    Welder midWeld(params.weldTolerance);
    struct SubQuad {
        std::uint32_t node[4];
        int panel;
    };
    std::vector<SubQuad> quads;
    quads.reserve(order.size() * static_cast<std::size_t>(n) * static_cast<std::size_t>(n));

    std::vector<std::uint32_t> grid(static_cast<std::size_t>(n + 1) * static_cast<std::size_t>(n + 1));
    for (std::size_t s : order) {
        const int panelIndex = candidate[s];
        const PlatePanel& p = structure.panels[static_cast<std::size_t>(panelIndex)];
        for (int i = 0; i <= n; ++i)
            for (int j = 0; j <= n; ++j) {
                const double u = static_cast<double>(i) / n, v = static_cast<double>(j) / n;
                const Vec3 point = p.corner[0] * ((1.0 - u) * (1.0 - v)) +
                                   p.corner[1] * (u * (1.0 - v)) + p.corner[2] * (u * v) +
                                   p.corner[3] * ((1.0 - u) * v);
                grid[static_cast<std::size_t>(i) * static_cast<std::size_t>(n + 1) +
                     static_cast<std::size_t>(j)] = midWeld.weld(point, mid);
            }
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) {
                const auto at = [&](int a, int b) {
                    return grid[static_cast<std::size_t>(a) * static_cast<std::size_t>(n + 1) +
                                static_cast<std::size_t>(b)];
                };
                quads.push_back({{at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1)},
                                 panelIndex});
            }
    }
    if (quads.empty()) {
        report("the zone meshed no elements");
        return patch;
    }

    // --- Nodal normals -----------------------------------------------------------
    //
    // One normal per mid-surface node, area weighted, and the element extruded
    // +/- t/2 along it. On flat plating every element comes out *exactly*
    // prismatic; on curved plating the spread of these normals across an element is
    // the offset ratio the ANS interpolation is charged for, and it is measured
    // below rather than assumed small.
    std::vector<Vec3> nodeNormal(mid.size(), Vec3{0, 0, 0});
    for (const SubQuad& q : quads) {
        Vec3 area = quadNormal(mid[q.node[0]], mid[q.node[1]], mid[q.node[2]], mid[q.node[3]]);
        if (dot(area, axis) < 0) area = -area;
        for (std::uint32_t node : q.node) nodeNormal[node] += area;
    }
    for (std::size_t i = 0; i < nodeNormal.size(); ++i) {
        nodeNormal[i] = normalize(nodeNormal[i]);
        if (length2(nodeNormal[i]) < 0.5) nodeNormal[i] = axis;  // a node with no area
    }

    // --- Node numbering ----------------------------------------------------------
    //
    // Sorted along the patch's long in-plane direction, with the short one running
    // fastest and the two faces of a mid node adjacent. That is the same reasoning
    // as `makePlateMesh`'s node order: it is what keeps the assembled bandwidth
    // proportional to the *short* side, and the static solver the refinement study
    // uses is banded.
    double loLong = 1e300, hiLong = -1e300, loShort = 1e300, hiShort = -1e300;
    for (const Vec3& point : mid) {
        const double a = dot(point - patch.centre, patch.right);
        const double b = dot(point - patch.centre, patch.up);
        loLong = std::min(loLong, a);
        hiLong = std::max(hiLong, a);
        loShort = std::min(loShort, b);
        hiShort = std::max(hiShort, b);
    }
    const bool rightIsLong = (hiLong - loLong) >= (hiShort - loShort);
    std::vector<std::uint32_t> rank(mid.size());
    {
        std::vector<std::uint32_t> byPosition(mid.size());
        for (std::uint32_t i = 0; i < mid.size(); ++i) byPosition[i] = i;
        const auto key = [&](std::uint32_t i) {
            const double a = dot(mid[i] - patch.centre, patch.right);
            const double b = dot(mid[i] - patch.centre, patch.up);
            const double major = rightIsLong ? a : b, minor = rightIsLong ? b : a;
            // Quantised, so that floating-point noise in two nominally identical
            // coordinates cannot interleave two columns.
            return std::pair<long long, long long>{
                static_cast<long long>(std::llround(major * 1e6)),
                static_cast<long long>(std::llround(minor * 1e6))};
        };
        std::sort(byPosition.begin(), byPosition.end(),
                  [&](std::uint32_t a, std::uint32_t b) { return key(a) < key(b); });
        for (std::uint32_t r = 0; r < byPosition.size(); ++r) rank[byPosition[r]] = r;
    }

    const double half = 0.5 * patch.thickness;
    patch.mesh.position.assign(mid.size() * 2 * 3, 0.0);
    patch.outerFace.assign(mid.size() * 2, 0u);
    for (std::uint32_t i = 0; i < mid.size(); ++i) {
        const Vec3 lower = mid[i] - nodeNormal[i] * half;
        const Vec3 upper = mid[i] + nodeNormal[i] * half;
        const std::size_t low = static_cast<std::size_t>(rank[i]) * 2;
        for (int k = 0; k < 3; ++k) {
            patch.mesh.position[low * 3 + static_cast<std::size_t>(k)] = lower[k];
            patch.mesh.position[(low + 1) * 3 + static_cast<std::size_t>(k)] = upper[k];
        }
        patch.outerFace[low + 1] = 1u;
    }

    // --- Elements ---------------------------------------------------------------

    patch.worstJacobian = 1e300;
    double worstSpread = 0;
    for (const SubQuad& q : quads) {
        const Vec3& a = mid[q.node[0]];
        const Vec3& b = mid[q.node[1]];
        const Vec3& c = mid[q.node[2]];
        const Vec3& d = mid[q.node[3]];
        const double area = quadArea(a, b, c, d);
        if (!(area > 0)) continue;

        // Node ordering is a contract: 0-3 on the -zeta face, counter-clockwise
        // seen from +zeta, with 4 above 0. zeta is the thickness direction, so a
        // quad wound the other way has to be reversed or the assumed strains would
        // cure locking *across* the plate instead of through it.
        std::uint32_t corner[4] = {q.node[0], q.node[1], q.node[2], q.node[3]};
        if (dot(quadNormal(a, b, c, d), axis) < 0) std::swap(corner[1], corner[3]);

        // The two faces of this element are the mid-surface quad offset by
        // +/- (t/2) n_a at each corner, so the top face's displacement from the
        // bottom one is t * n_a and its variation across the corners is the face
        // offset the ANS interpolation is charged for. Dividing by t leaves the
        // angle: offset/t is the spread of the nodal normals, in radians.
        Vec3 mean{0, 0, 0};
        for (std::uint32_t node : corner) mean += nodeNormal[node];
        mean = normalize(mean);
        double spread = 0;
        for (std::uint32_t node : corner)
            spread = std::max(spread, length(nodeNormal[node] - mean));
        worstSpread = std::max(worstSpread, spread);
        if (spread > params.normalSpreadWarning) ++patch.distortedElements;

        double edge[4];
        for (int e = 0; e < 4; ++e)
            edge[e] = length(mid[corner[static_cast<std::size_t>((e + 1) % 4)]] -
                             mid[corner[static_cast<std::size_t>(e)]]);
        const double shortest = *std::min_element(edge, edge + 4);
        const double longest = *std::max_element(edge, edge + 4);
        if (shortest > 0) patch.worstAspect = std::max(patch.worstAspect, longest / shortest);

        for (int k = 0; k < 4; ++k)
            patch.mesh.index.push_back(static_cast<std::uint32_t>(rank[corner[k]]) * 2);
        for (int k = 0; k < 4; ++k)
            patch.mesh.index.push_back(static_cast<std::uint32_t>(rank[corner[k]]) * 2 + 1);

        double nodes[kDof];
        patch.mesh.gather(patch.mesh.elementCount() - 1, patch.mesh.position, nodes);
        patch.worstJacobian = std::min(patch.worstJacobian, solidshell::smallestJacobian(nodes));

        patch.panelOf.push_back(q.panel);
        patch.elementArea.push_back(area);
        patch.area += area;
    }
    patch.mass = patch.area * patch.thickness * patch.material.density;
    patch.worstNormalSpread = worstSpread;
    patch.spuriousStiffness = 90.0 * worstSpread * worstSpread;
    if (worstSpread > params.normalSpreadWarning)
        report("the plating turns " + std::to_string(worstSpread) + " rad across " +
               std::to_string(patch.distortedElements) + " of " +
               std::to_string(patch.elementCount()) +
               " elements, so their faces are that far from parallel and the worst is about " +
               std::to_string(100.0 * patch.spuriousStiffness) +
               "% too stiff in bending (07-fem-spike-findings.md §6 limit 1). Refining"
               " MeshParams::subdivision will not help -- the panels are flat facets, so the"
               " turn is at a seam whatever the subdivision; the girth layout in Scantlings is"
               " what sets it");
    if (!(patch.worstJacobian > 0))
        report("an element came out inverted or degenerate; nothing computed on this patch"
               " means anything");

    // Distinct panels and how much of each was meshed, so a tear can be reported as
    // a fraction of the panel rather than as a single dead element.
    {
        std::map<int, double> covered;
        for (std::size_t e = 0; e < patch.panelOf.size(); ++e)
            covered[patch.panelOf[e]] += patch.elementArea[e];
        for (const auto& [index, area] : covered) {
            patch.panels.push_back(index);
            patch.panelArea.push_back(area);
        }
    }

    // --- Boundary ---------------------------------------------------------------
    //
    // An element edge used once is on the perimeter. Clamping both faces of its two
    // mid nodes is a genuine clamp; see the header §4 for why that, and what it
    // costs.
    patch.mesh.fixed.assign(patch.mesh.nodeCount() * 3, 0u);
    patch.mesh.prescribed.assign(patch.mesh.nodeCount() * 3, 0.0);

    // Which mid-surface pairs are joined by an element edge. The perimeter is the
    // edges used once; a stiffener's fibres may only span the edges used at all,
    // so that a member crossing a hole in the patch is broken into runs instead of
    // being bridged by a fibre through thin air.
    std::map<std::pair<std::uint32_t, std::uint32_t>, int> use;
    for (const SubQuad& q : quads)
        for (int e = 0; e < 4; ++e) {
            const std::uint32_t a = q.node[static_cast<std::size_t>(e)];
            const std::uint32_t b = q.node[static_cast<std::size_t>((e + 1) % 4)];
            if (a == b) continue;
            ++use[{std::min(a, b), std::max(a, b)}];
        }

    if (params.edge == Edge::Clamped) {
        int pinnedNodes = 0;
        std::vector<std::uint8_t> onEdge(mid.size(), 0u);
        for (const auto& [edge, count] : use) {
            if (count != 1) continue;
            onEdge[edge.first] = 1u;
            onEdge[edge.second] = 1u;
        }
        for (std::uint32_t i = 0; i < mid.size(); ++i) {
            if (!onEdge[i]) continue;
            ++pinnedNodes;
            for (int side = 0; side < 2; ++side)
                for (int k = 0; k < 3; ++k)
                    patch.mesh.pin(static_cast<std::size_t>(rank[i]) * 2 +
                                       static_cast<std::size_t>(side),
                                   k);
        }
        if (pinnedNodes == static_cast<int>(mid.size()))
            report("every node of the zone is on its clamped perimeter, so there is nothing"
                   " free to deform");
    }

    // --- Stiffeners -------------------------------------------------------------
    //
    // See the header §3. `RigidSupport` pins every plating node the stiffener runs
    // through -- the upper bound; `Modelled` builds the member out of eccentric
    // fibres tied to the plating by `constraint.hpp` and pins nothing. The nodes
    // are exact either way: a member is built along a panel seam, so the grid
    // points on that seam lie on the member's own segment.
    std::vector<std::uint8_t> onStiffener(mid.size(), 0u);
    if (params.stiffeners != Stiffeners::Ignored) {
        patch.stiffening.material = patch.material;
        const double reach = params.radius + 2.0;
        for (std::size_t memberIndex = 0; memberIndex < structure.members.size(); ++memberIndex) {
            const StructuralMember& member = structure.members[memberIndex];
            const Vec3 edge = member.b - member.a;
            const double lengthSquared = length2(edge);
            if (!(lengthSquared > 0)) continue;
            // Distance to the *segment*, not to its midpoint: a girder running the
            // length of the ship has its midpoint far from any zone and still passes
            // straight through one, and a midpoint filter would drop it.
            const double along = std::clamp(dot(impact - member.a, edge) / lengthSquared, 0.0, 1.0);
            if (length(impact - (member.a + edge * along)) > reach) continue;

            // The nodes this member passes through, in order along it. Ordered by
            // the parameter rather than by index, because the mesher numbers nodes
            // for bandwidth and a member may run either way through that numbering.
            std::vector<std::pair<double, std::uint32_t>> onMember;
            for (std::size_t i = 0; i < mid.size(); ++i) {
                const double t =
                    std::clamp(dot(mid[i] - member.a, edge) / lengthSquared, 0.0, 1.0);
                if (length(mid[i] - (member.a + edge * t)) > params.weldTolerance) continue;
                onMember.push_back({t, static_cast<std::uint32_t>(i)});
                onStiffener[i] = 1u;
            }
            if (params.stiffeners != Stiffeners::Modelled || onMember.size() < 2) continue;
            std::sort(onMember.begin(), onMember.end());

            // Which way the web rises against the pair direction. The tie measures
            // the eccentricity along the plating's own thickness direction, so a
            // web that is not perpendicular to the plate has no single offset and
            // is refused rather than projected -- a projected one would be a
            // plausible wrong second moment.
            const Vec3 rise = normalize(member.rise);
            double alignment = 0;
            for (const auto& [t, node] : onMember) alignment += dot(rise, nodeNormal[node]);
            alignment /= static_cast<double>(onMember.size());
            if (std::abs(alignment) < 0.9) {
                report("a stiffener's web rises " + std::to_string(alignment) +
                       " out of the plating's thickness direction, so it has no single"
                       " eccentricity and was left out of the zone");
                continue;
            }

            // Runs of stations actually joined by an element edge. A member that
            // leaves the patch and comes back is two runs, not one fibre spanning
            // the gap.
            constraint::SeamRun run;
            run.sign = alignment > 0 ? 1.0 : -1.0;
            // Which longitudinal these fibres are, so that a torn one can be named
            // back to `StructuralMesh::members` rather than reported as an anonymous
            // loss of steel. `promotion::reactionOf` is the reader.
            run.member = static_cast<int>(memberIndex);
            std::uint32_t previous = 0;
            bool have = false;
            // Counted once per *member*, not once per run: a member the patch
            // carries in two pieces is still one longitudinal, and a count that
            // said two would make `members` a property of the zone's shape.
            bool contributed = false;
            const auto flush = [&]() {
                if (run.bottom.size() >= 2 &&
                    constraint::addStiffener(run, member.profile, patch.thickness,
                                             patch.mesh.position, patch.stiffening) > 0)
                    contributed = true;
                run.bottom.clear();
                run.top.clear();
            };
            for (const auto& [t, node] : onMember) {
                if (have && use.find({std::min(previous, node), std::max(previous, node)}) ==
                                use.end())
                    flush();
                run.bottom.push_back(static_cast<std::uint32_t>(rank[node]) * 2);
                run.top.push_back(static_cast<std::uint32_t>(rank[node]) * 2 + 1);
                previous = node;
                have = true;
            }
            flush();
            if (contributed) ++patch.stiffening.members;
        }
        for (std::size_t i = 0; i < mid.size(); ++i) {
            if (!onStiffener[i]) continue;
            ++patch.stiffenerNodes;
            if (params.stiffeners != Stiffeners::RigidSupport) continue;
            for (int side = 0; side < 2; ++side)
                for (int k = 0; k < 3; ++k)
                    patch.mesh.pin(static_cast<std::size_t>(rank[i]) * 2 +
                                       static_cast<std::size_t>(side),
                                   k);
        }
        if (params.stiffeners == Stiffeners::Modelled && patch.stiffening.empty())
            report("no stiffener reached this zone, so Stiffeners::Modelled is the same"
                   " request as Stiffeners::Ignored");
    }

    std::size_t free = 0;
    for (std::size_t d = 0; d < patch.mesh.fixed.size(); ++d)
        if (!patch.mesh.fixed[d]) ++free;
    patch.freeFraction = patch.mesh.fixed.empty()
                             ? 0.0
                             : static_cast<double>(free) /
                                   static_cast<double>(patch.mesh.fixed.size());
    if (patch.freeFraction < 0.25)
        report("only " + std::to_string(100.0 * patch.freeFraction) +
               "% of the zone's degrees of freedom are free: the stiffener lines and the"
               " perimeter have eaten it, and a subdivision of at least 3 is needed for the"
               " plating between them to be able to deform");

    patch.platingTimestep = solidshell::criticalTimestep(patch.mesh, patch.material,
                                                         solidshell::Formulation::SolidShell);
    patch.criticalTimestep = patch.platingTimestep;

    // The fibres' share of the stability limit, taken rather than assumed. A tie
    // amplifies a fibre's stiffness by the square of its weight -- for the outer
    // fibre of a 200 mm flat bar on 12 mm plating that is a factor of 373 -- while
    // its mass arrives at the pair unamplified, so it is entirely possible for the
    // stiffener rather than the plating to set the step. Whether it does is a
    // property of the element size, and the answer is exact here because a rank-one
    // stiffness against a diagonal mass has a closed-form largest eigenvalue.
    if (!patch.stiffening.empty()) {
        std::vector<double> nodalMass(patch.nodeCount(), 0.0);
        for (std::size_t e = 0; e < patch.elementCount(); ++e) {
            double nodePosition[kDof], lumped[kNodes];
            patch.mesh.gather(e, patch.mesh.position, nodePosition);
            solidshell::elementMass(nodePosition, patch.material.density, lumped);
            for (int a = 0; a < kNodes; ++a)
                nodalMass[patch.mesh.index[e * kNodes + static_cast<std::size_t>(a)]] += lumped[a];
        }
        const constraint::RestFibers forms =
            constraint::restFibers(patch.stiffening, patch.mesh.position);
        constraint::lumpFiberMass(patch.stiffening, forms, patch.material.density, nodalMass);
        const double fibre =
            constraint::criticalTimestep(patch.stiffening, forms, patch.mesh.position, nodalMass,
                                         patch.material.youngsModulus);
        if (fibre > 0 && fibre < patch.criticalTimestep) patch.criticalTimestep = fibre;
    }
    return patch;
}

// --- Cost -----------------------------------------------------------------------

double estimatedCost(const Patch& patch, bool plastic) {
    if (patch.empty() || !(patch.criticalTimestep > 0)) return 0.0;
    const double perElement = plastic ? kPlasticMicroseconds : kElasticMicroseconds;
    const double steps = 1.0 / patch.criticalTimestep;
    return static_cast<double>(patch.elementCount()) * steps * perElement * 1e-6;
}

// --- The striking body ------------------------------------------------------------

double impactSpeed(double energy, double mass) {
    // Not `sqrt(2E/m)` guarded afterwards: a non-positive mass divides by zero and
    // a negative energy takes the root of a negative, and both come back as a NaN
    // that would then be a velocity nothing in the solver could reject. Zero is a
    // punch that does not move, which is a request a caller can see it made.
    if (!(mass > 0) || !(energy > 0)) return 0.0;
    return std::sqrt(2.0 * energy / mass);
}

// --- Solver ---------------------------------------------------------------------

Solver::Solver(const Patch& patch, const plasticity::Material& material, const SolveParams& params)
    : patch_(&patch), material_(material), params_(params) {
    const std::size_t nodes = patch.nodeCount();
    const std::size_t elements = patch.elementCount();

    rest_ = patch.mesh.position;
    position_ = rest_;
    // Before anything is built from `rest_`: the mass, the failure strain and the
    // Gauss volumes are all properties of the *undeformed* configuration, and a
    // pre-loaded patch's undeformed configuration is not the one it was meshed in.
    if (params_.preload.active()) applyPreload();
    velocity_.assign(nodes * 3, 0.0);
    force_.assign(nodes * 3, 0.0);
    mass_.assign(nodes, 0.0);
    elementForce_.assign(elements * kDof, 0.0);
    dissipation_.assign(elements, 0.0);
    pinned_.assign(nodes * 3, 0u);
    for (std::size_t d = 0; d < nodes * 3 && d < patch.mesh.fixed.size(); ++d)
        pinned_[d] = patch.mesh.fixed[d];

    if (params_.plastic) {
        elementStress_.assign(elements * kGauss * 6, 0.0);
        gaussVolume_.assign(elements * kGauss, 0.0);
        plastic_.resize(elements);
    } else {
        stiffness_.assign(elements * kDof * kDof, 0.0);
    }

    // The step-invariant half of every element, formed once. See
    // `SolveParams::cacheRestForms`: this is half of the per-step element kernel
    // and it is a function of `rest_`, which nothing after this point moves.
    result_.cachedRestForms = params_.cacheRestForms;
    if (result_.cachedRestForms) forms_.resize(elements);

    for (std::size_t e = 0; e < elements; ++e) {
        double nodePosition[kDof];
        patch.mesh.gather(e, rest_, nodePosition);
        // From `rest_`, not from `position_`. Without a `Preload` the two are the
        // same array and any mistake here is invisible; with one they differ by the
        // pre-strain, and forms built from the deformed configuration would silently
        // give every element the wrong B. Mutation testing found this exact edit
        // surviving the whole suite.
        if (result_.cachedRestForms)
            solidshell::computeRestForms(nodePosition, solidshell::Formulation::SolidShell,
                                         forms_[e]);
        double lumped[kNodes];
        solidshell::elementMass(nodePosition, patch.material.density, lumped);
        for (int a = 0; a < kNodes; ++a)
            mass_[patch.mesh.index[e * kNodes + static_cast<std::size_t>(a)]] += lumped[a];
        if (params_.plastic) {
            solidshell::initialisePlasticState(nodePosition, material_, plastic_[e]);
            solidshell::gaussVolumes(nodePosition, &gaussVolume_[e * kGauss]);
        } else {
            solidshell::elementStiffness(nodePosition, patch.material,
                                         solidshell::Formulation::SolidShell,
                                         &stiffness_[e * kDof * kDof]);
        }
    }

    // The stiffeners. `rest_`, not `patch.mesh.position`: a `Preload` has already
    // moved the rest configuration by this point, and a fibre whose rest length
    // came from the meshed geometry would be handed the pre-strain for free -- the
    // same mistake, in the same place, that mutation testing found in the element
    // forms above.
    if (!patch.stiffening.empty()) {
        fiberForms_ = constraint::restFibers(patch.stiffening, rest_);
        if (!fiberForms_.ok)
            result_.problems.push_back("a stiffener fibre came out with no length");
        constraint::lumpFiberMass(patch.stiffening, fiberForms_, patch.material.density, mass_);
        if (params_.plastic) fiber_.assign(patch.stiffening.fiberCount(), {});
    }

    // Node -> incident element corners, as CSR. Gathering rather than scattering
    // fixes the accumulation order, which is what makes the parallel answer
    // bit-identical to the serial one -- the same reason `fem.cpp` does it.
    adjacencyOffset_.assign(nodes + 1, 0u);
    for (std::size_t e = 0; e < elements; ++e)
        for (int a = 0; a < kNodes; ++a)
            ++adjacencyOffset_[patch.mesh.index[e * kNodes + static_cast<std::size_t>(a)] + 1];
    for (std::size_t i = 0; i < nodes; ++i) adjacencyOffset_[i + 1] += adjacencyOffset_[i];
    adjacencyEntry_.assign(elements * kNodes, 0u);
    {
        std::vector<std::uint32_t> cursor(adjacencyOffset_.begin(), adjacencyOffset_.end() - 1);
        for (std::size_t e = 0; e < elements; ++e)
            for (int a = 0; a < kNodes; ++a) {
                const std::uint32_t node =
                    patch.mesh.index[e * kNodes + static_cast<std::size_t>(a)];
                adjacencyEntry_[cursor[node]++] =
                    static_cast<std::uint32_t>(e) * kNodes + static_cast<std::uint32_t>(a);
            }
    }

    // The punch's footprint: outer-face nodes inside a rectangle in the patch's own
    // in-plane frame. Outer face only, so the plate may thin under it.
    const bool punchPresent =
        params_.indenter.halfLength > 0 && params_.indenter.halfWidth > 0;
    for (std::size_t i = 0; punchPresent && i < nodes; ++i) {
        if (i < patch.outerFace.size() && patch.outerFace[i] == 0u) continue;
        const Vec3 point{position_[i * 3], position_[i * 3 + 1], position_[i * 3 + 2]};
        const Vec3 offset = point - patch.centre;
        if (std::abs(dot(offset, patch.right)) > params_.indenter.halfLength) continue;
        if (std::abs(dot(offset, patch.up)) > params_.indenter.halfWidth) continue;
        if (pinned_[i * 3] || pinned_[i * 3 + 1] || pinned_[i * 3 + 2]) continue;
        driven_.push_back(static_cast<std::uint32_t>(i));
    }
    if (punchPresent && driven_.empty())
        result_.problems.push_back("the indenter's footprint contains no free node");

    // The striking body. See `zone.hpp` §6: the mass and the arrival speed are the
    // entry point, and the energy is what they carry rather than a third input that
    // could disagree with them.
    if (params_.indenter.drive == Drive::Inertial) {
        if (!(params_.indenter.mass > 0) || !(params_.indenter.speed > 0)) {
            // Refused rather than fallen back on. A striker with no mass or no
            // approach speed carries no energy, and running it as a prescribed punch
            // instead would hand back a travel the caller asked for an energy to
            // decide -- which is exactly the assumption the inertial drive exists to
            // remove. `impactSpeed` returns zero for an energy or a mass it cannot
            // use, so this is the shape a mis-piped collision arrives in and it has
            // to be loud.
            //
            // **The speed is refused here rather than left to arrest on the first
            // step, and that was a measurement.** A striker released at zero was
            // meant to stop immediately, because `arrested` is `!(speed > 0)`. It
            // does not: the internal force at the rest configuration is a rounding
            // residue rather than an identical zero, so the grip hands back a speed
            // around 1e-29 m/s and the run goes to `maxSteps` having travelled
            // 1e-33 m. Deciding a run's whole termination on an exact floating-point
            // zero is the trapdoor CLAUDE.md records `sectionElements` falling
            // through, and the fix is to refuse the input rather than to widen the
            // test into a tolerance nothing measures.
            result_.problems.push_back(
                "an inertial indenter needs a positive striking mass and approach speed;"
                " nothing was solved");
            done_ = true;
        } else if (driven_.empty()) {
            result_.problems.push_back(
                "an inertial indenter that touches nothing would never stop; nothing was solved");
            done_ = true;
        } else if (!(params_.indenter.stopAt > 0) && !(params_.duration > 0)) {
            // **Refused, and this is the one place the inertial drive really can run
            // on.** Its *speed* is bounded -- nothing does work on the punch, so
            // `1/2 m v^2 <= 1/2 m v_0^2` for the whole run -- but its *travel* is
            // not: a punch that has perforated the plating under it meets no force
            // at all and coasts at very nearly its arrival speed for ever. Measured:
            // twice the energy needed to tear the reference strip put the punch
            // 24.9 m past the plating and ended on `maxSteps`. So an inertial drive
            // must be bounded by something that is not a step count, and it is made
            // to say so at construction rather than after two million steps.
            result_.problems.push_back(
                "an inertial indenter needs a travel or a duration to bound it: a punch that"
                " perforates coasts, and maxSteps is a budget rather than a bound. Nothing was"
                " solved");
            done_ = true;
        } else {
            punchMass_ = params_.indenter.mass;
            punchSpeed_ = params_.indenter.speed;
            for (std::uint32_t node : driven_) drivenMass_ += mass_[node];
        }
        if (params_.indenter.rampTime > 0)
            result_.problems.push_back(
                "rampTime is a prescribed punch's statement about compliance and an inertial"
                " punch arrives travelling; it was ignored");
        result_.indenterMass = punchMass_;
        result_.indenterSpeed = punchSpeed_;
        result_.indenterEnergy = 0.5 * punchMass_ * punchSpeed_ * punchSpeed_;
        result_.indenterKinetic = result_.indenterEnergy;
    }

    // The driven perimeter. A pinned DOF whose prescribed value is zero is a
    // clamp and takes the cheap path; one with a value follows the surrounding
    // structure -- see `zone.hpp` §4 and `coupling.hpp`. Testing the value rather
    // than carrying a separate flag is not a shortcut: a DOF driven to exactly
    // zero *is* a clamp, in position, in velocity and in the work it does.
    for (std::size_t d = 0; d < nodes * 3 && d < patch.mesh.prescribed.size(); ++d)
        if (pinned_[d] && patch.mesh.prescribed[d] != 0.0)
            edgeDof_.push_back(static_cast<std::uint32_t>(d));
    result_.drivenEdgeDof = static_cast<int>(edgeDof_.size());
    if (!edgeDof_.empty()) edgeFree_.assign(nodes * 3, 0.0);

    result_.timestep = params_.timestep > 0
                           ? params_.timestep
                           : patch.criticalTimestep * (params_.timestepSafety / 0.9);
    if (!(result_.timestep > 0)) {
        result_.problems.push_back("no stable timestep: the patch has no usable element");
        done_ = true;
    }
    for (const std::string& problem : patch.problems) result_.problems.push_back(problem);

    // Prime the element state, so the stress the patch is *handed* exists before
    // the first step rather than appearing during it. Without this a pre-loaded
    // patch reports zero mean stress and zero strain energy right up until step
    // one, and the energy account then sees the pre-load's stored energy arrive
    // as work nobody did. Costs one force evaluation, once.
    //
    // Whatever plastic work the ship did to the patch before it was promoted is
    // the hull girder's, not the zone's, so the dissipation baseline is reset to
    // zero here: `SolveResult::dissipation` counts what the *indenter* spends.
    computeForces();
    result_.dissipation = 0.0;
    accumulateEnergy();
    result_.initialStrainEnergy = result_.strainEnergy;
    // The priming evaluation above is a promotion cost, like the stiffness
    // formation and the power iteration beside it, and counting it as a step would
    // put one element pass into a profile whose whole purpose is the *per-step*
    // ratio. Zeroed here so `profile` means "the run" and nothing else.
    result_.profile = SolveResult::Profile{};
}

// The pre-load, as an initial strain rather than an initial stress.
//
// The strain field of a uniaxial `sigma_xx(z) = stress + gradient * (z - ref)` is
//
//     eps_xx = e0 + k h,   eps_yy = eps_zz = -nu (e0 + k h),   h = z - ref
//
// with e0 = stress/E and k = gradient/E, and it is compatible: the displacement
// field below reproduces it exactly, which is the classical pure-bending solution
// of elasticity with a uniform axial part added. So the patch is handed a state
// that is *equilibrated* -- the interior nodal forces vanish identically and only
// the clamped edge carries reaction -- and a pre-loaded patch with no punch on it
// therefore does not move. `tests/test_promotion.cpp` asserts that, because a
// pre-load that quietly rings is worse than none.
//
// `rest = position - u(position)` rather than solving `position = rest + u(rest)`:
// the two differ at second order in the strain, which at 84 MPa in steel is one
// part in 2.4 million.
void Solver::applyPreload() {
    const Preload& pre = params_.preload;
    const double youngs = patch_->material.youngsModulus;
    if (!(youngs > 0)) {
        result_.problems.push_back("cannot pre-load a patch whose material has no stiffness");
        return;
    }
    const double nu = patch_->material.poissonRatio;
    const double e0 = pre.stress / youngs;
    const double k = pre.gradient / youngs;
    const Vec3 centre = patch_->centre;

    for (std::size_t node = 0; node < patch_->nodeCount(); ++node) {
        const double x = position_[node * 3] - centre.x;
        const double y = position_[node * 3 + 1] - centre.y;
        const double h = position_[node * 3 + 2] - pre.reference;
        const double strain = e0 + k * h;
        rest_[node * 3] -= strain * x;
        rest_[node * 3 + 1] -= -nu * strain * y;
        rest_[node * 3 + 2] -=
            -nu * (e0 * h + 0.5 * k * h * h) - 0.5 * k * x * x + 0.5 * nu * k * y * y;
    }
}

namespace {
// One clock read. Named so the profiling calls read as measurements rather than
// as bookkeeping.
inline std::chrono::steady_clock::time_point tick() {
    return std::chrono::steady_clock::now();
}
inline double since(std::chrono::steady_clock::time_point& mark) {
    const auto now = std::chrono::steady_clock::now();
    const double seconds = std::chrono::duration<double>(now - mark).count();
    mark = now;
    return seconds;
}
}  // namespace

void Solver::computeForces() {
    const std::size_t elements = patch_->elementCount();
    auto mark = tick();
    const auto body = [&](std::size_t e) {
        double rest[kDof], current[kDof];
        patch_->mesh.gather(e, rest_, rest);
        patch_->mesh.gather(e, position_, current);
        double* out = &elementForce_[e * kDof];
        if (!params_.plastic) {
            if (result_.cachedRestForms)
                solidshell::internalForce(forms_[e], &stiffness_[e * kDof * kDof], rest, current,
                                          out);
            else
                solidshell::internalForce(&stiffness_[e * kDof * kDof], rest, current, out);
            return;
        }
        if (plastic_[e].torn) {
            // Output-identical to letting the update run -- every point has failed, so
            // `plasticity::update` returns no stress and no dissipation -- and it
            // saves the whole 7 µs for the rest of a dead element's life. The same
            // argument as `solid_shell.cpp`'s pre-loop check for a degraded element.
            std::fill(out, out + kDof, 0.0);
            double* stress = elementStress_.data() + e * kGauss * 6;
            std::fill(stress, stress + kGauss * 6, 0.0);
            dissipation_[e] = 0.0;
            return;
        }
        const solidshell::PlasticUpdate update =
            result_.cachedRestForms
                ? solidshell::elementPlasticUpdate(forms_[e], rest, current, material_,
                                                   plastic_[e], out,
                                                   &elementStress_[e * kGauss * 6])
                : solidshell::elementPlasticUpdate(rest, current, material_,
                                                   solidshell::Formulation::SolidShell,
                                                   plastic_[e], out,
                                                   &elementStress_[e * kGauss * 6]);
        dissipation_[e] = update.dissipation;
    };

    if (params_.jobs != nullptr && elements >= 64) {
        params_.jobs->parallelFor(0, elements, 8, [&](std::size_t begin, std::size_t end) {
            for (std::size_t e = begin; e < end; ++e) body(e);
        });
    } else {
        for (std::size_t e = 0; e < elements; ++e) body(e);
    }
    result_.profile.element += since(mark);

    const std::size_t nodes = patch_->nodeCount();
    for (std::size_t node = 0; node < nodes; ++node) {
        double sum[3] = {0, 0, 0};
        for (std::uint32_t k = adjacencyOffset_[node]; k < adjacencyOffset_[node + 1]; ++k) {
            const std::uint32_t entry = adjacencyEntry_[k];
            const double* f = &elementForce_[(entry / kNodes) * kDof + (entry % kNodes) * 3];
            sum[0] += f[0];
            sum[1] += f[1];
            sum[2] += f[2];
        }
        force_[node * 3] = sum[0];
        force_[node * 3 + 1] = sum[1];
        force_[node * 3 + 2] = sum[2];
    }
    // The stiffeners, after the plating's gather and in a fixed order, so the
    // answer stays bit-identical whatever the worker count. There are a couple of
    // hundred fibres against tens of thousands of element-corner writes, so this is
    // not where a thread would go anyway.
    fiberEnergy_ = 0;
    if (!patch_->stiffening.empty()) {
        const constraint::FiberForces fibres =
            constraint::fiberForces(patch_->stiffening, fiberForms_, position_, material_,
                                    params_.plastic ? &fiber_ : nullptr, force_,
                                    params_.fiberFailure);
        fiberEnergy_ = fibres.strainEnergy;
        result_.dissipation += fibres.dissipation;
    }
    if (params_.plastic) {
        int torn = 0;
        for (std::size_t e = 0; e < elements; ++e) {
            result_.dissipation += dissipation_[e];
            if (plastic_[e].torn) ++torn;
        }
        result_.tornElements = torn;
    }
    result_.profile.gather += since(mark);
}

void Solver::accumulateEnergy() {
    double kinetic = 0;
    for (std::size_t node = 0; node < patch_->nodeCount(); ++node) {
        const double* v = &velocity_[node * 3];
        kinetic += 0.5 * mass_[node] * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    }
    result_.kinetic = kinetic;

    double strain = 0;
    if (params_.plastic) {
        // U = sum_gp w_gp * 1/2 sigma : C^-1 : sigma, which needs no strain and no B
        // matrix -- and is therefore right for an element whose total strain has a
        // plastic part in it. A torn point carries no stress, so its stored energy
        // leaves the account discontinuously, which is what element deletion is.
        const double mu = material_.shearModulus(), bulk = material_.bulkModulus();
        for (std::size_t e = 0; e < patch_->elementCount(); ++e)
            for (int gp = 0; gp < kGauss; ++gp) {
                const double* s = &elementStress_[e * kGauss * 6 + static_cast<std::size_t>(gp) * 6];
                const double pressure = (s[0] + s[1] + s[2]) / 3.0;
                double density = 0;
                for (int i = 0; i < 3; ++i)
                    density += 0.5 * s[i] * ((s[i] - pressure) / (2.0 * mu) + pressure / (3.0 * bulk));
                for (int i = 3; i < 6; ++i) density += 0.5 * s[i] * (s[i] / mu);
                strain += gaussVolume_[e * kGauss + static_cast<std::size_t>(gp)] * density;
            }
    } else {
        for (std::size_t e = 0; e < patch_->elementCount(); ++e) {
            double rest[kDof], current[kDof];
            patch_->mesh.gather(e, rest_, rest);
            patch_->mesh.gather(e, position_, current);
            double rotation[9];
            if (result_.cachedRestForms)
                solidshell::elementRotation(forms_[e], current, rotation);
            else
                solidshell::elementRotation(rest, current, rotation);
            double u[kDof];
            for (int a = 0; a < kNodes; ++a)
                for (int i = 0; i < 3; ++i) {
                    double s = 0;
                    for (int k = 0; k < 3; ++k) s += rotation[i * 3 + k] * current[a * 3 + k];
                    u[a * 3 + i] = s - rest[a * 3 + i];
                }
            const double* k = &stiffness_[e * kDof * kDof];
            for (int i = 0; i < kDof; ++i) {
                double s = 0;
                for (int j = 0; j < kDof; ++j) s += k[i * kDof + j] * u[j];
                strain += 0.5 * u[i] * s;
            }
        }
    }
    // The fibres' stored energy, which `computeForces` has already accumulated: it
    // is `sigma^2 A L / 2E` per fibre and needs no B matrix, exactly as the
    // plating's does not. Added here rather than reported separately because the
    // energy balance is a statement about the whole patch, and a stiffener that
    // stored energy outside the account would show up as the solver inventing it.
    result_.strainEnergy = strain + fiberEnergy_;
}

// How much of a driven perimeter's displacement is imposed at `time`. Smoothstep,
// so both ends of the ramp are at zero velocity -- see `SolveParams::edgeRamp`.
//
// The lower clamp is unreachable from the only call site, which passes
// `result_.time + dt` and `result_.time` starts at zero and only rises; mutation
// testing confirms that removing it changes nothing. It stays because this reads
// as a pure function of time, and a smoothstep that returned a negative fraction
// for a negative argument would drive the boundary the wrong way.
double Solver::edgeFraction(double time) const {
    if (!(params_.edgeRamp > 0)) return 1.0;
    const double s = std::min(1.0, std::max(0.0, time / params_.edgeRamp));
    return s * s * (3.0 - 2.0 * s);
}

bool Solver::step() {
    if (done_) return false;
    const double dt = result_.timestep;
    const Indenter& punch = params_.indenter;

    computeForces();
    auto mark = tick();

    for (std::size_t node = 0; node < patch_->nodeCount(); ++node)
        for (int k = 0; k < 3; ++k) {
            const std::size_t d = node * 3 + static_cast<std::size_t>(k);
            if (pinned_[d]) {
                // The velocity this DOF would have taken if nothing held it. It is
                // discarded on the next line, so a driven edge has to take its copy
                // here; what the boundary applies is the difference, exactly as the
                // punch's impulse is the difference for a node it grips.
                if (!edgeFree_.empty()) {
                    double freeVelocity = velocity_[d];
                    if (mass_[node] > 0) freeVelocity += force_[d] / mass_[node] * dt;
                    edgeFree_[d] = freeVelocity;
                }
                velocity_[d] = 0.0;
                continue;
            }
            if (mass_[node] > 0) velocity_[d] += force_[d] / mass_[node] * dt;
        }

    if (params_.damping != 1.0) {
        double before = 0;
        for (std::size_t node = 0; node < patch_->nodeCount(); ++node) {
            const double* v = &velocity_[node * 3];
            before += 0.5 * mass_[node] * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
        }
        for (double& v : velocity_) v *= params_.damping;
        result_.dampingLoss += before * (1.0 - params_.damping * params_.damping);
    }

    // The driven perimeter, after the damping and before the punch. After the
    // damping because the boundary's motion is imposed rather than computed and
    // scaling it would be scaling the surrounding ship's answer; before the punch
    // only for reading, since `driven_` excludes every pinned node and the two sets
    // cannot overlap.
    //
    // The target is taken from the **meshed** position rather than from `rest_`: a
    // `Preload` moves the rest configuration out from under the mesh, and a drive
    // measured from there would carry the pre-load's own displacement field into
    // the boundary condition. It is the same trap `restFibers` and the element
    // forms above each record.
    if (!edgeDof_.empty()) {
        const double fraction = edgeFraction(result_.time + dt);
        for (std::uint32_t dof : edgeDof_) {
            const std::size_t d = dof;
            const std::size_t node = d / 3;
            const double target =
                (patch_->mesh.position[d] + patch_->mesh.prescribed[d] * fraction - position_[d]) /
                dt;
            result_.boundaryWork += mass_[node] * (target - edgeFree_[d]) * target;
            velocity_[d] = target;
        }
    }

    // The punch. Its speed rises over `rampTime` and is constant after; the impulse
    // it has to apply to hold a node at that speed *is* the resisting force, and
    // taking it that way rather than as the internal force alone keeps the node's
    // own inertia in the account, which is what makes the energy balance close.
    double speed = punch.speed;
    if (punch.rampTime > 0 && result_.time < punch.rampTime)
        speed = punch.speed * (result_.time / punch.rampTime);
    // The striking body, when there is one. `punchMass_` is zero on the prescribed
    // drive, so everything below this block is the arithmetic that path always ran,
    // in the order it ran it -- which is what lets the tests assert bit-identity
    // rather than agreement.
    //
    // **The grip is a perfectly inelastic collision, and taking it that way is what
    // makes the infinite-mass limit exact.** The punch and the nodes it holds are
    // rigidly attached, so momentum along the travel is what survives the step:
    //
    //     v' = (m_punch v + sum_i m_i u_i) / (m_punch + sum_i m_i)
    //
    // As `m_punch` grows this tends to `v` and the impulse tends to
    // `sum_i m_i (v - u_i)`, which is the expression the prescribed punch has
    // always used -- so the two drives are one formulation and not two, and a heavy
    // enough inertial punch reproduces a prescribed run rather than resembling it.
    // A separate `m dv/dt = -F` with F taken from the previous step would be neither
    // momentum-conserving nor convergent to that limit.
    if (punchMass_ > 0) {
        double nodeMomentum = 0;
        for (std::uint32_t node : driven_) {
            double along = 0;
            for (int k = 0; k < 3; ++k) along += velocity_[node * 3 + static_cast<std::size_t>(k)] *
                                                 patch_->axis[k];
            // `along` is the +axis component and the punch travels along -axis, so
            // the node's speed *with* the punch is its negative.
            nodeMomentum -= mass_[node] * along;
        }
        punchSpeed_ = (punchMass_ * punchSpeed_ + nodeMomentum) / (punchMass_ + drivenMass_);
        speed = punchSpeed_;
        result_.indenterSpeed = punchSpeed_;
        result_.indenterKinetic = 0.5 * punchMass_ * punchSpeed_ * punchSpeed_;
    }
    double impulse = 0;
    for (std::uint32_t node : driven_) {
        double along = 0;
        for (int k = 0; k < 3; ++k) along += velocity_[node * 3 + static_cast<std::size_t>(k)] *
                                             patch_->axis[k];
        const double change = -speed - along;
        for (int k = 0; k < 3; ++k)
            velocity_[node * 3 + static_cast<std::size_t>(k)] += change * patch_->axis[k];
        impulse += mass_[node] * change;
    }
    result_.force = -impulse / dt;
    result_.peakForce = std::max(result_.peakForce, result_.force);
    result_.work += result_.force * speed * dt;

    for (std::size_t d = 0; d < position_.size(); ++d) position_[d] += velocity_[d] * dt;
    result_.penetration += speed * dt;
    result_.time += dt;
    ++result_.steps;
    result_.profile.integrate += since(mark);

    accumulateEnergy();
    result_.profile.energy += since(mark);

    if (params_.historyStride > 0 && result_.steps % params_.historyStride == 0)
        result_.history.push_back({result_.time, result_.penetration, result_.force, result_.work,
                                   result_.strainEnergy, result_.dissipation, result_.kinetic,
                                   // `speed` and not a test on the drive: the grip
                                   // above assigns `speed = punchSpeed_`, so the
                                   // conditional this used to carry was dead. A
                                   // mutant that removed it changed nothing, which
                                   // is how it was found.
                                   speed, result_.tornElements});

    // The striking body has stopped: it spent what it arrived with and the
    // penetration above is what that bought. Tested *before* the travel cap, so a
    // run that satisfies both is reported as ended by the energy -- which is the
    // truthful reading, since a striker with nothing left would not have gone
    // further whatever the cap said. And tested after the step, so the depth
    // reported is the deepest reached rather than the one before it.
    const bool arrested = punchMass_ > 0 && !(punchSpeed_ > 0);
    const bool reachedDepth = punch.stopAt > 0 && result_.penetration >= punch.stopAt;
    const bool reachedTime = params_.duration > 0 && result_.time >= params_.duration;
    if (arrested) {
        done_ = true;
        result_.completed = true;
        result_.indenterArrested = true;
    } else if (reachedDepth || reachedTime) {
        done_ = true;
        result_.completed = true;
        // A travel cap that ends an energy-driven run has decided the answer the
        // energy was meant to decide, and silently. That is the failure
        // `ImpactDamage::energyUnspent` guards on the membrane path, and it is
        // guarded here rather than left to be inferred from `indenterArrested`.
        if (punchMass_ > 0)
            result_.problems.push_back(
                std::string("the run ended on the ") + (reachedDepth ? "travel" : "duration") +
                " cap with the striker still carrying " +
                std::to_string(result_.indenterKinetic) +
                " J; the cap decided the penetration, not the energy");
    } else if (result_.steps >= params_.maxSteps) {
        done_ = true;
        result_.problems.push_back("stopped at maxSteps with the indenter at " +
                                   std::to_string(result_.penetration) + " m");
    }
    result_.profile.other += since(mark);
    return !done_;
}

std::vector<int> Solver::tornPanelsAt(double fraction) const {
    std::vector<int> torn;
    if (!params_.plastic) return torn;
    std::map<int, double> deleted;
    for (std::size_t e = 0; e < patch_->elementCount(); ++e)
        if (plastic_[e].torn) deleted[patch_->panelOf[e]] += patch_->elementArea[e];
    for (std::size_t i = 0; i < patch_->panels.size(); ++i) {
        auto found = deleted.find(patch_->panels[i]);
        if (found == deleted.end()) continue;
        if (found->second < fraction * patch_->panelArea[i]) continue;
        torn.push_back(patch_->panels[i]);
    }
    return torn;
}

void Solver::collectTorn() {
    result_.tornElements = 0;
    result_.yieldedElements = 0;
    result_.tornArea = 0;
    for (std::size_t e = 0; e < patch_->elementCount(); ++e) {
        if (!params_.plastic) continue;
        if (plastic_[e].torn) {
            ++result_.tornElements;
            result_.tornArea += patch_->elementArea[e];
        }
        for (int gp = 0; gp < kGauss; ++gp)
            if (plastic_[e].point[gp].equivalentPlasticStrain > 0) {
                ++result_.yieldedElements;
                break;
            }
    }
    // The stiffeners' half, taken from the *committed* fibre state for the same
    // reason the plating's is taken from `plastic_` here rather than kept from
    // `computeForces`: `adopt` runs a force evaluation on a throwaway copy to
    // recover the stress, and a fibre sitting just under its damage limit would be
    // tipped over it by the act of being read. Counting here counts what was
    // adopted.
    result_.tornFibers = 0;
    result_.tornFiberVolume = 0;
    for (std::size_t i = 0; i < fiber_.size(); ++i) {
        if (!fiber_[i].failed) continue;
        ++result_.tornFibers;
        if (i < fiberForms_.length.size() && i < patch_->stiffening.fiber.size())
            result_.tornFiberVolume += patch_->stiffening.fiber[i].area * fiberForms_.length[i];
    }
    result_.tornPanels = tornPanelsAt(params_.tearFraction);
    // The area `breach.hpp` will open is the whole panel, not the deleted part.
    result_.tornPanelArea = 0;
    for (int index : result_.tornPanels)
        for (std::size_t i = 0; i < patch_->panels.size(); ++i)
            if (patch_->panels[i] == index) result_.tornPanelArea += patch_->panelArea[i];
}

const SolveResult& Solver::run() {
    const auto begin = std::chrono::steady_clock::now();
    while (step()) {
    }
    const auto end = std::chrono::steady_clock::now();
    result_.wallSeconds = std::chrono::duration<double>(end - begin).count();
    const double elementSteps =
        static_cast<double>(result_.steps) * static_cast<double>(patch_->elementCount());
    result_.microsecondsPerElementStep =
        elementSteps > 0 ? result_.wallSeconds * 1e6 / elementSteps : 0.0;
    collectTorn();
    return result_;
}

void Solver::adopt(const std::vector<double>& position, const std::vector<double>& velocity,
                   const std::vector<solidshell::ElementPlasticState>& state, int steps,
                   double time, double penetration, double work, double dissipation) {
    if (position.size() == position_.size()) position_ = position;
    if (velocity.size() == velocity_.size()) velocity_ = velocity;
    if (params_.plastic && state.size() == plastic_.size()) plastic_ = state;
    // The accelerator drives a prescribed punch and knows nothing about a striking
    // mass, so a state adopted onto an inertial drive carries a punch speed that was
    // never decelerated. Said rather than absorbed, for the same reason the missing
    // stiffeners below are: the answer would look like a slightly deep zone.
    if (punchMass_ > 0)
        result_.problems.push_back(
            "a state was adopted onto an inertial indenter; the striker's speed is the one this"
            " solver last computed and not the one the other path ran");
    result_.steps = steps;
    result_.time = time;
    result_.penetration = penetration;
    result_.work = work;
    result_.dissipation = dissipation;
    // The stress the energy account needs is not carried across -- it is a function
    // of the state, so it is recomputed here through the same path every CPU step
    // uses. That also leaves `result_.tornElements` correct.
    //
    // **On a copy of the history, because `elementPlasticUpdate` commits.** Running
    // it to recover the stress would otherwise advance every integration point by a
    // further increment -- small, but a point sitting just under its damage limit
    // would be tipped over it by the act of being read, and the adopted state would
    // no longer be the state that was adopted. Found by writing the test that says
    // adopting a run's own output reproduces that run.
    const std::vector<solidshell::ElementPlasticState> committed = plastic_;
    const std::vector<constraint::FiberState> committedFibers = fiber_;
    computeForces();
    plastic_ = committed;
    fiber_ = committedFibers;
    result_.dissipation = dissipation;
    // The accelerator does not carry the stiffeners: `gpu::ZoneGpuSolver` uploads
    // the element arrays and nothing else, so a patch with fibres in it that was
    // stepped elsewhere has been stepped without them. Said rather than silently
    // absorbed, because the answer would look like a slightly soft zone.
    if (!patch_->stiffening.empty())
        result_.problems.push_back("this state was computed by a path that does not carry the"
                                   " zone's " + std::to_string(patch_->stiffening.fiberCount()) +
                                   " stiffener fibres, so they took no part in it");
    accumulateEnergy();
    collectTorn();
    done_ = true;
    result_.completed = true;
}

void Solver::translate(const Vec3& velocity) {
    for (std::size_t node = 0; node < patch_->nodeCount(); ++node)
        for (int k = 0; k < 3; ++k) velocity_[node * 3 + static_cast<std::size_t>(k)] += velocity[k];
    accumulateEnergy();
}

void Solver::meanStress(double out[6]) const {
    for (int i = 0; i < 6; ++i) out[i] = 0.0;
    double volume = 0;
    const std::size_t elements = patch_->elementCount();

    if (params_.plastic) {
        for (std::size_t e = 0; e < elements; ++e)
            for (int gp = 0; gp < kGauss; ++gp) {
                const double w = gaussVolume_[e * kGauss + static_cast<std::size_t>(gp)];
                const double* s = &elementStress_[e * kGauss * 6 + static_cast<std::size_t>(gp) * 6];
                volume += w;
                for (int i = 0; i < 6; ++i) out[i] += w * s[i];
            }
    } else {
        // The elastic path keeps no stress, so it is recovered from the
        // displacement -- which is the same route `elementStress` serves the patch
        // test by, and exact for the small strains a pre-load applies.
        for (std::size_t e = 0; e < elements; ++e) {
            double rest[kDof], current[kDof], displacement[kDof];
            patch_->mesh.gather(e, rest_, rest);
            patch_->mesh.gather(e, position_, current);
            for (int d = 0; d < kDof; ++d) displacement[d] = current[d] - rest[d];
            double stress[kGauss * 6], weight[kGauss];
            solidshell::elementStress(rest, displacement, patch_->material,
                                      solidshell::Formulation::SolidShell, stress);
            solidshell::gaussVolumes(rest, weight);
            for (int gp = 0; gp < kGauss; ++gp) {
                volume += weight[gp];
                for (int i = 0; i < 6; ++i) out[i] += weight[gp] * stress[gp * 6 + i];
            }
        }
    }
    if (volume > 0)
        for (int i = 0; i < 6; ++i) out[i] /= volume;
}

double Solver::largestDisplacement() const {
    double worst = 0;
    for (std::size_t node = 0; node < patch_->nodeCount(); ++node) {
        double d = 0;
        for (int k = 0; k < 3; ++k) {
            const double e = position_[node * 3 + static_cast<std::size_t>(k)] -
                             rest_[node * 3 + static_cast<std::size_t>(k)];
            d += e * e;
        }
        worst = std::max(worst, std::sqrt(d));
    }
    return worst;
}

// --- The whole chain -------------------------------------------------------------

ZoneDamage indent(const StructuralMesh& structure, const Vec3& impact, const MeshParams& mesh,
                  const SolveParams& solve, const plasticity::Material& material) {
    ZoneDamage damage;
    damage.patch = buildPatch(structure, impact, mesh);
    if (damage.patch.empty()) {
        damage.result.problems = damage.patch.problems;
        return damage;
    }
    Solver solver(damage.patch, material, solve);
    damage.result = solver.run();
    return damage;
}

}  // namespace sim::zone
