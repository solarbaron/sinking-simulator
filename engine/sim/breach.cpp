// SPDX-License-Identifier: MIT
#include "breach.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <set>
#include <unordered_map>
#include <utility>

namespace sim {
namespace {

constexpr double kInf = std::numeric_limits<double>::infinity();

void boundingBox(const TriMesh& m, Vec3& lo, Vec3& hi) {
    lo = {kInf, kInf, kInf};
    hi = {-kInf, -kInf, -kInf};
    for (const Vec3& v : m.verts) {
        lo = {std::min(lo.x, v.x), std::min(lo.y, v.y), std::min(lo.z, v.z)};
        hi = {std::max(hi.x, v.x), std::max(hi.y, v.y), std::max(hi.z, v.z)};
    }
}

// An empty mesh leaves lo > hi, so this is false for it and the winding sum is
// never asked about a mesh that has no surface.
bool insideBox(const Vec3& p, const Vec3& lo, const Vec3& hi) {
    return p.x >= lo.x && p.x <= hi.x && p.y >= lo.y && p.y <= hi.y && p.z >= lo.z && p.z <= hi.z;
}

// Locates points against the ship's spaces. The bounding boxes are the whole
// point: a probe near the shell is inside one compartment's box and outside the
// other fifteen, so a query costs one solid-angle sum over a few hundred
// triangles rather than sixteen.
class SpaceLocator {
public:
    explicit SpaceLocator(const Ship& ship) : ship_(&ship) {
        boundingBox(ship.hull, hullLo_, hullHi_);
        lo_.resize(ship.compartments.size());
        hi_.resize(ship.compartments.size());
        for (std::size_t i = 0; i < ship.compartments.size(); ++i)
            boundingBox(ship.compartments[i].mesh, lo_[i], hi_[i]);
    }

    int at(const Vec3& p) const {
        // Compartments first: they are the specific answer, and a point inside
        // one is inside the hull by construction.
        for (std::size_t i = 0; i < ship_->compartments.size(); ++i)
            if (inside(i, p)) return static_cast<int>(i);
        if (!insideBox(p, hullLo_, hullHi_)) return kSea;
        return std::abs(meshWindingNumber(ship_->hull, p)) > 0.5 ? kEnclosedVoid : kSea;
    }

    // A *second* compartment holding this point, or kSea when there is only one.
    //
    // Overlapping subdivision is not hypothetical: the reference ferry's forward
    // wing tanks lie entirely inside her forward holds, so a tear into one of
    // them is silently attributed to the other, and `Ship::validate()` does not
    // catch it because the subdivision still totals less than the hull. Which
    // compartment wins is arbitrary; saying that it *was* arbitrary is not.
    //
    // The claim has to survive a millimetre's displacement in all six directions
    // before it is made, because two compartments that merely *abut* also both
    // contain a point on the bulkhead between them -- meshWindingNumber() counts
    // a point lying in the plane of a face as inside, so a probe that lands on a
    // shared face reads as being in both. A region with volume in it survives the
    // displacement; a face does not. The check therefore fails towards missing a
    // real overlap rather than towards inventing one, which is the right way
    // round for a complaint about somebody's ship file.
    int overlapAt(const Vec3& p) const {
        int first = kSea, second = kSea;
        for (std::size_t i = 0; i < ship_->compartments.size(); ++i) {
            if (!inside(i, p)) continue;
            if (first == kSea) {
                first = static_cast<int>(i);
                continue;
            }
            second = static_cast<int>(i);
            break;
        }
        if (second == kSea) return kSea;
        for (int axis = 0; axis < 3; ++axis)
            for (const double sign : {-1.0, 1.0}) {
                Vec3 q = p;
                q[axis] += sign * 1e-3;
                if (!inside(static_cast<std::size_t>(first), q) ||
                    !inside(static_cast<std::size_t>(second), q))
                    return kSea;
            }
        return second;
    }

private:
    bool inside(std::size_t compartment, const Vec3& p) const {
        return insideBox(p, lo_[compartment], hi_[compartment]) &&
               std::abs(meshWindingNumber(ship_->compartments[compartment].mesh, p)) > 0.5;
    }

    const Ship* ship_;
    Vec3 hullLo_{}, hullHi_{};
    std::vector<Vec3> lo_, hi_;
};

// --- Welding panel corners ----------------------------------------------------
//
// Two panels are adjacent when they share an edge, and they share an edge when
// two of their corners are the same corner. Generated panels do hold their shared
// corners bit-identically -- `Section::at()` is a pure function of the girth
// fraction, so both panels either side of a seam evaluate it to the same bits --
// but the mirrored pair either side of the centreline meets at y = 0 and y = -0.0,
// and a hand-built mesh need not be exact at all. So corners are welded on a
// tolerance rather than compared.

struct Cell {
    std::int64_t x = 0, y = 0, z = 0;
    bool operator==(const Cell& o) const = default;
};

struct CellHash {
    std::size_t operator()(const Cell& c) const noexcept {
        std::uint64_t h = 1469598103934665603ull;
        for (std::int64_t v : {c.x, c.y, c.z}) {
            h ^= static_cast<std::uint64_t>(v);
            h *= 1099511628211ull;
        }
        return static_cast<std::size_t>(h);
    }
};

// Cells are one tolerance across, so two points within tolerance of each other
// are always in the same cell or in one of its 26 neighbours. Searching all 27 is
// what makes the weld independent of where the grid happens to fall -- bucketing
// by a single cell would split a pair that straddles a cell boundary, which is
// the failure that only shows up on some meshes and looks like a merging bug.
class Welder {
public:
    explicit Welder(double tolerance) : eps_(tolerance > 0 ? tolerance : 1e-12) {}

    int add(const Vec3& p) {
        const Cell home = cellOf(p);
        for (std::int64_t dx = -1; dx <= 1; ++dx)
            for (std::int64_t dy = -1; dy <= 1; ++dy)
                for (std::int64_t dz = -1; dz <= 1; ++dz) {
                    const auto it = cells_.find(Cell{home.x + dx, home.y + dy, home.z + dz});
                    if (it == cells_.end()) continue;
                    for (int id : it->second)
                        if (length2(points_[static_cast<std::size_t>(id)] - p) <= eps_ * eps_)
                            return id;
                }
        const int id = static_cast<int>(points_.size());
        points_.push_back(p);
        cells_[home].push_back(id);
        return id;
    }

private:
    Cell cellOf(const Vec3& p) const {
        return {static_cast<std::int64_t>(std::floor(p.x / eps_)),
                static_cast<std::int64_t>(std::floor(p.y / eps_)),
                static_cast<std::int64_t>(std::floor(p.z / eps_))};
    }

    double eps_;
    std::vector<Vec3> points_;
    std::unordered_map<Cell, std::vector<int>, CellHash> cells_;
};

// --- Connected components -----------------------------------------------------

class DisjointSet {
public:
    explicit DisjointSet(std::size_t n) : parent_(n) {
        for (std::size_t i = 0; i < n; ++i) parent_[i] = i;
    }
    std::size_t find(std::size_t i) {
        while (parent_[i] != i) {
            parent_[i] = parent_[parent_[i]];
            i = parent_[i];
        }
        return i;
    }
    void join(std::size_t a, std::size_t b) {
        const std::size_t ra = find(a), rb = find(b);
        if (ra != rb) parent_[ra] = rb;
    }

private:
    std::vector<std::size_t> parent_;
};

// The two spaces a panel separates, found by stepping out along its own normal
// until the sides disagree. See BreachParams::probeDistance for why this marches
// rather than probing once: a flat panel chording across a curved shell can have
// its centroid 0.15 m inside the hull on this very ferry, and a fixed probe short
// enough to stay inside a double bottom would report the shell plating as
// separating the double bottom from itself.
struct PanelSpaces {
    int a = kSea, b = kSea;
    Vec3 pointA{}, pointB{};  // where the answer was finally read, for a second opinion
};

PanelSpaces spacesEitherSide(const SpaceLocator& locator, const PlatePanel& panel,
                             const Vec3& normal, double step) {
    const Vec3 centre = panel.centroid();
    // As far as the panel itself reaches and no further: the surface a panel
    // stands for passes somewhere within the panel, so no point of it is further
    // from the centroid than its furthest corner. Stated as a corner distance
    // rather than as half a diagonal because a warped quad has two diagonals and
    // the choice between them would be arbitrary.
    double reach = 0;
    for (const Vec3& corner : panel.corner) reach = std::max(reach, length(corner - centre));
    PanelSpaces spaces;
    for (double d = step > 0 ? step : 1e-6;; d = std::min(2.0 * d, reach)) {
        spaces.pointA = centre + normal * d;
        spaces.pointB = centre - normal * d;
        spaces.a = locator.at(spaces.pointA);
        spaces.b = locator.at(spaces.pointB);
        if (spaces.a != spaces.b || d >= reach) break;
    }
    return spaces;
}

std::string spaceName(const Ship& ship, int space) {
    if (space == kSea) return "sea";
    if (space == kEnclosedVoid) return "unmodelled_void";
    const std::size_t i = static_cast<std::size_t>(space);
    if (i >= ship.compartments.size()) return "compartment_" + std::to_string(space);
    const std::string& name = ship.compartments[i].name;
    return name.empty() ? "compartment_" + std::to_string(space) : name;
}

}  // namespace

// ---------------------------------------------------------------------------
// Point location
// ---------------------------------------------------------------------------

double meshWindingNumber(const TriMesh& mesh, const Vec3& point) {
    // Van Oosterom & Strackee: the solid angle a triangle subtends at the origin
    // is 2*atan2(a.(b x c), |a||b||c| + (a.b)|c| + (b.c)|a| + (c.a)|b|). Summed
    // over a closed surface this is 4*pi inside and 0 outside, with no case
    // analysis anywhere -- the atan2 carries the sign and the quadrant.
    double solidAngle = 0;
    for (const Tri& t : mesh.tris) {
        const Vec3 a = mesh.verts[t.a] - point;
        const Vec3 b = mesh.verts[t.b] - point;
        const Vec3 c = mesh.verts[t.c] - point;
        const double la = length(a), lb = length(b), lc = length(c);
        const double numerator = dot(a, cross(b, c));
        const double denominator =
            la * lb * lc + dot(a, b) * lc + dot(b, c) * la + dot(c, a) * lb;
        // atan2(0, 0) is 0 rather than a NaN, which is the answer wanted when the
        // point sits on a vertex: that triangle subtends nothing measurable and
        // the query was ill-posed anyway.
        solidAngle += 2.0 * std::atan2(numerator, denominator);
    }
    return solidAngle / (4.0 * kPi);
}

int spaceAt(const Ship& ship, const Vec3& bodyPoint) {
    return SpaceLocator(ship).at(bodyPoint);
}

// ---------------------------------------------------------------------------
// Failed panels -> openings
// ---------------------------------------------------------------------------

std::vector<Opening> BreachSet::openings() const {
    std::vector<Opening> out;
    out.reserve(breaches.size());
    for (const Breach& b : breaches) out.push_back(b.opening);
    return out;
}

double BreachSet::totalArea() const {
    double total = 0;
    for (const Breach& b : breaches) total += b.opening.area;
    return total;
}

BreachSet breachesFromFailedPanels(const Ship& ship, const StructuralMesh& mesh,
                                   const std::vector<int>& failedPanels,
                                   const BreachParams& params) {
    BreachSet out;

    // --- 1. Clean the input ---------------------------------------------------
    //
    // Sorted and de-duplicated, because a panel listed twice would otherwise
    // contribute its area twice -- an opening that is quietly too big is exactly
    // the failure this whole file exists to get right.
    std::vector<int> panels;
    panels.reserve(failedPanels.size());
    for (int index : failedPanels) {
        if (index < 0 || static_cast<std::size_t>(index) >= mesh.panels.size()) {
            out.problems.push_back("failed panel " + std::to_string(index) +
                                   " is not in a mesh of " +
                                   std::to_string(mesh.panels.size()) + " panels");
            continue;
        }
        panels.push_back(index);
    }
    std::sort(panels.begin(), panels.end());
    for (std::size_t i = 1; i < panels.size(); ++i)
        if (panels[i] == panels[i - 1])
            out.problems.push_back("failed panel " + std::to_string(panels[i]) +
                                   " was listed more than once and was counted once");
    panels.erase(std::unique(panels.begin(), panels.end()), panels.end());
    if (panels.empty()) return out;

    // --- 2. Ask each panel which two spaces it separates -----------------------

    // Panels that produce no opening are tallied by reason rather than listed one
    // by one: ramming the ferry fails a patch of side plating that faces an
    // unmodelled void, and forty identical lines would bury the one that matters.
    std::map<std::string, std::pair<int, int>> skipped;  // reason -> (count, first panel)
    const auto skip = [&](const std::string& reason, int panel) {
        auto& tally = skipped[reason];
        if (tally.first == 0) tally.second = panel;
        ++tally.first;
    };
    // Remarks about the ship definition rather than about a panel, so they are
    // said once however many panels notice them.
    std::set<std::string> notes;
    const auto note = [&](const std::string& text) { notes.insert(text); };

    const SpaceLocator locator(ship);
    std::map<std::pair<int, int>, std::vector<std::size_t>> classes;
    std::vector<std::pair<int, int>> connects(panels.size(), {kSea, kSea});

    for (std::size_t i = 0; i < panels.size(); ++i) {
        const PlatePanel& panel = mesh.panels[static_cast<std::size_t>(panels[i])];
        const Vec3 normal = panel.normal();
        if (length2(normal) < 0.5) {
            // normalize() hands back a zero vector rather than a NaN, so a
            // degenerate quad would otherwise probe the same point twice and read
            // as a space opening onto itself -- true, but for the wrong reason.
            skip("are degenerate, so there is no plane to probe either side of", panels[i]);
            continue;
        }
        const PanelSpaces probe = spacesEitherSide(locator, panel, normal, params.probeDistance);
        for (const auto& [space, point] : {std::pair{probe.a, probe.pointA},
                                           std::pair{probe.b, probe.pointB}}) {
            if (space < 0) continue;
            const int overlap = locator.overlapAt(point);
            if (overlap != kSea)
                note("the subdivision puts " + spaceName(ship, space) + " and " +
                     spaceName(ship, overlap) + " in the same place, so a tear between them is "
                     "attributed to whichever was declared first");
        }
        int a = probe.a, b = probe.b;
        if (a > b) std::swap(a, b);  // kEnclosedVoid < kSea < any compartment

        if (a == b) {
            skip("separate " + spaceName(ship, a) + " from itself", panels[i]);
            continue;
        }
        if (a == kEnclosedVoid) {
            skip("open " + spaceName(ship, b) + " onto a part of the hull that no compartment "
                 "describes, so nothing floods through them", panels[i]);
            continue;
        }
        connects[i] = {a, b};
        classes[{a, b}].push_back(i);
    }

    // --- 3. Merge failures that share an edge ---------------------------------

    Welder welder(params.weldEpsilon);
    std::vector<std::array<int, 4>> corners(panels.size());
    for (std::size_t i = 0; i < panels.size(); ++i) {
        if (connects[i].first == connects[i].second) continue;  // skipped above
        const PlatePanel& panel = mesh.panels[static_cast<std::size_t>(panels[i])];
        for (int k = 0; k < 4; ++k) corners[i][static_cast<std::size_t>(k)] = welder.add(panel.corner[k]);
    }

    // `classes` is ordered by the pair of spaces, so the openings come out in a
    // fixed order however the caller listed the failures.
    for (const auto& [pair, members] : classes) {
        DisjointSet sets(members.size());
        std::map<std::pair<int, int>, std::vector<std::size_t>> edges;
        for (std::size_t m = 0; m < members.size(); ++m) {
            const std::array<int, 4>& corner = corners[members[m]];
            for (int k = 0; k < 4; ++k) {
                const int u = corner[static_cast<std::size_t>(k)];
                const int v = corner[static_cast<std::size_t>((k + 1) % 4)];
                if (u == v) continue;  // a collapsed side of a triangular quad
                edges[{std::min(u, v), std::max(u, v)}].push_back(m);
            }
        }
        // Sharing an *edge*, not merely two corners: two panels that touch only
        // at a corner appear in no common edge bucket and stay two openings.
        for (const auto& [edge, sharing] : edges) {
            (void)edge;
            for (std::size_t j = 1; j < sharing.size(); ++j) sets.join(sharing[0], sharing[j]);
        }

        // Components in order of their lowest member, which -- since `members` is
        // ascending -- is order of first appearance.
        std::map<std::size_t, std::size_t> componentOf;
        std::vector<std::vector<std::size_t>> components;
        for (std::size_t m = 0; m < members.size(); ++m) {
            const std::size_t root = sets.find(m);
            const auto [it, fresh] = componentOf.try_emplace(root, components.size());
            if (fresh) components.emplace_back();
            components[it->second].push_back(m);
        }

        for (std::size_t k = 0; k < components.size(); ++k) {
            double area = 0;
            Vec3 moment{};
            Breach breach;
            breach.panels.reserve(components[k].size());
            for (std::size_t m : components[k]) {
                const int index = panels[members[m]];
                const PlatePanel& panel = mesh.panels[static_cast<std::size_t>(index)];
                // Summed in ascending panel order and by nothing else: the
                // opening's area *is* the failed area, with no coefficient on it,
                // and the test asserts that to the last bit.
                const double panelArea = panel.area();
                area += panelArea;
                moment += panel.centroid() * panelArea;
                breach.panels.push_back(index);
            }

            Opening& opening = breach.opening;
            opening.name = params.namePrefix + "_" + spaceName(ship, pair.first) + "_" +
                           spaceName(ship, pair.second) + "_" + std::to_string(k + 1);
            opening.a = pair.first;
            opening.b = pair.second;
            // The failed region's own centroid, not any one panel's. The head at
            // the opening -- and so the flow through it -- is decided by this.
            opening.pos = area > 0 ? moment / area
                                   : mesh.panels[static_cast<std::size_t>(breach.panels[0])].centroid();
            opening.area = area;
            opening.dischargeCoeff = params.dischargeCoeff;
            opening.open = true;
            opening.kind = OpeningKind::Breach;
            out.breaches.push_back(std::move(breach));
        }
    }

    for (const std::string& text : notes) out.problems.push_back(text);
    for (const auto& [reason, tally] : skipped)
        out.problems.push_back(std::to_string(tally.first) + " failed panel(s) " + reason +
                               " (the first is panel " + std::to_string(tally.second) + ")");
    return out;
}

std::size_t applyBreaches(Ship& ship, const BreachSet& set) {
    ship.openings.reserve(ship.openings.size() + set.breaches.size());
    for (const Breach& b : set.breaches) ship.openings.push_back(b.opening);
    return set.breaches.size();
}

}  // namespace sim
