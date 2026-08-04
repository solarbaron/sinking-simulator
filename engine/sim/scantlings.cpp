// SPDX-License-Identifier: MIT
#include "scantlings.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <utility>

namespace sim {
namespace {

// --- Section properties -------------------------------------------------------

// One rectangle of a built-up section, in the profile's own axes: `n` runs along
// the web away from the plate, `m` across it.
struct Rect {
    double area = 0;
    double n = 0;      // centroid along the web
    double m = 0;      // centroid across the web
    double iOwnN = 0;  // second moment about its own centroid, integrating n^2
    double iOwnM = 0;  // ... integrating m^2
};

Rect rectangle(double alongN, double alongM, double n, double m) {
    Rect r;
    r.area = alongN * alongM;
    r.n = n;
    r.m = m;
    r.iOwnN = alongM * alongN * alongN * alongN / 12.0;
    r.iOwnM = alongN * alongM * alongM * alongM / 12.0;
    return r;
}

// --- Ray casting against the hull ---------------------------------------------
//
// The structure samples the hull rather than sharing its vertices (see the
// header), and the sample it needs is "how far out is the shell at this height,
// on this station". A ray cast answers that for any closed mesh, whatever built
// it -- a station grid, an offset table or an importer -- which is the point.

// Moller-Trumbore for the one ray direction this file ever uses: +y, from the
// centreline. Returns the distance to the intersection, or -1 for a miss.
//
// Triangles lying *in* the ray's plane -- the transom and stem caps, when a
// station lands exactly on them -- have a zero determinant and are rejected,
// which is the behaviour wanted: they are not shell.
double intersectY(const Vec3& v0, const Vec3& v1, const Vec3& v2, double x, double z) {
    constexpr double kEdge = 1e-9;
    const Vec3 e1 = v1 - v0;
    const Vec3 e2 = v2 - v0;
    // cross(dir, e2) with dir = (0, 1, 0).
    const Vec3 p{e2.z, 0.0, -e2.x};
    const double det = dot(e1, p);
    if (std::abs(det) < 1e-18) return -1.0;
    const double inv = 1.0 / det;
    const Vec3 t{x - v0.x, -v0.y, z - v0.z};
    const double u = dot(t, p) * inv;
    if (u < -kEdge || u > 1.0 + kEdge) return -1.0;
    const Vec3 q = cross(t, e1);
    const double v = q.y * inv;  // dot(dir, q) with dir = (0, 1, 0)
    if (v < -kEdge || u + v > 1.0 + kEdge) return -1.0;
    const double hit = dot(e2, q) * inv;
    return hit >= 0.0 ? hit : -1.0;
}

class HullProbe {
public:
    explicit HullProbe(const TriMesh& mesh) : mesh_(&mesh) {
        if (mesh.verts.empty()) return;
        lo_ = hi_ = mesh.verts[0];
        for (const Vec3& v : mesh.verts) {
            lo_ = {std::min(lo_.x, v.x), std::min(lo_.y, v.y), std::min(lo_.z, v.z)};
            hi_ = {std::max(hi_.x, v.x), std::max(hi_.y, v.y), std::max(hi_.z, v.z)};
        }
        span_.reserve(mesh.tris.size());
        for (const Tri& t : mesh.tris) {
            const double a = mesh.verts[t.a].x, b = mesh.verts[t.b].x, c = mesh.verts[t.c].x;
            span_.push_back({std::min(a, std::min(b, c)), std::max(a, std::max(b, c))});
        }
    }

    const Vec3& lo() const { return lo_; }
    const Vec3& hi() const { return hi_; }
    bool empty() const { return mesh_->tris.empty(); }

    // A station's worth of triangles, gathered once so that the hundred or so
    // height samples that follow test a fiftieth of the mesh each. Without this
    // the sampling is O(stations * samples * triangles) and shows up in the test
    // suite's wall time.
    std::vector<std::uint32_t> station(double x) const {
        std::vector<std::uint32_t> out;
        for (std::size_t i = 0; i < span_.size(); ++i)
            if (span_[i].first <= x && x <= span_[i].second)
                out.push_back(static_cast<std::uint32_t>(i));
        return out;
    }

    // Outermost shell offset at (x, z), or 0 where the ray misses entirely.
    // Outermost rather than nearest, so internal geometry -- if a hull ever
    // carries any -- cannot be mistaken for the shell.
    double halfBreadth(const std::vector<std::uint32_t>& candidates, double x, double z) const {
        double best = -1.0;
        for (std::uint32_t i : candidates) {
            const Tri& t = mesh_->tris[i];
            const double hit =
                intersectY(mesh_->verts[t.a], mesh_->verts[t.b], mesh_->verts[t.c], x, z);
            if (hit > best) best = hit;
        }
        return best > 0.0 ? best : 0.0;
    }

private:
    const TriMesh* mesh_;
    Vec3 lo_{}, hi_{};
    std::vector<std::pair<double, double>> span_;
};

// The girth curve of one transverse section: from the centreline at the keel,
// outboard across the flat of bottom, round the bilge and up to the sheer. Port
// side only -- the hull mesh format is a mirrored half-breadth table, so the
// starboard side is the mirror by construction rather than by assumption.
struct Section {
    double x = 0;
    std::vector<Vec3> point;
    std::vector<double> arc;  // cumulative length; arc.back() is the girth

    double girth() const { return arc.empty() ? 0.0 : arc.back(); }

    // Point at a fraction of the girth, and the inward normal of the plating
    // there. `normal` is what a stiffener's web rises along and what the
    // eccentricity of an offset beam is measured on.
    Vec3 at(double fraction, Vec3* normal = nullptr) const {
        if (point.empty()) return {};
        if (point.size() == 1) {
            if (normal) *normal = {0, 0, 1};
            return point[0];
        }
        const double target = std::clamp(fraction, 0.0, 1.0) * girth();
        std::size_t i = static_cast<std::size_t>(
            std::lower_bound(arc.begin(), arc.end(), target) - arc.begin());
        if (i == 0) i = 1;
        if (i >= point.size()) i = point.size() - 1;
        const double lo = arc[i - 1], hi = arc[i];
        const double u = hi > lo ? (target - lo) / (hi - lo) : 0.0;
        if (normal) {
            const Vec3 tangent = normalize(point[i] - point[i - 1]);
            *normal = normalize(cross(Vec3{1, 0, 0}, tangent));
        }
        return point[i - 1] + (point[i] - point[i - 1]) * u;
    }
};

Section sampleSection(const HullProbe& probe, double x, int samples) {
    Section s;
    s.x = x;
    if (probe.empty()) return s;

    const double zSpan = probe.hi().z - probe.lo().z;
    const double xSpan = probe.hi().x - probe.lo().x;
    if (zSpan <= 0.0) return s;

    // Rays are cast a whisker inside the bounding box in both x and z. On the
    // boundary they are coplanar with the transom, stem and deck caps, and a
    // coplanar triangle is a division by nearly zero rather than an answer. The
    // sample is reported at the station that was asked for, so panels still span
    // the full length; only the half-breadth comes from the nudged ray, and it is
    // wrong by dy/dx times a ten-millionth of the ship.
    const double nudgeX = 1e-7 * std::max(xSpan, 1e-9);
    const double nudgeZ = 1e-7 * zSpan;
    const double xq = std::clamp(x, probe.lo().x + nudgeX, probe.hi().x - nudgeX);
    const std::vector<std::uint32_t> candidates = probe.station(xq);

    const int n = std::max(3, samples);
    const double z0 = probe.lo().z + nudgeZ;
    const double z1 = probe.hi().z - nudgeZ;
    std::vector<Vec3> raw;
    raw.reserve(static_cast<std::size_t>(n));
    for (int k = 0; k < n; ++k) {
        const double z = z0 + (z1 - z0) * k / (n - 1);
        raw.push_back({x, probe.halfBreadth(candidates, xq, z), z});
    }

    // Trim heights where the section does not exist -- below a rising floor, or
    // above the sheer on a hull whose deck is not flat -- and close the outline
    // onto the centreline at the lowest height that does.
    const double tiny = 1e-9 * std::max(probe.hi().y - probe.lo().y, 1e-9);
    std::size_t first = 0;
    while (first < raw.size() && raw[first].y <= tiny) ++first;
    std::size_t last = raw.size();
    while (last > first && raw[last - 1].y <= tiny) --last;
    if (first >= last) return s;

    s.point.reserve(last - first + 1);
    s.point.push_back({x, 0.0, raw[first].z});
    for (std::size_t i = first; i < last; ++i) s.point.push_back(raw[i]);

    s.arc.resize(s.point.size());
    s.arc[0] = 0.0;
    for (std::size_t i = 1; i < s.point.size(); ++i)
        s.arc[i] = s.arc[i - 1] + length(s.point[i] - s.point[i - 1]);
    return s;
}

Vec3 mirrored(const Vec3& v) { return {v.x, -v.y, v.z}; }

void addPanel(StructuralMesh& mesh, const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d,
              double thickness, int material, PanelRole role, int source) {
    PlatePanel p;
    p.corner[0] = a;
    p.corner[1] = b;
    p.corner[2] = c;
    p.corner[3] = d;
    p.thickness = thickness;
    p.material = material;
    p.role = role;
    p.source = source;
    if (p.area() <= 1e-12) return;  // a band the hull has closed to nothing
    mesh.panels.push_back(p);
}

void addMember(StructuralMesh& mesh, const Vec3& a, const Vec3& b, const Vec3& rise,
               const StiffenerProfile& profile, double plateThickness, int material,
               MemberRole role) {
    StructuralMember m;
    m.a = a;
    m.b = b;
    m.rise = normalize(rise);
    m.profile = profile;
    m.attachedPlateThickness = plateThickness;
    m.material = material;
    m.role = role;
    if (m.length() <= 1e-9) return;
    mesh.members.push_back(m);
}

}  // namespace

// --- Materials and profiles ---------------------------------------------------

StructuralMaterial mildSteel() {
    StructuralMaterial m;
    m.name = "steel_a";
    m.yieldStrength = 235.0e6;
    return m;
}

StructuralMaterial ah36Steel() { return StructuralMaterial{}; }

StiffenerProfile flatBar(double webHeight, double webThickness) {
    StiffenerProfile p;
    p.kind = ProfileKind::FlatBar;
    p.webHeight = webHeight;
    p.webThickness = webThickness;
    return p;
}

StiffenerProfile tee(double webHeight, double webThickness, double flangeWidth,
                     double flangeThickness) {
    StiffenerProfile p;
    p.kind = ProfileKind::Tee;
    p.webHeight = webHeight;
    p.webThickness = webThickness;
    p.flangeWidth = flangeWidth;
    p.flangeThickness = flangeThickness;
    return p;
}

StiffenerProfile angle(double webHeight, double webThickness, double flangeWidth,
                       double flangeThickness) {
    StiffenerProfile p = tee(webHeight, webThickness, flangeWidth, flangeThickness);
    p.kind = ProfileKind::Angle;
    return p;
}

ProfileSection profileSection(const StiffenerProfile& profile) {
    ProfileSection s;
    const double hw = std::max(0.0, profile.webHeight);
    const double tw = std::max(0.0, profile.webThickness);
    const bool flanged = profile.kind != ProfileKind::FlatBar && profile.flangeWidth > 0 &&
                         profile.flangeThickness > 0;
    const double bf = flanged ? profile.flangeWidth : 0.0;
    const double tf = flanged ? profile.flangeThickness : 0.0;

    // An angle's flange hangs to one side of the web; a tee's straddles it. The
    // difference is invisible to bending in the plane of the web and is the whole
    // story about the axis at right angles to it, which is what decides tripping.
    const double flangeOffset = profile.kind == ProfileKind::Angle ? 0.5 * bf : 0.0;

    Rect parts[2];
    int count = 0;
    if (hw > 0 && tw > 0) parts[count++] = rectangle(hw, tw, 0.5 * hw, 0.0);
    if (flanged) parts[count++] = rectangle(tf, bf, hw + 0.5 * tf, flangeOffset);
    if (count == 0) return s;

    double area = 0, momentN = 0, momentM = 0;
    for (int i = 0; i < count; ++i) {
        area += parts[i].area;
        momentN += parts[i].area * parts[i].n;
        momentM += parts[i].area * parts[i].m;
    }
    if (area <= 0) return s;
    const double cn = momentN / area;
    const double cm = momentM / area;

    double iN = 0, iM = 0;
    for (int i = 0; i < count; ++i) {
        iN += parts[i].iOwnN + parts[i].area * (parts[i].n - cn) * (parts[i].n - cn);
        iM += parts[i].iOwnM + parts[i].area * (parts[i].m - cm) * (parts[i].m - cm);
    }

    s.area = area;
    s.centroid = cn;
    s.height = hw + tf;
    s.secondMoment = iN;
    s.secondMomentWeak = iM;
    // St Venant torsion of an open section built from thin rectangles.
    s.torsionConstant = (hw * tw * tw * tw + bf * tf * tf * tf) / 3.0;
    return s;
}

StiffenedSection stiffenedSection(const StiffenerProfile& profile, double plateThickness,
                                  double plateWidth) {
    StiffenedSection s;
    const double tp = std::max(0.0, plateThickness);
    const double bp = std::max(0.0, plateWidth);
    const ProfileSection ps = profileSection(profile);

    const double plateArea = bp * tp;
    const double plateOwn = bp * tp * tp * tp / 12.0;
    // The profile is welded to the plate's face, which is half a thickness from
    // the mid-surface everything here is measured against.
    const double profileCentroid = 0.5 * tp + ps.centroid;

    s.area = plateArea + ps.area;
    if (s.area <= 0) return s;
    s.neutralAxis = ps.area * profileCentroid / s.area;
    s.secondMoment = plateOwn + plateArea * s.neutralAxis * s.neutralAxis + ps.secondMoment +
                     ps.area * (profileCentroid - s.neutralAxis) * (profileCentroid - s.neutralAxis);
    s.height = 0.5 * tp + ps.height;

    const double toPlate = s.neutralAxis + 0.5 * tp;
    const double toStiffener = s.height - s.neutralAxis;
    if (toPlate > 1e-12) s.modulusPlate = s.secondMoment / toPlate;
    if (toStiffener > 1e-12) s.modulusStiffener = s.secondMoment / toStiffener;
    return s;
}

double smearedThickness(double plateThickness, const StiffenerProfile& profile, double spacing) {
    if (spacing <= 1e-9) return plateThickness;
    return plateThickness + profileSection(profile).area / spacing;
}

// --- Mesh accessors -----------------------------------------------------------

double PlatePanel::area() const {
    // Two triangles rather than a shoelace formula, so a warped quad -- which
    // every panel on a curved shell is -- still gets its true area.
    return 0.5 * length(cross(corner[1] - corner[0], corner[2] - corner[0])) +
           0.5 * length(cross(corner[2] - corner[0], corner[3] - corner[0]));
}

Vec3 PlatePanel::centroid() const {
    const double a1 = 0.5 * length(cross(corner[1] - corner[0], corner[2] - corner[0]));
    const double a2 = 0.5 * length(cross(corner[2] - corner[0], corner[3] - corner[0]));
    const Vec3 c1 = (corner[0] + corner[1] + corner[2]) / 3.0;
    const Vec3 c2 = (corner[0] + corner[2] + corner[3]) / 3.0;
    const double total = a1 + a2;
    if (total <= 1e-15) return corner[0];
    return (c1 * a1 + c2 * a2) / total;
}

Vec3 PlatePanel::normal() const {
    return normalize(cross(corner[2] - corner[0], corner[3] - corner[1]));
}

double StructuralMember::length() const { return sim::length(b - a); }

double StructuralMesh::plateArea() const {
    double total = 0;
    for (const PlatePanel& p : panels) total += p.area();
    return total;
}

double StructuralMesh::plateMass() const {
    double total = 0;
    for (const PlatePanel& p : panels) {
        const std::size_t m = static_cast<std::size_t>(p.material);
        const double density = m < materials.size() ? materials[m].density : 7850.0;
        total += p.area() * p.thickness * density;
    }
    return total;
}

double StructuralMesh::memberLength() const {
    double total = 0;
    for (const StructuralMember& s : members) total += s.length();
    return total;
}

double StructuralMesh::memberMass() const {
    double total = 0;
    for (const StructuralMember& s : members) {
        const std::size_t m = static_cast<std::size_t>(s.material);
        const double density = m < materials.size() ? materials[m].density : 7850.0;
        total += s.length() * profileSection(s.profile).area * density;
    }
    return total;
}

std::size_t StructuralMesh::panelCount(PanelRole role) const {
    std::size_t n = 0;
    for (const PlatePanel& p : panels)
        if (p.role == role) ++n;
    return n;
}

std::size_t StructuralMesh::memberCount(MemberRole role) const {
    std::size_t n = 0;
    for (const StructuralMember& s : members)
        if (s.role == role) ++n;
    return n;
}

// --- Generation ---------------------------------------------------------------

namespace {

// Where the strake seams and the longitudinals sit, as fractions of the section
// girth. **Global, not per station**: a longitudinal has to be one continuous
// member from end to end, so its girth *fraction* is what stays constant while
// the girth it lives on shrinks towards the ends. The consequence -- that the
// delivered spacing at the ends is smaller than amidships -- is what a real ship
// does too, up to the point where longitudinals are dropped instead, which this
// cannot express.
struct GirthLayout {
    std::vector<double> fraction;  // ascending, 0 .. 1
};

const ShellRegion* regionAt(const Scantlings& s, double x, double girthFraction) {
    const ShellRegion* found = nullptr;
    // Last match wins, so a thickened patch declared after the strake it sits in
    // overrides it rather than fighting it.
    for (const ShellRegion& r : s.shell)
        if (x >= r.xFrom && x <= r.xTo && girthFraction >= r.girthFrom - 1e-12 &&
            girthFraction <= r.girthTo + 1e-12)
            found = &r;
    return found;
}

GirthLayout buildGirthLayout(const Scantlings& s, double midshipGirth) {
    std::vector<double> breaks{0.0, 1.0};
    for (const ShellRegion& r : s.shell) {
        breaks.push_back(std::clamp(r.girthFrom, 0.0, 1.0));
        breaks.push_back(std::clamp(r.girthTo, 0.0, 1.0));
    }
    std::sort(breaks.begin(), breaks.end());
    breaks.erase(std::unique(breaks.begin(), breaks.end(),
                             [](double a, double b) { return std::abs(a - b) < 1e-9; }),
                 breaks.end());

    GirthLayout layout;
    layout.fraction.push_back(0.0);
    for (std::size_t k = 0; k + 1 < breaks.size(); ++k) {
        const double lo = breaks[k], hi = breaks[k + 1];
        const double mid = 0.5 * (lo + hi);

        // The tightest spacing any region asks for over this strake, because a
        // spacing is a maximum: exceeding it somewhere is a scantling failure,
        // undershooting it is only steel.
        double spacing = 0.0;
        for (const ShellRegion& r : s.shell) {
            if (mid < r.girthFrom - 1e-12 || mid > r.girthTo + 1e-12) continue;
            const double want = r.longitudinalSpacing > 0 ? r.longitudinalSpacing
                                                          : s.longitudinalSpacing;
            if (want > 0 && (spacing <= 0 || want < spacing)) spacing = want;
        }
        if (spacing <= 0) spacing = s.longitudinalSpacing;

        int bands = 1;
        if (spacing > 0 && midshipGirth > 0)
            bands = std::max(1, static_cast<int>(std::lround((hi - lo) * midshipGirth / spacing)));
        for (int i = 1; i <= bands; ++i) layout.fraction.push_back(lo + (hi - lo) * i / bands);
    }
    layout.fraction.back() = 1.0;
    return layout;
}

// Bay boundaries for a member that runs between two arbitrary stations: the frame
// stations it spans, closed at each end by its own. A deck or a girder that stops
// between frames is ordinary; forcing it onto a frame station would be a lie
// about where it stops.
std::vector<double> baysBetween(const std::vector<double>& frames, double xFrom, double xTo) {
    std::vector<double> edges{xFrom};
    for (double x : frames)
        if (x > xFrom + 1e-9 && x < xTo - 1e-9) edges.push_back(x);
    edges.push_back(xTo);
    return edges;
}

// Frame stations: the hull length divided into a whole number of bays. Laying
// them out from one end at the requested pitch instead leaves a ragged last bay,
// and a bay a millimetre long is a degenerate panel rather than a ship.
std::vector<double> frameStations(double xLo, double xHi, double spacing, int* bays) {
    const double length = xHi - xLo;
    int n = 1;
    if (spacing > 0 && length > 0) n = std::max(1, static_cast<int>(std::lround(length / spacing)));
    *bays = n;
    std::vector<double> out;
    out.reserve(static_cast<std::size_t>(n) + 1);
    for (int i = 0; i <= n; ++i) out.push_back(xLo + length * i / n);
    return out;
}


// The generation passes share a good deal of state -- the hull being sampled, the
// mesh being filled, and the running account of everything the description asked
// for that the hull could not give -- so they are members rather than a parameter
// list longer than most of them.
struct Builder {
    Builder(const TriMesh& hull, const Scantlings& description)
        : probe(hull), scantlings(description) {
        mesh.materials = description.materials;
        if (mesh.materials.empty()) mesh.materials.push_back(ah36Steel());
        xNudge = 1e-7 * (probe.hi().x - probe.lo().x);
        zNudge = 1e-7 * (probe.hi().z - probe.lo().z);
    }

    HullProbe probe;
    const Scantlings& scantlings;
    StructuralMesh mesh;
    std::vector<std::string> problems;
    std::vector<Section> sections;   // one per frame station
    GirthLayout layout;
    double xNudge = 0, zNudge = 0;

    void report(const std::string& what) { problems.push_back(what); }

    // A material index that is always in range, so a description naming one that
    // does not exist still builds -- validateScantlings() has already said so, and
    // refusing to produce geometry as well would hide the rest of the report.
    int material(int index) const {
        return std::clamp(index, 0, static_cast<int>(mesh.materials.size()) - 1);
    }

    // Half-breadth anywhere on the hull, with the station's triangle list cached
    // for the run of heights that follows. Every pass below walks heights inside a
    // station, so one cached entry turns an O(triangles) gather per query into one
    // per station -- measurably, the difference between the whole ferry structure
    // taking three milliseconds and taking a second.
    double breadthAt(double x, double z) {
        const double xq = std::clamp(x, probe.lo().x + xNudge, probe.hi().x - xNudge);
        const double zq = std::clamp(z, probe.lo().z + zNudge, probe.hi().z - zNudge);
        if (!(xq == cachedX_)) {
            cachedStation_ = probe.station(xq);
            cachedX_ = xq;
        }
        return probe.halfBreadth(cachedStation_, xq, zq);
    }

    // Frame stations, the section on each, and the girth layout they share.
    // False means the hull cannot carry structure at all and the passes below have
    // nothing to work on.
    bool layOutFrames();

    void buildShellAndFrames();
    void buildDecks();
    void buildBulkheads();
    void buildGirders();

private:
    double cachedX_ = std::numeric_limits<double>::quiet_NaN();
    std::vector<std::uint32_t> cachedStation_;
};

bool Builder::layOutFrames() {
    if (probe.empty()) {
        report("the hull mesh is empty, so there is no surface to build structure on");
        return false;
    }

    int bays = 1;
    mesh.frameStations = frameStations(probe.lo().x, probe.hi().x, scantlings.frameSpacing, &bays);
    mesh.frameSpacing = (probe.hi().x - probe.lo().x) / bays;
    if (scantlings.frameSpacing > 0 &&
        std::abs(mesh.frameSpacing - scantlings.frameSpacing) > 0.01 * scantlings.frameSpacing)
        report("frame spacing delivered as " + std::to_string(mesh.frameSpacing) +
               " m rather than the requested " + std::to_string(scantlings.frameSpacing) +
               " m, so that a whole number of bays spans the hull");

    sections.reserve(mesh.frameStations.size());
    for (double x : mesh.frameStations)
        sections.push_back(sampleSection(probe, x, scantlings.girthSamples));

    // The girth layout is set amidships and used everywhere, which is what makes a
    // longitudinal one continuous member from end to end.
    std::size_t midship = 0;
    for (std::size_t i = 1; i < sections.size(); ++i)
        if (std::abs(sections[i].x) < std::abs(sections[midship].x)) midship = i;
    const double midshipGirth = sections[midship].girth();
    if (midshipGirth <= 0) {
        report("the hull has no girth amidships, so no shell can be laid out on it");
        return false;
    }
    layout = buildGirthLayout(scantlings, midshipGirth);
    return true;
}

// Shell plating, its longitudinals and the transverse frames.
void Builder::buildShellAndFrames() {
    int uncovered = 0;
    for (std::size_t i = 0; i + 1 < sections.size(); ++i) {
        const Section& aft = sections[i];
        const Section& fwd = sections[i + 1];
        if (aft.point.empty() || fwd.point.empty()) continue;
        const double xMid = 0.5 * (aft.x + fwd.x);

        for (std::size_t j = 0; j + 1 < layout.fraction.size(); ++j) {
            const double g0 = layout.fraction[j], g1 = layout.fraction[j + 1];
            const ShellRegion* region = regionAt(scantlings, xMid, 0.5 * (g0 + g1));
            if (region == nullptr) {
                ++uncovered;
                continue;
            }
            const int source = static_cast<int>(region - scantlings.shell.data());
            const int mat = material(region->material);

            const Vec3 p0 = aft.at(g0), p1 = fwd.at(g0);
            const Vec3 p2 = fwd.at(g1), p3 = aft.at(g1);
            addPanel(mesh, p0, p1, p2, p3, region->thickness, mat, PanelRole::Shell, source);
            addPanel(mesh, mirrored(p3), mirrored(p2), mirrored(p1), mirrored(p0),
                     region->thickness, mat, PanelRole::Shell, source);

            // Each strake carries the longitudinal on its *upper* seam, so the
            // deck edge at girth 1 gets none -- that connection belongs to the
            // deck -- and the keel at girth 0 gets none either, where a centre
            // girder goes instead.
            if (!region->stiffened || j + 2 >= layout.fraction.size()) continue;
            Vec3 riseAft{}, riseFwd{};
            const Vec3 a = aft.at(g1, &riseAft);
            const Vec3 b = fwd.at(g1, &riseFwd);
            const Vec3 rise = normalize(riseAft + riseFwd);
            addMember(mesh, a, b, rise, region->longitudinal, region->thickness, mat,
                      MemberRole::Longitudinal);
            addMember(mesh, mirrored(a), mirrored(b), mirrored(rise), region->longitudinal,
                      region->thickness, mat, MemberRole::Longitudinal);
        }
    }
    if (uncovered > 0)
        report(std::to_string(uncovered) +
               " shell panels had no plating region covering them and were left unplated");

    // --- Transverse frames ----------------------------------------------------

    if (scantlings.framed) {
        for (const Section& section : sections) {
            if (section.point.empty()) continue;
            for (std::size_t j = 0; j + 1 < layout.fraction.size(); ++j) {
                const double g0 = layout.fraction[j], g1 = layout.fraction[j + 1];
                const ShellRegion* region = regionAt(scantlings, section.x, 0.5 * (g0 + g1));
                if (region == nullptr) continue;
                Vec3 rise{};
                const Vec3 a = section.at(g0);
                const Vec3 b = section.at(g1, &rise);
                addMember(mesh, a, b, rise, scantlings.frameProfile, region->thickness,
                          material(scantlings.frameMaterial), MemberRole::Frame);
                addMember(mesh, mirrored(a), mirrored(b), mirrored(rise), scantlings.frameProfile,
                          region->thickness, material(scantlings.frameMaterial), MemberRole::Frame);
            }
        }
    }
}

// Decks, their longitudinals and their beams.
void Builder::buildDecks() {
    for (std::size_t d = 0; d < scantlings.decks.size(); ++d) {
        const Deck& deck = scantlings.decks[d];
        if (deck.z < probe.lo().z - 1e-9 || deck.z > probe.hi().z + 1e-9) {
            report("deck " + deck.name + " at z = " + std::to_string(deck.z) +
                   " m is outside the hull and was skipped");
            continue;
        }
        const double xa = std::max(deck.xFrom, probe.lo().x);
        const double xb = std::min(deck.xTo, probe.hi().x);
        if (xb - xa <= 1e-9) {
            report("deck " + deck.name + " has no length inside the hull and was skipped");
            continue;
        }

        // Bay boundaries: the frame stations inside the deck, with the deck's own
        // ends closing the two part bays. A deck that stops between frames is
        // ordinary; a deck forced onto frame stations would be a lie about where
        // it stops.
        const std::vector<double> edges = baysBetween(mesh.frameStations, xa, xb);

        std::vector<double> halfBreadth(edges.size());
        for (std::size_t i = 0; i < edges.size(); ++i) halfBreadth[i] = breadthAt(edges[i], deck.z);
        const double widest = *std::max_element(halfBreadth.begin(), halfBreadth.end());
        if (widest <= 1e-9) {
            report("deck " + deck.name + " lies outside the shell everywhere and was skipped");
            continue;
        }

        // Deck longitudinals run at fixed |y|, not at a fixed fraction of the
        // breadth: a deck longitudinal is straight and stops where the deck
        // narrows past it. The outboard band of each bay is whatever is left.
        const double spacing = deck.longitudinalSpacing > 0 ? deck.longitudinalSpacing
                                                            : scantlings.longitudinalSpacing;
        const int bands = spacing > 0 ? std::max(1, static_cast<int>(std::ceil(widest / spacing)))
                                      : 1;
        const double pitch = spacing > 0 ? spacing : widest;
        const int mat = material(deck.material);
        const int source = static_cast<int>(d);

        for (std::size_t i = 0; i + 1 < edges.size(); ++i) {
            const double ea = halfBreadth[i], eb = halfBreadth[i + 1];
            for (int j = 0; j < bands; ++j) {
                const double y0 = j * pitch, y1 = (j + 1) * pitch;
                const double a0 = std::min(y0, ea), a1 = std::min(y1, ea);
                const double b0 = std::min(y0, eb), b1 = std::min(y1, eb);
                addPanel(mesh, {edges[i], a0, deck.z}, {edges[i + 1], b0, deck.z},
                         {edges[i + 1], b1, deck.z}, {edges[i], a1, deck.z}, deck.thickness, mat,
                         PanelRole::Deck, source);
                addPanel(mesh, {edges[i], -a0, deck.z}, {edges[i], -a1, deck.z},
                         {edges[i + 1], -b1, deck.z}, {edges[i + 1], -b0, deck.z}, deck.thickness,
                         mat, PanelRole::Deck, source);
            }

            // Stiffening hangs below the deck, which is where the eccentricity of
            // an offset beam points and half of what sets the deck's own section.
            const Vec3 down{0, 0, -1};
            if (deck.stiffened)
                for (int j = 1; j < bands; ++j) {
                    const double y = j * pitch;
                    if (y > ea - 1e-9 || y > eb - 1e-9) continue;
                    addMember(mesh, {edges[i], y, deck.z}, {edges[i + 1], y, deck.z}, down,
                              deck.longitudinal, deck.thickness, mat,
                              MemberRole::DeckLongitudinal);
                    addMember(mesh, {edges[i], -y, deck.z}, {edges[i + 1], -y, deck.z}, down,
                              deck.longitudinal, deck.thickness, mat,
                              MemberRole::DeckLongitudinal);
                }
        }

        if (deck.beamed)
            for (std::size_t i = 0; i < edges.size(); ++i) {
                const Vec3 down{0, 0, -1};
                addMember(mesh, {edges[i], 0, deck.z}, {edges[i], halfBreadth[i], deck.z}, down,
                          deck.beam, deck.thickness, mat, MemberRole::DeckBeam);
                addMember(mesh, {edges[i], 0, deck.z}, {edges[i], -halfBreadth[i], deck.z}, down,
                          deck.beam, deck.thickness, mat, MemberRole::DeckBeam);
            }
    }
}

// Watertight and structural bulkheads, transverse and longitudinal.
void Builder::buildBulkheads() {
    for (std::size_t b = 0; b < scantlings.bulkheads.size(); ++b) {
        const Bulkhead& bulkhead = scantlings.bulkheads[b];
        const double z0 = std::max(bulkhead.zFrom, probe.lo().z);
        const double z1 = std::min(bulkhead.zTo, probe.hi().z);
        if (z1 - z0 <= 1e-9) {
            report("bulkhead " + bulkhead.name + " has no height inside the hull and was skipped");
            continue;
        }
        const double spacing = bulkhead.stiffenerSpacing > 0 ? bulkhead.stiffenerSpacing
                                                             : scantlings.longitudinalSpacing;
        const int mat = material(bulkhead.material);
        const int source = static_cast<int>(b);
        const int zBands = spacing > 0 ? std::max(1, static_cast<int>(std::ceil((z1 - z0) / spacing)))
                                       : 1;

        if (bulkhead.transverse) {
            // A transverse bulkhead is stiffened vertically, so the plate panels
            // are bounded left and right by stiffeners. The horizontal seams are
            // a mesh subdivision and nothing else -- there is no stringer on them,
            // and this representation cannot put one there.
            std::vector<double> zs(static_cast<std::size_t>(zBands) + 1);
            std::vector<double> half(zs.size());
            for (std::size_t k = 0; k < zs.size(); ++k) {
                zs[k] = z0 + (z1 - z0) * static_cast<double>(k) / zBands;
                half[k] = breadthAt(bulkhead.position, zs[k]);
            }
            const double widest = *std::max_element(half.begin(), half.end());
            if (widest <= 1e-9) {
                report("bulkhead " + bulkhead.name + " lies outside the shell and was skipped");
                continue;
            }
            const int yBands =
                spacing > 0 ? std::max(1, static_cast<int>(std::ceil(widest / spacing))) : 1;
            const double pitch = spacing > 0 ? spacing : widest;

            for (int k = 0; k < zBands; ++k)
                for (int j = 0; j < yBands; ++j) {
                    const double y0 = j * pitch, y1 = (j + 1) * pitch;
                    const double a0 = std::min(y0, half[k]), a1 = std::min(y1, half[k]);
                    const double b0 = std::min(y0, half[k + 1]), b1 = std::min(y1, half[k + 1]);
                    const double x = bulkhead.position;
                    addPanel(mesh, {x, a0, zs[k]}, {x, a1, zs[k]}, {x, b1, zs[k + 1]},
                             {x, b0, zs[k + 1]}, bulkhead.thickness, mat, PanelRole::Bulkhead,
                             source);
                    addPanel(mesh, {x, -a0, zs[k]}, {x, -b0, zs[k + 1]}, {x, -b1, zs[k + 1]},
                             {x, -a1, zs[k]}, bulkhead.thickness, mat, PanelRole::Bulkhead, source);
                }

            if (bulkhead.stiffened)
                for (int k = 0; k < zBands; ++k)
                    for (int j = 1; j < yBands; ++j) {
                        const double y = j * pitch;
                        if (y > half[k] - 1e-9 || y > half[k + 1] - 1e-9) continue;
                        const double x = bulkhead.position;
                        addMember(mesh, {x, y, zs[k]}, {x, y, zs[k + 1]}, {1, 0, 0},
                                  bulkhead.stiffener, bulkhead.thickness, mat,
                                  MemberRole::BulkheadStiffener);
                        addMember(mesh, {x, -y, zs[k]}, {x, -y, zs[k + 1]}, {1, 0, 0},
                                  bulkhead.stiffener, bulkhead.thickness, mat,
                                  MemberRole::BulkheadStiffener);
                    }
            continue;
        }

        // A longitudinal bulkhead sits at a fixed |y| and is stiffened
        // longitudinally, which is why -- unlike a transverse one -- it works in
        // the hull girder. Where the hull has narrowed inside it the panel simply
        // does not exist; it is dropped and counted rather than left poking out
        // through the shell, which is silent and looks like extra steel.
        const double xa = std::max(bulkhead.xFrom, probe.lo().x);
        const double xb = std::min(bulkhead.xTo, probe.hi().x);
        if (xb - xa <= 1e-9) {
            report("longitudinal bulkhead " + bulkhead.name +
                   " has no length inside the hull and was skipped");
            continue;
        }
        const std::vector<double> edges = baysBetween(mesh.frameStations, xa, xb);

        // Whether the shell is outboard of the bulkhead at every grid corner,
        // tabulated with x outermost so the station cache is hit rather than
        // thrashed.
        const std::size_t zEdges = static_cast<std::size_t>(zBands) + 1;
        std::vector<char> fits(edges.size() * zEdges, 0);
        for (std::size_t i = 0; i < edges.size(); ++i)
            for (std::size_t k = 0; k < zEdges; ++k) {
                const double z = z0 + (z1 - z0) * static_cast<double>(k) / zBands;
                fits[i * zEdges + k] = breadthAt(edges[i], z) >= bulkhead.position ? 1 : 0;
            }

        int outside = 0;
        for (std::size_t i = 0; i + 1 < edges.size(); ++i)
            for (int k = 0; k < zBands; ++k) {
                const std::size_t k0 = static_cast<std::size_t>(k);
                const double za = z0 + (z1 - z0) * k / zBands;
                const double zb = z0 + (z1 - z0) * (k + 1) / zBands;
                const bool inside = fits[i * zEdges + k0] && fits[i * zEdges + k0 + 1] &&
                                    fits[(i + 1) * zEdges + k0] && fits[(i + 1) * zEdges + k0 + 1];
                if (!inside) {
                    ++outside;
                    continue;
                }
                const double y = bulkhead.position;
                addPanel(mesh, {edges[i], y, za}, {edges[i + 1], y, za}, {edges[i + 1], y, zb},
                         {edges[i], y, zb}, bulkhead.thickness, mat, PanelRole::Bulkhead, source);
                addPanel(mesh, {edges[i], -y, za}, {edges[i], -y, zb}, {edges[i + 1], -y, zb},
                         {edges[i + 1], -y, za}, bulkhead.thickness, mat, PanelRole::Bulkhead,
                         source);
                if (!bulkhead.stiffened || k + 1 >= zBands) continue;
                addMember(mesh, {edges[i], y, zb}, {edges[i + 1], y, zb}, {0, -1, 0},
                          bulkhead.stiffener, bulkhead.thickness, mat,
                          MemberRole::BulkheadStiffener);
                addMember(mesh, {edges[i], -y, zb}, {edges[i + 1], -y, zb}, {0, 1, 0},
                          bulkhead.stiffener, bulkhead.thickness, mat,
                          MemberRole::BulkheadStiffener);
            }
        if (outside > 0)
            report("longitudinal bulkhead " + bulkhead.name + " has " + std::to_string(outside) +
                   " panels outside the shell; shorten it or move it inboard");
    }
}

// Deep longitudinal members on named lines.
void Builder::buildGirders() {
    for (const Girder& girder : scantlings.girders) {
        const double xa = std::max(girder.xFrom, probe.lo().x);
        const double xb = std::min(girder.xTo, probe.hi().x);
        if (xb - xa <= 1e-9) {
            report("girder " + girder.name + " has no length inside the hull and was skipped");
            continue;
        }
        const std::vector<double> edges = baysBetween(mesh.frameStations, xa, xb);

        const int mat = material(girder.material);
        std::vector<char> fits(edges.size(), 0);
        for (std::size_t i = 0; i < edges.size(); ++i)
            fits[i] = breadthAt(edges[i], girder.z) >= std::abs(girder.y) ? 1 : 0;

        int outside = 0;
        for (std::size_t i = 0; i + 1 < edges.size(); ++i) {
            if (!fits[i] || !fits[i + 1]) {
                ++outside;
                continue;
            }
            addMember(mesh, {edges[i], girder.y, girder.z}, {edges[i + 1], girder.y, girder.z},
                      girder.rise, girder.profile, girder.attachedPlateThickness, mat,
                      MemberRole::Girder);
            if (girder.bothSides && std::abs(girder.y) > 1e-9)
                addMember(mesh, {edges[i], -girder.y, girder.z},
                          {edges[i + 1], -girder.y, girder.z}, mirrored(girder.rise),
                          girder.profile, girder.attachedPlateThickness, mat, MemberRole::Girder);
        }
        if (outside > 0)
            report("girder " + girder.name + " runs outside the shell over " +
                   std::to_string(outside) + " bays, which were dropped");
    }
}

}  // namespace

StructuralMesh makeStructuralMesh(const TriMesh& hull, const Scantlings& scantlings,
                                  std::vector<std::string>* problems) {
    Builder builder(hull, scantlings);
    // Whatever is wrong with the description on its own is said first, so one call
    // reports both what could not be described and what could not be built.
    for (const std::string& s : validateScantlings(scantlings)) builder.report(s);

    if (builder.layOutFrames()) {
        builder.buildShellAndFrames();
        builder.buildDecks();
        builder.buildBulkheads();
        builder.buildGirders();
    }

    if (problems)
        for (const std::string& s : builder.problems) problems->push_back(s);
    return std::move(builder.mesh);
}

// --- Hull girder ---------------------------------------------------------------

HullGirderSection hullGirderSection(const StructuralMesh& mesh, double x) {
    HullGirderSection s;
    s.x = x;

    double area = 0, moment = 0, second = 0;
    double zLo = 1e300, zHi = -1e300;
    constexpr double kFlat = 1e-9;

    // A cut that lands exactly on a frame station sits on the seam between two
    // bays, and both of them straddle it. Counting both doubles the whole
    // section -- area, second moment, everything except the neutral axis, which
    // is a ratio and so looks perfectly correct. The extent is therefore taken
    // half-open, [xLo, xHi), so a seam is served by the bay forward of it; only
    // at the forward end of the ship, where there is no such bay, does it close.
    double foreEnd = -1e300;
    for (const PlatePanel& p : mesh.panels)
        for (const Vec3& c : p.corner) foreEnd = std::max(foreEnd, c.x);
    for (const StructuralMember& m : mesh.members)
        foreEnd = std::max(foreEnd, std::max(m.a.x, m.b.x));
    const bool atForeEnd = x >= foreEnd - kFlat;
    const auto straddles = [&](double xLo, double xHi) {
        if (xHi - xLo <= kFlat) return false;
        if (x < xLo - kFlat) return false;
        return atForeEnd ? x <= xHi + kFlat : x < xHi - kFlat;
    };

    for (const PlatePanel& p : mesh.panels) {
        double xLo = p.corner[0].x, xHi = p.corner[0].x;
        for (const Vec3& c : p.corner) {
            xLo = std::min(xLo, c.x);
            xHi = std::max(xHi, c.x);
        }
        // Zero extent along x means the panel *is* the cut plane: a transverse
        // bulkhead, which carries no longitudinal stress. This is decided by
        // geometry rather than by the panel's label, so mis-tagging one cannot
        // put it into the hull girder.
        if (!straddles(xLo, xHi)) continue;

        // The cut through a convex quad is a segment; take the two crossings
        // furthest apart, which is robust when an edge lies in the plane.
        Vec3 crossing[4];
        int count = 0;
        for (int e = 0; e < 4 && count < 4; ++e) {
            const Vec3& a = p.corner[e];
            const Vec3& b = p.corner[(e + 1) % 4];
            const double da = a.x - x, db = b.x - x;
            if ((da <= 0 && db >= 0) || (da >= 0 && db <= 0)) {
                const double denominator = da - db;
                const double u = std::abs(denominator) > 1e-15 ? da / denominator : 0.0;
                crossing[count++] = a + (b - a) * u;
            }
        }
        if (count < 2) continue;
        int bestI = 0, bestJ = 1;
        double bestLength = -1;
        for (int i = 0; i < count; ++i)
            for (int j = i + 1; j < count; ++j) {
                const double l = length(crossing[i] - crossing[j]);
                if (l > bestLength) {
                    bestLength = l;
                    bestI = i;
                    bestJ = j;
                }
            }
        const Vec3 q0 = crossing[bestI], q1 = crossing[bestJ];
        const double cut = bestLength;
        if (cut <= 1e-12) continue;

        const double a = p.thickness * cut;
        const double zc = 0.5 * (q0.z + q1.z);
        const double dz = q1.z - q0.z;
        const double dy = q1.y - q0.y;
        // A rectangle of length `cut` and thickness t, inclined: the in-plane and
        // through-thickness terms rotate into the horizontal axis together.
        const double own =
            a * (dz * dz + p.thickness * p.thickness * (dy / cut) * (dy / cut)) / 12.0;

        area += a;
        moment += a * zc;
        second += own + a * zc * zc;
        const double halfThickness = 0.5 * p.thickness * std::abs(dy) / cut;
        zLo = std::min(zLo, std::min(q0.z, q1.z) - halfThickness);
        zHi = std::max(zHi, std::max(q0.z, q1.z) + halfThickness);
    }

    for (const StructuralMember& m : mesh.members) {
        const double xLo = std::min(m.a.x, m.b.x), xHi = std::max(m.a.x, m.b.x);
        if (!straddles(xLo, xHi)) continue;

        const double u = std::clamp((x - m.a.x) / (m.b.x - m.a.x), 0.0, 1.0);
        const Vec3 root = m.a + (m.b - m.a) * u;
        const ProfileSection ps = profileSection(m.profile);
        if (ps.area <= 0) continue;
        const Vec3 centre = root + m.rise * (0.5 * m.attachedPlateThickness + ps.centroid);

        // The profile's principal axes are along the web and across it; the hull
        // girder wants the horizontal one, so the two rotate in by the web's
        // direction cosines. A tee on the bottom shell contributes its strong
        // axis; the same tee on the side shell contributes its weak one, which is
        // most of why a side longitudinal is worth so much less than a bottom one.
        const double own = m.rise.z * m.rise.z * ps.secondMoment +
                           m.rise.y * m.rise.y * ps.secondMomentWeak;

        area += ps.area;
        moment += ps.area * centre.z;
        second += own + ps.area * centre.z * centre.z;
        const Vec3 tip = root + m.rise * (0.5 * m.attachedPlateThickness + ps.height);
        zLo = std::min(zLo, std::min(root.z, tip.z));
        zHi = std::max(zHi, std::max(root.z, tip.z));
    }

    if (area <= 0) return s;
    s.area = area;
    s.neutralAxis = moment / area;
    s.secondMoment = second - area * s.neutralAxis * s.neutralAxis;
    s.zKeel = zLo;
    s.zDeck = zHi;
    if (s.zDeck - s.neutralAxis > 1e-9) s.modulusDeck = s.secondMoment / (s.zDeck - s.neutralAxis);
    if (s.neutralAxis - s.zKeel > 1e-9) s.modulusKeel = s.secondMoment / (s.neutralAxis - s.zKeel);
    return s;
}

double ruleMinimumSectionModulus(double lengthPp, double beam, double blockCoefficient) {
    if (lengthPp <= 0 || beam <= 0) return 0.0;
    // IACS unified requirement S7: the block coefficient is not taken below 0.60,
    // because the requirement is a strength floor and a fine hull is not excused
    // from it.
    const double cb = std::max(0.60, blockCoefficient);
    double c = 10.75;
    if (lengthPp < 300.0) c = 10.75 - std::pow((300.0 - lengthPp) / 100.0, 1.5);
    else if (lengthPp > 350.0) c = 10.75 - std::pow((lengthPp - 350.0) / 150.0, 1.5);
    return c * lengthPp * lengthPp * beam * (cb + 0.7) * 1.0e-6;
}

// --- Validation -----------------------------------------------------------------

std::vector<std::string> validateScantlings(const Scantlings& s) {
    std::vector<std::string> problems;
    if (s.frameSpacing <= 0) problems.push_back("frame spacing is not positive");
    if (s.longitudinalSpacing <= 0) problems.push_back("longitudinal spacing is not positive");
    if (s.girthSamples < 8)
        problems.push_back("fewer than 8 girth samples will not resolve a section");
    if (s.materials.empty()) problems.push_back("no materials are defined");
    if (s.shell.empty()) problems.push_back("no shell plating is defined");

    const int materialCount = static_cast<int>(s.materials.size());
    const auto checkMaterial = [&](int index, const std::string& what) {
        if (index < 0 || index >= materialCount)
            problems.push_back(what + " names material " + std::to_string(index) +
                               ", which does not exist");
    };

    for (const ShellRegion& r : s.shell) {
        if (r.thickness <= 0)
            problems.push_back("shell region " + r.name + " has a non-positive thickness");
        if (r.girthTo <= r.girthFrom)
            problems.push_back("shell region " + r.name + " covers no girth");
        if (r.girthFrom < -1e-12 || r.girthTo > 1.0 + 1e-12)
            problems.push_back("shell region " + r.name + " leaves the girth range [0, 1]");
        checkMaterial(r.material, "shell region " + r.name);
    }

    // Every girth fraction must be covered at every station, or there is a hole in
    // the shell. The check runs over the grid of *declared* breakpoints in both x
    // and girth, because a region's girth range and its length range are
    // independent: a strake that is complete amidships and absent forward is a
    // gap that a girth-only check cannot see. It is exactly the arrangement the
    // reference ferry uses -- thinner side plating at the ends, declared as its
    // own region -- so getting this wrong would look like a working description.
    std::vector<double> girthBreaks{0.0, 1.0};
    std::vector<double> lengthBreaks;
    for (const ShellRegion& r : s.shell) {
        girthBreaks.push_back(r.girthFrom);
        girthBreaks.push_back(r.girthTo);
        lengthBreaks.push_back(r.xFrom);
        lengthBreaks.push_back(r.xTo);
    }
    std::sort(girthBreaks.begin(), girthBreaks.end());
    std::sort(lengthBreaks.begin(), lengthBreaks.end());
    int gaps = 0;
    for (std::size_t i = 0; i + 1 < girthBreaks.size(); ++i) {
        const double girth = 0.5 * (girthBreaks[i] + girthBreaks[i + 1]);
        if (girthBreaks[i + 1] - girthBreaks[i] < 1e-9 || girth < 0.0 || girth > 1.0) continue;
        for (std::size_t j = 0; j + 1 < lengthBreaks.size(); ++j) {
            if (lengthBreaks[j + 1] - lengthBreaks[j] < 1e-9) continue;
            const double x = 0.5 * (lengthBreaks[j] + lengthBreaks[j + 1]);
            bool covered = false;
            for (const ShellRegion& r : s.shell)
                if (girth >= r.girthFrom && girth <= r.girthTo && x >= r.xFrom && x <= r.xTo)
                    covered = true;
            if (covered) continue;
            if (++gaps <= 4)
                problems.push_back("no shell region covers girth fraction " +
                                   std::to_string(girth) + " at x = " + std::to_string(x));
        }
    }
    if (gaps > 4)
        problems.push_back(std::to_string(gaps - 4) + " further gaps in the shell plating");

    for (const Deck& d : s.decks) {
        if (d.thickness <= 0) problems.push_back("deck " + d.name + " has a non-positive thickness");
        if (d.xTo <= d.xFrom) problems.push_back("deck " + d.name + " has no length");
        checkMaterial(d.material, "deck " + d.name);
    }
    for (const Bulkhead& b : s.bulkheads) {
        if (b.thickness <= 0)
            problems.push_back("bulkhead " + b.name + " has a non-positive thickness");
        if (b.zTo <= b.zFrom) problems.push_back("bulkhead " + b.name + " has no height");
        if (!b.transverse && b.xTo <= b.xFrom)
            problems.push_back("longitudinal bulkhead " + b.name + " has no length");
        if (!b.transverse && b.position <= 0)
            problems.push_back("longitudinal bulkhead " + b.name +
                               " is on the centreline; give it a positive |y|");
        checkMaterial(b.material, "bulkhead " + b.name);
    }
    for (const Girder& g : s.girders) {
        if (g.xTo <= g.xFrom) problems.push_back("girder " + g.name + " has no length");
        if (length(g.rise) < 1e-9)
            problems.push_back("girder " + g.name + " has no direction for its web to rise in");
        checkMaterial(g.material, "girder " + g.name);
    }
    return problems;
}

// --- Reference ship ---------------------------------------------------------------

Scantlings ferryScantlings() {
    Scantlings s;
    s.materials = {ah36Steel(), mildSteel()};

    // Longitudinally framed with transverse web frames every four 600 mm frame
    // spaces, which is ordinary for the size. Ordinary frames between the webs
    // are not expressible here -- see docs/02-simulation.md section 3.
    s.frameSpacing = 2.40;
    s.frameProfile = tee(0.700, 0.011, 0.150, 0.014);
    s.frameMaterial = 0;
    s.longitudinalSpacing = 0.70;

    // Girth fractions: the ferry's midship girth runs about 3 m across the flat of
    // bottom, 5 m round the bilge to the turn at z = 4.2 m, then 10.8 m of wall
    // side to the sheer at 15 m. So the bottom is the first eighth, the bilge the
    // next seventh, and the side everything above -- with the sheer strake, which
    // is always thickened because it is the extreme fibre of the hull girder,
    // taking the top 4%.
    ShellRegion bottom;
    bottom.name = "bottom";
    bottom.girthFrom = 0.00;
    bottom.girthTo = 0.13;
    bottom.thickness = 0.0145;
    bottom.longitudinal = tee(0.200, 0.010, 0.090, 0.012);
    bottom.longitudinalSpacing = 0.70;

    ShellRegion bilge;
    bilge.name = "bilge";
    bilge.girthFrom = 0.13;
    bilge.girthTo = 0.28;
    bilge.thickness = 0.0155;
    bilge.longitudinal = tee(0.200, 0.010, 0.090, 0.012);
    bilge.longitudinalSpacing = 0.70;

    ShellRegion side;
    side.name = "side";
    side.girthFrom = 0.28;
    side.girthTo = 0.96;
    side.thickness = 0.0120;
    side.longitudinal = angle(0.160, 0.009, 0.070, 0.011);
    side.longitudinalSpacing = 0.70;

    ShellRegion sheer;
    sheer.name = "sheer_strake";
    sheer.girthFrom = 0.96;
    sheer.girthTo = 1.00;
    sheer.thickness = 0.0140;
    sheer.longitudinal = angle(0.160, 0.009, 0.070, 0.011);
    sheer.longitudinalSpacing = 0.70;

    // The ends carry less load and less plating; the midship two thirds carries
    // the hull girder. Declared after the strakes so that it overrides them.
    ShellRegion forwardEnd = side;
    forwardEnd.name = "side_forward";
    forwardEnd.xFrom = 44.0;
    forwardEnd.thickness = 0.0100;
    ShellRegion aftEnd = side;
    aftEnd.name = "side_aft";
    aftEnd.xTo = -44.0;
    aftEnd.thickness = 0.0100;

    s.shell = {bottom, bilge, side, sheer, forwardEnd, aftEnd};

    // Decks, on the same levels as the compartment plan in ships/ferry.ship: the
    // double bottom top at 1.8 m, the bulkhead deck at 7.0 m, the vehicle deck
    // head at 12.5 m and the weather deck at the sheer.
    Deck tankTop;
    tankTop.name = "tank_top";
    tankTop.z = 1.80;
    tankTop.xFrom = -44.0;
    tankTop.xTo = 44.0;
    tankTop.thickness = 0.0120;
    tankTop.longitudinal = tee(0.180, 0.009, 0.080, 0.011);
    tankTop.beam = tee(0.450, 0.010, 0.120, 0.012);

    Deck bulkheadDeck;
    bulkheadDeck.name = "bulkhead_deck";
    bulkheadDeck.z = 7.00;
    bulkheadDeck.xFrom = -60.0;
    bulkheadDeck.xTo = 60.0;
    bulkheadDeck.thickness = 0.0110;
    bulkheadDeck.longitudinal = tee(0.180, 0.009, 0.080, 0.011);
    bulkheadDeck.beam = tee(0.500, 0.010, 0.130, 0.012);

    Deck vehicleDeckHead;
    vehicleDeckHead.name = "vehicle_deck_head";
    vehicleDeckHead.z = 12.50;
    vehicleDeckHead.xFrom = -50.0;
    vehicleDeckHead.xTo = 50.0;
    vehicleDeckHead.thickness = 0.0080;
    vehicleDeckHead.longitudinal = flatBar(0.140, 0.008);
    vehicleDeckHead.beam = tee(0.400, 0.009, 0.110, 0.011);

    Deck weatherDeck;
    weatherDeck.name = "weather_deck";
    weatherDeck.z = 15.00;
    weatherDeck.xFrom = -40.0;
    weatherDeck.xTo = 34.0;
    weatherDeck.thickness = 0.0080;
    weatherDeck.material = 1;  // mild steel: the weather deck is not a strength member here
    weatherDeck.longitudinal = flatBar(0.140, 0.008);
    weatherDeck.beam = tee(0.400, 0.009, 0.110, 0.011);

    s.decks = {tankTop, bulkheadDeck, vehicleDeckHead, weatherDeck};

    // Watertight subdivision, on the same stations as ships/ferry.ship's
    // compartment boundaries.
    const double stations[] = {-44.0, -38.0, -8.0, 20.0, 44.0};
    for (double x : stations) {
        Bulkhead b;
        b.name = "wt_" + std::to_string(static_cast<int>(x));
        b.transverse = true;
        b.position = x;
        b.zFrom = 0.0;
        b.zTo = 7.0;
        b.thickness = 0.0095;
        b.stiffener = flatBar(0.180, 0.010);
        b.stiffenerSpacing = 0.70;
        s.bulkheads.push_back(b);
    }

    // The inboard boundary of the wing tanks. It cannot follow the compartment
    // plan's y = 8 m: the hull is narrower than that at the tank top over the
    // whole length, so a bulkhead there would hang outside the shell.
    Bulkhead wing;
    wing.name = "wing_bulkhead";
    wing.transverse = false;
    wing.position = 6.0;
    wing.zFrom = 1.80;
    wing.zTo = 7.00;
    wing.xFrom = -34.0;
    wing.xTo = 34.0;
    wing.thickness = 0.0095;
    wing.stiffener = angle(0.150, 0.009, 0.070, 0.010);
    wing.stiffenerSpacing = 0.70;
    s.bulkheads.push_back(wing);

    Girder centre;
    centre.name = "centre_girder";
    centre.y = 0.0;
    centre.z = 0.0;
    centre.xFrom = -50.0;
    centre.xTo = 50.0;
    centre.profile = flatBar(1.800, 0.012);
    centre.rise = {0, 0, 1};
    centre.attachedPlateThickness = 0.0145;
    s.girders.push_back(centre);

    Girder sideGirder;
    sideGirder.name = "side_girder";
    sideGirder.y = 2.0;
    sideGirder.z = 0.0;
    sideGirder.xFrom = -44.0;
    sideGirder.xTo = 44.0;
    sideGirder.bothSides = true;
    sideGirder.profile = flatBar(1.800, 0.010);
    sideGirder.rise = {0, 0, 1};
    sideGirder.attachedPlateThickness = 0.0145;
    s.girders.push_back(sideGirder);

    Girder deckGirder;
    deckGirder.name = "deck_girder";
    deckGirder.y = 5.0;
    deckGirder.z = 12.5;
    deckGirder.xFrom = -50.0;
    deckGirder.xTo = 50.0;
    deckGirder.bothSides = true;
    deckGirder.profile = tee(0.600, 0.010, 0.180, 0.014);
    deckGirder.rise = {0, 0, -1};
    deckGirder.attachedPlateThickness = 0.0080;
    s.girders.push_back(deckGirder);

    return s;
}

}  // namespace sim
