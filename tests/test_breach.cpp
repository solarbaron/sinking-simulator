// SPDX-License-Identifier: MIT
//
// Validation of failed structure -> holes in the flooding network.
//
// The milestone this serves is `docs/06-roadmap.md`'s: the hull tears, "and the
// resulting hole floods at a rate the hole's own area determines". So the
// load-bearing test here is not that an `Opening` comes out with plausible
// fields -- it is that the opening, handed to the real flooding solver, passes
// `Cd·A·√(2Δp/ρ)`, and that doubling the failed plating doubles the water.
//
// Three instruments, deliberately different from each other:
//
//   * **A rectangular barge**, where every panel is planar and every compartment
//     is a box, so area, centroid and head are exact rationals rather than
//     tessellations of something. Merging and connectivity are asserted here
//     against arithmetic anyone can do on paper.
//   * **The reference ferry**, whose hull is curved, whose compartments do not
//     reach the shell everywhere, and whose structural mesh has 8 900 panels laid
//     out by a generator that knows nothing about compartments. This is the only
//     thing that says the probe survives a real ship; it is also where the
//     three-way answer (compartment / sea / unmodelled void) earns its keep.
//   * **An authored twin.** The same ship is flooded twice, once through an
//     opening this file produced and once through an `Opening` written out by
//     hand with the same numbers, and the two water volumes are required to agree
//     to the last bit. That is the whole "indistinguishable to the solver" claim,
//     stated as an equality rather than as a hope.
#include "engine/sim/breach.hpp"

#include "engine/core/geometry.hpp"
#include "engine/sim/scantlings.hpp"
#include "game/prototype/ferry.hpp"
#include "harness.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

using namespace sim;
using testing::expectEqual;
using testing::expectNear;
using testing::expectTrue;

namespace {

// --- Panel construction -------------------------------------------------------
//
// Rectangles, wound the way scantlings.cpp winds them: corner order round the
// panel. Every one has area (span x span) and centroid at its own middle, both
// exactly, which is what makes the assertions below arithmetic rather than
// tolerance-chasing.

PlatePanel panelAtY(double y, double x0, double x1, double z0, double z1) {
    PlatePanel p;
    p.corner[0] = {x0, y, z0};
    p.corner[1] = {x1, y, z0};
    p.corner[2] = {x1, y, z1};
    p.corner[3] = {x0, y, z1};
    p.thickness = 0.012;
    p.role = PanelRole::Shell;
    return p;
}

PlatePanel panelAtX(double x, double y0, double y1, double z0, double z1) {
    PlatePanel p;
    p.corner[0] = {x, y0, z0};
    p.corner[1] = {x, y1, z0};
    p.corner[2] = {x, y1, z1};
    p.corner[3] = {x, y0, z1};
    p.thickness = 0.009;
    p.role = PanelRole::Bulkhead;
    return p;
}

PlatePanel panelAtZ(double z, double x0, double x1, double y0, double y1) {
    PlatePanel p;
    p.corner[0] = {x0, y0, z};
    p.corner[1] = {x1, y0, z};
    p.corner[2] = {x1, y1, z};
    p.corner[3] = {x0, y1, z};
    p.thickness = 0.010;
    p.role = PanelRole::Deck;
    return p;
}

// An L-shaped prism: a closed mesh whose notch is inside its bounding box and
// inside its convex hull but outside the solid. Nothing built from half-spaces
// can be this shape, which is the point of testing against it.
TriMesh lPrism() {
    const double footprint[6][2] = {{0, 0}, {4, 0}, {4, 1}, {1, 1}, {1, 4}, {0, 4}};
    TriMesh m;
    for (int k = 0; k < 2; ++k)
        for (const auto& p : footprint) m.verts.push_back({p[0], p[1], static_cast<double>(k)});
    const auto tri = [&](int a, int b, int c) {
        m.tris.push_back({static_cast<std::uint32_t>(a), static_cast<std::uint32_t>(b),
                          static_cast<std::uint32_t>(c)});
    };
    // Fanned from the inside corner, which is the one vertex every other vertex
    // of an L can see. The caps are wound against the sides so that each shared
    // edge is traversed once each way, which is what makes the mesh closed.
    for (int i = 1; i + 1 < 6; ++i) {
        tri(0, i + 1, i);                      // bottom, facing -z
        tri(6, 6 + i, 6 + i + 1);              // top, facing +z
    }
    for (int i = 0; i < 6; ++i) {              // sides
        const int j = (i + 1) % 6;
        tri(i, j, 6 + j);
        tri(i, 6 + j, 6 + i);
    }
    // Wind it outward by the sign of its own volume rather than by reasoning
    // about the fan: getting that wrong is exactly the mistake the winding number
    // is being asked to detect, so it must not be assumed here.
    if (integrate(m).volume < 0)
        for (Tri& t : m.tris) std::swap(t.b, t.c);
    return m;
}

Compartment box(const char* name, Vec3 lo, Vec3 hi) {
    Compartment c;
    c.name = name;
    c.mesh = makeBox(lo, hi);
    c.permeability = 1.0;
    c.ventedToAtmosphere = true;  // air is a separate physics; keep it out of the way
    return c;
}

// --- The barge ----------------------------------------------------------------
//
// 100 x 20 x 10 m, two holds meeting on a bulkhead at midship, both of them
// reaching the shell all round and stopping at z = 6 so that the space above them
// is inside the hull and inside nothing -- the unmodelled void the probe has to
// tell apart from the sea.

constexpr double kBargeLength = 100.0, kBargeBeam = 20.0, kBargeDepth = 10.0;
constexpr double kBargeDraft = 4.0;
constexpr double kHoldTop = 6.0;

// Panel indices in the barge's structural mesh. Named because every assertion
// below is about a particular piece of plating and "panel 7" says nothing.
enum PanelIndex {
    // A 4 x 2 grid of 2 x 1 m panels on the starboard shell of the aft hold,
    // x = -40..-32, z = 1..3. Column-major: kGrid + 2*column + row.
    kGrid = 0,
    kGridPanels = 8,
    kBulkhead = 8,     // on the midship bulkhead, starboard side
    kShellAtBulkhead,  // shell plating whose forward edge is the bulkhead's outboard edge
    kBottom,           // flat of bottom, aft hold
    kHoldTopPlate,     // the top of the aft hold: compartment below, unmodelled void above
    kWingBulkhead,     // a longitudinal bulkhead buried inside the aft hold
    kDegenerate,       // a collapsed quad, which has no plane to probe either side of
    // Two shell panels of *different* areas sharing a seam, 6 m2 under 3 m2. An
    // unweighted mean of their centres would put the hole 125 mm too high, so
    // this pair is the one that can tell the two rules apart.
    kWide,
    kNarrow,
    // A panel chording across the corner where the side meets the bottom, so its
    // own centroid is half a metre inside the hull and a short probe reads the
    // aft hold on both sides. This is the barge's version of what the ferry does
    // at the turn of the bilge, and the only panel here that needs the march.
    kChord,
    // Two triangular panels -- quads with a collapsed side, which is what the
    // generator produces where a deck runs out against the shell -- meeting at
    // the collapsed corner and nowhere else.
    kTriangleA,
    kTriangleB,
    // A seam whose two sides are a nanometre apart instead of identical, and
    // placed so the two corners fall in different cells of the weld grid.
    kSeamBelow,
    kSeamAbove,
    kPanelCount
};

struct Barge {
    Ship ship;
    StructuralMesh mesh;
    int aftHold = 0, fwdHold = 0;
};

Barge makeBarge() {
    Barge barge;
    Ship& s = barge.ship;
    s.hull = makeBox({-0.5 * kBargeLength, -0.5 * kBargeBeam, 0},
                     {0.5 * kBargeLength, 0.5 * kBargeBeam, kBargeDepth});
    s.deckEdgeZ = kHoldTop;
    s.compartments = {
        box("aft_hold", {-0.5 * kBargeLength, -0.5 * kBargeBeam, 0}, {0, 0.5 * kBargeBeam, kHoldTop}),
        box("fwd_hold", {0, -0.5 * kBargeBeam, 0},
            {0.5 * kBargeLength, 0.5 * kBargeBeam, kHoldTop}),
    };
    s.lightshipMass = kBargeLength * kBargeBeam * kBargeDraft * kRhoSeawater;
    s.lightshipCog = {0, 0, 3.0};
    s.gyradii = {7.0, 28.0, 28.0};
    s.initialise(0.0);
    barge.aftHold = s.findCompartment("aft_hold");
    barge.fwdHold = s.findCompartment("fwd_hold");

    const double y = -0.5 * kBargeBeam;
    std::vector<PlatePanel>& panels = barge.mesh.panels;
    for (int column = 0; column < 4; ++column)
        for (int row = 0; row < 2; ++row)
            panels.push_back(panelAtY(y, -40.0 + 2.0 * column, -38.0 + 2.0 * column, 1.0 + row,
                                      2.0 + row));
    panels.push_back(panelAtX(0.0, y, y + 2.0, 1.0, 3.0));
    panels.push_back(panelAtY(y, -2.0, 0.0, 1.0, 3.0));
    panels.push_back(panelAtZ(0.0, -20.0, -18.0, y, y + 2.0));
    panels.push_back(panelAtZ(kHoldTop, -20.0, -18.0, -1.0, 1.0));
    panels.push_back(panelAtX(-25.0, -5.0, -3.0, 1.0, 3.0));
    panels.push_back(panelAtY(y, -10.0, -10.0, 1.0, 1.0));  // collapsed
    panels.push_back(panelAtY(y, -16.0, -10.0, 1.0, 2.0));
    panels.push_back(panelAtY(y, -16.0, -10.0, 2.0, 2.5));

    // Chording the corner between the side shell and the flat of bottom.
    PlatePanel chord;
    chord.corner[0] = {-30.0, y, 1.0};
    chord.corner[1] = {-28.0, y, 1.0};
    chord.corner[2] = {-28.0, y + 1.0, 0.0};
    chord.corner[3] = {-30.0, y + 1.0, 0.0};
    chord.role = PanelRole::Shell;
    panels.push_back(chord);

    // Triangles: quads with corner[2] == corner[3] and corner[0] == corner[1]
    // respectively, meeting at (12, y, 3) and nowhere else.
    PlatePanel triangleA;
    triangleA.corner[0] = {10.0, y, 1.0};
    triangleA.corner[1] = {12.0, y, 1.0};
    triangleA.corner[2] = {12.0, y, 3.0};
    triangleA.corner[3] = {12.0, y, 3.0};
    triangleA.role = PanelRole::Shell;
    panels.push_back(triangleA);
    PlatePanel triangleB;
    triangleB.corner[0] = {12.0, y, 3.0};
    triangleB.corner[1] = {12.0, y, 3.0};
    triangleB.corner[2] = {14.0, y, 5.0};
    triangleB.corner[3] = {14.0, y, 3.0};
    triangleB.role = PanelRole::Shell;
    panels.push_back(triangleB);

    // The weld grid's cells are one tolerance across and start at the origin, so
    // a corner at a whole number of metres and the same corner a nanometre below
    // it fall in *different* cells on every axis at once while being a thousandth
    // of the tolerance apart. Offset in all three so that no one axis of the
    // weld's neighbourhood search can be dropped without this seam splitting.
    constexpr double kNudge = 1e-9;
    panels.push_back(panelAtY(y, -8.0, -4.0, 1.0, 2.0));
    PlatePanel above;
    above.corner[0] = {-8.0 - kNudge, y - kNudge, 2.0 - kNudge};
    above.corner[1] = {-4.0 - kNudge, y - kNudge, 2.0 - kNudge};
    above.corner[2] = {-4.0, y, 3.0};
    above.corner[3] = {-8.0, y, 3.0};
    above.role = PanelRole::Shell;
    panels.push_back(above);

    barge.mesh.materials = {ah36Steel()};
    return barge;
}

BreachSet breachOf(const Barge& barge, const std::vector<int>& failed) {
    return breachesFromFailedPanels(barge.ship, barge.mesh, failed);
}

const Breach* named(const BreachSet& set, int a, int b) {
    for (const Breach& breach : set.breaches)
        if ((breach.opening.a == a && breach.opening.b == b) ||
            (breach.opening.a == b && breach.opening.b == a))
            return &breach;
    return nullptr;
}

// --- The winding number, which everything else stands on ----------------------

void testWindingNumber() {
    const TriMesh cube = makeBox({-1, -1, -1}, {1, 1, 1});
    expectTrue("the test cube is a closed manifold", isClosedManifold(cube));

    // The origin is the adversarial case, not the easy one. makeBox triangulates
    // its +x face as (1,2,6) and (1,6,5), whose shared diagonal runs from
    // (1,-1,-1) to (1,1,1) and therefore passes exactly through (1,0,0); the -x
    // face's diagonal passes through (-1,0,0) for the same reason. An x-axis ray
    // cast from the origin leaves through an edge at both ends, so a parity count
    // over triangle hits sees four crossings and calls the centre of the cube
    // "outside". The solid-angle sum has no ray to place.
    expectNear("winding number at the centre of a cube", meshWindingNumber(cube, {0, 0, 0}), 1.0,
               1e-12);
    expectNear("winding number just inside a face", meshWindingNumber(cube, {0.999, 0, 0}), 1.0,
               1e-9);
    expectNear("winding number just outside a face", meshWindingNumber(cube, {1.001, 0, 0}), 0.0,
               1e-9);
    expectNear("winding number far away", meshWindingNumber(cube, {50, 40, 30}), 0.0, 1e-12);

    // Sign, not just magnitude: outward winding gives +1, and a mesh wound the
    // other way gives -1. `spaceAt` tests the magnitude precisely so that an
    // inverted compartment still locates rather than silently reading as empty
    // sea, and this is the assertion that says so deliberately.
    TriMesh inverted = cube;
    for (Tri& t : inverted.tris) std::swap(t.b, t.c);
    expectNear("winding number of an inside-out mesh", meshWindingNumber(inverted, {0, 0, 0}), -1.0,
               1e-12);

    // Nothing is asserted *on* the surface, deliberately. A coplanar triangle
    // makes the triple product cancel to a signed zero, and atan2 reads +pi off
    // +0 and -pi off -0, so the answer there lands on 1 or on 0 according to
    // which way the cancellation rounded -- measured as 1 on the +x, +y and +z
    // faces of a box and 0 on the other three, and 0.25 or 0.75 on an edge
    // depending on where the tessellation put its diagonal. That is why the probe
    // is kept off surfaces, and why the overlap check asks a second time a
    // millimetre away; see testOverlappingSubdivision().

    // A non-convex mesh, where "inside" cannot be decided by half-spaces.
    const TriMesh ell = lPrism();
    expectTrue("the L prism is a closed manifold", isClosedManifold(ell));
    expectNear("the L prism holds 7 m3", integrate(ell).volume, 7.0, 1e-12);
    expectNear("winding number in the L's long arm", meshWindingNumber(ell, {3.5, 0.5, 0.5}), 1.0,
               1e-9);
    expectNear("winding number in the L's tall arm", meshWindingNumber(ell, {0.5, 3.5, 0.5}), 1.0,
               1e-9);
    expectNear("winding number in the notch, which is inside the bounding box",
               meshWindingNumber(ell, {3.5, 3.5, 0.5}), 0.0, 1e-9);
}

// --- Point location: compartment, sea, and the space that is neither ----------

void testSpaceLocation() {
    const Barge barge = makeBarge();
    const Ship& s = barge.ship;

    expectEqual("a point in the aft hold locates the aft hold", spaceAt(s, {-25, 0, 3}),
                barge.aftHold);
    expectEqual("a point in the forward hold locates the forward hold", spaceAt(s, {25, 0, 3}),
                barge.fwdHold);
    expectEqual("a point outside the shell is the sea", spaceAt(s, {-25, -12, 3}), kSea);
    expectEqual("a point above the hull is the sea", spaceAt(s, {0, 0, 12}), kSea);

    // The one that matters. Above z = 6 is inside the hull and inside no
    // compartment. Calling it sea would open every tear up there straight to the
    // ocean; calling it a compartment would flood a space that does not exist.
    expectEqual("a point inside the hull but in no compartment is neither",
                spaceAt(s, {-25, 0, 8}), kEnclosedVoid);

    // Winding, not just magnitude. A compartment mesh wound inside out is still a
    // consistently wound closed manifold and its winding number is -1, so a test
    // written as `w > 0.5` would report the space as empty sea and open every
    // tear into it straight to the ocean. (`Ship::validate()` does catch this,
    // via a negative volume -- but by then the answer has already been used.)
    Ship inverted = barge.ship;
    for (Tri& t : inverted.compartments[static_cast<std::size_t>(barge.aftHold)].mesh.tris)
        std::swap(t.b, t.c);
    expectEqual("an inside-out compartment still locates", spaceAt(inverted, {-25, 0, 3}),
                barge.aftHold);

    // And the threshold sits at the midpoint rather than at 1, so a mesh with a
    // defect in it degrades instead of flipping. Removing one face of the aft
    // hold drops its winding number by about 9% at the point tested, which is a
    // long way from 0.5 and would be the wrong side of any threshold near 1.
    Ship holed = barge.ship;
    TriMesh& mesh = holed.compartments[static_cast<std::size_t>(barge.aftHold)].mesh;
    const double face = -0.5 * kBargeBeam;
    std::vector<Tri> kept;
    for (const Tri& t : mesh.tris) {
        const bool onFace = std::abs(mesh.verts[t.a].y - face) < 1e-12 &&
                            std::abs(mesh.verts[t.b].y - face) < 1e-12 &&
                            std::abs(mesh.verts[t.c].y - face) < 1e-12;
        if (!onFace) kept.push_back(t);
    }
    expectEqual("the aft hold's outboard face is two triangles",
                static_cast<long long>(mesh.tris.size() - kept.size()), 2);
    mesh.tris = kept;
    const double defective = meshWindingNumber(mesh, {-25, 0, 3});
    expectTrue("a compartment missing a face has a winding number short of one",
               defective < 0.97);
    expectTrue("but nowhere near a half", defective > 0.8);
    expectEqual("so it still locates", spaceAt(holed, {-25, 0, 3}), barge.aftHold);

    // On the real ship the same situation arises from the *subdivision*, not from
    // a deck: amidships the engine-room boxes stop at |y| = 8 m while the shell
    // is out past 8.6 m, so the plating over the engine room mostly faces a void.
    Ship ferry = game::buildFerry();
    ferry.initialise(0.0);
    expectEqual("the ferry's starboard engine room locates by name",
                spaceAt(ferry, {6, -4, 4}), ferry.findCompartment("engine_room_s"));
    expectEqual("outboard of the ferry's engine room is an unmodelled void",
                spaceAt(ferry, {6, -8.3, 4}), kEnclosedVoid);
    expectEqual("outboard of the ferry's shell is the sea", spaceAt(ferry, {6, -12, 4}), kSea);
}

// --- Connectivity -------------------------------------------------------------

void testConnectivity() {
    const Barge barge = makeBarge();

    // A shell panel over the aft hold: sea to that hold, and nothing else. The
    // identity is asserted, not the count alone -- a routine that returned "one
    // opening between the two holds" would pass a count check.
    const BreachSet shell = breachOf(barge, {kGrid});
    expectEqual("one shell panel makes one opening", static_cast<long long>(shell.breaches.size()),
                1);
    if (shell.breaches.size() == 1) {
        const Opening& o = shell.breaches[0].opening;
        expectEqual("a shell breach starts at the sea", o.a, kSea);
        expectEqual("a shell breach ends in the compartment behind the plate", o.b, barge.aftHold);
        expectTrue("a shell breach is a breach", o.kind == OpeningKind::Breach);
        expectTrue("a shell breach is open", o.open);
        expectNear("a shell breach carries the torn plate's discharge coefficient",
                   o.dischargeCoeff, kTornPlateDischarge, 0.0);
    }
    expectTrue("a clean shell failure has nothing to report", shell.problems.empty());

    // A bulkhead panel: the two compartments it separates, by identity.
    const BreachSet bulkhead = breachOf(barge, {kBulkhead});
    expectEqual("one bulkhead panel makes one opening",
                static_cast<long long>(bulkhead.breaches.size()), 1);
    if (bulkhead.breaches.size() == 1) {
        const Opening& o = bulkhead.breaches[0].opening;
        expectEqual("a bulkhead breach joins the aft hold", o.a, barge.aftHold);
        expectEqual("a bulkhead breach joins the forward hold", o.b, barge.fwdHold);
        expectTrue("a bulkhead breach does not touch the sea", o.a != kSea && o.b != kSea);
    }

    // The bottom shell is shell too, whatever its role says.
    const BreachSet bottom = breachOf(barge, {kBottom});
    expectEqual("the flat of bottom opens to the sea",
                bottom.breaches.empty() ? -99 : bottom.breaches[0].opening.a, kSea);

    // Role is not connectivity, in both directions. The top of the hold is
    // tagged `Deck` and separates a compartment from a space the ship does not
    // model; the longitudinal bulkhead is tagged `Bulkhead` and separates the aft
    // hold from itself. Both must open nothing, and both must say why.
    const BreachSet void_ = breachOf(barge, {kHoldTopPlate});
    expectEqual("a tear into an unmodelled void opens nothing",
                static_cast<long long>(void_.breaches.size()), 0);
    expectEqual("and says so", static_cast<long long>(void_.problems.size()), 1);
    expectTrue("naming the compartment it failed to connect",
               !void_.problems.empty() &&
                   void_.problems[0].find("aft_hold") != std::string::npos);

    const BreachSet internal = breachOf(barge, {kWingBulkhead});
    expectEqual("a bulkhead that is not a compartment boundary opens nothing",
                static_cast<long long>(internal.breaches.size()), 0);
    expectTrue("and says it separates a space from itself",
               internal.problems.size() == 1 &&
                   internal.problems[0].find("from itself") != std::string::npos);

    const BreachSet degenerate = breachOf(barge, {kDegenerate});
    expectEqual("a collapsed panel opens nothing",
                static_cast<long long>(degenerate.breaches.size()), 0);
    expectTrue("and is reported rather than read as a space opening onto itself",
               degenerate.problems.size() == 1 &&
                   degenerate.problems[0].find("degenerate") != std::string::npos);

    // Out of range and repeated indices are the two ways a caller gets the input
    // wrong. Neither may cost area, and neither may pass silently.
    const BreachSet doubled = breachOf(barge, {kGrid, kGrid, kGrid});
    expectEqual("a panel listed three times still makes one opening",
                static_cast<long long>(doubled.breaches.size()), 1);
    if (doubled.breaches.size() == 1)
        expectNear("with its own area, not three times it", doubled.breaches[0].opening.area,
                   barge.mesh.panels[kGrid].area(), 0.0);
    expectEqual("and the repeats are reported", static_cast<long long>(doubled.problems.size()), 2);

    const BreachSet outOfRange = breachOf(barge, {kGrid, kPanelCount, -3});
    expectEqual("out-of-range indices are dropped, not dereferenced",
                static_cast<long long>(outOfRange.breaches.size()), 1);
    expectEqual("and reported", static_cast<long long>(outOfRange.problems.size()), 2);
    // Named, because reading past the end of the panel array can produce a
    // *different* complaint -- the garbage it reads usually has no normal, which
    // reports as a degenerate panel and keeps the count at two.
    expectTrue("as out of range, not as something read out of bounds", [&] {
        for (const std::string& problem : outOfRange.problems)
            if (problem.find("is not in a mesh") == std::string::npos) return false;
        return true;
    }());

    // The marching probe. `kChord` cuts the corner where the side meets the
    // bottom, so its centroid is half a metre inside the hull: a single probe at
    // the default 50 mm reads the aft hold on *both* sides and would drop the
    // breach with nothing worse than a note. The first two assertions establish
    // that -- without them the third would pass on a fixed probe as well and
    // prove nothing.
    const PlatePanel& chord = barge.mesh.panels[kChord];
    const Vec3 near = chord.centroid() + chord.normal() * 0.05;
    const Vec3 far = chord.centroid() - chord.normal() * 0.05;
    expectEqual("a chording panel's near probe is inside the compartment",
                spaceAt(barge.ship, near), barge.aftHold);
    expectEqual("and so is its far probe, at the same distance",
                spaceAt(barge.ship, far), barge.aftHold);
    const BreachSet chorded = breachOf(barge, {kChord});
    expectEqual("but the probe marches out until the sides differ, so it opens one hole",
                static_cast<long long>(chorded.breaches.size()), 1);
    if (chorded.breaches.size() == 1) {
        expectEqual("to the sea", chorded.breaches[0].opening.a, kSea);
        expectEqual("from the aft hold", chorded.breaches[0].opening.b, barge.aftHold);
    }

    // An empty failure set is not an error and must not invent an opening.
    const BreachSet nothing = breachOf(barge, {});
    expectEqual("no failed panels, no openings", static_cast<long long>(nothing.breaches.size()),
                0);
    expectTrue("and nothing to report", nothing.problems.empty());
}

// --- Overlapping subdivision --------------------------------------------------

void testOverlappingSubdivision() {
    // Two compartments describing the same steel, which is what the reference
    // ferry does forward. A point out there is in both, so which one a tear joins
    // is decided by declaration order and nothing else; that is a fact about the
    // ship definition and it has to reach the caller.
    Ship s;
    s.hull = makeBox({-50, -10, 0}, {50, 10, 10});
    s.deckEdgeZ = 6.0;
    s.compartments = {box("hold", {-50, -10, 0}, {50, 10, 6}),
                      box("wing_tank", {-50, 4, 0}, {50, 10, 6})};
    s.lightshipMass = 100.0 * 20.0 * 4.0 * kRhoSeawater;
    s.lightshipCog = {0, 0, 3.0};
    s.gyradii = {7.0, 28.0, 28.0};
    s.initialise(0.0);
    expectEqual("the ship definition passes validate() all the same",
                static_cast<long long>(s.validate().size()), 0);
    expectEqual("a point in the overlap resolves to the compartment declared first",
                spaceAt(s, {0, 7, 3}), s.findCompartment("hold"));

    StructuralMesh mesh;
    mesh.panels.push_back(panelAtY(10.0, -4.0, -2.0, 1.0, 3.0));  // port shell, over the tank
    mesh.materials = {ah36Steel()};
    const BreachSet set = breachesFromFailedPanels(s, mesh, {0});
    expectEqual("a tear over the overlap still opens one hole",
                static_cast<long long>(set.breaches.size()), 1);
    expectTrue("and the ambiguity is reported rather than hidden, naming both spaces", [&] {
        for (const std::string& problem : set.problems)
            if (problem.find("in the same place") != std::string::npos &&
                problem.find("hold") != std::string::npos &&
                problem.find("wing_tank") != std::string::npos)
                return true;
        return false;
    }());

    // The same complaint, on the ship the engine is validated against. This is
    // not a hypothetical: `fwd_hold_p` is authored y = 0..20 m where
    // `wing_tank_fwd_p` is 8..20 m over the same length and height, so the tank
    // is wholly inside the hold and 217 m3 of her floods twice.
    Ship ferry = game::buildFerry();
    ferry.initialise(0.0);
    expectEqual("the ferry's forward wing tank is inside her forward hold",
                spaceAt(ferry, {32, 9, 4}), ferry.findCompartment("fwd_hold_p"));
    expectEqual("while the aft wing tank, authored on the same plan, is not",
                spaceAt(ferry, {-24, 9, 4}), ferry.findCompartment("wing_tank_aft_p"));
    expectTrue("and Ship::validate() does not notice, because the total still fits",
               ferry.validate().empty());

    // Compartments that merely abut are *not* an overlap, however emphatically
    // the winding number says a point on the bulkhead between them is inside
    // both. This is not hypothetical either: the ferry's forepeak and forward
    // holds meet on the plane x = 44 m and share no volume at all, and her
    // bulkhead deck at z = 7 m lays panels straddling it, whose probe therefore
    // lands exactly on the shared face.
    const int forepeak = ferry.findCompartment("forepeak");
    const int fwdHoldS = ferry.findCompartment("fwd_hold_s");
    const Vec3 onBulkhead{44.0, -6.0, 6.95};
    // These two are the guard: if the signed-zero cancellation ever stops making
    // both compartments claim this point, the check below goes vacuous, and it
    // should fail loudly rather than quietly stop testing anything.
    expectTrue("a point on the ferry's forepeak bulkhead reads as inside the forepeak",
               std::abs(meshWindingNumber(
                   ferry.compartments[static_cast<std::size_t>(forepeak)].mesh, onBulkhead)) > 0.5);
    expectTrue("and as inside the hold on the other side of it",
               std::abs(meshWindingNumber(
                   ferry.compartments[static_cast<std::size_t>(fwdHoldS)].mesh, onBulkhead)) > 0.5);
    StructuralMesh straddling;
    straddling.panels.push_back(panelAtZ(7.0, 43.0, 45.0, -6.5, -5.5));  // centroid x = 44 exactly
    straddling.materials = {ah36Steel()};
    expectNear("and the panel that finds it really is centred on the bulkhead",
               straddling.panels[0].centroid().x, 44.0, 1e-12);
    expectTrue("but abutting compartments are not reported as sharing a place", [&] {
        for (const std::string& problem : breachesFromFailedPanels(ferry, straddling, {0}).problems)
            if (problem.find("in the same place") != std::string::npos) return false;
        return true;
    }());

    // The same, on a bulkhead running fore and aft, because a check that only
    // stepped along x would clear the case above and still be wrong here.
    Ship split;
    split.hull = makeBox({-50, -10, 0}, {50, 10, 10});
    split.deckEdgeZ = 6.0;
    split.compartments = {box("port", {-50, 0, 0}, {50, 10, 6}),
                          box("starboard", {-50, -10, 0}, {50, 0, 6})};
    split.lightshipMass = 100.0 * 20.0 * 4.0 * kRhoSeawater;
    split.lightshipCog = {0, 0, 3.0};
    split.gyradii = {7.0, 28.0, 28.0};
    split.initialise(0.0);
    StructuralMesh keel;
    keel.panels.push_back(panelAtZ(0.0, -20.0, -18.0, -1.0, 1.0));  // centroid y = 0 exactly
    keel.materials = {ah36Steel()};
    expectNear("a keel panel centred on the centreline", keel.panels[0].centroid().y, 0.0, 1e-12);
    expectTrue("does not report the two sides as sharing a place", [&] {
        for (const std::string& problem : breachesFromFailedPanels(split, keel, {0}).problems)
            if (problem.find("in the same place") != std::string::npos) return false;
        return true;
    }());
    // ...and still opens to the sea, from whichever side won the coin toss.
    const BreachSet keelBreach = breachesFromFailedPanels(split, keel, {0});
    expectEqual("while still opening one hole", static_cast<long long>(keelBreach.breaches.size()),
                1);
    if (keelBreach.breaches.size() == 1)
        expectEqual("to the sea", keelBreach.breaches[0].opening.a, kSea);
}

// --- Area, centroid and merging ----------------------------------------------

void testMergingAndArea() {
    const Barge barge = makeBarge();
    const std::vector<PlatePanel>& p = barge.mesh.panels;

    // Two panels sharing a horizontal seam. One hole, and its area is the sum of
    // theirs to the last bit -- no coefficient, no effective-area fudge.
    const BreachSet pair = breachOf(barge, {kGrid + 0, kGrid + 1});
    expectEqual("two panels sharing an edge make one opening",
                static_cast<long long>(pair.breaches.size()), 1);
    if (pair.breaches.size() == 1) {
        expectNear("the merged area is exactly the sum of the failed areas",
                   pair.breaches[0].opening.area, p[kGrid + 0].area() + p[kGrid + 1].area(), 0.0);
        expectNear("which is 4 m2 of 2 x 1 m plating", pair.breaches[0].opening.area, 4.0, 1e-12);
        // The region's own centroid: two stacked 2 x 1 m panels at z = 1..3, so
        // the middle of the tear, not the middle of either panel.
        expectNear("the opening sits at the failed region's centroid, in x",
                   pair.breaches[0].opening.pos.x, -39.0, 1e-12);
        expectNear("the opening sits at the failed region's centroid, in z",
                   pair.breaches[0].opening.pos.z, 2.0, 1e-12);
        expectEqual("and it remembers both panels",
                    static_cast<long long>(pair.breaches[0].panels.size()), 2);
    }

    // A vertical seam merges the same way, so the merge is about adjacency rather
    // than about one axis happening to line up.
    expectEqual("two panels sharing a vertical seam make one opening",
                static_cast<long long>(breachOf(barge, {kGrid + 0, kGrid + 2}).breaches.size()), 1);

    // The adversarial case. Panels 0 and 3 are diagonal neighbours: they share
    // the single corner (-38, -10, 2) and nothing else. The pinch between them
    // has zero width, so no water crosses it, and a merged opening would be
    // centred on the one point in the region where there is no hole.
    const BreachSet corner = breachOf(barge, {kGrid + 0, kGrid + 3});
    expectEqual("panels touching only at a corner stay two openings",
                static_cast<long long>(corner.breaches.size()), 2);
    if (corner.breaches.size() == 2) {
        expectNear("with the areas kept apart", corner.breaches[0].opening.area, 2.0, 1e-12);
        expectNear("and the total area still conserved",
                   corner.breaches[0].opening.area + corner.breaches[1].opening.area,
                   p[kGrid + 0].area() + p[kGrid + 3].area(), 1e-15);
        expectTrue("and the two holes in different places",
                   std::abs(corner.breaches[0].opening.pos.x - corner.breaches[1].opening.pos.x) >
                       1.0);
    }

    // The same rule has to survive a *collapsed* side, which is what the
    // generator produces where a deck or a strake runs out against the shell:
    // both of these quads are really triangles, and they meet at the corner both
    // of them collapsed to. Treating a zero-length side as an edge would put them
    // in the same bucket and fuse two holes that touch at a point.
    expectTrue("the two triangles really do share a corner",
               length(p[kTriangleA].corner[2] - p[kTriangleB].corner[0]) < 1e-12);
    expectTrue("and really are quads with a side collapsed",
               length(p[kTriangleA].corner[2] - p[kTriangleA].corner[3]) < 1e-12 &&
                   length(p[kTriangleB].corner[0] - p[kTriangleB].corner[1]) < 1e-12);
    expectEqual("two triangles meeting at a collapsed corner stay two openings",
                static_cast<long long>(breachOf(barge, {kTriangleA, kTriangleB}).breaches.size()),
                2);

    // Corners a nanometre apart are the same corner. Generated panels share their
    // corners bit for bit, but a mesh that has been through a deformation solver
    // -- which is where failed panels will come from -- will not, and a weld that
    // bucketed on a single grid cell would split this seam because its two sides
    // fall either side of a cell boundary.
    const double seamGap = length(p[kSeamBelow].corner[3] - p[kSeamAbove].corner[0]);
    expectTrue("the seam's two sides are not identical", seamGap > 0.0);
    expectTrue("but are within the weld tolerance", seamGap < 1e-6);
    const BreachSet seam = breachOf(barge, {kSeamBelow, kSeamAbove});
    expectEqual("so the seam merges into one hole",
                static_cast<long long>(seam.breaches.size()), 1);
    if (seam.breaches.size() == 1)
        expectNear("carrying both plates' area", seam.breaches[0].opening.area,
                   p[kSeamBelow].area() + p[kSeamAbove].area(), 0.0);

    // ...but a corner contact is not a firewall: bridge the two with the panel
    // that shares an edge with each and the three become one tear.
    expectEqual("a staircase joined edge to edge is one tear",
                static_cast<long long>(
                    breachOf(barge, {kGrid + 0, kGrid + 1, kGrid + 3}).breaches.size()),
                1);

    // Failures at opposite ends of the ship. Same pair of spaces, still two
    // holes: a merge rule keyed on the space pair alone would wrongly fuse them.
    const BreachSet apart = breachOf(barge, {kGrid + 0, kShellAtBulkhead});
    expectEqual("failures far apart on the same compartment stay separate",
                static_cast<long long>(apart.breaches.size()), 2);
    if (apart.breaches.size() == 2)
        expectNear("and their areas add up to the failed area",
                   apart.breaches[0].opening.area + apart.breaches[1].opening.area,
                   p[kGrid + 0].area() + p[kShellAtBulkhead].area(), 1e-15);

    // Panels that *do* share an edge but join different pairs of spaces. The
    // shell plate at the bulkhead and the bulkhead's outboard panel meet along
    // the line x = 0, y = -10, z = 1..3. One joins the sea to the aft hold, the
    // other joins the two holds; an `Opening` has two ends, so merging them would
    // have to throw one of the four away.
    const BreachSet crossing = breachOf(barge, {kBulkhead, kShellAtBulkhead});
    expectEqual("a shared edge across a change of connectivity does not merge",
                static_cast<long long>(crossing.breaches.size()), 2);
    expectTrue("the sea reaches the aft hold", named(crossing, kSea, barge.aftHold) != nullptr);
    expectTrue("and the two holds reach each other",
               named(crossing, barge.aftHold, barge.fwdHold) != nullptr);

    // The whole grid: one tear of eight panels, sixteen square metres, centred on
    // the middle of the patch.
    std::vector<int> grid;
    for (int i = 0; i < kGridPanels; ++i) grid.push_back(kGrid + i);
    const BreachSet whole = breachOf(barge, grid);
    expectEqual("a torn plate is one hole, not eight",
                static_cast<long long>(whole.breaches.size()), 1);
    if (whole.breaches.size() == 1) {
        double sum = 0;
        for (int i : grid) sum += p[static_cast<std::size_t>(i)].area();
        expectNear("area is conserved across the whole tear", whole.breaches[0].opening.area, sum,
                   0.0);
        expectNear("and comes to 16 m2", whole.breaches[0].opening.area, 16.0, 1e-12);
        expectNear("the tear's centroid in x", whole.breaches[0].opening.pos.x, -36.0, 1e-12);
        expectNear("the tear's centroid in z", whole.breaches[0].opening.pos.z, 2.0, 1e-12);
        expectEqual("and every panel is accounted for",
                    static_cast<long long>(whole.breaches[0].panels.size()), kGridPanels);
    }

    // The answer must not depend on the order the failures were listed in: a
    // fracture model will hand them over in whatever order its elements failed.
    std::vector<int> shuffled{kGrid + 5, kGrid + 1, kGrid + 7, kGrid + 0,
                              kGrid + 4, kGrid + 3, kGrid + 6, kGrid + 2};
    const BreachSet reordered = breachOf(barge, shuffled);
    expectEqual("listing order does not change the number of openings",
                static_cast<long long>(reordered.breaches.size()),
                static_cast<long long>(whole.breaches.size()));
    if (reordered.breaches.size() == 1 && whole.breaches.size() == 1) {
        expectNear("nor the area", reordered.breaches[0].opening.area,
                   whole.breaches[0].opening.area, 0.0);
        expectNear("nor the position", reordered.breaches[0].opening.pos.x,
                   whole.breaches[0].opening.pos.x, 0.0);
        expectTrue("nor the name", reordered.breaches[0].opening.name ==
                                       whole.breaches[0].opening.name);
    }

    // A tear of two unequal plates: the merged position is the *area-weighted*
    // centroid, which is what sets the head and therefore the flow. 6 m2 centred
    // at z = 1.5 under 3 m2 centred at z = 2.25 gives 1.75 m; the unweighted mean
    // of the two panel centres is 1.875 m, so the two rules are 125 mm apart and
    // this assertion can tell them apart.
    const BreachSet weighted = breachOf(barge, {kWide, kNarrow});
    expectEqual("two unequal plates sharing a seam are one hole",
                static_cast<long long>(weighted.breaches.size()), 1);
    if (weighted.breaches.size() == 1) {
        expectNear("of 9 m2", weighted.breaches[0].opening.area,
                   p[kWide].area() + p[kNarrow].area(), 0.0);
        expectNear("sitting at the area-weighted centroid",
                   weighted.breaches[0].opening.pos.z, 1.75, 1e-12);
        expectTrue("which is not the mean of the panel centres",
                   std::abs(1.75 - 0.5 * (p[kWide].centroid().z + p[kNarrow].centroid().z)) > 0.1);
    }
}

// --- The flooding rate --------------------------------------------------------
//
// The milestone's actual claim. A barge big enough that the head at the hole
// barely moves while the hole is measured: 200 x 40 m of waterplane against a
// hold 20 m long, so thirty seconds of flooding sinks her 7 mm and raises the
// water inside by 70 mm, against an opening 2 m under and 1.5 m up.

constexpr double kRigLength = 200.0, kRigBeam = 40.0, kRigDepth = 12.0;
constexpr double kRigDraft = 4.0;
constexpr double kHoleTop = 2.5, kHoleBottom = 1.5;  // body frame; centre at z = 2
constexpr double kPanelWidth = 0.5;

struct FloodRig {
    Ship ship;
    StructuralMesh mesh;
};

// `panels` rectangles of 0.5 x 1.0 m side by side along the shell, all at the
// same height, so failing more of them changes the area and nothing else.
FloodRig makeFloodRig(int panels) {
    FloodRig rig;
    Ship& s = rig.ship;
    s.hull = makeBox({-0.5 * kRigLength, -0.5 * kRigBeam, 0},
                     {0.5 * kRigLength, 0.5 * kRigBeam, kRigDepth});
    s.deckEdgeZ = 10.0;
    s.compartments = {box("hold", {-10, -0.5 * kRigBeam, 0}, {10, 0.5 * kRigBeam, 10})};
    s.lightshipMass = kRigLength * kRigBeam * kRigDraft * kRhoSeawater;
    s.lightshipCog = {0, 0, 3.0};
    s.gyradii = {14.0, 58.0, 58.0};
    s.initialise(0.0);

    for (int i = 0; i < panels; ++i)
        rig.mesh.panels.push_back(panelAtY(-0.5 * kRigBeam, i * kPanelWidth, (i + 1) * kPanelWidth,
                                           kHoleBottom, kHoleTop));
    rig.mesh.materials = {ah36Steel()};
    return rig;
}

// Settle, then open the hole and integrate for `seconds`. Returns the water
// taken; `head` comes back as the depth of the opening below the sea at the
// moment it opened.
double floodThrough(FloodRig& rig, const std::vector<Opening>& openings, double seconds,
                    double* head) {
    const double dt = 0.01;
    for (int i = 0; i < 6000; ++i) rig.ship.step(dt, 0.0);  // 60 s to settle
    if (head != nullptr && !openings.empty()) {
        const Vec3 world = rig.ship.state.orientation.toMat3() * openings[0].pos +
                           rig.ship.state.position;
        *head = -world.z;  // the sea is at z = 0
    }
    for (const Opening& o : openings) rig.ship.openings.push_back(o);
    const int steps = static_cast<int>(seconds / dt);
    for (int i = 0; i < steps; ++i) rig.ship.step(dt, 0.0);
    return rig.ship.compartments[0].waterVolume;
}

void testFloodingRateFollowsArea() {
    const double seconds = 30.0;

    // --- One panel ------------------------------------------------------------
    FloodRig one = makeFloodRig(1);
    const BreachSet set1 = breachesFromFailedPanels(one.ship, one.mesh, {0});
    expectEqual("the rig's failed panel makes one opening",
                static_cast<long long>(set1.breaches.size()), 1);
    if (set1.breaches.size() != 1) return;
    const Opening opening1 = set1.breaches[0].opening;
    expectNear("of half a square metre", opening1.area, kPanelWidth * (kHoleTop - kHoleBottom),
               0.0);

    double head1 = 0;
    const double volume1 = floodThrough(one, {opening1}, seconds, &head1);

    // Torricelli through the orifice: q = Cd A sqrt(2 g h), integrated over the
    // interval. `docs/02-simulation.md` §1 and README both state the law; nothing
    // in this expression comes from the solver.
    expectNear("the rig settles with the hole 2 m under", head1, 2.0, 1e-3);
    const double rate1 = opening1.dischargeCoeff * opening1.area * std::sqrt(2.0 * kGravity * head1);
    // The head is not quite constant: she sinks 7 mm over the interval, which
    // lifts sqrt(h) by 0.18% at the end and about 0.09% on average. Anything
    // larger than that is a discrepancy rather than drift.
    expectNear("water taken matches Cd A sqrt(2 g h) over the interval", volume1,
               rate1 * seconds, 0.005 * rate1 * seconds);
    expectTrue("and the run was not vacuous", volume1 > 50.0);

    // --- Twice the failed plating --------------------------------------------
    FloodRig two = makeFloodRig(2);
    const BreachSet set2 = breachesFromFailedPanels(two.ship, two.mesh, {0, 1});
    expectEqual("two adjacent failed panels still make one opening",
                static_cast<long long>(set2.breaches.size()), 1);
    if (set2.breaches.size() != 1) return;
    const Opening opening2 = set2.breaches[0].opening;
    expectNear("of twice the area", opening2.area, 2.0 * opening1.area, 1e-15);
    expectNear("at the same depth", opening2.pos.z, opening1.pos.z, 1e-15);

    double head2 = 0;
    const double volume2 = floodThrough(two, {opening2}, seconds, &head2);
    const double rate2 = opening2.dischargeCoeff * opening2.area * std::sqrt(2.0 * kGravity * head2);
    expectNear("the wider hole also matches the orifice law", volume2, rate2 * seconds,
               0.005 * rate2 * seconds);

    // This is the milestone sentence, as a number: the hole's own area sets the
    // rate. The ratio runs a shade over two because the wider hole sinks the ship
    // faster and so raises its own head.
    expectNear("doubling the failed area doubles the water taken", volume2 / volume1, 2.0, 0.01);

    // --- Nothing fails, nothing floods ---------------------------------------
    FloodRig dry = makeFloodRig(1);
    const BreachSet none = breachesFromFailedPanels(dry.ship, dry.mesh, {});
    expectEqual("an empty failure set makes no openings",
                static_cast<long long>(none.breaches.size()), 0);
    expectEqual("and applying it adds nothing to the network",
                static_cast<long long>(applyBreaches(dry.ship, none)), 0);
    expectEqual("leaving the ship with no openings at all",
                static_cast<long long>(dry.ship.openings.size()), 0);
    const double volumeNone = floodThrough(dry, {}, seconds, nullptr);
    expectNear("and she stays exactly dry", volumeNone, 0.0, 0.0);
    // The guard: "stays dry" is worthless on a ship that cannot flood, so the
    // identical hull with one panel failed must take a great deal of water.
    expectTrue("on a ship that floods heavily when a panel does fail", volume1 > 50.0);
}

// --- Indistinguishable from an authored opening -------------------------------

void testAuthoredTwin() {
    FloodRig produced = makeFloodRig(2);
    const BreachSet set = breachesFromFailedPanels(produced.ship, produced.mesh, {0, 1});
    if (set.breaches.size() != 1) {
        expectTrue("the twin test needs exactly one opening", false);
        return;
    }

    // The same hole, written the way `game/prototype/ferry.cpp` writes one.
    Opening authored;
    authored.name = "breach_by_hand";
    authored.a = kSea;
    authored.b = 0;
    authored.pos = {0.5, -0.5 * kRigBeam, 2.0};
    authored.area = 1.0;
    authored.dischargeCoeff = kTornPlateDischarge;
    authored.kind = OpeningKind::Breach;
    authored.open = true;

    const Opening& grown = set.breaches[0].opening;
    expectNear("the produced opening is where the authored one is (x)", grown.pos.x, authored.pos.x,
               0.0);
    expectNear("the produced opening is where the authored one is (y)", grown.pos.y, authored.pos.y,
               0.0);
    expectNear("the produced opening is where the authored one is (z)", grown.pos.z, authored.pos.z,
               0.0);
    expectNear("and is the same size", grown.area, authored.area, 0.0);

    FloodRig byHand = makeFloodRig(2);
    double ignored = 0;
    const double volumeGrown = floodThrough(produced, {grown}, 20.0, &ignored);
    const double volumeAuthored = floodThrough(byHand, {authored}, 20.0, &ignored);
    // Bit-for-bit: the solver sees two identical structs, so anything other than
    // equality would mean the produced opening carries a field the authored one
    // does not.
    expectNear("a produced opening floods the ship exactly as an authored one does",
               volumeGrown, volumeAuthored, 0.0);
    expectTrue("with something to compare", volumeGrown > 10.0);

    // And the ship's own consistency check must not notice the difference.
    Ship ferry = game::buildFerry();
    ferry.initialise(0.0);
    const std::size_t before = ferry.validate().size();
    std::vector<std::string> meshProblems;
    const StructuralMesh mesh = makeStructuralMesh(ferry.hull, ferryScantlings(), &meshProblems);
    std::vector<int> failed;
    for (std::size_t i = 0; i < mesh.panels.size(); ++i)
        if (mesh.panels[i].role == PanelRole::Bulkhead) failed.push_back(static_cast<int>(i));
    const BreachSet ferrySet = breachesFromFailedPanels(ferry, mesh, failed);
    expectTrue("failing every bulkhead on the ferry opens something",
               !ferrySet.breaches.empty());
    applyBreaches(ferry, ferrySet);
    expectEqual("and the ship still validates", static_cast<long long>(ferry.validate().size()),
                static_cast<long long>(before));
}

// --- The reference ferry ------------------------------------------------------
//
// Built once: the structural mesh is 8 900 panels and probing all of them costs
// about 150 ms, which is worth paying once and not five times.

struct FerryDamage {
    Ship ship;
    StructuralMesh mesh;
    BreachSet everything;  // every panel failed: the complete damage map
};

const FerryDamage& ferryDamage() {
    static const FerryDamage damage = [] {
        FerryDamage d;
        d.ship = game::buildFerry();
        d.ship.initialise(0.0);
        std::vector<std::string> problems;
        d.mesh = makeStructuralMesh(d.ship.hull, ferryScantlings(), &problems);
        std::vector<int> all(d.mesh.panels.size());
        for (std::size_t i = 0; i < all.size(); ++i) all[i] = static_cast<int>(i);
        d.everything = breachesFromFailedPanels(d.ship, d.mesh, all);
        return d;
    }();
    return damage;
}

void testFerryConnectivity() {
    const FerryDamage& damage = ferryDamage();
    const Ship& ship = damage.ship;
    const int engineRoomS = ship.findCompartment("engine_room_s");
    const int aftHoldS = ship.findCompartment("aft_hold_s");
    const int accommodation = ship.findCompartment("accommodation");
    const int vehicleDeck = ship.findCompartment("vehicle_deck");

    // Every shell panel that reaches the starboard engine room, failed together,
    // is one hole between the sea and that engine room. The panels come from the
    // complete map rather than from an index written down here, because a panel
    // index is a property of the generator's loop order and pinning one would
    // make this fail the next time a strake is added.
    const Breach* erShell = named(damage.everything, kSea, engineRoomS);
    expectTrue("the ferry has shell plating against the starboard engine room",
               erShell != nullptr);
    if (erShell != nullptr) {
        const BreachSet alone = breachesFromFailedPanels(ship, damage.mesh, erShell->panels);
        expectEqual("failing it opens exactly one hole",
                    static_cast<long long>(alone.breaches.size()), 1);
        if (alone.breaches.size() == 1) {
            expectEqual("between the sea", alone.breaches[0].opening.a, kSea);
            expectEqual("and the starboard engine room", alone.breaches[0].opening.b, engineRoomS);
        }
        // One panel out of that patch, on its own, still joins the same two
        // spaces: the connectivity is a property of the plate, not of the patch.
        const BreachSet single =
            breachesFromFailedPanels(ship, damage.mesh, {erShell->panels.front()});
        expectEqual("and one plate of it opens the same two spaces",
                    static_cast<long long>(single.breaches.size()), 1);
        if (single.breaches.size() == 1) {
            expectEqual("sea", single.breaches[0].opening.a, kSea);
            expectEqual("to the starboard engine room", single.breaches[0].opening.b, engineRoomS);
        }
        expectTrue("every panel in the patch is shell plating", [&] {
            for (int i : erShell->panels)
                if (damage.mesh.panels[static_cast<std::size_t>(i)].role != PanelRole::Shell)
                    return false;
            return true;
        }());
    }

    // The watertight bulkhead at x = -8 m: engine room to aft hold, on the same
    // side. Getting this wrong by one compartment is the difference between
    // flooding the space next door and flooding the space across the ship.
    const Breach* bulkhead = named(damage.everything, engineRoomS, aftHoldS);
    expectTrue("the ferry's aft engine-room bulkhead separates two compartments",
               bulkhead != nullptr);
    if (bulkhead != nullptr) {
        expectNear("and stands at x = -8 m", bulkhead->opening.pos.x, -8.0, 1e-6);
        expectTrue("on the starboard side", bulkhead->opening.pos.y < 0.0);
        const BreachSet one =
            breachesFromFailedPanels(ship, damage.mesh, {bulkhead->panels.front()});
        expectEqual("one bulkhead plate makes one opening",
                    static_cast<long long>(one.breaches.size()), 1);
        if (one.breaches.size() == 1) {
            expectEqual("joining the starboard engine room", one.breaches[0].opening.a,
                        std::min(engineRoomS, aftHoldS));
            expectEqual("to the starboard aft hold", one.breaches[0].opening.b,
                        std::max(engineRoomS, aftHoldS));
        }
    }

    // Role really is not connectivity, and the weather deck proves it in the
    // direction that would embarrass a label: it is `PanelRole::Deck`, it has the
    // sky above it, and a hole in it opens the accommodation to the sea.
    int weatherPanel = -1;
    for (std::size_t i = 0; i < damage.mesh.panels.size() && weatherPanel < 0; ++i) {
        const PlatePanel& panel = damage.mesh.panels[i];
        if (panel.role != PanelRole::Deck) continue;
        const Vec3 c = panel.centroid();
        if (std::abs(c.z - 15.0) < 1e-6 && std::abs(c.x) < 10.0 && std::abs(c.y) < 4.0)
            weatherPanel = static_cast<int>(i);
    }
    expectTrue("the ferry has weather deck plating amidships", weatherPanel >= 0);
    if (weatherPanel >= 0) {
        const BreachSet deck = breachesFromFailedPanels(ship, damage.mesh, {weatherPanel});
        expectEqual("a hole in the weather deck makes one opening",
                    static_cast<long long>(deck.breaches.size()), 1);
        if (deck.breaches.size() == 1) {
            expectEqual("to the sea", deck.breaches[0].opening.a, kSea);
            expectEqual("from the accommodation", deck.breaches[0].opening.b, accommodation);
        }
    }
    // ...while the vehicle deck head, tagged the same way, joins two
    // compartments and touches no sea at all.
    expectTrue("the vehicle deck head joins two compartments",
               named(damage.everything, vehicleDeck, accommodation) != nullptr);

    // The wing bulkhead runs at |y| = 6 m, inside the holds rather than on a
    // compartment boundary, so tearing it opens nothing at all.
    int wingPanel = -1;
    for (std::size_t i = 0; i < damage.mesh.panels.size() && wingPanel < 0; ++i) {
        const PlatePanel& panel = damage.mesh.panels[i];
        if (panel.role != PanelRole::Bulkhead) continue;
        const Vec3 c = panel.centroid();
        if (std::abs(c.y - 6.0) < 1e-6 && c.x > 0 && c.x < 30 && c.z > 2.5 && c.z < 6.5)
            wingPanel = static_cast<int>(i);
    }
    expectTrue("the ferry has a longitudinal wing bulkhead at y = 6 m", wingPanel >= 0);
    if (wingPanel >= 0) {
        const BreachSet wing = breachesFromFailedPanels(ship, damage.mesh, {wingPanel});
        expectEqual("tearing a bulkhead that is not a compartment boundary opens nothing",
                    static_cast<long long>(wing.breaches.size()), 0);
        expectTrue("and it is reported rather than dropped",
                   wing.problems.size() == 1 &&
                       wing.problems[0].find("from itself") != std::string::npos);
    }

    // Amidships the engine-room boxes stop short of the shell, so most of the
    // plating over them faces a void. That must open nothing and say so -- the
    // dangerous alternative is calling the void "sea" and flooding a compartment
    // that the tear does not actually reach.
    int voidPanel = -1;
    for (std::size_t i = 0; i < damage.mesh.panels.size() && voidPanel < 0; ++i) {
        const PlatePanel& panel = damage.mesh.panels[i];
        if (panel.role != PanelRole::Shell) continue;
        const Vec3 c = panel.centroid();
        if (c.x > 0 && c.x < 12 && c.y < -8.2 && c.z > 3.5 && c.z < 5.0)
            voidPanel = static_cast<int>(i);
    }
    expectTrue("the ferry has side plating amidships outboard of the engine room",
               voidPanel >= 0);
    if (voidPanel >= 0) {
        const BreachSet gap = breachesFromFailedPanels(ship, damage.mesh, {voidPanel});
        expectEqual("plating over an unmodelled void opens nothing",
                    static_cast<long long>(gap.breaches.size()), 0);
        expectTrue("and names the void",
                   gap.problems.size() == 1 &&
                       gap.problems[0].find("no compartment describes") != std::string::npos);
    }

    // The reason the probe marches, on the geometry that forced it. A shell panel
    // is a flat chord across a curved surface, and where a girth band spans the
    // crease at the turn of the bilge its centroid ends up well inside the hull.
    // Find the worst such panel by asking the hull itself, then require that it
    // still opens to the sea.
    int buried = -1;
    double deepest = 0;
    for (std::size_t i = 0; i < damage.mesh.panels.size(); ++i) {
        const PlatePanel& panel = damage.mesh.panels[i];
        if (panel.role != PanelRole::Shell) continue;
        const Vec3 c = panel.centroid();
        const Vec3 n = panel.normal();
        // How far out along the normal both sides stay inside the hull.
        double inside = 0;
        for (double d = 0.02; d <= 0.40; d += 0.02) {
            if (std::abs(meshWindingNumber(ship.hull, c + n * d)) < 0.5) break;
            if (std::abs(meshWindingNumber(ship.hull, c - n * d)) < 0.5) break;
            inside = d;
        }
        if (inside > deepest) {
            deepest = inside;
            buried = static_cast<int>(i);
        }
    }
    expectTrue("the ferry has shell plating whose centroid is buried in the hull",
               deepest > 0.10);
    if (buried >= 0) {
        const BreachSet deep = breachesFromFailedPanels(ship, damage.mesh, {buried});
        expectEqual("and a fixed 50 mm probe would have dropped it",
                    spaceAt(ship, damage.mesh.panels[static_cast<std::size_t>(buried)].centroid() +
                                      damage.mesh.panels[static_cast<std::size_t>(buried)].normal() *
                                          0.05),
                    spaceAt(ship, damage.mesh.panels[static_cast<std::size_t>(buried)].centroid() -
                                      damage.mesh.panels[static_cast<std::size_t>(buried)].normal() *
                                          0.05));
        expectEqual("the marching probe opens it", static_cast<long long>(deep.breaches.size()), 1);
        if (deep.breaches.size() == 1)
            expectEqual("to the sea", deep.breaches[0].opening.a, kSea);
    }

    // The complete map, as a sanity check on the whole ship at once: every
    // opening joins two different spaces, neither of which is a void, both of
    // which exist, and its area is the sum of its panels'.
    const int compartments = static_cast<int>(ship.compartments.size());
    bool wellFormed = true, areaHolds = true;
    double totalMerged = 0, totalPanels = 0;
    for (const Breach& breach : damage.everything.breaches) {
        const Opening& o = breach.opening;
        if (o.a == o.b || o.a == kEnclosedVoid || o.b == kEnclosedVoid) wellFormed = false;
        if (o.a < kSea || o.a >= compartments || o.b < kSea || o.b >= compartments)
            wellFormed = false;
        if (o.area <= 0 || breach.panels.empty()) wellFormed = false;
        double sum = 0;
        for (int i : breach.panels) sum += damage.mesh.panels[static_cast<std::size_t>(i)].area();
        if (std::abs(sum - o.area) > 1e-9 * std::max(sum, 1.0)) areaHolds = false;
        totalMerged += o.area;
        totalPanels += sum;
    }
    expectTrue("every opening on the ferry joins two real, different spaces", wellFormed);
    expectTrue("and carries exactly its panels' area", areaHolds);
    expectNear("with no area created or lost in the merge", totalMerged, totalPanels,
               1e-9 * totalPanels);
    expectTrue("the whole-ship map is not empty", damage.everything.breaches.size() > 30);
}

// The damage map is worth keeping after a run: it is the one artefact that says
// what a ship's subdivision and its structure actually imply about each other,
// and re-deriving it by hand from a failing assertion is miserable.
void testFerryDamageInventory() {
    const FerryDamage& damage = ferryDamage();
    const std::string path = testing::scratchDir() + "ferry_damage_map.txt";
    std::FILE* file = std::fopen(path.c_str(), "w");
    expectTrue("the damage inventory can be written", file != nullptr);
    if (file == nullptr) return;
    for (const Breach& breach : damage.everything.breaches)
        std::fprintf(file, "%-44s %10.4f %4zu\n", breach.opening.name.c_str(), breach.opening.area,
                     breach.panels.size());
    for (const std::string& problem : damage.everything.problems)
        std::fprintf(file, "# %s\n", problem.c_str());
    std::fclose(file);

    // Read it back and check it against the set it came from, so the file is
    // evidence rather than decoration.
    file = std::fopen(path.c_str(), "r");
    expectTrue("and read back", file != nullptr);
    if (file == nullptr) return;
    char name[256];
    double area = 0;
    unsigned long panels = 0;
    int rows = 0;
    double areaSum = 0;
    unsigned long panelSum = 0;
    while (std::fscanf(file, "%255s %lf %lu", name, &area, &panels) == 3) {
        ++rows;
        areaSum += area;
        panelSum += panels;
    }
    std::fclose(file);
    expectEqual("every breach reached the inventory", rows,
                static_cast<long long>(damage.everything.breaches.size()));
    expectNear("with the total area intact", areaSum, damage.everything.totalArea(), 0.5);
    long long panelTotal = 0;
    for (const Breach& breach : damage.everything.breaches)
        panelTotal += static_cast<long long>(breach.panels.size());
    expectEqual("and every panel accounted for", static_cast<long long>(panelSum), panelTotal);
}

// --- The ferry, rammed --------------------------------------------------------

void testFerryFloodsThroughItsOwnStructure() {
    const FerryDamage& damage = ferryDamage();
    const int engineRoomS = damage.ship.findCompartment("engine_room_s");
    const Breach* patch = named(damage.everything, kSea, engineRoomS);
    if (patch == nullptr) {
        expectTrue("the ferry needs shell plating on the engine room to be rammed", false);
        return;
    }

    // Half the patch, so the hole is a plausible collision rather than the whole
    // side of the ship.
    std::vector<int> failed(patch->panels.begin(),
                            patch->panels.begin() + (patch->panels.size() + 1) / 2);

    // The ferry's own casualty is shut, so that anything that floods below got in
    // through structure rather than through an authored hole. Her engine-room air
    // pipe stays open: a sealed space stops flooding once its air balances the
    // head, which is real physics tested elsewhere and would only mask this.
    Ship control = game::buildFerry();
    control.initialise(0.0);
    for (Opening& o : control.openings) o.open = (o.name == "vent_er_s");
    Ship damaged = control;

    const BreachSet set = breachesFromFailedPanels(damaged, damage.mesh, failed);
    expectTrue("the collision opens at least one hole", !set.breaches.empty());
    expectTrue("of a few square metres", set.totalArea() > 1.0 && set.totalArea() < 40.0);
    const std::size_t hole = damaged.openings.size();
    expectEqual("added to the flooding network", static_cast<long long>(applyBreaches(damaged, set)),
                static_cast<long long>(set.breaches.size()));

    // The rate the hole's own area determines, at the instant it opens. The
    // engine room is dry and at atmospheric, so the pressure difference across
    // the tear is rho*g*h and nothing else, and the expected flow is the orifice
    // law evaluated on the ship's actual attitude. Everything on the right-hand
    // side comes from the breach and from where the ship is floating.
    const Opening& o = damaged.openings[hole];
    const Vec3 world = damaged.state.orientation.toMat3() * o.pos + damaged.state.position;
    const double head = -world.z;
    expectTrue("the tear is well below the waterline", head > 1.0);
    damaged.step(1e-4, 0.0);
    expectNear("the flow through a torn plate is Cd A sqrt(2 g h)",
               damaged.openings[hole].lastFlow,
               o.dischargeCoeff * o.area * std::sqrt(2.0 * kGravity * head),
               1e-6 * o.dischargeCoeff * o.area * std::sqrt(2.0 * kGravity * head));
    expectTrue("and it is water, not air", damaged.openings[hole].lastFlowWasWater);

    const double dt = 0.02;
    for (int i = 0; i < 3000; ++i) {  // 60 s
        control.step(dt, 0.0);
        damaged.step(dt, 0.0);
    }
    expectNear("with the authored casualty shut the ferry stays dry",
               control.totalFloodwaterMass(), 0.0, 1.0);
    const double flooded = damaged.compartments[static_cast<std::size_t>(engineRoomS)].waterVolume;
    expectTrue("but a hole in her own plating floods the engine room", flooded > 300.0);
    expectTrue("and she settles by the amount that water weighs",
               damaged.diagnostics(0.0).draftMidship > control.diagnostics(0.0).draftMidship + 0.05);
}

}  // namespace

void runBreachTests() {
    std::printf("\n--- structural failure to flooding openings ---\n");
    testWindingNumber();
    testSpaceLocation();
    testConnectivity();
    testOverlappingSubdivision();
    testMergingAndArea();
    testFloodingRateFollowsArea();
    testAuthoredTwin();
    testFerryConnectivity();
    testFerryDamageInventory();
    testFerryFloodsThroughItsOwnStructure();
}
