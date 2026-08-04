// SPDX-License-Identifier: MIT
//
// Validation of the ship definition format.
//
// Two questions, and the second is the one that matters.
//
// **Does the format describe the ferry?** `game/prototype/ferry.cpp` stays in the
// tree precisely so this can be asked: `ships/ferry.ship` is loaded and held
// against the compiled ship, hull triangle by hull triangle, compartment by
// compartment, and then both are floated and their hydrostatics compared. A
// format that dropped permeability, or the vented flag, or an opening's
// discharge coefficient would still produce a ship, and it would still float --
// it would just be a different ship. Comparing the *outputs* rather than the
// fields is what makes an omission visible.
//
// **Does it capture what matters?** Equal hydrostatics are a static claim. So
// both ships are then run through the damage-control scenario in lockstep and
// their trajectories compared, with a deliberately mis-copied file as the
// negative control -- because a comparison that cannot fail proves nothing, and
// this suite would otherwise be asserting that two ships did nothing identically.
//
// The malformed-input half is aimed at CLAUDE.md's `World::load`, which failed
// open and left a half-built world. The instrument that caught it -- every
// truncation of a valid file -- is pointed here at a compact test barge, which
// exercises every state of the parser at a thousandth of the ferry's parse cost.
#include "engine/sim/shipfile.hpp"
#include "game/prototype/ferry.hpp"
#include "harness.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using sim::Compartment;
using sim::Diagnostics;
using sim::Opening;
using sim::Pump;
using sim::Ship;
using sim::ShipDefinition;
using sim::Vec3;
using testing::expectEqual;
using testing::expectNear;
using testing::expectTrue;

namespace {

std::string label(const std::string& what, double value) {
    std::ostringstream text;
    text << what << " (" << value << ")";
    return text.str();
}

// --- Tolerances ---------------------------------------------------------------
//
// `ships/ferry.ship` states half-breadths to six decimals, which is a micron and
// four orders of magnitude finer than a shipyard measures. That rounding is the
// *only* difference between the two ships, and everything below is the measured
// consequence of it, tightened to the nearest round number:
//
//   hull vertices          5.0e-7 m   -- exactly the 1 um quantisation
//   displacement           3.3e-9 rel
//   compartment volume     2.7e-8 rel -- the clip amplifies; wing_tank_aft_s
//   GM, KB, GZ             4.1e-8 m
//   draft, freeboard       8.0e-15 m  -- the sinkage solve converges to the same root
//
// These are not the format's resolution. The format carries whatever a double
// carries: testFormatIsLosslessAtFullPrecision() writes seventeen significant
// figures and gets the mesh back bit for bit. They are this *file's* resolution,
// chosen for a table a human can read and diff.
constexpr double kVertexTol = 1e-6;      // m
constexpr double kVolumeRelTol = 1e-7;   // dimensionless
constexpr double kLengthTol = 1e-6;      // m: draft, GM, KB, GZ, freeboard
constexpr double kMassRelTol = 1e-8;     // dimensionless

std::string shipDir() {
#ifdef SHIPSIM_SHIP_DIR
    return std::string(SHIPSIM_SHIP_DIR) + "/";
#else
    return std::string("ships/");
#endif
}

std::string readFile(const std::string& path) {
    std::ifstream file(path, std::ios::binary);
    std::ostringstream contents;
    contents << file.rdbuf();
    return contents.str();
}

// A compact but complete ship: every block type, every required key, a hull that
// closes, two compartments that survive the clip, an opening to the sea and one
// between spaces, and a pump. Small enough that every byte truncation of it can
// be parsed inside a unit test.
const char* const kBarge = R"(# Test barge.
ship_format 1

ship barge
    deck_edge_z 3.0
    lightship_mass 1.2e6
    lightship_cog 0.0 0.0 2.0
    gyradii 3.0 8.0 8.0

hull
    waterlines 0.0 2.0 4.0
    station -10.0  4.0 5.0 5.0
    station   0.0  5.0 5.0 5.0
    station  10.0  4.0 5.0 5.0

compartment hold_fwd
    box 0 -6 0   12 6 3
    permeability 0.95
    vented false

compartment hold_aft
    box -12 -6 0   0 6 3
    permeability 0.90
    vented true

opening breach
    between sea hold_fwd
    at 4 -5 1
    area 0.5
    discharge 0.6
    kind breach
    open true

opening door
    between hold_fwd hold_aft
    at 0 0 0.5
    area 2.0
    discharge 0.7
    kind door
    open false

pump bilge
    drains hold_fwd
    capacity 0.05
    max_head 20
    on false
)";

// Replace the first occurrence of `from`, asserting there was one -- a scripted
// edit that silently matched nothing has removed a real assertion in this repo
// before, and a mutation test that failed to mutate would pass for free.
std::string mutate(const std::string& text, const std::string& from, const std::string& to) {
    const std::size_t at = text.find(from);
    expectTrue("mutation '" + from + "' has something to replace", at != std::string::npos);
    if (at == std::string::npos) return text;
    return text.substr(0, at) + to + text.substr(at + from.size());
}

// A definition that is *not* default-constructed, so "the loader left `out`
// alone" is distinguishable from "the loader cleared `out`".
ShipDefinition sentinel() {
    ShipDefinition marker;
    marker.name = "untouched";
    marker.ship.lightshipMass = 12345.0;
    marker.ship.deckEdgeZ = -99.0;
    return marker;
}

bool isSentinel(const ShipDefinition& d) {
    return d.name == "untouched" && d.ship.lightshipMass == 12345.0 && d.ship.deckEdgeZ == -99.0 &&
           d.ship.hull.tris.empty() && d.ship.compartments.empty() && d.ship.openings.empty() &&
           d.ship.pumps.empty();
}

// --- The load-bearing comparison ----------------------------------------------

struct Divergence {
    double vertex = 0;      // m
    double volumeRel = 0;   // dimensionless
    double massRel = 0;     // dimensionless
    double length = 0;      // m -- draft, GM, KB, freeboard, GZ
    bool   structural = false;  // counts, names, indices, flags, kinds
};

// Everything that is not a floating-point quantity: the two ships must agree
// exactly, and a mismatch here means the format lost a field rather than a digit.
Divergence compare(const Ship& a, const Ship& b) {
    Divergence d;
    const auto structuralUnless = [&](bool ok) { d.structural = d.structural || !ok; };

    structuralUnless(a.hull.verts.size() == b.hull.verts.size());
    structuralUnless(a.hull.tris.size() == b.hull.tris.size());
    if (a.hull.verts.size() == b.hull.verts.size())
        for (std::size_t i = 0; i < a.hull.verts.size(); ++i) {
            const Vec3 delta = a.hull.verts[i] - b.hull.verts[i];
            d.vertex = std::max({d.vertex, std::abs(delta.x), std::abs(delta.y),
                                 std::abs(delta.z)});
        }
    // Index-for-index, so the winding is compared too. A hull with the same
    // vertices and a flipped triangle integrates to a different displacement,
    // and CLAUDE.md's defect table has exactly that costing 40%.
    if (a.hull.tris.size() == b.hull.tris.size())
        for (std::size_t i = 0; i < a.hull.tris.size(); ++i)
            structuralUnless(a.hull.tris[i].a == b.hull.tris[i].a &&
                             a.hull.tris[i].b == b.hull.tris[i].b &&
                             a.hull.tris[i].c == b.hull.tris[i].c);

    structuralUnless(a.compartments.size() == b.compartments.size());
    if (a.compartments.size() == b.compartments.size())
        for (std::size_t i = 0; i < a.compartments.size(); ++i) {
            const Compartment& x = a.compartments[i];
            const Compartment& y = b.compartments[i];
            structuralUnless(x.name == y.name);
            structuralUnless(x.permeability == y.permeability);
            structuralUnless(x.ventedToAtmosphere == y.ventedToAtmosphere);
            if (x.grossVolume > 0)
                d.volumeRel = std::max(d.volumeRel,
                                       std::abs(x.grossVolume - y.grossVolume) / x.grossVolume);
        }

    structuralUnless(a.openings.size() == b.openings.size());
    if (a.openings.size() == b.openings.size())
        for (std::size_t i = 0; i < a.openings.size(); ++i) {
            const Opening& x = a.openings[i];
            const Opening& y = b.openings[i];
            structuralUnless(x.name == y.name && x.a == y.a && x.b == y.b);
            structuralUnless(x.area == y.area && x.dischargeCoeff == y.dischargeCoeff);
            structuralUnless(x.open == y.open && x.kind == y.kind);
            structuralUnless(x.pos.x == y.pos.x && x.pos.y == y.pos.y && x.pos.z == y.pos.z);
        }

    structuralUnless(a.pumps.size() == b.pumps.size());
    if (a.pumps.size() == b.pumps.size())
        for (std::size_t i = 0; i < a.pumps.size(); ++i) {
            const Pump& x = a.pumps[i];
            const Pump& y = b.pumps[i];
            structuralUnless(x.name == y.name && x.compartment == y.compartment);
            structuralUnless(x.capacity == y.capacity && x.maxHead == y.maxHead && x.on == y.on);
        }

    structuralUnless(a.deckEdgeZ == b.deckEdgeZ);
    structuralUnless(a.seaDensity == b.seaDensity);
    structuralUnless(a.lightshipCog.x == b.lightshipCog.x &&
                     a.lightshipCog.y == b.lightshipCog.y && a.lightshipCog.z == b.lightshipCog.z);
    structuralUnless(a.gyradii.x == b.gyradii.x && a.gyradii.y == b.gyradii.y &&
                     a.gyradii.z == b.gyradii.z);
    structuralUnless(a.zetaHeave == b.zetaHeave && a.zetaRoll == b.zetaRoll &&
                     a.zetaPitch == b.zetaPitch);
    structuralUnless(a.addedMassSurge == b.addedMassSurge && a.addedMassSway == b.addedMassSway &&
                     a.addedMassHeave == b.addedMassHeave);
    structuralUnless(a.addedInertiaRoll == b.addedInertiaRoll &&
                     a.addedInertiaPitch == b.addedInertiaPitch &&
                     a.addedInertiaYaw == b.addedInertiaYaw);
    d.massRel = std::abs(a.lightshipMass - b.lightshipMass) / a.lightshipMass;

    const Diagnostics da = a.diagnostics(0.0);
    const Diagnostics db = b.diagnostics(0.0);
    d.massRel = std::max(d.massRel, std::abs(da.displacementMass - db.displacementMass) /
                                        da.displacementMass);
    d.volumeRel = std::max(d.volumeRel,
                           std::abs(da.buoyantVolume - db.buoyantVolume) / da.buoyantVolume);
    d.volumeRel = std::max(d.volumeRel,
                           std::abs(da.waterplaneArea - db.waterplaneArea) / da.waterplaneArea);
    d.length = std::max({d.length, std::abs(da.draftMidship - db.draftMidship),
                         std::abs(da.gmTransverse - db.gmTransverse),
                         std::abs(da.centreOfBuoyancy.z - db.centreOfBuoyancy.z),
                         std::abs(da.freeboardMin - db.freeboardMin)});
    // The whole righting arm curve, not just the value at zero heel: GZ is what
    // decides whether the ship survives, and it is the integral of the shape of
    // the hull well above the waterline -- the part a lazy format would omit.
    for (int degrees = 0; degrees <= 60; degrees += 5) {
        const double heel = degrees * sim::kDegToRad;
        d.length = std::max(d.length,
                            std::abs(a.rightingArmAtHeel(heel, 0.0) -
                                     b.rightingArmAtHeel(heel, 0.0)));
    }
    return d;
}

// --- The damage-control scenario, as tests/../main.cpp schedules it -------------

struct Sample {
    double t = 0, draft = 0, heel = 0, trim = 0, gm = 0, floodwater = 0;
};

void setOpening(Ship& s, std::string_view name, bool open) {
    for (Opening& o : s.openings)
        if (o.name == name) o.open = open;
}
void setPump(Ship& s, std::string_view name, bool on) {
    for (Pump& p : s.pumps)
        if (p.name == name) p.on = on;
}

// `which` matches game/prototype/main.cpp's scenarios exactly, including the
// times: the point is to reproduce the runs `verify.sh` treats as validated
// behaviour, not to invent a new one.
void applyEvents(Ship& s, const std::string& which, double t, double dt) {
    const auto fires = [&](double when) { return t <= when && when < t + dt; };
    if (which == "none") return;
    if (fires(45.0)) setOpening(s, "wt_door_er", false);
    if (fires(60.0)) {
        setPump(s, "bilge_er_s", true);
        setPump(s, "bilge_er_p", true);
        setPump(s, "bilge_ah_s", true);
    }
    if (which != "full") return;
    if (fires(90.0)) {
        setOpening(s, "cf_valve_aft_p", true);
        setOpening(s, "cf_valve_fwd_p", true);
    }
    if (fires(120.0)) {
        setOpening(s, "downflood_ramp_s", false);
        setOpening(s, "downflood_port_s", false);
        setOpening(s, "downflood_port_p", false);
    }
}

// Steps `ship` in place, so the caller can inspect what the schedule did to it.
// That is not a convenience: it is how this suite proves the events fired at all,
// and it doubles as a check that the names in the file are the names the damage
// control board addresses.
std::vector<Sample> runScenario(Ship& ship, const std::string& which, double duration, double dt) {
    std::vector<Sample> track;
    double t = 0, nextReport = 0;
    while (t < duration) {
        applyEvents(ship, which, t, dt);
        if (t >= nextReport) {
            const Diagnostics d = ship.diagnostics(0.0);
            track.push_back({t, d.draftMidship, d.heelDeg, d.trimDeg, d.gmTransverse,
                             d.floodwaterMass});
            nextReport += 15.0;
        }
        ship.step(dt, 0.0);
        t += dt;
    }
    return track;
}

// Worst disagreement between two tracks, as one number. Draft and GM in metres,
// heel and trim in degrees, floodwater *relative* -- a few grams of water out of
// fourteen hundred tonnes is not a divergence, and dividing tonnes by an
// arbitrary constant to make it comparable would report it as the largest one.
double trackDivergence(const std::vector<Sample>& a, const std::vector<Sample>& b) {
    if (a.size() != b.size() || a.empty()) return 1e9;
    double worst = 0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        const double flood = std::max({a[i].floodwater, b[i].floodwater, 1.0});
        worst = std::max({worst, std::abs(a[i].draft - b[i].draft), std::abs(a[i].heel - b[i].heel),
                          std::abs(a[i].trim - b[i].trim), std::abs(a[i].gm - b[i].gm),
                          std::abs(a[i].floodwater - b[i].floodwater) / flood});
    }
    return worst;
}

bool openingIs(const Ship& s, std::string_view name, bool open) {
    for (const Opening& o : s.openings)
        if (o.name == name) return o.open == open;
    return false;
}
bool pumpIs(const Ship& s, std::string_view name, bool on) {
    for (const Pump& p : s.pumps)
        if (p.name == name) return p.on == on;
    return false;
}

// --- Tests ---------------------------------------------------------------------

void testFerryFileMatchesTheCompiledFerry() {
    ShipDefinition loaded;
    std::string error;
    const std::string path = shipDir() + "ferry.ship";
    const bool ok = sim::loadShipFile(path, loaded, error);
    expectTrue("ships/ferry.ship loads: " + error, ok);
    if (!ok) return;
    expectTrue("the file names the ship", loaded.name == "ferry_120m");

    Ship compiled = game::buildFerry();
    compiled.initialise(0.0);
    loaded.ship.initialise(0.0);

    // Guard against a vacuous comparison. Two empty ships agree perfectly.
    expectEqual("the ferry has its sixteen compartments",
                static_cast<long long>(loaded.ship.compartments.size()), 16);
    expectEqual("and its twenty-five openings",
                static_cast<long long>(loaded.ship.openings.size()), 25);
    expectEqual("and its four pumps", static_cast<long long>(loaded.ship.pumps.size()), 4);
    const Diagnostics intact = loaded.ship.diagnostics(0.0);
    expectNear("intact displacement is the ferry's, near 8984 t",
               intact.displacementMass / 1000.0, 8984.0, 5.0);
    expectNear("intact GM is in the ro-pax band", intact.gmTransverse, 2.0, 0.3);

    // The loader's own structural check must have nothing to say, and neither
    // must the simulation's -- an overlapping subdivision or a dangling index
    // would be `validate()`'s to find, and it finds neither.
    expectEqual("the loaded ferry passes Ship::validate()",
                static_cast<long long>(loaded.ship.validate().size()), 0);

    const Divergence d = compare(compiled, loaded.ship);
    expectTrue("every non-numeric field survives the round trip", !d.structural);
    expectTrue(label("hull vertices agree to 1 um", d.vertex), d.vertex < kVertexTol);
    expectTrue(label("compartment, buoyant and waterplane volumes agree", d.volumeRel),
               d.volumeRel < kVolumeRelTol);
    expectTrue(label("displacement agrees", d.massRel), d.massRel < kMassRelTol);
    expectTrue(label("draft, GM, KB, freeboard and the whole GZ curve agree to 1 um", d.length),
               d.length < kLengthTol);
}

// The tolerances above are only worth stating if they can fail. One station's
// half-breadth moved by a centimetre -- a hundredth of the smallest edit anyone
// would make on purpose -- has to break every one of them.
void testTheComparisonCanFail() {
    const std::string text = readFile(shipDir() + "ferry.ship");
    expectTrue("ferry.ship is readable as text", text.size() > 1000);
    // Station +25 at the 3.2 m waterline: unique in the file, and inside the
    // wetted envelope at the design draft.
    const std::string perturbed = mutate(text, "8.915221", "8.905221");

    ShipDefinition bad;
    std::string error;
    const bool ok = sim::parseShipFile(perturbed, "<perturbed>", bad, error);
    expectTrue("the perturbed file is still well-formed: " + error, ok);
    if (!ok) return;

    Ship compiled = game::buildFerry();
    compiled.initialise(0.0);
    bad.ship.initialise(0.0);
    const Divergence d = compare(compiled, bad.ship);
    expectTrue(label("a 1 cm offsets error moves the hull past the vertex tolerance", d.vertex),
               d.vertex > kVertexTol);
    expectTrue(label("and past the volume tolerance", d.volumeRel), d.volumeRel > kVolumeRelTol);
    expectTrue(label("and past the displacement tolerance", d.massRel), d.massRel > kMassRelTol);

    // And the other half of the comparison, which is the half that would go
    // unnoticed. A machinery space at 0.95 permeability instead of 0.85 changes
    // no vertex, no volume, no hydrostatic number and no intact GZ -- the ship
    // floats identically and floods twelve per cent faster. Only the field-by-
    // field check sees it, which is why there is one.
    ShipDefinition leaky;
    if (!sim::parseShipFile(mutate(text, "permeability 0.85", "permeability 0.95"), "<leaky>",
                            leaky, error)) {
        expectTrue("the permeability mutation is well-formed: " + error, false);
        return;
    }
    leaky.ship.initialise(0.0);
    const Divergence p = compare(compiled, leaky.ship);
    expectTrue("a wrong permeability is invisible to every geometric tolerance",
               p.vertex < kVertexTol && p.volumeRel < kVolumeRelTol && p.massRel < kMassRelTol &&
                   p.length < kLengthTol);
    expectTrue("but the field comparison catches it", p.structural);
}

// Equal hydrostatics is a static claim. This is the dynamic one: the same
// breach, the same doors closed at the same second, the same pumps started, and
// the two ships must trace the same trajectory.
void testFloodingScenariosMatch() {
    ShipDefinition loaded;
    std::string error;
    if (!sim::loadShipFile(shipDir() + "ferry.ship", loaded, error)) {
        expectTrue("ships/ferry.ship loads for the scenario run: " + error, false);
        return;
    }
    Ship compiled = game::buildFerry();
    compiled.initialise(0.0);
    loaded.ship.initialise(0.0);

    // 'full' is the scenario that touches everything: watertight door at t+45,
    // bilge pumps at t+60, counterflooding valves at t+90, vehicle deck secured
    // at t+120. 150 s therefore exercises the whole flow network, every opening
    // flag and every pump. dt is 0.04 s -- the top of the band main.cpp records
    // as converged to four significant figures -- because this runs on every
    // save. The 900 s runs of all three scenarios, which are what decide the
    // published outcomes, are `verify.sh full`'s job.
    constexpr double kDuration = 150.0;
    constexpr double kDt = 0.04;
    // Measured: the worst disagreement over these 150 s is 3.8e-7, in heel
    // degrees -- the 1 um offsets rounding entering heel as roughly its ratio to
    // the beam. It grows *linearly*, not exponentially: 1e-6 by 285 s on the
    // 'doors' run, which is why the 900 s scenarios in verify.sh still land on
    // the same verdict. The tolerance is that with a factor of twenty-five in
    // hand, and the label prints what was actually measured.
    constexpr double kTrackTol = 1e-5;

    const std::vector<Sample> a = runScenario(compiled, "full", kDuration, kDt);
    const std::vector<Sample> b = runScenario(loaded.ship, "full", kDuration, kDt);
    const double divergence = trackDivergence(a, b);
    expectTrue(label("the loaded ferry floods exactly as the compiled one does", divergence),
               divergence < kTrackTol);

    // Guards. Two ships that both sat still would track together perfectly, and
    // a schedule that silently addressed no one would leave both as 'none'.
    expectTrue(label("the scenario actually floods the ship", a.back().floodwater / 1000.0),
               a.back().floodwater > 500e3 && a.back().draft > a.front().draft + 0.05 &&
                   std::abs(a.back().heel) > 1.0);
    for (const Ship& s : {compiled, loaded.ship}) {
        expectTrue("the damage control schedule found the watertight door by name",
                   openingIs(s, "wt_door_er", false));
        expectTrue("... the bilge pumps", pumpIs(s, "bilge_er_s", true));
        expectTrue("... the counterflooding valves", openingIs(s, "cf_valve_aft_p", true));
        expectTrue("... and the vehicle deck shell doors",
                   openingIs(s, "downflood_ramp_s", false));
    }
}

// The negative control for the scenario comparison. A breach 4% larger is a
// difference no static field check would catch if `area` were mis-parsed by a
// digit, and the trajectory must reject it by orders of magnitude.
void testAMisCopiedFileDivergesInTheScenario() {
    const std::string text = readFile(shipDir() + "ferry.ship");
    ShipDefinition bad;
    std::string error;
    if (!sim::parseShipFile(mutate(text, "area 2.4", "area 2.5"), "<perturbed>", bad, error)) {
        expectTrue("the perturbed file is well-formed: " + error, false);
        return;
    }
    Ship compiled = game::buildFerry();
    compiled.initialise(0.0);
    bad.ship.initialise(0.0);

    const std::vector<Sample> a = runScenario(compiled, "full", 30.0, 0.04);
    const std::vector<Sample> b = runScenario(bad.ship, "full", 30.0, 0.04);
    expectTrue(label("a 4% larger breach diverges within 30 s", trackDivergence(a, b)),
               trackDivergence(a, b) > 1e-3);
}

// The ferry file's six decimals are an authoring choice, not the format's limit.
// Written at full precision the offsets come back bit for bit, which is what
// makes an exporter -- the ship editor docs 05 §3 commits to -- possible at all.
void testFormatIsLosslessAtFullPrecision() {
    const std::vector<double> waterlines{0.0, 1.0, 2.5, 4.0};
    std::vector<sim::Station> stations;
    std::ostringstream text;
    text.precision(17);
    text << "ship_format 1\nship exact\n deck_edge_z 3.0\n lightship_mass 1e6\n"
            " lightship_cog 0 0 2\n gyradii 3 8 8\nhull\n waterlines";
    for (double z : waterlines) text << " " << z;
    text << "\n";
    for (int i = 0; i < 5; ++i) {
        sim::Station s;
        s.x = -20.0 + 10.0 * i;
        text << " station " << s.x;
        for (double z : waterlines) {
            // Irrational on purpose: a value with no short decimal form is the
            // only kind that can expose a lossy round trip.
            const double halfBeam = 6.0 * std::pow(1.0 + z, 1.0 / 3.0) * std::sin(1.0 + 0.1 * i);
            s.halfBeam.push_back(halfBeam);
            text << " " << halfBeam;
        }
        text << "\n";
        stations.push_back(s);
    }
    // Vacuity guard: if these all had short decimal forms the test would prove
    // nothing about precision.
    std::ostringstream shortForm;
    shortForm << stations[0].halfBeam[1];
    expectTrue("the probe offsets have no short decimal form",
               std::strtod(shortForm.str().c_str(), nullptr) != stations[0].halfBeam[1]);

    ShipDefinition def;
    std::string error;
    const bool ok = sim::parseShipFile(text.str(), "<exact>", def, error);
    expectTrue("a full-precision file parses: " + error, ok);
    if (!ok) return;

    const sim::TriMesh reference = sim::makeHullFromStations(stations, waterlines);
    bool bitExact = reference.verts.size() == def.ship.hull.verts.size();
    for (std::size_t i = 0; bitExact && i < reference.verts.size(); ++i)
        bitExact = reference.verts[i].x == def.ship.hull.verts[i].x &&
                   reference.verts[i].y == def.ship.hull.verts[i].y &&
                   reference.verts[i].z == def.ship.hull.verts[i].z;
    expectTrue("seventeen significant figures round-trip bit for bit", bitExact);
}

// --- Failing closed ------------------------------------------------------------

void testEveryTruncationFailsClosed() {
    const std::string valid(kBarge);
    ShipDefinition whole;
    std::string error;
    expectTrue("the test barge is a valid ship: " + error,
               sim::parseShipFile(valid, "<barge>", whole, error));

    // Whether a prefix is legal is predictable without the loader: it is legal
    // exactly when it carries a stamped format, a complete ship block, a hull
    // with a waterline list and two stations, and no half-finished block after
    // them. Derived here rather than observed, so the expectation is independent
    // of the implementation it is checking.
    std::vector<std::size_t> lineEnds;
    for (std::size_t i = 0; i < valid.size(); ++i)
        if (valid[i] == '\n') lineEnds.push_back(i + 1);
    expectTrue("the barge has lines to cut at", lineEnds.size() > 20);

    long long accepted = 0, rejected = 0;
    bool asPredicted = true, nothingHalfBuilt = true;
    for (std::size_t cut : lineEnds) {
        const std::string prefix = valid.substr(0, cut);

        bool stamped = false, inShip = false, inHull = false, hasWaterlines = false;
        int shipKeys = 0, stations = 0, blockKeys = 0, wantKeys = 0;
        int compartments = 0, openings = 0, pumps = 0;
        std::size_t position = 0;
        while (position < prefix.size()) {
            const std::size_t end = prefix.find('\n', position);
            std::string line =
                prefix.substr(position, end == std::string::npos ? std::string::npos
                                                                 : end - position);
            position = end == std::string::npos ? prefix.size() : end + 1;
            const std::size_t hash = line.find('#');
            if (hash != std::string::npos) line = line.substr(0, hash);
            std::istringstream words(line);
            std::string first;
            if (!(words >> first)) continue;
            if (first == "ship_format") stamped = true;
            else if (first == "ship") { inShip = true; wantKeys = 0; }
            else if (first == "hull") { inHull = true; wantKeys = 0; }
            else if (first == "compartment") { ++compartments; blockKeys = 0; wantKeys = 3; }
            else if (first == "opening") { ++openings; blockKeys = 0; wantKeys = 6; }
            else if (first == "pump") { ++pumps; blockKeys = 0; wantKeys = 4; }
            else if (first == "waterlines") hasWaterlines = true;
            else if (first == "station") ++stations;
            else if (wantKeys > 0) ++blockKeys;
            else if (inShip && !inHull) ++shipKeys;
        }
        const bool shouldLoad = stamped && shipKeys == 4 && inHull && hasWaterlines &&
                                stations >= 2 && (wantKeys == 0 || blockKeys == wantKeys);

        ShipDefinition out = sentinel();
        std::string why;
        const bool loaded = sim::parseShipFile(prefix, "<truncated>", out, why);
        asPredicted = asPredicted && loaded == shouldLoad;
        if (loaded) {
            ++accepted;
            // Whatever survived must be the complete blocks it came from -- and
            // exactly those, with their values, not a name and a default.
            nothingHalfBuilt =
                nothingHalfBuilt && out.name == "barge" && out.ship.deckEdgeZ == 3.0 &&
                out.ship.lightshipMass == 1.2e6 && !out.ship.hull.tris.empty() &&
                out.ship.compartments.size() == static_cast<std::size_t>(compartments) &&
                out.ship.openings.size() == static_cast<std::size_t>(openings) &&
                out.ship.pumps.size() == static_cast<std::size_t>(pumps);
            if (!out.ship.compartments.empty())
                nothingHalfBuilt = nothingHalfBuilt &&
                                   out.ship.compartments[0].name == "hold_fwd" &&
                                   out.ship.compartments[0].permeability == 0.95 &&
                                   !out.ship.compartments[0].ventedToAtmosphere;
        } else {
            ++rejected;
            nothingHalfBuilt = nothingHalfBuilt && isSentinel(out);
            expectTrue("a rejected prefix says why", !why.empty());
        }
    }
    expectTrue("every line truncation loads exactly when it ends after a complete ship",
               asPredicted);
    expectTrue("no truncation produced a half-built ship", nothingHalfBuilt);
    // A loader that accepted everything, or nothing, satisfies one of those on
    // its own.
    expectTrue(label("some truncations were accepted", static_cast<double>(accepted)),
               accepted >= 5);
    expectTrue(label("and some were rejected", static_cast<double>(rejected)), rejected >= 10);

    // And every *byte* truncation, which can cut a number in half. The invariant
    // is weaker -- "0.9" is a legal truncation of "0.95" -- but the load-bearing
    // half is not: a failure must leave nothing behind, and a success must never
    // invent a compartment that was not declared.
    long long byteAccepted = 0;
    bool bytesSane = true;
    for (std::size_t cut = 0; cut <= valid.size(); ++cut) {
        ShipDefinition out = sentinel();
        std::string why;
        if (sim::parseShipFile(valid.substr(0, cut), "<bytes>", out, why)) {
            ++byteAccepted;
            bytesSane = bytesSane && out.name == "barge" && out.ship.compartments.size() <= 2 &&
                        out.ship.openings.size() <= 2 && out.ship.pumps.size() <= 1 &&
                        sim::isClosedManifold(out.ship.hull) &&
                        sim::integrate(out.ship.hull).volume > 0;
        } else {
            bytesSane = bytesSane && isSentinel(out);
        }
    }
    expectTrue("no byte truncation invents a space or survives a failure", bytesSane);
    expectTrue(label("byte truncations accepted", static_cast<double>(byteAccepted)),
               byteAccepted > 0 && byteAccepted < static_cast<long long>(valid.size()));
}

void testMalformedFilesAreRefusedByName() {
    const std::string barge(kBarge);
    struct Case {
        const char* what;
        std::string text;
        const char* expect;  // substring the diagnostic must carry
    };
    const std::vector<Case> cases = {
        {"an empty file", "", "ship_format"},
        {"a file of nothing but comments", "# a ship\n# honest\n", "ship_format"},
        {"a file with no format stamp", "ship barge\n", "must begin with"},
        {"a format from the future", mutate(barge, "ship_format 1", "ship_format 2"),
         "newer than this build"},
        {"a fractional format", mutate(barge, "ship_format 1", "ship_format 1.5"),
         "whole version number"},
        {"an unknown key", mutate(barge, "permeability 0.95", "permeabilty 0.95"), "unknown key"},
        {"a key outside any block", mutate(barge, "ship barge", "deck_edge_z 3.0"),
         "outside any block"},
        {"an opening onto a compartment that does not exist",
         mutate(barge, "between sea hold_fwd", "between sea hold_middle"), "does not exist"},
        {"a pump draining a compartment that does not exist",
         mutate(barge, "drains hold_fwd", "drains hold_middle"), "does not exist"},
        {"a pump draining the sea", mutate(barge, "drains hold_fwd", "drains sea"),
         "cannot draw from the sea"},
        {"a negative opening area", mutate(barge, "area 0.5", "area -0.5"), "must be within (0,"},
        {"a zero opening area", mutate(barge, "area 0.5", "area 0"), "must be within (0,"},
        {"a discharge coefficient above one", mutate(barge, "discharge 0.6", "discharge 1.4"),
         "must be within (0, 1]"},
        {"a permeability above one", mutate(barge, "permeability 0.95", "permeability 1.5"),
         "must be within (0, 1]"},
        {"stations out of order",
         mutate(barge, "station   0.0  5.0 5.0 5.0", "station -12.0  5.0 5.0 5.0"),
         "must ascend in x"},
        {"a repeated station", mutate(barge, "station   0.0  5.0", "station -10.0  5.0"),
         "must ascend in x"},
        {"waterlines out of order", mutate(barge, "waterlines 0.0 2.0 4.0", "waterlines 0.0 4.0 2.0"),
         "waterlines must ascend"},
        {"too few half-breadths", mutate(barge, "station   0.0  5.0 5.0 5.0", "station 0.0 5.0 5.0"),
         "one half-breadth per waterline"},
        {"too many half-breadths",
         mutate(barge, "station   0.0  5.0 5.0 5.0", "station 0.0 5.0 5.0 5.0 5.0"),
         "one half-breadth per waterline"},
        {"a negative half-breadth",
         mutate(barge, "station   0.0  5.0 5.0 5.0", "station 0.0 -5.0 5.0 5.0"),
         "cannot be negative"},
        {"stations before the waterline list",
         mutate(barge, "waterlines 0.0 2.0 4.0", "# waterlines come later"),
         "needs the waterline list"},
        {"a number where a name should be", mutate(barge, "compartment hold_fwd", "compartment 3.5"),
         "is a number, not a name"},
        {"a compartment called sea", mutate(barge, "compartment hold_fwd", "compartment sea"),
         "reserved name"},
        {"a duplicated compartment name", mutate(barge, "compartment hold_aft", "compartment hold_fwd"),
         "defined twice"},
        {"a duplicated opening name", mutate(barge, "opening door", "opening breach"),
         "defined twice"},
        {"a key set twice in one block",
         mutate(barge, "permeability 0.95", "permeability 0.95\n    permeability 0.5"),
         "set twice"},
        {"two ships in one file", mutate(barge, "hull\n", "ship other\nhull\n"),
         "exactly one ship"},
        {"a missing required key", mutate(barge, "    permeability 0.95\n", ""),
         "has no permeability"},
        {"a compartment that does not say whether it is vented",
         mutate(barge, "    vented false\n", ""), "vented"},
        {"an opening with no position", mutate(barge, "    at 4 -5 1\n", ""), "has no at"},
        {"a ship with no lightship weight", mutate(barge, "    lightship_mass 1.2e6\n", ""),
         "has neither"},
        {"a ship stating mass and draft both",
         mutate(barge, "lightship_mass 1.2e6", "lightship_mass 1.2e6\n    lightship_draft 2.0"),
         "would disagree"},
        {"a ship with no hull",
         "ship_format 1\nship hull_less\n deck_edge_z 3\n lightship_mass 1e6\n"
         " lightship_cog 0 0 2\n gyradii 3 8 8\n",
         "no 'hull' block"},
        {"a hull with one station",
         mutate(barge, "    station   0.0  5.0 5.0 5.0\n    station  10.0  4.0 5.0 5.0\n", ""),
         "at least two stations"},
        // Dropping just the `hull` header leaves the offsets stranded in the
        // ship block, and that is what the diagnostic must say -- naming the
        // first line that cannot be where it is, rather than the missing block.
        {"offsets outside a hull block", mutate(barge, "hull\n", "# where did the hull go\n"),
         "unknown key 'waterlines' in a ship block"},
        {"a box that misses the hull", mutate(barge, "box 0 -6 0   12 6 3", "box 40 -6 0  60 6 3"),
         "clips to nothing"},
        {"an inside-out box", mutate(barge, "box 0 -6 0   12 6 3", "box 12 -6 0  0 6 3"),
         "hi <= lo"},
        {"an opening onto itself", mutate(barge, "between hold_fwd hold_aft", "between hold_fwd hold_fwd"),
         "connects a space to itself"},
        {"a unit suffix on a number", mutate(barge, "area 0.5", "area 0.5m"),
         "is not a finite number"},
        {"a not-a-number", mutate(barge, "area 0.5", "area nan"), "is not a finite number"},
        {"an infinity", mutate(barge, "lightship_cog 0.0 0.0 2.0", "lightship_cog 0.0 0.0 inf"),
         "is not a finite number"},
        {"a boolean that is not one", mutate(barge, "open true", "open yes"),
         "is not true or false"},
        {"an opening kind nobody defined", mutate(barge, "kind breach", "kind window"),
         "is not one of breach"},
        {"a missing argument", mutate(barge, "area 0.5", "area"), "takes one number"},
        {"a surplus argument", mutate(barge, "deck_edge_z 3.0", "deck_edge_z 3.0 4.0"),
         "takes one number"},
        // A bare `station` has fewer tokens than the two the arity message
        // subtracts. The count has to clamp, not wrap around an unsigned zero.
        {"a station with nothing on it", mutate(barge, "station -10.0  4.0 5.0 5.0", "station"),
         "expected 3, got 0"},
        {"a single waterline", mutate(barge, "waterlines 0.0 2.0 4.0", "waterlines 0.0"),
         "at least two heights"},
        {"a pump that does not say what it drains", mutate(barge, "    drains hold_fwd\n", ""),
         "does not say what it drains"},
        {"an unknown key in a hull block", mutate(barge, "    waterlines", "    waterline"),
         "unknown key 'waterline' in a hull block"},
        {"an unknown key in a pump block", mutate(barge, "    max_head 20", "    maxhead 20"),
         "unknown key 'maxhead' in a pump block"},
        {"an unknown key in an opening block", mutate(barge, "    discharge 0.6", "    cd 0.6"),
         "unknown key 'cd' in an opening block"},
        // The weather deck is what `freeboardMin` -- and therefore the
        // prototype's verdict -- is measured to. A deck edge outside the
        // offsets table is a transposed digit, not a tall ship.
        {"a deck edge above the hull", mutate(barge, "deck_edge_z 3.0", "deck_edge_z 30.0"),
         "outside the hull"},
        {"a deck edge below the keel", mutate(barge, "deck_edge_z 3.0", "deck_edge_z -1.0"),
         "outside the hull"},
    };

    for (const Case& c : cases) {
        ShipDefinition out = sentinel();
        std::string why;
        const bool loaded = sim::parseShipFile(c.text, "barge.ship", out, why);
        expectTrue(std::string("refused: ") + c.what, !loaded);
        if (loaded) continue;
        expectTrue(std::string("... naming the problem, not just failing: ") + c.what + " -> " + why,
                   why.find(c.expect) != std::string::npos);
        // Every diagnostic carries the origin, and every one that can name a line
        // does. Only the whole-file complaints ("no hull block") legitimately
        // cannot.
        expectTrue(std::string("... and the origin: ") + c.what,
                   why.rfind("barge.ship", 0) == 0);
        expectTrue(std::string("... leaving nothing behind: ") + c.what, isSentinel(out));
    }
}

void testTheFilePathItselfFailsClosed() {
    ShipDefinition out = sentinel();
    std::string error;
    expectTrue("a path that does not exist is refused",
               !sim::loadShipFile(testing::scratchDir() + "no-such-ship.ship", out, error));
    expectTrue("... naming the path", error.find("no-such-ship.ship") != std::string::npos);
    expectTrue("... and leaving the definition alone", isSentinel(out));

    // A directory, which opens but does not read as a file on most systems --
    // and where it does read, it reads as garbage, which must also be refused.
    out = sentinel();
    expectTrue("a directory is not a ship", !sim::loadShipFile(testing::scratchDir(), out, error));
    expectTrue("... leaving the definition alone", isSentinel(out));

    const std::string path = testing::scratchDir() + "barge.ship";
    {
        std::ofstream file(path, std::ios::binary);
        file << kBarge;
    }
    out = sentinel();
    expectTrue("a ship written to disk loads back: " + error,
               sim::loadShipFile(path, out, error));
    expectTrue("... as the ship it was", out.name == "barge" && out.ship.compartments.size() == 2);
}

// Two ways to state the same lightship weight, and the format has to mean the
// same thing by both. `lightship_draft` is the one the ferry uses, because it
// keeps the hull form and the loading consistent when the offsets are edited.
void testLightshipDraftAndMassAgree() {
    ShipDefinition byMass, byDraft;
    std::string error;
    const std::string base(kBarge);
    expectTrue("the barge loads by mass: " + error,
               sim::parseShipFile(base, "<mass>", byMass, error));
    expectTrue("and by draft: " + error,
               sim::parseShipFile(mutate(base, "lightship_mass 1.2e6", "lightship_draft 2.0"),
                                  "<draft>", byDraft, error));

    // Independently derived: the barge's own hull, integrated below z = 2.
    const double volume =
        sim::integrateBelowPlane(byDraft.ship.hull, Vec3{0, 0, 1}, 2.0).volume;
    expectNear("lightship_draft is rho times the volume it displaces", byDraft.ship.lightshipMass,
               volume * sim::kRhoSeawater, 1e-6 * volume * sim::kRhoSeawater);
    expectTrue("which is not the mass the other file states",
               std::abs(byDraft.ship.lightshipMass - byMass.ship.lightshipMass) > 1.0);

    // And a draft that floats nothing is refused rather than producing a
    // weightless ship. The hull here starts at z = 6, so a 2 m draft is under
    // the keel -- the case a ship whose baseline is not at the origin can hit.
    ShipDefinition out = sentinel();
    const std::string dry = R"(ship_format 1
ship stilts
    deck_edge_z 9
    lightship_draft 2.0
    lightship_cog 0 0 8
    gyradii 3 8 8
hull
    waterlines 6.0 8.0 10.0
    station -10.0  4.0 5.0 5.0
    station  10.0  4.0 5.0 5.0
)";
    expectTrue("a lightship draft below the keel is refused",
               !sim::parseShipFile(dry, "<dry>", out, error));
    expectTrue("... naming what is wrong: " + error,
               error.find("no part of the hull in the water") != std::string::npos);
    expectTrue("... leaving nothing behind", isSentinel(out));

    // And so is one above the weather deck, which would float her at the light
    // condition with the whole envelope under.
    out = sentinel();
    expectTrue("a lightship draft above the hull is refused",
               !sim::parseShipFile(mutate(base, "lightship_mass 1.2e6", "lightship_draft 6.0"),
                                   "<drowned>", out, error));
    expectTrue("... leaving nothing behind", isSentinel(out));
}

// The optional keys are optional because the engine's defaults are honest
// placeholders, not because they do not matter. Stating them has to take.
void testOptionalKeysAreCarried() {
    ShipDefinition def;
    std::string error;
    std::string text(kBarge);
    text = mutate(text, "    gyradii 3.0 8.0 8.0\n",
                  "    gyradii 3.0 8.0 8.0\n"
                  "    sea_density 998.2\n"
                  "    damping 0.4 0.0 0.25\n"
                  "    added_mass 0.1 0.8 1.2\n"
                  "    added_inertia 0.3 1.1 0.7\n");
    expectTrue("a ship stating every optional key loads: " + error,
               sim::parseShipFile(text, "<optional>", def, error));
    expectNear("sea_density is carried", def.ship.seaDensity, sim::kRhoFresh, 0);
    expectNear("heave damping is carried", def.ship.zetaHeave, 0.4, 0);
    // Zero roll damping is legal and load-bearing: when a radiation model
    // supplies the damping, the modal stand-in must go to zero or it is counted
    // twice -- a defect this repo has already shipped once.
    expectNear("zero roll damping is legal", def.ship.zetaRoll, 0.0, 0);
    expectNear("sway added mass is carried", def.ship.addedMassSway, 0.8, 0);
    expectNear("yaw added inertia is carried", def.ship.addedInertiaYaw, 0.7, 0);

    // Omitted, they leave the engine's values exactly as they were.
    ShipDefinition plain;
    expectTrue("and a ship stating none of them loads: " + error,
               sim::parseShipFile(kBarge, "<plain>", plain, error));
    const Ship untouched;
    expectNear("sea density defaults to seawater", plain.ship.seaDensity, untouched.seaDensity, 0);
    expectNear("roll damping keeps the engine's default", plain.ship.zetaRoll,
               untouched.zetaRoll, 0);
    expectNear("heave added mass keeps the engine's default", plain.ship.addedMassHeave,
               untouched.addedMassHeave, 0);
}

// A mod may reach a compartment declared later in the file. Nothing in the
// format is order-dependent except the waterline list, which the station lines
// genuinely need first.
void testReferencesResolveInAnyOrder() {
    const std::string forward = R"(ship_format 1
ship reordered
    deck_edge_z 3.0
    lightship_mass 1.2e6
    lightship_cog 0 0 2
    gyradii 3 8 8
pump bilge
    drains hold
    capacity 0.05
    max_head 20
    on false
opening breach
    between sea hold
    at 4 -5 1
    area 0.5
    discharge 0.6
    kind breach
    open true
compartment hold
    box -12 -6 0   12 6 3
    permeability 0.95
    vented false
hull
    waterlines 0.0 2.0 4.0
    station -10.0  4.0 5.0 5.0
    station   0.0  5.0 5.0 5.0
    station  10.0  4.0 5.0 5.0
)";
    ShipDefinition def;
    std::string error;
    expectTrue("blocks resolve regardless of the order they appear in: " + error,
               sim::parseShipFile(forward, "<reordered>", def, error));
    expectEqual("the forward reference resolved", def.ship.openings.at(0).b, 0);
    expectEqual("and so did the pump's", def.ship.pumps.at(0).compartment, 0);
}

}  // namespace

void runShipFileTests() {
    std::printf("\n--- ship definition format ---\n");
    testFerryFileMatchesTheCompiledFerry();
    testTheComparisonCanFail();
    testFloodingScenariosMatch();
    testAMisCopiedFileDivergesInTheScenario();
    testFormatIsLosslessAtFullPrecision();
    testEveryTruncationFailsClosed();
    testMalformedFilesAreRefusedByName();
    testTheFilePathItselfFailsClosed();
    testLightshipDraftAndMassAgree();
    testOptionalKeysAreCarried();
    testReferencesResolveInAnyOrder();
}
