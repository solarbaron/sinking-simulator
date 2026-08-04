// SPDX-License-Identifier: MIT
#include "ferry.hpp"

#include <cmath>

namespace game {

using namespace sim;

namespace {

constexpr double kHalfBeam    = 10.0;   // moulded half breadth, m
constexpr double kDesignDraft = 5.5;    // m
constexpr double kDepth       = 15.0;   // m to the top of the accommodation block
// Ro-pax ferries carry a great deal of steel and passenger accommodation high up,
// and run intact GM in the 1.5-3 m band -- stiff enough to pass the criteria,
// nowhere near stiff enough to survive water on the vehicle deck.
constexpr double kKG          = 7.80;   // m, vertical centre of gravity, light condition

// Longitudinal fullness: parallel midbody between +/-24 m, a fine bow and a
// fuller stern. u is normalised half-length.
double longitudinalFullness(double x) {
    const double u = x / 60.0;
    if (u > 0.4) {
        const double t = (u - 0.4) / 0.6;
        return 1.0 - 0.88 * t * t;
    }
    if (u < -0.4) {
        const double t = (-u - 0.4) / 0.6;
        return 1.0 - 0.55 * t * t;
    }
    return 1.0;
}

// Vertical fullness: rise of floor from a narrow keel out to full beam at the
// turn of the bilge, 4 m up. Together with the above this gives Cb = 0.66 and
// Cwp = 0.86 at the design draft -- a fine-ended but wide-waterplane form, which
// is what makes a ferry's damaged stability fragile.
double verticalFullness(double z) {
    if (z >= 4.0) return 1.0;
    return 0.30 + 0.70 * std::pow(z / 4.0, 0.75);
}

TriMesh buildHull() {
    const std::vector<double> waterlines{0.0, 1.0, 1.8, 2.6, 3.2, 4.2,
                                         5.5, 7.0, 9.0, 11.0, 12.5, kDepth};
    std::vector<Station> stations;
    for (int i = 0; i <= 24; ++i) {
        Station s;
        s.x = -60.0 + 120.0 * i / 24.0;
        const double fx = longitudinalFullness(s.x);
        for (double z : waterlines)
            s.halfBeam.push_back(kHalfBeam * fx * verticalFullness(z));
        stations.push_back(s);
    }
    return makeHullFromStations(stations, waterlines);
}

// A compartment is the hull interior cut by its own bulkhead and deck planes, so
// its volume is whatever the hull form actually leaves between them. The planes
// below are deliberately set wider than the hull at the ends and low down; the
// clip trims them back to the shell rather than the author having to guess where
// the turn of the bilge is.
Compartment carve(const TriMesh& hull, std::string name, Vec3 lo, Vec3 hi,
                  double permeability, bool vented = false) {
    Compartment c;
    c.name = std::move(name);
    c.mesh = clipToBox(hull, lo, hi);
    c.permeability = permeability;
    c.ventedToAtmosphere = vented;
    return c;
}

Opening orifice(std::string name, int a, int b, Vec3 pos, double area, double cd,
                OpeningKind kind, bool open) {
    Opening o;
    o.name = std::move(name);
    o.a = a;
    o.b = b;
    o.pos = pos;
    o.area = area;
    o.dischargeCoeff = cd;
    o.kind = kind;
    o.open = open;
    return o;
}

}  // namespace

Ship buildFerry() {
    Ship ship;
    ship.hull = buildHull();
    ship.deckEdgeZ = 7.0;  // vehicle deck: the level that decides whether she lives

    // --- Subdivision -------------------------------------------------------
    // Lower deck spaces sit on the double bottom at 1.8 m and run to the
    // bulkhead deck at 7.0 m.
    const TriMesh& h = ship.hull;
    ship.compartments = {
        // Stops at the bulkhead deck: above 7 m the space forward is vehicle deck.
        carve(h, "forepeak",          {  44, -20, 1.8}, {  70,  20,  7.0}, 0.97),
        // Outboard limit 8 m, where the wing tanks start -- exactly as the aft
        // holds are drawn below. Authored at 20 m until the breach work found it:
        // that put the whole wing tank *inside* the hold, so 217 m3 of the ship
        // could flood twice. Ship::validate() could not see it, because a ship is
        // only ~89% subdivided and the total never came near the hull.
        carve(h, "fwd_hold_p",        {  20,   0, 1.8}, {  44,   8,  7.0}, 0.95),
        carve(h, "fwd_hold_s",        {  20,  -8, 1.8}, {  44,   0,  7.0}, 0.95),
        carve(h, "engine_room_p",     {  -8,   0, 1.8}, {  20,   8,  7.0}, 0.85),
        carve(h, "engine_room_s",     {  -8,  -8, 1.8}, {  20,   0,  7.0}, 0.85),
        carve(h, "aft_hold_p",        { -38,   0, 1.8}, {  -8,   8,  7.0}, 0.95),
        carve(h, "aft_hold_s",        { -38,  -8, 1.8}, {  -8,   0,  7.0}, 0.95),
        carve(h, "steering_gear",     { -70, -20, 1.8}, { -38,  20,  7.0}, 0.90),
        // Wing tanks now taper into the turn of the bilge instead of being a slab
        // hanging outside the shell, which is what they are on a real ship.
        carve(h, "wing_tank_fwd_p",   {  20,   8, 1.8}, {  44,  20,  7.0}, 0.98),
        carve(h, "wing_tank_fwd_s",   {  20, -20, 1.8}, {  44,  -8,  7.0}, 0.98),
        // The mid pair, flanking the machinery space, and they were missing.
        //
        // Forward and aft wing tanks were authored; these were not, so over the
        // 28 m of engine room the band between y = 8 and the shell at ~9.2 m
        // belonged to no compartment at all. Ramming her amidships tore 63 bays
        // of which **26 opened onto nothing** -- 41% of the hole led nowhere and
        // the ship simply did not flood through it. `breachesFromFailedPanels`
        // reported it every run ("open sea onto a part of the hull that no
        // compartment describes"); nothing was reading the warning.
        //
        // The sharp evidence that it was an omission rather than a design: a
        // strike at either quarter, where the fwd and aft tanks do reach the
        // shell, leaves **zero** panels unmatched. Machinery spaces on a real
        // ro-pax are flanked this way precisely so that side damage floods a void
        // instead of the engine room, which is what the subdivision rules are
        // written around and what the fwd and aft pairs already do here.
        carve(h, "wing_tank_mid_p",   {  -8,   8, 1.8}, {  20,  20,  7.0}, 0.98),
        carve(h, "wing_tank_mid_s",   {  -8, -20, 1.8}, {  20,  -8,  7.0}, 0.98),
        carve(h, "wing_tank_aft_p",   { -38,   8, 1.8}, {  -8,  20,  7.0}, 0.98),
        carve(h, "wing_tank_aft_s",   { -38, -20, 1.8}, {  -8,  -8,  7.0}, 0.98),
        carve(h, "double_bottom_fwd", {   4, -20, 0.0}, {  44,  20,  1.8}, 0.98),
        carve(h, "double_bottom_aft", { -44, -20, 0.0}, {   4,  20,  1.8}, 0.98),
        // One undivided space, 100 m long and 19 m wide, one deck above the
        // waterline. Nothing else on the ship can generate a free surface moment
        // remotely this large.
        carve(h, "vehicle_deck",      { -50, -20, 7.0}, {  50,  20, 12.5}, 0.90, true),
        carve(h, "accommodation",     { -40, -20,12.5}, {  34,  20, 15.0}, 0.95, true),
    };

    const int forepeak   = ship.findCompartment("forepeak");
    const int fwdHoldP   = ship.findCompartment("fwd_hold_p");
    const int fwdHoldS   = ship.findCompartment("fwd_hold_s");
    const int erP        = ship.findCompartment("engine_room_p");
    const int erS        = ship.findCompartment("engine_room_s");
    const int aftHoldP   = ship.findCompartment("aft_hold_p");
    const int aftHoldS   = ship.findCompartment("aft_hold_s");
    const int steering   = ship.findCompartment("steering_gear");
    const int wingFwdP   = ship.findCompartment("wing_tank_fwd_p");
    const int wingFwdS   = ship.findCompartment("wing_tank_fwd_s");
    const int wingMidP   = ship.findCompartment("wing_tank_mid_p");
    const int wingMidS   = ship.findCompartment("wing_tank_mid_s");
    const int wingAftP   = ship.findCompartment("wing_tank_aft_p");
    const int wingAftS   = ship.findCompartment("wing_tank_aft_s");
    const int dbFwd      = ship.findCompartment("double_bottom_fwd");
    const int dbAft      = ship.findCompartment("double_bottom_aft");
    const int vehDeck    = ship.findCompartment("vehicle_deck");

    // --- Flow network ------------------------------------------------------
    ship.openings = {
        // The casualty: a 2.4 m^2 tear in the starboard shell, 2.5 m below the
        // waterline, opening into the starboard engine room.
        //
        // It reaches the machinery space by going *through* the mid wing tank,
        // whose inboard bulkhead it opens as well -- the tank is only about 1.2 m
        // deep here, so anything that tears 2.4 m2 of shell has gone through both.
        // That is the scenario this ship was built to pose and it is deliberately
        // the bad case: a shallower strike would flood the wing tank alone, which
        // is what the tank is for and what makes the difference between a ship
        // that lolls and one that does not.
        orifice("breach_er_s", kSea, erS, {6, -9.0, 3.0}, 2.4, 0.62,
                OpeningKind::Breach, true),

        // Watertight door on the centreline, left open on passage as they so
        // often are. Closing it is the single highest-value action available.
        orifice("wt_door_er", erS, erP, {6, 0.0, 2.2}, 3.6, 0.75,
                OpeningKind::Door, true),

        // Escape trunks up to the vehicle deck. Downflooding paths once the
        // engine rooms fill above 7 m.
        orifice("escape_er_s", erS, vehDeck, {14, -4.0, 7.0}, 1.0, 0.70,
                OpeningKind::Hatch, false),
        orifice("escape_er_p", erP, vehDeck, {14,  4.0, 7.0}, 1.0, 0.70,
                OpeningKind::Hatch, false),

        // Vehicle deck side openings -- freeing ports and the shell doors around
        // the ramp, all sitting just above the deck at 7 m. Dry and irrelevant
        // right up until the ship settles far enough to put them under, at which
        // point they admit the sea onto a 100 x 19 m undivided deck and the free
        // surface moment does the rest.
        orifice("downflood_ramp_s", kSea, vehDeck, { 46, -9.2, 7.2}, 1.2, 0.60,
                OpeningKind::Breach, true),
        orifice("downflood_port_s", kSea, vehDeck, {-10, -9.6, 7.1}, 0.8, 0.60,
                OpeningKind::Breach, true),
        orifice("downflood_port_p", kSea, vehDeck, {-10,  9.6, 7.1}, 0.8, 0.60,
                OpeningKind::Breach, true),

        // An unsealed cable transit through the aft engine room bulkhead. Four
        // hundred square centimetres, and it will flood the aft hold anyway.
        orifice("cable_transit", erS, aftHoldS, {-8, -4.0, 2.0}, 0.04, 0.60,
                OpeningKind::Pipe, true),

        // Air escapes. Generously sized on manned spaces, deliberately mean on
        // the wing tanks, whose small air pipes will trap air and hold them
        // partly buoyant for a long time.
        orifice("vent_er_s",    erS,      kSea, { 10, -6.0, 12.5}, 0.50, 0.80, OpeningKind::Vent, true),
        orifice("vent_er_p",    erP,      kSea, { 10,  6.0, 12.5}, 0.50, 0.80, OpeningKind::Vent, true),
        orifice("vent_fh_p",    fwdHoldP, kSea, { 32,  6.0, 12.5}, 0.20, 0.80, OpeningKind::Vent, true),
        orifice("vent_fh_s",    fwdHoldS, kSea, { 32, -6.0, 12.5}, 0.20, 0.80, OpeningKind::Vent, true),
        orifice("vent_ah_p",    aftHoldP, kSea, {-24,  6.0, 12.5}, 0.20, 0.80, OpeningKind::Vent, true),
        orifice("vent_ah_s",    aftHoldS, kSea, {-24, -6.0, 12.5}, 0.20, 0.80, OpeningKind::Vent, true),
        orifice("vent_steer",   steering, kSea, {-46,  0.0, 12.5}, 0.15, 0.80, OpeningKind::Vent, true),
        orifice("vent_fpk",     forepeak, kSea, { 50,  0.0, 12.5}, 0.10, 0.80, OpeningKind::Vent, true),
        orifice("airpipe_wfp",  wingFwdP, kSea, { 32,  9.0, 12.5}, 0.02, 0.70, OpeningKind::Vent, true),
        orifice("airpipe_wfs",  wingFwdS, kSea, { 32, -9.0, 12.5}, 0.02, 0.70, OpeningKind::Vent, true),
        orifice("airpipe_wmp",  wingMidP, kSea, {  6,  9.0, 12.5}, 0.02, 0.70, OpeningKind::Vent, true),
        orifice("airpipe_wms",  wingMidS, kSea, {  6, -9.0, 12.5}, 0.02, 0.70, OpeningKind::Vent, true),
        orifice("airpipe_wap",  wingAftP, kSea, {-24,  9.0, 12.5}, 0.02, 0.70, OpeningKind::Vent, true),
        orifice("airpipe_was",  wingAftS, kSea, {-24, -9.0, 12.5}, 0.02, 0.70, OpeningKind::Vent, true),
        orifice("airpipe_dbf",  dbFwd,    kSea, { 24,  0.0, 12.5}, 0.02, 0.70, OpeningKind::Vent, true),
        orifice("airpipe_dba",  dbAft,    kSea, {-24,  0.0, 12.5}, 0.02, 0.70, OpeningKind::Vent, true),

        // Counterflooding: sea suction into the port wing tanks, to buy back
        // upright trim at the cost of freeboard.
        orifice("cf_valve_aft_p", kSea, wingAftP, {-24, 9.7, 3.4}, 0.40, 0.65,
                OpeningKind::Pipe, false),
        orifice("cf_valve_fwd_p", kSea, wingFwdP, { 32, 9.7, 3.4}, 0.40, 0.65,
                OpeningKind::Pipe, false),
        // Cross-flooding duct: the passive version, which equalises without a
        // valve to open but drops freeboard on both sides.
        orifice("crossflood_aft", wingAftS, wingAftP, {-24, 0.0, 3.4}, 0.30, 0.85,
                OpeningKind::Pipe, false),
    };

    // --- Pumps -------------------------------------------------------------
    // 216 m^3/h each. Against a 2.4 m^2 breach admitting some 38,000 m^3/h this
    // is not a rounding error away from useless -- it *is* useless, and the sim
    // should make that obvious.
    ship.pumps = {
        {"bilge_er_s", erS,  0.060, 25.0, false, 0.0},
        {"bilge_er_p", erP,  0.060, 25.0, false, 0.0},
        {"bilge_ah_s", aftHoldS, 0.060, 25.0, false, 0.0},
        {"ballast_wing_p", wingAftP, 0.100, 25.0, false, 0.0},
    };

    // --- Weights -----------------------------------------------------------
    // Light displacement is defined as whatever floats her at the design draft
    // with every compartment dry, so the hull form and the loading stay
    // consistent even if the offsets above are edited.
    const double designVolume =
        integrateBelowPlane(ship.hull, Vec3{0, 0, 1}, kDesignDraft).volume;
    ship.lightshipMass = designVolume * kRhoSeawater;
    ship.lightshipCog = {-1.5, 0.0, kKG};
    ship.gyradii = {7.0, 28.0, 29.0};

    return ship;
}

}  // namespace game
