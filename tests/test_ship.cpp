// SPDX-License-Identifier: MIT
//
// Counter-current exchange through a horizontal opening.
//
// A hatch in a deck with the sea standing on it passes water down and air up
// through the same hole at the same time, driven by the density difference and
// not by any net pressure difference at all -- so a model taking one dp at the
// orifice centre reports a hatch at rest while a tonne a second goes through it.
// That is the identical failure `fire.cpp` found in a doorway and fixed by
// integrating over its height; a hatch has no height, so `Ship` splits its area.
//
// Everything here is asserted against the algebra the split is derived from,
// never against the expression `ship.cpp` evaluates: the three equations that
// define the counterflow are solved independently below -- including once by
// bisection, which shares no arithmetic with the closed form at all -- and the
// nulls are asserted as *exact* zeros because they are exact by IEEE arithmetic
// rather than merely small.
#include "engine/core/geometry.hpp"
#include "engine/sim/ship.hpp"
#include "game/prototype/ferry.hpp"
#include "harness.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <tuple>
#include <vector>

using namespace sim;
using testing::expectEqual;
using testing::expectNear;
using testing::expectTrue;

namespace {

// The reference atmosphere, written as `ship.cpp` writes it so the two are the
// same double bit for bit. A test that rounded this would find the exact zeros
// below coming out at 1e-17 and would have to loosen to hide its own arithmetic.
constexpr double kRhoAmb = kPatm / (kRAir * kTAmbient);

// --- Fixtures ----------------------------------------------------------------

// A barge with one sealed space under an open hatch in its deck, and nothing else
// in the network. Small enough that every quantity in it is hand-checkable and
// large enough to float.
//
// The hatch is at the body origin *deliberately*. `R * pos` is then zero at every
// attitude, so the opening sits at the same world point however she is heeled and
// the tilt law below can be asserted with the two sides' pressures held genuinely
// fixed -- otherwise heeling moves the orifice, moves the head over it, and the
// cosine being measured is buried under a change in dp.
struct Glug {
    Ship ship;
    int  hold = 0;
};

Glug makeGlugBarge(double area, double cd, OpeningKind kind, double holdTopZ = -0.5) {
    Glug g;
    Ship& s = g.ship;
    s.hull = makeBox({-20, -8, -6}, {20, 8, 8});
    s.deckEdgeZ = 8.0;

    Compartment c;
    c.name = "hold";
    c.mesh = makeBox({-10, -6, -5}, {10, 6, holdTopZ});
    c.permeability = 1.0;
    s.compartments = {c};

    Opening o;
    o.name = "hatch";
    o.a = kSea;
    o.b = 0;
    o.pos = {0, 0, 0};
    o.area = area;
    o.dischargeCoeff = cd;
    o.kind = kind;
    s.openings = {o};

    // Ballasted to sit with the body origin a little under water, so the sea lies
    // on the hatch and the hold under it is dry: dense on light, which is the one
    // arrangement that drives an exchange.
    s.lightshipMass = 40.0 * 16.0 * 3.5 * kRhoSeawater;
    s.lightshipCog = {0, 0, -2.0};
    s.gyradii = {5.0, 12.0, 12.0};
    s.initialise(0.0);
    return g;
}

// The pressures `Ship::sideStateAt` presents at a hatch on the body origin, from
// public state only. The sea side is a head to the still-water level; the hold's
// gas is at ambient, so its buoyancy head is exactly zero and its pressure at the
// hatch is exactly `airPressure`.
double seaPressureAtOrigin(const Ship& s) {
    return kPatm + s.seaDensity * kGravity * (0.0 - s.state.position.z);
}

// --- The closed form, three ways ---------------------------------------------

// The buoyancy head across a hole of area `area`: the weight of a plug of the
// density difference half a hydraulic diameter deep. The area here is always the
// *authored* one -- a hole does not get smaller when the ship heels -- which is
// why every caller below passes it separately from the flow area.
double buoyancyHead(double area, double rhoUp, double rhoLo) {
    return (rhoUp - rhoLo) * kGravity * 0.5 * std::sqrt(4.0 * area / kPi);
}

// The single-orifice rate this engine gives every *other* opening, on the net head.
// The area split is solved against this, so it is what the two regimes have to
// agree on at the edge of the exchange window.
double singleOrifice(double flowArea, double cd, double rho, double dp) {
    return cd * flowArea * std::sqrt(2.0 * std::abs(dp) / rho);
}

// The counterflow as `ship.cpp` states it: two Torricelli streams sharing one hole
// and passing the net between them.
struct Streams { double down = 0, up = 0; };

Streams closedFormExchange(double flowArea, double dpB, double cd, double rhoUp, double rhoLo,
                           double dpNet) {
    const double vDown = std::sqrt(2.0 * (dpNet + dpB) / rhoUp);
    const double vUp = std::sqrt(2.0 * (dpB - dpNet) / rhoLo);
    const double cdA = cd * flowArea;
    const double qNet = (dpNet >= 0 ? 1.0 : -1.0) *
                        singleOrifice(flowArea, cd, dpNet >= 0 ? rhoUp : rhoLo, dpNet);
    const double aDown = (qNet + cdA * vUp) / (vDown + vUp);
    return {aDown * vDown, (cdA - aDown) * vUp};
}

// The same answer with **no shared arithmetic**: find the area split by bisection
// on the statement that defines it -- the two streams pass the net between them --
// and read the rates straight off. If the closed form has been mis-*derived*
// rather than mis-typed, this is what says so, because it solves the defining
// equation instead of evaluating its solution.
Streams bisectedExchange(double flowArea, double dpB, double cd, double rhoUp, double rhoLo,
                         double dpNet) {
    const double vDown = std::sqrt(2.0 * (dpNet + dpB) / rhoUp);
    const double vUp = std::sqrt(2.0 * (dpB - dpNet) / rhoLo);
    const double cdA = cd * flowArea;
    const double qNet = (dpNet >= 0 ? 1.0 : -1.0) *
                        singleOrifice(flowArea, cd, dpNet >= 0 ? rhoUp : rhoLo, dpNet);
    // phi is the share of the hole carrying the heavy fluid down. The residual is
    // monotone increasing in phi, so bisection is safe and exact to the last bit in
    // about sixty halvings; two hundred is free and leaves no doubt.
    double lo = 0.0, hi = 1.0;
    for (int i = 0; i < 200; ++i) {
        const double phi = 0.5 * (lo + hi);
        const double net = phi * cdA * vDown - (1.0 - phi) * cdA * vUp;
        if (net > qNet) hi = phi; else lo = phi;
    }
    const double phi = 0.5 * (lo + hi);
    return {phi * cdA * vDown, (1.0 - phi) * cdA * vUp};
}

// --- Tests -------------------------------------------------------------------

// The headline: a hatch under water passes both phases at once, at the rate the
// counterflow algebra gives, and the *defining* equations hold on the numbers the
// model published rather than on the ones it was given.
void testHorizontalExchangeMatchesTheCounterflowClosedForm() {
    Glug g = makeGlugBarge(1.0, 0.6, OpeningKind::Hatch);
    Ship& s = g.ship;

    // Pin the attitude and the sinkage so the step below is a statement about the
    // orifice and not about the rigid body underneath it.
    s.state.position.z = -0.30;
    s.state.velocity = {};
    s.state.angularVelocity = {};

    const double pSea = seaPressureAtOrigin(s);
    const double pHold = s.compartments[0].airPressure;
    const double dpNet = pSea - pHold;

    s.step(1e-4, 0.0);

    const Opening& o = s.openings[0];
    const double down = o.lastExchangeDown;
    const double up = o.lastExchangeUp;

    // Non-vacuous first: both streams really ran, and in opposite directions.
    // Without this every assertion below passes on a pair of zeros.
    expectTrue("water is falling through the hatch", down > 0.1);
    expectTrue("and air is rising through the same hatch at the same time", up > 0.1);
    expectTrue("the hold is still dry enough to be an air space",
               s.compartments[0].fillFraction() < 0.01);

    const double dpB = buoyancyHead(o.area, kRhoSeawater, kRhoAmb);
    const Streams want =
        closedFormExchange(o.area, dpB, o.dischargeCoeff, kRhoSeawater, kRhoAmb, dpNet);
    expectNear("the falling stream is the counterflow rate", down, want.down, 1e-12 * want.down);
    expectNear("and so is the rising one", up, want.up, 1e-12 * want.up);

    // And the same pair from the *defining equation* rather than from its solution.
    // These two share nothing but the inputs.
    const Streams bisected =
        bisectedExchange(o.area, dpB, o.dischargeCoeff, kRhoSeawater, kRhoAmb, dpNet);
    expectNear("bisecting the split gives the same falling stream", down, bisected.down,
               1e-12 * want.down);
    expectNear("and the same rising one", up, bisected.up, 1e-12 * want.up);

    // The two defining equations, re-derived from what was *published*. The area
    // each stream got is its own rate over its own Torricelli velocity; the two
    // must add up to the discharge-reduced area of the hole, and the difference of
    // the volumes must be the single-orifice rate on the net head -- which is what
    // makes this model and the one for every other opening the same model.
    const double vDown = std::sqrt(2.0 * (dpNet + dpB) / kRhoSeawater);
    const double vUp = std::sqrt(2.0 * (dpB - dpNet) / kRhoAmb);
    const double aDown = down / vDown;
    const double aUp = up / vUp;
    expectNear("the two streams share exactly one hole",
               (aDown + aUp) / (o.dischargeCoeff * o.area), 1.0, 1e-14);
    const double qNet = singleOrifice(o.area, o.dischargeCoeff, kRhoSeawater, dpNet);
    expectNear("and pass exactly the single-orifice net between them", down - up, qNet,
               1e-12 * qNet);

    // The masses follow the volumes at each side's own density, and nothing else.
    expectNear("the falling stream is water at the sea's own density",
               o.lastExchangeMassDown / down, kRhoSeawater, 1e-12 * kRhoSeawater);
    expectNear("the rising stream is the hold's own air", o.lastExchangeMassUp / up, kRhoAmb,
               1e-12 * kRhoAmb);
    expectTrue("the descending stream came from the sea", o.lastExchangeUpper == kSea);

    // The whole point, stated as the failure it replaces. A model reading one dp at
    // the orifice centre reports the *net*, and the net here is 46% of what falls
    // and none of what rises: the air going the other way is invisible to it
    // entirely, and it is the air leaving that lets the water in.
    std::printf("     one 1 m2 hatch under %.3g m of water: %.4g m3/s down, %.4g m3/s up,"
                " net %.4g -- a net-only model sees the net and misses %.0f%% of the traffic\n",
                dpNet / (kRhoSeawater * kGravity), down, up, down - up,
                100.0 * (1.0 - (down - up) / (down + up)));
    expectTrue("a net rate misses most of what crosses", (down - up) / (down + up) < 0.5);
}

// **The null cases, and they are exact.** A model that ran an exchange whenever it
// saw a hatch would pass every assertion above and be wrong on every one of these.
void testTheExchangeIsExactlyZeroWhereverItShouldNotRun() {
    // (a) Two gas spaces at one temperature. The density difference is +0.0 by IEEE
    //     arithmetic, not a small number, because `buoyancyDensity` is taken
    //     against a reference atmosphere -- so this is asserted at exactly zero.
    {
        Glug g = makeGlugBarge(1.0, 0.6, OpeningKind::Hatch);
        Ship& s = g.ship;
        s.state.position.z = 0.30;   // the hatch is now in air on both sides
        s.state.velocity = {};
        s.step(1e-4, 0.0);
        expectTrue("the hatch really is dry on both sides", !s.openings[0].lastFlowWasWater);
        expectNear("two spaces of air at one temperature exchange exactly nothing",
                   s.openings[0].lastExchangeDown, 0.0, 0.0);
        expectNear("in neither direction", s.openings[0].lastExchangeUp, 0.0, 0.0);
    }

    // (b) Light already lying on heavy: the stable arrangement, which is what a
    //     hatch above a flooded space is. Nothing should overturn.
    {
        Glug g = makeGlugBarge(1.0, 0.6, OpeningKind::Hatch);
        Ship& s = g.ship;
        s.state.position.z = 0.30;                 // air above the hatch
        s.compartments[0].waterVolume = s.compartments[0].floodableVolume() * 0.999;
        s.state.velocity = {};
        s.step(1e-4, 0.0);
        expectTrue("the hold really is full of water", s.compartments[0].fillFraction() > 0.9);
        expectNear("air lying on water overturns exactly nothing",
                   s.openings[0].lastExchangeDown, 0.0, 0.0);
        expectNear("in neither direction", s.openings[0].lastExchangeUp, 0.0, 0.0);
    }

    // (c) Water on water. Same density, so the same exact zero -- and a different
    //     path through the code from (a), because `isWater` is true on both sides.
    {
        Glug g = makeGlugBarge(1.0, 0.6, OpeningKind::Hatch);
        Ship& s = g.ship;
        s.state.position.z = -0.30;                // sea above the hatch
        s.compartments[0].waterVolume = s.compartments[0].floodableVolume() * 0.999;
        s.state.velocity = {};
        s.step(1e-4, 0.0);
        expectNear("water on water exchanges exactly nothing", s.openings[0].lastExchangeDown,
                   0.0, 0.0);
        expectNear("in neither direction", s.openings[0].lastExchangeUp, 0.0, 0.0);
    }

    // (d) The same geometry, the same fluids, the same head -- and the opening
    //     declared as anything but a hatch. This is the control that says the
    //     mechanism is keyed off the *kind* and not off the arrangement, and it is
    //     also the guarantee that nothing already in this engine moved: every
    //     opening on the ferry that is not shut is a breach, a door, a vent or a
    //     pipe, and all four of them take the old path character for character.
    const OpeningKind others[4] = {OpeningKind::Breach, OpeningKind::Door, OpeningKind::Vent,
                                   OpeningKind::Pipe};
    const char* names[4] = {"breach", "door", "vent", "pipe"};
    double vertical[4] = {0, 0, 0, 0};
    for (int k = 0; k < 4; ++k) {
        Glug g = makeGlugBarge(1.0, 0.6, others[k]);
        Ship& s = g.ship;
        s.state.position.z = -0.30;
        s.state.velocity = {};
        s.step(1e-4, 0.0);
        vertical[k] = s.openings[0].lastFlow;
        expectNear(std::string("a ") + names[k] + " under the same water exchanges nothing at all",
                   s.openings[0].lastExchangeDown, 0.0, 0.0);
        expectNear(std::string("nor a ") + names[k] + " in the other direction",
                   s.openings[0].lastExchangeUp, 0.0, 0.0);
        expectTrue(std::string("but the ") + names[k] + " still floods on the net head",
                   vertical[k] > 0);
    }
    // All four vertical kinds are the *same* single-orifice answer, which is what
    // "the old path is untouched" means when it is written as an assertion.
    for (int k = 1; k < 4; ++k)
        expectNear("every non-hatch kind takes the identical single-orifice path", vertical[k],
                   vertical[0], 0.0);
}

// **At rest the two streams are exactly equal volumes, and that is the state a
// net-only model calls a hatch at rest.**
//
// Engineered rather than waited for: the hold's air mass is set so its pressure is
// *bit-identical* to the head standing over the hatch, by nudging the mass through
// neighbouring doubles until the gas law returns the target exactly. A net head of
// even a millipascal would leave a residual net flow and turn this exact statement
// into an approximate one.
void testAtRestTheTwoStreamsAreExactlyEqualVolumes() {
    Ship s;
    s.hull = makeBox({-20, -8, -6}, {20, 8, 10});
    s.deckEdgeZ = 10.0;
    auto room = [&](const char* name, double z0, double z1, bool vented) {
        Compartment c;
        c.name = name;
        c.mesh = makeBox({-10, -6, z0}, {10, 6, z1});
        c.permeability = 1.0;
        c.ventedToAtmosphere = vented;
        return c;
    };
    s.compartments = {room("hold", -4, 0, false), room("deck", 0, 6, true)};
    s.compartments[1].waterVolume = 24.0;      // 0.2 m over a 10 x 12 m floor

    Opening o;
    o.name = "hatch";
    o.a = 0;
    o.b = 1;
    o.pos = {0, 0, 0};
    o.area = 0.9;
    o.dischargeCoeff = 0.65;
    o.kind = OpeningKind::Hatch;
    s.openings = {o};

    s.lightshipMass = 40.0 * 16.0 * 3.0 * kRhoSeawater;
    s.lightshipCog = {0, 0, -1.0};
    s.gyradii = {5.0, 12.0, 12.0};
    s.initialise(0.0);

    // One step first, and it is not a formality: `solveOffsetForVolume` bisects
    // from the offset it was last left at, so the deck's free surface -- and with
    // it the head this test is about to match to the last bit -- only settles on
    // its final value after the solve has been run from a converged seed.
    s.step(1e-9, 0.0);

    Compartment& hold = s.compartments[0];
    const Compartment& deck = s.compartments[1];
    const double zHatch = s.state.position.z;              // the hatch is on the origin
    const double target = kPatm + s.seaDensity * kGravity * (deck.surfaceWorldZ - zHatch);
    const double va = hold.airVolume();
    // Solve `m R T / va == target` in doubles rather than in reals. The closed-form
    // mass is within an ulp or two; walk to whichever neighbour reproduces `target`
    // bit for bit, so the net head below is +0.0 and not 1e-11.
    double m = target * va / (kRAir * hold.gasTemperature);
    for (int i = 0; i < 8 && m * kRAir * hold.gasTemperature / va != target; ++i)
        m = std::nextafter(m, m * kRAir * hold.gasTemperature / va > target ? 0.0 : 2.0 * m);
    hold.airMass = m;

    // **Exactly is not available here, and the reason is arithmetic rather than
    // modelling.** The gas law scales a mass of 1.3e3 kg into a pressure of 1e5 Pa
    // by about 86, so one ulp of mass moves the pressure by 2e-11 Pa while one ulp
    // of pressure is 1.5e-11: the reachable pressures are *sparser* than the
    // target grid, and the walk above lands on the nearest attainable neighbour.
    // What is asserted is that it is a neighbour -- and then the exact statement is
    // made against the residual head that leaves, rather than against a zero the
    // representation cannot hold.
    const double achieved = hold.airMass * kRAir * hold.gasTemperature / va;
    const double ulpP = std::nextafter(target, 2.0 * target) - target;
    expectTrue("the hold sits within one ulp of the head standing over the hatch",
               std::abs(achieved - target) <= 2.0 * ulpP);

    s.step(1e-4, 0.0);
    const Opening& ran = s.openings[0];
    expectTrue("water is falling into the hold", ran.lastExchangeDown > 0.1);
    expectTrue("and air is rising out of it", ran.lastExchangeUp > 0.1);
    expectTrue("the deck is the upper side", ran.lastExchangeUpper == 1);

    // The closed form at zero net head, where the two Torricelli streams reduce to
    // one expression: `Cd A sqrt(2 dp_b) / (sqrt(rho_up) + sqrt(rho_lo))`.
    //
    // **`rho_lo` is the hold's real gas density and not the reference atmosphere,
    // and this test asserted the wrong one first.** The two densities in
    // `SideState` are there precisely because they answer different questions: the
    // buoyancy head is a *difference*, taken against a reference atmosphere so that
    // two cold spaces give exactly +0.0; the Torricelli denominator is a *local*
    // density, `p/(R T)`, and this hold has been pressurised to 2 kPa above
    // atmosphere to put the net head at zero. Using `kRhoAmb` in both places is
    // wrong by 1.65e-4 here -- small, entirely systematic, and exactly the kind of
    // thing a loose tolerance would have hidden.
    const double rhoLo = achieved / (kRAir * hold.gasTemperature);
    const double dpB = buoyancyHead(o.area, kRhoSeawater, kRhoAmb);
    const double want = o.dischargeCoeff * o.area * std::sqrt(2.0 * dpB) /
                        (std::sqrt(kRhoSeawater) + std::sqrt(rhoLo));
    expectNear("the falling stream is the at-rest counterflow rate", ran.lastExchangeDown, want,
               1e-6 * want);

    // The net it carries is the single-orifice rate on the residual head, which is
    // the model's own defining equation and is exact. That is the statement "no net
    // volume" turns into once the representation is admitted: 9.9e-8 m3/s out of
    // 1.8, from 1.5e-11 Pa of head that a double cannot spend.
    const double residual = std::abs(achieved - target);
    const double qResidual = singleOrifice(o.area, o.dischargeCoeff, kRhoSeawater, residual);
    std::printf("     at rest: %.6g m3/s down, %.6g m3/s up, net %.3g m3/s from %.2g Pa of"
                " residual head a double cannot spend\n",
                ran.lastExchangeDown, ran.lastExchangeUp,
                ran.lastExchangeDown - ran.lastExchangeUp, residual);
    expectTrue("and it carries no net volume beyond that residual head",
               std::abs(ran.lastExchangeDown - ran.lastExchangeUp) <= qResidual + 1e-15);
    expectTrue("which is nine orders of magnitude under what crosses",
               std::abs(ran.lastExchangeDown - ran.lastExchangeUp) < 1e-7 * ran.lastExchangeDown);

    // Which is the whole failure being fixed, said as a number: the net a
    // single-dp model reports is zero, and tonnes a second are going through.
    expectTrue("a net-only model would have reported this hatch as still",
               ran.lastExchangeMassDown > 1000.0);
}

// **The two regimes have to be the same model, and at the edge of the window they
// are the same number.**
//
// Past `|dp_net| >= dp_b` the hole is flushed one way and the single-orifice solve
// takes over. Sweeping the net head across that edge is what caught the first
// version of this work: requiring the two *volumes* to be equal drove the exchange
// to zero as the edge was approached, so the model reported a hole passing nothing
// at all in a window under a millimetre of water wide, and then jumped to 2.3 m^3/s
// the instant the net took over. Solving the split against the net makes the two
// agree identically, and this is the test that says so.
// **And the cusp they meet at has a closed form of its own.** The approach is not
// linear: the rising stream's velocity goes as the square root of how far inside
// the window the net head is, so the gap closes as sqrt(eps) rather than as eps.
// Working the split out at the edge gives
//
//     q_exchange - q_single  =  Cd A (1 - 1/sqrt(2)) sqrt(2 eps / rho_lo)
//
// which is the same square root with the infinite derivative that makes `fire.cpp`
// split its bands *at* the neutral plane rather than quadrature through it.
// Asserting that, at four values of eps spanning six decades, is a far stronger
// statement than any single "the jump is small": it says the two regimes are the
// same model *and* names the rate at which they become it.
void testTheExchangeMeetsTheSingleOrificeSolveAtTheEdgeOfItsWindow() {
    const double dpB = buoyancyHead(1.0, kRhoSeawater, kRhoAmb);
    const double eps[4] = {1e-1, 1e-3, 1e-5, 1e-7};
    for (int k = 0; k < 4; ++k) {
        double flow[2] = {0, 0};
        bool exchanged[2] = {false, false};
        for (int side = 0; side < 2; ++side) {
            Glug g = makeGlugBarge(1.0, 0.6, OpeningKind::Hatch);
            Ship& s = g.ship;
            s.state.velocity = {};
            s.state.angularVelocity = {};
            // Sink her to the depth whose head is dp_b either side of the edge, so
            // one sample runs the exchange and the other the single-orifice solve.
            s.state.position.z =
                -(dpB + (side == 0 ? -eps[k] : eps[k])) / (s.seaDensity * kGravity);
            s.step(1e-6, 0.0);
            flow[side] = s.openings[0].lastFlow;   // + is sea -> hold, in both regimes
            exchanged[side] = s.openings[0].lastExchangeDown > 0;
        }
        // Non-vacuous: this pair really did straddle the two regimes.
        expectTrue("the inside sample ran the exchange", exchanged[0]);
        expectTrue("and the outside one ran the single-orifice solve", !exchanged[1]);
        expectTrue("both regimes are passing water", flow[0] > 1.0 && flow[1] > 1.0);

        const double want = 0.6 * 1.0 * (1.0 - 1.0 / std::sqrt(2.0)) *
                            std::sqrt(2.0 * eps[k] / kRhoAmb);
        const double got = flow[0] - flow[1];
        std::printf("     cusp at %6.0e Pa inside the edge: %.5g m3/s of gap against a leading"
                    " term of %.5g (%.2e of the flow)\n", eps[k], got, want, got / flow[1]);
        // **The tolerance goes as eps and not as a constant fraction**, because the
        // closed form above is the *leading* term of an expansion in sqrt(eps) and
        // what it drops is O(eps). Measured, the absolute residual is 5.7e-3 at
        // eps = 0.1 and 6.1e-5 at 1e-3 -- a factor of 93 for a factor of 100 in
        // eps, which is that O(eps) and not anything else. A single relative
        // tolerance loose enough to pass the first sample would be forty times
        // slacker than the last one needs.
        expectNear("the gap at the edge is Cd A (1 - 1/sqrt 2) sqrt(2 eps / rho_lo)", got, want,
                   0.08 * eps[k] + 1e-12);
    }
}

// **Trapped air arrests a breach. It does not arrest a hatch.**
//
// The behaviour the front page publishes -- a sealed space stops flooding when its
// air pressure balances the outside head, which is why upside-down hulls float for
// hours -- is a statement about *vertical* openings. Through a horizontal one the
// air has a way out the water is not blocking, so the arrest comes later and
// deeper, and the extra head it takes is exactly the buoyancy head of the hole.
//
// Two barges identical in every respect but the kind of the opening, each run to a
// steady state and then run the same distance again, so "arrested" is a *measured*
// statement about the second window and not an eyeballed one about the first.
//
// **The expectation this test started with was wrong, and the way it was wrong is
// worth keeping.** It predicted the hatch arresting too, just deeper -- at the head
// plus the hole's own buoyancy head, where the net flow reverses. It does not
// arrest at all, and the reason is that the pressure has no second equation to
// satisfy. Differentiating `p = m R T / (V0 - W)` gives
//
//     dp/dt = p (q_down - q_up) / V_air
//
// so the pressure is stationary exactly when the two *volumes* balance, which is
// where the net head is zero -- and at zero net head the exchange is at its
// strongest, not its weakest. The steady state is therefore water pouring in at the
// full at-rest counterflow rate with the pressure pinned to the outside head, and
// it persists until the space is full. A vertical opening has no such state: there
// the same condition, zero net head, is the condition for no flow at all.
void testTrappedAirArrestsABreachButNotAHatch() {
    auto settle = [](OpeningKind kind, int steps) {
        Glug g = makeGlugBarge(1.0, 0.6, kind);
        Ship& s = g.ship;
        for (int i = 0; i < steps; ++i) {
            s.state.position.z = -0.30;    // hold the head; only the opening may act
            s.state.velocity = {};
            s.state.angularVelocity = {};
            s.step(2e-3, 0.0);
        }
        return g;
    };
    const double head = kRhoSeawater * kGravity * 0.30;

    // --- the control: a vertical hole of the same size, in the same barge --------
    Glug b1 = settle(OpeningKind::Breach, 15000);   // 30 s
    Glug b2 = settle(OpeningKind::Breach, 30000);   // 60 s
    const double vBoyle =
        b1.ship.compartments[0].floodableVolume() * (1.0 - kPatm / (kPatm + head));
    expectNear("the breached barge stops where Boyle's law says it must",
               b2.ship.compartments[0].waterVolume, vBoyle, 0.02 * vBoyle);
    expectNear("its air stands at exactly the head over the hole",
               b2.ship.compartments[0].airPressure - kPatm, head, 0.02 * head);
    // Arresting is a statement about the *rate*, not about the total, and asserting
    // the total stopped moving would have been asserting convergence rather than
    // physics: measured, the breach is still creeping the last 9% of the way in at
    // 30 s. What has happened by then is that it decelerated by a factor of ten,
    // which is what "the trapped air is holding the sea out" means.
    const double bRate1 = b1.ship.compartments[0].waterVolume / 30.0;
    const double bRate2 = (b2.ship.compartments[0].waterVolume -
                           b1.ship.compartments[0].waterVolume) / 30.0;
    expectTrue("and it had all but stopped by thirty seconds", bRate2 < 0.15 * bRate1);

    // --- the hatch: same hole, same head, no arrest ------------------------------
    Glug h1 = settle(OpeningKind::Hatch, 15000);
    Glug h2 = settle(OpeningKind::Hatch, 30000);
    const double w1 = h1.ship.compartments[0].waterVolume;
    const double w2 = h2.ship.compartments[0].waterVolume;

    // The pressure it settles at is the outside head and nothing more -- not the
    // head plus the buoyancy head, which is where the *flow* would reverse.
    expectNear("the glugging space settles at exactly the head over the hole",
               h2.ship.compartments[0].airPressure - kPatm, head, 0.02 * head);

    // And it is still taking water at the at-rest counterflow rate, which is the
    // closed form for zero net head -- so the second thirty seconds admit as much
    // as a steady rate says they must.
    const double dpB = buoyancyHead(1.0, kRhoSeawater, kRhoAmb);
    const double steady = 0.6 * 1.0 * std::sqrt(2.0 * dpB) /
                          (std::sqrt(kRhoSeawater) + std::sqrt(kRhoAmb));
    std::printf("     the same 1 m2 hole over 60 s: a breach takes %.1f m3 and has decelerated"
                " %.0fx; a hatch takes %.1f m3 and is still running at %.3f m3/s against the"
                " at-rest closed form of %.3f\n",
                b2.ship.compartments[0].waterVolume, bRate1 / bRate2, w2, (w2 - w1) / 30.0,
                steady);
    expectNear("the hatch is still glugging at the at-rest rate after a minute",
               (w2 - w1) / 30.0, steady, 0.02 * steady);

    // Non-vacuous, and the finding: the two are not close, and the hatch has not
    // slowed down at all.
    expectTrue("a horizontal opening defeats the trapped air a vertical one holds",
               w2 > 3.0 * b2.ship.compartments[0].waterVolume);
}

// The exchange is a *vertical* counterflow, so what it has to work with is the
// opening projected on the horizontal. Heel her and that goes as the cosine of the
// deck's tilt, exactly -- and at ninety degrees a deck hatch is a scuttle in a wall
// with no height for a neutral plane to sit in, so it carries nothing.
//
// Asserted as a ratio at four angles rather than as four rates, because the hatch
// sits on the body origin: `R * pos` is zero at every attitude, both sides' pressures
// are then identical between the runs, and the only thing that can move the answer
// is the projection itself.
//
// **The reference is the deck normal the ship is actually holding, not
// `cos(heel)`, and the difference is not pedantry.** `Quat::fromAxisAngle` carries
// a right angle as a half-angle sine and cosine, so its matrix reports the deck
// normal's world-z component as 2.2e-16 where `std::cos(pi/2)` is 6.1e-17 -- three
// and a half times apart, both of them "zero". Asserting against `std::cos` at a
// tolerance loose enough to swallow that would swallow a broken projection too, so
// the law is asserted against the ship's own normal at *zero* tolerance and the
// cosine is checked separately at the angles a quaternion can resolve.
void testTheExchangeGoesAsTheCosineOfTheDeckSTilt() {
    const double angles[4] = {0.0, kPi / 6.0, kPi / 3.0, kPi / 2.0};
    double q[4] = {0, 0, 0, 0};
    double tilt[4] = {0, 0, 0, 0};
    for (int k = 0; k < 4; ++k) {
        Glug g = makeGlugBarge(1.0, 0.6, OpeningKind::Hatch);
        Ship& s = g.ship;
        s.state.orientation = Quat::fromAxisAngle(Vec3{1, 0, 0}, angles[k]);
        s.state.position.z = -0.30;
        s.state.velocity = {};
        s.state.angularVelocity = {};
        // Read off the ship before the step, because the step ends by integrating
        // the rigid body and a barge held at sixty degrees is under an enormous
        // moment: the attitude the orifice was solved at is this one and not the
        // one left behind afterwards.
        tilt[k] = (s.state.orientation.toMat3() * Vec3{0, 0, 1}).z;
        s.step(1e-4, 0.0);
        q[k] = s.openings[0].lastExchangeDown;
    }
    expectTrue("the upright case is not zero, so the ratios below mean something", q[0] > 0.1);

    // The law itself: rate over projected fraction is one constant across four
    // attitudes spanning sixteen orders of magnitude in the projection.
    for (int k = 1; k < 4; ++k)
        expectNear("the exchange is exactly the upright rate times the deck's own projection",
                   (q[k] / tilt[k]) / (q[0] / tilt[0]), 1.0, 1e-15);

    // And that projection is the cosine of the heel, at the angles where a
    // quaternion has the resolution to say so.
    for (int k = 1; k < 3; ++k)
        expectNear("the projection is the cosine of the heel", tilt[k], std::cos(angles[k]), 1e-15);

    // The null. A deck hatch heeled through a right angle is a scuttle in a wall,
    // and `Opening` has no height for a neutral plane to sit in, so it carries
    // nothing: 2.2e-16 of the upright rate, which is the quaternion's own idea of
    // zero and not a number this model invented.
    expectTrue("and at ninety degrees it has vanished", q[3] / q[0] < 1e-15);
    expectTrue("it vanished because the deck did, not for some other reason",
               tilt[3] < 1e-15);
}

// **Mass in equals mass out.** The barge above exchanges with the sea, so the sum
// inside is not conserved -- what is conserved is the sum inside less everything
// the openings reported carrying across the boundary, and the two must agree to
// round-off because each transfer is applied to both ends of the same edge.
//
// Run on a two-compartment fixture where the exchange is *internal*, so the ledger
// closes with nothing crossing a boundary at all and the assertion cannot be
// satisfied by a large flux cancelling a large error.
void testMassIsConservedAcrossTheExchange() {
    Ship s;
    s.hull = makeBox({-20, -8, -6}, {20, 8, 10});
    s.deckEdgeZ = 10.0;
    auto room = [&](const char* name, double z0, double z1) {
        Compartment c;
        c.name = name;
        c.mesh = makeBox({-10, -6, z0}, {10, 6, z1});
        c.permeability = 1.0;
        return c;
    };
    // Upper holds a shallow layer of water; lower is dry and sealed. The only edge
    // between them is the hatch on their shared deck.
    s.compartments = {room("lower", -4, 0), room("upper", 0, 6)};
    s.compartments[1].waterVolume = 12.0;      // 0.1 m over a 10 x 12 m floor

    Opening o;
    o.name = "hatch";
    o.a = 0;
    o.b = 1;
    o.pos = {0, 0, 0};
    o.area = 0.8;
    o.dischargeCoeff = 0.65;
    o.kind = OpeningKind::Hatch;
    s.openings = {o};

    s.lightshipMass = 40.0 * 16.0 * 3.0 * kRhoSeawater;
    s.lightshipCog = {0, 0, -1.0};
    s.gyradii = {5.0, 12.0, 12.0};
    s.initialise(0.0);

    auto total = [&] {
        double m = 0;
        for (const Compartment& c : s.compartments) m += c.waterVolume * s.seaDensity + c.airMass;
        return m;
    };
    const double m0 = total();
    double worst = 0.0;
    double movedWater = 0.0, movedGas = 0.0;
    const double dt = 1e-3;
    for (int i = 0; i < 1000; ++i) {
        s.step(dt, 0.0);
        movedWater += s.openings[0].lastExchangeMassDown * dt;
        movedGas += s.openings[0].lastExchangeMassUp * dt;
        worst = std::max(worst, std::abs(total() - m0));
    }
    // Non-vacuous from both ends: a real exchange ran, and it was two-way.
    expectTrue("water really fell through the hatch", movedWater > 100.0);
    expectTrue("and air really rose through the same hatch", movedGas > 0.1);
    expectTrue("the water that moved is a real share of the layer there was",
               movedWater > 0.05 * 12.0 * kRhoSeawater);

    // The assertion sits at what was measured and not at an order of magnitude
    // above it: a model that had lost the property would be out by the mass it
    // moved -- a tonne and a half through this hatch -- so any comfortable
    // tolerance would pass on one. The figure is printed rather than only written
    // in this comment, because nothing tests a comment.
    std::printf("     internal exchange: %.4g kg of water down, %.4g kg of air up,"
                " total mass held to %.3g of %.4g kg (%.2e)\n",
                movedWater, movedGas, worst, m0, worst / m0);
    expectNear("no mass is created or destroyed by the exchange", worst / m0, 0.0, 2.5e-15);
}

// The energy the two gas streams carry, when both of them are gas.
//
// A single `lastGasMassFlow` cannot describe two donors at two temperatures, which
// is why `lastGasEnthalpyFlow` exists. Assert that it does add up: a hot space
// under a cold one exchanges through the hatch between them, and the network's own
// energy has to move by exactly what the field says crossed.
void testTwoGasStreamsBalanceTheirOwnEnthalpy() {
    Ship s;
    s.hull = makeBox({-20, -8, -6}, {20, 8, 10});
    s.deckEdgeZ = 10.0;
    auto room = [&](const char* name, double z0, double z1) {
        Compartment c;
        c.name = name;
        c.mesh = makeBox({-10, -6, z0}, {10, 6, z1});
        c.permeability = 1.0;
        c.gasThermalTime = std::numeric_limits<double>::infinity();  // adiabatic
        return c;
    };
    // Cold, dense gas on top of hot, light gas: unstable, and it overturns.
    s.compartments = {room("hot", -4, 0), room("cold", 0, 6)};
    s.compartments[0].gasTemperature = 900.0;

    Opening o;
    o.name = "hatch";
    o.a = 0;
    o.b = 1;
    o.pos = {0, 0, 0};
    o.area = 0.5;
    o.dischargeCoeff = 0.6;
    o.kind = OpeningKind::Hatch;
    s.openings = {o};

    s.lightshipMass = 40.0 * 16.0 * 3.0 * kRhoSeawater;
    s.lightshipCog = {0, 0, -1.0};
    s.gyradii = {5.0, 12.0, 12.0};
    s.initialise(0.0);

    auto energy = [&] {
        double e = 0;
        for (const Compartment& c : s.compartments) e += c.airMass * c.gasTemperature;
        return e;
    };
    const double e0 = energy();
    double worst = 0.0, exchanged = 0.0;
    const double dt = 1e-3;
    for (int i = 0; i < 500; ++i) {
        s.step(dt, 0.0);
        exchanged += s.openings[0].lastExchangeMassUp * dt;
        worst = std::max(worst, std::abs(energy() - e0));
    }
    expectTrue("hot gas really did rise through the hatch", exchanged > 0.01);
    expectTrue("and the two sides really were at different temperatures",
               s.compartments[0].gasTemperature > 1.5 * s.compartments[1].gasTemperature);
    expectNear("no single donor temperature is offered for two streams",
               s.openings[0].lastGasDonorTemperature, 0.0, 0.0);
    // Both endpoints are inside the ledger, so every joule the hatch takes from one
    // it gives to the other and the total cannot move at all -- the enthalpy term
    // is antisymmetric across the opening whatever the clamps do.
    std::printf("     two gas streams: %.4g kg of hot gas up, network energy held to %.3g"
                " of %.4g kg K (%.2e)\n",
                exchanged, worst, e0, worst / e0);
    expectNear("an internal exchange moves no energy out of the network", worst / e0, 0.0, 2e-15);

    // And the field that had to exist for this to be checkable at all: the two
    // streams' enthalpies, summed as one signed rate, against the same two masses
    // and temperatures taken apart. A single donor temperature cannot say this.
    //
    // The temperatures are taken *before* the step that fills the field, because
    // the step ends by mixing what arrived into what was there.
    const double tLower = s.compartments[0].gasTemperature;   // rises: the light one
    const double tUpper = s.compartments[1].gasTemperature;   // falls: the dense one
    s.step(dt, 0.0);
    const Opening& ran = s.openings[0];
    // a is the lower space, so a -> b is upward and the rising stream carries the
    // positive sign.
    const double want =
        kGammaAir * (tLower * ran.lastExchangeMassUp - tUpper * ran.lastExchangeMassDown);
    const double stream = kGammaAir * tLower * ran.lastExchangeMassUp;
    expectTrue("both streams carried real enthalpy", stream > 1.0);

    // **The net is 4.5e-5 of either stream, and that is a closed form rather than
    // an accident.** Equal volumes of ideal gas at one pressure carry equal
    // enthalpy whatever their temperatures -- `gamma T rho V` is `gamma p V / R` --
    // so a volume-balanced exchange between two spaces at the same pressure moves
    // no enthalpy at all, and what survives is `gamma Q dp / R`. The two terms
    // above therefore cancel to four and a half decimal places by construction,
    // which is why this is asserted against the size of a *stream* and not against
    // the size of the difference: a relative tolerance on the difference would be
    // a tolerance of 1e-12 on a quantity known to 2e-12, and would fail on correct
    // code for reasons that have nothing to do with the model.
    std::printf("     enthalpy: net %.6g against %.6g in either stream (%.2e), residual %.2e\n",
                ran.lastGasEnthalpyFlow, stream, std::abs(want) / stream,
                std::abs(ran.lastGasEnthalpyFlow - want) / stream);
    expectNear("the published enthalpy rate is the two streams, taken apart",
               ran.lastGasEnthalpyFlow, want, 2.3e-16 * stream);
    expectTrue("and it is far smaller than either of them, as an ideal gas requires",
               std::abs(ran.lastGasEnthalpyFlow) < 0.2 * stream);
}

// --- The ferry ---------------------------------------------------------------

// **What this cost the ship, which is nothing, and it is asserted rather than
// asserted-about.** Every hatch the ferry has is an escape trunk from an engine
// room to the vehicle deck, and every one of them is shut in `ships/ferry.ship`
// and in all three scenarios. So no published figure can move, and the way to say
// that in a test is to run her and find no opening reporting an exchange.
//
// It is a null, so it needs the companion below or it proves nothing.
void testTheFerryAsShippedRunsNoExchangeAtAll() {
    Ship s = game::buildFerry();
    s.initialise(0.0);
    int hatches = 0;
    for (const Opening& o : s.openings)
        if (o.kind == OpeningKind::Hatch) ++hatches;
    expectTrue("the ferry has hatches for this to have been able to fire", hatches == 2);

    double worst = 0.0;
    for (int i = 0; i < 6000; ++i) {
        s.step(0.01, 0.0);
        for (const Opening& o : s.openings)
            worst = std::max(worst, std::max(o.lastExchangeDown, o.lastExchangeUp));
    }
    expectTrue("she flooded, so the run was not a no-op", s.totalFloodwaterMass() > 1e5);
    expectNear("no opening on the ferry as shipped exchanges anything", worst, 0.0, 0.0);
}

// The companion, and the size of the effect on her own compartments: leave an
// escape trunk open with water on the vehicle deck above it and a dry engine room
// below, which is the state the trunk exists to be a hazard in.
//
// The number this measures -- how much a 1 m^2 hatch passes at rest -- is the whole
// case for the mechanism, so it is asserted against the closed form and not merely
// observed to be non-zero.
void testAnOpenEscapeTrunkDrainsTheVehicleDeckIntoADryEngineRoom() {
    Ship s = game::buildFerry();
    // No hull damage: the only thing happening is the hatch.
    for (Opening& o : s.openings)
        if (o.kind == OpeningKind::Breach) o.open = false;
    for (Opening& o : s.openings)
        if (o.name == "escape_er_p") o.open = true;
    const int deck = s.findCompartment("vehicle_deck");
    const int erP = s.findCompartment("engine_room_p");
    expectTrue("the ferry has the two spaces this test uses", deck != kSea && erP != kSea);
    // Enough water on the deck to cover the trunk head at y = +4 m, and shallow
    // enough that its head stays under the buoyancy head -- past that the hole is
    // flushed one way and the exchange correctly stops.
    s.compartments[static_cast<std::size_t>(deck)].waterVolume = 120.0;
    s.initialise(0.0);

    int hatch = -1;
    for (std::size_t i = 0; i < s.openings.size(); ++i)
        if (s.openings[i].name == "escape_er_p") hatch = static_cast<int>(i);

    // The head, the attitude and the two pressures as they stand *before* the step,
    // because those are what the orifice was solved at: `step` re-levels every free
    // surface and integrates the rigid body afterwards, and reading them back from
    // the far side of that is reading a different ship.
    const Compartment& vd = s.compartments[static_cast<std::size_t>(deck)];
    const Compartment& er = s.compartments[static_cast<std::size_t>(erP)];
    const Mat3 R = s.state.orientation.toMat3();
    const Opening& o = s.openings[static_cast<std::size_t>(hatch)];
    const Vec3 wp = R * o.pos + s.state.position;
    const double pUp = vd.airPressure + s.seaDensity * kGravity * (vd.surfaceWorldZ - wp.z);
    const double pLo = er.airPressure;
    const double tilt = (R * Vec3{0, 0, 1}).z;

    s.step(1e-3, 0.0);
    expectTrue("the trunk is running", o.lastExchangeDown > 0.1);
    expectTrue("and it is draining the deck into the engine room", o.lastExchangeUpper == deck);

    // The buoyancy head comes off the authored area -- the hole is the size it was
    // built -- while the flow area is what still faces up at this heel. Passing the
    // two separately is the whole distinction, and getting it wrong would show up
    // here as a quarter power of the heel.
    const double dpB = buoyancyHead(o.area, kRhoSeawater, kRhoAmb);
    const Streams want = closedFormExchange(o.area * std::abs(tilt), dpB, o.dischargeCoeff,
                                            kRhoSeawater, kRhoAmb, pUp - pLo);
    expectNear("the trunk passes the counterflow rate for the head it is under",
               o.lastExchangeDown, want.down, 1e-12 * want.down);
    expectNear("and vents the engine room's air back up the same trunk", o.lastExchangeUp, want.up,
               1e-12 * want.up);

    // The size of it, which is the answer to "does this earn its place": a single
    // 1 m^2 trunk at rest moves more than a tonne of water a second, against a
    // model that reported the hatch as still.
    expectTrue("a 1 m2 escape trunk passes over a tonne a second at rest",
               o.lastExchangeMassDown > 1000.0);
}

}  // namespace

// **A ship with no mass must not step into a NaN, and must not be reported sunk.**
//
// `lightshipMass` defaults to 0, and `Ship::step` does not refuse -- so the
// effective mass was exactly zero on a dry hull, `fBody` was zero with it, and the
// linear acceleration was `0/0`. Nothing downstream reads that as a NaN: the
// orientation's normalisation is `if (l > 1e-300)`, false for a NaN, so the
// quaternion stays NaN for good; `heelTrimFromRotation` then gives a NaN heel, and
// `Diagnostics`'s `std::abs(d.heelDeg) < 90.0` is false for a NaN too. The reported
// answer was `afloat = false`. A misconfiguration arrived as a sinking.
//
// The rotational half of the same line has always been defended -- `inverse`
// returns a zero matrix on a singular inertia -- and this asserts the two halves
// now agree.
void testAMasslessShipStepsToFiniteStateRatherThanNaN() {
    Ship ship;
    ship.hull = game::buildFerry().hull;
    ship.initialise(0.0);
    expectNear("the fixture really has no mass, which is the whole point",
               ship.lightshipMass, 0.0, 0.0);

    for (int i = 0; i < 20; ++i) ship.step(0.02, Sea(0.0));

    const Vec3 p = ship.state.position;
    expectTrue("position stays finite",
               std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z));
    const Diagnostics d = ship.diagnostics(Sea(0.0));
    expectTrue("heel stays finite", std::isfinite(d.heelDeg));
    expectTrue("and trim", std::isfinite(d.trimDeg));
    // The point of the whole test: a NaN heel makes `afloat` false, so a finite
    // heel is what stops a missing mass from being reported as a sinking.
    expectTrue("so she is not reported sunk on account of a NaN", d.afloat);
}

// **A name the ship does not carry must be refusable, not answered with the sea.**
//
// `findCompartment`'s one-argument form answers `kSea` for a miss, and kSea is not
// a sentinel -- it is a real endpoint of the flow network, the open sea. Seventeen
// unchecked lookups in the ferry's builder were seventeen chances for a renamed or
// mistyped compartment to become a hole in the shell that then floods at a
// completely plausible rate. The refusable form is what makes a miss visible, and
// it refuses by leaving the caller's index alone rather than by writing an index
// of its own: there is no value it could write that some caller would not read as
// a space.
void testAMissedCompartmentLookupIsRefusedRatherThanAnsweredWithTheSea() {
    const Ship ship = game::buildFerry();

    // Against a vacuous test: the ferry really does carry one of these names and
    // really does not carry the other, and the two differ the way a typo differs.
    int found = kNoCompartment;
    expectTrue("the ferry carries engine_room_s", ship.findCompartment("engine_room_s", found));
    expectTrue("... at a compartment that exists",
               found >= 0 && found < static_cast<int>(ship.compartments.size()));

    constexpr int kUntouched = -12345;
    int index = kUntouched;
    expectTrue("a name it does not carry is refused",
               !ship.findCompartment("engine_room_starboard", index));
    expectEqual("... leaving the caller's index alone rather than writing the sea into it",
                index, kUntouched);

    // The one-argument form still answers kSea, which is the contract
    // `shipfile.cpp`'s resolve() refuses on -- and is exactly why nothing else may
    // use it without testing for kSea on the next line.
    expectEqual("the one-argument form still answers the open sea, as its one caller expects",
                ship.findCompartment("engine_room_starboard"), kSea);
    expectEqual("and the two forms agree wherever the name is there",
                ship.findCompartment("engine_room_s"), found);
}

// **The two answers lead to opposite outcomes, which is what makes them two
// sentinels rather than two spellings of one.**
//
// `Ship::validate()` has to let an opening with one end in the sea past -- the
// ferry's breach, its three downflooding openings, its sixteen vents and air pipes
// and its two counterflooding valves are all exactly that -- so a mistyped name
// resolved to kSea produces a definition that validates clean and floods, and 22
// of the ferry's 54 endpoints are legitimately indistinguishable from one. The
// identical opening built from
// `kNoCompartment` is named, by the opening that carries it. Asserted as a pair,
// because either half alone is consistent with the defect being present.
void testAnEndpointFromAMissedLookupIsNamedWhereASeaEndpointCannotBe() {
    // Initialised, because `validate()`'s overlap check reads the gross volumes
    // `initialise` integrates and calls every pair of empty boxes an overlap.
    Ship ferry = game::buildFerry();
    ferry.initialise(0.0);
    expectTrue("the ferry as shipped validates clean", ferry.validate().empty());

    int erS = kNoCompartment;
    expectTrue("the fixture resolves the space the cable transit leaves",
               ferry.findCompartment("engine_room_s", erS));

    // `cable_transit` as the ferry authors it, with its aft end unresolved: 0.04 m2
    // at z = 2.0 m, which is 3.5 m under the design waterline.
    Opening transit;
    transit.name = "cable_transit";
    transit.a = erS;
    transit.pos = {-8, -4.0, 2.0};
    transit.area = 0.04;
    transit.dischargeCoeff = 0.60;
    transit.kind = OpeningKind::Pipe;
    transit.open = true;

    Ship resolvedToSea = ferry;
    transit.b = kSea;
    resolvedToSea.openings.push_back(transit);
    expectTrue("a miss answered with the sea is a legal opening, and validates clean",
               resolvedToSea.validate().empty());

    Ship refused = ferry;
    transit.b = kNoCompartment;
    refused.openings.push_back(transit);
    const std::vector<std::string> problems = refused.validate();
    expectEqual("a miss answered with kNoCompartment is one problem, not none",
                static_cast<long long>(problems.size()), 1);
    expectTrue("... naming the opening that carries it",
               !problems.empty() && problems[0].find("cable_transit") != std::string::npos);

    // A pump has no sea to stand in for a miss -- its index is checked against the
    // range and nothing else -- so that half was never exposed. It is here to say
    // the openings half now agrees with the half that was already right.
    Ship pumped = ferry;
    pumped.pumps.push_back({"bilge_er_s", kNoCompartment, 0.060, 25.0, false, 0.0});
    expectEqual("and a pump drawing on a refused name is one problem too",
                static_cast<long long>(pumped.validate().size()), 1);
}

// **The refusal must not cost the ferry an opening.**
//
// A builder that refused to load, or that quietly dropped an endpoint it could not
// resolve, would be the original defect wearing a different face. The counts are
// the cheap half; the load-bearing half is that every endpoint in the network is
// either the sea or a compartment that exists, and none is the sentinel a miss
// would leave behind.
void testTheFerryBuildsItsWholeFlowNetworkFromNamesThatResolve() {
    Ship ferry = game::buildFerry();
    ferry.initialise(0.0);
    expectEqual("the ferry builds all 27 openings", static_cast<long long>(ferry.openings.size()),
                27);
    expectEqual("and all 4 pumps", static_cast<long long>(ferry.pumps.size()), 4);

    const int n = static_cast<int>(ferry.compartments.size());
    bool everyEndpointResolves = true;
    int seaEndpoints = 0, compartmentEndpoints = 0;
    for (const Opening& o : ferry.openings)
        for (const int e : {o.a, o.b}) {
            if (e == kSea) {
                ++seaEndpoints;
                continue;
            }
            everyEndpointResolves = everyEndpointResolves && e >= 0 && e < n;
            ++compartmentEndpoints;
        }
    for (const Pump& p : ferry.pumps)
        everyEndpointResolves = everyEndpointResolves && p.compartment >= 0 && p.compartment < n;
    expectTrue("every endpoint is the sea or a compartment that exists", everyEndpointResolves);

    // Both populations counted apart and pinned exactly, not merely as non-zero: a
    // network of pure sea endpoints would satisfy `everyEndpointResolves` on its
    // own, and so would one with none. The breach, three downflooding openings,
    // eight vents, eight air pipes and two counterflooding valves put 22 of the 54
    // endpoints in the sea; the other 32 are the ones a lookup has to get right.
    expectEqual("22 endpoints are the sea", static_cast<long long>(seaEndpoints), 22);
    expectEqual("and the other 32 came from a name", static_cast<long long>(compartmentEndpoints),
                32);
    expectEqual("which is every endpoint there is",
                static_cast<long long>(seaEndpoints + compartmentEndpoints),
                2 * static_cast<long long>(ferry.openings.size()));
    expectTrue("and the definition it builds validates clean", ferry.validate().empty());
}

// --- Endpoints the solver may not read ---------------------------------------

namespace {

// The whole compartment state after a fixed run, not a total. The claim under
// test is bit-identity, and a total hides a compensating pair.
struct FerryRun {
    std::vector<double> water, air, temperature;
    Opening opening;   // the one added below, carrying its last tick's diagnostics
    std::size_t problems = 0;
};

// The ferry as shipped, plus `cable_transit` -- the same 0.04 m2 pipe 3.5 m under
// the design waterline that the lookup tests use -- with its aft end set to
// whatever is being asked about, stepped for `ticks`.
FerryRun runFerryWithTransit(int endpointB, bool open, int ticks) {
    Ship s = game::buildFerry();
    s.initialise(0.0);
    int erS = kNoCompartment;
    expectTrue("the fixture resolves the space the cable transit leaves",
               s.findCompartment("engine_room_s", erS));

    Opening transit;
    transit.name = "cable_transit";
    transit.a = erS;
    transit.b = endpointB;
    transit.pos = {-8, -4.0, 2.0};
    transit.area = 0.04;
    transit.dischargeCoeff = 0.60;
    transit.kind = OpeningKind::Pipe;
    transit.open = open;
    s.openings.push_back(transit);

    FerryRun r;
    r.problems = s.validate().size();
    for (int i = 0; i < ticks; ++i) s.step(0.01, Sea(0.0));
    for (const Compartment& c : s.compartments) {
        r.water.push_back(c.waterVolume);
        r.air.push_back(c.airMass);
        r.temperature.push_back(c.gasTemperature);
    }
    r.opening = s.openings.back();
    return r;
}

// Exact equality, deliberately. `==` on a NaN is false, where `expectNear` with a
// zero tolerance would *pass* one -- `std::abs(NaN - want) > 0.0` is false. A run
// that read garbage doubles out of bounds is exactly the run that would arrive
// carrying NaNs, so the comparison has to be the one that catches them.
bool sameState(const FerryRun& a, const FerryRun& b) {
    if (a.water.size() != b.water.size()) return false;
    for (std::size_t i = 0; i < a.water.size(); ++i)
        if (a.water[i] != b.water[i] || a.air[i] != b.air[i] ||
            a.temperature[i] != b.temperature[i])
            return false;
    return true;
}

bool passedNothing(const Opening& o) {
    return o.lastFlow == 0.0 && o.lastGasMassFlow == 0.0 && o.lastGasEnthalpyFlow == 0.0 &&
           o.lastExchangeDown == 0.0 && o.lastExchangeUp == 0.0 && o.lastExchangeMassDown == 0.0 &&
           o.lastExchangeMassUp == 0.0;
}

}  // namespace

// **An endpoint that names no space is refused by the solver, not read.**
//
// `Ship::solveFlowNetwork` indexes `compartments[i]` for every endpoint that is
// not `kSea`, and kSea is the *only* negative it exempts. `kNoCompartment` (-3)
// and breach.hpp's `kEnclosedVoid` (-2) therefore went straight through as
// indices: `compartments[-3]` is an out-of-bounds read, and `dWater[-3] -= dv`
// sixty lines later an out-of-bounds *write* into a heap vector. Against the code
// before kNoCompartment existed, renaming one ferry compartment reached "double
// free or corruption (out)" and exit 134.
//
// `Ship::validate()` names such an endpoint and that was the entire defence.
// It is advisory -- nothing on the step path consults it -- and the network is
// not frozen when it runs: `applyBreaches` pushes openings into a ship that is
// already stepping, so the endpoints read at tick N are not the set validate()
// saw. Both halves are asserted here: validate() still names it, *and* the solver
// refuses it on its own.
//
// The assertion is not "it did not crash", which is what undefined behaviour is
// entitled to look like on a good day. It is that a refused endpoint is exactly a
// shut opening -- same floodwater, same air mass, same gas temperature in every
// compartment, bit for bit, after 400 ticks of a ferry that is flooding hard
// through her breach the whole time.
void testTheFlowSolverRefusesAnEndpointThatNamesNoSpace() {
    constexpr int kTicks = 400;
    const FerryRun shut = runFerryWithTransit(kSea, false, kTicks);
    expectEqual("the ferry plus a shut transit validates clean",
                static_cast<long long>(shut.problems), 0);

    // Against vacuity, and it is the load-bearing control: the identity below
    // would hold trivially if this opening could never move anything. The same
    // pipe with a *legal* endpoint floods the engine room it opens, and does so by
    // enough to be visible against a ship already flooding through her breach.
    const FerryRun sea = runFerryWithTransit(kSea, true, kTicks);
    expectEqual("the same transit opened to the sea also validates clean",
                static_cast<long long>(sea.problems), 0);
    expectTrue("an open cable transit is passing water on the last tick",
               sea.opening.lastFlow != 0.0);
    expectTrue("and it has changed the ship, so the comparison below has something to catch",
               !sameState(sea, shut));

    // Every way an index can fail to name a space. The two sentinels are the ones
    // this is about; the three out-of-range integers are here because the defect
    // was a missing *range* check and a check that only tested for -3 and -2 would
    // be a check written to the two bugs already found.
    const int n = static_cast<int>(shut.water.size());
    const struct { const char* what; int endpoint; } refusals[] = {
        {"kNoCompartment, the answer a name lookup refuses with", kNoCompartment},
        {"kEnclosedVoid, inside the hull and inside no compartment", -2},
        {"one past the last compartment", n},
        {"far past it", n + 97},
        {"an uninitialised-looking negative", -12345},
    };

    for (const auto& r : refusals) {
        const FerryRun run = runFerryWithTransit(r.endpoint, true, kTicks);
        expectEqual(std::string("validate() names an endpoint at ") + r.what,
                    static_cast<long long>(run.problems), 1);
        expectTrue(std::string("... and the solver passes nothing through it: ") + r.what,
                   passedNothing(run.opening));
        expectTrue(std::string("... leaving the ship bit-identical to the shut opening: ") + r.what,
                   sameState(run, shut));
    }
}

// **The pump half of the same loop, which indexes on an author-supplied int with
// no sentinel to exempt at all.**
//
// A pump has no sea to be mistaken for a space, so a missed name lookup could
// never reach it -- but `compartments[p.compartment]` and `dWater[p.compartment]`
// are the identical unchecked reads, and validate()'s range test on pumps is
// advisory in exactly the same way. Asserted alongside the openings half so the
// two cannot drift apart.
void testARunningPumpOnAnEndpointThatNamesNoSpaceMovesNothing() {
    constexpr int kTicks = 200;

    auto runWith = [&](int compartment, bool on) {
        Ship s = game::buildFerry();
        s.initialise(0.0);
        s.pumps.push_back({"bilge_nowhere", compartment, 0.060, 25.0, on, 0.0});
        const std::size_t problems = s.validate().size();
        for (int i = 0; i < kTicks; ++i) s.step(0.01, Sea(0.0));
        std::vector<double> water;
        for (const Compartment& c : s.compartments) water.push_back(c.waterVolume);
        return std::tuple{water, s.pumps.back().lastFlow, problems};
    };

    int erS = kNoCompartment;
    expectTrue("the fixture resolves a space a bilge pump could draw on",
               game::buildFerry().findCompartment("engine_room_s", erS));

    const auto [offWater, offFlow, offProblems] = runWith(erS, false);
    expectEqual("a legal pump left off is no problem", static_cast<long long>(offProblems), 0);

    // The control: the same pump, running, on a compartment that exists. It has to
    // move water or the nulls below say nothing.
    const auto [onWater, onFlow, onProblems] = runWith(erS, true);
    expectEqual("nor is one that is running", static_cast<long long>(onProblems), 0);
    expectTrue("a running bilge pump on a flooding engine room discharges", onFlow > 0.0);
    expectTrue("and it leaves less water behind than the same pump switched off",
               onWater != offWater);

    for (const int bad : {kNoCompartment, -2, -12345, static_cast<int>(offWater.size()) + 3}) {
        const auto [water, flow, problems] = runWith(bad, true);
        expectEqual("validate() names a pump on a compartment that does not exist",
                    static_cast<long long>(problems), 1);
        expectTrue("... and it discharges exactly nothing", flow == 0.0);
        expectTrue("... leaving the ship where the pump switched off leaves her",
                   water == offWater);
    }
}

void runShipTests() {
    std::printf("\n--- horizontal openings and buoyant exchange ---\n");
    testHorizontalExchangeMatchesTheCounterflowClosedForm();
    testTheExchangeIsExactlyZeroWhereverItShouldNotRun();
    testAtRestTheTwoStreamsAreExactlyEqualVolumes();
    testTheExchangeMeetsTheSingleOrificeSolveAtTheEdgeOfItsWindow();
    testTrappedAirArrestsABreachButNotAHatch();
    testTheExchangeGoesAsTheCosineOfTheDeckSTilt();
    testMassIsConservedAcrossTheExchange();
    testTwoGasStreamsBalanceTheirOwnEnthalpy();
    testTheFerryAsShippedRunsNoExchangeAtAll();
    testAnOpenEscapeTrunkDrainsTheVehicleDeckIntoADryEngineRoom();
    testAMasslessShipStepsToFiniteStateRatherThanNaN();

    std::printf("\n--- compartment lookup: a miss is not the open sea ---\n");
    testAMissedCompartmentLookupIsRefusedRatherThanAnsweredWithTheSea();
    testAnEndpointFromAMissedLookupIsNamedWhereASeaEndpointCannotBe();
    testTheFerryBuildsItsWholeFlowNetworkFromNamesThatResolve();

    std::printf("\n--- the flow solver refuses an endpoint it may not index ---\n");
    testTheFlowSolverRefusesAnEndpointThatNamesNoSpace();
    testARunningPumpOnAnEndpointThatNamesNoSpaceMovesNothing();
}
