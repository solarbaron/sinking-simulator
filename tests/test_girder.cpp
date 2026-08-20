// SPDX-License-Identifier: MIT
//
// Validation of the hull girder -- the ship as a beam.
//
// The file is split the way the code is, and for the same reason. The beam
// integration has exact answers and is asserted exactly: a uniform beam carries
// no load, a concentrated weight amidships gives a midship moment of exactly
// W*L/8, and a free-floating body must close to zero shear and zero moment at
// both perpendiculars. The *distributions* that feed it are assumptions -- a
// trapezoidal lightship curve is a construction, not a measurement -- so those
// are checked for the properties they must have (they integrate to the right
// total and the right centre) rather than against a value nobody can state.
#include "engine/sim/girder.hpp"
#include "engine/sim/scantlings.hpp"
#include "harness.hpp"

#include <cmath>
#include <cstdio>
#include <vector>

using namespace sim;
using testing::expectEqual;
using testing::expectNear;
using testing::expectTrue;

namespace {

std::vector<double> uniformStations(double lo, double hi, int n) {
    std::vector<double> x;
    for (int i = 0; i < n; ++i) x.push_back(lo + (hi - lo) * i / (n - 1));
    return x;
}

// --- The beam, against algebra -----------------------------------------------

// Weight exactly balancing buoyancy everywhere carries no shear and no moment.
// Trivial to state and the first thing a sign error breaks.
void testBalancedBeamCarriesNothing() {
    const std::vector<double> x = uniformStations(-60, 60, 41);
    const std::vector<double> flat(x.size(), 1.0e5);
    const HullGirder g = integrateGirder(x, flat, flat);

    double worstShear = 0, worstMoment = 0;
    for (const GirderStation& s : g.stations) {
        worstShear = std::max(worstShear, std::abs(s.shear));
        worstMoment = std::max(worstMoment, std::abs(s.moment));
    }
    expectTrue("a beam in balance carries no shear", worstShear < 1e-9);
    expectTrue("and no bending moment", worstMoment < 1e-9);
    expectNear("weight and buoyancy integrate to the same total", g.totalWeight, g.totalBuoyancy,
               1e-6);
}

// The textbook case, and the reason this file exists. A beam of length L with a
// weight W concentrated amidships, floating on uniform buoyancy W/L, carries a
// midship bending moment of exactly W*L/8. Every step of that is hand-checkable:
// the shear is +-W/2 at the quarter points and the moment is a parabola.
//
// Concentrated weight amidships *sags*, so the sign must be negative.
void testConcentratedWeightAmidshipsGivesWLOverEight() {
    const double length = 120.0, weight = 8.0e6;   // N
    const int n = 401;                             // odd, so a station lands on midship
    const std::vector<double> x = uniformStations(-0.5 * length, 0.5 * length, n);

    // A "point" weight is one station's worth: W spread over the width that
    // station owns, so the total is exactly W however fine the grid.
    const double spacing = length / (n - 1);
    std::vector<double> w(x.size(), 0.0);
    w[(n - 1) / 2] = weight / spacing;
    const std::vector<double> b(x.size(), weight / length);

    const HullGirder g = integrateGirder(x, w, b);

    expectNear("the concentrated weight integrates to W", g.totalWeight, weight, 1e-6 * weight);
    expectNear("buoyancy integrates to the same", g.totalBuoyancy, weight, 1e-9 * weight);

    const GirderStation& mid = g.stations[(n - 1) / 2];
    expectNear("midship bending moment is W L / 8", mid.moment, -weight * length / 8.0,
               0.01 * weight * length / 8.0);
    expectTrue("a weight amidships sags rather than hogs", g.maxMoment < 0);

    // Shear peaks at +-W/2 immediately either side of the weight.
    expectNear("shear reaches -W/2 just aft of the weight", g.stations[(n - 1) / 2 - 1].shear,
               -0.5 * weight, 0.02 * weight);
    expectNear("and +W/2 just forward of it", g.stations[(n - 1) / 2 + 1].shear, 0.5 * weight,
               0.02 * weight);

    // Free ends. This is the check that makes the whole calculation trustworthy.
    expectTrue("shear closes at the forward perpendicular", std::abs(g.shearClosure) < 1e-9);
    expectTrue("bending moment closes at the forward perpendicular",
               std::abs(g.momentClosure) < 1e-9);
}

// The mirror image: buoyancy concentrated amidships hogs. Same magnitude,
// opposite sign. Without this, a sign error that flipped *both* would pass the
// test above by matching a negated closed form.
void testConcentratedBuoyancyAmidshipsHogs() {
    const double length = 120.0, force = 8.0e6;
    const int n = 401;
    const std::vector<double> x = uniformStations(-0.5 * length, 0.5 * length, n);
    const double spacing = length / (n - 1);

    std::vector<double> b(x.size(), 0.0);
    b[(n - 1) / 2] = force / spacing;
    const std::vector<double> w(x.size(), force / length);

    const HullGirder g = integrateGirder(x, w, b);
    expectNear("midship moment is the mirror of the sagging case",
               g.stations[(n - 1) / 2].moment, force * length / 8.0,
               0.01 * force * length / 8.0);
    expectTrue("buoyancy amidships hogs", g.hogging());
    expectTrue("and closes at the far end", std::abs(g.momentClosure) < 1e-9);
}

// An out-of-balance ship must not quietly produce a bending moment curve. The
// closure residual is the instrument, and it has to actually respond.
void testImbalanceShowsInTheClosure() {
    const std::vector<double> x = uniformStations(-60, 60, 41);
    const std::vector<double> w(x.size(), 1.0e5);
    std::vector<double> b(x.size(), 1.0e5);
    for (double& v : b) v *= 1.10;   // 10% too buoyant: she is not floating here

    const HullGirder g = integrateGirder(x, w, b);
    // Against the W L / 8 reference the residual works out at 0.4 for a 10%
    // imbalance -- large next to the ~0 a balanced beam gives, and the threshold
    // is set from that arithmetic rather than from a round number.
    expectTrue("a 10% imbalance leaves a large residual at the far end",
               std::abs(g.momentClosure) > 0.2);
    // A 10% imbalance trips three guards at once -- shear closure, moment closure
    // and the buoyancy/weight total -- so "says so" does not say which, and any one
    // of the three could be deleted with the suite still green. They are three
    // different first moves for whoever reads the message: an integration bug at
    // the ends, a moment-arm bug, or a ship that is simply not floating here.
    const std::vector<std::string> said = validateGirder(g);
    expectTrue("and validate() says so", !said.empty());
    const auto mentions = [&](const char* s) {
        for (const std::string& p : said)
            if (p.find(s) != std::string::npos) return true;
        return false;
    };
    expectTrue("naming the imbalance, which is what the fixture built",
               mentions("not floating at this attitude"));
    expectTrue("and the moment residual it leaves at the far end",
               mentions("bending moment does not close"));

    // And the balanced case must be clean, or the check above is meaningless.
    const HullGirder ok = integrateGirder(x, w, std::vector<double>(x.size(), 1.0e5));
    expectTrue("a balanced beam raises nothing", validateGirder(ok).empty());
}

// --- The distributions, against their defining properties --------------------

// A 120 m x 18 m box barge floating at 5 m, scaled by `k` in every length. The
// scale factor exists for testBalanceConvergesTheSameWayAtEverySize() below:
// geometrically similar hulls are the same hydrostatic problem written in
// different units, so anything that behaves differently across `k` is carrying a
// dimensional mistake.
Ship scaledBarge(double lcgOffset, double k) {
    Ship s;
    std::vector<Station> stations;
    for (int i = 0; i <= 40; ++i) {
        Station q;
        q.x = k * (-60.0 + 120.0 * i / 40.0);
        q.halfBeam = {9.0 * k, 9.0 * k};
        stations.push_back(q);
    }
    s.hull = makeHullFromStations(stations, {0.0, 12.0 * k});
    s.deckEdgeZ = 12.0 * k;
    s.lightshipMass = k * k * k * 120.0 * 18.0 * 5.0 * kRhoSeawater;   // floats at 5k m
    s.lightshipCog = {lcgOffset * k, 0.0, 6.0 * k};
    s.gyradii = {6.0 * k, 30.0 * k, 30.0 * k};
    return s;
}

Ship barge(double lcgOffset) { return scaledBarge(lcgOffset, 1.0); }

// The trapezoidal lightship curve is an assumption, but it is an assumption with
// two exact properties: it integrates to the lightship weight and its centroid is
// the lightship LCG. Those are what it was constructed to satisfy, so they are
// what to check.
void testWeightDistributionCarriesItsOwnTotalAndCentre() {
    for (double lcg : {-6.0, 0.0, 4.5}) {
        Ship ship = barge(lcg);
        ship.initialise(0.0);
        const std::vector<double> x = girderStations(ship, 81);
        const std::vector<double> w = weightDistribution(ship, x);

        // Simpson, not trapezoidal. The distribution is linear so its integral is
        // exact either way, but the first moment integrates w*x, which is
        // quadratic -- and trapezoidal is not exact for a quadratic. The first
        // version of this check used it and reported the construction as 0.03%
        // wrong when the construction is exact and the *check* was not.
        double total = 0, moment = 0;
        for (std::size_t i = 2; i < x.size(); i += 2) {
            const double h = (x[i] - x[i - 2]) / 6.0;
            total += h * (w[i - 2] + 4.0 * w[i - 1] + w[i]);
            moment += h * (w[i - 2] * x[i - 2] + 4.0 * w[i - 1] * x[i - 1] + w[i] * x[i]);
        }
        const double want = ship.lightshipMass * kGravity;
        expectNear("the weight curve integrates to the lightship weight", total, want,
                   1e-9 * want);
        expectNear("and its centroid is the lightship LCG", moment / total, lcg, 1e-9);
    }
}

// Buoyancy comes from slab clipping, so its integral must be the displacement the
// hydrostatics already report -- computed a completely different way.
void testBuoyancyDistributionIntegratesToDisplacement() {
    Ship ship = barge(0.0);
    ship.initialise(0.0);
    const Sea still(0.0);

    const std::vector<double> x = girderStations(ship, 81);
    const std::vector<double> b = buoyancyDistribution(ship, still, x);
    double total = 0;
    for (std::size_t i = 1; i < x.size(); ++i)
        total += 0.5 * (x[i] - x[i - 1]) * (b[i] + b[i - 1]);

    const Diagnostics d = ship.diagnostics(still);
    const double buoyancy = ship.seaDensity * kGravity * d.buoyantVolume;
    expectTrue("the ship is actually floating", buoyancy > 1e6);
    expectNear("the buoyancy curve integrates to the displacement", total, buoyancy,
               0.02 * buoyancy);
}

// A prismatic barge, evenly loaded, is the one real ship with an exactly known
// answer: uniform weight against uniform buoyancy, so no bending anywhere. Any
// residual is the distribution machinery's own error, which makes this the
// calibration for everything that follows.
void testAnEvenlyLoadedBargeBarelyBends() {
    Ship ship = barge(0.0);
    ship.initialise(0.0);
    const HullGirder g = hullGirder(ship, Sea(0.0), 81);

    expectTrue("the calculation produced stations", g.stations.size() == 81);
    expectTrue("weight and buoyancy balance", std::abs(g.totalBuoyancy - g.totalWeight) <
                                                  0.02 * g.totalWeight);

    // Scale the residual moment against something meaningful: W*L/8 is what a
    // fully concentrated load would give, so this says "a fraction of a per cent
    // of the worst case".
    const double reference = g.totalWeight * (ship.hullHi.x - ship.hullLo.x) / 8.0;
    expectTrue("an evenly loaded prismatic barge carries almost no bending moment",
               std::abs(g.maxMoment) < 0.02 * reference);
    expectTrue("and the calculation closes", std::abs(g.momentClosure) < 0.05);
    expectTrue("with nothing to report", validateGirder(g).empty());
}

// The case the whole tier exists for. A wave crest amidships lifts the middle and
// the ship hogs; a trough drops it and she sags. Same ship, same water, opposite
// signs -- and neither is anything the code was told.
void testCrestHogsAndTroughSags() {
    Ship ship = barge(0.0);
    ship.initialise(0.0);

    // A wave twice the ship's length, so one crest spans her.
    const double length = ship.hullHi.x - ship.hullLo.x;
    const double omega = std::sqrt(kGravity * 2.0 * kPi / (2.0 * length));
    const WaveField field = WaveField::regular(2.5, omega, 0.0);

    Sea crest, trough;
    crest.waves = &field;
    trough.waves = &field;
    // Scan a full period for the moments when the surface amidships is highest
    // and lowest, rather than assuming where the crest is.
    const double period = 2.0 * kPi / omega;
    double best = -1e30, worst = 1e30;
    for (int i = 0; i < 720; ++i) {
        const double t = period * i / 720.0;
        const double eta = field.elevation(ship.state.position.x, 0.0, t);
        if (eta > best) { best = eta; crest.time = t; }
        if (eta < worst) { worst = eta; trough.time = t; }
    }
    expectTrue("the scan found a real crest and trough", best > 1.0 && worst < -1.0);

    const HullGirder onCrest = hullGirder(ship, crest, 81);
    const HullGirder inTrough = hullGirder(ship, trough, 81);

    expectTrue("a crest amidships hogs", onCrest.maxMoment > 0);
    expectTrue("a trough amidships sags", inTrough.maxMoment < 0);

    // Vacuity guard: both must be substantial, or "opposite signs" is a statement
    // about rounding. Compare against the still-water case, which is near zero.
    const double reference = onCrest.totalWeight * length / 8.0;
    expectTrue("and both are a real bending moment, not noise",
               onCrest.maxMoment > 0.02 * reference && inTrough.maxMoment < -0.02 * reference);

    // Both peaks belong near amidships, not at an end.
    expectTrue("the hogging peak is amidships", std::abs(onCrest.maxMomentX) < 0.25 * length);
    expectTrue("the sagging peak is amidships", std::abs(inTrough.maxMomentX) < 0.25 * length);
}

// Flooding a compartment is a concentrated weight, and a concentrated weight is
// what bends a hull girder. This is the coupling that makes the tier worth
// having: the flooding solver already knows where the water is.
void testFloodingBendsTheHull() {
    Ship dry = barge(0.0);
    dry.initialise(0.0);
    const HullGirder before = hullGirder(dry, Sea(0.0), 81);

    Ship wet = barge(0.0);
    Compartment c;
    c.name = "hold";
    c.mesh = clipToBox(wet.hull, {-15, -9, 0}, {15, 9, 6});
    c.permeability = 0.95;
    wet.compartments.push_back(c);
    wet.initialise(0.0);   // computes gross volumes, so floodableVolume() is real
    wet.compartments[0].waterVolume = 0.6 * wet.compartments[0].floodableVolume();
    wet.initialise(0.0);   // and again, to re-level the water and re-find the draft
    const HullGirder after = hullGirder(wet, Sea(0.0), 81);

    // **And carries exactly the water that is in it.** `weightDistribution` spreads
    // each compartment's water over its own extent, which the header states as the
    // contract, and `integrateGirder` reads the result trapezoidally -- so what the
    // girder ends up carrying is the sum of `perLength` times each station's own
    // slab, not `perLength` times the span. Asserted as a fraction because the
    // failure it catches is a fraction: one station spacing out of the compartment's
    // length.
    const double water = wet.compartments[0].waterVolume * wet.seaDensity * kGravity;
    const double delivered = after.totalWeight - before.totalWeight;
    std::printf("     floodwater: %.4e N in the hold, %.4e N delivered to the girder (%.2f%%)\n",
                water, delivered, 100.0 * delivered / water);
    expectNear("the girder carries exactly the water in the hold", delivered, water,
               1e-6 * water);

    expectTrue("the flooded case carries more weight", after.totalWeight > before.totalWeight);
    expectTrue("water amidships sags the hull", after.maxMoment < before.maxMoment);
    const double reference = after.totalWeight * (wet.hullHi.x - wet.hullLo.x) / 8.0;
    expectTrue("and by an amount worth calling a bending moment",
               std::abs(after.maxMoment) > 0.01 * reference);
}

// --- Balancing her on the wave -----------------------------------------------

// The balance drives two residuals to zero, and they are not the same kind of
// quantity: one is a net force in N, the other a net trimming moment in N m. A
// tolerance has the dimensions of the thing it bounds, so those two cannot share
// one -- and this loop gave them one, testing the moment against `1e-6 * weight`.
// Read as a lever that is asking for the centre of buoyancy within a fixed
// *micrometre* of the centre of gravity, on a ship of any size at all.
//
// The signature of a tolerance whose dimensions are wrong is that it is not
// scale free, and that is what this checks rather than a step count on one hull.
// The same barge at four sizes spanning 12 m to 1200 m is one hydrostatic problem
// written in four sets of units, so a criterion with the right dimensions has to
// stop it after the same number of Newton steps. Measured: nine steps at every
// size against `weight * length`, and 11, 13, 14 and 15 -- climbing with the
// hull, because a micrometre is a finer demand on a longer ship -- against
// `weight` alone.
//
// It is worth a real number: each step is three whole-hull clip-and-integrate
// passes, and `hullGirder()` takes this path on every Tier-0 review.
void testBalanceConvergesTheSameWayAtEverySize() {
    int first = -1;
    for (double k : {0.1, 1.0, 3.0, 10.0}) {
        Ship ship = scaledBarge(4.5, k);
        ship.initialise(0.0);
        const Sea still(0.0);

        int used = -1;
        expectTrue("the barge balances at every size", balanceOnWave(ship, still, 40, &used));
        // Vacuity guard: a ship that was already poised would report no steps at
        // every size, and the comparison below would be between four zeroes. The
        // +4.5 m LCG is there to make her trim.
        expectTrue("and balancing her is a real iteration, not a no-op", used >= 5);
        if (first < 0) first = used;
        expectEqual("the same problem takes the same number of steps at every size", used, first);
        expectTrue("which is the nine measured, not the forty-step budget", used <= 10);

        // What that criterion promises, read back from the balanced ship through
        // diagnostics() rather than from the integral the solver used: buoyancy
        // within 1e-6 of weight. Asserted at the tolerance itself and not at the
        // 7.5e-9 measured, because the loop can exit the moment *both* legs are
        // inside and it is the moment leg that binds here -- how far past its own
        // bound the force leg happens to be is luck, and 1e-6 is the promise.
        const Diagnostics d = ship.diagnostics(still);
        const double weight = d.displacementMass * kGravity;
        const double buoyancy = ship.seaDensity * d.buoyantVolume * kGravity;
        expectTrue("she is floating at all", buoyancy > 0);
        expectNear("buoyancy equals weight to the tolerance the loop stops on", buoyancy, weight,
                   1e-6 * weight);
    }
}

// The other half of scaling the moment tolerance by a length: there has to *be*
// a length. `hullLo`/`hullHi` are cached by Ship::initialise(), and on a ship
// that has never been initialised the tolerance would come out zero -- a
// criterion `std::abs(m0) < 0` that no residual can ever be under, so the loop
// would run its whole budget on any ship it did not otherwise give up on. She is
// refused instead, which is what girderStations() already does with the same ship.
void testAShipWithNoHullExtentIsRefused() {
    Ship raw = barge(0.0);   // deliberately not initialise()d
    expectTrue("the fixture really has no cached hull extent", !(raw.hullHi.x > raw.hullLo.x));

    int used = -1;
    expectTrue("a ship with no length to measure a moment against is refused",
               !balanceOnWave(raw, Sea(0.0), 40, &used));
    // And refused before spending anything, which is what says the refusal came
    // from the missing length rather than from the iteration failing later.
    expectEqual("and refused before taking a single Newton step", used, 0);
    expectTrue("the same ship girderStations() has always refused",
               girderStations(raw, 41).empty());
}

// The 2 x 2 Newton solve, and specifically its refusal to divide by a singular
// Jacobian. That refusal is what the balance falls back on when a hull leaves the
// water, and it can only work if it is measured against the size of the terms
// that cancelled: `a` is a waterplane stiffness in N/m and `d` a trim stiffness
// in N m/rad, so on a real hull their product is around 1e17 and *nothing* the
// problem produces comes near an absolute floor of 1e-30.
void testASingularBalanceJacobianIsRefusedOnItsOwnScale() {
    // A ferry's own numbers: 1.9e7 N/m of waterplane stiffness against 2.2e10
    // N m/rad of trim stiffness, coupled by the offset of the centre of flotation.
    const double a = 1.9e7, d = 2.2e10, b = -3.0e8, c = -3.0e8;
    const double force = 1.0e6, moment = 4.0e7;   // N and N m

    double dz = 0, dTrim = 0;
    expectTrue("a well conditioned Jacobian solves",
               solveBalanceStep(a, b, c, d, force, moment, dz, dTrim));
    // Checked by substitution rather than against a second Cramer's rule, which
    // would agree with a transposed one.
    expectNear("and the step satisfies the force row", a * dz + b * dTrim, force,
               1e-9 * std::abs(force));
    expectNear("and the moment row", c * dz + d * dTrim, moment, 1e-9 * std::abs(moment));

    // Now the same matrix made singular, at exactly the same scale: `b * c` is
    // built to sit 1e4 away from `a * d`, so every term is a ship's and the
    // determinant is nothing but what failed to cancel.
    const double bSing = 2.0e10, cSing = (a * d - 1.0e4) / bSing;
    const double det = a * d - bSing * cSing;
    expectNear("the singular case is built where it was meant to be", det, 1.0e4, 1.0e3);
    // The two halves of the point, and the first is why the second is needed:
    // an absolute floor calls a determinant of ten thousand perfectly healthy.
    expectTrue("an absolute floor of 1e-30 waves that determinant through",
               std::abs(det) > 1e-30);
    expectTrue("though it is a part in 4e13 of the products that made it",
               std::abs(det) < 1e-12 * (std::abs(a * d) + std::abs(bSing * cSing)));
    dz = dTrim = 0;
    expectTrue("so a Jacobian singular on its own scale is refused",
               !solveBalanceStep(a, bSing, cSing, d, force, moment, dz, dTrim));

    // And the degenerate end, which the relative form still has to catch: nothing
    // of the hull in the water, so every term is exactly zero and the scale it
    // would be compared against is zero too.
    expectTrue("a hull out of the water entirely is refused as well",
               !solveBalanceStep(0, 0, 0, 0, -1.0e8, 0, dz, dTrim));
}

// The trim residual is small **about the right point**, which is not the same
// claim as the residual being small.
//
// The loop solves for the attitude, and the centre of gravity's world position is
// a function of the attitude: the body-frame cog is fixed, but the ship rotates
// underneath it. Reading that world position once, at the orientation the ship
// entered with, and then driving the centre of buoyancy onto it is a balance
// against a stale target, and the staleness is `xg (cos t - 1) + zg sin t` -- to
// first order `KG * trim`, and independent of how well the loop converged.
//
// That is exactly the failure a convergence test cannot see, so this asserts the
// geometry instead: after balancing, the centre of buoyancy is under the centre
// of gravity **as she is now floating**. Measured against the frozen point the
// old code aimed at, the barge below stopped 0.11 m away -- 9.2e-4 of her length,
// 900x the tolerance the loop stops on, and `6.0 * 0.018554 = 0.1113` to four
// figures, which is what names the mechanism rather than merely sizing the error.
//
// Swept over size for the same reason the step-count test above is: geometric
// similarity makes the offset a fixed fraction of the length, so a defect that
// carries a stray absolute term shows up as a fraction that moves with `k`.
void testTheBalanceFollowsTheCentreOfGravityAsSheTrims() {
    double firstFraction = 0.0;
    bool haveFirst = false;
    for (double k : {0.1, 1.0, 3.0, 10.0}) {
        Ship entry = scaledBarge(4.5, k);
        entry.initialise(0.0);
        const Sea still(0.0);
        const double length = entry.hullHi.x - entry.hullLo.x;

        // The frozen point, computed exactly as the balance used to compute it:
        // the world x of the centre of gravity at the *entry* orientation.
        // `massProperties()` is orientation-free, so the body-frame cog this is
        // built from is the same one the balanced ship reports below.
        const Diagnostics before = entry.diagnostics(still);
        const Mat3 R0 = entry.state.orientation.toMat3();
        const double frozen = (R0 * before.centreOfGravity + entry.state.position).x;
        // Precondition for the closed form at the end: she enters upright, so the
        // frozen reading is the body-frame LCG and the drift below is the whole of
        // what trimming does to it.
        expectNear("the fixture enters upright", (R0 * Vec3{1, 0, 0}).z, 0.0, 1e-15);

        Ship ship = entry;
        expectTrue("the barge balances", balanceOnWave(ship, still));

        const Diagnostics d = ship.diagnostics(still);
        const Mat3 R = ship.state.orientation.toMat3();
        const double lcgLive = (R * d.centreOfGravity + ship.state.position).x;
        const double lcb = (R * d.centreOfBuoyancy + ship.state.position).x;
        const double trim = std::asin(-(R * Vec3{1, 0, 0}).z);

        // Vacuity guard, and the whole reason the fixture carries a +4.5 m LCG: on
        // a ship that does not trim the live point and the frozen one are the same
        // point, and every assertion below is satisfied by arithmetic.
        expectTrue("the offset LCG really trims her", std::abs(trim) > 0.015);

        // What the balance promises. 1e-6 of a length is the loop's own moment
        // tolerance read back as a distance, and it is asserted at the tolerance
        // rather than at the 7.4e-7 measured because the loop stops the moment
        // *both* legs are inside -- how far past its own bound the binding leg
        // happens to land is luck.
        expectNear("the centre of buoyancy sits under the live centre of gravity", lcb, lcgLive,
                   1e-6 * length);

        // And what aiming at the frozen point cost, closed form rather than
        // remembered.
        const double drift = lcgLive - frozen;
        expectNear("the frozen point is KG*trim away from the live one, exactly", drift,
                   d.centreOfGravity.x * (std::cos(trim) - 1.0) +
                       d.centreOfGravity.z * std::sin(trim),
                   1e-12 * length);
        expectTrue("which is hundreds of times the tolerance the loop converges to",
                   std::abs(drift) > 500.0 * 1e-6 * length);

        // Geometric similarity: one problem in four sets of units, so the offset
        // is one fraction of the length. A stray absolute term would not be.
        const double fraction = drift / length;
        if (!haveFirst) { firstFraction = fraction; haveFirst = true; }
        expectNear("and it is the same fraction of the length at every size", fraction,
                   firstFraction, 1e-9);
    }
}

// **The success the function reports has to be about both residuals**, because a
// ship converged in heave and still rotating is precisely the state `girder.hpp`
// says destroys the result: the shear and bending curves grow towards the far end
// instead of closing. The moment was computed on the way out and dropped.
//
// Starving the budget is how to build that state deliberately. `initialise()`
// leaves the barge at the draft her weight asks for and no trim, so the force
// residual starts near zero and stays near zero, while one Newton step is nowhere
// near enough for the +4.5 m LCG's 1.9 degrees of trim. Measured: 8.8e-5 of her
// weight out in heave -- eleven times *inside* the force leg -- against 1.13e-2 of
// `weight * length` out in moment, eleven times outside the moment one.
void testABalanceStillRotatingIsNotReportedAsSuccess() {
    const Sea still(0.0);
    Ship starved = barge(4.5);
    starved.initialise(0.0);
    const double length = starved.hullHi.x - starved.hullLo.x;

    int used = -1;
    expectTrue("one Newton step is not a balance, and is not reported as one",
               !balanceOnWave(starved, still, 1, &used));
    expectEqual("and it is the budget that ran out, not the step refusing", used, 1);

    // The leg the old return could see. Asserted explicitly, because it is the
    // reason the old return said yes -- and so the reason reverting the moment leg
    // turns this test red rather than leaving it green for some other cause.
    const Diagnostics d = starved.diagnostics(still);
    const double weight = d.displacementMass * kGravity;
    const double buoyancy = starved.seaDensity * d.buoyantVolume * kGravity;
    expectTrue("she is comfortably inside the force leg all the same",
               std::abs(buoyancy - weight) < 1e-3 * weight);

    // The leg it could not. Taken about the live centre of gravity from
    // `diagnostics()`, which shares no expression with the residual the solver
    // formed, so this is a second route to the same quantity rather than a copy.
    const Mat3 R = starved.state.orientation.toMat3();
    const double lcg = (R * d.centreOfGravity + starved.state.position).x;
    const double lcb = (R * d.centreOfBuoyancy + starved.state.position).x;
    expectTrue("while the trimming moment is an order of magnitude outside its own",
               std::abs(buoyancy * (lcb - lcg)) > 5e-3 * weight * length);

    // And the physical consequence, which is what makes the verdict worth having:
    // the curves the caller would have integrated off this ship.
    const HullGirder ruined = hullGirder(starved, still, 41, false);
    expectTrue("so the bending moment does not close at the forward perpendicular",
               std::abs(ruined.momentClosure) > 0.05);
    const std::vector<std::string> said = validateGirder(ruined);
    bool named = false;
    for (const std::string& p : said)
        if (p.find("bending moment does not close") != std::string::npos) named = true;
    expectTrue("and validateGirder() names that, not something else", named);

    // The control, without which "returns false" is satisfied by a function that
    // always returns false: the same ship, the same code, given the budget.
    Ship poised = barge(4.5);
    poised.initialise(0.0);
    expectTrue("the same ship with the default budget does balance",
               balanceOnWave(poised, still, 40, &used));
    expectTrue("using far less of it than the starved run was given", used > 1 && used < 20);
    const HullGirder good = hullGirder(poised, still, 41, false);
    expectTrue("and her bending moment closes", std::abs(good.momentClosure) < 0.01);
    expectTrue("with nothing for validateGirder() to report", validateGirder(good).empty());
}

// --- From moment to stress ---------------------------------------------------

// sigma = M / Z is arithmetic, so it is asserted as arithmetic. The part worth
// testing is the *sign*: hogging stretches the deck and compresses the keel, and
// a stress magnitude alone cannot tell those apart -- which matters because deck
// plating in compression buckles well below yield.
void testStressFollowsTheMomentAndItsSign() {
    const Scantlings sc = ferryScantlings();
    Ship ferry = barge(0.0);
    ferry.initialise(0.0);
    const StructuralMesh structure = makeStructuralMesh(ferry.hull, sc);
    expectTrue("the structural mesh has panels to work with", !structure.panels.empty());

    const std::vector<double> x = girderStations(ferry, 41);
    const std::size_t mid = x.size() / 2;

    // A pure hogging moment, imposed rather than derived, so the arithmetic has
    // nothing else in it.
    std::vector<double> w(x.size(), 0.0), b(x.size(), 0.0);
    const double spacing = x[1] - x[0];
    const double force = 1.0e7;
    b[mid] = force / spacing;
    for (double& v : w) v = force / (x.back() - x.front());
    const HullGirder hog = integrateGirder(x, w, b);
    expectTrue("the imposed case hogs", hog.hogging());

    const double yield = 235e6;   // mild shipbuilding steel, Pa
    const std::vector<GirderStress> stresses = girderStress(hog, structure, yield);
    expectTrue("stress was evaluated where there is structure", !stresses.empty());

    for (const GirderStress& s : stresses) {
        if (std::abs(s.moment) < 1e-6) continue;
        // **Against `M y / I`, not against `M / Z`.** `girderStress` sets
        // `s.stressDeck = station.moment / section.modulusDeck` and copies both of
        // those into `s.moment` and `s.modulusDeck`, so asserting
        // `s.stressDeck == s.moment / s.modulusDeck` was `x == x` -- it survived any
        // `hullGirderSection` whatever, including a `modulusDeck` out by a factor of
        // ten, and it ran about forty times a suite. The section is fetched again
        // here and the stress rebuilt from the fibre and the second moment, which is
        // the definition the modulus is a shorthand for.
        const HullGirderSection at = hullGirderSection(structure, s.x);
        expectNear("deck stress is M y over I", s.stressDeck,
                   s.moment * (at.zDeck - at.neutralAxis) / at.secondMoment,
                   1e-9 * std::abs(s.stressDeck));
        expectTrue("hogging puts the deck in tension and the keel in compression",
                   s.stressDeck > 0 && s.stressKeel < 0);
    }

    // Sagging must reverse both, or the sign is decoration.
    std::vector<double> w2 = b, b2 = w;   // swap: weight concentrated amidships
    const HullGirder sag = integrateGirder(x, w2, b2);
    expectTrue("the mirrored case sags", !sag.hogging());
    const std::vector<GirderStress> sagStress = girderStress(sag, structure, yield);
    bool reversed = true;
    for (const GirderStress& s : sagStress)
        if (std::abs(s.moment) > 1e-6 && !(s.stressDeck < 0 && s.stressKeel > 0)) reversed = false;
    expectTrue("sagging compresses the deck and stretches the keel", reversed);

    // Utilisation is the worst fibre over yield, and must scale with the moment.
    double atX = 0;
    const double u = worstUtilisation(stresses, &atX);
    expectTrue("utilisation is a real fraction", u > 0);

    HullGirder doubled = hog;
    for (GirderStation& st : doubled.stations) st.moment *= 2.0;
    expectNear("doubling the moment doubles the utilisation",
               worstUtilisation(girderStress(doubled, structure, yield)), 2.0 * u, 1e-9 * u);
    expectNear("and halving the yield strength doubles it too",
               worstUtilisation(girderStress(hog, structure, 0.5 * yield)), 2.0 * u, 1e-9 * u);
}

}  // namespace

void runGirderTests() {
    std::printf("\n--- hull girder ---\n");
    testBalancedBeamCarriesNothing();
    testConcentratedWeightAmidshipsGivesWLOverEight();
    testConcentratedBuoyancyAmidshipsHogs();
    testImbalanceShowsInTheClosure();
    testWeightDistributionCarriesItsOwnTotalAndCentre();
    testBuoyancyDistributionIntegratesToDisplacement();
    testAnEvenlyLoadedBargeBarelyBends();
    testCrestHogsAndTroughSags();
    testFloodingBendsTheHull();
    testBalanceConvergesTheSameWayAtEverySize();
    testAShipWithNoHullExtentIsRefused();
    testASingularBalanceJacobianIsRefusedOnItsOwnScale();
    testTheBalanceFollowsTheCentreOfGravityAsSheTrims();
    testABalanceStillRotatingIsNotReportedAsSuccess();
    testStressFollowsTheMomentAndItsSign();
}
