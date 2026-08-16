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

Ship barge(double lcgOffset) {
    Ship s;
    std::vector<Station> stations;
    for (int i = 0; i <= 40; ++i) {
        Station q;
        q.x = -60.0 + 120.0 * i / 40.0;
        q.halfBeam = {9.0, 9.0};
        stations.push_back(q);
    }
    s.hull = makeHullFromStations(stations, {0.0, 12.0});
    s.deckEdgeZ = 12.0;
    s.lightshipMass = 120.0 * 18.0 * 5.0 * kRhoSeawater;   // floats at 5 m
    s.lightshipCog = {lcgOffset, 0.0, 6.0};
    s.gyradii = {6.0, 30.0, 30.0};
    return s;
}

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
        expectNear("deck stress is M over Z_deck", s.stressDeck, s.moment / s.modulusDeck,
                   1e-9 * std::abs(s.moment / s.modulusDeck));
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
    testStressFollowsTheMomentAndItsSign();
}
