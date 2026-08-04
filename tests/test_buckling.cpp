// SPDX-License-Identifier: MIT
//
// Validation of the buckling checks.
//
// Plate and column buckling are among the very few things in structural
// mechanics with exact analytic answers, so almost nothing here is a tolerance
// against a plausible number. The buckling coefficient is a minimum of a convex
// expression and equals exactly 4 at every whole aspect ratio; the Euler column
// load is exact; the plasticity cap is continuous at its transition by
// construction, and a discontinuity there would be a real defect that no single
// evaluation could see.
#include "engine/sim/buckling.hpp"
#include "harness.hpp"

#include <cmath>
#include <cstdio>
#include <vector>

using namespace sim;
using testing::expectNear;
using testing::expectTrue;

namespace {

// --- The coefficient ---------------------------------------------------------

void testPlateCoefficientIsFourAtWholeAspectRatios() {
    // k = min over m of (m/alpha + alpha/m)^2. At alpha = m the two terms are
    // both 1, so k = 4 exactly -- not approximately, and at every whole ratio.
    for (double alpha : {1.0, 2.0, 3.0, 4.0, 7.0}) {
        expectNear("k is exactly 4 at aspect ratio " + std::to_string(alpha),
                   plateBucklingCoefficient(alpha, 1.0), 4.0, 1e-12);
    }
    // And it is symmetric in which side is which only through the caller: the
    // coefficient itself is a function of alpha, so 1/2 is not 2.
    expectTrue("the coefficient rises between whole ratios",
               plateBucklingCoefficient(std::sqrt(2.0), 1.0) > 4.0);

    // The crossover between one and two half-waves is at alpha = sqrt(2), where
    // both give (1/sqrt2 + sqrt2)^2 = 4.5. That is the worst case for a plate.
    expectNear("the crossover peak is 4.5", plateBucklingCoefficient(std::sqrt(2.0), 1.0), 4.5,
               1e-12);

    // A long plate tends to 4 from above, never below.
    double worst = 0;
    for (int i = 1; i <= 400; ++i) {
        const double k = plateBucklingCoefficient(1.0 + 0.05 * i, 1.0);
        expectTrue("k never falls below 4", k >= 4.0 - 1e-12);
        worst = std::max(worst, k);
    }
    expectTrue("and never exceeds the crossover peak", worst <= 4.5 + 1e-12);
    expectNear("a long plate tends to 4", plateBucklingCoefficient(50.0, 1.0), 4.0, 0.01);
}

// --- The plasticity cap ------------------------------------------------------

void testJohnsonOstenfeldIsContinuousAndBounded() {
    const double yield = 355e6;

    // Below half yield it does nothing at all.
    for (double s : {10e6, 100e6, 0.5 * yield})
        expectNear("elastic buckling below half yield is unchanged", johnsonOstenfeld(s, yield), s,
                   1e-9);

    // Continuity at the transition, stated as the closed form rather than as a
    // numerical limit. At sigma_cr = sigma_y / 2 the lower branch returns the
    // input and the upper returns sigma_y (1 - sigma_y / (4 * sigma_y/2)) =
    // sigma_y / 2: the same number exactly, which is what makes the correction
    // continuous by construction rather than by luck.
    //
    // The first version of this check compared two points 1 Pa either side of the
    // transition and demanded they agree to 1e-3 Pa. The function has slope 1
    // there, so they legitimately differ by 2 Pa -- that assertion was demanding a
    // constant, not continuity, and the code was right.
    expectNear("both branches meet exactly at half yield", johnsonOstenfeld(0.5 * yield, yield),
               0.5 * yield, 1e-9);
    const double eps = 1.0;
    expectNear("and the step across the transition is the slope, not a jump",
               johnsonOstenfeld(0.5 * yield + eps, yield) - johnsonOstenfeld(0.5 * yield - eps, yield),
               2.0 * eps, 0.05 * eps);

    // Bounded above by yield, and it approaches it. This is the whole point: an
    // infinitely stiff plate is still only as strong as its material.
    expectTrue("the cap never exceeds yield", johnsonOstenfeld(1e15, yield) <= yield);
    expectNear("and approaches it for a very stiff plate", johnsonOstenfeld(1e5 * yield, yield),
               yield, 1e-4 * yield);
    expectTrue("the correction is monotone",
               johnsonOstenfeld(2 * yield, yield) > johnsonOstenfeld(yield, yield));
}

// --- Plate buckling ----------------------------------------------------------

// The formula is arithmetic, so it is asserted as arithmetic -- against a value
// computed here from the constants, not read off the implementation.
void testPlateBucklingMatchesTheClosedForm() {
    const StructuralMaterial steel = ah36Steel();
    const double t = 0.012, b = 0.70, a = 2.40;

    const BucklingCheck c = plateBuckling(t, a, b, 50e6, steel);
    const double k = plateBucklingCoefficient(a / b, 1.0);
    const double want = k * kPi * kPi * steel.youngsModulus /
                        (12.0 * (1.0 - steel.poissonRatio * steel.poissonRatio)) * (t / b) * (t / b);
    expectNear("the elastic stress is k pi^2 E t^2 / (12 (1-nu^2) b^2)", c.elasticStress, want,
               1e-6 * want);
    expectNear("and the coefficient is the one for this aspect ratio", c.coefficient, k, 1e-12);

    // The short side is b whichever way round the caller names the sides.
    const BucklingCheck swapped = plateBuckling(t, b, a, 50e6, steel);
    expectNear("naming the sides the other way round gives the same answer", swapped.elasticStress,
               c.elasticStress, 1e-9 * c.elasticStress);

    // Thickness squared, and spacing squared the other way. Both are exact.
    const BucklingCheck thicker = plateBuckling(2 * t, a, b, 50e6, steel);
    expectNear("doubling the thickness quadruples the elastic stress", thicker.elasticStress,
               4.0 * c.elasticStress, 1e-6 * c.elasticStress);
    const BucklingCheck wider = plateBuckling(t, a, 2 * b, 50e6, steel);
    expectTrue("and doubling the stiffener spacing weakens it",
               wider.elasticStress < 0.5 * c.elasticStress);
}

// The reason the module exists: thin ship plating buckles a long way below yield,
// so a stress check alone passes a panel that is already gone.
void testThinPlatingBucklesFarBelowYield() {
    const StructuralMaterial steel = ah36Steel();
    // 10 mm plating on 800 mm spacing -- ordinary ship deck structure.
    const BucklingCheck thin = plateBuckling(0.010, 2.40, 0.80, 0.0, steel);
    expectTrue("ordinary deck plating buckles below yield",
               thin.criticalStress < 0.6 * steel.yieldStrength);

    // A stress that a yield check calls safe can still buckle the panel. That is
    // the whole claim, so it is asserted directly rather than implied.
    const double applied = 0.45 * steel.yieldStrength;
    const BucklingCheck loaded = plateBuckling(0.010, 2.40, 0.80, applied, steel);
    expectTrue("a panel at 45% of yield is nonetheless overloaded in buckling",
               loaded.utilisation > 1.0);

    // Thick plating on close spacing does not: the cap has to bite, or the module
    // would report every heavy structure as infinitely strong.
    const BucklingCheck thick = plateBuckling(0.040, 2.40, 0.60, 0.0, steel);
    expectTrue("heavy plating is capped at yield rather than running away",
               thick.criticalStress <= steel.yieldStrength);
    expectTrue("and is genuinely strong", thick.criticalStress > 0.9 * steel.yieldStrength);
}

// --- Column buckling ---------------------------------------------------------

void testColumnBucklingIsEuler() {
    const StructuralMaterial steel = ah36Steel();
    const StiffenedSection section = stiffenedSection(flatBar(0.200, 0.010), 0.012, 0.70);
    expectTrue("the stiffened section is real", section.area > 0 && section.secondMoment > 0);

    const double span = 2.40;
    const BucklingCheck c = columnBuckling(section, span, 0.0, steel);
    const double want = kPi * kPi * steel.youngsModulus * section.secondMoment /
                        (section.area * span * span);
    expectNear("the Euler stress is pi^2 E I / (A L^2)", c.elasticStress, want, 1e-9 * want);

    // Inverse square in the span, exactly.
    const BucklingCheck longer = columnBuckling(section, 2 * span, 0.0, steel);
    expectNear("doubling the span quarters the elastic stress", longer.elasticStress,
               0.25 * c.elasticStress, 1e-9 * c.elasticStress);

    // A stiffener between frames is short and stocky, so the cap should bite --
    // if it does not, the check is reporting an eigenvalue rather than a strength.
    expectTrue("a stiffener over a frame spacing is limited by yield, not by Euler",
               c.criticalStress < c.elasticStress);
}

// --- Hull girder scale -------------------------------------------------------

// Buckling is checked on the *compressed* fibre, and which fibre that is follows
// the sign of the moment. Getting it backwards checks the deck when the keel is
// the one folding.
void testGirderBucklingFollowsTheCompressedFibre() {
    const Scantlings sc = ferryScantlings();

    Ship ship;
    std::vector<Station> stations;
    for (int i = 0; i <= 40; ++i) {
        Station q;
        q.x = -60.0 + 120.0 * i / 40.0;
        q.halfBeam = {9.0, 9.0};
        stations.push_back(q);
    }
    ship.hull = makeHullFromStations(stations, {0.0, 12.0});
    ship.deckEdgeZ = 12.0;
    ship.lightshipMass = 120.0 * 18.0 * 5.0 * kRhoSeawater;
    ship.lightshipCog = {0, 0, 6.0};
    ship.gyradii = {6, 30, 30};
    ship.initialise(0.0);

    const StructuralMesh structure = makeStructuralMesh(ship.hull, sc);
    const std::vector<double> x = girderStations(ship, 41);
    const std::size_t mid = x.size() / 2;
    const double spacing = x[1] - x[0];
    const double force = 2.0e7;

    // Hogging: buoyancy concentrated amidships. Compresses the keel.
    std::vector<double> w(x.size(), force / (x.back() - x.front())), b(x.size(), 0.0);
    b[mid] = force / spacing;
    const HullGirder hog = integrateGirder(x, w, b);
    expectTrue("the hogging case really hogs", hog.hogging());

    const auto hogStress = girderStress(hog, structure, ah36Steel().yieldStrength);
    const auto hogBuckle = girderBuckling(hogStress, structure, sc);
    expectTrue("hogging produces buckling checks", !hogBuckle.empty());
    bool keelChecked = true;
    for (const GirderBuckling& g : hogBuckle)
        if (g.deckInCompression) keelChecked = false;
    expectTrue("hogging compresses the keel, so the keel is what is checked", keelChecked);

    // Sagging is the mirror, and must check the deck.
    const HullGirder sag = integrateGirder(x, b, w);
    expectTrue("the sagging case really sags", !sag.hogging());
    const auto sagBuckle =
        girderBuckling(girderStress(sag, structure, ah36Steel().yieldStrength), structure, sc);
    expectTrue("sagging produces checks", !sagBuckle.empty());
    bool deckChecked = true;
    for (const GirderBuckling& g : sagBuckle)
        if (!g.deckInCompression) deckChecked = false;
    expectTrue("sagging compresses the deck, so the deck is what is checked", deckChecked);

    // Utilisation must scale with the load, and a ship in pure tension must have
    // nothing to report -- the guard against a check that always fires.
    double atX = 0;
    const double u = worstBucklingUtilisation(hogBuckle, &atX);
    expectTrue("something is loaded", u > 0);

    HullGirder doubled = hog;
    for (GirderStation& s : doubled.stations) s.moment *= 2.0;
    const auto twice = girderBuckling(girderStress(doubled, structure, ah36Steel().yieldStrength),
                                      structure, sc);
    expectNear("doubling the moment doubles the buckling utilisation",
               worstBucklingUtilisation(twice), 2.0 * u, 1e-6 * u);

    HullGirder none = hog;
    for (GirderStation& s : none.stations) s.moment = 0.0;
    expectTrue("a hull carrying no moment has nothing to buckle",
               girderBuckling(girderStress(none, structure, ah36Steel().yieldStrength), structure,
                              sc)
                   .empty());
}

}  // namespace

void runBucklingTests() {
    std::printf("\n--- buckling ---\n");
    testPlateCoefficientIsFourAtWholeAspectRatios();
    testJohnsonOstenfeldIsContinuousAndBounded();
    testPlateBucklingMatchesTheClosedForm();
    testThinPlatingBucklesFarBelowYield();
    testColumnBucklingIsEuler();
    testGirderBucklingFollowsTheCompressedFibre();
}
