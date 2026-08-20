// SPDX-License-Identifier: MIT
//
// Validation of the rigid-plastic indentation model.
//
// Every quantity here has a closed form and three of them are inverses of each
// other, which is the strongest structure a test file can be handed: strain and
// penetration invert exactly, energy and penetration invert exactly, and the
// energy is the force integrated over the penetration. A round trip through a
// pair of inverses catches an error that a single evaluation against a
// hand-computed number would not, because the hand-computed number and the code
// can be wrong the same way and an inverse cannot.
#include "engine/sim/indentation.hpp"
#include "engine/sim/ship.hpp"
#include "harness.hpp"

#include <cmath>
#include <cstdio>
#include <vector>

using namespace sim;
using testing::expectNear;
using testing::expectTrue;

namespace {

IndentedPanel referencePanel() {
    IndentedPanel p;
    p.span = 2.40;
    p.thickness = 0.012;
    p.contactWidth = 2.00;
    p.yieldStrength = 355.0e6;
    p.failureStrain = 0.20;
    return p;
}

// The bounds below are in ulps rather than in metres, because every closed form
// here is now written without a cancelling subtraction and its error is therefore
// a fixed handful of roundings at every depth -- see the two round trips.
constexpr double kUlp = 2.220446049250313e-16;

// --- The closed forms, and their inverses ------------------------------------

void testStrainAndPenetrationInvertExactly() {
    const double span = 2.40;

    // The geometry: a tent of depth d over a half-span h stretches the plate from
    // h to sqrt(h^2 + d^2), so the strain is that ratio minus one. Assert the
    // geometry directly at a point where it is trivial -- d = h gives sqrt(2) - 1.
    expectNear("a penetration equal to the half-span strains the plate by sqrt(2) - 1",
               membraneStrain(span, 0.5 * span), std::sqrt(2.0) - 1.0, 1e-12);
    expectNear("no penetration is no strain", membraneStrain(span, 0.0), 0.0, 1e-15);

    // Round trip, over a wide range. This is the check that cannot be satisfied
    // by two matching mistakes.
    //
    // Geometric, and down to a micron, because a *linear* sweep cannot see this
    // pair's failure mode. The version that stood here ran d = 0.005 i L for
    // i = 1..200: every one of its points is within a factor of 200 of the top,
    // and both closed forms are a difference of two nearly equal lengths whose
    // conditioning goes as (h/d)^2. Its worst point was its first, 12 mm, where the
    // round trip was out by 1.4e-14 against a tolerance of 1e-12 L = 2.4e-12 --
    // 172x of headroom, and the headroom *grew* over the rest of the sweep. Run
    // against the same code, the sweep below is out by 3.7e-5 relative at its
    // small end, 1.7e11 ulps and eleven orders past the bound.
    //
    // Relative rather than absolute, for the same reason: an absolute bound is
    // slackest in exactly the corner where the arithmetic is hardest, and a round
    // trip through a pair of exact inverses is honestly measured in ulps of the
    // answer rather than in metres.
    double worst = 0, worstAt = 0;
    for (int i = 0; i <= 600; ++i) {
        const double d = 1e-6 * span * std::pow(1e6, i / 600.0);   // 2.4 um to 2.4 m
        const double back = penetrationForStrain(span, membraneStrain(span, d));
        if (std::abs(back - d) / d > worst) {
            worst = std::abs(back - d) / d;
            worstAt = d;
        }
    }
    std::printf("     strain round trip over six decades: worst %.2f ulp at d = %.3g m\n",
                worst / kUlp, worstAt);
    // 8 ulps, against a measured 1.905 (4.2e-16, at d = 0.53 m). It is allowed to
    // be this tight because neither direction cancels any more, so the error is the
    // rounding of four operations and does not grow as the dent shrinks. The 4x
    // margin is for the compiler and is more than this one needs: 1.905 is
    // bit-identical at -O0, -O1, -O2, -O3 and -Os, with -ffp-contract=fast and with
    // -march=native, which spans everything `full` and the sanitizers compile.
    expectTrue("strain and penetration invert to within a few ulps at every depth",
               worst < 8.0 * kUlp);
    expectTrue("and the sweep reached the depths where that is hard", worstAt > 0);

    // Small penetrations: the strain goes as 2 (d/L)^2, so halving the
    // penetration quarters the strain. A model that got the geometry linear
    // instead would still pass a single-point check.
    const double a = membraneStrain(span, 0.02 * span);
    const double b = membraneStrain(span, 0.01 * span);
    expectNear("membrane strain is quadratic in penetration for small dents", a / b, 4.0, 0.02);

    // The same limit asserted where it is hard rather than where it is easy, and
    // against the series to two terms so the tolerance can be arithmetic rather
    // than physical: eps = r^2/2 (1 - r^2/4 + r^4/8 - ...), and at r = 2e-6 the
    // first neglected term is 2e-24, nine orders below an ulp. Measured 0.00 ulp
    // against that; the cancelling form is out by 2.2e-5.
    const double tiny = 1e-6 * span;
    const double r = tiny / (0.5 * span);
    const double series = 0.5 * r * r * (1.0 - 0.25 * r * r);
    expectNear("a micron-deep dent strains the plate by r^2/2 (1 - r^2/4)",
               membraneStrain(span, tiny) / series, 1.0, 8.0 * kUlp);
    // And the quartering at that depth, to the same standard: the ratio is
    // 4 (1 - r^2) / (1 - r^2/4) = 4 (1 - 0.75 r^2) to O(r^4). Measured 1.00 ulp,
    // against 2.8e-5 for the form this replaced.
    expectNear("and halving a micron-deep dent still quarters it",
               membraneStrain(span, 2.0 * tiny) / membraneStrain(span, tiny) /
                   (4.0 * (1.0 - 0.75 * r * r)),
               1.0, 8.0 * kUlp);
}

void testEnergyIsTheForceIntegratedOverThePenetration() {
    const IndentedPanel p = referencePanel();

    // Quadrature of the force against the closed-form energy: two completely
    // separate routes to the same number.
    const double target = 0.35;
    const int steps = 200000;
    double quadrature = 0;
    for (int i = 0; i < steps; ++i) {
        const double d = target * (i + 0.5) / steps;
        quadrature += indentationForce(p, d) * (target / steps);
    }
    const double closedForm = indentationEnergy(p, target);
    expectNear("the closed-form energy is the force integrated", closedForm, quadrature,
               1e-6 * closedForm);

    // And energy inverts back to penetration -- geometrically down to a micron and
    // relatively, for the reason set out in the round trip above. This is the
    // worse-conditioned of the two pairs when written directly: `sqrt(h^2 + d^2)
    // - h` on the way out and `(E/scale + h)^2 - h^2` on the way back both
    // difference two O(h^2) numbers whose true difference is O(d^2). The sweep that
    // stood here -- d = 0.003 i, from 3 mm -- began four decades above the trouble
    // and asserted 1e-9 m absolute, which 3 mm meets with six orders to spare.
    double worst = 0, worstAt = 0;
    const double tearing = penetrationForStrain(p.span, p.failureStrain);
    for (int i = 0; i <= 600; ++i) {
        const double d = 1e-6 * p.span * std::pow(1e5, i / 600.0);   // 2.4 um to 0.24 m
        if (d >= tearing) continue;
        const double back = penetrationForEnergy(p, indentationEnergy(p, d));
        if (std::abs(back - d) / d > worst) {
            worst = std::abs(back - d) / d;
            worstAt = d;
        }
    }
    std::printf("     energy round trip over five decades: worst %.2f ulp at d = %.3g m\n",
                worst / kUlp, worstAt);
    // 8 ulps against a measured 1.070 (2.4e-16), bit-identical across the same
    // seven optimisation settings, and tight for the same reason. The cancelling
    // forms reach 2.6e-5 relative at the small end of this sweep, eleven orders
    // outside it.
    expectTrue("energy and penetration invert to within a few ulps at every depth",
               worst < 8.0 * kUlp);
    expectTrue("and the sweep reached the depths where that is hard", worstAt > 0);

    // The small-dent limit against its closed form. For d << h the tent is a
    // parabola and the energy is sigma_y t w d^2 / h -- E = scale h eps -- so the
    // same two-term series carries it, and the same 2e-24 remainder lets the
    // tolerance be arithmetic. Measured -0.50 ulp; written as the subtraction it is
    // out by 3.3e-5, because sqrt(h^2 + d^2) and h agree to eleven digits there.
    const double tiny = 1e-6 * p.span;
    const double h = 0.5 * p.span;
    const double r = tiny / h;
    const double scale = 2.0 * p.yieldStrength * p.thickness * p.contactWidth;
    const double eSeries = scale * h * 0.5 * r * r * (1.0 - 0.25 * r * r);
    expectNear("a micron-deep dent costs sigma_y t w d^2 / h", indentationEnergy(p, tiny) / eSeries,
               1.0, 8.0 * kUlp);
    // ...and its inverse, which is the one term this file had no check on at all
    // below 3 mm: d = sqrt(u (u + 2h)) = sqrt(2 h u) (1 + u/(4h) - ...) with
    // u = E/scale. Measured 0.00 ulp; the cancelling form is out by 5.1e-6.
    const double u = indentationEnergy(p, tiny) / scale;
    expectNear("and that energy puts the dent back at sqrt(2 h E / scale)",
               penetrationForEnergy(p, indentationEnergy(p, tiny)) /
                   (std::sqrt(2.0 * h * u) * (1.0 + 0.25 * u / h)),
               1.0, 8.0 * kUlp);
}

// The scalings are exact and each isolates one term, so a factor misplaced
// between them shows up here and nowhere else.
void testEnergyScalesWithTheThingsItShould() {
    const IndentedPanel base = referencePanel();
    const double d = 0.25;
    const double e = indentationEnergy(base, d);

    IndentedPanel thicker = base;
    thicker.thickness *= 2.0;
    expectNear("energy is linear in thickness", indentationEnergy(thicker, d), 2.0 * e, 1e-9 * e);

    IndentedPanel wider = base;
    wider.contactWidth *= 3.0;
    expectNear("and linear in the struck width", indentationEnergy(wider, d), 3.0 * e, 1e-9 * e);

    IndentedPanel stronger = base;
    stronger.yieldStrength *= 1.5;
    expectNear("and linear in yield strength", indentationEnergy(stronger, d), 1.5 * e, 1e-9 * e);

    // Span is the one that is *not* linear: a longer bay strains less for the
    // same dent, so it absorbs less. Assert the direction and that it is real.
    IndentedPanel longer = base;
    longer.span *= 2.0;
    expectTrue("a longer bay absorbs less energy for the same dent",
               indentationEnergy(longer, d) < 0.7 * e);
}

// --- Tearing -----------------------------------------------------------------

void testTearingHappensAtTheFailureStrainAndNotBefore() {
    const IndentedPanel p = referencePanel();
    const double tearing = penetrationForStrain(p.span, p.failureStrain);
    expectTrue("the tearing penetration is a real distance", tearing > 0.1 && tearing < p.span);

    expectNear("the strain at the tearing penetration is the failure strain",
               membraneStrain(p.span, tearing), p.failureStrain, 1e-12);
    expectTrue("just short of it the panel is intact", !indentAt(p, 0.999 * tearing).torn);
    expectTrue("just past it the panel has torn", indentAt(p, 1.001 * tearing).torn);

    // Energy to tear, against the closed form rather than against its own body.
    //
    // **This asserted `energyToTear(p) == indentationEnergy(p, tearing)`, which
    // is `energyToTear`'s one-line definition** (`indentation.cpp`) evaluated on
    // the same arguments in the same order: bit-exactly zero residual against a
    // 1e-9 tolerance, unfailable for any implementation that stores what it
    // divided.
    //
    // The header states an independent form — `sigma_y t A eps_f`, reached only
    // through the failure-strain regularisation — and the whole "5%, not ten
    // times" correction rests on the span cancelling out of it. That is the
    // claim worth pinning, and it was the one thing here nothing checked.
    const double energy = energyToTear(p);
    const double closedForm =
        p.yieldStrength * p.thickness * (p.span * p.contactWidth) * p.failureStrain;
    std::printf("     energy to tear %.0f J, closed form sigma_y t A eps_f %.0f J\n",
                energy, closedForm);
    expectNear("the energy to tear is sigma_y t A eps_f", energy, closedForm,
               1e-9 * closedForm);

    // **And the span really does cancel.** Holding the struck area fixed while
    // moving the span is what makes the hole a 5% question rather than a tenfold
    // one; if this ever stops holding, every published hole area moves with it.
    IndentedPanel shortSpan = p;
    shortSpan.span = 0.70;
    shortSpan.contactWidth = (p.span * p.contactWidth) / shortSpan.span;  // same area
    const double shortEnergy = energyToTear(shortSpan);
    expectNear("the same struck area tears at the same energy whatever the span",
               shortEnergy, energy, 1e-9 * energy);
    // The guard: the two spans really are different, or the line above compares
    // a panel with itself.
    expectTrue("and the two spans genuinely differ",
               std::abs(shortSpan.span - p.span) > 1.0);
    expectTrue("more energy than that does not drive it deeper by this model",
               penetrationForEnergy(p, 10.0 * energy) <= tearing + 1e-12);

    // A tougher steel tears later and absorbs more. Both, or the failure strain
    // is decorative.
    IndentedPanel tough = p;
    tough.failureStrain = 0.35;
    expectTrue("a more ductile plate tears at a deeper dent",
               penetrationForStrain(tough.span, tough.failureStrain) > tearing);
    expectTrue("and absorbs more energy doing it", energyToTear(tough) > energy);
}

// The number that decides whether a collision floods anything: how much energy a
// bay can swallow before it lets go. It should be a credible size -- megajoules
// for a ship, not kilojoules and not gigajoules.
void testEnergyToTearIsAPlausibleSizeForShipStructure() {
    const IndentedPanel p = referencePanel();
    const double energy = energyToTear(p);
    expectTrue("a 12 mm bay absorbs megajoules, not kilojoules", energy > 1.0e5);
    expectTrue("and not gigajoules", energy < 1.0e8);

    // A ship of 8984 t at 6 m/s carries about 160 MJ. The number of bays that
    // would absorb it is a sanity check on the whole model: it should be tens,
    // not one and not ten thousand.
    const double impact = 0.5 * 8.984e6 * 6.0 * 6.0;
    const double bays = impact / energy;
    std::printf("     a 12 mm bay absorbs %.3e J; a 6 m/s ferry is %.0f of them\n", energy, bays);
    // **The upper half of this used to be unreachable.** `energy > 1.0e5` two lines
    // above forces `bays < 1617`, so `bays < 5000` could never fire and only the
    // lower half carried information. The label says "tens to hundreds" and the
    // comment above says "not ten thousand"; 5000 was neither. Bounded at 1000 now,
    // which is inside what the energy assertion already implies and is therefore a
    // statement rather than a restatement.
    expectTrue("a ferry at 6 m/s is worth tens to hundreds of bays, not one",
               bays > 10.0 && bays < 1000.0);
}

// --- Against a real ship -----------------------------------------------------

void testImpactDamageGrowsWithEnergyAndIsLocal() {
    Ship ferry;
    std::vector<Station> stations;
    for (int i = 0; i <= 40; ++i) {
        Station q;
        q.x = -60.0 + 120.0 * i / 40.0;
        q.halfBeam = {9.0, 9.0};
        stations.push_back(q);
    }
    ferry.hull = makeHullFromStations(stations, {0.0, 12.0});
    ferry.deckEdgeZ = 12.0;
    ferry.lightshipMass = 120.0 * 18.0 * 5.0 * kRhoSeawater;
    ferry.lightshipCog = {0, 0, 6.0};
    ferry.gyradii = {6, 30, 30};
    ferry.initialise(0.0);

    const Scantlings sc = ferryScantlings();
    const StructuralMesh structure = makeStructuralMesh(ferry.hull, sc);
    const Vec3 impact{0.0, -9.0, 4.0};   // starboard side, amidships, below the waterline

    const ImpactDamage nudge = impactDamage(structure, impact, 6.0, 2.0e5, sc);
    const ImpactDamage light = impactDamage(structure, impact, 6.0, 2.0e6, sc);
    const ImpactDamage heavy = impactDamage(structure, impact, 6.0, 2.0e8, sc);

    expectTrue("the strike found panels", !light.panels.empty());
    expectTrue("a strike too weak to tear a single bay dents it and stops",
               nudge.torn.empty() && nudge.penetration > 0);
    expectTrue("a light strike tears something", !light.torn.empty());
    expectTrue("a heavy strike tears more", heavy.torn.size() > light.torn.size());
    expectTrue("and reaches further into her", heavy.panels.size() > light.panels.size());
    expectTrue("the hole grows with the strike", heavy.tornArea > 2.0 * light.tornArea);

    // Penetration saturates once a bay lets go, and that is by design rather than
    // by accident: past tearing the membrane model has nothing further to say, so
    // depth stops being the thing that grows and *extent* becomes it. Asserting it
    // keeps that a decision rather than a surprise.
    // Not exactly equal: each panel tears at its own depth, because the failure
    // strain is regularised on that panel's own thickness and a heavier strike
    // reaches thinner plating further out. Saturation is a bound, not an identity,
    // and asserting the identity was this test's mistake.
    expectNear("penetration saturates at the tearing depth", heavy.penetration,
               light.penetration, 0.02 * light.penetration);

    // A radius too tight silently reinstates the very defect the outward march
    // removed, so the truncation is reported rather than hidden.
    const ImpactDamage cramped = impactDamage(structure, impact, 2.0, 2.0e8, sc);
    expectTrue("a tight radius leaves energy unspent and says so", cramped.energyUnspent > 0);
    expectTrue("and it therefore tears less than a generous one",
               cramped.tornArea < heavy.tornArea);
    expectTrue("a strike that fits inside its radius spends everything",
               light.energyUnspent < 1e-6 * 2.0e6);

    // Local: nothing outside the contact radius is touched, and a strike
    // elsewhere finds different panels. Without this the model could be damaging
    // the whole ship and passing every check above.
    const ImpactDamage elsewhere =
        impactDamage(structure, {40.0, -9.0, 4.0}, 3.0, 2.0e8, sc);
    expectTrue("a strike forward finds panels too", !elsewhere.panels.empty());
    bool disjoint = true;
    for (int a : heavy.panels)
        for (int b : elsewhere.panels)
            if (a == b) disjoint = false;
    expectTrue("and they are not the same panels", disjoint);

    // Energy accounting: what the panels absorbed cannot exceed what was
    // delivered. Once they tear the model stops absorbing, so it is a bound and
    // not an equality.
    expectTrue("absorbed energy does not exceed the energy delivered",
               heavy.energyAbsorbed <= 2.0e8 * (1.0 + 1e-9));
    expectNear("and a light strike absorbs essentially all of it", light.energyAbsorbed, 2.0e6,
               0.02 * 2.0e6);

    // Nothing at all for no energy.
    expectTrue("no energy, no damage",
               impactDamage(structure, impact, 3.0, 0.0, sc).panels.empty());
}

void testValidateCatchesModelAbuse() {
    expectTrue("an ordinary bay raises nothing", validateIndentation(referencePanel()).empty());

    IndentedPanel stubby = referencePanel();
    stubby.span = 0.1;   // span/thickness of 8: a block, not a membrane
    expectTrue("a span only a few thicknesses long is refused",
               !validateIndentation(stubby).empty());

    IndentedPanel wide = referencePanel();
    wide.contactWidth = 100.0;
    expectTrue("a contact far wider than the span is refused",
               !validateIndentation(wide).empty());
}

}  // namespace

void runIndentationTests() {
    std::printf("\n--- impact indentation ---\n");
    testStrainAndPenetrationInvertExactly();
    testEnergyIsTheForceIntegratedOverThePenetration();
    testEnergyScalesWithTheThingsItShould();
    testTearingHappensAtTheFailureStrainAndNotBefore();
    testEnergyToTearIsAPlausibleSizeForShipStructure();
    testImpactDamageGrowsWithEnergyAndIsLocal();
    testValidateCatchesModelAbuse();
}
