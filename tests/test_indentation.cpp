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
    double worst = 0;
    for (int i = 1; i <= 200; ++i) {
        const double d = 0.005 * i * span;
        const double back = penetrationForStrain(span, membraneStrain(span, d));
        worst = std::max(worst, std::abs(back - d));
    }
    expectTrue("strain and penetration invert to machine precision", worst < 1e-12 * span);

    // Small penetrations: the strain goes as 2 (d/L)^2, so halving the
    // penetration quarters the strain. A model that got the geometry linear
    // instead would still pass a single-point check.
    const double a = membraneStrain(span, 0.02 * span);
    const double b = membraneStrain(span, 0.01 * span);
    expectNear("membrane strain is quadratic in penetration for small dents", a / b, 4.0, 0.02);
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

    // And energy inverts back to penetration.
    double worst = 0;
    for (int i = 1; i <= 100; ++i) {
        const double d = 0.003 * i;
        const double back = penetrationForEnergy(p, indentationEnergy(p, d));
        if (d < penetrationForStrain(p.span, p.failureStrain))
            worst = std::max(worst, std::abs(back - d));
    }
    expectTrue("energy and penetration invert to machine precision", worst < 1e-9);
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
    expectTrue("a ferry at 6 m/s is worth tens to hundreds of bays, not one",
               bays > 10.0 && bays < 5000.0);
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
