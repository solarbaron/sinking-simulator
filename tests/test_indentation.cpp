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

// --- The ship the strike tests are driven against ------------------------------

// Starboard side, amidships, below the waterline.
const Vec3 kImpact{0.0, -9.0, 4.0};

// The reference ferry's structure. Three tests need the same hull, the same
// scantlings and the same impact point, and one of them has to find the struck bay
// *without* asking `impactDamage` which one it was -- so the fixture is built here
// rather than three times over.
StructuralMesh ferryStructure() {
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
    return makeStructuralMesh(ferry.hull, ferryScantlings());
}

// How far a reported panel is from the strike -- the same distance `impactDamage`
// sorts and bounds on, recomputed from the mesh rather than read back out of it.
double reachOf(const StructuralMesh& structure, int index) {
    return length(structure.panels[static_cast<std::size_t>(index)].centroid() - kImpact);
}

// The bay under the impact point, found by scanning the mesh. Independent of the
// strike, which is the point: the panel model `impactDamage` builds can only be
// checked against a panel someone else picked.
int nearestShellPanel(const StructuralMesh& structure) {
    int nearest = -1;
    for (std::size_t i = 0; i < structure.panels.size(); ++i) {
        if (structure.panels[i].role != PanelRole::Shell) continue;
        if (nearest < 0 || reachOf(structure, static_cast<int>(i)) < reachOf(structure, nearest))
            nearest = static_cast<int>(i);
    }
    return nearest;
}

double yieldStrengthOf(const StructuralMesh& structure, const PlatePanel& panel) {
    return panel.material >= 0 && panel.material < static_cast<int>(structure.materials.size())
               ? structure.materials[static_cast<std::size_t>(panel.material)].yieldStrength
               : ah36Steel().yieldStrength;
}

void testImpactDamageGrowsWithEnergyAndIsLocal() {
    const Scantlings sc = ferryScantlings();
    const StructuralMesh structure = ferryStructure();
    const Vec3 impact = kImpact;

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

// What the strike actually hands the panel model, against numbers derived here.
//
// **Every energy assertion in this file is blind to the span by construction.** A
// bay tears at `sigma_y t (span x contactWidth) eps_f` and `impactDamage` sets
// `contactWidth = area / span`, so the span cancels out of the product and a bay
// that spans the wrong way tears at all but the same energy. The one penetration
// assertion above is a ratio between two strikes, where it cancels again. Nothing
// here ever compared a per-panel model input against an independently derived
// value, only against another strike -- so `span = max(...)` in place of `min`, and
// a struck width that multiplies by the span where it should divide, both passed
// the whole file.
//
// They are not small. The header records what the zone FEM cost to settle: the
// hole moves 5%, because the failure strain is nearly flat over this range, but
// **penetration and resisting force move by 3.4x**, in the direction of reporting
// the hull far softer than it is. Both witnesses are below.
void testTheStrikeBuildsThePanelModelItClaimsTo() {
    const Scantlings sc = ferryScantlings();
    const StructuralMesh structure = ferryStructure();

    // The two spacings the mesh was built with. A plate spans the *short* way
    // between its supports, so on this longitudinally framed side the span is the
    // 0.70 m longitudinal spacing and not the 2.40 m frame spacing; both are named
    // so what follows rests on which of them is smaller rather than on a constant.
    expectNear("the ferry is framed transversely at 2.40 m", structure.frameSpacing, 2.40, 1e-12);
    expectNear("and longitudinally at 0.70 m", sc.longitudinalSpacing, 0.70, 1e-12);
    const double span = std::min(structure.frameSpacing, sc.longitudinalSpacing);

    // The bay under the impact, picked out of the mesh rather than read back out of
    // the answer.
    const int nearest = nearestShellPanel(structure);
    expectTrue("the impact point has a bay under it", nearest >= 0);
    const PlatePanel& bay = structure.panels[static_cast<std::size_t>(nearest)];
    const double sigmaY = yieldStrengthOf(structure, bay);

    // --- witness one: the penetration of a strike too weak to tear anything ------
    //
    // Below tearing the failure strain plays no part, so the answer is the membrane
    // inverse alone. With `w = A/L` and `h = L/2`,
    //
    //     scale = 2 sigma_y t w = 2 sigma_y t A / L
    //     u     = E / scale     = E L / (2 sigma_y t A) = c L,  c = E / (2 sigma_y t A)
    //     d     = sqrt(u (u + 2h)) = sqrt(c^2 L^2 + c L^2)     = L sqrt(c^2 + c)
    //
    // -- **the penetration is exactly proportional to the span.** The two routes the
    // span takes into the answer, through the half-span and through the struck width
    // that divides by it, pull the same way instead of cancelling, which is why this
    // is the term the mistake shows up in and the energy is not. 2.40 / 0.70 = 3.43
    // is the header's 3.4x, and it is the whole of the difference between 0.686 m of
    // reported denting and the true 0.205.
    //
    // On this bay: A = 1.7136 m2, t = 12.0 mm, sigma_y = 355 MPa, E = 200 kJ, so
    // c = 200e3 / (2 x 355e6 x 0.012 x 1.7136) = 1.36988e-2 and
    // d = 0.70 sqrt(c^2 + c) = 82.488 mm. Spanning the long way gives 282.82 mm;
    // multiplying by the span where the model divides gives 118.67 mm.
    const double energy = 2.0e5;
    const ImpactDamage nudge = impactDamage(structure, kImpact, 6.0, energy, sc);
    expectTrue("a 200 kJ nudge dents one bay and stops in it",
               nudge.torn.empty() && nudge.panels.size() == 1 && nudge.panels[0] == nearest);
    const double c = energy / (2.0 * sigmaY * bay.thickness * bay.area());
    const double expected = span * std::sqrt(c * c + c);
    std::printf("     a 200 kJ nudge dents the struck bay %.3f mm; derived %.3f mm\n",
                1e3 * nudge.penetration, 1e3 * expected);
    expectNear("the strike dents the bay by L sqrt(c^2 + c), c = E / (2 sigma_y t A)",
               nudge.penetration, expected, 1e-9 * expected);

    // --- witness two: what a torn patch costs -----------------------------------
    //
    // The other half of the same finding, and the half the header's "5%, not
    // tenfold" rests on. A bay tears at `sigma_y t (L w) eps_f` exactly -- at the
    // tearing depth sqrt(h^2 + d_f^2) = h (1 + eps_f), so the energy collapses to
    // that -- and `L w` is the bay's own area, so the span reaches the answer only
    // through the failure-strain regularisation.
    //
    // Asserted on the real ship, against an area, a yield strength and a failure
    // strain written out here. A strike that runs out of *hull* rather than out of
    // energy tears every bay it can reach and nothing else, so what it absorbed is
    // the sum of their tearing energies -- a derived total, not a bound.
    const double radius = 2.0;
    const ImpactDamage cramped = impactDamage(structure, kImpact, radius, 2.0e8, sc);
    const plasticity::Material steel = plasticity::shipSteel();
    double toTear = 0;
    long long reachable = 0;
    for (const PlatePanel& p : structure.panels) {
        if (p.role != PanelRole::Shell) continue;
        if (length(p.centroid() - kImpact) > radius) continue;
        ++reachable;
        // The regularisation written out rather than called: an element no larger
        // than the plate is thick contains the neck, a bigger one sees it diluted by
        // t/l, and `l` is the span -- the length the membrane strain is smeared
        // over. Here t/L = 0.012/0.70 = 1.7143e-2 and eps_f = 0.159564; on the frame
        // spacing it would be 5.0e-3 and 0.151670. That 5.2% is the entire cost of
        // the span error on the hole, against 3.4x on the depth of it.
        const double share = std::min(1.0, p.thickness / span);
        const double failureStrain =
            steel.failure.uniformStrain +
            (steel.failure.fractureStrain - steel.failure.uniformStrain) * share;
        toTear += yieldStrengthOf(structure, p) * p.thickness * p.area() * failureStrain;
    }
    expectTrue("the cramped strike tore every bay it could reach and still had energy left",
               cramped.torn.size() == static_cast<std::size_t>(reachable) &&
                   cramped.energyUnspent > 0);
    // 8 bays of 1.7136 m2 12 mm AH36: 8 x 355e6 x 0.012 x 1.7136 x 0.159564 =
    // 9.31848e6 J. Spanning the long way makes it 8.85748e6; multiplying by the span
    // instead of dividing makes it 0.49 of it.
    std::printf("     %lld bays inside a %.1f m reach tear for %.5e J; the strike absorbed %.5e\n",
                reachable, radius, toTear, cramped.energyAbsorbed);
    expectNear("a torn patch costs sigma_y t A eps_f, bay by bay", cramped.energyAbsorbed, toTear,
               1e-9 * toTear);
}

// The march itself: how far it may go, and in what order it goes.
void testTheMarchIsOutwardAndBoundedByTheRadius() {
    const Scantlings sc = ferryScantlings();
    const StructuralMesh structure = ferryStructure();

    // --- the radius is a bound on the reach --------------------------------------
    //
    // **The only locality check this file had compares two strikes 40 m apart**,
    // which a 6 m reach and a 12 m reach separate equally well: nothing asserted
    // that a bay the strike *reached* was inside the radius at all, so doubling the
    // radius was survived. Asserted on a strike that runs out of hull rather than
    // out of energy, because then the reached set is the whole reachable set and the
    // count is two-sided -- no bay outside, and every bay inside.
    const double radius = 2.0;
    const ImpactDamage cramped = impactDamage(structure, kImpact, radius, 2.0e8, sc);
    expectTrue("the strike has energy to spare, so it stops for want of hull and not of energy",
               cramped.energyUnspent > 0);

    long long reachable = 0;
    for (const PlatePanel& p : structure.panels)
        if (p.role == PanelRole::Shell && length(p.centroid() - kImpact) <= radius) ++reachable;

    double farthest = 0;
    bool inside = true;
    for (int index : cramped.panels) {
        const double reach = reachOf(structure, index);
        farthest = std::max(farthest, reach);
        if (reach > radius) inside = false;
    }
    std::printf("     a %.1f m reach holds %lld bays; the strike touched %zu, out to %.4f m\n",
                radius, reachable, cramped.panels.size(), farthest);
    expectTrue("no bay outside the radius is touched", inside);
    expectEqual("and every bay inside it is", static_cast<long long>(cramped.panels.size()),
                reachable);
    // "Inside the radius" is satisfied vacuously by a patch that stops well short of
    // it, so say that the patch pushes against the bound: 1.62 m of the 2.00 asked
    // for, the outermost bay centroid the mesh has inside that reach.
    expectTrue("and the patch reaches most of the way out to the bound", farthest > 0.7 * radius);

    // --- and the energy is spent outward -----------------------------------------
    //
    // The outward march was itself the fix for a measured defect -- energy shared
    // over a fixed patch by area made the hole a property of the contact radius
    // rather than of the collision -- and every assertion it arrived with is a
    // count, an area or a total. **None of them names an order**, so spending the
    // energy on the farthest bay first passed all of them.
    const ImpactDamage heavy = impactDamage(structure, kImpact, 6.0, 2.0e8, sc);
    expectTrue("the heavy strike reaches enough bays for an order to exist",
               heavy.panels.size() > 10);
    bool outward = true;
    double previous = -1.0;
    for (int index : heavy.panels) {
        const double reach = reachOf(structure, index);
        if (reach < previous) outward = false;
        previous = reach;
    }
    expectTrue("the bays are reported nearest first, which is the order they are spent in",
               outward);

    // The same statement where it decides a number rather than an ordering: a strike
    // too weak to tear anything stops in the bay *nearest* the impact, and the bays
    // within reach are not alike, so which one it picks is visible in the answer.
    const ImpactDamage nudge = impactDamage(structure, kImpact, 6.0, 2.0e5, sc);
    expectTrue("a strike that stops in one bay stops in the nearest one",
               nudge.panels.size() == 1 && nudge.panels[0] == nearestShellPanel(structure));
}

// **Which problem, not whether there was one.**
//
// The slenderness fixture here was `span = 0.1` with the reference panel's 2.00 m
// contact width left on it, and that trips *both* refusals: span/thickness of 8,
// and a contact twenty times the span where the other rule allows four. So
// `!validateIndentation(stubby).empty()` had been passing for the wrong reason --
// it asserted that something came back, and deleting the slenderness guard outright
// left it green, which is exactly what a mutation sweep found. Each fixture below
// breaks one rule and is asserted against that rule by name; the panel that breaks
// both is kept, and asserted to be told about both.
void testValidateCatchesModelAbuse() {
    expectTrue("an ordinary bay raises nothing", validateIndentation(referencePanel()).empty());

    IndentedPanel stubby = referencePanel();
    stubby.span = 0.10;          // span/thickness of 8.3: a block, not a membrane
    stubby.contactWidth = 0.30;  // and inside 4 x span, so the contact rule stays quiet
    const std::vector<std::string> stubbyProblems = validateIndentation(stubby);
    expectEqual("a span only a few thicknesses long is refused, on its own account",
                static_cast<long long>(stubbyProblems.size()), 1);
    expectTrue("and what it is told is that the bending it omits is no longer negligible",
               stubbyProblems.size() == 1 &&
                   stubbyProblems[0].find("span/thickness below 20") != std::string::npos);

    IndentedPanel wide = referencePanel();
    wide.contactWidth = 100.0;   // 42 x the span; span/thickness is a healthy 200
    const std::vector<std::string> wideProblems = validateIndentation(wide);
    expectEqual("a contact far wider than the span is refused, on its own account",
                static_cast<long long>(wideProblems.size()), 1);
    expectTrue("and what it is told is that the surrounding structure would share the load",
               wideProblems.size() == 1 &&
                   wideProblems[0].find("contact much wider than the span") != std::string::npos);

    IndentedPanel both = stubby;
    both.contactWidth = 2.00;    // the fixture as it stood: stubby *and* over-wide
    expectEqual("a bay that breaks both rules is told about both",
                static_cast<long long>(validateIndentation(both).size()), 2);
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
    testTheStrikeBuildsThePanelModelItClaimsTo();
    testTheMarchIsOutwardAndBoundedByTheRadius();
    testValidateCatchesModelAbuse();
}
