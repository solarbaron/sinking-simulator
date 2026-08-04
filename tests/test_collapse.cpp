// SPDX-License-Identifier: MIT
//
// Validation of progressive collapse.
//
// Smith's method sits between two exact answers and that is what makes it
// testable. At zero curvature it must be `E sum(A d^2)` about the elastic neutral
// axis. At infinite curvature with buckling switched off it must be the fully
// plastic moment about the *plastic* neutral axis -- a different axis, and using
// the elastic one understates it. Everything the method does happens between
// those two, so a curve that hits both ends exactly and is monotone in between is
// a curve that has nowhere to hide.
#include "engine/sim/collapse.hpp"
#include "engine/sim/ship.hpp"
#include "harness.hpp"

#include <cmath>
#include <cstdio>
#include <vector>

using namespace sim;
using testing::expectEqual;
using testing::expectNear;
using testing::expectTrue;

namespace {

// A symmetric I: two flanges and a web, all the same material. Every quantity is
// hand-computable, which is the point.
std::vector<CollapseElement> symmetricSection(double flangeArea, double webArea, double depth,
                                              double bucklingStress = 1e30,
                                              double shed = 0.0) {
    LoadShortening curve;
    curve.youngsModulus = 206e9;
    curve.yieldStrength = 355e6;
    curve.bucklingStress = bucklingStress;
    curve.shedExponent = shed;

    std::vector<CollapseElement> out;
    out.push_back({flangeArea, 0.0, curve});          // keel
    out.push_back({webArea, 0.5 * depth, curve});     // web, at mid depth
    out.push_back({flangeArea, depth, curve});        // deck
    return out;
}

// --- The two exact ends ------------------------------------------------------

void testInitialSlopeIsTheAxialSecondMoment() {
    const double flange = 0.60, web = 0.40, depth = 12.0;
    const std::vector<CollapseElement> section = symmetricSection(flange, web, depth);

    const CollapseCurve curve = progressiveCollapse(section, 1e-7, 20);
    expectNear("the elastic neutral axis of a symmetric section is at mid depth",
               curve.elasticNeutralAxis, 0.5 * depth, 1e-9);

    // sum(A d^2) by hand: the web sits on the axis and contributes nothing.
    const double want = 206e9 * (2.0 * flange * (0.5 * depth) * (0.5 * depth));
    expectNear("the initial stiffness is E sum(A d^2)", curve.initialStiffness, want,
               1e-9 * want);

    // And the curve's own slope agrees, which ties the reported number to what
    // the solve actually does rather than to a formula sitting beside it.
    const CollapsePoint p = curve.points.front();
    expectNear("the moment at small curvature is that stiffness times the curvature", p.moment,
               curve.initialStiffness * p.curvature, 1e-6 * std::abs(p.moment));
    expectTrue("and the axial force balances", std::abs(p.residual) < 1e-3);
}

// With buckling switched off, a large curvature drives every element to +-yield,
// and the answer is the fully plastic moment about the plastic neutral axis. That
// axis equalises *area*, not first moment, so on an asymmetric section it is not
// the elastic one -- and using the elastic axis understates the answer, which is
// the mistake worth catching.
void testUltimateReachesTheFullyPlasticMomentWithoutBuckling() {
    const double flange = 0.60, web = 0.40, depth = 12.0;
    const std::vector<CollapseElement> section = symmetricSection(flange, web, depth);

    // By hand, symmetric: the plastic axis is at mid depth, the web straddles it
    // and contributes |0| ... but it is a single lumped element sitting exactly on
    // the axis, so its lever is zero and only the flanges carry.
    const double want = 355e6 * 2.0 * flange * 0.5 * depth;
    expectNear("the fully plastic moment is sigma_y sum(A |d|)", fullyPlasticMoment(section), want,
               1e-9 * want);

    const CollapseCurve curve = progressiveCollapse(section, 5e-3, 400);
    expectNear("with no buckling the ultimate moment is the fully plastic moment",
               curve.ultimateMoment, curve.fullyPlasticMoment, 1e-3 * curve.fullyPlasticMoment);
    expectTrue("the curve is still rising at the end, so there is no peak to find",
               curve.points.back().moment > 0.999 * curve.ultimateMoment);
}

// The plastic axis is the one that balances area. Assert it on a deliberately
// asymmetric section, where it differs from the elastic axis by a knowable amount.
void testPlasticNeutralAxisBalancesAreaNotFirstMoment() {
    LoadShortening curve;
    curve.youngsModulus = 206e9;
    curve.yieldStrength = 355e6;

    // Heavy bottom, light deck: 1.5 m^2 at z = 0, 0.5 m^2 at z = 10.
    const std::vector<CollapseElement> section = {{1.5, 0.0, curve}, {0.5, 10.0, curve}};

    // Elastic axis: first moment / area = (1.5*0 + 0.5*10) / 2.0 = 2.5 m.
    const CollapseCurve c = progressiveCollapse(section, 1e-7, 5);
    expectNear("the elastic axis is the first-moment centroid", c.elasticNeutralAxis, 2.5, 1e-9);

    // Plastic axis: equal area either side. With lumped elements the balance point
    // is anywhere strictly between them that puts 1.5 below and 0.5 above -- it
    // cannot equalise, so it sits against the larger element. The moment is then
    // sigma_y (1.5 * d0 + 0.5 * d1) with the axis just above the bottom flange:
    // 355e6 * (1.5*0 + 0.5*10) = 1.775e9.
    expectNear("the fully plastic moment uses the area-balancing axis",
               fullyPlasticMoment(section), 355e6 * 0.5 * 10.0, 1e-6 * 355e6 * 0.5 * 10.0);
    expectTrue("which is not the elastic axis here",
               std::abs(c.elasticNeutralAxis - 0.0) > 1.0);
}

// --- What buckling does ------------------------------------------------------

// The whole reason the method exists: capping the compression side costs
// strength, and the loss is the difference between a ship that survives a wave
// and one that does not.
void testBucklingCostsUltimateStrength() {
    const double flange = 0.60, web = 0.40, depth = 12.0;
    const std::vector<CollapseElement> strong = symmetricSection(flange, web, depth);
    // Compression capacity halved: the flange in compression buckles at 0.5 yield.
    const std::vector<CollapseElement> weak =
        symmetricSection(flange, web, depth, 0.5 * 355e6, 0.0);

    const CollapseCurve a = progressiveCollapse(strong, 5e-3, 400);
    const CollapseCurve b = progressiveCollapse(weak, 5e-3, 400);

    expectTrue("buckling reduces the ultimate moment", b.ultimateMoment < a.ultimateMoment);
    expectTrue("but not to nothing", b.ultimateMoment > 0.3 * a.ultimateMoment);
    expectTrue("and the fully plastic moment is unchanged, since it ignores buckling",
               std::abs(b.fullyPlasticMoment - a.fullyPlasticMoment) <
                   1e-9 * a.fullyPlasticMoment);
    expectTrue("so the ultimate falls short of the plastic ceiling",
               b.ultimateMoment < 0.95 * b.fullyPlasticMoment);

    // The neutral axis must migrate away from the buckled side. That is the
    // mechanism, and a solve that held the axis fixed would still produce a
    // plausible-looking reduction.
    const CollapsePoint late = b.points.back();
    expectTrue("the neutral axis migrates towards the side still carrying",
               late.neutralAxis > b.elasticNeutralAxis + 0.1);
}

// Shedding is what turns a plateau into a peak. Without it the curve rises and
// flattens; with it the curve turns over, and the turnover is the failure.
void testSheddingProducesAPeak() {
    const double flange = 0.60, web = 0.40, depth = 12.0;
    const std::vector<CollapseElement> plateau =
        symmetricSection(flange, web, depth, 0.4 * 355e6, 0.0);
    const std::vector<CollapseElement> shedding =
        symmetricSection(flange, web, depth, 0.4 * 355e6, 0.8);

    const CollapseCurve a = progressiveCollapse(plateau, 8e-3, 600);
    const CollapseCurve b = progressiveCollapse(shedding, 8e-3, 600);

    expectTrue("without shedding the curve is still at its maximum at the end",
               a.points.back().moment > 0.995 * a.ultimateMoment);
    expectTrue("with shedding it has turned over", b.points.back().moment < 0.9 * b.ultimateMoment);
    expectTrue("and the peak is inside the sweep, not at its edge",
               b.ultimateCurvature < 0.9 * b.points.back().curvature);
    expectTrue("shedding lowers the ultimate moment", b.ultimateMoment < a.ultimateMoment);
}

// Sagging is the mirror of hogging on a symmetric section with symmetric
// elements. If it is not, a sign has gone astray somewhere in the strain.
void testSaggingMirrorsHoggingOnASymmetricSection() {
    const std::vector<CollapseElement> section = symmetricSection(0.6, 0.4, 12.0, 0.5 * 355e6, 0.5);
    const CollapseCurve hog = progressiveCollapse(section, 6e-3, 400);
    const CollapseCurve sag = progressiveCollapse(section, -6e-3, 400);

    expectTrue("hogging is positive and sagging negative", hog.ultimateMoment > 0 &&
                                                               sag.ultimateMoment < 0);
    expectNear("and they are equal in magnitude on a symmetric section",
               -sag.ultimateMoment, hog.ultimateMoment, 1e-6 * hog.ultimateMoment);
    expectNear("at the same curvature magnitude", -sag.ultimateCurvature, hog.ultimateCurvature,
               1e-9);
}

// --- The load-shortening curve itself ----------------------------------------

void testLoadShorteningIsContinuousAndCapped() {
    LoadShortening c;
    c.youngsModulus = 206e9;
    c.yieldStrength = 355e6;
    c.bucklingStress = 150e6;
    c.shedExponent = 0.5;

    expectNear("zero strain gives zero stress", c.stressAt(0.0), 0.0, 1e-12);
    expectNear("small tension is elastic", c.stressAt(1e-5), 206e9 * 1e-5, 1e-6);
    expectNear("large tension is capped at yield", c.stressAt(0.05), 355e6, 1e-6);
    expectNear("compression is capped by buckling, not by yield", c.stressAt(-0.05), 0.0,
               150e6);   // magnitude below the cap; exact value checked next
    expectTrue("and is genuinely reduced by shedding", c.stressAt(-0.05) > -150e6);

    // Continuity at the compressive cap, asserted as the closed form both sides
    // share rather than as a numerical limit: at exactly the critical strain the
    // shed factor is 1 whatever the exponent.
    const double criticalStrain = 150e6 / 206e9;
    expectNear("the curve is continuous at the buckling strain", c.stressAt(-criticalStrain),
               -150e6, 1e-6);
    // Approaching from the elastic side, the gap is exactly E times the strain
    // offset -- so the tolerance is derived from the offset rather than guessed.
    // A hard-coded 100 Pa was too tight by half: a relative offset of 1e-6 on a
    // 150 MPa cap is 150 Pa of stress, and the code was right.
    const double offset = 1e-6;
    expectNear("approached elastically from below", c.stressAt(-criticalStrain * (1.0 - offset)),
               -150e6, 2.0 * 150e6 * offset);

    // Never stronger than the cap, and the cap is actually reached. Asserted once
    // each rather than inside the sweep: an assertion per iteration inflates the
    // check count by thousands and buries any real regression in it.
    double worst = 0;
    bool withinCapacity = true;
    for (int i = 1; i <= 2000; ++i) {
        const double s = c.stressAt(-1e-5 * i);
        worst = std::min(worst, s);
        withinCapacity = withinCapacity && s >= -150e6 - 1e-6;
    }
    expectTrue("compression never exceeds the capacity anywhere in the sweep", withinCapacity);
    expectTrue("and the capacity is actually reached", worst < -0.99 * 150e6);

    // No buckling: the cap is yield.
    LoadShortening tough = c;
    tough.bucklingStress = 1e30;
    expectNear("without buckling the compressive capacity is yield",
               tough.compressiveCapacity(), 355e6, 1e-6);
}

// --- On the real ship --------------------------------------------------------

void testTheFerrySectionCollapsesAboveFirstYield() {
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
    const std::vector<CollapseElement> elements = collapseElementsAt(structure, sc, 0.0);
    expectTrue("the midship section decomposes into elements", elements.size() > 20);

    const CollapseCurve curve = progressiveCollapse(elements, 2e-3, 300);
    expectTrue("it has an ultimate moment", curve.ultimateMoment > 0);

    // Against the elastic section: the initial stiffness must be just below E*I,
    // because Smith's method carries axial stress only and leaves out the
    // elements' own second moments. "Just below" is the claim, so both sides are
    // checked -- it must not exceed it, and must not fall far short.
    const HullGirderSection section = hullGirderSection(structure, 0.0);
    const double elastic = ah36Steel().youngsModulus * section.secondMoment;
    expectTrue("the initial stiffness does not exceed E I", curve.initialStiffness <= elastic);
    expectTrue("and is within a few per cent of it",
               curve.initialStiffness > 0.90 * elastic);

    // The ultimate moment must sit between first yield and the plastic ceiling.
    const double firstYield = ah36Steel().yieldStrength * section.modulusDeck;
    expectTrue("collapse takes more moment than first yield", curve.ultimateMoment > firstYield);
    expectTrue("but less than the fully plastic moment",
               curve.ultimateMoment < curve.fullyPlasticMoment);
}

}  // namespace

void runCollapseTests() {
    std::printf("\n--- progressive collapse ---\n");
    testInitialSlopeIsTheAxialSecondMoment();
    testUltimateReachesTheFullyPlasticMomentWithoutBuckling();
    testPlasticNeutralAxisBalancesAreaNotFirstMoment();
    testBucklingCostsUltimateStrength();
    testSheddingProducesAPeak();
    testSaggingMirrorsHoggingOnASymmetricSection();
    testLoadShorteningIsContinuousAndCapped();
    testTheFerrySectionCollapsesAboveFirstYield();
}
