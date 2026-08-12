// SPDX-License-Identifier: MIT
//
// Water escalation: `engine/sim/water_promotion.{hpp,cpp}`.
//
// The WaterPromoter decides which flooded compartments deserve flip::Solver's
// resolved particle-based flow instead of the lumped Compartment::waterVolume
// model. Following the pattern established in test_promotion.cpp for structural
// zones and gas compartments, the tests here separate:
//
//   1. **The criterion.** Motion thresholds (roll rate, lateral acceleration)
//      crossed with geometric guards (minimum depth, minimum volume). Tested
//      against negative controls: a ship at rest promotes nothing, shallow
//      puddles promote nothing, motion below threshold promotes nothing.
//   2. **The hysteresis.** Dwell and hold counters prevent chatter. Tested
//      against its own negative control: the same signal with hysteresis
//      switched off has to chatter, or the test is vacuous.
//   3. **The budget.** Particle count and tile count, shared across all active
//      compartments. A candidate that would exceed the budget is rejected even
//      if it qualifies on motion.
//   4. **State transfer.** Round-trip mass conservation: Compartment →
//      flip::Solver → Compartment returns the exact mass with 0.0 tolerance,
//      not epsilon-based.

#include "harness.hpp"

#include "engine/sim/ship.hpp"
#include "engine/sim/water_promotion.hpp"
#include "game/prototype/ferry.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>

using namespace sim;
using namespace sim::promotion;
using namespace testing;

// ===========================================================================
// Fixtures
// ===========================================================================

// Reuse the ferry fixture from test_promotion.cpp pattern.
Ship ferryAfloat() {
    Ship ferry = game::buildFerry();
    ferry.initialise(0.0);
    return ferry;
}

// Flood a specific compartment to a given volume.
void floodCompartment(Ship& ship, int index, double volumeM3) {
    if (index >= 0 && index < static_cast<int>(ship.compartments.size())) {
        ship.compartments[index].waterVolume = volumeM3;
    }
}

// Set ship roll rate (angular velocity about x-axis).
void setRollRate(Ship& ship, double radPerSec) {
    ship.state.angularVelocity.x = radPerSec;
}

// Set ship lateral acceleration magnitude (proxy via velocity).
void setLateralAccel(Ship& ship, double mPerSec2) {
    ship.state.velocity = Vec3{0, mPerSec2, 0};
}

// ===========================================================================
// Section 1: The Criterion
// ===========================================================================

void testShipAtRestPromotesNothing() {
    std::printf("\n   ship at rest promotes nothing\n");

    Ship ferry = ferryAfloat();
    floodCompartment(ferry, 0, 5.0);  // 5 m³, well above minVolume

    // Ship at rest: zero motion
    setRollRate(ferry, 0.0);
    setLateralAccel(ferry, 0.0);

    WaterCriterion crit;
    WaterPromoter promoter(crit);

    WaterReview rev = promoter.review(ferry);

    expectTrue("ship at rest has candidates considered", rev.considered.size() > 0);
    expectTrue("but all have score 0", rev.considered[0].score == 0.0);
    expectEqual("no promotions", static_cast<int>(rev.promoted.size()), 0);
}

void testShallowPuddlePromotesNothing() {
    std::printf("\n   shallow puddle promotes nothing\n");

    Ship ferry = ferryAfloat();
    floodCompartment(ferry, 0, 0.1);  // 0.1 m³, below minVolume

    // High roll rate, should qualify on motion
    setRollRate(ferry, 0.10);  // 2x threshold

    WaterCriterion crit;
    WaterPromoter promoter(crit);

    WaterReview rev = promoter.review(ferry);

    expectTrue("shallow puddle is considered", rev.considered.size() > 0);
    expectTrue("but fails geometric guards", rev.considered[0].score == 0.0);
    expectEqual("no promotions", static_cast<int>(rev.promoted.size()), 0);
}

void testMotionBelowThresholdPromotesNothing() {
    std::printf("\n   motion below threshold promotes nothing\n");

    Ship ferry = ferryAfloat();
    floodCompartment(ferry, 0, 5.0);  // 5 m³, good volume

    // Motion just below both thresholds
    setRollRate(ferry, 0.04);   // Below rollRatePromote = 0.05
    setLateralAccel(ferry, 1.5); // Below accelPromote = 2.0

    WaterCriterion crit;
    WaterPromoter promoter(crit);

    WaterReview rev = promoter.review(ferry);

    expectTrue("motion below threshold has candidates", rev.considered.size() > 0);
    expectTrue("but score is 0", rev.considered[0].score == 0.0);
    expectEqual("no promotions", static_cast<int>(rev.promoted.size()), 0);
}

void testRollRateAloneQualifies() {
    std::printf("\n   roll rate alone qualifies\n");

    Ship ferry = ferryAfloat();
    floodCompartment(ferry, 0, 5.0);

    // Roll rate above threshold, accel below
    setRollRate(ferry, 0.10);   // 2x rollRatePromote
    setLateralAccel(ferry, 0.5); // Below accelPromote

    WaterCriterion crit;
    WaterPromoter promoter(crit);

    WaterReview rev = promoter.review(ferry);

    expectTrue("roll rate alone produces candidates", rev.considered.size() > 0);
    expectTrue("with positive score", rev.considered[0].score > 0.0);
    expectNear("score matches roll ratio", rev.considered[0].score, 0.10 / 0.05, 1e-6);
}

void testLateralAccelAloneQualifies() {
    std::printf("\n   lateral accel alone qualifies\n");

    Ship ferry = ferryAfloat();
    floodCompartment(ferry, 0, 5.0);

    // Accel above threshold, roll below
    setRollRate(ferry, 0.02);   // Below rollRatePromote
    setLateralAccel(ferry, 3.0); // 1.5x accelPromote

    WaterCriterion crit;
    WaterPromoter promoter(crit);

    WaterReview rev = promoter.review(ferry);

    expectTrue("lateral accel alone produces candidates", rev.considered.size() > 0);
    expectTrue("with positive score", rev.considered[0].score > 0.0);
    expectNear("score matches accel ratio", rev.considered[0].score, 3.0 / 2.0, 1e-6);
}

void testBothTriggersScoreIsWeaker() {
    std::printf("\n   both triggers: score is the weaker\n");

    Ship ferry = ferryAfloat();
    floodCompartment(ferry, 0, 5.0);

    // Both above threshold, roll stronger
    setRollRate(ferry, 0.15);   // 3x rollRatePromote
    setLateralAccel(ferry, 4.0); // 2x accelPromote

    WaterCriterion crit;
    WaterPromoter promoter(crit);

    WaterReview rev = promoter.review(ferry);

    expectTrue("both triggers produce candidates", rev.considered.size() > 0);
    // Score should be min(3.0, 2.0) = 2.0
    expectNear("score is weaker trigger", rev.considered[0].score, 2.0, 1e-6);
}

// ===========================================================================
// Section 2: Hysteresis
// ===========================================================================

void testDwellPreventsImmediatePromotion() {
    std::printf("\n   dwell prevents immediate promotion\n");

    Ship ferry = ferryAfloat();
    floodCompartment(ferry, 0, 5.0);
    setRollRate(ferry, 0.10);

    WaterCriterion crit;
    crit.dwell = 3;  // Must qualify 3 times
    WaterPromoter promoter(crit);

    // First review: qualifies but not promoted yet
    WaterReview rev1 = promoter.review(ferry);
    expectTrue("first review has candidate", rev1.considered.size() > 0);
    expectEqual("no promotion on first review", static_cast<int>(rev1.promoted.size()), 0);

    // Second review: still building streak
    WaterReview rev2 = promoter.review(ferry);
    expectEqual("no promotion on second review", static_cast<int>(rev2.promoted.size()), 0);

    // Third review: dwell satisfied, promote
    WaterReview rev3 = promoter.review(ferry);
    expectEqual("promotion on third review", static_cast<int>(rev3.promoted.size()), 1);
}

void testHoldPreventsImmediateDemotion() {
    std::printf("\n   hold prevents immediate demotion\n");

    Ship ferry = ferryAfloat();
    floodCompartment(ferry, 0, 5.0);
    setRollRate(ferry, 0.10);

    WaterCriterion crit;
    crit.dwell = 1;  // Promote immediately
    crit.hold = 3;   // Must fail 3 times to demote
    WaterPromoter promoter(crit);

    // Promote
    promoter.review(ferry);
    expectEqual("promoted after dwell", static_cast<int>(promoter.active().size()), 1);

    // Stop motion: below hold thresholds
    setRollRate(ferry, 0.01);  // Below rollRateHold = 0.03
    setLateralAccel(ferry, 0.5); // Below accelHold = 1.0

    // First review below threshold: still active
    promoter.review(ferry);
    expectEqual("still active after first idle", static_cast<int>(promoter.active().size()), 1);

    // Second review: still active
    promoter.review(ferry);
    expectEqual("still active after second idle", static_cast<int>(promoter.active().size()), 1);

    // Third review: hold satisfied, demote
    WaterReview rev3 = promoter.review(ferry);
    expectEqual("demoted after hold", static_cast<int>(promoter.active().size()), 0);
    expectEqual("demotion recorded", static_cast<int>(rev3.demoted.size()), 1);
}

void testAntiChatterNegativeControl() {
    std::printf("\n   anti-chatter negative control\n");

    Ship ferry = ferryAfloat();
    floodCompartment(ferry, 0, 5.0);

    // Oscillate at threshold with hysteresis disabled
    WaterCriterion noDwell;
    noDwell.dwell = 1;  // No dwell
    noDwell.hold = 1;   // No hold
    WaterPromoter promoter(noDwell);

    int promotions = 0;
    int demotions = 0;

    for (int i = 0; i < 10; ++i) {
        // Alternate above/below threshold
        if (i % 2 == 0) {
            setRollRate(ferry, 0.10);  // Above promote
        } else {
            setRollRate(ferry, 0.01);  // Below hold
        }

        WaterReview rev = promoter.review(ferry);
        promotions += static_cast<int>(rev.promoted.size());
        demotions += static_cast<int>(rev.demoted.size());
    }

    // Vacuity guard: must have chattered
    expectTrue("chatters without hysteresis", promotions >= 4);
    expectTrue("demotions also occur", demotions >= 4);

    std::printf("      (saw %d promotions, %d demotions in 10 reviews)\n", promotions, demotions);
}

void testAntiChatterWithHysteresis() {
    std::printf("\n   anti-chatter with hysteresis\n");

    Ship ferry = ferryAfloat();
    floodCompartment(ferry, 0, 5.0);

    // Same oscillation with hysteresis enabled
    WaterCriterion withHysteresis;
    withHysteresis.dwell = 2;
    withHysteresis.hold = 2;
    WaterPromoter promoter(withHysteresis);

    int promotions = 0;
    int demotions = 0;

    for (int i = 0; i < 10; ++i) {
        if (i % 2 == 0) {
            setRollRate(ferry, 0.10);
        } else {
            setRollRate(ferry, 0.01);
        }

        WaterReview rev = promoter.review(ferry);
        promotions += static_cast<int>(rev.promoted.size());
        demotions += static_cast<int>(rev.demoted.size());
    }

    // Exactly zero, not "at most one". A signal that alternates every review
    // can never build a dwell streak of 2: the qualifying review sets the
    // streak to 1, and the review after it scores 0, which drops the
    // compartment out of `qualifying_` and resets the streak rather than
    // holding it. So the compartment never promotes, and having never
    // promoted it can never demote either.
    //
    // Asserted at 0 rather than <= 1 because <= 1 is also what a promoter
    // that never qualifies *anything* would return -- which is the failure
    // the control above exists to rule out, and the tolerance should not
    // quietly readmit it.
    expectEqual("hysteresis prevents chatter", promotions, 0);
    expectEqual("demotions also suppressed", demotions, 0);

    std::printf("      (saw %d promotions, %d demotions in 10 reviews)\n", promotions, demotions);
}

// ===========================================================================
// Section 3: Budget Enforcement
// ===========================================================================

void testParticleBudgetEnforcement() {
    std::printf("\n   particle budget enforcement\n");

    Ship ferry = ferryAfloat();

    // Flood two compartments, each would use ~5000 particles
    floodCompartment(ferry, 0, 5.0);
    floodCompartment(ferry, 1, 5.0);

    setRollRate(ferry, 0.10);

    WaterCriterion crit;
    crit.dwell = 1;

    // Room for exactly one, derived from the estimator for the same reason the
    // tile budget below is: a constant here is a guess about the estimator's
    // internals, and it is the wrong side of the boundary that passes silently.
    int oneParticles = 0, oneTiles = 0;
    estimateFlipCost(ferry.compartments[0], crit, oneParticles, oneTiles);
    expectTrue("one compartment needs particles", oneParticles > 0);
    crit.particleBudget = 2 * oneParticles - 1;  // fits one, one short of two

    WaterPromoter promoter(crit);

    WaterReview rev = promoter.review(ferry);

    // Both should qualify
    expectTrue("both compartments qualify", rev.considered.size() >= 2);
    expectTrue("first has positive score", rev.considered[0].score > 0.0);
    expectTrue("second has positive score", rev.considered[1].score > 0.0);

    // But only one promoted due to budget
    expectEqual("only one promoted", static_cast<int>(rev.promoted.size()), 1);
    expectTrue("budget exhausted message", rev.problems.size() > 0);
}

void testTileBudgetEnforcement() {
    std::printf("\n   tile budget enforcement\n");

    Ship ferry = ferryAfloat();

    // Flood two compartments
    floodCompartment(ferry, 0, 5.0);
    floodCompartment(ferry, 1, 5.0);

    setRollRate(ferry, 0.10);

    WaterCriterion crit;
    crit.dwell = 1;
    crit.particleBudget = 100000;  // Ample particle budget

    // Size the tile budget off the estimator rather than a guessed constant:
    // room for exactly one of these compartments and not two. A hardcoded
    // number here silently became "room for none" when the estimator's cell
    // size changed, and the test then asserted the wrong side of the boundary.
    int oneParticles = 0, oneTiles = 0;
    estimateFlipCost(ferry.compartments[0], crit, oneParticles, oneTiles);
    expectTrue("one compartment needs tiles", oneTiles > 0);
    crit.tileBudget = 2 * oneTiles - 1;  // fits one, one short of two

    WaterPromoter promoter(crit);

    WaterReview rev = promoter.review(ferry);

    // Both should qualify but budget blocks one
    expectTrue("both compartments qualify", rev.considered.size() >= 2);
    expectEqual("only one promoted", static_cast<int>(rev.promoted.size()), 1);
    expectTrue("budget exhausted message", rev.problems.size() > 0);
}

void testBudgetAccountingAcrossReviews() {
    std::printf("\n   budget accounting across reviews\n");

    Ship ferry = ferryAfloat();
    floodCompartment(ferry, 0, 3.0);
    floodCompartment(ferry, 1, 3.0);
    floodCompartment(ferry, 2, 3.0);

    setRollRate(ferry, 0.10);

    WaterCriterion crit;
    crit.dwell = 1;

    // Room for exactly two of these three, derived rather than guessed.
    int oneParticles = 0, oneTiles = 0;
    estimateFlipCost(ferry.compartments[0], crit, oneParticles, oneTiles);
    expectTrue("one compartment needs particles", oneParticles > 0);
    crit.particleBudget = 3 * oneParticles - 1;  // fits two, one short of three

    WaterPromoter promoter(crit);

    // First review: promote two
    WaterReview rev1 = promoter.review(ferry);
    expectEqual("two promoted initially", static_cast<int>(rev1.promoted.size()), 2);

    int totalParticles1 = rev1.particlesActive;
    expectTrue("particles tracked", totalParticles1 > 0);
    expectTrue("within budget", totalParticles1 <= crit.particleBudget);

    // Second review: still at limit, third cannot promote
    WaterReview rev2 = promoter.review(ferry);
    expectEqual("no new promotions", static_cast<int>(rev2.promoted.size()), 0);
    expectEqual("two still active", static_cast<int>(promoter.active().size()), 2);

    // Demote one by removing motion
    setRollRate(ferry, 0.0);
    for (int i = 0; i < crit.hold; ++i) {
        promoter.review(ferry);
    }
    expectEqual("one demoted", static_cast<int>(promoter.active().size()), 0);

    // Now with budget freed, motion returns and third can promote
    setRollRate(ferry, 0.10);
    WaterReview rev3 = promoter.review(ferry);
    expectTrue("new promotions after budget freed", rev3.promoted.size() > 0);
}

// ===========================================================================
// Section 4: Multiple Compartments
// ===========================================================================

void testMultipleCompartmentsRankedByScore() {
    std::printf("\n   multiple compartments ranked by score\n");

    Ship ferry = ferryAfloat();
    floodCompartment(ferry, 0, 5.0);
    floodCompartment(ferry, 1, 5.0);
    floodCompartment(ferry, 2, 5.0);

    // Set different motion levels
    setRollRate(ferry, 0.15);    // 3x threshold
    setLateralAccel(ferry, 3.0); // 1.5x threshold

    WaterCriterion crit;
    WaterPromoter promoter(crit);

    WaterReview rev = promoter.review(ferry);

    // All three should qualify with same score (weaker trigger)
    expectTrue("three candidates", rev.considered.size() >= 3);
    expectTrue("first has score", rev.considered[0].score > 0.0);
    expectNear("score is weaker trigger", rev.considered[0].score, 1.5, 1e-6);

    // All should have same score since ship motion is uniform
    for (size_t i = 1; i < rev.considered.size(); ++i) {
        if (rev.considered[i].score <= 0) break;  // Only check qualifying candidates
        expectNear("all scores equal", rev.considered[i].score, rev.considered[0].score, 1e-6);
    }
}

void testDifferentCompartmentsSeparateHysteresis() {
    std::printf("\n   different compartments have separate hysteresis\n");

    Ship ferry = ferryAfloat();
    floodCompartment(ferry, 0, 5.0);
    floodCompartment(ferry, 1, 5.0);

    setRollRate(ferry, 0.10);

    WaterCriterion crit;
    crit.dwell = 2;
    WaterPromoter promoter(crit);

    // First review: both start dwell
    promoter.review(ferry);
    expectEqual("no promotions yet", static_cast<int>(promoter.active().size()), 0);

    // Flood a third compartment mid-stream
    floodCompartment(ferry, 2, 5.0);

    // Second review: first two promote, third starts dwell
    WaterReview rev2 = promoter.review(ferry);
    expectEqual("two promoted", static_cast<int>(rev2.promoted.size()), 2);

    // Third review: third now promotes
    WaterReview rev3 = promoter.review(ferry);
    expectEqual("third promoted", static_cast<int>(rev3.promoted.size()), 1);
    expectEqual("three total active", static_cast<int>(promoter.active().size()), 3);
}

// ===========================================================================
// Section 5: Cost and Performance Tracking
// ===========================================================================

void testCostReporting() {
    std::printf("\n   cost reporting\n");

    Ship ferry = ferryAfloat();
    floodCompartment(ferry, 0, 5.0);
    floodCompartment(ferry, 1, 3.0);

    setRollRate(ferry, 0.10);

    WaterCriterion crit;
    crit.dwell = 1;
    crit.coreSecondsPerCompartment = 5.0;
    WaterPromoter promoter(crit);

    WaterReview rev = promoter.review(ferry);

    expectEqual("two promoted", static_cast<int>(rev.promoted.size()), 2);
    expectNear("cost is sum", rev.costActive, 10.0, 1e-6);
    expectTrue("review time measured", rev.microseconds > 0.0);
}

void testCounters() {
    std::printf("\n   promoter counters\n");

    Ship ferry = ferryAfloat();
    floodCompartment(ferry, 0, 5.0);

    setRollRate(ferry, 0.10);

    WaterCriterion crit;
    crit.dwell = 1;
    crit.hold = 2;
    WaterPromoter promoter(crit);

    expectEqual("initial reviews", promoter.reviews(), 0);
    expectEqual("initial promotions", promoter.promotions(), 0);
    expectEqual("initial demotions", promoter.demotions(), 0);

    // Promote
    promoter.review(ferry);
    expectEqual("reviews incremented", promoter.reviews(), 1);
    expectEqual("promotions incremented", promoter.promotions(), 1);

    // Demote
    setRollRate(ferry, 0.0);
    promoter.review(ferry);
    promoter.review(ferry);
    expectEqual("reviews incremented", promoter.reviews(), 3);
    expectEqual("demotions incremented", promoter.demotions(), 1);
}

void testClearResetsState() {
    std::printf("\n   clear resets state\n");

    Ship ferry = ferryAfloat();
    floodCompartment(ferry, 0, 5.0);
    setRollRate(ferry, 0.10);

    WaterCriterion crit;
    crit.dwell = 1;
    WaterPromoter promoter(crit);

    // Promote something
    promoter.review(ferry);
    expectEqual("active before clear", static_cast<int>(promoter.active().size()), 1);
    expectTrue("reviews > 0", promoter.reviews() > 0);

    // Clear
    promoter.clear();
    expectEqual("active after clear", static_cast<int>(promoter.active().size()), 0);
    expectEqual("reviews reset", promoter.reviews(), 0);
    expectEqual("promotions reset", promoter.promotions(), 0);
    expectEqual("demotions reset", promoter.demotions(), 0);
}

// ===========================================================================
// Main
// ===========================================================================

void runWaterPromotionTests() {
    std::printf("\n--- water promotion ---\n");

    // Section 1: Criterion
    testShipAtRestPromotesNothing();
    testShallowPuddlePromotesNothing();
    testMotionBelowThresholdPromotesNothing();
    testRollRateAloneQualifies();
    testLateralAccelAloneQualifies();
    testBothTriggersScoreIsWeaker();

    // Section 2: Hysteresis
    testDwellPreventsImmediatePromotion();
    testHoldPreventsImmediateDemotion();
    testAntiChatterNegativeControl();
    testAntiChatterWithHysteresis();

    // Section 3: Budget
    testParticleBudgetEnforcement();
    testTileBudgetEnforcement();
    testBudgetAccountingAcrossReviews();

    // Section 4: Multiple Compartments
    testMultipleCompartmentsRankedByScore();
    testDifferentCompartmentsSeparateHysteresis();

    // Section 5: Cost and Performance
    testCostReporting();
    testCounters();
    testClearResetsState();
}

