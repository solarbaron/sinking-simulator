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

// Lateral acceleration is a *difference* between two reviews, not a state the
// ship carries, so producing one takes two samples and the time between them.
// This drives a promoter through the first review at rest and returns the dt the
// caller should pass to the second, at which point the velocity below has
// arrived and the acceleration is `mPerSec2` exactly.
//
// **The helper this replaced set `state.velocity = {0, a, 0}` and called it an
// acceleration**, which is how every test in this file passed against an
// implementation that compared a speed in m/s against a threshold in m/s^2.
// A test that shares the code's misconception cannot see it: the units error was
// in the harness and in the thing being tested, agreeing with each other.
constexpr double kAccelDt = 0.5;   // s between the two reviews

void primeAccel(WaterPromoter& promoter, Ship& ship) {
    ship.state.velocity = Vec3{};
    promoter.review(ship, kAccelDt);
}

void setLateralAccel(Ship& ship, double mPerSec2) {
    // Velocity after kAccelDt of that acceleration, from rest.
    ship.state.velocity = Vec3{0, mPerSec2 * kAccelDt, 0};
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

// Roll is rotation about the ship's own axis, and the criterion is read at
// attitudes where that axis is nowhere near the world's. A ship heeled 90
// degrees has her bow axis along world x still -- heel rotates *about* it -- so
// the case that separates the two frames is **pitch**: at 90 degrees of trim the
// bow points at the sky, and a rotation about world x is then pure yaw as the
// ship experiences it, with no roll in it at all.
//
// Against the implementation this replaced -- `abs(angularVelocity.x)` -- this
// reads a full 0.2 rad/s of roll rate on a ship that is not rolling.
void testRollRateIsBodyFrame() {
    std::printf("\n   roll rate is read in the body frame\n");

    Ship ferry = ferryAfloat();
    floodCompartment(ferry, 0, 5.0);

    // Bow up: rotate 90 degrees about the transverse (y) axis.
    ferry.state.orientation = Quat::fromAxisAngle(Vec3{0, 1, 0}, kPi / 2.0);
    ferry.state.angularVelocity = Vec3{0.2, 0, 0};   // world x, 4x the threshold

    WaterCriterion crit;
    WaterPromoter promoter(crit);
    WaterReview rev = promoter.review(ferry);

    expectTrue("the ship is considered", rev.considered.size() > 0);
    // Her bow axis is world +z here, and the angular velocity is along world x,
    // so the roll component is exactly zero.
    expectNear("rotation about world x is not roll when the bow points up",
               rev.considered[0].rollRate, 0.0, 1e-12);
    expectTrue("so nothing qualifies on it", rev.considered[0].score == 0.0);

    // The control: the *same* rate about the axis that really is her bow reads
    // in full. Without this the test above would pass on an implementation that
    // returned zero unconditionally.
    ferry.state.angularVelocity = Vec3{0, 0, 0.2};   // world z == her bow now
    WaterPromoter upright(crit);
    WaterReview rolling = upright.review(ferry);
    expectNear("but rotation about her own bow axis is",
               rolling.considered[0].rollRate, 0.2, 1e-12);
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

    WaterCriterion crit;
    WaterPromoter promoter(crit);
    primeAccel(promoter, ferry);
    setLateralAccel(ferry, 3.0); // 1.5x accelPromote

    WaterReview rev = promoter.review(ferry, kAccelDt);

    expectTrue("lateral accel alone produces candidates", rev.considered.size() > 0);
    expectTrue("with positive score", rev.considered[0].score > 0.0);
    expectNear("score matches accel ratio", rev.considered[0].score, 3.0 / 2.0, 1e-6);

    // **The negative control for the units, and it is the whole point of the
    // rewrite.** The same ship reviewed with no time between samples has no
    // acceleration, so a criterion reading a genuine difference scores 0 while
    // one reading `length(velocity)` scores 1.5/... whatever dt it was given.
    // Against the implementation this replaced, this assertion fails.
    WaterPromoter noTime(crit);
    noTime.review(ferry);                       // seeds, dt defaults to 0
    WaterReview flat = noTime.review(ferry, 0); // no elapsed time
    expectTrue("with no elapsed time there is no acceleration to read",
               flat.considered[0].score == 0.0);
}

void testBothTriggersScoreIsWeaker() {
    std::printf("\n   both triggers: score is the weaker\n");

    Ship ferry = ferryAfloat();
    floodCompartment(ferry, 0, 5.0);

    // Both above threshold, roll stronger
    setRollRate(ferry, 0.15);   // 3x rollRatePromote

    WaterCriterion crit;
    WaterPromoter promoter(crit);
    primeAccel(promoter, ferry);
    setLateralAccel(ferry, 4.0); // 2x accelPromote

    WaterReview rev = promoter.review(ferry, kAccelDt);

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

// The two budgets state the same volume, and that is an invariant rather than a
// coincidence: `estimateFlipCost` derives both from the wetted volume, so a
// budget pair that disagrees makes the larger of the two unreachable. It
// disagreed by 6.25x once -- 2000 tiles admitting 16 m3 against 100 000
// particles admitting 100 -- and the consequence was that no compartment over
// 16 m3 could be promoted at all, on a ship whose compartments run to 1232.
//
// Asserted against the arithmetic rather than against the constants, so it fails
// when someone edits one number and not the other, which is exactly how it broke.
// A copied ship has no promoted water, and its promoter must agree.
//
// `rightingArmAtHeel`, `diagnostics` and `girder.cpp` all take a ship copy to ask
// a hydrostatic question, and none of them wants a live particle field -- they
// want the scalar `waterVolume` every one of their readers understands. So the
// copy drops `activeWaterFields_`, which was already true, and it must clear
// `waterPromoter_` **with** it: copying the promoter while dropping the fields
// leaves a copy whose promoter believes compartments are active and whose map is
// empty, and the next review on that copy would demote compartments whose state
// had never been transferred anywhere.
//
// The two halves are one invariant and this asserts them together, because the
// version that copied the promoter and dropped the fields passed every test in
// this file.
// `minDepth` is declared and deliberately not enforced, and this is the
// assertion that says so out loud.
//
// The only depth available is the mean over the compartment's whole bounding-box
// footprint, which on a ship whose compartments are tens of metres a side reads
// ~0.02 m for water that stands well over half a metre against a bulkhead under
// heel. Gating on it filtered *every* compartment aboard and returned an empty
// candidate list from every review -- so the field stays as documentation of the
// threshold the real free-surface calculation will be compared against, and
// `minVolume` is the only geometric guard.
//
// **An unenforced parameter that silently starts being enforced is a trap**, and
// this one is primed: the header invites re-enabling it. This test names the
// current state so the change is visible rather than mysterious.
//
// **What it cannot do is be the thing that catches it.** Enforcing `minDepth` was
// tried as a deliberate mutation, and the suite did not report a failure -- it
// *hung*, output stopping mid-line in a later suite, which on a first reading
// looked like a clean pass with zero failures. That is the characteristic kill in
// this codebase and the reason a mutation harness here needs a time bound rather
// than a pass/fail read; see CLAUDE.md. So the assertions below document the
// state, and the hang is what would actually announce the change. Both are worth
// having and neither substitutes for the other.
void testMinDepthIsNotEnforced() {
    std::printf("\n   minDepth is declared and deliberately not enforced\n");

    Ship ferry = ferryAfloat();
    floodCompartment(ferry, 0, 5.0);
    setRollRate(ferry, 0.15);

    WaterCriterion crit;
    WaterPromoter promoter(crit);
    const WaterReview rev = promoter.review(ferry);

    expectTrue("the compartment is considered", rev.considered.size() > 0);

    // The guard that makes this non-vacuous: the depth really is far below the
    // threshold, so a criterion that enforced it would certainly reject this.
    expectTrue("and its mean depth is well under minDepth",
               rev.considered[0].depth < crit.minDepth);
    std::printf("      (depth %.4f m against minDepth %.2f m)\n",
                rev.considered[0].depth, crit.minDepth);

    // ...and it qualifies anyway.
    expectTrue("yet it qualifies, because only minVolume gates promotion",
               rev.considered[0].score > 0.0);
}

void testCopiedShipHasNoPromotedWater() {
    std::printf("\n   a copied ship has no promoted water, promoter included\n");

    Ship ferry = ferryAfloat();
    floodCompartment(ferry, 0, 5.0);

    // Drive the original to a genuine promotion so there is something to lose.
    WaterCriterion crit;
    WaterPromoter promoter(crit);
    setRollRate(ferry, 0.15);
    promoter.review(ferry);
    const WaterReview rev = promoter.review(ferry);

    expectTrue("the original promoted something to copy away",
               promoter.active().size() > 0);
    expectEqual("and it promoted this review", static_cast<int>(rev.promoted.size()),
                static_cast<int>(promoter.active().size()));

    // The ship's own promoter is the one the copy carries; drive it too.
    ferry.waterPromoter_.review(ferry);
    ferry.waterPromoter_.review(ferry);
    expectTrue("the ship's own promoter has active compartments",
               ferry.waterPromoter_.active().size() > 0);

    const Ship copy = ferry;

    expectEqual("the copy has no FLIP fields",
                static_cast<int>(copy.activeWaterFields_.size()), 0);
    expectEqual("and its promoter agrees that nothing is active",
                static_cast<int>(copy.waterPromoter_.active().size()), 0);
    expectEqual("and its review count is reset, not inherited",
                copy.waterPromoter_.reviews(), 0);

    // The control: the *original* still has its promotions. Without this the
    // assertions above would pass on a copy constructor that cleared the source.
    expectTrue("while the original still has its own",
               ferry.waterPromoter_.active().size() > 0);
}

void testBudgetsAdmitTheSameVolume() {
    std::printf("\n   the two budgets admit the same volume\n");

    WaterCriterion crit;

    // What each budget admits, at the cost model the promoter actually uses.
    const double h = 0.05;                       // m, estimateFlipCost's cell size
    const double particlesPerM3 = 1000.0;
    const double tilesPerM3 = 1.0 / (64.0 * h * h * h);   // 4x4x4 cells per tile

    const double m3FromParticles = crit.particleBudget / particlesPerM3;
    const double m3FromTiles = crit.tileBudget / tilesPerM3;

    std::printf("      particles admit %.1f m3, tiles admit %.1f m3\n",
                m3FromParticles, m3FromTiles);

    expectNear("the particle and tile budgets admit the same volume",
               m3FromParticles, m3FromTiles, 1e-9);

    // And the guard against a vacuous version of the test above: if either
    // budget were zero they would agree trivially and admit nothing.
    expectTrue("and it is not zero", m3FromParticles > 1.0);

    // The per-compartment ceiling is the one the decision is really made on, and
    // it must not exceed what the shared budget can hold -- a compartment that
    // passes the volume gate and is then refused by the budget is the silent
    // rejection this pair replaced.
    expectTrue("the per-compartment ceiling fits inside the shared budget",
               crit.maxVolumePerCompartment <= m3FromParticles);
}

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
// Section 4: State Transfer
// ===========================================================================

void testRoundTripMassConservation() {
    std::printf("\n   round-trip mass conservation\n");

    Ship ferry = ferryAfloat();
    Compartment& forepeak = ferry.compartments[0];

    // Set known water volume
    const double initialVolume = 5.0;  // m³
    forepeak.waterVolume = initialVolume;

    WaterCriterion crit;

    // Promote to FLIP
    auto field = promoteWater(forepeak, ferry, crit);
    expectTrue("field created", field != nullptr);

    // Demote back to compartment
    demoteWater(forepeak, std::move(field));

    // Mass conservation: exact, not epsilon-based. The round trip is
    // Compartment::waterVolume → mass = volume * density → particles with
    // setTotalMass() → totalMass() → volume = mass / density. The density
    // cancels, setTotalMass() is exact (§5), and totalMass() is compensated,
    // so the volume returned is the volume sent, bit for bit.
    expectNear("waterVolume conserved exactly", forepeak.waterVolume, initialVolume, 0.0);
}

void testEmptyCompartmentReturnsNull() {
    std::printf("\n   empty compartment returns nullptr\n");

    Ship ferry = ferryAfloat();
    Compartment& forepeak = ferry.compartments[0];

    // Empty compartment
    forepeak.waterVolume = 0.0;

    WaterCriterion crit;

    auto field = promoteWater(forepeak, ferry, crit);
    expectTrue("empty compartment returns nullptr", field == nullptr);
}

void testCentroidPreservation() {
    std::printf("\n   centroid preservation\n");

    Ship ferry = ferryAfloat();
    Compartment& forepeak = ferry.compartments[0];

    forepeak.waterVolume = 5.0;

    WaterCriterion crit;

    auto field = promoteWater(forepeak, ferry, crit);
    expectTrue("field created", field != nullptr);

    // Demote and check centroid is reasonable
    demoteWater(forepeak, std::move(field));

    // Centroid should be within compartment bounds
    expectTrue("centroid x in bounds",
               forepeak.waterCentroid.x >= forepeak.bboxLo.x &&
               forepeak.waterCentroid.x <= forepeak.bboxHi.x);
    expectTrue("centroid y in bounds",
               forepeak.waterCentroid.y >= forepeak.bboxLo.y &&
               forepeak.waterCentroid.y <= forepeak.bboxHi.y);
    expectTrue("centroid z in bounds",
               forepeak.waterCentroid.z >= forepeak.bboxLo.z &&
               forepeak.waterCentroid.z <= forepeak.bboxHi.z);

    std::printf("      (centroid at [%.3f, %.3f, %.3f] m)\n",
                forepeak.waterCentroid.x, forepeak.waterCentroid.y, forepeak.waterCentroid.z);
}

void testSurfaceOffsetMatchesVolume() {
    std::printf("\n   surface offset matches volume\n");

    Ship ferry = ferryAfloat();
    Compartment& forepeak = ferry.compartments[0];

    const double initialVolume = 5.0;
    forepeak.waterVolume = initialVolume;

    WaterCriterion crit;

    auto field = promoteWater(forepeak, ferry, crit);
    expectTrue("field created", field != nullptr);

    // Record the grid plan area before demoting
    const double gridPlanArea = field->grid.planArea();

    // Demote and check surface offset
    demoteWater(forepeak, std::move(field));

    // Surface offset should equal volume / gridPlanArea, where gridPlanArea is
    // discretized to whole cells. The grid is constructed with ceil() so its
    // area can differ from the bbox by up to h^2 per axis.
    const double expectedOffset = forepeak.waterVolume / gridPlanArea;

    // Tolerance accounts for grid discretization: h = 0.05 m, so one cell edge
    // is a ~0.2% change in a ~25 m² compartment footprint, giving ~0.2% depth
    // error. Use 1e-6 as conservative bound on the arithmetic error alone.
    expectNear("surfaceOffset matches volume/planArea",
               forepeak.surfaceOffset, expectedOffset, 1e-6);

    std::printf("      (surfaceOffset %.6f m, grid planArea %.3f m²)\n",
                forepeak.surfaceOffset, gridPlanArea);
}

// ===========================================================================
// Section 5: Multiple Compartments
// ===========================================================================

void testMultipleCompartmentsRankedByScore() {
    std::printf("\n   multiple compartments ranked by score\n");

    Ship ferry = ferryAfloat();
    floodCompartment(ferry, 0, 5.0);
    floodCompartment(ferry, 1, 5.0);
    floodCompartment(ferry, 2, 5.0);

    // Set different motion levels
    setRollRate(ferry, 0.15);    // 3x threshold

    WaterCriterion crit;
    WaterPromoter promoter(crit);
    // The acceleration is a difference between reviews, so it takes two of them
    // and the dt between. Reviewing once leaves it at zero and the score falls
    // back to the roll ratio -- which is a correct answer to a different
    // question, and would make the 1.5 below unreachable.
    primeAccel(promoter, ferry);
    setLateralAccel(ferry, 3.0); // 1.5x threshold

    WaterReview rev = promoter.review(ferry, kAccelDt);

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
// Section 6: Cost and Performance Tracking
// ===========================================================================

void testCostReporting() {
    std::printf("\n   cost reporting\n");

    Ship ferry = ferryAfloat();
    floodCompartment(ferry, 0, 5.0);
    floodCompartment(ferry, 1, 3.0);

    setRollRate(ferry, 0.10);

    WaterCriterion crit;
    crit.dwell = 1;
    WaterPromoter promoter(crit);

    WaterReview rev = promoter.review(ferry);

    expectEqual("two promoted", static_cast<int>(rev.promoted.size()), 2);
    expectNear("cost is the sum over active compartments",
               rev.costActive, 2.0 * crit.coreSecondsPerCompartment, 1e-6);
    expectTrue("review time measured", rev.microseconds > 0.0);

    // **What this test does not assert, stated so the gap is visible.**
    //
    // The two compartments hold 5 m³ and 3 m³ and are billed identically,
    // because `cand.cost = criterion.coreSecondsPerCompartment` is a
    // per-compartment constant. `water_probe --cost` measured the truth as
    // 27.9 core-s/sim-s at 1 m³ and 3030 at 100 -- it scales with the water,
    // steeply -- so the cost model is not merely mis-valued but the wrong shape,
    // and this test asserts additivity over a term that should not be uniform.
    //
    // It is left asserting what the code does rather than what is true, because
    // changing `cost` to scale with volume is a change to the promoter's
    // behaviour and not to a test. The earlier version of this test pinned
    // `coreSecondsPerCompartment = 5.0` -- the refuted value -- which made it
    // pass regardless of the field it was ostensibly testing; it now reads the
    // criterion's own default so that a change to that default is at least
    // visible here.
    expectNear("the two sizes are billed the same, which is the known defect",
               rev.promoted[0].cost, rev.promoted[1].cost, 0.0);
    expectTrue("and the compartments really are different sizes, or the line above proves nothing",
               std::abs(rev.promoted[0].volume - rev.promoted[1].volume) > 1.0);
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
    testRollRateIsBodyFrame();
    testRollRateAloneQualifies();
    testLateralAccelAloneQualifies();
    testBothTriggersScoreIsWeaker();

    // Section 2: Hysteresis
    testDwellPreventsImmediatePromotion();
    testHoldPreventsImmediateDemotion();
    testAntiChatterNegativeControl();
    testAntiChatterWithHysteresis();

    // Section 3: Budget
    testMinDepthIsNotEnforced();
    testCopiedShipHasNoPromotedWater();
    testBudgetsAdmitTheSameVolume();
    testParticleBudgetEnforcement();
    testTileBudgetEnforcement();
    testBudgetAccountingAcrossReviews();

    // Section 4: State Transfer
    testRoundTripMassConservation();
    testEmptyCompartmentReturnsNull();
    testCentroidPreservation();
    testSurfaceOffsetMatchesVolume();

    // Section 5: Multiple Compartments
    testMultipleCompartmentsRankedByScore();
    testDifferentCompartmentsSeparateHysteresis();

    // Section 6: Cost and Performance
    testCostReporting();
    testCounters();
    testClearResetsState();
}

