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
#include "engine/sim/flip.hpp"
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
// A compartment that is pumped out is demoted immediately, without waiting for
// the dwell.
//
// `!volumeHolds` is the second of the two demotion conditions and no test had
// ever taken it: every demotion in the suite came from `idleReviews >= hold`,
// i.e. the ship going quiet. Water leaving the compartment is the other way a
// promotion ends, and it is the one a bilge pump causes.
//
// The distinction matters because the two have different timing. Going quiet is
// deliberately slow -- `hold` reviews of hysteresis, so a ship that rolls once
// more does not thrash. Losing the water is immediate: there is nothing left to
// resolve, and holding a FLIP field open for three more reviews over a compartment
// with 0.1 m³ in it is pure cost. This asserts the immediacy, not just the fact.
void testDrainingACompartmentDemotesItAtOnce() {
    std::printf("\n   a compartment that drains is demoted at once, not after the hold\n");

    Ship ferry = ferryAfloat();
    floodCompartment(ferry, 0, 5.0);
    setRollRate(ferry, 0.15);          // well above promote and hold alike

    WaterCriterion crit;
    crit.dwell = 1;                    // promote on the first qualifying review
    WaterPromoter promoter(crit);

    promoter.review(ferry);
    expectEqual("the compartment is active to begin with",
                static_cast<int>(promoter.active().size()), 1);

    // Pump it out, and keep the ship rolling hard so the *motion* still holds --
    // otherwise this would demote through the idle path and prove nothing.
    floodCompartment(ferry, 0, 0.1);   // below minVolume
    const WaterReview rev = promoter.review(ferry);

    expectEqual("draining it demotes it on the very next review",
                static_cast<int>(rev.demoted.size()), 1);
    expectEqual("and it is no longer active",
                static_cast<int>(promoter.active().size()), 0);

    // The control that makes the timing claim mean something: the ship is still
    // rolling above the hold threshold, so the idle path cannot be what demoted
    // it, and `hold` is 3 -- more reviews than have happened.
    expectTrue("while the motion was still above the hold threshold",
               std::abs(ferry.state.angularVelocity.x) >= crit.rollRateHold);
    expectTrue("and fewer reviews have passed than the hold would need",
               promoter.reviews() < crit.hold);
}

// Lateral acceleration alone can hold a compartment active.
//
// `accelHolds` is read on every review and had never been the deciding term:
// every hold and anti-chatter test in this file drives roll only, so `accel` was
// 0 throughout and the branch was reached but never *taken*. A criterion that had
// dropped the acceleration term from the hold test entirely would pass all of
// them.
void testAccelerationAloneHoldsACompartment() {
    std::printf("\n   lateral acceleration alone can hold a compartment active\n");

    Ship ferry = ferryAfloat();
    floodCompartment(ferry, 0, 5.0);

    WaterCriterion crit;
    crit.dwell = 1;
    WaterPromoter promoter(crit);

    // Promote on roll.
    setRollRate(ferry, 0.15);
    primeAccel(promoter, ferry);
    promoter.review(ferry, kAccelDt);
    expectEqual("promoted on roll", static_cast<int>(promoter.active().size()), 1);

    // Now stop rolling, but keep throwing her sideways: roll is below its hold
    // threshold and acceleration is above its own, so the compartment must stay.
    //
    // **A sustained acceleration is a velocity that keeps climbing**, not one
    // that is set and reset. The first version of this test set the velocity and
    // then zeroed it before the next review, so successive samples differed by
    // +v then -v and half the reviews saw the ship decelerating to a stop -- the
    // compartment demoted, and the failure was the test's own, not the code's.
    // That is the same misconception the units fix in this file was about,
    // arriving from the other side.
    setRollRate(ferry, 0.0);
    double vy = 0.0;
    const double accel = 2.0 * crit.accelHold;
    for (int i = 0; i < crit.hold + 2; ++i) {
        vy += accel * kAccelDt;                     // still accelerating
        ferry.state.velocity = Vec3{0, vy, 0};
        promoter.review(ferry, kAccelDt);
    }

    expectEqual("acceleration alone keeps it active past the hold count",
                static_cast<int>(promoter.active().size()), 1);
    expectEqual("with no demotions at all", promoter.demotions(), 0);

    // The control: with neither trigger it goes, and within the hold count. Same
    // ship, same promoter, so the only thing that changed is the acceleration.
    for (int i = 0; i < crit.hold + 1; ++i) {
        ferry.state.velocity = Vec3{};
        promoter.review(ferry, kAccelDt);
    }
    expectEqual("but with the acceleration gone as well, it demotes",
                static_cast<int>(promoter.active().size()), 0);
}

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

    // **The same invariant through copy *assignment*, which is a different
    // function and was covered by nothing.** `Ship::operator=(const Ship&)` is
    // referenced by zero objects in the whole build -- it carries its own copy of
    // the promoter/field-map pairing (a `clear()` beside a delete loop) and could
    // have copied the promoter while dropping the fields, which is precisely the
    // defect the constructor test above exists to prevent, and every test would
    // still have passed. Half the special members were dead code holding an
    // invariant asserted only on the other half.
    Ship assigned = ferryAfloat();
    assigned = ferry;

    expectEqual("the assigned-to ship has no FLIP fields either",
                static_cast<int>(assigned.activeWaterFields_.size()), 0);
    expectEqual("and its promoter agrees that nothing is active",
                static_cast<int>(assigned.waterPromoter_.active().size()), 0);
    expectEqual("and its review count is reset too",
                assigned.waterPromoter_.reviews(), 0);
    expectTrue("and the source is still intact after being assigned from",
               ferry.waterPromoter_.active().size() > 0);

    // **And through the move constructor**, which is the one member that
    // *transfers* the promoter and the field map rather than clearing them --
    // and which was referenced by zero objects anywhere: every candidate site is
    // copy-elided, so the transfer at ship.cpp had never executed once.
    //
    // `std::move` is explicit rather than relying on a return value, because
    // elision is exactly what kept this untested.
    const int activeBeforeMove = static_cast<int>(ferry.waterPromoter_.active().size());
    const int reviewsBeforeMove = ferry.waterPromoter_.reviews();
    expectTrue("the ship about to be moved has something worth transferring",
               activeBeforeMove > 0 && reviewsBeforeMove > 0);

    Ship moved = std::move(ferry);

    expectEqual("a move transfers the promoter's active list rather than clearing it",
                static_cast<int>(moved.waterPromoter_.active().size()), activeBeforeMove);
    expectEqual("and its review count",
                moved.waterPromoter_.reviews(), reviewsBeforeMove);
    // The source is left valid and empty: `activeWaterFields_` was moved from, so
    // the destructor must not double-delete. Asserting the map rather than the
    // promoter, because that is the one holding raw pointers.
    expectEqual("and leaves the source's field map empty",
                static_cast<int>(ferry.activeWaterFields_.size()), 0);  // NOLINT: moved-from, deliberately
}

// **What the estimator bills against what the promoter actually allocates.**
//
// These deliberately differ, and the ratio between them is a documented number
// that nothing checked. `WaterCriterion` says it outright: `estimateFlipCost`
// bills 1000 particles per m3 and `flip::seedBox` at 2^3 per cell on a 0.05 m
// grid puts `8/h^3` = 64 000 there, so the estimate is a *budget denomination*
// that is self-consistent with `tileBudget`, and any sentence about bytes has to
// use the other number.
//
// A prose ratio between two live computations is the shape this repo keeps
// finding wrong -- the zone's elastic and plastic per-element costs sat 5.6x from
// each other behind exactly such a sentence. So it is measured here. The excess
// over 64 000 is the wetted box rounding onto the grid: `promoteWater` takes the
// compartment's full bounding box in x and y and a depth of `volume / planArea`,
// and `seedBox` fills whole cells, so the seeded box is a little larger than the
// water it stands for.
void testTheEstimateAndTheAllocationDifferByTheDocumentedFactor() {
    std::printf("\n   the estimate against the allocation\n");

    // Volumes above the seeding's resolution limit -- below it `promoteWater`
    // refuses and there is no allocation to compare against; see
    // `testTheRoundTripConservesAtEveryVolume`.
    for (double volume : {3.0, 5.0, 12.0}) {
        Ship ferry = ferryAfloat();
        Compartment& forepeak = ferry.compartments[0];
        forepeak.waterVolume = volume;

        WaterCriterion crit;
        int predictedParticles = 0, predictedTiles = 0;
        estimateFlipCost(forepeak, crit, predictedParticles, predictedTiles);

        auto field = promoteWater(forepeak, ferry, crit);
        expectTrue("the compartment promotes", field != nullptr);
        const double actual = static_cast<double>(field->particles.size());
        const double perM3 = actual / volume;

        std::printf("      %5.1f m3: billed %6d, allocated %7.0f -- %6.0f/m3, %.1fx the estimate\n",
                    volume, predictedParticles, actual, perM3,
                    actual / std::max(1.0, static_cast<double>(predictedParticles)));

        // **Not asserted as a rate per cubic metre, because it is not one.**
        // Measured here: 123 947/m3 at 3 m3, 74 368 at 5, 61 973 at 12 -- the
        // count is the same 371 840 across 3, 5 and 8 m3 and then doubles. It
        // quantises to *particle planes*: `seedBox` lays planes every `h/perAxis`
        // = 25 mm and keeps those under the surface, so the count steps as the
        // wetted depth crosses each plane and is flat between them. The nominal
        // `8/h^3` = 64 000 is what a box filling whole cells would give, and no
        // real compartment does. Pinning it would be a figure that holds at the
        // volumes tried and fails at the ones that were not -- the shape this
        // suite has been unpicking elsewhere.
        //
        // What is stable is the order of the gap, which is the documented claim.
        expectTrue("the promoter allocates orders more than the budget bills",
                   actual > 20.0 * predictedParticles);
        expectTrue("and the rate is somewhere near the nominal 8 per cell",
                   perM3 > 40000.0 && perM3 < 200000.0);
        // The estimator's own rate, so this reads the documented 64x rather than
        // re-typing 1000 and comparing the test's copy of the model with itself.
        const double billedPerM3 = predictedParticles / volume;
        expectNear("the estimator bills the 1000 per m3 the budget is denominated in",
                   billedPerM3, 1000.0, 1e-9);

        demoteWater(forepeak, std::move(field));
    }
}

// **Mass across the promote/demote round trip, at volumes that span the seeding's
// own resolution limit rather than at one that clears it.**
//
// The round trip is asserted elsewhere at 5.0 m3 with a tolerance of exactly zero,
// which is the right assertion at a volume chosen before anyone knew the limit
// existed. Below it the water did not come back short -- it came back as *nothing*.
// `seedBox` puts its lowest plane of particles a quarter of a cell above the floor,
// 12.5 mm at h = 0.05, and filters every position against the box top, so a wetted
// depth under that seeds no particle at all; `setTotalMass` then had nothing to
// distribute and `demoteWater` read the total back as zero.
//
// Measured on this forepeak, 232 m2 of plan area, so the limit falls at 2.9 m3:
// 2.8 m3 in gave 0.0000 m3 out and 3.0 m3 round-tripped exactly.
void testTheRoundTripConservesAtEveryVolume() {
    std::printf("\n   the round trip across the seeding's resolution limit\n");

    int promoted = 0, refused = 0;
    for (double v : {0.5, 1.0, 2.0, 2.8, 2.9, 3.0, 5.0, 12.0}) {
        Ship ferry = ferryAfloat();
        Compartment& c = ferry.compartments[0];
        c.waterVolume = v;
        WaterCriterion crit;

        auto field = promoteWater(c, ferry, crit);
        if (!field) {
            // Refused. The compartment keeps its water, which is the whole point:
            // a tier that cannot represent this depth must leave it quiescent.
            ++refused;
            expectNear("a refused compartment keeps every drop", c.waterVolume, v, 0.0);
            continue;
        }
        ++promoted;
        expectTrue("a promoted field is never empty", !field->particles.empty());
        demoteWater(c, std::move(field));
        expectNear("and a promoted one round-trips exactly", c.waterVolume, v, 0.0);
    }

    // **`demoteWater`'s own guard, which the refusal above hides.** With
    // `promoteWater` refusing, no empty field ever reaches `demoteWater` from
    // inside this module -- verified by mutation: deleting its guard leaves the
    // suite green. But it takes ownership of a field it did not create, and
    // reading `totalMass()` off an empty one gives zero and writes that back over
    // whatever the compartment held. So it is handed one directly.
    {
        Ship ferry = ferryAfloat();
        Compartment& c = ferry.compartments[0];
        c.waterVolume = 7.5;
        auto barren = std::make_unique<flip::Field>();
        demoteWater(c, std::move(barren));
        expectNear("a field with no particles does not empty the compartment",
                   c.waterVolume, 7.5, 0.0);
    }

    // Vacuity on both sides: a guard that refused everything would conserve
    // perfectly and promote nothing, and one that refused nothing is the bug.
    expectTrue("some volumes promote", promoted > 0);
    expectTrue("and some are refused, or the limit is not being exercised", refused > 0);
    std::printf("      %d promoted, %d refused, none lost\n", promoted, refused);
}

void testBudgetsAdmitTheSameVolume() {
    std::printf("\n   the two budgets admit the same volume\n");

    WaterCriterion crit;

    // **Asked of `estimateFlipCost` rather than re-typed from it.** The first
    // version of this test wrote out `particlesPerM3 = 1000.0` and `h = 0.05` by
    // hand, which made it blind to exactly the edit it existed to catch: change
    // the density or the cell size in `estimateFlipCost` and the two budgets
    // silently stop agreeing while this test goes on comparing its own copy of
    // the model against itself. It caught edits to the two budget *fields* only.
    //
    // A probe compartment of known volume turns the estimator into the two rates.
    Compartment probe;
    probe.waterVolume = 10.0;                 // m³, large enough that the tile
                                              // count is not floored at 1
    int probeParticles = 0, probeTiles = 0;
    estimateFlipCost(probe, crit, probeParticles, probeTiles);

    const double particlesPerM3 = probeParticles / probe.waterVolume;
    const double tilesPerM3 = probeTiles / probe.waterVolume;

    expectTrue("the estimator reports a non-zero particle rate", particlesPerM3 > 0);
    expectTrue("and a non-zero tile rate", tilesPerM3 > 0);

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

// **A compartment the budget cannot afford keeps its dwell.** That is the
// semantics this tier settled on -- dwell counts qualification, not
// qualification-and-affordability -- and `WaterCriterion::dwell` carries the
// evidence for choosing it over the other reading. This is the test that holds the
// code to it.
//
// Three tied 5 m3 compartments and room for exactly one. Compartment 0 promotes;
// compartments 1 and 2 qualify on every review and are refused on every review, and
// the question is what being refused does to their hysteresis. It used to destroy
// it: the promotion loop stopped at the first compartment the budget could not
// admit, so 1 and 2 were never written back to `qualifying_` and lost their streaks
// on the spot. They were not deferred, they were restarted, and restarted by a
// decision taken about compartment 0.
//
// What that left was an oscillation rather than a standstill, which matters for how
// this test is built. Each of them rebuilt to a streak of two on every even review,
// was refused there and dropped, and stood at one again on every odd review -- so
// whether a starved compartment ever promoted at all came down to the *parity* of
// the review the budget happened to free on. That is not a property of the ship, of
// the water, or of the criterion.
//
// The last third is the negative control, and it is the part that separates the two
// designs rather than merely observing that something eventually promotes. Free the
// budget and ask *when* compartment 1 arrives. With its dwell retained it promotes
// on that very review, having become promotable nine reviews earlier. With the
// streak oscillating it stands at one on the review chosen here and promotes
// nothing. An assertion that it "promotes eventually" passes under both, and so
// does one that frees the budget on an even review -- see the note at the control
// itself.
void testABudgetStarvedCompartmentKeepsItsDwell() {
    std::printf("\n   a budget-starved compartment keeps its dwell\n");

    Ship ferry = ferryAfloat();
    floodCompartment(ferry, 0, 5.0);
    floodCompartment(ferry, 1, 5.0);
    floodCompartment(ferry, 2, 5.0);
    setRollRate(ferry, 0.15);

    WaterCriterion crit;
    crit.dwell = 2;   // the default, written out because the arithmetic below is in it

    // Room for exactly one, sized off the estimator rather than off a constant, for
    // the reason `testParticleBudgetEnforcement` gives: a hardcoded number here is a
    // guess about the estimator's internals and lands on the wrong side of the
    // boundary in silence when they change.
    int oneParticles = 0, oneTiles = 0;
    estimateFlipCost(ferry.compartments[0], crit, oneParticles, oneTiles);
    expectTrue("one compartment needs particles and tiles",
               oneParticles > 0 && oneTiles > 0);
    crit.particleBudget = 2 * oneParticles - 1;   // fits one, one short of two
    crit.tileBudget = 2 * oneTiles - 1;

    WaterPromoter promoter(crit);

    // Read through a lookup rather than by index: before the fix `starved` is empty
    // on every review, and a test that indexed into it would crash instead of
    // reporting, which is the wrong kind of red.
    const auto refusedFor = [](const WaterReview& r, int compartment) {
        for (const WaterStarved& st : r.starved)
            if (st.compartment == compartment) return st.refusedReviews;
        return 0;
    };
    const auto qualifiedFor = [](const WaterReview& r, int compartment) {
        for (const WaterStarved& st : r.starved)
            if (st.compartment == compartment) return st.qualifyingReviews;
        return 0;
    };

    // Review 1: all three qualify and none has served its dwell, so nothing is
    // starved. Being short of dwell is the hysteresis working and must not be
    // reported as a budget refusal -- telling those two apart is the whole point of
    // the channel.
    const WaterReview first = promoter.review(ferry);
    expectEqual("nothing promotes on the first review",
                static_cast<long long>(first.promoted.size()), 0);
    expectEqual("and nothing is starved, because nothing is promotable yet",
                static_cast<long long>(first.starved.size()), 0);

    // Review 2: all three have served their dwell and the budget admits one.
    const WaterReview second = promoter.review(ferry);
    expectEqual("one promotes on the second review",
                static_cast<long long>(second.promoted.size()), 1);
    expectEqual("and it is compartment 0, the tie-break order",
                second.promoted.empty() ? -1 : second.promoted[0].compartment, 0);
    expectEqual("the other two are starved and named",
                static_cast<long long>(second.starved.size()), 2);
    expectEqual("compartment 1 has been refused once", refusedFor(second, 1), 1);
    expectEqual("compartment 2 with it", refusedFor(second, 2), 1);

    // And the candidate list says so on the candidate itself, which is where this
    // tier puts a refusal -- the answer it already gives for a compartment turned
    // away on volume or on grid spacing.
    bool starvedCandidateSaysWhy = false;
    for (const WaterCandidate& c : second.considered)
        if (c.compartment == 1 && c.why.find("budget refused it") != std::string::npos)
            starvedCandidateSaysWhy = true;
    expectTrue("and the starved candidate carries the reason in its own why",
               starvedCandidateSaysWhy);

    // Reviews 3 to 10. They qualify every time and are refused every time, and the
    // refusal count is the measurement: a count stuck at one is exactly what a
    // destroyed streak looks like from outside, which is why it is asserted to
    // climb rather than merely to be positive.
    bool bothNamedEveryReview = true, refusalCountClimbs = true, nothingElsePromoted = true;
    int previousRefused = refusedFor(second, 1);
    int qualifyingAtTen = 0;
    for (int review = 3; review <= 10; ++review) {
        const WaterReview r = promoter.review(ferry);
        if (!r.promoted.empty()) nothingElsePromoted = false;
        if (r.starved.size() != 2 || refusedFor(r, 1) == 0 || refusedFor(r, 2) == 0)
            bothNamedEveryReview = false;
        if (refusedFor(r, 1) != previousRefused + 1) refusalCountClimbs = false;
        previousRefused = refusedFor(r, 1);
        qualifyingAtTen = qualifiedFor(r, 1);
    }
    expectTrue("both are named starved on every review of the eight",
               bothNamedEveryReview);
    expectTrue("nothing else promotes while the budget is full", nothingElsePromoted);
    expectTrue("and the refusal count climbs on every one of them", refusalCountClimbs);
    expectEqual("so by review 10 compartment 1 has been refused nine reviews running",
                previousRefused, 9);
    // The streak underneath it, because the refusal count is derived from it and a
    // reader handed one is entitled to the other.
    expectEqual("on a streak of ten qualifying reviews", qualifyingAtTen, 10);

    // The control. Free the budget the only way a caller can from outside: drain
    // the promoted compartment, which demotes it at once on the volume guard
    // (`testDrainingACompartmentDemotesItAtOnce`). The criterion is fixed for the
    // promoter's lifetime so the budget cannot be raised mid-run -- and draining is
    // the more honest scenario anyway, being what happens on a ship when a
    // compartment is pumped out.
    //
    // **This is review 11, and the odd count is deliberate.** The streak the old
    // code left behind was not stuck at zero, it *oscillated*: compartment 1 built
    // to two on every even review, was refused there and dropped, and was back at
    // one on every odd review. So whether a starved compartment ever promoted
    // depended on the parity of the review the budget happened to free on, and a
    // control run at review 12 would have passed against the defect for that reason
    // alone. Freeing on an odd review is where the two designs disagree: with the
    // streak kept, compartment 1 has eleven reviews of dwell and promotes; with it
    // oscillating it stands at one, short of `dwell = 2`, and promotes nothing.
    floodCompartment(ferry, 0, 0.0);
    const WaterReview freed = promoter.review(ferry);

    expectEqual("draining the promoted compartment demotes it",
                static_cast<long long>(freed.demoted.size()), 1);
    expectEqual("and the starved compartment promotes on that very review, not two later",
                static_cast<long long>(freed.promoted.size()), 1);
    expectEqual("which is compartment 1, the next in the tie-break order",
                freed.promoted.empty() ? -1 : freed.promoted[0].compartment, 1);
    expectEqual("while compartment 2 is starved still, one place behind it",
                refusedFor(freed, 2), 10);

    // And it is not silent about it, which is the standard the gas tier set after
    // the same defect was found there: a budget that refuses a qualifying
    // compartment says so, in a message carrying the word `water_probe` counts.
    bool saidSoOutLoud = false;
    for (const std::string& problem : freed.problems)
        if (problem.find("budget") != std::string::npos &&
            problem.find("dwell is retained") != std::string::npos)
            saidSoOutLoud = true;
    expectTrue("and the review says out loud that one is still being refused",
               saidSoOutLoud);
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

    // **And handing that nullptr straight to `demoteWater` must be a no-op**,
    // which is the other half of the same contract and was covered by nothing:
    // the guard `if (!field) return;` sat one line away from the test that
    // already had a `nullptr` in hand and threw it away. A promoter that
    // demotes a compartment it never managed to promote takes exactly this
    // path, and without the guard it dereferences null.
    forepeak.waterVolume = 7.25;      // a value that must survive untouched
    forepeak.surfaceOffset = 1.5;
    demoteWater(forepeak, nullptr);

    expectNear("demoting a null field leaves the volume exactly alone",
               forepeak.waterVolume, 7.25, 0.0);
    expectNear("and the surface offset with it",
               forepeak.surfaceOffset, 1.5, 0.0);
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

    // **Both operands rebuilt from the compartment, not read back from the
    // demote path.** This asserted `forepeak.waterVolume / gridPlanArea` where
    // the volume had just been *written* by `demoteWater` and the area came from
    // the engine's own grid -- the demote path's arithmetic checked against the
    // demote path's own outputs, which cannot fail however wrong either is.
    //
    // The grid is built with `ceil((hi - lo) / h)` cells per axis at the cell
    // size `promoteWater` uses, so its plan area is reconstructible from the
    // compartment's bounding box alone. And the volume to divide is the one the
    // test *set*, which also closes the round trip: a demote that lost mass now
    // fails here as well as in the mass test.
    const double h = 0.05;   // m, `promoteWater`'s cell size
    const double nx = std::ceil((forepeak.bboxHi.x - forepeak.bboxLo.x) / h);
    const double ny = std::ceil((forepeak.bboxHi.y - forepeak.bboxLo.y) / h);
    const double rebuiltPlanArea = nx * ny * h * h;

    expectNear("the grid's plan area is the compartment's box, discretised",
               gridPlanArea, rebuiltPlanArea, 1e-9);

    const double expectedOffset = initialVolume / rebuiltPlanArea;

    // 1e-6 m on an offset of ~0.02 m is a bound on the arithmetic, not on the
    // discretisation -- the cell count is reconstructed exactly above, so there
    // is no discretisation error left to absorb.
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

    // **The tie is exact, and this loop used to assert it at 1e-6 and then assert
    // nothing whatever about the order that came out of it** -- a test that named
    // the hazard and did not check it. `rollRate` and `accel` are properties of
    // the *ship*, computed once above `waterCandidates`' compartment loop, so
    // every qualifying compartment divides the same two doubles by the same two
    // thresholds and takes `std::min` of the same pair. The scores are the same
    // *bits*, and zero is therefore the tolerance the measurement supports.
    //
    // Asserted at zero rather than at 1e-6 for the reason a loose bound is nearly
    // a vacuous one: a 1e-6 also passes on scores that genuinely differ, and that
    // is precisely the case in which the order below would be settled by the
    // score comparison instead of by the tie-break these tests are aimed at.
    int qualifying = 0;
    for (const WaterCandidate& c : rev.considered) {
        if (c.score <= 0) continue;
        ++qualifying;
        expectTrue("every qualifying score is bit-identical to the first",
                   c.score == rev.considered[0].score);
    }
    expectEqual("and exactly the three flooded compartments qualify", qualifying, 3);

    // Which leaves the whole order to the tie-break -- and the order is what
    // `review()` spends the budget along, so it is a decision and not a
    // presentation detail. See `testTiedCandidatesAreOrderedByCompartmentIndex`.
    expectEqual("the tie is broken by compartment index: first is 0",
                rev.considered[0].compartment, 0);
    expectEqual("second is 1", rev.considered[1].compartment, 1);
    expectEqual("third is 2", rev.considered[2].compartment, 2);
}

// The order `waterCandidates` hands back is a decision, not a presentation
// detail: `WaterPromoter::review` walks `considered` from the front, promotes
// until the particle or tile budget is spent, and then `break`s. It picks the
// compartments that get FLIP water.
//
// And on this tier that sort is almost nothing but ties. The two motion
// quantities are properties of the ship, computed once above the compartment
// loop, so `score` takes one of two values over the whole list: one
// bit-identical number on every qualifying compartment, and zero on all the rest.
//
// **Against the implementation this replaced -- `return a.score > b.score` and
// no second leg -- the resulting order was undefined.** `std::sort` is not
// stable, and the ferry's 18 compartments are past libstdc++'s 16-element
// insertion-sort threshold, so introsort partitions the range about a
// median-of-three pivot and the tied runs come back permuted. Measured on this
// exact fixture, on GCC 16.1.1 / libstdc++, the old comparator returned
//
//     1 2 0 17 16 15 14 13 12 11 10 9 8 7 6 5 4 3
//
// where the fixed one returns 0 1 2 3 ... 17. That permutation is a function of
// the compartment count and of the standard library's internals; nothing about
// the ship is in it, and it decided which compartments were resolved --
// `testTheBudgetSpendsTheTieBreakOrder` below is the same fixture with room for
// one promotion, and it promoted compartment **1** for exactly this reason.
void testTiedCandidatesAreOrderedByCompartmentIndex() {
    std::printf("\n   tied candidates are ordered by compartment index\n");

    Ship ferry = ferryAfloat();
    floodCompartment(ferry, 0, 5.0);
    floodCompartment(ferry, 1, 5.0);
    floodCompartment(ferry, 2, 5.0);
    setRollRate(ferry, 0.15);   // 3x threshold, and the only trigger that fires

    WaterCriterion crit;
    WaterPromoter promoter(crit);
    const WaterReview rev = promoter.review(ferry);

    // The non-vacuity guard, and it is about the *length* of the list rather
    // than its contents. libstdc++ hands any range of 16 or fewer elements
    // straight to an insertion sort, which is stable by construction -- so on a
    // fifteen-compartment fixture the old comparator would also have returned
    // index order and this test would prove nothing at all. The ferry has 18.
    expectTrue("the fixture is long enough that introsort really partitions it",
               rev.considered.size() > 16);

    // The expected order is closed-form rather than observed. The three flooded
    // compartments are the ship's first three and they alone score; every other
    // compartment is dry and scores exactly zero. Score-descending and
    // index-ascending therefore agree, and the list must come back as the
    // identity permutation.
    bool inIndexOrder = true;
    bool scoresDescend = true;
    for (size_t i = 0; i < rev.considered.size(); ++i) {
        if (rev.considered[i].compartment != static_cast<int>(i)) inIndexOrder = false;
        if (i > 0 && rev.considered[i].score > rev.considered[i - 1].score) scoresDescend = false;
    }
    if (!inIndexOrder) {
        std::printf("      order was:");
        for (const WaterCandidate& c : rev.considered) std::printf(" %d", c.compartment);
        std::printf("\n");
    }
    expectTrue("scores are non-increasing along the list", scoresDescend);
    expectTrue("and every tie is broken by ascending compartment index", inIndexOrder);
}

// Score is still the primary key; the tie-break only separates equals.
//
// Flooding compartments 2 and 5 -- and neither end of the ship -- separates the
// two orderings. Ranked by score the pair leads; ranked by index alone
// compartments 0 and 1 would. A comparator that had lost the score leg, or had
// its two legs the wrong way round, passes the test above and fails here.
//
// The zero-score remainder is asserted as well, and that is not padding: with
// only two compartments scoring, the old comparator happened to leave *those
// two* in index order on this fixture and scrambled the fifteen behind them.
// Checking the head alone would have been a green test over a broken sort.
void testScoreOutranksTheTieBreak() {
    std::printf("\n   score outranks the tie-break\n");

    Ship ferry = ferryAfloat();
    floodCompartment(ferry, 2, 5.0);
    floodCompartment(ferry, 5, 5.0);
    setRollRate(ferry, 0.15);

    WaterCriterion crit;
    WaterPromoter promoter(crit);
    const WaterReview rev = promoter.review(ferry);

    expectTrue("the whole ship is considered", rev.considered.size() > 16);
    expectEqual("the flooded pair leads, the lower index of the two first",
                rev.considered[0].compartment, 2);
    expectEqual("then the other one", rev.considered[1].compartment, 5);
    expectTrue("both on a positive score", rev.considered[1].score > 0.0);
    expectTrue("and the third candidate scores zero", rev.considered[2].score == 0.0);

    // The zero-score remainder follows in index order, with the two that were
    // ranked out of it skipped. Closed-form for the same reason as above: every
    // one of these scores is exactly zero, so only the tie-break orders them.
    bool restInIndexOrder = true;
    int expected = 0;
    for (size_t i = 2; i < rev.considered.size(); ++i) {
        while (expected == 2 || expected == 5) ++expected;
        if (rev.considered[i].compartment != expected) restInIndexOrder = false;
        ++expected;
    }
    if (!restInIndexOrder) {
        std::printf("      order was:");
        for (const WaterCandidate& c : rev.considered) std::printf(" %d", c.compartment);
        std::printf("\n");
    }
    expectTrue("and the zero-score remainder is in index order too", restInIndexOrder);
}

// The consequence, which is why any of this matters. `review()` promotes down
// `considered` until the budget is exhausted and then stops, so the tie-break is
// what chooses the compartments that get resolved. Three tied compartments and
// room for one: the lowest index wins, and it wins on every run and every build
// rather than on whichever one introsort happened to move to the front.
void testTheBudgetSpendsTheTieBreakOrder() {
    std::printf("\n   the budget is spent in tie-break order\n");

    Ship ferry = ferryAfloat();
    floodCompartment(ferry, 0, 5.0);
    floodCompartment(ferry, 1, 5.0);
    floodCompartment(ferry, 2, 5.0);
    setRollRate(ferry, 0.15);

    WaterCriterion crit;
    crit.dwell = 1;

    // Room for exactly one, sized off the estimator rather than off a constant,
    // for the reason `testParticleBudgetEnforcement` gives: a hardcoded number
    // here is a guess about the estimator's internals, and it lands on the wrong
    // side of the boundary in silence when they change.
    int oneParticles = 0, oneTiles = 0;
    estimateFlipCost(ferry.compartments[0], crit, oneParticles, oneTiles);
    expectTrue("one compartment needs particles and tiles",
               oneParticles > 0 && oneTiles > 0);
    crit.particleBudget = 2 * oneParticles - 1;   // fits one, one short of two
    crit.tileBudget = 2 * oneTiles - 1;

    WaterPromoter promoter(crit);
    const WaterReview rev = promoter.review(ferry);

    std::printf("      qualifying order:");
    for (const WaterCandidate& c : rev.considered)
        if (c.score > 0) std::printf(" %d", c.compartment);
    std::printf("\n");

    expectEqual("only one of the three fits", static_cast<int>(rev.promoted.size()), 1);
    expectTrue("and the budget said so out loud", rev.problems.size() > 0);

    // Read through a guard rather than indexed directly: if the budget ever
    // admits none, this reports the wrong compartment instead of walking off the
    // end of the vector, and the check count is the same either way.
    const int promoted = rev.promoted.empty() ? -1 : rev.promoted[0].compartment;
    expectEqual("and it is the lowest-index qualifying compartment", promoted, 0);
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

// **A threshold of zero refuses, rather than dividing by it.**
//
// The score is `rate / rollRatePromote`, and a zero threshold does not make
// everything qualify -- it makes the score `+inf` for a moving ship and a NaN for
// one at rest, which is `computeLateralAccel`'s documented answer when `dt <= 0`.
// A NaN then goes into the `a.score > b.score` comparator the ranking uses, and a
// comparator that is not a strict weak ordering is undefined behaviour rather than
// a wrong number. Both sibling tiers guard the identical construct
// (`promotion.cpp:271` and `:830`); this one did not.
void testAZeroThresholdRefusesRatherThanDividingByIt() {
    std::printf("\n   a zero promote threshold is refused\n");

    Ship ferry = ferryAfloat();
    setRollRate(ferry, 0.30);
    ferry.compartments[0].waterVolume = 8.0;

    WaterCriterion sane;
    const std::vector<WaterCandidate> found = waterCandidates(ferry, sane, Vec3{}, 0.02);
    expectTrue("the ship qualifies on a real threshold, or this proves nothing",
               !found.empty());

    // Not empty -- this list is what was *considered*, and a compartment dropped
    // from it is a silent refusal. Scored zero, with the reason, and above all
    // finite: the failure being guarded against is a NaN score reaching the sort.
    const auto refused = [&](WaterCriterion c) {
        const std::vector<WaterCandidate> got = waterCandidates(ferry, c, Vec3{}, 0.02);
        expectTrue("the compartment is still considered", !got.empty());
        bool named = false;
        for (const WaterCandidate& w : got) {
            // The failure being guarded against is a NaN reaching the sort, so this
            // is the assertion that matters and it applies to every entry.
            expectTrue("every score is finite", std::isfinite(w.score));
            expectNear("and none of them scored", w.score, 0.0, 0.0);
            if (w.why.find("threshold is not positive") != std::string::npos) named = true;
        }
        // Not every entry: a compartment too large to resolve, or below the motion
        // thresholds, is refused earlier and carries its own reason. What has to be
        // true is that the ones reaching the divide say why they did not.
        expectTrue("and the threshold is named by the ones that reached it", named);
    };
    WaterCriterion noRoll = sane;
    noRoll.rollRatePromote = 0.0;
    refused(noRoll);
    WaterCriterion noAccel = sane;
    noAccel.accelPromote = 0.0;
    refused(noAccel);
}

// **The tile half of `estimateFlipCost` was never compared with anything.** The
// particle half is checked against the solver's real seeding a few tests up;
// `predictedTiles` was computed there and dropped on the floor, and the only tile
// assertions in this file were `tilesPerM3 > 0` and `oneTiles > 0`.
//
// `water_promotion.hpp` records a disagreement in prose -- "by the solver's own
// count rather than the estimator's 12 500, 17 019 tiles" -- and calls the estimate
// "self-consistent" because the budget is denominated in it. The 1.36x there is
// measured on a **cube** of water, by `water_probe --cost`. On a real ship
// compartment it is not 1.36x and it is not a factor at all:
//
//     forepeak,  3 m3:  billed 375, allocated 12 300  -- 32.8x
//     forepeak, 12 m3:  billed 1500, allocated 12 300 --  8.2x
//
// **The allocated count is the same at both volumes.** The estimator bills
// `V / (64 h^3)`, strictly proportional to the water. The solver allocates 4x4x4
// tiles over the compartment's *footprint*, and the ferry's compartments are wide
// and shallow -- 3 m3 over the forepeak's 232 m2 is 13 mm deep and 12 m3 is 52 mm,
// so both are one tile layer thick and the count is the plan area, four times over.
//
// So the estimate is not merely low. It is a function of the wrong variable, which
// is the same defect this tier was already caught by once when `coreSecondsPerCompartment`
// turned out to be a per-compartment constant where the truth scales with the water.
// Here it is the mirror image: a per-volume term where the truth is an area.
//
// The consequence is the budget. `tileBudget = 12500` nominally admits 100 m3; a
// single forepeak holding 3 m3 already allocates 12 300 real tiles -- 98% of that
// capacity for 3% of the volume. What the budget means depends on the shape of the
// compartment it is spent on, which is the property a budget is supposed not to have.
void testTheTileEstimatorIsAFunctionOfTheWrongVariable() {
    const auto measure = [](double volume, int& billed, double& depth) {
        Ship ferry = ferryAfloat();
        Compartment& forepeak = ferry.compartments[0];
        forepeak.waterVolume = volume;
        const Vec3 span = forepeak.bboxHi - forepeak.bboxLo;
        depth = volume / std::max(1e-9, span.x * span.y);
        WaterCriterion crit;
        int particles = 0;
        estimateFlipCost(forepeak, crit, particles, billed);
        auto field = promoteWater(forepeak, ferry, crit);
        if (!field) return 0;
        // One step builds the sparse structure; `tiles()` is empty before it.
        flip::Params params;
        flip::Solver solver;
        flip::Account account;
        solver.step(*field, 0.001, params, account);
        return solver.tiles();
    };

    int smallBilled = 0, largeBilled = 0;
    double smallDepth = 0, largeDepth = 0;
    const int smallReal = measure(3.0, smallBilled, smallDepth);
    const int largeReal = measure(12.0, largeBilled, largeDepth);
    std::printf("     forepeak tiles, billed against allocated: 3 m3 %d vs %d (%.1fx,"
                " %.0f mm deep), 12 m3 %d vs %d (%.1fx, %.0f mm)\n",
                smallBilled, smallReal, static_cast<double>(smallReal) / smallBilled,
                1000.0 * smallDepth, largeBilled, largeReal,
                static_cast<double>(largeReal) / largeBilled, 1000.0 * largeDepth);

    // Guards first: a refused promotion would make every ratio zero, and the
    // forepeak refuses below about 2.9 m3 because the water is thinner than the
    // lowest particle plane -- which is why 3 m3 is the small case and not 1.
    expectTrue("both volumes promoted", smallReal > 0 && largeReal > 0);
    expectTrue("and both were billed for something", smallBilled > 0 && largeBilled > 0);
    // The billed side is exactly proportional to volume, which is the thing being
    // contrasted. Asserted so the contrast cannot quietly become a comparison of
    // two quantities that both vary.
    expectEqual("the estimator bills strictly in proportion to the water",
                static_cast<long long>(4 * smallBilled), static_cast<long long>(largeBilled));
    // And the allocated side does not move at all, because both depths are inside
    // one 4x4x4 tile at h = 0.05, i.e. 0.2 m.
    expectTrue("both fills are thinner than one tile", largeDepth < 0.2);
    expectEqual("so the solver allocates the same tiles for four times the water",
                static_cast<long long>(smallReal), static_cast<long long>(largeReal));
    expectTrue("and the estimate is short by a large factor at both, not by 1.36x",
               smallReal > 8 * smallBilled && largeReal > 4 * largeBilled);
}

void runWaterPromotionTests() {
    std::printf("\n--- water promotion ---\n");

    // Section 1: Criterion
    testShipAtRestPromotesNothing();
    testAZeroThresholdRefusesRatherThanDividingByIt();
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
    testDrainingACompartmentDemotesItAtOnce();
    testAccelerationAloneHoldsACompartment();
    testMinDepthIsNotEnforced();
    testCopiedShipHasNoPromotedWater();
    testBudgetsAdmitTheSameVolume();
    testTheRoundTripConservesAtEveryVolume();
    testTheEstimateAndTheAllocationDifferByTheDocumentedFactor();
    testParticleBudgetEnforcement();
    testTileBudgetEnforcement();
    testBudgetAccountingAcrossReviews();
    testABudgetStarvedCompartmentKeepsItsDwell();

    // Section 4: State Transfer
    testRoundTripMassConservation();
    testEmptyCompartmentReturnsNull();
    testCentroidPreservation();
    testSurfaceOffsetMatchesVolume();

    // Section 5: Multiple Compartments
    testMultipleCompartmentsRankedByScore();
    testTiedCandidatesAreOrderedByCompartmentIndex();
    testScoreOutranksTheTieBreak();
    testTheBudgetSpendsTheTieBreakOrder();
    testDifferentCompartmentsSeparateHysteresis();

    // Section 6: Cost and Performance
    testCostReporting();
    testCounters();
    testClearResetsState();
    testTheTileEstimatorIsAFunctionOfTheWrongVariable();
}

