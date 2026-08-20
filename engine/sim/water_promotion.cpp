// SPDX-License-Identifier: MIT

#include "water_promotion.hpp"
#include "ship.hpp"
#include "flip.hpp"
#include "engine/core/math.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>

namespace sim {
namespace promotion {

// --- Helpers --------------------------------------------------------------

// Roll is rotation about the ship's *own* longitudinal axis, and
// `state.angularVelocity` is a world vector. Taking its x component treats the
// two frames as the same, which they are only while the ship is upright -- and
// the whole point of this criterion is the case where she is not. A ferry lolled
// to 58 degrees has her bow axis nowhere near world x.
//
// The same confusion, one frame further out, is already in CLAUDE.md's table:
// `state.velocity.x` read as "speed" when it is a world vector.
static double computeRollRate(const Ship& ship) {
    const Mat3 R = ship.state.orientation.toMat3();   // body -> world
    // The body x axis expressed in world coordinates; the roll rate is the
    // component of the world angular velocity along it.
    const Vec3 bowAxis{R.m[0], R.m[3], R.m[6]};
    return std::abs(dot(ship.state.angularVelocity, bowAxis));
}

// Lateral acceleration, differenced from the velocity the caller last saw.
//
// **This was `length(state.velocity)` compared against a threshold in m/s^2** --
// a speed standing in for an acceleration, so the units on the two sides of the
// comparison never matched. `water_probe`'s beam-sea control read 2.26 "m/s2" at
// 2 m of wave, which was 2.26 m/s of orbital speed: the criterion was firing on
// how fast the ship was moving rather than on how hard she was being thrown.
//
// A ship drifting steadily at 3 m/s has no lateral acceleration at all and would
// have promoted every compartment aboard; a ship snapped sideways from rest by a
// collision has a large one and, at the instant it matters most, a small speed.
// The sign of that error is the dangerous one -- it fires when nothing is
// happening and stays quiet when something is.
//
// `dt <= 0` means the caller cannot say how much time passed, which is the case
// for every test that reviews a static ship. There is no acceleration to be had
// from a single sample, so it reports zero rather than inventing one.
static double computeLateralAccel(const Ship& ship, const Vec3& previousVelocity,
                                  double dt) {
    if (dt <= 0) return 0;
    const Vec3 dv = ship.state.velocity - previousVelocity;
    const Mat3 R = ship.state.orientation.toMat3();
    // Body y (to port) in world coordinates: the transverse direction, which is
    // the one water sloshes along and the one a ro-pax capsizes about.
    const Vec3 portAxis{R.m[1], R.m[4], R.m[7]};
    return std::abs(dot(dv, portAxis)) / dt;
}

// Mean depth if the water were spread flat over the compartment's whole
// bounding-box footprint. That is a *floor*, not the depth anywhere in
// particular: the ferry's compartments are tens of metres on a side, so 5 m3
// in the forepeak reads 0.0215 m however deep the water actually stands
// against a bulkhead once she heels. `criterion.minDepth` of 0.5 m is
// therefore not a quantity this function can be compared against as a
// promotion gate -- doing so silently filtered *every* compartment on the
// ship, and the review came back with an empty candidate list. Only
// `minVolume` gates promotion until there is a real free-surface calculation
// here; see the TODO below.
static double computeWaterDepth(const Compartment& comp) {
    if (comp.waterVolume <= 0) return 0;

    // Rough estimate: depth = volume / horizontal cross-section area.
    // Compute bounding box from mesh vertices.
    // TODO: proper free-surface calculation from TriMesh.
    if (comp.mesh.verts.empty()) return 0;

    Vec3 min = comp.mesh.verts[0];
    Vec3 max = comp.mesh.verts[0];
    for (const auto& v : comp.mesh.verts) {
        min.x = std::min(min.x, v.x);
        min.y = std::min(min.y, v.y);
        min.z = std::min(min.z, v.z);
        max.x = std::max(max.x, v.x);
        max.y = std::max(max.y, v.y);
        max.z = std::max(max.z, v.z);
    }

    double footprint = (max.x - min.x) * (max.y - min.y);
    if (footprint <= 0) return 0;

    return comp.waterVolume / footprint;
}

// --- 1. Candidate generation ----------------------------------------------

void estimateFlipCost(const Compartment& comp, const WaterCriterion& /*criterion*/,
                      int& particles, int& tiles) {
    // Heuristic: ~1000 particles per m³ for reasonable resolution.
    //
    // Rounded on the same terms as the tile count below. A multiplication by
    // 1000 happens to be exact for the volumes tried here where the tile
    // division is not, so this one loses nothing today -- which is a fact about
    // the arithmetic and not a property worth depending on.
    particles = static_cast<int>(std::llround(comp.waterVolume * 1000.0));

    // Tiles: 4×4×4 cells each, cell size from criterion.solver default grid.
    // We don't have a grid yet, so use a reasonable default cell size.
    // FLIP typically uses h = 0.05 m (5 cm cells) as the default.
    // Wetted volume / (64 * h³) gives approximate tile count.
    //
    // **Rounded rather than truncated, and the difference is not cosmetic.**
    // `10.0 / (64 * 0.05³)` is 1249.9999999999998 in binary floating point, not
    // 1250, so a `static_cast<int>` threw a whole tile away on an exact-looking
    // input -- and the budget derived from it then admitted 100.08 m³ where the
    // particle budget admitted 100.00. That 0.08% disagreement is the same defect
    // shape as the 6.25× one this pair already had, arriving through a cast
    // instead of through a constant. `testBudgetsAdmitTheSameVolume` only saw it
    // once it asked `estimateFlipCost` for the rate instead of re-typing it.
    double h = 0.05;  // m, typical FLIP cell edge length
    double cellVolume = h * h * h;
    double tilesExact = comp.waterVolume / (64.0 * cellVolume);
    tiles = std::max(1, static_cast<int>(std::llround(tilesExact)));
}

std::vector<WaterCandidate> waterCandidates(const Ship& ship, const WaterCriterion& criterion,
                                            const Vec3& previousVelocity, double dt) {
    std::vector<WaterCandidate> candidates;

    double rollRate = computeRollRate(ship);
    double accel = computeLateralAccel(ship, previousVelocity, dt);

    for (size_t i = 0; i < ship.compartments.size(); ++i) {
        const auto& comp = ship.compartments[i];

        // Skip dry or nearly-dry compartments
        if (comp.waterVolume < criterion.minVolume) {
            WaterCandidate cand;
            cand.compartment = static_cast<int>(i);
            cand.name = comp.name;
            cand.rollRate = rollRate;
            cand.accel = accel;
            cand.depth = 0;
            cand.volume = comp.waterVolume;
            cand.score = 0;
            cand.why = "insufficient volume";
            candidates.push_back(cand);
            continue;
        }

        WaterCandidate cand;
        double depth = computeWaterDepth(comp);
        cand.compartment = static_cast<int>(i);
        cand.name = comp.name;
        cand.rollRate = rollRate;
        cand.accel = accel;
        cand.depth = depth;
        cand.volume = comp.waterVolume;

        estimateFlipCost(comp, criterion, cand.particles, cand.tiles);
        cand.cost = criterion.coreSecondsPerCompartment;

        // Too large to resolve at this grid spacing, and that is a property of
        // the compartment rather than of how busy the budget happens to be. It
        // gets its own reason: a vehicle deck refused because it cannot be
        // afforded is a finding about this tier, and it read as a transient
        // "budget exhausted" for as long as the two were the same branch.
        if (comp.waterVolume > criterion.maxVolumePerCompartment) {
            cand.score = 0;
            cand.why = "too large to resolve at this grid spacing";
            candidates.push_back(cand);
            continue;
        }

        // Check motion thresholds
        bool rollQualifies = rollRate >= criterion.rollRatePromote;
        bool accelQualifies = accel >= criterion.accelPromote;

        if (!rollQualifies && !accelQualifies) {
            cand.why = "below motion thresholds";
            cand.score = 0;
            candidates.push_back(cand);  // Include in considered with score 0
            continue;
        }

        // Score: how far past threshold is the **weaker** trigger?
        // Both above → score is min of the two ratios.
        // One above → score is that ratio.
        // **Both thresholds refused before they are divided by, and reported rather
        // than dropped.** A threshold of zero does not make everything qualify: it
        // makes the score `+inf` for a moving ship and a NaN for one at rest, which
        // is `computeLateralAccel`'s own answer when `dt <= 0`. A NaN then reaches
        // the `a.score > b.score` comparator below, and a comparator that is not a
        // strict weak ordering is undefined behaviour rather than a wrong number.
        // Both sibling tiers guard the identical construct (`promotion.cpp:271`
        // and `:830`).
        //
        // Kept in the list with a reason, not `continue`d past. This list is what
        // a caller reads to see what was *considered*, and a compartment that
        // vanishes from it is the silent refusal this tier was already caught
        // doing once, 154 times over.
        if (!(criterion.rollRatePromote > 0) || !(criterion.accelPromote > 0)) {
            cand.score = 0;
            cand.why = "a promotion threshold is not positive, so no score can be formed";
            candidates.push_back(cand);
            continue;
        }
        double rollScore = rollRate / criterion.rollRatePromote;
        double accelScore = accel / criterion.accelPromote;

        if (rollQualifies && accelQualifies) {
            cand.score = std::min(rollScore, accelScore);
        } else if (rollQualifies) {
            cand.score = rollScore;
        } else {
            cand.score = accelScore;
        }

        cand.why = "qualifies";
        candidates.push_back(cand);
    }

    // Most urgent first, and **ties broken by compartment index** -- the same
    // second leg both sibling tiers carry (`promotion.cpp:332` and `:349` for the
    // structural zones, `:871` for the gas ones), for the reason written down
    // above the first of them: "Deterministic ties, by panel index, so two runs
    // of the same load promote the same patches."
    //
    // On this tier a tie is not the unlucky case, it is the only case. `rollRate`
    // and `accel` are ship-level quantities, computed once above the compartment
    // loop, so `rollScore` and `accelScore` are the *same two doubles* for every
    // compartment and `std::min` of them is one bit-identical `score` on every
    // qualifying compartment, with zero on all the rest. A comparator on `score`
    // alone therefore ties across the whole list.
    //
    // `std::sort` is not stable, so what came back was introsort's permutation
    // rather than the ship's order. Measured on the ferry -- 18 compartments,
    // just past libstdc++'s 16-element insertion-sort threshold, so the range is
    // really partitioned -- three flooded holds at indices 0, 1, 2 came back as
    // `1 2 0`, and the rest of the list reversed. The order was a function of the
    // compartment count and of which standard library built it, and of nothing
    // else; add a compartment to the ferry and it moves.
    //
    // That order is load-bearing. `WaterPromoter::review` walks `considered` from
    // the front, promotes until the particle or tile budget is spent and then
    // stops, so *which* compartments get FLIP water was being decided by a
    // tie-break that did not exist -- on the same fixture, with room for one
    // promotion, it resolved compartment 1 and left compartment 0 quiescent.
    // `testTiedCandidatesAreOrderedByCompartmentIndex` and
    // `testTheBudgetSpendsTheTieBreakOrder` pin both halves of that.
    std::sort(candidates.begin(), candidates.end(),
              [](const WaterCandidate& a, const WaterCandidate& b) {
                  if (a.score != b.score) return a.score > b.score;
                  return a.compartment < b.compartment;
              });

    return candidates;
}

// --- 2. State machine -----------------------------------------------------

WaterPromoter::WaterPromoter(WaterCriterion criterion)
    : criterion_(std::move(criterion)) {}

WaterReview WaterPromoter::review(const Ship& ship, double dt) {
    auto t0 = std::chrono::steady_clock::now();

    // The first review has no earlier velocity to difference against, so it
    // reports no acceleration however large dt is. Promoting on the strength of
    // a difference taken from a default-constructed zero would read the ship's
    // whole speed as one review's worth of acceleration -- the same units error
    // this replaced, arriving once per promoter rather than once per review.
    const double accelDt = havePreviousVelocity_ ? dt : 0.0;

    WaterReview result;
    result.considered = waterCandidates(ship, criterion_, previousVelocity_, accelDt);

    reviews_++;

    // --- Step 1: Check active compartments for demotion ---

    double rollRate = computeRollRate(ship);
    double accel = computeLateralAccel(ship, previousVelocity_, accelDt);

    previousVelocity_ = ship.state.velocity;
    havePreviousVelocity_ = true;

    std::vector<WaterActive> stillActive;

    for (auto& act : active_) {
        // Update current state
        const auto& comp = ship.compartments[act.compartment];
        act.rollRate = rollRate;
        act.accel = accel;
        act.depth = computeWaterDepth(comp);
        act.volume = comp.waterVolume;

        // Check hold thresholds (lower than promote)
        bool rollHolds = rollRate >= criterion_.rollRateHold;
        bool accelHolds = accel >= criterion_.accelHold;
        bool volumeHolds = act.volume >= criterion_.minVolume;

        if (!rollHolds && !accelHolds) {
            act.idleReviews++;
        } else {
            act.idleReviews = 0;
        }

        // Demote if idle too long OR volume drops below minimum
        if (act.idleReviews >= criterion_.hold || !volumeHolds) {
            result.demoted.push_back(act);
            demotions_++;
        } else {
            stillActive.push_back(act);
        }
    }

    active_ = std::move(stillActive);

    // --- Step 2: Dwell, which counts qualification and nothing else ---

    // **Its own pass, ahead of the budget, and that separation is the whole fix.**
    // `qualifying_` holds the consecutive reviews a compartment has *qualified*
    // for, which is what its own comment has always said it held. It used to be
    // built inside the promotion loop below, and that loop stops at the first
    // candidate the budget cannot admit -- so every candidate ranked below the stop
    // was never written back. Its streak did not pause there: it went to zero and
    // restarted at one on the next review. A compartment held under a full budget
    // could therefore never reach `dwell` however long it qualified. It was not
    // deferred, it was repeatedly restarted, and restarted by a budget decision
    // taken about a *different* compartment ranked above it.
    //
    // Measured on three tied 5 m3 compartments with room for one and `dwell = 2`.
    // Compartment 0 promotes on review 2. Compartments 1 and 2 then *oscillated*:
    // each rebuilt to a streak of two on every even review, was refused there and
    // dropped in the same breath, and stood at one again on every odd review --
    // sixteen reviews of continuous qualification and not one promotion. The nastier
    // consequence is the parity. Whether a starved compartment ever promoted at all
    // depended on whether the budget happened to free on an even review or an odd
    // one, which is not a property of the ship, of the water, or of the criterion.
    //
    // With the streak kept they stand at nine consecutive refusals by review 10 and
    // the first of them promotes on the very review that frees the budget, whichever
    // review that is. `testABudgetStarvedCompartmentKeepsItsDwell` pins both halves
    // and frees the budget on an odd review on purpose, because an even one passes
    // against the defect.
    //
    // That is the gas tier's defect running backwards. There a refusal was silent
    // and 154 of them read as hysteresis working correctly; here the hysteresis is
    // silently *reset* by a decision that says nothing about whether the
    // compartment qualifies. Both siblings build the streak in a separate pass over
    // the whole ranked list before any budget is consulted (`promotion.cpp:441` for
    // the structural zones, `:967` for the gas ones) -- the structure the header
    // already claimed this tier followed "exactly".
    //
    // An already-active compartment leaves the list rather than accumulating,
    // which is this tier's existing behaviour and is kept as it was. It reads as a
    // divergence from the siblings, which have no such skip, and it is not one in
    // outcome: every route out of `active_` here -- motion below the hold
    // thresholds for `hold` consecutive reviews, or volume below `minVolume` --
    // also drops the compartment out of the candidate list, so anything demoted has
    // already had its streak cleared by ceasing to qualify. The skip saves carrying
    // a count that nothing can read; it does not change when a compartment
    // re-promotes.
    std::vector<std::pair<int, int>> stillQualifying;
    for (const WaterCandidate& cand : result.considered) {
        if (cand.score <= 0) continue;  // didn't qualify this review

        const bool alreadyActive =
            std::any_of(active_.begin(), active_.end(), [&](const WaterActive& a) {
                return a.compartment == cand.compartment;
            });
        if (alreadyActive) continue;

        const auto it = std::find_if(qualifying_.begin(), qualifying_.end(),
                                     [&](const auto& p) { return p.first == cand.compartment; });
        stillQualifying.push_back(
            {cand.compartment, (it != qualifying_.end()) ? it->second + 1 : 1});
    }
    qualifying_ = std::move(stillQualifying);

    // --- Step 3: Promotion, down the ranked order, until the budget is spent ---

    // `max(1, dwell)`, the way both siblings write it (`promotion.cpp:475` and
    // `:989`): a streak starts at one, so zero and one are the same criterion, and
    // a negative dwell is not a licence to promote before qualifying. The
    // starvation report below measures against the same value, so what "promotable"
    // means cannot drift between the decision and the report.
    const int dwellNeeded = std::max(1, criterion_.dwell);

    int totalParticles = 0, totalTiles = 0;
    for (const auto& a : active_) {
        totalParticles += a.particles;
        totalTiles += a.tiles;
    }

    // The budget is spent from the front of the ranked list and promotion stops at
    // the first compartment that does not fit; it is not repacked with a smaller
    // one from further down, because the order is the priority order and skipping
    // ahead would resolve a 1 m3 trickle in place of a compartment that just
    // missed. `testTheBudgetSpendsTheTieBreakOrder` pins that.
    //
    // Stopping *promoting* is not the same as stopping *looking*, and the
    // difference is what makes the starvation visible. The scan runs to the end of
    // the list either way, so every compartment that had served its dwell and got
    // nothing can be named -- not only the first one the budget happened to test.
    // The loop used to `break`, which is why this was invisible and why the streaks
    // below it were being destroyed.
    bool budgetSpent = false;
    for (WaterCandidate& cand : result.considered) {
        if (cand.score <= 0) continue;

        const auto it = std::find_if(qualifying_.begin(), qualifying_.end(),
                                     [&](const auto& p) { return p.first == cand.compartment; });
        if (it == qualifying_.end()) continue;   // already active, so not a candidate
        const int streak = it->second;
        if (streak < dwellNeeded) continue;      // still building its dwell

        if (!budgetSpent && totalParticles + cand.particles <= criterion_.particleBudget &&
            totalTiles + cand.tiles <= criterion_.tileBudget) {
            WaterActive newActive;
            newActive.compartment = cand.compartment;
            newActive.name = cand.name;
            newActive.rollRate = cand.rollRate;
            newActive.accel = cand.accel;
            newActive.depth = cand.depth;
            newActive.volume = cand.volume;
            newActive.score = cand.score;
            newActive.particles = cand.particles;
            newActive.tiles = cand.tiles;
            newActive.cost = cand.cost;
            newActive.promotedAtReview = reviews_;
            newActive.idleReviews = 0;

            active_.push_back(newActive);
            result.promoted.push_back(newActive);
            promotions_++;
            totalParticles += cand.particles;
            totalTiles += cand.tiles;
            continue;
        }

        // Starved: it qualified, it served its dwell, and it got nothing. How long
        // that has been true is the number worth carrying -- one review is the
        // budget doing its job, and a count that climbs on every review is a
        // compartment that will never be resolved while the ones above it hold
        // their water.
        budgetSpent = true;
        WaterStarved starved;
        starved.compartment = cand.compartment;
        starved.name = cand.name;
        starved.qualifyingReviews = streak;
        starved.refusedReviews = streak - dwellNeeded + 1;
        result.starved.push_back(starved);

        cand.why = "qualifies, and the particle or tile budget refused it on " +
                   std::to_string(starved.refusedReviews) +
                   " consecutive review(s); its dwell is retained";
    }

    // **The refusal gets a channel, on this tier's own terms and on the gas tier's
    // both.** `considered` is where a rejected candidate says why -- that is
    // already the answer here for one refused on volume or on grid spacing -- and
    // `problems` carries the count, which is how `gasCandidates` reports its
    // refusals after the same defect was found there. The message this replaced,
    // "particle or tile budget exhausted", named no compartment and said nothing
    // about how long: it could not tell a busy review from starvation, and those
    // are the two things a reader of it has to tell apart.
    //
    // The word "budget" is load-bearing rather than decorative: `water_probe`
    // counts the reviews that refused by looking for it in these strings.
    if (!result.starved.empty()) {
        const WaterStarved* longest = &result.starved.front();
        for (const WaterStarved& s : result.starved)
            if (s.refusedReviews > longest->refusedReviews) longest = &s;
        result.problems.push_back(
            "the particle or tile budget refused " + std::to_string(result.starved.size()) +
            " qualifying compartment(s); nothing already promoted was demoted for them, and their"
            " dwell is retained, so each promotes on the first review that admits it (longest: '" +
            longest->name + "', refused on " + std::to_string(longest->refusedReviews) +
            " consecutive review(s))");
    }

    // --- Step 4: Compute totals ---

    for (const auto& a : active_) {
        result.particlesActive += a.particles;
        result.tilesActive += a.tiles;
        result.costActive += a.cost;
    }

    auto t1 = std::chrono::steady_clock::now();
    result.microseconds = std::chrono::duration<double, std::micro>(t1 - t0).count();

    return result;
}

void WaterPromoter::clear() {
    active_.clear();
    qualifying_.clear();
    previousVelocity_ = Vec3{};
    havePreviousVelocity_ = false;
    reviews_ = 0;
    promotions_ = 0;
    demotions_ = 0;
}

// --- 3. State transfer (stubs) --------------------------------------------
//
// These will be implemented in ship.cpp where Ship internals are accessible.
std::unique_ptr<flip::Field> promoteWater(const Compartment& comp, const Ship& /*ship*/,
                                          const WaterCriterion& /*criterion*/) {
    if (comp.waterVolume <= 0) return nullptr;

    // Create FLIP field with compartment bounding box
    auto field = std::make_unique<flip::Field>();
    // Grid cell size: use FLIP default 0.05 m (5 cm cells).
    const double h = 0.05;

    // Compute wetted depth from volume and plan area, the same depth
    // `computeWaterDepth()` reports.
    //
    // **Not "matching `estimateFlipCost`'s assumption of ~1000 particles per m³",
    // which this comment used to claim and which is not what happens.** The
    // seeding below is 2³ per cell, `8/h³` = 64 000 per m³. `WaterCriterion`
    // documents that 64× deliberately -- the estimate is what the *budget* is
    // denominated in and is self-consistent with `tileBudget` -- so the number is
    // not wrong, but a comment here asserting the two agree is.
    const double planArea = (comp.bboxHi.x - comp.bboxLo.x) * (comp.bboxHi.y - comp.bboxLo.y);
    const double wettedDepth = (planArea > 0) ? comp.waterVolume / planArea : 0;

    // **Water too shallow for the lattice to hold is refused, not silently
    // dropped.** `seedBox` puts its lowest plane of particles at `0.5/perAxis`
    // of a cell above the floor -- 0.25 h, or 12.5 mm at h = 0.05 -- and filters
    // every position against the box top. A wetted depth under that admits *no*
    // particle at all, and what followed was not an empty solve but a lost
    // compartment: `setTotalMass` had nothing to distribute the mass over,
    // `demoteWater` read `totalMass()` back as zero, and the water was gone.
    //
    // Measured on the ferry's forepeak, 232 m² of plan area: 2.8 m³ went in and
    // 0.0000 m³ came out, silently, while 3.0 m³ round-tripped exactly. The
    // round-trip test used 5.0 m³. Mass conservation here is exact or it is
    // nothing, so this returns `nullptr` and leaves the compartment quiescent --
    // the same answer it gives for a compartment with no water in it, and one
    // the callers already handle.
    const int particlesPerAxis = 2;
    const double lowestParticlePlane = 0.5 * h / static_cast<double>(particlesPerAxis);
    if (!(wettedDepth > lowestParticlePlane)) return nullptr;

    // Grid constructor: need lo[3], hi[3], and h. Vec3 doesn't have .data().
    double lo[3] = {comp.bboxLo.x, comp.bboxLo.y, comp.bboxLo.z};
    double hi[3] = {comp.bboxHi.x, comp.bboxHi.y, comp.bboxLo.z + wettedDepth};

    // Compute cell count from bounding box and h
    int nx = std::max(1, static_cast<int>(std::ceil((hi[0] - lo[0]) / h)));
    int ny = std::max(1, static_cast<int>(std::ceil((hi[1] - lo[1]) / h)));
    int nz = std::max(1, static_cast<int>(std::ceil((hi[2] - lo[2]) / h)));

    field->grid.h = h;
    field->grid.lo[0] = lo[0];
    field->grid.lo[1] = lo[1];
    field->grid.lo[2] = lo[2];
    field->grid.n[0] = nx;
    field->grid.n[1] = ny;
    field->grid.n[2] = nz;

    // Seed particles: 2^3 lattice per cell, the count the refusal above is derived
    // from -- one constant, so the guard cannot drift from what it guards.
    flip::seedBox(*field, lo, hi, particlesPerAxis, kRhoSeawater);

    // Set exact mass from compartment water volume (m³ → kg)
    const double mass = comp.waterVolume * kRhoSeawater;
    flip::setTotalMass(*field, mass);

    // Initial velocity: particles at rest in body frame
    // TODO: add ship motion contribution when Ship::step() wiring is in place

    return field;
}

void demoteWater(Compartment& comp, std::unique_ptr<flip::Field> field) {
    if (!field) return;
    // **A field with no particles is not a compartment with no water.** Reading
    // `totalMass()` off an empty field gives zero and writing that back destroys
    // whatever the compartment was holding. `promoteWater` no longer produces such
    // a field, but the guard belongs on this side too: this function takes
    // ownership of a field it did not create, and the failure it prevents is the
    // silent loss of the one quantity this tier conserves exactly.
    if (field->particles.empty()) return;

    // Read FLIP state back into compartment quiescent representation.

    // Mass: exact from particle total
    const double mass = field->totalMass();
    comp.waterVolume = mass / kRhoSeawater;

    // Centroid: mass-weighted position
    double cog[3];
    flip::centroid(*field, cog);
    comp.waterCentroid = Vec3{cog[0], cog[1], cog[2]};

    // Surface offset: quiescent level above floor
    // Plan area from mesh horizontal cross-section (TODO: proper from TriMesh)
    // For now use grid plan area as conservative estimate
    const double planArea = field->grid.planArea();
    const double floorZ = comp.bboxLo.z;
    const double level = flip::quiescentLevel(*field, kRhoSeawater, planArea, floorZ);
    comp.surfaceOffset = level - floorZ;

    // field is automatically destroyed when unique_ptr goes out of scope
}

} // namespace promotion
} // namespace sim
