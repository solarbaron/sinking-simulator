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
    particles = static_cast<int>(comp.waterVolume * 1000.0);

    // Tiles: 4×4×4 cells each, cell size from criterion.solver default grid.
    // We don't have a grid yet, so use a reasonable default cell size.
    // FLIP typically uses h = 0.05 m (5 cm cells) as the default.
    // Wetted volume / (64 * h³) gives approximate tile count.
    double h = 0.05;  // m, typical FLIP cell edge length
    double cellVolume = h * h * h;
    double tilesExact = comp.waterVolume / (64.0 * cellVolume);
    tiles = std::max(1, static_cast<int>(tilesExact));
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

    // Sort by score descending: most-urgent compartments first
    std::sort(candidates.begin(), candidates.end(),
              [](const WaterCandidate& a, const WaterCandidate& b) {
                  return a.score > b.score;
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

    // --- Step 2: Check candidates for promotion ---

    // Update qualifying_ map: compartments that qualify this review
    std::vector<std::pair<int, int>> stillQualifying;

    for (const auto& cand : result.considered) {
        if (cand.score <= 0) continue;  // didn't qualify

        // Is it already active?
        bool alreadyActive = std::any_of(active_.begin(), active_.end(),
                                         [&](const WaterActive& a) {
                                             return a.compartment == cand.compartment;
                                         });
        if (alreadyActive) continue;

        // Find in qualifying_ list
        auto it = std::find_if(qualifying_.begin(), qualifying_.end(),
                               [&](const auto& p) { return p.first == cand.compartment; });

        int streak = (it != qualifying_.end()) ? it->second + 1 : 1;

        if (streak >= criterion_.dwell) {
            // Promote! But check budget first.
            int totalParticles = 0, totalTiles = 0;
            for (const auto& a : active_) {
                totalParticles += a.particles;
                totalTiles += a.tiles;
            }

            if (totalParticles + cand.particles <= criterion_.particleBudget &&
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
            } else {
                // Budget exhausted, stop promoting
                result.problems.push_back("particle or tile budget exhausted");
                break;
            }
        } else {
            // Still building dwell streak
            stillQualifying.push_back({cand.compartment, streak});
        }
    }

    qualifying_ = std::move(stillQualifying);

    // --- Step 3: Compute totals ---

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

    // Compute wetted depth from volume and plan area, matching computeWaterDepth()
    // and estimateFlipCost()'s assumption of ~1000 particles per m³.
    const double planArea = (comp.bboxHi.x - comp.bboxLo.x) * (comp.bboxHi.y - comp.bboxLo.y);
    const double wettedDepth = (planArea > 0) ? comp.waterVolume / planArea : 0;

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

    // Seed particles: 2^3 lattice per cell
    const int particlesPerCell = 2;
    flip::seedBox(*field, lo, hi, particlesPerCell, kRhoSeawater);

    // Set exact mass from compartment water volume (m³ → kg)
    const double mass = comp.waterVolume * kRhoSeawater;
    flip::setTotalMass(*field, mass);

    // Initial velocity: particles at rest in body frame
    // TODO: add ship motion contribution when Ship::step() wiring is in place

    return field;
}

void demoteWater(Compartment& comp, std::unique_ptr<flip::Field> field) {
    if (!field) return;

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
