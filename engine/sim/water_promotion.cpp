// SPDX-License-Identifier: MIT

#include "water_promotion.hpp"
#include "ship.hpp"
#include "engine/core/math.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>

namespace sim {
namespace promotion {

// --- Helpers --------------------------------------------------------------

static double computeRollRate(const Ship& ship) {
    // Roll is rotation about the longitudinal (x) axis in body frame.
    // Ship::state.angularVelocity is in world frame, so we need to transform
    // it to body frame to extract the roll component.
    // For now, assume small angles and take the x component magnitude.
    return std::abs(ship.state.angularVelocity.x);
}

static double computeLateralAccel(const Ship& ship) {
    // Lateral acceleration is transverse (y-axis in body frame).
    // For now, use world-frame velocity magnitude as a proxy.
    // TODO: proper body-frame transformation and time derivative.
    return length(ship.state.velocity);
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

std::vector<WaterCandidate> waterCandidates(const Ship& ship, const WaterCriterion& criterion) {
    std::vector<WaterCandidate> candidates;

    double rollRate = computeRollRate(ship);
    double accel = computeLateralAccel(ship);

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

WaterReview WaterPromoter::review(const Ship& ship) {
    auto t0 = std::chrono::steady_clock::now();

    WaterReview result;
    result.considered = waterCandidates(ship, criterion_);

    reviews_++;

    // --- Step 1: Check active compartments for demotion ---

    double rollRate = computeRollRate(ship);
    double accel = computeLateralAccel(ship);

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
    reviews_ = 0;
    promotions_ = 0;
    demotions_ = 0;
}

// --- 3. State transfer (stubs) --------------------------------------------
//
// These will be implemented in ship.cpp where Ship internals are accessible.
// For now, placeholders that compile.

flip::Solver* promoteWater(const Compartment& /*comp*/, const Ship& /*ship*/,
                           const WaterCriterion& /*criterion*/) {
    // TODO: implement in ship.cpp
    // - Create flip::Solver with criterion.solver params
    // - Get compartment AABB from comp.mesh.bounds()
    // - flip::seedBox(solver, aabb)
    // - flip::setTotalMass(solver, comp.waterVolume * 1000.0)  // kg
    // - Set initial velocity from ship motion at comp.waterCentroid
    return nullptr;
}

void demoteWater(Compartment& /*comp*/, flip::Solver* /*solver*/) {
    // TODO: implement in ship.cpp
    // - Read solver->totalMass() → comp.waterVolume (kg → m³)
    // - Read flip::quiescentLevel(solver) → comp.surfaceOffset
    // - Update comp.waterCentroid from particle centroid
    // - delete solver
}

} // namespace promotion
} // namespace sim
