// SPDX-License-Identifier: MIT
//
// Water promotion: quiescent ↔ dynamic escalation for compartment flooding.
//
// Follows the pattern established in `promotion.hpp` for structural zones (§1-5)
// and gas compartments (§6). A `WaterPromoter` decides which flooded compartments
// deserve `flip::Solver`'s resolved particle-based flow instead of the lumped
// `Compartment::waterVolume` model.
//
// **Why this exists:** A quiescent model treats water as a scalar volume with a
// flat free surface. That is right when the ship is at rest and wrong when she
// rolls: sloshing, green water, and the dynamic free-surface forces that capsized
// ships in beam seas are invisible to a volume. FLIP gives the answer, at 1-10
// core-seconds per simulated second per active compartment, so the decision is
// when to pay that cost.
//
// **The criterion (§2):** Ship motion (roll rate, lateral acceleration) crossed
// with geometric guards (minimum depth, minimum volume). Hysteresis and dwell
// prevent chatter, following `GasPromoter`'s structure exactly.
//
// **State transfer (§3):**
//   - Escalation: `Compartment::waterVolume` → `flip::Solver` particles via
//     `flip::seedBox` and `flip::setTotalMass`, exact mass conservation.
//   - Demotion: `flip::Solver::totalMass()` → `Compartment::waterVolume`,
//     `flip::quiescentLevel()` gives the surface height.
//
// **Budget (§4):** Particle count and tile count, shared across all active
// compartments. Unlike structural elements (per-zone) and gas cells (per-compartment),
// particles are the memory bottleneck and must be globally bounded.
//
// **Not here yet:**
//   - Integration into `Ship::step()` — that wiring is in `ship.cpp`.
//   - FLIP solver ownership — `Ship` will hold `std::map<int, flip::Solver*>`.
//   - Boundary conditions — FLIP sees AABB initially, TriMesh in Phase 2.
//   - Openings — quiescent network handles inter-compartment flow for now.

#pragma once

#include "flip.hpp"
#include "engine/core/math.hpp"
#include <string>
#include <vector>

namespace sim {

// Forward declarations to avoid circular dependencies
class Ship;
struct Compartment;

namespace promotion {

// --- 1. The criterion -----------------------------------------------------

struct WaterCriterion {
    // Ship motion thresholds: roll rate and lateral acceleration.
    // Roll rate triggers sloshing at the compartment's natural frequency;
    // lateral acceleration captures transient events (collision, hard rudder).
    double rollRatePromote = 0.05;      // rad/s — ~3°/s, excites typical sloshing
    double rollRateHold = 0.03;         // rad/s
    double accelPromote = 2.0;          // m/s² — ~0.2g, distinguishes casualty response
    double accelHold = 1.0;             // m/s²

    // Geometric guards: shallow puddles can't slosh, tiny volumes aren't worth it.
    //
    // **`minDepth` is not enforced yet.** The only depth available is the mean
    // over the compartment's whole bounding-box footprint, and on a ship whose
    // compartments are 24 x 8 m that reads ~0.02 m for volumes that stand well
    // over half a metre against a bulkhead under heel. Gating on it filtered
    // every compartment aboard and returned an empty candidate list from every
    // review. It stays here as the threshold the real free-surface calculation
    // will be compared against (`computeWaterDepth`'s TODO); until then
    // `minVolume` is the only geometric guard, on promotion and on hold alike.
    double minDepth = 0.5;              // m — restoring force goes as depth
    double minVolume = 1.0;             // m³ — FLIP has fixed overhead

    // Hysteresis: consecutive reviews a candidate must qualify for (dwell) or
    // an active compartment must fail (hold) before state changes.
    int dwell = 2;
    int hold = 3;

    // Budget: shared across all active compartments.
    // Particles determine memory (Vec3 x, Vec3 v, Mat3 C, double m = ~80 bytes).
    // Tiles scale with wetted volume (4×4×4 cells, sparse hash map).
    int particleBudget = 100000;        // total, all compartments
    int tileBudget = 2000;              // 4×4×4 tiles

    // Core-seconds per simulated second per compartment, for cost reporting.
    // Measured rather than assumed, on the same terms as structural/gas tiers.
    double coreSecondsPerCompartment = 5.0;  // estimate, will be measured

    // FLIP solver parameters for all promoted compartments.
    flip::Params solver;
};

struct WaterCandidate {
    int compartment = -1;               // index into Ship::compartments
    std::string name;
    double rollRate = 0;                // rad/s, abs(angularVelocity.x) in roll axis
    double accel = 0;                   // m/s², lateral acceleration magnitude
    double depth = 0;                   // m, current water depth
    double volume = 0;                  // m³, Compartment::waterVolume

    // How far past its own threshold the **weaker** of the two triggers is,
    // following `GasCandidate::score`. A score of 1.0 means both just cleared.
    double score = 0;

    int particles = 0;                  // estimated particle count for this compartment
    int tiles = 0;                      // estimated tile count
    double cost = 0;                    // core-seconds per simulated second
    std::string why;                    // human-readable reason (qualifies / rejected)
};

// Every tracked compartment that qualifies, ranked, before the budget is applied.
// Exposed so a caller can see what was considered and rejected (same pattern as
// `promotion::candidates()` and `promotion::gasCandidates()`).
std::vector<WaterCandidate> waterCandidates(const Ship& ship, const WaterCriterion& criterion);

// --- 2. The state machine -------------------------------------------------

struct WaterActive {
    int compartment = -1;
    std::string name;
    double rollRate = 0, accel = 0, depth = 0, volume = 0;
    double score = 0;
    int particles = 0, tiles = 0;
    double cost = 0;
    int promotedAtReview = 0;           // review number when this was promoted
    int idleReviews = 0;                // consecutive reviews below hold threshold
};

struct WaterReview {
    std::vector<WaterCandidate> considered;  // ranked, before budget
    std::vector<WaterActive> promoted;       // newly promoted this review
    std::vector<WaterActive> demoted;        // dropped this review
    int particlesActive = 0;                 // across all active compartments
    int tilesActive = 0;
    double costActive = 0;                   // core-seconds/sim-second total
    double microseconds = 0;                 // wall time for the review itself
    std::vector<std::string> problems;
};

// The decision, and nothing else: it does not build a FLIP solver, step one, or
// own one. That is `Ship`'s job, on the same terms as `GasPromoter` leaves the
// LES grid to `les.hpp` and the calling code.
class WaterPromoter {
public:
    explicit WaterPromoter(WaterCriterion criterion = {});

    // Review the ship's compartments and decide which deserve FLIP.
    // Does not mutate the ship, does not step solvers, does not own them.
    WaterReview review(const Ship& ship);

    const std::vector<WaterActive>& active() const { return active_; }
    const WaterCriterion& criterion() const { return criterion_; }

    int reviews() const { return reviews_; }
    int promotions() const { return promotions_; }  // cumulative
    int demotions() const { return demotions_; }    // cumulative

    void clear();

private:
    WaterCriterion criterion_;
    std::vector<WaterActive> active_;
    // Panel index -> consecutive reviews qualifying, ascending (for dwell).
    // Follows `GasPromoter::qualifying_` pattern exactly.
    std::vector<std::pair<int, int>> qualifying_;
    int reviews_ = 0, promotions_ = 0, demotions_ = 0;
};

// --- 3. State transfer helpers --------------------------------------------
//
// These are declared here but will be implemented in `ship.cpp` where they have
// access to both `Ship` internals and `flip::Solver`. They are the interface
// between the promoter (which decides) and the ship (which owns the solvers).

// Estimate particle and tile counts for a compartment, for budget purposes.
// Uses heuristic: ~1000 particles/m³, tiles from wetted volume / 64 cells per tile.
void estimateFlipCost(const Compartment& comp, const WaterCriterion& criterion,
                      int& particles, int& tiles);

// Create a FLIP solver for a compartment, seeded with its exact water mass.
// Returns a new solver owned by the caller (Ship will store in its activeWater_ map).
// Initial velocity set from ship motion at compartment centroid.
flip::Solver* promoteWater(const Compartment& comp, const Ship& ship,
                           const WaterCriterion& criterion);

// Read FLIP state back into compartment quiescent representation.
// Deletes the solver (caller must remove from activeWater_ map).
// Mass conservation is exact (not tolerance-based).
void demoteWater(Compartment& comp, flip::Solver* solver);

} // namespace promotion
} // namespace sim
