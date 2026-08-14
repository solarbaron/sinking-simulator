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
#include <memory>
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
    //
    // **Both are read in the body frame**, which is not a detail. `state`'s
    // angular velocity and velocity are *world* vectors, and a ship lolled to 58
    // degrees has her own axes nowhere near the world's; a criterion that takes
    // the world x component as "roll" is measuring something else entirely at
    // exactly the attitudes this is for. See `computeRollRate`.
    //
    // The acceleration is differenced between reviews rather than taken from a
    // speed -- it was a speed once, compared against a threshold in m/s^2, and
    // the two sides of that comparison were never in the same units.
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
    //
    // **These two are not independent, and setting them as though they were made
    // one of them unreachable.** `estimateFlipCost` puts 1000 particles and
    // `1/(64 h³)` = 125 tiles in every cubic metre of water, so each budget
    // states a volume: 100 000 particles is 100 m³ and 2000 tiles is **16 m³**.
    // The tile budget therefore binds 6.25× sooner and the particle budget can
    // never be reached at all -- while the comment above calls particles the
    // memory bottleneck, which at 115 bytes a cell they are not: a tile is 7.4 kB
    // and a particle is 80 bytes, so tiles are 92% of the footprint.
    //
    // The consequence was not a tuning matter. Any compartment over 16 m³ was
    // refused outright, however hard she rolled, so the tier could not reach a
    // single compartment that matters: `water_probe` seeded the vehicle deck --
    // 462 m³, the largest free surface a ro-pax has and the reason this tier
    // exists -- and watched it be refused 154 times while a 1 m³ trickle in a
    // forward hold was promoted instead. The refusals were silent, because
    // nothing read `WaterReview::problems`.
    //
    // The old ceiling was ~97 MB and the vehicle deck would have been 2.8 GB, so
    // the *number* was defensible and the way it was expressed was not. Both are
    // kept, in agreement, and `maxVolumePerCompartment` states the real limit in
    // the units the decision is actually made in.
    //
    // **Those megabytes are at the solver's real seeding, not the estimator's.**
    // `estimateFlipCost` bills 1000 particles per m³; `flip::seedBox` at 2³ per
    // cell puts `8/h³` = 64 000 there, so a memory figure taken off the estimate
    // is 64× light. The estimate is what the *budget* is denominated in and it is
    // self-consistent, but any sentence about bytes has to use the other number.
    int particleBudget = 100000;        // total, all compartments == 100 m³
    int tileBudget = 12500;             // 4×4×4 tiles, == the same 100 m³

    // The largest single compartment worth promoting, m³. A compartment past this
    // is not refused silently as "over budget" -- it is rejected with a reason,
    // because a vehicle deck that cannot be afforded is a finding about this tier
    // and not a transient budget condition.
    //
    // 100 m³ at h = 0.05 is 12 500 tiles and 6.4 M particles — ~604 MB, which is
    // the largest single body of water worth resolving on one core.
    //
    // A larger compartment does **not** simply need a coarser grid. The vehicle
    // deck floods shallow and wide — 462 m³ over 1868 m² is 0.247 m deep — so the
    // h = 0.20 that would make it affordable puts 1.2 cells across the whole
    // depth, and `flip_probe`'s sloshing study wants ~2 cells of amplitude before
    // a voxelised surface has any restoring force. It would run, conserve mass
    // exactly, and not slosh. That compartment wants a depth-averaged model
    // rather than a bigger budget; see the roadmap's Phase 5 entry.
    double maxVolumePerCompartment = 100.0;

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
//
// `previousVelocity` and `dt` are what the lateral acceleration is differenced
// from. `dt <= 0` reports zero acceleration rather than inventing one from a
// single sample, which is the honest answer for a caller reviewing a static ship.
std::vector<WaterCandidate> waterCandidates(const Ship& ship, const WaterCriterion& criterion,
                                            const Vec3& previousVelocity = {}, double dt = 0);

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
    //
    // `dt` is the model time since the previous review, and it is what the
    // lateral acceleration is differenced over. A caller that passes nothing gets
    // a criterion on roll rate alone -- which is the right answer for a static
    // ship, where there is no acceleration to measure, and not a silent
    // degradation: `WaterCandidate::accel` reads 0 and says so.
    WaterReview review(const Ship& ship, double dt = 0);

    const std::vector<WaterActive>& active() const { return active_; }
    const WaterCriterion& criterion() const { return criterion_; }

    int reviews() const { return reviews_; }
    int promotions() const { return promotions_; }  // cumulative
    int demotions() const { return demotions_; }    // cumulative

    void clear();

private:
    WaterCriterion criterion_;
    std::vector<WaterActive> active_;
    // The velocity at the previous review, which is the other half of the
    // acceleration difference. Seeded on the first review, where there is no
    // earlier sample and therefore no acceleration to report.
    Vec3 previousVelocity_{};
    bool havePreviousVelocity_ = false;
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

// Create a FLIP field for a compartment, seeded with its exact water mass.
// Returns a new field owned by the caller (Ship will store in its activeWater_ map).
// Initial velocity set from ship motion at compartment centroid.
std::unique_ptr<flip::Field> promoteWater(const Compartment& comp, const Ship& ship,
                                          const WaterCriterion& criterion);

// Read FLIP state back into compartment quiescent representation.
// Takes ownership of the field and destroys it after reading state.
// Mass conservation is exact (not tolerance-based).
void demoteWater(Compartment& comp, std::unique_ptr<flip::Field> field);

} // namespace promotion
} // namespace sim
