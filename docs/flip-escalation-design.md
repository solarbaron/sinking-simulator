# FLIP Escalation Design

## Status: DRAFT — awaiting agent findings

## Overview

Quiescent ↔ dynamic water escalation following the `promotion.hpp` pattern established for structural zones and gas compartments.

## Existing Infrastructure (Phase 5, completed)

- ✅ `flip::Solver` — sparse APIC solver, validated against closed forms
- ✅ `flip::seedBox` + `flip::setTotalMass` — exact mass seeding
- ✅ `flip::quiescentLevel` — mass → level conversion (round-trip identity)
- ✅ `Compartment::waterVolume` — quiescent water state

## Missing Pieces

1. **WaterPromoter** — criterion + state machine (following GasPromoter pattern)
2. **Wiring into Ship::step()** — where escalation checks run
3. **State transfer** — Compartment ↔ flip::Solver bidirectional
4. **Boundary conditions** — how FLIP sees compartment geometry
5. **Tests** — round-trip conservation, integration scenarios

## Design: WaterPromoter

Following `promotion.hpp` §6 structure:

```cpp
struct WaterCriterion {
    // Ship motion thresholds
    double rollRatePromote = 0.05;      // rad/s
    double rollRateHold = 0.03;
    double accelPromote = 2.0;          // m/s² lateral
    double accelHold = 1.0;
    
    // Water state thresholds
    double minDepth = 0.5;              // m — too shallow to slosh
    double minVolume = 1.0;             // m³ — not worth the cost
    
    // Hysteresis
    int dwell = 2;                      // reviews to qualify
    int hold = 3;                       // reviews to demote
    
    // Budget
    int particleBudget = 100000;        // total across all compartments
    int tileBudget = 2000;              // 4×4×4 tiles
    
    flip::Params solver;
};

struct WaterCandidate {
    int compartment = -1;
    std::string name;
    double rollRate = 0;                // rad/s
    double accel = 0;                   // m/s²
    double depth = 0;                   // m
    double volume = 0;                  // m³
    double score = 0;                   // weaker trigger excess
    int particles = 0;                  // estimated
    int tiles = 0;                      // estimated
    double cost = 0;                    // core-seconds/simulated-second
    std::string why;
};

class WaterPromoter {
public:
    explicit WaterPromoter(WaterCriterion criterion = {});
    
    WaterReview review(const Ship& ship);
    
    const std::vector<WaterActive>& active() const;
    // ... same structure as GasPromoter
};
```

## Criterion Rationale

### Ship Motion Triggers

**Roll rate** — captures resonant sloshing. Natural frequency ω = √(g k tanh(kd)) for water depth d. For typical compartment (10m wide, 2m deep): ω ≈ 1.4 rad/s, period ≈ 4.5s. Ship roll at 0.05 rad/s (3°/s) excites this.

**Lateral acceleration** — captures transient events (collision, hard rudder). 2.0 m/s² ≈ 0.2g distinguishes normal maneuvering from casualty response.

### Geometric Guards

**Minimum depth** — shallow puddles can't slosh (restoring force goes as depth). 0.5m threshold.

**Minimum volume** — cost floor. FLIP has fixed overhead (tile structure, CG solve). Below 1 m³, not worth it.

### Budget Strategy

Unlike structural/gas promoters, water has **particle count** as the primary cost driver:
- Particles determine memory (Vec3 x, Vec3 v, Mat3 C, double m per particle)
- Tiles scale with wetted volume (sparse, but CG solve cost)

Budget is **shared across all active compartments** — same as gas cells, opposite of structural elements (which are per-zone).

## State Transfer

### Escalation: Compartment → FLIP

```cpp
flip::Solver* promote(Compartment& comp, const WaterCriterion& crit) {
    auto* solver = new flip::Solver(crit.solver);
    
    // Geometry from compartment mesh
    AABB bounds = comp.mesh.bounds();
    solver->setSolid(bounds);  // TODO: actual mesh boundary
    
    // Seed particles with exact mass
    flip::seedBox(solver, bounds.lo, bounds.hi, crit.solver.dx);
    flip::setTotalMass(solver, comp.waterVolume * kRhoWater);
    
    // Initial velocity from ship motion
    Vec3 v_ship = ship.velocity + cross(ship.angularVelocity, comp.waterCentroid);
    solver->setInitialVelocity(v_ship);
    
    return solver;
}
```

### Demotion: FLIP → Compartment

```cpp
void demote(Compartment& comp, const flip::Solver* solver) {
    // Exact mass conservation (particle sum)
    double mass = solver->totalMass();
    comp.waterVolume = mass / kRhoWater;
    
    // Quiescent level from mass
    double z = flip::quiescentLevel(comp, mass);
    comp.surfaceOffset = z;  // body frame
    
    // Centroid from FLIP state
    Vec3 centroid = solver->centerOfMass();
    comp.waterCentroid = centroid;
    
    // Velocity discarded (quiescent assumption)
    
    delete solver;
}
```

### Conservation Assertion

```cpp
// In tests:
double m0 = comp.waterVolume * kRhoWater;
flip::Solver* s = promote(comp, criterion);
// ... run solver ...
demote(comp, s);
double m1 = comp.waterVolume * kRhoWater;
expectNear(m1, m0, 0.0, 0.0);  // exact, not tolerance
```

## Integration into Ship::step()

Following `Ship::solveFloodingNetwork()` pattern:

```cpp
void Ship::step(double dt, const Sea& sea) {
    // ... existing steps ...
    
    solveFloodingNetwork(dt, sea);  // quiescent transfers
    
    // Review escalation (not every tick — once per dwell period)
    if (ticks_ % waterPromoter_.criterion().dwell == 0) {
        WaterReview review = waterPromoter_.review(*this);
        applyWaterReview(review);
    }
    
    // Step active solvers
    for (auto& active : waterPromoter_.active()) {
        Compartment& comp = compartments[active.compartment];
        flip::Solver* solver = activeWater_[active.compartment];
        
        // Boundary conditions from ship acceleration
        Vec3 a_ship = /* ... from RK4 */;
        solver->setGravity(Vec3{0, 0, -kGravity} + a_ship);
        
        solver->step(dt);
        
        // Read back for coupling
        comp.waterVolume = solver->totalMass() / kRhoWater;
        comp.waterCentroid = solver->centerOfMass();
    }
    
    updateFreeSurfaces();
    // ... existing steps ...
}
```

## Boundary Conditions

**Ship motion in FLIP frame:**
- Gravity vector includes ship acceleration: `g_eff = g_world - a_ship`
- In ship body frame, gravity rotates with ship attitude
- FLIP works in inertial frame → transform or adapt formulation

**Compartment geometry:**
- Phase 1: AABB solid boundary (existing `flip::setSolid`)
- Phase 2: Actual TriMesh boundary (new `flip::setSolidMesh`)

**Openings:**
- FLIP doesn't see openings initially
- Quiescent network handles inter-compartment flow
- Future: source/sink terms in FLIP for breach inflow

## Testing Strategy

### Unit Tests (test_flip.cpp extensions)

1. **Round-trip identity**
   ```cpp
   testQuiescentToFlipAndBack() {
       // Seed 2.7183 m³ → FLIP → demote → 2.7183 m³ exactly
   }
   ```

2. **Mass conservation through motion**
   ```cpp
   testFlipConservesMassUnderShipMotion() {
       // Roll ship, step FLIP, mass unchanged to 0.0 tolerance
   }
   ```

3. **Criterion coverage**
   ```cpp
   testWaterPromoterDoesNotChatter() {
       // Oscillate at threshold, count promotions (should be 1)
   }
   ```

### Integration Tests (new test_ship_flip.cpp)

1. **Flooding compartment escalates**
2. **Calm water demotes**
3. **Budget enforced**
4. **Multiple compartments**

## Performance Targets

From `docs/06-roadmap.md`:
- Tier-0 beam: 0.10 core-s/sim-s (whole ship)
- Tier-1 reduction: 0.35–0.41 core-s/sim-s (one patch)
- Tier-2 FEM: 1155 core-s/sim-s (same patch)

FLIP cost (from flip_probe measurements):
- ~1000 particles/m³ seeding density
- 10 m³ compartment = 10k particles
- CG solve dominates (24 iterations for machine precision)
- Estimate: **1–10 core-s/sim-s per active compartment**

Budget implications:
- 5 compartments at once = 5–50 core-s/sim-s
- Compare to 10 Hz target = 0.1 sim-s per tick
- Affordable if promotions are rare and brief

## Open Questions (awaiting agent findings)

1. How does Ship currently detect motion (roll rate, accel)? Or derive from state?
2. What's the review cadence for structural/gas promoters?
3. Is there a Compartment index → flip::Solver* map already, or create one?
4. How are compartment meshes stored? TriMesh vs AABB?
5. Threading: can FLIP solvers run parallel per compartment?

## Next Steps

1. ✅ Understand promotion.hpp pattern (done)
2. 🔄 Wait for agent findings on Ship integration points
3. Draft WaterPromoter class skeleton
4. Wire into Ship::step() (no-op initially)
5. Implement state transfer functions
6. Add round-trip test
7. Implement criterion and budget
8. Integration tests
9. Validate against milestone scenario

---

**Estimated complexity: SIMPLE-MEDIUM** (1-2 weeks)
- Pattern established (GasPromoter template)
- Arithmetic proven (seedBox, setTotalMass, quiescentLevel)
- Main work is wiring and testing
