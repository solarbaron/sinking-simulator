# 05 — Ship Data, Modding, and Validation

## Why this is not an afterthought

A ship simulator with hardcoded ships has as many ships as someone had time to
hardcode. The current prototype's ferry is ~200 lines of C++ (`game/prototype/ferry.cpp`)
and it took real effort to get its compartment boxes to stay inside its own hull.
That does not scale to a fleet, and it certainly does not scale to a fleet built
by other people.

Equally: a simulator nobody has checked against reality is a physics-flavoured
toy. Validation is in this document, next to the data, because they are the same
problem — you cannot validate against a real ship without being able to *describe*
a real ship.

---

## 1. Ship definition format

Declarative, text, version-controlled, human-diffable. Reference structure:

```
ships/mv-example/
  ship.toml            identity, dimensions, loading conditions
  hull.offsets         station/waterline table, or a reference to hull.step
  structure.toml       scantlings: plate thickness, stiffener spacing, materials
  subdivision.toml     bulkhead and deck planes; compartments derived by CSG
  openings.toml        doors, hatches, vents, valves, pipes, with positions and areas
  systems/             electrical, bilge, ballast, fire main, fuel — as graphs
  machinery.toml       engines, propellers, rudders, steering gear
  loading/             loading conditions: cargo, ballast, consumables, KG, trim
  hydro/               precomputed BEM coefficients, windage tables (build artefacts)
  art/                 meshes, materials, sounds
```

**Compartments are derived, not authored.** You place bulkheads and decks as
planes or surfaces; the toolchain intersects them with the hull interior and
produces the compartment polyhedra. This is the fix for the prototype's biggest
data problem, it guarantees the compartments tile the hull without gaps or
overlaps, and it means changing the hull form does not silently invalidate every
compartment volume.

**Everything derived is cached and content-hashed.** BEM coefficient solves take
hours; windage sweeps take longer; FEM Craig–Bampton reduction takes minutes.
These are build artefacts keyed on the hash of their inputs, computed once, and
distributed with the ship.

## 2. Importing real ships

Real hull forms come as:

- **IGES / STEP** — the CAD interchange formats naval architects actually use.
  Import via OpenCASCADE, tessellate to the watertight mesh the sim needs.
- **Offset tables** — the classical station/waterline half-breadth table. Already
  supported by `makeHullFromStations()`.
- **DXF / lines plans** — traced body plans, for historical vessels where that is
  all that survives.
- **Point clouds / photogrammetry** — for scanned wrecks and preserved ships.

Public hull forms available for validation and as shipped content: DTMB 5415
(naval combatant), KCS (container ship), KVLCC2 (tanker), S175, JBC, Duisburg Test
Case, Wigley and Series 60 parametric hulls, and the various IMO/HARDER damage
stability test cases.

## 3. Editor

In-engine, because a ship editor that cannot run the ship is only half a tool.

- Hull: NURBS surface editing plus the offsets table view, with instantaneous
  hydrostatics feedback (displacement, LCB, waterplane, GZ curve) as you drag.
- Subdivision: place bulkheads, watch the compartments and the damage stability
  indices update live. The probabilistic damage stability index `A` from
  SOLAS Ch. II-1 is computable in seconds and is the single most informative
  number about a subdivision arrangement.
- Openings and systems: graph editors over the compartment topology.
- Scenario editor: initial conditions, sea state, damage events, scripted or
  free-form.
- Hot reload throughout — edit a bulkhead, see the ship's fate change.

## 4. Scripting and mods

- **WASM** for behavioural mods (systems logic, scenario scripts, AI, UI). Chosen
  over Lua for sandboxing, performance, and language choice. Mods cannot corrupt
  the sim state; they interact through a declared API.
- **Data mods** need no code at all — a new ship is a directory.
- Deterministic constraint: mods that touch simulation state run inside the
  deterministic boundary and are subject to the same rules (no wall clock, no
  unordered iteration, engine RNG only).

---

## 5. Validation

Every subsystem gets a validation target before it is called done. Not "it looks
right" — a number, compared against a number someone else measured.

### Already in place

`tests/test_core.cpp`, 29 checks, all against closed-form answers:

| Check | Reference |
|---|---|
| Closed-mesh volume and centroid | exact algebra on a box |
| Plane-clipped volume, axis-aligned and tilted | corner tetrahedron of a unit cube = 1/6 exactly |
| Volume→plane-offset solve round trip | self-consistency to 1e-6 relative |
| Floating draft of a homogeneous box | Archimedes |
| Metacentric height of a box barge | KB + BM − KG with BM = L·B³/12V |
| Free-surface loss of GM | ρ·i/Δ, matched within 15% |
| Trapped air arresting flooding | Boyle's law, pV conserved within 2% |
| Water mass conservation across a network | flow through the breach = water held |

### Planned, by subsystem

**Hydrostatics and damage stability**
- Published stability booklets for real vessels (displacement, KM, GZ curves at
  the stated loading conditions) — the ground truth naval architects themselves
  use.
- IMO/SOLAS Ch. II-1 probabilistic damage stability: compute the attained index
  `A` for a known ship and compare against the certified value.
- The **HARDER** project database of damage stability model tests.

**Seakeeping**
- ITTC benchmark cases; response amplitude operators (RAOs) for S175, DTMB 5415
  and KCS against published model-basin data across headings and frequencies.
- Roll decay tests: compare simulated decay envelope against measured, which is
  the direct check on the Ikeda damping model.
- Parametric roll: reproduce the known onset conditions for a container ship in
  head seas. A model that produces parametric roll for the right reason at the
  right frequency is a model that got the nonlinear restoring right.

**Structures**
- Standard FEA benchmarks first (NAFEMS suite) for the linear solver.
- Plate and stiffened panel collapse tests against published experimental
  load-shortening curves.
- Hull girder ultimate strength against progressive collapse analyses for known
  ships.
- Impact and tearing against published ship collision experiments and the
  large-scale grounding tests in the literature.
- Cross-check the Craig–Bampton reduced model against the full model it was
  reduced from — this one is free and should run in CI.

**Fire and thermal**
- The NIST/FDS verification and validation suite compartment-fire cases.
- ISO 9705 room corner test.
- Steel temperature rise in a standard fire resistance test against the
  Eurocode 3 curves.

**Fluids**
- Sloshing model tests (there is extensive published data driven by LNG carrier
  design).
- Dam-break benchmarks for the FLIP solver.

**The integration test that matters**

Reconstruct real casualties from the official investigation reports and check
whether the simulation kills the ship the same way and on roughly the same
timeline:

- **Herald of Free Enterprise** (1987) — bow doors open, water on the vehicle
  deck, capsize in ~90 seconds. Tests exactly the mechanism the current prototype
  already reproduces.
- **Estonia** (1994) — bow visor failure, progressive vehicle deck flooding,
  rapid list then capsize. Tests structure→flooding coupling.
- **Costa Concordia** (2012) — grounding, long raking damage, multiple
  compartments, slow list over hours, eventual capsize onto a shoal. Tests
  grounding, progressive flooding across many compartments, long timescales.
- **Sewol** (2014) — cargo shift under a turn with inadequate GM. Tests cargo
  shift and marginal intact stability.
- **El Faro** (2015) — flooding plus loss of propulsion in a hurricane. Tests
  seakeeping, machinery, and flooding together.

These are not scenarios to ship as entertainment. They are the most thoroughly
documented full-scale experiments in the field, and reproducing them is the
strongest available evidence that the coupled model is right.

### Continuous validation

All of the above runs in CI, not as a one-time exercise. Numerical results are
recorded per commit and regressions in *accuracy* are treated exactly like
regressions in *correctness* — the whole point of the project is a number that can
be trusted, and a number that silently drifts is worse than no number.
