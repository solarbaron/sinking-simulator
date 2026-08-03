# 06 — Roadmap

## Scale, honestly

This is a multi-year project. The subsystems described in `02-simulation.md` are
individually the subject of entire commercial products — potential-flow seakeeping
codes, explicit FEA solvers, and fire CFD are each sold as standalone software by
companies with staff. Building credible versions of all of them, coupling them,
and wrapping them in a game is not a six-month effort by any staffing level.

Estimates below are in **engineer-months** for someone competent in the relevant
domain, assuming the prior phases are done. They are not calendar time and they
are not padded for the discovery that always happens. Multiply by your own
confidence factor.

The order is chosen so that **something is playable at the end of every phase**,
and so that each phase's output validates the previous one.

---

## Phase 0 — Foundation ✅ *complete*

Numerical core, flooding, air, damaged stability, validation harness.

- Closed-mesh volume integration with plane clipping
- Compartments, two-phase orifice network, compressible trapped air
- 6-DOF rigid body with hydrostatics, added mass, measured damping
- Damaged GZ curves with floodwater re-levelling
- Mesh boolean (clip, weld, cap by ear clipping) so compartments are carved out
  of the hull form rather than authored as boxes
- Watertightness checking and load-time ship definition validation
- 116 closed-form validation checks
- A 120 m ferry that lolls over, capsizes, or survives depending on what you do
- Explicit co-rotational tet FEM, CPU reference plus a Vulkan compute back-end,
  validated against beam theory and benchmarked on the target GPU

**Deliverables:** `./build/shipsim`, `./build/shipsim_tests`, `./build/fem_spike`.
All run today.

---

## Phase 1 — Engine skeleton — *~8 em*

Everything in `01-architecture.md` that the current single-threaded prototype does
without.

- ~~Job system~~ **done** — work-stealing with task helping, deterministic
  reductions, verified under ThreadSanitizer (`01-architecture.md` §1)
- ~~**Arena allocators**~~ **done** — per-frame, per-lane bump allocators with
  ASan poisoning of the unused region (`01-architecture.md` §5). These are what
  make the Chase-Lev handle scheme below safe without a reclamation scheme.
- ~~**Job throughput and grain-scaling benchmark**~~ **done** — `tools/job_bench`.
  Dispatch costs 17–23 ns uncontended; efficiency plateaus at ~2 µs chunks; the
  penalty for bad grain is ~40×. Numbers and caveats in `01-architecture.md` §1.
- ~~**Grain auto-tuning**~~ **done** — `parallelForAuto()` probes for cost and
  targets 10 µs chunks, clamped to 2–64 chunks per lane. Matches or beats the
  best hand-picked grain without being told the element cost. Restricted to
  `parallelFor`; `parallelReduce` keeps an explicit grain because a
  timing-derived one would break bit-identical reductions
  (`01-architecture.md` §1).
- ~~**Chase-Lev revisit**~~ **cancelled on evidence.** Dispatch is 0.2% of a
  10 µs chunk, and the sweep shows no dispatch-limited regime at plateau grains,
  so a faster queue has nothing to recover. Reopen only if a profile shows queue
  contention. This is what the benchmark was for.
- ~~**Archetype ECS**~~ **done** — SoA chunk storage, generational entity
  handles, archetype transitions, chunked queries (`01-architecture.md` §2)
- Reflection and serialisation, save/load — also the fix for component ids
  currently being first-use-order rather than stable across builds
- Multi-rate scheduler with time dilation — continuation-style, so the steady
  state never calls `wait()`; this is what makes fibers unnecessary rather than
  merely deferred (`01-architecture.md` §1)
- Vulkan device, render graph, bindless setup, a debug renderer
- **Milestone:** watch the current ferry scenario in 3D, in real time, with a
  free camera and a wireframe of the flooding compartments.

This is the least glamorous phase and skipping it is how projects like this die at
month 30.

## Phase 2 — Sea and ship — *~12 em*

- Wave field (spectral, FFT, shared between physics and rendering)
- BEM offline pipeline + Cummins state-space radiation at runtime
- Nonlinear Froude–Krylov on the instantaneous wetted surface
- Ikeda viscous roll damping
- Propulsion, rudder, MMG manoeuvring
- Ocean rendering, hull rendering, basic materials
- **Milestone:** drive a ship in a real seaway. Validated against published RAOs.

At the end of Phase 2 there is a ship simulator, without the damage.

## Phase 3 — Structure — *~18 em*

The longest and highest-risk phase.

- Structural mesh generation from scantlings
- Craig–Bampton offline reduction; Tier-0 beam and Tier-1 reduced models
- Solid-shell elements for plating; explicit tet FEM for genuinely 3D regions
  (see `07-fem-spike-findings.md` §4 for why this split is not optional)
- Co-rotational elasticity and J2 plasticity
- Adaptive zone promotion/demotion and interface coupling
- Ductile damage and mesh-splitting fracture
- GPU element solver
- FEM → flooding coupling (a tear becomes an orifice)
- Deformation and tear rendering
- **Milestone:** ram the ferry. The hull deforms, tears where the stress says it
  should, and the resulting hole floods at a rate the hole's own area determines.

This is the phase the whole concept is named for. If it works, everything after
it is addition; if it does not, the project is a very good flooding simulator and
should be honest about that.

## Phase 4 — Fire and heat — *~10 em*

- Multi-zone compartment fire, species transport through the opening network
- Implicit thermal FEM on the structural mesh
- Temperature-dependent material strength (fire → FEM coupling)
- Suppression systems, and their effect on stability
- LES promotion for the local compartment
- Volumetric fire and smoke rendering
- **Milestone:** an engine room fire that heats a bulkhead until it fails under
  the head of water behind it, and the flooding spreads. Three subsystems, none
  of which know about each other, producing one consequence.

## Phase 5 — Fluids — *~9 em*

- Sparse FLIP/APIC solver for interior water
- Quiescent ↔ dynamic escalation with exact mass conservation
- SPH spray, jets from breaches, rain
- Sloshing, green water on deck
- Screen-space fluid rendering
- **Milestone:** stand in a flooding compartment as it fills around you, in a
  rolling ship, and have the water behave.

## Phase 6 — Crew and multiplayer — *~12 em*

- Avatar controller in a 6-DOF accelerating reference frame
- Interaction: doors, valves, pumps, extinguishers, tools
- Authoritative server, replication, prediction (`04-multiplayer.md`)
- NPC damage control parties and passenger evacuation
- **Milestone:** four people fight a flooding casualty together and lose an
  argument about whether to close the door.

## Phase 7 — Content pipeline — *~10 em*

- Ship definition format, CSG subdivision, importers (STEP/IGES/offsets)
- In-engine ship editor with live hydrostatics and SOLAS index
- Scenario editor
- WASM mod runtime
- **Milestone:** someone who is not you ships a ship.

## Phase 8 — Breadth — *ongoing*

The fleet and the long tail: cargo dynamics, ice and icing, grounding, collision
between two deformable ships, submarines, sailing vessels, machinery depth,
electrical and hydraulic networks, historical vessels, VR polish, audio
propagation.

This phase has no end and that is correct — it is where a simulator stops being a
tech demo and becomes a *ship* simulator covering "all kinds of ships and sizes
and utility".

---

## Critical path and risk

```
Phase 0 ✅ ──▶ Phase 1 ──┬──▶ Phase 2 ──┬──▶ Phase 3 ──┬──▶ Phase 4
                         │              │              └──▶ Phase 5
                         │              └──▶ Phase 6
                         └──▶ Phase 7 (can run in parallel from Phase 2)
```

**Highest risks, in order:**

1. **Phase 3 performance.** Adaptive FEM at interactive rates is the central
   technical bet. **Now measured** (`07-fem-spike-findings.md`): the GPU is fast
   enough, but only with the right element — uniform linear tets at
   plate-resolving resolution would have been 10³–10⁴× slower than real time.
   Residual risk moves to solid-shell element implementation, which is
   better-understood work than an open performance question. The tier structure
   still means a failure here degrades to Tier-1-only (elastic deformation,
   scripted-threshold tearing) rather than killing the project.
2. **Coupling stability.** Partitioned multiphysics can go unstable in ways
   neither solver does alone. Mitigation: every coupling gets an energy-balance
   check in CI; the flooding↔rigid-body predictor-corrector already in the code
   is the template.
3. **Scope.** The document you are reading describes more work than most
   commercial engines contain. Every phase must end in something playable, and
   if a phase overruns badly the correct response is to ship what exists rather
   than to press on.
4. **Validation data access.** Some of the best model-basin data is proprietary.
   Mitigation: the public benchmark hulls listed in `05-data-modding-validation.md`
   are enough to validate everything on the critical path.

## What to do next

**Both pre-Phase-1 items are done.**

1. ~~Compartment CSG~~ — done, and it found a hull winding bug that had been
   overstating displacement by 40%. See `README.md`.
2. ~~GPU tet-FEM spike~~ — done. Full results in
   **[07 — FEM spike findings](07-fem-spike-findings.md)**. Headline: the GPU
   delivers 450–670 M element-updates/s, the formulation validates against beam
   theory, and **uniform linear tetrahedra are ruled out for plating** — Tier 2
   needs solid-shell elements where the structure is thin. That change is now
   folded into `02-simulation.md` §3 and into the Phase 3 line item above.

So Phase 1 is next, with two lessons carried forward from Phase 0:

- **Every geometric representation gets a validity check that runs at load.** The
  winding bug cost nothing to fix and would have cost a great deal to find later,
  because the simulation kept producing believable numbers the whole time it was
  wrong.
- **Measure the expensive bet before building on it.** The spike cost a day and
  changed an architectural decision inside an 18-engineer-month phase, before any
  of that phase existed to be rewritten.
