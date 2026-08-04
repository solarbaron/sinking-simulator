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

## Phase 1 — Engine skeleton ✅ *complete*

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
- ~~**Reflection and serialisation**~~ **done** — macro-based type description,
  schema-tolerant save format, and stable name-derived component ids so
  archetype ordering no longer depends on registration order
  (`01-architecture.md` §2). Save/load of a whole World still to come. across builds
- ~~**Multi-rate scheduler with time dilation**~~ **done** — per-system rates in
  simulation time, dilation bands that gate solvers a fast-forward cannot
  afford, catch-up ceilings, dependency levels run in parallel, and an integer
  nanosecond clock so the schedule never drifts (`01-architecture.md` §3)
- ~~**World save/load**~~ **done** — whole-world persistence with surviving
  entity handles, unknown-component tolerance and a byte-identical round trip
  (`01-architecture.md` §2)
- ~~**Offscreen verification harness**~~ **done** — dependency-free PNG codec,
  externally validated (`01-architecture.md` §4)
- ~~**Vulkan device**~~ **done** — shared instance/device/queue/allocation,
  graduated from the FEM spike, degrades cleanly with no GPU
  (`01-architecture.md` §4)
- ~~**Offscreen debug renderer**~~ **done** — colour+depth target, render pass,
  indexed triangle pipeline, PNG readback, all assertions closed-form
  (`01-architecture.md` §4)
- Render graph and bindless setup — deferred past the Phase 1 milestone, which
  the push-constant path already satisfies
- ✅ **Milestone met:** `tools/ferry_view` renders the ferry casualty in 3D with
  an orbiting camera and the compartments tinted by fill fraction, driven by the
  same `sim::Ship` and `core::Scheduler` the headless scenarios use. A cutaway
  rather than a wireframe — `fillModeNonSolid` is a device feature this build
  does not request, and clipping the starboard half away with the existing
  `sim::clipByPlane` shows the interior better anyway.

This is the least glamorous phase and skipping it is how projects like this die at
month 30.

### What Phase 1 actually taught

Every subsystem shipped green on its functional tests while still containing a
real defect. In each case a *different instrument* found it:

| Defect | Found by |
|---|---|
| Job records recycled while still queued | zero-worker configuration |
| Chase-Lev slot-reuse data race | ThreadSanitizer |
| Auto-grain clamp off by one | ThreadSanitizer as a *slow* operating point, not as a race detector |
| Arena aligning offsets instead of addresses | sweeping alignments 1–256 rather than the usual few |
| Arena overruns invisible to ASan | manual poisoning, verified with a deliberate negative control |
| Scheduler drifting over long runs | a 999-vs-1000 count that could have been "fixed" with a tolerance |
| `World::load` failing open, leaving a half-built world | feeding it every truncation of a valid save |

The pattern is consistent enough to be a rule: **a green functional test is
evidence the code does what you thought of, not that it is correct.** The
sanitizer builds, the adversarial orderings (far surface drawn second, destroy
from the front, load every truncation) and the closed-form expectations are what
actually found things.

Three test expectations were themselves wrong — a projected quad area, a
bow-on camera angle, and a mirror check that would have passed on two blank
frames. Each was fixed by deriving the expected value from geometry rather than
by loosening the assertion, which is the habit that makes the rest of the suite
worth trusting.

## Phase 2 — Sea and ship — *~12 em*

- ~~**Wave field**~~ **done** — JONSWAP/PM directional spectra, equal-energy
  binning, analytic evaluation with orbital kinematics (`02-simulation.md` §2)
- ~~**Buoyancy under a wavy free surface**~~ **done** — `integrateBelowSurface()`,
  the routine nonlinear Froude-Krylov needs (`01-architecture.md`)
- ~~**Ikeda viscous roll damping**~~ **done** — all five components, validated by
  roll decay and scale invariance (`02-simulation.md` §2)
- ~~**Propulsion, rudder, MMG manoeuvring**~~ **done** — blade-element propeller
  over four quadrants, Fujii rudder, KVLCC2 manoeuvring set (§7)
- ~~**Nonlinear Froude–Krylov**~~ **done** — `Sea` carries an optional `WaveField`;
  `Ship` integrates buoyancy against the actual surface and openings read the
  local elevation, so a breach under a crest floods faster than one under a trough
- ~~**Response amplitude operators**~~ **done** — `engine/sim/rao.{hpp,cpp}`
  measures RAOs the way a basin does: regular wave, discard the transient, fit a
  harmonic. Both asymptotes and *both* Froude–Krylov sinc zeros come out
  emergent (`02-simulation.md` §2). This is the instrument the milestone is
  defined against, so it had to exist before the thing it measures.
- ~~**Cummins state-space radiation at runtime**~~ **done, strip theory not BEM** —
  `engine/sim/radiation.{hpp,cpp}` and `Ship::attachRadiation()`. The offline BEM
  pipeline this line originally promised is *not* here and this is not a
  substitute for it: sectional coefficients come from an exact 2D close-fit source
  solve over Lewis forms, integrated by strip theory. `RadiationTable` is the only
  thing the Cummins machinery consumes, so replacing strip theory with a real
  panel code later changes how that table is filled and nothing else.
  Stations are derived from the hull mesh itself, so any ship asset gets a table
  (`02-simulation.md` §2)
- ~~**Ocean rendering**~~ **done** — a displaced grid driven by `sim::WaveField`
  itself, asserted against `elevation()` through the whole render path
  (`03-renderer-audio.md`)
- ~~**Hull rendering and basic materials**~~ **done** — `engine/gpu/hull.{hpp,cpp}`
  and `engine/gpu/material.{hpp,cpp}`. A lit solid on one pipeline with the sea, so
  the ship shares its depth buffer rather than being composited beside it; an
  analytic metallic-roughness BRDF; and a material set that is a **text file**, so a
  mod adds a surface without recompiling (`03-renderer-audio.md`,
  `05-data-modding-validation.md` §4). 0.26 ms of GPU time at 1080p for a whole ship
  and sea on the target card
- ~~**Propulsion coupled into the seakeeping path**~~ **done** — `Ship::propulsion`
  applies the MMG hull, propeller and rudder forces to the 6-DOF body, so a ship
  makes way under its own power. The **encounter frequency is emergent**: a hull
  translating through `A cos(kx − ωt + φ)` meets waves at `|ω − kU|` with nothing
  imposing it, measured against `encounterFrequency()` to 0.1%. `measureRaoAt()`
  runs the ship up to its own speed and steers to hold heading, because the
  reference hull is directionally unstable and departs into a turn after about
  1900 s of straight running (`02-simulation.md` §2)
- ~~**Ikeda coupled into the seakeeping path, and the roll added-mass frame**~~
  **done** — `Ship::attachRollDamping()` replaces `zetaRoll` with B44 re-evaluated
  every tick at the amplitude, frequency and speed the ship is actually at (52 ns
  a call against a 16 µs tick, so caching it would only buy a stale answer). The
  free-decay log decrement matches `2πζ/√(1−ζ²)` to 0.25% cycle by cycle, and the
  decrement *changes* through a decay exactly as an amplitude-dependent
  coefficient says it must, where the stand-in it replaced held it constant.
  `A_inf`'s rotational block is now referred from the baseline origin to the
  centre of gravity — a factor of 3.5 in roll on a ferry-like hull, with the sign
  pinned by a semicircular section whose roll added mass is analytically zero.
  Along the way: the angular stiffness scaling every modal damper was being
  measured about the wrong point, so `zetaRoll = 0.08` was delivering 0.144
  (`02-simulation.md` §2)
- **Milestone:** drive a ship in a real seaway. Validated against published RAOs.
  **Half met.** *Drive a ship in a real seaway* — yes: a hull under power, in a
  directional spectral sea, responding through nonlinear Froude–Krylov with
  radiation memory, drawn in the sea it is actually responding to, holding its
  heading. That works today.

  *Validated against published RAOs* — not yet, and what remains is now an
  **input** problem rather than a capability one. Comparison needs (a) a real
  hull's offsets, so the mesh is the benchmark ship and not an approximation of
  it, and (b) the published curves themselves. The S-175 containership is the
  obvious target; its RAOs are tabulated at Fn ≈ 0.275, which the machinery can
  now reach in principle but which the reference barge's propeller cannot. Sizing
  a powerplant for a real hull and sourcing the benchmark data are the two tasks
  left, and neither is a physics gap.

  `tools/seaway_view` is the visible form of the half that works: a named ship
  built from principal particulars, under power at a chosen Froude number, in a
  directional spectral sea, drawn in the water it is responding to. **Fn 0.275 is
  reachable** — measured, the S-175 form settles at Fn 0.081 at 1.4 rev/s, 0.241
  at 4.0 and 0.277 at 4.6 — so the benchmark *condition* is no longer the
  obstacle either. It runs in `verify.sh full`, because it is the only thing that
  exercises hull generation, radiation, propulsion, the wave field and both
  renderers against each other rather than one at a time.

  Three things a look at its output showed that no unit test would have:

  - ~~**The ocean patch edge is visible against the sky.**~~ **Fixed by the
    geometric cascade** (`03-renderer-audio.md`, "The cascade"). A finite grid
    ends, and at 1050 m it ended inside the frame; a bigger uniform grid was never
    the answer, because resolution is set by the *shortest* component and widening
    the patch squares the cost. Rings that halve their resolution each time they
    double their extent make reach exponential in the level count instead: the sea
    now reaches 67 km — the distance at which its edge is one pixel from the
    horizon for a 69 m eye at 720p — for **46 ms** of displacement against the old
    patch's 55 ms at 525 m. What it did *not* fix is the near field, which still
    costs what it always did and is still the FFT's job.
  - **The area curve has no parallel middle body.** `f(u) = 1 − (1−e)|u|ⁿ` is
    smooth everywhere, so a fine hull comes out canoe-like where a real
    containership has a long constant midship. It matters for appearance and for
    sectional-area-driven quantities; it does not much affect Cb, Cp or LCB,
    which is why the coefficient tests are all green. A three-parameter curve with
    an explicit parallel middle body is the fix.
  - **Irregular frequencies get worse as hulls get fuller.** The S-175 needs 8 of
    600 station-frequency solves repaired; the KVLCC2, with a nearly rectangular
    midship at Cm 0.998, needs 76 of 600. That is 12.7%, which is the region
    `radiation.hpp` warns about, and it is an argument for the extended integral
    equation rather than the interpolation currently standing in for it.

**Performance, now measured rather than projected.** The wave-field query is
essentially 100% of a wavy tick. Evaluating the surface once per hull *vertex*
rather than once per triangle *corner* removed a 6× redundancy for no change in
any answer, taking a 128-component sea from 52% to 28% of a 100 Hz budget and a
576-component sea from 225% to 119%. A vectorised sincos is still worth ~4× and
is the next step for large spectra — it is simply no longer the *first* step,
which is what extrapolating from a per-component figure had suggested.

At the end of Phase 2 there is a ship simulator, without the damage.

## Phase 3 — Structure — *~18 em*

The longest and highest-risk phase.

- ~~Structural mesh generation from scantlings~~ — done: `engine/sim/scantlings.hpp`,
  measured in `02-simulation.md` §3. Stiffeners are discrete rather than smeared,
  and the structural mesh is independent of the hydrodynamic hull mesh; both
  decisions are recorded there with what they rule out.
- ~~**Tier-0 beam**~~ **done** — `engine/sim/girder.{hpp,cpp}` and
  `collapse.{hpp,cpp}`: the hull girder as a free beam balanced on the wave, first
  yield, buckling, and Smith's-method progressive collapse swept along the length
  (`02-simulation.md` §3). Craig–Bampton reduction and the Tier-1 reduced model are
  **not** done, and are the piece between this and Tier-2
- ~~**Solid-shell elements for plating**~~ **done** —
  `engine/sim/solid_shell.{hpp,cpp}`. Explicit tet FEM for genuinely 3D regions
  already existed from the Phase 0 spike (`engine/sim/fem.{hpp,cpp}`). **What does
  not exist is a consumer**: nothing builds elements over a `StructuralMesh` and
  solves, so the elements are validated and unused
  (see `07-fem-spike-findings.md` §4 for why this split is not optional)
- ~~Co-rotational elasticity and J2 plasticity~~ — done: `engine/sim/solid_shell.hpp`
  and `engine/sim/plasticity.hpp`, measured in `02-simulation.md` §3. Radial return
  with isotropic hardening (kinematic available, defaulted off for want of a
  measurement); rate dependence deliberately deferred, with the reason recorded
- Adaptive zone promotion/demotion and interface coupling — **the largest thing
  outstanding.** It is a cost problem as much as an engineering one: measured, an
  elastoplastic solid-shell element costs 7.3 µs against 273 ns elastic, so a
  200 m² collision zone is ~2 hours of wall time per simulated second on 24
  threads. Promotion has to be rare and small, which is the whole reason the tier
  structure exists
- Ductile damage — done, `engine/sim/plasticity.hpp`: equivalent plastic strain to
  failure, regularised against the element's own size and against triaxiality.
  **Mesh-splitting fracture is not**: a failed integration point is deleted, and the
  maximum principal direction it returns — the plane a tear would open on — has no
  consumer yet
- GPU element solver — not started for solid-shell. The Phase 0 spike has a
  Vulkan compute back-end for tets (`07-fem-spike-findings.md`), which is the
  pattern to follow
- ~~FEM → flooding coupling (a tear becomes an orifice)~~ — done:
  `engine/sim/breach.hpp`, measured in `02-simulation.md` §3. Failed panels become
  merged `Opening`s whose area, position and connectivity come from the structure;
  what supplies the failed panels is still to come
- ~~Ship-to-ship contact detection and rigid-body response~~ — done:
  `engine/sim/collision.hpp`, measured in `02-simulation.md` §2. Two hulls
  interpenetrate as an exact solid; the force follows the overlap volume and is
  reported as a patch, a pressure and an energy rather than as a velocity change,
  which is the load case the FEM-active zone needs. The rigid half of the
  milestone; what it does *not* do is deform or tear
- Deformation and tear rendering — not started. `engine/gpu/hull.{hpp,cpp}` draws
  an undeformed hull from `sim::Ship`; nothing feeds it a deformed one
- **Milestone:** ram the ferry. The hull deforms, tears where the stress says it
  should, and the resulting hole floods at a rate the hole's own area determines.

  **Met, and `tools/ram_view` is it.** A 5175 t hull strikes the 8984 t ferry
  abeam; the penetration volume gives a force and the energy that went into the
  meeting; that energy is spent outward through her plating until bays tear; the
  torn bays become openings in the flooding network; and she floods at
  `Cd·A·√(2Δp/ρ)` through them. Not a step of that chain is reimplemented in the
  tool, and it runs in `verify.sh full` because it is the only thing that
  exercises collision, indentation, breach and flooding *against each other*.

  At 6 m/s: 1.14 s in contact, 149 MN peak, 74 MJ absorbed, 66 bays torn,
  112.8 m² of hole, 66.75 m² of it reaching a compartment, 2210 t of water and
  she is lost with GM −1.62 m.

  Two things fall out of it rather than being put in. **Beyond a threshold the
  outcome stops caring how big the hole is** — from 1.5 to 6 m/s the hole grows
  from 3.4 to 113 m² and the floodwater barely moves, because a small breach fills
  the compartment behind it inside 900 s just as a large one does. Damage
  stability is decided by *which* compartments open, which is what the subdivision
  rules are written around. And **where she is struck decides how she dies**:
  amidships she takes 2100 t and lolls 6°, at the quarter she takes 8100 t and goes
  over 26° the other way.

  She is lost in every case tried, and that is consistent rather than suspicious:
  `ram_view` applies no damage control at all, which is Phase 0's `none` scenario,
  and `none` loses her too. The scenarios that let her live are the ones where
  somebody acts.

  **What the milestone does not yet include.** The striking ship is rigid, so
  every joule goes into tearing the struck plating — `indentation.hpp` records that
  this is tight at low energy and loose above a few tens of megajoules. The
  deformation is computed by a membrane model rather than by the solid-shell
  elements and plasticity that now exist; wiring those in is the Tier-2 zone
  solver, and it is the largest thing still outstanding in Phase 3.

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

**Pulled forward:** the reading half of the ship definition format landed early,
because the ferry being C++ was blocking anything that needed a second ship.
`engine/sim/shipfile.{hpp,cpp}` builds a whole `sim::Ship` from text,
`ships/ferry.ship` is the prototype ferry in it, and `./shipsim --ship=<path>`
runs one. It matches `game::buildFerry()` to 1e-8 relative on displacement and
reaches the same 900 s outcomes in all three flooding scenarios. Still to do here:
the exporter, CSG subdivision from a shared plane set rather than a box per
compartment, and the importers. See `05-data-modding-validation.md` §1 for what
the format deliberately cannot yet say.

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
