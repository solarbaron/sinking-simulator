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
- 200 395 closed-form validation checks
- A 120 m ferry that lolls over or capsizes depending on what you do — and, since
  GM stopped being sampled at a fixed ±0.03 rad, does not survive any of the three
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
  (`02-simulation.md` §3)
- ~~**Solid-shell elements for plating**~~ **done** —
  `engine/sim/solid_shell.{hpp,cpp}`. Explicit tet FEM for genuinely 3D regions
  already existed from the Phase 0 spike (`engine/sim/fem.{hpp,cpp}`)
  (see `07-fem-spike-findings.md` §4 for why this split is not optional)
- ~~**A consumer for them**~~ **done** — `engine/sim/zone.{hpp,cpp}`, measured in
  `02-simulation.md` §3, run at ship scale by `tools/zone_probe` in `verify.sh
  full`. `StructuralMesh` + a load → solid-shell elements → an explicit solve →
  which panels tore, as indices `breachesFromFailedPanels()` takes unchanged. A
  three-metre zone on the ferry's side is 224 elements and a 0.037 s ram through it
  is 4.5 s of wall time on 23 workers against 15.4 s on one, bit-identical either
  way — the speedup saturates at 3.4× because each step is a barrier, not because
  the element kernel does not thread. Its
  binding limit is geometric rather than numerical: the elements are exactly
  prismatic on flat plating and 319% too stiff in bending across this hull's
  shoulder, and the cure is a finer girth layout in `Scantlings` rather than a finer
  zone. **Stiffeners were not meshed** — there is no way to attach a web to a
  solid-shell plate without a multi-point constraint — so the zone offered the two
  bounds that leaves and published the bracket
- ~~**The multi-point constraint**~~ **done** — `engine/sim/constraint.{hpp,cpp}`,
  measured in `02-simulation.md` §3 under *Eccentric stiffeners*. A solid-shell has
  no rotational degree of freedom and does not need one: its two nodes through the
  thickness *are* the rotation of the cross-section, so a member at offset `e` is
  the exact linear combination `u_bottom + ((e + t/2)/t)(u_top − u_bottom)` —
  exact for a finite rotation, measured at 3.6 × 10⁻¹⁵ m after 0.70 rad.

  `zone::Stiffeners::Modelled` builds a member out of axial fibres tied that way,
  two-point Gauss through each rectangle of the profile so the section's area,
  neutral axis and second moment come out of the finite element model **equal to
  `scantlings::stiffenedSection` to 2 × 10⁻¹⁰**, checked by sweeping the bending
  axis and reading all three off the parabola. It lands inside the published
  bracket (5.54 MN against 4.90 and 7.72), and because the member is *condensed*
  it adds no degrees of freedom and therefore no zero-energy mode — the stiffened
  patch has the same six as the bare one, measured on the whole spectrum.

  Two things it does not do, and one cost. It cannot represent tripping: the
  fibres contribute exactly zero to it and the plating alone restrains it, at the
  closed-form `16 D / b`, so the formulation *over*-restrains where the hinge
  leaves free. It carries no weak-axis second moment.
  And it re-introduces an in-plane length scale into a stable step that was
  thickness-governed: free at the ferry's resolution, 2.4× at a four-times finer
  one. **What it leaves ready**: the Tier-1/Tier-2 interface coupling below, which
  needs the same tie to drive a zone's boundary from a reduced model, and a
  whole-ship section mesh that has webs in it
- ~~**Craig–Bampton reduction**~~ **done** — `engine/sim/reduction.{hpp,cpp}`,
  measured in `02-simulation.md` §3. The missing middle: nothing could give a
  structural answer for a whole hold or the region between two bulkheads, because
  Tier 0 is a beam and Tier 2 is a patch. Boundary DOF kept exactly, the interior
  carried by constraint modes plus a handful of fixed-interface normal modes, and
  a small dense mass and stiffness pair out of it.

  Four properties are identities rather than tolerances and are asserted as such:
  zero modes reproduces static condensation and static condensation is **exact at
  the interface for any load** (2 × 10⁻¹⁰ m of a 0.31 m deflection against an
  independent solve); a free-free substructure keeps exactly six zero eigenvalues;
  the reduced frequencies come down **from above**, monotonically, because a
  reduction can only stiffen; and the reduced pair is `TᵀKT` and `TᵀMT`, formed
  the long way in the tests. Symmetric eigensolvers written rather than taken —
  dense Householder/QL, and subspace iteration whose mode count is verified by a
  **Sturm sequence** because a skipped mode is otherwise silent.

  **Three published figures were wrong and are corrected there**: the cost saving
  is 1200× rather than 10⁻⁵; "the static interface response improves with mode
  count" is false, it starts exact and the modes buy the *interior* and the
  dynamics; and the standard "cut off at twice the band of interest" buys 0.6%
  inside the 10 Hz hull-girder band, not four figures, because the cutoff is a
  frequency of the fixed-interface spectrum and the band is a frequency of the
  assembled one.

  Cost on the same plating: **0.10 core-seconds per simulated second** for Tier 0
  over the whole ship, **0.35–0.41** for this patch at Tier 1, **490** for the
  same patch at Tier 2. What it cannot do is yield, tear, buckle, contact or
  rotate — it is linear by construction — so `checkValidity` exists to say when the
  region has to be promoted, and it under-predicts a concentration so the warning
  is late rather than early. It reads the attached members' own stress as well as
  the elements', and reports the two halves apart: judging a stiffened region by
  its plating alone put the utilisation 11% low, in the unsafe direction

- ~~Co-rotational elasticity and J2 plasticity~~ — done: `engine/sim/solid_shell.hpp`
  and `engine/sim/plasticity.hpp`, measured in `02-simulation.md` §3. Radial return
  with isotropic hardening (kinematic available, defaulted off for want of a
  measurement); rate dependence deliberately deferred, with the reason recorded
- ~~Adaptive zone promotion/demotion and interface coupling~~ — **done for
  Tier 0 ↔ Tier 2**: `engine/sim/promotion.{hpp,cpp}`, measured in
  `02-simulation.md` §3. Three things that were three separate holes: a criterion
  that decides which patches deserve Tier 2 under an element budget, a
  `SectionReduction` that turns what tore back into a `StructuralMesh` Tier-0
  already knows how to read, and a `zone::Preload` that hands the zone the
  girder's own `M (z − z_na) / I` on the way in.

  It is a cost problem as much as an engineering one, and the cost is **linear in
  the number of zones** — measured, two zones take 2.00× one. So the criterion is
  built around *not* firing: a station qualifies only when it is past an absolute
  threshold **and** standing above the ship's own **median** for that trigger, so
  a hull at 0.9 of capacity everywhere promotes nothing (0 candidates with the
  guard, 31 without). A contact patch is the other trigger and is absolute,
  against the struck bay's own `4 σ_y (t/span)²` — 402 kPa on her side, bracketed
  at ±20%. Chatter is held off by hysteresis *and* a dwell, each tested against
  its own negative control: 1 promotion against 5 without the band, 0 against 20
  without the dwell.

  Two costs, three orders apart: the **decision** is 7 µs over 8900 panels and the
  **Tier-0 answer it reads** is 167 ms, of which 137 ms is the Smith sweep — so it
  is reviewed at a cadence and explicitly not every tick.

  The coupling to **Tier 0** still goes to a *section* rather than to retained
  interface DOF, and that is right for a beam: `hullGirderSection` reads a
  thickness and nothing else. The coupling to **Tier 1** now exists beside it —
  `engine/sim/coupling.{hpp,cpp}`, `02-simulation.md` §3 — and is exact: a zone's
  perimeter is driven from a Craig–Bampton model of the plating round it and a
  torn zone goes back as a mesh with the dead elements deleted. Against the same
  plate meshed and solved in one piece, the coupled zone reproduces the monolithic
  punch reaction exactly and its displacement field to **1.06e-15 m of a 2.0e-4 m
  peak**, at 0, 4 and 12 modes alike, where the clamped zone it replaces is **74%**
  out and 433× too stiff. The *mesher* that was missing now exists —
  `engine/sim/section.{hpp,cpp}`, `docs/02-simulation.md` §3 — and cuts a region
  between two transverse planes: a hold of the ferry is 2 068 elements and agrees
  with `hullGirderSection` to 0.3% on area and second moment once the girders it
  cannot attach are accounted for. It cannot **weld** a junction — `makeStructuralMesh`
  shares no corner between two panel roles and a solid-shell node pair carries one
  thickness direction where a corner has two — so it **ties** one instead:
  `solidshell::Mpc`, an eight-master constraint over the face the free edge lands in,
  eliminated rather than penalised. The ferry's hold goes from seven components to
  one, its `GJ` up 44.9% and its lowest fixed-interface frequency from 0.7785 Hz to
  2.3026 Hz — where untied it was the decks' own frequency to four figures, the shell
  contributing nothing. It costs the assembled band 146 → 1 520.
  The section
  reduction used to carry plating only, so a collision that opened fourteen bays left
  their longitudinals at full strength — the un-conservative direction, and the last
  named gap in Phase 3. **It is closed**: the fibres carry damage
  (`constraint.hpp` §2b, `02-simulation.md` §3 under *Eccentric stiffeners*), and
  `promotion::reactionOf`/`reduce` consume it, naming the member back to
  `StructuralMesh::members` and scaling its web and flange thicknesses — which moves
  area, `I_own` and the Steiner term together and leaves the centroid alone, the
  exact analogue of a plate's thickness. **What it was worth**, on the ferry, for a
  ram opening 125.6 m² of her side and the 163 m of longitudinal in it: the section
  area lost goes 6.861% → **8.455%**, the second moment 5.429% → **6.590%**, the
  hogging ultimate moment 5.711% → **7.068%** and the sagging one 11.846% →
  **17.507%**. So about a fifth of what a collision takes out of her hull girder was
  invisible, and about a third of it in sagging — stiffener loss bites hardest in
  compression, because a panel that has lost its stiffener buckles far earlier than
  one that has merely thinned. At the zone, against the same solve with the control
  `SolveParams::fiberFailure` off, the old model left the longitudinal at
  `ε_p = 0.513` — **2.83× its own failure strain** — still at full section, carrying
  2.23× the force and costing the ram **25.3% more energy to open the same hole**.

  The criterion is the plating's with two different arguments, and establishing them
  was the work: the averaging length is the **fibre's own rest length** (not the
  seam's — the tie extrapolates, so on a curved patch they differ), and the neck
  width is the **profile rectangle's thickness** (not the plate's). A bar's
  triaxiality is `±1/3` exactly, so the Rice–Tracey multiplier is **exactly 1** in
  tension and **infinite** in compression — a squashed stiffener never tears. That
  closed form is used rather than `plasticity::triaxiality` of a Voigt stress
  because `vonMises` of a uniaxial stress is *not* identically `|σ|` (measured: one
  ULP low in four of eight sampled magnitudes) and the compression branch is a `<=`
  on the cutoff, where one ULP is the difference between never failing and failing
  at 37% of the tensile strain

  **The "200 m² is two core-hours per simulated second" figure needs its element
  size beside it or it says nothing.** It assumes 50 mm elements; the same 200 m²
  at the 300 mm the zone actually uses is 2 200 elements and four minutes. What is
  bounded is the element *count*, and area and resolution trade as `area / h²` —
  the step does not care how big an element is in plane. The event is short too: a
  6 m/s bow reaching 0.22 m into her side is 0.037 s, not one second

- **A Tier-0 defect the coupling found, and fixed.** `longitudinalStrength` sized
  its progressive-collapse sweep from `firstYieldCurvature`, which is set by the
  *weakest* element in a section — safe only while none is anomalously weak, which
  stops being true exactly when it matters, because a damaged bay's critical
  stress falls as `t²`. Thinning forty of the ferry's side panels to an eighth
  dropped first yield 25× and the reported ultimate moment to 1.26e8 N m against a
  true 1.89e9; taking the same plating away *entirely* reported 1.87e9, so the
  hull girder got fifteen times stronger when material was removed.
  `collapseCurve()` now extends the sweep only when the peak lands on its last
  point, so an intact section's answer is **bit-identical** to before — her worst
  strength margin is 4.2366 either way
- Ductile damage — done, `engine/sim/plasticity.hpp`: equivalent plastic strain to
  failure, regularised against the element's own size and against triaxiality.
  **Mesh-splitting fracture is not**: a failed integration point is deleted, and the
  maximum principal direction it returns — the plane a tear would open on — has no
  consumer yet
- **GPU element solver for the solid-shell** — built, re-mapped, and **now the
  faster path on throughput while still not being usable, for a different reason**.
  `engine/gpu/zone_gpu.{hpp,cpp}` plus `solidshell_forces.comp`,
  `solidshell_forces_wg.comp` and `solidshell_integrate.comp`, following the tet
  back end's pattern: forces to per-element slots, a CSR gather in a fixed order,
  every substep in one command buffer. EAS internal variables and eight Gauss points
  of plastic history are device-resident and the seven enhanced parameters are
  condensed in the shader each step. The full account is
  `07-fem-spike-findings.md` §8; four things belong here.

  **The profile came first and it changed the job.** Element evaluation is 98.5% of
  a Tier-2 solve on one worker, so the pattern does apply — but half of that was
  `computeForms` rebuilding each element's strain-displacement matrices from its
  *rest* configuration every step, which an explicit solve never moves. Hoisting
  them into a `solidshell::RestForms` is **2.0× on the CPU for bit-identical
  answers**, at every zone size from 192 to 17 800 elements. The tet has always had
  this — it uploads `restInverse` and `restVolume` — and the solid-shell simply
  never grew the equivalent. The per-element cost *was* measured and was right; what
  was never asked is which part of it depended on the state being advanced.

  **One invocation per element was 0.23–0.68× the CPU; one workgroup per element is
  1.26–2.43×.** The first mapping got *worse* past 3 000 elements, which was
  diagnosed as register spilling from the shape of the curve alone. The driver's own
  statistics now confirm it — **1 936 bytes per thread of spill, 484 floats, against
  the tet's zero** — and the remap takes it to 96 bytes and halves the register
  count. The degradation past 3 000 elements is gone and the curve rises
  monotonically. **The tet's mapping does not carry over**: a linear tet has twelve
  DOF and no history and fits in registers, which is the whole reason it reaches
  450 M element-updates/s.

  **Float is still not sufficient, and the reason is the torn count.** Two defects
  were found and fixed — absolute ship coordinates in float lose the displacement to
  cancellation (solve about the patch centroid: 190× better), and Kaa is equilibrated
  before it is factored. `computeRestForms` now also normalises the enhanced modes,
  so κ(Kaa) is a constant 3.50 rather than (h/t)⁴; **measured by A/B, that changes
  nothing for this kernel**, because the shader's equilibration was already
  addressing the same conditioning. What remains: at 768 and 3 072 elements the float
  kernel tears 41 and 248 where the double reference tears **32 and 162**, while the
  negative control — the double solver on a mesh jittered by 2 × 10⁻⁷ m — tears
  exactly 32 and 162. Plastic dissipation runs 27–34% high where the control moves
  0.06–0.8%. **A previous claim here that it "tears 60 elements where the double
  reference tears none" does not reproduce under any configuration and has been
  withdrawn** — see §8; the GPU's peak damage on that case is 0.41 against the CPU's
  0.65, so it is further from tearing, not nearer.

  **`alpha` in double was the next step, it has been taken, and it is a negative.**
  §8 item 3 named it as the only part shown to need the digits. Five kernels are now
  compiled from one source, differing only in how much of the enhanced block is fp64
  (`--eas=float|tight|solve|condense|newton`). Over 5 505 steps they land between
  **40 and 44** torn at 768 elements and between **204 and 247** at 3 072, against a
  reference of 32 and 162 that the negative control hits exactly — a spread about half
  the gap, whose sign reverses between the two sizes, and which sorts by the *stopping
  rule* each kernel carries rather than by its arithmetic. Holding the stopping rule
  fixed, the whole block in fp64 moves alpha by 1/810 of the amount float is already
  wrong by, because Kaa's inputs are a float tangent and a float stress. It costs
  **5–15×** on the kernel, turning 1.26–2.43× into **0.09–0.47×**. **The list of cheap
  things to try is now empty**, and the error is
  where §8 said it was: per-step float in `u`, the return map and the integrator
- ~~FEM → flooding coupling (a tear becomes an orifice)~~ — done:
  `engine/sim/breach.hpp`, measured in `02-simulation.md` §3. Failed panels become
  merged `Opening`s whose area, position and connectivity come from the structure,
  and `zone.hpp` now supplies them from elements rather than from a script. What
  `breach.hpp` cannot yet take is a torn *area*, so a partly torn bay is rounded to
  a whole panel — see `SolveParams::tearFraction`
- ~~Ship-to-ship contact detection and rigid-body response~~ — done:
  `engine/sim/collision.hpp`, measured in `02-simulation.md` §2. Two hulls
  interpenetrate as an exact solid; the force follows the overlap volume and is
  reported as a patch, a pressure and an energy rather than as a velocity change,
  which is the load case the FEM-active zone needs. The rigid half of the
  milestone; what it does *not* do is deform or tear
- ~~Deformation and tear rendering~~ — done: `engine/gpu/damage.{hpp,cpp}` and
  `SceneMesh::appendShip`'s damaged overload, measured in `03-renderer-audio.md`.
  A damaged ship draws as damaged. `zone::Solver`'s displaced nodes become a
  piecewise-linear field over the patch's own elements — *interpolating*, so the
  drawn surface passes through every node to 2.8e-17 m — the hull is refined where
  the damage is because 5 m plating cannot show a 0.2 m dent, torn panels are
  **removed** so the pixel behind them shows the compartment or the sea, and the
  plating round a hole takes a material named in a `.materials` file.

  Two things are worth carrying forward. **Joining a deformed patch to an
  undeformed hull is the ocean cascade's problem in a triangle mesh** and it is
  solved the same way — a split decided per *edge* from its endpoints alone, an
  interned midpoint, and transition templates so a neighbour uses it — with the
  same instrument: an edge census, zero unmatched edges, and a control that leaks
  840. And **damage costs nothing until there is damage**: an undamaged hull comes
  back bit-identically, so the damaged and undamaged paths render byte-for-byte
  identical frames, which is asserted rather than intended.

  What it does not do is the plan's compute-shader skinning off a resident node
  buffer. `buildDamagedHull` is a 2 ms CPU rebuild that happens when the damage
  changes, not a per-frame step, and that is the right shape until something
  promotes and solves a zone *live* — which is the item above this one
- **Milestone:** ram the ferry. The hull deforms, tears where the stress says it
  should, and the resulting hole floods at a rate the hole's own area determines.

  **Met, and `tools/ram_view` is it.** A 5175 t hull strikes the 8984 t ferry
  abeam; the penetration volume gives a force and the energy that went into the
  meeting; that energy is spent outward through her plating until bays tear; the
  torn bays become openings in the flooding network; and she floods at
  `Cd·A·√(2Δp/ρ)` through them. Not a step of that chain is reimplemented in the
  tool, and it runs in `verify.sh full` because it is the only thing that
  exercises collision, indentation, breach and flooding *against each other*.

  At 6 m/s: 1.14 s in contact, 149 MN peak, 74 MJ absorbed, 0.205 m into her,
  63 bays torn, 107.7 m² of hole reaching four compartments, and she capsizes —
  19° of heel by 150 s, through 90° before 300 s, and lying at 165° with 16 613 t
  in her at 900 s.

  Those figures were first published with the indentation model spanning the
  frames at 2.40 m where the plating actually spans the longitudinals at 0.70 m.
  The zone FEM found it, having no span in it at all. Worth being exact about what
  it cost, because the natural reading of that finding is too broad: the energy to
  tear a bay is `σ_y·t·area·ε_f`, in which **the span cancels**, reaching the
  answer only through the failure-strain regularisation — and that is nearly flat
  here, 0.1517 against 0.1596. So the torn count moved from 66 to 63 and the 3 m/s
  case did not move at all. What was wrong by 3.4× was **force and penetration**:
  2.96 MN against 10.35, and 0.686 m of denting against 0.205, both in the
  direction of reporting the hull as far softer than it is.

  Two things fall out of it rather than being put in.

  **There is a saturation regime, and it is narrower than this section used to
  claim.** From 3 to 5.5 m/s the hole grows from 24 to about 95 m², four times
  over, and the floodwater moves 3%: 7182 t against 7415 t. A breach past a
  certain size fills the compartments behind it inside 900 s no faster than a
  smaller one does, so damage stability is decided by *which* compartments open —
  which is what the subdivision rules are written around.

  It is bounded on both sides, and the earlier draft of this paragraph had the
  lower bound wrong for a reason worth recording. It claimed saturation all the
  way down to 1.5 m/s, over a 32× range of hole size. That was an artefact of a
  defect: the mid wing tanks did not exist, so the extra hole area amidships
  opened onto no compartment at all and could not flood anything by construction.
  With the tanks in place, a strike at 1.5 m/s opens 3.4 m² and takes 2254 t,
  lolling her 9.4°, while 3 m/s takes 7182 t and puts her at 60.6° — the low end
  does not saturate, it simply had nowhere to put the water.

  **And she goes over.** Above 5.75 m/s the amidships strike stops being a loll
  and becomes a capsize: at 5.5 m/s she lolls to 63° with 7415 t, at 5.75 she
  rolls through to 167° with 16 425 t, and beyond that the inverted result is
  steady to four figures (16 613 t at 6 m/s, 16 627 at 6.5). The threshold is
  sharp, and 6 m/s — the speed this section publishes — sits about 4% above it.
  What carries her past the loll is the ram's own angular impulse: driven from
  rest with identical damage she settles at 44°, and the difference is roll energy
  the collision put in, not flooding.

  So **where she is struck decides how she dies**: amidships she goes over
  entirely, at the quarter she takes 7323 t and lolls to 48° the other way without
  inverting.

  Those figures are checked rather than quoted. `scripts/check-figures.sh` re-runs
  `ram_view` at each operating point and fails naming the line to update. It
  brackets the capsize rather than trusting a tonnage, because either side of that
  threshold a single number looks plausible on its own while describing a
  different casualty. The previous set went stale unnoticed twice — once when a
  wing tank authored *inside* a hold stopped flooding twice over, and once when
  the mid wing tanks were added — so the numbers in this section are now tested
  from outside it.

  She is lost in every case tried, and that is consistent rather than suspicious:
  `ram_view` applies no damage control at all, which is Phase 0's `none` scenario,
  and `none` loses her too. The scenarios that let her live are the ones where
  somebody acts.

  **And it draws.** `ram_view --frames=N --out=DIR` writes the flooding sequence
  with the damage in it: plating dished in, the torn bays cut out as holes with the
  ferry's own compartment meshes visible through them, and exposed metal round
  their edges. 27 631 triangles at 1280 × 720, 0.15 ms of GPU. The dent's *shape*
  is the membrane model's own tent kinematics, stated in the tool rather than
  implied — `indentation.hpp` reports a depth and a torn set and not a surface —
  and `gpu::HullDamage::addZone` takes the zone's displaced nodes instead when a
  run can afford the solve.

  Drawing it surfaced one thing about the ship that nothing else had: a
  `Compartment::mesh` is `clipToBox(hull, ...)`, so its outboard face is not near
  the shell, it **is** the shell, to the last bit. Drawn together the two z-fight
  over her whole side. The tool insets its copies by three per cent; the general
  answer is that an interior wants its own inboard surface rather than a copy of
  the shell, and that arrives with the compartment rendering in Phase 5.

  **What the milestone does not yet include.** The striking ship is rigid, so
  every joule goes into tearing the struck plating — `indentation.hpp` records that
  this is tight at low energy and loose above a few tens of megajoules. The
  deformation is computed by a membrane model rather than by the solid-shell
  elements and plasticity that now exist. Every piece needed to replace it is now
  here — `promotion.hpp` sites a zone from the contact patch, hands it the
  girder's stress, and takes its damage back — and `tools/zone_probe` runs that
  chain end to end on her side.

  ~~What `ram_view` still lacks is the *energy* bookkeeping: the zone takes a
  prescribed punch travel where a collision delivers a number of joules, and
  closing that needs the striking body's mass, which is `collision.hpp`'s to
  supply.~~ **The zone takes the joules now** — `zone.hpp` §6,
  `zone::Drive::Inertial`, and `zone::impactSpeed` turns `collision.hpp`'s own
  `ImpulseSolution::effectiveMass` and `energyLost` into the arrival speed that
  carries them. `zone_probe` drives the same collision both ways and prints the two
  answers against the membrane model's closed-form inverse.

  **Nothing in the figures above moves, and that is the point rather than a
  let-off.** `ram_view` spends its energy through `indentation.hpp`, which has been
  energy-driven since it was written; the gap was that the *Tier-2* zone had no
  such entry point, so the expensive model could not be asked the question the
  cheap one already answered. What the two now say about the ferry's own side, on
  the 2.755 MJ a 0.22 m punch costs: the membrane model on the longitudinal span
  tears the bay outright, and the FEM reaches the same 0.22 m with 11% of the
  striker's energy still unspent, tearing 42 elements where a punch driven at a
  constant 6 m/s tears 80. Replacing the tool's membrane with a zone is a
  core-minutes decision and is not made here; what is made is that the question can
  now be put to it.

  With the zone solver, the Tier-0 coupling and the Craig–Bampton reduction done,
  the largest thing outstanding in Phase 3 was **the Tier-1 coupling** — three
  pieces, none of which is the reduction itself, and all three are now built. A
  ~~**mesher** that can produce a hold-sized substructure, since the only one that
  exists is `zone::buildPatch`: plating only, within a radius, stopping at
  thickness seams~~ — **done**, `engine/sim/section.{hpp,cpp}`, with the junctions
  between panel roles **tied** rather than left open. An ~~**assembly** pass that
  matches two substructures' interfaces and stacks their modal blocks, which is
  what makes component mode synthesis a synthesis~~ — **done**, and it takes any
  number of components rather than two: `section::buildChain` cuts a length at N+1
  frame stations, reduces each piece once and assembles them, reproducing the same
  length in one piece to 1e-10 in `EA`, `EI` and `GJ`. And a ~~**two-way zone
  interface**, so a promoted zone is driven by the structure around it instead of
  clamped on a boundary that cannot move~~ — **done**,
  `engine/sim/coupling.{hpp,cpp}`.

  ~~What is left is **reach, not machinery**~~ — **the reach is done, and what it
  cost was a diagnosis rather than a rewrite.** Two-bay sections used to work over
  62.4 m of a 120 m ship in five disconnected islands — the longest unbroken run
  being 26.4 m, x = −7.2…19.2, which is the hold every figure here was measured on —
  with `buildSection` reporting an inverted element outside that. **Nothing was
  inverted**: every one of
  those refusals is a *collapsed* hexahedron — the triangular prism that one of the
  ferry's 166 degenerate `PlatePanel`s extrudes to — whose Jacobian is exactly zero
  at the closed edge and sound wherever the element is integrated, and
  `smallestJacobian` samples the corners. `section.hpp` §7 has the measurement; the
  reach is now the whole 120 m, and `./scripts/verify.sh full` runs
  `section_probe --scan=2` so it stays that way.

  ~~What is left is that a cut plane unties the junctions on it, which costs a
  five-section chain 8.9% of its torsional stiffness~~ — **done**, and the premise
  was half wrong: a prescribed degree of freedom cannot also be derived, but an
  *interior* cut plane is **shared** rather than prescribed. The junctions on it are
  now tied to the line the other surface draws in the same plane — every master a
  boundary DOF both sections have — and the constraint is applied to the assembled
  model rather than inside a section, where condensation being exact makes before and
  after the same number. The five-section chain's 8.9% becomes 0.17% and the box
  girder's 19.2% at eight sections becomes 1.6e-11; `docs/02-simulation.md`
  §*The in-plane line tie* has the measurements, and the old table is the negative
  control `SectionParams::interfaceTies = false` still reproduces. It was invisible
  in `EA` and `EI` at every stage, which is exactly the trap that section is written
  around.

This is the phase the whole concept is named for. If it works, everything after
it is addition; if it does not, the project is a very good flooding simulator and
should be honest about that.

## Phase 4 — Fire and heat — *~10 em*

- ✅ Multi-zone compartment fire, species transport through the opening network —
  `engine/sim/fire.{hpp,cpp}`. A hot upper layer over a cool lower one per
  compartment, Heskestad's plume entraining the one into the other, a design
  fire as the input, and gas crossing the openings `ship.hpp` already carries.
  Layer state is mass and *internal energy*, which makes the closure exact:
  both layers share one pressure, so `p = (γ−1)U/V` and therefore
  `V_u/V = U_u/U` — the volume split **is** the internal-energy split, with no
  root find. Validated against Heskestad's two entrainment branches meeting at
  the flame tip, the classical doorway integral and its
  `(T_∞/T_h)^(1/3)` neutral plane, MQH, Thomas, and an account that closes to
  1e-14 of scale on the ferry under a 4 MW fire.

  **Two corrections to what the flooding network suggested would carry over.**

  *A doorway is not one orifice.* `Ship::solveFlowNetwork` takes a single Δp at
  the orifice centre. With a hot layer on one side that is not merely
  inaccurate: the steady state of a doorway is hot gas out of the top and cool
  air into the bottom in **equal** mass flows, so a model that sees only the net
  is at rest and moves nothing, and the compartment behind it heats without
  bound. The vent integral therefore runs over height, split at each layer
  interface and at the neutral plane, closed form in each band. Measured on the
  ferry: 2.4 kg/s exchanged through the watertight door, 0.004 kg/s through a
  wing tank's air pipe, and the shortfall a single Δp reports is not a constant
  factor — it runs from 14% to a factor of nine depending on where the orifice
  centre happens to fall relative to the neutral plane.

  *The equalisation clamp cannot be repaired.* Clamping a transfer to the mass
  that equalises the two pressures starves a burning compartment, because it
  ignores the pressure the fire adds over the same step — 500 kW in an ISO room
  sat 800 Pa above atmosphere and vented a twentieth of what it should. Adding
  the source term to the target fixes that and then rings, because scaling one
  direction of a bidirectional flow to hit a net target is not continuous in the
  state: the steady state became a limit cycle between −3 Pa and −590 Pa. The
  clamp is replaced by what it approximates — the compartment pressures are
  solved **implicitly**, Gauss-Seidel over a monotone 1-D bisection per
  compartment. A 1.6 m² doorway relaxes a 100 Pa imbalance in 0.65 ms; unlike
  conduction, the gas really is sub-millisecond and explicit was never viable.

  Not here, and named rather than hidden: buoyancy-driven exchange through
  *horizontal* openings (a hatch in a deck has no Δp at all), the stack effect
  of the 5.5 m of air pipe between a deckhead and its gooseneck, a heeled
  compartment's layer geometry, and the `p dV` term a compartment that floods
  while it burns would need. `Ship`'s own gas is isothermal at `kTAmbient` by
  construction, so the write-back is a pressure proxy until `ship.cpp` carries a
  gas temperature.
- ✅ Implicit thermal FEM on the structural mesh — `engine/sim/thermal.{hpp,cpp}`.
  Backward Euler on the same solid-shell mesh and the same `BandedSpd`
  factorisation the statics use, one scalar unknown per node; Dirichlet, flux and
  convective boundaries; EN 1993-1-2 carbon steel with `k(T)` and `c(T)` live,
  closed by Picard on a **secant** heat capacity. Validated against the
  semi-infinite `erf`, the steady plate, the Kirchhoff transform for the
  nonlinear steady operator, and an energy account that closes to 5e-15 of the
  enthalpy moved.

  **A correction, because this item's own justification was wrong.** "Implicit"
  is usually argued for from an explicit stability limit of `rho c h^2 / 2k`
  said to be *milliseconds* on ship plating. It is **seconds** — 4.66 s for the
  ferry's 12 mm AH36, because steel's diffusivity is 1.5e-5 m²/s and 12 mm is
  not a small length. Milliseconds would need 0.3 mm elements. Explicit
  conduction is therefore *viable* on the unrefined structural mesh, and the
  real case for implicit is the `h²`: four elements through the same plate —
  which is what a through-thickness gradient needs, and a through-thickness
  gradient is what bows a plate — takes the limit to 0.29 s, and a 1.5 mm
  surface layer to 0.073 s. `thermal.hpp` carries the measurement and
  `explicitLimit()` reproduces it as a closed form.
- ✅ Temperature-dependent material strength (fire → FEM coupling) —
  `thermal::carbonSteelReduction` and `thermal::atTemperature`. EN 1993-1-2 §3.2
  Table 3.1 with the standard's own linear interpolation: `k_y,θ` on the
  effective yield, `k_p,θ` on the proportional limit, `k_E,θ` on Young's modulus.
  `atTemperature` returns a reduced `StructuralMaterial` or `plasticity::Material`
  by value, so a hot ship is the same ship with a reduced material and nothing
  downstream changed. The standard's four-branch σ(ε) curve, ellipse included, is
  there as the reference the reduced model is measured against.

  **The temperature-field question is settled by measurement, and the answer is
  per element.** Every element-level entry point in `solid_shell.hpp` and every
  material lookup in `buckling`, `collapse` and `indentation` already takes one
  material per element, so per-element temperature needed **no interface change at
  all**; building the reduced material costs 11.5 ns against
  `elementPlasticUpdate`'s 5.65 µs, 0.20%. Per-Gauss-point would need a signature
  change, and does not earn it: a 12 mm plate against a 900 °C compartment through
  a 200 W/(m²·K) film has a spread of **19.1 K across the whole plate** at its
  worst and under 1 K after half an hour, because the Biot number is 0.070. The
  worst spread in `k_y` that implies is 0.039.

  **A structural consequence that is not the yield factor.** The ferry's midship
  section loses **17.6% of its ultimate sagging moment at 400 °C, where `k_y` is
  still exactly 1** — plate buckling goes with `E`, and `k_E` is 0.70 there. At
  600 °C it keeps 0.376 against `k_y = 0.47`. A single 0.7 × 2.4 m × 12 mm panel
  changes regime outright: squash-governed cold (Johnson-Ostenfeld biting at
  213.6 MPa), elastic-buckling-governed at 600 °C (69.1 MPa, 0.323 of cold). The
  fully plastic moment does scale by exactly `k_y`, and is the closed-form anchor.

  **The return map needed nothing.** At any fixed temperature the scaled curve is
  still monotonically hardening — `k_y` multiplies `σ_y` and its slope alike — so
  radial return still lands on the unique root and step-independence still holds
  exactly. A yield surface that *shrinks* between steps is handled by the map it
  already has: the stored stress starts outside the new surface and is returned to
  it, which is stress relaxation.

  **What was missing was elongation, not strength, and it is four times larger.**
  Done, and the item below records what it turned out to be. Creep is *not*
  missing in the same sense: §3.2.1's curves carry it implicitly for 2–50 K/min,
  and both fires measured here (43 K/min under ISO 834, 114 K/min post-flashover)
  sit at or above that band, where the implicit treatment is conservative. A fire
  that stabilises and soaks is the case to revisit.
- ✅ **Thermal elongation and restraint (EN 1993-1-2 §3.4.1.1).**
  `thermal::carbonSteelElongation` and the eigenstrain path through
  `solidshell::elementStress`, `elementPlasticUpdate` and `thermalLoad`.

  **Where the strain enters, settled.** A thermal strain is an eigenstrain, not a
  load: `σ = C(ε_total − ε*)`, and the subtraction goes in the constitutive law —
  in both `elementStress` and `elementPlasticUpdate`, immediately before the
  material sees the strain. The equivalent nodal force `∫BᵀCε* dV` is a
  *consequence* of that and not an alternative to it; it is needed only where the
  constitutive law has been linearised away, which is the static solve.
  `elementPlasticUpdate` needs no load term at all — it already returns
  `−∫Bᵀσ dV` with the real `σ`. Doing only one of the two produces the same
  2.68 GPa fiction with opposite signs, on a bar under no load at all. A freely
  expanding body carries **exactly zero** stress, asserted at zero.

  **The 164.6 °C anchor reproduces: 164.630 °C**, `k_y` exactly 1 and `k_E`
  0.9354 — all expansion, no weakening — and it is right for E = 206 GPa and only
  for 206 GPa. At the 210 GPa `carbonSteelStress` defaults to it is 161.546 °C.

  **But it is the right answer to the wrong question, and buckling is why.** A
  restrained member goes into *compression*, and against `buckling.hpp`'s own
  checks on this ferry's own scantlings the plating goes long first: the 8 mm
  vehicle deck head in its 0.70 × 2.40 m bay at **59.0 °C — a 39 K rise**, the
  12 mm side shell at 103.1 °C, the 14.5 mm bottom at 121.1 °C, and only then the
  side longitudinal as an Euler column at 149.8 °C. Run through
  `collapseElementsAt` over the whole midship section, the ferry loses her first
  panel at **59.0 °C with no wave, no cargo and no bending moment at all**, at a
  temperature where `k_y` is exactly 1 and the panel has lost 0.0% of its own
  capacity to heat. Thermal expansion does not weaken the ship; it uses her up.

  **Two closed forms carry that.** On the elastic branch `k_E` **cancels exactly** —
  the restraint stress is `k_E E ε*` and the elastic buckling stress is `k_E E`
  times a pure geometry — so the buckling temperature is a function of geometry
  and the elongation curve *alone*, independent of E, of the grade and of the
  reduction factor. Asserted by changing all three. For a column that reads
  `ε*(T) = π²/λ²`, and equating it to the restrained yield strain gives a
  **critical slenderness of 73.19**: the ferry's longitudinals sit at 34–45, which
  is exactly why they pre-empt yield only through the plasticity cap and the
  plating pre-empts it by a factor of three in temperature rise.

  **What "restrained" means on a ship, because a fully restrained bar is a
  laboratory object.** A uniformly heated region expands freely and carries
  nothing; an axially heated chain with free ends is statically determinate and
  carries nothing either, however the temperature varies along it. Stress comes
  from a *transverse* gradient — a hot strake beside a cold one, sharing a seam.
  For two parallel strips tied at both ends, hot fraction *f*, the fully
  restrained figure is the *f* → 0 limit, and **the limit is approached fast**:
  181.5 °C at *f* = 0.1 against 164.6 at *f* = 0, and 313.6 °C only when half the
  section is hot. One compartment of a ferry is nearer a tenth than a half, so
  this is not a laboratory result.

  **A correction to the standard's own curve.** EN 1993-1-2 §3.4.1.1 is *not*
  continuous at both ends of the 750–860 °C phase change. It is exact at 860 —
  `2e-5·860 − 6.2e-3 == 1.1e-2` to the last bit in binary as well as in decimal —
  and it steps **down by exactly 8.4e-6 at 750 °C**, 0.076% of the elongation and
  1.73 MPa of restrained stress. Asserted rather than smoothed away.
- ✅ Suppression systems, and their effect on stability — `fire::Drencher`,
  `fire::Scupper`, `fire::scupperFlow`, `fire::sprayMassFlow`. A mass flow of
  water into a compartment, a share of it evaporating in the hot layer (which is
  where the cooling comes from) and the rest landing on the deck (which is where
  the stability cost comes from), plus the freeing ports that are supposed to get
  it back out — with the blocked case and the submerged case, because those are
  the two ways a port stops working.

  **Suppression water needed no path of its own, and establishing that was half
  the item.** `ship.cpp` contains no free-surface *correction* anywhere: it
  carries floodwater as real mass at the real centroid of the real water body and
  re-levels it at every attitude it is asked about, so the effect is emergent. A
  cubic metre added to `Compartment::waterVolume` therefore arrives with its
  mass, its centroid, its free surface and its cost to GM already right.
  Measured, on a box barge where the second moment is arithmetic and on the
  ferry's own vehicle deck where it is not: `GM_solid − GM_liquid = ρ μ I / Δ` to
  **1.4e-6 relative** on the box, which is the truncation of the finite-difference
  GM itself, and 2e-4 on the ferry, where `I` is itself a mesh integration.
  The permeability is in it and is not optional — dropping it overstates
  the loss by 11% on that deck. `test_core.cpp` already checked this at 15%; the
  slack there is the two barges not having quite the same KG, not the model.

  **The headline number is not GM, and that is the finding.** The ferry's
  undivided 100 × 19 m vehicle deck has a free-surface moment of 5.19e7 kg·m, so
  **ten tonnes of water takes her 2.00 m of GM negative** and everything after
  that sits at the same −3.7 m whether the deck is holding thirty tonnes or two
  thousand. GM saturates; the **angle of loll** is what tracks the water. Four
  hours of a SOLAS ro-ro drencher (5 L/(min·m²)) over one 20 m section:

  | Freeing ports | Water on deck | GM | Loll | Drained |
  |---|---|---|---|---|
  | clear, 6 × 1 m | **36 t** (steady) | −3.77 m | **0.8°** | 397 t |
  | blocked | **433 t** | −3.68 m | **9.7°** | 0 t |

  A twelve-fold difference in water moves GM by 0.09 m and the loll by a factor
  of twelve. The ports hold the deck at the weir's own equilibrium depth,
  `h = (q / ((2/3) C_d b √(2g)))^(2/3)`, verified as a fixed point to 1e-3.

  **A caveat on a published quantity — since repaired.**
  `Diagnostics::gmTransverse` was finite-differenced at a fixed ±0.03 rad, which
  is the right question for a ship at a finite angle and the wrong one for a
  shallow layer. With 50 t on the vehicle deck — a 2.9 cm layer, which spans the
  deck only out to 0.0031 rad — it read **+0.59 m where the initial GM is
  −3.77 m**, because the water has pulled off the high side long before ±0.03.
  The pocketing is real, and it is why she lolls rather than capsizing; but a
  positive GM where the truth is −3.77 m is the unsafe direction on the one
  number every stability judgement here keys off, so it was repaired rather than
  left standing. `Ship::diagnostics()` now halves the sampling angle until the
  slope stops moving and publishes the angle it settled at; the write-up, the
  wedge closed form and what the linear region actually costs are in
  `02-simulation.md` §1.

  **One published result moved, and it is a verdict rather than a figure.** GM
  feeds no physics — it is read, never integrated — so no tonnage, heel, draft or
  torn-panel count anywhere moves, and the 69 figures `scripts/check-figures.sh`
  gates are unchanged to the digit. What moved is `shipsim --scenario=full`,
  whose outcome line is decided by the sign of GM. Every figure below is measured
  at the end of the gate's own 900 s runs:

  | scenario | water on the deck | layer | pockets at | GM at ±0.03 | converged GM | verdict |
  |---|---|---|---|---|---|---|
  | `none`  | 119 t  | 6.9 cm | 7.4e-3 rad | −0.62 m | **−3.01 m** | LOST, unchanged |
  | `doors` | 3842 t | 2.23 m | 0.234 rad  | −2.51 m | **−2.51 m** | LOST, unchanged |
  | `full`  | 11 t   | 6.4 mm | 6.8e-4 rad | **+1.37 m** | **−3.30 m** | **SURVIVED → LOST** |

  `full` is the *successful* damage-control response — every action taken, 1442 t
  aboard, 7.2° of heel — and it has been reporting "SURVIVED but the deck edge is
  under" on the strength of a 6 mm puddle sampled at forty-four times the angle
  it spans the deck to. `doors` does not move at all, its water being deep enough
  that ±0.03 rad was inside the linear region all along, which is the control
  that says this is pocketing and not a global change of scale.

  **And a fixed smaller angle would not have fixed it.** ±0.001 rad — thirty
  times finer, and comfortably inside the 50 t case that started this — is still
  outside the linear region of a 6 mm layer, and reports −2.91 m for `full`
  against the converged −3.30 m. There is no fixed angle: the layer sets it.

  **The evaporated fraction, which is the coefficient nothing here measures.**
  Per kilogram, cooling runs 582 kJ at `e = 0.1` to 2387 kJ at 0.9 — a factor of
  **4.1**, not the factor of nine the latent heats alone suggest, because the
  sensible heat to saturation is a floor that no fraction goes below and is a
  third of the total at 0.3. The water landing runs the other way and by more
  (9×). **And on a real ship neither matters much**, because the fraction is a
  ceiling and not a rate: what evaporates is bounded by the power available to
  boil it, and the ferry's 20 MW fire converts **3.3%** of a 31 kg/s flow, not
  30%. The stability answer is set by the fire's power and the port area.

  **A correction found by measuring rather than by thinking, twice.** The first
  scupper compared `Compartment::surfaceOffset` against a sill *height* — but
  that field is the water plane's offset along the **body** up-vector, not a
  height, and this ship lolls within a minute of the drencher starting. The ports
  silently stopped draining at the moment they were needed and the run looked
  like a system with no freeing ports at all. A sill *position* and
  `surfaceWorldZ − sillWorldZ` is exact at any attitude and buys the physics the
  height could not: on a lolled deck only the low-side ports run, and on the
  ferry — a fifth of a degree of heel and trim apiece under a 2 cm layer —
  **two of six** carry the whole flow while the other four stand dry. Which is
  the real argument for distributing freeing ports along a deck edge, and it is
  emergent rather than assumed.

  The second: the cooling was first a constant per-kilogram sink with the
  *engaged mass* capped at what the layer could supply, which is the entrainment
  cap's idiom and is wrong here, because this cap binds at the **steady state**
  rather than on a transient — a drencher routinely out-absorbs the fire it is
  aimed at, by ten times on the ISO-room fixture. The applied rate became
  `available / dt` and the model stopped converging under refinement. It is now
  an exact relaxation towards the water's own temperature — the boundary loss's
  `expm1` treatment — which has a real equilibrium `Q_fire = Q_spray(T)` and is
  the same at any step.

  Not here, and named rather than hidden: the **steam mass** is not added to the
  gas, because `fire.hpp`'s closure needs one `γ` for both layers and injecting
  steam as air would be a 60% error in its own contribution to pressure — so the
  expansion that lets a sprinkler push smoke out of a compartment is missing;
  the design fire is still an input, so water cools the gas but **does not put
  the fire out**; a drenched compartment's gas volume is fixed at `attach()`, so
  the water on the deck does not shrink the gas space; and boundary cooling
  needed no new mechanism at all — it is `wallConductance` and `wallTemperature`,
  which `GasCompartment` already carries.
- ✅ **Suppression by gas: CO₂ and inert-gas total flooding** — `fire::AgentSpecies`,
  `fire::AgentSystem`, `fire::agentMassForFraction`, `fire::exposureAt`,
  `GasCompartment::settlingVelocity`, `DesignFire::limitingOxygen`. A second gas
  species with its own `R` and `γ`, carried in both layers, separating under
  gravity, displacing the oxygen the fire breathes, and leaking back out through
  the ship's own openings.

  **The closure generalises exactly rather than being approximated.** Per layer
  `p V_k = (m_air R_air + m_agent R_agent) T_k = (γ_k−1) U_k`, so
  `p = [(γ_u−1)U_u + (γ_l−1)U_l]/V` and the volume split is still the
  *pressure-energy* split with no root find. `Layer::excessEnergy` carries the
  whole difference — the extra internal energy a pure-air layer would need to
  press as hard — and it is **exactly 0.0** without an agent, so every figure
  already published is bit-identical. Asserted as bit-identity, not as closeness.

  **The agent stratifies downward, and that is the hazard.** CO₂ is M = 44 against
  air's 29. The discharge is delivered *mixed* (in proportion to the layers'
  moles: a total-flooding jet is what NFPA 12's two-minute limit exists to
  produce), and every bit of stratification is then a buoyant separation flux with
  a single coefficient, `settlingVelocity`. Two one-way streams, not a diffusive
  exchange, because a diffusive exchange drives towards uniformity and gravity
  drives towards segregation: the agent falls out of the layer it does not belong
  in, and the carrier gas is floated out of the layer it does, in proportion to
  the agent fraction there. Both are `expm1` relaxations, both are exactly zero
  without an agent, and the fixed point of the pair is the fully stratified state.
  **Both limits are closed forms and both are asserted:**

  | | concentration at the deck | O₂ at the deck | at the deckhead |
  |---|---|---|---|
  | `w = 0`, perfectly mixed | 40.0% | 12.57% | the same |
  | `w → ∞`, perfectly stratified | **100%** | **0%** | nearly air |

  In the stratified limit the blanket is pure agent under pure air at one
  temperature, so it occupies the agent's own mole fraction of the height: **2.40 m
  of a 6 m machinery space at 40%**, asserted as a length. The direction is a
  consequence of the molar mass and not a hard-coded drift — the same run with
  IG-100 (nitrogen, M = 28) puts the blanket under the *deckhead* and leaves a
  person on the deck breathing all of the oxygen.

  **The mass is not free, and the pressure is what loses you the concentration.**
  40% by volume means the space is holding `1/(1−y)` times the moles it started
  with, so a sealed 720 m³ machinery space at its design concentration sits
  **67.6 kPa above atmospheric** — seven metres of water head on every boundary —
  on 893.4 kg of CO₂. That is the term that drives the agent back out, and where
  the opening is decides everything. Thirty minutes after a 120 s discharge, with
  one 0.8 × 2.0 m door and nothing else different:

  | | agent lost | concentration held |
  |---|---|---|
  | sealed | **exactly 0 kg** | 40.0% |
  | door at the deck (sill 0 m) | **892.8 kg of 893.4** | **0.04%** |
  | the same door at the deckhead (sill 4 m) | 220.1 kg | **46.8%** |

  A door at the deck loses 99.9% of the charge. The *identical* door four metres
  higher loses a quarter of it and leaves the space **richer than it was aimed
  at**, because what a high vent sheds is the air standing over the blanket. A
  single Δp at the orifice centre cannot tell those two holes apart; the height
  integral the doorway needed for hot gas is what makes this a mechanism rather
  than an assumption, which is the second time in this file that has been the case.

  **Extinguishing is displacement.** The air the fire breathes is the air its plume
  entrains, so the two layers are weighted by `Plume::entrainment` — the
  correlation already there, and no new coefficient — and the design curve is
  scaled by a linear ramp between the fuel's limiting oxygen concentration and
  ambient. `X_O2 = 0.2095(1−y)`, so a fuel at `L = 0.13` is out at **y = 37.95%**
  and nowhere else, which is why the marine design concentration is 40 and not 30.
  The control beside it: the same 3 MW fire in the same space **with no agent is
  still burning ten minutes later**, at full oxygen availability and more than
  200 K over ambient, while the flooded one is out and cooling.

  **The phase change is modelled, and the number is large.** CO₂ leaves the bottle
  as a liquid and arrives as cold vapour and dry-ice snow at 194.65 K; the snow
  sublimes in the compartment and the compartment pays for it, at **constant
  volume** — `L_sub − R T_sub`, not `L_sub`, and using the constant-pressure
  enthalpy instead over-states the sink 2.5-fold. Per kilogram the gas receives
  **−6.31 kJ**: negative, so the discharge is a net internal-energy sink even
  counting the agent's own gas. On the gas alone that is **145.2 K and 85.1 kPa —
  *below* atmospheric despite a tonne of added gas**. Against a 25 W/(m²·K)
  boundary at ambient the same space comes all the way back to 288.15 K and
  168.9 kPa. Both are asserted. The gas-only chill is a **ceiling and not a
  prediction**: everything the space *contains* is missing as a thermal mass, and
  this model has only `wallConductance` against `wallTemperature`.

  **What it does to the people in the space, as data.** At 40% CO₂ the oxygen is
  12.57% — below the 19.5% entry limit, above the 10% that takes consciousness,
  and well above the 6% that kills. **On oxygen alone the atmosphere is
  survivable. It is lethal anyway**, because 40% CO₂ is four times CO₂'s own
  lethal concentration and CO₂ is not merely an asphyxiant. Stratified, the deck is
  inside a blanket at ~100% and there is no oxygen there at all while the gas at
  the deckhead is still nearly air — so the space's *mean* concentration describes
  neither. IG-541 at the same 40% leaves the same 12.57% oxygen and 3.2% CO₂,
  under the 4% IDLH: not lethal, and still not an atmosphere anyone may enter.
  That contrast is the whole reason the blend exists, and here it is emergent from
  one species record rather than asserted.

  **Two corrections, both found by measuring rather than by thinking.**

  *Freezing `γ` in the `p dV` split is wrong and it shows.* With a second species
  each layer's own `γ` moves as its composition does, and the obvious
  generalisation — the same elimination with `a` and `b` held constant — leaves a
  counter-current exchange between two layers at **one** temperature at two
  temperatures. Nothing had been heated. The corrected split carries
  `K_k = T_k(dR_k − a_k dC_k)` and is **exact** for the isothermal case, which is
  what turns the blanket's depth from a 1% claim into a machine-precision one; the
  spurious split was 5 K and 1% of the interface height.

  *Delivering the discharge by volume under-doses the layer the fire is in.*
  "Perfectly mixed" has to mean one *concentration*, and a concentration is a mole
  fraction; splitting by volume gives one partial density, which at one pressure is
  a mole fraction proportional to `T`. Measured on a 3 MW fire: the volume split
  leaves the cool lower layer at **26.8%** while the hot upper one is at 41.6%, and
  the fire — which breathes the cool layer — goes on burning through a flood that
  has already reached its design concentration. Splitting by moles is uniform
  exactly, at any temperature ratio, and identical whenever the layers are at one
  temperature.

  **On the ferry**, a 1215 kg bank into the starboard engine room with the
  watertight door standing open: 18% of the charge leaves the ship, the fired space
  reaches 31.4% and the *other* engine room reaches 24.3% — so flooding one space
  half-floods its neighbour and **neither reaches the design concentration**. And
  the finding that ties the two halves together: stratification makes a low leak
  *worse*, not better. The same discharge into a space with a door at the deck
  holds 39.3% when perfectly mixed and **17.6%** when it separates, because the
  blanket forms exactly where the hole is and drains out of it — and the fire that
  the mixed case extinguishes goes on to release 5.6× as much energy.

  **Mutation-tested, and the survivors were the finding.** 50 substitutions plus
  four deliberate controls, one at a time, each rebuilt and run against the whole
  suite. The per-mutant bound is the substep controller's own arithmetic floor
  rather than a wall clock — `testTheSubstepControllerStaysNearItsFloorWhileFlooding`
  asserts one substep per `maxSubstep` of model time, so a collapsed controller
  fails an assertion in seconds instead of running for hours — with a generous
  wall-clock backstop behind it, scaled off the box's own measured clean run.
  **First pass 40/50; after closing the ten survivors, MUTAFTER.** One kill was a
  **hang** and not a failure, which is this codebase's characteristic kill and the
  reason the bound is not optional. All four controls behaved on both passes.

  Three of the ten survivors were the same shape: **a test that only ever asked
  the question of a pure-air layer, or asked it over a step too short to see the
  answer.** The volume split is a ratio and only its numerator carries the upper
  layer's own excess, so every fixture with air on top was blind to dropping it.
  The boundary loss and the spray cooling are both `C(1 − e^{−r dt/C})`, which
  tends to `r dt` **whatever `C` is** — so every fixture that stepped short enough
  was blind to which heat capacity they used, which is the *same* blindness this
  repo already records for a `c_p` where a `c_v` belonged. And the
  discharge-by-volume error above survived because the entrainment weighting
  masked it: the fire sees a mixture of both layers, and the mixture was near
  enough right.

  **The harness needed correcting too.** Five suites at once on one GPU makes the
  render tests fail intermittently, and a mutation harness scores that as a kill.
  A false kill inflates the rate and hides a real gap, which is the one direction
  a kill rate cannot afford, so Vulkan is compiled out of the mutant builds — none
  of these substitutions is anywhere near it.

  Not here, and named rather than hidden: oxygen is **diluted and never consumed**,
  so a sealed compartment fire that would vitiate itself out does not; a space
  holding a hot smoke layer, cool air *and* an agent blanket has three strata and
  two zones to put them in, so the blanket and the cool air share one; the agent's
  own radiative absorption is not modelled; a re-condensing agent (the adiabatic
  case above ends below CO₂'s own sublimation point) has no solid phase to go to;
  and **`les::demote` writes a `GasCompartment`'s layers without touching their
  agent mass**, so a promoted-and-demoted compartment loses the agent it was
  holding — the reader side is clamped so the state stays consistent, but the mass
  is gone. High-expansion foam and dry powder are still ahead, and neither is a
  gas.
- ✅ **LES promotion for the local compartment** — `engine/sim/les.{hpp,cpp}`,
  `promotion::GasPromoter`, `promotion::gasCandidates`, `les::promote` and
  `les::demote`. One burning compartment's two well-mixed layers replaced by a
  coarse resolved flow on a regular staggered grid, and handed back as two layers
  when it is done; every other compartment stays two-zone throughout.

  **The formulation is low-Mach finite volume, and the case for it is that it is
  the two-zone model's own closure resolved rather than a second model beside it.**
  At one thermodynamic pressure the internal energy density is uniform, so the
  total is `p V / (γ−1)` — which is exactly the identity `fire.hpp` calls "the
  single load-bearing algebraic fact in this file", `V_u/V = U_u/U`. The energy
  therefore crosses the fidelity boundary as **one scalar with nothing to
  interpolate**. Mass is the field and temperature is derived (`ρ = m/V_c`,
  `T = U_c/(m c_v)`), so the equation of state holds identically at every cell at
  every step instead of drifting, which is what the obvious
  variable-density-incompressible arrangement does. Heat enters the *velocity*
  field through the dilatation `div u = (γ−1)(q‴ − Q/V)/(γp)`, whose integral over
  a sealed box is `dp/dt = (γ−1)Q/V` — `fire.hpp`'s own sealed closed form,
  recovered rather than imposed, and reproduced to 1e-9 against the design fire's
  analytic release. The subgrid closure is constant-coefficient Smagorinsky: the
  cheapest thing that is a closure rather than an absence of one, with one
  constant, and switching it off is measured — on a 9 × 7 × 5 m box at 1 m cells,
  150 kW over a 1.4 m pool for 30 s, peak speed 2.63 m/s at `Cs = 0` against
  2.47 m/s at 0.20, so it damps and it is not decoration. A fully compressible
  solve was rejected on the step it would demand: `dx/c` is 1.5 ms on a 0.5 m cell
  for a flow at Mach 0.003.

  **Conservation across the boundary is the acceptance test and it comes out
  exact.** Promote a 10 × 8 × 5 m compartment at 1 m cells carrying a 1.7 m layer at
  500 K, and demote it again: the interface returns **bit-exact** (8.9e-16 m at the
  depths where it is not), the layer masses to **6e-16 relative**, the energy split
  to 2e-16, and the products to the layer they were in. Through 60 s of a 200 kW
  pool fire in the same box the mass residual is
  **1.2e-16** of scale and the energy residual **1.9e-15**, asserted at every step
  rather than at the end — a model can close its books over a run while breaking
  them in the middle.

  **Two corrections found by measuring rather than by thinking.**

  *A demotion cannot integrate the cells it promoted.* The obvious reduction sums
  the cell masses above the interface, and it comes back **9% wrong on a field it
  has just promoted** — because promotion has to smear the interface across the
  cell row it falls inside, and no weighted sum of a mixed row pulls the two
  layers apart again. The equivalent two-layer profile *can*: `ρ_u(H−z_i) +
  ρ_l(z_i−z_f) = M/A` is one linear equation in `z_i`, closed form, and it reads
  only the top row, the bottom row and the total — the three things the smear
  leaves untouched. On that 1.7 m layer at 1 m cells, 9e-2 became 6e-16. The price
  is named and measured: a layer **thinner than one cell** is not representable,
  and the reduction faithfully reports the layer the grid does hold — a 0.2 m layer
  under a 1 m cell comes back 1 m deep at 315 K instead of 500 K, and between one
  cell and one cell short of the floor the round trip is exact.

  *The arithmetic was the limit before the physics was.* A naive sum over the 400
  cells of that box carries `√n ε`, and the two-layer reduction divides
  `M/A − ρ_u H` — a difference of numbers three times its own size — by a density
  gap. That noise was **larger than the round-trip error it was measuring**:
  8.6e-14 m on the interface. Compensated summation took it to the last bit. A
  test cannot be tighter than its own arithmetic.

  **Both directions of the two-zone comparison, because a model that always agrees
  has bought nothing.** Where the zone model is right — a well-mixed compartment
  cooling to its boundary — the resolved model reproduces the lumped exponential
  `exp(−t hA/(Mc_v))` to **1.59% of the 111.9 K the box loses**, worst case, and
  to 0.002 K after five time constants; the two-zone model reproduces the same
  closed form to 1e-4, so all three agree. Where it is wrong, 150 kW over a 1.4 m
  pool for 60 s in a 9 × 7 × 5 m box at 1 m cells puts **45.0 K across the
  deckhead** the zone model says is one temperature, raises a plume at
  **+2.21 m/s** over the fire with a return flow of −0.36 m/s in the corner, and a
  ceiling jet running outboard at ±0.314 m/s on both sides. The control that makes
  that a statement about the fire being local rather than about the model being
  noisy: **the same 150 kW for the same 60 s spread over the whole floor** gives
  8.5 K and no plume (+0.08 m/s), and both fires leave the room at the same mean
  temperature to 0.7 K.

  **The criterion is Alpert's own ceiling jet, and the heat release cancels out of
  it.** The two branches of the correlation either side of `r/H = 0.18` have the
  ratio `(16.9/5.38)(r/H)^(2/3)` — no `Q` in it at all — so the geometric trigger
  is a property of the compartment and a bigger fire cannot switch it on, which is
  what makes the second trigger (the layer standing above the cool gas) do
  independent work. The two branches meet at the crossover to **0.1427%**,
  derived from the published coefficients rather than remembered, exactly as
  Heskestad's two entrainment branches meet to 1.78%. The threshold is anchored on
  the conventional `L/H = 3` limit of a zone model, which for a square room is a
  spread of 5.16; at `spreadPromote = 5.0` a 6 × 5 × 4 m tank (3.09) and a
  12 × 10 × 5 m machinery space (4.23) stay two-zone while a 20 × 8 × 5 m hold
  (5.24) and the ferry's 100 × 19 m vehicle deck (14.75) do not. Hysteresis and
  dwell are `Promoter`'s, tested against the same negative controls: without the
  band a fire hovering at the threshold is built and thrown away six times in
  twelve reviews, with it once.

  **Cost is *not* linear in the resolution, and that is the difference from the
  structural tier.** A Tier-2 zone is exactly `1.7 × elements` core-seconds per
  simulated second because its step is fixed. Halving a cell here multiplies the
  cells by eight *and* halves the advective step *and* lengthens the pressure
  solve, so the measured cost **per cell per simulated second rises** rather than
  staying flat: on one core, over 20 simulated seconds of a 300 kW fire in a
  12 × 10 × 5 m space, 1.6e-6 at 168 cells (1.5 m) against 1.2e-5 at 1456 (0.75 m).
  `cellBudget` therefore bounds the arithmetic honestly and the wall clock only
  loosely.

  **Four of the tests interrogate single cells rather than totals**, because the
  useful defect here is one that cancels when asked globally. The box's wall area
  is right even when every cell's share of it is wrong, so over one 0.1 ms step in
  a 9 × 7 × 5 m box at 1 m cells the corner, edge, face and interior cells are
  asked separately and their losses stand at **3 : 2 : 1 : 0** to 5e-6. One heated
  cell warms over one 1 ms step by `T(γ−1)q δt/(γpV_c)(1 + (γ−1)/N)` to 1.6e-4,
  whose large-box limit is `q δt/(m c_p)` — heating at constant pressure, which is
  what a cell in a low-Mach model is under even though the box as a whole is at
  constant volume.

  **And a fourth that only mutation testing asked for.** Every one of those is
  taken at a step so short that the boundary exponential is linear, and a linear
  boundary term does not care what heat capacity is in it: `C(1 − e^(−r δt/C))`
  tends to `r δt` whatever `C` is. Swapping `c_p` for `c_v` in it therefore
  survived the whole suite. The time constant is only visible at a step comparable
  with it, so the test now takes **one 1 s step at h = 1000 W/(m²·K)**, where a
  1 m cell's own `m c_p/(hA)` is about a second: 15 616 kJ leaves the box against
  the constant-pressure closed form's 15 616 kJ, where constant volume would have
  removed 12 616 kJ — a 19% gap, which is the vacuity guard as much as the
  measurement.

  And a centred fire in a 10.5 × 9 × 4.5 m box — 11 × 9 × 5 cells, odd on both
  plan axes and no two spacings equal — leaves a field that is
  **bit-identical to its own mirror** in both plan axes, with the velocities
  antisymmetric to the bit: red-black relaxation reads only the other colour
  within a half sweep, the Laplacian and the flux gather are accumulated axis pair
  by axis pair, and floating-point addition is commutative — so a reflection maps
  every sum onto itself term for term. Anything less than exact equality there is
  an index reading the wrong axis.

  The quiescent control is exact for the same reason: 120 s of a compartment
  uniformly at **350 K** — not ambient, so a buoyancy referred to a fixed ambient
  density could not hide — leaves every cell mass bit-identical and every face
  velocity at exactly 0.0 over 480 substeps, and the same 9 × 7 × 5 m field at 1 m
  cells against a 288.15 K boundary at 30 W/(m²·K) stratifies by 1.85 K and moves
  at 0.25 m/s.

  Not here, and named rather than hidden: **vents** — a resolved compartment is
  sealed for the duration of its promotion, because handing `fire.hpp`'s vent
  integral a field means deciding what a doorway does with a ceiling jet arriving
  at an angle, which is a separate question; **kinetic energy**, which is 2.4e-6
  of the internal energy at 1 m/s and is not booked; **combustion**, still a design
  fire, though the *entrainment* is now computed rather than correlated;
  radiation, soot optics and a second species; a **heeled** compartment, on the
  same terms as the two-zone interface; and the **plan shape**, which
  `fire::GasCompartment` does not carry at all — the grid is the rectangle its area
  and perimeter imply, falling back to a square of the same area when no rectangle
  has both, and saying so.
- ✅ Volumetric fire and smoke rendering — `engine/gpu/smoke.{hpp,cpp}`,
  `engine/gpu/smoke_gpu.cpp`, `engine/gpu/shaders/smoke.{vert,frag}`,
  `tools/smoke_view`. Two homogeneous emitting, absorbing slabs per compartment on
  the prism the zone model solves on, composited over the lit solid with
  premultiplied alpha. Full detail and every figure in
  `03-renderer-audio.md`'s "Fire and smoke — what exists".

  **There is no raymarch, and that is the physics rather than an optimisation.**
  The plan for this item said "volumetric raymarching against the gas solver's
  density and temperature fields". A two-zone model has no fields: it has two
  masses, two internal energies, two species loadings and an interface height per
  compartment. Each layer is exactly uniform, so the transfer integral is closed
  form — `B(1 − e^(−kd)) + e^(−kd) L_bg` per layer — and the path lengths come
  from an analytic segment-against-half-space intersection. Extinction and
  emissivity are the same exponential, so a clear compartment composites to the
  background *bit for bit*, which is asserted against `HullRenderer`'s own frame.

  **What was refused.** No plume, because the model carries an entrainment *rate*
  and not a shape; no flame, because it carries a mean height and no flame
  temperature; no ceiling jet, no horizontal structure and no gradient inside a
  layer, because a zone is well mixed by definition; and the interface is drawn
  sharp because the model says it is a plane. Each of those would have been
  structure the simulation does not have.

  **Three results that were not the expected ones.** *This fire has no fire in
  it*: a grey layer first puts one byte of red on the screen at 834 K and the
  ferry's 4 MW machinery fire peaks at 531 K, so what a two-zone fire looks like at
  that power is smoke — the glow is a result, and `smoke_view` asserts its absence
  rather than tuning it into view. *The layer goes optically black in about a
  minute* — optical depth across the engine room passes 10 before t = 100 s and
  reaches 511 — so the only visible information a two-zone model has is the
  interface height. And *the layer does not descend monotonically*: it reaches
  2.96 m at t ≈ 300 s and recovers to 3.13 m as the room reaches its vented steady
  state, which the first version of the test asserted away.

  Closed forms on pixels, all agreeing to **one least-significant bit**:
  Beer–Lambert on a slab of known thickness, its square when the thickness
  doubles, the medium stopping at whatever is solid, every pixel against an
  independent statement of the same integral in double, and the interface
  bracketed by two `sim::clipToPixel` projections 2 px apart. **0.007 ms** for the
  volumetric pass at 960 × 540 on the GTX 1070 Ti, against 0.056 ms for the lit
  solid under it.

  The pass first shipped culling the wrong face, and **every closed-form check
  passed at 1 LSB anyway** — from outside a volume both senses give one fragment
  per pixel and the same colour. Only a camera *inside* the medium can tell them
  apart, and that check was off by 152.
- ✅ **Milestone:** an engine room fire that heats a bulkhead until it fails under
  the head of water behind it, and the flooding spreads. Three subsystems, none
  of which know about each other, producing one consequence.
  `fire::wallExchange`, `fire::filmCoefficient`, `thermal::HeatedMember`,
  `thermal::beamColumnMagnifier`, `thermal::twoStripStress`,
  `thermal::temperatureForElongation`, `thermal::Solver::setFilm`, and
  `tools/bulkhead_probe`, which runs the whole chain on the reference ferry.

  **Every piece existed and nothing joined them.** Before this, no file outside
  `thermal.{hpp,cpp}` and its own test named `thermal::` at all. The four links,
  and what each needed:

  1. *Fire → steel.* **The fire already had the loss term; what it did not have
     was a wall.** `GasCompartment::wallConductance` and `wallTemperature` have
     relaxed each layer over its own wetted area since the suppression work — but
     the conductance was MQH's lumped `h_k`, standing in for the whole path, and
     the temperature was pinned at ambient for as long as the fire ran.
     `wallExchange` replaces both: the conductance becomes the **gas-side film**,
     because the wall's own resistance is now solved rather than lumped, and the
     temperature becomes the steel's own surface. The same `expm1` relaxation
     that was standing in for a boundary now *is* one.
  2. *Steel → load.* `thermal::HeatedMember`: restrained expansion as an axial
     compression, the head as a lateral moment integrated off the ship's own free
     surface, joined by the exact beam-column magnifier.
  3. *Failure → hole.* `breachesFromFailedPanels`, unchanged.
  4. *Hole → flooding.* `Ship::step`, with nothing added.

  **Radiation is in the coupling and it is not an approximation.** `thermal.hpp`
  leaves radiation out because it is nonlinear where conduction is not. It does
  not have to be: `σ(T_g⁴ − T_s⁴) = σ(T_g²+T_s²)(T_g+T_s)(T_g−T_s)` is an
  *identity*, so a film coefficient of `εσ(T_g²+T_s²)(T_g+T_s)` delivers the
  Stefan–Boltzmann flux exactly — 112.5 W/(m²·K) for 800 °C gas against 227 °C
  steel, of which 87.5 is radiation, and the two forms agree to **1.0e-13** over
  280–1400 K. Where they *disagree* the factored form is the right one: near
  equilibrium `T_g⁴ − T_s⁴` is a cancellation and the factored form has none, so
  a nanokelvin apart the difference form has lost six digits. At `T_g == T_s` the
  flux is exactly zero, which is what the cold control rests on.

  **The result, on the ferry's own x = −8 bulkhead** — 9.5 mm plating on 180 × 10
  flat bars at 700 mm, spanning the tank top to the bulkhead deck, 4 MW machinery
  fire in the port engine room, both aft holds 45% flooded:

  | | fire alone | head alone | both |
  |---|---|---|---|
  | member peak | 252.2 °C | 20 °C | 174.8 °C |
  | worst utilisation | 0.934 | 0.274 | **1.001** |
  | failed | no | no | **at 1935 s** |

  It fails at a member temperature of **151.6 °C**, where `k_y` is exactly 1 and
  the steel has lost none of its strength. At that state, each cause alone:
  restrained expansion **0.513** of the member's capacity, the head **0.234**
  unmagnified and **0.488** magnified — the axial load is what magnifies it, by
  **2.086**. **A purely additive check reads 0.747 and has not failed**, which is
  the whole of what the coupling buys and is what two subsystems that do not know
  about each other would have produced. One panel goes, 0.490 m², and 307.7 m³
  reaches the machinery spaces through it; four compartments end wet and she
  stays afloat at 6.2° of heel.

  **The milestone's sentence reproduces, and the honest answer is that it does so
  over a band.** The restraint on a heated member is set by the stiffness of the
  structure at its ends — a bulkhead deck and a tank top — and that is in *none*
  of the three subsystems, so it is an input the model cannot derive. The tool
  measures the band rather than asserting a value: the sentence is true for
  **0.2194 ≤ r < 0.2505**, a factor of 1.142. Below it nothing fells the bulkhead;
  at or above it the fire alone does and the head is not what fails it.

  **Why the band is that narrow is itself the finding: the water that loads the
  bulkhead also cools it.** A dry hold behind leaves the member **77.4 K hotter**
  (252.2 °C against 174.8 °C), so the two controls are not one change apart. The
  same-state decomposition above is the sharper statement for exactly that reason.

  **Two controls had to be repaired before either was a control at all.**
  *The fire has to end.* Steel asymptotes to the gas temperature it stands in and
  the restrained-buckling limit is a fixed temperature, so a fire that burns
  indefinitely fells the bulkhead on its own — measured at 2760 s — and the
  sentence would have been true only over a window in *time* the tool itself had
  chosen. A design fire that grows, burns and decays gives the steel a peak.
  *And the ship has a second leak path.* The ferry carries an unsealed 0.04 m²
  cable transit through this very bulkhead; left open it floods the machinery
  spaces, downfloods the vehicle deck and lolls her, which fells the bulkhead in
  the **water-only control at 3645 s with nothing burning anywhere**. Real, and a
  different experiment, so the milestone run seals it and `--cable-transit`
  reopens it.

  **What the plating does, and it is not a hole.** Under the horizontal
  compression the hot band of the bulkhead takes from the cold band below it —
  `twoStripStress` at the hot fraction the fire's own interface height gives —
  the 0.70 × 0.70 m panels between the stiffeners buckle at t = 1360 s. A buckled
  panel goes out of plane, sheds its in-plane load and stays watertight. It is
  reported and it is deliberately not fed to `breach`.

  **Three defects the chain found in its own first version**, each of which
  produced a plausible number: filling a hold *after* `Ship::initialise` leaves
  the air mass from when it was empty, so the trapped gas reads 184 kPa — 8.3 m
  of head that is not water — and the bulkhead failed in the first step of every
  run, controls included. One film per gas layer misbooks **8%** of the exchange,
  because the radiative coefficient goes as `T_s²` and the foot of a bulkhead is
  held near ambient by the water behind while its head is at 500 K; banding the
  films by element row takes that to **0.102%** of a 261.7 kW peak exchange on
  the ferry, and to 0.000% on the unit fixture — and it needs no re-preparing,
  which is what lets the steel's energy account span the run and close to
  **7e-13** of the 0.1534 GJ it moved. And a
  restraint window measured off a run that *applied* its own damage comes back
  equal to whatever restraint that run was given, because the failure floods the
  compartment and relieves the head — so the bound is measured on a fourth pass
  that evaluates the failure and does not open it.

  Not here, and named rather than hidden: no redistribution, so a member that has
  gone hands nothing to its neighbours; first yield of the extreme fibre rather
  than a plastic mechanism, so the check is conservative by the shape factor;
  pin-ended, so real end fixity fails later; and the axial restraint is uniform
  along the member where the temperature is not — the member is reduced to one
  *equivalent uniform temperature*, the inverse of the elongation curve at the
  mean of its own elongations, which is the right integral but still one number.

## Phase 5 — Fluids — *~9 em*

- ✅ Sparse FLIP/APIC solver for interior water — `engine/sim/flip.{hpp,cpp}`.
  Particles carrying mass and an affine velocity, a staggered pressure grid stored
  as 4×4×4 tiles in a hash map, a Jacobi-preconditioned conjugate-gradient
  projection with a voxelised free surface, and RK2 advection through the
  extrapolated grid field. `tools/flip_probe` runs the studies; `tests/test_flip.cpp`
  owns the closed forms.

  **It is APIC, and the measurement rather than the citation says why.** All three
  transfers are reachable from one code path (`Params::affine`, `Params::flipBlend`)
  so the choice can be re-derived. On a 0.8 m cube in solid-body rotation, PIC
  keeps **89.2%** of its angular momentum through ten particle→grid→particle
  transfers and APIC **98.2%** — 1.08% against 0.18% per transfer, **6.1× less**.
  In a sloshing tank the ordering reverses: FLIP carries **87.6%** particle-borne
  velocity the grid never sees, APIC **1.3%**, and PIC damps the very mode it is
  meant to measure, crossing the centreline three times in five seconds where APIC
  crosses eight. On the period itself APIC comes out **+4.45%** against PIC's
  +9.16% and FLIP's +9.81%. So APIC buys FLIP's conservation at a sixty-seventh of
  FLIP's noise, and is twice as accurate as either on the one closed form
  available.
  The kernel is the **quadratic B-spline** because only there is APIC's `D_p` a
  constant — `h²/4 I`, asserted directly to 3.2e-16 relative — where trilinear
  weights make it *singular* whenever a particle lands on a node.

  **What is asserted, and at what.** A hydrostatic column does not move at all:
  after 25 steps **not one particle coordinate has changed**, because the residual
  velocity is 6.2e-16 m/s and `x + v dt` at that size is a no-op — against 12.3 mm
  and 0.49 m/s of unprojected free fall. Its pressure is `ρ g h (K − k)` **cell by
  cell** to 2.6e-15 relative. The projection reproduces a discrete Helmholtz
  decomposition — a MAC curl plus a discrete gradient half again its size — giving
  the divergence-free part back to **8.8e-16** relative per face and the potential
  back as the pressure to **9.4e-16** per cell, on a 7×5×6 grid so an index that
  reads the wrong axis cannot survive. A block in free fall follows its
  integrator's own closed form `−g dt² N(N+1)/2` to 1.1e-15, with **exactly zero**
  sideways. Mass is not conserved to a tolerance, it is *exact*: `expectNear(...,
  0.0, 0.0)` on the residual, with the particle count asserted as an integer
  beside it, through a dam break that clamps 18 610 particles against walls.

  **Sparsity is tested for what it is for.** An empty 400³-cell room — 64 million
  cells, 7.4 GB dense at the 115 bytes a cell this structure costs — allocates
  **zero tiles and zero bytes**. The same water in a 20³ room and a 400³ room
  allocates the same 48 tiles and produces **bit-identical** particle positions
  after fifty steps, which is the statement that the room's extent does not enter
  the arithmetic at all. Arrival is the case sparse solvers fail at, so a single
  particle is walked 6 m across six tile seams in 1 mm steps with the invariant
  checked at every position: the face weights of each component sum to the particle
  mass, to 1.8e-15 kg of 7. A halo one tile too small drops part of that sum and
  nothing else notices.

  **And the tile count is asserted against one derived by hand**, not against
  another run: a 6×6×4-cell block dilated by the 3-cell halo is tiles 0…3, 0…3,
  1…4, which is 64. That check is here because `Solver::tiles()` shipped returning
  `tileKey_.size()` — three coordinates per tile — and reported **192**. Every
  other tile assertion in the suite compares one count with another, equal in two
  rooms or changed over a fall, and a factor of three is invisible to all of them.
  Same shape as the half-bandwidth `reduction.cpp` reported and the suite believed.

  **Sloshing, against `2π/√(g k tanh k d)`.** The first mode comes out **+5.66%**
  on six cells of depth and **+3.97%** on nine, converging on about +3%, which is
  the nonlinear correction at `a/d = 1/3` and does not go away with the grid. The
  amplitude has to be a *couple of cells*: a voxelised free surface cannot
  represent a sub-cell tilt and therefore has no restoring force to give, and the
  first version of that test seeded a sixth of a cell, did not slosh at all, and
  reported a period of zero.

  **Shared with `les.cpp`, and deliberately not.** Shared: the staggered grid, the
  all-Neumann compatibility discipline (a fully submerged compartment determines
  its pressure only up to a constant, so the mean is removed — the same problem
  `les.cpp` meets from the heating side), compensated summation on every conserved
  total, and a substep controller that *publishes* when its budget binds. Not
  shared: the linear solver. `les.cpp` uses red-black SOR because its answer must
  commute with a reflection; this needs a residual at machine precision instead,
  and on the same 24-cell column SOR at 1.7 takes **1 170 sweeps** to reach 1e-13
  and does not reach 1e-15 in 200 000, where Jacobi-CG reaches **exact zero in 24
  iterations**. Also not shared: the transported variable — mass is on the
  particles here, which is why its conservation is an identity rather than a
  residual.

  **Not here, named rather than hidden:** solids are one axis-aligned box, so no
  stair, girder or heeled deck; the free surface is voxelised, so first order and
  no ghost-fluid pressure; no volume correction, so FLIP's usual particle
  clumping is unopposed; no viscosity, surface tension or two-phase air; no moving
  solids; single threaded. And it is **not wired into `sim::Ship`** — that is the
  next item and doing it here would have made this untestable.
- ✅ Quiescent ↔ dynamic escalation with exact mass conservation —
  `engine/sim/water_promotion.{hpp,cpp}` following the `StructuralPromoter` /
  `GasPromoter` pattern. `WaterPromoter` is a motion-gated state machine: a
  compartment qualifies when roll rate OR lateral acceleration exceed thresholds,
  builds a dwell streak over consecutive reviews, promotes when the streak reaches
  the dwell count and particle/tile budgets allow, and demotes when motion drops
  below hold thresholds for `hold` reviews. `promoteWater()` creates a
  `flip::Field` seeded through the wetted depth (not full bbox height —
  forepeak at 5 m³ seeds 745k particles, not 77M), and `demoteWater()` reads mass,
  centroid and quiescent level back. Round-trip mass is **exact** (0.0 tolerance,
  not epsilon), tested in Section 4 of `tests/test_water_promotion.cpp`.
  
  **Seeding only wetted depth matters.** The first version seeded `bboxLo.z` to
  `bboxHi.z` (full compartment height) and generated 77 million particles for 5 m³
  of water in a 232 m² forepeak — 100× too many. Corrected to seed through
  `bboxLo.z + waterVolume/planArea`, the same compartment gets 745k particles at
  ~1000/m³, which is the target density. A measurement probe (`/tmp/wp_probe.cpp`)
  found it; the fix is one line changing `hi[2]`.
  
  **Wired into Ship, not yet active.** `Ship::waterPromoter_` and
  `Ship::activeWaterFields_` exist with proper copy/move semantics (copy leaves
  fields empty since girder calculations copy ships without active FLIP; move
  transfers ownership and clears source). `Ship::step()` carries TODO markers for
  periodic `review()` calls, promotion/demotion handling, and FLIP field stepping.
  The infrastructure is in place; activation is the next step.

  **A capsizing ferry is a slow event, and that is the measurement that says
  activation is safe.** `tools/water_probe` runs the criterion alongside the real
  `doors` scenario without promoting anything, and the ship that lolls to 58.8°
  with 5591 t aboard peaks at **0.0042 rad/s of roll rate — 12× below the 0.05
  threshold** — and promotes nothing in 892 reviews. Progressive flooding is
  minutes of quasi-static heeling, not the sloshing this tier is for. So turning
  the promoter on inside `Ship::step()` moves no gated scenario figure, which is a
  measurement rather than an argument, and the front page's three outcomes are
  safe from it.

  **A promoter that fires nowhere would be broken rather than safe**, so the probe
  carries its own control: the same ferry in a 2 m regular beam sea at her own
  9.92 s roll period reaches **0.1525 rad/s, 3× over the threshold, and promotes**.
  The pair is the finding — silent through a capsize, live in a seaway.

  **The two budgets were set independently and are not independent, which made
  the whole tier unreachable.** `estimateFlipCost` puts 1000 particles and
  `1/(64h³)` = 125 tiles in every m³ of water, so each budget states a *volume*:
  100 000 particles is 100 m³, and 2000 tiles was **16 m³**. The tile budget
  therefore bound 6.25× sooner and the particle budget could never be reached at
  all.

  **The memory share took three attempts to state, and the sequence is the
  lesson.** This entry first said the header's "particles are the memory
  bottleneck" was wrong and put tiles at 92% of the footprint — a byte claim
  computed at the estimator's 1000 particles/m³, in the same breath as
  establishing that byte claims must use the solver's real 64 000. Corrected to
  15%, still off, because `flip::Particle` was documented as "~80 bytes" in a
  comment that itemises sixteen doubles. At the measured 128 B a cubic metre is
  920 kB of tiles against **8.19 MB** of particles: **tiles are 10%**, and the
  line originally being corrected was right all along.

  The consequence was not a tuning matter. **Any compartment over 16 m³ was
  refused outright**, however hard she rolled, on a ship whose compartments run to
  1232 m³. The probe seeded the vehicle deck — 462 m³, the largest free surface a
  ro-pax has and the reason this tier exists at all — and watched it refused 154
  times while a 1 m³ trickle in a forward hold was promoted instead. **The
  refusals were silent**, because nothing read `WaterReview::problems`; a
  criterion qualifying 2525 times and promoting once looks exactly like
  hysteresis working until something prints the reason.

  The old ceiling was ~146 MB and the deck it refused would have been 4.2 GB, so
  the number was defensible and the way it was expressed was not. **Those
  megabytes are at the solver's real seeding and not the estimator's**:
  `estimateFlipCost` bills 1000 particles per m³ while `flip::seedBox` at 2³ per
  cell puts 64 000 there, so a memory figure read off the budget arithmetic is
  **8.7×** light — including, at first, the ones in this entry, which said 64×.
  The particle *count* is 64× light; tiles are geometric and do not move with the
  seeding, so the total is not, and the two are not interchangeable. The budget is denominated
  in the estimate and is self-consistent; any sentence about bytes has to use the
  other number. The budgets now agree at 100 m³ — measured to
  1.00× rather than asserted — and a compartment past `maxVolumePerCompartment`
  is rejected with *its own reason* rather than as a transient budget condition,
  because a vehicle deck that cannot be afforded is a finding about this tier and
  not a passing shortage. With that, the beam-sea control promotes **7**
  compartments — both engine rooms, three wing tanks and two holds — where it
  promoted one. A vehicle deck still needs a coarser grid than h = 0.05, which is
  the next question this tier has to answer.

  `testBudgetsAdmitTheSameVolume` asserts the agreement against the arithmetic
  rather than against the constants, so it fails when one number is edited and
  not the other — which is exactly how it broke.

  **And then the cost was measured, which settles activation on its own.**
  `coreSecondsPerCompartment` was 5.0, marked *"estimate, will be measured"*, and
  it is the number the tier's affordability rests on. `water_probe --cost` steps a
  real `flip::Solver` at h = 0.05:

  | volume | particles | tiles | ms/step | core-s/sim-s | × the estimate |
  |---:|---:|---:|---:|---:|---:|
  | 1 m³ | 65 650 | 525 | 279 | **27.9** | 5.6× |
  | 5 m³ | 325 424 | 1 560 | 1 462 | **146.2** | 29× |
  | 20 m³ | 1 276 292 | 4 104 | 5 616 | **561.6** | 112× |
  | 100 m³ | 6 351 696 | 17 019 | 30 304 | **3030.4** | 606× |

  The estimate was not merely low, it was the **wrong shape**: a per-compartment
  constant where the true cost scales with the water. So the seven compartments
  the beam-sea control promotes would cost **~1000 core-seconds per simulated
  second** if each held only 5 m³, against the 1.0 that realtime needs.

  **It cannot be stepped around.** At the production dt = 0.02 the substeps double
  from 2 to 4 and the wall time doubles with them, leaving core-s/sim-s within 4%
  of the dt = 0.01 figures: the CFL condition sets the work, not the timestep. The
  cause is the same 64× that made the memory figures wrong — `flip::seedBox` at 2³
  per cell puts 64 000 particles in a cubic metre at this `h`, and every one is
  scanned per substep.

  So **activation is blocked on cell size, not on wiring**, and that is the third
  independent finding pointing at the same place: the deck needs a coarser grid
  the physics will not allow, the budget admits volumes the clock cannot afford,
  and the estimate that hid both was never run. `Ship::step()` keeps its TODO
  markers deliberately — wiring a tier that costs 1000× realtime would be a
  working feature nobody could switch on, and the honest state is a measured
  ceiling with the integration left until there is something affordable to
  integrate.

  **And the vehicle deck cannot be bought with a coarser grid, which is a
  finding rather than a budget.** The obvious answer to "462 m³ costs too much at
  h = 0.05" is to coarsen until it fits: the deck needs 57 750 tiles and 29.6
  million particles at 2³ seeding per cell — **4.2 GB**, not the 484 MB a
  1000-particles/m³ estimate suggests — and at h = 0.20 the same water is 902
  tiles and 66 MB, which is affordable outright. But the deck floods *shallow
  and wide*: 462 m³ over its own 1868.4 m² is **0.247 m deep**, so h = 0.20 puts
  **1.2 cells across the entire depth** where `flip_probe`'s own sloshing study
  requires about two cells of *amplitude* before a voxelised free surface has any
  restoring force at all. Coarsening to afford the memory destroys the physics
  the tier exists for, and it would do it quietly — the run would complete, the
  mass would still be exact, and the deck would simply not slosh.

  So the vehicle deck is not a promotion candidate at any cell size this solver
  currently offers, and `maxVolumePerCompartment` says so as a limit rather than
  leaving it to be discovered as a silent refusal. What it needs is a different
  model — a shallow-water or depth-averaged treatment, where the small dimension
  is integrated out rather than resolved — and that is a Phase 5 item in its own
  right, not a parameter of this one.

  **This paragraph used to end by saying the tier "can serve the deep narrow ones
  — engine rooms, holds and wing tanks", and the measurement says otherwise.** At
  a fixed 20 m³, `water_probe --cost` gives **665** core-s/sim-s for a 1 × 1 × 20 m
  tank, **586** for the 2:1 hold and **535** for a 0.25 m-deep deck: the deep
  narrow shape is the *most* expensive of the three, not the least. Cost follows
  the water's own extent, and a tall column touches more tiles per cubic metre
  than a flat sheet does. So the consolation prize this entry awarded itself was
  written without measuring — inside an entry whose entire subject is a cost that
  had never been measured. **There is no compartment shape at this cell size that
  this tier can afford**, and the honest statement of scope is that it serves
  nothing on this ship today.

  **The control found two frame errors that every unit test had agreed with.**
  `computeLateralAccel` returned `length(state.velocity)` — a *speed* in m/s —
  compared against a threshold labelled m/s². The beam-sea run read 2.26 "m/s²"
  that was 2.26 m/s of orbital velocity; differencing the velocity properly gives
  0.72 m/s², so the old figure was 3× too large. The sign of that error is the
  dangerous one: a ship drifting steadily at 3 m/s has no lateral acceleration at
  all and would have promoted every compartment aboard, while a ship snapped
  sideways from rest by a collision has a large one and a small speed. And
  `computeRollRate` took `angularVelocity.x` — a **world** component — as roll,
  which is true only while she is upright, at exactly the attitudes this criterion
  is not for. Both are now read in the body frame off `state.orientation`.

  This is the third time this repository has confused a world vector for a body
  one (`state.velocity.x` read as speed, the roll stiffness about the wrong point),
  and the first two are already in CLAUDE.md's table. **The tests could not have
  caught it: the harness set `state.velocity = {0, a, 0}` and called it an
  acceleration**, so the fixture and the code shared one misconception and agreed
  with each other. `setLateralAccel` now produces a real difference over a real dt,
  and two new tests fail against the old implementation — verified by reverting it,
  which takes the suite to 7 failures.
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
   scripted-threshold tearing) rather than killing the project — and Tier 1 now
   exists and is measured at 1200× cheaper than Tier 2 on the same plating, so
   that fallback is a real one rather than a plan.
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
