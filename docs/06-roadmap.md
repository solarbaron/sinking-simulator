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
- 195 086 closed-form validation checks
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
  is 2800× rather than 10⁻⁵; "the static interface response improves with mode
  count" is false, it starts exact and the modes buy the *interior* and the
  dynamics; and the standard "cut off at twice the band of interest" buys 0.6%
  inside the 10 Hz hull-girder band, not four figures, because the cutoff is a
  frequency of the fixed-interface spectrum and the band is a frequency of the
  assembled one.

  Cost on the same plating: **0.10 core-seconds per simulated second** for Tier 0
  over the whole ship, **0.35–0.41** for this patch at Tier 1, **1155** for the
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
- Suppression by gas: CO₂ and inert-gas total flooding (oxygen displacement, and
  the sealing requirement that makes it fail with a door open), high-expansion
  foam, dry powder. All of them need a second gas species with its own `R` and
  `γ`, which is the same change the steam mass above wants.
- LES promotion for the local compartment
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
   scripted-threshold tearing) rather than killing the project — and Tier 1 now
   exists and is measured at 2800× cheaper than Tier 2 on the same plating, so
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
