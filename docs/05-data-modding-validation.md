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

Declarative, text, version-controlled, human-diffable.

### What exists now

`engine/sim/shipfile.{hpp,cpp}` reads a **single text file** that builds a whole
`sim::Ship`, and `ships/ferry.ship` is the prototype ferry expressed in it.
`game/prototype/ferry.cpp` stays in the tree as the reference the file is checked
against, not as the way ships are authored.

The grammar is whitespace-separated tokens, one directive per line, `#` to end of
line for comments — deliberately the same lexical idiom as
`engine/gpu/materials/*.materials`, because this repo already has one text format
and a second style would be a second thing to learn. Blocks are `ship`, `hull`,
`compartment`, `opening` and `pump`; the header of `shipfile.hpp` is the grammar
reference. Abridged:

```
ship_format 1

ship ferry_120m
    deck_edge_z 7.0
    lightship_draft 5.5            # or lightship_mass, in kg
    lightship_cog -1.5 0.0 7.80
    gyradii        7.0 28.0 29.0

hull
    waterlines  0.000  1.000  1.800 ...     # ascending, m
    station -60.0  1.350000  2.463693 ...   # x, then one half-breadth per waterline

compartment engine_room_s
    box  -8 -8 1.8   20 0 7.0      # six bulkhead and deck planes
    permeability 0.85
    vented false

opening breach_er_s
    between sea engine_room_s
    at 6 -9.0 3.0
    area 2.4
    discharge 0.62
    kind breach
    open true

pump bilge_er_s
    drains engine_room_s
    capacity 0.060
    max_head 25.0
    on false
```

SI throughout. No third-party parser: this repo has a hand-written PNG codec
rather than take libpng, and a ship file does not justify a JSON or TOML
dependency either.

**It is the same ship, and that is measured rather than asserted.**
`tests/test_shipfile.cpp` loads
`ships/ferry.ship` and holds it against `game::buildFerry()` — hull triangle by
hull triangle, every compartment volume, every opening field, then both floated
and their displacement, draft, GM, KB, waterplane area, freeboard and the whole
0–60° GZ curve compared. The file's six-decimal offsets table is the only
difference between the two ships and it costs **3.3e-9 relative on displacement,
2.7e-8 relative on compartment volume, 4.1e-8 m on GM and the GZ curve**. The
format itself is lossless — seventeen significant figures round-trip bit for bit;
six decimals is a micron, and an authoring choice about what a human can diff.
Both ships are then run through the damage-control scenario in lockstep, and
`scripts/verify.sh full` runs all three 900 s flooding scenarios from the file and
requires the same verdicts as the compiled ferry.

**It fails closed, and says why.** Everything parses into a scratch definition and
is committed only when the whole file has parsed, every name has resolved and
every compartment has come back from the clip with geometry in it. Diagnostics
carry `<origin>:<line>: <reason>`. CLAUDE.md records a `World::load` that failed
open and left a half-built world; the instrument that caught it — every truncation
of a valid file — is pointed at this loader too, at both line and byte
granularity.

**Unknown keys are rejected, not ignored,** which is the opposite of the usual mod
advice and is deliberate. `permeabilty 0.85` on a machinery space is a typo that
leaves the space at a *plausible* default and nothing downstream looks wrong; the
ship just floods more than the author wrote. Forward compatibility is bought
instead with the `ship_format` stamp: a file from a newer build is refused **by
name** rather than silently missing the field that changed the ship's meaning.
Data a mod carries that the simulation does not read belongs in a sibling file —
a ship is a directory, and `art/`, `systems/` and `loading/` are where that goes.

Everything a key can get wrong is caught at load: a value out of range, a
non-finite number, a unit suffix on a number, a number where a name should be, a
name that never resolves, a station list out of order, the wrong count of
half-breadths, a duplicate name, a key set twice, a box that clips to nothing, a
hull that does not close, and a deck edge or a lightship draft outside the hull
the offsets describe — the last two because they show up only as a freeboard or a
displacement nobody has anything to check against.

Every key is required unless the format documents a default, and the only
documented defaults are `sea_density`, `damping`, `added_mass` and
`added_inertia` — engine-wide physics placeholders rather than statements about
this ship. Permeability and the vented flag in particular are *not* optional: the
defaults for both are plausible, and a plausible wrong number is the failure mode
this repo keeps finding.

### What the format cannot yet express

Named so that nobody discovers them by finding a field silently missing:

- **Compartments are authored, not derived.** Each one is six axis-aligned planes
  clipped against the hull. That is `ferry.cpp`'s model as data — it is *not* the
  derived subdivision described below, and it carries the same weakness: nothing
  guarantees the spaces tile the hull. `Ship::validate()` catches overlap; gaps
  are indistinguishable from unmodelled structure. Raked and curved bulkheads,
  L-shaped spaces, and any union or subtraction are simply not sayable.
- **The hull is an offsets table and nothing else.** No STEP/IGES, no DXF, no
  point cloud, and no reference to an external offsets file — the table is inline.
- **One file, not a directory.** No scantlings, no systems graphs, no machinery or
  propulsion (`Ship::propulsion` is reachable only from code), no loading
  conditions, no precomputed hydrodynamics, no art. Radiation is attached by
  calling `Ship::attachRadiation`, not by a key.
- **No initial condition.** Every ship loads intact and upright. Water already in
  a compartment, an initial list, a starting draft other than equilibrium: none of
  them can be stated.
- **No scenarios.** Damage events, sea state and the damage-control schedule are
  still C++ in `game/prototype/main.cpp`.
- **No exporter.** The engine reads the format and does not write it, so the
  editor of §3 has half its round trip. This is the next thing to build, and the
  bit-exact round-trip test exists so that it can be.
- Names are single tokens with no spaces, and a ship has no human-readable title.
- Permeability is one scalar per space; there is no distribution and no
  distinction between a space full of machinery and a space full of cars.

### Where it is going

The reference structure below is still the target:

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

Whether it becomes a directory of TOML files or a directory of files in the
format above is not settled. What the shipped loader has already shown is that a
whole `sim::Ship` fits in one file readably, and that a parser worth trusting is a
few hundred lines — so the argument for a dependency is weaker than it looked.

**Compartments are derived, not authored.** You place bulkheads and decks as
planes or surfaces; the toolchain intersects them with the hull interior and
produces the compartment polyhedra. This is the fix for the prototype's biggest
data problem, it guarantees the compartments tile the hull without gaps or
overlaps, and it means changing the hull form does not silently invalidate every
compartment volume.

`ships/ferry.ship` does **not** do this yet — it states one box per compartment,
which is the prototype's model written down rather than fixed. The half-step it
does take is that a box is clipped to the hull rather than trusted, so the wing
tanks taper into the turn of the bilge and a box that misses the hull is refused
at load instead of becoming a space that never floods. The remaining half is a
shared plane set and a partition of the interior between them.

**Everything derived is cached and content-hashed.** BEM coefficient solves take
hours; windage sweeps take longer; FEM Craig–Bampton reduction takes 0.1–5 s for
a patch of side shell and would take minutes for a whole ship.
These are build artefacts keyed on the hash of their inputs, computed once, and
distributed with the ship.

## 2. Importing real ships

Real hull forms come as:

- **IGES / STEP** — the CAD interchange formats naval architects actually use.
  Import via OpenCASCADE, tessellate to the watertight mesh the sim needs.
- **Offset tables** — the classical station/waterline half-breadth table. Already
  supported by `makeHullFromStations()`.
- **Principal particulars** — **implemented**, `engine/sim/hullform.{hpp,cpp}`.
  See below: the case where no offsets exist at all, which is most ships.
- **DXF / lines plans** — traced body plans, for historical vessels where that is
  all that survives.
- **Point clouds / photogrammetry** — for scanned wrecks and preserved ships.

Public hull forms available for validation and as shipped content: DTMB 5415
(naval combatant), KCS (container ship), KVLCC2 (tanker), S175, JBC, Duisburg Test
Case, Wigley and Series 60 parametric hulls, and the various IMO/HARDER damage
stability test cases.

### From principal particulars alone — **implemented**

Offsets are proprietary for almost every real hull and published for almost none,
while *every* reference gives length, beam, draft and a handful of form
coefficients. `makeHullFromParticulars()` builds a hull that **measures** as the
requested ship.

The construction is the classical one. `Cp = Cb / Cm` separates how full the
midship section is from how full the ship is along its length. The midship
section is a rectangle with a radiused bilge, and the radius follows in closed
form from Cm, since each bilge removes `r²(1 − π/4)`. The sectional area curve is
`f(u) = 1 − (1 − e)|u|ⁿ` with one exponent and one end value per end; the mean of
the exponents sets Cp and their difference sets LCB, both through closed forms,
so the pair is *solved* rather than searched. Every station is then the midship
section scaled in breadth, which makes the sectional areas exact by construction
and leaves only tessellation between the request and the result.

**It is not the real hull, and that distinction is load-bearing.** Two ships with
identical coefficients can have visibly different bodies and measurably different
seakeeping. For stability, flooding and manoeuvring — dominated by volume,
waterplane and their distribution — this is a reasonable stand-in. For
**validation against published RAOs it is not**, because the comparison would
then be against a hull that is not the benchmark ship. That needs real offsets,
and any such claim should say which table it used.

**Measured accuracy**, S-175 form, block-coefficient error against what was asked:

| stations | 11 waterlines | 21 | 41 |
|---|---|---|---|
| 21 | 0.195% | 0.104% | 0.180% |
| 41 | 0.351% | 0.053% | **0.022%** |
| 161 | 0.392% | 0.094% | 0.019% |

The error is dominated by **waterline** count, not station count, because the
waterlines are what resolve the bilge arc — going 41→161 stations moves it from
0.022% to 0.019% where going 11→41 waterlines moves it sixteenfold. The default
therefore spends its triangles on waterlines.

LCB is analytic off the area curve, so refining *waterlines* leaves it alone:
worst 4.0 × 10⁻⁵ of Lpp over that sweep. Refining **stations** is refining the
area curve itself, and it moves — 2.3 × 10⁻⁴ over all nine meshes. This used to
claim 6 × 10⁻⁵ "at every resolution", which was true of the sweep the test ran and
not of the table beside it.

> **Correction, and none of these nine numbers had a producer.** The table read
> 0.52 / 0.23 / 0.15, 0.43 / 0.13 / 0.06, 0.40 / 0.10 / 0.025. Grepping every one
> of them across `tests/`, `tools/` and `engine/` returned nothing:
> `test_hullform.cpp` swept the same waterline counts, computed the middle row as
> an absolute residual, asserted only that refining improves things, and printed
> none of it. It now prints all nine and `check-figures.sh` gates them.
>
> The measured values are two to five times smaller than what was published, and
> the most likely reason is in the paragraph below this one: the stations were
> moved off uniform spacing and concentrated in the tapers. The table appears to
> predate that change. It is offered as the likeliest explanation rather than a
> demonstrated one — nothing recorded which build produced the old numbers, which
> is the whole argument for gating a figure instead of measuring it once.
>
> **One row is not monotone and that is not a rounding artefact.** At 21 stations
> the error goes 0.195% → 0.104% → 0.180% as waterlines refine, so more waterlines
> is briefly *worse*. At 41 and 161 stations it falls monotonically. A coarse
> station set and a fine waterline set disagree about where the bilge is, and the
> volume error changes sign somewhere in between; the default is nowhere near that
> corner, but the table should not be read as "finer is always better".

**Two defects this turned up, both silent.** The first version solved the area
curve *without* the transom and stem end values, then rescaled the exponents to
recover Cp — which leaves LCB wherever the rescale put it. The error tracked
`transomFraction` exactly: 0.0003 of Lpp at a cruiser stern, 0.021 at a wide
transom, which is a metre and a half of LCB on a frigate. And
`bilgeRadiusForMidshipCoefficient()` clamps, because a radius cannot exceed the
draft or the half-beam — so a fine enough Cm on a shallow hull is unreachable and
was being quietly rounded up, meeting Cp and missing Cb, which is the confusing
way round. Both are now reported through the `problems` out-parameter.

**Parallel middle body.** Real ships have a span of constant midship section — a
VLCC is close to half parallel — and the first version of the area curve had
none, so a fine hull came out canoe-like. Nothing about Cb, Cp or LCB can see the
difference, which is why every coefficient test passed; a rendered frame showed
it immediately. The curve is now flat over |u| ≤ p and tapers over the rest, with
both integrals still closed forms.

That change also moved where the stations belong. Uniform spacing spends them
where the section never changes and starves the ends where it changes fastest —
and the taper is *steeper* with a parallel body, because the same fall-off is
squeezed into (1 − p) of the length. Measured on the S-175 form with p = 0.40,
uniform spacing at 41 stations put Cb 0.76% high and needed 161 stations to reach
0.04%. Concentrating them in the tapers, where a constant section needs only two,
brings 41 stations to 0.44%.

`kvlcc2Particulars()` and `s175Particulars()` are supplied, the former matching
the `HullParams` already in `propulsion.cpp` so the two describe one ship.

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

The reader half of this landed with the ship definition format; the writer has
not. An editor needs to emit the format as faithfully as it consumes it, which is
why `tests/test_shipfile.cpp` establishes that full-precision offsets survive the
round trip bit for bit — the exporter has a target to hit before it is written.

## 4. Scripting and mods

- **WASM** for behavioural mods (systems logic, scenario scripts, AI, UI). Chosen
  over Lua for sandboxing, performance, and language choice. Mods cannot corrupt
  the sim state; they interact through a declared API.
- **Data mods** need no code at all — a new ship is a directory. The first half of
  this is real: `./shipsim --ship=path/to/ship` runs any file the format can
  describe, and nothing recompiles. A bad file stops the program with a line
  number rather than falling back to a built-in ship, because a mod that silently
  does not load is worse than one that refuses to.
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
| Free-surface loss of GM | ρ·i/Δ, matched within 15% — the slack is the two barges not having quite the same KG, not the model. `test_fire.cpp` remakes the comparison with the centroid held identical and gets **1.4e-6** on a box barge, and 2e-4 on the ferry's vehicle deck where the second moment is itself a mesh integration |
| Trapped air arresting flooding | Boyle's law, pV conserved within 2% |
| Water mass conservation across a network | flow through the breach = water held |

`tests/test_shipfile.cpp`, on the ship definition format, against the compiled
ferry it has to reproduce:

| Check | Reference |
|---|---|
| Hull mesh from `ships/ferry.ship` | `game::buildFerry()`, vertex and index for index, to 1e-6 m |
| Displacement and waterplane | the same ship, to 1e-8 relative |
| Compartment volumes | the same ship, to 1e-7 relative |
| Draft, GM, KB, freeboard, GZ at 13 heel angles | the same ship, to 1e-6 m |
| 150 s of the 'full' damage-control scenario | the same ship, tracks to 4e-7 |
| 900 s outcomes for all three scenarios | `scripts/verify.sh full`, identical verdicts |
| Full-precision offsets round trip | bit-exact against `makeHullFromStations` |
| Every line and byte truncation of a valid file | accepted exactly when it ends after a complete ship; nothing half built |

Each of the first four carries a negative control: a station offset moved by a
centimetre must break the static tolerances, and a breach 4% larger must break the
trajectory one. A comparison that cannot fail is not a measurement.

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
- ~~Cross-check the Craig–Bampton reduced model against the full model it was
  reduced from — this one is free and should run in CI.~~ **Done**, and it is
  free: `tests/test_reduction.cpp` builds the full model's own dense spectrum and
  asserts that the reduced frequencies bound it **from above** and fall
  monotonically as modes are added, which is stronger than closeness because a
  reduction can only stiffen. Static condensation is checked against
  `solidshell::solveStatic` on the same problem, to 2 × 10⁻¹⁰ m of a 0.31 m
  deflection. See `02-simulation.md` §3.

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
