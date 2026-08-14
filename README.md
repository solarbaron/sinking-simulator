# shipsim

A ship simulator where the ship is actually simulated.

Not a vehicle with a health bar and a sinking animation — a hull with material
properties, a structure that yields and tears under real stress, compartments
that flood through real orifices at real rates, air that compresses and finds its
way out, fires that heat steel until it loses strength, and a rigid body whose
stability is whatever those things add up to. If the ship capsizes it is because
the righting arm went negative, and you can go and read why.

**Status:** slice 1 (progressive flooding and damaged stability) is implemented,
validated and running headless. Everything else is specified in `docs/`.

---

## Build and run

Needs a C++23 compiler and CMake ≥ 3.24. No external dependencies for the current
slice.

```sh
cmake -S . -B build -G Ninja
ninja -C build
./build/shipsim_tests                       # 200380 validation checks against closed-form answers
./build/shipsim --scenario=none             # 120 m ferry, holed, nobody does anything
./build/shipsim --scenario=doors            # close the watertight door
./build/shipsim --scenario=full             # full damage control response
./build/shipsim --scenario=full --csv=run.csv
./build/shipsim --ship=ships/ferry.ship     # the same ferry, read from a text file
./build/fem_spike                           # explicit tet FEM: validation + GPU benchmark
```

`fem_spike` needs Vulkan and `glslc`; without them CMake skips that target and
everything else still builds.

## What slice 1 already does

A 120 m ro-pax ferry — 8984 t, Cb 0.66, intact GM 2.00 m, 18 compartments carved
out of the hull form over four levels — is holed by a 2.4 m² tear in the
starboard shell 2.5 m below the waterline. From there the simulation takes over.
Three runs, same ship, same damage, different decisions, each `--duration=1800`:

| Scenario | Action taken | Outcome |
|---|---|---|
| `none` | nothing | GM negative at t+690 s, lolls to **53° by t+1800 s**, 6311 t aboard at the end |
| `doors` | close the watertight door at t+45 s | **capsizes at t+930 s** |
| `full` | door, pumps, early counterflood, secure the vehicle deck | **lost** — GM −3.23 m under 6 mm of water on the vehicle deck, 9.3° list, 1556 t |

**The duration is part of the figure**, which is why it is written above the
table: `full` is still flooding at the end, so a run of a different length
publishes different numbers and neither is stale. `docs/06-roadmap.md` quotes the
same run at the gate's own t+900 s, where she is at −3.30 m, 7.2° and 1442 t.

The middle row is not a bug and it is the reason this is worth building. Closing
the door stops the two engine rooms cross-equalising, so the asymmetric
floodwater becomes a list; the list pushes the starboard vehicle-deck openings
under; water spreads across a 100 × 19 m undivided deck — 3950 t of it by the
end — and the free surface moment finishes her. Symmetric flooding drowns her
more slowly than asymmetric flooding rolls her. Nobody scripted that; it falls
out of the integrals.

The bottom row used to read "survives", and it changed on the day GM stopped
being finite-differenced at a fixed ±0.03 rad. Every tonnage, list and draft in
it is the same; what moved is the number that judges them. Eleven tonnes of water
on that deck is a **6.6 mm** layer, and it pulls off the high side at 7.07e-4 rad
— so a slope taken at forty-two times that angle sees **6%** of that deck's free
surface and reports +1.38 m where the metacentric height is −3.23 m. The full
damage-control response buys her a great deal and it does not buy her positive
stability. See `docs/02-simulation.md` §1.

Every figure in that sentence now comes from `shipsim --scenario=full
--duration=1800 --gm-detail` and is gated. The angle used to read 6.8e-4 rad,
computed from the deck's **nominal** 100 × 19 m; the deck the ship actually
carries is 1868.4 m² with a **mean** breadth of 18.684 m, and a bounding box
would have said 20.00 m — the deck at its widest, which is the one width no part
of the calculation wants. That is a 4% error in the angle and it changes nothing
here, but it is the kind that survives precisely because it looks derived.

The same physics declines to let you cheat. In an earlier run the counterflooding
attempt failed, and it failed for the right reason: by the time the valves were
opened the ship had listed far enough to lift the port sea suctions clear of the
water. Counterflooding is a thing you do *early* or not at all, and the
simulation will not tell you that in a tooltip.

Trapped air is doing real work too. In the `doors` run, the final compartment
table shows spaces the sea never directly reached:

```
compartment            gross m3   fill %    P kPa   water t
fwd_hold_s                  965     20.6    127.7       194
wing_tank_fwd_s             108     32.5    150.1        35
aft_hold_s                 1232     25.9    123.2       311
```

Each vent goes under at its own angle, because each sits at its own half-breadth:
the starboard wing tank's air pipe at about **32°** of heel and the forward
hold's vent at about **42°**. Past that the sea starts down the vent — and stops,
because the air it has to displace can only leave by the same submerged pipe.
Each one settles at a fill fraction set by where its air pressure balances the
outside head, 1.2 to 1.5 atmospheres. That behaviour was never written down
anywhere; it is two lines of orifice physics meeting Boyle's law. The aft hold is
in the table for a different reason — the sea reaches it through an unsealed
cable transit from the engine room within the first minute — but its air is
trapped by the same submerged vent.

## How the flooding model works

Three deliberate refusals of the standard shortcuts:

**Floodwater is mass, not a correction term.** Naval architecture normally handles
loose water with a tabulated free-surface correction subtracted from GM. Here the
water inside each compartment is re-levelled against gravity every tick, its true
volume and centroid computed by clipping the compartment mesh with the free
surface plane, and carried as mass at that centroid. Free surface effect is then
not modelled at all — it simply happens, because a few hundred tonnes of water
slides to the low side. The test suite confirms the emergent loss matches the
classical ρ·μ·i/Δ to **1.4e-6** on a box barge and 2e-4 on the ferry's own
vehicle deck, with the frozen and the free ship held at the same displacement and
the same centre of gravity so that nothing but the water's freedom to move
differs between them. (An older check quotes 15%; its two barges do not have
quite the same KG.)

**Air is a real compressible species.** Every compartment tracks air mass and
pressure. A sealed space stops flooding when its air pressure balances the
external head — the reason a capsized hull floats for hours. `testTrappedAirArrests-
Flooding` checks the result against Boyle's law. Vents, air pipes and open hatches
are edges in the same network the water uses.

**One network, two phases.** Breaches, doors, hatches, vents, pipes and cross-flood
ducts are all orifices. Flow is `Cd·A·√(2Δp/ρ)`. What moves through an opening is
whatever is sitting against it on the high-pressure side — which is why the same
hole vents air when it is above the internal waterline and admits sea when it is
below, with no special-casing anywhere.

The geometric core underneath is one routine: given a closed triangle mesh and a
plane, return the volume and centroid below the plane. Buoyancy asks it of the
hull against the sea; floodwater asks it of a compartment against its own free
surface. It uses the divergence theorem about a reference point *on* the cutting
plane, so cap tetrahedra vanish identically and no cap geometry is ever built —
correct for any mesh topology, including planes that cut a compartment into
several separate loops.

**Compartments are carved, not authored.** Each one is the hull interior clipped
by its own bulkhead and deck planes, using a mesh boolean that welds the cut
edges, chains them into loops and caps them by ear clipping. So a wing tank
tapers into the turn of the bilge and a forepeak narrows into the stem, instead
of being a box that pokes out through the shell. `testSubdivisionTiles` cuts a
hull into a grid of cells and checks the volumes sum back to the whole.

That work also turned up a bug that had been silently poisoning every number in
the project: `makeHullFromStations` wound its keel and deck strips opposite to
the shell plating, so the hull was never a closed, consistently-wound manifold.
Volume integrals on such a mesh do not fail — they return plausible, wrong
answers that depend on where you put the reference point. The hull was overstating
displacement by 40%. There is now an `isClosedManifold` check and `Ship::validate()`
runs it at load, because this class of error is invisible until you go looking.

## The FEM spike

Phase 3 — deformation, tearing, structural failure — is the largest and riskiest
part of the plan, so it got probed before anything was built on it. There is now
an explicit co-rotational tetrahedral FEM with a CPU reference and a Vulkan
compute back-end, and it has been measured on the target GPU.

- **The formulation validates.** A cantilever under its own weight converges to
  Euler-Bernoulli theory under uniform refinement: 63.6% error at 2 elements
  through the thickness, 32.2% at 4, 10.8% at 8.
- **The GPU kernel is correct.** One step from identical deformed state agrees
  with the CPU to 2 × 10⁻⁵ relative.
- **It is not reproducible over long runs.** That agreement decays to 5 × 10⁻³
  over a thousand steps, because the GPU contracts multiply-add pairs differently
  and an explicit scheme at the CFL limit amplifies it. Not a bug — a property of
  explicit dynamics, and it confirms the decision to replicate damage *events*
  rather than field state across the network.
- **Throughput: 450–670 M element-updates/s** on a GTX 1070 Ti. About 100× one
  CPU core, or 4× the whole 24-thread CPU.

And the finding that changed the plan: **uniform linear tetrahedra are the wrong
element for ship plating.** They lock in bending, so they need many elements
through the thickness — and the explicit stability limit is set by the smallest
element dimension, so those same elements collapse the timestep. Cost scales as
h⁻⁴ and the two constraints close on each other. A realistic 20 m collision zone
at plate-resolving tet resolution is 10³–10⁴× slower than real time.

The fix is solid-shell elements for plating, with tets kept for genuinely
three-dimensional regions — and that element now exists
(`engine/sim/solid_shell.{hpp,cpp}`): an eight-node hex with assumed natural
strains and seven enhanced strain parameters, one element through the plate. On
the same mesh at the slenderness of real plating it lands within 0.4% of the
closed form where a plain hex is 1 400× too stiff and a linear tet 3 800× too
stiff, and a square metre of 20 mm plating costs 37 core-seconds per simulated
second instead of 4.1 million. It also corrected the spike: the stable timestep is
set by the plate thickness after all, not by the in-plane element size, because a
solid-shell keeps its through-thickness stretch. Full numbers, limits and what
mutation testing found in
[docs/07 — FEM spike findings](docs/07-fem-spike-findings.md).

That element now has a consumer (`engine/sim/zone.{hpp,cpp}`): a patch of a
ship's own plating meshed into solid-shells around an impact, driven by a rigid
punch, solved explicitly with plasticity and ductile failure, and reporting which
plating panels tore — straight into the flooding network. A three-metre zone on the
reference ferry is 224 elements and a 0.037 s ram through it is 4.5 s of wall time on 23
threads, bit-identical however many of them run it. Its binding limit turned out to be
geometric rather than numerical: the elements are *exactly* prismatic on flat
plating and 319% too stiff in bending across this hull's shoulder, and the cure is a
finer girth layout in the scantlings rather than a finer mesh. Run against
`indentation.hpp`'s rigid-plastic membrane model on the same bay the two agree to a
factor of two, and the factor is predicted by the hardening curve and the
plane-strain constraint the membrane model does not have — which is also what
identified an order-of-magnitude error in which spacing that model was being handed
as its span.

And it now decides for itself when to exist (`engine/sim/promotion.{hpp,cpp}`).
The beam model is on everywhere and a zone is promoted only where the response is
*local* — a station standing above the ship's own median utilisation, or a contact
patch past the struck bay's plastic collapse pressure — because a Tier-2 element
costs 4.0 core-seconds per simulated second and that is linear in the number of
zones, so a criterion that fires everywhere is unaffordable rather than merely
noisy. A promoted zone is handed the girder's own stress to start from, which on
the ferry's side turns out to make her *stiffer* rather than weaker, and hands back
what it lost as a thinner ship that `girder.hpp`, `buckling.hpp` and `collapse.hpp`
read unchanged. Wiring the two together found a real defect in the third: a
progressive-collapse sweep sized from first yield never reached its peak on a
damaged section, and reported a hull girder that grew fifteen times stronger when
material was removed from it.

Between the beam and the zone there was, until recently, nothing — so a whole
hold, a superstructure, or the region between two bulkheads had no structural
answer at all. `engine/sim/reduction.{hpp,cpp}` is that middle tier:
**Craig–Bampton component mode synthesis**, keeping every interface degree of
freedom exactly and carrying the interior on a handful of fixed-interface modes,
with the symmetric eigensolvers written out rather than taken from a library. It
is unusually well supplied with exact answers and they are asserted as identities:
static condensation is exact at the interface for *any* load (2 × 10⁻¹⁰ m of a
0.31 m deflection against an independent solve), a free component keeps exactly
six zero eigenvalues, and the reduced frequencies come down **from above**,
monotonically, because a reduction can only stiffen. On the same patch of the
ferry's side it costs 0.4 core-seconds per simulated second against Tier 2's 1155
— 1200× — and what it buys with that is nothing nonlinear at all: it cannot yield,
tear, buckle or rotate, so it also reports when it has stopped being valid and the
region has to be promoted. Three published figures in the plan turned out to be
optimistic and are corrected there.

The middle tier now *drives* the fine one (`engine/sim/coupling.{hpp,cpp}`). A
zone's perimeter used to be clamped — the plating outside it could not move,
whatever happened inside — and it now follows a reduced model of the structure
round it, with a torn zone handed back to that model as a mesh with the dead
elements deleted so the surroundings feel the damage. Against the same plate meshed
and solved in one piece, the coupled zone reproduces the monolithic answer to
**10⁻¹⁵ m of a 2 × 10⁻⁴ m peak**, at zero modes and at twelve alike, where the
clamped zone it replaces is 74% out and 433× too stiff. That it does not merely
*approach* the monolithic answer is the point: static condensation is exact at an
interface, so the coupling is an identity and is tested as one — and the natural
convergence study, sweeping mode count for an improvement, would have measured
nothing at all.

The spike cost a day and moved an 18-engineer-month bet from assumed to measured,
identifying one architectural change before any of that phase existed to rewrite.

## Where this is going

`docs/` contains the full plan. In short:

- **[01 — Architecture](docs/01-architecture.md)** — job system, ECS, multi-rate
  solver scheduling, determinism.
- **[02 — Simulation](docs/02-simulation.md)** — the physics roadmap: seakeeping,
  adaptive tetrahedral FEM with plasticity and tearing, fire and thermodynamics,
  particle fluids, aerodynamics, machinery, ice, cargo.
- **[03 — Renderer](docs/03-renderer-audio.md)** — Vulkan 1.3 targeting Pascal-class
  hardware, ocean, volumetrics, interior water.
- **[04 — Multiplayer](docs/04-multiplayer.md)** — authoritative crew co-op.
- **[05 — Data & validation](docs/05-data-modding-validation.md)** — ship format,
  editor, and how each subsystem gets checked against real data.
- **[06 — Roadmap](docs/06-roadmap.md)** — phases, dependencies, honest costs.
- **[07 — FEM spike findings](docs/07-fem-spike-findings.md)** — what the
  measurements said and what changed because of them.

## Working on it

```sh
./scripts/install-git-hooks.sh   # once: conflict markers, stray binaries, worktrees
./scripts/verify.sh              # quick    build + the test suite     ~110–150 s
./scripts/verify.sh full         # + rebuild, GPU, scenarios, figures     ~2400 s
./scripts/verify.sh sanitize     # + AddressSanitizer and ThreadSanitizer
```

Those are measured, on a box that was busy with other work — which it usually is,
and `quick` came out at 112 s and 147 s an hour apart — and they grow with the
suite, so re-measure rather than trusting them. Most of `full` is the suite
again: six repeat runs for flakiness, and `scripts/check-figures.sh` re-running
every tool this repository quotes a number from, including the ones on this page.

**Warnings are failures.** The build is `-Wall -Wextra -Wpedantic` and has been
warning-clean since Phase 0, so one warning means the signal is decaying. Note
that an incremental build cannot see a warning in a file it did not recompile,
which is why `full` compiles everything from scratch in a throwaway directory.

`CLAUDE.md` carries the conventions, the testing philosophy, the settled
architecture decisions, and a table of every defect that shipped green on its
functional tests along with the instrument that actually caught it. That table is
the most useful page in the repo.

## Repository layout

```
engine/core/     math, closed-mesh volume integration, mesh booleans
engine/sim/      compartments, orifice network, 6-DOF rigid body, tet FEM
engine/gpu/      Vulkan compute back-end and shaders for the FEM
game/prototype/  the reference ferry in C++ and the headless scenario driver
ships/           ships as data: the same ferry, in the text definition format
tools/fem_spike/ FEM validation and GPU benchmark
tests/           closed-form validation of the numerical core
docs/            the plan
```

## Conventions

SI units throughout, no exceptions. Body frame is +x forward, +y to port, +z up,
origin at midship on the baseline. Heel is positive starboard-down, trim positive
bow-down. Angles are radians everywhere except at display boundaries.
