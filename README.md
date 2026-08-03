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
./build/shipsim_tests                       # 29 validation checks against closed-form answers
./build/shipsim --scenario=none             # 120 m ferry, holed, nobody does anything
./build/shipsim --scenario=doors            # close the watertight door
./build/shipsim --scenario=full             # full damage control response
./build/shipsim --scenario=full --csv=run.csv
```

## What slice 1 already does

A 120 m ro-pax ferry, 16 compartments over four levels, is holed by a 2.4 m²
tear in the starboard shell 2.5 m below the waterline. From there the simulation
takes over. Three runs, same ship, same damage, different decisions:

| Scenario | Action taken | Outcome |
|---|---|---|
| `none` | nothing | founders at **t+555 s**, GM −2.6 m, 3251 t of floodwater |
| `doors` | close the watertight door at t+45 s | **capsizes at t+330 s** — *sooner* |
| `full` | door, pumps, early counterflood, secure the vehicle deck | **survives** — steady 6.3° list, GM +1.7 m |

The middle row is not a bug and it is the reason this is worth building. Closing
the door stops the two engine rooms cross-equalising, so ~700 t of asymmetric
floodwater becomes a list; the list pushes the starboard vehicle-deck openings
under; water spreads across a 100 × 19 m undivided deck; the free surface moment
destroys GM and she goes over. Symmetric flooding drowns her more slowly than
asymmetric flooding rolls her. Nobody scripted that — it falls out of the
integrals.

The same physics also declines to let you cheat. In an earlier run the
counterflooding attempt failed, and it failed for the right reason: by the time
the valves were opened the ship had listed far enough to lift the port sea
suctions clear of the water. Counterflooding is a thing you do *early* or not at
all, and the simulation will not tell you that in a tooltip.

## How the flooding model works

Three deliberate refusals of the standard shortcuts:

**Floodwater is mass, not a correction term.** Naval architecture normally handles
loose water with a tabulated free-surface correction subtracted from GM. Here the
water inside each compartment is re-levelled against gravity every tick, its true
volume and centroid computed by clipping the compartment mesh with the free
surface plane, and carried as mass at that centroid. Free surface effect is then
not modelled at all — it simply happens, because a few hundred tonnes of water
slides to the low side. The test suite confirms the emergent loss matches the
classical ρ·i/Δ to within 15%.

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

## Repository layout

```
engine/core/     math, closed-mesh volume integration
engine/sim/      compartments, orifice network, 6-DOF rigid body
game/prototype/  the ferry definition and the headless scenario driver
tests/           closed-form validation of the numerical core
docs/            the plan
```

## Conventions

SI units throughout, no exceptions. Body frame is +x forward, +y to port, +z up,
origin at midship on the baseline. Heel is positive starboard-down, trim positive
bow-down. Angles are radians everywhere except at display boundaries.
