# 04 — Multiplayer

Target: **co-op crew.** Several players aboard one ship, each with a station and a
job — bridge, engine room, damage control party, fire team. Optionally several
ships in one world (convoy, collision scenarios, salvage).

This is the architectural decision with the longest reach, which is why it is
settled before the systems it constrains get written.

---

## The choice: authoritative server, not lockstep

Deterministic lockstep is the obvious fit for a heavy simulation — send inputs
only, everyone computes the same world. It is rejected here for one reason: **the
simulation cost is not uniform across clients.** A player standing in a flooding
engine room needs a promoted particle fluid solver and an active FEM zone. A
player on the bridge does not. Lockstep would force every client to simulate every
other client's expensive local detail, and the frame budget does not survive that.

So: **one authoritative simulation**, running on a dedicated server or on a host
client, with clients receiving state.

Determinism is still built (see `01-architecture.md` §4) — it is what makes
client-side prediction reconcile cleanly, what makes replays exact, and what makes
a 10⁶-variable simulation debuggable at all. It is just not used for lockstep.

## What replicates, and how

The state falls into tiers with wildly different costs, and each gets a different
treatment.

| State | Size | Method | Rate |
|---|---|---|---|
| Ship rigid body (pos, orientation, velocities) | ~50 B | snapshot + interpolation | 20 Hz |
| Flooding network (water volume, air mass per compartment) | ~2 KB for 200 compartments | delta-compressed snapshot | 10 Hz |
| Opening states (doors, valves, breaches) | bitfield + areas | reliable events | on change |
| Machinery / systems state | few KB | delta snapshot | 10 Hz |
| Player avatars | ~40 B each | snapshot + prediction | 30 Hz |
| Fire / gas fields | large | **not replicated** — see below | — |
| FEM node state | very large | **not replicated** — see below | — |
| Particle fluids | very large | **not replicated** — see below | — |

The insight that makes this affordable: **the coarse state is tiny.** The entire
flooding condition of a 200-compartment ship — the thing that determines whether
she lives — is a few kilobytes. It is the *fields* that are enormous, and fields
do not need to be replicated because they are reproducible from their boundary
conditions.

### Fields are re-derived, not replicated

- **FEM**: the server replicates *damage events* — this element tore, this zone
  yielded, this is the resulting hole area and position. Clients run their own
  local FEM for visual deformation, seeded by those events. Client and server may
  differ by a millimetre in the exact shape of a dent; nobody can tell, and the
  hole area that governs flooding came from the server.
- **Fire and gas**: the server owns the authoritative per-compartment layer state
  (temperature, layer height, species concentrations) — small. Clients run local
  volumetric detail for rendering, driven by those boundary values.
- **Particle fluids**: entirely client-local. The server owns the water *volume*
  in each compartment; the client decides where the individual particles are.
  Two players in the same flooding room see slightly different splashes and the
  same rising water level.

This is the standard "replicate the summary, simulate the detail" pattern, and
this simulation happens to have an unusually clean split between the two.

## Prediction and reconciliation

- **Player avatars**: standard client-side prediction with server reconciliation.
  Complicated here by the fact that the floor is accelerating — the avatar's
  reference frame is the ship, so prediction runs in ship-local coordinates and
  the ship's own motion is interpolated separately. This is the same trick as
  moving platforms, applied to a platform with six degrees of freedom.
- **Ship motion**: not predicted from player input (there is no meaningful input
  latency on a rudder order — the ship takes 30 seconds to respond anyway).
  Interpolated between snapshots with a small buffer.
- **Interactions** — closing a door, opening a valve, starting a pump: predicted
  optimistically on the client for immediate feedback, corrected if the server
  disagrees. Doors are the common case and a 100 ms mispredicted door is
  acceptable; a mispredicted flooding rate is not, which is why the *consequence*
  is never predicted, only the *animation*.

## Late join and reconnection

A ship 20 minutes into a casualty has a large amount of accumulated state.
Late-joining clients receive a full baseline snapshot (compartment states, opening
states, damage event log, system states) and re-derive their local fields from it.
The damage event log is capped and compacted — a torn zone's history collapses to
its current geometry.

## Transport

UDP with a custom reliability layer: unreliable-sequenced for snapshots, reliable-
ordered for events, with delta compression against the last acknowledged baseline
per client. Snapshot deltas are generated from the reflection data described in
`01-architecture.md` §2, so a new replicated field costs an annotation and nothing
else.

Bandwidth target: **< 40 kbit/s down per client** in steady state, spiking during
active damage. This is achievable precisely because of the summary/detail split.

## Server cost

The server runs the authoritative sim without rendering, but *with* the expensive
solvers — it must, because the FEM decides hole sizes and the fire decides
bulkhead strength. Mitigation: the server runs solvers at the coarsest fidelity
that preserves the outcome (Tier-1 FEM plus promoted zones only where damage is
actually occurring, multi-zone fire only, no particle fluids at all). A dedicated
server for one 8-player ship should fit in ~4 cores.

## Anti-cheat

Co-op, so the threat model is mild. The authoritative server means a client cannot
fabricate ship state. Clients do get full knowledge of the ship's condition, which
is fine — there is no hidden information worth protecting in a co-op damage
control scenario.

## What this rules out

Being honest about the cost of this choice:

- **No 100-player ships.** The bandwidth is fine; the interaction contention and
  server sim cost are not. Design target is 2–16.
- **Rewind-based lag compensation is not used.** Rewinding a coupled multiphysics
  simulation is not practical. Nothing in this game needs sub-100 ms hit
  registration, so this costs nothing.
- **The host cannot be a potato.** Listen-server play puts the full authoritative
  sim on one player's machine. Dedicated server is the recommended configuration
  and the one that gets tested.
