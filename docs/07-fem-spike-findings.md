# 07 — FEM Spike: Results

The Phase 3 de-risking probe from `06-roadmap.md`, run on the target hardware
(GTX 1070 Ti, Pascal, 8 GB). Source in `tools/fem_spike/`, `engine/sim/fem.cpp`
and `engine/gpu/`. Reproduce with `./build/fem_spike`.

The question was whether adaptive full-3D tetrahedral FEM on a ship is a plan or
a wish. The answer is **both, depending on the element** — and that is a more
useful answer than either alone.

---

## 1. The formulation is correct

Explicit co-rotational linear tets, validated against a cantilever under its own
weight, where Euler-Bernoulli gives a tip deflection of qL⁴/8EI:

| Elements through thickness | Tets | Deflection | Error vs theory |
|---|---|---|---|
| 2 | 480 | 0.0801 mm | 63.6% |
| 4 | 3 840 | 0.1491 mm | 32.2% |
| 8 | 30 720 | 0.1962 mm | 10.8% |

Theory: 0.2199 mm. Refinement is uniform in all three directions, so element
aspect ratio is held constant and the study is not confounded by it.

Converging, from below, at roughly first order. That is textbook behaviour for
linear tetrahedra and it confirms the formulation. It also, immediately, is the
problem — see §4.

## 2. The GPU kernel computes the same thing as the CPU

Both paths run identical float arithmetic in the same order; the CSR gather
exists precisely so that nodal force accumulation is not at the mercy of warp
scheduling. Starting from a deformed, twisted state so the polar decomposition
is genuinely exercised:

| Steps | RMS velocity difference | Max position difference |
|---|---|---|
| 1 | 1.96 × 10⁻⁵ | 1.5 × 10⁻⁸ m |
| 10 | 1.45 × 10⁻⁴ | 3.0 × 10⁻⁸ m |
| 100 | 1.21 × 10⁻³ | 2.7 × 10⁻⁷ m |
| 1 000 | 4.61 × 10⁻³ | 6.0 × 10⁻⁷ m |

**One step agrees to 2 × 10⁻⁵ relative.** That is float rounding across a few
hundred operations with cancellation — the kernel is right.

**The growth from there is the finding.** The GPU is free to contract
multiply-add pairs differently, and an explicit scheme at the CFL limit amplifies
that seed difference in its highest-frequency modes, which sit at the stability
boundary. So:

> **The GPU FEM path is not bit-reproducible against the CPU, and will not be
> reproducible across drivers or vendors either.**

This is not fixable by being careful; it is a property of explicit dynamics.
It confirms, empirically rather than by assertion, the decision in
`04-multiplayer.md` to replicate damage *events* and let each client run its own
local FEM for visuals. Any design that assumed client and server FEM stay in
lockstep would have failed here.

It also bounds the determinism claim in `01-architecture.md` §4: determinism
holds for the CPU simulation path, not for GPU-accelerated field solvers.

*(A diagnostic worth recording: comparing peak nodal velocities reads over 200%
difference and looks alarming. It is meaningless. Switching gravity on at t = 0
is an impulse that excites modes oscillating with a period of ~4 steps, where an
explicit scheme has no accuracy at all. Sweeping the GPU step count from 198 to
202 does not reduce the difference, ruling out a step offset. Use integral or
low-frequency quantities to compare solvers, never the max norm on velocity.)*

## 3. Throughput

Steel plate, varying mesh density, 2000 steps per measurement, GPU timestamps:

| Tets | Nodes | ms/step | M element-updates/s | Wall ms per simulated second |
|---|---|---|---|---|
| 960 | 315 | 0.0083 | 116 | 3 423 |
| 11 520 | 2 745 | 0.0173 | 666 | 14 305 |
| 46 080 | 10 285 | 0.0786 | 587 | 65 005 |
| 230 400 | 45 225 | 0.5253 | 439 | 869 490 |
| 491 520 | 95 337 | 1.0461 | 470 | 1 731 384 |

**Peak ≈ 670 M element-updates/s; ≈ 450–470 M sustained at working set sizes.**
Small meshes are launch-overhead bound; the largest are memory bound on the node
gather.

For scale, the same kernel on one CPU core manages **4.5–6.7 M element-updates/s**.
The GPU is worth about 100 cores — or, put less flatteringly, about **4× the whole
24-thread CPU**. Useful, but not a silver bullet, and a well-threaded CPU fallback
is a viable degradation path rather than a joke.

The last column is the one that matters. It is dominated not by element count but
by the CFL limit: halving element size doubles the element count in each of three
directions *and* halves the timestep, so cost scales as h⁻⁴.

## 4. What this actually means for Phase 3

Working the numbers forward at 500 M element-updates/s, for steel
(c = 5172 m/s, dt = 0.5·h/c):

| Element size | Elements affordable at 1/50 slow motion | at 1/100 |
|---|---|---|
| 10 mm | ~24 000 | ~48 000 |
| 25 mm | ~60 000 | ~120 000 |
| 50 mm | ~120 000 | ~240 000 |

A realistic collision zone — 20 m of side shell, 10 m deep — at 10 mm resolution
with 2 elements through 20 mm plating is on the order of **2–4 million elements**.
That is 10³–10⁴× slower than real time. Not interactive, not even at heavy time
dilation.

At 50 mm elements the same zone is ~160 000 elements and runs at roughly 60× slower
than real time, which *is* usable with the time dilation the engine already plans
for. But 50 mm elements cannot resolve bending in 20 mm plate — and §1 shows
exactly how badly linear tets behave when under-resolved through the thickness:
**63% error at 2 elements through.**

So the two constraints close on each other, and the conclusion is not about
throughput at all:

> **Linear tetrahedra are the wrong element for ship plating.** They need many
> elements through the thickness to avoid locking, and that resolution is the one
> thing the CFL limit makes unaffordable.

### Revision to the plan

Tier 2 should not be uniformly tetrahedral. It should be:

- **Solid-shell or assumed-strain (EAS/ANS) shell elements for plating** — one
  element through the thickness, no locking, and a timestep set by the *in-plane*
  element size rather than the plate thickness. That last point is worth
  emphasising: it is a 5–10× timestep win on top of the element count win,
  because plate thickness is what currently drives the CFL limit.
- **Tetrahedra retained where the geometry is genuinely three-dimensional** —
  castings, engine seats, thick brackets, and the crush zone itself once plating
  has folded and is no longer thin.
- **Adaptive promotion tet-ward** as a shell element crumples past the point where
  shell kinematics hold.

This is more work than a uniform tet solver, and it is the difference between a
20 m damage zone being affordable and not.

### Optimisation headroom, separately

The kernel is unoptimised. The polar decomposition alone does four 3×3 matrix
inversions per element per step. An analytic or quaternion-based polar
decomposition is roughly 3× cheaper, and the node gather would benefit from
shared-memory staging. Expect **2–3×** from that work — real, but an order of
magnitude short of what the element change buys.

Hardware scaling is worth about the same as the element change: a 4090-class card
is 10–15× a 1070 Ti in FP32 throughput.

## 5. Verdict

Phase 3 is viable, at ~18 engineer-months as estimated, **provided the element
technology changes**. Specifically:

- ✅ Explicit GPU FEM at useful throughput on modest hardware — confirmed.
- ✅ Adaptive zone promotion with a reduced global model — unchanged, and now
  more clearly necessary.
- ✅ Damage events rather than field replication for multiplayer — confirmed
  empirically.
- ⚠️ Uniform linear tetrahedra at plate-resolving resolution — **ruled out.**
  Replace with solid-shell elements for plating.
- ⚠️ Determinism does not extend to GPU field solvers. Bound the claim.

The spike cost a day and moved a 18-engineer-month bet from "assumed" to
"measured, with one architectural change identified before any of it was built".
That is what it was for.
