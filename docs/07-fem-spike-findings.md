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

  > **Correction, from building it — see §6.** The timestep claim is wrong as
  > written, and the mechanism matters more than the number. A *solid-shell* keeps
  > its through-thickness stretch degree of freedom — deliberately, because a crush
  > zone needs the plate to be able to thin — so its highest frequency is the
  > thickness dilatational mode and the stable step is `t / c_p` **however large
  > the element is in plane**. Measured: `dt·c_p/t = 0.999` at in-plane sizes from
  > 5t to 50t, i.e. flat. The step is *not* set by the in-plane size. The 5–10×
  > figure survives by accident, because it is the right comparison against a tet
  > mesh with 8 elements through the thickness (measured 6.7×), but anyone sizing
  > elements from the stated mechanism would over-predict the affordable step by
  > the in-plane aspect ratio, which is 25× at 50 mm elements on 20 mm plate. A
  > timestep governed by in-plane size needs a classical shell element with no
  > thickness stretch, or selective mass scaling on the thickness mode.
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
  Replace with solid-shell elements for plating. **Done — §6.**
- ⚠️ Determinism does not extend to GPU field solvers. Bound the claim.

The spike cost a day and moved a 18-engineer-month bet from "assumed" to
"measured, with one architectural change identified before any of it was built".
That is what it was for.

---

## 6. The solid-shell element — **implemented**

`engine/sim/solid_shell.{hpp,cpp}`, checked by `tests/test_solid_shell.cpp`. The
element §4 specified, built and measured.

### What was chosen, and what it costs

An eight-node hexahedron with three translational degrees of freedom per node —
no rotations, no director — so it stacks against tetrahedra at an interface for
free, which is what the "promotion from shell to tet as an element crumples" plan
in `02-simulation.md` §3 needs. Three cures, all parameter-free:

| Locking mode | Cure | What it buys, measured |
|---|---|---|
| Transverse shear | **ANS** (Dvorkin–Bathe), covariant shear sampled at the mid-edge points | the difference between 0.018 and 0.996 of the closed form at L/t = 100 |
| Thickness (Poisson) | **EAS**, 3 modes carrying `E_ζζ = ζ, ξζ, ηζ` | 0.816 → 1.000: exactly the ratio of the plane-stress modulus to the oedometer one, 1.225 |
| In-plane shear of a distorted element | **EAS**, 4 Simo–Rifai modes | 0.579 → 0.983 in in-plane bending |
| Curvature thickness | **ANS** (Betsch–Stein), `E_ζζ` sampled at the four in-plane corners | free on a flat element; needed on a warped one |

Seven enhanced parameters, condensed out at element level. **Reduced integration
with hourglass control was rejected**: it is cheaper and would pass the patch test
too, but the hourglass stiffness is a tuned coefficient with no physical value and
no measurement that sets it — too small and the mesh develops zero-energy modes
under exactly the loads this solver exists to compute, too large and it stiffens
the bending mode it is meant to leave alone. ANS and EAS have no free parameters,
so there is nothing to tune wrongly and nothing to re-tune per material.

The element stiffness has **exactly six zero eigenvalues and eighteen positive
ones** on a cube, on a plate element and on one distorted in all three directions
— so the enhanced modes have not bought a spurious mechanism, which is the failure
a deflection test cannot see.

### Cost, one core, measured

| | |
|---|---|
| Element stiffness formation | **21.1 µs** — once, when a zone is promoted |
| The same for a plain hex | 13.4 µs, so the assumed strains cost **1.57×** |
| Internal force per step | **267 ns** |
| `fem.cpp` linear tet internal force per step | 129 ns, so one solid-shell = **2.1 tets** |
| Explicit step, 20 mm plate | 3.25 µs (`t / c_p`, in-plane size irrelevant) |

Per element the solid-shell is twice a tet. That is not the comparison that
matters. Per **square metre of 20 mm plating per simulated second**, single core:

| | elements/m² | steps/s | cost |
|---|---|---|---|
| Linear tets at 2.5 mm (8 through the thickness, the spike's 11%-error mesh) | 7.7 M | 4.1 M | **4.1 × 10⁶ s** |
| Solid-shell at 50 mm, one through the thickness | 400 | 3.3 × 10⁵ | **37 s** |

**A factor of 1.1 × 10⁵.** A 20 m × 10 m collision zone is 200 m² of shell, so one
simulated second costs ~2 core-hours as tets and ~2 core-minutes as solid-shells —
about 5 minutes of wall time per simulated second on the 24-thread CPU, or well
under a minute on the GPU at the throughput §3 measured. That is inside the time
dilation the engine already plans for, and it is what makes the 20 m damage zone
affordable rather than theoretical.

Stiffness formation at 21 µs means promoting a 10⁵-element zone costs ~2 s on one
core, ~0.1 s threaded. That is a hitch at promotion, not a per-frame cost, and it
is the number to watch if zones are promoted during play rather than at impact.

### What it was checked against

Closed forms throughout, not eyeballed output. 207 assertions:

- **The patch test**, exactly, on a distorted patch: displacement to 1.1 × 10⁻¹⁵
  and constant stress to 2 × 10⁻¹³. Distortion in-plane is arbitrary.
- **Rigid body motion** carries no force, including finite rotations to 3 rad; and
  the stronger statement, **frame indifference** — a 25% stretch rotated 1.7 rad
  gives the rotated force to 8 × 10⁻¹³ of it.
- **Cylindrical bending** against `PL³/3Db`, converging at **order 2.1–2.4** and
  reaching the closed form to 1 × 10⁻⁵.
- **Plates**: the Navier double series (summed, not quoted) to 0.04%; the clamped
  square plate at coefficient 0.001265 against Timoshenko's tabulated 0.00126;
  supported and clamped strips at `5qL⁴/384D` and `qL⁴/384D` to 0.07%.
- **σ_zz relaxes to machine zero** through a bent plate, which is the whole point
  of the enhanced thickness strain and is invisible to any deflection test.

### Locking, demonstrated on one mesh

The same 8 × 2 × 1 mesh, the same load, the same solver, ratio to the closed form:

| L/t | solid-shell | ANS without EAS | plain 8-node hex | linear tets |
|---|---|---|---|---|
| 10 | 1.003 | 0.820 | 0.569 | 0.199 |
| 50 | 0.996 | 0.813 | 0.067 | 0.024 |
| 100 | 0.996 | 0.813 | 0.018 | 0.006 |
| 500 | 0.996 | 0.813 | **0.0007** | **0.0003** |

At the slenderness of real plating a plain hex is 1 400× too stiff and a linear
tet 3 800× too stiff, on a mesh where the solid-shell is within 0.4%. The middle
column is the honest accounting of what each cure is for: ANS alone fixes shear
locking and leaves the 22.5% thickness-locking penalty untouched.

### Limits — measured, not guessed

1. **The element must be close to prismatic through its thickness.** The ANS
   interpolation is exact only when the top and bottom faces are parallel and the
   thickness direction is straight. Where they are not, it introduces a parasitic
   strain, and bending squares it into a spurious stiffness: measured at
   **≈ 90 × (offset/t)²** excess stiffness, so a 5% face offset costs 24% and a 10%
   offset costs 90%. It is a consistency error and not an inconsistency —
   refinement at a fixed element *shape* converges — but on a fixed mesh the
   constant is large. **Practical rule: keep the thickness direction within a few
   degrees of the surface normal, and change plate thickness at a seam rather than
   across an element.**
2. **The warped patch test fails, in proportion to the warp**: 0.7% stress error at
   2% warp, 17% at 20%. A plain hex passes it exactly, which locates the fault in
   the ANS interpolation rather than in the mesh. On a smooth hull the warp per
   element falls as O(h), so consistency is recovered under refinement.
3. **Trapezoidal distortion in-plane costs about a third** (0.983 → 0.644 in
   in-plane bending); parallelogram distortion costs 5%. This is the known
   behaviour of the Simo–Rifai enhanced modes, which are exact for a constant
   Jacobian.
4. **The timestep is thickness-governed, not in-plane-governed** — the correction
   to §4 above.
5. **Stress converges one order below displacement** (measured order 1.0 against
   `Mz/I`), which is ordinary for a trilinear element but matters for a yield or
   fracture criterion evaluated on it. Budget a finer mesh for stress than for
   deflection, or recover stress by patch averaging.
6. **No GPU path, no plasticity, no fracture, no `StructuralMesh` consumer yet.**
   The element is elastic and co-rotational: it is the element technology, not the
   Tier-2 solver. A float GPU kernel can be derived from it the way `fem_gpu.cpp`
   was derived from `fem.cpp`, and will inherit §2's reproducibility bound.
7. **It is in `double`, where `fem.cpp` is in `float`.** The identities that
   establish an element is correct — the patch test, rigid-body invariance, the
   rank of the stiffness — are exact, and in float their noise floor sits at 1e-6,
   the same order as the defects they exist to catch.

### What mutation testing found

52 mutants, each a single plausible edit. The first pass killed 41 and **the 11
survivors were all real holes**, each in a part of the element no test reached:

- Four of the seven enhanced parameters — the in-plane Simo–Rifai set — were
  **never exercised**, because nothing loaded the element in its own plane. Moving
  one of them onto the wrong strain component changed nothing.
- The enhanced parameters are recovered from the displacement inside
  `elementStress`, and **discarding them entirely, or recovering them with the
  wrong sign, passed everything** — the condensed stiffness is identical either
  way, so no deflection test can see it. Only `σ_zz` can.
- The polar decomposition could be **truncated to one iteration** and pass, because
  the rigid-body test rotated an *undeformed* element, where the deformation
  gradient is already a rotation and Higham's iteration returns it untouched.
- The lumped mass could be split evenly instead of by row sum, and the pressure
  load could ignore its shape functions, because **every mesh under test was
  rectangular**, where the two agree exactly.
- The inverted-element guard could look at one corner instead of eight.
- The internal force's **sign** was untested: it was only ever asked for zero.

Adding tests for those — in-plane bending on trapezoids, the through-thickness
stress state, frame indifference of a 25%-stretched element, mass and pressure
resultants on distorted geometry, a folded element, and the sign and equilibrium
of the internal force — took the suite from 152 to 207 checks and killed **51 of
52**. The single survivor is genuinely equivalent: loosening the polar
decomposition's convergence tolerance from 1e-16 to 1e-6 moves the frame
indifference residual from 8.3e-13 to 5.9e-13, which is to say it moves it around
inside the rounding floor, because Newton convergence is quadratic.

Two defects in the code were found this way or by the assertions it prompted, both
of which would have shipped green on a plausible suite: a **duplicated node in the
mesh generator**, caught immediately by the Jacobian guard, which is the argument
for having the guard; and `solveStatic` **failing open** on an inverted element
when every degree of freedom happened to be prescribed, so there was no work to do
and it reported success.

Two more were in the *tests*, and both were caught by the habit of asking whether
the expectation or the design is wrong before loosening anything: a uniaxial
**strain** compared against a uniaxial **stress** modulus — out by exactly
`(1-ν)/((1+ν)(1-2ν))`, 34.6%, which is what identified it — and a convergence order
computed from two numbers that were both already at the noise floor, reporting a
confident 4.03 that meant nothing.
