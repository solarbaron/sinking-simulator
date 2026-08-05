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
simulated second costs **~8 × 10⁸ core-seconds as tets and ~7 400 core-seconds as
solid-shells** — about 5 minutes of wall time per simulated second on the
24-thread CPU, or well under a minute on the GPU at the throughput §3 measured.
That is inside the time dilation the engine already plans for, and it is what
makes the 20 m damage zone affordable rather than theoretical.

> **Correction.** That sentence previously read "~2 core-hours as tets and
> ~2 core-minutes as solid-shells", and **both figures were wrong** — 200 × 37 s is
> 7 400 core-seconds, which is two core-*hours*, and 200 × 4.1 × 10⁶ s is 26
> core-*years*, not two core-hours. The wall-time figure beside them was right all
> along (7 400 / 24 = 5.1 minutes), which is why the error survived: the number
> anyone would sanity-check was correct and the two feeding it were not. Found while
> costing the elastoplastic path against it.

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
6. **No GPU path and no `StructuralMesh` consumer yet.** A float GPU kernel can be
   derived from this the way `fem_gpu.cpp` was derived from `fem.cpp`, and will
   inherit §2's reproducibility bound. **Plasticity and ductile failure are no
   longer missing** — see §7 and `02-simulation.md` §3 — but the element is still
   the element technology rather than the Tier-2 solver.
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

---

## 7. Plasticity and ductile failure — **implemented**

`engine/sim/plasticity.{hpp,cpp}`, hooked into the element by
`solidshell::elementPlasticUpdate`, checked by `tests/test_plasticity.cpp`.
`02-simulation.md` §3 records the formulation, the failure criterion, the cost
table and the limits. This section records what the *instruments* found, because
that is the part that does not fit in a design document.

### The step-independence test is the one that pays

Backward Euler on the flow rule is exact for radial loading: the consistency
condition fixes the plastic flow from the total deviatoric strain and from nothing
else. So one step and ten thousand steps to the same final strain must give the
same stress, and they do, to 10⁻¹². A forward-Euler update does not have that
property and no accuracy test at a fixed step size would tell the two apart.

It is paired with a **negative control** — a deliberately non-proportional path,
which must come out step-dependent — because the first thing a step-independence
test does when the plasticity is broken is pass.

### What mutation testing found

68 mutants, each a single plausible edit. The first pass killed 58; **six of the
ten survivors were real holes**, and one of the tests written to close them found a
defect in the code:

- **`elasticStress` was reachable by no test at all.** Dropping its deviatoric
  split — returning `K tr(ε) + 2μ ε_ii` instead of `K tr(ε) + 2μ(ε_ii − tr/3)` —
  passed everything, on a public entry point. It is now tied both to `elasticModuli`
  computed by a different route and, bit for bit, to the elastic branch of the
  return map.
- **The element's size could be taken off its bottom face**, because every element
  under test was prismatic and on a prism the bottom face and the mid-surface have
  the same area. A tapered element — 40 mm at the bottom, 60 mm at the top —
  distinguishes them, and its volume is the exact integral of `(a + (b−a)s)²`, which
  2×2×2 Gauss reproduces exactly because that is a quadratic.
- **`initialisePlasticState` could decline to clear the history it was
  initialising**, because every caller happened to hand it a freshly constructed
  state. That is exactly the condition under which a re-promoted zone inherits the
  damage of the last collision.
- **An element could call itself torn on its first dead integration point.** Every
  tearing test until then strained the element uniformly, where all eight points
  die on the same step and "any" and "all" agree.
- **The failure plane's normal was only ever asked for on an axis-aligned pull**,
  where the stress tensor is already diagonal, the Jacobi sweep has nothing to do
  and the eigenvectors come back as the identity however badly they are
  accumulated. Freezing them entirely passed. The same tear, rotated 0.9 rad about
  an arbitrary axis, kills it.
- The suite went from 1 007 to 1 045 checks and from 58 to 65 kills.

**The defect that came out of it.** Writing the partial-failure test exposed
something no assertion had been aimed at: with four of eight integration points
gone, the enhanced-strain problem `∫Gᵀσ dV = 0` loses rank, the Newton stops
converging, α wanders, and the surviving points were driven to a triaxiality
*below the damage cutoff* — their damage froze at 0.78 while their plastic strain
ran on from 0.49 to 0.89. **The element stopped tearing.** An element now drops its
enhanced modes the moment any point fails, and re-runs the step in which the point
died so that nothing ill-posed is committed.

**The three survivors are argued equivalent, and one of them measured so.**
Loosening the return map's scalar-Newton step tolerance from 1e-15 to 1e-6 leaves
the worst consistency residual **bit-identical** at 8.88 × 10⁻¹⁵, because Newton
convergence is quadratic and the corrected iterate is already past anything the
suite resolves — the same argument, and the same measurement, as the polar
decomposition's tolerance in §6. The geometric midpoint in the Swift fit is
conditioning rather than correctness: 200 arithmetic bisections reach the same
root. And the pre-loop check that skips the enhanced modes for an
already-degraded element is output-identical to letting the in-step retry catch
it — it saves a wasted Newton pass per step for the rest of a torn element's life,
which is a cost difference and not a behaviour one.

### A flaky test, caught by the gate and not by any run of the suite

The cost test asserted that the elastoplastic element path costs more than the
elastic one. True by 35% on an idle machine, and 60 sequential runs never saw it
fail — but **25 of 48 runs fail with sixteen copies of the suite running at once**,
because under contention the two measurements converge and the inequality becomes
noise. It cost `verify.sh sanitize` one of its six repeat runs, which is the only
instrument that would have found it; it had also, earlier, reported one mutant
killed that was not.

`test_solid_shell.cpp` had already written the rule down in its own cost test —
*a timing assertion tight enough to be interesting is a flaky test on a shared
machine* — and this is what ignoring it looks like. Replaced with bounds two orders
loose, which fire only if something has gone structurally wrong. **The numbers are
for printing, not for asserting.**

That the failure needed *contention* to appear is the reusable part: a suite with
no threads, no GPU and no file I/O in it still had a nondeterminism, and the
operating point that exposed it was load. It is the same lesson as the auto-grain
clamp that only ThreadSanitizer's slowness reached.

### Two test-harness defects, both of the same shape

Both were bisections walking off a function that stops being monotone at failure,
and both produced a confidently wrong physical number rather than a crash:

- The helper that drives a point in **uniaxial stress** bisects the lateral strain
  to null the transverse stress. A torn probe returns *zero* transverse stress,
  which reads as "not enough lateral stretch", so the bracket walked to its far
  end and produced a near-hydrostatic state at η = 15.9 — where the failure strain
  collapses to 4 × 10⁻⁶. A bar reported tearing at a seventh of its failure strain.
  Fixed by searching with failure switched off, which is legitimate because damage
  does not touch the stress until the point tears.
- The load-control loop above it had the same fault: past the failure strain the
  axial stress drops to zero, the "still too low" branch keeps pushing, and the
  bracket runs away. A bar reported tearing at **93% strain**. Fixed by counting
  failure as overshoot.

### And one where the test used a mesh the element is known to be wrong on

The plastic patch test — uniform deformation gradient, so the enhanced parameters
must stay at exactly zero and every Gauss point must carry the point law's stress —
failed by 1.6% and looked like a constitutive bug. It was §6 limit 1: the element
was distorted by moving one face's node alone, making it *warped*, and the warped
patch test fails in proportion to the warp. Distorted prismatically instead — the
same offset on a node and the one above it — it is exact.

---

## 8. A GPU element solver for the solid-shell — **built, and it does not pay**

`engine/gpu/zone_gpu.{hpp,cpp}`, `engine/gpu/shaders/solidshell_forces.comp` and
`solidshell_integrate.comp`, driven by `tools/zone_gpu_probe`. The Phase 3 item
`06-roadmap.md` named, following the tet back-end as instructed.

The kernel exists, it is correct on the closed forms available to it, and **it is
between 1.5× and 4× slower than the 24-thread CPU at every zone size that is
affordable to solve.** The float answer also stops tracking the double one long
before a run finishes. Both are measurements and both are reported here rather
than being worked around, because the second one decides whether this element
belongs on a GPU at all.

The most useful thing that came out of the work is none of that: it is what the
profiling found on the way in.

### The profile first, and it changed the job

**Amdahl's law decides an accelerator, and the kernel's peak throughput does
not.** So before any Vulkan, `zone::Solver` was instrumented per phase and run on
the ferry's own side patch — 192 elements, 450 nodes, a 2 m punch driven 0.048 m
in 6 608 steps, `tools/zone_probe --radius=2.5`:

| phase | one worker | 23 workers |
|---|---|---|
| element evaluation | **98.5%** | 93.9% |
| CSR nodal gather | 0.3% | 1.4% |
| integration | 0.2% | 0.8% |
| energy accounting | 1.0% | 3.9% |

So the tet pattern does apply: this is an element-dominated explicit loop and not
a direct solve. **But half of that 98.5% was work that did not have to happen at
all.**

`elementPlasticUpdate` began by calling `computeForms(rest, …)`, rebuilding the
element's strain-displacement matrices, its enhanced-strain interpolation, its
Gauss weights and its rest Jacobian — every step, for every element, from the
*rest* configuration, which an explicit solve never moves. Hoisting them into a
`solidshell::RestForms` formed once at promotion:

| | uncached | cached | |
|---|---|---|---|
| 192 elements, 1 worker | 5.48 s | **2.73 s** | 2.01× |
| 192 elements, 23 workers | 1.48 s | **0.90 s** | 1.64× |
| 17 800 elements, 23 workers | 21.12 s | **13.03 s** | 1.62× |

**The two answers are bit-identical** — same arithmetic, same order, on the same
numbers — asserted on every reported quantity and on every node position rather
than compared to a tolerance, in `tests/test_solid_shell.cpp` and
`tests/test_zone.cpp`.

Two things are worth carrying forward.

**The tet has always had this.** `fem.cpp` uploads `TetMesh::restInverse` and
`restVolume` and `tet_forces.comp` reads them out of a buffer. The solid-shell
simply never grew the equivalent, and nobody noticed because the per-element cost
*was* measured — 7.3 µs, recorded in `02-simulation.md` §3, and correct. What was
never asked is which part of that 7.3 µs depended on the state being advanced. **A
cost model built from a correct total can still point at the wrong optimisation.**

**The obvious memory objection does not bite, and that was measured rather than
argued.** The cache is 12.0 kB per element against the 4.6 kB per-element
stiffness whose L3 cliff `02-simulation.md` §3 puts at ~6 500 elements, so the
expectation was a crossover at ~2 500. There is none: the cache still wins 1.62×
at 17 800 elements and 214 MB, four times past this machine's last-level cache,
because the arithmetic it removes is 2.2 µs per element per step while streaming
12 kB at the saturated 29 GB/s that section measures is 0.41 µs. Five to one in
the cache's favour even entirely from DRAM.

> **A correction to how this was first measured, because the mistake is the
> reusable part.** Two instruments — a standalone micro-benchmark of
> `computeForms`, and a `steady_clock` pair wrapped round the call *in situ* —
> both said it was 4.2 µs of a 4.3 µs element update, i.e. **97%**, and the
> in-situ pair contradicted itself on arithmetic (its two sub-timers summed to
> more than the phase containing them). The A/B on the real run says 51%. The
> micro-benchmark measured a loop-invariant call on a hot cache and the in-situ
> timer measured a call the compiler could no longer inline. **Only the A/B on the
> real run was right**, which is `CLAUDE.md`'s "timing the real tick instead of
> extrapolating" arriving a second time. The 97% figure was believed for about an
> hour and would have gone into this document.

### What the kernel does, and what lives on the device

One invocation per element, forces to a per-element slot, a companion node kernel
gathering them over the CSR adjacency in a fixed order — the tet's structure
exactly, and for the same reason: no float atomics, and a summation order fixed by
the mesh rather than by warp scheduling.

Per element the device holds 1 505 floats of `RestForms`, 24 force slots and 129
floats of state: eight integration points of plastic strain, back stress,
accumulated plastic strain and damage, **plus the seven enhanced parameters**. The
EAS variables are per-element internal state and they are condensed inside the
element every step, by a Newton on `∫GᵀσdV = 0` with a 7×7 Cholesky in the shader
— the CPU's algorithm, in float. They are read as a warm start and written back,
and never cross the bus. Plasticity is **in scope** and is the reason the kernel is
the size it is; it is also the only version worth building, since the elastic
element is 0.27 µs against the elastoplastic 7.3.

The RestForms hoist is what made the kernel tractable at all: the ANS/EAS geometry
pipeline never has to exist in GLSL, which is the part of this element hardest to
verify.

### Throughput, end to end on the real patch

Not element-updates per second. Wall time for a whole run of the same patch, the
same steps, the same punch — CPU on 23 threads in double, GPU on a 1070 Ti in
float:

| elements | steps | CPU wall | GPU wall | GPU / CPU | GPU µs per element-step |
|---|---|---|---|---|---|
| 192 | 5 505 | 0.70 s | 2.79 s | **0.25×** | 2.6 |
| 768 | 5 505 | 1.80 s | 4.35 s | **0.41×** | 1.03 |
| 3 072 | 5 505 | 6.19 s | 9.11 s | **0.68×** | 0.54 |
| 8 192 | 1 500 | 6.63 s | 27.53 s | **0.24×** | 2.24 |
| 16 384 | 1 000 | 9.29 s | 40.30 s | **0.23×** | 2.46 |

It improves with occupancy to 3 072 elements and then **gets worse again**, which
is the shape that identifies the cause. One thread per element needs about five
hundred floats of thread-private, dynamically indexed state — two copies of eight
integration points of history, a 6×6 algorithmic tangent, a 7×7 Kaa and its
factor, and four 24-vectors. Pascal has 255 registers per thread, so all of it
spills to local memory, which is global memory with an L2 in front; past ~3 000
elements the spill working set leaves the 2 MB L2 and every access is a DRAM
round trip.

> **The tet's mapping does not carry over, and that is the finding.** A linear tet
> is twelve degrees of freedom and no history: its whole state fits in registers,
> which is why `tet_forces.comp` reaches 450–670 M element-updates/s. A solid-shell
> with EAS and eight points of plastic history is two orders of magnitude more
> live state. "One invocation per element" is not the neutral choice it looks
> like — it is the choice that decides whether the kernel runs out of registers.

The design to try next is **one workgroup per element**: eight or thirty-two
threads cooperating over the Gauss points, with the history, the tangent and Kaa
in shared memory rather than in registers. That is a different kernel, not a
tuning pass, and it should be costed before it is written — which is what this
section is for.

### Precision — float is not enough, and where it fails is specific

Two defects and one unresolved gap, in the order they were found.

**1. Absolute ship coordinates in float are fatal, and the fix is free.** The
element works on `u = Rᵀx − X`, a difference of node positions. A side patch on
this ferry sits at (0, −9.9, 8), so |x| ≈ 13 m, where float's 24-bit mantissa
resolves 1 × 10⁻⁶ m — and the displacements of the first thousand steps are
10⁻⁵ m. The difference is then computed from two numbers agreeing in six of their
seven digits. Measured: node positions 4.3 × 10⁻⁷ m apart after a **single** step.
Solving about the patch centroid puts every coordinate inside the zone radius and
costs nothing, because the shift is common to `rest` and `position`; it improved
the work agreement at 1 000 steps by **190×** (5.9 × 10⁻² to 3.1 × 10⁻⁴).

**2. The enhanced-strain condensation has no correct digits in float.** Kaa is
`∫GᵀCG dV`, and the enhanced modes are not commensurate: the three thickness modes
carry `E_ζζ` so their columns of G are scaled by the Voigt transform's `1/t²` —
2.8 × 10⁴ for 12 mm plating — while the in-plane modes are scaled by `1/h²`, of
order ten. Kaa inherits the square of the ratio, so **κ(Kaa) ≈ (h/t)⁴ ≈ 10⁷** on a
perfectly ordinary 0.6 m × 12 mm plate element. Float has seven digits; the
symptom is not noise but an alpha with none of them. Measured at 100 steps, where
the two node fields still agreed to 1.4 × 10⁻⁶ m and no point had yielded, the
GPU's enhanced parameters were **exactly zero** — the residual gate fired on the
first iteration and left them at their warm start — while the CPU's were
3.9 × 10⁻⁸.

Jacobi equilibration of Kaa in the shader (solve for `Dα` with `D = √diag Kaa`,
exact in exact arithmetic, and within a factor of seven of the best diagonal
scaling by van der Sluis) recovers the peak: 3.935 × 10⁻⁸ against the CPU's
3.911 × 10⁻⁸, 0.6%. It does **not** recover every element — the worst element's
alpha is still out by more than its own magnitude.

**3. The run-scale divergence is unresolved, and the negative control is what
says so.** Over a full 5 505-step run the float kernel ends 28% high in plastic
dissipation and **tears 60 elements where the double reference tears none**:

| after 5 505 steps | CPU double | GPU float | relative |
|---|---|---|---|
| work in | 1.1579 MJ | 1.1610 MJ | 2.7 × 10⁻³ |
| plastic dissipation | 1.1114 MJ | 1.4214 MJ | **2.8 × 10⁻¹** |
| torn elements | 0 | **60** | — |
| worst node position | — | 1.88 × 10⁻² m | 39% of the 0.048 m travelled |
| peak equivalent plastic strain | 0.0625 | out by 0.0329 | 53% |

The obvious reading is §2's: an explicit scheme at the CFL limit amplifies float
rounding in the modes at its stability boundary, so of course the two part
company. **That reading is wrong here, and the control that shows it is cheap.**
Perturbing the *double* solver's mesh by 2 × 10⁻⁷ m — the size of the float
representation error — and running the identical 5 505 steps moves the dissipation
by 9.6 × 10⁻⁷ relative, the peak plastic strain by 3.2 × 10⁻⁵, and tears nothing.
So this problem does **not** amplify a perturbation of that size by anything like
five orders of magnitude, and the divergence is not a seeded chaotic one. It is
per-step float error, injected 5 505 times, and it is **biased** — the float run
always dissipates *more*, because noise acts as an imperfection that seeds strain
localisation earlier.

Whether what remains is float being genuinely insufficient or a residual defect in
the kernel is **not settled**, and this section says so rather than picking the
flattering answer. What is settled is that the kernel is not usable for the
question a zone exists to answer: its output is *which panels tore*, and it tears
sixty where the reference tears none.

### What mutation testing found

Eighteen mutants, each a single plausible edit, across the cached element forms,
the zone's use of them, the host's buffer layout and the two shaders. The first
pass killed six; **every one of the survivors was a hole in the suite except two**,
and one of them was a defect sitting in the tree.

- **A defect, and it was live.** `zone::Solver`'s constructor could build the
  cached forms from `position_` instead of `rest_` and pass the entire suite.
  Without a `Preload` the two are the same array, so the mistake is invisible; with
  one they differ by the pre-strain and every element gets the wrong B. The mutant
  was left applied in the tree by a harness accident for the better part of an hour
  and nothing noticed. The bit-identity test now runs a second time **with a
  pre-load on**, which is the only condition that separates the two arrays, and
  asserts the pre-load moved the answer so the case cannot become the case above it.
- **The device suite asked only whether the plate moved and yielded**, so five
  shader mutants survived it: a wrong `sqrt(2/3)`, a lost engineering-shear factor
  of two in the yield norm, a transposed rest Jacobian, a host-side state offset,
  and skipping the Kaa equilibration. A short GPU-against-CPU comparison on the
  integral quantities kills four of them.
- **The transposed rest Jacobian needed a sheared mesh.** On the axis-aligned
  rectangular plate the fixture used, the rest Jacobian is diagonal and a transpose
  is not a change — the same trap `CLAUDE.md` records for the mass lumping and the
  pressure load. The fixture is sheared now, in y only so the element stays
  prismatic through its thickness.
- **`elementRotation` off refused forms** returned a matrix of zeros rather than the
  identity, and nothing asked.
- **`Solver::adopt` shipped with no test at all**, on the caller's own path, in a
  change whose headline feature was well covered. Writing one found that it
  *advanced* the state it was adopting: recovering the stress runs
  `elementPlasticUpdate`, which commits, so a point sitting just under its damage
  limit would be tipped over it by the act of being read. It runs on a copy now.

The suite went from 6 002 to 6 062 checks and from 6 to 15 kills of 18.

**Two survivors, and both are findings rather than gaps.**

`kRoot23` can be wrong by 0.4% — a 0.8% error in the yield radius — and no test
sees it, because the float path's own agreement with the double reference is
5.8 × 10⁻² on the plastic dissipation after four hundred steps of a barely-yielded
plate. **The instrument is less precise than the error it would have to catch.**
That is the same statement as the precision section above, arrived at from the
other direction, and it is the clearest single argument that this kernel cannot be
validated to the standard the CPU element was.

Keeping the enhanced modes on a degraded element is *behaviourally* equivalent
here, though not provably so: the equilibration's positive-diagonal check refuses
to factor a rank-deficient Kaa and leaves alpha at the warm start, which the step
of the first failure already zeroed. It is argued equivalent rather than shown to
be, and the argument depends on a guard added for an unrelated reason — which is
worth knowing if that guard is ever removed.

Two deliberate controls behaved as controls should: dropping the *cached* call in
the elastic path in favour of the uncached one survives, which is the bit-identity
claim being true rather than a hole; and dropping the equilibration is killed,
which says the fix earns its place.

### What to do about it

1. **Keep the `RestForms` hoist.** It is 2× on the CPU, free, and bit-identical.
   It is the one part of this work that is unambiguously worth having.
2. **Fix the conditioning in the element, not in the shader.** The enhanced modes
   are a basis and their scaling is a free choice; normalising each column of G by
   its own natural scale makes Kaa O(1) and is exact, so it costs the double CPU
   path nothing and would remove the largest single obstacle to a float one. That
   is a change to `solid_shell.cpp` and it should be made and validated against the
   existing 207 assertions *before* any more GPU work.

   **Two things decide whether it is safe, and neither is the scaling itself.**
   `Kua` has to be scaled with `Kaa` or the condensation is no longer the same
   operator — scaling column *j* of G by *s* and α*ⱼ* by 1/*s* leaves `G α`
   unchanged, which is why it is exact, but only if every place α appears agrees on
   *s*. And **α is persistent per-element state**: `ElementPlasticState` carries the
   enhanced parameters between steps, so a scale that is recomputed differently on
   any step silently reinterprets the history already stored.

   Both fall out if the scale is taken from the **rest** configuration and cached —
   which is exactly where `RestForms` now lives. The `(h/t)⁴` conditioning is a
   property of the element's shape, not of its current deformation, so a rest-based
   scale is constant for the element's life, consistent for the stored α by
   construction, and free. Deriving it from the *current* configuration would be
   the version that looks equivalent and is not.
3. **If a float kernel is still wanted after that, keep alpha in double.** The EAS
   block is 7×7. Even at Pascal's 1/32 fp64 rate it is a small share of the kernel,
   and it is the only part that has been shown to need the digits.
4. **Re-map the kernel to a workgroup per element** before measuring throughput
   again. The current numbers measure register spilling, not the element.

Until 2 and 4 are done, **the CPU is the faster and the more trustworthy path for
Tier 2**, and the honest statement of this item's status is that the pattern
applies, the kernel was built to it, and the element does not fit the mapping.
