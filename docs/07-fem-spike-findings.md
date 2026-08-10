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

> **The paragraph above is the only thing in §1 or §2 that `fem_spike` cannot be
> asked to produce.** Every other cell of both tables comes back on
> `./build/fem_spike` exactly as printed, to every digit. The peak-velocity
> comparison and the 198–202 step sweep were scaffolding and were never checked in,
> so "over 200%" and "~4 steps" rest on the authority of the session that ran them
> and cannot be re-derived from this tree.

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

> **The two count columns re-derive; the three clock columns never have, and an
> attempt was refused by the box.** `./build/fem_spike` returns 960/315,
> 11 520/2 745, 46 080/10 285, 230 400/45 225 and 491 520/95 337 exactly, so the
> mesh sizes §4 works forward from are sound. The timings were re-run with an
> unrelated desktop application holding the device at 74–100% for the whole session:
> the 960-tet row came back within 1% (0.0082 ms/step, 116.9 M) and the two largest
> read **4.5× slow**. That is the signature of a busy GPU rather than of a moved
> figure — a competing consumer costs a long kernel far more than a launch-bound
> one, which is the mechanism §8 measures directly — so **nothing in this table is
> treated as re-measured**, and nothing here is gated: `scripts/check-figures.sh`
> reaches §8 and no other section of this file. The single-core 4.5–6.7 M range is
> worse off still: `fem_spike` benchmarks one mesh size on the CPU, so even an idle
> box would return one point of a range that has two.

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
  > the element is in plane**. Measured: `dt·c_p/t` runs 0.992 at 5t to 0.9999 at
  > 50t, i.e. flat. The step is *not* set by the in-plane size. The 5–10×
  > figure survives by accident, because it is the right comparison against a tet
  > mesh with 8 elements through the thickness (measured 6.7×), but anyone sizing
  > elements from the stated mechanism would over-predict the affordable step by
  > the in-plane aspect ratio, which is 2.5× at 50 mm elements on 20 mm plate. A
  > timestep governed by in-plane size needs a classical shell element with no
  > thickness stretch, or selective mass scaling on the thickness mode.
  >
  > > **Both numbers in the sentence above were wrong, and one by a factor of ten.**
  > > It read `dt·c_p/t = 0.999` "from 5t to 50t" and an aspect ratio of "25× at
  > > 50 mm elements on 20 mm plate". 50 mm on 20 mm plate is an aspect ratio of
  > > **2.5**, not 25, and the over-prediction is the same 2.5: the stability sweep
  > > in `tests/test_solid_shell.cpp` reports `dt·c_p/h = 0.3867` at t = 20 mm,
  > > h = 50 mm, so believing the in-plane size sets the step over-predicts it by
  > > 1/0.3867 = **2.6×**. The flatness figure was quoted from the wide end of its
  > > own sweep: `dt·c_p/t` is 0.9924 at 5t and only reaches 0.999 at 10t. Re-derived
  > > by running `./build/shipsim_tests` and reading the *explicit stability limit*
  > > table under `--- solid-shell elements ---`, which is where the 6.7× beside it
  > > comes from and always did.
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

All five rows come off the *measured cost, one core* block `./build/shipsim_tests`
prints under `--- solid-shell elements ---`, and each is the **minimum over five
passes** rather than one reading, because contention can only add time:

| | |
|---|---|
| Element stiffness formation | **21.9 µs** — once, when a zone is promoted |
| The same for a plain hex | 13.5 µs, so the assumed strains cost **1.62×** |
| Internal force per step | **293 ns** |
| `fem.cpp` linear tet internal force per step | 129 ns, so one solid-shell = **2.3 tets** |
| Explicit step, 20 mm plate | 3.25 µs (`t / c_p`, in-plane size irrelevant) |

> **Two of those had drifted, and the tet row is what says so rather than the box
> being busy.** The table read 21.1 µs, 13.4 µs, **267 ns**, 129 ns and 1.57× /
> 2.1 tets. The solid-shell's own two figures are 4% and 10% above what was
> published; the plain hex is within 1% and **the tet is exactly 129 ns**. All four
> are timed in the same loop, in the same run, on the same box, so a load that
> inflated the solid-shell would have inflated the tet beside it — and across five
> passes the tet's minimum lands on 129 while the internal force never once came
> near 267 (293, 295, 303, 307, 309). That is a **control**, not an argument: the
> re-measurement was taken while the box was 22–38% busy with sibling agents, and
> the tet row is the only reason the two moved rows can be read as the element
> changing rather than the machine. *Why* the element got 10% slower was **not
> chased and is not explained here**: `internalForce` is a polar decomposition and
> a 24 × 24 matvec on a stiffness formed outside the timed loop, so it touches
> neither the enhanced-mode normalisation nor the constraint path that arrived
> after these figures were taken.

Per element the solid-shell is twice a tet. That is not the comparison that
matters. Per **square metre of 20 mm plating per simulated second**, single core:

| | elements/m² | steps/s | cost |
|---|---|---|---|
| Linear tets at 2.5 mm (8 through the thickness, the spike's 11%-error mesh) | 7.7 M | 4.1 M | **4.1 × 10⁶ s** |
| Solid-shell at 50 mm, one through the thickness | 400 | 3.4 × 10⁵ | **40 s** |

**A factor of 1.0 × 10⁵.** A 20 m × 10 m collision zone is 200 m² of shell, so one
simulated second costs **~8 × 10⁸ core-seconds as tets and ~8 000 core-seconds as
solid-shells** — about five and a half minutes of wall time per simulated second on
the 24-thread CPU, or well under a minute on the GPU at the throughput §3 measured.
That is inside the time dilation the engine already plans for, and it is what
makes the 20 m damage zone affordable rather than theoretical.

> **Correction.** That sentence previously read "~2 core-hours as tets and
> ~2 core-minutes as solid-shells", and **both figures were wrong** — 200 × 40 s is
> 8 000 core-seconds, which is two core-*hours*, and 200 × 4.1 × 10⁶ s is 26
> core-*years*, not two core-hours. The wall-time figure beside them was right all
> along (8 000 / 24 = 5.6 minutes), which is why the error survived: the number
> anyone would sanity-check was correct and the two feeding it were not. Found while
> costing the elastoplastic path against it.
>
> > The tet column and the two derived from it are unchanged on re-measurement
> > (4.07 × 10⁶ s, minimum of five). The solid-shell column moved with the 267 → 293 ns
> > above it: 37 → 40 s, 1.1 × 10⁵ → 1.0 × 10⁵, 7 400 → 8 000 core-seconds. The
> > `steps/s` column was independently wrong: the test forms it from
> > `criticalTimestep(..., 0.9)`, i.e. 1/(0.9 × 3.25 µs) = **3.4 × 10⁵**, and 3.3 × 10⁵
> > never reproduced the 37 s printed beside it.

Stiffness formation at 22 µs means promoting a 10⁵-element zone costs ~2 s on one
core, ~0.1 s threaded. That is a hitch at promotion, not a per-frame cost, and it
is the number to watch if zones are promoted during play rather than at impact.

### What it was checked against

Closed forms throughout, not eyeballed output. **428 assertions** — counted by
running `runSolidShellTests()` on its own against `testing::checkCount()`. It was
207 when the mutation pass below finished and has roughly doubled since; the number
here had never been re-counted:

- **The patch test**, exactly, on a distorted patch: displacement to 1.0 × 10⁻¹⁵
  and constant stress to 1.7 × 10⁻¹³. Distortion in-plane is arbitrary.
- **Rigid body motion** carries no force, including finite rotations to 3 rad; and
  the stronger statement, **frame indifference** — a 25% stretch rotated 1.7 rad
  gives the rotated force to 4.7 × 10⁻¹³ of it.
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

## 8. A GPU element solver for the solid-shell — **re-mapped: it now wins on throughput and still loses on precision**

`engine/gpu/zone_gpu.{hpp,cpp}`, `engine/gpu/shaders/solidshell_forces.comp`,
`solidshell_forces_wg.comp` and `solidshell_integrate.comp`, driven by
`tools/zone_gpu_probe`. The Phase 3 item `06-roadmap.md` named, following the tet
back-end as instructed.

> **This section was rewritten after a re-measurement, and three of its claims did
> not survive.** The first version reported the kernel as 0.23–0.68× the CPU and
> "not paying", diagnosed register spilling as the cause, prescribed one workgroup
> per element as the fix, and reported that float tears 60 elements where double
> tears none. The re-measurement says: **the diagnosis was right and the fix works
> — 1.26–2.43× against the CPU, which is 2.1× to 10.7× off the kernel time.** The
> tearing figure
> **does not reproduce under any configuration** and appears never to have come
> from the run the rest of that table came from. And the conditioning fix that item
> 2 below prescribed, which was expected to be what rescued float, **changes
> nothing measurable for this kernel** — for a reason worth knowing. Each is
> separated by measurement below.
>
> Three documents carried the 60 figure — this one, `06-roadmap.md` and
> `02-simulation.md` — each quoting the previous rather than the tool. That is
> `CLAUDE.md`'s "nothing at all tests a comment", for the second time.

> **And the last item on this section's own list has since been closed, as a negative.**
> "Keep alpha in double" was item 3 of *What to do about it* and the only thing left to
> try. It is built — five kernels from one source, `--eas=float|tight|solve|condense|
> newton` — and it does not close the gap: the five land between 40 and 44 torn elements
> at 768 where the reference is 32 and between 204 and 247 at 3 072 where it is 162, the
> spread across the precisions is about half the gap they were meant to close, and its
> sign reverses between the two sizes. What the five *do* sort by is the stopping rule
> they carry, not their arithmetic. It costs 5–15× on the kernel, which undoes the
> remap. *Alpha in double — measured, and it is not what was missing*, below, has the
> tables, the mechanism, and the two claims in item 3 that turned out to be wrong.

> **The throughput numbers here have since been re-measured on a box checked idle,
> and they hold.** 1.26–2.43× re-measures as **1.27–2.44×** over six interleaved
> repeats, and 0.23–0.68× as 0.22–0.64×; every cell is within 6% of what it
> replaces. The range is quoted below as 1.26–2.43× where other documents quote it,
> because a re-measurement agreeing to 1% is a confirmation and not a new figure.
> What that pass *did* overturn is a correction: five cells this section had marked
> as "not reproducing" reproduce exactly, and the arrows are withdrawn — see *The
> correction to the correction*. And the fp64 cost table, which was taken against an
> unrelated GPU consumer and said so, has been re-taken idle; it moved.

The kernel exists, it is correct on the closed forms available to it, and after the
remap it is **1.26× to 2.43× the 23-worker CPU** end to end on the real patch. The
float answer still stops tracking the double one long before a run finishes, and
that — not throughput — is what stands between this kernel and use. It is now known
*not* to be the enhanced block's precision.

The most useful thing that came out of the original work is neither: it is what the
profiling found on the way in.

### The profile first, and it changed the job

**Amdahl's law decides an accelerator, and the kernel's peak throughput does
not.** So before any Vulkan, `zone::Solver` was instrumented per phase and run on
the ferry's own side patch — 192 elements, 450 nodes, a 2 m punch driven 0.060 m
in 6 608 steps, `tools/zone_probe --radius=2.5 --depth=0.06 --forms-cache=never`:

| phase | one worker | 23 workers |
|---|---|---|
| element evaluation | **98.5%** | 93.9% |
| CSR nodal gather | 0.3% | 1.4% |
| integration | 0.2% | 0.8% |
| energy accounting | 1.0% | 3.9% |

> **The command above used to read `tools/zone_probe --radius=2.5`, and that is a
> different experiment.** `--depth` defaults to **0.45 m**, which is 42 396 steps
> and half the patch torn, and it profiles 92.1 / 1.8 / 1.2 / 4.9 on 23 workers —
> neither column of this table. Two parameters were missing beside the number:
> the depth, and `--forms-cache=never`. The second matters because **this is the
> pre-hoist profile** — it is the thing the next paragraph is about — and the
> switch that turns the hoist off did not exist when it was taken. Supplied both,
> the one-worker column comes back exactly: **98.5–98.6 / 0.3 / 0.2 / 0.9** on
> three passes, and the 23-worker column comes back 93.3–93.8 / 1.1 / 0.8 /
> 4.3–4.7 on a box 31–36% busy, where the serial energy phase is the one a
> competing load inflates. The travel was also wrong: 0.048 m is the last *history
> sample* the run prints, and the punch is driven 0.060 m. This is the same fault
> as *The correction to the correction* below — a parameter of the experiment not
> written next to the number — in the section that records it.

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

> **What re-measures here is the ratio, and only the one-worker one.** On the
> commands above with `--forms-cache=never` against the default, the one-worker
> pair is 6.57 s / 3.29 s = **2.00×** against the published 2.01× — the absolute
> times are 20% high on a box 27–36% busy, and the ratio is untouched by that
> because both halves are one thread. The 23-worker pair comes out 2.43 s / 1.36 s
> = 1.79× against 1.64×, which is a 23-way split on a box with a third of its
> cores in other hands and is not evidence of anything. **The third row cannot be
> re-run at all**: nothing records its step count, no `--radius`/`--sub` is given,
> and the configuration that produces exactly 17 800 elements had to be searched
> for — it is `zone_probe --force --radius=11 --sub=10`. At the depth the other two
> rows use that is 7 001 steps and 106 s / 154 s, and the published 13.03 s implies
> a run roughly eight times shorter. The row is retained as taken and flagged as
> unreproducible, not corrected, because there is nothing to correct it *to*.

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

### The remap, and what it did

`solidshell_forces_wg.comp` is the same element under **one workgroup of 32 threads
per element** instead of one thread per element. `gpu::Mapping` selects between
them, both are built, both are tested, and `zone_gpu_probe --mapping=` runs either,
so the comparison below is one command apart rather than a rebuild apart.

The split is not uniform, because the element's work is not:

| phase | threads | share of a Newton iteration |
|---|---|---|
| load, `u = Rᵀx − X`, final `f = −R f_int` | 24, one per DOF | small |
| strain at the Gauss points | 48 tasks over 32 threads | ~10% |
| **return map** | **8, one per Gauss point** | ~20% |
| **Kaa = ∫GᵀCG dV** | 49 tasks over 32 threads | **~70%** |
| equilibrate, 7×7 Cholesky, solve | thread 0 alone | ~1% |

The one that mattered is the return map: a point's plastic history now lives in
**its own thread's** registers rather than all eight copies living in every
thread's. Kaa is where the arithmetic is (8 × 49 × 36 multiply-adds against the
return map's ~4 000) and it parallelises to 49 independent entries, so it was the
phase that had to spread wide. The 6×6 tangents cross from the return-map thread to
the Kaa threads through shared memory, which is the whole trade.

**The per-element force slots and the fixed-order CSR gather are untouched.** The
same 24 floats land in the same slot in the same order and `solidshell_integrate.comp`
is unchanged, so the two mappings are interchangeable behind an identical buffer —
which is also what makes it possible to run them against each other.

### The compiler's own account, which had not been asked

§8's first version diagnosed register spilling from the **shape of the throughput
curve**. `VK_KHR_pipeline_executable_properties` reports it directly, and
`zone_gpu_probe --stats` now does. `node_integrate.comp` is queried first as a
calibration, because this driver reports `Local Memory Size` with a large fixed
offset and the raw number means nothing on its own:

| shader | registers/thread | spill over calibration |
|---|---|---|
| `node_integrate.comp` (calibration) | 32 | 0 B |
| `tet_forces.comp` † | 40 | **0 B** |
| `solidshell_forces.comp` — one invocation | 128 | **1 936 B = 484 floats** |
| `solidshell_forces_wg.comp` — one workgroup | 64 | **96 B = 24 floats** |

Re-run and identical, every row: these come from the driver's compiler and not from
a clock, so they are the one part of this section that a busy box cannot move.

> † **except the tet row, which `--stats` does not print.** It queries the two
> solid-shell mappings, the calibration and the four fp64 variants;
> `tet_forces.comp` is not among them, so that cell cannot be re-derived from the
> tool that produced the rest of the table and is carried on the authority of the
> run that first took it. Flagged rather than fixed, because the tet is not what
> this section is about — but it is exactly one cell of one table with a different
> provenance from its neighbours, which is the thing this section has now been
> caught by twice.

**484 floats.** The first version guessed "about five hundred floats of
thread-private state" from counting the arrays by hand and was right to three
figures, and the tet's zero is that section's claim about the tet checked rather
than repeated. The remap takes the spill down **20×** and halves the register
count; derived from the register count and Pascal's 65 536-register file, warp
occupancy goes from 16 of 64 warps per SM to 32.

### Throughput, end to end on the real patch

Not element-updates per second. Wall time for a whole run of the same patch, the
same steps, the same punch — CPU on 23 threads in double, GPU on a 1070 Ti in
float. Both mappings, same binary, one sweep, **six interleaved repeats of every
cell**, `tools/zone_gpu_probe/sweep.sh`:

| elements | steps | CPU wall | invocation | vs CPU | workgroup | vs CPU | workgroup ÷ invocation |
|---|---|---|---|---|---|---|---|
| 192 | 5 505 | 0.67 s (0.67–0.70) | 2.71 s (2.71–2.83) | **0.25×** | 0.43 s (0.43–0.47) | **1.56×** | 6.3× |
| 768 | 5 505 | 1.62 s (1.62–2.76) | 4.10 s (4.10–4.42) | **0.40×** | 1.28 s (1.28–1.30) | **1.27×** | 3.2× |
| 3 072 | 5 505 | 5.63 s (5.63–10.28) | 8.83 s (8.83–9.09) | **0.64×** | 4.19 s (4.19–4.39) | **1.34×** | 2.1× |
| 8 192 | 1 500 | 6.17 s (6.17–7.12) | 26.02 s (26.02–28.19) | **0.24×** | 2.77 s (2.77–2.88) | **2.23×** | 9.4× |
| 16 384 | 1 000 | 8.84 s (8.84–9.55) | 40.98 s (40.98–42.28) | **0.22×** | 3.62 s (3.62–6.45) | **2.44×** | 11.3× |

Each cell is the **minimum over its repeats**, with the full range beside it.
Minimum rather than mean, because contention can only *add* time: an inflated run
is a run with someone else's work in it, and averaging it in publishes a number
that is partly the neighbour's. The ratios are formed from the two minima. The
last column is how much faster the remap is than the mapping it replaced, which
is the quantity the prose calls "off the kernel time".

**The published range survives the re-measurement: 1.26–2.43× against 1.27–2.44×
here**, and the invocation mapping's 0.23–0.68× against 0.22–0.64×. Every number
in the previous version of this table is within 6% of its re-measurement and most
are within 2%.

> ### The box was idle, and that was established rather than assumed
>
> The previous version of this table asserted idleness in a callout. The fp64 work
> below then published a table that was *not* taken on an idle box and said so, and
> the two claims are indistinguishable from the outside — which is the reason this
> is now instrumented instead of stated.
>
> **Three instruments, because they fail differently.**
>
> - **Before every run**, in the gap where nothing of ours is executing, the CPU's
>   busy fraction from a `/proc/stat` delta and the GPU's utilisation and clock from
>   the driver. Across all sixty rows the worst reading was **9.7% of the CPU and
>   13% of the GPU**, against a desktop floor of about 2%. Deliberately not the load
>   average: a 1-minute average still carries the sweep's own previous run, so it
>   would flag every row after the first.
> - **During every run**, everything the machine burned minus everything this probe
>   burned, over the run's own wall clock, so our own 23 workers cancel. The before
>   check cannot see a neighbour that starts two seconds later, and one did; this
>   column was added *because* of what the third instrument found in the sweep
>   above, and it catches the same events live in the fp64 sweep below.
> - **The repeats themselves**, which are the sharpest of the three because they
>   measure the thing that matters rather than a proxy for it. **The GPU column
>   repeats to within 2–9% at every size bar one run**, where the fp64 session below
>   — with an unrelated consumer on the device — saw the same configuration move
>   **35%**. A kernel time that repeats to 2% cannot have been taken against a
>   competing load. The exception is named in the next paragraph rather than dropped.
>
> **The GPU was idle. The CPU was not, and the difference matters.** Two sibling
> agents were gating this repository throughout, and their `shipsim_tests` runs land
> on the CPU reference: four of the sixty rows carry a CPU wall 15–83% above the
> minimum for their configuration, and one carries a GPU wall 78% above it, because
> that suite has device tests in it too. Every one of those rows read **under 3%
> busy** on its before-check. That is why the table quotes minima and ranges rather
> than a mean of six, and why the CPU column's range is wide at 768 and 3 072 while
> the GPU column's is not.
>
> The one figure to distrust is unchanged and for the reason previously given: at
> 192 elements six workgroups is a handful of warps, so it is latency-bound and
> anything else on the device shows in it directly.
>
> > **A later audit could not repeat this table, and the way it failed is the
> > callout's own mechanism observed again.** With a desktop application holding the
> > device at 74–100% for an entire session, single runs of the same commands gave
> > 0.79× / 1.32× / 1.02× / 1.02× / 1.22× on the workgroup mapping against the
> > 1.56 / 1.27 / 1.34 / 2.23 / 2.44 above — while the **invocation** mapping at
> > 3 072, the longest kernel in the sweep, came back at **0.62× against its
> > published 0.64×**. The one measurement that survived the contention is the one
> > made of the most work, which is what *the contended session understated the
> > cost* says below, and it is why none of the rest was written into the table.
> > Every **deterministic** cell of §8 — every torn count, dissipation, alpha and
> > register count — was re-derived in that session and holds; nothing on a clock
> > was.

**The degradation past 3 000 elements is gone**, and that is the load-bearing part.
The old curve improved to 3 072 and then fell off a cliff — 0.64× to 0.24× to 0.22×
— which is a spill working set leaving the 2 MB L2. The new one does not: from 768
up it rises with size, 1.27× → 1.34× → 2.23× → 2.44×, which is what a kernel does
when it is bound by occupancy. The diagnosis was right.

**It is not monotone across the whole range, and the 192-element point is the
exception rather than the trend.** 1.56× there sits above the 1.27× at 768, because
192 elements is six workgroups short of filling one SM's worth of the device while
also being the size at which the CPU's own 23-way split has the least to work
with — both sides are inefficient and the ratio says little. The four larger sizes
are the ones the shape should be read off.

> **The tet's mapping does not carry over, and that was the finding.** A linear tet
> is twelve degrees of freedom and no history: its whole state fits in registers —
> measured above at zero spill — which is why `tet_forces.comp` reaches 450–670 M
> element-updates/s. A solid-shell with EAS and eight points of plastic history is
> two orders of magnitude more live state. "One invocation per element" is not the
> neutral choice it looks like: it is the choice that decides whether the kernel
> runs out of registers.

**What the remap did not fix is the precision**, and every quantity below is
unchanged by it to four significant figures. On the unit-scale fixture the two
mappings differ by 2.8 × 10⁻⁵ on the plastic dissipation where either differs from
the double CPU by 6.8 × 10⁻² — three orders of magnitude, which is what says the
remap changed the mapping and not the computation. Speed was never the reason this
kernel could not be used.

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

**3. The run-scale divergence is real, it is 27–34% on the dissipation, and it
moves the torn count by a quarter to a half.** Over a full 5 505-step run:

| after 5 505 steps, 192 elements | CPU double | GPU float | relative |
|---|---|---|---|
| work in | 1.1579 MJ | 1.1610 MJ | 2.7 × 10⁻³ |
| plastic dissipation | 1.1114 MJ | 1.4214 MJ | **2.8 × 10⁻¹** |
| torn elements | 0 | **0** | — |
| peak damage (a point fails at 1.0) | 0.6527 | 0.4094 | 37% |
| worst node position | — | 1.88 × 10⁻² m | 39% of the 0.048 m travelled |
| peak equivalent plastic strain | 0.0625 | out by 0.0329 | 53% |

> ### The correction: "tears 60 elements where the reference tears none" does not reproduce
>
> Every other cell of that table reproduces to four significant figures on today's
> tree. The torn count does not: it is **0 against 0**, and the GPU's peak damage
> is 0.41 against the CPU's 0.65 — the float run is *further* from tearing than the
> double one, by a factor of 1.6, with **zero failed Gauss points on either side**.
> For 60 elements to tear, 480 points would have to reach 1.0.
>
> It is not a configuration difference. Every combination of the three float
> accommodations §8 describes was run, and none of them tears anything:
>
> | enhanced-mode normalisation | Kaa equilibration | centroid shift | GPU torn | GPU peak damage |
> |---|---|---|---|---|
> | on | on | on (as shipped) | 0 | 0.4094 |
> | **off** | on | on | 0 | 0.4095 |
> | on | **off** | on | 0 | 0.4094 |
> | **off** | **off** | on | 0 | 0.4095 |
> | on | on | **off** | 0 | 0.4303 |
> | **off** | on | **off** | 0 | 0.4303 |
>
> A run that tore 60 of 192 elements could not have produced the same dissipation
> to five digits as one that tore none, so the figure cannot have come from the run
> the rest of that table came from. `zone_gpu.cpp` separately records "by 5 500
> steps the GPU had torn 58 elements" for the *pre-centroid-shift* kernel, and 58
> and 60 are close enough that the most likely account is a pre-fix number carried
> forward into a post-fix table. It is recorded here as unreproduced rather than
> explained away, and the probe now prints the peak damage on both sides so that
> "0 torn against 0 torn" is a margin rather than two zeros that might be zero for
> unrelated reasons.

**The conditioning fix changes nothing for this kernel, and the reason is worth
knowing.** `computeRestForms` now normalises the enhanced modes and κ(Kaa) is a
constant 3.50 rather than (h/t)⁴ (item 2 of *What to do about it*, below). The
expectation was that this would be what rescued the float path. Measured by A/B on
the identical 5 505-step run — the normalisation switched off in `computeForms` and
nothing else changed:

| 192 elements, 5 505 steps | normalisation on | normalisation off |
|---|---|---|
| plastic dissipation | 1.4214 MJ | 1.4214 MJ |
| worst node position | 1.878 × 10⁻² m | 1.878 × 10⁻² m |
| peak eps_p error | 3.293 × 10⁻² | 3.294 × 10⁻² |
| torn | 0 | 0 |

**Nothing moves.** The reason is that the shader was *already* equilibrating Kaa
with a Jacobi scaling before factoring it, and Jacobi equilibration and column
normalisation address the same conditioning — so the second fix had nothing left to
do. Turning the shader's equilibration off as well leaves the answer identical
again. The element-level normalisation is still the better fix (it is free, it is
exact, and it fixes every path including the CPU's and any future one that does not
equilibrate), but **it is not what stands between this kernel and float**, and the
brief that predicted it would was measuring the wrong obstacle.

**So where float actually fails is the tearing, and the negative control now says
so at sizes where tearing is live.** 192 elements never tears on either path, which
is why the original run could say nothing about it. At 768 and 3 072 elements the
*double* reference does tear, and that is where the comparison bites:

**All of it at a fixed 5 505 steps**, which is the whole of what was ever wrong with
this table — see the correction below.

| | 768 elements | 3 072 elements |
|---|---|---|
| CPU double, torn | **32** | **162** |
| GPU float, torn (workgroup) | **40** | **247** |
| GPU float, torn (invocation) | **44** | **213** |
| **control: double, mesh jittered 2 × 10⁻⁷ m, torn** | **32** | **162** |
| CPU dissipation | **1.5194** MJ | **1.2910** MJ |
| GPU dissipation, relative | **2.57** × 10⁻¹ | **3.38** × 10⁻¹ |
| control dissipation, relative | **7.69** × 10⁻³ | **5.64** × 10⁻⁴ |

> ### The correction to the correction: this table was right, and the arrows were the error
>
> An earlier revision annotated five of these cells with arrows — 40 → 41, 247 → 248,
> 213 → 241, 1.5194 → 1.5205 MJ, 1.2910 → 1.2963 MJ — on the grounds that the tool no
> longer produced the figures printed, and concluded that "five more cells of this
> table did not come from the run the rest of it did". **Every cell above is what the
> tool produces today**, on six repeats of a single sweep, identical on every one.
> The arrows are withdrawn.
>
> **What produced the other set of numbers was eight extra steps at 768, and forty at
> 3 072.** The probe derives
> its step count from the punch *depth* unless `--steps` is given, and the critical
> timestep falls with element size — 1.8163 µs at 192 elements, 1.8138 at 768, 1.8033
> at 3 072 — so the derived count is 5 505, 5 513 and 5 545. This table was taken at a
> **fixed** 5 505; the re-run that produced the arrows let the count be derived. Run
> today's tree at the derived count and the "corrections" come back exactly, all five:
>
> | 768, derived 5 513 steps | 3 072, derived 5 545 steps |
> |---|---|
> | CPU dissipation **1.5205** MJ, torn 32 / **41** | CPU dissipation **1.2963** MJ, torn 162 / **248** workgroup, **241** invocation |
>
> **So neither table had cells from two different runs.** Each was internally
> consistent; they were two different experiments, and the second was read as a
> re-measurement of the first.
>
> **The reusable part is how convincing the wrong diagnosis was.** That re-run
> established that its figures were not noise (identical on a repeat), not the
> reduction order (identical at 1, 4, 12 and 23 workers), and not the tree — it built
> `776d15e` from a `git archive` and got its own numbers back rather than the
> published ones. All three are true and none of them could see the fault, because
> **every one of those controls varied something other than the command line, and the
> command line was what differed.** `6b928ba`'s own binary, rebuilt here from a `git
> archive` and given `--steps=5505`, produces 1.5194 MJ and 40 torn — the figures its
> commit message says it could not reproduce.
>
> A step count is a parameter of the experiment and it was not recorded next to the
> numbers. It is now, in this table's own heading and in every row of the CSV
> `sweep.sh` writes.

**The control tears exactly the reference's count at both sizes and the float
kernel is 25% and 52% over it.** That is the argument, and it is much stronger than
the one it replaces: a geometric perturbation the size of float's own
representation error, applied to the double solver, does not move the torn count by
a single element, while float moves it by a quarter to a half. At 192 elements the
same control moves the dissipation by 9.6 × 10⁻⁷ and the peak plastic strain by
3.2 × 10⁻⁵ — reproducing the original control exactly.

So the divergence is not a seeded chaotic one. It is per-step float error, injected
thousands of times, and it is **biased**: the float run always dissipates *more*,
because noise acts as an imperfection that seeds strain localisation earlier — and
earlier localisation is precisely what turns into extra torn elements.

What is settled is that the kernel is still not usable for the question a zone
exists to answer. Its output is *which panels tore*, and it tears a quarter to a
half too many. What is **no longer** true is that it is slow.

### Alpha in double — measured, and it is not what was missing

Item 3 of *What to do about it* was the last thing left on this list: "keep alpha in
double — the EAS block is 7×7, and it is the only part that has been shown to need the
digits." It has been built and measured, and **it does not close the gap.** What it
does close is a different question that turned out to be the interesting one.

`solidshell_forces_wg.comp` now compiles into **five kernels from one source**, selected
by `gpu::EasPrecision` and `zone_gpu_probe --eas=`. One source matters: any difference
between two of them is the precision and nothing else, and at level 0 the SPIR-V is
**byte-identical** to the kernel this work started from, so the float column below is
the same kernel §8 has been quoting throughout.

| `--eas=` | what is in fp64 | stopping rule |
|---|---|---|
| `float` | nothing — as shipped | 1e-9·σ_y V, 12 iterations |
| `tight` | nothing — **the control** | the CPU's 1e-16, 40 iterations |
| `solve` | the Jacobi equilibration, the 7×7 Cholesky, both substitutions | shipped |
| `condense` | + Kaa and the residual, over an fp64 copy of G and the weights | shipped |
| `newton` | + alpha, its correction, its persistent per-element state | the CPU's |

`tight` is there because `newton` changes the tolerance and the arithmetic together.
Without a kernel that changes only the tolerance there is no way to say which of the two
moved anything — and it turns out to be entirely the tolerance.

**The fp64 copy of G is not optional at `condense` and above.** A double accumulation of
`∫GᵀCG` over a G that was rounded to float is an operator known to seven digits however
it is summed, so a level that widened the arithmetic without widening the operator would
have measured nothing and reported it as a result.

#### First, the float baseline

Same probe, same five sizes, workgroup mapping, one sweep of six repeats, and the
step count **fixed** at each size rather than derived — 5 505, 5 505, 5 505, 1 500,
1 000, which is what *Throughput* above runs and what *Precision* above records:

| | 192 | 768 | 3 072 | 8 192 | 16 384 |
|---|---|---|---|---|---|
| CPU double, torn | 0 | **32** | **162** | 0 | 0 |
| GPU float, torn | 0 | **40** | **247** | 0 | 0 |
| control: double, jittered 2e-7 m | — | **32** | **162** | — | — |
| GPU dissipation, relative | 2.79e-1 | 2.57e-1 | 3.38e-1 | 9.32e-2 | 2.28e-2 |
| control dissipation, relative | — | 7.69e-3 | 5.64e-4 | — | — |

Every cell here is identical on all six repeats — both paths are deterministic, and
the sweep asserts that by recording each repeat rather than the first. At 192
elements every figure reproduces §8's original to four significant figures:
dissipation 1.1114 against 1.4214 MJ, peak damage 0.6527 against 0.4094, zero torn
on both sides. The reference's 32 and 162 reproduce, and so does the negative
control's exact 32 and 162 — which is the load-bearing part of the argument and it
is intact.

> **An earlier revision of this subsection was titled "and two more cells of the table
> above do not reproduce", and it was wrong.** It compared a derived step count
> against a fixed one; the arrows it added to *Precision* are withdrawn there, with
> the mechanism. Nothing in the conclusion moved then and nothing moves now — float
> is 25% and 52% over a reference that a float-sized geometric jitter reproduces
> exactly — but two revisions of this document have now published a "figure that does
> not reproduce" which did reproduce, under the parameters it was taken at. **The
> first cost a real figure** (the 60 torn elements, correctly withdrawn); the second
> cost five figures that were correct. A withdrawal is a claim like any other and
> wants the same standard of evidence as the thing it withdraws.

#### Then the ladder, and nothing on it converges to the reference

**Re-taken at a fixed 5 505 steps**, like everything else in this section — the first
version of this table let the count be derived and so is not comparable with the
reference rows above it:

| torn elements | stopping rule | 768 | 3 072 |
|---|---|---|---|
| CPU double, the reference | 1e-16 | **32** | **162** |
| control: double, mesh jittered 2e-7 m | 1e-16 | **32** | **162** |
| `float`, as shipped | 1e-9 | 40 | 247 |
| `solve` — + the 7×7 solve in fp64 | 1e-9 | 44 | 243 |
| `condense` — + the condensation in fp64 | 1e-9 | 44 | 243 |
| `tight` — float, the CPU's stopping rule | 1e-16 | 40 | **205** |
| `newton` — the **whole block** in fp64 | 1e-16 | 44 | **204** |

**Sorted by stopping rule rather than by precision, because that is what the numbers
sort by.** At 3 072 the three kernels on the shipped 1e-9 gate land on 243–247 and the
two on the CPU's 1e-16 land on 204–205, whatever their arithmetic: the gate is worth
**42 elements** and everything fp64 does inside a gate is worth **four**. At 768 the
gate is worth nothing at all and fp64 is worth four, in the wrong direction.

**And `tight` and `newton` agree to one element at 3 072** — 205 against 204 — where the
derived-step table this was drawn from had them 36 apart (205 against 241). That is the
section's own central claim about this ladder, *"the whole of the fp64 block's effect on
this kernel is the tolerance it was bundled with"*, and fixing the step count is what
made it hold at the size where tearing is busiest.

> **It does not hold at 768, and a previous revision of this paragraph said it did.**
> The sentence read "`tight` and `newton` now agree to one element **at both sizes** —
> 40/44 and 205/204", quoting as agreement a pair four elements apart in the table
> immediately above it, one paragraph after that same four had been called out as
> "everything fp64 does inside a gate". Both figures re-measure exactly —
> `zone_gpu_probe --radius=2.5 --sub=8 --steps=5505 --mapping=workgroup --eas=tight`
> gives 40 and `--eas=newton` gives 44 — so this was never a drifted number, it was a
> claim that contradicted its own evidence on the page. The honest statement of the
> finding is the one the table already makes: what the numbers sort by is the stopping
> rule, and inside a stopping rule fp64 is worth about four elements either way.

**Read the spread, not any one row.** Five kernels that differ only in the precision of
the enhanced block land between 40 and 44 at 768 and between 204 and 247 at 3 072, while
the reference is 32 and 162 and the negative control hits both exactly. The variation
*among the precisions* is about **half** the gap they were supposed to close at both
sizes, and its sign reverses with size: fp64 makes it worse at 768 and better at 3 072.
A quantity that responds to a variable does not change which way it responds when the
mesh is refined.

Over the full five sizes the fp64 block changes the torn count nowhere it was zero.
**Re-taken at the fixed step counts this column claims** — the previous version of
this table gave 5 505 / 5 505 in its own `steps` column while carrying figures from
runs of 5 513 and 5 545, which is the same defect as the arrows above and in the one
place it was self-evidently wrong:

| elements | steps | CPU torn | float torn | fp64 torn | float diss. rel. | fp64 diss. rel. |
|---|---|---|---|---|---|---|
| 192 | 5 505 | 0 | 0 | 0 | 2.79e-1 | 2.79e-1 |
| 768 | 5 505 | **32** | 40 | 44 | 2.57e-1 | 2.69e-1 |
| 3 072 | 5 505 | **162** | 247 | 204 | 3.38e-1 | 3.01e-1 |
| 8 192 | 1 500 | 0 | 0 | 0 | 9.32e-2 | 9.32e-2 |
| 16 384 | 1 000 | 0 | 0 | 0 | 2.28e-2 | 2.28e-2 |

#### Why it cannot help, which was visible before it was built and is now measured

Kaa is `∫GᵀCG dV` and the residual is `∫Gᵀσ dV`. **C and σ are produced by a float return
map into float shared memory.** Widening the arithmetic that consumes them cannot recover
digits they never had, so the enhanced block's accuracy is capped at float's ~1e-7 on its
own inputs however it is compiled. That cap was not binding when §8 wrote item 3, because
κ(Kaa) was then (h/t)⁴ ≈ 10⁷ and the block genuinely had no correct digits; item 2's
normalisation has since put κ(Kaa) at a constant **3.50**, and a well-conditioned 7×7
solve on float inputs already delivers everything float can deliver.

Measured directly, on the same 400-step unit-scale fixture the two mappings are compared
on: holding the stopping rule fixed, the **entire** enhanced block in fp64 moves alpha by
**4.6 × 10⁻⁹** where the float kernel is already **3.8 × 10⁻⁶** away from the double CPU.
A factor of **810**. `tests/test_zone_gpu.cpp` asserts that ratio at 100, so it fails if
fp64 ever does start mattering here.

#### The one thing that did move, and §8 had the mechanism wrong

Item 2 recorded: "measured at 100 steps … the GPU's enhanced parameters were **exactly
zero** — the residual gate fired on the first iteration and left them at their warm
start". The zero is real and reproduces. **It is not the residual gate.** `kResidualTol`
is already the CPU's `1e-12·σ_y V`; what differs is `kEnhancedWorkTol`, **1e-9 against
the CPU's 1e-16** — and the shader tests it *after* computing the correction and *before*
applying it, so the element returns alpha at its warm start of zero. Measured on the
192-element patch, and on the unit-scale fixture in the suite:

| | CPU double | `float` | `solve` | `condense` | `tight` | `newton` |
|---|---|---|---|---|---|---|
| peak \|alpha\|, 192 elements at 50 steps | 7.671e-9 | **0** | **0** | **0** | 2.010e-8 | 2.010e-8 |
| peak \|alpha\|, unit fixture, 120 light steps | 5.511e-9 | **0** | **0** | **0** | 5.206e-9 | 5.206e-9 |

Three things fall out of that table and each is a separate finding.

- **fp64 in the solve and in the condensation changes nothing whatsoever** — they return
  the identical bit-zero, because the correction they compute to sixteen digits is
  discarded at the same gate the float kernel discards it at.
- **`tight` and `newton` are indistinguishable *on alpha*.** Every digit of the two rows
  above agrees. They are not indistinguishable on the torn count — 40 against 44 at 768,
  205 against 204 at 3 072 — so the fp64 arithmetic is worth about four elements and the
  tolerance is worth the rest. The whole of the fp64 block's *large* effect on this
  kernel is the tolerance it was bundled with.
- **The tolerance is worth something on its own.** In float, with the CPU's rule, alpha
  goes from wrong by 100% of the reference (it is zero) to wrong by 20%, and the peak
  lands within 5.5%. It is also the only variant that moves the torn count meaningfully:
  **247 → 205 at 3 072**, which takes it from 52% over the reference to 27%. It does
  **nothing at 768** (40 → 40, 25% over and still 25% over), it costs **4.5×** on the
  kernel, and 27% over is still not an answer — so the shipped gate stays at 1e-9 and
  this is recorded rather than adopted. It is the one thing on this page that a future
  attempt should start from, and it is not a precision question at all.

  > **Three of those numbers were wrong, and all three the same way.** The bullet read
  > "248 → 205 … 53% over the reference to 27% … not reproduced at 768 (41 → 40, 28%
  > over to 25%) … costs 3.5×". 248 and 41 are the **derived**-step float counts (5 545
  > and 5 513 steps); 205 and 40 are the **fixed** 5 505-step `tight` counts, which is
  > what every table in this subsection is taken at. Comparing across the two is the
  > exact defect the *correction to the correction* above is about, in the prose of the
  > section that documents it. At a fixed 5 505 the shipped float kernel tears 247 and
  > 40, so the tight rule is worth 42 elements at 3 072 and **nothing at all** at 768 —
  > which is what the ladder table two screens up already said. The 3.5× is the
  > figure the cost table below carries in its *"as first published, under
  > contention"* column; the current one, in its left-hand column, is 4.5×.

**A caution on comparing alpha's *value* across item 2.** Normalising the enhanced modes
rescales alpha by construction — column *j* of G by *s* and α*ⱼ* by 1/*s*, which is what
makes the fix exact — so the numbers are not commensurable either side of it and §8's
"the CPU's were 3.9e-8" should not be compared with the 3.445e-8 measured today at the
same 100 steps. What *is* commensurable is whether the device's alpha is zero, and there
item 2 did move something: at 100 steps the device now returns 6.202e-8 rather than
nothing. The bit-zero reproduces at 50 steps and below. That crossover shift is the one
thing item 2 changed on this path, and it is a change in where the work gate stops
firing, not in how many digits the block has.

#### The cost, which is not small

`zone_gpu_probe --stats`, so the compiler's account rather than the curve's:

| kernel | registers/thread | spill over calibration | shared |
|---|---|---|---|
| `solidshell_forces_wg.comp` | 64 | 96 B | 2 832 B |
| + fp64 7×7 solve | 72 | 128 B | 2 928 B |
| + fp64 condensation | 80 | 128 B | 2 928 B |
| + fp64 alpha and Newton | 80 | 128 B | 2 960 B |

Derived from the register count and Pascal's 65 536-register file, warp occupancy falls
from 32 warps per SM to 25. That is the cheap part. The expensive part is the arithmetic
and the extra Newton iterations the tighter rule buys:

Kernel time at 768 elements, two interleaved passes of all five, relative to `float`,
minima:

| kernel | × `float` | *as first published, under contention* |
|---|---|---|
| `solve` | **×2.4** | ×1.8 |
| `condense` | **×4.4** | ×2.9 |
| `tight` — float, the CPU's stopping rule | **×4.5** | ×3.5 |
| `newton` | **×14.7** | ×6.5 |

> **All four were understated, and `newton` by a factor of 2.3.** This pass was itself
> taken while a sibling agent was building — 25–70% of the box in other hands, recorded
> per row — so it is quoted for *ratios* and not for absolute times: all five kernels
> ran in one interleaved pass and were equally exposed. The check that this is
> permissible is `newton`, which was also measured on the idle sweep above: **×14.7
> here against ×14.4 there**, 2% apart. A ratio taken across five kernels in one pass
> survives contention; the ratio that does not is the one taken across two *sessions*,
> which is what produced the ×6.5.

And `newton` against `float` across the whole size range. **Re-measured on a box
checked idle** — the first version of this row was taken against an unrelated GPU
consumer at 100% utilisation and said so; this is the re-run it asked for, four
interleaved repeats per cell, minima quoted:

| elements | 192 | 768 | 3 072 | 8 192 | 16 384 |
|---|---|---|---|---|---|
| `newton` ÷ `float`, kernel time | **×14.5** | **×14.4** | **×10.3** | **×5.3** | **×5.2** |
| *the same row, taken under contention* | ×10.4 | ×8.5 | ×9.9 | ×5.3 | ×6.0 |
| float mapping, against the CPU | 1.56× | 1.27× | 1.34× | 2.23× | 2.44× |
| **measured for the fp64 kernel** | **0.11×** | **0.09×** | **0.13×** | **0.42×** | **0.47×** |

> **The contended session understated the cost, and the mechanism is worth keeping.**
> The three larger sizes agree to within 13%; **768 elements is out by 70%** and 192
> by 39%. That is not noise in one direction by luck. A competing consumer inflates a
> *short* kernel proportionally more than a long one — latency it can hide behind its
> own work is latency the short kernel is made of — so it inflates the `float`
> denominator far more than the `newton` numerator and the ratio comes out too small.
> The previous revision noticed the second half of this ("`newton` is the steadiest of
> the five, because at six times the work it is the least sensitive") and drew the
> reassuring conclusion from it. The correct conclusion is the opposite: **if only one
> side of a ratio is insensitive to contention, the ratio is biased, not robust.**
>
> The last row is now **measured** rather than *implied*: `newton`'s own wall time
> against the CPU minimum for the same size and step count, out of the same sweeps,
> rather than the float figure scaled by the kernel ratio. It is worse than the
> implied 0.14–0.43× at the small end and slightly better at the large one.

**Five to fifteen times, and it is enough to undo the remap.** The fp64 block puts the
kernel at **0.09–0.47× the 23-worker CPU** — lower at the small end than the
one-invocation mapping's 0.22–0.64×, and never better than half the CPU anywhere.

§8 item 3's estimate — "even at Pascal's 1/32 fp64 rate it is a small share of the
kernel" — is wrong twice over: the enhanced block is ~70% of a Newton iteration by
§8's own phase table in *The remap, and what it did*, not a small share; and the
stopping rule that has to come with it multiplies the iteration count on top.

#### Where the error is instead, which the same runs place directly

The two zeros above do more than explain themselves. **At 50 steps and below the shipped
kernel's enhanced modes are switched off entirely** — alpha is bit-zero on every element,
so the enhanced block contributes exactly nothing to the answer. And at 100 steps, still
before a single point has yielded, the device's node field is already **1.41 × 10⁻⁶ m
away from the CPU's, 5.8 × 10⁻² of the 2.4 × 10⁻⁵ m the punch has travelled**.

So the divergence is well established in a regime where the enhanced block is provably
not participating. It cannot be the enhanced block. What is left is `u = Rᵀx − X`, the
polar decomposition, the return map and the integrator — per-step float error injected
thousands of times, and biased, exactly as *Precision* above describes.

**So the item closes as a negative.** Alpha in double was the last thing on this list, it
has been tried, and the digits the enhanced block was missing were not what stands
between this kernel and its answer. The enhanced block was simply not the part that
needed the digits.

### What mutation testing found on the fp64 block

> **The re-measurement pass that produced the tables above ran no mutants, and that is
> deliberate rather than an omission.** It changed no C++ and no GLSL: what it added is
> `tools/zone_gpu_probe/sweep.sh`, a driver that runs the whole sweep in one invocation
> and writes a CSV with a run id, a commit and the step count on every row, plus five
> deterministic cells of this section in `scripts/check-figures.sh`. There is nothing
> for a mutant of the solver to be detected by that was not already there. The level-0
> SPIR-V is unchanged and that was checked rather than assumed — `776d15e`'s shader
> source and today's compile to byte-identical output under the same `glslc -O`, 36 480
> bytes, `cmp`-clean, and identical to the build artifact the probe loads.

Forty-one mutants — the twenty-six from the remap, brought up to date, plus twelve of
the fp64 block and its host wiring, and three deliberate controls.
**29 of 38 real mutants killed; 3 of 3 controls survive.** `tools/zone_gpu_probe/mutate.py`.

> **The harness now mutates a copy of the tree outside the repository.** Two harnesses
> here have left a mutant applied in the working tree after being killed before their
> cleanup ran, and a `finally` does not survive a kill. A copy cannot have that failure
> mode at all, and it also makes the harness safe to run while something else is
> building the real tree — which is how this pass and the gate ran together.

**Four of the survivors are the same finding as the measurements above, arrived at from
the other end**, and they are why this pass is worth reading next to the tables.

- **Dropping the fp64 warm start survives. So does never writing the fp64 alpha back.**
  Two independent mutants of the persistent fp64 alpha buffer — the *literal* reading of
  "keep alpha in double" — and nothing observes either. With the CPU's stopping rule the
  enhanced Newton converges to the same alpha from a cold start as from a warm one, so
  carrying it in double between steps buys nothing that any quantity this repo measures
  can see.
- **Making the fp64 condensation read its Gauss weights from the wrong offset survives**,
  and this is the sharpest statement of *why* fp64 in the condensation cannot help.
  **Kaa is the enhanced Newton's Jacobian, not its residual.** The root is set by
  `∫Gᵀσ dV`; Kaa only decides the path to it. A wrong Kaa changes how many iterations
  the Newton takes and not what it converges to — so widening Kaa to fp64 cannot move
  the converged alpha either. The mutant found the argument the measurement had only
  the shape of.
- **Cutting the tight variant's iteration cap from 40 back to 12 survives**, which says
  the Newton reaches the CPU's 1e-16 work gate inside twelve iterations. The cap is not
  what the tight rule costs; the extra iterations below it are.

**Two survivors are hazards no test on this machine can see, and both are the same
shape as §8's barriers.**

- **Three of the four `barrier()` deletions survive**, exactly as recorded before: the
  1070 Ti's subgroup is 32 and so is the workgroup. Reproduced rather than repeated.
- **Never requesting `shaderFloat64` of the device survives.** This driver runs a shader
  carrying the `Float64` capability on a device created without the feature enabled.
  That is undefined behaviour the spec forbids and this driver permits, so the request
  is correct by the specification and by nothing the suite can observe. An edit removing
  it will pass everything here and fail on a stricter driver.

**One survivor is §8's own reversal, reproduced.** Replacing the shader's Jacobi
equilibration of Kaa with the identity survives, as it has since `computeRestForms`
normalised the enhanced modes.

> **And one control was killed, which is the third time in this repo that a control
> behaving unexpectedly has been a defect in the *harness*.** Rewriting `mutate.py` for
> this pass mistyped the "write `gl_WorkGroupSize.x` as the literal 32" control: the
> literal went into the loop's *bound* rather than its *stride*, so it ran eight
> iterations past `kDof` and wrote off the end of a shared array. Not an equivalent edit
> at all, and the suite killed it in nineteen checks. Nothing but the control's own
> recorded expectation would have said so — a mutant that dies looks like a suite doing
> its job. Fixed, and it survives.
>
> A second control was killed once and survives on a re-run, and the reason is worth
> recording separately because it is not about this work: it failed on
> `test_promotion.cpp`'s **"two zones cost about two zones"**, a wall-clock ratio
> asserted between 1.4 and 3.0. Under an unrelated GPU consumer and forty mutants'
> worth of load, that ratio can leave its band. The gate's six repeat runs did not
> catch it and neither did the other forty mutants. It is a pre-existing flaky
> assertion, not a consequence of anything here, and it is recorded rather than
> quietly re-run.

### What mutation testing found on the remap

Twenty-eight mutants of the new kernel and its host — `tools/zone_gpu_probe/mutate.py`,
which is checked in because the result below is only readable next to the mutants
that produced it. Twenty-six real ones and two deliberate controls.

**21 of 26 killed; both controls survived as controls should.**

**The test that does the killing is the one that did not exist before: the two
mappings against each other.** Two float kernels running the same arithmetic agree
to 2.8 × 10⁻⁵ on the plastic dissipation, where either against the double CPU is
6.8 × 10⁻². That is three orders of magnitude of headroom, and it is where
`kRoot23 wrong by 0.4%` now dies — **the mutant this section previously recorded as
an unkillable survivor**, on the grounds that "the instrument is less precise than
the error it would have to catch". A second float kernel is a more precise
instrument, and it was available for the cost of building the mapping anyway.

**Two survivors were holes in the new tests, and both are fixed.**

- **Collapsing the host's shader selection so both mappings load the same file
  survived everything.** The cross-mapping comparison would then have been one
  kernel against itself — perfect agreement, every tolerance passed, a completely
  vacuous test. This is `CLAUDE.md`'s "two solvers that both did nothing agree
  perfectly" in its most literal form, and the vacuity guards that were there
  (the patch moved, it yielded, alpha was live) could not see it because all of
  them were true. `ZoneGpuSolver::elementShader()` exists so the guard can be an
  equality on what was actually loaded.
- **`elementOut`'s second float per element was written by both shaders and read by
  nobody.** The torn flag the host uses comes from `plastic[kStateTorn]`; this was a
  parallel copy that had drifted out of the readback path. Setting it to the wrong
  value survived because nothing consumes it. It is gone, and the buffer is one
  float per element.

**Three survivors are the same finding: on this device the barriers are
unobservable.** `solidshell_forces_wg.comp` declares a 32-thread workgroup and the
1070 Ti's subgroup is 32, so a workgroup is exactly one subgroup running in
lockstep — measured, via `zone_gpu_probe --stats`, not assumed. Of four mutants
that each delete one `barrier()`, three survive the whole suite. The fourth dies,
and the reason is worth stating precisely because the obvious one is wrong: not
that the hardware failed to synchronise, but that `barrier()` is also a *memory*
barrier and without it the compiler kept a shared value in a register. **So which
of the four dies is a fact about this driver's optimiser, not about the code.** The
barriers are required by the Vulkan memory model and become load-bearing on any
device with independent thread scheduling; they are correct by the uniformity
argument in the shader's header and by nothing else, and an edit that removes one
will pass every test in this repository.

**One survivor is a fact about the element, and it is the clearest single piece of
evidence for the conditioning result above.** Replacing the shader's Jacobi
equilibration of Kaa with the identity now **survives** — where §8's original
mutation pass recorded it as *killed*, "which says the fix earns its place". The
element-level normalisation has since made it redundant, exactly as the A/B on the
full run says. Two instruments, one conclusion, arrived at from opposite ends.

**One survivor is a tolerance nothing tests.** Dropping `yieldStrength` from the
enhanced Newton's convergence scale makes the gate about 10⁸ times tighter; the
Newton runs longer and converges to the same answer, so nothing sees it. That is
a real gap rather than a false alarm — the scale is what makes the criterion
dimensionless — but it is a gap in what can be observed from outside the element,
and no assertion on an integral quantity will close it.

> **Two defects in the harness, and both were caught by the same thing — a control
> behaving unexpectedly — rather than by reading it.** It rebuilt only the shaders
> for a shader mutant, which is right in isolation and wrong in sequence: after a
> mutant that edited the *host*, restoring the file on disk does not relink the test
> binary, so the next mutant ran against the previous one's host object. And an
> aborted run left a mutant applied in the working tree, which is `CLAUDE.md`'s
> "the mutant was left applied in the tree by a harness accident" happening for the
> second time. The restore is in a `finally` now and the build is unconditional.
> **A harness that scores its own controls is what makes both of these visible;
> without them the first shows up as two controls surviving, which is what they
> were supposed to do.**

### What mutation testing found the first time round

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

> **Superseded, and the way out was not more precision but a second instrument.**
> This mutant dies now. The remap produced a *second float kernel*, and two float
> kernels agree with each other three orders of magnitude better than either agrees
> with double — so the comparison that could not see 0.4% against the CPU sees it
> easily against the other mapping. The argument above was right about the
> instrument and wrong to conclude the kernel could not be validated: what it
> needed was a differential reference at its own precision, not a better absolute
> one.

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

> **The second of those has since reversed, and the reversal is the point.**
> Dropping the equilibration now *survives*, because `computeRestForms` normalises
> the enhanced modes and κ(Kaa) no longer depends on the element's shape. The fix
> earned its place when it was written and has stopped earning it since, which is
> a thing a mutation score can tell you and a passing suite cannot.

### What to do about it

1. **Keep the `RestForms` hoist.** It is 2× on the CPU, free, and bit-identical.
   It is the one part of this work that is unambiguously worth having.
2. ~~**Fix the conditioning in the element, not in the shader.**~~ **Done.**
   `computeRestForms` normalises each column of G by its own weighted L2 norm.
   Measured on this element, κ(Kaa) against slenderness:

   | h/t | before | after |
   |---|---|---|
   | 5  | 2.19e3 | 3.50 |
   | 10 | 3.50e4 | 3.50 |
   | 30 | 2.84e6 | 3.50 |
   | 60 | 4.54e7 | 3.50 |

   Before, that is **(h/t)⁴ to three figures**, which is the diagnosis §8 reached
   from the float failure. After, it is constant in the geometry, and the residual
   3.50 is exactly the material's own anisotropy — 277 GPa against 79.2 GPa — so
   the geometric spread is entirely gone rather than merely reduced. No material is
   needed to do it, which is why it fits in the rest forms.

   Exactness is asserted as an identity, not inferred from the suite still passing.
   The test re-derives the condensation independently from the public `RestForms`
   and its own constitutive matrix, then repeats it with the enhanced basis rescaled
   across twelve orders of magnitude: Kua → Kua S and Kaa → S Kaa S, so
   Kua S (S Kaa S)⁻¹ S Kuaᵀ collapses back exactly. Wrecking κ to **2.86e23** moves
   the element by 1.91e-06 of 1.03e12 — a relative 1.9e-18. Both vacuity guards are
   there: that the rescaling really did ruin Kaa, and that the independent
   condensation reproduces the shipped element.

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
3. ~~**Keep alpha in double — this is now the only thing left to try.** The EAS block
   is 7×7. Even at Pascal's 1/32 fp64 rate it is a small share of the kernel, and it
   is the only part that has been shown to need the digits.~~ **Done, and it is a
   negative.** *Alpha in double — measured, and it is not what was missing*, above, has
   the tables. The short version: five kernels from one source, differing only in how
   much of the enhanced block is fp64, land between 40 and 44 torn at 768 elements and
   between 204 and 247 at 3 072, against a reference of 32 and 162 that the negative
   control reproduces exactly — a spread about half the gap, whose sign reverses
   between the two sizes, and which sorts by the stopping rule each kernel carries
   rather than by its arithmetic. Holding the stopping rule fixed, the whole block in
   fp64 moves alpha by
   1/810 of the amount the float kernel is already wrong by, because Kaa's inputs are a
   float tangent and a float stress and widening what consumes them recovers nothing.
   It costs **5–15× on the kernel**, which turns 1.26–2.43× against the CPU into
   **0.09–0.47×** — re-measured on an idle box, where the first pass was taken
   against a GPU consumer at 100% and reported 5–10× and 0.14–0.43×. Contention
   inflates a short kernel more than a long one, so it inflated the `float`
   denominator and understated the cost, by 70% at 768 elements. Two claims in this
   item were wrong: the block is ~70% of a Newton
   iteration rather than "a small share", and it was not "the only part shown to need
   the digits" — that was true when κ(Kaa) was (h/t)⁴ and stopped being true when
   item 2 made it 3.50.
4. ~~**Re-map the kernel to a workgroup per element** before measuring throughput
   again. The current numbers measure register spilling, not the element.~~
   **Done**, and the diagnosis held: `solidshell_forces_wg.comp`, 1.26–2.43× against
   the CPU where the invocation mapping was 0.23–0.68×, and 484 floats of spill per
   thread down to 24. Every precision figure is unchanged by it. **Re-measured on a
   box checked idle, six interleaved repeats: 1.27–2.44× and 0.22–0.64×**, so the
   pair holds to 1% and the spill figures — which come from the driver's compiler
   rather than a clock — are identical.

**The status of this item has changed twice, and the second change is that the list is
now empty.** Throughput was the visible problem and it is solved. Precision was the one
that decided whether the element belongs on a GPU at all, and it is **unchanged**: the
float kernel tears a quarter to a half too many elements, which is exactly the output a
zone exists to produce. Item 3 was the remaining hypothesis about *why*, and it is now
disproved rather than untried — the enhanced block is well conditioned, float already
resolves it, and putting it in fp64 costs 5–15× to move the answer by less than the
noise between five kernels that all claim to compute it.

So the CPU remains the *trustworthy* path for Tier 2, it is no longer the faster one,
and the two statements should still not be made together as though one implied the
other. What has changed is that there is no longer a cheap experiment outstanding.
**Anything further would have to widen the parts that actually carry the error** — `u`,
the return map, the integrator — which is the whole element in fp64 on a device that
runs it at 1/32, and the measurement above is the reason not to start.
