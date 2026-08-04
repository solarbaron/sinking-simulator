# 02 — Simulation

The physics plan, in dependency order. Section 1 is implemented; the rest is the
build.

---

## 1. Hydrostatics, flooding and air — **implemented**

See `README.md` for behaviour and `engine/sim/ship.cpp` for the code. Summary of
what exists and what it still needs.

**Done:** closed-mesh volume/centroid integration under an arbitrary plane;
internal free surfaces re-levelled per tick; floodwater as real mass at its real
centroid (free-surface effect emergent, validated against ρ·i/Δ); compressible
trapped air per compartment with isothermal state and pressure-limited transfer;
a single two-phase orifice network covering breaches, doors, hatches, vents,
pipes and cross-flood ducts; bilge and ballast pumps with head-dependent output;
6-DOF rigid body with measured hydrostatic stiffness, added mass, quadratic drag
and modal damping; damaged GZ curves computed by forced-heel sweep with the
floodwater free to re-level.

**Next, in order:**

1. ~~**Compartment geometry from CSG rather than authored boxes.**~~ **Done.**
   `clipByPlane` / `clipToBox` carve each compartment out of the hull interior:
   Sutherland-Hodgman per triangle, cut edges welded on a spatial hash, chained
   into loops and capped by ear clipping. Validated by cutting a hull into a grid
   and checking the cell volumes sum to the whole. Still to do: exact predicates
   (Shewchuk) for robustness on degenerate input, nested cap loops (a plane
   through a hollow mast), and non-planar bulkheads.
2. **Permeability from contents, not a constant.** Currently a scalar per
   compartment. Should be derived from the actual cargo, machinery and outfit
   volumes placed in the space, and should change as cargo shifts or burns.
3. **Sloshing.** The quasi-static horizontal free surface is right for slow
   flooding and wrong for a half-full tank in a seaway. Escalation path in §6:
   the lumped model stays until a compartment's excitation exceeds a threshold,
   then it is promoted to a particle solver and demoted when it settles.
4. **Downflooding angle bookkeeping.** Track the heel at which each opening
   immerses, and surface it — this is the number that decides survival and it is
   currently only implicit.
5. **Progressive flooding through the structure.** ~~Requires §3: once elements
   tear, new orifices appear automatically.~~ **The coupling is done** —
   `engine/sim/breach.{hpp,cpp}`, written up in §3 below: given a structural mesh
   and the panels that failed, it returns `Opening`s that drop into
   `Ship::openings` and are indistinguishable to the solver from authored ones.
   What is still missing is the other half — a fracture model to decide *which*
   panels fail. Today the failure set is supplied by the caller.

---

## 2. Seakeeping and hydrodynamics

The ship currently floats in still water. This adds the sea. The wave field is
built; the loads it drives are not.

### Wave field — **implemented**

`engine/sim/waves.{hpp,cpp}`. A directional spectrum discretised into components
and evaluated analytically:

```
eta(x, y, t) = sum_i a_i cos(k_i . x - omega_i t + phi_i)
```

`SeaState` asks for the sea in the terms a forecast gives it — Hs, Tp, mean
direction, spreading exponent, peak enhancement — and a seed. `WaveField` turns
that into N × M components; `elevation()` and `kinematics()` evaluate them
directly at whatever point is asked for. A `WaveField` can be built from several
superposed `SeaState`s, which is how a swell train sits under a wind sea; the
variances add and each system draws from its own phase stream.

**JONSWAP normalised by quadrature, not by the usual approximation.** The shape
is scaled so that its zeroth moment is exactly Hs²/16 whatever γ is. The
customary DNV shortcut for that factor, `1 − 0.287 ln γ`, is a fit and is wrong by
a few parts in a thousand (worst observed 7.9 × 10⁻³ at γ = 7); integrating the
shape costs one quadrature at construction and is exact. γ = 1 recovers
Pierson–Moskowitz, and the normalisation then comes out as 1.0 to 8 × 10⁻¹⁶,
which is a useful check that the quadrature is not lying.

**Frequencies are spaced by equal energy, not uniformly.** The spectrum is cut
into N intervals each carrying exactly m0/N of the variance, and each interval is
represented by its own energy *centroid* frequency. The outermost intervals are
open — the first runs down to zero and the last up to infinity — so there is no
band truncation anywhere. What follows:

- **No energy is lost, in the tail or anywhere else.** Measured Hs round trip
  (`4·sqrt(m0)` against the Hs asked for) is exact to 5 × 10⁻¹⁵ relative at
  N = 1, 8, 48, 128 and 512 — floating-point summation error and nothing else.
  Uniform spacing over a truncated band loses the tail quietly, and a sea missing
  its short waves is too smooth at exactly the scale of the openings that flood a
  ship.
- **Components sit where the energy is.** At Hs = 3 m, Tp = 9 s, γ = 3.3 and
  N = 64, the interval containing the peak is 0.0035 rad/s wide while the
  open-topped one starts at 2.7ω_p and runs to infinity; a uniform grid would
  spend most of its components resolving a tail that carries 1.5% of the
  variance.
- **The first moment is preserved exactly**, because a centroid is what makes
  `Σ ω̄_i E_i = ∫ ω S dω` an identity rather than an approximation. The mean
  period comes out at the textbook Pierson–Moskowitz ratio T1 = 0.77177 Tp to
  2 × 10⁻¹² relative.
- **The second moment is not**, and this is the price. One frequency cannot
  represent the spread of frequencies inside its own interval, so m2 is biased
  low and the zero-crossing period comes out long: +1.9% at N = 12, +0.93% at
  N = 48, +0.33% at N = 384, halving as N doubles. It is a discretisation bias,
  not an error, and the tests assert its sign as well as its size.

Pierson–Moskowitz has an exact cumulative energy distribution,
`exp(−1.25 (ω_p/ω)⁴)`, so the equal-energy interval edges have a closed form to
check the whole quantile machinery against; the numerically tabulated path that
JONSWAP needs reproduces it to 6 × 10⁻¹⁵ relative.

**Directions are discretised the same way** — M intervals of equal directional
energy, each at its centroid — over a `cos^(2s)((θ−θ₀)/2)` Longuet-Higgins
spreading function, normalised through Γ functions so it integrates to exactly 1
over the circle. That form covers the whole circle; the `cos^(2s)(θ−θ₀)`
alternative has to be truncated at ±90° and renormalised, and a spreading
function that does not normalise makes the sea the wrong height in a way that
only shows up in directional seas, which is all of them. The interval edges and
centroids are explicitly antisymmetrised about the mean direction rather than
left to cancel numerically, so a single-direction sea is exactly long-crested
instead of long-crested to 10⁻¹⁷. Spreading is currently constant with frequency;
real seas spread more at high frequency, which is a refinement, not a bug.

**Kinematics** are deep-water linear theory: ω² = gk, orbital velocity
`(a ω e^{kz} cos ψ, a ω e^{kz} sin ψ)` and the local (Eulerian) acceleration that
goes with it, which is what Froude–Krylov and Morison want. **z is clamped to ≤ 0
before the exponential.** A nonlinear Froude–Krylov integration over the
*instantaneous* wetted surface will ask for kinematics above the still-water
plane, inside a crest, and `exp(kz)` there grows without bound — this is the
standard way that integration explodes in a steep sea. Wheeler stretching is the
refinement when the nonlinear loads arrive.

**Determinism.** Phases come from Threefry-2x64-20 keyed on (seed, sea-state
index) and *indexed by component*, never advanced as a stream — so components can
be generated in any order or on any thread and the sea is bit-identical. The
residual risk is the one §4 of `01-architecture.md` already names: `cos`, `exp`
and `pow` here are still libm's, and move to the in-tree implementations when
those land.

**Cost, measured (one core, `-O2`, glibc `libm`).** Construction is ~2 ms per sea state,
one-off. Evaluation is ~10–12 ns per component for elevation and ~20 ns for full
kinematics, so a 576-component sea costs ~6 µs per elevation query.

The extrapolation from that figure — "a 2000-panel hull at 100 Hz is 23 ms/tick
against a 10 ms budget, so the work is a vectorised sincos" — was measured once
the Froude–Krylov coupling existed, and **the number was right while the
prescription was wrong.** The whole-tick figures for the 120 m ferry (1196
triangles, 600 vertices) came out at 22.5 ms with the default 48×12 spectrum,
and the wave field accounted for essentially 100% of it. But the first fix was
not vectorisation:

| | queries per buoyancy integration | ms/tick, 128 components | ms/tick, 576 |
|---|---|---|---|
| height evaluated per triangle corner | 3948 | 5.19 | 22.47 |
| height evaluated per **vertex** | 960 | 2.76 | 11.87 |

Euler's formula gives `tris ≈ 2·verts` on a closed mesh, so the corner form asked
the same question about six times over. Hoisting it is arithmetically identical —
same values, same subtraction, bit-identical results — and worth 2.1× on its own.
Two further properties fell out of measuring rather than assuming:

- **Vertical plating is free.** The integrand `F = (0, 0, z − h)` points along z,
  so a panel with a horizontal normal carries no flux and never reaches the
  quadrature. On this hull that is 45% of the triangles — a ship's sides, exactly
  where a naive cost model would put the panels.
- Fitting tick cost against component count gives ~2000 elevation queries per
  step, against 960 for the buoyancy integral. The rest is compartment
  free-surface queries, which is where the next reduction is if one is needed.

A vectorised sincos is still worth roughly 4× on top and is the right next step
for large spectra; it is no longer the *first* step, and 128 components now fit
in 28% of the budget without it.

**How many components? — and a methodological lesson.** The default is 48 × 12 =
576. Two experiments were run to find out whether that is necessary, and *both
were dominated by their own noise*, which is the more useful result.

The first measured wave-train recurrence by autocorrelation. It ranked 96
components worse than 64, every "worst lag" landed near the half-record limit,
and the ranking was not monotone — the signature of an extreme-value statistic
over an estimator whose variance grows with lag, not of a real effect.

The second ran the ferry in a seaway and compared RMS motions across component
counts, using **seed-to-seed spread as the yardstick**. That framing was right
and it is what exposed the problem: over 160 s records with three seeds, seed
noise (heave 0.18–0.41 m) is as large as the variation across component counts
(0.53–0.75 m). Only one finding survives it, and it is consistent in sign across
every comparison:

| | roll, RMS deg | pitch, RMS deg |
|---|---|---|
| long-crested (M = 1) | 1.65 – 2.33 | 2.04 – 2.49 |
| spread (M = 8) | 2.49 – 2.98 | 1.61 – 1.73 |

**Directional resolution matters and frequency resolution, above about 16 bins,
does not.** A long-crested sea under-predicts roll by 30–40% and over-predicts
pitch, because spreading moves energy off the head-on components that drive pitch
and onto the oblique ones that drive roll. Equal-energy binning already
reproduces Hs to within 2% at N = 8, so the frequency axis was never the
constraint.

The first version of that study also computed heave RMS about zero rather than
about the mean, so it reported the ship's static draft (6.15 m) rather than its
motion. It looked beautifully converged across every configuration, for the
obvious reason.

The lesson worth keeping: **a statistical instrument answered neither question,
while the deterministic one answered both.** An RAO needs a single regular wave,
has no random phase and therefore no seed noise, and is checkable against closed
forms. The response spectrum then follows as |RAO(ω)|²·S(ω) without simulating a
random sea at all. Component count is consequently a *recurrence and visual
variety* question — how long before the sea visibly repeats — not an accuracy
one.

**Still to build here:** the FFT-on-a-tiled-grid path for rendering (physics
keeps the analytic sum — interpolating a rendering grid loses exactly the detail
that decides whether a hatch immerses); second-order Stokes correction for steep
seas; and a breaking criterion driving spray/foam emission.

### Wave field — what the tests are pointed at

`tests/test_waves.cpp`. Any sum of cosines looks like the sea, so nothing is
eyeballed: every assertion is a closed form fixed before the code ran — the
analytic integral of S, the dispersion relation, the Pierson–Moskowitz period
ratios from Γ functions, the exact Airy wave, `e^{-π}` at half a wavelength down,
`E[cos Δ] = s/(s+1)` for the spreading.

Two of them earned their place by catching something:

- **The variance-convergence test was wrong before the code was.** It sampled
  elevation on a regular grid stepping 101 m through a sea whose dominant
  wavelength was 126 m, so the dominant components were aliased rather than
  averaged and the sampled variance came out 1% high — systematically, every run.
  The fix was pseudorandom sample points, not a looser tolerance. It now takes
  300 000 samples for a relative standard error of `sqrt(2/M)` = 0.26%, measured
  at 0.22% RMS over 30 field and sampling seeds, and asserts 1%.
- **A one-component sea must be the exact Airy wave**, `a cos(kx − ωt + φ)`, with
  the kinematics that go with it. A spectrum averages sign and convention errors
  into something that still looks like a sea; a single component cannot.

The suite was checked by mutation: seventeen plausible-but-wrong implementations
— g inverted in the dispersion relation, the spreading left unnormalised, the
depth decay with the wrong sign, `+ωt` instead of `−ωt`, the missing factor of
two in the amplitude, one round of RNG mixing instead of twenty, the band
truncated at 3ω_p — and each must fail at least one check. Two survived the first
pass and both were real gaps: placing components at interval *midpoints* instead
of energy centroids was invisible until the mean-period tolerance was tightened
from 0.1% to 10⁻⁸ (which is what the construction actually guarantees), and
swapping JONSWAP's σ = 0.07/0.09 branches changed nothing measurable at all until
a test was added that divides S by the bare Pierson–Moskowitz shape and demands
the peak enhancement back. Both of those are shape errors that leave Hs, Tp and
the total energy correct, which is precisely why they needed to be hunted rather
than waited for.

### Radiation and diffraction

The honest way to do wave loads on a large ship is potential flow.

- **Radiation — implemented by strip theory**, `engine/sim/radiation.{hpp,cpp}`.
  See the next section for what it is, what it measures, and where it stops.
- **Diffraction — still to build.** The exciting force from the incident wave
  scattering off a stationary hull. Strip theory can supply it by the same
  sectional route (Haskind relations tie it to the radiation potential already
  computed, which is the cheap way in), and nothing of it exists yet.
- **The offline BEM is still the plan for the real thing:** solve the
  boundary-element problem over the whole hull for a grid of frequencies and
  headings with NEMOH or Capytaine, and ship the coefficient tables as ship
  assets. `RadiationTable` is deliberately the *only* thing the Cummins machinery
  consumes, so that swap changes how the table is filled in and nothing else.
- **Nonlinear Froude–Krylov and restoring** — **implemented**. `Sea` carries an
  optional `WaveField`; when present, `Ship` transforms the hull into world
  coordinates and integrates buoyancy against the actual surface with
  `integrateBelowSurface()`, rather than clipping by a plane at a mean waterline.
  Openings read the *local* surface too, so a breach under a crest floods faster
  than the same breach under a trough — a large part of why damage in a seaway is
  worse than damage alongside. Flat water keeps the old single-plane path, which
  is exact and much cheaper, and `Sea` is implicitly constructible from a `double`
  so still-water callers are unchanged.

  Three behaviours are asserted, and the middle one is the load-bearing check:

  | Check | Why it matters |
  |---|---|
  | A zero-amplitude wave field reproduces still water | drives the entire wavy path against an answer the flat path computes a completely different way |
  | A wave ten times the ship length lifts it by the full amplitude | the ship contours a long wave |
  | A wave a fifth of the ship length barely moves it | crests and troughs cancel along the hull |

  Still-water quantities stay still-water on purpose: GZ and GM are defined about
  a mean waterline, and quoting a "GM in waves" would invent a figure naval
  architecture does not have.

  **The hull mesh must resolve the wavelength.** This turned up as a failing test
  and is worth stating loudly, because the failure mode is not blurring. A panel
  spanning several wavelengths samples the surface at three quadrature points and
  reports whatever they say, so an under-tessellated hull *invents or destroys*
  displacement — measured at ±6% on a 60 m barge under a 12 m wave, with the sign
  depending only on wave phase. That error would ride silently through every
  seakeeping result. Two panels along the length already gave the exact answer in
  that case, but the rule to build to is several panels per wavelength.

### Does any of it change the damage-control answer? — **measured**

Worth asking directly, because Phase 0's ferry scenarios are the project's
centrepiece and everything since has moved the physics under them. Two changes in
particular: the roll stiffness defect meant every ship had been rolling at 1.8×
its configured damping since Phase 0, and Ikeda's B44 is strongly
amplitude-dependent where the linear stand-in is not — so the two differ most at
exactly the heel angles a lolling ship lives at.

`shipsim --bilge-keels=40` attaches Ikeda instead of the stand-in. Comparing the
`doors` capsize trajectory over 900 s:

| | worst difference | as % of that column's range |
|---|---|---|
| heel | 0.076° | 0.13% |
| GM | 0.007 m | 0.13% |
| floodwater | 1.13 t | 0.02% |
| displacement | 1.2 t | 0.02% |
| draft | 0.0006 m | 0.03% |

All three verdicts are unchanged. (GZ moves 35% of *its* range, which sounds
alarming and is not: GZ is near zero throughout a loll, so its whole range here is
7 mm and the difference is 2.4 mm.)

**The conclusion is that the damage-control outcome is governed by statics, not by
roll damping**, and that is what it should be. A loll is a *static* instability —
GM goes negative and the ship settles at a new equilibrium — so how fast it gets
there barely matters. Which of the two damping models is used, and the 1.8× error
that was in both of them, move the answer by about a tenth of a per cent.

That is a reassuring result rather than a boring one. It says the Phase 0
conclusions survived having their dynamics substantially rewritten underneath
them, and it says where roll damping *will* matter: a ship rolling in a seaway
near resonance, not a ship lolling in flat water.

### Response amplitude operators — **implemented**

`engine/sim/rao.{hpp,cpp}`. An RAO is the transfer function between a regular
wave and a ship's motion, and it is the quantity model basins publish — so it is
the one place this simulator can be checked against the outside world rather than
against itself. The measurement is made the way a basin makes it: run the ship in
a regular wave, discard the transient, fit a harmonic to what is left. Nothing
linearises the ship or reads a coefficient out of the solver.

Supporting this, `WaveField` gained an explicit-component constructor and
`WaveField::regular(A, ω, direction, phase)`. The constructor **recomputes**
wavenumber and the direction unit vector from ω and direction rather than
trusting them: a caller-supplied wavenumber contradicting ω²/g would give a wave
that propagates at a forbidden speed while looking entirely normal.

The response is extracted by least squares against `[1, cos ωt, sin ωt]`, not a
DFT bin, because a record holding a non-integer number of cycles leaks — and
leaks smoothly enough to pass for physics. Two exactness claims are kept apart in
the tests: least squares recovers the frequency it is *fitting* exactly on any
window, but rejecting a *different* frequency needs a whole number of cycles.
Measured, a window short by one sample in twenty thousand moves the amplitude in
the fifth decimal place.

**What the barge sweep shows.** A 60 m × 16 m box at 4 m draft, head seas:

| ω | λ/L | heave | pitch | heave phase |
|---|---|---|---|---|
| 0.25 | 16.4 | 1.03 | 1.07 | −9° |
| 0.70 | 2.1 | 0.89 | 1.32 | +9° |
| 1.00 | 0.97 | **0.04** | 0.55 | −71° |
| 1.45 | 0.49 | **0.02** | 0.08 | −130° |
| 2.20 | 0.21 | 0.02 | 0.00 | — |

Three features are emergent rather than coded, and each is a closed form:

- **Both asymptotes.** Long waves are ridden one-for-one and short waves ignored.
  The 3% excess at ω = 0.25 is not error: it is the dynamic magnification of an
  oscillator below resonance, `1/(1 − (ω/ω_n)²) = 1.026` for a heave natural
  frequency of √(g/T) ≈ 1.57 rad/s.
- **Pitch lags the surface elevation by 96°**, because pitch follows the wave
  *slope*, which is the derivative of elevation — 90° out of phase.
- **The Froude–Krylov sinc zeros.** A constant-section box integrates the wave
  pressure along its length as `sin(kL/2)/(kL/2)`, which vanishes at `kL = 2π` and
  `4π` — wavelengths of exactly L and L/2, where the crest over one half of the
  ship cancels the trough over the other. Both notches appear, at 61.6 m and
  29.3 m against a 60 m ship. This is now the strongest test in the seakeeping
  suite: no coefficient can be tuned to produce it, and if the surface
  integration ever stops resolving the hull, or the wave's spatial phase is
  dropped, the notches are the first thing to go while both asymptotes still pass.

**Not yet valid:** forward speed. `encounterFrequency()` is implemented and tested
for sign and magnitude (ω_e = ω − ω²U cos µ / g), but `RaoSettings::forwardSpeed`
only sets the frequency the response is fitted at — it does not propel the hull.
Zero-speed RAOs are the validated case. Comparison against published RAOs for a
real hull waits on the radiation coefficients below, since without added mass the
heave natural frequency is too high.
### Radiation hydrodynamics — **implemented** (strip theory)

`engine/sim/radiation.{hpp,cpp}`. A ship carries water with it and radiates waves
away, and the force from that depends on the *history* of the motion. Cummins
(1962):

```
(M + A_inf) x'' + integral_0^t K(t - tau) x'(tau) dtau + C x = F_ext
```

Three pieces are implemented: sectional coefficients per station, the
Cummins/Ogilvie machinery that turns `A(ω)` and `B(ω)` into the memory kernel
`K(t)`, and a state-space approximation of the convolution so the memory costs a
matrix-vector product per tick instead of an integral over the past.

**No boundary-element solver was used and none is pretended.** NEMOH and
Capytaine are not available here. What replaces them is strip theory with a real
two-dimensional solve per station — not a chart lookup and not a regression.

#### The 2D section problem

Section shape is a **Lewis form**: the image of the unit circle under
`z = M(ζ + a₁/ζ + a₃/ζ³)`, determined exactly by beam, draft and sectional area
coefficient, which is what a hull table gives you. The area coefficient reduces
to a quadratic in `1 + a₃`; **only the larger root is the section you asked for**
— the smaller one lands on `a₃ = −1/2` at the circle and produces a
plausible-looking shape that is not it. `validateLewisSection()` reports both
ways the family can fail: an area coefficient outside it, and a mapping that
folds over itself (a folded contour still panels and still solves, and the answer
is meaningless).

The radiation problem on that contour is solved by a **close-fit source
distribution** in the manner of Frank (1967), with the deep-water free-surface
Green's function

```
G = ln r − ln r' − 2 Re{J} + 2πi e^{ν(y+η)} cos(ν(x−ξ)),
J = PV ∫₀^∞ e^{k(y+η)} e^{ik(x−ξ)}/(k−ν) dk = e^{Z}[E₁(Z) + iπ sgn(x−ξ)]
```

with `Z = ν((y+η) + i(x−ξ))`. Writing `J` closed-form through the complex
exponential integral is what makes this affordable — the alternative is a
principal-value quadrature per panel pair per frequency. `E₁` is evaluated as the
product `e^z E₁(z)`, which stays of order `1/z` where `E₁` alone overflows.

**Im{J} does not satisfy the free-surface condition; only Re{J} does.** The
outgoing behaviour has to come from the regular standing wave. Getting that wrong
produces a purely *real* source distribution and therefore identically zero
damping everywhere — a failure that looks like a plausible added-mass calculation
with the damping switched off, and which the near-field/far-field energy check
below catches immediately.

Port/starboard symmetry is exploited exactly: heave is solved with symmetric
source pairs and sway/roll with antisymmetric ones, so a symmetric section cannot
produce an asymmetric answer through rounding, and the cross-plane couplings are
exactly zero rather than small.

**The infinite-frequency added mass is solved separately and exactly.** As
`ω → ∞` the free-surface condition degenerates to `φ = 0`, the wave terms vanish
identically, and `G = ln r − ln r'`. That matters twice: the panel method would
resolve a very high frequency badly, and it gives `A_inf` by a route that shares
no code with the Ogilvie relation that also produces it.

#### What was measured

| Quantity | Measured |
|---|---|
| Heaving semicircle, `A₃₃(ω→∞)` | 1.0022 × ρπa²/2, against the exact 1 |
| Heaving semicircle, added-mass minimum | **0.5983** at `ω√(a/g) = 0.894`, against Ursell's published ≈0.60 near 0.9 |
| Near-field vs far-field damping | 4.5 × 10⁻³ at 40 panels, 2.2 × 10⁻³ at 80 |
| Reciprocity `A₂₄` vs `A₄₂` | 3.9 × 10⁻⁴ at 80 panels, first order in panel size |
| Ogilvie round trip `B → K → B` | worst 6.3 × 10⁻³ of peak `B₃₃` |
| `A_inf`: rigid-lid solve vs Ogilvie | 0.57% apart; the Ogilvie value varies 2.3% across `ω` |
| Closed-form transforms | 2.8 × 10⁻⁷ (B), 4.4 × 10⁻⁶ (A) |
| Memory decay | `K₃₃` falls to 1% of peak at **20.3 s**, 0.1% at 56.6 s |
| State space, 6 states | 7.6% relative RMS, **2.8% of K(0)** peak error |
| Prony on a planted 4-pole signal | 3.0 × 10⁻⁹ relative RMS |
| Geometric scaling by 2.5 | added mass exact to 7.5 × 10⁻¹⁵, damping to 1.9 × 10⁻¹² |
| Runtime cost | 0.53 µs per tick for 13 state-space models |

The **published comparison is Ursell (1949)**, the heaving semi-immersed circular
cylinder, whose added-mass coefficient dips to about 0.60 of `ρπa²/2` near
`ω√(a/g) ≈ 0.9` and returns to 1.0 at high frequency. Both the depth of the
minimum and its location are asserted, and the high-frequency limit is asserted
against the exact 1.0 — which follows from an image argument with no
hydrodynamics in it at all: with `φ = 0` on the free surface the section and its
negative image are a full circle translating in unbounded fluid.

The 2D solver was cross-checked during development against an **independent
Ursell-type multipole expansion** written for the purpose: Richardson-extrapolated
they agree to four figures (0.98745 vs 0.98749 at `νa = 0.2`). That multipole
method is not in the tree — it diverges for box-like sections, because its
Laurent series in `1/z` only converges outside a circle that a full midship
section does not fit inside. That is the reason the close-fit method is the one
that shipped.

#### Irregular frequencies — the defect that had to be found, not waited for

A source distribution over a closed contour solves the exterior problem correctly
**except at the eigenfrequencies of the interior Dirichlet problem**, where the
system is near-singular. On the 25 × 6.5 m midship section these sit at
ω = 1.33, 1.94, 2.48, 2.94 rad/s — squarely inside the seakeeping band — and they
match the closed form `ν_m = (mπ/B) coth(mπT/B)` for odd `m`. They produce
**negative damping**, a ship extracting energy from still water, and they do
**not** refine away: doubling the panels narrows the spike and leaves it there,
because they are a property of the geometry rather than of the discretisation.

They were found by the near-field/far-field energy check, not by any functional
test — every coefficient looked like a coefficient. That check compares damping
from integrating pressure over the hull against damping from the amplitude of the
wave radiated to infinity: the same number by conservation of energy, computed
along completely different paths. It reads under 0.01 in the clean band and runs
to 60 on a spike, so the signal is unambiguous.

**What is implemented is detection and repair**, not removal: solves whose energy
residual exceeds 0.2, or whose diagonal damping is negative, are rejected and
interpolated over from the neighbours that passed, and `RadiationTable` reports
`repairedSolves` out of `totalSolves` so a hull whose entire grid is being
patched cannot pass for a clean one. On the reference ferry that is **13 of 180**
section solves. The proper cure is an extended integral equation with unknowns on
an interior lid; it was attempted and abandoned — imposing `φ = 0` on a lid
without giving the lid its own unknowns over-constrains the exterior solution
too, because a source-only formulation's interior field is not zero, and the
variants with lid unknowns did not converge in the time available. It remains the
right fix.

#### Cummins and Ogilvie

The Ogilvie relations are exact identities, which makes them the strongest check
available on everything upstream:

```
K(t)     = (2/π) ∫₀^∞ B(ω) cos(ωt) dω
B(ω)     =       ∫₀^∞ K(t) cos(ωt) dt
A(ω)     = A_inf − (1/ω) ∫₀^∞ K(t) sin(ωt) dt
```

Both legs treat their input as **piecewise linear between samples and integrate
each interval against cos or sin exactly**. That is not a refinement: the
integrands oscillate arbitrarily fast at large `t`, so any rule that samples them
would need a mesh refined without limit, while the closed form is as accurate at
`t = 100 s` as at `t = 0`.

Because `A(ω)` and `B(ω)` are a Kramers–Kronig pair, `A_inf` recovered from
`A(ω)` and `K(t)` must not depend on `ω` — and must agree with the rigid-lid
panel solve, which shares no code with the transform. It does, to 0.57%, with a
2.3% spread across the band. A frequency table that was internally inconsistent
could not do that.

#### The state-space approximation

`μ(t) = ∫₀^t K(t−τ) v(τ) dτ` is replaced by `ẋ = Ax + Bv, μ = Cx` with
`C e^{At} B` fitted to `K`. `A` is block diagonal — one 2×2 block per complex pole
pair, one 1×1 per real pole — so the model is an explicit set of damped modes.
Poles come from a **least-squares linear-prediction (Prony) fit**; unstable roots
are reflected inside the unit circle, because an unstable radiation model does not
merely lose accuracy, it diverges.

Prony is exact for a signal that really is a sum of damped sinusoids, and the
tests plant one and require it back to 3 × 10⁻⁹ before the method is turned loose
on a hull. On the real `K₃₃`, six states give 7.6% relative RMS and 2.8% of `K(0)`
peak error — honest but not impressive; the published figures for this technique
are better, and the gap is the long low-frequency tail that a truncated `B(ω)`
grid leaves behind. Fitting on a coarser sampling of `K` helps (2.9% relative RMS
at `dt = 0.2 s`), which is a clue that the residual is high-frequency structure
rather than shape.

The runtime update is the **exact zero-order-hold discretisation** of each block,
in closed form. It is unconditionally stable and does not care whether the tick
matches the interval `K` was fitted at — an explicit integrator would go unstable
on exactly the stiff fast-decaying modes the fit produces. Measured at 0.53 µs
per tick for 13 models.

`RadiationForce::memoryForce()` returns `μ`; `A_inf` is **not** applied there, it
belongs on the left-hand side added to the rigid-body mass matrix, and is exposed
for that.

#### Validity limits — where this stops being trustworthy

1. **Strip theory.** The flow at a station is taken as two-dimensional, which
   needs the section small against both the length and the wavelength. A beamy
   ship cannot satisfy both comfortably; `validateRadiationHull()` flags a band
   whose shortest wavelength is under a beam, and a hull whose beam exceeds a
   quarter of its length.
2. **No forward-speed correction.** Everything here is zero-speed. The
   speed-dependent terms in strip theory (the `U/ω` corrections to the coupled
   heave–pitch and sway–yaw–roll coefficients, and the transom terms) are not
   implemented. At Froude numbers above about 0.2 this matters.
3. **Surge is identically zero.** A strip has no longitudinal radiation problem
   at all; surge added mass is entirely a three-dimensional end effect. The
   matrix reports zero rather than inventing a number.
4. **Lewis forms are a three-parameter family.** Bulbous bows, transom sterns,
   hard chines and any section with a real bulb are outside it. The mapping will
   return the closest attainable form and `validateLewisSection()` will say so.
5. **Irregular frequencies are repaired, not removed** — see above.
6. **First-order convergence in panel count.** Flat constant-strength panels with
   midpoint collocation. Halving the panel size halves the error; there is no
   cheaper accuracy to be had without higher-order panels.
7. **Deep water, no current, no viscosity.** Roll radiation damping from potential
   flow is real but small; the viscous part is the next section and is where roll
   damping actually comes from.
8. **The frequency grid is truncated.** `B(ω)` is taken to zero above the last
   grid point, and that truncation is what limits the Ogilvie round trip to
   6 × 10⁻³ rather than machine precision.

#### Cost

One section solve is O(n²) influence coefficients and O(n³) factorisation, the
former dominating. Measured on the 25 × 6.5 m midship section at ω = 0.8: **1.8–2.0
ms at 24 panels per half section and 4.7–5.8 ms at 40** — the spread is the
machine, not the method. A full table for the reference ferry — 9 stations × 20
frequencies × 24 panels — is **420–490 ms**, one-off at ship load. Runtime is
**0.53 µs per tick** for 13 state-space models, which is what the whole
state-space apparatus is for.

The cost is not uniform across sections: the inner loop is the complex
exponential integral, and its power series takes more terms for a deep section
than a shallow one, so a fine end station is cheaper than midships by more than
its panel count suggests. The quadrature over each panel drops from 8 points to 4
when the panel is well separated and the radiated waves are long, worth a factor
of two overall.

A hull whose stations are geometrically similar could share solves outright — the
2D problem depends only on `ν` times a length — and nothing exploits that yet.
That is the obvious next factor if table construction ever becomes the
bottleneck.

### Radiation coupled into the ship — **implemented, opt-in**

`Ship::attachRadiation(waterlineZ)` builds a table from the hull's *own* sections
and attaches it. `radiationHullFromMesh()` derives the stations by clipping the
hull into slabs and reading each one — sectional area is the slab's volume over
its thickness, so the same integrator that gets displacement right also gets the
stations right. On a box barge that is exact: beam and draft to zero error, area
coefficient to one ulp.

It is opt-in, and the scalar coefficients stay the default. The flooding
scenarios predate radiation and are validated against their own behaviour;
arriving physics should not silently rewrite a validated result.

**What the hull's own numbers said about the guesses.** For the 60 × 16 × 4 m
barge, as multiples of displaced mass:

| | guessed | measured |
|---|---|---|
| heave | 1.10 | **1.885** |
| sway | 0.90 | **0.190** |

Both wrong, in opposite directions, and both explicable: a wide shallow section
drags a great deal of water vertically and almost none edgewise, where it is
nearly a thin plate.

**Two errors were made wiring this up and both are worth recording.**

*The fit window.* `attachRadiation()` deliberately takes **no timestep**. The
first version passed the simulation `dt` straight into the state-space fit, which
sampled K(t) over 1.3 s of a 20 s decay. A fit that sees only the first 6% of the
retardation function sees something nearly constant, places its poles near zero,
and produces a model that *integrates* velocity instead of damping it — the ship
reached NaN in five steps. The fit is now sampled from `memoryDecayTime()`, which
is a property of the hull. `RadiationForce::step()` is an exact zero-order hold
and genuinely does not care what dt it is driven at, which is what lets the two
be separated.

*Damping counted twice.* The modal damping coefficients **were** radiation — a
lumped stand-in for it. Adding the real thing alongside them cost 27% of
mid-frequency heave. With radiation attached, the sway, heave and pitch modal
dampers are switched off. Roll keeps its stand-in *unless Ikeda is also attached*,
which is now the case that closes the loop: roll radiation damping is genuinely
small and the mechanism that matters is viscous, so deleting roll's stand-in
before Ikeda landed would have left the mode that most needs damping with almost
none. Surge keeps its damper too, since strip theory contributes no surge
radiation at all. Quadratic drag is untouched: it is viscous and separate.

That last point generalises. `A_inf[0][0]` is a **structural zero, not a measured
zero** — a strip has no longitudinal radiation problem, and surge added mass is
entirely a three-dimensional end effect. Taking `A_inf` wholesale would have
deleted a real term while looking like an improvement.

**Validation.** The load-bearing check is that the time-domain model agrees with
the frequency-domain table it was built from. Cummins splits radiation into an
instantaneous `A_inf` plus a memory convolution, and the memory carries the whole
frequency dependence — so a free decay settling at ω_d must behave as though its
added mass were `A(ω_d)`, *not* `A_inf`. Measured: the decay implies 1.170 × M
against the table's 1.258 × M at that frequency, 6.9% apart, where `A_inf` is
1.885 × M and would be 61% away. The test carries that gap as an explicit guard,
because a 15% tolerance on the wrong quantity would otherwise look like agreement.

Attaching radiation also produces a heave resonance the ship did not have before,
and leaves the Froude–Krylov sinc zeros exactly where they were — correct, since
radiation changes the response and not the excitation.

**Cost.** Building the table is 2.1 s at 9 stations and 4.8 s at 21; the runtime
memory term is free — under a microsecond a tick for seven state-space models,
below the noise floor of a tick measurement. That is the intended shape: an
expensive coefficient table computed once, and a convolution cheap enough to run
at 100 Hz.

**The reference-point transfer — was a documented gap, now closed.** Strip theory
assembles `A_inf` about the body-frame origin, midship on the baseline, because
that is where `stripTheoryTable()` puts it. The inertia it is added to is about
the **centre of gravity**, and rotational added mass is no more origin-independent
than a moment of inertia is. `transferAddedMass()` moves a 6×6 by the congruence
`A' = TᵀAT` with

```
T = [ I  d~ ]     d~ a = d x a,   d = new reference point - old
    [ 0  I  ]
```

which is the parallel-axis theorem in 6×6 clothing — feed a rigid body's own mass
matrix in and the rotational block comes out as `I + m(d dᵀ − |d|² I)`. Written
out for roll:

```
A44' = A44 + 2 dz A24 + dz^2 A22 + dy^2 A33
```

**It is not a correction, it is the difference between two different numbers.** On
the 60 × 16 × 4 m barge with KG 5 m, `A_inf` roll falls from 5.70 × 10⁷ to
3.48 × 10⁷ kg m² — 39% — and on a 120 × 22 × 6 m ferry-like hull with KG 7.2 m
from 3.82 × 10⁸ to 1.10 × 10⁸, a factor of **3.5**. Against that ship's own dry
roll inertia of 7.0 × 10⁸ that is the difference between adding 54% and adding
16%.

*Establishing the sign, rather than assuming it.* The literature convention was
not trusted; the code's own was measured, using a section whose answer is
analytic. **Rotate a semicircle about the centre of its own circle and no point of
the contour moves normal to itself**, so the radiation potential is identically
zero and `a44 = a24 = 0` exactly, at every frequency. Put that centre on the
waterline and everything strip theory reports for the roll mode about the baseline
is *pure transfer*: measured at `a44 = T² a22` and `a24 = −T a22` to six figures,
with the raw sectional values coming back at 4 × 10⁻²⁶ and 6 × 10⁻¹³ against an
`a22` of 1.5 × 10⁴ — machine zero. Transferring back up by the draft has to
annihilate the roll block, and does. **The opposite sign does not give a small
error there; it gives `4 T² a22`**, so the check discriminates by a factor of
four rather than by a tolerance. Independently, the height at which sway and roll
decouple, `−A24/A22`, comes out at 5.7 m above the baseline on the barge and
10.8 m on the ferry-like hull — near the waterline in both cases, where a
section's added-mass roll centre belongs. The other sign would put it under
the keel.

A second, corroborating measurement: `A44` about the baseline swings by a factor
of 2.7 across the seakeeping band on the barge (1.03 × 10⁸ down to 3.8 × 10⁷),
while `A44` about the centre of gravity at 5 m — close to that decoupling height —
is flat to 4% over the same band. Almost all of the apparent frequency dependence
of roll added inertia about an arbitrary point is the sway added mass riding a
lever. That is a measurement, not a proof of sign on its own: it fixes the
*distance* to the decoupling point and the semicircle fixes which side of the
baseline it lies on.

**Two things the earlier note got wrong**, both found while doing this:

- *Pitch is not unaffected.* Zero surge added mass kills the **vertical** part of
  the pitch transfer, which is what that note was reaching for, but not the
  **longitudinal** part: `A55' = A55 + 2 dx A35 + dx² A33`, which is exactly
  `∫(x − dx)² a33 dx` and is the same identity strip theory assembles `A55` from
  in the first place. Measured on the 120 m hull with the cog 6 m aft of midship,
  `A55` moves 4.9% and `A66` 3.7%. Heave really is unaffected: the translational
  block is reference-point independent, full stop, which is why the heave
  validation above never noticed any of this.
- *The memory force had the same mismatch, in the damping rather than the
  inertia.* `RadiationForce::memoryForce()` returns its moment about the origin
  the table was assembled about; `integrateRigidBody` was applying it as a moment
  about the centre of gravity. The correction is `M_G = M_O − d × F`, the same one
  the propulsion forces already carried. Removing it again shifts the measured
  roll added inertia by 17%.

**What is still approximated.** Only the *diagonal* of the transferred `A_inf` is
used: the integrator inverts a 3×3 inertia and divides by a per-axis mass, so
there is nowhere for the sway–roll coupling `A_24` to go. That is a limitation of
this integrator's shape, not of the physics. Doing the transfer first shrinks the
term being dropped by **8×** on the barge (−3.96 × 10⁶ to −4.80 × 10⁵), because
the centre of gravity sits close to the height at which sway and roll decouple —
so the transfer improves the approximation twice over.

### Propulsion coupled into the ship — **implemented, opt-in**

`Ship::propulsion` holds an optional `Manoeuvring`. Ship owns the motion: the
manoeuvring state is overwritten from the rigid body every tick and only ever
read back as forces, so the two cannot disagree about where the ship is. The
horizontal-plane forces `X`, `Y` and the yaw moment `N` are applied in the body
frame; `X` and `Y` act at midship while the integrator takes moments about the
centre of gravity, so the offset contributes a yaw moment of its own — dropping
it would make the ship turn about the wrong point in a way that still looks like
turning.

**Stand-ins removed, for the third time in this file.** A drag coefficient of
0.10 on the bow's projected area was standing in for hull resistance, and 1.00 on
the side for cross-flow drag. The MMG polynomial computes both properly (`R_0'`,
`Y_v` and its higher terms), so with a manoeuvring set attached those go, and the
yaw damper with them since `N_r` is what *it* was standing in for. Heave keeps its
drag: nothing in a horizontal-plane model speaks to it. Conversely, surge added
mass now comes from `m_x'` — filling exactly the hole strip theory leaves, since
`A_inf[0][0]` is a structural zero.

**One overlap is not resolved and is not hidden.** With both radiation and
propulsion attached, sway is damped by radiation (wave-making) *and* by `Y_v`
(lifting and cross-flow). These are different physical mechanisms that dominate
in different frequency bands — radiation at wave frequencies, `Y_v` at
manoeuvring frequencies — so adding them double-counts in the middle. At zero
speed there is no overlap at all, because `evaluateHull()` returns only
resistance when the speed is zero. The proper treatment is a frequency crossover
between the seakeeping and manoeuvring models, which is a research topic rather
than an oversight to be tidied away.

**The encounter frequency is emergent, and this is the point of the whole
exercise.** Nothing imposes it. The surface is `A cos(kx − ωt + φ)` and the ship's
x is moving, so a hull under way meets waves at `|ω − kU|` purely because it
translates through a spatially varying field. Measured on the barge at three
frequencies, against `encounterFrequency()` at the speed the ship actually
achieved:

| ω | U, m/s | predicted ω_e | measured ω_e | error |
|---|---|---|---|---|
| 0.55 | 2.02 | 0.6123 | 0.6125 | 0.0% |
| 0.70 | 2.31 | 0.8152 | 0.8160 | 0.1% |
| 0.85 | 2.13 | 1.0068 | 1.0075 | 0.1% |

The response frequency there is measured *independently*, by scanning the
harmonic fit for the frequency leaving least unexplained; reading
`encounterOmega` back out of `measureRaoAt()` would only confirm that a formula
had been evaluated. A first attempt reported 3–6% error and that was entirely the
measurement: an exponential drift filter with a 10 s time constant sitting on a
9 s wave is a high-pass filter on the signal.

**Two things this turned up about running a ship at all.**

*Acceleration takes minutes.* A 3.9 × 10⁶ kg hull reaches 61% of its final speed
after 300 s and settles by about 1800 s. A test asking whether it had converged
at 400 s was asking the wrong question of correct behaviour.

*The hull is directionally unstable*, which is unsurprising with KVLCC2
derivatives — that hull is the standard unstable example. Left alone it runs dead
straight for about 1900 s and then departs into a *steady* turn: surge 1.0515,
sway 0.4969, both constant, heading rotating uniformly. `measureRaoAt()` therefore
steers, with proportional-on-heading and derivative-on-yaw-rate gains in
`RaoSettings`. A ship that has quietly turned out of head seas is no longer
measuring the RAO that was asked for, while still producing a number.

That departure was first read as chaos, because the speed trace was oscillating
between +1.6 and −1.2 m/s. It was not: `RigidState::velocity` is a **world**
vector, and its x component is the ship's speed only while the ship still points
along world x. The surge speed was constant the whole time. `rao.cpp` and the
tests now both take speed along the ship's own bow.

**Forward speed in `RaoSettings` is now measured, not asserted.**
`accelerateSeconds` runs the ship up in still water first, and the encounter
frequency is re-derived from the mean surge over the recording window, because
added resistance in waves means the speed under way is not the speed it settled
at in flat water. `forwardSpeed` survives only for hulls with no machinery
modelled — a towed model — and is ignored outright when propulsion is attached.

One consequence worth stating: the reference wave is sampled **at the ship**, not
at the origin. A moving hull meets the wave at ω_e while the surface at a fixed
point still oscillates at ω, so fitting the origin's history at ω_e would return
a small meaningless amplitude and every RAO would be normalised by it.

### Viscous roll damping — **implemented**

Potential flow gives no roll damping at all — it is entirely viscous and it is the
difference between a ship that rolls 8° and one that rolls 35°.

`engine/sim/roll_damping.{hpp,cpp}` implements Ikeda's method and returns the
equivalent linear coefficient B44 (N·m·s/rad) together with its components, so
which mechanism is doing the work is always visible. Inputs are hull form
(Lpp, B, d, Cb, Cm, roll-axis height above keel, bilge radius), bilge keel
dimensions, and the operating point — roll amplitude, roll frequency, forward
speed.

**Empirical basis, component by component.** All of it is fitted to model
experiments; none of it is derivable.

| Component | Source | Behaviour |
|---|---|---|
| Friction | Kato (1958) + Tamiya's speed correction, via ITTC 7.5-02-07-04.5 (2.13)–(2.17) | Amplitude *cancels exactly*; ∝ √ω; the only component with a Reynolds number in it |
| Eddy | Simplified Ikeda regression, Kawahara, Maekawa & Ikeda, STAB 2009 eq. (31), fitted to Ikeda's own sectional method | ∝ ω·φ_a exactly (a linearised quadratic moment); collapses with forward speed per ITTC (2.21) |
| Lift | Ikeda et al. (1978), ITTC (2.10)–(2.11) | Exactly zero at zero speed, exactly linear in speed, amplitude independent |
| Bilge keel | Ikeda's sectional pressure model, ITTC (2.24)–(2.32), split into normal-force and hull-pressure parts | Grows with amplitude, but to a **finite non-zero limit** as φ_a → 0, because C_D rises as the Keulegan–Carpenter number falls |
| Wave (radiation) | **Not computed** — an input coefficient | Belongs to the BEM pipeline above; 5–30% of the total per ITTC |

Nondimensionalisation follows ITTC (2.4): `B44hat = B44 / (ρ ∇ B²) · √(B/2g)`.
Note the *multiplication* by √(B/2g); dividing by it leaves a quantity with units
of s⁻² and destroys the Froude scale invariance that makes B44hat useful.

**Measured on the reference hull** (170 × 25 × 6.5 m ro-pax, Cb 0.55, Cm 0.98,
KG 11.0 m, 34 × 0.6 m bilge keels; T_roll 13.9 s at GM 2.0 m): viscous B44hat is
0.0081 bare and 0.0439 with keels at ω̂ = 1.1 and 20°, i.e. 6.3% of critical
damping at the natural roll frequency and 10° — which is where `Ship::zetaRoll`
was guessed at 0.08. Bilge keels are 80% of the total, the eddy component
essentially all of the rest at zero speed, and friction 1.6% at full scale
against 8.2% for a 2 m model of the same hull (ITTC quotes 1–3% and 8–10%).

**Validity limits — where this stops being trustworthy.** `validateRollDamping()`
reports the first group; the rest cannot be detected from the parameters at all
and are recorded here because that is the only place they can live.

1. *Regression domain* (Kawahara et al., their own stated limits): 0.5 ≤ Cb ≤
   0.85, 2.5 ≤ B/d ≤ 4.5, −1.5 ≤ OG/d ≤ 0.2, 0.9 ≤ Cm ≤ 0.99, and for keels
   0.01 ≤ b_BK/B ≤ 0.06, 0.05 ≤ l_BK/Lpp ≤ 0.4.
2. *Hull form.* The fits come from a methodical series of conventional
   displacement cargo hulls. Kawahara et al. state explicitly that accuracy
   "might decrease remarkably" for buttock-flow sterns — which is exactly the
   modern large passenger ship and pure car carrier. Hard-chine hulls need the
   separate chine formulation (ITTC §2.3); barges, planing hulls and multihulls
   each have their own. None of those are implemented.
3. *Section idealisation for bilge keels.* Vertical side, horizontal bottom,
   quarter-circle bilge, keel at the 45° point and normal to the shell. A slender
   high-speed hull with a small bilge radius gets the moment levers wrong.
4. *High centre of gravity.* OG/d below −1.5 is extrapolation, and a high roll
   axis was the specific case that broke the earlier version of the simple
   formula. A top-heavy ro-pax at −0.7 is inside the range but near the edge of
   the data.
5. *Appendages other than bilge keels are absent*: no skeg, rudder, stabiliser
   fin, shaft bracket or bossing. ITTC gives a skeg formulation (§2.2.5.2) that
   is not implemented.
6. *Equivalent linearisation.* B44 is only valid at the amplitude, frequency and
   speed it was asked for; it must be re-evaluated as the operating point moves.
   Very large rolls where the bilge keel emerges from the water need Bassler's
   piecewise treatment and are not modelled.
7. *Scale.* The friction component alone has a Reynolds effect, so a B44hat
   measured on a model must not be reused at full scale without correcting it.
8. Deep, calm, unrestricted water. No shallow-water correction; no forward-speed
   effect on the bilge keel component (ITTC judges it small and it is neglected).

**Next**: per-ship tuning against roll decay tests — a per-component scale factor
identified from a measured decay curve, which is how this is done in practice and
what turns a ±25% method into a ±10% one for a specific ship. Then the skeg and
rudder components.

### Viscous roll damping coupled into the ship — **implemented, opt-in**

`Ship::attachRollDamping(waterlineZ, bilgeKeelLength, bilgeKeelBreadth)` reads the
hull form off the mesh and puts B44 into the rigid-body roll equation **in place
of** `Ship::zetaRoll`.

**Stand-ins removed, for the fourth time in this file**, and this one had a trap
in it: `zetaRoll = 0.08` was picked to sit where Ikeda puts a ro-pax with bilge
keels — 6.3% of critical, measured above — so leaving both in would not have
looked wrong anywhere. It would simply have doubled the damping of the mode that
decides whether a damaged ship capsizes. The test that catches it is the log
decrement, and it catches it by a factor of four.

Opt-in, like radiation and propulsion: the flooding scenarios predate all of it
and attach neither this nor radiation, so they keep the scalar coefficients.

They are not, however, *untouched* — and the difference between "the outcome line
is the same" and "nothing changed" is worth being precise about. The stiffness
correction below applies to every ship, so their roll damping moved by
`sqrt(KM/GM)` ≈ 1.8×. All three outcome lines are identical; the traces are not
bit-identical, and comparing them properly rather than grepping the verdict is
what says whether that matters. On the `doors` capsize the two runs track to
**0.1° of heel out of 52°** at t = 780 s and 7 t of floodwater out of 4894 — and
the numbers that decide survival, the damaged GZ curve and the per-compartment
final state, agree to the last printed digit. A slow loll driven by flooding is
not sensitive to how hard roll is damped, which is the reassuring answer rather
than the assumed one.

**Where the hull form comes from — derived, not authored.**
`rollDampingHullFromMesh()` reads Lpp, beam and draft off the wetted body, Cb from
the volume the hydrostatics already integrate, and Cm from the largest sectional
area found by clipping the hull into slabs — the same route
`radiationHullFromMesh()` takes to its stations, and for the same reason: a hull
that displaces what it should then also has a block coefficient that says so, and
there is no second set of main dimensions to drift out of step with the shape.

Three inputs are deliberately *not* derived, and the boundary is the point:

| Input | Where it comes from | Why |
|---|---|---|
| Lpp, B, d, Cb, Cm | the mesh | they are properties of the shape, and the shape is the mesh |
| bilge keels | the caller | a watertight envelope has no appendages in it — there is nothing to measure |
| bilge radius | left at the "estimate from Cm and B/d" sentinel | Ikeda's keel model idealises the section as vertical side, flat bottom, quarter-circle bilge; a radius measured off a shape that is not that is not the radius the formulae want |
| roll axis | the live centre of gravity, every tick | it is *loading*, not form, and a flooding ship's loading moves |

**The operating point — evaluated every tick, and why that is the cheap answer.**
B44 is an equivalent *linear* coefficient, valid only at the amplitude, frequency
and speed it was asked for, all three of which move. Measured, one `rollDamping()`
call is **52 ns against a 16 µs tick — 0.33%**, so caching it buys nothing but a
stale answer and the question settles itself. The three inputs come from
quantities the integrator already holds:

- **Frequency: the natural roll frequency** `sqrt(C44 / (I44 + A44))`, from the
  measured hydrostatic stiffness and the effective inertia. Roll is lightly damped
  and sharply resonant, so a ship's roll sits at its own frequency whatever drives
  it; and Ikeda's frequency dependence is only `sqrt(omega)` for friction and
  linear for everything else, so a frequency wrong by a fifth is a B44 wrong by
  less than the method's own ±25%. **Measured against the frequency the ship
  actually rolled at in a free decay: 0.098% apart.** An estimator fitted to
  recent motion would need filter state, would lag every change of loading, and
  would have nothing to say at all about a ship that is not currently rolling.
- **Amplitude: the phase-plane radius** `sqrt(phi² + (phidot/omega_n)²)`, which for
  a lightly damped oscillator *is* the envelope of the roll being executed —
  exactly the amplitude equivalent linearisation is defined at, instantaneous, and
  with no history to keep. Using the heel alone instead reports zero damping every
  time the ship passes through upright, which is where it is rolling fastest;
  substituted deliberately, it breaks the log decrement by 12%.
- **Speed: surge along the ship's own bow**, not `state.velocity.x`, which is a
  world vector — the mistake recorded further up this file.

**`waveDamping` stays zero when radiation is attached.** B44W is the radiation
share of roll damping, and the Cummins memory convolution is already applying it;
supplying both is the same double count that cost 27% of mid-frequency heave.
`attachRollDamping()` therefore sets it to zero explicitly rather than by
omission, and `validateRollDamping()` still reports it as missing for the
no-radiation case, which is the honest thing for it to say.

**Validation.** The reference hull for this is a 120 × 22 m ferry-like form at 6 m
draft — Cb 0.742, Cm 0.910, B/d 3.67, KG 7.2 m (OG/d −0.20), GM 2.09 m, roll
period 11.94 s, bilge keels 40 × 1.0 m. It is not the box barge the rest of the
seakeeping suite uses, deliberately: Ikeda's regressions are fitted over
0.5 ≤ Cb ≤ 0.85 and 0.9 ≤ Cm ≤ 0.99 and **a box is 1.0 on both**, so a barge's
eddy coefficient is pure extrapolation and its default bilge radius is exactly
zero. The tests assert the hull is inside every published bound rather than
assuming it.

| Check | Measured |
|---|---|
| Log decrement against `2πζ/√(1−ζ²)` from B44, per cycle, 3° release | worst **0.25%** over 14 cycles |
| Same, at dt = 0.02 / 0.01 / 0.005 s | 0.24% / 0.25% / 0.25% — not a timestep artefact |
| Frequency Ikeda was asked at vs the frequency the ship rolled at | 0.098% |
| 8° release, 15 cycles: bilge keels / bare hull / no damping at all | 6.30° → 0.58° / 7.92° → 6.99° / 8.00° → **8.13°** |
| 20° release: decrement, first cycle to last | 0.316 → 0.139, monotone |
| Same decay on the `zetaRoll` stand-in | 0.4996 → 0.5043 against the closed form's 0.50427 |

Damping ratio at the natural roll frequency, straight from the coefficient:

| roll amplitude | ζ with bilge keels | ζ bare | bilge-keel share |
|---|---|---|---|
| 2.5° | 0.0259 | 0.00062 | 97.6% |
| 5° | 0.0323 | 0.00105 | 96.7% |
| 10° | 0.0470 | 0.00193 | 95.9% |
| 20° | 0.0810 | 0.00368 | 95.5% |

Between 2.6% and 8.1% of critical with keels, against 0.06–0.37% bare: a factor of
**22 to 42**, which is why the bilge-keel test needs no tolerance to be
convincing. The keels dominate more here than the 50–80% Kawahara reports,
because this hull's roll axis sits close to the waterline and the eddy component
nearly vanishes there; the ro-pax in `test_roll_damping.cpp`, with OG/d = −0.69,
is the fixture that checks the component split against ITTC's bands.

Every ingredient of the first row comes from somewhere else: the stiffness from
GM (a forced-heel GZ sweep, not the finite difference the damper is scaled by),
the effective inertia from the period of a decay with damping switched off
entirely, B44 from `roll_damping.cpp` called directly, and the decrement from the
simulation. The last two rows of the first table are the pair that matters — **a
linear damper run through the same integrator, the same hull and the same peak
finder holds its decrement constant to 0.9%**, so the factor of 2.3 the nonlinear
run shows over the same amplitudes is the damper and not the measurement.

The "no damping at all" column is the control, and it *grows* by 0.2% over
fourteen cycles rather than decaying. That is the ceiling on what the integrator
itself contributes, and it is what makes the bare-hull decay — 12% over the same
record — a measurement of Ikeda rather than of numerical dissipation.

**A defect this turned up, which every functional test had passed.** The measured
angular stiffnesses that scale the modal dampers were being finite-differenced by
rotating the ship about the **body origin** while taking the moment about the
**centre of gravity**. That leaves the weight with a moment arm, so the stiffness
came out as `ρg∇·KM` instead of `ρg∇·GM`: too large by `ρg∇·KG`. Against a
longitudinal metacentric height of a couple of hundred metres that is nothing, and
pitch never noticed — against a transverse GM of a couple of metres and a KG of
5–7 m it is larger than the quantity being measured. The inflation factor is
`sqrt(KM/GM)`: measured on a 120 × 18 m hull at GM 2.28 m and KG 5.0 m, a nominal
`zetaRoll = 0.08` was delivering a log decrement of 0.90825, i.e. **0.1436 of
critical** against the 0.0800 on the label, where `0.08·sqrt((2.28+5.0)/2.28)`
predicts 0.1430. The fix is to rotate about the centre of gravity, after which the
same run gives 0.50433 against the closed form's 0.50427 — four decimal places on
a quantity that was 80% out.

It shipped green because `zetaRoll` is a dimensionless fudge factor: an inflated
stiffness only means the damper is stronger than its label, and nothing checked
the label. It was found by **timing a free decay and comparing the logarithmic
decrement against the ζ the coefficient claims** — the same instrument that
validates Ikeda, pointed at the thing Ikeda replaces. It matters here beyond the
label, because the natural roll frequency Ikeda is evaluated at is built from that
same stiffness and would have been 2.1× too high.

**Limits.** Everything in the validity list above still applies, plus three that
belong to the coupling rather than to the method:

1. *The operating-point frequency is the natural one*, so a ship driven hard well
   off resonance — a following sea overtaking slowly — is evaluated at the wrong
   frequency by whatever the detuning is. Ikeda is at worst linear in ω, so the
   error is bounded by the detuning and is not sharp.
2. *The natural frequency itself uses `A_inf`, not `A(ω_n)`*, because that is what
   the mass matrix carries; the self-consistent value would need an iteration. On
   the barge the two added masses are 0.4% apart and it does not matter at all; on
   the ferry-like hull `A_inf` is 21% below `A(ω_d)`, which is dry inertia plus
   16% against dry plus 20% — a 1.7% error in ω and so in B44.
3. *Hull form is snapshotted at `attachRollDamping()`* and only the roll axis
   tracks the ship afterwards, so a hull that floods until its draft and Cb have
   visibly changed is being damped on its intact form.

All three are refinements with a clear route, not gaps in the coupling.

### Hull-to-hull contact — **implemented**

`engine/sim/collision.hpp` / `.cpp`. Everything downstream of an impact existed
before this — a structural mesh, buckling, a solid-shell element, failed panels
that become flooding openings — and everything upstream did too. Two `sim::Ship`s
passed through each other without noticing. This is the middle: whether two hulls
interfere, *where* and *how much*, and what force that puts on both of them. It is
the rigid-body half only; there is no fracture and no deformation here.

#### Detection: a penetration volume, not a contact manifold

The usual answer is a list of contact points with normals and depths. It is the
wrong one here for three reasons. A ship hull is **not convex** — bulbous bow,
transom tuck, flare above the knuckle — and a convex decomposition that gets one
of those wrong invents contact where there is none. The point *count* is a
tessellation artefact, so a penalty force summed over it scales with how finely
the shell was wound rather than with the ship. And the consumer is a structural
model that asks which panels failed: "42 MN over 18 m² centred here" is a load
case, "nine points with depths" is not.

So detection computes the **interpenetration solid** A ∩ B — volume, centroid,
principal extents, and the two surface patches that bound it — exactly, for
arbitrary closed meshes, using only machinery the engine already had.

The device is the signed-cone identity that already makes `integrate()` correct
for any apex, stated for the *indicator* rather than for the volume: for a closed,
outward-wound mesh and any apex `o`,

    1_B(x) = Σ_j sign_j · 1_{tet_j}(x),   tet_j = (o, a_j, b_j, c_j)

pointwise almost everywhere. That turns "is this point inside a general closed
mesh" into a sum over **convex** pieces, and a tetrahedron is an intersection of
four half-spaces — the one thing this engine's geometry is built to clip against.
Applied to a triangle T of ∂A it gives `area(T ∩ B) = Σ_j sign_j area(T ∩ tet_j)`,
and every term is a convex polygon clip: no welding, no cap stitching, no ear
clipping, no mesh boolean anywhere. Summing the divergence-theorem flux of those
polygons over ∂A ∩ B and ∂B ∩ A gives the overlap's volume, centroid and second
moments, because ∂(A ∩ B) is precisely those two pieces — and the same sweep hands
back each patch's area, centroid and area-weighted normal for free.

Cost is O(triangles of A × triangles of B), so both meshes are first cut to the
overlap region with `clipToBox` — legitimate because A ∩ B ⊆ box implies
A ∩ B = (A ∩ box) ∩ (B ∩ box) — and the box is tightened against the clipped
meshes' own bounds a couple of times. The box is then padded outward, which costs
nothing and buys something worth having: A ∩ B is then *strictly* inside it, so the
cut faces `clipToBox` introduces lie outside the overlap and contribute exactly
zero to both the volume and the patches. No provenance tracking is needed.

**Depth and patch come from the overlap's own second moments.** A uniform box of
side lengths (d₁, d₂, d₃) has covariance eigenvalues dᵢ²/12, so the
eigen-decomposition of the overlap's covariance inverts to an *equivalent box*:
exactly for a box, and sensibly for anything lens-shaped. The smallest extent is
the penetration depth, volume over it is the projected patch area, and the
corresponding eigenvector is a contact normal that owes nothing to any
tessellation. The normal actually reported is the better-conditioned one — the
area-weighted mean of the two patches' *surface* normals — with the principal axis
as fallback and cross-check.

**Measured.** Two boxes: volume, centroid, extents, both patch areas, both patch
centroids and both patch normals exact to 1e-9. A square and the same square at
45° intersect in the regular octagon of area 2(√2−1)s², to 1e-9. Two tessellated
spheres reach the analytic lens volume at **second order** — measured error ratios
3.73, 3.92 and 3.98 per halving of facet size, against the 4 that owes. A box
inside a hull returns the box's own volume and the box's whole surface as patch.
Detection costs **5.4 ms** on a bow-into-side contact between two 464-triangle
hulls and **11 ms** at 3272 triangles; two hulls that are not touching cost
nothing at all, because the bounding boxes reject first. That is not cheap against
a 100 Hz budget, but a collision lasts of order a second, and the O(n·m) sweep is
the obvious thing to index if it ever needs to be faster.

#### Two defects the closed forms found

Both shipped green on everything else and both are worth recording.

**`Mat3{}` is the identity, not zero.** `Mat3` carries a default member
initialiser of `{1,0,0,0,1,0,0,0,1}`, so a value-initialised accumulator starts at
I. The second-moment integral then reported every solid's second moment one too
large on the diagonal — for a 2 × 4 × 6 box, a 6% error in the shortest extent,
which reads as a slightly fat box rather than as a broken integral. Found by
asserting the closed-form covariance d²/12, not by any volume test: the volume was
right throughout.

**Coincident faces were counted twice.** ∂(A ∩ B) is (∂A inside B) plus (∂B inside
A), and those two sets are disjoint *except* where a face of one hull lies exactly
in a face of the other — which both sweeps then claim, doubling that face's flux.
It is not an exotic case: two boxes meeting face to face share four side planes,
and the overlap volume comes out **5/3 too large** with no other symptom. Two
squares at 45° share their end planes and come out 4/3 too large — larger than the
whole of one of the solids, which is how it was caught. Shared faces now take half
weight in *both* sweeps, which counts them once and keeps the routine symmetric in
its arguments; faces coincident with *opposing* normals are two hulls touching
from outside, bound no overlap, and take weight zero. A cheap standing guard was
added with it — the overlap can never exceed either solid — and it is what would
have caught this on the first run instead of on a closed form aimed at something
else.

#### Response: a compliant force, not an impulse

At 6 m/s closing, two 20,000 t ships are **not** an impulsive problem. Take a
contact stiffness that puts the peak force where ship-collision measurements put
it and the contact half-period is π√(m_red/K) ≈ 0.5 s — fifty to a hundred
simulation steps. An impulse collapses all of that into one step and hands the
structural model a velocity change, when what it needs is a force history: a rise
time, a peak, a patch that grows and moves, and an energy that accumulates. So the
primary model is a **penalty** contact,

    F = stiffness · volume · (1 + dissipation · approachRate),  clamped ≥ 0

with Coulomb friction, regularised, on the tangential component. Force
proportional to the overlap *volume* rather than to a depth is the hydroelastic
form and the one that is mesh-independent: no contact count appears in it, it is
continuous as the patch grows, and `stiffness` is an honest physical quantity — a
contact pressure per metre of penetration, in Pa/m. The dissipation term is
Hunt–Crossley, proportional to the penetration, so the force is zero at first
touch and at separation and never goes tensile, which a linear dashpot does at
both ends.

The resultant acts at the **overlap solid's centroid**, which for a pressure
proportional to local penetration is the centre of pressure — the depth-weighted
centroid of the patch, and not the patch's area centroid. On the box fixture the
two are 0.3 m apart; on a ship-length lever that is a moment error of the same
order as the moment the contact is trying to produce.

The impulsive solver is still there. It is the closed form the penalty model's
conservation is checked against, it is the right model for a low-speed nudge that
would otherwise be resolved in one step, and its energy figure
`(1 − e²)/2 · m_eff · u²` is the *target* the penalty force's work integral has to
reach. Both apply equal and opposite loads at one shared point, so both conserve
linear and angular momentum identically.

**Hunt–Crossley's small-dissipation law is asymptotic, and reading it as a
calibration is how a default gets chosen twice as elastic as intended.** The law
says e = 1 − (2/3)·dissipation·u. Measured on the flat-block fixture at 6 m/s the
slope is 3.92 against the predicted 4.00 at 0.005 s/m — and has fallen to 2.05 by
0.25 s/m, where the linear extrapolation says e = 0 and the truth is e = 0.49. The
default was set from the extrapolation first and is now set from the measurement:
**0.8 s/m**, delivering e = 0.20 at 6 m/s and absorbing 96% of the closing energy,
which is where ship-collision energy ratios sit. It is speed dependent by
construction — 0.65 at 1 m/s, 0.12 at 10 m/s — and that is the right direction,
since a gentle touch does less plastic damage than a ram.

#### Coupling: nothing here replaces a lumped stand-in

Three times in `ship.cpp` — radiation, Ikeda, the MMG polynomial — coupling a real
model in has meant *deleting* a fraction-of-critical damper that was secretly
doing the same job, because adding the real thing alongside its own stand-in damps
the ship twice. **Contact does not have that problem**, and the reason is worth
stating rather than leaving to be rediscovered: every damping, added-mass and drag
term in `integrateRigidBody()` is a *fluid* mechanism. Radiation, viscous roll,
cross-flow drag, hull resistance — all of them act on a ship that is nowhere near
another ship, none of them was ever a proxy for hull-to-hull contact, and none of
them switches off while contact is happening. `Ship` therefore loses nothing and
gains one thing: an `externalForce` / `externalMoment` accumulator in the world
frame, taken about the centre of gravity, applied after the damping loop and
**cleared** once consumed.

It is divided by the *effective* mass, added mass included, which is right — a
struck ship accelerates with the water it has to shove aside — and which means the
two ships in a collision do not conserve momentum between themselves alone. The
difference is in the water. The conservation laws are asserted against bare rigid
bodies, where they hold exactly.

The one double count that does exist is small and in the other direction: while
the hulls overlap, both ships claim buoyancy from the same displaced water. A ram
deep enough to matter — 10 m³ of overlap — credits 10 t of extra buoyancy against
20,000 t of ship, below the noise of everything else, and removing it would mean
subtracting the overlap from each hull's own hydrostatic integral for no
observable change.

#### What a ram looks like in numbers

Two S-175-like hulls, one at rest, the other striking her port side 30 m forward
of midship at 6 m/s, both floating free with hydrostatics, damping and drag:

| | |
|---|---|
| closing speed | 6.0 m/s |
| contact duration | 1.62 s (162 steps at dt = 10 ms) |
| peak normal force | 323 MN |
| penetration at peak | 0.398 m |
| projected patch at peak | 11.4 m² |
| mean contact pressure at peak | 28.3 MPa |
| energy absorbed | 233 MJ |
| patch location, struck ship's frame | x = +30.0 m, y = +9.3 m, z = +8.3 m |

The patch lands within 5 cm of the station aimed at, on the port side, inside the
hull's depth — which is the number the structural model consumes.

#### What this cannot represent

- **Two contact regions at once.** A bow and a quarter touching simultaneously
  produce one volume, one centroid and one force, applied *between* them, where
  nothing is touching. The volume is right and the moment is wrong. This is the one
  thing a contact-point list does better. `normalAgreement` and `patchSeparation`
  are published so the case announces itself: two regions a ship's length apart put
  the two patch centroids a ship's length apart, where a single region puts them
  within its own penetration depth.
- **Deep or engulfing penetration.** Once the overlap stops being a thin lens the
  equivalent box is no longer flat and its smallest extent stops being a depth. In
  a real ship that state arrives long after the shell has failed, but nothing here
  detects it.
- **The structure.** The contact is rigid. A real 6 m/s ram spends most of its
  energy crushing bow and side at roughly constant force over metres of
  penetration; a rigid contact reaches a much higher peak force over a much smaller
  penetration for the same absorbed energy. **Peak force is an over-estimate and
  penetration an under-estimate**; only the energy and the patch location survive
  the approximation intact.
- **Hydrodynamic interaction.** The water squeezed out from between two closing
  hulls resists, and the added mass of a hull alongside another is not that of a
  hull alone. Neither is modelled; both make a real collision slightly softer.
- **Coincident surfaces**, in one direction. Faces coincident with the *same*
  outward normal are exact — that is the fix above, and a shared face is split
  evenly between the two reported patches, which is a choice rather than a
  measurement. Faces coincident with *opposing* normals are exact whenever the
  coincident tetrahedron is the only one claiming them, which covers two convex
  hulls laid flush; on a non-convex solid a face on its own skin can also be
  claimed by a tetrahedron it is not coplanar with, and no per-face rule undoes
  that. Measured on a box wedged flush into an L-prism's notch while overlapping
  its other arm: 2.10 m³ reported against 1.60 m³. Deciding membership per point
  on a coincident face is the coincident-face problem every mesh boolean has and
  it is not solved here, so **that case is reported in `problems`** instead of
  being presented as a measurement. It needs exact coplanarity, to a part in 10⁹
  of the contact region, which two independently placed hulls do not produce.
- **`ContactHistory::work`** is a rectangle rule on the contact power at the start
  of each step, so it is **first order in dt** — 1.5% of the closing energy at
  dt = 1 ms, halving exactly with the step. It is a diagnostic; nothing in the
  dynamics reads it back.

### The hard cases

- **Green water on deck**: waves over the bulwark, water loose on deck. Handled by
  the same escalation as sloshing (§6): the deck becomes a shallow-water solver
  when wet.
- **Slamming**: bow flare and bottom slam impacts, Wagner/von Kármán added-mass
  impact theory to get the pressure pulse, delivered to the FEM as a load.
- **Parametric roll**: emerges naturally from nonlinear restoring in head seas
  once the above is in place. A good validation target precisely because it is
  hard to get by accident.
- **Broaching and surf-riding** in following seas.

---

## 3. Structure — adaptive FEM

### Scantlings → structural mesh — **implemented**

`engine/sim/scantlings.hpp` / `.cpp`. The geometry foundation: a ship described
the way a ship is actually specified — shell plating thickness by strake and
region, frame spacing and section, longitudinal spacing and section, deck and
bulkhead plating, girders — turned into plating panels and stiffener line
elements against a hull form.

**Two decisions everything downstream inherits.**

**Stiffeners are discrete line elements, not smeared into an equivalent plate.**
The smeared alternative (`t + A/s`, computed by `smearedThickness()` and kept for
the Tier-0 beam) is a third of the elements and exactly right for axial
stiffness. It is rejected because it cannot represent stiffened panel collapse:
collapse is the plate buckling in half-waves *between* stiffeners, and a
homogenised panel has no between. The missing thing is a length scale, not a
value, so no choice of `t_eq` recovers it. Measured cost of the smearing: for a
200 × 10 flat bar on 12 mm plating at 700 mm spacing, the panel's second moment
falls from 2.49 × 10⁻⁵ m⁴ to 1.91 × 10⁻⁷ — **a factor of 130** — because the
stiffener's area is moved into the plate, where its lever arm is zero. The two
panels have identical area and identical axial stiffness, which is exactly why an
axial check would not notice. What discrete costs: roughly 3× the elements, and every
stiffener is an eccentric beam whose offset from the plate mid-surface has to be
carried explicitly (`StructuralMember::rise`) or the section modulus is wrong by
that same Steiner term.

**The structural mesh is independent of the hydrodynamic hull mesh.** They are
refined for different reasons and neither refinement is negotiable: the hull
mesh's resolution comes from a *volume* convergence study (21 waterlines, 41
stations — a 3 m station spacing on a 120 m ship) and is paid for every tick by
the buoyancy integral; the structural mesh's comes from the frame spacing, which
is a property of the ship. Sharing forces either refining the hull mesh to frame
spacing and paying for it in the hydrostatic inner loop forever, or snapping the
frames to whatever stations the hydrostatics chose. So the hull is a *reference
surface* only, sampled by ray casting — which also means the same scantlings work
against a hull from `makeHullFromParticulars`, from a `.ship` offset table or
from an importer, with no shared topology. What that rules out: the two meshes
cannot deform together for free. When the FEM moves the structure the
hydrodynamic hull has to be updated by projection, an extra mapping with its own
error, where a shared vertex would have moved once.

Supporting conventions: panel corners lie on the plate's **mid-surface** (taken
as the hull surface itself — the 6 mm difference from shipbuilding moulded lines
is what a shell element wants); **girth** runs from the centreline keel outboard,
round the bilge and up to the deck edge, and regions and stiffener seams are
*fractions* of it, so a longitudinal is one continuous member while the girth it
lives on shrinks towards the ends.

**Measured, on the 120 m reference ferry** (`ferryScantlings()`, 2.4 m web
frames, 0.70 m longitudinals, four decks, six bulkheads, three girder lines):

| | |
|---|---|
| Generation | **3.0 ms** — once, at load |
| Panels | 8 900 (3 100 shell, 4 172 deck, 1 628 bulkhead) |
| Members | 11 846 (3 000 longitudinals, 3 162 frames, 3 772 deck longitudinals, 332 deck beams, 1 378 bulkhead stiffeners, 202 girder segments) |
| Memory | 2.46 MB (120 B/panel, 128 B/member) |
| `hullGirderSection()` | 0.09 ms per cut |
| At 0.6 m ordinary framing | 15.1 ms, 31 786 panels, 44 162 members, 9.03 MB |
| Steel | plating 1 084 t + stiffening 631 t = **1 716 t** |
| Midship section | A = 1.80 m², neutral axis 6.71 m above baseline, I = 46.2 m⁴, Z<sub>deck</sub> = 5.57 m³, Z<sub>keel</sub> = 6.88 m³ |

The section modulus is against an IACS unified-requirement minimum of 3.27 m³ for
L = 120 m, B = 20 m, C<sub>b</sub> = 0.66 (`ruleMinimumSectionModulus()`), so the
arrangement clears the rule by 1.7× at the deck. That is on the generous side and
the reason is the hull, not the scantlings: this ferry's offsets carry full
breadth to 15 m, giving L/D = 8 where a real ro-pax runs 10–12, so the girder is
unusually deep.

**Is 1 716 t credible for a 120 m ro-pax?** It is low, and knowably so. The
generator builds plating and primary stiffening and nothing else — no brackets,
no double-bottom floors, no pillars, no foundations, no stem or stern frame, no
weld metal, no minor structure. Those are conventionally another 15–30% of hull
steel, which puts the honest equivalent at 2 000–2 200 t against a real hull steel
weight of roughly 2 500–3 500 t for the type. The internal proportions are right
where they can be checked: 0.58 t of stiffening per tonne of plating against a
shipbuilding norm of 0.4–0.7, and a mean plating thickness of 11.4 mm. The
comparison against the ferry's own 8 984 t displacement at the design draft is
weaker than it looks, because that vessel's lightship is heavy for its size by
construction.

**What the representation cannot express.** Recorded because each is a decision
deferred, not an oversight:

- **Bulb flats** — the commonest longitudinal on a real ship. IACS handles them by
  converting to an equivalent angle with a published dimensional transformation,
  and a transformation reproduced from memory would be a plausible wrong number in
  the one place the whole section modulus depends on it. Add it with the rule text
  in hand.
- **Ordinary frames between web frames.** `frameSpacing` is one pitch. The ferry
  is described with 2.4 m webs and no 600 mm ordinary frames between them.
- **Longitudinals terminating towards the ends.** Seams are constant girth
  *fractions*, so longitudinals converge instead of being dropped in groups.
- **Horizontal stringers on bulkheads.** A transverse bulkhead's horizontal panel
  seams are a mesh subdivision with no member on them.
- **Stem, transom and stern-frame plating.** The shell runs keel to sheer only;
  the hull mesh's end caps carry no structure.
- **The product of inertia of an unsymmetrical section.** `secondMomentWeak` is
  about the geometric axis parallel to the web. For the hull girder this cancels
  identically wherever a web is vertical or horizontal and is confined to
  stiffeners on the turn of the bilge; it would matter for lateral-torsional
  buckling of an angle.
- **Asymmetric hulls.** The starboard shell is the mirror of the port, which the
  hull mesh format enforces anyway.
- **Node connectivity.** Panels and members are independent geometric records;
  merging coincident corners into shared FE nodes is the next step, not this one.
- **Authoring.** Scantlings are built in C++, as hulls were before
  `ships/ferry.ship`. A `structure` block in the ship file is Phase 7's.

**How it is checked** (`tests/test_scantlings.cpp`). Three instruments do the
work, and the first two each found a defect the third would not have:

- A **box hull**, where the shell is `L·B + 2·L·D` exactly and every panel is
  planar, so tiling and steel weight are equalities rather than estimates. This
  caught `hullGirderSection()` counting *both* bays either side of a frame
  station: area and second moment came out doubled while the neutral axis — a
  ratio — looked perfectly correct at every station. The cut interval is now
  half-open.
- **Quadrature over the section's own width profile**, which reaches the same
  second moment without ever forming a parallel-axis term.
- **Refinement**: the shell panels are flat chords across a curved surface, so
  they converge to the hull's own area *from below* — 2.1 × 10⁻³ at the reference
  spacings, halving to 5.7 × 10⁻⁴ by the third refinement — and section properties
  at a tapered station converge under frame refinement.

**Mutation testing killed 36 of 37 mutants**, and the eight it did not kill on
the first pass were all real holes in the suite: the hull girder taking a
stiffener's strong axis where its web lay flat; starboard webs not mirrored; the
frame count floored instead of rounded; panel area doubled from one triangle
(exact for a rectangle, wrong for every trapezoid); a later plating region failing
to override an earlier one; deck and bulkhead panels not clipped to the shell; a
girder built outside the hull; and longitudinals shifted onto the wrong strake
seam, which leaves one on the centreline keel while weight and count say nothing.
The single survivor is genuinely equivalent: a station-and-waterline hull mesh has
a single-valued half-breadth, so nearest and outermost hit are the same ray.

### Failed structure → flooding openings — **implemented**

`engine/sim/breach.{hpp,cpp}`. The half of the Phase 3 milestone that is not the
FEM: given a `StructuralMesh` and the panels that failed, produce the holes that
failure implies in the flooding network. It contains no fracture model and does
not pretend to — the failure set is the caller's answer, eventually the FEM's.

```cpp
BreachSet breachesFromFailedPanels(const Ship&, const StructuralMesh&,
                                   const std::vector<int>& failedPanels,
                                   const BreachParams& = {});
std::size_t applyBreaches(Ship&, const BreachSet&);   // straight into Ship::openings
int    spaceAt(const Ship&, const Vec3& bodyPoint);   // compartment / kSea / kEnclosedVoid
double meshWindingNumber(const TriMesh&, const Vec3&);
```

**Connectivity is read off the geometry, never off the label.** Which two spaces
a hole joins is the one thing that must not be wrong — a shell panel opens a
compartment to the sea, a bulkhead panel opens one compartment to another, and
confusing them is the difference between a breach and a door. `PanelRole` looks
like the answer and is not, and the reference ferry breaks it three ways: her
wing bulkhead runs at |y| = 6 m *inside* the holds, so tearing it opens a space
to itself; her weather deck is `PanelRole::Deck` with the sky above it, so a hole
there is open to the sea; and amidships her engine-room boxes stop at |y| = 8 m
while the shell is out past 8.6 m, so most of the plating over the engine room
faces nothing the subdivision describes. Each failed panel is therefore probed
either side of its own plane, and the probe returns three answers — a
compartment, `kSea`, or `kEnclosedVoid` for inside the hull but inside no
compartment. The third is not an error: it is a hole into a space the ship
definition does not model, and it opens nothing and says so.

Point location is a **summed solid angle** (Van Oosterom–Strackee), not a ray
cast. Compartments are carved from the hull on axis-aligned bulkhead and deck
planes, so their edges lie along exactly the axes a ray would be written along,
and a ray through an edge is counted twice or not at all. The solid angle sum has
no ray to place. On a *surface* its answer is decided by the sign of a
floating-point zero — a coplanar triangle cancels the triple product to ±0 and
`atan2` reads ±π off it, giving 1 on the +x, +y and +z faces of a box and 0 on
the other three — which is why probes are kept off surfaces and why the overlap
check below asks twice.

**The probe marches.** A flat panel chords across a curved shell, so its centroid
is not on the surface it stands for: **measured at 0.15 m** on this ferry, where
a girth band spans the crease at the turn of the bilge. A single fixed probe
cannot serve — one small enough to stay inside a 1.8 m double bottom is smaller
than that and reads the same space on both sides, which drops the breach. So the
probe doubles outward until the two sides disagree and stops at the first
disagreement, which is the nearest boundary and therefore this panel's own
surface, giving up at the distance from the centroid to the furthest corner.
Measured stability against the starting step, over all 8 900 panels of the ferry:

| first step | openings | torn area |
|---|---|---|
| 20 mm | 55 | 10 136 m² |
| **50 mm** (default) | **61** | **10 179 m²** |
| 100 mm | 65 | 10 182 m² |

**A torn plate is one hole, not forty.** Failed panels that share an *edge* and
join the same pair of spaces merge into one opening whose area is the sum of
theirs and whose position is the region's area-weighted centroid. Both qualifiers
carry weight. *Same pair*: a shell panel and the bulkhead panel it lands against
share an edge and are not the same hole, because an `Opening` has exactly two
ends. *An edge, not a corner*: two panels touching only at a corner stay two
openings, because the pinch between them has zero width and the merged centroid
would sit at the one point in the region where there is no hole. It costs nothing
in total flow — `Cd·A·√(2Δp/ρ)` is linear in area — and it keeps the head right
where the pieces sit at different depths, which is the case where splitting is
the better quadrature of `∫√h dA` and merging is the approximation.

**Discharge coefficient: 0.60.** A clean sharp-edged orifice runs 0.61–0.62 (the
ferry's authored breach uses 0.62) and a rounded or ducted entry 0.8–0.95, so the
only question is which side of the sharp-edged value a tear sits. It sits at or
below it: fractured plating petals, and an edge protruding into the flow is
re-entrant, whose contraction coefficient is 0.5 in the Borda limit.
Damage-stability practice takes 0.6 for a damage opening. 0.60 is therefore the
top of the honest range rather than the middle, so the error is towards flooding
too slowly. It is a constant because nothing upstream can refine it: what would
is the tear's aspect ratio and which way the plating folded, and a set of failed
panels records neither.

**Does the water follow the area?** That is the milestone's actual claim, so it
is measured rather than asserted. A barge with 200 × 40 m of waterplane, a hold
20 m long and a hole 2 m under, against `Cd·A·√(2gh)·t` over 30 s:

| failed plating | water taken | orifice law | ratio |
|---|---|---|---|
| 0.5 m² | 56.416 m³ | 56.368 m³ | 1.00084 |
| 1.0 m² | 112.926 m³ | 112.736 m³ | 1.00168 |
| 2.0 m² | 226.231 m³ | 225.472 m³ | 1.00337 |

**Doubling the failed plating multiplies the water by 2.0017.** The excess is not
noise and is not a discrepancy: the wider hole sinks the ship faster and so
raises its own head, and the residual scales exactly linearly with area, as a
first-order sinkage term must. Handing the flooding solver an opening this file
produced and one written out by hand with the same numbers gives water volumes
that agree to the last bit.

**What it found in the ship it was pointed at.** The probe reports the reference
ferry's forward wing tanks as lying *entirely inside* her forward holds —
`fwd_hold_p` is authored y = 0…20 m where `wing_tank_fwd_p` is 8…20 m over the
same length and height — so 217 m³ of her is described twice and a tear out there
is attributed to whichever compartment was declared first. `Ship::validate()`
does not catch it, because the subdivision still totals 25 259 m³ against a
28 273 m³ hull and its overlap test only fires past 100.1%. The aft wing tanks,
authored on the same plan, are correct. Reported rather than fixed: the flooding
scenarios are validated against the ship as it stands.

Distinguishing a real overlap from two compartments that merely *abut* needs care
for the signed-zero reason above — both claim a point on the bulkhead between
them, and the ferry's bulkhead deck lays a panel whose probe lands exactly on the
forepeak boundary at x = 44 m. So the claim is only made when it survives a
millimetre's displacement in all six directions: a region with volume in it does,
a face does not. The check therefore fails towards missing a real overlap rather
than towards inventing one.

**What this cannot yet express:**

- **Which panels fail.** The whole fracture half. Everything here is consequence.
- **Partial failure.** A panel is torn or intact; there is no fraction of a plate
  gone, and no separate treatment of a stiffener failing without its plating.
- **A tear's shape.** The opening is an area at a point. A long vertical slit and
  a compact hole of the same area are the same orifice here, though the slit's
  flow should be the height integral of `√h` rather than `√h` at the centroid —
  which is exactly the error the merge introduces and the reason corner-touching
  pieces are deliberately left separate.
- **Deformation.** The panels are read where the scantlings put them. A crushed
  bow's plating has moved, and neither the structural mesh nor the compartment
  meshes follow it yet.
- **Air paths.** A tear vents air as readily as it admits water and the network
  already handles that, but structure between compartments that fails *without*
  becoming a flooding path — a buckled but unsplit bulkhead — has no
  representation.
- **Progressive failure.** One call, one failure set. A hole that grows is a new
  call, and nothing dedupes it against the openings already added.

**How it is checked** (`tests/test_breach.cpp`). A rectangular barge where every
area, centroid and head is an exact rational; the reference ferry, whose hull is
curved and whose compartments do not reach the shell everywhere; and an authored
twin flooded alongside a produced opening and required to agree bit for bit.
**Mutation testing killed 47 of 49 mutants.** Five of the survivors on the first
pass were real holes and are now closed: the marching probe was never exercised
at all (no test had a panel whose centroid was off its own surface); a collapsed
quad side was treated as an edge, which fuses two triangles meeting at a point;
the weld's neighbourhood search could be reduced to a single cell; an inside-out
compartment mesh read as empty sea; and an out-of-range panel index survived
because reading past the end of the array produced a *different* complaint that
kept the problem count right. The two remaining survivors are equivalent: the
bounding box is a pure pre-filter that the winding number overrules, and taking
the probe's reach from one corner instead of the furthest agrees on every quad
whose corners are equidistant from its centroid, which is every one the generator
makes.

### Adaptive tetrahedral FEM

### Tier-0: the hull girder — **implemented**

`engine/sim/girder.{hpp,cpp}`. The ship as a free beam, which is the cheapest
structural answer worth having and the one naval architecture has computed by
hand for a century. Weight and buoyancy cancel in total and in moment — she
floats, and floats level — but not station by station, and what is left over
bends the hull.

```
q(x) = w(x) − b(x)     V(x) = ∫q     M(x) = ∫V      hogging positive
```

Weight *minus* buoyancy, and that ordering is the whole sign convention: with the
other one a hogging ship comes out negative. Getting it backwards produces a
plausible curve of the right magnitude that names every failure as its opposite.

**The free ends are the instrument.** Nothing holds a floating ship up at the
perpendiculars, so V and M must both be zero there — guaranteed by the two
equilibrium conditions and by nothing else. A residual is therefore a direct
measure of every error upstream. Measured on the ferry: shear closes to 4 × 10⁻⁵
of W/2 and moment to 1.8 × 10⁻³ of WL/8.

**Balancing on the wave is not optional.** The first version computed the girder
with the ship left at her still-water attitude, and the result was not merely
inaccurate — the shear and moment curves grew monotonically to the forward
perpendicular and never closed, because an unbalanced ship carries a net force
and a net moment that swamp the wave-induced bending. `balanceOnWave()` solves
for the sinkage and trim that make buoyancy equal weight *and* put B under G:
two residuals, two unknowns, Newton with a numerical Jacobian, one whole-hull
integral per residual rather than one per station.

**Measured on the 120 m ferry**, 8984 t, under a standard wave of her own length
and L/20 height:

| condition | bending moment | |
|---|---|---|
| still water | +1.57 × 10⁸ N·m | hogging |
| crest amidships | +4.70 × 10⁸ N·m | hogging |
| trough amidships | −1.65 × 10⁸ N·m | sagging |

The wave-induced component is therefore about ±3.2 × 10⁸ N·m, against an
IACS-style rule-of-thumb sagging moment `0.11 C L² B (Cb + 0.7)` of
3.59 × 10⁸ N·m — **agreement within 12% with an independent industry estimate**.
That is a sanity check rather than a validation: the rule is written for a
specific design wave and load case, so what it establishes is the right order of
magnitude arrived at down a completely different road.

**From moment to stress.** `girderStress()` divides the moment curve by the
section modulus from `hullGirderSection()` — the whole of classical longitudinal
strength, `sigma = M / Z`. The sign carries more information than the magnitude:
hogging stretches the deck and compresses the keel, sagging reverses both, and
deck plating *in compression* buckles well below yield. A stress magnitude alone
cannot see that failure.

Measured on the ferry, with the reference scantlings (1716 t of steel, midship
I = 46.2 m⁴, Z_deck = 5.57 m³ against an IACS minimum of 3.26 m³):

| condition | M, N·m | deck | keel | utilisation |
|---|---|---|---|---|
| still water | +9.35 × 10⁷ | +28.9 MPa | −18.6 MPa | 0.12 |
| crest amidships | +4.70 × 10⁸ | +84.4 MPa | −68.4 MPa | **0.36** |
| trough amidships | −1.65 × 10⁸ | −29.6 MPa | +24.0 MPa | 0.13 |

Hogging dominates because she already hogs in still water, so the wave crest adds
to an existing moment while the trough subtracts from it.

**That 0.36 is the cross-check worth having.** This hull carries 1.7× the rule
minimum section modulus, so a rule-minimum ship in the same condition would run at
0.36 × 1.7 = 0.61 of yield — squarely inside the 0.6–0.75 band real ships are
designed to. The moment and the section were computed by completely separate
routes, neither aimed at that number.

### Progressive collapse — **implemented**

`engine/sim/collapse.{hpp,cpp}`. First yield is not strength. A hull girder does
not fail when its worst panel reaches its limit: that panel sheds load onto its
neighbours, the neutral axis migrates towards the side still carrying, and the
section goes on taking moment until enough of it has gone that the total starts
to fall. **Smith's method** — impose a curvature, let every element answer with
its own load-shortening curve, move the neutral axis until axial forces balance,
sum moments, step and repeat — is what classification societies use, and the peak
of the resulting curve is the ultimate strength.

It is testable because it sits between two exact answers. At zero curvature the
slope must be `E·Σ(A d²)` about the elastic neutral axis. With buckling switched
off, a large curvature must give the fully plastic moment about the **plastic**
neutral axis — the one that balances *area*, not first moment, which on an
asymmetric section is a different height and using the elastic one understates
the answer.

The elements' own second moments are deliberately absent from that slope. Smith's
method carries axial stress alone — a strip of plating resists by stretching, not
by bending about its own mid-thickness — so the initial stiffness is slightly
below `E·I` as `hullGirderSection()` reports it. Measured on the ferry: 90–100%
of it, the gap being the plating's own inertia.

**What it says about the ferry** (369 elements at midship, A = 1.80 m², I = 46.2 m⁴):

| | moment | as a fraction of first yield |
|---|---|---|
| first yield (deck governs both ways) | 1.978 × 10⁹ N·m | 1.00 |
| fully plastic | 2.573 × 10⁹ N·m | 1.30 |
| **ultimate hogging** | 1.987 × 10⁹ N·m | **1.00** |
| **ultimate sagging** | 1.288 × 10⁹ N·m | **0.65** |

**In sagging she fails at 65% of first yield.** The deck buckles in compression
long before any fibre reaches yield stress, so a first-yield check overstates her
sagging strength by a factor of **1.54** — it would clear a moment half again
larger than the one that actually collapses her. Hogging shows no such gap,
because it compresses the keel, which is heavy, closely framed structure that
yields rather than buckles.

That asymmetry — sagging ultimate at 65% of hogging ultimate — is the standard
finding for a ship with a light, widely stiffened deck, and it is emergent here:
nothing in the method knows which fibre is which, only which one is in
compression.

Against the standard wave of her own length the margins are 4.2× hogging and 7.8×
sagging, which is comfortable and consistent with her carrying 1.7× the rule
minimum section modulus.

**Where does she break?** `longitudinalStrength()` runs the collapse calculation at
every station of a bending moment curve and reports the ratio, because "how
strong is midship" and "where does she fail" are different questions and the
ferry separates them.

Her scantlings already taper — `side_forward` and `side_aft` drop the side plating
from 12 mm to 10 mm beyond ±44 m, and the element count falls at the ends as decks
and bulkheads run out. Measured, her ultimate moment runs from 1.18 × 10⁹ N·m
amidships down to 6.8 × 10⁸ at the ends, **a factor of 1.74**. So her weakest
sections are her ends. Her worst *margin* is nonetheless amidships, because the
ends carry almost no bending moment. A calculation reporting either number alone
would be answering the other question.

Sweeping the whole ship costs 142 ms hogging and 239 ms sagging over 41 stations.
Getting there needed one fix: the neutral-axis bisection ran a fixed 200
iterations, and a double carries 53 bits, so a bracket nine section-depths wide is
exhausted in about sixty halvings. Two thirds of the solve was refining digits
that do not exist. Stopping on bracket width instead is **4× faster** — 70.7 ms to
17.3 ms for a 300-step curve — with identical answers to every figure.

**The load-shortening curve is the model, and it is one number.** Past the
compressive cap, `σ = σ_c (ε_c/ε)ⁿ`, continuous at the cap by construction, with
`n ≈ 0.45` for plating and zero for stiffeners. The published curves are per
failure mode and far more elaborate; this family is chosen so the *shape* is
right and the assumption is visible rather than buried. `n = 0` gives a perfect
plateau and an upper bound on strength — the curve then never turns over, and
`ultimateMoment` is the plateau rather than a peak.

### Buckling — **implemented**

`engine/sim/buckling.{hpp,cpp}`. Dividing a moment by a section modulus is the
right answer for the fibre in *tension*. For the one in compression it is not even
the right question: a plate panel under edge compression folds at a stress that
can be a small fraction of yield.

Two modes, checked separately because they fail at different stresses and are
cured by different things — **plate buckling** between stiffeners, going as
`k π²E/(12(1−ν²)) · (t/b)²`, and **stiffener column buckling** between frames, as
`π²EI/(AL²)`. The buckling coefficient is a minimum of `(m/α + α/m)²` over integer
half-wave counts: exactly 4 at every whole aspect ratio, peaking at 4.5 at the
crossovers, tending to 4 for a long panel. That is why "4" is the number quoted
for ship plating, and why a square panel is no weaker than a long one.

Elastic buckling stress is not strength — above about half yield the material goes
plastic before the instability arrives, and the elastic formula runs away to
values the plate cannot reach. The Johnson–Ostenfeld cap is continuous at the
transition **by construction**: at `σ_cr = σ_y/2` the lower branch returns its
input and the upper returns `σ_y(1 − σ_y/(4·σ_y/2)) = σ_y/2`, the same number
exactly.

**What it changes.** The ferry, same conditions as the table above but against
AH36's 355 MPa:

| condition | yield utilisation | compressed fibre | buckling utilisation |
|---|---|---|---|
| still water | 0.08 | keel | 0.09 |
| crest, hogging | 0.24 | keel | 0.26 |
| trough, sagging | **0.08** | deck | **0.30** |

Hogging barely moves — the keel is thick, closely framed structure. Sagging moves
by a factor of **3.7**: a condition that reads as 8% of yield is at 30% of what
the deck can actually carry, because the deck is the thinner, more widely
stiffened structure and it is furthest from the neutral axis. A yield check alone
would rank sagging as the benign case. It is the dangerous one.

**Limits, and they are not small.** This is elastic buckling of an ideal flat
panel: no initial distortion, no welding residual stress, no lateral pressure, no
interaction between the two modes, simply-supported edges throughout. A real rule
check applies knock-down factors for every one of those. What is here answers "is
this panel in the dangerous region", which is a question the hull girder result
could not ask at all.

**What this tier still cannot do.** It is a beam: it knows nothing about where
stress goes *within* a section, nothing about shear lag, nothing about local
loads, and nothing about what happens *after* a panel buckles — post-buckling
strength, load shedding to the stiffeners, and the progressive collapse that
follows are what the FEM tiers are for.

The lightship weight curve is a trapezoidal fit matched to total weight and LCG —
the Prohaska/Biles construction, and an *assumption*. Floodwater is not assumed:
each compartment's water is spread over that compartment's own extent, because
the flooding solver already knows where it is. A flooded hold amidships sags the
hull, and that coupling is asserted.

The chosen fidelity target is full 3D tetrahedral finite elements. The physical
obstacle is arithmetic: a 300 m hull with 20 mm plating needs ~10 mm elements to
resolve bending through the plate thickness, which is on the order of 10¹¹
elements for the whole ship. That is not a hardware problem to be solved by
waiting.

The resolution is **spatial adaptivity with substructuring**, which is how crash
and blast FEA has always handled the same conflict:

### Three-tier model

**Tier 0 — global girder, always on.** The hull as a Timoshenko beam with
section properties derived from the real structural arrangement. ~200 DOF.
Gives hull girder bending, hogging/sagging in waves, shear, torsion. Runs at
100 Hz for free. Its input exists: `hullGirderSection()` returns area, neutral
axis, second moment and section moduli at any cut, in 0.09 ms, taking membership
from geometry rather than from element labels so nothing athwartships can leak
into it.

**Tier 1 — reduced 3D, always on.** The full 3D shell/solid model of the entire
ship, condensed offline by **Craig–Bampton component mode synthesis**: retain the
interface DOF plus the lowest few hundred fixed-interface normal modes, discard
the rest. Reproduces the linear-elastic response of the full model to within a
few percent for the frequency range that matters, at 10⁻⁵ of the cost. Everything
away from damage stays here forever.

**Tier 2 — full nonlinear tet FEM, adaptive.** Around an impact, a fire, or a
growing crack, a region is *promoted*: the reduced model is replaced by genuine
3D tetrahedra at full resolution, coupled to the surrounding Tier-1 model through
the retained interface DOF. Inside that region the physics is uncompromised —
real stress tensors, real plasticity, real fracture. Budget 10⁵–10⁶ elements
across all active zones.

This is not a compromise on the answer where the answer matters. A collision
loads a 20 m stretch of side shell; the response of the bow 200 m away is linear
and the reduced model gets it right. The nonlinearity is local, so the expensive
model is local.

### Constitutive model

- **Elasticity**: co-rotational formulation to handle large rotations with small
  strains; St. Venant–Kirchhoff or Neo-Hookean where strains get large.
- **Plasticity**: J2 (von Mises) flow with combined isotropic and kinematic
  hardening; radial return mapping. Rate dependence via Johnson–Cook, because
  steel is markedly stronger under impact loading and ignoring that overestimates
  damage.
- **Damage and fracture**: ductile damage accumulation (Johnson–Cook damage or
  Gurson–Tvergaard–Needleman for void growth, which is more physical for the
  thick sections that matter). Failure surface reached → element splits along the
  maximum principal stress plane. Element deletion is the cheap fallback but
  loses mass and looks wrong; the plan is mesh splitting with remeshing of the
  affected neighbourhood.
- **Buckling**: falls out of the geometrically nonlinear formulation provided the
  mesh resolves the plate panels and the initial imperfections are seeded.
  Stiffened panel collapse is the failure mode that actually breaks ships' backs
  and it must not be scripted.
- **Welds and joints**: modelled as cohesive-zone interface elements with their
  own (lower) strength. Structures fail at connections far more often than in the
  middle of a plate.

### Plasticity and ductile failure — **implemented**

`engine/sim/plasticity.{hpp,cpp}`, checked by `tests/test_plasticity.cpp`, hooked
into the element by `solidshell::elementPlasticUpdate`. This is the half of the
list above that makes the solid-shell element do more than spring back.

**J2 flow, backward Euler, solved rather than stepped.** The return map is exact
for radial loading: the consistency condition `2μ(E − P) = √(2/3)·σ_y(√(2/3)·P)`
determines the plastic flow from the total deviatoric strain and from nothing
else, so **one step and ten thousand steps to the same final strain give the same
stress to 10⁻¹²** — asserted for both hardening curves, with a deliberately
non-proportional path as the negative control, because that property is what
separates a return map from a forward-Euler update and a test that cannot see the
difference is measuring nothing.

**Hardening is isotropic, with linear kinematic available and defaulted off.**
Forward, the two are indistinguishable — growing the yield surface and moving it
give the same monotonic curve, which is exactly why a suite that only pulls
cannot tell them apart. They part on a reversal, by closed form: the elastic span
is `2σ_y(ε_p)` for isotropic hardening and **exactly `2σ_y0` for kinematic**,
whatever the prestrain. Kinematic is off by default because nothing measures it —
`StructuralMaterial` carries E, ν, ρ and a yield strength and no hardening curve
at all — and because it costs six doubles per integration point, 380 MB at the
Tier-2 budget, to buy a Bauschinger effect a monotonic ram does not exercise.
**Rate dependence is deliberately absent**: it would make the map viscoplastic and
destroy the step-size independence above, which is the best instrument this file
has. Until it exists the flow stress is quasi-static, so the model under-predicts
the resisting force and over-predicts the penetration.

**The hardening curve is fitted, not tabulated.** Swift, `σ_y = K(ε₀+ε_p)ⁿ`, from
the three numbers a tensile test reports — yield strength, true stress at the
ultimate load, true uniform strain — with the fit exact at both points and its
Considère necking strain equal to the uniform strain that went in. The
construction that looks natural is wrong and is kept as a guard: taking `ε₀ =
σ_y/E`, the elastic strain at yield, puts the curve **21% high at 0.2 plastic
strain**, because a power law forced through the yield point of a steel that has a
plateau is far too steep.

#### Element size is the failure criterion, not a detail

An element does not tear at a material constant. Everything after Considère
happens inside a neck whose width is set by the plate thickness, so an element of
in-plane size *l* containing that neck reads an average — the uniform part, plus
the necking part diluted by *t/l*:

```
eps_f(l) = eps_uniform + (eps_fracture - eps_uniform) * min(t/l, 1)
```

Both constants are measurable and neither is fitted. `eps_uniform` is Considère's
necking strain, computed from the hardening curve — so the failure model and the
flow curve cannot drift apart, and a test asserts they have not. `eps_fracture` is
the local true strain at fracture, `ln(A₀/A_f)` from the reduction of area. The
clamp at `l = t` is not a fudge either: an element no larger than the plate is
thick resolves the neck itself and should see the local fracture strain, which is
what the expression returns there. **So it interpolates between two measured end
points rather than between two fitted ones.** What is mesh-invariant is the
necking elongation `(eps_f(l) − eps_g)·l = eps_e·t`, asserted as an identity; and
the formula is checked against its own derivation by averaging a rectangular neck
of width *t* over an element of length *l* and demanding the answer back.

Triaxiality is the other half of "depending on the strain state": Rice–Tracey void
growth gives `eps_f ∝ exp(−3η/2)`, referred to η = 1/3 where a tensile test
measures the constants, so in-plane biaxial tension (η = 2/3) fails at exp(−1/2) =
0.607 of uniaxial. Below a cutoff at η = −1/3 voids close and no damage
accumulates at all, which also bounds the multiplier by *e*. Damage is
**accumulated** — `D += Δε_p/ε_f(η)`, failing at D = 1 — not compared, so a path
that wanders through several stress states spends the right fraction of its life
in each.

For AH36 on 20 mm plating that gives:

| element | l/t | ε_f uniaxial | ε_f biaxial |
|---|---|---|---|
| 20 mm | 1.0 | 0.799 | 0.484 |
| 50 mm | 2.5 | 0.408 | 0.248 |
| 100 mm | 5.0 | 0.278 | 0.169 |
| 200 mm | 10 | 0.213 | 0.129 |
| 400 mm | 20 | 0.181 | 0.110 |

which brackets the 0.15–0.35 that ship-collision practice uses at the mesh
densities it uses. The element measures its own size — `elementSize()` returns
`sqrt(mid-surface area)` and `volume/area`, so a *sheared* element reports the
perpendicular distance between its faces rather than the length of its slanted
edge — and `initialisePlasticState` resolves the failure strain from it. **The
failure strain is a property of the element, not of the material.**

**The limit, measured rather than asserted away.** This regularises the strain one
element must reach to tear. It does *not* make the elongation of a uniformly
strained gauge length mesh-independent, and nothing without a softening mechanism
could: with hardening and no softening the plastic strain in a bar never
localises, so every element reaches its own failure strain at nearly the same
load. Measured on a 200 mm bar with a 5% thinned point, elongation before the
first tear: **43.3 mm as one element, 65.9 mm as four** — a factor of 1.52. What
the criterion does deliver, and what a structural tear needs, is that the tear
starts at the right place and at the right strain for the mesh: the thinned point
is the one that tears in both meshes, each at its own ε_f (0.213 and 0.408).

#### Plasticity inside an EAS element

The seven enhanced parameters stop having a closed form the moment the material
is nonlinear. `α = −K_aa⁻¹K_ua^T u` is a statement about a linear element; with
plasticity, α is whatever satisfies `r(α) = ∫Gᵀσ(Bu + Gα) dV = 0`, seven nonlinear
equations solved by Newton on the **algorithmic** tangent — which is why the
return map returns one, and why that tangent is checked against a central finite
difference of the return map itself rather than transcribed and hoped for.
Skipping the solve is not a small error: the enhanced thickness modes are the
entire reason σ_zz relaxes through a bent plate, and a plate that cannot thin does
not yield where a real one does.

Two things had to be got right that are easy to get wrong:

- **‖r‖ is not a scale-free convergence measure here.** The enhanced thickness
  modes carry `E_ζζ`, so their columns of G are scaled by the Voigt transform's
  1/t² and K_aa inherits its square. Measured on a bent element, a residual of
  10⁻³ corresponds to an error in α of 6 × 10⁻¹⁸ against an α of 2.7 × 10⁻⁶ —
  twelve significant digits. Chasing ‖r‖ below its own floor costs forty
  iterations and moves nothing; the forty-iteration answer and the four-iteration
  answer agree to every digit printed. Convergence is taken on `|δ·r|`, the work
  the correction would do, against the element's yield energy `σ_y·V`.
- **A dead integration point makes the enhanced problem ill-posed**, and the
  consequence is not a wobble. K_aa loses rank as the tangent at the dead points
  goes to zero; measured on a plate torn under an in-plane strain gradient, with
  four of eight points gone the survivors were driven to a triaxiality below the
  damage cutoff, their damage froze at 0.78 while their plastic strain ran on from
  0.49 to 0.89, and **the element never finished tearing**. So an element drops its
  enhanced modes the moment any point fails and finishes its life as the ANS hex,
  re-running the step in which the point died so nothing ill-posed is committed.

#### Cost, one core, measured

| | |
|---|---|
| Return map, elastic point | **35 ns** |
| Return map, plastic point (Swift, Newton on the consistency condition) | **180 ns**, 5.2× |
| The algorithmic tangent on top of that | +9% |
| Enhanced-strain Newton, warm start | 4 iterations (5 cold) |
| Solid-shell elastic internal force, the reference | 273 ns |
| Solid-shell elastoplastic update, nothing yielding | 5.4 µs, **19.6×** |
| Solid-shell elastoplastic update, marching and yielding | 7.3 µs, **26.9×** |
| Per-point state | 120 B; 10⁶ elements × 8 points = **916 MB** |

**35% of that 7.3 µs is re-forming B on the rest geometry, which never changes.**
It is loop-invariant and hoistable, at 12 kB per element — against the 4.6 kB per
element the elastic path already spends on a condensed stiffness. That is the
first optimisation to reach for and it is a memory decision, so it belongs with
the Tier-2 solver rather than here.

The consequence for the budget in `07-fem-spike-findings.md` §6 is real and should
not be glossed: a 200 m² collision zone at 50 mm elements costs ~7 400 core-seconds
per simulated second elastic, and **~2 × 10⁵ core-seconds elastoplastic** — five
minutes of wall time on the 24-thread CPU against about two hours. Plasticity is
what makes the crush zone expensive, not the element.

#### What it cannot do yet

1. **Small strain in a co-rotational frame.** The strain measure is additive, so
   at 0.2 it differs from the logarithmic one by about 10%, and the failure
   strains above are quoted as true strains. Large rotations are exact; large
   *strains* are not.
2. **Element deletion, not splitting.** A failed point carries no stress at all,
   discontinuously, which an explicit scheme feels as a small shock. The maximum
   principal direction at the moment of failure is returned — it is the plane a
   tear would open on — but nothing consumes it yet.
3. **No rate dependence, no temperature, no Gurson.** See above for why the first
   is deferred deliberately rather than forgotten.
4. **No `StructuralMesh` consumer.** `plasticity::shipSteel()` exists because the
   material database this section plans does not; its elastic constants are
   asserted against `ah36Steel()` so the two cannot drift, but the hardening curve
   and the fracture strain cannot yet come out of a ship file.

### Element technology — revised after measurement, then built

The original plan was uniform linear tetrahedra throughout Tier 2. **The spike in
`07-fem-spike-findings.md` ruled that out**, and the reason is worth stating
because it is not the obvious one.

Linear tets lock in bending: measured error against beam theory is 63% at two
elements through the thickness, 32% at four, 11% at eight. Ship plating is thin,
so avoiding that error needs many elements through 20 mm of steel — and the
explicit stability limit is set by the *smallest* element dimension, so those
same elements collapse the timestep. Cost scales as h⁻⁴. The two constraints
close on each other and leave no workable resolution.

So Tier 2 is mixed:

- **Solid-shell / assumed-strain (EAS, ANS) elements for plating.** One element
  through the thickness, no locking.
- **Tetrahedra where the geometry really is three-dimensional**: castings, engine
  seats, thick brackets, and the crush zone once plating has folded and shell
  kinematics no longer apply.
- **Promotion from shell to tet** as an element crumples past that point.

### The solid-shell element — **implemented**

`engine/sim/solid_shell.{hpp,cpp}`, validated in `tests/test_solid_shell.cpp`,
recorded in full in `07-fem-spike-findings.md` §6.

An eight-node hexahedron with three translations per node — no rotational degrees
of freedom, no director, so it stacks against tetrahedra at an interface for free,
which is what the promotion path above needs. Locking is cured by **assumed
natural strain** for the transverse shear (Dvorkin–Bathe) and the thickness strain
(Betsch–Stein), plus **seven enhanced assumed strain** parameters (Simo–Rifai)
condensed out at element level. Both are parameter-free. Reduced integration with
hourglass control was rejected: it is cheaper, but its hourglass stiffness is a
tuned coefficient that no measurement sets, and getting it wrong either invents
zero-energy modes or stiffens the bending it exists to preserve.

**Measured, on the same 8 × 2 × 1 mesh, as a ratio to the closed form:**

| L/t | solid-shell | ANS without EAS | plain 8-node hex | linear tets |
|---|---|---|---|---|
| 10 | 1.003 | 0.820 | 0.569 | 0.199 |
| 100 | 0.996 | 0.813 | 0.018 | 0.006 |
| 500 | 0.996 | 0.813 | 0.0007 | 0.0003 |

At the slenderness of real plating the plain hex is 1 400× too stiff and the
linear tet 3 800× too stiff, where the solid-shell is within 0.4%. The middle
column is the accounting: assumed strains alone fix shear locking and leave the
thickness-locking penalty of 1.225 — which is exactly the ratio of the plane-stress
modulus to the oedometer one, so the enhanced modes buy a closed form rather than
a tuning.

**Cost, one core:** 21.1 µs to form an element stiffness (once, at promotion),
267 ns for an internal force step against 129 ns for `fem.cpp`'s linear tet. Per
element that is 2×; per square metre of 20 mm plating per simulated second it is
**37 s against 4.1 × 10⁶ s**, a factor of 1.1 × 10⁵, because the tet mesh pays in
element count *and* in timestep and the solid-shell pays in neither.

**Two corrections this produced.** The timestep claim above was wrong: a
solid-shell keeps its through-thickness stretch mode — deliberately, a crush zone
needs the plate to thin — so the stable step is `t/c_p` **regardless of in-plane
size** (measured flat to 0.1% from 5t to 50t). The win over a bending-resolving tet
mesh is real and measured at 6.7×, but it comes from the thickness element being
the whole plate rather than an eighth of it, not from in-plane sizing. And the
element's binding limit is not accuracy but **geometry**: the assumed strains are
exact only for an element prismatic through its thickness, and non-parallel faces
cost ≈ 90 × (offset/t)² in spurious stiffness. Keep the thickness direction near
the surface normal and change plate thickness at a seam, not across an element.

### Solver

Explicit central-difference time integration inside Tier 2 (standard for
impact/fracture, no global stiffness matrix, trivially parallel). Tier 1 is
implicit and cheap. Lumped mass matrix.

GPU: the element loop is a good compute-shader workload — gather nodal state,
compute deformation gradient, stress, internal force, write per-element forces,
then gather into nodes through a CSR adjacency. **Measured on a GTX 1070 Ti:
450–670 M element-updates/s**, roughly 100× a single CPU core and 4× the whole
24-thread CPU. Nodal forces are gathered rather than scattered, which avoids
float atomics entirely and fixes the accumulation order.

### Slow damage

Separate from the fast solver, on a very long timescale:

- **Fatigue**: rainflow counting on the Tier-0/Tier-1 stress history, Miner's rule
  against S-N curves per detail class. A ship that has worked hard in heavy
  weather for years has cracks where the hot spots are.
- **Corrosion**: thickness diminution by zone, driven by coating condition,
  ballast/cargo/atmosphere exposure, and time. Directly reduces section modulus
  and therefore the loads at which everything above triggers.

### Material database

Per material: density, E, ν, yield and ultimate strength, hardening curve,
fracture strain vs triaxiality, Johnson–Cook rate and thermal coefficients,
thermal conductivity, specific heat, expansion, melting point, and — critically —
**temperature-dependent strength reduction curves** (Eurocode 3 for structural
steel, which loses roughly half its yield at 550 °C and nearly all of it at
800 °C). Coverage: mild steel, higher-tensile grades AH/DH/EH 32/36/40, stainless,
aluminium 5083/5383/6082 (which loses strength at *200* °C — the reason aluminium
superstructures are a fire problem), GRP and sandwich laminates, timber, ferro-
cement, and HY-80/100 for naval hulls.

---

## 4. Fire, heat and gas

Fire matters because of what it does to the structure and the atmosphere, not
because of the flames.

### Combustion

Two-tier, like the FEM:

- **Baseline: multi-zone model.** Each compartment carries an upper hot layer and
  a lower cool layer with an interface height — the classical CFAST formulation.
  Cheap, well validated, correct for smoke filling and layer descent, which is
  what kills people.
- **Promoted: LES combustion** in the compartment of interest. Eddy-dissipation or
  a flamelet model on a coarse (10–20 cm) grid. This is FDS-class physics at
  game resolution, used only where the player is.

Fuel is real: each compartment has an inventory of combustibles with heat of
combustion, ignition temperature, pyrolysis rate and oxygen demand. Fire goes out
when it runs out of either, and an under-ventilated fire produces carbon monoxide
instead of carbon dioxide — which is the actual hazard.

### Gas transport

Species-resolved (O₂, CO₂, CO, unburnt fuel, water vapour, soot) advected through
the same opening network the water uses. A door is a door. Buoyancy-driven flow
through vertical openings uses the standard two-way orifice formulation (hot gas
out the top, cool air in the bottom, with a neutral plane in between).

### Heat

- Conduction through bulkheads and decks by an implicit FEM thermal solve on the
  structural mesh — the same mesh, so the temperature field maps directly onto the
  strength reduction in §3.
- Radiation between hot surfaces and to flame volumes (view factors precomputed
  per compartment).
- Convection to gas layers.

The coupling that makes this worth the effort: **a fire in a machinery space heats
a bulkhead, the bulkhead loses strength, the bulkhead fails under hydrostatic
load from the flooded space next door, and the flooding spreads.** Every step of
that is modelled by a different subsystem and none of them know about the others.

### Suppression

Water spray and deluge (droplet evaporation cooling gas and wetting fuel), CO₂ and
inert gas total flooding (oxygen displacement, with the compartment sealing
requirement that makes it fail if a door is open), high-expansion foam, dry
powder, and boundary cooling. All of them add water or gas mass to compartments,
which is to say all of them affect stability. Firefighting has sunk ships.

---

## 5. Aerodynamics and wind

- **Wind field**: mean profile with a logarithmic boundary layer over the sea
  surface, plus a gust spectrum (NPD or Harris), plus local disturbance around
  the superstructure.
- **Windage**: offline panel-method or RANS sweep over heading angles produces
  force and moment coefficient tables per ship; runtime interpolates. Heeling
  moment from wind is a required input to the IMO weather criterion and to any
  honest capsize model.
- **Superstructure flow**: funnel exhaust dispersion, helicopter deck turbulence,
  and the recirculation zones that decide where smoke goes.
- **Sailing vessels**: for anything with sails, a lifting-line or vortex-lattice
  solver over the sail plan with real angles of attack, sail trim, twist,
  reefing, and the coupled heel/leeway equilibrium. A square-rigger and a modern
  sloop are different aerodynamic problems and both should work.
- **Ventilation systems**: fans, ducts and dampers as another set of edges in the
  gas network, which is how smoke gets somewhere it should not be.

---

## 6. Free-surface fluids

The lumped compartment model is right for slow flooding and wrong for violent
water. The plan is not to replace it but to **escalate**.

Each compartment carries a state: `Quiescent` (lumped, analytic free surface) or
`Dynamic` (particles). Promotion triggers on lateral acceleration, fill fraction
in the sloshing-sensitive band, a nearby impact, or player proximity. Demotion
happens when kinetic energy drops below a threshold for a sustained period, at
which point the particle mass and momentum are integrated back into a lumped
level. Mass is conserved exactly across both transitions by construction.

Solver choice for the dynamic state: **FLIP/APIC on a sparse grid** for interior
water, because it handles the free surface, the pressure projection and the
coupling to compartment boundaries well, and because sparse grids match the
geometry (water is in a few rooms, not everywhere). SPH is the alternative and is
better for spray and jets; the current expectation is FLIP for volumes and a
separate SPH/particle system for spray, jets from breaches, and rain.

Uses: sloshing in partly filled tanks (a genuine stability hazard, and the design
driver for LNG carriers), green water on deck, water moving between rooms as the
ship rolls, the jet from a hull breach, bilge water, and the visual of a
compartment filling around you.

---

## 7. Propulsion, manoeuvring and machinery — **partly implemented**

### What exists — `engine/sim/propulsion.hpp` / `.cpp`

Propeller, rudder and MMG hull, identified separately so damage to one degrades
one term. The module produces body-frame surge/sway/yaw contributions in the
same convention as `engine/sim/ship.hpp` (+x forward, +y to port, +z up; yaw
about +z, so positive yaw rate and positive rudder angle both mean "to port").
It also carries a self-contained 3-DOF horizontal-plane integrator
(`Manoeuvring`) so a turning circle can be measured without dragging in flooding
and hydrostatics; coupling into the 6-DOF rigid body is a separate change.

**Propeller.** Open-water K_T and K_Q from a single equivalent blade section at
0.7R, evaluated over all four quadrants of (advance speed, shaft speed) as one
expression rather than four cases. Wake fraction (Va = u(1−w)) and thrust
deduction ((1−t)T) are explicit. The four-quadrant behaviour is structural, not
tuned: astern rotation in ahead flow gives astern thrust because the section
angle of attack says so, so a crash stop cannot produce forward thrust.

The choice of a blade-element form over the Wageningen B-series regression is
deliberate. The regression is 39 terms for K_T and 47 for 10 K_Q, is
first-quadrant only, says nothing about astern rotation, and a single
mis-transcribed coefficient produces a curve that looks entirely plausible and
is wrong. The blade-element form instead makes the invariants algebraic:
`eta = J·K_T/(2π·K_Q) = tan(β)·(L cos β − D sin β)/(L sin β + D cos β)`, which is
identically 1 at zero section drag and strictly below 1 for any positive drag —
so the efficiency ceiling is a property of the algebra, not of the coefficients.

**Rudder.** Fujii's normal force `F_N = ½ρA_R f_α U_R² sin α_R` with
`f_α = 6.13Λ/(Λ+2.25)`, the inflow speed and angle taken inside the propeller
race, and the hull-interaction factors (1−t_R), (1+a_H) and the effective lever
x_R + a_H·x_H applied to surge, sway and yaw respectively. Past the stall angle
the normal force falls away from the linear extrapolation towards the
finite-aspect-ratio broadside value instead of growing without bound. The race
term is written as `u_P(1−κ) + κ√(u_P² + 8K_T n²D²/π)` rather than the usual
`u_P(1 + κ(√(1 + 8K_T/(πJ²)) − 1))`; the two are identical for J > 0 and only
the first is finite at bollard pull, which is the condition a tug lives in.

**Hull.** The MMG standard-method polynomial (surge with R_0, X_vv, X_vr, X_rr,
X_vvvv; sway and yaw with the full cubic set), written in the form where the
reference speed U divides out term by term — the raw non-dimensional form has
v′⁴ and r′³ terms that are singular at U = 0 and this one is bounded there. The
published derivatives are identified in a frame with y and r to starboard; they
carry into this frame's port-positive convention unchanged because every term is
odd under the joint flip (v, r, Y, N) → (−v, −r, −Y, −N), which is asserted in
`tests/test_propulsion.cpp` rather than assumed.

**Measured behaviour** (KVLCC2 defaults, 15.5 kn approach, 111 rpm):

| quantity | model | published / expected |
|---|---|---|
| bollard K_T, P/D = 1.0, A_E/A_0 = 0.70 | 0.347 | ≈ 0.35 (B4-70) |
| bollard 10 K_Q, same | 0.493 | ≈ 0.48 (B4-70) |
| zero-thrust advance ratio | 0.849 | ≈ 0.85 (B4-70) |
| peak open-water efficiency | 0.672 at J = 0.710 | ≈ 0.66 at J ≈ 0.68 |
| steady turning radius, 35° rudder | 1.13 L (360 m) | ≈ 1.1 L |
| drift angle at 35° | 19.4° | ≈ 17–20° |
| speed retained at 35° | 41 % | ≈ 40 % |
| steady turning radius, 20° / 10° | 1.76 L / 2.73 L | — |

### Which coefficients are real and which are placeholders

**Published, and used as published.** Fujii's `f_α = 6.13Λ/(Λ+2.25)`. The MMG
hull, propeller-hull and rudder-hull interaction structure and its KVLCC2
coefficient set (R_0′, the X/Y/N derivative polynomial, m_x′, m_y′, J_z′, w_P,
t_P, t_R, a_H, x_H′, l_R′, ε, κ) and the KVLCC2 principal particulars. **These
were transcribed from the MMG standard-method literature and have not been
checked against a primary source in this worktree** — the turning-circle and
drift-angle agreement above is evidence the set is at least self-consistent, but
treat any individual number as unverified until someone opens the paper.

**Calibration, not measurement.** Every free constant of the blade-section
model: `pitchEffectiveness` 0.853, `sectionLiftSlope` 4.0,
`sectionNormalForceMax` 1.4, `sectionDragCoeff` 0.018, `solidity` 0.194,
`solidityExponent` 0.2. They were fitted so the first-quadrant curve reproduces
the published B4-70 magnitudes in the table above; the fit is good but it is a
fit, and the propeller is not a specific propeller. The rudder stall constants
(`stallAngle` 30°, `postStallWidth` 20°, `postStallDrop` 0.25,
`broadsideCoeff` 1.20) are engineering judgement in the right range, not
identified data. `flowStraightening` is a single value (0.5) where the KVLCC2
identification lists 0.395 and 0.640 for the two signs of the rudder drift
angle; that asymmetry comes from single-screw propeller rotation, is not
otherwise modelled here, and using it would make port and starboard turns
differ for reasons the rest of the model cannot explain.

**Per-ship identification is the eventual answer.** These defaults exist so the
module has something to be tested against, not because a ferry is a VLCC.

### Validity limits

- **Deep water, no waves, no current, no bank.** Shallow-water, bank and
  ship-to-ship interaction are all unmodelled and all matter in the places
  ships actually manoeuvre.
- **Horizontal plane only.** No heel-induced hydrodynamic terms, so the heel a
  ship takes in a hard turn (and the sway/yaw coupling that follows) is absent.
- **Astern thrust is exactly the mirror of ahead thrust** at zero advance. Real
  propellers make roughly 70 % astern because the sections are not symmetric.
  The model's symmetry is a clean, testable idealisation and is optimistic.
- **Crash-stop distances are optimistic** for a second reason: with no engine
  model the shaft reverses instantly and delivers full torque astern. The
  measured head reach of 3.0 L is a fraction of a real VLCC's 8–15 L. This is
  the prime-mover item below, not a propeller problem.
- **No cavitation, no ventilation, no propeller emergence.** The propeller is
  assumed fully submerged in solid water at all times.
- **The blade-section model is a one-radius approximation.** It has no blade
  count dependence at all, and only a weak, calibrated dependence on blade area
  ratio; radial load distribution, skew and rake do not appear.
- **The MMG polynomial is a fit within its identification range.** Non-dimensional
  yaw rate r′ = rL/U is clamped at 1.5 (`kMaxYawRatio`) before the polynomial is
  evaluated. A hard-over turn reaches about 0.89, so the clamp is inactive in
  normal manoeuvring — it exists so a near-stationary spinning ship produces
  something bounded rather than a confident extrapolation. Low-speed and
  pure-rotation manoeuvring needs a different model.
- **The rudder race model assumes the rudder is downstream of the propeller.**
  Going astern the augmentation is switched off and the rudder sees the bare
  wake, which is the right sign but not a real astern-flow steering model.

### Still to build

- **Propellers**: cavitation inception and its thrust breakdown, ventilation
  when the propeller emerges in a seaway, per-ship open-water curves (either
  measured tables or a verified B-series regression) instead of one calibrated
  section model.
- **Manoeuvring**: shallow water effects, bank effect, ship-to-ship interaction,
  current, and the heel coupling. Coupling the module's forces into the 6-DOF
  rigid body in `engine/sim/ship.cpp`.
- **Steering**: hydraulic steering gear with real rates and failure modes,
  azimuth thrusters, Voith-Schneider, waterjets, bow and stern thrusters with
  their speed-dependent effectiveness loss.
- **Prime movers**: diesel engine model with turbocharger lag, fuel rack, governor,
  cooling and lubrication circuits, and the failure modes that follow when those
  circuits are damaged or flooded. Steam plant (boiler, turbine, condenser) for
  older and naval vessels. Gas turbines. Diesel-electric and hybrid drivetrains
  with a real electrical load-flow model.
- **Ship systems as networks**: electrical, hydraulic, fuel, fresh water, bilge,
  ballast, fire main, compressed air — each a graph of sources, sinks, valves and
  pipes. Damage severs edges. Loss of the fire main because a pump room flooded
  is a real and common cascade, and it should be modelled as a cascade rather
  than as a scripted consequence.

---

## 8. Cargo, ice and the rest

- **Cargo**: bulk cargo as a granular material that can shift and liquefy (the
  mechanism behind a long list of bulk carrier losses); containers with real
  lashing loads that part in heavy rolling; RoRo vehicles with tyre friction and
  lashings; liquid cargo with sloshing; heavy lift with the crane-induced
  stability problem.
- **Ice**: hull-ice interaction loads for icebreaking, ice-strengthened scantlings,
  and **topside icing** — spray freezing on the superstructure adds high weight
  and has capsized fishing vessels within hours.
- **Grounding**: seabed contact with soil mechanics for the reaction, hull raking
  damage as a moving FEM load, and the stability problem of being partly supported
  by the ground.
- **Collision**: two *deformable* ships, both FEM-active in the contact zone. The
  rigid-body half is done — §2 "Hull-to-hull contact" — so what is left here is
  the deformation: the contact patch and force history it reports are already the
  load case a FEM-active zone would be driven by.
- **Submarines**: pressure hull with depth-dependent loading and collapse depth,
  main ballast and trim tanks, the fact that submerged stability is a different
  problem (no waterplane, so BM = 0 and only BG matters).
- **Crew and passengers**: damage control parties that take real time to reach a
  door, firefighting teams with air supply limits, and passenger evacuation on a
  social-force model over a ship that is listing — walking speed as a function of
  deck angle is a well-measured quantity and it collapses beyond about 20°.

---

## Which of these are load-bearing

If time runs short, the ordering that preserves the most of the concept:

1. Flooding + air + stability **(done)**
2. Waves and seakeeping — without a seaway, nothing else has excitation
3. Structural FEM — the deformation and tearing the concept is named for
4. Fire and thermal — the second casualty type, and the best coupling story
5. Fluids — mostly a fidelity and presentation upgrade over the lumped model
6. Everything else — breadth, and breadth is what makes it a *ship* simulator
   rather than a *sinking* simulator
