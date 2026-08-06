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
modelled — a towed model — and does not reach the reported answer when propulsion
is attached. It is not entirely unused there, which an earlier draft of this
paragraph claimed: it seeds the encounter frequency that sizes the settle and
record windows, and those are counted in response periods. Since
`accelerateSeconds` defaults to **zero**, the default powered sweep accelerates
through its own recording window with that window sized from a speed the ship
never held. Set `accelerateSeconds` whenever the prototype has machinery.

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

**The 233 MJ is more than the two rigid ships carried into it**, and that is not
an error. Their rigid closing energy is 223 MJ; the extra comes from the water
each hull has to accelerate, which the contact removes momentum from as surely as
it does from the steel. The striker decelerates along its own bow and carries
surge added mass; the struck ship is driven along her beam and carries sway added
mass, which is nearly her displacement again, so the pair's true closing energy is
301 MJ. A bound written against the rigid masses passes this run by half a
percent and would read as a model creating energy the first time anything moved.



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
that same Steiner term. `constraint.hpp` is what carries it in the FEM, and
`stiffenedSection()` here is the independent route its section properties are
checked against — the 130× above reproduces exactly (130.25, at identical area).

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

**How it is checked** (`tests/test_scantlings.cpp`). Four instruments do the
work, and the first, second and fourth each found a defect the others would not
have:

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
- **A sweep along the length, with the two element populations counted apart.**
  The half-open interval above fixed the double count and introduced the opposite
  one, which stood for as long as it did because a single cut has nothing to be
  compared with. `sectionElements()` decides membership with a 10⁻⁹ m tolerance
  either side of a seam and then extracted the cut through each panel *exactly*;
  a plane admitted by the first test could therefore lie a hair outside the panel
  itself, and then no edge changed sign across it and the panel was dropped
  without a word. One unit in the last place was enough: the ferry is 120 m over
  50 bays, so station 33 is `-60 + 120·33/50 = 19.200000000000003`, and a cut
  asked for at the 19.2 a drawing carries fell 3.6 × 10⁻¹⁵ m aft of the bay that
  owns the seam. **Named that way, eleven of the ferry's 51 stations lost every
  plate panel** — 1.80 m² became 0.43 m², and which eleven was decided by nothing
  but the direction the division rounded — while every stiffener survived,
  because the member branch clamps its interpolation parameter into the member
  and was never asked the same question. What came back was 23.8% of the area, which is exactly the
  stiffeners' share, with a plausible neutral axis and a plausible taper either
  side of it. The plate cut is now clamped into the panel it belongs to, so the
  two branches snap to a seam the same way; and the section is asserted over a
  fine sweep of the whole length rather than at points, with plate and stiffener
  counts asserted separately, because a total cannot tell "the plating vanished"
  from "the taper moved".

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

A second pass, of **31 mutants aimed at `sectionElements()` alone**, killed 27 —
two of them only by a segfault, with no `FAIL` line, which a harness grepping for
one would have scored survivors. It found two more real holes: the guard that
rejects a degenerate cut could be raised from a picometre to a centimetre and eat
a 5 mm gunwale strake unnoticed, and `SectionElement::width` was written and never
read anywhere in the engine, so any value at all would do. Both are now asserted
against the box's closed-form girth. The four remaining survivors are equivalent:
`xHi - (xHi - xLo)` is bit-identical to `xLo` on every panel of the ferry;
`makeStructuralMesh` never puts a member forward of all plating, so the forward
end is set by the panels either way; a single crossing is already rejected by the
degenerate-cut guard before it is used; and dropping the member branch's clamp
moves a root by under a nanometre, which is the difference between that branch and
the plate branch failing catastrophically at the same seam.

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

**Tier 1 — reduced 3D. The reduction is implemented; the whole-ship model is
not.** A 3D shell model of a region, condensed by **Craig–Bampton component mode
synthesis**: retain the interface DOF plus the lowest few fixed-interface normal
modes, discard the rest — `engine/sim/reduction.{hpp,cpp}`, below. Two figures in
the plan this paragraph used to state are now measured and both were optimistic.
The cost saving on the same patch of the ferry's side is **2800×**, not 10⁻⁵, and
"a few percent for the frequency range that matters" needs the frequency range
named: the standard cutoff at twice the band of interest buys **0.6% inside the
10 Hz hull-girder band and 8% up to 20 Hz**, and five or six times the band is
what buys a part in ten thousand. The pass that assembles two substructures now
exists (`assemble`, §below), and there is now something worth joining — a zone and
the plating round it, §*Tier-1 to Tier-2 coupling*. It can also now see a
**stiffener**, which it could not until recently and which matters more than a
patch-scale figure suggests: a longitudinal has no nodes and no elements of its
own, so a substructure that only read the mesh reduced stiffened plating as bare
plating, and a hull girder is mostly longitudinals. What does not exist is a mesher
that can produce a whole-*ship* substructure, so "everything away from damage stays
here forever" is still the plan rather than the code.

**Tier 2 — full nonlinear tet FEM, adaptive.** Around an impact, a fire, or a
growing crack, a region is *promoted*: the reduced model is replaced by genuine
3D tetrahedra at full resolution, coupled to the surrounding Tier-1 model through
the retained interface DOF. Inside that region the physics is uncompromised —
real stress tensors, real plasticity, real fracture. Budget 10⁵–10⁶ elements
across all active zones. **The interface coupling is implemented**
(`coupling.{hpp,cpp}`) and is exact for a linear zone; what it cannot yet carry
back to Tier 1 is plastic softening short of a tear, because a reduction is linear
and the Tier-2 explicit solver forms no tangent operator.

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

### The Tier-2 zone — **implemented**

**How big a zone can usefully be, measured — and it is a cliff, not a slope.**
A zone stores a 24 × 24 stiffness matrix per element: **4.6 kB each**. This machine
has 30 MB of L3, so the working set crosses it at about 6500 elements, and the
parallel speedup does not degrade gracefully across that line — it falls off it:

| elements | working set | speedup on 27 lanes | effective bandwidth |
|---|---|---|---|
| 2 000 | 9 MB | 13.5× | 275 GB/s |
| 6 000 | 28 MB | 13.4× | 260 GB/s |
| 8 000 | 37 MB | 10.5× | 193 GB/s |
| **12 000** | **55 MB** | **2.7×** | **53 GB/s** |
| 20 000 | 92 MB | 1.6× | 29 GB/s |

Below the line the matrices stay resident and the solve scales. Above it every
element streams from DRAM and the bandwidth saturates at ~29 GB/s however many
threads are pointed at it. **So a big zone does not merely cost more, it scales
worse**, and the practical ceiling is 8000–10 000 elements rather than anything
the element count alone would suggest.

That also corrects the diagnosis this section was first written with. The zone's
own measurement — 3.4× on a 224-element patch — was attributed to the
synchronisation barrier, and "the barrier is what to attack, not the kernel".
Measured separately, `JobSystem::parallelFor` delivers **7.8× at 224 elements with
47 µs of overhead** for compute-bound work, so the barrier is not what is holding
it down. A kernel with the *same memory footprint* as a real element update gives
about 4× at that size, which is what is actually being seen.

Two caveats on that, both real. The 224-element point is small enough that a
single step is tens of microseconds and the measurement is noisy — repeated runs
gave 1.8× and 4.0×. And the probe is a synthetic matvec with the right footprint,
not the elastoplastic kernel, so it settles the *bandwidth* question and only
bounds the barrier one.

The per-element stiffness store is what puts the ceiling there, so the question is
what to do about it — and the obvious answer is wrong. **Recomputing instead of
storing costs 133× the bandwidth it saves**: forming an element stiffness is
21.1 µs (`07-fem-spike-findings.md`) while streaming its 4.6 kB at the saturated
29 GB/s is 0.16 µs. This paragraph first suggested exactly that trade without
costing it, which is the error this file keeps recording other people making.

What *is* available is that the matrix is **symmetric** — 300 independent entries
of 576 — so storing a triangle halves the working set and moves the knee from
~6500 elements to ~13 000 for no arithmetic at all. Beyond that, the honest answer
is to keep zones under the line rather than to make the line move.


`engine/sim/zone.{hpp,cpp}`, checked by `tests/test_zone.cpp`, run at ship scale by
`tools/zone_probe` in `verify.sh full`. The consumer the element and the return map
did not have: **`StructuralMesh` + a load → solid-shell elements → an explicit
solve → which panels tore**, as indices `breachesFromFailedPanels()` takes
unchanged.

```cpp
Patch  buildPatch(const StructuralMesh&, const Vec3& impact, const MeshParams& = {});
class  Solver { bool step(); const SolveResult& run(); ... };
double estimatedCost(const Patch&, bool plastic = true);   // core-seconds per simulated second
ZoneDamage indent(const StructuralMesh&, const Vec3&, const MeshParams&, const SolveParams&, ...);
```

#### Cost decides the design, so it is measured before the solve and after

The stable step is `t / c_p` — thickness governed, flat in the in-plane element
size — so for 12 mm plating it is 1.8 µs and 5.5 × 10⁵ steps buy one simulated
second. At 7.3 µs per elastoplastic element that is

```
core-seconds per simulated second = 4.0 x elementCount
```

and **the zone is bounded by element count, not by area**. Since the step does not
care about the in-plane size, area and resolution trade as `elements = area / h²`:
200 m² is 80 000 elements at 50 mm and 2 200 at 300 mm — two hours of wall time on
24 threads against four minutes. *Quoting an affordable area without its in-plane
size says nothing*, which is the correction this made to the figure in
`07-fem-spike-findings.md` §6.

The other half is that the event is short. A 6 m/s bow reaching 0.22 m into the
ferry's side is 0.037 s, not one second. Measured, on the reference ferry, a
three-metre zone at four elements across each 0.70 m bay:

| | |
|---|---|
| Zone | 14 panels → **224 elements**, 522 nodes, 24.0 m² of 12 mm plating |
| Promotion — meshing, and a power iteration per element for the step | **7 ms**, once |
| Predicted | 900 core-seconds per simulated second |
| Delivered, 23 workers | 21 290 steps, **4.5 s of wall time**, 0.94 µs/element/step |
| The same run, one worker | 15.4 s, 3.24 µs/element/step |

The element loop is dispatched through `core::JobSystem` and the nodal accumulation
is a CSR gather in a fixed order, so **the threaded answer is bit-identical to the
serial one at any worker count** — asserted, not hoped for, because that is the
property replays and multiplayer rest on.

> **Half of that was work that did not have to happen, and profiling found it while
> costing a GPU port.** `elementPlasticUpdate` began by calling `computeForms`,
> rebuilding the element's strain-displacement matrices, its enhanced-strain
> interpolation, its Gauss weights and its rest Jacobian — every step, for every
> element, from the **rest** configuration, which an explicit solve never moves.
> `solidshell::RestForms` forms them once at promotion instead, and
> `zone::SolveParams::cacheRestForms` is on by default:
>
> | | rebuilt each step | cached | |
> |---|---|---|---|
> | 192 elements, 1 worker | 5.48 s | **2.73 s** | 2.01× |
> | 192 elements, 23 workers | 1.48 s | **0.90 s** | 1.64× |
> | 17 800 elements, 23 workers | 21.12 s | **13.03 s** | 1.62× |
>
> **The two answers are bit-identical**, asserted on every reported quantity and
> every node position rather than compared to a tolerance, and driven past first
> tear so the element-deletion path is compared too. The per-element elastoplastic
> cost therefore falls from 7.3 µs to **3.0 µs**, and `estimatedCost` is
> `1.7 × elementCount` core-seconds per simulated second for 12 mm plating rather
> than 4.0 — see `07-fem-spike-findings.md` §8, which also records that the first
> two instruments used to size this both reported 97% where the A/B on the real run
> says 51%.
>
> `fem.cpp`'s tetrahedron has had this since the spike — it uploads `restInverse`
> and `restVolume` and the shader reads them from a buffer. The solid-shell never
> grew the equivalent, and the reason it went unnoticed is worth stating: the
> per-element cost **was** measured, and the measurement was right. What was never
> asked is which part of it depended on the state being advanced. *A cost model
> built from a correct total can still point at the wrong optimisation.*

**Threading saturates at about 3.4×, and the reason is the step and not the
element.** Measured on the same zone: 15.4 s at one worker, 6.9 at four, 5.2 at
eight, 4.9 at sixteen, 4.5 at twenty-three. Fitting `s + p/N` puts the serial and
barrier cost at **188 µs per step against 542 µs of element work** — twenty-one
thousand barriers in a four-second run, each dispatching twenty-eight chunks. A
larger zone does better (736 elements reaches 4.3×) but not much. What is left to
attack is the per-step barrier and the per-step energy accounting, not the element
kernel, and any GPU path will meet the same structure.

#### Meshing: the element's geometry limit is the mesher's problem

`07-fem-spike-findings.md` §6 limit 1 is binding: the ANS interpolation is exact
only for an element prismatic through its thickness, and non-parallel faces cost
≈ `90 × (offset/t)²`.

**One normal per mid-surface node**, area-weighted over the elements around it,
and the element extruded ± t/2 along it. On flat plating the offset is then
*identically zero* — not small, zero, and the tests assert the identity. On curved
plating the spread of nodal normals across an element **is** the offset ratio, and
it is reported rather than incurred silently. Extruding each element along its own
face normal instead is exactly prismatic everywhere and was rejected: the elements
then share no nodes on a curve and the patch falls into loose plates.

Measured on the ferry, and the result is a constraint on where a zone may go:

| zone | offset/t | excess bending stiffness |
|---|---|---|
| flat of side, z = 9 m, r = 2 m | **0.0000** | 0% |
| flat of side, z = 8 m, r = 4 m | 0.019 | 3% |
| flat of bottom | **0.0000** | 0% |
| across the shoulder at z ≈ 4.2 m | 0.188 | **319%** |

The shoulder is a 43° facet between two girth bands where this hull reaches full
breadth. **Refining `MeshParams::subdivision` does not help**, and that is the
useful part: a panel is a flat facet, so all the turning is at the seam whatever
the subdivision. The cure is a finer girth layout in `Scantlings` — halving the
band halves the offset and quarters the penalty, asserted against the closed form
`offset = facet angle / 4` on a cylinder.

**Thickness seams stop a zone** rather than being averaged across, because a node
between a 12 mm and a 16 mm strake has no single position and splitting the
difference puts a taper inside both elements. The truncation is reported; so is the
thickness used, if a caller overrides it.

The mesher itself is validated against a closed form the same way the element was:
a patch meshed from panels, loaded by `uniformPressureLoad` and solved by
`solveStatic`, reproduces **Timoshenko's clamped square plate** — 0.00126 q a⁴/D —
converging −3.7%, −0.6%, +0.3% at 4×4, 8×8 and 16×16 elements.

#### Boundaries, and the stiffeners that are not there

The patch is cut out of a ship, so its edge is a lie either way. **Clamped**, and
the price is paid in zone size rather than in a tuned spring stiffness: the
boundary error is a Saint-Venant effect that decays away from the edge, so growing
the radius makes it go away and the convergence is measurable.

**Stiffeners were not meshed**, and the reason was that the only element in the
inventory is the solid-shell hex and there is no way to attach a web to a plate
with it that is not wrong: a web sharing one node row along the seam is a hinge
with a zero-energy tripping mode; a web widened to the plate's element size has its
strong-axis stiffness wrong by `(h/t_web)³`, 3 000× here; and smearing is what
`§3` above rejects with a factor of 130. What was needed is a multi-point constraint
tying an eccentric beam to a shell — the same machinery Tier-1/Tier-2 interface
coupling needs. It now exists: `engine/sim/constraint.{hpp,cpp}`, and
`Stiffeners::Modelled` builds the member out of it. See *Eccentric stiffeners*
below.

**Leaving them out is not the neutral choice, and measuring it is what showed
that.** With no supports the plating spans from one clamped zone edge to the other,
so the span it uses is the *zone radius* — a meshing parameter. That is exactly the
defect `indentation.hpp` records in its own history, where the size of the hole came
out a property of the contact radius rather than of the collision. So the two
honest readings are offered as the two bounds:

| | ferry side, 2 m punch, at 0.078 m | against a membrane on the real 0.70 m span |
|---|---|---|
| `Stiffeners::Ignored` | 6.6 MN at 0.35 m (spanning the whole zone) | ~5× too soft |
| `Stiffeners::RigidSupport` (default) | **18.9 MN** at 0.078 m | 10.6 MN — the expected ×1.8 |

`RigidSupport` pins every plating node a stiffener line runs through, which is what
a member far stiffer than its plating does in the limit; `Ignored` is the lower
bound; the truth is between and the bracket is published rather than a point
estimate. The default needs `subdivision ≥ 3` — the stiffener lines *are* the panel
seams, so at 2 a bay is left with one free node — and `Patch::freeFraction` says so
when it is not met.

> Those membrane figures are **per bay times the number of bays the punch covers**,
> and re-deriving them turned up how easily that is dropped. `indentationForce` on
> a 0.70 m span of 12 mm plate under a 2 m contact is 3.71 MN at 0.078 m and
> 12.05 MN at 0.35 m; the 2 m punch spans `2.0 / 0.70 = 2.86` bays, and
> `3.71 × 2.86 = 10.6` and `12.05 × 2.86 = 34.4`. Both published numbers
> reproduce. What does **not** reproduce is the "seven times too soft" in
> `zone::Stiffeners::Ignored`'s own comment: `6.6 → 34` is 5.15, which is what this
> section and `zone.hpp` §3 both say. The comment was the outlier and is corrected.

#### Eccentric stiffeners — **implemented**

`engine/sim/constraint.{hpp,cpp}`, checked by `tests/test_constraint.cpp`.

**The formulation, and why it needs no rotational degree of freedom.** A
solid-shell keeps `kDof == 24` — three translations per node and no rotation
anywhere — and carries bending as *differential displacement between its two
faces*. The rotation of a plate's cross-section is therefore already in the model,
as the pair of nodes through the thickness, and a point rigidly attached at
through-thickness offset `e` is an exact linear function of that pair:

```
u(e) = u_bottom + ((e + t/2)/t) · (u_top − u_bottom)
```

At `e = 0.1 m` on 12 mm plating the weight is 8.83, so the tie is an
*extrapolation* far outside the element — which is what an eccentric member needs.
It is exact for a **finite** rigid rotation, because a rotation is linear and
commutes with the interpolation: measured at 3.6 × 10⁻¹⁵ m after 0.70 rad and 2.9 m
of travel, which is under one ulp of the amplified operand.

**The member is a set of axial fibres**, two-node bars along the stiffener line,
each at its own offset and each tied that way. "Plane sections remain plane" is not
assumed; it is *imposed by the tie*, which is the physical statement that the web
is welded to the plate. The stations are **two-point Gauss through each rectangle
of the profile**, which makes the area, the first moment and the second moment of
each rectangle exact — so the fibres carry the profile's own `I_own` as well as its
Steiner term. One station per rectangle would lose `I_own`, and that is 25% of a
200×10 flat bar's second moment about the weld line.

**Against `scantlings::stiffenedSection`, two independent routes to one number.**
Impose pure bending about an axis `z₀` and the stored energy is
`E κ² (I_NA + A (z₀ − z_NA)²) L / 2` — a parabola in `z₀` whose minimum is the
neutral axis, whose value there is the second moment, and whose curvature is the
area. Sweeping `z₀` therefore checks all three of `stiffenedSection`'s outputs at
once, and an offset that is wrong moves the vertex. For a 200×10 flat bar on 12 mm
plating at 700 mm spacing:

| | finite element | `stiffenedSection` |
|---|---|---|
| stored energy, five axes from −0.05 m to +0.15 m | — | ratio **1.000000000** at every one |
| neutral axis, from the parabola's vertex | −0.02038462 m | −0.02038462 m |
| area, from the parabola's curvature | 1.04000000 × 10⁻² m² | 1.04000000 × 10⁻² m² |

Both routes are exact, so this is an identity and is asserted at 10⁻⁸ having been
measured at 2 × 10⁻¹⁰. It runs through a *linear* assembly on purpose: the same
comparison through `zone::Solver` picks up the co-rotational frame's own O(θ²)
error, which at `κ = 10⁻³` over a two-metre patch is **1.6%** and does not converge
away with the mesh — that is measured too, in `test_zone.cpp`, at `κ = 10⁻⁵` where
it is 10⁻⁴ and the comparison is about the formulation instead.

The negative controls are what make that non-vacuous. Tying the same fibres to the
mid-surface instead — zero eccentricity, which is what a wrong weight produces —
leaves the panel **exactly** the bare plate again, Steiner term and `I_own`
together, which is `§3`'s smearing argument written as an energy. And one fibre per
rectangle instead of two comes out 22.8% soft, which is precisely `I_own`.

**Tripping is not a zero-energy mode, and the reason is structural.** The member is
*condensed*: it has no degrees of freedom of its own, so it cannot add a mechanism.
Measured on a free-free stiffened patch, the whole spectrum: the bare plating has
exactly six zero eigenvalues and so does the stiffened one, with the threshold
scaled off the first elastic eigenvalue (1.17 × 10⁵ N/m) rather than fixed —
this repo has twice been bitten by a constant rigid-body cutoff. Tipping the web by
10⁻³ rad costs 29.2 J, and the tie carries the eccentric kinematics exactly: every
fibre moves laterally by its own offset times the rotation, to **0.000 × 10⁰ m**.

**And the honest half.** A bar has stiffness along its own axis and none across it,
so the fibres contribute *exactly zero* to that tipping motion — asserted as an
identity, not as "small". What restrains tripping here is the plating alone, and
the number is a closed form: a strip of width `b` with the seam down the middle and
both far edges clamped resists a seam rotation at `16 D / b` per unit length, met to
0.04% at 4 × 32 elements. So this formulation **over**-restrains tripping where the
hinge leaves it free — the web is forced to follow the plate's cross-section
exactly. That is the opposite error, and lateral-torsional buckling of a stiffener
stays `buckling.hpp`'s question.

**The stiffener re-introduces an in-plane length scale into the stable step, and
that is a cost decision.** `§1`'s "thickness-governed, flat in the in-plane element
size" is a property of the *plating*. A fibre's is not: the tie amplifies its
stiffness by the square of the weight — 373× for the outer fibre of a 200 mm bar on
12 mm plating — while its mass arrives at the pair unamplified, and `EA/L` grows as
the elements shrink. Measured on a 2.0 × 0.7 m strip:

| subdivision | plating step | with fibres | ratio |
|---|---|---|---|
| 2 (0.25 × 0.175 m) | 1.816 µs | 1.816 µs | 1.000 |
| 4 | 1.812 µs | 1.622 µs | 1.118 |
| 8 (0.0625 × 0.044 m) | 1.797 µs | 0.738 µs | **2.436** |

At the resolution the reference ferry uses the stiffener is free; past it, it sets
the step and the zone costs proportionally more. The figure is **exact rather than
bounded**: a rank-one stiffness against a diagonal mass has the closed-form largest
eigenvalue `k Σ vᵢ²/mᵢ`, checked against a power iteration on the same block.

**The mass is lumped equally over the pair, and the consistent alternative is a
trap.** `TᵀMT` preserves the total for any weight — the row sums are `(1−w) m` and
`w m` — but for `w = 8.83` those are `−7.83 m` and `+8.83 m`: a **negative nodal
mass**, which an explicit scheme cannot integrate. There is no lumping of a
positive mass over two nodes six millimetres apart that reproduces the first moment
of a mass a hundred millimetres away; the first moment is what the extrapolation is
for. So the total is exact and the stiffener's rotary inertia about the seam is
given up, and a test asserts every nodal mass is positive so the code cannot
quietly go back.

**Where it lands.** `Stiffeners::Modelled` on the flat reference strip, 50 mm into
a 200 mm half-span:

| | force | work |
|---|---|---|
| `Ignored` | 4.90 MN | 159 kJ |
| `Modelled` | **5.54 MN** | **179 kJ** |
| `RigidSupport` | 7.72 MN | 241 kJ |

Inside the bracket on both, which is the thing only a solve can say — a tie with
the wrong sign would put the member on the other side of the plate and could come
out either way. `RigidSupport` stays the default: it is what every figure in this
section was taken with, and it is the setting that makes the comparison against
`indentation.hpp` an equal one.

**What the arithmetic costs is nothing; what the step costs is everything.**
`fiberForces` is **20 ns per fibre per step**, measured over 20 000 evaluations of
a 240-fibre set. At the 0.42 fibres per element the ferry's resolution delivers
that is **8.4 ns per element per step**, or 0.27% of the 3.1 µs elastoplastic
element — invisible in an end-to-end A/B, which came back at +1.0% and −0.1% on two
runs of the same case. So the only cost that matters is the stable step in the
table above, and it is a cost of *resolution* rather than of stiffeners.

**What it still cannot do**, stated rather than discovered: no tripping (above);
no weak-axis second moment, because every fibre sits on the stiffener line, so a
tee with a wide flange loses more than a flat bar does; and **the stiffener does
not tear** — a fibre yields on its own uniaxial flow curve and hardens, but it has
no damage variable, so a zone whose plating has torn away from under a longitudinal
still has the longitudinal. That last one is what `promotion::reduce` is waiting
on, and it is a failure criterion rather than more constraint machinery. The GPU
path does not carry fibres either, and `Solver::adopt` says so rather than
absorbing it.

#### What mutation testing found here

33 mutants, each a single plausible edit to `constraint.cpp`, to the new code in
`zone.cpp`, or to `solveStatic`'s extra-block path. The first pass killed 27. **All
six survivors but one were real holes**, and closing them added five checks:

- **A fibre's plastic strain was stored unsigned.** `state.plasticStrain +=
  increment` instead of `direction * increment` passes every tensile test, because
  the sign only matters to the *next* trial stress. A stiffener on the far side of
  a dent is in compression, and the error is invisible until the load turns round.
  The suite now yields a fibre at −1% strain and releases it: it has to come back
  in tension at `E |ε_p|`, capped by the flow stress it has hardened to.
- **A run of exactly two stations built nothing.** A member crossing a single
  element of the patch is a member; the off-by-one was invisible because no zone
  large enough to mesh is that small. Checked on `addStiffener` directly.
- **Fibres bridged a gap in the plating.** The check that a fibre may only span an
  element *edge* survived its own removal, because nothing in the fixtures had a
  member the patch carries in two pieces. There is now a U of panels whose two arms
  both carry the seam and whose spine does not: 1.4 m of fibres over a 2.1 m
  member, and a bridged gap would show as 2.1.
- **`zone::Solver` never lumped the fibres' mass.** `lumpFiberMass` was tested;
  the solver *calling* it was not, which is a different statement. A zone whose
  stiffeners weighed nothing would accelerate too easily and no force or energy
  would say so.
- **A `DofBlock` reaching outside the element bandwidth was dropped in silence** —
  `BandedSpd::add` discards an out-of-band entry without a word. It does not bite
  the stiffener, whose four nodes are already banded, but `DofBlock` is the general
  shape an interface coupling has and it would bite that.

  **And the first version of that test passed against the mutant**, which is worth
  recording on its own. A spring whose off-diagonal is dropped while its diagonal
  survives is two enormous springs to earth: both nodes go to zero, they agree
  perfectly, and "the block tied them" is satisfied by a coupling that is gone. The
  check that works is *what sets the tied deflection* — the plate's compliance
  (5.4 × 10⁻² m) or the spring's stiffness (`load/k` = 2 × 10⁻⁸ m). A factor of
  2.7 million between them.

The one survivor is argued **equivalent**: the return map's first estimate of the
plastic increment can drop the hardening modulus, `f/E` instead of `f/(E+H)`,
because it is only the Newton start and one step lands on the root exactly for a
linear curve. Kept as the control, because a mutation harness that reports
everything killed is reporting nothing.

One mutant is worth a note on *cost*: shrinking the fibre timestep by dropping its
square root made the suite take twenty minutes instead of one, because every zone
solve then ran 30 000 times as many steps. It is killed by the reference-resolution
assertion, but a mutation sweep over an explicit solver needs a per-mutant timeout
or one bad edit eats the budget.

#### Against the membrane model, and one thing it settles

`indentation.hpp` on the same bay is the load-bearing comparison, and it is set up
so that the membrane model's own idealisation holds: a long strip, a rigid line
punch across it, boundaries held. Two models built from different physics:

| at 0.24 m into a 0.8 m span of 20 mm plate | force | energy |
|---|---|---|
| rigid-plastic membrane | 11.7 MN | 1.51 MJ |
| solid-shell FEM | 24.8 MN | 3.00 MJ |
| ratio | 2.12 | 1.99 |

**And the factor of two is predicted, not accepted.** The membrane model carries
the tension at σ_y flat and in one direction. The FEM carries it at the hardening
flow stress averaged over the path — 509 MPa against 355 — and under the
plane-strain constraint the clamped side edges impose, which raises the yield stress
to 2/√3 of it. Correcting the membrane model by those two closed forms brings the
disagreement to **20%**, and what is left is bending and the strip's own ends. The
run is quasi-static (kinetic energy 0.3% of the work) so the comparison against a
rate-independent model is fair.

They part company on *where* a bay tears: the FEM tears at 0.055 m where the
membrane model says 0.092 m, earlier by the predictable direction, because the
membrane model spreads one strain over the whole leg where the FEM has the punch
edge and the clamped support concentrating it.

> **What this settles, and it was a correction to `indentation.hpp`.** `impactDamage()`
> took the membrane span as the **frame** spacing — 2.4 m on this ship — and the
> struck width as the longitudinal spacing. For a longitudinally framed side that is
> the wrong way round: the plating spans between *longitudinals*, 0.70 m. The FEM has
> no span in it at all, only plating and supports, so it is the instrument that can
> say. On the ferry's own side at 0.078 m of penetration it resists at **18.9 MN**,
> against 10.6 MN for a membrane on the 0.70 m span and **1.10 MN** for one on the
> 2.4 m span. The short span is right, and it is now what the code uses.
>
> **But "a factor of ten" was the wrong thing to conclude from it,** and this draft
> said so before the fix was measured. Force differs by ten; the *hole* does not,
> because the hole is set by energy per unit struck area and the span very nearly
> cancels out of it. The energy to tear a bay is `σ_y·t·A·ε_f`, reaching the span
> only through the failure-strain regularisation, which is flat here: **0.646 MJ/m²
> on the long span against 0.680 on the short, 5%.** Running the pre-fix energies
> through the corrected span moves the published holes from 27 / 107 / 234 / 415 m²
> to 25.6 / 102.0 / 222.4 / 392.8 — five per cent, not ten times. What *was* wrong
> by 3.4× is force and penetration, 0.686 m of denting against 0.205, in the
> direction of reporting the hull far softer than it is.
>
> Reasoning from a force ratio to an area ratio is the mistake, and it survived
> three documents because each copy quoted the last rather than the measurement.

#### The energy account, and where it stops closing

Work in = strain energy + plastic dissipation + kinetic + damping, checked every
run. Measured on a 32-element bay: **−0.15%** elastic, **+0.43%** plastic,
**+0.45%** plastic with damping (which is in the account, or every balance above
would be satisfied by a solver that quietly removed energy).

Two limits came out of it and neither was assumed:

- **The residual is not an integrator error.** Quartering the step leaves it where
  it was. It grows with the *rotation* the elements carry — +0.9% at 0.01 rad,
  −0.03% at 0.05, −3.7% at 0.15 — which is the co-rotational frame's additive small
  strain measure, `§3`'s "what it cannot do yet" item 1, measured as an energy error
  rather than as a strain one.
- **Element deletion breaks it outright.** Past the first tear the account runs
  −70%: the deleted elements' stored energy vanishes, and dropping the enhanced
  modes at the step a point dies commits a strain jump that reads as dissipation
  nobody paid for. The balance is a check on the *solver*, and it is only meaningful
  before anything tears.

A third tie is worth more than either: **the plastic path with a material that
cannot yield reproduces the elastic path node for node**, to 8 × 10⁻¹⁶ m of 4 mm
travelled. Different force routine, different strain-energy expression, different
state — and the same trajectory.

#### What mutation testing found

56 mutants, each a single plausible edit. The first pass killed 43, and **every one
of the thirteen real survivors was a hole in the suite** — plus one deliberate
no-op, kept as a control, because a mutation harness that reports everything killed
is reporting nothing.

- **Three of them were in the weld**, and they share a shape: every mesh under test
  had its duplicate points either bit-identical or centimetres apart, so nothing
  probed the band the tolerance actually decides. Forgetting to square the tolerance,
  halving the bucket, and dropping the neighbour probe all passed. The fix in the
  tests is a gap of *one and a half tolerances* — the band a hundred-tolerance gap
  cannot see — and a shared corner placed one ulp across a bucket boundary, which is
  the real case the mirrored starboard panels produce.
- **The fix in the code came from the same place.** The probe was `±1` cell with the
  cell at twice the tolerance, which is correct and *silently coupled*: shrink the
  cell and the weld stops being a distance without anything saying so. The reach is
  now derived from `tolerance / cell`, which makes the cell size a free choice again
  and the coupling impossible to reintroduce.
- **Which way is out was never asked.** Every synthetic mesh happened to be wound so
  that the outward test never fired, so both dropping it and ignoring the caller's
  own direction passed. Panels wound the other way now have to produce an identical
  patch, and the ferry's own axes are asserted on both sides and on her bottom.
- **Area weighting of the nodal normals was invisible** because every mesh was
  uniform, and **the trapezoid was missing**: taking a quad's area from one doubled
  triangle is exact on a parallelogram, which every synthetic panel was. The same
  trap `test_scantlings.cpp` records finding in `PlatePanel::area()`.
- **A second real defect**: the stiffener proximity filter measured from a member's
  *midpoint*, so a girder running the length of a ship — midpoint two hundred metres
  away, passing straight through the zone — would have been dropped. It measures to
  the segment now. The mutant that put the midpoint back needed a test with a real
  radius in it before it would die, which is its own lesson: a test that sets
  `radius = 1000` has switched off every filter downstream of it.
- The spurious-stiffness penalty could be made linear in the offset instead of
  quadratic, and `estimatedCost` could forget the step count entirely, because both
  were printed and never asserted on.

Adding those took the suite from 164 checks to 179 and killed 52 of 56. **Four are
argued equivalent**: relabelling a quad's corners to reverse its winding gives the
same element (the formulation is symmetric in ξ and η); skipping the update for an
element every point of which has failed is output-identical to running it, because
the return map returns no stress there; dropping the member proximity filter is a
cost difference; and halving the weld's cell size is now compensated by the derived
reach, which is exactly what that change was for.

#### What it does not do yet

1. ~~**No coupling to Tier 1.**~~ **Done, and linear on the way back** —
   `coupling.{hpp,cpp}`, §*Tier-1 to Tier-2 coupling*. The patch's perimeter is
   driven from a Craig–Bampton model of the plating round it instead of being
   clamped, and a torn zone goes back to that model as a mesh with the dead
   elements deleted, so the surroundings feel the damage. What still does not go
   back is **plastic softening short of a tear**: a reduction is linear by
   construction and this solver forms no tangent operator, so between first yield
   and first tear the surroundings are told the zone is stiffer than it is — in the
   unsafe direction. Coupling to Tier 0 is unchanged and remains
   `promotion.{hpp,cpp}`'s.
2. **The indenter is kinematic and rigid** — a prescribed rectangular punch, no
   contact search, no friction, no release, and the striking body does not crush. A
   prescribed motion cannot run away, which is what makes it testable; a delivered
   *energy* needs the striking body's mass and is what `collision.hpp` would supply.
3. **Element deletion, not splitting**, so the hole is the deleted area — reported
   as whole panels because that is what `breachesFromFailedPanels` consumes.
   `tearFraction` decides at what share of a panel it counts, and neither end of it
   is right: a first dead element over-states a slit by the subdivision squared and
   requiring all of them under-states a tear that has crossed the bay. The real fix
   is a breach interface that takes an area.
4. **A GPU path exists, is now faster than the CPU, and still cannot be used** —
   and there is no rate dependence in the material, so the resistance is
   under-predicted by the 10–30% steel gains at collision strain rates.

   `engine/gpu/zone_gpu.{hpp,cpp}` is the Vulkan back-end for this solver, built to
   `fem_gpu.cpp`'s pattern with the EAS variables and the plastic history resident
   on the device. One invocation per element ran at **0.23–0.68× the 24-thread CPU**
   and got worse past 3 000 elements, because that mapping needs ~500 floats of
   dynamically indexed private state and Pascal spills it — 1 936 bytes per thread,
   confirmed from the driver's own pipeline statistics rather than inferred.
   `solidshell_forces_wg.comp` re-maps it to **one workgroup per element**, which
   spills 96 bytes and runs at **1.26–2.43×**, rising monotonically with size.

   What stops it being used is precision, not speed. At 768 and 3 072 elements the
   float kernel tears **40 and 247 elements against the double reference's 32 and
   162**, and plastic dissipation runs 26–34% high — while the negative control, the
   same double solver on a mesh jittered by 2 × 10⁻⁷ m, tears exactly 32 and 162 and
   moves the dissipation by 0.06–0.8%. **An earlier claim here that it "tears 60
   elements where the CPU tears none" does not reproduce and has been withdrawn.**
   The enhanced modes have since been normalised in `solid_shell.cpp`, and measured
   by A/B that changes nothing on this path, because the shader was already
   equilibrating Kaa. `07-fem-spike-findings.md` §8 has the measurements; the next
   thing to try is keeping `alpha` in double.

### Adaptive zone promotion — **implemented**

**End to end on the ferry, which nothing had run.** The pieces were unit-tested
separately; this is the whole path — a real ship, a real seaway, Tier-0, the
criterion, and a zone:

| condition | worst yield | worst buckling | worst collapse | promoted |
|---|---|---|---|---|
| still water | 0.079 | 0.088 | 0.079 | none |
| crest, hogging | 0.237 | 0.264 | 0.236 | none |
| trough, sagging | 0.083 | **0.298** | 0.128 | none |

Nothing promotes, and that is the right answer rather than a broken one: her worst
utilisation anywhere is 0.298 against a buckling threshold of 0.80. An intact
ro-pax in a standard wave does not need a finite-element model of her side, and a
criterion that said otherwise would be promoting for the sake of it.

The contact trigger fires exactly where the plating says it should. Over a bay of
12 mm plating spanning 0.713 m, `platingCollapsePressure` gives **402 kPa**, and:

| contact force | pressure | promotes |
|---|---|---|
| 1.0 MN | 221 kPa | no |
| **2.0 MN** | **442 kPa** | **yes, 464 elements** |
| 20.0 MN | 4421 kPa | yes, 464 elements |

Cost separates as designed: Tier-0 is 180–530 ms (the Smith sweep dominating) and
the decision that reads it is **4–5 µs**, five orders apart, which is why the
criterion is cheap to run often and the Tier-0 answer it consumes is not. A
promoted zone at 464 elements sits comfortably inside the ~8000-element ceiling
the cache cliff imposes.

**Two ways to drive this API wrongly, both found by doing it.** A contact patch
whose centre is *inboard of the plating* finds no panels — the ferry's shell at
midship is at y ≈ 9.16 m, and aiming at 8.5 m with a 1.2 m radius reaches nothing.
And a single `review()` on a fresh `Promoter` never promotes anything, because the
dwell requires a candidate to qualify on consecutive reviews; promotion lands on
review 2. Both look like refusals and neither is.


`engine/sim/promotion.{hpp,cpp}`, checked by `tests/test_promotion.cpp`, run at
ship scale by `tools/zone_probe` in `verify.sh full`. The piece that makes the
tiers one system: **which patches deserve Tier 2, what a solved zone hands back to
Tier 0, and what the girder hands it on the way in.**

```cpp
TierZero        tierZero(const Ship&, const Sea&, const StructuralMesh&, const Scantlings&, ...);
class           Promoter { Review review(const StructuralMesh&, const TierZero&, contacts); ... };
PreloadCheck    preloadFor(const HullGirder&, const StructuralMesh&, const zone::Patch&, ...);
SectionReduction reactionOf(const StructuralMesh&, const zone::Patch&, const zone::Solver&);
StructuralMesh  reduce(const StructuralMesh&, const SectionReduction&);
```

A zone couples to Tier 0 through a *section* rather than through retained
interface DOF. That is cruder than the plan above and it is the honest thing
available for the **beam**; everything below is arranged so that inserting
Craig–Bampton replaces the coupling and not the criterion. That has now happened
alongside rather than instead: `coupling.{hpp,cpp}` drives a zone's edge from a
reduced model of the plating round it (§*Tier-1 to Tier-2 coupling*), and this
section is unchanged by it, exactly as intended. The two answer different
questions — Tier 0 wants a thinner ship for a section modulus, Tier 1 wants a
displacement field at an interface — so they coexist rather than one replacing the
other, and §*Tier-1 to Tier-2 coupling* records why the *thinner ship* shape does
not transfer.

#### The criterion: what Tier-2 adds, not what is worst

The cost is `4.0 × elementCount` core-seconds per simulated second and it is
**linear in the number of zones** — measured, not assumed: two zones take
2.00× one zone of the same size. So "promote where utilisation is high" is not
usable. A ship at 0.9 of her buckling capacity everywhere is a *well* designed
one — making utilisation uniform is what scantling design is for — and promoting
everywhere is a ship-shaped Tier-2 model, three orders past affordable.

What Tier-2 adds over Tier-0 is an answer where the response is **local and
nonlinear**. So a station qualifies only when it is both past an absolute
threshold for its trigger *and* standing at least `localExcess` above the
**median** over the ship's own stations for that same trigger:

| trigger | promote | hold | why there |
|---|---|---|---|
| yield, `σ/σ_y` | 0.90 | 0.80 | highest: a section modulus is not an approximation to something better |
| plate/column buckling | 0.80 | 0.70 | lowest: past critical the problem is nonlinear and Tier-0 only reports a ratio |
| collapse, `M/M_ult` | 0.80 | 0.70 | where load shedding along the length starts to matter |

The median is the flat-ship guard and it is deliberately robust rather than a
mean. Measured on synthetic profiles: flat at 0.90 gives **0 candidates with the
guard and 31 without**; a hull worked to 0.88 amidships with light ends has median
0.88 but mean 0.72, and only the median still says "flat".

A **contact patch** is the other trigger and it is not treated the same way: a
beam cannot represent a local load at all, so there is no background to stand
above and the test is absolute — mean contact pressure against the struck bay's
own plastic collapse pressure `p_L = 4 σ_y (t/span)²`, the three-hinge mechanism
of a clamped strip. `span` is the **shorter side of the struck panel itself**,
which is the one definition that cannot be got the wrong way round — and getting
it the wrong way round is precisely the defect `indentation.hpp` records. On the
ferry's 12 mm side plating over a 0.713 m bay that is **402 kPa**, and the
criterion is bracketed against it: 80% of the hinging force promotes nothing and
120% promotes one zone.

Candidates are ranked by `utilisation / its own promote threshold` — not by raw
utilisation, because a buckling utilisation and a yield utilisation are not the
same currency — then thinned so no two zones sit closer than `2 × radius`, then
accepted until the element budget is spent.

**What it misses**, and none of it is small:

1. **Everything Tier-0 cannot see**: no torsion, no shear lag, no racking, no
   local pressure. A load that does not change `M(x)` and is not handed in as a
   contact patch is invisible — slamming, sloshing, a dropped weight, a fire.
2. **It is one-dimensional.** A station says *where along*, never where around the
   girth. The site is the panel nearest the extreme fibre on the centreline, and a
   beam cannot say which **side** is in trouble; only a contact patch can.
3. **A uniform overload promotes nothing**, by construction. That is a decision —
   a ship uniformly past capacity has a girder answer and `collapse.hpp` gives it
   — and it is the one most worth revisiting.
4. **It is a snapshot at a cadence.** Two costs, three orders apart and both
   measured on the ferry: the **decision** is 7 µs over 8900 panels, which is
   tick-cheap; the **Tier-0 answer it reads** is 167 ms, of which 137 ms is the
   Smith's-method sweep. So promotion is emphatically *not* reviewed every tick,
   and whipping and springing live below the cadence.
5. **It has no memory** — corrosion, fatigue and previous damage enter only once
   something has changed the structural mesh, which `reduce()` below does.
6. **The element count is estimated, not meshed**: panels within the radius times
   `subdivision²`, because `buildPatch` is 7 ms and the criterion runs over every
   station. It over-states wherever a fold or a seam truncates the patch, which is
   the safe direction for a budget.
7. **Nothing is evicted.** A collision arriving on a full budget is refused and
   reported, not traded against a zone already running, because stopping a zone
   mid-solve throws away its plastic history.

#### Chatter, and the two mechanisms against it

Both, because they catch different things, and each tested against its own
negative control — the same signal with the mechanism off, which has to chatter or
the test proves nothing:

| signal | mechanism | with | without |
|---|---|---|---|
| 0.82 ± 0.03 across a 0.80 threshold | hysteresis (hold at 0.70) | **1** promotion, 0 demotions | 5 promotions |
| 0.95/0.55 alternating every review | dwell of 2 | **0** promotions | 20 promotions, 20 demotions |

Hysteresis kills any oscillation that fits inside the band whatever its frequency;
the dwell kills one *wider* than the band provided it is faster than the dwell.
Neither subsumes the other. Guards both ways: a load that stays up is promoted on
exactly the review the dwell is satisfied, and a zone whose load has gone survives
`hold − 1` quiet reviews and is dropped on the `hold`-th.

#### The pre-load, and a correction to the obvious argument

`preloadFor()` reads the moment and the section at the patch's station and returns
`σ_xx(z) = M (z − z_na) / I` as a `zone::Preload`. It is imposed as an initial
**strain** — the rest configuration is the meshed one with the exact elasticity
displacement field for that stress state taken back out — so the stress at step
zero comes out of the same validated constitutive path as everything else and no
new code goes inside the element. Measured: a patch asked for 84 MPa carries
84.05 MPa of `σ_xx` and under 0.03 MPa of everything else, stores `σ²V/2E` to
within 0.1%, and **does not move**: over 20 ms — many fundamental periods — the
worst node travels 4.8 nm and the patch holds 1.5 × 10⁻⁵ J of kinetic energy
against 439 J stored. The bending form is equilibrated too, which is a separate
statement: its displacement field carries a curvature term in `x²` and a field
missing it would spring the moment it was let go.

It is refused where it would not be traction-free: a panel whose normal leans `φ`
out of the athwartships plane is left with `σ sin²φ` unbalanced on its own face,
so a transverse bulkhead — which carries no hull girder stress, which is why
`hullGirderSection` leaves it out — gets no pre-load and is told why.

**How much it matters is not the 24% the ratio suggests, and measuring it is what
showed that.** Three findings, each from a different measurement:

- **The capacity it spends is exact.** With no punch at all, a patch yields when
  the pre-load reaches σ_y and not before — measured between 0.99 and 1.01 of
  355 MPa, both signs. So a zone handed the ferry's 84 MPa starts with **76.3%**
  of the uniaxial capacity an unloaded one claims. That part of the argument holds.
- **Under a punch it moves *yield onset*, and by a lot — but only in the regime
  the criterion actually promotes in.** A patch pre-loaded to 0.9 σ_y first yields
  at **0.059** of the penetration an unloaded one needs. At the ferry's 84 MPa the
  effect is a few per cent, and at 84 MPa nothing would have promoted a
  girder-triggered zone in the first place: her worst buckling utilisation on a
  3 m crest is 0.26.
- **It barely moves *tearing*, and that is the correction.** The pre-strain is
  elastic — 0.9 σ_y is 1.6 × 10⁻³ — against a regularised failure strain around
  0.15, two orders larger. Measured: the first element lets go at **0.9942** of
  the unloaded penetration. "A pre-loaded zone fails earlier" is true of yielding
  and false of tearing.

**And at ship scale it makes the zone *stronger*, not weaker.** `tools/zone_probe
--no-preload` exists so this is a measurement rather than an argument. Driving a
2 m punch into the ferry's own side at z = 8 m, on the 3 m crest that puts
13.1 MPa of hogging tension through that plating:

| at 0.078 m of penetration | resisting force |
|---|---|
| zone told it starts unstressed | **18.90 MN** — the figure this file published |
| zone handed the girder's 13.1 MPa | **20.25 MN**, +7.1% |

The reason is geometry. The girder's stress runs along **x**; the membrane stress
a punch raises in longitudinally framed side plating runs across the bay, which is
**vertical**, because the plating spans between longitudinals. The two are
perpendicular, and von Mises subtracts the product of perpendicular stresses
rather than adding them: `σ_vm² = σ_x² − σ_xσ_z + σ_z²`, so a *tensile* girder
stress raises the transverse stress at which the plating yields — from 355 to
389 MPa at 84 MPa of pre-load, closed form. **So on the ferry's side above the
neutral axis in hogging, ignoring the pre-load under-states her resistance.** It
over-states it below the neutral axis, and in sagging, where the same plating is
in compression; the sign is a property of where the patch is and which way she is
bending, and it is not available from the magnitude alone.

First yield under a punch is a *bending* event at the clamp, which presents both
signs of surface stress at once, so both signs of pre-load bring first yield
forward. A membrane argument about the sign of *that* is not available and quoting
one would be wrong.

One more thing the tool found rather than a test: a pre-load is not damage, and
neither is numerical dust. The reduction below reported twice as many damaged
panels under a pre-load as without one, every extra one intact to six figures,
until it was floored at a part in ten thousand of the plate — a micron on 12 mm,
below the tolerance the strake was rolled to.

#### The reaction back: damage is a thinner ship

`reactionOf()` measures, per panel the zone meshed, the fraction of its thickness
still working: a torn element counts zero; an intact one counts the thickness it
**actually has now**, taken as volume over mid-surface area from its deformed
geometry — measured, not modelled from a plastic strain and an assumed Poisson
ratio; and the part of a panel the zone did not mesh counts as intact, because
nothing looked at it.

`reduce()` folds that into the `StructuralMesh` as a thinner ship, and **nothing
in Tier-0 changes at all** — `hullGirderSection`, `girderStress`,
`girderBuckling`, `collapseElementsAt` and `longitudinalStrength` already read a
thickness. Measured, thinning forty side panels near the sheer at midship:

| effectiveness | area (m²) | I (m⁴) | Z_deck (m³) | M_ult hog | M_ult sag |
|---|---|---|---|---|---|
| 1.00 | 1.8013 | 46.205 | 5.573 | 1.987e9 | −1.288e9 |
| 0.50 | 1.7524 | 44.986 | 5.348 | 1.934e9 | −1.177e9 |
| 0.00 | 1.7034 | 43.714 | 5.118 | 1.875e9 | −1.169e9 |

Losing that plating costs **5.6%** of the hogging ultimate moment and **9.2%** of
the sagging one. **Two things are not monotone and both are measured rather than
assumed:**

- **The section modulus at the *undamaged* fibre can rise.** Damage just above the
  neutral axis pulls the axis down and the keel's lever shrinks faster than the
  second moment does: measured, the axis moves 6.713 → 6.681 m and `modulusKeel`
  goes 6.875 → 6.890 m³, *upwards*. So a far-fibre section modulus is not a
  conservative reading of damage; the ultimate moment is, and that is what the
  conservatism claim is made against.
- **The sagging ultimate moment does not fall at every step**, rising by up to
  0.35% over part of the range, because a thinner panel buckles *earlier* and
  therefore sheds earlier and Smith's method lets the neutral axis migrate in
  response. It is in the model rather than in the quadrature — eightfold the
  curvature steps leaves it where it is — which is why the claim is "never more
  than intact" rather than step-to-step monotone.

**What the reduction does not carry**, in the un-conservative direction: the zone
meshes plating only, so the stiffeners running through a torn panel are left at
full strength. A collision that opens fourteen bays has certainly destroyed the
longitudinals in them. It needs the multi-point constraint that would let the zone
mesh a web at all. The other is that a **dented** panel is far weaker in
compression than a merely thinner one; `dentedCompressiveCapacity()` computes that
by Perry–Robertson from the measured deviation (`η = 6 w₀/t`, exact at `η = 0`)
and it is deliberately *not* folded in, because it is a compression-only knockdown
and an effective thickness is not. Folding it in needs a `collapse.hpp` that takes
a per-element imperfection.

#### A Tier-0 defect this coupling found

`longitudinalStrength` sized its progressive-collapse sweep from
`firstYieldCurvature`, which is set by the **weakest** element in the section.
That is safe only while no element is anomalously weak — and it stops being safe
exactly when it matters, because a damaged bay's critical stress falls as `t²`.

Measured on the ferry with those forty side panels at an eighth of their
thickness: first yield falls **25×**, and a sweep to six times it reports an
ultimate moment of **1.26 × 10⁸ N m** against a true **1.89 × 10⁹**. The signature
is unmistakable — taking the same plating away *entirely* reported
**1.87 × 10⁹**, so the hull girder got fifteen times stronger when material was
removed. Worse, before an unrelated epsilon was fixed a fully torn panel came out
at an effectiveness of 1e-16 rather than zero, leaving a plate 10⁻¹⁸ m thick that
`collapseElementsAt` kept, and the ferry's ultimate moment was reported as
10⁻²¹ N m.

The fix is `collapseCurve()`: size from `firstYieldCurvature` exactly as before,
and extend **only if the peak lands on the last point of the sweep**, jumping to
`extremeFibreYieldCurvature` — the classical first-yield curvature, ignoring
buckling, which no single weak element can move. An intact section never enters
that branch and its answer is **bit-identical** to the old sizing, which is what
makes it a fix rather than a re-tuning; the ferry's worst strength margin is
4.2366 before and after.

#### What mutation testing found

61 mutants, each a single plausible edit to the criterion, the pre-load field, the
reaction or the corrected sweep. The first pass killed 37 of 59 and **every one of
the 21 real survivors was a hole in the suite**; after the fixes below, 58 of 61
die. The three that live are one deliberate no-op, kept as a control because a
harness that reports everything killed is reporting nothing, and two argued
equivalent.

- **A flat profile cannot tell a median from a mean from a minimum.** Three
  mutants of the background statistic all survived, because the only profile under
  test was flat and all three agree on a flat profile. The fix is a *skewed* one —
  a hull worked to 0.88 amidships with light ends, where the median says "flat"
  and the mean does not — and a ship already past capacity everywhere with one
  station 0.15 worse, which separates a difference from a ratio.
- **Which fibre a zone lands on was never asked.** Flipping the collapse fibre,
  forcing yield to the deck and ignoring `deckInCompression` all passed: nothing
  looked at where the candidate actually was. Hogging now has to put the zone at
  the keel and sagging at the deck, on a section whose two fibres are 15 m apart.
- **A criterion sampled once is a criterion untested.** The contact trigger was
  probed at 20 MN and 0.2 MN, two orders either side of a 2.8 MN threshold, so
  using the wrong power of the contact radius passed. It is bracketed at ±20% now.
- **Ranking and dedup were invisible** because no test had two candidates with
  different triggers. Normalising the score by each trigger's own threshold is the
  difference between promoting a station at 0.85 of a 0.80 buckling limit and one
  at 0.92 of a 0.90 yield limit, and nothing had asked.
- **Two real defects in the tests' reach rather than in the code**: the demotion
  side of the hysteresis was never exercised (the chatter test never let anything
  promote), and the element budget was only ever tested against an empty one, so
  forgetting what was already running survived.
- **The pre-load's `x²` bending term** could be dropped, because every pre-load
  test used a uniform stress. A patch standing on its edge under a gradient, run
  forward 20 ms, now has to stay still — and would not, since the field would no
  longer be compatible with the strain it claims.
- **A reduction measured against the nominal plate thickness** survived every flat
  fixture, because on flat plating an element's volume over its mid-surface area
  *is* the nominal thickness. It is not on a curved one, where the extrusion
  follows nodal normals that disagree, so a bilge zone would have reported section
  lost the moment it was promoted. A cylinder patch, unsolved, now has to report
  nothing at all.
- Two survivors are argued **equivalent given the mesher**: the "part nobody
  looked at" term and the zero snap in `reactionOf` are both dead code while
  `buildPatch` meshes whole panels or none. That invariant is now asserted, so if
  the mesher ever starts clipping the tests will say so rather than the arithmetic
  quietly reappearing.

### Tier 1: Craig–Bampton reduction — **implemented**

`engine/sim/reduction.{hpp,cpp}`, checked by `tests/test_reduction.cpp`. The
missing middle: Tier 0 is the whole ship as a beam and knows nothing about where
stress goes within a section; Tier 2 is solid-shell elements over a patch at
`4.0 × elementCount` core-seconds per simulated second. Nothing could give a
structural answer for a whole hold, a superstructure, or the region between two
bulkheads. This is what fills that gap.

```cpp
class  Substructure { ... };                              // K, M, the partition, K_ii factored
Reduction craigBampton(const Substructure&, const ReduceParams& = {});
std::vector<double> recover(const Substructure&, const Reduction&, const std::vector<double>&);
bool   staticSolve(const Reduction&, load, held, state, std::string* problem = nullptr);
Validity checkValidity(const Substructure&, const Reduction&, const std::vector<double>&);
double Substructure::memberStress(std::size_t, const std::vector<double>&) const;
constraint::AttachedForms constraint::attachedForms(...);   // blocks and stress forms together
Eigenpairs symmetricEigen(const std::vector<double>&, int n);     // no third-party dependency
Eigenpairs generalisedEigen(const std::vector<double>&, const std::vector<double>&, int n);
```

Partition the degrees of freedom into **boundary** (the interface the substructure
shares with the rest of the ship) and **interior**. Keep every boundary DOF
exactly; represent the interior by the static constraint modes
`Ψ = −K_ii⁻¹ K_ib` plus the lowest `m` fixed-interface normal modes `Φ`. With
`u = T x`, `x = [u_b ; q]`, `T = [[Ψ, Φ], [I, 0]]`, the reduced pair is

```
K_r = [[K_bb − K_bi K_ii⁻¹ K_ib,  0     ]      M_r = [[M_bb + ΨᵀM_iiΨ,  ΨᵀM_iiΦ]
       [0,                        Λ     ]]            [ΦᵀM_iiΨ,         ΦᵀM_iiΦ]]
```

The `(b,q)` stiffness block is *identically* zero, not small:
`K_biΦ + ΨᵀK_iiΦ = K_biΦ − K_biK_ii⁻¹K_iiΦ = 0`. It and the modal block `Λ` are
stored as the closed forms they are, because computing `ΦᵀK_iiΦ` would be a
noisier route to an eigenvalue already in hand. The modal *mass* block is the
identity by the mass normalisation and is nevertheless computed, because it is
cheap and because asserting on it is what says the eigenvectors really came back
mass-normalised. The tests check all of it by forming `TᵀKT` and `TᵀMT` the long
way.

#### The properties are identities, and they are asserted as such

- **Zero modes is Guyan static condensation, and static condensation is exact at
  the interface** — for any load, at any mode count. Against
  `solidshell::solveStatic` on the same problem the two agree to **2 × 10⁻¹⁰ m of
  a 0.31 m** deflection, and **8 × 10⁻¹⁰ m of 0.027 m** on a real patch of the
  ferry's curved side. What limits it is the conditioning of two independent
  solves, not the reduction.
- **This corrects the usual statement of that property.** "The static interface
  response improves with mode count" is false: it starts exact and stays there.
  Measured with the load moved *into* the interior, the interface error is
  2.5 × 10⁻¹¹ m at zero modes and the same 2.5 × 10⁻¹¹ m at thirty-two. What the
  modes buy is the **interior** recovery — 9.7 × 10⁻⁴ m of error at zero modes,
  1.2 × 10⁻⁸ m at thirty-two — and dynamics.
- **Rigid body modes survive exactly.** A rigid interface motion is the exact
  static solution of the interior with no interior load, so `Ψ u_b` *is* the rigid
  interior field — reproduced to 4 × 10⁻¹⁰ m of one metre. A free-free
  substructure keeps exactly six zero eigenvalues, at 2 × 10⁻¹⁰ of the first
  elastic one. Because the mass goes through the same `T`, a rigid translation of
  the reduced model weighs what the substructure weighs and a rigid rotation
  carries the full model's own rotary inertia.
- **Reduced frequencies come down from above, monotonically.** Rayleigh–Ritz on a
  subspace of the full space, and the subspaces nest, so a reduction can only
  stiffen. The tests assert the *direction* — a reduction whose error had the
  wrong sign would pass any "close enough" test and fail this one.
- **Symmetry and definiteness carry over — but not the same definiteness for both
  matrices.** `M_r` is positive *definite*, because `T` has full column rank and
  `M` is a positive diagonal. `K_r` is positive **semi**-definite and for a
  free-free substructure exactly six-fold singular: it being positive definite
  would mean the rigid body modes had been lost. Both are symmetric to the last
  bit because the analytically symmetric blocks are computed on one triangle and
  mirrored — which makes the symmetry *check* vacuous on its own, so the test that
  carries the content is the independent `TᵀKT`.

#### Whether to reduce the mass too: yes, and the alternative is measurable

Both matrices go through the same `T`, and the two things one might do instead are
wrong in ways that can be shown rather than argued. Keeping `M_bb` alone — the
"Guyan mass" a static condensation leaves — hands the model **12.5%** of the
substructure's inertia on the test plate, so every frequency comes out high by the
square root of the ratio. And projecting `K` through `T` while getting `M` from
somewhere else loses the from-above property above, which is a statement about one
subspace and only true when both matrices are projected onto it.

The mass itself is **row-sum lumped**, matching Tier 2: the reduced region and the
zone that replaces it then agree about inertia at the moment of promotion, a
diagonal `M_ii` makes the fixed-interface eigenproblem a standard symmetric one,
`M_bi` vanishes identically, and total mass and rigid translation are conserved
exactly. The price is a slightly low rotary inertia in absolute terms, which
touches nothing asserted here because every comparison is against the *same* full
finite element model.

There is a sharper statement available and it is what mutation testing forced:
row-sum lumping is `m_a = ρ ∫N_a dV`, and since `Σ N_a x_a` is the isoparametric
map, `Σ m_a x_a = ρ ∫x dV` **exactly** — the lumped masses sit at the meshed
volume's own centre of gravity. On the ferry's curved, tapered plating that holds
to **2 × 10⁻¹⁴ m**, where splitting each element's mass evenly is out by
2 × 10⁻⁵ m.

#### The eigensolvers, written rather than taken

No third-party dependencies, for the same reason this repo has a hand-written PNG
codec. Two, because the problems are different shapes:

- **Dense**: Householder tridiagonalisation with the transform accumulated, then
  implicit QL with Wilkinson shifts. Checked against the second-difference
  operator, whose spectrum `2 − 2cos(kπ/(n+1))` and sine eigenvectors are exact:
  eigenvalues to 10⁻¹³, eigenvectors to 10⁻¹², orthonormality to 10⁻¹³, and a
  triply repeated eigenvalue resolved with an orthonormal eigenspace.
- **Subspace iteration** for the fixed-interface modes, because only the lowest
  `m` of several thousand are wanted. It iterates `X ← K_ii⁻¹M_iiX` on a block of
  `q > m` vectors, using `solidshell::BandedSpd` for the factorisation.

**Subspace iteration converges to the lowest modes; it does not prove it found
them,** and a skipped mode is the silent failure here — the reduced model would
still be symmetric, positive definite, and bound the full model from above, and
simply be missing a mode. So the count is verified by a **Sturm sequence check**:
an LDLᵀ factorisation of `K_ii − σM_ii` has exactly as many negative pivots as
there are eigenvalues below `σ` (Sylvester's law of inertia), and `σ` between the
last mode kept and the first discarded must give exactly `m`. The same inertia
count makes a frequency cutoff an *exact* request answered before any eigenvector
is computed.

Two defects the tests found in the first version of that machinery, both of the
kind that produce plausible numbers rather than a crash:

- **The block goes numerically rank deficient.** Repeated multiplication by
  `K_ii⁻¹M_ii` drives every column towards the lowest modes; on a 210-DOF interior
  with a block of 200 the projected mass matrix stopped being positive definite,
  the projected eigenproblem was refused, and *the whole reduction quietly fell
  back to Guyan* with one line in `problems` to say so. The cure is to
  re-orthogonalise the block in the mass inner product each iteration — two passes
  of classical Gram–Schmidt — after which the projected mass matrix is the identity
  by construction and there is no Cholesky left to fail. It is not free: it is what
  takes the reduction of a 47-mode ferry patch from 0.33 s to 1.5 s.
- **A failed eigensolve returned no pairs and the reduction read them anyway**,
  which is a crash where the right answer is a Guyan reduction and a stated reason.

#### How the interface is chosen, and what it costs

The interface is every node the substructure shares with anything outside it, and
that is a property of the cut rather than of the mesh, so the caller supplies it.
`nodesNearPlanes` covers a region cut out between two bulkheads; `nodesPinned`
turns a `zone::buildPatch` patch into a substructure whose interface is the
clamped perimeter it already has. `HexMesh::fixed` is deliberately **ignored** —
a substructure is a free component and what holds it is a constraint on the
*reduced* model — and `problems()` says so, because a dropped constraint that is
not reported is indistinguishable from one that was honoured.

**The interface must restrain the interior, and the factorisation does not catch
it when it does not.** `K_ii` is singular if the interface leaves a rigid body
mode free, but a mechanism gives an exactly zero pivot in exact arithmetic and a
tiny *positive* one in floating point, so `BandedSpd::factor` succeeds and the
solve returns nonsense. It is a geometric precondition, so it is checked
geometrically — at least three non-collinear nodes — before anything is factored,
and the tolerance is scaled by the interface's own extent so that an interface
collinear to within a rounding is refused too.

**Cost is quadratic in the interface and linear in the interior**: `Ψ` is
`n_i × n_b` dense and `ΨᵀM_iiΨ` is `O(n_i n_b²)`. That is the argument for cutting
at a bulkhead — a transverse section is a thin ring of nodes — rather than around
a patch, whose perimeter is a large fraction of it.

#### Reverse Cuthill–McKee is not unconditionally better, measured

The interior is solved banded and banded storage is `n(band+1)`, so the numbering
the mesher happened to emit decides the memory and the factorisation time. RCM
removes that dependence — except that with a minimum-degree start it came out
**59 against 41** on the test plate and **83 against 53** on one ferry patch,
*worse* than the numbering it replaced, while beating it 89 to 173 and 137 to 341
on others. Two changes followed: the start is the George–Liu pseudo-peripheral
node, and the result is compared against leaving the numbering alone and the
narrower of the two kept. Comparing is free — a bandwidth is one pass over the
adjacency — and it makes the ordering incapable of being a regression:

| substructure | interior DOF | as meshed | RCM alone | chosen | banded store |
|---|---|---|---|---|---|
| test plate | 210 | 41 | 65 | **41** | 0.07 MB |
| ferry patch, r = 2.5 m | 1038 | 173 | 71 | **71** | 0.60 MB |
| ferry patch, r = 3 m | 1374 | 53 | 71 | **53** | 0.59 MB |
| ferry patch, r = 5 m | 4182 | 341 | 119 | **119** | 4.0 MB |
| ferry patch, r = 8 m | 8742 | 437 | 167 | **167** | 11.7 MB |

The **reversal** in the name is measured to buy exactly nothing here: plain
Cuthill–McKee gives the identical bandwidth on all five. What reversal reduces is
the *profile*, and constant-band storage does not exploit a profile. It is kept
because it costs one line and a skyline solver would want it, not because it was
observed to help.

#### How many modes buy what accuracy — measured, and the rule needs care

On a 272-element patch of the ferry's own side (1902 DOF, 528 of them interface)
with half the interface free, against a dense generalised eigensolve of the same
full model:

| modes | the fixed-interface cutoff it is | worst error, modes < 10 Hz | worst error, modes < 20 Hz | reduce |
|---|---|---|---|---|
| 0 (Guyan) | — | 3.9 × 10⁻² | 4.9 × 10⁻¹ | 0.11 s |
| 2 | 16 Hz | 6.0 × 10⁻³ | 8.3 × 10⁻² | 0.12 s |
| 5 | 24 Hz | 5.0 × 10⁻³ | 1.5 × 10⁻² | 0.19 s |
| 10 | 34 Hz | 3.8 × 10⁻³ | 5.3 × 10⁻³ | 0.21 s |
| 20 | 56 Hz | 5.7 × 10⁻⁴ | 1.9 × 10⁻³ | 0.50 s |
| 47 | 116 Hz | 8.8 × 10⁻⁵ | 2.0 × 10⁻⁴ | 1.50 s |
| 100 | 301 Hz | 3.7 × 10⁻⁶ | 8.3 × 10⁻⁶ | 5.20 s |

Every one of those is *above* the full model's, and monotonically decreasing.

**The trap is that the cutoff is a frequency of the fixed-interface problem and
the band of interest is a frequency of the assembled one, and they are not the
same spectrum.** With part of the interface free the assembled model is far softer
than the fixed-interface one — first assembled mode 2.97 Hz against a first
fixed-interface mode at 14.7 Hz — so the standard "cut off at twice the highest
frequency of interest" keeps **two** modes here, not the dozens the assembled
frequency count would suggest. Two modes is 0.6% inside the 10 Hz hull-girder band
and 8% up to 20 Hz. **Five or six times the band, not twice, is what buys a part
in ten thousand**, and the header says so beside the default rather than leaving
"a few percent" to be discovered.

The other end of the same measurement is what decides when to call this at all: a
**small stiff substructure needs no modes**, because its first fixed-interface
frequency is already above the band and Guyan condensation is exactly right at the
interface however soft the load. `firstFixedFrequency` is computed and reported
even when zero modes are kept, so "no modes were needed" is a number.

#### Cost against the tiers either side, on the same plating

| tier | what it covers | one-off | core-seconds per simulated second |
|---|---|---|---|
| **Tier 0** | the whole ship, as a beam | — | **0.10** (`hullGirder` 0.995 ms at 100 Hz) |
| **Tier 1** | this patch, linear | 0.11 s to reduce (1.5 s at 47 modes) + 0.01 s to factor | **0.35–0.41** (1 ms implicit steps, 528–575 DOF) |
| **Tier 2** | this patch, nonlinear | 8 ms to mesh | **1155** (`4.0 × elementCount`, 1.72 µs explicit steps) |

**Tier 1 is 2800× cheaper than Tier 2 on the same plating** and about four times
the cost of the beam that covers the entire ship. The explicit step is 1.72 µs and
thickness-governed, so Tier 2 takes 580 000 steps where a linear implicit model
takes a thousand — that ratio, not the matrix size, is where the three orders of
magnitude come from.

**The interface sets the running cost, not the mode count.** Going from 0 to 47
modes moves a step from 350 µs to 409 µs, 17%, while the 528 interface DOF are 92%
of the reduced model either way.

#### What a reduced region cannot do, and what the caller must do instead

**A reduced model is linear by construction and that cannot be worked around
inside it.** `K` is formed once from the undeformed geometry and `Ψ`, `Φ` are a
fixed subspace derived from it. So it **cannot yield** — no stress state to
return-map, no plastic history, no path dependence; it **cannot tear, buckle or
contact**; it **cannot rotate**, because the stiffness is small-strain in the
global frame, unlike the co-rotational Tier-2 element; and it is **stale the
moment anything changes**, because thinning or tearing changes `K` and the
reduction has to be rebuilt.

The tests assert that plainly rather than describing it: a load a thousand times
larger gives a stress a thousand times larger, straight through 355 MPa without a
word.

What the caller must do is **promote the region to Tier 2** — `promotion.hpp`
already owns that decision and `zone::buildPatch` already meshes the result — and
the reduced model's job at that point is to stop answering and hand over its
interface displacements as the zone's boundary condition. `checkValidity` is the
trigger: it recovers the interior displacement, evaluates stresses through the
same `solidshell::elementStress` Tier 2 uses **and the attached members' through
the same rank-one form their stiffness came from**, and reports the peak von
Mises utilisation over both halves.

`Validity` carries the two halves apart — `platingVonMises` and `memberVonMises`
— as well as the governing `utilisation` over them. One number is what a
promotion trigger has to compare, because promotion is a decision about the
region and a region is no more linear than its most utilised part; two is what
says *which* half, because a plating at 0.6 under a longitudinal at 1.1 is a
section shedding load off a failed stiffener and the reverse is a panel buckling
between intact frames. They are not two readings of one stress: on the stiffened
cantilever the plating's peak is transverse bending one element *off* the seam —
75 mm, at the loaded end — where σ_yy is 61.3 MPa against an axial 18.9, while
the member's is pure axial 65.2 MPa at the far fibre, on the seam. The
von-Mises-to-axial ratio at the plating's own governing point is **3.08**.

A fibre is an axial bar, so its stress tensor has one non-zero entry and its von
Mises is `|σ|` *identically* — the two halves are therefore the same equivalent
stress computed the same way, and comparing them against one yield is exact. What
is given up is the fibre *model*, not the comparison: a bar carries no transverse
stress at all, and that omission has **no fixed sign** — a transverse tension of
half the axial lowers von Mises to 0.866 |σ|, a transverse compression of equal
size raises it to √3 |σ|.

**And the warning is late, not early**, which is the direction that matters and is
measured rather than asserted: a truncated basis cannot represent a
concentration, so the recovered peak comes out **under**-predicted — 9.25 MPa
against the full model's 11.08 MPa with Guyan alone, converging to 11.01 MPa with
64 modes. `checkValidity` is a trigger for promotion, not a strength check.

#### Stiffness and mass the elements do not carry

A `HexMesh` is elements, and **not everything structural is an element**. A
stiffener is represented as axial fibres tied to the plating and condensed onto
the plating's own DOF (`constraint.hpp`), so it has nodes nowhere and elements
nowhere, and a substructure that only looked at the mesh reduced a stiffened patch
as **bare plating** — measured at 3e-16 relative to an element-only assembly,
while one 200 × 10 flat bar across a 0.6 m patch is worth 7.9% of its displacement
field. Longitudinals carry a large share of a hull girder, so any whole-ship model
built on that would have been wrong in the unsafe direction.

`reduction::Attachment` is a `std::vector<solidshell::DofBlock>` and a nodal mass
array — exactly what `constraint::stiffnessBlocks` and `constraint::lumpFiberMass`
already produce, and exactly what `solidshell::solveStatic` already takes, so the
static reference a test compares against is the *same* description of the same
physics rather than a second one that could disagree. The old constructor
delegates to the new one with an empty `Attachment`, so "no attachment is what
this always did" is one code path rather than a promise — asserted bit for bit
anyway, because it is the negative control for everything here.

**Three things had to be right, and each fails silently otherwise.**

- **The sparsity pattern.** The CSR is built from element node adjacency and a
  block ties DOF no single element need share. The lookup is a binary search that
  on a missing column returns the slot of the *next* one, so an unchecked scatter
  lands a stiffener's stiffness on a neighbouring DOF — a plausible field, quietly
  wrong, which is worse than dropping it. The blocks are folded into the adjacency
  *before* the pattern is sized, and every scatter checks the slot it found is the
  slot it asked for.
- **The bandwidth.** `BandedSpd::add` drops anything outside its band without a
  word. The band here is taken from the assembled pattern rather than estimated
  from the elements, so folding the blocks in fixes this too — and the interior
  renumbering sees them, which it would not if the blocks were scattered
  afterwards. Measured: fibres between adjacent mesh stations sit inside one
  element and change the band not at all (41 → 41), while a stiffener modelled
  every second station takes it 41 → 65. The banded fill counts what fell outside
  the band and refuses the substructure if anything did, so "the pattern covers it
  by construction" is a measurement.
- **Mass.** A block is stiffness alone. Without the steel the model is stiffer and
  no heavier and every frequency comes out high, which is most of what this tier is
  for. It is lumped onto the same diagonal as the elements' — which is also what
  `zone::Solver` does, so the two tiers agree about inertia at the moment of
  promotion. A substructure given stiffness and no mass says so in `problems()`.

**Validated against closed forms, not against itself.** Prescribe
`u_x = (ε + κz)x` and nothing else: every fibre lies along x and the tie puts the
point at offset `e` at `(ε + κe)x` exactly, so what the attachment adds to
`uᵀKu` is `E L (ε²A + 2εκS + κ²I)` — the profile's area, first moment and second
moment about the plate mid-surface. Those three come from
`scantlings::profileSection`, computed from rectangle dimensions and owing nothing
to fibres; the second is cross-checked against `stiffenedSection` by the parallel
axis theorem and the first against `smearedThickness`. **Measured: 2.2e-13
relative.** Both `ε` and `κ` are non-zero because with either alone the cross term
`S` drops out, and `S` is the only quantity that knows which side of the plate the
web is on. Statically, the reduced stiffened model reproduces
`solidshell::solveStatic` with the same blocks to **5e-10** of the peak on a
stiffened plate and to **4e-14** on the stiffened patch a coupling drives, in both
cases against a bare-plating answer 7.9%–99% away.

**What it gives up, measured.** The fibre mass is split equally over the
through-thickness pair rather than condensed, because `TᵀMT` is *negative* on one
node of any eccentric tie (`constraint.hpp`), so about the seam the model carries
6.8e-4 kg m² for a 200 × 10 bar where its own eccentricity carries 2.1e-1 — a
factor of **312**, all of it. Tier 1 is implicit and could integrate a negative
diagonal where Tier 2 could not, but taking the consistent mass here and the
lumped one there would make the two tiers disagree about the inertia of the same
steel.

`checkValidity` used to walk the elements alone, so an attached member's own
stress was not in the utilisation: a stiffened region was judged by its plating,
and the stiffener makes the plating *less* utilised, so the promotion trigger
moved later on both counts. Measured on the stiffened cantilever, 58.2 MPa
reported against 65.2 MPa in the member — utilisation **0.164 against a true
0.184**, 11% low and low in the *unsafe* direction. It is closed by
`Attachment::stress`, which `constraint::attachedForms` produces alongside the
blocks; a `DofBlock` on its own cannot supply it, because a rank-one `s v vᵀ`
pins down the strain energy `s (v·u)²` and nothing that separates the stress from
the volume it acts in. A substructure given stiffness with no stress form says so
in `problems()` rather than reporting the low number silently.

**One result worth having, because it is not the obvious one.** On the test plate
the 200 × 10 bar raises the first frequency by **307%** and its own mass moves
that same frequency by **0.015%** — because a member that stiff makes its own line
a *node* of the first mode, which becomes plate bending beside a seam that does not
move. The mass is not negligible; it is 13.7% on the worst of the first twelve
modes, and every one of the twelve falls when it is added, none by more than the
worst nodal mass ratio allows. On a 60 × 6 member — one the panel can actually bend
with — the first frequency lands within **4.5%** of `√(I_stiffened/I_plate ÷
m_ratio)` from `stiffenedSection`, where a smeared panel of identical area and mass
is out by **2.5×**. That is the factor-of-130 Steiner argument `scantlings.hpp` §1
makes, arriving as a frequency.

#### What it does not do yet

1. ~~**There is no whole-ship Tier-1 model.**~~ **Done** — `section::buildChain`,
   §*A ship: a chain of sections* below. `section.{hpp,cpp}` cuts a region between
   two transverse planes and hands it over with its longitudinals attached and its
   panel roles *tied*, and a chain of those, each reduced once and assembled, is a
   model of the whole length: measured against the same length in one piece at
   1e-10 in `EA`, `EI` and `GJ` on the box girder and 4e-10 on the ferry.

   ~~What is left is the mesher's reach.~~ **Also done, and the diagnosis was worth
   more than the fix.** Two-bay sections used to work over 62.4 m of a 120 m ship in
   five disconnected islands, the longest unbroken run being 26.4 m, and
   `buildSection` refused the rest "with an inverted element". **Nothing was
   inverted.** Every refusal on this ship was a
   *collapsed* hexahedron: `makeStructuralMesh` emits 166 degenerate `PlatePanel`s
   (quads with two coincident corners — a bulkhead running into the keel, a deck
   strake running out at the stem), each of which extrudes to a triangular prism
   whose Jacobian is exactly zero on the closed edge and sound at every point the
   element is integrated at. `solidshell::smallestJacobian` samples the *corners*,
   which is right for a general hex and wrong for this one. See `section.hpp` §7:
   the reach is now **120.0 m of 120.0 m**, 49 of 49 two-bay windows mesh, reduce
   and solve, and the whole hull meshes as one 8 900-element component.

   A cut plane is still an interface, so it unties the junctions on it — see the
   section below for what that costs, which is the one thing about a chain that is
   worse than the monolith rather than equal to it.
2. ~~**Nothing assembles two substructures.**~~ **Done** — `matchBoundaries()` and
   `assemble()`. Two components meeting at an interface have the same displacement
   there, so their shared boundary DOF *are* the same unknown and coupling is
   scatter-add of both reduced pairs into one, shared boundary landing on the same
   row and column, modal blocks stacked. Exact: nothing is approximated at the
   join, and whatever error the assembled model carries came from truncating each
   component's modes.

   Validated against the structure it claims to reproduce — a plate split down the
   middle by element, each half reduced independently, against the same plate
   meshed in one piece whose free-free spectrum is formed densely from the operator
   and owes nothing to the reduction. Worst of the first four elastic modes:
   **16 at zero modes a side, 4.6e-3 at four, 2.6e-5 at twelve**, always from above.
   The zero-mode figure is the useful one — Guyan is *exact* for statics and
   1600% wrong for dynamics, which is the whole reason the normal modes are there.

   Two controls, because a wrong assembly still produces a full set of plausible
   frequencies. Halves that were never coupled have **twelve** rigid body modes
   rather than six, each floating free of the other; without that the test passes
   on an assembly that joined nothing. And the split must be by *element*, asserted
   through total mass — splitting by **node** would count the interface mass twice
   and pull every frequency down while looking entirely reasonable.

   Matching requires the **axis** to agree and not merely the position: a
   coincident node would otherwise couple x to y, and the assembled model would be
   wrong in the way that is hardest to see. It is now also checked *after* the fact:
   a `Joint` is data a caller can build by hand, so `Assembly::axisDisagreements`
   counts assembled rows that merged two different axes and says so.

   ~~It joins **two** components and only two.~~ **Closed** — `assemble(parts,
   joints)` takes any number, §*A ship: a chain of sections* below. The prediction
   recorded here was that the fix was to carry the boundary DOF identity through an
   assembly so the same position match could be run against it. **That was half
   right, and the half that was wrong is the interesting half.** Carrying the
   identity is right and `Assembly::boundaryPoint` carries it — a caller needs it to
   find an assembled row *geometrically*, which is what prescribing a plane-sections
   field on a chain's end cuts is. But it is not what generalises the assembly: what
   stopped at two was the bookkeeping, and N components settle it by a **union-find
   over every pairwise map**, which needs no identity at all because every map is
   still expressed in the substructures' boundary DOF the caller already has. Both
   routes are built and tested against each other; they produce the same matrix to
   the last bit, and folding one component at a time costs O(N³n_b²) against
   O(N²n_b²) — measured at ×5.1 on a chain of eight.
3. ~~**A zone's edge is still clamped.**~~ **Done** — `coupling.{hpp,cpp}`,
   §*Tier-1 to Tier-2 coupling* below. What remains is that a reduced model cannot
   represent a zone that has *yielded* without tearing, which is item 6 there.
4. ~~**Stiffeners are invisible**, inherited from `zone.hpp` §3 — the substructure
   is as good as the mesh it is given.~~ **Done** — `reduction::Attachment`,
   §*Stiffness and mass the elements do not carry* above. A stiffener is
   *condensed*, not meshed, and the two consequences of that have had different
   fates: its own stress **is** now in `checkValidity`'s utilisation, through
   `Attachment::stress`, while its rotary inertia about the seam is still not in
   `M` and cannot be without a negative nodal mass. The substructure is still as
   good as what it is given — it just can now be given more than a mesh.
5. **No damping.** `M_r` and `K_r` are what comes out; a modal damping ratio or a
   Rayleigh pair is a caller's to add, and there is no measurement here that sets
   one.

#### What mutation testing found

80 mutants, each a single plausible edit. The first pass ran 68 of them and killed
52; of the sixteen that lived, one was the deliberate no-op kept as a control —
because a harness that reports everything killed is reporting nothing — one was a
mutant of mine that turned out to be a no-op by accident, and **every one of the
remaining fourteen was either a hole in the suite or a defect in the code**. Two
were defects: the rank-deficient subspace block and reverse Cuthill–McKee being
worse than the numbering it replaced, both described above. The suite went from
103 checks to 134, and the final pass kills 75 of 80 with four argued equivalent.

- **The peak of a bent plate is very nearly a uniaxial stress at a free surface**,
  so dropping the shear terms from von Mises changed nothing, and neither did
  looking at one Gauss point instead of eight. The fix is a plate turned 45° —
  the stress comes back in the *global* frame, so a stress along the plate is
  equal parts `σ_xx`, `σ_yy` and `σ_xy` and dropping the shear halves it — under
  an **eccentric** in-plane pull, which puts `4P/A` on the face it acts on and
  `−2P/A` on the other so the `ζ = −1` Gauss points carry less than half the peak.
  Both guards are asserted, not assumed.
- **Row-sum lumping and volume/8 are identical on a uniform brick mesh**, so
  nothing distinguished them. The centre-of-gravity identity above does, on the
  ferry's own distorted plating, by nine orders of magnitude.
- **A Cholesky that accepts a zero pivot looks like one that refuses it** unless
  the zero is on the *last* diagonal: anywhere else the row below divides by it,
  produces a NaN, and the very next pivot test rejects the matrix anyway.
- **A wrong Wilkinson shift converges to the same answer** and only costs sweeps —
  54 as written and 64 with the shift halved, where dropping it altogether does not
  converge at all. The sweep count is deterministic, so it is asserted where a wall
  clock could not be.
- **A node count is not a second test for a degenerate interface.** Two nodes are
  collinear; the count survives only in the message. And *nearly* collinear is the
  case that matters, because exactly collinear is what a synthetic mesh produces
  and a real one never does.
- **A refusal is not enough without its reason**: a test that only checked
  `problems()` was non-empty passed a Cholesky failure being reported as a QL
  failure — and then passed again on `find("collinear")` matching inside
  "non-**collinear**" in a different message.

Four are argued **equivalent**, and saying why is the useful part. Counting
non-positive rather than negative pivots in the inertia count cannot differ,
because an exactly zero pivot is already perturbed to a positive tiny before the
sign test reaches it. Dropping the reversal in RCM gives the identical bandwidth
on all five test meshes, for the reason above. The guard that refuses an interior
renumbering which did not cover every node is unreachable while the renumbering
is correct — but it is what converts "only the first component is visited" from a
silent half model into a reported one, so two mutants share one guard and the
other one dies. And the branch that replaces a collapsed column in the subspace
block is defence the suite never triggers, now that the orthogonalisation stops
columns collapsing in the first place.

A fifth is worth recording because it was killed by *deleting* code rather than by
adding a test: the interface check tested a node count and then collinearity, and
changing the count threshold from three to two altered nothing, because two nodes
are collinear and one is and none is. The count was dead. It survives only in the
message, where it does tell a caller something the word "collinear" does not.

**A later pass over the member half of `checkValidity`** ran 38 mutants and killed
25 on the first attempt; every one of the thirteen that lived was a hole, and
three of them are worth carrying:

- **A uniform prescribed field makes every segment of a seam identical**, so a
  recovery that read the *neighbouring* segment's degrees of freedom returned
  exactly the right number — and the test that checked all sixteen fibres against
  a closed form passed. The field is now `(ε + κz)·g(x)` with `g` quadratic, which
  gives each segment its own elongation and keeps the form closed because
  `(x_b² − x_a²)/(x_b − x_a)` is `x_a + x_b`.
- **The worst fibre under bending is never the first and never ties**, so the
  sweep's lower bound and its tie-break were both unreachable from the sixteen-
  fibre seam. They are reached by a substructure carrying *one* member, and by one
  carrying the same member twice — which is what two identical longitudinals
  either side of a symmetric section do for real, and which is what makes
  "the first of the tied members" a specification rather than an accident.
- **Every fibre on a flat plate lies in a plane**, so the last four terms of the
  twelve-term stress form are identically zero and dropping one is a no-op. A
  diagonal brace — one end tied to the −ζ face, the other to the +ζ face of the
  next column — reaches them, checked against the change in distance between its
  two tied points rather than against the condensed form the stress comes from.

Two of the thirty-eight were killed **only by SIGSEGV, with no `FAIL` line at
all**, so a harness that counted failures rather than exit status would have
scored them survivors. The harness counts signals, exit codes and a per-mutant
timeout, and it compiles from a copy outside the repository so a killed run
cannot leave a mutant applied in the tree. The final pass kills 38 of 38.

### Tier-1 section mesher — **implemented**

`engine/sim/section.{hpp,cpp}`, checked by `tests/test_section.cpp` — which is where
every box-girder figure below is an assertion rather than a reading — and measured at
ship scale by `tools/section_probe`, which is where every ferry figure below comes
from and how it is re-checked rather than re-quoted. **A region of a real ship between two
transverse planes, meshed as solid-shell elements with its longitudinals attached,
with the two cut sections as the interface a `reduction::Substructure` is built on.**
It is what item 1 above was waiting for, and it changes two things that were
believed about the problem.

#### What it found before it built anything

**`makeStructuralMesh` produces three topologically disjoint panel sets.** Of the
reference ferry's 9 390 distinct panel corners, the number shared between a `Shell`
panel and a `Deck` panel is **zero**. Shell-to-bulkhead is zero, deck-to-bulkhead is
zero, and two *different* bulkheads share none either. The three roles are laid out
on three independent grids — the shell on girth fractions of each station, a deck on
fixed |y| lines clipped to the hull, a bulkhead on its own — so a deck edge lands in
the *middle* of a shell panel, missing the nearest shell corner by up to **0.31 m**.

So a section of this ship arrives as several disconnected surfaces however it is
meshed, and that is a property of the input rather than of the meshing. It is why
`Section` counts its connected components, says which of them reach the interface,
and measures how much free edge is lying on another surface's plating without being
joined to it — 370 m of it on an eleven-bay hold, at a worst gap of 9 mm, which is
the deck edge sitting on the shell to within the hull tessellation.

**The second half of the junction answer is a formulation limit and would bite on a
conforming mesh too.** A solid-shell carries its thickness as *geometry*, so two
plates meeting at an angle have two thickness directions and a shared node pair
would have to point between them. Measured on a box girder whose four plates really
do share corners: welding the right-angled corners extrudes those nodes along the
45° mean normal, thins the plating towards the corner by `cos 45`, and costs **9.4%
of EA and 8.9% of EI**, both in the unsafe direction. The mesher therefore refuses to
weld across a fold and reports the junction instead.

#### The validation that proves nothing, and the two that do

**A section's EA, its neutral axis and its EI are insensitive to whether the
junctions are welded**, and that is a statement about the test rather than about the
ship. Prescribe both cut sections to a plane-sections field and every
longitudinally continuous strip carries `σ = E ε` whatever it is attached to,
because the ends alone already say what its strain is. On the box, the section whose
four corners are **cut** reproduces `EA` to a relative **6e-13** and `EI` to **5e-7**,
against −9.4% and −8.9% for the welded one. The unjoined mesh is the *more* accurate
of the two for this question, so a mesher that welded nothing would score best on a
section-properties test.

Two things do see the junctions:

- **Torsion.** A closed cell carries torque by Bredt's `4A²/∮(ds/t)`; an open one
  does not. On the box the welded section gives **0.966** of Bredt and the cut one
  **0.083** — a factor of 11.6 at `L = 8 m`. The factor depends on length and that
  matters: **1.7 at 2 m, 11.6 at 8 m, 166 at 32 m**, because a short open section is
  held by the warping restraint of its own end planes rather than by torsion. A
  one-bay section would have shown almost nothing.
- **The lowest fixed-interface frequency**, which is the sharper instrument on a real
  ship. On the ferry's hold between x = −7.2 m and 19.2 m the whole *unjoined*
  section's first fixed-interface mode is **0.7785 Hz**; the decks *on their own* are
  0.7785 Hz to four figures and the shell on its own is **1.6999 Hz**. The softest
  thing in the section is a 26 m deck held on two edges instead of four, and adding
  the shell it ought to be welded to changes it by nothing at all. Tied — see the
  junction tie below — the same section is **2.3026 Hz**, above both pieces. (Every
  figure is bracketed by `Substructure::eigenvaluesBelow`, an exact inertia count,
  because the subspace iteration reports that it did not converge and a frequency it
  produced is not evidence on its own.)

  > **Correction.** This paragraph and `reduction.hpp` both read **3.4600 Hz** for
  > the shell alone until `tools/section_probe` — the program this section says
  > produced the figure — was re-run. It gives 1.6999 Hz, and no combination of
  > subdivision or member setting reproduces 3.46. Nothing tests a comment; this is
  > the fourth time that has cost this repository a published number.

#### Against Tier 0, and what closes the gap

`hullGirderSection` reaches `A`, `z_na` and `I` by summing `A`, `Az` and `Az²` over a
transverse cut, sharing no line of code with the mesher, the element or the solver.
On the ferry's hold at `subdivision = 1`:

| | A (m²) | z_na (m) | I (m⁴) |
|---|---|---|---|
| Tier 0, `hullGirderSection` | 1.80133 | 6.71317 | 46.2047 |
| Tier 0 less the members the mesher cannot attach | 1.72557 | 6.8534 | 43.746 |
| the section | **1.73266** | **6.86502** | **43.868** |
| the section, bare plating | 1.37058 | 6.99161 | 35.650 |

> **Correction.** The third row read 1.73122 / 6.85797 / 43.875 and the tool it names
> has never printed that: on master it gives 1.73394 / 6.86054 / 43.91346, and with
> the halo 1.73266 / 6.86502 / 43.86796. Only the bare-plating row reproduced. See
> the resolution table's correction above; the same drift, the same cause.

The mesher's shortfall against Tier 0 is an **accounting**, not a discrepancy.
`Section::attachedMemberArea + missedMemberArea` equals `sectionElements`' stiffener
area to 1e-9. What is missed is the three girders, which sit off the longitudinal
spacing and so pass through no node of a mesh whose nodes are panel seams, plus —
since the halo — four members lying wholly in the bay against the strake seam at
x = 19.2, whose runs are now two stations of different thickness and therefore no run
at all: 0.07576 m² against 0.07464 before. Correcting for them leaves **+0.41% on the
area** and **+0.28% on the second moment** — the middle row subtracts the girders'
second moment about the
section's own neutral axis; carrying the shift of that axis through as well gives
43.70 and +0.40%, which is the size of the residual either way. And even that is
explained: the frames and deck beams restrain the section's Poisson contraction,
which a beam model has no way to represent and which is worth +0.52% of the area,
measured by omitting them.

The **negative control** lands where the section properties say it must: bare
plating is 23.9% short on `EA` where Tier 0 says the stiffeners are 23.8% of the
section, and 22.8% short on `EI`. That is `reduction.hpp` §8's warning — one flat bar
was 7.9% on one patch — arriving at ship scale.

#### Resolution: the Tier-2 ceiling was never the binding constraint

`applyBeamLoad` prescribes the interface and relaxes the interior, so what it
reports **is** `½ uᵦᵀ K_r uᵦ` for the Guyan reduction of the section — the static
condensation identity, asserted against a real `craigBampton` rather than taken on
trust. Refining it is therefore refining the reduced answer:

| subdivision | elements | A_eff (m²) | z_na (m) | I_eff (m⁴) | GJ (N m²) | solve |
|---|---|---|---|---|---|---|
| 1 | 2 068 | 1.72945 | 6.86238 | 43.8169 | 3.6164e12 | 0.45 s |
| 2 | 8 272 | 1.73081 | 6.85933 | 43.8598 | 3.5998e12 | 3.22 s |
| 3 | 18 612 | 1.73126 | 6.85828 | 43.8728 | 3.5947e12 | 13.23 s |
| 4 | 33 088 | 1.73148 | 6.85777 | 43.8787 | 3.5923e12 | 37.56 s |

> **Correction.** The four rows above are what `section_probe --sweep=4 --no-reduce`
> prints. The ones they replace were not: before the halo the same program gave
> 1.73075 / 6.85784 / 43.86283 at subdivision 1 against a published 1.73122 /
> 6.85797 / 43.8749, a table that had drifted from the program it names and that no
> gate re-runs — `scripts/check-figures.sh` covers `ram_view` and not this. The halo
> then moved it again, by less than the drift: the hold's forward plane at x = 19.2
> is a strake seam, and a section that can now see it stops four member runs there.

**One element per panel is converged.** From subdivision 1 to 4 the area moves
0.12%, the neutral axis 4.6 mm and the second moment 0.14%, while the element
count moves sixteenfold. The one quantity with plate bending in it — the torsional
stiffness — moves 0.67%, five times as much, which is the same membrane/bending
split the junction measurements make. All three of the first move about half again
as much as they did before the halo, for a reason worth stating: a thickness seam
costs a **fixed** two dropped fibre segments per member per seam, so refining the
mesh dilutes it and a coarse mesh is where it is dearest.

So a hold-sized Tier-1 section is **2 068 elements**, a quarter of the ~8 000 where
Tier 2's per-element stiffness store leaves L3. That ceiling is a property of an
explicit solve touching every element every step and it does not apply here. What
does bind is the **interface**: a transverse cut of this ferry is 195 panel corners
and therefore 2 340 boundary DOF at subdivision 1 — the whole of the reduced model
at zero modes and 93% of it with the 178 the cutoff asks for — and
`ΨᵀMΨ` is `O(n_i n_b²)` — so refining costs sixteen times per doubling and buys
0.05%.

**`ReduceParams`' default is wrong at this size and says nothing about it.** With the
hold's lowest fixed-interface frequency at 0.78 Hz, the 20 Hz cutoff asks for **178
modes**; that takes 275 s and the subspace iteration does not converge in its 60
iterations. Guyan alone is 6.2 s and is exactly right at the interface for any static
load. A substructure softer than the ship it belongs to is a sign that the section is
not yet a component, not a reason to keep more modes.

#### Thickness seams: taper rather than stop

`zone.hpp` §2 stops a patch at a plate thickness seam. A section cannot: the ferry's
shell crosses four seams per side between keel and sheer, and stopping at each would
deliver eight loose girth bands. So each mid-surface node carries the area-weighted
mean thickness of the sub-quads round it and the taper lands inside one element
either side of the seam — `dt/t = 0.131` at worst on this ship, which
`07-fem-spike-findings.md` §6 limit 1 prices at 154% excess *plate bending* stiffness.

That sounds fatal and is not, for the same reason the junctions do not matter to
`EA`: the excess is on the plating's own bending about its own mid-surface, which is
0.0175 m⁴ of the ferry's 46.2 m⁴ — **0.038%**. Measured directly on a box whose flanges
step from 10 mm to 16 mm, `ThicknessSeam::Split` (every element exactly prismatic,
section disconnected at the strake) and `ThicknessSeam::Taper` differ in `EA` by
**4.5e-7**, at a taper the spike's rule prices at 610%.

**What the taper does cost is stiffener steel, and mutation testing is what found
it.** `constraint::addStiffener` turns one `plateThickness` into one tie weight per
fibre, `(e + t/2)/t`, and applies it to a node pair whose real separation is the
*local* nodal thickness — so a member run crossing a seam, where that thickness is a
mean of two strakes, puts its fibres at `e·t_local/t_run`. Measured on the ferry:
**47 mm** out on a 700 mm frame, a quarter of its Steiner term in the wrong place,
and invisible in `EA`, in `EI` and in the mass. A run therefore stops at a thickness
change, and the seam station — whose neighbours both differ from it — is left in a run
of one and dropped. That costs **5.2% of the section's member steel** (146.8 t against
139.1 t) and **none** of its longitudinally effective member area, because a
longitudinal runs along the ship at constant thickness and it is the athwartships
members that cross strakes. Missing steel is visible in a mass; a misplaced
eccentricity is not, which is the direction to fail in. The fix that costs neither is
a per-station thickness in `constraint::SeamRun`, and it is that file's change rather
than this one's.

#### Three smaller things worth carrying

- **Cut on a frame station, not at a bulkhead.** The obvious place to cut a ship is
  at a watertight bulkhead, and on this ship it is the wrong place: the bulkheads are
  at x = −44, −38, −8, 20 and 44 and the frames are multiples of 2.4 m, so a bulkhead
  plane passes through **188 panels**. They are dropped and counted, which is why the
  28 m bulkhead cut ends up with *less* plating per metre than the 26.4 m frame cut.
- **The interface has to be chosen on the mid-surface.**
  `reduction::nodesNearPlanes` at anything like its 1e-9 default keeps only the node
  of each through-thickness pair that happened to land on the plane, and the plating's
  normal leans out of the transverse direction wherever the hull is not
  parallel-sided — up to `t/2`, 8 mm on this ship. Half an interface is not a cut, it
  is a hinge. Measured and asserted on a section forward of amidships, where the
  node-based rule returns strictly fewer nodes than the mid-surface one.
- **The node numbering is a hundredfold, not a rounding.** `solidshell::solveStatic`
  numbers its free DOF in the mesh's own order and has no renumbering pass, unlike
  `reduction::Substructure`. Which ordering is right is a property of the section — a
  hold has twelve stations along x and hundreds round the girth, a slender box is the
  other way round — so three candidates are built and the narrowest kept, including
  `reduction::bandwidthReducingOrder`, which is exposed for this. Getting it wrong
  cost a 6 240-DOF solve **9.53 s instead of 0.06**, and took the ferry hold's
  half-bandwidth from 146 to 1 382 — 0.14 s of banded factorisation against 5.3.

#### What mutation testing found

40 mutants, each a single plausible edit, run against the section suite alone with a
per-mutant timeout. The first pass killed 27 and every one of the nine survivors was
a real hole or an equivalent mutant; **39 are killed now and the one that survives is
argued equivalent.** The ones worth recording:

- **A tie weight computed against the wrong thickness** — the 47 mm defect above.
  It was found by a mutant that swapped the meshed thickness for the member's own
  nominal one, which changed *nothing measurable*, which is what made it worth
  looking at. The invariant that kills it now is that a fibre's tied point sits at
  the offset the fibre records, from its own pair's mid-surface: exact to **1.8e-13 m**
  where it had been 4.7e-2.
- **The weld radius and the bucket probe mask each other.** Comparing the squared
  distance against the tolerance rather than against its square makes the weld radius
  a millimetre — and the mutant survived a test that nudged a panel by ten microns,
  because the bucket probe never looks that far. It dies at **three** microns, inside
  the probe's reach and outside the tolerance, where only the distance test can
  refuse it.
- **Reversing whole plates proves nothing about orientation.** A test that flipped
  every other *panel* in the list flipped whole surfaces, so their normals still all
  agreed and a mutant that dropped the orientation walk entirely survived. It dies
  when alternate *bays within* each surface are flipped, which is what
  `makeStructuralMesh` actually does when it mirrors the starboard side.
- **`spanningComponents` and `floatingComponents` are not the same test.** A
  component reaching one cut plane and stopping is restrained, so nothing refuses it,
  and it carries none of the section's axial or bending stiffness. Only a case with
  such a component tells the two apart.
- **The junction census must exclude a surface's own plating**, or it degenerates
  into "this edge is free" and the ferry's 370 m means nothing. The control is a
  single flat plate, whose two long edges are free and whose junction length must be
  zero.
- The one **equivalent** mutant halves the welder's bucket size. It cannot matter,
  because the probe's reach is computed from the tolerance and the cell rather than
  hard-coded — which is the fix `zone.cpp` made after exactly this mutant survived
  there.

#### What it does not do

1. **It cannot weld a junction** — it *ties* one instead. See the section below.
   What remains untied is the junction nodes on a cut plane, 60 m of 370 on an
   eleven-bay hold, because an interface degree of freedom is prescribed and cannot
   also be derived.
2. **It takes one material.** A section spanning two is meshed entirely as the first
   and says so; on this ship the weather deck's mild steel differs from AH36 in yield
   alone, so nothing the reduction reads is affected.
2b. **A member lying on the forward cut plane goes to the next section.** A frame or
   a deck beam has no extent along x, so it sits on a station, and a station is a cut
   plane for the sections on both sides of it — taken by both, a chain carries it
   twice. The rule is half-open at the aft plane, the same one the transverse
   *plating* already used, and `Section::membersOnForwardPlane` counts what it left
   behind (70 on a two-bay ferry section). A section standing alone is therefore short
   one ring of frames, which is worth nothing to the hull girder — an athwartships
   member carries no longitudinal stress and `hullGirderSection` drops it too — and
   worth 10.1% of a two-section chain's stiffener mass, which is what the rule exists
   for.
3. **A member is only picked up where it runs along a panel seam.** A girder
   positioned off the longitudinal spacing is invisible to a mesh built from panels,
   which is 0.075 m² and 4% of this ferry's section area. It is reported rather than
   quietly absorbed.
4. **A member touching two surfaces is refused rather than projected.** The
   eccentricity is measured along the plating's thickness direction, and a member
   reaching a corner picks up node pairs whose directions are 90° apart. Twelve of
   this ship's members in a two-bay section; the alternative is a plausible wrong
   second moment.

#### Reaching the whole ship: the "inverted element" was a triangle

This section was written, validated and mutation-tested on **a quarter of this ship**,
and nobody noticed until the reach was asked for directly. Of the 49 two-bay windows,
21 meshed, reduced and solved; the rest reported "an element came out inverted or
degenerate" and were refused by `applyBeamLoad` and `Substructure`.

**And they were never a range.** The 21 were five islands — x = −43.2…−38.4,
−36…−31.2, −28.8…−9.6, −7.2…19.2 and 21.6…28.8 — so 62.4 m of 120 m worked in total
while the longest *unbroken* run was 26.4 m: exactly the eleven-bay hold every figure
above was measured on. The two windows containing the bulkhead at x = −8 failed while
sitting in the middle of the best island, which is the clue that the cause was never
"the ends are curved".

**The guess was curvature approaching the plate thickness. The measurement says
otherwise: there is not one negative Jacobian anywhere on this ship.** Over all 49
two-bay windows every non-positive determinant is exactly zero, and every one of
them is on a *collapsed* element:

- `makeStructuralMesh` emits **degenerate `PlatePanel`s** — quads with two coincident
  corners, i.e. triangles. 90 bulkhead panels and 76 deck panels of 8 900; **no shell
  panels at all**. A bulkhead grid running into the centreline or the keel; a deck
  laid on fixed |y| lines clipped to a hull that narrows past them. 39.1 m² of
  12 802.9 — **0.305%** of her plating.
- Extruding one gives a **collapsed hexahedron**: a triangular prism written in eight
  nodes, two of them the same node. A covariant base vector vanishes on the closed
  edge, so `det J` is exactly zero *there and only there*.
- `solidshell::smallestJacobian` samples the eight **corners**. The element is
  integrated at the centre and the 2×2×2 Gauss points, and at every one of those it is
  sound: the collapsed elements' worst Gauss determinant on the ferry is 5.6e-6
  against 2.7e-5 for the worst *sound* element. No closer to singular than the mesh
  already was.

`solidshell::ElementShape` separates the two. A section is refused when a corner that
nothing coincides with has gone non-positive, **or** when the quadrature has — not
when a wedge's closed edge reads zero. Both halves are load-bearing and both have a
negative control in `tests/test_solid_shell.cpp`, including a wedge folded through its
thickness at one corner, which the quadrature alone would accept.

**What a wedge is worth**, measured on the box girder with *every* panel triangulated
— the worst case, against 1.9% of elements on the ferry:

| subdivision | EA/exact | EI/exact | GJ/Bredt, tied | (quads, tied) |
|---|---|---|---|---|
| 1 | 1.00000000 | 1.1370088 | 1.2015 | 1.0986 |
| 2 | 1.00000000 | 1.0068191 | 1.0516 | 1.0512 |
| 3 | 1.00000000 | 1.0021759 | 1.0315 | 1.0367 |
| 4 | 1.00000000 | 1.0011108 | 1.0126 | 1.0297 |

`EA` is exact at every refinement — a collapsed element still passes the patch test —
and `EI`'s 14% converges away twentyfold on the first refinement, so it is a
discretisation error and not the formulation. The test asserts the *convergence*,
because a tolerance on the coarse value would pass on an element that was simply
wrong by a constant.

**The reach**, from `tools/section_probe --scan=2`, which `verify.sh full` now runs:

| | before | after |
|---|---|---|
| two-bay windows that mesh, reduce and solve | 21 / 49 | **49 / 49** |
| of those, a single connected piece | 21 / 49 | 46 / 49 |
| union of hull inside a working window | 62.4 m (52.0%) | **120.0 m (100%)** |
| longest *unbroken* run of working windows | 26.4 m, x = −7.2…19.2 | **120.0 m, bow to stern** |
| the whole 120 m as one section | refused | 8 900 elements, **one component** |

**The old range was never the contiguous 43 m this document used to quote.** The
working windows were an archipelago: x = −43.2…−38.4 and −36…−31.2 worked, so did
21.6…28.8, and the two windows containing the bulkhead at x = −8 did *not*, even
though they sit in the middle of the range. What a chain needs is the unbroken run,
and that was 26.4 m — eleven bays, which is exactly the hold every measurement in
this section was made on.

One window regressed, and it is worth naming rather than netting off: x = −43.2…−38.4
used to come out in one piece and now comes out in two, because `junctionWeightLimit`
refuses a tie there that was closing a 22.4 mm gap with a weight of 2.10. It still
meshes, reduces and solves.

Two things it did not fix, both now reported rather than latent:

1. **A tie can take more of a slave's mass than the slave has.** The junction search
   admits a master whose mid-surface is within `junctionTolerance` — 25 mm, an
   *absolute* figure — while the through-thickness split is relative to the master's
   own thickness. On 10 mm plating a slave 21.5 mm off the mid-surface is admitted and
   asks for a weight of 2.65, putting −1.65 of the slave's steel on one master face;
   `Substructure` then refuses the **whole section** for a non-positive nodal mass.
   `SectionParams::junctionWeightLimit` refuses the junction instead — one full share
   is the limit, and on this ship that separates cleanly, every junction being at or
   below 1.885 or at 2.10 and above. It costs three windows their single-component
   status, which is why the two rows above differ. The fix that needs no limit is a
   `junctionTolerance` scaled off the plating, which changes the junction *census* as
   well as the tie.
2. **`hullGirderSection` sampled exactly on a frame station loses 76% of the ship**,
   and it is the Tier-0 reference this whole section is compared against. Measured:
   0.42932 m² at x = 19.2, 21.6 and 24.0 against 1.80133 at 18.8 and 19.6 — the
   plating either side of the plane fails the half-open `straddles` test in
   `sectionElements` on a floating-point knife edge and only the longitudinals
   survive. Not every station does it (16.8 is fine), which is what says it is
   representation rather than geometry. `section_probe --scan` now samples Tier 0 at
   the centre of a bay and reports both areas rather than their ratio; fixing it
   belongs with `girder.hpp` and moves published figures.
3. **The thickness-seam rule cost five times more at the ends than amidships —
   *fixed*, see §The halo.** `nodeThickness` was an area-weighted mean over the
   sub-quads *inside* the section, so a station where a strake steps handed one
   thickness to the section aft of it, another to the one forward, and the mean to a
   section spanning it. A member run stops at a thickness change, so a spanning
   section stopped runs its two halves never saw and dropped the seam node's run of
   one. Cutting a window in two conserved the stiffener steel **exactly** amidships,
   and did not at x = −24…−19.2 or at the bow shoulder. The direction is the spanning
   section's: the one piece carried 1.9% less than its two halves at x = −24…−19.2 —
   inside the range that always worked, so this is older than the reach — and
   **25.2%** less at x = 40.8…45.6. It is the nodal-thickness twin of the
   nodal-normal problem in §*A ship: a chain of sections*, and it had the same fix: a
   halo.

### The junction tie — **implemented**

`solidshell::Mpc` and `solidshell::DofExpansion` in `engine/sim/solid_shell.{hpp,cpp}`,
carried by `reduction::Attachment::constrained`, built by `section::buildSection`.
Checked by `tests/test_solid_shell.cpp`, `tests/test_reduction.cpp` and
`tests/test_section.cpp`. **The decks are attached to the shell.**

A free-edge node lying on another surface is tied to the point of that surface it
lands in — bilinear in the master face's two in-plane coordinates, and through its
thickness by the weight `constraint.hpp` §1 derives:

```
u[slave] = Σ_a N_a(ξ, η) · [ (1−w) u[bottom_a] + w u[top_a] ]
```

Eight masters, three constraints per extruded node, and the slave keeps **no unknown
of its own**: `solveStatic` and `Substructure` both scatter through the
transformation, so the matrix they factor *is* `TᵀKT`. No penalty stiffness, no
Lagrange multiplier, no tuned parameter — which matters, because a stiffness no
measurement sets is precisely what `solid_shell.hpp` rejects hourglass control for.

#### Why it could not be a `DofBlock`, which is the load-bearing finding

`constraint.hpp` already eliminates an eccentric stiffener's fibre with a `Tie` and
hands the result over as a `DofBlock` of extra stiffness. The obvious move was to do
the same here, and it does not work. **A fibre's endpoints are not mesh nodes**: they
have no rows of their own, so `TᵀKT` over the masters is the whole of what the fibre
contributes and adding a block is exact. A junction ties a node that *is* a mesh
node, with elements of its own — eliminating it means rewriting rows that already
exist, and no amount of added stiffness redirects a row. So the constraint had to
become something the assembler understands rather than something the assembler is
handed.

#### It is not a weld, and the difference is measurable

A weld merges two node pairs into one, which has one thickness direction where two
plates at an angle need two, and pays for it in steel. A tie moves no node: the same
nodes in the same places, the same elements, the same mass, still four surfaces on
the box — and one connected component. On the box girder, at one element per panel:

| corners | components | `EA` | `EI` | `GJ` / Bredt | band |
|---|---|---|---|---|---|
| cut | 4 | 1.000000 | 1.000000 | 0.083 | 62 |
| welded | 1 | 0.905788 | 0.911115 | 0.966 | 104 |
| **tied** | **1** | **1.000000** | **1.002158** | **1.099** | **236** |

The tied mesh is the only one that is *both* exact in `EA` and closed in torsion.
`1.099` is not an overshoot of Bredt: a closed cell does not stop having the open
section's own `Σ s t³/3`, and `1.000 + 0.083` is what the two together predict.

`EI` costs 0.216% at one element per panel, and it is a **consistency error rather
than the formulation**: at a butt corner the tied node sits half a plate thickness
past the end of the other plate's mid-surface, so the tie extrapolates a bilinear
approximation to a quadratic transverse contraction. Refining gives 2.16e-3, 2.85e-4,
8.62e-5 — twenty-five fold over a threefold refinement. Clamping the tie onto the
face instead would trade that for an `O(t)` error refinement cannot reach, which is
why the overshoot is bounded rather than forbidden.

#### On the ferry, where the fixed-interface frequency is the instrument

The hold between x = −7.2 m and 19.2 m, one element per panel:

| | components | free edge joined | `GJ` | first fixed mode |
|---|---|---|---|---|
| untied | 7 | 0 of 370 m | 3.6164e12 | 0.7785 Hz |
| tied | **1** | **310 of 370 m** | **5.2387e12** | **2.3026 Hz** |

0.7785 Hz is the decks' own frequency to four figures — adding the shell changed
nothing, because in the model it was not attached. 2.3026 Hz is above the decks'
0.7785 *and* the shell's 1.6999, which is what a joined structure does. `GJ` moves
+44.9%.

`A_eff` moves +0.19% and `I_eff` +0.12%, the same local consistency error the box
shows. **Neither could have told you the junction was open** — which is the whole
point of the section above.

#### What it cost the band, measured rather than assumed

A tie couples nodes that no element edge joins, and the node numbering has to know:

| ferry hold | DOF half-bandwidth | banded static solve |
|---|---|---|
| untied | 146 | 0.16 s |
| tied, ordering scored on the element graph | 10 769 | (effectively dense) |
| tied, ordering scored on the tied graph | 1 520 | 5.34 s |

The mesher compares three candidate orderings and now scores them on the *tied*
graph. That is worth a factor of seven. What is left is real: joining a section's
decks to its shell closes the cross-section into a tube with internal webs, and a
wavefront that used to sweep one flat sheet has to cross a deck. `Substructure`'s
interior band goes 125 → 455 and its Craig-Bampton 6.0 s → 57.7 s; at subdivision 2
the banded static solve goes 1.1 s → 149 s, which is why `tools/section_probe` runs
its resolution sweep untied.

#### What it refuses

- **Chains.** A node whose master is itself a slave is left untied and counted.
  Composing a chain silently is a modelling error that assembles.
- **Interface degrees of freedom.** A reduction keeps its interface *exactly* and
  `applyBeamLoad` prescribes it; a degree of freedom cannot be both prescribed and
  derived. `Substructure` refuses a constrained boundary DOF outright.
- **A face overshoot past `junctionOvershoot`.** The weights are a partition of unity
  but a `d` overshoot puts `−d/2` on a master, and a negative weight is a negative
  share of the slave's steel in `TᵀMT`. `Substructure`'s positive-nodal-mass check is
  the backstop. Unlike the eccentric fibre there is nothing to give up: a fibre
  *extrapolates* and must abandon its first moment, while a junction interpolates and
  keeps it.
- **A degree of freedom both pinned and tied.** Two claims on one unknown, refused by
  `solveStatic`.

#### What mutation testing said

The recurring shape in this repository is an error that cancels when asked globally,
and a junction tie has exactly that shape — so the tests that carry the weight ask
about the tie *alone*:

- The constraint set is checked as an algebraic object before anything is solved
  through it: `Σ w = 1` to 1.1e-16 (a rigid translation is reproduced exactly) and
  `Σ w X_master = X_slave` to 8.9e-16 m (so is a rigid rotation). A weight dropped,
  mis-scaled or attached to the wrong axis dies here and nowhere else — an energy
  comparison would report it as a small stiffness change.
- The constrained patch test: a plate whose boundary is prescribed to a linear field,
  with an interior node eliminated in favour of the two either side. A linear field
  satisfies that constraint identically, so the answer must be **exactly** the
  unconstrained one — 6.5e-18 m of 3.0e-3. Its vacuity guard is a constraint the exact
  field does *not* satisfy, which must move the answer.
- The substructure's interior must be in equilibrium on the field `solveStatic` found
  through the same constraints — two independent assemblies of one operator — while
  the unconstrained substructure is not, on the same field.
- The band the mesher reports is rebuilt from the mesh and the constraints in the
  test, because a reported band that does not cover the assembly is a decoration.

Thirty-two mutants, thirty killed. The two that survive are equivalent on every
input this repository has and are recorded as such in `section.hpp` §5: no geometry
here produces a **chain**, and the projection's Gauss-Newton lands in one step on
every face here because they are all flat enough.

**The run also found a defect older than any of this.** `Substructure::totalMass()`
and `stiffnessTimes()` read past the end of arrays the constructor never sized when
the substructure *refused* — reachable from an inverted element long before a
constraint existed, and reported to the caller who skips `ready()` as a read of
address zero. It surfaced as a segmentation fault three tests downstream of the
mutant that provoked it, which is also how the mutation harness was found to be
wrong: its first version counted `FAIL` lines only, and a crash prints none, so it
scored **eight false survivors**. A mutation runner has to look at the exit code.

### A ship: a chain of sections — **implemented**

`reduction::assemble(parts, joints)` and `section::buildChain` in
`engine/sim/{reduction,section}.{hpp,cpp}`, checked by `tests/test_reduction.cpp` and
`tests/test_section.cpp` — where every box-girder figure below is an assertion — and
measured at ship scale by `tools/section_probe --chain=N`, which is where every ferry
figure comes from and how it is re-checked rather than re-quoted. **The whole length,
cut into N pieces, each reduced once, assembled into one model.**

```cpp
struct Component { const Substructure*; const Reduction*; };
struct Joint     { int a, b; InterfaceMap map; };
std::vector<Joint> matchComponents(const std::vector<Component>&, double tolerance);
std::vector<Joint> matchNeighbours(const std::vector<Component>&, double tolerance);
Assembly assemble(const std::vector<Component>&, const std::vector<Joint>&);
Assembly assemble(const Assembly&, const Component&, const InterfaceMap&);  // the fold
int      assembledComponents(const Assembly&);

Chain           buildChain(const StructuralMesh&, const ChainParams&);
BeamResponse    applyBeamLoad(const Chain&, const BeamLoad&);
TorsionResponse applyTwist(const Chain&, double twist, double reference);
```

#### The generalisation is a union-find, not the identity

The note this replaced predicted that carrying each assembled boundary row's identity
through would be the one change needed, "so that the same position match can be run
against it". Half of that held. The identity **is** worth carrying —
`Assembly::boundaryPoint` does, and a caller needs it to name an assembled row
geometrically, which is what prescribing a plane-sections field on a chain's two end
cuts is. But it is not what generalises the assembly. What stopped at two components
was the *bookkeeping*: which assembled row a component's boundary DOF belongs on. Two
components settle that with one `InterfaceMap`; N settle it with a **union-find over
every pairwise map**, because a DOF shared by A-and-B and by B-and-C is one row and
no pairwise map says so on its own. Every map is still expressed in the
substructures' boundary DOF, which the caller has, so that route needs no identity at
all.

Both routes are built, because "one is better" is a measurement rather than an
argument. The fold — carry the identity, match against the assembly, add one
component — **produces the same matrices and the same `from` maps to the last bit**,
and re-scatters a matrix that grows every fold: O(N³n_b²) against O(N²n_b²). Measured
on the box girder at ×0.9 for one component, ×1.7 for two, ×3.2 for four and ×5.1 to
×6.0 for eight — the spread is the timer, the growth is the shape. It is kept because
a caller that acquires components one at a time has no N to assemble at once, and
because a route that agrees bit for bit is a stronger test of the index arithmetic
than any tolerance.

#### It is exact at zero modes, and that is a sharper test than a convergence

The assembled boundary is the union of the cut planes; Guyan condensation is exact
there for any load with no interior load; the scatter-add approximates nothing. So a
chain and the same length in one piece are the *same static model*, not a converging
approximation to it — provided they are the same structure, which means junctions off
on both, for the reason below.

| like for like, untied | chain / one piece |
|---|---|
| box girder, 4 sections, `EA` | −9.9e-13 |
| box girder, 4 sections, `EI` | +1.4e-10 |
| box girder, 4 sections, `GJ` | −1.3e-10 |
| ferry x = −7.2…2.4, 2 sections, `EA` / z_na / `EI` / `GJ` | +2.5e-10 / +1.9e-10 / +4.8e-10 / +3.9e-12 |
| ferry x = −7.2…16.8, 5 sections, all four | ≤ 6.3e-10 |

Dynamically it converges from above, as one component does: the chain's lowest
frequency with both end planes held is 0.844705 Hz at zero modes and 0.843943 at two,
against the monolith's own fixed-interface 0.843941 Hz computed through none of the
assembly.

#### What a cut plane costs, which `EA` cannot see

**A junction node on a cut plane is an interface degree of freedom, and one of those
is prescribed rather than derived — so it cannot also be tied.** Cutting a length into
N pieces turns N−1 interior stations into interfaces and unties every junction on
them. The previous note hoped an interior cut plane in a chain might stop being an
interface; it does not, and the effect runs the other way:

| box girder, tied | 1 section | 2 | 4 | 8 |
|---|---|---|---|---|
| junction edge joined, of 64.0 m | 58.0 | 52.0 | 40.0 | 16.0 |
| `EA` against the closed form | 1.2e-12 | 8.5e-14 | 1.5e-12 | 1.1e-12 |
| `GJ` against one piece | −2.6e-12 | −3.29% | −9.24% | **−19.20%** |
| `GJ` / Bredt | 1.0986 | 1.0625 | 0.9972 | 0.8876 |

On the ferry the same shape: x = −7.2…2.4 as two sections ties 9.6 m of 134.4 against
the monolith's 72.0 and loses **3.5%** of `GJ`; x = −7.2…16.8 as five ties 28.8 m of
336.0 and loses **8.9%**. `EA` moves 2.5e-4 and `EI` 4.5e-5 in the same runs.

**That row of `EA` figures is the whole reason this section is written the way it
is.** A chain of eight sections that ties nothing to anything reproduces the closed
form to twelve figures, because prescribing plane sections at the ends makes every
longitudinal strip carry σ = Eε whatever it is joined to. Any validation of an
assembly that stops at `EA` and `EI` scores a do-nothing implementation best. It is
the same trap §*The validation that proves nothing* records for the mesher, and it
caught this work too — the first version of the ferry comparison agreed to 1e-10 in
`EA` while `GJ` was 33% out.

The fix that would close it is a tie whose masters all lie *in* the cut plane, which
both sections would then derive identically from shared boundary DOF. That is a line
tie through the plating thickness rather than the bilinear face tie the junction work
built, so it is new machinery rather than a parameter. Until it exists the rule is:
**cut a ship into as few pieces as the interface cost allows.**

#### The interface was not coincident by construction — *fixed*, see §The halo

Two sections cut on the same frame station agree on every *mid-surface* point of it —
they come from the same panel corners. They did not agree on the mesh nodes, because a
node is the mid-surface offset by t/2 along the **nodal normal** and `buildSection`
averaged that normal over the sub-quads inside the section: aft of the plane for one,
forward for the other. Prismatic hull, same vector, exact coincidence. Shoulder, not:

| plane x | −45.6 | −33.6 | −21.6 | −9.6 | 0 | 12.0 | 24.0 | 36.0 | 48.0 |
|---|---|---|---|---|---|---|---|---|---|
| worst gap, m | 6.4e-4 | 3.5e-4 | 3.4e-6 | 0 | 0 | 0 | 1.9e-4 | 7.1e-4 | 1.0e-3 |
| nodes past 1e-9 | 116 | 112 | 112 | 0 | 0 | 0 | 112 | 116 | 120 |

That row was measured when the mesher refused two-bay sections outside about
x = −26.4…19.2, so most of it was unreachable and only **x = −21.6** could be
checked. `section.hpp` §7 removed that limit — every plane in the row is reachable
now — and the figures below are still the ones taken at x = −21.6, because that is
where they were measured and re-measuring them elsewhere is a different claim.
The gap there is 3.4 µm: three orders of magnitude below
the plating and three above `matchBoundaries`' 1e-9 default. At the default **336 of
that plane's 1 170 boundary DOF found no partner**, and a chain assembled out of them,
solved, and was torn along 29% of the cut. `Chain::unmatched` counts them and
`ready()` refuses. At 1e-5 the whole plane matches and the chain reproduces the
monolith to 3.1e-7 in `EA`, 3.0e-7 in `EI` and 1.1e-6 in `GJ` — so joining two nodes
3.4 µm apart is a legitimate answer, and the reason to set the tolerance was that the
default silently did something worse.

**The whole row is now zero and the tolerance is not needed at all.** The table is
kept because `SectionParams::halo = false` still reproduces every figure in it, and it
is the negative control the fix below is measured against.

### The halo — **implemented**

`SectionParams::halo`, in `engine/sim/section.cpp`. Checked by
`tests/test_section.cpp` and by `tools/section_probe --invariance=2`, which is in the
`full` gate.

**A node is the mid-surface offset by t/2 along the nodal normal, and both of those
were averages over the sub-quads inside the section** — so both were statements about
the *cut*, and the two defects above are one defect. The fix averages both over one
bay of plating **beyond** each cut plane and meshes only what is inside.

A bay is not a distance, and it is not a corner test either: the halo is every panel
whose extent along x reaches the window, and the **weld** decides what that is worth.
A halo panel's grid points go through the same `Welder` as the section's own, so one
that lands on a node of an inside sub-quad joins that node's averages and one that
does not contributes to nodes nobody keeps and is compacted away. Nothing has to be
predicted, and there is no second notion of "the same point" — the trapdoor
`sectionElements` fell through. The bound is exact, which is why it is allowed to be a
bound: a panel reaching a node of an inside sub-quad has a grid point at that node's
x, which lies between the planes.

> **This was a corner test first, and mutation testing is what said it should not be.**
> "Shares a welded corner with a panel inside" is the same set only when every node is
> a panel corner; at `subdivision > 1` a node can sit part way along an edge, where a
> panel butting mid-edge reaches it and shares no corner. Three mutants survived on
> that predicate — widening it, an off-by-one in it, a second weld tolerance inside it
> — all equivalent on a conforming hull and none of them equivalent in general.
> Deleting it is cheaper than assuming it, in code and in time.

It is **exact rather than close**. The halo's panels join the candidate list in
ascending panel index, so a node reached by the same panels from either side of a cut
accumulates the same terms in the same order and comes out at the same double. The
interface test is `worstGap == 0`, not `worstGap < tol`.

| ferry, 4.8 m windows either side of a station | before | after |
|---|---|---|
| stations where a node on the plane moves | 32 of 47 | **0 of 47** |
| the furthest it moved | 2.5e−3 m | **0** |
| unmatched boundary DOF at x = −21.6, at 1e-9 | 336 / 1 170 | **0 / 1 170** |
| stations whose stiffener steel depends on the cut | 32 of 47 | **0 of 47** |
| the worst of that | 12.2% | 3.9e−14 |
| stiffener steel, whole ship in 1-bay windows | 517 629 kg | 501 263 kg |
| …in 2-bay windows | 501 427 kg | 501 263 kg |
| …in 5-bay windows | 501 206 kg | 501 263 kg |
| …as one 120 m piece | 501 263 kg | 501 263 kg |

And the solved consequence on the same plane, from `section_probe --chain=2
--from=-26.4 --to=-16.8` — a chain of two against the same length in one piece, each
at the tolerance it can actually use:

| | before, at 1e-5 | after, at the 1e-9 default |
|---|---|---|
| `EA` | +3.081e−7 | **+3.217e−10** |
| z_na | +2.323e−7 | **+5.363e−10** |
| `EI` | +3.013e−7 | **+4.966e−9** |
| `GJ` | −1.051e−6 | **+4.161e−12** |

The tolerance was never free: merging two DOF 3.4 µm apart is a legitimate answer and
costs about what the gap is worth. Not needing it is worth three orders in `EA` and
six in `GJ` — and at the default the same chain used to be refused outright.

**The stiffener block above is what says the halo is right rather than merely
consistent.** Cut the ship four ways and it gives one answer, and that answer is the one the whole ship
in one piece gives — the only one of the four that never had a cut to depend on.
Before, a 1-bay partition carried 3.3% *more* steel than the monolith, because a
section that cannot see across its own cut does not know there is a strake seam there
and keeps a run the ship does not have. The steel the bow shoulder appears to lose —
18 270 kg over two 2.4 m windows becoming 13 669 — is that steel, and it was never
there. What remains is the thickness-seam rule's own cost, unchanged and now uniform:
501 263 kg of the members' own 631 451 kg, **79.4%**, whose fix is still the
per-station thickness in `constraint::SeamRun`.

**Blast radius: one deletion.** The halo's sub-quads are dropped, and the nodes only
they reached compacted out, immediately after the two averages are formed and before
anything else in `buildSection` runs — so the free-edge count, the junction search,
the surface census, the component walk and the member runs all see the sub-quad list
they always saw.

**What it does not move.** Amidships the halo's sub-quads carry the same normal and
thickness as the ones inside, so a node moves 1.7e−18 m — a few ulps of t/2, fourteen
orders below the millimetre it removes at the ends. The reach stays 49 of 49 windows
and 46 in one piece; the hold's first fixed-interface mode (0.7785 untied, 2.3026
tied), DOF half-bandwidth (146 and 1 520), component count (7 and 1) and 309.6 m of
tied edge are unchanged to every digit, and its `GJ` moves 3e−5. A chain of four still
reproduces the same length in one piece to 9.9e−13 in `EA`, 1.4e−10 in `EI` and
1.3e−10 in `GJ` on the box, all three unmoved. Meshing the whole ship costs 0.44 s
against 0.45 s, and 50 one-bay windows 0.15 s against 0.11 s.

**What it does move.** The hold's forward plane at x = 19.2 is a strake seam, so a
hold-length section can now see one, and four members lying wholly in the bay against
it end up in runs of two stations of different thickness — no run at all. `A_eff` goes
1.73394 → 1.73266 and `I_eff` 43.913 → 43.868, both about a tenth of a percent. That
is the *consistent* answer rather than a worse one: it is what the same plating gives
in any other window containing x = 19.2. A section's plate *mass* is likewise no
longer exactly Σ A_p t_p over the panels it owns, because the thickness taper now
reaches the cut plane too: a window at a strake step hands its neighbour parts per
million, and the sum over any partition of the ship is unchanged at 1 084 106 kg.

#### Two things a chain got wrong that a total would never have shown

**A ring of frames counted twice.** A member with no extent along x — a frame, a deck
beam — lies *on* a station, and a station is a cut plane for the sections on both
sides of it. Both built it. Two four-bay sections of the ferry came to **10.1% more
stiffener mass** than the same length in one piece: 608 fibres, and their stiffness
and mass land on shared interface DOF where nothing but a total would see it. The
transverse *plating* already had the half-open rule that fixes it — a bulkhead is
taken at the aft plane and not the forward one — and the members simply never got it.
`tests/test_section.cpp` now asserts additivity: element for element, 4.4e-16 in
plate mass, 1.5e-14 in stiffener mass, 2 400 + 2 400 = 4 800 fibres, 460 + 460 = 920
members. Nothing in the suite noticed a 22% change in a section's own `memberMass`
when the rule landed, which is what a missing identity looks like.

**Two component counts, and neither answers the other.**
`reduction::assembledComponents` works on the reduced model, where a component's pair
is dense over its own DOF whether the mesh behind it is one piece or seven — so it
sees a chain that failed to join at a cut plane and is blind to a section whose decks
are not tied to its shell. `Section::components` is the other way round. `Chain`
carries both, and `Chain::components` joins them: a mesh component of one section and
a mesh component of the next are the same piece of ship when they share an assembled
boundary row. **That is the one that must be 1.** It is also what
`applyBeamLoad(Chain)` restrains — three rigid-body restraints per piece, twelve on
the untied box and twenty-one on an untied ferry section, because prescribing `u_x`
alone leaves each piece free in y, z and rotation about x.

#### What it costs, measured

The interface is the whole cost and it is set by the cross-section, not by the length:
a transverse cut of the ferry is **1 170 boundary DOF** at `subdivision = 1`, so a
chain of N is a *dense* (N+1) × 1 170 reduced model.

| ferry x = −7.2…16.8, ten bays | 5 sections | one piece |
|---|---|---|
| assembled / boundary DOF | 7 020 | 2 340 |
| reduce, untied / tied | 1.62 s / 9.13 s | (57.7 s tied, eleven bays) |
| assemble | 0.23 s | — |
| one static solve, untied / tied | 24.4 s dense | 0.13 s / 3.60 s banded |

So cutting a ship into more pieces makes every piece cheaper to reduce and the
assembly quadratically dearer to *solve*, and that is the trade a caller is making
whether or not it is stated. It is the right trade only because a reduction is built
once and used many times; for a single static answer the monolithic banded solve wins
by two orders of magnitude.

#### What mutation testing found

47 mutants, each a single plausible edit to `reduction.cpp` or `section.cpp`, run
against the reduction, section, coupling, constraint and solid-shell suites with a
per-mutant timeout and the verdict taken from the **exit code** rather than from the
`FAIL` lines — the harness that counted lines scored eight false survivors here once
already, on a mutant that segfaulted three tests downstream of the one it broke.

The first pass killed 33 of 46 and left 12; the final pass kills **42 of 47** with
four argued equivalent and one a deliberate no-op control. Of the twelve survivors
**seven were holes in the suite**, and they cluster in one place:

- **`matchBoundaries` was never asked about on its own.** Deleting the axis check,
  deleting the "this DOF of B is taken" mark and deleting the `break` that keeps the
  *first* match all survived. They survive **together**: a substructure's boundary
  DOF come out ordered by node and then by axis, so a greedy first-match walks the
  two lists in step and lands on the right partner for the wrong reason, and each of
  the three mutants is masked by the other two. Two coincident nodes on the same axis
  — an unwelded seam — separates them, and so does a list whose axes are in the other
  order. Asking the primitive directly, on lists built here, kills all three.
- **A refusal is not enough without its reason**, for the second time in this
  document. Deleting the "this assembly carries no identity" guard from
  `matchBoundaries(Assembly, Substructure)` leaves the match running against an empty
  list, which raises the *other* complaint — "these share no boundary DOF" — so a
  test that only asked whether `problems` was non-empty passed on it.
- **An interface map shorter than the boundary it claims to describe** was accepted.
  That is what a component whose reduction fell back to empty produces, and it joins
  the DOF the map reaches and leaves the rest as two unknowns.
- **`Assembly::worstMergedGap` was never non-zero**, because nothing else in the
  suite merges two DOF that are not exactly coincident. Replacing `max` with `min`
  survived. The ferry shoulder at a 1e-5 tolerance is the case with something in it.
- **Two routes mean two copies of the same guard**, and testing one tests one. The
  fold's axis check survived everything the N-way one killed. Closing it also turned
  up a real asymmetry: the fold *counted* a crossed axis and said nothing about it,
  where the N-way route complained.

The remaining survivors are argued **equivalent**, and saying why is more useful than
a score:

- **The bounds check in `boundaryIdentity`** (`3n+2 >= size` against `>`) can only
  differ when `3n+2 == 3N` exactly, which has no integer solution.
- **Transposing the scatter's read** cannot matter: `Reduction::stiffness` and `mass`
  are symmetric to the last bit by construction — the analytically symmetric blocks
  are computed on one triangle and mirrored — which is stated in `reduction.hpp` and
  is what this mutant measures.
- **Either root choice in the union-find** is a valid merge.
- **`matchComponents` in place of `matchNeighbours`** in `buildChain` is equivalent in
  *result* on a chain and not in cost: two sections two bays apart share no boundary
  DOF, so the extra pairs produce no joints and only a 1 170 × 1 170 distance sweep
  each. It is kept as the deliberate difference between the two calls.
- One is the **no-op control**, because a harness that reports everything killed is
  reporting nothing.

### Tier-1 to Tier-2 coupling — **implemented**

`engine/sim/coupling.{hpp,cpp}`, checked by `tests/test_coupling.cpp`. The piece
that made the three tiers one system in the middle as well as at the bottom: **a
zone's perimeter follows the reduced model round it instead of being clamped, and
a torn zone goes back into that model so the surroundings feel it.**

```cpp
Coupling    couple(const Substructure& surroundings, const Reduction&,
                   const Substructure& zone, const Reduction&, double tolerance = 1e-9);
bool        prescribedStaticSolve(const Assembly&, load, held, prescribed, state, problem);
EdgeDrive   edgeDrive(const Coupling&, const Substructure& zone, assembledState);
std::size_t applyEdgeDrive(const EdgeDrive&, HexMesh&);
DamagedMesh withoutTornElements(const zone::Patch&, const zone::Solver&);
```

#### It is exact, not close, and that is the whole design of the tests

The coupling is `matchBoundaries` + `assemble` + a static solve + writing the
interface displacement into `HexMesh::prescribed`. Nothing is approximated at any
step: Guyan static condensation is exact at the interface for *any* load at *any*
mode count, assembly is scatter-add, and a shared boundary DOF **is** one unknown.
So a coupled zone does not approach the monolithic answer — it reproduces it.

Measured on a 1.2 m square plate, 8 × 8 solid-shell elements, clamped rim, with a
gripped punch pushed 200 µm into the middle quarter, against the same plate meshed
and solved in one piece by `solveStatic`:

| | punch reaction | the zone's field vs the whole plate's |
|---|---|---|
| monolithic (the reference) | 1 446 N | — |
| zone clamped, as before | 625 671 N (**433×**) | 1.47e-4 m out, 74% of peak |
| zone free | 0 N exactly | — |
| **coupled**, 0 / 4 / 12 modes | 1 446 N, identical | **1.06e-15 m** of a 2.0e-4 m peak |

Three things in that table are worth more than the agreement.

**The free bound is exactly zero, and as a bracket it is vacuous.** A patch with
nothing holding its edge, handed the same displacement at every gripped punch node,
answers with a rigid translation: no strain, no reaction. So "the coupled force
lies between free and clamped" is satisfied by any positive number and proves
nothing. The bracket that carries content is the **edge-following ratio** — how
much of the punch's travel the perimeter takes up — where both ends are exact
closed forms rather than measurements: a clamped edge follows **0**, a free one
follows **1**, and the plate's own answer is **0.737**. The coupled zone reproduces
that to 1e-9.

**Mode count buys nothing statically, and the test asserts that rather than a
convergence.** This contradicts the natural expectation and it is the same
property `reduction.hpp` §1 records: the reduction starts exact at the interface
and stays there. A test that swept modes looking for improvement would have been
measuring nothing — and a coupling that had quietly *lost* the property would
still improve with modes and still look reasonable, which is why the assertion is
invariance. The sweep carries its own guard that the modes were really kept
(the assembled model grows by 2 × modes), or "adding modes changed nothing" would
be a statement about a sweep that never added any.

**The negative control is the model this replaces.** The clamped zone, measured
with the same instrument, misses the monolithic field by 74% of its peak. Without
that the 1e-15 agreement could have been measuring a well-converged mesh.

The whole chain runs through the **explicit** Tier-2 solver too, not only the
direct one: driven perimeter, prescribed punch, smoothstep ramp, damped to rest —
it settles onto the monolithic field within **5.5e-5 of the peak**, where the
clamped edge it replaces is 1.3e4 times further out.

#### Direction 2: what a torn zone hands back, and why not a thinner ship

`promotion.hpp` reports Tier-2 damage to Tier 0 as an **effective thickness**, and
for a beam that is exactly right — `hullGirderSection` and everything downstream
read a thickness and nothing else. **It is the wrong shape for Tier 1, and the
reason is geometric.** A solid-shell carries its plate thickness in the *positions
of its nodes*: mid-surface at zero, faces at ±t/2. Thinning a zone therefore moves
every node it shares with the structure round it. Measured: a 20% knockdown on
12 mm plating puts every interface node **1.200 mm** out of place — exactly
(t − t′)/2, a closed form — and `matchBoundaries` at the 1e-9 tolerance a coupling
uses then finds **zero** shared DOF. A thinned zone is not a modified component; it
is a different component that no longer fits.

What does transfer, exactly and with no remeshing, is **element deletion**. A torn
element carries nothing, which is a statement about the stiffness assembly rather
than about geometry, so `withoutTornElements` removes it, drops any node left with
no element at all — an orphan is a zero stiffness row, and `reduction.hpp` §3
records that the banded factorisation does *not* reliably catch one — and hands
back a mesh whose surviving interface nodes are where they always were. Measured
on the same plate with a corner element deleted, which orphans exactly its own
corner: **90 of the 96 perimeter DOF still shared, worst gap 5.6e-17 m** — against
zero shared and 1.2 mm for the thickness knockdown.

And it is felt outside. The slit softens the plate by **15.1%**, the surrounding
plating moves **8.7 µm** — 4.4% of the peak — and the recovered surrounding field
matches the monolithic slit plate to **6.6e-15 m**. The control is the field a
one-way coupling would still be reporting: the *intact* surroundings, which are
wrong by 10⁹ times that tolerance.

#### What it does not do yet

1. **Plastic softening short of a tear does not go back.** A reduction is linear
   (§*Tier 1*, what a reduced region cannot do) and the Tier-2 explicit solver
   forms no tangent operator, so between first yield and first tear the
   surroundings are told the zone is stiffer than it is — the unsafe direction.
   Closing it needs a reduction built from a tangent the solver does not form.
2. **A stiffener has to be handed to Tier 1 explicitly, and a caller who does not
   still gets the bare plating.** A zone under `Stiffeners::Modelled` carries its
   members as `constraint::Stiffening` — fibres condensed onto the plating's own
   DOF, with no nodes and no elements of their own — so a `reduction::Substructure`
   built from the mesh and a material alone is **the bare plating, to the last
   bit** (measured: its operator matches an element-only assembly to 3e-16
   relative). One 200 × 10 flat bar across a 0.6 m patch is worth **7.9%** of its
   displacement field, so that was never a rounding, and it is why the fix below
   exists.

   `reduction::Attachment` (`reduction.hpp` §8) closes it: the same
   `solidshell::DofBlock` list `constraint::stiffnessBlocks` already produces, plus
   the nodal mass `constraint::lumpFiberMass` already lumps, folded into the one
   CSR and the one lumped diagonal the elements go into — so the sparsity pattern,
   the bandwidth, the interior renumbering and the reduction carry the member with
   nothing added downstream. Measured on a stiffened plate: the attachment adds
   exactly the profile's area, first moment and second moment about the plate
   mid-surface, to **2.2e-13** of the closed form, checked against
   `scantlings::profileSection` and cross-checked against `stiffenedSection`; the
   reduced model reproduces `solveStatic`'s own stiffened field to **4e-14** of it,
   where the bare plating is twelve orders further away. What it does *not* do is
   find the fibres by itself — the caller assembles the `Attachment`, because a
   `Substructure` is built by the caller and a coupling is handed two that already
   exist.

   Two things that fold into that and are worth having written down. **The
   stiffener's rotary inertia about the seam is given up**, because the fibre mass
   is split equally over the through-thickness pair rather than condensed —
   `constraint.hpp` records why (`T^T M T` is negative on one node of any eccentric
   tie), and the price is measured: about the seam, the model carries 6.8e-4 kg m²
   for that member where its own eccentricity would carry 2.1e-1, a factor of
   **312**. And **a stiffener stiff enough is a node line of the first mode**: on
   the plate above, the 200 × 10 bar raises the first frequency by 307% but its own
   mass moves it by 0.015%, because the first mode has become plate bending beside
   a seam that does not move. The fibre mass is worth 13.7% on the worst of the
   first twelve, and on a 60 × 6 member — one the panel can actually bend with —
   the first frequency lands within 4.5% of what `stiffenedSection` predicts, where
   a smeared panel of identical area and mass is out by 2.5×.

   Separately, and unchanged: `constraint.hpp` gives a fibre no damage variable and
   never deletes it, so "this longitudinal is gone" does not exist to be read —
   `promotion.hpp` records the same gap for Tier 0. That one is a **failure
   criterion for the fibres**, not coupling machinery, and it stays blocked until
   the fibres have one. It has a second edge now: `withoutTornElements` renumbers,
   and an `Attachment` built against the original numbering does not survive that,
   so a torn stiffened zone needs its fibres rebuilt against the damaged mesh.
3. **It is a static coupling.** The assembled solve is `K x = f`; there is no
   co-simulation, no shared time integration and no attempt to run a reduced
   surrounding forward alongside an explicit zone. For a collision that is the
   right first answer — the surroundings are quasi-static over a 0.13 s event —
   but a whipping or slamming load is not, and nothing here would say so.
4. **It couples two components and only two**, inherited from `assemble`. A zone
   between two holds is three.
5. **The interfaces must be node-coincident.** The split is by element, so the two
   meshes share nodes and the coupling is an identity rather than an
   interpolation. Two meshes at different resolutions either side of a cut would
   need `constraint::Tie` — which exists — and a matching pass that finds a face
   rather than a node. Nothing here needs it and nothing here has it.
6. **A tear that reaches the perimeter is handled but not measured.** The zone
   then loses interface nodes, the surroundings keep boundary DOF the zone no
   longer shares, and `couple` counts them rather than refusing — which is
   physically right, the plating really did tear away. There is no measurement of
   what that does to the answer.

#### What mutation testing found

35 mutants, each a single plausible edit. 34 killed; one is argued equivalent.
Nothing here was a defect in the code — which is worth saying plainly, because the
usual result in this repo is that it was — but **eight were holes in the suite**,
and seven of those were guards nothing reached:

- **Three whole preconditions were dead weight.** A reduction that is not the
  zone's, an assembled state of the wrong length, and the boundary conditions
  carried across mesh surgery each had a guard, and removing the guard entirely
  passed the whole suite. The middle one is not cosmetic: every index in a
  `Coupling` is in range for the assembly it was built from, so a short state is a
  read past the end of the *caller's* vector and a plausible drive on the other
  side of it.
- **A wrong answer can be collapsed by `std::unique`.** `carriedInterface` drops
  an interface node the tear removed; a mutant that mapped it to whichever node
  now sits at its place in the ordering survived, because asked over the whole
  perimeter the wrong answer lands on a node that is in the list anyway and the
  duplicate is removed. It dies when the orphan is asked about on its own. This is
  the same shape as the defect the previous agent recorded — two springs to earth
  that both went to zero and agreed perfectly.
- **`applyEdgeDrive` never un-pinning** was untested because every caller in the
  suite cleared the boundary conditions first. It matters exactly where a partly
  matched interface is: the unmatched part has to stay clamped rather than go
  quietly free.
- **A repeated prescribed DOF** takes the last value rather than the sum, and
  nothing named one twice.

The one argued equivalent is the lower clamp in the boundary ramp's smoothstep:
the only call site passes `result_.time + dt`, and `result_.time` starts at zero
and only rises, so the clamp cannot bind. It stays, with the reasoning written
where it is, because `edgeFraction` reads as a pure function of time.

Two mutants are worth recording for what they did *not* need:

- **Taking the drive's target from the rest configuration instead of the meshed
  one** is invisible without a `Preload`, because without one the two arrays are
  equal. The same edit in the same file has survived a whole suite before — it is
  in `zone.cpp`'s own comments twice. The test that kills it drives a pre-loaded
  patch in *x*, where the rest configuration has moved 1.22e-4 m, a third of the
  travel.
- **A linear ramp instead of the smoothstep** is killed by the energy account
  rather than by any displacement: a velocity step at the end of the ramp rings
  the patch, and the account notices before the field does.

### Solver

Explicit central-difference time integration inside Tier 2 (standard for
impact/fracture, no global stiffness matrix, trivially parallel). Tier 1 is
implicit and cheap — measured at 0.35–0.41 core-seconds per simulated second at
1 ms steps against Tier 2's 1155 on the same plating, because a linear model is
unconditionally stable and takes a thousand steps where the explicit one takes
580 000. Lumped mass matrix in both, so the two agree about inertia at the moment
a region is promoted.

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
