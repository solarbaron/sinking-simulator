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
5. **Progressive flooding through the structure.** Requires §3: once elements
   tear, new orifices appear automatically.

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
dampers are switched off. Roll keeps its stand-in, because roll radiation damping
is genuinely small and the mechanism that matters is viscous; deleting it before
Ikeda is wired into this integrator would leave the mode that most needs damping
with almost none. Surge keeps its damper too, since strip theory contributes no
surge radiation at all. Quadratic drag is untouched: it is viscous and separate.

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

**Limits, beyond the ones the radiation module already lists.** Only the
*diagonal* of `A_inf` is used: the integrator inverts a 3×3 inertia and divides by
a per-axis mass, so there is nowhere for the sway–roll coupling `A_24` to go. That
coupling is real and matters for roll. And `A_inf`'s rotational block is
referenced to the body-frame origin on the baseline, while the inertia it is
added to is about the centre of gravity — a frame mismatch that affects roll and
must be resolved when Ikeda damping is wired in and roll is validated properly.
Heave is unaffected (translation-invariant) and pitch is unaffected given strip
theory's zero surge added mass.

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
rudder components, and feeding B44 into the rigid-body roll equation in place of
`Ship::zetaRoll`.

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

## 3. Structure — adaptive tetrahedral FEM

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
100 Hz for free.

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

### Element technology — revised after measurement

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
  through the thickness, no locking, and the timestep is governed by in-plane
  size instead of plate thickness — worth 5–10× on the step alone, on top of the
  element count reduction.
- **Tetrahedra where the geometry really is three-dimensional**: castings, engine
  seats, thick brackets, and the crush zone once plating has folded and shell
  kinematics no longer apply.
- **Promotion from shell to tet** as an element crumples past that point.

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
- **Collision**: two deformable ships, both FEM-active in the contact zone.
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
