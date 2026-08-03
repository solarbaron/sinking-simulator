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
one-off. Evaluation is ~11 ns per component for elevation and ~20 ns for full
kinematics, so a 576-component sea costs 6.3 µs per elevation query and 11.4 µs
per kinematics query. That is affordable for hundreds of query points per tick
and **not** for thousands: a 2000-panel hull wanting kinematics at 100 Hz is
23 ms/tick on one core, against a 10 ms budget. The naive `std::cos` loop is the
whole of it, so the work is a vectorised sincos over the component array, and a
per-tick phase recurrence where the query point is not moving. Worth doing before
the Froude–Krylov integration, not after.

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

The honest way to do wave loads on a large ship is potential flow. Plan:

- **Offline:** solve the boundary-element (panel) problem over the hull for a grid
  of frequencies and headings, producing added mass A(ω), radiation damping B(ω),
  and diffraction/Froude–Krylov exciting forces. Use NEMOH or Capytaine as the
  reference solver; ship the coefficient tables as ship assets.
- **Runtime:** Cummins impulse-response formulation. The radiation force has
  memory — it depends on the history of motion — so the convolution integral is
  approximated by a fitted state-space model (4–8 states per DOF, identified
  offline). This is the standard marine-simulation approach and it is cheap at
  runtime.
- **Nonlinear Froude–Krylov and restoring**: integrated over the *instantaneous*
  wetted surface each tick rather than the mean position. The engine already has
  exactly the routine this needs. This is what captures a ship's behaviour in
  large waves, where linear theory quietly stops being true.

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

## 7. Propulsion, manoeuvring and machinery

- **Propellers**: open-water curves (K_T, K_Q vs J) per propeller, wake fraction
  and thrust deduction from the hull, four-quadrant data so astern and crash-stop
  work, cavitation inception and its thrust breakdown, and ventilation when the
  propeller emerges in a seaway.
- **Manoeuvring**: MMG-style modular model — hull, propeller, rudder terms
  separately identified — rather than a monolithic derivative set, so damage to
  one part degrades one term. Shallow water effects, bank effect, ship-to-ship
  interaction, and current.
- **Steering**: rudder with stall, hydraulic steering gear with real rates and
  failure modes, azimuth thrusters, Voith-Schneider, waterjets, bow and stern
  thrusters with their speed-dependent effectiveness loss.
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
