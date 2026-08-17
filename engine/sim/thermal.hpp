// SPDX-License-Identifier: MIT
//
// Implicit heat conduction on the structural mesh -- the first item of
// `docs/06-roadmap.md` Phase 4 that everything else in that phase needs.
//
// The milestone Phase 4 is aimed at is "an engine room fire that heats a
// bulkhead until it fails under the head of water behind it". Three things stand
// between here and there: a compartment fire that says what the gas does, a
// conduction solve that says what the steel does, and a strength model that says
// what the steel then carries. This file is the middle one **and, since the
// milestone landed, the third as well** -- `carbonSteelReduction`,
// `carbonSteelElongation` and `HeatedMember` below are what a heated member
// carries, and `fire::wallExchange` is what hands the boundary across. There is
// still no combustion here, no species transport and no radiation view factors:
// radiation crosses the boundary as a film coefficient formed on the *fire* side,
// where the exchange `sigma (T_g^4 - T_s^4)` factors exactly into
// `h (T_g - T_s)`, so the conduction operator stays linear and nothing here has
// to know about it.
//
// --- Why implicit, and a correction to the usual argument ----------------------
//
// The usual argument is that explicit conduction, stable only below
// `rho c h^2 / (2 k)`, is hopeless on ship plating because that limit is
// *milliseconds*. **It is not, and the figure is wrong by three orders of
// magnitude.** Steel's thermal diffusivity is
// `alpha = 53.3 / (7850 * 440) = 1.545e-5 m^2/s`, so for the ferry's 12 mm
// plating the limit is `0.012^2 / (2 * 1.545e-5) = 4.66 s`. Seconds, not
// milliseconds. `explicitLimit()` returns exactly that -- to machine precision,
// because for a rectangular element the largest generalised eigenvalue of
// `(K, C_lumped)` is exactly the jump through the thinnest direction, which
// `tests/test_thermal.cpp` asserts as a closed form rather than a tolerance.
//
// Milliseconds would need 0.3 mm elements. **Steel is a poor conductor by the
// standards of this estimate**: it is 1/25 of copper's diffusivity, and a 12 mm
// plate equilibrates through its own thickness in about nine seconds.
//
// So explicit conduction is *viable* on the unrefined structural mesh, and it is
// worth saying so plainly rather than repeating a figure that makes the decision
// look easier than it is. Three things make implicit the right choice anyway,
// and the first two are measured in the test file:
//
//   * **The limit falls as `h^2` and the coupling will refine.** Four elements
//     through the same 12 mm plate -- which is what a through-thickness
//     temperature gradient needs to be resolved, and a gradient through the
//     thickness is what bows a plate -- takes it to 0.29 s. A 1.5 mm surface
//     layer resolving a fire's thermal boundary layer takes it to 0.073 s.
//     Implicit pays nothing for either.
//   * **The step should be set by the fire, not by the mesh.** A compartment fire
//     is a twenty-minute event whose boundary condition is smooth over tens of
//     seconds. Backward Euler runs it in as many steps as the *physics* needs:
//     the measured error at 80x the explicit limit on the semi-infinite problem
//     -- twelve steps over a boundary condition that is a step discontinuity --
//     is 10.0 K in an 880 K jump, 1.1%. Backward Euler is first order and that is
//     what first order costs; it is stable at any step, which is the point.
//   * **The factorisation is reused.** At a fixed step `C/dt + K` is factored
//     once and every later step is two triangular solves -- `n b` rather than
//     `n b^2`. That is what makes a long implicit step cheap enough that trading
//     it against a hundred explicit ones is a win rather than a wash.
//
// The `h^2` in the explicit limit is the whole argument, and it is an argument
// about where this is *going*, not about where the mesh is today.
//
// --- Temperature is in KELVIN --------------------------------------------------
//
// Every temperature crossing this interface is absolute, per CLAUDE.md's "SI
// everywhere". That is not a free choice dressed up as a convention: Phase 4's
// next item after this one is radiation, `sigma (T_gas^4 - T_steel^4)`, which is
// meaningless in Celsius, and a module that took Celsius here and Kelvin there
// would produce a plausible wrong number in exactly the way this repo keeps
// finding. The material curves below are *published* in Celsius and convert on
// the way in, in one place, which is `kCelsius`.
//
// --- What this reuses, and one thing it does not -------------------------------
//
// The mesh, the banded solver and the Cuthill-McKee ordering are
// `solid_shell.hpp`'s and `reduction.hpp`'s. `C/dt + K` is symmetric positive
// definite for any positive step -- C is a mass matrix and K is positive
// semi-definite -- so `solidshell::BandedSpd` is exactly the right factorisation
// and a second solver would be a second place to be wrong.
//
// **`solidshell::RestForms` is not reusable here, and it is worth saying why.**
// It looks like it should be: it caches the Gauss weights and the Jacobians a
// conduction operator wants. The weights it does cache are the same object, and
// `gaussVolumes()` is the check that this file's are right. But its `b` is a
// 6 x 24 *strain*-displacement matrix for a vector field, and every one of the
// assumed-strain cures that make the solid-shell an element rather than a locked
// hex has already been applied to it: row 2 is Betsch-Stein sampled at the
// in-plane corners, rows 4 and 5 are Dvorkin-Bathe sampled at the mid-edges.
// Reading `grad N_a` out of it would be reading a matrix for something it does
// not mean. And it should not be cured anyway -- locking is a property of a
// constrained vector field, and the scalar Laplacian has no shear to lock.
// `Forms` below is therefore its own object: the plain Cartesian shape gradient,
// the shape functions, and the same `det J`.
//
// Body frame and SI units per CLAUDE.md.
#pragma once

#include "plasticity.hpp"    // Material, for the strength model at temperature
#include "scantlings.hpp"    // StructuralMaterial
#include "solid_shell.hpp"   // HexMesh, BandedSpd, kNodes/kGauss/kDof

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace sim::thermal {

using solidshell::kDof;
using solidshell::kGauss;
using solidshell::kNodes;

// 0 degrees Celsius, for the material curves and for a caller writing a fire
// temperature down the way a fire test reports it.
inline constexpr double kCelsius = 273.15;

// --- Carbon steel, EN 1993-1-2:2005 --------------------------------------------
//
// The one published, dimensionally complete thermal description of structural
// steel over a fire's range, and the one every structural fire calculation in
// Europe is done against. Sources, by clause:
//
//   * conductivity, 3.4.1.3: `lambda = 54 - 3.33e-2 theta` W/(m K) over
//     20 <= theta < 800 C, and a flat 27.3 above it to 1200 C.
//   * specific heat, 3.4.1.2: a cubic to 600 C, then two hyperbolae, then a flat
//     650 J/(kg K) from 900 C to 1200 C.
//   * density, 3.4.1.4: 7850 kg/m^3, taken as independent of temperature. It is
//     not quite -- steel expands about 1.4% by volume at 700 C -- but the code
//     says to treat it as constant and mixing a constant density with a
//     temperature-dependent everything else is the kind of half-applied model
//     that reads as a physical effect.
//
// **Both matter over a fire's range and the second matters enormously.**
// Conductivity falls 36%, from 53.3 to 34.0 W/(m K), between 20 C and 600 C:
// steel gets *worse* at spreading heat exactly as it starts to lose strength.
// Specific heat is the larger effect and it is not monotone -- the ferrite to
// austenite transition puts a spike at 735 C where `c` reaches
// **5000 J/(kg K), eleven times its room-temperature value**. That spike is a
// latent heat wearing a heat capacity's clothes, and a solve that misses it
// gallops a plate through 735 C in one step that should have taken eleven.
//
// So it is modelled, and `Problem::temperatureDependent` switches it on. What
// makes it affordable, and what makes the energy account still close to machine
// precision across the spike, is that the capacity is formed from a **secant**
// and not a tangent -- see `Solver` below.
//
// Outside 20..1200 C the standard says nothing and these clamp to the end
// values, which keeps `c` integrable and `k` positive rather than extrapolating a
// cubic into a region where it turns over.
double carbonSteelConductivity(double kelvin);  // W/(m K)
double carbonSteelSpecificHeat(double kelvin);  // J/(kg K)

// Specific enthalpy above 20 C: the exact integral of `carbonSteelSpecificHeat`,
// piece by piece, with the constants of integration chosen to make it continuous.
//
// This exists because **`c` is not what conserves energy once `c` depends on
// temperature** -- `h` is. A step from `T_a` to `T_b` moves `rho (h(T_b) -
// h(T_a))` joules per cubic metre, and that is only `rho c dT` when `c` is
// constant over the interval. Across the 735 C spike it is not: a 30 K step
// through it carries about six times what the endpoint `c` values would suggest.
double carbonSteelEnthalpy(double kelvin);  // J/kg

// --- Strength at temperature, EN 1993-1-2:2005 §3.2 -----------------------------
//
// The other half of the same standard, and the half a fire is actually *for*:
// conduction says what the steel's temperature is, and this says what the steel
// then carries. Three factors, all dimensionless, all 1 at 20 C, tabulated at
// hundred-degree stations from 20 C to 1200 C with **linear interpolation between
// them** -- the interpolation is the standard's, not a smoothing applied here, so
// the kinks are load-bearing and are asserted rather than rounded off.
//
//   * `effectiveYield`, `k_y,theta`: the strength at 2% strain. This is the factor
//     every member-resistance calculation in EN 1993-1-2 uses, and it is the one
//     that makes 600 C the number fire engineers quote -- 0.47 there, so a little
//     over half the strength is gone.
//   * `proportionalLimit`, `k_p,theta`: where the response stops being linear.
//     **It falls far faster than the yield does**, and that difference is the whole
//     shape of the hot stress-strain curve: 0.18 at 600 C against the yield's 0.47,
//     a factor of 2.6. A model that scales one number and calls it "the yield" is
//     making a choice, and `carbonSteelStress` below is what measures its size.
//   * `youngsModulus`, `k_E,theta`: 0.31 at 600 C. **Also below the yield factor**,
//     from 500 C up, and that is not a detail -- an elastic buckling stress is
//     proportional to E and a squash load is proportional to f_y, so a structure
//     that fails by instability loses strength *faster* than `k_y` says. Measured
//     on the reference ferry at 600 C: the midship section keeps 0.356 of its
//     ultimate sagging moment against `k_y` = 0.47, and one 0.7 x 2.4 m panel of
//     12 mm plating keeps 0.323. At **400 C**, where `k_y` is still exactly 1, the
//     section has already lost 21.2% -- all of it through `k_E`.
//
//     (These read 0.376 and 17.6% while the stiffener column check used the
//     transverse frame's profile: a frame is stiff enough that `min(yield,
//     buckling)` took the yield every time, so only the plate panels carried the
//     `k_E` sensitivity. Corrected with the profile, and this header was the site
//     that correction did not reach.)
//
// **What is not here, and what it costs.** Thermal elongation used to be the first
// item on this list and is now `carbonSteelElongation` below, where the estimate
// made here -- yield on expansion alone at 164.6 C -- is reproduced to 164.630 C
// against the element rather than against a hand solve. What that estimate got
// wrong is not the number but the failure mode; see there.
//
// **Creep, §3.2.4, is absent and it does not matter here -- with a number.**
// EN 1993-1-2 §3.2.1 states that the curves above already carry creep implicitly,
// for heating rates between **2 and 50 K/min**. Both fires measured in
// `tests/test_thermal.cpp` sit at or above that band: the ISO 834 standard curve
// takes this ferry's 12 mm plating from 20 C to 618 C in 14 minutes, **43 K/min**,
// which is inside it; a post-flashover compartment held at 900 C does it in five,
// **114 K/min**, which is above it and is the direction in which the implicit
// creep is *conservative*. So for a fire that grows, creep is already in the
// numbers. What is outside them is a fire that stabilises and holds -- a
// twenty-minute soak at temperature under sustained load -- and that is the case
// to revisit, not the growth phase.
//
// Outside 20..1200 C these clamp, matching `carbonSteelConductivity` and for the
// same reason. Below 20 C they clamp to 1: steel is *stronger* cold, and returning
// the room-temperature strength there is the conservative direction as well as the
// only one the standard supports.
struct SteelReduction {
    double effectiveYield = 1.0;     // k_y,theta -- f_y,theta / f_y
    double proportionalLimit = 1.0;  // k_p,theta -- f_p,theta / f_y
    double youngsModulus = 1.0;      // k_E,theta -- E_a,theta / E_a
};

SteelReduction carbonSteelReduction(double kelvin);

// The three separately, for a caller that wants one. Each is exactly the matching
// member of `carbonSteelReduction`, which the tests assert rather than assume.
double carbonSteelYieldFactor(double kelvin);
double carbonSteelProportionalFactor(double kelvin);
double carbonSteelModulusFactor(double kelvin);

// The standard's own stress-strain relation, §3.2.1 Table 3.1 and Figure 3.1, at
// engineering strain `strain` and temperature `kelvin`. Four branches:
//
//   * linear at `k_E E` up to `f_p = k_p f_y`;
//   * an **ellipse** from there to `f_y,theta = k_y f_y` at exactly 2% strain,
//     constructed to be tangent to the line at one end and horizontal at the
//     other -- both of those are exact and both are asserted;
//   * a plateau at `f_y,theta` to 15%;
//   * a straight fall to zero at 20%.
//
// It is here as the **reference the reduced material is measured against**, not as
// something the solver evaluates. `plasticity::Material` is a J2 flow curve with
// one yield strength, and `atTemperature` below scales it by `k_y`; that model is
// exact at and beyond 2% strain and is *over-strong* below it, because the real
// curve has already left the straight line at `k_p f_y`.
//
// **The size of that, measured, and it is not where it looks like it should be.**
// The gap is bounded by `(k_y - k_p) f_y` by construction, and the temperature that
// maximises it is not the hottest one -- it is **400 C**, where `k_y` is still
// exactly 1 and `k_p` has already fallen to 0.42. Measured there: 0.388 f_y, at
// 0.24% strain. At 600 C it is 0.189 f_y and at 700 C 0.095 f_y, because by then
// both factors are small and their difference is smaller. At 20 C it is exactly
// zero -- `k_p == k_y == 1`, the ellipse degenerates, and the standard's own curve
// is elastic-perfectly-plastic. Past 2% strain, which is where a collapse, a
// squash load or a plastic hinge lives, the two agree identically at every
// temperature.
//
// Negative strain returns the negation of the positive branch: the standard writes
// the curve in tension and says the compressive response is the same, which for a
// hot member is the interesting direction.
double carbonSteelStress(double strain, double kelvin, double yieldStrength = 355.0e6,
                         double youngsModulus = 210.0e9);

// --- Thermal elongation, EN 1993-1-2:2005 §3.4.1.1 --------------------------------
//
// **The largest single thing a fire does to structure, and it is not a weakening.**
// The relative elongation from 20 C, `dl/l`, in three branches:
//
//     20 <= theta <  750   1.2e-5 theta + 0.4e-8 theta^2 - 2.416e-4
//     750 <= theta <= 860   1.1e-2                       (the phase change)
//     860 <  theta <= 1200  2e-5 theta - 6.2e-3
//
// The plateau is the ferrite-to-austenite transition: face-centred cubic packs
// more densely than body-centred cubic, and over that 110 K the lattice
// contraction cancels the thermal expansion exactly. It is the same transition the
// 735 C spike in `carbonSteelSpecificHeat` is -- one physical event appearing in
// two of the standard's clauses, and a model that carried one without the other
// would be describing half a phase change.
//
// **The standard's own curve is discontinuous at 750 C, by 8.4e-6, and it is
// asserted rather than smoothed.** The polynomial reaches 1.100840e-2 at 750 and
// the plateau is 1.1e-2 flat, so the elongation steps *down* on entering the phase
// change. At 860 the two branches agree to the last bit -- `2e-5 * 860 - 6.2e-3` is
// exactly 1.1e-2 in binary as well as in decimal. So the curve is continuous at one
// end of the plateau and not at the other, which is a fact about EN 1993-1-2 and
// not about this file. It is 0.076% of the elongation there and 1.73 MPa of
// restrained stress; interpolating it away would be inventing a curve the standard
// does not publish, and `tests/test_thermal.cpp` asserts the jump's exact size so
// that a later smoothing has to argue with a test.
//
// **`carbonSteelElongation(860 C + kCelsius)` still does not return 1.1e-2**, and
// the reason is the one this file already records at `carbonSteelReduction`: no
// double satisfies `k - 273.15 == 860`. The round trip lands 1.1e-13 K high, over
// the branch boundary, so the linear expression is evaluated a hair past 860 and
// the answer is 2.3e-15 out. That is 1e-13 K of *representation* against the 21 K
// it would take to matter, and it is asserted rather than hidden behind a
// tolerance, because a caller who compares against 1.1e-2 with `==` should find
// the property written down.
//
// --- What it is: an eigenstrain -------------------------------------------------
//
// A thermal strain is not a stress and it is not a load. It is strain the material
// carries for free, subtracted from the total before the constitutive law is
// shown anything:
//
//     sigma = C (eps_total - eps_th)
//
// `solid_shell.hpp` says where that subtraction goes and why the equivalent nodal
// force is a consequence of it rather than an alternative to it. What belongs here
// is the consequence: **heating a body generates no stress. Preventing it from
// expanding does.**
//
//   * A free bar, uniformly heated, carries **exactly zero** stress at any
//     temperature. Asserted at zero and not at a tolerance, because that is the
//     assertion that says the eigenstrain was subtracted and not added.
//   * A bar restrained against expanding carries `sigma = -k_E(T) E eps_th(T)`,
//     compression, and reaches `k_y(T) f_y` at **164.630 C** for this repo's
//     `ah36Steel()` -- a 144.6 K rise, `k_E` = 0.9354, `k_y` exactly 1. The whole
//     of the effect is expansion and none of it is strength loss, which is the
//     figure `restrainedYieldTemperature` returns and the element reproduces.
//   * **164.6 C is right for E = 206 GPa and only for E = 206 GPa.** At the
//     210 GPa that `carbonSteelStress` takes as its default it is 161.546 C. These
//     take the modulus from the material rather than from a constant for that
//     reason.
//
// --- What "restrained" means on a ship, which is the part that decides ------------
//
// A fully restrained bar is a laboratory object. A bulkhead stiffener is held by
// structure that is itself hot and itself expanding, and **a uniformly heated
// region expands freely and generates no stress at all** -- so the question is
// never "how hot" on its own.
//
// It is also not answered by a length: an axially heated *chain* of members with
// free ends is statically determinate and carries **zero** stress however the
// temperature varies along it. The bar just gets longer. What generates stress is
// a **transverse** gradient -- a hot strake beside a cold one, sharing a seam, the
// hot one wanting to grow along the ship and the cold one holding it -- which is
// exactly the shape of a fire in one compartment of many.
//
// For two parallel strips tied at both ends and free overall, hot fraction `f` of
// the area, the hot strip carries
//
//     sigma = -eps_th / (1/E_hot + (f / (1 - f)) / E_cold)
//
// and the fully restrained answer is its `f -> 0` limit. **The limit is approached
// fast, and that is why this is not a laboratory result**: with a tenth of the
// section hot the hot strip yields at 181.5 C against the fully restrained 164.6,
// and with a fiftieth at 167.7. Half the section hot is the case that is genuinely
// different -- 313.6 C -- and half a ship's midship section is not what one
// compartment burns.
//
// --- And buckling arrives first ---------------------------------------------------
//
// **A restrained member that cannot expand goes into compression, and compression
// is not a yield problem.** Against `buckling.hpp`'s own checks on this ferry's own
// scantlings, with the restraint stress as the applied compression:
//
//   * the 8 mm plating of the vehicle deck head, in its 0.70 x 2.40 m bay, buckles
//     at **59.0 C -- a 39 K rise**. A hot shower would do it, if the deck were
//     truly restrained.
//   * the 12 mm side plating, same bay, at **103.1 C**;
//   * the 14.5 mm bottom plating at 121.1 C;
//   * the side longitudinal with its attached plating, as an Euler column over the
//     2.40 m frame bay, at 149.8 C -- only 9% earlier than yield, because it is a
//     stocky column and Johnson-Ostenfeld caps it near `f_y`.
//
// So the ordering is plating, then stiffener, then yield, and the plating goes at
// **a third of the temperature rise** the yield figure quotes. 164.6 C is the right
// answer to the wrong question.
//
// **Two closed forms make that more than a table.** In the elastic branch the
// modulus reduction factor **cancels exactly**: the restraint stress is
// `k_E E eps_th` and the elastic buckling stress is `k_E E` times a pure geometry,
// so the temperature at which they meet is a function of geometry and the
// elongation curve **alone** -- independent of E, of `k_E`, and of the steel grade.
// For an Euler column that reads `eps_th(T) = pi^2 / lambda^2` with `lambda` the
// slenderness, and equating it to the 164.630 C yield strain gives a **critical
// slenderness of 73.19**: above it a fully restrained column buckles first, below
// it it would yield first were it not for the plasticity cap. The ferry's
// longitudinals sit at lambda = 34 to 45, which is why they pre-empt yield only
// narrowly and the plating pre-empts it by a mile.
//
// Limits, inherited from `buckling.hpp` and worth naming because they all point the
// same way: these are ideal flat panels under *uniaxial* compression. Real
// restrained heating is biaxial, real panels carry initial distortion and welding
// residual stress, and every one of those makes the real temperature **lower**.

// Relative elongation `dl/l` from 20 C, dimensionless. Clamped outside 20..1200 C,
// matching `carbonSteelConductivity` and for the same reason.
double carbonSteelElongation(double kelvin);

// The eigenstrain between two temperatures: `elongation(kelvin) -
// elongation(reference)`. This -- not `carbonSteelElongation` -- is what an element
// subtracts, because a structure is stress-free at the temperature it was built
// at and that need not be 20 C.
double thermalStrain(double kelvin, double referenceKelvin = kCelsius + 20.0);

// The same as the six Voigt components `solid_shell.hpp` wants: `{e, e, e, 0, 0, 0}`.
// A thermal expansion of an isotropic material is pure dilatation and carries no
// shear at all, which is why a J2 yield surface never sees it.
void thermalEigenstrain(double kelvin, double referenceKelvin, double out[6]);

// Stress in a member held to exactly its reference length: `-k_E(T) E eps_th`,
// **negative in compression**, the same sign convention as `elementStress`.
//
// `material` is the material at 20 C -- the reduction factor is applied here, so
// handing this the output of `atTemperature` would apply `k_E` twice.
double restrainedStress(const StructuralMaterial& material, double kelvin,
                        double referenceKelvin = kCelsius + 20.0);

// The temperature at which `|restrainedStress|` first reaches `k_y(T) f_y`: the
// 164.630 C anchor for `ah36Steel()` from 20 C. Kelvin.
//
// Zero -- not a temperature -- when the reference is already at or above 1200 C,
// because there is nothing above it to search.
//
// **1200 C is the answer for a degenerate material, and that is not "never".**
// EN 1993-1-2 takes both `k_y` and `k_E` to *exactly* zero there, so a member at
// the top of the range has no strength left and is at its limit whatever stress it
// is carrying -- a material with no stiffness at all comes back 1200 C rather than
// zero. Both of these searches therefore always find a crossing, which mutation
// testing established rather than the other way round, and
// `tests/test_thermal.cpp` asserts the two factors that make it so.
double restrainedYieldTemperature(const StructuralMaterial& material,
                                  double referenceKelvin = kCelsius + 20.0);

// The temperature at which `|restrainedStress|` first reaches the member's own
// buckling capacity, `johnsonOstenfeld(k_E(T) elasticStress, k_y(T) f_y)`.
//
// `elasticStress` is the member's **elastic** buckling stress at 20 C -- exactly
// `BucklingCheck::elasticStress` out of `plateBuckling` or `columnBuckling` on the
// unheated material -- so one routine serves a plate and a column and there is no
// second place for the reduction to be applied. Zero if it never does.
double restrainedBucklingTemperature(double elasticStress, const StructuralMaterial& material,
                                     double referenceKelvin = kCelsius + 20.0);

// --- Restraint that is not total, and a temperature that is not uniform ----------
//
// Two closed forms that turn the paragraphs above into functions, because the
// milestone this file exists for -- a bulkhead heated on one side and loaded on the
// other -- has neither a fully restrained member nor a uniform one.
//
// **The two-strip stress.** `hotFraction` is `f` in the header note above: two
// parallel strips tied at both ends and free overall, the hot one carrying
//
//     sigma = -eps_th / (1/E_hot + (f / (1 - f)) / E_cold)
//
// written here as `-k_E E eps_th / (1 + (f/(1-f)) k_E)` so that `f = 0` returns
// `restrainedStress` **bit for bit** rather than to a rounding -- the two
// expressions are algebraically the same and `1/(1/E)` is not `E` in binary. `f = 1`
// is a body with nothing cold left to hold it and returns exactly zero.
//
// The ratio to the fully restrained answer is `1 / (1 + (f/(1-f)) k_E(T))`, which is
// the *only* place the modulus survives: it is a stiffness ratio, and the hot strip
// is the softer of the two.
double twoStripStress(const StructuralMaterial& material, double kelvin, double hotFraction,
                      double referenceKelvin = kCelsius + 20.0);

// **The equivalent uniform temperature.** A member's restraint force is set by the
// elongation it is prevented from making, which is the *integral* of `eps_th` along
// it and not `eps_th` at its mean temperature -- the elongation curve is quadratic,
// so those differ. `temperatureForElongation` inverts `carbonSteelElongation`, so a
// caller averages the elongations of its own temperature field and asks what one
// temperature would have produced it. Everything downstream is then uniform, which
// is what keeps the `k_E` cancellation below exact.
//
// The **lowest** such temperature, and the plateau makes that a real distinction
// rather than a formality. EN 1993-1-2's curve is flat at 1.1e-2 over 750-860 C, so
// every temperature in that band inverts to one number -- and that number is
// **749.533 C, not 750**, because the curve *steps down* by 8.4e-6 on entering the
// plateau (see the note above) and the polynomial branch therefore reaches 1.1e-2 a
// little before the branch boundary. 860 C is the exception and for the reason
// recorded above: `carbonSteelElongation(860 C)` lands 2.3e-15 *over* the plateau
// value, so it inverts to itself. Below 750 C, which is where a member that still
// carries anything lives, the curve is strictly increasing and the round trip closes
// to a rounding of the temperature itself.
//
// Clamped to 20 C below zero elongation and to 1200 C above the curve's top, on the
// same terms as `carbonSteelElongation` itself.
double temperatureForElongation(double elongation);

// --- A heated member that is also carrying a load --------------------------------
//
// The Phase 4 milestone in one function: "an engine room fire that heats a bulkhead
// until it fails **under the head of water behind it**". The header above establishes
// what heat alone does to a restrained member -- compression, and buckling before
// yield. This is what happens when that member is *also* being bent, which is the
// only case in which the milestone's own sentence can be true.
//
// --- Why it is not a sum -----------------------------------------------------------
//
// A member under axial compression `N` and a lateral load deflects, and the axial
// load then acts through that deflection. For a pin-ended beam-column under a uniform
// lateral load the mid-span moment is **exactly**
//
//     M = M_0 * 2 (sec u - 1) / u^2,     u = (pi/2) sqrt(N / N_E)
//
// (Timoshenko, *Theory of Elastic Stability* §1.11), with `M_0` the first-order
// moment and `N_E` the Euler load. It is 1 at `u = 0` -- exactly, and its series
// there is `1 + 5u^2/12 + ...`, because `(sec u - 1)/u^2` is 0/0 at the origin
// and evaluating it as written loses every digit to cancellation nearby.
//
// **The implementation does not use that series**, and this said it did.
// `beamColumnMagnifier` rewrites the whole thing as `sinc(u/2)^2 / cos u`, which
// is the same number exactly with no subtraction in it, and `thermal.cpp` says so
// in its own words: "No series expansion and therefore no truncation term to
// size." Two comments on one function, describing two implementations, one of
// which had been gone long enough that nobody could say when. Nothing tests a
// comment; this one was wrong in the direction that invites a reader to go
// looking for a truncation error that does not exist.
// It diverges at `N = N_E`. **That divergence is the whole mechanism**: the fire
// supplies `N`, the head of water supplies `M_0`, and neither on its own is what
// fells the member.
//
// --- `k_E` cancels here too --------------------------------------------------------
//
// `u` depends on `N/N_E`, and both are proportional to the reduced modulus: the
// restraint stress is `k_E E eps_th` and the Euler stress is `k_E` times the cold
// one. So **the magnification is a function of the elongation curve and the member's
// geometry alone** -- independent of `k_E`, of `E` and of the steel grade -- exactly
// as the restrained buckling temperature is. Asserted by changing all three.
//
// --- The check ----------------------------------------------------------------------
//
// EN 1993-1-2 §4.2.3.5's shape, with this repo's own pieces rather than the
// standard's tabulated interaction factors -- a table reproduced from memory is a
// plausible wrong number, and the exact magnifier above is a better object than a
// code fit anyway:
//
//     utilisation = |sigma_N| / sigma_column(T)  +  Psi(u) |sigma_M| / (k_y(T) f_y)
//
// with `sigma_column(T) = johnsonOstenfeld(k_E(T) sigma_E, k_y(T) f_y)`, which is the
// same capacity `restrainedBucklingTemperature` uses and reduces to it exactly when
// there is no lateral load at all.
//
// **What it is not.** First yield of the extreme fibre, not a plastic mechanism, so
// it is conservative by the shape factor -- 1.5 for a rectangle, less for a stiffener
// with its plating. Pin-ended, so a member with real end fixity fails later. Elastic
// magnification, so it says nothing about what happens after the first hinge. And no
// redistribution: a member that has gone hands its load to its neighbours, and
// nothing here notices.
struct HeatedMember {
    // Section modulus to the fibre the lateral load puts in compression, m^3.
    double modulus = 0;
    // First-order moment from the lateral load, N m. The caller integrates its own
    // pressure distribution: a bulkhead's is hydrostatic over part of its span and
    // zero above, which is not `q L^2 / 8`.
    double lateralMoment = 0;
    // Elastic column buckling stress at 20 C, Pa -- `BucklingCheck::elasticStress`
    // out of `columnBuckling` on the unheated material, exactly as
    // `restrainedBucklingTemperature` takes it.
    //
    // **Zero is a member with no buckling capacity, not one that cannot buckle**, so
    // a default-constructed `HeatedMember` fails the instant it is heated above its
    // reference temperature. That is deliberate and it is the conservative
    // direction: the alternative reading -- an unfilled field meaning "infinitely
    // stiff" -- would let a caller who forgot it get a member that never buckles at
    // all, which is the failure this repo keeps finding under the name "fails open".
    double eulerStress = 0;
    // The share of the fully restrained stress the member actually carries, 1 for a
    // laboratory bar. **The model cannot derive this**: it is set by the stiffness of
    // the structure at the member's ends, which is a bulkhead deck and a tank top and
    // is in neither the fire model nor the conduction one. `twoStripStress` gives it
    // for the one case that *is* derivable -- a hot strip held by a cold one beside
    // it, in the same plane.
    double restraint = 1.0;
};

// Which term reached one first.
enum class MemberLimit { None, Column, Interaction };

struct MemberState {
    double axialStress = 0;       // Pa, magnitude, compression
    double eulerStress = 0;       // Pa, at temperature: k_E times the cold value
    double axialOverEuler = 0;    // the magnifier's argument
    double magnifier = 1;         // Psi(u)
    double bendingStress = 0;     // Pa, magnified
    double columnCapacity = 0;    // Pa, Johnson-Ostenfeld at temperature
    double yieldCapacity = 0;     // Pa, k_y f_y
    double utilisation = 0;       // the interaction; >= 1 has failed
    // The same sum with the magnifier forced to 1. Published because the difference
    // between the two **is** the coupling: a purely additive check is what two
    // subsystems that do not know about each other would produce, and the gap is the
    // measure of what the third one buys.
    double additiveUtilisation = 0;
    MemberLimit limit = MemberLimit::None;
};

// `2 (sec u - 1) / u^2` at `u = (pi/2) sqrt(x)`, `x = N / N_E`. Exactly 1 at x = 0,
// infinite at x >= 1, and refused (returning 1) for x < 0 -- a member in tension is
// stiffened by its axial load and this formula is not the one for that.
double beamColumnMagnifier(double axialOverEuler);

MemberState memberState(const HeatedMember& member, const StructuralMaterial& material,
                        double kelvin, double referenceKelvin = kCelsius + 20.0);

// The temperature at which `memberState(...).utilisation` first reaches 1, by the
// same scan-and-bisect `restrainedYieldTemperature` uses and over the same range.
// Zero when it never does.
//
// With `lateralMoment` zero this is `restrainedBucklingTemperature` on the same
// member, **to the bit**: the identity that says the interaction was not bolted on
// beside the axial check but contains it. It is not obvious that it should be --
// the two searches bisect the same crossing on differently rounded margins,
// `sigma/sigma_c - 1` against `sigma - sigma_c`, which are free to disagree about
// the sign of a probe that lands between them -- so it is measured over a sweep of
// fifty elastic stresses from 5 MPa to 2 GPa rather than assumed, and every one of
// the fifty agrees exactly.
double memberFailureTemperature(const HeatedMember& member, const StructuralMaterial& material,
                                double referenceKelvin = kCelsius + 20.0);

// --- Materials at temperature ----------------------------------------------------
//
// The coupling itself, and it is deliberately two free functions returning
// *values* rather than a temperature field stored inside a material. A material
// with a temperature in it would have to be kept in step with a thermal solve that
// runs on its own clock, and the repo already records what a cache that quietly
// goes stale costs. A reduced material is a value: build it, use it, throw it away.
//
// **This is also the whole answer to "where does temperature enter".** Every
// element-level entry point in `solid_shell.hpp` already takes one material *per
// call* -- `elementStiffness`, `elementStress`, `elementPlasticUpdate`,
// `criticalTimestep` -- so a per-element temperature needs no interface change at
// all, only a caller that builds a reduced material per element. The same is true
// one tier down: `buckling.hpp`, `collapse.hpp` and `indentation.hpp` all look
// their material up through `SectionElement::material`, an index into
// `StructuralMesh::materials`, so a mesh with a temperature *field* is a mesh with
// one material entry per distinct temperature and nothing else changes.
//
// What that does *not* reach without a signature change is per-Gauss-point
// temperature: `elementPlasticUpdate` loops the eight points against one material,
// and `elementStiffness` forms one elasticity matrix outside the loop.
// `gaussTemperature` below exists so the size of that gap is a measurement rather
// than an argument, and the measurement says **per element is enough**:
//
//   * **Cost.** `atTemperature` on a `plasticity::Material` is 11.5 ns against
//     `elementPlasticUpdate`'s 5.65 us -- **0.20%** of one element update, so
//     rebuilding the material every element every step is free. Per Gauss point
//     would be eight times that, 1.6%, and still cheap; cost is not what decides
//     this.
//   * **What the gradient actually is.** A 12 mm plate with a post-flashover
//     compartment on one face -- gas at 900 C through an effective film of
//     200 W/(m^2 K), which is EN's 25 for convection plus about 175 for the
//     radiation this file does not carry -- and adiabatic on the other, solved by
//     `Solver` with `temperatureDependent`: the spread **across the whole plate**
//     peaks at **19.1 K**, eight seconds in, and is under 1 K by half an hour. The
//     corresponding spread in `k_y` never exceeds **0.039**. Steel plating is not
//     thermally thick -- the Biot number is `h t / k(600 C) = 0.070` -- and a body
//     at `Bi << 1` has no through-thickness gradient to resolve. A ramped fire is
//     gentler still: the ISO 834 curve on the same plate peaks at 7.2 K. Both
//     figures bound any single element's spread whatever `nz` is, because both are
//     the spread of the whole plate.
//   * **What the error would be if the gradient were large.** The quadrature error
//     from reducing at the mean instead of averaging the reductions is at most
//     `|delta slope| * dT / 8` at a kink of Table 3.1, and the sharpest kink is at
//     400 C. Measured: 0.0064 in `k_y` for a 20 K spread across one element,
//     0.032 for 100 K.
//
// So the through-thickness case is settled by the physics rather than by the
// budget. The case per-Gauss-point *would* buy something is an **in-plane** one --
// the edge of a fire crossing a 2.4 m frame bay -- and there the honest answer is
// that a boundary condition the mesh does not resolve is a meshing problem, which
// is the same statement this file already makes about refinement.

// `youngsModulus` scaled by `k_E,theta` and `yieldStrength` by `k_y,theta`. Density
// is left alone, per §3.4.1.4 and for the same reason `Problem` leaves it alone.
// The thermal properties on this struct are the 20 C ones by construction and are
// **not** re-evaluated here -- a caller who wants `k(T)` and `c(T)` has
// `Problem::temperatureDependent`, and quietly moving them would give a solve two
// different conductivities depending on which door it came in by.
StructuralMaterial atTemperature(const StructuralMaterial& material, double kelvin);

// The same for the flow model: `youngsModulus` by `k_E,theta` and **the whole flow
// curve** by `k_y,theta` -- yield strength, hardening modulus, Swift strength
// coefficient and kinematic modulus alike.
//
// Scaling the whole curve rather than only its intercept is what keeps the model
// consistent, and it has an exact consequence: Considere's necking strain is
// `d sigma_y / d eps_p = sigma_y`, and multiplying both sides by the same positive
// number does not move the root, so `plasticity::uniformElongation` is **invariant
// under temperature** -- bit-exactly for a Swift curve, whose necking strain
// `n - eps_0` is built from two fields this does not touch, and to one unit in the
// last place for a linear one, where `(k sigma_y0)/(k H)` is not always
// `sigma_y0/H`. `Failure` is therefore untouched too, which is the
// honest answer rather than a convenient one: EN 1993-1-2 tabulates no ductility,
// hot steel is *more* ductile rather than less, and inventing a temperature
// dependence for the fracture strain would be a plausible wrong number in the one
// place the tearing criterion depends on it.
plasticity::Material atTemperature(const plasticity::Material& material, double kelvin);

// --- One element ---------------------------------------------------------------
//
// The conduction analogue of `solidshell::RestForms`: everything the element
// derives from its geometry, which nothing in a conduction solve ever moves.
// Same eight-node hexahedron, same 2x2x2 rule, same node ordering contract --
// though unlike the mechanical element the ordering is *not* load-bearing here,
// because a scalar Laplacian has no thickness direction to get wrong. A rotated
// element gives a bit-identical answer up to the summation order, and
// `tests/test_thermal.cpp` asserts it.
struct Forms {
    double gradient[kGauss][kNodes][3];  // dN_a/dx_i, Cartesian
    double shape[kGauss][kNodes];        // N_a
    double weight[kGauss];               // 2x2x2 Gauss weight (one) times det J
    bool ok = false;                     // false on an inverted or degenerate element
};

// Fill `out` from an element's node positions. False -- and `out.ok` false -- when
// `det J` is not positive at one of the eight integration points, which is the
// same condition `solidshell::computeRestForms` refuses on.
//
// A **collapsed** hexahedron -- the wedge a degenerate plate panel extrudes to,
// 166 of them on the reference ferry -- is accepted, for the reason
// `solid_shell.hpp` gives at `ElementShape`: its `det J` is exactly zero at the
// coincident corners and only there, and positive at every point the element is
// integrated at.
bool computeForms(const double nodes[kDof], Forms& out);

// --- From a nodal temperature field to a material ---------------------------------
//
// `Solver::temperature()` is per node; a material is per element or per Gauss
// point. These two are the whole of the bridge.

// A nodal field interpolated to the eight Gauss points. Geometry-independent --
// the trilinear shape functions at +-1/sqrt(3) are the same numbers for every
// element -- so it takes no `Forms`, and it uses **the same eight points in the
// same order** as `solidshell`'s integration rule. That is a contract between two
// files rather than a coincidence, so `tests/test_thermal.cpp` asserts it by
// interpolating a linear field and comparing against `solidshell::gaussVolumes`'
// own weighting rather than by reading either file's constants.
void gaussTemperature(const double nodal[kNodes], double out[kGauss]);

// The volume-weighted mean of a nodal field over the element: `integral T dV / V`,
// by the same 2x2x2 rule everything else here uses. This is the temperature a
// per-element material is built at.
//
// It is a *volume* average and not the mean of the eight nodes, and on a
// distorted element those differ for exactly the reason `elementMass` is row-sum
// lumped rather than volume/8.
double elementTemperature(const Forms& forms, const double nodal[kNodes]);

// Per-element mean temperature over a whole mesh, from a nodal field. Empty, and
// `false`, if the field is not one value per node or an element will not form.
bool elementTemperatures(const solidshell::HexMesh& mesh, const std::vector<double>& nodal,
                         std::vector<double>& out);

// Conductance, 8x8 row-major: `integral grad N_a . k grad N_b dV`, with `k` given
// per Gauss point so a temperature-dependent conductivity costs nothing extra.
// Symmetric, positive semi-definite, and **exactly singular on the constant
// field**: `K 1 = 0`, because `sum_a grad N_a = 0` is a partition of unity. That
// identity is not decoration -- it is what makes the energy account below close,
// and `tests/test_thermal.cpp` asserts it directly.
void conductance(const Forms& forms, const double conductivity[kGauss],
                 double out[kNodes * kNodes]);

// Capacity, 8x8 row-major: `integral rho c N_a N_b dV`, with the *volumetric*
// heat capacity `rho c` in J/(m^3 K) given per Gauss point.
//
// `lumped` replaces it by the diagonal of its own **row sums**, which is the same
// lumping `solidshell::elementMass` uses and for the same reason: volume/8 is
// exact only on a parallelepiped. Row-sum lumping leaves `1^T C dT` **exactly
// unchanged** for any `dT`, by the symmetry of C -- so the two capacities differ
// in where the heat sits and never in how much there is, and the energy account
// is identical under both. That is asserted rather than argued.
//
// Consistent is the default. Lumped exists because a consistent capacity gives
// backward Euler a physically impossible *undershoot* ahead of a steep front
// unless `dt >= rho c h^2 / (6 k)`, and a fire's surface is exactly a steep
// front. **Measured**, on 5 mm elements whose criterion is 0.333 s: at a tenth
// of it the interior dips **14.6 K below its own initial temperature**, and at
// three times it the dip is 7e-13 K. Lumped does not dip at either step -- it is
// unconditionally monotone -- and, by the row-sum identity above, moves exactly
// the same energy. So the choice costs nothing in the account and only decides
// whether a caller who has to take short steps sees a wiggle.
void capacity(const Forms& forms, const double volumetricCapacity[kGauss], bool lumped,
              double out[kNodes * kNodes]);

// Largest stable *explicit* step for one element: `2 / lambda_max`, with
// `lambda_max` the largest eigenvalue of `K x = lambda C x` on the lumped
// capacity, from `reduction::generalisedEigen`.
//
// Exposed because the case for solving this implicitly at all is a number, and a
// number should be measured on the element rather than estimated from a nominal
// thickness -- the same argument `solidshell::criticalTimestep` makes, which got
// the mechanical estimate wrong for the same reason. Here it corrects a figure by
// a factor of a thousand: see the header note above.
//
// For a rectangular element this comes out at exactly `h_min^2 / (2 alpha)` --
// the thin-direction jump is the highest mode, and the full 3D checkerboard is
// three times *lower*. The test file asserts that identity, which is also what
// caught a power iteration whose start vector was orthogonal to the mode it was
// looking for.
double explicitLimit(const double nodes[kDof], const StructuralMaterial& material);
double explicitLimit(const solidshell::HexMesh& mesh, const StructuralMaterial& material);

// --- The mesh boundary ---------------------------------------------------------

// A face of the mesh carried by exactly one element. Faces are numbered in
// `solid_shell.hpp`'s node ordering: 0 is the zeta = -1 face, 1 the zeta = +1
// face, then the four sides.
//
// This is here rather than in a test because a boundary condition has to name a
// surface, and naming it by element-and-face-index is how a caller ends up
// applying a fire to the inside of a bulkhead. A caller filters these by
// `centroid` and `normal` -- geometry it can check -- and hands the survivors
// back.
//
// `normal` is outward: the winding of the face table is not consistently
// oriented, so it is fixed here by comparing the face centroid against the
// element centroid rather than by trusting the table.
struct BoundaryFace {
    std::uint32_t element = 0;
    std::uint8_t  face = 0;
    std::uint32_t node[4]{};
    Vec3   centroid{};
    Vec3   normal{};   // unit, outward
    double area = 0;   // m^2, by the same 2x2 rule the surface integral uses
};

// Every exterior face, in element then face order. A face whose area is zero --
// the degenerate side of a collapsed hexahedron -- is dropped, because it carries
// no heat and its normal is not defined.
std::vector<BoundaryFace> boundaryFaces(const solidshell::HexMesh& mesh);

// A boundary condition of the second and third kinds over a set of faces:
//
//     q_in = flux + coefficient * (ambient - T)      W/m^2, positive INTO the solid
//
// The two are one struct because they are one surface integral and separating
// them would be two places to get `dA` wrong. A pure flux sets `coefficient` to
// zero; a pure convective surface sets `flux` to zero. Radiation is *not* here --
// it is a separate roadmap item and it is nonlinear in a way conduction is not.
//
// Sign: `flux` positive heats the solid, whichever way the face points. That is
// deliberate. The alternative -- resolving the caller's flux against the outward
// normal -- would make the answer depend on a winding the caller cannot see, and
// `normal` is published above for a caller who does want to resolve one.
struct Film {
    std::vector<BoundaryFace> face;
    double flux = 0.0;         // W/m^2 into the solid
    double coefficient = 0.0;  // W/(m^2 K)
    double ambient = 293.15;   // K
};

// --- The problem ---------------------------------------------------------------

struct Problem {
    // Borrowed, not owned, and it must outlive the `Solver` built from it. The
    // solve never moves a node, so a mesh that deforms elsewhere invalidates the
    // factorisation rather than being tracked.
    const solidshell::HexMesh* mesh = nullptr;
    StructuralMaterial material{};

    // Take `k` and `c` from the EN 1993-1-2 curves at the local temperature
    // instead of from `material`. `material.density` is used either way.
    bool temperatureDependent = false;
    bool lumpedCapacity = false;

    // Per node. `prescribed` non-zero holds that node at `prescribedValue`, in
    // kelvin, exactly -- moved to the right-hand side the way
    // `solidshell::solveStatic` moves a prescribed displacement, not by a penalty.
    // Empty means nothing is prescribed, which is legal for a transient and
    // singular for a steady solve unless a film holds it.
    std::vector<std::uint8_t> prescribed;
    std::vector<double>       prescribedValue;

    // W/m^3, one per element. Empty means none. This is where a coupled fire
    // would *not* put its heat -- a fire heats a surface -- but it is what a
    // manufactured-solution test needs and it is one line.
    std::vector<double> volumetricSource;

    std::vector<Film> film;
};

// --- The energy account ---------------------------------------------------------
//
// The analogue of `zone::SolveResult`'s work-against-strain-energy balance, and
// the check most likely to catch a wrong quadrature weight. Every figure is in
// joules and accumulated from the start of the run.
//
// With constant properties it closes to **machine precision** rather than to the
// integrator's order -- measured at 5e-15 of the enthalpy moved over sixty steps
// of three different sizes -- and that is a property of the formulation rather
// than of the step size. Summing the discrete system with `1^T`:
//
//     1^T C (T1 - T0)/dt  +  1^T K T1  =  1^T f + 1^T r
//
// and `1^T K = 0` exactly, so the enthalpy rate is the sum of the applied and
// reaction heat with nothing left over -- at any `dt`, converged or not. The one
// thing that can break it is `C` being wrong, which is exactly what the check is
// for.
//
// **With temperature-dependent properties it closes to the Picard tolerance
// instead**, and the difference is worth being exact about rather than rounding
// up to "machine precision". The secant capacity below makes `1^T C dT`
// identically the enthalpy change *of the state the system was solved for*; an
// iteration stopped at 1e-8 K has not quite reached that state, and the leftover
// is the whole of the residual. Measured: 1.4e-12 of the heat supplied at a 1e-6 K
// tolerance and 1.6e-15 at 1e-11 K, which is the scaling a wrong quadrature would
// not have.
struct Account {
    double enthalpy = 0;        // J, above 20 C, in the mesh now
    double enthalpyChange = 0;  // J, since the solver was prepared
    double prescribedHeat = 0;  // J in through nodes held at a temperature
    double filmHeat = 0;        // J in through flux and convective surfaces
    double sourceHeat = 0;      // J in from volumetric sources

    // enthalpyChange - (prescribedHeat + filmHeat + sourceHeat). Joules, and it
    // should be at the rounding of the largest term.
    double residual() const {
        return enthalpyChange - prescribedHeat - filmHeat - sourceHeat;
    }

    // W. The largest row of `A T1 - b` over the *free* nodes after the last solve
    // -- which is zero if and only if the band actually held every term that was
    // scattered into it.
    //
    // It is reported rather than assumed because `BandedSpd::add` **silently
    // drops** an entry outside its band, so a bandwidth computed one element short
    // reads as a slightly soft answer and not as an error. That is the shape of
    // the defect `CLAUDE.md` records against `reduction`'s node ordering, arriving
    // by a different door.
    double equilibriumResidual = 0;
};

// --- The solver ------------------------------------------------------------------
//
//     (C/dt + K) T1 = C/dt T0 + f
//
// Backward Euler: unconditionally stable, unconditionally monotone with a lumped
// capacity, and first order in time. First order is a real cost and it is
// deliberate -- Crank-Nicolson is second order and oscillates on exactly the step
// change in surface temperature a fire is, and the roadmap's milestone is a
// bulkhead under a fire and not a smooth ramp. The order is *asserted* in the
// tests rather than assumed, because a scheme that quietly came out first order
// when it was meant to be second is a defect a tolerance would hide.
//
// --- The factorisation is reused ---------------------------------------------
//
// `C/dt + K` does not depend on the state when the properties do not, so with a
// fixed step it is factored **once** and every subsequent step is a pair of
// triangular solves. `factorisations()` reports how many times it was actually
// factored, so a caller can tell whether it is paying `n b^2` per step or `n b`.
// One scalar unknown per node against the mechanical solver's three makes the
// band a third as wide and the factorisation a ninth of the cost.
//
// --- Temperature-dependent properties, and the secant ---------------------------
//
// With `Problem::temperatureDependent` the system is nonlinear and is closed by
// Picard iteration: re-evaluate `k` and `c` at the current iterate, re-assemble,
// re-solve, until the temperature stops moving. `iterations()` reports the last
// step's count.
//
// The capacity used is the **secant**
//
//     c = (h(T1) - h(T0)) / (T1 - T0)
//
// at each Gauss point, and not `c(T1)` or `c(T0)` or `c` at the midpoint. Two
// things follow, and they are the reason for it:
//
//   * `1^T C (T1 - T0)` is then *identically* the enthalpy change, so the energy
//     account above still closes -- to the Picard tolerance, which is the only
//     thing left between the accepted state and the one the last system was solved
//     for. A tangent `c` leaves a residual of order `dt` per step instead, which
//     accumulates and which no tolerance on the iteration would shrink.
//   * The 735 C spike is integrated rather than sampled. A tangent evaluated at
//     either end of a step that crosses the spike sees `c ~ 700` where the true
//     average over the interval is several thousand, so the plate walks through
//     the phase change at six times the rate it should. The secant cannot miss it,
//     because it *is* the average.
//
// A tangent would also be the wrong object for the same reason `plasticity.hpp`'s
// softening wanted a secant: the quantity being conserved is an integral, and the
// derivative at a point is not it.
class Solver {
public:
    // Build the forms, choose a node numbering and take the initial state.
    // Returns false with a reason on an inverted element, a malformed problem, or
    // a mesh that has no nodes.
    bool prepare(const Problem& problem, double uniformTemperature, std::string* why = nullptr);
    bool prepare(const Problem& problem, const std::vector<double>& initialTemperature,
                 std::string* why = nullptr);

    // Advance one step. False on a step that is not positive, a system that will
    // not factor, or a Picard iteration that did not converge -- the state is
    // still advanced in the last case, the same choice
    // `solidshell::elementPlasticUpdate` makes, because a caller that has to keep
    // stepping is better served by a slightly wrong state it is told about.
    bool step(double timestep, std::string* why = nullptr);

    // Drop the capacity and solve `K T = f` directly. Singular, and refused, when
    // nothing holds the temperature -- no prescribed node and no convective film.
    bool solveSteady(std::string* why = nullptr);

    // Move one film's scalars without re-preparing. `index` is into
    // `Problem::film` as it was handed to `prepare`; the face list, and therefore
    // every surface integral cached from it, is untouched.
    //
    // **This is what a coupled fire needs and nothing else here provides.** A gas
    // temperature moves every coupling step while the bulkhead it stands against
    // does not, so re-preparing would rebuild the forms, renumber the mesh and
    // reset the account -- and an account that restarts every step cannot say
    // whether the run conserved anything. False on an index that is not there,
    // rather than growing the vector: a film with no faces heats nothing, so a
    // caller who mis-indexed would otherwise see silence.
    //
    // The factorisation is dropped, because `h` is in the matrix and `h * ambient`
    // is in the load. That costs nothing under `temperatureDependent`, where Picard
    // refactors every iteration anyway.
    bool setFilm(std::size_t index, double flux, double coefficient, double ambient);

    const std::vector<double>& temperature() const { return temperature_; }
    const Account& account() const { return account_; }

    std::size_t halfBandwidth() const { return band_; }
    std::size_t freeNodes() const { return free_; }
    int factorisations() const { return factorisations_; }
    int iterations() const { return iterations_; }
    double time() const { return time_; }

    // Picard controls. `tolerance` is in kelvin on the largest nodal change
    // between iterations.
    void setPicard(double tolerance, int maximumIterations) {
        picardTolerance_ = tolerance;
        picardLimit_ = maximumIterations;
    }

private:
    // One face of the boundary a `Film` acts on, with its surface integrals taken
    // once. `mass` is `integral N_a N_b dA` and `load` is `integral N_a dA`; both
    // are geometry, so a convective boundary costs one 4x4 scatter per step and no
    // quadrature at all.
    struct FilmFace {
        std::uint32_t node[4]{};
        std::uint32_t film = 0;
        double mass[16]{};
        double load[4]{};
    };

    void number();
    void refreshProperties(const std::vector<double>& evaluateAt,
                           const std::vector<double>& previous, bool secant);
    bool buildAndSolve(double inverseStep, bool useCapacity, bool refactor, std::string* why);
    void residual(double inverseStep, bool useCapacity, std::vector<double>& out) const;
    double enthalpyOf(const std::vector<double>& field) const;
    double specificEnthalpy(double kelvin) const;

    Problem problem_{};
    std::vector<Forms> forms_;
    std::vector<double> conductance_;  // 64 per element
    std::vector<double> capacity_;     // 64 per element
    std::vector<FilmFace> filmFace_;
    std::vector<double> load_;         // full length: volumetric source + film supply
    std::vector<double> temperature_;
    std::vector<double> previous_;
    std::vector<std::ptrdiff_t> map_;   // node -> free slot, or -1
    std::vector<double> rhs_;
    std::vector<double> work_;
    solidshell::BandedSpd system_{0, 0};
    Account account_{};

    std::size_t nodes_ = 0, elements_ = 0, free_ = 0, band_ = 0;
    double time_ = 0, factoredStep_ = 0, enthalpyStart_ = 0;
    double picardTolerance_ = 1e-8;
    int picardLimit_ = 40, factorisations_ = 0, iterations_ = 0;
    bool factored_ = false, propertiesFresh_ = false, ready_ = false;
};

}  // namespace sim::thermal
