// SPDX-License-Identifier: MIT
//
// J2 (von Mises) plasticity by radial return, and the ductile failure criterion
// that decides an element has torn -- the constitutive half of
// `docs/02-simulation.md` section 3. `solid_shell.{hpp,cpp}` is the element that
// deforms; without this it springs back, and nothing decides which panels fail.
//
// Nothing here knows what an element is. It is a strain-driven point law: hand it
// the total strain at an integration point and the history that point has
// accumulated, and it returns the stress and advances the history. That is the
// interface an explicit solver wants, and it is also why almost every assertion
// in `tests/test_plasticity.cpp` is an exact identity rather than a comparison
// against a reference run.
//
// --- Why radial return, and what it buys ---------------------------------------
//
// The return map is **backward Euler on the flow rule**, solved exactly rather
// than stepped. For radial loading -- a deviatoric strain path that keeps its
// direction, which is what uniaxial tension, pure shear and every proportional
// load are -- the consistency condition at the end of a step,
//
//     2 mu (E - P) = sqrt(2/3) sigma_y(sqrt(2/3) P)
//
// determines the plastic flow P from the total deviatoric strain E and from
// nothing else. **The answer therefore does not depend on how many steps were
// taken to get there**, not approximately but to rounding, and that holds for a
// nonlinear hardening curve as well because the scalar consistency equation is
// solved to convergence inside every step. A forward-Euler update has no such
// property: it drifts off the yield surface and the drift is proportional to the
// step. That difference is the single strongest test available on this file, and
// `testStepIndependence` spends it.
//
// --- Hardening: isotropic, and kinematic as well -------------------------------
//
// Isotropic hardening grows the yield surface; kinematic hardening translates it.
// Both are here (linear Prager kinematic), because the second is exactly testable
// in a way the first is not: under **pure kinematic hardening the elastic span on
// a reversal is exactly 2 sigma_y0**, whatever the prestrain, while under pure
// isotropic hardening it is exactly 2 sigma_y(eps_p). Those are two closed forms
// that differ, so a suite that checks both cannot be satisfied by an
// implementation that has quietly mislabelled one as the other -- and the flow
// direction is taken from the *relative* stress sigma - alpha, which is a place
// an implementation can be wrong without any monotonic test noticing.
//
// It is defaulted **off** (`kinematicModulus = 0`), for two reasons. There is no
// measurement that sets it: `StructuralMaterial` carries density, E, nu and a
// yield strength and nothing else, so a kinematic modulus would be a number
// invented here. And it costs six doubles of state per integration point, which
// at the Tier-2 budget of 10^5-10^6 elements is 38-380 MB of memory whose only
// purpose is a Bauschinger effect that a predominantly monotonic ram does not
// exercise. Ship collision practice uses isotropic hardening; kinematic hardening
// earns its keep in cyclic ratcheting and fatigue, which is `Slow damage`, not
// this. The path exists and is tested, so switching it on is a material
// parameter rather than a rewrite.
//
// **Rate dependence is deliberately absent.** Steel is 10-30% stronger at
// collision strain rates and `docs/02-simulation.md` section 3 plans Johnson-Cook
// for it. Adding it makes the map viscoplastic -- the flow stress becomes a
// function of dgamma/dt -- and that destroys the step-size independence above,
// which is this file's best instrument. It deserves its own session and its own
// closed form (a fixed strain rate must reproduce the rate multiplier exactly).
// Until then the flow stress is the quasi-static one, so the model
// **under-predicts the resisting force and over-predicts the penetration**.
//
// SI units throughout, per CLAUDE.md.
#pragma once

#include <limits>

namespace sim::plasticity {

// --- Conventions ---------------------------------------------------------------
//
// Six-component Voigt vectors ordered **[xx, yy, zz, xy, yz, zx]**, matching
// `solid_shell.cpp` so a strain goes straight from `B u` into `update`.
//
// **Stress entries are tensor components; strain entries 3..5 are engineering
// shears** (gamma_xy = 2 eps_xy). The asymmetry is not an oversight: it is what
// makes plastic work a plain dot product, `sigma : d eps = sum_i s[i] de[i]`, with
// no factors of two at the call site. The price is that every *norm* has to
// remember it --
//
//     ||dev sigma||^2 = sum_{i<3} s_i^2 + 2 sum_{i>=3} s_i^2
//
// -- and dropping that 2 gives a yield surface that is no longer isotropic while
// still passing every uniaxial test. `testRotationalInvariance` is aimed at
// exactly that.
inline constexpr int kVoigt = 6;

// Passed as the failure strain when failure is not wanted: no damage accumulates
// and `State::failed` never becomes true. Plasticity tests pass it so that they
// measure the flow rule alone.
inline constexpr double kNeverFails = std::numeric_limits<double>::infinity();

// --- The hardening curve -------------------------------------------------------
//
// `Linear` exists because its return map is *affine* in dgamma, so Newton lands on
// the root in one step and the closed form for the whole uniaxial response is a
// line -- which is what makes the elementary checks exact instead of tolerant.
// `Swift` exists because real steel is not linear, and because it hands over a
// closed form the linear law cannot: **Considere's criterion**, dsigma/deps =
// sigma, puts the onset of necking at eps_p = n - eps_0 exactly. That is where
// `Failure::uniformStrain` comes from, so the mesh-regularised failure strain
// below is derived from the hardening curve rather than tabulated beside it.
//
// Voce saturation hardening is not here. It is the other curve steel is commonly
// fitted with, and adding it is a `flowStress`/`flowSlope` pair plus a numerically
// solved Considere point -- but the numerical Considere point is the reason to
// wait: it would be the first number in this file that is not a closed form.
enum class Hardening { Linear, Swift };

struct FlowCurve {
    Hardening kind = Hardening::Linear;

    double yieldStrength = 355.0e6;  // sigma_y(0), Pa. Both laws return this at eps_p = 0.

    // Linear: sigma_y = sigma_y0 + H eps_p.
    double hardeningModulus = 0.0;  // H, Pa

    // Swift: sigma_y = K (eps_0 + eps_p)^n.
    double strengthCoefficient = 0.0;  // K, Pa
    double referenceStrain = 0.0;      // eps_0
    double hardeningExponent = 0.0;    // n

    // Prager kinematic hardening: the back stress grows as (2/3) H_kin d eps_p.
    // Zero by default -- see the header note.
    double kinematicModulus = 0.0;  // H_kin, Pa
};

FlowCurve linearHardening(double yieldStrength, double hardeningModulus);

// Swift from its own parameters. `referenceStrain` is the offset that stands in
// for the yield plateau; it is *not* sigma_y/E. Taking it as the elastic strain at
// yield -- the construction that looks natural -- puts the whole curve 30% high at
// 0.2 plastic strain, because a Hollomon power law forced through the yield point
// of a steel that has a plateau is far too steep. Use `swiftFromTensile` unless
// you have K, eps_0 and n from a fit.
FlowCurve swiftHardening(double yieldStrength, double referenceStrain, double exponent);

// Swift fitted to the three numbers a tensile test reports: the yield strength,
// the **true** stress at the ultimate load, and the **true** uniform strain
// (eps_u = ln(1 + e_u), sigma_u_true = sigma_u_eng (1 + e_u)). The fit is exact:
// the returned curve passes through both points and its Considere necking strain
// is eps_u to rounding.
//
// Requires sigma_u > sigma_y exp(eps_u), which every structural steel satisfies
// comfortably; below that no Swift curve has both the strength ratio and the
// necking strain, and the fit degrades to the closest one that has the ratio.
FlowCurve swiftFromTensile(double yieldStrength, double ultimateTrueStress,
                           double uniformTrueStrain);

double flowStress(const FlowCurve& curve, double equivalentPlasticStrain);  // sigma_y, Pa
double flowSlope(const FlowCurve& curve, double equivalentPlasticStrain);   // d sigma_y / d eps_p

// Considere: the plastic strain at which the hardening rate falls to the flow
// stress and deformation stops being uniform. Zero if the curve never satisfies
// it (a linear curve with H <= sigma_y0 necks at the first increment).
double uniformElongation(const FlowCurve& curve);

// --- Ductile failure, and the element size that is not a detail -----------------
//
// A steel element does not tear at a material constant. Everything after the
// Considere point happens inside a neck whose width is set by the plate thickness,
// not by the mesh, so an element of in-plane size l that contains that neck reads
// an *average* strain -- the uniform part, plus the necking part diluted by t/l:
//
//     eps_f(l) = eps_uniform + (eps_fracture - eps_uniform) * min(t/l, 1)
//
// Both constants are properties with a definition and a measurement. `eps_uniform`
// is Considere's necking strain, which `uniformElongation` computes from the
// hardening curve. `eps_fracture` is the local true strain at fracture, ln(A0/Af)
// from the reduction of area. The clamp at l = t is not a fudge: an element no
// larger than the plate is thick resolves the neck itself and should see the local
// fracture strain, which is exactly what the formula returns there. So the
// expression interpolates between two measurable end points rather than between
// two fitted ones.
//
// What is mesh-invariant is the **necking elongation**, (eps_f(l) - eps_g) * l =
// eps_e * t, and that identity is asserted rather than described. The limit is
// stated in `docs/02-simulation.md` section 3 and is real: this regularises a tear
// that localises into one row of elements. It does not make the elongation of a
// uniformly strained gauge length mesh-independent, because with hardening and no
// softening the strain never localises in the first place.
//
// Triaxiality is the second half of "depending on the strain state". Voids grow at
// a rate proportional to exp(3 sigma_m / 2 sigma_eq) -- Rice-Tracey -- so the
// failure strain falls exponentially with the triaxiality eta = sigma_m/sigma_eq.
// The reference is eta = 1/3, uniaxial tension, which is where a tensile test
// measures the two constants above; in-plane biaxial tension is eta = 2/3 and
// fails at exp(-1/2) = 0.607 of it. Below `cutoffTriaxiality` voids close rather
// than grow and no damage accumulates at all, which also bounds the multiplier:
// with the defaults it never exceeds e.
struct Failure {
    double uniformStrain = 0.0;   // eps_g, from Considere. Mesh independent.
    double fractureStrain = 0.0;  // eps_f local, ln(A0/Af) from the reduction of area.
    double triaxialitySensitivity = 1.5;      // Rice-Tracey
    double referenceTriaxiality = 1.0 / 3.0;  // uniaxial tension
    double cutoffTriaxiality = -1.0 / 3.0;    // below this, no damage
};

// The failure strain this element size sees. `elementLength` is the in-plane
// characteristic length and `thickness` the plate thickness; `solidshell::
// elementSize` measures both off the element's own geometry, which is the point --
// the failure strain is a property of the element, not of the material.
double regularisedFailureStrain(const Failure& failure, double elementLength,
                                double thickness);

// The Rice-Tracey multiplier on the failure strain. Infinite below the cutoff,
// meaning no damage accumulates at all.
double triaxialityFactor(const Failure& failure, double triaxiality);

// --- Material ------------------------------------------------------------------
//
// Separate from `scantlings.hpp`'s `StructuralMaterial`, which carries density, E,
// nu and a yield strength and no hardening curve, no fracture strain and no
// reduction of area. That is the material database `docs/02-simulation.md`
// section 3 plans and does not yet have; until it exists the plastic properties
// cannot come out of a ship file. `shipSteel()` and `StructuralMaterial`'s elastic
// constants are asserted against each other so they cannot drift apart silently.
struct Material {
    double youngsModulus = 206.0e9;  // Pa
    double poissonRatio = 0.30;
    FlowCurve flow;
    Failure failure;

    double shearModulus() const { return youngsModulus / (2.0 * (1.0 + poissonRatio)); }
    double bulkModulus() const { return youngsModulus / (3.0 * (1.0 - 2.0 * poissonRatio)); }
};

// AH36 higher-tensile shipbuilding steel. Every constant has a provenance:
//   E, nu, sigma_y   from `scantlings.hpp`'s ah36Steel()
//   sigma_uts        490 MPa, the bottom of AH36's specified 490-620 MPa band
//   e_u              0.16 engineering uniform elongation
//   reduction of area 0.55
// The last two are representative rather than specified -- a grade certificate
// gives total elongation, not uniform elongation, and reports reduction of area
// only sometimes -- so they are the numbers to replace first when real coupon data
// arrives. They move the failure strain and nothing else.
Material shipSteel();

// --- What a yielded point has left, as an elastic modulus ------------------------
//
// A linear model of a region that has flowed needs one number from this file: how
// much stiffness is left. There are two candidates and they are **not** two
// approximations to the same thing -- they answer different questions, and
// `coupling.hpp` §5 measures which question a reduced model is asking.
//
// **Secant.** The modulus that reproduces the *total* strain a point is carrying
// from the stress it is carrying, sigma = E_s eps. In uniaxial terms the total
// strain at flow is eps = sigma_y(eps_p)/E + eps_p, so
//
//     1/E_s = 1/E + eps_p / sigma_y(eps_p)
//
// exactly. Plastic flow is deviatoric -- the return map never touches the trace --
// so the bulk modulus is untouched and the whole of the softening is in the shear
// modulus. Using 1/E = 1/(9K) + 1/(3G), the pair (K unchanged, G_s) delivers that
// E_s identically when
//
//     1/G_s = 1/G + 3 eps_p / sigma_y(eps_p)
//
// which is what `secantShearModulus` returns. That identity is the test: an
// isotropic pair built from (K, G_s) has to come back with exactly the uniaxial
// secant, not nearly.
//
// At `eps_p == 0` it returns G **exactly**, and by an early return rather than by
// the arithmetic: `1/(1/G)` is not `G` for every double, only for most of them.
// That matters because a caller builds a stiffness correction out of the
// difference, and a difference that is only nearly zero on an unyielded element is
// a coupling that has stopped being exact where it used to be. It is one line, and
// a mutant that deleted it survived the whole suite until the test below swept
// moduli the reciprocal round trip loses.
//
// **Tangent.** The modulus that relates an *increment* of stress to an increment
// of strain, 1/E_t = 1/E + 1/H' with H' = `flowSlope`, and by the same identity
//
//     1/G_t = 1/G + 3 / H'
//
// It is the modulus of a point that **is** flowing. An elastic point's tangent is
// G, and G is *not* the limit of this as eps_p goes to zero -- H' is finite there,
// so G_t drops by a finite step the instant a point flows at all, where G_s leaves
// G continuously. Which points are flowing is the caller's to know; this returns
// the flowing value at any argument, including zero, so that the size of that step
// can be read off directly.
// For AH36 at the first increment of flow that step is a factor of 34 -- G_t/G =
// 0.0296 against G_s/G = 0.999999, measured in `tests/test_plasticity.cpp`. A model
// that has to reproduce a *total* displacement under a *total* load and reaches for
// the tangent is therefore not making a small error at small plastic strain; it is
// making its largest error there.
//
// Both clamp at zero: a perfectly plastic curve has H' = 0 and no tangent
// stiffness at all, and a curve with no strength has no secant one.
double secantShearModulus(const Material& material, double equivalentPlasticStrain);
double tangentShearModulus(const Material& material, double equivalentPlasticStrain);

// The uniaxial partner of the two above, exposed because it is the closed form the
// isotropic pair is checked against rather than a second route to the same answer.
double secantYoungsModulus(const Material& material, double equivalentPlasticStrain);

// (E, nu) of the isotropic pair with these moduli -- the inverse of
// `Material::bulkModulus` and `shearModulus`, and how a softened shear modulus
// becomes something `StructuralMaterial` and therefore `solidshell::
// elementStiffness` can carry. Not a member of `Material` because the pair a
// caller softens belongs to the *elastic* material its stiffness was assembled
// from, which is a `StructuralMaterial` and has no flow curve.
void isotropicFromBulkShear(double bulkModulus, double shearModulus, double* youngsModulus,
                            double* poissonRatio);

// --- State ---------------------------------------------------------------------
//
// One per integration point. Trivially copyable on purpose: probing the law with a
// candidate strain, which the tests and any equilibrium iteration both need, is
// done by copying the state and throwing the copy away.
struct State {
    double plasticStrain[kVoigt] = {};  // engineering shear convention; deviatoric
    double backStress[kVoigt] = {};     // Pa, deviatoric
    double equivalentPlasticStrain = 0.0;
    double damage = 0.0;  // 0..1, monotone. 1 means torn.
    bool failed = false;
};

// What one increment did. Everything here is diagnostic except `failedNow`, which
// is how a caller learns a hole has opened.
struct Increment {
    bool yielded = false;
    bool failedNow = false;
    double plasticMultiplier = 0.0;         // dgamma
    double equivalentPlasticStrainIncrement = 0.0;  // sqrt(2/3) dgamma
    double dissipation = 0.0;   // J/m^3, (sigma - alpha) : d eps_p. Never negative.
    double plasticWork = 0.0;   // J/m^3, sigma : d eps_p. May be negative on a
                                // reversal under kinematic hardening -- the back
                                // stress is giving stored energy back.
    double vonMises = 0.0;      // Pa, of the returned stress
    double triaxiality = 0.0;   // sigma_m / sigma_eq of the returned stress
    // Evaluations of the consistency condition. One for a perfectly plastic curve,
    // where the no-hardening bracket *is* the answer; exactly two for a linear one,
    // where g is affine so a single Newton step lands on the root and the second
    // evaluation only confirms it; a handful for Swift.
    int iterations = 0;
    double failureNormal[3] = {0.0, 0.0, 0.0};  // max principal direction, when failedNow
};

// The whole law. Advances `state` from its stored history to `strain`, writes the
// stress, and optionally the **algorithmic** (consistent) tangent -- 6x6 row-major,
// the derivative of this very map, which is what makes an element-level
// equilibrium iteration converge quadratically. The tangent is checked against a
// finite difference of `update` itself, so it cannot be a plausible-looking wrong
// formula.
//
// `failureStrain` is the mesh-regularised value from `regularisedFailureStrain`,
// passed rather than looked up because it belongs to the element and not to the
// material. Use `kNeverFails` to exercise the flow rule alone.
//
// Once `state.failed` is set the point carries **no stress at all** and no further
// history: failure is irreversible and unloading does not heal it. That is element
// deletion, the cheap fallback `docs/02-simulation.md` section 3 names; the stress
// drops discontinuously, which an explicit scheme feels as a small shock.
Increment update(const Material& material, double failureStrain, const double strain[kVoigt],
                 State& state, double stress[kVoigt], double tangent[kVoigt * kVoigt] = nullptr);

// --- Exposed for the tests, and useful on their own ----------------------------

void elasticModuli(const Material& material, double c[kVoigt * kVoigt]);  // row-major 6x6
void elasticStress(const Material& material, const double strain[kVoigt], double stress[kVoigt]);

void deviator(const double stress[kVoigt], double out[kVoigt]);
double deviatorNorm(const double stress[kVoigt]);  // ||dev sigma||, the tensor norm
double vonMises(const double stress[kVoigt]);      // sqrt(3/2) ||dev sigma||
double meanStress(const double stress[kVoigt]);    // tr(sigma)/3
double triaxiality(const double stress[kVoigt]);   // mean / von Mises, 0 if the latter is 0

// Unit eigenvector of the largest principal stress, sign-fixed so the largest
// component is positive. This is the plane a tear would open on, which is what a
// splitting fracture model needs and what element deletion throws away.
void maxPrincipalDirection(const double stress[kVoigt], double out[3]);

}  // namespace sim::plasticity
