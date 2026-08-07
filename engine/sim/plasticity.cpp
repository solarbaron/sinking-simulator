// SPDX-License-Identifier: MIT
#include "plasticity.hpp"

#include <algorithm>
#include <cmath>

namespace sim::plasticity {
namespace {

// sqrt(2/3) turns a plastic multiplier into an equivalent plastic strain and a
// yield strength into a deviatoric radius. It appears often enough, and wrongly
// often enough, to be named once.
const double kRoot23 = std::sqrt(2.0 / 3.0);

// The tensor inner product of two Voigt vectors where one carries engineering
// shears and the other carries tensor components. sigma : d eps is a plain dot
// product precisely because of that pairing -- see the header.
double dot(const double a[kVoigt], const double b[kVoigt]) {
    double s = 0.0;
    for (int i = 0; i < kVoigt; ++i) s += a[i] * b[i];
    return s;
}

// Squared Frobenius norm of a *tensor-component* Voigt vector. The off-diagonals
// appear twice in the tensor, hence the 2.
double tensorNormSquared(const double v[kVoigt]) {
    return v[0] * v[0] + v[1] * v[1] + v[2] * v[2] +
           2.0 * (v[3] * v[3] + v[4] * v[4] + v[5] * v[5]);
}

// Jacobi eigen-decomposition of the symmetric 3x3 a stress Voigt vector encodes.
// Small and unconditionally convergent; the alternative closed form has a
// degenerate branch at repeated eigenvalues, which is exactly the case an
// equibiaxial plate sits in.
void symmetricEigen(const double stress[kVoigt], double values[3], double vectors[3][3]) {
    double a[3][3] = {{stress[0], stress[3], stress[5]},
                      {stress[3], stress[1], stress[4]},
                      {stress[5], stress[4], stress[2]}};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) vectors[i][j] = (i == j) ? 1.0 : 0.0;

    for (int sweep = 0; sweep < 64; ++sweep) {
        double off = a[0][1] * a[0][1] + a[0][2] * a[0][2] + a[1][2] * a[1][2];
        if (off <= 0.0) break;
        for (int p = 0; p < 2; ++p)
            for (int q = p + 1; q < 3; ++q) {
                const double apq = a[p][q];
                if (apq == 0.0) continue;
                const double theta = (a[q][q] - a[p][p]) / (2.0 * apq);
                const double t = (theta >= 0.0 ? 1.0 : -1.0) /
                                 (std::abs(theta) + std::sqrt(theta * theta + 1.0));
                const double c = 1.0 / std::sqrt(t * t + 1.0), s = t * c;
                for (int k = 0; k < 3; ++k) {
                    const double akp = a[k][p], akq = a[k][q];
                    a[k][p] = c * akp - s * akq;
                    a[k][q] = s * akp + c * akq;
                }
                for (int k = 0; k < 3; ++k) {
                    const double apk = a[p][k], aqk = a[q][k];
                    a[p][k] = c * apk - s * aqk;
                    a[q][k] = s * apk + c * aqk;
                }
                for (int k = 0; k < 3; ++k) {
                    const double vkp = vectors[k][p], vkq = vectors[k][q];
                    vectors[k][p] = c * vkp - s * vkq;
                    vectors[k][q] = s * vkp + c * vkq;
                }
            }
    }
    for (int i = 0; i < 3; ++i) values[i] = a[i][i];
}

}  // namespace

// --- Hardening curves ----------------------------------------------------------

FlowCurve linearHardening(double yieldStrength, double hardeningModulus) {
    FlowCurve curve;
    curve.kind = Hardening::Linear;
    curve.yieldStrength = yieldStrength;
    curve.hardeningModulus = hardeningModulus;
    return curve;
}

FlowCurve swiftHardening(double yieldStrength, double referenceStrain, double exponent) {
    FlowCurve curve;
    curve.kind = Hardening::Swift;
    curve.yieldStrength = yieldStrength;
    curve.referenceStrain = referenceStrain;
    curve.hardeningExponent = exponent;
    // K is fixed by sigma_y(0) = sigma_y0, so the curve starts where the elastic
    // solution left off and there is no step at first yield.
    curve.strengthCoefficient = yieldStrength / std::pow(referenceStrain, exponent);
    return curve;
}

FlowCurve swiftFromTensile(double yieldStrength, double ultimateTrueStress,
                           double uniformTrueStrain) {
    // With n = eps_u + eps_0 forced by Considere, the strength ratio is
    //     sigma_u / sigma_y = (n / eps_0)^n,
    // whose logarithm falls monotonically from infinity (eps_0 -> 0) to eps_u
    // (eps_0 -> infinity). One bisection, no derivatives, no starting guess.
    const double target = std::log(ultimateTrueStress / yieldStrength);
    const auto ratio = [&](double eps0) {
        const double n = uniformTrueStrain + eps0;
        return n * std::log(n / eps0);
    };
    double lo = 1e-9, hi = 1e3;
    for (int i = 0; i < 200; ++i) {
        // Geometric midpoint, because eps_0 is bracketed across twelve decades and
        // this converges on its *order* first. It is conditioning, not correctness:
        // mutation testing confirmed that 200 arithmetic bisections reach the same
        // root, since 1e3/2^200 is below anything double precision can distinguish.
        const double mid = std::sqrt(lo * hi);
        if (ratio(mid) > target)
            lo = mid;
        else
            hi = mid;
    }
    const double eps0 = std::sqrt(lo * hi);
    return swiftHardening(yieldStrength, eps0, uniformTrueStrain + eps0);
}

double flowStress(const FlowCurve& curve, double equivalentPlasticStrain) {
    const double p = std::max(0.0, equivalentPlasticStrain);
    switch (curve.kind) {
        case Hardening::Linear: return curve.yieldStrength + curve.hardeningModulus * p;
        case Hardening::Swift: break;
    }
    return curve.strengthCoefficient *
           std::pow(curve.referenceStrain + p, curve.hardeningExponent);
}

double flowSlope(const FlowCurve& curve, double equivalentPlasticStrain) {
    const double p = std::max(0.0, equivalentPlasticStrain);
    switch (curve.kind) {
        case Hardening::Linear: return curve.hardeningModulus;
        case Hardening::Swift: break;
    }
    return curve.strengthCoefficient * curve.hardeningExponent *
           std::pow(curve.referenceStrain + p, curve.hardeningExponent - 1.0);
}

double uniformElongation(const FlowCurve& curve) {
    switch (curve.kind) {
        case Hardening::Linear:
            // sigma_y0 + H eps = H  =>  eps = 1 - sigma_y0/H, and nothing at all if
            // the curve never hardens as fast as it is strong.
            if (!(curve.hardeningModulus > curve.yieldStrength)) return 0.0;
            return 1.0 - curve.yieldStrength / curve.hardeningModulus;
        case Hardening::Swift: break;
    }
    return std::max(0.0, curve.hardeningExponent - curve.referenceStrain);
}

// --- Failure -------------------------------------------------------------------

double regularisedFailureStrain(const Failure& failure, double elementLength,
                                double thickness) {
    // An element no larger than the plate is thick contains the neck rather than
    // averaging over it, so it sees the local fracture strain; a bigger one sees
    // the neck diluted by t/l. Degenerate geometry falls to the unregularised end
    // rather than to zero, because reporting "fails instantly" for a bad Jacobian
    // would be a failure-open.
    double share = 1.0;
    if (elementLength > 0.0 && thickness > 0.0)
        share = std::min(1.0, thickness / elementLength);
    return failure.uniformStrain + (failure.fractureStrain - failure.uniformStrain) * share;
}

double triaxialityFactor(const Failure& failure, double eta) {
    if (eta <= failure.cutoffTriaxiality) return std::numeric_limits<double>::infinity();
    return std::exp(-failure.triaxialitySensitivity * (eta - failure.referenceTriaxiality));
}

// --- Material ------------------------------------------------------------------

Material shipSteel() {
    Material material;
    material.youngsModulus = 206.0e9;
    material.poissonRatio = 0.30;

    // Engineering 490 MPa at 0.16 elongation, converted to true measures.
    const double engineeringUltimate = 490.0e6, engineeringUniform = 0.16;
    material.flow = swiftFromTensile(355.0e6,
                                     engineeringUltimate * (1.0 + engineeringUniform),
                                     std::log(1.0 + engineeringUniform));

    // The two ends of the regularisation, both measurable: Considere's necking
    // strain and the true strain at fracture from a 55% reduction of area.
    material.failure.uniformStrain = uniformElongation(material.flow);
    material.failure.fractureStrain = std::log(1.0 / (1.0 - 0.55));
    return material;
}

// --- What a yielded point has left ---------------------------------------------

double secantShearModulus(const Material& material, double equivalentPlasticStrain) {
    const double shear = material.shearModulus();
    const double plastic = std::max(0.0, equivalentPlasticStrain);
    // **The early return is what makes this exact, and the arithmetic below is
    // not.** `1/(1/G) == G` for most doubles -- including AH36's shear modulus,
    // which is why deleting this line survived the first round of mutation testing
    // -- but not for all of them, and one unit in the last place turns a caller's
    // knockdown ratio from 1 into 0.999...89 and puts a block of near-zeros where
    // there should be no block at all. See the header, and the modulus sweep in
    // `tests/test_plasticity.cpp` that now carries a case it fails.
    if (plastic == 0.0) return shear;
    const double strength = flowStress(material.flow, plastic);
    if (!(strength > 0.0)) return 0.0;
    const double compliance = 1.0 / shear + 3.0 * plastic / strength;
    return compliance > 0.0 ? 1.0 / compliance : 0.0;
}

double tangentShearModulus(const Material& material, double equivalentPlasticStrain) {
    const double shear = material.shearModulus();
    const double slope = flowSlope(material.flow, std::max(0.0, equivalentPlasticStrain));
    // Perfect plasticity: no incremental stiffness at all. The limit of the
    // expression below, taken rather than divided by.
    if (!(slope > 0.0)) return 0.0;
    return 1.0 / (1.0 / shear + 3.0 / slope);
}

double secantYoungsModulus(const Material& material, double equivalentPlasticStrain) {
    const double plastic = std::max(0.0, equivalentPlasticStrain);
    if (plastic == 0.0) return material.youngsModulus;
    const double strength = flowStress(material.flow, plastic);
    if (!(strength > 0.0)) return 0.0;
    return 1.0 / (1.0 / material.youngsModulus + plastic / strength);
}

void isotropicFromBulkShear(double bulk, double shearModulus, double* youngsModulus,
                            double* poissonRatio) {
    const double denominator = 3.0 * bulk + shearModulus;
    if (youngsModulus)
        *youngsModulus = denominator != 0.0 ? 9.0 * bulk * shearModulus / denominator : 0.0;
    if (poissonRatio)
        *poissonRatio =
            denominator != 0.0 ? (3.0 * bulk - 2.0 * shearModulus) / (2.0 * denominator) : 0.5;
}

// --- Stress algebra ------------------------------------------------------------

void elasticModuli(const Material& material, double c[kVoigt * kVoigt]) {
    const double e = material.youngsModulus, nu = material.poissonRatio;
    const double lambda = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    const double mu = material.shearModulus();
    for (int i = 0; i < kVoigt * kVoigt; ++i) c[i] = 0.0;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) c[i * kVoigt + j] = lambda;
        c[i * kVoigt + i] += 2.0 * mu;
        c[(3 + i) * kVoigt + (3 + i)] = mu;
    }
}

// Split into a pressure and a deviator rather than written with Lame's lambda, so
// that this and the elastic branch of `update` are the same arithmetic in the same
// order and agree bit for bit. Against `elasticModuli`, which does use lambda, they
// agree only to rounding -- which is itself a check worth having.
void elasticStress(const Material& material, const double strain[kVoigt], double stress[kVoigt]) {
    const double mu = material.shearModulus();
    const double volumetric = strain[0] + strain[1] + strain[2];
    const double pressure = material.bulkModulus() * volumetric;
    for (int i = 0; i < 3; ++i) {
        stress[i] = pressure + 2.0 * mu * (strain[i] - volumetric / 3.0);
        stress[3 + i] = mu * strain[3 + i];
    }
}

void deviator(const double stress[kVoigt], double out[kVoigt]) {
    const double mean = (stress[0] + stress[1] + stress[2]) / 3.0;
    for (int i = 0; i < 3; ++i) {
        out[i] = stress[i] - mean;
        out[3 + i] = stress[3 + i];
    }
}

double deviatorNorm(const double stress[kVoigt]) {
    double s[kVoigt];
    deviator(stress, s);
    return std::sqrt(tensorNormSquared(s));
}

double vonMises(const double stress[kVoigt]) { return std::sqrt(1.5) * deviatorNorm(stress); }

double meanStress(const double stress[kVoigt]) {
    return (stress[0] + stress[1] + stress[2]) / 3.0;
}

double triaxiality(const double stress[kVoigt]) {
    const double q = vonMises(stress);
    return q > 0.0 ? meanStress(stress) / q : 0.0;
}

void maxPrincipalDirection(const double stress[kVoigt], double out[3]) {
    double values[3], vectors[3][3];
    symmetricEigen(stress, values, vectors);
    int best = 0;
    for (int i = 1; i < 3; ++i)
        if (values[i] > values[best]) best = i;
    for (int i = 0; i < 3; ++i) out[i] = vectors[i][best];

    double norm = std::sqrt(out[0] * out[0] + out[1] * out[1] + out[2] * out[2]);
    if (!(norm > 0.0)) {
        out[0] = 1.0;
        out[1] = out[2] = 0.0;
        return;
    }
    // Sign is arbitrary for an eigenvector; fix it so two callers agree.
    int largest = 0;
    for (int i = 1; i < 3; ++i)
        if (std::abs(out[i]) > std::abs(out[largest])) largest = i;
    if (out[largest] < 0.0) norm = -norm;
    for (int i = 0; i < 3; ++i) out[i] /= norm;
}

// --- The return map ------------------------------------------------------------

Increment update(const Material& material, double failureStrain, const double strain[kVoigt],
                 State& state, double stress[kVoigt], double tangent[kVoigt * kVoigt]) {
    Increment increment;

    // A torn point carries nothing and learns nothing. Checked before anything
    // else so that "unloading does not heal it" is structural rather than a
    // consequence of the arithmetic below.
    if (state.failed) {
        for (int i = 0; i < kVoigt; ++i) stress[i] = 0.0;
        if (tangent != nullptr)
            for (int i = 0; i < kVoigt * kVoigt; ++i) tangent[i] = 0.0;
        return increment;
    }

    const double mu = material.shearModulus();
    const double kappa = material.bulkModulus();

    // Elastic predictor. The volumetric response is elastic always -- plastic flow
    // is deviatoric -- so the pressure is final at this point and takes no part in
    // the return.
    double elastic[kVoigt];
    for (int i = 0; i < kVoigt; ++i) elastic[i] = strain[i] - state.plasticStrain[i];
    const double volumetric = elastic[0] + elastic[1] + elastic[2];
    const double pressure = kappa * volumetric;

    double trial[kVoigt];  // deviatoric trial stress
    for (int i = 0; i < 3; ++i) {
        trial[i] = 2.0 * mu * (elastic[i] - volumetric / 3.0);
        trial[3 + i] = mu * elastic[3 + i];  // sigma_xy = mu gamma_xy
    }

    double relative[kVoigt];  // trial minus back stress: what the yield surface sees
    for (int i = 0; i < kVoigt; ++i) relative[i] = trial[i] - state.backStress[i];
    const double relativeNorm = std::sqrt(tensorNormSquared(relative));

    const double yieldRadius = kRoot23 * flowStress(material.flow, state.equivalentPlasticStrain);
    const double excess = relativeNorm - yieldRadius;

    if (!(excess > 0.0)) {
        for (int i = 0; i < 3; ++i) {
            stress[i] = pressure + trial[i];
            stress[3 + i] = trial[3 + i];
        }
        if (tangent != nullptr) elasticModuli(material, tangent);
        // From the stress, not from the relative stress: with a back stress the two
        // differ and the caller wants the real invariant.
        increment.vonMises = vonMises(stress);
        increment.triaxiality = triaxiality(stress);
        return increment;
    }

    // --- Consistency: solve for the plastic multiplier ---------------------------
    //
    //   g(dg) = ||xi_trial|| - (2 mu + (2/3) H_kin) dg
    //                        - sqrt(2/3) sigma_y(eps_p + sqrt(2/3) dg) = 0
    //
    // g is strictly decreasing, g(0) = excess > 0, and dropping the hardening term
    // gives an upper bound, so the root is bracketed before the first step. It is
    // entered from that upper bound: for a **perfectly plastic** curve the bound is
    // the answer and the loop costs one evaluation; for a **linear** curve g is
    // affine, so one Newton step lands on the root exactly and the second
    // evaluation only confirms it; for Swift g is convex, so Newton approaches
    // monotonically from above. The bisection fallback is for the curve nobody has
    // added yet.
    const double kinematic = material.flow.kinematicModulus;
    const double linearPart = 2.0 * mu + (2.0 / 3.0) * kinematic;
    double lo = 0.0, hi = excess / linearPart;
    double gamma = hi;
    int iterations = 0;
    for (; iterations < 50; ++iterations) {
        const double accumulated = state.equivalentPlasticStrain + kRoot23 * gamma;
        const double residual =
            relativeNorm - linearPart * gamma - kRoot23 * flowStress(material.flow, accumulated);
        if (residual > 0.0)
            lo = gamma;
        else
            hi = gamma;
        if (std::abs(residual) <= 1e-14 * yieldRadius) break;

        const double slope = -linearPart - (2.0 / 3.0) * flowSlope(material.flow, accumulated);
        double next = gamma - residual / slope;
        if (!(next > lo && next < hi)) next = 0.5 * (lo + hi);
        // Converge on the **step** as well as on the residual. g is a difference of
        // terms of size ||xi_trial||, so its own rounding floor grows with the
        // elastic predictor and can sit above the tolerance above -- at 12% plastic
        // strain the predictor is 80x the yield radius and the floor rises with it.
        // A correction this small cannot move the stress, so grinding on the
        // residual there buys nothing and costs an unbounded iteration count. With
        // this, a linear curve always costs exactly two evaluations.
        if (std::abs(next - gamma) <= 1e-15 * gamma) {
            gamma = next;
            break;
        }
        gamma = next;
    }

    // The flow direction, unit in the tensor norm. Taken from the *relative* stress
    // -- with a back stress present, taking it from the stress instead is a wrong
    // implementation that no monotonic test can see.
    double direction[kVoigt];
    for (int i = 0; i < kVoigt; ++i) direction[i] = relative[i] / relativeNorm;

    double plasticIncrement[kVoigt];  // engineering shear convention, like the state
    for (int i = 0; i < 3; ++i) {
        plasticIncrement[i] = gamma * direction[i];
        plasticIncrement[3 + i] = 2.0 * gamma * direction[3 + i];
    }

    for (int i = 0; i < 3; ++i) {
        stress[i] = pressure + trial[i] - 2.0 * mu * gamma * direction[i];
        stress[3 + i] = trial[3 + i] - 2.0 * mu * gamma * direction[3 + i];
    }

    double back[kVoigt];
    for (int i = 0; i < kVoigt; ++i)
        back[i] = state.backStress[i] + (2.0 / 3.0) * kinematic * gamma * direction[i];

    double effective[kVoigt];
    for (int i = 0; i < kVoigt; ++i) effective[i] = stress[i] - back[i];

    increment.yielded = true;
    increment.iterations = iterations + 1;
    increment.plasticMultiplier = gamma;
    increment.equivalentPlasticStrainIncrement = kRoot23 * gamma;
    // Both by contraction rather than from the closed form sigma_y * d eps_p, so
    // the test that they agree is a test and not a tautology.
    increment.plasticWork = dot(stress, plasticIncrement);
    increment.dissipation = dot(effective, plasticIncrement);

    for (int i = 0; i < kVoigt; ++i) {
        state.plasticStrain[i] += plasticIncrement[i];
        state.backStress[i] = back[i];
    }
    state.equivalentPlasticStrain += increment.equivalentPlasticStrainIncrement;

    increment.vonMises = vonMises(stress);
    increment.triaxiality = triaxiality(stress);

    // --- Damage ------------------------------------------------------------------
    //
    // Accumulated rather than compared, so a path that wanders through several
    // stress states spends the right fraction of its life in each. A constant
    // failure strain and a monotonic path make this exactly eps_p / eps_f, which
    // is the elementary criterion and is what the failure tests assert.
    const double critical = failureStrain * triaxialityFactor(material.failure,
                                                              increment.triaxiality);
    if (critical > 0.0 && std::isfinite(critical))
        state.damage += increment.equivalentPlasticStrainIncrement / critical;

    if (state.damage >= 1.0) {
        state.damage = 1.0;
        state.failed = true;
        increment.failedNow = true;
        maxPrincipalDirection(stress, increment.failureNormal);
        for (int i = 0; i < kVoigt; ++i) stress[i] = 0.0;
        if (tangent != nullptr)
            for (int i = 0; i < kVoigt * kVoigt; ++i) tangent[i] = 0.0;
        return increment;
    }

    if (tangent == nullptr) return increment;

    // --- Algorithmic tangent ------------------------------------------------------
    //
    //   C = kappa 1(x)1 + 2 mu theta I_dev - 2 mu thetaBar n(x)n
    //
    // theta = 1 - 2 mu dg / ||xi_trial|| is the softening the finite step already
    // committed to; thetaBar carries the hardening. This is the derivative of the
    // map above, not of the continuum flow rule, and the difference is exactly the
    // theta term -- which vanishes as the step does. Checked against a central
    // finite difference of `update` in tests/test_plasticity.cpp.
    const double theta = 1.0 - 2.0 * mu * gamma / relativeNorm;
    const double plasticModulus =
        flowSlope(material.flow, state.equivalentPlasticStrain) + kinematic;
    const double thetaBar = 1.0 / (1.0 + plasticModulus / (3.0 * mu)) - (1.0 - theta);

    for (int i = 0; i < kVoigt * kVoigt; ++i) tangent[i] = 0.0;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j)
            tangent[i * kVoigt + j] = kappa + 2.0 * mu * theta * ((i == j ? 1.0 : 0.0) - 1.0 / 3.0);
        tangent[(3 + i) * kVoigt + (3 + i)] = mu * theta;
    }
    for (int i = 0; i < kVoigt; ++i)
        for (int j = 0; j < kVoigt; ++j)
            tangent[i * kVoigt + j] -= 2.0 * mu * thetaBar * direction[i] * direction[j];
    return increment;
}

}  // namespace sim::plasticity
