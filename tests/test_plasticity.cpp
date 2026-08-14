// SPDX-License-Identifier: MIT
//
// Validation of J2 plasticity, the return map, and the ductile failure criterion.
//
// Plasticity is unusually well supplied with exact answers, and this file spends
// them rather than comparing against a stored curve:
//
//   * **uniaxial tension** has a closed form for the whole response -- the axial
//     stress is the hardening curve evaluated at the plastic strain, and the
//     lateral strain is -nu sigma/E - eps_p/2, because the elastic part contracts
//     by Poisson's ratio and the plastic part is incompressible. Both are asserted,
//     the transverse stress that closed form implies must come back at machine
//     zero, and the same points are reached a second time by bisecting the return
//     map itself, using none of the closed form;
//   * **the return map is exact for radial loading.** One step to a final strain
//     and ten thousand steps to the same final strain must give the same stress.
//     That is a property backward Euler has and a forward-Euler update does not,
//     and it holds for a nonlinear hardening curve too because the scalar
//     consistency equation is solved to convergence inside each step. It is the
//     strongest single discriminator available here, so it is also given a
//     **negative control**: a deliberately non-proportional path *must* come out
//     step dependent, or the test above is measuring nothing;
//   * **the yield surface is a cylinder about the hydrostatic axis.** Pressure of
//     any magnitude produces no plastic strain at all, and superposing an arbitrary
//     wandering pressure on a deviatoric path must not move the deviatoric answer.
//     An implementation that got the deviator wrong still looks right in uniaxial
//     tension;
//   * **isotropy.** The same strain tensor, rotated, must give the rotated stress.
//     The Voigt norm carries a factor of two on its shear entries and dropping it
//     produces a yield surface that is not round -- invisible to every axis-aligned
//     load in this file, and caught by this one alone;
//   * **consistency, incompressibility and dissipation** are equalities, not
//     tolerances: after the return the stress sits on the yield surface, the
//     plastic strain increment has zero trace, and the dissipation computed by
//     tensor contraction equals sigma_y d eps_p;
//   * **the algorithmic tangent** is checked against a central finite difference of
//     the return map itself, so it cannot be a plausible-looking wrong formula.
//
// For the element hook the tie that matters is that a material which never yields
// must reproduce, to rounding, exactly what the validated elastic element already
// computes -- and that the enhanced-strain residual, which is the *definition* of
// the seven EAS parameters once the material is nonlinear, comes back at its floor.
//
// **Six of these tests exist because mutation testing said the rest were not
// enough**, and two of them found defects rather than gaps. 68 single-edit mutants;
// the first pass killed 58, and of the ten survivors six were real holes:
// `elasticStress` was reachable by no test at all and could drop its deviatoric
// split; the element's size could be taken off its bottom face, because every
// element under test was prismatic; `initialisePlasticState` could decline to clear
// the history it was initialising; an element could call itself torn on its *first*
// dead integration point, because every tearing test until then strained the
// element uniformly, where all eight die on the same step; and the failure plane's
// normal was only ever asked for on an axis-aligned pull, where the stress tensor
// is already diagonal and the eigenvectors come back as the identity however badly
// they are computed. Writing the partial-failure test then exposed the defect: with
// four of eight points gone the enhanced-strain problem loses rank, the surviving
// points were driven below the damage cutoff, and the element **stopped tearing** --
// its damage frozen at 0.78 while its plastic strain ran on from 0.49 to 0.89.
// `docs/07-fem-spike-findings.md` section 7 records the full account.
#include "engine/sim/plasticity.hpp"
#include "engine/sim/scantlings.hpp"
#include "engine/sim/solid_shell.hpp"
#include "harness.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <numbers>
#include <random>
#include <string>
#include <vector>

using namespace sim;
using namespace sim::plasticity;
using testing::expectEqual;
using testing::expectNear;
using testing::expectTrue;

namespace {

constexpr int kV = kVoigt;

// C++ has no compound literals, and a strain written out over four lines stops
// being readable. Converts to `const double*`, which is what `update` takes.
struct V6 {
    double v[kVoigt];
    explicit V6(double a = 0, double b = 0, double c = 0, double d = 0, double e = 0,
                double f = 0)
        : v{a, b, c, d, e, f} {}
    operator const double*() const { return v; }  // NOLINT: that is the point
};

Material steel() { return shipSteel(); }

// Linear hardening makes every closed form in the uniaxial tests a line, so those
// assertions are exact rather than tolerant.
Material linearSteel(double hardening = 2.0e9) {
    Material material;
    material.youngsModulus = 206.0e9;
    material.poissonRatio = 0.30;
    material.flow = linearHardening(355.0e6, hardening);
    material.failure.uniformStrain = uniformElongation(material.flow);
    material.failure.fractureStrain = 0.8;
    return material;
}

double maxAbs(const double* v, int n) {
    double s = 0.0;
    for (int i = 0; i < n; ++i) s = std::max(s, std::abs(v[i]));
    return s;
}

// Voigt <-> tensor, for the rotation tests. `engineering` says the 3..5 entries are
// twice the tensor component, which strains are and stresses are not.
void toTensor(const double voigt[kV], bool engineering, double m[3][3]) {
    const double half = engineering ? 0.5 : 1.0;
    m[0][0] = voigt[0];
    m[1][1] = voigt[1];
    m[2][2] = voigt[2];
    m[0][1] = m[1][0] = half * voigt[3];
    m[1][2] = m[2][1] = half * voigt[4];
    m[2][0] = m[0][2] = half * voigt[5];
}

void fromTensor(const double m[3][3], bool engineering, double voigt[kV]) {
    const double twice = engineering ? 2.0 : 1.0;
    voigt[0] = m[0][0];
    voigt[1] = m[1][1];
    voigt[2] = m[2][2];
    voigt[3] = twice * m[0][1];
    voigt[4] = twice * m[1][2];
    voigt[5] = twice * m[2][0];
}

void rotateVoigt(const double r[3][3], const double voigt[kV], bool engineering, double out[kV]) {
    double m[3][3], t[3][3], rotated[3][3];
    toTensor(voigt, engineering, m);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            double s = 0.0;
            for (int k = 0; k < 3; ++k) s += m[i][k] * r[j][k];
            t[i][j] = s;
        }
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            double s = 0.0;
            for (int k = 0; k < 3; ++k) s += r[i][k] * t[k][j];
            rotated[i][j] = s;
        }
    fromTensor(rotated, engineering, out);
}

// Rotation about an axis, Rodrigues. Used to demand that the yield surface is
// round, which no axis-aligned load can establish.
void axisRotation(const double axis[3], double angle, double out[3][3]) {
    const double n = std::sqrt(axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]);
    const double a = axis[0] / n, b = axis[1] / n, c = axis[2] / n;
    const double s = std::sin(angle), t = 1.0 - std::cos(angle), co = std::cos(angle);
    out[0][0] = co + a * a * t;
    out[0][1] = a * b * t - c * s;
    out[0][2] = a * c * t + b * s;
    out[1][0] = b * a * t + c * s;
    out[1][1] = co + b * b * t;
    out[1][2] = b * c * t - a * s;
    out[2][0] = c * a * t - b * s;
    out[2][1] = c * b * t + a * s;
    out[2][2] = co + c * c * t;
}

// One point driven in **uniaxial stress**: given the axial strain, the lateral
// strain that leaves the transverse stress at zero is found by bisection through
// the return map itself, not from the closed form. That is the independent route
// the closed form is checked against.
struct UniaxialProbe {
    double axialStress = 0.0;
    double lateralStrain = 0.0;
    double transverseStress = 0.0;
    State state;
    Increment increment;
};

UniaxialProbe uniaxialPull(const Material& material, double failureStrain, const State& before,
                           double axialStrain, int bisections = 80) {
    UniaxialProbe result;
    result.state = before;
    if (before.failed) return result;

    const double span = std::abs(axialStrain) + 1e-9;
    double lo = -span, hi = span;
    double stress[kV];
    for (int iteration = 0; iteration < bisections; ++iteration) {
        const double mid = 0.5 * (lo + hi);
        State probe = before;
        const double strain[kV] = {axialStrain, mid, mid, 0.0, 0.0, 0.0};
        // **kNeverFails, deliberately.** Which lateral strain nulls the transverse
        // stress is a purely mechanical question, and damage does not touch the
        // stress until the point tears, so the search gets the same answer either
        // way -- but only this way is it a search. With the real failure strain a
        // torn probe returns *zero* transverse stress, which reads as "not enough
        // lateral stretch", and the bracket walks to its far end: the first version
        // of this helper did that, produced a near-hydrostatic state at eta = 15.9,
        // and made a bar tear at a seventh of its failure strain. The commit below
        // uses the real failure strain.
        update(material, kNeverFails, strain, probe, stress);
        // Transverse stress rises with lateral stretch, so the bracket closes on
        // the root from a single monotone comparison.
        if (stress[1] > 0.0)
            hi = mid;
        else
            lo = mid;
    }
    result.lateralStrain = 0.5 * (lo + hi);
    const double strain[kV] = {axialStrain, result.lateralStrain, result.lateralStrain, 0.0, 0.0,
                               0.0};
    result.state = before;
    result.increment = update(material, failureStrain, strain, result.state, stress);
    result.axialStress = stress[0];
    result.transverseStress = std::max(std::abs(stress[1]), std::abs(stress[2]));
    return result;
}

// --- the hardening curves ------------------------------------------------------

void testFlowCurves() {
    const FlowCurve linear = linearHardening(355.0e6, 2.0e9);
    expectNear("linear curve starts at the yield strength", flowStress(linear, 0.0), 355.0e6, 1e-6);
    expectNear("linear curve is a line", flowStress(linear, 0.05), 355.0e6 + 0.1e9, 1e-3);
    expectNear("linear slope is H", flowSlope(linear, 0.37), 2.0e9, 1e-6);

    // Considere for a line: sigma_y0 + H eps = H, so eps = 1 - sigma_y0/H exactly.
    expectNear("linear Considere strain", uniformElongation(linear), 1.0 - 355.0e6 / 2.0e9, 1e-15);
    expectTrue("a curve that hardens no faster than it is strong never deforms uniformly",
               uniformElongation(linearHardening(355.0e6, 300.0e6)) == 0.0);

    // Swift, fitted to a tensile test. The fit is exact at both ends and its
    // Considere point is the uniform strain that went in -- three round trips, and
    // between them they pin K, eps_0 and n.
    const double yield = 355.0e6, ultimate = 568.4e6, uniform = 0.1484;
    const FlowCurve swift = swiftFromTensile(yield, ultimate, uniform);
    expectNear("Swift fit reproduces the yield strength", flowStress(swift, 0.0), yield,
               1e-6 * yield);
    expectNear("Swift fit reproduces the ultimate true stress", flowStress(swift, uniform),
               ultimate, 1e-6 * ultimate);
    expectNear("Swift fit reproduces the uniform strain", uniformElongation(swift), uniform, 1e-9);

    // Considere, found independently: the strain where the slope crosses the
    // stress. Bisection on the curve, nothing to do with the formula above.
    double lo = 1e-6, hi = 2.0;
    for (int i = 0; i < 200; ++i) {
        const double mid = 0.5 * (lo + hi);
        if (flowSlope(swift, mid) > flowStress(swift, mid))
            lo = mid;
        else
            hi = mid;
    }
    expectNear("Considere found by bisection agrees with uniformElongation()", 0.5 * (lo + hi),
               uniformElongation(swift), 1e-9);

    // The slope is the derivative of the stress, checked rather than assumed.
    for (double p : {0.001, 0.02, 0.15, 0.4}) {
        const double h = 1e-7;
        const double fd = (flowStress(swift, p + h) - flowStress(swift, p - h)) / (2.0 * h);
        expectNear("Swift slope is the derivative at eps_p = " + std::to_string(p),
                   flowSlope(swift, p), fd, 1e-5 * std::abs(fd));
    }

    // The construction that looks natural and is wrong: taking eps_0 as the elastic
    // strain at yield puts the whole curve far too high, because a power law forced
    // through the yield point of a steel with a plateau is much too steep. Recorded
    // as a guard so nobody "simplifies" swiftFromTensile back into it.
    const FlowCurve naive = swiftHardening(yield, yield / 206.0e9, uniform);
    expectTrue("eps_0 = sigma_y/E overstates the flow stress at 0.2 by more than 10%",
               flowStress(naive, 0.2) > 1.1 * flowStress(swift, 0.2));

    // The plastic material and the structural one must not drift apart.
    const Material material = steel();
    const StructuralMaterial structural = ah36Steel();
    expectNear("shipSteel E matches ah36Steel", material.youngsModulus, structural.youngsModulus,
               1.0);
    expectNear("shipSteel nu matches ah36Steel", material.poissonRatio, structural.poissonRatio,
               1e-15);
    expectNear("shipSteel yield matches ah36Steel", flowStress(material.flow, 0.0),
               structural.yieldStrength, 1.0);
    expectNear("the failure model's uniform strain is the curve's Considere point",
               material.failure.uniformStrain, uniformElongation(material.flow), 1e-12);

    // **The pairing above covers `ah36Steel` and nothing else, and the ship has
    // two materials.** `zone::Solver` reads its elastic constants, mass and
    // critical timestep from `Patch::material` -- a `StructuralMaterial` taken
    // from `structure.materials[panel.material]` -- while its *return map* runs
    // on a separately supplied `plasticity::Material`, which every caller
    // defaults to `shipSteel()`. Two sources for one patch's steel.
    //
    // For index 0 they agree, which is what the three assertions above pin. For
    // index 1 they do not: `mildSteel()` yields at 235 MPa against `shipSteel()`'s
    // 355, so a zone on the ferry's weather deck (`scantlings.cpp` sets
    // `weatherDeck.material = 1`) would integrate 235 MPa plating that does not
    // yield until 355 -- **51% too strong**, and chosen for promotion precisely
    // because Tier 0 read it as the weaker material.
    //
    // This asserts the disagreement rather than the agreement, because the fix is
    // to thread the right plastic material through `zone::indent` and that is a
    // change to the solver's interface, not to a constant. Until then the hazard
    // is real and this is what says so out loud: if someone makes the two agree,
    // or adds a third material, this fails and points at the thread that has to
    // be pulled.
    const StructuralMaterial mild = mildSteel();
    expectTrue("the ship's second material is genuinely weaker than the plastic default",
               mild.yieldStrength < flowStress(material.flow, 0.0));
    const double ratio = flowStress(material.flow, 0.0) / mild.yieldStrength;
    std::printf("     the plastic default is %.2fx the yield of the ship's mild steel"
                " (%.0f vs %.0f MPa) -- a zone on a mild-steel panel integrates the\n"
                "     wrong one, and `zone::Solver` has no way to tell\n",
                ratio, flowStress(material.flow, 0.0) / 1e6, mild.yieldStrength / 1e6);
    expectNear("and the gap is the documented 1.51x", ratio, 355.0 / 235.0, 0.01);
}

// --- uniaxial tension against the closed form ----------------------------------

void testUniaxialTension() {
    for (int kind = 0; kind < 2; ++kind) {
        const Material material = kind == 0 ? linearSteel() : steel();
        const std::string label = kind == 0 ? "linear" : "Swift";
        const double e = material.youngsModulus, nu = material.poissonRatio;
        const double yield = flowStress(material.flow, 0.0);

        // Below yield the response is Hooke's law, and yield arrives exactly at
        // eps = sigma_y/E -- not one increment early or late. The lateral strain is
        // -nu times the axial one, so the state is exactly uniaxial and the
        // boundary is exactly E eps = sigma_y.
        for (int side = 0; side < 2; ++side) {
            const double axial = (yield / e) * (side == 0 ? 1.0 - 1e-5 : 1.0 + 1e-5);
            State state;
            double stress[kV];
            const Increment increment =
                update(material, kNeverFails, V6(axial, -nu * axial, -nu * axial), state, stress);
            if (side == 0) {
                expectTrue(label + ": just below sigma_y/E is elastic", !increment.yielded);
                expectNear(label + ": elastic axial stress is E eps", stress[0], e * axial,
                           1e-9 * yield);
                expectTrue(label + ": elastic leaves no plastic strain",
                           state.equivalentPlasticStrain == 0.0);
                expectNear(label + ": elastic transverse stress is zero",
                           maxAbs(stress + 1, 2), 0.0, 1e-8 * yield);
            } else {
                expectTrue(label + ": just above sigma_y/E has yielded", increment.yielded);
            }
        }

        // The whole plastic branch, driven from the closed form. Feeding the strain
        // that theory says produces a given plastic strain, the map must return the
        // hardening curve at that plastic strain and *nothing* transverse.
        for (double plastic : {0.001, 0.01, 0.05, 0.12}) {
            const double sigma = flowStress(material.flow, plastic);
            const double axial = sigma / e + plastic;
            const double lateral = -nu * sigma / e - 0.5 * plastic;
            State state;
            double stress[kV];
            const Increment increment =
                update(material, kNeverFails, V6(axial, lateral, lateral), state, stress);

            const std::string at = label + " at eps_p = " + std::to_string(plastic);
            expectNear(at + ": axial stress is the hardening curve", stress[0], sigma,
                       1e-9 * sigma);
            expectNear(at + ": equivalent plastic strain", state.equivalentPlasticStrain, plastic,
                       1e-11 * plastic);
            expectNear(at + ": transverse stress is machine zero", maxAbs(stress + 1, 2), 0.0,
                       1e-7 * sigma);
            expectNear(at + ": no shear appears", maxAbs(stress + 3, 3), 0.0, 1e-9 * sigma);
            // Plastic incompressibility, seen as a lateral strain: the plastic part
            // contracts by exactly half, the elastic part by nu. From a virgin
            // state, so the trace is the increment's own and there is no
            // cancellation between two large accumulated numbers.
            expectNear(at + ": plastic lateral strain is -eps_p/2", state.plasticStrain[1],
                       -0.5 * plastic, 1e-11 * plastic);
            expectNear(at + ": the plastic strain increment has zero trace",
                       state.plasticStrain[0] + state.plasticStrain[1] + state.plasticStrain[2],
                       0.0, 1e-16 * maxAbs(state.plasticStrain, kV));
            // A linear curve makes the consistency equation affine, so one Newton
            // step lands on the root exactly and the second only confirms it.
            if (kind == 0)
                expectEqual(at + ": a linear curve costs two evaluations", increment.iterations,
                            2);
            else
                expectTrue(at + ": Swift converges in a handful of evaluations",
                           increment.iterations <= 8);
        }

        // The same points reached the other way round: prescribe the axial strain
        // and let a bisection through the return map find the lateral strain that
        // nulls the transverse stress. None of the closed form is used.
        for (double axial : {0.005, 0.02, 0.08}) {
            State virgin;
            const UniaxialProbe probe = uniaxialPull(material, kNeverFails, virgin, axial);
            const double plastic = probe.state.equivalentPlasticStrain;
            const std::string at = label + " pulled to eps = " + std::to_string(axial);
            expectNear(at + ": stress is the hardening curve", probe.axialStress,
                       flowStress(material.flow, plastic), 1e-8 * probe.axialStress);
            expectNear(at + ": axial strain splits into elastic and plastic",
                       probe.axialStress / e + plastic, axial, 1e-10 * axial);
            expectNear(at + ": lateral strain is -nu sigma/E - eps_p/2", probe.lateralStrain,
                       -nu * probe.axialStress / e - 0.5 * plastic, 1e-9 * axial);
            expectTrue(at + ": the plastic strain is not a rounding artefact",
                       plastic > 0.5 * axial);
        }
    }
}

// `elasticStress` is a public entry point and was, until mutation testing said so,
// reachable by no test at all: dropping the deviatoric split from it -- returning
// K tr(eps) + 2 mu eps_ii instead of K tr(eps) + 2 mu (eps_ii - tr/3) -- passed
// everything. Tied here to the moduli, computed by a different route, and to the
// elastic branch of the return map, which must be the *same* arithmetic in the same
// order and therefore agree bit for bit.
void testElasticStressAgreesWithTheModuli() {
    const Material material = steel();
    double c[kV * kV];
    elasticModuli(material, c);
    const double strain[kV] = {3.0e-4, -1.1e-4, 7.0e-5, 2.2e-4, -1.5e-4, 9.0e-5};
    expectTrue("the elastic-stress test has a volumetric part to get wrong",
               std::abs(strain[0] + strain[1] + strain[2]) > 1e-5);

    double want[kV];
    for (int i = 0; i < kV; ++i) {
        double s = 0.0;
        for (int j = 0; j < kV; ++j) s += c[i * kV + j] * strain[j];
        want[i] = s;
    }
    double direct[kV];
    elasticStress(material, strain, direct);
    double worst = 0.0;
    for (int i = 0; i < kV; ++i) worst = std::max(worst, std::abs(direct[i] - want[i]));
    expectNear("elasticStress agrees with elasticModuli times the strain", worst, 0.0,
               1e-12 * maxAbs(want, kV));

    State virgin;
    double mapped[kV];
    const Increment increment = update(material, kNeverFails, strain, virgin, mapped);
    expectTrue("the elastic-stress test stays elastic", !increment.yielded);
    for (int i = 0; i < kV; ++i)
        expectTrue("the return map's elastic branch is elasticStress, bit for bit",
                   mapped[i] == direct[i]);
}

// Backward Euler is idempotent: after the return the stress is *on* the surface, so
// asking the same question again must produce the same stress and no further flow.
// A return map that stopped its scalar Newton early would leave the stress a little
// outside and the second call would quietly finish the job.
void testTheReturnIsIdempotent() {
    for (int kind = 0; kind < 2; ++kind) {
        const Material material = kind == 0 ? linearSteel() : steel();
        const std::string label = kind == 0 ? "linear" : "Swift";
        const double strain[kV] = {0.04, -0.011, -0.013, 0.006, -0.004, 0.009};

        State state;
        double first[kV], second[kV];
        update(material, kNeverFails, strain, state, first);
        const double plastic = state.equivalentPlasticStrain;
        expectTrue(label + ": the idempotence test yielded", plastic > 0.01);

        const Increment again = update(material, kNeverFails, strain, state, second);
        // The boolean `yielded` is deliberately *not* asserted false. The trial
        // state sits exactly on the surface, so whether the second call sees an
        // excess of +1e-8 Pa or -1e-8 Pa is floating-point equality between two
        // independently recomputed quantities, and no answer to that is more
        // correct than the other. Putting a tolerance into the yield check to make
        // the flag behave would buy a tidier boolean with a tuned parameter that no
        // measurement sets -- the same trade `solid_shell.hpp` rejects for
        // hourglass stiffness. What is substantive is that no *flow* happens, and
        // that is asserted at the rounding floor.
        expectNear(label + ": repeating the same strain adds no plastic strain",
                   state.equivalentPlasticStrain, plastic, 1e-15 * plastic);
        expectTrue(label + ": and no meaningful plastic multiplier",
                   again.plasticMultiplier <= 1e-14 * plastic);
        double worst = 0.0;
        for (int i = 0; i < kV; ++i) worst = std::max(worst, std::abs(second[i] - first[i]));
        expectNear(label + ": and returns the same stress", worst, 0.0, 1e-13 * maxAbs(first, kV));
    }
}

// --- the property radial return has and forward Euler does not ------------------

void testStepIndependence() {
    for (int kind = 0; kind < 2; ++kind) {
        const Material material = kind == 0 ? linearSteel() : steel();
        const std::string label = kind == 0 ? "linear" : "Swift";

        // A radial deviatoric path: the components stay in fixed proportion, so the
        // deviatoric direction never moves.
        const double target[kV] = {0.06, -0.018, -0.018, 0.0, 0.0, 0.0};

        double reference[kV] = {};
        double referencePlastic = 0.0;
        for (int steps : {1, 7, 100, 10000}) {
            State state;
            double stress[kV] = {};
            for (int s = 1; s <= steps; ++s) {
                const double f = static_cast<double>(s) / steps;
                double strain[kV];
                for (int i = 0; i < kV; ++i) strain[i] = f * target[i];
                update(material, kNeverFails, strain, state, stress);
            }
            if (steps == 1) {
                for (int i = 0; i < kV; ++i) reference[i] = stress[i];
                referencePlastic = state.equivalentPlasticStrain;
                expectTrue(label + ": the radial path is well into plasticity",
                           referencePlastic > 0.04);
                continue;
            }
            double difference = 0.0;
            for (int i = 0; i < kV; ++i)
                difference = std::max(difference, std::abs(stress[i] - reference[i]));
            expectNear(label + ": " + std::to_string(steps) + " steps give the one-step stress",
                       difference, 0.0, 1e-12 * maxAbs(reference, kV));
            expectNear(label + ": " + std::to_string(steps) +
                           " steps give the one-step plastic strain",
                       state.equivalentPlasticStrain, referencePlastic, 1e-12 * referencePlastic);
        }

        // **Negative control.** Backward Euler is exact for a radial path and only
        // for a radial path. A path that turns a corner in deviatoric space must
        // come out step dependent, or the test above is measuring nothing.
        const auto walk = [&](int perLeg) {
            State state;
            double stress[kV] = {};
            const double legs[2][kV] = {{0.05, -0.015, -0.015, 0, 0, 0},
                                        {0.0, 0.05, -0.05, 0.02, 0, 0}};
            double from[kV] = {};
            for (const auto& leg : legs) {
                for (int s = 1; s <= perLeg; ++s) {
                    const double f = static_cast<double>(s) / perLeg;
                    double strain[kV];
                    for (int i = 0; i < kV; ++i) strain[i] = from[i] + f * (leg[i] - from[i]);
                    update(material, kNeverFails, strain, state, stress);
                }
                for (int i = 0; i < kV; ++i) from[i] = leg[i];
            }
            return state.equivalentPlasticStrain;
        };
        const double coarse = walk(1), fine = walk(400);
        expectTrue(label + ": a non-proportional path IS step dependent",
                   std::abs(coarse - fine) > 1e-4 * fine);
    }
}

// --- the yield surface is a cylinder --------------------------------------------

void testHydrostaticStatesDoNotYield() {
    const Material material = steel();
    const double kappa = material.bulkModulus();

    // Volumetric strains far past anything the deviator could survive, both signs.
    // Nothing plastic may happen at all: not a small amount, none.
    for (double volumetric : {0.05, -0.05, 0.5, -0.5}) {
        State state;
        double stress[kV];
        const double third = volumetric / 3.0;
        const Increment increment =
            update(material, kNeverFails, V6(third, third, third), state, stress);
        const std::string at = "hydrostatic strain " + std::to_string(volumetric);
        expectTrue(at + ": does not yield", !increment.yielded);
        expectTrue(at + ": leaves the equivalent plastic strain at exactly zero",
                   state.equivalentPlasticStrain == 0.0);
        expectTrue(at + ": leaves the plastic strain tensor at exactly zero",
                   maxAbs(state.plasticStrain, kV) == 0.0);
        expectTrue(at + ": accumulates exactly no damage", state.damage == 0.0);
        for (int i = 0; i < 3; ++i)
            expectNear(at + ": normal stress is K tr(eps)", stress[i], kappa * volumetric,
                       1e-9 * std::abs(kappa * volumetric));
        expectNear(at + ": carries no shear", maxAbs(stress + 3, 3), 0.0, 0.0);
        expectNear(at + ": von Mises stress is zero", vonMises(stress), 0.0,
                   1e-9 * std::abs(kappa * volumetric));
    }

    // Stronger: an arbitrary wandering pressure superposed on a radial deviatoric
    // path must not move the deviatoric answer. This is what separates a correct
    // deviator from one that has quietly kept a volumetric term.
    const double deviatoricTarget[kV] = {0.04, -0.005, -0.035, 0.02, -0.03, 0.016};
    const auto run = [&](bool withPressure, double outStress[kV]) {
        State state;
        for (int i = 0; i < kV; ++i) outStress[i] = 0.0;
        for (int s = 1; s <= 200; ++s) {
            const double f = static_cast<double>(s) / 200.0;
            const double wobble =
                withPressure ? 0.02 * std::sin(11.0 * f) + 0.015 * std::cos(3.0 * f) : 0.0;
            double strain[kV];
            for (int i = 0; i < kV; ++i)
                strain[i] = f * deviatoricTarget[i] + (i < 3 ? wobble : 0.0);
            update(material, kNeverFails, strain, state, outStress);
        }
        return state;
    };
    double plainStress[kV], pressedStress[kV];
    const State plain = run(false, plainStress);
    const State pressed = run(true, pressedStress);

    expectTrue("the superposed pressure test is not vacuous", plain.equivalentPlasticStrain > 0.02);
    expectNear("a wandering pressure leaves the equivalent plastic strain unchanged",
               pressed.equivalentPlasticStrain, plain.equivalentPlasticStrain,
               1e-13 * plain.equivalentPlasticStrain);
    double plainDeviator[kV], pressedDeviator[kV];
    deviator(plainStress, plainDeviator);
    deviator(pressedStress, pressedDeviator);
    double worst = 0.0;
    for (int i = 0; i < kV; ++i)
        worst = std::max(worst, std::abs(plainDeviator[i] - pressedDeviator[i]));
    expectNear("a wandering pressure leaves the deviatoric stress unchanged", worst, 0.0,
               1e-11 * maxAbs(plainDeviator, kV));
    // ... and it did change the pressure, so the comparison had something to compare.
    expectTrue("the wandering pressure actually moved the mean stress",
               std::abs(meanStress(pressedStress) - meanStress(plainStress)) > 1.0e9);
}

// --- consistency, incompressibility, dissipation --------------------------------

void testConsistencyAndIncompressibility() {
    std::mt19937 rng(20260803u);
    std::uniform_real_distribution<double> noise(-3e-4, 3e-4);

    for (int kinematic = 0; kinematic < 2; ++kinematic) {
        Material material = steel();
        if (kinematic != 0) material.flow.kinematicModulus = 1.5e9;
        const std::string label = kinematic == 0 ? "isotropic" : "combined";

        // A cyclic path with noise on it: reversals are where a mislabelled back
        // stress and a sign error in the dissipation both show up, and a monotonic
        // pull would never visit them.
        const double shape[kV] = {1.0, -0.3, -0.3, 0.4, -0.2, 0.15};

        State state;
        double stress[kV] = {};
        int plasticSteps = 0;
        double worstConsistency = 0.0, worstTrace = 0.0, worstDissipation = 0.0;
        double totalDissipation = 0.0;

        for (int step = 1; step <= 400; ++step) {
            double strain[kV];
            const double phase = 2.0 * std::numbers::pi * step / 160.0;
            for (int i = 0; i < kV; ++i)
                strain[i] = 0.05 * shape[i] * std::sin(phase) + noise(rng);
            const double before = state.equivalentPlasticStrain;
            const Increment increment = update(material, kNeverFails, strain, state, stress);
            if (!increment.yielded) continue;
            ++plasticSteps;

            // After the return the stress sits **on** the surface: the norm of the
            // relative deviatoric stress is exactly the yield radius.
            double relative[kV];
            deviator(stress, relative);
            for (int i = 0; i < kV; ++i) relative[i] -= state.backStress[i];
            const double radius = std::sqrt(
                relative[0] * relative[0] + relative[1] * relative[1] + relative[2] * relative[2] +
                2.0 * (relative[3] * relative[3] + relative[4] * relative[4] +
                       relative[5] * relative[5]));
            const double want =
                std::sqrt(2.0 / 3.0) * flowStress(material.flow, state.equivalentPlasticStrain);
            worstConsistency = std::max(worstConsistency, std::abs(radius - want) / want);

            // The plastic strain stays deviatoric, exactly. Measured on the
            // *accumulated* tensor against its own magnitude, not on the increment:
            // differencing two accumulated numbers of size 0.05 to recover an
            // increment of size 1e-5 loses eleven digits to cancellation and would
            // be measuring the test's arithmetic rather than the code's. The
            // increment's own trace is checked without cancellation in
            // testUniaxialTension, from a virgin state.
            const double trace =
                state.plasticStrain[0] + state.plasticStrain[1] + state.plasticStrain[2];
            // Against the *accumulated* flow, which is what the drift is
            // proportional to: each increment contributes a trace of order
            // eps * dgamma, so the sum grows with the accumulated plastic strain
            // and not with the current tensor, which a reversal can take back
            // through zero.
            worstTrace = std::max(worstTrace,
                                  std::abs(trace) / state.equivalentPlasticStrain);

            // Dissipation, computed by contraction inside the map, against the
            // closed form sigma_y d eps_p.
            const double closed = flowStress(material.flow, state.equivalentPlasticStrain) *
                                  (state.equivalentPlasticStrain - before);
            worstDissipation = std::max(
                worstDissipation, std::abs(increment.dissipation - closed) / std::max(closed, 1.0));
            expectTrue(label + ": dissipation is never negative", increment.dissipation >= 0.0);
            totalDissipation += increment.dissipation;
        }

        std::printf("     %-9s cyclic path: %d plastic steps, accumulated eps_p = %.3f\n",
                    label.c_str(), plasticSteps, state.equivalentPlasticStrain);
        expectTrue(label + ": the cyclic path yielded, repeatedly", plasticSteps > 100);
        std::printf("     %-9s worst consistency %.2e, worst deviatoric drift %.2e\n",
                    label.c_str(), worstConsistency, worstTrace);
        expectNear(label + ": the stress sits on the yield surface", worstConsistency, 0.0, 1e-14);
        expectNear(label + ": the accumulated plastic strain stays deviatoric", worstTrace, 0.0,
                   1e-14);
        expectNear(label + ": dissipation equals sigma_y d eps_p", worstDissipation, 0.0, 1e-12);
        expectTrue(label + ": total dissipation is positive", totalDissipation > 0.0);
    }
}

void testUnloadingIsElastic() {
    const Material material = linearSteel();
    double c[kV * kV];
    elasticModuli(material, c);

    State state;
    double stress[kV];
    const double loaded[kV] = {0.03, -0.009, -0.009, 0.004, 0.0, 0.0};
    update(material, kNeverFails, loaded, state, stress);
    expectTrue("the unloading test starts from a plastic state",
               state.equivalentPlasticStrain > 0.01);

    const State frozen = state;
    double before[kV];
    for (int i = 0; i < kV; ++i) before[i] = stress[i];

    // Take the strain back. The response must be the *original* elastic modulus --
    // unloading does not see the plastic tangent -- and nothing may be added to the
    // plastic strain.
    double unloaded[kV];
    for (int i = 0; i < kV; ++i) unloaded[i] = 0.999 * loaded[i];
    const Increment increment = update(material, kNeverFails, unloaded, state, stress);
    expectTrue("unloading does not yield", !increment.yielded);
    expectTrue("unloading leaves the equivalent plastic strain untouched",
               state.equivalentPlasticStrain == frozen.equivalentPlasticStrain);
    for (int i = 0; i < kV; ++i)
        expectTrue("unloading leaves the plastic strain tensor untouched",
                   state.plasticStrain[i] == frozen.plasticStrain[i]);

    double expected[kV];
    for (int i = 0; i < kV; ++i) {
        double s = 0.0;
        for (int j = 0; j < kV; ++j) s += c[i * kV + j] * (unloaded[j] - loaded[j]);
        expected[i] = s;
    }
    double worst = 0.0;
    for (int i = 0; i < kV; ++i)
        worst = std::max(worst, std::abs((stress[i] - before[i]) - expected[i]));
    expectNear("unloading follows the original elastic modulus", worst, 0.0,
               1e-9 * maxAbs(expected, kV));
    expectTrue("the unloading step was large enough to measure",
               maxAbs(expected, kV) > 1.0e6);
}

// --- isotropy of the yield surface ----------------------------------------------

void testRotationalInvariance() {
    const Material material = steel();
    const double axis[3] = {0.37, -0.81, 0.45};
    double r[3][3];
    axisRotation(axis, 1.1, r);

    // The same strain path, in two frames. The stress must be the rotated stress
    // and the equivalent plastic strain -- an invariant -- must be identical.
    const double target[kV] = {0.05, -0.01, -0.02, 0.03, -0.015, 0.02};

    State plain, rotated;
    double plainStress[kV] = {}, rotatedStress[kV] = {};
    for (int s = 1; s <= 40; ++s) {
        const double f = static_cast<double>(s) / 40.0;
        double strain[kV], turned[kV];
        for (int i = 0; i < kV; ++i) strain[i] = f * target[i];
        rotateVoigt(r, strain, true, turned);
        update(material, kNeverFails, strain, plain, plainStress);
        update(material, kNeverFails, turned, rotated, rotatedStress);
    }
    expectTrue("the rotation test is not vacuous", plain.equivalentPlasticStrain > 0.02);

    double want[kV];
    rotateVoigt(r, plainStress, false, want);
    double worst = 0.0;
    for (int i = 0; i < kV; ++i) worst = std::max(worst, std::abs(rotatedStress[i] - want[i]));
    expectNear("a rotated strain gives the rotated stress", worst, 0.0, 1e-9 * maxAbs(want, kV));
    expectNear("a rotated strain gives the same equivalent plastic strain",
               rotated.equivalentPlasticStrain, plain.equivalentPlasticStrain,
               1e-12 * plain.equivalentPlasticStrain);
    expectTrue("the rotation actually moved the stress components",
               std::abs(rotatedStress[0] - plainStress[0]) > 0.05 * maxAbs(plainStress, kV));
}

void testPureShearYieldsAtYieldOverRootThree() {
    const Material material = linearSteel();
    const double mu = material.shearModulus();
    const double yield = flowStress(material.flow, 0.0);
    // tau_y = sigma_y / sqrt(3), so gamma_y = sigma_y / (sqrt(3) mu). Checked on all
    // three shear planes: mixing up two Voigt shear slots survives a test that only
    // ever loads one of them.
    const double gammaYield = yield / (std::sqrt(3.0) * mu);
    for (int component = 3; component < 6; ++component) {
        const std::string at = "shear component " + std::to_string(component);
        State below, above;
        double stress[kV];
        double strain[kV] = {};
        strain[component] = gammaYield * (1.0 - 1e-9);
        const Increment cold = update(material, kNeverFails, strain, below, stress);
        expectTrue(at + " is elastic below tau_y", !cold.yielded);
        expectNear(at + " reaches sigma_y at tau_y", vonMises(stress), yield, 1e-6 * yield);
        expectNear(at + " carries mu gamma", stress[component], mu * strain[component],
                   1e-9 * yield);
        expectNear(at + " carries no normal stress in pure shear", maxAbs(stress, 3), 0.0,
                   1e-9 * yield);

        strain[component] = gammaYield * (1.0 + 1e-6);
        const Increment hot = update(material, kNeverFails, strain, above, stress);
        expectTrue(at + " yields above tau_y", hot.yielded);
    }
}

// --- kinematic hardening, and the one thing that distinguishes it ---------------

void testBauschingerEffect() {
    const double hardening = 3.0e9;
    Material isotropic = linearSteel(hardening);
    Material kinematic = linearSteel(0.0);
    kinematic.flow.kinematicModulus = hardening;
    const double yield = flowStress(isotropic.flow, 0.0);

    // Forward, the two are *indistinguishable*: growing the surface and moving it
    // give the same flow stress on a monotonic path. That is the point -- and it is
    // why a suite that only pulls cannot tell them apart.
    State a, b;
    const double pull = 0.02;
    UniaxialProbe forwardIsotropic, forwardKinematic;
    for (int s = 1; s <= 40; ++s) {
        const double axial = pull * s / 40.0;
        forwardIsotropic = uniaxialPull(isotropic, kNeverFails, a, axial);
        a = forwardIsotropic.state;
        forwardKinematic = uniaxialPull(kinematic, kNeverFails, b, axial);
        b = forwardKinematic.state;
    }
    expectNear("forward, isotropic and kinematic hardening are the same curve",
               forwardKinematic.axialStress, forwardIsotropic.axialStress,
               1e-7 * forwardIsotropic.axialStress);
    const double plastic = a.equivalentPlasticStrain;
    expectTrue("the Bauschinger test reached real plastic strain", plastic > 0.015);
    expectNear("the two hardening laws accumulate the same plastic strain",
               b.equivalentPlasticStrain, plastic, 1e-9 * plastic);

    // Reversed, they part company by a closed form. The reverse yield stress is
    // found by bisecting on the axial strain at which plastic flow restarts.
    const auto reverseYieldStress = [&](const Material& material, const State& state) {
        double lo = pull - 0.02, hi = pull;
        for (int i = 0; i < 100; ++i) {
            const double mid = 0.5 * (lo + hi);
            if (uniaxialPull(material, kNeverFails, state, mid).increment.yielded)
                lo = mid;
            else
                hi = mid;
        }
        return uniaxialPull(material, kNeverFails, state, lo).axialStress;
    };
    const double isotropicReverse = reverseYieldStress(isotropic, a);
    const double kinematicReverse = reverseYieldStress(kinematic, b);
    const double isotropicSpan = forwardIsotropic.axialStress - isotropicReverse;
    const double kinematicSpan = forwardKinematic.axialStress - kinematicReverse;

    // The elastic span on reversal is 2 sigma_y(eps_p) for isotropic hardening and
    // exactly 2 sigma_y0 for kinematic, whatever the prestrain -- the whole
    // Bauschinger effect in two closed forms.
    expectNear("isotropic hardening: the elastic span on reversal is 2 sigma_y(eps_p)",
               isotropicSpan, 2.0 * flowStress(isotropic.flow, plastic), 1e-3 * yield);
    expectNear("kinematic hardening: the elastic span on reversal is exactly 2 sigma_y0",
               kinematicSpan, 2.0 * yield, 1e-3 * yield);
    expectTrue("the two spans genuinely differ",
               isotropicSpan - kinematicSpan > 0.5 * hardening * plastic);
}

// --- the algorithmic tangent ----------------------------------------------------

void testAlgorithmicTangent() {
    Material material = steel();
    material.flow.kinematicModulus = 8.0e8;

    struct Case {
        const char* label;
        double strain[kV];
    };
    const Case cases[] = {
        {"elastic", {0.0008, -0.0002, -0.0002, 0.0, 0.0, 0.0}},
        {"uniaxial plastic", {0.03, -0.009, -0.009, 0.0, 0.0, 0.0}},
        {"shear plastic", {0.0, 0.0, 0.0, 0.02, 0.01, -0.006}},
        {"general plastic", {0.02, -0.004, -0.011, 0.007, -0.009, 0.013}},
    };

    double elastic[kV * kV];
    elasticModuli(material, elastic);

    for (const Case& item : cases) {
        // Every probe restarts from the same history, so the finite difference is
        // of this map and not of a moving one.
        State origin;
        double warm[kV];
        update(material, kNeverFails, V6(0.004, -0.0012, -0.0012), origin, warm);

        State centre = origin;
        double stress[kV], tangent[kV * kV];
        update(material, kNeverFails, item.strain, centre, stress, tangent);

        const double h = 1e-7;
        const double scale = maxAbs(tangent, kV * kV);
        double worst = 0.0;
        for (int j = 0; j < kV; ++j) {
            double plus[kV], minus[kV], up[kV], down[kV];
            for (int i = 0; i < kV; ++i) plus[i] = minus[i] = item.strain[i];
            plus[j] += h;
            minus[j] -= h;
            State pa = origin, pb = origin;
            update(material, kNeverFails, plus, pa, up);
            update(material, kNeverFails, minus, pb, down);
            for (int i = 0; i < kV; ++i)
                worst =
                    std::max(worst, std::abs(tangent[i * kV + j] - (up[i] - down[i]) / (2.0 * h)));
        }
        expectNear(std::string("algorithmic tangent matches a finite difference: ") + item.label,
                   worst, 0.0, 1e-5 * scale);

        double asymmetry = 0.0;
        for (int i = 0; i < kV; ++i)
            for (int j = 0; j < kV; ++j)
                asymmetry =
                    std::max(asymmetry, std::abs(tangent[i * kV + j] - tangent[j * kV + i]));
        expectNear(std::string("the tangent is symmetric: ") + item.label, asymmetry, 0.0,
                   1e-12 * scale);
    }

    // The elastic branch returns the elastic moduli, not something close.
    State virgin;
    double stress[kV], tangent[kV * kV];
    update(material, kNeverFails, V6(1e-5), virgin, stress, tangent);
    double worst = 0.0;
    for (int i = 0; i < kV * kV; ++i) worst = std::max(worst, std::abs(tangent[i] - elastic[i]));
    expectNear("the elastic branch returns the elastic moduli exactly", worst, 0.0, 0.0);

    // A plastic tangent must be softer than the elastic one, which is the whole
    // reason for having it.
    State plastic;
    update(material, kNeverFails, V6(0.03, -0.009, -0.009), plastic, stress, tangent);
    expectTrue("the plastic tangent is softer than the elastic one",
               tangent[0] < 0.9 * elastic[0]);
}

// --- failure --------------------------------------------------------------------

void testFailureStrainRegularisation() {
    const Failure failure = steel().failure;
    const double thickness = 0.020;

    expectNear("an element one thickness across sees the local fracture strain",
               regularisedFailureStrain(failure, thickness, thickness), failure.fractureStrain,
               1e-15);
    expectNear("an element smaller than the thickness is not extrapolated past it",
               regularisedFailureStrain(failure, 0.5 * thickness, thickness),
               failure.fractureStrain, 1e-15);
    expectTrue("a very large element sees only the uniform strain",
               std::abs(regularisedFailureStrain(failure, 100.0, thickness) -
                        failure.uniformStrain) < 1e-3);

    double previous = 1e30;
    for (double length : {0.02, 0.03, 0.05, 0.1, 0.2, 0.5}) {
        const double strain = regularisedFailureStrain(failure, length, thickness);
        expectTrue("the failure strain falls with element size at l = " + std::to_string(length),
                   strain < previous);
        previous = strain;

        // The identity the regularisation exists to enforce: the *necking*
        // elongation an element must accumulate before it tears is a property of the
        // plate, not of the mesh.
        expectNear("necking elongation is mesh invariant at l = " + std::to_string(length),
                   (strain - failure.uniformStrain) * length,
                   (failure.fractureStrain - failure.uniformStrain) * thickness, 1e-15);
    }

    // The derivation, checked against itself. Take a strain field that is a uniform
    // background plus a neck of width t at the local fracture strain, average it
    // over an element of length l, and the answer must be the regularised failure
    // strain. That is what the formula *is*, and it is not a fit.
    for (double length : {0.02, 0.04, 0.08, 0.16}) {
        constexpr int kSamples = 400001;
        double sum = 0.0;
        for (int i = 0; i < kSamples; ++i) {
            const double x = length * ((i + 0.5) / kSamples - 0.5);  // centred on the neck
            sum += std::abs(x) <= 0.5 * thickness ? failure.fractureStrain : failure.uniformStrain;
        }
        expectNear("a neck of width t averaged over an element of length " +
                       std::to_string(length),
                   sum / kSamples, regularisedFailureStrain(failure, length, thickness), 1e-5);
        // ... and a *constant* failure strain does not reproduce that average, so
        // the check above is discriminating rather than trivially satisfied.
        if (length > 0.03)
            expectTrue("a mesh-independent failure strain would not match that average",
                       std::abs(sum / kSamples - failure.fractureStrain) > 0.1);
    }

    // Triaxiality: 1 at the reference, exp(-1/2) in equibiaxial tension, and no
    // damage at all below the cutoff -- which also bounds the multiplier by e.
    expectNear("the triaxiality factor is 1 where the constants were measured",
               triaxialityFactor(failure, failure.referenceTriaxiality), 1.0, 1e-15);
    expectNear("equibiaxial tension fails at exp(-1/2) of uniaxial",
               triaxialityFactor(failure, 2.0 / 3.0), std::exp(-0.5), 1e-15);
    expectTrue("below the cutoff no damage accumulates at all",
               std::isinf(triaxialityFactor(failure, failure.cutoffTriaxiality - 1e-9)));
    expectTrue("the multiplier is bounded by e",
               triaxialityFactor(failure, failure.cutoffTriaxiality + 1e-12) <=
                   std::exp(1.0) + 1e-9);
}

void testFailureIsReachedAndIsIrreversible() {
    const Material material = steel();
    const double thickness = 0.020, length = 0.050;
    const double critical = regularisedFailureStrain(material.failure, length, thickness);

    // Uniaxial tension sits exactly at the reference triaxiality, so the damage is
    // eps_p / eps_f and the threshold is an equality.
    const auto pullTo = [&](double plastic, State& state) {
        const double sigma = flowStress(material.flow, plastic);
        const double axial = sigma / material.youngsModulus + plastic;
        const double lateral =
            -material.poissonRatio * sigma / material.youngsModulus - 0.5 * plastic;
        double stress[kV];
        return update(material, critical, V6(axial, lateral, lateral), state, stress);
    };

    State under;
    pullTo(0.995 * critical, under);
    expectTrue("pulled to just under the failure strain, it has not failed", !under.failed);
    expectNear("damage is the plastic strain over the failure strain", under.damage, 0.995, 1e-9);

    State over;
    const Increment tearing = pullTo(1.005 * critical, over);
    expectTrue("pulled past the failure strain, it has failed", over.failed);
    expectTrue("and it says so in the increment", tearing.failedNow);
    expectNear("the damage saturates at one", over.damage, 1.0, 0.0);
    // The tear opens on the plane normal to the maximum principal stress, which for
    // a uniaxial pull is the pull axis.
    expectNear("the failure plane normal is the pull axis", std::abs(tearing.failureNormal[0]),
               1.0, 1e-12);
    expectNear("the failure plane normal has no transverse component",
               maxAbs(tearing.failureNormal + 1, 2), 0.0, 1e-12);

    // The same tear, in a **rotated** frame. Until mutation testing said so, the
    // failure normal was only ever asked for on an axis-aligned pull, where the
    // stress tensor is already diagonal, the Jacobi sweep has nothing to do and the
    // eigenvectors come back as the identity however badly they are accumulated.
    // Freezing them entirely passed.
    {
        const double axis[3] = {0.31, 0.62, -0.72};
        double r[3][3];
        axisRotation(axis, 0.9, r);
        const double plastic = 1.005 * critical;
        const double sigma = flowStress(material.flow, plastic);
        const double axial = sigma / material.youngsModulus + plastic;
        const double lateral =
            -material.poissonRatio * sigma / material.youngsModulus - 0.5 * plastic;
        const double aligned[kV] = {axial, lateral, lateral, 0, 0, 0};
        double turned[kV];
        rotateVoigt(r, aligned, true, turned);
        expectTrue("the rotated tear is genuinely off-axis", maxAbs(turned + 3, 3) > 0.05 * axial);

        State state;
        double stress[kV];
        const Increment increment = update(material, critical, turned, state, stress);
        expectTrue("a rotated pull tears too", increment.failedNow);
        // The pull axis, rotated. Compared through the absolute dot product, since
        // an eigenvector's sign carries no information.
        double dot = 0.0;
        for (int i = 0; i < 3; ++i) dot += increment.failureNormal[i] * r[i][0];
        expectNear("the failure plane normal follows the rotation", std::abs(dot), 1.0, 1e-9);
    }

    // Irreversible. Unloading, reloading and pushing the other way all leave it
    // failed and stressless.
    State failed = over;
    const State frozen = failed;
    for (double axial : {0.0, -0.05, 0.5, 0.001}) {
        double stress[kV];
        const Increment increment = update(material, critical, V6(axial), failed, stress);
        expectTrue("a failed point stays failed at eps = " + std::to_string(axial), failed.failed);
        expectNear("a failed point carries no stress at eps = " + std::to_string(axial),
                   maxAbs(stress, kV), 0.0, 0.0);
        expectTrue("a failed point accumulates no more plastic strain", !increment.yielded);
        expectTrue("a failed point's history is frozen",
                   failed.equivalentPlasticStrain == frozen.equivalentPlasticStrain &&
                       failed.damage == frozen.damage);
    }

    // Compression does not open voids, so below the cutoff triaxiality no damage
    // accumulates at all -- **exactly** none, not a little. Driven here with a
    // superposed pressure so the state sits clearly below the cutoff rather than
    // on it: uniaxial compression is at eta = -1/3, which *is* the cutoff, and
    // whether it lands a rounding error above or below the boundary is not
    // something to build an assertion on. What holds robustly for uniaxial
    // compression is the bounded rate below.
    {
        State state;
        double stress[kV];
        double strain[kV] = {};
        for (int s = 1; s <= 200; ++s) {
            const double f = 0.0005 * s;
            strain[0] = -2.0 * f;  // deviatoric compression ...
            strain[1] = strain[2] = f;
            strain[0] -= 0.004;  // ... plus a hydrostatic squeeze, eta well below -1/3
            strain[1] -= 0.004;
            strain[2] -= 0.004;
            update(material, critical, strain, state, stress);
        }
        expectTrue("triaxial compression really did flow", state.equivalentPlasticStrain > 0.05);
        expectTrue("and it is genuinely below the cutoff", triaxiality(stress) < -0.4);
        expectTrue("below the cutoff, exactly no damage accumulates", state.damage == 0.0);
        expectTrue("and therefore it never fails", !state.failed);
    }

    // And at the cutoff itself the multiplier is bounded by e, so uniaxial
    // compression damages at most 1/e as fast as uniaxial tension at the same
    // plastic strain -- a closed form, and the statement that survives the
    // knife edge.
    {
        const auto damageAfter = [&](double sign) {
            State state;
            double stress[kV];
            for (int s = 1; s <= 100; ++s) {
                const double plastic = 0.001 * s;
                const double sigma = sign * flowStress(material.flow, plastic);
                const double axial = sigma / material.youngsModulus + sign * plastic;
                const double lateral = -material.poissonRatio * sigma / material.youngsModulus -
                                       0.5 * sign * plastic;
                update(material, critical, V6(axial, lateral, lateral), state, stress);
            }
            return state.damage;
        };
        const double tension = damageAfter(+1.0), compression = damageAfter(-1.0);
        expectTrue("the tension leg accumulated damage", tension > 0.1);
        expectTrue("compression damages at most 1/e as fast as tension",
                   compression <= tension / std::exp(1.0) + 1e-12);
    }

    // Damage is *accumulated*, not compared: a path that spends part of its life at
    // a high triaxiality must spend the corresponding share of its life. Summed here
    // independently and compared with what the state carries.
    {
        State state;
        double stress[kV];
        double expected = 0.0;
        double strain[kV] = {};
        double lowest = 1e30, highest = -1e30;
        for (int s = 1; s <= 300; ++s) {
            // A deviatoric ramp with a growing hydrostatic tension on it, so the
            // triaxiality sweeps a real range instead of sitting still. The
            // volumetric part has to be small: at 3% it would put the mean stress
            // at 5 GPa against a 0.6 GPa flow stress, eta above 8, and the failure
            // strain at 4e-6 -- which is what the first version of this test did,
            // and it tore at a plastic strain of 1e-5.
            const double f = static_cast<double>(s) / 300.0;
            const double volumetric = 0.0016 * f;
            strain[0] = 0.05 * f + volumetric;
            strain[1] = strain[2] = -0.025 * f + volumetric;
            const double before = state.equivalentPlasticStrain;
            const Increment increment = update(material, critical, strain, state, stress);
            if (increment.yielded) {
                lowest = std::min(lowest, increment.triaxiality);
                highest = std::max(highest, increment.triaxiality);
                expected += (state.equivalentPlasticStrain - before) /
                            (critical * triaxialityFactor(material.failure, increment.triaxiality));
            }
        }
        expectTrue("the accumulating-damage path did yield", state.equivalentPlasticStrain > 0.02);
        expectTrue("and did not tear, so the damage was not clamped", !state.failed);
        expectTrue("and the triaxiality really moved along it", highest - lowest > 0.15);
        expectNear("damage is the path integral of d eps_p / eps_f(eta)", state.damage, expected,
                   1e-9 * std::max(expected, 1e-9));
    }
}

// --- what the regularisation can and cannot deliver -----------------------------

void testLocalisationAndItsLimit() {
    // Points in series carrying the same load, the first 5% thinner -- the
    // imperfection that starts a neck. Under load control the thin one carries more
    // stress, yields further, and reaches its failure strain first. That much the
    // criterion does deliver, and it is what picks the site of a tear.
    const Material material = steel();
    const double thickness = 0.020;

    struct Chain {
        double elongation = 0.0;
        double peakPlastic = 0.0;
        int firstToFail = -1;
    };
    const auto pull = [&](int count, double totalLength) {
        const double length = totalLength / count;
        const double critical = regularisedFailureStrain(material.failure, length, thickness);
        std::vector<State> points(static_cast<std::size_t>(count));
        std::vector<double> strainOf(static_cast<std::size_t>(count), 0.0);
        Chain chain;
        const double peak = 1.5 * flowStress(material.flow, material.failure.fractureStrain);
        for (int step = 1; step <= 300 && chain.firstToFail < 0; ++step) {
            const double load = peak * step / 300.0;
            for (int p = 0; p < count; ++p) {
                const std::size_t index = static_cast<std::size_t>(p);
                const double target = p == 0 ? load / 0.95 : load;  // the thinned point

                // Bisect for the strain that carries the load, counting **failure
                // as overshoot**. Without that the search is chasing a
                // non-monotone function: past the failure strain the stress drops
                // to zero, the "stress is still too low" branch keeps pushing, and
                // the bracket runs away to the end of its range. The first version
                // of this test did exactly that and reported a bar tearing at 93%
                // strain.
                double lo = strainOf[index], hi = 2.0;
                for (int i = 0; i < 45; ++i) {
                    const double mid = 0.5 * (lo + hi);
                    const UniaxialProbe probe =
                        uniaxialPull(material, critical, points[index], mid, 50);
                    if (probe.state.failed || probe.axialStress >= target)
                        hi = mid;
                    else
                        lo = mid;
                }
                const UniaxialProbe probe =
                    uniaxialPull(material, critical, points[index], hi, 50);
                points[index] = probe.state;
                strainOf[index] = hi;
                if (probe.state.failed && chain.firstToFail < 0) chain.firstToFail = p;
            }
            chain.elongation = 0.0;
            chain.peakPlastic = 0.0;
            for (int p = 0; p < count; ++p) {
                const std::size_t index = static_cast<std::size_t>(p);
                chain.elongation += strainOf[index] * length;
                chain.peakPlastic =
                    std::max(chain.peakPlastic, points[index].equivalentPlasticStrain);
            }
        }
        return chain;
    };

    const Chain coarse = pull(1, 0.2);
    const Chain fine = pull(4, 0.2);
    expectEqual("the thinned point is the one that tears, 1-element mesh", coarse.firstToFail, 0);
    expectEqual("the thinned point is the one that tears, 4-element mesh", fine.firstToFail, 0);
    // Each mesh tears at *its own* failure strain, which is the whole point: the
    // number the element compares against came from the element's size.
    expectNear("the 200 mm element tore at its own failure strain", coarse.peakPlastic,
               regularisedFailureStrain(material.failure, 0.2, thickness), 0.005);
    expectNear("the 50 mm element tore at its own failure strain", fine.peakPlastic,
               regularisedFailureStrain(material.failure, 0.05, thickness), 0.005);
    expectTrue("and those are genuinely different numbers",
               fine.peakPlastic - coarse.peakPlastic > 0.15);

    // And the limit, measured rather than asserted away. With hardening and no
    // softening the plastic strain in a bar does not localise, so the elongation of
    // the whole bar before the first tear stays mesh dependent. The criterion
    // regularises the strain one element must reach to tear -- which is what governs
    // a structural tear, where a crack runs element by element and each deletion is
    // itself the softening event. It does not turn a uniformly strained gauge length
    // into a mesh-independent one, and nothing without a softening mechanism could.
    std::printf("     elongation of a 200 mm bar before the first tear\n");
    std::printf("       1 element  of 200 mm: %6.2f mm\n", 1000.0 * coarse.elongation);
    std::printf("       4 elements of  50 mm: %6.2f mm   (ratio %.2f)\n", 1000.0 * fine.elongation,
                fine.elongation / coarse.elongation);
    expectTrue("both meshes tore at all", coarse.elongation > 0.0 && fine.elongation > 0.0);
}

// --- the element hook -----------------------------------------------------------

// A single solid-shell element, centred on the origin so that a bending field
// leaves the polar rotation at exactly the identity.
void plateElement(double lx, double ly, double thickness, double nodes[solidshell::kDof]) {
    const double corner[4][2] = {{-0.5 * lx, -0.5 * ly},
                                 {+0.5 * lx, -0.5 * ly},
                                 {+0.5 * lx, +0.5 * ly},
                                 {-0.5 * lx, +0.5 * ly}};
    for (int a = 0; a < solidshell::kNodes; ++a) {
        nodes[a * 3 + 0] = corner[a % 4][0];
        nodes[a * 3 + 1] = corner[a % 4][1];
        nodes[a * 3 + 2] = (a < 4 ? -0.5 : 0.5) * thickness;
    }
}

void testElementSize() {
    double nodes[solidshell::kDof];
    plateElement(0.05, 0.05, 0.02, nodes);
    double inPlane = 0.0, thickness = 0.0;
    solidshell::elementSize(nodes, &inPlane, &thickness);
    expectNear("a square element reports its in-plane size", inPlane, 0.05, 1e-15);
    expectNear("and its thickness", thickness, 0.02, 1e-15);

    plateElement(0.06, 0.04, 0.012, nodes);
    solidshell::elementSize(nodes, &inPlane, &thickness);
    expectNear("a rectangular element reports sqrt(area)", inPlane, std::sqrt(0.06 * 0.04), 1e-15);
    expectNear("and its thickness", thickness, 0.012, 1e-15);

    // Shear the top face sideways. The perpendicular thickness has not changed; the
    // slanted edge is longer, and reporting *that* would overstate how much neck the
    // element averages over.
    plateElement(0.05, 0.05, 0.02, nodes);
    for (int a = 4; a < 8; ++a) nodes[a * 3 + 0] += 0.02;
    solidshell::elementSize(nodes, &inPlane, &thickness);
    expectNear("a sheared element reports the perpendicular thickness", thickness, 0.02, 1e-14);
    expectTrue("which is shorter than its slanted edge",
               thickness < std::sqrt(0.02 * 0.02 + 0.02 * 0.02));

    // A **tapered** element, where the mid-surface is not either face. Every
    // element under test was prismatic until mutation testing pointed out that
    // taking the size off the bottom face alone passed the lot -- on a prism the
    // bottom face and the mid-surface have the same area, so the whole averaging
    // step was unexercised. Sides 40 mm at the bottom and 60 mm at the top: the
    // mid-surface is 50 mm, and the volume is the exact integral of
    // (a + (b-a) s)^2, which 2x2x2 Gauss reproduces exactly because that is a
    // quadratic.
    plateElement(0.04, 0.04, 0.02, nodes);
    for (int a = 4; a < 8; ++a)
        for (int i = 0; i < 2; ++i) nodes[a * 3 + i] *= 1.5;  // 40 mm -> 60 mm
    solidshell::elementSize(nodes, &inPlane, &thickness);
    const double frustum = 0.02 * (0.04 * 0.04 + 0.04 * 0.06 + 0.06 * 0.06) / 3.0;
    expectNear("a tapered element reports its mid-surface size", inPlane, 0.05, 1e-15);
    expectNear("and volume over mid-surface area as its thickness", thickness, frustum / 0.0025,
               1e-15);
    expectTrue("which is not what either face alone would say",
               std::abs(thickness - 0.02) > 2e-4);

    // The failure strain follows the element, not the material.
    const Material material = steel();
    solidshell::ElementPlasticState state;
    plateElement(0.05, 0.05, 0.02, nodes);
    solidshell::initialisePlasticState(nodes, material, state);
    expectNear("the element resolves its own failure strain", state.failureStrain,
               regularisedFailureStrain(material.failure, 0.05, 0.02), 1e-15);
    solidshell::ElementPlasticState coarse;
    plateElement(0.20, 0.20, 0.02, nodes);
    solidshell::initialisePlasticState(nodes, material, coarse);
    expectTrue("a coarser element of the same plate tears at a lower strain",
               coarse.failureStrain < state.failureStrain);

    std::printf("     failure strain vs element size, 20 mm plate\n");
    for (double length : {0.02, 0.05, 0.10, 0.20, 0.40}) {
        const double base = regularisedFailureStrain(material.failure, length, 0.02);
        std::printf("       l = %5.0f mm  l/t = %4.1f   uniaxial %.3f   biaxial %.3f\n",
                    1000.0 * length, length / 0.02, base,
                    base * triaxialityFactor(material.failure, 2.0 / 3.0));
    }
}

void testElementReducesToTheElasticElement() {
    // A material that cannot yield must reproduce, to rounding, what the validated
    // elastic element already computes. That is the tie between this path and the
    // one with 207 assertions behind it.
    Material material = steel();
    material.flow = linearHardening(1.0e13, 0.0);
    const StructuralMaterial structural = ah36Steel();

    double nodes[solidshell::kDof];
    plateElement(0.05, 0.04, 0.02, nodes);
    // Distort it: a patch of rectangles hides errors a real mesh does not.
    nodes[3] += 0.004;
    nodes[7 * 3 + 1] -= 0.003;
    nodes[5 * 3 + 2] += 0.001;

    double current[solidshell::kDof];
    std::mt19937 rng(4242u);
    std::uniform_real_distribution<double> jiggle(-2e-4, 2e-4);
    for (int i = 0; i < solidshell::kDof; ++i) current[i] = nodes[i] + jiggle(rng);

    for (solidshell::Formulation form :
         {solidshell::Formulation::Displacement, solidshell::Formulation::AssumedNaturalStrain,
          solidshell::Formulation::SolidShell}) {
        double stiffness[solidshell::kDof * solidshell::kDof];
        solidshell::elementStiffness(nodes, structural, form, stiffness);
        double elastic[solidshell::kDof];
        solidshell::internalForce(stiffness, nodes, current, elastic);

        solidshell::ElementPlasticState state;
        solidshell::initialisePlasticState(nodes, material, state);
        double plasticForce[solidshell::kDof];
        const solidshell::PlasticUpdate result =
            solidshell::elementPlasticUpdate(nodes, current, material, form, state, plasticForce);

        const std::string label = solidshell::name(form);
        expectTrue(label + ": the enhanced-strain solve converged", result.converged);
        expectEqual(label + ": nothing yielded", result.yieldedPoints, 0);
        double worst = 0.0;
        for (int i = 0; i < solidshell::kDof; ++i)
            worst = std::max(worst, std::abs(plasticForce[i] - elastic[i]));
        expectTrue(label + ": the elastic force is not trivially zero",
                   maxAbs(elastic, solidshell::kDof) > 1.0);
        expectNear(label + ": the plastic path reproduces the elastic internal force", worst, 0.0,
                   1e-9 * maxAbs(elastic, solidshell::kDof));
    }
}

void testElementUniformDeformationIsAPatchTest() {
    // A uniform deformation gradient produces a uniform stress, and every enhanced
    // mode is odd in a natural coordinate, so their residual is identically zero and
    // the seven parameters stay at exactly zero -- in the plastic range as much as
    // the elastic one. That is the plastic patch test.
    const Material material = steel();
    double nodes[solidshell::kDof];
    plateElement(0.05, 0.05, 0.02, nodes);
    // Distorted **in plane and prismatically**: the same offset on a node and the
    // one above it, so the top and bottom faces stay parallel. That distinction is
    // not cosmetic. `docs/07-fem-spike-findings.md` section 6 limit 2 records that
    // the ANS interpolation is exact only for a prismatic element and that the
    // *warped* patch test fails in proportion to the warp -- so moving one face's
    // node alone, which the first version of this test did, measures that known
    // limit and not the constitutive update. It came out at 1.6% and looked like a
    // plasticity bug.
    const double shift[4][2] = {{0.006, 0.0}, {0.0, 0.003}, {-0.004, 0.002}, {0.0, -0.005}};
    for (int a = 0; a < 4; ++a)
        for (int i = 0; i < 2; ++i) {
            nodes[a * 3 + i] += shift[a][i];
            nodes[(a + 4) * 3 + i] += shift[a][i];
        }

    const double gradient[3] = {0.03, -0.008, -0.012};
    double current[solidshell::kDof];
    for (int a = 0; a < solidshell::kNodes; ++a)
        for (int i = 0; i < 3; ++i)
            current[a * 3 + i] = nodes[a * 3 + i] * (1.0 + gradient[i]);

    solidshell::ElementPlasticState state;
    solidshell::initialisePlasticState(nodes, material, state);
    double force[solidshell::kDof], stress[solidshell::kGauss * 6];
    const solidshell::PlasticUpdate result = solidshell::elementPlasticUpdate(
        nodes, current, material, solidshell::Formulation::SolidShell, state, force, stress);

    expectTrue("uniform deformation: converged", result.converged);
    expectEqual("uniform deformation: every point yielded", result.yieldedPoints,
                solidshell::kGauss);
    expectNear("uniform deformation: the enhanced parameters stay at exactly zero",
               maxAbs(state.enhanced, solidshell::kEas), 0.0, 0.0);

    // The point law, asked the same question directly.
    State reference;
    double want[kV];
    const double strain[kV] = {gradient[0], gradient[1], gradient[2], 0, 0, 0};
    update(material, state.failureStrain, strain, reference, want);
    double worst = 0.0, worstPlastic = 0.0;
    for (int gp = 0; gp < solidshell::kGauss; ++gp) {
        for (int i = 0; i < kV; ++i)
            worst = std::max(worst, std::abs(stress[gp * 6 + i] - want[i]));
        worstPlastic =
            std::max(worstPlastic, std::abs(state.point[gp].equivalentPlasticStrain -
                                            reference.equivalentPlasticStrain));
    }
    expectNear("uniform deformation: every Gauss point carries the point law's stress", worst, 0.0,
               1e-8 * maxAbs(want, kV));
    expectNear("uniform deformation: and its plastic strain", worstPlastic, 0.0, 1e-12);

    // Equilibrium: a self-equilibrated stress state applies no net force.
    double sum[3] = {0, 0, 0};
    for (int a = 0; a < solidshell::kNodes; ++a)
        for (int i = 0; i < 3; ++i) sum[i] += force[a * 3 + i];
    expectNear("uniform deformation: the nodal forces sum to zero", maxAbs(sum, 3), 0.0,
               1e-9 * maxAbs(force, solidshell::kDof));
}

void testElementBendingIntoYield() {
    // Bending is where the enhanced parameters have work to do: the through-thickness
    // modes are the only reason sigma_zz relaxes, and once the material yields they
    // are no longer given by a closed form.
    const Material material = steel();
    double nodes[solidshell::kDof];
    plateElement(0.05, 0.05, 0.02, nodes);

    solidshell::ElementPlasticState state;
    solidshell::initialisePlasticState(nodes, material, state);

    const auto bend = [&](double k, double out[solidshell::kDof]) {
        for (int a = 0; a < solidshell::kNodes; ++a) {
            const double x = nodes[a * 3 + 0], z = nodes[a * 3 + 2];
            out[a * 3 + 0] = x + k * x * z;
            out[a * 3 + 1] = nodes[a * 3 + 1];
            out[a * 3 + 2] = z - 0.5 * k * x * x;
        }
    };

    double current[solidshell::kDof];
    solidshell::PlasticUpdate result;
    double force[solidshell::kDof], stress[solidshell::kGauss * 6];
    const double curvature = 3.0;  // 1/m; surface strain is curvature * t/2 = 3%
    for (int step = 1; step <= 20; ++step) {
        bend(curvature * step / 20.0, current);
        result = solidshell::elementPlasticUpdate(nodes, current, material,
                                                  solidshell::Formulation::SolidShell, state,
                                                  force, stress);
    }

    expectTrue("bending: the enhanced-strain Newton converged", result.converged);
    expectTrue("bending: it took more than one iteration, so it did something",
               result.iterations > 1);
    expectTrue("bending: the element yielded", result.yieldedPoints > 0);
    expectTrue("bending: the enhanced parameters are genuinely non-zero",
               maxAbs(state.enhanced, solidshell::kEas) > 1e-6);

    // r = int G^T sigma dV = 0 is the *definition* of the enhanced parameters, so
    // the converged state has to satisfy it -- but ||r|| is the wrong yardstick,
    // because the enhanced thickness modes are scaled by 1/t^2 and ||r|| inherits
    // that. The scale-free statement is that the Newton correction does no work
    // against the element's yield energy, and that is asserted at 1e-16.
    const double volume = 0.05 * 0.05 * 0.02;
    const double yieldEnergy = flowStress(material.flow, 0.0) * volume;
    std::printf("     bending into yield: ||r|| = %.2e sigma_y V, correction work = %.2e\n",
                result.enhancedResidual / yieldEnergy, result.enhancedWork / yieldEnergy);
    expectTrue("bending: the enhanced-strain correction does no work",
               result.enhancedWork <= 1e-16 * yieldEnergy);

    // sigma_zz still relaxes. It cannot go to machine zero as it does in the elastic
    // element: through a partly yielded section eps_zz has to be *nonlinear* in z --
    // plastic incompressibility ties it to the plastic part of eps_xx, which is zero
    // in the elastic core -- and the three enhanced thickness modes can only carry a
    // linear variation. What is left is a genuine limit of the element, measured
    // here rather than tolerated silently.
    double worstNormal = 0.0, largest = 0.0;
    for (int gp = 0; gp < solidshell::kGauss; ++gp) {
        worstNormal = std::max(worstNormal, std::abs(stress[gp * 6 + 2]));
        largest = std::max(largest, std::abs(stress[gp * 6 + 0]));
    }
    std::printf("     bending into yield: max |sigma_zz| / max |sigma_xx| = %.4f\n",
                worstNormal / largest);
    expectTrue("bending: sigma_zz stays small against sigma_xx in the plastic range",
               worstNormal < 0.1 * largest);

    // Equilibrium of the internal force, in force and in moment. The moment is the
    // one that needs the stress to be symmetric.
    double sum[3] = {0, 0, 0}, moment[3] = {0, 0, 0};
    for (int a = 0; a < solidshell::kNodes; ++a) {
        for (int i = 0; i < 3; ++i) sum[i] += force[a * 3 + i];
        const double* x = &current[a * 3];
        const double* f = &force[a * 3];
        moment[0] += x[1] * f[2] - x[2] * f[1];
        moment[1] += x[2] * f[0] - x[0] * f[2];
        moment[2] += x[0] * f[1] - x[1] * f[0];
    }
    const double scale = maxAbs(force, solidshell::kDof);
    expectNear("bending: the nodal forces sum to zero", maxAbs(sum, 3), 0.0, 1e-9 * scale);
    expectNear("bending: and carry no net moment", maxAbs(moment, 3), 0.0, 1e-10 * scale);

    // The warm start is a warm start and nothing more: from zero the Newton reaches
    // the same answer, at the cost of iterations.
    solidshell::ElementPlasticState cold = state, warm = state;
    for (int k = 0; k < solidshell::kEas; ++k) cold.enhanced[k] = 0.0;
    double coldForce[solidshell::kDof], warmForce[solidshell::kDof];
    bend(1.05 * curvature, current);
    const solidshell::PlasticUpdate coldResult = solidshell::elementPlasticUpdate(
        nodes, current, material, solidshell::Formulation::SolidShell, cold, coldForce);
    const solidshell::PlasticUpdate warmResult = solidshell::elementPlasticUpdate(
        nodes, current, material, solidshell::Formulation::SolidShell, warm, warmForce);
    double worst = 0.0;
    for (int i = 0; i < solidshell::kDof; ++i)
        worst = std::max(worst, std::abs(coldForce[i] - warmForce[i]));
    expectNear("a cold start reaches the same force as a warm one", worst, 0.0,
               1e-9 * maxAbs(warmForce, solidshell::kDof));
    expectTrue("and the warm start is worth having", warmResult.iterations <= coldResult.iterations);
    std::printf("     enhanced-strain Newton: %d iterations warm, %d cold\n",
                warmResult.iterations, coldResult.iterations);
}

void testElementFrameIndifference() {
    // The plastic history lives in the co-rotated frame, so rotating a deformed
    // element must leave every scalar of the state alone and rotate the force.
    const Material material = steel();
    double nodes[solidshell::kDof];
    plateElement(0.05, 0.05, 0.02, nodes);

    double stretched[solidshell::kDof];
    for (int a = 0; a < solidshell::kNodes; ++a) {
        stretched[a * 3 + 0] = nodes[a * 3 + 0] * 1.03;
        stretched[a * 3 + 1] = nodes[a * 3 + 1] * 0.99;
        stretched[a * 3 + 2] = nodes[a * 3 + 2] * 0.985;
    }

    solidshell::ElementPlasticState plain;
    solidshell::initialisePlasticState(nodes, material, plain);
    double plainForce[solidshell::kDof];
    solidshell::elementPlasticUpdate(nodes, stretched, material,
                                     solidshell::Formulation::SolidShell, plain, plainForce);
    expectTrue("the frame-indifference element actually yielded",
               plain.point[0].equivalentPlasticStrain > 0.005);

    const double axis[3] = {0.2, 0.5, -0.84};
    double r[3][3];
    axisRotation(axis, 1.7, r);
    double turned[solidshell::kDof];
    for (int a = 0; a < solidshell::kNodes; ++a)
        for (int i = 0; i < 3; ++i) {
            double s = 0.0;
            for (int k = 0; k < 3; ++k) s += r[i][k] * stretched[a * 3 + k];
            turned[a * 3 + i] = s;
        }

    solidshell::ElementPlasticState rotated;
    solidshell::initialisePlasticState(nodes, material, rotated);
    double rotatedForce[solidshell::kDof];
    solidshell::elementPlasticUpdate(nodes, turned, material, solidshell::Formulation::SolidShell,
                                     rotated, rotatedForce);

    double worstState = 0.0;
    for (int gp = 0; gp < solidshell::kGauss; ++gp)
        worstState = std::max(worstState, std::abs(rotated.point[gp].equivalentPlasticStrain -
                                                   plain.point[gp].equivalentPlasticStrain));
    expectNear("a rotated element accumulates the same plastic strain", worstState, 0.0,
               1e-12 * plain.point[0].equivalentPlasticStrain);

    double worstForce = 0.0;
    for (int a = 0; a < solidshell::kNodes; ++a)
        for (int i = 0; i < 3; ++i) {
            double s = 0.0;
            for (int k = 0; k < 3; ++k) s += r[i][k] * plainForce[a * 3 + k];
            worstForce = std::max(worstForce, std::abs(rotatedForce[a * 3 + i] - s));
        }
    expectNear("and its internal force is the rotated force", worstForce, 0.0,
               1e-9 * maxAbs(plainForce, solidshell::kDof));
}

void testElementTears() {
    const Material material = steel();
    double nodes[solidshell::kDof];
    plateElement(0.05, 0.05, 0.02, nodes);

    solidshell::ElementPlasticState state;
    solidshell::initialisePlasticState(nodes, material, state);

    double current[solidshell::kDof], force[solidshell::kDof];
    solidshell::PlasticUpdate result;
    double peakForce = 0.0;
    int tornAtStep = -1;
    for (int step = 1; step <= 400; ++step) {
        const double stretch = 0.003 * step;
        for (int a = 0; a < solidshell::kNodes; ++a) {
            current[a * 3 + 0] = nodes[a * 3 + 0] * (1.0 + stretch);
            current[a * 3 + 1] = nodes[a * 3 + 1] * (1.0 - 0.5 * stretch);
            current[a * 3 + 2] = nodes[a * 3 + 2] * (1.0 - 0.5 * stretch);
        }
        result = solidshell::elementPlasticUpdate(nodes, current, material,
                                                  solidshell::Formulation::SolidShell, state,
                                                  force);
        peakForce = std::max(peakForce, maxAbs(force, solidshell::kDof));
        if (state.torn && tornAtStep < 0) tornAtStep = step;
    }
    expectTrue("stretched far enough, the element tears", state.torn);
    expectTrue("and it carried real force before it did", peakForce > 1.0e5);
    expectNear("a torn element carries no force at all", maxAbs(force, solidshell::kDof), 0.0, 0.0);
    expectEqual("every integration point failed", result.failedPoints, solidshell::kGauss);
    expectTrue("the tear happened after real plastic flow", tornAtStep > 20);

    // Irreversible at the element level too.
    for (int i = 0; i < solidshell::kDof; ++i) current[i] = nodes[i];
    solidshell::elementPlasticUpdate(nodes, current, material, solidshell::Formulation::SolidShell,
                                     state, force);
    expectTrue("unloading a torn element does not heal it", state.torn);
    expectNear("and it still carries nothing", maxAbs(force, solidshell::kDof), 0.0, 0.0);

    std::printf("     a 50 mm element of 20 mm plate tore at step %d, eps_f = %.3f\n", tornAtStep,
                state.failureStrain);

    // Re-initialising must *clear* the history. Not clearing it survived the whole
    // suite, because every caller here happened to hand over a freshly constructed
    // state -- which is exactly the condition under which a re-promoted zone would
    // inherit the damage of the last collision.
    solidshell::initialisePlasticState(nodes, material, state);
    expectTrue("re-initialising clears the torn flag", !state.torn);
    expectTrue("re-initialising clears the enhanced parameters",
               maxAbs(state.enhanced, solidshell::kEas) == 0.0);
    for (int gp = 0; gp < solidshell::kGauss; ++gp)
        expectTrue("re-initialising clears every integration point",
                   !state.point[gp].failed && state.point[gp].damage == 0.0 &&
                       state.point[gp].equivalentPlasticStrain == 0.0 &&
                       maxAbs(state.point[gp].plasticStrain, kV) == 0.0);
}

void testPartialFailureIsNotATear() {
    // An element is torn when it has nothing left to carry, and that is *not* the
    // same as having lost a point. Under a strain gradient the far side reaches its
    // failure strain first, and calling that a tear would delete an element that is
    // still carrying most of its load. Mutation testing found this: `torn` could be
    // "any point failed" and nothing noticed, because every tearing test until now
    // strained the element uniformly, where all eight points fail on the same step.
    const Material material = steel();
    double nodes[solidshell::kDof];
    plateElement(0.05, 0.05, 0.02, nodes);

    solidshell::ElementPlasticState state;
    solidshell::initialisePlasticState(nodes, material, state);

    double current[solidshell::kDof], force[solidshell::kDof];
    int partialStep = -1, partialCount = 0;
    double partialForce = 0.0;
    for (int step = 1; step <= 1000; ++step) {
        // A membrane stretch with an in-plane gradient on it, so the +xi points
        // strain harder than the -xi ones and reach the failure strain first.
        const double mean = 0.0015 * step;
        for (int a = 0; a < solidshell::kNodes; ++a) {
            const double x = nodes[a * 3 + 0];
            // u_x = x * mean * (1 + 12 x), so eps_xx = mean (1 + 24 x): +-60% across
            // the element. The factor in the *displacement* is half the one in the
            // strain, and using the strain's factor here put the far side into
            // compression, where it never fails at all.
            const double stretch = mean * (1.0 + 12.0 * x);
            current[a * 3 + 0] = x * (1.0 + stretch);
            current[a * 3 + 1] = nodes[a * 3 + 1] * (1.0 - 0.5 * stretch);
            current[a * 3 + 2] = nodes[a * 3 + 2] * (1.0 - 0.5 * stretch);
        }
        const solidshell::PlasticUpdate result = solidshell::elementPlasticUpdate(
            nodes, current, material, solidshell::Formulation::SolidShell, state, force);
        if (result.failedPoints > 0 && result.failedPoints < solidshell::kGauss &&
            partialStep < 0) {
            partialStep = step;
            partialCount = result.failedPoints;
            partialForce = maxAbs(force, solidshell::kDof);
            expectTrue("a partly failed element is not torn", !state.torn);
            // ... and it has dropped its enhanced strains, which is what lets the
            // surviving points go on accumulating damage instead of stalling.
            expectTrue("a partly failed element has dropped its enhanced strains",
                       maxAbs(state.enhanced, solidshell::kEas) == 0.0);
            expectTrue("and its enhanced-strain solve is well posed again", result.converged);
        }
        if (state.torn) break;
    }
    expectTrue("the gradient did produce a partly failed element", partialStep > 0);
    expectTrue("and it still carried force", partialForce > 1.0e4);
    std::printf("     partial failure at step %d with %d of %d points, still carrying %.2e N\n",
                partialStep, partialCount, solidshell::kGauss, partialForce);

    // That the element eventually tears *completely* is asserted in
    // `testElementTears`, on a uniform stretch, and not here. To drive every point
    // of a 50 mm element past its failure strain through an in-plane gradient the
    // leading corner has to reach nearly 200% stretch, which collapses the element
    // to a fortieth of its thickness -- a geometry a small-strain co-rotational
    // formulation does not represent, so whatever it reported there would not be a
    // statement about the failure criterion.
}

// --- cost -----------------------------------------------------------------------

void testCost() {
    std::printf("\n   measured cost, one core\n");
    const Material material = steel();
    const StructuralMaterial structural = ah36Steel();

    const auto time = [](const char* label, int repeats, auto&& body) {
        const auto begin = std::chrono::steady_clock::now();
        body(repeats);
        const auto end = std::chrono::steady_clock::now();
        const double nanoseconds =
            std::chrono::duration<double, std::nano>(end - begin).count() / repeats;
        std::printf("     %-52s %8.1f ns\n", label, nanoseconds);
        return nanoseconds;
    };

    double sink = 0.0;
    const double elasticPoint = time("return map, elastic point", 2000000, [&](int n) {
        for (int i = 0; i < n; ++i) {
            State state;
            double stress[kV];
            update(material, kNeverFails, V6(1e-4, -3e-5, -3e-5), state, stress);
            sink += stress[0];
        }
    });
    const double plasticPoint =
        time("return map, plastic point (Swift, Newton)", 2000000, [&](int n) {
            for (int i = 0; i < n; ++i) {
                State state;
                double stress[kV];
                update(material, kNeverFails, V6(0.02, -0.006, -0.006), state, stress);
                sink += stress[0];
            }
        });
    const double withTangent =
        time("return map, plastic point with the tangent", 2000000, [&](int n) {
            for (int i = 0; i < n; ++i) {
                State state;
                double stress[kV], tangent[kV * kV];
                update(material, kNeverFails, V6(0.02, -0.006, -0.006), state, stress, tangent);
                sink += stress[0] + tangent[0];
            }
        });

    double nodes[solidshell::kDof];
    plateElement(0.05, 0.05, 0.02, nodes);
    double stiffness[solidshell::kDof * solidshell::kDof];
    solidshell::elementStiffness(nodes, structural, solidshell::Formulation::SolidShell, stiffness);
    double current[solidshell::kDof];
    for (int i = 0; i < solidshell::kDof; ++i) current[i] = nodes[i] * 1.0002;
    double force[solidshell::kDof];

    const double elasticElement =
        time("solid-shell elastic internal force (the reference)", 500000, [&](int n) {
            for (int i = 0; i < n; ++i) {
                solidshell::internalForce(stiffness, nodes, current, force);
                sink += force[0];
            }
        });

    solidshell::ElementPlasticState fresh;
    solidshell::initialisePlasticState(nodes, material, fresh);
    const double elasticPath =
        time("solid-shell elastoplastic update, nothing yielding", 200000, [&](int n) {
            for (int i = 0; i < n; ++i) {
                solidshell::ElementPlasticState state = fresh;
                solidshell::elementPlasticUpdate(nodes, current, material,
                                                 solidshell::Formulation::SolidShell, state,
                                                 force);
                sink += force[0];
            }
        });

    // The case an explicit solver actually runs: one element carried forward, its
    // enhanced parameters warm, taking small plastic increments. Timing a *cold*
    // start on a large jump instead -- which the obvious benchmark does -- charges
    // the Newton four or five iterations it would never need in a march, and
    // overstates the plastic path by about half. Restarted every 200 increments so
    // the strain stays inside the failure strain.
    double yielding[solidshell::kDof];
    const double plasticPath =
        time("solid-shell elastoplastic update, marching and yielding", 100000, [&](int n) {
            solidshell::ElementPlasticState state = fresh;
            for (int i = 0; i < n; ++i) {
                const int within = i % 200;
                if (within == 0) state = fresh;
                const double stretch = 0.0001 * (within + 1);
                for (int a = 0; a < solidshell::kNodes; ++a) {
                    yielding[a * 3 + 0] = nodes[a * 3 + 0] * (1.0 + stretch);
                    yielding[a * 3 + 1] = nodes[a * 3 + 1] * (1.0 - 0.5 * stretch);
                    yielding[a * 3 + 2] = nodes[a * 3 + 2] * (1.0 - 0.5 * stretch);
                }
                solidshell::elementPlasticUpdate(nodes, yielding, material,
                                                 solidshell::Formulation::SolidShell, state,
                                                 force);
                sink += force[0];
            }
        });

    // How much of that is re-forming B on the rest geometry, which never changes.
    // `elementSize` is computeForms plus a little arithmetic, so it isolates the
    // term almost exactly.
    const double forms = time("of which: re-forming B (elementSize proxy)", 200000, [&](int n) {
        for (int i = 0; i < n; ++i) {
            double a = 0.0, b = 0.0;
            solidshell::elementSize(nodes, &a, &b);
            sink += a + b;
        }
    });

    std::printf("     per integration point: elastic %.0f ns, plastic %.0f ns (%.1fx); "
                "the tangent adds %.0f%%\n",
                elasticPoint, plasticPoint, plasticPoint / elasticPoint,
                100.0 * (withTangent / plasticPoint - 1.0));
    std::printf("     per element: elastic reference %.0f ns, elastoplastic %.0f ns elastic-path "
                "/ %.0f ns marching (%.1fx / %.1fx)\n",
                elasticElement, elasticPath, plasticPath, elasticPath / elasticElement,
                plasticPath / elasticElement);
    std::printf("     of the %.0f ns, %.0f ns (%.0f%%) is re-forming B on the rest geometry, "
                "which is loop invariant\n",
                plasticPath, forms, 100.0 * forms / plasticPath);
    std::printf("     per-point state: %zu bytes; 10^6 elements x 8 points = %.0f MB\n",
                sizeof(State), static_cast<double>(sizeof(State)) * 8.0e6 / (1024.0 * 1024.0));

    expectTrue("the cost measurement did work that could not be optimised away",
               std::isfinite(sink) && sink != 0.0);
    // **Deliberately loose, and the loose version replaced a tight one that was a
    // real flaky test.** This asserted `plasticPath > elasticPath` -- true by 35% on
    // an idle machine and measured here at 25 failures in 48 runs with sixteen
    // copies of this suite running at once, because under contention the two
    // converge and the inequality becomes noise. It cost the gate a repeat run, and
    // it had already reported one mutant killed that was not. `test_solid_shell.cpp`
    // says the same thing in its own cost test and says it first: a timing
    // assertion tight enough to be interesting is a flaky test on a shared machine.
    // Two orders of headroom each, so these fire only if something has gone
    // structurally wrong; the numbers themselves are printed, which is what they
    // are for.
    expectTrue("a return-map evaluation is under 10 microseconds", plasticPoint < 1.0e4);
    expectTrue("an elastoplastic element update is under a millisecond",
               plasticPath < 1.0e6 && elasticPath < 1.0e6);
}

// --- What a yielded point has left ------------------------------------------------
//
// Two closed forms and a discontinuity. The secant is what `coupling.hpp` §5 hands
// a reduced model, so the property that matters most is not its value at any one
// strain but that it leaves G **exactly** at zero plastic strain: a coupling that
// is exact on an unyielded zone stops being exact the moment this returns
// `G * (1 - 1e-17)`.
void testWhatAYieldedPointHasLeft() {
    std::printf("\n--- the modulus a yielded point has left ---\n");
    const Material steel = shipSteel();
    const double shear = steel.shearModulus(), bulk = steel.bulkModulus();

    // 1. Zero plastic strain returns the elastic modulus to the last bit, and it has
    //    to be swept rather than asserted once.
    //
    //    **`1/(1/x)` is not the identity in floating point, but it is the identity
    //    for most doubles**, and AH36's shear modulus is one of them. So the version
    //    of this that checked ship steel alone passed against an implementation with
    //    the exact-zero shortcut removed -- mutation testing found exactly that, and
    //    it is the shortcut the whole negative control in `test_coupling.cpp` rests
    //    on: one unit in the last place makes the knockdown ratio 0.999...89 instead
    //    of 1, which emits a block of near-zeros where there should be no block.
    //    206 GPa survives the round trip; **127 GPa does not**, in either modulus,
    //    which is why it is in the sweep.
    int swept = 0, luckyRoundTrips = 0;
    for (double youngs : {69.0e9, 110.0e9, 127.0e9, 200.0e9, 206.0e9, 310.0e9}) {
        Material other = steel;
        other.youngsModulus = youngs;
        const double otherShear = other.shearModulus();
        expectTrue("an unyielded point's secant shear modulus is G, bit for bit",
                   secantShearModulus(other, 0.0) == otherShear);
        expectTrue("and so is its secant Young's modulus",
                   secantYoungsModulus(other, 0.0) == youngs);
        expectTrue("a negative plastic strain is clamped to zero rather than extrapolated",
                   secantShearModulus(other, -0.1) == otherShear);
        if (1.0 / (1.0 / otherShear) == otherShear && 1.0 / (1.0 / youngs) == youngs)
            ++luckyRoundTrips;
        ++swept;
    }
    // The guard that stops the sweep being vacuous: at least one of those moduli has
    // to be one the reciprocal round trip does *not* survive, or every case above is
    // satisfied by arithmetic that happens to be exact and the shortcut is untested.
    std::printf("  %d moduli swept, %d of which survive 1/(1/x) by luck\n", swept,
                luckyRoundTrips);
    expectTrue("the sweep contains a modulus the reciprocal round trip loses",
               luckyRoundTrips < swept);

    // 2. The isotropic pair built from (K, G_s) *is* the uniaxial secant. That is the
    //    identity 1/E = 1/(9K) + 1/(3G) carried through the plastic term, and it is
    //    what makes "keep the bulk modulus, soften the shear one" the right split
    //    rather than a plausible one.
    double worstPair = 0, worstDirect = 0;
    for (double plastic : {1.0e-5, 1.0e-4, 1.0e-3, 1.0e-2, 5.0e-2, 0.15}) {
        double youngs = 0, poisson = 0;
        isotropicFromBulkShear(bulk, secantShearModulus(steel, plastic), &youngs, &poisson);
        const double closed = secantYoungsModulus(steel, plastic);
        // ...and the closed form itself, from the stress-strain point: sigma_y over
        // the total strain that produced it. Nothing here reuses the expression
        // under test.
        const double strength = flowStress(steel.flow, plastic);
        const double total = strength / steel.youngsModulus + plastic;
        worstPair = std::max(worstPair, std::fabs(youngs - closed) / closed);
        worstDirect = std::max(worstDirect, std::fabs(closed - strength / total) / closed);
        expectTrue("the secant Poisson ratio moves towards a half and never past it",
                   poisson > steel.poissonRatio && poisson < 0.5);
    }
    std::printf("  (K, G_s) reproduces the uniaxial secant to %.2e relative;\n"
                "  the uniaxial secant is sigma_y/eps to %.2e\n",
                worstPair, worstDirect);
    expectTrue("the isotropic pair is the uniaxial secant, to rounding", worstPair < 1e-15);
    expectTrue("and the uniaxial secant is the stress over the total strain",
               worstDirect < 1e-15);

    // 3. Monotone, and bounded below by nothing physical -- it falls forever, which
    //    is why a caller with a real element has to notice when it reaches zero.
    double previous = shear;
    for (double plastic : {1.0e-4, 1.0e-3, 1.0e-2, 0.1, 1.0}) {
        const double got = secantShearModulus(steel, plastic);
        expectTrue("the secant modulus falls with plastic strain", got < previous && got > 0.0);
        previous = got;
    }

    // 4. **The tangent is discontinuous at first yield and the secant is not.** This
    //    is the whole reason `coupling.hpp` §5's control loses: an elastic point has
    //    modulus G, and the tangent's limit as the plastic strain goes to zero is
    //    not G but a finite fraction of it.
    const double tangentAtOnset = tangentShearModulus(steel, 0.0);
    const double secantJustPast = secantShearModulus(steel, 1.0e-9);
    std::printf("  at the first increment of flow: G_t/G = %.4f, G_s/G = %.9f\n",
                tangentAtOnset / shear, secantJustPast / shear);
    expectTrue("the tangent drops by a finite step at first yield",
               tangentAtOnset < 0.1 * shear && tangentAtOnset > 0.0);
    expectTrue("the secant leaves G continuously", secantJustPast > 0.999999 * shear);
    expectTrue("so the two disagree by more than an order of magnitude where the softening "
               "is smallest",
               secantJustPast > 10.0 * tangentAtOnset);

    // 5. Perfect plasticity has no tangent stiffness at all, taken as the limit
    //    rather than divided by. The secant on the same curve is still finite.
    Material flat = steel;
    flat.flow = linearHardening(355.0e6, 0.0);
    expectTrue("a perfectly plastic curve has no tangent shear stiffness",
               tangentShearModulus(flat, 0.01) == 0.0);
    expectTrue("but it still has a secant one", secantShearModulus(flat, 0.01) > 0.0);

    // 6. On a linear curve both are closed forms of their own, which is the check
    //    that the Swift arithmetic above is not being marked against itself.
    Material linear = steel;
    const double hardening = 2.0e9;
    linear.flow = linearHardening(355.0e6, hardening);
    const double e = linear.youngsModulus;
    const double wantTangent = 1.0 / (1.0 / linear.shearModulus() + 3.0 / hardening);
    expectNear("the linear-curve tangent shear modulus is 1/(1/G + 3/H)",
               tangentShearModulus(linear, 0.03), wantTangent, 1e-9 * wantTangent);
    const double plastic = 0.03, sigma = 355.0e6 + hardening * plastic;
    expectNear("and the linear-curve secant Young's modulus is sigma/(sigma/E + eps_p)",
               secantYoungsModulus(linear, plastic), sigma / (sigma / e + plastic),
               1e-9 * sigma / (sigma / e + plastic));

    // 7. `isotropicFromBulkShear` on the *elastic* pair is the identity, to rounding
    //    -- the round trip, which is why `coupling.hpp` §5 emits no block at all for
    //    an unyielded element rather than a block of zeros.
    double youngs = 0, poisson = 0;
    isotropicFromBulkShear(bulk, shear, &youngs, &poisson);
    std::printf("  the elastic round trip is off by %.2e in E and %.2e in nu, which is why an "
                "unyielded element gets no block\n",
                std::fabs(youngs - steel.youngsModulus) / steel.youngsModulus,
                std::fabs(poisson - steel.poissonRatio));
    expectTrue("the round trip through (K, G) recovers E and nu to rounding",
               std::fabs(youngs - steel.youngsModulus) < 1e-12 * steel.youngsModulus &&
                   std::fabs(poisson - steel.poissonRatio) < 1e-14);
}

}  // namespace

void runPlasticityTests() {
    std::printf("\n--- plasticity and ductile failure ---\n");
    testFlowCurves();
    testUniaxialTension();
    testElasticStressAgreesWithTheModuli();
    testTheReturnIsIdempotent();
    testStepIndependence();
    testHydrostaticStatesDoNotYield();
    testConsistencyAndIncompressibility();
    testUnloadingIsElastic();
    testRotationalInvariance();
    testPureShearYieldsAtYieldOverRootThree();
    testBauschingerEffect();
    testAlgorithmicTangent();
    testFailureStrainRegularisation();
    testFailureIsReachedAndIsIrreversible();
    testLocalisationAndItsLimit();
    testElementSize();
    testElementReducesToTheElasticElement();
    testElementUniformDeformationIsAPatchTest();
    testElementBendingIntoYield();
    testElementFrameIndifference();
    testElementTears();
    testPartialFailureIsNotATear();
    testWhatAYieldedPointHasLeft();
    testCost();
}
