// SPDX-License-Identifier: MIT
//
// Validation of the radiation hydrodynamics.
//
// The Ogilvie relations are exact identities, so most of this file is aimed at
// them rather than at the coefficients themselves. Three instruments do the
// work, and they are independent of each other:
//
//   * **Closed forms.** K(t) = e^{-at} has B(omega) = a/(a^2 + omega^2) and
//     A(omega) = A_inf - 1/(a^2 + omega^2) exactly; a heaving semicircle has an
//     infinite-frequency added mass of exactly rho pi a^2 / 2 by an image
//     argument that involves no hydrodynamics at all; a circular section rolling
//     about its own centre displaces nothing and must have identically zero roll
//     coefficients; a prismatic hull's strip integrals are algebra.
//   * **Two routes to the same number.** Radiation damping computed by
//     integrating pressure over the hull, against the same damping computed from
//     the amplitude of the wave radiated to infinity. A_inf from a rigid-lid
//     panel solve, against A_inf recovered from A(omega) and K(t) through
//     Ogilvie. Neither pair shares any code path beyond the section geometry.
//   * **Froude scaling.** Scale the hull by lambda and the frequencies by
//     1/sqrt(lambda) and every nondimensional coefficient must be *identical*,
//     because the 2D problem depends only on nu times a length. That is what
//     catches a power of a length in the wrong place, which is how the Ikeda
//     nondimensionalisation was found to be wrong in this repo before.
//
// Published comparison: Ursell (1949), the heaving semi-immersed circular
// cylinder, whose added-mass coefficient dips to about 0.60 of rho pi a^2 / 2
// near omega sqrt(a/g) = 0.9 and returns to 1.0 at high frequency. Both the
// value of the minimum and its location are asserted, and the high-frequency
// limit is asserted against the exact 1.0.
#include "engine/sim/radiation.hpp"
#include "harness.hpp"

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

using namespace sim;
using testing::expectEqual;
using testing::expectNear;
using testing::expectTrue;

namespace {

// --- Test geometry -----------------------------------------------------------

// A semi-immersed circular cylinder of radius a, as a Lewis form. sigma = pi/4
// is exactly the semicircle's area coefficient, so a1 and a3 must come out zero.
LewisSection semicircle(double a) { return lewisSection(2.0 * a, a, kPi / 4.0); }

// 170 m ro-pax, matching the hull tests/test_roll_damping.cpp uses: beam 25 m,
// draft 6.5 m, full amidships and fine at the ends. Nine stations and twenty
// frequencies keep the whole suite inside its time budget; the convergence
// behaviour is measured separately rather than assumed.
RadiationHull ferry(double scale = 1.0) {
    RadiationHull hull;
    hull.draft = 6.5 * scale;
    hull.density = kRhoSeawater;
    hull.panelsPerHalfSection = 24;
    const int stations = 9;
    for (int i = 0; i < stations; ++i) {
        const double u = -1.0 + 2.0 * i / (stations - 1);   // -1 aft, +1 forward
        RadiationStation s;
        s.x = 85.0 * scale * u;
        // Finer forward than aft, so the heave-pitch and sway-yaw couplings are
        // genuinely non-zero: a fore-and-aft symmetric hull makes them numerical
        // dust and any test on them vacuous.
        s.beam = 25.0 * scale * std::max(0.25, std::pow(1.0 - u * u, 0.35)) *
                 (1.0 - 0.18 * u);
        s.draft = 6.5 * scale;
        s.areaCoefficient = 0.98 - 0.35 * u * u;
        hull.stations.push_back(s);
    }
    return hull;
}

// A prismatic barge: every station identical, so the strip integrals collapse to
// algebra and the length integration can be checked against a closed form.
RadiationHull barge(int stations, double length) {
    RadiationHull hull;
    hull.draft = 4.0;
    hull.density = 1000.0;
    hull.panelsPerHalfSection = 20;
    for (int i = 0; i < stations; ++i) {
        RadiationStation s;
        s.x = -0.5 * length + length * i / (stations - 1);
        s.beam = 12.0;
        s.draft = 4.0;
        s.areaCoefficient = 0.90;
        hull.stations.push_back(s);
    }
    return hull;
}

// Built once: a full table costs a few hundred milliseconds and every test below
// wants the same one.
const RadiationTable& ferryTable() {
    static const RadiationHull hull = ferry();
    static const std::vector<double> grid = radiationFrequencyGrid(0.12, 4.0, 20);
    static const RadiationTable table = stripTheoryTable(hull, grid);
    return table;
}

double peakDamping(const RadiationTable& table, int i, int j) {
    double peak = 0;
    for (int f = 0; f < table.size(); ++f)
        peak = std::max(peak, std::abs(table.damping[static_cast<std::size_t>(f)]
                                                    [static_cast<std::size_t>(i)]
                                                    [static_cast<std::size_t>(j)]));
    return peak;
}

// --- Lewis forms -------------------------------------------------------------

void testLewisFormsReproduceTheirInputs() {
    struct Case { double beam, draft, sigma; };
    const Case cases[] = {{25.0, 6.5, 0.98}, {12.0, 4.0, 0.90}, {6.0, 6.0, 0.70},
                          {30.0, 5.0, 0.85}, {8.0, 10.0, 0.75}};
    for (const Case& c : cases) {
        const LewisSection s = lewisSection(c.beam, c.draft, c.sigma);
        // Beam and draft come out of the mapping at theta = pi/2 and theta = 0.
        expectNear("Lewis half-beam at the waterline", s.point(0.5 * kPi).x, 0.5 * c.beam,
                   1e-12 * c.beam);
        expectNear("Lewis draft at the keel", -s.point(0.0).y, c.draft, 1e-12 * c.draft);
        // The analytic area of the mapping must equal what was asked for, which
        // is the whole point of solving the quadratic rather than guessing a3.
        expectNear("Lewis sectional area", s.sectionArea(), c.sigma * c.beam * c.draft,
                   1e-10 * c.sigma * c.beam * c.draft);
        expectTrue("this Lewis form is representable", validateLewisSection(s).empty());
    }

    // sigma = pi/4 at B = 2T is the semicircle exactly, and the mapping has to
    // find that rather than the other root of the quadratic -- which lands on
    // a3 = -1/2 and produces a plausible-looking section that is not a circle.
    const LewisSection circle = semicircle(1.0);
    expectNear("the semicircle has a1 = 0", circle.a1, 0.0, 1e-14);
    expectNear("the semicircle has a3 = 0", circle.a3, 0.0, 1e-14);
    expectNear("the semicircle's scale is its radius", circle.scale, 1.0, 1e-14);
    for (int i = 0; i <= 8; ++i) {
        const double theta = 0.5 * kPi * i / 8.0;
        const Vec3 p = circle.point(theta);
        expectNear("the semicircle's contour is at radius a", std::hypot(p.x, p.y), 1.0, 1e-14);
    }

    // Outside the family, and folded over: both must be reported, not silently
    // approximated.
    // A Lewis form may bulge outside the B x T rectangle, so sigma > 1 is not
    // itself unattainable; for B/2T = 2 the family runs out at about 1.23.
    // **"Both must be reported" is what the comment says and what the assertions
    // did not check.** `folded` forces `a3` past univalence *after* the fit, which
    // also throws the sectional area off what was asked for -- so the area guard
    // fires as well, and deleting the univalence check entirely leaves the caller
    // told that sigma = 0.90 is unattainable for this beam/draft ratio. That is a
    // real and different failure, the one `tooFull` is here to produce, and a
    // reader who believed it would go looking for the Lewis family's limits
    // instead of a contour that crosses itself.
    const auto mentions = [](const std::vector<std::string>& v, const char* s) {
        for (const std::string& p : v)
            if (p.find(s) != std::string::npos) return true;
        return false;
    };

    const std::vector<std::string> tooFullSaid = validateLewisSection(lewisSection(20.0, 5.0, 1.45));
    expectTrue("an unattainable area coefficient is reported", !tooFullSaid.empty());
    expectTrue("naming the family it fell outside of",
               mentions(tooFullSaid, "outside the Lewis family"));

    LewisSection folded = lewisSection(20.0, 5.0, 0.90);
    folded.a3 = 0.6;   // forced past univalence
    const std::vector<std::string> foldedSaid = validateLewisSection(folded);
    expectTrue("a folded Lewis form is reported", !foldedSaid.empty());
    expectTrue("and the fold is named, not just the area it threw off",
               mentions(foldedSaid, "not univalent"));
}

// --- The 2D section solver, against Ursell -----------------------------------

void testHeavingCircleAgainstUrsell() {
    const LewisSection circle = semicircle(1.0);
    const double density = 1000.0;
    const double reference = density * kPi / 2.0;   // rho pi a^2 / 2 with a = 1

    // The omega -> infinity limit is exact and needs no reference at all: with
    // phi = 0 on the free surface the section plus its negative image is a full
    // circle translating in unbounded fluid, whose added mass is rho pi a^2, and
    // half of that belongs to the submerged half.
    const SectionCoefficients rigid = sectionCoefficientsInfiniteFrequency(circle, density, 160);
    expectNear("semicircle added mass at infinite frequency is rho pi a^2 / 2",
               rigid.addedMass[kSectionHeave][kSectionHeave] / reference, 1.0, 0.01);
    expectNear("there is no damping at infinite frequency",
               rigid.damping[kSectionHeave][kSectionHeave], 0.0, 0.0);

    // Ursell's curve: a minimum near 0.60 at omega sqrt(a/g) close to 0.9.
    double lowest = 1e30, lowestAt = 0;
    for (int i = 0; i <= 16; ++i) {
        const double nuA = 0.4 + 0.05 * i;
        const SectionCoefficients c =
            sectionCoefficients(circle, std::sqrt(nuA * kGravity), density, 80);
        const double ratio = c.addedMass[kSectionHeave][kSectionHeave] / reference;
        if (ratio < lowest) { lowest = ratio; lowestAt = nuA; }
    }
    expectNear("Ursell's added-mass minimum for a heaving semicircle", lowest, 0.60, 0.02);
    expectNear("and where it sits, in omega sqrt(a/g)", std::sqrt(lowestAt), 0.9, 0.08);

    // Energy: damping from the pressure integral against damping from the wave
    // radiated to infinity. Independent routes, so the agreement measures the
    // discretisation and nothing else. Refining the panels must improve it, and
    // the improvement must be first order -- flat constant-strength panels with
    // midpoint collocation cannot do better.
    double coarse = 0, fine = 0;
    for (double nuA : {0.2, 0.5, 1.0, 2.0}) {
        const double omega = std::sqrt(nuA * kGravity);
        coarse = std::max(coarse, sectionCoefficients(circle, omega, density, 40).energyResidual);
        fine = std::max(fine, sectionCoefficients(circle, omega, density, 80).energyResidual);
    }
    expectTrue("the near-field/far-field energy balance holds to better than 1%", coarse < 0.01);
    expectTrue("and doubling the panels roughly halves the error",
               fine < 0.7 * coarse && fine > 0.3 * coarse);
    // Reciprocity makes A_ij = A_ji exactly in potential flow. Constant sources
    // with midpoint collocation do not preserve it, so the departure is a pure
    // discretisation error and has to vanish with the panels -- a far stronger
    // statement than any fixed tolerance, and it needs a section whose coupling
    // is not zero by symmetry to say anything at all.
    const LewisSection wedge = lewisSection(12.0, 4.0, 0.90);
    double previousReciprocity = 0;
    for (int panels : {20, 40, 80}) {
        const SectionCoefficients c = sectionCoefficients(wedge, 0.8, 1000.0, panels);
        const double mean = 0.5 * (c.addedMass[kSectionSway][kSectionRoll] +
                                   c.addedMass[kSectionRoll][kSectionSway]);
        expectTrue("the sway-roll coupling is not zero, so reciprocity says something",
                   std::abs(mean) > 0);
        const double error = std::abs(c.addedMass[kSectionSway][kSectionRoll] -
                                      c.addedMass[kSectionRoll][kSectionSway]) / std::abs(mean);
        if (previousReciprocity > 0)
            expectTrue("reciprocity error falls at least as fast as the panel size",
                       error < 0.65 * previousReciprocity);
        previousReciprocity = error;
    }
    expectTrue("and is under 0.1% by 80 panels", previousReciprocity < 1e-3);
    std::printf("     reciprocity (A24 vs A42) %.2e at 80 panels, first order in panel size\n",
                previousReciprocity);

    std::printf("     semicircle: Ca(inf) %.4f (exact 1), min Ca %.4f at omega sqrt(a/g) %.3f;"
                " energy residual %.2e at 40 panels, %.2e at 80\n",
                rigid.addedMass[kSectionHeave][kSectionHeave] / reference, lowest,
                std::sqrt(lowestAt), coarse, fine);
}

// A circular section rolling about its own centre pushes no water anywhere: the
// boundary velocity is everywhere tangential. Both roll coefficients must be
// *identically* zero, not small, and the sway-roll coupling with them. Nothing
// about the solver knows this, so it is a real check on the generalised normals.
void testCircleRollingAboutItsCentreDisplacesNothing() {
    const LewisSection circle = semicircle(3.0);
    const SectionCoefficients c = sectionCoefficients(circle, 1.2, kRhoSeawater, 40);
    expectTrue("the same section has non-trivial heave and sway, so this is not vacuous",
               c.addedMass[kSectionHeave][kSectionHeave] > 0 &&
                   c.addedMass[kSectionSway][kSectionSway] > 0);
    const double scale = c.addedMass[kSectionSway][kSectionSway] * 9.0;   // times a^2
    expectNear("a circle rolling about its centre has no roll added inertia",
               c.addedMass[kSectionRoll][kSectionRoll], 0.0, 1e-16 * scale);
    expectNear("nor any roll damping", c.damping[kSectionRoll][kSectionRoll], 0.0, 1e-16 * scale);
    expectNear("nor any sway-roll coupling", c.addedMass[kSectionSway][kSectionRoll], 0.0,
               1e-16 * scale);
}

// --- Strip integration -------------------------------------------------------

// A prismatic hull's strip integrals are algebra: every station has the same
// sectional coefficient a, so A33 = a * sum(w) = a * L exactly, A55 = a * sum(w
// x^2), and A35 = -a * sum(w x) = 0 by symmetry. The trapezium weights are
// re-derived here rather than taken from the engine.
void testPrismaticHullIntegratesAlongItsLength() {
    const double length = 60.0;
    const int stations = 9;
    const RadiationHull hull = barge(stations, length);
    const std::vector<double> grid{0.8};
    const RadiationTable table = stripTheoryTable(hull, grid);

    const LewisSection section = lewisSection(12.0, 4.0, 0.90);
    SectionCoefficients c = sectionCoefficients(section, 0.8, 1000.0, 20);
    // The hull matrices are about the baseline, the section about the waterline.
    const double a33 = c.addedMass[kSectionHeave][kSectionHeave];
    const double b33 = c.damping[kSectionHeave][kSectionHeave];

    double sumW = 0, sumWx = 0, sumWxx = 0;
    const double h = length / (stations - 1);
    for (int i = 0; i < stations; ++i) {
        const double x = -0.5 * length + h * i;
        const double w = (i == 0 || i == stations - 1) ? 0.5 * h : h;
        sumW += w;
        sumWx += w * x;
        sumWxx += w * x * x;
    }
    expectNear("the trapezium weights sum to the length", sumW, length, 1e-12 * length);
    expectNear("and their first moment vanishes about midship", sumWx, 0.0, 1e-12 * length * length);

    expectTrue("the barge has non-trivial sectional coefficients", a33 > 0 && b33 > 0);
    expectNear("A33 is the sectional value times the length", table.addedMass[0][2][2],
               a33 * sumW, 1e-9 * a33 * sumW);
    expectNear("B33 likewise", table.damping[0][2][2], b33 * sumW, 1e-9 * b33 * sumW);
    expectNear("A55 is the second moment of the sectional heave added mass",
               table.addedMass[0][4][4], a33 * sumWxx, 1e-9 * a33 * sumWxx);
    expectNear("A35 vanishes on a fore-and-aft symmetric hull", table.addedMass[0][2][4], 0.0,
               1e-9 * a33 * sumWxx / length);
    expectNear("A66 is the second moment of the sectional sway added mass",
               table.addedMass[0][5][5],
               c.addedMass[kSectionSway][kSectionSway] * sumWxx,
               1e-9 * c.addedMass[kSectionSway][kSectionSway] * sumWxx);
    // Strip theory has nothing to say about surge, and says so rather than
    // inventing a number.
    expectNear("surge added mass is left at zero by strip theory", table.addedMass[0][0][0], 0.0,
               0.0);

    // Roll about the baseline is *not* roll about the waterline: it picks up the
    // draft lever on the sway mode. Re-derived here from the section values.
    const double t = hull.draft;
    // The panel solve does not return an exactly symmetric coupling, and the
    // transfer symmetrises it. That asymmetry is a discretisation error and is
    // asserted small rather than assumed away.
    const double coupling =
        0.5 * (c.addedMass[kSectionSway][kSectionRoll] + c.addedMass[kSectionRoll][kSectionSway]);
    expectTrue("the sway-roll coupling is symmetric to better than 0.5% at 20 panels",
               std::abs(c.addedMass[kSectionSway][kSectionRoll] -
                        c.addedMass[kSectionRoll][kSectionSway]) < 5e-3 * std::abs(coupling));
    const double rollAboutBaseline = c.addedMass[kSectionRoll][kSectionRoll] -
                                     2.0 * t * coupling +
                                     t * t * c.addedMass[kSectionSway][kSectionSway];
    expectNear("A44 is transferred from the waterline to the baseline", table.addedMass[0][3][3],
               rollAboutBaseline * sumW, 1e-9 * std::abs(rollAboutBaseline * sumW));
    expectTrue("the transfer actually moved something",
               std::abs(rollAboutBaseline - c.addedMass[kSectionRoll][kSectionRoll]) >
                   0.1 * std::abs(rollAboutBaseline));
}

// --- Froude scaling ----------------------------------------------------------
//
// The 2D radiation problem depends on the section only through nu times a
// length, so scaling the hull by lambda and the frequencies by 1/sqrt(lambda)
// must leave every nondimensional coefficient bit-for-bit the same and every
// dimensional one scaled by a fixed power of lambda. A single misplaced power of
// a length shows up here and nowhere else.
void testGeometricScalingIsExact() {
    const double lambda = 2.5;
    const std::vector<double> small = radiationFrequencyGrid(0.2, 3.0, 8);
    std::vector<double> large;
    for (double omega : small) large.push_back(omega / std::sqrt(lambda));

    const RadiationTable a = stripTheoryTable(ferry(1.0), small);
    const RadiationTable b = stripTheoryTable(ferry(lambda), large);
    expectEqual("both tables have the same number of frequencies", a.size(), b.size());

    struct Entry { const char* name; int i, j; double massPower, dampingPower; };
    const Entry entries[] = {
        {"heave", 2, 2, 3.0, 2.5},
        {"sway", 1, 1, 3.0, 2.5},
        {"roll", 3, 3, 5.0, 4.5},
        {"pitch", 4, 4, 5.0, 4.5},
        {"yaw", 5, 5, 5.0, 4.5},
        {"sway-roll", 1, 3, 4.0, 3.5},
        {"heave-pitch", 2, 4, 4.0, 3.5},
    };
    double worstMass = 0, worstDamping = 0;
    for (const Entry& e : entries) {
        const double massFactor = std::pow(lambda, e.massPower);
        const double dampingFactor = std::pow(lambda, e.dampingPower);
        bool sawSomething = false;
        for (int f = 0; f < a.size(); ++f) {
            const std::size_t k = static_cast<std::size_t>(f);
            const std::size_t i = static_cast<std::size_t>(e.i), j = static_cast<std::size_t>(e.j);
            const double smallMass = a.addedMass[k][i][j];
            const double largeMass = b.addedMass[k][i][j];
            const double smallDamping = a.damping[k][i][j];
            const double largeDamping = b.damping[k][i][j];
            // Only entries that carry real weight: a ratio taken on numerical
            // dust says nothing, and the "is non-zero somewhere" guard below is
            // what keeps this from passing vacuously.
            if (std::abs(smallMass) > 1e-6 * std::abs(a.addedMass[k][i][i])) {
                sawSomething = true;
                worstMass = std::max(worstMass,
                                     std::abs(largeMass / (smallMass * massFactor) - 1.0));
            }
            if (std::abs(smallDamping) > 1e-12 * peakDamping(a, e.i, e.j))
                worstDamping = std::max(
                    worstDamping, std::abs(largeDamping / (smallDamping * dampingFactor) - 1.0));
        }
        expectTrue(std::string("scaling test: ") + e.name + " is non-zero somewhere", sawSomething);
    }
    expectTrue("added mass scales as lambda^n exactly", worstMass < 1e-9);
    expectTrue("damping scales as lambda^(n-0.5) exactly", worstDamping < 1e-9);

    const double infRatio = b.addedMassInfinite[2][2] / a.addedMassInfinite[2][2];
    expectNear("A_inf scales as lambda^3", infRatio, std::pow(lambda, 3.0),
               1e-9 * std::pow(lambda, 3.0));

    // Negative control. If the frequency is *not* Froude scaled with the hull,
    // the nondimensional coefficients must change -- otherwise the test above
    // would pass against an implementation that ignored frequency entirely.
    const RadiationTable c = stripTheoryTable(ferry(lambda), small);
    expectTrue("the same hull at unscaled frequencies is a different answer",
               std::abs(c.addedMass[0][2][2] / (a.addedMass[0][2][2] * std::pow(lambda, 3.0)) -
                        1.0) > 0.05);
    std::printf("     geometric scaling by %.1f: added mass to %.1e, damping to %.1e relative\n",
                lambda, worstMass, worstDamping);
}

// --- Ogilvie, against closed forms -------------------------------------------

// K(t) = exp(-a t) has B(omega) = a / (a^2 + omega^2) and
// A(omega) = A_inf - 1 / (a^2 + omega^2), both exactly. K(t) = exp(-a t)
// cos(b t) has B(omega) = (a/2) [1/(a^2 + (omega-b)^2) + 1/(a^2 + (omega+b)^2)].
// Neither involves the hull, so this is a direct test of the transforms.
void testOgilvieTransformsAgainstClosedForms() {
    const double dt = 0.002;
    const int count = 20000;   // 40 s, by which exp(-0.4 t) is 1e-7
    const double decay = 0.4;
    const double frequency = 1.3;

    std::vector<double> plain(static_cast<std::size_t>(count));
    std::vector<double> oscillating(static_cast<std::size_t>(count));
    for (int n = 0; n < count; ++n) {
        const double t = n * dt;
        plain[static_cast<std::size_t>(n)] = std::exp(-decay * t);
        oscillating[static_cast<std::size_t>(n)] = std::exp(-decay * t) * std::cos(frequency * t);
    }

    double worstPlain = 0, worstOscillating = 0, worstMass = 0;
    for (double omega : {0.05, 0.2, 0.5, 1.0, 1.3, 2.0, 4.0}) {
        const double wantPlain = decay / (decay * decay + omega * omega);
        worstPlain = std::max(worstPlain,
                              std::abs(dampingFromRetardation(plain, dt, omega) - wantPlain));
        const double wantOscillating =
            0.5 * decay *
            (1.0 / (decay * decay + (omega - frequency) * (omega - frequency)) +
             1.0 / (decay * decay + (omega + frequency) * (omega + frequency)));
        worstOscillating =
            std::max(worstOscillating,
                     std::abs(dampingFromRetardation(oscillating, dt, omega) - wantOscillating));
        const double wantMass = 5.0 - 1.0 / (decay * decay + omega * omega);
        worstMass = std::max(worstMass,
                             std::abs(addedMassFromRetardation(plain, dt, omega, 5.0) - wantMass));
    }
    expectTrue("B(omega) from an exponential K matches a/(a^2+omega^2)", worstPlain < 2e-6);
    expectTrue("B(omega) from a damped cosine K matches its closed form",
               worstOscillating < 2e-6);
    // The sine transform divides by omega, so the same tail truncation that
    // costs 1e-7 in B costs 20 times that in A at the bottom of the band.
    expectTrue("A(omega) from an exponential K matches A_inf - 1/(a^2+omega^2)",
               worstMass < 1e-5);
    // Guard against a vacuous comparison: the quantity being matched has to move.
    expectTrue("B(omega) actually varies over the band tested",
               dampingFromRetardation(plain, dt, 0.05) >
                   20.0 * dampingFromRetardation(plain, dt, 4.0));
    std::printf("     closed-form transforms: worst error %.1e (B), %.1e (B oscillating),"
                " %.1e (A)\n", worstPlain, worstOscillating, worstMass);
}

// --- Ogilvie, on the real hull -----------------------------------------------

void testRetardationRoundTripsThroughOgilvie() {
    const RadiationTable& table = ferryTable();
    const double dt = 0.05;
    const int count = 1200;   // 60 s
    const std::vector<double> k = retardationFunction(table, 2, 2, dt, count);
    expectTrue("K(0) is positive and large", k[0] > 0);

    const double peak = peakDamping(table, 2, 2);
    expectTrue("the damping curve is not flat", peak > 0);

    double worst = 0, worstAt = 0;
    for (int f = 0; f < table.size(); ++f) {
        const double got = dampingFromRetardation(k, dt, table.omega[static_cast<std::size_t>(f)]);
        const double want = table.damping[static_cast<std::size_t>(f)][2][2];
        if (std::abs(got - want) > worst) {
            worst = std::abs(got - want);
            worstAt = table.omega[static_cast<std::size_t>(f)];
        }
    }
    expectTrue("B -> K -> B round-trips to better than 1% of the peak", worst < 0.01 * peak);
    expectTrue("the round-trip check spanned the grid", worstAt > 0);

    // A ship forgets. K must decay, and the window has to be long enough that
    // the tail is genuinely gone rather than merely off the end.
    const double toOnePercent = memoryDecayTime(k, dt, 0.01);
    const double toATenth = memoryDecayTime(k, dt, 0.001);
    expectTrue("K decays to 1% of its peak well inside the window",
               toOnePercent > 0 && toOnePercent < 0.6 * count * dt);
    expectTrue("and keeps decaying after that", toATenth > toOnePercent);
    expectTrue("the tail of the window is small", std::abs(k[static_cast<std::size_t>(count - 1)]) <
                                                      0.01 * std::abs(k[0]));

    // Ogilvie couples A and B: A_inf recovered from A(omega) and K(t) cannot
    // depend on omega, and must agree with the rigid-lid panel solve, which
    // shares no code with the transform.
    double sum = 0, spread = 0;
    int used = 0;
    std::vector<double> recovered;
    for (int f = 0; f < table.size(); ++f) {
        const double omega = table.omega[static_cast<std::size_t>(f)];
        if (omega < 0.3 || omega > 2.5) continue;
        const double value = infiniteAddedMassFromRetardation(
            k, dt, omega, table.addedMass[static_cast<std::size_t>(f)][2][2]);
        recovered.push_back(value);
        sum += value;
        ++used;
    }
    expectTrue("the A_inf recovery used a range of frequencies", used >= 6);
    const double mean = sum / used;
    for (double value : recovered) spread = std::max(spread, std::abs(value / mean - 1.0));
    expectTrue("A_inf recovered through Ogilvie does not depend on omega", spread < 0.05);
    expectNear("and agrees with the rigid-lid panel solve", mean / table.addedMassInfinite[2][2],
               1.0, 0.03);

    std::printf("     Ogilvie: B->K->B worst %.2e of peak; A_inf %.4e (Ogilvie) vs %.4e"
                " (rigid lid), %.2f%% apart, %.2f%% spread over omega\n",
                worst / peak, mean, table.addedMassInfinite[2][2],
                100.0 * std::abs(mean / table.addedMassInfinite[2][2] - 1.0), 100.0 * spread);
    std::printf("     memory: K falls to 1%% of its peak at %.1f s and 0.1%% at %.1f s\n",
                toOnePercent, toATenth);
}

// --- What damping is allowed to do -------------------------------------------

void testDampingIsPhysical() {
    const RadiationTable& table = ferryTable();
    expectTrue("the table was built", table.size() > 0);

    int negatives = 0;
    double worstNegative = 0;
    for (int f = 0; f < table.size(); ++f)
        for (int d = 1; d < 6; ++d) {
            const double b = table.damping[static_cast<std::size_t>(f)][static_cast<std::size_t>(d)]
                                          [static_cast<std::size_t>(d)];
            if (b < 0) {
                ++negatives;
                worstNegative = std::min(worstNegative, b / std::max(1.0, peakDamping(table, d, d)));
            }
        }
    expectEqual("no diagonal damping is negative anywhere in the table", negatives, 0);
    if (negatives) std::printf("     worst negative: %.3e of peak\n", worstNegative);

    // B must vanish at both ends of the spectrum: a very long wave carries no
    // energy away and a very short one does not reach the hull. Asserted as a
    // fraction of the peak, with the peak required to be interior so the test is
    // not satisfied by a monotone curve.
    for (int d : {1, 2, 3}) {
        const double peak = peakDamping(table, d, d);
        expectTrue("this mode radiates something", peak > 0);
        const std::size_t last = static_cast<std::size_t>(table.size() - 1);
        expectTrue("B -> 0 at low frequency", table.damping[0][static_cast<std::size_t>(d)]
                                                           [static_cast<std::size_t>(d)] <
                                                  0.6 * peak);
        // Roll radiates over a much wider band than heave or sway, so the top of
        // a band chosen for heave still has real roll damping in it. What must
        // hold for every mode is that the curve is past its peak and still
        // falling; the tighter bound is asserted for heave alone below.
        const double top = table.damping[last][static_cast<std::size_t>(d)]
                                        [static_cast<std::size_t>(d)];
        const double nextToTop = table.damping[last - 1][static_cast<std::size_t>(d)]
                                              [static_cast<std::size_t>(d)];
        expectTrue("B is past its peak at the top of the band", top < 0.3 * peak);
        expectTrue("and still falling there", top < nextToTop);
        bool interior = false;
        for (int f = 1; f + 1 < table.size(); ++f)
            if (table.damping[static_cast<std::size_t>(f)][static_cast<std::size_t>(d)]
                             [static_cast<std::size_t>(d)] >= peak)
                interior = true;
        expectTrue("the damping peak is interior, so the curve really does rise and fall",
                   interior);
    }

    expectTrue("heave damping in particular has all but gone at the top of the band",
               table.damping[static_cast<std::size_t>(table.size() - 1)][2][2] <
                   0.05 * peakDamping(table, 2, 2));

    // Added mass tends to A_inf from below at high frequency for heave.
    const std::size_t last = static_cast<std::size_t>(table.size() - 1);
    expectNear("A(omega) approaches A_inf at the top of the band",
               table.addedMass[last][2][2] / table.addedMassInfinite[2][2], 1.0, 0.10);
    expectTrue("and is far from it at low frequency",
               table.addedMass[0][2][2] > 1.5 * table.addedMassInfinite[2][2]);

    // The symmetric ship must not couple the vertical and lateral planes at all.
    for (int f = 0; f < table.size(); ++f) {
        const std::size_t k = static_cast<std::size_t>(f);
        const double scale = table.addedMass[k][2][2];
        expectNear("heave does not couple to sway", table.addedMass[k][2][1], 0.0, 1e-12 * scale);
        expectNear("heave does not couple to roll", table.addedMass[k][2][3], 0.0, 1e-12 * scale);
    }

    // dampingAt() interpolates and, outside the grid, returns zero -- which is
    // the physically right extrapolation rather than a convenience: radiation
    // damping really does vanish at both ends of the spectrum.
    const std::size_t mid = static_cast<std::size_t>(table.size() / 2);
    expectNear("dampingAt reproduces a grid point exactly",
               table.dampingAt(2, 2, table.omega[mid]), table.damping[mid][2][2],
               1e-12 * table.damping[mid][2][2]);
    const double between = 0.5 * (table.omega[mid] + table.omega[mid + 1]);
    const double linear = 0.5 * (table.damping[mid][2][2] + table.damping[mid + 1][2][2]);
    expectTrue("the two bracketing values differ, so interpolating means something",
               std::abs(table.damping[mid][2][2] - table.damping[mid + 1][2][2]) >
                   0.01 * table.damping[mid][2][2]);
    expectNear("and interpolates linearly between them", table.dampingAt(2, 2, between), linear,
               1e-9 * std::abs(linear));
    expectNear("B is zero below the grid", table.dampingAt(2, 2, 0.5 * table.omega.front()), 0.0,
               0.0);
    expectNear("B is zero above the grid", table.dampingAt(2, 2, 2.0 * table.omega.back()), 0.0,
               0.0);

    // **The two grid points the assertions above step over.** "Zero outside the
    // grid" is the right extrapolation and `omega.front()` and `omega.back()` are
    // not outside it -- they are the grid. A caller sampling the table on its own
    // frequencies, which is the natural way to build a retardation kernel, hits
    // both ends every time.
    const double firstB = table.damping.front()[2][2];
    const double lastB = table.damping.back()[2][2];
    expectTrue("the end points carry damping worth reporting, or this proves nothing",
               std::abs(firstB) > 0 && std::abs(lastB) > 0);
    expectNear("dampingAt reproduces the first grid point", table.dampingAt(2, 2, table.omega.front()),
               firstB, 1e-12 * std::abs(firstB));
    expectNear("and the last", table.dampingAt(2, 2, table.omega.back()), lastB,
               1e-12 * std::abs(lastB));

    std::printf("     table: %d frequencies, worst kept energy residual %.2e,"
                " %d of %d section solves rejected as irregular frequencies\n",
                table.size(), table.worstEnergyResidual, table.repairedSolves, table.totalSolves);
}

// --- State space -------------------------------------------------------------

// Prony is exact for a signal that really is a sum of damped sinusoids. Before
// it is trusted on a hull it has to reproduce one it cannot have guessed.
void testStateSpaceRecoversKnownPoles() {
    const double dt = 0.02;
    const int count = 1500;
    struct Mode { double decay, frequency, cosine, sine; };
    const Mode wanted[] = {{0.35, 1.1, 2.0, -0.5}, {0.9, 2.7, -0.75, 0.3}};
    std::vector<double> k(static_cast<std::size_t>(count));
    for (int n = 0; n < count; ++n) {
        const double t = n * dt;
        double value = 0;
        for (const Mode& m : wanted)
            value += std::exp(-m.decay * t) *
                     (m.cosine * std::cos(m.frequency * t) + m.sine * std::sin(m.frequency * t));
        k[static_cast<std::size_t>(n)] = value;
    }

    const StateSpaceFit fit = fitStateSpace(k, dt, 4);
    expectTrue("the fit converged", fit.converged);
    expectEqual("four states for two complex pairs", fit.model.stateCount(), 4);
    expectTrue("the model is stable", fit.model.stable());
    expectTrue("a synthetic sum of damped sinusoids is recovered to 1e-6 relative",
               fit.relativeRms < 1e-6);

    // The poles themselves, not just the curve through them.
    for (const Mode& m : wanted) {
        double best = 1e30;
        for (const RadiationStateSpace::Mode& got : fit.model.modes)
            best = std::min(best, std::hypot(got.decay - m.decay, got.frequency - m.frequency));
        expectTrue("each planted pole is recovered", best < 1e-6);
    }
    // Vacuity guard: the signal has to be something a constant could not match.
    expectTrue("the synthetic K is not trivially small at late times",
               std::abs(k[0]) > 0 && std::abs(k[static_cast<std::size_t>(count / 4)]) > 0);
    std::printf("     Prony on a planted 4-pole signal: relative rms %.2e\n", fit.relativeRms);
}

void testStateSpaceApproximatesTheHull() {
    const RadiationTable& table = ferryTable();
    const double dt = 0.1;
    const int count = 400;   // 40 s
    const std::vector<double> k = retardationFunction(table, 2, 2, dt, count);

    double previous = 1e30;
    for (int order : {4, 6, 8}) {
        const StateSpaceFit fit = fitStateSpace(k, dt, order);
        expectTrue("the hull fit converged", fit.converged);
        expectTrue("the hull fit is stable", fit.model.stable());
        expectTrue("more states do not make the fit worse", fit.relativeRms < previous * 1.2);
        previous = fit.relativeRms;
        if (order == 6) {
            expectTrue("six states reproduce K to better than 10% relative rms",
                       fit.relativeRms < 0.10);
            expectTrue("and to better than 5% of K(0) anywhere",
                       fit.peakError < 0.05 * std::abs(k[0]));
            std::printf("     state space on the ferry's K33: order 6, %d states,"
                        " relative rms %.2e, peak error %.2f%% of K(0)\n",
                        fit.model.stateCount(), fit.relativeRms,
                        100.0 * fit.peakError / std::abs(k[0]));
        }
        // The impulse response is the thing that actually gets used, so check it
        // directly rather than trusting the residual that produced it.
        double worst = 0;
        for (int n = 0; n < count; ++n)
            worst = std::max(worst, std::abs(fit.model.impulseResponse(n * dt) -
                                             k[static_cast<std::size_t>(n)]));
        expectNear("impulseResponse() agrees with the fit's own residual", worst, fit.peakError,
                   1e-9 * std::abs(k[0]));
    }
}

// The evaluator has to reproduce K when driven with an impulse: that is what a
// convolution *is*, and it exercises the zero-order-hold update, the state
// layout and the output row together.
void testRadiationForceReproducesTheConvolution() {
    const RadiationTable& table = ferryTable();
    const double dt = 0.1;
    const int count = 400;
    RadiationForce force(table, dt, count, 6);
    expectTrue("some entries got state-space models", force.modelCount() > 0);
    expectTrue("A_inf is passed through unchanged",
               force.addedMassInfinite()[2][2] == table.addedMassInfinite[2][2]);

    const std::vector<double> k = retardationFunction(table, 2, 2, dt, count);
    const StateSpaceFit fit = fitStateSpace(k, dt, 6);

    // A velocity of 1/dt held for one step is a unit impulse in displacement, so
    // the memory force afterwards traces the impulse response.
    force.reset();
    std::array<double, 6> velocity{};
    velocity[2] = 1.0 / dt;
    force.step(velocity, dt);
    velocity[2] = 0.0;
    // Zero-order hold: one step at 1/dt is a rectangle of unit area, not a Dirac
    // impulse, so what comes back is K averaged over the step. Comparing against
    // the midpoint is the second-order-accurate way to say that, and it is what
    // makes this a check on the update rather than on the sampling.
    double worst = 0;
    for (int n = 1; n < 200; ++n) {
        worst = std::max(worst, std::abs(force.memoryForce()[2] -
                                         fit.model.impulseResponse((n - 0.5) * dt)));
        force.step(velocity, dt);
    }
    expectTrue("the evaluator reproduces its own model's impulse response",
               worst < 0.01 * std::abs(k[0]));

    // Against the real K, which is the number that matters.
    force.reset();
    velocity[2] = 1.0 / dt;
    force.step(velocity, dt);
    velocity[2] = 0.0;
    double worstAgainstK = 0;
    for (int n = 1; n < 200; ++n) {
        worstAgainstK = std::max(
            worstAgainstK,
            std::abs(force.memoryForce()[2] -
                     0.5 * (k[static_cast<std::size_t>(n)] + k[static_cast<std::size_t>(n - 1)])));
        force.step(velocity, dt);
    }
    expectTrue("and tracks the true retardation function to 5% of K(0)",
               worstAgainstK < 0.05 * std::abs(k[0]));

    // Still water, still ship: no memory, no force. And the memory must not
    // survive a reset.
    force.reset();
    const std::array<double, 6> zero{};
    for (int n = 0; n < 50; ++n) force.step(zero, dt);
    for (int d = 0; d < 6; ++d)
        expectNear("a ship that has never moved feels no radiation force",
                   force.memoryForce()[static_cast<std::size_t>(d)], 0.0, 0.0);

    // A steady velocity must produce a *bounded* force -- the integral of K is
    // finite because K decays -- and it must settle rather than grow.
    force.reset();
    std::array<double, 6> steady{};
    steady[2] = 1.0;
    for (int n = 0; n < 300; ++n) force.step(steady, dt);
    const double settled = force.memoryForce()[2];
    for (int n = 0; n < 300; ++n) force.step(steady, dt);
    const double later = force.memoryForce()[2];
    expectTrue("a steady heave velocity gives a bounded memory force",
               std::isfinite(settled) && std::abs(settled) > 0);
    expectNear("and it has settled", later, settled, 0.02 * std::abs(settled));
    // The settled value is the integral of K, and that integral is B(omega -> 0),
    // which radiation damping makes exactly zero: a hull towed at constant
    // vertical speed radiates nothing in the limit. So K has zero mean, and the
    // settled force must be small against the scale K sets -- not against
    // itself, which is what makes this a statement about the physics rather than
    // about two truncations agreeing.
    double absoluteArea = 0;
    for (std::size_t n = 0; n + 1 < k.size(); ++n)
        absoluteArea += 0.5 * dt * (std::abs(k[n]) + std::abs(k[n + 1]));
    expectTrue("the area under |K| is a real scale", absoluteArea > 0);
    expectTrue("K integrates to nearly zero, because B(omega -> 0) is zero",
               std::abs(settled) < 0.2 * absoluteArea);
    std::printf("     evaluator: %d models, impulse response within %.2f%% of K(0);"
                " steady-state force %.3e is %.1f%% of the area under |K|\n",
                force.modelCount(), 100.0 * worstAgainstK / std::abs(k[0]), settled,
                100.0 * std::abs(settled) / absoluteArea);
}

// --- Validity ----------------------------------------------------------------

void testValidityIsReported() {
    const RadiationHull hull = ferry();
    const std::vector<double> sane = radiationFrequencyGrid(0.2, 1.5, 10);
    const std::vector<std::string> clean = validateRadiationHull(hull, sane);
    for (const std::string& s : clean) std::printf("     unexpected: %s\n", s.c_str());
    expectTrue("the reference ferry inside a sane band is not flagged", clean.empty());

    expectTrue("a band whose shortest wave is under a beam is flagged",
               !validateRadiationHull(hull, radiationFrequencyGrid(0.2, 6.0, 10)).empty());

    RadiationHull stubby = hull;
    for (RadiationStation& s : stubby.stations) s.x *= 0.05;   // B/L far past slender
    expectTrue("a hull too beamy for strip theory is flagged",
               !validateRadiationHull(stubby, sane).empty());

    RadiationHull coarse = hull;
    coarse.panelsPerHalfSection = 4;
    expectTrue("too few panels is flagged", !validateRadiationHull(coarse, sane).empty());

    RadiationHull empty;
    expectTrue("an empty hull is rejected outright", !validateRadiationHull(empty, sane).empty());
    expectEqual("and produces an empty table rather than a NaN",
                stripTheoryTable(empty, sane).size(), 0);

    RadiationHull offset = hull;
    for (RadiationStation& s : offset.stations) s.x += 200.0;   // all forward of midship
    expectTrue("stations that do not straddle midship are flagged",
               !validateRadiationHull(offset, sane).empty());

    // Degenerate inputs must return zeros, not NaNs.
    const SectionCoefficients nothing = sectionCoefficients(lewisSection(0, 0, 0), 1.0, 1000.0, 20);
    expectNear("a degenerate section produces no added mass", nothing.addedMass[1][1], 0.0, 0.0);
    const std::vector<double> emptyK = retardationFunction(RadiationTable{}, 2, 2, 0.1, 10);
    expectEqual("an empty table produces no retardation function",
                static_cast<long long>(emptyK.size()), 0);
    expectNear("memoryDecayTime of nothing is zero", memoryDecayTime({}, 0.1), 0.0, 0.0);
    expectTrue("fitting nothing does not converge", !fitStateSpace({}, 0.1, 4).converged);
}

}  // namespace

void runRadiationTests() {
    std::printf("\n--- radiation (strip theory, Cummins/Ogilvie) ---\n");
    testLewisFormsReproduceTheirInputs();
    testHeavingCircleAgainstUrsell();
    testCircleRollingAboutItsCentreDisplacesNothing();
    testPrismaticHullIntegratesAlongItsLength();
    testGeometricScalingIsExact();
    testOgilvieTransformsAgainstClosedForms();
    testRetardationRoundTripsThroughOgilvie();
    testDampingIsPhysical();
    testStateSpaceRecoversKnownPoles();
    testStateSpaceApproximatesTheHull();
    testRadiationForceReproducesTheConvolution();
    testValidityIsReported();
}
