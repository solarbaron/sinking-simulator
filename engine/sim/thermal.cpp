// SPDX-License-Identifier: MIT
#include "thermal.hpp"

#include "buckling.hpp"   // johnsonOstenfeld, for the restrained-buckling temperature
#include "reduction.hpp"  // bandwidthReducingOrder

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <limits>
#include <map>

namespace sim::thermal {
namespace {

// The same node ordering `solid_shell.hpp` fixes. Repeated here rather than
// reached for, because it is a *contract* stated in that header and not a private
// implementation detail -- and the conduction element does not care about it in
// the way the mechanical one does, so a shared symbol would suggest a coupling
// that is not there.
constexpr double kXi[kNodes]  = {-1, +1, +1, -1, -1, +1, +1, -1};
constexpr double kEta[kNodes] = {-1, -1, +1, +1, -1, -1, +1, +1};
constexpr double kZta[kNodes] = {-1, -1, -1, -1, +1, +1, +1, +1};

// The six faces, each four nodes. Windings are *not* consistently outward -- face
// 0 winds counter-clockwise seen from +zeta and so points into the element -- and
// `boundaryFaces` fixes that against the element centroid rather than by
// transcribing a corrected table, which is where this sort of code acquires a
// flipped face that only shows on one of the six.
constexpr int kFaces[6][4] = {{0, 1, 2, 3}, {4, 5, 6, 7}, {0, 1, 5, 4},
                              {3, 2, 6, 7}, {0, 3, 7, 4}, {1, 2, 6, 5}};

constexpr double kCorner[4][2] = {{-1, -1}, {+1, -1}, {+1, +1}, {-1, +1}};

double determinant3(const double m[3][3]) {
    return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
           m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
           m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
}

double invert3(const double m[3][3], double inv[3][3]) {
    const double det = determinant3(m);
    if (det == 0.0) return 0.0;
    const double d = 1.0 / det;
    inv[0][0] = (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * d;
    inv[0][1] = (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * d;
    inv[0][2] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * d;
    inv[1][0] = (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * d;
    inv[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * d;
    inv[1][2] = (m[0][2] * m[1][0] - m[0][0] * m[1][2]) * d;
    inv[2][0] = (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * d;
    inv[2][1] = (m[0][1] * m[2][0] - m[0][0] * m[2][1]) * d;
    inv[2][2] = (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * d;
    return det;
}

// EN 1993-1-2 is written in Celsius; every curve below converts once, here.
double celsius(double kelvin) { return kelvin - kCelsius; }

}  // namespace

// --- Material curves -------------------------------------------------------------

double carbonSteelConductivity(double kelvin) {
    const double t = std::clamp(celsius(kelvin), 20.0, 1200.0);
    // 3.4.1.3. The standard's own two pieces do not quite meet -- the linear one
    // reaches 27.36 at 800 C against the flat 27.3 above it, a 0.2% step. It is
    // reproduced as written rather than smoothed, because a value smoothed here
    // would no longer be the one a fire engineer's hand calculation uses.
    if (t >= 800.0) return 27.3;
    return 54.0 - 3.33e-2 * t;
}

double carbonSteelSpecificHeat(double kelvin) {
    const double t = std::clamp(celsius(kelvin), 20.0, 1200.0);
    // 3.4.1.2. The two hyperbolic pieces both reach 5000 J/(kg K) at 735 C, which
    // is the ferrite-austenite transition standing in for a latent heat.
    if (t < 600.0) return 425.0 + 0.773 * t - 1.69e-3 * t * t + 2.22e-6 * t * t * t;
    if (t < 735.0) return 666.0 + 13002.0 / (738.0 - t);
    if (t < 900.0) return 545.0 + 17820.0 / (t - 731.0);
    return 650.0;
}

double carbonSteelEnthalpy(double kelvin) {
    // The exact integral of the curve above, piece by piece, with the constants
    // chosen so `h` is continuous. `h(20 C) = 0`.
    //
    // Continuity has to be imposed rather than inherited: the standard's pieces
    // are themselves discontinuous by about 0.3 J/(kg K) at 600 C and 0.4 at
    // 900 C, so integrating each piece from its own lower limit and adding would
    // leave `h` with steps in it and the energy account with a residual that grew
    // with the number of times a plate crossed 600 C.
    const auto polynomial = [](double t) {
        return 425.0 * t + 0.3865 * t * t - (1.69e-3 / 3.0) * t * t * t +
               (2.22e-6 / 4.0) * t * t * t * t;
    };
    const auto first = [](double t) { return 666.0 * t - 13002.0 * std::log(738.0 - t); };
    const auto second = [](double t) { return 545.0 * t + 17820.0 * std::log(t - 731.0); };

    // Junction values, evaluated once so the four pieces cannot drift apart.
    static const double kAt600 = polynomial(600.0) - polynomial(20.0);
    static const double kAt735 = kAt600 + first(735.0) - first(600.0);
    static const double kAt900 = kAt735 + second(900.0) - second(735.0);

    const double raw = celsius(kelvin);
    const double t = std::clamp(raw, 20.0, 1200.0);
    double h;
    if (t < 600.0) {
        h = polynomial(t) - polynomial(20.0);
    } else if (t < 735.0) {
        h = kAt600 + first(t) - first(600.0);
    } else if (t < 900.0) {
        h = kAt735 + second(t) - second(735.0);
    } else {
        h = kAt900 + 650.0 * (t - 900.0);
    }
    // Outside the standard's range `c` is clamped, so `h` continues linearly at
    // the end value. Anything else would make `h` non-monotone and the secant
    // capacity negative, which turns an unconditionally stable scheme into one
    // that is not.
    if (raw < 20.0) h += carbonSteelSpecificHeat(kelvin) * (raw - 20.0);
    if (raw > 1200.0) h += 650.0 * (raw - 1200.0);
    return h;
}

// --- Strength at temperature -------------------------------------------------------

namespace {

// EN 1993-1-2:2005 Table 3.1, transcribed. Thirteen stations at hundred-degree
// intervals from 20 C; the standard prescribes linear interpolation between them,
// so this is the whole model and not a fit to it.
//
// It is *not* exported. A test that asserted the interpolant against this array
// would be asserting the code against itself -- move a row and both sides move
// together -- so `tests/test_thermal.cpp` carries its own transcription of the
// published table and compares against that. Two independent transcriptions is
// the point.
constexpr int kReductionStations = 13;
constexpr double kReductionCelsius[kReductionStations] = {20,  100, 200, 300,  400,  500, 600,
                                                          700, 800, 900, 1000, 1100, 1200};
constexpr double kYieldFactor[kReductionStations] = {1.000, 1.000, 1.000, 1.000, 1.000,
                                                     0.780, 0.470, 0.230, 0.110, 0.060,
                                                     0.040, 0.020, 0.000};
constexpr double kProportionalFactor[kReductionStations] = {1.0000, 1.0000, 0.8070, 0.6130,
                                                            0.4200, 0.3600, 0.1800, 0.0750,
                                                            0.0500, 0.0375, 0.0250, 0.0125,
                                                            0.0000};
constexpr double kModulusFactor[kReductionStations] = {1.0000, 1.0000, 0.9000, 0.8000,
                                                       0.7000, 0.6000, 0.3100, 0.1300,
                                                       0.0900, 0.0675, 0.0450, 0.0225,
                                                       0.0000};

// Linear interpolation on the station list, clamped at both ends.
//
// **The tabulated points come back exactly.** The `t == station` early exit is
// what says so, and it is worth being honest about what it currently buys:
// **nothing, on this table.** `s` is exactly 1 at the upper station and exactly 0
// at the lower, so the general expression is `a + (b - a) * 1` and `a + (b - a) *
// 0`; the second is `a` for every finite `a` by IEEE, and the first happens to be
// `b` for **all thirty-six adjacent pairs of Table 3.1**. That was assumed here
// first -- with a worked example claiming 0.0675 and 0.045 were a unit in the last
// place apart, which they are not -- and mutation testing found the branch
// unobservable and the claim invented.
//
// It stays, because `a + (b - a) == b` is a property of *these numbers* and not of
// the method, and a table row added later has no reason to keep it.
// `tests/test_thermal.cpp` now asserts the property over every pair directly, so
// the day it stops holding is caught there rather than in a caller.
double interpolate(const double value[kReductionStations], double t) {
    if (t <= kReductionCelsius[0]) return value[0];
    if (t >= kReductionCelsius[kReductionStations - 1]) return value[kReductionStations - 1];
    for (int i = 1; i < kReductionStations; ++i) {
        if (t > kReductionCelsius[i]) continue;
        if (t == kReductionCelsius[i]) return value[i];
        const double lo = kReductionCelsius[i - 1], hi = kReductionCelsius[i];
        const double s = (t - lo) / (hi - lo);
        return value[i - 1] + (value[i] - value[i - 1]) * s;
    }
    return value[kReductionStations - 1];
}

}  // namespace

SteelReduction carbonSteelReduction(double kelvin) {
    const double t = celsius(kelvin);
    SteelReduction out;
    out.effectiveYield = interpolate(kYieldFactor, t);
    out.proportionalLimit = interpolate(kProportionalFactor, t);
    out.youngsModulus = interpolate(kModulusFactor, t);
    return out;
}

double carbonSteelYieldFactor(double kelvin) {
    return interpolate(kYieldFactor, celsius(kelvin));
}
double carbonSteelProportionalFactor(double kelvin) {
    return interpolate(kProportionalFactor, celsius(kelvin));
}
double carbonSteelModulusFactor(double kelvin) {
    return interpolate(kModulusFactor, celsius(kelvin));
}

double carbonSteelStress(double strain, double kelvin, double yieldStrength,
                         double youngsModulus) {
    // The standard writes the curve in tension and takes compression as its
    // mirror. Doing the reflection here rather than in four branches keeps the
    // ellipse's construction in one place.
    const double sign = strain < 0.0 ? -1.0 : 1.0;
    const double eps = std::abs(strain);

    const SteelReduction k = carbonSteelReduction(kelvin);
    const double fy = k.effectiveYield * yieldStrength;
    const double fp = k.proportionalLimit * yieldStrength;
    const double ea = k.youngsModulus * youngsModulus;

    // The four strain landmarks. Only the first depends on temperature; 2%, 15%
    // and 20% are fixed by the standard.
    const double epsY = 0.02, epsT = 0.15, epsU = 0.20;
    if (!(ea > 0.0)) return 0.0;             // 1200 C: no stiffness and no strength
    const double epsP = fp / ea;

    if (eps <= epsP) return sign * ea * eps;
    if (eps >= epsU) return 0.0;
    if (eps >= epsT) return sign * fy * (1.0 - (eps - epsT) / (epsU - epsT));
    if (eps >= epsY) return sign * fy;

    // The elliptical transition. `c`, `a` and `b` are the standard's own symbols;
    // the construction makes the ellipse tangent to the elastic line at `epsP` and
    // horizontal at `epsY`, which is why it is an ellipse and not a spline.
    //
    // Degenerate when `fy == fp` -- every temperature at or below 100 C, where both
    // factors are 1 -- and there the curve is elastic-perfectly-plastic with no
    // transition at all. `c` is then zero over a zero denominator, so the case is
    // taken rather than divided.
    const double gap = fy - fp;
    const double span = epsY - epsP;
    if (!(gap > 0.0) || !(span > 0.0)) return sign * fy;

    const double denominator = span * ea - 2.0 * gap;
    if (!(denominator > 0.0)) return sign * fy;
    const double c = gap * gap / denominator;
    const double aSquared = span * (span + c / ea);
    const double bSquared = c * span * ea + c * c;
    const double reach = epsY - eps;
    const double inside = aSquared - reach * reach;
    if (!(inside > 0.0) || !(aSquared > 0.0)) return sign * fy;
    return sign * (fp - c + std::sqrt(bSquared / aSquared) * std::sqrt(inside));
}

// --- Thermal elongation -------------------------------------------------------------

double carbonSteelElongation(double kelvin) {
    const double t = celsius(kelvin);
    // Below 20 C the standard says nothing. Clamping to zero rather than running
    // the polynomial backwards is the same choice the conductivity makes, and it
    // has a consequence worth stating: a structure *cooled* below its reference
    // temperature contracts and goes into tension, which this does not model. A
    // fire is not that problem and inventing an extrapolation for it here would
    // put an unvalidated branch on the path of every heated element.
    if (t <= 20.0) return 0.0;
    if (t < 750.0) return 1.2e-5 * t + 0.4e-8 * t * t - 2.416e-4;
    if (t <= 860.0) return 1.1e-2;
    if (t <= 1200.0) return 2.0e-5 * t - 6.2e-3;
    return 2.0e-5 * 1200.0 - 6.2e-3;
}

double thermalStrain(double kelvin, double referenceKelvin) {
    return carbonSteelElongation(kelvin) - carbonSteelElongation(referenceKelvin);
}

void thermalEigenstrain(double kelvin, double referenceKelvin, double out[6]) {
    const double e = thermalStrain(kelvin, referenceKelvin);
    out[0] = out[1] = out[2] = e;
    out[3] = out[4] = out[5] = 0.0;
}

double restrainedStress(const StructuralMaterial& material, double kelvin,
                        double referenceKelvin) {
    return -carbonSteelModulusFactor(kelvin) * material.youngsModulus *
           thermalStrain(kelvin, referenceKelvin);
}

namespace {

// First temperature above `referenceKelvin` at which `margin` changes sign, by a
// scan and then a bisection. Kelvin, or zero when it never does.
//
// A scan is used rather than a root find from one end because the margin is
// **piecewise** in temperature -- kinked at every hundred-degree station of
// Table 3.1, and genuinely discontinuous at the 750 C step in the elongation and
// at the Johnson-Ostenfeld transition -- so a method that assumes smoothness would
// converge on whichever root it happened to bracket. The half-kelvin grid is
// finer than any feature of either curve; the tabulated stations are 100 K apart
// and the two discontinuities are single points, so no crossing can hide between
// two samples of a function that is monotone between kinks.
template <typename Margin>
double firstCrossing(const Margin& margin, double referenceKelvin) {
    const double top = kCelsius + 1200.0;
    if (!(referenceKelvin < top)) return 0.0;
    double lo = referenceKelvin, fLo = margin(lo);
    if (fLo >= 0.0) return lo;
    for (double t = referenceKelvin + 0.5; t <= top + 0.25; t += 0.5) {
        const double hi = std::min(t, top);
        const double fHi = margin(hi);
        if (fHi >= 0.0) {
            // 200 bisections is far past what a double can distinguish over this
            // range; it costs microseconds and removes the question of whether the
            // count was tuned to the answer.
            double a = lo, b = hi;
            for (int i = 0; i < 200; ++i) {
                const double m = 0.5 * (a + b);
                if (margin(m) >= 0.0) b = m;
                else a = m;
            }
            return 0.5 * (a + b);
        }
        lo = hi;
        fLo = fHi;
        if (hi >= top) break;
    }
    return 0.0;
}

}  // namespace

double restrainedYieldTemperature(const StructuralMaterial& material, double referenceKelvin) {
    return firstCrossing(
        [&](double k) {
            return std::abs(restrainedStress(material, k, referenceKelvin)) -
                   carbonSteelYieldFactor(k) * material.yieldStrength;
        },
        referenceKelvin);
}

double restrainedBucklingTemperature(double elasticStress, const StructuralMaterial& material,
                                     double referenceKelvin) {
    if (!(elasticStress > 0.0)) return 0.0;
    return firstCrossing(
        [&](double k) {
            const SteelReduction r = carbonSteelReduction(k);
            return std::abs(restrainedStress(material, k, referenceKelvin)) -
                   johnsonOstenfeld(r.youngsModulus * elasticStress,
                                    r.effectiveYield * material.yieldStrength);
        },
        referenceKelvin);
}

double twoStripStress(const StructuralMaterial& material, double kelvin, double hotFraction,
                      double referenceKelvin) {
    const double f = std::clamp(hotFraction, 0.0, 1.0);
    if (f >= 1.0) return 0.0;   // nothing cold left to hold it: free expansion
    // `-k_E E eps / (1 + (f/(1-f)) k_E)`, and not the reciprocal-sum form the header
    // quotes it in, so that f = 0 divides by a literal 1.0 and returns
    // `restrainedStress` bit for bit.
    return restrainedStress(material, kelvin, referenceKelvin) /
           (1.0 + (f / (1.0 - f)) * carbonSteelModulusFactor(kelvin));
}

double temperatureForElongation(double elongation) {
    const double lo = kCelsius + 20.0, hi = kCelsius + 1200.0;
    if (!(elongation > 0.0)) return lo;
    if (elongation >= carbonSteelElongation(hi)) return hi;
    // The curve is monotone non-decreasing, so a bisection on "has it got there yet"
    // converges on the *lowest* temperature that has -- which is the bottom of the
    // 750-860 C plateau for any elongation the plateau holds, and the exact answer
    // everywhere else. 200 halvings is past what a double distinguishes over 1180 K,
    // the same count and the same reasoning as `firstCrossing`.
    double a = lo, b = hi;
    for (int i = 0; i < 200; ++i) {
        const double m = 0.5 * (a + b);
        if (carbonSteelElongation(m) >= elongation) b = m;
        else a = m;
    }
    return 0.5 * (a + b);
}

double beamColumnMagnifier(double axialOverEuler) {
    if (!(axialOverEuler > 0.0)) return 1.0;
    if (axialOverEuler >= 1.0) return std::numeric_limits<double>::infinity();
    // `2 (sec u - 1) / u^2` is 0/0 at the origin and loses most of its digits to
    // cancellation anywhere near it. The identities `sec u - 1 = (1 - cos u)/cos u`
    // and `1 - cos u = 2 sin^2(u/2)` turn it into
    //
    //     Psi = (sin(u/2) / (u/2))^2 / cos u
    //
    // which is the same number exactly and has no subtraction in it at all: the
    // sinc factor is 1.0 to the bit for small u, so the limit is approached rather
    // than computed. No series expansion and therefore no truncation term to size.
    const double u = 0.5 * kPi * std::sqrt(axialOverEuler);
    const double sinc = std::sin(0.5 * u) / (0.5 * u);
    return sinc * sinc / std::cos(u);
}

MemberState memberState(const HeatedMember& member, const StructuralMaterial& material,
                        double kelvin, double referenceKelvin) {
    const SteelReduction r = carbonSteelReduction(kelvin);
    const double restraint = std::max(member.restraint, 0.0);

    MemberState s;
    s.axialStress = std::abs(restraint * restrainedStress(material, kelvin, referenceKelvin));
    s.eulerStress = r.youngsModulus * std::max(member.eulerStress, 0.0);
    s.yieldCapacity = r.effectiveYield * material.yieldStrength;
    s.columnCapacity = johnsonOstenfeld(s.eulerStress, s.yieldCapacity);

    // A member with no stiffness left has already gone, whatever it is carrying: at
    // 1200 C EN 1993-1-2 takes both factors to exactly zero, which is the case
    // `restrainedBucklingTemperature` records as always finding a crossing.
    const double infinite = std::numeric_limits<double>::infinity();
    s.axialOverEuler = s.eulerStress > 0.0 ? s.axialStress / s.eulerStress
                                           : (s.axialStress > 0.0 ? infinite : 0.0);
    s.magnifier = beamColumnMagnifier(s.axialOverEuler);

    // First-order bending stress. Guarded against `0 * infinity`: a member carrying
    // no lateral load at all carries none however far past its Euler load it is,
    // and that case is the one the axial identity below rests on.
    const double first = member.modulus > 0.0 ? std::abs(member.lateralMoment) / member.modulus
                                              : (member.lateralMoment != 0.0 ? infinite : 0.0);
    s.bendingStress = first > 0.0 ? first * s.magnifier : 0.0;

    const double axialTerm = s.columnCapacity > 0.0
                                 ? s.axialStress / s.columnCapacity
                                 : (s.axialStress > 0.0 ? infinite : 0.0);
    const auto bend = [&](double stress) {
        return s.yieldCapacity > 0.0 ? stress / s.yieldCapacity : (stress > 0.0 ? infinite : 0.0);
    };
    s.utilisation = axialTerm + bend(s.bendingStress);
    s.additiveUtilisation = axialTerm + bend(first);
    if (s.utilisation >= 1.0)
        s.limit = axialTerm >= 1.0 ? MemberLimit::Column : MemberLimit::Interaction;
    return s;
}

double memberFailureTemperature(const HeatedMember& member, const StructuralMaterial& material,
                               double referenceKelvin) {
    return firstCrossing(
        [&](double k) { return memberState(member, material, k, referenceKelvin).utilisation - 1.0; },
        referenceKelvin);
}

StructuralMaterial atTemperature(const StructuralMaterial& material, double kelvin) {
    const SteelReduction k = carbonSteelReduction(kelvin);
    StructuralMaterial out = material;
    out.youngsModulus *= k.youngsModulus;
    out.yieldStrength *= k.effectiveYield;
    return out;
}

plasticity::Material atTemperature(const plasticity::Material& material, double kelvin) {
    const SteelReduction k = carbonSteelReduction(kelvin);
    plasticity::Material out = material;
    out.youngsModulus *= k.youngsModulus;
    // The whole curve, so that `flowStress` scales by `k_y` at every plastic
    // strain and Considere's root does not move. Which of these are live depends
    // on `flow.kind`, and scaling all of them is what keeps a curve that is
    // switched from Linear to Swift after being reduced still reduced.
    out.flow.yieldStrength *= k.effectiveYield;
    out.flow.hardeningModulus *= k.effectiveYield;
    out.flow.strengthCoefficient *= k.effectiveYield;
    out.flow.kinematicModulus *= k.effectiveYield;
    return out;
}

// --- Element forms ---------------------------------------------------------------

bool computeForms(const double nodes[kDof], Forms& out) {
    out.ok = true;
    const double q = 1.0 / std::sqrt(3.0);
    for (int gp = 0; gp < kGauss; ++gp) {
        const double xi = (gp & 1) ? q : -q;
        const double eta = (gp & 2) ? q : -q;
        const double zta = (gp & 4) ? q : -q;

        double dN[kNodes][3];
        for (int a = 0; a < kNodes; ++a) {
            const double x = 1.0 + xi * kXi[a];
            const double y = 1.0 + eta * kEta[a];
            const double z = 1.0 + zta * kZta[a];
            out.shape[gp][a] = 0.125 * x * y * z;
            dN[a][0] = 0.125 * kXi[a] * y * z;
            dN[a][1] = 0.125 * kEta[a] * x * z;
            dN[a][2] = 0.125 * kZta[a] * x * y;
        }

        // jac[i][k] = dx_i / dxi_k, the same layout solid_shell.cpp uses.
        double jac[3][3];
        for (int i = 0; i < 3; ++i)
            for (int k = 0; k < 3; ++k) {
                double s = 0.0;
                for (int a = 0; a < kNodes; ++a)
                    s += dN[a][k] * nodes[static_cast<std::size_t>(a) * 3 + static_cast<std::size_t>(i)];
                jac[i][k] = s;
            }
        double inv[3][3];
        const double det = invert3(jac, inv);
        if (!(det > 0.0)) {
            out.ok = false;
            return false;
        }
        out.weight[gp] = det;  // 2x2x2 Gauss weights are all one

        // dN_a/dx_i = dN_a/dxi_k * dxi_k/dx_i, and inv[k][i] is dxi_k/dx_i.
        for (int a = 0; a < kNodes; ++a)
            for (int i = 0; i < 3; ++i) {
                double s = 0.0;
                for (int k = 0; k < 3; ++k) s += dN[a][k] * inv[k][i];
                out.gradient[gp][a][i] = s;
            }
    }
    return true;
}

void gaussTemperature(const double nodal[kNodes], double out[kGauss]) {
    const double q = 1.0 / std::sqrt(3.0);
    for (int gp = 0; gp < kGauss; ++gp) {
        // The same bit pattern `computeForms` uses, which is the same one
        // `solidshell` uses: xi from bit 0, eta from bit 1, zeta from bit 2.
        const double xi = (gp & 1) ? q : -q;
        const double eta = (gp & 2) ? q : -q;
        const double zta = (gp & 4) ? q : -q;
        double sum = 0.0;
        for (int a = 0; a < kNodes; ++a)
            sum += 0.125 * (1.0 + xi * kXi[a]) * (1.0 + eta * kEta[a]) * (1.0 + zta * kZta[a]) *
                   nodal[a];
        out[gp] = sum;
    }
}

double elementTemperature(const Forms& forms, const double nodal[kNodes]) {
    if (!forms.ok) return 0.0;
    double gauss[kGauss];
    gaussTemperature(nodal, gauss);
    double weighted = 0.0, volume = 0.0;
    for (int gp = 0; gp < kGauss; ++gp) {
        weighted += forms.weight[gp] * gauss[gp];
        volume += forms.weight[gp];
    }
    return volume > 0.0 ? weighted / volume : 0.0;
}

bool elementTemperatures(const solidshell::HexMesh& mesh, const std::vector<double>& nodal,
                         std::vector<double>& out) {
    out.clear();
    if (nodal.size() != mesh.nodeCount()) return false;
    out.resize(mesh.elementCount(), 0.0);
    bool ok = true;
    for (std::size_t e = 0; e < mesh.elementCount(); ++e) {
        double nodes[kDof], value[kNodes];
        for (int a = 0; a < kNodes; ++a) {
            const std::size_t n = mesh.index[e * kNodes + static_cast<std::size_t>(a)];
            value[a] = nodal[n];
            for (int i = 0; i < 3; ++i) nodes[a * 3 + i] = mesh.position[n * 3 + static_cast<std::size_t>(i)];
        }
        Forms forms;
        if (!computeForms(nodes, forms)) {
            // The `continue` is belt-and-braces: `elementTemperature` refuses a
            // `Forms` that did not build and returns zero, which is the fill value
            // this element already has. Mutation testing scored deleting it as an
            // equivalent mutant for exactly that reason, and it is kept because
            // two readers of "this element has no answer" agreeing is the point.
            ok = false;
            continue;
        }
        out[e] = elementTemperature(forms, value);
    }
    return ok;
}

void conductance(const Forms& forms, const double conductivity[kGauss],
                 double out[kNodes * kNodes]) {
    std::fill(out, out + kNodes * kNodes, 0.0);
    if (!forms.ok) return;
    for (int gp = 0; gp < kGauss; ++gp) {
        const double w = forms.weight[gp] * conductivity[gp];
        for (int a = 0; a < kNodes; ++a)
            for (int b = 0; b < kNodes; ++b) {
                double s = 0.0;
                for (int i = 0; i < 3; ++i) s += forms.gradient[gp][a][i] * forms.gradient[gp][b][i];
                out[a * kNodes + b] += w * s;
            }
    }
}

void capacity(const Forms& forms, const double volumetricCapacity[kGauss], bool lumped,
              double out[kNodes * kNodes]) {
    std::fill(out, out + kNodes * kNodes, 0.0);
    if (!forms.ok) return;
    for (int gp = 0; gp < kGauss; ++gp) {
        const double w = forms.weight[gp] * volumetricCapacity[gp];
        for (int a = 0; a < kNodes; ++a)
            for (int b = 0; b < kNodes; ++b)
                out[a * kNodes + b] += w * forms.shape[gp][a] * forms.shape[gp][b];
    }
    if (!lumped) return;
    double row[kNodes];
    for (int a = 0; a < kNodes; ++a) {
        double s = 0.0;
        for (int b = 0; b < kNodes; ++b) s += out[a * kNodes + b];
        row[a] = s;
    }
    std::fill(out, out + kNodes * kNodes, 0.0);
    for (int a = 0; a < kNodes; ++a) out[a * kNodes + a] = row[a];
}

// --- Explicit stability, for the record -------------------------------------------

double explicitLimit(const double nodes[kDof], const StructuralMaterial& material) {
    Forms forms;
    if (!computeForms(nodes, forms)) return 0.0;
    double k[kGauss], rc[kGauss];
    for (int gp = 0; gp < kGauss; ++gp) {
        k[gp] = material.conductivity;
        rc[gp] = material.density * material.specificHeat;
    }
    double ke[kNodes * kNodes], ce[kNodes * kNodes];
    conductance(forms, k, ke);
    capacity(forms, rc, /*lumped=*/true, ce);

    // Largest eigenvalue of `K x = lambda C x`; the forward-Euler limit is
    // `2/lambda_max`.
    //
    // **By a full eigensolve and not by power iteration**, and that is a
    // correction rather than a preference. Power iteration from a checkerboard
    // start is the obvious implementation and it is *wrong on exactly the elements
    // this file exists for*: the highest mode of a plate element is the jump
    // through its thin direction, and a checkerboard start -- with or without a
    // linear perturbation -- is orthogonal to that mode by symmetry, so the
    // iteration converges to a lower eigenvalue and reports a stable step that is
    // three times too long. It was caught by asserting the result against the
    // element's own `rho c h^2 / 2k`, which is the only reason it is not still
    // here. `reduction::generalisedEigen` has no start vector to get wrong, and a
    // third eigensolver in this repository would be a third place to be wrong.
    std::vector<double> a(kNodes * kNodes), b(kNodes * kNodes);
    for (int i = 0; i < kNodes * kNodes; ++i) {
        a[static_cast<std::size_t>(i)] = ke[i];
        b[static_cast<std::size_t>(i)] = ce[i];
    }
    const reduction::Eigenpairs pairs = reduction::generalisedEigen(a, b, kNodes);
    if (!pairs.converged || pairs.count == 0) return 0.0;
    const double lambda = pairs.value[static_cast<std::size_t>(pairs.count - 1)];
    return lambda > 0.0 ? 2.0 / lambda : 0.0;
}

double explicitLimit(const solidshell::HexMesh& mesh, const StructuralMaterial& material) {
    double smallest = 0.0;
    bool any = false;
    for (std::size_t e = 0; e < mesh.elementCount(); ++e) {
        double nodes[kDof];
        mesh.gather(e, mesh.position, nodes);
        const double dt = explicitLimit(nodes, material);
        if (!(dt > 0.0)) continue;
        smallest = any ? std::min(smallest, dt) : dt;
        any = true;
    }
    return any ? smallest : 0.0;
}

// --- The boundary ------------------------------------------------------------------

namespace {

// The 2x2 rule over one face: shape functions, the area vector, and the area.
struct FaceQuadrature {
    double shape[4][4];   // [gauss point][corner]
    double weight[4];     // |x_s x x_t| at that point, the 2x2 weight being one
    Vec3   normal[4];     // unit, in the face table's own winding
};

bool faceQuadrature(const solidshell::HexMesh& mesh, const std::uint32_t node[4],
                    FaceQuadrature& out) {
    const double q = 1.0 / std::sqrt(3.0);
    double total = 0.0;
    for (int g = 0; g < 4; ++g) {
        const double s = (g & 1) ? q : -q;
        const double t = (g & 2) ? q : -q;
        double ds[4], dt[4];
        for (int c = 0; c < 4; ++c) {
            out.shape[g][c] = 0.25 * (1.0 + s * kCorner[c][0]) * (1.0 + t * kCorner[c][1]);
            ds[c] = 0.25 * kCorner[c][0] * (1.0 + t * kCorner[c][1]);
            dt[c] = 0.25 * kCorner[c][1] * (1.0 + s * kCorner[c][0]);
        }
        Vec3 xs{}, xt{};
        for (int c = 0; c < 4; ++c) {
            const std::size_t n = node[c] * 3;
            const Vec3 p{mesh.position[n], mesh.position[n + 1], mesh.position[n + 2]};
            xs += p * ds[c];
            xt += p * dt[c];
        }
        const Vec3 area = cross(xs, xt);
        out.weight[g] = length(area);
        out.normal[g] = normalize(area);
        total += out.weight[g];
    }
    return total > 0.0;
}

}  // namespace

std::vector<BoundaryFace> boundaryFaces(const solidshell::HexMesh& mesh) {
    std::vector<BoundaryFace> out;
    // A face carried by two elements is interior. Keyed on the sorted node set, so
    // the two elements' opposite windings still match -- the same key
    // `solidshell::uniformPressureLoad` uses, for the same reason.
    std::map<std::array<std::uint32_t, 4>, int> shared;
    for (std::size_t e = 0; e < mesh.elementCount(); ++e)
        for (const auto& face : kFaces) {
            std::array<std::uint32_t, 4> key{};
            for (int i = 0; i < 4; ++i)
                key[static_cast<std::size_t>(i)] =
                    mesh.index[e * kNodes + static_cast<std::size_t>(face[i])];
            std::sort(key.begin(), key.end());
            ++shared[key];
        }

    for (std::size_t e = 0; e < mesh.elementCount(); ++e) {
        Vec3 elementCentre{};
        for (int a = 0; a < kNodes; ++a) {
            const std::size_t n = mesh.index[e * kNodes + static_cast<std::size_t>(a)] * 3;
            elementCentre += Vec3{mesh.position[n], mesh.position[n + 1], mesh.position[n + 2]};
        }
        elementCentre *= 1.0 / kNodes;

        for (int f = 0; f < 6; ++f) {
            BoundaryFace bf;
            std::array<std::uint32_t, 4> key{};
            for (int i = 0; i < 4; ++i) {
                bf.node[i] = mesh.index[e * kNodes + static_cast<std::size_t>(kFaces[f][i])];
                key[static_cast<std::size_t>(i)] = bf.node[i];
            }
            std::sort(key.begin(), key.end());
            if (shared[key] != 1) continue;

            FaceQuadrature fq;
            if (!faceQuadrature(mesh, bf.node, fq)) continue;  // collapsed to a line

            bf.element = static_cast<std::uint32_t>(e);
            bf.face = static_cast<std::uint8_t>(f);
            for (int g = 0; g < 4; ++g) {
                bf.area += fq.weight[g];
                bf.normal += fq.normal[g] * fq.weight[g];
            }
            if (!(bf.area > 0.0)) continue;
            for (int i = 0; i < 4; ++i) {
                const std::size_t n = bf.node[i] * 3;
                bf.centroid += Vec3{mesh.position[n], mesh.position[n + 1], mesh.position[n + 2]};
            }
            bf.centroid *= 0.25;
            bf.normal = normalize(bf.normal);
            // The table's winding is not consistently outward; the element's own
            // centroid says which way is out and cannot be transcribed wrongly.
            if (dot(bf.normal, bf.centroid - elementCentre) < 0.0) bf.normal = -bf.normal;
            out.push_back(bf);
        }
    }
    return out;
}

// --- The solver ---------------------------------------------------------------------

namespace {
// The datum every enthalpy here is measured above: 20 C, which is where
// EN 1993-1-2's tables start.
constexpr double kDatum = kCelsius + 20.0;
}  // namespace

double Solver::specificEnthalpy(double kelvin) const {
    if (problem_.temperatureDependent) return carbonSteelEnthalpy(kelvin);
    return problem_.material.specificHeat * (kelvin - kDatum);
}

bool Solver::prepare(const Problem& problem, double uniformTemperature, std::string* why) {
    const std::size_t n = problem.mesh != nullptr ? problem.mesh->nodeCount() : 0;
    return prepare(problem, std::vector<double>(n, uniformTemperature), why);
}

bool Solver::prepare(const Problem& setup, const std::vector<double>& initialTemperature,
                     std::string* why) {
    ready_ = false;
    factored_ = false;
    propertiesFresh_ = false;
    factorisations_ = 0;
    iterations_ = 0;
    time_ = 0;
    account_ = Account{};
    problem_ = setup;

    if (problem_.mesh == nullptr || problem_.mesh->nodeCount() == 0) {
        if (why) *why = "no mesh";
        return false;
    }
    const solidshell::HexMesh& mesh = *problem_.mesh;
    nodes_ = mesh.nodeCount();
    elements_ = mesh.elementCount();
    if (elements_ == 0) {
        if (why) *why = "mesh has no elements";
        return false;
    }
    if (initialTemperature.size() != nodes_) {
        if (why)
            *why = "initial temperature has " + std::to_string(initialTemperature.size()) +
                   " values for " + std::to_string(nodes_) + " nodes";
        return false;
    }

    // A zero conductivity or a zero heat capacity is a solve in which nothing
    // conducts and nothing stores, and every test of it passes trivially. Refuse
    // it rather than return a field of zeros that reads as a converged answer.
    if (!(problem_.material.density > 0.0)) {
        if (why) *why = "density is not positive";
        return false;
    }
    if (!problem_.temperatureDependent &&
        (!(problem_.material.conductivity > 0.0) || !(problem_.material.specificHeat > 0.0))) {
        if (why) *why = "conductivity and specific heat must both be positive";
        return false;
    }
    if (!problem_.prescribed.empty() &&
        (problem_.prescribed.size() != nodes_ || problem_.prescribedValue.size() != nodes_)) {
        if (why) *why = "prescribed temperatures are not one per node";
        return false;
    }
    if (!problem_.volumetricSource.empty() && problem_.volumetricSource.size() != elements_) {
        if (why) *why = "volumetric source is not one per element";
        return false;
    }

    // A node index past the end of the mesh reads out of bounds in every loop
    // below, so it is caught here rather than by a sanitizer on someone else's
    // machine. `HexMesh` carries no invariant that says otherwise.
    for (std::size_t i = 0; i < mesh.index.size(); ++i)
        if (mesh.index[i] >= nodes_) {
            if (why)
                *why = "element " + std::to_string(i / kNodes) + " names node " +
                       std::to_string(mesh.index[i]) + " of " + std::to_string(nodes_);
            return false;
        }

    forms_.assign(elements_, Forms{});
    for (std::size_t e = 0; e < elements_; ++e) {
        double nodes[kDof];
        mesh.gather(e, mesh.position, nodes);
        if (!computeForms(nodes, forms_[e])) {
            if (why) *why = "element " + std::to_string(e) + " is inverted or degenerate";
            return false;
        }
    }
    conductance_.assign(elements_ * kNodes * kNodes, 0.0);
    capacity_.assign(elements_ * kNodes * kNodes, 0.0);

    // The film surfaces, integrated once.
    filmFace_.clear();
    for (std::size_t f = 0; f < problem_.film.size(); ++f) {
        for (const BoundaryFace& bf : problem_.film[f].face) {
            if (bf.element >= elements_) {
                if (why) *why = "film face names element " + std::to_string(bf.element);
                return false;
            }
            if (bf.face >= 6) {
                if (why)
                    *why = "film face names face " + std::to_string(bf.face) + " of a hexahedron";
                return false;
            }
            FilmFace ff;
            ff.film = static_cast<std::uint32_t>(f);
            for (int i = 0; i < 4; ++i) {
                ff.node[i] = mesh.index[bf.element * kNodes +
                                        static_cast<std::size_t>(kFaces[bf.face][i])];
                if (ff.node[i] != bf.node[i]) {
                    if (why)
                        *why = "film face " + std::to_string(bf.element) + "/" +
                               std::to_string(bf.face) + " does not match the mesh";
                    return false;
                }
            }
            FaceQuadrature fq;
            if (!faceQuadrature(mesh, ff.node, fq)) continue;  // zero area carries no heat
            for (int g = 0; g < 4; ++g)
                for (int a = 0; a < 4; ++a) {
                    ff.load[a] += fq.weight[g] * fq.shape[g][a];
                    for (int b = 0; b < 4; ++b)
                        ff.mass[a * 4 + b] += fq.weight[g] * fq.shape[g][a] * fq.shape[g][b];
                }
            filmFace_.push_back(ff);
        }
    }

    temperature_ = initialTemperature;
    for (std::size_t n = 0; n < nodes_; ++n)
        if (!problem_.prescribed.empty() && problem_.prescribed[n])
            temperature_[n] = problem_.prescribedValue[n];
    previous_.assign(nodes_, 0.0);
    load_.assign(nodes_, 0.0);
    work_.assign(nodes_, 0.0);

    number();

    enthalpyStart_ = enthalpyOf(temperature_);
    account_.enthalpy = enthalpyStart_;
    ready_ = true;
    return true;
}

// --- Numbering ---------------------------------------------------------------------
//
// One scalar unknown per free node. Two candidate orderings are built and the
// narrower kept, which is exactly what `section.cpp` and `reduction.cpp` both do
// with this function and for the reason `reduction.hpp` records: Cuthill-McKee is
// **not** unconditionally better than a mesh's own numbering, and measuring costs
// one pass over the adjacency.
//
// The bandwidth is computed on the *free* numbering, not the node numbering. That
// distinction is the whole of the defect `CLAUDE.md` records against the section
// mesher -- an order chosen on the element graph while the solver assembled a
// constrained one -- and a bandwidth that is too small here does not fail loudly,
// it silently drops terms in `BandedSpd::add`.
void Solver::number() {
    const solidshell::HexMesh& mesh = *problem_.mesh;

    // Every node of an element reaches every other. Built by pushing and then
    // uniquing rather than by a visited mark, because a mark has to be cleared
    // against what the *previous* element already added and getting that wrong
    // duplicates edges silently -- which Cuthill-McKee then weights by.
    std::vector<std::vector<std::uint32_t>> adjacency(nodes_);
    for (std::size_t e = 0; e < elements_; ++e)
        for (int a = 0; a < kNodes; ++a) {
            const std::uint32_t na = mesh.index[e * kNodes + static_cast<std::size_t>(a)];
            for (int b = 0; b < kNodes; ++b) {
                const std::uint32_t nb = mesh.index[e * kNodes + static_cast<std::size_t>(b)];
                if (na != nb) adjacency[na].push_back(nb);
            }
        }
    for (std::vector<std::uint32_t>& row : adjacency) {
        std::sort(row.begin(), row.end());
        row.erase(std::unique(row.begin(), row.end()), row.end());
    }

    // rank[node] -> position. `bandwidthReducingOrder` returns the inverse.
    std::vector<std::uint32_t> identity(nodes_), cuthill(nodes_);
    for (std::uint32_t n = 0; n < nodes_; ++n) identity[n] = n;
    {
        const std::vector<std::uint32_t> order = reduction::bandwidthReducingOrder(adjacency);
        if (order.size() == nodes_)
            for (std::uint32_t r = 0; r < order.size(); ++r) cuthill[order[r]] = r;
        else
            cuthill = identity;
    }

    const auto isFree = [&](std::size_t n) {
        return problem_.prescribed.empty() || !problem_.prescribed[n];
    };
    // Free slots in ascending rank, then the band that numbering delivers.
    const auto measure = [&](const std::vector<std::uint32_t>& rank,
                             std::vector<std::ptrdiff_t>& map) {
        std::vector<std::uint32_t> byRank(nodes_);
        for (std::uint32_t n = 0; n < nodes_; ++n) byRank[rank[n]] = n;
        map.assign(nodes_, -1);
        std::size_t slot = 0;
        for (std::uint32_t r = 0; r < nodes_; ++r)
            if (isFree(byRank[r])) map[byRank[r]] = static_cast<std::ptrdiff_t>(slot++);
        std::size_t width = 0;
        for (std::size_t e = 0; e < elements_; ++e) {
            std::ptrdiff_t lo = -1, hi = -1;
            for (int a = 0; a < kNodes; ++a) {
                const std::ptrdiff_t d = map[mesh.index[e * kNodes + static_cast<std::size_t>(a)]];
                if (d < 0) continue;
                if (lo < 0 || d < lo) lo = d;
                if (d > hi) hi = d;
            }
            if (lo >= 0) width = std::max(width, static_cast<std::size_t>(hi - lo));
        }
        // A film couples only the four nodes of one face, which are four nodes of
        // one element, so the element sweep already covers it. Asserting that
        // rather than assuming it costs nothing and would catch a film face wired
        // to the wrong element.
        for (const FilmFace& ff : filmFace_) {
            std::ptrdiff_t lo = -1, hi = -1;
            for (int a = 0; a < 4; ++a) {
                const std::ptrdiff_t d = map[ff.node[a]];
                if (d < 0) continue;
                if (lo < 0 || d < lo) lo = d;
                if (d > hi) hi = d;
            }
            if (lo >= 0) width = std::max(width, static_cast<std::size_t>(hi - lo));
        }
        return slot == 0 ? std::size_t{0} : width;
    };

    std::vector<std::ptrdiff_t> mapIdentity, mapCuthill;
    const std::size_t bandIdentity = measure(identity, mapIdentity);
    const std::size_t bandCuthill = measure(cuthill, mapCuthill);
    if (bandCuthill < bandIdentity) {
        band_ = bandCuthill;
        map_ = mapCuthill;
    } else {
        band_ = bandIdentity;
        map_ = mapIdentity;
    }
    free_ = 0;
    for (std::size_t n = 0; n < nodes_; ++n)
        if (map_[n] >= 0) ++free_;
    rhs_.assign(free_, 0.0);
}

// --- Properties and element matrices -------------------------------------------------

void Solver::refreshProperties(const std::vector<double>& evaluateAt,
                               const std::vector<double>& previous, bool secant) {
    const solidshell::HexMesh& mesh = *problem_.mesh;
    const double rho = problem_.material.density;
    for (std::size_t e = 0; e < elements_; ++e) {
        double k[kGauss], rc[kGauss];
        if (!problem_.temperatureDependent) {
            for (int gp = 0; gp < kGauss; ++gp) {
                k[gp] = problem_.material.conductivity;
                rc[gp] = rho * problem_.material.specificHeat;
            }
        } else {
            double now[kNodes], was[kNodes];
            for (int a = 0; a < kNodes; ++a) {
                const std::uint32_t n = mesh.index[e * kNodes + static_cast<std::size_t>(a)];
                now[a] = evaluateAt[n];
                was[a] = previous[n];
            }
            for (int gp = 0; gp < kGauss; ++gp) {
                double tNow = 0.0, tWas = 0.0;
                for (int a = 0; a < kNodes; ++a) {
                    tNow += forms_[e].shape[gp][a] * now[a];
                    tWas += forms_[e].shape[gp][a] * was[a];
                }
                k[gp] = carbonSteelConductivity(tNow);
                double c;
                if (!secant) {
                    c = carbonSteelSpecificHeat(tNow);
                } else {
                    const double d = tNow - tWas;
                    // The secant, and its limit where the step has not moved. The
                    // cancellation in the numerator is not a hazard here -- the
                    // energy account differences the same two enthalpies, so
                    // whatever `h` rounds to, the account cancels it identically.
                    c = std::abs(d) > 1e-9
                            ? (carbonSteelEnthalpy(tNow) - carbonSteelEnthalpy(tWas)) / d
                            : carbonSteelSpecificHeat(tNow);
                }
                rc[gp] = rho * c;
            }
        }
        conductance(forms_[e], k, &conductance_[e * kNodes * kNodes]);
        capacity(forms_[e], rc, problem_.lumpedCapacity, &capacity_[e * kNodes * kNodes]);
    }

    // The applied load: volumetric sources and whatever the films supply that does
    // not depend on the temperature. The convective `-h T` part is a matrix and
    // lives in the assembly.
    std::fill(load_.begin(), load_.end(), 0.0);
    if (!problem_.volumetricSource.empty())
        for (std::size_t e = 0; e < elements_; ++e) {
            const double qv = problem_.volumetricSource[e];
            if (qv == 0.0) continue;
            for (int gp = 0; gp < kGauss; ++gp)
                for (int a = 0; a < kNodes; ++a)
                    load_[mesh.index[e * kNodes + static_cast<std::size_t>(a)]] +=
                        qv * forms_[e].weight[gp] * forms_[e].shape[gp][a];
        }
    for (const FilmFace& ff : filmFace_) {
        const Film& film = problem_.film[ff.film];
        const double supply = film.flux + film.coefficient * film.ambient;
        if (supply == 0.0) continue;
        for (int a = 0; a < 4; ++a) load_[ff.node[a]] += supply * ff.load[a];
    }
}

// --- Assembly and solve ----------------------------------------------------------------

bool Solver::buildAndSolve(double inverseStep, bool useCapacity, bool refactor, std::string* why) {
    const solidshell::HexMesh& mesh = *problem_.mesh;
    if (free_ == 0) return true;  // everything prescribed: nothing to solve

    if (refactor) system_ = solidshell::BandedSpd(free_, band_);
    std::fill(rhs_.begin(), rhs_.end(), 0.0);

    for (std::size_t n = 0; n < nodes_; ++n)
        if (map_[n] >= 0) rhs_[static_cast<std::size_t>(map_[n])] = load_[n];

    const auto scatter = [&](std::size_t row, std::size_t column, double value) {
        const std::ptrdiff_t r = map_[row];
        if (r < 0) return;
        const std::ptrdiff_t c = map_[column];
        if (c >= 0) {
            if (refactor)
                system_.add(static_cast<std::size_t>(r), static_cast<std::size_t>(c), value);
        } else {
            // A prescribed node moves to the right-hand side exactly, the same way
            // `solidshell::solveStatic` moves a prescribed displacement -- not by a
            // penalty, whose size would then be a tolerance the steady patch test
            // had to live with.
            rhs_[static_cast<std::size_t>(r)] -= value * temperature_[column];
        }
    };

    for (std::size_t e = 0; e < elements_; ++e) {
        const double* ke = &conductance_[e * kNodes * kNodes];
        const double* ce = &capacity_[e * kNodes * kNodes];
        std::uint32_t node[kNodes];
        for (int a = 0; a < kNodes; ++a)
            node[a] = mesh.index[e * kNodes + static_cast<std::size_t>(a)];
        for (int a = 0; a < kNodes; ++a) {
            if (map_[node[a]] < 0) continue;
            const std::size_t r = static_cast<std::size_t>(map_[node[a]]);
            for (int b = 0; b < kNodes; ++b) {
                // C/dt T0 is a load over *every* column, prescribed or not, which
                // is why it is added here and not folded into the scatter.
                if (useCapacity)
                    rhs_[r] += ce[a * kNodes + b] * inverseStep * previous_[node[b]];
                scatter(node[a], node[b],
                        ke[a * kNodes + b] + (useCapacity ? ce[a * kNodes + b] * inverseStep : 0.0));
            }
        }
    }
    for (const FilmFace& ff : filmFace_) {
        const double h = problem_.film[ff.film].coefficient;
        if (h == 0.0) continue;
        for (int a = 0; a < 4; ++a)
            for (int b = 0; b < 4; ++b) scatter(ff.node[a], ff.node[b], h * ff.mass[a * 4 + b]);
    }

    if (refactor) {
        ++factorisations_;
        if (!system_.factor()) {
            if (why)
                *why = useCapacity
                           ? "the system is not positive definite"
                           : "steady conduction is singular: nothing holds the temperature";
            factored_ = false;
            return false;
        }
        factored_ = true;
    }
    system_.solve(rhs_);
    for (std::size_t n = 0; n < nodes_; ++n)
        if (map_[n] >= 0) temperature_[n] = rhs_[static_cast<std::size_t>(map_[n])];
    return true;
}

// `A T1 - b` over every node: zero at a free one, and the heat that had to be
// supplied at a prescribed one.
void Solver::residual(double inverseStep, bool useCapacity, std::vector<double>& out) const {
    const solidshell::HexMesh& mesh = *problem_.mesh;
    out.assign(nodes_, 0.0);
    for (std::size_t e = 0; e < elements_; ++e) {
        const double* ke = &conductance_[e * kNodes * kNodes];
        const double* ce = &capacity_[e * kNodes * kNodes];
        std::uint32_t node[kNodes];
        for (int a = 0; a < kNodes; ++a)
            node[a] = mesh.index[e * kNodes + static_cast<std::size_t>(a)];
        for (int a = 0; a < kNodes; ++a)
            for (int b = 0; b < kNodes; ++b) {
                out[node[a]] += ke[a * kNodes + b] * temperature_[node[b]];
                if (useCapacity)
                    out[node[a]] += ce[a * kNodes + b] * inverseStep *
                                    (temperature_[node[b]] - previous_[node[b]]);
            }
    }
    for (const FilmFace& ff : filmFace_) {
        const double h = problem_.film[ff.film].coefficient;
        if (h == 0.0) continue;
        for (int a = 0; a < 4; ++a)
            for (int b = 0; b < 4; ++b)
                out[ff.node[a]] += h * ff.mass[a * 4 + b] * temperature_[ff.node[b]];
    }
    for (std::size_t n = 0; n < nodes_; ++n) out[n] -= load_[n];
}

double Solver::enthalpyOf(const std::vector<double>& field) const {
    const solidshell::HexMesh& mesh = *problem_.mesh;
    const double rho = problem_.material.density;
    double total = 0.0;
    for (std::size_t e = 0; e < elements_; ++e) {
        double t[kNodes];
        for (int a = 0; a < kNodes; ++a)
            t[a] = field[mesh.index[e * kNodes + static_cast<std::size_t>(a)]];
        for (int gp = 0; gp < kGauss; ++gp) {
            double tgp = 0.0;
            for (int a = 0; a < kNodes; ++a) tgp += forms_[e].shape[gp][a] * t[a];
            total += forms_[e].weight[gp] * rho * specificEnthalpy(tgp);
        }
    }
    return total;
}

bool Solver::step(double timestep, std::string* why) {
    if (!ready_) {
        if (why) *why = "solver was not prepared";
        return false;
    }
    if (!(timestep > 0.0)) {
        if (why) *why = "timestep is not positive";
        return false;
    }
    previous_ = temperature_;
    const double invDt = 1.0 / timestep;

    bool converged = true;
    iterations_ = 0;
    if (!problem_.temperatureDependent) {
        // Constant properties: the element matrices and the applied load are the
        // same every step, and so is `C/dt + K` at a fixed step. Rebuilding either
        // would be pure waste -- the same argument `solidshell::RestForms` makes,
        // and here it removes the *factorisation* and not merely the forms.
        if (!propertiesFresh_) {
            refreshProperties(temperature_, previous_, /*secant=*/true);
            propertiesFresh_ = true;
        }
        const bool refactor = !factored_ || factoredStep_ != timestep;
        if (!buildAndSolve(invDt, true, refactor, why)) return false;
        factoredStep_ = timestep;
        iterations_ = 1;
    } else {
        // Picard. Properties are re-evaluated at the current iterate and the system
        // is re-formed and re-factored, which is why `factorisations()` is worth
        // reporting: a nonlinear step costs one factorisation per iteration where a
        // linear one costs none at all after the first.
        std::vector<double> before;
        converged = false;
        for (int it = 0; it < picardLimit_; ++it) {
            before = temperature_;
            refreshProperties(temperature_, previous_, /*secant=*/true);
            if (!buildAndSolve(invDt, true, /*refactor=*/true, why)) return false;
            ++iterations_;
            double move = 0.0;
            for (std::size_t n = 0; n < nodes_; ++n)
                move = std::max(move, std::abs(temperature_[n] - before[n]));
            if (move <= picardTolerance_) {
                converged = true;
                break;
            }
        }
        // The properties, and therefore `C`, must match the temperature that was
        // finally accepted, or `1^T C dT` is the enthalpy change of a state that
        // was never reached and the account carries a residual that looks like a
        // quadrature error.
        refreshProperties(temperature_, previous_, /*secant=*/true);
        factored_ = false;
        propertiesFresh_ = false;
    }

    residual(invDt, true, work_);
    double worstFree = 0.0, reaction = 0.0;
    for (std::size_t n = 0; n < nodes_; ++n) {
        if (map_[n] >= 0)
            worstFree = std::max(worstFree, std::abs(work_[n]));
        else
            reaction += work_[n];
    }
    account_.equilibriumResidual = worstFree;
    account_.prescribedHeat += reaction * timestep;

    double filmPower = 0.0, sourcePower = 0.0;
    for (const FilmFace& ff : filmFace_) {
        const Film& film = problem_.film[ff.film];
        for (int a = 0; a < 4; ++a) {
            filmPower += (film.flux + film.coefficient * film.ambient) * ff.load[a];
            if (film.coefficient == 0.0) continue;
            for (int b = 0; b < 4; ++b)
                filmPower -= film.coefficient * ff.mass[a * 4 + b] * temperature_[ff.node[b]];
        }
    }
    if (!problem_.volumetricSource.empty())
        for (std::size_t e = 0; e < elements_; ++e)
            for (int gp = 0; gp < kGauss; ++gp)
                sourcePower += problem_.volumetricSource[e] * forms_[e].weight[gp];
    account_.filmHeat += filmPower * timestep;
    account_.sourceHeat += sourcePower * timestep;

    account_.enthalpy = enthalpyOf(temperature_);
    account_.enthalpyChange = account_.enthalpy - enthalpyStart_;
    time_ += timestep;

    if (!converged) {
        if (why)
            *why = "Picard did not converge in " + std::to_string(picardLimit_) + " iterations";
        return false;
    }
    return true;
}

bool Solver::setFilm(std::size_t index, double flux, double coefficient, double ambient) {
    if (index >= problem_.film.size()) return false;
    Film& film = problem_.film[index];
    film.flux = flux;
    film.coefficient = coefficient;
    film.ambient = ambient;
    // `coefficient` is in the matrix and `coefficient * ambient` is in the load, so
    // both have to be rebuilt. Neither the forms nor the numbering depend on either.
    factored_ = false;
    propertiesFresh_ = false;
    return true;
}

bool Solver::solveSteady(std::string* why) {
    if (!ready_) {
        if (why) *why = "solver was not prepared";
        return false;
    }
    previous_ = temperature_;
    if (!problem_.temperatureDependent) {
        refreshProperties(temperature_, previous_, /*secant=*/false);
        propertiesFresh_ = true;
        if (!buildAndSolve(0.0, false, /*refactor=*/true, why)) return false;
        iterations_ = 1;
    } else {
        std::vector<double> before;
        bool converged = false;
        iterations_ = 0;
        for (int it = 0; it < picardLimit_; ++it) {
            before = temperature_;
            refreshProperties(temperature_, previous_, /*secant=*/false);
            if (!buildAndSolve(0.0, false, /*refactor=*/true, why)) return false;
            ++iterations_;
            double move = 0.0;
            for (std::size_t n = 0; n < nodes_; ++n)
                move = std::max(move, std::abs(temperature_[n] - before[n]));
            if (move <= picardTolerance_) {
                converged = true;
                break;
            }
        }
        refreshProperties(temperature_, previous_, /*secant=*/false);
        if (!converged) {
            if (why)
                *why = "Picard did not converge in " + std::to_string(picardLimit_) + " iterations";
            return false;
        }
    }
    // The steady system is `K` alone, so a transient step must build its own.
    factored_ = false;
    propertiesFresh_ = false;

    residual(0.0, false, work_);
    double worstFree = 0.0;
    for (std::size_t n = 0; n < nodes_; ++n)
        if (map_[n] >= 0) worstFree = std::max(worstFree, std::abs(work_[n]));
    account_.equilibriumResidual = worstFree;
    account_.enthalpy = enthalpyOf(temperature_);
    account_.enthalpyChange = account_.enthalpy - enthalpyStart_;
    return true;
}

}  // namespace sim::thermal
