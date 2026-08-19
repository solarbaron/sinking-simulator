// SPDX-License-Identifier: MIT
#include "roll_damping.hpp"

#include <algorithm>
#include <cmath>
#include <string>

namespace sim {
namespace {

double sqr(double x)  { return x * x; }
double cube(double x) { return x * x * x; }

// Everything the component formulae share, worked out once.
struct Geometry {
    double L = 0, B = 0, d = 0, cb = 0, cm = 0;
    double og = 0;      // OG: still-water level down to the roll axis, + downward
    double ogOverD = 0; // OG/d
    double h0 = 0;      // half-beam/draft ratio, B/(2d)
    double bilgeR = 0;  // bilge radius, m
    double volume = 0;  // displacement volume, m^3
    bool   valid = false;
};

Geometry derive(const RollDampingHull& h) {
    Geometry g;
    if (h.lengthPp <= 0 || h.beam <= 0 || h.draft <= 0 || h.blockCoeff <= 0 ||
        h.midshipCoeff <= 0 || h.seaDensity <= 0)
        return g;
    g.L = h.lengthPp;
    g.B = h.beam;
    g.d = h.draft;
    g.cb = h.blockCoeff;
    g.cm = h.midshipCoeff;
    g.og = h.draft - h.rollAxisAboveKeel;
    g.ogOverD = g.og / g.d;
    g.h0 = g.B / (2.0 * g.d);
    g.bilgeR = h.bilgeRadiusOrDefault();
    g.volume = h.displacementVolume();
    g.valid = true;
    return g;
}

// --- Friction: Kato (1958), ITTC (2.13)-(2.17) ------------------------------
//
// Kato's frictional coefficient Cf = 1.328 (3.22 rf^2 phi_a^2 / (T nu))^(-1/2)
// is inversely proportional to the roll amplitude, and B44F0 carries a factor
// phi_a, so the amplitude cancels *exactly*:
//
//   B44F0 = 4/(3 pi) rho Sf rf^3 phi_a omega Cf
//         = 4/(3 pi) rho Sf rf^2 omega (1.328 / sqrt(3.22)) sqrt(T nu)
//
// and with T = 2 pi / omega the frequency dependence collapses to sqrt(omega).
// Written in the cancelled form the result is finite and correct at phi_a = 0,
// which is where the naive form would divide by zero. Friction really is
// amplitude-independent; that is a property of the model, not an approximation.
double frictionDamping(const RollDampingHull& h, const Geometry& g, double omega, double speed) {
    if (omega <= 0) return 0.0;
    // Wetted girth. ITTC (2.16) prints 1.7 d and Kawahara's eq. (23) prints
    // 1.75 d for the same quantity; the sources genuinely disagree. It is a 3%
    // difference in a term that is 1-3% of a full-scale ship's roll damping, so
    // the ITTC value is used and the discrepancy noted rather than tested for.
    const double girth = 1.7 * g.d + g.cb * g.B;
    const double sf = g.L * girth;                                   // ITTC (2.16)
    const double rf = ((0.887 + 0.145 * g.cb) * girth - 2.0 * g.og)  // ITTC (2.15)
                      / kPi;
    if (rf <= 0) return 0.0;

    const double kKato = 1.328 / std::sqrt(3.22);  // 0.74005...
    const double b0 = 4.0 / (3.0 * kPi) * h.seaDensity * sf * sqr(rf) * kKato *
                      std::sqrt(2.0 * kPi * h.kinematicViscosity * omega);

    // Tamiya's forward-speed correction, ITTC (2.17).
    return b0 * (1.0 + 4.1 * std::abs(speed) / (omega * g.L));
}

// --- Eddy: simplified Ikeda, Kawahara et al. (2009) eq. (31) ----------------
//
// The regression returns C_R for the whole ship; the nondimensional eddy damping
// is then
//   B44E0hat = 4 Lpp d^4 omegahat phi_a / (3 pi volume B^2) * C_R
//            = 4 omegahat phi_a / (3 pi Cb (B/d)^3) * C_R
// (the two forms are identical once volume = Cb Lpp B d is substituted).
double eddyCoefficientCR(double x1, double x2, double x3, double x4) {
    const double x2_2 = sqr(x2), x2_3 = x2_2 * x2, x2_4 = x2_3 * x2;

    const double aE = (-0.0182 * x2 + 0.0155) * cube(x1 - 1.8) - 79.414 * x2_4 +
                      215.695 * x2_3 - 215.883 * x2_2 + 93.894 * x2 - 14.848;
    const double bE1 = (-0.2 * x1 + 1.6) * (3.98 * x2 - 5.1525) * x4 *
                       ((0.9717 * x2_2 - 1.55 * x2 + 0.723) * x4 + (0.04567 * x2 + 0.9408));
    const double bE2 = (0.25 * x4 + 0.95) * x4 - 219.2 * x2_3 + 443.7 * x2_2 - 283.3 * x2 + 59.6;
    const double bE3 = (46.5 - 15.0 * x1) * x2 + 11.2 * x1 - 28.6;

    return aE * std::exp(bE1 + bE2 * std::pow(x3, bE3));
}

double eddyDamping(const RollDampingHull& h, const Geometry& g, double amplitude, double omega,
                   double speed) {
    if (amplitude <= 0 || omega <= 0) return 0.0;  // strictly proportional to both
    const double x1 = g.B / g.d, x2 = g.cb, x3 = g.cm, x4 = g.ogOverD;
    const double cr = eddyCoefficientCR(x1, x2, x3, x4);
    if (!(cr > 0)) return 0.0;  // the polynomial fit can go negative off its domain

    const double omegaHat = omega * std::sqrt(g.B / (2.0 * kGravity));
    const double hat0 = 4.0 * omegaHat * amplitude / (3.0 * kPi * x2 * cube(x1)) * cr;
    const double b0 = hat0 * h.nondimensionalScale();

    // ITTC (2.21): the eddy component collapses with forward speed as the bilge
    // vortices are swept away. K = omega L / U, factor = (0.04 K)^2/(1+(0.04 K)^2),
    // rearranged so that U = 0 needs no special case.
    const double k = 0.04 * omega * g.L;
    return b0 * sqr(k) / (sqr(k) + sqr(speed));
}

// --- Hull lift: ITTC (2.10)-(2.11) ------------------------------------------
//
// A hull moving ahead while rolling meets the flow at an angle of attack and
// generates lift, which opposes the roll. Strictly zero at zero speed and linear
// in speed. Independent of roll amplitude -- it is a linear mechanism.
double liftDamping(const RollDampingHull& h, const Geometry& g, double speed) {
    const double u = std::abs(speed);
    if (u <= 0) return 0.0;

    double kappa = 0.0;
    if (g.cm > 0.97)      kappa = 0.3;
    else if (g.cm > 0.92) kappa = 0.1;

    const double kN = 2.0 * kPi * g.d / g.L + kappa * (4.1 * g.B / g.L - 0.045);
    const double l0 = 0.3 * g.d;
    const double lR = 0.5 * g.d;
    const double shape = 1.0 + 1.4 * g.og / lR + 0.7 * sqr(g.og) / (l0 * lR);

    return 0.5 * h.seaDensity * u * g.L * g.d * kN * l0 * lR * shape;
}

// --- Bilge keels: ITTC (2.24)-(2.32) ----------------------------------------
//
// Two mechanisms from the same vortices: the normal force on the keel itself,
// and the pressure field the keel induces on the shell plating ahead of and
// behind it. The section is idealised as a vertical side, a horizontal bottom
// and a quarter-circle bilge, with the keel at the middle of the arc and normal
// to it -- Ikeda's own assumption, and the reason a chine hull needs a different
// method.
//
// The drag coefficient CD = 22.5 b_BK/(pi l f phi_a) + 2.4 diverges as phi_a
// tends to zero, and it multiplies phi_a. The product is finite, so both
// components are evaluated as phi_a * CD rather than as CD, which keeps them
// exact and well-defined at zero amplitude. That limit is not zero: a bilge keel
// still damps an infinitesimal roll, because the Keulegan-Carpenter number falls
// with the amplitude and the drag coefficient rises to meet it. Only the
// *increase* with amplitude vanishes.
struct BilgeKeel {
    double normal = 0;
    double hull = 0;
};

BilgeKeel bilgeKeelDamping(const RollDampingHull& h, const Geometry& g, double amplitude,
                           double omega) {
    BilgeKeel out;
    // No keel, no vortex, no damping. The formulae below are empirical fits for
    // a real keel and do not tend to zero as its breadth does (Cp+ is a constant
    // 1.2), so the absence of a keel is handled here rather than extrapolated.
    if (h.bilgeKeelBreadth <= 0 || h.bilgeKeelLength <= 0 || omega <= 0 || amplitude < 0)
        return out;

    const double bBK = h.bilgeKeelBreadth;
    const double lBK = h.bilgeKeelLength;
    const double m1 = g.bilgeR / g.d;
    const double m2 = g.ogOverD;
    const double m3 = 1.0 - m1 - m2;
    const double m4 = g.h0 - m1;

    // Distance from the roll axis to the tip of the keel, ITTC (2.31). The keel
    // sits at the 45-degree point of the bilge arc, which is (1 - sqrt2/2) R in
    // from both the side and the bottom.
    const double inset = (1.0 - std::sqrt(2.0) / 2.0) * m1;
    const double l = g.d * std::sqrt(sqr(g.h0 - inset) + sqr(1.0 - m2 - inset));
    if (l <= 0) return out;

    // Local flow speed-up around the turn of the bilge, ITTC (2.26).
    const double f = 1.0 + 0.3 * std::exp(-160.0 * (1.0 - g.cm));

    // phi_a * CD and phi_a * Cp-, both finite at phi_a = 0. ITTC (2.24), (2.29).
    const double kcTerm = 22.5 * bBK / (kPi * l * f);
    const double phiCd = kcTerm + 2.4 * amplitude;
    const double phiCpMinus = -(kcTerm + 1.2 * amplitude);
    const double phiCpPlus = 1.2 * amplitude;

    out.normal = 8.0 / (3.0 * kPi) * h.seaDensity * cube(l) * omega * bBK * sqr(f) * phiCd * lBK;

    // Length of the negative-pressure region behind the keel, ITTC (2.27).
    const double s0 = 0.3 * kPi * l * f * amplitude + 1.95 * bBK;
    const double arc = 0.25 * kPi * g.bilgeR;

    double m7 = 0.0, m8 = 0.0;
    if (s0 > arc) {
        m7 = s0 / g.d - 0.25 * kPi * m1;
        m8 = m7 + 0.414 * m1;
    } else {
        m7 = 0.0;
        m8 = g.bilgeR > 0 ? 1.414 * m1 * (1.0 - std::cos(s0 / g.bilgeR)) : 0.0;
    }

    const double denom = (g.h0 - 0.215 * m1) * (1.0 - 0.215 * m1);
    if (std::abs(denom) < 1e-12) return out;
    const double m5 = (0.414 * g.h0 + 0.0651 * sqr(m1) - (0.382 * g.h0 + 0.0106) * m1) / denom;
    const double m6 = (0.414 * g.h0 + 0.0651 * sqr(m1) - (0.382 + 0.0106 * g.h0) * m1) / denom;

    // ITTC (2.30): the girthwise integral of the induced pressure times its
    // moment lever, split into the part multiplying Cp- and the part multiplying
    // Cp+. Both A0 and B0 are pure section geometry.
    //
    // The two signs here -- the leading sign of the m1 term in B0, and the minus
    // on A0 * Cp- -- were pinned by agreement with Kawahara's independent
    // regression, which the wrong choice misses by 15% and the right one by 3%.
    // tests/test_roll_damping.cpp keeps that comparison, at a tolerance tight
    // enough to catch a future flip.
    //
    // **`cube(m2)`, not `sqr(m2)`.** Each term of A0 and B0 is a girth times a
    // lever, both nondimensionalised on d, so each must be degree 2 in the m's --
    // and the elided `1`s are `d/d`, degree 1, not degree 0. Term two is degree 3
    // over degree 1; term three is `m1` times a degree-1 lever times a degree-0
    // ratio. Term one was degree 2 over degree 1, the only term in either integral
    // that did not scale like the rest. That argument needs no source document.
    //
    // Which cube it should be does. The regression settles it: swept over Ikeda's
    // own validated `OG/d` range and nothing else, the sectional model misses
    // Kawahara by
    //
    //     sqr(m2)    3.6% at the fixture, rising monotonically to 15.6% at -1.5
    //     cube(m2)   4.9% worst anywhere on the range, at -0.25
    //     cube(m4)   9-15% everywhere, worst 23.7%
    //
    // `cube(m4)` is what a first-principles reading gives if term one is the flat
    // bottom's integral -- pressure rising linearly with girth from the centreline,
    // lever the horizontal offset, which lands on `m4^3/(3(H0 - 0.215 m1))` exactly
    // the way term two's `[(1-m1)^3/3 - m2 (1-m1)^2/2]/(1 - 0.215 m1)` falls out of
    // the side's. It is refuted: 9-15% is not a regression's scatter. So the
    // reading is wrong, term one is not the flat bottom, and this comment does not
    // claim to know what it is.
    //
    // **What the old test could not see.** It swept omega and amplitude, nine
    // points, and held the hull at `OG/d = -0.6923`. Neither omega nor amplitude
    // distinguishes the two readings -- `m2` contains neither -- and -0.69 is
    // within 0.05 of where `sqr` and `cube` cross. A nine-point grid that varies
    // only the axes the defect is constant along is a one-point grid.
    const double a0 = (m3 + m4) * m8 - sqr(m7);
    const double b0 = cube(m2) / (3.0 * (g.h0 - 0.215 * m1)) +
                      sqr(1.0 - m1) * (2.0 * m3 - m2) / (6.0 * (1.0 - 0.215 * m1)) +
                      m1 * (m3 * m5 + m4 * m6);
    const double integral = sqr(g.d) * (-a0 * phiCpMinus + b0 * phiCpPlus);

    out.hull = 4.0 / (3.0 * kPi) * h.seaDensity * sqr(l) * sqr(f) * omega * integral * lBK;

    if (out.normal < 0) out.normal = 0;
    if (out.hull < 0) out.hull = 0;
    return out;
}

}  // namespace

double RollDampingHull::displacementVolume() const {
    return blockCoeff * lengthPp * beam * draft;
}

// ITTC (2.32). A quarter-circle bilge of this radius reproduces the midship
// section area coefficient of a hull with the given B/d, clamped so it can never
// exceed the draft or the half-beam.
double RollDampingHull::bilgeRadiusOrDefault() const {
    if (bilgeRadius >= 0) return bilgeRadius;
    if (draft <= 0 || beam <= 0 || midshipCoeff <= 0 || midshipCoeff >= 1.0) return 0.0;
    const double h0 = beam / (2.0 * draft);
    double r = 2.0 * draft * std::sqrt(h0 * (midshipCoeff - 1.0) / (kPi - 4.0));
    r = std::min(r, h0 >= 1.0 ? draft : 0.5 * beam);
    return std::max(0.0, r);
}

double RollDampingHull::nondimensionalScale() const {
    const double v = displacementVolume();
    if (v <= 0 || beam <= 0) return 0.0;
    return seaDensity * v * sqr(beam) / std::sqrt(beam / (2.0 * kGravity));
}

RollDamping rollDamping(const RollDampingHull& hull, const RollDampingCondition& condition) {
    RollDamping out;
    const Geometry g = derive(hull);
    if (!g.valid) return out;

    const double amplitude = std::max(0.0, condition.rollAmplitude);
    const double omega = condition.rollFrequency;
    const double speed = condition.forwardSpeed;

    out.friction = frictionDamping(hull, g, omega, speed);
    out.eddy = eddyDamping(hull, g, amplitude, omega, speed);
    out.lift = liftDamping(hull, g, speed);
    const BilgeKeel bk = bilgeKeelDamping(hull, g, amplitude, omega);
    out.bilgeKeelNormal = bk.normal;
    out.bilgeKeelHull = bk.hull;
    out.wave = hull.waveDamping;

    out.total = out.friction + out.eddy + out.lift + out.bilgeKeelNormal + out.bilgeKeelHull +
                out.wave;
    const double scale = hull.nondimensionalScale();
    out.totalHat = scale > 0 ? out.total / scale : 0.0;
    return out;
}

RollDampingHull rollDampingHullFromMesh(const TriMesh& hull, double waterlineZ,
                                        double bilgeKeelLength, double bilgeKeelBreadth,
                                        double density) {
    RollDampingHull out;
    out.seaDensity = density;
    out.bilgeKeelLength = bilgeKeelLength;
    out.bilgeKeelBreadth = bilgeKeelBreadth;
    if (hull.verts.empty()) return out;

    Vec3 lo = hull.verts[0], hi = hull.verts[0];
    for (const Vec3& v : hull.verts) {
        lo = {std::min(lo.x, v.x), std::min(lo.y, v.y), std::min(lo.z, v.z)};
        hi = {std::max(hi.x, v.x), std::max(hi.y, v.y), std::max(hi.z, v.z)};
    }
    if (!(waterlineZ > lo.z)) return out;

    // Clip generously wide in the directions that are not being cut, so the box
    // trims the hull only at the waterline and the shell does the rest.
    const double wide = (hi.x - lo.x) + (hi.y - lo.y) + (hi.z - lo.z) + 1.0;
    const Vec3 boxLo{lo.x - wide, lo.y - wide, lo.z - wide};
    const TriMesh wetted = clipToBox(hull, boxLo, {hi.x + wide, hi.y + wide, waterlineZ});
    if (wetted.verts.empty()) return out;

    Vec3 wlo = wetted.verts[0], whi = wetted.verts[0];
    for (const Vec3& v : wetted.verts) {
        wlo = {std::min(wlo.x, v.x), std::min(wlo.y, v.y), std::min(wlo.z, v.z)};
        whi = {std::max(whi.x, v.x), std::max(whi.y, v.y), std::max(whi.z, v.z)};
    }
    // Beam is taken over the whole wetted body rather than at the waterline
    // exactly. For any hull without tumblehome those are the same number, and for
    // one with it the wider figure is what the bilge-keel lever arms are built on.
    out.lengthPp = whi.x - wlo.x;
    out.beam = whi.y - wlo.y;
    out.draft = waterlineZ - wlo.z;
    if (out.lengthPp <= 0 || out.beam <= 0 || out.draft <= 0) return out;

    const double volume = integrate(wetted).volume;
    if (volume <= 0) return out;
    out.blockCoeff = volume / (out.lengthPp * out.beam * out.draft);

    // Cm is the *midship* section coefficient, and midship on a real hull is
    // wherever the section is largest -- not necessarily x = 0, and certainly not
    // the mean of the stations. Scanning for the maximum is what makes this right
    // for a hull with its parallel midbody off centre.
    constexpr int kSlabs = 41;
    const double thickness = out.lengthPp / kSlabs;
    double largestArea = 0;
    for (int i = 0; i < kSlabs; ++i) {
        const double x = wlo.x + out.lengthPp * (i + 0.5) / kSlabs;
        const TriMesh slab = clipToBox(hull, {x - 0.5 * thickness, boxLo.y, boxLo.z},
                                       {x + 0.5 * thickness, hi.y + wide, waterlineZ});
        if (slab.verts.empty()) continue;
        largestArea = std::max(largestArea, integrate(slab).volume / thickness);
    }
    out.midshipCoeff = largestArea / (out.beam * out.draft);
    return out;
}

std::vector<std::string> validateRollDamping(const RollDampingHull& hull,
                                             const RollDampingCondition& condition) {
    std::vector<std::string> problems;
    auto note = [&](const std::string& s) { problems.push_back(s); };
    auto range = [&](const char* what, double value, double lo, double hi) {
        if (value < lo || value > hi)
            note(std::string(what) + " is " + std::to_string(value) + ", outside Ikeda's range [" +
                 std::to_string(lo) + ", " + std::to_string(hi) + "]");
    };

    if (hull.lengthPp <= 0) note("length between perpendiculars is not positive");
    if (hull.beam <= 0) note("beam is not positive");
    if (hull.draft <= 0) note("draft is not positive");
    if (hull.seaDensity <= 0) note("sea density is not positive");
    if (hull.kinematicViscosity <= 0) note("kinematic viscosity is not positive");
    if (condition.rollFrequency <= 0) note("roll frequency is not positive");
    if (condition.rollAmplitude < 0) note("roll amplitude is negative");
    if (!problems.empty()) return problems;

    // The regression domain of the simplified method (Kawahara et al. 2009).
    range("block coefficient", hull.blockCoeff, 0.5, 0.85);
    range("beam/draft ratio", hull.beam / hull.draft, 2.5, 4.5);
    range("midship coefficient", hull.midshipCoeff, 0.9, 0.99);
    range("OG/draft (roll axis below the waterline is positive)",
          (hull.draft - hull.rollAxisAboveKeel) / hull.draft, -1.5, 0.2);

    if (hull.bilgeKeelBreadth > 0 && hull.bilgeKeelLength > 0) {
        range("bilge keel breadth / beam", hull.bilgeKeelBreadth / hull.beam, 0.01, 0.06);
        range("bilge keel length / Lpp", hull.bilgeKeelLength / hull.lengthPp, 0.05, 0.4);
    }
    // **Kawahara's C_R collapses to zero at the top of its own Cb range.** The
    // quartic in Cb inside `aE` has a root at about 0.844, so a hull at Cb = 0.85
    // -- which the range check above admits, because 0.85 is the published bound --
    // comes back with `C_R <= 0`, and `eddyDamping` then returns exactly zero
    // through its `if (!(cr > 0))` guard. For a bare hull the eddy term is over 90%
    // of the viscous total, so that is the quietest failure this method has: the
    // dominant component silently absent on a hull the validator just passed.
    //
    // Found by sweeping the four regression inputs across the declared box, which
    // nothing had ever done -- every point at Cb = 0.85 was dead and every point
    // below 0.84 was healthy. Reported rather than clamped: the bound is Kawahara's
    // and moving it would be inventing coverage the fit does not have.
    if (hull.draft > 0 && hull.beam > 0) {
        const double ogOverD = (hull.draft - hull.rollAxisAboveKeel) / hull.draft;
        if (!(eddyCoefficientCR(hull.beam / hull.draft, hull.blockCoeff, hull.midshipCoeff,
                                ogOverD) > 0))
            note("the eddy regression is non-positive at this hull's form, so the eddy "
                 "component -- most of the viscous damping on a bare hull -- comes out "
                 "as exactly zero; the fit's Cb quartic has a root near 0.844");
    }

    if (hull.midshipCoeff > 0.99)
        note("midship coefficient exceeds 0.99, past the last entry of Ikeda's lift-slope table");
    if (hull.waveDamping <= 0)
        note("wave (radiation) damping is not supplied, so the total omits 5-30% of the real "
             "damping -- see docs/02-simulation.md");

    // Ikeda's fits come from conventional cargo-hull sections. Anything far from
    // that -- a hard chine, a very flat barge, a buttock-flow stern -- needs a
    // different method entirely, and there is no parameter that reveals it, so
    // the best that can be done is to say where the numbers came from.
    return problems;
}

}  // namespace sim
