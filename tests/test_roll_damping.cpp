// SPDX-License-Identifier: MIT
//
// Validation of Ikeda's viscous roll damping.
//
// An empirical method cannot be checked against first principles -- that is what
// makes it empirical -- so it is checked against the things that must hold no
// matter what the fitted constants are:
//
//   * exact algebraic properties of each component (the eddy moment is quadratic
//     in roll rate, so its linear equivalent is exactly proportional to omega
//     and to roll amplitude; friction's amplitude cancels; lift is exactly zero
//     at zero speed and exactly linear in speed);
//   * Froude scale invariance of the nondimensional coefficient, which every
//     component except friction must satisfy exactly, and which friction must
//     violate in a specific way (it is the only component with a Reynolds number
//     in it);
//   * agreement between two independently published formulations of the same
//     component -- the bilge-keel damping is computed here from Ikeda's sectional
//     pressure model and compared against Kawahara's regression, which was fitted
//     to it. Two routes to the same number is the strongest check available;
//   * and finally the closed form that says what the coefficient *means*: feed
//     B44 into a linear 1-DOF roll equation, integrate a free decay, and the
//     logarithmic decrement must be 2 pi zeta / sqrt(1 - zeta^2).
#include "engine/sim/hullform.hpp"
#include "engine/sim/roll_damping.hpp"
#include "harness.hpp"

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

using namespace sim;
using testing::expectNear;
using testing::expectTrue;

namespace {

// The test hull: a 170 m ro-pax. Cb 0.55 and Cm 0.98 are conventional for the
// type, KG 11 m puts the roll axis 4.5 m above the waterline (OG/d = -0.69) as a
// high-sided vehicle ferry has, and the bilge keels are 34 m x 0.6 m -- 0.20 Lpp
// by 0.024 B. Every parameter sits inside the published validity range of the
// method, which testValidityIsReported() confirms rather than assumes.
RollDampingHull roPax() {
    RollDampingHull h;
    h.lengthPp = 170.0;
    h.beam = 25.0;
    h.draft = 6.5;
    h.blockCoeff = 0.55;
    h.midshipCoeff = 0.98;
    h.rollAxisAboveKeel = 11.0;
    h.bilgeKeelLength = 34.0;
    h.bilgeKeelBreadth = 0.6;
    return h;
}

RollDampingHull bareRoPax() {
    RollDampingHull h = roPax();
    h.bilgeKeelLength = 0;
    h.bilgeKeelBreadth = 0;
    return h;
}

RollDampingCondition at(double amplitudeDeg, double omega, double speed = 0.0) {
    RollDampingCondition c;
    c.rollAmplitude = amplitudeDeg * kDegToRad;
    c.rollFrequency = omega;
    c.forwardSpeed = speed;
    return c;
}

// The natural roll frequency of the test hull, so that the coefficients are
// quoted where they are actually used. k_roll = 0.35 B is the usual estimate for
// a ro-pax, the added roll inertia is 25% (matching Ship::addedInertiaRoll), and
// GM = 2.0 m is a normal loaded condition for the type.
constexpr double kGm = 2.0;
constexpr double kAddedInertiaRoll = 0.25;

double rollInertia(const RollDampingHull& h) {
    const double mass = h.seaDensity * h.displacementVolume();
    return mass * (0.35 * h.beam) * (0.35 * h.beam) * (1.0 + kAddedInertiaRoll);
}

double rollStiffness(const RollDampingHull& h) {
    return h.seaDensity * kGravity * h.displacementVolume() * kGm;
}

double naturalRollFrequency(const RollDampingHull& h) {
    return std::sqrt(rollStiffness(h) / rollInertia(h));
}

// --- Nondimensionalisation ---------------------------------------------------

// B44hat = B44 / (rho * volume * B^2) * sqrt(B / (2 g)). Recomputed here from
// the definition rather than from the engine's own helper, so a unit slip in the
// helper cannot hide behind itself.
double nondimensional(const RollDampingHull& h, double b44) {
    const double volume = h.blockCoeff * h.lengthPp * h.beam * h.draft;
    return b44 / (h.seaDensity * volume * h.beam * h.beam) *
           std::sqrt(h.beam / (2.0 * kGravity));
}

void testNondimensionalisationMatchesItsDefinition() {
    const RollDampingHull h = roPax();
    const RollDamping d = rollDamping(h, at(10.0, naturalRollFrequency(h), 8.0));
    expectTrue("the test hull has non-trivial damping", d.total > 0);
    expectNear("B44hat matches its definition", d.totalHat, nondimensional(h, d.total),
               1e-12 * std::abs(d.totalHat));
    expectNear("the scale factor inverts the nondimensionalisation",
               d.totalHat * h.nondimensionalScale(), d.total, 1e-9 * d.total);
}

// Geometrically scale the hull by lambda and Froude-scale the motion with it
// (omega / sqrt(lambda), speed * sqrt(lambda)). Every component whose physics is
// inertia-dominated must then have an identical B44hat and a B44 exactly
// lambda^4.5 larger. Friction must not: it carries a Reynolds number and the
// viscosity does not scale, which is the whole reason model-scale roll decay
// tests need a friction correction. Getting the nondimensionalisation wrong by
// any power of a length shows up here immediately.
void testFroudeScaleInvariance() {
    const double lambda = 3.0;
    const RollDampingHull small = roPax();
    RollDampingHull large = small;
    large.lengthPp *= lambda;
    large.beam *= lambda;
    large.draft *= lambda;
    large.rollAxisAboveKeel *= lambda;
    large.bilgeKeelLength *= lambda;
    large.bilgeKeelBreadth *= lambda;

    const double omega = naturalRollFrequency(small);
    const double speed = 8.0;
    const RollDamping a = rollDamping(small, at(12.0, omega, speed));
    const RollDamping b =
        rollDamping(large, at(12.0, omega / std::sqrt(lambda), speed * std::sqrt(lambda)));

    const double expected = std::pow(lambda, 4.5);
    struct Case { const char* name; double small, large; };
    const Case inviscidlyScaled[] = {
        {"eddy", a.eddy, b.eddy},
        {"lift", a.lift, b.lift},
        {"bilge keel (normal force)", a.bilgeKeelNormal, b.bilgeKeelNormal},
        {"bilge keel (hull pressure)", a.bilgeKeelHull, b.bilgeKeelHull},
    };
    for (const Case& c : inviscidlyScaled) {
        expectTrue(std::string("scaling test: ") + c.name + " is non-zero", c.small > 0);
        expectNear(std::string(c.name) + " B44hat is Froude scale invariant",
                   nondimensional(large, c.large), nondimensional(small, c.small),
                   1e-11 * nondimensional(small, c.small));
        expectNear(std::string(c.name) + " B44 scales as lambda^4.5", c.large / c.small, expected,
                   1e-9 * expected);
    }

    // Friction: B44F is proportional to Sf * rf^2 * sqrt(omega), i.e. lambda^3.75,
    // so its nondimensional value falls as lambda^-0.75. Stated exactly, because
    // an approximate statement here would also pass for a wrong exponent.
    expectNear("friction B44 scales as lambda^3.75", b.friction / a.friction,
               std::pow(lambda, 3.75), 1e-9 * std::pow(lambda, 3.75));
    expectNear("friction B44hat falls as lambda^-0.75",
               nondimensional(large, b.friction) / nondimensional(small, a.friction),
               std::pow(lambda, -0.75), 1e-9);
    expectTrue("friction is relatively weaker on the larger hull",
               nondimensional(large, b.friction) < nondimensional(small, a.friction));

    // The total is a mixture, so it must sit strictly between the two behaviours.
    expectTrue("total B44hat is not scale invariant once friction is included",
               nondimensional(large, b.total) < nondimensional(small, a.total));
}

// --- Amplitude dependence ----------------------------------------------------

// The failure this guards against is a bilge-keel or eddy term that quietly
// keeps a constant floor at zero amplitude, or a friction term that vanishes
// with it. Both are invisible at 10 degrees and both are wrong.
void testAmplitudeDependence() {
    const RollDampingHull h = roPax();
    const double omega = naturalRollFrequency(h);

    const RollDamping zero = rollDamping(h, at(0.0, omega, 6.0));
    const RollDamping five = rollDamping(h, at(5.0, omega, 6.0));

    // Eddy damping is exactly proportional to roll amplitude: no vortex is shed
    // by a hull that is not moving.
    expectNear("eddy damping is exactly zero at zero roll amplitude", zero.eddy, 0.0, 0.0);
    const RollDamping ten = rollDamping(h, at(10.0, omega, 6.0));
    expectNear("eddy damping is exactly proportional to roll amplitude", ten.eddy,
               2.0 * five.eddy, 1e-9 * five.eddy);

    // Friction's amplitude cancels algebraically (Kato's Cf goes as 1/phi_a
    // against a leading phi_a), so it is not merely close at two amplitudes --
    // it is identical, and non-zero at zero amplitude.
    expectTrue("friction damping is non-zero at zero roll amplitude", zero.friction > 0);
    expectNear("friction damping is amplitude independent", ten.friction, zero.friction,
               1e-12 * zero.friction);
    expectTrue("lift damping is non-zero at zero roll amplitude", zero.lift > 0);
    expectNear("lift damping is amplitude independent", ten.lift, zero.lift, 1e-12 * zero.lift);

    // Bilge keels are the subtle case. Their drag coefficient rises as the
    // Keulegan-Carpenter number falls, exactly cancelling the leading phi_a, so
    // the equivalent linear damping tends to a *finite* non-zero limit as the
    // amplitude goes to zero -- a bilge keel still damps a small roll. What must
    // vanish is the amplitude-dependent part, and it does.
    // Stated per component, because the two parts share the term that produces
    // the limit and a sum would hide the loss of it from either one.
    expectTrue("bilge keel normal force has a finite non-zero zero-amplitude limit",
               zero.bilgeKeelNormal > 0);
    expectTrue("bilge keel hull pressure has a finite non-zero zero-amplitude limit",
               zero.bilgeKeelHull > 0);
    expectTrue("bilge keel damping grows above that limit with amplitude",
               ten.bilgeKeel() > zero.bilgeKeel() * 1.05);

    // The normal-force part is exactly affine in amplitude: CD = 22.5/Kc + 2.4
    // contributes a constant and a term proportional to phi_a, and nothing else
    // in it depends on amplitude at all.
    const RollDamping twenty = rollDamping(h, at(20.0, omega, 6.0));
    expectNear("bilge keel normal force is exactly affine in roll amplitude",
               twenty.bilgeKeelNormal - zero.bilgeKeelNormal,
               2.0 * (ten.bilgeKeelNormal - zero.bilgeKeelNormal),
               1e-9 * (ten.bilgeKeelNormal - zero.bilgeKeelNormal));

    // Monotone across a sweep, for every amplitude-dependent component.
    double prevEddy = -1, prevBkN = -1, prevBkH = -1;
    for (double deg = 0.0; deg <= 30.0 + 1e-9; deg += 1.0) {
        const RollDamping d = rollDamping(h, at(deg, omega, 6.0));
        expectTrue("eddy damping increases with roll amplitude", d.eddy > prevEddy);
        expectTrue("bilge keel normal force increases with roll amplitude",
                   d.bilgeKeelNormal > prevBkN);
        expectTrue("bilge keel hull pressure increases with roll amplitude",
                   d.bilgeKeelHull > prevBkH);
        prevEddy = d.eddy;
        prevBkN = d.bilgeKeelNormal;
        prevBkH = d.bilgeKeelHull;
    }
}

// The bilge-keel normal force is the one component that can be derived rather
// than fitted, so it is worth deriving. A pair of flat plates of area b_BK*l_BK
// at radius l from the roll axis, moving at l*phidot into a flow sped up by the
// factor f at the turn of the bilge, feels a quadratic drag moment
//   M = 2 * (1/2 rho CD b_BK l_BK (f l phidot)^2) * l = rho CD b_BK l_BK f^2 l^3
//       * phidot |phidot|
// and ITTC (3.6) linearises B2 phidot|phidot| to (8/3pi) omega phi_a B2. So
//   B44BKN = 8/(3 pi) rho CD b_BK l_BK f^2 l^3 omega phi_a
// which is ITTC (2.25) exactly. Only CD itself is empirical. Checking this pins
// the 8/3pi, the cube of the lever, the factor of two for the pair, and the
// square on f -- none of which the regression comparison can separate out.
void testBilgeKeelNormalForceIsLinearisedPlateDrag() {
    const RollDampingHull h = roPax();
    const double omega = naturalRollFrequency(h);

    // The section idealisation the method is built on, recomputed here from its
    // description: vertical side, horizontal bottom, quarter-circle bilge of
    // radius R, keel at the 45-degree point of the arc and normal to the shell.
    const double h0 = h.beam / (2.0 * h.draft);
    const double ogOverD = (h.draft - h.rollAxisAboveKeel) / h.draft;
    const double r = 2.0 * h.draft * std::sqrt(h0 * (h.midshipCoeff - 1.0) / (kPi - 4.0));
    const double inset = (1.0 - std::sqrt(2.0) / 2.0) * r / h.draft;
    const double lever = h.draft * std::sqrt((h0 - inset) * (h0 - inset) +
                                             (1.0 - ogOverD - inset) * (1.0 - ogOverD - inset));
    const double f = 1.0 + 0.3 * std::exp(-160.0 * (1.0 - h.midshipCoeff));
    expectTrue("the assumed bilge radius is a sane fraction of the draft",
               r > 0.1 * h.draft && r < h.draft);

    for (double deg : {3.0, 12.0, 28.0}) {
        const double phi = deg * kDegToRad;
        const double cd = 22.5 * h.bilgeKeelBreadth / (kPi * lever * f * phi) + 2.4;
        const double want = 8.0 / (3.0 * kPi) * h.seaDensity * cd * h.bilgeKeelBreadth *
                            h.bilgeKeelLength * f * f * lever * lever * lever * omega * phi;
        expectNear("bilge keel normal force is linearised plate drag at " +
                       std::to_string(static_cast<int>(deg)) + " deg",
                   rollDamping(h, at(deg, omega)).bilgeKeelNormal, want, 1e-9 * want);
    }
}

// The eddy moment is physically B2 * phidot * |phidot|. Equivalent linearisation
// of that (ITTC 3.6) gives Be = 8/(3 pi) * omega * phi_a * B2, so B44E divided by
// (omega * phi_a) must be one constant over the whole operating envelope. If it
// is not, the linearisation has been applied to the wrong power of something.
void testEddyIsALinearisedQuadraticMoment() {
    const RollDampingHull h = bareRoPax();
    double reference = -1;
    for (double omega : {0.3, 0.5236, 0.8, 1.2})
        for (double deg : {2.0, 7.5, 15.0, 25.0}) {
            const RollDamping d = rollDamping(h, at(deg, omega));
            const double b2 = d.eddy / (omega * deg * kDegToRad);
            if (reference < 0) reference = b2;
            expectNear("eddy damping / (omega * phi_a) is one constant", b2, reference,
                       1e-9 * reference);
        }
    expectTrue("that constant is positive", reference > 0);
}

// --- Bilge keels -------------------------------------------------------------

void testBilgeKeelPresence() {
    const RollDampingHull h = roPax();
    const double omega = naturalRollFrequency(h);
    const RollDamping with = rollDamping(h, at(12.0, omega, 6.0));
    const RollDamping without = rollDamping(bareRoPax(), at(12.0, omega, 6.0));

    expectNear("no bilge keels means exactly zero normal-force damping",
               without.bilgeKeelNormal, 0.0, 0.0);
    expectNear("no bilge keels means exactly zero hull-pressure damping", without.bilgeKeelHull,
               0.0, 0.0);
    expectTrue("bilge keels add damping", with.total > without.total);

    // Removing only one dimension must also zero it -- a keel of zero breadth or
    // zero length is not a keel.
    RollDampingHull noBreadth = h;
    noBreadth.bilgeKeelBreadth = 0;
    RollDampingHull noLength = h;
    noLength.bilgeKeelLength = 0;
    expectNear("zero bilge keel breadth gives zero bilge keel damping",
               rollDamping(noBreadth, at(12.0, omega)).bilgeKeel(), 0.0, 0.0);
    expectNear("zero bilge keel length gives zero bilge keel damping",
               rollDamping(noLength, at(12.0, omega)).bilgeKeel(), 0.0, 0.0);

    // Everything else must be untouched by the keels: the bare-hull components
    // are functions of the hull form alone.
    expectNear("bilge keels do not change the friction component", with.friction,
               without.friction, 1e-12 * without.friction);
    expectNear("bilge keels do not change the eddy component", with.eddy, without.eddy,
               1e-12 * without.eddy);

    // Larger keels, more damping, monotonically in both dimensions.
    double prev = 0;
    for (double breadth = 0.3; breadth <= 1.2 + 1e-9; breadth += 0.1) {
        RollDampingHull wider = h;
        wider.bilgeKeelBreadth = breadth;
        const double bk = rollDamping(wider, at(12.0, omega)).bilgeKeel();
        expectTrue("wider bilge keels damp more", bk > prev);
        prev = bk;
    }
    prev = 0;
    for (double length = 10.0; length <= 60.0 + 1e-9; length += 5.0) {
        RollDampingHull longer = h;
        longer.bilgeKeelLength = length;
        const double bk = rollDamping(longer, at(12.0, omega)).bilgeKeel();
        expectTrue("longer bilge keels damp more", bk > prev);
        prev = bk;
    }
}

// Kawahara, Maekawa & Ikeda (STAB 2009) eq. (38): a regression fitted to the
// output of Ikeda's sectional bilge-keel model. The engine implements the
// sectional model directly (ITTC 2.24-2.32); this is the independent second
// opinion. Reimplemented here from the paper rather than shared with the engine,
// which is the point -- a shared helper would agree with itself.
double simplifiedIkedaBilgeKeelHat(const RollDampingHull& h, double amplitudeDeg, double omega) {
    const double x1 = h.beam / h.draft;
    const double x2 = h.blockCoeff;
    const double x3 = h.midshipCoeff;
    const double x4 = (h.draft - h.rollAxisAboveKeel) / h.draft;
    const double x5 = omega * std::sqrt(h.beam / (2.0 * kGravity));
    const double x6 = amplitudeDeg;
    const double x7 = h.bilgeKeelBreadth / h.beam;
    const double x8 = h.bilgeKeelLength / h.lengthPp;

    const double f1 = (-0.3651 * x2 + 0.3907) * (x1 - 2.83) * (x1 - 2.83) - 2.21 * x2 + 2.632;
    const double f2 = 0.00255 * x6 * x6 + 0.122 * x6 + 0.4794;
    const double f3 = (-0.8913 * x7 * x7 - 0.0733 * x7) * x8 * x8 +
                      (5.2857 * x7 * x7 - 0.01185 * x7 + 0.00189) * x8;
    const double aBk = f1 * f2 * f3;
    const double bBk1 =
        (5.0 * x7 + 0.3 * x1 - 0.2 * x8 + 0.00125 * x6 * x6 - 0.0425 * x6 - 1.86) * x4;
    const double bBk2 = -15.0 * x7 + 1.2 * x2 - 0.1 * x1 - 0.0657 * x4 * x4 + 0.0586 * x4 + 1.6164;
    const double bBk3 = 2.5 * x4 + 15.75;
    return aBk * std::exp(bBk1 + bBk2 * std::pow(x3, bBk3)) * x5;
}

void testBilgeKeelAgreesWithTheSimplifiedRegression() {
    // **Swept over OG/d, because omega and amplitude cannot see the error this
    // comparison exists to catch.** The previous version of this test held the
    // hull at the roPax's own `OG/d = -0.6923` and swept omega and amplitude --
    // nine points, all of them on one hull. B0's first term is a function of
    // `m2 = OG/d` alone, so a wrong exponent on it is *constant* along both swept
    // axes, and -0.69 happens to sit within 0.05 of where the wrong reading and
    // the right one cross. The grid was nine points wide and one point deep, and
    // it passed at 3.6% against a `sqr` where a `cube` belongs.
    //
    // Swept, that same defect reads 15.6%. Measured worst over this grid:
    //
    //     OG/d      -1.50  -1.25  -1.00  -0.75  -0.50  -0.25   0.00   0.20
    //     as shipped  3.7%   1.9%   2.2%   3.1%   4.2%   4.9%   4.6%   4.0%
    //     with sqr   15.6%  11.3%   7.6%   4.3%   2.8%   4.6%   4.6%   4.2%
    //
    // The range is Ikeda's own, and it is taken from `validateRollDamping` rather
    // than written down twice: every hull below is asserted to be inside it and
    // one just outside is asserted to be refused, so narrowing the validator
    // without narrowing this sweep is a failure rather than a silent gap.
    const double kOgLo = -1.5, kOgHi = 0.2;
    auto atOg = [](double ogOverD) {
        RollDampingHull h = roPax();
        h.rollAxisAboveKeel = h.draft - ogOverD * h.draft;
        return h;
    };
    auto complainsAboutOg = [](const RollDampingHull& hull) {
        for (const std::string& p : validateRollDamping(hull, at(10.0, 0.5236)))
            if (p.find("OG/draft") != std::string::npos) return true;
        return false;
    };
    expectTrue("a hull just below Ikeda's OG/d range is refused", complainsAboutOg(atOg(-1.6)));
    expectTrue("and one just above it", complainsAboutOg(atOg(0.3)));

    double worst = 0;
    for (double ogOverD : {kOgLo, -1.25, -1.0, -0.75, -0.5, -0.25, 0.0, kOgHi}) {
        const RollDampingHull h = atOg(ogOverD);
        expectTrue("the swept hull is inside Ikeda's OG/d range", !complainsAboutOg(h));
        for (double omega : {0.4, 0.5236, 0.7})
            for (double deg : {5.0, 10.0, 20.0}) {
                const double got = nondimensional(h, rollDamping(h, at(deg, omega)).bilgeKeel());
                const double want = simplifiedIkedaBilgeKeelHat(h, deg, omega);
                expectTrue("the regression predicts positive bilge keel damping", want > 0);
                // 6% is a band the regression's own scatter fits inside while a
                // wrong sign, a wrong exponent or a factor of two anywhere in the
                // sectional pressure model does not -- 4.9% is the worst this
                // grid reaches, so there is one point of margin and not more.
                // Loosening it to the ~10% the regression's authors quote would
                // let a sign error through, which was checked by flipping one and
                // watching this test go green; it would also let `sqr(m2)` back.
                expectNear("sectional bilge keel damping matches Kawahara's regression at " +
                               std::to_string(static_cast<int>(deg)) + " deg, OG/d " +
                               std::to_string(ogOverD),
                           got, want, 0.06 * want);
                const double err = std::abs(got - want) / want;
                if (err > worst) worst = err;
            }
    }
    // Vacuity: two routes that agreed to the bit would be one route. The
    // regression is a fit to the sectional model's *output*, so some daylight is
    // the evidence they are independent implementations.
    expectTrue("the sweep found hulls the two routes do not agree exactly on", worst > 0.02);
}

// --- Forward speed -----------------------------------------------------------

void testForwardSpeedEffects() {
    const RollDampingHull h = roPax();
    const double omega = naturalRollFrequency(h);

    const RollDamping stopped = rollDamping(h, at(10.0, omega, 0.0));
    expectNear("lift damping is exactly zero at zero forward speed", stopped.lift, 0.0, 0.0);

    // Lift is a linear mechanism: exactly proportional to speed.
    const RollDamping slow = rollDamping(h, at(10.0, omega, 5.0));
    const RollDamping fast = rollDamping(h, at(10.0, omega, 10.0));
    expectTrue("lift damping is non-zero under way", slow.lift > 0);
    expectNear("lift damping is exactly proportional to speed", fast.lift, 2.0 * slow.lift,
               1e-9 * slow.lift);

    double prevLift = -1, prevFriction = -1, prevEddy = 1e300;
    for (double u = 0.0; u <= 12.0 + 1e-9; u += 0.5) {
        const RollDamping d = rollDamping(h, at(10.0, omega, u));
        expectTrue("lift damping increases with speed", d.lift > prevLift);
        expectTrue("friction damping increases with speed", d.friction > prevFriction);
        // ITTC (2.21): the bilge vortices are swept downstream, so the eddy
        // component collapses with speed. Monotone, and it must reach zero in
        // the limit rather than levelling off at a floor.
        expectTrue("eddy damping falls with speed", d.eddy < prevEddy);
        prevLift = d.lift;
        prevFriction = d.friction;
        prevEddy = d.eddy;
    }
    expectTrue("eddy damping is nearly gone at high speed",
               rollDamping(h, at(10.0, omega, 400.0)).eddy < 0.01 * stopped.eddy);

    // The bilge keel component is taken as speed independent (ITTC 2.2.5.1), and
    // must be exactly so rather than accidentally drifting.
    expectNear("bilge keel damping is speed independent", fast.bilgeKeel(), stopped.bilgeKeel(),
               1e-12 * stopped.bilgeKeel());
}

// --- What the coefficient means ---------------------------------------------

// The proof that B44 is a linear damping coefficient in N m s/rad and not
// something a factor of two or a radian-per-degree away from one: put it into
//   (I + A) phiddot + B44 phidot + C phi = 0
// integrate a free decay, and compare the measured logarithmic decrement against
// the analytic value for a linear damped oscillator. Nothing about Ikeda is
// being tested here -- only that the number he produces means what it claims.
void testRollDecayLogarithmicDecrement() {
    const RollDampingHull h = roPax();
    const double inertia = rollInertia(h);
    const double stiffness = rollStiffness(h);
    const double omegaN = std::sqrt(stiffness / inertia);
    const double b44 = rollDamping(h, at(10.0, omegaN, 6.0)).total;

    const double zeta = b44 / (2.0 * std::sqrt(inertia * stiffness));
    expectTrue("the damping ratio is light, as roll always is", zeta > 0.01 && zeta < 0.3);
    const double wanted = 2.0 * kPi * zeta / std::sqrt(1.0 - zeta * zeta);

    // RK4 on the linear oscillator. B44 is held fixed, which is exactly the
    // linearisation the coefficient represents.
    const double dt = 1e-4;
    double phi = 10.0 * kDegToRad, rate = 0.0;
    auto accel = [&](double p, double v) { return -(b44 * v + stiffness * p) / inertia; };

    std::vector<double> peaks;
    double y0 = phi, y1 = phi;
    const int steps = static_cast<int>(8.0 * 2.0 * kPi / omegaN / dt);
    for (int i = 0; i < steps; ++i) {
        const double k1p = rate,                  k1v = accel(phi, rate);
        const double k2p = rate + 0.5 * dt * k1v, k2v = accel(phi + 0.5 * dt * k1p, k2p);
        const double k3p = rate + 0.5 * dt * k2v, k3v = accel(phi + 0.5 * dt * k2p, k3p);
        const double k4p = rate + dt * k3v,       k4v = accel(phi + dt * k3p, k4p);
        phi += dt / 6.0 * (k1p + 2 * k2p + 2 * k3p + k4p);
        rate += dt / 6.0 * (k1v + 2 * k2v + 2 * k3v + k4v);

        // Parabolic interpolation of the maximum through three samples; sampling
        // the peak at the nearest step would bias every peak downward and quietly
        // inflate the measured decrement.
        if (y1 > y0 && y1 > phi) {
            const double denom = y0 - 2.0 * y1 + phi;
            const double shift = denom != 0 ? 0.5 * (y0 - phi) / denom : 0.0;
            peaks.push_back(y1 - 0.25 * (y0 - phi) * shift);
        }
        y0 = y1;
        y1 = phi;
    }

    expectTrue("the free decay produced several peaks", peaks.size() >= 5);
    for (std::size_t i = 1; i < peaks.size(); ++i) {
        expectTrue("the decay is monotone", peaks[i] < peaks[i - 1]);
        expectNear("logarithmic decrement matches 2 pi zeta / sqrt(1 - zeta^2)",
                   std::log(peaks[i - 1] / peaks[i]), wanted, 1e-6);
    }

    // And the same statement the other way round: the decay envelope is
    // exp(-zeta * omega_n * t), so after n periods the amplitude is down by
    // exp(-n * delta). Checked against the last peak rather than the first pair.
    const double n = static_cast<double>(peaks.size() - 1);
    expectNear("the decay envelope over the whole record matches the same zeta",
               peaks.back(), peaks.front() * std::exp(-n * wanted), 1e-9 * peaks.front());
}

// --- Magnitudes --------------------------------------------------------------

// A method can pass every structural test above and still be out by an order of
// magnitude. This is the check on absolute size, stated two ways: as B44hat at
// the operating point the literature quotes it at, and as a fraction of critical
// damping at the hull's own natural roll frequency, which is the number that
// decides how far a ship rolls.
void testMagnitudesAreShipLike() {
    const RollDampingHull bare = bareRoPax();
    const RollDampingHull keeled = roPax();

    // Ikeda's published figures are drawn at omegahat = omega sqrt(B/2g) around
    // 1.1 and roll amplitudes around 20 degrees; quoted at any other operating
    // point B44hat is a different number, because it is a strong function of
    // both. This is the comparison point for the literature values.
    //
    // Those literature values are for the *total* B44, radiation included, and
    // this engine computes only the viscous part -- the wave component is an
    // input from the BEM pipeline that does not exist yet. ITTC 7.5-02-07-04.5
    // puts the wave component at 5-30% of the total for a conventional hull, so
    // the viscous-only share of a total in the customary 0.01-0.05 band is
    // 0.70-0.95 of it, i.e. 0.007-0.0475. The bands below are that arithmetic
    // and nothing else; they are not fitted to the answer.
    const double omegaLit = 1.1 / std::sqrt(keeled.beam / (2.0 * kGravity));
    const double bareHat = rollDamping(bare, at(20.0, omegaLit)).totalHat;
    const double keeledHat = rollDamping(keeled, at(20.0, omegaLit)).totalHat;
    expectTrue("bare hull viscous B44hat is 0.70-0.95 of the 0.01-0.05 total band",
               bareHat > 0.01 * 0.70 && bareHat < 0.05 * 0.95);
    expectTrue("keeled viscous B44hat is 0.70-0.95 of the 0.03-0.10 total band",
               keeledHat > 0.03 * 0.70 && keeledHat < 0.10 * 0.95);
    expectTrue("bilge keels are the dominant contribution", keeledHat > 2.0 * bareHat);

    // The same statement made directly: supplying a wave component at the middle
    // of ITTC's documented share must land the totals inside the literature
    // bands as published.
    RollDampingHull bareWithWave = bare, keeledWithWave = keeled;
    bareWithWave.waveDamping = 0.25 / 0.75 * rollDamping(bare, at(20.0, omegaLit)).total;
    keeledWithWave.waveDamping = 0.25 / 0.75 * rollDamping(keeled, at(20.0, omegaLit)).total;
    expectTrue("bare hull total B44hat is in the 0.01-0.05 band with a 25% wave share",
               rollDamping(bareWithWave, at(20.0, omegaLit)).totalHat > 0.01 &&
                   rollDamping(bareWithWave, at(20.0, omegaLit)).totalHat < 0.05);
    expectTrue("keeled total B44hat is in the 0.03-0.10 band with a 25% wave share",
               rollDamping(keeledWithWave, at(20.0, omegaLit)).totalHat > 0.03 &&
                   rollDamping(keeledWithWave, at(20.0, omegaLit)).totalHat < 0.10);

    // At the hull's own natural roll frequency the coefficient is smaller,
    // because omegahat there is only about 0.6. The physically meaningful
    // statement is the damping ratio, and a ro-pax with bilge keels is
    // conventionally 5-10% of critical -- which is where Ship::zetaRoll sits.
    const double omegaN = naturalRollFrequency(keeled);
    const double critical = 2.0 * std::sqrt(rollInertia(keeled) * rollStiffness(keeled));
    const double zetaBare = rollDamping(bare, at(10.0, omegaN)).total / critical;
    const double zetaKeeled = rollDamping(keeled, at(10.0, omegaN)).total / critical;
    expectTrue("bare hull damping is 1-4% of critical at the natural roll frequency",
               zetaBare > 0.01 && zetaBare < 0.04);
    expectTrue("with bilge keels it is 5-12% of critical",
               zetaKeeled > 0.05 && zetaKeeled < 0.12);

    // Component ranking. Bilge keels are 50-80% of the total for a hull that has
    // them (Kawahara et al. 4.4), and the eddy component is what is left of a
    // bare hull once friction is discounted.
    const RollDamping d = rollDamping(keeled, at(10.0, omegaN, 6.0));
    expectTrue("bilge keels are 50-80% of the total for a keeled hull",
               d.bilgeKeel() > 0.5 * d.total && d.bilgeKeel() < 0.9 * d.total);
    const RollDamping bareAtRest = rollDamping(bare, at(10.0, omegaN));
    expectTrue("eddy damping dominates a bare hull at rest",
               bareAtRest.eddy > 0.9 * bareAtRest.total);

    // The absolute check on the Reynolds scaling, not just its exponent. ITTC
    // puts friction at 1-3% of the roll damping of a full-scale ship and 8-10%
    // for a 2 m model of the same hull -- the only component with a scale effect
    // at all. Reproducing both ends from one formula is a real constraint.
    // The bands here are widened by the missing radiation component, which would
    // enlarge both totals and so reduce both shares.
    expectTrue("friction is 1-3% of a full-scale bare hull's damping",
               bareAtRest.friction > 0.01 * bareAtRest.total &&
                   bareAtRest.friction < 0.03 * bareAtRest.total);

    const double lambda = 2.0 / keeled.lengthPp;  // a 2 m towing-tank model
    RollDampingHull model = keeled;
    model.lengthPp *= lambda;
    model.beam *= lambda;
    model.draft *= lambda;
    model.rollAxisAboveKeel *= lambda;
    model.bilgeKeelLength *= lambda;
    model.bilgeKeelBreadth *= lambda;
    model.seaDensity = kRhoFresh;  // a towing tank is fresh water
    const RollDamping m = rollDamping(model, at(10.0, omegaN / std::sqrt(lambda)));
    expectTrue("friction is 5-15% for a 2 m model, bracketing ITTC's quoted 8-10%",
               m.friction > 0.05 * m.total && m.friction < 0.15 * m.total);
    expectTrue("the model's friction share is an order of magnitude above the ship's",
               m.friction / m.total > 8.0 * (d.friction / d.total));
    std::printf("     friction share: ship %.1f%% of bare-hull damping, 2 m model %.1f%%"
                " of total\n",
                100 * bareAtRest.friction / bareAtRest.total, 100 * m.friction / m.total);

    std::printf("     ro-pax 170x25x6.5, Cb 0.55, Cm 0.98, KG 11.0, keels 34.0 x 0.60 m\n");
    std::printf("     omega_n %.4f rad/s (T %.1f s), omegahat %.3f;  bare %.5f  keeled %.5f"
                "  at omegahat 1.1 / 20 deg\n",
                omegaN, 2.0 * kPi / omegaN, omegaN * std::sqrt(keeled.beam / (2.0 * kGravity)),
                bareHat, keeledHat);
    for (double speed : {0.0, 6.0})
        for (double deg : {2.5, 10.0, 25.0}) {
            const RollDamping r = rollDamping(keeled, at(deg, omegaN, speed));
            std::printf("     %4.1f m/s %4.1f deg: B44hat %.5f  zeta %.3f  fric %4.1f%%"
                        "  eddy %4.1f%%  lift %4.1f%%  keels %4.1f%%\n",
                        speed, deg, r.totalHat, r.total / critical, 100 * r.friction / r.total,
                        100 * r.eddy / r.total, 100 * r.lift / r.total,
                        100 * r.bilgeKeel() / r.total);
        }
}

// --- Validity ----------------------------------------------------------------

// Ikeda's method has a documented domain, and silently extrapolating past it is
// how an empirical model produces confident nonsense. The test hull must be
// inside it; a hull outside must say so.
void testValidityIsReported() {
    RollDampingHull h = roPax();
    h.waveDamping = 1.0e7;  // otherwise the missing radiation term is itself a note
    const RollDampingCondition c = at(10.0, naturalRollFrequency(h), 6.0);
    const std::vector<std::string> clean = validateRollDamping(h, c);
    for (const std::string& s : clean) std::printf("     unexpected: %s\n", s.c_str());
    expectTrue("the ro-pax test hull is inside Ikeda's validity range", clean.empty());

    expectTrue("a missing wave component is reported",
               !validateRollDamping(roPax(), c).empty());

    RollDampingHull barge = h;
    barge.draft = 3.0;   // B/d = 8.3, far outside 2.5-4.5
    barge.blockCoeff = 0.92;
    expectTrue("a flat full barge is reported as out of range",
               validateRollDamping(barge, c).size() >= 2);

    RollDampingHull topHeavy = h;
    topHeavy.rollAxisAboveKeel = 30.0;  // OG/d = -3.6
    expectTrue("a roll axis far above the waterline is reported as out of range",
               !validateRollDamping(topHeavy, c).empty());

    RollDampingHull broken;
    expectTrue("an empty hull is rejected outright", !validateRollDamping(broken, c).empty());
    expectNear("an empty hull produces no damping rather than a NaN",
               rollDamping(broken, c).total, 0.0, 0.0);
}

// Two implementations of the bilge radius, and nothing had ever set them side by
// side. `hullform.cpp`'s `bilgeRadiusForMidshipCoefficient` is the one with a test
// -- a round trip against its own inverse. `RollDampingHull::bilgeRadiusOrDefault`
// is the one the ship actually runs: `rollDampingHullFromMesh` leaves `bilgeRadius`
// at its `-1` sentinel, so the derived branch is what production takes, and
// `derive()` feeds the answer straight into `Geometry::bilgeR`, on which both
// bilge-keel components are built. It had zero mentions anywhere under `tests/`.
//
// They agree, and by algebra rather than by coincidence:
//
//     hullform      r = sqrt((1-Cm) B d / (2 (1 - pi/4)))  =  sqrt(2 B d (1-Cm) / (4-pi))
//     roll damping  r = 2d sqrt(h0 (Cm-1) / (pi-4)),  h0 = B/2d
//                     = sqrt(4 d^2 B (1-Cm) / (2 d (4-pi)))  =  sqrt(2 B d (1-Cm) / (4-pi))
//
// and the clamps coincide too, which is the part that looks like a difference and
// is not: `h0 >= 1` is `B >= 2d` is `0.5B >= d`, so `h0 >= 1 ? d : 0.5B` *is*
// `min(0.5B, d)`, written the other way round.
//
// None of that is written down anywhere, and neither routine's header mentions the
// other. This is the repo's own rule applied where it had not been: when a second
// way to compute something exists, assert they agree rather than trusting that they
// do.
void testTheTwoBilgeRadiiAgree() {
    double worst = 0, largest = 0;
    int clamped = 0, free_ = 0;
    for (double beam : {8.0, 14.0, 20.0, 32.0}) {
        for (double draft : {2.0, 4.0, 6.5, 11.0}) {
            for (double cm : {0.55, 0.70, 0.85, 0.95, 0.99}) {
                RollDampingHull hull;
                hull.beam = beam;
                hull.draft = draft;
                hull.midshipCoeff = cm;
                // Left at the -1 sentinel deliberately: that is the branch the ship
                // takes, and pinning the other one would test nothing.
                const double derived = hull.bilgeRadiusOrDefault();
                const double independent = bilgeRadiusForMidshipCoefficient(beam, draft, cm);
                worst = std::max(worst, std::fabs(derived - independent));
                largest = std::max(largest, derived);
                // Did the section-folding clamp bind? Both apply it, and a sweep
                // that never reached it would leave the halves of it untested --
                // which is where the two are written most differently.
                if (derived >= std::min(0.5 * beam, draft) - 1e-12) ++clamped; else ++free_;
            }
        }
    }
    std::printf("     the two bilge radii: worst disagreement %.3e m over 80 hulls,"
                " largest radius %.3f m\n", worst, largest);
    std::printf("     the fold clamp bound on %d of them and not on %d\n", clamped, free_);

    // Asserted at what was measured. They are the same expression, so the only
    // difference available is the order of the floating-point operations.
    expectTrue("the roll-damping default and the hull-form estimator agree",
               worst <= 1e-12 * std::max(1.0, largest));
    // Guards: a sweep that clamped everywhere would compare two `min` calls and
    // nothing else, and one that never clamped would leave the clamp -- the half
    // written differently in the two files -- unexercised.
    expectTrue("and the sweep drove both sides of the fold clamp", clamped > 0 && free_ > 0);
    expectTrue("and there was a real radius to compare", largest > 1.0);

    // The divergence, stated rather than merely known. `bilgeRadiusOrDefault`
    // refuses a midship coefficient outside (0, 1) and returns zero;
    // `bilgeRadiusForMidshipCoefficient` clamps the coefficient into range and
    // answers anyway. Both are defensible and no live caller passes such a hull --
    // `RollDampingHull` is built from a mesh, and `validateRollDamping` rejects
    // these before Ikeda sees them -- but the next caller to try one should find
    // out here rather than from a zero that looks like a very fine bilge.
    RollDampingHull square;
    square.beam = 20.0;
    square.draft = 6.0;
    square.midshipCoeff = 1.0;  // a box section: no bilge radius at all
    expectNear("a box section has no bilge radius, both ways",
               square.bilgeRadiusOrDefault(), 0.0, 1e-12);
    expectNear("and so it does by the hull-form route",
               bilgeRadiusForMidshipCoefficient(20.0, 6.0, 1.0), 0.0, 1e-12);
    square.midshipCoeff = 0.0;  // nonphysical: the section has no area at all
    expectNear("at a nonphysical Cm the roll-damping route refuses and returns zero",
               square.bilgeRadiusOrDefault(), 0.0, 1e-12);
    expectTrue("where the hull-form route clamps the coefficient and answers",
               bilgeRadiusForMidshipCoefficient(20.0, 6.0, 0.0) > 1.0);
}

}  // namespace

void runRollDampingTests() {
    std::printf("\n--- viscous roll damping (Ikeda) ---\n");
    testNondimensionalisationMatchesItsDefinition();
    testFroudeScaleInvariance();
    testAmplitudeDependence();
    testBilgeKeelNormalForceIsLinearisedPlateDrag();
    testEddyIsALinearisedQuadraticMoment();
    testBilgeKeelPresence();
    testTheTwoBilgeRadiiAgree();
    testBilgeKeelAgreesWithTheSimplifiedRegression();
    testForwardSpeedEffects();
    testRollDecayLogarithmicDecrement();
    testMagnitudesAreShipLike();
    testValidityIsReported();
}
