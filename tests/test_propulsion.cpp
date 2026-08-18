// SPDX-License-Identifier: MIT
//
// Validation of propulsion, steering and manoeuvring.
//
// A propeller model that produces a plausible-looking thrust curve can still be
// wrong by a factor of D, by a factor of 2 pi, or in the sign of astern thrust,
// and none of those show up as anything but "the ship feels a bit fast". So
// nothing here is eyeballed. Every check is either a dimensional identity that
// must hold exactly (T = rho n^2 D^4 K_T under independent scaling of n, D and
// rho), a physical bound that cannot be exceeded without a factor error
// (open-water efficiency <= 1), an exact symmetry (a port turn is the mirror of
// a starboard turn), or a closed form recomputed here from the parameters.
#include "engine/sim/propulsion.hpp"
#include "harness.hpp"

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

using namespace sim;
using testing::expectNear;
using testing::expectTrue;

namespace {

// --- Propeller: dimensional identities --------------------------------------

// K_T and K_Q are defined by T = rho n^2 D^4 K_T and Q = rho n^2 D^5 K_Q. The
// two code paths -- the J-parameterised coefficient and the four-quadrant
// evaluation -- must agree, and holding J fixed while scaling n, D and rho
// independently must scale thrust and torque by exactly the powers above. This
// is the check that catches a ship with ten times too much thrust.
void testDimensionalIdentities() {
    PropellerParams base;
    base.diameter = 5.0;
    base.wakeFraction = 0.0;  // so surge speed is the advance speed

    const double j = 0.5;
    const double n = 2.0;
    const double rho = kRhoSeawater;
    const double va = j * n * base.diameter;

    const PropellerState s = evaluatePropeller(base, va, n, rho);
    expectNear("advance ratio round trip", s.advanceRatio, j, 1e-12);

    const double kt = propellerThrustCoefficient(base, j);
    const double kq = propellerTorqueCoefficient(base, j);
    const double d4 = std::pow(base.diameter, 4.0);
    expectNear("T = rho n^2 D^4 K_T", s.thrust, rho * n * n * d4 * kt, 1e-6 * std::abs(s.thrust));
    expectNear("Q = rho n^2 D^5 K_Q", s.torque,
               rho * n * n * d4 * base.diameter * kq, 1e-6 * std::abs(s.torque));
    expectNear("reported K_T matches the J-parameterised curve", s.thrustCoefficient, kt, 1e-12);
    expectNear("reported K_Q matches the J-parameterised curve", s.torqueCoefficient, kq, 1e-12);

    // Scale n by 3 at fixed J: thrust x 9, torque x 9, power x 27.
    const PropellerState fastShaft = evaluatePropeller(base, 3.0 * va, 3.0 * n, rho);
    expectNear("thrust scales as n^2", fastShaft.thrust, 9.0 * s.thrust, 1e-6 * 9.0 * s.thrust);
    expectNear("torque scales as n^2", fastShaft.torque, 9.0 * s.torque, 1e-6 * 9.0 * s.torque);
    expectNear("delivered power scales as n^3", fastShaft.deliveredPower,
               27.0 * s.deliveredPower, 1e-6 * 27.0 * std::abs(s.deliveredPower));

    // Scale D by 2 at fixed J and n: thrust x 16, torque x 32.
    PropellerParams big = base;
    big.diameter = 2.0 * base.diameter;
    const PropellerState bigProp = evaluatePropeller(big, 2.0 * va, n, rho);
    expectNear("thrust scales as D^4", bigProp.thrust, 16.0 * s.thrust, 1e-6 * 16.0 * s.thrust);
    expectNear("torque scales as D^5", bigProp.torque, 32.0 * s.torque, 1e-6 * 32.0 * s.torque);
    expectNear("K_T is invariant under a pure size change", bigProp.thrustCoefficient, kt, 1e-12);
    expectNear("K_Q is invariant under a pure size change", bigProp.torqueCoefficient, kq, 1e-12);

    // Scale rho: linear in both, coefficients untouched.
    const PropellerState fresh = evaluatePropeller(base, va, n, kRhoFresh);
    const double ratio = kRhoFresh / kRhoSeawater;
    expectNear("thrust is linear in density", fresh.thrust, ratio * s.thrust,
               1e-9 * std::abs(s.thrust));
    expectNear("torque is linear in density", fresh.torque, ratio * s.torque,
               1e-9 * std::abs(s.torque));
    expectNear("K_T is invariant under a density change", fresh.thrustCoefficient, kt, 1e-12);

    // Wake fraction and thrust deduction, by definition.
    PropellerParams wakeful = base;
    wakeful.wakeFraction = 0.30;
    wakeful.thrustDeduction = 0.18;
    const PropellerState w = evaluatePropeller(wakeful, 10.0, n, rho);
    expectNear("Va = u (1 - w)", w.advanceSpeed, 10.0 * 0.70, 1e-12);
    expectNear("hull feels (1 - t) T", w.thrustOnHull, 0.82 * w.thrust,
               1e-9 * std::abs(w.thrust));
}

// Bisect for the advance ratio at which K_T changes sign.
double zeroThrustAdvanceRatio(const PropellerParams& p) {
    double lo = 0.0, hi = 3.0;
    for (int i = 0; i < 200; ++i) {
        const double mid = 0.5 * (lo + hi);
        if (propellerThrustCoefficient(p, mid) > 0.0) lo = mid; else hi = mid;
    }
    return 0.5 * (lo + hi);
}

// Bollard pull is the most thrust a given shaft speed can make, K_T must fall
// away monotonically from it, and it must reach zero at a finite advance ratio
// -- otherwise the ship has no top speed.
void testBollardPullAndMonotonicity() {
    const PropellerParams p;  // P/D = 1.0, A_E/A_0 = 0.70

    const double bollard = propellerThrustCoefficient(p, 0.0);
    expectTrue("bollard K_T is positive", bollard > 0.0);

    double previous = bollard;
    bool monotone = true, everExceeded = false;
    for (int i = 1; i <= 300; ++i) {
        const double j = i * 0.005;
        const double kt = propellerThrustCoefficient(p, j);
        if (kt >= previous) monotone = false;
        if (kt > bollard) everExceeded = true;
        previous = kt;
    }
    expectTrue("K_T decreases monotonically over 0 <= J <= 1.5", monotone);
    expectTrue("no J beats bollard pull", !everExceeded);

    const double j0 = zeroThrustAdvanceRatio(p);
    expectTrue("K_T crosses zero at a finite J", j0 > 0.4 && j0 < 1.5);
    expectNear("K_T is zero at the crossing", propellerThrustCoefficient(p, j0), 0.0, 1e-9);

    // Calibration targets: the published Wageningen B4-70 curve at P/D = 1.0.
    // These are the numbers the free constants of the blade-element model were
    // fitted to, so this is a regression fence on the calibration, not an
    // independent validation.
    expectNear("bollard K_T matches the B4-70 magnitude", bollard, 0.35, 0.03);
    expectNear("bollard 10 K_Q matches the B4-70 magnitude",
               10.0 * propellerTorqueCoefficient(p, 0.0), 0.49, 0.04);
    expectNear("zero-thrust J matches the B4-70 value", j0, 0.85, 0.05);

    // **Printed as well as bracketed.** The three assertions above are against the
    // *published* B4-70 values with the band the calibration was fitted to; the
    // model's own figures -- which is what section 7's table publishes -- were not
    // written anywhere a reader or a gate could see them. A bracket around a
    // reference is not a check that the model still produces what the document
    // says it produces.
    std::printf("     bollard K_T %.3f, 10 K_Q %.3f, zero thrust at J %.3f\n", bollard,
                10.0 * propellerTorqueCoefficient(p, 0.0), j0);

    // More pitch must buy more top-end and more bollard torque.
    PropellerParams coarse = p;
    coarse.pitchRatio = 1.3;
    expectTrue("coarser pitch pushes the zero-thrust J out",
               zeroThrustAdvanceRatio(coarse) > j0 + 0.2);
    expectTrue("coarser pitch costs more torque at bollard",
               propellerTorqueCoefficient(coarse, 0.0) > propellerTorqueCoefficient(p, 0.0));
}

// eta = J K_T / (2 pi K_Q) is bounded above by 1 for any real propeller: it is
// useful work out over shaft work in. A model that exceeds it has a sign or a
// factor wrong somewhere, and nothing else about the curve would look odd.
void testOpenWaterEfficiency() {
    const PropellerParams p;
    const double j0 = zeroThrustAdvanceRatio(p);

    expectNear("efficiency is zero at bollard pull", propellerOpenWaterEfficiency(p, 0.0), 0.0,
               1e-15);
    expectNear("efficiency is zero again where K_T crosses zero",
               propellerOpenWaterEfficiency(p, j0), 0.0, 1e-9);

    double peak = 0.0, peakJ = 0.0, worstIdentity = 0.0;
    for (int i = 1; i < 1000; ++i) {
        const double j = j0 * i / 1000.0;
        const double eta = propellerOpenWaterEfficiency(p, j);
        if (eta > peak) { peak = eta; peakJ = j; }
        // The identity that guarantees the bound: eta = J K_T / (2 pi K_Q).
        const double kt = propellerThrustCoefficient(p, j);
        const double kq = propellerTorqueCoefficient(p, j);
        worstIdentity =
            std::max(worstIdentity, std::abs(eta - j * kt / (2.0 * kPi * kq)));
    }
    expectNear("efficiency equals J K_T / (2 pi K_Q) everywhere", worstIdentity, 0.0, 1e-12);
    expectTrue("efficiency never exceeds 1", peak < 1.0);
    expectTrue("efficiency peaks strictly inside the working range",
               peakJ > 0.1 * j0 && peakJ < 0.98 * j0);
    expectTrue("peak efficiency is physically plausible", peak > 0.55 && peak < 0.80);
    std::printf("     open-water peak eta = %.4f at J = %.3f (zero thrust at J = %.3f)\n", peak,
                peakJ, j0);
}

// Delivered power is 2 pi n Q, thrust power is T Va, and their ratio is the
// efficiency. Asserting all three separately is what makes a stray 2 pi visible.
void testPowerConsistency() {
    PropellerParams p;
    p.wakeFraction = 0.0;
    const double n = 1.7, rho = kRhoSeawater;

    for (double j : {0.1, 0.3, 0.5, 0.7}) {
        const double va = j * n * p.diameter;
        const PropellerState s = evaluatePropeller(p, va, n, rho);
        const std::string at = " @ J=" + std::to_string(j);
        expectNear("delivered power = 2 pi n Q" + at, s.deliveredPower,
                   2.0 * kPi * n * s.torque, 1e-9 * std::abs(s.deliveredPower));
        expectNear("thrust power = T Va" + at, s.thrustPower, s.thrust * va,
                   1e-9 * std::abs(s.thrustPower));
        expectNear("power ratio is the open-water efficiency" + at,
                   s.thrustPower / s.deliveredPower, propellerOpenWaterEfficiency(p, j), 1e-12);
        expectNear("reported efficiency agrees" + at, s.efficiency,
                   propellerOpenWaterEfficiency(p, j), 1e-12);
        expectTrue("delivered power exceeds thrust power" + at,
                   s.deliveredPower > s.thrustPower);
    }
}

// Four quadrants. The one that matters operationally is astern rotation with
// ahead flow: if that produced forward thrust, a crash stop would accelerate.
void testFourQuadrantSigns() {
    PropellerParams p;
    p.wakeFraction = 0.0;

    const PropellerState dead = evaluatePropeller(p, 0.0, 0.0);
    expectNear("stopped shaft, dead in the water: no thrust", dead.thrust, 0.0, 0.0);
    expectNear("stopped shaft, dead in the water: no torque", dead.torque, 0.0, 0.0);

    // Zero advance: reversing rotation reverses thrust and torque exactly. The
    // blade section sees the mirrored angle of attack, so this is an algebraic
    // identity of the model rather than a tuned result.
    for (double n : {0.5, 1.0, 2.5}) {
        const PropellerState ahead = evaluatePropeller(p, 0.0, n);
        const PropellerState astern = evaluatePropeller(p, 0.0, -n);
        expectNear("reversing rotation reverses thrust at zero advance", astern.thrust,
                   -ahead.thrust, 1e-9 * std::abs(ahead.thrust));
        expectNear("reversing rotation reverses torque at zero advance", astern.torque,
                   -ahead.torque, 1e-9 * std::abs(ahead.torque));
        expectTrue("ahead rotation at zero advance pushes forward", ahead.thrust > 0.0);
    }

    // Crash stop: making way ahead, shaft reversed. Thrust must be astern at
    // every speed and every astern shaft speed, with no window of positive
    // thrust hiding between quadrants.
    bool everForward = false, everDriving = false;
    for (int i = 1; i <= 40; ++i) {
        const double u = i * 0.5;  // 0.5 .. 20 m/s
        for (double n : {-0.25, -1.0, -2.0, -4.0}) {
            const PropellerState s = evaluatePropeller(p, u, n);
            if (s.thrust > 0.0) everForward = true;
            if (s.torque > 0.0) everDriving = true;
        }
    }
    expectTrue("a crash stop never produces forward thrust", !everForward);
    expectTrue("astern rotation always needs astern torque", !everDriving);

    // Backing down with the shaft turning ahead (the other off-design quadrant)
    // must still push the ship forward.
    for (int i = 1; i <= 20; ++i)
        expectTrue("ahead rotation in astern flow still thrusts ahead",
                   evaluatePropeller(p, -i * 0.5, 2.0).thrust > 0.0);

    // A stopped shaft in a stream is a drag disc, not a neutral object: it must
    // resist, never assist. Modelling this as exactly zero -- which the bare
    // T = rho n^2 D^4 K_T form does -- under-predicts crash-stop deceleration,
    // so the four-quadrant model deliberately does not.
    const PropellerState locked = evaluatePropeller(p, 8.0, 0.0);
    expectTrue("a locked propeller in ahead flow drags", locked.thrust < 0.0);
    expectTrue("locked-rotor drag is a sane fraction of bollard thrust",
               std::abs(locked.thrust) < evaluatePropeller(p, 0.0, 2.0).thrust);
    expectNear("locked-rotor drag scales as speed squared",
               evaluatePropeller(p, 16.0, 0.0).thrust, 4.0 * locked.thrust,
               1e-9 * 4.0 * std::abs(locked.thrust));
}

// --- Rudder -----------------------------------------------------------------

double fujiiSlope(const RudderParams& r) { return 6.13 * r.aspectRatio / (r.aspectRatio + 2.25); }

void testRudderNormalCoefficient() {
    const RudderParams r = kvlcc2Rudder();
    const double fAlpha = fujiiSlope(r);

    expectNear("no angle, no normal force", rudderNormalCoefficient(r, 0.0), 0.0, 0.0);
    expectNear("C_N is odd in the angle of attack", rudderNormalCoefficient(r, -0.2),
               -rudderNormalCoefficient(r, 0.2), 1e-15);
    expectNear("flow straight along the chord from astern makes no force",
               rudderNormalCoefficient(r, kPi), 0.0, 1e-12);

    // Below stall the model is exactly Fujii's f_alpha sin(alpha), so it is
    // linear in the angle to within the sine's own curvature.
    for (double deg : {1.0, 5.0, 10.0, 20.0}) {
        const double a = deg * kDegToRad;
        expectNear("Fujii law holds below stall @ " + std::to_string(deg) + " deg",
                   rudderNormalCoefficient(r, a), fAlpha * std::sin(a), 1e-15);
    }
    const double small = 0.5 * kDegToRad;
    expectNear("C_N is linear in alpha at small angles",
               rudderNormalCoefficient(r, 2.0 * small) / rudderNormalCoefficient(r, small), 2.0,
               1e-4);

    // Continuity across the stall angle: the two branches must meet.
    const double eps = 1e-7;
    expectNear("C_N is continuous at the stall angle",
               rudderNormalCoefficient(r, r.stallAngle + eps),
               rudderNormalCoefficient(r, r.stallAngle - eps), 1e-6);

    // Past stall the force falls away and never recovers past its stall value.
    const double atStall = rudderNormalCoefficient(r, r.stallAngle);
    expectTrue("C_N drops immediately past the stall angle",
               rudderNormalCoefficient(r, r.stallAngle + 10.0 * kDegToRad) < 0.9 * atStall);
    double peak = 0.0, peakAngle = 0.0;
    bool everBeatTheLinearLaw = false;
    for (int i = 0; i <= 1800; ++i) {
        const double a = i * 0.1 * kDegToRad;
        const double cn = rudderNormalCoefficient(r, a);
        if (cn > peak) { peak = cn; peakAngle = a; }
        if (a > r.stallAngle && cn >= fAlpha * std::sin(a)) everBeatTheLinearLaw = true;
    }
    expectTrue("C_N never exceeds the linear Fujii extrapolation past stall",
               !everBeatTheLinearLaw);
    expectNear("the peak of C_N is at the stall angle", peakAngle, r.stallAngle, 0.2 * kDegToRad);
    expectNear("the peak of C_N is the Fujii value there", peak, atStall, 1e-12);
}

// Closed form for the rudder with the shaft stopped and the ship on a straight
// course: the inflow is simply the wake-reduced surge speed times epsilon.
void testRudderForcesAgainstClosedForm() {
    const RudderParams r = kvlcc2Rudder();
    const PropellerParams p = kvlcc2Propeller();
    const PropellerState stopped = evaluatePropeller(p, 6.0, 0.0);

    const double delta = 20.0 * kDegToRad;
    const RudderState s = evaluateRudder(r, p, stopped, 6.0, 0.0, 0.0, delta);

    const double expectedInflow = r.raceFactor * 6.0 * (1.0 - p.wakeFraction);
    expectNear("inflow speed with the shaft stopped", s.inflowSpeed, expectedInflow, 1e-12);
    expectNear("angle of attack is the rudder angle on a straight course", s.angleOfAttack, delta,
               1e-15);

    const double fn = 0.5 * kRhoSeawater * r.area * expectedInflow * expectedInflow *
                      fujiiSlope(r) * std::sin(delta);
    expectNear("normal force from Fujii's formula", s.normalForce, fn, 1e-9 * fn);
    expectNear("surge force carries (1 - t_R)", s.surgeForce,
               -(1.0 - r.steeringResistance) * fn * std::sin(delta), 1e-9 * fn);
    expectNear("sway force carries (1 + a_H)", s.swayForce,
               -(1.0 + r.hullLiftFactor) * fn * std::cos(delta), 1e-9 * fn);
    expectNear("yaw moment acts at x_R + a_H x_H", s.yawMoment,
               (r.x + r.hullLiftFactor * r.hullLiftLever) * (-fn * std::cos(delta)),
               1e-9 * fn * std::abs(r.x));

    // Signs: a rudder put over to port pushes the stern to starboard and swings
    // the bow to port, and costs surge speed doing it.
    expectTrue("port rudder pushes the stern to starboard", s.swayForce < 0.0);
    expectTrue("port rudder yaws the bow to port", s.yawMoment > 0.0);
    expectTrue("rudder drag opposes surge", s.surgeForce < 0.0);
    const RudderState mirrored = evaluateRudder(r, p, stopped, 6.0, 0.0, 0.0, -delta);
    expectNear("starboard rudder is the exact mirror", mirrored.yawMoment, -s.yawMoment, 1e-12);
    expectNear("zero rudder makes no side force",
               evaluateRudder(r, p, stopped, 6.0, 0.0, 0.0, 0.0).swayForce, 0.0, 0.0);

    // Quadratic in inflow speed: doubling the speed quadruples the force.
    const RudderState fast = evaluateRudder(r, p, evaluatePropeller(p, 12.0, 0.0), 12.0, 0.0, 0.0,
                                            delta);
    expectNear("rudder force grows as the square of the inflow speed", fast.normalForce,
               4.0 * s.normalForce, 1e-9 * 4.0 * fn);

    // The propeller race must accelerate the flow over the rudder, which is why
    // a ship steers on the engine and not on speed alone.
    const PropellerState turning = evaluatePropeller(p, 6.0, 1.5);
    const RudderState inRace = evaluateRudder(r, p, turning, 6.0, 0.0, 0.0, delta);
    expectTrue("the propeller race speeds up the rudder inflow",
               inRace.inflowSpeed > 1.2 * s.inflowSpeed);
    expectTrue("and so increases the rudder force", inRace.normalForce > s.normalForce);

    // Bollard condition: dead in the water with the shaft turning, the rudder
    // still bites. The alternative J-based race formula divides by zero here.
    const PropellerState bollard = evaluatePropeller(p, 0.0, 1.5);
    const RudderState atRest = evaluateRudder(r, p, bollard, 0.0, 0.0, 0.0, delta);
    expectTrue("the rudder works at bollard pull", atRest.normalForce > 0.0);
    expectTrue("bollard rudder inflow is finite", std::isfinite(atRest.inflowSpeed));

    // Stall reaches the ship, not just the coefficient: a rudder put hard over
    // past the stall angle makes less turning moment than one at the stall.
    const RudderState atStall = evaluateRudder(r, p, stopped, 6.0, 0.0, 0.0, r.stallAngle);
    const RudderState overStalled =
        evaluateRudder(r, p, stopped, 6.0, 0.0, 0.0, r.stallAngle + 15.0 * kDegToRad);
    expectTrue("past the stall angle the yaw moment falls",
               overStalled.yawMoment < atStall.yawMoment);
    expectTrue("the stalled state is reported", overStalled.stalled && !atStall.stalled);
}

// --- MMG hull ---------------------------------------------------------------

void testHullSymmetryAndClosedForm() {
    const HullParams h = kvlcc2Hull();
    const double rho = kRhoSeawater;

    // Straight ahead, no drift, no yaw: nothing but resistance. A non-zero sway
    // force or yaw moment here means an asymmetry has crept into the
    // derivatives, and a ship that will not steer straight.
    for (double u : {1.0, 5.0, 9.0}) {
        const PlanarForces f = evaluateHull(h, u, 0.0, 0.0, rho);
        expectNear("no drift, no sway force @ u=" + std::to_string(u), f.sway, 0.0, 0.0);
        expectNear("no drift, no yaw moment @ u=" + std::to_string(u), f.yawMoment, 0.0, 0.0);
        expectNear("resistance is the closed form @ u=" + std::to_string(u), f.surge,
                   -0.5 * rho * h.length * h.draft * h.resistance * u * u,
                   1e-9 * std::abs(f.surge));
    }
    expectTrue("resistance opposes sternway too", evaluateHull(h, -5.0, 0.0, 0.0, rho).surge > 0.0);

    // Mirror invariance. The published MMG derivatives are identified in a frame
    // with y and r to starboard; they carry over to this frame's port-positive
    // convention only because every term is odd under the joint flip. Assert it
    // rather than trust it.
    for (double v : {-2.0, -0.5, 0.5, 2.0})
        for (double r : {-0.01, -0.002, 0.003, 0.02}) {
            const PlanarForces a = evaluateHull(h, 7.0, v, r, rho);
            const PlanarForces b = evaluateHull(h, 7.0, -v, -r, rho);
            const double scale = std::abs(a.sway) + std::abs(a.yawMoment) + 1.0;
            expectTrue("surge is even under the port/starboard mirror",
                       std::abs(a.surge - b.surge) < 1e-9 * scale);
            expectTrue("sway is odd under the port/starboard mirror",
                       std::abs(a.sway + b.sway) < 1e-9 * scale);
            expectTrue("yaw is odd under the port/starboard mirror",
                       std::abs(a.yawMoment + b.yawMoment) < 1e-6 * scale);
        }

    // Sway and yaw damping must oppose the motion, or the hull is unstable in
    // the most basic sense.
    const PlanarForces drifting = evaluateHull(h, 7.0, 1.0, 0.0, rho);
    expectTrue("drift to port makes a sway force to starboard", drifting.sway < 0.0);
    expectTrue("drift also makes drag", drifting.surge < evaluateHull(h, 7.0, 0.0, 0.0, rho).surge);
    const PlanarForces yawing = evaluateHull(h, 7.0, 0.0, 0.01, rho);
    expectTrue("yaw rate is damped", yawing.yawMoment < 0.0);

    // Nothing may blow up when the reference speed goes to zero -- the raw
    // non-dimensional polynomial has v'^4 and r'^3 terms that do.
    for (double u : {0.0, 1e-9, 1e-4})
        for (double r : {0.0, 0.05}) {
            const PlanarForces f = evaluateHull(h, u, 0.0, r, rho);
            expectTrue("hull forces stay finite at vanishing speed",
                       std::isfinite(f.surge) && std::isfinite(f.sway) &&
                           std::isfinite(f.yawMoment));
        }
}

// --- Turning circle ---------------------------------------------------------

struct TurnResult {
    double approachSpeed = 0;
    double steadySpeed = 0;
    double yawRate = 0;
    double driftAngle = 0;
    double radius = 0;         // from U / |r|
    double geometricRadius = 0;  // from three points on the track
    double yawRatio = 0;         // r' = r L / U
    double yawRateDrift = 0;     // change in r over the last minute of the run
};

// Run to a steady approach speed, then put the rudder over and run to a steady
// turn. The radius is measured two independent ways: from the kinematics
// (U / |r|) and from the circumradius of three widely spaced points on the
// actual track. They agree only if the turn really is steady and circular.
//
// The turn is held for 2400 s because it approaches its steady state
// exponentially with a time constant of about two minutes: at 600 s the yaw rate
// is still moving by 3e-3 per minute and the two radius measurements disagree in
// the fourth digit. That is the model converging, not the model wrong, and
// loosening the agreement tolerance instead of waiting would have hidden it.
TurnResult turningCircle(double rudderDeg, double revs) {
    Manoeuvring ship = kvlcc2();
    ship.revsPerSecond = revs;
    ship.state.surgeSpeed = 8.0;

    const double dt = 0.05;
    for (int i = 0; i < 12000; ++i) ship.step(dt);  // 600 s straight

    TurnResult out;
    out.approachSpeed = ship.state.surgeSpeed;
    expectNear("the approach run stays on a straight course (sway)", ship.state.swaySpeed, 0.0,
               0.0);
    expectNear("the approach run stays on a straight course (yaw)", ship.state.yawRate, 0.0, 0.0);

    ship.rudderAngle = rudderDeg * kDegToRad;
    for (int i = 0; i < 48000; ++i) ship.step(dt);  // 2400 s of turning

    const double rBefore = ship.state.yawRate;
    ship.state.x = 0;
    ship.state.y = 0;
    std::vector<double> px, py;
    for (int i = 0; i < 1200; ++i) {  // one more minute, sampling the track
        ship.step(dt);
        if (i % 400 == 0) { px.push_back(ship.state.x); py.push_back(ship.state.y); }
    }
    out.yawRateDrift = ship.state.yawRate - rBefore;

    const double u = ship.state.surgeSpeed, v = ship.state.swaySpeed;
    out.steadySpeed = std::sqrt(u * u + v * v);
    out.yawRate = ship.state.yawRate;
    out.driftAngle = std::atan2(v, u);
    out.radius = out.steadySpeed / std::abs(out.yawRate);
    out.yawRatio = out.yawRate * kvlcc2Hull().length / out.steadySpeed;

    // Circumradius of the three sampled track points: R = abc / (4 * area).
    const double ax = px[1] - px[0], ay = py[1] - py[0];
    const double bx = px[2] - px[1], by = py[2] - py[1];
    const double cx = px[2] - px[0], cy = py[2] - py[0];
    const double area = 0.5 * std::abs(ax * cy - ay * cx);
    out.geometricRadius = std::sqrt(ax * ax + ay * ay) * std::sqrt(bx * bx + by * by) *
                          std::sqrt(cx * cx + cy * cy) / (4.0 * area);
    return out;
}

void testTurningCircle() {
    const double length = kvlcc2Hull().length;
    const double revs = 1.86;  // ~ the KVLCC2 design shaft speed

    const TurnResult hard = turningCircle(35.0, revs);
    std::printf("     approach %.2f m/s (%.1f kn); 35 deg turn: U=%.2f m/s r=%.5f rad/s "
                "drift=%.1f deg R=%.0f m = %.2f L\n",
                hard.approachSpeed, hard.approachSpeed / 0.5144, hard.steadySpeed, hard.yawRate,
                hard.driftAngle * kRadToDeg, hard.radius, hard.radius / length);

    expectTrue("the approach speed is a plausible service speed",
               hard.approachSpeed > 7.0 && hard.approachSpeed < 9.0);

    // The turn must settle rather than diverge or wind in.
    expectTrue("the turn reaches a steady yaw rate",
               std::abs(hard.yawRateDrift) < 1e-6 * std::abs(hard.yawRate));
    expectTrue("the turn is to port for port rudder", hard.yawRate > 0.0);
    expectTrue("the ship drifts outward in the turn", hard.driftAngle < 0.0);
    // KVLCC2's measured drift angle at hard-over is a little under 20 degrees.
    // It is the characteristic the smaller sway derivatives actually control --
    // the radius barely moves when Y_r' is wrong, and the drift angle does.
    expectTrue("the drift angle at hard-over is a plausible 14 to 22 degrees",
               -hard.driftAngle > 14.0 * kDegToRad && -hard.driftAngle < 22.0 * kDegToRad);

    // A loaded tanker at hard-over turns in one to two ship lengths of radius;
    // KVLCC2's measured steady turning diameter is a shade over two lengths.
    expectTrue("the steady turning radius is a plausible few ship lengths",
               hard.radius > 0.8 * length && hard.radius < 3.0 * length);
    expectNear("kinematic and track-fitted radii agree", hard.geometricRadius, hard.radius,
               1e-5 * hard.radius);

    // The polynomial limiter must not be what produced that answer.
    expectTrue("the yaw-rate limiter is not active in a hard-over turn",
               std::abs(hard.yawRatio) < 0.9 * kMaxYawRatio);

    // More rudder, tighter turn, and more speed lost doing it.
    const TurnResult moderate = turningCircle(20.0, revs);
    const TurnResult gentle = turningCircle(10.0, revs);
    expectTrue("more rudder gives a tighter turn",
               hard.radius < moderate.radius && moderate.radius < gentle.radius);
    expectTrue("more rudder costs more speed",
               hard.steadySpeed < moderate.steadySpeed && moderate.steadySpeed < gentle.steadySpeed);
    expectTrue("a hard-over turn costs a large fraction of the approach speed",
               hard.steadySpeed < 0.7 * hard.approachSpeed);
    std::printf("     20 deg: R = %.2f L, U = %.2f m/s;  10 deg: R = %.2f L, U = %.2f m/s\n",
                moderate.radius / length, moderate.steadySpeed, gentle.radius / length,
                gentle.steadySpeed);

    // A starboard turn is the exact mirror of the port turn. Any drift here is
    // an asymmetry in the assembled model, not in the hull alone.
    const TurnResult starboard = turningCircle(-35.0, revs);
    expectNear("a starboard turn mirrors the port turn: yaw rate", starboard.yawRate,
               -hard.yawRate, 1e-10 * std::abs(hard.yawRate));
    expectNear("a starboard turn mirrors the port turn: radius", starboard.radius, hard.radius,
               1e-10 * hard.radius);
    expectNear("a starboard turn mirrors the port turn: speed", starboard.steadySpeed,
               hard.steadySpeed, 1e-10 * hard.steadySpeed);
}

// The four-quadrant propeller earns its keep here: from a steady approach, order
// full astern and the ship must slow down, with the propeller pulling backwards
// the whole way. The head reach is reported but not asserted against a real
// ship's -- with no engine model yet the shaft reverses instantly and delivers
// full torque astern, so the distance is optimistic by a large factor. See
// docs/02-simulation.md section 7.
void testCrashStop() {
    Manoeuvring ship = kvlcc2();
    ship.revsPerSecond = 1.86;
    ship.state.surgeSpeed = 8.0;
    const double dt = 0.05;
    for (int i = 0; i < 12000; ++i) ship.step(dt);
    const double approach = ship.state.surgeSpeed;

    ship.revsPerSecond = -1.86;
    ship.state.x = 0;
    bool everAccelerated = false, everPushedAhead = false, stopped = false;
    double previous = approach, headReach = 0.0, stopTime = 0.0;
    for (int i = 0; i < 24000; ++i) {  // 1200 s
        ship.step(dt);
        if (previous > 0.0 && ship.state.surgeSpeed > previous) everAccelerated = true;
        if (ship.state.surgeSpeed > 0.0 && ship.propellerState().thrust > 0.0)
            everPushedAhead = true;
        if (!stopped && ship.state.surgeSpeed <= 0.0) {
            stopped = true;
            headReach = ship.state.x;
            stopTime = (i + 1) * dt;
        }
        previous = ship.state.surgeSpeed;
    }
    expectTrue("a crash stop never speeds the ship up while it still has headway",
               !everAccelerated);
    expectTrue("the propeller never pushes ahead during a crash stop", !everPushedAhead);
    expectTrue("the ship stops within the run", stopped);
    expectTrue("the head reach is ahead of where the order was given", headReach > 0.0);
    expectTrue("the ship ends up making sternway", ship.state.surgeSpeed < 0.0);
    expectNear("the ship stays on a straight course through the crash stop", ship.state.yawRate,
               0.0, 0.0);
    expectNear("no rudder means no sway through the crash stop", ship.state.swaySpeed, 0.0, 0.0);
    std::printf("     crash stop from %.2f m/s: head reach %.0f m (%.2f L) in %.0f s\n", approach,
                headReach, headReach / kvlcc2Hull().length, stopTime);
}

}  // namespace

void runPropulsionTests() {
    std::printf("\n--- propulsion, steering and manoeuvring ---\n");
    testDimensionalIdentities();
    testBollardPullAndMonotonicity();
    testOpenWaterEfficiency();
    testPowerConsistency();
    testFourQuadrantSigns();
    testRudderNormalCoefficient();
    testRudderForcesAgainstClosedForm();
    testHullSymmetryAndClosedForm();
    testTurningCircle();
    testCrashStop();
}
