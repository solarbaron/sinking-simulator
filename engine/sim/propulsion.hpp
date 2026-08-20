// SPDX-License-Identifier: MIT
//
// Propulsion, steering and manoeuvring: propeller, rudder and MMG hull, kept as
// three separately identified terms rather than one derivative set.
//
// The modularity is the point. A monolithic set of manoeuvring derivatives is
// identified for one intact ship at one loading; the moment a propeller is
// damaged, a rudder jams or a shaft stops, the whole set is invalid. Here the
// hull, the propeller and the rudder each produce their own surge/sway/yaw
// contribution, so losing one degrades one term and the rest keep meaning what
// they meant.
//
// Frame and sign conventions (matching engine/core/math.hpp):
//
//   * body frame is +x forward (bow), +y to port, +z up;
//   * the yaw rate r is about +z, so positive r swings the bow to port;
//   * the sway speed v is positive to port;
//   * the rudder angle delta is positive when the trailing edge swings to port,
//     which is the deflection that turns the bow to port -- i.e. positive delta
//     produces positive yaw rate, consistent with the right-hand rule about +z;
//   * propeller revolutions n are positive for ahead rotation. Negative n is
//     astern rotation and is handled by the same code path, not a special case.
//
// The published MMG literature uses the opposite lateral convention (y and r to
// starboard, N positive to starboard). Every term of the hull polynomial below
// is odd under the simultaneous flip (v, r, Y, N) -> (-v, -r, -Y, -N), because a
// port/starboard-symmetric hull cannot be anything else, so the published
// coefficients transfer into this frame unchanged. That invariance is asserted
// in tests/test_propulsion.cpp rather than assumed.
#pragma once

#include "../core/math.hpp"

namespace sim {

// ---------------------------------------------------------------------------
// Propeller
// ---------------------------------------------------------------------------
//
// Open-water characteristics come from a single equivalent blade section at
// 0.7R, evaluated over all four quadrants of (advance speed, rotation rate).
//
// Why not the Wageningen B-series regression: that regression is 39 terms for
// K_T and 47 for 10 K_Q, it is first-quadrant only, and a single mis-transcribed
// coefficient produces a curve that looks entirely plausible and is wrong. It
// also says nothing about astern rotation, which is precisely the case a crash
// stop needs. The blade-element form below instead makes the invariants that
// matter structural rather than incidental:
//
//   * thrust reverses exactly when rotation reverses at zero advance, because
//     the section angle of attack does;
//   * open-water efficiency is provably <= 1 (see below), so a sign or factor
//     error cannot hide behind a plausible-looking curve;
//   * astern rotation with ahead flow gives astern thrust by construction, so a
//     crash stop cannot produce forward thrust.
//
// With the resultant inflow angle at 0.7R
//     beta = atan2(Va, 0.7 pi n D)
// and the section lift and drag resolved into axial and tangential components,
//     eta = J K_T / (2 pi K_Q) = tan(beta) * (L cos b - D sin b) / (L sin b + D cos b)
// which is identically 1 for zero section drag and strictly less than 1 for any
// positive drag. The Betz-style ceiling is therefore a property of the algebra,
// not of the coefficients.
//
// The free constants (solidity, section lift slope, saturation and drag, pitch
// effectiveness) are calibrated so that the first-quadrant curve of a
// P/D = 1.0, A_E/A_0 = 0.70 four-bladed propeller reproduces the published
// Wageningen B4-70 magnitudes: K_T(0) ~ 0.35, 10 K_Q(0) ~ 0.49, zero thrust at
// J ~ 0.85, peak open-water efficiency ~ 0.67 near J ~ 0.71. They are a fit, not
// measurements; see docs/02-simulation.md section 7.
struct PropellerParams {
    double diameter = 5.0;          // D, m
    double pitchRatio = 1.0;        // P/D at 0.7R
    double bladeAreaRatio = 0.70;   // A_E / A_0
    double wakeFraction = 0.25;     // w, so Va = u (1 - w)
    double thrustDeduction = 0.20;  // t, so the hull feels (1 - t) T

    // --- Blade-section model constants (calibration, not measurement) -------
    // Hydrodynamic pitch is less than geometric pitch because the induced
    // velocity turns the inflow; this factor carries that reduction and sets
    // where K_T crosses zero.
    double pitchEffectiveness = 0.853;
    // Lift-curve slope of the equivalent section, reduced from 2 pi by finite
    // span and cascade interference.
    double sectionLiftSlope = 4.0;
    // Saturation of the section normal force: propeller blades in cascade stall
    // gently, so the normal force rolls off smoothly instead of breaking down.
    double sectionNormalForceMax = 1.4;
    double sectionDragCoeff = 0.018;  // profile drag at zero incidence
    // Effective solidity, sigma = solidity * (A_E/A_0)^solidityExponent. The
    // exponent is well below 1 because added blade area also raises the induced
    // velocity, which is why real K_T depends only weakly on A_E/A_0.
    double solidity = 0.194;
    double solidityExponent = 0.2;
    double representativeRadius = 0.7;  // r/R of the equivalent section
};

struct PropellerState {
    double revsPerSecond = 0;      // n, as supplied; negative is astern rotation
    double advanceSpeed = 0;       // Va = u (1 - w), m/s
    // The three shaft-normalised quantities are reported only while the shaft,
    // rather than the inflow, sets the water speed at the blade -- specifically
    // while (0.7 pi n D)^2 > 0.1 (Va^2 + (0.7 pi n D)^2), which is |J| < 6.6 and
    // caps |K_T| at 4.7. Outside that they are zero, which covers a locked shaft
    // and also a shaft crawling through zero during a reversal: there they are
    // not merely undefined, they diverge as 1/n^2, and a reversal walks straight
    // through the divergence. See evaluatePropeller for the derivation. Thrust,
    // torque and power below are unaffected and stay valid everywhere.
    double advanceRatio = 0;       // J = Va / (n D)
    double thrustCoefficient = 0;  // K_T = T / (rho n^2 D^4)
    double torqueCoefficient = 0;  // K_Q = Q / (rho n^2 D^5)
    double thrust = 0;             // T, N, along +x
    double thrustOnHull = 0;       // (1 - t) T, N -- what the hull actually feels
    double torque = 0;             // Q, N m; positive means the shaft is driven
    double deliveredPower = 0;     // 2 pi n Q, W
    double thrustPower = 0;        // T Va, W
    double efficiency = 0;         // open-water eta; zero outside the working range
};

// Advance speed at the propeller disc.
double advanceSpeed(const PropellerParams& propeller, double surgeSpeed);

// Open-water coefficients as functions of J. Defined for ahead rotation (n > 0)
// and any sign of J; for astern rotation use evaluatePropeller, which is not
// parameterisable by J alone.
double propellerThrustCoefficient(const PropellerParams& propeller, double advanceRatio);
double propellerTorqueCoefficient(const PropellerParams& propeller, double advanceRatio);

// eta = J K_T / (2 pi K_Q), and zero outside the working range (where either
// coefficient is non-positive and the ratio stops meaning "useful work out over
// work in").
double propellerOpenWaterEfficiency(const PropellerParams& propeller, double advanceRatio);

// Full four-quadrant evaluation. surgeSpeed is the ship's speed through the
// water along +x; revsPerSecond may be either sign.
PropellerState evaluatePropeller(const PropellerParams& propeller, double surgeSpeed,
                                 double revsPerSecond, double density = kRhoSeawater);

// ---------------------------------------------------------------------------
// Rudder
// ---------------------------------------------------------------------------
//
// Normal force from Fujii's formula,
//     F_N = 0.5 rho A_R f_alpha U_R^2 sin(alpha_R),   f_alpha = 6.13 L / (L + 2.25)
// (L here is the rudder's geometric aspect ratio), with the inflow speed and
// angle taken inside the propeller race, and a stall model past the critical
// angle so the linear law is never extrapolated to the stops.
struct RudderParams {
    double area = 20.0;         // A_R, m^2
    double aspectRatio = 1.827; // geometric aspect ratio, drives Fujii's f_alpha
    double span = 6.0;          // H_R, m; the race covers D/H_R of it
    double x = -50.0;           // longitudinal position, body frame (negative = aft)

    // Stall. Below stallAngle the Fujii law holds; past it the normal force
    // falls away from the linear extrapolation towards the finite-aspect-ratio
    // flat-plate value, with a dip in between. A ship's rudder does not simply
    // stop working past 35 degrees, but it stops gaining.
    double stallAngle = 30.0 * kDegToRad;
    double postStallWidth = 20.0 * kDegToRad;
    double postStallDrop = 0.25;
    double broadsideCoeff = 1.20;  // normal-force coefficient at 90 deg incidence

    // --- Hull interaction --------------------------------------------------
    double flowStraightening = 0.5;   // gamma_R
    double flowStraighteningLever = -227.2;  // l_R, m (negative = aft of midship)
    double raceFactor = 1.09;         // epsilon: wake ratio hull/propeller
    double raceContraction = 0.5;     // kappa: fraction of the race the rudder sees
    double steeringResistance = 0.387;  // t_R, so surge gets (1 - t_R) of the force
    double hullLiftFactor = 0.312;      // a_H, extra hull lift induced by the rudder
    double hullLiftLever = -148.5;      // x_H, m, where that extra lift acts
};

struct RudderState {
    double inflowSpeed = 0;        // U_R, m/s
    double axialInflow = 0;        // u_R, m/s (inside the race)
    double lateralInflow = 0;      // v_R, m/s (after flow straightening)
    double angleOfAttack = 0;      // alpha_R, rad, wrapped to (-pi, pi]
    double normalCoefficient = 0;  // C_N(alpha_R)
    double normalForce = 0;        // F_N, N; positive pushes the stern to starboard
    double surgeForce = 0;         // X_R, N, hull-interaction factor applied
    double swayForce = 0;          // Y_R, N, ditto
    double yawMoment = 0;          // N_R, N m about +z, ditto
    bool   stalled = false;
};

// C_N(alpha): Fujii's f_alpha sin(alpha) below the stall angle, rolling off
// past it. Odd in alpha, and zero at 0 and +-pi.
double rudderNormalCoefficient(const RudderParams& rudder, double angleOfAttack);

// The rudder in the propeller race. u, v, r are the ship's surge, sway and yaw
// rate in the body frame; rudderAngle is delta.
RudderState evaluateRudder(const RudderParams& rudder, const PropellerParams& propeller,
                           const PropellerState& propellerState, double u, double v, double r,
                           double rudderAngle, double density = kRhoSeawater);

// ---------------------------------------------------------------------------
// MMG hull
// ---------------------------------------------------------------------------
//
// Yasukawa & Yoshimura's MMG standard-method polynomial, in the cancelled form
// where the reference speed U divides out of each term. Writing it that way is
// not cosmetic: the raw non-dimensional form has v'^4 and r'^3 terms that are
// singular at U = 0, and this form is bounded there.
struct HullParams {
    double length = 320.0;  // L between perpendiculars, m
    double beam = 58.0;     // B, m
    double draft = 20.8;    // d, m
    double blockCoefficient = 0.810;
    double xCog = 11.2;     // x_G, m, positive forward of midship

    // Added mass and added yaw inertia, non-dimensionalised on 0.5 rho L^2 d
    // and 0.5 rho L^4 d respectively.
    double addedMassSurge = 0.022;   // m_x'
    double addedMassSway  = 0.223;   // m_y'
    double addedInertiaYaw = 0.011;  // J_z'
    double yawGyradiusRatio = 0.25;  // k_zz / L, about the centre of gravity

    double resistance = 0.022;  // R_0'

    double surgeVV = -0.040;    // X_vv'
    double surgeVR = 0.002;     // X_vr'
    double surgeRR = 0.011;     // X_rr'
    double surgeVVVV = 0.771;   // X_vvvv'

    double swayV = -0.315;      // Y_v'
    double swayR = 0.083;       // Y_r'
    double swayVVV = -1.607;    // Y_vvv'
    double swayVVR = 0.379;     // Y_vvr'
    double swayVRR = -0.391;    // Y_vrr'
    double swayRRR = 0.008;     // Y_rrr'

    double yawV = -0.137;       // N_v'
    double yawR = -0.049;       // N_r'
    double yawVVV = -0.030;     // N_vvv'
    double yawVVR = -0.294;     // N_vvr'
    double yawVRR = 0.055;      // N_vrr'
    double yawRRR = -0.013;     // N_rrr'

    double mass(double density = kRhoSeawater) const {
        return density * length * beam * draft * blockCoefficient;
    }
};

// Non-dimensional yaw rate r' = r L / U is clamped to this before the
// polynomial is evaluated. A hard-over turn reaches about 0.9; beyond ~1.5 the
// polynomial fit is extrapolating far outside its identification range and the
// clamp keeps a near-stationary, spinning ship from producing nonsense rather
// than pretending the fit is still valid there.
inline constexpr double kMaxYawRatio = 1.5;

struct PlanarForces {
    double surge = 0;      // X, N
    double sway = 0;       // Y, N
    double yawMoment = 0;  // N, N m about +z, at midship
};

PlanarForces evaluateHull(const HullParams& hull, double u, double v, double r,
                          double density = kRhoSeawater);

// ---------------------------------------------------------------------------
// Assembled 3-DOF manoeuvring model
// ---------------------------------------------------------------------------
//
// Surge, sway and yaw only. This is deliberately not the ship's rigid body --
// engine/sim/ship.cpp owns that -- but a self-contained horizontal-plane model
// so the manoeuvring terms can be exercised, and a turning circle measured,
// without dragging in flooding and hydrostatics.
struct ManoeuvringState {
    double surgeSpeed = 0;  // u, m/s
    double swaySpeed = 0;   // v, m/s (positive to port)
    double yawRate = 0;     // r, rad/s (positive swings the bow to port)
    double x = 0;           // world position, m
    double y = 0;
    double heading = 0;     // rad, from +x world, counter-clockwise
};

struct Manoeuvring {
    HullParams hull;
    PropellerParams propeller;
    RudderParams rudder;
    double density = kRhoSeawater;

    double revsPerSecond = 0;  // commanded shaft speed
    double rudderAngle = 0;    // delta, rad

    ManoeuvringState state;

    // Body-frame totals, hull + propeller + rudder.
    PlanarForces forces() const;

    // Individual contributions, for instrumentation and for damage models that
    // need to scale one term.
    PropellerState propellerState() const;
    RudderState rudderState() const;

    void step(double dt);
};

// Reference ship: the KVLCC2 tanker, the SIMMAN benchmark hull the MMG standard
// method was published against. See docs/02-simulation.md section 7 for which of
// these numbers are published and which are placeholders.
HullParams kvlcc2Hull();
PropellerParams kvlcc2Propeller();
RudderParams kvlcc2Rudder();
Manoeuvring kvlcc2();

}  // namespace sim
