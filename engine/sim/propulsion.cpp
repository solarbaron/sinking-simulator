// SPDX-License-Identifier: MIT
#include "propulsion.hpp"

#include <algorithm>
#include <cmath>

namespace sim {
namespace {

// Wrap an angle into (-pi, pi]. Everything downstream is a sine or cosine of it,
// so this matters only for reporting and for the |alpha| used by the stall
// model -- but a stall model fed 350 degrees instead of -10 is exactly the kind
// of quiet wrongness this project is trying to avoid.
double wrapAngle(double a) { return std::atan2(std::sin(a), std::cos(a)); }

// Axial and tangential force of the equivalent blade section at 0.7R, already
// divided by the disc dynamic pressure 0.5 rho V_R^2 (pi/4) D^2, and with the
// tangential one already carrying its moment arm r_e / D = 0.5 * (r/R).
//
// beta is the resultant inflow angle at that radius, atan2(Va, 0.7 pi n D). It
// spans the full circle, so the four quadrants are one expression rather than
// four cases.
struct SectionForces {
    double axial = 0;       // -> C_T*
    double tangential = 0;  // -> C_Q*
};

SectionForces bladeSection(const PropellerParams& p, double beta) {
    const double theta =
        std::atan(p.pitchEffectiveness * p.pitchRatio / (p.representativeRadius * kPi));
    const double alpha = wrapAngle(theta - beta);

    // Normal force saturating towards a flat-plate value: linear at small
    // incidence, rolling over as the section stalls. tanh keeps it monotone in
    // sin(alpha), which is what makes K_T monotone in J.
    const double sa = std::sin(alpha);
    const double cn = p.sectionNormalForceMax *
                      std::tanh(p.sectionLiftSlope * sa / p.sectionNormalForceMax);
    const double lift = cn * std::cos(alpha);
    const double drag = cn * sa + p.sectionDragCoeff;

    const double sigma = p.solidity * std::pow(p.bladeAreaRatio, p.solidityExponent);
    const double sb = std::sin(beta), cb = std::cos(beta);
    SectionForces f;
    f.axial = sigma * (lift * cb - drag * sb);
    f.tangential = 0.5 * p.representativeRadius * sigma * (lift * sb + drag * cb);
    return f;
}

// K_T and K_Q share the same factor: with V_R^2 = n^2 D^2 (J^2 + (0.7 pi)^2),
//     T = C_T* * 0.5 rho V_R^2 (pi/4) D^2 = C_T* (pi/8)(J^2 + (0.7 pi)^2) rho n^2 D^4
// so the bracket below is exactly the conversion from disc coefficient to
// K-coefficient. Any factor error here shows up as a ship with the wrong thrust
// by a constant, which is why the tests check T = rho n^2 D^4 K_T directly.
double coefficientScale(const PropellerParams& p, double advanceRatio) {
    const double vt = p.representativeRadius * kPi;
    return (kPi / 8.0) * (advanceRatio * advanceRatio + vt * vt);
}

}  // namespace

double advanceSpeed(const PropellerParams& propeller, double surgeSpeed) {
    return surgeSpeed * (1.0 - propeller.wakeFraction);
}

double propellerThrustCoefficient(const PropellerParams& p, double advanceRatio) {
    const double beta = std::atan2(advanceRatio, p.representativeRadius * kPi);
    return coefficientScale(p, advanceRatio) * bladeSection(p, beta).axial;
}

double propellerTorqueCoefficient(const PropellerParams& p, double advanceRatio) {
    const double beta = std::atan2(advanceRatio, p.representativeRadius * kPi);
    return coefficientScale(p, advanceRatio) * bladeSection(p, beta).tangential;
}

double propellerOpenWaterEfficiency(const PropellerParams& p, double advanceRatio) {
    const double kt = propellerThrustCoefficient(p, advanceRatio);
    const double kq = propellerTorqueCoefficient(p, advanceRatio);
    if (advanceRatio <= 0.0 || kt <= 0.0 || kq <= 0.0) return 0.0;
    return advanceRatio * kt / (2.0 * kPi * kq);
}

PropellerState evaluatePropeller(const PropellerParams& p, double surgeSpeed,
                                 double revsPerSecond, double density) {
    PropellerState s;
    s.revsPerSecond = revsPerSecond;
    s.advanceSpeed = advanceSpeed(p, surgeSpeed);

    const double tangentialSpeed = p.representativeRadius * kPi * revsPerSecond * p.diameter;
    const double resultantSq = s.advanceSpeed * s.advanceSpeed + tangentialSpeed * tangentialSpeed;
    if (resultantSq < 1e-18) return s;  // dead in the water with the shaft stopped

    const double beta = std::atan2(s.advanceSpeed, tangentialSpeed);
    const SectionForces f = bladeSection(p, beta);

    const double discPressure = 0.5 * density * resultantSq * (kPi / 4.0) * p.diameter * p.diameter;
    s.thrust = f.axial * discPressure;
    s.torque = f.tangential * discPressure * p.diameter;
    s.thrustOnHull = (1.0 - p.thrustDeduction) * s.thrust;
    s.deliveredPower = 2.0 * kPi * revsPerSecond * s.torque;
    s.thrustPower = s.thrust * s.advanceSpeed;

    // J, K_T and K_Q are the *shaft-normalised* coefficients -- J = Va/(nD),
    // K_T = T/(rho n^2 D^4), K_Q = Q/(rho n^2 D^5) -- so they mean something only
    // while the shaft, and not the inflow, sets the resultant the blade sees.
    // That is the condition this guard has to express, and it is a statement
    // about a *ratio*, never about n on its own.
    //
    // It used to read `std::abs(revsPerSecond) > 1e-9`, which is not scale-free
    // for the same reason ||r|| is not scale-free in solid_shell.cpp's EAS loop:
    // 1e-9 is a bare number laid against a frequency in Hz, with no reference
    // frequency anywhere in the expression to divide it by. Whether a shaft is
    // "turning" is not a property of n. Measured on the KVLCC2 propeller
    // (D = 9.86 m) at 6 m/s of headway, with the shaft passing through zero --
    // which every reversal does, continuously, not as a contrived case:
    //
    //     n = 1e-4 Hz   old guard passes   J = 4.0e3   K_T = -1.5e6
    //     n = 1e-6 Hz   old guard passes   J = 4.0e5   K_T = -1.5e10
    //     n = 1e-8 Hz   old guard passes   J = 4.0e7   K_T = -1.5e14
    //
    // against quantities propulsion.hpp documents as J ~ 0.85 and K_T ~ 0.1-0.4.
    // It caught n at or within a nanohertz of exactly zero and nothing else, so
    // the infinities it exists to keep off an instrument panel walked past it.
    //
    // The scale-free form. tangentialSpeed = 0.7 pi n D and resultantSq =
    // Va^2 + tangentialSpeed^2 are computed six lines up in m/s and m^2/s^2, so
    // their ratio is dimensionless, lies in [0, 1] for any ship, propeller or
    // unit system, and is exactly "how much of the water the blade meets is put
    // there by the blade": 1 at bollard pull, 0 for a locked shaft. It is also J
    // in disguise, since J = Va/(nD) = 0.7 pi Va / tangentialSpeed, so with
    // Vt^2 > eps (Va^2 + Vt^2) the admitted band is |J| < 0.7 pi sqrt((1-eps)/eps).
    //
    // Which fixes eps by arithmetic rather than by taste. coefficientScale above
    // is (pi/8)(J^2 + (0.7 pi)^2), and at the threshold J^2 + (0.7 pi)^2 reduces
    // to (0.7 pi)^2 / eps, so the largest coefficient this block can ever report
    // is exactly
    //
    //     |K_T| <= max|C_T*| (pi/8)(0.7 pi)^2 / eps = 0.245 * 1.899 / eps
    //
    // max|C_T*| = 0.245 being the disc coefficient's own ceiling, measured by
    // sweeping the section model over every beta (it is bounded by solidity 0.181
    // times a normal force saturating at 1.4). eps = 0.1 gives |J| < 6.6 and
    // |K_T| <= 4.7, about ten times the documented ceiling of 0.4; eps = 1e-6,
    // the loosest form that is still dimensionally coherent, would give |J| < 2199
    // and |K_T| <= 4.7e5, which is no better than what it replaces. No threshold
    // removes the divergence -- it only chooses where to cut it -- and this one
    // cuts while the number is still the right order of magnitude. Past |J| ~ 7
    // the propeller is a drag disc and not a propulsor; thrust, torque and power
    // are computed outside this block and stay finite and correct there.
    //
    // Squared rather than `std::abs(n) * D > 1e-3 * std::sqrt(resultantSq)`: same
    // criterion, no sqrt, and nothing to say about resultantSq == 0, which the
    // early return above has already taken.
    constexpr double kMinShaftShareOfResultant = 0.1;
    if (tangentialSpeed * tangentialSpeed > kMinShaftShareOfResultant * resultantSq) {
        const double nd = revsPerSecond * p.diameter;
        s.advanceRatio = s.advanceSpeed / nd;
        const double d4 = p.diameter * p.diameter * p.diameter * p.diameter;
        s.thrustCoefficient = s.thrust / (density * revsPerSecond * revsPerSecond * d4);
        s.torqueCoefficient = s.torque / (density * revsPerSecond * revsPerSecond * d4 * p.diameter);
        if (s.thrustPower > 0.0 && s.deliveredPower > 0.0)
            s.efficiency = s.thrustPower / s.deliveredPower;
    }
    return s;
}

// ---------------------------------------------------------------------------
// Rudder
// ---------------------------------------------------------------------------

namespace {

// The coefficient multiplying sin(alpha) in the rudder normal force: Fujii's
// f_alpha below the stall angle, and past it a value decaying from f_alpha
// towards the broadside flat-plate coefficient with a dip superimposed on the
// way down. Keeping the sin(alpha) factor outside the piecewise definition is
// what makes the whole curve continuous at the stall angle and exactly zero at
// zero and +-pi -- a piecewise fit gets both wrong by a little and never says so.
double rudderNormalSlope(const RudderParams& rudder, double absAlpha) {
    const double fAlpha = 6.13 * rudder.aspectRatio / (rudder.aspectRatio + 2.25);
    if (absAlpha <= rudder.stallAngle) return fAlpha;
    const double x = (absAlpha - rudder.stallAngle) / std::max(rudder.postStallWidth, 1e-9);
    return rudder.broadsideCoeff + (fAlpha - rudder.broadsideCoeff) * std::exp(-x) -
           rudder.postStallDrop * fAlpha * x * std::exp(1.0 - x);
}

}  // namespace

double rudderNormalCoefficient(const RudderParams& rudder, double angleOfAttack) {
    const double alpha = wrapAngle(angleOfAttack);
    return rudderNormalSlope(rudder, std::abs(alpha)) * std::sin(alpha);
}

RudderState evaluateRudder(const RudderParams& rudder, const PropellerParams& propeller,
                           const PropellerState& propellerState, double u, double v, double r,
                           double rudderAngle, double density) {
    RudderState s;

    // Axial inflow. Inside the race the accelerated slipstream is
    //   u_race = u_P (1 - kappa) + kappa sqrt(u_P^2 + 8 K_T n^2 D^2 / pi),
    // which is the usual u_P (1 + kappa (sqrt(1 + 8 K_T / (pi J^2)) - 1)) with
    // the u_P divided into the root. The two are identical for J > 0 and only
    // the second is finite at bollard pull, where J -> 0 with 8 K_T / (pi J^2)
    // -> infinity. Writing it the first way is a division by zero at exactly the
    // condition a tug spends its life in.
    const double uP = u * (1.0 - propeller.wakeFraction);
    const double n = propellerState.revsPerSecond;
    double uRace = uP;
    if (uP >= 0.0 && n > 0.0 && propellerState.thrustCoefficient > 0.0) {
        const double nd = n * propeller.diameter;
        const double inner =
            uP * uP + (8.0 / kPi) * propellerState.thrustCoefficient * nd * nd;
        uRace = uP * (1.0 - rudder.raceContraction) +
                rudder.raceContraction * std::sqrt(std::max(inner, 0.0));
    }
    // Only the D/H_R of the rudder span that lies in the race sees the
    // accelerated flow; the rest sees the hull wake.
    const double covered = std::clamp(propeller.diameter / std::max(rudder.span, 1e-9), 0.0, 1.0);
    const double blended = std::sqrt(covered * uRace * uRace + (1.0 - covered) * uP * uP);
    s.axialInflow = rudder.raceFactor * (uRace < 0.0 ? -blended : blended);

    // Lateral inflow: the sway velocity at the rudder's effective lever, reduced
    // by flow straightening along the hull.
    s.lateralInflow = rudder.flowStraightening * (v + rudder.flowStraighteningLever * r);

    s.inflowSpeed = std::sqrt(s.axialInflow * s.axialInflow + s.lateralInflow * s.lateralInflow);
    if (s.inflowSpeed < 1e-9) return s;

    // Angle between the rudder chord and the flow. With the chord pointing
    // forward as (cos d, -sin d), the water arriving at (-u_R, -v_R) meets it at
    // delta + atan2(v_R, u_R) -- which for reversed flow (u_R < 0, going astern)
    // lands near +-pi and correctly reverses the force.
    s.angleOfAttack = wrapAngle(rudderAngle + std::atan2(s.lateralInflow, s.axialInflow));
    s.stalled = std::abs(s.angleOfAttack) > rudder.stallAngle;

    // sin(alpha_R) written out as sin(delta) cos(phi) + cos(delta) sin(phi) with
    // cos(phi) = u_R / U_R and sin(phi) = v_R / U_R, rather than round-tripped
    // through atan2 and back through sin. Algebraically identical; numerically
    // it is exactly zero where it should be. Astern with no drift and no rudder,
    // the round trip lands on sin(pi) = 1.2e-16 instead of 0, which is a phantom
    // yaw moment on a ship that is supposed to be running dead straight.
    const double sinAlpha = (std::sin(rudderAngle) * s.axialInflow +
                             std::cos(rudderAngle) * s.lateralInflow) / s.inflowSpeed;
    s.normalCoefficient = rudderNormalSlope(rudder, std::abs(s.angleOfAttack)) * sinAlpha;

    s.normalForce = 0.5 * density * rudder.area * s.inflowSpeed * s.inflowSpeed *
                    s.normalCoefficient;

    // The force acts along the chord normal, which for positive delta (trailing
    // edge to port) points to port; the reaction on the ship is to starboard.
    const double fx = -s.normalForce * std::sin(rudderAngle);
    const double fy = -s.normalForce * std::cos(rudderAngle);

    // Hull interaction: the rudder's own drag is partly cancelled by the change
    // in hull resistance (t_R), the hull carries extra lift induced by the
    // rudder (a_H) and that extra lift acts at x_H, not at the rudder.
    s.surgeForce = (1.0 - rudder.steeringResistance) * fx;
    s.swayForce = (1.0 + rudder.hullLiftFactor) * fy;
    s.yawMoment = (rudder.x + rudder.hullLiftFactor * rudder.hullLiftLever) * fy;
    return s;
}

// ---------------------------------------------------------------------------
// MMG hull
// ---------------------------------------------------------------------------

PlanarForces evaluateHull(const HullParams& h, double u, double v, double r, double density) {
    PlanarForces f;
    const double L = h.length;
    const double speed = std::sqrt(u * u + v * v);
    const double q = 0.5 * density * L * h.draft;

    // Resistance always opposes surge, so it carries u|u| rather than U^2.
    f.surge = -q * h.resistance * u * std::abs(u);

    if (speed < 1e-9) return f;  // no drift, no yaw: only resistance survives

    // Written with U divided out term by term. v' = v/U and r' = rL/U, so e.g.
    // U^2 X_vv' v'^2 = X_vv' v^2 exactly, with no cancellation error and no
    // singularity as U -> 0.
    const double rl = std::clamp(r * L / speed, -kMaxYawRatio, kMaxYawRatio) * speed;
    const double vv = v * v, rr = rl * rl;

    f.surge += q * (h.surgeVV * vv + h.surgeVR * v * rl + h.surgeRR * rr +
                    h.surgeVVVV * vv * vv / (speed * speed));
    f.sway = q * (h.swayV * v * speed + h.swayR * rl * speed + h.swayVVV * vv * v / speed +
                  h.swayVVR * vv * rl / speed + h.swayVRR * v * rr / speed +
                  h.swayRRR * rr * rl / speed);
    f.yawMoment = q * L *
                  (h.yawV * v * speed + h.yawR * rl * speed + h.yawVVV * vv * v / speed +
                   h.yawVVR * vv * rl / speed + h.yawVRR * v * rr / speed +
                   h.yawRRR * rr * rl / speed);
    return f;
}

// ---------------------------------------------------------------------------
// Assembled model
// ---------------------------------------------------------------------------

PropellerState Manoeuvring::propellerState() const {
    return evaluatePropeller(propeller, state.surgeSpeed, revsPerSecond, density);
}

RudderState Manoeuvring::rudderState() const {
    return evaluateRudder(rudder, propeller, propellerState(), state.surgeSpeed, state.swaySpeed,
                          state.yawRate, rudderAngle, density);
}

PlanarForces Manoeuvring::forces() const {
    const PropellerState prop =
        evaluatePropeller(propeller, state.surgeSpeed, revsPerSecond, density);
    const RudderState rud = evaluateRudder(rudder, propeller, prop, state.surgeSpeed,
                                           state.swaySpeed, state.yawRate, rudderAngle, density);
    PlanarForces f = evaluateHull(hull, state.surgeSpeed, state.swaySpeed, state.yawRate, density);
    f.surge += prop.thrustOnHull + rud.surgeForce;
    f.sway += rud.swayForce;
    f.yawMoment += rud.yawMoment;
    return f;
}

void Manoeuvring::step(double dt) {
    const PlanarForces f = forces();

    const double L = hull.length;
    const double m = hull.mass(density);
    const double refMass = 0.5 * density * L * L * hull.draft;
    const double mx = hull.addedMassSurge * refMass;
    const double my = hull.addedMassSway * refMass;
    const double jz = hull.addedInertiaYaw * (0.5 * density * L * L * L * L * hull.draft);
    const double izz = m * (hull.yawGyradiusRatio * L) * (hull.yawGyradiusRatio * L);
    const double xg = hull.xCog;

    const double u = state.surgeSpeed, v = state.swaySpeed, r = state.yawRate;

    // Sway and yaw are coupled through the longitudinal offset of the centre of
    // gravity, so they are solved together rather than sequentially.
    const double a11 = m + my, a12 = xg * m;
    const double a21 = xg * m, a22 = izz + xg * xg * m + jz;
    const double b1 = f.sway - (m + mx) * u * r;
    const double b2 = f.yawMoment - xg * m * u * r;
    const double det = a11 * a22 - a12 * a21;

    const double surgeAccel = (f.surge + (m + my) * v * r + xg * m * r * r) / (m + mx);
    const double swayAccel = (b1 * a22 - a12 * b2) / det;
    const double yawAccel = (a11 * b2 - a21 * b1) / det;

    state.surgeSpeed += surgeAccel * dt;
    state.swaySpeed += swayAccel * dt;
    state.yawRate += yawAccel * dt;
    state.heading += state.yawRate * dt;

    const double ch = std::cos(state.heading), sh = std::sin(state.heading);
    state.x += (state.surgeSpeed * ch - state.swaySpeed * sh) * dt;
    state.y += (state.surgeSpeed * sh + state.swaySpeed * ch) * dt;
}

// ---------------------------------------------------------------------------
// Reference ship
// ---------------------------------------------------------------------------

HullParams kvlcc2Hull() { return HullParams{}; }

PropellerParams kvlcc2Propeller() {
    PropellerParams p;
    p.diameter = 9.86;
    p.pitchRatio = 0.721;
    p.bladeAreaRatio = 0.70;
    p.wakeFraction = 0.35;
    p.thrustDeduction = 0.220;
    return p;
}

RudderParams kvlcc2Rudder() {
    const HullParams h = kvlcc2Hull();
    RudderParams rd;
    rd.area = h.length * h.draft / 54.45;
    rd.aspectRatio = 1.827;
    rd.span = 15.8;
    rd.x = -0.5 * h.length;
    rd.flowStraighteningLever = -0.710 * h.length;
    rd.hullLiftLever = -0.464 * h.length;
    return rd;
}

Manoeuvring kvlcc2() {
    Manoeuvring m;
    m.hull = kvlcc2Hull();
    m.propeller = kvlcc2Propeller();
    m.rudder = kvlcc2Rudder();
    return m;
}

}  // namespace sim
