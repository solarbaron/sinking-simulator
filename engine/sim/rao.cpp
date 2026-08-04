// SPDX-License-Identifier: MIT
#include "rao.hpp"

#include "waves.hpp"

#include <cmath>

namespace sim {
namespace {

// Symmetric 3x3 solve for the normal equations, with a singularity check rather
// than a tolerance. The system goes singular exactly when omega is zero, where
// the cosine column is the constant column; refusing there is better than
// returning a large amplitude and a meaningless phase.
bool solve3x3(const double a[3][3], const double b[3], double x[3]) {
    const double det = a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) -
                       a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0]) +
                       a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
    // Scale the determinant against the matrix's own magnitude so the test means
    // "singular", not "small numbers".
    double scale = 0;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) scale = std::max(scale, std::abs(a[i][j]));
    if (scale <= 0 || std::abs(det) < 1e-12 * scale * scale * scale) return false;

    for (int col = 0; col < 3; ++col) {
        double m[3][3];
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) m[i][j] = (j == col) ? b[i] : a[i][j];
        const double d = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
                         m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
                         m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
        x[col] = d / det;
    }
    return true;
}

// Speed along the ship's own bow, not along world +x. RigidState::velocity is a
// world vector, so its x component is the ship's speed only while the ship
// happens to point that way -- and a hull that has turned reports a surge speed
// that oscillates with its heading while the real one is dead constant. That
// mistake produced a "chaotic" speed trace that was a perfectly steady turn.
double surgeSpeed(const Ship& ship) {
    const Mat3 R = ship.state.orientation.toMat3();
    return dot(ship.state.velocity, R * Vec3{1, 0, 0});
}

// Hold the initial heading with the rudder.
//
// Not a convenience. Real manoeuvring derivatives often describe a
// directionally *unstable* hull -- the reference set here runs straight for
// about 1900 s and then departs into a steady turn -- and a ship that has turned
// out of head seas has stopped measuring the thing that was asked for, quietly
// and while still producing a plausible number. Proportional on heading error,
// derivative on yaw rate; positive rudder swings the bow to port, so correcting
// a bow-to-port error needs negative rudder.
void holdHeading(Ship& ship, double headingGain, double rateGain) {
    if (!ship.propulsion.has_value()) return;
    const Mat3 R = ship.state.orientation.toMat3();
    const Vec3 forward = R * Vec3{1, 0, 0};
    const double heading = std::atan2(forward.y, forward.x);
    const double demand = -headingGain * heading - rateGain * ship.state.angularVelocity.z;
    ship.propulsion->rudderAngle = std::clamp(demand, -35.0 * kPi / 180.0, 35.0 * kPi / 180.0);
}

}  // namespace

HarmonicFit fitHarmonic(const std::vector<double>& samples, double dt, double omega) {
    HarmonicFit fit;
    const std::size_t n = samples.size();
    if (n == 0) return fit;

    double mean = 0;
    for (double v : samples) mean += v;
    mean /= static_cast<double>(n);
    fit.mean = mean;
    if (n < 3 || omega <= 0.0 || dt <= 0.0) return fit;

    // Normal equations for the basis [1, cos(wt), sin(wt)].
    double a[3][3] = {}, b[3] = {};
    for (std::size_t i = 0; i < n; ++i) {
        const double t = static_cast<double>(i) * dt;
        const double basis[3] = {1.0, std::cos(omega * t), std::sin(omega * t)};
        for (int r = 0; r < 3; ++r) {
            for (int c = 0; c < 3; ++c) a[r][c] += basis[r] * basis[c];
            b[r] += basis[r] * samples[i];
        }
    }

    double x[3];
    if (!solve3x3(a, b, x)) return fit;

    fit.mean = x[0];
    // y = m + A cos(wt) + B sin(wt) = m + amplitude cos(wt + phase), so
    // A = amplitude cos(phase) and B = -amplitude sin(phase).
    fit.amplitude = std::hypot(x[1], x[2]);
    fit.phase = std::atan2(-x[2], x[1]);

    double sumSquares = 0;
    for (std::size_t i = 0; i < n; ++i) {
        const double t = static_cast<double>(i) * dt;
        const double model = fit.mean + fit.amplitude * std::cos(omega * t + fit.phase);
        const double e = samples[i] - model;
        sumSquares += e * e;
    }
    fit.residual = std::sqrt(sumSquares / static_cast<double>(n));
    return fit;
}

double encounterFrequency(double omega, double forwardSpeed, double heading) {
    // omega_e = omega - k U cos(mu), k = omega^2 / g. Heading 0 means the wave
    // travels the way the ship points, so it overtakes and the encounter
    // frequency drops; head seas raise it.
    return omega - omega * omega * forwardSpeed * std::cos(heading) / kGravity;
}

RaoPoint measureRaoAt(const Ship& prototype, double omega, const RaoSettings& settings) {
    RaoPoint point;
    point.omega = omega;
    point.waveAmplitude = settings.waveAmplitude;
    if (omega <= 0.0 || settings.timestep <= 0.0) return point;

    const double k = deepWaterWavenumber(omega);
    point.waveLength = 2.0 * kPi / k;

    Ship ship = prototype;
    ship.initialise(0.0);

    const double dt = settings.timestep;

    // Let a powered ship reach its own speed in still water before the wave
    // arrives. A speed the hull cannot actually hold would put the fit at a
    // frequency the ship never met, which is why this is run rather than asked
    // for: with propulsion attached, settings.forwardSpeed is ignored entirely.
    double speed = settings.forwardSpeed;
    if (ship.propulsion.has_value() && settings.accelerateSeconds > 0.0) {
        const Sea still(0.0);
        const auto steps = static_cast<long long>(settings.accelerateSeconds / dt);
        for (long long i = 0; i < steps; ++i) {
            holdHeading(ship, settings.headingGain, settings.rateGain);
            ship.step(dt, still);
        }
        speed = surgeSpeed(ship);
    }

    point.encounterOmega = encounterFrequency(omega, speed, settings.heading);
    double omegaResponse = std::abs(point.encounterOmega);
    if (omegaResponse <= 0.0) return point;

    const WaveField field = WaveField::regular(settings.waveAmplitude, omega, settings.heading);
    Sea sea;
    sea.waves = &field;
    const double responsePeriod = 2.0 * kPi / omegaResponse;
    const auto settleSteps =
        static_cast<long long>(settings.settleCycles * responsePeriod / dt);
    const auto recordSteps =
        static_cast<long long>(settings.recordCycles * responsePeriod / dt);
    if (recordSteps < 8) return point;

    std::vector<double> heave, pitch, roll, wave;
    heave.reserve(static_cast<std::size_t>(recordSteps));
    pitch.reserve(static_cast<std::size_t>(recordSteps));
    roll.reserve(static_cast<std::size_t>(recordSteps));
    wave.reserve(static_cast<std::size_t>(recordSteps));
    double speedSum = 0;

    for (long long i = 0; i < settleSteps + recordSteps; ++i) {
        sea.time = static_cast<double>(i) * dt;
        holdHeading(ship, settings.headingGain, settings.rateGain);
        ship.step(dt, sea);
        if (i < settleSteps) continue;

        // Heel and trim straight from the orientation, not via diagnostics().
        // Diagnostics repeats the whole wavy buoyancy integration and copies the
        // hull to do it -- around 960 wave-field queries, which is the entire
        // cost of a tick -- to deliver, among much else, two angles that are a
        // pure function of the rotation matrix. Sampling every step would have
        // doubled the price of a sweep for nothing.
        double heelRad = 0, trimRad = 0;
        heelTrimFromRotation(ship.state.orientation.toMat3(), heelRad, trimRad);
        heave.push_back(ship.state.position.z);
        pitch.push_back(trimRad);
        roll.push_back(heelRad);
        speedSum += surgeSpeed(ship);
        // The driving wave, sampled **at the ship** over the very same window.
        // Fitting it rather than deriving its phase analytically means the
        // reported motion phases cannot pick up a sign error from the
        // book-keeping of when the window started.
        //
        // At the ship, not at the origin: a moving hull meets the wave at the
        // encounter frequency, while the surface at a fixed point still
        // oscillates at omega. Fitting the origin's history at omega_e would
        // return a small, meaningless amplitude and every RAO would be
        // normalised by it. For a ship at rest the two are the same sample.
        wave.push_back(field.elevation(ship.state.position.x, ship.state.position.y, sea.time));
    }

    // Added resistance in waves slows a ship, so the speed over the record is not
    // the speed it settled at in flat water. Re-derive the encounter frequency
    // from what actually happened and fit at that.
    point.forwardSpeed = recordSteps > 0 ? speedSum / static_cast<double>(recordSteps) : speed;
    if (ship.propulsion.has_value()) {
        point.encounterOmega = encounterFrequency(omega, point.forwardSpeed, settings.heading);
        omegaResponse = std::abs(point.encounterOmega);
        if (omegaResponse <= 0.0) return point;
    }

    const HarmonicFit waveFit = fitHarmonic(wave, dt, omegaResponse);
    const HarmonicFit heaveFit = fitHarmonic(heave, dt, omegaResponse);
    const HarmonicFit pitchFit = fitHarmonic(pitch, dt, omegaResponse);
    const HarmonicFit rollFit = fitHarmonic(roll, dt, omegaResponse);

    // Normalise against the wave the ship actually saw, not the one requested.
    // They differ if the record is short or the response is not at omega, and
    // using the requested amplitude would quietly absorb that error into the RAO.
    const double amplitude = waveFit.amplitude > 0 ? waveFit.amplitude : settings.waveAmplitude;
    const double slope = k * amplitude;

    point.heave = heaveFit.amplitude / amplitude;
    point.pitch = slope > 0 ? pitchFit.amplitude / slope : 0.0;
    point.roll = slope > 0 ? rollFit.amplitude / slope : 0.0;

    const auto wrap = [](double a) {
        while (a > kPi) a -= 2.0 * kPi;
        while (a <= -kPi) a += 2.0 * kPi;
        return a;
    };
    point.heavePhase = wrap(heaveFit.phase - waveFit.phase);
    point.pitchPhase = wrap(pitchFit.phase - waveFit.phase);
    point.rollPhase = wrap(rollFit.phase - waveFit.phase);

    // How much of the motion the single harmonic could not explain. A linear
    // response is a pure sinusoid at the encounter frequency; anything else means
    // the answer is a number, but not an RAO.
    const auto nonlinearity = [](const HarmonicFit& f) {
        return f.amplitude > 0 ? f.residual / f.amplitude : 0.0;
    };
    point.heaveNonlinearity = nonlinearity(heaveFit);
    point.pitchNonlinearity = nonlinearity(pitchFit);
    point.rollNonlinearity = nonlinearity(rollFit);
    return point;
}

std::vector<RaoPoint> measureRao(const Ship& prototype, const std::vector<double>& omegas,
                                 const RaoSettings& settings) {
    std::vector<RaoPoint> points;
    points.reserve(omegas.size());
    for (double omega : omegas) points.push_back(measureRaoAt(prototype, omega, settings));
    return points;
}

}  // namespace sim
