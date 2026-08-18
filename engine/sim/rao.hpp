// SPDX-License-Identifier: MIT
// Response amplitude operators, measured from the simulator rather than assumed.
//
// An RAO is the transfer function between a regular wave and a ship's motion:
// heave in metres per metre of wave amplitude, pitch and roll in radians per
// radian of wave slope. It is the quantity model basins publish, so it is the
// one place where this simulator can be checked against the outside world
// instead of against itself.
//
// The measurement is deliberately made the way a basin makes it: run the ship in
// a regular wave, throw away the transient, and fit a harmonic to what is left.
// Nothing here linearises the ship or reads a coefficient out of the solver --
// if the simulation is wrong, the RAO is wrong, which is the entire point.
#pragma once

#include "ship.hpp"

#include <limits>
#include <vector>

namespace sim {

// A least-squares fit of `mean + amplitude * cos(omega * t + phase)` to a sampled
// signal, with t = 0 at the first sample.
//
// Least squares rather than a plain DFT bin because the record is not guaranteed
// to hold a whole number of cycles: a truncated window leaks a real signal into
// a wrong amplitude, and it does so smoothly enough to look like physics.
struct HarmonicFit {
    double amplitude = 0;  // always >= 0
    double phase = 0;      // rad, in (-pi, pi]
    double mean = 0;       // the drift the wave rides on
    double residual = 0;   // RMS of what the fit could not explain

    // **False means the normal equations were refused**, and the four numbers above
    // are the defaults rather than a measurement. Without it a refused fit is
    // indistinguishable from a real one that found no motion: both report an
    // amplitude and a residual of zero, and every consumer that asked
    // `amplitude > 0` read the first as the second. `mean` is the exception and is
    // filled either way, because a mean needs no solve.
    bool fitted = false;
};

HarmonicFit fitHarmonic(const std::vector<double>& samples, double dt, double omega);

// One frequency of a sweep.
struct RaoPoint {
    double omega = 0;           // wave frequency, rad/s
    double encounterOmega = 0;  // what the ship actually feels, rad/s
    double forwardSpeed = 0;    // m/s, measured over the record if under power
    double waveLength = 0;      // m
    double waveAmplitude = 0;   // m, as driven

    double heave = 0;  // m / m
    double pitch = 0;  // rad / rad of wave slope (k * amplitude)
    double roll = 0;   // rad / rad of wave slope

    double heavePhase = 0;  // rad, relative to the wave crest at the origin
    double pitchPhase = 0;
    double rollPhase = 0;

    // Fraction of the recorded motion the single harmonic failed to explain.
    // Large values mean the response is not linear at this amplitude -- the
    // number is still reported, but it is no longer an RAO.
    //
    // **NaN means there is no answer**, which is a different thing from linear and
    // used to be reported as 0.0 -- the most reassuring value in the range, from
    // the one input that means the question was not answered. A threshold test on
    // this field fails for a NaN, which is the safe direction and was not the old
    // one.
    //
    // Two paths reach it and both used to read 0.0. A *fit* that was refused, and
    // a *measurement* that was: `measureRaoAt` returns a default-constructed point
    // when the record is too short to fit at all, so the default itself has to be
    // the honest value rather than the reassuring one.
    double heaveNonlinearity = std::numeric_limits<double>::quiet_NaN();
    double pitchNonlinearity = std::numeric_limits<double>::quiet_NaN();
    double rollNonlinearity = std::numeric_limits<double>::quiet_NaN();
};

struct RaoSettings {
    // Small enough to stay linear, large enough to sit well above the noise a
    // discrete hull integration produces. 0.5 m on a 120 m ship is both.
    double waveAmplitude = 0.5;
    double timestep = 0.02;      // s
    double heading = 0.0;        // rad; 0 = following, pi = head seas

    // Still-water running before the wave is switched on, so a ship with
    // propulsion reaches its own speed rather than being told one.
    //
    // **Set this when the prototype carries propulsion.** It defaults to zero, and
    // at zero a powered ship accelerates *through* her own measurement window: the
    // reported RAO is still fitted at the speed she averaged over the record, so
    // the answer is not wrong, but the window was sized from a speed she never
    // held and the surge is not steady across it.
    double accelerateSeconds = 0.0;

    // Speed to *assume* when the ship has no propulsion attached -- a towed
    // model, or a hull whose machinery is not modelled.
    //
    // When `prototype` does carry propulsion this does not reach the reported
    // answer: the encounter frequency is re-derived from the mean surge over the
    // record, because a requested speed the hull cannot reach would put the fit at
    // a frequency the ship never met. It is not quite unused, though -- it sizes
    // the settle and record windows, which are counted in response periods, so a
    // badly wrong value still costs accuracy through the fit's window length.
    double forwardSpeed = 0.0;   // m/s
    int settleCycles = 15;       // transient discarded before recording
    int recordCycles = 10;       // fitted window

    // Autopilot gains, used only when the ship has propulsion. A hull with real
    // manoeuvring derivatives is frequently directionally unstable, and one that
    // has wandered out of the heading it was given is no longer measuring the
    // RAO that was asked for.
    double headingGain = 3.0;    // rad of rudder per rad of heading error
    double rateGain = 30.0;      // rad of rudder per rad/s of yaw rate
};

// Encounter frequency for a ship making way. omega_e = omega - omega^2 U cos(mu) / g,
// where mu is the angle between the ship's heading and the wave's direction of
// travel. At zero speed this is just omega, but getting the sign wrong is a
// classic error that only shows up in following seas.
double encounterFrequency(double omega, double forwardSpeed, double heading);

// Run `prototype` in a regular wave at each frequency and measure its response.
// The ship is copied per frequency, so the sweep is order-independent and each
// point starts from the same state.
std::vector<RaoPoint> measureRao(const Ship& prototype, const std::vector<double>& omegas,
                                 const RaoSettings& settings = {});

// The single-frequency form, for a targeted check.
RaoPoint measureRaoAt(const Ship& prototype, double omega, const RaoSettings& settings = {});

}  // namespace sim
