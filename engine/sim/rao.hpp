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
    double heaveNonlinearity = 0;
    double pitchNonlinearity = 0;
    double rollNonlinearity = 0;
};

struct RaoSettings {
    // Small enough to stay linear, large enough to sit well above the noise a
    // discrete hull integration produces. 0.5 m on a 120 m ship is both.
    double waveAmplitude = 0.5;
    double timestep = 0.02;      // s
    double heading = 0.0;        // rad; 0 = following, pi = head seas

    // Still-water running before the wave is switched on, so a ship with
    // propulsion reaches its own speed rather than being told one.
    double accelerateSeconds = 0.0;

    // Speed to *assume* when the ship has no propulsion attached -- a towed
    // model, or a hull whose machinery is not modelled. When `prototype` does
    // carry propulsion this is ignored and the speed the ship actually achieves
    // is measured instead, because a requested speed the hull cannot reach would
    // put the fit at a frequency the ship never met.
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
