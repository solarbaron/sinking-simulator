// SPDX-License-Identifier: MIT
//
// Validation of the response-amplitude-operator machinery.
//
// This file has two halves, and they fail in different ways. The harmonic fit is
// pure arithmetic with exact answers, so it is asserted exactly. The RAO itself
// is the output of the whole simulator, so it is asserted against the two
// asymptotes every seakeeping text agrees on and no implementation gets right by
// accident: a ship in a wave far longer than itself rides the surface (heave
// ratio 1, pitch equal to the wave slope), and a ship in a wave far shorter than
// itself barely moves.
//
// The asymptotes matter more than a mid-frequency number would, because they are
// the two places where the answer does not depend on added mass, damping, or any
// coefficient anyone could tune.
#include "engine/sim/rao.hpp"
#include "harness.hpp"

#include <cmath>
#include <cstdio>
#include <vector>

using namespace sim;
using testing::expectNear;
using testing::expectTrue;

namespace {

// A plain box barge: 60 m x 16 m, floating at 4 m. Tessellated along its length,
// because a single long panel under a short wave invents displacement -- see
// testHullMustResolveTheWavelength in test_core.cpp.
Ship makeBarge(int lengthDivisions = 32) {
    Ship s;
    std::vector<Station> stations;
    for (int i = 0; i <= lengthDivisions; ++i) {
        Station station;
        station.x = -30.0 + 60.0 * i / lengthDivisions;
        station.halfBeam = {8.0, 8.0};
        stations.push_back(station);
    }
    s.hull = makeHullFromStations(stations, {0.0, 12.0});
    s.deckEdgeZ = 12.0;
    s.lightshipMass = 60.0 * 16.0 * 4.0 * kRhoSeawater;
    s.lightshipCog = {0, 0, 5.0};
    s.gyradii = {5.0, 16.0, 16.0};
    return s;
}

std::vector<double> sampleCosine(double mean, double amplitude, double omega, double phase,
                                 double dt, int count) {
    std::vector<double> out;
    out.reserve(static_cast<std::size_t>(count));
    for (int i = 0; i < count; ++i)
        out.push_back(mean + amplitude * std::cos(omega * (i * dt) + phase));
    return out;
}

// --- The fit -----------------------------------------------------------------

void testHarmonicFitIsExactOnAPureSinusoid() {
    const double dt = 0.01, omega = 0.8, amplitude = 2.0, phase = 0.6, mean = 3.5;
    const std::vector<double> signal = sampleCosine(mean, amplitude, omega, phase, dt, 3000);
    const HarmonicFit fit = fitHarmonic(signal, dt, omega);

    expectNear("amplitude is recovered exactly", fit.amplitude, amplitude, 1e-9);
    expectNear("phase is recovered exactly", fit.phase, phase, 1e-9);
    expectNear("the mean is separated from the oscillation", fit.mean, mean, 1e-9);
    expectTrue("a pure sinusoid leaves no residual", fit.residual < 1e-9);
}

// The reason this is a least-squares fit and not a DFT bin. A record holding a
// non-integer number of cycles leaks in a DFT, and it leaks *smoothly* -- the
// amplitude comes out a few percent wrong and looks like a physical result.
void testHarmonicFitDoesNotLeakOnAPartialCycle() {
    const double dt = 0.01, omega = 0.8, amplitude = 2.0, phase = -1.1;
    const double period = 2.0 * kPi / omega;
    // 4.37 cycles: deliberately not a whole number, and not close to one.
    const auto count = static_cast<int>(4.37 * period / dt);
    const std::vector<double> signal = sampleCosine(0.0, amplitude, omega, phase, dt, count);

    // Guard against the test being vacuous: confirm the window really is partial.
    const double cycles = count * dt / period;
    expectTrue("the record genuinely holds a partial cycle",
               std::abs(cycles - std::round(cycles)) > 0.2);

    const HarmonicFit fit = fitHarmonic(signal, dt, omega);
    expectNear("a partial window still recovers the amplitude", fit.amplitude, amplitude, 1e-9);
    expectNear("and the phase", fit.phase, phase, 1e-9);
}

// A second harmonic is exactly what a nonlinear response adds, so the fit must
// ignore it in the amplitude and report it in the residual. For a cosine of
// amplitude A2 over a whole number of its own cycles, the RMS is A2/sqrt(2).
//
// Two different exactness claims live here and it is worth keeping them apart.
// Least squares recovers the frequency it is *fitting* exactly on any window --
// that is the previous test. Rejecting a *different* frequency needs the two to
// be orthogonal over the window, which needs a whole number of cycles. So the
// sample count is chosen to make 20 cycles land exactly, rather than truncated
// to whatever `int()` gives.
void testHarmonicFitReportsNonlinearityAsResidual() {
    const double omega = 1.0, a1 = 1.5, a2 = 0.4;
    const double period = 2.0 * kPi / omega;
    const int samplesPerCycle = 1000, cycles = 20;
    const double dt = period / samplesPerCycle;
    const int count = cycles * samplesPerCycle;

    std::vector<double> signal;
    signal.reserve(static_cast<std::size_t>(count));
    for (int i = 0; i < count; ++i) {
        const double t = i * dt;
        signal.push_back(a1 * std::cos(omega * t) + a2 * std::cos(2.0 * omega * t + 0.3));
    }

    const HarmonicFit fit = fitHarmonic(signal, dt, omega);
    expectNear("the fundamental is unaffected by a second harmonic", fit.amplitude, a1, 1e-9);
    expectNear("the harmonic lands in the residual at its RMS", fit.residual,
               a2 / std::sqrt(2.0), 1e-9);

    // And the leakage when the window is *not* a whole number of cycles is real,
    // small, and worth having a number for: a record short by a thousandth of a
    // cycle moves the fundamental in the fifth decimal place. That is why the RAO
    // sweep sizes its window from the response period rather than the wave's.
    std::vector<double> clipped(signal.begin(), signal.end() - 1);
    const HarmonicFit partial = fitHarmonic(clipped, dt, omega);
    const double leakage = std::abs(partial.amplitude - a1) / a1;
    expectTrue("a partial window leaks a little", leakage > 1e-9);
    expectTrue("but only a little", leakage < 1e-4);
}

void testHarmonicFitRefusesDegenerateInput() {
    const std::vector<double> flat(500, 7.25);
    const HarmonicFit constantFit = fitHarmonic(flat, 0.01, 0.9);
    expectTrue("a constant signal has no oscillation", constantFit.amplitude < 1e-9);
    expectNear("and its mean is itself", constantFit.mean, 7.25, 1e-12);

    // At omega = 0 the cosine column *is* the constant column, so the system is
    // singular. Refusing beats returning a large amplitude with a free phase.
    const std::vector<double> signal = sampleCosine(1.0, 2.0, 0.5, 0.0, 0.01, 500);
    const HarmonicFit zeroFit = fitHarmonic(signal, 0.01, 0.0);
    expectTrue("a zero frequency yields no amplitude", zeroFit.amplitude == 0.0);
    expectTrue("an empty record is handled", fitHarmonic({}, 0.01, 1.0).amplitude == 0.0);
}

// --- Encounter frequency -----------------------------------------------------

void testEncounterFrequencySigns() {
    const double omega = 0.9, speed = 8.0;
    expectNear("at rest, encounter frequency is the wave frequency",
               encounterFrequency(omega, 0.0, 0.0), omega, 1e-12);

    const double following = encounterFrequency(omega, speed, 0.0);
    const double head = encounterFrequency(omega, speed, kPi);
    expectTrue("overtaking waves are met less often", following < omega);
    expectTrue("head seas are met more often", head > omega);
    expectNear("the two are symmetric about the wave frequency", following + head, 2.0 * omega,
               1e-12);
    expectNear("beam seas are unaffected by speed", encounterFrequency(omega, speed, kPi / 2.0),
               omega, 1e-12);

    // The magnitude is omega^2 U / g, not omega U / g -- an easy substitution to
    // make and impossible to see in a plot.
    expectNear("the shift is k U, with k from the dispersion relation", omega - following,
               omega * omega * speed / kGravity, 1e-12);
}

// --- The asymptotes ----------------------------------------------------------

// In a wave far longer than the ship, the hull is a float on a slowly tilting
// surface: it rises the full wave height and tilts to the full wave slope. Both
// ratios go to 1, and neither depends on any hydrodynamic coefficient.
void testLongWavesAreRiddenExactly() {
    const Ship barge = makeBarge();
    RaoSettings settings;
    settings.waveAmplitude = 0.4;
    settings.heading = kPi;  // head seas, so pitch is excited
    settings.settleCycles = 12;
    settings.recordCycles = 8;

    // lambda = 2 pi g / omega^2. At omega = 0.32 that is ~600 m against a 60 m
    // ship, so the residual is order (L/lambda)^2, a per cent or so.
    const RaoPoint point = measureRaoAt(barge, 0.32, settings);
    expectTrue("the test wave really is long compared with the ship",
               point.waveLength > 8.0 * 60.0);

    expectNear("heave follows a very long wave one for one", point.heave, 1.0, 0.08);
    expectNear("pitch follows the slope of a very long wave", point.pitch, 1.0, 0.15);
    expectTrue("the response is a clean harmonic, so calling it an RAO is honest",
               point.heaveNonlinearity < 0.15);
}

// The other end: a wave much shorter than the ship averages out along the hull
// and the ship ignores it. Without this, a broken RAO that simply returned 1
// everywhere would pass the long-wave check.
void testShortWavesAreIgnored() {
    const Ship barge = makeBarge(64);
    RaoSettings settings;
    settings.waveAmplitude = 0.25;
    settings.heading = kPi;
    settings.settleCycles = 25;
    settings.recordCycles = 12;

    const RaoPoint point = measureRaoAt(barge, 2.2, settings);
    expectTrue("the test wave really is short compared with the ship",
               point.waveLength < 60.0 / 4.0);
    expectTrue("a wave much shorter than the ship barely lifts it", point.heave < 0.25);
}

// The strongest check available here, because the answer is a closed form, it is
// emergent rather than coded, and no coefficient can be tuned to produce it.
//
// A box barge has constant cross-section, so the Froude-Krylov heave force is the
// pressure integrated along the length: proportional to the integral of cos(k x)
// over [-L/2, L/2], which is sin(kL/2)/(kL/2). That vanishes at kL = 2 pi and
// kL = 4 pi -- wavelengths of exactly L and L/2, where the crest over one half of
// the ship is cancelled by the trough over the other. A real ship shows the same
// notch, blurred by its shape; a box shows it sharply.
//
// If the surface integration ever silently stops resolving the hull, or the
// wave's spatial phase is dropped, these notches are the first thing to go, and
// nothing about the asymptotes would notice.
void testHeaveExcitationCancelsWhenTheWavelengthMatchesTheShip() {
    const Ship barge = makeBarge(64);
    RaoSettings settings;
    settings.waveAmplitude = 0.3;
    settings.heading = kPi;
    settings.settleCycles = 18;
    settings.recordCycles = 12;

    const double length = 60.0;
    // omega for a chosen wavelength, from omega^2 = g k.
    const auto omegaFor = [](double wavelength) {
        return std::sqrt(kGravity * 2.0 * kPi / wavelength);
    };

    const RaoPoint atL = measureRaoAt(barge, omegaFor(length), settings);
    const RaoPoint atHalfL = measureRaoAt(barge, omegaFor(length / 2.0), settings);
    // Two flanking frequencies where the sinc is nowhere near a zero.
    const RaoPoint below = measureRaoAt(barge, omegaFor(length * 1.6), settings);
    const RaoPoint between = measureRaoAt(barge, omegaFor(length * 0.68), settings);

    expectTrue("the wavelength really equals the ship length",
               std::abs(atL.waveLength - length) < 0.5);

    // Vacuity guard: the flanking points must be a real response, or "smaller
    // than the neighbours" would be a statement about noise.
    expectTrue("the flanking frequencies produce a substantial heave",
               below.heave > 0.2 && between.heave > 0.05);

    expectTrue("heave collapses when the wavelength equals the ship length",
               atL.heave < 0.25 * below.heave);
    expectTrue("and again at half the ship length, the second sinc zero",
               atHalfL.heave < 0.5 * between.heave);
}

// A sweep has to be monotone in the right places and it has to be reproducible.
void testSweepIsOrderedAndDeterministic() {
    const Ship barge = makeBarge();
    RaoSettings settings;
    settings.waveAmplitude = 0.4;
    settings.heading = kPi;
    settings.settleCycles = 10;
    settings.recordCycles = 6;

    const std::vector<double> omegas{0.35, 0.55, 1.90};
    const std::vector<RaoPoint> first = measureRao(barge, omegas, settings);
    const std::vector<RaoPoint> second = measureRao(barge, omegas, settings);

    expectTrue("the sweep returns a point per frequency", first.size() == omegas.size());
    for (std::size_t i = 0; i < first.size(); ++i)
        expectTrue("re-running the sweep gives bit-identical heave",
                   first[i].heave == second[i].heave);

    expectTrue("wavelength falls as frequency rises",
               first[0].waveLength > first[1].waveLength &&
                   first[1].waveLength > first[2].waveLength);
    expectTrue("the long wave is ridden more than the short one",
               first[0].heave > first[2].heave);
}

}  // namespace

void runRaoTests() {
    std::printf("\n--- response amplitude operators ---\n");
    testHarmonicFitIsExactOnAPureSinusoid();
    testHarmonicFitDoesNotLeakOnAPartialCycle();
    testHarmonicFitReportsNonlinearityAsResidual();
    testHarmonicFitRefusesDegenerateInput();
    testEncounterFrequencySigns();
    testLongWavesAreRiddenExactly();
    testShortWavesAreIgnored();
    testHeaveExcitationCancelsWhenTheWavelengthMatchesTheShip();
    testSweepIsOrderedAndDeterministic();
}
