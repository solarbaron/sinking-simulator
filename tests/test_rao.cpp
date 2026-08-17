// SPDX-License-Identifier: MIT
//
// Validation of the response-amplitude-operator machinery and of the radiation
// coupling that feeds it.
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

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

using namespace sim;
using testing::expectEqual;
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


// --- Radiation coupling ------------------------------------------------------

// Stations are read off the hull mesh by clipping it into slabs, so a box barge
// has an exact answer: every section is a rectangle of the full beam and draft,
// and its area coefficient is exactly 1.
void testStationsFromAMeshAreExactForABox() {
    const Ship box = makeBarge(32);
    const RadiationHull sections = radiationHullFromMesh(box.hull, 4.0, 11);

    expectEqual("every requested station was produced",
                static_cast<long long>(sections.stations.size()), 11LL);
    expectNear("draft comes from the mesh", sections.draft, 4.0, 1e-12);

    double worstBeam = 0, worstDraft = 0, worstSigma = 0;
    for (const RadiationStation& st : sections.stations) {
        worstBeam = std::max(worstBeam, std::abs(st.beam - 16.0));
        worstDraft = std::max(worstDraft, std::abs(st.draft - 4.0));
        worstSigma = std::max(worstSigma, std::abs(st.areaCoefficient - 1.0));
    }
    expectTrue("a box has full beam at every station", worstBeam < 1e-9);
    expectTrue("and full draft", worstDraft < 1e-9);
    // Sectional area over B*T. Anything but 1 means the slab volume, the beam or
    // the draft disagree with each other.
    expectTrue("and a rectangular area coefficient of exactly 1", worstSigma < 1e-12);

    // Stations must be inset from the ends: a slab at the stem collects no volume
    // and its area coefficient is whatever the rounding says.
    expectTrue("stations are inset from the extreme ends",
               sections.stations.front().x > -30.0 && sections.stations.back().x < 30.0);
}

// A barge with radiation attached, built once: the 2D solve over every station
// and frequency costs seconds, and this file needs it three times.
const Ship& radiatingBarge() {
    static const Ship ship = [] {
        Ship s = makeBarge(32);
        s.initialise(0.0);
        s.attachRadiation(4.0, 9);
        return s;
    }();
    return ship;
}

// The regression this exists for: an earlier version fitted the retardation
// function over 1.3 s of a 20 s decay, which places the state-space poles near
// zero. The model then *integrates* velocity instead of damping it, and the ship
// reached NaN in five steps. A free decay that stays bounded and shrinks is the
// direct statement that the memory term is a damper.
void testRadiationMemoryDampsRatherThanDiverges() {
    Ship ship = radiatingBarge();
    expectTrue("radiation is attached", ship.radiation.has_value());
    expectTrue("some entries got state-space models", ship.radiation->modelCount() > 0);

    const Sea still(0.0);
    const double datum = ship.state.position.z;
    const double release = 0.5;
    ship.state.position.z += release;

    double earlyPeak = 0, latePeak = 0, worst = 0;
    bool finite = true;
    for (int i = 1; i <= 6000; ++i) {          // 30 s at dt = 0.005
        ship.step(0.005, still);
        const double z = std::abs(ship.state.position.z - datum);
        finite = finite && std::isfinite(z);   // asserted once, not 6000 times
        worst = std::max(worst, z);
        if (i > 400 && i <= 1600) earlyPeak = std::max(earlyPeak, z);
        if (i > 4800) latePeak = std::max(latePeak, z);
    }
    expectTrue("heave stays finite for the whole run", finite);

    // A diverging model blows past the release amplitude immediately; a damper
    // cannot exceed it by more than the integrator's own overshoot.
    expectTrue("the motion never exceeds what it was released from", worst < 1.2 * release);
    expectTrue("the oscillation is real before it decays", earlyPeak > 0.05 * release);
    expectTrue("and has decayed by the end", latePeak < 0.5 * earlyPeak);
}

// The strongest check on the coupling: the time-domain Cummins model must agree
// with the frequency-domain table it was built from.
//
// Cummins splits radiation into an instantaneous A_inf plus a memory
// convolution, and the memory carries the whole frequency dependence. So a free
// decay settling at omega_d must behave as though its added mass were A(omega_d)
// from the table -- *not* A_inf. Those differ by 60% here, which is what makes
// this a test rather than a coincidence: getting A_inf into the mass matrix and
// dropping the memory term entirely would fail it badly.
void testFreeDecayAgreesWithTheFrequencyDomainTable() {
    Ship ship = radiatingBarge();
    const RadiationHull sections = radiationHullFromMesh(ship.hull, 4.0, 9);
    const RadiationTable table = stripTheoryTable(sections, radiationFrequencyGrid(0.2, 2.5, 40));

    const double mass = 60.0 * 16.0 * 4.0 * kRhoSeawater;
    const double stiffness = kRhoSeawater * kGravity * 60.0 * 16.0;

    const Sea still(0.0);
    const double datum = ship.state.position.z;
    ship.state.position.z += 0.5;

    // Upward zero crossings of the heave deviation, linearly interpolated.
    std::vector<double> crossings;
    double previous = 0.5;
    const double dt = 0.005;
    for (int i = 1; i <= 8000; ++i) {
        ship.step(dt, still);
        const double z = ship.state.position.z - datum;
        if (previous < 0 && z >= 0)
            crossings.push_back(i * dt - dt * z / (z - previous));
        previous = z;
    }
    expectTrue("the decay produced enough cycles to time", crossings.size() >= 4);
    if (crossings.size() < 4) return;

    double total = 0;
    for (std::size_t i = 2; i < crossings.size(); ++i) total += crossings[i] - crossings[i - 1];
    const double period = total / static_cast<double>(crossings.size() - 2);
    const double omegaDamped = 2.0 * kPi / period;

    // What added mass that period implies, against what the table says at that
    // frequency.
    const double implied = stiffness * (period / (2.0 * kPi)) * (period / (2.0 * kPi)) - mass;
    double fromTable = 0;
    for (int i = 0; i + 1 < table.size(); ++i)
        if (table.omega[i] <= omegaDamped && omegaDamped <= table.omega[i + 1]) {
            const double f = (omegaDamped - table.omega[i]) /
                             (table.omega[i + 1] - table.omega[i]);
            fromTable = table.addedMass[static_cast<std::size_t>(i)][2][2] * (1 - f) +
                        table.addedMass[static_cast<std::size_t>(i) + 1][2][2] * f;
        }
    expectTrue("the decay frequency lies inside the table", fromTable > 0);

    expectNear("the free decay sees A(omega_d) from the table", implied, fromTable,
               0.15 * fromTable);

    // The guard that makes the above meaningful. If the memory term were missing
    // and only A_inf reached the mass matrix, `implied` would sit at A_inf, and a
    // 15% tolerance on the wrong quantity would still look like agreement unless
    // the two are known to be far apart.
    const double aInf = table.addedMassInfinite[2][2];
    expectTrue("A_inf and A(omega_d) are far enough apart for this to discriminate",
               std::abs(aInf - fromTable) > 0.4 * fromTable);
    expectTrue("and the decay follows A(omega_d), not A_inf",
               std::abs(implied - fromTable) < std::abs(implied - aInf));
}


// --- Propulsion coupling -----------------------------------------------------

// A manoeuvring set scaled to the barge. The non-dimensional derivatives stay
// KVLCC2's: wrong in detail for a box, right in order of magnitude, and not what
// is under test here -- these tests ask whether the coupling moves the ship and
// whether the sea it meets is the right one, not whether a barge has a tanker's
// hydrodynamics.
Ship poweredBarge(double revsPerSecond) {
    Ship s = makeBarge(32);
    Manoeuvring m = kvlcc2();
    m.hull.length = 60.0;
    m.hull.beam = 16.0;
    m.hull.draft = 4.0;
    m.hull.blockCoefficient = 1.0;
    m.hull.xCog = 0.0;
    m.propeller.diameter = 2.5;
    m.rudder.area = 3.0;
    m.rudder.span = 2.2;
    m.rudder.x = -28.0;
    m.rudder.flowStraighteningLever = -42.6;   // scaled from -227.2 by 60/320
    m.rudder.hullLiftLever = -27.8;            // scaled from -148.5
    m.revsPerSecond = revsPerSecond;
    s.propulsion = m;
    s.initialise(0.0);
    return s;
}

// Speed along the ship's own bow. RigidState::velocity is a *world* vector, so
// its x component is the ship's speed only while the ship still points that way.
// Reading it directly made a perfectly steady turn look like a chaotic speed
// trace, because the world-frame component oscillates with heading while the
// real surge speed is dead constant.
double surgeOf(const Ship& ship) {
    const Mat3 R = ship.state.orientation.toMat3();
    return dot(ship.state.velocity, R * Vec3{1, 0, 0});
}

double runToSpeed(Ship& ship, double seconds, double dt = 0.02) {
    const Sea still(0.0);
    const auto steps = static_cast<long long>(seconds / dt);
    for (long long i = 0; i < steps; ++i) ship.step(dt, still);
    return surgeOf(ship);
}

void testPropulsionDrivesTheShip() {
    double previous = 0;
    for (double revs : {0.5, 1.0, 1.5}) {
        Ship ship = poweredBarge(revs);
        const double speed = runToSpeed(ship, 400.0);
        expectTrue("the ship makes way ahead", speed > previous + 0.05);
        previous = speed;
    }
    expectTrue("and reaches a useful speed, not a crawl", previous > 0.5);

    // Steady state means thrust balances resistance, so the speed stops changing.
    //
    // The window matters and was got wrong first time: a 3.9e6 kg hull has an
    // acceleration time constant of *minutes*. Measured, this one is at 61% of
    // its final speed after 300 s and within 0.3% after 1200 s, so asking
    // whether it had converged at 400 s was asking the wrong question of correct
    // behaviour.
    Ship steady = poweredBarge(1.5);
    const double first = runToSpeed(steady, 1200.0);
    const double second = runToSpeed(steady, 300.0);
    expectTrue("speed settles rather than climbing without limit",
               std::abs(second - first) < 0.02 * first);

    // Astern rotation must drive the ship astern, not merely slow it. This is the
    // case a crash stop needs and the one a first-quadrant propeller model gets
    // wrong while looking entirely plausible ahead.
    Ship astern = poweredBarge(-1.5);
    expectTrue("astern revolutions drive the ship astern", runToSpeed(astern, 400.0) < -0.1);
}

// The rudder convention is stated in propulsion.hpp: positive delta swings the
// bow to port, which is positive yaw about +z. A sign error here turns the ship
// the wrong way, which looks like steering and is not.
void testRudderTurnsTowardsTheCommandedSide() {
    Ship port = poweredBarge(1.5);
    runToSpeed(port, 300.0);
    Ship starboard = port;
    port.propulsion->rudderAngle = 20.0 * kPi / 180.0;
    starboard.propulsion->rudderAngle = -20.0 * kPi / 180.0;

    const Sea still(0.0);
    for (int i = 0; i < 3000; ++i) {         // 60 s of helm
        port.step(0.02, still);
        starboard.step(0.02, still);
    }
    expectTrue("positive rudder swings the bow to port", port.state.angularVelocity.z > 1e-4);
    expectTrue("negative rudder swings it to starboard",
               starboard.state.angularVelocity.z < -1e-4);
    // Port/starboard symmetry: a symmetric hull must turn equally hard each way.
    expectNear("the two turns are mirror images", port.state.angularVelocity.z,
               -starboard.state.angularVelocity.z,
               0.05 * std::abs(port.state.angularVelocity.z));
}

// The milestone check, and the reason forward speed had to exist at all.
//
// Nothing imposes the encounter frequency. The surface is cos(k x - omega t + phi)
// and the ship's x is moving, so a hull under way meets waves at |omega - k U|
// purely because it translates through a spatially varying field. If that does
// not fall out, comparing against published RAOs -- which are all tabulated
// against encounter frequency at a Froude number -- would be meaningless.
//
// The response frequency is measured *independently*, by scanning the harmonic
// fit for the frequency that leaves least unexplained. Reading encounterOmega
// back out of measureRaoAt() would only confirm that a formula was evaluated.
void testEncounterFrequencyEmergesFromMovingThroughTheWave() {
    Ship ship = poweredBarge(2.0);
    runToSpeed(ship, 400.0);

    const double omega = 0.7;
    const double dt = 0.02;
    const WaveField field = WaveField::regular(0.3, omega, kPi);   // head seas
    Sea sea;
    sea.waves = &field;

    std::vector<double> heave;
    double speedSum = 0;
    int samples = 0;
    for (int i = 1; i <= 15000; ++i) {       // 300 s, first third discarded
        sea.time = i * dt;
        ship.step(dt, sea);
        if (i > 5000) {
            heave.push_back(ship.state.position.z);
            speedSum += surgeOf(ship);
            ++samples;
        }
    }
    const double speed = speedSum / samples;
    expectTrue("the ship was genuinely under way", speed > 0.5);

    double bestOmega = 0, bestResidual = 1e30;
    for (double w = 0.3; w <= 1.6; w += 0.001) {
        const HarmonicFit fit = fitHarmonic(heave, dt, w);
        if (fit.residual < bestResidual) {
            bestResidual = fit.residual;
            bestOmega = w;
        }
    }

    const double predicted = encounterFrequency(omega, speed, kPi);

    // Vacuity guard, and the whole point: at rest omega_e is omega, so a test
    // that merely found "some frequency near omega" would pass on a ship that
    // never moved. The shift has to be far larger than the search resolution.
    expectTrue("the Doppler shift is large enough to be the thing measured",
               predicted - omega > 0.05);
    expectNear("a moving ship meets waves at the encounter frequency", bestOmega, predicted,
               0.02 * predicted);
    expectTrue("and not at the wave's own frequency",
               std::abs(bestOmega - omega) > 0.5 * (predicted - omega));
}

// measureRaoAt must run the ship up to speed rather than take one on trust,
// report what it actually achieved, and keep it pointed where it was aimed.
//
// The acceleration window is deliberately past the point where this hull departs
// into a turn unsteered -- measured at about 1900 s. That makes the check
// discriminating rather than decorative: without the autopilot the ship ends up
// in a steady turn at 1.05 m/s of surge instead of running straight at 1.99, and
// every RAO measured from it would be at the wrong speed, the wrong heading, and
// entirely plausible-looking.
void testRaoMeasuresTheSpeedItActuallyReached() {
    const Ship driven = poweredBarge(1.5);
    RaoSettings settings;
    settings.waveAmplitude = 0.3;
    settings.heading = kPi;
    settings.accelerateSeconds = 2200.0;
    settings.forwardSpeed = 99.0;   // absurd, and must be ignored
    settings.settleCycles = 8;
    settings.recordCycles = 6;

    const RaoPoint point = measureRaoAt(driven, 0.7, settings);
    expectTrue("the measured speed is the ship's own, not the requested one",
               point.forwardSpeed > 0.5 && point.forwardSpeed < 20.0);
    expectTrue("the ship was held straight rather than left to fall into a turn",
               point.forwardSpeed > 1.7);
    expectTrue("head seas raise the encounter frequency", point.encounterOmega > 0.7);
    expectNear("and it follows from the speed that was measured", point.encounterOmega,
               encounterFrequency(0.7, point.forwardSpeed, kPi), 1e-9);
    expectTrue("the response is still a usable RAO", point.heave > 0.05 && point.heave < 3.0);
}

// =============================================================================
// --- Viscous roll damping, and the frame the roll added inertia lives in ------
//
// Two things are under test here and they fail in different ways.
//
// The *frame transfer* is arithmetic with exact answers, so it is asserted
// exactly: strip theory assembles A_inf about the body-frame origin -- midship
// on the baseline -- while the integrator adds it to an inertia about the centre
// of gravity, and the two differ by the sway-roll coupling twice over. The sign
// of that transfer is the whole question, and it is pinned by a section whose
// roll added mass is analytically zero rather than by a convention read out of a
// textbook.
//
// The *damping* is Ikeda's empirical method driving the real integrator, so it is
// asserted against the closed form that says what a linear damping coefficient
// means: log decrement 2 pi zeta / sqrt(1 - zeta^2). What makes that a test and
// not a tautology is that B44 is strongly amplitude-dependent, so the decrement
// has to change through a single decay in the way the coefficient says it will --
// and a linear damper, run through the same code, has to hold it constant.
// =============================================================================

// A ship-shaped hull for the roll work: 120 x 22 m at 6 m draft, with rise of
// floor and fine ends.
//
// The box barge above cannot be used. Ikeda's regressions are fitted over
// 0.5 <= Cb <= 0.85 and 0.9 <= Cm <= 0.99 and a box is 1.0 on both, so its eddy
// coefficient would be pure extrapolation and its default bilge radius exactly
// zero. This hull comes out at Cb 0.742, Cm 0.910, B/d 3.67 and OG/d -0.20,
// inside every published bound -- which testRollDampingFormComesOffTheHull
// asserts rather than assumes.
constexpr double kRollWaterline = 6.0;
constexpr double kRollKeelLength = 40.0;
constexpr double kRollKeelBreadth = 1.0;
constexpr double kRollKg = 7.2;          // puts the roll axis 1.2 m above the waterline

Ship rollShip() {
    const std::vector<double> waterlines{0.0, 0.8, 1.6, 2.6, 3.6, 4.8, 6.0, 8.0, 11.0, 14.0};
    std::vector<Station> stations;
    for (int i = 0; i <= 40; ++i) {
        Station s;
        s.x = -60.0 + 120.0 * i / 40.0;
        const double u = std::abs(s.x) / 60.0;
        const double fx = u <= 0.35 ? 1.0 : 1.0 - 0.85 * std::pow((u - 0.35) / 0.65, 2.0);
        for (double z : waterlines) {
            const double fz = z >= 3.0 ? 1.0 : 0.55 + 0.45 * std::pow(z / 3.0, 0.6);
            s.halfBeam.push_back(11.0 * fx * fz);
        }
        stations.push_back(s);
    }
    Ship ship;
    ship.hull = makeHullFromStations(stations, waterlines);
    ship.deckEdgeZ = 11.0;
    const RollDampingHull form = rollDampingHullFromMesh(ship.hull, kRollWaterline);
    ship.lightshipMass = form.blockCoeff * form.lengthPp * form.beam * form.draft * kRhoSeawater;
    ship.lightshipCog = {0, 0, kRollKg};
    ship.gyradii = {0.35 * form.beam, 0.25 * form.lengthPp, 0.26 * form.lengthPp};
    ship.initialise(0.0);
    return ship;
}

struct RollDecay {
    std::vector<double> peaks;   // rad, successive maxima of the heel angle
    std::vector<double> times;   // s

    double period() const {
        return times.size() >= 2 ? (times.back() - times.front()) /
                                       static_cast<double>(times.size() - 1)
                                 : 0.0;
    }
    double decrement(std::size_t i) const { return std::log(peaks[i - 1] / peaks[i]); }
};

// Release the ship at a heel angle and record the successive maxima.
//
// dt is 0.02 s -- 600 steps per roll period. Checked against 0.01 and 0.005: the
// measured decrement moves in the fourth decimal place, so nothing below is a
// statement about the integrator's step size.
RollDecay rollDecay(Ship& ship, double heelDeg, double seconds, double dt = 0.02) {
    const Sea still(0.0);
    ship.state.orientation = Quat::fromAxisAngle(Vec3{1, 0, 0}, heelDeg * kDegToRad);
    RollDecay out;
    double y0 = 0, y1 = 0;
    const auto steps = static_cast<int>(seconds / dt);
    for (int i = 0; i < steps; ++i) {
        ship.step(dt, still);
        double heel = 0, trim = 0;
        heelTrimFromRotation(ship.state.orientation.toMat3(), heel, trim);
        // Parabolic interpolation through three samples, as in
        // test_roll_damping.cpp: taking the nearest step instead biases every
        // peak downward and quietly inflates the decrement.
        if (i >= 2 && y1 > y0 && y1 > heel) {
            const double denom = y0 - 2.0 * y1 + heel;
            const double shift = denom != 0 ? 0.5 * (y0 - heel) / denom : 0.0;
            out.peaks.push_back(y1 - 0.25 * (y0 - heel) * shift);
            out.times.push_back((i - 1 + shift) * dt);
        }
        y0 = y1;
        y1 = heel;
    }
    return out;
}

// --- The frame transfer ------------------------------------------------------

// transferAddedMass() is the parallel-axis theorem wearing a 6x6, so a rigid
// body's own mass matrix must come back out of it obeying the theorem exactly.
// Nothing hydrodynamic is involved; this is the arithmetic, checked before the
// physics leans on it.
void testAddedMassTransferObeysTheParallelAxisTheorem() {
    const double m = 1234.5;
    const Vec3 offset{2.0, -3.0, 5.0};
    const double ixx = 700.0, iyy = 900.0, izz = 1100.0;

    Matrix6 aboutCentre{};
    for (std::size_t i = 0; i < 3; ++i) aboutCentre[i][i] = m;
    aboutCentre[3][3] = ixx;
    aboutCentre[4][4] = iyy;
    aboutCentre[5][5] = izz;

    const Matrix6 moved = transferAddedMass(aboutCentre, offset);
    const double d2 = dot(offset, offset);
    expectNear("Ixx picks up m (dy^2 + dz^2)", moved[3][3], ixx + m * (d2 - offset.x * offset.x),
               1e-9);
    expectNear("Iyy likewise", moved[4][4], iyy + m * (d2 - offset.y * offset.y), 1e-9);
    expectNear("Izz likewise", moved[5][5], izz + m * (d2 - offset.z * offset.z), 1e-9);
    expectNear("and the products of inertia appear", moved[3][4], -m * offset.x * offset.y, 1e-9);

    // Translational added mass is reference-point independent. This is why heave
    // was never affected by the mismatch the transfer exists to fix.
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j)
            expectTrue("the translational block is untouched by the transfer",
                       std::abs(moved[i][j] - aboutCentre[i][j]) < 1e-9);

    // Moving back must land where it started, or the transform is not a group
    // action and one of the two directions is wrong.
    const Matrix6 back = transferAddedMass(moved, Vec3{0, 0, 0} - offset);
    double worst = 0;
    for (std::size_t i = 0; i < 6; ++i)
        for (std::size_t j = 0; j < 6; ++j)
            worst = std::max(worst, std::abs(back[i][j] - aboutCentre[i][j]));
    expectTrue("transferring there and back is the identity", worst < 1e-9);
}

// **The check that pins the sign**, and the reason no textbook convention had to
// be trusted.
//
// Rotate a semicircle about the centre of its own circle and no point of the
// contour moves normal to itself: the radiation potential is identically zero, so
// a44 and a24 about that centre are exactly zero at every frequency. Put the
// centre on the waterline and everything strip theory then reports for the roll
// mode about the baseline is *pure transfer* -- a44 = T^2 a22 and a24 = -T a22 --
// and transferring it back up by the draft has to annihilate it.
//
// The opposite sign does not give a small error there. It gives 4 T^2 a22.
void testSemicircularSectionsPinTheTransferSign() {
    const double radius = 3.0;
    RadiationHull hull;
    hull.draft = radius;
    hull.panelsPerHalfSection = 64;
    const double length = 30.0;
    for (int i = 0; i < 9; ++i) {
        RadiationStation st;
        st.x = -0.5 * length + length * i / 8.0;
        st.beam = 2.0 * radius;
        st.draft = radius;
        st.areaCoefficient = kPi / 4.0;   // exactly a semicircle
        hull.stations.push_back(st);
    }
    const RadiationTable table = stripTheoryTable(hull, {0.8});
    const Matrix6& baseline = table.addedMass[0];

    expectTrue("the semicircular hull has real sway added mass", baseline[1][1] > 0);
    expectNear("roll about the baseline is the sway added mass on a draft lever",
               baseline[3][3], radius * radius * baseline[1][1], 1e-6 * baseline[3][3]);
    expectNear("and the coupling is that lever once", baseline[1][3], -radius * baseline[1][1],
               1e-6 * std::abs(baseline[1][3]));

    // Vacuity guard: if the roll block were already small there would be nothing
    // for the transfer to annihilate.
    const double scale = baseline[3][3];
    expectTrue("there is a substantial roll block to cancel",
               scale > 0.5 * radius * radius * baseline[1][1]);

    const Matrix6 aboutCentre = transferAddedMass(baseline, Vec3{0, 0, radius});
    expectTrue("transferring to the circle's centre annihilates the roll added mass",
               std::abs(aboutCentre[3][3]) < 1e-6 * scale);
    expectTrue("and the sway-roll coupling with it",
               std::abs(aboutCentre[1][3]) < 1e-6 * radius * baseline[1][1]);

    // The wrong sign, spelled out, so the tolerance above is known to be
    // discriminating rather than merely tight.
    const double wrongSign =
        baseline[3][3] - 2.0 * radius * baseline[1][3] + radius * radius * baseline[1][1];
    expectTrue("the opposite sign would leave four drafts' worth of it behind",
               wrongSign > 3.5 * scale);
}

// What the transfer is worth on a real hull, and why only roll cares.
void testRollAddedInertiaIsReferredToTheCentreOfGravity() {
    const Ship& ship = radiatingBarge();
    const Vec3 cog = ship.diagnostics(Sea(0.0)).centreOfGravity;
    expectTrue("the barge's centre of gravity is well above the baseline", cog.z > 3.0);
    expectTrue("the barge is carrying a radiation model", ship.radiation.has_value());

    // Straight off the attached model, which is the matrix the integrator is
    // actually adding to the inertia -- no need to rebuild the table to ask what
    // A_inf is, and this way there is no second table to disagree with the first.
    const Matrix6& origin = ship.radiation->addedMassInfinite();
    const Matrix6 aboutCog = transferAddedMass(origin, cog);

    expectTrue("the transfer is not a rounding correction in roll",
               aboutCog[3][3] < 0.7 * origin[3][3]);
    expectTrue("but it is still positive, as an added inertia must be", aboutCog[3][3] > 0);
    expectNear("heave added mass does not move with the reference point", aboutCog[2][2],
               origin[2][2], 1e-9 * origin[2][2]);
    // Pitch moves only through the *longitudinal* offset of the cog, because
    // strip theory's surge added mass is a structural zero; this barge's cog is
    // on the midship section, so it does not move at all.
    expectNear("pitch is untouched for a cog on the midship section", aboutCog[4][4],
               origin[4][4], 1e-9 * origin[4][4]);

    // Referring the matrix to the cog also shrinks the coupling the integrator
    // has no way to carry, because the cog sits near the height at which sway and
    // roll decouple. That is a real improvement in the approximation, not a
    // cosmetic one.
    expectTrue("the sway-roll coupling the integrator must drop shrinks with the transfer",
               std::abs(aboutCog[1][3]) < 0.5 * std::abs(origin[1][3]));

    // The height at which sway and roll decouple is a property of the sections,
    // and it has to land on the ship rather than under her keel. This is the same
    // statement the semicircle makes, on a hull with no closed form.
    const double rollCentre = -origin[1][3] / origin[1][1];
    expectTrue("the added-mass roll centre sits near the waterline, not below the keel",
               rollCentre > 0.0 && rollCentre < 12.0);
}

// The dynamic statement, and the roll counterpart of
// testFreeDecayAgreesWithTheFrequencyDomainTable: a roll free decay settling at
// omega_d must behave as though its added roll inertia were A44(omega_d) *about
// the centre of gravity*. About the baseline origin that number is three times
// larger, so this discriminates between the two frames rather than merely
// confirming that some added inertia is present.
void testRollDecaySeesAddedInertiaAboutTheCentreOfGravity() {
    Ship ship = radiatingBarge();
    // The comparison table only has to answer "what is A44 near the decay
    // frequency", so it is built on a grid bracketing it rather than over the
    // whole seakeeping band -- same hull, same panel count, eight frequencies
    // instead of forty, and finer where it is read. A 2D panel solve is seconds
    // of work at full width and this suite runs eight times in the gate.
    const RadiationHull sections = radiationHullFromMesh(ship.hull, 4.0, 9);
    const RadiationTable table = stripTheoryTable(sections, radiationFrequencyGrid(0.55, 1.2, 8));

    const Diagnostics diag = ship.diagnostics(Sea(0.0));
    const double mass = 60.0 * 16.0 * 4.0 * kRhoSeawater;
    const double stiffness = mass * kGravity * diag.gmTransverse;
    const double dryInertia = mass * 5.0 * 5.0;   // gyradius kxx = 5 m
    expectTrue("the barge is stable enough to roll about", diag.gmTransverse > 1.0);

    const Sea still(0.0);
    ship.state.orientation = Quat::fromAxisAngle(Vec3{1, 0, 0}, 4.0 * kDegToRad);
    std::vector<double> crossings;
    double previous = 4.0 * kDegToRad;
    const double dt = 0.01;
    for (int i = 1; i <= 12000; ++i) {
        ship.step(dt, still);
        double heel = 0, trim = 0;
        heelTrimFromRotation(ship.state.orientation.toMat3(), heel, trim);
        if (previous < 0 && heel >= 0)
            crossings.push_back(i * dt - dt * heel / (heel - previous));
        previous = heel;
    }
    expectTrue("the roll decay produced enough cycles to time", crossings.size() >= 4);
    if (crossings.size() < 4) return;

    double total = 0;
    for (std::size_t i = 2; i < crossings.size(); ++i) total += crossings[i] - crossings[i - 1];
    const double omegaDamped =
        2.0 * kPi * static_cast<double>(crossings.size() - 2) / total;
    const double implied = stiffness / (omegaDamped * omegaDamped) - dryInertia;

    double aboutCog = 0, aboutOrigin = 0;
    for (int i = 0; i + 1 < table.size(); ++i) {
        const auto lo = static_cast<std::size_t>(i);
        if (table.omega[lo] > omegaDamped || omegaDamped > table.omega[lo + 1]) continue;
        const double f = (omegaDamped - table.omega[lo]) / (table.omega[lo + 1] - table.omega[lo]);
        aboutOrigin = table.addedMass[lo][3][3] * (1 - f) + table.addedMass[lo + 1][3][3] * f;
        aboutCog = transferAddedMass(table.addedMass[lo], diag.centreOfGravity)[3][3] * (1 - f) +
                   transferAddedMass(table.addedMass[lo + 1], diag.centreOfGravity)[3][3] * f;
    }
    expectTrue("the roll decay frequency lies inside the table", aboutCog > 0);

    expectNear("the decay sees A44(omega_d) referred to the centre of gravity", implied, aboutCog,
               0.15 * aboutCog);
    // The guard that makes the tolerance mean something.
    expectTrue("the two frames are far enough apart for this to discriminate",
               std::abs(aboutOrigin - aboutCog) > 0.5 * aboutCog);
    expectTrue("and the decay follows the centre-of-gravity value, not the baseline one",
               std::abs(implied - aboutCog) < std::abs(implied - aboutOrigin));
}

// --- Ikeda in the integrator -------------------------------------------------

// The hull form Ikeda is evaluated on is read off the mesh, so it must agree with
// the mesh. A box barge already checks the exact case in
// testStationsFromAMeshAreExactForABox; what is checked here is that a real hull
// form lands inside the domain the regressions were fitted over, because outside
// it the eddy polynomial returns whatever it returns.
// `waveDamping` must stay zero while a radiation model is attached, because the
// memory convolution already applies the radiation share of B44 -- and supplying
// both is the documented double count that cost 27% of mid-frequency heave when
// radiation first landed.
//
// **That rule was an imperative in a header and was enforced in exactly one
// place**: `attachRollDamping()` sets the field to zero. But `rollDampingForm` is
// a public optional and `waveDamping` a public field, so a caller who builds the
// form by hand, or assigns to it after attaching, re-creates the defect in
// silence -- the total flows straight into the roll damping coefficient.
//
// `validateRollDamping()` cannot catch it: it sees the form alone and has no way
// to know a `RadiationForce` exists, which is why it reports the *correct*
// configuration as a problem and every caller filters that string out. Only
// `Ship::validate()` can see both, and now does.
void testRollWaveDampingAndRadiationAreNotBothApplied() {
    Ship ship = rollShip();
    ship.attachRollDamping(kRollWaterline, kRollKeelLength, kRollKeelBreadth);
    ship.attachRadiation(kRollWaterline, 9);

    expectTrue("attachRollDamping leaves waveDamping at zero",
               ship.rollDampingForm && ship.rollDampingForm->waveDamping == 0.0);
    const auto clean = ship.validate();
    const bool cleanComplains =
        std::any_of(clean.begin(), clean.end(), [](const std::string& s) {
            return s.find("damped twice in roll") != std::string::npos;
        });
    expectTrue("so a correctly configured ship says nothing about it", !cleanComplains);

    // Now do what the header forbids, the way a caller actually could.
    ship.rollDampingForm->waveDamping = 5.0e7;
    const auto doubled = ship.validate();
    const bool caught =
        std::any_of(doubled.begin(), doubled.end(), [](const std::string& s) {
            return s.find("damped twice in roll") != std::string::npos;
        });
    expectTrue("but supplying both is reported", caught);

    // The guard: without radiation the same wave damping is legitimate -- it is
    // the 5-30% of B44 a caller with no radiation model is *supposed* to supply,
    // so the check must key on the pair and not on the field alone.
    Ship noRadiation = rollShip();
    noRadiation.attachRollDamping(kRollWaterline, kRollKeelLength, kRollKeelBreadth);
    noRadiation.rollDampingForm->waveDamping = 5.0e7;
    const auto lone = noRadiation.validate();
    const bool falsePositive =
        std::any_of(lone.begin(), lone.end(), [](const std::string& s) {
            return s.find("damped twice in roll") != std::string::npos;
        });
    expectTrue("while wave damping without radiation is not a problem", !falsePositive);
}

void testRollDampingFormComesOffTheHull() {
    Ship ship = rollShip();
    const RollDampingHull form =
        ship.attachRollDamping(kRollWaterline, kRollKeelLength, kRollKeelBreadth);

    expectNear("length is the wetted length of the mesh", form.lengthPp, 120.0, 1e-9);
    expectNear("beam is the widest the wetted body gets", form.beam, 22.0, 1e-9);
    expectNear("draft is the waterline down to the keel", form.draft, kRollWaterline, 1e-12);
    expectTrue("the block coefficient is a real hull's, not a box's",
               form.blockCoeff > 0.5 && form.blockCoeff < 0.85);
    expectTrue("and so is the midship coefficient",
               form.midshipCoeff > 0.9 && form.midshipCoeff < 0.99);
    // Cb <= Cm always: the midship section is the largest, so the prismatic
    // coefficient Cb/Cm cannot exceed 1.
    expectTrue("the block coefficient does not exceed the midship coefficient",
               form.blockCoeff <= form.midshipCoeff);

    // Displacement, derived two ways: from the coefficients the form reports, and
    // from the hydrostatic integral the ship actually floats on.
    const double fromForm =
        form.blockCoeff * form.lengthPp * form.beam * form.draft * kRhoSeawater;
    expectNear("the form's displacement is the hull's displacement", fromForm,
               ship.diagnostics(Sea(0.0)).displacementMass, 1e-6 * fromForm);

    // The roll axis is loading, not form, and comes from the live centre of
    // gravity rather than from anything in the mesh.
    expectNear("the roll axis is the centre of gravity above the keel", form.rollAxisAboveKeel,
               ship.diagnostics(Sea(0.0)).centreOfGravity.z, 1e-9);

    // Everything the method complains about, other than the wave component it
    // deliberately does not compute. A hull outside the regression domain would
    // make every coefficient below an extrapolation.
    RollDampingCondition at;
    at.rollAmplitude = 10.0 * kDegToRad;
    at.rollFrequency = 0.53;
    int unexpected = 0;
    for (const std::string& problem : validateRollDamping(form, at))
        if (problem.find("wave") == std::string::npos) ++unexpected;
    expectEqual("the test hull sits inside every published bound of the method",
                static_cast<long long>(unexpected), 0LL);

    // waveDamping is left at zero on purpose: with a radiation model attached the
    // memory convolution already applies it, and supplying both is the double
    // count that cost 27% of mid-frequency heave when radiation first landed.
    expectTrue("the wave component is left to the radiation model", form.waveDamping == 0.0);
}

// **The headline check.** Feed Ikeda's B44 to the real rigid-body integrator,
// release the ship at a heel, and require the logarithmic decrement of every
// cycle to be the 2 pi zeta / sqrt(1 - zeta^2) implied by the coefficient at that
// cycle's amplitude.
//
// Every ingredient comes from a different place, which is what stops this being a
// restatement of the integrator's own arithmetic:
//
//   * the restoring stiffness from GM, which diagnostics() gets from a forced-heel
//     GZ sweep and not from the finite difference the damper is scaled by;
//   * the effective roll inertia from the period of a decay with the damping
//     switched off entirely;
//   * B44 from roll_damping.cpp, called directly;
//   * the decrement from the simulation.
//
// The release is deliberately small. At 20 degrees this hull's restoring is
// visibly nonlinear -- the roll period runs 1.3% short -- and the closed form is
// a linear-restoring statement, so a large release would be testing GZ curvature
// under the name of damping.
void testRollFreeDecayMatchesTheIkedaCoefficient() {
    Ship undamped = rollShip();
    undamped.zetaRoll = 0.0;   // no modal stand-in, and no Ikeda attached
    const RollDecay natural = rollDecay(undamped, 3.0, 180.0);
    expectTrue("the undamped decay produced enough peaks to time", natural.peaks.size() >= 8);
    if (natural.peaks.size() < 8) return;

    // Vacuity guard, and the control for everything below: a ship with no roll
    // damping of any kind must not decay. Measured drift is +0.2% over fourteen
    // cycles, which is the ceiling on what the integrator itself contributes.
    expectTrue("with no roll damping at all the ship does not decay",
               std::abs(natural.peaks.back() / natural.peaks.front() - 1.0) < 0.02);

    const Diagnostics diag = undamped.diagnostics(Sea(0.0));
    const double stiffness = diag.displacementMass * kGravity * diag.gmTransverse;
    const double omegaNatural = 2.0 * kPi / natural.period();
    const double inertia = stiffness / (omegaNatural * omegaNatural);
    expectTrue("the stiffness and inertia are ship-sized", stiffness > 0 && inertia > 0);

    Ship ship = rollShip();
    const RollDampingHull form =
        ship.attachRollDamping(kRollWaterline, kRollKeelLength, kRollKeelBreadth);
    const RollDecay decay = rollDecay(ship, 3.0, 180.0);
    expectTrue("the damped decay produced enough peaks", decay.peaks.size() >= 8);
    if (decay.peaks.size() < 8) return;

    // The integrator picks its own operating point. Roll being sharply resonant
    // is the argument for using the natural frequency, and this is the check on
    // it: the frequency Ikeda was asked at has to be the frequency the ship
    // actually rolled at.
    expectNear("Ikeda was evaluated at the frequency the ship rolls at",
               ship.rollCondition.rollFrequency, omegaNatural, 0.01 * omegaNatural);
    expectTrue("and at an amplitude that followed the decay down",
               ship.rollCondition.rollAmplitude > 0 &&
                   ship.rollCondition.rollAmplitude < 0.5 * decay.peaks.front());

    double smallest = 1e30, largest = 0;
    for (std::size_t i = 1; i < decay.peaks.size(); ++i) {
        RollDampingCondition at;
        // The amplitude over the cycle, which for a slowly decaying envelope is
        // the mean of its endpoints. B44 is affine in amplitude, so the mean of
        // the amplitude is the amplitude of the mean.
        at.rollAmplitude = 0.5 * (decay.peaks[i - 1] + decay.peaks[i]);
        at.rollFrequency = omegaNatural;
        const double zeta =
            rollDamping(form, at).total / (2.0 * std::sqrt(inertia * stiffness));
        const double wanted = 2.0 * kPi * zeta / std::sqrt(1.0 - zeta * zeta);
        smallest = std::min(smallest, wanted);
        largest = std::max(largest, wanted);
        expectTrue("the damping ratio is light, as roll always is", zeta > 0.005 && zeta < 0.2);
        expectNear("logarithmic decrement matches 2 pi zeta / sqrt(1 - zeta^2)",
                   decay.decrement(i), wanted, 0.02 * wanted);
    }

    // Guard against the whole comparison being a constant against a constant. If
    // B44 did not move with amplitude, a linear damper would pass this test and
    // the amplitude handling would be untested.
    expectTrue("the predicted decrement genuinely varies across the record",
               largest > 1.2 * smallest);
}

// Bilge keels are the strongest check available on the coupling because the
// statement is monotone, physical and untunable: adding them can only take energy
// out faster. Nothing about their size or the fitted constants is asserted.
void testBilgeKeelsShortenTheRollDecay() {
    Ship keeled = rollShip();
    const RollDampingHull withKeels =
        keeled.attachRollDamping(kRollWaterline, kRollKeelLength, kRollKeelBreadth);
    Ship bare = rollShip();
    const RollDampingHull withoutKeels = bare.attachRollDamping(kRollWaterline);
    Ship none = rollShip();
    none.zetaRoll = 0.0;

    // At the coefficient level first, so a failure downstream can be told apart
    // from a failure here.
    RollDampingCondition at;
    at.rollAmplitude = 5.0 * kDegToRad;
    at.rollFrequency = 0.53;
    const RollDamping keeledB = rollDamping(withKeels, at);
    const RollDamping bareB = rollDamping(withoutKeels, at);
    expectTrue("bilge keels add damping", keeledB.total > bareB.total);
    expectTrue("and they are the dominant contribution once fitted",
               keeledB.bilgeKeel() > 0.5 * keeledB.total);
    expectTrue("a hull without them reports exactly none", bareB.bilgeKeel() == 0.0);

    const RollDecay withDecay = rollDecay(keeled, 8.0, 180.0);
    const RollDecay withoutDecay = rollDecay(bare, 8.0, 180.0);
    const RollDecay noneDecay = rollDecay(none, 8.0, 180.0);
    const std::size_t cycles =
        std::min({withDecay.peaks.size(), withoutDecay.peaks.size(), noneDecay.peaks.size()});
    expectTrue("all three decays ran long enough to compare", cycles >= 8);
    if (cycles < 8) return;

    // **The damping-ratio table `docs/02-simulation.md` publishes, printed.**
    // It was a hand-made table with no producer anywhere in the tree: grepping
    // its own figures -- 0.0259, 0.0323, 0.0470, 0.0810, 0.00062, 0.00193 --
    // across `tests/`, `tools/` and `engine/` returned nothing for any of them, so
    // re-deriving it meant writing a driver against the library. That is the exact
    // condition `check-figures.sh:676-682` describes about the GM-detail figures
    // before `--gm-detail` existed, and the answer that worked there was to teach
    // the binary to print what the document claims. Everything this needs was
    // already sitting in this function.
    //
    // `none` carries no damping of any kind, so its decay period is the natural
    // one rather than a damped one, and the critical coefficient follows from the
    // stiffness the same way `testRollFreeDecayMatchesTheIkedaCoefficient` takes
    // it: `2 sqrt(C I)` with `I = C / omega_n^2`, which is `2 C / omega_n`.
    const double omegaNatural = 2.0 * kPi / noneDecay.period();
    const Diagnostics rollDiag = none.diagnostics(Sea(0.0));
    const double rollStiffness = rollDiag.displacementMass * kGravity * rollDiag.gmTransverse;
    const double criticalRoll = 2.0 * rollStiffness / omegaNatural;
    expectTrue("the critical roll coefficient is ship-sized", criticalRoll > 0);
    std::printf("     roll damping ratio at omega_n %.4f rad/s, OG/d %.3f\n", omegaNatural,
                (withKeels.draft - withKeels.rollAxisAboveKeel) / withKeels.draft);
    std::printf("     %9s %12s %12s %12s\n", "amplitude", "zeta keels", "zeta bare", "keel share");
    for (double deg : {2.5, 5.0, 10.0, 20.0}) {
        RollDampingCondition sample;
        sample.rollAmplitude = deg * kDegToRad;
        sample.rollFrequency = omegaNatural;
        const RollDamping k = rollDamping(withKeels, sample);
        const RollDamping b = rollDamping(withoutKeels, sample);
        std::printf("     %8.1f %12.4f %12.5f %11.1f%%\n", deg, k.total / criticalRoll,
                    b.total / criticalRoll, 100.0 * k.bilgeKeel() / k.total);
    }

    const double kept = withDecay.peaks[cycles - 1] / withDecay.peaks.front();
    const double bareKept = withoutDecay.peaks[cycles - 1] / withoutDecay.peaks.front();
    const double noneKept = noneDecay.peaks[cycles - 1] / noneDecay.peaks.front();

    // Vacuity guard: the bare hull must itself be damped by Ikeda rather than by
    // the integrator, or "keels beat bare" would be a comparison against noise.
    expectTrue("the bare hull's decay is Ikeda's and not the integrator's",
               bareKept < 0.9 * noneKept);
    expectTrue("bilge keels shorten the decay a great deal", kept < 0.25 * bareKept);
    expectTrue("and every cycle of the keeled decay is shorter than the bare one's",
               withDecay.decrement(1) > withoutDecay.decrement(1) &&
                   withDecay.decrement(cycles - 1) > withoutDecay.decrement(cycles - 1));
}

// Roll damping is strongly amplitude-dependent, which is the one property that
// cannot survive being replaced by any constant coefficient. A large decay must
// therefore show its decrement *falling* as the amplitude falls -- and the same
// ship run on the fraction-of-critical stand-in, through the same integrator and
// the same peak finder, must hold it constant.
//
// That pairing is what makes this a test of the physics rather than of the
// measurement: both runs share every source of error except the damper.
void testLargeRollDecayHasAnAmplitudeDependentDecrement() {
    Ship ikeda = rollShip();
    ikeda.attachRollDamping(kRollWaterline, kRollKeelLength, kRollKeelBreadth);
    Ship linear = rollShip();   // zetaRoll = 0.08, the lumped stand-in

    const RollDecay nonlinearDecay = rollDecay(ikeda, 20.0, 180.0);
    const RollDecay linearDecay = rollDecay(linear, 20.0, 180.0);
    expectTrue("both decays produced enough cycles",
               nonlinearDecay.peaks.size() >= 8 && linearDecay.peaks.size() >= 8);
    if (nonlinearDecay.peaks.size() < 8 || linearDecay.peaks.size() < 8) return;

    // The release really was large, and really did decay a long way -- otherwise
    // "the decrement changed" would be a statement about a narrow amplitude band.
    expectTrue("the decay spans a wide range of amplitudes",
               nonlinearDecay.peaks.front() > 8.0 * nonlinearDecay.peaks.back());

    const std::size_t last = nonlinearDecay.peaks.size() - 1;
    expectTrue("the decrement falls as the roll decays",
               nonlinearDecay.decrement(1) > 1.5 * nonlinearDecay.decrement(last));
    for (std::size_t i = 2; i <= last; ++i)
        expectTrue("and falls monotonically, as B44 does with amplitude",
                   nonlinearDecay.decrement(i) < nonlinearDecay.decrement(i - 1));

    // The control. A constant fraction of critical gives a constant decrement, so
    // whatever the nonlinear run is showing is the damper and not the hull, the
    // peak finder or the restoring curve.
    const std::size_t lastLinear = linearDecay.peaks.size() - 1;
    const double spread =
        std::abs(linearDecay.decrement(1) / linearDecay.decrement(lastLinear) - 1.0);
    expectTrue("a linear damper holds its decrement constant through the same decay",
               spread < 0.03);

    // And the stand-in means what its name says: zetaRoll = 0.08 delivers 0.08 of
    // critical. That is a statement about the *stiffness* the damper is scaled by,
    // and it did not hold until the finite difference was made to rotate the ship
    // about her centre of gravity rather than about the body origin.
    const double zeta = 0.08;
    expectNear("zetaRoll delivers the fraction of critical it claims",
               linearDecay.decrement(lastLinear), 2.0 * kPi * zeta / std::sqrt(1.0 - zeta * zeta),
               0.02);
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
    testStationsFromAMeshAreExactForABox();
    testRadiationMemoryDampsRatherThanDiverges();
    testFreeDecayAgreesWithTheFrequencyDomainTable();
    testPropulsionDrivesTheShip();
    testRudderTurnsTowardsTheCommandedSide();
    testEncounterFrequencyEmergesFromMovingThroughTheWave();
    testRaoMeasuresTheSpeedItActuallyReached();

    // --- viscous roll damping and the roll added-inertia frame ---
    std::printf("\n--- roll damping and added inertia coupled into the ship ---\n");
    testAddedMassTransferObeysTheParallelAxisTheorem();
    testSemicircularSectionsPinTheTransferSign();
    testRollAddedInertiaIsReferredToTheCentreOfGravity();
    testRollDecaySeesAddedInertiaAboutTheCentreOfGravity();
    testRollWaveDampingAndRadiationAreNotBothApplied();
    testRollDampingFormComesOffTheHull();
    testRollFreeDecayMatchesTheIkedaCoefficient();
    testBilgeKeelsShortenTheRollDecay();
    testLargeRollDecayHasAnAmplitudeDependentDecrement();
}
