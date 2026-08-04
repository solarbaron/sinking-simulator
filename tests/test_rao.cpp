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

#include <cmath>
#include <cstdio>
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
}
