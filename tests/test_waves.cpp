// SPDX-License-Identifier: MIT
//
// Validation of the spectral wave field.
//
// A wave field is the easiest thing in a simulator to get plausibly wrong: any
// sum of cosines looks like the sea. So nothing here is eyeballed. Every check is
// against a closed form worked out beforehand -- the analytic integral of the
// spectrum, the deep-water dispersion relation, textbook period ratios for
// Pierson-Moskowitz, the exact Airy wave, the exp(k z) decay -- and the checks
// are aimed at the failures that stay silent:
//
//   * a spreading function that does not normalise (the sea is then the wrong
//     height, and only in a directional sea, which is every real one);
//   * g on the wrong side of the dispersion relation (waves of the right shape
//     travelling at the wrong speed);
//   * frequency spacing that loses energy in the tail (a sea that is too smooth
//     exactly at the scale of a hatch coaming).
#include "engine/sim/waves.hpp"
#include "harness.hpp"

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

using namespace sim;
using testing::expectEqual;
using testing::expectNear;
using testing::expectTrue;

namespace {

// Independent composite Simpson. Deliberately not the engine's quadrature: a
// test that reuses the code under test's integrator cannot catch a bad
// integrator.
template <class F>
double integrate(F f, double lo, double hi, int panels) {
    const double h = (hi - lo) / panels;
    double sum = f(lo) + f(hi);
    for (int i = 1; i < panels; ++i) sum += (i & 1 ? 4.0 : 2.0) * f(lo + h * i);
    return sum * h / 3.0;
}

std::string label(const char* what, double value) {
    char buffer[128];
    std::snprintf(buffer, sizeof buffer, "%s %.4g", what, value);
    return buffer;
}

// --- The analytic spectrum ---------------------------------------------------

// The single fact the whole model rests on: whatever gamma is, the spectrum must
// integrate to Hs^2/16, because that is the definition of Hs.
void testSpectrumIntegratesToTheRequestedHs() {
    for (double gamma : {1.0, 2.0, 3.3, 5.0}) {
        SeaState sea;
        sea.significantHeight = 3.5;
        sea.peakPeriod = 11.0;
        sea.peakEnhancement = gamma;
        sea.shape = gamma == 1.0 ? SpectrumShape::PiersonMoskowitz : SpectrumShape::Jonswap;
        const Spectrum spectrum(sea);
        const auto density = [&](double omega) { return spectrum.density(omega); };
        // Fine over the energetic band, coarse over the omega^-5 tail. Above
        // 500 rad/s the closed-form tail integral leaves 1.4e-12 of the total,
        // which is why the tolerance below is 1e-9 and not tighter.
        const double m0 =
            integrate(density, 1e-3, 5.0, 200000) + integrate(density, 5.0, 500.0, 20000);
        expectNear(label("integral of S(omega) is Hs^2/16 at gamma", gamma), m0,
                   3.5 * 3.5 / 16.0, 1e-9);
    }
}

// The peak-enhancement normalisation. Pierson-Moskowitz needs none, and the
// JONSWAP value must agree with the DNV closed-form approximation
// 1 - 0.287 ln(gamma), which is derived independently of anything here.
void testPeakEnhancementNormalisation() {
    SeaState pm;
    pm.shape = SpectrumShape::PiersonMoskowitz;
    expectNear("Pierson-Moskowitz needs no normalisation", Spectrum(pm).normalisation(), 1.0,
               1e-12);

    for (double gamma : {1.5, 3.3, 6.0}) {
        SeaState sea;
        sea.peakEnhancement = gamma;
        expectNear(label("normalisation matches DNV 1 - 0.287 ln gamma at", gamma),
                   Spectrum(sea).normalisation(), 1.0 - 0.287 * std::log(gamma), 0.01);
    }
}

// The JONSWAP peak enhancement is deliberately asymmetric -- sigma is 0.07 below
// the peak and 0.09 above it, so the spectrum is sharper on the long-wave side.
// Swapping the two branches leaves Hs, Tp and the total energy untouched and only
// changes the shape, which is exactly the kind of error that survives every other
// check here. So it gets its own closed form: dividing S by the bare
// Pierson-Moskowitz shape, computed independently below, must leave
// gamma^exp(-(x-1)^2 / (2 sigma^2)) and nothing else.
void testPeakEnhancementShapeIsAsymmetric() {
    SeaState sea;
    sea.significantHeight = 3.0;
    sea.peakPeriod = 10.0;
    sea.peakEnhancement = 3.3;
    const Spectrum spectrum(sea);
    const double omegaPeak = spectrum.peakFrequency();
    const double scale =
        5.0 / 16.0 * sea.significantHeight * sea.significantHeight * spectrum.normalisation() /
        omegaPeak;
    const auto bare = [](double x) {
        return std::pow(x, -5.0) * std::exp(-1.25 * std::pow(x, -4.0));
    };

    for (double delta : {0.03, 0.07, 0.15, 0.3}) {
        for (int side = -1; side <= 1; side += 2) {
            const double x = 1.0 + side * delta;
            const double sigma = side < 0 ? 0.07 : 0.09;
            const double want = std::pow(3.3, std::exp(-0.5 * delta * delta / (sigma * sigma)));
            const double got = spectrum.density(x * omegaPeak) / (scale * bare(x));
            expectNear(label(side < 0 ? "enhancement below the peak at delta"
                                      : "enhancement above the peak at delta",
                             delta),
                       got, want, 1e-12 * want);
        }
    }
    expectTrue("the enhancement is broader above the peak than below",
               spectrum.density(1.1 * omegaPeak) / bare(1.1) >
                   spectrum.density(0.9 * omegaPeak) / bare(0.9));
}

// S must peak at omega_p, for JONSWAP as well as Pierson-Moskowitz: both factors
// of the JONSWAP product are maximised there, so the product is too.
void testSpectrumPeaksAtTheRequestedPeriod() {
    for (double gamma : {1.0, 3.3}) {
        SeaState sea;
        sea.peakPeriod = 12.0;
        sea.peakEnhancement = gamma;
        sea.shape = gamma == 1.0 ? SpectrumShape::PiersonMoskowitz : SpectrumShape::Jonswap;
        const Spectrum spectrum(sea);
        const double omegaPeak = 2.0 * kPi / 12.0;
        expectNear(label("peak frequency is 2 pi / Tp at gamma", gamma), spectrum.peakFrequency(),
                   omegaPeak, 1e-15);
        const double atPeak = spectrum.density(omegaPeak);
        bool dominates = true;
        for (int i = 1; i <= 400; ++i) {
            const double offset = 0.005 * i;
            dominates = dominates && spectrum.density(omegaPeak - offset) < atPeak &&
                        spectrum.density(omegaPeak + offset) < atPeak;
        }
        expectTrue(label("S is maximised at omega_p at gamma", gamma), dominates);
    }
}

// --- Dispersion --------------------------------------------------------------

void testDeepWaterDispersion() {
    SeaState sea;
    sea.frequencyCount = 24;
    sea.directionCount = 5;
    const WaveField field(sea);
    expectEqual("component count is N x M",
                static_cast<long long>(field.components().size()), 24 * 5);

    bool dispersion = true, phase = true, group = true;
    for (const WaveComponent& c : field.components()) {
        dispersion = dispersion && std::abs(c.omega * c.omega - kGravity * c.wavenumber) <
                                       1e-12 * c.omega * c.omega;
        const double speed = c.omega / c.wavenumber;
        phase = phase && std::abs(speed - phaseSpeed(c.omega)) < 1e-12 * speed;
        group = group && std::abs(groupSpeed(c.omega) - 0.5 * speed) < 1e-12 * speed;
    }
    expectTrue("every component satisfies omega^2 = g k", dispersion);
    expectTrue("phase speed is omega / k", phase);
    expectTrue("group speed is half the phase speed in deep water", group);

    // And the helpers themselves, against algebra rather than against each other.
    expectNear("wavenumber of a 10 s wave", deepWaterWavenumber(2.0 * kPi / 10.0),
               (2.0 * kPi / 10.0) * (2.0 * kPi / 10.0) / kGravity, 1e-15);
    expectNear("a 10 s deep-water wave is 156.13 m long",
               2.0 * kPi / deepWaterWavenumber(2.0 * kPi / 10.0),
               kGravity * 100.0 / (2.0 * kPi), 1e-9);
}

// --- The discretisation ------------------------------------------------------

// The single most important check in this file. Equal-energy bins that run to
// zero and infinity lose nothing, so the round trip is exact, not approximate.
void testSignificantHeightRoundTrip() {
    for (double hs : {0.5, 3.0, 14.0}) {
        for (int frequencies : {1, 8, 48, 128}) {
            SeaState sea;
            sea.significantHeight = hs;
            sea.frequencyCount = frequencies;
            sea.directionCount = 7;
            const WaveField field(sea);
            expectNear(label("m0 of the discretisation is Hs^2/16 at Hs", hs) + " x " +
                           std::to_string(frequencies),
                       field.zerothMoment(), hs * hs / 16.0, 1e-12 * hs * hs);
            expectNear(label("4 sqrt(m0) recovers Hs at Hs", hs) + " x " +
                           std::to_string(frequencies),
                       field.significantHeight(), hs, 1e-12 * hs);
        }
    }
}

// The discretised m0 must also match an independent numerical integral of the
// analytic S(omega). This is the check that a truncated band would fail: it
// compares the sum of the components against the continuous spectrum they are
// supposed to represent, not against the construction that produced them.
void testDiscretisationMatchesTheAnalyticIntegral() {
    SeaState sea;
    sea.significantHeight = 4.2;
    sea.peakPeriod = 10.5;
    sea.frequencyCount = 64;
    sea.directionCount = 9;
    const WaveField field(sea);
    const Spectrum spectrum(sea);
    const auto density = [&](double omega) { return spectrum.density(omega); };
    const double analytic =
        integrate(density, 1e-3, 5.0, 200000) + integrate(density, 5.0, 500.0, 20000);
    expectNear("discrete m0 equals the numerical integral of S", field.zerothMoment(), analytic,
               1e-8);
}

// Equal-energy bins have a closed-form inverse for Pierson-Moskowitz, because
// its cumulative energy distribution is exactly exp(-1.25 (omega_p/omega)^4).
// Checking the bin edges against it validates the whole quantile machinery,
// including the numerically tabulated part that the JONSWAP case relies on.
void testFrequencyBinEdgesAgainstTheClosedForm() {
    SeaState sea;
    sea.shape = SpectrumShape::PiersonMoskowitz;
    sea.peakPeriod = 9.0;
    sea.frequencyCount = 64;
    sea.directionCount = 1;
    const WaveField field(sea);
    const double omegaPeak = 2.0 * kPi / 9.0;
    const auto& bins = field.frequencyBins();
    expectEqual("one bin per frequency", static_cast<long long>(bins.size()), 64);
    expectTrue("the first bin reaches down to zero frequency", bins.front().omegaLow == 0.0);
    expectTrue("the last bin reaches up to infinity", std::isinf(bins.back().omegaHigh));

    double worst = 0.0;
    for (int i = 1; i < 64; ++i) {
        const double want = omegaPeak * std::pow(1.25 / std::log(64.0 / i), 0.25);
        worst = std::max(worst, std::abs(bins[static_cast<std::size_t>(i)].omegaLow - want) / want);
    }
    expectTrue("Pierson-Moskowitz bin edges match exp(-1.25 (wp/w)^4) to 1e-12", worst < 1e-12);

    bool contiguous = true, equalEnergy = true;
    for (int i = 0; i < 64; ++i) {
        const FrequencyBin& bin = bins[static_cast<std::size_t>(i)];
        contiguous = contiguous && bin.omega > bin.omegaLow && bin.omega < bin.omegaHigh;
        equalEnergy =
            equalEnergy && std::abs(bin.energy - bins[0].energy) < 1e-15 * bins[0].energy;
    }
    expectTrue("every bin's centroid lies inside the bin", contiguous);
    expectTrue("every bin carries the same energy", equalEnergy);
}

// Each bin must actually contain the energy it claims: the analytic integral of
// S over the bin's own frequency range, computed here by the test's quadrature,
// must equal the energy the component was given. This is the assertion aimed
// squarely at a discretisation that loses energy in the tail -- it compares the
// components against the continuous spectrum interval by interval, so a band that
// stops short, or bins whose widths do not tile the spectrum, shows up as a
// mismatch in the bin that was shortchanged rather than as a small error in a
// total that other bins could mask.
void testEachBinHoldsTheEnergyItClaims() {
    SeaState sea;
    sea.significantHeight = 3.0;
    sea.peakPeriod = 9.0;
    sea.peakEnhancement = 3.3;
    sea.frequencyCount = 32;
    sea.directionCount = 1;
    const WaveField field(sea);
    const Spectrum spectrum(sea);
    const auto density = [&](double omega) { return spectrum.density(omega); };

    const auto& bins = field.frequencyBins();
    double worst = 0.0, total = 0.0;
    for (std::size_t i = 0; i + 1 < bins.size(); ++i) {
        const double got = integrate(density, std::max(bins[i].omegaLow, 1e-3), bins[i].omegaHigh,
                                     4000);
        total += got;
        worst = std::max(worst, std::abs(got - bins[i].energy) / bins[i].energy);
    }
    expectTrue("every closed bin holds exactly the energy its component carries", worst < 1e-9);
    // The open-topped bin, integrated out to 2000 rad/s: a 6 mm wave at 3 ms
    // period. Whatever is above that is 5e-15 of m0 by the closed-form tail
    // integral, so the whole spectrum is accounted for.
    const double tail = integrate(density, bins.back().omegaLow, 2000.0, 400000);
    total += tail;
    expectNear("the open-topped bin holds the whole remaining tail", tail, bins.back().energy,
               1e-9 * bins.back().energy);
    expectNear("the bins tile the spectrum with nothing left over", total,
               field.zerothMoment(), 1e-9);
}

// The discretised spectrum -- bin energy over bin width -- must peak at the
// requested Tp. With equal-energy bins the narrowest bin is the most energetic
// one, which is a different route to the same answer than the analytic peak.
void testDiscretisedSpectrumPeaksAtTp() {
    for (double gamma : {1.0, 3.3}) {
        SeaState sea;
        sea.peakPeriod = 13.0;
        sea.peakEnhancement = gamma;
        sea.shape = gamma == 1.0 ? SpectrumShape::PiersonMoskowitz : SpectrumShape::Jonswap;
        sea.frequencyCount = 64;
        sea.directionCount = 1;
        const WaveField field(sea);
        const double omegaPeak = 2.0 * kPi / 13.0;

        const auto& bins = field.frequencyBins();
        std::size_t best = 0;
        double bestDensity = -1.0;
        for (std::size_t i = 0; i < bins.size(); ++i) {
            const double width = bins[i].omegaHigh - bins[i].omegaLow;
            if (!std::isfinite(width) || width <= 0.0) continue;
            const double d = bins[i].energy / width;
            if (d > bestDensity) { bestDensity = d; best = i; }
        }
        expectTrue(label("the densest bin brackets omega_p at gamma", gamma),
                   bins[best].omegaLow <= omegaPeak && omegaPeak <= bins[best].omegaHigh);
        expectNear(label("the densest bin's period is Tp at gamma", gamma),
                   2.0 * kPi / bins[best].omega, 13.0, 0.13);
    }
}

// Textbook Pierson-Moskowitz period ratios: T1 = 0.7718 Tp and T2 = 0.7104 Tp,
// both derived from Gamma functions and neither used anywhere in the engine.
//
// T1 is preserved exactly by the discretisation, because each bin sits at its
// energy *centroid* and centroids sum to the first moment. T2 is not: a single
// frequency cannot represent the spread of frequencies inside its own bin, so
// the discrete m2 is biased low and T2 comes out slightly long. The bias is
// asserted with its own tolerance rather than hidden behind a loose one on both.
void testMeanPeriodsAgainstGammaFunctions() {
    SeaState sea;
    sea.shape = SpectrumShape::PiersonMoskowitz;
    sea.peakPeriod = 11.0;
    sea.frequencyCount = 96;
    sea.directionCount = 1;
    const WaveField field(sea);

    // m1/m0 = (1/4) 1.25^(-3/4) Gamma(3/4) / (1/5), in units of omega_p.
    const double meanRatio =
        1.0 / (0.25 * std::pow(1.25, -0.75) * std::tgamma(0.75) / 0.2);
    const double crossRatio =
        1.0 / std::sqrt(0.25 * std::pow(1.25, -0.5) * std::tgamma(0.5) / 0.2);
    // T1 is exact, not approximate, so the tolerance is 1e-8 and not a percent.
    // A discretisation that placed each component at its bin's midpoint instead
    // of its centroid would still look right to within a percent, and this is the
    // assertion that says it is not right.
    expectNear("PM mean period T1 is 0.7718 Tp", field.meanPeriod(), meanRatio * 11.0,
               1e-8 * 11.0);
    // The bias is +0.66% at N = 96 and halves as N doubles; 1.5% is the bound.
    expectNear("PM zero-crossing period T2 is 0.7104 Tp", field.zeroCrossingPeriod(),
               crossRatio * 11.0, 0.015 * 11.0);
    expectTrue("T2 is biased long, not short, by the centroid discretisation",
               field.zeroCrossingPeriod() >= crossRatio * 11.0);
}

// --- Directional spreading ---------------------------------------------------

void testSpreadingNormalisation() {
    const double mean = 0.7;
    for (double s : {0.0, 0.5, 1.0, 2.0, 4.0, 10.0}) {
        const auto d = [&](double theta) { return directionalSpreading(theta, mean, s); };
        expectNear(label("cos^2s spreading integrates to 1 at s", s),
                   integrate(d, mean - kPi, mean + kPi, 200000), 1.0, 1e-9);
    }

    // Closed forms that do not go through lgamma: s = 1 gives (1 + cos d)/(2 pi)
    // and s = 2 gives a peak of 4/(3 pi).
    expectNear("s = 1 peak density is 1/pi", directionalSpreading(mean, mean, 1.0), 1.0 / kPi,
               1e-14);
    expectNear("s = 1 density at 90 degrees is 1/(2 pi)",
               directionalSpreading(mean + 0.5 * kPi, mean, 1.0), 0.5 / kPi, 1e-14);
    expectNear("s = 2 peak density is 4/(3 pi)", directionalSpreading(mean, mean, 2.0),
               4.0 / (3.0 * kPi), 1e-14);
    expectNear("s = 0 is isotropic", directionalSpreading(mean + 2.0, mean, 0.0), 0.5 / kPi,
               1e-14);
    expectNear("spreading is symmetric about the mean direction",
               directionalSpreading(mean + 1.1, mean, 3.0),
               directionalSpreading(mean - 1.1, mean, 3.0), 1e-15);
    expectNear("spreading wraps around the circle", directionalSpreading(mean + 1.1, mean, 3.0),
               directionalSpreading(mean + 1.1 + 4.0 * kPi, mean, 3.0), 1e-15);
}

// The first circular moment of cos^(2s)(d/2) spreading is exactly s/(s+1) --
// another closed form the engine never uses. Checked analytically, then on the
// discretised bins, which is what actually reaches the components.
void testSpreadingFirstMoment() {
    for (double s : {1.0, 4.0, 20.0}) {
        const auto d = [&](double delta) {
            return directionalSpreading(delta, 0.0, s) * std::cos(delta);
        };
        expectNear(label("analytic E[cos delta] is s/(s+1) at s", s),
                   integrate(d, -kPi, kPi, 200000), s / (s + 1.0), 1e-9);

        SeaState sea;
        sea.spreadingExponent = s;
        sea.directionCount = 64;
        sea.frequencyCount = 1;
        const WaveField field(sea);
        double sum = 0.0;
        for (const DirectionBin& bin : field.directionBins())
            sum += bin.energy * std::cos(bin.direction);
        expectNear(label("discrete E[cos delta] is s/(s+1) at s", s), sum, s / (s + 1.0), 2e-3);
    }
}

void testDirectionBinsPartitionTheEnergy() {
    SeaState sea;
    sea.spreadingExponent = 6.0;
    sea.directionCount = 16;
    sea.frequencyCount = 4;
    sea.meanDirection = 1.3;
    const WaveField field(sea);
    const auto& bins = field.directionBins();
    double energy = 0.0;
    bool ordered = true;
    for (std::size_t i = 0; i < bins.size(); ++i) {
        energy += bins[i].energy;
        ordered = ordered && bins[i].directionLow < bins[i].direction &&
                  bins[i].direction < bins[i].directionHigh;
        if (i > 0) ordered = ordered && bins[i - 1].directionHigh == bins[i].directionLow;
    }
    expectNear("direction bins carry the whole energy", energy, 1.0, 1e-14);
    expectTrue("direction bins are contiguous and centroid-ordered", ordered);
    expectNear("the spread is symmetric about the mean direction",
               0.5 * (bins.front().direction + bins.back().direction), 1.3, 1e-15);
}

// A very large s must be long-crested: elevation varies along the direction of
// travel and, to tolerance, not across it.
void testLargeSpreadingExponentIsLongCrested() {
    SeaState sea;
    sea.significantHeight = 3.0;
    sea.peakPeriod = 9.0;
    sea.spreadingExponent = 5.0e5;
    sea.directionCount = 8;
    sea.frequencyCount = 32;
    sea.meanDirection = 0.0;  // travelling along +x
    const WaveField field(sea);

    double acrossSpread = 0.0, alongSpread = 0.0;
    const double reference = field.elevation(0.0, 0.0, 0.0);
    double alongMin = reference, alongMax = reference;
    for (int i = -40; i <= 40; ++i) {
        acrossSpread = std::max(acrossSpread, std::abs(field.elevation(0.0, i * 2.5, 0.0) - reference));
        const double along = field.elevation(i * 2.5, 0.0, 0.0);
        alongMin = std::min(alongMin, along);
        alongMax = std::max(alongMax, along);
    }
    alongSpread = alongMax - alongMin;
    expectTrue("a long-crested sea varies along the direction of travel", alongSpread > 1.0);
    expectTrue("a long-crested sea is flat across it over 200 m", acrossSpread < 0.02);
    expectTrue("the crest-wise variation is two orders of magnitude smaller",
               acrossSpread * 100.0 < alongSpread);
}

// --- The exact Airy wave -----------------------------------------------------

// One frequency, one direction: the field must be exactly a cos(k x - omega t +
// phi), with the kinematics that go with it. This is the test that catches the
// sign and convention errors a spectrum would average away.
void testSingleComponentIsAnExactAiryWave() {
    SeaState sea;
    sea.significantHeight = 2.0;
    sea.peakPeriod = 8.0;
    sea.frequencyCount = 1;
    sea.directionCount = 1;
    sea.meanDirection = 0.0;
    const WaveField field(sea);
    expectEqual("N = 1, M = 1 gives one component",
                static_cast<long long>(field.components().size()), 1);

    const WaveComponent c = field.components()[0];
    expectNear("the single amplitude carries all the variance", c.amplitude,
               std::sqrt(2.0) * 2.0 / 4.0, 1e-15);
    expectNear("the single direction is exactly the mean direction", c.direction, 0.0, 0.0);
    expectTrue("the phase is inside [0, 2 pi)", c.phase >= 0.0 && c.phase < 2.0 * kPi);

    const double k = c.wavenumber, w = c.omega, a = c.amplitude, p = c.phase;
    bool exact = true, velocityExact = true, accelerationExact = true;
    for (double t : {0.0, 0.37, 4.2, -1.9}) {
        for (double x : {0.0, 3.5, -21.0, 140.0}) {
            const double psi = k * x - w * t + p;
            exact = exact && std::abs(field.elevation(x, 0.0, t) - a * std::cos(psi)) < 1e-14;
            const WaveKinematics kin = field.kinematics({x, 0.0, 0.0}, t);
            velocityExact = velocityExact &&
                            std::abs(kin.velocity.x - a * w * std::cos(psi)) < 1e-14 &&
                            std::abs(kin.velocity.y) < 1e-300 &&
                            std::abs(kin.velocity.z - a * w * std::sin(psi)) < 1e-14;
            accelerationExact = accelerationExact &&
                                std::abs(kin.acceleration.x - a * w * w * std::sin(psi)) < 1e-13 &&
                                std::abs(kin.acceleration.z + a * w * w * std::cos(psi)) < 1e-13;
        }
    }
    expectTrue("elevation is exactly a cos(k x - omega t + phi)", exact);
    expectTrue("surface orbital velocity is (a w cos psi, 0, a w sin psi)", velocityExact);
    expectTrue("surface orbital acceleration is (a w^2 sin psi, 0, -a w^2 cos psi)",
               accelerationExact);

    // The wave must travel in +x at the phase speed: a crest followed for t
    // seconds is still a crest.
    const double speed = phaseSpeed(w);
    bool travels = true;
    for (double t : {0.5, 3.0, 17.0})
        travels = travels && std::abs(field.elevation(speed * t, 0.0, t) -
                                      field.elevation(0.0, 0.0, 0.0)) < 1e-12;
    expectTrue("the wave travels along +x at the phase speed", travels);

    // w = d eta / dt at the surface, checked against a central difference of the
    // elevation rather than against the same algebra.
    const double h = 1e-5;
    const double dEta = (field.elevation(7.0, 0.0, 1.0 + h) - field.elevation(7.0, 0.0, 1.0 - h)) /
                        (2.0 * h);
    expectNear("vertical orbital velocity is d eta / dt at z = 0",
               field.velocity({7.0, 0.0, 0.0}, 1.0).z, dEta, 1e-6);
    const double dU = (field.velocity({7.0, 0.0, -3.0}, 1.0 + h).x -
                       field.velocity({7.0, 0.0, -3.0}, 1.0 - h).x) /
                      (2.0 * h);
    expectNear("horizontal acceleration is du/dt", field.acceleration({7.0, 0.0, -3.0}, 1.0).x, dU,
               1e-6);

    // At z = 0 a single component satisfies u = omega eta exactly.
    expectNear("u equals omega times eta for one component at the surface",
               field.velocity({11.0, 0.0, 0.0}, 2.5).x, w * field.elevation(11.0, 0.0, 2.5),
               1e-14);
}

// A wave travelling along +y must vary with y and not with x -- the other half
// of the convention check.
void testDirectionRotatesTheWave() {
    SeaState sea;
    sea.frequencyCount = 1;
    sea.directionCount = 1;
    sea.meanDirection = 0.5 * kPi;
    const WaveField field(sea);
    const WaveComponent c = field.components()[0];
    expectNear("dirX is cos(direction)", c.dirX, std::cos(0.5 * kPi), 1e-15);
    expectNear("dirY is sin(direction)", c.dirY, 1.0, 1e-15);

    bool exact = true;
    for (double y : {0.0, 12.0, -55.0}) {
        const double psi = c.wavenumber * y - c.omega * 3.0 + c.phase;
        exact = exact && std::abs(field.elevation(0.0, y, 3.0) - c.amplitude * std::cos(psi)) < 1e-14;
    }
    expectTrue("a wave heading +y is a function of y alone", exact);
    expectNear("moving along x does not change a +y wave", field.elevation(500.0, 4.0, 3.0),
               field.elevation(-500.0, 4.0, 3.0), 1e-12);
}

// --- Depth decay -------------------------------------------------------------

void testOrbitalVelocityDecaysExponentially() {
    SeaState sea;
    sea.significantHeight = 2.5;
    sea.peakPeriod = 7.0;
    sea.frequencyCount = 1;
    sea.directionCount = 1;
    const WaveField field(sea);
    const double k = field.components()[0].wavenumber;

    // The ratio of orbital speed at two depths is exactly exp(k (z1 - z2)).
    for (double z1 : {-1.0, -5.0, -12.0}) {
        const double z2 = z1 - 4.0;
        const Vec3 upper = field.velocity({3.0, 0.0, z1}, 0.9);
        const Vec3 lower = field.velocity({3.0, 0.0, z2}, 0.9);
        expectNear(label("orbital speed ratio is exp(k dz) at z", z1),
                   length(lower) / length(upper), std::exp(k * (z2 - z1)), 1e-12);
    }
    // A depth of half a wavelength is the textbook e^-pi = 4.3% of the surface.
    const double halfWave = -kPi / k;
    expectNear("half a wavelength down leaves e^-pi of the motion",
               length(field.velocity({3.0, 0.0, halfWave}, 0.9)) /
                   length(field.velocity({3.0, 0.0, 0.0}, 0.9)),
               std::exp(-kPi), 1e-12);

    // Above the still-water plane the exponential is clamped: a Froude-Krylov
    // integration over an instantaneous wetted surface will ask for z > 0 inside
    // a crest, and exp(k z) there grows without bound.
    const Vec3 above = field.velocity({3.0, 0.0, 6.0}, 0.9);
    const Vec3 at = field.velocity({3.0, 0.0, 0.0}, 0.9);
    expectNear("kinematics above still water are clamped to z = 0", length(above), length(at),
               1e-15);
}

// --- Statistics --------------------------------------------------------------

// The variance of the surface over many points and times must converge to m0.
// This is a different statement from "the amplitudes sum to m0", which is true by
// construction: it is the check that the phases really are spread and the
// components really do superpose incoherently.
//
// Sample count and tolerance, stated rather than tuned. The surface is Gaussian
// to a good approximation, so the variance estimator over M independent samples
// has a relative standard error of sqrt(2/M); at M = 300 000 that is 0.26%.
// Measured over 30 independent field and sampling seeds the RMS relative error
// came out at 0.22% and the worst at 0.50%, which is the expected distribution.
// The tolerance below is 1%, just under four standard errors.
//
// The sample points are pseudorandom rather than a lattice. A regular grid is
// the trap here: the first version of this test walked x in steps of 101 m
// through a sea whose dominant wavelength was 126 m, so the dominant components
// were aliased instead of averaged and the variance came out 1% high --
// systematically, not by chance. Random points have no resonance to hit.
void testSurfaceVarianceConvergesToM0() {
    SeaState sea;
    sea.significantHeight = 3.0;
    sea.peakPeriod = 9.0;
    sea.spreadingExponent = 8.0;
    sea.frequencyCount = 24;
    sea.directionCount = 6;
    sea.seed = 0xA1B2C3D4;
    const WaveField field(sea);

    // A plain 64-bit LCG, local to the test: the sample points must be
    // reproducible and must not come from the generator under test.
    std::uint64_t state = 0x9E3779B97F4A7C15ull;
    const auto uniform = [&state] {
        state = state * 6364136223846793005ull + 1442695040888963407ull;
        return static_cast<double>(state >> 11) * 0x1p-53;
    };

    constexpr long long kSamples = 300000;
    double sum = 0.0, sumSquares = 0.0;
    for (long long i = 0; i < kSamples; ++i) {
        // 20 km square, 3000 s: about 60 dominant wavelengths and 300 dominant
        // periods, so the domain average is itself m0 to well inside the
        // sampling error.
        const double x = uniform() * 20000.0;
        const double y = uniform() * 20000.0;
        const double t = uniform() * 3000.0;
        const double e = field.elevation(x, y, t);
        sum += e;
        sumSquares += e * e;
    }
    const double mean = sum / static_cast<double>(kSamples);
    const double variance = sumSquares / static_cast<double>(kSamples) - mean * mean;
    const double m0 = field.zerothMoment();
    expectNear("sampled surface variance converges to m0", variance, m0, 0.01 * m0);
    // The mean has standard error sqrt(m0/M) = 1.4 mm; 10 mm is seven of those.
    expectNear("sampled surface mean is zero", mean, 0.0, 0.01);
    // And the sampled Hs, which is what a wave buoy would report. Half the
    // relative tolerance of the variance, because of the square root.
    expectNear("sampled 4 sqrt(variance) recovers Hs", 4.0 * std::sqrt(variance), 3.0, 0.016);
}

// --- Determinism -------------------------------------------------------------

void testSeedDeterminism() {
    SeaState sea;
    sea.frequencyCount = 32;
    sea.directionCount = 8;
    sea.seed = 12345;

    const WaveField a(sea);
    const WaveField b(sea);
    SeaState other = sea;
    other.seed = 12346;
    const WaveField c(other);

    bool identical = true, differs = false;
    for (int i = 0; i < 200; ++i) {
        const double x = 3.1 * i, y = -1.7 * i, t = 0.21 * i;
        identical = identical && a.elevation(x, y, t) == b.elevation(x, y, t);
        differs = differs || a.elevation(x, y, t) != c.elevation(x, y, t);
    }
    expectTrue("the same seed gives bit-identical elevation", identical);
    expectTrue("a different seed gives a different sea", differs);

    bool samePhases = true, differentPhases = false;
    for (std::size_t i = 0; i < a.components().size(); ++i) {
        samePhases = samePhases && a.components()[i].phase == b.components()[i].phase;
        differentPhases = differentPhases || a.components()[i].phase != c.components()[i].phase;
        // Frequencies and directions are deterministic regardless of seed.
        samePhases = samePhases && a.components()[i].omega == c.components()[i].omega;
    }
    expectTrue("phases and geometry are reproduced exactly", samePhases);
    expectTrue("the seed changes the phases", differentPhases);

    // Phases must fill the circle rather than clustering: a generator seeded per
    // component but badly mixed would give a visibly banded sea.
    double phaseSum = 0.0, phaseSquares = 0.0;
    for (const WaveComponent& component : a.components()) {
        phaseSum += component.phase;
        phaseSquares += component.phase * component.phase;
    }
    const double n = static_cast<double>(a.components().size());
    const double phaseMean = phaseSum / n;
    expectNear("phases average pi", phaseMean, kPi, 0.25);
    expectNear("phase variance is (2 pi)^2 / 12", phaseSquares / n - phaseMean * phaseMean,
               4.0 * kPi * kPi / 12.0, 0.6);
}

// --- Degenerate and composite fields -----------------------------------------

void testDegenerateFields() {
    const WaveField still;
    expectTrue("a default wave field has no components", still.components().empty());
    expectNear("still water is flat", still.elevation(12.0, -4.0, 3.0), 0.0, 0.0);
    expectNear("still water has no orbital motion", length(still.velocity({1, 2, -3}, 4.0)), 0.0,
               0.0);
    expectNear("still water has zero significant height", still.significantHeight(), 0.0, 0.0);

    SeaState calm;
    calm.significantHeight = 0.0;
    const WaveField flat(calm);
    expectNear("Hs = 0 gives a flat sea", flat.elevation(30.0, 20.0, 10.0), 0.0, 0.0);
    expectNear("Hs = 0 gives zero m0", flat.zerothMoment(), 0.0, 0.0);

    SeaState degenerate;
    degenerate.frequencyCount = 0;
    degenerate.directionCount = -3;
    const WaveField clamped(degenerate);
    expectEqual("a degenerate component count is clamped to one",
                static_cast<long long>(clamped.components().size()), 1);

    // A non-positive peak period has no peak frequency to place anything at, so
    // the sea state contributes nothing rather than a NaN.
    SeaState noPeriod;
    noPeriod.peakPeriod = 0.0;
    const WaveField empty(noPeriod);
    expectEqual("a zero peak period contributes no components",
                static_cast<long long>(empty.components().size()), 0);
    expectNear("a zero peak period leaves the sea flat", empty.elevation(5.0, 6.0, 7.0), 0.0, 0.0);
    expectNear("a zero peak period has zero density", Spectrum(noPeriod).density(1.0), 0.0, 0.0);
}

// Wind sea plus swell: variances add, because the two systems are independent.
void testSuperposedSeaStates() {
    SeaState wind;
    wind.significantHeight = 2.5;
    wind.peakPeriod = 7.0;
    wind.frequencyCount = 24;
    wind.directionCount = 6;
    wind.seed = 11;

    SeaState swell;
    swell.significantHeight = 1.5;
    swell.peakPeriod = 14.0;
    swell.meanDirection = 1.2;
    swell.spreadingExponent = 60.0;
    swell.frequencyCount = 16;
    swell.directionCount = 4;
    swell.seed = 11;  // same seed: the streams must still differ

    const WaveField combined(std::vector<SeaState>{wind, swell});
    expectEqual("components of both systems are present",
                static_cast<long long>(combined.components().size()), 24 * 6 + 16 * 4);
    const double want = (2.5 * 2.5 + 1.5 * 1.5) / 16.0;
    expectNear("superposed variances add", combined.zerothMoment(), want, 1e-12 * want);
    expectNear("the combined Hs is the root sum of squares", combined.significantHeight(),
               std::sqrt(2.5 * 2.5 + 1.5 * 1.5), 1e-11);

    // Same seed, different stream: the swell must not be a copy of the wind sea's
    // phase sequence.
    const WaveField windOnly(wind);
    bool sharedPhase = false;
    for (std::size_t i = 0; i < 64; ++i)
        sharedPhase =
            sharedPhase || combined.components()[24 * 6 + i].phase == windOnly.components()[i].phase;
    expectTrue("a second sea state draws from its own phase stream", !sharedPhase);
}

}  // namespace

// --- Regular waves -----------------------------------------------------------

// A monochromatic train is the test signal for a response amplitude operator, so
// it has to be exactly the wave it claims to be, not approximately.
void testRegularWaveIsTheExactAiryForm() {
    const double amplitude = 1.75, omega = 0.7, phase = 0.4;
    const WaveField field = WaveField::regular(amplitude, omega, 0.0, phase);
    expectEqual("a regular wave has exactly one component",
                static_cast<long long>(field.components().size()), 1LL);

    const double k = omega * omega / kGravity;
    expectNear("the wavenumber comes from the dispersion relation",
               field.components()[0].wavenumber, k, 1e-12);

    double worst = 0;
    for (double t : {0.0, 1.3, 7.9, 40.0})
        for (double x : {-120.0, -13.0, 0.0, 55.5, 310.0}) {
            const double want = amplitude * std::cos(k * x - omega * t + phase);
            worst = std::max(worst, std::abs(field.elevation(x, 0.0, t) - want));
        }
    expectTrue("elevation is A cos(kx - wt + phase) everywhere", worst < 1e-12);

    // A wave along +x must not vary along y at all -- not "to tolerance".
    double crossVariation = 0;
    for (double y : {-500.0, -20.0, 20.0, 500.0})
        crossVariation = std::max(crossVariation,
                                  std::abs(field.elevation(30.0, y, 2.0) -
                                           field.elevation(30.0, 0.0, 2.0)));
    expectTrue("a long-crested wave along +x is constant along y", crossVariation == 0.0);

    // m0 = A^2/2 for a single component, so Hs = 2 sqrt(2) A. This is worth
    // asserting because it is the one statistic the spectral path and the
    // explicit path must agree on.
    expectNear("significant height of a regular wave is 2 sqrt(2) A",
               field.significantHeight(), 2.0 * std::sqrt(2.0) * amplitude, 1e-12);
}

// The crest must travel at the phase speed g/omega. A wave of the right shape
// moving at the wrong speed is the failure this file was written to catch, and a
// regular wave makes it checkable exactly rather than statistically.
void testRegularWaveCrestTravelsAtThePhaseSpeed() {
    const double omega = 0.9;
    const WaveField field = WaveField::regular(2.0, omega, 0.0, 0.0);
    const double c = kGravity / omega;

    // The crest at t = 0 sits at x = 0; at time t it must sit at x = c t.
    double worst = 0;
    for (double t : {0.0, 3.0, 11.0, 25.0}) {
        const double crestX = c * t;
        worst = std::max(worst, std::abs(field.elevation(crestX, 0.0, t) - 2.0));
    }
    expectTrue("the crest stays a crest when followed at g/omega", worst < 1e-9);

    // And is demonstrably not a crest when followed at the group speed, which is
    // half of it -- otherwise the check above would pass for a standing wave.
    expectTrue("following at the group speed instead does not track the crest",
               std::abs(field.elevation(0.5 * c * 25.0, 0.0, 25.0) - 2.0) > 0.5);
}

// The explicit-component constructor recomputes the derived fields. Hand it a
// wavenumber that contradicts omega and check it is overruled: a component that
// kept it would propagate at a forbidden speed and look entirely normal.
void testExplicitComponentsCannotBreakDispersion() {
    WaveComponent bogus;
    bogus.amplitude = 1.0;
    bogus.omega = 1.2;
    bogus.direction = kPi / 3.0;
    bogus.wavenumber = 99.0;   // nonsense
    bogus.dirX = 0.0;          // inconsistent with direction
    bogus.dirY = 0.0;
    const WaveField field(std::vector<WaveComponent>{bogus});

    const WaveComponent& c = field.components()[0];
    expectNear("a contradictory wavenumber is recomputed from omega", c.wavenumber,
               1.2 * 1.2 / kGravity, 1e-12);
    expectNear("the direction unit vector is recomputed too", c.dirX, std::cos(kPi / 3.0), 1e-12);
    expectNear("and its y part", c.dirY, std::sin(kPi / 3.0), 1e-12);

    // A zeroed direction vector would have made elevation independent of
    // position, so this also proves the recomputation reached the evaluator.
    expectTrue("the field actually varies in space after correction",
               std::abs(field.elevation(20.0, 5.0, 0.0) - field.elevation(0.0, 0.0, 0.0)) > 1e-6);
}

void runWaveTests() {
    std::printf("\n--- spectral wave field ---\n");
    testSpectrumIntegratesToTheRequestedHs();
    testPeakEnhancementNormalisation();
    testPeakEnhancementShapeIsAsymmetric();
    testSpectrumPeaksAtTheRequestedPeriod();
    testDeepWaterDispersion();
    testSignificantHeightRoundTrip();
    testDiscretisationMatchesTheAnalyticIntegral();
    testFrequencyBinEdgesAgainstTheClosedForm();
    testEachBinHoldsTheEnergyItClaims();
    testDiscretisedSpectrumPeaksAtTp();
    testMeanPeriodsAgainstGammaFunctions();
    testSpreadingNormalisation();
    testSpreadingFirstMoment();
    testDirectionBinsPartitionTheEnergy();
    testLargeSpreadingExponentIsLongCrested();
    testSingleComponentIsAnExactAiryWave();
    testDirectionRotatesTheWave();
    testOrbitalVelocityDecaysExponentially();
    testSurfaceVarianceConvergesToM0();
    testSeedDeterminism();
    testDegenerateFields();
    testSuperposedSeaStates();
    testRegularWaveIsTheExactAiryForm();
    testRegularWaveCrestTravelsAtThePhaseSpeed();
    testExplicitComponentsCannotBreakDispersion();
}
