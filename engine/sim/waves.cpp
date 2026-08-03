// SPDX-License-Identifier: MIT
#include "waves.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

namespace sim {
namespace {

// --- Counter-based RNG ------------------------------------------------------
//
// Threefry-2x64-20, as `01-architecture.md` section 4 requires: stateless,
// explicitly keyed, indexed by a counter rather than advanced by a stream. The
// phase of component i is a pure function of (seed, stream, i), so components can
// be generated in any order, on any thread, and a replay of the same seed
// reproduces the sea bit for bit.

constexpr std::uint64_t kThreefryParity = 0x1BD11BDAA9FC1A22ull;
constexpr int kThreefryRotation[8] = {16, 42, 12, 31, 16, 32, 24, 21};

std::uint64_t rotateLeft(std::uint64_t value, int bits) {
    return (value << bits) | (value >> (64 - bits));
}

std::array<std::uint64_t, 2> threefry2x64(std::uint64_t key0, std::uint64_t key1,
                                          std::uint64_t counter0, std::uint64_t counter1) {
    const std::uint64_t schedule[3] = {key0, key1, kThreefryParity ^ key0 ^ key1};
    std::uint64_t x0 = counter0 + schedule[0];
    std::uint64_t x1 = counter1 + schedule[1];
    for (int round = 0; round < 20; ++round) {
        x0 += x1;
        x1 = rotateLeft(x1, kThreefryRotation[round & 7]);
        x1 ^= x0;
        if ((round & 3) == 3) {
            const int inject = round / 4 + 1;
            x0 += schedule[inject % 3];
            x1 += schedule[(inject + 1) % 3] + static_cast<std::uint64_t>(inject);
        }
    }
    return {x0, x1};
}

// Uniform phase in [0, 2 pi). The top 53 bits are used, which is the whole
// mantissa; the low bits of a Threefry output are as good as the high ones, but
// taking the high ones keeps this identical to the usual convention.
double phaseFrom(std::uint64_t seed, std::uint64_t stream, std::uint64_t index) {
    const auto bits = threefry2x64(seed, stream, index, 0);
    const double unit = static_cast<double>(bits[0] >> 11) * 0x1p-53;
    return unit * 2.0 * kPi;
}

// --- Quadrature -------------------------------------------------------------

// Panels per bin when integrating a bin's energy and centroid. The integrands
// are smooth over an interval carrying 1/N of the variance, so this is far more
// than accuracy needs and it is paid once, at construction.
constexpr int kBinPanels = 256;

// Composite Simpson with an explicit midpoint per panel, so the rule is exact
// for cubics on every panel and the cumulative form below can report a value at
// every node rather than every second node.
template <class F>
double simpson(F f, double lo, double hi, int panels) {
    if (!(hi > lo) || panels < 1) return 0.0;
    const double h = (hi - lo) / panels;
    double sum = 0.0;
    for (int i = 0; i < panels; ++i) {
        const double a = lo + h * i;
        const double b = (i + 1 == panels) ? hi : lo + h * (i + 1);
        sum += (b - a) / 6.0 * (f(a) + 4.0 * f(0.5 * (a + b)) + f(b));
    }
    return sum;
}

// --- JONSWAP shape ----------------------------------------------------------
//
// Everything is done in the reduced frequency x = omega / omega_p, where the
// shape is
//
//     phi(x) = x^-5 exp(-1.25 x^-4) gamma^r(x),
//     r(x)   = exp(-(x - 1)^2 / (2 sigma^2)),   sigma = 0.07 (x <= 1), 0.09 (x > 1)
//
// and is independent of Hs, Tp and the units. Two facts make the discretisation
// cheap and accurate:
//
//   * the peak-enhancement factor is 1 to within 1e-14 outside [0.40, 1.90], so
//     outside that window the shape is the bare Pierson-Moskowitz form;
//   * that bare form integrates in closed form,
//     integral of phi from 0 to X = (1/5) exp(-1.25 X^-4),
//     which is also the Pierson-Moskowitz cumulative distribution of energy.
//
// So both tails are exact and only the middle needs quadrature.

constexpr double kEnhancementLow = 0.40;
constexpr double kEnhancementHigh = 1.90;

double enhancement(double x, double gamma) {
    if (gamma == 1.0 || x < kEnhancementLow || x > kEnhancementHigh) return 1.0;
    const double sigma = x <= 1.0 ? 0.07 : 0.09;
    const double d = (x - 1.0) / sigma;
    return std::pow(gamma, std::exp(-0.5 * d * d));
}

double shape(double x, double gamma) {
    // Below x = 0.05 the exponential is exp(-2e5) and the shape is identically
    // zero in double precision; the guard also keeps x^-5 away from overflow.
    if (x < 0.05) return 0.0;
    const double x2 = x * x;
    const double x4 = x2 * x2;
    return std::exp(-1.25 / x4) / (x4 * x) * enhancement(x, gamma);
}

// Closed-form integral of the bare shape from 0 to x, and from x to infinity.
double bareShapeBelow(double x) {
    if (x < 0.05) return 0.0;
    const double x2 = x * x;
    return 0.2 * std::exp(-1.25 / (x2 * x2));
}
double bareShapeAbove(double x) { return 0.2 - bareShapeBelow(x); }

// Integral of phi over (0, infinity). Equals 1/5 exactly when gamma = 1.
double shapeIntegral(double gamma) {
    const auto phi = [gamma](double x) { return shape(x, gamma); };
    // Split at the peak: sigma steps there, so the second derivative does too.
    return bareShapeBelow(kEnhancementLow) +
           simpson(phi, kEnhancementLow, 1.0, 4000) +
           simpson(phi, 1.0, kEnhancementHigh, 4000) + bareShapeAbove(kEnhancementHigh);
}

// Cumulative energy distribution of the shape, invertible. Closed form outside
// the enhancement window, tabulated inside it.
class ShapeCdf {
public:
    explicit ShapeCdf(double gamma) : gamma_(gamma) {
        const auto phi = [gamma](double x) { return shape(x, gamma); };
        below_ = bareShapeBelow(kEnhancementLow);
        node_.reserve(2 * kPanels + 1);
        cumulative_.reserve(2 * kPanels + 1);
        double running = 0.0;
        node_.push_back(kEnhancementLow);
        cumulative_.push_back(0.0);
        const double edges[3] = {kEnhancementLow, 1.0, kEnhancementHigh};
        for (int half = 0; half < 2; ++half) {
            const double h = (edges[half + 1] - edges[half]) / kPanels;
            for (int i = 0; i < kPanels; ++i) {
                const double a = node_.back();
                const double b = (i + 1 == kPanels) ? edges[half + 1] : edges[half] + h * (i + 1);
                running += (b - a) / 6.0 * (phi(a) + 4.0 * phi(0.5 * (a + b)) + phi(b));
                node_.push_back(b);
                cumulative_.push_back(running);
            }
        }
        total_ = below_ + running + bareShapeAbove(kEnhancementHigh);
    }

    // x such that the fraction of energy below it is p.
    double inverse(double p) const {
        if (p <= 0.0) return 0.0;
        if (p >= 1.0) return std::numeric_limits<double>::infinity();
        const double target = p * total_;
        if (target <= below_) {
            // (1/5) exp(-1.25 x^-4) = target
            const double ln = std::log(5.0 * target);
            return std::pow(1.25 / -ln, 0.25);
        }
        if (target >= below_ + cumulative_.back()) {
            // (1/5) (1 - exp(-1.25 x^-4)) = total - target
            const double inner = 1.0 - 5.0 * (total_ - target);
            if (inner <= 0.0) return std::numeric_limits<double>::infinity();
            return std::pow(1.25 / -std::log(inner), 0.25);
        }
        const double want = target - below_;
        std::size_t lo = 0, hi = cumulative_.size() - 1;
        while (hi - lo > 1) {
            const std::size_t mid = (lo + hi) / 2;
            if (cumulative_[mid] <= want) lo = mid; else hi = mid;
        }
        const double xa = node_[lo], xb = node_[hi];
        const double ca = cumulative_[lo], cb = cumulative_[hi];
        double x = xa + (xb - xa) * (want - ca) / (cb - ca);
        const auto phi = [g = gamma_](double v) { return shape(v, g); };
        for (int iteration = 0; iteration < 4; ++iteration) {
            const double slope = phi(x);
            if (!(slope > 0.0)) break;
            x -= (ca + simpson(phi, xa, x, 2) - want) / slope;
            x = std::clamp(x, xa, xb);
        }
        return x;
    }

private:
    static constexpr int kPanels = 2000;
    double gamma_ = 1.0;
    double below_ = 0.0;
    double total_ = 0.2;
    std::vector<double> node_;
    std::vector<double> cumulative_;
};

// Energy and first moment of the shape over [xLow, xHigh], xHigh possibly
// infinite. The infinite piece is integrated in y = 1/x, where phi dx becomes
// y^3 exp(-1.25 y^4) dy and x phi dx becomes y^2 exp(-1.25 y^4) dy -- both
// bounded and smooth at y = 0, whereas x^-4 in the original variable has an
// unbounded derivative at the far end and Simpson would lose most of its order.
void shapeEnergyAndMoment(double xLow, double xHigh, double gamma, double& energy,
                          double& moment) {
    const auto phi = [gamma](double x) { return shape(x, gamma); };
    const auto xPhi = [gamma](double x) { return x * shape(x, gamma); };
    energy = 0.0;
    moment = 0.0;
    if (std::isfinite(xHigh)) {
        energy = simpson(phi, xLow, xHigh, kBinPanels);
        moment = simpson(xPhi, xLow, xHigh, kBinPanels);
        return;
    }
    const double split = std::max(xLow, 4.0);
    if (split > xLow) {
        energy += simpson(phi, xLow, split, kBinPanels);
        moment += simpson(xPhi, xLow, split, kBinPanels);
    }
    const auto energyInY = [gamma](double y) {
        if (y <= 0.0) return 0.0;
        const double y2 = y * y;
        return y2 * y * std::exp(-1.25 * y2 * y2) * enhancement(1.0 / y, gamma);
    };
    const auto momentInY = [gamma](double y) {
        if (y <= 0.0) return 0.0;
        const double y2 = y * y;
        return y2 * std::exp(-1.25 * y2 * y2) * enhancement(1.0 / y, gamma);
    };
    energy += simpson(energyInY, 0.0, 1.0 / split, kBinPanels);
    moment += simpson(momentInY, 0.0, 1.0 / split, kBinPanels);
}

// --- cos^(2s) spreading -----------------------------------------------------
//
// Written in the half-angle t = (theta - theta0) / 2, where the density is
// cos^(2s)(t) on [-pi/2, pi/2]. That is the Longuet-Higgins form: it covers the
// whole circle, unlike cos^(2s)(theta - theta0) which has to be truncated at
// +-pi/2 and then renormalised.

double spreadDensity(double t, double s) {
    const double c = std::cos(t);
    if (c <= 0.0) return s > 0.0 ? 0.0 : 1.0;
    return std::exp(2.0 * s * std::log(c));
}

// Half-width outside which the density is below exp(-80) of its peak, clamped to
// the half-angle range. A uniform grid over the full range would miss the peak
// entirely for large s: at s = 5e5 the density is narrower than 0.01 rad.
double spreadHalfWidth(double s) {
    if (!(s > 0.0)) return 0.5 * kPi;
    const double c = std::exp(-40.0 / s);
    return std::acos(std::clamp(c, 0.0, 1.0));
}

}  // namespace

// --- Dispersion -------------------------------------------------------------

double deepWaterWavenumber(double omega) {
    return omega > 0.0 ? omega * omega / kGravity : 0.0;
}
double phaseSpeed(double omega) { return omega > 0.0 ? kGravity / omega : 0.0; }
double groupSpeed(double omega) { return 0.5 * phaseSpeed(omega); }

double directionalSpreading(double theta, double meanDirection, double s) {
    const double s0 = std::max(s, 0.0);
    const double delta = std::remainder(theta - meanDirection, 2.0 * kPi);
    // C(s) = Gamma(s+1) / (2 sqrt(pi) Gamma(s+1/2)), which is what makes the
    // density integrate to one over the circle.
    const double norm = std::exp(std::lgamma(s0 + 1.0) - std::lgamma(s0 + 0.5)) /
                        (2.0 * std::sqrt(kPi));
    return norm * spreadDensity(0.5 * delta, s0);
}

// --- Spectrum ---------------------------------------------------------------

Spectrum::Spectrum(const SeaState& sea)
    : omegaPeak_(sea.peakPeriod > 0.0 ? 2.0 * kPi / sea.peakPeriod : 0.0),
      gamma_(sea.shape == SpectrumShape::PiersonMoskowitz ? 1.0
                                                          : std::max(sea.peakEnhancement, 1e-6)),
      significantHeight_(std::max(sea.significantHeight, 0.0)) {
    normalisation_ = 0.2 / shapeIntegral(gamma_);
}

double Spectrum::density(double omega) const {
    if (!(omega > 0.0) || !(omegaPeak_ > 0.0)) return 0.0;
    const double hs = significantHeight_;
    return 5.0 / 16.0 * hs * hs * normalisation_ * shape(omega / omegaPeak_, gamma_) / omegaPeak_;
}

double Spectrum::zerothMoment() const {
    return significantHeight_ * significantHeight_ / 16.0;
}

// --- WaveField --------------------------------------------------------------

WaveField::WaveField(const SeaState& sea) { addSeaState(sea, 0); }

WaveField::WaveField(const std::vector<SeaState>& seas) {
    for (std::size_t i = 0; i < seas.size(); ++i) addSeaState(seas[i], i);
}

WaveField::WaveField(std::vector<WaveComponent> components)
    : components_(std::move(components)) {
    // Recompute rather than trust: wavenumber and the direction unit vector are
    // functions of omega and direction, and letting a caller set them
    // independently would allow a wave that propagates at a speed the dispersion
    // relation forbids -- visually and statistically indistinguishable from a
    // correct one.
    for (WaveComponent& c : components_) {
        c.wavenumber = deepWaterWavenumber(c.omega);
        c.dirX = std::cos(c.direction);
        c.dirY = std::sin(c.direction);
    }
}

WaveField WaveField::regular(double amplitude, double omega, double direction, double phase) {
    WaveComponent c;
    c.amplitude = amplitude;
    c.omega = omega;
    c.direction = direction;
    c.phase = phase;
    return WaveField(std::vector<WaveComponent>{c});
}

void WaveField::addSeaState(const SeaState& sea, std::uint64_t stream) {
    const Spectrum spectrum(sea);
    const int frequencies = std::max(sea.frequencyCount, 1);
    const int directions = std::max(sea.directionCount, 1);
    const double omegaPeak = spectrum.peakFrequency();
    const double gamma = spectrum.peakEnhancement();
    const double m0 = spectrum.zerothMoment();
    if (!(omegaPeak > 0.0)) return;

    // Frequency: equal-energy intervals, each placed at its energy centroid.
    const ShapeCdf cdf(gamma);
    std::vector<double> edge(static_cast<std::size_t>(frequencies) + 1, 0.0);
    for (int i = 0; i <= frequencies; ++i)
        edge[static_cast<std::size_t>(i)] =
            cdf.inverse(static_cast<double>(i) / frequencies);

    const std::size_t firstFrequencyBin = frequencyBins_.size();
    for (int i = 0; i < frequencies; ++i) {
        const double xLow = edge[static_cast<std::size_t>(i)];
        const double xHigh = edge[static_cast<std::size_t>(i) + 1];
        double energy = 0.0, moment = 0.0;
        shapeEnergyAndMoment(xLow, xHigh, gamma, energy, moment);
        FrequencyBin bin;
        bin.omegaLow = xLow * omegaPeak;
        bin.omegaHigh = xHigh * omegaPeak;
        bin.omega = (energy > 0.0 ? moment / energy : 0.5 * (xLow + xHigh)) * omegaPeak;
        bin.energy = m0 / frequencies;
        frequencyBins_.push_back(bin);
    }

    // Direction: the same construction on the spreading function.
    const double s = std::max(sea.spreadingExponent, 0.0);
    const double halfWidth = spreadHalfWidth(s);
    constexpr int kSpreadPanels = 4000;
    std::vector<double> spreadNode(kSpreadPanels + 1, 0.0);
    std::vector<double> spreadCumulative(kSpreadPanels + 1, 0.0);
    {
        const double h = 2.0 * halfWidth / kSpreadPanels;
        double running = 0.0;
        spreadNode[0] = -halfWidth;
        for (int i = 0; i < kSpreadPanels; ++i) {
            const double a = spreadNode[static_cast<std::size_t>(i)];
            const double b = (i + 1 == kSpreadPanels) ? halfWidth : -halfWidth + h * (i + 1);
            running += (b - a) / 6.0 *
                       (spreadDensity(a, s) + 4.0 * spreadDensity(0.5 * (a + b), s) +
                        spreadDensity(b, s));
            spreadNode[static_cast<std::size_t>(i) + 1] = b;
            spreadCumulative[static_cast<std::size_t>(i) + 1] = running;
        }
    }
    const double spreadTotal = spreadCumulative.back();

    const auto spreadInverse = [&](double p) {
        if (p <= 0.0) return -halfWidth;
        if (p >= 1.0) return halfWidth;
        const double want = p * spreadTotal;
        std::size_t lo = 0, hi = spreadCumulative.size() - 1;
        while (hi - lo > 1) {
            const std::size_t mid = (lo + hi) / 2;
            if (spreadCumulative[mid] <= want) lo = mid; else hi = mid;
        }
        const double ta = spreadNode[lo], tb = spreadNode[hi];
        const double ca = spreadCumulative[lo], cb = spreadCumulative[hi];
        double t = ta + (tb - ta) * (want - ca) / (cb - ca);
        for (int iteration = 0; iteration < 4; ++iteration) {
            const double slope = spreadDensity(t, s);
            if (!(slope > 0.0)) break;
            t -= (ca + simpson([s](double v) { return spreadDensity(v, s); }, ta, t, 2) - want) /
                 slope;
            t = std::clamp(t, ta, tb);
        }
        return t;
    };

    std::vector<double> spreadEdge(static_cast<std::size_t>(directions) + 1, 0.0);
    for (int j = 0; j <= directions; ++j)
        spreadEdge[static_cast<std::size_t>(j)] =
            spreadInverse(static_cast<double>(j) / directions);
    // The spreading is symmetric about the mean direction, so the edges must be
    // too. Enforcing it rather than hoping for it is what makes a single-
    // direction sea exactly long-crested instead of long-crested to 1e-17, and
    // it costs nothing.
    for (int j = 0; j <= directions / 2; ++j) {
        const auto a = static_cast<std::size_t>(j);
        const auto b = static_cast<std::size_t>(directions - j);
        const double halfDifference = 0.5 * (spreadEdge[a] - spreadEdge[b]);
        spreadEdge[a] = halfDifference;
        spreadEdge[b] = -halfDifference;
    }

    const std::size_t firstDirectionBin = directionBins_.size();
    std::vector<double> centre(static_cast<std::size_t>(directions), 0.0);
    for (int j = 0; j < directions; ++j) {
        const double tLow = spreadEdge[static_cast<std::size_t>(j)];
        const double tHigh = spreadEdge[static_cast<std::size_t>(j) + 1];
        const double energy =
            simpson([s](double v) { return spreadDensity(v, s); }, tLow, tHigh, kBinPanels);
        const double moment =
            simpson([s](double v) { return v * spreadDensity(v, s); }, tLow, tHigh, kBinPanels);
        centre[static_cast<std::size_t>(j)] =
            energy > 0.0 ? moment / energy : 0.5 * (tLow + tHigh);
    }
    for (int j = 0; j <= (directions - 1) / 2; ++j) {
        const auto a = static_cast<std::size_t>(j);
        const auto b = static_cast<std::size_t>(directions - 1 - j);
        const double halfDifference = 0.5 * (centre[a] - centre[b]);
        centre[a] = halfDifference;
        centre[b] = -halfDifference;
    }
    for (int j = 0; j < directions; ++j) {
        DirectionBin bin;
        // t is the half-angle, so the direction offset is twice it.
        bin.directionLow = sea.meanDirection + 2.0 * spreadEdge[static_cast<std::size_t>(j)];
        bin.directionHigh = sea.meanDirection + 2.0 * spreadEdge[static_cast<std::size_t>(j) + 1];
        bin.direction = sea.meanDirection + 2.0 * centre[static_cast<std::size_t>(j)];
        bin.energy = 1.0 / directions;
        directionBins_.push_back(bin);
    }

    // Components: one per (frequency, direction) pair, each carrying m0 / (N M)
    // of the variance, so the amplitude is the same for all of them.
    const double amplitude = std::sqrt(2.0 * m0 / (frequencies * directions));
    components_.reserve(components_.size() +
                        static_cast<std::size_t>(frequencies) * directions);
    for (int i = 0; i < frequencies; ++i) {
        const FrequencyBin& fbin = frequencyBins_[firstFrequencyBin + static_cast<std::size_t>(i)];
        for (int j = 0; j < directions; ++j) {
            const DirectionBin& dbin =
                directionBins_[firstDirectionBin + static_cast<std::size_t>(j)];
            WaveComponent component;
            component.amplitude = amplitude;
            component.omega = fbin.omega;
            component.wavenumber = deepWaterWavenumber(fbin.omega);
            component.direction = dbin.direction;
            component.dirX = std::cos(dbin.direction);
            component.dirY = std::sin(dbin.direction);
            component.phase = phaseFrom(
                sea.seed, stream, static_cast<std::uint64_t>(i) * static_cast<std::uint64_t>(
                                      directions) + static_cast<std::uint64_t>(j));
            components_.push_back(component);
        }
    }
}

double WaveField::elevation(double x, double y, double t) const {
    double sum = 0.0;
    for (const WaveComponent& c : components_) {
        const double psi = c.wavenumber * (c.dirX * x + c.dirY * y) - c.omega * t + c.phase;
        sum += c.amplitude * std::cos(psi);
    }
    return sum;
}

WaveKinematics WaveField::kinematics(const Vec3& position, double t) const {
    WaveKinematics out;
    const double z = std::min(position.z, 0.0);
    for (const WaveComponent& c : components_) {
        const double psi =
            c.wavenumber * (c.dirX * position.x + c.dirY * position.y) - c.omega * t + c.phase;
        const double cosPsi = std::cos(psi);
        const double sinPsi = std::sin(psi);
        out.elevation += c.amplitude * cosPsi;
        const double decay = std::exp(c.wavenumber * z);
        const double along = c.amplitude * c.omega * decay * cosPsi;  // horizontal, along dir
        const double up = c.amplitude * c.omega * decay * sinPsi;
        out.velocity += Vec3{c.dirX * along, c.dirY * along, up};
        const double alongDot = c.amplitude * c.omega * c.omega * decay * sinPsi;
        const double upDot = -c.amplitude * c.omega * c.omega * decay * cosPsi;
        out.acceleration += Vec3{c.dirX * alongDot, c.dirY * alongDot, upDot};
    }
    return out;
}

Vec3 WaveField::velocity(const Vec3& position, double t) const {
    return kinematics(position, t).velocity;
}
Vec3 WaveField::acceleration(const Vec3& position, double t) const {
    return kinematics(position, t).acceleration;
}

double WaveField::spectralMoment(int order) const {
    double sum = 0.0;
    for (const WaveComponent& c : components_) {
        double weight = 0.5 * c.amplitude * c.amplitude;
        for (int i = 0; i < order; ++i) weight *= c.omega;
        sum += weight;
    }
    return sum;
}

double WaveField::significantHeight() const { return 4.0 * std::sqrt(spectralMoment(0)); }

double WaveField::meanPeriod() const {
    const double m0 = spectralMoment(0), m1 = spectralMoment(1);
    return m1 > 0.0 ? 2.0 * kPi * m0 / m1 : 0.0;
}

double WaveField::zeroCrossingPeriod() const {
    const double m0 = spectralMoment(0), m2 = spectralMoment(2);
    return m2 > 0.0 ? 2.0 * kPi * std::sqrt(m0 / m2) : 0.0;
}

}  // namespace sim
