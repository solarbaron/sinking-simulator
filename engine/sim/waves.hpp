// SPDX-License-Identifier: MIT
//
// Directional wave spectrum, discretised into components and evaluated
// analytically.
//
// The sea is a sum of long-crested Airy waves:
//
//     eta(x, y, t) = sum_i a_i cos(k_i . x - omega_i t + phi_i)
//
// with amplitudes drawn from a JONSWAP (or Pierson-Moskowitz) spectrum spread
// over direction by a cos^(2s) law, and phases from a seeded counter-based
// generator. It is evaluated directly at whatever point is asked for rather than
// interpolated off a grid, because the points that matter -- a hatch coaming, a
// breach, a hull panel -- are exactly where a grid loses the detail that decides
// whether they immerse.
//
// **Frequency discretisation is equal-energy, not uniform.** The band is cut
// into N intervals each carrying exactly m0/N of the variance, and each interval
// is represented by its own energy centroid frequency. Three properties follow,
// and they are the reason for the choice:
//
//   * No energy is lost anywhere, including in the tail. The outermost intervals
//     run to zero and to infinity, so there is no truncation: the discrete
//     zeroth moment equals the analytic one to machine precision and
//     Hs = 4 sqrt(m0) round-trips exactly. Uniform spacing over a truncated band
//     loses the tail quietly, and a lost tail is a sea that is too smooth at the
//     scale of the openings that flood a ship.
//   * Components cluster where the energy is, near the peak, instead of being
//     spent resolving a tail that carries almost none.
//   * Because each interval is represented by its energy *centroid*, the first
//     spectral moment is preserved exactly too, so the mean period is right and
//     not merely close.
//
// The cost of the choice is that the second moment is biased slightly low -- a
// centroid cannot represent the spread of frequencies inside its own interval --
// so the zero-crossing period comes out a percent or two long at N = 48. That is
// stated in the tests rather than hidden.
//
// Directions are discretised the same way: M intervals of equal directional
// energy, each represented by its centroid. Kinematics are deep-water linear
// theory with exp(k z) depth decay.
#pragma once

#include "../core/math.hpp"

#include <cstdint>
#include <vector>

namespace sim {

enum class SpectrumShape {
    Jonswap,            // fetch-limited, peak-enhanced
    PiersonMoskowitz,   // fully developed; identical to JONSWAP with gamma = 1
};

// One wave system. Superpose several of these for a wind sea plus swell trains.
struct SeaState {
    double significantHeight = 3.0;   // Hs, m
    double peakPeriod = 9.0;          // Tp, s
    double meanDirection = 0.0;       // direction of *propagation*, rad, world frame
    double spreadingExponent = 10.0;  // s in cos^(2s)(theta/2); larger is longer-crested
    double peakEnhancement = 3.3;     // gamma; ignored for Pierson-Moskowitz
    SpectrumShape shape = SpectrumShape::Jonswap;
    int frequencyCount = 48;          // N
    int directionCount = 12;          // M
    std::uint64_t seed = 0x5eaf00d;   // phases are a pure function of this
};

// Deep-water linear dispersion. omega^2 = g k, so everything else follows.
double deepWaterWavenumber(double omega);
double phaseSpeed(double omega);   // c  = omega / k = g / omega
double groupSpeed(double omega);   // cg = c / 2 in deep water

// Normalised cos^(2s) spreading about `meanDirection`, integrating to exactly 1
// over any 2*pi interval of theta. s = 0 is isotropic; s -> infinity is
// long-crested.
double directionalSpreading(double theta, double meanDirection, double s);

// The analytic spectrum S(omega), in m^2 s/rad, normalised so that its zeroth
// moment is exactly Hs^2/16 whatever gamma is. The usual DNV shortcut for that
// normalisation (1 - 0.287 ln gamma) is an approximation; this integrates the
// shape instead, which is exact and costs one quadrature at construction.
class Spectrum {
public:
    explicit Spectrum(const SeaState& sea);

    double density(double omega) const;             // S(omega)
    double peakFrequency() const { return omegaPeak_; }
    double peakEnhancement() const { return gamma_; }
    // The factor that makes the shape integrate to Hs^2/16. Exactly 1 for
    // gamma = 1, and close to the DNV approximation for gamma = 3.3.
    double normalisation() const { return normalisation_; }
    double zerothMoment() const;                    // Hs^2 / 16, by construction

private:
    double omegaPeak_ = 0;
    double gamma_ = 1;
    double significantHeight_ = 0;
    double normalisation_ = 1;
};

// An equal-energy frequency interval. `omega` is the interval's energy centroid,
// which is the frequency the component is placed at. The outermost intervals are
// open: omegaLow is 0 for the first and omegaHigh is infinity for the last.
struct FrequencyBin {
    double omegaLow = 0;
    double omegaHigh = 0;
    double omega = 0;    // energy centroid
    double energy = 0;   // m^2; equal for every bin by construction
};

// An equal-energy direction interval, in the same shape.
struct DirectionBin {
    double directionLow = 0;
    double directionHigh = 0;
    double direction = 0;  // energy centroid
    double energy = 0;     // fraction of the sea state's energy; 1/M by construction
};

struct WaveComponent {
    double amplitude = 0;   // m
    double omega = 0;       // rad/s
    double wavenumber = 0;  // rad/m, = omega^2 / g
    double phase = 0;       // rad
    double direction = 0;   // rad, direction of propagation
    double dirX = 1;        // cos(direction)
    double dirY = 0;        // sin(direction)
};

struct WaveKinematics {
    double elevation = 0;  // surface elevation at (x, y); independent of z
    Vec3 velocity{};       // orbital velocity at the sample point
    Vec3 acceleration{};   // local (Eulerian) acceleration; the convective term
                           // is second order and linear theory drops it
};

class WaveField {
public:
    // Still water.
    WaveField() = default;

    explicit WaveField(const SeaState& sea);
    // Several superposed wave systems -- a wind sea plus swell. Each gets its own
    // phase stream, so adding a swell train does not disturb the wind sea.
    explicit WaveField(const std::vector<SeaState>& seas);

    double elevation(double x, double y, double t) const;
    // Kinematics at a point. z is measured from the still-water plane and is
    // clamped to <= 0 before the exp(k z) decay: above still water the
    // exponential grows without bound, which is the classic way a Froude-Krylov
    // integration explodes in a steep sea. Wheeler stretching is the refinement
    // when nonlinear kinematics arrive.
    WaveKinematics kinematics(const Vec3& position, double t) const;
    Vec3 velocity(const Vec3& position, double t) const;
    Vec3 acceleration(const Vec3& position, double t) const;

    const std::vector<WaveComponent>& components() const { return components_; }
    // Bins of every sea state, appended in the order the sea states were given.
    const std::vector<FrequencyBin>& frequencyBins() const { return frequencyBins_; }
    const std::vector<DirectionBin>& directionBins() const { return directionBins_; }

    // Moments of the *discretised* spectrum: sum omega^order a^2 / 2.
    double spectralMoment(int order) const;
    double zerothMoment() const { return spectralMoment(0); }
    double significantHeight() const;    // 4 sqrt(m0)
    double meanPeriod() const;           // T1 = 2 pi m0 / m1
    double zeroCrossingPeriod() const;   // T2 = 2 pi sqrt(m0 / m2)

private:
    void addSeaState(const SeaState& sea, std::uint64_t stream);

    std::vector<WaveComponent> components_;
    std::vector<FrequencyBin> frequencyBins_;
    std::vector<DirectionBin> directionBins_;
};

}  // namespace sim
