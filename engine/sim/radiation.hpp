// SPDX-License-Identifier: MIT
//
// Radiation hydrodynamics: the part of the fluid force that depends on the
// *history* of the ship's motion.
//
// A ship pushing water aside carries some of it along and radiates the rest away
// as waves. The carried water is added mass and acts instantaneously; the
// radiated waves travel out, and while they are still nearby they push back. So
// the force at time t depends on where the ship has been, not only on where it is
// now. Cummins (1962) writes that as
//
//     (M + A_inf) xddot + integral_0^t K(t - tau) xdot(tau) dtau + C x = F_ext
//
// where A_inf is the infinite-frequency added mass and K is the retardation
// (impulse-response) function. Ogilvie (1964) ties K to the frequency-domain
// coefficients any seakeeping method produces:
//
//     K(t)     = (2/pi) integral_0^inf B(omega) cos(omega t) domega
//     B(omega) =        integral_0^inf K(t)     cos(omega t) dt
//     A(omega) = A_inf - (1/omega) integral_0^inf K(t) sin(omega t) dt
//
// Those are exact identities, not approximations, which makes them the strongest
// check available on everything upstream of them: A(omega) and B(omega) are a
// Kramers-Kronig pair, so a frequency-domain table that is internally
// inconsistent cannot satisfy them. tests/test_radiation.cpp exploits that.
//
// **What this module is, and is not.** docs/02-simulation.md plans an offline
// boundary-element solve (NEMOH/Capytaine) over the whole hull, with the
// coefficient tables shipped as ship assets. That solver is not here and this is
// not a substitute for it. What is here is:
//
//   * **Strip theory.** The hull is cut into stations; each station's
//     two-dimensional radiation problem is solved; the results are integrated
//     along the length. Valid for slender bodies at moderate frequencies, and
//     silent about everything three-dimensional -- see the validity list at the
//     bottom of this header and in docs/02-simulation.md.
//   * **An exact 2D solver per station.** Not a chart lookup and not a
//     regression: a close-fit source distribution over the section contour with
//     the deep-water free-surface Green's function, in the manner of Frank
//     (1967). The Green's function is evaluated in closed form through the
//     complex exponential integral rather than by quadrature.
//   * **The Cummins/Ogilvie machinery and a state-space approximation of the
//     convolution**, which is what makes the memory affordable at 100 Hz.
//
// Section shape comes from Lewis-form conformal mapping, so a station is
// described by three numbers -- beam, draft, sectional area coefficient -- which
// is what a hull table gives you.
//
// Conventions. SI throughout. Body frame as elsewhere in the engine: +x forward,
// +y to port, +z up, origin at midship on the baseline. Degrees of freedom are
// indexed 0..5 = surge, sway, heave, roll, pitch, yaw. The 2D section solver
// works in its own plane with x transverse and y vertically *up* from the still
// waterline, because that is the frame the free-surface condition is written in;
// SectionCoefficients therefore reports about the waterline centreline and
// stripTheoryTable() does the transfer to the baseline origin.
#pragma once

#include "../core/math.hpp"

#include <array>
#include <string>
#include <vector>

namespace sim {

// --- Section geometry --------------------------------------------------------

// A Lewis form: the image of the unit circle under
//
//     z = scale * (zeta + a1/zeta + a3/zeta^3)
//
// which, on the circle, is the section
//
//     halfBreadth(theta) =  scale * ((1 + a1) sin theta - a3 sin 3 theta)
//     depth(theta)       =  scale * ((1 - a1) cos theta + a3 cos 3 theta)
//
// with theta = 0 at the keel and pi/2 at the waterline. Three numbers -- beam,
// draft and sectional area coefficient -- determine a1, a3 and the scale
// exactly, which is why Lewis forms are the standard idealisation in strip
// theory: a hull offset table is usually not available and a station's area
// almost always is.
struct LewisSection {
    double beam = 0;              // B, m, at the waterline
    double draft = 0;             // T, m
    double areaCoefficient = 0;   // sigma = sectional area / (B T)

    double scale = 0;             // the mapping's scale factor
    double a1 = 0;
    double a3 = 0;

    double sectionArea() const;                     // m^2, from the mapping
    Vec3 point(double theta) const;                 // (transverse, vertical, 0), vertical <= 0
    Vec3 tangent(double theta) const;               // d(point)/d(theta), unnormalised
};

// Solve the Lewis mapping for the given beam, draft and area coefficient.
// Returns a section whose analytic area matches sigma * B * T to machine
// precision when sigma is attainable, and the closest attainable form when it is
// not; validateLewisSection() reports the latter case rather than hiding it.
LewisSection lewisSection(double beam, double draft, double areaCoefficient);

// Every way this section sits outside what a Lewis form can represent. Empty
// means the three-parameter family reproduces the requested section exactly.
std::vector<std::string> validateLewisSection(const LewisSection& section);

// --- Two-dimensional radiation coefficients ----------------------------------

// The three modes a 2D section has. Indices into SectionCoefficients.
enum SectionMode { kSectionSway = 0, kSectionHeave = 1, kSectionRoll = 2 };

// Sectional coefficients, per unit length: added mass in kg/m (sway, heave),
// kg m^2/m for roll, kg m/m for the sway-roll coupling; damping the same
// divided by seconds. Taken about the *waterline centreline*, which is where the
// 2D problem is naturally posed.
struct SectionCoefficients {
    double addedMass[3][3]{};
    double damping[3][3]{};

    // |B_near / B_far - 1| for the heave mode: the damping computed by
    // integrating pressure over the hull, against the same damping computed from
    // the amplitude of the wave radiated to infinity. Those are the same number
    // by conservation of energy, computed along completely different paths, so
    // this is a direct measure of how well the discretisation is doing. It is
    // reported rather than asserted because the caller chose the panel count.
    //
    // **Check it.** It reads under 0.01 in the useful band and runs to 60 at an
    // irregular frequency, where the source formulation returns nonsense --
    // including negative damping. stripTheoryTable() uses it to reject and
    // interpolate over those; a caller going straight to sectionCoefficients()
    // has to do the same. See the note in radiation.cpp and
    // docs/02-simulation.md.
    double energyResidual = 0;
};

// One 2D radiation solve. `panels` is the number of panels on *each* half of the
// contour; port/starboard symmetry is exploited exactly, so the sway and roll
// modes are solved with antisymmetric source pairs and heave with symmetric
// ones, and a section symmetric about its centreline can never produce an
// asymmetric answer through rounding.
//
// Cost is O(panels^3) for the factorisation and O(panels^2) for the influence
// coefficients, the latter dominating at the panel counts that are useful here.
SectionCoefficients sectionCoefficients(const LewisSection& section, double omega,
                                        double density, int panels = 48);

// The omega -> infinity limit, where the free-surface condition degenerates to
// phi = 0 and the problem becomes a rigid-lid one with a closed-form Green's
// function. Solved exactly rather than by taking a large omega, which the panel
// method would resolve badly. Damping is identically zero in this limit.
SectionCoefficients sectionCoefficientsInfiniteFrequency(const LewisSection& section,
                                                         double density, int panels = 48);

// --- Hull --------------------------------------------------------------------

struct RadiationStation {
    double x = 0;                 // m, body frame, +forward from midship
    double beam = 0;              // m
    double draft = 0;             // m
    double areaCoefficient = 0;   // sigma
};

// Stations must be sorted by x and there must be at least two. `draft` is the
// still-water draft, and station drafts are measured from the same waterline;
// the roll axis for the assembled matrices is the body-frame origin, i.e. on the
// baseline at midship, not at the waterline.
struct RadiationHull {
    std::vector<RadiationStation> stations;
    double draft = 0;                   // m, still-water draft (baseline to waterline)
    double density = kRhoSeawater;      // kg/m^3
    int panelsPerHalfSection = 48;
};

// --- Frequency-domain table --------------------------------------------------

using Matrix6 = std::array<std::array<double, 6>, 6>;

// A(omega) and B(omega) on a frequency grid, plus the infinite-frequency added
// mass. This is the object an offline BEM run would produce, and the Cummins
// machinery below consumes nothing else -- so replacing strip theory with a real
// panel code later changes only how this table is filled in.
struct RadiationTable {
    std::vector<double> omega;        // rad/s, strictly increasing
    std::vector<Matrix6> addedMass;   // A(omega_k)
    std::vector<Matrix6> damping;     // B(omega_k)
    Matrix6 addedMassInfinite{};      // A_inf

    // Worst per-station energy residual over the frequencies that were *kept* --
    // see SectionCoefficients::energyResidual.
    double worstEnergyResidual = 0;
    // How many (station, frequency) solves were rejected as irregular
    // frequencies and interpolated over. Non-zero is normal and small; a large
    // fraction of the grid means the table should not be trusted.
    int repairedSolves = 0;
    int totalSolves = 0;

    int size() const { return static_cast<int>(omega.size()); }
    // B_ij(omega), linearly interpolated in omega and zero outside the grid.
    // Zero outside is the physically right extrapolation: radiation damping
    // vanishes at both ends of the spectrum.
    double dampingAt(int i, int j, double omegaQuery) const;
};

// Strip theory. Each station's 2D problem is solved at every frequency and the
// sectional coefficients are integrated along the length by the trapezium rule
// over the station positions:
//
//   A22 = int a22 dx        A33 = int a33 dx        A44 = int a44 dx
//   A24 = int a24 dx        A26 = int x a22 dx      A35 = -int x a33 dx
//   A46 = int x a24 dx      A55 = int x^2 a33 dx    A66 = int x^2 a22 dx
//
// and symmetrically. Surge is left at zero: a strip has no longitudinal
// radiation problem at all, and surge added mass is entirely a three-dimensional
// end effect.
RadiationTable stripTheoryTable(const RadiationHull& hull, const std::vector<double>& omega);

// A frequency grid that resolves a hull of this length: geometric spacing
// between the given bounds, which puts points where B(omega) has structure
// instead of spending them on a flat tail.
std::vector<double> radiationFrequencyGrid(double omegaMin, double omegaMax, int count);

// --- Cummins / Ogilvie -------------------------------------------------------

// The retardation function K_ij(t) = (2/pi) int_0^inf B_ij(omega) cos(omega t)
// domega, sampled at t = 0, dt, 2 dt, ... (count samples). The omega integral is
// taken over the table's grid by Simpson's rule on a refined mesh, so the
// result is limited by the table's own resolution and by where it was truncated,
// not by the quadrature.
std::vector<double> retardationFunction(const RadiationTable& table, int i, int j, double dt,
                                        int count);

// The inverse leg of the Ogilvie pair: B(omega) = int_0^inf K(t) cos(omega t) dt,
// by Simpson over the sampled K. Round-tripping B -> K -> B is an identity and
// the test suite requires it to hold.
double dampingFromRetardation(const std::vector<double>& k, double dt, double omega);

// A(omega) = A_inf - (1/omega) int_0^inf K(t) sin(omega t) dt. Equivalently, and
// this is how the test suite uses it, A_inf recovered from a *known* A(omega) is
// independent of omega -- which couples the added mass and the damping and so
// checks both at once.
double addedMassFromRetardation(const std::vector<double>& k, double dt, double omega,
                                double addedMassInfinite);
double infiniteAddedMassFromRetardation(const std::vector<double>& k, double dt, double omega,
                                        double addedMassAtOmega);

// The time by which |K(t)| has fallen to `fraction` of its peak and stays there.
// This is the ship's memory: convolutions can be truncated here, and the
// state-space model must reproduce at least this much of K.
double memoryDecayTime(const std::vector<double>& k, double dt, double fraction = 0.01);

// --- State-space approximation of the convolution ----------------------------

// mu(t) = int_0^t K(t - tau) v(tau) dtau is replaced by
//
//     xdot = A x + B v,     mu = C x
//
// with C exp(A t) B fitted to K(t). A is block diagonal: one 2x2 block per
// complex pole pair and one 1x1 block per real pole, so `order` is the number of
// states and is even when every pole is complex. Poles come from a
// least-squares linear-prediction (Prony) fit to the sampled K, which is exact
// for a K that really is a sum of damped sinusoids -- the test suite checks that
// against a synthetic one before trusting it on a real hull.
struct RadiationStateSpace {
    struct Mode {
        double decay = 0;      // sigma, s^-1; the pole is -sigma +- i frequency
        double frequency = 0;  // rad/s; zero for a real pole
        double c0 = 0;         // output row
        double c1 = 0;         // second output component; unused for a real pole
    };
    std::vector<Mode> modes;

    int stateCount() const;
    double impulseResponse(double t) const;   // C exp(A t) B
    bool stable() const;                      // every pole strictly in the left half-plane
};

struct StateSpaceFit {
    RadiationStateSpace model;
    double rmsError = 0;        // over the samples used for the fit
    double peakError = 0;
    double relativeRms = 0;     // rmsError / rms(K)
    bool converged = false;
};

// Fit `order` states to the sampled retardation function. Orders of 4-8 are what
// the literature recommends and what this was tuned against; the fit quality is
// returned rather than assumed.
StateSpaceFit fitStateSpace(const std::vector<double>& k, double dt, int order);

// --- Runtime evaluator -------------------------------------------------------

// What a ship calls each tick. Holds one state-space model per (i, j) entry with
// non-negligible damping, advances them with the ship's velocity, and returns
// the memory force. The infinite-frequency added mass is *not* applied here: it
// belongs on the left-hand side of the equation of motion, added to the rigid
// body's mass matrix, and is exposed for that purpose.
class RadiationForce {
public:
    RadiationForce() = default;
    // Builds state-space models for every entry whose damping is not negligible
    // against the largest entry in the table. `order` states per entry.
    RadiationForce(const RadiationTable& table, double dt, int samples, int order);

    // Advance the memory states by dt with the body-frame velocity vector
    // (surge, sway, heave, roll rate, pitch rate, yaw rate). The update is the
    // exact zero-order-hold discretisation of the block-diagonal system, so it
    // is unconditionally stable and does not care whether dt matches the dt the
    // model was fitted at.
    void step(const std::array<double, 6>& velocity, double dt);
    void reset();

    // mu_i = sum_j int K_ij(t - tau) v_j(tau) dtau, the radiation memory force
    // (N for i = 0..2, N m for i = 3..5). Subtract it from the applied force.
    std::array<double, 6> memoryForce() const;

    const Matrix6& addedMassInfinite() const { return addedMassInfinite_; }
    // Fit quality of the entry with the largest damping, for reporting.
    double worstRelativeRms() const { return worstRelativeRms_; }
    int modelCount() const { return static_cast<int>(entries_.size()); }

private:
    struct Entry {
        int i = 0, j = 0;
        RadiationStateSpace model;
        std::vector<double> state;   // 2 per complex mode, 1 per real mode
    };
    std::vector<Entry> entries_;
    Matrix6 addedMassInfinite_{};
    double worstRelativeRms_ = 0;
};

// Every way this hull and frequency grid sit outside what strip theory and this
// implementation can honestly claim. Advisory, like Ship::validate().
std::vector<std::string> validateRadiationHull(const RadiationHull& hull,
                                               const std::vector<double>& omega);

}  // namespace sim
