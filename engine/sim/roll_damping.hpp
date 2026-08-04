// SPDX-License-Identifier: MIT
//
// Viscous roll damping by Ikeda's method.
//
// Potential flow gives a ship no roll damping worth the name. Every other mode
// radiates waves strongly and is heavily damped by pressure forces an ideal-flow
// solver reproduces; roll radiates almost nothing, so what stops a ship rolling
// is viscosity -- skin friction, vortices shed from the bilge and from the bilge
// keels, and the lift on a hull that is also going forwards. Ignore it and a ship
// that should roll 8 degrees rolls 35, which is the difference between riding out
// a beam sea and lying on your side in it.
//
// There is no honest first-principles route to those numbers at simulation cost,
// so this is Ikeda's empirical decomposition, the same one every seakeeping code
// uses. Each component is a fit to model experiments, summed, and reported as an
// equivalent *linear* coefficient B44 (N m s/rad) obtained by matching the energy
// dissipated per cycle. Two consequences follow and both matter:
//
//   * B44 is only meaningful at the amplitude and frequency it was asked for.
//     The eddy and bilge-keel moments are physically quadratic in roll rate, so
//     their linear equivalents grow with roll amplitude. Feeding a coefficient
//     computed at 5 degrees into a 25-degree roll underestimates the damping by
//     a factor of several.
//   * The wave (radiation) component is *not* computed here. It is potential
//     flow, it belongs to the boundary-element pipeline (docs/02-simulation.md),
//     and until that exists it is an input: set RollDampingHull::waveDamping to
//     B44W in SI units, or leave it zero and accept an underestimate of roughly
//     5-30% of the total.
//
// Sources, and what was taken from each:
//   * ITTC Recommended Procedure 7.5-02-07-04.5, "Numerical Estimation of Roll
//     Damping" (2011/2017) -- the friction component (Kato 1958, with Tamiya's
//     speed correction), the hull lift component, the bilge-keel normal-force
//     and hull-pressure components, the eddy speed reduction, and the
//     nondimensionalisation.
//   * Kawahara, Maekawa & Ikeda, "A Simple Prediction Formula of Roll Damping of
//     Conventional Cargo Ships on the Basis of Ikeda's Method and Its
//     Limitation", STAB 2009 -- the eddy component, as the regression fitted to
//     Ikeda's own sectional method over a methodical hull series. The full
//     sectional eddy calculation needs Lewis-form conformal mapping coefficients
//     per station and a hull offset table; the regression needs six hull
//     numbers, reproduces the sectional method to within about 10%, and comes
//     with a published range of validity that validateRollDamping() enforces.
//
// Conventions. SI throughout, angles in radians unless a name says otherwise.
// The roll axis height is given as KR, metres above the keel, because that is
// what the rest of the engine knows; Ikeda's OG (still-water level down to the
// roll axis, positive downward) is derived as draft - KR internally.
#pragma once

#include "../core/geometry.hpp"
#include "../core/math.hpp"

#include <string>
#include <vector>

namespace sim {

// Hull form and appendages. These are the parameters Ikeda's method is written
// in: a ship is reduced to its main dimensions plus a midship section shape.
struct RollDampingHull {
    double lengthPp = 0;         // Lpp, m
    double beam = 0;             // B, m, moulded at the waterline
    double draft = 0;            // d, m, mean moulded draft
    double blockCoeff = 0;       // Cb
    double midshipCoeff = 0;     // Cm, midship section area / (B * d)
    double rollAxisAboveKeel = 0;  // KR, m; normally KG

    // Bilge radius, m. Negative means "estimate it from Cm and B/d" -- see
    // bilgeRadiusOrDefault(). Only the bilge-keel component uses it.
    double bilgeRadius = -1;

    // Bilge keels, taken as a symmetric pair. Either dimension zero means the
    // ship has none, and both bilge-keel components are then exactly zero.
    double bilgeKeelLength = 0;   // l_BK, m, streamwise length of one keel
    double bilgeKeelBreadth = 0;  // b_BK, m, projection normal to the shell

    double seaDensity = kRhoSeawater;         // kg/m^3
    double kinematicViscosity = 1.14e-6;      // m^2/s, seawater at 15 C

    // Wave (radiation) roll damping at the frequency of interest, N m s/rad.
    // Not computed here -- see the file header. Passed straight through into the
    // total so callers with a BEM table are not forced to add it themselves.
    double waveDamping = 0;

    double displacementVolume() const;   // m^3, Cb * Lpp * B * d
    double bilgeRadiusOrDefault() const; // m

    // B44 = B44hat * nondimensionalScale(). ITTC (2.4):
    //   B44hat = B44 / (rho * volume * B^2) * sqrt(B / (2 g))
    // so the scale is rho * volume * B^2 / sqrt(B / (2 g)), in kg m^2 / s.
    double nondimensionalScale() const;
};

// The operating point. B44 is a function of all three and quoting it without
// them is meaningless.
struct RollDampingCondition {
    double rollAmplitude = 0;   // phi_a, rad
    double rollFrequency = 0;   // omega, rad/s
    double forwardSpeed = 0;    // U, m/s
};

// Equivalent linear roll damping coefficients, N m s/rad. Components are kept
// separately because which one dominates is the whole diagnostic value of the
// method: friction at model scale, eddy for a bare hull, bilge keels for
// anything with them, lift only when the ship is moving.
struct RollDamping {
    double friction = 0;          // B44F, Kato + Tamiya
    double eddy = 0;              // B44E, bilge vortex shedding
    double lift = 0;              // B44L, zero at zero forward speed
    double bilgeKeelNormal = 0;   // B44BKN, normal force on the keel itself
    double bilgeKeelHull = 0;     // B44BKH, pressure the keel induces on the shell
    double wave = 0;              // B44W, the input coefficient, passed through
    double total = 0;             // sum of the above
    double totalHat = 0;          // nondimensional, ITTC (2.4)

    double bilgeKeel() const { return bilgeKeelNormal + bilgeKeelHull; }
};

RollDamping rollDamping(const RollDampingHull& hull, const RollDampingCondition& condition);

// Derive the hull-form half of a RollDampingHull from an actual hull mesh, so a
// ship asset does not have to carry a second, separately-authored set of main
// dimensions that can drift away from the shape it actually has. Lpp and beam
// are read off the wetted body, draft from the waterline down to the keel, Cb
// from the volume the same integrator gives the hydrostatics, and Cm from the
// largest sectional area found by clipping the hull into slabs -- the same route
// radiationHullFromMesh() takes to its stations, and for the same reason: a hull
// that displaces what it should then also has a block coefficient that says so.
//
// **What is deliberately not derived.** Bilge keels are appendages, and a
// watertight envelope has no appendages in it -- there is nothing in the mesh to
// measure, so they are arguments. The bilge radius is left at its "estimate from
// Cm and B/d" sentinel rather than fitted to the mesh, because Ikeda's bilge-keel
// model idealises the section as a vertical side, a horizontal bottom and a
// quarter-circle bilge, and a radius measured off a shape that is not that is not
// the radius the formulae want. The roll axis is left at zero: it is a property
// of loading rather than of form, it moves as a ship floods, and Ship sets it
// from the live centre of gravity every tick.
RollDampingHull rollDampingHullFromMesh(const TriMesh& hull, double waterlineZ,
                                        double bilgeKeelLength = 0, double bilgeKeelBreadth = 0,
                                        double density = kRhoSeawater);

// Ikeda's method is a fit, and a fit has a domain. This reports, in plain words,
// every way the given hull and operating point sit outside the range over which
// the constituent formulae were validated. An empty result means the answer is
// as trustworthy as the method gets; a non-empty one does not mean the answer is
// garbage, only that nobody checked. It is deliberately advisory -- the caller
// decides, exactly as with Ship::validate().
std::vector<std::string> validateRollDamping(const RollDampingHull& hull,
                                             const RollDampingCondition& condition);

}  // namespace sim
