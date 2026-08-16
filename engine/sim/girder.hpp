// SPDX-License-Identifier: MIT
// The hull girder: the ship as a beam, and whether she breaks her back.
//
// This is Tier-0 of the structural model in `docs/02-simulation.md` §3 -- the
// cheapest structural answer that is worth having, and the one classical naval
// architecture has been computing by hand for a century. A ship afloat is a free
// beam loaded by two distributions that almost cancel: her weight acting down,
// and buoyancy acting up. Almost, because they cancel in *total* and in *moment*
// -- she floats, and she floats level -- but not station by station. Machinery is
// concentrated aft, cargo amidships, and the ends are fine and carry little
// buoyancy. What is left over bends the hull.
//
//     q(x) = w(x) - b(x)          load per unit length, N/m
//     V(x) = integral of q        shear force, N
//     M(x) = integral of V        bending moment, N m
//
// Weight *minus* buoyancy, and that ordering is the whole of the sign
// convention. With buoyancy minus weight, a hogging ship comes out negative:
// excess buoyancy amidships makes the shear negative over the after half, and
// its integral reaches a minimum at midship. Naval architecture calls hogging
// positive, so the load is written the other way round and dM/dx = V still
// holds. Getting this backwards produces a plausible curve of the right
// magnitude that names every failure as its opposite.
//
// **The free ends are the check.** Nothing holds a floating ship up at the
// perpendiculars, so V and M must both be exactly zero there. That is guaranteed
// by the two equilibrium conditions and by nothing else: the integral of q
// vanishes because she floats, and the integral of x q vanishes because she
// floats level. A residual at the far end is therefore a direct measure of every
// error upstream -- an unbalanced ship, a mis-distributed weight, a quadrature
// that does not close. It is reported rather than hidden.
//
// Sign convention: **hogging positive**. A ship supported amidships -- on a wave
// crest, or simply with fine ends -- arches, putting the deck in tension; that is
// hogging and M > 0. A ship supported at the ends, in a trough, sags and M < 0.
// The two failures look nothing alike and are told apart by that sign.
//
// What this tier cannot do: it is a beam, so it knows nothing about where the
// stress goes within a section, nothing about buckling, nothing about shear lag,
// and nothing about local loads. It answers "is the hull girder overloaded", and
// that is the question that decides whether a ship in a seaway survives being
// long.
#pragma once

#include "scantlings.hpp"
#include "ship.hpp"

#include <string>
#include <vector>

namespace sim {

struct GirderStation {
    double x = 0;              // m, body frame
    double weightPerLength = 0;    // N/m, downward magnitude
    double buoyancyPerLength = 0;  // N/m, upward magnitude
    double loadPerLength = 0;      // N/m, weight - buoyancy
    double shear = 0;              // N
    double moment = 0;             // N m, hogging positive
};

struct HullGirder {
    std::vector<GirderStation> stations;

    double maxShear = 0;          // N, largest magnitude
    double maxShearX = 0;
    double maxMoment = 0;         // N m, largest magnitude, signed
    double maxMomentX = 0;

    // What is left at the after and forward perpendiculars, where both must be
    // zero. Non-zero means the ship was not in equilibrium, or the distributions
    // do not integrate to the same total.
    //
    // **Normalised on `W/2` and `W L / 8`, not on the curve's own peak** -- which
    // is what this comment used to say and what `girder.cpp` rejects by name:
    // "Normalising by the peak of the same quantity is degenerate", because a
    // balanced ship carrying almost no moment has a peak that *is* the residual
    // and reports 1.0. The message strings were corrected when that was found and
    // this declaration was not, which is the half a reader meets first.
    //
    // The difference is not small and it is optimistic. On the suite's own 10%
    // imbalance fixture these print 0.20 and 0.40 while the residual really is
    // 100% of the peak in both -- so a reader taking them as fractions of the peak
    // is out by 5x and 2.5x, in the direction of calling an unbalanced ship
    // balanced. `validateGirder`'s `> 0.05` gate is 5% of `W/2`, on the same terms.
    double shearClosure = 0;
    double momentClosure = 0;

    double totalWeight = 0;     // N
    double totalBuoyancy = 0;   // N

    bool hogging() const { return maxMoment > 0; }
};

// Integrate a load distribution into shear and moment.
//
// Separated from everything that produces a distribution, because this half has
// exact answers and the other half does not. `x` must be ascending. Trapezoidal
// in both integrals, which is exact for a piecewise-linear load and is what makes
// the closed-form checks in tests/test_girder.cpp meaningful.
HullGirder integrateGirder(const std::vector<double>& x,
                           const std::vector<double>& weightPerLength,
                           const std::vector<double>& buoyancyPerLength);

// Buoyancy per unit length along the hull, by slab clipping against the sea
// surface -- the same route `radiationHullFromMesh()` and `measureHull()` take,
// so a hull that displaces what it should also has a buoyancy curve that
// integrates to what it should.
//
// Works against a wavy `Sea`, which is the case that matters: a ship poised on a
// crest amidships is the standard hogging condition, and one in a trough is the
// standard sagging condition.
std::vector<double> buoyancyDistribution(const Ship& ship, const Sea& sea,
                                         const std::vector<double>& x);

// Weight per unit length along the hull.
//
// Lightship weight is distributed by the classical trapezoidal fit: a uniform
// part plus a linear part, with the two chosen so the distribution integrates to
// the lightship mass and to its longitudinal centre of gravity. That is the
// Prohaska/Biles construction, and it is an *assumption* -- a real ship's weight
// curve is stepped, and getting it from the arrangement is what Tier-1 is for.
// Floodwater is not assumed: each compartment's water is spread over that
// compartment's own longitudinal extent, because the simulator knows where it is.
std::vector<double> weightDistribution(const Ship& ship, const std::vector<double>& x);

// Poise the ship on the sea she is given: solve for the sinkage and trim that
// make buoyancy equal weight *and* put the centre of buoyancy under the centre
// of gravity.
//
// This is the "balancing on the wave" step of the classical standard-wave
// calculation, and skipping it does not merely add error -- it destroys the
// result. An unbalanced ship has a net force and a net moment, so the shear and
// bending moment curves grow monotonically towards the far end instead of
// closing, and the wave-induced bending is buried under an imbalance that is not
// physical at all. Measured on a barge under a two-ship-length wave, omitting it
// gave a moment curve that peaked at the forward perpendicular and never returned
// to zero.
//
// Returns false if it cannot balance -- a ship that cannot float on this wave has
// no hull girder answer to give.
bool balanceOnWave(Ship& ship, const Sea& sea, int iterations = 40);

// The whole calculation for a ship as she floats. `balance` poises her on the sea
// first; turn it off only when the caller has already done so.
HullGirder hullGirder(const Ship& ship, const Sea& sea, int stationCount = 41,
                      bool balance = true);

// Stations spread over the hull's length, ends included.
std::vector<double> girderStations(const Ship& ship, int count);

// --- From moment to stress -----------------------------------------------------
//
// The bending moment alone does not say whether a ship survives; a moment is only
// dangerous relative to the section carrying it. Dividing by the section modulus
// is the whole of classical longitudinal strength:
//
//     sigma = M / Z
//
// Hogging puts the deck in tension and the keel in compression, because the hull
// arches and the top fibres stretch. Sagging reverses both. Telling those apart
// matters more than the magnitude does: deck plating in compression buckles at a
// stress well below yield, and that is the failure a stress magnitude alone
// cannot see.
struct GirderStress {
    double x = 0;
    double moment = 0;          // N m, hogging positive
    double modulusDeck = 0;     // m^3
    double modulusKeel = 0;     // m^3
    double stressDeck = 0;      // Pa, tension positive
    double stressKeel = 0;      // Pa, tension positive
    double utilisation = 0;     // worst |stress| over the yield strength
};

// Combine a bending moment curve with the section properties of a structural
// mesh. Stations with no structure at that x are skipped rather than reported as
// infinitely strong.
std::vector<GirderStress> girderStress(const HullGirder& girder, const StructuralMesh& structure,
                                       double yieldStrength);

// Worst utilisation anywhere along the length, and where. Zero if nothing could
// be evaluated.
double worstUtilisation(const std::vector<GirderStress>& stresses, double* atX = nullptr);

// Every way this calculation is outside what a beam model can honestly claim.
std::vector<std::string> validateGirder(const HullGirder& girder);

}  // namespace sim
