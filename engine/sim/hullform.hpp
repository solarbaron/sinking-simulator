// SPDX-License-Identifier: MIT
// Hull generation from principal particulars, and measurement back out of them.
//
// Every reference describes a ship the same way: length, beam, draft, and a
// handful of form coefficients. Offsets -- the table of half-breadths that
// actually defines the shape -- are proprietary for almost every real hull and
// published for almost none. So a simulator that can only be given offsets can
// only be given ships nobody has the offsets for.
//
// This generates a hull that *measures* as the requested ship: the same
// displacement, the same block and prismatic coefficients, the same longitudinal
// centre of buoyancy. It is emphatically **not** the real hull. Two ships with
// identical coefficients can have visibly different bodies and measurably
// different seakeeping, so:
//
//   * for stability, flooding and manoeuvring, where the answer is dominated by
//     volume, waterplane and their distribution, this is a reasonable stand-in;
//   * for validation against published RAOs it is **not**, because the
//     comparison would then be against a hull that is not the benchmark ship.
//     Use real offsets for that, and say which you used.
//
// The construction is the classical one. Cp = Cb / Cm splits "how full is the
// midship section" from "how full is the ship along its length". The midship
// section is a rectangle with a radiused bilge, whose radius follows in closed
// form from Cm. The sectional area curve is a two-parameter family solved for Cp
// and LCB together. Every station is then the midship section scaled in breadth
// by the area curve, which makes the sectional areas exact by construction and
// leaves only the tessellation between the request and the result.
#pragma once

#include "../core/geometry.hpp"

#include <string>
#include <vector>

namespace sim {

struct HullParticulars {
    double lengthPp = 100.0;   // L between perpendiculars, m
    double beam = 16.0;        // B, moulded, m
    double draft = 6.0;        // T, design, m
    double depth = 12.0;       // D, moulded, to the weather deck, m

    double blockCoefficient = 0.650;    // Cb = displaced volume / (L B T)
    double midshipCoefficient = 0.980;  // Cm = midship section area / (B T)

    // Longitudinal centre of buoyancy, as a fraction of Lpp forward of midship.
    // Positive is forward. Real ships sit within roughly +-0.03; full hulls carry
    // it forward, fine fast hulls aft.
    double lcbFraction = 0.0;

    // Immersed transom breadth as a fraction of B. Zero gives a cruiser stern
    // that closes to a point; most modern commercial and naval hulls do not.
    double transomFraction = 0.0;

    // A station of exactly zero breadth makes degenerate triangles at the stem,
    // so the ends close to this fraction of B instead. Small enough to be a stem
    // bar, large enough to be a triangle.
    double stemFraction = 0.01;

    // Length of constant midship section, as a fraction of Lpp. Full hulls carry
    // a long one -- a VLCC is nearly half parallel -- and fine fast hulls carry
    // almost none. Zero reproduces the earlier single-curve behaviour exactly.
    double parallelMiddleBodyFraction = 0.0;

    // Measured convergence on the S-175 form: the block-coefficient error is
    // dominated by *waterline* count, because that is what resolves the bilge
    // arc. At 41 stations it runs 0.351% with 11 waterlines, 0.053% with 21 and
    // 0.022% with 41; going from 41 to 161 stations at fixed waterlines moves it
    // only to 0.019%. So the default spends its triangles on waterlines.
    //
    // These are the same nine measurements `docs/05` publishes, and both are now
    // re-derived from `test_hullform.cpp`'s own stdout by the figure gate -- they
    // read 0.43 / 0.13 / 0.06 here and in that table for as long as neither was
    // produced by anything.
    int stationCount = 41;     // odd, so one lands exactly on midship
    int waterlineCount = 21;
};

// Measured back out of a mesh, by the same integrator the hydrostatics use.
// These are what a generated hull should be checked against -- not the numbers
// that were requested, which is the difference between a test and a tautology.
struct HullCoefficients {
    double displacedVolume = 0;   // m^3 at the given draft
    double blockCoefficient = 0;
    double prismaticCoefficient = 0;
    double midshipCoefficient = 0;
    double waterplaneCoefficient = 0;
    double lcbFraction = 0;       // fraction of Lpp forward of midship
    double lengthPp = 0, beam = 0, draft = 0;
};

// The bilge radius that yields the requested midship coefficient, in closed
// form: a rectangle B x T loses r^2 (1 - pi/4) at each bilge, so
// Cm = 1 - 2 r^2 (1 - pi/4) / (B T).
double bilgeRadiusForMidshipCoefficient(double beam, double draft, double midshipCoefficient);

// Half-breadth of the midship section at height z above the baseline.
double midshipHalfBreadth(double beam, double draft, double bilgeRadius, double z);

// The sectional area curve, as a fraction of the midship section area, at
// station position u in [-1, 1] (u = -1 aft perpendicular, +1 forward).
//
// Two exponents, one per end: their mean sets Cp and their difference sets LCB,
// and both relationships are closed forms, so the pair is solved rather than
// searched. See hullform.cpp.
struct AreaCurve {
    double aftExponent = 1.0;
    double forwardExponent = 1.0;
    double transomFraction = 0.0;
    double stemFraction = 0.01;

    // Half-extent of the parallel middle body, as a fraction of Lpp/2. Real
    // ships have one -- a tanker's runs most of its length -- and a curve
    // without it comes out canoe-like: smooth everywhere, tapering from
    // amidships in both directions. It barely moves Cb, Cp or LCB, which is
    // exactly why a coefficient test cannot see its absence.
    double parallelMiddleBody = 0.0;

    double operator()(double u) const;
    double prismaticCoefficient() const;  // analytic, not quadratured
    double lcbFraction() const;           // fraction of Lpp forward of midship
};

// Solve the area-curve exponents for a target Cp and LCB. Both targets are
// achievable over a wide but not unlimited range; `problems` records any target
// that had to be clamped rather than silently missing it.
AreaCurve solveAreaCurve(double prismaticCoefficient, double lcbFraction,
                         double transomFraction = 0.0, double stemFraction = 0.01,
                         double parallelMiddleBody = 0.0,
                         std::vector<std::string>* problems = nullptr);

// Midship coefficient actually delivered by a given bilge radius -- the inverse
// of bilgeRadiusForMidshipCoefficient(). Worth having separately because that
// function clamps: a radius cannot exceed the draft or the half-beam, so a low
// enough Cm is simply unreachable and the caller has to be told rather than
// handed a hull that is fuller than it asked for.
double midshipCoefficientForBilgeRadius(double beam, double draft, double bilgeRadius);

// Generate the hull.
TriMesh makeHullFromParticulars(const HullParticulars& particulars,
                                std::vector<std::string>* problems = nullptr);

// Measure a hull at a given draft. `lengthPp` and `beam` are taken from the mesh
// extents unless supplied as positive values.
HullCoefficients measureHull(const TriMesh& hull, double draft, double lengthPp = 0,
                             double beam = 0);

// Every way these particulars are outside what this construction can honestly
// produce. Advisory, like Ship::validate(). Empty means the request is ordinary.
std::vector<std::string> validateParticulars(const HullParticulars& particulars);

// --- Reference ships ---------------------------------------------------------
//
// Principal particulars from the open literature for hulls that seakeeping and
// manoeuvring work is routinely published against. **These are coefficients, not
// offsets**: the hull that comes out measures like the ship and is not the ship.
// Anything claiming to validate against that ship's published data needs its
// real offset table, and should say so.
HullParticulars kvlcc2Particulars();   // SIMMAN benchmark tanker
HullParticulars s175Particulars();     // ITTC benchmark containership

}  // namespace sim
