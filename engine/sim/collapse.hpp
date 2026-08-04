// SPDX-License-Identifier: MIT
// Progressive collapse: the moment at which the hull girder actually fails.
//
// `girderStress()` reports first yield and `buckling.hpp` reports first
// instability. Neither is the strength of the ship. A hull girder does not fail
// when its worst panel reaches its limit -- that panel sheds load onto its
// neighbours, the neutral axis migrates towards the side still carrying, and the
// section goes on taking moment until enough of it has gone that the total
// starts to fall. The peak of that curve is the ultimate strength, and the
// difference between it and first yield is the margin a ship actually has.
//
// Smith's method, which is what classification societies use:
//
//   1. Cut the section into elements -- strips of plating and stiffeners.
//   2. Impose a curvature. Every element's strain follows from its distance to
//      the neutral axis, and nothing else.
//   3. Each element answers with a stress from its own load-shortening curve:
//      elastic, then capped by yield in tension and by whichever of yield or
//      buckling comes first in compression, then shedding.
//   4. Move the neutral axis until the axial forces sum to zero. They must: no
//      one is pulling on the ends of the ship.
//   5. Sum the moments. Step the curvature and repeat.
//
// The curve that comes out has three landmarks worth naming: the initial slope,
// which is `E` times the section's `sum(A d^2)`; first yield; and the peak.
//
// **The elements' own second moments are deliberately absent from the slope.**
// Smith's method treats each element as carrying axial stress alone -- a strip of
// plating resisting the hull girder by stretching, not by bending about its own
// mid-thickness. That is right, and it means the initial slope is *slightly*
// below `E I` as `hullGirderSection()` reports it. The gap is the elements' own
// inertia, which for ship plating is a fraction of a per cent; `collapseCurve`
// reports both so the difference is visible rather than mysterious.
#pragma once

#include "scantlings.hpp"

#include <vector>

namespace sim {

// How one element answers a strain. Tension is positive throughout.
//
// Elastic to the first cap; capped at yield in tension; capped in compression by
// whichever of yield or buckling arrives first; then shedding, because a buckled
// plate does not hold its critical load, it folds and gives some of it back.
//
// `shedExponent` is a model, and a crude one: `sigma = sigma_c (eps_c/eps)^n`
// past the cap, continuous at the cap by construction. Zero gives a perfect
// plateau -- no shedding, an upper bound on strength. Around 0.4-0.5 is
// representative of stiffened ship plating. The published load-shortening curves
// are per failure mode and considerably more elaborate; this family is chosen so
// the *shape* is right and the assumption is one number a caller can see.
struct LoadShortening {
    double youngsModulus = 206.0e9;
    double yieldStrength = 355.0e6;
    // Critical compressive stress. At or above yield the element yields before it
    // buckles, which is what a heavy bottom structure does.
    double bucklingStress = 1e30;
    double shedExponent = 0.0;

    double stressAt(double strain) const;
    // Compressive capacity actually available: min(yield, buckling).
    double compressiveCapacity() const;
};

// An element of the section, with the curve it answers by.
struct CollapseElement {
    double area = 0;     // m^2
    double height = 0;   // m above the baseline
    LoadShortening curve;
};

struct CollapsePoint {
    double curvature = 0;    // 1/m, positive hogging
    double neutralAxis = 0;  // m above the baseline
    double moment = 0;       // N m, hogging positive
    double residual = 0;     // N, the axial force the solve could not balance
};

struct CollapseCurve {
    std::vector<CollapsePoint> points;

    double ultimateMoment = 0;     // N m, the peak, signed
    double ultimateCurvature = 0;  // 1/m
    double initialStiffness = 0;   // N m^2, the slope at the origin: E * sum(A d^2)
    double elasticNeutralAxis = 0; // m
    // The fully plastic moment of the same section, ignoring buckling entirely:
    // every element at +-yield about the plastic neutral axis. The ceiling that
    // no amount of stiffening can pass, and the number the ultimate moment is
    // worth quoting against.
    double fullyPlasticMoment = 0;
};

// Sweep a curvature range and return the curve. `maxCurvature` is signed: give it
// a negative value for the sagging branch. Steps are uniform.
CollapseCurve progressiveCollapse(const std::vector<CollapseElement>& elements,
                                  double maxCurvature, int steps = 200);

// The neutral axis that balances axial force at a given curvature, and the moment
// there. Exposed because the migration of the neutral axis is the mechanism, and
// a caller may want to watch it.
CollapsePoint collapseAt(const std::vector<CollapseElement>& elements, double curvature);

// Fully plastic moment: every element at +-yield about the plastic neutral axis,
// which is the height that puts equal *area* either side rather than equal first
// moment. Closed form, and the ceiling on `ultimateMoment`.
double fullyPlasticMoment(const std::vector<CollapseElement>& elements);

// Build elements from a structural mesh at a station, giving each the buckling
// stress its own plating and spacing imply. `shedExponent` is passed through.
std::vector<CollapseElement> collapseElementsAt(const StructuralMesh& structure,
                                                const Scantlings& scantlings, double x,
                                                double shedExponent = 0.45);

}  // namespace sim
