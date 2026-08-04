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

#include "girder.hpp"
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
//
// It sweeps exactly as far as it is told and nothing checks that the peak is
// inside -- use `collapseCurve()` below unless the range is the point of the
// call, because **sizing the sweep is the part that is easy to get wrong.**
CollapseCurve progressiveCollapse(const std::vector<CollapseElement>& elements,
                                  double maxCurvature, int steps = 200);

// The collapse curve, swept far enough that its peak is genuinely a peak.
// `direction` is the sign of the moment the section has to carry.
//
// Starts at six times `firstYieldCurvature`, exactly as before, and extends only
// when the peak lands on the last point of the sweep -- which is the signature of
// a curve that was still rising. See `extremeFibreYieldCurvature` for the case
// that makes that happen and for what it cost when nothing checked.
CollapseCurve collapseCurve(const std::vector<CollapseElement>& elements, double direction,
                            int steps = 150);

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

// Curvature at which the first element reaches its own limit -- yield, or
// buckling if that comes first. This is *first yield of the section*, and it is
// set by the **weakest** element in it.
double firstYieldCurvature(const std::vector<CollapseElement>& elements);

// Curvature at which the extreme fibre would reach **yield**, ignoring buckling
// entirely. The classical first-yield curvature, and the one no single weak
// element can move: yield strengths across a ship vary by a factor of one and a
// half where a plate's critical stress varies as `t^2`.
//
// It exists because sizing a collapse sweep from `firstYieldCurvature` is only
// safe when no element is anomalously weak, and **that assumption breaks exactly
// when it matters**. A wasted strake -- or a bay a Tier-2 zone has just reported
// as damaged, which is what `promotion.hpp` does -- buckles at a couple of MPa,
// which puts first yield five orders of magnitude below the curvature at which
// the section actually collapses. A sweep sized from it never reaches the peak
// and reports a ship that has lost one bay as having lost nearly all her
// strength. Measured on the reference ferry: thinning forty side panels to an
// eighth took the reported ultimate moment from 1.85e9 to 8.4e7 N m, and taking
// them away *entirely* put it back to 1.7e9 -- a strength that improves when
// material is removed, which is the signature of a truncated sweep and not of
// any physics.
//
// `longitudinalStrength` now sizes from `firstYieldCurvature` as before and
// extends only if the peak lands on the last point of the sweep, so a section
// whose peak was already inside is unaffected to the last bit.
double extremeFibreYieldCurvature(const std::vector<CollapseElement>& elements);

// --- Along the length ---------------------------------------------------------

struct StrengthStation {
    double x = 0;
    double appliedMoment = 0;   // N m, signed, hogging positive
    double ultimateMoment = 0;  // N m, signed, in the direction the moment acts
    double fullyPlastic = 0;    // N m, magnitude
    double margin = 0;          // |ultimate / applied|; zero where nothing is applied
};

// Ultimate strength at every station of a bending moment curve, against the
// moment actually there.
//
// Midship is where the moment usually peaks and where the section is usually
// strongest, and those two facts point in opposite directions -- so the station
// that fails first is not reliably either. This walks both and reports the ratio,
// which is the only quantity that decides anything.
//
// The ultimate is computed in the direction the moment acts at that station,
// because a hull is not equally strong both ways: sagging compresses the deck,
// which is the thin, widely stiffened structure, and on this ferry that costs a
// third of the strength.
std::vector<StrengthStation> longitudinalStrength(const HullGirder& girder,
                                                  const StructuralMesh& structure,
                                                  const Scantlings& scantlings,
                                                  double shedExponent = 0.45,
                                                  int curvatureSteps = 150);

// The smallest margin anywhere, and where. Zero if nothing could be evaluated.
double worstStrengthMargin(const std::vector<StrengthStation>& stations, double* atX = nullptr);

}  // namespace sim
