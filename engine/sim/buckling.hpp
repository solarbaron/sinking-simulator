// SPDX-License-Identifier: MIT
// Buckling: the failure a stress magnitude cannot see.
//
// `girderStress()` divides a bending moment by a section modulus and reports how
// close the extreme fibre is to yield. That is the right answer for the fibre in
// *tension*. For the one in compression it is not even the right question: a
// plate panel under edge compression goes unstable and folds at a stress that can
// be a small fraction of yield, and the thinner it is relative to its width, the
// smaller that fraction gets. A ship whose deck reports 40% of yield may already
// have lost the deck.
//
// Two modes matter at the hull-girder scale and they are checked separately
// because they fail at different stresses and are cured by different things:
//
//   * **Plate buckling** between stiffeners. The panel bounded by two
//     longitudinals and two frames buckles as a plate,
//
//         sigma_cr = k pi^2 E / (12 (1 - nu^2)) * (t / b)^2
//
//     with b the *shorter* side. Cured by closer stiffener spacing or thicker
//     plate, and the (t/b)^2 says which is cheaper.
//   * **Stiffener column buckling** between frames. The longitudinal, with its
//     attached strip of plating, buckles as an Euler column,
//
//         sigma_cr = pi^2 E I / (A L^2)
//
//     Cured by a deeper stiffener or closer frames.
//
// **Elastic buckling stress is not strength.** Above about half yield the
// material starts to go plastic before the elastic instability arrives, and the
// elastic formula runs away to values the plate cannot reach -- an infinitely
// thick plate does not have infinite strength. The Johnson-Ostenfeld correction
// caps it, and is continuous at the transition by construction:
//
//     sigma_cr <= sigma_y / 2   ->   unchanged
//     sigma_cr >  sigma_y / 2   ->   sigma_y (1 - sigma_y / (4 sigma_cr))
//
// which is what every classification rule uses and what makes the answer a
// strength rather than an eigenvalue.
//
// Limits, and they are not small. This is elastic buckling of an ideal flat
// panel: no initial distortion, no welding residual stress, no lateral pressure,
// no interaction between the two modes, and simply-supported edges assumed
// throughout. Real panels carry all of those and a real rule check applies
// knock-down factors for them. What is here answers "is this panel in the
// dangerous region", which is the question the hull girder result could not ask
// at all.
#pragma once

#include "girder.hpp"
#include "scantlings.hpp"

#include <vector>

namespace sim {

struct BucklingCheck {
    double elasticStress = 0;    // Pa, the eigenvalue
    double criticalStress = 0;   // Pa, after the plasticity correction
    double appliedStress = 0;    // Pa, compression positive
    double utilisation = 0;      // applied / critical; >= 1 is a buckled panel
    double coefficient = 0;      // k, for a plate; unused for a column
};

// Buckling coefficient for a rectangular plate under uniaxial edge compression,
// simply supported on all four edges:
//
//     k = min over integer m of (m / alpha + alpha / m)^2,   alpha = a / b
//
// where a is along the load and b across it. Exactly 4 whenever alpha is a whole
// number, peaking at 4.5 at the crossovers, and tending to 4 for a long panel --
// which is why "4" is the number quoted for ship plating, and why a *square*
// panel is no weaker than a long one.
double plateBucklingCoefficient(double loadedLength, double width);

// Plate buckling between stiffeners. `width` is the stiffener spacing and
// `loadedLength` the frame spacing, or the other way round; the routine takes the
// shorter as b itself, because a plate does not care which of its sides you call
// which.
BucklingCheck plateBuckling(double thickness, double loadedLength, double width,
                            double appliedCompression, const StructuralMaterial& material);

// Euler column buckling of a stiffener with its attached plating, over a span of
// `length` between frames. `section` is the combined section from
// `stiffenedSection()`.
BucklingCheck columnBuckling(const StiffenedSection& section, double length,
                             double appliedCompression, const StructuralMaterial& material);

// The plasticity cap, exposed because it is worth testing on its own and worth
// reusing: it applies to any elastic instability, not only these two.
double johnsonOstenfeld(double elasticStress, double yieldStrength);

// --- Hull girder scale --------------------------------------------------------

struct GirderBuckling {
    double x = 0;
    double compressiveStress = 0;   // Pa, whichever fibre is in compression
    bool deckInCompression = false; // false means the keel is
    BucklingCheck plate;
    BucklingCheck column;
    double utilisation = 0;         // the worse of the two
};

// Check every station of a bending moment curve against the panels that carry it.
// Stations in pure tension are skipped: there is nothing to buckle.
//
// The compressive fibre follows the sign of the moment. Hogging compresses the
// keel, sagging compresses the deck -- and sagging is the dangerous one on most
// ships, because the deck is the thinner, more widely stiffened structure and it
// is furthest from the neutral axis.
std::vector<GirderBuckling> girderBuckling(const std::vector<GirderStress>& stresses,
                                           const StructuralMesh& structure,
                                           const Scantlings& scantlings);

// Worst utilisation anywhere, and where.
double worstBucklingUtilisation(const std::vector<GirderBuckling>& checks, double* atX = nullptr);

}  // namespace sim
