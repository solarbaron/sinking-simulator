// SPDX-License-Identifier: MIT
#include "buckling.hpp"

#include <algorithm>
#include <cmath>

namespace sim {
namespace {

// pi^2 E / (12 (1 - nu^2)): the plate flexural constant every buckling formula
// carries. Written once so a stray factor cannot appear in only one of them.
double plateConstant(const StructuralMaterial& m) {
    return kPi * kPi * m.youngsModulus / (12.0 * (1.0 - m.poissonRatio * m.poissonRatio));
}

}  // namespace

double plateBucklingCoefficient(double loadedLength, double width) {
    if (!(loadedLength > 0) || !(width > 0)) return 0.0;
    const double alpha = loadedLength / width;
    // The plate buckles into whichever number of half-waves is cheapest, so the
    // coefficient is a minimum over m and not a formula in alpha. Searching a few
    // m either side of alpha is exact: (m/a + a/m)^2 is convex in m.
    const int centre = std::max(1, static_cast<int>(std::round(alpha)));
    double best = 0;
    for (int m = std::max(1, centre - 2); m <= centre + 2; ++m) {
        const double k = (m / alpha + alpha / m) * (m / alpha + alpha / m);
        if (best == 0.0 || k < best) best = k;
    }
    return best;
}

double johnsonOstenfeld(double elasticStress, double yieldStrength) {
    if (!(yieldStrength > 0)) return 0.0;
    if (elasticStress <= 0.5 * yieldStrength) return elasticStress;
    return yieldStrength * (1.0 - yieldStrength / (4.0 * elasticStress));
}

BucklingCheck plateBuckling(double thickness, double loadedLength, double width,
                            double appliedCompression, const StructuralMaterial& material) {
    BucklingCheck c;
    c.appliedStress = appliedCompression;
    if (!(thickness > 0) || !(loadedLength > 0) || !(width > 0)) return c;

    // A plate does not know which of its sides the caller called which: b is the
    // short one and a the long one, whatever they were named.
    const double b = std::min(loadedLength, width);
    const double a = std::max(loadedLength, width);

    c.coefficient = plateBucklingCoefficient(a, b);
    c.elasticStress = c.coefficient * plateConstant(material) * (thickness / b) * (thickness / b);
    c.criticalStress = johnsonOstenfeld(c.elasticStress, material.yieldStrength);
    c.utilisation = c.criticalStress > 0 ? std::max(0.0, appliedCompression) / c.criticalStress : 0.0;
    return c;
}

BucklingCheck columnBuckling(const StiffenedSection& section, double length,
                             double appliedCompression, const StructuralMaterial& material) {
    BucklingCheck c;
    c.appliedStress = appliedCompression;
    if (!(section.area > 0) || !(section.secondMoment > 0) || !(length > 0)) return c;

    c.elasticStress = kPi * kPi * material.youngsModulus * section.secondMoment /
                      (section.area * length * length);
    c.criticalStress = johnsonOstenfeld(c.elasticStress, material.yieldStrength);
    c.utilisation = c.criticalStress > 0 ? std::max(0.0, appliedCompression) / c.criticalStress : 0.0;
    return c;
}

std::vector<GirderBuckling> girderBuckling(const std::vector<GirderStress>& stresses,
                                           const StructuralMesh& structure,
                                           const Scantlings& scantlings,
                                           double yieldStrength) {
    std::vector<GirderBuckling> out;
    out.reserve(stresses.size());

    const double frameSpacing = structure.frameSpacing > 0 ? structure.frameSpacing
                                                           : scantlings.frameSpacing;
    const double stiffenerSpacing = scantlings.longitudinalSpacing;
    StructuralMaterial material =
        structure.materials.empty() ? ah36Steel() : structure.materials.front();
    // An explicit yield strength overrides the plating's own, and only the yield
    // strength -- the modulus and Poisson ratio still come from the steel, because
    // what a caller is choosing here is which of a mixed ship's yields to be a
    // ratio to, not a different material.
    if (yieldStrength > 0) material.yieldStrength = yieldStrength;

    for (const GirderStress& s : stresses) {
        // Only the compressed fibre can buckle. Hogging compresses the keel,
        // sagging the deck; a station with no moment has neither.
        const bool deckCompressed = s.stressDeck < 0;
        const double compression =
            std::max(deckCompressed ? -s.stressDeck : -s.stressKeel, 0.0);
        if (compression <= 0) continue;

        // The plating that carries it: the thickness at the compressed fibre,
        // taken from the panels actually present at this station rather than from
        // the scantling description, so a locally thickened strake counts.
        double thickness = 0;
        double best = -1e30;
        for (const PlatePanel& p : structure.panels) {
            const Vec3 c = p.centroid();
            if (std::abs(c.x - s.x) > 0.5 * frameSpacing) continue;
            // Deck fibre is the highest panel present, keel fibre the lowest.
            const double score = deckCompressed ? c.z : -c.z;
            if (score > best) {
                best = score;
                thickness = p.thickness;
            }
        }
        if (!(thickness > 0)) continue;

        GirderBuckling g;
        g.x = s.x;
        g.compressiveStress = compression;
        g.deckInCompression = deckCompressed;
        g.plate = plateBuckling(thickness, frameSpacing, stiffenerSpacing, compression, material);

        const StiffenedSection combined =
            stiffenedSection(scantlings.frameProfile, thickness, stiffenerSpacing);
        g.column = columnBuckling(combined, frameSpacing, compression, material);

        g.utilisation = std::max(g.plate.utilisation, g.column.utilisation);
        out.push_back(g);
    }
    return out;
}

double worstBucklingUtilisation(const std::vector<GirderBuckling>& checks, double* atX) {
    double worst = 0;
    for (const GirderBuckling& c : checks)
        if (c.utilisation > worst) {
            worst = c.utilisation;
            if (atX) *atX = c.x;
        }
    return worst;
}

}  // namespace sim
