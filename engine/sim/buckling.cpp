// SPDX-License-Identifier: MIT
#include "buckling.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace sim {
namespace {

// pi^2 E / (12 (1 - nu^2)): the plate flexural constant every buckling formula
// carries. Written once so a stray factor cannot appear in only one of them.
double plateConstant(const StructuralMaterial& m) {
    return kPi * kPi * m.youngsModulus / (12.0 * (1.0 - m.poissonRatio * m.poissonRatio));
}

// Applied over critical, and **a NaN on either side must not come out as zero**.
//
// Both callers wrote this inline as
// `critical > 0 ? std::max(0.0, applied) / critical : 0.0`, which swallows a NaN
// twice over and independently. `std::max(0.0, NaN)` is `0.0 < NaN ? NaN : 0.0`,
// and `0.0 < NaN` is false, so it returns **0.0** -- the argument order decides,
// and `std::max(NaN, 0.0)` would have propagated. Separately, `NaN > 0` is false
// and takes the else, which is also 0.0. So a panel whose applied stress or whose
// capacity is not a number reported **zero utilisation**: the single most
// reassuring answer a collapse check can give, produced by the one input meaning
// the check could not be made at all.
//
// The two early returns in each caller are already written `!(x > 0)`, which is
// the NaN-safe form and sends a NaN to the "nothing computed" return. The idiom
// was in the file. It was the way out that had not been given it.
//
// A zero `critical` is a different thing and keeps its zero: that is the caller
// having been refused, and it is reported by `criticalStress` being zero beside it.
double utilisationOf(double applied, double critical) {
    if (std::isnan(applied) || std::isnan(critical))
        return std::numeric_limits<double>::quiet_NaN();
    if (!(critical > 0)) return 0.0;
    return std::max(0.0, applied) / critical;
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

    // **`b` is the side across the load, and that is not always the shorter one.**
    // This took `min` and `max`, on the argument -- written into `buckling.hpp` in
    // those words -- that "a plate does not care which of its sides you call
    // which". It does. `sigma_cr = k pi^2 D / (b^2 t)` is only the Timoshenko form
    // when `b` is the dimension perpendicular to the load, because that is the
    // dimension the half-waves span; `k` is measured against `a/b` on the same
    // convention. The header said so itself eight lines further down -- "a is
    // along the load and b across it" -- so the file carried both definitions and
    // the code followed the one that is wrong.
    //
    // It matters in one direction only, and it is the bad one. Taking the shorter
    // side always finds a `k` near 4; the true `k` for a panel loaded across its
    // long dimension is larger, but it divides by a `b` that is larger squared,
    // and the square wins. A 0.70 x 2.40 panel loaded along the 0.70:
    //
    //     as written   k = 4.072 on b = 0.70  ->  222.8 MPa
    //     as meant     k = 13.84 on b = 2.40  ->   64.4 MPa
    //
    // a factor of 3.46 unconservative, on a check whose whole purpose is to say
    // that a stress well under yield has already lost the panel.
    //
    // Nothing in the tree reaches it today: every caller passes the frame spacing
    // as the loaded length and the longitudinal spacing as the width -- 2.40
    // against 0.70 on this ship -- so `min` picked the width and got the right
    // answer for the wrong reason. A transversely framed ship, where the frames
    // are the close spacing and the girders the wide one, walks straight into it.
    const double b = width;
    const double a = loadedLength;

    c.coefficient = plateBucklingCoefficient(a, b);
    c.elasticStress = c.coefficient * plateConstant(material) * (thickness / b) * (thickness / b);
    c.criticalStress = johnsonOstenfeld(c.elasticStress, material.yieldStrength);
    c.utilisation = utilisationOf(appliedCompression, c.criticalStress);
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
    c.utilisation = utilisationOf(appliedCompression, c.criticalStress);
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

        // The longitudinal at that same fibre, for the same reason the thickness
        // above comes from the panels present. This passed
        // `scantlings.frameProfile`, the transverse web frame, to a column whose
        // span is the frame spacing -- and `buckling.hpp` calls the column mode
        // "The longitudinal, with its attached strip of plating".
        //
        // A longitudinal runs *along* x, so it is present at this station when the
        // station lies between its ends -- not when its centroid is within half a
        // bay, which is the panels' test and would exclude nearly all of them.
        const StiffenerProfile* longitudinal = nullptr;
        double attached = 0;
        double bestMember = -1e30;
        for (const StructuralMember& m : structure.members) {
            if (m.role != MemberRole::Longitudinal) continue;
            if (s.x < std::min(m.a.x, m.b.x) || s.x > std::max(m.a.x, m.b.x)) continue;
            const double z = 0.5 * (m.a.z + m.b.z);
            const double score = deckCompressed ? z : -z;
            if (score > bestMember) {
                bestMember = score;
                longitudinal = &m.profile;
                attached = m.attachedPlateThickness;
            }
        }
        // No longitudinal at this fibre is not a column of zero strength -- it is a
        // check that cannot be made. `g.column` stays default, so `g.utilisation`
        // below is the plate mode alone rather than a knockdown from a member that
        // is not there.
        if (longitudinal != nullptr) {
            const StiffenedSection combined =
                stiffenedSection(*longitudinal, attached, stiffenerSpacing);
            g.column = columnBuckling(combined, frameSpacing, compression, material);
        }

        g.utilisation = std::max(g.plate.utilisation, g.column.utilisation);
        out.push_back(g);
    }
    return out;
}

double worstBucklingUtilisation(const std::vector<GirderBuckling>& checks, double* atX) {
    double worst = 0;
    // **A station that could not be checked is worse than any station that could.**
    // `c.utilisation > worst` is false for a NaN, so a maximum fold drops exactly
    // the stations whose check failed and returns the worst of the ones that
    // succeeded -- making a NaN one level down invisible one level up, which is
    // how a zero utilisation would have reached a caller even after
    // `utilisationOf` stopped manufacturing one. Once a NaN is seen it stays: no
    // later comparison against it can be true, so the loop stops rather than
    // pretending to keep looking.
    for (const GirderBuckling& c : checks) {
        if (std::isnan(worst)) break;
        if (std::isnan(c.utilisation) || c.utilisation > worst) {
            worst = c.utilisation;
            if (atX) *atX = c.x;
        }
    }
    return worst;
}

}  // namespace sim
