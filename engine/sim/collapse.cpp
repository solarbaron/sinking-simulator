// SPDX-License-Identifier: MIT
#include "collapse.hpp"

#include "buckling.hpp"

#include <algorithm>
#include <cmath>

namespace sim {

double LoadShortening::compressiveCapacity() const {
    return std::min(yieldStrength, bucklingStress);
}

double LoadShortening::stressAt(double strain) const {
    if (youngsModulus <= 0) return 0.0;
    if (strain >= 0.0) return std::min(youngsModulus * strain, yieldStrength);

    const double capacity = compressiveCapacity();
    const double criticalStrain = capacity / youngsModulus;   // magnitude
    const double magnitude = -strain;
    if (magnitude <= criticalStrain) return -youngsModulus * magnitude;
    if (shedExponent <= 0.0) return -capacity;
    // Continuous at the cap by construction: at magnitude == criticalStrain the
    // ratio is 1 and the expression is exactly -capacity, whatever the exponent.
    return -capacity * std::pow(criticalStrain / magnitude, shedExponent);
}

CollapsePoint collapseAt(const std::vector<CollapseElement>& elements, double curvature) {
    CollapsePoint p;
    p.curvature = curvature;
    if (elements.empty()) return p;

    double zLo = elements.front().height, zHi = elements.front().height;
    for (const CollapseElement& e : elements) {
        zLo = std::min(zLo, e.height);
        zHi = std::max(zHi, e.height);
    }

    // Net axial force as a function of where the neutral axis sits. Raising the
    // axis lowers the strain everywhere (for positive curvature), so the force is
    // monotone in it and bisection cannot get lost. That monotonicity is why no
    // Newton iteration is needed and why a solve that fails to converge means the
    // section has no capacity at all, rather than a bad initial guess.
    const auto force = [&](double axis) {
        double sum = 0;
        for (const CollapseElement& e : elements)
            sum += e.area * e.curve.stressAt(curvature * (e.height - axis));
        return sum;
    };

    // Widen the bracket generously: at large curvature the balancing axis can sit
    // outside the section entirely.
    const double span = std::max(zHi - zLo, 1e-6);
    double lo = zLo - 4.0 * span, hi = zHi + 4.0 * span;
    if (curvature < 0.0) std::swap(lo, hi);   // force is increasing in the axis instead

    // Bisect to machine precision and stop, rather than to a fixed count. A
    // double carries 53 bits, so a bracket nine section-depths wide is exhausted
    // in about sixty halvings; the 200 this started with spent two thirds of the
    // solve refining digits that do not exist. The sweep below calls this once
    // per curvature step per station, so that waste was the whole cost.
    const double tolerance = 1e-12 * span;
    for (int i = 0; i < 100 && std::abs(hi - lo) > tolerance; ++i) {
        const double mid = 0.5 * (lo + hi);
        (force(mid) > 0.0 ? lo : hi) = mid;
    }
    p.neutralAxis = 0.5 * (lo + hi);
    p.residual = force(p.neutralAxis);

    for (const CollapseElement& e : elements) {
        const double stress = e.curve.stressAt(curvature * (e.height - p.neutralAxis));
        p.moment += e.area * stress * (e.height - p.neutralAxis);
    }
    return p;
}

double fullyPlasticMoment(const std::vector<CollapseElement>& elements) {
    if (elements.empty()) return 0.0;
    // The plastic neutral axis puts equal *area* either side -- equal force, since
    // every element is at the same magnitude of stress -- rather than equal first
    // moment, which is the elastic axis. On a section whose material is uniform
    // the two differ, and using the elastic one here understates the plastic
    // moment.
    double zLo = elements.front().height, zHi = elements.front().height;
    for (const CollapseElement& e : elements) {
        zLo = std::min(zLo, e.height);
        zHi = std::max(zHi, e.height);
    }
    const auto imbalance = [&](double axis) {
        double sum = 0;
        for (const CollapseElement& e : elements)
            sum += e.area * e.curve.yieldStrength * (e.height > axis ? 1.0 : -1.0);
        return sum;
    };
    double lo = zLo, hi = zHi;
    const double plasticTolerance = 1e-12 * std::max(zHi - zLo, 1e-6);
    for (int i = 0; i < 100 && std::abs(hi - lo) > plasticTolerance; ++i) {
        const double mid = 0.5 * (lo + hi);
        (imbalance(mid) > 0.0 ? lo : hi) = mid;
    }
    const double axis = 0.5 * (lo + hi);

    double moment = 0;
    for (const CollapseElement& e : elements)
        moment += e.area * e.curve.yieldStrength * std::abs(e.height - axis);
    return moment;
}

CollapseCurve progressiveCollapse(const std::vector<CollapseElement>& elements,
                                  double maxCurvature, int steps) {
    CollapseCurve curve;
    if (elements.empty() || steps < 2 || maxCurvature == 0.0) return curve;

    double area = 0, moment = 0;
    for (const CollapseElement& e : elements) {
        area += e.area;
        moment += e.area * e.height;
    }
    if (area <= 0) return curve;
    curve.elasticNeutralAxis = moment / area;

    // Smith's method carries axial stress only, so the slope is E sum(A d^2) --
    // the elements' own second moments are absent by construction, not by
    // oversight. See the header.
    for (const CollapseElement& e : elements) {
        const double d = e.height - curve.elasticNeutralAxis;
        curve.initialStiffness += e.curve.youngsModulus * e.area * d * d;
    }
    curve.fullyPlasticMoment = fullyPlasticMoment(elements);

    curve.points.reserve(static_cast<std::size_t>(steps));
    for (int i = 0; i < steps; ++i) {
        const double k = maxCurvature * (i + 1) / steps;
        const CollapsePoint p = collapseAt(elements, k);
        curve.points.push_back(p);
        if (std::abs(p.moment) > std::abs(curve.ultimateMoment)) {
            curve.ultimateMoment = p.moment;
            curve.ultimateCurvature = p.curvature;
        }
    }
    return curve;
}

double firstYieldCurvature(const std::vector<CollapseElement>& elements) {
    if (elements.empty()) return 0.0;
    double area = 0, moment = 0;
    for (const CollapseElement& e : elements) {
        area += e.area;
        moment += e.area * e.height;
    }
    if (area <= 0) return 0.0;
    const double axis = moment / area;

    // The fibre that yields first is the one with the largest strain demand per
    // unit curvature relative to its own capacity -- which is not simply the
    // furthest one, because a panel that buckles at half yield gets there in half
    // the strain.
    double smallest = 0;
    for (const CollapseElement& e : elements) {
        const double lever = std::abs(e.height - axis);
        if (lever <= 0 || e.curve.youngsModulus <= 0) continue;
        const double limit = std::min(e.curve.yieldStrength, e.curve.compressiveCapacity());
        const double k = limit / (e.curve.youngsModulus * lever);
        if (smallest == 0.0 || k < smallest) smallest = k;
    }
    return smallest;
}

std::vector<StrengthStation> longitudinalStrength(const HullGirder& girder,
                                                  const StructuralMesh& structure,
                                                  const Scantlings& scantlings,
                                                  double shedExponent, int curvatureSteps) {
    std::vector<StrengthStation> out;
    out.reserve(girder.stations.size());

    for (const GirderStation& g : girder.stations) {
        const std::vector<CollapseElement> elements =
            collapseElementsAt(structure, scantlings, g.x, shedExponent);
        if (elements.empty()) continue;

        StrengthStation s;
        s.x = g.x;
        s.appliedMoment = g.moment;

        // Sweep far enough to be certain the peak is inside. Collapse arrives at
        // a small multiple of first yield; six is comfortable and still cheap.
        const double yieldCurvature = firstYieldCurvature(elements);
        if (!(yieldCurvature > 0)) continue;
        // In the direction the moment actually acts. A hull is not equally strong
        // both ways and evaluating the wrong branch would report the stronger one.
        const double sign = g.moment >= 0 ? 1.0 : -1.0;
        const CollapseCurve curve =
            progressiveCollapse(elements, sign * 6.0 * yieldCurvature, curvatureSteps);

        s.ultimateMoment = curve.ultimateMoment;
        s.fullyPlastic = curve.fullyPlasticMoment;
        s.margin = std::abs(g.moment) > 1.0
                       ? std::abs(curve.ultimateMoment / g.moment)
                       : 0.0;
        out.push_back(s);
    }
    return out;
}

double worstStrengthMargin(const std::vector<StrengthStation>& stations, double* atX) {
    double worst = 0;
    for (const StrengthStation& s : stations) {
        if (s.margin <= 0) continue;           // nothing applied here to be worst at
        if (worst == 0.0 || s.margin < worst) {
            worst = s.margin;
            if (atX) *atX = s.x;
        }
    }
    return worst;
}

std::vector<CollapseElement> collapseElementsAt(const StructuralMesh& structure,
                                                const Scantlings& scantlings, double x,
                                                double shedExponent) {
    std::vector<CollapseElement> out;
    const double frameSpacing =
        structure.frameSpacing > 0 ? structure.frameSpacing : scantlings.frameSpacing;
    const double stiffenerSpacing = scantlings.longitudinalSpacing;

    for (const SectionElement& e : sectionElements(structure, x)) {
        if (e.area <= 0) continue;
        const StructuralMaterial material =
            e.material >= 0 && e.material < static_cast<int>(structure.materials.size())
                ? structure.materials[static_cast<std::size_t>(e.material)]
                : ah36Steel();

        CollapseElement c;
        c.area = e.area;
        c.height = e.height;
        c.curve.youngsModulus = material.youngsModulus;
        c.curve.yieldStrength = material.yieldStrength;
        c.curve.shedExponent = shedExponent;

        if (e.stiffener) {
            // A stiffener is a column between frames, and it is stocky enough that
            // it usually yields first. It does not shed the way a plate panel
            // does, so it holds its capacity.
            const StiffenedSection combined =
                stiffenedSection(scantlings.frameProfile, e.thickness, stiffenerSpacing);
            c.curve.bucklingStress =
                columnBuckling(combined, frameSpacing, 0.0, material).criticalStress;
            c.curve.shedExponent = 0.0;
        } else {
            c.curve.bucklingStress =
                plateBuckling(e.thickness, frameSpacing, stiffenerSpacing, 0.0, material)
                    .criticalStress;
        }
        out.push_back(c);
    }
    return out;
}

}  // namespace sim
