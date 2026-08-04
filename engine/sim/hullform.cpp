// SPDX-License-Identifier: MIT
#include "hullform.hpp"

#include <algorithm>
#include <cmath>

namespace sim {
namespace {

// A rectangle loses this fraction of r^2 at each radiused bilge: the corner
// square minus the quarter circle.
const double kBilgeLoss = 1.0 - kPi / 4.0;

// Bisection on a monotone function. Used twice below, and used rather than
// Newton because both targets are monotone but their derivatives are not worth
// writing down, and a solve that cannot diverge is worth more here than one that
// converges in fewer steps.
double bisect(double lo, double hi, double target, const auto& f, int iterations = 80) {
    const bool rising = f(hi) > f(lo);
    for (int i = 0; i < iterations; ++i) {
        const double mid = 0.5 * (lo + hi);
        const bool high = rising ? (f(mid) > target) : (f(mid) < target);
        (high ? hi : lo) = mid;
    }
    return 0.5 * (lo + hi);
}

}  // namespace

double bilgeRadiusForMidshipCoefficient(double beam, double draft, double midshipCoefficient) {
    if (beam <= 0 || draft <= 0) return 0.0;
    const double lost = (1.0 - std::clamp(midshipCoefficient, 0.0, 1.0)) * beam * draft;
    const double radius = std::sqrt(std::max(0.0, lost / (2.0 * kBilgeLoss)));
    // A bilge radius cannot exceed the half-beam or the draft without the
    // section folding through itself.
    return std::min(radius, std::min(0.5 * beam, draft));
}

double midshipHalfBreadth(double beam, double draft, double bilgeRadius, double z) {
    const double half = 0.5 * beam;
    if (z <= 0.0) return std::max(0.0, half - bilgeRadius);
    if (z >= bilgeRadius || bilgeRadius <= 0.0) return half;
    (void)draft;
    const double dz = bilgeRadius - z;
    return half - bilgeRadius + std::sqrt(std::max(0.0, bilgeRadius * bilgeRadius - dz * dz));
}

// --- Sectional area curve ----------------------------------------------------
//
// f(u) = 1 - (1 - e) |u|^n, one exponent and one end value per end. At u = 0 the
// section is the midship section; at |u| = 1 it is the end value.

namespace {

// One half of the curve, integrated over [0, 1]: p of parallel middle body at
// f = 1, then a taper to the end value e over the remaining (1 - p).
//
// With t = (u - p) / (1 - p) the taper is the same 1 - (1-e) t^n as before, so
//   area   = p + (1-p) [1 - (1-e)/(n+1)]
//   moment = p^2/2 + (1-p) { p [1 - (1-e)/(n+1)] + (1-p) [1/2 - (1-e)/(n+2)] }
// both exact, which keeps the solve analytic rather than quadratured.
struct HalfIntegrals {
    double area = 0;
    double moment = 0;
};

HalfIntegrals halfIntegrals(double n, double end, double p) {
    const double body = 1.0 - end;
    const double taper = 1.0 - body / (n + 1.0);
    HalfIntegrals h;
    h.area = p + (1.0 - p) * taper;
    h.moment = 0.5 * p * p +
               (1.0 - p) * (p * taper + (1.0 - p) * (0.5 - body / (n + 2.0)));
    return h;
}

}  // namespace

double AreaCurve::operator()(double u) const {
    const double clamped = std::clamp(u, -1.0, 1.0);
    const double magnitude = std::abs(clamped);
    const double p = std::clamp(parallelMiddleBody, 0.0, 0.98);
    if (magnitude <= p) return 1.0;   // constant midship section
    const double end = clamped < 0.0 ? transomFraction : stemFraction;
    const double n = clamped < 0.0 ? aftExponent : forwardExponent;
    const double t = (magnitude - p) / (1.0 - p);
    return 1.0 - (1.0 - end) * std::pow(t, n);
}

double AreaCurve::prismaticCoefficient() const {
    const double p = std::clamp(parallelMiddleBody, 0.0, 0.98);
    return 0.5 * (halfIntegrals(aftExponent, transomFraction, p).area +
                  halfIntegrals(forwardExponent, stemFraction, p).area);
}

double AreaCurve::lcbFraction() const {
    // The aft half contributes with the opposite sign because u is negative there.
    const double p = std::clamp(parallelMiddleBody, 0.0, 0.98);
    const double moment = halfIntegrals(aftExponent, transomFraction, p).moment -
                          halfIntegrals(forwardExponent, stemFraction, p).moment;
    const double area = 2.0 * prismaticCoefficient();
    if (area <= 0.0) return 0.0;
    // moment/area is in units of u, i.e. of half-lengths; halve again for Lpp.
    return -moment / area * 0.5;
}

double midshipCoefficientForBilgeRadius(double beam, double draft, double bilgeRadius) {
    if (beam <= 0 || draft <= 0) return 0.0;
    return 1.0 - 2.0 * bilgeRadius * bilgeRadius * kBilgeLoss / (beam * draft);
}

AreaCurve solveAreaCurve(double prismaticCoefficient, double lcbFraction,
                         double transomFraction, double stemFraction,
                         double parallelMiddleBody, std::vector<std::string>* problems) {
    AreaCurve curve;
    curve.parallelMiddleBody = std::clamp(parallelMiddleBody, 0.0, 0.98);
    // The end values shape both Cp and LCB, so they go in *before* the solve.
    // Setting them afterwards and rescaling the exponents to recover Cp -- which
    // is what the first version did -- leaves LCB wherever the rescale put it.
    // The error tracked transomFraction exactly: 0.0003 of Lpp at a cruiser
    // stern, 0.021 at a wide transom, which is a metre and a half on a frigate.
    curve.transomFraction = std::clamp(transomFraction, 0.0, 1.0);
    curve.stemFraction = std::clamp(stemFraction, 1e-4, 1.0);
    const double targetCp = std::clamp(prismaticCoefficient, 0.30, 0.995);
    if (problems && targetCp != prismaticCoefficient)
        problems->push_back("prismatic coefficient clamped to [0.30, 0.995]");

    // Inner solve: given a fore-aft skew, find the common exponent that hits Cp.
    // Cp rises monotonically with the exponent -- a larger exponent keeps the
    // section full further towards the ends.
    const auto withSkew = [&](double skew) {
        AreaCurve c = curve;   // carries the end values into every trial
        const auto cpFor = [&](double n) {
            c.aftExponent = n * std::exp(skew);
            c.forwardExponent = n * std::exp(-skew);
            return c.prismaticCoefficient();
        };
        const double n = bisect(0.05, 200.0, targetCp, cpFor);
        c.aftExponent = n * std::exp(skew);
        c.forwardExponent = n * std::exp(-skew);
        return c;
    };

    // Outer solve: skew for LCB. Positive skew fattens the aft body, which moves
    // the centre of buoyancy aft, so LCB falls as skew rises.
    const auto lcbFor = [&](double skew) { return withSkew(skew).lcbFraction(); };
    const double loLcb = lcbFor(1.6), hiLcb = lcbFor(-1.6);
    const double target = std::clamp(lcbFraction, std::min(loLcb, hiLcb), std::max(loLcb, hiLcb));
    if (problems && std::abs(target - lcbFraction) > 1e-9)
        problems->push_back("LCB target unreachable at this Cp; clamped to " +
                            std::to_string(target));

    return withSkew(bisect(-1.6, 1.6, target, lcbFor));
}

// --- Generation --------------------------------------------------------------

TriMesh makeHullFromParticulars(const HullParticulars& p, std::vector<std::string>* problems) {
    if (problems)
        for (const std::string& s : validateParticulars(p)) problems->push_back(s);

    const double bilge = bilgeRadiusForMidshipCoefficient(p.beam, p.draft, p.midshipCoefficient);

    // A bilge radius is capped by the draft and the half-beam, so a section can
    // only be so fine. Say so rather than quietly building a fuller ship: the
    // request is then met in Cp and missed in Cb, which is the confusing way
    // round.
    const double achievableCm = midshipCoefficientForBilgeRadius(p.beam, p.draft, bilge);
    if (problems && achievableCm > p.midshipCoefficient + 1e-6)
        problems->push_back("midship coefficient " + std::to_string(p.midshipCoefficient) +
                            " needs a bilge radius larger than the draft or half-beam; the "
                            "finest achievable here is " + std::to_string(achievableCm));

    // Cp is asked of the section this hull will actually have, not of the one
    // that was requested and clamped away.
    const double cp = p.blockCoefficient / std::max(achievableCm, 1e-6);
    const AreaCurve curve = solveAreaCurve(cp, p.lcbFraction, p.transomFraction, p.stemFraction,
                                          p.parallelMiddleBodyFraction, problems);

    // Waterlines: clustered towards the baseline and towards the design draft,
    // because the bilge radius lives at one end and the waterplane at the other,
    // and a uniform ladder resolves neither. Above the draft the section is
    // wall-sided, so a couple of levels carry it to the deck.
    const int below = std::max(4, p.waterlineCount * 3 / 4);
    const int above = std::max(2, p.waterlineCount - below);
    std::vector<double> waterlines;
    waterlines.reserve(static_cast<std::size_t>(below + above));
    for (int k = 0; k < below; ++k) {
        const double t = static_cast<double>(k) / (below - 1);
        waterlines.push_back(p.draft * 0.5 * (1.0 - std::cos(kPi * t)));
    }
    for (int k = 1; k <= above; ++k)
        waterlines.push_back(p.draft + (p.depth - p.draft) * k / above);

    // Station positions: uniform without a parallel middle body, and concentrated
    // in the tapers with one.
    //
    // Uniform spacing spends stations where the section never changes and starves
    // the ends where it changes fastest, and the taper is *steeper* with a
    // parallel body because the same fall-off is squeezed into (1 - p) of the
    // length. Measured on the S-175 form with a 0.40 parallel body, uniform
    // spacing at 41 stations put the block coefficient 0.76% high; it needed 161
    // to come back to 0.04%. A section that is constant needs two stations to
    // describe it exactly, so the rest belong in the tapers.
    const int stations = std::max(5, p.stationCount | 1);
    const double parallel = std::clamp(p.parallelMiddleBodyFraction, 0.0, 0.98);
    std::vector<double> positions;
    positions.reserve(static_cast<std::size_t>(stations));
    if (parallel <= 1e-6) {
        for (int i = 0; i < stations; ++i)
            positions.push_back(-1.0 + 2.0 * i / (stations - 1));
    } else {
        // A handful across the flat span, purely so the mesh has faces there;
        // the geometry is exact with two.
        constexpr int kMiddle = 5;
        const int perTaper = std::max(3, (stations - kMiddle) / 2);
        for (int i = 0; i < perTaper; ++i)
            positions.push_back(-1.0 + (1.0 - parallel) * i / (perTaper - 1));
        for (int i = 1; i < kMiddle - 1; ++i)
            positions.push_back(-parallel + 2.0 * parallel * i / (kMiddle - 1));
        for (int i = 0; i < perTaper; ++i)
            positions.push_back(parallel + (1.0 - parallel) * i / (perTaper - 1));
    }

    std::vector<Station> out;
    out.reserve(positions.size());
    for (double u : positions) {
        const double scale = curve(u);
        Station station;
        station.x = 0.5 * p.lengthPp * u;
        station.halfBeam.reserve(waterlines.size());
        // Every station is the midship section scaled in breadth, which makes the
        // sectional area curve exact by construction: area scales with breadth
        // when the shape does not change.
        for (double z : waterlines)
            station.halfBeam.push_back(scale * midshipHalfBreadth(p.beam, p.draft, bilge, z));
        out.push_back(station);
    }
    return makeHullFromStations(out, waterlines);
}

// --- Measurement -------------------------------------------------------------

HullCoefficients measureHull(const TriMesh& hull, double draft, double lengthPp, double beam) {
    HullCoefficients c;
    if (hull.verts.empty() || draft <= 0) return c;

    Vec3 lo = hull.verts[0], hi = hull.verts[0];
    for (const Vec3& v : hull.verts) {
        lo = {std::min(lo.x, v.x), std::min(lo.y, v.y), std::min(lo.z, v.z)};
        hi = {std::max(hi.x, v.x), std::max(hi.y, v.y), std::max(hi.z, v.z)};
    }
    c.lengthPp = lengthPp > 0 ? lengthPp : hi.x - lo.x;
    c.beam = beam > 0 ? beam : hi.y - lo.y;
    c.draft = draft;
    if (c.lengthPp <= 0 || c.beam <= 0) return c;

    const double waterline = lo.z + draft;
    const VolumeIntegral submerged = integrateBelowPlane(hull, {0, 0, 1}, waterline);
    c.displacedVolume = submerged.volume;
    c.blockCoefficient = submerged.volume / (c.lengthPp * c.beam * draft);
    c.lcbFraction = submerged.volume > 0 ? submerged.centroid.x / c.lengthPp : 0.0;

    // Waterplane area by finite difference of the volume, which is what it is:
    // dV/dz at the waterline. The same trick the hydrostatics use.
    const double h = std::min(0.05, 0.02 * draft);
    const double dV = integrateBelowPlane(hull, {0, 0, 1}, waterline + h).volume -
                      integrateBelowPlane(hull, {0, 0, 1}, waterline - h).volume;
    c.waterplaneCoefficient = dV / (2.0 * h) / (c.lengthPp * c.beam);

    // Midship section area from a thin slab, as radiationHullFromMesh does.
    const double slab = 0.005 * c.lengthPp;
    const double wide = (hi.y - lo.y) + 1.0;
    const TriMesh midship = clipToBox(hull, {-0.5 * slab, lo.y - wide, lo.z - wide},
                                      {0.5 * slab, hi.y + wide, waterline});
    const double area = integrate(midship).volume / slab;
    c.midshipCoefficient = area / (c.beam * draft);
    c.prismaticCoefficient =
        c.midshipCoefficient > 0 ? c.blockCoefficient / c.midshipCoefficient : 0.0;
    return c;
}

std::vector<std::string> validateParticulars(const HullParticulars& p) {
    std::vector<std::string> problems;
    if (p.lengthPp <= 0) problems.push_back("length is not positive");
    if (p.beam <= 0) problems.push_back("beam is not positive");
    if (p.draft <= 0) problems.push_back("draft is not positive");
    if (p.depth <= p.draft) problems.push_back("moulded depth is not above the design draft");
    if (p.blockCoefficient <= 0.2 || p.blockCoefficient >= 1.0)
        problems.push_back("block coefficient outside (0.2, 1.0)");
    if (p.midshipCoefficient <= 0.4 || p.midshipCoefficient > 1.0)
        problems.push_back("midship coefficient outside (0.4, 1.0]");
    if (p.blockCoefficient > p.midshipCoefficient)
        problems.push_back("Cb exceeds Cm, which makes Cp > 1 -- the ship would be fuller "
                           "along its length than its own midship section");
    if (std::abs(p.lcbFraction) > 0.10)
        problems.push_back("LCB more than 10% of Lpp from midship is not a ship");
    if (p.lengthPp / p.beam < 2.0)
        problems.push_back("length/beam below 2 is a pontoon, and strip theory will not apply");
    if (p.stationCount < 11)
        problems.push_back("fewer than 11 stations will not resolve the ends");
    return problems;
}

// --- Reference ships ---------------------------------------------------------

HullParticulars kvlcc2Particulars() {
    // Matches the HullParams in propulsion.cpp, so the two describe one ship.
    HullParticulars p;
    p.lengthPp = 320.0;
    p.beam = 58.0;
    p.draft = 20.8;
    p.depth = 30.0;
    p.blockCoefficient = 0.8098;
    p.midshipCoefficient = 0.998;
    p.lcbFraction = 0.035;
    p.transomFraction = 0.0;
    // A VLCC is close to half parallel. Typical for the type rather than a
    // published figure for this hull, and it barely moves the coefficients.
    p.parallelMiddleBodyFraction = 0.45;
    p.stationCount = 41;
    return p;
}

HullParticulars s175Particulars() {
    HullParticulars p;
    p.lengthPp = 175.0;
    p.beam = 25.4;
    p.draft = 9.5;
    p.depth = 15.4;
    p.blockCoefficient = 0.572;
    p.midshipCoefficient = 0.98;
    p.lcbFraction = -0.02;
    p.transomFraction = 0.10;
    p.parallelMiddleBodyFraction = 0.18;
    p.stationCount = 41;
    return p;
}

}  // namespace sim
