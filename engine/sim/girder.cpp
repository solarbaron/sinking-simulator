// SPDX-License-Identifier: MIT
#include "girder.hpp"

#include <algorithm>
#include <cmath>

namespace sim {

HullGirder integrateGirder(const std::vector<double>& x,
                           const std::vector<double>& weightPerLength,
                           const std::vector<double>& buoyancyPerLength) {
    HullGirder g;
    const std::size_t n = x.size();
    if (n < 2 || weightPerLength.size() != n || buoyancyPerLength.size() != n) return g;

    g.stations.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
        GirderStation& s = g.stations[i];
        s.x = x[i];
        s.weightPerLength = weightPerLength[i];
        s.buoyancyPerLength = buoyancyPerLength[i];
        s.loadPerLength = weightPerLength[i] - buoyancyPerLength[i];
    }

    // Shear, then moment, both trapezoidal. Exact for a piecewise-linear load,
    // which is what the station values define.
    g.stations[0].shear = 0.0;
    for (std::size_t i = 1; i < n; ++i) {
        const double dx = x[i] - x[i - 1];
        g.stations[i].shear = g.stations[i - 1].shear +
                              0.5 * dx * (g.stations[i].loadPerLength +
                                          g.stations[i - 1].loadPerLength);
    }
    g.stations[0].moment = 0.0;
    for (std::size_t i = 1; i < n; ++i) {
        const double dx = x[i] - x[i - 1];
        g.stations[i].moment =
            g.stations[i - 1].moment + 0.5 * dx * (g.stations[i].shear + g.stations[i - 1].shear);
    }

    for (const GirderStation& s : g.stations) {
        const double dx = std::abs(s.x - g.stations.front().x);
        (void)dx;
        if (std::abs(s.shear) > std::abs(g.maxShear)) {
            g.maxShear = s.shear;
            g.maxShearX = s.x;
        }
        if (std::abs(s.moment) > std::abs(g.maxMoment)) {
            g.maxMoment = s.moment;
            g.maxMomentX = s.x;
        }
        g.totalWeight += s.weightPerLength;
        g.totalBuoyancy += s.buoyancyPerLength;
    }
    // Totals as integrals, not as sums of station values.
    g.totalWeight = 0;
    g.totalBuoyancy = 0;
    for (std::size_t i = 1; i < n; ++i) {
        const double dx = x[i] - x[i - 1];
        g.totalWeight += 0.5 * dx * (weightPerLength[i] + weightPerLength[i - 1]);
        g.totalBuoyancy += 0.5 * dx * (buoyancyPerLength[i] + buoyancyPerLength[i - 1]);
    }

    // The free-end conditions. Both are zero at the after perpendicular by
    // construction, so the forward end carries the whole residual.
    //
    // Scaled against the ship, not against the curve's own peak. Normalising by
    // the peak of the same quantity is degenerate: a correctly balanced ship
    // carrying essentially no bending moment has a peak that *is* the residual,
    // and the ratio comes out at 1.0 — reporting a perfect calculation as a total
    // failure. The reference here is W L / 8, the moment a fully concentrated
    // load would produce, and W / 2, the shear it would produce; those are the
    // largest values the quantities can physically reach.
    const double length = x.back() - x.front();
    const double shearScale = std::max(0.5 * g.totalWeight, 1e-9);
    const double momentScale = std::max(g.totalWeight * length / 8.0, 1e-9);
    g.shearClosure = g.stations.back().shear / shearScale;
    g.momentClosure = g.stations.back().moment / momentScale;
    return g;
}

std::vector<double> girderStations(const Ship& ship, int count) {
    std::vector<double> x;
    const int n = std::max(3, count);
    const double lo = ship.hullLo.x, hi = ship.hullHi.x;
    if (!(hi > lo)) return x;
    x.reserve(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) x.push_back(lo + (hi - lo) * i / (n - 1));
    return x;
}

std::vector<double> buoyancyDistribution(const Ship& ship, const Sea& sea,
                                         const std::vector<double>& x) {
    std::vector<double> out(x.size(), 0.0);
    if (x.size() < 2) return out;

    const Mat3 R = ship.state.orientation.toMat3();
    TriMesh world = ship.hull;
    for (Vec3& v : world.verts) v = R * v + ship.state.position;

    Vec3 lo = world.verts[0], hi = world.verts[0];
    for (const Vec3& v : world.verts) {
        lo = {std::min(lo.x, v.x), std::min(lo.y, v.y), std::min(lo.z, v.z)};
        hi = {std::max(hi.x, v.x), std::max(hi.y, v.y), std::max(hi.z, v.z)};
    }
    const double wide = (hi.y - lo.y) + (hi.z - lo.z) + 1.0;

    for (std::size_t i = 0; i < x.size(); ++i) {
        // Slab width: the interval this station owns, half-width at the ends.
        const double left = i == 0 ? x[0] : 0.5 * (x[i - 1] + x[i]);
        const double right = i + 1 == x.size() ? x.back() : 0.5 * (x[i] + x[i + 1]);
        const double width = right - left;
        if (width <= 0) continue;

        // Station positions are body-frame; the slab has to be cut in the frame
        // the mesh is now in. For the small angles a girder calculation is valid
        // at, the difference between a body-frame plane and a world-frame one is
        // second order -- but the *offset* is not, so the slab follows the hull.
        const Vec3 forward = R * Vec3{1, 0, 0};
        const double leftOffset = dot(forward, ship.state.position) + left;
        const double rightOffset = dot(forward, ship.state.position) + right;

        TriMesh slab = world;
        slab = clipByPlane(slab, forward * -1.0, -leftOffset);
        slab = clipByPlane(slab, forward, rightOffset);
        if (slab.verts.empty()) continue;

        const VolumeIntegral wet =
            sea.flat() ? integrateBelowPlane(slab, {0, 0, 1}, sea.level)
                       : integrateBelowSurface(slab, [&](double px, double py) {
                             return sea.heightAt(px, py);
                         });
        (void)wide;
        out[i] = ship.seaDensity * kGravity * std::max(0.0, wet.volume) / width;
    }
    return out;
}

std::vector<double> weightDistribution(const Ship& ship, const std::vector<double>& x) {
    std::vector<double> out(x.size(), 0.0);
    if (x.size() < 2) return out;
    const double lo = x.front(), hi = x.back();
    const double length = hi - lo;
    if (!(length > 0)) return out;

    // Lightship: a uniform part plus a linear part, chosen so the distribution
    // integrates to the lightship mass and to its LCG. Writing the linear part
    // about the mid-length makes the two conditions independent:
    //   w(x) = W/L + k (x - xm),  integral = W,  first moment = k L^3 / 12
    // so k follows directly from the required LCG offset.
    const double weight = ship.lightshipMass * kGravity;
    const double middle = 0.5 * (lo + hi);
    const double lcgOffset = ship.lightshipCog.x - middle;
    const double slope = 12.0 * weight * lcgOffset / (length * length * length);
    for (std::size_t i = 0; i < x.size(); ++i)
        out[i] = weight / length + slope * (x[i] - middle);

    // A trapezoid that goes negative is telling you the LCG is further off than a
    // linear distribution can represent. Clamping is wrong -- it silently loses
    // weight -- so the negative part is kept and validateGirder() reports it.

    // Floodwater is not assumed. Each compartment's water is spread over that
    // compartment's own extent, because the simulator knows exactly where it is.
    for (const Compartment& c : ship.compartments) {
        const double mass = c.waterVolume * ship.seaDensity;
        if (mass <= 0) continue;
        const double a = std::max(c.bboxLo.x, lo), b = std::min(c.bboxHi.x, hi);
        const double span = b - a;
        if (!(span > 0)) continue;
        const double perLength = mass * kGravity / span;
        // **By the overlap of each station's own slab, not by whether the station
        // falls inside the compartment.** `integrateGirder` reads this array
        // trapezoidally, which weights station `i` by the slab
        // `[0.5(x[i-1]+x[i]), 0.5(x[i]+x[i+1])]` -- half-width at the ends. Giving
        // every station in `[a, b]` the full `perLength` therefore delivers
        // `M h (W/S)` rather than `W`, and both end stations carry a whole slab
        // where they own half of one.
        //
        // Measured on the suite's own flooding fixture -- an 81-station barge, so
        // `h = 1.5`, with a hold from -15 to +15 whose bulkheads land exactly on
        // stations -- the girder carried **105.00%** of the water in the hold. On
        // the ferry at the default 41 stations no bulkhead coincides and the sign
        // varies by compartment: the forepeak came out +12.5%, the engine room
        // -3.6%. `validateGirder` compares that weight against the true displaced
        // buoyancy at a 2% threshold, so a correctly floating flooded ship could be
        // reported as not floating.
        //
        // The construction below is the one `buoyancyDistribution` already uses
        // fifty lines above, for the same reason. It also fixes the other tail: a
        // compartment shorter than one slab used to contribute nothing at all.
        for (std::size_t i = 0; i < x.size(); ++i) {
            const double left = i == 0 ? x[0] : 0.5 * (x[i - 1] + x[i]);
            const double right = i + 1 == x.size() ? x.back() : 0.5 * (x[i] + x[i + 1]);
            const double width = right - left;
            if (width <= 0) continue;
            const double overlap = std::min(right, b) - std::max(left, a);
            if (overlap <= 0) continue;
            out[i] += perLength * overlap / width;
        }
    }
    return out;
}

bool solveBalanceStep(double a, double b, double c, double d, double force, double moment,
                      double& dz, double& dTrim) {
    const double det = a * d - b * c;
    // The size the determinant would have had if nothing had cancelled. Written
    // as a negated `>` so that a Jacobian which is entirely zero -- no hull in the
    // water at all, every term exactly 0, and `scale` zero with it -- and a NaN
    // one both refuse rather than divide. 1e-12 is four decades above double
    // round-off, so it fires only once the two products agree to within a few
    // thousand ULP, which is past the point where a Jacobian finite-differenced
    // from clipped volumes carries any digits at all.
    const double scale = std::abs(a * d) + std::abs(b * c);
    if (!(std::abs(det) > 1e-12 * scale)) return false;
    dz = (force * d - b * moment) / det;
    dTrim = (a * moment - force * c) / det;
    // Checked before the caller applies it, so a nonsense correction never
    // reaches the ship.
    return std::isfinite(dz) && std::isfinite(dTrim);
}

bool balanceOnWave(Ship& ship, const Sea& sea, int iterations, int* iterationsUsed) {
    if (iterationsUsed) *iterationsUsed = 0;
    const Diagnostics d = ship.diagnostics(sea);
    const double weight = d.displacementMass * kGravity;
    if (!(weight > 0)) return false;
    // A moment residual has to be measured against a lever, and the hull's own
    // length is the only one the problem offers. A ship with no cached hull
    // extent has neither a length nor any stations -- girderStations() refuses
    // the same ship -- so there is nothing here to balance.
    const double length = ship.hullHi.x - ship.hullLo.x;
    if (!(length > 0)) return false;
    // Longitudinal centre of gravity, in world coordinates along the ship.
    const Mat3 R0 = ship.state.orientation.toMat3();
    const double lcgWorld = (R0 * d.centreOfGravity + ship.state.position).x;

    // Two residuals -- net vertical force and net trimming moment -- against two
    // unknowns, sinkage and trim. Newton with a numerical Jacobian; each residual
    // evaluation is one whole-hull integral, so this is a handful of them rather
    // than one per station.
    const auto residuals = [&](const Ship& s, double& force, double& moment) {
        const Mat3 R = s.state.orientation.toMat3();
        TriMesh world = s.hull;
        for (Vec3& v : world.verts) v = R * v + s.state.position;
        const VolumeIntegral wet =
            sea.flat() ? integrateBelowPlane(world, {0, 0, 1}, sea.level)
                       : integrateBelowSurface(world, [&](double px, double py) {
                             return sea.heightAt(px, py);
                         });
        const double buoyancy = s.seaDensity * kGravity * wet.volume;
        force = buoyancy - weight;
        // Trimming moment about the centre of gravity.
        moment = buoyancy * (wet.centroid.x - lcgWorld);
    };

    const auto shifted = [&](double dz, double dTrim) {
        Ship s = ship;
        s.state.position.z += dz;
        if (dTrim != 0.0)
            s.state.orientation = Quat::fromAxisAngle({0, 1, 0}, dTrim) * s.state.orientation;
        return s;
    };

    // Converged when both residuals are small, each against its own scale --
    // which cannot be one scale, because they are not one quantity. `f0` is a
    // force in N and `m0` a moment in N m, and testing the moment against
    // `1e-6 * weight` asked for the centre of buoyancy within a fixed
    // *micrometre* of the centre of gravity, on a ship of any size at all.
    //
    // What that costs was measured rather than assumed, and it is not what it
    // looks like: the criterion is reachable, because Newton drives the residual
    // to round-off and the loop does stop. What it is not is scale free. The same
    // barge at 12 m, 120 m, 360 m and 1200 m stops after 11, 13, 14 and 15 steps
    // against `weight`, climbing with the hull because a micrometre is a finer
    // demand on a longer ship, and after 9 steps at every one of those sizes
    // against `weight * length`. The ferry on a 3 m crest goes from 12 steps to
    // 10. Each step is three whole-hull clip-and-integrate passes, and
    // hullGirder() takes this path on every Tier-0 review.
    //
    // Against `weight * length` the two legs are the same relative tolerance, and
    // the moment one reads as "the centre of buoyancy is within 1e-6 of a length
    // of the centre of gravity".
    const double forceTolerance = 1e-6 * weight;
    const double momentTolerance = 1e-6 * weight * length;

    double dz = 0.0, dTrim = 0.0;
    for (int i = 0; i < iterations; ++i) {
        double f0 = 0, m0 = 0;
        residuals(shifted(dz, dTrim), f0, m0);
        if (std::abs(f0) < forceTolerance && std::abs(m0) < momentTolerance) break;
        if (iterationsUsed) ++*iterationsUsed;

        const double hz = 0.01, ht = 1e-4;
        double f1 = 0, m1 = 0, f2 = 0, m2 = 0;
        residuals(shifted(dz + hz, dTrim), f1, m1);
        residuals(shifted(dz, dTrim + ht), f2, m2);

        const double a = (f1 - f0) / hz, b = (f2 - f0) / ht;
        const double c = (m1 - m0) / hz, d2 = (m2 - m0) / ht;
        double sinkage = 0, trim = 0;
        if (!solveBalanceStep(a, b, c, d2, f0, m0, sinkage, trim)) return false;
        // Damped, because a hull leaving the water makes the Jacobian lie.
        dz -= 0.7 * sinkage;
        dTrim -= 0.7 * trim;
    }

    ship = shifted(dz, dTrim);
    double f = 0, m = 0;
    residuals(ship, f, m);
    return std::abs(f) < 1e-3 * weight;
}

HullGirder hullGirder(const Ship& ship, const Sea& sea, int stationCount, bool balance) {
    Ship poised = ship;
    if (balance) balanceOnWave(poised, sea);
    const std::vector<double> x = girderStations(poised, stationCount);
    return integrateGirder(x, weightDistribution(poised, x),
                           buoyancyDistribution(poised, sea, x));
}

std::vector<GirderStress> girderStress(const HullGirder& girder,
                                       const StructuralMesh& structure, double yieldStrength) {
    std::vector<GirderStress> out;
    if (!(yieldStrength > 0)) return out;
    out.reserve(girder.stations.size());
    for (const GirderStation& station : girder.stations) {
        const HullGirderSection section = hullGirderSection(structure, station.x);
        // No structure here -- the ends, past the last frame. Reporting zero
        // stress would read as "safe" rather than "not asked".
        if (!(section.modulusDeck > 0) || !(section.modulusKeel > 0)) continue;

        GirderStress s;
        s.x = station.x;
        s.moment = station.moment;
        s.modulusDeck = section.modulusDeck;
        s.modulusKeel = section.modulusKeel;
        // Hogging arches the hull, stretching the deck and compressing the keel.
        s.stressDeck = station.moment / section.modulusDeck;
        s.stressKeel = -station.moment / section.modulusKeel;
        s.utilisation = std::max(std::abs(s.stressDeck), std::abs(s.stressKeel)) / yieldStrength;
        out.push_back(s);
    }
    return out;
}

double worstUtilisation(const std::vector<GirderStress>& stresses, double* atX) {
    double worst = 0;
    for (const GirderStress& s : stresses)
        if (s.utilisation > worst) {
            worst = s.utilisation;
            if (atX) *atX = s.x;
        }
    return worst;
}

std::vector<std::string> validateGirder(const HullGirder& g) {
    std::vector<std::string> problems;
    if (g.stations.size() < 3) {
        problems.push_back("fewer than three stations is not a beam");
        return problems;
    }
    // **"of peak" was what these said, and the block above rejects that
    // normalisation by name.** `shearClosure` is the residual over `W/2` and
    // `momentClosure` over `W L / 8` -- deliberately *not* over the curve's own
    // peak, because a balanced ship's peak is its residual and the ratio would
    // read 1.0 on a perfect calculation. On the suite's own 10% imbalance
    // fixture the residual really is 100% of the peak in both quantities while
    // these print 20% and 40%, so the label was wrong by 5x and 2.5x in the one
    // case the tests construct. The numbers were always right; the units in the
    // sentence were not.
    if (std::abs(g.shearClosure) > 0.05)
        problems.push_back("shear does not close at the forward perpendicular (" +
                           std::to_string(100.0 * g.shearClosure) +
                           "% of W/2) -- the ship is not in equilibrium");
    if (std::abs(g.momentClosure) > 0.05)
        problems.push_back("bending moment does not close at the forward perpendicular (" +
                           std::to_string(100.0 * g.momentClosure) + "% of W L / 8)");
    const double imbalance = std::abs(g.totalBuoyancy - g.totalWeight) /
                             std::max(g.totalWeight, 1.0);
    if (imbalance > 0.02)
        problems.push_back("buoyancy and weight differ by " + std::to_string(100.0 * imbalance) +
                           "% -- the ship is not floating at this attitude");
    for (const GirderStation& s : g.stations)
        if (s.weightPerLength < 0) {
            problems.push_back("the trapezoidal lightship distribution goes negative near x = " +
                               std::to_string(s.x) +
                               "; this LCG is further off than a linear weight curve can carry");
            break;
        }
    return problems;
}

}  // namespace sim
