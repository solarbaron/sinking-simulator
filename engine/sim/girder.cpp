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
        for (std::size_t i = 0; i < x.size(); ++i)
            if (x[i] >= a && x[i] <= b) out[i] += perLength;
    }
    return out;
}

bool balanceOnWave(Ship& ship, const Sea& sea, int iterations) {
    const Diagnostics d = ship.diagnostics(sea);
    const double weight = d.displacementMass * kGravity;
    if (!(weight > 0)) return false;
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

    double dz = 0.0, dTrim = 0.0;
    for (int i = 0; i < iterations; ++i) {
        double f0 = 0, m0 = 0;
        residuals(shifted(dz, dTrim), f0, m0);
        if (std::abs(f0) < 1e-6 * weight && std::abs(m0) < 1e-6 * weight) break;

        const double hz = 0.01, ht = 1e-4;
        double f1 = 0, m1 = 0, f2 = 0, m2 = 0;
        residuals(shifted(dz + hz, dTrim), f1, m1);
        residuals(shifted(dz, dTrim + ht), f2, m2);

        const double a = (f1 - f0) / hz, b = (f2 - f0) / ht;
        const double c = (m1 - m0) / hz, d2 = (m2 - m0) / ht;
        const double det = a * d2 - b * c;
        if (std::abs(det) < 1e-30) return false;
        // Damped, because a hull leaving the water makes the Jacobian lie.
        dz -= 0.7 * (f0 * d2 - b * m0) / det;
        dTrim -= 0.7 * (a * m0 - f0 * c) / det;
        if (!std::isfinite(dz) || !std::isfinite(dTrim)) return false;
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

std::vector<std::string> validateGirder(const HullGirder& g) {
    std::vector<std::string> problems;
    if (g.stations.size() < 3) {
        problems.push_back("fewer than three stations is not a beam");
        return problems;
    }
    if (std::abs(g.shearClosure) > 0.05)
        problems.push_back("shear does not close at the forward perpendicular (" +
                           std::to_string(100.0 * g.shearClosure) +
                           "% of peak) -- the ship is not in equilibrium");
    if (std::abs(g.momentClosure) > 0.05)
        problems.push_back("bending moment does not close at the forward perpendicular (" +
                           std::to_string(100.0 * g.momentClosure) + "% of peak)");
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
