// SPDX-License-Identifier: MIT
//
// The optics of a two-zone fire, with no device anywhere in it.
//
// Every formula `shaders/smoke.frag` evaluates is written out here as well, in
// double, and `tests/test_smoke_render.cpp` predicts pixels from this file rather
// than asking the shader what it thinks. That is the arrangement `material.hpp`
// and `hull.frag` already have, and its value is that two implementations written
// from the same formula disagree loudly when one of them drifts.
#include "smoke.hpp"

#include <algorithm>
#include <cmath>

namespace gpu {
namespace {

// 8-point Gauss-Legendre on [-1, 1]. Enough for Planck over a decade of
// wavelength; the refinement study in tests/test_smoke.cpp is what says so.
constexpr double kGaussNode[8] = {-0.9602898564975363, -0.7966664774136267,
                                  -0.5255324099163290, -0.1834346424956498,
                                  0.1834346424956498,  0.5255324099163290,
                                  0.7966664774136267,  0.9602898564975363};
constexpr double kGaussWeight[8] = {0.1012285362903763, 0.2223810344533745,
                                    0.3137066458778873, 0.3626837833783620,
                                    0.3626837833783620, 0.3137066458778873,
                                    0.2223810344533745, 0.1012285362903763};

// A direction component that a reciprocal can be taken of. The slab test wants
// `1 / d`, and a ray exactly parallel to a face has `d = 0`; substituting a tiny
// finite number gives an intersection at +-1e30, which the min/max below reject
// or accept correctly, and -- unlike an infinity -- can never meet a zero and
// produce a NaN. **The shader does the same substitution with the same constant**,
// so the two agree on the parallel case rather than merely both being plausible.
constexpr double kParallel = 1e-30;

double safeInverse(double d) {
    return 1.0 / (std::abs(d) < kParallel ? (d < 0 ? -kParallel : kParallel) : d);
}

}  // namespace

// --- Blackbody ----------------------------------------------------------------

double stefanBoltzmann() {
    const double pi5 = sim::kPi * sim::kPi * sim::kPi * sim::kPi * sim::kPi;
    const double k4 = kBoltzmann * kBoltzmann * kBoltzmann * kBoltzmann;
    const double h3 = kPlanck * kPlanck * kPlanck;
    return 2.0 * pi5 * k4 / (15.0 * kLightSpeed * kLightSpeed * h3);
}

double planckRadiance(double wavelength, double temperature) {
    if (wavelength <= 0 || temperature <= 0) return 0.0;
    const double x = kPlanck * kLightSpeed / (wavelength * kBoltzmann * temperature);
    // expm1 rather than exp(x) - 1: as `x` goes to zero the subtraction cancels and
    // the relative error grows like eps/x, so the Rayleigh-Jeans limit
    // `B -> 2 c k T / lambda^4` is lost long before the limit is reached.
    //
    // **The obvious justification for this is wrong and was written here first.**
    // The whole-spectrum integral does *not* care: measured over 1 nm to 1 m at
    // 300 K, expm1 gives -6.661e-15 of sigma T^4 / pi and exp(x) - 1 gives
    // -6.772e-15, because the region where the cancellation bites carries almost
    // none of the power. What does care is the limit itself, which
    // `tests/test_smoke.cpp` asserts as a convergence rather than at a point --
    // and nothing tests a comment, which is how this one survived being false.
    const double denominator = std::expm1(x);
    if (!(denominator > 0.0)) return 0.0;
    const double l2 = wavelength * wavelength;
    return 2.0 * kPlanck * kLightSpeed * kLightSpeed / (l2 * l2 * wavelength) / denominator;
}

double blackbodyRadiance(double temperature) {
    if (temperature <= 0) return 0.0;
    const double t2 = temperature * temperature;
    return stefanBoltzmann() * t2 * t2 / sim::kPi;
}

double blackbodyBandRadiance(double lo, double hi, double temperature, int intervals) {
    if (!(hi > lo) || lo <= 0 || temperature <= 0 || intervals < 1) return 0.0;
    // Uniform in log wavelength: the substitution lambda = e^u turns
    // `integral B dlambda` into `integral B(e^u) e^u du`, which spreads the nodes
    // over the decades rather than piling them at the long end.
    const double uLo = std::log(lo), uHi = std::log(hi);
    const double step = (uHi - uLo) / intervals;
    double total = 0.0;
    for (int i = 0; i < intervals; ++i) {
        const double a = uLo + step * i;
        const double half = 0.5 * step;
        const double mid = a + half;
        for (int g = 0; g < 8; ++g) {
            const double u = mid + half * kGaussNode[g];
            const double lambda = std::exp(u);
            total += kGaussWeight[g] * half * planckRadiance(lambda, temperature) * lambda;
        }
    }
    return total;
}

void blackbodyBands(double temperature, double out[3]) {
    // **Four intervals, not the sixty-four the wide-range default uses.** A display
    // band is 100 nm of a smooth function and an 8-point Gauss rule nearly
    // resolves it outright. Measured against a 4096-interval reference: four
    // intervals agree to 5e-13 everywhere from 300 K to 2500 K, and two would do
    // as well *except* in the blue band at 300 K, where the rule sits far down the
    // Wien tail and slips to 1e-8 -- on a value of 4.5e-35 W/(m^2 sr), so it could
    // not matter, which is exactly the kind of reasoning that later turns out to
    // have been about the wrong band. Four costs nothing and needs no argument.
    //
    // It is worth measuring at all because this is the only part of the file that
    // runs per layer per frame: at sixty-four it cost 40 us a call, which is 1.4 ms
    // a frame over a ship's compartments, and a millisecond spent on a colour that
    // is already right in the fourteenth digit is a millisecond spent on nothing.
    for (int c = 0; c < 3; ++c)
        out[c] = blackbodyBandRadiance(kBandEdges[2 - c], kBandEdges[3 - c], temperature,
                                       kDisplayBandIntervals);
}

void emissiveColour(double temperature, const SmokeShading& shading, double out[3]) {
    out[0] = out[1] = out[2] = 0.0;
    double reference[3];
    blackbodyBands(shading.referenceTemperature, reference);
    const double peak = std::max({reference[0], reference[1], reference[2]});
    if (!(peak > 0.0)) return;
    double bands[3];
    blackbodyBands(temperature, bands);
    for (int c = 0; c < 3; ++c) out[c] = bands[c] / peak;
}

// --- The medium ---------------------------------------------------------------

sim::Vec3 SmokeVolume::centreWorld() const {
    return rotation * ((lo + hi) * 0.5) + translation;
}

SmokeLayer layerFrom(const sim::fire::Layer& layer, double volume, const SmokeShading& shading) {
    SmokeLayer out;
    if (volume > 0) {
        // Concentration from the model's own layer volume, not from the drawn box:
        // kg of tracer per m^3 is a statement the zone model makes, and the box is
        // only where the path length comes from.
        const double concentration = std::max(layer.products, 0.0) / volume;
        out.extinction = shading.massExtinction * concentration;
    }
    emissiveColour(layer.temperature(), shading, out.emission);
    return out;
}

std::vector<SmokeVolume> volumesFromFire(const sim::fire::Model& model, const sim::Ship& ship,
                                         const SmokeShading& shading) {
    std::vector<SmokeVolume> volumes;
    volumes.reserve(model.gas.size());
    const sim::Mat3 rotation = ship.state.orientation.toMat3();

    for (const sim::fire::GasCompartment& g : model.gas) {
        if (!(g.floorArea > 0) || !(g.ceilingZ > g.floorZ)) continue;

        // Plan extents: the compartment's own aspect ratio in plan, scaled so the
        // area is the model's `floorArea` exactly. A gas space with no ship
        // compartment behind it has no aspect ratio to borrow and gets a square.
        double aspect = 1.0;
        double cx = 0.0, cy = 0.0;
        if (g.shipCompartment >= 0 &&
            g.shipCompartment < static_cast<int>(ship.compartments.size())) {
            const sim::Compartment& c =
                ship.compartments[static_cast<std::size_t>(g.shipCompartment)];
            const double bx = c.bboxHi.x - c.bboxLo.x;
            const double by = c.bboxHi.y - c.bboxLo.y;
            if (bx > 0 && by > 0) aspect = bx / by;
            cx = 0.5 * (c.bboxLo.x + c.bboxHi.x);
            cy = 0.5 * (c.bboxLo.y + c.bboxHi.y);
        }
        const double sy = std::sqrt(g.floorArea / aspect);
        const double sx = aspect * sy;

        SmokeVolume v;
        v.name = g.name;
        v.lo = {cx - 0.5 * sx, cy - 0.5 * sy, g.floorZ};
        v.hi = {cx + 0.5 * sx, cy + 0.5 * sy, g.ceilingZ};
        // Clamped into the prism rather than trusted. `interfaceZ()` is
        // `ceilingZ - upperVolume / floorArea` and a layer that has grown past the
        // whole space would put the plane under the floor, where the two path
        // lengths stop meaning what they are named.
        v.interfaceZ = std::clamp(g.interfaceZ(), g.floorZ, g.ceilingZ);
        const double upperVolume = std::clamp(g.upperVolume(), 0.0, g.gasVolume);
        v.upper = layerFrom(g.upper, upperVolume, shading);
        v.lower = layerFrom(g.lower, std::max(g.gasVolume - upperVolume, 0.0), shading);
        v.rotation = rotation;
        v.translation = ship.state.position;
        volumes.push_back(std::move(v));
    }
    return volumes;
}

// --- The transfer integral ------------------------------------------------------

double transmittance(double extinction, double distance) {
    const double tau = extinction * distance;
    if (!(tau > 0.0)) return 1.0;
    return std::exp(-tau);
}

bool intersectVolume(const SmokeVolume& volume, const sim::Vec3& origin,
                     const sim::Vec3& direction, double maxDistance, double& tEnter,
                     double& tExit) {
    // Into the body frame. The transform is a rotation and `direction` is a unit
    // vector, so the ray parameter is metres in both frames and no rescaling is
    // needed -- which is the reason the prism is stored in the body frame at all.
    const sim::Mat3 bodyFromWorld = volume.rotation.transposed();
    const sim::Vec3 ob = bodyFromWorld * (origin - volume.translation);
    const sim::Vec3 db = bodyFromWorld * direction;

    double t0 = 0.0, t1 = maxDistance;
    for (int axis = 0; axis < 3; ++axis) {
        const double inverse = safeInverse(db[axis]);
        double entry = (volume.lo[axis] - ob[axis]) * inverse;
        double exit = (volume.hi[axis] - ob[axis]) * inverse;
        if (entry > exit) std::swap(entry, exit);
        t0 = std::max(t0, entry);
        t1 = std::min(t1, exit);
    }
    tEnter = t0;
    tExit = t1;
    return t1 > t0;
}

LayerPath layerPath(const SmokeVolume& volume, const sim::Vec3& origin,
                    const sim::Vec3& direction, double tEnter, double tExit) {
    LayerPath path;
    const double span = tExit - tEnter;
    if (!(span > 0)) return path;

    const sim::Mat3 bodyFromWorld = volume.rotation.transposed();
    const sim::Vec3 ob = bodyFromWorld * (origin - volume.translation);
    const sim::Vec3 db = bodyFromWorld * direction;

    if (std::abs(db.z) < kParallel) {
        // Horizontal in the body frame: the whole segment sits in one layer, and
        // which one is decided at the entry height. `upperFirst` is not read.
        const bool inUpper = ob.z >= volume.interfaceZ;
        path.upper = inUpper ? span : 0.0;
        path.lower = inUpper ? 0.0 : span;
        return path;
    }

    const double cross = (volume.interfaceZ - ob.z) / db.z;
    if (db.z < 0) {
        // Descending: above the interface first, so the upper layer comes first.
        path.upperFirst = true;
        path.upper = std::clamp(std::min(tExit, cross), tEnter, tExit) - tEnter;
    } else {
        path.upperFirst = false;
        path.upper = tExit - std::clamp(std::max(tEnter, cross), tEnter, tExit);
    }
    path.upper = std::clamp(path.upper, 0.0, span);
    path.lower = span - path.upper;
    return path;
}

void compositeOver(const SmokeVolume& volume, const LayerPath& path, const double background[3],
                   double out[3]) {
    const double tUpper = transmittance(volume.upper.extinction, path.upper);
    const double tLower = transmittance(volume.lower.extinction, path.lower);
    // Kirchhoff: the emissivity of a slab is the same `1 - exp(-k d)` that its
    // transmittance is one minus. A layer with no soot in it does not glow, at any
    // temperature.
    const double whole = tUpper * tLower;
    for (int c = 0; c < 3; ++c) {
        const double eUpper = volume.upper.emission[c] * (1.0 - tUpper);
        const double eLower = volume.lower.emission[c] * (1.0 - tLower);
        const double source = path.upperFirst ? eUpper + tUpper * eLower
                                              : eLower + tLower * eUpper;
        out[c] = source + whole * background[c];
    }
}


// --- The camera ----------------------------------------------------------------

bool depthBasisFrom(const float mvp[16], DepthBasis& out) {
    if (mvp == nullptr) return false;
    // Column-major, m[column * 4 + row], matching sim::Mat4 and GLSL.
    //
    // For `perspective * lookAt`, with `s = dot(forward, p - eye)` the distance
    // along the view axis:
    //
    //     row 3 = ( forward, -dot(forward, eye) )        so  clip.w = s
    //     row 2 = ( -a forward,  a dot(forward, eye) + b )   so  clip.z = b - a s
    //
    // and therefore `depth = clip.z / clip.w = -a + b / s`, which inverts to
    // `s = b / (depth + a)`. Both constants are already in the matrix.
    const double w[3] = {mvp[0 * 4 + 3], mvp[1 * 4 + 3], mvp[2 * 4 + 3]};
    const double w3 = mvp[3 * 4 + 3];
    const double z[3] = {mvp[0 * 4 + 2], mvp[1 * 4 + 2], mvp[2 * 4 + 2]};
    const double z3 = mvp[3 * 4 + 2];

    const sim::Vec3 forward{w[0], w[1], w[2]};
    const double norm = sim::length(forward);
    // An orthographic projection has row 3 equal to (0, 0, 0, 1) and no centre of
    // projection; it is refused rather than approximated. The bracket is loose
    // because the row is a unit vector for any `perspective * lookAt` and zero for
    // anything else -- there is no third case to be careful about.
    if (!(norm > 0.5) || !(norm < 2.0)) return false;
    out.forward = forward / norm;

    // `a` from the largest component of row 3 rather than from a fixed one: a
    // camera looking down a world axis genuinely has zeros in the other two, and
    // dividing by one of those is how a camera-derived constant becomes an
    // infinity that still renders something.
    int best = 0;
    for (int i = 1; i < 3; ++i)
        if (std::abs(w[i]) > std::abs(w[best])) best = i;
    out.a = -z[best] / w[best];
    out.b = z3 + out.a * w3;
    return std::abs(out.b) > 0.0;
}

}  // namespace gpu
