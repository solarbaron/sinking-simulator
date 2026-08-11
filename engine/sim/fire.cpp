// SPDX-License-Identifier: MIT
#include "fire.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

namespace sim::fire {
namespace {

constexpr double kInf = std::numeric_limits<double>::infinity();

// Below this a layer is treated as absent rather than as a very cold gas: a
// temperature of U/(m c_v) with m at machine epsilon is a number, not a
// temperature.
constexpr double kMassFloor = 1e-9;      // kg

// Pressure differences smaller than this do not move gas. The same idea as
// `Ship::solveFlowNetwork`'s 1e-3 Pa noise floor, and the same reason: below it
// the orifice law is resolving round-off. Tighter than the ship's because the
// buoyancy pressures that drive a fire are two orders of magnitude smaller than
// the water heads that drive flooding -- a doorway with a 300 K layer behind it
// is working on about 6 Pa per metre of height.
constexpr double kPressureFloor = 1e-9;  // Pa

// Incoming gas joins the receiving compartment's upper layer only if it is
// warmer than the lower layer by at least this much. Without the margin, ambient
// air entering a compartment whose lower layer is *exactly* at ambient would test
// as buoyant and be deposited under the deckhead, which is precisely backwards.
constexpr double kDepositionMargin = 0.5;   // K

// Entrainment is tapered to zero over the last this-much of the compartment's
// volume, as the interface reaches the floor.
//
// This is physics rather than a numerical guard: a plume rising through a lower
// layer that has been consumed has no cool air left to entrain, and the fire is
// then burning inside the hot layer. Without the taper the plume drains the last
// of the lower layer in one step and leaves a residue of mass with no energy --
// a layer at absolute zero, an infinite density, and a NaN one vent integral
// later.
constexpr double kEntrainmentTaperFraction = 0.02;

// The band above saturation over which a spray goes from wetting the deck to
// boiling away, K.
//
// A smoothing of a boundary that is not sharp in reality: drops entering an
// unsaturated gas evaporate well below the nominal boiling point, and drops that
// are already hot evaporate above it. Switching the latent term on at exactly
// 373.15 K instead would put a **discontinuity in the sink** into a term the
// substep controller subdivides against, at the one temperature a drenched layer
// actually settles at -- the ferry's deck sits within this band, not away from
// it.
//
// The width is a modelling choice and it moves a *temperature*, not a mass: the
// evaporating share at equilibrium is fixed by the energy balance
// `Q_fire = flow (cp dT + e L s)`, so a wider band puts the same steam into the
// layer at a higher temperature. `tests/test_fire.cpp` measures that the water
// landing on the deck barely moves with it, which is what the stability half of
// this file depends on.
constexpr double kEvaporationBand = 50.0;   // K

// Sweeps of the implicit pressure solve, which is blocked by *vent* -- see
// `substep`. Each vent's own pair is solved exactly, so a compartment with one
// vent needs one sweep; the sweeps are only there for the coupling between
// several vents on the same compartment, which is genuinely the weak one.
constexpr int kPressureSweeps = 24;

// Watts to kilowatts, for the plume and enclosure correlations, all of which are
// published in kW. One place, one direction.
constexpr double kWattsPerKilowatt = 1000.0;

double cube(double x) { return x * x * x; }

// `x^(3/2)` for x >= 0.
double pow32(double x) { return x <= 0 ? 0.0 : x * std::sqrt(x); }

}  // namespace

// ---------------------------------------------------------------------------
// Plume
// ---------------------------------------------------------------------------

double Plume::flameHeight() const {
    if (heatRelease <= 0) return 0.0;
    const double qkw = heatRelease / kWattsPerKilowatt;
    return 0.235 * std::pow(qkw, 0.4) - 1.02 * diameter;
}

double Plume::virtualOrigin() const {
    if (heatRelease <= 0) return 0.0;
    const double qkw = heatRelease / kWattsPerKilowatt;
    return 0.083 * std::pow(qkw, 0.4) - 1.02 * diameter;
}

double Plume::entrainment(double height) const {
    if (height <= 0 || heatRelease <= 0) return 0.0;
    const double qc = convectiveHeatRelease() / kWattsPerKilowatt;   // kW
    if (qc <= 0) return 0.0;

    const double lf = flameHeight();
    if (lf > 0 && height <= lf) {
        // Inside the flame the plume is still accelerating and the far-field
        // power law does not hold. Heskestad's flame-region result is linear in
        // height and anchored at the flame tip.
        return 0.0056 * qc * (height / lf);
    }

    // Far field, measured from the virtual origin. `z0` is usually negative -- a
    // wide fire behaves like a point source below the pan -- but a small, tall
    // one can put it above the base, so the clamp is not decoration.
    const double z = std::max(height - virtualOrigin(), 0.0);
    return 0.071 * std::cbrt(qc) * std::pow(z, 5.0 / 3.0) + 0.0018 * qc;
}

// ---------------------------------------------------------------------------
// Design fire
// ---------------------------------------------------------------------------

double DesignFire::heatRelease(double t) const {
    if (t <= 0 || peakHeatRelease <= 0) return 0.0;
    const double tGrow =
        growthCoefficient > 0 ? std::sqrt(peakHeatRelease / growthCoefficient) : 0.0;
    if (t < tGrow) return growthCoefficient * t * t;

    const double tSteadyEnd = tGrow + steadyDuration;
    if (t <= tSteadyEnd) return peakHeatRelease;
    // No decay ramp means the fire burns at its peak indefinitely, which is the
    // steady case every enclosure correlation is written for.
    if (decayDuration <= 0) return peakHeatRelease;
    const double into = t - tSteadyEnd;
    if (into >= decayDuration) return 0.0;
    return peakHeatRelease * (1.0 - into / decayDuration);
}

double DesignFire::totalEnergy() const {
    if (peakHeatRelease <= 0) return 0.0;
    if (decayDuration <= 0) return kInf;   // never goes out
    const double tGrow =
        growthCoefficient > 0 ? std::sqrt(peakHeatRelease / growthCoefficient) : 0.0;
    return growthCoefficient * cube(tGrow) / 3.0 + peakHeatRelease * steadyDuration +
           0.5 * peakHeatRelease * decayDuration;
}

// ---------------------------------------------------------------------------
// Layers
// ---------------------------------------------------------------------------

double Layer::temperature() const {
    // Both guards matter, and the second one is not symmetry. A layer that has
    // been emptied down to a residue of mass while its energy has gone to zero
    // reports T = 0, which becomes an infinite density in `sideOf` and a NaN in
    // the next vent integral. `substep` tapers entrainment so the state should
    // not arise; this is the second line of defence, and the tests assert it is
    // never reached on a realistic case.
    if (mass <= kMassFloor || energy <= 0) return kTAmbient;
    return energy / (mass * kCvAir);
}

double Layer::productFraction() const {
    if (mass <= kMassFloor) return 0.0;
    return products / mass;
}

double GasCompartment::pressure() const {
    if (gasVolume <= 0) return kPatm;
    return (kGammaAir - 1.0) * totalEnergy() / gasVolume;
}

double GasCompartment::gaugeAtFloor() const {
    return pressure() - (kPatm - kRhoAmbient * kGravity * floorZ);
}

double GasCompartment::upperVolume() const {
    const double u = totalEnergy();
    if (u <= 0) return 0.0;
    // The exact closure: the volume split *is* the internal-energy split.
    return gasVolume * (upper.energy / u);
}

double GasCompartment::interfaceZ() const {
    if (floorArea <= 0) return ceilingZ;
    const double z = ceilingZ - upperVolume() / floorArea;
    return std::clamp(z, floorZ, ceilingZ);
}

void GasCompartment::fillAmbient(double seedFraction) {
    const double f = std::clamp(seedFraction, 1e-12, 0.5);
    // The absolute pressure that makes the gauge pressure exactly zero, so a
    // compartment nobody has lit produces exactly zero vent flow.
    const double pTarget = kPatm - kRhoAmbient * kGravity * floorZ;
    const double u = pTarget * gasVolume / (kGammaAir - 1.0);
    upper.energy = f * u;
    lower.energy = u - upper.energy;
    upper.mass = upper.energy / (kCvAir * kTAmbient);
    lower.mass = lower.energy / (kCvAir * kTAmbient);
    upper.products = 0;
    lower.products = 0;
}

// ---------------------------------------------------------------------------
// Vent geometry
// ---------------------------------------------------------------------------

VentShape ventShapeFor(const Opening& o) {
    VentShape s;
    if (o.area <= 0) return s;
    const double side = std::sqrt(o.area);
    switch (o.kind) {
        case OpeningKind::Door:
            // A door is tall. 2.0 m is the standard clear height; below about
            // 4 m^2 that leaves a plausible width, and above it the door would
            // have to be absurdly wide, so fall back to square.
            if (o.area <= 4.0) {
                s.height = 2.0;
                s.width = o.area / 2.0;
            } else {
                s.height = side;
                s.width = side;
            }
            break;
        case OpeningKind::Hatch:
            s.horizontal = true;
            s.height = 0.0;
            s.width = side;
            break;
        case OpeningKind::Breach:
        case OpeningKind::Vent:
        case OpeningKind::Pipe:
            s.height = side;
            s.width = side;
            break;
    }
    return s;
}

double VentSide::gaugeAt(double z) const {
    // gauge(z) = gauge(floor) - integral of (rho_gas - rho_ambient) g dz.
    // Subtracting the ambient column is what makes still air give exactly zero at
    // every height, rather than two large numbers that nearly cancel.
    const double zi = std::max(interfaceZ, floorZ);
    if (z <= zi) return gaugeAtFloor - (rhoLower - kRhoAmbient) * kGravity * (z - floorZ);
    return gaugeAtFloor - (rhoLower - kRhoAmbient) * kGravity * (zi - floorZ) -
           (rhoUpper - kRhoAmbient) * kGravity * (z - zi);
}

VentSide ambientSide() {
    VentSide s;
    s.gaugeAtFloor = 0;
    // A finite floor reference, because `gaugeAt` would otherwise evaluate
    // `0.0 * infinity`; both densities are ambient, so the difference against the
    // ambient column is exactly 0.0 at any finite reference anyway.
    //
    // The interface is at *positive* infinity, which makes the whole outside one
    // cool lower layer. That is the meaningful reading -- outside air is not
    // anybody's smoke layer -- and it is what makes `VentResult::fromUpperA/B`
    // say what it claims to: gas drawn in from outdoors did not come out of an
    // upper layer, and reporting that it did would be a diagnostic that lies.
    s.floorZ = 0;
    s.interfaceZ = kInf;
    s.rhoLower = s.rhoUpper = kRhoAmbient;
    s.tLower = s.tUpper = kTAmbient;
    s.yLower = s.yUpper = 0;
    return s;
}

namespace {

// One constant-sign, constant-density band of the vent integral.
void accumulateBand(const Vent& v, const VentSide& a, const VentSide& b, double z1, double z2,
                    double dp1, double dp2, VentResult& r) {
    if (z2 <= z1) return;
    const double u1 = std::abs(dp1), u2 = std::abs(dp2);
    if (std::max(u1, u2) < kPressureFloor) return;

    // The band's midpoint, and **any point in the band would name the same layer**:
    // the span was split at every interface before this was called, so no interface
    // lies strictly inside `[z1, z2]` and `densityAt` cannot change across it.
    // Measured rather than assumed -- sampling at `z1` instead is bit-identical
    // over twenty thousand random two-layer configurations. The midpoint is kept
    // because it says what is meant.
    const double zm = 0.5 * (z1 + z2);
    const bool aToB = (dp1 + dp2) > 0;

    // dp is linear in z over the band, so the integral of sqrt|dp| is closed
    // form. Quadrature would converge at half order against the square root's
    // infinite derivative at the neutral plane, which is why the band was split
    // there in the first place.
    const double slope = (dp2 - dp1) / (z2 - z1);
    const double integral = std::abs(slope) > kPressureFloor
                                ? (2.0 / (3.0 * std::abs(slope))) * std::abs(pow32(u2) - pow32(u1))
                                : std::sqrt(u1) * (z2 - z1);

    const VentSide& donor = aToB ? a : b;
    const double rho = donor.densityAt(zm);
    const double mdot = v.dischargeCoeff * v.width * std::sqrt(2.0 * rho) * integral;
    if (mdot <= 0) return;

    const bool fromUpper = zm >= donor.interfaceZ;
    const double temp = donor.temperatureAt(zm);
    if (aToB) {
        r.massAToB += mdot;
        r.enthalpyAToB += mdot * kCpAir * temp;
        if (fromUpper) r.fromUpperA += mdot;
    } else {
        r.massBToA += mdot;
        r.enthalpyBToA += mdot * kCpAir * temp;
        if (fromUpper) r.fromUpperB += mdot;
    }
}

}  // namespace

VentResult ventMassFlow(const Vent& v, const VentSide& a, const VentSide& b) {
    VentResult r;
    if (!v.open || v.blockedByWater) return r;

    if (v.horizontal) {
        // No height to integrate over. A single pressure difference is all this
        // geometry offers -- and see fire.hpp for what it therefore misses.
        if (v.area <= 0) return r;
        const double z = v.sillZ;
        const double dp = a.gaugeAt(z) - b.gaugeAt(z);
        if (std::abs(dp) < kPressureFloor) return r;
        const VentSide& donor = dp > 0 ? a : b;
        const double rho = donor.densityAt(z);
        const double mdot = v.dischargeCoeff * v.area * std::sqrt(2.0 * rho * std::abs(dp));
        const bool fromUpper = z >= donor.interfaceZ;
        if (dp > 0) {
            r.massAToB = mdot;
            r.enthalpyAToB = mdot * kCpAir * donor.temperatureAt(z);
            if (fromUpper) r.fromUpperA = mdot;
        } else {
            r.massBToA = mdot;
            r.enthalpyBToA = mdot * kCpAir * donor.temperatureAt(z);
            if (fromUpper) r.fromUpperB = mdot;
        }
        return r;
    }

    if (v.width <= 0 || v.soffitZ <= v.sillZ) return r;

    // Break the span at each side's layer interface: within a band both densities
    // are constant, so dp is linear and the neutral plane is exactly locatable.
    //
    // A **fixed** four-element sort, with the unused slots parked on the soffit
    // and the duplicates falling out as zero-width bands the loop already skips.
    // The natural version fills the first `n` slots and sorts `[begin, begin+n)`,
    // and at -O3 GCC cannot bound `n`: it reasons the range could run to
    // `SIZE_MAX/8` and reports five `-Warray-bounds` inside `std::sort`. Nothing
    // below -O3 sees them, which is the entire reason `verify.sh full` compiles
    // the engine there -- the same class of thing that sat on `reduction.cpp`'s
    // triangular solves for as long as they existed.
    std::array<double, 4> pts{v.sillZ, v.soffitZ, v.soffitZ, v.soffitZ};
    if (a.interfaceZ > v.sillZ && a.interfaceZ < v.soffitZ) pts[2] = a.interfaceZ;
    if (b.interfaceZ > v.sillZ && b.interfaceZ < v.soffitZ) pts[3] = b.interfaceZ;
    std::sort(pts.begin(), pts.end());

    for (std::size_t i = 0; i + 1 < pts.size(); ++i) {
        const double z1 = pts[i];
        const double z2 = pts[i + 1];
        if (z2 <= z1) continue;
        const double dp1 = a.gaugeAt(z1) - b.gaugeAt(z1);
        const double dp2 = a.gaugeAt(z2) - b.gaugeAt(z2);

        if (dp1 * dp2 < 0) {
            // The neutral plane: the elevation at which the flow reverses. Split
            // exactly on it. This is the whole reason a doorway is not one orifice
            // with one pressure difference.
            const double zn = z1 + (z2 - z1) * dp1 / (dp1 - dp2);
            r.neutralPlaneZ = zn;
            r.bidirectional = true;
            accumulateBand(v, a, b, z1, zn, dp1, 0.0, r);
            accumulateBand(v, a, b, zn, z2, 0.0, dp2, r);
        } else {
            accumulateBand(v, a, b, z1, z2, dp1, dp2, r);
        }
    }
    return r;
}

// ---------------------------------------------------------------------------
// Account
// ---------------------------------------------------------------------------

// Normalised by the largest term in the account, **including the state itself**.
//
// The fluxes alone are not enough. A model whose compartments only exchange with
// each other has every boundary term at exactly zero, and normalising by those
// would report its own last-bit round-off as a whole-number fraction: 1.1e-8 J
// of drift in 34 MJ of gas came out as 1.1e-8 "of scale", which reads like a
// failure and is a hundred million times better than one.
double Account::energyResidualFraction() const {
    // `suppressionCooling` belongs in the scale for exactly the reason the other
    // fluxes do, and leaving it out would have been invisible on any run with a
    // fire in it: a drencher on a compartment nobody lit removes real joules
    // while `heatReleased` stays at zero, and the residual would then be
    // normalised by the state instead of by the largest thing that moved.
    //
    // Mutation testing confirms that "invisible": deleting this line survives the
    // whole suite, because every fixture that runs a drencher also has a fire,
    // and `heatReleased` dominates the maximum in all of them. It stays because
    // the case it covers is a real one -- a boundary-cooled or already-hot space
    // with no design fire in it -- not because a test currently reaches it.
    const double scale = std::max({std::abs(heatReleased), std::abs(enthalpyIn),
                                   std::abs(enthalpyOut), std::abs(wallLoss),
                                   std::abs(suppressionCooling),
                                   std::abs(energy - initialEnergy), std::abs(energy), 1.0});
    return energyResidual() / scale;
}

double Account::massResidualFraction() const {
    const double scale = std::max({std::abs(massIn), std::abs(massOut),
                                   std::abs(mass - initialMass), std::abs(mass), 1e-9});
    return massResidual() / scale;
}

// ---------------------------------------------------------------------------
// Suppression
// ---------------------------------------------------------------------------

double sprayMassFlow(double area, double litresPerMinutePerSquareMetre) {
    if (area <= 0 || litresPerMinutePerSquareMetre <= 0) return 0.0;
    // Litres of fresh water per minute to kilograms per second. Fresh, not sea:
    // the density here converts a *rule's* unit, and the rules are written in
    // litres. What lands on the deck is put into the ship at the ship's own
    // `seaDensity`, in `applyTo`, because that is the density `Ship` weighs
    // floodwater at and a second one would make the mass depend on which file
    // asked.
    return area * litresPerMinutePerSquareMetre * 1e-3 * kRhoFresh / 60.0;
}

double scupperFlow(double insideHead, double outsideHead) {
    const double a = std::max(insideHead, 0.0);
    const double c = std::max(outsideHead, 0.0);
    if (a == c) return 0.0;
    // Antisymmetric by construction rather than by two branches that have to be
    // kept in step: a port with the sea above it is the same port run backwards.
    const double sign = a > c ? 1.0 : -1.0;
    const double hi = std::max(a, c), lo = std::min(a, c);
    const double drop = hi - lo;
    const double root2g = std::sqrt(2.0 * kGravity);
    // The drowned band, of height `lo`, at the constant head `drop`; then the
    // free band, over which the head falls linearly to zero at the surface.
    return sign * root2g * (lo * std::sqrt(drop) + (2.0 / 3.0) * pow32(drop));
}

// ---------------------------------------------------------------------------
// The boundary the structure sees
// ---------------------------------------------------------------------------

double filmCoefficient(double gasKelvin, double surfaceKelvin, const BoundaryFilm& params) {
    const double g = gasKelvin, s = surfaceKelvin;
    // Factored, not `(g^4 - s^4)/(g - s)`: the quotient form is 0/0 when the gas and
    // the surface are at the same temperature, which is the state every coupled run
    // starts in and the state the exact control has to survive.
    return params.convective +
           params.emissivity * kStefanBoltzmann * (g * g + s * s) * (g + s);
}

WallExchange wallExchange(const GasCompartment& gas,
                          const std::vector<thermal::BoundaryFace>& face,
                          const std::vector<double>& surfaceKelvin, const BoundaryFilm& params,
                          double bandHeight) {
    WallExchange out;
    if (face.size() != surfaceKelvin.size() || face.empty()) return out;

    const double zi = gas.interfaceZ();
    const double tUpper = gas.upper.temperature(), tLower = gas.lower.temperature();

    // The band a face belongs to, as a pure function of its own centroid, so that the
    // same face set always produces the same films in the same order. `zBase` is the
    // lowest face centroid rather than the compartment floor: a caller who meshes
    // part of a bulkhead should get bands over the part it meshed.
    double zBase = face[0].centroid.z;
    for (const thermal::BoundaryFace& f : face)
        if (f.area > 0) zBase = std::min(zBase, f.centroid.z);
    const bool banded = bandHeight > 0;
    std::size_t bands = 2;
    if (banded) {
        double zTop = zBase;
        for (const thermal::BoundaryFace& f : face)
            if (f.area > 0) zTop = std::max(zTop, f.centroid.z);
        bands = static_cast<std::size_t>((zTop - zBase) / bandHeight) + 1;
    }
    const auto bandOf = [&](const thermal::BoundaryFace& f) -> std::size_t {
        if (!banded) return f.centroid.z >= zi ? 0u : 1u;
        const auto b = static_cast<std::size_t>((f.centroid.z - zBase) / bandHeight);
        return std::min(b, bands - 1);
    };

    out.film.assign(bands, thermal::Film{});
    out.area.assign(bands, 0.0);
    out.surface.assign(bands, 0.0);

    // Two passes: the first sorts the faces and takes each band's area-weighted mean
    // surface temperature, the second forms the coefficient at that mean. A single
    // pass cannot, because the coefficient is a function of the mean and the mean is
    // not known until every face in the band has been seen.
    for (std::size_t i = 0; i < face.size(); ++i) {
        const thermal::BoundaryFace& f = face[i];
        if (!(f.area > 0)) continue;
        const std::size_t b = bandOf(f);
        out.area[b] += f.area;
        out.surface[b] += f.area * surfaceKelvin[i];
        out.film[b].face.push_back(f);
    }
    for (std::size_t b = 0; b < bands; ++b) {
        if (out.area[b] <= 0) continue;
        out.surface[b] /= out.area[b];
        // A band's own gas temperature: the layer its faces sit in on average. With
        // no banding this is exact, because the split *is* the interface; with bands
        // it quantises the interface to the band height, which is the price of a
        // film membership that does not move.
        double meanZ = 0;
        for (const thermal::BoundaryFace& f : out.film[b].face) meanZ += f.area * f.centroid.z;
        meanZ /= out.area[b];
        out.film[b].ambient = (banded ? meanZ >= zi : b == 0) ? tUpper : tLower;
        out.film[b].coefficient = filmCoefficient(out.film[b].ambient, out.surface[b], params);
    }

    for (std::size_t i = 0; i < face.size(); ++i) {
        const thermal::BoundaryFace& f = face[i];
        if (!(f.area > 0)) continue;
        const std::size_t b = bandOf(f);
        // Per face, though `sum_i a_i h (T_g - T_s,i)` is identically
        // `A_band h (T_g - <T_s>)` -- one coefficient per band is what makes them
        // the same number, and the two forms agree to 6e-16 relative. It is
        // written per face because `exactHeat` beside it *cannot* be, and the
        // difference between the two is the whole point of the pair.
        const double tg = out.film[b].ambient, ts = surfaceKelvin[i];
        out.heat += f.area * out.film[b].coefficient * (tg - ts);
        out.exactHeat += f.area * (params.convective * (tg - ts) +
                                   params.emissivity * kStefanBoltzmann *
                                       (tg * tg * tg * tg - ts * ts * ts * ts));
    }
    out.linearisationError = out.heat - out.exactHeat;

    for (std::size_t b = 0; b < bands; ++b) {
        out.totalArea += out.area[b];
        out.wallTemperature += out.area[b] * out.surface[b];
        out.wallConductance += out.area[b] * out.film[b].coefficient;
    }
    if (out.totalArea > 0) {
        out.wallTemperature /= out.totalArea;
        out.wallConductance /= out.totalArea;
    }
    return out;
}

// ---------------------------------------------------------------------------
// Model internals
// ---------------------------------------------------------------------------

namespace {

// The per-substep rate accumulators, one entry per gas compartment, plus the
// boundary crossings the account wants. Everything here is a *rate*: nothing is
// multiplied by the step until the step has been accepted, which is what lets a
// trial step be abandoned without leaving a mark on the account.
struct Deltas {
    std::vector<double> dmU, dmL, dEU, dEL, dsU, dsL;
    // Water the drenchers landed less what the freeing ports took back, kg/s.
    // Signed, and per compartment, because `applyTo` writes it per compartment.
    std::vector<double> dWater;
    double heat = 0, radiative = 0, wall = 0, entrained = 0;
    double suppression = 0;
    double waterDelivered = 0, waterEvaporated = 0, waterDrained = 0;
    double enthalpyIn = 0, enthalpyOut = 0;
    double massIn = 0, massOut = 0;
    double productsGenerated = 0, productsOut = 0;
    explicit Deltas(std::size_t n)
        : dmU(n, 0.0), dmL(n, 0.0), dEU(n, 0.0), dEL(n, 0.0), dsU(n, 0.0), dsL(n, 0.0),
          dWater(n, 0.0) {}
};

// The gas column against one side of a vent.
VentSide sideOf(const GasCompartment& g) {
    VentSide s;
    s.gaugeAtFloor = g.gaugeAtFloor();
    s.floorZ = g.floorZ;
    s.interfaceZ = g.interfaceZ();
    s.tUpper = g.upper.temperature();
    s.tLower = g.lower.temperature();
    // Densities for the *buoyancy* term come off `kPatm`, not off the
    // compartment's own pressure.
    //
    // That is deliberate and it is the one place where this file separates the
    // thermodynamic pressure from the hydrostatic one. Density varies by a factor
    // of four across a fire's temperature range and by one part in a hundred
    // thousand across its pressure range, so taking the reference pressure loses
    // nothing physical -- and it buys an exact property that the alternative
    // destroys: at `T = kTAmbient` this returns `kRhoAmbient` **exactly**, so the
    // gauge profile of a compartment full of still ambient air is exactly zero at
    // every height and produces exactly zero flow.
    //
    // Using `g.pressure()` here instead leaves the layer density 2.6e-4 kg/m^3
    // off ambient, which is a spurious buoyancy of a hundredth of a pascal per
    // metre. That is invisible next to a fire and fatal to the control that says
    // an unlit fire changes nothing: it moved 74 of the ship's state doubles.
    s.rhoUpper = kPatm / (kRAir * s.tUpper);
    s.rhoLower = kPatm / (kRAir * s.tLower);
    s.yUpper = g.upper.productFraction();
    s.yLower = g.lower.productFraction();
    return s;
}

// Is there water sitting against this opening?
//
// A deliberate second reading of the rule `Ship::sideStateAt` applies, because
// that one is private and returns a pressure this model does not want. Only the
// phase question is duplicated, and `tests/test_fire.cpp` checks the two agree
// across the ferry's whole opening list rather than trusting that they do.
bool waterAgainst(const Ship& ship, const Sea& sea, int shipCompartment, const Vec3& bodyPos) {
    const Vec3 worldPos = ship.state.orientation.toMat3() * bodyPos + ship.state.position;
    if (shipCompartment == kSea) return worldPos.z < sea.heightAt(worldPos.x, worldPos.y);
    const Compartment& c = ship.compartments[static_cast<std::size_t>(shipCompartment)];
    return c.waterVolume > 1e-9 && worldPos.z < c.surfaceWorldZ;
}

// The pdV-consistent split of a compartment's energy input between its layers.
//
// With `V_u = V U_u / U` and both layers at one pressure, each does `p dV_k/dt`
// of work on the other; substituting the closure gives
// `dU_u/dt = [E_u + (gamma-1) f E] / gamma`. The lower layer's rate is taken as
// the remainder rather than from its own formula, so the two sum to `E` in
// floating point and not merely in algebra.
//
// **That last property is real and it is not what makes the account close**, which
// is worth correcting rather than repeating: this comment used to claim it was.
// Replacing the remainder with the algebraically equal `[E_l + (gamma-1)(1-f)E] /
// gamma` was measured on the sealed room, on a wall-loss-only compartment and on
// the ferry, and the energy residual **does not move** -- 7.8e-16 of the heat
// released either way, and 4.7e-16 for the *substituted* form on one of them. The
// account's floor is set by `E` itself, which is a sum of two nearly cancelling
// enthalpy rates, not by the split of it. The remainder stays because it is the
// cheaper expression and because summing to `E` exactly is a property worth
// having; no test here can tell the two apart, and that is a statement about the
// size of the effect rather than about the tests.
void layerSplit(double energyUpper, double energyTotal, double eUpper, double eLower,
                double& dUpper, double& dLower) {
    const double e = eUpper + eLower;
    const double f = energyTotal > 0 ? energyUpper / energyTotal : 0.0;
    dUpper = (eUpper + (kGammaAir - 1.0) * f * e) / kGammaAir;
    dLower = e - dUpper;
}

// Move one directed stream, resolved into the portion that left the donor's upper
// layer and the portion that left its lower one. Those portions carry different
// temperatures and different product loadings, which is exactly why they are kept
// apart rather than averaged into one enthalpy: a doorway can be venting 600 K
// smoke over the top while drawing 288 K air under the same lintel, and one mean
// temperature for that stream would be a temperature nothing in the compartment
// has.
void applyStream(int fromGas, int toGas, const VentSide& donor, const VentSide& recv,
                 double massUpper, double massLower, Deltas& d) {
    for (int part = 0; part < 2; ++part) {
        const bool upper = part == 0;
        const double mass = upper ? massUpper : massLower;
        if (mass <= 0) continue;
        const double temp = upper ? donor.tUpper : donor.tLower;
        const double yield = upper ? donor.yUpper : donor.yLower;
        const double h = mass * kCpAir * temp;
        const double s = mass * yield;

        if (fromGas == kSea) {
            d.massIn += mass;
            d.enthalpyIn += h;
        } else {
            const std::size_t i = static_cast<std::size_t>(fromGas);
            if (upper) {
                d.dmU[i] -= mass;
                d.dEU[i] -= h;
                d.dsU[i] -= s;
            } else {
                d.dmL[i] -= mass;
                d.dEL[i] -= h;
                d.dsL[i] -= s;
            }
        }

        if (toGas == kSea) {
            d.massOut += mass;
            d.enthalpyOut += h;
            d.productsOut += s;
        } else {
            const std::size_t i = static_cast<std::size_t>(toGas);
            // Buoyancy decides the destination, not which layer it left: gas that
            // came out of a hot layer next door arrives buoyant and runs along the
            // deckhead, and cool air arrives dense and sinks whatever door it came
            // through.
            const bool joinUpper = temp > recv.tLower + kDepositionMargin;
            if (joinUpper) {
                d.dmU[i] += mass;
                d.dEU[i] += h;
                d.dsU[i] += s;
            } else {
                d.dmL[i] += mass;
                d.dEL[i] += h;
                d.dsL[i] += s;
            }
        }
    }
}


}  // namespace

// ---------------------------------------------------------------------------
// Model
// ---------------------------------------------------------------------------

int Model::gasIndexOf(int shipCompartment) const {
    if (shipCompartment == kSea) return kSea;
    for (std::size_t i = 0; i < gas.size(); ++i)
        if (gas[i].shipCompartment == shipCompartment) return static_cast<int>(i);
    return -1;
}

int Model::findGas(std::string_view name) const {
    for (std::size_t i = 0; i < gas.size(); ++i)
        if (gas[i].name == name) return static_cast<int>(i);
    return -1;
}

void Model::attach(const Ship& ship, const std::vector<int>& shipCompartments) {
    gas.clear();
    vents.clear();

    for (int idx : shipCompartments) {
        if (idx < 0 || idx >= static_cast<int>(ship.compartments.size())) continue;
        const Compartment& c = ship.compartments[static_cast<std::size_t>(idx)];
        GasCompartment g;
        g.shipCompartment = idx;
        g.name = c.name;
        // The gas floor is the floodwater surface when there is one and the bottom
        // of the space when there is not. A half-flooded engine room has half the
        // gas space and a much shorter distance for the layer to fall.
        g.floorZ = c.bboxLo.z;
        if (c.waterVolume > 1e-9) g.floorZ = std::max(g.floorZ, c.surfaceOffset);
        g.ceilingZ = c.bboxHi.z;
        g.gasVolume = std::max(c.airVolume(), 0.0);
        const double height = std::max(g.ceilingZ - g.floorZ, 1e-6);
        // Prismatic equivalent: the footprint that, at the space's own height,
        // holds the gas it actually has. A ship compartment is not a box, it is a
        // hull clipped to one, so the bounding box would over-state the floor area
        // by the turn of the bilge and the interface would descend too slowly by
        // exactly that ratio.
        g.floorArea = g.gasVolume / height;
        g.perimeter = 2.0 * ((c.bboxHi.x - c.bboxLo.x) + (c.bboxHi.y - c.bboxLo.y));
        g.fillAmbient();
        gas.push_back(std::move(g));
    }

    for (std::size_t i = 0; i < ship.openings.size(); ++i) {
        const Opening& o = ship.openings[i];
        const int ga = gasIndexOf(o.a);
        const int gb = gasIndexOf(o.b);
        // One end tracked and the other a compartment this model knows nothing
        // about: dropped on purpose. Treating an untracked compartment as the
        // atmosphere would put an infinite reservoir of cool air behind a
        // bulkhead, which is worse than having no path at all.
        const bool aOk = (ga >= 0) || (o.a == kSea);
        const bool bOk = (gb >= 0) || (o.b == kSea);
        if (!aOk || !bOk) continue;
        if (ga < 0 && gb < 0) continue;

        const VentShape shape = ventShapeFor(o);
        Vent v;
        v.opening = static_cast<int>(i);
        v.name = o.name;
        v.a = ga >= 0 ? ga : kSea;
        v.b = gb >= 0 ? gb : kSea;
        v.width = shape.width;
        v.area = o.area;
        v.horizontal = shape.horizontal;
        v.sillZ = o.pos.z - 0.5 * shape.height;
        v.soffitZ = o.pos.z + 0.5 * shape.height;
        v.dischargeCoeff = o.dischargeCoeff;
        v.open = o.open;
        vents.push_back(std::move(v));
    }

    resetAccount();
    time = 0;
}

void Model::resetAccount() {
    appliedMass_.assign(gas.size(), 0.0);
    pendingWater_.assign(gas.size(), 0.0);
    for (std::size_t i = 0; i < gas.size(); ++i)
        appliedMass_[i] = gas[i].pressure() * gas[i].gasVolume / (kRAir * kTAmbient);

    account = Account{};
    for (const GasCompartment& g : gas) {
        account.initialEnergy += g.totalEnergy();
        account.initialMass += g.totalMass();
        account.products += g.upper.products + g.lower.products;
    }
    account.energy = account.initialEnergy;
    account.mass = account.initialMass;
}

std::vector<std::string> Model::validate() const {
    std::vector<std::string> problems;
    for (const GasCompartment& g : gas) {
        if (g.gasVolume <= 0) problems.push_back("gas space '" + g.name + "' has no volume");
        if (g.ceilingZ <= g.floorZ)
            problems.push_back("gas space '" + g.name + "' has a ceiling at or below its floor");
        if (g.floorArea <= 0) problems.push_back("gas space '" + g.name + "' has no floor area");
    }
    const int n = static_cast<int>(gas.size());
    for (const Vent& v : vents) {
        if (v.a != kSea && (v.a < 0 || v.a >= n))
            problems.push_back("vent '" + v.name + "' references a gas space that does not exist");
        if (v.b != kSea && (v.b < 0 || v.b >= n))
            problems.push_back("vent '" + v.name + "' references a gas space that does not exist");
        if (v.a == v.b) problems.push_back("vent '" + v.name + "' connects a space to itself");
        if (!v.horizontal && v.soffitZ <= v.sillZ)
            problems.push_back("vent '" + v.name + "' has no height to integrate over");
    }
    for (const DesignFire& f : fires) {
        if (f.compartment < 0 || f.compartment >= n)
            problems.push_back("fire '" + f.name + "' is in a gas space that does not exist");
        if (f.diameter <= 0) problems.push_back("fire '" + f.name + "' has non-positive diameter");
        if (f.convectiveFraction <= 0 || f.convectiveFraction > 1)
            problems.push_back("fire '" + f.name + "' has a convective fraction outside (0, 1]");
    }
    for (const Drencher& s : drenchers) {
        if (s.gasCompartment < 0 || s.gasCompartment >= n)
            problems.push_back("drencher '" + s.name + "' is in a gas space that does not exist");
        if (s.evaporatedFraction < 0 || s.evaporatedFraction > 1)
            problems.push_back("drencher '" + s.name +
                               "' has an evaporated fraction outside [0, 1]");
        if (s.flow < 0) problems.push_back("drencher '" + s.name + "' has a negative flow");
        // A nozzle hotter than saturation removes no sensible heat and the model
        // would quietly deliver water that cools nothing.
        if (s.waterTemperature >= kTSaturation)
            problems.push_back("drencher '" + s.name + "' supplies water at or above boiling");
    }
    for (const Scupper& sc : scuppers) {
        if (sc.gasCompartment < 0 || sc.gasCompartment >= n) {
            problems.push_back("scupper '" + sc.name + "' is in a gas space that does not exist");
            continue;
        }
        // A scupper drains a ship compartment's water, so a gas space that is not
        // attached to one has nothing for it to drain. Silent otherwise: the port
        // would simply never run, which looks exactly like a blocked one.
        if (gas[static_cast<std::size_t>(sc.gasCompartment)].shipCompartment == kSea)
            problems.push_back("scupper '" + sc.name +
                               "' drains a gas space with no ship compartment behind it");
        if (sc.width <= 0) problems.push_back("scupper '" + sc.name + "' has no width");
        // A sill below the deck it drains would sit permanently under water and
        // pass the whole compartment overboard.
        const GasCompartment& g = gas[static_cast<std::size_t>(sc.gasCompartment)];
        if (sc.sillPos.z < g.floorZ - 1e-9)
            problems.push_back("scupper '" + sc.name + "' has its sill below the deck it drains");
    }
    return problems;
}

double Model::totalEntrainment(double atTime) const {
    double total = 0;
    for (const DesignFire& f : fires) {
        if (f.compartment < 0 || f.compartment >= static_cast<int>(gas.size())) continue;
        const GasCompartment& g = gas[static_cast<std::size_t>(f.compartment)];
        const Plume p{f.heatRelease(atTime), f.diameter, f.convectiveFraction};
        total += p.entrainment(g.interfaceZ() - f.baseZ);
    }
    return total;
}

StepResult Model::step(double dt, const Ship& ship, const Sea& sea) {
    StepResult out;
    out.time = time;
    if (dt <= 0) return out;

    double remaining = dt;
    double h = std::min(dt, maxSubstep);
    int budget = maxSubsteps;
    while (remaining > 1e-12 * dt && budget-- > 0) {
        h = std::min(h, remaining);
        if (substep(h, ship, sea, out)) {
            remaining -= h;
            ++out.substeps;
            // Creep back up, so one transient does not pin the step for the whole
            // run.
            //
            // 1.5 rather than 2.0, and **the reason once given for that does not
            // survive being measured.** This said doubling oscillates against the
            // rejection test and 1.5 settles; the ISO room at 500 kW takes 2485
            // substeps over 600 s at 1.5 and 2491 at 2.0, at 2 MW it takes 2533 and
            // 2525 -- 0.3% either way, and the ferry's engine rooms and the tall box
            // sit at the arithmetic floor of four per second under both. The factor
            // is a free choice between two that behave the same, not a tuned one.
            // `tests/test_fire.cpp` asserts the floor rather than the factor.
            h = std::min(h * 1.5, maxSubstep);
        } else {
            h *= 0.5;
        }
    }
    out.time = time;
    for (const DesignFire& f : fires) out.heatRelease += f.heatRelease(time);
    return out;
}

bool Model::substep(double dt, const Ship& ship, const Sea& sea, StepResult& out) {
    const std::size_t n = gas.size();
    // A caller that pushed a compartment on after `attach()` without resetting
    // the account would otherwise index past the end of this. `appliedMass_` has
    // the same shape and the same exposure; resizing here rather than asserting
    // keeps a hand-assembled model -- which is how every unit fixture in
    // `tests/test_fire.cpp` is built -- from being a special case.
    if (pendingWater_.size() != n) pendingWater_.resize(n, 0.0);
    Deltas d(n);
    // The non-vent energy rate per compartment: what the fire and the boundary do
    // to its pressure over this step whatever the vents carry. The implicit
    // pressure solve below needs it as the constant term of each equation.
    std::vector<double> sourceRate(n, 0.0);

    // --- The fires and their plumes ----------------------------------------
    for (const DesignFire& f : fires) {
        if (f.compartment < 0 || f.compartment >= static_cast<int>(n)) continue;
        const std::size_t i = static_cast<std::size_t>(f.compartment);
        GasCompartment& g = gas[i];
        // At the *midpoint* of the substep, not its start. The heat release is a
        // known function of time rather than a state variable, so integrating it
        // exactly costs nothing and is worth more than it looks: evaluated at the
        // left end, the first substep of a fire that starts at t = 0 releases
        // nothing at all, and the sealed compartment's closed form came out
        // exactly one substep of heat short -- 50 kJ, forever, on every run. The
        // midpoint rule is exact for a steady fire and second order on a t^2
        // growth curve.
        const double q = f.heatRelease(time + 0.5 * dt);
        if (q <= 0) continue;

        d.heat += q;
        d.radiative += q * f.radiativeLossFraction;
        const double toGas = q * (1.0 - f.radiativeLossFraction);
        sourceRate[i] += toGas;

        // Entrainment is evaluated at the interface: the plume drags air out of
        // the lower layer over exactly the height it rises through it. This is the
        // term that drives the layer down, and it dominates the fire's own thermal
        // expansion by about two orders of magnitude.
        const Plume plume{q, f.diameter, f.convectiveFraction};
        double mp = plume.entrainment(g.interfaceZ() - f.baseZ);
        // Taper as the lower layer is consumed: there is nothing left down there
        // to entrain. See `kEntrainmentTaperFraction`.
        const double lowerVolume = std::max(g.gasVolume - g.upperVolume(), 0.0);
        mp *= std::clamp(lowerVolume / (kEntrainmentTaperFraction * g.gasVolume), 0.0, 1.0);
        // And never more than half of what the layer holds in one step, so the
        // explicit update cannot take it negative whatever the taper does.
        mp = std::min(mp, 0.5 * std::max(g.lower.mass, 0.0) / dt);

        d.entrained += mp;
        const double tl = g.lower.temperature();
        const double yl = g.lower.productFraction();
        d.dmL[i] -= mp;
        d.dmU[i] += mp;
        d.dEL[i] -= mp * kCpAir * tl;
        d.dEU[i] += mp * kCpAir * tl + toGas;
        d.dsL[i] -= mp * yl;
        d.dsU[i] += mp * yl;

        // Products are a tracer riding the gas, not part of its mass: a 1 MW fire
        // burns about 0.05 kg/s of fuel against several kg/s through the vent, so
        // adding it to the gas mass would be noise dressed as physics -- and it
        // would put a term in the mass account with nothing on the other side.
        const double ms = f.productYield * q;
        d.dsU[i] += ms;
        d.productsGenerated += ms;
    }

    // --- Suppression: what the drenchers take out ---------------------------
    //
    // Before the vents for the same reason the boundary loss is: this is a
    // non-vent energy rate and the implicit pressure solve needs it as the
    // constant term of its equation. It is not a small one. A ro-ro drencher
    // sized to SOLAS puts 5 L/(min m^2) over the space, which over one 20 m
    // section of the ferry's vehicle deck is 31.1 kg/s and a nominal cooling
    // demand of 32.1 MW at the default evaporated fraction -- half as much again
    // as the 20 MW lorry fire it is aimed at. Leaving it out of the pressure
    // equation would let the vents carry what a compartment tens of megawatts
    // hotter would have carried.
    for (Drencher& s : drenchers) {
        s.lastCooling = s.lastEvaporated = s.lastToDeck = 0;
        if (!s.on || s.flow <= 0) continue;
        if (s.gasCompartment < 0 || s.gasCompartment >= static_cast<int>(n)) continue;
        const std::size_t i = static_cast<std::size_t>(s.gasCompartment);
        GasCompartment& g = gas[i];
        const double e = std::clamp(s.evaporatedFraction, 0.0, 1.0);

        // What one kilogram of spray takes out of the layer it falls through.
        //
        //   * **Sensible**: the drops leave at the layer's own temperature, or at
        //     saturation if the layer is hotter than that -- water cannot be
        //     heated past boiling and cannot be heated past the gas heating it.
        //   * **Latent**: `e L` on the evaporating share, ramped in over
        //     `kEvaporationBand` above saturation rather than switched on at it.
        //
        // Both are functions of the layer temperature and that is the whole point.
        // The first version made `perKg` a constant of the drencher and capped the
        // *engaged mass* at what the layer could supply, which is the same idiom
        // the entrainment cap uses -- and it was wrong here for a reason the
        // entrainment cap does not have. Entrainment's cap binds on a transient;
        // this one binds at the **steady state**, because a drencher routinely
        // out-absorbs the fire it is aimed at. A cap that binds every step forever
        // makes the applied rate `available / dt`, which is a function of the
        // substep and not of the physics: the model stopped converging under
        // refinement.
        //
        // Measured, on the fixture that shows it worst -- 20 kg/s into the ISO
        // room's 2 MW fire, a nominal demand ten times the release. Under the cap
        // the layer sat pinned at the water temperature and the applied rate was
        // whatever the step made it. Under the relaxation it settles at a steady
        // 333.4 K with 9.3% of the demand applied, and refining the substep 16x
        // moves the water landing on the deck by 0.6%.
        //
        // Written this way the sink vanishes as the layer reaches the water's own
        // temperature, so there is a real equilibrium -- `Q_fire = Q_spray(T)` --
        // and it is the same at any step size.
        const double tu = g.upper.temperature();
        const double excess = tu - s.waterTemperature;
        const double share = std::clamp((tu - kTSaturation) / kEvaporationBand, 0.0, 1.0);
        const double perKg =
            kCpWater * std::max(std::min(tu, kTSaturation) - s.waterTemperature, 0.0) +
            e * kLatentHeat * share;
        const double demand = s.flow * perKg;

        // Applied **exactly, not explicitly**, as an equivalent relaxation of the
        // layer towards the temperature of the water -- the same treatment, and
        // for the same reason, as the boundary loss below. The layer a spray falls
        // into can be arbitrarily thin: this model's seed layer is half a
        // millimetre and holds 0.12 kg, against which 30 kg/s of water is a time
        // constant of microseconds. Explicitly that is unusable at any step the
        // rest of the model wants to take.
        //
        // `-expm1(-x)` lies in [0, 1] and is below `x`, so this can neither take
        // the layer past the water's temperature at any dt, nor remove more than
        // the demand; and it tends to the demand as dt shrinks, which is what
        // makes the scheme converge.
        const double capacity = g.upper.mass * kCvAir;
        double cooling = 0;
        if (demand > 0 && excess > 0 && capacity > 0 && dt > 0) {
            const double kEff = demand / excess;   // W/K, the equivalent coefficient
            cooling = capacity * excess * -std::expm1(-kEff * dt / capacity) / dt;
        }

        // The steam leaves in the same proportion as the heat: a step that could
        // only take a fraction of the demand only boiled that fraction away.
        const double evaporated =
            demand > 0 ? e * share * s.flow * (cooling / demand) : 0.0;
        // Everything that did not leave as steam lands. This is the term the
        // stability half of the roadmap item is about, and it is the *whole* flow
        // when the compartment is cold.
        const double toDeck = s.flow - evaporated;

        s.lastCooling = cooling;
        s.lastEvaporated = evaporated;
        s.lastToDeck = toDeck;

        d.dEU[i] -= cooling;
        d.suppression += cooling;
        sourceRate[i] -= cooling;
        d.waterDelivered += s.flow;
        d.waterEvaporated += evaporated;
        d.dWater[i] += toDeck;
    }

    // --- Suppression: what the freeing ports get back out -------------------
    //
    // No gas term at all, so this sits outside the pressure solve: a scupper
    // moves water between the deck and the sea and the model's only interest in
    // it is the water. It is still accumulated as a rate into `Deltas` rather
    // than applied here, so that a rejected substep leaves no trace of it.
    for (Scupper& sc : scuppers) {
        sc.lastFlow = sc.lastInsideHead = sc.lastOutsideHead = 0;
        if (sc.blocked || sc.width <= 0 || sc.dischargeCoeff <= 0) continue;
        if (sc.gasCompartment < 0 || sc.gasCompartment >= static_cast<int>(n)) continue;
        const std::size_t i = static_cast<std::size_t>(sc.gasCompartment);
        const int shipIdx = gas[i].shipCompartment;
        if (shipIdx == kSea || shipIdx >= static_cast<int>(ship.compartments.size())) continue;
        const Compartment& c = ship.compartments[static_cast<std::size_t>(shipIdx)];

        // Both heads are true vertical depths in world coordinates, taken at the
        // port's own position -- the same comparison `Ship::sideStateAt` makes
        // when it decides whether an opening has water against it. On a lolled
        // ship the low-side ports carry a head the high-side ones do not, and
        // that asymmetry is the physics, not an artefact.
        const Vec3 sillWorld = ship.state.orientation.toMat3() * sc.sillPos + ship.state.position;
        // A dry compartment parks its free surface a metre below its own floor, so
        // `surfaceWorldZ` only means anything when there is water to have a
        // surface.
        const double insideHead =
            c.waterVolume > 1e-9 ? std::max(c.surfaceWorldZ - sillWorld.z, 0.0) : 0.0;
        const double outsideHead =
            std::max(sea.heightAt(sillWorld.x, sillWorld.y) - sillWorld.z, 0.0);
        sc.lastInsideHead = insideHead;
        sc.lastOutsideHead = outsideHead;

        double q = sc.dischargeCoeff * sc.width * scupperFlow(insideHead, outsideHead);
        if (q == 0) continue;

        // The ship's own water level is a step behind -- `step()` reads the ship
        // const and may take many substeps against one snapshot of it -- so the
        // rate is capped by what is actually there to move. At any sane tick the
        // correction is nothing: a 2 m port under 0.1 m of head passes 0.056 m^3/s
        // against a deck holding hundreds. It is here for the pathological caller,
        // and because a port that drains water the ship does not have would show
        // up as a broken water account rather than as anything visible.
        const double held = c.waterVolume + pendingWater_[i];
        if (q > 0) q = std::min(q, std::max(held, 0.0) / std::max(dt, 1e-12));
        else q = -std::min(-q, std::max(c.floodableVolume() - held, 0.0) / std::max(dt, 1e-12));

        sc.lastFlow = q;
        const double kg = q * ship.seaDensity;
        d.dWater[i] -= kg;
        d.waterDrained += kg;
    }

    // --- Boundary heat loss -------------------------------------------------
    //
    // Before the vents, not after: the implicit pressure solve needs every
    // non-vent energy rate as the constant term of its equation.
    for (std::size_t i = 0; i < n; ++i) {
        GasCompartment& g = gas[i];
        if (g.wallConductance <= 0) continue;
        const double zi = g.interfaceZ();
        // Ceiling plus the walls the hot layer wets; floor plus the walls the cool
        // one does. Two zones exist precisely so that this is not one area at one
        // temperature -- charging the whole enclosure at the upper-layer
        // temperature, which is what MQH's own `A_T` does, costs about a quarter
        // of the steady temperature rise.
        const double aUpper = g.floorArea + g.perimeter * (g.ceilingZ - zi);
        const double aLower = g.floorArea + g.perimeter * (zi - g.floorZ);
        // **Exactly, not explicitly.** In isolation this term is a linear
        // relaxation towards the boundary temperature, `m c_v dT/dt = -h A dT`,
        // whose solution over a step is an exponential; taking the average rate
        // from that solution costs one `expm1` and is unconditionally stable.
        //
        // Explicit was tried and it is the stiffest term in the model, by a wide
        // margin and from an unexpected direction. A layer of *negligible
        // thickness* still wets the whole deckhead: the seed upper layer this
        // model starts with is half a millimetre thick and holds 0.12 kg, while
        // the ceiling it touches is 188 m^2. Its time constant is milliseconds.
        // Explicitly, it oscillated -- 288 K, 300 K, 177 K, 1260 K on successive
        // substeps -- and each excursion through zero energy tripped the
        // non-negativity clamp, which is the only thing in this file that can
        // break the energy account. It did, at 0.09%.
        auto relax = [&](double area, double temp, double mass) {
            const double capacity = mass * kCvAir;
            const double rate = g.wallConductance * area;
            if (capacity <= 0 || rate <= 0) return 0.0;
            const double removed = capacity * (temp - g.wallTemperature) *
                                   -std::expm1(-rate * dt / capacity);
            return removed / dt;
        };
        const double qu = relax(aUpper, g.upper.temperature(), g.upper.mass);
        const double ql = relax(aLower, g.lower.temperature(), g.lower.mass);
        d.dEU[i] -= qu;
        d.dEL[i] -= ql;
        d.wall += qu + ql;
        sourceRate[i] -= qu + ql;
    }

    // --- The opening network, solved implicitly in the pressures ------------
    //
    // `Ship::solveFlowNetwork` moves air explicitly and clamps the transfer to
    // the mass that equalises the two pressures, because gas is the stiffest
    // thing in that network and an explicit step that ignores the target state
    // blows past equilibrium and rings.
    //
    // **That clamp does not carry over, and it was tried.** It fails twice:
    //
    //   * *Equal pressures is not equilibrium* once one side has a hot layer. A
    //     doorway with equal floor pressures is exchanging kilograms a second in
    //     both directions at once. The invariant is zero *net* mass flow.
    //   * *A fire adds energy during the step.* Clamping the transfer to the
    //     imbalance that exists at the start of the step ignores the pressure the
    //     fire will add over the same dt. Measured: a 500 kW fire in an ISO room
    //     held 800 Pa above atmosphere, venting a twentieth of what it should,
    //     with the smoke layer driven onto the floor.
    //
    // Repairing the target by adding the source term fixes the pressure runaway
    // and then rings anyway, because scaling one direction of a bidirectional
    // flow to hit a net target is not a continuous function of the state: the
    // measured steady state was a limit cycle alternating between -3 Pa and
    // -590 Pa on successive substeps.
    //
    // So the clamp is not repaired here, it is replaced by the thing it was
    // approximating: **the compartment pressures are solved implicitly.** Freeze
    // the slow variables -- layer temperatures, the interface, the product
    // loadings -- and find the end-of-step gauge pressures that satisfy
    //
    //     x_i = gauge_i + (gamma-1) dt E_i(x) / V_i
    //
    // where `E_i` is the net energy rate into compartment i. `E_i` is monotone
    // decreasing in `x_i` -- raising a compartment's pressure can only increase
    // what leaves and reduce what enters -- so each compartment's equation is a
    // safe 1-D bisection, and Gauss-Seidel over the few compartments closes the
    // coupling. A sealed compartment has no `x` dependence at all and the
    // solution is exact in one evaluation, which is what keeps
    // `p = p0 + (gamma-1) E / V` exact for the sealed case.
    //
    // Backward Euler in the fast variable, forward in the slow ones: stable at
    // any step, and the step is then set by accuracy rather than by the 0.65 ms
    // an explicit doorway would demand.
    std::vector<VentSide> sides(n);
    for (std::size_t i = 0; i < n; ++i) sides[i] = sideOf(gas[i]);
    const VentSide outside = ambientSide();

    // The vents that are actually in play this step, already clipped to the gas
    // space on both sides. Collected once: the pressure solve evaluates them
    // tens of times each.
    std::vector<Vent> active;
    std::vector<int> activeOwner;
    for (std::size_t vi = 0; vi < vents.size(); ++vi) {
        Vent& v = vents[vi];
        v.massAToB = v.massBToA = 0;
        v.bidirectional = false;
        v.blockedByWater = false;
        v.activeSillZ = v.activeSoffitZ = 0;
        if (!v.open) continue;

        // Water against either side means this opening belongs to the flooding
        // solve this tick, not to the gas model. Splitting them by phase is the
        // same rule `Ship::solveFlowNetwork` already applies, and it is what stops
        // the two models moving the same mass twice.
        if (v.opening >= 0) {
            const Opening& o = ship.openings[static_cast<std::size_t>(v.opening)];
            const int sca = v.a == kSea ? kSea : gas[static_cast<std::size_t>(v.a)].shipCompartment;
            const int scb = v.b == kSea ? kSea : gas[static_cast<std::size_t>(v.b)].shipCompartment;
            if (waterAgainst(ship, sea, sca, o.pos) || waterAgainst(ship, sea, scb, o.pos)) {
                v.blockedByWater = true;
                continue;
            }
        }

        // Bring the span into the gas space on both sides.
        //
        // **Shifted before it is trimmed**, and that is not a nicety. The ferry
        // authors its air escapes at the gooseneck -- `vent_er_s` sits at
        // z = 12.5 m, on the weather deck -- while the engine room it drains tops
        // out at 7.0 m. `Opening::pos` for a pipe is the *outboard* end, because
        // that is the end whose height decides whether the sea is over it, which
        // is all the flooding solve needs. Trimming alone would delete every vent
        // on the ship and leave the machinery spaces sealed.
        //
        // So the vent is slid, at its own height, to the nearest end of the space
        // it opens into: a pipe from a compartment terminates at the deckhead.
        // **What this does not model is the pipe itself.** Five and a half metres
        // of duct between the deckhead and the gooseneck is a chimney, and its
        // stack effect on hot gas is a real driver that is simply not here.
        double lo = -kInf, hi = kInf;
        for (int side = 0; side < 2; ++side) {
            const int gi = side == 0 ? v.a : v.b;
            if (gi == kSea) continue;
            const GasCompartment& g = gas[static_cast<std::size_t>(gi)];
            lo = std::max(lo, g.floorZ);
            hi = std::min(hi, g.ceilingZ);
        }
        if (hi <= lo) continue;

        Vent clipped = v;
        if (clipped.soffitZ > hi) {
            const double shift = clipped.soffitZ - hi;
            clipped.soffitZ -= shift;
            clipped.sillZ -= shift;
        }
        if (clipped.sillZ < lo) {
            const double shift = lo - clipped.sillZ;
            clipped.sillZ += shift;
            clipped.soffitZ += shift;
        }
        clipped.sillZ = std::max(clipped.sillZ, lo);
        clipped.soffitZ = std::min(clipped.soffitZ, hi);
        if (!clipped.horizontal && clipped.soffitZ <= clipped.sillZ) continue;
        v.activeSillZ = clipped.sillZ;
        v.activeSoffitZ = clipped.soffitZ;
        active.push_back(clipped);
        activeOwner.push_back(static_cast<int>(vi));
    }

    // The unknowns are the end-of-step gauge pressures. Start each compartment
    // where its *own* sources alone would put it -- which is already the exact
    // answer for a sealed space, and is what keeps `p0 + (gamma-1)E/V` exact.
    std::vector<double> x(n, 0.0), cap(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        cap[i] = (kGammaAir - 1.0) * dt / std::max(gas[i].gasVolume, 1e-12);
        x[i] = sides[i].gaugeAtFloor + cap[i] * sourceRate[i];
    }

    auto sideAt = [&](int gi, double gauge) -> VentSide {
        if (gi == kSea) return outside;
        VentSide s = sides[static_cast<std::size_t>(gi)];
        s.gaugeAtFloor = gauge;
        return s;
    };

    // **Blocked by vent, not by compartment**, and that distinction is the whole
    // difference between this converging and not.
    //
    // Sweeping compartment by compartment -- solve A holding B, then B holding A
    // -- is the obvious arrangement and it fails outright whenever a vent is
    // large relative to the spaces it joins. Measured on two 60 m^3 boxes sharing
    // a 2 m^2 doorway with a 200 Pa imbalance: each sweep closed **0.037 Pa** of
    // it, because the orifice law is so flat near zero pressure difference that
    // the iteration is sublinear, not merely slow. Twelve sweeps left the boxes
    // 119 Pa apart and the test that found it saw two compartments that would not
    // equalise in forty seconds through a doorway.
    //
    // Solving one *vent* at a time fixes it, because the unknown is then the
    // energy the vent carries and both of its endpoints move together: the block
    // is exactly the pair, and a single vent is solved in one sweep. `tv[k]` is
    // what vent k is currently contributing, so removing and re-adding it keeps
    // every compartment's pressure the sum of its sources and its vents.
    std::vector<double> tv(active.size(), 0.0);   // W into side A of each vent

    for (int sweep = 0; sweep < kPressureSweeps; ++sweep) {
        double moved = 0;
        for (std::size_t k = 0; k < active.size(); ++k) {
            const Vent& v = active[k];
            const double ca = v.a == kSea ? 0.0 : cap[static_cast<std::size_t>(v.a)];
            const double cb = v.b == kSea ? 0.0 : cap[static_cast<std::size_t>(v.b)];
            if (ca + cb <= 0) continue;

            // Both endpoints with this vent's own contribution taken back out.
            const double xa0 =
                v.a == kSea ? 0.0 : x[static_cast<std::size_t>(v.a)] - ca * tv[k];
            const double xb0 =
                v.b == kSea ? 0.0 : x[static_cast<std::size_t>(v.b)] + cb * tv[k];

            // `t - carried(t)` is monotone increasing: raising the vent's delivery
            // to A raises A's pressure and lowers B's, which can only reduce what
            // the vent then carries. Bisection is therefore safe and needs no
            // derivative of an orifice law whose derivative is infinite at the
            // neutral plane.
            auto residual = [&](double t) {
                const VentResult r =
                    ventMassFlow(v, sideAt(v.a, xa0 + ca * t), sideAt(v.b, xb0 - cb * t));
                return t - (r.enthalpyBToA - r.enthalpyAToB);
            };

            const double r0 = residual(0.0);
            double lo = 0.0, hi = 0.0, t = 0.0;
            if (r0 != 0) {
                double span = std::max(std::abs(r0), 1e-9);
                int guard = 0;
                if (r0 < 0) {
                    hi = span;
                    while (residual(hi) < 0 && guard++ < 200) { span *= 2.0; hi = span; }
                } else {
                    lo = -span;
                    while (residual(lo) > 0 && guard++ < 200) { span *= 2.0; lo = -span; }
                }
                if (guard >= 200) out.pressureSolveCapped = true;
                for (int it = 0; it < 200; ++it) {
                    const double mid = 0.5 * (lo + hi);
                    if (hi - lo <= 1e-9 * std::max(1.0, std::abs(mid))) break;
                    if (residual(mid) < 0) lo = mid; else hi = mid;
                }
                t = 0.5 * (lo + hi);
            }

            moved = std::max(moved, std::abs(t - tv[k]) * (ca + cb));
            if (v.a != kSea) x[static_cast<std::size_t>(v.a)] = xa0 + ca * t;
            if (v.b != kSea) x[static_cast<std::size_t>(v.b)] = xb0 - cb * t;
            tv[k] = t;
        }
        out.pressureSweeps = sweep + 1;
        if (moved < 1e-9) break;
    }

    // Evaluate each vent once at the solved pressures and move what it carries.
    for (std::size_t k = 0; k < active.size(); ++k) {
        const Vent& cv = active[k];
        const VentSide sa = sideAt(cv.a, cv.a == kSea ? 0.0 : x[static_cast<std::size_t>(cv.a)]);
        const VentSide sb = sideAt(cv.b, cv.b == kSea ? 0.0 : x[static_cast<std::size_t>(cv.b)]);
        const VentResult res = ventMassFlow(cv, sa, sb);

        Vent& v = vents[static_cast<std::size_t>(activeOwner[k])];
        v.massAToB = res.massAToB;
        v.massBToA = res.massBToA;
        v.bidirectional = res.bidirectional;
        v.neutralPlaneZ = res.neutralPlaneZ;
        if (res.massAToB <= 0 && res.massBToA <= 0) continue;

        applyStream(cv.a, cv.b, sa, sb, res.fromUpperA, res.massAToB - res.fromUpperA, d);
        applyStream(cv.b, cv.a, sb, sa, res.fromUpperB, res.massBToA - res.fromUpperB, d);
    }

    // --- Accept or reject ---------------------------------------------------
    //
    // Nothing above has touched the state or the account, so a rejection here is
    // free and leaves no trace. Below 1 ns there is nothing left to subdivide and
    // the step is taken regardless, which is what stops the halving from spinning
    // on a genuinely singular configuration.
    if (dt > 1e-9) {
        for (std::size_t i = 0; i < n; ++i) {
            const GasCompartment& g = gas[i];
            const double mScale = std::max(g.totalMass(), kMassFloor);
            const double eScale = std::max(g.totalEnergy(), 1.0);
            const double dm = std::abs(d.dmU[i]) + std::abs(d.dmL[i]);
            const double de = std::abs(d.dEU[i]) + std::abs(d.dEL[i]);
            if (dm * dt > maxRelativeChange * mScale) return false;
            if (de * dt > maxRelativeChange * eScale) return false;
            // A layer taken negative would be caught by the clamps below, and
            // those clamps are the only thing here that can put a hole in the
            // energy account -- they conjure whatever they clamp away. Rejecting
            // the step instead keeps the account exact.
            double dU = 0, dL = 0;
            layerSplit(g.upper.energy, g.totalEnergy(), d.dEU[i], d.dEL[i], dU, dL);
            if (g.upper.energy + dU * dt < 0 || g.lower.energy + dL * dt < 0) return false;
            if (g.upper.mass + d.dmU[i] * dt < 0 || g.lower.mass + d.dmL[i] * dt < 0) return false;
        }
    }

    // --- Commit -------------------------------------------------------------
    for (std::size_t i = 0; i < n; ++i) {
        GasCompartment& g = gas[i];
        double dU = 0, dL = 0;
        layerSplit(g.upper.energy, g.totalEnergy(), d.dEU[i], d.dEL[i], dU, dL);
        g.upper.energy = std::max(g.upper.energy + dU * dt, 0.0);
        g.lower.energy = std::max(g.lower.energy + dL * dt, 0.0);
        g.upper.mass = std::max(g.upper.mass + d.dmU[i] * dt, 0.0);
        g.lower.mass = std::max(g.lower.mass + d.dmL[i] * dt, 0.0);
        g.upper.products = std::max(g.upper.products + d.dsU[i] * dt, 0.0);
        g.lower.products = std::max(g.lower.products + d.dsL[i] * dt, 0.0);
    }

    // The water the drenchers landed and the ports did not take back, converted
    // to volume at the density `Ship` weighs floodwater with -- so that the mass
    // this model says it delivered is the mass the ship's displacement gains.
    for (std::size_t i = 0; i < n; ++i)
        if (d.dWater[i] != 0.0) pendingWater_[i] += d.dWater[i] * dt / ship.seaDensity;

    account.heatReleased += d.heat * dt;
    account.radiativeLoss += d.radiative * dt;
    account.wallLoss += d.wall * dt;
    account.suppressionCooling += d.suppression * dt;
    account.waterDelivered += d.waterDelivered * dt;
    account.waterEvaporated += d.waterEvaporated * dt;
    account.waterDrained += d.waterDrained * dt;
    account.enthalpyIn += d.enthalpyIn * dt;
    account.enthalpyOut += d.enthalpyOut * dt;
    account.massIn += d.massIn * dt;
    account.massOut += d.massOut * dt;
    account.productsGenerated += d.productsGenerated * dt;
    account.productsOut += d.productsOut * dt;
    out.entrainment = d.entrained;
    out.suppressionCooling = d.suppression;

    time += dt;
    account.energy = 0;
    account.mass = 0;
    account.products = 0;
    for (const GasCompartment& g : gas) {
        account.energy += g.totalEnergy();
        account.mass += g.totalMass();
        account.products += g.upper.products + g.lower.products;
    }
    return true;
}

void Model::applyTo(Ship& ship) {
    if (pendingWater_.size() != gas.size()) pendingWater_.resize(gas.size(), 0.0);
    for (std::size_t i = 0; i < gas.size(); ++i) {
        const GasCompartment& g = gas[i];
        if (g.shipCompartment == kSea) continue;
        Compartment& c = ship.compartments[static_cast<std::size_t>(g.shipCompartment)];
        // The mass that reproduces this model's pressure under the ship's own
        // isothermal formula. See fire.hpp for why this is a pressure proxy rather
        // than a mass, and why the write is a delta.
        const double mEq = g.pressure() * g.gasVolume / (kRAir * kTAmbient);
        c.airMass += mEq - appliedMass_[i];
        appliedMass_[i] = mEq;

        // Suppression water, handed over to the flooding solve and forgotten.
        // From here the ship owns it: it re-levels against gravity, it moves the
        // centre of gravity, it makes a free surface, and it may leave again
        // through the ship's own openings. Nothing in this file models any of
        // that, because `Ship` already does and doing it twice would be two
        // answers to one question.
        //
        // **Guarded rather than written unconditionally.** `x += 0.0` is not
        // quite the identity -- it turns a negative zero positive -- and
        // `Compartment::waterVolume` is clamped by a `std::clamp` that returns
        // its argument unchanged when it compares equal to the bound, so a
        // negative zero is reachable. The exact control this file is under says
        // *bit*-identical, so the store is skipped rather than made harmless.
        //
        // Mutation testing says the suite **cannot currently tell**: replacing
        // this condition with `true` survives every test, because no fixture ever
        // has a compartment at a negative zero or outside its own bounds when
        // nothing is owed. It stays anyway. The condition costs one comparison,
        // and the alternative is a write and a clamp on every tracked compartment
        // on every call for a model that is doing nothing -- which is the shape of
        // thing the exact control exists to forbid, whether or not today's
        // fixtures can catch it.
        //
        // What will not fit stays **owed** rather than being dropped. A
        // compartment flooded solid still has a drencher running into it, and
        // discarding the overflow here would put a hole in the water account of
        // exactly the kind the layer clamps put in the energy one. Kept, the
        // invariant `written + owed == waterDelivered - evaporated - drained`
        // holds for the model's whole life, which is what
        // `tests/test_fire.cpp` asserts.
        if (pendingWater_[i] != 0.0) {
            const double want = c.waterVolume + pendingWater_[i];
            const double got = std::clamp(want, 0.0, c.floodableVolume());
            c.waterVolume = got;
            pendingWater_[i] = want - got;
        }
    }
}

// ---------------------------------------------------------------------------
// Enclosure correlations
// ---------------------------------------------------------------------------

double mqhTemperatureRise(double heatRelease, double ventArea, double ventHeight,
                          double wallConductance, double wallArea) {
    if (heatRelease <= 0 || ventArea <= 0 || ventHeight <= 0 || wallConductance <= 0 ||
        wallArea <= 0)
        return 0.0;
    const double qkw = heatRelease / kWattsPerKilowatt;
    const double hk = wallConductance / kWattsPerKilowatt;   // kW/(m^2 K)
    const double denom = ventArea * std::sqrt(ventHeight) * hk * wallArea;
    return 6.85 * std::cbrt(qkw * qkw / denom);
}

double thomasFlashoverPower(double wallArea, double ventArea, double ventHeight) {
    if (wallArea <= 0 || ventArea <= 0 || ventHeight <= 0) return 0.0;
    return kWattsPerKilowatt * (7.8 * wallArea + 378.0 * ventArea * std::sqrt(ventHeight));
}

}  // namespace sim::fire
