// SPDX-License-Identifier: MIT
#include "ship.hpp"

#include <cmath>
#include <limits>
#include <string_view>

namespace sim {
namespace {

constexpr double kInf = std::numeric_limits<double>::infinity();

// Density of the reference atmosphere: the air outside the hull, and the air in
// every compartment that has not been given a reason to be anything else.
//
// **Written as this expression and not as a number.** Everything below that
// carries a gas temperature subtracts this from `kPatm / (kRAir * T)`, and at
// T == kTAmbient the two operands are then the same double bit for bit, so the
// difference is exactly +0.0 and the whole buoyancy term drops out by IEEE
// arithmetic rather than by being small. That is what makes a cold ship
// bit-identical to the model that had no temperature in it.
constexpr double kRhoAirAmbient = kPatm / (kRAir * kTAmbient);

// Gas density at `kelvin`, at the reference pressure.
//
// Deliberately *not* p/(R T) with the compartment's own pressure. This value is
// used only for the buoyancy head between two spaces, and there it is the
// thermal part of the density difference that matters -- a fire takes the gas
// from 288 K to 800 K and halves its density, while pressurising a compartment
// enough to flood it moves the density by around a percent. Taking the pressure
// term as well would make the head non-zero for two cold compartments at
// different pressures, which is a real effect the flooding model has never had
// and is not the effect being added here.
inline double gasBuoyancyDensity(double kelvin) { return kPatm / (kRAir * kelvin); }

// The head a compartment's gas carries between the datum its pressure is quoted
// at -- the middle of the gas space -- and some other height in the same space.
//
// The gas has weight, and until now this model did not admit it: `airPressure`
// was returned at every height in the compartment, so two spaces joined by a high
// opening and a low one saw the *same* pressure difference at both and moved gas
// the same way through each. There is no neutral plane in that, and therefore no
// buoyant exchange -- which is the entire mechanism by which smoke leaves a
// burning compartment. Adding a temperature without adding this would have bought
// a hotter gas that still went nowhere.
//
// The head is taken against the reference atmosphere rather than absolutely.
// Absolutely, two *cold* compartments whose gas centroids sit at different
// heights would exchange air through a shared opening for ever -- a real
// stratification, but one this flooding model has never carried and not the
// effect being added. Against the reference it is the thermal buoyancy that
// survives, and at ambient the coefficient is exactly +0.0, so a cold ship gets
// `airPressure` back unchanged, bit for bit.
inline double gasBuoyancyHead(const Compartment& c, double worldZ) {
    return -(gasBuoyancyDensity(c.gasTemperature) - kRhoAirAmbient) * kGravity *
           (worldZ - c.gasCentroidWorldZ);
}

// The diameter of the circle with the same area as an orifice.
//
// The buoyancy exchange below needs a *length* and this is the only one an
// `Opening` carries: it has an area and a centre and no height. Taken from the
// authored area rather than from the projection onto the horizontal, because it
// stands for the physical size of the hole -- the scale of the overturning that
// the heavy fluid has to organise itself into to get through -- and a hole does
// not get smaller when the ship heels. What heeling takes away is the area
// available to a vertical counterflow, which is a separate factor.
inline double hydraulicDiameter(double area) { return std::sqrt(4.0 * area / kPi); }

// Ship "up" expressed in the body frame. Every free surface -- the sea outside and
// the floodwater inside -- is a plane with this normal.
inline Vec3 bodyFrameUp(const Mat3& R) { return R.transposed() * Vec3{0, 0, 1}; }

// Highest and lowest points of an axis-aligned box in the `up` direction. Exact
// for the box-carved compartments this engine builds, and a bound for anything
// else.
inline double boxTopAlong(const Vec3& up, const Vec3& lo, const Vec3& hi) {
    return (up.x > 0 ? up.x * hi.x : up.x * lo.x) + (up.y > 0 ? up.y * hi.y : up.y * lo.y) +
           (up.z > 0 ? up.z * hi.z : up.z * lo.z);
}
inline double boxBottomAlong(const Vec3& up, const Vec3& lo, const Vec3& hi) {
    return (up.x > 0 ? up.x * lo.x : up.x * hi.x) + (up.y > 0 ? up.y * lo.y : up.y * hi.y) +
           (up.z > 0 ? up.z * lo.z : up.z * hi.z);
}

void boundingBox(const TriMesh& m, Vec3& lo, Vec3& hi) {
    lo = {kInf, kInf, kInf};
    hi = {-kInf, -kInf, -kInf};
    for (const Vec3& v : m.verts) {
        lo = {std::min(lo.x, v.x), std::min(lo.y, v.y), std::min(lo.z, v.z)};
        hi = {std::max(hi.x, v.x), std::max(hi.y, v.y), std::max(hi.z, v.z)};
    }
}

// Inertia of a homogeneous box of mass m and side lengths e, about its own centre.
Mat3 boxInertia(double m, const Vec3& e) {
    Mat3 I = Mat3::zero();
    I(0, 0) = m * (e.y * e.y + e.z * e.z) / 12.0;
    I(1, 1) = m * (e.x * e.x + e.z * e.z) / 12.0;
    I(2, 2) = m * (e.x * e.x + e.y * e.y) / 12.0;
    return I;
}

}  // namespace

// ---------------------------------------------------------------------------
// Setup
// ---------------------------------------------------------------------------

Ship::~Ship() {
    for (auto& [idx, field] : activeWaterFields_) {
        delete field;
    }
}

Ship::Ship(const Ship& other)
    : hull(other.hull)
    , compartments(other.compartments)
    , openings(other.openings)
    , pumps(other.pumps)
    , lightshipMass(other.lightshipMass)
    , lightshipCog(other.lightshipCog)
    , gyradii(other.gyradii)
    , seaDensity(other.seaDensity)
    , deckEdgeZ(other.deckEdgeZ)
    , hullLo(other.hullLo)
    , hullHi(other.hullHi)
    , zetaHeave(other.zetaHeave)
    , zetaRoll(other.zetaRoll)
    , zetaPitch(other.zetaPitch)
    , addedMassSurge(other.addedMassSurge)
    , addedMassSway(other.addedMassSway)
    , addedMassHeave(other.addedMassHeave)
    , addedInertiaRoll(other.addedInertiaRoll)
    , addedInertiaPitch(other.addedInertiaPitch)
    , addedInertiaYaw(other.addedInertiaYaw)
    , radiation(other.radiation)
    , propulsion(other.propulsion)
    , rollDampingForm(other.rollDampingForm)
    , rollCondition(other.rollCondition)
    , rollDampingApplied(other.rollDampingApplied)
    , externalForce(other.externalForce)
    , externalMoment(other.externalMoment)
    , state(other.state)
    , cachedWaterplaneArea_(other.cachedWaterplaneArea_)
    , cachedKRoll_(other.cachedKRoll_)
    , cachedKPitch_(other.cachedKPitch_)
    , stabilityRefreshCounter_(other.stabilityRefreshCounter_)
{
    // `activeWaterFields_` is left empty and `waterPromoter_` is left
    // default-constructed, deliberately and together. See ship.hpp: a copy is a
    // ship with no promoted water, and copying the promoter while dropping the
    // fields would leave the two disagreeing about what is active.
}

Ship& Ship::operator=(const Ship& other) {
    if (this != &other) {
        // Clean up existing fields
        for (auto& [idx, field] : activeWaterFields_) {
            delete field;
        }
        activeWaterFields_.clear();

        hull = other.hull;
        compartments = other.compartments;
        openings = other.openings;
        pumps = other.pumps;
        lightshipMass = other.lightshipMass;
        lightshipCog = other.lightshipCog;
        gyradii = other.gyradii;
        seaDensity = other.seaDensity;
        deckEdgeZ = other.deckEdgeZ;
        hullLo = other.hullLo;
        hullHi = other.hullHi;
        zetaHeave = other.zetaHeave;
        zetaRoll = other.zetaRoll;
        zetaPitch = other.zetaPitch;
        addedMassSurge = other.addedMassSurge;
        addedMassSway = other.addedMassSway;
        addedMassHeave = other.addedMassHeave;
        addedInertiaRoll = other.addedInertiaRoll;
        addedInertiaPitch = other.addedInertiaPitch;
        addedInertiaYaw = other.addedInertiaYaw;
        radiation = other.radiation;
        propulsion = other.propulsion;
        rollDampingForm = other.rollDampingForm;
        rollCondition = other.rollCondition;
        rollDampingApplied = other.rollDampingApplied;
        externalForce = other.externalForce;
        externalMoment = other.externalMoment;
        state = other.state;
        // Cleared, not copied, and paired with the empty field map above: see
        // the copy constructor and ship.hpp.
        waterPromoter_.clear();
        cachedWaterplaneArea_ = other.cachedWaterplaneArea_;
        cachedKRoll_ = other.cachedKRoll_;
        cachedKPitch_ = other.cachedKPitch_;
        stabilityRefreshCounter_ = other.stabilityRefreshCounter_;
    }
    return *this;
}

Ship::Ship(Ship&& other) noexcept
    : hull(std::move(other.hull))
    , compartments(std::move(other.compartments))
    , openings(std::move(other.openings))
    , pumps(std::move(other.pumps))
    , lightshipMass(other.lightshipMass)
    , lightshipCog(other.lightshipCog)
    , gyradii(other.gyradii)
    , seaDensity(other.seaDensity)
    , deckEdgeZ(other.deckEdgeZ)
    , hullLo(other.hullLo)
    , hullHi(other.hullHi)
    , zetaHeave(other.zetaHeave)
    , zetaRoll(other.zetaRoll)
    , zetaPitch(other.zetaPitch)
    , addedMassSurge(other.addedMassSurge)
    , addedMassSway(other.addedMassSway)
    , addedMassHeave(other.addedMassHeave)
    , addedInertiaRoll(other.addedInertiaRoll)
    , addedInertiaPitch(other.addedInertiaPitch)
    , addedInertiaYaw(other.addedInertiaYaw)
    , radiation(std::move(other.radiation))
    , propulsion(std::move(other.propulsion))
    , rollDampingForm(other.rollDampingForm)
    , rollCondition(other.rollCondition)
    , rollDampingApplied(other.rollDampingApplied)
    , externalForce(other.externalForce)
    , externalMoment(other.externalMoment)
    , state(other.state)
    , waterPromoter_(std::move(other.waterPromoter_))
    , activeWaterFields_(std::move(other.activeWaterFields_))
    , cachedWaterplaneArea_(other.cachedWaterplaneArea_)
    , cachedKRoll_(other.cachedKRoll_)
    , cachedKPitch_(other.cachedKPitch_)
    , stabilityRefreshCounter_(other.stabilityRefreshCounter_)
{
}

Ship& Ship::operator=(Ship&& other) noexcept {
    if (this != &other) {
        // Clean up existing fields
        for (auto& [idx, field] : activeWaterFields_) {
            delete field;
        }

        hull = std::move(other.hull);
        compartments = std::move(other.compartments);
        openings = std::move(other.openings);
        pumps = std::move(other.pumps);
        lightshipMass = other.lightshipMass;
        lightshipCog = other.lightshipCog;
        gyradii = other.gyradii;
        seaDensity = other.seaDensity;
        deckEdgeZ = other.deckEdgeZ;
        hullLo = other.hullLo;
        hullHi = other.hullHi;
        zetaHeave = other.zetaHeave;
        zetaRoll = other.zetaRoll;
        zetaPitch = other.zetaPitch;
        addedMassSurge = other.addedMassSurge;
        addedMassSway = other.addedMassSway;
        addedMassHeave = other.addedMassHeave;
        addedInertiaRoll = other.addedInertiaRoll;
        addedInertiaPitch = other.addedInertiaPitch;
        addedInertiaYaw = other.addedInertiaYaw;
        radiation = std::move(other.radiation);
        propulsion = std::move(other.propulsion);
        rollDampingForm = other.rollDampingForm;
        rollCondition = other.rollCondition;
        rollDampingApplied = other.rollDampingApplied;
        externalForce = other.externalForce;
        externalMoment = other.externalMoment;
        state = other.state;
        waterPromoter_ = std::move(other.waterPromoter_);
        activeWaterFields_ = std::move(other.activeWaterFields_);
        cachedWaterplaneArea_ = other.cachedWaterplaneArea_;
        cachedKRoll_ = other.cachedKRoll_;
        cachedKPitch_ = other.cachedKPitch_;
        stabilityRefreshCounter_ = other.stabilityRefreshCounter_;
    }
    return *this;
}

void Ship::initialise(const Sea& sea) {
    boundingBox(hull, hullLo, hullHi);

    for (Compartment& c : compartments) {
        c.grossVolume = integrate(c.mesh).volume;
        boundingBox(c.mesh, c.bboxLo, c.bboxHi);
        c.waterVolume = std::clamp(c.waterVolume, 0.0, c.floodableVolume());
        // Start dry compartments full of air at atmospheric pressure, at whatever
        // temperature the caller asked for -- a compartment that opens already
        // hot holds correspondingly less air.
        c.airMass = kPatm * std::max(c.airVolume(), 0.0) / (kRAir * c.gasTemperature);
        c.airPressure = kPatm;
        c.lastAirVolume = std::max(c.airVolume(), 0.0);
    }

    // Drop the hull to its floating waterline before the first tick, so the sim
    // does not open with a large transient heave.
    state.position.z = 0;
    updateInternalFreeSurfaces(sea);
    const MassProperties mp = massProperties();
    state.position.z = equilibriumDraftAt(state.orientation, sea, mp.mass);
    updateInternalFreeSurfaces(sea);
}

std::vector<std::string> Ship::validate() const {
    std::vector<std::string> problems;
    auto note = [&](const std::string& s) { problems.push_back(s); };

    if (!isClosedManifold(hull)) note("hull is not a closed, consistently wound manifold");
    const double hullVolume = integrate(hull).volume;
    if (hullVolume <= 0) note("hull volume is not positive");

    double subdivided = 0;
    for (const Compartment& c : compartments) {
        if (c.mesh.tris.empty()) { note("compartment '" + c.name + "' is empty"); continue; }
        if (!isClosedManifold(c.mesh))
            note("compartment '" + c.name + "' is not a closed manifold");
        if (c.grossVolume <= 0)
            note("compartment '" + c.name + "' has non-positive volume");
        if (c.permeability <= 0 || c.permeability > 1)
            note("compartment '" + c.name + "' has permeability outside (0, 1]");
        subdivided += c.grossVolume;
    }
    // Compartments that overlap would push the total past the hull; the converse
    // (a shortfall) is normal and just means unmodelled voids and structure.
    if (subdivided > hullVolume * 1.001)
        note("compartments total " + std::to_string(subdivided) + " m3, more than the hull's " +
             std::to_string(hullVolume) + " m3 -- the subdivision overlaps");

    // The total is a weak instrument, and it missed a real one. A ship is
    // typically 85-90% subdivided, so one compartment can sit *wholly inside*
    // another and the sum still comes nowhere near the hull -- which is exactly
    // what the reference ferry did with a wing tank inside a hold, 217 m3 of
    // floodwater with nowhere physical to go. Pairwise is what catches it.
    //
    // Bounding boxes, not meshes: abutting compartments share a face and
    // intersect in zero volume, while a nested one intersects in its own. That
    // distinction is the whole check, and it is exact for the box-carved
    // compartments this engine builds. Advisory, because an L-shaped space could
    // overlap in box and not in fact.
    for (std::size_t i = 0; i < compartments.size(); ++i)
        for (std::size_t j = i + 1; j < compartments.size(); ++j) {
            const Compartment& a = compartments[i];
            const Compartment& b = compartments[j];
            Vec3 lo{std::max(a.bboxLo.x, b.bboxLo.x), std::max(a.bboxLo.y, b.bboxLo.y),
                    std::max(a.bboxLo.z, b.bboxLo.z)};
            Vec3 hi{std::min(a.bboxHi.x, b.bboxHi.x), std::min(a.bboxHi.y, b.bboxHi.y),
                    std::min(a.bboxHi.z, b.bboxHi.z)};
            const Vec3 e{hi.x - lo.x, hi.y - lo.y, hi.z - lo.z};
            if (e.x <= 0 || e.y <= 0 || e.z <= 0) continue;
            const double shared = e.x * e.y * e.z;
            const double smaller = std::min(a.grossVolume, b.grossVolume);
            if (smaller > 0 && shared > 0.01 * smaller)
                note("compartments '" + a.name + "' and '" + b.name + "' overlap by about " +
                     std::to_string(shared) + " m3 (" +
                     std::to_string(100.0 * shared / smaller) + "% of the smaller) -- that volume "
                     "would flood twice");
        }

    for (const Opening& o : openings) {
        const int n = static_cast<int>(compartments.size());
        if (o.a != kSea && (o.a < 0 || o.a >= n))
            note("opening '" + o.name + "' references a compartment that does not exist");
        if (o.b != kSea && (o.b < 0 || o.b >= n))
            note("opening '" + o.name + "' references a compartment that does not exist");
        if (o.a == o.b) note("opening '" + o.name + "' connects a space to itself");
        if (o.area < 0) note("opening '" + o.name + "' has negative area");
    }
    for (const Pump& p : pumps)
        if (p.compartment < 0 || p.compartment >= static_cast<int>(compartments.size()))
            note("pump '" + p.name + "' references a compartment that does not exist");

    if (lightshipMass <= 0) note("lightship mass is not positive");
    return problems;
}

int Ship::findCompartment(std::string_view name) const {
    for (std::size_t i = 0; i < compartments.size(); ++i)
        if (compartments[i].name == name) return static_cast<int>(i);
    return kSea;
}

double Ship::totalFloodwaterMass() const {
    double m = 0;
    for (const Compartment& c : compartments) m += c.waterVolume * seaDensity;
    return m;
}

// ---------------------------------------------------------------------------
// Free surfaces
// ---------------------------------------------------------------------------

void Ship::updateInternalFreeSurfaces(const Sea& sea) {
    (void)sea;
    const Mat3 R = state.orientation.toMat3();
    const Vec3 up = bodyFrameUp(R);

    for (Compartment& c : compartments) {
        c.waterVolume = std::clamp(c.waterVolume, 0.0, c.floodableVolume());

        if (c.waterVolume <= 1e-9) {
            // Dry. Most compartments on a ship are, most of the time, and the
            // solve below is the single most expensive thing in the tick -- so
            // this early-out is worth more than any micro-optimisation in it.
            // Park the surface below the deepest point of the space, so an
            // opening down there still reads as being in air.
            double deepest = kInf;
            for (const Vec3& v : c.mesh.verts) deepest = std::min(deepest, dot(up, v));
            c.surfaceOffset = deepest - 1.0;
            c.surfaceWorldZ = c.surfaceOffset + state.position.z;
            c.waterCentroid = (c.bboxLo + c.bboxHi) * 0.5;
        } else {
            // Permeability: the liquid occupies a geometric region larger than its
            // own volume, because part of that region is structure and cargo it
            // flows around.
            const double geometricRegion = c.waterVolume / std::max(c.permeability, 1e-6);
            const PlaneSweep sweep(c.mesh, up);
            c.surfaceOffset =
                sweep.solveOffsetForVolume(geometricRegion, c.grossVolume, c.surfaceOffset);
            c.surfaceWorldZ = c.surfaceOffset + state.position.z;

            const VolumeIntegral vi = sweep.below(c.surfaceOffset);
            c.waterCentroid = vi.volume > 1e-9 ? vi.centroid
                                               : (c.bboxLo + c.bboxHi) * 0.5;
        }

        // The datum the gas pressure is quoted at: the middle of the gas space,
        // between whatever the water surface is standing at and the deckhead.
        // For a box compartment that is the centroid of the gas region, which is
        // the height at which the well-mixed value m R T / V is the true local
        // pressure of a column in hydrostatic balance.
        // A dry compartment parks its free surface a metre *below* its own floor
        // so that an opening down there still reads as being in air, so the
        // surface offset cannot be used as the bottom of the gas space directly.
        const double top = boxTopAlong(up, c.bboxLo, c.bboxHi);
        const double floorOffset = boxBottomAlong(up, c.bboxLo, c.bboxHi);
        const double gasBottom = std::clamp(c.surfaceOffset, floorOffset, top);
        c.gasCentroidWorldZ = 0.5 * (gasBottom + top) + state.position.z;

        // The work the rising water does on the trapped air.
        //
        // At the default gasThermalTime of zero there is none to account for: the
        // structure takes it as fast as it is done, the gas stays at ambient, and
        // Boyle's law is the pressure law -- which is the model this file has
        // always had and the one the flooding scenarios are validated against. The
        // assignment below is exact, so every expression that reads the
        // temperature afterwards is character for character what it was.
        //
        // Given a relaxation time, the same compression is adiabatic on the way
        // in and the gas heats: at constant mass T V^(gamma-1) is conserved, which
        // is the isentropic relation T/T0 = (p/p0)^((gamma-1)/gamma) written the
        // way this loop can actually evaluate it. Mass changes are *not* handled
        // here -- solveFlowNetwork() has already taken the enthalpy they carried,
        // and taking a volume term and a mass term in one expression would be
        // taking the same joules twice.
        const double vaRaw = std::max(c.airVolume(), 0.0);
        if (c.gasThermalTime <= 0) {
            c.gasTemperature = kTAmbient;
        } else if (vaRaw > 1e-12 && c.lastAirVolume > 1e-12 && !c.ventedToAtmosphere) {
            c.gasTemperature *= std::pow(c.lastAirVolume / vaRaw, kGammaAir - 1.0);
        }
        c.lastAirVolume = vaRaw;

        if (c.ventedToAtmosphere) {
            // Open to the atmosphere, so the pressure is fixed and the mass is
            // whatever fills the space at this temperature. Air *arriving* from
            // outside arrives at ambient and cools whatever is in here; air
            // leaving takes its own temperature with it and changes nothing. Left
            // out, a vented space could be filled from the atmosphere without ever
            // losing the heat that filling dilutes -- energy from nowhere.
            const double mNew = kPatm * std::max(c.airVolume(), 0.0) / (kRAir * c.gasTemperature);
            if (mNew > c.airMass && mNew > 0)
                c.gasTemperature += (kTAmbient - c.gasTemperature) * ((mNew - c.airMass) / mNew);
            c.airPressure = kPatm;
            c.airMass = kPatm * std::max(c.airVolume(), 0.0) / (kRAir * c.gasTemperature);
        } else {
            const double va = std::max(c.airVolume(), 1e-3 * std::max(c.grossVolume, 1.0));
            c.airPressure = c.airMass * kRAir * c.gasTemperature / va;
        }
    }
}

// Heat out of the gas and into the steel around it.
//
// Exact exponential rather than an explicit Euler step, because the relaxation
// time a fire wants is seconds while the flooding solve steps at milliseconds,
// and an explicit step would be the one place in this file where a stiff term is
// integrated by the least stable method available. `dt / tau` of any size is
// stable and monotone here, and tau = infinity is the adiabatic limit rather than
// a division by zero.
void Ship::relaxGasToStructure(double dt) {
    for (Compartment& c : compartments) {
        // A zero thermal time is not relaxed here but pinned outright where the
        // pressure is computed, which is the one place that has to see an exact
        // kTAmbient. Repeating the pin here as well was dead code: it wrote a
        // value that updateInternalFreeSurfaces overwrote with the same value a
        // few lines later, and nothing between the two could read it.
        if (c.gasThermalTime <= 0 || !std::isfinite(c.gasThermalTime) || dt <= 0) continue;
        c.gasTemperature =
            kTAmbient + (c.gasTemperature - kTAmbient) * std::exp(-dt / c.gasThermalTime);
    }
}

Ship::SideState Ship::sideStateAt(int idx, const Vec3& worldPos, const Sea& sea) const {
    if (idx == kSea) {
        // The head at an opening is measured to the *local* free surface, not to
        // the still-water level: a breach under a crest floods faster than the
        // same breach under a trough, and that difference is a large part of why
        // damage in a seaway is worse than damage alongside.
        const double surface = sea.heightAt(worldPos.x, worldPos.y);
        if (worldPos.z < surface)
            return {kPatm + seaDensity * kGravity * (surface - worldPos.z), true, seaDensity,
                    seaDensity};
        return {kPatm, false, kRhoAirAmbient, kRhoAirAmbient};
    }
    const Compartment& c = compartments[static_cast<std::size_t>(idx)];
    if (c.waterVolume > 1e-9 && worldPos.z < c.surfaceWorldZ)
        return {c.airPressure + gasBuoyancyHead(c, c.surfaceWorldZ) +
                    seaDensity * kGravity * (c.surfaceWorldZ - worldPos.z),
                true, seaDensity, seaDensity};
    const double p = c.airPressure + gasBuoyancyHead(c, worldPos.z);
    return {p, false, p / (kRAir * c.gasTemperature), gasBuoyancyDensity(c.gasTemperature)};
}

// A hatch is in a deck, so its normal is the body's own +z and the sea is above
// it: a hatch to open water is a weather-deck hatch. Which of `a` and `b` is
// physically the upper one is therefore *not* the authored order -- it is decided
// in the world frame, every tick, because heeling her past ninety degrees puts the
// sea over what used to be a deckhead and the exchange has to reverse with her.
Ship::HorizontalSides Ship::horizontalSidesOf(const Opening& o, const Mat3& R) const {
    HorizontalSides h;
    if (o.kind != OpeningKind::Hatch || o.area <= 0) return h;

    // The cosine of the angle between the deck and the horizon. Heel her and the
    // hole stops facing up; at ninety degrees a deck hatch is a scuttle in a wall,
    // and since `Opening` carries no height there is no neutral plane for the
    // doorway integral to find, so zero is the honest limit for this struct rather
    // than a modelling shortcut. The *sign* is what says which way is down.
    const double tilt = (R * Vec3{0, 0, 1}).z;
    if (tilt == 0) return h;

    const auto heightOf = [&](int i) {
        if (i == kSea) return kInf;
        const Compartment& c = compartments[static_cast<std::size_t>(i)];
        return 0.5 * (c.bboxLo.z + c.bboxHi.z);
    };
    const double ha = heightOf(o.a);
    const double hb = heightOf(o.b);
    // Two spaces whose centres sit at the same height are not above and below each
    // other, so nothing here can say which way the heavy fluid would fall. That
    // also catches a hatch authored between the sea and the sea.
    if (!(ha != hb)) return h;

    const bool aAbove = (ha > hb) == (tilt > 0);
    h.valid = true;
    h.upper = aAbove ? o.a : o.b;
    h.lower = aAbove ? o.b : o.a;
    h.area = o.area * std::abs(tilt);
    return h;
}

// ---------------------------------------------------------------------------
// Orifice network
// ---------------------------------------------------------------------------

void Ship::solveFlowNetwork(double dt, const Sea& sea) {
    const Mat3 R = state.orientation.toMat3();
    const auto n = compartments.size();

    // Deltas are accumulated against the start-of-tick state so that opening order
    // does not bias the result -- an explicit Jacobi sweep over the network.
    //
    // `dMassT` is the third of these and it is the reason the *stored* state is a
    // temperature while the *accumulated* one is an energy. Temperatures do not
    // add: two openings delivering gas at 900 K and 300 K into the same space do
    // not deliver 1200 K, and a Jacobi sweep is a sum. m*T is the gas's internal
    // energy in units of cv, it is additive, and dividing it by the mass that
    // arrived with it is the mixing rule -- exact, for a gas with constant cv, and
    // with no root find anywhere.
    std::vector<double> dWater(n, 0.0), dAir(n, 0.0), dMassT(n, 0.0);

    auto waterAvailable = [&](int i) {
        return i == kSea ? kInf : compartments[i].waterVolume + dWater[i];
    };
    auto waterCapacity = [&](int i) {
        if (i == kSea) return kInf;
        return compartments[i].floodableVolume() - (compartments[i].waterVolume + dWater[i]);
    };
    auto airAvailable = [&](int i) {
        return i == kSea ? kInf : compartments[i].airMass + dAir[i];
    };

    auto gasTemperatureOf = [&](int i) {
        return i == kSea ? kTAmbient : compartments[i].gasTemperature;
    };

    for (Opening& o : openings) {
        o.lastFlow = 0;
        o.lastGasMassFlow = 0;
        o.lastGasDonorTemperature = 0;
        o.lastGasEnthalpyFlow = 0;
        o.lastExchangeDown = 0;
        o.lastExchangeUp = 0;
        o.lastExchangeMassDown = 0;
        o.lastExchangeMassUp = 0;
        o.lastExchangeUpper = kSea;
        if (!o.open || o.area <= 0) continue;

        const Vec3 worldPos = R * o.pos + state.position;
        const SideState sa = sideStateAt(o.a, worldPos, sea);
        const SideState sb = sideStateAt(o.b, worldPos, sea);

        const double dp = sa.pressure - sb.pressure;

        // --- Counter-current exchange through a horizontal opening ------------
        //
        // Tested *before* the noise floor below, because the case this exists for
        // is precisely the one the noise floor throws away: dp at rest at the
        // orifice and a tonne a second going through it anyway.
        //
        // A doorway has a height, so `fire.cpp` integrates Torricelli up it and
        // splits the band at the neutral plane. A hatch has none -- there is one
        // elevation and one dp -- so what is split instead is the *area*: the
        // heavy fluid falls through part of the hole while the light one rises
        // through the rest. The buoyancy head drives the two in opposite
        // directions, so it *adds* to each one's own driving pressure:
        //
        //     dp_b = (rho_up - rho_lo) g D/2,   D = sqrt(4A/pi)
        //     v_dn = sqrt(2 (dp_net + dp_b) / rho_up)    the heavy fluid, falling
        //     v_up = sqrt(2 (dp_b - dp_net) / rho_lo)    the light one, rising
        //
        // **What sets the split between them is the net.** Two linear equations,
        // and no free parameter in them:
        //
        //     A_dn + A_up             = Cd A      they share one hole
        //     A_dn v_dn - A_up v_up   = q_net     and pass the net between them
        //
        // where `q_net` is the single-orifice Torricelli on the net head over the
        // whole area -- exactly what this file gives every other opening.
        //
        // **That second equation is the one that took a second try, and the first
        // attempt failed in a way worth recording.** Requiring the two *volumes* to
        // be equal instead -- neither space changes volume, so what falls has to be
        // let past by what rises -- is right at rest and wrong everywhere else. It
        // drives the exchange to zero as the net head approaches dp_b, so the model
        // reported a hole passing **nothing at all** in a window 0.7 mm of water
        // deep, and then jumped to 2.3 m^3/s the instant the net took over. Solved
        // against the net instead, the two regimes agree at the edge *identically*:
        // at dp_net = dp_b, v_up is zero, A_dn comes out at Cd A / sqrt(2), and
        // A_dn v_dn is q_net to the last bit.
        //
        // **D/2 is the one modelling choice in it, and it is a choice of length
        // scale rather than a fitted coefficient.** The head that drives the
        // overturning is the weight of a plug of fluid about half a hole across --
        // an opening's own size is the only length in the problem, which is why
        // every published correlation for this (Epstein, Cooper) reports a rate
        // going as `sqrt(g') D^{5/2}`, and the form above reproduces that exponent
        // exactly. In the ship's own limit, water over air with no net head, it
        // collapses to `Q = Cd A sqrt(g D)` each way. The discharge coefficient is
        // the opening's own, so an author who wants a different constant has one to
        // turn.
        //
        // Two limits, and both are exact rather than approached:
        //
        //   * `dp_b = 0` -- equal densities, or light already lying on heavy. The
        //     branch below is simply not taken, and the single-orifice code beneath
        //     runs character for character as it always did, which is what keeps a
        //     cold ship bit-identical. Two cold compartments give a density
        //     difference of exactly +0.0 by IEEE arithmetic and not merely a small
        //     one.
        //   * `dp_net = 0` -- the pure exchange, and the whole point. `q_net`
        //     vanishes, the areas fall out as `A_dn : A_up = v_up : v_dn`, and the
        //     two volumes are then equal: `Q = Cd A sqrt(2 dp_b) / (sqrt(rho_up) +
        //     sqrt(rho_lo))` each way, with no net whatever. That is precisely the
        //     state a model reading one dp at the orifice centre calls a hatch at
        //     rest.
        //
        // Past `|dp_net| >= dp_b` the hole is flushed one way and the exchange has
        // stopped of its own accord -- Cooper's critical pressure difference,
        // falling out of the two driving pressures rather than being imposed.
        const HorizontalSides hs = horizontalSidesOf(o, R);
        if (hs.valid) {
            const SideState& su = hs.upper == o.a ? sa : sb;
            const SideState& sl = hs.upper == o.a ? sb : sa;
            const double drho = su.buoyancyDensity - sl.buoyancyDensity;
            const double dpB = drho > 0 ? drho * kGravity * 0.5 * hydraulicDiameter(o.area) : 0.0;
            const double dpDown = su.pressure - sl.pressure;
            if (dpB > std::abs(dpDown) && su.flowDensity > 0 && sl.flowDensity > 0) {
                const double vDown = std::sqrt(2.0 * (dpDown + dpB) / su.flowDensity);
                const double vUp = std::sqrt(2.0 * (dpB - dpDown) / sl.flowDensity);
                const double cdA = o.dischargeCoeff * hs.area;
                // The net the pair passes, at the density of whichever side is
                // donating it -- the same expression, on the same head, that the
                // single-orifice solve below would apply to this opening.
                const double rhoNet = dpDown >= 0 ? su.flowDensity : sl.flowDensity;
                const double qNet = (dpDown >= 0 ? 1.0 : -1.0) * cdA *
                                    std::sqrt(2.0 * std::abs(dpDown) / rhoNet);
                // Both areas come out non-negative for every net head inside the
                // window, so there is nothing to clamp here: `cdA v_up` always
                // dominates a negative `q_net` and `cdA v_dn` a positive one,
                // because each stream carries dp_b on top of the net that is trying
                // to cancel it.
                const double aDown = (qNet + cdA * vUp) / (vDown + vUp);
                const double aUp = cdA - aDown;
                double dvDown = aDown * vDown * dt;
                double dvUp = aUp * vUp * dt;

                // **Clamped by one common factor, not stream by stream.** Scaling
                // them apart would change the net the pair passes, which is the
                // quantity the split was solved against; scaling them together
                // leaves the net proportional and leaves the equal volumes at rest
                // exactly equal.
                const auto share = [&](int from, int to, const SideState& s, double dv) {
                    if (dv <= 0) return 1.0;
                    const double room =
                        s.isWater ? std::min(std::max(waterAvailable(from), 0.0),
                                             std::max(waterCapacity(to), 0.0))
                                  : std::max(airAvailable(from), 0.0) / s.flowDensity;
                    return room < dv ? room / dv : 1.0;
                };
                const double scale = std::min(share(hs.upper, hs.lower, su, dvDown),
                                              share(hs.lower, hs.upper, sl, dvUp));
                dvDown *= scale;
                dvUp *= scale;

                // One stream. Returns the mass it carried, so the diagnostics
                // below do not have to re-derive which density applied.
                const auto move = [&](int from, int to, const SideState& s, double dv) {
                    const double dm = s.flowDensity * dv;
                    if (s.isWater) {
                        if (from != kSea) dWater[from] -= dv;
                        if (to   != kSea) dWater[to]   += dv;
                        return dm;
                    }
                    const double dE = kGammaAir * gasTemperatureOf(from) * dm;
                    if (from != kSea) { dAir[from] -= dm; dMassT[from] -= dE; }
                    if (to   != kSea) { dAir[to]   += dm; dMassT[to]   += dE; }
                    return dm;
                };
                const double mDown = move(hs.upper, hs.lower, su, dvDown);
                const double mUp = move(hs.lower, hs.upper, sl, dvUp);
                const double perSecond = dt > 0 ? 1.0 / dt : 0.0;
                o.lastExchangeDown = dvDown * perSecond;
                o.lastExchangeUp = dvUp * perSecond;
                o.lastExchangeMassDown = mDown * perSecond;
                o.lastExchangeMassUp = mUp * perSecond;
                o.lastExchangeUpper = hs.upper;

                // The same two streams as a signed a -> b net, for consumers that
                // only ever wanted one number. Water crosses one way at most -- the
                // rising stream is never the denser phase -- so `lastFlow` still
                // says everything there is to say about it.
                const double sign = hs.upper == o.a ? 1.0 : -1.0;
                o.lastFlowWasWater = su.isWater;
                if (su.isWater) o.lastFlow = sign * o.lastExchangeDown;
                const double tUp = gasTemperatureOf(hs.upper);
                const double tLo = gasTemperatureOf(hs.lower);
                if (!su.isWater) {
                    o.lastGasMassFlow += sign * o.lastExchangeMassDown;
                    o.lastGasEnthalpyFlow += sign * kGammaAir * tUp * o.lastExchangeMassDown;
                }
                if (!sl.isWater) {
                    o.lastGasMassFlow -= sign * o.lastExchangeMassUp;
                    o.lastGasEnthalpyFlow -= sign * kGammaAir * tLo * o.lastExchangeMassUp;
                }
                // Left at zero when both streams were gas: two donors at two
                // temperatures, and a single donor temperature would be a
                // diagnostic that lies. `lastGasEnthalpyFlow` is the one that still
                // adds up.
                o.lastGasDonorTemperature = su.isWater && !sl.isWater ? tLo : 0.0;
                continue;
            }
        }

        if (std::abs(dp) < 1e-3) continue;  // below the noise floor of the solve

        const bool aIsDonor = dp > 0;
        const int donor = aIsDonor ? o.a : o.b;
        const int recv  = aIsDonor ? o.b : o.a;
        // What actually moves is whatever is sitting against the orifice on the
        // high-pressure side. This one line is why a hole above the internal
        // waterline vents air, and the same hole below it admits water.
        const bool water = aIsDonor ? sa.isWater : sb.isWater;

        const double pDonor = std::max(sa.pressure, sb.pressure);
        // Hot gas is thin gas, and a thin gas leaves faster through the same hole
        // for the same head: Torricelli goes as 1/sqrt(rho). This is the density
        // the *donor* is at, which is now a question with an answer.
        const double tDonor = gasTemperatureOf(donor);
        const double rho = water ? seaDensity : pDonor / (kRAir * tDonor);

        // Torricelli / sharp-edged orifice.
        const double q = o.dischargeCoeff * o.area * std::sqrt(2.0 * std::abs(dp) / rho);

        if (water) {
            double dv = q * dt;
            dv = std::min(dv, std::max(waterAvailable(donor), 0.0));
            dv = std::min(dv, std::max(waterCapacity(recv), 0.0));
            if (donor != kSea) dWater[donor] -= dv;
            if (recv  != kSea) dWater[recv]  += dv;
            o.lastFlow = (aIsDonor ? 1.0 : -1.0) * (dt > 0 ? dv / dt : 0.0);
        } else {
            double dm = rho * q * dt;
            dm = std::min(dm, std::max(airAvailable(donor), 0.0));

            // Air is by far the stiffest part of this network: a nearly-full
            // compartment has a tiny gas volume, so an explicit step that ignores
            // the target state will blow straight past equilibrium and ring. Clamp
            // the transfer to the mass that equalises the two pressures, which is
            // the implicit answer this step is trying to approximate anyway.
            //
            // **The clamp survives the temperature intact, and it stays a closed
            // form.** The worry is real -- the mass that equalises the pressures
            // now depends on the temperature the receiving side ends up at, and
            // that temperature depends on how much mass arrives -- but the
            // circularity resolves algebraically rather than needing an implicit
            // solve. Pressure is (gamma-1)U/V, so it is *linear in the energy
            // delivered*, and the energy delivered is linear in dm: the (m + dm)
            // in the mixing rule cancels against the m + dm in the gas law and
            // never has to be inverted. What each branch below solves is one
            // linear equation.
            //
            // The two branches are not two approximations of one thing, they are
            // the clamp for the two energy equations a compartment can be running.
            // A compartment at gasThermalTime = 0 is isothermal by construction --
            // arriving gas does not change its temperature, because the structure
            // takes the heat -- so the mass that equalises it is Boyle's, exactly
            // as before. One with a relaxation time keeps what arrives, including
            // the flow work, and equalises at a different mass.
            // **Both clamps are stated at the compartment's own pressure datum,
            // not at the opening.** `sa` and `sb` are pressures *at the hole*,
            // while the mass a compartment holds sets its pressure at the middle
            // of its gas space, and the two differ by that side's own buoyancy
            // head. Equalising the wrong pair silently pins the transfer to zero:
            // a hot space whose gas is 25 Pa lighter at a high doorway reads, at
            // its own datum, as already equal to the cold space next door, and
            // sends nothing at all through a door the physics says should be
            // pouring smoke. Subtracting each side's own head puts the target back
            // on the datum the mass is measured against. Zero at ambient, so a
            // cold ship is untouched.
            const double pDonorSide = aIsDonor ? sa.pressure : sb.pressure;
            const double pRecvSide  = aIsDonor ? sb.pressure : sa.pressure;
            if (recv != kSea) {
                const Compartment& rc = compartments[recv];
                const double va = std::max(rc.airVolume(), 1e-3 * std::max(rc.grossVolume, 1.0));
                const double pTarget = pDonorSide - gasBuoyancyHead(rc, worldPos.z);
                if (rc.gasThermalTime <= 0) {
                    const double mEq = pTarget * va / (kRAir * rc.gasTemperature);
                    dm = std::min(dm, std::max(0.0, mEq - (rc.airMass + dAir[recv])));
                } else {
                    // (m T)_new = (m T) + gamma T_donor dm, and p = R (m T) / V.
                    const double mtEq = pTarget * va / kRAir;
                    const double mt = rc.airMass * rc.gasTemperature + dMassT[recv];
                    dm = std::min(dm, std::max(0.0, (mtEq - mt) / (kGammaAir * tDonor)));
                }
            }
            if (donor != kSea) {
                const Compartment& dc = compartments[donor];
                const double va = std::max(dc.airVolume(), 1e-3 * std::max(dc.grossVolume, 1.0));
                const double pTarget = pRecvSide - gasBuoyancyHead(dc, worldPos.z);
                if (dc.gasThermalTime <= 0) {
                    const double mEq = pTarget * va / (kRAir * dc.gasTemperature);
                    dm = std::min(dm, std::max(0.0, (dc.airMass + dAir[donor]) - mEq));
                } else {
                    // The same equation from the other end. It also keeps the
                    // donor's energy positive: the bound is (m T - p V/R)/(gamma T),
                    // which is at most m/gamma, so the blowdown term below can
                    // never take more energy out than the gas had.
                    const double mtEq = pTarget * va / kRAir;
                    const double mt = dc.airMass * dc.gasTemperature + dMassT[donor];
                    dm = std::min(dm, std::max(0.0, (mt - mtEq) / (kGammaAir * tDonor)));
                }
            }

            // What crosses the hole is mass and the energy riding on it, and the
            // energy is *enthalpy*, not internal energy: the donor has to do the
            // work of pushing the gas out through the orifice, and the receiver
            // collects that work as well as the heat. Taking cv T here instead of
            // cp T would leave the donor too warm and the receiver too cool by
            // exactly the flow work, and would lose the two textbook limits this
            // network is made of -- a vessel blowing down cools isentropically,
            // and one being charged from a reservoir ends up hotter than the
            // reservoir. It is antisymmetric across the opening, so the sum over
            // the network is conserved to roundoff whatever the clamps do.
            const double dEnergy = kGammaAir * tDonor * dm;
            if (donor != kSea) { dAir[donor] -= dm; dMassT[donor] -= dEnergy; }
            if (recv  != kSea) { dAir[recv]  += dm; dMassT[recv]  += dEnergy; }
            o.lastFlow = (aIsDonor ? 1.0 : -1.0) * (dt > 0 ? (rho > 0 ? dm / rho / dt : 0.0) : 0.0);
            o.lastGasMassFlow = (aIsDonor ? 1.0 : -1.0) * (dt > 0 ? dm / dt : 0.0);
            o.lastGasDonorTemperature = tDonor;
            o.lastGasEnthalpyFlow = kGammaAir * tDonor * o.lastGasMassFlow;
        }
        o.lastFlowWasWater = water;
    }

    for (Pump& p : pumps) {
        p.lastFlow = 0;
        if (!p.on) continue;
        Compartment& c = compartments[static_cast<std::size_t>(p.compartment)];
        // Centrifugal pumps lose output as discharge head rises; overboard discharge
        // means lifting from the compartment's water surface to the sea surface.
        const double head = std::max(0.0, sea.level - c.surfaceWorldZ);
        const double efficiency = std::clamp(1.0 - head / std::max(p.maxHead, 1e-6), 0.0, 1.0);
        const double q = p.capacity * efficiency;
        const double dv = std::min(q * dt, std::max(c.waterVolume + dWater[p.compartment], 0.0));
        dWater[p.compartment] -= dv;
        p.lastFlow = dt > 0 ? dv / dt : 0.0;
    }

    for (std::size_t i = 0; i < n; ++i) {
        Compartment& c = compartments[i];
        c.waterVolume =
            std::clamp(c.waterVolume + dWater[i], 0.0, c.floodableVolume());
        const double mNew = std::max(0.0, c.airMass + dAir[i]);
        // Back out of energy into the temperature that is actually stored. The
        // guard is not a tolerance: a compartment can be flooded solid, and a gas
        // with no mass has no temperature to speak of. Keeping the last one it had
        // is harmless -- nothing reads it while the mass is zero -- and inventing
        // a zero would put an absolute temperature of 0 K into the pressure law.
        const double mt = c.airMass * c.gasTemperature + dMassT[i];
        if (mNew > 1e-12 && mt > 0) c.gasTemperature = mt / mNew;
        c.airMass = mNew;
    }
}

// ---------------------------------------------------------------------------
// Mass properties
// ---------------------------------------------------------------------------

Ship::MassProperties Ship::massProperties() const {
    MassProperties mp;
    mp.mass = lightshipMass;
    Vec3 moment = lightshipCog * lightshipMass;

    for (const Compartment& c : compartments) {
        const double m = c.waterVolume * seaDensity;
        if (m <= 0) continue;
        mp.mass += m;
        moment += c.waterCentroid * m;
    }
    mp.cog = mp.mass > 0 ? moment / mp.mass : lightshipCog;

    Mat3 I = Mat3::zero();
    I(0, 0) = lightshipMass * gyradii.x * gyradii.x;
    I(1, 1) = lightshipMass * gyradii.y * gyradii.y;
    I(2, 2) = lightshipMass * gyradii.z * gyradii.z;
    I = I + pointInertia(lightshipMass, lightshipCog - mp.cog);

    for (const Compartment& c : compartments) {
        const double m = c.waterVolume * seaDensity;
        if (m <= 0) continue;
        // The water body's own inertia, approximated as a box filling the
        // compartment footprint to the current level. Omitting it would make a
        // half-flooded hold roll like a point mass and badly under-damp the ship.
        const Vec3 ext = c.bboxHi - c.bboxLo;
        const Vec3 waterExt{ext.x, ext.y, std::max(ext.z * c.fillFraction(), 1e-3)};
        I = I + boxInertia(m, waterExt) + pointInertia(m, c.waterCentroid - mp.cog);
    }
    mp.inertiaAboutCog = I;
    return mp;
}

// ---------------------------------------------------------------------------
// Rigid body
// ---------------------------------------------------------------------------

double Ship::equilibriumDraftAt(const Quat& orientation, const Sea& sea,
                                double targetMass) const {
    const Vec3 up = bodyFrameUp(orientation.toMat3());
    const double offset =
        solvePlaneOffsetForVolume(hull, up, targetMass / seaDensity, -kInf, kInf);
    // world_z(x) = dot(up, x) + position.z, so the surface sits at offset + position.z.
    return sea.level - offset;
}

void Ship::integrateRigidBody(double dt, const Sea& sea) {
    const Mat3 R = state.orientation.toMat3();
    const Vec3 up = bodyFrameUp(R);
    const double planeOffset = sea.level - state.position.z;

    const PlaneSweep hullSweep(hull, up);
    const MassProperties mp = massProperties();
    const Vec3 cogWorld = R * mp.cog + state.position;

    // Buoyancy. Flat water is clipped by a single plane in the body frame, which
    // is exact and cheap. A wave field is not a plane and is a function of world
    // x and y, so the hull has to be carried into world coordinates and
    // integrated against the surface itself -- this is the nonlinear
    // Froude-Krylov restoring force, evaluated over the instantaneous wetted
    // surface rather than about a mean waterline.
    VolumeIntegral sub;
    Vec3 cbWorld;
    if (sea.flat()) {
        sub = hullSweep.below(planeOffset);
        cbWorld = R * sub.centroid + state.position;
    } else {
        worldHullScratch_ = hull;
        for (Vec3& v : worldHullScratch_.verts) v = R * v + state.position;
        sub = integrateBelowSurface(worldHullScratch_,
                                    [&](double x, double y) { return sea.heightAt(x, y); });
        cbWorld = sub.centroid;  // already world, no transform needed
    }

    Vec3 force{0, 0, -mp.mass * kGravity};
    Vec3 torque{0, 0, 0};  // about cogWorld

    Vec3 buoyancy{0, 0, 0};
    if (sub.volume > 0) {
        buoyancy = Vec3{0, 0, seaDensity * kGravity * sub.volume};
        force += buoyancy;
        torque += cross(cbWorld - cogWorld, buoyancy);
    }

    // --- Damping ------------------------------------------------------------
    // Stiffnesses are measured, not assumed: the waterplane area and the angular
    // restoring rates come from finite-differencing the same hydrostatic integral
    // that produced the buoyancy force. Refreshed periodically because they change
    // far more slowly than the ship moves.
    if (stabilityRefreshCounter_ % 16 == 0) {
        const double h = 0.05;
        cachedWaterplaneArea_ = (hullSweep.below(planeOffset + h).volume -
                                 hullSweep.below(planeOffset - h).volume) / (2 * h);

        // **The perturbation has to rotate the ship about its centre of gravity.**
        // The damping this feeds is applied to a torque about the cog, so the
        // stiffness paired with it must be about the cog too. Rotating about the
        // body origin instead -- which is what leaving `state.position` alone
        // does -- moves the centre of gravity out from under the point the moment
        // is taken about, and the weight then contributes a moment of rho g V
        // times KG that no righting arm has in it.
        //
        // That term is nothing next to the longitudinal metacentric height, so
        // pitch never noticed; against a transverse GM of 2.3 m and a KG of 5 m
        // it was larger than the quantity being measured, and a nominal
        // zetaRoll of 0.08 was delivering 0.144 of critical -- found by timing a
        // free decay's logarithmic decrement, not by any check on the stiffness.
        auto angularStiffness = [&](const Vec3& axisWorld) {
            const double eps = 2e-3;
            const Quat q2 = Quat::fromAxisAngle(axisWorld, eps) * state.orientation;
            const Mat3 R2 = q2.toMat3();
            const Vec3 position2 = cogWorld - R2 * mp.cog;
            const double planeOffset2 = sea.level - position2.z;
            const VolumeIntegral s2 = PlaneSweep(hull, bodyFrameUp(R2)).below(planeOffset2);
            const Vec3 fb2{0, 0, seaDensity * kGravity * s2.volume};
            const Vec3 cb2 = R2 * s2.centroid + position2;
            const Vec3 t2 = cross(cb2 - cogWorld, fb2);
            return std::max(0.0, -(dot(t2, axisWorld) - dot(torque, axisWorld)) / eps);
        };
        const Vec3 xw = R * Vec3{1, 0, 0};
        const Vec3 yw = R * Vec3{0, 1, 0};
        cachedKRoll_  = angularStiffness(xw);
        cachedKPitch_ = angularStiffness(yw);
    }
    ++stabilityRefreshCounter_;

    const double dispMass = seaDensity * std::max(sub.volume, 0.0);
    Vec3 mEff{mp.mass + addedMassSurge * dispMass,
              mp.mass + addedMassSway  * dispMass,
              mp.mass + addedMassHeave * dispMass};

    Mat3 Ieff = mp.inertiaAboutCog;
    Ieff(0, 0) *= (1.0 + addedInertiaRoll);
    Ieff(1, 1) *= (1.0 + addedInertiaPitch);
    Ieff(2, 2) *= (1.0 + addedInertiaYaw);

    // With a radiation model, the guessed coefficients give way to this hull's
    // own infinite-frequency added mass.
    //
    // **The transfer is not optional.** stripTheoryTable() assembles A_inf about
    // the body-frame origin -- midship on the baseline -- while everything below
    // takes moments about the centre of gravity, and rotational added mass is no
    // more origin-independent than a moment of inertia is. Referring it to the
    // cog costs one congruence, A' = T^T A T, and it is not a small correction:
    // on the reference barge A_44 falls by 39%, because a section's added-mass
    // roll centre sits near the waterline and the baseline is a long lever away
    // from it. Heave comes through untouched (translational added mass is
    // reference-point independent) and pitch moves only through the longitudinal
    // offset of the cog, since strip theory's surge added mass is zero.
    //
    // Only the diagonal of the result is taken: the integrator below inverts a
    // 3x3 inertia and divides by a per-axis mass, so there is nowhere for the
    // sway-roll coupling A_24 to go. That is an approximation this integrator's
    // shape forces rather than one the physics justifies -- but referring the
    // matrix to the cog first shrinks the term being dropped by roughly eight
    // times on the barge, because the cog sits near the point where sway and roll
    // decouple. See docs/02-simulation.md.
    if (radiation.has_value()) {
        const Matrix6 aInf = transferAddedMass(radiation->addedMassInfinite(), mp.cog);
        // Surge keeps its coefficient. Strip theory does not decline to answer
        // for surge by oversight -- a strip has no longitudinal radiation problem
        // at all, and surge added mass is entirely a three-dimensional end
        // effect -- so A_inf[0][0] is a structural zero, not a measurement of
        // zero. Taking it would delete a real term rather than improve it.
        mEff.y = mp.mass + aInf[1][1];
        mEff.z = mp.mass + aInf[2][2];
        Ieff(0, 0) = mp.inertiaAboutCog(0, 0) + aInf[3][3];
        Ieff(1, 1) = mp.inertiaAboutCog(1, 1) + aInf[4][4];
        Ieff(2, 2) = mp.inertiaAboutCog(2, 2) + aInf[5][5];
    }

    // Surge added mass, the one entry strip theory structurally cannot supply,
    // comes from the manoeuvring set when there is one: m_x' non-dimensionalised
    // on 0.5 rho L^2 d. This is exactly the hole documented above being filled by
    // the only model here that has an answer for it.
    if (propulsion.has_value()) {
        const HullParams& mmg = propulsion->hull;
        const double reference = 0.5 * seaDensity * mmg.length * mmg.length * mmg.draft;
        mEff.x = mp.mass + mmg.addedMassSurge * reference;
    }

    const Vec3 vBody = R.transposed() * state.velocity;
    const Vec3 wBody = R.transposed() * state.angularVelocity;

    const double kHeave = seaDensity * kGravity * std::max(cachedWaterplaneArea_, 1.0);
    Vec3 cLin{2 * zetaHeave * std::sqrt(kHeave * mEff.x) * 0.05,
              2 * zetaHeave * std::sqrt(kHeave * mEff.y) * 0.30,
              2 * zetaHeave * std::sqrt(kHeave * mEff.z)};
    Vec3 cAng{2 * zetaRoll  * std::sqrt(cachedKRoll_  * Ieff(0, 0)),
              2 * zetaPitch * std::sqrt(cachedKPitch_ * Ieff(1, 1)),
              0.02 * Ieff(2, 2)};

    // With radiation attached, the modal dampers in the modes it covers have to
    // go, or the ship is damped twice.
    //
    // That is not a subtlety, it is what these coefficients *were*: a lumped
    // stand-in for radiation, which is the dominant damping in heave, sway and
    // pitch. Leaving them in alongside the real thing cost 27% of mid-frequency
    // heave in the first version of this coupling -- a wrong answer that looks
    // entirely reasonable, since a ship that moves too little in waves reads as
    // stiff rather than as broken.
    //
    // Yaw keeps its damper, and so does surge -- strip theory contributes no
    // surge radiation at all. Roll keeps this one only until Ikeda is attached
    // below: roll radiation damping is genuinely small and the mechanism that
    // matters is viscous, so a ship with radiation but no viscous roll model
    // would otherwise have the one mode that most needs damping running on almost
    // none. Quadratic drag is untouched throughout: it is viscous and separate.
    if (radiation.has_value()) {
        cLin.y = 0.0;
        cLin.z = 0.0;
        cAng.y = 0.0;
    }

    // Viscous roll damping, Ikeda's method, in place of roll's fraction-of-
    // critical stand-in -- the fourth time in this file that coupling a real
    // model in has meant deleting a lumped one that was secretly doing its job.
    // zetaRoll was fitted to sit where Ikeda puts a ro-pax with bilge keels (8%
    // of critical against a measured 6.3%), so leaving both in would not look
    // wrong; it would just double the damping of the mode that decides whether a
    // damaged ship capsizes.
    //
    // **The operating point.** B44 is an equivalent *linear* coefficient and is
    // only valid at the amplitude, frequency and speed it was asked for, all
    // three of which move every tick. Evaluated per tick rather than cached,
    // because one call is 52 ns against a 16 us tick -- 0.3%, so a cache would
    // buy nothing but a stale answer. The three inputs come from quantities this
    // integrator already has in hand:
    //
    //   * frequency: the natural roll frequency sqrt(C44 / (I44 + A44)), from the
    //     measured hydrostatic stiffness and the effective inertia two lines up.
    //     Roll is lightly damped and sharply resonant, so a ship's roll sits at
    //     its natural frequency whatever is driving it; and Ikeda's frequency
    //     dependence is only sqrt(omega) for friction and linear for the rest, so
    //     a frequency wrong by a fifth is a B44 wrong by less than the method's
    //     own +/-25%. Estimating it from recent motion instead would need filter
    //     state, would lag every change of loading, and would have nothing to say
    //     at all about a ship that is not currently rolling.
    //   * amplitude: the phase-plane radius sqrt(phi^2 + (phidot/omega_n)^2),
    //     which for a lightly damped oscillator *is* the envelope of the roll it
    //     is presently executing -- exactly the amplitude equivalent
    //     linearisation is defined at, available instantaneously and with no
    //     history to keep. Reading the heel alone would report zero damping every
    //     time the ship passed through upright, which is where it is rolling
    //     fastest.
    //   * speed: surge along the ship's own bow. state.velocity is a *world*
    //     vector and its x component is the ship's speed only while the ship
    //     still points along world x.
    if (rollDampingForm.has_value()) {
        // The roll axis is loading, not form: it moves as the ship floods, so it
        // is taken from the live centre of gravity rather than snapshotted.
        rollDampingForm->rollAxisAboveKeel = mp.cog.z - hullLo.z;
        rollDampingForm->seaDensity = seaDensity;

        const double omegaRoll = std::sqrt(cachedKRoll_ / std::max(Ieff(0, 0), 1e-9));
        double heel = 0, trim = 0;
        heelTrimFromRotation(R, heel, trim);

        rollCondition.rollFrequency = omegaRoll;
        rollCondition.rollAmplitude =
            omegaRoll > 0 ? std::hypot(heel, wBody.x / omegaRoll) : std::abs(heel);
        rollCondition.forwardSpeed = vBody.x;
        rollDampingApplied = rollDamping(*rollDampingForm, rollCondition);

        // Zero means the form is not usable -- rollDamping() reports every
        // component as zero for a hull it cannot make sense of -- and the stand-in
        // stays. It is not a guard against a capsizing ship: when the
        // finite-difference roll stiffness clamps to zero the stand-in is
        // 2 zeta sqrt(0 * I), which is also zero, so there is nothing to fall back
        // *to*. What the ship keeps in that state is its quadratic drag.
        if (rollDampingApplied.total > 0) cAng.x = rollDampingApplied.total;
    }

    // Quadratic drag on the projected areas of the box around the hull.
    const Vec3 ext = hullHi - hullLo;
    const double draft = std::max(0.1, sea.level - (state.position.z + hullLo.z));
    const Vec3 projArea{ext.y * draft, ext.x * draft, ext.x * ext.y};
    Vec3 dragCoeff{0.10, 1.00, 1.50};

    // The same argument as radiation, one deck down. A drag coefficient of 0.10
    // on the bow's projected area is a stand-in for hull resistance, and 1.00 on
    // the side is a stand-in for cross-flow drag. The MMG polynomial computes
    // both properly -- R_0' for resistance, Y_v and its higher terms for the
    // lateral force -- so with a manoeuvring set attached the stand-ins go, and
    // the yaw damper with them, since N_r is what that was standing in for.
    // Heave keeps its drag: nothing in a horizontal-plane model speaks to it.
    if (propulsion.has_value()) {
        dragCoeff.x = 0.0;
        dragCoeff.y = 0.0;
        cLin.x = 0.0;
        cLin.y = 0.0;
        cAng.z = 0.0;
    }

    Vec3 fBody = R.transposed() * force;
    Vec3 tBody = R.transposed() * torque;
    for (int i = 0; i < 3; ++i) {
        const double v = vBody[i];
        fBody[i] -= cLin[i] * v
                  + 0.5 * seaDensity * dragCoeff[i] * projArea[i] * std::abs(v) * v;
        tBody[i] -= cAng[i] * wBody[i];
    }

    // Propeller, rudder and the MMG hull polynomial, in the horizontal plane.
    //
    // Ship owns the motion: the manoeuvring state is overwritten from the rigid
    // body here and only read back as forces, so the two cannot disagree about
    // where the ship is or how fast it is going.
    if (propulsion.has_value()) {
        propulsion->density = seaDensity;
        propulsion->state.surgeSpeed = vBody.x;
        propulsion->state.swaySpeed = vBody.y;
        propulsion->state.yawRate = wBody.z;

        const PlanarForces planar = propulsion->forces();
        const Vec3 horizontal{planar.surge, planar.sway, 0.0};
        fBody += horizontal;
        // X and Y act at midship -- the body origin -- while this integrator
        // takes moments about the centre of gravity, so the offset contributes a
        // yaw moment of its own. On a hull with x_G forward of midship that term
        // is not small, and dropping it would make the ship turn about the wrong
        // point in a way that still looks like turning.
        tBody.z += planar.yawMoment + cross(Vec3{} - mp.cog, horizontal).z;
    }

    // The wave-memory force: waves this ship radiated earlier, still nearby and
    // still pushing back. Advanced with the body-frame velocity and subtracted,
    // per the sign convention in radiation.hpp.
    //
    // The velocity handed over is the body origin's, which is the point the
    // table is referenced to, so no transfer is owed on the way in. The moment
    // that comes back is about that same origin, and tBody is about the centre of
    // gravity -- so the memory *force* contributes a moment of its own through
    // the offset, exactly as the propulsion forces do above. Dropping it would
    // put the whole radiation memory on the wrong lever in roll.
    if (radiation.has_value()) {
        radiation->step({vBody.x, vBody.y, vBody.z, wBody.x, wBody.y, wBody.z}, dt);
        const std::array<double, 6> memory = radiation->memoryForce();
        const Vec3 memoryForce{memory[0], memory[1], memory[2]};
        const Vec3 memoryMoment = Vec3{memory[3], memory[4], memory[5]} -
                                  cross(mp.cog, memoryForce);
        for (int i = 0; i < 3; ++i) {
            fBody[i] -= memoryForce[i];
            tBody[i] -= memoryMoment[i];
        }
    }

    // Contact with another hull, and anything else outside this ship's own
    // physics. Added after the damping loop because a contact force is not
    // something the sea damps, and cleared once consumed.
    //
    // It is divided by the *effective* mass below, added mass included, which is
    // right -- a struck ship accelerates with the water it has to shove aside --
    // and which means the two ships in a collision do not conserve momentum
    // between themselves alone. The difference is in the water. See
    // collision.hpp; the conservation laws are asserted against the bare rigid
    // bodies, where they hold exactly.
    fBody += R.transposed() * externalForce;
    tBody += R.transposed() * externalMoment;
    externalForce = {};
    externalMoment = {};

    const Vec3 aBody{fBody.x / mEff.x, fBody.y / mEff.y, fBody.z / mEff.z};
    const Vec3 alphaBody = inverse(Ieff) * (tBody - cross(wBody, Ieff * wBody));

    // Integrate about the centre of gravity, then map back to the body origin.
    const Vec3 rc = R * mp.cog;
    Vec3 vCog = state.velocity + cross(state.angularVelocity, rc);
    vCog += (R * aBody) * dt;
    state.angularVelocity += (R * alphaBody) * dt;
    state.velocity = vCog - cross(state.angularVelocity, rc);
    state.position += state.velocity * dt;
    state.orientation = state.orientation.integrated(state.angularVelocity, dt);
}

RadiationTable Ship::attachRadiation(double waterlineZ, int stationCount, int stateOrder) {
    const RadiationHull sections = radiationHullFromMesh(hull, waterlineZ, stationCount, seaDensity);
    // The band a ship actually meets. Below 0.2 rad/s waves are swell the hull
    // simply follows; above 2.5 they are ripples against a ship's length, and the
    // 2D solve gets expensive exactly where the answer stops mattering.
    const std::vector<double> omega = radiationFrequencyGrid(0.2, 2.5, 40);
    const RadiationTable table = stripTheoryTable(sections, omega);

    // Sample the fit from the hull's own memory rather than from a guess. Ask K
    // for the heave entry over a generous window, find where it has decayed, then
    // fit over twice that -- the model has to reproduce the whole tail, because
    // truncating it is precisely what turns a damper into an integrator.
    constexpr double kFitStep = 0.1;   // s; K is smooth on this scale
    const std::vector<double> probe = retardationFunction(table, 2, 2, kFitStep, 1200);
    const double decay = memoryDecayTime(probe, kFitStep, 0.01);
    const int samples = std::clamp(static_cast<int>(2.0 * decay / kFitStep), 200, 1200);

    radiation.emplace(table, kFitStep, samples, stateOrder);
    return table;
}

RollDampingHull Ship::attachRollDamping(double waterlineZ, double bilgeKeelLength,
                                        double bilgeKeelBreadth) {
    RollDampingHull form = rollDampingHullFromMesh(hull, waterlineZ, bilgeKeelLength,
                                                   bilgeKeelBreadth, seaDensity);
    // The roll axis is loading rather than form; seed it from the ship as it
    // stands so a caller who never steps still gets a usable coefficient.
    // integrateRigidBody() refreshes it every tick thereafter.
    form.rollAxisAboveKeel = massProperties().cog.z - hullLo.z;

    // Left at zero deliberately, and not because it is unknown. B44W is the
    // radiation share of roll damping; when a RadiationForce is attached the
    // memory convolution is already applying it, so adding it here would be the
    // same double count that cost 27% of mid-frequency heave the first time
    // radiation was wired in. A caller with no radiation model who wants the
    // 5-30% ITTC puts on the wave component can set it themselves -- and
    // validateRollDamping() says so when it is missing.
    form.waveDamping = 0.0;

    rollDampingForm = form;
    return form;
}

void Ship::step(double dt, const Sea& sea) {
    updateInternalFreeSurfaces(sea);
    solveFlowNetwork(dt, sea);
    relaxGasToStructure(dt);
    updateInternalFreeSurfaces(sea);

    // TODO: Water promotion review (periodic, like structural/gas reviews would be)
    // - Call waterPromoter_.review(*this) at appropriate cadence
    // - Handle promoted list: create flip::Field via promotion::promoteWater, store in activeWaterFields_
    // - Handle demoted list: read state back via promotion::demoteWater, delete field, remove from map
    // - Step active FLIP fields: for each in activeWaterFields_, create flip::Solver, call solver.step(dt)
    // - Update compartment waterVolume from FLIP state when active

    integrateRigidBody(dt, sea);
}

// ---------------------------------------------------------------------------
// Stability analysis
// ---------------------------------------------------------------------------

double Ship::rightingArmAtHeel(double heelRad, const Sea& sea) const {
    // GZ depends only on hull form and mass distribution, not on where the sea
    // surface happens to sit in world coordinates.
    (void)sea;
    // Positive heel puts starboard (-y) down.
    const Quat q = Quat::fromAxisAngle(Vec3{1, 0, 0}, heelRad);
    const Mat3 R = q.toMat3();
    const Vec3 up = bodyFrameUp(R);

    // Re-level the floodwater at this attitude. Water volumes are conserved; only
    // their centroids move. This is the entire mechanism of free-surface effect.
    double mass = lightshipMass;
    Vec3 moment = lightshipCog * lightshipMass;
    for (const Compartment& c : compartments) {
        const double m = c.waterVolume * seaDensity;
        if (m <= 0) continue;
        const double region = c.waterVolume / std::max(c.permeability, 1e-6);
        const double off = solvePlaneOffsetForVolume(c.mesh, up, region, -kInf, kInf);
        const VolumeIntegral vi = integrateBelowPlane(c.mesh, up, off);
        mass += m;
        moment += (vi.volume > 1e-9 ? vi.centroid : (c.bboxLo + c.bboxHi) * 0.5) * m;
    }
    const Vec3 cogLocal = moment / mass;

    // Sink to the draft that balances this displacement, then read off B.
    const double offset = solvePlaneOffsetForVolume(hull, up, mass / seaDensity, -kInf, kInf);
    const VolumeIntegral sub = integrateBelowPlane(hull, up, offset);
    if (sub.volume <= 1e-9) return 0.0;

    // GZ is the transverse separation of B and G in the world horizontal.
    // B to starboard of G is a righting couple, hence the sign.
    const Vec3 bWorld = R * sub.centroid;
    const Vec3 gWorld = R * cogLocal;
    return -(bWorld.y - gWorld.y);
}

Diagnostics Ship::diagnostics(const Sea& sea) const {
    Diagnostics d;
    const Mat3 R = state.orientation.toMat3();
    const Vec3 up = bodyFrameUp(R);
    const double planeOffset = sea.level - state.position.z;

    // Displacement follows the actual free surface, so a ship on a crest reads a
    // larger displaced volume than the same ship in a trough. The stability
    // figures below deliberately stay still-water: GZ and GM are defined about a
    // mean waterline, and quoting a "GM in waves" would be inventing a quantity
    // naval architecture does not have.
    VolumeIntegral sub;
    if (sea.flat()) {
        sub = integrateBelowPlane(hull, up, planeOffset);
    } else {
        TriMesh worldHull = hull;
        for (Vec3& v : worldHull.verts) v = R * v + state.position;
        const VolumeIntegral world = integrateBelowSurface(
            worldHull, [&](double x, double y) { return sea.heightAt(x, y); });
        sub.volume = world.volume;
        // Report the centre of buoyancy in the body frame, as the flat path does.
        sub.centroid = R.transposed() * (world.centroid - state.position);
    }
    const MassProperties mp = massProperties();

    d.buoyantVolume = sub.volume;
    d.displacementMass = mp.mass;
    d.floodwaterMass = totalFloodwaterMass();
    d.centreOfGravity = mp.cog;
    d.centreOfBuoyancy = sub.centroid;
    heelTrimFromRotation(R, d.heelDeg, d.trimDeg);
    d.heelDeg *= kRadToDeg;
    d.trimDeg *= kRadToDeg;
    d.draftMidship = sea.level - (state.position.z + hullLo.z);

    const double h = 0.05;
    d.waterplaneArea = (integrateBelowPlane(hull, up, planeOffset + h).volume -
                        integrateBelowPlane(hull, up, planeOffset - h).volume) / (2 * h);

    double heelRad = 0, trimRad = 0;
    heelTrimFromRotation(R, heelRad, trimRad);
    d.gzRighting = rightingArmAtHeel(heelRad, sea);

    // GM is the slope of the righting arm at the origin, and a finite difference
    // only delivers that while the arm is linear across the angle it is taken at.
    //
    // **A shallow layer stops being linear almost immediately.** Tilt a compartment
    // holding a layer of mean depth h across a breadth b, and the surface stays a
    // full-breadth plane only while tan(heel) <= 2h/b; past that the water has
    // pulled off the high side, the surface is narrower than the compartment and
    // the free-surface moment collapses as `3/tau - 2/tau^1.5` in
    // `tau = tan(heel) b / 2h`. That closed form is exact for a box compartment and
    // holds on the ferry's own tapered vehicle deck to about 1%.
    //
    // 0.03 rad -- what this line used to be, unconditionally -- is ten times the
    // pocketing angle of a 2.9 cm layer on that deck, and it delivers 24% of the
    // free-surface effect. **The ship then reads +0.59 m where her initial GM is
    // -3.77 m**, which is the unsafe direction on the one number every stability
    // judgement here keys off.
    //
    // So the angle is measured rather than assumed: halve it until the slope stops
    // moving. The secant slope is monotone in the sampling angle and constant below
    // the pocketing angle, so "it stopped moving" is exactly "we are inside the
    // linear region" and the angle it settles at is a lower bound on where that
    // region ends. This is the only route that does not need a free-surface
    // *correction* term, a compartment idealised as a box, or a count of how many
    // compartments are wet -- none of which this model has anywhere else.
    //
    // Conditioning is not the trade it looks like. `rightingArmAtHeel` is
    // deterministic and its round-off cancels almost exactly across the difference:
    // measured, the ferry's GZ(0) is 6e-18 m intact and 5e-12 m with a free
    // surface, and GZ(eps) - GZ(-eps) carries about 5e-16 m of noise on the ferry
    // and on a box barge alike. The slope error is then ~5e-16/(2 eps) -- 2.7e-10 m
    // at the floor below, under three thousandths of the tolerance. There is no
    // sampling angle here that is wrong in both directions at once.
    const auto slopeAt = [&](double e) {
        return (rightingArmAtHeel(e, sea) - rightingArmAtHeel(-e, sea)) / (2 * e);
    };

    // 0.03 rad and no refinement at all when nothing is wet, on the same test
    // rightingArmAtHeel() itself uses to decide a compartment contributes. **That
    // is deliberate and it is a guarantee, not an optimisation**: an intact ship
    // has no free surface to pocket, so every figure this repository publishes
    // about one comes out bit-identical. What the intact ship is left with is the
    // hull's own O(eps^2) truncation -- +0.13% on the ferry, +0.34% on the box
    // barge -- which is a different and much smaller error, and moving it would
    // move every gated figure taken on an intact ship for no gain.
    bool anyFreeSurface = false;
    for (const Compartment& c : compartments)
        if (c.waterVolume * seaDensity > 0) { anyFreeSurface = true; break; }

    double eps = 0.03;
    double gm = slopeAt(eps);
    if (anyFreeSurface) {
        // **The bound is the round-off limit, not a policy about what is worth
        // reporting**, and getting that wrong is the one thing mutation testing
        // caught here. It was 15 halvings -- 9.2e-7 rad -- on the reasoning that
        // anything finer is a safety rail nothing reaches. That reasoning was
        // false: a *control* that raised the bound was killed, and the sweep
        // behind it shows why. Between 10 and 100 kg of water on the ferry's
        // vehicle deck the answer is perfectly well posed and 15 halvings could
        // not reach it -- 100 kg came out right and was flagged unreliable, and
        // 10 kg came out at -3.25 m against a true -3.78.
        //
        // 24 halvings bottom out at 1.79e-9 rad, where the noise above costs
        // 5e-16/(2 eps) = 1.4e-7 m of slope: seven times under the tolerance
        // floor, and the last power of two that stays there. Deeper starts
        // chasing round-off.
        //
        // The tolerance is what the slope has to stop moving by: a micrometre of
        // GM, or 1e-4 of it, whichever is larger.
        d.gmSlopeConverged = false;
        for (int i = 0; i < 24; ++i) {
            const double halved = 0.5 * eps;
            const double gmHalved = slopeAt(halved);
            if (std::abs(gmHalved - gm) <= std::max(1e-6, 1e-4 * std::abs(gm))) {
                // Report the finer value with the coarser angle: the two agree, so
                // linearity is established out to `eps`, and `gmHalved` is the
                // better-converged of the two numbers that establish it.
                gm = gmHalved;
                d.gmSlopeConverged = true;
                break;
            }
            eps = halved;
            gm = gmHalved;
            // Counted where the angle actually moves, so the count and the angle
            // cannot disagree: the pass that *stops* the search leaves both alone.
            ++d.gmHalvings;
        }
    }
    d.gmTransverse = gm;
    d.gmSampledAtRad = eps;

    // Lowest point of the weather deck edge relative to the sea.
    double minFreeboard = kInf;
    for (const Vec3& v : hull.verts) {
        if (v.z < deckEdgeZ - 1e-6) continue;
        const double wz = dot(up, v) + state.position.z;
        minFreeboard = std::min(minFreeboard, wz - sea.level);
    }
    d.freeboardMin = std::isfinite(minFreeboard) ? minFreeboard : 0.0;

    const double hullVolume = integrate(hull).volume;
    d.afloat = d.buoyantVolume < 0.995 * hullVolume && std::abs(d.heelDeg) < 90.0;
    return d;
}

// The convergence flag is read *before* the sign and not alongside it, because a
// slope that is not a metacentric height has no sign worth reading. See
// StabilityJudgement in ship.hpp for why the answer is "no reading" rather than a
// conservative "negative".
StabilityJudgement judgeStability(const Diagnostics& d) {
    if (!d.gmSlopeConverged) return StabilityJudgement::Unresolved;
    return d.gmTransverse < 0 ? StabilityJudgement::Negative : StabilityJudgement::Positive;
}

double FreeSurfaceLayer::momentFractionAt(double heelRad) const {
    if (depth <= 0 || breadth <= 0) return 0.0;
    const double tau = std::abs(std::tan(heelRad)) * breadth / (2.0 * depth);
    if (tau <= 1.0) return 1.0;
    return 3.0 / tau - 2.0 / (tau * std::sqrt(tau));
}

FreeSurfaceLayer largestFreeSurface(const Ship& s) {
    // The floor's plan area, from a 5 cm slab standing on it. One-sided *upward*
    // from `bboxLo.z`, so it is not measured through a slab half outside the space
    // -- and taken at the floor rather than at the water surface because that is
    // where a layer thin enough to matter actually sits.
    constexpr double kSlab = 0.05;

    FreeSurfaceLayer best;
    for (std::size_t i = 0; i < s.compartments.size(); ++i) {
        const Compartment& c = s.compartments[i];
        if (c.waterVolume * s.seaDensity <= 0) continue;
        const double area =
            integrate(clipToBox(c.mesh, c.bboxLo,
                                {c.bboxHi.x, c.bboxHi.y, c.bboxLo.z + kSlab})).volume / kSlab;
        const double length = c.bboxHi.x - c.bboxLo.x;
        if (area <= best.planArea || area <= 0 || length <= 0) continue;
        best.compartment = static_cast<int>(i);
        best.planArea = area;
        // Over the *geometric* region the water occupies, `V/mu`, which is what
        // sets the depth and hence the angle -- the same quantity the free-surface
        // moment in rightingArmAtHeel() is taken over.
        best.depth = c.waterVolume / std::max(c.permeability, 1e-6) / area;
        best.breadth = area / length;
        best.pocketingRad = std::atan(2.0 * best.depth / best.breadth);
    }
    return best;
}

}  // namespace sim
