// SPDX-License-Identifier: MIT
#include "ship.hpp"

#include <limits>
#include <string_view>

namespace sim {
namespace {

constexpr double kInf = std::numeric_limits<double>::infinity();

// Ship "up" expressed in the body frame. Every free surface -- the sea outside and
// the floodwater inside -- is a plane with this normal.
inline Vec3 bodyFrameUp(const Mat3& R) { return R.transposed() * Vec3{0, 0, 1}; }

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

void Ship::initialise(double seaLevel) {
    boundingBox(hull, hullLo, hullHi);

    for (Compartment& c : compartments) {
        c.grossVolume = integrate(c.mesh).volume;
        boundingBox(c.mesh, c.bboxLo, c.bboxHi);
        c.waterVolume = std::clamp(c.waterVolume, 0.0, c.floodableVolume());
        // Start dry compartments full of air at atmospheric pressure.
        c.airMass = kPatm * std::max(c.airVolume(), 0.0) / (kRAir * kTAmbient);
        c.airPressure = kPatm;
    }

    // Drop the hull to its floating waterline before the first tick, so the sim
    // does not open with a large transient heave.
    state.position.z = 0;
    updateInternalFreeSurfaces(seaLevel);
    const MassProperties mp = massProperties();
    state.position.z = equilibriumDraftAt(state.orientation, seaLevel, mp.mass);
    updateInternalFreeSurfaces(seaLevel);
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

void Ship::updateInternalFreeSurfaces(double seaLevel) {
    (void)seaLevel;
    const Mat3 R = state.orientation.toMat3();
    const Vec3 up = bodyFrameUp(R);

    for (Compartment& c : compartments) {
        c.waterVolume = std::clamp(c.waterVolume, 0.0, c.floodableVolume());

        // Permeability: the liquid occupies a geometric region larger than its own
        // volume, because part of that region is structure and cargo it flows around.
        const double geometricRegion = c.waterVolume / std::max(c.permeability, 1e-6);
        c.surfaceOffset = solvePlaneOffsetForVolume(c.mesh, up, geometricRegion, -kInf, kInf);
        c.surfaceWorldZ = c.surfaceOffset + state.position.z;

        const VolumeIntegral vi = integrateBelowPlane(c.mesh, up, c.surfaceOffset);
        c.waterCentroid = vi.volume > 1e-9 ? vi.centroid
                                           : (c.bboxLo + c.bboxHi) * 0.5;

        // Isothermal, not adiabatic: a few tonnes of steel bulkhead is an enormous
        // heat sink next to the air in a compartment, so trapped air tracks ambient
        // temperature far faster than a compartment floods.
        if (c.ventedToAtmosphere) {
            c.airPressure = kPatm;
            c.airMass = kPatm * std::max(c.airVolume(), 0.0) / (kRAir * kTAmbient);
        } else {
            const double va = std::max(c.airVolume(), 1e-3 * std::max(c.grossVolume, 1.0));
            c.airPressure = c.airMass * kRAir * kTAmbient / va;
        }
    }
}

Ship::SideState Ship::sideStateAt(int idx, const Vec3& worldPos, double seaLevel) const {
    if (idx == kSea) {
        if (worldPos.z < seaLevel)
            return {kPatm + seaDensity * kGravity * (seaLevel - worldPos.z), true};
        return {kPatm, false};
    }
    const Compartment& c = compartments[static_cast<std::size_t>(idx)];
    if (c.waterVolume > 1e-9 && worldPos.z < c.surfaceWorldZ)
        return {c.airPressure + seaDensity * kGravity * (c.surfaceWorldZ - worldPos.z), true};
    return {c.airPressure, false};
}

// ---------------------------------------------------------------------------
// Orifice network
// ---------------------------------------------------------------------------

void Ship::solveFlowNetwork(double dt, double seaLevel) {
    const Mat3 R = state.orientation.toMat3();
    const auto n = compartments.size();

    // Deltas are accumulated against the start-of-tick state so that opening order
    // does not bias the result -- an explicit Jacobi sweep over the network.
    std::vector<double> dWater(n, 0.0), dAir(n, 0.0);

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

    for (Opening& o : openings) {
        o.lastFlow = 0;
        if (!o.open || o.area <= 0) continue;

        const Vec3 worldPos = R * o.pos + state.position;
        const SideState sa = sideStateAt(o.a, worldPos, seaLevel);
        const SideState sb = sideStateAt(o.b, worldPos, seaLevel);

        const double dp = sa.pressure - sb.pressure;
        if (std::abs(dp) < 1e-3) continue;  // below the noise floor of the solve

        const bool aIsDonor = dp > 0;
        const int donor = aIsDonor ? o.a : o.b;
        const int recv  = aIsDonor ? o.b : o.a;
        // What actually moves is whatever is sitting against the orifice on the
        // high-pressure side. This one line is why a hole above the internal
        // waterline vents air, and the same hole below it admits water.
        const bool water = aIsDonor ? sa.isWater : sb.isWater;

        const double pDonor = std::max(sa.pressure, sb.pressure);
        const double rho = water ? seaDensity : pDonor / (kRAir * kTAmbient);

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
            const double pDonorSide = aIsDonor ? sa.pressure : sb.pressure;
            const double pRecvSide  = aIsDonor ? sb.pressure : sa.pressure;
            if (recv != kSea) {
                const Compartment& rc = compartments[recv];
                const double va = std::max(rc.airVolume(), 1e-3 * std::max(rc.grossVolume, 1.0));
                const double mEq = pDonorSide * va / (kRAir * kTAmbient);
                dm = std::min(dm, std::max(0.0, mEq - (rc.airMass + dAir[recv])));
            }
            if (donor != kSea) {
                const Compartment& dc = compartments[donor];
                const double va = std::max(dc.airVolume(), 1e-3 * std::max(dc.grossVolume, 1.0));
                const double mEq = pRecvSide * va / (kRAir * kTAmbient);
                dm = std::min(dm, std::max(0.0, (dc.airMass + dAir[donor]) - mEq));
            }

            if (donor != kSea) dAir[donor] -= dm;
            if (recv  != kSea) dAir[recv]  += dm;
            o.lastFlow = (aIsDonor ? 1.0 : -1.0) * (dt > 0 ? (rho > 0 ? dm / rho / dt : 0.0) : 0.0);
        }
        o.lastFlowWasWater = water;
    }

    for (Pump& p : pumps) {
        p.lastFlow = 0;
        if (!p.on) continue;
        Compartment& c = compartments[static_cast<std::size_t>(p.compartment)];
        // Centrifugal pumps lose output as discharge head rises; overboard discharge
        // means lifting from the compartment's water surface to the sea surface.
        const double head = std::max(0.0, seaLevel - c.surfaceWorldZ);
        const double efficiency = std::clamp(1.0 - head / std::max(p.maxHead, 1e-6), 0.0, 1.0);
        const double q = p.capacity * efficiency;
        const double dv = std::min(q * dt, std::max(c.waterVolume + dWater[p.compartment], 0.0));
        dWater[p.compartment] -= dv;
        p.lastFlow = dt > 0 ? dv / dt : 0.0;
    }

    for (std::size_t i = 0; i < n; ++i) {
        compartments[i].waterVolume =
            std::clamp(compartments[i].waterVolume + dWater[i], 0.0,
                       compartments[i].floodableVolume());
        compartments[i].airMass = std::max(0.0, compartments[i].airMass + dAir[i]);
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

double Ship::equilibriumDraftAt(const Quat& orientation, double seaLevel,
                                double targetMass) const {
    const Vec3 up = bodyFrameUp(orientation.toMat3());
    const double offset =
        solvePlaneOffsetForVolume(hull, up, targetMass / seaDensity, -kInf, kInf);
    // world_z(x) = dot(up, x) + position.z, so the surface sits at offset + position.z.
    return seaLevel - offset;
}

void Ship::integrateRigidBody(double dt, double seaLevel) {
    const Mat3 R = state.orientation.toMat3();
    const Vec3 up = bodyFrameUp(R);
    const double planeOffset = seaLevel - state.position.z;

    const VolumeIntegral sub = integrateBelowPlane(hull, up, planeOffset);
    const MassProperties mp = massProperties();
    const Vec3 cogWorld = R * mp.cog + state.position;

    Vec3 force{0, 0, -mp.mass * kGravity};
    Vec3 torque{0, 0, 0};  // about cogWorld

    Vec3 buoyancy{0, 0, 0};
    if (sub.volume > 0) {
        buoyancy = Vec3{0, 0, seaDensity * kGravity * sub.volume};
        const Vec3 cbWorld = R * sub.centroid + state.position;
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
        cachedWaterplaneArea_ =
            (integrateBelowPlane(hull, up, planeOffset + h).volume -
             integrateBelowPlane(hull, up, planeOffset - h).volume) / (2 * h);

        auto angularStiffness = [&](const Vec3& axisWorld) {
            const double eps = 2e-3;
            const Quat q2 = Quat::fromAxisAngle(axisWorld, eps) * state.orientation;
            const Mat3 R2 = q2.toMat3();
            const VolumeIntegral s2 = integrateBelowPlane(hull, bodyFrameUp(R2), planeOffset);
            const Vec3 fb2{0, 0, seaDensity * kGravity * s2.volume};
            const Vec3 cb2 = R2 * s2.centroid + state.position;
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
    const Vec3 mEff{mp.mass + addedMassSurge * dispMass,
                    mp.mass + addedMassSway  * dispMass,
                    mp.mass + addedMassHeave * dispMass};

    Mat3 Ieff = mp.inertiaAboutCog;
    Ieff(0, 0) *= (1.0 + addedInertiaRoll);
    Ieff(1, 1) *= (1.0 + addedInertiaPitch);
    Ieff(2, 2) *= (1.0 + addedInertiaYaw);

    const Vec3 vBody = R.transposed() * state.velocity;
    const Vec3 wBody = R.transposed() * state.angularVelocity;

    const double kHeave = seaDensity * kGravity * std::max(cachedWaterplaneArea_, 1.0);
    const Vec3 cLin{2 * zetaHeave * std::sqrt(kHeave * mEff.x) * 0.05,
                    2 * zetaHeave * std::sqrt(kHeave * mEff.y) * 0.30,
                    2 * zetaHeave * std::sqrt(kHeave * mEff.z)};
    const Vec3 cAng{2 * zetaRoll  * std::sqrt(cachedKRoll_  * Ieff(0, 0)),
                    2 * zetaPitch * std::sqrt(cachedKPitch_ * Ieff(1, 1)),
                    0.02 * Ieff(2, 2)};

    // Quadratic drag on the projected areas of the box around the hull.
    const Vec3 ext = hullHi - hullLo;
    const double draft = std::max(0.1, seaLevel - (state.position.z + hullLo.z));
    const Vec3 projArea{ext.y * draft, ext.x * draft, ext.x * ext.y};
    const Vec3 dragCoeff{0.10, 1.00, 1.50};

    Vec3 fBody = R.transposed() * force;
    Vec3 tBody = R.transposed() * torque;
    for (int i = 0; i < 3; ++i) {
        const double v = vBody[i];
        fBody[i] -= cLin[i] * v
                  + 0.5 * seaDensity * dragCoeff[i] * projArea[i] * std::abs(v) * v;
        tBody[i] -= cAng[i] * wBody[i];
    }

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

void Ship::step(double dt, double seaLevel) {
    updateInternalFreeSurfaces(seaLevel);
    solveFlowNetwork(dt, seaLevel);
    updateInternalFreeSurfaces(seaLevel);
    integrateRigidBody(dt, seaLevel);
}

// ---------------------------------------------------------------------------
// Stability analysis
// ---------------------------------------------------------------------------

double Ship::rightingArmAtHeel(double heelRad, double seaLevel) const {
    // GZ depends only on hull form and mass distribution, not on where the sea
    // surface happens to sit in world coordinates.
    (void)seaLevel;
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

Diagnostics Ship::diagnostics(double seaLevel) const {
    Diagnostics d;
    const Mat3 R = state.orientation.toMat3();
    const Vec3 up = bodyFrameUp(R);
    const double planeOffset = seaLevel - state.position.z;

    const VolumeIntegral sub = integrateBelowPlane(hull, up, planeOffset);
    const MassProperties mp = massProperties();

    d.buoyantVolume = sub.volume;
    d.displacementMass = mp.mass;
    d.floodwaterMass = totalFloodwaterMass();
    d.centreOfGravity = mp.cog;
    d.centreOfBuoyancy = sub.centroid;
    heelTrimFromRotation(R, d.heelDeg, d.trimDeg);
    d.heelDeg *= kRadToDeg;
    d.trimDeg *= kRadToDeg;
    d.draftMidship = seaLevel - (state.position.z + hullLo.z);

    const double h = 0.05;
    d.waterplaneArea = (integrateBelowPlane(hull, up, planeOffset + h).volume -
                        integrateBelowPlane(hull, up, planeOffset - h).volume) / (2 * h);

    double heelRad = 0, trimRad = 0;
    heelTrimFromRotation(R, heelRad, trimRad);
    d.gzRighting = rightingArmAtHeel(heelRad, seaLevel);
    const double eps = 0.03;
    d.gmTransverse = (rightingArmAtHeel(eps, seaLevel) - rightingArmAtHeel(-eps, seaLevel))
                     / (2 * eps);

    // Lowest point of the weather deck edge relative to the sea.
    double minFreeboard = kInf;
    for (const Vec3& v : hull.verts) {
        if (v.z < deckEdgeZ - 1e-6) continue;
        const double wz = dot(up, v) + state.position.z;
        minFreeboard = std::min(minFreeboard, wz - seaLevel);
    }
    d.freeboardMin = std::isfinite(minFreeboard) ? minFreeboard : 0.0;

    const double hullVolume = integrate(hull).volume;
    d.afloat = d.buoyantVolume < 0.995 * hullVolume && std::abs(d.heelDeg) < 90.0;
    return d;
}

}  // namespace sim
