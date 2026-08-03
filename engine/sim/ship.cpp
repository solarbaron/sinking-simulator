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

void Ship::initialise(const Sea& sea) {
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

Ship::SideState Ship::sideStateAt(int idx, const Vec3& worldPos, const Sea& sea) const {
    if (idx == kSea) {
        // The head at an opening is measured to the *local* free surface, not to
        // the still-water level: a breach under a crest floods faster than the
        // same breach under a trough, and that difference is a large part of why
        // damage in a seaway is worse than damage alongside.
        const double surface = sea.heightAt(worldPos.x, worldPos.y);
        if (worldPos.z < surface)
            return {kPatm + seaDensity * kGravity * (surface - worldPos.z), true};
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

void Ship::solveFlowNetwork(double dt, const Sea& sea) {
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
        const SideState sa = sideStateAt(o.a, worldPos, sea);
        const SideState sb = sideStateAt(o.b, worldPos, sea);

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
        const double head = std::max(0.0, sea.level - c.surfaceWorldZ);
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

        auto angularStiffness = [&](const Vec3& axisWorld) {
            const double eps = 2e-3;
            const Quat q2 = Quat::fromAxisAngle(axisWorld, eps) * state.orientation;
            const Mat3 R2 = q2.toMat3();
            const VolumeIntegral s2 = PlaneSweep(hull, bodyFrameUp(R2)).below(planeOffset);
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
    Vec3 mEff{mp.mass + addedMassSurge * dispMass,
              mp.mass + addedMassSway  * dispMass,
              mp.mass + addedMassHeave * dispMass};

    Mat3 Ieff = mp.inertiaAboutCog;
    Ieff(0, 0) *= (1.0 + addedInertiaRoll);
    Ieff(1, 1) *= (1.0 + addedInertiaPitch);
    Ieff(2, 2) *= (1.0 + addedInertiaYaw);

    // With a radiation model, the guessed coefficients give way to this hull's
    // own infinite-frequency added mass. Only the diagonal is taken: the
    // integrator below inverts a 3x3 inertia and divides by a per-axis mass, so
    // there is nowhere for the sway-roll coupling A_24 to go. That coupling is
    // real and matters for roll, and taking the diagonal is an approximation this
    // integrator's shape forces rather than one the physics justifies -- see
    // docs/02-simulation.md.
    if (radiation.has_value()) {
        const Matrix6& aInf = radiation->addedMassInfinite();
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
    // Roll and yaw keep theirs. Roll radiation damping is genuinely small and the
    // mechanism that matters is viscous -- eddies and bilge keels, which
    // roll_damping.{hpp,cpp} computes by Ikeda's method and which is not yet
    // wired into this integrator; deleting roll's stand-in before that lands
    // would leave the one mode that most needs damping with almost none. Surge
    // keeps its damper too, since strip theory contributes no surge radiation at
    // all. Quadratic drag is untouched throughout: it is viscous and separate.
    if (radiation.has_value()) {
        cLin.y = 0.0;
        cLin.z = 0.0;
        cAng.y = 0.0;
    }

    // Quadratic drag on the projected areas of the box around the hull.
    const Vec3 ext = hullHi - hullLo;
    const double draft = std::max(0.1, sea.level - (state.position.z + hullLo.z));
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

    // The wave-memory force: waves this ship radiated earlier, still nearby and
    // still pushing back. Advanced with the body-frame velocity and subtracted,
    // per the sign convention in radiation.hpp.
    if (radiation.has_value()) {
        radiation->step({vBody.x, vBody.y, vBody.z, wBody.x, wBody.y, wBody.z}, dt);
        const std::array<double, 6> memory = radiation->memoryForce();
        for (int i = 0; i < 3; ++i) {
            fBody[i] -= memory[i];
            tBody[i] -= memory[i + 3];
        }
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

void Ship::step(double dt, const Sea& sea) {
    updateInternalFreeSurfaces(sea);
    solveFlowNetwork(dt, sea);
    updateInternalFreeSurfaces(sea);
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
    const double eps = 0.03;
    d.gmTransverse = (rightingArmAtHeel(eps, sea) - rightingArmAtHeel(-eps, sea))
                     / (2 * eps);

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

}  // namespace sim
