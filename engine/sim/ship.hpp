// SPDX-License-Identifier: MIT
//
// Progressive flooding with trapped, compressible air.
//
// The model deliberately avoids the naval-architecture shortcuts (added-weight vs
// lost-buoyancy, tabulated free-surface corrections). Instead:
//
//   * floodwater is carried as real mass at the real centroid of the real water
//     body inside each compartment, re-levelled against gravity every tick. Free
//     surface effect is then not a correction term -- it is what happens when a
//     few hundred tonnes of water slides to the low side of a compartment.
//   * air is a compressible species with its own mass, pressure and flow paths.
//     A sealed compartment stops flooding when its air pressure balances the
//     external head, which is why upside-down hulls float for hours.
//   * every hole, door, hatch, vent and pipe is an orifice in one network. Water
//     and air move through the same edges; which phase moves depends only on
//     what is sitting at the opening on the high-pressure side.
#pragma once

#include "../core/geometry.hpp"
#include "waves.hpp"

#include <string>
#include <vector>

namespace sim {

// Sentinel compartment index meaning "the sea / open atmosphere outside the hull".
inline constexpr int kSea = -1;

// The water the ship floats on.
//
// Implicitly constructible from a still-water level, so every existing call site
// that passed a bare `double` keeps working and keeps taking the fast path. That
// matters: flat water is clipped by a single plane, which is exact and cheap,
// while a wave field costs a height evaluation per quadrature point.
struct Sea {
    Sea() = default;
    Sea(double stillWaterLevel) : level(stillWaterLevel) {}  // NOLINT: implicit on purpose

    double level = 0.0;               // still-water level, m
    const WaveField* waves = nullptr; // null means flat water
    double time = 0.0;                // s, for evaluating the wave field

    bool flat() const { return waves == nullptr; }
    double heightAt(double x, double y) const {
        return waves == nullptr ? level : level + waves->elevation(x, y, time);
    }
};

struct Compartment {
    std::string name;
    TriMesh mesh;                  // ship body frame
    double permeability = 0.95;    // void fraction; machinery spaces are much lower

    // State.
    double waterVolume = 0;        // m^3 of liquid water held
    double airMass = 0;            // kg of gas trapped
    bool   ventedToAtmosphere = false;  // open compartments never pressurise

    // Cached at initialise().
    double grossVolume = 0;        // m^3, geometric
    Vec3   bboxLo{}, bboxHi{};     // body frame, for the water body's own inertia

    // Derived each tick.
    double surfaceOffset = 0;      // plane offset of the internal free surface, body frame
    double surfaceWorldZ = 0;      // world height of that surface
    Vec3   waterCentroid{};        // body frame
    double airPressure = kPatm;    // Pa

    double floodableVolume() const { return grossVolume * permeability; }
    double airVolume()       const { return floodableVolume() - waterVolume; }
    double fillFraction()    const {
        const double f = floodableVolume();
        return f > 1e-9 ? waterVolume / f : 0.0;
    }
};

enum class OpeningKind { Breach, Door, Hatch, Vent, Pipe };

// An orifice between two spaces. Either endpoint may be kSea.
struct Opening {
    std::string name;
    int    a = kSea;
    int    b = kSea;
    Vec3   pos{};                  // body frame; the orifice centre
    double area = 0;               // m^2
    double dischargeCoeff = 0.6;   // sharp-edged hole ~0.6, smooth duct ~0.9
    bool   open = true;
    OpeningKind kind = OpeningKind::Breach;

    // Diagnostic, filled each tick.
    double lastFlow = 0;           // m^3/s, positive means a -> b
    bool   lastFlowWasWater = false;
};

// A bilge or ballast pump moving water out of a compartment, overboard.
struct Pump {
    std::string name;
    int    compartment = 0;
    double capacity = 0;           // m^3/s at zero head
    double maxHead = 20.0;         // m; output falls off as head approaches this
    bool   on = false;
    double lastFlow = 0;
};

struct RigidState {
    Vec3 position{};               // world position of the body origin
    Quat orientation{};            // body -> world
    Vec3 velocity{};               // world velocity of the body origin
    Vec3 angularVelocity{};        // world, rad/s
};

// Everything worth putting on an instrument panel.
struct Diagnostics {
    double displacementMass = 0;   // kg, total
    double floodwaterMass = 0;     // kg
    double buoyantVolume = 0;      // m^3 submerged
    double draftMidship = 0;       // m
    double heelDeg = 0;            // + starboard down
    double trimDeg = 0;            // + bow down
    double freeboardMin = 0;       // m, lowest point of the weather deck edge above sea
    double waterplaneArea = 0;     // m^2
    double gmTransverse = 0;       // m, effective (includes free-surface effect)
    double gzRighting = 0;         // m, righting arm at the current heel
    Vec3   centreOfGravity{};      // body frame
    Vec3   centreOfBuoyancy{};     // body frame
    bool   afloat = true;
};

class Ship {
public:
    // --- Configuration ---
    TriMesh hull;                        // watertight envelope, body frame
    std::vector<Compartment> compartments;
    std::vector<Opening> openings;
    std::vector<Pump> pumps;

    double lightshipMass = 0;            // kg, structure + machinery + outfit
    Vec3   lightshipCog{};               // body frame
    Vec3   gyradii{0, 0, 0};             // kxx, kyy, kzz in m, about lightshipCog

    double seaDensity = kRhoSeawater;
    double deckEdgeZ = 0;                // body-frame height of the weather deck edge
    Vec3   hullLo{}, hullHi{};           // hull bounding box, cached at initialise()

    // Damping as a fraction of critical, per mode. Roll is the lightly damped one
    // and is what makes a flooding ship feel alive rather than syrupy.
    double zetaHeave = 0.35;
    double zetaRoll  = 0.08;
    double zetaPitch = 0.30;

    // Added-mass coefficients as multiples of the displaced mass / inertia.
    double addedMassSurge = 0.05, addedMassSway = 0.9, addedMassHeave = 1.1;
    double addedInertiaRoll = 0.25, addedInertiaPitch = 1.0, addedInertiaYaw = 0.6;

    RigidState state;

    // --- Lifecycle ---

    // Caches gross volumes and puts the ship at its floating equilibrium.
    void initialise(const Sea& sea);

    // Advance the flooding network and the rigid body by dt seconds.
    void step(double dt, const Sea& sea);

    Diagnostics diagnostics(const Sea& sea) const;

    // Righting arm at a forced heel angle, with the ship free to sink to its
    // equilibrium draft and the floodwater free to re-level. Sweeping this gives
    // a true damaged-condition GZ curve.
    double rightingArmAtHeel(double heelRad, const Sea& sea) const;

    // Structural sanity of the ship definition, checked once after initialise().
    // Returns a human-readable problem per entry; empty means the definition is
    // self-consistent. Bad subdivision does not crash anything -- it quietly
    // produces wrong volumes and therefore a wrong ship -- so it is worth an
    // explicit check rather than trust.
    std::vector<std::string> validate() const;

    // Convenience accessors.
    int findCompartment(std::string_view name) const;
    double totalFloodwaterMass() const;

private:
    struct MassProperties {
        double mass = 0;
        Vec3   cog{};              // body frame
        Mat3   inertiaAboutCog{};  // body frame
    };

    void updateInternalFreeSurfaces(const Sea& sea);
    void solveFlowNetwork(double dt, const Sea& sea);
    void integrateRigidBody(double dt, const Sea& sea);

    MassProperties massProperties() const;

    // Pressure and phase presented by one side of an opening.
    struct SideState {
        double pressure = kPatm;
        bool   isWater = false;
    };
    SideState sideStateAt(int compartmentIndex, const Vec3& worldPos, const Sea& sea) const;

    // Sinkage/attitude-only hydrostatic solve, used by the GZ sweep.
    double equilibriumDraftAt(const Quat& orientation, const Sea& sea, double targetMass) const;

    // Reused each tick when the sea is wavy: the hull must be integrated in world
    // coordinates because the surface is a function of world x and y, and
    // reallocating a 600-vertex mesh per tick is pure waste.
    mutable TriMesh worldHullScratch_;

    double cachedWaterplaneArea_ = 0;
    double cachedKRoll_ = 0;
    double cachedKPitch_ = 0;
    int    stabilityRefreshCounter_ = 0;
};

}  // namespace sim
