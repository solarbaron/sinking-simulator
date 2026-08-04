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
#include "propulsion.hpp"
#include "radiation.hpp"
#include "roll_damping.hpp"
#include "waves.hpp"

#include <optional>
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
    //
    // Every one of these is a lumped stand-in for a mechanism that is now
    // modelled somewhere: heave/sway/pitch for radiation, yaw for N_r, roll for
    // Ikeda's viscous components. They are the default because the flooding
    // scenarios predate all of it and are validated against their own behaviour;
    // each is switched off in integrateRigidBody() the moment the real model
    // arrives, which is not a tidy-up but the entire point -- adding the real
    // thing alongside its own stand-in damps the ship twice.
    double zetaHeave = 0.35;
    double zetaRoll  = 0.08;
    double zetaPitch = 0.30;

    // Added-mass coefficients as multiples of the displaced mass / inertia.
    //
    // These are plausible textbook figures, not this hull's. They are used only
    // while no radiation model is attached; attachRadiation() replaces them with
    // A_inf computed from the hull's own sections. Keeping both means the
    // flooding scenarios, which predate radiation and are validated against
    // their own behaviour, are untouched by its arrival.
    double addedMassSurge = 0.05, addedMassSway = 0.9, addedMassHeave = 1.1;
    double addedInertiaRoll = 0.25, addedInertiaPitch = 1.0, addedInertiaYaw = 0.6;

    // Frequency-dependent radiation, when the ship has been given it. Held by
    // value, not shared: the memory states are part of *this* ship's history, and
    // two copies riding the same states would be two ships pretending to be one.
    std::optional<RadiationForce> radiation;

    // Propulsion, steering and horizontal-plane hydrodynamics. Optional, like
    // radiation: a drifting casualty has none, and that is the case Phase 0 was
    // built around.
    //
    // `Manoeuvring` carries its own ManoeuvringState, and this ship does *not*
    // let it integrate. Ship owns the motion; the state is overwritten from the
    // rigid body every tick and only ever read back as forces. Anything else
    // would be two integrators disagreeing about where the ship is.
    std::optional<Manoeuvring> propulsion;

    // Viscous roll damping by Ikeda's method, when the ship has been given it.
    // Potential flow gives roll almost no damping at all, so this is where roll
    // damping actually comes from, and it replaces zetaRoll outright rather than
    // adding to it.
    //
    // Held as the hull-form parameters rather than as a coefficient, because
    // Ikeda's B44 is a function of the operating point -- roll amplitude, roll
    // frequency, forward speed -- and the operating point changes every tick.
    // integrateRigidBody() evaluates it per tick at 52 ns against a 16 us tick,
    // which is cheap enough that caching it would only buy a stale answer.
    //
    // `rollAxisAboveKeel` is overwritten from the live centre of gravity each
    // tick, on the same terms as the manoeuvring state: hull *form* is fixed and
    // is snapshotted by attachRollDamping(), but the roll axis is a property of
    // loading and a flooding ship's loading moves.
    //
    // **waveDamping must stay zero while `radiation` is attached.** It is the
    // radiation share of B44, and with a radiation model the memory convolution
    // is already applying it; supplying both is the same double count that cost
    // 27% of mid-frequency heave when radiation first landed.
    std::optional<RollDampingHull> rollDampingForm;

    // Diagnostic, filled each tick when rollDampingForm is set: the operating
    // point Ikeda was asked at and every component of the answer. Which mechanism
    // is doing the work is the whole diagnostic value of the method, so it is
    // published rather than consumed silently.
    RollDampingCondition rollCondition{};
    RollDamping rollDampingApplied{};

    // Force and moment from something that is not the sea: today, contact with
    // another ship's hull. World frame; the moment is taken about the centre of
    // gravity, which is the point integrateRigidBody() takes every other moment
    // about. Added to by whoever is applying the load, consumed and **cleared** by
    // the next step(), so a caller that stops pushing gets a ship that stops being
    // pushed rather than one that coasts on a stale force.
    //
    // Unlike radiation, Ikeda and the MMG polynomial, this replaces nothing.
    // Every damping and added-mass term above is a *fluid* mechanism that acts on
    // a ship with no other ship anywhere near it, so none of them was ever a
    // lumped stand-in for contact and none of them is switched off while contact
    // is happening. See engine/sim/collision.hpp §4.
    Vec3 externalForce{};
    Vec3 externalMoment{};

    RigidState state;

    // --- Lifecycle ---

    // Caches gross volumes and puts the ship at its floating equilibrium.
    void initialise(const Sea& sea);

    // Advance the flooding network and the rigid body by dt seconds.
    void step(double dt, const Sea& sea);

    // Build a radiation model from this hull's own sections and attach it, so
    // added mass and the wave-memory force stop being guesses. `waterlineZ` is
    // the still-water plane in body coordinates; call after initialise().
    //
    // Deliberately takes no timestep. The state-space fit is sampled from the
    // *ship's memory*, which is a property of the hull -- tens of seconds -- and
    // has nothing to do with how finely the simulation is stepped. An earlier
    // version passed the simulation dt straight through and fitted K(t) over
    // 1.3 s of a 20 s decay: the fit then sees a nearly constant K, places its
    // poles near zero, and the resulting model integrates velocity instead of
    // damping it. The ship diverged to NaN in five steps. RadiationForce::step()
    // is an exact zero-order hold and genuinely does not care what dt it is
    // driven at, which is what makes this separation possible.
    //
    // Returns the table it built, because worstEnergyResidual and repairedSolves
    // are the numbers that say whether to trust it, and swallowing them here
    // would be the wrong kind of convenience.
    RadiationTable attachRadiation(double waterlineZ, int stationCount = 21, int stateOrder = 6);

    // Read this hull's form off its own mesh and attach Ikeda's viscous roll
    // damping, so the roll mode stops running on a guessed fraction of critical.
    // `waterlineZ` is the still-water plane in body coordinates; call after
    // initialise(). Bilge keels are arguments because a watertight envelope has
    // no appendages in it -- see rollDampingHullFromMesh().
    //
    // Returns the form it derived, because validateRollDamping() against it is
    // what says whether this hull is inside the domain the regressions were
    // fitted over, and swallowing that would be the wrong kind of convenience.
    RollDampingHull attachRollDamping(double waterlineZ, double bilgeKeelLength = 0,
                                      double bilgeKeelBreadth = 0);

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

    // Mass, centre of gravity and inertia as the rigid body sees them, floodwater
    // included. Public because the contact solver needs exactly the numbers the
    // integrator uses, and a second implementation of them would be a second
    // answer that agrees until the day it does not.
    struct MassProperties {
        double mass = 0;
        Vec3   cog{};              // body frame
        Mat3   inertiaAboutCog{};  // body frame
    };
    MassProperties massProperties() const;

private:
    void updateInternalFreeSurfaces(const Sea& sea);
    void solveFlowNetwork(double dt, const Sea& sea);
    void integrateRigidBody(double dt, const Sea& sea);

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
