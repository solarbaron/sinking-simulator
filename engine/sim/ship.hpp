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
//     external head, which is why upside-down hulls float for hours -- and it
//     does *not* stop when the opening is a horizontal one, because there the air
//     has a way out the water is not blocking. See below.
//   * every hole, door, hatch, vent and pipe is an orifice in one network. Water
//     and air move through the same edges; which phase moves depends only on
//     what is sitting at the opening on the high-pressure side -- except through a
//     horizontal one, where both phases move at once and in opposite directions.
//
// **A horizontal opening is not a small vertical one.** Through a hatch in a deck
// with the sea standing on it, water falls while air rises through the *same*
// hole, driven by the density difference and not by any net pressure difference
// at all. A model that takes one dp at the orifice centre is at rest in exactly
// that case and moves nothing, which is the identical failure the fire work found
// in a doorway: the steady state of a doorway is hot gas out of the top and cool
// air into the bottom in equal mass flows, and a net-only model sees a still one.
// `fire.cpp` fixed the doorway by integrating over its height. A hatch has no
// height to integrate over, so `Ship::solveFlowNetwork` splits its *area* instead
// -- see the derivation there.
#pragma once

#include "../core/geometry.hpp"
#include "propulsion.hpp"
#include "radiation.hpp"
#include "roll_damping.hpp"
#include "water_promotion.hpp"
#include "waves.hpp"

#include <map>
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
    double gasTemperature = kTAmbient;  // K, well mixed over the gas space
    bool   ventedToAtmosphere = false;  // open compartments never pressurise

    // How fast this compartment's gas gives its heat up to the structure around
    // it, as a relaxation time towards kTAmbient. **Zero is the default and it
    // means instantaneous**, which reproduces the isothermal gas this model
    // carried before temperature was a state at all: a few tonnes of steel
    // bulkhead is an enormous heat sink next to the air in a compartment, so
    // trapped air tracks ambient far faster than a compartment floods, and
    // Boyle's law is then the right pressure law for flooding.
    //
    // That default is not laziness, it is what the flooding scenarios were
    // validated under -- `testTrappedAirArrestsFlooding` asserts pV is conserved
    // through a compression to 2%, which is a statement that the work of
    // compression goes into the steel and not into the gas. It is also what makes
    // a cold ship bit-identical to the model without temperature: at zero the
    // temperature is pinned to kTAmbient by exact assignment, and every
    // expression below that reads it recovers its old form character for
    // character.
    //
    // Set it positive and the gas becomes a real thermodynamic system: pdV work
    // heats it under compression, flow carries enthalpy rather than just mass,
    // and density -- hence buoyancy through a vertical opening -- follows the
    // temperature. Infinity is the adiabatic limit and is where the isentropic
    // closed form T/T0 = (p/p0)^((gamma-1)/gamma) is exact. A compartment with a
    // fire in it wants this; a cold one has no use for it.
    //
    // This is the same idiom as `radiation` and `rollDampingForm`: the real model
    // arrives alongside the lumped stand-in and replaces it only where it is
    // asked for, so nothing that was validated against the stand-in moves.
    double gasThermalTime = 0;     // s

    // Cached at initialise().
    double grossVolume = 0;        // m^3, geometric
    Vec3   bboxLo{}, bboxHi{};     // body frame, for the water body's own inertia

    // Derived each tick.
    double surfaceOffset = 0;      // plane offset of the internal free surface, body frame
    double surfaceWorldZ = 0;      // world height of that surface
    Vec3   waterCentroid{};        // body frame
    double airPressure = kPatm;    // Pa

    // World height of the middle of the gas space -- the free surface when there
    // is water, the bottom of the space when there is not, up to the deckhead.
    // This is the datum `airPressure` is quoted at: a well-mixed gas is only
    // uniform in pressure to the extent its own weight is negligible, and the
    // whole mechanism of smoke spread is that for a hot gas it is not. Openings
    // above this height see less than airPressure, below it more.
    double gasCentroidWorldZ = 0;

    // Gas volume at the end of the previous free-surface update, so the pdV work
    // of the water rising against the trapped air can be taken between one update
    // and the next. Only read when gasThermalTime is non-zero.
    double lastAirVolume = 0;      // m^3

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

    // The gas side of the same flow, which a volumetric rate cannot express on
    // its own: `lastFlow` is quoted at the donor's density, and the donor's
    // density is now a function of its temperature. Mass and the temperature it
    // arrives at are what a species or energy account downstream actually needs,
    // and deriving them back out of a volume would mean re-deriving the donor.
    // Zero on a tick that moved water, or nothing.
    double lastGasMassFlow = 0;         // kg/s, positive means a -> b
    double lastGasDonorTemperature = 0; // K, the temperature the gas left at

    // The enthalpy riding on that mass, in the kg*K units the network's own energy
    // accumulator uses -- `gamma * T_donor * mdot`, positive a -> b.
    //
    // For every opening that moves one stream this is exactly
    // `kGammaAir * lastGasDonorTemperature * lastGasMassFlow` and is redundant. It
    // exists for the one that moves *two*: a horizontal opening exchanging gas
    // both ways has two donors at two temperatures, and no single donor
    // temperature describes it -- `lastGasDonorTemperature` is then zero and this
    // is the only field that still adds up. An energy account should read this.
    double lastGasEnthalpyFlow = 0;     // kg*K/s, positive means a -> b

    // --- What a net rate cannot say ------------------------------------------
    //
    // Through a horizontal opening the heavy fluid falls while the light one rises
    // through the same hole at once -- in *equal* volumes when nothing else is
    // pushing, which is the case the mechanism exists for and is exactly the case a
    // net rate calls a hatch at rest. Reporting only `lastFlow` would report a
    // hatch doing nothing while a tonne a second went through it. These four say
    // what actually crossed, and their difference is the net `lastFlow` carries.
    //
    // All zero on any opening that is not exchanging, which is every opening that
    // is not a `Hatch` and every hatch whose stratification is stable or whose net
    // head has already flushed it one way.
    double lastExchangeDown = 0;        // m^3/s falling, out of the upper side
    double lastExchangeUp = 0;          // m^3/s rising, out of the lower side
    double lastExchangeMassDown = 0;    // kg/s falling
    double lastExchangeMassUp = 0;      // kg/s rising
    int    lastExchangeUpper = kSea;    // the endpoint the falling stream left, when one ran
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

    // The half-angle `gmTransverse` was central-differenced over. A metacentric
    // height is the slope of the righting arm *at the origin*, so this is the
    // statement that makes the number one: the arm was measured to be linear
    // across +/- this heel, and nothing is claimed beyond it.
    //
    // It is 0.03 rad on a ship with no free surface and shrinks -- by a factor of
    // sixteen on the ferry with fifty tonnes on her vehicle deck -- when a shallow
    // layer pockets. See Ship::diagnostics().
    double gmSampledAtRad = 0;

    // How many halvings it took to get there, so that `gmSampledAtRad` can be read
    // as a *statement about the refinement* and not just as an angle: zero on a
    // ship with no free surface, and 24 -- the bound -- on one where the search
    // ran to the floor. It is exactly `log2(0.03 / gmSampledAtRad)`, and it is
    // reported rather than left to be derived because the two ways of saying it
    // must not be able to drift apart.
    int gmHalvings = 0;

    // False when that refinement ran out of room: the free surface pockets at an
    // angle the righting arm cannot be resolved across, so `gmTransverse` is the
    // best slope available rather than a metacentric height, and a consumer that
    // keys a stability judgement off GM needs to tell "negative" from
    // "unanswerable".
    //
    // A free-surface moment does not depend on how much water there is, only on
    // the shape of its surface -- while the angle that surface survives to does.
    // So the initial GM is *discontinuous* in the amount of loose water, and the
    // ferry's vehicle deck straddles that discontinuity within a factor of two:
    // half a litre reports her own **+2.00 m** over +/-3.8e-3 rad, because a layer
    // that thin is not there at any angle anything can measure; one litre reports
    // no usable window at all and sets this false. Both are honest and the flag is
    // what tells them apart.
    bool gmSlopeConverged = true;

    double gzRighting = 0;         // m, righting arm at the current heel
    Vec3   centreOfGravity{};      // body frame
    Vec3   centreOfBuoyancy{};     // body frame
    bool   afloat = true;
};

// What transverse stability a set of diagnostics will support being asked about.
//
// **`Unresolved` is not a degree of danger.** It is the statement that
// `gmTransverse` is not a metacentric height, so *neither* sign is available from
// it -- and that is a different thing from a GM that came out near zero. The
// discontinuity documented on `gmSlopeConverged` is what makes it different: on
// the ferry's vehicle deck half a litre reports +2.00 m and one tonne reports
// -3.78 m, both converged, and the litre between them resolves to nothing at all.
// A consumer handed that litre and reading only the sign gets whatever the
// refinement happened to be holding, which is how `--scenario=full` published
// SURVIVED off a 6 mm puddle.
//
// **So an unresolved GM must not be answered with a conservative "negative"
// either**, which is the tempting reading and is wrong on the measurement. The
// condition arises when the layer is thin enough that its free surface is not
// there at any angle the arm can be sampled across: microns, in every case
// observed. In `--scenario=none` it is three steps inside a 0.12 s window at
// t+679.4, while the first 1.2 to 4.8 kg of water reach the vehicle deck -- a ship
// in no more danger there than one step earlier. Defaulting to "negative" would
// cry wolf precisely where the water is negligible, and would say nothing new
// where it is not: the 50 t layer that started all of this converges perfectly
// well and is genuinely -3.77 m. The honest answer is that the instrument has no
// reading, and the caller has to decide what to do about that rather than be
// handed a sign it cannot support.
enum class StabilityJudgement {
    Unresolved,   // gmSlopeConverged is false: no metacentric height exists to judge
    Negative,     // GM < 0: she has no initial stability and will loll or worse
    Positive,     // GM >= 0
};

// The one place that rule lives. Two tools drew a survival verdict from their own
// copy of `gmTransverse < 0`, neither reachable from the test suite, and both
// scored a flooded ship exactly as confidently as a dry one.
StabilityJudgement judgeStability(const Diagnostics& d);

// The shallow-layer geometry of one flooded compartment: what the wedge closed
// form needs, measured off the compartment's own mesh rather than assumed.
//
// **Reported, never used.** Nothing in the flooding solution reads any of this --
// the free surface is a re-levelled water body and not a correction term, which is
// the whole point of the model. It exists because the *angle a metacentric height
// stops meaning anything above* is a published claim about this ship, and until
// now the only way to re-derive it was to write C++ against the library.
struct FreeSurfaceLayer {
    int    compartment = -1;  // index into Ship::compartments, or -1 if none is wet
    double planArea = 0;      // m^2, of the compartment floor
    double depth = 0;         // m, mean over the geometric region V/mu
    double breadth = 0;       // m, plan area / length: the mean breadth, not the maximum
    double pocketingRad = 0;  // atan(2h/b): above this the surface has left the high side

    // The fraction of the linear theory's free-surface moment that survives at
    // `heel`. Below the pocketing angle the surface spans the compartment and this
    // is 1; above it the water is a triangular prism of the same area whose base is
    // `b/sqrt(tau)` wide, and the moment collapses as `3/tau - 2/tau^1.5` in
    // `tau = tan(heel) b / 2h`. Exact for a box, and about 1% on the ferry's own
    // tapered vehicle deck.
    double momentFractionAt(double heelRad) const;
};

class Ship;

// The wettest free surface aboard, by floor area -- which on any ro-pax is the
// vehicle deck by an order of magnitude, and is the surface that decides her.
// `compartment` is -1 when nothing is wet, on the same test `Ship::diagnostics()`
// uses to decide whether to refine the sampling angle at all.
FreeSurfaceLayer largestFreeSurface(const Ship& s);

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

    // --- Water promotion ---
    // Decided which flooded compartments deserve flip::Solver's resolved flow instead
    // of the lumped waterVolume model. Follows the same pattern as structural/gas
    // promotion: the promoter decides, Ship owns the fields and steps them.
    promotion::WaterPromoter waterPromoter_;
    std::map<int, flip::Field*> activeWaterFields_;  // compartment index -> FLIP field

    // --- Lifecycle ---

    Ship() = default;
    ~Ship();

    // Copy/move: activeWaterFields_ contains raw pointers that must be deep-copied.
    // Move leaves the source empty (nullptr fields cleared).
    Ship(const Ship& other);
    Ship& operator=(const Ship& other);
    Ship(Ship&& other) noexcept;
    Ship& operator=(Ship&& other) noexcept;

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

    // Heat leaking from the gas into the structure around it, once per step.
    // Separate from the free-surface update because that runs twice a step and a
    // relaxation applied twice would relax at twice the rate it was asked for.
    void relaxGasToStructure(double dt);

    // Pressure and phase presented by one side of an opening.
    struct SideState {
        double pressure = kPatm;
        bool   isWater = false;

        // The density of that fluid, twice over, because the two uses of a density
        // here are different questions.
        //
        // `flowDensity` is the real local density `p / (R T)`: what Torricelli
        // divides by and what turns a volumetric rate into a mass one.
        //
        // `buoyancyDensity` is only ever read as a *difference across the opening*,
        // and there it is deliberately `kPatm / (R T)` -- the same choice
        // `gasBuoyancyHead` makes and for the same reason. Against a reference
        // atmosphere the difference between two gas spaces at ambient temperature
        // is exactly +0.0 whatever their pressures, so a cold ship gets no buoyant
        // exchange at all rather than a small one that would have to be argued
        // about. Two compartments at different *pressures* stratifying through a
        // hatch is a real effect and a different one; see §1 of docs/02-simulation.
        double flowDensity = kPatm / (kRAir * kTAmbient);
        double buoyancyDensity = kPatm / (kRAir * kTAmbient);
    };
    SideState sideStateAt(int compartmentIndex, const Vec3& worldPos, const Sea& sea) const;

    // Which end of a horizontal opening is the upper one, and how much of its area
    // still faces up. Invalid for every opening that is not a horizontal one.
    struct HorizontalSides {
        bool   valid = false;
        int    upper = kSea;
        int    lower = kSea;
        double area = 0;           // m^2, the opening projected on the horizontal
    };
    HorizontalSides horizontalSidesOf(const Opening& o, const Mat3& R) const;

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
