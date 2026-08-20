// SPDX-License-Identifier: MIT
//
// A two-zone compartment fire that pushes hot gas through the ship's own opening
// network -- `docs/06-roadmap.md` Phase 4's "multi-zone compartment fire, species
// transport through the opening network".
//
// This is very nearly *only* the gas. There is no combustion chemistry here (the
// heat release rate is an input, a design fire), no soot particle model and no
// radiation view factors. What this file answers is: given a fire of a known power
// in a known space, how hot does the gas get, how far down does the smoke layer
// come, and what crosses the doors and vents.
//
// **The one thing it does reach outside itself for is the boundary**, and that is
// the Phase 4 milestone: `wallExchange` at the bottom of this file turns a
// two-zone gas into the `thermal::Film`s a conduction solve wants, and takes the
// steel's own surface temperature back. It is a handful of lines and it needed no
// new mechanism on this side at all -- see the note there. Strength stays out:
// what a hot bulkhead then carries is `thermal::HeatedMember`'s business and this
// file does not know what a bulkhead is.
//
// --- Why this belongs next to the flooding network -----------------------------
//
// `ship.hpp` already carries every hole, door, hatch, vent and pipe as an orifice
// in one network, and that network is already two-phase: `Ship::solveFlowNetwork`
// moves water *or* air through each opening depending on what is sitting against
// it on the high-pressure side. A fire is the same network with buoyancy, a heat
// release and a hot upper layer. So this file does not build a second network --
// it reads `Ship::openings` and works out what the gas does through the openings
// it is given.
//
// The coupling is deliberately **one way by default**: `step()` takes the ship by
// const reference and writes nothing. See `Model::applyTo` for the write path and
// for why it is opt-in.
//
// --- The model, and what it assumes --------------------------------------------
//
// Each tracked compartment holds two well-mixed layers: a hot upper one under the
// deckhead and a cool lower one on the floor, separated by a horizontal interface
// that descends as the fire's plume pumps air out of the lower layer into the
// upper. That is the standard two-zone description (Quintiere; CFAST; Karlsson &
// Quintiere ch. 6) and it is the right first answer -- a field model is a later
// roadmap item and is explicitly named as such.
//
// State per layer is **mass and internal energy**, not mass and temperature. That
// choice is what makes the closure exact rather than iterative:
//
//   * `T_k = U_k / (m_k c_v)`
//   * the two layers share one pressure, so `p = (gamma - 1) (U_u + U_l) / V`
//   * and therefore `V_u / V = U_u / (U_u + U_l)` **exactly** -- the volume split
//     is the internal-energy split, with no iteration and no root find. That
//     identity is asserted directly in `tests/test_fire.cpp`; it is the single
//     load-bearing algebraic fact in this file.
//
// The layers exchange volume at the common pressure, so each does `p dV_k/dt` of
// work on the other. Substituting `V_u = V U_u / U` into the open-system energy
// equation of each layer closes the system in one line (see `layerSplit()` in
// fire.cpp): with `E_k` the non-work energy input to layer k and `E = E_u + E_l`,
//
//     dU_u/dt = [E_u + (gamma - 1) f E] / gamma,   f = U_u / U
//
// and the two rates sum to `E`, so total energy is conserved to machine precision
// whatever the split does. That is the property the account below checks.
//
// **Compressible, not isobaric.** Many zone models fix `p = p_atm` and let the
// expansion leave through the vent. That is fine for a room with a door and
// wrong for a sealed space, and this ship has plenty of nearly sealed spaces --
// a wing tank vents through a 0.02 m^2 air pipe. Carrying the pressure makes the
// sealed case exact: with no vent and no wall loss, `p(t) = p_0 + (gamma-1) E/V`
// to machine precision, which `tests/test_fire.cpp` asserts as a closed form.
//
// **Temperature is in KELVIN**, everywhere, on the same terms as `thermal.hpp`:
// the next thing that reads a gas temperature is a radiation term in `T^4`, which
// is meaningless in Celsius. `thermal::kCelsius` is the one conversion point.
//
// --- The vent flow, and why a single pressure difference is not enough ----------
//
// `Ship::solveFlowNetwork` takes **one** pressure difference, at the orifice
// centre, and that is right for a small hole. It is not right for a doorway with
// a hot layer on one side, and the failure is not a matter of accuracy:
//
//   A steady fire in a room with an open door has hot gas leaving through the top
//   of the doorway and cool air entering through the bottom, in nearly equal mass
//   flows. The *net* mass flow is close to zero while the *enthalpy* flow is
//   megawatts. A single-pressure-difference treatment sees the near-zero net,
//   moves almost nothing, and the compartment then heats without bound because
//   nothing carries the energy out. The two-layer vent integral is not a
//   refinement of that model; it is the difference between having ventilation and
//   not having it.
//
// So `ventMassFlow()` integrates `Cd sqrt(2 rho |dp|)` over the vent's height,
// splitting the span at each layer interface and at the neutral plane, where
// `dp` changes sign. Within each resulting band both densities are constant, so
// `dp(z)` is linear and the integral is **closed form** -- `(2/3|b|) |u^{3/2}|`
// between the ends -- with no quadrature and no convergence question. The
// square root has an infinite derivative at the neutral plane, which is exactly
// why splitting there rather than quadraturing through it matters.
//
// **`Opening` carries an area and a centre but no height**, and a height is what
// the integral needs. `ventShapeFor()` derives one from the opening's kind: a
// door is 2 m tall and as wide as its area requires, everything else is taken as
// square. Every field of the result is overridable per opening, because a derived
// aspect ratio is a guess and an authored one is not. The ferry's openings span
// 0.02 m^2 air pipes to a 3.6 m^2 watertight door, and
// `tests/test_fire.cpp` measures what the neutral plane is worth across that
// range rather than asserting that it matters.
//
// **Horizontal openings are the acknowledged gap.** A hatch in a deck is not a
// neutral-plane problem: buoyancy drives an unstable exchange across it
// (Epstein/Cooper) that has no pressure difference at all. Those are treated here
// as pressure-driven only, so a fire under a closed-to-flow-but-open hatch will
// not shed smoke upward the way it should. Named rather than hidden.
//
// --- The second species, and why it goes *down* ---------------------------------
//
// A total-flooding agent -- CO2, or one of the inert-gas blends -- is a second gas
// with its own molar mass, its own `R` and its own `gamma`, and every closure above
// generalises to carry it **exactly**:
//
//     p V_k = (m_air R_air + m_agent R_agent) T_k          (Dalton, per layer)
//     T_k   = U_k / (m_air c_v,air + m_agent c_v,agent)
//   =>  p V_k = (gamma_k - 1) U_k  with  gamma_k - 1 = R_mix,k / c_v,mix,k
//   =>  p     = [ (gamma_u-1) U_u + (gamma_l-1) U_l ] / V
//   and V_u/V = (gamma_u-1) U_u / [ (gamma_u-1) U_u + (gamma_l-1) U_l ].
//
// So the volume split is still the *pressure-energy* split, still with no root
// find. `Layer::excessEnergy` carries the whole of the difference: it is the extra
// internal energy a pure-air layer would need to press as hard as this one does,
// it is proportional to the agent mass in the layer and therefore **exactly 0.0**
// when there is none, and every formula above is then bit-for-bit the one that was
// here before agents existed. That is asserted, not hoped for.
//
// **The interesting part is that the agent stratifies the wrong way up.** CO2 is
// M = 44 against air's 29, so it is half again as dense and it sinks -- to the
// deck, which is where the people are, and out of any opening near the deck, which
// is where the concentration goes. Nothing about the hot-layer machinery gives that
// for free: two well-mixed zones with the agent dropped into them would report one
// uniform concentration and would have thrown away the entire hazard.
//
// The mechanism is a **buoyant separation flux** across the interface
// (`GasCompartment::settlingVelocity`), and it is two one-way streams rather than a
// diffusive exchange, because a diffusive exchange drives towards *uniformity* and
// what gravity actually does is drive towards segregation:
//
//   * the agent in the layer it does not belong in falls (or rises) into the other
//     one, at a rate `agent_there / (thickness / w)`;
//   * the carrier gas in the layer the agent does belong in is floated out, at a
//     rate `carrier_there / (thickness / w)` times the agent's own volume fraction
//     there -- so a compartment with no agent in it has **both** streams at exactly
//     zero and nothing changes.
//
// Both are applied as exact relaxations (`expm1`, as the wall loss and the drencher
// are), so a layer of half a millimetre cannot be over-drained at any step. The
// fixed point of the pair is *all* the agent in one layer and *all* the carrier in
// the other, and that state is exactly the fully stratified one: the blanket's
// height is then the agent's own partial volume, `H y`, which
// `tests/test_fire.cpp` asserts as a closed form. The other limit, `w = 0`, is the
// perfectly mixed one, and it is asserted against the ideal-gas arithmetic. The
// value of `w` between them is a modelling choice and is documented as one.
//
// **Extinguishing is displacement, not cooling.** The air the fire breathes is the
// air its plume entrains, so the two layers are weighted by `Plume::entrainment`
// itself -- the correlation already here, and not a new coefficient -- and the heat
// release the gas receives is the design curve times a linear availability ramp
// between the fuel's limiting oxygen concentration and ambient. Weighting by the
// plume also makes the oxygen a *continuous* function of the interface height,
// which matters because an agent blanket puts that interface across the fire's own
// base. The design curve is still the input, exactly as `Drencher` leaves it --
// what changed is that the model now reports how much of it the atmosphere allows,
// on the same terms as `StepResult::suppressionCooling` reports applied rather than
// demanded cooling.
//
// **Named rather than hidden**, in this half as in the others: oxygen is *diluted*
// here and never *consumed*, so a sealed compartment fire that would vitiate itself
// out does not; a compartment holding a hot smoke layer, cool air and an agent
// blanket has three strata and two zones to put them in, so the blanket and the
// cool air share one; the agent's own radiative absorption is not modelled; and
// `les::demote` writes a `GasCompartment`'s layers without touching their agent
// mass, so a promoted-and-demoted compartment loses the agent it was holding.
//
// Body frame and SI units per CLAUDE.md. The interface is horizontal in the
// *body* frame, which is exact for an upright ship and approximate for a heeled
// one; a heeled compartment's layer geometry is a real problem and is not solved
// here.
#pragma once

#include "ship.hpp"
#include "thermal.hpp"   // Film, BoundaryFace: the boundary the structure sees

#include <string>
#include <vector>

namespace sim::fire {

// Air's caloric constants, derived from `kRAir` and `kGammaAir` rather than
// quoted, so that `c_p - c_v = R` and `c_p / c_v = gamma` hold *exactly*. The
// pressure closure `p = (gamma-1) U / V` is only true if they do, and a quoted
// 1005 J/(kg K) against this repo's R and gamma would break it in the fourth
// digit -- enough to show up in an energy account asserted at machine precision.
inline constexpr double kCvAir = kRAir / (kGammaAir - 1.0);   // 717.63 J/(kg K)
inline constexpr double kCpAir = kCvAir + kRAir;              // 1004.68 J/(kg K)

// Ambient air density at `kPatm` and `kTAmbient`.
inline constexpr double kRhoAmbient = kPatm / (kRAir * kTAmbient);   // 1.2253 kg/m^3

// ---------------------------------------------------------------------------
// The flooding agent
// ---------------------------------------------------------------------------

// The molar gas constant, `N_A k_B`, which the 2019 SI defines **exactly**: both
// factors are exact defining constants, so this is not a measurement and carries
// no uncertainty. `tests/test_fire.cpp` re-derives it from them rather than
// quoting it, on the same terms as `kStefanBoltzmann`.
inline constexpr double kUniversalGas = 8.31446261815324;   // J/(mol K)

// Air's molar mass **as this repo's own `kRAir` implies it**, 28.965 g/mol.
// Derived rather than quoted: every mole fraction below is a ratio of `R`s, and a
// quoted 28.9647 against a `kRAir` that means 28.9650 would put a part in ten
// thousand of drift into a concentration this file asserts as a closed form.
inline constexpr double kMolarMassAir = kUniversalGas / kRAir;   // kg/mol

// Oxygen in dry air, by volume. The one number every asphyxiation threshold below
// is measured against.
inline constexpr double kOxygenFractionAir = 0.2095;

// What a diluted atmosphere does to a person, by volume of oxygen. The standard
// occupational-hygiene ladder: 19.5% is the minimum for entry, below about 10%
// consciousness goes, below 6% death follows in minutes.
inline constexpr double kOxygenEntryLimit    = 0.195;
inline constexpr double kOxygenIncapacitating = 0.10;
inline constexpr double kOxygenLethal         = 0.06;

// One flooding agent, as a gas.
//
// **Trivially copyable and stringless on purpose.** A `VentSide` carries one of
// these and `substep` builds a `VentSide` per compartment per pressure-solve
// evaluation -- tens of times a step -- so a `std::string` in here would be an
// allocation on the hot path. The name lives on `AgentSystem`, which is authored
// once.
struct AgentSpecies {
    double molarMass = 0.0440095;   // kg/mol; CO2 by default
    double gamma = 1.289;           // CO2 at 300 K

    // Sublimation enthalpy at one atmosphere, J/kg, and the temperature it happens
    // at. Zero for an agent that is stored and delivered as a gas. See
    // `AgentSystem::solidFraction` for what they are for.
    double sublimationHeat = 5.71e5;   // J/kg
    double sublimationTemperature = 194.65;   // K, CO2's 1 atm sublimation point

    // Toxicity by volume fraction, which is the whole difference between a CO2
    // system and an inert-gas one and is the reason this is a species property
    // rather than a constant. CO2 is not merely an asphyxiant: 7-10% brings
    // unconsciousness in minutes and above 10% it convulses and kills, at
    // concentrations where the *oxygen* it leaves would still be breathable. A
    // true inert gas has no such action at all, so both are 1.0 for one -- a
    // fraction no atmosphere can exceed.
    double incapacitatingFraction = 0.07;
    double lethalFraction = 0.10;

    double gasConstant() const { return kUniversalGas / molarMass; }
    double cv() const { return gasConstant() / (gamma - 1.0); }
    // `c_p - c_v = R` and `c_p / c_v = gamma` hold exactly, for the same reason
    // `kCpAir` is derived rather than quoted.
    double cp() const { return cv() + gasConstant(); }
    // Denser than air at the same pressure and temperature, so it settles
    // downward. `R = R_u / M`, so a *smaller* R is a *heavier* gas.
    bool heavierThanAir() const { return gasConstant() < kRAir; }
};

// The agents the marine rules name. Molar masses are the IUPAC 2021 values; the
// blends are the mole-weighted mean of their components.
//
//   * **CO2** -- SOLAS FSS Code ch. 5. Heavy, and lethal at its own design
//     concentration.
//   * **IG-01 (argon)**, **IG-100 (nitrogen)**, **IG-55 (50/50 N2/Ar)** and
//     **IG-541 (52% N2, 40% Ar, 8% CO2)** -- ISO 14520. Carried because nitrogen
//     is *lighter* than air and argon is heavier, so the direction a flooding
//     agent stratifies is a consequence of its molar mass in this model rather
//     than a hard-coded downward drift. `tests/test_fire.cpp` checks that by
//     changing the species and nothing else.
inline constexpr AgentSpecies kCarbonDioxide{};
inline constexpr AgentSpecies kNitrogen{.molarMass = 0.0280134,
                                        .gamma = 1.4,
                                        .sublimationHeat = 0.0,
                                        .sublimationTemperature = 77.36,
                                        .incapacitatingFraction = 1.0,
                                        .lethalFraction = 1.0};
inline constexpr AgentSpecies kArgon{.molarMass = 0.039948,
                                     .gamma = 1.667,
                                     .sublimationHeat = 0.0,
                                     .sublimationTemperature = 87.30,
                                     .incapacitatingFraction = 1.0,
                                     .lethalFraction = 1.0};
// **The blends' `gamma` is derived from their components, not quoted.**
//
// Both of these carried a hard-coded gamma -- 1.52 and 1.51 -- two lines under a
// comment saying the blends are the mole-weighted mean of their components, and
// beside a `molarMass` that really is one. For an ideal-gas mixture the *molar*
// c_v is what mole-weights, so gamma follows from
//
//     1 / (gamma_mix - 1) = sum_i x_i / (gamma_i - 1)
//
// which gives 1.5001 for IG-55 and 1.4594 for IG-541 against the 1.52 and 1.51
// that were written here. The second is the one that matters: `cv()` is
// `R / (gamma - 1)`, so a gamma 3.5% high is a c_v **9.9% low**, and c_v is what
// `Layer::heatCapacity` and `Layer::excessEnergy` are built on -- the latter
// feeding `GasCompartment::pressure()` and the `layerSplit` volume closure.
//
// The self-contained argument needs no reference book: adding 8% CO2, whose
// gamma of 1.289 is the lowest in this file, must pull the blend's gamma *down*
// from IG-55's. The old pair moved it down by 0.01; the components say 0.041.
//
// Nothing caught it because the only functional test of either blend discharges
// them isothermally, where the wrong c_v appears in the injected energy and in
// the receiving heat capacity and cancels exactly.
inline constexpr double mixtureGamma(double xa, double ga, double xb, double gb,
                                     double xc = 0.0, double gc = 1.4) {
    return 1.0 + 1.0 / (xa / (ga - 1.0) + xb / (gb - 1.0) +
                        (xc > 0.0 ? xc / (gc - 1.0) : 0.0));
}

inline constexpr AgentSpecies kIG55{.molarMass = 0.5 * 0.0280134 + 0.5 * 0.039948,
                                    .gamma = mixtureGamma(0.5, kNitrogen.gamma, 0.5, kArgon.gamma),
                                    .sublimationHeat = 0.0,
                                    .sublimationTemperature = 0.0,
                                    .incapacitatingFraction = 1.0,
                                    .lethalFraction = 1.0};
inline constexpr AgentSpecies kIG541{
    .molarMass = 0.52 * 0.0280134 + 0.40 * 0.039948 + 0.08 * 0.0440095,
    .gamma = mixtureGamma(0.52, kNitrogen.gamma, 0.40, kArgon.gamma,
                          0.08, kCarbonDioxide.gamma),
    .sublimationHeat = 0.0,
    .sublimationTemperature = 0.0,
    // 8% CO2 in the blend, so the atmosphere's CO2 reaches 8% of the agent
    // fraction. At IG-541's 40% design concentration that is 3.2% CO2 -- which is
    // the point of the blend: it is *deliberately* below the 4% IDLH while an
    // equivalent CO2 flood is at 40%.
    .incapacitatingFraction = 0.07 / 0.08,
    .lethalFraction = 0.10 / 0.08};

// t-squared growth coefficients, W/s^2. The canonical NFPA 72 / SFPE set: the
// fire reaches 1055 kW (1000 Btu/s) in 600, 300, 150 and 75 seconds.
inline constexpr double kGrowthSlow      = 1.055e6 / (600.0 * 600.0);  // 2.931
inline constexpr double kGrowthMedium    = 1.055e6 / (300.0 * 300.0);  // 11.72
inline constexpr double kGrowthFast      = 1.055e6 / (150.0 * 150.0);  // 46.89
inline constexpr double kGrowthUltrafast = 1.055e6 / (75.0 * 75.0);    // 187.6

// ---------------------------------------------------------------------------
// The plume
// ---------------------------------------------------------------------------

// Heskestad's axisymmetric fire plume, as published in the SFPE Handbook.
//
// **The correlations are written in kW and metres and return kg/s.** Every
// interface on this struct is in SI watts; the conversion happens once, inside,
// and is the single place a factor of a thousand could hide. That is not
// pedantry: this repo has already published a figure that was wrong by exactly
// that factor, in a comment nothing tested.
struct Plume {
    double heatRelease = 0;          // W, total
    double diameter = 1.0;           // m, equivalent circular fire base
    double convectiveFraction = 0.7; // the share that drives the buoyant plume

    double convectiveHeatRelease() const { return convectiveFraction * heatRelease; }

    // Mean flame height above the fire base, `L = 0.235 Q^(2/5) - 1.02 D`.
    // Negative for a wide, weak fire, which means there is no coherent flame; the
    // entrainment below then uses the far-field branch throughout.
    double flameHeight() const;

    // Virtual origin, `z0 = 0.083 Q^(2/5) - 1.02 D`, relative to the fire base.
    // Routinely negative -- a wide fire behaves like a point source *below* the
    // pan -- and that sign is load bearing.
    double virtualOrigin() const;

    // Entrained mass flow at `height` above the fire base, kg/s. Two branches:
    //
    //   * in the flame, `z <= L`:  `m = 0.0056 Qc (z/L)`
    //   * above it:                `m = 0.071 Qc^(1/3) (z-z0)^(5/3) + 0.0018 Qc`
    //
    // They are independent formulae and they agree at `z = L` to **1.78%** --
    // measured, not remembered, and `tests/test_fire.cpp` derives that ratio from
    // the published coefficients and asserts it to 1e-9. It is a constant: both
    // the flame height and the virtual origin carry the same `Q^(2/5)`, so the
    // agreement does not depend on the fire's size or its diameter at all.
    double entrainment(double height) const;
};

// ---------------------------------------------------------------------------
// The design fire
// ---------------------------------------------------------------------------

// A heat release rate as a function of time: growth, steady, decay. The boundary
// of this file's scope -- what burns and how fast is an input, not a result.
struct DesignFire {
    std::string name;
    int    compartment = 0;      // index into Model::gas
    double baseZ = 0;            // body frame; the height of the fire base
    double diameter = 1.0;       // m

    double growthCoefficient = 0;   // W/s^2; Q = alpha t^2 while growing
    double peakHeatRelease = 0;     // W; the cap the growth phase reaches
    double steadyDuration = 0;      // s at the peak, after growth, before decay
    double decayDuration = 0;       // s of linear decay to zero; 0 means never

    double convectiveFraction = 0.7;
    // The share of the heat release that leaves the gas without heating it --
    // flame radiation absorbed by surfaces and lost through the vent.
    //
    // **Defaults to zero, and that is deliberate.** The obvious value is the
    // textbook 0.3-0.35 radiative fraction, but the boundary loss in this model
    // is a lumped `h_k A (T_g - T_wall)` whose coefficient is the one MQH
    // correlates against, and MQH is fitted using the *total* heat release.
    // Subtracting a radiative fraction here as well double counts the same loss
    // and drops the steady upper-layer temperature about 25% below MQH. Set it
    // non-zero only alongside a boundary model that does not already carry it.
    double radiativeLossFraction = 0.0;

    // Combustion products generated per joule released, kg/J. The default is a
    // 0.05 kg/kg soot-and-products yield on a fuel of 20 MJ/kg heat of
    // combustion; it is a passive scalar and nothing in this file reads it back
    // into the physics, which is what makes it safe to carry at a guessed yield.
    double productYield = 0.05 / 20.0e6;

    // The oxygen volume fraction below which this fuel will not burn, and the
    // whole of how a flooding agent puts a fire out here.
    //
    // The heat release the gas receives is the design curve times a linear ramp
    // from this to ambient -- `(X_O2 - LOC) / (0.2095 - LOC)`, clamped -- taken on
    // the **lower** layer, because that is what the plume entrains. Linear rather
    // than a cliff at the limit for the reason `kEvaporationBand` in fire.cpp
    // exists: a discontinuity in the source term is a discontinuity the substep
    // controller has to subdivide against, at the one concentration a flooded
    // space actually settles near.
    //
    // 0.13 is a middling value for the flaming combustion of ordinary liquid and
    // solid fuels; NFPA 12's own minimum design concentrations imply 0.128 for
    // heptane and 0.149 for methane. Whatever it is, it is a **fuel** property
    // this file takes as an input, exactly as the heat release curve is.
    //
    // What this does *not* model, and it matters: the fire is never allowed to
    // vitiate its own air. There is no oxygen state here, only dilution by the
    // agent, so a sealed compartment fire that would go out on its own does not.
    double limitingOxygen = 0.13;

    double heatRelease(double time) const;   // W

    // The share of `heatRelease` an atmosphere at `oxygenFraction` allows, in
    // [0, 1]. Exactly 1.0 at `kOxygenFractionAir`, which is what keeps a model
    // with no agent in it bit-identical.
    double oxygenAvailability(double oxygenFraction) const;

    // The time at which the fire has released `energy` joules, or infinity if it
    // never does. Used by the tests to state a conservation account over a run
    // whose released energy is known independently of the integrator.
    double totalEnergy() const;              // J over the whole design curve
};

// ---------------------------------------------------------------------------
// Gas state
// ---------------------------------------------------------------------------

// One well-mixed layer. Mass and *internal energy*, not mass and temperature --
// see the header comment; the pressure closure depends on it.
struct Layer {
    double mass = 0;       // kg, air *and* agent together
    double energy = 0;     // J, internal (m c_v T)
    double products = 0;   // kg of combustion products carried

    // Flooding agent, kg, and it is **part of `mass`** -- unlike `products`, which
    // is a passive tracer riding on top of it. That difference is not cosmetic: a
    // machinery space at its design concentration is a *third* agent by mass, so
    // an agent kept outside the mass would leave the density, the vent flow and
    // the mass account all wrong by that much, while a soot yield of 0.05 kg/kg
    // on a 20 MJ/kg fuel is parts per million and safely a tracer.
    double agent = 0;

    // Every mixture query needs to know what the agent *is*, and a `Layer` does
    // not carry it -- `GasCompartment::agentSpecies` does.
    //
    // **The default is exact wherever it is reachable and wrong nowhere else.**
    // With `agent == 0` every expression below multiplies the species terms by
    // exactly zero, so a caller with no agent gets the identical double whatever
    // it passes; and nothing inside `fire.cpp` ever takes the default, because a
    // `GasCompartment` always has its own species to hand. The default exists for
    // the callers outside this file that read `upper.temperature()` on gas nobody
    // has flooded.
    double temperature(const AgentSpecies& s = kCarbonDioxide) const;   // K
    double productFraction() const;          // kg/kg
    double agentMassFraction() const;        // kg/kg of the layer's own mass

    // The agent's share **by volume**, which for an ideal-gas mixture is its mole
    // fraction and its partial-pressure fraction alike -- and is the unit every
    // suppression rule and every toxicity threshold is written in.
    double agentFraction(const AgentSpecies& s = kCarbonDioxide) const;

    // `sum_i m_i c_v,i`, J/K, and `sum_i m_i R_i`, J/K. The two sums every mixture
    // property here is a ratio of.
    double heatCapacity(const AgentSpecies& s = kCarbonDioxide) const;
    double gasConstantSum(const AgentSpecies& s = kCarbonDioxide) const;

    // The extra internal energy a layer of *pure air* would need to exert the same
    // pressure in the same volume: `U [(R_mix/c_v,mix) / (gamma_air - 1) - 1]`.
    //
    // This is the one number the whole second species costs the closure. With it,
    // `p = (gamma_air - 1)(U + x_u + x_l)/V` and `V_u/V = (U_u + x_u)/(U + x_u + x_l)`
    // -- the same two lines as before, with one term added. Written out it is
    //
    //     x = U * agent * (c_v,air R_agent - c_v,agent R_air)
    //             / (c_v,mix_sum * R_air)
    //
    // whose numerator carries `agent` as a *factor*, so it is **exactly 0.0** for a
    // layer with no agent in it -- not small, zero -- and `y + 0.0 == y`. That is
    // what keeps every pre-existing figure bit-identical, and `tests/test_fire.cpp`
    // asserts it as bit-identity rather than as closeness.
    double excessEnergy(const AgentSpecies& s = kCarbonDioxide) const;
};

// The gas in one compartment: two layers, a box of known height and footprint,
// and a lumped boundary conductance.
struct GasCompartment {
    int    shipCompartment = kSea;   // index into Ship::compartments, kSea if standalone
    std::string name;

    // Geometry. `attach()` fills these from the ship's own compartment mesh; a
    // test may set them exactly. The space is treated as a prism of `floorArea`
    // between `floorZ` and `ceilingZ`, which is what a two-zone model needs and
    // is what makes the interface height a single number.
    double floorZ = 0, ceilingZ = 0;   // body frame, m
    double floorArea = 0;              // m^2
    double perimeter = 0;              // m, for the wall area either side of the interface
    double gasVolume = 0;              // m^3 available to gas

    // Lumped boundary conductance, W/(m^2 K), applied over the area each layer
    // actually wets. This is MQH's `h_k`: `sqrt(k rho c / t)` while the boundary
    // is thermally thick, `k / delta` once the heat has been through it.
    //
    // Applying it over the *wetted* area rather than the whole enclosure is the
    // point of having two zones at all, and it matters: charging the full
    // enclosure area at the upper-layer temperature, which is what MQH's own
    // `A_T` does, under-predicts the steady temperature by about a quarter.
    double wallConductance = 30.0;
    double wallTemperature = kTAmbient;   // K

    // What has been, or could be, discharged in here. One species per compartment:
    // two different agents mixed in one space would need a mixture rule on every
    // layer, and `validate()` refuses a vent between compartments carrying
    // different ones rather than silently averaging them.
    AgentSpecies agentSpecies{};

    // How fast the agent separates from the carrier gas by buoyancy, m/s.
    //
    // **This is the coefficient nothing here measures**, and it is exactly the
    // `Drencher::evaporatedFraction` situation: the two limits are both closed
    // forms this file asserts, and the value between them is a choice. At zero the
    // compartment is perfectly mixed and the concentration is the ideal-gas
    // arithmetic everywhere; in the limit the compartment is perfectly stratified,
    // the blanket is pure agent and its depth is the agent's own partial volume.
    // Both extremes are real: a discharge in progress is violently stirred by its
    // own jets -- NFPA 12 requires two minutes for a surface fire *because* rapid
    // discharge mixes -- and a space left alone afterwards settles.
    //
    // The default gives a time constant of `thickness / w`: about five minutes to
    // separate a 6 m machinery space, which is the order the holding-time
    // requirements are written at. It is a guess and it is labelled one;
    // `tests/test_fire.cpp` measures what it is worth rather than asserting it.
    double settlingVelocity = 0.02;   // m/s

    Layer upper, lower;

    double pressure() const;        // Pa, absolute, at the floor
    double gaugeAtFloor() const;    // Pa above the ambient hydrostatic profile
    double upperVolume() const;     // m^3
    double lowerVolume() const;     // m^3
    double interfaceZ() const;      // body frame, m
    double totalEnergy() const { return upper.energy + lower.energy; }
    double totalMass() const { return upper.mass + lower.mass; }
    double totalAgent() const { return upper.agent + lower.agent; }
    // The agent's share of the whole space by volume, which is the number a
    // suppression rule states a design concentration in.
    double agentFraction() const;

    // Fill both layers with still ambient air. The upper layer gets a token
    // `seedFraction` of the volume so that its temperature is defined before the
    // fire has made one; with both layers at ambient the interface then sits that
    // fraction below the deckhead and the gauge pressure is exactly zero.
    void fillAmbient(double seedFraction = 1e-4);
};

// ---------------------------------------------------------------------------
// Vents
// ---------------------------------------------------------------------------

// The height and width the vent integral needs, which `Opening` does not carry.
struct VentShape {
    double height = 0;      // m, sill to soffit
    double width = 0;       // m
    bool   horizontal = false;
};

// A default shape for an opening, from its kind and area:
//
//   * `Door`   -- 2.0 m tall, or square if its area cannot make a 2 m door.
//   * `Hatch`  -- horizontal: a hatch is in a deck.
//   * anything else -- square. A breach is a tear of no particular shape and an
//     air pipe is a round hole; neither has an aspect ratio worth inventing.
//
// A guess, and labelled as one. `Model::attach` uses it and `Model::vents` is
// public so a caller with real drawings can overwrite the result.
VentShape ventShapeFor(const Opening& o);

// One opening as the gas model sees it.
struct Vent {
    int    opening = -1;             // index into Ship::openings, -1 if synthetic
    std::string name;
    int    a = kSea, b = kSea;       // indices into Model::gas, or kSea for outside air
    double sillZ = 0, soffitZ = 0;   // body frame
    double width = 0;                // m
    double area = 0;                 // m^2, used for horizontal vents
    double dischargeCoeff = 0.6;
    bool   horizontal = false;
    bool   open = true;

    // Diagnostic, filled every step. Both are non-negative: a vent spanning the
    // interface runs in *both* directions at once, and reporting only the net
    // would throw away the entire mechanism.
    double massAToB = 0;             // kg/s
    double massBToA = 0;             // kg/s
    double neutralPlaneZ = 0;        // body frame; valid only when `bidirectional`
    bool   bidirectional = false;
    bool   blockedByWater = false;   // water against it: the flooding solve owns it
    // The span the last step actually integrated over, after the vent was slid
    // and trimmed into the gas space on both sides. Published because that
    // adjustment is the difference between the ferry's air escapes existing and
    // not -- they are authored 5.5 m above the compartment they drain -- and a
    // silent geometric correction is one nothing can check.
    double activeSillZ = 0, activeSoffitZ = 0;

    double netMass() const { return massAToB - massBToA; }
};

// The gas column against one side of a vent, as the integral sees it: two
// constant-density bands and a reference pressure.
//
// Pressures here are **gauge against the ambient hydrostatic profile**
// `p_atm(z) = kPatm - rho_inf g z`, not absolute. That is what makes a
// compartment of still ambient air produce *exactly* zero flow rather than
// something that cancels to a rounding error, and it is what the zero-heat-release
// control in `tests/test_fire.cpp` rests on.
struct VentSide {
    double gaugeAtFloor = 0;     // Pa at `floorZ`
    double floorZ = 0;           // m
    double interfaceZ = 0;       // m; at or below `floorZ` means "upper layer only"
    double rhoLower = kRhoAmbient, rhoUpper = kRhoAmbient;   // kg/m^3
    double tLower = kTAmbient, tUpper = kTAmbient;           // K
    double yLower = 0, yUpper = 0;                           // product mass fractions
    double aLower = 0, aUpper = 0;                           // agent mass fractions
    AgentSpecies agentSpecies{};                             // what `aLower`/`aUpper` are

    double densityAt(double z) const { return z >= interfaceZ ? rhoUpper : rhoLower; }
    double temperatureAt(double z) const { return z >= interfaceZ ? tUpper : tLower; }
    double productFractionAt(double z) const { return z >= interfaceZ ? yUpper : yLower; }
    double agentFractionAt(double z) const { return z >= interfaceZ ? aUpper : aLower; }
    // The mixture's `c_p`, J/(kg K), which is what an enthalpy flux wants. Exactly
    // `kCpAir` where there is no agent.
    double heatCapacityOf(bool upper) const;
    double heatCapacityAt(double z) const { return heatCapacityOf(z >= interfaceZ); }
    // The temperature at which *pure air* would be as dense as this gas is: `T`
    // times `R_mix/R_air`. It is what decides which layer an incoming stream joins,
    // because a doorway can deliver gas that is warm and still sinks -- half again
    // as heavy is worth 150 K of buoyancy, and a rule written on temperature alone
    // would file a CO2 blanket under the deckhead. Exactly `T` without an agent, so
    // the pure-air deposition rule is unchanged to the bit.
    //
    // **Two callers and one definition**, deliberately: the vent integral asks by
    // height and the deposition rule asks by layer, and an earlier draft answered
    // the second by writing the formula out a second time inside `applyStream`.
    // That is the "two answers that agree until the day they do not" this file has
    // already deleted once, and it also left the by-height form with no caller and
    // therefore no test.
    double buoyancyTemperatureOf(bool upper) const;
    double buoyancyTemperatureAt(double z) const {
        return buoyancyTemperatureOf(z >= interfaceZ);
    }
    // Gauge pressure at `z`: the ambient profile is already subtracted, so the
    // integrand is the *buoyancy* difference and still air gives zero.
    double gaugeAt(double z) const;
};

// Outside air: gauge pressure zero at every height, by construction.
VentSide ambientSide();

// What crosses a vent, resolved into the layer it leaves and the layer it joins.
struct VentResult {
    double massAToB = 0;         // kg/s
    double massBToA = 0;         // kg/s
    double enthalpyAToB = 0;     // W
    double enthalpyBToA = 0;     // W
    // Of `massAToB`, how much left A's *upper* layer; likewise for B. The rest
    // came from the lower one. Deposition on the receiving side is by buoyancy,
    // not by which layer it left.
    //
    // The stream's species loading is deliberately **not** here. It would be the
    // donor's product fraction times the mass, which is a number the caller
    // already has both halves of, and a second copy of it inside this struct was
    // written on every band and read by nothing at all -- the whole of Phase 2's
    // "two answers that agree until the day they do not" in miniature. Mutation
    // testing found it: every mutant of those two lines survived, because nothing
    // could see them.
    double fromUpperA = 0, fromUpperB = 0;
    double neutralPlaneZ = 0;
    bool   bidirectional = false;
};

// The vent integral. Splits `[sillZ, soffitZ]` at each side's layer interface and
// at every neutral plane, then integrates each band in closed form. Positive
// `massAToB` means gas leaving side A.
VentResult ventMassFlow(const Vent& v, const VentSide& a, const VentSide& b);

// ---------------------------------------------------------------------------
// Suppression
// ---------------------------------------------------------------------------
//
// `docs/06-roadmap.md` Phase 4's "suppression systems, **and their effect on
// stability**", and the second half of that phrase is the whole item. Water put
// into a compartment to fight a fire is water in a compartment. Firefighting has
// capsized ships that the fire would not have.
//
// **Suppression water is floodwater, and it needed no new path.** `Ship` carries
// floodwater as real mass at the real centroid of the real water body and
// re-levels it against gravity every tick; `rightingArmAtHeel` re-levels it again
// at each heel it is asked about. There is no free-surface *correction term*
// anywhere in `ship.cpp` -- no `rho i / disp`, no tabulated allowance -- because
// the effect is what the re-levelling does. So a cubic metre added to
// `Compartment::waterVolume` arrives with its mass, its centroid, its free
// surface and its cost to GM already correct, and `applyTo` writes exactly that.
// `tests/test_fire.cpp` measures the emergent effect against the closed form it
// is supposed to reproduce and finds `GM_solid - GM_liquid = rho mu I / disp` to
// **1.4e-6** relative on a box barge, where `I` is `b^3 l / 12` and nothing is
// approximated, and to 2e-4 on the ferry's own vehicle deck, where `I` is itself
// a numerical integration over the compartment's mesh. The `mu` is not optional:
// `rightingArmAtHeel` levels the *geometric* region `waterVolume / mu` while the
// mass is `waterVolume * rho`, so on the ferry's 0.90 deck leaving it out
// over-states the loss by 11%.

// Water's caloric constants at one atmosphere. Quoted, not derived: unlike the
// air constants above, nothing here has to satisfy an identity.
inline constexpr double kCpWater      = 4186.0;    // J/(kg K), liquid
inline constexpr double kLatentHeat   = 2.257e6;   // J/kg of vaporisation at 373.15 K
inline constexpr double kTSaturation  = 373.15;    // K, boiling at kPatm

// A spray density in litres per minute per square metre -- the unit every
// suppression rule is written in, SOLAS II-2 Reg. 20's 5 L/(min m^2) for a
// ro-ro space among them -- as a mass flow in kg/s over `area`.
//
// One line, and it exists for the reason `Plume`'s kilowatts do: this is the
// single place a factor of sixty or a factor of a thousand could hide, and this
// repo has already published a figure that was wrong by exactly such a factor in
// a comment that nothing tested.
double sprayMassFlow(double area, double litresPerMinutePerSquareMetre);

// A drencher, sprinkler or watermist head set: a mass flow of water into one
// compartment, part of which evaporates in the hot layer and part of which does
// not.
//
// **The split is the entire mechanism.** What evaporates is where the cooling
// comes from -- 2.257 MJ/kg, so a few kg/s is megawatts, comparable with the
// fire itself. What does not evaporate lands on the deck, and on a ship that is
// where the stability cost comes from. The same kilogram cannot do both.
struct Drencher {
    std::string name;
    int    gasCompartment = 0;         // index into Model::gas
    double flow = 0;                   // kg/s of water delivered
    double waterTemperature = kTAmbient;   // K, as it leaves the nozzle
    bool   on = false;

    // The share of the *engaged* water that evaporates in the hot layer.
    //
    // **This is the coefficient nothing here measures**, and it is the one the
    // answer turns on. The literature spreads it over most of the unit interval:
    // a fine mist in a hot post-flashover layer converts nearly all of it, while
    // a ro-ro drencher throwing coarse drops through a 500 K layer onto vehicles
    // converts very little and puts the rest on the deck. 0.3 is a deliberately
    // middling default and it is a guess.
    //
    // What `tests/test_fire.cpp` establishes instead of asserting a value is the
    // *sensitivity*, and it has two halves that point opposite ways.
    //
    // Per kilogram, cooling runs 582 kJ at 0.1 to 2387 kJ at 0.9: a factor of
    // **4.1**, not the factor of nine the latent heats alone suggest, because the
    // sensible heat of raising the water to saturation is a floor no fraction can
    // go below and is a third of the total at 0.3. The water landing runs the
    // other way and by *more* -- the deck collects `1 - e`, so 0.1 and 0.9 differ
    // by **9x**. On the arithmetic alone the coefficient matters more to the
    // stability answer than to the cooling one.
    //
    // **And then on a real ship it barely matters at all**, which is the finding
    // worth having. This fraction is a ceiling and not a rate: what actually
    // evaporates is bounded by the power available to boil it. The ferry's 20 MW
    // lorry fire against a 31 kg/s drencher converts **3.3%** of the flow over
    // four hours, not 30%, because most of the fire's energy goes into raising
    // 447 tonnes of water to saturation and out through the vents. So the water
    // on the deck -- the whole stability half of this item -- is set by the fire's
    // power and the freeing ports rather than by this number: on the ISO-room
    // fixture, where the same over-supply holds, moving it from 0.1 to 0.9 moves
    // the water landing by **2%** against the 9x the arithmetic alone predicts.
    double evaporatedFraction = 0.3;

    // Diagnostic, filled every step.
    double lastCooling = 0;      // W taken out of the upper layer
    double lastEvaporated = 0;   // kg/s that left as steam
    double lastToDeck = 0;       // kg/s that landed
};

// A freeing port, scupper or drain: the thing that is supposed to get the water
// back out, and the reason a drencher does not simply fill the ship.
//
// Modelled as a **rectangular notch of width `width` with its sill at `sillZ` and
// no top**, which is what a freeing port in a bulwark is and, more to the point,
// is what can be written down without inventing a height. `Opening` carries an
// area and a centre and no height at all; deriving one is what `ventShapeFor`
// has to do and it is a guess there. A notch with no top needs no such guess,
// because a layer a drencher makes is centimetres deep and never reaches a
// soffit.
//
// The discharge is one integral of Torricelli's velocity over the wetted band,
// split at the outside water level exactly as `ventMassFlow` splits at a layer
// interface, and closed form in each band -- see `scupperFlow`.
struct Scupper {
    std::string name;
    // Indexed into `Model::gas` rather than into `Ship::compartments`, on the same
    // terms as `DesignFire::compartment`, so that every entry in this file speaks
    // one index space. The cost is that a deck with freeing ports has to be a
    // tracked gas space even if nothing is burning in it; `validate()` says so
    // rather than leaving a port that silently never runs, which is
    // indistinguishable from a blocked one.
    int    gasCompartment = 0;

    // The middle of the sill, in the body frame, exactly as `Opening::pos` is the
    // orifice centre.
    //
    // **A position and not a height, and that was found by measuring rather than
    // by thinking.** The first version carried a `sillZ` and compared it against
    // `Compartment::surfaceOffset`, which reads like a height and is not one: it
    // is the offset of the water's plane along the *body* up-vector, and the
    // moment this ship takes up an angle of loll the two differ by the beam times
    // the sine of it. The ferry lolls within a minute of the drencher starting,
    // so the ports silently stopped draining at the exact moment they were needed
    // and the run looked like a suppression system with no freeing ports at all.
    //
    // With a position the head is `surfaceWorldZ - sillWorldZ`, which is the same
    // comparison `Ship::sideStateAt` makes and is exact at any attitude. It also
    // buys the physics the height could not express: on a lolled ship the low-side
    // ports run hard and the high-side ones stand dry.
    Vec3   sillPos{};
    double width = 0;               // m of clear opening
    double dischargeCoeff = 0.6;
    // Blocked by a vehicle parked against it, by debris, by a lashing, or by
    // nobody having opened it. The case that turns a survivable drencher run into
    // a capsize, so it is a field rather than a caller deleting the entry.
    bool   blocked = false;

    double lastFlow = 0;            // m^3/s, positive means overboard
    // The two heads the integral was taken between, m above the sill. Published
    // because "the port did not run" has two quite different causes -- no water
    // inside, or the sea standing over it -- and they want different actions.
    double lastInsideHead = 0, lastOutsideHead = 0;
};

// Volumetric discharge through a notch of unit `Cd * width`, m^3/s per metre of
// width, given the inside and outside water heads above the sill. Both heads are
// clamped at zero: water below the sill cannot reach it.
//
// Positive means outward. The integral is Torricelli's `sqrt(2 g d)` over the
// wetted height, split at the outside level:
//
//   * below it the notch is drowned and the driving head is the constant
//     difference `a - c`, contributing `c sqrt(2 g (a - c))`;
//   * above it the discharge is free and the head runs linearly to zero at the
//     inside surface, contributing `(2/3) (a - c)^(3/2) sqrt(2 g)`.
//
// With a dry outside this is the classical free weir `(2/3) Cd b sqrt(2g) h^(3/2)`
// exactly; with the two levels equal it is exactly zero; and with the outside
// higher it runs backwards and the port admits the sea, which is what a freeing
// port does once the ship has settled far enough to put it under. One formula for
// all three, so there is no regime boundary to get wrong.
double scupperFlow(double insideHead, double outsideHead);

// ---------------------------------------------------------------------------
// Suppression by gas
// ---------------------------------------------------------------------------
//
// `docs/06-roadmap.md` Phase 4's "CO2 and inert-gas total flooding (oxygen
// displacement, and the sealing requirement that makes it fail with a door
// open)".
//
// **Two of those three clauses are the easy half.** Getting a design
// concentration into a space is arithmetic -- `agentMassForFraction` below is the
// closed form and a bank of bottles either holds that much or does not. Holding
// it is the problem, and it is the problem *because* of the first clause: the
// agent is a mass, hundreds of kilograms of it, and putting mass into a sealed
// steel box raises its pressure, and a raised pressure drives gas out through
// every opening the box has. The concentration is then lost through the ship's
// own opening network, at a rate that depends on **where** the openings are,
// because a heavy agent lies on the deck and a low door drains the blanket while a
// high one vents air. A single pressure difference at an orifice centre cannot
// tell those two apart. The vent integral can, and that is the second time in this
// file that the height resolution turns out to be the difference between a
// mechanism existing and not.

// One bank of bottles and its nozzles.
struct AgentSystem {
    std::string name;
    int    gasCompartment = 0;        // index into Model::gas
    double flow = 0;                  // kg/s at the nozzles while discharging
    // What is in the bottles, kg. **Finite by default in the sense that matters**:
    // zero means an unlimited supply, which is a test fixture and not a ship. A
    // real bank runs out, and running out is what turns "we reached the design
    // concentration" into "we reached it and then watched it leak away".
    double charge = 0;
    bool   on = false;

    // Where the discharge is delivered, between the two layers.
    //
    // **In proportion to their volume, and that is deliberate rather than lazy.** A
    // total-flooding discharge is a high-momentum jet: NFPA 12 requires it inside
    // two minutes for a surface fire *because* speed is what mixes the space. So
    // the agent arrives mixed, and every bit of stratification this model produces
    // is then the work of `GasCompartment::settlingVelocity` and can be switched
    // off by setting it to zero -- which is what makes the mixed limit a control
    // rather than an assumption. Delivering into the layer the agent "belongs" in
    // would have hidden the separation inside the source term, and would have
    // driven 30 kg/s into a seed layer of 0.12 kg on the first substep of every
    // run.

    // The state the agent arrives in, and what the compartment pays for it.
    //
    // **CO2 discharge is a phase change and it is modelled, not waved away.** The
    // agent leaves the bottle as a liquid at storage pressure and flashes across
    // the nozzle; at one atmosphere CO2 has no liquid phase at all, so what comes
    // out is cold vapour and a snow of dry ice, both at the 1 atm sublimation point
    // of 194.65 K. The snow then sublimes *in the compartment*, and the compartment
    // pays for it. Per kilogram discharged the gas therefore receives
    //
    //     u = c_v,agent T_discharge - solidFraction (L_sub - R_agent T_sub)
    //
    // -- an internal energy, not an enthalpy, because the expansion work was done
    // by the storage pressure outside this model's boundary and the compartment is
    // rigid. `L_sub - R T` rather than `L_sub` because a closed rigid vessel books
    // energy at constant volume and the sublimation enthalpy is a constant-pressure
    // quantity; using the enthalpy here would over-state the sink by 6.9%.
    //
    // The number it comes to for CO2 is **-6.28 kJ per kilogram** -- negative, so
    // the discharge is a net internal-energy *sink* even counting the agent's own
    // internal energy. That is not a rounding: flooding a 1000 m^3 space to 35%
    // takes a tonne of CO2 and 194 MJ out of a gas whose whole heat capacity is
    // 0.88 MJ/K, and on the gas alone that is 127 K of chill and a pressure that
    // ends up *below* atmospheric. It is why the boundary term is not optional in
    // an agent run: the steel and everything bolted to it supply that heat, and
    // this model has them only as `wallConductance` against `wallTemperature`.
    // Everything the space *contains* is missing, so the transient chill here is a
    // ceiling rather than a prediction, and it is named rather than hidden.
    //
    // Setting `solidFraction` to zero and `dischargeTemperature` to `kTAmbient`
    // gives the pure-displacement case, which is what the concentration closed
    // forms are asserted on -- a mole fraction does not care what temperature the
    // gas is at, which is the reason that separation is clean.
    double dischargeTemperature = 194.65;   // K
    // The share of the discharged mass that arrives as solid. NFPA 12's figure for
    // CO2 flashed from storage to atmospheric.
    double solidFraction = 0.25;

    // Diagnostic, filled every step.
    double lastFlow = 0;        // kg/s actually leaving the bottles
    double lastCooling = 0;     // W the discharge takes out of the gas (may be < 0)
    double delivered = 0;       // kg out of the bottles over the system's life
};

// The mass of agent that would take `g` to `volumeFraction` by volume, kg, over
// and above what it is already carrying. Negative if it is already above it.
//
// The closed form, and it is pure ideal-gas arithmetic:
// `n_agent = y/(1-y) n_carrier`, and a mole count is `m R / R_universal`, so
//
//     m_agent(target) = y/(1-y) * (m_carrier R_air) / R_agent.
//
// Deliberately **not** the design quantity a rule would give you: SOLAS FSS ch. 5
// sizes a machinery-space bank on 40% of the gross volume with the free-gas volume
// taken at 0.56 m^3/kg and no credit for anything, which carries a margin for
// leakage during discharge. This is the sealed arithmetic the margin is measured
// from.
double agentMassForFraction(const GasCompartment& g, double volumeFraction);

// What the atmosphere in a compartment does to a person standing in it, at one
// height. Data, and the thresholds beside it, rather than a verdict.
//
// This exists because the roadmap item is about a damage-control *decision*. A
// machinery space at CO2's design concentration is not a space with a suppressed
// fire in it, it is a space that kills anyone still inside, and a simulator that
// models the suppression without modelling that is answering the easier question.
// The two numbers that make the point are `agentFraction` and `oxygenFraction`,
// and they say different things: at 35% CO2 the oxygen is still 13.6%, which is
// impaired but survivable -- and the CO2 is three and a half times its own lethal
// concentration. The agent kills first, and by its own action.
struct Exposure {
    double z = 0;                     // body frame, where it was sampled
    double temperature = kTAmbient;   // K
    double agentFraction = 0;         // by volume
    double oxygenFraction = kOxygenFractionAir;   // by volume
    double productFraction = 0;       // kg/kg of smoke in the gas
    bool   inUpperLayer = false;

    // Why, so a caller can report the cause and not only the outcome. Both are
    // set when both apply.
    bool   oxygenIncapacitating = false, oxygenLethal = false;
    bool   agentIncapacitating = false, agentLethal = false;

    bool incapacitating() const { return oxygenIncapacitating || agentIncapacitating; }
    bool lethal() const { return oxygenLethal || agentLethal; }
};

// Sample the atmosphere at a height in the body frame. Heights outside
// `[floorZ, ceilingZ]` are clamped into the space, because a breathing zone taken
// off a deck the compartment does not reach should not silently report ambient.
Exposure exposureAt(const GasCompartment& g, double z);

// A standing person's breathing zone above the deck, m. 1.5 m is the height the
// occupational-exposure literature samples at, and it is well inside the blanket a
// machinery-space CO2 flood makes.
inline constexpr double kBreathingZone = 1.5;

// ---------------------------------------------------------------------------
// The boundary the structure sees
// ---------------------------------------------------------------------------
//
// `docs/06-roadmap.md`'s Phase 4 milestone -- "an engine room fire that heats a
// bulkhead until it fails under the head of water behind it" -- and the one link
// of it that belongs in this file: what boundary condition a two-zone gas imposes
// on a `thermal::Problem`.
//
// **The fire already had the loss term. What it did not have was a wall.**
// `GasCompartment::wallConductance` and `wallTemperature` have been here since the
// suppression work, and `step()` relaxes each layer towards `wallTemperature` over
// its own wetted area exactly. But `wallConductance` was MQH's lumped `h_k` --
// `sqrt(k rho c / t)`, a *wall* conductance standing in for the whole path -- and
// `wallTemperature` was a constant, pinned at ambient for as long as the fire ran.
// A compartment could burn for an hour against a boundary that never got warm.
//
// With a conduction solve on the other side, both halves change and neither needs
// a new mechanism:
//
//   * `wallTemperature` becomes the steel's own surface temperature, which the
//     solve produces;
//   * `wallConductance` becomes the **gas-side film**, because the wall's own
//     resistance is no longer being lumped into it -- it is being solved.
//
// So the same `expm1` relaxation that was standing in for a boundary now *is* the
// boundary, and the energy it takes out of the gas is the energy the film puts
// into the steel. `wallExchange` below produces both ends of that from one set of
// temperatures, which is what stops them from being two answers.
//
// **Radiation is here and it is not an approximation.** `thermal.hpp` says plainly
// that radiation is not in the conduction solve, because it is nonlinear in a way
// conduction is not. It does not have to be:
//
//     sigma (T_g^4 - T_s^4) = sigma (T_g^2 + T_s^2)(T_g + T_s) (T_g - T_s)
//
// is an **identity**, so a film coefficient of `eps sigma (T_g^2+T_s^2)(T_g+T_s)`
// delivers the Stefan-Boltzmann flux exactly at the two temperatures it was formed
// at. Re-formed every coupling step, it is exact at every one of them; what it is
// not is exact *within* a step, and the size of that is a measurement rather than
// an argument -- `WallExchange::linearisationError` reports it. This is the same
// arrangement `thermal.hpp` already describes in prose at `gaussTemperature`,
// where a 900 C compartment is quoted through "an effective film of 200 W/(m^2 K),
// which is EN's 25 for convection plus about 175 for the radiation this file does
// not carry".

// Convection and the emissivity of the gas-to-steel exchange.
//
// 25 W/(m^2 K) is EN 1991-1-2 §3.1's coefficient on the fire-exposed side of a
// standard-fire boundary. 0.7 is the resultant emissivity, and it is the carbon
// steel figure -- EN 1993-1-2 §4.2.5.1, the material part, rather than the
// actions part quoted above it.
//
// **This used to derive the 0.7 from "the code's 0.8 surface emissivity ... after
// the configuration factor", and there is no 0.8 anywhere in this repository.**
// The only emissivity in the tree is this one. So the parenthetical explained a
// constant in terms of a second constant that does not exist, and 0.7/0.8 is
// 0.875, which is not a configuration factor either standard names -- EN's
// default is 1.0. The value is unchanged; what is removed is a derivation that
// could not be checked and would have misled anyone who tried.
//
// The clause attribution above is from the standards' structure and has not been
// re-read against the documents themselves; if the two parts are ever confirmed,
// this comment is where the confirmation belongs.
struct BoundaryFilm {
    double convective = 25.0;   // W/(m^2 K)
    double emissivity = 0.7;    // dimensionless, resultant
};

// Stefan-Boltzmann.
inline constexpr double kStefanBoltzmann = 5.670374419e-8;   // W/(m^2 K^4)

// `h_c + eps sigma (T_g^2 + T_s^2)(T_g + T_s)`, W/(m^2 K). Multiplying it by
// `T_g - T_s` reproduces `h_c (T_g - T_s) + eps sigma (T_g^4 - T_s^4)`.
//
// **Algebraically exactly, and numerically better than the thing it reproduces.**
// The two agree to 1.0e-13 relative over 280-1400 K wherever the two temperatures
// differ by a kelvin or more, and where they differ by less it is `T_g^4 - T_s^4`
// that is losing the digits, to the cancellation the factored form has none of.
// That is not a curiosity: a coupled solve spends its time *near* equilibrium, which
// is exactly where the difference form is worst, and at `T_g == T_s` this returns a
// flux of exactly zero rather than a rounding of one -- which is what the
// zero-heat-release control rests on.
double filmCoefficient(double gasKelvin, double surfaceKelvin,
                       const BoundaryFilm& params = {});

// The two-zone boundary condition, as a `thermal::Problem` wants it, plus the two
// numbers the gas wants back.
//
// **The split is by height, and that is the reason the model has two zones at all.**
// A bulkhead standing in a compartment with a smoke layer over it has hot gas on its
// upper part and cool air on its lower, and a single film at a mean temperature
// would heat the whole bulkhead evenly -- which is precisely the thing the head of
// water behind it is not doing, and the milestone turns on the two profiles being
// different shapes. Faces are assigned by their own centroid height against
// `GasCompartment::interfaceZ()`, in the **body** frame, which is where the
// interface is defined.
//
// **`bandHeight` exists because one coefficient over a whole layer is not good
// enough, and the size of that was measured rather than assumed.** A `Film` carries
// one coefficient for its whole face set, so a layer treated as one film has its
// radiative coefficient formed at the layer's *mean* surface temperature -- and the
// radiative coefficient goes as `T_s^2`, so on a bulkhead whose upper part is at
// 500 K and whose foot is held near ambient by the water behind it, the two ends of
// one layer want coefficients differing by half. Measured on the ferry's own
// bulkhead: 8.8 kW of a 110 kW exchange, 8%, at the worst step. Banding by height
// fixes it for nothing, because surface temperature on a bulkhead varies with height
// and almost not at all across it.
//
// Positive `bandHeight` therefore cuts the faces into horizontal bands of that
// height and gives each its own film. **The band index is a pure function of the
// face's own centroid**, so a caller that hands in the same face set every step gets
// the same films in the same order every step -- which is what makes
// `thermal::Solver::setFilm` usable at all, and what stops a coupled run from having
// to re-prepare (and so reset its own energy account) every time the smoke layer
// moves. Zero or negative gives exactly two films, split at the interface.
struct WallExchange {
    // Ready for `Problem::film`, in band order: with `bandHeight` positive, band 0
    // is the lowest. Empty bands are kept, so the indexing does not depend on which
    // heights happened to have faces at them.
    std::vector<thermal::Film> film;
    std::vector<double> area;      // m^2 per film
    std::vector<double> surface;   // K, the area mean each film's coefficient was formed at

    // What the gas should now carry, for the caller to write onto the
    // `GasCompartment`. Area-weighted over every face handed in.
    //
    // **This extrapolates**, and it is the honest weakness of the coupling: the
    // fire charges its boundary loss over the compartment's whole wetted enclosure,
    // and the faces handed in are one bulkhead of it. Writing these says "the rest
    // of the enclosure behaves like the part that is meshed". It is still a strict
    // improvement on a constant `h_k` against a boundary pinned at ambient, and it
    // is named rather than hidden.
    double wallTemperature = 0;    // K
    double wallConductance = 0;    // W/(m^2 K)
    double totalArea = 0;          // m^2

    // W into the solid at the temperatures handed in, by the films as built.
    double heat = 0;
    // The same, face by face, with each face's own surface temperature and the
    // exact `sigma (T_g^4 - T_s^4)`. `heat - exactHeat` is zero when every face in a
    // band is at the same temperature, which is what makes it a measure of the
    // *spread within a band* rather than of the linearisation, and what makes it
    // fall as `bandHeight` does.
    double linearisationError = 0;   // W
    double exactHeat = 0;            // W
};

// Build it. `surfaceKelvin` is one temperature per entry of `face`, which is what
// a caller gets by averaging `thermal::Solver::temperature()` over each face's four
// nodes. A mismatched or empty pair is refused: the result is empty, which
// `Problem` will treat as no boundary at all, rather than a film at a silently
// invented temperature.
WallExchange wallExchange(const GasCompartment& gas,
                          const std::vector<thermal::BoundaryFace>& face,
                          const std::vector<double>& surfaceKelvin,
                          const BoundaryFilm& params = {}, double bandHeight = 0.0);

// ---------------------------------------------------------------------------
// The account
// ---------------------------------------------------------------------------

// Where the mass, the energy and the products went. Same discipline as
// `zone::SolveResult`'s energy account and `Ship`'s flooding: a model that cannot
// say what crossed its own boundary is a model that cannot be checked.
//
// `zone.hpp` records its residual closing on unit fixtures and running to -102%
// on a real tearing case, which is exactly why `tests/test_fire.cpp` asserts this
// one on the ferry's engine rooms and not only on a box.
struct Account {
    double heatReleased = 0;     // J the design fires put into the gas
    double radiativeLoss = 0;    // J the design fires took straight back out
    double wallLoss = 0;         // J conducted into the boundaries
    double suppressionCooling = 0;  // J the drenchers took out of the gas
    double enthalpyIn = 0;       // J in across the model's boundary (outside air)
    double enthalpyOut = 0;      // J out across it
    double massIn = 0;           // kg in across it
    double massOut = 0;          // kg out across it
    double productsGenerated = 0;   // kg
    double productsOut = 0;         // kg across the boundary

    // The agent books, kept apart from the air's for the reason the water's are:
    // a species that leaks is a species whose account has to close on its own, or
    // "the concentration fell" and "the model lost some" are the same observation.
    // `agentDischarged` is a mass *source* inside the boundary and therefore
    // appears in the mass residual too.
    double agentDischarged = 0;  // kg out of the bottles and into the gas
    double agentEnergy = 0;      // J the discharge added to the gas; negative for CO2
    double agentIn = 0;          // kg of agent in across the model's boundary
    double agentOut = 0;         // kg of agent out across it
    double initialAgent = 0, agent = 0;   // kg at attach(), and now

    // The water books, kept on the same terms as the gas ones and for the same
    // reason. All three are *crossings*, not holdings: what is standing on the
    // deck at any moment is the ship's own `Compartment::waterVolume`, because
    // `applyTo` hands the water over and the flooding solve owns it from then on
    // -- it may drain through the ship's own openings, or be pumped, or run to
    // the low side and out. A second copy of that holding kept here would be the
    // "two answers that agree until the day they do not" this file has already
    // deleted once.
    double waterDelivered = 0;   // kg out of the nozzles
    double waterEvaporated = 0;  // kg of it that left as steam
    double waterDrained = 0;     // kg over the freeing ports, negative if they admitted

    double initialEnergy = 0, initialMass = 0;   // at attach()
    double energy = 0, mass = 0, products = 0;   // now

    double energyResidual() const {
        return (heatReleased - radiativeLoss + enthalpyIn - enthalpyOut - wallLoss -
                suppressionCooling + agentEnergy) -
               (energy - initialEnergy);
    }
    double massResidual() const {
        return (massIn + agentDischarged - massOut) - (mass - initialMass);
    }
    double productsResidual() const { return (productsGenerated - productsOut) - products; }
    double agentResidual() const {
        return (agentDischarged + agentIn - agentOut) - (agent - initialAgent);
    }

    // What suppression handed to the flooding network, kg. Delivered, less what
    // evaporated, less what the freeing ports took back overboard -- which is
    // exactly what `applyTo` has written into `Compartment::waterVolume` over the
    // model's life, and `tests/test_fire.cpp` asserts the two agree to the bit.
    double waterToShip() const { return waterDelivered - waterEvaporated - waterDrained; }

    // The residuals above are absolute. These are the fractions worth quoting,
    // normalised by the largest single term rather than by the net -- a net that
    // is nearly zero would make any residual look enormous, which is one way a
    // conservation check tells you nothing.
    double energyResidualFraction() const;
    double massResidualFraction() const;
};

// What one call to `step()` did, for a caller that wants to watch rather than
// poll the state.
struct StepResult {
    double time = 0;             // s, model time after the step
    double heatRelease = 0;      // W, summed over the design fires now
    // Cooling as **applied** on the last internal step, after the cap that stops
    // a drencher taking the layer below the temperature of its own water. Same
    // discipline as `entrainment` below: the demand is
    // `flow * (cp dT + e L)` and reporting that instead would be reporting a
    // number the model did not use. On the ferry the two differ by two orders of
    // magnitude, because a deck drencher can absorb far more than any fire in the
    // space is releasing.
    double suppressionCooling = 0;  // W
    // Entrainment as **applied** on the last internal step -- after the taper
    // that shuts the plume off as the cool layer is consumed, and after the cap
    // that stops it draining more than the layer holds. `Model::totalEntrainment`
    // is the raw correlation; the difference between the two is the whole of what
    // happens when the smoke layer reaches the floor, so reporting only the
    // correlation would be reporting a number the model did not use.
    double entrainment = 0;      // kg/s
    // The agent systems' discharge on the last internal step, kg/s, after the cap
    // that stops a bank delivering more than it holds -- applied, not demanded, on
    // the same terms as everything else here.
    double agentDischarge = 0;   // kg/s
    // The lowest oxygen availability any design fire saw on the last internal step
    // (1 in clean air, 0 once the space is inert). Published because "the fire went
    // out" and "the fire's curve ended" are the same heat release and different
    // events.
    double oxygenAvailability = 1.0;
    int    substeps = 0;         // how many internal steps the accuracy cap took
    int    pressureSweeps = 0;   // Gauss-Seidel sweeps the last substep's solve took
    // The pressure solve failed to bracket a root in 80 doublings. Should never
    // happen -- the residual is monotone and unbounded in both directions -- and
    // is published rather than swallowed so that "should never" can be asserted.
    bool   pressureSolveCapped = false;
    // `maxSubsteps` ran out before the whole `dt` was taken, so the gas is **short
    // of the time it was asked for** and otherwise indistinguishable from a step
    // that took all of it. Reachable rather than theoretical: a rejected trial
    // spends a budget slot and halves `h` *without* advancing `remaining`, so the
    // budget bounds trials and not committed time, and a defect that drives the
    // rejection test into a corner does not fail -- it under-advances quietly.
    // `les::StepResult::incomplete` was added for a loop of exactly this shape and
    // named this file when it was.
    //
    // The exposure is concrete. `tools/bulkhead_probe` steps the gas, the
    // conduction solve and the ship by one `coupling` each and prints all three
    // under a `t` it computes itself as `(step + 1) * coupling`; a short gas step
    // puts three different model times on one row of a trace that stays perfectly
    // smooth.
    bool   incomplete = false;
};

// ---------------------------------------------------------------------------
// The model
// ---------------------------------------------------------------------------

class Model {
public:
    std::vector<GasCompartment> gas;
    std::vector<DesignFire>     fires;
    std::vector<Vent>           vents;
    std::vector<Drencher>       drenchers;
    std::vector<Scupper>        scuppers;
    std::vector<AgentSystem>    agents;
    Account                     account;

    double time = 0;             // s of model time
    // Largest internal step. The vent flow is stiff near equilibrium -- the
    // orifice law's `sqrt` has an infinite derivative there -- so `step()`
    // subdivides rather than trusting the caller's tick.
    double maxSubstep = 0.25;    // s
    // A trial substep is **rejected and halved** if it would move any
    // compartment's mass or energy by more than this fraction. Rejected before
    // anything is committed, including the account, so a rejected step leaves no
    // trace and the committed rates were always evaluated at the committed step.
    double maxRelativeChange = 0.02;
    int    maxSubsteps = 100000;

    // Build the gas model over `shipCompartments` of `ship`, deriving each
    // tracked compartment's box from its mesh and every vent from the openings
    // that touch it. Openings between a tracked compartment and an untracked one
    // are dropped -- an untracked space has no gas state to exchange with, and
    // pretending it is the atmosphere would invent an infinite reservoir behind a
    // bulkhead.
    void attach(const Ship& ship, const std::vector<int>& shipCompartments);

    // Index into `gas` for a ship compartment, or -1.
    int gasIndexOf(int shipCompartment) const;
    int findGas(std::string_view name) const;

    // Zero the account and take the current gas state as its baseline. `attach()`
    // does this; a model assembled by hand -- a bare compartment in a test, a
    // fire scenario with no ship behind it -- has to say so, and an account whose
    // baseline is zero reports the whole initial internal energy as a residual.
    void resetAccount();

    // Advance the gas by `dt`. Reads `ship` for what is against each opening --
    // an opening with water on either side belongs to the flooding solve and is
    // skipped -- and **writes nothing to it**.
    StepResult step(double dt, const Ship& ship, const Sea& sea);

    // The write path back into the flooding network, and it is opt-in.
    //
    // `Ship`'s gas is isothermal at `kTAmbient` by construction:
    // `updateInternalFreeSurfaces` recomputes `airPressure = airMass R T_amb / V`
    // every tick, so there is no way to hand it a hot compartment. What it
    // actually *uses* the gas for is the pressure at an orifice, so this writes
    // the mass that reproduces this model's pressure under the ship's own
    // formula: `m_eq = p_fire V / (R T_amb)`.
    //
    // That makes the flooding network see the right pressure and makes
    // `Compartment::airMass` a pressure proxy in a burning compartment. It is a
    // real distortion of the ship's books and it is the reason this is a separate
    // call rather than something `step()` does. Writing the *delta* since the
    // last call, rather than the absolute value, is what keeps a zero-heat-release
    // run bit-identical: the delta is then exactly 0.0 and `x += 0.0` changes
    // nothing.
    //
    // **This is also the path suppression water takes**, and it is the whole of
    // it: the net of what the drenchers landed and the freeing ports took back
    // goes into `Compartment::waterVolume`, and every stability consequence
    // follows from the ship's own solve with nothing else added. The water half
    // is written as a pending delta that this call consumes, so a model with no
    // drenchers writes nothing at all -- not `+= 0.0`, but no store, which is
    // what keeps the exact control exact even for a compartment whose water
    // volume is a negative zero.
    void applyTo(Ship& ship);

    // Water this model owes the ship, m^3, indexed like `gas`. Published because
    // it is the number `applyTo` is about to write and because a caller stepping
    // the gas faster than the ship needs to see it accumulate.
    const std::vector<double>& pendingWater() const { return pendingWater_; }

    // Total plume entrainment now, kg/s. Published because "the layer descended"
    // is not evidence that entrainment did it -- the fire's own expansion pushes
    // the interface down too, about two orders of magnitude slower.
    double totalEntrainment(double atTime) const;

    // Problems with the model definition, in the same spirit as `Ship::validate`:
    // a bad definition does not crash, it quietly produces a wrong answer.
    std::vector<std::string> validate() const;

private:
    // Attempts one internal step. Returns false, having changed nothing at all,
    // when `dt` would breach `maxRelativeChange`.
    bool substep(double dt, const Ship& ship, const Sea& sea, StepResult& out);

    // The air mass this model last wrote into each ship compartment, so `applyTo`
    // can write a delta. Indexed like `gas`.
    std::vector<double> appliedMass_;
    // Water the drenchers have landed and the freeing ports have not yet taken
    // back, m^3, waiting for `applyTo`. Signed: a submerged freeing port admits.
    std::vector<double> pendingWater_;
};

// The correlations the tests check against, exposed because a caller wants them
// for a design check and because a correlation kept private is a correlation
// nobody can disagree with.

// McCaffrey-Quintiere-Harkleroad: the steady upper-layer temperature rise of a
// well-ventilated compartment fire, K.
//
//   `dT = 6.85 (Q^2 / (Ao sqrt(Ho) hk At))^(1/3)`
//
// Arguments in SI: `heatRelease` W, `ventArea` m^2, `ventHeight` m,
// `wallConductance` W/(m^2 K), `wallArea` m^2. The published form is in kW and
// kW/(m^2 K); the conversion is here, once.
double mqhTemperatureRise(double heatRelease, double ventArea, double ventHeight,
                          double wallConductance, double wallArea);

// Thomas's flashover correlation, W: `Q = 7.8 At + 378 Ao sqrt(Ho)` in kW. An
// independent fit to the same phenomenon as MQH inverted at a 500 K rise, which
// is what makes it worth carrying -- two correlations that agree are evidence,
// one is a definition.
double thomasFlashoverPower(double wallArea, double ventArea, double ventHeight);

}  // namespace sim::fire
