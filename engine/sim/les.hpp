// SPDX-License-Identifier: MIT
//
// A resolved-flow model for **one** compartment -- `docs/06-roadmap.md` Phase 4's
// "LES promotion for the local compartment", and the fidelity tier above
// `fire.hpp`'s two zones.
//
// `fire.hpp` gives every tracked compartment two well-mixed control volumes: a hot
// upper layer and a cool lower one. That is the right first answer and it is cheap.
// What it cannot do is *compute* anything the two-zone description asserts -- a
// plume, a ceiling jet, a horizontal temperature gradient under a deckhead, or a
// pocket of cool air behind a stack of containers. Each of those is a field, and a
// zone model has no fields. This file gives one compartment fields, on a coarse
// regular grid, and hands the answer back as two zones when it is done.
//
// --- 1. What this is, and what it deliberately is not --------------------------
//
// **It is a coarse low-Mach finite-volume solve on a regular staggered grid.** Not
// a research LES: there is no body-fitted mesh, no dynamic procedure, no wall
// model, no combustion chemistry and no radiative transport. The subgrid closure
// is constant-coefficient Smagorinsky, which is the cheapest thing that is a
// closure at all rather than an absence of one -- it supplies the grid-scale
// dissipation a buoyant plume on a 0.5 m mesh has nowhere else to get, and it has
// exactly one constant. A dynamic Smagorinsky needs a test filter and a plane
// average, which on a compartment with six walls and no homogeneous direction is
// more machinery than this earns.
//
// **Low Mach, and that is the load-bearing choice.** A fire drives the gas at
// about 1 m/s against a sound speed of 340; a fully compressible solve would spend
// its whole budget on acoustics it does not want, at a step of `dx/c` -- 1.5 ms on
// a 0.5 m cell. So the pressure is split the way every fire code splits it: one
// **thermodynamic** pressure for the whole box, uniform in space, and a dynamic
// pressure that appears only in the momentum equation and carries no thermodynamics
// at all.
//
// That split is not merely convenient here. It is *the two-zone model's own
// closure*, resolved:
//
//   * At one uniform pressure, `u = p / (gamma - 1)` per unit volume everywhere, so
//     the **internal energy density is uniform** and the total is `p V / (gamma-1)`.
//   * `fire.hpp` states the same fact as `V_u / V = U_u / U` and calls it "the
//     single load-bearing algebraic fact in this file". It is the same identity.
//   * Therefore the resolved field's total energy is one scalar, exactly as the
//     two-zone model's is, and **energy transfers across the fidelity boundary with
//     nothing to interpolate**. Whatever the flow does, `E` is `E`.
//
// The consequence for the fields is that **mass is the field and temperature is
// derived**, not the other way round: `rho_c = m_c / V_c` and
// `T_c = U_c / (m_c c_v)` with `U_c = E V_c / V`. There is no separate temperature
// equation to drift out of step with the density, which is the failure mode of the
// obvious variable-density-incompressible arrangement -- there `rho` and `T` are
// advected apart and the equation of state stops holding after a few hundred steps.
// Here it holds identically, at every cell, at every step, by construction.
//
// Heat therefore enters the *velocity* field rather than a temperature field. The
// continuity equation at uniform pressure is
//
//     div u = (gamma-1)/(gamma p) * q'''  -  1/(gamma p) * dp/dt
//
// and integrating that over a sealed box gives `dp/dt = (gamma-1) Q / V` -- which
// is `fire.hpp`'s own sealed closed form, `p(t) = p_0 + (gamma-1) E / V`,
// recovered rather than imposed. The local source is then the *deviation* of the
// heating from its box mean, which sums to exactly zero over the box and so is a
// compatible right-hand side for the all-Neumann projection below. A cell that is
// heated expands, pushes mass out, and thereby gets lighter and hotter: the
// buoyancy is a consequence of the mass leaving, not a separate model of it.
//
// **Boussinesq in the inertia, not in the buoyancy.** The projection is
// constant-coefficient -- the momentum equation divides by a single reference
// density -- while the buoyant acceleration carries the full `(rho_ref - rho)/rho_ref`
// with no linearisation in it. That is the usual compromise and it costs an O(1)
// error in the *inertia* of a hot parcel, which a coarse grid does not resolve
// anyway, in exchange for a Poisson operator with constant coefficients.
//
// **The reference density is a linearisation point, not a free constant**, and
// getting that wrong is easy: `g (rho_ref - rho) / rho_ref` is `g - g rho / rho_ref`,
// so the reference sets the *scale* of the buoyancy as well as its offset. Only the
// offset is absorbed by the projection -- a spatially uniform body force in a closed
// box is exactly balanced by the pressure, which is why a box uniformly at 350 K
// does not move even though its density is nowhere near ambient. The scale is a real
// Boussinesq sensitivity and it is measured rather than argued: moving the reference
// from the box's own mid-density to ambient, about 35%, moves the peak velocity one
// step later by **0.27%**.
//
// Zero takes the box's mid-density, `(rho_min + rho_max)/2`. Two reasons, and the
// second is the load-bearing one: it centres the linearisation on the gas actually
// present, and it is **exactly** the cell density when the box is uniform, so the
// buoyant acceleration on every face of a quiescent compartment is exactly `0.0`
// rather than a rounding of it. The box *mean* would be the obvious choice and
// cannot have that property -- `N rho` is not in general a representable double, so
// even a compensated sum of `N` identical cells need not return exactly `N` times
// one of them. That is the same discipline `fire.hpp`'s `VentSide::gaugeAt` follows,
// for the same reason: an exact control is worth more than a tidy formula.
//
// --- 2. What crosses the fidelity boundary, and what is conserved ---------------
//
// `promote()` lays a two-zone compartment onto the grid; `demote()` reads the grid
// back as two zones. Both are exact statements about mass, energy and products,
// and neither has a tuning constant in it.
//
//   * **Down**: each layer's mass and products are spread uniformly by volume over
//     the region it occupies, split at `GasCompartment::interfaceZ()`; a cell the
//     interface passes through gets a volume-weighted share of both. The energy
//     needs no distribution at all -- it is one scalar and it is copied.
//   * **Up**: the two-layer equivalent of a resolved field is the profile that has
//     the field's own ceiling and floor densities and the field's own total mass:
//
//         rho_u (H - z_i) + rho_l (z_i - z_f) = M / A
//
//     one linear equation in `z_i`, closed form, no iteration and no fitted
//     constant. On a field that came *from* two zones it returns the interface it
//     came from **exactly**, which is the algebra: the bracket that multiplies
//     `(z_i0 - z_i)` is `m_u/(H - z_i0) - m_l/(z_i0 - z_f)`, and that is
//     `A(rho_u - rho_l)`, nonzero unless the box is uniform.
//   * The upper layer's mass is then `rho_u A (H - z_i)` -- **reconstructed from
//     that profile, not integrated over the cells** -- and the lower is taken as
//     the remainder, so the two halves add to the total in floating point and not
//     merely in algebra. That is the same trick `fire.cpp`'s `layerSplit` uses, and
//     the reconstruction is not a stylistic preference: promotion has to smear the
//     interface across the cell row it falls inside, and integrating the cell
//     masses above `z_i` then comes back **9% wrong on a field it has just
//     promoted**, because no weighted sum of a mixed row can pull the two layers
//     back apart. The reconstruction can, because the smear leaves the top row, the
//     bottom row and the total untouched and those are the only three numbers it
//     reads. Measured: 9e-2 against 1e-16.
//   * Products follow the same route through the *mass fraction* the top row
//     carries, `s_u = y_u m_u`, for the same reason and with the same exactness.
//   * A **uniform** box determines no interface at all: `rho_u == rho_l` and the
//     equation above is `0 = 0`. That is not a numerical difficulty, it is the
//     physics -- every split of a uniform gas is the same gas -- so the interface
//     is *held* at its last determined value. `fire::GasCompartment::fillAmbient`
//     seeds a token upper layer for the same reason, and holding is what makes the
//     round trip exact on a compartment nobody has lit.
//
// The round-trip error is therefore geometric rounding and nothing else, and it is
// measured rather than assumed: `tests/test_promotion.cpp` reports it and asserts
// at what it measured.
//
// --- 3. What is not here, named rather than hidden -----------------------------
//
//  1. **Vents.** The resolved compartment is sealed for the duration of its
//     promotion. `fire.hpp`'s vent integral wants a gauge pressure profile against
//     a *layer* structure, and handing it a field means deciding what a doorway
//     does with a ceiling jet arriving at an angle -- a real question and a
//     separate one. A compartment is promoted for what happens *inside* it; the
//     ventilation answer is still the two-zone one, and the honest consequence is
//     that a promotion should be short next to the compartment's own ventilation
//     time constant. `demote()` hands the state back so the vents can have it.
//  2. **Kinetic energy.** The account is internal energy only. At 1 m/s and 1.2
//     kg/m^3 the kinetic energy density is 0.6 J/m^3 against 253 kJ/m^3 of internal
//     energy -- 2.4e-6 -- so the work the buoyancy does is not booked. It is named
//     because "the energy account closes" would otherwise be a claim about a
//     quantity that was never in it.
//  3. **Combustion.** The heat release is a design fire, exactly as in `fire.hpp`.
//     What *is* resolved is the entrainment: `fire::Plume`'s Heskestad correlation
//     is not used here at all, because the plume entrains by advecting. That is the
//     one thing this file computes that the two-zone model can only correlate.
//  4. **Radiation, soot optics, and a second species.** Products ride as a passive
//     tracer, as they do in `fire.hpp`.
//  5. **A heeled compartment.** The grid is axis-aligned in the body frame, so
//     gravity is `-z`. Same limitation, same place, as the two-zone interface.
//  6. **The plan shape.** `fire::GasCompartment` carries a footprint *area* and a
//     *perimeter* and no shape at all, so the grid is the rectangle that has both
//     -- the roots of `t^2 - (P/2) t + A = 0`. A real ship compartment's perimeter
//     is its bounding box's while its area is the prismatic equivalent, so those
//     roots are often complex; then a square of side `sqrt(A)`, and `problems`
//     says so. The area is always right, because the area is what the interface
//     height and the wall exchange are computed from.
//
// SI units, body frame, temperature in kelvin, per CLAUDE.md.
#pragma once

#include "fire.hpp"

#include <string>
#include <vector>

namespace sim::les {

// ---------------------------------------------------------------------------
// The grid
// ---------------------------------------------------------------------------

// A uniform axis-aligned box of cells with velocities on the faces (MAC/staggered).
// Staggered rather than collocated because a collocated projection on a regular
// grid has a checkerboard null space, and a checkerboard is exactly the shape a
// plume's pressure field would put into it.
struct Grid {
    int    n[3]{0, 0, 0};    // cells along x, y, z
    double h[3]{0, 0, 0};    // m, cell edge along each axis
    double lo[3]{0, 0, 0};   // m, the low corner in the body frame

    int cells() const { return n[0] * n[1] * n[2]; }
    bool empty() const { return cells() <= 0; }
    int cell(int i, int j, int k) const { return (k * n[1] + j) * n[0] + i; }
    double cellVolume() const { return h[0] * h[1] * h[2]; }
    // Cells times cell volume, **not** the product of the three side lengths. The
    // two differ in the last bits, and the divergence source below is only a
    // compatible right-hand side if the volume the mean is taken over is the one
    // the cells add up to.
    double volume() const { return cellVolume() * cells(); }
    // Each factor widened before the multiply, as `flip.hpp` does for the same
    // quantity. `n[0] * n[1]` formed the product in `int` first, which is safe only
    // because `gridFor` clamps each axis to 512 -- and `Field` is a public struct
    // whose `Grid` a caller can fill in directly.
    double planArea() const {
        return h[0] * h[1] * static_cast<double>(n[0]) * static_cast<double>(n[1]);
    }
    double centre(int axis, int index) const {
        return lo[axis] + (static_cast<double>(index) + 0.5) * h[axis];
    }

    // The staggered component along `axis` has one more face than cells along that
    // axis and exactly as many as cells along the others.
    int faceCount(int axis, int along) const { return n[along] + (along == axis ? 1 : 0); }
    int faces(int axis) const {
        return faceCount(axis, 0) * faceCount(axis, 1) * faceCount(axis, 2);
    }
    int face(int axis, int i, int j, int k) const {
        return (k * faceCount(axis, 1) + j) * faceCount(axis, 0) + i;
    }
    double faceArea(int axis) const { return cellVolume() / h[axis]; }
};

// The plan rectangle a footprint area and a perimeter imply: `L + W = P/2` and
// `L W = A`, so the sides are the roots of `t^2 - (P/2) t + A = 0`. Returns false,
// and a square of side `sqrt(A)`, when the discriminant is negative -- which is the
// ordinary case on a ship, where the perimeter is the bounding box's and the area
// is the prismatic equivalent, so no rectangle has both.
bool planRectangle(double floorArea, double perimeter, double* lengthX, double* lengthY);

struct Params {
    // Target cell edge. The grid takes the nearest whole number of cells on each
    // axis, so the cells the box actually gets are near-cubic but not cubic.
    double cellSize = 0.5;      // m
    int    maxCells = 12000;    // the budget; the cell size grows to fit it
    int    minPerAxis = 3;      // fewer than three cells across resolves nothing

    // Constant-coefficient Smagorinsky: `nu_t = (Cs Delta)^2 |S|`, `Delta` the cube
    // root of the cell volume. 0.20 is the classic value for free shear; a fire
    // plume is a free shear flow and there is no wall layer resolved here to want
    // anything else.
    double smagorinsky = 0.20;
    double viscosity = 1.5e-5;  // m^2/s, air's own; a floor under the eddy value

    // The advective Courant number a substep is allowed to reach before it is
    // rejected and halved, and the largest substep whatever the velocities say.
    double courant = 0.40;
    double maxSubstep = 0.25;         // s
    // A substep is rejected, having changed nothing, if it would move any cell's
    // mass by more than this fraction of what the cell holds. Same discipline as
    // `fire::Model::maxRelativeChange`, and here it is what actually sets the step
    // at ignition, when a cell in the flame is expanding through its own mass in
    // under a second.
    double maxRelativeChange = 0.05;
    int    maxSubsteps = 200000;

    // The projection. Red-black SOR: order-independent within a colour, which is
    // what makes the answer a pure function of the field rather than of the sweep
    // order, and therefore what lets a mirror-symmetry test be an exact statement.
    int    projectionSweeps = 400;
    double projectionTolerance = 1e-7;   // relative, on the divergence residual
    double relaxation = 1.7;             // SOR factor

    // kg/m^3. Zero takes the box's own mid-density; see the header for what the
    // reference is and is not, and for what moving it costs.
    double buoyancyReference = 0;
};

// ---------------------------------------------------------------------------
// The field
// ---------------------------------------------------------------------------

// A heat release with a place and a size. `fire::DesignFire` carries a `baseZ` and
// a diameter but no horizontal position, because a zone model has nowhere to put
// one; here the offset is from the compartment's plan centre.
struct HeatSource {
    double offsetX = 0, offsetY = 0;   // m from the plan centre
    double baseZ = 0;                  // m, body frame
    double diameter = 1.0;             // m
    double flameHeight = 0;            // m above the base over which the heat is released
    double power = 0;                  // W into the gas, already net of any radiative loss
    double productRate = 0;            // kg/s of combustion products
};

// One resolved compartment. Everything a `fire::GasCompartment` holds, resolved,
// plus the velocity field it does not have.
struct Field {
    Grid grid;
    double floorZ = 0, ceilingZ = 0;   // m, body frame; the grid spans them in z
    // The boundary, on exactly the terms `fire::GasCompartment` states it: a
    // gas-side film and the temperature of the steel behind it.
    double wallConductance = 30.0;     // W/(m^2 K)
    double wallTemperature = kTAmbient;

    // **The whole thermodynamic state is this one scalar**, because the internal
    // energy density is uniform at a uniform pressure. See the header.
    double energy = 0;                 // J, internal, over the whole box
    std::vector<double> mass;          // kg per cell
    std::vector<double> products;      // kg per cell
    std::vector<double> velocity[3];   // m/s on the faces of each component

    // The two-layer interface this field last determined, held while it does not
    // determine one. Diagnostic while resolved; the thing `demote()` splits on.
    double interfaceZ = 0;             // m
    double time = 0;                   // s of model time

    std::vector<std::string> problems;

    bool empty() const { return grid.empty() || mass.empty(); }
    double pressure() const;                  // Pa, uniform over the box
    double cellEnergy() const;                // J, the same in every cell
    double totalMass() const;                 // kg
    double totalProducts() const;             // kg
    double density(int cell) const;           // kg/m^3
    double temperature(int cell) const;       // K
    double meanTemperature() const;           // K, mass weighted: E / (M c_v)
    // The plan-mean density of a horizontal cell row, kg/m^3.
    double rowDensity(int k) const;
};

// ---------------------------------------------------------------------------
// The account
// ---------------------------------------------------------------------------

// Where the mass, the energy and the products went, on the same terms as
// `fire::Account`.
//
// **The energy line is weaker than it looks and saying so is the point.** The
// total internal energy *is* a state variable here, so "the energy account closes"
// is very nearly a statement about one `+=`. What is not by construction, and what
// the tests actually assert, is that the state matches the design fire's own
// analytic release through the sealed closed form `p = p_0 + (gamma-1) E / V`, and
// that promotion and demotion move it without loss. The **mass** line is the one
// that checks the transport: every face flux is applied to the two cells it joins
// with opposite signs, so a sealed box holds its mass to summation error, and an
// off-by-one in a neighbour index does not.
struct Account {
    double heatReleased = 0;       // J the sources put into the gas
    double wallLoss = 0;           // J the boundary took out
    double productsGenerated = 0;  // kg

    double initialEnergy = 0, initialMass = 0;
    double energy = 0, mass = 0, products = 0;

    double energyResidual() const {
        return (heatReleased - wallLoss) - (energy - initialEnergy);
    }
    double massResidual() const { return mass - initialMass; }
    double productsResidual() const { return productsGenerated - products; }

    // Normalised by the largest term including the state, for the reason
    // `fire::Account` gives: a run whose fluxes are all zero would otherwise report
    // its own last-bit drift as a whole-number fraction.
    double energyResidualFraction() const;
    double massResidualFraction() const;
};

struct StepResult {
    double time = 0;              // s of model time after the step
    int    substeps = 0;
    int    rejections = 0;        // trial substeps thrown away, having changed nothing
    int    projectionSweeps = 0;  // on the last accepted substep
    double projectionResidual = 0;   // 1/s, what the last projection left
    bool   projectionCapped = false; // it ran out of sweeps before its tolerance
    // The substep budget ran out before the whole `dt` was taken, so the field is
    // **short of the time it was asked for**. Published rather than swallowed for
    // the reason `fire::StepResult::pressureSolveCapped` is: "should never happen"
    // is only worth saying if something can assert it.
    bool   incomplete = false;
    double courant = 0;           // the worst advective Courant number committed
    double peakSpeed = 0;         // m/s
    double heatRelease = 0;       // W, summed over the sources now
};

// ---------------------------------------------------------------------------
// Across the boundary
// ---------------------------------------------------------------------------

// The grid a compartment gets: its own heights in z, the plan rectangle its area
// and perimeter imply in x and y, and a cell count under `Params::maxCells`.
Grid gridFor(const fire::GasCompartment& gas, const Params& params = {});

// Cells a promotion here would cost, without building it -- the budget's unit, on
// the same terms as `promotion::estimateElements`.
int estimateCells(const fire::GasCompartment& gas, const Params& params = {});

// Two zones down onto the grid. Exact in mass, energy and products.
Field promote(const fire::GasCompartment& gas, const Params& params = {});

// The two-layer interface a field implies: the height at which a profile carrying
// the field's own top-row and bottom-row densities holds the field's own mass.
// Returns `held` unchanged when the field is uniform enough that the equation
// degenerates -- see the header; a uniform gas has no interface to find.
double equivalentInterface(const Field& field, double held);

// The grid back up into two zones, writing mass, energy and products onto `gas`
// and leaving its geometry and its boundary alone. Exact in all three: the upper
// layer is summed and the lower is the remainder.
void demote(const Field& field, fire::GasCompartment& gas);

// The heat sources one gas compartment's design fires impose at `atTime`, placed at
// the plan centre because a `fire::DesignFire` carries no horizontal position, and
// spread over the flame height `fire::Plume` reports rather than over the pan --
// which is where the heat is actually released, and the difference is a factor of
// several in the source density a cell at the fire base has to absorb.
std::vector<HeatSource> sourcesFor(const fire::Model& model, int gasCompartment, double atTime);

// Zero the account and take the field as its baseline.
void resetAccount(const Field& field, Account& account);

// Advance the field by `dt` against a fixed set of sources. Subdivides internally,
// rejecting and halving a trial substep that would move any cell too far -- and a
// rejected substep leaves no trace at all, in the field or in the account.
StepResult step(Field& field, double dt, const std::vector<HeatSource>& sources,
                const Params& params, Account& account);

// Problems with a field's definition, in the spirit of `fire::Model::validate`.
std::vector<std::string> validate(const Field& field);

// ---------------------------------------------------------------------------
// Alpert's ceiling jet -- the correlation the promotion criterion reads
// ---------------------------------------------------------------------------

// Alpert (1972), the maximum ceiling-jet excess temperature under an unconfined
// ceiling `height` above the fire base, at radius `radius` from the plume axis:
//
//   r/H <= 0.18:  dT = 16.9 Q^(2/3) / H^(5/3)
//   r/H >  0.18:  dT = 5.38 (Q/r)^(2/3) / H
//
// Published in kW and metres; the conversion is inside, once, for the reason
// `fire::Plume`'s is.
double alpertCeilingJetRise(double heatRelease, double radius, double height);

// The **ratio** of the two: how much hotter the jet is over the fire than it is at
// `radius`. The heat release cancels out of it exactly --
//
//     16.9 Q^(2/3) H^(-5/3)  /  [5.38 Q^(2/3) r^(-2/3) H^(-1)]  =  (16.9/5.38) (r/H)^(2/3)
//
// -- which is what makes it a statement about the *compartment* and not about the
// fire in it, and therefore usable as a promotion criterion that a bigger fire does
// not simply switch on. The two branches meet at `r/H = 0.18` to within **0.1427%**
// -- derived from the published coefficients and measured, not remembered, exactly
// as `fire::Plume`'s two entrainment branches meet to 1.78% -- so the ratio passes
// through one at the crossover and this function is continuous there to that much.
double ceilingJetSpread(double radius, double height);

inline constexpr double kAlpertNear = 16.9;        // K on kW^(2/3) m^(-5/3)
inline constexpr double kAlpertFar = 5.38;         // K on kW^(2/3) m^(-2/3) m^(-1)
inline constexpr double kAlpertCrossover = 0.18;   // r/H at which the branches swap
// The far coefficient's radius exponent read `-1/3` and is `-2/3`: the correlation
// above is `5.38 (Q/r)^(2/3) / H`, and `les.cpp`'s `cbrt(qkw*qkw/(r*r))` is the
// same. It is checkable without either -- two branches of one correlation both turn
// `kW^(2/3)` times a length group into kelvin, so their length groups must have the
// same dimension, and `m^(-2/3) m^(-1)` is `m^(-5/3)`, which is what the near
// coefficient carries. `-1/3` would have made them different units.

}  // namespace sim::les
