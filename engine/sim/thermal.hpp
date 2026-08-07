// SPDX-License-Identifier: MIT
//
// Implicit heat conduction on the structural mesh -- the first item of
// `docs/06-roadmap.md` Phase 4 that everything else in that phase needs.
//
// The milestone Phase 4 is aimed at is "an engine room fire that heats a
// bulkhead until it fails under load". Three things stand between here and
// there: a compartment fire that says what the gas does, a conduction solve that
// says what the steel does, and a strength model that says what the steel then
// carries. This file is the middle one, and it is deliberately *only* the middle
// one -- there is no combustion here, no species transport, no radiation view
// factors and no coupling back into `solid_shell.hpp`'s constitutive law. Each of
// those is a separate roadmap item with its own validation, and one conduction
// solve that is right is worth more than three that are sketched.
//
// --- Why implicit, and a correction to the usual argument ----------------------
//
// The usual argument is that explicit conduction, stable only below
// `rho c h^2 / (2 k)`, is hopeless on ship plating because that limit is
// *milliseconds*. **It is not, and the figure is wrong by three orders of
// magnitude.** Steel's thermal diffusivity is
// `alpha = 53.3 / (7850 * 440) = 1.545e-5 m^2/s`, so for the ferry's 12 mm
// plating the limit is `0.012^2 / (2 * 1.545e-5) = 4.66 s`. Seconds, not
// milliseconds. `explicitLimit()` returns exactly that -- to machine precision,
// because for a rectangular element the largest generalised eigenvalue of
// `(K, C_lumped)` is exactly the jump through the thinnest direction, which
// `tests/test_thermal.cpp` asserts as a closed form rather than a tolerance.
//
// Milliseconds would need 0.3 mm elements. **Steel is a poor conductor by the
// standards of this estimate**: it is 1/25 of copper's diffusivity, and a 12 mm
// plate equilibrates through its own thickness in about nine seconds.
//
// So explicit conduction is *viable* on the unrefined structural mesh, and it is
// worth saying so plainly rather than repeating a figure that makes the decision
// look easier than it is. Three things make implicit the right choice anyway,
// and the first two are measured in the test file:
//
//   * **The limit falls as `h^2` and the coupling will refine.** Four elements
//     through the same 12 mm plate -- which is what a through-thickness
//     temperature gradient needs to be resolved, and a gradient through the
//     thickness is what bows a plate -- takes it to 0.29 s. A 1.5 mm surface
//     layer resolving a fire's thermal boundary layer takes it to 0.073 s.
//     Implicit pays nothing for either.
//   * **The step should be set by the fire, not by the mesh.** A compartment fire
//     is a twenty-minute event whose boundary condition is smooth over tens of
//     seconds. Backward Euler runs it in as many steps as the *physics* needs:
//     the measured error at 80x the explicit limit on the semi-infinite problem
//     -- twelve steps over a boundary condition that is a step discontinuity --
//     is 10.0 K in an 880 K jump, 1.1%. Backward Euler is first order and that is
//     what first order costs; it is stable at any step, which is the point.
//   * **The factorisation is reused.** At a fixed step `C/dt + K` is factored
//     once and every later step is two triangular solves -- `n b` rather than
//     `n b^2`. That is what makes a long implicit step cheap enough that trading
//     it against a hundred explicit ones is a win rather than a wash.
//
// The `h^2` in the explicit limit is the whole argument, and it is an argument
// about where this is *going*, not about where the mesh is today.
//
// --- Temperature is in KELVIN --------------------------------------------------
//
// Every temperature crossing this interface is absolute, per CLAUDE.md's "SI
// everywhere". That is not a free choice dressed up as a convention: Phase 4's
// next item after this one is radiation, `sigma (T_gas^4 - T_steel^4)`, which is
// meaningless in Celsius, and a module that took Celsius here and Kelvin there
// would produce a plausible wrong number in exactly the way this repo keeps
// finding. The material curves below are *published* in Celsius and convert on
// the way in, in one place, which is `kCelsius`.
//
// --- What this reuses, and one thing it does not -------------------------------
//
// The mesh, the banded solver and the Cuthill-McKee ordering are
// `solid_shell.hpp`'s and `reduction.hpp`'s. `C/dt + K` is symmetric positive
// definite for any positive step -- C is a mass matrix and K is positive
// semi-definite -- so `solidshell::BandedSpd` is exactly the right factorisation
// and a second solver would be a second place to be wrong.
//
// **`solidshell::RestForms` is not reusable here, and it is worth saying why.**
// It looks like it should be: it caches the Gauss weights and the Jacobians a
// conduction operator wants. The weights it does cache are the same object, and
// `gaussVolumes()` is the check that this file's are right. But its `b` is a
// 6 x 24 *strain*-displacement matrix for a vector field, and every one of the
// assumed-strain cures that make the solid-shell an element rather than a locked
// hex has already been applied to it: row 2 is Betsch-Stein sampled at the
// in-plane corners, rows 4 and 5 are Dvorkin-Bathe sampled at the mid-edges.
// Reading `grad N_a` out of it would be reading a matrix for something it does
// not mean. And it should not be cured anyway -- locking is a property of a
// constrained vector field, and the scalar Laplacian has no shear to lock.
// `Forms` below is therefore its own object: the plain Cartesian shape gradient,
// the shape functions, and the same `det J`.
//
// Body frame and SI units per CLAUDE.md.
#pragma once

#include "scantlings.hpp"    // StructuralMaterial
#include "solid_shell.hpp"   // HexMesh, BandedSpd, kNodes/kGauss/kDof

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace sim::thermal {

using solidshell::kDof;
using solidshell::kGauss;
using solidshell::kNodes;

// 0 degrees Celsius, for the material curves and for a caller writing a fire
// temperature down the way a fire test reports it.
inline constexpr double kCelsius = 273.15;

// --- Carbon steel, EN 1993-1-2:2005 --------------------------------------------
//
// The one published, dimensionally complete thermal description of structural
// steel over a fire's range, and the one every structural fire calculation in
// Europe is done against. Sources, by clause:
//
//   * conductivity, 3.4.1.3: `lambda = 54 - 3.33e-2 theta` W/(m K) over
//     20 <= theta < 800 C, and a flat 27.3 above it to 1200 C.
//   * specific heat, 3.4.1.2: a cubic to 600 C, then two hyperbolae, then a flat
//     650 J/(kg K) from 900 C to 1200 C.
//   * density, 3.4.1.4: 7850 kg/m^3, taken as independent of temperature. It is
//     not quite -- steel expands about 1.4% by volume at 700 C -- but the code
//     says to treat it as constant and mixing a constant density with a
//     temperature-dependent everything else is the kind of half-applied model
//     that reads as a physical effect.
//
// **Both matter over a fire's range and the second matters enormously.**
// Conductivity falls 36%, from 53.3 to 34.0 W/(m K), between 20 C and 600 C:
// steel gets *worse* at spreading heat exactly as it starts to lose strength.
// Specific heat is the larger effect and it is not monotone -- the ferrite to
// austenite transition puts a spike at 735 C where `c` reaches
// **5000 J/(kg K), eleven times its room-temperature value**. That spike is a
// latent heat wearing a heat capacity's clothes, and a solve that misses it
// gallops a plate through 735 C in one step that should have taken eleven.
//
// So it is modelled, and `Problem::temperatureDependent` switches it on. What
// makes it affordable, and what makes the energy account still close to machine
// precision across the spike, is that the capacity is formed from a **secant**
// and not a tangent -- see `Solver` below.
//
// Outside 20..1200 C the standard says nothing and these clamp to the end
// values, which keeps `c` integrable and `k` positive rather than extrapolating a
// cubic into a region where it turns over.
double carbonSteelConductivity(double kelvin);  // W/(m K)
double carbonSteelSpecificHeat(double kelvin);  // J/(kg K)

// Specific enthalpy above 20 C: the exact integral of `carbonSteelSpecificHeat`,
// piece by piece, with the constants of integration chosen to make it continuous.
//
// This exists because **`c` is not what conserves energy once `c` depends on
// temperature** -- `h` is. A step from `T_a` to `T_b` moves `rho (h(T_b) -
// h(T_a))` joules per cubic metre, and that is only `rho c dT` when `c` is
// constant over the interval. Across the 735 C spike it is not: a 30 K step
// through it carries about six times what the endpoint `c` values would suggest.
double carbonSteelEnthalpy(double kelvin);  // J/kg

// --- One element ---------------------------------------------------------------
//
// The conduction analogue of `solidshell::RestForms`: everything the element
// derives from its geometry, which nothing in a conduction solve ever moves.
// Same eight-node hexahedron, same 2x2x2 rule, same node ordering contract --
// though unlike the mechanical element the ordering is *not* load-bearing here,
// because a scalar Laplacian has no thickness direction to get wrong. A rotated
// element gives a bit-identical answer up to the summation order, and
// `tests/test_thermal.cpp` asserts it.
struct Forms {
    double gradient[kGauss][kNodes][3];  // dN_a/dx_i, Cartesian
    double shape[kGauss][kNodes];        // N_a
    double weight[kGauss];               // 2x2x2 Gauss weight (one) times det J
    bool ok = false;                     // false on an inverted or degenerate element
};

// Fill `out` from an element's node positions. False -- and `out.ok` false -- when
// `det J` is not positive at one of the eight integration points, which is the
// same condition `solidshell::computeRestForms` refuses on.
//
// A **collapsed** hexahedron -- the wedge a degenerate plate panel extrudes to,
// 166 of them on the reference ferry -- is accepted, for the reason
// `solid_shell.hpp` gives at `ElementShape`: its `det J` is exactly zero at the
// coincident corners and only there, and positive at every point the element is
// integrated at.
bool computeForms(const double nodes[kDof], Forms& out);

// Conductance, 8x8 row-major: `integral grad N_a . k grad N_b dV`, with `k` given
// per Gauss point so a temperature-dependent conductivity costs nothing extra.
// Symmetric, positive semi-definite, and **exactly singular on the constant
// field**: `K 1 = 0`, because `sum_a grad N_a = 0` is a partition of unity. That
// identity is not decoration -- it is what makes the energy account below close,
// and `tests/test_thermal.cpp` asserts it directly.
void conductance(const Forms& forms, const double conductivity[kGauss],
                 double out[kNodes * kNodes]);

// Capacity, 8x8 row-major: `integral rho c N_a N_b dV`, with the *volumetric*
// heat capacity `rho c` in J/(m^3 K) given per Gauss point.
//
// `lumped` replaces it by the diagonal of its own **row sums**, which is the same
// lumping `solidshell::elementMass` uses and for the same reason: volume/8 is
// exact only on a parallelepiped. Row-sum lumping leaves `1^T C dT` **exactly
// unchanged** for any `dT`, by the symmetry of C -- so the two capacities differ
// in where the heat sits and never in how much there is, and the energy account
// is identical under both. That is asserted rather than argued.
//
// Consistent is the default. Lumped exists because a consistent capacity gives
// backward Euler a physically impossible *undershoot* ahead of a steep front
// unless `dt >= rho c h^2 / (6 k)`, and a fire's surface is exactly a steep
// front. **Measured**, on 5 mm elements whose criterion is 0.333 s: at a tenth
// of it the interior dips **14.6 K below its own initial temperature**, and at
// three times it the dip is 7e-13 K. Lumped does not dip at either step -- it is
// unconditionally monotone -- and, by the row-sum identity above, moves exactly
// the same energy. So the choice costs nothing in the account and only decides
// whether a caller who has to take short steps sees a wiggle.
void capacity(const Forms& forms, const double volumetricCapacity[kGauss], bool lumped,
              double out[kNodes * kNodes]);

// Largest stable *explicit* step for one element: `2 / lambda_max`, with
// `lambda_max` the largest eigenvalue of `K x = lambda C x` on the lumped
// capacity, from `reduction::generalisedEigen`.
//
// Exposed because the case for solving this implicitly at all is a number, and a
// number should be measured on the element rather than estimated from a nominal
// thickness -- the same argument `solidshell::criticalTimestep` makes, which got
// the mechanical estimate wrong for the same reason. Here it corrects a figure by
// a factor of a thousand: see the header note above.
//
// For a rectangular element this comes out at exactly `h_min^2 / (2 alpha)` --
// the thin-direction jump is the highest mode, and the full 3D checkerboard is
// three times *lower*. The test file asserts that identity, which is also what
// caught a power iteration whose start vector was orthogonal to the mode it was
// looking for.
double explicitLimit(const double nodes[kDof], const StructuralMaterial& material);
double explicitLimit(const solidshell::HexMesh& mesh, const StructuralMaterial& material);

// --- The mesh boundary ---------------------------------------------------------

// A face of the mesh carried by exactly one element. Faces are numbered in
// `solid_shell.hpp`'s node ordering: 0 is the zeta = -1 face, 1 the zeta = +1
// face, then the four sides.
//
// This is here rather than in a test because a boundary condition has to name a
// surface, and naming it by element-and-face-index is how a caller ends up
// applying a fire to the inside of a bulkhead. A caller filters these by
// `centroid` and `normal` -- geometry it can check -- and hands the survivors
// back.
//
// `normal` is outward: the winding of the face table is not consistently
// oriented, so it is fixed here by comparing the face centroid against the
// element centroid rather than by trusting the table.
struct BoundaryFace {
    std::uint32_t element = 0;
    std::uint8_t  face = 0;
    std::uint32_t node[4]{};
    Vec3   centroid{};
    Vec3   normal{};   // unit, outward
    double area = 0;   // m^2, by the same 2x2 rule the surface integral uses
};

// Every exterior face, in element then face order. A face whose area is zero --
// the degenerate side of a collapsed hexahedron -- is dropped, because it carries
// no heat and its normal is not defined.
std::vector<BoundaryFace> boundaryFaces(const solidshell::HexMesh& mesh);

// A boundary condition of the second and third kinds over a set of faces:
//
//     q_in = flux + coefficient * (ambient - T)      W/m^2, positive INTO the solid
//
// The two are one struct because they are one surface integral and separating
// them would be two places to get `dA` wrong. A pure flux sets `coefficient` to
// zero; a pure convective surface sets `flux` to zero. Radiation is *not* here --
// it is a separate roadmap item and it is nonlinear in a way conduction is not.
//
// Sign: `flux` positive heats the solid, whichever way the face points. That is
// deliberate. The alternative -- resolving the caller's flux against the outward
// normal -- would make the answer depend on a winding the caller cannot see, and
// `normal` is published above for a caller who does want to resolve one.
struct Film {
    std::vector<BoundaryFace> face;
    double flux = 0.0;         // W/m^2 into the solid
    double coefficient = 0.0;  // W/(m^2 K)
    double ambient = 293.15;   // K
};

// --- The problem ---------------------------------------------------------------

struct Problem {
    // Borrowed, not owned, and it must outlive the `Solver` built from it. The
    // solve never moves a node, so a mesh that deforms elsewhere invalidates the
    // factorisation rather than being tracked.
    const solidshell::HexMesh* mesh = nullptr;
    StructuralMaterial material{};

    // Take `k` and `c` from the EN 1993-1-2 curves at the local temperature
    // instead of from `material`. `material.density` is used either way.
    bool temperatureDependent = false;
    bool lumpedCapacity = false;

    // Per node. `prescribed` non-zero holds that node at `prescribedValue`, in
    // kelvin, exactly -- moved to the right-hand side the way
    // `solidshell::solveStatic` moves a prescribed displacement, not by a penalty.
    // Empty means nothing is prescribed, which is legal for a transient and
    // singular for a steady solve unless a film holds it.
    std::vector<std::uint8_t> prescribed;
    std::vector<double>       prescribedValue;

    // W/m^3, one per element. Empty means none. This is where a coupled fire
    // would *not* put its heat -- a fire heats a surface -- but it is what a
    // manufactured-solution test needs and it is one line.
    std::vector<double> volumetricSource;

    std::vector<Film> film;
};

// --- The energy account ---------------------------------------------------------
//
// The analogue of `zone::SolveResult`'s work-against-strain-energy balance, and
// the check most likely to catch a wrong quadrature weight. Every figure is in
// joules and accumulated from the start of the run.
//
// With constant properties it closes to **machine precision** rather than to the
// integrator's order -- measured at 5e-15 of the enthalpy moved over sixty steps
// of three different sizes -- and that is a property of the formulation rather
// than of the step size. Summing the discrete system with `1^T`:
//
//     1^T C (T1 - T0)/dt  +  1^T K T1  =  1^T f + 1^T r
//
// and `1^T K = 0` exactly, so the enthalpy rate is the sum of the applied and
// reaction heat with nothing left over -- at any `dt`, converged or not. The one
// thing that can break it is `C` being wrong, which is exactly what the check is
// for.
//
// **With temperature-dependent properties it closes to the Picard tolerance
// instead**, and the difference is worth being exact about rather than rounding
// up to "machine precision". The secant capacity below makes `1^T C dT`
// identically the enthalpy change *of the state the system was solved for*; an
// iteration stopped at 1e-8 K has not quite reached that state, and the leftover
// is the whole of the residual. Measured: 1.4e-12 of the heat supplied at a 1e-6 K
// tolerance and 1.6e-15 at 1e-11 K, which is the scaling a wrong quadrature would
// not have.
struct Account {
    double enthalpy = 0;        // J, above 20 C, in the mesh now
    double enthalpyChange = 0;  // J, since the solver was prepared
    double prescribedHeat = 0;  // J in through nodes held at a temperature
    double filmHeat = 0;        // J in through flux and convective surfaces
    double sourceHeat = 0;      // J in from volumetric sources

    // enthalpyChange - (prescribedHeat + filmHeat + sourceHeat). Joules, and it
    // should be at the rounding of the largest term.
    double residual() const {
        return enthalpyChange - prescribedHeat - filmHeat - sourceHeat;
    }

    // W. The largest row of `A T1 - b` over the *free* nodes after the last solve
    // -- which is zero if and only if the band actually held every term that was
    // scattered into it.
    //
    // It is reported rather than assumed because `BandedSpd::add` **silently
    // drops** an entry outside its band, so a bandwidth computed one element short
    // reads as a slightly soft answer and not as an error. That is the shape of
    // the defect `CLAUDE.md` records against `reduction`'s node ordering, arriving
    // by a different door.
    double equilibriumResidual = 0;
};

// --- The solver ------------------------------------------------------------------
//
//     (C/dt + K) T1 = C/dt T0 + f
//
// Backward Euler: unconditionally stable, unconditionally monotone with a lumped
// capacity, and first order in time. First order is a real cost and it is
// deliberate -- Crank-Nicolson is second order and oscillates on exactly the step
// change in surface temperature a fire is, and the roadmap's milestone is a
// bulkhead under a fire and not a smooth ramp. The order is *asserted* in the
// tests rather than assumed, because a scheme that quietly came out first order
// when it was meant to be second is a defect a tolerance would hide.
//
// --- The factorisation is reused ---------------------------------------------
//
// `C/dt + K` does not depend on the state when the properties do not, so with a
// fixed step it is factored **once** and every subsequent step is a pair of
// triangular solves. `factorisations()` reports how many times it was actually
// factored, so a caller can tell whether it is paying `n b^2` per step or `n b`.
// One scalar unknown per node against the mechanical solver's three makes the
// band a third as wide and the factorisation a ninth of the cost.
//
// --- Temperature-dependent properties, and the secant ---------------------------
//
// With `Problem::temperatureDependent` the system is nonlinear and is closed by
// Picard iteration: re-evaluate `k` and `c` at the current iterate, re-assemble,
// re-solve, until the temperature stops moving. `iterations()` reports the last
// step's count.
//
// The capacity used is the **secant**
//
//     c = (h(T1) - h(T0)) / (T1 - T0)
//
// at each Gauss point, and not `c(T1)` or `c(T0)` or `c` at the midpoint. Two
// things follow, and they are the reason for it:
//
//   * `1^T C (T1 - T0)` is then *identically* the enthalpy change, so the energy
//     account above still closes -- to the Picard tolerance, which is the only
//     thing left between the accepted state and the one the last system was solved
//     for. A tangent `c` leaves a residual of order `dt` per step instead, which
//     accumulates and which no tolerance on the iteration would shrink.
//   * The 735 C spike is integrated rather than sampled. A tangent evaluated at
//     either end of a step that crosses the spike sees `c ~ 700` where the true
//     average over the interval is several thousand, so the plate walks through
//     the phase change at six times the rate it should. The secant cannot miss it,
//     because it *is* the average.
//
// A tangent would also be the wrong object for the same reason `plasticity.hpp`'s
// softening wanted a secant: the quantity being conserved is an integral, and the
// derivative at a point is not it.
class Solver {
public:
    // Build the forms, choose a node numbering and take the initial state.
    // Returns false with a reason on an inverted element, a malformed problem, or
    // a mesh that has no nodes.
    bool prepare(const Problem& problem, double uniformTemperature, std::string* why = nullptr);
    bool prepare(const Problem& problem, const std::vector<double>& initialTemperature,
                 std::string* why = nullptr);

    // Advance one step. False on a step that is not positive, a system that will
    // not factor, or a Picard iteration that did not converge -- the state is
    // still advanced in the last case, the same choice
    // `solidshell::elementPlasticUpdate` makes, because a caller that has to keep
    // stepping is better served by a slightly wrong state it is told about.
    bool step(double timestep, std::string* why = nullptr);

    // Drop the capacity and solve `K T = f` directly. Singular, and refused, when
    // nothing holds the temperature -- no prescribed node and no convective film.
    bool solveSteady(std::string* why = nullptr);

    const std::vector<double>& temperature() const { return temperature_; }
    const Account& account() const { return account_; }

    std::size_t halfBandwidth() const { return band_; }
    std::size_t freeNodes() const { return free_; }
    int factorisations() const { return factorisations_; }
    int iterations() const { return iterations_; }
    double time() const { return time_; }

    // Picard controls. `tolerance` is in kelvin on the largest nodal change
    // between iterations.
    void setPicard(double tolerance, int maximumIterations) {
        picardTolerance_ = tolerance;
        picardLimit_ = maximumIterations;
    }

private:
    // One face of the boundary a `Film` acts on, with its surface integrals taken
    // once. `mass` is `integral N_a N_b dA` and `load` is `integral N_a dA`; both
    // are geometry, so a convective boundary costs one 4x4 scatter per step and no
    // quadrature at all.
    struct FilmFace {
        std::uint32_t node[4]{};
        std::uint32_t film = 0;
        double mass[16]{};
        double load[4]{};
    };

    void number();
    void refreshProperties(const std::vector<double>& evaluateAt,
                           const std::vector<double>& previous, bool secant);
    bool buildAndSolve(double inverseStep, bool useCapacity, bool refactor, std::string* why);
    void residual(double inverseStep, bool useCapacity, std::vector<double>& out) const;
    double enthalpyOf(const std::vector<double>& field) const;
    double specificEnthalpy(double kelvin) const;

    Problem problem_{};
    std::vector<Forms> forms_;
    std::vector<double> conductance_;  // 64 per element
    std::vector<double> capacity_;     // 64 per element
    std::vector<FilmFace> filmFace_;
    std::vector<double> load_;         // full length: volumetric source + film supply
    std::vector<double> temperature_;
    std::vector<double> previous_;
    std::vector<std::ptrdiff_t> map_;   // node -> free slot, or -1
    std::vector<double> rhs_;
    std::vector<double> work_;
    solidshell::BandedSpd system_{0, 0};
    Account account_{};

    std::size_t nodes_ = 0, elements_ = 0, free_ = 0, band_ = 0;
    double time_ = 0, factoredStep_ = 0, enthalpyStart_ = 0;
    double picardTolerance_ = 1e-8;
    int picardLimit_ = 40, factorisations_ = 0, iterations_ = 0;
    bool factored_ = false, propertiesFresh_ = false, ready_ = false;
};

}  // namespace sim::thermal
