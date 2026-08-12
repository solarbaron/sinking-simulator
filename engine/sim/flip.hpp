// SPDX-License-Identifier: MIT
//
// A **sparse FLIP/APIC solver for interior water** -- `docs/06-roadmap.md` Phase 5's
// first item, and the fidelity tier above the well-mixed water volume
// `sim::Ship`'s compartments carry today.
//
// A flooded compartment in `ship.hpp` is one number: `Compartment::waterVolume`,
// with a flat free surface at whatever height that volume implies. That is the
// right first answer, it is what the stability solve wants, and it costs nothing.
// What it cannot do is *move*: it has no velocity field, so it cannot slosh, it
// cannot pile up against a bulkhead in a roll, it cannot run through a door as a
// jet, and it cannot tell you what any of that does to the ship. Each of those is
// a field. This file gives one region of water fields, on a sparse staggered grid
// with particles, and hands the answer back as a mass and a level when it is done.
//
// **This file is the solver and nothing else.** It does not render, it is not
// wired into `sim::Ship`'s flooding model, and it does not contain the
// quiescent <-> dynamic escalation that is Phase 5's *second* item. What it does
// carry, deliberately, is everything that escalation needs to be exact -- see §5.
//
// --- 1. Which method, and the measurement that chose it ------------------------
//
// PIC, FLIP and APIC differ only in what a particle carries across the grid and
// back. PIC replaces the particle's velocity with the interpolated grid velocity;
// it is unconditionally stable and it is a low-pass filter applied once per step,
// so a sloshing tank damps out. FLIP adds the interpolated *change* to the
// particle's own velocity, so it dissipates almost nothing -- and keeps a
// per-particle velocity component the grid never sees, which shows up as noise.
// APIC (Jiang et al. 2015) carries an affine velocity field per particle, `C_p`,
// which is exactly the information a single velocity loses; it is filtered like
// PIC and conserves angular momentum like FLIP.
//
// **This solver is APIC**, `Params::flipBlend = 0` and `Params::affine = true`, and
// the other two are reachable from the same code path so that the choice can be
// measured rather than cited. What was measured, in `tests/test_flip.cpp` and
// printed in full by `tools/flip_probe`:
//
//   * **Angular momentum through the grid round trip**, on a 0.8 m cube of water
//     in solid-body rotation at 2 rad/s on a 0.05 m grid, transfer only. PIC keeps
//     **98.65%** of it after one round trip and **89.23%** after ten -- 1.08% per
//     trip. APIC keeps **99.29%** and **98.24%** -- 0.18% per trip, **6.1 times
//     less**. The residue is not the method: it is the free surface, where the
//     transfer discards the momentum it deposits on the air side (§3) and the
//     extrapolation puts back the fluid's value instead. FLIP loses exactly none,
//     because with no grid dynamics its round trip is the identity -- which is
//     also the control that caught a real defect in the ordering of `saveGrid`.
//   * **Particle noise**, the RMS of `v_p - interp(v_grid)` as a fraction of the
//     RMS particle speed, after five seconds of a sloshing tank. FLIP reaches
//     **87.6%**; APIC holds **1.31%**; PIC is 0.00% by construction, since it *is*
//     the interpolated field. Sixty-seven times less noise than FLIP.
//   * **The sloshing period itself**, same tank, against `g k tanh(k d)`. APIC
//     comes out **+4.45%**, PIC **+9.16%** and FLIP **+9.81%** -- and PIC crosses
//     the centreline three times in five seconds where APIC crosses eight, because
//     PIC damps the mode it is meant to be measuring.
//
// So APIC is twice as accurate as either on the one closed form available, at a
// sixty-seventh of FLIP's noise, which is why it is the default and why the blend
// is left at zero rather than at the 0.95 a pure-FLIP solver would need.
//
// The kernel is the **quadratic B-spline**, and that is not a free choice either.
// APIC needs the second moment of its own kernel,
//
//     D_p = sum_n w_np (x_n - x_p) (x_n - x_p)^T,
//
// and only for a quadratic (or higher) B-spline is that a *constant*: `h^2/4 I`,
// independent of where the particle sits. With trilinear weights `D_p` is
// diagonal but position dependent -- `h^2 f(1-f)` per axis -- and therefore
// **singular** whenever a particle lands on a grid node, which is the one place a
// solver must not have a special case. The constant is asserted directly
// (`tests/test_flip.cpp`, over 200 000 offsets spanning forty cells, worst
// deviation 2.3e-16 of `h^2/4`, alongside `sum w = 1` to 2.3e-16 and the first
// moment to 6.2e-17) rather than taken on trust, because it is the whole reason
// `C_p` can be formed by a multiply instead of a matrix inverse.
//
// --- 2. Sparse, and what "sparse" has to mean to be worth it -------------------
//
// The grid is stored as **tiles of 4x4x4 cells in a hash map**, and a tile exists
// only when some cell within `kHalo = 3` of a cell containing particles falls
// inside it. Three consequences, each of which is a test rather than a claim:
//
//   1. **An empty compartment costs nothing.** No particles, no tiles, no bytes,
//      no work: `Solver::bytes()` is 0 and `step()` returns having done nothing,
//      whatever the domain is. A 400x400x400-cell compartment -- 64 million cells,
//      7.4 GB dense at the 115 bytes a cell this structure costs -- with no water
//      in it allocates zero.
//   2. **Cost follows the water, not the room.** The same body of water in a 20^3
//      domain and in a 400^3 domain allocates the *same* tiles and produces
//      **bit-identical** particle positions after fifty steps. That is a stronger
//      statement than "about the same cost": it says the domain extent does not
//      enter the arithmetic at all. And the tile count is checked against one
//      derived by hand from the water's own extent rather than against another run
//      -- `tiles()` shipped returning `tileKey_.size()`, three coordinates per
//      tile, and every count-against-count assertion in the suite was happy.
//   3. **Arrival is the hard part and it is where sparse solvers fail.** Water
//      moving into a region that had no storage a moment ago must find storage at
//      the moment it needs it -- and the moment it needs it is one transfer
//      *before* a particle enters the region, because the kernel reaches two cells
//      ahead of the particle. The invariant that catches a missing tile is
//      exact and cheap: the face weights of each component sum to the total
//      particle mass, since the kernel is a partition of unity. A halo one tile
//      too small drops mass out of that sum and nothing else notices.
//
// The halo is 3 cells rather than 2 because the transfer reaches `+2` along a
// component's own axis (§1), the velocity extrapolation that fills the air side of
// the surface needs one layer beyond that, and the RK2 advection samples the field
// at a midpoint up to a Courant number away. Two would work for the transfer alone
// and would be wrong for the advection, silently, only for fast water.
//
// --- 3. What the projection is, and what it shares with `les.cpp` --------------
//
// `les.cpp` -- the low-Mach compartment fire that landed a week before this -- is
// the closest prior art in the repo and the conventions below are deliberately
// its:
//
//   * A **MAC/staggered grid**: pressure at cell centres, velocity on faces.
//     Collocated is cheaper to index and has a checkerboard null space, which is
//     exactly the shape the pressure field would put into it.
//   * The **all-Neumann compatibility discipline**. A sealed box determines its
//     pressure only up to a constant and only accepts a right-hand side that sums
//     to zero; `les.cpp` gets that by taking the heating's deviation from its own
//     mean, and this file gets it by removing the mean of the residual on a fully
//     submerged region. A tank with no free surface in it is not an edge case
//     here, it is *the* case for a flooded compartment, so it is solved rather
//     than excluded.
//   * **Compensated summation on every conserved total**, for the reason
//     `les.cpp`'s `Accumulator` gives: a test cannot be tighter than its own
//     arithmetic, and every assertion this file is under is a conservation
//     assertion.
//   * A **substep controller that publishes when its budget ran out** rather than
//     silently under-advancing. `les::StepResult::incomplete` exists because a
//     mutation drove its rejection test into a corner and the suite *hung*; the
//     same shape and the same publication are here, and `Params::maxSubsteps` and
//     `Params::projectionIterations` are hard bounds, so no input can make a step
//     take unbounded time. That is what lets a mutation sweep bound a mutant by
//     arithmetic rather than by a wall clock.
//
// What is deliberately **not** shared:
//
//   * **The linear solver.** `les.cpp` uses red-black SOR, because its answer must
//     commute with a reflection of the grid and a colour sweep is order
//     independent within a colour. This file uses **Jacobi-preconditioned
//     conjugate gradients** instead, because the sharpest test of a water
//     projection is that a hydrostatic column does not move *at all*, and that is
//     a statement about the residual reaching machine precision rather than about
//     symmetry. Measured on that column's own Poisson problem, written out
//     standalone in `tools/flip_probe --solver`: red-black SOR at `les.cpp`'s own
//     relaxation of 1.7 needs **1 170 sweeps** to reach a 1e-13 residual on 24
//     cells and **does not reach 1e-15 in 200 000**; Jacobi-preconditioned CG
//     reaches exactly zero in **24 iterations**. At 64 cells it is 8 221 sweeps
//     against 64 iterations. Conjugate gradients cost a global dot product, so the
//     answer is no longer bit-exact under a reflection; this file therefore makes
//     no mirror-symmetry claim, and `les.cpp` keeps SOR for the claim it does make.
//   * **The variable being transported.** `les.cpp` advects mass on the grid and
//     derives temperature; here the mass is *on the particles* and the grid never
//     owns it. That is why mass conservation here is exact in a stronger sense
//     than a flux-form scheme can be -- see §4.
//   * **Anything thermodynamic.** Water is incompressible and at one density.
//
// --- 4. Mass, and why "exact" is the right word here ---------------------------
//
// A finite-volume scheme conserves mass to *summation error*, because every face
// flux is added to one cell and subtracted from another. A particle scheme does
// better: a particle's mass is never modified after it is created, so the total is
// the same sum of the same numbers at every step, and `Account::massResidual()` is
// **exactly 0.0**, not 1e-16. The tests assert `== 0.0` and the harness's
// `expectNear` is given a tolerance of zero, because anything else would be a
// tolerance on an identity.
//
// There are exactly two ways to lose mass here and both are closed:
//
//   * A particle leaving the domain. Advection clamps a particle back inside the
//     solid box and zeroes the velocity component that drove it out;
//     `Account::clamped` counts every time that happens, so "no mass left" is a
//     claim with a witness rather than an absence of evidence.
//   * A particle being deleted. Nothing here deletes particles.
//
// `Account::particles` is therefore constant, and it is asserted as an integer
// rather than inferred from the mass, because two compensating errors in a mass
// sum are a thing that happens and two compensating errors in a count are not.
//
// --- 5. What Phase 5's *second* item will need, and what is here for it --------
//
// "Quiescent <-> dynamic escalation with exact mass conservation" is the same
// problem `promotion.hpp` solves for structure and `les::promote`/`les::demote`
// solve for gas: state has to cross a fidelity boundary without a tuning constant
// and without losing the conserved quantity. That machinery is *not* built here,
// on purpose -- it needs `sim::Ship`, and wiring it now would make this
// untestable. What is here is the pair of exact statements it will stand on:
//
//   * **Down**, `seedBox` fills a box with a lattice of particles and
//     `setTotalMass` gives them a total that is *exactly* the mass asked for. The
//     last particle takes the remainder, so `N` equal shares that do not sum to
//     the whole in floating point still sum to the whole. That is the same trick
//     `fire.cpp`'s `layerSplit` and `les.cpp`'s two-layer reduction use, and it is
//     the reason the round trip below can be an identity rather than a tolerance.
//   * **Up**, `quiescentLevel` returns the still-water level a body of water
//     implies over a given plan area: `z_floor + M / (rho A)`. One closed form, no
//     iteration, no fitted constant, and it is the *same* number
//     `Compartment::waterVolume` already means. Seeding 2.7183 m and reading it
//     back returns **exactly 2.7183 m** -- 0.0 m of error, not a rounding of one --
//     because the mass that went in is bit for bit the mass `M/(rho A)` divides.
//
// Neither of those is escalation. Both are the arithmetic escalation has to be
// built on, and having them here means the second item is a wiring job rather than
// a redesign.
//
// --- 6. What is not here, named rather than hidden -----------------------------
//
//  1. **Solids are one axis-aligned box.** A compartment with a stair, a girder or
//     a heeled deck in it needs a solid *mask* and a fractional face area. The face
//     classification below is already per face, so that is an addition rather than
//     a rewrite, but it is not written and a curved boundary would be wrong.
//  2. **The free surface is voxelised.** An air cell is at `p = 0` at its centre,
//     so the surface is resolved to the cell and the pressure is first order there.
//     A ghost-fluid free surface would fix it and would also cost the hydrostatic
//     column its exactness, which is the test this file is built around.
//  3. **No volume correction.** FLIP-family solvers drift in *volume* even while
//     conserving mass exactly, because nothing stops particles clumping. There is
//     no density projection, no resampling and no particle merge/split here. §4's
//     mass claim is a claim about mass; the volume claim is only that the
//     projection keeps the grid divergence-free, and that is asserted separately.
//  4. **No viscosity, no surface tension, no two-phase air.** The air is a void at
//     zero pressure.
//  5. **No moving solids.** A rolling ship moves the compartment, which means a
//     time-dependent frame acceleration and a moving boundary. `Params::gravity`
//     is a vector so the first half is expressible; the second is not.
//  6. **Single threaded.** The job system is not used. A sparse tile grid is the
//     right shape for it and nothing here forecloses it, but a solver whose
//     correctness is still being established should not be establishing its
//     thread-safety at the same time.
//
// SI units, body frame (+x bow, +y port, +z up), per CLAUDE.md. Gravity is
// `-z` by default and lives in `Params` rather than in the code.
#pragma once

#include "../core/math.hpp"

#include <cstddef>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

namespace sim::flip {

// ---------------------------------------------------------------------------
// The grid
// ---------------------------------------------------------------------------

// Cells per tile edge. Four is small enough that a puddle costs a few tiles and
// large enough that the hash is not the whole cost; the halo in §2 is 3 cells, so
// a larger tile would not reduce the number of *neighbour* tiles a fluid tile
// pulls in, only the number of tiles.
inline constexpr int kTile = 4;
inline constexpr int kTileCells = kTile * kTile * kTile;

// Cells of storage kept around the water. See §2 for why it is 3 and not 2.
inline constexpr int kHalo = 3;

// A uniform cubic grid. Cubic rather than the three independent spacings
// `les::Grid` carries, because APIC's `D_p = h^2/4 I` is a scalar only when the
// cell is a cube -- an anisotropic cell makes it a diagonal matrix and every
// affine term a per-axis one, which is a real generalisation and not this one.
//
// The domain is `n[0] x n[1] x n[2]` cells of solid-bounded water. Everything
// outside is solid; the sparse structure is free to allocate tiles out there,
// because the wall face of the outermost cell is stored on the cell beyond it.
struct Grid {
    double h = 0.05;        // m, cell edge
    double lo[3]{0, 0, 0};  // m, body frame, the low corner of cell (0,0,0)
    int    n[3]{0, 0, 0};   // cells along x, y, z

    bool empty() const { return n[0] <= 0 || n[1] <= 0 || n[2] <= 0; }
    double cellVolume() const { return h * h * h; }
    // Cells times cell volume rather than the product of the three side lengths,
    // for the reason `les::Grid::volume` gives: the two differ in the last bits.
    double volume() const {
        return cellVolume() * static_cast<double>(n[0]) * static_cast<double>(n[1]) *
               static_cast<double>(n[2]);
    }
    double planArea() const { return h * h * static_cast<double>(n[0]) * static_cast<double>(n[1]); }
    double hi(int axis) const { return lo[axis] + h * static_cast<double>(n[axis]); }

    // The cell a coordinate falls in, unclamped: outside the domain this is
    // negative or past `n`, which is exactly what the solid test wants.
    int cellOf(int axis, double x) const;
    bool inside(int i, int j, int k) const {
        return i >= 0 && j >= 0 && k >= 0 && i < n[0] && j < n[1] && k < n[2];
    }
    // The centre of cell (i, j, k) along `axis`.
    double centre(int axis, int index) const {
        return lo[axis] + (static_cast<double>(index) + 0.5) * h;
    }
    // The position of the low face of cell `index` along `axis` -- where the
    // staggered velocity component `axis` lives.
    double facePosition(int axis, int index) const {
        return lo[axis] + static_cast<double>(index) * h;
    }
};

// ---------------------------------------------------------------------------
// The particles
// ---------------------------------------------------------------------------

// One material point. `affine` is APIC's `C_p`, row major: `affine[3 r + c]` is
// `d v_r / d x_c`, so the velocity the particle would hand a node at `x` is
// `velocity + C (x - position)`.
//
// Plain data, trivially copyable, never destructed -- the same rule ECS components
// follow, because a particle array is the one thing here that gets large.
struct Particle {
    double position[3]{0, 0, 0};
    double velocity[3]{0, 0, 0};
    double affine[9]{0, 0, 0, 0, 0, 0, 0, 0, 0};
    double mass = 0;   // kg
};

// One region of water: the grid it lives in and the particles that are it.
struct Field {
    Grid grid;
    std::vector<Particle> particles;
    double time = 0;   // s of model time

    std::vector<std::string> problems;

    bool empty() const { return grid.empty() || particles.empty(); }
    // kg, compensated. The sum of the same numbers at every step -- see §4.
    double totalMass() const;
    // m^3 the mass occupies at `density`; the *material* volume, not the volume
    // the particles' bounding cells happen to add up to.
    double volumeAt(double density) const;
};

// ---------------------------------------------------------------------------
// Parameters
// ---------------------------------------------------------------------------

struct Params {
    double density = kRhoSeawater;   // kg/m^3, uniform
    double gravity[3]{0, 0, -kGravity};   // m/s^2, body frame

    // The transfer. `affine` off and `flipBlend` 0 is PIC; off and 1 is FLIP; on
    // and 0 is APIC. See §1 for the measurement behind the default.
    bool   affine = true;
    double flipBlend = 0.0;

    // The advective Courant number a substep is allowed to reach. Unlike a
    // finite-volume advection this is not a stability limit -- the particles carry
    // the transport, so a large step is inaccurate rather than unstable -- but the
    // extrapolated halo is `kHalo` cells deep and a particle must not leave it.
    double courant = 0.9;
    double maxSubstep = 0.005;   // s
    // A hard bound, published when it binds. See §3: the reason this exists at all
    // is that a mutation which collapses a substep controller must *fail* rather
    // than run forever.
    int    maxSubsteps = 4096;

    // The projection. Jacobi-preconditioned CG; the tolerance is relative to the
    // right-hand side with an absolute floor, exactly as `les.cpp`'s is.
    // 1e-13 is where the hydrostatic column bottoms out: measured over a 24-cell
    // column, 1e-12 leaves 3.8e-15 m/s of residual velocity, 1e-13 leaves 5.5e-16,
    // and 1e-15 leaves 5.8e-16 for three more iterations. Tightening past the knee
    // buys nothing and loosening to 1e-8 costs five orders of magnitude.
    int    projectionIterations = 400;
    double projectionTolerance = 1e-13;

    // Layers of velocity extrapolated into the air side of the surface. Must be at
    // least 2 for the transfer alone (§2); 3 covers the RK2 midpoint as well.
    int    extrapolationDepth = 3;

    // How far inside a wall a clamped particle is put, as a fraction of `h`.
    // Non-zero so that a particle on a wall is unambiguously in the cell beside
    // it rather than in the solid cell behind it -- the same "on the boundary"
    // discipline `section.cpp` needed and did not have.
    double wallMargin = 1e-3;
};

// ---------------------------------------------------------------------------
// The account
// ---------------------------------------------------------------------------

// Where the mass went, on the same terms as `les::Account` -- except that here the
// answer is an identity rather than a residual, and §4 says why.
struct Account {
    double    initialMass = 0;
    double    mass = 0;
    long long initialParticles = 0;
    long long particles = 0;
    long long clamped = 0;   // times a particle was pushed back off a wall

    double    massResidual() const { return mass - initialMass; }
    // Normalised by the state, for the reason `les::Account` gives: an empty field
    // would otherwise report its own zero as a whole-number fraction.
    double    massResidualFraction() const;
    long long particleResidual() const { return particles - initialParticles; }
};

struct StepResult {
    double time = 0;         // s of model time after the step
    int    substeps = 0;
    int    projectionIterations = 0;   // on the last substep
    double projectionResidual = 0;     // relative, what the last projection left
    bool   projectionCapped = false;   // it ran out of iterations
    // The substep budget ran out before the whole `dt` was taken, so the field is
    // short of the time it was asked for. Published rather than swallowed, for
    // `les::StepResult::incomplete`'s reason.
    bool   incomplete = false;
    double courant = 0;       // the worst advective Courant number committed
    double peakSpeed = 0;     // m/s
    double maxDivergence = 0; // 1/s, the worst cell after the last projection
    int    tiles = 0;
    int    fluidCells = 0;
    long long clamped = 0;
};

// ---------------------------------------------------------------------------
// Seeding, and the two exact statements Phase 5's second item stands on (§5)
// ---------------------------------------------------------------------------

// Fill the intersection of `[lo, hi)` with the domain with a lattice of
// `perAxis^3` particles per cell, at rest, each carrying `density` times its share
// of the cell volume. Particles are *appended*, so two calls make two bodies.
//
// The lattice is deterministic and centred: the `m`th of `perAxis` along an axis
// sits at `(m + 0.5) / perAxis` of the cell. No jitter -- a jittered seed is
// irreproducible across a promotion, and the escalation in §5 has to be able to
// rebuild the same field twice.
void seedBox(Field& field, const double lo[3], const double hi[3], int perAxis,
             double density);

// Scale every particle's mass so the total is **exactly** `mass`. The last
// particle takes the remainder against a compensated sum of the others, so the
// total is the number asked for and not a rounding of it. Does nothing, and says
// so, if there are no particles.
bool setTotalMass(Field& field, double mass);

// The still-water level a body of water implies over `planArea`, measured from
// `floorZ`: `floorZ + M / (rho A)`. The closed form `Compartment::waterVolume`
// already means; §5.
double quiescentLevel(const Field& field, double density, double planArea, double floorZ);

// Diagnostics that are also the state an escalation would hand across.
double kineticEnergy(const Field& field);            // J
void   linearMomentum(const Field& field, double out[3]);   // kg m/s
void   angularMomentum(const Field& field, const double about[3], double out[3]);
void   centroid(const Field& field, double out[3]);   // m, mass weighted

// Problems with a field's definition, in the spirit of `fire::Model::validate`.
std::vector<std::string> validate(const Field& field, const Params& params = {});

// ---------------------------------------------------------------------------
// The transfer kernel, exposed because its identities are what the tests assert
// ---------------------------------------------------------------------------

// Quadratic B-spline weights and node offsets for a coordinate `s` measured in
// cells from node zero of a component's own lattice. Fills `w[0..2]` with the
// weights of nodes `base`, `base+1`, `base+2` and `offset[0..2]` with
// `node - s` in cells; returns `base`.
//
// The two identities that make this the kernel: `sum w = 1` and
// `sum w offset^2 = 1/4`, both independent of `s`. Neither is an approximation
// and both are asserted directly.
int splineWeights(double s, double w[3], double offset[3]);

// `sum_n w_n (x_n - x_p)^2` in metres squared, for a particle at fractional
// position `s` in cells. Exactly `h^2/4`; exposed so a test can say so.
double kernelSecondMoment(double s, double h);

// ---------------------------------------------------------------------------
// The solver
// ---------------------------------------------------------------------------

enum class Cell : std::uint8_t { Solid = 0, Air = 1, Fluid = 2 };

// A face is Solid when either cell beside it is; Fluid when neither is solid and
// at least one holds particles; Air otherwise. Only Fluid faces are unknowns of
// the projection; Air faces are filled by extrapolation and Solid faces are zero.
enum class Face : std::uint8_t { Solid = 0, Air = 1, Fluid = 2 };

class Solver {
public:
    // Advance `field` by `dt`, subdividing to `Params::courant`. The sparse
    // structure, the grid velocities and the pressure are left in place
    // afterwards, so a caller -- or a test -- can interrogate them.
    StepResult step(Field& field, double dt, const Params& params, Account& account);

    // The pieces, in the order `step` runs them. Public because the sharpest tests
    // here are about one stage in isolation: a transfer that is a partition of
    // unity, a projection that reproduces a discrete Helmholtz decomposition, an
    // extrapolation that fills exactly the faces the transfer will read.
    void rebuild(const Field& field, const Params& params);
    void transferToGrid(const Field& field, const Params& params);
    // Take the current grid as the "before" state the FLIP difference is measured
    // against. It is called **after** the first extrapolation and not at the end of
    // the transfer, and that ordering is load bearing: a face on the air side of
    // the surface is zero when the transfer leaves it and non-zero once it has been
    // extrapolated, so saving too early makes FLIP's difference the whole
    // extrapolated velocity instead of the change in it. That defect does not fail
    // quietly -- it drives a rotating body's kinetic energy up twelvefold in ten
    // transfers -- but it is invisible at `flipBlend = 0`, which is the default.
    void saveGrid();
    void addBodyForce(double dt, const Params& params);
    void project(double dt, const Params& params);
    void extrapolate(const Params& params);
    void transferToParticles(Field& field, const Params& params) const;
    void advect(Field& field, double dt, const Params& params, Account& account) const;

    // --- queries ------------------------------------------------------------
    // Three coordinates per tile, so the count is a third of the vector. Reported
    // as `tileKey_.size()` at first, and nothing caught it: every tile assertion
    // in the suite compared one count against another -- equal in two rooms,
    // changed over a fall -- and a factor of three is invisible to all of those.
    // `tests/test_flip.cpp` now derives the count from the water's own tile
    // bounding box instead, which is the difference between checking a number and
    // checking the number a solver reports.
    int    tiles() const { return static_cast<int>(tileKey_.size() / 3); }
    int    fluidCells() const { return static_cast<int>(fluid_.size()); }
    // Bytes the sparse structure holds. Zero when there is no water, which is §2's
    // first claim and is asserted as `== 0`.
    std::size_t bytes() const;
    bool   allocated(int i, int j, int k) const { return index(i, j, k) >= 0; }
    Cell   cellAt(int i, int j, int k) const;
    Face   faceAt(int axis, int i, int j, int k) const;
    // The staggered velocity on the low face of cell (i, j, k) along `axis`.
    double faceVelocity(int axis, int i, int j, int k) const;
    // The transferred mass on that face: `sum_p w_np m_p`. Summed over every face
    // of one component this is the total particle mass exactly -- §2's arrival
    // invariant.
    double faceMass(int axis, int i, int j, int k) const;
    double totalFaceMass(int axis) const;
    double pressureAt(int i, int j, int k) const;    // Pa
    // The divergence of the *current* face velocities in one cell, 1/s. Asked per
    // cell rather than as a norm, because an error that cancels globally is this
    // repo's characteristic defect shape.
    double divergenceAt(int i, int j, int k) const;
    double maxDivergence() const;
    int    lastIterations() const { return iterations_; }
    double lastResidual() const { return residual_; }
    bool   lastCapped() const { return capped_; }
    // True when the region has no free surface at all, so the pressure is
    // determined only up to a constant. Not an error -- see §3.
    bool   singular() const { return singular_; }
    // Particles found in a cell outside the solid box by the last `rebuild`.
    // Nothing here creates one; a non-zero count is a field that was seeded
    // without being validated, and it is a count rather than a silent drop.
    int    outsideParticles() const { return outside_; }

    // Interpolate the current grid velocity at a point, with the same kernel the
    // transfer uses. What advection samples.
    void   sampleVelocity(const double x[3], double out[3]) const;

    // Write one face velocity, returning false if that face has no storage. Two
    // callers are foreseen and one exists: a test that hands the projection a
    // field with a known decomposition, and -- Phase 5's third item -- a jet
    // entering through a breach, which is a prescribed velocity on a face and
    // nothing else.
    bool   setFaceVelocity(int axis, int i, int j, int k, double value);

private:
    int    index(int i, int j, int k) const;   // flat cell index, or -1
    int    tileIndex(int ti, int tj, int tk) const;
    void   buildNeighbours();
    void   classifyFaces();
    void   buildFluidList();
    double dotFluid(const std::vector<double>& a, const std::vector<double>& b) const;
    // The twenty-seven tiles around one tile, gathered once so the transfer's
    // inner loop is arithmetic rather than hash probes. `void*` because the
    // gathered block is an implementation detail with no business in the header.
    void   gatherHalo(int tile, void* out) const;
    int    haloIndex(const void* halo, int i, int j, int k) const;

    Grid   grid_;
    std::unordered_map<std::uint64_t, int> map_;
    std::vector<std::int32_t> tileKey_;   // 3 per tile

    std::vector<std::uint8_t> cell_;
    std::vector<std::int32_t> nb_;        // 6 per cell: -x +x -y +y -z +z, or -1
    std::vector<double>       vel_[3], old_[3], wgt_[3];
    std::vector<std::uint8_t> face_[3];
    std::vector<double>       pressure_;

    // The projection runs over a compact list of fluid cells rather than over the
    // tiles, so its cost is the water and not the storage around it.
    std::vector<std::int32_t> fluid_;      // flat cell index per fluid cell
    std::vector<std::int32_t> fluidOf_;    // flat cell index -> fluid slot, or -1
    // 6 per fluid cell: the neighbour's fluid slot, -1 for air (Dirichlet zero) or
    // -2 for solid (the term drops out of both the matrix and its diagonal).
    std::vector<std::int32_t> fnb_;
    std::vector<double>       diag_, rhs_, sol_, res_, dir_, mat_, pre_;

    // Particles gathered behind the tile they occupy, so the transfer pays the
    // hash once per tile instead of eighty-one times per particle.
    std::vector<std::int32_t> order_, tileStart_, particleTile_;
    std::vector<std::uint8_t> mark_, valid_, filled_;
    std::vector<double>       scratch_;

    int    iterations_ = 0;
    double residual_ = 0;
    bool   capped_ = false;
    bool   singular_ = false;
    int    outside_ = 0;
};

// Zero the account and take the field as its baseline, on `les::resetAccount`'s
// terms.
void resetAccount(const Field& field, Account& account);

}  // namespace sim::flip
