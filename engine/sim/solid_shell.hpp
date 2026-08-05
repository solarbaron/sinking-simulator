// SPDX-License-Identifier: MIT
//
// Solid-shell elements -- the element technology `docs/07-fem-spike-findings.md`
// §4 specified after the spike ruled out uniform linear tetrahedra for plating.
//
// The spike measured a linear tet cantilever at 63% error with two elements
// through the thickness, 32% at four, 11% at eight. Ship plating is 10-25 mm, so
// resolving bending with tets means 2-3 mm elements, and the explicit stability
// limit is set by the *smallest* element dimension. Cost scales as h^-4. The two
// constraints close on each other and leave no workable resolution.
//
// A solid-shell is an eight-node hexahedron that carries bending correctly with
// **one element through the thickness**. It keeps three translational degrees of
// freedom per node -- no rotational DOF, no director, no drilling problem -- so it
// stacks against tetrahedra at an interface for free, which is what
// `docs/02-simulation.md` §3's "promotion from shell to tet as an element crumples"
// needs. What it costs is that the cures for locking have to be built in
// explicitly, because a plain trilinear hex at plate aspect ratios is *worse* than
// the tet it replaces.
//
// --- The three locking modes, and what is done about each ---------------------
//
// **Transverse shear locking.** Under cylindrical bending a trilinear hex carries
// u_x = A*xi*zeta, whose exact bending strain is right but which also produces a
// parasitic transverse shear gamma_xz proportional to xi. Full 2x2x2 integration
// sees that shear at every point and charges energy for it; the spurious energy
// scales as (in-plane size / thickness)^2, so a 10 mm plate with 50 mm elements is
// stiffened by a factor of ~25. Cure: **assumed natural strain** (Dvorkin-Bathe),
// sampling the covariant shear E_xi,zeta only at xi = 0 (mid-edge points), where
// the parasitic term is identically zero. Exact, parameter-free, no tuning.
//
// **Thickness (Poisson) locking.** With one element through the thickness the
// transverse normal strain E_zeta,zeta is constant in zeta. Bending needs it to be
// *linear* in zeta, because sigma_zz = 0 on both faces requires
// eps_zz = -nu/(1-nu) (eps_xx + eps_yy) and eps_xx is linear in zeta. Forced to a
// constant, the element behaves as plane strain through the thickness rather than
// plane stress. Measured on a strip in cylindrical bending: **22.5% too stiff**,
// which is exactly the ratio of the plane-stress modulus E/(1-nu^2) to the
// oedometer modulus E(1-nu)/((1+nu)(1-2nu)) -- a closed form, not a fudge, and
// tests/test_solid_shell.cpp asserts the two against each other. Cure: one
// **enhanced assumed strain** parameter
// carrying E_zeta,zeta = zeta * alpha, plus two more (xi*zeta, eta*zeta) so the
// curvature may vary across the element.
//
// **In-plane (membrane) locking of a distorted element.** A non-parallelogram hex
// under in-plane bending develops parasitic in-plane shear. Cure: the Simo-Rifai
// four-mode in-plane EAS set. Together with the three above that is seven enhanced
// parameters, condensed out at element level.
//
// A fourth, **curvature thickness locking**, appears once the mid-surface is
// curved -- which every shell plate on a ship is. E_zeta,zeta is biquadratic in
// (xi, eta) for a warped element; ANS sampling at the four in-plane corners and
// bilinear interpolation (Betsch-Stein) removes it. It costs nothing on a flat
// element, where the interpolation is already exact.
//
// **Reduced integration with hourglass control was rejected.** It is cheaper, and
// Flanagan-Belytschko hourglass vectors are orthogonal to linear fields so it
// would pass the patch test too. It is rejected because the hourglass stiffness is
// a tuned coefficient with no physical value: too small and the mesh develops
// zero-energy modes under the very loads this solver exists to compute; too large
// and it stiffens the physical bending mode it is supposed to leave alone. There
// is no measurement that sets it. ANS and EAS have no free parameters at all, and
// the price -- see the cost table in `07-fem-spike-findings.md` §6 -- is affordable.
//
// --- What this is in double, when fem.cpp is in float -------------------------
//
// `fem.hpp` uses flat float arrays because its CPU reference and its GLSL compute
// kernel have to operate on bit-identical std430 layouts. The things that
// establish an element is *correct* -- the patch test, rigid-body invariance, the
// rank of the element stiffness -- are exact identities, and in float their noise
// floor sits at 1e-6, the same order as the defects they exist to catch. So the
// element is formulated in double.
//
// **There is now a float GPU path -- `engine/gpu/zone_gpu.{hpp,cpp}` -- and the
// measurement it produced belongs here rather than only there.** Float is *not*
// sufficient for this element as formulated, and the reason is a property of the
// formulation and not of the GPU: `Kaa = int G^T C G dV`, the enhanced-strain
// block condensed out of every element, mixes modes that are not commensurate.
// The three thickness modes carry E_zeta,zeta so their columns of G are scaled by
// the Voigt transform's `1/t^2` -- 2.8e4 for 12 mm plating -- while the in-plane
// modes are scaled by `1/h^2`, of order ten. Kaa inherits the square of the ratio,
// so **kappa(Kaa) ~ (h/t)^4 ~ 1e7 on perfectly ordinary plating**, and float's
// seven digits leave none in alpha. Equilibrating Kaa before factoring recovers
// most of it and is what the shader does; the better fix is to normalise the
// enhanced modes here, which is a free basis change and would cost the double path
// nothing. `07-fem-spike-findings.md` §8 has the measurements.
//
// Body frame and SI units per CLAUDE.md.
#pragma once

#include "plasticity.hpp"
#include "scantlings.hpp"  // StructuralMaterial

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace sim::solidshell {

inline constexpr int kNodes = 8;   // nodes per element
inline constexpr int kDof = 24;    // 3 translations each
inline constexpr int kEas = 7;     // enhanced assumed strain parameters
inline constexpr int kGauss = 8;   // 2x2x2 Gauss points

// **Node ordering is load-bearing.** Nodes 0-3 are the zeta = -1 face and 4-7 the
// zeta = +1 face, each counter-clockwise seen from +zeta, with node 4 above node 0.
// zeta is therefore the **thickness direction**, and every assumed-strain
// modification in this file is written against that. An element handed its nodes
// in some other order is not wrong so much as no longer a shell: the ANS sampling
// would cure locking across the plate instead of through it. `smallestJacobian()`
// catches an inverted element; it cannot catch a merely rotated one, so the
// ordering is a contract.
//
// Natural coordinates of the nodes:
//   0(-1,-1,-1) 1(+1,-1,-1) 2(+1,+1,-1) 3(-1,+1,-1)
//   4(-1,-1,+1) 5(+1,-1,+1) 6(+1,+1,+1) 7(-1,+1,+1)

// Which cures are switched on. `Displacement` and `AssumedNaturalStrain` exist so
// the tests can attribute each share of the locking to the cure that removes it --
// a locking claim demonstrated on the same mesh, rather than asserted.
enum class Formulation {
    Displacement,          // plain trilinear hex, 2x2x2 Gauss. Locks. The control.
    AssumedNaturalStrain,  // + ANS transverse shear and thickness strain. No EAS.
    SolidShell,            // + 7-parameter EAS. This is the element.
};
const char* name(Formulation form);

// --- The step-invariant half of an element --------------------------------------
//
// **Everything an element derives from its rest configuration, computed once.**
//
// The strain-displacement matrices, the enhanced-strain interpolation, the Gauss
// weights and the rest Jacobian are functions of the *undeformed* geometry alone.
// An explicit solve never moves the rest configuration, so recomputing them every
// step is pure waste -- and it is not a small share of the bill. Measured by
// running the ferry's own side patch both ways, a punch into 192 elements for
// 6 608 steps: **2.01x end to end on one worker** (5.48 s against 2.73 s) and
// 1.64x on 23, for a **bit-identical** answer. Per element the elastoplastic
// update falls from 7.18 µs to 3.06 µs.
//
// It is worth being explicit about why this went unnoticed, because the shape of
// the mistake is reusable: the per-element cost was *measured* (7.3 µs, recorded
// in `02-simulation.md` §3) and the measurement was correct. What was never asked
// is which part of that 7.3 µs depended on the state being advanced. A cost model
// built from a correct total can still point at the wrong optimisation.
//
// **And the first two attempts to size it were both wrong in the same direction.**
// A standalone micro-benchmark and an in-situ `steady_clock` pair both reported
// 4.2 µs of a 4.3 µs update -- 97% -- and the in-situ pair contradicted itself
// arithmetically, its two sub-timers summing to more than the phase containing
// them. The A/B on the real run says 51%. `07-fem-spike-findings.md` §8 records it,
// because "time the real tick rather than extrapolating" arriving a second time is
// worth writing down.
//
// `fem.cpp`'s tetrahedron has had this all along: `TetMesh::restInverse` and
// `restVolume` are exactly this object for a linear tet, and `tet_forces.comp`
// reads them out of a buffer rather than rebuilding them. The solid-shell simply
// never grew the equivalent.
//
// **Size looked like the trade and measurement says it is not.** 12.0 kB per
// element in double, against the 4.6 kB per-element stiffness `02-simulation.md`
// §3 measures an L3 cliff for at ~6 500 elements -- so a crossover was expected at
// ~2 500. There is none. The cache still wins **1.62x at 17 800 elements and
// 214 MB**, four times past this machine's last-level cache, because the
// arithmetic it removes is 2.2 µs per element per step while streaming 12 kB at
// the saturated 29 GB/s that section measures is 0.41 µs. Five to one in the
// cache's favour even entirely from DRAM.
struct RestForms {
    double b[kGauss][6][kDof];    // strain-displacement, Cartesian, after ANS
    double g[kGauss][6][kEas];    // enhanced-strain interpolation, Cartesian
    double weight[kGauss];        // 2x2x2 Gauss weight times det J
    // (dX/dxi)^-1 at the element centre, row-major. The polar decomposition needs
    // it every step and it is the same matrix every step; this is the solid-shell's
    // `TetMesh::restInverse`.
    double restJacobianInverse[9];
    int  easCount = 0;
    bool ok = false;
};

// Fill `out` from an element's rest node positions. False -- and `out.ok` false --
// on an inverted or degenerate element, the same condition `elementStiffness`
// declines to work on.
bool computeRestForms(const double rest[kDof], Formulation form, RestForms& out);

// --- One element --------------------------------------------------------------
//
// `nodes` is always 8 nodes x 3 coordinates, node-major, in the ordering above.

// Element stiffness, 24x24 row-major, in the global frame. Symmetric to rounding;
// the enhanced parameters are condensed out here, so this is the whole element.
void elementStiffness(const double nodes[kDof], const StructuralMaterial& material,
                      Formulation form, double out[kDof * kDof]);

// Cartesian stress at the eight Gauss points, six Voigt components each in the
// order [xx, yy, zz, xy, yz, zx]. The enhanced parameters are recovered from the
// displacement first, so this is the stress the element actually carries -- which
// is what a patch test has to look at.
void elementStress(const double nodes[kDof], const double displacement[kDof],
                   const StructuralMaterial& material, Formulation form,
                   double out[kGauss * 6]);

// Row-sum lumped nodal mass: integral of N_a over the element, times density. Row
// sum rather than volume/8 because a distorted hex does not divide evenly and the
// row sum is the lumping that still moves rigid translation exactly.
void elementMass(const double nodes[kDof], double density, double out[kNodes]);

// The volume each of the eight Gauss points carries: the 2x2x2 weight times det J,
// summing to the element volume. Exposed because an energy balance over an
// elastoplastic element is a quadrature -- strain energy is `sum_gp w_gp * sigma :
// C^-1 : sigma` -- and volume/8 is exact only on a parallelepiped, so a solver that
// assumed it would report a distorted element's stored energy wrongly and in a way
// no deflection test can see.
void gaussVolumes(const double nodes[kDof], double out[kGauss]);

// Rotation carried by the element: the polar factor of the deformation gradient at
// the element centre. Column-major 3x3, matching fem.cpp's M3 layout.
void elementRotation(const double rest[kDof], const double current[kDof], double out[9]);
// The same rotation from a pre-computed `RestForms`, which spares the rest
// element's Jacobian and its inverse. **Bit-identical** to the call above -- the
// cached matrix is the one that call would have formed -- and `tests/test_solid_shell.cpp`
// asserts the identity rather than a tolerance, because anything looser would hide
// a stale or transposed cache.
void elementRotation(const RestForms& forms, const double current[kDof], double out[9]);

// Co-rotational internal force. `out` is the force the element applies **to its
// nodes**, the same sign convention as `fem::tetForces`, so a caller accumulates
// it directly. f = -R K (R^T x - X): exactly zero for any rigid body motion,
// including a finite rotation, which is the whole reason the co-rotational form is
// used with a linear material.
void internalForce(const double stiffness[kDof * kDof], const double rest[kDof],
                   const double current[kDof], double out[kDof]);
// The same, off a pre-computed `RestForms`. Bit-identical to the call above.
void internalForce(const RestForms& forms, const double stiffness[kDof * kDof],
                   const double rest[kDof], const double current[kDof], double out[kDof]);

// Largest stable explicit step for this element on its own: safety * 2/omega_max,
// with omega_max the largest eigenvalue of the lumped-mass generalised problem,
// by power iteration.
//
// This is computed rather than estimated from a length because the estimate is
// exactly what `07-fem-spike-findings.md` §4 got wrong. A solid-shell keeps its
// through-thickness stretch mode, so its highest frequency is thickness-governed
// however large the element is in plane -- see the measurement in §6.
double criticalTimestep(const double nodes[kDof], const StructuralMaterial& material,
                        Formulation form, double safety = 0.9);

// Smallest Jacobian determinant over the eight corners. Non-positive means the
// element is inverted or degenerate and nothing computed from it means anything.
double smallestJacobian(const double nodes[kDof]);

// --- Plasticity and tearing ---------------------------------------------------
//
// `plasticity.{hpp,cpp}` is the point law; this is the element that feeds it. The
// two things the element has to supply that a point law cannot are its **own
// size** -- because the failure strain is not a material constant, it depends on
// how much of a neck one element averages over -- and an **equilibrium condition
// for the enhanced strains**, because the seven EAS parameters are no longer given
// by a closed form once the material stops being linear.

// The element's characteristic in-plane length and its thickness, which
// `plasticity::regularisedFailureStrain` resolves the failure strain against.
// In-plane is sqrt(mid-surface area) and thickness is volume/area -- the
// perpendicular distance between the faces, so a *sheared* element reports the
// thickness it actually has rather than the length of its slanted edge.
void elementSize(const double nodes[kDof], double* inPlane, double* thickness);

// Everything one element remembers. Eight integration points, each with its own
// plastic history, plus the enhanced parameters kept as a warm start for the
// Newton below -- they are a function of the displacement and the history, not
// independent state, and starting from zero reaches the same answer at the cost of
// a few more iterations, which the tests check.
struct ElementPlasticState {
    plasticity::State point[kGauss];
    double enhanced[kEas] = {};
    double failureStrain = 0.0;  // resolved from this element's own geometry
    bool torn = false;           // every integration point has failed
};

// Sets `failureStrain` from the element's geometry and clears the history.
void initialisePlasticState(const double nodes[kDof], const plasticity::Material& material,
                            ElementPlasticState& state);

struct PlasticUpdate {
    bool converged = false;
    int iterations = 0;     // enhanced-parameter Newton iterations
    int yieldedPoints = 0;  // of kGauss
    int failedPoints = 0;
    double dissipation = 0.0;  // J over the element volume, this increment

    // ||int G^T sigma dV||. A diagnostic, **not** the convergence measure: the
    // enhanced thickness modes are scaled by 1/t^2, so for 20 mm plate a residual
    // of 1e-3 already means alpha is right to twelve digits. See the note in
    // `elementPlasticUpdate`.
    double enhancedResidual = 0.0;
    // |delta . r| of the last Newton correction, in joules -- the scale-free
    // measure that actually decides convergence, against the element's yield
    // energy sigma_y * V.
    double enhancedWork = 0.0;
};

// Elastoplastic internal force, and the history advanced to `current`.
//
// The enhanced strains are what makes this more than a loop over eight points. In
// the elastic element `alpha = -Kaa^-1 Kua^T u` closes in one line because
// everything is linear. Once the material yields, alpha is whatever satisfies
//
//     r(alpha) = integral G^T sigma(B u + G alpha) dV = 0
//
// -- the statement that the enhanced modes carry no stress -- and that is seven
// nonlinear equations solved here by Newton on the **algorithmic** tangent, which
// is why `plasticity::update` returns one. Skipping the solve and holding alpha at
// its elastic value would not be a small error: the enhanced thickness modes are
// the entire reason sigma_zz relaxes to zero through a bent plate, and a plate that
// cannot thin does not yield where a real one does.
//
// `force` is what the element applies **to its nodes**, the same sign convention
// as `internalForce` and `fem::tetForces`. With a material that never yields this
// returns exactly what `internalForce` returns for the condensed elastic
// stiffness, to rounding -- which is the tie between this path and the validated
// one.
//
// **An element that has lost an integration point drops its enhanced strains and
// finishes its life as the ANS hex**, and `state.enhanced` comes back zeroed to say
// so. `r(alpha) = 0` is a statement about a continuum, and an element with a dead
// point in it is not one: Kaa loses rank, and the measured consequence was not a
// wobble but an element that *stopped tearing* -- its surviving points driven below
// the damage cutoff, damage frozen at 0.78 while plastic strain ran on from 0.49 to
// 0.89. The step in which the point dies is re-run without the enhanced modes, so
// nothing computed on a rank-deficient system is ever committed.
//
// The history is committed even when the Newton did not converge, because a caller
// that has to keep stepping is better served by a slightly wrong state it is told
// about than by a state that silently stopped advancing. Check `converged`.
PlasticUpdate elementPlasticUpdate(const double rest[kDof], const double current[kDof],
                                  const plasticity::Material& material, Formulation form,
                                  ElementPlasticState& state, double force[kDof],
                                  double stress[kGauss * 6] = nullptr);

// The same update off a pre-computed `RestForms`, which is what an explicit solver
// wants: the forms are 97% of the cost above and none of it depends on `current`.
// **Bit-identical** to the call above -- same arithmetic, same order, on the same
// numbers -- which is the property that makes the cache safe to switch on and off,
// and which `tests/test_solid_shell.cpp` asserts on every output including the
// committed plastic history.
PlasticUpdate elementPlasticUpdate(const RestForms& forms, const double rest[kDof],
                                  const double current[kDof],
                                  const plasticity::Material& material,
                                  ElementPlasticState& state, double force[kDof],
                                  double stress[kGauss * 6] = nullptr);

// --- Meshes -------------------------------------------------------------------

struct HexMesh {
    std::vector<double>        position;    // 3 per node
    std::vector<std::uint32_t> index;       // 8 per element
    std::vector<std::uint8_t>  fixed;       // 3 per node; non-zero pins that DOF
    std::vector<double>        prescribed;  // 3 per node; the value a pinned DOF takes

    std::size_t nodeCount() const { return position.size() / 3; }
    std::size_t elementCount() const { return index.size() / 8; }

    void pin(std::size_t node, int axis, double value = 0.0);
    // Gather one element's 24 values out of a per-node array.
    void gather(std::size_t element, const std::vector<double>& nodal, double out[kDof]) const;
};

// A plate lx by ly by `thickness`, nx by ny elements in plane and nz through it.
// nz = 1 is the case the formulation exists for.
//
// Nodes are numbered with the thickness index fastest, then y, then x. That is not
// cosmetic: it is what keeps the assembled bandwidth at 3*(2*(ny+1)+2) instead of
// 3*((nx+1)*(ny+1)+...), which is the difference between a banded solve costing
// n*b^2 = 2e7 and 2e9 for the meshes the convergence study runs.
HexMesh makePlateMesh(double lx, double ly, double thickness, int nx, int ny, int nz = 1);

// Smallest stable explicit step over the whole mesh.
double criticalTimestep(const HexMesh& mesh, const StructuralMaterial& material,
                        Formulation form, double safety = 0.9);

// --- Static solution ----------------------------------------------------------
//
// A direct solver is here rather than in a test because validating an element
// against closed forms *is* a static problem: settling an explicit scheme to
// equilibrium instead would put the answer at the mercy of the damping and the
// step count, and a plate at a/t = 100 needs 10^5 steps to settle. It is also the
// only way to get a patch test to machine precision.

// Symmetric positive-definite banded direct solver, lower band stored.
//
// Exposed because the element's own validation has to build a *different* element
// on the *same* mesh and solve it the same way -- the linear tetrahedron the spike
// rejected, and the plain hex -- and a locking comparison down two different
// solver paths would prove nothing.
class BandedSpd {
public:
    BandedSpd(std::size_t dofCount, std::size_t halfBandwidth);

    // Accumulates into the lower triangle. Entries with row < column are ignored,
    // so a caller may loop over a full symmetric element matrix and each term
    // lands exactly once.
    void add(std::size_t row, std::size_t column, double value);

    // Cholesky in place. False means a non-positive pivot: the system is singular
    // or indefinite, which for a stiffness matrix means the boundary conditions
    // left a rigid body mode.
    bool factor();
    void solve(std::vector<double>& rightHandSide) const;  // in place, after factor()

    std::size_t dofCount() const { return n_; }
    std::size_t halfBandwidth() const { return b_; }

private:
    std::size_t n_, b_;
    std::vector<double> a_;  // a_[i * (b_ + 1) + (i - j)] is A(i, j), j <= i
};

// Stiffness the hex mesh does not carry, as a dense symmetric block over an
// arbitrary set of global degrees of freedom (3 * node + axis). It exists because
// `constraint.hpp` ties an eccentric stiffener to the plating and condenses it
// onto the shell's own DOF, and a second assembly path for that would be a second
// place the two could disagree.
//
// It used to say here that "an interface spring coupling two substructures is the
// same shape". It is, and **that is not what the interface coupling turned out to
// want**: `coupling.{hpp,cpp}` joins Tier 1 to Tier 2 with no extra stiffness at
// all, because a mesh split by element shares its interface nodes exactly and two
// components meeting there have the same displacement rather than a spring between
// them. A spring would need a stiffness, and a stiffness no measurement sets is
// precisely what this file rejects hourglass control for.
struct DofBlock {
    std::vector<std::uint32_t> dof;
    std::vector<double> stiffness;  // dof.size()^2, row-major, symmetric
};

// --- Multi-point constraints ----------------------------------------------------
//
// **A `DofBlock` can only add stiffness, and that is why a junction between two
// plates is not one.** `constraint.hpp` eliminates an eccentric stiffener's fibre
// with a `Tie` and hands the result over as a block, and that works because the
// fibre's endpoints are *not mesh nodes*: they have no rows of their own, so
// `T^T K T` over the masters is the whole of what the fibre contributes. A
// junction ties a node that **is** a mesh node, with elements of its own, and
// eliminating it means changing rows that already exist rather than adding new
// ones. No block can do that -- see `section.hpp` §5.
//
// So the constraint is carried as data and applied by the assembler. One slave
// degree of freedom, written as a weighted sum of masters:
//
//     u[slave] = sum_a weight[a] * u[master[a]]
//
// The slave keeps no unknown of its own; `solveStatic` and
// `reduction::Substructure` both scatter through the transformation, so the system
// they factor **is** `T^T K T` and nothing is enforced by a penalty. A weight set
// that sums to one reproduces every rigid body translation exactly, and one whose
// masters interpolate a point reproduces a rigid rotation exactly as well -- the
// same argument `constraint.hpp` §1 makes for its two-master case, which is this
// with `master.size() == 2`.
struct Mpc {
    std::uint32_t slave = 0;
    std::vector<std::uint32_t> master;
    std::vector<double> weight;  // same length as `master`
};

// What every degree of freedom stands for once a set of `Mpc`s is applied: itself,
// or the masters it was eliminated in favour of. Built once and shared by both
// assemblers, so there is one place that decides what a constrained scatter means.
//
// **Chained constraints are refused, not resolved.** A master that is itself a
// slave has a well-defined expansion and composing it silently is exactly the kind
// of thing that turns a mesher bug into a plausible answer -- so is a degree of
// freedom constrained twice. Both come back as `ok() == false` with a reason.
class DofExpansion {
public:
    struct Term {
        std::uint32_t dof;
        double weight;
    };

    DofExpansion() = default;
    DofExpansion(std::size_t dofCount, const std::vector<Mpc>& constrained);

    bool ok() const { return ok_; }
    const std::string& problem() const { return problem_; }
    // True when nothing is constrained, in which case every expansion is the
    // identity and both assemblers reduce to what they always did.
    bool empty() const { return eliminated_ == 0; }
    std::size_t eliminatedCount() const { return eliminated_; }
    std::size_t dofCount() const { return isSlave_.size(); }
    bool eliminated(std::uint32_t dof) const { return isSlave_[dof] != 0; }

    const Term* begin(std::uint32_t dof) const { return term_.data() + start_[dof]; }
    const Term* end(std::uint32_t dof) const { return term_.data() + start_[dof + 1]; }

    // Fill every eliminated degree of freedom of `values` from its masters. What
    // turns a solution over the free degrees of freedom back into one over the
    // mesh's own.
    void recover(std::vector<double>& values) const;

private:
    std::vector<std::size_t> start_;
    std::vector<Term> term_;
    std::vector<std::uint8_t> isSlave_;
    std::size_t eliminated_ = 0;
    bool ok_ = true;
    std::string problem_;
};

// Assemble and solve K u = f. Pinned DOF take their prescribed values and their
// contribution moves to the right-hand side, so a non-zero prescribed
// displacement -- which is what a patch test is -- is handled exactly rather than
// by a penalty. Returns false, with a reason in `problem`, on an inverted element
// or a singular system.
bool solveStatic(const HexMesh& mesh, const StructuralMaterial& material, Formulation form,
                 const std::vector<double>& load, std::vector<double>& displacement,
                 std::string* problem = nullptr);

// The same, with extra stiffness assembled alongside the elements. The blocks are
// counted in the bandwidth, so a block coupling distant DOF widens the band rather
// than being dropped on the floor.
bool solveStatic(const HexMesh& mesh, const StructuralMaterial& material, Formulation form,
                 const std::vector<DofBlock>& extra, const std::vector<double>& load,
                 std::vector<double>& displacement, std::string* problem = nullptr);

// The same, with multi-point constraints. Every scatter -- element, block, load and
// prescribed term alike -- goes through the expansion, so the system factored is
// `T^T K T` exactly and a slave degree of freedom is filled in afterwards from its
// masters rather than solved for.
//
// **A constrained degree of freedom that the mesh also pins is refused.** The two
// say different things about the same unknown and there is no reading of "both"
// that is not a guess; a caller who wants the boundary to win should not have
// constrained it.
bool solveStatic(const HexMesh& mesh, const StructuralMaterial& material, Formulation form,
                 const std::vector<DofBlock>& extra, const std::vector<Mpc>& constrained,
                 const std::vector<double>& load, std::vector<double>& displacement,
                 std::string* problem = nullptr);

// Consistent nodal loads for a uniform pressure acting on every +zeta face that
// lies on the boundary of the mesh -- the top surface of a plate. Positive
// pressure pushes *into* the element, so a plate loaded from above deflects in -z.
std::vector<double> uniformPressureLoad(const HexMesh& mesh, double pressure);

}  // namespace sim::solidshell
