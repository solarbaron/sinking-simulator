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
// kernel have to operate on bit-identical std430 layouts. There is no GPU kernel
// for this element yet, and the things that establish an element is *correct* --
// the patch test, rigid-body invariance, the rank of the element stiffness -- are
// exact identities. In float their noise floor sits at 1e-6, which is the same
// order as the defects they exist to catch. So the element is formulated in
// double; a float GPU path can be derived from it the way `fem_gpu.cpp` was
// derived from `fem.cpp`, and will have its own reproducibility bound (§2).
//
// Body frame and SI units per CLAUDE.md.
#pragma once

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

// Rotation carried by the element: the polar factor of the deformation gradient at
// the element centre. Column-major 3x3, matching fem.cpp's M3 layout.
void elementRotation(const double rest[kDof], const double current[kDof], double out[9]);

// Co-rotational internal force. `out` is the force the element applies **to its
// nodes**, the same sign convention as `fem::tetForces`, so a caller accumulates
// it directly. f = -R K (R^T x - X): exactly zero for any rigid body motion,
// including a finite rotation, which is the whole reason the co-rotational form is
// used with a linear material.
void internalForce(const double stiffness[kDof * kDof], const double rest[kDof],
                   const double current[kDof], double out[kDof]);

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

// Assemble and solve K u = f. Pinned DOF take their prescribed values and their
// contribution moves to the right-hand side, so a non-zero prescribed
// displacement -- which is what a patch test is -- is handled exactly rather than
// by a penalty. Returns false, with a reason in `problem`, on an inverted element
// or a singular system.
bool solveStatic(const HexMesh& mesh, const StructuralMaterial& material, Formulation form,
                 const std::vector<double>& load, std::vector<double>& displacement,
                 std::string* problem = nullptr);

// Consistent nodal loads for a uniform pressure acting on every +zeta face that
// lies on the boundary of the mesh -- the top surface of a plate. Positive
// pressure pushes *into* the element, so a plate loaded from above deflects in -z.
std::vector<double> uniformPressureLoad(const HexMesh& mesh, double pressure);

}  // namespace sim::solidshell
