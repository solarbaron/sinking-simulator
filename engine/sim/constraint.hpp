// SPDX-License-Identifier: MIT
//
// Multi-point constraints, and the eccentric stiffener they were built for.
//
// `zone.hpp` §3 records that a stiffener cannot be attached to solid-shell plating
// with the element inventory as it stood, and that what was needed is "a
// multi-point constraint tying an eccentric beam to the shell, which is the same
// machinery Tier-1/Tier-2 interface coupling needs and does not exist yet". This
// is that machinery, and §2 below is the stiffener built on it.
//
// **The half of that sentence about interface coupling turned out to be wrong, and
// it is worth correcting where it was written rather than only where it was
// disproved.** `coupling.{hpp,cpp}` couples Tier 1 to Tier 2 and uses none of this
// file. A structure split by *element* has interface nodes that are literally the
// same points, so the two components' shared DOF **are** one unknown and the
// coupling is scatter-add -- exact, and with nothing to interpolate. A tie is what
// a *non-matching* interface needs: two meshes at different resolutions either side
// of a cut, where a node on one side lands inside a face on the other. That case is
// real and will arrive with a whole-ship mesher; it had simply been assumed to be
// this case.
//
// --- 1. The tie, and why no rotational degree of freedom is needed -------------
//
// A solid-shell has `kNodes == 8` and `kDof == 24`: **three translations per node
// and no rotation anywhere**. It carries bending as *differential displacement
// between its two faces* -- that is the whole reason it stacks against a
// tetrahedron for free. So the rotation of a plate's cross-section is not a
// missing degree of freedom that has to be invented; it is already in the model,
// as the pair of nodes through the thickness.
//
// A point rigidly attached to the plating at a through-thickness offset `e` from
// the mid-surface is therefore an **exact linear function of that pair**:
//
//     u(e) = u_bottom + ((e + t/2) / t) * (u_top - u_bottom)
//
// -- one weight, two masters, no rotational degree of freedom anywhere, and no
// penalty stiffness or Lagrange multiplier. `e` is not confined to the plate: at
// `e = 0.1 m` on 12 mm plating the weight is 8.83 and the tie is an
// **extrapolation** far outside the element, which is exactly what an eccentric
// stiffener needs and is also where the surprises are (see the mass note below).
//
// **It is exact for a finite rigid rotation**, not just for small motions. A
// rotation R is linear, so `(1-w) R x_bottom + w R x_top = R ((1-w) x_bottom +
// w x_top) = R x(e)`, and a tied point therefore rides a rotating patch without
// straining. That is the same property the co-rotational element formulation
// exists to have, and `tests/test_constraint.cpp` asserts it as an identity.
//
// **The tie is condensed, not enforced.** The tied point has no degrees of freedom
// of its own at all: its force is distributed to the masters as `f_master +=
// w * f_tied` and its stiffness as `T^T K T`. The consequence is worth being
// explicit about, because it is the whole point of the exercise: **a formulation
// that eliminates the attached member's degrees of freedom cannot introduce a
// zero-energy mode.** The hinge `zone.hpp` §3 rejects has one because a web
// sharing one node row keeps an independent rotation about the seam; here there
// is no independent anything, and the null space of the stiffened system is
// exactly the null space of the plating it is tied to. Measured, not argued --
// see `tests/test_constraint.cpp`.
//
// **Where the tie costs something, and it is the mass.** The consistent condensed
// mass is `T^T M T`, whose row sums for one tied point are `(1-w) m` and `w m`.
// Those sum to `m` -- total mass is preserved for any weight -- but for `w = 8.83`
// they are `-7.83 m` and `+8.83 m`: **a negative nodal mass**, which an explicit
// scheme cannot integrate. There is no lumping of a positive mass over two nodes
// six millimetres apart that reproduces the first moment of a mass a hundred
// millimetres away; the first moment is what the extrapolation is for. So the
// stiffener's mass is split **equally** over the pair, which keeps the total
// exactly and gives up the stiffener's rotary inertia about the seam. It is
// recorded here rather than buried because the alternative -- the consistent
// lumping -- looks obviously right, and `tests/test_constraint.cpp` asserts that
// every nodal mass comes out positive so that nothing can quietly go back to it.
//
// **What the tie costs in arithmetic is nothing.** `fiberForces` is 20 ns per
// fibre per step, measured over 20 000 evaluations of a 240-fibre set; at the
// 0.42 fibres per element the reference ferry's resolution delivers that is 8 ns
// per element per step against a 3.1 µs elastoplastic element, and an end-to-end
// A/B cannot see it. The cost that is real is the **stable step** -- see
// `fiberFrequencySquared`.
//
// --- 2. The stiffener as axial fibres ------------------------------------------
//
// `scantlings.hpp` §1 keeps stiffeners as discrete line elements and says why:
// smearing loses the Steiner term `A d^2` about the combined neutral axis, which
// on a 200x10 flat bar on 12 mm plating at 700 mm spacing is a factor of **130**
// in the panel's second moment at identical area. So the eccentricity is the whole
// point, and a formulation that gets the offset wrong will look plausible and be
// badly wrong.
//
// The member is represented as a set of **axial fibres**: two-node bars running
// along the stiffener line, each at its own through-thickness offset, each tied to
// the plating by the constraint above. The "plane sections remain plane"
// hypothesis a beam element would assume is not assumed here -- it is *imposed by
// the tie*, which is the physical statement that the web is welded to the plate
// and does not shear off it.
//
// **The fibre stations are two-point Gauss through each rectangle of the profile,
// and that makes the section properties exact rather than approximate.** The
// energy density of an axial fibre is `E (eps0 + e kappa)^2 / 2`, quadratic in the
// offset, and two-point Gauss integrates a cubic exactly. So for a web of height
// `h` split at `h/2 +/- h/(2 sqrt 3)`,
//
//     sum A_i = A        sum A_i y_i = A h / 2      sum A_i y_i^2 = t_w h^3 / 3
//
// -- area, first moment and second moment all exact, which means the fibres carry
// the profile's **own** second moment as well as its Steiner term. Dropping to one
// fibre per rectangle would lose `I_own`, and for that flat bar `I_own` is 27% of
// the stiffened panel's total. `tests/test_constraint.cpp` checks the three
// moments against `scantlings::profileSection` rather than against this file's own
// arithmetic.
//
// **What this cannot do, stated rather than discovered.**
//
//   * **Tripping.** A bar has stiffness along its own axis and none across it, so
//     the fibres contribute *exactly zero* to a motion that tips the web about the
//     seam -- measured, and it is zero to rounding rather than small. What
//     restrains tripping here is the plating's own rotational stiffness, through
//     the tie, and the web is forced to follow the plate's cross-section exactly.
//     That is the opposite error from the hinge: the hinge leaves tripping free,
//     this over-restrains it. Lateral-torsional buckling of a stiffener remains
//     `buckling.hpp`'s question.
//   * **The weak-axis second moment.** Every fibre sits on the stiffener line, so
//     the profile's `secondMomentWeak` is not represented -- only its area is, in
//     bending of the panel about the plate normal. For a flat bar that is a 10 mm
//     web against a 200 mm one and the omission is small; for a tee with a wide
//     flange it is not.
//   * **The stiffener does not tear.** A fibre yields and hardens; it has no
//     damage variable and is never deleted. A zone whose plating has torn away
//     from under a stiffener still has the stiffener.
//
// SI units, body frame per CLAUDE.md.
#pragma once

#include "plasticity.hpp"
#include "scantlings.hpp"
#include "solid_shell.hpp"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace sim::constraint {

// --- 1. The tie ----------------------------------------------------------------

// One point whose translation is a fixed linear combination of a solid-shell's
// through-thickness node pair:
//
//     u = (1 - weight) * u[bottom] + weight * u[top]
//
// `bottom` is the node on the -zeta face and `top` the node on the +zeta face --
// in a `zone::Patch` that is the pair `(2m, 2m+1)`, with `Patch::outerFace`
// naming the second. A weight in [0, 1] interpolates inside the plate; outside it
// the tie extrapolates, which is what an eccentric member needs.
struct Tie {
    std::uint32_t bottom = 0;
    std::uint32_t top = 0;
    double weight = 0.5;
};

// The weight that puts a tied point `offset` metres from the mid-surface, measured
// along the direction the pair runs (bottom -> top), on plating `thickness` thick.
// Zero thickness has no pair direction and no answer; it returns 0.5, the
// mid-surface, which is the only weight that is right for *both* nodes when they
// coincide.
double tieWeight(double offset, double thickness);

// Where the tied point is, given a nodal coordinate or displacement array with
// three doubles per node.
Vec3 tiedPoint(const Tie& tie, const std::vector<double>& nodal);

// Distribute a force applied at the tied point onto the pair: the transpose of the
// interpolation, which is what makes the condensed system symmetric and the
// virtual work identical.
void applyTiedForce(const Tie& tie, const Vec3& force, std::vector<double>& nodal);

// --- 2. The eccentric stiffener -------------------------------------------------

// One axial fibre: a two-node bar whose ends are tied to the plating.
struct Fiber {
    Tie end[2];
    double area = 0;    // m^2 of the profile this fibre stands for
    double offset = 0;  // m from the plate mid-surface, along the pair direction
};

// What one fibre remembers. Uniaxial isotropic hardening; no kinematic hardening
// and no damage -- see the header note on what this cannot do.
struct FiberState {
    double plasticStrain = 0;            // signed, so a reversal unloads elastically
    double equivalentPlasticStrain = 0;  // monotone, drives the flow curve
};

// The fibre stations of one profile: offsets from the plate mid-surface and the
// area each stands for. Two-point Gauss per rectangle, so at most four for a tee.
//
// `sign` is `+1` when the web rises along the pair direction (bottom -> top) and
// `-1` when it rises the other way, which is the usual case for shell plating
// whose outward normal points away from the ship.
struct ProfileFibers {
    int count = 0;
    double offset[4] = {};  // m from the plate mid-surface, signed
    double area[4] = {};    // m^2
};
ProfileFibers profileFibers(const StiffenerProfile& profile, double plateThickness, double sign);

// Every fibre attached to one mesh, plus what they are made of.
struct Stiffening {
    std::vector<Fiber> fiber;
    StructuralMaterial material;

    double mass = 0;      // kg of stiffener steel
    double length = 0;    // m of stiffener line covered
    // How many `StructuralMember`s contributed a fibre. One per member, not one
    // per run: a member the patch carries in two disconnected pieces is still one
    // longitudinal, and a count that said two would be a property of the zone's
    // shape rather than of the ship's.
    int members = 0;
    std::vector<std::string> problems;

    std::size_t fiberCount() const { return fiber.size(); }
    bool empty() const { return fiber.empty(); }
};

// The step-invariant half of a fibre set: the rest length of every fibre, formed
// once from the rest configuration. The same object, for the same reason, as
// `solidshell::RestForms` -- an explicit solve never moves the rest configuration,
// and a fibre's rest length is two tie evaluations and a square root.
//
// **It is taken from the rest array and not from the mesh**, because a
// `zone::Preload` moves the rest configuration out from under the mesh, and a
// fibre whose rest length came from the meshed geometry would be handed the
// pre-strain as a free lunch.
struct RestFibers {
    std::vector<double> length;  // m, per fibre
    bool ok = false;             // false if any fibre came out with no length
};
RestFibers restFibers(const Stiffening& stiffening, const std::vector<double>& rest);

// What one force evaluation produced.
struct FiberForces {
    double strainEnergy = 0;  // J, recoverable
    double dissipation = 0;   // J, plastic, this increment
    int yielded = 0;          // fibres that yielded this increment
};

// Co-rotational axial force of every fibre, **added into** `force` through the
// ties. `f = N n` on the first end and `-N n` on the second, with `n` the fibre's
// *current* direction and `N = sigma A`; the strain is `(l - L) / L`, which is
// exactly zero for any rigid body motion including a finite rotation, so a fibre
// on a patch that has rotated carries no force. That is the same statement
// `solidshell::internalForce` makes, and the fibres have to make it too or a
// rotating zone would grow force out of nothing.
//
// `state` may be null, which solves the elastic bar and is the path the geometric
// tests use.
FiberForces fiberForces(const Stiffening& stiffening, const RestFibers& forms,
                        const std::vector<double>& current, const plasticity::Material& material,
                        std::vector<FiberState>* state, std::vector<double>& force);

// Add the fibres' mass to a per-node lumped mass array, **split equally over each
// pair**. Not `T^T M T`: see the header. The total is exact; the first moment
// about the mid-surface is not represented, so the stiffener's rotary inertia
// about the seam is given up.
void lumpFiberMass(const Stiffening& stiffening, const RestFibers& forms, double density,
                   std::vector<double>& nodalMass);

// --- 3. Stiffness, for the static path and for the step -------------------------

// The condensed stiffness of one fibre. A bar is rank one -- `K = (EA/L) d d^T` --
// and the tie is linear, so the condensed 12 x 12 block is rank one too:
//
//     K_condensed = scale * v v^T,   scale = EA/L
//
// over the four master nodes (two ties, two nodes each). Stored as the vector
// rather than the block because that is both smaller and what makes the
// eigenvalue below exact.
struct FiberStiffness {
    std::uint32_t dof[12] = {};  // global DOF, 3 * node + axis
    double vector[12] = {};
    double scale = 0;  // EA / L
};
FiberStiffness fiberStiffness(const Fiber& fiber, const std::vector<double>& rest,
                              double restLength, double youngsModulus);

// The largest eigenvalue of `M^-1 K` for one fibre, in rad^2/s^2. **Exact, not a
// bound**: a rank-one stiffness `k v v^T` against a diagonal mass has the single
// non-zero eigenvalue `k * sum_i v_i^2 / m_i`, so no power iteration is needed and
// none of the usual Gershgorin slack is paid.
//
// This is why the stiffener's effect on the stable step is knowable before the
// run: the tie amplifies a fibre's stiffness by the square of its weight -- 373x
// for the outer fibre of a 200 mm bar on 12 mm plating -- while its mass arrives
// unamplified, and whether that beats the plate's own thickness-governed frequency
// is a measurement rather than a guess.
double fiberFrequencySquared(const FiberStiffness& stiffness,
                             const std::vector<double>& nodalMass);

// Every fibre's condensed stiffness as a block `solidshell::solveStatic` takes, so
// the static path and the explicit one are the same physics assembled twice rather
// than two formulations. One block per fibre: 12 DOF, rank one.
std::vector<solidshell::DofBlock> stiffnessBlocks(const Stiffening& stiffening,
                                                  const std::vector<double>& rest,
                                                  const RestFibers& forms, double youngsModulus);

// Smallest stable explicit step the fibres allow, given the assembled nodal mass:
// `safety * 2 / omega_max` over every fibre. Infinite -- returned as zero -- when
// there are no fibres, so a caller takes the plating's own step.
double criticalTimestep(const Stiffening& stiffening, const RestFibers& forms,
                        const std::vector<double>& rest, const std::vector<double>& nodalMass,
                        double youngsModulus, double safety = 0.9);

// --- 4. Building the fibres over a plating mesh ---------------------------------

// A run of shell node pairs along one stiffener line, in order, together with the
// direction the web rises. `pair[i]` is the *bottom* node index; the top node is
// `pair[i] + 1` for a `zone::Patch`, so the caller supplies both explicitly rather
// than relying on that.
struct SeamRun {
    std::vector<std::uint32_t> bottom, top;
    double sign = -1.0;  // +1 if the web rises bottom -> top, -1 the other way
};

// Fibres for one profile along one run of pairs. Appends to `out`; returns how
// many fibres were added. A run shorter than two stations adds nothing.
std::size_t addStiffener(const SeamRun& run, const StiffenerProfile& profile,
                         double plateThickness, const std::vector<double>& rest,
                         Stiffening& out);

}  // namespace sim::constraint
