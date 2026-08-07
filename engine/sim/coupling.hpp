// SPDX-License-Identifier: MIT
//
// Tier-1 to Tier-2 coupling: a zone's edge driven by the reduced model around it,
// and what a torn zone hands back.
//
// `zone.hpp` §5 item 1 records the hole this closes, and `promotion.hpp`'s opening
// note records the same thing from the other side: a zone is coupled to
// **Tier 0** through a section, where the three-tier plan wants it coupled to
// Tier 1 through retained interface DOF. The reduction that makes that possible
// has existed since `reduction.{hpp,cpp}` -- boundary DOF kept exactly, in global
// numbering -- and so has the synthesis that joins two of them. What did not exist
// is the piece in this file: something that takes the interface displacement out of
// an assembled reduced model and *puts it on the zone's perimeter*, and something
// that puts the zone's own stiffness back into the model so the surroundings feel
// it.
//
// Both directions are here, and the second is not an extra. A zone whose edge is
// driven but whose softness nothing outside can see is still one-way coupling
// wearing a better boundary condition.
//
// --- 1. What "coupled" is, exactly ---------------------------------------------
//
// A zone and its surroundings are two components sharing an interface, which is
// precisely the object `reduction::assemble` already builds. So the coupling is:
//
//   1. Cut the structure by **element**, into a surrounding region and the zone.
//      By element, never by node: `reduction.hpp` §synthesis says why -- a node
//      split would double the interface mass.
//   2. Reduce each with `craigBampton`. The surroundings are Tier 1 proper; the
//      zone is reduced *only to be assembled*, and its own answer never comes from
//      the reduction.
//   3. `matchBoundaries` + `assemble` + a static solve, which gives the
//      displacement of every shared interface DOF.
//   4. Put that displacement on the zone's perimeter as a **prescribed** boundary
//      condition, and solve the zone -- statically here, or with the Tier-2
//      explicit solver, which now honours `HexMesh::prescribed` (`zone.hpp`
//      §4).
//
// **In the linear case this is not an approximation and the tests assert it as an
// identity.** `reduction.hpp` §1 property 1 is that Guyan static condensation is
// exact at the interface for *any* load at *any* mode count; assembly is
// scatter-add and approximates nothing; so the interface displacement the
// assembled model reports is the monolithic model's own, to the conditioning of
// two independent solves. Driving the zone's full model with it reproduces the
// monolithic field inside the zone, exactly.
//
// The consequence is worth stating plainly because it contradicts the natural
// expectation: **statically, mode count buys nothing here.** A coupled answer at
// zero modes is not an approximation that improves -- it starts exact and stays
// there, and `tests/test_coupling.cpp` asserts the agreement across mode counts
// rather than a convergence. What modes buy is dynamics and interior recovery
// under an interior load, which is what §5 of `reduction.hpp` measures.
//
// --- 2. Prescribing on a reduced model -----------------------------------------
//
// `reduction::assembledStaticSolve` holds DOF at **zero**, and a coupling needs
// them held at a value: a punch is a prescribed motion, and so is a hull girder
// pulling on a cut. `prescribedStaticSolve` does it by elimination rather than by
// penalty --
//
//     K_ff x_f = f_f - (K x_p)_f
//
// -- which is a right-hand side and nothing else, so it reuses the existing solve
// unchanged instead of introducing a second factorisation path or a penalty
// stiffness that no measurement sets. That is the same argument `zone.hpp` §4
// makes against a tuned elastic edge and `solid_shell.hpp` makes against hourglass
// control.
//
// The DOF a caller wants to prescribe has to *be* in the reduced model, which
// means it has to be an interface DOF of its component. That is not a workaround:
// `reduction.hpp` §3 defines the interface as every node the substructure shares
// with anything outside it, and a punch is outside.
//
// --- 3. What the zone hands back, and why it is not a thinner ship -------------
//
// `promotion.hpp` §5 folds Tier-2 damage back into Tier 0 as an **effective
// thickness**, and for Tier 0 that is exactly right: `hullGirderSection`,
// `girderBuckling` and `collapseElementsAt` all read a thickness and nothing else,
// so a thinner ship needs no new model anywhere.
//
// **It is the wrong shape for Tier 1, and the reason is geometric rather than a
// matter of taste.** A Tier-1 component is a `solidshell::HexMesh`, and a
// solid-shell carries its plate thickness in the *positions of its nodes* -- the
// mid-surface at zero and the faces at +/- t/2. Thinning a zone therefore moves
// every one of its nodes, including the interface nodes it shares with the
// surroundings, so the two components stop being coincident and
// `matchBoundaries` matches nothing at all. A thickness knockdown is not a
// modification of a component; it is a different component that no longer fits.
// `tests/test_coupling.cpp` measures that rather than asserting it, because it is
// the kind of claim that sounds like an excuse.
//
// What Tier 1 can take, exactly and with no remeshing, is **element deletion**.
// A torn element carries nothing, which is a statement about the stiffness
// assembly and not about geometry, so `withoutTornElements` removes it and hands
// back a mesh whose surviving interface nodes are where they always were. The
// surrounding reduced model then feels a genuinely softer zone, because it is
// assembled against a genuinely softer component.
//
// **A zone that has *yielded* without tearing is softer than its elastic self, and
// that used to be missing too.** It is §5 now. The note that stood here said
// closing it needed "a reduction built from a tangent operator the Tier-2 explicit
// solver does not form", and both halves of that were wrong: a static coupling
// wants a **secant**, not a tangent, and the secant of a J2 point is a closed form
// in the equivalent plastic strain the solver has carried all along. The
// measurement is in §5 and in `tests/test_coupling.cpp`; what the note got right is
// the direction, which is measured at **1.95x** the true perimeter travel a fifth
// of the way to a tear and 2.27x at a half.
//
// --- 4. The stiffeners, and what a coupling has to hand over --------------------
//
// A zone under `zone::Stiffeners::Modelled` carries its members as
// `constraint::Stiffening` -- fibres condensed onto the plating's own DOF, with no
// nodes and no elements of their own. A `reduction::Substructure` built from the
// mesh and a material therefore reduced a stiffened patch as **the bare plating**,
// and one 200x10 flat bar across a 0.6 m patch is worth 7.9% of its displacement
// field, so that was not a rounding.
//
// **It is closed, and the shape of the fix is worth knowing because the coupling
// has to do something about it.** `reduction::Attachment` (`reduction.hpp` §8)
// takes exactly the `solidshell::DofBlock` list `constraint::stiffnessBlocks`
// already produces and exactly the nodal mass `constraint::lumpFiberMass` already
// lumps. The substructure folds both into the one CSR and the one lumped diagonal
// it assembles the elements into, so the bandwidth, the interior renumbering and
// the reduction all carry the stiffener with nothing added downstream.
//
// **The caller has to hand it over.** A substructure built from a stiffened patch
// with no `Attachment` is still the bare plating, to the last bit -- that is
// asserted in `tests/test_coupling.cpp` as the negative control, next to the
// measurement that the same patch *with* one reproduces `solveStatic`'s own
// stiffened field to 4e-14 of it. Nothing in `couple()` reaches into a
// `zone::Patch` to find its `stiffening`, because a `Substructure` is built by the
// caller and a coupling is handed two that already exist.
//
// Separately, and unchanged: `constraint.hpp` gives a fibre no damage variable and
// never deletes it, so "this longitudinal is gone" does not exist to be read --
// `promotion.hpp` §5 records the same gap for Tier 0. A zone whose plating has
// torn out from under a longitudinal reports the plating and keeps the member.
// That one is a failure criterion for the fibres, not coupling machinery. It has
// grown a second edge now that Tier 1 can see the member: `withoutTornElements`
// hands back a mesh with elements removed and renumbered, and an `Attachment` built
// against the *original* numbering does not survive that -- it has to be rebuilt
// from the fibres against the damaged mesh, which nothing here does yet.
//
// --- 5. Where the exactness stops ------------------------------------------------
//
// The interface displacement is exact for **any** load, including one applied
// inside a component -- that is `reduction.hpp` §1 property 1 in full, and it is
// what the zone's own answer rests on, because the zone is re-solved in its full
// model rather than recovered from its reduction.
//
// What is *not* exact at zero modes is a component's recovered **interior** field
// under an interior load: `recover` gives `Psi u_b + Phi q`, and with no modes the
// particular solution is missing. The surroundings are recovered that way, so a
// pressure on the surrounding plating is the case where mode count starts to
// matter even statically. Every measurement in `tests/test_coupling.cpp` loads
// through prescribed interface DOF alone, which is why the recovered surrounding
// field comes back at 7e-15 m rather than at a truncation error.
//
// SI units, body frame per CLAUDE.md.
#pragma once

#include "reduction.hpp"
#include "solid_shell.hpp"
#include "zone.hpp"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace sim::coupling {

// --- 1. A boundary condition on a reduced model ---------------------------------

// One assembled DOF held at a value rather than at zero.
//
// **This moved.** It only ever touched an `Assembly`, and a chain of Tier-1 sections
// needs the same thing from the other side of this file -- a ship is loaded by
// prescribing a plane-sections field on its two end cuts -- so it lives in
// `reduction.hpp` next to the solve it eliminates into. The names here are the same
// names and mean the same thing; §2 below is still where the reasoning is written
// down.
using Prescribed = reduction::Prescribed;

// K x = f with `held` at zero and `prescribed` at their values. See §2: the
// prescribed values are eliminated into the right-hand side, so this is
// `reduction::assembledStaticSolve` with a different load and a longer held list,
// and there is exactly one dense factorisation in the engine either way.
//
// A DOF listed in both is prescribed: `held` is the special case value zero, and
// resolving the clash the other way would silently discard the value.
bool prescribedStaticSolve(const reduction::Assembly& assembly,
                           const std::vector<double>& load,
                           const std::vector<std::uint32_t>& held,
                           const std::vector<Prescribed>& prescribed,
                           std::vector<double>& state, std::string* problem = nullptr);

// --- 2. Matching a zone to the structure round it --------------------------------

// A surrounding region and a zone, matched at the interface and assembled.
//
// Thin over `reduction::matchBoundaries` and `reduction::assemble`, which do the
// work. What it adds is the bookkeeping a caller would otherwise repeat -- which
// assembled DOF each component's boundary landed on, and which of the zone's
// boundary DOF the surroundings actually touch.
struct Coupling {
    reduction::InterfaceMap map;  // surroundings -> zone, in boundary DOF indices
    reduction::Assembly assembly;

    // Assembled index of each component's boundary DOF: `assembly.fromA` and
    // `fromB` restricted to the boundary, named for what they are so a caller does
    // not have to remember which component went in first.
    std::vector<int> surroundDof, zoneDof;
    // Per zone boundary DOF: non-zero where the surroundings share it. **Not every
    // zone boundary DOF is shared, and that is not an error** -- a zone whose
    // interface deliberately carries a punch, or any other DOF something outside
    // the reduced model drives, has boundary DOF nothing in the surroundings
    // touches. What *is* an error is a perimeter DOF that goes unmatched, and only
    // the caller knows which of its boundary DOF are perimeter.
    std::vector<std::uint8_t> zoneShared;

    std::size_t sharedDof = 0;     // how many DOF the two components hold in common
    std::size_t zoneUnshared = 0;  // zone boundary DOF the surroundings do not touch
    double worstGap = 0;           // m, the furthest apart a matched pair actually was
    std::vector<std::string> problems;

    bool ready() const { return !assembly.empty() && sharedDof > 0; }
};

Coupling couple(const reduction::Substructure& surroundings,
                const reduction::Reduction& reducedSurroundings,
                const reduction::Substructure& zone, const reduction::Reduction& reducedZone,
                double tolerance = 1e-9);

// --- 3. Driving the zone's edge ---------------------------------------------------

// The displacement the surrounding model puts on the zone's own nodes, in the
// zone mesh's global DOF numbering (3 * node + axis).
//
// **Shared DOF only.** A zone boundary DOF the surroundings do not touch is left
// undriven, because nothing outside the zone said anything about it -- driving it
// at the zero it happens to hold would be inventing a clamp.
struct EdgeDrive {
    std::vector<double> displacement;   // 3 per zone node; zero away from the interface
    std::vector<std::uint8_t> driven;   // 3 per zone node; non-zero where imposed
    std::size_t count = 0;              // DOF actually driven
    double largest = 0;                 // m, the biggest of them
    std::vector<std::string> problems;

    bool empty() const { return count == 0; }
};

EdgeDrive edgeDrive(const Coupling& coupling, const reduction::Substructure& zone,
                    const std::vector<double>& assembledState);

// Pin a mesh's driven DOF to the drive: `fixed` set and `prescribed` written.
//
// It **pins**, rather than assuming the caller already did. A `zone::buildPatch`
// patch arrives with its perimeter pinned at zero and this replaces the zero; a
// bare `makePlateMesh` arrives with nothing pinned and this is what holds it. What
// it never does is *un*-pin: a DOF the mesh pins and the drive does not name keeps
// the boundary condition it had, so applying a drive to a partly matched interface
// leaves the unmatched part clamped rather than quietly setting it free.
//
// Returns how many DOF it wrote, which a caller should compare against
// `EdgeDrive::count` -- they differ only if the drive names a node the mesh does
// not have.
std::size_t applyEdgeDrive(const EdgeDrive& drive, solidshell::HexMesh& mesh);

// --- 4. What a torn zone hands back ------------------------------------------------

// A mesh with elements removed, renumbered so that nothing refers to what is gone.
//
// Nodes left with no element at all are **dropped**, and that is not tidying: an
// orphan node's stiffness row is identically zero, so a substructure containing
// one has a singular interior and `reduction.hpp` §3 records that the banded
// factorisation does not reliably catch it -- an exactly zero pivot comes back as
// a tiny positive one and the solve returns nonsense.
struct DamagedMesh {
    solidshell::HexMesh mesh;
    std::vector<std::uint32_t> element;  // per kept element, its index in the original
    std::vector<std::uint32_t> node;     // per kept node, its index in the original

    std::size_t removedElements = 0;
    std::size_t orphanedNodes = 0;  // dropped because every element on them went
    std::vector<std::string> problems;
};

// `removed` is one byte per element, non-zero to delete.
DamagedMesh withoutElements(const solidshell::HexMesh& mesh,
                            const std::vector<std::uint8_t>& removed);

// The same, reading the tear out of a solved zone. An element `zone::Solver`
// reports torn carries no stress and no stiffness, which is exactly the element
// this leaves out.
DamagedMesh withoutTornElements(const zone::Patch& patch, const zone::Solver& solver);

// The nodes of a damaged mesh that were interface nodes of the original, in the
// damaged mesh's own numbering. A tear that reaches the perimeter removes interface
// nodes, and the surrounding model then has boundary DOF the zone no longer shares
// -- which is the physically right answer and is why `couple` counts them rather
// than refusing.
std::vector<std::uint32_t> carriedInterface(const DamagedMesh& damaged,
                                            const std::vector<std::uint32_t>& interfaceNodes);

// --- 5. What a *yielded* zone hands back --------------------------------------------
//
// §3 said this could not be done without "a reduction built from a tangent operator
// the Tier-2 explicit solver does not form". **That framing was wrong, and the
// measurement that says so is in `tests/test_coupling.cpp`.** Two things were wrong
// with it and they point in opposite directions.
//
// **A tangent is not what a static coupling wants.** The assembled solve is
// `K x = f` on a *total* load, so the operator that makes it right is the one for
// which `K u = f_int(u)` at the state the zone is actually in -- a **secant**. A
// tangent answers `K du = df`, which is a different question; used for a total
// solve it over-softens, and it over-softens worst exactly where the softening is
// smallest, because `dsigma_y/deps_p` drops by a finite step at the first
// increment of flow while the secant leaves the elastic modulus continuously
// (`plasticity.hpp`, "what a yielded point has left"). Both are built here, the
// tangent as the control. Measured: just past first yield the secant leaves 0.37%
// of the surroundings' own displacement, the tangent 31.3%, and reporting nothing
// at all 1.14% -- so the tangent is 85x the secant there and **27x worse than
// changing nothing**. Deep in the band it is 3.4x the secant. It over-softens
// throughout, where the model it was meant to replace over-stiffens.
//
// **And nothing had to be *formed* that the solver does not already have.** The
// secant modulus of a J2 point is a closed form in its equivalent plastic strain
// alone, which `zone::Solver` has carried per integration point since it existed.
// What was missing was not an operator but the arithmetic that turns that number
// into element stiffness.
//
// --- How it reaches the reduced model, with nothing added to `reduction` ---------
//
// A `reduction::Attachment` is scatter-added into the same CSR the elements go
// into, and nothing anywhere requires a block to be positive. So a **negative**
// block is a stiffness *correction*: for each softened element,
//
//     dK_e = K_e(E_s, nu_s) - K_e(E, nu)
//
// on the element's own 24 DOF, so that the assembly the substructure forms is
// exactly `K_e(E_s, nu_s)`. It names the degrees of freedom the element already
// couples, so the sparsity, the bandwidth and the interior renumbering are
// untouched, and `reduction.{hpp,cpp}` needed no change at all. `solid_shell.hpp`
// says a `DofBlock` "can only add stiffness"; that is a statement about what a
// *constraint* can be expressed as, not a sign convention, and §5 of this file is
// the case that reads it the other way.
//
// The plate's bulk modulus is kept and the whole knockdown is in the shear
// modulus, because J2 flow is deviatoric and never touches the trace.
//
// --- Where it stops, and how it composes with a tear ---------------------------
//
// **Torn elements are skipped, and that is not an omission.** A torn element's
// secant stiffness is not small, it is *absent*, and handing back a block equal to
// `-K_e` leaves the assembly with rows no element supports -- which is the
// singular interior §4 removes an orphan node to avoid, and which
// `reduction.hpp` §3 records the banded factorisation does not reliably catch.
// Tearing goes back through `withoutTornElements` and softening goes back through
// this; they are different mechanisms and they compose in that order. `torn` counts
// what was skipped so a caller that has not deleted them can see it.
//
// **A secant is a statement about the state it was measured at.** It is exact where
// the point is at yield, which is every point of a monotonically loaded zone. A
// point that has yielded and then unloaded is reported with the stiffness it had at
// its peak rather than the elastic stiffness it has now, because the secant is
// built from `sigma_y(eps_p)` and not from the stress the point currently carries
// -- which is the bounded, monotone choice, where the ratio of the two would run to
// zero as an unloaded point's stress does.
//
// **And it is a linearisation, so it does not iterate itself.** The zone's state
// sets the softening, the softening sets the interface displacement, and the
// interface displacement sets the zone's state. `tests/test_coupling.cpp` measures
// that loop from a single run, by keeping what each pass recovered: on the plate
// there, three passes leave 9.22% and six leave 8.55%, and the sixth pass moves the
// field by 0.018% of the peak against the third-to-sixth 1.27% -- a factor of 72, so
// the loop has closed and the 8.55% is the isotropic knockdown rather than the
// iteration. Nothing here runs the loop, because how many passes a caller can afford
// is a caller's decision --
// and a `Softening` is cheap beside it, one element stiffness pair per yielded
// element against a whole zone re-solve.

// Which modulus a softened element is given. `Secant` is the answer; `Tangent` is
// the control the docs used to prescribe, kept because a claim that one of them is
// wrong is worth being able to re-measure rather than re-read.
enum class Modulus { Secant, Tangent };

struct Softening {
    // Ready to hand to `reduction::Substructure`. The mass is a zero per node --
    // stated rather than left empty, because plastic flow moves no steel and the
    // substructure would otherwise report an attachment that is "stiffer but no
    // heavier", which is the opposite of what this is. The stress forms are zeros
    // parallel to the blocks: a softening correction is not a member and carries no
    // member stress, and `reduction::checkValidity` takes a maximum over them, so a
    // zero is inert there. It does count in `Substructure::attachedMembers`.
    reduction::Attachment attachment;

    std::vector<std::uint32_t> element;  // per block, which element it corrects
    std::size_t softened = 0;            // elements that had flowed
    std::size_t torn = 0;                // elements skipped because they are gone

    // The knockdown, as G_s/G. `worstRatio` is the softest single element and
    // `meanRatio` is the volume-weighted average over the whole patch, torn
    // elements counted at zero -- so an intact elastic zone reports 1.0 for both and
    // a caller can read "how stale is my reduced model" off `meanRatio` without
    // reducing or assembling anything. That is the cheap honest answer for a caller
    // that would rather be told its model is stale than have it silently corrected:
    // the ratios come out of the plastic state and the Gauss volumes alone.
    double worstRatio = 1.0;
    double meanRatio = 1.0;
    double peakPlasticStrain = 0.0;

    std::vector<std::string> problems;

    bool empty() const { return softened == 0; }
};

// The correction for a mesh whose elements carry `state`. `elastic` is the material
// the substructure is being built with -- the correction is referenced to it, so
// the two must be the same material or the block corrects a stiffness nobody
// assembled. `material` is where the flow curve comes from.
//
// **An element with no plastic strain gets no block at all**, rather than a block
// of zeros, so a zone that has not yielded builds *bit for bit* the substructure it
// built before this existed. That is the negative control the whole of §5 rests on
// and it is a property of the code rather than of a tolerance.
Softening softening(const solidshell::HexMesh& mesh, const StructuralMaterial& elastic,
                    const plasticity::Material& material,
                    const std::vector<solidshell::ElementPlasticState>& state,
                    Modulus modulus = Modulus::Secant,
                    solidshell::Formulation form = solidshell::Formulation::SolidShell);

// The same, reading the state out of a solved zone. The mesh is the patch's, so a
// caller that has already deleted torn elements should use the call above against
// the damaged mesh and the states it kept -- this one is for a zone that has
// yielded and not torn, which is the case §5 exists for.
Softening softening(const zone::Patch& patch, const zone::Solver& solver,
                    const plasticity::Material& material, Modulus modulus = Modulus::Secant);

}  // namespace sim::coupling
