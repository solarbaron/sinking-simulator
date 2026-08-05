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
// **What it cannot carry, and this is a real limit rather than a rounding.** A
// zone that has *yielded* without tearing is softer than its elastic self and
// nothing here says so: the reduction is linear by construction
// (`reduction.hpp` §6) and there is no tangent stiffness to hand it. Between
// first yield and first tear, direction 2 under-reports the softening -- always
// in the stiff direction, which is the unsafe one. Closing that means a reduction
// built from a tangent operator the Tier-2 explicit solver does not form.
//
// --- 4. What Tier 1 cannot see at all: the stiffeners ---------------------------
//
// This is larger than it looks and it is not about damage. A
// `reduction::Substructure` is built from a `solidshell::HexMesh` and a material.
// A zone under `zone::Stiffeners::Modelled` carries its members as
// `constraint::Stiffening` -- fibres condensed onto the plating's own DOF, with no
// nodes and no elements of their own. **So the substructure a coupling builds from
// a stiffened patch is the bare plating**, and `tests/test_coupling.cpp` measures
// exactly that: the substructure's operator matches an element-only assembly to
// 3e-16 relative, while one 200x10 flat bar across a 0.6 m patch is worth 7.9% of
// its displacement field.
//
// The fix is named rather than done, because it belongs in `reduction.cpp` and
// wants its own validation against a monolithic stiffened plate:
// `constraint::stiffnessBlocks` already produces precisely the
// `solidshell::DofBlock` list `solidshell::solveStatic` takes, and what is missing
// is a `Substructure` that accepts extra blocks and scatters them into the same
// CSR it assembles the elements into.
//
// Separately, and unchanged: `constraint.hpp` gives a fibre no damage variable and
// never deletes it, so "this longitudinal is gone" does not exist to be read --
// `promotion.hpp` §5 records the same gap for Tier 0. A zone whose plating has
// torn out from under a longitudinal reports the plating and keeps the member.
// That one is a failure criterion for the fibres, not coupling machinery.
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
struct Prescribed {
    std::uint32_t dof = 0;  // index into the assembled model
    double value = 0;       // m
};

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

}  // namespace sim::coupling
