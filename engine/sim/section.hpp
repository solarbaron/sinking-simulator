// SPDX-License-Identifier: MIT
//
// Cut a ship into a piece worth reducing: the plating between two transverse
// planes, meshed as solid-shell elements, with its longitudinals attached, and the
// two cut sections as the interface a `reduction::Substructure` is built on.
//
// `zone::buildPatch` is the only other thing that turns ship structure into a
// `solidshell::HexMesh`, and it is centred on an impact point, holds one panel
// role, one thickness and one material, and stops at every seam that is not all
// three. That is right for a Tier-2 patch and it is exactly wrong for a Tier-1
// component, which is a *region* of the ship rather than a neighbourhood of an
// event. `reduction.hpp` §3 names the region it wants -- "between two bulkheads,
// cut at the two transverse sections it is attached through" -- and
// `nodesNearPlanes` has been sitting there waiting for a mesher to hand it one.
//
// --- 1. What this found out first, because it changes what the file can claim ---
//
// **`makeStructuralMesh` produces three topologically disjoint panel sets, and no
// mesher can weld them.** Measured on the reference ferry: of 9 390 distinct panel
// corners, the number shared between a `Shell` panel and a `Deck` panel is
// **zero**; shell-to-bulkhead is zero; deck-to-bulkhead is zero; and even two
// *different* bulkheads share none. The three roles are laid out on three
// independent grids -- the shell on girth fractions of each station, a deck on
// fixed |y| lines clipped to the hull, a bulkhead on its own -- so a deck edge
// lands in the *middle* of a shell panel, missing the nearest shell corner by up
// to 0.31 m on this ship.
//
// So a section of a real ship comes out as several disconnected surfaces, and that
// is a property of the input, not of the meshing. This file therefore does two
// things it would not otherwise do: it **counts the connected components** and
// says which of them reach the interface, and it **measures how close a free edge
// comes to another surface's plating**, which is what turns "the deck is not
// attached" from an inference into a number. `Section::junctionEdges` is that
// number.
//
// **The second half of the junction answer is a formulation limit rather than an
// input defect, and it would bite even on a conforming mesh.** A solid-shell
// carries its thickness as *geometry*: a node pair straddles the mid-surface along
// the surface normal, and every assumed-strain cure in `solid_shell.hpp` is
// written against that direction. Two plates meeting at an angle have two
// different thickness directions, so a shared node pair would have to point
// somewhere between them -- at a right-angled corner, 45 degrees out of both,
// which `07-fem-spike-findings.md` §6 limit 1 prices at `90 x (offset/t)^2`. So
// this file **refuses to weld across a fold**: two sub-quads share a node only
// when their normals agree to `SectionParams::foldLimit`. A junction is left open
// and reported rather than closed wrongly and not.
//
// **The price of welding it anyway is measured, on the box in the tests, and it is
// not the spurious bending stiffness that looks like the danger.** Extruding a
// corner node along the 45-degree mean normal makes the plating *thin* towards the
// corner by `cos 45`, so the section simply loses steel: on a 2 x 1 m box at one
// element per 0.5 m, welding the four corners costs **9.4% of `EA` and 8.9% of
// `EI`**, both in the unsafe direction and both a function of how many elements
// touch a corner rather than of the ship. Leaving the junction open costs neither
// -- see §2.
//
// What would close one honestly is a tie between two plates' node pairs, or a
// small element filling the corner. `constraint.hpp`'s own header already names
// the first as the missing case -- "a node on one side lands inside a face on the
// other ... that case is real and will arrive with a whole-ship mesher". It has
// arrived. It is not built here.
//
// --- 2. What the junctions cost, which is not what it looks like ----------------
//
// **For the quantity Tier 0 offers an independent answer to -- the hull girder's
// EA, its neutral axis and its EI -- an open junction costs exactly nothing, and
// that is a statement about the test rather than about the ship.** Cut a section
// at two transverse planes and prescribe the axial displacement of both cut
// sections to a plane-sections field, `u_x = eps x` or `u_x = theta x (z - z_na)`:
// every longitudinally continuous strip of plating then carries `sigma = E eps` or
// `E theta (z - z_na)` *whatever it is attached to*, because the ends alone already
// say what its strain is. The section's axial and bending stiffness come out right
// on a mesh whose plates are not joined at all -- measured on the box at a relative
// error of **6e-13 in `EA`** with its four corners cut, against 9.4% with them
// welded. The unjoined mesh is the *more* accurate one for this question.
//
// So **any validation of a section mesher that stops at `EA` and `EI` has proved
// nothing about the junctions**, and two things that do are carried here instead:
//
//   * **Torsion.** A closed cell carries torque by Bredt shear flow at
//     `4 A_enclosed^2 / integral(ds/t)` and an open one does not. On the box,
//     welded gives 0.966 of Bredt and cut gives 0.083 -- a factor of **11.6** at
//     `L = 8 m`. The factor is a function of length and that is worth knowing:
//     1.7 at 2 m, 11.6 at 8 m, 166 at 32 m, because a *short* open section is held
//     by the warping restraint of its own end planes rather than by torsion, and
//     that restraint decays. A one-bay section would have shown almost nothing.
//   * **The lowest fixed-interface frequency**, which is what Tier 1 is for and
//     which turns out to be the sharper instrument on a real ship. On the
//     reference ferry's hold between x = -7.2 m and x = 19.2 m the whole section's
//     first fixed-interface mode is **0.7785 Hz**; the decks *on their own* are
//     0.7785 Hz, to four figures, and the shell on its own is 3.4600 Hz. The
//     softest thing in the section is a 26 m deck held on two edges instead of
//     four, and adding the shell it should be welded to changes it by nothing at
//     all. That is the junctions' cost stated as a number.
//
// --- 3. Thickness seams: taper, do not stop -------------------------------------
//
// `zone.hpp` §2 stops the patch at a plate thickness seam, because a node between
// a 12 mm and a 15.5 mm strake has no single position and averaging puts a taper
// *inside* both elements. A section cannot do that: the ferry's shell crosses four
// thickness seams per side between the keel and the sheer, and a mesher that
// stopped at each would deliver eight loose girth bands instead of a section.
//
// So the seam is **tapered** across one element either side of it: each
// mid-surface node carries the area-weighted mean thickness of the sub-quads
// around it, and the element is extruded +/- t_node/2. `Section::worstTaper` is
// the resulting `dt / t` per element and `Section::taperStiffness` converts it by
// the same `90 x (offset/t)^2` the spike measured, so a section reports that its
// seam elements are several times too stiff in plate bending rather than quietly
// being so.
//
// **Several times too stiff sounds fatal and is not, and the reason is worth
// stating because it is the same reason §2 gives.** The excess is on the plate's
// bending stiffness about its own mid-surface. On the reference ferry's midship
// section that term -- every panel's own `I` about its own centroid -- is
// 0.0175 m^4 of a total 46.2 m^4, **0.038%**, so even 600% of it is 0.2% of the
// section. It matters for local plate bending and therefore for the
// fixed-interface frequencies; it does not matter for the hull girder.
//
// `ThicknessSeam::Split` is the control that measures that rather than asserting
// it: it splits the nodes at a seam so every element is exactly prismatic, at the
// cost of a section no longer joined across strakes. On a box whose flanges step
// from 10 mm to 16 mm the taper reaches `dt/t = 0.26`, which the spike's rule
// prices at 610% excess plate bending, and the two answers for `EA` differ by
// **4.5e-7**. That is the whole of what the taper costs a membrane quantity.
//
// **What the taper does cost, and it was not the plate bending, is stiffener
// steel.** `constraint::addStiffener` turns one `plateThickness` into one tie
// weight per fibre, `(e + t/2) / t`, and applies it to a node pair whose real
// separation is the *local* nodal thickness -- so a member run crossing a seam,
// where that thickness is a mean of two strakes, puts its fibres at
// `e * t_local / t_run`. Measured on the reference ferry: **47 mm** out on a 700 mm
// frame, a quarter of its Steiner term in the wrong place, and invisible in every
// aggregate. A run therefore **stops at a thickness change** the way a patch does,
// and the seam node -- whose two neighbours both differ from it -- is left in a run
// of one and dropped. That costs 5.2% of the section's member steel on the ferry
// (146.8 t against 139.1 t) and **none** of its longitudinally effective member
// area, because a longitudinal runs along the ship at constant thickness and it is
// the athwartships members that cross strakes. Missing steel is visible in a mass;
// a misplaced eccentricity is not, which is why this is the direction to fail in.
//
// The fix that costs neither is a per-station thickness in `constraint::SeamRun`,
// so that a tie weight is computed against the pair it is applied to. It is not
// built here: `constraint.{hpp,cpp}` is validated and mutation-tested against its
// current interface, and widening it is that file's change rather than this one's.
//
// --- 4. Resolution: a Tier-1 section is not a Tier-2 patch ----------------------
//
// `reduction.hpp` records a hard ceiling around 8 000 elements for Tier 2, where
// the per-element stiffness store leaves L3. **That ceiling is not this file's.**
// It is a property of an explicit solve that touches every element every step;
// a reduction touches each element once, and what it costs afterwards is set by
// the interface, not by the element count (`reduction.hpp` §7).
//
// What sets the resolution here is convergence of the *reduced* answer, and it was
// measured on the reference ferry's hold between x = -7.2 m and x = 19.2 m, eleven
// frame bays. `applyBeamLoad` prescribes the interface and relaxes the interior, so
// what it reports **is** `0.5 u_b^T K_r u_b` for the Guyan reduction of the section
// -- `reduction.hpp` property 1 -- and refining it is refining the reduced answer:
//
//     subdivision  elements  A_eff m^2   z_na m    I_eff m^4   GJ N m^2
//     ------------------------------------------------------------------
//         1          2 068    1.73122    6.85797    43.8749    3.6165e12
//         2          8 272    1.73222    6.85626    43.9042    3.5999e12
//         3         18 612    1.73255    6.85567    43.9127    3.5948e12
//         4         33 088    1.73271    6.85535    43.9165    3.5924e12
//
// **One element per panel is converged**: from `subdivision = 1` to 4 the area
// moves 0.086%, the neutral axis 2.6 mm and the second moment 0.095%, while the
// element count moves sixteenfold. The torsional stiffness -- the one quantity here
// with plate bending in it -- moves 0.67%, seven times as much, which is the same
// split §2 makes between membrane and bending questions.
//
// So a hold-sized Tier-1 section is **2 068 elements**, a quarter of the Tier-2
// ceiling, and the ceiling was never the binding constraint. What binds is the
// interface: a transverse cut of this ferry is 195 panel corners and therefore 780
// nodes and 2 340 boundary DOF at `subdivision = 1` -- the whole of the reduced
// model at zero modes, and 93% of it with the 178 the default cutoff asks for --
// and `Psi^T M Psi` is `O(n_i n_b^2)`. Refining costs sixteen times per doubling, not
// four -- and it buys 0.05%.
//
// **One warning about `ReduceParams` at this size, because the default is wrong
// here and says nothing about it.** The hold's lowest fixed-interface frequency is
// 0.78 Hz (§2), so the 20 Hz default cutoff asks for **178 modes**; that takes 275 s
// and the subspace iteration does not converge in its 60 iterations. Guyan alone --
// `modes = 0` -- is 6.2 s and is exactly right at the interface for any static load.
// A substructure softer than the ship it is part of is a sign that the section is
// not yet a component, not a reason to keep more modes.
//
// SI units, body frame per CLAUDE.md.
#pragma once

#include "constraint.hpp"
#include "reduction.hpp"
#include "scantlings.hpp"
#include "solid_shell.hpp"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace sim::section {

// What happens where the plating changes thickness. See §3.
enum class ThicknessSeam {
    // One element either side of the seam is tapered from one thickness to the
    // other. The section stays joined; the seam elements are too stiff in plate
    // bending by `Section::taperStiffness`.
    Taper,
    // The nodes at the seam are split, so every element is exactly prismatic and
    // the section comes apart into constant-thickness bands. The control §3 uses
    // to measure what the taper costs, and never the right answer for a section
    // that has to carry shear across a strake boundary.
    Split,
};

struct SectionParams {
    // The two transverse cut planes, m in the body frame. `xFrom < xTo`.
    //
    // **Cut on a frame station.** Panels are bounded by frames, so a cut there
    // passes between panels and takes nothing with it; the reference ferry's
    // watertight bulkheads are at x = -44, -38, -8, 20 and 44, none of which is a
    // multiple of the 2.4 m frame spacing, and cutting at one of them puts the
    // plane through 188 panels. Those panels are **excluded** and counted in
    // `Section::straddlingPanels`, which leaves a hole in the section rather than
    // a short bay -- so a cut that is not on a seam is a request this file
    // declines to answer quietly.
    double xFrom = 0.0, xTo = 0.0;

    // Elements across each panel edge, both ways. An integer and the same both
    // ways for the same reason `zone::MeshParams::subdivision` is: a shared edge
    // is then divided identically by the panels either side of it whichever way
    // round each numbers its corners, so the mesh conforms with no seam-matching
    // pass.
    //
    // One, not four. A Tier-1 section is reduced, so its mesh serves stiffness and
    // inertia rather than tearing, and §4's measurement says one element per
    // stiffener bay is where the section's own stiffness already is.
    int subdivision = 1;

    // Which plating to take. All three by default: the deck plating alone is 44%
    // of the reference ferry's midship second moment and the shell 32%, so a
    // section missing either is not a hull girder.
    bool shell = true, deck = true, bulkhead = true;

    // Build the stiffener fibres and the `reduction::Attachment` that carries
    // them. Off is the negative control §2 of `reduction.hpp` describes: bare
    // plating, measurably softer by what `scantlings::hullGirderSection` says the
    // members are worth.
    bool members = true;

    ThicknessSeam thicknessSeam = ThicknessSeam::Taper;

    // Two sub-quads share a node when their positions agree to this and their
    // normals agree to `foldLimit`. Matches `zone::MeshParams::weldTolerance`.
    double weldTolerance = 1e-6;  // m

    // Beyond this fold two sub-quads are different surfaces and are **not**
    // welded, however close their nodes are. See §1: a shared node pair has one
    // thickness direction and two plates at an angle have two.
    double foldLimit = 0.785;  // rad, 45 degrees

    // How near a cut plane a mid-surface node must be to be on it. The mesh's
    // nodes are the mid-surface offset by +/- t/2 along the surface normal, so a
    // node of a panel whose normal leans out of the transverse plane is up to t/2
    // off it -- 8 mm on this ship. The interface is therefore chosen on the
    // **mid-surface** and both of each pair's extruded nodes taken, which is exact;
    // `reduction::nodesNearPlanes` at its 1e-9 default would silently keep one of
    // each pair and leave the interface with half its nodes.
    double planeTolerance = 1e-6;  // m

    // Reported, not enforced, exactly as in `zone::MeshParams`.
    double normalSpreadWarning = 0.05;  // rad
    // Above this `dt / t` across one element the taper of §3 is worth a word.
    double taperWarning = 0.05;

    // How far from another surface's plating a free edge has to be before it is
    // *not* an unwelded junction. A deck edge clipped to the hull sits on the
    // shell to rounding, so anything within a plate thickness or so of another
    // surface is a joint that this file could not make. Scaled off the plating
    // rather than absolute would be better and needs a thickness this does not
    // have until the mesh exists; 25 mm is more than the thickest plate on the
    // reference ship, whose bilge strake is 15.5 mm.
    double junctionTolerance = 0.025;  // m
};

// A meshed section: the elements, what they came from, how well they joined up,
// and everything `reduction::Substructure` needs to consume it.
struct Section {
    solidshell::HexMesh mesh;

    // The members, as eccentric fibres tied to the plating by `constraint.hpp`,
    // and the same thing packaged as `reduction.hpp` §8 takes it. Both are filled;
    // `attachment` is derived from `stiffening` and is what a `Substructure` wants.
    constraint::Stiffening stiffening;
    reduction::Attachment attachment;

    // The interface: every node on either cut plane, ascending and unique. This is
    // what `reduction::Substructure`'s third argument takes.
    std::vector<std::uint32_t> interfaceNodes;
    // The two planes separately, because a section is loaded by moving one
    // relative to the other and a caller needs to know which is which.
    std::vector<std::uint32_t> aftNodes, forwardNodes;

    std::vector<int>          panelOf;      // per element -> StructuralMesh::panels
    std::vector<double>       elementArea;  // per element, m^2 of mid-surface
    std::vector<int>          panels;       // the distinct panels covered, ascending
    // Per node: 1 on the +normal face, 0 on the -normal face, as `zone::Patch`.
    std::vector<std::uint8_t> outerFace;
    // Per node, m: the plate thickness the element extrusion used there. A seam
    // node under `ThicknessSeam::Taper` carries the mean of the strakes it joins.
    std::vector<double>       nodeThickness;

    StructuralMaterial material;
    double xFrom = 0, xTo = 0;
    double length() const { return xTo - xFrom; }

    double area = 0;        // m^2 of mid-surface meshed
    double plateMass = 0;   // kg of plating
    double memberMass = 0;  // kg of stiffener
    double mass() const { return plateMass + memberMass; }

    // --- Geometry quality, the same instruments `zone::Patch` carries -----------
    double worstNormalSpread = 0;  // rad
    double spuriousStiffness = 0;  // 90 * spread^2
    int    distortedElements = 0;
    double worstJacobian = 0;      // <= 0 is fatal
    double worstAspect = 0;
    // The DOF half-bandwidth the node numbering delivers. `solidshell::solveStatic`
    // numbers its free degrees of freedom in the mesh's own order and has no
    // renumbering pass, so this is what a static solve of the section costs -- and it
    // is not a rounding: the ferry hold measures 146 here against 1 382 before
    // `reduction::bandwidthReducingOrder` was one of the candidates, which is 5.3 s
    // of banded factorisation against 0.14.
    std::size_t halfBandwidth = 0;

    // --- Thickness seams, see §3 -----------------------------------------------
    int    taperedElements = 0;    // elements whose corners do not share a thickness
    double worstTaper = 0;         // max dt / t_mean over an element
    double taperStiffness = 0;     // 90 * worstTaper^2, the excess plate bending

    // --- Topology, see §1. This is what the file exists to report ---------------
    int    surfaces = 0;    // maximal sets of sub-quads welded to one another
    int    components = 0;  // connected components of the assembled element mesh
    // Components touching neither cut plane. Each is a mechanism in `K_ii` and
    // `reduction::Substructure` will *not* catch it -- its interface check is
    // geometric and a floating component leaves a tiny positive pivot rather than
    // a zero one, which is the failure `reduction.hpp` §3 records measuring.
    int    floatingComponents = 0;
    // Components reaching *both* cut planes. Only these carry any of the section's
    // axial or bending stiffness; one held at a single end is a cantilever hanging
    // off the interface and contributes nothing a beam would recognise.
    int    spanningComponents = 0;
    // Per node, which component it belongs to, or -1 for a node no element reached.
    // Exposed because every rigid body motion is per component, so a caller
    // restraining a disconnected section has to restrain each piece -- see
    // `applyBeamLoad`, where getting this wrong made the whole factorisation fail.
    std::vector<int> componentOf;
    // Element edges used once that are not on a cut plane: the section's own free
    // boundary, plus every junction that did not weld.
    double freeEdgeLength = 0;  // m
    // How much of that runs within `SectionParams::junctionTolerance` of another
    // surface's plating -- an edge that is sitting on a plate it is not joined to.
    // Zero is a section whose only free edges are genuinely free.
    double junctionEdges = 0;   // m
    double worstJunctionGap = 0;  // m, the furthest such an edge was from it

    int straddlingPanels = 0;   // panels a cut plane passed through, and dropped
    int membersAttached = 0;    // members that contributed fibres
    int membersRefused = 0;     // members whose web is not along the plate normal
    int membersMissed = 0;      // members that lie on no run of mesh nodes
    // Times a member's run of stations was broken because the plating under it
    // changes thickness. `constraint::addStiffener` turns one `plateThickness` into
    // one tie weight per fibre, and a tie weight is only right for the pair
    // separation it was computed against, so a run that crossed a strake seam would
    // put its fibres at `e * t_local / t_run` instead of at `e` -- 47 mm out on the
    // reference ferry's frames. A stiffener is therefore discontinuous at a
    // thickness seam, which loses one segment of it per seam crossed.
    int memberRunsSplitByThickness = 0;
    // The **longitudinally effective** member area an average transverse cut of this
    // section sees, m^2 -- each member's profile area weighted by the fraction of the
    // section's length it spans, so a member covering one bay of eleven counts for a
    // eleventh. Only members with an extent along x are counted, because an
    // athwartships one carries no longitudinal stress and `hullGirderSection` drops
    // it for the same reason.
    //
    // The two sum to what `sectionElements` reports for the stiffeners, so a section
    // coming out short against `hullGirderSection` is an **accounting** rather than a
    // discrepancy. On the reference ferry the shortfall is the three girders, which
    // sit off the longitudinal spacing and therefore pass through no node of a mesh
    // built from panel seams.
    double attachedMemberArea = 0;
    double missedMemberArea = 0;

    std::vector<std::string> problems;

    std::size_t elementCount() const { return mesh.elementCount(); }
    std::size_t nodeCount() const { return mesh.nodeCount(); }
    bool empty() const { return mesh.elementCount() == 0; }
};

// Mesh the section. Always returns something: an impossible request yields an
// empty section and a full account in `problems`, in the same spirit as
// `makeStructuralMesh` and `zone::buildPatch`.
Section buildSection(const StructuralMesh& structure, const SectionParams& params);

// --- Loading a section like a beam ----------------------------------------------
//
// A section's whole claim on Tier 0 is that it carries the same `EA`, neutral axis
// and `EI` that `hullGirderSection` computes from the same scantlings by a
// completely different route -- summing `A`, `A z` and `A z^2` over a transverse
// cut. Asking a mesh that question is a static problem with a prescribed boundary,
// and getting it right is fiddly enough in three places (a test, a tool, a caller)
// that it lives here once.
//
// The field is plane sections: `u_x = strain * x + curvature * x * (z - reference)`
// on both cut planes, with `u_y` and `u_z` left **free** there. Free is the whole
// point -- pinning them would suppress the Poisson contraction and report a
// section stiffer than it is by roughly `1 / (1 - nu^2)`, and on a section shorter
// than it is deep that error is not local, it is everything.
struct BeamLoad {
    double strain = 0;     // axial strain, dimensionless
    double curvature = 0;  // 1/m, hogging positive: tension above `reference`
    double reference = 0;  // m above the baseline that the curvature turns about
};

// What a prescribed beam field costs the section, and therefore what its section
// properties are.
//
// `axialForce` and `bendingMoment` are the resultants on the **forward** cut
// plane, taken from the reaction there. With a pure `strain` they give
// `EA = axialForce / strain` and, if `reference` is the true neutral axis, a
// bending moment of zero -- so the neutral axis is *found* rather than assumed,
// and comparing it against `HullGirderSection::neutralAxis` is a third
// independent number rather than a restatement of the first two.
struct BeamResponse {
    bool ok = false;
    double axialForce = 0;     // N through the section
    double bendingMoment = 0;  // N m about `BeamLoad::reference`, hogging positive
    double strainEnergy = 0;   // J stored, plating and fibres together

    double axialStiffness = 0;    // N, `axialForce / strain`: this is `EA`
    double bendingStiffness = 0;  // N m^2, `bendingMoment / curvature`: this is `EI`
    double peakDisplacement = 0;  // m

    // `max |K u|` over the degrees of freedom nothing holds. Zero to the solver's
    // conditioning on a solve that converged, so it measures the solve rather than
    // decorating it -- and it is the only thing here that would notice a stiffness
    // assembled one way and a reaction taken another.
    double residual = 0;  // N
    // What the three rigid-body restraints carry. **Exactly zero** is the right
    // answer and not a tolerance: they select one of a family of zero-energy
    // motions rather than resisting anything, so a reading here is a defect.
    double restraintReaction = 0;  // N

    std::string problem;
};

// Solve it. `material` is the section's own unless a caller wants to sweep E.
//
// The interface is prescribed and the interior is free, so this is exactly the
// static condensation `reduction.hpp` property 1 calls exact: the strain energy it
// reports **is** `0.5 u_b^T K_r u_b` for the Guyan reduction of this section, and
// `tests/test_section.cpp` asserts that identity against a real `craigBampton`
// rather than taking it on the header's word.
BeamResponse applyBeamLoad(const Section& section, const StructuralMaterial& material,
                           const BeamLoad& load);

// The same, with the section twisted: the forward plane is rotated by `twist`
// radians about the x axis through (`y`, `z`) = (0, `reference`) and the aft plane
// held, with **every** interface degree of freedom prescribed -- a rigid disc,
// because a twist is a statement about the section's shape and leaving it free
// would let the section shear instead of turning.
//
// `torsionalStiffness` is `torque * length / twist`, which is `GJ` for a section
// long enough that warping is free. It is here because it is the cheapest question
// that can tell a welded section from an unwelded one: a closed cell carries
// torsion by Bredt shear flow at `4 A^2 / integral(ds/t)` and an open one by
// `sum s t^3 / 3`, and the two differ by two to three orders of magnitude. See §2.
struct TorsionResponse {
    bool ok = false;
    double torque = 0;               // N m about x
    double strainEnergy = 0;         // J
    double torsionalStiffness = 0;   // N m^2, `GJ`
    std::string problem;
};

TorsionResponse applyTwist(const Section& section, const StructuralMaterial& material,
                           double twist, double reference = 0.0);

}  // namespace sim::section
