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
// What closes one honestly is a tie between two plates' node pairs, and **that is
// now built -- see §5**, which is where the numbers this section quotes for an
// unjoined section stop being what the file delivers and start being its negative
// control.
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
//     reference ferry's hold between x = -7.2 m and x = 19.2 m the whole *unjoined*
//     section's first fixed-interface mode is **0.7785 Hz**; the decks *on their
//     own* are 0.7785 Hz, to four figures, and the shell on its own is 1.6999 Hz.
//     The softest thing in the section is a 26 m deck held on two edges instead of
//     four, and adding the shell it should be welded to changes it by nothing at
//     all. That is the junctions' cost stated as a number, and §5 is what moves it:
//     tied, the same section is 2.3026 Hz.
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
//         1          2 068    1.72945    6.86238    43.8169    3.6164e12
//         2          8 272    1.73081    6.85933    43.8598    3.5998e12
//         3         18 612    1.73126    6.85828    43.8728    3.5947e12
//         4         33 088    1.73148    6.85777    43.8787    3.5923e12
//
// **The four columns above are what `tools/section_probe --sweep=4 --no-reduce`
// prints, and the ones they replace were not.** Before the halo the same tool gave
// 1.73075, 6.85784, 43.86283 at `subdivision = 1` against a published 1.73122,
// 6.85797, 43.8749 -- a table that had drifted from the program it names and that no
// gate re-ran, which is the failure CLAUDE.md records four times over. The halo moved
// it again, by rather less than the drift: the hold's forward plane at x = 19.2 is a
// strake seam, and a section that can now see it stops four member runs there.
//
// **One element per panel is converged**: from `subdivision = 1` to 4 the area
// moves 0.12%, the neutral axis 4.6 mm and the second moment 0.14%, while the
// element count moves sixteenfold. The torsional stiffness -- the one quantity here
// with plate bending in it -- moves 0.67%, five times as much, which is the same
// split §2 makes between membrane and bending questions. The first three move about
// half again as much as they did before the halo, and for a reason worth stating: a
// thickness seam costs a **fixed** two dropped fibre segments per member per seam
// (§3), so refining the mesh dilutes it, and a coarse mesh is where it is dearest.
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
// --- 5. The junction tie: joining the plates without welding them ----------------
//
// A free-edge node lying on another surface is **tied** to the point of that
// surface it lands in. Three constraints per extruded node, one per axis, over the
// eight mesh nodes of the master face:
//
//     u[slave] = sum_a N_a(xi, eta) * [ (1-w) u[bottom_a] + w u[top_a] ]
//
// -- bilinear in the master face's two in-plane coordinates, and through its
// thickness by exactly the weight `constraint.hpp` §1 derives. The slave keeps no
// unknown of its own: `solidshell::solveStatic` and `reduction::Substructure` both
// scatter through the transformation, so the matrix they factor **is** `T^T K T`
// and there is no penalty stiffness anywhere.
//
// **The tie is not a weld and the difference is the whole point.** A weld merges
// two node pairs into one, which has one thickness direction where two plates at an
// angle need two, and pays for it in steel -- §1's 9.4% of `EA`. A tie moves no
// node and merges nothing: the mesh it produces has the same nodes in the same
// places, the same elements, the same mass, still four surfaces on the box, and one
// connected component. `surfaces` counts what shares nodes; `components` counts what
// is joined; before this the two were the same question.
//
// **A `DofBlock` cannot express it, and finding out why is the load-bearing part.**
// `constraint.hpp` eliminates a stiffener fibre with a `Tie` and hands the result to
// `reduction::Attachment::stiffness`, and that works because a fibre's endpoints are
// *not mesh nodes*: they have no rows, so `T^T K T` over the masters is the whole of
// what the fibre contributes and adding a block is exact. A junction ties a node
// that **is** a mesh node, with elements of its own. Eliminating it means rewriting
// rows that already exist, and no amount of added stiffness redirects a row. So the
// constraint is carried as data -- `solidshell::Mpc`, in
// `reduction::Attachment::constrained` -- and applied by the assembler.
//
// --- What it does, on the box, where every quantity has a closed form -------------
//
//     corners      components   EA        EI        GJ / Bredt
//     ---------------------------------------------------------
//     cut               4       1.000000  1.000000    0.083
//     welded            1       0.905788  0.911115    0.966
//     tied              1       1.000000  1.002158    1.099
//
// The tied mesh is the only one that is **both** exact in `EA` and closed in
// torsion. `1.099` is not an overshoot of Bredt: a closed cell does not stop having
// the open section's own `sum s t^3 / 3`, and `1.000 + 0.083` is what the two
// together predict.
//
// `EI` costs 0.216% at one element per panel and that is a **consistency error, not
// the formulation**: at a butt corner the tied node sits half a plate thickness past
// the end of the other plate's mid-surface, so the tie extrapolates a bilinear
// approximation to a quadratic transverse contraction. It falls 2.16e-3, 2.85e-4,
// 8.62e-5 over a threefold refinement. Clamping the tie onto the face instead would
// trade it for an `O(t)` error that refinement cannot reach, which is why the
// overshoot is bounded rather than forbidden -- `SectionParams::junctionOvershoot`.
//
// --- And on the ferry, where the fixed-interface frequency is the instrument ------
//
// The hold between x = -7.2 m and x = 19.2 m, at one element per panel:
//
//     untied: 7 components, 369.6 m of free edge lying on plating it is not joined
//             to, first fixed-interface mode **0.7785 Hz** -- which is the decks'
//             own 0.7785 Hz to four figures, because the shell contributes nothing.
//     tied:   1 component, 309.6 m of that edge joined, first fixed-interface mode
//             **2.3026 Hz** -- above the decks' 0.7785 *and* the shell's 1.6999,
//             which is what a joined structure does.
//
// `GJ` goes 3.6164e12 to 5.2387e12 N m^2, +44.9%. `A_eff` moves +0.19% and `I_eff`
// +0.12%, which is the same local consistency error the box shows and is why §2's
// warning stands: **neither of those two could have told you the junction was
// open.**
//
// The 60 m that stays open is the nodes on a cut plane, and it is open because a
// master *face* spans x: half its nodes are interior to the section, and a boundary
// degree of freedom written as a function of an interior one is not one a reduction
// keeps exactly -- `reduction::Substructure` refuses it. On an eleven-bay section that
// is two stations of twelve.
//
// **§9 closes it, for a chain, with a tie to the master's line rather than its face.**
// A lone section's cut planes stay open and should: they are what a load is prescribed
// on, and a prescribed degree of freedom cannot also be derived. The 60 m above is
// therefore still what this section delivers, and it is exactly
// `planeTiedEdgesAft + planeTiedEdgesForward` -- 28.8 + 31.2 -- which is what a
// neighbour on either side would join.
//
// --- What it costs, which is the band ---------------------------------------------
//
// A tie couples nodes that no element edge joins, and the numbering has to know:
//
//     ferry hold        DOF half-bandwidth   banded solve
//     ------------------------------------------------------
//     untied                    146              0.16 s
//     tied, ordering blind    10 769            (dense)
//     tied, ordering aware     1 520              5.34 s
//
// The ordering is chosen by comparing three candidates, and it is scored on the
// *tied* graph -- element edges plus tie edges -- because an ordering scored on the
// sub-quads alone is scored on a graph the solver does not have. That is worth a
// factor of seven. What is left is real: joining a section's decks to its shell
// closes the cross-section into a tube with internal webs, and a wavefront that used
// to sweep one flat sheet now has to cross a deck. `Substructure`'s interior band
// goes 125 to 455 and its Craig-Bampton 6.0 s to 57.7 s. At subdivision 2 the banded
// static solve goes 1.1 s to 149 s, which is why `tools/section_probe` runs its
// resolution sweep untied.
//
// --- What the tie cannot do ------------------------------------------------------
//
//   * **Chains are refused, not composed.** A node whose master is itself a slave --
//     three surfaces meeting at a point -- is left untied and counted, because a
//     chain resolved silently is a modelling error that assembles.
//   * **The weights are a partition of unity but not always positive.** A face
//     overshoot `d` puts `-d/2` on a master, and a real junction's 9 mm gap makes the
//     through-thickness weight 1.69 rather than 0.5. The condensed mass is `T^T M T`
//     row-summed, so a negative weight is a negative share of the slave's steel;
//     `reduction::Substructure` refuses a non-positive nodal mass and is the backstop.
//     Unlike `constraint.hpp`'s eccentric fibre there is nothing to give up here --
//     the fibre *extrapolates* and has to abandon its first moment, while a junction
//     interpolates and keeps it.
//   * **The moment transfer is only as good as the master element.** A plate's
//     rotation about the junction line is carried into the other plate's *in-plane*
//     displacement gradient, which for a bilinear element is the element's own
//     average slope. Membrane continuity -- the shear flow torsion needs -- is
//     transmitted exactly; the edge moment is transmitted to `O(h)`.
//
// --- What mutation testing left standing ------------------------------------------
//
// Thirty-two single edits across this file, `solid_shell.cpp` and `reduction.cpp`;
// thirty are killed. **The two that are not are equivalent on every input this
// repository has, and saying which is more useful than a score:**
//
//   * **Deleting the chain refusal changes nothing**, because no input here produces
//     a chain: the "already a master" skip absorbs the second side of every junction,
//     on the box, on a box with a deck through it, on one with a bracket at the
//     corner, and on the ferry. The check is kept because the *consequence* of a
//     chain is not local -- `solidshell::DofExpansion` refuses the whole set and the
//     section stops reducing -- and a mesher that would rather say "this junction is
//     open" than emit one is the honest half of that. Both refusals are tested
//     directly, in `tests/test_solid_shell.cpp` and `tests/test_reduction.cpp`.
//   * **Running the projection for one iteration instead of twenty-four changes
//     nothing**, to the last digit of `|sum w X - X_slave|` on the ferry's own warped
//     plating. Gauss-Newton from the face centre lands in one step on a face that is
//     flat, and every sub-quad here is flat enough. The loop is kept because a
//     genuinely warped master face is a property of a *finer* mesher, not of this
//     one, and the cost is twenty-three iterations of a 2x2 solve per junction node.
//
// The run also found something older than any of this: `reduction::Substructure`'s
// `totalMass()` and `stiffnessTimes()` read past the end of arrays the constructor
// never sized when the substructure **refused**, which was reachable from an
// inverted element long before a constraint existed. It surfaced as a segmentation
// fault three tests downstream of the mutant that provoked it -- which is why the
// mutation harness has to look at the *exit code* and not only at the failure lines,
// something its first version got wrong and scored eight false survivors on.
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

    // Average the nodal normal and the nodal thickness over a bay of plating
    // **beyond** each cut plane as well, and mesh only what is inside. See §8: it is
    // what makes a node's position a property of the ship rather than of the window,
    // and it is what lets a chain match its interface at `matchBoundaries`' own
    // default instead of at a tolerance chosen to cover the disagreement.
    //
    // Off is the negative control §8's measurements are made against, and it is what
    // this file did before the halo existed. It is never the right answer for a
    // section that is going to be one link of a chain.
    bool halo = true;

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

    // Tie every free-edge node that lands on another surface to that surface. See
    // §5. Off is the negative control the junction measurements are made against,
    // and it is what this file did before the tie existed.
    bool junctions = true;

    // Build the **in-plane line ties** of §9 as well: the junction nodes that lie on a
    // cut plane, tied to the line the other surface draws on that same plane rather
    // than to a face that leaves it. A section does not apply them -- its interface is
    // prescribed -- it carries them in `Section::planeTies` and `buildChain` puts them
    // on the assembled model, where an interior plane is shared instead.
    //
    // Off is the negative control every figure in §9 is measured against, and it is
    // what this file did before: `Section::junctionsOnInterface` then counts every
    // junction node on a cut plane instead of only the ones no line reaches. It has no
    // effect at all on a section solved on its own.
    bool interfaceTies = true;

    // How far outside a master face, in the face's own natural coordinates, a tied
    // point may land. Zero would refuse every L-junction: two plates butting at a
    // corner have each other's mid-surface *ending* at the corner line, so a node
    // half a plate thickness the other side of it is outside the last face by
    // `t / (2 h)` -- 2% of a half-metre element on 10 mm plating. It is small and
    // it is not zero, and clamping it instead would tie a node at `z` to the field
    // at `z + t/2`, which loses exactly the lever arm `EI` is made of.
    //
    // The bound matters because a tie that extrapolates has a **negative weight**,
    // and a negative weight means a negative share of the slave's mass on that
    // master. The worst negative weight of a bilinear face at overshoot `d` is
    // `-d/2`, so this is also the guarantee that the condensed mass stays positive
    // -- `reduction::Substructure` checks it and would refuse the section rather
    // than integrate it.
    double junctionOvershoot = 0.25;

    // The same guarantee in the *other* direction, which `junctionOvershoot` does
    // not cover and which is what actually bites on a real hull.
    //
    // A tie splits the slave through the master's thickness by
    // `constraint::tieWeight`, `(e + t/2) / t` for a slave `e` off the master's
    // mid-surface. Inside the plate that is a weight in [0, 1]; outside it one of
    // the two goes negative, and a real junction is always outside -- a deck edge
    // sits on the shell's *outer* face, plus whatever gap two plates of different
    // thickness leave. The reference ferry's ordinary deck-to-shell junctions run
    // at `w = 1.69`, so `-0.69` on the inner face, and that is not a tolerance to
    // be tightened: it is where the steel is.
    //
    // What is not survivable is a master face asked to give up **more of the
    // slave's steel than the slave has**. `reduction::Substructure` row-sums
    // `T^T M T` and refuses a non-positive nodal mass, so past one full share the
    // tie does not make the section slightly wrong, it makes it unusable -- and it
    // refuses the whole section rather than the one junction. One share is
    // therefore the limit, and a tie past it is left open and counted in
    // `Section::junctionsThroughThickness`, the same way an over-reaching face
    // coordinate is. Measured on the reference ferry: every junction on the ship
    // is inside it but six, at x = 45.6..50.4 m, where a deck lands 2.15
    // thicknesses off a bilge strake's mid-surface and asks for `w = 2.65`.
    double junctionWeightLimit = 1.0;
};

// One in-plane line tie: a junction node on a cut plane written as a linear function
// of the four mesh nodes of the segment the other surface draws on that same plane.
// See §9.
//
// **Every degree of freedom in it lies on `plane`**, which is the whole point and the
// only thing that distinguishes it from a `§5` tie. That makes it a relation among the
// section's *boundary* unknowns, so it is meaningless inside a `reduction::Substructure`
// -- which keeps its interface exactly and refuses a boundary degree of freedom that is
// a function of others -- and exact on the assembled model of a chain, where the plane
// is shared by the two sections that meet on it and prescribed by neither.
struct PlaneTie {
    // Which cut plane: 1 the aft, 2 the forward. From the mesher's own membership test
    // on the **mid-surface**, not from a position test on the node -- a node whose
    // plating leans out of the transverse plane is up to `t/2` off it, which is the
    // distinction §6 note 1 records costing a day.
    std::uint8_t plane = 0;
    // In the section's own mesh degrees of freedom. `buildChain` maps them through
    // `reduction::Substructure::boundaryDof` and `Assembly::from`.
    solidshell::Mpc mpc;
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
    // Smallest `det J` over the centre and the 2x2x2 rule of every element -- the
    // places `solidshell` integrates and the only places it requires positivity.
    // <= 0 is fatal. **Not** the smallest nodal determinant, which is zero by
    // construction on the collapsed elements of §7 and says nothing about them.
    double worstJacobian = 0;
    // Elements with two coincident nodes: the wedge a degenerate plate panel
    // extrudes to. Sound, integrated, and counted rather than refused -- see §7.
    int    collapsedElements = 0;
    // Elements that are genuinely folded. Non-zero means the section is worthless.
    int    invertedElements = 0;
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

    // --- The junction tie, see §5 ------------------------------------------------
    //
    // How much of `junctionEdges` is joined after all: an edge both of whose nodes
    // were tied into the surface they were lying on. `junctionEdges - tiedEdges` is
    // what is still carrying no shear.
    double tiedEdges = 0;  // m
    int    junctionTies = 0;      // mid-surface nodes tied, so `2 x` mesh nodes
    // Free-edge nodes on another surface that could **not** be tied, and why. A
    // node whose master would itself be a slave is refused rather than chained: see
    // `solidshell::DofExpansion`. So is one on a cut plane -- the interface is
    // prescribed and a prescribed degree of freedom cannot also be a function of
    // others -- and one whose master face is further outside than
    // `SectionParams::junctionOvershoot`.
    int    junctionsChained = 0, junctionsOnInterface = 0, junctionsOutsideFace = 0;
    // And one whose through-thickness split would hand a master face more than
    // `SectionParams::junctionWeightLimit` shares of the slave's mass, negative.
    int    junctionsThroughThickness = 0;
    // The largest overshoot outside a master face, in that face's own natural
    // coordinates, and the largest through-thickness weight any tie used. Both
    // bound how negative a master's share of the slave's mass can be, which is the
    // one way this formulation can produce something an integrator cannot use.
    double worstJunctionOvershoot = 0;
    double worstJunctionWeight = 0;

    // --- The in-plane line tie, see §9 ---------------------------------------------
    //
    // What a cut plane's ring of junctions comes to when it is not left open: the
    // constraints, in this section's mesh degrees of freedom, and what applying them
    // would join. A section on its own applies none of them -- `buildChain` does, and
    // only at the interior planes, because an outermost plane is still prescribed.
    std::vector<PlaneTie> planeTies;
    int planeTieNodes = 0;  // mid-surface nodes tied this way, so `6 x` constraints
    // How much more of `junctionEdges` is joined once the ties on the aft plane, on
    // the forward plane, or on both at once are applied. `Both` is an edge reaching
    // from one cut to the other -- a section one element long -- which needs both
    // planes to be interior before it is joined at all, so it is counted apart rather
    // than credited to either.
    double planeTiedEdgesAft = 0, planeTiedEdgesForward = 0, planeTiedEdgesBoth = 0;
    // Junction nodes on a cut plane a line tie could not take, and why. `Unreached` is
    // **no one segment on this plane reaching both of the slave's extruded nodes**:
    // either no other surface draws one within `junctionTolerance`, which is the honest
    // answer where the other plate stops short of the cut, or the slave straddles the
    // node between two of them and each half of its pair falls the other way. Tying the
    // two halves to two different segments would split the slave against two plates.
    int planeTiesUnreached = 0, planeTiesOutsideLine = 0;
    int planeTiesThroughThickness = 0, planeTiesChained = 0;
    // The same two bounds §5 reports, over the line ties, and one a face tie has no
    // analogue of.
    //
    // **`worstPlaneTieSlip` is the whole of what a line approximates and a face does
    // not.** Projecting onto a face leaves a residual along the master's own normal
    // and nothing else, and the through-thickness weight is exactly what that residual
    // becomes -- so a face tie drops nothing. Projecting onto a line leaves that, plus
    // whatever ran along the master's surface *across* the line; on a cut plane that
    // direction is x, so it is the slave's own normal leaning fore or aft, bounded by
    // `t/2`. Zero on prismatic plating, m elsewhere.
    double worstPlaneTieOvershoot = 0, worstPlaneTieWeight = 0, worstPlaneTieSlip = 0;

    int straddlingPanels = 0;   // panels a cut plane passed through, and dropped
    // Panels outside the section that reached one of its nodes, were averaged into
    // that node's normal and thickness, and were then thrown away -- the halo of §8.
    // **Reached, not looked at**: the halo is drawn from every panel whose extent
    // touches the window and the weld decides which of them land on anything, so this
    // is the count that says the halo did something. Zero for a section that contains
    // the whole of every surface it touches -- and an invariance test on a window with
    // no halo is a test of nothing, which is why it is reported at all.
    int haloPanels = 0;
    // Members with no extent along x -- frames, deck beams -- lying on the forward
    // cut plane, left to the section forward of it. A station is a cut plane for the
    // sections on **both** sides of it, so a member sitting on one belongs to
    // whichever claims it and not to both: taken by both, a chain carries that
    // frame's steel and stiffness twice, which is 10.1% of the stiffener mass of two
    // four-bay sections of the reference ferry and is invisible in any total but the
    // one it inflates. The rule is the half-open one `SectionParams::xFrom` already
    // describes for transverse plating.
    int membersOnForwardPlane = 0;
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
    // The same two, as first and second moments about the **baseline** -- `A z` and
    // `A z^2 + I_own`, weighted by the same length fraction, with the profile's own
    // second moment rotated in by the web's direction cosines exactly as
    // `sectionElements` does it.
    //
    // **An area alone cannot correct an `I` comparison, and reporting only the area
    // made it look as though the two tiers disagreed by 5% amidships when they agree
    // to 0.4%.** The girders this mesh cannot reach are 2.5 m^4 of the ferry's 46.2 --
    // 5.3% of her second moment against 4.4% of her area, because they sit low in a
    // double bottom and a second moment is a lever arm squared. Subtracting the area
    // and not the moment is subtracting the wrong one of the two.
    //
    // About the baseline rather than about a neutral axis, because the axis a caller
    // wants to compare on is `hullGirderSection`'s and not this section's: a moment
    // about z = 0 shifts to any axis by the parallel axis theorem and a moment about
    // some other axis does not shift back without the area it was taken with.
    double attachedMemberFirstMoment = 0, missedMemberFirstMoment = 0;    // m^3
    double attachedMemberSecondMoment = 0, missedMemberSecondMoment = 0;  // m^4

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
// `displacementOut`, when given, receives the solved field over the section's own
// mesh degrees of freedom -- with the junction ties' eliminated rows already filled
// in from their masters, because a zero left at one of those is the hole
// `reduction::recover` shipped with and it reads as an enormous artificial gradient
// in every element that touches the node. It is what `axialStress` below consumes.
BeamResponse applyBeamLoad(const Section& section, const StructuralMaterial& material,
                           const BeamLoad& load, std::vector<double>* displacementOut = nullptr);

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

// --- 10. The stress the section carries, which is where the two tiers can differ ---
//
// `BeamResponse` reports **resultants** -- a force and a moment on the forward plane
// -- and those are exactly the quantities `hullGirderSection` reaches by summing over
// a transverse cut. That is what makes them a real comparison, and it is also what
// makes them a limited one: two models can agree on `EA`, on the neutral axis and on
// `EI` to ten figures and still disagree everywhere about the stress, because a
// resultant is an integral and an integral does not see how its integrand is spread.
//
// A beam idealisation can produce exactly one distribution:
//
//     sigma_xx = N / A + M (z - z_na) / I
//
// -- linear in `z`, constant at a given `z` across the whole breadth, and with no
// memory of how the load reached the section. A reduced 3D model produces whatever
// the plating does, and the three ways it differs -- shear lag across a wide deck,
// the warping of a section that is not closed, and the transverse restraint a ring
// of frames adds -- are the whole of what this tier is for. **None of the three is
// visible in `EA` or `EI`**, which is the same warning §2 makes about the junctions,
// arriving at the other end of the model.
//
// So the stress is exposed where the element carries it: the eight Gauss points of
// every element a transverse plane passes through, each with the volume it
// integrates. **The volume is the weight and not a decoration** -- a deck element
// spans metres of this ferry and a bilge strake spans centimetres, so an unweighted
// mean over Gauss points is an average over the mesher's resolution rather than over
// the steel.
struct StressSample {
    Vec3 at{};           // m, body frame: where the element integrates
    double volume = 0;   // m^3 this point carries
    double sigmaXX = 0;  // Pa, longitudinal, tension positive
};

// Every Gauss point of every element the transverse plane `x` passes through, by the
// same half-open rule `sectionElements` uses -- an element counts when its own nodes
// straddle `x` as `[xLo, xHi)`, so a plane on a frame station is served by the bay
// forward of it and not by both. Counting both would double the steel at every seam
// while leaving the *fit* below looking perfectly correct, which is the failure mode
// `scantlings.hpp` records for the same convention.
//
// `displacement` is what `applyBeamLoad` solved for; hand it back the same section
// and the same material, or this is a stress from one model read off another.
std::vector<StressSample> axialStress(const Section& section, const StructuralMaterial& material,
                                      const std::vector<double>& displacement, double x);

// The beam a stress field would have to be, and how much of it is not one.
//
// A volume-weighted least squares of `sigma = axial + gradient * (z - about)`, which
// is the only shape Tier 0 can produce. **`residualRms` is the number with no Tier-0
// counterpart**: it is identically zero for a beam, and on a section it is everything
// the plane-sections assumption threw away. It is reported next to `peak` because a
// residual means nothing except against the field it is a residual of.
struct BeamFit {
    bool ok = false;
    std::size_t samples = 0;
    double volume = 0;         // m^3 the fit was taken over
    double axial = 0;          // Pa at z = `about`
    double gradient = 0;       // Pa/m
    double neutralAxis = 0;    // m where the fitted line crosses zero; `about` if flat
    double residualRms = 0;    // Pa, volume weighted
    double residualWorst = 0;  // Pa
    double peak = 0;           // Pa, largest |sigma| in the field, signed
    Vec3   peakAt{};
    double zLo = 0, zHi = 0;   // m, the extreme fibres the samples reached
};
BeamFit fitBeam(const std::vector<StressSample>& samples, double about);

// The extreme fibre a section modulus is taken at, as the stress actually found
// there. `band` is the fraction of the sampled depth counted as "the deck" or "the
// keel", because a fibre is a line and a mesh has elements.
//
// **`worst / mean` is shear lag stated as a number.** A beam says the two are equal:
// every point at one `z` carries one stress. A wide deck does not -- the stress
// crowds towards the sides, which is where the shell can feed it in -- and that ratio
// is a thing Tier 0 has no way to report at all.
struct FibreStress {
    bool ok = false;
    std::size_t samples = 0;
    double z = 0;       // m, volume-weighted mean height of the band
    double mean = 0;    // Pa
    double worst = 0;   // Pa, largest magnitude, signed
    double volume = 0;  // m^3
};
FibreStress fibreStress(const std::vector<StressSample>& samples, bool deck, double band = 0.02);

// --- 7. Reaching the whole ship: the collapsed element ----------------------------
//
// **Everything above was measured on a quarter of this ship, and the reason was one
// line in a validity check.** Of the 49 two-bay windows of this ferry, 21 meshed,
// reduced and solved; the rest reported "an element came out inverted or degenerate"
// and were refused by `applyBeamLoad` and `reduction::Substructure`.
//
// **And they were not a range.** The 21 were five islands -- x = -43.2..-38.4,
// -36..-31.2, -28.8..-9.6, -7.2..19.2 and 21.6..28.8 -- so 62.4 m of 120 m worked in
// total while the longest *unbroken* run was **26.4 m**, x = -7.2 to 19.2, which is
// exactly the eleven-bay hold every figure in the sections above was measured on. A
// chain needs the unbroken run and got a fifth of the ship. The two windows
// containing the bulkhead at x = -8 failed while sitting in the middle of the best
// island, which is the clue that the cause was never "the ends are curved".
//
// **Nothing was ever inverted.** Over all 49 two-bay windows of this ship there is
// not one negative Jacobian anywhere, and the worst determinant over the sound
// elements is 2.7e-5 -- at the stem and the stern as much as amidships. The
// mechanism is not curvature approaching the plate thickness, not warp, not taper:
//
//   * `makeStructuralMesh` emits **degenerate `PlatePanel`s** -- quads with two
//     coincident corners, which are triangles. 90 bulkhead panels and 76 deck panels
//     of the ferry's 8 900, no shell panels at all: a bulkhead grid running into the
//     centreline or the keel, and a deck laid on fixed |y| lines clipped to a hull
//     that narrows past them. 39.1 m^2 of 12 802.9, **0.305%** of her plating.
//   * Extruding one gives a **collapsed hexahedron** -- a triangular prism written
//     in eight nodes, two of them the same node. One covariant base vector vanishes
//     on the closed edge, so `det J` is *exactly* zero there.
//   * `solidshell::smallestJacobian` samples the eight **corners**. That is the right
//     test for a general hex and the wrong one for this: the element is integrated at
//     the centre and the 2x2x2 Gauss points, and at every one of those it is sound.
//     On the ferry the collapsed elements' worst Gauss determinant is 5.6e-6 against
//     2.7e-5 for the worst *sound* element -- no closer to singular than the mesh
//     already was.
//
// So the fix is a classification rather than a loosening: `solidshell::ElementShape`
// separates *collapsed* from *inverted*, and a section is refused only when a corner
// that nothing coincides with has gone non-positive, or when the quadrature has.
// `tests/test_solid_shell.cpp` carries three negative controls for that, including a
// wedge folded through its thickness at one corner -- which the quadrature alone
// would accept and the corner rule refuses.
//
// **What it cost, and what a wedge is worth.** A collapsed hex is a stiffer
// approximation than the quad it replaces, the way a constant-strain triangle is
// stiffer than a bilinear quad. Measured on the box with **every** panel triangulated
// -- the worst case, against 1.9% of elements on the ferry:
//
//     subdivision   EA/exact      EI/exact    GJ/Bredt tied    (quads, tied)
//     -----------------------------------------------------------------------
//         1        1.00000000    1.1370088       1.2015           1.0986
//         2        1.00000000    1.0068191       1.0516           1.0512
//         3        1.00000000    1.0021759       1.0315           1.0367
//         4        1.00000000    1.0011108       1.0126           1.0297
//
// `EA` is **exact at every refinement**, because a collapsed element still passes the
// patch test and `EA` is the integral of a constant stress; `EI` converges away at
// twenty fold on the first refinement. It is a discretisation error and not the
// formulation, which is why the test asserts the convergence rather than the number.
// On the ferry the wedges are 0.305% of the plating, so even the coarse 14% is
// 0.04% of a section's bending stiffness.
//
// --- What that buys, measured by `tools/section_probe --scan=2` --------------------
//
//                                                        before          after
//     ---------------------------------------------------------------------------
//     windows that mesh, reduce and solve                21 / 49        49 / 49
//     of those, a single connected piece                 21 / 49        46 / 49
//     union of hull inside a working window              62.4 m (52%)   120.0 m (100%)
//     longest unbroken run of working windows            26.4 m         120.0 m
//     the whole 120 m as one section                     refused        8 900 elements,
//                                                                       one component
//     the whole 120 m as a chain of five                 refused        one piece, 5 058
//                                                                       boundary DOF
//
// **One window regressed and netting it off would be dishonest**: x = -43.2..-38.4
// used to come out in one piece and now comes out in two, because
// `junctionWeightLimit` below refuses a tie there that was closing a 22.4 mm gap. It
// still meshes, reduces and solves; it is one of the three that separate the two
// counts above.
//
// --- What it did *not* fix, and one thing it exposed --------------------------------
//
//   * **A tie can still take more of a slave's mass than the slave has.** The
//     junction search accepts a master whose mid-surface is within
//     `junctionTolerance` -- 25 mm, absolute -- while the through-thickness split is
//     relative to the master's own thickness. On 10 mm plating a slave 21.5 mm off
//     the mid-surface is admitted and asks for a weight of 2.65, which puts -1.65 of
//     the slave's steel on one master face and makes `reduction::Substructure` refuse
//     the **whole section** for a non-positive nodal mass. That is now
//     `SectionParams::junctionWeightLimit`, which refuses the junction instead. It
//     costs three of the 49 windows their single-component status -- those junctions
//     were closing a 21 mm gap -- and it is why the two counts above differ. The
//     fix that would need no limit is a `junctionTolerance` scaled off the plating,
//     which this file already says it wants and which changes the junction *census*
//     as well as the tie, so it is a separate piece of work.
//   * **`hullGirderSection` sampled exactly on a frame station lost 76% of the
//     ship, and does not any more.** It read 0.42932 m^2 at x = 19.2, 21.6 and
//     24.0 against 1.80133 at 19.6 and 18.8 -- the plating either side of the plane
//     failing the half-open `straddles` test in `sectionElements` on a
//     floating-point knife edge, leaving only the longitudinals. Not every station
//     (16.8 was fine), which is what said it was representation rather than
//     geometry.
//
//     `scantlings.cpp` carries the fix, with this measurement quoted back at it:
//     `straddles` is tolerant by `kFlat` while the crossing search is exact, so the
//     plane could be admitted while lying a hair outside the panel -- the ferry's
//     station 33 is 19.200000000000003, so a cut asked for at 19.2 sat 3.6e-15 aft
//     of the bay that owned it. The cut is clamped to the panel's own edge now.
//     Re-measured: 1.80133 at 19.2, identical to 18.8 and 16.8 to five decimals,
//     and 1.80129, 1.80113, 1.80083 at 19.6, 21.6 and 24.0 -- the gentle taper the
//     ship actually has.
//
//     This paragraph said "fixing it belongs with `girder.hpp`, and it moves
//     published figures" for as long as the fix had already landed elsewhere.
//     `tools/section_probe --scan` still samples Tier 0 at the centre of a bay,
//     which is now a choice rather than a workaround.
//   * **The thickness-seam rule costs five times more at the ends than amidships,
//     and how much depends on where the section was cut.** `nodeThickness` was an
//     area-weighted mean over the sub-quads *inside* the section, so a station where
//     a strake steps carried one thickness to the section aft of it, another to the
//     one forward of it, and the mean to a section spanning it. A member run stops at
//     a thickness change (§3), so a spanning section stopped runs its two halves never
//     saw and dropped the seam node's run of one. Measured: cutting a window in two
//     conserved the stiffener steel **exactly** amidships, and at x = -24 .. -19.2 --
//     inside the range that always worked, so this is older than the reach -- and at
//     the bow shoulder it did not. **The direction is the spanning section's**: the
//     one piece carries 1.9% less than its two halves at x = -24 .. -19.2 and
//     **25.2%** less at x = 40.8 .. 45.6, because the one piece is the one that sees
//     the seam. It is the nodal-thickness twin of §6 note 1's nodal normal and it has
//     the same fix: a halo, which is §8 and is now built.
//
// --- What mutation testing said about §7 -------------------------------------------
//
// Nineteen single edits across `solid_shell.cpp`, `section.cpp` and `reduction.cpp`.
// The first run killed twelve, and **the seven survivors were worth more than the
// score**, because five of them said the same thing: a predicate is only tested by an
// input whose verdict it changes, and every element the suite fed `elementShape` was
// one it was already happy with.
//
//   * Dropping the quadrature test from `integrable`, dropping the centre sample,
//     dropping a Gauss point, and sampling one zeta level instead of two **all
//     survived**. Nothing had ever handed it an element that is sound at all eight
//     corners and folded inside. One now does -- found by searching a quarter-unit
//     lattice, because no closed form hands such a thing over -- along with a sweep of
//     two hundred elements checked against an independent evaluation of the same nine
//     determinants, and an element deliberately pinched at its centre, which is where
//     no random element's minimum ever landed.
//   * `Section::invertedElements` never firing survived, because every section in the
//     suite has none. The case now tested is the one the check exists for: 1.00 m of
//     plating over a 0.10 m fold, where the extrusion turns through itself. At 0.15 m
//     the same fold leaves every *Gauss point* positive and still inverts two corners,
//     which is the control for the other half.
//   * Counting a wedge's apex twice in the nodal averages survived. A collapsed
//     sub-quad names its apex in two of its four corners, and nothing measured a nodal
//     thickness where that mattered; a two-panel fixture with different thicknesses
//     either side now does, against the closed form.
//
// With those tests the score is eighteen of nineteen. **The one that is left is
// equivalent on every input this ship has, and saying which is more useful than the
// number:** replacing `min(w, 1 - w)` with `min(w, w)` in the weight limit changes
// nothing, because every over-weight tie on this hull has `w` itself negative --
// -1.649 at x = 45.6..50.4, -1.104 at -43.2..-38.4 -- so `w` is already the smaller of
// the pair. Separating them needs a junction whose slave lies past the master's *other*
// face, which this hull does not have.
//
// Two mutants were killed **only by a segmentation fault with no failure line at
// all** -- reverting the corner rule in `elementShape` and in `solveStatic`, both of
// which reach `reduction::Substructure`'s refusal path. A harness that counted `FAIL`
// lines would have scored them as survivors, which is the mistake §5's run records
// making the first time round.
//
// --- 8. The halo: a section is a property of the ship, not of the window -----------
//
// **A node is the mid-surface offset by `t/2` along the nodal normal, and both of
// those were averages over the sub-quads inside the section.** So both were
// statements about the *cut*, and the same physical station came out in a different
// place depending on which window contained it. §6 note 1 records one half of that --
// a nodal normal leaning aft for one section and forward for the other, 1.0 mm at the
// bow and 336 of a plane's 1 170 boundary DOF finding no partner at
// `matchBoundaries`' default. §7's third bullet records the other -- a nodal
// thickness that made a member run break in a spanning section and not in either
// half, 25.2% of the stiffener steel at the bow shoulder. They are one defect.
//
// The fix is to average both over a **halo**: one bay of plating beyond each cut
// plane, taken into the two nodal averages and then thrown away, so that the section
// meshes only what is inside but places its nodes knowing what is outside.
//
// **A bay is not a distance, and it is not a corner test either.** The halo is every
// panel whose extent along x reaches the window, and the **weld** decides what that
// is worth: a halo panel's grid points go through the same `Welder` as the section's
// own, so one that lands on a node of an inside sub-quad joins that node's averages,
// and one that does not contributes only to nodes nobody keeps and is compacted away.
// Nothing has to be predicted and there is no second notion of "the same point" --
// which is the trapdoor CLAUDE.md records `sectionElements` falling through.
//
// The bound is exact, which is why it is allowed to be a bound: a panel reaching a
// node of an inside sub-quad has a grid point at that node's x, which is between the
// planes, so its own extent must reach them.
//
// **This was a corner test first, and mutation testing is what said it should not
// be.** "Shares a welded corner with a panel inside" is the same set only when every
// node is a panel corner; at `subdivision > 1` a node can sit part way along an edge,
// where a panel butting mid-edge reaches it and shares no corner at all. Three
// mutants survived on that predicate -- widening it, an off-by-one in it, a second
// weld tolerance inside it -- all equivalent on a conforming hull and none of them
// equivalent in general. Deleting it is cheaper than assuming it, in code and in
// time: 50 one-bay windows of the ferry mesh in 0.15 s rather than 0.19.
//
// **It is exact rather than close, and that is worth the one line it costs.** The
// halo's panels join `candidate` in ascending panel index, so a node reached by the
// same panels from either side of a cut has its normal and its thickness accumulated
// from the same terms in the same order and comes out at **the same double**. The
// interface test is therefore `worstGap == 0`, not `worstGap < tol`; and
// `tools/section_probe --invariance=2` sweeps 47 stations of the ferry and finds not
// one node moved. A tolerance there would have been a place for the next defect to
// hide.
//
//     ferry, 4.8 m windows either side of a station        before        after
//     ---------------------------------------------------------------------------
//     stations where a node on the plane moves            32 of 47      0 of 47
//     the furthest it moved                               2.5e-3 m      0
//     unmatched boundary DOF at x = -21.6, at 1e-9        336 / 1170    0 / 1170
//     stations whose stiffener steel depends on the cut   32 of 47      0 of 47
//     the worst of that                                   12.2%         3.9e-14
//     stiffener steel, whole ship in 1-bay windows        517 629 kg    501 263 kg
//                                          2-bay          501 427 kg    501 263 kg
//                                          5-bay          501 206 kg    501 263 kg
//                                     one 120 m piece     501 263 kg    501 263 kg
//
// And the solved consequence, on the same plane, `tools/section_probe --chain=2
// --from=-26.4 --to=-16.8`: a chain of two against the same length in one piece, at
// the tolerance each can actually use.
//
//                    before, at 1e-5      after, at the 1e-9 default
//     ---------------------------------------------------------------
//     EA               +3.081e-7                  +3.217e-10
//     z_na             +2.323e-7                  +5.363e-10
//     EI               +3.013e-7                  +4.966e-9
//     GJ               -1.051e-6                  +4.161e-12
//
// The tolerance was not free: merging two DOF 3.4 µm apart is a legitimate answer and
// costs about what the gap is worth. Not needing it is worth three orders in `EA` and
// six in `GJ`, and at the default the same chain used to be refused outright.
//
// **The stiffener block above is the one that says the halo is right rather than
// merely consistent.** Cut the ship four ways and the halo gives one answer; it is the
// answer the whole ship in one piece gives, which is the only one of the four that
// never had a cut to depend on. Without the halo a 1-bay partition carried 3.3%
// *more* steel than the monolith and a 5-bay partition slightly less, because a
// section that cannot see across its own cut does not know there is a strake seam
// there and keeps a run the ship does not have. The steel the shoulder appears to
// lose -- 18 270 kg over two 2.4 m windows becoming 13 669 -- is that steel, and it
// was never there.
//
// What remains is §3's own cost, unchanged and now uniform: 501 263 kg of the
// members' own 631 451 kg, **79.4%**, of which the thickness-seam rule is part and
// `membersMissed` and `membersRefused` are the rest. The fix for that half is still
// the per-station thickness in `constraint::SeamRun` §3 names, and the halo has made
// it worth more rather than less -- it is now the only reason a run breaks at a seam
// at all, and it costs the same wherever the section was cut.
//
// **The blast radius is one deletion.** The halo's sub-quads are dropped, and the
// nodes only they reached are compacted out, immediately after the two averages are
// formed and before anything else in `buildSection` runs. The free-edge count, the
// junction search, the surface census, the component walk and the member runs
// therefore see exactly the sub-quad list they always saw -- which is what the
// warning in §6 was about, answered by being right in one place instead of six.
// `Section::surfaces` is counted after the compaction for the same reason, so a
// surface reached only by halo plating is not one the section is made of.
//
// **What it does not move.** Amidships the halo's sub-quads carry the same normal and
// the same thickness as the ones inside, so the average is the same vector formed
// from twice as many terms and a node moves by 1.7e-18 m -- a few ulps of `t/2`, and
// fourteen orders below the millimetre it removes at the ends. The reach is unchanged
// at 49 of 49 windows and 46 of them in one piece; the ferry hold's first
// fixed-interface mode (0.7785 untied, 2.3026 tied), its DOF half-bandwidth (146 and
// 1 520), its component count (7 and 1) and its 309.6 m of tied edge are unchanged to
// every digit, and its `GJ` moves 3e-5 -- 5.2388e12 to 5.2387e12 tied. A chain of
// four reproduces the same length in one piece exactly as well as it did: EA 9.9e-13,
// EI 1.4e-10, GJ 1.3e-10 on the box, all three unmoved.
//
// **What it does move, and it is not a rounding.** The hold's forward plane at
// x = 19.2 is a strake seam, so a hold-length section can now see one, and four
// members lying wholly in the bay against it are left in runs of two stations of
// different thickness -- no run at all. `A_eff` goes 1.73394 to 1.73266 and `I_eff`
// 43.913 to 43.868, both about a tenth of a percent, and `missedMemberArea` 0.07464
// to 0.07576 m^2. That is the *consistent* answer rather than a worse one -- it is
// what the same plating gives in any other window containing x = 19.2 -- and the
// resolution table in §4 has been re-measured against it. A section's plate mass is
// likewise no longer exactly `sum A_p t_p` over the panels it owns, because §3's
// taper now reaches the cut plane: a window at a strake step hands its neighbour what
// it takes, parts per million, and the sum over any partition of the ship is
// unchanged at 1 084 106 kg.
//
// Meshing the whole ship costs 0.44 s against 0.45 s, and 50 one-bay windows of it
// 0.15 s against 0.11 s: one more panel row through the grid, and nothing else.
//
// `SectionParams::halo` turns it off, and off is the negative control every figure
// above is measured against. It is what this file did before, it still reproduces
// every number in §6 note 1's table, and the tests and the probe both run it
// alongside so that "the two agree" cannot be reported by a hull that was prismatic
// all along.
//
// --- 6. A ship: a chain of sections ----------------------------------------------
//
// One section is a component and a component is not a ship. A ship is the whole
// length cut at N+1 transverse planes, each piece reduced once, and the pieces
// assembled -- `reduction.hpp` §Assembling many, whose union-find over the pairwise
// interface maps is the part that stops at two components without it.
//
// Consecutive sections meet at a cut plane and share every degree of freedom on it,
// so **the assembly is exact at zero modes**: the reduced boundary is the union of
// the cut planes, static condensation is exact there for any load with no interior
// load (`reduction.hpp` property 1), and the scatter-add approximates nothing. A
// chain of N sections and the same length meshed as one section therefore agree on
// `EA`, the neutral axis, `EI` and `GJ` to the conditioning of the solves and not to
// a truncation -- which is a much sharper instrument than a convergence would be,
// and it is what `tests/test_section.cpp` asserts. Modes buy dynamics, and the
// chain's fixed-interface spectrum does converge on the monolithic one from above.
//
// --- What a chain is not, measured on the reference ferry -------------------------
//
// **1. The interface nodes were not coincident by construction, and the parallel
// middle body is what hid it. §8 is the fix and this is what it was fixing.** Two
// sections cut on the same frame station agree on every *mid-surface* point of that
// station -- they come from the same panel corners. They did not agree on the mesh
// nodes, because a node is the mid-surface offset by `t/2` along the **nodal
// normal**, and this file averaged that normal over the sub-quads *inside the
// section*: aft of the plane for one, forward of it for the other. Where the hull is
// prismatic those are the same vector and the nodes land on top of each other
// exactly. Where it is not they did not:
//
//     plane x       -45.6  -33.6  -21.6   -9.6    0.0   12.0   24.0   36.0   48.0
//     worst gap m   6.4e-4 3.5e-4 3.4e-6    0      0      0   1.9e-4 7.1e-4 1.0e-3
//     nodes > 1e-9    116    112    112      0      0      0    112    116    120
//
// That row was measured when most of it was unreachable -- two-bay sections then only
// reduced over 62.4 m of this ship in five islands, and only **x = -21.6** of the nine
// planes above was inside one. §7 removed that limit and every plane in the row can be
// cut now; the figures below are still the ones taken at x = -21.6, because that is
// where they were measured and quoting them anywhere else would be a different claim.
// The gap there is 3.4e-6 m -- three orders
// of magnitude below the plating and three above `matchBoundaries`' default. At that
// default **336 of the plane's 1 170 boundary DOF found no partner**, and a chain
// assembled out of them, solved, and was torn along 29% of the cut. So `Chain` counts
// what did not match rather than trusting the tolerance.
//
// **What the tolerance costs was measured rather than argued.** At 1e-5 the whole
// plane matches and the chain reproduces the same length in one piece to 3.1e-7 in
// `EA`, 3.0e-7 in `EI` and 1.1e-6 in `GJ` -- so joining two nodes 3.4 µm apart is a
// legitimate answer and the reason to set the tolerance was that the default silently
// did something worse.
//
// **The whole row is now zero and the tolerance is not needed at all** --
// `SectionParams::halo`, §8. The table above is kept because it is the negative
// control the fix is measured against and `halo = false` still reproduces every
// figure in it.
//
// One warning that cost a day: **which boundary DOF are on a cut plane is a question
// for the mesher and not for a position test.** `Section::aftNodes` and
// `forwardNodes` are chosen on the mid-surface and take both nodes of each pair,
// because a node whose plating leans out of the transverse plane is up to `t/2` off
// it -- 8 mm here. Selecting the end planes by `x == station` instead drops exactly
// the plating that is not square to the cut, and a chain built that way came out 33%
// short in `EA` on a shoulder while reporting nothing wrong at all. It is the same
// mistake in reverse that hid the 336 above behind a count of 24.
//
// **2. Every interior cut plane used to open a ring of junctions the monolith closes,
// and §9 is what closed it.** A junction node on a cut plane gets no *face* tie --
// §5's master face spans x, so half its nodes are interior to the section, and a
// boundary degree of freedom written as a function of an interior one is not one a
// reduction keeps exactly. Cutting a length into N pieces turns N-1 interior stations
// into interfaces, so a chain used to tie N-1 stations fewer than the monolith does.
// Measured on x = -7.2 to 2.4 m of the ferry, four bays: in one piece 36 ties and
// 72.0 m of junction edge joined; as two two-bay sections, 12 + 12 ties and
// 4.8 + 4.8 m. On the box in the tests it cost **19.2% of `GJ` at N = 8** while `EA`
// stayed exact at 1e-12 at every N, which is why the measurement had to be torsion.
//
// The fix is a tie whose masters all lie *in* the cut plane, which both sections
// derive identically from shared boundary DOF -- a line tie through the plating
// thickness rather than the bilinear face tie §5 builds. It is §9, it is new
// machinery rather than a parameter, and the table above is now the negative control
// `SectionParams::interfaceTies = false` reproduces.
//
// --- What it costs ----------------------------------------------------------------
//
// The interface is the whole cost and it is set by the cross-section, not by the
// length: a transverse cut of the ferry is 1 170 boundary DOF at `subdivision = 1`,
// so a chain of N sections is a **dense** `(N+1) x 1 170` reduced model -- 3 510 for
// two, 5 850 for four -- and the assembled solve is a dense Cholesky at that size.
// Cutting a ship into more pieces makes every piece cheaper to reduce and the
// assembly quadratically dearer to solve, which is the trade a caller is making
// whether or not it is stated. `tools/section_probe --chain=N` measures both halves.
//
// --- 9. The in-plane line tie: what an interior cut plane no longer costs ----------
//
// **The first question was whether this is possible at all, because a prescribed
// degree of freedom cannot also be derived and that is why the ties were dropped.
// The answer is that the premise was half right: an interior cut plane is not
// prescribed, it is *shared*.** The two are different things and §6 note 2 ran them
// together. What a reduction keeps exactly is the interface; what a *load* prescribes
// is the chain's two outermost planes. An interior plane is neither -- it is a set of
// unknowns that two reduced models both write into -- and a relation among those
// unknowns is something the assembled model can carry.
//
// So the impossibility is real but narrower than it looked, and it lands in exactly
// one place: **a `Substructure` cannot carry this constraint and an `Assembly` can.**
//
//   * §5's master is a *face*, and a sub-quad face spans x. Half its nodes are one
//     station inside the section, so a slave on the cut plane would be a function of
//     an interior unknown -- which is not something a reduction keeps exactly, and
//     `reduction::Substructure` refuses it in as many words rather than dropping it.
//   * The master that works is the **line** that surface draws on the same cut plane:
//     the sub-quad edge both of whose ends are on the plane. Two mid-surface nodes,
//     four mesh nodes with the through-thickness split, and every one of them on the
//     plane. `Section::planeTies` is that, and `SectionParams::interfaceTies` is the
//     negative control.
//   * A section therefore **carries** the constraint and does not apply it. Its own
//     two planes are what a load is prescribed on, so applying one there really would
//     be a second claim on one unknown. `buildChain` applies them at the **interior**
//     planes only, on `Assembly::stiffness` and `Assembly::mass`, as `T^T K T`.
//
// **Doing it after the condensation rather than before is exact, not an
// approximation, and that is the load-bearing step.** Guyan condensation is exact for
// any prescribed boundary displacement (`reduction.hpp` property 1), so the reduced
// energy is `0.5 u_b^T K_r u_b` for *every* `u_b`. Restricting `u_b` to the subspace
// the tie defines and relaxing the interior is therefore the same number as applying
// the tie first and condensing after -- the constraint touches boundary rows only, so
// it commutes with eliminating the interior. Nothing is truncated and no penalty is
// introduced, exactly as in §5.
//
// --- What it does, on the box, where every quantity has a closed form ---------------
//
//     sections                        1          2          4          8
//     -------------------------------------------------------------------------
//     junction edge joined, of 64.0 m
//       cut planes open              58.0       52.0       40.0       16.0
//       line ties                    58.0       58.0       58.0       58.0
//     GJ against the same length in one piece
//       cut planes open           -2.6e-12    -3.289%    -9.236%   -19.204%
//       line ties                 -2.6e-12   +9.3e-12   +2.9e-11   -1.6e-11
//     EA against the closed form
//       cut planes open            1.2e-12    8.0e-14    1.5e-12    1.1e-12
//       line ties                  1.2e-12    1.2e-13    1.5e-12    1.4e-12
//
// **58.0 m at every N is the whole claim, and it is the monolith's own figure.** The
// 6.0 m that never joins is the two outermost planes, which one piece leaves open as
// well. `GJ` comes back to the conditioning of two independent solves -- 1e-11, four
// orders below the 3.3% the first cut used to cost -- and it does so at N = 8, where
// every station but two is an interface.
//
// **`EA` is the row that says why none of this could have been validated on `EA`.**
// It is exact to 1e-12 in *both* columns at every N, so an implementation that did
// nothing at all would have scored best on it. §2 says this about a section and it is
// no less true of a chain: prescribing plane sections at the ends makes every
// longitudinal strip carry `sigma = E eps` whatever it is joined to.
//
// --- And on the ferry, where it is a factor of fifty rather than an identity --------
//
//     ferry, tied                             junction edge joined   GJ vs one piece
//     ------------------------------------------------------------------------------
//     x = -7.2..2.4, 2 sections, planes open   9.6 m of 134.4           -3.334%
//     x = -7.2..2.4, 2 sections, line ties     72.0 m (= one piece)     -0.095%
//     x = -7.2..16.8, 5 sections, planes open  28.8 m of 336.0          -8.939%
//     x = -7.2..16.8, 5 sections, line ties    273.6 m (= one piece)    -0.174%
//
// **It joins exactly what the monolith joins and it does not carry exactly what the
// monolith carries, and the reason is worth stating rather than calling it
// conditioning.** On the box the two are the *same constraint*: a slave whose own
// normal is square to the cut lands on the master face's edge, where the bilinear shape
// functions collapse onto the two nodes the line tie uses, so the line tie is the face
// tie restricted rather than an approximation to it. On a real hull the slave lands a
// little off that edge, so the monolith's face tie and the chain's line tie are two
// slightly different constraints on the same junction -- and the residual is what the
// line drops. 0.1% and 0.17% against 3.3% and 8.9%.
//
// The lowest fixed-interface frequency is the third instrument and it separates the
// two by their **sign**. A tied chain is the monolith with each piece reduced, so its
// Rayleigh quotient is an upper bound and it comes down onto the monolith's own
// frequency from above -- 12.44833 Hz against 12.44824 at four sections and six modes.
// A chain whose cut planes are open is a *different, softer* structure and falls
// *through* it, to 12.37090; no number of modes brings it back, because there is
// nothing there to converge to.
//
// --- What the line approximates, which is one number ------------------------------
//
// Projecting onto a face leaves a residual along the master's own normal and nothing
// else, and the through-thickness weight is exactly what that residual becomes -- so
// a face tie drops nothing. Projecting onto a line leaves that, **plus** whatever ran
// along the master's surface across the line. On a cut plane that direction is x, so
// what is dropped is the slave's own nodal normal leaning fore or aft, bounded above
// by `t/2` and measured on the reference ferry at **1.9e-5 m** -- four hundred times
// below the bound, because a transverse cut of a ship is very nearly square to the
// plating it passes through. `Section::worstPlaneTieSlip` is that number and it is
// reported rather than assumed; a mesher run on plating raked hard against its own
// cut planes would see it grow, and would say so.
//
// Everything else is §5's, because it is the same construction: the weights are a
// partition of unity so a rigid translation is exact, the masters interpolate the
// slave's own point so a rigid rotation is too, the overshoot outside the segment is
// bounded by `SectionParams::junctionOvershoot`, the through-thickness weight by
// `junctionWeightLimit`, and a master that is itself a slave is refused rather than
// composed. On the ferry the line ties run at the same `w = 1.69` the face ties do,
// which is the same 9 mm gap seen from the same side.
//
// --- The two sides of a cut derive the same tie, and that is checked ---------------
//
// A constraint on an interior plane is one relation among shared unknowns, so the two
// sections either side of it must agree about it or the ship is tied to itself twice.
// They do, **to the last bit**, and it is §8 that makes them: the plane's mid-surface
// nodes, nodal normals and nodal thicknesses are properties of the ship rather than of
// the window, so both sections start from the same doubles. The one thing that is not
// geometry -- which of several equidistant segments wins -- is broken by the segment's
// own endpoint positions rather than by a node index, because node indices are a
// property of the window and positions are not.
//
// `buildChain` derives every interior plane twice and compares before adopting one.
// `Chain::planeTiesDisagreeing` is zero on the box at N = 2, 4 and 8 and on the ferry,
// and `worstPlaneTieDisagreement` is **0.0 rather than small**. A chain built on a
// plane the two sides disagreed about applies nothing there and says which plane.
//
// **Exact equality and not a tolerance, and the negative control is what says so.**
// Turn the halo off and cut the ferry's stern shoulder at x = -21.6: the two sides pick
// the *same* masters -- those come from panel corners and never depended on the cut --
// and differ only in the weights, by **9.63e-8**. A structural comparison finds nothing
// and a 1e-6 tolerance finds nothing; comparing against zero finds it. It is the same
// argument §8 makes for testing `worstGap == 0`, and the same place the next defect
// would have hidden.
//
// --- What it costs ----------------------------------------------------------------
//
// Nothing measurable. A tie is `6` constraints per node over `4` masters and the fold
// is `O(n)` per constraint against a dense assembly that is `O(n^2)` to build and
// `O(n^3)` to factor: the box at N = 8 has 28 tied nodes and the fold is under a
// millisecond of a solve that is seconds. There is no band to widen -- the assembled
// model is dense already -- which is the one way this is *cheaper* than §5, whose ties
// cost the ferry hold a half-bandwidth of 146 against 1 520.
//
// The one thing a caller has to know is that **`Chain::planeTieDof` must be held**. An
// eliminated row is left isolated with a unit diagonal, so a solve that leaves it free
// returns zero there, and `assembledFrequencies` would return one spurious eigenvalue
// of 1 rad/s per tie -- 0.159 Hz, underneath every elastic mode a ship has. All three
// of `applyBeamLoad`, `applyTwist` and `chainFrequencies` hold them and fill them back
// in from their masters afterwards.
//
// --- What mutation testing said, which is mostly what it demanded be built ---------
//
// Thirty-seven single edits across this file's implementation; twenty-nine are killed.
// **The first pass killed eighteen, and the sixteen it left standing were worth more
// than the score, because eleven of them said the same thing: the box's corners are a
// butt joint, and a butt joint exercises almost none of a tie.** Two mid-surfaces
// meeting on a corner line have no offset along either's normal, so the
// through-thickness weight is 1/2 whatever the plating is, both halves of the pair
// weigh the same, and swapping them is a no-op. Every junction figure in §5 and §9 was
// measured on exactly that. Four fixtures now exist because of it:
//
//   * **A box with a deck laid inboard of its side plating**, stopping 12 mm short of
//     its mid-surface at a height that lands in the *middle* of a side element, over
//     plating that steps from 10 to 20 mm at mid-height. The junction is then at
//     `w = 1.4639` -- a closed form in the gap and the two thicknesses -- so the split
//     is not a half, the two halves of the pair are worth different things, and the
//     master's thickness has to be interpolated along the segment rather than taken
//     from either end. It kills four mutants that no box corner can.
//   * **The same deck at exactly mid-height**, where it straddles the node between two
//     side elements and each half of its pair falls the other way. That is the input
//     the one-segment rule exists for and nothing else here produces it.
//   * **The box wound the other way at every other bay**, which is what
//     `makeStructuralMesh` does when it mirrors the starboard side. The sub-quads
//     either side of a cut then name the segment they share in opposite directions,
//     which is what ordering it by position is for.
//   * **The shoulder chain with the halo off**, where the two sides of the cut pick the
//     same masters and differ only in the weights, by 9.63e-8.
//
// **The eight that are left are equivalent on every input this repository has, and
// saying which is more useful than the score:**
//
//   * **Accepting a master segment on *either* cut plane** changes nothing, because a
//     candidate must already be within `junctionTolerance` -- 25 mm -- of the slave,
//     and the two planes of a section that reduces are metres apart.
//   * **Dropping the geometric tie-break between equidistant segments** changes
//     nothing, because nothing here produces two candidates equal in both overshoot
//     and distance; a slave that straddles a node produces two that differ and is
//     refused by the rule above instead.
//   * **Deleting the line-tie chain refusal** changes nothing, for the same reason §5
//     records for its own: the "already a master" skip absorbs the second side of every
//     junction here. `solidshell::DofExpansion` is the backstop and `buildChain` runs
//     it over the assembled set.
//   * **Leaving the eliminated rows free in `applyTwist`** changes nothing, because the
//     fold already isolates them: holding a row and isolating it are two statements of
//     the same thing and either alone is enough. Deleting the isolation is killed --
//     the matrix is published as `T^T K T` and a caller reads it -- and deleting the
//     hold from `chainFrequencies`, where the isolation is *not* enough, is killed too.
//   * **Three separate edits to the restraint guard** -- letting the anchor, or its
//     partner, or the lever land on an eliminated row -- each change nothing, while
//     **disabling the guard as a whole is killed**. A slave's three axes are eliminated
//     together, so any one of the three tests catches the same node and the other two
//     are redundant *given it*. The guard is what stopped an eight-section box chain
//     going singular in `applyBeamLoad`, which is how it came to be written.
//   * **Not filling the eliminated rows back in** changes nothing, because the only
//     thing that reads them is `BeamResponse::peakDisplacement` and every load here has
//     its peak on a prescribed end-plane row. It is kept because a zero left at an
//     eliminated degree of freedom is precisely the hole `reduction::recover` shipped
//     with -- a displacement field that read 850 788 MPa against a true 427 -- and the
//     next reader of that field would find it.
//   * **Leaving the line ties out of `Chain::components`** changes nothing, and the
//     reason is close to a theorem: a junction line runs along x, so a surface tied at
//     a cut plane is tied again at a station inside the section unless the section is
//     one element long, which `reduction::Substructure` refuses. Separating them needs
//     a junction lying wholly *in* a cut plane -- a bulkhead on a frame station -- which
//     this ship does not have, because none of its bulkheads is on one.
//
// One mutant was killed **only by a segmentation fault with no failure line**: tying
// the outermost planes as well, which puts a constraint on a prescribed row and takes
// `reduction::assembledStaticSolve` through a partition it cannot make. A harness that
// counted `FAIL` lines would have scored it a survivor, which is the mistake §5's run
// records making and §7's records making again.

struct ChainParams {
    // The cut planes, ascending: N+1 of them make N sections. **Frame stations**, for
    // the reason `SectionParams::xFrom` gives -- a plane anywhere else passes through
    // panels, which are then dropped from both neighbours and leave a gap in the
    // chain rather than a short bay.
    std::vector<double> station;

    // What each section is meshed with. `xFrom` and `xTo` are overwritten per
    // section; everything else is used as given.
    SectionParams section;

    // How each section is reduced. `modes = 0` is Guyan, which is exact at the
    // interface for any static load and is what the static comparisons here use.
    reduction::ReduceParams reduce;

    // How far apart two boundary DOF on a shared cut plane may be and still be the
    // same unknown. The default is `matchBoundaries`' own, and with
    // `SectionParams::halo` on it is not a tolerance at all: two sections cut on the
    // same station put that station's nodes at the same double, so the whole plane
    // matches at any positive value. It was 3.4 µm out at x = -21.6 and 1.0 mm at the
    // bow before -- see note 1 and §8 -- which is what `Chain::unmatched` exists to
    // say rather than leave to be assumed, and it still says it under `halo = false`.
    double matchTolerance = 1e-9;
};

// A ship built from sections: the pieces, their reductions, and the one assembled
// model.
//
// `substructure` holds a reference to `section[i].mesh`, so neither vector may be
// reallocated once the chain is built. Both are sized once here and never grown,
// which is why this is a struct a caller receives rather than one it fills in.
struct Chain {
    std::vector<Section> section;                      // in order, aft to forward
    std::vector<reduction::Substructure> substructure;  // one per section
    std::vector<reduction::Reduction> reduced;
    reduction::Assembly assembly;

    // The assembled boundary DOF on the two outermost cut planes, which is what a
    // caller prescribes to load the chain like a beam. Everything else -- every
    // interior cut plane, every modal coordinate -- is free.
    std::vector<std::uint32_t> aftDof, forwardDof;

    // Per interior cut plane, N-1 of them: how many boundary DOF the two sections
    // either side of it turned out to share, and how many of the aft section's
    // forward-plane DOF found no partner. **`unmatched` is the number that matters**;
    // a chain assembles and solves with it non-zero, and what it has then built is a
    // ship torn along part of a cut.
    std::vector<std::size_t> shared, unmatched;
    double worstGap = 0;  // m, the furthest apart a matched pair actually was

    // **Two component counts, because there are two questions and neither answers
    // the other.** `reduction::assembledComponents` works on the reduced model, where
    // a component's pair is dense over its own DOF whether the mesh behind it is one
    // piece or seven -- so it sees a chain that failed to join at a cut plane and is
    // blind to a section whose decks are not tied to its shell. `Section::components`
    // is the other way round. `components` here joins the two: a mesh component of
    // one section and a mesh component of the next are the same piece of ship when
    // they share an assembled boundary row. **That is the one that must be 1.**
    int components = 0;
    int reducedComponents = 0;
    // Which piece each `(section, mesh component)` belongs to, flattened by
    // `componentBase[i] + c`. Exposed for the same reason `Section::componentOf` is:
    // every rigid body motion is per piece, so anything restraining a chain has to
    // restrain each of them -- `applyBeamLoad` below does, and an untied chain of a
    // real ship has seven.
    std::vector<int> pieceOf;
    std::vector<std::size_t> componentBase;
    // Junction ties the cut planes cost, against what the same length in one piece
    // would have tied: see note 2. `tiedEdges` is summed over the sections **and over
    // the interior planes the in-plane line ties of §9 close**, which is what makes it
    // reach the monolith's own figure instead of falling by one station per cut.
    double tiedEdges = 0, junctionEdges = 0;

    // --- The in-plane line ties, §9 --------------------------------------------------
    //
    // The ring of junctions on each interior cut plane, applied to the assembled model:
    // `assembly.stiffness` and `assembly.mass` are `T^T K T` and `T^T M T` by the time
    // a caller sees them, and the eliminated rows are left isolated.
    //
    // **`planeTieDof` has to be held.** An eliminated row carries a unit diagonal and
    // nothing else, so a solve that leaves it free returns a zero displacement there
    // and -- for `assembledFrequencies` -- one spurious eigenvalue of 1 rad/s per row,
    // which would land underneath every elastic mode of a ship. All three of
    // `applyBeamLoad`, `applyTwist` and `chainFrequencies` below hold them and then
    // fill them back in from their masters.
    std::vector<solidshell::Mpc> planeTies;    // in assembled degrees of freedom
    std::vector<std::uint32_t> planeTieDof;    // the rows they eliminate, ascending
    int planeTieNodes = 0;
    // Interior planes whose two sections did **not** produce the same constraint. Zero
    // is the claim §9 rests on -- with the halo a station's nodes, normals and
    // thicknesses are the ship's rather than the window's, so the two sections derive
    // the same tie from the same doubles -- and it is checked rather than assumed,
    // because a chain built on a plane the two sides disagree about is a ship tied to
    // itself twice in two different ways.
    int planeTiesDisagreeing = 0;
    double worstPlaneTieDisagreement = 0;

    double meshSeconds = 0, reduceSeconds = 0, assembleSeconds = 0;
    std::size_t boundaryDof() const { return static_cast<std::size_t>(assembly.boundary); }
    std::vector<std::string> problems;

    double xFrom = 0, xTo = 0;
    double length() const { return xTo - xFrom; }

    // Every section meshed, every reduction ready, and nothing unmatched.
    bool ready() const;
};

Chain buildChain(const StructuralMesh& structure, const ChainParams& params);

// `applyBeamLoad` and `applyTwist`, on the assembled model instead of one section.
//
// The same field on the same two planes: `u_x = strain * x + curvature * x * (z -
// reference)` on the chain's outermost cuts with `u_y` and `u_z` free, and for the
// twist a rigid disc with every degree of freedom prescribed. What the reduced model
// cannot do is restrain a rigid body mode on an *interior* node, because it has none
// -- so the three restraints that take out the translations in y and z and the
// rotation about x go on boundary DOF of the outermost planes, chosen the same way
// `applyBeamLoad` chooses them and reported through `restraintReaction` the same way.
//
// A chain in more than one piece is refused rather than solved: the second piece is a
// mechanism and the factorisation would report a rigid body mode without saying which
// of the two it belonged to.
//
// `stateOut` receives the solved assembled state -- boundary rows and modal
// coordinates together, with the interior planes' eliminated rows filled in from
// their masters. It is what `sectionDisplacement` expands back onto a section's own
// mesh, and it is the only route from a whole-ship model to a stress.
BeamResponse applyBeamLoad(const Chain& chain, const BeamLoad& load,
                           std::vector<double>* stateOut = nullptr);
TorsionResponse applyTwist(const Chain& chain, double twist, double reference = 0.0);

// One section's own mesh displacement field, out of a solved assembled state.
//
// `reduction::componentState` picks that section's reduced coordinates out of the
// assembly and `reduction::recover` expands them through `u = T x`, so what comes
// back is the *exact* Guyan interior for the boundary the chain settled on -- not an
// interpolation of it, and not a second solve. Empty if the state is not this
// chain's, for the same reason `componentState` refuses one: every index would
// otherwise be in range and a short state would come back as a plausible field
// quietly missing its modal content.
//
// **This is what makes a whole-ship Tier-1 model able to answer a stress question at
// all.** Everything above it reports resultants, and a resultant is what Tier 0
// already had.
std::vector<double> sectionDisplacement(const Chain& chain, std::size_t index,
                                        const std::vector<double>& state);

// The chain's lowest natural frequencies with both end planes held, rad/s ascending.
//
// This is the same boundary condition `reduction::Substructure::fixedInterfaceModes`
// applies to the same length in one piece, so the two are directly comparable and the
// monolithic answer owes nothing to any of this. It is an upper bound and it falls as
// modes are added, which is the property `reduction.hpp` §1 item 3 is about.
//
// **It is the instrument that sees whether the sections are joined.** A chain whose
// interfaces did not match leaves its interior sections held at nothing at all, so
// its spectrum opens with six near-zero rigid body modes per loose section rather
// than with an elastic mode. A chain of sections whose *plating* is not tied reads
// the decks' own frequency instead, which is the failure §2 of this file is about.
// Neither is visible in `EA` or `EI`.
std::vector<double> chainFrequencies(const Chain& chain);

// --- 11. The whole ship, and what she does not agree with Tier 0 about --------------
//
// Every piece above existed and nobody had put them together. The mesher reaches all
// 120 m (§7), the halo makes a station's nodes a property of the ship rather than of
// the cut (§8), and the in-plane line tie stops an interior plane costing torsion
// (§9) -- and the model those three exist to make had never been built, nor asked
// whether it agrees with the beam it is supposed to refine. `tools/section_probe
// --whole`, `--profile` and `--wave` are that, and the disagreements are the part
// worth reading.
//
// **All timings below were taken with the machine also running a mutation harness.**
// They are the shape of the cost, not a benchmark.
//
// --- What the whole ship costs, both ways ---------------------------------------
//
//     the reference ferry, subdivision 1, -60 .. 60 m
//     ------------------------------------------------------------------------------
//     one piece    mesh 0.44 s: 8 900 elements, 18 780 nodes, **one** component
//                  1 084 106 kg of plate + 501 263 kg of member -- the same two
//                  figures §8's partition table publishes, to the kilogram
//                  1 408.6 m of junction edge, 1 337.6 m of it tied
//                  DOF half-bandwidth **5 384**; three banded solves **1 044 s**
//     chain of 5   build 221 s (mesh 0.2, reduce 221, assemble 0.1)
//                  5 058 assembled DOF, 390 MB dense, 0 unmatched, worst gap **0.0 m**
//                  three dense solves **9.0 s**
//     chain of 10  build 93 s (mesh 0.1, reduce 92, assemble 0.4)
//                  9 828 assembled DOF, 1 474 MB dense, dense solve **69 s**
//
// **The band is why one piece is not a model.** `solveStatic` factors at `n b^2`: the
// eleven-bay hold is 1 520 tied and the whole ship is 5 384 over 56 340 unknowns --
// 1.6e12 flops against 5.5e10, thirty times -- and it measures 1 044 s against 5.3. A
// tie joins nodes no element edge joins, so the ordering has to carry the whole
// cross-section *and* every junction across fifty bays. It is not a property of this
// ordering. The monolithic solve is therefore behind `--whole-solve` and is in no gate.
//
// **And more pieces is not simply better.** Cutting finer makes every reduction
// cheaper -- 221 s at five, 93 s at ten, because a Guyan condensation is
// `O(n_i b_i^2)` and both shrink -- and makes the assembled model quadratically larger
// to hold and cubically dearer to factor: 9 s at five, 69 s at ten, 390 MB against
// 1.5 GB. The minimum is somewhere in between and it is set by the cross-section, not
// by the length.
//
// --- The chain against the monolith at ship scale, which is not 1e-10 ---------------
//
//     whole ship        one piece      chain of 5     chain of 10    5 against 10
//     -------------------------------------------------------------------------------
//     EA  N            1.47220e11      1.47255e11     1.47342e11      +5.9e-4
//     z_na m              4.97345         4.97283        4.97306      +4.6e-5
//     EI  N m^2        2.38603e12      2.38679e12     2.38676e12      -1.5e-5
//     GJ  N m^2        1.25970e12      1.26030e12     1.26000e12      -2.4e-4
//
// **§6 says a chain and the same length in one piece agree "to the conditioning of the
// solves and not to a truncation". That is still true and the number is 3e-4, not the
// 1e-10 the box reaches, and the reason is §9's rather than the solver's.** Nothing is
// truncated: `unmatched` is zero, `worstGap` is exactly 0.0 m, `planeTiesDisagreeing`
// is zero, and both chains tie 1 340.8 m of the monolith's 1 337.6 m. What differs is
// that on the box a cut plane's *line* tie and the monolith's *face* tie are the same
// constraint -- the slave lands on the face's own edge, where the bilinear shape
// functions collapse onto the two nodes the line uses -- and on a real hull it lands a
// little off it. §9 measures that trade at 0.095% of `GJ` for a two-section chain of
// the hold; five more cut planes in a 120 m ship come to 0.024%, in the same direction:
// the finer chain is the softer one. **Quoting the box's 1e-10 as a ship-scale figure
// would be this file's own recorded failure mode, so the ship-scale figure is measured
// and sits next to it.**
//
// The gate can afford **two partitions against each other** and not the 1 044 s
// monolith. It is the same claim: an interior cut plane is *shared* rather than
// prescribed, so cutting one ship two ways has to give one model.
//
// --- Against Tier 0 along the length, where two of three differences are accounting --
//
// `--profile=2` tiles the hull with two-bay windows and asks both tiers for `A`, the
// neutral axis and `I`. Amidships they agree to **+0.356% in area and +0.347% in second
// moment**; at the ends Tier 1 reads **80% low**. Three separate things are in that and
// only the third is physics.
//
// **1. An area alone cannot correct an `I` comparison.** `missedMemberArea` was the
// whole of the accounting, and the girders this mesh cannot reach -- three of them, off
// the longitudinal spacing, so no node of a mesh built from panel seams lands on one --
// are 4.4% of the ferry's area and **5.3% of her second moment**, because they sit low
// in a double bottom and a second moment is a lever arm squared. Subtracting the area
// and not the moment left the two tiers looking 5.0% apart amidships where they agree
// to 0.35%. `missedMemberFirstMoment` and `missedMemberSecondMoment` are the fix, and
// `tests/test_section.cpp` asserts all three against `sectionElements`.
//
// **2. Where plane sections is asserted matters as much as what is compared.** Tier 0
// gives a section at a *station* and asserts plane sections at every one of them; a
// Tier-1 window asserts it at two and lets the plating do as it likes in between. The
// experiment holds the steel fixed and moves only the planes -- one eight-bay window
// cut into k pieces combined in series, as a fraction of Tier 0's own answer for the
// same length:
//
//     planes apart      2.4 m    4.8 m    9.6 m   19.2 m
//     ---------------------------------------------------
//     stern            0.9847   0.5161   0.5421   0.6275
//     amidships        1.0059   1.0055   1.0075   1.0078
//     bow              0.9303   0.3153   0.3324   0.3951
//
// Amidships nothing moves. At the ends the answer halves as soon as the planes are more
// than one bay apart, and **the first hypothesis -- section-level averaging -- is
// measured and rejected**: the harmonic mean of a window's own station properties is
// within one per cent of the arithmetic mean everywhere, because the section *total*
// barely changes over two bays even at the stem. What changes is which steel is
// continuous from one plane to the other: a girth band closing to nothing, a strake
// turning to meet the stem and spending its length athwartships, where it carries no
// longitudinal stress and Tier 0 still counts the full cut through it.
//
// It is not monotone in between, and that was predicted wrongly here before it was
// measured: eight bays cut into four pieces and into two are cut at *different frames*,
// so which station a plane lands on matters as well as how far apart they are. The
// claim is carried by the two ends of the sweep.
//
// **At the ends the FEM is the better answer and the beam is not imprecise, it is
// answering a different question.** No station-by-station section property can say that
// the material at one station is not the material at the next.
//
// **3. What is left amidships is the frames, and there the FEM is right.**
// `hullGirderSection` drops every member with no extent along x -- an athwartships
// member carries no longitudinal stress -- so Tier 0 scores a structure with no frames
// in it *identically*, which is checked rather than assumed. Tier 1 does not: a strip
// that cannot contract in y and z carries more than `E eps` for the same strain.
// Measured on the two-bay midship window, **+0.3561% with the frames and -0.0828%
// without, so they are worth +0.4389%** -- and without them the two tiers agree to a
// tenth of a per cent, which is what says the rest of the accounting is right. The
// effect grows with section length (0.31% at one bay, 0.44% at two, 0.48% at four)
// because a cut plane is free in y and z and a longer section has proportionally less
// of itself next to one.
//
// **Two figures `tests/test_section.cpp` carried for that were the wrong section's.**
// "0.44% measured" and "0.52% ... measured by omitting them" sit above a *two-bay*
// test and are the **hold**'s numbers -- and the hold's numbers from before §8's halo.
// Re-measured today: the hold gives +0.411% and +0.508% and the two-bay window gives
// +0.356% and +0.439%. So one of the two was near enough to be believed and neither
// belonged to the section the test builds. Both are corrected in place and both now
// come out of a run rather than out of a comment.
//
// `docs/02-simulation.md` §*Against Tier 0* carries **+0.41%**, **+0.28%** and
// **+0.52%** for the hold and all three reproduce: +0.411%, +0.279% and +0.508%.
//
// --- The hull-girder response: the first time the two tiers were asked one question --
//
// `--wave` poises the ferry on a crest of her own length and applies **the load** to
// the Tier-1 model rather than the answer: `girder.hpp`'s weight minus buoyancy per
// station, spread over the elements in that station's slab in proportion to their
// volume, so the resultant per station -- which is the whole of what sets `V(x)` and
// `M(x)` -- is exact and the local distribution deliberately is not. Handing the model
// Tier 0's own bending moment would assume most of what is being compared.
//
// Three things had to be right and each is a place a load can vanish:
//
//   * **A junction tie's slave has no row.** `reduction::reduceLoad` reads the boundary
//     and interior partitions, and an eliminated degree of freedom is in neither, so a
//     load left on one is **silently dropped** -- 0.84 MN of it on this ship. It belongs
//     to the masters by the transpose of the constraint, exactly as the reaction does in
//     `stiffnessTimes`. The same fold is needed again for the assembly's own line ties.
//     `reduceLoad` is not changed here; the caller folds first, and the hole is worth
//     knowing about before the next caller finds it.
//   * **Six restraints, not three.** `applyBeamLoad` prescribes both end planes and has
//     three motions left; a ship floating free has six. They are statically determinate,
//     so on a balanced load they carry nothing -- measured at 3.8e3 N against a largest
//     applied 1.2e7 -- and that reading is the end-to-end check that the load balanced.
//   * **Guyan is exact at the interface for a load applied inside it**
//     (`reduction.hpp` property 1), so every cut plane's displacement is the monolith's.
//     What zero modes does not buy is the *interior* recovery under an interior load, so
//     the stress inside a section is not read off the reduced model: the chain's own
//     exact interface displacement is prescribed on that section's mesh and it is solved
//     directly with its own share of the load. That is exact and costs one banded
//     factorisation per section.
//
// --- The stress, which is the half a beam cannot produce ----------------------------
//
// At the peak-moment station, driven by the wave load through the whole-ship model:
//
//     x = 6.0 m, M = 4.573e8 N m       Tier 0          Tier 1
//     ---------------------------------------------------------------
//     deck fibre                       82.06 MPa      118.42 MPa
//       the deck's own mean            (one number)    83.26 MPa
//     keel fibre                      -66.52 MPa      -88.49 MPa
//     neutral axis                      6.7132 m        6.7961 m
//     what a beam cannot carry          0               8.06 MPa rms
//
// **The mean over the deck is the beam's answer to 1.5%, and the worst is 42% above
// it.** Both halves matter and neither alone would do: the mean agreeing is what says
// the moment is arriving, and the worst not agreeing is what says the field is not a
// beam. The peak sits at `(5.31, 10.00, 14.81)` -- the deck edge at the ship's side,
// which is where the shell can feed stress in and is exactly where shear lag puts it.
//
// `Section::fibreStress` reports `worst / mean` and it is **1.42** here, rising to 1.50
// at x = 18 -- taken only over the sections carrying at least half the peak moment,
// because towards the ends the deck's own mean passes through zero and the ratio runs
// off to seventy while saying nothing. On the prismatic box in the tests the same ratio
// is 1.006, which is the sampling band's own depth and is the floor this instrument can
// resolve.
//
// **It needed the real load to see.** Shear lag is driven by `dM/dx`, and a section
// handed a constant moment -- which is all `applyBeamLoad` can prescribe -- has none of
// it to show. That is why the wave *load* is applied and not the wave's moment.
//
// --- And what she does, where the obvious reading is the wrong one -------------------
//
// **A cut plane's mean `u_z` is not the ship's deflection and using it cost this
// measurement a false answer before it was checked.** The load is applied to every node
// of a slab in proportion to its steel, which puts the right resultant on each station
// and buoyancy on the deck as well as on the shell, so the plating deflects locally
// under it and the mean carries that too. Measured: the second difference of the mean
// `u_z` is **3.8 times** the curvature the stress field carries at the same station,
// which is not something bending can do.
//
// The bending is in `u_x`. A beam's axial displacement is `-z dw/dx`, so a least squares
// of `u_x` against `z` over a cut plane's own rows gives the slope, and a panel bulging
// between frames contributes nothing to it. `--wave` reports that curve and the mean
// `u_z` beside it, because the gap between them *is* the local response and hiding it
// would be worse than publishing it.
//
// Both curves are defined only up to a heave and a trim -- she is floating -- so the
// best-fit line comes off each before they are compared, and the difference is run
// against **two** beams: one with Tier 0's `I` and one with each chain section's own
// measured `EI`, so that "the ends are softer" is separated from "a beam has no shear"
// rather than argued about. Measured on a peak bending deflection of 69.8 mm: the two
// differ by **9.5 mm over the middle two thirds and 33.1 mm at worst**, at the forward
// perpendicular where the section is finest and least like a beam, and putting Tier 1's
// own `EI` into the same integration takes the worst to 25.7 mm. So about a quarter of
// the difference is `EI(x)`. **Shear and the ends are not separated further and this
// brief did not close that.**
//
// --- What mutation testing said about all of it ------------------------------------
//
// Twenty-eight single edits to the code above, compiled from a copy outside the
// repository, with a per-mutant timeout and the verdict taken from the **exit code** as
// well as the `FAIL` lines. None of these twenty-eight needed the exit code -- no
// segmentation fault and no timeout -- and it is kept because §5, §7 and §9 each record
// a run where it was the only thing that scored a mutant correctly.
//
// The first pass killed twenty, and **six of the eight survivors said two things**:
//
//   * **Nothing tested that a Gauss point's weight is a volume.** A mutant giving every
//     point a weight of 1.0 survived the whole file. One bay of the box carries
//     `2(B+H) t L/n` of steel and the weights now have to sum to it, exactly.
//   * **Every fit was taken about the neutral axis**, where the axial term is zero and
//     the cross term `s1 t0` in the normal equations goes with it -- so a sign error
//     there, and a reflected neutral axis, were both invisible. Both are now fitted
//     about the keel as well, on the same samples.
//   * And every load was **hogging**, so a peak taken as the largest *signed* stress
//     rather than the largest magnitude was never wrong and a worst fibre reported as a
//     magnitude never lost a sign. A mirrored sample set and a keel-in-compression
//     assertion close both. A **cut with no depth** was never handed to `fibreStress`
//     either, because nothing this mesher builds produces one; it is constructed.
//
// **The two that are left are equivalent on every input this repository has, and saying
// which is more useful than the score.** Sign-flipping `xi` in the shape function that
// *places* a Gauss point moves the point along the ship and pairs it with another
// point's stress -- and every load `applyBeamLoad` can prescribe gives a stress that is
// uniform along the ship, so nothing sees it. Separating it needs a longitudinal stress
// gradient inside one element, which needs a distributed load. And an off-by-one on the
// section index in `sectionDisplacement` is caught by the next of the three bounds
// tests, which is the same redundancy §9 records for its restraint guard.
//
// Twenty-six of twenty-eight.

}  // namespace sim::section
