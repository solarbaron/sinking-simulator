// SPDX-License-Identifier: MIT
//
// The Tier-2 zone: `StructuralMesh` + a load -> solid-shell elements -> an
// explicit solve -> which panels tore.
//
// `solid_shell.{hpp,cpp}` is a validated locking-free element and
// `plasticity.{hpp,cpp}` a validated return map with ductile failure. Until this
// file existed **neither had a consumer**: nothing built elements over a ship's
// structure and solved. `breach.{hpp,cpp}` sits on the other side waiting for a
// list of failed panel indices. This is the bridge, and its output is exactly
// what `breachesFromFailedPanels` already takes.
//
// It replaces `indentation.{hpp,cpp}` -- rigid-plastic membrane stretching of one
// bay, with no bending, no stiffeners and no relief from the bays around it --
// for the cases that can afford it. It does not delete it: the membrane model
// costs microseconds where this costs core-minutes, and `tests/test_zone.cpp`
// runs the two against each other precisely because two models built from
// different physics agreeing to a factor is worth more than either agreeing with
// itself.
//
// --- 1. Cost is a design input, not a footnote ---------------------------------
//
// Measured (`02-simulation.md` §3): an elastoplastic solid-shell element costs
// **3.1 µs** per step against 273 ns elastic, and the stable step is `t / c_p` --
// **thickness-governed, flat in the in-plane element size**. So for 12 mm plating
// the step is 1.8 µs, 5.5 x 10^5 steps per simulated second, and
//
//     core-seconds per simulated second = 1.7 x elementCount
//
// That 3.1 µs was **7.3** until the step-invariant element forms stopped being
// rebuilt every step -- see `SolveParams::cacheRestForms`, which is the single
// largest cost decision in this file and was found by profiling rather than by
// reading. Every figure below that predates it is 2.4x pessimistic per element.
//
// That number is the whole design. Two consequences follow and both shape this
// file:
//
//   * **The zone must be small, and it is bounded by element *count*, not by
//     area.** Since the step does not care about the in-plane size, area and
//     resolution trade against each other as `elements = area / h^2`. A 200 m^2
//     zone is 80 000 elements at 50 mm and 2 200 at 300 mm -- two hours of wall
//     time on 24 threads against four minutes. Quoting an affordable *area*
//     without its in-plane size says nothing.
//   * **The event is short.** A 6 m/s bow reaching 0.8 m into a side is 0.13 s,
//     not one second. Cost scales with the simulated duration, and the duration is
//     set by the physics rather than chosen.
//
// Put together, at the resolution `MeshParams::subdivision = 4` delivers on the
// reference ferry (0.6 m x 0.175 m elements, four across the 0.70 m between
// longitudinals), a fourteen-panel zone is 224 elements at **380** core-seconds
// per simulated second, and driving a punch 0.23 m into her takes 0.04 s of
// simulated time -- **4.5 s of wall time on 23 workers**, against 15.4 s on one;
//
// **That figure read 900 until it was checked against the sentence it is in.**
// 900 is `224 x 5.5e5 x 7.3 us` -- the pre-`cacheRestForms` cost, the one §1
// above declares superseded and says every figure below it is 2.4x pessimistic.
// The measurement eight words later refutes it without leaving the line: 15.4
// core-seconds over 0.04 s of simulated time is **385**, which is the live 3.1 us
// to within 1%. A prediction and its own measurement, adjacent and 2.36x apart.
// the speedup
// saturates at 3.4x because a step is a barrier and the zone is small, not because
// the element loop fails to thread. Measured, not extrapolated. The same
// 24 m^2 at 50 mm elements is 9 600 elements and forty times that. `Solver`
// measures and reports µs/element/step so the cost is visible rather than
// discovered, and `estimatedCost()` predicts it before the solve starts.
//
// --- 2. Meshing: the geometry limit is the mesher's problem ---------------------
//
// `07-fem-spike-findings.md` §6 limit 1 is the binding one: the ANS interpolation
// is exact only for an element **prismatic through its thickness**, and
// non-parallel faces cost about `90 x (offset/t)^2` in spurious stiffness -- a 5%
// offset is 24% too stiff. The rule it states is "keep the thickness direction
// within a few degrees of the surface normal, and change plate thickness at a
// seam rather than across an element", and both halves land here.
//
// **Thickness direction.** Each mid-surface node carries one normal, the
// area-weighted average of the elements around it, and the element is extruded
// +/- t/2 along it. Where the plating is flat every element is then *exactly*
// prismatic -- offset identically zero, not merely small. Where it is curved the
// nodal normals disagree across an element by the surface's turning angle, and
// that angle **is** the offset ratio: for nodal normals spread by an angle theta,
// `offset/t = theta` to first order. `Patch::worstNormalSpread` measures it and
// `Patch::spuriousStiffness` converts it, so a zone laid over a hard chine reports
// that it is 20% too stiff instead of quietly being so.
//
// The alternative -- extruding each element along its own face normal, which is
// exactly prismatic everywhere -- was rejected because the elements then share no
// nodes at all on a curved surface and the patch falls apart into loose plates.
//
// **Thickness seams.** A node between a 12 mm and a 16 mm strake has no single
// position: the two elements would want it 2 mm apart. Averaging puts a taper
// *inside* both elements, which is the case the limit forbids. So the patch
// **stops at a thickness seam** -- the flood fill will not cross one -- and the
// seam becomes a zone boundary, which is the "change thickness at a seam" the rule
// asks for. `Patch::problems` says when a zone was truncated that way, because a
// truncation that is not reported is indistinguishable from a small ship.
//
// --- 3. Stiffeners: two bounds, and now the member itself ----------------------
//
// The zone was plating only, and the reason is that **the only element in the
// inventory was the solid-shell hex, and there is no way to attach a web to a plate
// with it that is not wrong.**
//
//   * A web sharing one row of nodes along the seam is a **hinge**. It carries the
//     axial and strong-axis bending the stiffener is for, and it has a zero-energy
//     rotation about the seam -- the tripping mode -- which an explicit scheme
//     turns into ringing rather than into an answer.
//   * A web widened to the plate's element size to get two rows of shared nodes
//     has its strong-axis stiffness wrong by `(h_element / t_web)^3`, a factor of
//     3 000 at the sizes here.
//   * Smearing into the plate is the thing `scantlings.hpp` §1 rejects with a
//     measurement: the panel's second moment falls by **130x** at identical area.
//
// What was needed is a multi-point constraint tying an eccentric beam to the
// shell, and that is now `constraint.{hpp,cpp}`: the solid-shell's two nodes
// through the thickness *are* the rotation of the plate's cross-section, so a
// member at through-thickness offset `e` is an exact linear function of that pair
// and needs no rotational degree of freedom anywhere. `Stiffeners::Modelled`
// builds the member out of axial fibres tied that way. What it buys and what it
// still cannot do are in `constraint.hpp`'s own header; the short version is that
// it carries the section's area, neutral axis and second moment exactly -- checked
// against `scantlings::stiffenedSection` -- and that it *over*-restrains tripping
// where the hinge leaves it free.
//
// **Leaving them out entirely is not the neutral choice, and measuring it is what
// showed that.** On the ferry's side under a 2 m punch at 0.35 m of penetration,
// a zone with no stiffeners resists at 6.6 MN where a membrane model spanning the
// real 0.70 m longitudinal spacing says 34 MN -- a factor of five, because with no
// supports the plating spans from one clamped zone edge to the other. The span it
// used was therefore the **zone radius**, a meshing parameter: exactly the failure
// `indentation.hpp` records in its own history, where the size of the hole came out
// a property of the contact radius rather than of the collision.
//
// So `Stiffeners` offers the two bounds the omission leaves, and the member that
// now sits between them:
//
//   * `Ignored` -- the lower bound. Soft by the factor above; conservative for a
//     damage-stability question and useless for a survivability one.
//   * `RigidSupport` -- the upper bound, and still the default. Every plating node
//     a stiffener runs through is pinned, which is what a member far stiffer than
//     the plate it carries does in the limit. Measured on the same zone, the
//     resisting force comes out 23-30 MN against the membrane model's 28 MN on the
//     same 0.70 m span, so the plating between supports is behaving as the membrane
//     model says it should.
//   * `Modelled` -- the member itself, eccentric fibres tied to the plating. It has
//     to land *inside* the bracket, and `tests/test_zone.cpp` asserts that rather
//     than trusting it: a stiffener that came out stiffer than a rigid support or
//     softer than no support at all would be a formulation error that no single
//     run could see.
//
// `RigidSupport` remains the default because it is the setting that makes the two
// models comparable at all -- `indentation.hpp` assumes one bay with its boundaries
// held, which is precisely a rigid support on every stiffener line -- and because
// it is what every published figure in `02-simulation.md` §3 was taken with.
//
// --- 4. Boundary conditions --------------------------------------------------
//
// The patch is cut out of a ship, so its edge is a lie either way. Free lets it
// translate away under the first load; clamped is too stiff, because the plating
// outside the zone really does pull in. Clamped is the default, and **the price is
// paid in zone size rather than in a tuned spring stiffness**: the boundary error
// is a Saint-Venant effect that decays away from the edge, so growing the radius
// makes it go away and the convergence is measurable. `tests/test_zone.cpp`
// measures it. An elastic edge with a stiffness standing for the surrounding
// structure would be softer and closer, and it would also be a coefficient no
// measurement sets -- which is the argument `solid_shell.hpp` uses against
// hourglass control.
//
// The consequence to keep in mind: **a zone one bay wide is a clamped bay**, which
// is the same idealisation `indentation.hpp` makes, and a zone five bays wide is
// not. That is deliberate -- it is what makes the two models comparable at radius
// one and different at radius five.
//
// **There is now a third option, and it is neither a bound nor a coefficient.**
// `HexMesh::prescribed` says what value a pinned DOF takes, and the solver below
// honours it: a perimeter DOF whose prescribed value is non-zero is *driven* there
// rather than held at the meshed position. `coupling.hpp` fills that array from a
// Tier-1 reduced model of the structure round the patch, so the edge follows what
// the rest of the ship is doing. Clamped is then the special case "the surroundings
// said zero", which is what makes the two paths one path.
//
// The static path already read `prescribed`; the explicit one ignored it, silently,
// which meant a boundary condition set on a patch was honoured or dropped depending
// on which solver saw it. Both read it now.
//
// --- 5. What this does not do yet ---------------------------------------------
//
//  1. **Coupling to Tier 1 is static, and it is a linearisation on the way back.**
//     `coupling.{hpp,cpp}` now drives the perimeter from a Craig-Bampton model of
//     the structure round the patch, hands a torn zone back to that model as a
//     mesh with the dead elements removed, and hands a *yielded* one back as a
//     secant knockdown built from the equivalent plastic strain this solver
//     already carries per integration point (`coupling.hpp` §5). This item used to
//     say the second of those could not be done without a tangent operator; at 45%
//     of the way to a tear it measured at 66% of the surroundings' own
//     displacement, a secant closes it to 9%, and a tangent would have made it
//     *worse* than doing nothing at small plastic strain. What is left is that the knockdown is isotropic where a J2
//     secant is not, and that the zone's state and the interface displacement set
//     each other -- so the loop is staggered and a caller decides how many passes
//     it can afford. Coupling to Tier 0 is unchanged and remains the path
//     `promotion.{hpp,cpp}` owns.
//  2. **The indenter is rigid, and it no longer has to be kinematic.** Nodes inside
//     a rectangular footprint are driven; there is still no contact search, no
//     friction, no release, and the striking body still does not crush. What has
//     changed is where the motion comes from: this item used to say a delivered
//     *energy* "would need the striking body's mass and is what `collision.hpp`
//     would eventually supply", and that is now `Drive::Inertial` -- see §6. It
//     supplies neither more nor less than `collision.hpp` already computes.
//  3. **Element deletion, not splitting.** A torn element carries no stress and
//     keeps its mass. The hole is therefore the deleted area, which this reports as
//     whole panels because that is what `breachesFromFailedPanels` consumes -- see
//     `SolveParams::tearFraction`.
//  4. **The GPU path is now the faster one and still cannot be trusted** --
//     `engine/gpu/zone_gpu.{hpp,cpp}`, whose state comes back through
//     `Solver::adopt` below. Re-mapped to one workgroup per element it runs at
//     1.26-2.43x this solver on 23 workers, where one invocation per element ran at
//     0.23-0.68x. What stops it being used is precision: in float it tears a
//     quarter to a half more elements than this double path does, and the torn set
//     is the answer a zone exists to produce. `07-fem-spike-findings.md` §8.
//     This item said "no GPU path" while the same file named `gpu::ZoneGpuSolver`
//     fifty lines further down, and then said "slower" after it had stopped being.
//  5. **No rate dependence in the material**, so the resistance is under-predicted
//     by the 10-30% steel gains at collision strain rates.
//
// --- 6. A collision delivers joules; this used to consume a travel --------------
//
// `collision.hpp` ends an event holding two numbers -- `ImpulseSolution::
// effectiveMass` and an energy `(1 - e^2)/2 m_eff u^2` its penalty force's work
// integral is aimed at -- and the zone had nowhere to put either. It took a punch
// speed and a `stopAt`, so the depth of the hole was an *assumption* rather than a
// consequence, which is the same failure `indentation.hpp` records against its own
// first version. `Drive::Inertial` closes it.
//
// **The entry point is a mass and a velocity, not an energy, and that is a
// decision rather than a convenience.** Three shapes were available:
//
//   * *An energy.* Rejected. This is explicit dynamics: an energy is not a
//     boundary condition, it is a first integral of one, and the zone would have
//     to invent a mass to get a motion back out of it. Two collisions carrying
//     identical joules at different masses do not do identical damage here and
//     the solver can tell them apart -- the plating ahead of a punch sees `rho c v`,
//     280 MPa at 6 m/s, so how fast the joules arrive is part of the answer. A mass
//     invented inside this file to hide that would be exactly the coefficient no
//     measurement sets that §4 rejects for the boundary.
//   * *A force history*, played in from `collision.hpp`'s penalty contact. Rejected,
//     and it is the shape that would really have cost the testability: **a
//     prescribed force is the drive that runs away.** If the structure fails it
//     keeps pushing at the same newtons into nothing and accelerates without
//     bound. It is also a force computed against a *rigid* hull, which §5 of
//     `collision.hpp` says over-states the peak and under-states the penetration,
//     so it is the wrong history for a structure that deforms.
//   * *A mass and a velocity.* Taken. It is what `collision.hpp` has, it reduces to
//     the energy exactly -- `impactSpeed(energyLost, effectiveMass)` gives the
//     `u sqrt(1 - e^2)` whose kinetic energy *is* that energy -- and it carries the
//     deceleration law as well, which an energy alone cannot.
//
// **Its speed cannot run away; its travel can, and the difference is the whole
// design.** Nothing here does work on the punch: the plating grips it, and a
// stretched plate pulls back along the travel rather than pushing forward. So the
// total energy of (punch + patch) is non-increasing -- internal forces are
// conservative, the grip and the return map and `damping` are dissipative, and a
// clamped edge does no work -- and therefore `1/2 m v^2 <= 1/2 m v_0^2` for the
// whole run. The punch can never be going faster than it arrived, whatever the
// structure does, and `tests/test_zone.cpp` asserts that on a run that tears
// because tearing is the case the bound exists for.
//
// **That bound is on the speed and not on the distance, and assuming otherwise is
// a mistake this file made first.** A punch that has perforated the plating under
// it meets no force at all and coasts: twice the energy needed to tear the
// reference strip put it 24.9 m past the plating and ended the run on `maxSteps`.
// Nothing is wrong there -- the zone simply has nothing further to say, because
// the rest of the energy goes into structure this patch does not contain, and the
// punch cannot re-engage what it is not already touching (there is no contact
// search; §5 item 2) -- but a solver that discovers it after two million steps is
// no use. So **an inertial drive is refused unless `Indenter::stopAt` or
// `SolveParams::duration` bounds it**, and the cap's meaning changes: on a
// prescribed punch it *is* the answer, and here it is the edge of what this patch
// can be asked about. A run that ends on it leaves `indenterKinetic` unspent,
// which is `ImpactDamage::energyUnspent` on the membrane path and means exactly
// the same thing -- the hole would have been bigger and the caller is looking at a
// bound rather than a result.
//
// **Four requests are refused outright rather than answered wrongly**, all of them
// with `problems` and no steps taken: no striking mass, no approach speed, nothing
// to bound the run, and a footprint that grips nothing. Falling back on the
// prescribed drive for any of them would hand the caller exactly the assumed
// penetration this entry point exists to remove -- and the last is not
// hypothetical, since a punch that grips nothing is never decelerated by anything
// and would coast to the cap reporting a perforation it never made.
//
// **Two things widen the speed bound and both are named rather than assumed away.**
// A `Preload` hands the patch `initialStrainEnergy` it did not pay for, and a
// driven perimeter keeps putting `boundaryWork` in; either can be released into the
// punch, so the honest statement on such a patch is
// `1/2 m v^2 <= 1/2 m v_0^2 + initialStrainEnergy + boundaryWork`. Both are in
// `SolveResult` already, so the wider bound is checkable and not merely stated.
//
// The step's own arithmetic gives a second, exact statement that needs no bound at
// all: the grip conserves momentum, so `mass * (v_before - v_after)` is precisely
// `force * timestep` every step. That is a closed form to floating point rather
// than to an order, and it is what the tests assert the grip against.
//
// **What is given up.** Three things, and none of them is the speed bound above:
//
//   * **The cost is not knowable before the run, and it is not bounded by the
//     prescribed run's either.** A prescribed punch takes `stopAt / speed` seconds
//     and `estimatedCost()` turns that into wall time. An inertial one takes as
//     long as it takes, and the tail is the expensive part: the step is fixed by
//     the material, so a striker that has nearly stopped needs many steps per
//     centimetre. On the ferry's own patch the same 0.22 m cost 45 970 steps by
//     energy against 21 290 by travel, and a striker that perforated and then
//     coasted at 0.5 m/s took 216 000. **`stopAt` bounds the reach and `duration`
//     bounds the cost; they are not the same bound and a run that can perforate
//     wants both.**
//   * **The answer at a given travel becomes path dependent.** The punch is fast at
//     first and slow at the end, where a prescribed one is neither, so the two
//     drives do not pass through the same states. They agree wherever the run is
//     quasi-static and they are entitled to differ where it is not; how far apart
//     they land is measured in `tests/test_zone.cpp` rather than assumed small.
//   * **Rebound is still not modelled.** The punch grips, so at arrest the plate
//     would spring back and drag it out; the run ends there instead. The reported
//     penetration is therefore the deepest the punch got, which is the number a
//     hole is read from, and not a residual dent.
//
// **Where it errs, and in which direction.** The punch is rigid, so every joule
// goes into the struck plating. That is the same bias `indentation.hpp` records
// against its own high-energy holes, for the same reason, and it makes this an
// **upper bound on the struck side's damage** -- tight when the striking bow is far
// stronger than what it hits, loose when the two are comparable, which for two
// ship sides is the usual case. Reducing the delivered energy by a share to stand
// for the bow is *not* offered, because that share is a coefficient nothing here
// measures. What has changed is that the honest fix is now representable at all:
// the punch is a body with a mass, so giving it a crushing characteristic is one
// force law between it and the plating. It could never have been put on a
// prescribed motion, which has no force to share.
//
// **This does not duplicate `indentation.hpp`.** That model has been energy-driven
// from the start and stays: it spends a budget *outward across panels*, bay after
// bay, with no dynamics in it at all, and costs microseconds. This spends one
// striking body's kinetic energy *through a single solve*, and the two answer
// different questions -- "how far does the damage reach" against "what does this
// patch do". `tools/zone_probe` now runs them against each other on energy as well
// as on force.
//
// SI units, body frame per CLAUDE.md.
#pragma once

#include "constraint.hpp"
#include "plasticity.hpp"
#include "scantlings.hpp"
#include "solid_shell.hpp"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace core {
class JobSystem;
}

namespace sim::zone {

// --- Meshing -------------------------------------------------------------------

// What holds the patch's edge. See §4 -- `Clamped` is the answer; `Free` exists so
// the tests can ask the questions only an unconstrained patch answers (a rigid
// translation carries no force; nothing at rest starts moving).
enum class Edge { Clamped, Free };

// What the stiffeners *inside* the zone do. See §3: they are not meshed, and the
// two settings are the two bounds that fact leaves.
enum class Stiffeners {
    // Not there at all. The plating then spans from one clamped zone edge to the
    // other, which is soft by a large factor -- measured on the ferry's side at
    // 0.35 m of penetration under a 2 m punch, **five times** too soft against a
    // membrane model spanning the real 0.70 m longitudinal spacing (6.6 MN against
    // 34). This comment said *seven* times until the membrane figure was
    // re-derived: `indentationForce` on that span under a 2 m contact is 12.05 MN
    // per bay and the punch covers 2.0/0.70 = 2.86 of them, so 34.4 MN and a ratio
    // of 5.15 -- which is what §3 and `02-simulation.md` §3 both already said.
    // Nothing tests a comment. Worse, the
    // span it does use is the *zone radius*, so the answer depends on a meshing
    // parameter. That is the failure `indentation.hpp` records in its own history,
    // where the hole came out a property of the contact radius rather than of the
    // collision.
    Ignored,
    // Every plating node lying on a stiffener line is pinned. A stiffener really
    // is far stiffer than the plating it carries, so this is the *upper* bound and
    // `Ignored` is the lower one; between them they bracket the answer, which is a
    // more honest thing to publish than either. It is also exactly the boundary
    // `indentation.hpp` assumes -- one bay, boundaries held -- so the two models
    // are comparing the same idealisation.
    //
    // It needs `subdivision` of at least 2, and really 3 or 4: the stiffener lines
    // are the panel seams, so at subdivision 2 a panel is left with one free node
    // and cannot deform. `Patch::problems` says so when it happens.
    RigidSupport,
    // The member itself: axial fibres at the profile's own through-thickness Gauss
    // stations, each tied to the plating by `constraint.hpp`'s multi-point
    // constraint. Nothing is pinned, so the plating between stiffeners spans the
    // real distance *and* the supports deflect, which is the whole point -- the
    // two bounds above differ by a factor of two to five on the reference ferry and
    // the answer is somewhere in there.
    //
    // It is not the default. `RigidSupport` is what every figure `02-simulation.md`
    // §3 publishes was taken with, and changing the default would silently move
    // them; a caller that wants the member asks for it.
    Modelled,
};

struct MeshParams {
    double radius = 4.0;         // m from the impact point, on panel centroids
    // Elements across each panel edge, both ways. **An integer, and the same in
    // both directions, on purpose**: a shared edge is then divided identically by
    // the panels either side of it whichever way round each of them numbers its
    // corners, so the mesh is conforming without a seam-matching pass. A per-edge
    // count taken from a target element size is not conforming, because a quad's
    // opposite edges have different lengths.
    //
    // Four, because the stiffener lines are the panel seams: at 2 a panel is left
    // with a single free node in the middle and cannot deform at all. Four leaves
    // a 3 x 3 free grid per bay, which is 0.6 m x 0.175 m on the reference ferry.
    int subdivision = 4;
    PanelRole role = PanelRole::Shell;
    Edge edge = Edge::Clamped;
    Stiffeners stiffeners = Stiffeners::RigidSupport;

    // Two panels are neighbours when they share an edge to this tolerance. Matches
    // `BreachParams::weldEpsilon` and the CSG's weld.
    double weldTolerance = 1e-6;  // m

    // The flood fill refuses a neighbour folded further than this from the panel it
    // came from, so a zone cannot turn through a chine, round the keel and come up
    // the other side. It is a fold *per panel*, not a total, so a gently curved
    // bilge is still followed.
    double foldLimit = 0.785;  // rad

    // Refuse to cross a plate thickness seam. See §2: crossing one puts a taper
    // inside an element, which is the geometry the element is known to be wrong on.
    //
    // **Setting it false does not mesh the second thickness** -- it meshes the
    // whole patch at the struck panel's thickness and says so, because one patch
    // has one extrusion distance. It exists so a caller can trade a known
    // thickness error against a truncated zone knowingly, which is a different
    // thing from being given one without being told.
    bool singleThickness = true;

    // Reported, not enforced: above this spread of nodal normals across one
    // element the patch is being laid over curvature the element cannot carry.
    // 0.05 rad is ~3 degrees and ~22% spurious bending stiffness.
    double normalSpreadWarning = 0.05;  // rad

    // Which way is out. Zero derives it, by taking the struck panel's normal away
    // from the structure's own centroid -- right for shell plating, arbitrary for a
    // bulkhead, where "out" is a choice rather than a fact and the caller should
    // make it.
    Vec3 outward{0, 0, 0};
};

// A meshed patch of plating: the elements, what they came from, and how good the
// geometry turned out to be.
struct Patch {
    solidshell::HexMesh mesh;

    std::vector<int>    panelOf;      // per element -> index into StructuralMesh::panels
    std::vector<double> elementArea;  // per element, m^2 of mid-surface
    std::vector<int>    panels;       // the distinct panels covered, ascending
    std::vector<double> panelArea;    // per entry of `panels`, m^2 meshed of it
    // Per node: 1 on the +axis face, 0 on the -axis face. The indenter touches the
    // outer face only, and the pair (2m, 2m+1) is one mid-surface point.
    std::vector<std::uint8_t> outerFace;

    StructuralMaterial material;
    double thickness = 0;   // m, one value: the patch stops at a thickness seam
    double area = 0;        // m^2 of mid-surface meshed
    double mass = 0;        // kg of steel in the patch

    // The outward normal at the struck panel. The indenter travels along -axis;
    // `right` and `up` complete a right-handed in-plane frame, and the indenter's
    // footprint is rectangular in them.
    Vec3 axis{}, right{}, up{}, centre{};
    int  struckPanel = -1;

    // Geometry quality, measured rather than assumed. See §2.
    double worstNormalSpread = 0;   // rad, over elements: max |n_node - n_element|
    double spuriousStiffness = 0;   // 90 * spread^2, the fraction of excess bending stiffness
    // How many elements are past `MeshParams::normalSpreadWarning`. The worst
    // alone is misleading: on the reference ferry a zone crossing the shoulder at
    // z = 4.2 m has three bad element rows out of thirty, and whether that matters
    // depends on whether they are where the punch is.
    int distortedElements = 0;
    double worstJacobian = 0;       // smallest det J over the patch; <= 0 is fatal
    double worstAspect = 0;         // longest / shortest mid-surface edge

    // How much of the patch is free to move: nodes not pinned by the perimeter or
    // by a stiffener line, over the total. A patch whose stiffener lines have eaten
    // its interior cannot deform, and the number says so before the solve does.
    double freeFraction = 0;
    int stiffenerNodes = 0;   // pinned because a stiffener runs through them

    // The members themselves, under `Stiffeners::Modelled` and empty otherwise.
    // Fibres tied to the plating by `constraint.hpp`; they carry no degrees of
    // freedom of their own, so they change the patch's stiffness and mass without
    // changing its node count.
    constraint::Stiffening stiffening;

    // Smallest stable explicit step over the patch, computed here rather than in
    // the solver because it is what makes the cost knowable *before* the run: it
    // is 0.5 ms of power iteration per element, which is a promotion cost like the
    // stiffness formation beside it and not a per-step one.
    //
    // **The fibres are in it.** A tie amplifies a fibre's stiffness by the square
    // of its weight while its mass arrives unamplified, so whether a stiffener
    // shortens the step is a measurement rather than an assumption; it is taken
    // here, exactly, because a rank-one stiffness against a diagonal mass has a
    // closed-form largest eigenvalue.
    double criticalTimestep = 0;    // s
    // What the plating alone would have allowed, so the fibres' share of the cost
    // is visible instead of being discovered as a slower run.
    double platingTimestep = 0;     // s

    std::vector<std::string> problems;

    std::size_t elementCount() const { return mesh.elementCount(); }
    std::size_t nodeCount() const { return mesh.nodeCount(); }
    bool empty() const { return mesh.elementCount() == 0; }
};

// Mesh a patch of plating round `impact`. Always returns something: an impossible
// request yields an empty patch and a full account in `problems`, in the same
// spirit as `makeStructuralMesh`.
Patch buildPatch(const StructuralMesh& structure, const Vec3& impact,
                 const MeshParams& params = {});

// --- Driving it ----------------------------------------------------------------

// Where the punch's motion comes from. See §6 -- the two are the same rigid
// gripping punch and differ only in what sets its speed.
enum class Drive {
    // A travel is prescribed: the punch moves at `speed` whatever the plating
    // does, as if infinitely heavy. The original, and still the default, because
    // it is the drive a *convergence* question wants -- two meshes compared at one
    // penetration are comparable, where two meshes compared at one energy stop at
    // two different penetrations and the difference has both effects in it.
    Prescribed,
    // A striking body of `mass` arrives at `speed` and the plating decelerates it.
    // The run ends when it stops, having spent `1/2 mass speed^2` on the patch --
    // so the travel is an output and the energy is the input, which is the way
    // round a collision actually delivers it.
    Inertial,
};

// A rigid flat punch, rectangular in the patch's own in-plane frame, travelling
// along -axis. Nodes on the **outer face only** are driven, so the plate is free
// to thin under the punch -- which is the whole reason the element keeps a
// thickness stretch mode. It grips: a driven node is held to the punch rather than
// released, which is right for a monotonic push and wrong for rebound.
// A footprint with a non-positive half-extent is **no punch at all**, which is
// how a patch is asked to carry no load. Leaving the footprint and setting the
// speed to zero is a different thing -- a stationary punch still holds the nodes
// it touches, and holding them is exactly what stops an unloaded patch from
// translating freely.
struct Indenter {
    double halfLength = 1.0;  // m along `Patch::right`
    double halfWidth = 1.0;   // m along `Patch::up`
    // m/s. Constant under `Prescribed`; the speed the striking body *arrives* at
    // under `Inertial`, after which it is the solver's to decide.
    //
    // Zero is a *stationary punch* under `Prescribed` -- it still holds the nodes
    // it touches, which is what stops an unloaded patch translating freely -- and
    // is **refused** under `Inertial`, where it is a body carrying no energy. It is
    // refused rather than allowed to stop on the first step, because that would
    // make a run's termination turn on an exact floating-point zero and it does not
    // reach one: the internal force at the rest configuration is a rounding residue,
    // and a striker released at zero picks up ~1e-29 m/s from it and runs to
    // `maxSteps` having travelled 1e-33 m. Measured, not supposed.
    double speed = 6.0;
    // m of penetration at which the run ends. Under `Inertial` it is a **cap and
    // not the answer** -- the energy is meant to end the run -- so a run that hits
    // it says so in `SolveResult::problems`, for the same reason
    // `ImpactDamage::energyUnspent` exists. Zero removes the cap; under `Inertial`
    // that is only allowed when `SolveParams::duration` bounds the run instead,
    // because a perforating punch coasts and `maxSteps` is a budget rather than a
    // bound. See §6.
    double stopAt = 0.5;
    // Seconds over which the speed rises linearly from zero. Zero starts the punch
    // at full speed, which is a genuine step in velocity: the plate ahead of it
    // sees rho*c*v, 280 MPa at 6 m/s, and rings for as long as the run lasts. A
    // ramp is a statement about the striking body's own compliance rather than a
    // numerical fudge, and the energy account holds either way.
    //
    // **`Prescribed` only.** An inertial punch's speed is its own state, and a body
    // with a mass genuinely does arrive travelling; a ramp on it would be a
    // statement about a compliance the mass has now made representable properly.
    // Setting one anyway is reported rather than quietly honoured or dropped.
    double rampTime = 0.0;    // s

    Drive drive = Drive::Prescribed;
    // kg, `Inertial` only, and it is the *effective* mass of the striking pair at
    // the contact rather than either ship's -- `collision.hpp`'s
    // `ImpulseSolution::effectiveMass`. A published damage table in this repo once
    // used one whole ship's mass as though all of it reached the plating; the
    // reduced mass is the quantity that does not make that mistake.
    double mass = 0.0;

    // There was a `kineticEnergy()` here, `1/2 mass speed^2`, and it is gone. It
    // had no consumer -- `SolveResult::indenterEnergy` is the same quantity, taken
    // from the state the solver actually ran with -- and mutation testing found it
    // by dropping the half and surviving the entire suite. An accessor nothing
    // calls is an accessor nothing tests; the table in CLAUDE.md already records
    // that exact shape shipping in a well-tested commit.
};

// The speed at which a body of `mass` carries `energy`: sqrt(2E/m). What to set
// `Indenter::speed` to when the collision is known by its joules rather than by a
// closing speed -- `collision.hpp` names both halves, so
//
//     indenter.mass  = solution.effectiveMass;
//     indenter.speed = zone::impactSpeed(solution.energyLost, solution.effectiveMass);
//
// makes `1/2 m v^2` exactly the `(1 - e^2)/2 m_eff u^2` that collision took out of
// the two ships' motion. Zero for a non-positive mass or energy, because there is
// no such body; the caller gets a punch that does not move rather than an
// infinity.
double impactSpeed(double energy, double mass);

// The stress the hull girder is *already* carrying through this patch when it is
// promoted, and the thing §5 item 1 records as missing.
//
// A patch of side shell in a hogging ship is not unstressed. On the reference
// ferry poised on a 3 m crest the deck carries 84 MPa against a 355 MPa yield
// (`02-simulation.md` §3), so a zone told it starts from zero starts with a
// quarter of a capacity it does not have. Uniaxial along the ship's x and linear
// in height, which is exactly and only what a beam model knows:
//
//     sigma_xx(z) = stress + gradient * (z - reference)
//
// **It is imposed as a strain, not as a stress.** The rest configuration handed
// to the elements is the meshed one with the corresponding displacement field
// taken back out, so the stress at step zero comes out of the same validated
// constitutive path as everything else and no new code goes inside the element.
// The field is the exact elasticity solution for that stress state -- see
// `applyPreload` -- so it is compatible and equilibrated, and a clamped patch
// carrying it and nothing else does not move.
//
// **Where it is not right.** The state is uniaxial, so it is traction-free only
// on a surface whose normal is perpendicular to x. That is every shell, deck and
// longitudinal-bulkhead panel on a parallel-body ship, and it is *not* a
// transverse bulkhead or a steeply raked stem: a panel whose normal leans by
// `phi` out of the athwartships plane is left with `sigma sin^2(phi)` of
// unbalanced traction on its own face. `promotion.hpp` measures that angle and
// declines to pre-load a patch that is too far out, because a pre-load nobody
// checked is worse than none.
struct Preload {
    double stress = 0;     // Pa at z = `reference`, tension positive
    double gradient = 0;   // Pa per m of height
    double reference = 0;  // m above the baseline; the girder's neutral axis

    bool active() const { return stress != 0.0 || gradient != 0.0; }
    double at(double z) const { return stress + gradient * (z - reference); }
};

struct SolveParams {
    Indenter indenter;

    // What the hull girder is already putting through this patch. Default is
    // none, which is the old behaviour and is right for a patch nothing is
    // bending.
    Preload preload;

    // False solves the co-rotational *elastic* element -- 273 ns against 3.1 µs,
    // and the path the geometric tests use, where plasticity would only add noise
    // to an identity.
    bool plastic = true;

    // Whether the stiffeners' fibres can tear -- `constraint.hpp` §2b. True is the
    // model. **False is the control**, and it is kept for one reason: it is the
    // un-conservative model that shipped before the fibres carried damage, and the
    // size of the error it makes is a measurement that has to stay re-runnable
    // rather than re-readable. It reaches the fibres and nothing else: the plating's
    // own tearing is untouched, so an A/B isolates the stiffeners exactly. Same
    // spirit as `coupling::Modulus::Tangent`, which is kept as the control for the
    // modulus the docs used to prescribe.
    bool fiberFailure = true;

    // Seconds over which a driven perimeter -- `Patch::mesh.prescribed`, which
    // `coupling.hpp` writes -- rises from the meshed position to its full value.
    // Zero imposes it on the first step.
    //
    // The ramp is a **smoothstep**, not the linear one `Indenter::rampTime` uses,
    // and the difference is not cosmetic. An indenter's linear ramp ends in a real
    // velocity discontinuity because a striking body genuinely arrives travelling;
    // a boundary drive is a statement about where the surrounding structure *is*,
    // so a velocity step at either end of it is purely numerical and rings the
    // patch for the rest of the run. `3s^2 - 2s^3` starts and ends at zero
    // velocity, so nothing is excited that the drive did not ask for.
    double edgeRamp = 0.0;  // s

    double damping = 1.0;   // velocity multiplier per step; 1.0 is none
    double timestep = 0.0;  // s; 0 takes `solidshell::criticalTimestep` over the patch
    double timestepSafety = 0.9;
    double duration = 0.0;  // s; 0 runs until the indenter reaches `stopAt`
    int    maxSteps = 4000000;

    // A panel is reported torn once this fraction of its meshed area has been
    // deleted. **The rule matters and neither end of it is right**: `breach.hpp`
    // takes panel indices and opens the whole panel, so reporting a panel on its
    // first dead element over-states a one-element-wide slit by the subdivision
    // squared, and requiring every element under-states a tear that has crossed the
    // bay. Half is the honest middle; the real fix is a breach interface that takes
    // an area, and that is a change to `breach.hpp` rather than to this.
    double tearFraction = 0.5;

    int historyStride = 0;   // steps between `SolveResult::history` samples; 0 = none

    // Optional. The element loop is embarrassingly parallel -- each element writes
    // only its own force slot and the nodal gather is by CSR, so the answer is
    // **bit-identical** to the serial one whatever the worker count, which
    // `tests/test_zone.cpp` asserts.
    core::JobSystem* jobs = nullptr;

    // Keep each element's `solidshell::RestForms` for the life of the solve
    // instead of rebuilding it every step.
    //
    // **This is the largest single cost decision in the file**, and it was found by
    // profiling rather than by reading: those matrices are a function of the rest
    // configuration alone, nothing here moves the rest configuration, and forming
    // them was half the element kernel. On is **2.0x faster end to end** on the
    // ferry's own patch, serial, and the two answers are **bit-identical** --
    // asserted in `tests/test_zone.cpp`, because a cache that is merely close is a
    // cache that has silently changed the physics.
    //
    // The obvious objection is memory: 12.0 kB per element, against the 4.6 kB
    // per-element stiffness whose L3 cliff `02-simulation.md` §3 measures at
    // ~6 500 elements. **It does not bite, and that was measured rather than
    // argued.** The cache still wins 1.62x at 17 800 elements and 214 MB, four
    // times past this machine's last-level cache, because the arithmetic it
    // removes is 2.2 µs per element per step while streaming 12 kB at the
    // saturated 29 GB/s that section measures is 0.41 µs. There is no crossover to
    // find; the ratio is five to one in the cache's favour even entirely from
    // DRAM.
    //
    // False exists as the control the bit-identity assertion needs, and because a
    // caller meshing a zone far past the sizes anyone can afford to solve should
    // be able to decline 200 MB.
    bool cacheRestForms = true;
};

struct Sample {
    double time = 0, penetration = 0, force = 0;
    double work = 0, strainEnergy = 0, dissipation = 0, kinetic = 0;
    // The punch's own speed. Constant on the prescribed drive and the whole story
    // on the inertial one, where it is what turns a travel history into an energy
    // history -- see §6.
    double indenterSpeed = 0;
    int tornElements = 0;
};

struct SolveResult {
    bool   completed = false;   // reached `stopAt` or `duration` rather than `maxSteps`
    int    steps = 0;
    double time = 0;            // s simulated
    double timestep = 0;        // s

    double penetration = 0;     // m the indenter has travelled
    double force = 0;           // N resisting, now
    double peakForce = 0;       // N

    // The energy account. `work` is what the indenter put in; the rest is where it
    // went. They balance to the integrator's order, which is the coupling check
    // `06-roadmap.md` asks every coupling for.
    double work = 0;               // J
    double strainEnergy = 0;       // J, recoverable
    double dissipation = 0;        // J, plastic
    double kinetic = 0;            // J
    double dampingLoss = 0;        // J removed by `SolveParams::damping`
    // What a `Preload` put in before the run started. Zero without one, so the
    // balance below is unchanged for every case that has none -- but a
    // pre-stressed patch stores 84 MPa's worth of energy the indenter never paid
    // for, and an account that did not subtract it would report the solver
    // inventing 5 kJ on its first step.
    double initialStrainEnergy = 0;  // J
    // What a driven perimeter put in, by the same reckoning as the indenter's
    // `work`: the impulse needed to hold each driven DOF on the boundary's motion,
    // times that motion. Kept separate from `work` because `work` answers "what did
    // the indenter spend" and a coupled zone has a second agent on it; zero without
    // a drive, so the balance below is unchanged for every case that has none.
    double boundaryWork = 0;         // J
    int drivenEdgeDof = 0;           // DOF the perimeter drive actually moves
    // **How far this is from closing at ship scale, which was never written down.**
    // The unit fixtures close it to a fraction of a percent, and that is what the
    // suite asserts. On the ferry's own patch it is **-9.5% before anything tears
    // and -102% once it does** -- dissipation 4.638 MJ against work 2.755 MJ, so
    // the dissipation alone exceeds the work. The mechanism is already diagnosed
    // (the co-rotational small-strain measure, which is why the residual does not
    // move with the timestep); what was missing is its *size* on a real run, and a
    // reader who had only seen the unit figures would reasonably take this as a
    // primary check on a tearing solve. It is not one.
    //
    // `indenterResidual()` below is, on exactly those runs: 0.02-0.16% where this
    // is 25-102% open, because it involves only the punch's kinematics and the
    // force actually reported, and it converges with the step where this does not.
    double energyResidual() const {
        return (work + boundaryWork) -
               ((strainEnergy - initialStrainEnergy) + dissipation + kinetic + dampingLoss);
    }

    // --- The striking body's half of the same account, `Drive::Inertial` only ----
    //
    // The account above says where the energy the indenter spent *went*. These say
    // where it came *from*, which a prescribed punch cannot answer -- it has an
    // unbounded supply by construction. All zero on the prescribed drive, so
    // nothing that has none of these changes.
    double indenterMass = 0;      // kg, the striking body's effective mass
    double indenterSpeed = 0;     // m/s, its travel speed now
    double indenterEnergy = 0;    // J, 1/2 m v0^2 -- what it arrived with
    double indenterKinetic = 0;   // J, 1/2 m v^2  -- what it still has
    // True when the run ended because the striker stopped, which is the *intended*
    // ending: the energy was spent and the penetration below is what it bought.
    // False on a run that hit `stopAt`, `duration` or `maxSteps` first, and those
    // are travel caps rather than answers.
    bool indenterArrested = false;

    // What the striker lost, against what the indenter spent. They differ by the
    // grip: bringing the punch and the plating it holds to a common speed each step
    // is a perfectly inelastic collision, and the punch's own share of that loss is
    // `1/2 m (dv)^2` summed over the run -- second order in the step, so this is
    // **first order in dt overall** and goes to zero under refinement. That makes
    // it a convergence check with a rate rather than a tolerance, which is more
    // than `energyResidual()` can offer: that one is limited by the co-rotational
    // strain measure and does *not* move with the step.
    double indenterResidual() const { return (indenterEnergy - indenterKinetic) - work; }

    int yieldedElements = 0;
    int tornElements = 0;
    double tornArea = 0;             // m^2 of deleted element
    std::vector<int> tornPanels;     // ascending indices into StructuralMesh::panels
    double tornPanelArea = 0;        // m^2 of those panels, as `breach.hpp` will see it

    // The stiffeners' half of the same account. A fibre fails on axial damage --
    // `constraint.hpp` §2b -- and a failed one carries no stress and no stiffness for
    // the rest of the run. Reported beside the plating's rather than folded into it
    // because the two fail at very different strains for a reason that is *not* their
    // failure strains: those are within 12% of each other. It is Rice-Tracey. A
    // fibre's triaxiality is exactly the reference, so its multiplier is exactly 1;
    // a plate element under a punch is not, and measured on the reference strip the
    // first plate point tears at 7.8% of its own regularised failure strain while the
    // first fibre tears at 100.0% of its.
    int tornFibers = 0;
    double tornFiberVolume = 0;      // m^3 of stiffener steel that has gone

    // Cost, measured. The figures are for printing and for `estimatedCost`, never
    // for asserting on: `test_plasticity.cpp` records what a tight timing assertion
    // costs on a shared machine.
    double wallSeconds = 0;
    double microsecondsPerElementStep = 0;
    // Whether the step-invariant element forms were kept. Reported because it is a
    // 3.6x cost decision and `SolveParams::RestFormsCache::Auto` takes it without
    // being asked.
    bool cachedRestForms = false;

    // Where the wall time went, per phase, summed over every step.
    //
    // **This exists because Amdahl's law decides whether an accelerator is worth
    // building, and the element kernel's own throughput does not.** A GPU element
    // solver moves `element` and nothing else unless the rest moves with it, so
    // the ratio between these five is the ceiling on what one can buy. They are
    // taken with `steady_clock` around each phase; the clock costs ~25 ns against
    // phases of hundreds of microseconds, so the instrument does not disturb what
    // it measures. `other` is the residue -- the punch, the history sampling and
    // the loop itself -- and is reported rather than dropped so the parts sum to
    // the whole.
    struct Profile {
        double element = 0;    // s in the per-element force/plastic kernel
        double gather = 0;     // s in the CSR nodal gather and the per-element reduction
        double integrate = 0;  // s advancing velocity and position, and the punch
        double energy = 0;     // s in `accumulateEnergy` -- a reduction over every Gauss point
        double other = 0;      // s in the rest of `step()`
        double total() const { return element + gather + integrate + energy + other; }
    } profile;

    std::vector<Sample> history;
    std::vector<std::string> problems;
};

// What a solve will cost before it is run, in **core-seconds of simulated
// second**, from the element count, the patch's own critical timestep and the
// measured per-element cost. Exposed because a zone that cannot be afforded should
// be refused at promotion rather than discovered at the end of a run.
double estimatedCost(const Patch& patch, bool plastic = true);

// The explicit solver. Central-difference (symplectic Euler) time integration on a
// lumped mass, which is what `02-simulation.md` §3 specifies for Tier 2 and what
// `fem.cpp` already does for tets; the element loop and the CSR nodal gather are
// the same shape as `fem.cpp`'s for the same reason -- a fixed accumulation order
// is what makes the answer independent of the worker count.
class Solver {
public:
    Solver(const Patch& patch, const plasticity::Material& material,
           const SolveParams& params = {});

    // One step. False once the run is over.
    bool step();
    // Steps until the indenter reaches `stopAt`, `duration` elapses, or `maxSteps`.
    const SolveResult& run();

    const SolveResult& result() const { return result_; }
    double timestep() const { return result_.timestep; }

    // Current nodal state, 3 per node, in the patch's node numbering.
    const std::vector<double>& position() const { return position_; }
    const std::vector<double>& velocity() const { return velocity_; }
    const std::vector<double>& nodalMass() const { return mass_; }
    // Per-element plastic history, for a test that needs to see inside.
    const std::vector<solidshell::ElementPlasticState>& elementState() const { return plastic_; }
    // Per-fibre plastic history and the fibres' rest lengths, for the same reason.
    const std::vector<constraint::FiberState>& fiberState() const { return fiber_; }
    const constraint::RestFibers& restFibers() const { return fiberForms_; }
    // Strain energy stored in the stiffeners alone, J. `SolveResult::strainEnergy`
    // is the plating and this together; separating them is what lets a test compare
    // the stiffener's contribution against a closed form without the plating's own
    // discretisation error in the way.
    double fiberStrainEnergy() const { return fiberEnergy_; }

    // Displacement of the patch's most-displaced node from the **rest**
    // configuration, m. A patch under no load must keep this at zero, which is a
    // statement no energy total can make.
    //
    // A `Preload` moves the rest configuration, so a pre-loaded patch reports the
    // pre-load's own displacement here from the moment it is built and never
    // returns to zero. "Has it moved since it was promoted" is `position()`
    // against the `position()` it started with, and that is what a test of a
    // pre-loaded patch has to ask.
    double largestDisplacement() const;

    // Volume-averaged Cauchy stress over the patch, Voigt [xx, yy, zz, xy, yz, zx]
    // in the body frame, Pa. A `Preload` is a statement about exactly this and
    // nothing else in the result makes it visible: an energy total cannot tell a
    // pre-stressed patch from an unstressed one, and neither can a displacement.
    // Torn elements are averaged in as the zero they carry.
    void meanStress(double out[6]) const;

    // The stress-free configuration the elements are measuring strain from. It is
    // the meshed geometry unless a `Preload` moved it -- which is the whole of how
    // a pre-load is applied, so it is worth being able to look at.
    const std::vector<double>& rest() const { return rest_; }

    // The panels at least `fraction` of whose meshed area has been deleted, as
    // ascending indices into `StructuralMesh::panels`. `SolveResult::tornPanels`
    // is this at `SolveParams::tearFraction`; it is exposed separately because the
    // threshold is a judgement and a caller -- or a test -- is entitled to ask what
    // the other answers would have been without paying for the solve again.
    std::vector<int> tornPanelsAt(double fraction) const;

    // --- What an accelerator needs ---------------------------------------------
    //
    // The constructor lumps the mass, builds the CSR node->element adjacency, finds
    // the nodes the punch holds and forms every element's `RestForms`. A GPU back
    // end needs all four, and **reimplementing them is exactly where two paths stop
    // being the same solver** -- a second mass lumping or a second adjacency order
    // would be a difference nobody could see in a force. So they are exposed and
    // `gpu::ZoneGpuSolver` uploads these very arrays.
    const std::vector<std::uint32_t>& adjacencyOffset() const { return adjacencyOffset_; }
    const std::vector<std::uint32_t>& adjacencyEntry() const { return adjacencyEntry_; }
    const std::vector<std::uint32_t>& drivenNodes() const { return driven_; }
    const std::vector<std::uint8_t>&  pinnedDof() const { return pinned_; }
    // Empty when `SolveParams::cacheRestForms` is false.
    const std::vector<solidshell::RestForms>& restForms() const { return forms_; }

    // Take a state computed elsewhere -- a GPU run -- and continue as if this
    // solver had produced it, so that the energy account, the tearing rules and the
    // panel reporting are the ones already validated rather than a second copy.
    // `steps` and `time` advance by what the other path ran.
    void adopt(const std::vector<double>& position, const std::vector<double>& velocity,
               const std::vector<solidshell::ElementPlasticState>& state, int steps, double time,
               double penetration, double work, double dissipation);

    // Give every node the same velocity. The only way to ask a patch for a rigid
    // body motion, which must carry no force at all however far it has travelled --
    // the statement the co-rotational formulation exists to make, and one that a
    // deflection test cannot.
    void translate(const Vec3& velocity);

private:
    void applyPreload();
    void computeForces();
    void accumulateEnergy();
    void collectTorn();

    const Patch* patch_;
    plasticity::Material material_;
    SolveParams params_;
    SolveResult result_;

    std::vector<double> rest_, position_, velocity_, mass_, force_;
    std::vector<double> elementForce_;   // 24 per element
    std::vector<double> elementStress_;  // 48 per element, the last stress it carried
    std::vector<double> gaussVolume_;    // 8 per element
    std::vector<double> stiffness_;      // 576 per element, elastic path only
    // Per element, when `SolveParams::restFormsCache` leaves it on. 12 kB each,
    // which is the whole of the trade -- see `SolveParams::restFormsCache`.
    std::vector<solidshell::RestForms> forms_;
    std::vector<solidshell::ElementPlasticState> plastic_;
    std::vector<double> dissipation_;    // per element, this step
    // The stiffeners. `fiberForms_` is the step-invariant half, taken from `rest_`
    // after any `Preload` has moved it.
    constraint::RestFibers fiberForms_;
    std::vector<constraint::FiberState> fiber_;
    double fiberEnergy_ = 0;

    std::vector<std::uint32_t> adjacencyOffset_, adjacencyEntry_;  // node -> element corner
    std::vector<std::uint32_t> driven_;   // node indices the indenter holds
    std::vector<std::uint8_t>  pinned_;   // per node
    // The driven perimeter. `edgeDof_` is the pinned DOF with a non-zero prescribed
    // value; `edgeFree_` is the velocity each pinned DOF *would* have taken this
    // step, which is what the boundary's reaction -- and therefore its work -- is
    // measured against. Both empty when nothing is driven, which is the test the
    // step loop branches on.
    std::vector<std::uint32_t> edgeDof_;
    std::vector<double> edgeFree_;
    double edgeFraction(double time) const;

    // The striking body, `Drive::Inertial` only. `punchMass_` stays at zero on the
    // prescribed drive and **is the branch the step loop tests**, so that path
    // executes the same arithmetic in the same order it always did -- a negative
    // control the tests assert bit-for-bit rather than to a tolerance.
    // `drivenMass_` is the plating the punch grips, which moves with it and is
    // therefore part of the assembly's inertia, not part of what resists.
    double punchMass_ = 0;    // kg
    double drivenMass_ = 0;   // kg
    double punchSpeed_ = 0;   // m/s
    bool done_ = false;
};

// --- The whole chain -----------------------------------------------------------

// Mesh a zone round `impact`, drive an indenter into it, and report what tore --
// the FEM answer to the question `sim::impactDamage` answers with a membrane, and
// in the same shape, so `breachesFromFailedPanels(ship, structure, damage.torn())`
// consumes it unchanged.
struct ZoneDamage {
    Patch patch;
    SolveResult result;
    const std::vector<int>& torn() const { return result.tornPanels; }
};

ZoneDamage indent(const StructuralMesh& structure, const Vec3& impact,
                  const MeshParams& mesh, const SolveParams& solve,
                  const plasticity::Material& material = plasticity::shipSteel());

}  // namespace sim::zone
