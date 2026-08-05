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
// **7.3 µs** per step against 273 ns elastic, and the stable step is `t / c_p` --
// **thickness-governed, flat in the in-plane element size**. So for 12 mm plating
// the step is 1.8 µs, 5.5 x 10^5 steps per simulated second, and
//
//     core-seconds per simulated second = 4.0 x elementCount
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
// longitudinals), a fourteen-panel zone is 224 elements at 900 core-seconds per
// simulated second, and driving a punch 0.23 m into her takes 0.04 s of simulated
// time -- **4.5 s of wall time on 23 workers**, against 15.4 s on one; the speedup
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
// --- 3. Stiffeners are not meshed, and the answer is a bracket -----------------
//
// The zone is plating only, and the reason is that **the only element in the
// inventory is the solid-shell hex, and there is no way to attach a web to a plate
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
// What is actually needed is a multi-point constraint tying an eccentric beam to
// the shell, which is the same machinery Tier-1/Tier-2 interface coupling needs
// and does not exist yet.
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
// So `Stiffeners` offers the two bounds the omission leaves, and the honest answer
// is the bracket rather than either end:
//
//   * `Ignored` -- the lower bound. Soft by the factor above; conservative for a
//     damage-stability question and useless for a survivability one.
//   * `RigidSupport` -- the upper bound, and the default. Every plating node a
//     stiffener runs through is pinned, which is what a member far stiffer than the
//     plate it carries does in the limit. Measured on the same zone, the resisting
//     force comes out 23-30 MN against the membrane model's 28 MN on the same
//     0.70 m span, so the plating between supports is behaving as the membrane
//     model says it should.
//
// It is also the setting that makes the two models comparable at all, because
// `indentation.hpp` assumes one bay with its boundaries held -- which is precisely
// a rigid support on every stiffener line.
//
// --- 4. Boundary conditions --------------------------------------------------
//
// The patch is cut out of a ship, so its edge is a lie either way. Free lets it
// translate away under the first load; clamped is too stiff, because the plating
// outside the zone really does pull in. Clamped is chosen, and **the price is paid
// in zone size rather than in a tuned spring stiffness**: the boundary error is a
// Saint-Venant effect that decays away from the edge, so growing the radius makes
// it go away and the convergence is measurable. `tests/test_zone.cpp` measures it.
// An elastic edge with a stiffness standing for the surrounding structure would be
// softer and closer, and it would also be a coefficient no measurement sets --
// which is the argument `solid_shell.hpp` uses against hourglass control.
//
// The consequence to keep in mind: **a zone one bay wide is a clamped bay**, which
// is the same idealisation `indentation.hpp` makes, and a zone five bays wide is
// not. That is deliberate -- it is what makes the two models comparable at radius
// one and different at radius five.
//
// --- 5. What this does not do yet ---------------------------------------------
//
//  1. **No coupling to Tier 1.** What does exist is coupling to Tier 0, in
//     `promotion.{hpp,cpp}`: it decides when a patch deserves a zone, hands the
//     zone the girder's stress as a `Preload` below, and turns what tore back into
//     a section the beam can read. The boundary is still fixed -- the patch's edge
//     is clamped and nothing outside it moves in response to what happens inside --
//     so the coupling is one way per solve.
//
//     The **reduction** the two-way version needs now exists:
//     `reduction.{hpp,cpp}` will take this patch's own mesh, treat its clamped
//     perimeter as an interface (`reduction::nodesPinned`) and hand back a reduced
//     pair whose boundary DOF are kept exactly. What is still missing is the
//     *coupling*: something that drives those interface DOF from the surrounding
//     structure instead of pinning them, and a matching pass that ties two
//     substructures' interfaces together. Until that exists the edge is clamped,
//     and `zone.cpp` is unchanged by the reduction's arrival.
//  2. **The indenter is kinematic and rigid.** Nodes inside a rectangular footprint
//     are driven at a prescribed velocity; there is no contact search, no friction,
//     no release, and the striking body does not crush. A prescribed motion cannot
//     run away, which is what makes it testable; a delivered *energy* would need
//     the striking body's mass and is what `collision.hpp` would eventually supply.
//  3. **Element deletion, not splitting.** A torn element carries no stress and
//     keeps its mass. The hole is therefore the deleted area, which this reports as
//     whole panels because that is what `breachesFromFailedPanels` consumes -- see
//     `SolveParams::tearFraction`.
//  4. **No GPU path**, and no rate dependence in the material, so the resistance is
//     under-predicted by the 10-30% steel gains at collision strain rates.
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
    // 0.35 m of penetration under a 2 m punch, **seven times** too soft against a
    // membrane model spanning the real 0.70 m longitudinal spacing. Worse, the
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

    // Smallest stable explicit step over the patch, computed here rather than in
    // the solver because it is what makes the cost knowable *before* the run: it
    // is 0.5 ms of power iteration per element, which is a promotion cost like the
    // stiffness formation beside it and not a per-step one.
    double criticalTimestep = 0;    // s

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

// A rigid flat punch, rectangular in the patch's own in-plane frame, moving at a
// constant speed along -axis. Nodes on the **outer face only** are driven, so the
// plate is free to thin under the punch -- which is the whole reason the element
// keeps a thickness stretch mode. It grips: a driven node is held to the punch
// rather than released, which is right for a monotonic push and wrong for
// rebound.
// A footprint with a non-positive half-extent is **no punch at all**, which is
// how a patch is asked to carry no load. Leaving the footprint and setting the
// speed to zero is a different thing -- a stationary punch still holds the nodes
// it touches, and holding them is exactly what stops an unloaded patch from
// translating freely.
struct Indenter {
    double halfLength = 1.0;  // m along `Patch::right`
    double halfWidth = 1.0;   // m along `Patch::up`
    double speed = 6.0;       // m/s
    double stopAt = 0.5;      // m of penetration; the run ends here
    // Seconds over which the speed rises linearly from zero. Zero starts the punch
    // at full speed, which is a genuine step in velocity: the plate ahead of it
    // sees rho*c*v, 280 MPa at 6 m/s, and rings for as long as the run lasts. A
    // ramp is a statement about the striking body's own compliance rather than a
    // numerical fudge, and the energy account holds either way.
    double rampTime = 0.0;    // s
};

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

    // False solves the co-rotational *elastic* element -- 273 ns against 7.3 µs,
    // and the path the geometric tests use, where plasticity would only add noise
    // to an identity.
    bool plastic = true;

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
};

struct Sample {
    double time = 0, penetration = 0, force = 0;
    double work = 0, strainEnergy = 0, dissipation = 0, kinetic = 0;
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
    double energyResidual() const {
        return work - ((strainEnergy - initialStrainEnergy) + dissipation + kinetic + dampingLoss);
    }

    int yieldedElements = 0;
    int tornElements = 0;
    double tornArea = 0;             // m^2 of deleted element
    std::vector<int> tornPanels;     // ascending indices into StructuralMesh::panels
    double tornPanelArea = 0;        // m^2 of those panels, as `breach.hpp` will see it

    // Cost, measured. The figures are for printing and for `estimatedCost`, never
    // for asserting on: `test_plasticity.cpp` records what a tight timing assertion
    // costs on a shared machine.
    double wallSeconds = 0;
    double microsecondsPerElementStep = 0;

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
    std::vector<solidshell::ElementPlasticState> plastic_;
    std::vector<double> dissipation_;    // per element, this step

    std::vector<std::uint32_t> adjacencyOffset_, adjacencyEntry_;  // node -> element corner
    std::vector<std::uint32_t> driven_;   // node indices the indenter holds
    std::vector<std::uint8_t>  pinned_;   // per node
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
