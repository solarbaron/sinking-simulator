// SPDX-License-Identifier: MIT
//
// Scantlings -> structural mesh: the geometry foundation of Phase 3.
//
// A ship is not specified as a mesh. It is specified as *scantlings* -- shell
// plating thickness by strake and region, frame spacing and section, longitudinal
// stiffener spacing and section, deck and bulkhead plating, girders. That is the
// description a classification society approves and a yard cuts steel from, and
// it is the description this file turns into elements.
//
// Two decisions are made here that everything downstream inherits.
//
// --- 1. Stiffeners are discrete line elements, not smeared into the plate ------
//
// The cheap alternative is a *smeared* (homogenised) plate: replace the plating
// plus its stiffeners with a plate of equivalent thickness t + A/s, where A is
// the stiffener area and s its spacing. It is a standard technique, it is a third
// of the elements, and it is exactly right for axial stiffness.
//
// It is rejected because of what it cannot represent:
//
//   * **Stiffened panel collapse.** `docs/02-simulation.md` §3 says plainly that
//     this is the failure mode that breaks ships' backs and it must not be
//     scripted. Collapse is the plate buckling in half-waves *between* the
//     stiffeners while the stiffeners themselves stay straight, then the
//     stiffener tripping. A homogenised panel has no between; it has one
//     stiffness everywhere and buckles as a slab, at the wrong load, in the wrong
//     mode. There is no tuning of t_eq that fixes this, because the missing thing
//     is a length scale, not a value.
//   * **Bending stiffness of the panel itself.** A stiffener's contribution to
//     the panel's own second moment is dominated by the Steiner term A*d^2 about
//     the combined neutral axis, and d is the stiffener's height. Smearing puts
//     that area *in the plate*, where d = 0. Measured, for a 200x10 flat bar on
//     12 mm plating at 700 mm spacing: the panel's second moment falls from
//     2.49e-5 m^4 to 1.91e-7, **a factor of 130**, at identical area and
//     identical axial stiffness -- which is why an axial check does not notice.
//   * **A stiffener failing separately from its plate.** Cracks start at
//     stiffener ends and bracket toes. A smeared panel has no stiffener end.
//
// What the discrete choice costs: roughly 3x the element count for a stiffened
// shell, and every stiffener is an eccentric beam -- its neutral axis is offset
// from the plate's mid-surface, so the coupling to the shell has to carry that
// offset or the section modulus is wrong by the Steiner term above. Each member
// therefore records the direction its web rises from the plate.
//
// The smeared thickness is *not* thrown away -- `smearedThickness()` computes it,
// because the Tier-0 girder beam and any coarse global model legitimately want
// it. It is a derived quantity here rather than the representation.
//
// --- 2. The structural mesh is independent of the hydrodynamic hull mesh ------
//
// It would be convenient to share vertices with the `TriMesh` the hydrostatics
// integrate. They are not shared, because the two meshes are refined for
// different reasons and neither refinement is negotiable:
//
//   * The hull mesh's resolution is chosen by *volume* convergence.
//     `hullform.hpp` records the measurement: the block coefficient is dominated
//     by waterline count, so the default spends triangles on waterlines (21) and
//     is content with 41 stations -- a 3 m station spacing on a 120 m ship. Its
//     cost is paid every tick, by the buoyancy integral, and `docs/02` §2 records
//     that surface queries against it were once the entire tick.
//   * The structural mesh's resolution is set by the *frame spacing*, which is a
//     property of the ship, not of a solver. 2.4 m web frames on that same ship
//     is 51 transverse stations, and ordinary frames at 600 mm would be 201.
//
// Sharing forces one of two bad outcomes: refine the hull mesh to frame spacing
// and pay for it in the hydrostatic inner loop forever, or snap the frames to
// whatever stations the hydrostatics chose -- which makes frame spacing a
// consequence of a tessellation study instead of a statement about the ship.
//
// So the hull mesh is used as a *reference surface only*. The structure samples
// it by ray casting, which means the same scantlings work against a hull from
// `makeHullFromParticulars`, from a `.ship` offset table, or from an importer,
// with no shared topology and no assumption that the hull is a station grid.
//
// What that rules out: the two meshes cannot deform together for free. When Phase
// 3's FEM moves the structure, the hydrodynamic hull has to be updated by
// projection from the structural mesh, and that mapping is an extra step with its
// own error -- against the alternative, where a shared vertex moves once. It also
// means the two can drift: editing a hull silently invalidates a structural mesh
// built from it, so the mesh is regenerated rather than cached across an edit.
//
// --- Conventions --------------------------------------------------------------
//
// **Panel corners lie on the plate's mid-surface**, which is taken to be the hull
// mesh surface itself. Shipbuilding moulded lines put the surface at the outer
// face of the plating instead; the difference is t/2, six millimetres on a twenty
// metre beam, and the mid-surface convention is what a shell element wants. It is
// applied to decks and bulkheads too, so there is one rule rather than three.
//
// **Girth** runs from the centreline at the keel, outboard across the flat of
// bottom, round the bilge and up the side to the deck edge -- the way a shell
// expansion is drawn. Regions and stiffener positions are given as *fractions* of
// it, so that a longitudinal is continuous along the ship even though the girth
// it lives on shrinks towards the ends.
//
// SI throughout. Body frame per CLAUDE.md: +x forward, +y to port, +z up, origin
// at midship on the baseline, so a z coordinate is a height above the baseline.
//
// Scantlings are not yet authored in `shipfile.cpp`; they are built in C++, the
// way hulls were before `ships/ferry.ship` existed. `ferryScantlings()` is the
// reference arrangement, keyed to the same 120 m ro-pax the rest of the engine is
// validated against.
#pragma once

#include "../core/geometry.hpp"

#include <string>
#include <vector>

namespace sim {

// --- Materials ---------------------------------------------------------------

// The subset of `docs/02-simulation.md` §3's material database that geometry
// needs. Strength and thermal properties arrive with the constitutive model;
// putting placeholders here would be a plausible wrong number, which is the
// failure mode this repo keeps finding.
// **Every value here is at 20 C**, and `thermal::atTemperature` is what turns one
// into the material a hot member has: `youngsModulus` scaled by EN 1993-1-2 §3.2's
// `k_E,theta` and `yieldStrength` by `k_y,theta`. It returns a value rather than
// storing a temperature, so a mesh with a temperature field is a mesh with one
// entry per distinct temperature in `StructuralMesh::materials` and every consumer
// -- `buckling`, `collapse`, `indentation`, `solid_shell` -- reaches it through the
// material index it already uses. Nothing here needs to know.
struct StructuralMaterial {
    std::string name = "steel_ah36";
    double density = 7850.0;         // kg/m^3
    double youngsModulus = 206.0e9;  // Pa
    double poissonRatio = 0.30;
    double yieldStrength = 355.0e6;  // Pa

    // Thermal, at 20 C, from EN 1993-1-2:2005 §3.4.1.3 and §3.4.1.2 evaluated
    // there: `54 - 3.33e-2 * 20` and the cubic at 20 C. They arrive here with
    // `thermal.hpp`, which is the constitutive model this header was waiting for.
    //
    // **These are the room-temperature values and a fire does not stay there.**
    // Conductivity falls 36% by 600 C and specific heat has an eleven-fold spike
    // at 735 C, so `thermal::Problem::temperatureDependent` takes both from the
    // published curves instead of from here. What these two are for is a solve
    // over a range narrow enough that constant properties are the right model --
    // and for the closed forms that validate the operator, every one of which
    // assumes a constant diffusivity.
    double conductivity = 53.334;   // W/(m K)
    double specificHeat = 439.80176;  // J/(kg K)
};

StructuralMaterial mildSteel();   // grade A, 235 MPa
StructuralMaterial ah36Steel();   // higher-tensile, 355 MPa

// --- Stiffener profiles ------------------------------------------------------

// Bulb flats -- by far the most common longitudinal on a real ship -- are
// deliberately absent. IACS treats them by converting to an equivalent angle
// with a published dimensional transformation, and a transformation reproduced
// from memory would be a plausible wrong number in the one place the whole
// section modulus depends on it. Add it with the rule text in hand.
enum class ProfileKind { FlatBar, Angle, Tee };

// Dimensions of a rolled or built-up stiffener. The web rises from the plate;
// the flange, where there is one, sits across its far end.
struct StiffenerProfile {
    ProfileKind kind = ProfileKind::FlatBar;
    double webHeight = 0.200;        // m
    double webThickness = 0.010;     // m
    double flangeWidth = 0.0;        // m, ignored for a flat bar
    double flangeThickness = 0.0;    // m
};

StiffenerProfile flatBar(double webHeight, double webThickness);
StiffenerProfile tee(double webHeight, double webThickness, double flangeWidth,
                     double flangeThickness);
StiffenerProfile angle(double webHeight, double webThickness, double flangeWidth,
                       double flangeThickness);

// Properties of the profile alone, about its own centroid, measured from the
// *attachment line* -- the face of the plate the web is welded to. Independent of
// the plate, so a profile can be tabulated once.
//
// `secondMomentWeak` is taken about the geometric axis parallel to the web, not
// about a principal axis. For a flat bar and a tee those are the same thing; an
// **angle is unsymmetrical** and also has a product of inertia, which is not
// carried here. That matters less than it sounds for the hull girder: a stiffener
// whose web is vertical or horizontal contributes through one axis only and the
// product term drops out identically, so the omission is confined to stiffeners
// sitting on the turn of the bilge. It would matter for lateral-torsional
// buckling of an angle, which is a Phase 3 solver's problem and not this file's.
struct ProfileSection {
    double area = 0;              // m^2
    double centroid = 0;          // m from the attachment line to the centroid
    double height = 0;            // m from the attachment line to the far fibre
    double secondMoment = 0;      // m^4, bending in the plane of the web
    double secondMomentWeak = 0;  // m^4, bending about the web's own plane
    double torsionConstant = 0;   // m^4, St Venant J of the open section
};
ProfileSection profileSection(const StiffenerProfile& profile);

// The stiffener working with its attached plating, which is the combination that
// actually carries load. Distances are measured from the **plate mid-surface**,
// positive towards the stiffener.
//
// `plateWidth` is the effective breadth of plating, which is the caller's
// business: it is the stiffener spacing for a short panel and rather less for a
// long one, because shear lag stops the far edge of the plate from working.
// Nothing here guesses it.
struct StiffenedSection {
    double area = 0;              // m^2
    double neutralAxis = 0;       // m from the plate mid-surface
    double secondMoment = 0;      // m^4 about the combined neutral axis
    double modulusPlate = 0;      // m^3, to the outer fibre of the plating
    double modulusStiffener = 0;  // m^3, to the far fibre of the stiffener
    double height = 0;            // m, plate mid-surface to the far fibre
};
StiffenedSection stiffenedSection(const StiffenerProfile& profile, double plateThickness,
                                  double plateWidth);

// Plating thickness that has the same axial stiffness and mass per unit width as
// the stiffened panel: t + A/s. This is what a homogenised model uses, and it is
// exactly what the header rejects as a *representation* -- but the Tier-0 girder
// beam wants it, so it is computed rather than hidden.
double smearedThickness(double plateThickness, const StiffenerProfile& profile, double spacing);

// --- The scantling description -----------------------------------------------

// A band of shell plating: a range of the length and a range of the section
// girth, with its own thickness and its own longitudinals. Bottom, bilge, side
// and sheer strake are four of these; a thickened region round a breach or a
// reduced one towards the ends are more.
//
// Girth fractions run 0 at the centreline keel to 1 at the deck edge. Regions
// must cover [0, 1] at every station, and `validateScantlings()` reports any
// (girth, length) cell they leave bare. Overlapping is legal and is how a
// thickened patch is written: the **last** region declared that covers a point
// wins, so a patch goes on top of the strake it sits in rather than having to be
// cut out of it.
struct ShellRegion {
    std::string name;
    double xFrom = -1e9, xTo = 1e9;        // m, body frame
    double girthFrom = 0.0, girthTo = 1.0; // fraction of the section girth
    double thickness = 0.012;              // m
    int material = 0;                      // index into Scantlings::materials

    StiffenerProfile longitudinal;         // the stiffeners in this band
    // Nominal girth spacing, measured amidships. Zero takes
    // Scantlings::longitudinalSpacing. The delivered spacing is this rounded to a
    // whole number of bands across the region, exactly as a yard would.
    double longitudinalSpacing = 0.0;
    bool stiffened = true;                 // false leaves the band unstiffened
};

// A horizontal deck at a fixed height, spanning a range of the length and
// clipped outboard to the hull. Longitudinals run at fixed |y| from the
// centreline out -- *not* at constant girth fraction, because a deck longitudinal
// is straight and terminates where the deck narrows past it.
struct Deck {
    std::string name;
    double z = 0;                     // m above the baseline
    double xFrom = -1e9, xTo = 1e9;   // m
    double thickness = 0.010;         // m
    int material = 0;

    StiffenerProfile longitudinal;
    double longitudinalSpacing = 0.0; // m; zero takes Scantlings::longitudinalSpacing
    bool stiffened = true;

    StiffenerProfile beam;            // transverse deck beams, one per frame
    bool beamed = true;
};

// A watertight or structural bulkhead. Transverse ones sit at a station and are
// stiffened vertically; longitudinal ones sit at a fixed |y| -- both sides -- and
// are stiffened longitudinally, which is the usual arrangement and is what
// decides whether the bulkhead is effective in the hull girder.
struct Bulkhead {
    std::string name;
    bool transverse = true;
    double position = 0;              // m: x for a transverse bulkhead, |y| for a longitudinal
    double zFrom = 0, zTo = 0;        // m above the baseline
    // Longitudinal bulkheads only: a longitudinal bulkhead runs between two
    // transverse ones, which is a decision rather than something to derive from
    // the hull. Ignored when `transverse`, where the extent is the hull itself.
    double xFrom = -1e9, xTo = 1e9;   // m
    double thickness = 0.009;         // m
    int material = 0;

    StiffenerProfile stiffener;
    double stiffenerSpacing = 0.0;    // m; zero takes Scantlings::longitudinalSpacing
    bool stiffened = true;
};

// A deep longitudinal member on a named line: the centre girder on the keel, side
// girders in the double bottom, a deck girder under a vehicle deck. Positioned
// explicitly rather than derived, because a girder is a decision rather than a
// spacing.
struct Girder {
    std::string name;
    double y = 0, z = 0;              // m, the attachment line
    double xFrom = -1e9, xTo = 1e9;   // m
    bool bothSides = false;           // mirror to -y as well
    StiffenerProfile profile;
    Vec3 rise{0, 0, 1};               // direction the web rises from the plating
    double attachedPlateThickness = 0.012;  // m, the plating it works with
    int material = 0;
};

struct Scantlings {
    // Nominal transverse frame spacing. Delivered as the hull length divided by a
    // whole number of bays, so the ends land on frames instead of on a short bay;
    // `problems` reports the difference if it is more than a percent.
    double frameSpacing = 2.40;       // m
    StiffenerProfile frameProfile;    // transverse frames on the shell
    int frameMaterial = 0;
    bool framed = true;

    double longitudinalSpacing = 0.70;  // m, default for anything that omits one

    // How finely the hull section is sampled when its girth is measured. The
    // panels chord across this, so it bounds how well curvature is followed.
    int girthSamples = 192;

    std::vector<ShellRegion> shell;
    std::vector<Deck> decks;
    std::vector<Bulkhead> bulkheads;
    std::vector<Girder> girders;
    std::vector<StructuralMaterial> materials{ah36Steel()};
};

// --- The structural mesh ------------------------------------------------------

enum class PanelRole { Shell, Deck, Bulkhead };
enum class MemberRole { Longitudinal, Frame, DeckLongitudinal, DeckBeam, BulkheadStiffener, Girder };

// A quadrilateral plating panel bounded by two frames and two longitudinals.
// Corners are in order round the panel and lie on the plate's mid-surface.
struct PlatePanel {
    Vec3 corner[4]{};
    double thickness = 0;
    int material = 0;
    PanelRole role = PanelRole::Shell;
    int source = -1;      // index into the ShellRegion / Deck / Bulkhead it came from

    double area() const;      // m^2, as two triangles, so a warped quad is still exact
    Vec3 centroid() const;    // area-weighted
    Vec3 normal() const;      // unit, from the quad's diagonals
};

// A stiffener as a line element with its section. `rise` is the direction the web
// grows from the plating: the eccentricity that a beam offset from a shell has to
// carry, and the thing smearing throws away.
struct StructuralMember {
    Vec3 a{}, b{};
    Vec3 rise{0, 0, 1};
    StiffenerProfile profile;
    double attachedPlateThickness = 0;
    int material = 0;
    MemberRole role = MemberRole::Longitudinal;

    double length() const;
};

struct StructuralMesh {
    std::vector<PlatePanel> panels;
    std::vector<StructuralMember> members;
    std::vector<StructuralMaterial> materials;

    // As built rather than as requested -- see Scantlings::frameSpacing.
    double frameSpacing = 0;
    std::vector<double> frameStations;   // m, ascending

    double plateArea() const;    // m^2
    double plateMass() const;    // kg
    double memberLength() const; // m
    double memberMass() const;   // kg
    double steelMass() const { return plateMass() + memberMass(); }

    std::size_t panelCount(PanelRole role) const;
    std::size_t memberCount(MemberRole role) const;
};

// Generate the structure. `hull` is a reference surface only: it is sampled, not
// shared. Anything the description asks for that the hull cannot provide -- a
// deck above the sheer, a bulkhead outside the hull, a region that leaves a gap
// in the girth -- is reported in `problems` rather than quietly dropped, along
// with everything `validateScantlings()` would have said about the description on
// its own, so one call is enough to know whether the answer is trustworthy.
//
// It always returns *something*: an unbuildable description yields a thin or
// empty mesh and a full account of why, in the same spirit as Ship::validate().
StructuralMesh makeStructuralMesh(const TriMesh& hull, const Scantlings& scantlings,
                                  std::vector<std::string>* problems = nullptr);

// Every way this description is not a ship. Advisory, like Ship::validate().
std::vector<std::string> validateScantlings(const Scantlings& scantlings);

// --- Hull girder ---------------------------------------------------------------

// The section properties of the hull girder at a transverse cut: what Tier 0
// needs, and the number a midship section drawing exists to produce.
//
// Membership is decided **geometrically, not by label**: a panel or member counts
// only if it has a non-zero extent along x and straddles the cut. Frames, deck
// beams and transverse bulkheads have zero x extent and so contribute nothing,
// which is correct -- an athwartships member carries no longitudinal stress --
// and cannot be got wrong by mis-tagging an element.
struct HullGirderSection {
    double x = 0;
    double area = 0;            // m^2 of longitudinally effective material
    double neutralAxis = 0;     // m above the baseline
    double secondMoment = 0;    // m^4 about the horizontal axis through the neutral axis
    double zDeck = 0, zKeel = 0;   // m, the extreme fibres actually present
    double modulusDeck = 0;     // m^3
    double modulusKeel = 0;     // m^3
};
// One longitudinally-effective piece of a transverse cut: a strip of plating, or
// a stiffener. `hullGirderSection()` is the sum of these, and progressive
// collapse (`collapse.hpp`) needs them individually, because the whole method is
// that each one goes non-linear at its own strain.
struct SectionElement {
    double area = 0;             // m^2
    double height = 0;           // m above the baseline, of its own centroid
    double ownSecondMoment = 0;  // m^4 about its own horizontal centroidal axis
    double zLo = 0, zHi = 0;     // m, its extreme fibres
    double thickness = 0;        // m, plating thickness or stiffener web thickness
    double width = 0;            // m, girth spanned; zero for a stiffener
    // A stiffener's own profile and the plating it is welded to, kept because the
    // column-buckling check needs both and `thickness` above is neither: for a
    // stiffener it is the *web* thickness. Empty and zero for a plate element.
    StiffenerProfile profile{};
    double attachedPlateThickness = 0;  // m
    int material = 0;
    bool stiffener = false;
};

// Decompose a transverse cut into its elements. The extent is half-open,
// [xLo, xHi), so a cut landing on a frame seam is served by the bay forward of it
// and not by both -- counting both doubles area and second moment while leaving
// the neutral axis, a ratio, looking perfectly correct.
//
// The seam is fuzzy by a nanometre and `x` snaps onto it, so a caller need not
// reproduce the arithmetic the stations were laid out with: `120*33/50` lands one
// unit in the last place above the 19.2 a drawing carries, and both name the same
// station. The section is continuous across a seam, and it is asked for at one.
std::vector<SectionElement> sectionElements(const StructuralMesh& mesh, double x);

HullGirderSection hullGirderSection(const StructuralMesh& mesh, double x);

// The midship section modulus a classification society will not go below, from
// the IACS unified requirement: Z = C L^2 B (Cb + 0.7) x 1e-6 m^3, with
// C = 10.75 - ((300 - L)/100)^1.5 over 90 m <= L <= 300 m. It is the one
// externally published number that says whether a set of scantlings is a ship or
// a guess, which is why it is here rather than in a test.
double ruleMinimumSectionModulus(double lengthPp, double beam, double blockCoefficient);

// --- Reference ship ------------------------------------------------------------

// The 120 m ro-pax the rest of the engine is validated against: `ships/ferry.ship`
// and `game/prototype/ferry.cpp`. Longitudinally framed with 2.4 m web frames,
// bulkheads on the same stations as that file's compartment boundaries, and a
// double bottom to 1.8 m. Scantlings are typical for the size rather than taken
// from any real ship's drawing -- there is no more a published midship section for
// a generic ferry than there is a published offset table.
Scantlings ferryScantlings();

}  // namespace sim
