// SPDX-License-Identifier: MIT
//
// Failed structure -> holes in the flooding network.
//
// `docs/02-simulation.md` §1 lists "progressive flooding through the structure --
// once elements tear, new orifices appear automatically" as the last thing the
// flooding model is missing, and `docs/06-roadmap.md` makes it the Phase 3
// milestone: the hull tears where the stress says it should, "and the resulting
// hole floods at a rate the hole's own area determines". This file is that last
// clause. Give it a structural mesh and the panels that failed; it gives back
// `Opening`s that drop straight into `Ship::openings` and are, to the flooding
// solver, indistinguishable from ones a human authored.
//
// **It contains no fracture model and does not pretend to.** Which panels failed
// is the caller's answer -- eventually the FEM's, today a test's or a scripted
// scenario's. Everything here is about the consequence.
//
// --- 1. Connectivity is read off the geometry, never off the label -------------
//
// The one question that must not be got wrong is *which two spaces the hole
// joins*: a shell panel opens a compartment to the sea, a bulkhead panel opens
// one compartment to another, and confusing them is the difference between a
// breach and a door. `PanelRole` looks like the answer and is not:
//
//   * A `PanelRole::Bulkhead` panel is not necessarily a compartment boundary.
//     The reference ferry's wing bulkhead runs at |y| = 6 m, *inside* the engine
//     rooms and holds, because the hull is too narrow at the tank top to carry it
//     on the compartment plan's y = 8 m. Tearing it opens a space to itself,
//     which is to say it opens nothing, and `Ship::validate()` would call the
//     resulting orifice a fault.
//   * A `PanelRole::Deck` panel at the weather deck has the sky above it. A hole
//     there is open to the sea in every sense the flooding network cares about.
//   * A `PanelRole::Shell` panel is not necessarily against a compartment either.
//     Amidships on the ferry the engine-room boxes stop at |y| = 8 m while the
//     shell is out at 8.6 m, so most of the side plating over the engine room
//     faces an unmodelled void, and a tear there floods nothing at all.
//
// So each failed panel is probed a short distance either side of its own plane,
// and whatever space is found there is what it connects. `hullGirderSection()`
// takes its membership from geometry for the same reason and states it in the
// same words: a decision that can be got wrong by mis-tagging an element will
// eventually be got wrong by mis-tagging an element.
//
// The probe distinguishes three answers -- a compartment, `kSea` outside the
// hull, and `kEnclosedVoid` inside the hull but inside no compartment. The third
// is the interesting one: it is not an error, it is a hole into a space this ship
// definition does not model, and the honest response is to open nothing and say
// so rather than to invent either a flooded compartment or a leak to the sea.
//
// A fourth case has no honest answer at all: **overlapping compartments**. The
// reference ferry's forward wing tanks lie entirely inside her forward holds --
// `fwd_hold_p` spans y = 0..20 m while `wing_tank_fwd_p` spans 8..20 m over the
// same length and height -- so 217 m3 of her is described twice, and a point out
// there is genuinely in two compartments. `Ship::validate()` does not catch it,
// because the subdivision still totals only 89% of the hull and a shortfall is
// normal. The lowest-numbered compartment wins, which is arbitrary, so the
// overlap is reported: an arbitrary answer that says so is a different thing from
// an arbitrary answer that does not.
//
// --- 2. A torn plate is one hole, not forty ------------------------------------
//
// Failed panels that share an *edge* and join the same pair of spaces merge into
// a single opening whose area is the sum of theirs and whose position is the
// area-weighted centroid of the region. Two things follow from the qualifiers.
//
// **Same pair of spaces.** A shell panel and the bulkhead panel it lands against
// share an edge and are not the same hole: one joins the sea to a compartment,
// the other joins two compartments, and an `Opening` has exactly two ends.
// Merging across the pair would silently pick one of them.
//
// **Share an edge, not a corner.** Two panels touching only at a corner stay two
// openings. The pinch between them has zero width, so no water crosses it; and
// the merged centroid would sit at the corner, which is the one point in the
// region where there is no hole. It costs nothing in total flow -- `Cd·A·√(2Δp/ρ)`
// is linear in area, so N holes at one head pass exactly what one hole of the
// same total area passes -- and it keeps the head right when the pieces sit at
// different depths, where splitting is the *better* quadrature of `∫√h dA` and
// merging is the approximation.
//
// --- 3. Discharge coefficient ---------------------------------------------------
//
// `kTornPlateDischarge` = 0.60. A clean sharp-edged orifice runs 0.61-0.62 (the
// ferry's authored breach uses 0.62) and a rounded or ducted entry 0.8-0.95, so
// the question is only which side of the sharp-edged value a tear sits. It sits
// at or below it: fractured plating petals, and an edge that protrudes into the
// flow is re-entrant, whose contraction coefficient is 0.5 in the Borda limit.
// Damage-stability practice takes 0.6 for a damage opening, and flooding model
// tests report 0.6-0.75 for openings that are clean rectangles cut in a
// bulkhead -- which a tear is not. 0.60 is therefore the top of the honest range
// rather than the middle of it, chosen so that the failure is towards flooding
// too slowly rather than too fast.
//
// It is a constant because nothing upstream can yet refine it. What would refine
// it is the tear's aspect ratio and which way the plating folded, and a set of
// failed panels records neither.
//
// --- 4. Which law the hole floods under, from the way the region faces ---------
//
// `Opening::kind` is not a label. It is the switch on the flooding model:
// `Ship::horizontalSidesOf` admits `OpeningKind::Hatch` and nothing else to the
// counter-current branch, and `ship.hpp` opens by arguing why that branch has to
// exist -- through a hole in a deck with the sea standing on it, water falls
// while air rises through the *same* hole, driven by the density difference and
// by no net pressure difference at all, so "a net-only model sees a still one".
// Every torn region here used to come out `Breach`, which is also the zero
// enumerator, so "nobody chose" and "chose `Breach`" were the same bytes and a
// hole torn through a deck was handed the single-dp vertical orifice law. That
// failure is not a refusal and not a crash; it is a plausible small net flow
// where a tonne a second crosses in two directions, which is the shape of thing
// this whole file exists to get right.
//
// **The crossover needs no constant.** A piece of plating whose normal stands at
// theta from the vertical presents `A cos(theta)` to anything crossing it
// downwards and `A sin(theta)` to anything crossing it sideways, and those are
// equal at 45 degrees. So the region's own two projections are summed over its
// panels and compared:
//
//     sum A_i |n_z,i|  >  sum A_i |n_xy,i|   ->  Hatch, else Breach
//
// That is a crossover rather than a tuned threshold: it is the one angle at
// which neither law is the better description, so it is the only place the
// answer can change without someone having chosen a number.
//
// **A region exactly on it comes out `Breach`**, the comparison being strict.
// The cost of landing on the wrong side there is bounded and it is not
// symmetric. `Opening` carries no normal -- `horizontalSidesOf` takes a hatch's
// to be the body's own +z, because a hatch is in a deck -- so a region called
// `Hatch` at theta exchanges over `A` where only `A cos(theta)` of it faces up,
// 1.41x too much at the crossover; a region called `Breach` exchanges over
// nothing at all. Nothing between the two is available: `kind` is the only
// switch there is, and there is no third law for it to select.
//
// **`PanelRole` is not the criterion, because it does not carry orientation.**
// `PanelRole::Shell` covers the flat of bottom, which is as horizontal as any
// deck -- and is a deck with the sea standing on it the moment she rolls past
// ninety degrees, which is exactly the case `ship.hpp` says the exchange has to
// reverse for -- while the same role covers the side plating, which is a wall.
// `Deck` and `Bulkhead` are generated at constant z and constant x or |y| today,
// so for those two the label and the geometry happen to agree; `Shell` is the
// counterexample and it is most of the ship's watertight envelope. This is §1's
// argument arriving at the same answer from the other end: a decision that can be
// got wrong by mis-tagging an element will eventually be got wrong by mis-tagging
// an element.
//
// A hatch produced this way keeps the `breach_` name prefix. The name records
// where the hole came from; the kind records how it floods, and they are
// different questions. `fire.cpp`'s `ventShapeFor` reads the same enumerator and
// gives a `Hatch` no height for its doorway integral, which is the right answer
// for a hole in a deck for the same reason.
//
// **And it puts a price on §2's merge rule that a breach does not have.** The
// claim there -- that splitting a torn region into N holes costs nothing in total
// flow -- is a statement about `Cd·A·√(2Δp/ρ)`, which is linear in area. The
// horizontal law is not: its driving head is `(ρ_up − ρ_lo) g D/2` with
// `D = √(4A/π)`, so the exchange goes as `A^{5/4}` -- the `D^{5/2}` that
// `ship.hpp` cites Epstein and Cooper for -- and N equal pieces of one deck hole
// pass `N^{-1/4}` of what the whole hole passes, 16% less at N = 2. Merging on a
// shared edge is therefore load-bearing for a hatch where it was only tidy for a
// breach, and the corner-touching split §2 keeps is now a modelling choice with a
// cost rather than a free one. It is kept: the pinch between two panels meeting
// at a corner still has zero width and still passes nothing, and the alternative
// is to merge regions that are not one hole.
//
// SI throughout, body frame per CLAUDE.md.
#pragma once

#include "scantlings.hpp"
#include "ship.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace sim {

// Discharge coefficient of a torn plate edge. See §3 above.
inline constexpr double kTornPlateDischarge = 0.60;

// Sentinel meaning "inside the hull, but inside no compartment" -- structure,
// voids, and anything the subdivision simply does not describe. Distinct from
// `kSea`, which is outside the hull altogether, because the two lead to opposite
// answers and only one of them floods anything.
inline constexpr int kEnclosedVoid = -2;

// The signed number of times `mesh` wraps around `point`, by summed solid angle:
// 1 inside a closed outward-wound mesh, 0 outside.
//
// Ray casting is the obvious alternative and is the one that fails quietly. A ray
// through an edge or a vertex is counted twice or not at all, and compartment
// meshes are carved from the hull on axis-aligned bulkhead and deck planes, so
// their edges lie along exactly the axes a ray would be written along. The solid
// angle sum has no preferred direction and no branch to get wrong.
//
// **On the surface the answer is decided by the sign of a floating-point zero**,
// and is not to be relied on. A triangle the point is coplanar with makes the
// triple product cancel to ±0, and `atan2(+0, negative)` is +π where
// `atan2(-0, negative)` is −π, so that face contributes ±2π and the total lands
// on 1 or on 0 according to which way the cancellation rounded. Measured on a
// box: a point on a face whose outward normal is +x, +y or +z reads 1, and on the
// opposite three reads 0. The consequence that matters is that two compartments
// which merely *abut* can both claim a point on the bulkhead between them, which
// is why the overlap check below asks a second time a millimetre away. Keep the
// query off surfaces.
double meshWindingNumber(const TriMesh& mesh, const Vec3& point);

// Which space a body-frame point lies in: a compartment index, `kSea`, or
// `kEnclosedVoid`. O(triangles in the ship); for repeated queries, build the
// breach set in one call rather than looping over this.
int spaceAt(const Ship& ship, const Vec3& bodyPoint);

struct BreachParams {
    // The first distance either side of a panel at which to look for the spaces
    // it separates, and the smallest. A flat panel *chords* across a curved
    // shell, so its centroid is not on the surface it stands for: measured at
    // 0.15 m on the reference ferry, where a girth band spans the crease at the
    // turn of the bilge and the panel cuts the corner. A single fixed probe
    // cannot serve, because one small enough to stay inside a double bottom is
    // smaller than that and reads the same space on both sides -- which drops the
    // breach silently, the worst available failure. So the probe marches: it
    // doubles until the two sides disagree, and stops at the *first*
    // disagreement, which is the nearest boundary and therefore the surface this
    // panel represents. It gives up at the distance from the panel's centroid to
    // its own furthest corner, on the grounds that the surface a panel stands for
    // passes within the panel's own footprint and anything further away is some
    // other panel's business.
    double probeDistance = 0.05;   // m

    double dischargeCoeff = kTornPlateDischarge;

    // Openings are named `<prefix>_<space a>_<space b>_<n>`.
    std::string namePrefix = "breach";

    // Corners closer together than this are the same corner, which is what makes
    // two panels adjacent. Matches the weld tolerance the CSG uses.
    double weldEpsilon = 1e-6;     // m
};

// One hole, and the failed panels that made it.
struct Breach {
    Opening opening;            // ready for Ship::openings
    std::vector<int> panels;    // ascending indices into StructuralMesh::panels
};

struct BreachSet {
    std::vector<Breach> breaches;

    // Every reason a failed panel produced no opening, and every complaint about
    // the input. Advisory, in the spirit of `Ship::validate()`: a tear into a
    // space the ship does not model is a fact about the ship definition, not a
    // failure of this routine, and swallowing it would leave a hole that floods
    // nothing looking exactly like no hole at all.
    std::vector<std::string> problems;

    std::vector<Opening> openings() const;
    double totalArea() const;    // m^2
};

// The openings a set of failed panels implies.
//
// `failedPanels` are indices into `mesh.panels`, in any order; duplicates and
// out-of-range entries are reported and ignored rather than counted twice. The
// mesh needs no relationship to `ship.hull` beyond sharing its body frame -- it
// is probed against the ship's own compartments, not assumed to match them.
BreachSet breachesFromFailedPanels(const Ship& ship, const StructuralMesh& mesh,
                                   const std::vector<int>& failedPanels,
                                   const BreachParams& params = {});

// Add every breach to the ship's flooding network. Returns how many were added.
std::size_t applyBreaches(Ship& ship, const BreachSet& set);

}  // namespace sim
