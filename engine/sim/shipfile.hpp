// SPDX-License-Identifier: MIT
//
// The ship definition format: a `sim::Ship` as text rather than as C++.
//
// `docs/05-data-modding-validation.md` §1 commits to ships being declarative,
// version-controlled and human-diffable, and §4 to a data mod needing no code at
// all. `game/prototype/ferry.cpp` is the thing that has to stop being necessary:
// it builds one ship, in a function, and a fleet cannot be built that way. This
// header is the first half of the fix -- the reader. `ships/ferry.ship` is the
// same ferry as data, and `tests/test_shipfile.cpp` holds the two side by side.
//
// **The grammar**, whitespace-separated tokens, one directive per line, `#` to
// end of line is a comment. Deliberately the same lexical idiom as
// `engine/gpu/material.cpp`, which is the text format this repo already has; a
// second style would be a second thing to learn for no gain.
//
//     ship_format 1                    must be the first directive in the file
//
//     ship <name>
//         deck_edge_z    <m>           weather deck edge, within the hull's own
//                                      z range -- it is what freeboard is measured to
//         lightship_mass <kg>          exactly one of mass or draft
//         lightship_draft <m>          mass = rho * (hull volume below this plane),
//                                      so an edited hull moves the loading with it
//         lightship_cog  <x> <y> <z>   m, body frame
//         gyradii        <kxx> <kyy> <kzz>   m, about the lightship CoG
//         sea_density    <kg/m3>       optional, default 1025
//         damping        <heave> <roll> <pitch>    optional, fraction of critical
//         added_mass     <surge> <sway> <heave>    optional, x displaced mass
//         added_inertia  <roll> <pitch> <yaw>      optional, x displaced inertia
//
//     hull
//         waterlines <z0> <z1> ...     m, ascending, at least two
//         station <x> <hb0> <hb1> ...  m; one half-breadth per waterline,
//                                      stations in ascending x, at least two
//
//     compartment <name>
//         box <loX> <loY> <loZ> <hiX> <hiY> <hiZ>   m, body frame
//         permeability <fraction>      void fraction, within (0, 1]
//         vented <true|false>          true means it never pressurises
//
//     opening <name>
//         between <a> <b>              compartment names, or `sea`
//         at <x> <y> <z>               m, body frame, the orifice centre
//         area <m2>
//         discharge <cd>               within (0, 1]
//         kind <breach|door|hatch|vent|pipe>
//         open <true|false>
//
//     pump <name>
//         drains <name>                the compartment it empties overboard
//         capacity <m3/s>              at zero head
//         max_head <m>
//         on <true|false>
//
// `drains` rather than `compartment` so that `compartment` means exactly one
// thing wherever it appears -- it opens a block. A word that opens a block in one
// place and sets a field in another makes the compartment declared straight after
// a pump unparseable, which is a grammar that only works until someone reorders
// their file.
//
// SI throughout, as CLAUDE.md requires. Nothing in the format is an angle, so
// the radians rule never comes up; if a raked bulkhead ever arrives it arrives in
// radians.
//
// **Every key is required except where a default is documented above,** and the
// only documented defaults are `sea_density`, `damping`, `added_mass` and
// `added_inertia` -- engine-wide physics placeholders rather than claims about
// this ship. `permeability` and `vented` in particular are not optional: both
// have plausible defaults, and a plausible wrong number is the failure mode this
// repo keeps finding. `kind` is required for the same reason -- a watertight door
// defaulting to "breach" is a lie about the ship, even though nothing in the
// physics reads it yet.
//
// **A compartment is carved, not authored.** `box` is six bulkhead and deck
// planes, and the compartment is the hull interior between them -- `clipToBox`,
// exactly as `ferry.cpp` does it. So a wing tank tapers into the turn of the
// bilge and a forepeak narrows into the stem without the author knowing where
// either is. A box that misses the hull entirely produces nothing and is
// rejected, because that is the prototype's oldest authoring mistake and it is
// silent: an empty compartment simply never floods.
//
// **Unknown keys are rejected, not ignored.** The brief for a mod format usually
// argues the other way -- old loaders meet new files -- but consider what a
// tolerated key costs here. `permeabilty 0.85` on a machinery space is a typo
// that leaves the space at the 0.95 default, and 0.95 is a *plausible* number, so
// nothing downstream looks wrong; the ship just floods 12% more than the author
// wrote. That is precisely the silent failure CLAUDE.md's defect table is made
// of. Forward compatibility is bought instead with `ship_format`: a file
// declaring a version this build does not know is refused **by name**, so an old
// loader meeting a new file says so rather than quietly dropping the field that
// changed the ship's meaning. Data a mod wants to carry that the *simulation*
// does not read belongs in a sibling file -- docs §1 already makes a ship a
// directory, and `art/`, `systems/` and `loading/` are where that goes.
//
// **The load path fails closed.** Everything is parsed into a scratch definition
// and `out` is assigned only once the whole file has parsed, every name has
// resolved and every compartment has come back from the clip with geometry. An
// error names the origin, the line and what was wrong. CLAUDE.md records a
// `World::load` that failed open and left a half-built world; the instrument that
// caught it -- every truncation of a valid file -- is pointed at this loader too,
// in `tests/test_shipfile.cpp`.
#pragma once

#include "ship.hpp"

#include <string>
#include <string_view>

namespace sim {

// The format version this build understands. Bumped only when a change would
// make an older loader misread a newer file; additive keys that an older loader
// would reject anyway still need it, which is the point of refusing rather than
// ignoring.
inline constexpr int kShipFormatVersion = 1;

struct ShipDefinition {
    std::string name;  // the identity token from the `ship` line
    Ship ship;
};

// Read `path`. On success `out` holds the definition; on any failure `out` is
// untouched and `error` carries "<path>:<line>: <reason>".
//
// The loader checks what it can decide locally: a closed, positively-oriented
// hull, ranges on every physical quantity, and that no compartment clipped away
// to nothing. `Ship::validate()` is still the deeper check -- overlapping
// subdivision and dangling indices are its job -- and callers should still run it
// after `initialise()`.
bool loadShipFile(const std::string& path, ShipDefinition& out, std::string& error);

// Same, from memory. `origin` only labels the diagnostics.
bool parseShipFile(std::string_view text, const std::string& origin, ShipDefinition& out,
                   std::string& error);

}  // namespace sim
