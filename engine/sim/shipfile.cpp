// SPDX-License-Identifier: MIT
#include "shipfile.hpp"

#include "../core/geometry.hpp"

#include <charconv>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <system_error>
#include <vector>

namespace sim {
namespace {

// Whitespace-separated tokens of one line, with everything from '#' dropped.
void tokenise(std::string_view line, std::vector<std::string_view>& out) {
    out.clear();
    const std::size_t comment = line.find('#');
    if (comment != std::string_view::npos) line = line.substr(0, comment);
    std::size_t i = 0;
    while (i < line.size()) {
        while (i < line.size() && (line[i] == ' ' || line[i] == '\t' || line[i] == '\r')) ++i;
        if (i >= line.size()) break;
        const std::size_t start = i;
        while (i < line.size() && line[i] != ' ' && line[i] != '\t' && line[i] != '\r') ++i;
        out.push_back(line.substr(start, i - start));
    }
}

// std::from_chars rather than strtod: the decimal point must not depend on the
// host's locale, or the same ship file means different things on two machines.
// The whole token has to be consumed, so "5.5m" is rejected rather than read as
// 5.5 -- a loader that accepts a typo is a loader that fails open.
//
// Non-finite is rejected here rather than by each caller's range check. No
// coordinate, area or mass in a ship is legitimately infinite, from_chars does
// accept "inf" and "nan", and a NaN in the offsets table poisons every
// hydrostatic integral downstream without ever raising anything.
bool parseNumber(std::string_view token, double& out) {
    const char* begin = token.data();
    const char* end = begin + token.size();
    const std::from_chars_result result = std::from_chars(begin, end, out);
    return result.ec == std::errc{} && result.ptr == end && std::isfinite(out);
}

std::string at(const std::string& origin, int line, const std::string& what) {
    std::ostringstream message;
    message << origin << ":" << line << ": " << what;
    return message.str();
}

std::string number(double v) {
    std::ostringstream text;
    text << v;
    return text.str();
}

enum class Block { None, Ship, Hull, Compartment, Opening, Pump };

// Everything lands in these first. Nothing touches the caller's ShipDefinition
// until the whole file has parsed, every reference has resolved and every
// compartment has come back from the clip with geometry in it.
struct PendingShip {
    std::string name;
    int declaredAt = 0;
    double deckEdgeZ = 0;
    double lightshipMass = 0;
    double lightshipDraft = 0;
    Vec3 cog{};
    Vec3 gyradii{};
    double seaDensity = kRhoSeawater;
    double zetaHeave = 0, zetaRoll = 0, zetaPitch = 0;
    Vec3 addedMass{};
    Vec3 addedInertia{};
    bool hasDeckEdgeZ = false, hasMass = false, hasDraft = false;
    bool hasCog = false, hasGyradii = false, hasSeaDensity = false;
    bool hasDamping = false, hasAddedMass = false, hasAddedInertia = false;
};

struct PendingCompartment {
    std::string name;
    int declaredAt = 0, boxAt = 0;
    Vec3 lo{}, hi{};
    double permeability = 0;
    bool vented = false;
    bool hasBox = false, hasPermeability = false, hasVented = false;
};

struct PendingOpening {
    std::string name, a, b;
    int declaredAt = 0, betweenAt = 0;
    Vec3 pos{};
    double area = 0, discharge = 0;
    OpeningKind kind = OpeningKind::Breach;
    bool open = false;
    bool hasBetween = false, hasAt = false, hasArea = false;
    bool hasDischarge = false, hasKind = false, hasOpen = false;
};

struct PendingPump {
    std::string name, compartment;
    int declaredAt = 0, compartmentAt = 0;
    double capacity = 0, maxHead = 0;
    bool on = false;
    bool hasCompartment = false, hasCapacity = false, hasMaxHead = false, hasOn = false;
};

}  // namespace

bool parseShipFile(std::string_view text, const std::string& origin, ShipDefinition& out,
                   std::string& error) {
    PendingShip shipBlock;
    std::vector<double> waterlines;
    std::vector<Station> stations;
    std::vector<PendingCompartment> compartments;
    std::vector<PendingOpening> openings;
    std::vector<PendingPump> pumps;

    bool sawFormat = false, sawShip = false, sawHull = false, sawWaterlines = false;
    int hullAt = 0;
    Block block = Block::None;

    std::vector<std::string_view> tokens;
    int lineNumber = 0;
    std::size_t cursor = 0;

    while (cursor <= text.size()) {
        const std::size_t newline = text.find('\n', cursor);
        const std::string_view line =
            text.substr(cursor, newline == std::string_view::npos ? std::string_view::npos
                                                                  : newline - cursor);
        cursor = newline == std::string_view::npos ? text.size() + 1 : newline + 1;
        ++lineNumber;

        tokenise(line, tokens);
        if (tokens.empty()) continue;
        const std::string key(tokens[0]);

        // --- helpers, all of which report and return false -------------------

        const auto fail = [&](const std::string& why) {
            error = at(origin, lineNumber, why);
            return false;
        };
        const auto arity = [&](std::size_t want, const char* shape) {
            if (tokens.size() == want) return true;
            return fail("'" + key + "' takes " + shape);
        };
        // Setting a key twice in one block is a typo, not an override: the second
        // value wins in silence and the first line reads as though it took effect.
        const auto once = [&](bool& flag) {
            if (flag) return fail("'" + key + "' is set twice in this block");
            flag = true;
            return true;
        };
        const auto num = [&](std::size_t index, double& value) {
            if (parseNumber(tokens[index], value)) return true;
            return fail("'" + std::string(tokens[index]) + "' is not a finite number");
        };
        const auto vec = [&](std::size_t index, Vec3& value) {
            return num(index, value.x) && num(index + 1, value.y) && num(index + 2, value.z);
        };
        // (0, hi]. For the quantities where zero is a broken definition rather
        // than a value: an opening with no area, a space with no permeability, a
        // radius of gyration of zero -- which is an axis with no inertia and
        // therefore unbounded angular acceleration.
        const auto positive = [&](double value, double hi) {
            if (value > 0 && value <= hi) return true;
            return fail("'" + key + "' must be within (0, " + number(hi) + "], got " +
                        number(value));
        };
        // [0, hi]. For the coefficients a ship may legitimately zero out. Roll
        // damping is the case that decides this: when a radiation model supplies
        // the damping, the modal stand-in has to go to zero or it is counted
        // twice, which is a defect this repo has already shipped once.
        const auto nonNegative = [&](double value, double hi) {
            if (value >= 0 && value <= hi) return true;
            return fail("'" + key + "' must be within [0, " + number(hi) + "], got " +
                        number(value));
        };
        const auto flag = [&](std::size_t index, bool& value) {
            if (tokens[index] == "true") { value = true; return true; }
            if (tokens[index] == "false") { value = false; return true; }
            return fail("'" + std::string(tokens[index]) + "' is not true or false");
        };
        // A name that parses as a number is a shifted field, not a name --
        // `pump 0.060` is someone's capacity landing where the name goes. Left
        // alone it produces a pump nothing will ever reference by name, and the
        // damage-control scenarios address everything by name.
        const auto declare = [&](std::size_t index, std::string& value) {
            double ignored = 0;
            if (parseNumber(tokens[index], ignored))
                return fail("'" + std::string(tokens[index]) + "' is a number, not a name");
            value = std::string(tokens[index]);
            return true;
        };

        // --- the format stamp, which must come first -------------------------

        if (key == "ship_format") {
            if (sawFormat) return fail("'ship_format' is declared twice");
            if (!arity(2, "one version number")) return false;
            double version = 0;
            if (!num(1, version)) return false;
            if (!(version >= 1) || !(version <= 1e6) || version != std::floor(version))
                return fail("'ship_format' takes a positive whole version number");
            if (version > kShipFormatVersion)
                return fail("ship_format " + number(version) +
                            " is newer than this build, which understands " +
                            number(kShipFormatVersion));
            sawFormat = true;
            continue;
        }
        if (!sawFormat)
            return fail("a ship file must begin with 'ship_format " +
                        number(kShipFormatVersion) + "', not '" + key + "'");

        // --- block headers ---------------------------------------------------

        if (key == "ship") {
            if (sawShip) return fail("a file describes exactly one ship");
            if (!arity(2, "one name")) return false;
            if (!declare(1, shipBlock.name)) return false;
            shipBlock.declaredAt = lineNumber;
            sawShip = true;
            block = Block::Ship;
            continue;
        }
        if (key == "hull") {
            if (sawHull) return fail("a file describes exactly one hull");
            if (!arity(1, "no arguments")) return false;
            sawHull = true;
            hullAt = lineNumber;
            block = Block::Hull;
            continue;
        }
        if (key == "compartment") {
            if (!arity(2, "one name")) return false;
            PendingCompartment c;
            if (!declare(1, c.name)) return false;
            if (c.name == "sea")
                return fail("'sea' is the reserved name for the water outside the hull");
            for (const PendingCompartment& earlier : compartments)
                if (earlier.name == c.name)
                    return fail("compartment '" + c.name + "' is defined twice");
            c.declaredAt = lineNumber;
            compartments.push_back(c);
            block = Block::Compartment;
            continue;
        }
        if (key == "opening") {
            if (!arity(2, "one name")) return false;
            PendingOpening o;
            if (!declare(1, o.name)) return false;
            for (const PendingOpening& earlier : openings)
                if (earlier.name == o.name)
                    return fail("opening '" + o.name + "' is defined twice");
            o.declaredAt = lineNumber;
            openings.push_back(o);
            block = Block::Opening;
            continue;
        }
        if (key == "pump") {
            if (!arity(2, "one name")) return false;
            PendingPump p;
            if (!declare(1, p.name)) return false;
            for (const PendingPump& earlier : pumps)
                if (earlier.name == p.name) return fail("pump '" + p.name + "' is defined twice");
            p.declaredAt = lineNumber;
            pumps.push_back(p);
            block = Block::Pump;
            continue;
        }

        if (block == Block::None) return fail("'" + key + "' outside any block");

        // --- keys, by block --------------------------------------------------

        if (block == Block::Ship) {
            if (key == "deck_edge_z") {
                if (!once(shipBlock.hasDeckEdgeZ) || !arity(2, "one number") ||
                    !num(1, shipBlock.deckEdgeZ))
                    return false;
            } else if (key == "lightship_mass") {
                if (!once(shipBlock.hasMass) || !arity(2, "one number in kg") ||
                    !num(1, shipBlock.lightshipMass) ||
                    !positive(shipBlock.lightshipMass, 1e12))
                    return false;
            } else if (key == "lightship_draft") {
                if (!once(shipBlock.hasDraft) || !arity(2, "one number in m") ||
                    !num(1, shipBlock.lightshipDraft) ||
                    !positive(shipBlock.lightshipDraft, 1e4))
                    return false;
            } else if (key == "lightship_cog") {
                if (!once(shipBlock.hasCog) || !arity(4, "three numbers") ||
                    !vec(1, shipBlock.cog))
                    return false;
            } else if (key == "gyradii") {
                if (!once(shipBlock.hasGyradii) || !arity(4, "three numbers") ||
                    !vec(1, shipBlock.gyradii))
                    return false;
                for (int i = 0; i < 3; ++i)
                    if (!positive(shipBlock.gyradii[i], 1e4)) return false;
            } else if (key == "sea_density") {
                if (!once(shipBlock.hasSeaDensity) || !arity(2, "one number in kg/m3") ||
                    !num(1, shipBlock.seaDensity) || !positive(shipBlock.seaDensity, 1e5))
                    return false;
            } else if (key == "damping") {
                if (!once(shipBlock.hasDamping) || !arity(4, "three numbers") ||
                    !num(1, shipBlock.zetaHeave) || !num(2, shipBlock.zetaRoll) ||
                    !num(3, shipBlock.zetaPitch))
                    return false;
                for (double z : {shipBlock.zetaHeave, shipBlock.zetaRoll, shipBlock.zetaPitch})
                    if (!nonNegative(z, 1.0)) return false;
            } else if (key == "added_mass") {
                if (!once(shipBlock.hasAddedMass) || !arity(4, "three numbers") ||
                    !vec(1, shipBlock.addedMass))
                    return false;
                for (int i = 0; i < 3; ++i)
                    if (!nonNegative(shipBlock.addedMass[i], 1e3)) return false;
            } else if (key == "added_inertia") {
                if (!once(shipBlock.hasAddedInertia) || !arity(4, "three numbers") ||
                    !vec(1, shipBlock.addedInertia))
                    return false;
                for (int i = 0; i < 3; ++i)
                    if (!nonNegative(shipBlock.addedInertia[i], 1e3)) return false;
            } else {
                return fail("unknown key '" + key + "' in a ship block");
            }
            continue;
        }

        if (block == Block::Hull) {
            if (key == "waterlines") {
                if (!once(sawWaterlines)) return false;
                if (tokens.size() < 3) return fail("'waterlines' takes at least two heights");
                for (std::size_t i = 1; i < tokens.size(); ++i) {
                    double z = 0;
                    if (!num(i, z)) return false;
                    // Strictly ascending: a repeated or reversed waterline gives
                    // makeHullFromStations a degenerate or inside-out strake, and
                    // the mesh it returns is still closed -- so nothing downstream
                    // complains, the displacement is simply wrong.
                    if (!waterlines.empty() && z <= waterlines.back())
                        return fail("waterlines must ascend: " + number(z) + " follows " +
                                    number(waterlines.back()));
                    waterlines.push_back(z);
                }
            } else if (key == "station") {
                if (!sawWaterlines)
                    return fail("'station' needs the waterline list, which must come first");
                // The count is checked before anything is read, because
                // makeHullFromStations silently substitutes a half-breadth of
                // zero for a missing one: a station one column short becomes a
                // knife edge at the top waterline and the mesh still closes.
                if (tokens.size() != waterlines.size() + 2) {
                    // tokens.size() is at least 1 here; the subtraction is on
                    // std::size_t, so clamp rather than wrap on a bare 'station'.
                    const std::size_t given = tokens.size() > 2 ? tokens.size() - 2 : 0;
                    return fail("'station' takes an x followed by one half-breadth per waterline"
                                " -- expected " + number(static_cast<double>(waterlines.size())) +
                                ", got " + number(static_cast<double>(given)));
                }
                Station s;
                if (!num(1, s.x)) return false;
                if (!stations.empty() && s.x <= stations.back().x)
                    return fail("stations must ascend in x: " + number(s.x) + " follows " +
                                number(stations.back().x));
                for (std::size_t i = 0; i < waterlines.size(); ++i) {
                    double hb = 0;
                    if (!num(i + 2, hb)) return false;
                    if (hb < 0) return fail("half-breadths cannot be negative, got " + number(hb));
                    s.halfBeam.push_back(hb);
                }
                stations.push_back(s);
            } else {
                return fail("unknown key '" + key + "' in a hull block");
            }
            continue;
        }

        if (block == Block::Compartment) {
            PendingCompartment& c = compartments.back();
            if (key == "box") {
                if (!once(c.hasBox) || !arity(7, "six numbers: lo x y z then hi x y z") ||
                    !vec(1, c.lo) || !vec(4, c.hi))
                    return false;
                for (int i = 0; i < 3; ++i)
                    if (c.hi[i] <= c.lo[i])
                        return fail("'box' has hi <= lo on the " + std::string(1, "xyz"[i]) +
                                    " axis");
                c.boxAt = lineNumber;
            } else if (key == "permeability") {
                if (!once(c.hasPermeability) || !arity(2, "one fraction") ||
                    !num(1, c.permeability) || !positive(c.permeability, 1.0))
                    return false;
            } else if (key == "vented") {
                if (!once(c.hasVented) || !arity(2, "true or false") || !flag(1, c.vented))
                    return false;
            } else {
                return fail("unknown key '" + key + "' in a compartment block");
            }
            continue;
        }

        if (block == Block::Opening) {
            PendingOpening& o = openings.back();
            if (key == "between") {
                if (!once(o.hasBetween) || !arity(3, "two space names") || !declare(1, o.a) ||
                    !declare(2, o.b))
                    return false;
                if (o.a == o.b) return fail("'between' connects a space to itself");
                o.betweenAt = lineNumber;
            } else if (key == "at") {
                if (!once(o.hasAt) || !arity(4, "three numbers") || !vec(1, o.pos)) return false;
            } else if (key == "area") {
                // Strictly positive. Negative is nonsense the flow solver would
                // run with; zero is an opening that does nothing, and an author
                // who wants that has `open false`, which says so.
                if (!once(o.hasArea) || !arity(2, "one number in m2") || !num(1, o.area) ||
                    !positive(o.area, 1e4))
                    return false;
            } else if (key == "discharge") {
                if (!once(o.hasDischarge) || !arity(2, "one coefficient") || !num(1, o.discharge) ||
                    !positive(o.discharge, 1.0))
                    return false;
            } else if (key == "kind") {
                if (!once(o.hasKind) || !arity(2, "one of breach, door, hatch, vent, pipe"))
                    return false;
                if (tokens[1] == "breach") o.kind = OpeningKind::Breach;
                else if (tokens[1] == "door") o.kind = OpeningKind::Door;
                else if (tokens[1] == "hatch") o.kind = OpeningKind::Hatch;
                else if (tokens[1] == "vent") o.kind = OpeningKind::Vent;
                else if (tokens[1] == "pipe") o.kind = OpeningKind::Pipe;
                else return fail("'" + std::string(tokens[1]) +
                                 "' is not one of breach, door, hatch, vent, pipe");
            } else if (key == "open") {
                if (!once(o.hasOpen) || !arity(2, "true or false") || !flag(1, o.open))
                    return false;
            } else {
                return fail("unknown key '" + key + "' in an opening block");
            }
            continue;
        }

        // Block::Pump. The space a pump empties is `drains`, not `compartment`,
        // so that `compartment` means one thing everywhere: it opens a block. A
        // key that opens a block in one context and sets a field in another makes
        // the compartment declared immediately after a pump unparseable.
        PendingPump& p = pumps.back();
        if (key == "drains") {
            if (!once(p.hasCompartment) || !arity(2, "one compartment name") ||
                !declare(1, p.compartment))
                return false;
            p.compartmentAt = lineNumber;
        } else if (key == "capacity") {
            if (!once(p.hasCapacity) || !arity(2, "one number in m3/s") || !num(1, p.capacity) ||
                !nonNegative(p.capacity, 1e4))
                return false;
        } else if (key == "max_head") {
            if (!once(p.hasMaxHead) || !arity(2, "one number in m") || !num(1, p.maxHead) ||
                !positive(p.maxHead, 1e4))
                return false;
        } else if (key == "on") {
            if (!once(p.hasOn) || !arity(2, "true or false") || !flag(1, p.on)) return false;
        } else {
            return fail("unknown key '" + key + "' in a pump block");
        }
    }

    // --- completeness ---------------------------------------------------------
    //
    // Checked here rather than as each block closes, so the diagnostic points at
    // the block header -- which is where the author has to go anyway -- and so a
    // file truncated mid-block is refused for the same reason as one that simply
    // forgot a key.
    const auto missing = [&](int line, const std::string& what) {
        error = at(origin, line, what);
        return false;
    };

    if (!sawFormat) {
        error = origin + ": empty, or not a ship file -- no 'ship_format' directive";
        return false;
    }
    if (!sawShip) {
        error = origin + ": no 'ship' block";
        return false;
    }
    if (!sawHull) {
        error = origin + ": no 'hull' block";
        return false;
    }

    if (!shipBlock.hasDeckEdgeZ)
        return missing(shipBlock.declaredAt, "ship '" + shipBlock.name + "' has no deck_edge_z");
    if (!shipBlock.hasCog)
        return missing(shipBlock.declaredAt, "ship '" + shipBlock.name + "' has no lightship_cog");
    if (!shipBlock.hasGyradii)
        return missing(shipBlock.declaredAt, "ship '" + shipBlock.name + "' has no gyradii");
    if (shipBlock.hasMass && shipBlock.hasDraft)
        return missing(shipBlock.declaredAt, "ship '" + shipBlock.name +
                                                 "' states both lightship_mass and"
                                                 " lightship_draft; they would disagree");
    if (!shipBlock.hasMass && !shipBlock.hasDraft)
        return missing(shipBlock.declaredAt, "ship '" + shipBlock.name +
                                                 "' has neither lightship_mass nor"
                                                 " lightship_draft");

    if (!sawWaterlines) return missing(hullAt, "hull has no waterline list");
    if (stations.size() < 2) return missing(hullAt, "hull needs at least two stations");

    for (const PendingCompartment& c : compartments) {
        if (!c.hasBox)
            return missing(c.declaredAt, "compartment '" + c.name + "' has no box");
        if (!c.hasPermeability)
            return missing(c.declaredAt, "compartment '" + c.name + "' has no permeability");
        if (!c.hasVented)
            return missing(c.declaredAt, "compartment '" + c.name + "' does not say whether it is"
                                         " vented; a sealed space traps air and floods far more"
                                         " slowly, so this cannot be left to a default");
    }
    for (const PendingOpening& o : openings) {
        if (!o.hasBetween) return missing(o.declaredAt, "opening '" + o.name + "' has no between");
        if (!o.hasAt) return missing(o.declaredAt, "opening '" + o.name + "' has no at");
        if (!o.hasArea) return missing(o.declaredAt, "opening '" + o.name + "' has no area");
        if (!o.hasDischarge)
            return missing(o.declaredAt, "opening '" + o.name + "' has no discharge");
        if (!o.hasKind) return missing(o.declaredAt, "opening '" + o.name + "' has no kind");
        if (!o.hasOpen) return missing(o.declaredAt, "opening '" + o.name + "' has no open");
    }
    for (const PendingPump& p : pumps) {
        if (!p.hasCompartment)
            return missing(p.declaredAt, "pump '" + p.name + "' does not say what it drains");
        if (!p.hasCapacity) return missing(p.declaredAt, "pump '" + p.name + "' has no capacity");
        if (!p.hasMaxHead) return missing(p.declaredAt, "pump '" + p.name + "' has no max_head");
        if (!p.hasOn) return missing(p.declaredAt, "pump '" + p.name + "' has no on");
    }

    // --- geometry -------------------------------------------------------------

    ShipDefinition scratch;
    scratch.name = shipBlock.name;
    Ship& ship = scratch.ship;

    ship.hull = makeHullFromStations(stations, waterlines);
    // The manifold check is here because CLAUDE.md's defect table has it finding
    // a hull wound inconsistently that was 40% heavy on displacement, and every
    // volume integral in this engine returns nonsense in silence on a mesh that
    // fails it. Offsets that describe a self-intersecting form get caught here,
    // at load, with a line number, rather than as a wrong number in a report.
    if (!isClosedManifold(ship.hull))
        return missing(hullAt, "the offsets do not produce a closed, consistently wound hull");
    const double hullVolume = integrate(ship.hull).volume;
    if (!(hullVolume > 0)) return missing(hullAt, "the offsets enclose no volume");

    // The offsets table *is* the envelope, from keel to weather deck, so a deck
    // edge outside its range is not a tall ship -- it is a transposed digit or a
    // unit mistake, and it would show up only as a freeboard nobody could check.
    // `freeboardMin` is what the prototype's verdict turns on.
    const double keelZ = waterlines.front();
    const double topZ = waterlines.back();
    if (shipBlock.deckEdgeZ < keelZ || shipBlock.deckEdgeZ > topZ)
        return missing(shipBlock.declaredAt,
                       "deck_edge_z " + number(shipBlock.deckEdgeZ) +
                           " m is outside the hull, which runs from " + number(keelZ) + " to " +
                           number(topZ) + " m");

    ship.deckEdgeZ = shipBlock.deckEdgeZ;
    ship.lightshipCog = shipBlock.cog;
    ship.gyradii = shipBlock.gyradii;
    if (shipBlock.hasSeaDensity) ship.seaDensity = shipBlock.seaDensity;
    if (shipBlock.hasDamping) {
        ship.zetaHeave = shipBlock.zetaHeave;
        ship.zetaRoll = shipBlock.zetaRoll;
        ship.zetaPitch = shipBlock.zetaPitch;
    }
    if (shipBlock.hasAddedMass) {
        ship.addedMassSurge = shipBlock.addedMass.x;
        ship.addedMassSway = shipBlock.addedMass.y;
        ship.addedMassHeave = shipBlock.addedMass.z;
    }
    if (shipBlock.hasAddedInertia) {
        ship.addedInertiaRoll = shipBlock.addedInertia.x;
        ship.addedInertiaPitch = shipBlock.addedInertia.y;
        ship.addedInertiaYaw = shipBlock.addedInertia.z;
    }

    // `lightship_draft` keeps the hull form and the loading consistent: light
    // displacement is *defined* as whatever floats her at that draft with every
    // space dry, so editing the offsets moves the mass with them. Stating a mass
    // instead pins it, which is what a real stability booklet does.
    //
    // A draft outside the hull's own z range is refused from both ends: below it
    // there is nothing in the water and the ship weighs nothing, above it the
    // whole envelope is submerged at the light condition, which is not a ship.
    if (shipBlock.hasDraft &&
        (shipBlock.lightshipDraft <= keelZ || shipBlock.lightshipDraft > topZ))
        return missing(shipBlock.declaredAt,
                       "lightship_draft " + number(shipBlock.lightshipDraft) +
                           " m puts no part of the hull in the water, or all of it -- the hull"
                           " runs from " + number(keelZ) + " to " + number(topZ) + " m");
    ship.lightshipMass = shipBlock.hasMass
                             ? shipBlock.lightshipMass
                             : integrateBelowPlane(ship.hull, Vec3{0, 0, 1},
                                                   shipBlock.lightshipDraft)
                                       .volume *
                                   ship.seaDensity;
    if (!(ship.lightshipMass > 0))
        return missing(shipBlock.declaredAt, "the ship's lightship mass works out non-positive");

    ship.compartments.reserve(compartments.size());
    for (const PendingCompartment& pc : compartments) {
        Compartment c;
        c.name = pc.name;
        c.mesh = clipToBox(ship.hull, pc.lo, pc.hi);
        c.permeability = pc.permeability;
        c.ventedToAtmosphere = pc.vented;
        // A box that misses the hull clips away to nothing. That is the
        // prototype's oldest authoring mistake and it is completely silent: the
        // compartment exists, has zero volume, never floods, and the ship simply
        // survives damage it should not have.
        if (c.mesh.tris.empty())
            return missing(pc.boxAt, "compartment '" + pc.name +
                                         "' clips to nothing -- its box does not"
                                         " overlap the hull");
        if (!isClosedManifold(c.mesh))
            return missing(pc.boxAt, "compartment '" + pc.name +
                                         "' clips to a mesh that is not a closed manifold");
        if (!(integrate(c.mesh).volume > 0))
            return missing(pc.boxAt, "compartment '" + pc.name + "' clips to zero volume");
        ship.compartments.push_back(std::move(c));
    }

    // --- references -----------------------------------------------------------
    //
    // Resolved after the whole file, so a mod may declare an opening before the
    // compartment it reaches. A name that never turns up is an error naming the
    // line that used it -- `Ship::validate()` would catch a dangling *index*, but
    // by then the name is gone and there is nothing useful left to say.
    const auto resolve = [&](const std::string& spaceName, int line, const std::string& what,
                             int& index) {
        if (spaceName == "sea") {
            index = kSea;
            return true;
        }
        index = ship.findCompartment(spaceName);
        if (index == kSea) {
            error = at(origin, line, what + " names compartment '" + spaceName +
                                         "', which does not exist");
            return false;
        }
        return true;
    };

    ship.openings.reserve(openings.size());
    for (const PendingOpening& po : openings) {
        Opening o;
        o.name = po.name;
        if (!resolve(po.a, po.betweenAt, "opening '" + po.name + "'", o.a)) return false;
        if (!resolve(po.b, po.betweenAt, "opening '" + po.name + "'", o.b)) return false;
        o.pos = po.pos;
        o.area = po.area;
        o.dischargeCoeff = po.discharge;
        o.kind = po.kind;
        o.open = po.open;
        ship.openings.push_back(std::move(o));
    }

    ship.pumps.reserve(pumps.size());
    for (const PendingPump& pp : pumps) {
        Pump p;
        p.name = pp.name;
        if (!resolve(pp.compartment, pp.compartmentAt, "pump '" + pp.name + "'", p.compartment))
            return false;
        if (p.compartment == kSea)
            return missing(pp.compartmentAt,
                           "pump '" + pp.name + "' cannot draw from the sea; it pumps a"
                                                " compartment overboard");
        p.capacity = pp.capacity;
        p.maxHead = pp.maxHead;
        p.on = pp.on;
        ship.pumps.push_back(std::move(p));
    }

    // Committed. Nothing above this line has touched the caller's definition.
    out = std::move(scratch);
    return true;
}

bool loadShipFile(const std::string& path, ShipDefinition& out, std::string& error) {
    std::ifstream file(path, std::ios::binary);
    if (!file) {
        error = path + ": cannot open";
        return false;
    }
    std::ostringstream contents;
    contents << file.rdbuf();
    if (file.bad()) {
        error = path + ": read failed";
        return false;
    }
    return parseShipFile(contents.str(), path, out, error);
}

}  // namespace sim
