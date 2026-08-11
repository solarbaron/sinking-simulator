// SPDX-License-Identifier: MIT
// Slice-1 test vessel: a 120 m ro-pax ferry, subdivided into 18 watertight
// compartments over four levels -- 0-1.8 m double bottom, 1.8-7.0 m holds and
// machinery, 7.0-12.5 m vehicle deck, 12.5-15.0 m above -- with a vehicle deck
// that has no subdivision at all.
//
// It said 14 over two decks until the mid wing tanks were authored, which is the
// same gap that let 41% of a ram amidships tear open onto no compartment at all.
// The vehicle deck is the point: it is the single largest free surface a ferry can
// have, and it is why ro-pax casualties capsize instead of settling.
#pragma once

#include "../../engine/sim/ship.hpp"

namespace game {

sim::Ship buildFerry();

// --- The verdict a run ends on ------------------------------------------------
//
// Two tools draw a survival verdict from one `sim::Diagnostics` -- `shipsim` in a
// sentence, `ram_view` in a word -- and each did it with its own copy of
// `gmTransverse < 0`. Neither copy was reachable from `shipsim_tests`, so the
// verdict was the one thing in the stability chain nothing tested, and it read a
// GM that the ship itself had already flagged as unusable with exactly the
// confidence of one it had not.
//
// They live here rather than in `sim` because a verdict is prose about a
// scenario, not physics; and here rather than in either tool because this is the
// only translation unit `shipsim`, `ram_view` and `shipsim_tests` all compile.
//
// **`UNDETERMINED` is a refusal, not a third degree of danger.** When
// `sim::judgeStability` comes back `Unresolved` there is no metacentric height to
// be negative or positive, so no survival verdict is stated at all; what the
// reader gets instead is the reason, and the observations that do not depend on
// GM -- heel, trim, draft, floodwater, and whether the deck edge is under -- on
// the lines that follow.
//
// Each returns one of a small fixed set of strings and never formats a number
// into one: `scripts/verify.sh` compares the compiled ferry's outcome with
// `ships/ferry.ship`'s as exact strings, and a verdict carrying a float would
// make that a comparison of two round-offs.

// The full sentence, for the final state of a `shipsim` run.
const char* floodingOutcome(const sim::Diagnostics& d);

// The same judgement in one word: LOST, SURVIVED or UNDETERMINED.
const char* stabilityWord(const sim::Diagnostics& d);

}  // namespace game
