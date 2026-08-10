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

}  // namespace game
