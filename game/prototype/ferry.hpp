// SPDX-License-Identifier: MIT
// Slice-1 test vessel: a 120 m ro-pax ferry, subdivided into 14 watertight
// compartments over two decks, with a vehicle deck that has no subdivision at all.
// The vehicle deck is the point: it is the single largest free surface a ferry can
// have, and it is why ro-pax casualties capsize instead of settling.
#pragma once

#include "../../engine/sim/ship.hpp"

namespace game {

sim::Ship buildFerry();

}  // namespace game
