// SPDX-License-Identifier: MIT
//
// What `GasPromoter` decides on the real ferry with a real fire in her.
//
// **Why this exists.** `GasPromoter` is exercised by `tests/test_promotion.cpp`
// and by nothing else: `grep -rl GasPromoter engine tools tests` returns the
// header, its own translation unit, and that one test file. Every compartment it
// has ever judged was constructed by a test to have the property under test.
//
// That is precisely the gap that hid every defect in the *water* tier. The unit
// tests there were thorough and green while the criterion compared a speed
// against an acceleration threshold, read roll off a world axis, and carried a
// budget pair that made the whole tier unreachable. None of it was visible until
// the criterion was run against a real ship, because a fixture built to produce
// a given number cannot disagree with the code that reads it.
//
// So this is the water tier's `water_probe`, pointed at the gas tier. It decides
// nothing and changes nothing: it drives `fire::Model` on the ferry, reviews it,
// and reports what the criterion saw -- including the things the tests cannot
// see, which are the resolution a promoted compartment would actually get and
// whether the two cell limits are in agreement.

#include "../../engine/sim/fire.hpp"
#include "../../engine/sim/les.hpp"
#include "../../engine/sim/promotion.hpp"
#include "../../engine/sim/ship.hpp"
#include "../../game/prototype/ferry.hpp"

#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

using namespace sim;

int main(int argc, char** argv) {
    double duration = 600.0;
    double power = 4.0e6;        // W
    double reviewEvery = 10.0;   // s of model time

    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a.starts_with("--duration=")) duration = std::atof(a.c_str() + 11);
        else if (a.starts_with("--power=")) power = std::atof(a.c_str() + 8) * 1.0e6;
        else if (a.starts_with("--review=")) reviewEvery = std::atof(a.c_str() + 9);
        // Without this the loop had no `else`, so a mistyped flag was silently
        // dropped and the run proceeded on its defaults -- and this tool's output
        // is compared against published figures, so the quiet outcome is a
        // different experiment handed back under the right name.
        else {
            std::printf("unknown option %s\n", a.c_str());
            return 2;
        }
    }

    Ship ship = game::buildFerry();
    const Sea sea;
    ship.initialise(sea);

    fire::Model model;
    model.attach(ship, {ship.findCompartment("engine_room_s"),
                        ship.findCompartment("engine_room_p"),
                        ship.findCompartment("vehicle_deck")});

    fire::DesignFire blaze;
    blaze.name = "machinery";
    blaze.compartment = model.findGas("engine_room_s");
    blaze.baseZ = 2.5;
    blaze.diameter = 2.5;
    blaze.growthCoefficient = fire::kGrowthFast;
    blaze.peakHeatRelease = power;
    blaze.steadyDuration = duration * 2.0;
    model.fires.push_back(blaze);

    for (const std::string& problem : model.validate())
        std::fprintf(stderr, "  fire model: %s\n", problem.c_str());

    promotion::GasPromoter promoter;
    const auto& crit = promoter.criterion();

    std::printf("--- what the gas promoter sees on a burning ferry ---\n");
    std::printf("    %.0f MW in engine_room_s, %.0f s, reviewed every %.0f s\n",
                power / 1e6, duration, reviewEvery);
    std::printf("    promote at spread >= %g or rise >= %g K, dwell %d, hold %d\n",
                crit.spreadPromote, crit.risePromote, crit.dwell, crit.hold);
    std::printf("    cellBudget %d across all compartments; grid.maxCells %d each\n\n",
                crit.cellBudget, crit.grid.maxCells);

    // **The two cell limits are set independently and are denominated in the same
    // unit**, which is the shape that made the water tier's budgets disagree by
    // 6.25x. Here they are equal, and equal is its own problem: one compartment
    // may consume the entire global budget, so the second compartment to qualify
    // is refused however small it is. Reported rather than asserted, because
    // whether that is wrong depends on how many compartments are meant to run at
    // once -- but it cannot be reasoned about while nobody has printed it.
    if (crit.cellBudget == crit.grid.maxCells)
        std::printf("    note: the global budget equals the per-compartment ceiling,\n"
                    "          so one compartment can consume all of it\n\n");

    double maxSpread = 0, maxRise = 0;
    int reviews = 0, everQualified = 0, budgetRefusals = 0;
    double firstPromotionAt = -1;
    int peakActive = 0, peakCells = 0;
    std::vector<std::pair<std::string, int>> promotedEver;

    double t = 0, nextReview = 0;
    while (t < duration) {
        const double step = std::min(2.0, duration - t);
        model.step(step, ship, sea);
        t += step;

        if (t < nextReview) continue;
        nextReview = t + reviewEvery;

        const promotion::GasReview r = promoter.review(model);
        reviews++;

        for (const auto& c : r.considered) {
            if (c.spread > maxSpread) maxSpread = c.spread;
            if (c.rise > maxRise) maxRise = c.rise;
            if (c.score > 0) everQualified++;
        }
        for (const auto& p : r.promoted) {
            if (firstPromotionAt < 0) firstPromotionAt = t;
            bool seen = false;
            for (auto& e : promotedEver)
                if (e.first == p.name) { e.second++; seen = true; }
            if (!seen) promotedEver.push_back({p.name, 1});
        }
        // The refusal nobody reads. In the water tier this was silent for 154
        // reviews and looked exactly like hysteresis working.
        for (const auto& problem : r.problems)
            if (problem.find("budget") != std::string::npos) budgetRefusals++;

        if (static_cast<int>(promoter.active().size()) > peakActive)
            peakActive = static_cast<int>(promoter.active().size());
        if (r.cellsActive > peakCells) peakCells = r.cellsActive;
    }

    std::printf("    %d reviews over %.0f s\n\n", reviews, duration);
    std::printf("    peak spread seen         %10.3f      (threshold %g)\n",
                maxSpread, crit.spreadPromote);
    std::printf("    peak layer rise seen     %10.1f K    (threshold %g)\n",
                maxRise, crit.risePromote);
    std::printf("    candidate-reviews qualifying          %8d\n", everQualified);
    std::printf("    promotions                            %8d\n", promoter.promotions());
    std::printf("    demotions                             %8d\n", promoter.demotions());
    if (firstPromotionAt >= 0)
        std::printf("    first promotion at t = %.0f s\n", firstPromotionAt);
    std::printf("    peak compartments active              %8d\n", peakActive);
    std::printf("    peak cells active                     %8d   (budget %d)\n",
                peakCells, crit.cellBudget);
    std::printf("    reviews that refused on budget        %8d\n", budgetRefusals);

    if (!promotedEver.empty()) {
        std::printf("\n    compartments promoted at least once:\n");
        for (const auto& e : promotedEver)
            std::printf("      %-24s %d time(s)\n", e.first.c_str(), e.second);
    }

    // **What resolution a promoted compartment actually gets**, which is the
    // question the water tier answered too late. `les::gridFor` coarsens until the
    // compartment fits `grid.maxCells`, so a large space is not refused -- it is
    // silently resolved at whatever cell size makes it fit, and a cell size
    // coarser than the physics needs produces a run that completes, conserves
    // mass, and models nothing. The vehicle deck is the case to watch, exactly as
    // it was for water.
    std::printf("\n    what a resolved compartment would actually get:\n");
    std::printf("      %-24s %10s %10s %12s\n", "compartment", "cells", "cell (m)", "of budget");
    for (const auto& gas : model.gas) {
        const les::Grid grid = les::gridFor(gas, crit.grid);
        const int cells = grid.n[0] * grid.n[1] * grid.n[2];
        std::printf("      %-24s %10d %10.2f %11.0f%%\n", gas.name.c_str(), cells,
                    grid.h[0], 100.0 * cells / std::max(1, crit.cellBudget));
    }

    std::printf("\n    verdict: ");
    if (promoter.promotions() == 0)
        std::printf("nothing promoted -- the criterion never fired on a real fire\n");
    else
        std::printf("%d promotions on a real ferry\n", promoter.promotions());
    return 0;
}
