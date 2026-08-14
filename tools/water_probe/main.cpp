// SPDX-License-Identifier: MIT
//
// What `WaterPromoter` would actually decide on a real flooding ferry, measured
// before any of it is wired into `Ship::step()`.
//
// **Why this is a tool and not a unit test.** The unit tests drive the promoter
// with a synthetic ship at a chosen roll rate, which answers "does the state
// machine work" and not "what does it do here". Those are different questions and
// the second one decides whether activating it is safe: the flooding scenarios
// carry gated published figures, and a promoter that fires during them changes
// numbers the front page publishes. This runs the promoter alongside the real
// scenario and reports what it saw, without changing a single compartment.
//
// It exists because of the entry in CLAUDE.md's table about extrapolated figures
// pointing at the wrong fix. The particle-budget arithmetic says a compartment of
// this size costs so many particles; what it cannot say is whether the criterion
// ever fires, how often, or on what.

#include "../../engine/sim/ship.hpp"
#include "../../engine/sim/water_promotion.hpp"
#include "../../engine/sim/waves.hpp"
#include "../../game/prototype/ferry.hpp"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

using namespace sim;

namespace {

void setOpening(Ship& s, std::string_view name, bool open) {
    for (Opening& o : s.openings)
        if (o.name == name) o.open = open;
}

void setPump(Ship& s, std::string_view name, bool on) {
    for (Pump& p : s.pumps)
        if (p.name == name) p.on = on;
}

// The 'doors' scenario, which is the one the README publishes a capsize for and
// therefore the one where a promoter firing would be most expensive.
void applyDamageControl(Ship& s, double t, double dt) {
    auto fires = [&](double at) { return t >= at && t - dt < at; };
    if (fires(45.0))  setOpening(s, "wt_door_er", false);
    if (fires(60.0)) {
        setPump(s, "bilge_er_s", true);
        setPump(s, "bilge_er_p", true);
        setPump(s, "bilge_ah_s", true);
    }
}

}  // namespace

// What a promoted compartment actually costs to step, in core-seconds per
// simulated second.
//
// `coreSecondsPerCompartment` was 5.0 and marked "estimate, will be measured".
// It is the number the whole tier's affordability rests on and it had never been
// run; this is what ran it. The measured answer is 27.9 at 1 m³ rising to 3030 at
// 100 m³, so the estimate was low by 5.6x to 606x and -- worse -- was a
// per-compartment *constant* where the truth scales with the water.
//
// **This comment described the pre-measurement state for one commit after the
// measurement existed**, which is the "each document quoted the previous one"
// failure in the instrument that did the measuring. It says what is true now.
//
// The compartment is a box of the given volume at the ferry's own proportions
// rather than a cube: a hold is wide and shallow, and the cell count -- which is
// what the projection actually costs -- follows the shape and not just the
// volume.
void measureStepCost(double dt, double h) {
    std::printf("--- what a promoted compartment costs to step ---\n");
    std::printf("    dt = %g s, h = %g m, %s\n\n", dt, h,
                "cost is core-seconds per simulated second");
    std::printf("      %8s %9s %8s %10s %12s %10s\n",
                "volume", "particles", "tiles", "ms/step", "core-s/sim-s", "substeps");

    for (double vol : {1.0, 5.0, 20.0, 50.0, 100.0}) {
        flip::Field field;
        // A compartment 2:1 in plan and shallow, which is a hold rather than a
        // cube: volume = 2b * b * d with d = b/4 gives b = (2 vol)^(1/3).
        const double b = std::cbrt(2.0 * vol);
        const double lo[3] = {0, 0, 0};
        const double hi[3] = {2.0 * b, b, vol / (2.0 * b * b)};

        field.grid.h = h;
        field.grid.lo[0] = lo[0]; field.grid.lo[1] = lo[1]; field.grid.lo[2] = lo[2];
        for (int a = 0; a < 3; ++a)
            field.grid.n[a] = std::max(1, static_cast<int>(std::ceil((hi[a] - lo[a]) / h)));

        flip::seedBox(field, lo, hi, 2, kRhoSeawater);
        flip::setTotalMass(field, vol * kRhoSeawater);

        flip::Params params;
        flip::Solver solver;
        flip::Account account;

        // One step to build the sparse structure, so the timing measures steady
        // state rather than first-touch allocation.
        solver.step(field, dt, params, account);

        const int steps = 5;
        const auto t0 = std::chrono::steady_clock::now();
        int substeps = 0;
        for (int i = 0; i < steps; ++i)
            substeps += solver.step(field, dt, params, account).substeps;
        const auto t1 = std::chrono::steady_clock::now();

        const double msPerStep =
            std::chrono::duration<double, std::milli>(t1 - t0).count() / steps;
        // Core-seconds of work per second of simulated time: a step covers dt of
        // model time, so the ratio is the wall time over dt.
        const double coreSecondsPerSimSecond = (msPerStep / 1000.0) / dt;

        std::printf("      %6.0f m3 %9d %8d %9.2f %12.2f %10.1f\n",
                    vol, static_cast<int>(field.particles.size()), solver.tiles(),
                    msPerStep, coreSecondsPerSimSecond,
                    static_cast<double>(substeps) / steps);
    }
    // **The cost per simulated second is flat in dt, and that is the point.**
    // Doubling dt to the production 0.02 doubles the substeps (2 -> 4) and the
    // wall time with them, leaving core-s/sim-s within 4%: the CFL condition sets
    // the work, not the timestep. So this is not a cost that can be stepped
    // around, and a caller who wants it cheaper has to change `h`.
    std::printf("\n    Cost per simulated second is flat in dt -- doubling dt doubles the\n"
                "    substeps and the wall time together, because the CFL condition sets\n"
                "    the work. `coreSecondsPerCompartment` defaults to the 1 m3 figure.\n");

    // **The same volume in a different shape, because the claim under test is
    // that this tier can serve "deep narrow compartments".** Every row above is
    // one aspect ratio (2:1 in plan, shallow), and a cost model that only ever
    // saw one shape cannot distinguish "cost follows volume" from "cost follows
    // the shape I happened to pick". A tall narrow tank and a shallow wide deck
    // of equal volume hold the same water and are not the same problem.
    std::printf("\n      the same 20 m3 in three shapes -- volume alone does not set the cost\n");
    std::printf("      %-22s %9s %8s %10s %12s\n", "shape", "particles", "tiles",
                "ms/step", "core-s/sim-s");
    struct Shape { const char* name; double x, y, z; };
    const double v = 20.0;
    const Shape shapes[] = {
        // A wing tank: tall, narrow, deep. What the tier claims to serve.
        {"deep narrow (1x1x20)",  std::sqrt(v / 20.0), std::sqrt(v / 20.0), 20.0},
        // The probe's own default, for comparison against the table above.
        {"hold (2:1, shallow)",   2.0 * std::cbrt(2.0 * v) / 2.0 * 2.0 / 2.0, 0, 0},
        // A vehicle deck: shallow and very wide. What it cannot.
        {"deck-like (0.25 deep)", 0, 0, 0.25},
    };
    for (const Shape& s : shapes) {
        double ex = s.x, ey = s.y, ez = s.z;
        if (std::strcmp(s.name, "hold (2:1, shallow)") == 0) {
            const double b = std::cbrt(2.0 * v);
            ex = 2.0 * b; ey = b; ez = v / (2.0 * b * b);
        } else if (std::strcmp(s.name, "deck-like (0.25 deep)") == 0) {
            // 2:1 in plan at a fixed 0.25 m depth.
            const double area = v / 0.25;
            ey = std::sqrt(area / 2.0); ex = 2.0 * ey; ez = 0.25;
        }

        flip::Field field;
        const double lo[3] = {0, 0, 0};
        const double hi[3] = {ex, ey, ez};
        field.grid.h = h;
        field.grid.lo[0] = 0; field.grid.lo[1] = 0; field.grid.lo[2] = 0;
        for (int a = 0; a < 3; ++a)
            field.grid.n[a] = std::max(1, static_cast<int>(std::ceil((hi[a] - lo[a]) / h)));
        flip::seedBox(field, lo, hi, 2, kRhoSeawater);
        flip::setTotalMass(field, v * kRhoSeawater);

        flip::Params params;
        flip::Solver solver;
        flip::Account account;
        solver.step(field, dt, params, account);
        const auto t0 = std::chrono::steady_clock::now();
        for (int i = 0; i < 3; ++i) solver.step(field, dt, params, account);
        const auto t1 = std::chrono::steady_clock::now();
        const double msPerStep =
            std::chrono::duration<double, std::milli>(t1 - t0).count() / 3.0;

        std::printf("      %-22s %9d %8d %10.2f %12.2f\n", s.name,
                    static_cast<int>(field.particles.size()), solver.tiles(),
                    msPerStep, (msPerStep / 1000.0) / dt);
    }
}

int main(int argc, char** argv) {
    double duration = 900.0;
    double dt = 0.01;
    double reviewEvery = 1.0;   // s of model time between reviews
    // The beam-sea control. A promoter that fires nowhere is broken rather than
    // safe, so the flooding run needs a counterpart that *should* promote: a
    // regular beam wave at the ferry's own roll period, which is the one forcing
    // guaranteed to produce large roll rates on an intact ship.
    double waveAmplitude = 0;   // m; 0 means flat water, the flooding scenario
    double wavePeriod = 0;      // s; 0 means take the ship's own roll period
    bool costOnly = false;      // --cost: measure the step cost and nothing else
    // The cell size the cost is measured at. `promoteWater` uses 0.05 and that is
    // where the tier is unaffordable, so the question this answers is whether any
    // coarser h is affordable *and* still resolves anything -- see --cost.
    double cellSize = 0.05;     // m

    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a.starts_with("--duration=")) duration = std::atof(a.c_str() + 11);
        else if (a.starts_with("--dt=")) dt = std::atof(a.c_str() + 5);
        else if (a.starts_with("--review=")) reviewEvery = std::atof(a.c_str() + 9);
        else if (a.starts_with("--wave=")) waveAmplitude = std::atof(a.c_str() + 7);
        else if (a.starts_with("--period=")) wavePeriod = std::atof(a.c_str() + 9);
        else if (a == "--cost") costOnly = true;
        else if (a.starts_with("--h=")) cellSize = std::atof(a.c_str() + 4);
    }
    // After the whole argument list, not during it: `--cost` acted immediately
    // once, so `--cost --dt=0.02` measured at the default 0.01 and printed a
    // header saying so. An option whose meaning depends on where it appears in
    // the line is a trap.
    if (costOnly) { measureStepCost(dt, cellSize); return 0; }

    Ship ship = game::buildFerry();
    ship.initialise(0.0);

    // Beam seas: the wave runs along +y, across the ship, which is the heading
    // that puts the most roll into a hull and the one a ro-pax capsizes in.
    WaveField waves;
    Sea sea(0.0);
    if (waveAmplitude > 0) {
        if (wavePeriod <= 0) {
            // Her own roll period, so the forcing is at resonance rather than at
            // an arbitrary frequency: T = 2*pi*sqrt(k_roll^2 / (g GM)).
            const Diagnostics d0 = ship.diagnostics(0.0);
            const double kRoll = ship.gyradii.x;
            wavePeriod = 2.0 * kPi * std::sqrt(kRoll * kRoll / (kGravity * d0.gmTransverse));
        }
        const double omega = 2.0 * kPi / wavePeriod;
        waves = WaveField::regular(waveAmplitude, omega, kPi / 2.0);
        sea.waves = &waves;

        // Water to promote. An intact ship has none, and a criterion gated on
        // `minVolume` cannot fire on an empty compartment however hard she rolls
        // -- so the control would come back "no promotions" for a reason that has
        // nothing to do with the motion it is meant to be testing.
        const int idx = ship.findCompartment("vehicle_deck");
        if (idx >= 0) {
            ship.compartments[idx].waterVolume =
                0.05 * ship.compartments[idx].floodableVolume();
            std::printf("    seeded %s with %.0f m3 so there is something to promote\n",
                        ship.compartments[idx].name.c_str(),
                        ship.compartments[idx].waterVolume);
        }
    }

    promotion::WaterPromoter promoter;
    const auto& crit = promoter.criterion();

    std::printf("--- what the water promoter sees on a flooding ferry ---\n");
    if (waveAmplitude > 0)
        std::printf("    beam sea: %.2f m amplitude at T = %.2f s (her roll period),"
                    " %.0f s at dt = %g\n", waveAmplitude, wavePeriod, duration, dt);
    else
        std::printf("    scenario 'doors', %.0f s at dt = %g, reviewed every %.2f s\n",
                    duration, dt, reviewEvery);
    std::printf("    promote at rollRate >= %g rad/s or accel >= %g m/s2,"
                " dwell %d, hold %d\n",
                crit.rollRatePromote, crit.accelPromote, crit.dwell, crit.hold);
    std::printf("    minVolume %g m3, particle budget %d, tile budget %d\n\n",
                crit.minVolume, crit.particleBudget, crit.tileBudget);

    // The extremes of what the criterion is reading, so a criterion that never
    // fires can be told apart from one that fires on everything. A promoter that
    // reports zero promotions is only meaningful next to the range of its input.
    double maxRoll = 0, maxAccel = 0;
    double maxRollAt = 0, maxAccelAt = 0;
    int reviews = 0, everQualified = 0;
    double firstPromotionAt = -1;
    int peakActive = 0, peakParticles = 0, peakTiles = 0;
    // **A review that refuses to promote says so in `problems`, and nothing was
    // reading it.** A criterion that qualifies 2525 times and promotes once is
    // either dwelling or being refused, and those are opposite findings: the
    // first is the hysteresis working, the second is a budget too small to run
    // the thing it is budgeting for. Without this the two are indistinguishable
    // from the outside, which is the silent-failure shape this repo keeps
    // finding.
    int budgetRefusals = 0;
    // Every compartment that was ever promoted, and how many reviews it held.
    std::vector<std::pair<std::string, int>> promotedEver;

    double t = 0, nextReview = 0;
    while (t < duration) {
        // The beam-sea control is an *intact* ship being rolled; damage control
        // belongs only to the flooding scenario.
        if (waveAmplitude <= 0) applyDamageControl(ship, t, dt);
        sea.time = t;
        ship.step(dt, sea);
        t += dt;

        if (t < nextReview) continue;
        nextReview = t + reviewEvery;

        const promotion::WaterReview r = promoter.review(ship, reviewEvery);
        reviews++;

        for (const auto& c : r.considered) {
            if (c.rollRate > maxRoll) { maxRoll = c.rollRate; maxRollAt = t; }
            if (c.accel > maxAccel)   { maxAccel = c.accel;   maxAccelAt = t; }
            if (c.score > 0) everQualified++;
        }
        for (const auto& p : r.promoted) {
            if (firstPromotionAt < 0) firstPromotionAt = t;
            bool seen = false;
            for (auto& e : promotedEver)
                if (e.first == p.name) { e.second++; seen = true; }
            if (!seen) promotedEver.push_back({p.name, 1});
        }
        if (static_cast<int>(promoter.active().size()) > peakActive)
            peakActive = static_cast<int>(promoter.active().size());
        if (r.particlesActive > peakParticles) peakParticles = r.particlesActive;
        if (r.tilesActive > peakTiles) peakTiles = r.tilesActive;
        for (const auto& p : r.problems)
            if (p.find("budget") != std::string::npos) budgetRefusals++;
    }

    const Diagnostics d = ship.diagnostics(sea);
    std::printf("    the run itself: heel %.1f deg, %.0f t of floodwater, GM %.2f m\n",
                d.heelDeg, d.floodwaterMass / 1000.0, d.gmTransverse);
    std::printf("    %d reviews over %.0f s\n\n", reviews, duration);

    std::printf("    peak roll rate seen      %10.4f rad/s at t = %.0f s"
                "   (threshold %g)\n", maxRoll, maxRollAt, crit.rollRatePromote);
    std::printf("    peak 'accel' seen        %10.4f m/s2  at t = %.0f s"
                "   (threshold %g)\n", maxAccel, maxAccelAt, crit.accelPromote);
    std::printf("    candidate-reviews qualifying          %8d\n", everQualified);
    std::printf("    promotions                            %8d\n", promoter.promotions());
    std::printf("    demotions                             %8d\n", promoter.demotions());
    if (firstPromotionAt >= 0)
        std::printf("    first promotion at t = %.0f s\n", firstPromotionAt);
    std::printf("    peak compartments active              %8d\n", peakActive);
    std::printf("    peak particles active                 %8d   (budget %d)\n",
                peakParticles, crit.particleBudget);
    std::printf("    peak tiles active                     %8d   (budget %d)\n",
                peakTiles, crit.tileBudget);
    std::printf("    reviews that refused on budget        %8d\n", budgetRefusals);

    // **Which budget binds, at the cost model the promoter actually uses.** The
    // two are set independently and are not independent: `estimateFlipCost` puts
    // 1000 particles and 1/(64 h^3) = 125 tiles in every cubic metre, so the
    // volume each budget admits is fixed and one of them is always the smaller.
    // A budget that can never bind is not a limit, it is a number in a header --
    // and the particle count is the one the comment calls "the memory
    // bottleneck".
    {
        const double h = 0.05;
        const double m3PerParticleBudget = crit.particleBudget / 1000.0;
        const double m3PerTileBudget = crit.tileBudget * (64.0 * h * h * h);
        std::printf("\n    the particle budget admits %.1f m3; the tile budget admits %.1f m3\n",
                    m3PerParticleBudget, m3PerTileBudget);
        std::printf("    so %s binds first, by %.2fx -- the other cannot ever be reached\n",
                    m3PerTileBudget < m3PerParticleBudget ? "the TILE budget" : "the PARTICLE budget",
                    m3PerParticleBudget > m3PerTileBudget
                        ? m3PerParticleBudget / m3PerTileBudget
                        : m3PerTileBudget / m3PerParticleBudget);
    }

    if (!promotedEver.empty()) {
        std::printf("\n    compartments promoted at least once:\n");
        for (const auto& e : promotedEver)
            std::printf("      %-24s %d time(s)\n", e.first.c_str(), e.second);
    }

    // **The verdict this tool exists to give.** Activating the promoter inside
    // `Ship::step()` is safe for the published figures exactly when it never
    // fires during them, and that is a measurement rather than an argument.
    std::printf("\n    verdict: ");
    if (promoter.promotions() == 0)
        std::printf("nothing promoted -- activating this changes no published figure\n");
    else
        std::printf("%d promotions -- activation WOULD move the gated scenario figures\n",
                    promoter.promotions());
    return 0;
}
