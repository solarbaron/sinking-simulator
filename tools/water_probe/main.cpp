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

#include <cstdio>
#include <cstdlib>
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

    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a.starts_with("--duration=")) duration = std::atof(a.c_str() + 11);
        else if (a.starts_with("--dt=")) dt = std::atof(a.c_str() + 5);
        else if (a.starts_with("--review=")) reviewEvery = std::atof(a.c_str() + 9);
        else if (a.starts_with("--wave=")) waveAmplitude = std::atof(a.c_str() + 7);
        else if (a.starts_with("--period=")) wavePeriod = std::atof(a.c_str() + 9);
    }

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
    int peakActive = 0, peakParticles = 0;
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
