// SPDX-License-Identifier: MIT
//
// Slice 1, headless. Hole the ferry and watch what the physics does about it.
//
//   ./shipsim [--scenario=none|doors|full] [--duration=900] [--dt=0.01] [--csv=path]
//
// There is no renderer here on purpose. Everything the eventual game shows on a
// damage-control board is already decided by this loop; getting it right in a
// terminal first is much cheaper than getting it right behind a Vulkan swapchain.
#include "ferry.hpp"

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <string>
#include <vector>

using namespace sim;

namespace {

struct Event {
    double time;
    std::string label;
    std::function<void(Ship&)> action;
    bool fired = false;
};

void setOpening(Ship& s, std::string_view name, bool open) {
    for (Opening& o : s.openings)
        if (o.name == name) o.open = open;
}

void setPump(Ship& s, std::string_view name, bool on) {
    for (Pump& p : s.pumps)
        if (p.name == name) p.on = on;
}

std::vector<Event> buildScenario(const std::string& which) {
    std::vector<Event> events;
    if (which == "none") return events;

    // Realistic reaction times: alarm, muster, then someone physically at the
    // door. Forty-five seconds is optimistic.
    events.push_back({45.0, "close watertight door between engine rooms",
                      [](Ship& s) { setOpening(s, "wt_door_er", false); }});
    events.push_back({60.0, "start bilge pumps (ER port/stbd, aft hold stbd)",
                      [](Ship& s) {
                          setPump(s, "bilge_er_s", true);
                          setPump(s, "bilge_er_p", true);
                          setPump(s, "bilge_ah_s", true);
                      }});

    if (which == "full") {
        // Counterflood early. Once a list develops, the port sea suctions lift
        // clear of the water and the option disappears -- which the 'doors'
        // scenario demonstrates by losing it.
        events.push_back({90.0, "open port wing tank sea suctions (counterflood)",
                          [](Ship& s) {
                              setOpening(s, "cf_valve_aft_p", true);
                              setOpening(s, "cf_valve_fwd_p", true);
                          }});
        // The decisive action on any ro-pax: keep the sea off the vehicle deck.
        // Nothing else on the damage control board comes close in value.
        events.push_back({120.0, "secure vehicle deck shell doors and ramp",
                          [](Ship& s) {
                              setOpening(s, "downflood_ramp_s", false);
                              setOpening(s, "downflood_port_s", false);
                              setOpening(s, "downflood_port_p", false);
                          }});
        events.push_back({300.0, "shut counterflooding valves",
                          [](Ship& s) {
                              setOpening(s, "cf_valve_aft_p", false);
                              setOpening(s, "cf_valve_fwd_p", false);
                          }});
    }
    return events;
}

void printGzCurve(const Ship& s, const char* title) {
    std::printf("\n  %s\n  heel deg :", title);
    for (int d = 0; d <= 60; d += 5) std::printf("%7d", d);
    std::printf("\n  GZ     m :");
    for (int d = 0; d <= 60; d += 5)
        std::printf("%7.3f", s.rightingArmAtHeel(d * kDegToRad, 0.0));
    std::printf("\n");
}

void printCompartments(const Ship& s) {
    std::printf("  %-20s %10s %8s %8s %9s\n", "compartment", "gross m3", "fill %", "P kPa", "water t");
    for (const Compartment& c : s.compartments) {
        if (c.waterVolume < 1.0 && c.fillFraction() < 1e-4) continue;
        std::printf("  %-20s %10.0f %8.1f %8.1f %9.0f\n", c.name.c_str(), c.grossVolume,
                    100.0 * c.fillFraction(), c.airPressure / 1000.0,
                    c.waterVolume * kRhoSeawater / 1000.0);
    }
}

}  // namespace

int main(int argc, char** argv) {
    std::string scenario = "doors";
    double duration = 900.0;
    // The flooding solution is converged well below this: heel, trim, draft and
    // floodwater agree to four significant figures from dt = 0.005 to dt = 0.04.
    // The system is stiff in restoring force but heavily damped, and the orifice
    // flows change far more slowly than the step, so 0.02 s is a safety margin
    // rather than a requirement.
    double dt = 0.02;
    std::string csvPath;

    for (int i = 1; i < argc; ++i) {
        const std::string_view a = argv[i];
        if (a.starts_with("--scenario=")) scenario = std::string(a.substr(11));
        else if (a.starts_with("--duration=")) duration = std::atof(argv[i] + 11);
        else if (a.starts_with("--dt=")) dt = std::atof(argv[i] + 5);
        else if (a.starts_with("--csv=")) csvPath = std::string(a.substr(6));
        else { std::fprintf(stderr, "unknown argument: %s\n", argv[i]); return 2; }
    }

    Ship ship = game::buildFerry();
    const double seaLevel = 0.0;
    ship.initialise(seaLevel);

    for (const std::string& problem : ship.validate())
        std::fprintf(stderr, "  ship definition: %s\n", problem.c_str());

    const Diagnostics intact = ship.diagnostics(seaLevel);
    std::printf("=== 120 m ro-pax ferry, intact condition ===\n");
    std::printf("  displacement      %10.0f t\n", intact.displacementMass / 1000.0);
    std::printf("  draft (midship)   %10.2f m\n", intact.draftMidship);
    const double lengthOverall = ship.hullHi.x - ship.hullLo.x;
    const double beam = ship.hullHi.y - ship.hullLo.y;
    std::printf("  L x B x D         %6.0f x %.0f x %.1f m\n", lengthOverall, beam,
                ship.hullHi.z - ship.hullLo.z);
    std::printf("  block coeff Cb    %10.3f\n",
                intact.buoyantVolume / (lengthOverall * beam * intact.draftMidship));
    std::printf("  waterplane Cwp    %10.3f\n", intact.waterplaneArea / (lengthOverall * beam));
    std::printf("  waterplane area   %10.0f m2\n", intact.waterplaneArea);
    std::printf("  KG                %10.2f m\n", intact.centreOfGravity.z);
    std::printf("  KB                %10.2f m\n", intact.centreOfBuoyancy.z);
    std::printf("  GM (transverse)   %10.2f m\n", intact.gmTransverse);
    std::printf("  freeboard to veh. deck %5.2f m\n", intact.freeboardMin);
    printGzCurve(ship, "Intact righting arm curve");

    auto events = buildScenario(scenario);
    std::printf("\n=== Scenario '%s': 2.4 m2 breach, starboard engine room, 2.5 m below waterline ===\n",
                scenario.c_str());
    if (events.empty()) std::printf("  (no damage control action taken)\n");
    for (const Event& e : events) std::printf("  t+%-6.0f %s\n", e.time, e.label.c_str());

    std::ofstream csv;
    if (!csvPath.empty()) {
        csv.open(csvPath);
        csv << "t,draft,heel_deg,trim_deg,gm_m,gz_m,floodwater_t,freeboard_m,displacement_t\n";
    }

    std::printf("\n  %8s %8s %8s %8s %8s %8s %11s %10s\n",
                "t s", "draft m", "heel", "trim", "GM m", "GZ m", "flood t", "freeb m");

    double t = 0.0;
    double nextReport = 0.0;
    const char* outcome = "still afloat at end of run";

    while (t < duration) {
        for (Event& e : events) {
            if (!e.fired && t >= e.time) {
                e.fired = true;
                e.action(ship);
                std::printf("  -- t+%.0fs  %s\n", t, e.label.c_str());
            }
        }

        if (t >= nextReport) {
            const Diagnostics d = ship.diagnostics(seaLevel);
            std::printf("  %8.0f %8.2f %8.2f %8.2f %8.2f %8.3f %11.0f %10.2f\n",
                        t, d.draftMidship, d.heelDeg, d.trimDeg, d.gmTransverse,
                        d.gzRighting, d.floodwaterMass / 1000.0, d.freeboardMin);
            if (csv)
                csv << t << ',' << d.draftMidship << ',' << d.heelDeg << ',' << d.trimDeg
                    << ',' << d.gmTransverse << ',' << d.gzRighting << ','
                    << d.floodwaterMass / 1000.0 << ',' << d.freeboardMin << ','
                    << d.displacementMass / 1000.0 << '\n';

            if (std::abs(d.heelDeg) > 60.0) { outcome = "CAPSIZED"; break; }
            if (!d.afloat)                  { outcome = "FOUNDERED"; break; }
            nextReport += 15.0;
        }

        ship.step(dt, seaLevel);
        t += dt;
    }

    const Diagnostics fin = ship.diagnostics(seaLevel);
    // "Still afloat" is not a verdict. A ship lying at 55 degrees with negative GM
    // and water still coming in has already been lost; it just has not finished.
    if (std::string_view(outcome).starts_with("still")) {
        if (fin.gmTransverse < 0 && std::abs(fin.heelDeg) > 20.0)
            outcome = "LOST - lolled over with negative GM, flooding continuing";
        else if (fin.gmTransverse < 0)
            outcome = "LOST - negative GM, loll imminent";
        else if (fin.freeboardMin < 0)
            outcome = "SURVIVED but the deck edge is under; no margin left";
        else
            outcome = "SURVIVED - positive GM, deck edge dry";
    }
    std::printf("\n=== Outcome at t+%.0fs: %s ===\n", t, outcome);
    std::printf("  heel %.1f deg, trim %.1f deg, draft %.2f m, %.0f t of floodwater aboard\n",
                fin.heelDeg, fin.trimDeg, fin.draftMidship, fin.floodwaterMass / 1000.0);
    std::printf("  effective GM %.2f m (intact was %.2f m)\n\n",
                fin.gmTransverse, intact.gmTransverse);
    printCompartments(ship);
    printGzCurve(ship, "Damaged righting arm curve");
    return 0;
}
