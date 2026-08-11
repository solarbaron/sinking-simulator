// SPDX-License-Identifier: MIT
//
// Slice 1, headless. Hole the ferry and watch what the physics does about it.
//
//   ./shipsim [--scenario=none|doors|full] [--duration=900] [--dt=0.01] [--csv=path]
//              [--ship=path] [--gm-detail]
//
// `--gm-detail` prints how the metacentric height at the end of the run was
// arrived at -- the angle, the halvings, whether it converged, the layer that set
// the angle, and the slope at the fixed +/-0.03 rad this used to be taken at.
// One quantity per line with stable labels, because the numbers the front page
// publishes about that repair were otherwise re-derivable only by writing C++.
//
// There is no renderer here on purpose. Everything the eventual game shows on a
// damage-control board is already decided by this loop; getting it right in a
// terminal first is much cheaper than getting it right behind a Vulkan swapchain.
//
// `--ship` loads a ship definition file instead of calling buildFerry(). It is
// how a data mod is run, and it is how `verify.sh` checks that `ships/ferry.ship`
// reaches the same outcomes over 900 s as the compiled ferry -- the unit suite
// can only afford the first 150 s of that.
#include "ferry.hpp"

#include "../../engine/sim/shipfile.hpp"

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <string>
#include <utility>
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

// --- `--gm-detail` -------------------------------------------------------------
//
// The metacentric height's own working, one quantity per line, so that the claims
// made about it are re-derivable without writing C++ against the library. Every
// line is `gm-detail: <label> <value>`; the labels are an interface and do not
// move.
//
// **The last line is a control and it is deliberately the wrong answer.**
// `gm_at_fixed_0.03rad_m` is the current code asked for the slope at the angle
// this repository used to take unconditionally -- not a remembered number from
// before the change, which would be a historical note rather than a measurement.
// It is what makes "a fixed sample is wrong by 4.7 m on this ship" reproducible,
// and it is here rather than on `Diagnostics` for exactly the reason it is
// interesting: an engine that returns a GM taken at an angle nobody checked is
// how this went wrong the first time. A caller has to ask for it by name, on a
// tool, with `--gm-detail` on the command line.
//
// `pockets_at_rad` and `fixed_sample_sees_frac` come from `sim::largestFreeSurface`
// -- the wedge closed form evaluated on the layer's own measured geometry -- and
// are what say *why* a fixed angle is wrong rather than merely that it is.
void printGmDetail(const Ship& s, const Diagnostics& d, const Sea& sea) {
    constexpr double kFixed = 0.03;
    const double atFixed =
        (s.rightingArmAtHeel(kFixed, sea) - s.rightingArmAtHeel(-kFixed, sea)) / (2 * kFixed);
    const FreeSurfaceLayer layer = largestFreeSurface(s);

    std::printf("\ngm-detail: gm_converged_m %.6f\n", d.gmTransverse);
    std::printf("gm-detail: gm_sampled_at_rad %.6e\n", d.gmSampledAtRad);
    std::printf("gm-detail: gm_halvings %d\n", d.gmHalvings);
    std::printf("gm-detail: gm_converged %s\n", d.gmSlopeConverged ? "yes" : "no");
    std::printf("gm-detail: gm_at_fixed_0.03rad_m %.6f\n", atFixed);
    std::printf("gm-detail: layer_compartment %s\n",
                layer.compartment < 0
                    ? "none"
                    : s.compartments[static_cast<std::size_t>(layer.compartment)].name.c_str());
    std::printf("gm-detail: layer_water_t %.6f\n",
                layer.compartment < 0 ? 0.0
                                      : s.compartments[static_cast<std::size_t>(layer.compartment)]
                                                .waterVolume * s.seaDensity / 1000.0);
    std::printf("gm-detail: layer_plan_area_m2 %.3f\n", layer.planArea);
    std::printf("gm-detail: layer_depth_m %.6e\n", layer.depth);
    std::printf("gm-detail: layer_breadth_m %.3f\n", layer.breadth);
    std::printf("gm-detail: pockets_at_rad %.6e\n", layer.pocketingRad);
    std::printf("gm-detail: fixed_sample_sees_frac %.6f\n", layer.momentFractionAt(kFixed));
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
    std::string shipPath;
    // Bilge keel length, m. Zero leaves the ship on the linear stand-in; any
    // positive value attaches Ikeda's viscous roll damping instead. A ro-pax of
    // this size carries keels around a third of her length.
    double bilgeKeelLength = 0.0;
    // Print the metacentric height's own working at the end of the run. Off by
    // default because the answer is `gmTransverse` and everything here is about
    // *how it was arrived at* -- see printGmDetail().
    bool gmDetail = false;

    for (int i = 1; i < argc; ++i) {
        const std::string_view a = argv[i];
        if (a.starts_with("--scenario=")) scenario = std::string(a.substr(11));
        else if (a.starts_with("--duration=")) duration = std::atof(argv[i] + 11);
        else if (a.starts_with("--dt=")) dt = std::atof(argv[i] + 5);
        else if (a.starts_with("--csv=")) csvPath = std::string(a.substr(6));
        else if (a.starts_with("--ship=")) shipPath = std::string(a.substr(7));
        else if (a.starts_with("--bilge-keels=")) bilgeKeelLength = std::atof(argv[i] + 14);
        else if (a == "--gm-detail") gmDetail = true;
        else { std::fprintf(stderr, "unknown argument: %s\n", argv[i]); return 2; }
    }

    Ship ship;
    if (shipPath.empty()) {
        ship = game::buildFerry();
    } else {
        // Fails closed, and says why. A ship that would not load must not be
        // silently replaced by the compiled one -- that is how a broken mod ships.
        ShipDefinition definition;
        std::string error;
        if (!loadShipFile(shipPath, definition, error)) {
            std::fprintf(stderr, "%s\n", error.c_str());
            return 2;
        }
        std::printf("  loaded ship '%s' from %s\n", definition.name.c_str(), shipPath.c_str());
        ship = std::move(definition.ship);
    }
    const double seaLevel = 0.0;
    ship.initialise(seaLevel);

    // Real viscous roll damping, on request. Roll is the mode that decides
    // whether a damaged ship lolls or lies over, and Ikeda's B44 is strongly
    // amplitude-dependent where the linear stand-in is not -- so the two differ
    // most exactly where this scenario spends its time.
    if (bilgeKeelLength > 0.0) {
        const RollDampingHull form = ship.attachRollDamping(5.5, bilgeKeelLength, 0.9);
        for (const std::string& problem : validateRollDamping(form, {}))
            std::fprintf(stderr, "  roll damping: %s\n", problem.c_str());
    }

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
    // An intact ship has no free surface to pocket and is sampled at the full
    // 0.03 rad, so this is a metacentric height on every ship this prototype
    // builds. It is not guaranteed to be one on a ship loaded from a file with
    // water already in a compartment, and an unresolved GM must not reach the
    // reader looking like a resolved one anywhere it is printed.
    std::printf("  GM (transverse)   %10.2f m%s\n", intact.gmTransverse,
                intact.gmSlopeConverged ? "" : "   (unresolved: no linear region)");
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
        // `gm_resolved` rides with `gm_m` and not somewhere else, because a plot
        // of the GM column alone is exactly the reading this whole exercise is
        // about: the column is a metacentric height on the rows where that is 1
        // and the last thing a bisection was holding on the rows where it is 0.
        csv << "t,draft,heel_deg,trim_deg,gm_m,gm_resolved,gz_m,floodwater_t,freeboard_m,"
               "displacement_t\n";
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
            // **The note is a line of its own, and it begins with a word.** The
            // obvious place for it is the end of the row above, and that would be a
            // trapdoor: `scripts/check-figures.sh` reads this table with
            // `awk 'NF == 8 && $1 ~ /^[0-9]+$/'`, so a marker appended to a row
            // would silently *remove that row from the gate* -- exactly on the
            // steps where GM is least trustworthy. A separate line whose first
            // field is not a number is invisible to those readers by construction.
            if (!d.gmSlopeConverged)
                std::printf("      note: the GM above is not a metacentric height -- the arm is"
                            " not linear at any angle down to +/-%.2e rad\n", d.gmSampledAtRad);
            if (csv)
                csv << t << ',' << d.draftMidship << ',' << d.heelDeg << ',' << d.trimDeg
                    << ',' << d.gmTransverse << ',' << (d.gmSlopeConverged ? 1 : 0) << ','
                    << d.gzRighting << ',' << d.floodwaterMass / 1000.0 << ','
                    << d.freeboardMin << ',' << d.displacementMass / 1000.0 << '\n';

            if (std::abs(d.heelDeg) > 60.0) { outcome = "CAPSIZED"; break; }
            if (!d.afloat)                  { outcome = "FOUNDERED"; break; }
            nextReport += 15.0;
        }

        ship.step(dt, seaLevel);
        t += dt;
    }

    const Diagnostics fin = ship.diagnostics(seaLevel);
    // The verdict itself is in `game::floodingOutcome`, where a test can reach it.
    if (std::string_view(outcome).starts_with("still")) outcome = game::floodingOutcome(fin);
    std::printf("\n=== Outcome at t+%.0fs: %s ===\n", t, outcome);
    std::printf("  heel %.1f deg, trim %.1f deg, draft %.2f m, %.0f t of floodwater aboard\n",
                fin.heelDeg, fin.trimDeg, fin.draftMidship, fin.floodwaterMass / 1000.0);
    // The GM line says which of the two things it is. A slope the refinement gave
    // up on is still worth printing -- it is the only number there is -- but it is
    // printed as what it is, next to the angle it bottomed out at, and it is not
    // called a metacentric height.
    if (fin.gmSlopeConverged)
        std::printf("  effective GM %.2f m (intact was %.2f m)\n\n",
                    fin.gmTransverse, intact.gmTransverse);
    else
        std::printf("  effective GM unresolved: %.2f m is the best slope available, and it is\n"
                    "  still moving at +/-%.2e rad (intact GM was %.2f m)\n\n",
                    fin.gmTransverse, fin.gmSampledAtRad, intact.gmTransverse);
    if (gmDetail) printGmDetail(ship, fin, seaLevel);
    printCompartments(ship);
    printGzCurve(ship, "Damaged righting arm curve");
    return 0;
}
