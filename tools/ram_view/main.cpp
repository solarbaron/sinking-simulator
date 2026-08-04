// SPDX-License-Identifier: MIT
//
// The Phase 3 milestone, in one act: **ram the ferry.**
//
// `docs/06-roadmap.md` states it as "the hull deforms, tears where the stress
// says it should, and the resulting hole floods at a rate the hole's own area
// determines." Every piece of that exists and each was built and validated
// separately. Nothing had ever run them as one event, which is where the
// interesting failures live.
//
// The chain, and not a step of it is reimplemented here:
//
//   collision   two closed hulls overlap; the penetration volume gives a force,
//               a patch, and the energy that went into the meeting
//   indentation that energy spent outward panel by panel; membrane stretching
//               says how deep it got and which bays tore
//   breach      torn panels become openings in the flooding network, merged,
//               placed at their own centroid, connected by geometry
//   ship        the openings flood at Cd A sqrt(2 dp / rho), the free surfaces
//               re-level, the trapped air compresses, and she lists or does not
//
// The striking ship is rigid here. She is not: a real bow crushes and takes a
// large share of the energy, and `indentation.hpp` records what leaving that out
// costs. This over-states the damage, knowingly, and the run prints enough to see
// by how much.
//
// **No damage control is applied.** Nobody closes a door, starts a pump or
// counterfloods, so this is Phase 0's `none` scenario with the breach computed
// rather than authored -- and `none` loses her too. She is lost at every speed and
// every aiming point tried, which is consistent rather than suspicious: the
// scenarios that let her live are the ones where somebody acts. What does respond
// is *how* she is lost. Struck amidships she takes 2100 t and lolls 6 degrees;
// struck at the quarter she takes 8100 t and goes over 26 the other way.
//
// It also shows where the outcome stops caring about the strike. From 1.5 to
// 6 m/s the hole grows from 3.4 to 113 m2 and the floodwater barely moves,
// because a 3.4 m2 breach fills the compartment behind it inside 900 s just as a
// 113 m2 one does. Beyond a threshold, damage stability is decided by *which*
// compartments are opened and not by how big the hole is -- which is exactly what
// the subdivision rules are written around, and it falls out here rather than
// being put in.
//
//   ./ram_view [--speed=M_PER_S] [--aim=X_METRES] [--reach=METRES]
//              [--duration=SECONDS] [--out=DIR] [--frames=N]
#include "engine/sim/breach.hpp"
#include "engine/sim/collision.hpp"
#include "engine/sim/hullform.hpp"
#include "engine/sim/indentation.hpp"
#include "engine/sim/scantlings.hpp"
#include "game/prototype/ferry.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

namespace {

struct Options {
    double speed = 6.0;       // m/s closing
    double aim = 0.0;         // m along the ferry, where the bow lands
    double reach = 12.0;      // m, how far damage may propagate from the impact
    double duration = 900.0;  // s of flooding after the strike
    std::string out;
    int frames = 0;
};

bool parse(int argc, char** argv, Options& o) {
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        const auto value = [&](const char* key) -> const char* {
            const std::string prefix = std::string("--") + key + "=";
            return a.rfind(prefix, 0) == 0 ? a.c_str() + prefix.size() : nullptr;
        };
        if (const char* v = value("speed")) o.speed = std::atof(v);
        else if (const char* v = value("aim")) o.aim = std::atof(v);
        else if (const char* v = value("reach")) o.reach = std::atof(v);
        else if (const char* v = value("duration")) o.duration = std::atof(v);
        else if (const char* v = value("out")) o.out = v;
        else if (const char* v = value("frames")) o.frames = std::atoi(v);
        else {
            std::printf("unknown argument: %s\n", a.c_str());
            return false;
        }
    }
    return true;
}

// A striking ship: a smaller hull, built from particulars like any other.
sim::Ship striker() {
    sim::HullParticulars p;
    p.lengthPp = 90.0;
    p.beam = 15.0;
    p.draft = 5.5;
    p.depth = 11.0;
    p.blockCoefficient = 0.68;
    p.midshipCoefficient = 0.97;
    p.parallelMiddleBodyFraction = 0.25;
    p.stationCount = 41;

    sim::Ship s;
    s.hull = sim::makeHullFromParticulars(p);
    s.deckEdgeZ = p.depth;
    s.lightshipMass = p.blockCoefficient * p.lengthPp * p.beam * p.draft * sim::kRhoSeawater;
    s.lightshipCog = {0.0, 0.0, 0.55 * p.depth};
    s.gyradii = {0.35 * p.beam, 0.25 * p.lengthPp, 0.25 * p.lengthPp};
    return s;
}

}  // namespace

int main(int argc, char** argv) {
    Options options;
    if (!parse(argc, argv, options)) return 2;

    sim::Ship ferry = game::buildFerry();
    ferry.initialise(0.0);
    const sim::Scantlings scantlings = sim::ferryScantlings();
    const sim::StructuralMesh structure = sim::makeStructuralMesh(ferry.hull, scantlings);

    sim::Ship bow = striker();
    bow.initialise(0.0);

    const sim::Diagnostics before = ferry.diagnostics(0.0);
    std::printf("=== ram the ferry ===\n");
    std::printf("struck : 120 m ro-pax, %.0f t, GM %.2f m, %zu structural panels\n",
                before.displacementMass / 1000.0, before.gmTransverse, structure.panels.size());
    std::printf("striker: 90 m hull, %.0f t, closing at %.1f m/s, aimed at x = %+.0f m\n",
                bow.diagnostics(0.0).displacementMass / 1000.0, options.speed, options.aim);

    // Place the striker abeam, on the ferry's starboard side, closing.
    const double standoff = 0.5 * (ferry.hullHi.y - ferry.hullLo.y) + 45.0;
    bow.state.position = {options.aim, -standoff, bow.state.position.z};
    bow.state.orientation = sim::Quat::fromAxisAngle({0, 0, 1}, sim::kPi / 2.0);
    bow.state.velocity = {0.0, options.speed, 0.0};

    // --- 1. The strike -------------------------------------------------------
    const sim::Sea sea(0.0);
    const double dt = 0.01;
    sim::ContactMaterial material;
    sim::ContactHistory history;

    bool touched = false;
    double peakForce = 0;
    for (int i = 0; i < 4000 && (!touched || history.duration > 0); ++i) {
        const sim::HullContact contact = sim::applyContact(ferry, bow, material, dt, &history);
        if (contact.touching) {
            touched = true;
            const sim::ContactBody a = sim::contactBodyOf(ferry), b = sim::contactBodyOf(bow);
            peakForce = std::max(peakForce, length(sim::contactLoad(contact, a, b, material).force));
        } else if (touched) {
            break;   // through and clear
        }
        ferry.step(dt, sea);
        bow.step(dt, sea);
    }

    if (!touched) {
        std::printf("\nno contact -- she was missed\n");
        return 0;
    }
    std::printf("\ncontact: %.2f s, peak %.0f MN, %.0f MJ taken out of the two ships\n",
                history.duration, peakForce / 1e6, history.work / 1e6);

    // --- 2. What that energy did to the plating ------------------------------
    const sim::Vec3 impact = history.loadAtPeak.point - ferry.state.position;
    const sim::ImpactDamage damage =
        sim::impactDamage(structure, impact, options.reach, history.work, scantlings);
    std::printf("damage : %.3f m into her over %zu bays, %zu torn, %.1f m2 of hole\n",
                damage.penetration, damage.panels.size(), damage.torn.size(), damage.tornArea);
    if (damage.energyUnspent > 0)
        std::printf("         %.0f MJ unspent -- reach bounded the answer, the hull did not\n",
                    damage.energyUnspent / 1e6);

    // --- 3. The hole -------------------------------------------------------
    const sim::BreachSet breaches = sim::breachesFromFailedPanels(ferry, structure, damage.torn);
    double openArea = 0;
    for (const sim::Breach& b : breaches.breaches) openArea += b.opening.area;
    std::printf("breach : %zu opening(s), %.2f m2 reaching a compartment\n",
                breaches.breaches.size(), openArea);
    for (const std::string& problem : breaches.problems)
        std::printf("         ! %s\n", problem.c_str());
    sim::applyBreaches(ferry, breaches);

    // --- 4. Whether she lives ------------------------------------------------
    std::printf("\n%8s %10s %8s %8s %8s\n", "t (s)", "flood t", "heel", "trim", "GM");
    const auto steps = static_cast<int>(options.duration / dt);
    for (int i = 1; i <= steps; ++i) {
        ferry.step(dt, sea);
        if (i % (steps / 6) == 0 || i == steps) {
            const sim::Diagnostics d = ferry.diagnostics(sea);
            std::printf("%8.0f %10.0f %7.2f° %7.2f° %8.2f\n", i * dt,
                        d.floodwaterMass / 1000.0, d.heelDeg, d.trimDeg, d.gmTransverse);
        }
    }

    const sim::Diagnostics after = ferry.diagnostics(sea);
    std::printf("\noutcome: %s -- %.0f t of water, heel %.1f deg, GM %.2f m\n",
                (!after.afloat || after.gmTransverse < 0) ? "LOST" : "SURVIVED",
                after.floodwaterMass / 1000.0, after.heelDeg, after.gmTransverse);
    std::printf("ok\n");
    return 0;
}
