// SPDX-License-Identifier: MIT
//
// The Phase 2 milestone, made visible: a real ship, under its own power, in a
// real seaway.
//
// Every subsystem Phase 2 built meets here, and nothing is reimplemented for
// display:
//
//   * the hull comes from principal particulars (`hullform.hpp`), so it is a
//     named ship rather than a shape somebody drew;
//   * the sea is a directional JONSWAP spectrum (`waves.hpp`), and the *same*
//     `WaveField` drives the physics and the picture -- if those two differed,
//     the wave under the bow would not be the wave the ship responded to;
//   * buoyancy is integrated over the instantaneous wetted surface against that
//     surface, not against a mean waterline;
//   * radiation gives frequency-dependent added mass and the wave-memory force;
//   * propeller, rudder and MMG hull drive and steer it, with an autopilot on
//     the heading because the reference derivatives describe a directionally
//     unstable hull;
//   * hull and sea go through one pipeline and one depth buffer, so the ship is
//     *in* the water rather than composited next to it.
//
// The camera rides alongside at a fixed offset in the ship's own frame, because
// a world-fixed camera loses a ship making way in about twenty seconds.
//
//   ./seaway_view [--out=DIR] [--frames=N] [--hs=METRES] [--tp=SECONDS]
//                 [--ship=kvlcc2|s175] [--heading=DEGREES] [--revs=PER_SECOND]
#include "engine/core/math.hpp"
#include "engine/core/png.hpp"
#include "engine/gpu/hull.hpp"
#include "engine/gpu/ocean.hpp"
#include "engine/sim/hullform.hpp"
#include "engine/sim/ship.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

namespace {

struct Options {
    std::string out = ".";
    std::string ship = "s175";
    int frames = 48;
    double significantHeight = 4.0;
    double peakPeriod = 9.0;
    double headingDeg = 180.0;   // 180 = head seas
    // Chosen so the S-175 settles near Fn 0.275, which is where its published
    // RAOs are tabulated. Measured: 1.4 rps gives Fn 0.081, 4.0 gives 0.241,
    // 6.0 gives 0.361 -- the machinery reaches the benchmark condition, it just
    // has to be asked to.
    double revs = 4.6;
};

bool parse(int argc, char** argv, Options& o) {
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        const auto value = [&](const char* key) -> const char* {
            const std::string prefix = std::string("--") + key + "=";
            return a.rfind(prefix, 0) == 0 ? a.c_str() + prefix.size() : nullptr;
        };
        if (const char* v = value("out")) o.out = v;
        else if (const char* v = value("ship")) o.ship = v;
        else if (const char* v = value("frames")) o.frames = std::atoi(v);
        else if (const char* v = value("hs")) o.significantHeight = std::atof(v);
        else if (const char* v = value("tp")) o.peakPeriod = std::atof(v);
        else if (const char* v = value("heading")) o.headingDeg = std::atof(v);
        else if (const char* v = value("revs")) o.revs = std::atof(v);
        else {
            std::printf("unknown argument: %s\n", a.c_str());
            return false;
        }
    }
    return true;
}

// A ship built entirely from its principal particulars, with a machinery fit
// scaled from the KVLCC2 reference set by length. The manoeuvring derivatives
// are the tanker's whatever hull is asked for -- right in order of magnitude,
// wrong in detail, and stated here rather than implied.
sim::Ship buildShip(const sim::HullParticulars& p, double revs) {
    sim::Ship ship;
    std::vector<std::string> problems;
    ship.hull = sim::makeHullFromParticulars(p, &problems);
    for (const std::string& s : problems) std::printf("  note: %s\n", s.c_str());

    ship.deckEdgeZ = p.depth;
    // Float at the design draft: the hull is generated to displace Cb*L*B*T there.
    ship.lightshipMass = p.blockCoefficient * p.lengthPp * p.beam * p.draft * sim::kRhoSeawater;
    ship.lightshipCog = {0.0, 0.0, 0.55 * p.depth};
    ship.gyradii = {0.35 * p.beam, 0.25 * p.lengthPp, 0.25 * p.lengthPp};

    sim::Manoeuvring m = sim::kvlcc2();
    const double scale = p.lengthPp / 320.0;
    m.hull.length = p.lengthPp;
    m.hull.beam = p.beam;
    m.hull.draft = p.draft;
    m.hull.blockCoefficient = p.blockCoefficient;
    m.hull.xCog = p.lcbFraction * p.lengthPp;
    m.propeller.diameter *= scale;
    m.rudder.area *= scale * scale;
    m.rudder.span *= scale;
    m.rudder.x *= scale;
    m.rudder.flowStraighteningLever *= scale;
    m.rudder.hullLiftLever *= scale;
    m.revsPerSecond = revs;
    ship.propulsion = m;
    return ship;
}

// The autopilot the RAO harness uses, applied from outside. Duplicated rather
// than exported because rao.cpp keeps it internal, and a tool is the wrong
// reason to widen an API.
void holdHeading(sim::Ship& ship) {
    if (!ship.propulsion.has_value()) return;
    const sim::Mat3 R = ship.state.orientation.toMat3();
    const sim::Vec3 forward = R * sim::Vec3{1, 0, 0};
    const double heading = std::atan2(forward.y, forward.x);
    const double demand = -3.0 * heading - 30.0 * ship.state.angularVelocity.z;
    ship.propulsion->rudderAngle =
        std::clamp(demand, -35.0 * sim::kDegToRad, 35.0 * sim::kDegToRad);
}

}  // namespace

int main(int argc, char** argv) {
    Options options;
    if (!parse(argc, argv, options)) return 2;

    const sim::HullParticulars particulars =
        options.ship == "kvlcc2" ? sim::kvlcc2Particulars() : sim::s175Particulars();
    std::printf("ship: %s, Lpp %.1f m, B %.1f m, T %.1f m, Cb %.3f\n", options.ship.c_str(),
                particulars.lengthPp, particulars.beam, particulars.draft,
                particulars.blockCoefficient);

    sim::Ship ship = buildShip(particulars, options.revs);
    ship.initialise(0.0);
    const sim::HullCoefficients measured =
        sim::measureHull(ship.hull, particulars.draft, particulars.lengthPp, particulars.beam);
    std::printf("  measured Cb %.4f, displacement %.0f m3, Cwp %.3f\n", measured.blockCoefficient,
                measured.displacedVolume, measured.waterplaneCoefficient);

    const sim::RadiationTable table = ship.attachRadiation(particulars.draft, 15);
    std::printf("  radiation: %d frequencies, %d of %d solves repaired, worst residual %.3f\n",
                table.size(), table.repairedSolves, table.totalSolves, table.worstEnergyResidual);

    sim::SeaState state;
    state.significantHeight = options.significantHeight;
    state.peakPeriod = options.peakPeriod;
    // The sea travels the way the ship is *met* from: heading 180 is head seas.
    state.meanDirection = options.headingDeg * sim::kDegToRad;
    state.frequencyCount = 16;
    state.directionCount = 8;
    const sim::WaveField field(state);
    std::printf("  sea: Hs %.2f m, Tp %.1f s, %zu components, shortest wave %.1f m\n",
                field.significantHeight(), options.peakPeriod, field.components().size(),
                gpu::shortestWavelength(field));

    sim::Sea sea;
    sea.waves = &field;

    // Run the ship up to speed in still water before the sea arrives, or the
    // first frames show a stationary ship being hit by a fully developed storm.
    const double dt = 0.02;
    const sim::Sea still(0.0);
    for (int i = 0; i < 30000; ++i) {   // 600 s
        holdHeading(ship);
        ship.step(dt, still);
    }
    const sim::Mat3 R0 = ship.state.orientation.toMat3();
    std::printf("  settled at %.2f m/s (Fn %.3f)\n",
                dot(ship.state.velocity, R0 * sim::Vec3{1, 0, 0}),
                dot(ship.state.velocity, R0 * sim::Vec3{1, 0, 0}) /
                    std::sqrt(sim::kGravity * particulars.lengthPp));

    gpu::Device device;
    std::string error;
    if (!device.create(error)) {
        std::printf("no usable GPU (%s) -- running headless\n", error.c_str());
    }

    gpu::MaterialLibrary library;
    if (device.valid() &&
        !library.load(std::string(SHIPSIM_MATERIAL_DIR) + "/marine.materials", error)) {
        std::printf("materials: %s\n", error.c_str());
        return 1;
    }

    constexpr std::uint32_t kWidth = 1280, kHeight = 720;
    gpu::HullRenderer renderer;
    if (device.valid() &&
        !renderer.create(device, kWidth, kHeight, SHIPSIM_SHADER_DIR, error)) {
        std::printf("renderer: %s\n", error.c_str());
        return 1;
    }

    // Resolve the *shortest* component, not the dominant one: sixteen cells
    // across the dominant wavelength still invents a quarter of Hs, because the
    // spectrum carries waves an order of magnitude shorter than its peak.
    const double patch = 3.0 * particulars.lengthPp;
    gpu::OceanGrid grid;
    grid.halfExtent = patch;
    grid.resolution =
        std::min(512, gpu::oceanResolutionFor(patch, gpu::shortestWavelength(field), 8.0));
    std::printf("  ocean grid: %d cells over %.0f m\n", grid.resolution, 2 * patch);

    gpu::OceanSurface surface;
    gpu::SceneMesh scene;
    gpu::HullPaint paint;
    paint.waterlineZ = particulars.draft;
    gpu::SceneView view;

    const int stepsPerFrame = static_cast<int>(0.25 / dt);
    double worstHeel = 0, worstTrim = 0;
    for (int frame = 0; frame < options.frames; ++frame) {
        for (int i = 0; i < stepsPerFrame; ++i) {
            sea.time += dt;
            holdHeading(ship);
            ship.step(dt, sea);
        }

        const sim::Mat3 R = ship.state.orientation.toMat3();
        double heel = 0, trim = 0;
        sim::heelTrimFromRotation(R, heel, trim);
        worstHeel = std::max(worstHeel, std::abs(heel * sim::kRadToDeg));
        worstTrim = std::max(worstTrim, std::abs(trim * sim::kRadToDeg));

        if (!device.valid()) continue;

        // The ocean patch follows the ship, so the sea is always drawn where the
        // ship is rather than where it started.
        grid.centreX = ship.state.position.x;
        grid.centreY = ship.state.position.y;
        surface.build(field, grid, sea.time);

        scene.clear();
        if (!scene.appendShip(ship, paint, library, gpu::HullShading{}, error)) {
            std::printf("scene: %s\n", error.c_str());
            return 1;
        }
        scene.appendOcean(surface, library.find("sea_water"));

        // Ride alongside and slightly above, in the ship's frame.
        const sim::Vec3 offset = R * sim::Vec3{-0.9 * particulars.lengthPp,
                                               1.1 * particulars.lengthPp,
                                               0.45 * particulars.lengthPp};
        const sim::Vec3 eye = ship.state.position + offset;
        const sim::Mat4 mvp =
            sim::perspective(50.0 * sim::kDegToRad, double(kWidth) / kHeight, 1.0,
                             12.0 * particulars.lengthPp) *
            sim::lookAt(eye, ship.state.position, {0, 0, 1});
        view.eye[0] = float(eye.x);
        view.eye[1] = float(eye.y);
        view.eye[2] = float(eye.z);

        float mvpf[16];
        for (int i = 0; i < 16; ++i) mvpf[i] = float(mvp.m[i]);
        const float clear[4] = {0.55f, 0.68f, 0.82f, 1.0f};

        core::Image image;
        if (!renderer.render(mvpf, view, scene, library, clear, image)) {
            std::printf("render failed on frame %d\n", frame);
            return 1;
        }
        char path[512];
        std::snprintf(path, sizeof(path), "%s/seaway_%03d.png", options.out.c_str(), frame);
        if (!core::writePng(path, image)) {
            std::printf("could not write %s\n", path);
            return 1;
        }
    }

    const sim::Mat3 R = ship.state.orientation.toMat3();
    std::printf("  %d frames, worst heel %.2f deg, worst trim %.2f deg, speed %.2f m/s\n",
                options.frames, worstHeel, worstTrim,
                dot(ship.state.velocity, R * sim::Vec3{1, 0, 0}));
    std::printf("ok\n");
    return 0;
}
