// SPDX-License-Identifier: MIT
//
// The Phase 4 picture: an engine-room fire, drawn from the gas the two-zone model
// actually solved for.
//
// Every figure `docs/03-renderer-audio.md`'s fire section publishes comes out of
// this tool, so that none of them can go stale the way a quoted number does. It
// drives `sim::fire::Model` on the reference ferry, turns each tracked gas space
// into the prism the model solves on with `gpu::volumesFromFire`, and renders the
// pair with `gpu::SmokeRenderer` -- nothing about the physics is duplicated or
// re-tuned for display.
//
// **What it does not draw**, and the reason is in `engine/gpu/smoke.hpp`: no
// plume, no flame, no ceiling jet, no turbulent detail. A two-zone model carries
// an entrainment *rate* and a mean flame height, not a shape, and inventing one
// would put structure in the picture that the simulation does not have.
//
//   ./smoke_view [--out=DIR] [--frames=N] [--duration=SECONDS] [--power=WATTS]
#include "engine/core/geometry.hpp"
#include "engine/core/math.hpp"
#include "engine/core/png.hpp"
#include "engine/gpu/hull.hpp"
#include "engine/gpu/material.hpp"
#include "engine/gpu/smoke.hpp"
#include "engine/sim/fire.hpp"
#include "engine/sim/ship.hpp"
#include "game/prototype/ferry.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

namespace {

constexpr std::uint32_t kWidth = 960;
constexpr std::uint32_t kHeight = 540;

int failures = 0;

void check(const char* what, bool ok, const std::string& detail = {}) {
    std::printf("  [%s] %-62s %s\n", ok ? "PASS" : "FAIL", what, detail.c_str());
    if (!ok) ++failures;
}

std::string number(double value, int digits = 3) {
    char buffer[64];
    std::snprintf(buffer, sizeof buffer, "%.*g", digits, value);
    return buffer;
}

sim::Mat4 camera(double angleRadians, double distance, double height, const sim::Vec3& target,
                 sim::Vec3& eyeOut) {
    eyeOut = {target.x + distance * std::cos(angleRadians),
              target.y + distance * std::sin(angleRadians), height};
    return sim::perspective(45.0 * sim::kDegToRad, static_cast<double>(kWidth) / kHeight, 1.0,
                            2000.0) *
           sim::lookAt(eyeOut, target, {0, 0, 1});
}

}  // namespace

int main(int argc, char** argv) {
    std::string outputDirectory = ".";
    int frameCount = 8;
    double duration = 600.0;
    double power = 4.0e6;

    for (int i = 1; i < argc; ++i) {
        const std::string_view argument = argv[i];
        if (argument.starts_with("--out=")) outputDirectory = std::string(argument.substr(6));
        else if (argument.starts_with("--frames=")) frameCount = std::atoi(argv[i] + 9);
        else if (argument.starts_with("--duration=")) duration = std::atof(argv[i] + 11);
        else if (argument.starts_with("--power=")) power = std::atof(argv[i] + 8);
        else { std::fprintf(stderr, "unknown argument: %s\n", argv[i]); return 2; }
    }

    std::printf("shipsim - volumetric fire and smoke from a two-zone model\n");

    gpu::Device device;
    std::string error;
    if (!device.create(error)) {
        std::printf("  no usable GPU (%s)\n", error.c_str());
        return 0;  // skip, do not fail
    }
    gpu::SmokeRenderer renderer;
    if (!renderer.create(device, kWidth, kHeight, SHIPSIM_SHADER_DIR, error)) {
        std::printf("  renderer unavailable (%s)\n", error.c_str());
        return 0;
    }
    gpu::MaterialLibrary library;
    if (!library.load(std::string(SHIPSIM_MATERIAL_DIR) + "/marine.materials", error)) {
        std::printf("  material set unavailable (%s)\n", error.c_str());
        return 0;
    }
    std::printf("  %s, %ux%u\n", device.deviceName().c_str(), kWidth, kHeight);

    // --- the casualty ---------------------------------------------------------
    sim::Ship ship = game::buildFerry();
    const sim::Sea sea;
    ship.initialise(sea);

    sim::fire::Model model;
    model.attach(ship, {ship.findCompartment("engine_room_s"),
                        ship.findCompartment("engine_room_p")});
    sim::fire::DesignFire fire;
    fire.name = "machinery";
    fire.compartment = model.findGas("engine_room_s");
    fire.baseZ = 2.5;
    fire.diameter = 2.5;
    fire.growthCoefficient = sim::fire::kGrowthFast;
    fire.peakHeatRelease = power;
    fire.steadyDuration = duration * 2.0;
    model.fires.push_back(fire);
    check("the fire model is self-consistent", model.validate().empty());

    const int burning = model.findGas("engine_room_s");
    const int neighbour = model.findGas("engine_room_p");
    check("both engine rooms are tracked", burning >= 0 && neighbour >= 0);
    if (failures != 0) return 1;

    const gpu::SmokeShading shading;

    // The scene: her **port** half, drawn as a backdrop, with the camera to
    // starboard -- so the starboard engine room's gas sits in the space the cut
    // removed and nothing is between it and the eye.
    //
    // `clipByPlane` caps its cut, so the section reads as solid steel rather than
    // as an open box; the port engine room's gas is therefore behind that cap and
    // correctly invisible, which is the depth buffer doing its job. The
    // compartment meshes are not drawn for the same reason: their own near walls
    // would stop every ray at the bulkhead.
    const sim::TriMesh cut = sim::clipByPlane(ship.hull, {0, -1, 0}, 0.0);
    gpu::SceneMesh mesh;
    const int steel = library.find("painted_steel_topside");
    check("the hull material is in the library", steel >= 0);
    mesh.appendMesh(cut, ship.state.orientation.toMat3(), ship.state.position,
                    static_cast<std::uint32_t>(std::max(steel, 0)), gpu::HullShading{});

    // Checked on the same terms as the two `findGas` lookups twenty lines above,
    // which this one sat between and did not follow.
    const int roomIndex = ship.findCompartment("engine_room_s");
    check("the burning compartment is on the ship", roomIndex >= 0);
    if (failures != 0) return 1;
    const sim::Compartment& room = ship.compartments[static_cast<std::size_t>(roomIndex)];
    const sim::Vec3 target{0.5 * (room.bboxLo.x + room.bboxHi.x), 0.0,
                           0.5 * (room.bboxLo.z + room.bboxHi.z)};

    // Overcast daylight. **Not night, and the reason is a measurement**: the
    // temperature at which a grey layer first puts a single byte on the screen is
    // printed below, and this fire's upper layer never reaches it. There is no
    // glow to make room for, so lighting the ship to see the smoke against is the
    // honest arrangement.
    gpu::SceneView view;
    view.sunDirection[0] = 0.25f;
    view.sunDirection[1] = -0.55f;
    view.sunDirection[2] = 0.80f;
    view.skyColour[0] = 0.34f;
    view.skyColour[1] = 0.39f;
    view.skyColour[2] = 0.46f;
    const float clear[4] = {0.42f, 0.48f, 0.56f, 1.0f};
    const int clearByte[3] = {static_cast<int>(std::lround(clear[0] * 255.0f)),
                              static_cast<int>(std::lround(clear[1] * 255.0f)),
                              static_cast<int>(std::lround(clear[2] * 255.0f))};

    // The temperature at which the upper layer first contributes one byte of red.
    // Bisected on `gpu::emissiveColour`, so it moves if the display mapping does
    // and is never a number written down beside one.
    double glowLo = 300.0, glowHi = 3000.0;
    for (int i = 0; i < 200; ++i) {
        const double mid = 0.5 * (glowLo + glowHi);
        double emission[3];
        gpu::emissiveColour(mid, shading, emission);
        (emission[0] * 255.0 >= 0.5 ? glowHi : glowLo) = mid;
    }
    const double glowTemperature = 0.5 * (glowLo + glowHi);

    // **A second threshold, because the check below does not test the first one.**
    // `glowTemperature` is where emission first registers *at all* -- one byte of
    // red. The pixel counter asks something much stronger: `red > blue + 8`, red
    // *dominance*, against a sky-blue clear colour whose blue byte is 143. Between
    // the two a layer is emitting and is not red-dominant, and the tool used to
    // assert a glow it had never claimed it would draw -- so it exited non-zero at
    // 10 MW (876 K, above 834, zero red-dominant pixels) for no reason but its own
    // disagreement with itself. Bisected on the same margin the counter uses, so
    // the two move together if the display mapping changes.
    double domLo = 300.0, domHi = 3000.0;
    for (int i = 0; i < 200; ++i) {
        const double mid = 0.5 * (domLo + domHi);
        double emission[3];
        gpu::emissiveColour(mid, shading, emission);
        (emission[0] * 255.0 >= emission[2] * 255.0 + 8.0 ? domHi : domLo) = mid;
    }
    const double dominantTemperature = 0.5 * (domLo + domHi);

    // --- the identity that says the picture is the model's ---------------------
    std::printf("\n1. The drawn prism is the prism the model solved on\n");
    {
        const std::vector<gpu::SmokeVolume> volumes = gpu::volumesFromFire(model, ship, shading);
        bool volumeExact = true, areaExact = true;
        for (std::size_t i = 0; i < volumes.size(); ++i) {
            const double drawn = volumes[i].planArea() * (volumes[i].hi.z - volumes[i].lo.z);
            if (std::abs(drawn - model.gas[i].gasVolume) > model.gas[i].gasVolume * 1e-12)
                volumeExact = false;
            if (std::abs(volumes[i].planArea() - model.gas[i].floorArea) >
                model.gas[i].floorArea * 1e-12)
                areaExact = false;
        }
        check("the drawn plan area is the model's floor area", areaExact);
        check("so the drawn gas volume is the model's gas volume, to machine precision",
              volumeExact);
        const double bbox = (room.bboxHi.x - room.bboxLo.x) * (room.bboxHi.y - room.bboxLo.y);
        check("and it is smaller than the bounding box, as the model says it must be",
              volumes[static_cast<std::size_t>(burning)].planArea() < bbox,
              number(100.0 * volumes[static_cast<std::size_t>(burning)].planArea() / bbox, 4) +
                  "% of the bounding box in plan");
    }

    // --- zero smoke is the unsmoked image -------------------------------------
    std::printf("\n2. An unlit fire is exactly no fire\n");
    {
        sim::Vec3 eye;
        float matrix[16];
        camera(-0.85, 58.0, 13.0, target, eye).toFloats(matrix);
        view.eye[0] = static_cast<float>(eye.x);
        view.eye[1] = static_cast<float>(eye.y);
        view.eye[2] = static_cast<float>(eye.z);

        gpu::HullRenderer hull;
        core::Image reference, empty, cold;
        bool ok = hull.create(device, kWidth, kHeight, SHIPSIM_SHADER_DIR, error) &&
                  hull.render(matrix, view, mesh, library, clear, reference) &&
                  renderer.render(matrix, view, mesh, library, {}, clear, empty);
        check("the lit pass renders both ways", ok);

        // A cold ship: every compartment at ambient with no products in it, so
        // every extinction is zero and every emission is Planck at 288 K, which is
        // nothing the display can show.
        sim::fire::Model cool;
        cool.attach(ship, {ship.findCompartment("engine_room_s"),
                           ship.findCompartment("engine_room_p")});
        const std::vector<gpu::SmokeVolume> coldVolumes =
            gpu::volumesFromFire(cool, ship, shading);
        ok = ok && renderer.render(matrix, view, mesh, library, coldVolumes, clear, cold);

        check("an empty volume list is the hull renderer's frame, bit for bit",
              ok && empty.rgba == reference.rgba);
        check("and so is a ship whose compartments hold nothing but ambient air",
              ok && cold.rgba == reference.rgba);
        // Guard against a vacuous identity: two blank frames also match.
        std::size_t drawn = 0;
        for (std::size_t i = 0; ok && i + 3 < reference.rgba.size(); i += 4)
            if (std::abs(static_cast<int>(reference.rgba[i]) - clearByte[0]) > 6 ||
                std::abs(static_cast<int>(reference.rgba[i + 1]) - clearByte[1]) > 6)
                ++drawn;
        check("the identity was checked on a frame with a ship in it",
              drawn > 20000, std::to_string(drawn) + " pixels of ship");
    }

    // --- the fire, frame by frame ---------------------------------------------
    std::printf("\n3. The casualty\n");
    std::printf("     %6s %8s %8s %8s %9s %9s %9s %8s %10s %14s\n", "frame", "t s", "Q MW",
                "T_u K", "z_i m", "k 1/m", "tau", "S m", "emit R", "file");

    double smokeMs = 0, opaqueMs = 0;
    double firstInterface = 0, lastInterface = 0, lowestInterface = 1e30;
    double firstOpacity = 0, lastOpacity = 0, peakTemperature = 0, peakEmission = 0;
    core::Image image;
    std::vector<std::size_t> litPixels, darkPixels;

    for (int frame = 0; frame < frameCount; ++frame) {
        const double targetTime = duration * frame / std::max(1, frameCount - 1);
        while (model.time < targetTime - 1e-9)
            model.step(std::min(2.0, targetTime - model.time), ship, sea);

        const std::vector<gpu::SmokeVolume> volumes = gpu::volumesFromFire(model, ship, shading);
        const gpu::SmokeVolume& hot = volumes[static_cast<std::size_t>(burning)];
        const sim::fire::GasCompartment& gas = model.gas[static_cast<std::size_t>(burning)];

        if (frame == 0) { firstInterface = hot.interfaceZ; firstOpacity = hot.upper.extinction; }
        lowestInterface = std::min(lowestInterface, hot.interfaceZ);
        lastInterface = hot.interfaceZ;
        lastOpacity = hot.upper.extinction;
        peakTemperature = std::max(peakTemperature, gas.upper.temperature());
        peakEmission = std::max(peakEmission, hot.upper.emission[0]);

        // Optical depth across the room at the deckhead, and the visibility that
        // implies: Jin's `S = K / k` with K = 3 for a reflecting sign. Reported
        // because it is the number that says whether the picture is smoke or fog.
        const double across = hot.hi.x - hot.lo.x;
        const double tau = hot.upper.extinction * across;
        const double visibility = hot.upper.extinction > 0 ? 3.0 / hot.upper.extinction : 1e30;

        sim::Vec3 eye;
        // A short sweep across the starboard quarter, level enough that the layer
        // interface reads as the horizontal plane it is.
        const double angle = -1.15 + 0.6 * frame / std::max(1, frameCount - 1);
        float matrix[16];
        camera(angle, 58.0, 13.0, target, eye).toFloats(matrix);
        view.eye[0] = static_cast<float>(eye.x);
        view.eye[1] = static_cast<float>(eye.y);
        view.eye[2] = static_cast<float>(eye.z);

        if (!renderer.render(matrix, view, mesh, library, volumes, clear, image)) {
            check("frame rendered", false);
            break;
        }
        smokeMs = std::max(smokeMs, renderer.lastFrame().smokeGpuSeconds * 1e3);
        opaqueMs = std::max(opaqueMs, renderer.lastFrame().opaqueGpuSeconds * 1e3);

        // Two counts, because the two things a two-zone layer can do to a frame
        // are opposite: emission adds red, and extinction takes everything away.
        std::size_t lit = 0, dark = 0;
        for (std::size_t i = 0; i + 3 < image.rgba.size(); i += 4) {
            if (image.rgba[i] > image.rgba[i + 2] + 8) ++lit;
            if (image.rgba[i] < 12 && image.rgba[i + 1] < 12 && image.rgba[i + 2] < 12) ++dark;
        }
        litPixels.push_back(lit);
        darkPixels.push_back(dark);

        char name[64];
        std::snprintf(name, sizeof(name), "smoke_%02d.png", frame);
        // Checked for the reason `ferry_view` now states: the figure gate reads
        // this tool's stdout and never opens the file it names.
        if (!core::writePng(outputDirectory + "/" + name, image)) {
            std::printf("could not write %s into %s\n", name, outputDirectory.c_str());
            return 1;
        }
        std::printf("     %6d %8.0f %8.2f %8.0f %9.2f %9.3f %9.1f %8s %10.2e %14s\n", frame,
                    model.time, fire.heatRelease(model.time) * 1e-6, gas.upper.temperature(),
                    hot.interfaceZ, hot.upper.extinction, tau,
                    (visibility > 1e6 ? std::string("clear") : number(visibility, 3)).c_str(),
                    hot.upper.emission[0], name);
    }

    std::printf("\n4. What it costs and what it shows\n");
    // **The layer descends and then recovers**, so the assertion is about the
    // descent and the recovery separately. It reaches its lowest point as the fire
    // stops growing and then rises a little as the room settles into its vented
    // steady state and the hot layer loses mass through the door. Asserting a
    // monotone descent would be asserting something the model does not do.
    check("the smoke layer descended", lowestInterface < firstInterface - 0.5,
          number(firstInterface, 4) + " m to " + number(lowestInterface, 4) + " m");
    check("and recovered only slightly once the room reached its steady state",
          lastInterface >= lowestInterface && lastInterface < lowestInterface + 0.5,
          "back to " + number(lastInterface, 4) + " m");
    check("the layer thickened optically", lastOpacity > firstOpacity * 10.0,
          "k " + number(firstOpacity) + " to " + number(lastOpacity) + " 1/m, visibility 3/k = " +
              number(3.0 / lastOpacity, 3) + " m");
    // **Both counts, because the layer can reach the frame either way.** Fifty lines
    // above, the counting loop says so outright -- "the two things a two-zone layer
    // can do to a frame are opposite: emission adds red, and extinction takes
    // everything away" -- and then this check looked only at extinction. A layer hot
    // enough to glow is not black, so at 12 MW it blacked out zero pixels, lit 16 309,
    // and was reported as never having reached the frame. The tool exited non-zero on
    // a picture that was right, which is the expensive kind of wrong: it makes a
    // correct render indistinguishable from a broken one.
    const auto grew = [](const std::vector<std::size_t>& v) -> std::size_t {
        return v.empty() || v.back() <= v.front() ? 0 : v.back() - v.front();
    };
    check("and the smoke reached the frame",
          !darkPixels.empty() && grew(darkPixels) + grew(litPixels) > 3000,
          std::to_string(grew(darkPixels)) + " pixels blacked out and " +
              std::to_string(grew(litPixels)) + " lit, against 3000 either way");

    // **Whether there is a glow at all is a result, not a setting.** The layer
    // emits Planck's spectrum at its own temperature, so it is a light source only
    // above the threshold printed here -- and a 4 MW machinery fire in 1150 m^3 of
    // engine room does not get there. What a two-zone fire looks like, at this
    // power, is smoke.
    std::printf("     a grey layer first puts one byte of red on the screen at %.0f K and"
                " first goes red-dominant at %.0f K;\n     this one peaked at %.0f K"
                " (emission %.2e)\n",
                glowTemperature, dominantTemperature, peakTemperature, peakEmission);
    if (peakTemperature < glowTemperature) {
        check("the layer never got hot enough to be a light source, and is not drawn as one",
              !litPixels.empty() && litPixels.back() == 0 && peakEmission * 255.0 < 0.5,
              "no red-dominant pixels, as Planck requires");
    } else if (peakTemperature < dominantTemperature) {
        // **The band the tool used to be wrong in, and it is a real state rather
        // than a gap in the argument.** The layer is emitting -- so the first check
        // would be false -- and it is not red-dominant, so the third would fail.
        // A 10 MW fire in this room lands here. What Planck requires of it is
        // exactly this: emission on the screen, no dominance.
        check("the layer emits but is not yet red-dominant, which is what Planck requires here",
              !litPixels.empty() && litPixels.back() == 0 && peakEmission * 255.0 >= 0.5,
              number(peakEmission * 255.0, 3) + " codes of red, " +
                  std::to_string(litPixels.empty() ? 0 : litPixels.back()) +
                  " red-dominant pixels; dominance needs " +
                  number(dominantTemperature, 4) + " K");
    } else {
        check("the layer got hot enough to glow, and does",
              !litPixels.empty() && litPixels.back() > litPixels.front(),
              std::to_string(litPixels.empty() ? 0 : litPixels.back()) + " red-dominant pixels");
    }

    check("the volumetric pass costs a small part of a frame",
          smokeMs > 0.0 && smokeMs < 16.0,
          "volumetric " + number(smokeMs, 3) + " ms, lit " + number(opaqueMs, 3) + " ms at " +
              std::to_string(kWidth) + "x" + std::to_string(kHeight));

    std::printf("\n   Not drawn, and named rather than hidden: the plume (a two-zone model\n"
                "   carries an entrainment rate, not a shape), the flame (a mean height and no\n"
                "   flame temperature), the ceiling jet, and any structure inside a layer --\n"
                "   a zone is well mixed by definition. See engine/gpu/smoke.hpp.\n");

    std::printf("\n%s\n", failures == 0 ? "ok" : "CHECKS FAILED");
    return failures == 0 ? 0 : 1;
}
