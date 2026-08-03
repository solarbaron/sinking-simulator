// SPDX-License-Identifier: MIT
//
// The Phase 1 milestone: the ferry casualty in 3D.
//
// Drives the same `sim::Ship` the headless scenarios use, through the same
// `core::Scheduler` the engine will use for everything else, and renders it with
// `gpu::OffscreenRenderer`. Nothing about the physics is duplicated or
// reimplemented for display -- the compartment tints are read straight off
// `Compartment::fillFraction()` as the flooding solver leaves them.
//
// The hull is drawn as a **cutaway**: the starboard half is clipped away with
// `sim::clipByPlane`, the same mesh boolean that carves the compartments in the
// first place. A solid hull would simply hide everything worth looking at, and a
// wireframe would need `fillModeNonSolid`, which is a device feature this build
// does not request.
//
//   ./ferry_view [--out=DIR] [--frames=N] [--duration=SECONDS]
#include "engine/core/geometry.hpp"
#include "engine/core/math.hpp"
#include "engine/core/png.hpp"
#include "engine/core/scheduler.hpp"
#include "engine/gpu/offscreen.hpp"
#include "game/prototype/ferry.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

using gpu::MeshVertex;

namespace {

constexpr std::uint32_t kWidth = 960;
constexpr std::uint32_t kHeight = 540;

int failures = 0;

void check(const char* what, bool ok, const std::string& detail = {}) {
    std::printf("  [%s] %-58s %s\n", ok ? "PASS" : "FAIL", what, detail.c_str());
    if (!ok) ++failures;
}

struct Geometry {
    std::vector<MeshVertex> vertices;
    std::vector<std::uint32_t> indices;

    void append(const sim::TriMesh& mesh, float r, float g, float b) {
        const auto base = static_cast<std::uint32_t>(vertices.size());
        for (const sim::Vec3& v : mesh.verts)
            vertices.push_back({{static_cast<float>(v.x), static_cast<float>(v.y),
                                 static_cast<float>(v.z)},
                                {r, g, b}});
        for (const sim::Tri& t : mesh.tris) {
            indices.push_back(base + t.a);
            indices.push_back(base + t.b);
            indices.push_back(base + t.c);
        }
    }
};

// Dry compartments read as pale steel; flooded ones as deep water. Interpolating
// on fill fraction is what makes the flooding legible without any extra state.
void floodColour(double fill, float& r, float& g, float& b) {
    const auto t = static_cast<float>(std::clamp(fill, 0.0, 1.0));
    r = 0.62f * (1.0f - t) + 0.05f * t;
    g = 0.64f * (1.0f - t) + 0.32f * t;
    b = 0.68f * (1.0f - t) + 0.85f * t;
}

// Ship body frame to world, from the rigid-body state the solver maintains.
sim::Mat4 shipModelMatrix(const sim::Ship& ship) {
    const sim::Mat3 rotation = ship.state.orientation.toMat3();
    sim::Mat4 model;
    for (int r = 0; r < 3; ++r)
        for (int c = 0; c < 3; ++c) model(r, c) = rotation(r, c);
    model(0, 3) = ship.state.position.x;
    model(1, 3) = ship.state.position.y;
    model(2, 3) = ship.state.position.z;
    return model;
}

Geometry buildScene(const sim::Ship& ship, bool cutaway) {
    Geometry geometry;
    // Keep the port half so the camera can look in from starboard. clipByPlane
    // caps the cut, so the shell reads as a solid section rather than a hole.
    const sim::TriMesh hull =
        cutaway ? sim::clipByPlane(ship.hull, {0, -1, 0}, 0.0) : ship.hull;
    geometry.append(hull, 0.30f, 0.33f, 0.38f);

    for (const sim::Compartment& compartment : ship.compartments) {
        float r = 0, g = 0, b = 0;
        floodColour(compartment.fillFraction(), r, g, b);
        geometry.append(compartment.mesh, r, g, b);
    }
    return geometry;
}

sim::Mat4 orbitCamera(double angleRadians, double distance, double height, sim::Vec3& eyeOut) {
    eyeOut = {distance * std::cos(angleRadians), distance * std::sin(angleRadians), height};
    const sim::Mat4 view = sim::lookAt(eyeOut, {0, 0, 0}, {0, 0, 1});
    const sim::Mat4 projection = sim::perspective(
        50.0 * sim::kDegToRad, static_cast<double>(kWidth) / kHeight, 1.0, 2000.0);
    return projection * view;
}

std::size_t countNonBackground(const core::Image& image, const float clear[3]) {
    const auto cr = static_cast<int>(clear[0] * 255.0f + 0.5f);
    const auto cg = static_cast<int>(clear[1] * 255.0f + 0.5f);
    const auto cb = static_cast<int>(clear[2] * 255.0f + 0.5f);
    std::size_t count = 0;
    for (std::uint32_t y = 0; y < image.height; ++y)
        for (std::uint32_t x = 0; x < image.width; ++x) {
            const std::uint8_t* p = image.pixel(x, y);
            if (std::abs(p[0] - cr) > 6 || std::abs(p[1] - cg) > 6 || std::abs(p[2] - cb) > 6)
                ++count;
        }
    return count;
}

}  // namespace

int main(int argc, char** argv) {
    std::string outputDirectory = ".";
    int frameCount = 12;
    double duration = 900.0;

    for (int i = 1; i < argc; ++i) {
        const std::string_view argument = argv[i];
        if (argument.starts_with("--out=")) outputDirectory = std::string(argument.substr(6));
        else if (argument.starts_with("--frames=")) frameCount = std::atoi(argv[i] + 9);
        else if (argument.starts_with("--duration=")) duration = std::atof(argv[i] + 11);
        else { std::fprintf(stderr, "unknown argument: %s\n", argv[i]); return 2; }
    }

    std::printf("shipsim - ferry casualty in 3D (Phase 1 milestone)\n");

    gpu::Device device;
    std::string error;
    if (!device.create(error)) {
        std::printf("  no usable GPU (%s)\n", error.c_str());
        return 0;  // skip, do not fail
    }
    gpu::OffscreenRenderer renderer;
    if (!renderer.create(device, kWidth, kHeight, SHIPSIM_SHADER_DIR, error)) {
        std::printf("  renderer unavailable (%s)\n", error.c_str());
        return 0;
    }
    std::printf("  %s, %ux%u\n", device.deviceName().c_str(), kWidth, kHeight);

    const float clear[4] = {0.05f, 0.07f, 0.10f, 1.0f};
    sim::Ship ship = game::buildFerry();
    ship.initialise(0.0);

    // --- closed-form checks --------------------------------------------------
    std::printf("\n1. Geometry and camera\n");
    {
        // Beam-on, not bow-on. Viewed from dead ahead a 120 m ship is foreshortened
        // to almost nothing and legitimately covers about 1% of the frame -- the
        // first version of this check used angle 0 and failed for that reason,
        // which was the camera choice being wrong rather than the renderer.
        sim::Vec3 eye;
        const sim::Mat4 mvp =
            orbitCamera(sim::kPi * 0.5, 110.0, 25.0, eye) * shipModelMatrix(ship);

        // The bow and stern are known points in the body frame; the camera maths
        // says exactly where they land, and the render must agree.
        double clip[4], bowX = 0, bowY = 0, sternX = 0, sternY = 0;
        mvp.transform({ship.hullHi.x, 0, 0}, clip);
        const bool bowVisible = sim::clipToPixel(clip, kWidth, kHeight, bowX, bowY);
        mvp.transform({ship.hullLo.x, 0, 0}, clip);
        const bool sternVisible = sim::clipToPixel(clip, kWidth, kHeight, sternX, sternY);
        check("bow and stern both project in front of the camera", bowVisible && sternVisible);
        check("the ship spans a useful part of the frame",
              std::abs(bowX - sternX) > kWidth * 0.4 && std::abs(bowX - sternX) < kWidth,
              std::to_string(static_cast<int>(std::abs(bowX - sternX))) + " px of " +
                  std::to_string(kWidth));

        // Projected bounding box of the hull's extreme points. A silhouette can
        // never exceed its own bounding box, and a ship profile fills most of
        // one -- so both directions are a real constraint rather than a guess.
        double boxMinX = 1e18, boxMaxX = -1e18, boxMinY = 1e18, boxMaxY = -1e18;
        for (int corner = 0; corner < 8; ++corner) {
            const sim::Vec3 point{corner & 1 ? ship.hullHi.x : ship.hullLo.x,
                                  corner & 2 ? ship.hullHi.y : ship.hullLo.y,
                                  corner & 4 ? ship.hullHi.z : ship.hullLo.z};
            double px = 0, py = 0;
            mvp.transform(point, clip);
            if (!sim::clipToPixel(clip, kWidth, kHeight, px, py)) continue;
            boxMinX = std::min(boxMinX, px); boxMaxX = std::max(boxMaxX, px);
            boxMinY = std::min(boxMinY, py); boxMaxY = std::max(boxMaxY, py);
        }
        const double boxArea = (boxMaxX - boxMinX) * (boxMaxY - boxMinY);

        Geometry scene = buildScene(ship, false);
        float matrix[16];
        mvp.toFloats(matrix);
        core::Image image;
        gpu::OffscreenRenderer::Draw draw{scene.vertices.data(), scene.vertices.size(),
                                          scene.indices.data(), scene.indices.size()};
        check("the intact hull renders", renderer.render(matrix, draw, clear, image));

        const std::size_t covered = countNonBackground(image, clear);
        check("the silhouette fits inside its projected bounding box",
              static_cast<double>(covered) <= boxArea * 1.02,
              std::to_string(covered) + " px vs box " + std::to_string(static_cast<int>(boxArea)));
        check("the silhouette fills most of that box, as a hull profile should",
              static_cast<double>(covered) > boxArea * 0.35,
              std::to_string(static_cast<int>(100.0 * covered / boxArea)) + "% of box");

        // Port and starboard views of a symmetric intact ship must mirror. A
        // handedness bug in the camera or the winding shows up here and nowhere
        // else -- a single viewpoint looks perfectly fine either way.
        sim::Vec3 portEye, starboardEye;
        const sim::Mat4 portMvp =
            orbitCamera(sim::kPi * 0.5, 110.0, 0.0, portEye) * shipModelMatrix(ship);
        const sim::Mat4 starboardMvp =
            orbitCamera(-sim::kPi * 0.5, 110.0, 0.0, starboardEye) * shipModelMatrix(ship);

        core::Image portImage, starboardImage;
        float portMatrix[16], starboardMatrix[16];
        portMvp.toFloats(portMatrix);
        starboardMvp.toFloats(starboardMatrix);
        renderer.render(portMatrix, draw, clear, portImage);
        renderer.render(starboardMatrix, draw, clear, starboardImage);

        std::size_t mismatched = 0;
        for (std::uint32_t y = 0; y < kHeight; ++y)
            for (std::uint32_t x = 0; x < kWidth; ++x) {
                const std::uint8_t* a = portImage.pixel(x, y);
                const std::uint8_t* b = starboardImage.pixel(kWidth - 1 - x, y);
                if (std::abs(a[0] - b[0]) > 8 || std::abs(a[1] - b[1]) > 8 ||
                    std::abs(a[2] - b[2]) > 8)
                    ++mismatched;
            }
        const double mismatchFraction =
            static_cast<double>(mismatched) / (double(kWidth) * kHeight);
        // Guard against a vacuous pass: two blank frames also mirror perfectly.
        check("the mirror comparison had something to compare",
              countNonBackground(portImage, clear) > 20000);
        check("port and starboard views of a symmetric hull mirror each other",
              mismatchFraction < 0.02,
              std::to_string(static_cast<int>(mismatchFraction * 1000) / 10.0) + "% differ");
    }

    // --- flooding must be visible -------------------------------------------
    std::printf("\n2. Flooding is legible\n");
    {
        sim::Vec3 eye;
        const sim::Mat4 camera = orbitCamera(-0.9, 160.0, 40.0, eye);
        float matrix[16];

        Geometry dry = buildScene(ship, true);
        (camera * shipModelMatrix(ship)).toFloats(matrix);
        core::Image dryImage;
        gpu::OffscreenRenderer::Draw dryDraw{dry.vertices.data(), dry.vertices.size(),
                                             dry.indices.data(), dry.indices.size()};
        check("the dry cutaway renders", renderer.render(matrix, dryDraw, clear, dryImage));

        // Flood every compartment and re-render from the identical camera. If the
        // tint were not actually driven by fill fraction the two frames would be
        // identical, which is the failure this catches.
        sim::Ship flooded = game::buildFerry();
        flooded.initialise(0.0);
        for (sim::Compartment& compartment : flooded.compartments)
            compartment.waterVolume = compartment.floodableVolume();

        Geometry wet = buildScene(flooded, true);
        core::Image wetImage;
        gpu::OffscreenRenderer::Draw wetDraw{wet.vertices.data(), wet.vertices.size(),
                                             wet.indices.data(), wet.indices.size()};
        renderer.render(matrix, wetDraw, clear, wetImage);

        std::size_t differing = 0;
        for (std::uint32_t y = 0; y < kHeight; ++y)
            for (std::uint32_t x = 0; x < kWidth; ++x) {
                const std::uint8_t* a = dryImage.pixel(x, y);
                const std::uint8_t* b = wetImage.pixel(x, y);
                if (std::abs(a[2] - b[2]) > 20) ++differing;
            }
        check("a fully flooded ship renders differently from a dry one", differing > 10000,
              std::to_string(differing) + " pixels changed");

        core::writePng(outputDirectory + "/ferry_dry.png", dryImage);
        core::writePng(outputDirectory + "/ferry_flooded.png", wetImage);
    }

    // --- the casualty, frame by frame ---------------------------------------
    std::printf("\n3. Rendering the casualty\n");
    core::Scheduler scheduler;
    const double dt = 0.02;
    scheduler.add({"ship", 1.0 / dt, 0.0, 1e30, 4096},
                  [&](double step) { ship.step(step, 0.0); });

    std::printf("     %6s %8s %8s %10s %12s\n", "frame", "t s", "heel", "flood t", "file");
    for (int frame = 0; frame < frameCount; ++frame) {
        const double target = duration * frame / std::max(1, frameCount - 1);
        while (scheduler.simulationTime() < target) scheduler.advance(0.5);

        Geometry scene = buildScene(ship, true);
        sim::Vec3 eye;
        // Sweep the quarter from astern to ahead on the starboard side, where the
        // cutaway faces the camera. A full orbit would spend half its frames
        // bow-on, where a 120 m ship is legitimately a sliver a few pixels wide,
        // and the other half looking at the uncut port shell.
        const double angle = -2.35 + 1.70 * frame / std::max(1, frameCount - 1);
        const sim::Mat4 mvp = orbitCamera(angle, 115.0, 32.0, eye) * shipModelMatrix(ship);
        float matrix[16];
        mvp.toFloats(matrix);

        core::Image image;
        gpu::OffscreenRenderer::Draw draw{scene.vertices.data(), scene.vertices.size(),
                                          scene.indices.data(), scene.indices.size()};
        if (!renderer.render(matrix, draw, clear, image)) {
            check("frame rendered", false);
            break;
        }

        char name[64];
        std::snprintf(name, sizeof(name), "ferry_%02d.png", frame);
        core::writePng(outputDirectory + "/" + name, image);

        const sim::Diagnostics diagnostics = ship.diagnostics(0.0);
        std::printf("     %6d %8.0f %8.2f %10.0f %12s\n", frame, scheduler.simulationTime(),
                    diagnostics.heelDeg, diagnostics.floodwaterMass / 1000.0, name);
    }
    check("the casualty produced a frame sequence", frameCount > 0);

    std::printf("\n%s\n", failures == 0 ? "all checks passed" : "CHECKS FAILED");
    return failures == 0 ? 0 : 1;
}
