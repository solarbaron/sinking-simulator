// SPDX-License-Identifier: MIT
//
// Validation of the volumetric pass.
//
// Nothing here is eyeballed. A participating medium is the single most inviting
// place in a renderer to accept "it looks like smoke", and the two-zone model is
// the single least forgiving thing to draw badly, because the numbers it produces
// are few enough that a wrong one has nowhere to hide. So every assertion below
// is a pixel value predicted before the render:
//
//   * **Beer-Lambert on a slab.** A medium bounded by two planes perpendicular to
//     x transmits `exp(-k t / |dir.x|)` along a ray, and the thickness and the
//     direction are both known before anything is drawn. The prediction uses no
//     code from `smoke.cpp` at all -- it is one line of arithmetic in this file.
//   * **Doubling the path squares the transmittance.** The same camera against a
//     slab of twice the thickness, held against `exp(-k d)` squared -- which is
//     the property that catches an integration carrying a stray term while every
//     single-thickness prediction still agrees.
//   * **The interface lands between two points `sim::clipToPixel` predicts** -- the
//     interface on the prism's near face, below which no ray reaches the hot layer
//     at all, and on its far face, above which every ray is in it for the whole
//     box. Seen through a box of finite depth those are two pixels apart at 70 m,
//     and the width between them is the box's depth over the viewing distance
//     rather than any softness in the model.
//   * **Zero smoke is the unsmoked image, bit for bit** -- against
//     `HullRenderer`'s own output, which also proves the lit pass here is that
//     pass rather than a copy of it that drifted.
//   * **The whole medium, per pixel, against `engine/gpu/smoke.cpp`**, which
//     states the same transfer integral in double from the formula rather than
//     from the shader.
//
// Skipped rather than failed when there is no GPU.
#include "engine/core/geometry.hpp"
#include "engine/core/math.hpp"
#include "engine/core/png.hpp"
#include "engine/gpu/hull.hpp"
#include "engine/gpu/material.hpp"
#include "engine/gpu/smoke.hpp"
#include "engine/sim/fire.hpp"
#include "engine/sim/ship.hpp"
#include "game/prototype/ferry.hpp"
#include "harness.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

using core::Image;
using testing::expectEqual;
using testing::expectNear;
using testing::expectTrue;

namespace {

constexpr std::uint32_t kWidth = 512;
constexpr std::uint32_t kHeight = 384;

// Chosen so each channel is an exact byte: 0.8 * 255 = 204, 0.6 * 255 = 153,
// 0.4 * 255 = 102. The background a prediction is held against is then the number
// this file wrote rather than the number the UNORM conversion happened to round
// it to.
constexpr float kClear[4] = {0.8f, 0.6f, 0.4f, 1.0f};
constexpr double kClearByte[3] = {204.0, 153.0, 102.0};

bool announced = false;

std::string label(const char* what, double value) {
    char buffer[192];
    std::snprintf(buffer, sizeof buffer, "%s %.6g", what, value);
    return buffer;
}

bool setup(gpu::Device& device, gpu::SmokeRenderer& renderer) {
    std::string error;
    if (!device.create(error)) {
        if (!announced) {
            std::printf("     no usable GPU (%s) - volumetric checks skipped\n", error.c_str());
            announced = true;
        }
        return false;
    }
    if (!renderer.create(device, kWidth, kHeight, SHIPSIM_SHADER_DIR, error)) {
        if (!announced) {
            std::printf("     smoke renderer unavailable (%s) - checks skipped\n", error.c_str());
            announced = true;
        }
        return false;
    }
    return true;
}

// A camera that can also say what ray goes through a pixel, built from first
// principles rather than by inverting the matrix -- so the ray and the matrix are
// two derivations of the same camera, and `testTheRaysAgreeWithTheProjection`
// checks they agree before anything else trusts either.
struct Camera {
    sim::Vec3 eye, forward, right, up;
    double f = 1.0, aspect = 1.0;
    sim::Mat4 mvp;

    sim::Vec3 rayFor(double px, double py) const {
        const double ndcX = px / kWidth * 2.0 - 1.0;
        const double ndcY = py / kHeight * 2.0 - 1.0;
        // Vulkan clip space: y points *down*, so a positive ndcY is a ray below
        // the axis and the world offset carries a minus sign. Getting this the
        // other way up renders a vertically mirrored medium over an unmirrored
        // ship, which is exactly the plausible-looking wrong image CLAUDE.md warns
        // about.
        return normalize(right * (ndcX * aspect / f) + up * (-ndcY / f) + forward);
    }
};

Camera cameraAt(const sim::Vec3& eye, const sim::Vec3& target, double fovYRadians, double nearPlane,
                double farPlane) {
    Camera camera;
    camera.eye = eye;
    camera.forward = normalize(target - eye);
    const sim::Vec3 side = cross(camera.forward, sim::Vec3{0, 0, 1});
    camera.right = normalize(side);
    camera.up = cross(camera.right, camera.forward);
    camera.f = 1.0 / std::tan(fovYRadians * 0.5);
    camera.aspect = static_cast<double>(kWidth) / kHeight;
    camera.mvp = sim::perspective(fovYRadians, camera.aspect, nearPlane, farPlane) *
                 sim::lookAt(eye, target, {0, 0, 1});
    return camera;
}

gpu::SceneView viewFor(const Camera& camera) {
    gpu::SceneView view;
    view.eye[0] = static_cast<float>(camera.eye.x);
    view.eye[1] = static_cast<float>(camera.eye.y);
    view.eye[2] = static_cast<float>(camera.eye.z);
    return view;
}

// A slab of gas bounded by two planes perpendicular to x, wide enough in y and z
// that every ray in the frame enters and leaves through those two planes. The
// path length is then `thickness / |dir.x|` with no geometry code in it.
gpu::SmokeVolume slab(double halfThickness, double extinction) {
    gpu::SmokeVolume v;
    v.name = "slab";
    v.lo = {-halfThickness, -4000.0, -4000.0};
    v.hi = {halfThickness, 4000.0, 4000.0};
    // The whole slab is the upper layer: the interface at the floor puts every
    // point of it at or above the interface.
    v.interfaceZ = -4000.0;
    v.upper.extinction = extinction;
    return v;
}

int predictedByte(double background, double transmittance, double source) {
    const double value = source + transmittance * background / 255.0;
    return static_cast<int>(std::lround(std::clamp(value, 0.0, 1.0) * 255.0));
}

bool renderMedium(gpu::SmokeRenderer& renderer, const Camera& camera,
                  const std::vector<gpu::SmokeVolume>& volumes, Image& out) {
    const gpu::SceneMesh empty;
    const gpu::MaterialLibrary library;
    float matrix[16];
    camera.mvp.toFloats(matrix);
    return renderer.render(matrix, viewFor(camera), empty, library, volumes, kClear, out);
}

// --- the ray construction, checked against the matrix it will be compared with -

void testTheRaysAgreeWithTheProjection() {
    const Camera camera = cameraAt({70.0, 6.0, 9.0}, {0.0, 0.0, 5.0},
                                   42.0 * sim::kDegToRad, 0.5, 500.0);
    int checked = 0;
    for (std::uint32_t y = 24; y < kHeight; y += 61)
        for (std::uint32_t x = 17; x < kWidth; x += 79) {
            const sim::Vec3 direction = camera.rayFor(x + 0.5, y + 0.5);
            const sim::Vec3 point = camera.eye + direction * 37.0;
            double clip[4], px = 0, py = 0;
            camera.mvp.transform(point, clip);
            expectTrue("a point on the ray projects", sim::clipToPixel(clip, kWidth, kHeight, px, py));
            expectNear("and lands on the pixel the ray came from, in x", px, x + 0.5, 1e-6);
            expectNear("and in y", py, y + 0.5, 1e-6);
            ++checked;
        }
    expectTrue("the ray check swept a useful number of pixels", checked >= 30);
}

// --- Beer-Lambert --------------------------------------------------------------

void testAUniformSlabTransmitsExpMinusKD() {
    gpu::Device device;
    gpu::SmokeRenderer renderer;
    if (!setup(device, renderer)) return;

    // Looking along -x from well outside the slab, tilted enough that |dir.x|
    // varies across the frame and the path length is a range rather than a number.
    const Camera camera = cameraAt({80.0, 0.0, 0.0}, {0.0, 0.0, 0.0},
                                   50.0 * sim::kDegToRad, 0.5, 500.0);
    const double halfThickness = 5.0;
    const double k = 0.15;

    Image thin, thick;
    expectTrue("the thin slab renders",
               renderMedium(renderer, camera, {slab(halfThickness, k)}, thin));
    expectTrue("the thick slab renders",
               renderMedium(renderer, camera, {slab(2.0 * halfThickness, k)}, thick));

    int worst = 0, worstSquare = 0, samples = 0;
    double minTransmittance = 1e30, maxTransmittance = -1e30;
    for (std::uint32_t y = 8; y < kHeight - 8; y += 13)
        for (std::uint32_t x = 8; x < kWidth - 8; x += 11) {
            const sim::Vec3 direction = camera.rayFor(x + 0.5, y + 0.5);
            // The closed form, in one line and with no geometry code behind it.
            const double path = 2.0 * halfThickness / std::abs(direction.x);
            const double transmitted = std::exp(-k * path);
            minTransmittance = std::min(minTransmittance, transmitted);
            maxTransmittance = std::max(maxTransmittance, transmitted);

            const std::uint8_t* thinPixel = thin.pixel(x, y);
            const std::uint8_t* thickPixel = thick.pixel(x, y);
            for (int c = 0; c < 3; ++c) {
                worst = std::max(worst, std::abs(static_cast<int>(thinPixel[c]) -
                                                 predictedByte(kClearByte[c], transmitted, 0.0)));
                // Twice the thickness, so the square of the transmittance: the
                // closed form again, and the property that catches an integration
                // with a stray term in it.
                worstSquare = std::max(
                    worstSquare, std::abs(static_cast<int>(thickPixel[c]) -
                                          predictedByte(kClearByte[c], transmitted * transmitted,
                                                        0.0)));
            }
            ++samples;
        }

    std::printf("     Beer-Lambert over %d pixels: transmittance %.4f..%.4f, worst error %d LSB"
                " (doubled path %d LSB)\n",
                samples, minTransmittance, maxTransmittance, worst, worstSquare);
    // One least-significant bit is the whole budget: the shader works in float,
    // the blend is exact, and the store is a round. Anything looser would pass on
    // a medium that had lost a per cent of its optical depth.
    expectTrue(label("a uniform slab transmits exp(-k d) to, LSB", worst), worst <= 1);
    expectTrue(label("and doubling the path squares it to, LSB", worstSquare), worstSquare <= 1);
    // Vacuity: a blank frame, or a medium that did nothing, would pass a pair of
    // predictions that were also blank.
    // A medium that ignored the path length would give one transmittance
    // everywhere, and both predictions would then be one number. The ratio is what
    // says the sweep asked the question.
    expectTrue(label("the sweep covered a real range of transmittance, max over min",
                     maxTransmittance / minTransmittance),
               maxTransmittance / minTransmittance > 1.4 && minTransmittance < 0.25);
    expectTrue("the sweep was not one pixel", samples > 1000);
}

// --- bit identity ---------------------------------------------------------------

void testZeroSmokeIsTheUnsmokedImage() {
    gpu::Device device;
    gpu::SmokeRenderer renderer;
    if (!setup(device, renderer)) return;
    gpu::HullRenderer hull;
    std::string error;
    if (!hull.create(device, kWidth, kHeight, SHIPSIM_SHADER_DIR, error)) {
        std::printf("     hull renderer unavailable (%s) - identity check skipped\n",
                    error.c_str());
        return;
    }
    gpu::MaterialLibrary library;
    expectTrue("the shipped material set loads",
               library.load(std::string(SHIPSIM_MATERIAL_DIR) + "/marine.materials", error));

    sim::Ship ship = game::buildFerry();
    ship.initialise(sim::Sea{});
    gpu::SceneMesh mesh;
    gpu::HullPaint paint;
    const gpu::HullShading shading;
    expectTrue("the ferry meshes", mesh.appendShip(ship, paint, library, shading, error));

    const Camera camera = cameraAt({40.0, -110.0, 35.0}, {0.0, 0.0, 6.0},
                                   45.0 * sim::kDegToRad, 1.0, 600.0);
    float matrix[16];
    camera.mvp.toFloats(matrix);
    const gpu::SceneView view = viewFor(camera);

    Image reference, empty, transparent;
    expectTrue("the hull renderer draws her",
               hull.render(matrix, view, mesh, library, kClear, reference));
    expectTrue("and so does the smoke renderer with nothing to composite",
               renderer.render(matrix, view, mesh, library, {}, kClear, empty));

    // A medium that is present, drawn, rasterised and blended -- and completely
    // clear. This is the stronger of the two: the pass runs, and the arithmetic
    // still has to come out as the identity.
    gpu::SmokeVolume clear = slab(20.0, 0.0);
    clear.upper.emission[0] = clear.upper.emission[1] = clear.upper.emission[2] = 4.0;
    clear.lower.emission[0] = 4.0;
    expectTrue("and with a clear but glowing volume over her",
               renderer.render(matrix, view, mesh, library, {clear}, kClear, transparent));

    expectTrue("an empty volume list is the hull renderer's frame, bit for bit",
               empty.rgba == reference.rgba);
    expectTrue("and so is a volume of zero extinction, however hot",
               transparent.rgba == reference.rgba);

    // Vacuity: two blank frames are also bit-identical. CLAUDE.md records a mirror
    // test in this repository that would have passed on exactly that.
    std::size_t drawn = 0;
    for (std::uint32_t y = 0; y < kHeight; ++y)
        for (std::uint32_t x = 0; x < kWidth; ++x) {
            const std::uint8_t* p = reference.pixel(x, y);
            if (std::abs(p[0] - 204) > 4 || std::abs(p[1] - 153) > 4 || std::abs(p[2] - 102) > 4)
                ++drawn;
        }
    expectTrue(label("the identity was checked on a frame with a ship in it, pixels",
                     static_cast<double>(drawn)),
               drawn > 20000);

    // And the same medium with extinction in it must *not* be the same frame, or
    // the two above prove nothing about the pass having run.
    gpu::SmokeVolume opaque = slab(20.0, 0.05);
    Image smoky;
    expectTrue("a medium with extinction renders",
               renderer.render(matrix, view, mesh, library, {opaque}, kClear, smoky));
    expectTrue("and changes the frame", smoky.rgba != reference.rgba);
}

// --- the interface ---------------------------------------------------------------

// The layer interface is a plane of no thickness at a height the model reports,
// so where it lands in the image is arithmetic rather than opinion.
//
// **It is not a single row, and pretending it were would be the wrong test.** Seen
// through a box of finite depth, the boundary between "this ray sees hot gas" and
// "this ray does not" is bracketed by two projected points: the interface on the
// box's *near* face, below which no ray reaches the hot layer at all, and the
// interface on its *far* face, above which every ray is in the hot layer for the
// whole box. Between them is a ramp a few pixels wide, and its width is the box's
// depth over the viewing distance -- geometry, not softness. So both brackets are
// projected with `sim::clipToPixel` and the frame is held against both.
void testTheInterfaceLandsWhereTheCameraSaysItDoes() {
    gpu::Device device;
    gpu::SmokeRenderer renderer;
    if (!setup(device, renderer)) return;

    const double interfaceZ = 6.5;
    gpu::SmokeVolume box;
    box.name = "compartment";
    box.lo = {-6.0, -60.0, 1.0};   // wide in y so every sampled ray leaves by the far face
    box.hi = {6.0, 60.0, 40.0};
    box.interfaceZ = interfaceZ;
    box.upper.extinction = 0.30;   // exp(-0.3 * 12) = 0.027: opaque across the box
    box.lower.extinction = 0.0;    // exactly clear, so "below" is exactly the clear colour

    // Level, and above the interface, so every ray that reaches the hot layer does
    // so through the near face. Level also means the ray's z slope depends on the
    // row alone and not on the column -- a horizontal plane projects to a
    // horizontal line -- which is what lets one pair of rows be checked against
    // every column below.
    const Camera camera = cameraAt({70.0, 0.0, 8.0}, {0.0, 0.0, 8.0},
                                   40.0 * sim::kDegToRad, 0.5, 500.0);

    Image image;
    expectTrue("the stratified box renders", renderMedium(renderer, camera, {box}, image));

    double clip[4], nearX = 0, nearY = 0, farX = 0, farY = 0;
    camera.mvp.transform({box.hi.x, 0.0, interfaceZ}, clip);
    expectTrue("the interface on the near face projects into the frame",
               sim::clipToPixel(clip, kWidth, kHeight, nearX, nearY));
    camera.mvp.transform({box.lo.x, 0.0, interfaceZ}, clip);
    expectTrue("and so does the interface on the far face",
               sim::clipToPixel(clip, kWidth, kHeight, farX, farY));
    const int nearRow = static_cast<int>(std::floor(nearY));
    const int farRow = static_cast<int>(std::floor(farY));

    std::printf("     interface: near face row %d, far face row %d (%d px of ramp)\n", nearRow,
                farRow, nearRow - farRow);
    // Rows increase downward in Vulkan, and the near face's interface is the lower
    // of the two: a ray that just grazes it is steeper than one that reaches the
    // interface only at the back of the box.
    expectTrue("the near face's interface is below the far face's", nearRow > farRow);
    expectTrue(label("and the interface is sharp: the ramp between them is, px",
                     nearRow - farRow),
               nearRow - farRow <= 4);

    int columns = 0;
    bool clearBelow = true, opaqueAbove = true, monotone = true;
    int brightestAbove = 0, darkestBelow = 255;
    for (std::uint32_t x = 200; x < 320; x += 3) {
        const std::uint8_t* below = image.pixel(x, static_cast<std::uint32_t>(nearRow + 1));
        const std::uint8_t* above = image.pixel(x, static_cast<std::uint32_t>(farRow - 1));
        if (!(below[0] == 204 && below[1] == 153 && below[2] == 102)) clearBelow = false;
        darkestBelow = std::min(darkestBelow, static_cast<int>(below[0]));
        brightestAbove = std::max(brightestAbove, static_cast<int>(above[0]));
        if (above[0] >= 30) opaqueAbove = false;
        // Between the brackets the attenuation only ever grows going up, because
        // the path in the hot layer only ever grows. A medium that had lost the
        // interface would be flat here instead.
        for (int y = nearRow + 1; y > farRow - 1; --y) {
            const int here = image.pixel(x, static_cast<std::uint32_t>(y))[0];
            const int up = image.pixel(x, static_cast<std::uint32_t>(y - 1))[0];
            if (up > here) monotone = false;
        }
        ++columns;
    }

    expectTrue("the scan covered a band of columns", columns > 30);
    expectEqual(label("below the near face's interface every column is exactly the clear"
                      " colour, darkest", darkestBelow),
                darkestBelow, 204);
    expectTrue("and it is exactly the clear colour in all three channels", clearBelow);
    expectTrue(label("above the far face's interface the hot layer is opaque, brightest",
                     brightestAbove),
               opaqueAbove);
    expectTrue("and the attenuation grows monotonically between the brackets", monotone);

    // Vacuity: the frame has to contain both regimes in quantity, or a wholly
    // opaque or wholly clear image would satisfy every check above.
    std::size_t smoky = 0, untouched = 0;
    for (std::uint32_t y = 0; y < kHeight; ++y)
        for (std::uint32_t x = 0; x < kWidth; ++x) {
            const std::uint8_t* p = image.pixel(x, y);
            if (p[0] == 204 && p[1] == 153 && p[2] == 102) ++untouched;
            else ++smoky;
        }
    expectTrue(label("the frame contained smoke, pixels", static_cast<double>(smoky)),
               smoky > 20000);
    expectTrue(label("and clear air, pixels", static_cast<double>(untouched)),
               untouched > 20000);
}

// --- against the model, per pixel ---------------------------------------------

// The whole medium, held against `engine/gpu/smoke.cpp` -- which states the same
// transfer integral in double, from the formula rather than from the shader. Two
// layers with different extinctions and different emissions, on a heeled ship so
// the body-frame transform is exercised rather than assumed to be the identity,
// and from two viewpoints so that both orders of meeting the layers are in the
// sweep.
void testEveryPixelAgreesWithTheModel() {
    gpu::Device device;
    gpu::SmokeRenderer renderer;
    if (!setup(device, renderer)) return;

    gpu::SmokeVolume box;
    box.name = "compartment";
    box.lo = {-9.0, -7.0, 0.0};
    box.hi = {9.0, 7.0, 9.0};
    box.interfaceZ = 5.0;
    box.upper.extinction = 0.22;
    box.upper.emission[0] = 0.55;
    box.upper.emission[1] = 0.18;
    box.upper.emission[2] = 0.03;
    box.lower.extinction = 0.045;
    box.lower.emission[0] = 0.02;
    box.rotation = sim::Quat::fromAxisAngle({1, 0, 0}, 18.0 * sim::kDegToRad).toMat3();
    box.translation = {1.5, -2.0, 0.75};

    // **Two cameras, one looking up through the medium and one looking down.**
    // The order the ray meets the layers in is the sign of its body-frame z
    // component, and a sweep that only ever descends -- or, as the first version of
    // this test did, only ever ascends once the 18 degrees of heel are taken out --
    // exercises one branch of that and leaves the other free. Mutation testing
    // found it: forcing `upperFirst` false unconditionally left the whole suite
    // green. The counts below are what stops that happening again.
    const Camera cameras[2] = {
        // Low and close, so a great many rays graze the interface rather than
        // plunging through it. Ascending in the body frame.
        cameraAt({34.0, 24.0, 8.5}, {1.0, -1.0, 5.0}, 46.0 * sim::kDegToRad, 0.5, 400.0),
        // High and steep, so the rays descend through the interface instead.
        cameraAt({20.0, 15.0, 42.0}, {1.0, -1.0, 4.0}, 46.0 * sim::kDegToRad, 0.5, 400.0)};

    int worst = 0, hits = 0, bothLayers = 0, hotFirst = 0, coolFirst = 0;
    double minValue = 1e30, maxValue = -1e30;
    for (const Camera& camera : cameras) {
        Image image;
        expectTrue("the two-layer medium renders", renderMedium(renderer, camera, {box}, image));

        gpu::DepthBasis basis;
        float matrix[16];
        camera.mvp.toFloats(matrix);
        expectTrue("the camera yields a depth basis", gpu::depthBasisFrom(matrix, basis));
        // Nothing opaque is drawn, so every ray runs to the far plane -- which is
        // the distance a cleared depth of 1 inverts to.
        const double farPlane = basis.b / (1.0 + basis.a);

        for (std::uint32_t y = 2; y < kHeight - 2; y += 2)
            for (std::uint32_t x = 2; x < kWidth - 2; x += 2) {
                const sim::Vec3 direction = camera.rayFor(x + 0.5, y + 0.5);
                double tEnter = 0, tExit = 0;
                if (!gpu::intersectVolume(box, camera.eye, direction, farPlane, tEnter, tExit))
                    continue;
                const gpu::LayerPath path =
                    gpu::layerPath(box, camera.eye, direction, tEnter, tExit);
                if (path.upper > 1e-3 && path.lower > 1e-3) {
                    ++bothLayers;
                    ++(path.upperFirst ? hotFirst : coolFirst);
                }

                const double background[3] = {kClearByte[0] / 255.0, kClearByte[1] / 255.0,
                                              kClearByte[2] / 255.0};
                double predicted[3];
                gpu::compositeOver(box, path, background, predicted);

                const std::uint8_t* pixel = image.pixel(x, y);
                for (int c = 0; c < 3; ++c) {
                    const int want =
                        static_cast<int>(std::lround(std::clamp(predicted[c], 0.0, 1.0) * 255.0));
                    worst = std::max(worst, std::abs(static_cast<int>(pixel[c]) - want));
                    minValue = std::min(minValue, predicted[c]);
                    maxValue = std::max(maxValue, predicted[c]);
                }
                ++hits;
            }
    }

    std::printf("     %d pixels inside the medium, %d with both layers on the ray"
                " (%d hot first, %d cool first); worst disagreement with the model %d LSB\n",
                hits, bothLayers, hotFirst, coolFirst, worst);
    expectTrue(label("the shader agrees with engine/gpu/smoke.cpp to, LSB", worst), worst <= 1);
    expectTrue("the sweep found the medium at all", hits > 2000);
    // Vacuity: a ray that only ever saw one layer would not test the split, and a
    // flat prediction would not test anything.
    expectTrue(label("and rays that crossed the interface, count",
                     static_cast<double>(bothLayers)),
               bothLayers > 800 && bothLayers > hits / 20);
    // And both orders have to be in the sweep, or half the composite is untested
    // while every pixel of it agrees.
    expectTrue(label("rays that met the hot layer first, count", static_cast<double>(hotFirst)),
               hotFirst > 400);
    expectTrue(label("rays that met the cool layer first, count", static_cast<double>(coolFirst)),
               coolFirst > 400);
    expectTrue(label("the predicted values spanned a real range, min", minValue),
               maxValue - minValue > 0.2);
}

// --- occlusion -----------------------------------------------------------------

// The ray stops at the hull. A depth *test* would drop the whole fragment where a
// bulkhead is nearer than the far face of the box; what has to happen is that the
// integral stops at the bulkhead, which is the reason the depth buffer is sampled
// rather than tested against.
void testTheMediumStopsAtWhateverIsSolid() {
    gpu::Device device;
    gpu::SmokeRenderer renderer;
    if (!setup(device, renderer)) return;
    std::string error;
    gpu::MaterialLibrary library;
    expectTrue("the shipped material set loads",
               library.load(std::string(SHIPSIM_MATERIAL_DIR) + "/marine.materials", error));
    const int material = library.find("painted_steel_topside");
    expectTrue("the topside material is in it", material >= 0);

    // A plate at x = 0, inside a slab that spans x in [-5, 5]: half the medium is
    // in front of it and half behind.
    sim::TriMesh plate;
    plate.verts = {{0.0, -20.0, -20.0}, {0.0, 20.0, -20.0}, {0.0, 20.0, 20.0}, {0.0, -20.0, 20.0}};
    plate.tris = {{0, 1, 2}, {0, 2, 3}};
    gpu::SceneMesh mesh;
    mesh.appendMesh(plate, sim::Mat3::identity(), {0, 0, 0}, static_cast<std::uint32_t>(material),
                    gpu::HullShading{});

    const Camera camera = cameraAt({80.0, 0.0, 0.0}, {0.0, 0.0, 0.0},
                                   50.0 * sim::kDegToRad, 0.5, 500.0);
    float matrix[16];
    camera.mvp.toFloats(matrix);
    const gpu::SceneView view = viewFor(camera);

    const double k = 0.25;
    Image plain, smoky;
    expectTrue("the plate renders on its own",
               renderer.render(matrix, view, mesh, library, {slab(5.0, 0.0)}, kClear, plain));
    expectTrue("and with a medium over it",
               renderer.render(matrix, view, mesh, library, {slab(5.0, k)}, kClear, smoky));

    int worst = 0, samples = 0;
    double minPath = 1e30, maxPath = -1e30;
    for (std::uint32_t y = 100; y < kHeight - 100; y += 5)
        for (std::uint32_t x = 100; x < kWidth - 100; x += 5) {
            const std::uint8_t* before = plain.pixel(x, y);
            // Only where the plate was actually drawn; the shaded value is taken
            // from the frame rather than predicted from the BRDF, because what is
            // under test here is the integration limit and not the shading.
            if (before[0] == 204 && before[1] == 153 && before[2] == 102) continue;
            const sim::Vec3 direction = camera.rayFor(x + 0.5, y + 0.5);
            // Closed form again: the medium in front of the plate runs from x = 5
            // to x = 0, so the path is 5 / |dir.x|. Nothing behind the plate may
            // contribute at all.
            const double path = 5.0 / std::abs(direction.x);
            minPath = std::min(minPath, path);
            maxPath = std::max(maxPath, path);
            const double transmitted = std::exp(-k * path);
            const std::uint8_t* after = smoky.pixel(x, y);
            for (int c = 0; c < 3; ++c)
                worst = std::max(worst, std::abs(static_cast<int>(after[c]) -
                                                 predictedByte(before[c], transmitted, 0.0)));
            ++samples;
        }

    std::printf("     occlusion: %d plate pixels, path in front of it %.3f..%.3f m,"
                " worst error %d LSB\n", samples, minPath, maxPath, worst);
    expectTrue("the plate covered a useful part of the frame", samples > 500);
    expectTrue(label("the medium stops at the plate to, LSB", worst), worst <= 1);
    // Vacuity: if the medium had contributed nothing at all, `worst` would still be
    // small only if the prediction were also nothing.
    expectTrue(label("the medium in front of the plate was optically thick enough to see,"
                     " transmittance", std::exp(-k * maxPath)),
               std::exp(-k * maxPath) < 0.35);
}

// --- ordering -------------------------------------------------------------------

// Two separated media compose in the order the ray meets them, and emission is
// what makes that matter: the nearer one's glow is attenuated by less. A renderer
// that drew them in list order would agree with this on one ordering and not the
// other, which is why both are rendered.
void testTheSortPutsTheNearerVolumeInFront() {
    gpu::Device device;
    gpu::SmokeRenderer renderer;
    if (!setup(device, renderer)) return;

    gpu::SmokeVolume nearVolume = slab(4.0, 0.20);
    nearVolume.lo.x = 16.0;
    nearVolume.hi.x = 24.0;
    nearVolume.upper.emission[0] = 0.60;
    gpu::SmokeVolume farVolume = slab(4.0, 0.20);
    farVolume.lo.x = -24.0;
    farVolume.hi.x = -16.0;
    farVolume.upper.emission[1] = 0.60;

    const Camera camera = cameraAt({80.0, 0.0, 0.0}, {0.0, 0.0, 0.0},
                                   30.0 * sim::kDegToRad, 0.5, 500.0);
    Image listOrder, reversed;
    expectTrue("near first renders", renderMedium(renderer, camera, {nearVolume, farVolume}, listOrder));
    expectTrue("far first renders", renderMedium(renderer, camera, {farVolume, nearVolume}, reversed));
    // The sort is inside the renderer, so the order the caller happened to use
    // cannot change the frame.
    expectTrue("the caller's list order does not change the frame",
               listOrder.rgba == reversed.rgba);

    gpu::DepthBasis basis;
    float matrix[16];
    camera.mvp.toFloats(matrix);
    expectTrue("the camera yields a depth basis", gpu::depthBasisFrom(matrix, basis));
    const double farPlane = basis.b / (1.0 + basis.a);

    int worst = 0, samples = 0;
    for (std::uint32_t y = 120; y < kHeight - 120; y += 5)
        for (std::uint32_t x = 160; x < kWidth - 160; x += 5) {
            const sim::Vec3 direction = camera.rayFor(x + 0.5, y + 0.5);
            double predicted[3] = {kClearByte[0] / 255.0, kClearByte[1] / 255.0,
                                   kClearByte[2] / 255.0};
            // Composited back to front by hand: the far one over the background,
            // then the near one over that.
            for (const gpu::SmokeVolume* v : {&farVolume, &nearVolume}) {
                double tEnter = 0, tExit = 0;
                if (!gpu::intersectVolume(*v, camera.eye, direction, farPlane, tEnter, tExit))
                    continue;
                double next[3];
                gpu::compositeOver(*v, gpu::layerPath(*v, camera.eye, direction, tEnter, tExit),
                                   predicted, next);
                for (int c = 0; c < 3; ++c) predicted[c] = next[c];
            }
            const std::uint8_t* pixel = listOrder.pixel(x, y);
            for (int c = 0; c < 3; ++c)
                worst = std::max(worst,
                                 std::abs(static_cast<int>(pixel[c]) -
                                          static_cast<int>(std::lround(
                                              std::clamp(predicted[c], 0.0, 1.0) * 255.0))));
            ++samples;
        }
    expectTrue("the ordering sweep had pixels", samples > 500);
    expectTrue(label("two media compose back to front to, LSB", worst), worst <= 1);

    // And the wrong order really would be a different frame: swap which one glows
    // and the two are no longer interchangeable.
    gpu::SmokeVolume swappedNear = nearVolume, swappedFar = farVolume;
    std::swap(swappedNear.lo.x, swappedFar.lo.x);
    std::swap(swappedNear.hi.x, swappedFar.hi.x);
    Image swapped;
    expectTrue("the swapped pair renders",
               renderMedium(renderer, camera, {swappedNear, swappedFar}, swapped));
    expectTrue("and moving the glowing medium to the back changes the frame",
               swapped.rgba != listOrder.rgba);
}

// --- inside the compartment ------------------------------------------------------

// The camera in the smoke, which is where a damage-control party is. The near
// faces of the box are behind the eye, so this only works because the far faces
// are the ones kept -- and the ray still has to start at the eye rather than at
// the box.
void testTheCameraCanBeInsideTheMedium() {
    gpu::Device device;
    gpu::SmokeRenderer renderer;
    if (!setup(device, renderer)) return;

    gpu::SmokeVolume box;
    box.lo = {-12.0, -9.0, 0.0};
    box.hi = {12.0, 9.0, 8.0};
    box.interfaceZ = 4.0;
    box.upper.extinction = 0.35;
    box.upper.emission[0] = 0.4;
    box.lower.extinction = 0.05;

    const Camera camera = cameraAt({6.0, 0.0, 2.0}, {-12.0, 0.0, 2.5},
                                   60.0 * sim::kDegToRad, 0.2, 300.0);
    Image image;
    expectTrue("the medium renders from inside it", renderMedium(renderer, camera, {box}, image));

    gpu::DepthBasis basis;
    float matrix[16];
    camera.mvp.toFloats(matrix);
    expectTrue("the camera yields a depth basis", gpu::depthBasisFrom(matrix, basis));
    const double farPlane = basis.b / (1.0 + basis.a);

    int worst = 0, samples = 0, insideStart = 0;
    for (std::uint32_t y = 6; y < kHeight - 6; y += 11)
        for (std::uint32_t x = 6; x < kWidth - 6; x += 9) {
            const sim::Vec3 direction = camera.rayFor(x + 0.5, y + 0.5);
            double tEnter = 0, tExit = 0;
            if (!gpu::intersectVolume(box, camera.eye, direction, farPlane, tEnter, tExit))
                continue;
            if (tEnter <= 1e-9) ++insideStart;
            const double background[3] = {kClearByte[0] / 255.0, kClearByte[1] / 255.0,
                                          kClearByte[2] / 255.0};
            double predicted[3];
            gpu::compositeOver(box, gpu::layerPath(box, camera.eye, direction, tEnter, tExit),
                               background, predicted);
            const std::uint8_t* pixel = image.pixel(x, y);
            for (int c = 0; c < 3; ++c)
                worst = std::max(worst,
                                 std::abs(static_cast<int>(pixel[c]) -
                                          static_cast<int>(std::lround(
                                              std::clamp(predicted[c], 0.0, 1.0) * 255.0))));
            ++samples;
        }
    expectTrue("the inside sweep had pixels", samples > 500);
    // Every ray must start at the eye, not at a box face -- otherwise the camera
    // is not inside anything and this test is the outside one again.
    expectTrue(label("and every one of them started at the eye, count",
                     static_cast<double>(insideStart)),
               insideStart == samples);
    expectTrue(label("the medium is right from inside it to, LSB", worst), worst <= 1);
}

// --- the ferry -------------------------------------------------------------------

// The Phase 4 milestone's picture, and what it costs. Not a pixel assertion: the
// numbers above are the assertions, and this is the frame a person can look at
// plus the price of drawing it on the target card.
void testTheFerryFireCosts() {
    gpu::Device device;
    gpu::SmokeRenderer renderer;
    if (!setup(device, renderer)) return;
    std::string error;
    gpu::MaterialLibrary library;
    if (!library.load(std::string(SHIPSIM_MATERIAL_DIR) + "/marine.materials", error)) return;

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
    fire.peakHeatRelease = 4.0e6;
    fire.steadyDuration = 900.0;
    model.fires.push_back(fire);
    for (int i = 0; i < 420; ++i) model.step(1.0, ship, sea);

    gpu::SmokeShading shading;
    const std::vector<gpu::SmokeVolume> volumes = gpu::volumesFromFire(model, ship, shading);
    expectTrue("the ferry fire produced volumes to draw", volumes.size() == 2);

    // Cut her open so the gas is visible, exactly as ferry_view does.
    const sim::TriMesh cut = sim::clipByPlane(ship.hull, {0, -1, 0}, 0.0);
    gpu::SceneMesh mesh;
    mesh.appendMesh(cut, ship.state.orientation.toMat3(), ship.state.position,
                    static_cast<std::uint32_t>(std::max(library.find("bare_steel"), 0)),
                    gpu::HullShading{});

    const Camera camera = cameraAt({55.0, -95.0, 26.0}, {10.0, 0.0, 6.0},
                                   45.0 * sim::kDegToRad, 1.0, 600.0);
    float matrix[16];
    camera.mvp.toFloats(matrix);
    gpu::SceneView view = viewFor(camera);
    view.skyColour[0] = 0.08f;
    view.skyColour[1] = 0.09f;
    view.skyColour[2] = 0.11f;
    const float night[4] = {0.02f, 0.025f, 0.035f, 1.0f};

    gpu::SmokeFrameCost best;
    double bestTotal = 1e30;
    Image image;
    for (int repeat = 0; repeat < 8; ++repeat) {
        if (!renderer.render(matrix, view, mesh, library, volumes, night, image)) return;
        const double total = renderer.lastFrame().opaqueGpuSeconds +
                             renderer.lastFrame().smokeGpuSeconds;
        if (total < bestTotal) {
            bestTotal = total;
            best = renderer.lastFrame();
        }
    }
    std::printf("     %s: %zu verts, %zu tris + %zu volumes | lit %.3f ms  volumetric %.3f ms"
                "  upload %.2f  submit %.2f  readback %.2f  total %.2f ms\n",
                device.deviceName().c_str(), best.vertices, best.triangles, best.volumes,
                best.opaqueGpuSeconds * 1e3, best.smokeGpuSeconds * 1e3,
                best.uploadSeconds * 1e3, best.submitSeconds * 1e3,
                best.readbackSeconds * 1e3, best.totalSeconds * 1e3);
    expectTrue("both passes were timed",
               best.opaqueGpuSeconds > 0.0 && best.smokeGpuSeconds > 0.0);
    // Every field of the cost record is read here, on purpose. A published timing
    // struct with a field nobody prints is a field nobody would notice going wrong
    // -- the same shape as the species loading `fire.hpp` deleted for being written
    // on every band and read by nothing.
    expectTrue("and the wall-clock account covers the passes it contains",
               best.totalSeconds >= best.uploadSeconds + best.submitSeconds +
                                        best.readbackSeconds &&
                   best.submitSeconds > 0.0 && best.readbackSeconds > 0.0 &&
                   best.uploadSeconds > 0.0);
    // The design statement: an analytic two-slab medium is not what a frame costs.
    // If this ever stops holding, the answer is not to march fewer samples -- there
    // are no samples -- but to look at what changed.
    expectTrue(label("the volumetric pass costs well under a 60 Hz frame, ms",
                     best.smokeGpuSeconds * 1e3),
               best.smokeGpuSeconds < 0.016);

    const std::string path = testing::scratchDir() + "smoke_engine_room.png";
    expectTrue("a representative frame is written: " + path, core::writePng(path, image));
}

}  // namespace

void runSmokeRenderTests() {
    std::printf("\n--- volumetric fire and smoke: the pass ---\n");
    testTheRaysAgreeWithTheProjection();
    testAUniformSlabTransmitsExpMinusKD();
    testZeroSmokeIsTheUnsmokedImage();
    testTheInterfaceLandsWhereTheCameraSaysItDoes();
    testEveryPixelAgreesWithTheModel();
    testTheMediumStopsAtWhateverIsSolid();
    testTheSortPutsTheNearerVolumeInFront();
    testTheCameraCanBeInsideTheMedium();
    testTheFerryFireCosts();
}
