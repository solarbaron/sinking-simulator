// SPDX-License-Identifier: MIT
//
// Validation of the offscreen renderer.
//
// Every assertion is a closed form worked out before the render, in the style of
// tests/test_camera.cpp. Rendering tests decay into "looks plausible" unless the
// expected answer is written down in advance, and a plausible-looking image is
// exactly what a broken depth buffer or a mirrored Y axis produces.
//
// Skipped rather than failed when there is no GPU.
#include "engine/core/math.hpp"
#include "engine/core/png.hpp"
#include "engine/gpu/offscreen.hpp"
#include "harness.hpp"

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

using core::Image;
using gpu::MeshVertex;
using testing::expectEqual;
using testing::expectTrue;

namespace {

constexpr std::uint32_t kWidth = 256;
constexpr std::uint32_t kHeight = 256;
const std::string kScratch =
    "/tmp/claude-1000/-home-solarbaron-git/49c48569-121e-4a97-b437-25941d58fd05/scratchpad/";

bool announced = false;

// Sets up device and renderer, or reports why not. Returns false to skip.
bool setup(gpu::Device& device, gpu::OffscreenRenderer& renderer) {
    std::string error;
    if (!device.create(error)) {
        if (!announced) {
            std::printf("     no usable GPU (%s) - render checks skipped\n", error.c_str());
            announced = true;
        }
        return false;
    }
    if (!renderer.create(device, kWidth, kHeight, SHIPSIM_SHADER_DIR, error)) {
        if (!announced) {
            std::printf("     renderer unavailable (%s) - render checks skipped\n", error.c_str());
            announced = true;
        }
        return false;
    }
    return true;
}

// Round-trips through a PNG on disk, so the harness the rest of the renderer
// will rely on is exercised on every render rather than assumed.
bool renderThroughPng(gpu::OffscreenRenderer& renderer, const sim::Mat4& mvp,
                      const gpu::OffscreenRenderer::Draw& draw, const float clear[4],
                      const std::string& name, Image& out) {
    float matrix[16];
    mvp.toFloats(matrix);
    Image rendered;
    if (!renderer.render(matrix, draw, clear, rendered)) return false;
    if (!core::writePng(kScratch + name, rendered)) return false;
    return core::readPng(kScratch + name, out);
}

std::size_t countPixels(const Image& image, std::uint8_t r, std::uint8_t g, std::uint8_t b,
                        int tolerance = 2) {
    std::size_t count = 0;
    for (std::uint32_t y = 0; y < image.height; ++y)
        for (std::uint32_t x = 0; x < image.width; ++x) {
            const std::uint8_t* p = image.pixel(x, y);
            if (std::abs(p[0] - r) <= tolerance && std::abs(p[1] - g) <= tolerance &&
                std::abs(p[2] - b) <= tolerance)
                ++count;
        }
    return count;
}

// An empty draw must clear the whole viewport and nothing else -- exactly
// width * height pixels, not approximately.
void testClearCoversExactlyTheViewport() {
    gpu::Device device;
    gpu::OffscreenRenderer renderer;
    if (!setup(device, renderer)) return;

    const float clear[4] = {0.25f, 0.5f, 0.75f, 1.0f};
    Image image;
    expectTrue("an empty draw renders",
               renderThroughPng(renderer, sim::Mat4::identity(), {}, clear, "clear.png", image));
    expectEqual("the image is the requested size", static_cast<long long>(image.width), kWidth);
    expectEqual("the image is the requested height", static_cast<long long>(image.height), kHeight);

    // 0.25, 0.5, 0.75 in UNORM8 are 64, 128, 191.
    expectEqual("every pixel is the clear colour",
                static_cast<long long>(countPixels(image, 64, 128, 191)),
                static_cast<long long>(kWidth) * kHeight);
}

// A triangle covering exactly half the viewport must colour within a pixel row
// of half the pixels. Anything wildly off means the projection or the viewport
// disagrees with what the camera maths predicted.
void testTriangleCoversItsPredictedArea() {
    gpu::Device device;
    gpu::OffscreenRenderer renderer;
    if (!setup(device, renderer)) return;

    // A right triangle over the whole clip-space square covers exactly half of it.
    const std::vector<MeshVertex> vertices = {
        {{-1.0f, -1.0f, 0.5f}, {1.0f, 0.0f, 0.0f}},
        {{ 1.0f, -1.0f, 0.5f}, {1.0f, 0.0f, 0.0f}},
        {{-1.0f,  1.0f, 0.5f}, {1.0f, 0.0f, 0.0f}},
    };
    const std::vector<std::uint32_t> indices = {0, 1, 2};

    gpu::OffscreenRenderer::Draw draw{vertices.data(), vertices.size(), indices.data(),
                                      indices.size()};
    const float clear[4] = {0.0f, 0.0f, 0.0f, 1.0f};
    Image image;
    expectTrue("the triangle renders",
               renderThroughPng(renderer, sim::Mat4::identity(), draw, clear, "triangle.png",
                                image));

    const std::size_t red = countPixels(image, 255, 0, 0);
    const auto total = static_cast<double>(kWidth) * kHeight;
    const double fraction = static_cast<double>(red) / total;
    // Half, within one row of pixels' worth of rasterisation edge effects.
    const double tolerance = static_cast<double>(kWidth) / total + 0.005;
    expectTrue("the triangle covers half the viewport, within a pixel row",
               std::abs(fraction - 0.5) < tolerance);
    expectTrue("the covered corner is inside the triangle",
               countPixels(image, 255, 0, 0) > 0 &&
                   image.pixel(4, 4)[0] > 200 && image.pixel(4, 4)[1] < 50);
    expectTrue("the opposite corner is background",
               image.pixel(kWidth - 5, kHeight - 5)[0] < 50);
}

// The test that actually proves depth testing works. Drawing near-then-far gives
// the right image even with the depth buffer disabled, so the far triangle is
// drawn *second* and must still lose.
void testNearerGeometryOccludesRegardlessOfDrawOrder() {
    gpu::Device device;
    gpu::OffscreenRenderer renderer;
    if (!setup(device, renderer)) return;

    // Two full-screen quads at different depths, green nearer than blue.
    auto quad = [](float z, float r, float g, float b) {
        return std::vector<MeshVertex>{
            {{-1.0f, -1.0f, z}, {r, g, b}}, {{1.0f, -1.0f, z}, {r, g, b}},
            {{ 1.0f,  1.0f, z}, {r, g, b}}, {{-1.0f, 1.0f, z}, {r, g, b}},
        };
    };

    std::vector<MeshVertex> vertices = quad(0.3f, 0.0f, 1.0f, 0.0f);  // near, green
    const std::vector<MeshVertex> far = quad(0.7f, 0.0f, 0.0f, 1.0f);  // far, blue
    vertices.insert(vertices.end(), far.begin(), far.end());

    // Near quad first, far quad second: with depth testing off the blue would
    // overwrite the green everywhere.
    const std::vector<std::uint32_t> indices = {0, 1, 2, 0, 2, 3, 4, 5, 6, 4, 6, 7};

    gpu::OffscreenRenderer::Draw draw{vertices.data(), vertices.size(), indices.data(),
                                      indices.size()};
    const float clear[4] = {0.0f, 0.0f, 0.0f, 1.0f};
    Image image;
    expectTrue("the depth scene renders",
               renderThroughPng(renderer, sim::Mat4::identity(), draw, clear, "depth.png", image));

    const std::size_t green = countPixels(image, 0, 255, 0);
    const std::size_t blue = countPixels(image, 0, 0, 255);
    expectEqual("the nearer surface owns every pixel", static_cast<long long>(green),
                static_cast<long long>(kWidth) * kHeight);
    expectEqual("the farther surface is completely hidden", static_cast<long long>(blue), 0);

    const std::uint8_t* centre = image.pixel(kWidth / 2, kHeight / 2);
    expectTrue("the centre pixel is the near colour",
               centre[1] > 200 && centre[2] < 50);
}

// A world point must land on the pixel sim::clipToPixel predicts -- the same
// closed form the camera tests established, now checked through the GPU.
void testWorldPointLandsWhereTheCameraPredicts() {
    gpu::Device device;
    gpu::OffscreenRenderer renderer;
    if (!setup(device, renderer)) return;

    const sim::Mat4 view = sim::lookAt({0, 0, 10}, {0, 0, 0}, {0, 1, 0});
    const sim::Mat4 projection =
        sim::perspective(90.0 * sim::kDegToRad, 1.0, 0.1, 100.0);
    const sim::Mat4 mvp = projection * view;

    // A small quad centred on a deliberately off-centre world point, so a
    // mirrored axis puts it in the wrong quadrant rather than back on itself.
    // Sized from the closed form rather than picked: at 90 degrees vertical
    // field of view the visible height at ten units of depth is twenty world
    // units, so a quad of this half-extent covers a predictable pixel area. The
    // first version used half = 0.35, which works out to about 81 pixels -- and
    // the test demanded more than 100, so it failed on arithmetic rather than on
    // anything being wrong with the renderer.
    const sim::Vec3 centre{3.0, 2.0, 0.0};
    const float half = 1.0f;
    const double visibleWorldHeight = 2.0 * 10.0 * std::tan(45.0 * sim::kDegToRad);
    const double expectedSidePixels = (2.0 * half / visibleWorldHeight) * kHeight;
    const double expectedArea = expectedSidePixels * expectedSidePixels;
    const std::vector<MeshVertex> vertices = {
        {{float(centre.x) - half, float(centre.y) - half, 0.0f}, {1.0f, 1.0f, 0.0f}},
        {{float(centre.x) + half, float(centre.y) - half, 0.0f}, {1.0f, 1.0f, 0.0f}},
        {{float(centre.x) + half, float(centre.y) + half, 0.0f}, {1.0f, 1.0f, 0.0f}},
        {{float(centre.x) - half, float(centre.y) + half, 0.0f}, {1.0f, 1.0f, 0.0f}},
    };
    const std::vector<std::uint32_t> indices = {0, 1, 2, 0, 2, 3};

    gpu::OffscreenRenderer::Draw draw{vertices.data(), vertices.size(), indices.data(),
                                      indices.size()};
    const float clear[4] = {0.0f, 0.0f, 0.0f, 1.0f};
    Image image;
    expectTrue("the positioned quad renders",
               renderThroughPng(renderer, mvp, draw, clear, "position.png", image));

    // Where the camera maths says it should be.
    double clip[4], predictedX = 0, predictedY = 0;
    mvp.transform(centre, clip);
    expectTrue("the point projects in front of the eye",
               sim::clipToPixel(clip, kWidth, kHeight, predictedX, predictedY));

    // Centroid of the drawn pixels, which is where it actually is.
    double sumX = 0, sumY = 0;
    std::size_t count = 0;
    for (std::uint32_t y = 0; y < image.height; ++y)
        for (std::uint32_t x = 0; x < image.width; ++x) {
            const std::uint8_t* p = image.pixel(x, y);
            if (p[0] > 200 && p[1] > 200 && p[2] < 60) {
                sumX += x;
                sumY += y;
                ++count;
            }
        }
    // Closed form: the quad's pixel area follows from its world size and the
    // frustum, so assert against that rather than against "some".
    expectTrue("the rendered area matches the projected size",
               std::abs(static_cast<double>(count) - expectedArea) < expectedArea * 0.15);
    if (count == 0) return;

    const double actualX = sumX / static_cast<double>(count);
    const double actualY = sumY / static_cast<double>(count);
    // Within two pixels: the centroid of a rasterised quad sits half a pixel off
    // its geometric centre by construction.
    expectTrue("the rendered position matches the camera prediction in x",
               std::abs(actualX - predictedX) < 2.0);
    expectTrue("the rendered position matches the camera prediction in y",
               std::abs(actualY - predictedY) < 2.0);

    // And it must be up and to the right of centre, which pins the Y flip.
    expectTrue("a point above and right of the aim lands above and right",
               actualX > kWidth * 0.5 && actualY < kHeight * 0.5);
}

// Rendering the same scene twice must give identical bytes; the renderer is
// reused across frames and must not accumulate state.
void testRenderingIsRepeatable() {
    gpu::Device device;
    gpu::OffscreenRenderer renderer;
    if (!setup(device, renderer)) return;

    const std::vector<MeshVertex> vertices = {
        {{-0.5f, -0.5f, 0.5f}, {0.2f, 0.4f, 0.6f}},
        {{ 0.5f, -0.5f, 0.5f}, {0.2f, 0.4f, 0.6f}},
        {{ 0.0f,  0.5f, 0.5f}, {0.2f, 0.4f, 0.6f}},
    };
    const std::vector<std::uint32_t> indices = {0, 1, 2};
    gpu::OffscreenRenderer::Draw draw{vertices.data(), vertices.size(), indices.data(),
                                      indices.size()};
    const float clear[4] = {0.1f, 0.1f, 0.1f, 1.0f};

    Image first, second;
    expectTrue("first render",
               renderThroughPng(renderer, sim::Mat4::identity(), draw, clear, "repeat_a.png",
                                first));
    expectTrue("second render",
               renderThroughPng(renderer, sim::Mat4::identity(), draw, clear, "repeat_b.png",
                                second));
    expectTrue("rendering the same scene twice is byte-identical", first.rgba == second.rgba);

    // And an empty draw afterwards must clear completely, proving the previous
    // frame's colour attachment is not leaking through.
    Image cleared;
    expectTrue("a subsequent empty draw renders",
               renderThroughPng(renderer, sim::Mat4::identity(), {}, clear, "repeat_c.png",
                                cleared));
    expectEqual("the previous frame does not survive the clear",
                static_cast<long long>(countPixels(cleared, 26, 26, 26)),
                static_cast<long long>(kWidth) * kHeight);
}

}  // namespace

void runOffscreenTests() {
    std::printf("\n--- offscreen renderer ---\n");
    testClearCoversExactlyTheViewport();
    testTriangleCoversItsPredictedArea();
    testNearerGeometryOccludesRegardlessOfDrawOrder();
    testWorldPointLandsWhereTheCameraPredicts();
    testRenderingIsRepeatable();
}
