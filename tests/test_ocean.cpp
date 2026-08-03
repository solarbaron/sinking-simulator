// SPDX-License-Identifier: MIT
//
// Validation of the sea surface.
//
// Any grid of displaced quads looks like the sea, exactly as any sum of cosines
// does -- `tests/test_waves.cpp` says so about the field and it is doubly true
// once the field is a picture. So nothing here is eyeballed. Every assertion is a
// closed form fixed before the render: a flat sea's single exact colour, the
// shoelace area of a projected quadrilateral, the sign flip of a wave train after
// half a period, and the linear-interpolation error `a (1 - cos(k h / 2))`, which
// is an identity rather than a bound.
//
// The load-bearing one is `testRenderedSurfaceAgreesWithTheWaveField`. It takes a
// known world point, asks `sim::clipToPixel` where the camera puts it, reads the
// elevation the renderer actually drew at that pixel, and holds it against
// `sim::WaveField::elevation()` -- with a tolerance derived per sample from the
// quantisation, the local slope and the metres a pixel covers, not chosen. That
// assertion is the whole premise of `docs/02-simulation.md`: the wave under the
// bow has to be the wave the ship responds to. It carries a negative control,
// because a comparison that cannot fail proves nothing: the same check against
// the field one second later must reject nearly every sample.
//
// Skipped rather than failed when there is no GPU.
#include "engine/core/math.hpp"
#include "engine/core/png.hpp"
#include "engine/gpu/ocean.hpp"
#include "engine/sim/waves.hpp"
#include "harness.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

using core::Image;
using testing::expectEqual;
using testing::expectNear;
using testing::expectTrue;

namespace {

constexpr std::uint32_t kWidth = 512;
constexpr std::uint32_t kHeight = 512;

bool announced = false;

// The PNG round trip is part of every check here, so the suite must not fail
// merely because one machine's scratch directory is another machine's nothing.
// testing::scratchDir() is the shared answer; it returns a trailing separator.
std::string outputDirectory() {
    std::string dir = testing::scratchDir();
    if (!dir.empty() && dir.back() == '/') dir.pop_back();
    return dir;
}

bool setup(gpu::Device& device, gpu::OceanRenderer& renderer) {
    std::string error;
    if (!device.create(error)) {
        if (!announced) {
            std::printf("     no usable GPU (%s) - ocean render checks skipped\n", error.c_str());
            announced = true;
        }
        return false;
    }
    if (!renderer.create(device, kWidth, kHeight, SHIPSIM_SHADER_DIR, error)) {
        if (!announced) {
            std::printf("     ocean renderer unavailable (%s) - render checks skipped\n",
                        error.c_str());
            announced = true;
        }
        return false;
    }
    return true;
}

bool renderThroughPng(gpu::OceanRenderer& renderer, const sim::Mat4& mvp,
                      const gpu::OceanSurface& surface, const gpu::OceanView& view,
                      const float clear[4], const std::string& name, Image& out) {
    float matrix[16];
    mvp.toFloats(matrix);
    Image rendered;
    if (!renderer.render(matrix, surface, view, clear, rendered)) return false;
    const std::string path = outputDirectory() + "/" + name;
    if (!core::writePng(path, rendered)) return false;
    return core::readPng(path, out);
}

// A single wave train, which is the only sea whose surface has a closed form.
//
// `sim::WaveField::regular()` now exists and was written partly for this, but
// this helper deliberately keeps the SeaState form. A one-bin spectrum places its
// component at the *energy centroid*, not at the peak: Hs 3 m / Tp 9 s gives
// omega 0.837 rather than 2 pi / Tp = 0.698, so a 88 m wave rather than a 126 m
// one. Swapping to regular() would therefore either change the wave every
// resolution check here is tuned against, or require the centroid frequency as a
// magic constant. Both are worse than deriving it. tests/test_waves.cpp already
// asserts a single component is exactly a cos(k x - omega t + phi).
sim::WaveField monochromaticSea(double significantHeight, double peakPeriod, double direction) {
    sim::SeaState sea;
    sea.significantHeight = significantHeight;
    sea.peakPeriod = peakPeriod;
    sea.meanDirection = direction;
    sea.frequencyCount = 1;
    sea.directionCount = 1;
    return sim::WaveField(sea);
}

sim::WaveField spectralSea(double significantHeight, double peakPeriod, int frequencies,
                           int directions, std::uint64_t seed) {
    sim::SeaState sea;
    sea.significantHeight = significantHeight;
    sea.peakPeriod = peakPeriod;
    sea.frequencyCount = frequencies;
    sea.directionCount = directions;
    sea.seed = seed;
    return sim::WaveField(sea);
}

// Straight down, with world +x as screen up. World x then maps to the image's
// vertical axis and world y to its horizontal one *independently of elevation*,
// because a point directly under the eye moves only in the axis it varies in.
// That is what turns "long-crested along +x" into a statement about image rows.
sim::Mat4 overheadCamera(double height, double fovDegrees) {
    const sim::Mat4 view = sim::lookAt({0, 0, height}, {0, 0, 0}, {1, 0, 0});
    const sim::Mat4 projection = sim::perspective(
        fovDegrees * sim::kDegToRad, static_cast<double>(kWidth) / kHeight, 1.0, 4.0 * height);
    return projection * view;
}

// An oblique view, for the frame a human would look at.
sim::Mat4 obliqueCamera(double distance, double height, double fovDegrees) {
    const sim::Mat4 view =
        sim::lookAt({-distance, -0.45 * distance, height}, {0, 0, 0}, {0, 0, 1});
    const sim::Mat4 projection = sim::perspective(
        fovDegrees * sim::kDegToRad, static_cast<double>(kWidth) / kHeight, 1.0,
        6.0 * distance);
    return projection * view;
}

bool sameColour(const std::uint8_t* a, const std::uint8_t* b, int tolerance) {
    return std::abs(int{a[0]} - int{b[0]}) <= tolerance &&
           std::abs(int{a[1]} - int{b[1]}) <= tolerance &&
           std::abs(int{a[2]} - int{b[2]}) <= tolerance;
}

std::string label(const char* what, double value) {
    char buffer[160];
    std::snprintf(buffer, sizeof buffer, "%s %.4g", what, value);
    return buffer;
}

// --- The calm sea ------------------------------------------------------------

// Hs = 0 keeps every component and gives each of them zero amplitude, so this
// drives the whole displacement loop rather than skipping it, and the answer is
// still exactly zero. A test against an answer that is exactly statable is worth
// more than one against a plausible one.
void testCalmSeaIsExactlyFlat() {
    sim::SeaState calm;
    calm.significantHeight = 0.0;
    calm.frequencyCount = 48;
    calm.directionCount = 12;
    const sim::WaveField field(calm);
    expectEqual("a calm sea still carries all its components",
                static_cast<long long>(field.components().size()), 48 * 12);

    const gpu::OceanGrid grid{0.0, 0.0, 60.0, 32};
    gpu::OceanSurface surface;
    surface.build(field, grid, 7.5);

    expectEqual("the grid has (resolution + 1)^2 vertices",
                static_cast<long long>(surface.vertices().size()), 33 * 33);
    // Six triangle corners reference every interior vertex. Displacing per corner
    // rather than per unique vertex would multiply the cost by six for exactly
    // the same surface -- the redundancy the physics tick was caught paying, and
    // the reason this ratio is asserted rather than assumed.
    expectEqual("the grid has six indices per cell",
                static_cast<long long>(surface.indices().size()), 32 * 32 * 6);

    bool flat = true, upright = true, positioned = true;
    const double h = grid.cellSize();
    for (int j = 0; j <= grid.resolution; ++j)
        for (int i = 0; i <= grid.resolution; ++i) {
            const gpu::OceanVertex& v = surface.vertices()[surface.vertexIndex(i, j)];
            flat = flat && v.position[2] == 0.0f;
            upright =
                upright && v.normal[0] == 0.0f && v.normal[1] == 0.0f && v.normal[2] == 1.0f;
            positioned = positioned &&
                         std::abs(v.position[0] - (-grid.halfExtent + h * i)) < 1e-4 &&
                         std::abs(v.position[1] - (-grid.halfExtent + h * j)) < 1e-4;
        }
    expectTrue("every vertex of a calm sea sits exactly on z = 0", flat);
    expectTrue("every normal of a calm sea is exactly +z", upright);
    expectTrue("vertices sit on the grid the patch describes", positioned);

    // And still water, which has no components at all rather than zero ones.
    gpu::OceanSurface still;
    still.build(sim::WaveField{}, grid, 3.0);
    bool stillFlat = true;
    for (const gpu::OceanVertex& v : still.vertices()) stillFlat = stillFlat && v.position[2] == 0.0f;
    expectTrue("still water is flat too", stillFlat);
}

// A flat sea under view-independent shading has exactly one colour, and covers
// exactly the quadrilateral its four corners project to. Both are closed forms,
// and between them they drive the entire path -- vertex format, normals, push
// constants, shader, render pass, readback and PNG round trip -- against answers
// stated in advance.
void testCalmSeaRendersAFlatPlane() {
    gpu::Device device;
    gpu::OceanRenderer renderer;
    if (!setup(device, renderer)) return;

    sim::SeaState calm;
    calm.significantHeight = 0.0;
    const sim::WaveField field(calm);
    const gpu::OceanGrid grid{0.0, 0.0, 60.0, 48};
    gpu::OceanSurface surface;
    surface.build(field, grid, 11.0);

    const gpu::OceanView view;
    const float clear[4] = {0.55f, 0.72f, 0.91f, 1.0f};
    const sim::Mat4 mvp = overheadCamera(200.0, 60.0);

    Image image;
    expectTrue("the calm sea renders",
               renderThroughPng(renderer, mvp, surface, view, clear, "ocean_calm.png", image));
    if (!image.valid() || image.width != kWidth) return;

    // n = (0, 0, 1) everywhere, so sky = 1 and ndotl = sun.z, giving
    //     colour = waterColour * (ambient + sunStrength * sun.z)
    // stored as round(255 c). Written out here rather than asked of the renderer:
    // a test that gets the expected answer from the code under test cannot catch
    // a wrong answer.
    const double factor = static_cast<double>(view.ambient) +
                          static_cast<double>(view.sunStrength) * view.sunDirection[2];
    int expected[3];
    for (int c = 0; c < 3; ++c) {
        const double linear =
            std::clamp(static_cast<double>(view.waterColour[c]) * factor, 0.0, 1.0);
        expected[c] = static_cast<int>(std::lround(linear * 255.0));
    }

    // Where the patch lands, from the camera maths alone. A projective map takes
    // a planar quadrilateral to a quadrilateral, so the shoelace area of the four
    // projected corners is the covered area exactly, not approximately.
    const double e = grid.halfExtent;
    const sim::Vec3 worldCorner[4] = {{-e, -e, 0}, {e, -e, 0}, {e, e, 0}, {-e, e, 0}};
    double corner[4][2];
    bool inside = true;
    for (int i = 0; i < 4; ++i) {
        double clip[4];
        mvp.transform(worldCorner[i], clip);
        inside = inside && sim::clipToPixel(clip, kWidth, kHeight, corner[i][0], corner[i][1]);
        inside = inside && corner[i][0] > 1.0 && corner[i][0] < kWidth - 1.0 &&
                 corner[i][1] > 1.0 && corner[i][1] < kHeight - 1.0;
    }
    expectTrue("the whole patch projects inside the viewport", inside);

    double twiceArea = 0.0, perimeter = 0.0;
    for (int i = 0; i < 4; ++i) {
        const int j = (i + 1) & 3;
        twiceArea += corner[i][0] * corner[j][1] - corner[j][0] * corner[i][1];
        perimeter += std::hypot(corner[j][0] - corner[i][0], corner[j][1] - corner[i][1]);
    }
    const double projectedArea = std::abs(twiceArea) * 0.5;

    // Pixel (0, 0) is outside the patch by the check above, so it is the
    // background, measured rather than assumed.
    const std::uint8_t* background = image.pixel(0, 0);
    std::size_t sea = 0, uniform = 0;
    int seen[3] = {-1, -1, -1};
    for (std::uint32_t y = 0; y < image.height; ++y)
        for (std::uint32_t x = 0; x < image.width; ++x) {
            const std::uint8_t* p = image.pixel(x, y);
            if (sameColour(p, background, 2)) continue;
            if (seen[0] < 0) { seen[0] = p[0]; seen[1] = p[1]; seen[2] = p[2]; }
            ++sea;
            if (p[0] == seen[0] && p[1] == seen[1] && p[2] == seen[2]) ++uniform;
        }

    // Guard against a vacuous pass: an empty frame is also perfectly uniform.
    expectTrue("the calm sea covered a large part of the frame",
               static_cast<double>(sea) > 0.15 * kWidth * kHeight);
    expectEqual("a flat sea is exactly one colour, to the last bit",
                static_cast<long long>(uniform), static_cast<long long>(sea));
    for (int c = 0; c < 3; ++c)
        expectTrue(label("the flat sea's colour matches the closed form in channel", c),
                   std::abs(seen[c] - expected[c]) <= 1);
    // Only the boundary row of pixels is in doubt, so the perimeter is the bound.
    expectTrue(label("the flat sea covers exactly its projected quadrilateral, out by",
                     static_cast<double>(sea) - projectedArea),
               std::abs(static_cast<double>(sea) - projectedArea) < perimeter);
}

// --- The grid is the wave field ----------------------------------------------

// The constraint the whole task turns on, checked at the geometry level before
// any pixel is involved: every vertex of the render grid carries the elevation
// `WaveField::elevation()` reports at that vertex's own world position.
void testGridElevationMatchesTheWaveField() {
    const double t = 13.7;
    const sim::WaveField field = spectralSea(3.0, 9.0, 48, 12, 0x5eaf00d);
    // Off-centre, so a patch that quietly assumed the origin shows up.
    const gpu::OceanGrid grid{120.0, -80.0, 100.0, 96};
    gpu::OceanSurface surface;
    surface.build(field, grid, t);

    double worst = 0.0, lowest = 1e30, highest = -1e30, worstPosition = 0.0;
    double alongX = 0.0, alongY = 0.0;
    const double h = grid.cellSize();
    for (int j = 0; j <= grid.resolution; ++j)
        for (int i = 0; i <= grid.resolution; ++i) {
            const gpu::OceanVertex& v = surface.vertices()[surface.vertexIndex(i, j)];
            // At the grid point the vertex stands for, not at the float32 position
            // it stores. The first version of this compared against
            // elevation(double(v.position[0]), ...) and failed at 1.6e-6 m -- the
            // builder evaluates the field at the *unrounded* grid point, so that
            // comparison was measuring the position quantisation times the local
            // slope, which is a different and much larger quantity than the one
            // the check is about. Position is asserted separately, below.
            const double x = grid.centreX - grid.halfExtent + h * i;
            const double y = grid.centreY - grid.halfExtent + h * j;
            const double want = field.elevation(x, y, t);
            worst = std::max(worst, std::abs(static_cast<double>(v.position[2]) - want));
            worstPosition = std::max({worstPosition, std::abs(v.position[0] - x),
                                      std::abs(v.position[1] - y)});
            lowest = std::min(lowest, want);
            highest = std::max(highest, want);
        }
    for (int i = 0; i <= grid.resolution; ++i) {
        const double a = surface.vertices()[surface.vertexIndex(i, 0)].position[2];
        const double b = surface.vertices()[surface.vertexIndex(0, i)].position[2];
        alongX = std::max(alongX, std::abs(a - surface.vertices()[surface.vertexIndex(0, 0)]
                                                   .position[2]));
        alongY = std::max(alongY, std::abs(b - surface.vertices()[surface.vertexIndex(0, 0)]
                                                   .position[2]));
    }

    // float32 storage is the whole of the error. The row recurrence agrees with
    // direct evaluation to 2e-14 m, and an elevation of a few metres held in
    // float32 rounds at about 5e-7 m.
    expectTrue(label("every grid vertex carries the wave field's own elevation, worst", worst),
               worst < 1e-6);
    // A patch reaching 220 m from the origin rounds to 220 * 2^-24 = 1.3e-5 m.
    expectTrue(label("every vertex stands on the grid point it represents, worst", worstPosition),
               worstPosition < 2e-5);
    // Guards: agreement is worth nothing if there is nothing to agree about, and
    // a surface that varied only along x would satisfy the first check just as
    // well on a short-crested sea that should vary in both.
    expectTrue("the surface is a sea, not a plane", highest - lowest > 2.0);
    expectTrue("a short-crested sea varies along x", alongX > 0.3);
    expectTrue("a short-crested sea varies along y", alongY > 0.3);
}

// Normals come from the spectrum's analytic slope. Checked against a central
// difference of `elevation()`, which is a different computation from the sin()
// the builder accumulates, with the tolerance derived from the difference's own
// truncation error rather than picked.
void testNormalsAreTheAnalyticSurfaceSlope() {
    const double t = 4.25;
    const sim::WaveField field = spectralSea(2.5, 8.0, 16, 5, 0x1234);
    const gpu::OceanGrid grid{0.0, 0.0, 80.0, 64};
    gpu::OceanSurface surface;
    surface.build(field, grid, t);

    // A central difference of step g is out by (g^2 / 6) eta''', and eta''' is at
    // most the sum of a k^3 over the components.
    const double step = 0.01;
    double thirdDerivative = 0.0;
    for (const sim::WaveComponent& c : field.components())
        thirdDerivative += c.amplitude * c.wavenumber * c.wavenumber * c.wavenumber;
    const double tolerance = step * step / 6.0 * thirdDerivative + 1e-6;

    double worst = 0.0, steepest = 0.0;
    const double h = grid.cellSize();
    for (int j = 0; j <= grid.resolution; j += 5)
        for (int i = 0; i <= grid.resolution; i += 5) {
            const gpu::OceanVertex& v = surface.vertices()[surface.vertexIndex(i, j)];
            // The unrounded grid point, for the same reason as in the elevation
            // check above: the builder evaluated the slope there.
            const double x = -grid.halfExtent + h * i;
            const double y = -grid.halfExtent + h * j;
            const double gx = (field.elevation(x + step, y, t) - field.elevation(x - step, y, t)) /
                              (2.0 * step);
            const double gy = (field.elevation(x, y + step, t) - field.elevation(x, y - step, t)) /
                              (2.0 * step);
            const double inverse = 1.0 / std::sqrt(gx * gx + gy * gy + 1.0);
            const double want[3] = {-gx * inverse, -gy * inverse, inverse};
            for (int c = 0; c < 3; ++c)
                worst = std::max(worst, std::abs(static_cast<double>(v.normal[c]) - want[c]));
            steepest = std::max(steepest, std::hypot(gx, gy));
        }
    expectTrue(label("normals are the analytic surface slope, worst", worst), worst < tolerance);
    // Guard: normals of (0, 0, 1) everywhere would also agree with a flat field.
    expectTrue(label("the normals actually tilt, steepest slope", steepest), steepest > 0.05);
}

// --- The load-bearing check ---------------------------------------------------

// A world point on the surface, projected with sim::clipToPixel, must be drawn
// at the elevation sim::WaveField::elevation() reports for it.
void testRenderedSurfaceAgreesWithTheWaveField() {
    gpu::Device device;
    gpu::OceanRenderer renderer;
    if (!setup(device, renderer)) return;

    const double t = 21.5;
    const sim::WaveField field = spectralSea(3.0, 9.0, 32, 8, 0xC0FFEE);
    // 140 m across at 0.625 m cells. Sized so that neither of the two things the
    // comparison cannot see past -- the grid's own interpolation error and the
    // sub-pixel offset of the pixel centre -- is larger than the millimetres the
    // elevation channel resolves to.
    const gpu::OceanGrid grid{0.0, 0.0, 70.0, 224};
    gpu::OceanSurface surface;
    surface.build(field, grid, t);

    gpu::OceanView view;
    view.shading = gpu::OceanShading::Elevation;
    // Nothing in a 3 m sea reaches 8 m, so nothing saturates; 16 m over 16 bits
    // quantises at 0.24 mm.
    view.elevationMin = -8.0f;
    view.elevationSpan = 16.0f;

    // Straight down from 280 m at 30 degrees: the patch fills the frame, and a
    // near-vertical view means no crest can hide the water behind it, so the ray
    // through a surface point's own pixel really does land on that point.
    const sim::Mat4 mvp = overheadCamera(280.0, 30.0);
    const float clear[4] = {0.0f, 0.0f, 0.0f, 1.0f};  // blue 0, so the surface tag is unambiguous

    Image image;
    expectTrue("the elevation channel renders",
               renderThroughPng(renderer, mvp, surface, view, clear, "ocean_elevation.png", image));
    if (!image.valid() || image.width != kWidth) return;

    const double quantisation = 0.5 * static_cast<double>(view.elevationSpan) / 65535.0;
    const double h = grid.cellSize();
    const double step = 0.05;

    // The piecewise-linear surface a *correct* grid of this spacing would have,
    // evaluated straight from the wave field. The ray meets a mesh, not the
    // analytic surface, so the gap between the two belongs in the tolerance -- but
    // it must be computed from the field, never from the artefact under test.
    //
    // The first version measured it with `surface.sampleElevation()`, and that is
    // a defect rather than a shortcut: a surface with a systematic elevation error
    // inflates its own tolerance by exactly the amount it is wrong by. A mutation
    // test settled it -- scaling every elevation by 1.02 was caught by four other
    // assertions here and sailed straight through this one.
    const auto meshFromField = [&](double x, double y) {
        const double fx = std::clamp((x + grid.halfExtent) / h, 0.0, grid.resolution - 1e-9);
        const double fy = std::clamp((y + grid.halfExtent) / h, 0.0, grid.resolution - 1e-9);
        const int i = static_cast<int>(std::floor(fx));
        const int j = static_cast<int>(std::floor(fy));
        const double u = fx - i, v = fy - j;
        const double cornerX = -grid.halfExtent + i * h;
        const double cornerY = -grid.halfExtent + j * h;
        const double z00 = field.elevation(cornerX, cornerY, t);
        const double z11 = field.elevation(cornerX + h, cornerY + h, t);
        // Same v00-v11 diagonal, and therefore the same two triangles, that
        // OceanSurface builds and sampleElevation interpolates over.
        if (v <= u) {
            const double z10 = field.elevation(cornerX + h, cornerY, t);
            return z00 + u * (z10 - z00) + v * (z11 - z10);
        }
        const double z01 = field.elevation(cornerX, cornerY + h, t);
        return z00 + v * (z01 - z00) + u * (z11 - z01);
    };

    long long samples = 0, exceeded = 0, staleRejected = 0;
    double worstError = 0.0, sumError = 0.0, sumSquares = 0.0;
    double worstBound = 0.0, worstRatio = 0.0, lowest = 1e30, highest = -1e30;

    for (int j = 6; j < grid.resolution; j += 11)
        for (int i = 6; i < grid.resolution; i += 11) {
            const gpu::OceanVertex& v = surface.vertices()[surface.vertexIndex(i, j)];
            const sim::Vec3 point{v.position[0], v.position[1], v.position[2]};
            double clip[4], px = 0, py = 0;
            mvp.transform(point, clip);
            if (!sim::clipToPixel(clip, kWidth, kHeight, px, py)) continue;
            if (px < 1.0 || py < 1.0 || px > kWidth - 2.0 || py > kHeight - 2.0) continue;

            const auto column = static_cast<std::uint32_t>(px);
            const auto row = static_cast<std::uint32_t>(py);
            double rendered = 0.0;
            if (!gpu::decodeOceanElevation(view, image.pixel(column, row), rendered)) continue;

            // The Jacobian of the projection here, from two short probes through
            // the same matrix. Columns are d(pixel)/dx and d(pixel)/dy.
            double clipX[4], clipY[4], ax = 0, ay = 0, bx = 0, by = 0;
            mvp.transform({point.x + step, point.y, point.z}, clipX);
            mvp.transform({point.x, point.y + step, point.z}, clipY);
            if (!sim::clipToPixel(clipX, kWidth, kHeight, ax, ay)) continue;
            if (!sim::clipToPixel(clipY, kWidth, kHeight, bx, by)) continue;
            const double j00 = (ax - px) / step, j01 = (bx - px) / step;
            const double j10 = (ay - py) / step, j11 = (by - py) / step;
            const double determinant = j00 * j11 - j01 * j10;
            if (std::abs(determinant) < 1e-9) continue;

            // The pixel *centre*, not the point's own projection, is what the
            // rasteriser shaded. Rather than widen the tolerance to swallow that
            // offset, solve for the world point on the analytic surface whose
            // projection *is* the pixel centre, and ask the field about that. What
            // is left over is the thing the test is actually about.
            //
            // Newton rather than one linear step: the surface height changes as
            // the solution moves, and that shifts the projection again. Ignoring
            // it left four samples of four hundred outside their bound -- a real
            // 2 mm effect, not noise, and the sort of residual that would have
            // been quietly absorbed by rounding the tolerance up.
            double dx = 0.0, dy = 0.0;
            for (int iteration = 0; iteration < 4; ++iteration) {
                const double z = field.elevation(point.x + dx, point.y + dy, t);
                double clipHere[4], qx = 0, qy = 0;
                mvp.transform({point.x + dx, point.y + dy, z}, clipHere);
                if (!sim::clipToPixel(clipHere, kWidth, kHeight, qx, qy)) break;
                const double ex = (column + 0.5) - qx, ey = (row + 0.5) - qy;
                dx += (j11 * ex - j01 * ey) / determinant;
                dy += (-j10 * ex + j00 * ey) / determinant;
            }
            const double hitX = point.x + dx, hitY = point.y + dy;
            const double analytic = field.elevation(hitX, hitY, t);
            const double error = rendered - analytic;

            // What is left: the 0.24 mm quantisation and the exact interpolation
            // error a correct grid of this spacing has at this very point.
            //
            // Times three. The height error shifts where the ray meets the mesh and
            // the slope turns that shift back into height, which accounts for about
            // 1.1; the rest is that the interpolation error vanishes at the
            // vertices, so evaluating it a few millimetres from where the ray
            // actually lands changes it by a large *relative* amount while staying
            // in millimetres absolutely. Measured worst ratio over these samples is
            // 1.57 and the run prints it, so the margin is visible rather than
            // assumed -- and a 2% elevation error, which puts 47 mm on the board,
            // is still rejected by an order of magnitude.
            const double interpolation = std::abs(meshFromField(hitX, hitY) - analytic);
            const double bound = quantisation + 3.0 * interpolation;

            ++samples;
            worstError = std::max(worstError, std::abs(error));
            worstBound = std::max(worstBound, bound);
            if (interpolation > 1e-9)
                worstRatio = std::max(worstRatio, (std::abs(error) - quantisation) / interpolation);
            if (std::abs(error) > bound) ++exceeded;
            sumError += error;
            sumSquares += error * error;
            lowest = std::min(lowest, analytic);
            highest = std::max(highest, analytic);

            // Negative control. One second is a small, entirely plausible error --
            // a renderer a frame or two behind the physics -- and the same
            // comparison must reject it. Without this the tolerance above could be
            // loose enough to accept anything.
            if (std::abs(rendered - field.elevation(point.x + dx, point.y + dy, t + 1.0)) > bound)
                ++staleRejected;
        }

    expectTrue("the check had samples to work with", samples > 200);
    if (samples == 0) return;
    const double mean = sumError / static_cast<double>(samples);
    const double rms = std::sqrt(sumSquares / static_cast<double>(samples));

    std::printf("     rendered vs WaveField::elevation: %lld samples, mean %+.5f m, rms %.5f m,"
                " worst %.5f m against a worst bound of %.5f m (ratio needed %.2f)\n",
                samples, mean, rms, worstError, worstBound, worstRatio);

    expectEqual("every rendered point agrees with WaveField::elevation inside its own bound",
                exceeded, 0);
    // A residual error is as likely one way as the other; a sign flip, a stale
    // time, a scale factor or a shifted patch would not be.
    expectNear("the rendered surface carries no systematic offset", mean, 0.0, 0.002);
    // Guards against a vacuous pass.
    expectTrue(label("the sampled surface varies by", highest - lowest), highest - lowest > 2.0);
    expectTrue("the agreement is three orders tighter than the variation being measured",
               rms * 1000.0 < highest - lowest);
    expectTrue(label("a one-second-stale field is rejected, fraction",
                     static_cast<double>(staleRejected) / static_cast<double>(samples)),
               staleRejected > samples * 8 / 10);
}

// --- Determinism and time -----------------------------------------------------

void testRenderingIsRepeatable() {
    gpu::Device device;
    gpu::OceanRenderer renderer;
    if (!setup(device, renderer)) return;

    const sim::WaveField field = spectralSea(2.0, 7.0, 24, 6, 99);
    const gpu::OceanGrid grid{0.0, 0.0, 90.0, 64};
    gpu::OceanSurface first, second;
    first.build(field, grid, 6.25);
    second.build(field, grid, 6.25);

    bool identicalGeometry = first.vertices().size() == second.vertices().size();
    for (std::size_t i = 0; identicalGeometry && i < first.vertices().size(); ++i)
        identicalGeometry =
            std::memcmp(&first.vertices()[i], &second.vertices()[i], sizeof(gpu::OceanVertex)) == 0;
    expectTrue("building the same sea twice gives identical geometry", identicalGeometry);

    const gpu::OceanView view;
    const float clear[4] = {0.05f, 0.09f, 0.14f, 1.0f};
    const sim::Mat4 mvp = obliqueCamera(160.0, 40.0, 45.0);

    Image a, b;
    expectTrue("first render",
               renderThroughPng(renderer, mvp, first, view, clear, "ocean_repeat_a.png", a));
    expectTrue("second render",
               renderThroughPng(renderer, mvp, second, view, clear, "ocean_repeat_b.png", b));
    expectTrue("rendering the same sea at the same time twice is byte-identical", a.rgba == b.rgba);

    // Guard: two blank frames are also byte-identical.
    std::size_t drawn = 0;
    for (std::uint32_t y = 0; y < a.height; ++y)
        for (std::uint32_t x = 0; x < a.width; ++x)
            if (!sameColour(a.pixel(x, y), a.pixel(0, 0), 2)) ++drawn;
    expectTrue("the repeated frames had a sea in them",
               static_cast<double>(drawn) > 0.10 * kWidth * kHeight);

    // An unbuilt surface afterwards must clear completely, so no geometry from
    // the previous frame survives in the buffers.
    gpu::OceanSurface empty;
    Image cleared;
    expectTrue("an empty surface renders",
               renderThroughPng(renderer, mvp, empty, view, clear, "ocean_empty.png", cleared));
    std::size_t background = 0;
    for (std::uint32_t y = 0; y < cleared.height; ++y)
        for (std::uint32_t x = 0; x < cleared.width; ++x)
            if (sameColour(cleared.pixel(x, y), cleared.pixel(0, 0), 0)) ++background;
    expectEqual("the previous frame does not survive the clear",
                static_cast<long long>(background), static_cast<long long>(kWidth) * kHeight);
}

// Time is not frozen: a wave train half a period on is exactly its own negative,
// and a spectral sea two seconds on redraws almost every pixel.
void testAdvancingTimeChangesTheSea() {
    const sim::WaveField train = monochromaticSea(3.0, 9.0, 0.0);
    const sim::WaveComponent c = train.components()[0];
    const double period = 2.0 * sim::kPi / c.omega;
    const gpu::OceanGrid grid{0.0, 0.0, 120.0, 64};

    gpu::OceanSurface now, later;
    now.build(train, grid, 3.0);
    later.build(train, grid, 3.0 + 0.5 * period);

    double worst = 0.0, largest = 0.0;
    for (std::size_t i = 0; i < now.vertices().size(); ++i) {
        const double a = now.vertices()[i].position[2];
        const double b = later.vertices()[i].position[2];
        worst = std::max(worst, std::abs(a + b));
        largest = std::max(largest, std::abs(a));
    }
    // cos(psi - pi) = -cos(psi), exactly.
    expectTrue(label("half a period on, a wave train is its own negative, worst", worst),
               worst < 1e-6);
    // Guard: 0 == -0 would satisfy that on a flat sea.
    expectTrue(label("the wave train has an amplitude, largest", largest), largest > 0.5);

    gpu::Device device;
    gpu::OceanRenderer renderer;
    if (!setup(device, renderer)) return;

    const sim::WaveField field = spectralSea(3.0, 9.0, 24, 6, 7);
    gpu::OceanSurface early, late;
    early.build(field, grid, 0.0);
    late.build(field, grid, 2.0);

    const gpu::OceanView view;
    const float clear[4] = {0.05f, 0.09f, 0.14f, 1.0f};
    const sim::Mat4 mvp = obliqueCamera(180.0, 55.0, 45.0);

    Image before, after;
    expectTrue("the sea renders at t = 0",
               renderThroughPng(renderer, mvp, early, view, clear, "ocean_t0.png", before));
    expectTrue("the sea renders at t = 2",
               renderThroughPng(renderer, mvp, late, view, clear, "ocean_t2.png", after));

    std::size_t sea = 0, changed = 0;
    for (std::uint32_t y = 0; y < before.height; ++y)
        for (std::uint32_t x = 0; x < before.width; ++x) {
            if (sameColour(before.pixel(x, y), before.pixel(0, 0), 2)) continue;
            ++sea;
            if (!sameColour(before.pixel(x, y), after.pixel(x, y), 1)) ++changed;
        }
    expectTrue("there was a sea to compare", static_cast<double>(sea) > 0.10 * kWidth * kHeight);
    expectTrue(label("advancing time redraws the sea, fraction changed",
                     static_cast<double>(changed) / static_cast<double>(std::max<std::size_t>(sea, 1))),
               static_cast<double>(changed) > 0.6 * static_cast<double>(sea));
}

// --- Long-crested -------------------------------------------------------------

// A wave train travelling along +x varies along x and not along y. With one
// direction the field's dirY is exactly zero, so the geometry claim is exact
// rather than a tolerance; the image claim is a tolerance, because it is about
// the rasteriser.
void testLongCrestedSeaVariesAlongXOnly() {
    const sim::WaveField train = monochromaticSea(4.0, 6.0, 0.0);
    const sim::WaveComponent c = train.components()[0];
    const double amplitude = c.amplitude;
    const double wavelength = 2.0 * sim::kPi / c.wavenumber;
    expectNear("the wave train travels along +x", c.dirY, 0.0, 0.0);

    // 400 m across, more than a wavelength in every direction the camera sees.
    const gpu::OceanGrid grid{0.0, 0.0, 200.0, 256};
    gpu::OceanSurface surface;
    surface.build(train, grid, 2.5);

    bool crestwiseIdentical = true, normalFlatInY = true;
    double alongCrest = 0.0;
    for (int i = 0; i <= grid.resolution; ++i) {
        const gpu::OceanVertex& reference = surface.vertices()[surface.vertexIndex(i, 0)];
        for (int j = 1; j <= grid.resolution; ++j) {
            const gpu::OceanVertex& v = surface.vertices()[surface.vertexIndex(i, j)];
            crestwiseIdentical = crestwiseIdentical && v.position[2] == reference.position[2];
            normalFlatInY = normalFlatInY && v.normal[1] == 0.0f;
        }
        alongCrest = std::max(alongCrest,
                              std::abs(static_cast<double>(reference.position[2]) -
                                       static_cast<double>(surface.vertices()[surface.vertexIndex(0, 0)]
                                                               .position[2])));
    }
    expectTrue("a +x wave train is bit-identical along every crest", crestwiseIdentical);
    expectTrue("its normals have no y component at all", normalFlatInY);
    expectTrue(label("and it does vary along x, by", alongCrest), alongCrest > amplitude);

    gpu::Device device;
    gpu::OceanRenderer renderer;
    if (!setup(device, renderer)) return;

    // Overhead with +x as screen up, and a patch wider than the frustum, so every
    // pixel is sea and every image *row* is one world x.
    gpu::OceanView view;
    view.shading = gpu::OceanShading::Elevation;
    const float clear[4] = {0.0f, 0.0f, 0.0f, 1.0f};
    const sim::Mat4 mvp = overheadCamera(200.0, 60.0);

    Image image;
    expectTrue("the long-crested sea renders",
               renderThroughPng(renderer, mvp, surface, view, clear, "ocean_longcrested.png",
                                image));
    if (!image.valid() || image.width != kWidth) return;

    std::size_t constantRows = 0, coveredRows = 0;
    int worstRowSpread = 0;
    double lowest = 1e30, highest = -1e30;
    for (std::uint32_t y = 0; y < image.height; ++y) {
        int lo[3] = {255, 255, 255}, hi[3] = {0, 0, 0};
        bool covered = true;
        for (std::uint32_t x = 0; x < image.width; ++x) {
            const std::uint8_t* p = image.pixel(x, y);
            if (p[2] < 128) { covered = false; break; }
            for (int ch = 0; ch < 3; ++ch) {
                lo[ch] = std::min(lo[ch], int{p[ch]});
                hi[ch] = std::max(hi[ch], int{p[ch]});
            }
        }
        if (!covered) continue;
        ++coveredRows;
        const int spread = std::max(hi[0] - lo[0], hi[1] - lo[1]);
        worstRowSpread = std::max(worstRowSpread, spread);
        if (spread == 0) ++constantRows;
        double elevation = 0.0;
        if (gpu::decodeOceanElevation(view, image.pixel(image.width / 2, y), elevation)) {
            lowest = std::min(lowest, elevation);
            highest = std::max(highest, elevation);
        }
    }

    std::printf("     long-crested: %zu/%u rows covered, %zu bit-constant, worst row spread %d\n",
                coveredRows, kHeight, constantRows, worstRowSpread);

    expectEqual("the patch covers every row of the frame", static_cast<long long>(coveredRows),
                static_cast<long long>(kHeight));
    // World y maps to the image's horizontal axis under this camera, so a row is
    // a line of constant x. Two codes of slack: the geometry is constant along the
    // row exactly, and what is left is the rasteriser's own interpolation.
    expectTrue(label("every image row is constant across the crests, worst spread",
                     worstRowSpread),
               worstRowSpread <= 2);
    // Closed form for the other axis: the frame spans several wavelengths, so the
    // decoded elevation must run from -a to +a. Tolerance is the 0.24 mm
    // quantisation plus the mesh's a (1 - cos(k h / 2)) at the crest, which for
    // this grid is 1.4 cm.
    const double meshBound =
        amplitude * (1.0 - std::cos(0.5 * c.wavenumber * grid.cellSize())) + 1e-3;
    expectTrue(label("the frame spans more than a wavelength, wavelength", wavelength),
               wavelength < 200.0);
    expectNear("the crest reaches the wave train's amplitude", highest, amplitude, meshBound);
    expectNear("and the trough reaches its negative", lowest, -amplitude, meshBound);
}

// The lighting answers to the waves, not to a constant.
//
// The first version of this counted pixels more than three codes brighter than
// flat water and demanded 5% of them, and failed -- not because the shading was
// wrong but because the assertion had no closed form behind it and no idea how
// much brightness a 3.5 m sea is actually entitled to span. The right question is
// answerable exactly: run every vertex normal through the shading expression the
// fragment shader documents and the range the image may show follows, so the test
// asserts that range instead of a guessed fraction.
void testLightingRespondsToTheWaves() {
    gpu::Device device;
    gpu::OceanRenderer renderer;
    if (!setup(device, renderer)) return;

    const double t = 17.0;
    const sim::WaveField field = spectralSea(3.5, 9.0, 32, 8, 0xBEEF);
    // Inside the 107 m the 400 m / 30 degree frustum sees at z = 0, so every
    // vertex of the patch contributes to the frame and the predicted range is
    // over the same surface the image shows.
    const gpu::OceanGrid grid{0.0, 0.0, 100.0, 160};
    gpu::OceanSurface wavy;
    wavy.build(field, grid, t);

    const gpu::OceanView view;
    const float clear[4] = {0.55f, 0.72f, 0.91f, 1.0f};

    // colour = waterColour * (ambient * (0.5 + 0.5 n.z) + sunStrength * max(n . sun, 0)),
    // stored as round(255 c). Written out here, not asked of the renderer.
    const auto blueOf = [&](const float n[3]) {
        const double ndotl = std::max(static_cast<double>(n[0]) * view.sunDirection[0] +
                                          static_cast<double>(n[1]) * view.sunDirection[1] +
                                          static_cast<double>(n[2]) * view.sunDirection[2],
                                      0.0);
        const double sky = 0.5 + 0.5 * static_cast<double>(n[2]);
        const double lit = std::clamp(static_cast<double>(view.waterColour[2]) *
                                          (static_cast<double>(view.ambient) * sky +
                                           static_cast<double>(view.sunStrength) * ndotl),
                                      0.0, 1.0);
        return static_cast<int>(std::lround(lit * 255.0));
    };
    const float flatNormal[3] = {0.0f, 0.0f, 1.0f};
    const int flatBlue = blueOf(flatNormal);

    int predictedLow = 255, predictedHigh = 0;
    for (const gpu::OceanVertex& v : wavy.vertices()) {
        const int blue = blueOf(v.normal);
        predictedLow = std::min(predictedLow, blue);
        predictedHigh = std::max(predictedHigh, blue);
    }

    const sim::Mat4 mvp = overheadCamera(400.0, 30.0);
    Image image;
    expectTrue("the lit sea renders",
               renderThroughPng(renderer, mvp, wavy, view, clear, "ocean_lit.png", image));
    if (!image.valid() || image.width != kWidth) return;

    int renderedLow = 255, renderedHigh = 0;
    std::size_t sea = 0;
    for (std::uint32_t y = 0; y < image.height; ++y)
        for (std::uint32_t x = 0; x < image.width; ++x) {
            const std::uint8_t* p = image.pixel(x, y);
            if (sameColour(p, image.pixel(0, 0), 2)) continue;  // background
            ++sea;
            renderedLow = std::min(renderedLow, int{p[2]});
            renderedHigh = std::max(renderedHigh, int{p[2]});
        }

    std::printf("     shading: flat %d, vertex normals predict [%d, %d], rendered [%d, %d]\n",
                flatBlue, predictedLow, predictedHigh, renderedLow, renderedHigh);

    expectTrue("there is a sea to shade", static_cast<double>(sea) > 0.4 * kWidth * kHeight);
    // Guard: matching a predicted range proves nothing if the range is one value.
    expectTrue(label("the wave normals predict a real spread of brightness, codes",
                     predictedHigh - predictedLow),
               predictedHigh - predictedLow >= 8);
    expectTrue("waves catch more light than flat water somewhere", renderedHigh > flatBlue);
    expectTrue("and less somewhere else", renderedLow < flatBlue);
    // The shading is linear in the normal, and the interpolated normal is a convex
    // combination renormalised by a factor within 1e-4 of one for neighbours this
    // close, so a fragment cannot be meaningfully outside the vertex range. Two
    // codes of slack covers that and the UNORM8 rounding.
    expectTrue("no fragment is lit outside what the wave normals allow",
               renderedHigh <= predictedHigh + 2 && renderedLow >= predictedLow - 2);
    // And it fills that range, which is what says the shading is driven by these
    // normals rather than by something that merely varies.
    expectTrue("the rendered brightness fills the range the normals predict",
               renderedHigh - renderedLow >= (predictedHigh - predictedLow) * 3 / 4);

    // The frame a person would actually look at, kept because a renderer whose
    // only validated view is straight down is a renderer nobody has looked at.
    // A plane seen obliquely is foreshortened into a band, so the coverage bar is
    // the tenth of the frame the other oblique checks use, not the half an
    // overhead view would give.
    Image oblique;
    expectTrue("the sea renders from an oblique camera",
               renderThroughPng(renderer, obliqueCamera(150.0, 35.0, 55.0), wavy, view, clear,
                                "ocean_shaded.png", oblique));
    std::size_t obliqueSea = 0;
    for (std::uint32_t y = 0; y < oblique.height; ++y)
        for (std::uint32_t x = 0; x < oblique.width; ++x)
            if (!sameColour(oblique.pixel(x, y), oblique.pixel(0, 0), 2)) ++obliqueSea;
    std::printf("     oblique view covers %.0f%% of the frame\n",
                100.0 * static_cast<double>(obliqueSea) / (double{kWidth} * kHeight));
    expectTrue("the oblique view has a sea in it",
               static_cast<double>(obliqueSea) > 0.10 * kWidth * kHeight);
}

// --- Resolution ---------------------------------------------------------------

// A cell centre lies on the shared diagonal, where the mesh height is the mean of
// the two diagonal corners. For eta = a cos(k x + phi) that mean is
// a cos(psi_mid) cos(k h / 2), so the error there is exactly
//
//     a (1 - cos(k h / 2)) cos(psi_mid)
//
// -- an identity, not a bound. It is asserted for a grid that resolves the wave
// and for one that does not, so the failure mode is demonstrated rather than
// assumed.
void testGridMustResolveTheWavelength() {
    const sim::WaveField train = monochromaticSea(3.0, 9.0, 0.0);
    const sim::WaveComponent c = train.components()[0];
    const double wavelength = 2.0 * sim::kPi / c.wavenumber;
    const double amplitude = c.amplitude;
    expectTrue("the wave train has an amplitude worth measuring", amplitude > 0.5);
    expectNear("dominantWavelength finds the only wave there is",
               gpu::dominantWavelength(train), wavelength, 1e-9 * wavelength);
    expectNear("shortestWavelength finds it too", gpu::shortestWavelength(train), wavelength,
               1e-9 * wavelength);

    double relative[2] = {0.0, 0.0};
    const double cellsPerWavelength[2] = {16.0, 4.0};
    for (int variant = 0; variant < 2; ++variant) {
        const int cells = 64;
        const double h = wavelength / cellsPerWavelength[variant];
        const gpu::OceanGrid grid{0.0, 0.0, 0.5 * cells * h, cells};
        expectEqual(label("oceanResolutionFor asks for this grid at cells per wavelength",
                          cellsPerWavelength[variant]),
                    gpu::oceanResolutionFor(grid.halfExtent, wavelength,
                                            cellsPerWavelength[variant]),
                    cells);

        gpu::OceanSurface surface;
        surface.build(train, grid, 0.0);

        double worst = 0.0, largestCosine = 0.0;
        for (int j = 0; j < cells; ++j)
            for (int i = 0; i < cells; ++i) {
                const double x = -grid.halfExtent + (i + 0.5) * h;
                const double y = -grid.halfExtent + (j + 0.5) * h;
                double mesh = 0.0;
                if (!surface.sampleElevation(x, y, mesh)) continue;
                worst = std::max(worst, std::abs(mesh - train.elevation(x, y, 0.0)));
                largestCosine = std::max(
                    largestCosine, std::abs(std::cos(c.wavenumber * x + c.phase)));
            }
        const double predicted =
            amplitude * (1.0 - std::cos(0.5 * c.wavenumber * h)) * largestCosine;
        expectNear(label("the interpolation error is a (1 - cos(k h / 2)) cos(psi) at cells per"
                         " wavelength",
                         cellsPerWavelength[variant]),
                   worst, predicted, 2e-6);
        relative[variant] = worst / amplitude;
    }

    std::printf("     cell-centre error: %.2f%% of amplitude at 16 cells per wavelength,"
                " %.1f%% at 4\n", 100.0 * relative[0], 100.0 * relative[1]);

    // 1 - cos(pi/16) = 1.92%, and the midpoint phases are 22.5 degrees apart so
    // the largest cosine is at least cos(11.25 deg) = 0.981.
    expectTrue(label("16 cells per wavelength reproduces the surface to 2%, at", relative[0]),
               relative[0] < 0.02);
    // 4 cells per wavelength puts the midpoint phases a quarter turn apart, so the
    // largest cosine is at least 1/sqrt(2) and the error is at least
    // (1 - cos(45 deg)) / sqrt(2) = 20.7% of the amplitude. Twenty per cent of the
    // wave height, from a grid that sounds like it ought to be enough.
    expectTrue(label("4 cells per wavelength does not, at", relative[1]), relative[1] > 0.20);
    expectTrue("the coarse grid is more than ten times worse", relative[1] > 10.0 * relative[0]);
}

// The finding that costs the most to ignore: resolving the *dominant* wavelength
// is not the criterion. A spectral sea carries components an order of magnitude
// shorter, and they are what the grid has to resolve.
void testResolvingTheDominantWavelengthIsNotEnough() {
    const double t = 5.0;
    const sim::WaveField field = spectralSea(3.0, 9.0, 48, 12, 0x5eaf00d);
    const double dominant = gpu::dominantWavelength(field);
    const double shortest = gpu::shortestWavelength(field);
    expectNear("the dominant wavelength is g Tp^2 / 2 pi to a percent", dominant,
               sim::kGravity * 81.0 / (2.0 * sim::kPi), 0.01 * 126.4);
    expectTrue(label("the shortest component is an order of magnitude shorter, at", shortest),
               shortest < 0.15 * dominant);

    const double halfExtent = 100.0;
    const auto measure = [&](int cells) {
        const gpu::OceanGrid grid{0.0, 0.0, halfExtent, cells};
        gpu::OceanSurface surface;
        surface.build(field, grid, t);
        const double h = grid.cellSize();
        double worst = 0.0;
        for (int j = 0; j < cells; ++j)
            for (int i = 0; i < cells; ++i) {
                const double x = -halfExtent + (i + 0.5) * h;
                const double y = -halfExtent + (j + 0.5) * h;
                double mesh = 0.0;
                if (!surface.sampleElevation(x, y, mesh)) continue;
                worst = std::max(worst, std::abs(mesh - field.elevation(x, y, t)));
            }
        return worst;
    };

    const int byDominant = gpu::oceanResolutionFor(halfExtent, dominant, 16.0);
    const int byShortest = gpu::oceanResolutionFor(halfExtent, shortest, 8.0);
    const double errorByDominant = measure(byDominant);
    const double errorByShortest = measure(byShortest);

    std::printf("     Hs 3 m Tp 9 s: dominant %.1f m, shortest %.2f m;"
                " %d cells -> %.3f m error, %d cells -> %.3f m\n",
                dominant, shortest, byDominant, errorByDominant, byShortest, errorByShortest);

    expectTrue("resolving the shortest component needs a far finer grid", byShortest > 4 * byDominant);
    // Sixteen cells across the dominant wavelength still invents a large fraction
    // of the wave height, because the short components are below its Nyquist.
    expectTrue(label("16 cells per dominant wavelength is still badly out, m", errorByDominant),
               errorByDominant > 0.10 * 3.0);
    expectTrue(label("8 cells per shortest component is not, m", errorByShortest),
               errorByShortest < 0.03 * 3.0);
    // The closed-form worst case: every component at its own extremum at once.
    double bound = 0.0;
    const double h = 2.0 * halfExtent / byShortest;
    for (const sim::WaveComponent& component : field.components())
        bound += component.amplitude * (1.0 - std::cos(0.5 * component.wavenumber * h));
    expectTrue("the measured error respects the sum of the per-component bounds",
               errorByShortest <= bound);
}

// --- Cost ---------------------------------------------------------------------

// Measured, not extrapolated, and on the shipped path rather than a model of it.
// The numbers land in docs/03-renderer-audio.md.
void testDisplacementCost() {
    const sim::WaveField field = spectralSea(3.0, 9.0, 48, 12, 0x5eaf00d);
    const auto components = static_cast<double>(field.components().size());
    std::printf("     CPU displacement, %g components:\n", components);

    for (int cells : {64, 128, 256}) {
        const gpu::OceanGrid grid{0.0, 0.0, 200.0, cells};
        gpu::OceanSurface surface;
        double best = 1e30;
        for (int repeat = 0; repeat < 2; ++repeat) {
            surface.build(field, grid, 1.0 + repeat);
            best = std::min(best, surface.buildSeconds());
        }
        const double vertices = static_cast<double>(surface.vertices().size());
        std::printf("       %4d cells  %7.0f vertices  %7.2f ms  %5.2f ns per vertex-component\n",
                    cells, vertices, best * 1e3, best * 1e9 / (vertices * components));
        // The structural fact behind the number: one evaluation per unique grid
        // vertex, shared by the six triangle corners that reference it.
        expectEqual(label("one vertex per grid point at cells", cells),
                    static_cast<long long>(surface.vertices().size()),
                    static_cast<long long>(cells + 1) * (cells + 1));
        expectEqual(label("six indices per cell at cells", cells),
                    static_cast<long long>(surface.indices().size()),
                    static_cast<long long>(cells) * cells * 6);
    }
}

}  // namespace

void runOceanTests() {
    std::printf("\n--- ocean surface ---\n");
    testCalmSeaIsExactlyFlat();
    testCalmSeaRendersAFlatPlane();
    testGridElevationMatchesTheWaveField();
    testNormalsAreTheAnalyticSurfaceSlope();
    testRenderedSurfaceAgreesWithTheWaveField();
    testRenderingIsRepeatable();
    testAdvancingTimeChangesTheSea();
    testLongCrestedSeaVariesAlongXOnly();
    testLightingRespondsToTheWaves();
    testGridMustResolveTheWavelength();
    testResolvingTheDominantWavelengthIsNotEnough();
    testDisplacementCost();
}
