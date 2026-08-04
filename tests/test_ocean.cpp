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

// What the comparison below found, aggregated over its samples.
struct ElevationAgreement {
    long long samples = 0, exceeded = 0, staleRejected = 0;
    double mean = 0.0, rms = 0.0, worstError = 0.0, worstBound = 0.0, worstRatio = 0.0;
    double lowest = 1e30, highest = -1e30;
};

// Holds a rendered frame's elevation channel against `WaveField::elevation()` at
// each of `probes`, with the tolerance derived per sample rather than chosen.
//
// Factored out because the cascade has to pass exactly this check in its near
// field, and the tolerance construction is the part that must not be reinvented:
// it was wrong once (see the note in meshFromField below) and the way it was
// wrong is not visible from the outside.
//
// `meshFromField` is the piecewise-linear surface a *correct* mesh of the
// relevant spacing would have, computed from the field. It must never come from
// the artefact under test.
template <typename MeshFromField>
ElevationAgreement compareRenderedElevation(const Image& image, const sim::Mat4& mvp,
                                            const gpu::OceanView& view,
                                            const sim::WaveField& field, double t,
                                            const std::vector<sim::Vec3>& probes,
                                            MeshFromField meshFromField) {
    ElevationAgreement out;
    const double quantisation = 0.5 * static_cast<double>(view.elevationSpan) / 65535.0;
    const double step = 0.05;
    double sumError = 0.0, sumSquares = 0.0;

    for (const sim::Vec3& point : probes) {
        double clip[4], px = 0, py = 0;
        mvp.transform(point, clip);
        if (!sim::clipToPixel(clip, image.width, image.height, px, py)) continue;
        if (px < 1.0 || py < 1.0 || px > image.width - 2.0 || py > image.height - 2.0) continue;

        const auto column = static_cast<std::uint32_t>(px);
        const auto row = static_cast<std::uint32_t>(py);
        double rendered = 0.0;
        if (!gpu::decodeOceanElevation(view, image.pixel(column, row), rendered)) continue;

        // The Jacobian of the projection here, from two short probes through the
        // same matrix. Columns are d(pixel)/dx and d(pixel)/dy.
        double clipX[4], clipY[4], ax = 0, ay = 0, bx = 0, by = 0;
        mvp.transform({point.x + step, point.y, point.z}, clipX);
        mvp.transform({point.x, point.y + step, point.z}, clipY);
        if (!sim::clipToPixel(clipX, image.width, image.height, ax, ay)) continue;
        if (!sim::clipToPixel(clipY, image.width, image.height, bx, by)) continue;
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
        // Newton rather than one linear step: the surface height changes as the
        // solution moves, and that shifts the projection again. Ignoring it left
        // four samples of four hundred outside their bound -- a real 2 mm effect,
        // not noise, and the sort of residual that would have been quietly
        // absorbed by rounding the tolerance up.
        double dx = 0.0, dy = 0.0;
        for (int iteration = 0; iteration < 4; ++iteration) {
            const double z = field.elevation(point.x + dx, point.y + dy, t);
            double clipHere[4], qx = 0, qy = 0;
            mvp.transform({point.x + dx, point.y + dy, z}, clipHere);
            if (!sim::clipToPixel(clipHere, image.width, image.height, qx, qy)) break;
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
        // 1.1; the rest is that the interpolation error vanishes at the vertices,
        // so evaluating it a few millimetres from where the ray actually lands
        // changes it by a large *relative* amount while staying in millimetres
        // absolutely. Measured worst ratio is printed, so the margin is visible
        // rather than assumed -- and a 2% elevation error, which puts 47 mm on the
        // board, is still rejected by an order of magnitude.
        const double interpolation = std::abs(meshFromField(hitX, hitY) - analytic);
        const double bound = quantisation + 3.0 * interpolation;

        ++out.samples;
        out.worstError = std::max(out.worstError, std::abs(error));
        out.worstBound = std::max(out.worstBound, bound);
        if (interpolation > 1e-9)
            out.worstRatio =
                std::max(out.worstRatio, (std::abs(error) - quantisation) / interpolation);
        if (std::abs(error) > bound) ++out.exceeded;
        sumError += error;
        sumSquares += error * error;
        out.lowest = std::min(out.lowest, analytic);
        out.highest = std::max(out.highest, analytic);

        // Negative control. One second is a small, entirely plausible error -- a
        // renderer a frame or two behind the physics -- and the same comparison
        // must reject it. Without this the tolerance above could be loose enough
        // to accept anything.
        if (std::abs(rendered - field.elevation(hitX, hitY, t + 1.0)) > bound) ++out.staleRejected;
    }

    if (out.samples > 0) {
        out.mean = sumError / static_cast<double>(out.samples);
        out.rms = std::sqrt(sumSquares / static_cast<double>(out.samples));
    }
    return out;
}

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

    const double h = grid.cellSize();

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

    std::vector<sim::Vec3> probes;
    for (int j = 6; j < grid.resolution; j += 11)
        for (int i = 6; i < grid.resolution; i += 11) {
            const gpu::OceanVertex& v = surface.vertices()[surface.vertexIndex(i, j)];
            probes.push_back({v.position[0], v.position[1], v.position[2]});
        }

    const ElevationAgreement found =
        compareRenderedElevation(image, mvp, view, field, t, probes, meshFromField);

    expectTrue("the check had samples to work with", found.samples > 200);
    if (found.samples == 0) return;

    std::printf("     rendered vs WaveField::elevation: %lld samples, mean %+.5f m, rms %.5f m,"
                " worst %.5f m against a worst bound of %.5f m (ratio needed %.2f)\n",
                found.samples, found.mean, found.rms, found.worstError, found.worstBound,
                found.worstRatio);

    expectEqual("every rendered point agrees with WaveField::elevation inside its own bound",
                found.exceeded, 0);
    // A residual error is as likely one way as the other; a sign flip, a stale
    // time, a scale factor or a shifted patch would not be.
    expectNear("the rendered surface carries no systematic offset", found.mean, 0.0, 0.002);
    // Guards against a vacuous pass.
    expectTrue(label("the sampled surface varies by", found.highest - found.lowest),
               found.highest - found.lowest > 2.0);
    expectTrue("the agreement is three orders tighter than the variation being measured",
               found.rms * 1000.0 < found.highest - found.lowest);
    expectTrue(label("a one-second-stale field is rejected, fraction",
                     static_cast<double>(found.staleRejected) / static_cast<double>(found.samples)),
               found.staleRejected > found.samples * 8 / 10);
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

// --- The cascade ---------------------------------------------------------------

// Every directed edge of a triangle soup, counted. An interior edge of a
// watertight surface appears exactly twice, once in each direction; an edge on
// the outer boundary appears once. **A T-junction is exactly an interior edge
// that appears once**, so this is the crack question asked combinatorially --
// no camera, no pixels, no tolerance -- and it is the same instrument that caught
// the hull wound inconsistently (CLAUDE.md).
struct EdgeCensus {
    long long interior = 0;      // seen twice, opposite directions
    long long unmatched = 0;     // seen once: a boundary edge or a crack
    long long malformed = 0;     // seen twice the same way, or more than twice
    long long onOuterEdge = 0;   // of the unmatched, those on the cascade's outer square
    long long degenerate = 0;    // an edge whose two ends are the same vertex
};

EdgeCensus censusEdges(const std::vector<gpu::OceanVertex>& vertices,
                       const std::vector<std::uint32_t>& indices, double outerHalfExtent,
                       double centreX, double centreY) {
    EdgeCensus census;
    // key = (low << 32) | high, value = (count, net direction).
    std::vector<std::pair<std::uint64_t, int>> edges;
    edges.reserve(indices.size());
    for (std::size_t t = 0; t + 2 < indices.size(); t += 3)
        for (int e = 0; e < 3; ++e) {
            const std::uint32_t a = indices[t + static_cast<std::size_t>(e)];
            const std::uint32_t b = indices[t + static_cast<std::size_t>((e + 1) % 3)];
            if (a == b) { ++census.degenerate; continue; }
            const std::uint64_t low = std::min(a, b), high = std::max(a, b);
            edges.emplace_back((low << 32) | high, a < b ? 1 : -1);
        }
    std::sort(edges.begin(), edges.end(),
              [](const auto& a, const auto& b) { return a.first < b.first; });

    const auto onOuter = [&](std::uint32_t v) {
        const gpu::OceanVertex& p = vertices[v];
        const double radius = std::max(std::abs(static_cast<double>(p.position[0]) - centreX),
                                       std::abs(static_cast<double>(p.position[1]) - centreY));
        return radius > outerHalfExtent * (1.0 - 1e-9);
    };

    std::size_t i = 0;
    while (i < edges.size()) {
        std::size_t j = i;
        int net = 0;
        while (j < edges.size() && edges[j].first == edges[i].first) net += edges[j++].second;
        const std::size_t count = j - i;
        if (count == 1) {
            ++census.unmatched;
            const auto low = static_cast<std::uint32_t>(edges[i].first >> 32);
            const auto high = static_cast<std::uint32_t>(edges[i].first & 0xffffffffu);
            if (onOuter(low) && onOuter(high)) ++census.onOuterEdge;
        } else if (count == 2 && net == 0) {
            ++census.interior;
        } else {
            ++census.malformed;
        }
        i = j;
    }
    return census;
}

// The structural claim, before any picture: the cascade is one watertight sheet
// whose only boundary is its outer square, and its vertex and triangle counts are
// the closed forms the ring construction predicts.
void testCascadeIsOneWatertightSheet() {
    const sim::WaveField field = spectralSea(3.0, 9.0, 24, 6, 0xCA5CADE);
    gpu::OceanCascade cascade;
    cascade.centreX = 40.0;   // off-centre, so a cascade that assumed the origin shows up
    cascade.centreY = -25.0;
    cascade.innerHalfExtent = 64.0;
    cascade.resolution = 32;
    cascade.levels = 4;

    gpu::OceanCascadeSurface surface;
    surface.build(field, cascade, 3.5);

    const double n = cascade.resolution;
    const auto levels = static_cast<double>(cascade.levels);
    // Level 0 is a full grid; a ring is the level's square less the quarter its
    // predecessor already covers, stored without the hole or its boundary.
    expectEqual("the cascade has (n + 1)^2 + (levels - 1)(3n^2/4 + n) vertices",
                static_cast<long long>(surface.vertices().size()),
                static_cast<long long>((n + 1) * (n + 1) + (levels - 1) * (0.75 * n * n + n)));
    // Two triangles a cell, except along a seam, where 2n cells carry three.
    expectEqual("and 2n^2 + (levels - 1)(3n^2/2 + 2n) triangles",
                static_cast<long long>(surface.indices().size() / 3),
                static_cast<long long>(2 * n * n + (levels - 1) * (1.5 * n * n + 2 * n)));
    expectEqual("the index buffer is whole triangles",
                static_cast<long long>(surface.indices().size() % 3), 0);
    // And per level, so the reported breakdown is the geometry rather than a
    // second derivation of it that could drift from the first.
    long long reported = 0;
    for (int level = 0; level < cascade.levels; ++level) {
        const gpu::OceanCascadeSurface::Level& info =
            surface.levels()[static_cast<std::size_t>(level)];
        expectEqual(label("level's triangle count at half extent", info.halfExtent),
                    static_cast<long long>(info.triangles),
                    level == 0 ? static_cast<long long>(2 * n * n)
                               : static_cast<long long>(1.5 * n * n + 2 * n));
        reported += static_cast<long long>(info.triangles);
        expectEqual(label("level's vertex count at half extent", info.halfExtent),
                    static_cast<long long>(info.vertices),
                    level == 0 ? static_cast<long long>((n + 1) * (n + 1))
                               : static_cast<long long>(0.75 * n * n + n));
    }
    expectEqual("the levels account for every triangle in the buffer", reported,
                static_cast<long long>(surface.indices().size() / 3));

    const EdgeCensus census = censusEdges(surface.vertices(), surface.indices(), cascade.reach(),
                                          cascade.centreX, cascade.centreY);
    std::printf("     cascade edges: %lld interior, %lld unmatched (%lld on the outer square),"
                " %lld malformed\n",
                census.interior, census.unmatched, census.onOuterEdge, census.malformed);

    expectEqual("no triangle of the cascade is degenerate", census.degenerate, 0);
    expectEqual("every edge is used once or twice and never twice the same way",
                census.malformed, 0);
    // 4n edges around the outermost level's square, and nothing else may be open.
    expectEqual("the only unmatched edges are the outer boundary", census.unmatched,
                static_cast<long long>(4 * n));
    expectEqual("and every one of them really is on the outer square", census.onOuterEdge,
                census.unmatched);
    // Guard against a vacuous pass: a mesh with no interior edges would satisfy
    // the first two trivially.
    expectTrue(label("there were interior edges to check, count",
                     static_cast<double>(census.interior)),
               census.interior > 1000);

    // The negative control. Without the stitch every seam leaves T-junctions:
    // each of the 2n coarse cells along a seam has one edge spanning two of the
    // finer level's, so three edges go unmatched per seam cell.
    gpu::OceanCascade cracked = cascade;
    cracked.stitchSeams = false;
    gpu::OceanCascadeSurface unstitched;
    unstitched.build(field, cracked, 3.5);
    const EdgeCensus broken = censusEdges(unstitched.vertices(), unstitched.indices(),
                                          cascade.reach(), cascade.centreX, cascade.centreY);
    std::printf("     without the stitch: %lld unmatched edges (%lld on the outer square)\n",
                broken.unmatched, broken.onOuterEdge);
    expectEqual("an unstitched cascade leaves 6n T-junction edges at each of its seams",
                broken.unmatched,
                static_cast<long long>(4 * n + 6 * n * (levels - 1)));
    expectTrue("which is what the census is able to see", broken.unmatched > census.unmatched);
}

// The ring construction needs the cell count to be a multiple of four, and a
// caller has no reason to know that. So the cascade rounds -- and it has to round
// *up*, because `resolution` is the answer to a correctness question and handing
// back a coarser grid than was asked for is exactly the silent degradation
// `docs/03-renderer-audio.md` says resolution is not allowed to be.
//
// Swept over every count from one to forty rather than the three anyone writes by
// hand, and each one is asked whether it is still one watertight sheet: an odd
// count, a count of one, and a count that is two mod four all take different paths
// through the hole arithmetic.
void testEveryResolutionRoundsUpAndStaysWatertight() {
    const sim::WaveField field = spectralSea(2.0, 7.0, 8, 4, 0x4444);
    bool roundedUp = true, multipleOfFour = true, watertight = true, counted = true;
    long long swept = 0;
    for (int asked = 1; asked <= 40; ++asked) {
        gpu::OceanCascade cascade;
        cascade.innerHalfExtent = 30.0;
        cascade.resolution = asked;
        cascade.levels = 3;
        const int used = cascade.cellsPerSide();
        roundedUp = roundedUp && used >= asked && used >= 4;
        multipleOfFour = multipleOfFour && used % 4 == 0;
        // The smallest multiple of four that is at least `asked`, and no larger.
        roundedUp = roundedUp && used - 4 < std::max(asked, 4);

        gpu::OceanCascadeSurface surface;
        surface.build(field, cascade, 0.5);
        const EdgeCensus census = censusEdges(surface.vertices(), surface.indices(),
                                              cascade.reach(), 0.0, 0.0);
        watertight = watertight && census.malformed == 0 && census.degenerate == 0 &&
                     census.unmatched == 4LL * used && census.onOuterEdge == census.unmatched;
        const auto n = static_cast<double>(used);
        counted = counted &&
                  surface.vertices().size() ==
                      static_cast<std::size_t>((n + 1) * (n + 1) + 2 * (0.75 * n * n + n));
        ++swept;
    }
    expectEqual("every resolution from 1 to 40 was swept", swept, 40);
    expectTrue("the cascade never quietly hands back a coarser grid than asked for", roundedUp);
    expectTrue("and always a multiple of four", multipleOfFour);
    expectTrue("every one of them is one watertight sheet", watertight);
    expectTrue("and has the vertex count the layout predicts", counted);
}

// The seam vertices are not two evaluations that happen to agree -- they are one
// entry in one array, addressed from both sides.
void testCascadeSeamVerticesAreShared() {
    const sim::WaveField field = spectralSea(3.0, 9.0, 24, 6, 0xC0FFEE);
    gpu::OceanCascade cascade;
    cascade.innerHalfExtent = 50.0;
    cascade.resolution = 32;
    cascade.levels = 3;
    gpu::OceanCascadeSurface surface;
    surface.build(field, cascade, 1.25);

    const int n = cascade.resolution;
    bool shared = true, positioned = true;
    long long checked = 0;
    for (int level = 1; level < cascade.levels; ++level)
        for (int j = -n / 4; j <= n / 4; ++j)
            for (int i = -n / 4; i <= n / 4; ++i) {
                if (std::max(std::abs(i), std::abs(j)) != n / 4) continue;  // the hole's boundary
                shared = shared && surface.vertexIndexAt(level, i, j) ==
                                       surface.vertexIndexAt(level - 1, 2 * i, 2 * j);
                // And the two levels agree about where that vertex is in the
                // world, which is the claim that makes sharing legitimate.
                const gpu::OceanVertex& v =
                    surface.vertices()[surface.vertexIndexAt(level, i, j)];
                positioned = positioned &&
                             std::abs(v.position[0] - i * cascade.cellSize(level)) < 1e-3 &&
                             std::abs(v.position[1] - j * cascade.cellSize(level)) < 1e-3;
                ++checked;
            }
    expectTrue(label("a seam vertex is one vertex, addressed from both levels, checked",
                     static_cast<double>(checked)),
               shared);
    expectTrue("and it stands where both levels say it stands", positioned);
    expectTrue("there were seam vertices to check", checked > 100);

    // Every vertex is referenced by at least one triangle: a level that stored
    // its hole would still pass the counts above and waste a quarter of the
    // displacement on water nobody draws.
    std::vector<bool> used(surface.vertices().size(), false);
    for (std::uint32_t index : surface.indices()) used[index] = true;
    const auto orphans = static_cast<long long>(std::count(used.begin(), used.end(), false));
    expectEqual("no vertex of the cascade is displaced and then never drawn", orphans, 0);
}

// Both builders keep their index buffer across calls and rebuild it only when the
// topology changes, which is worth doing and is a standing way to be wrong: a
// cache key that misses one of the things it depends on hands back the previous
// frame's triangles over this frame's vertices, and the result is still a
// plausible sea. So every reuse is held against a surface built from nothing.
//
// This exists because a mutation that dropped `stitchSeams` from the cascade's
// cache key was caught by no assertion in this file: every other test builds a
// fresh object for each configuration, so the cache was never asked to notice a
// change.
void testRebuildingCarriesNothingStale() {
    const sim::WaveField field = spectralSea(3.0, 8.0, 16, 6, 0xB0A7);

    gpu::OceanCascade a;
    a.innerHalfExtent = 40.0;
    a.resolution = 16;
    a.levels = 2;
    a.stitchSeams = true;

    // Each of these differs from `a` in exactly one thing the index buffer
    // depends on, and each is visited after a build of `a`.
    std::vector<gpu::OceanCascade> variants;
    gpu::OceanCascade differentLevels = a; differentLevels.levels = 4;
    gpu::OceanCascade differentCells = a; differentCells.resolution = 24;
    gpu::OceanCascade unstitched = a; unstitched.stitchSeams = false;
    gpu::OceanCascade moved = a; moved.centreX = 500.0;
    variants = {differentLevels, differentCells, unstitched, moved, a};

    gpu::OceanCascadeSurface reused;
    long long identical = 0;
    for (const gpu::OceanCascade& variant : variants) {
        reused.build(field, a, 2.0);
        reused.build(field, variant, 5.5);
        gpu::OceanCascadeSurface fresh;
        fresh.build(field, variant, 5.5);
        const bool sameIndices = reused.indices() == fresh.indices();
        const bool sameVertices =
            reused.vertices().size() == fresh.vertices().size() &&
            std::memcmp(reused.vertices().data(), fresh.vertices().data(),
                        fresh.vertices().size() * sizeof(gpu::OceanVertex)) == 0;
        expectTrue(label("a rebuilt cascade is the cascade, at levels",
                         static_cast<double>(variant.levels)),
                   sameIndices && sameVertices);
        if (sameIndices && sameVertices) ++identical;
    }
    expectEqual("every variant was compared", identical,
                static_cast<long long>(variants.size()));

    // Guard: the comparison has to be able to fail. The variants must not all
    // produce the same geometry as `a`.
    gpu::OceanCascadeSurface base;
    base.build(field, a, 5.5);
    long long distinct = 0;
    for (const gpu::OceanCascade& variant : variants) {
        gpu::OceanCascadeSurface fresh;
        fresh.build(field, variant, 5.5);
        if (fresh.indices() != base.indices() ||
            fresh.vertices().size() != base.vertices().size() ||
            std::memcmp(fresh.vertices().data(), base.vertices().data(),
                        base.vertices().size() * sizeof(gpu::OceanVertex)) != 0)
            ++distinct;
    }
    expectEqual("and four of the five variants really are different geometry", distinct,
                static_cast<long long>(variants.size()) - 1);

    // The uniform patch caches the same way and gets the same question.
    gpu::OceanGrid coarse{0.0, 0.0, 40.0, 8};
    gpu::OceanGrid fine{0.0, 0.0, 40.0, 20};
    gpu::OceanSurface recycled, once;
    recycled.build(field, coarse, 1.0);
    recycled.build(field, fine, 1.0);
    once.build(field, fine, 1.0);
    expectTrue("a rebuilt patch is the patch", recycled.indices() == once.indices() &&
                                                   recycled.vertices().size() ==
                                                       once.vertices().size());
}

// The seam has to be geometrically invisible as well as combinatorially closed:
// the mesh's own elevation, approached from the fine side and from the coarse
// side, must be the same number.
void testCascadeElevationAgreesAcrossALevelBoundary() {
    const double t = 8.75;
    // Short waves, so the surface is genuinely curved at the seam and a straight
    // coarse edge would be visibly wrong there.
    const sim::WaveField field = spectralSea(4.0, 6.0, 32, 8, 0x5EA5EA);
    gpu::OceanCascade cascade;
    cascade.innerHalfExtent = 60.0;
    cascade.resolution = 64;
    cascade.levels = 3;
    // Nothing dropped anywhere, so a step at the seam could only be geometry.
    cascade.minimumCellsPerWavelength = 0.5;

    gpu::OceanCascadeSurface surface;
    surface.build(field, cascade, t);
    gpu::OceanCascade cracked = cascade;
    cracked.stitchSeams = false;
    gpu::OceanCascadeSurface unstitched;
    unstitched.build(field, cracked, t);

    const double seam = cascade.innerHalfExtent;         // the level 0 / level 1 boundary
    const double coarse = cascade.cellSize(1);           // the coarse cell the seam is split into
    const double epsilon = 1e-4;
    // The two samples straddle the seam by 2 * epsilon, so a perfectly continuous
    // surface still reports that much times the steepest slope the field can
    // reach. Derived rather than chosen: sum a k over the components.
    double steepest = 0.0;
    for (const sim::WaveComponent& c : field.components()) steepest += c.amplitude * c.wavenumber;
    const double continuous = 2.0 * epsilon * steepest;

    double worstStep = 0.0, worstCrack = 0.0, largestSag = 0.0;
    long long compared = 0;
    // Walk the midpoints of the coarse cells along the seam. Those are exactly
    // the points a T-junction gets wrong: the fine level has a vertex there and
    // an unstitched coarse cell interpolates straight past it.
    for (int cell = -cascade.resolution / 4; cell < cascade.resolution / 4; ++cell) {
        const double x = (cell + 0.5) * coarse;
        double inside = 0.0, outside = 0.0;
        if (!surface.sampleElevation(x, seam - epsilon, inside)) continue;
        if (!surface.sampleElevation(x, seam + epsilon, outside)) continue;
        worstStep = std::max(worstStep, std::abs(inside - outside));

        double brokenInside = 0.0, brokenOutside = 0.0;
        if (unstitched.sampleElevation(x, seam - epsilon, brokenInside) &&
            unstitched.sampleElevation(x, seam + epsilon, brokenOutside))
            worstCrack = std::max(worstCrack, std::abs(brokenInside - brokenOutside));

        // How far the fine level's midpoint is off the straight coarse edge --
        // the sagitta, which is the size of the gap the stitch closes. Taken
        // from the field, not from the mesh.
        const double a = field.elevation(x - 0.5 * coarse, seam, t);
        const double b = field.elevation(x + 0.5 * coarse, seam, t);
        largestSag = std::max(largestSag, std::abs(field.elevation(x, seam, t) - 0.5 * (a + b)));
        ++compared;
    }

    std::printf("     seam: worst step %.2e m stitched (a continuous surface owes %.2e m),"
                " %.4f m unstitched, sagitta up to %.4f m over %lld coarse cells\n",
                worstStep, continuous, worstCrack, largestSag, compared);

    expectTrue("there were seam crossings to compare", compared > 10);
    // All that is allowed is the 2 * epsilon the samples are apart, times the
    // slope. There is no room in that for a step.
    expectTrue(label("the stitched seam has no step, worst", worstStep), worstStep <= continuous);
    // Guard against a vacuous pass: a flat sea has no step either. The gap the
    // stitch is closing has to be large compared with what is left.
    expectTrue(label("the surface really does bend across the seam, sagitta", largestSag),
               largestSag > 0.02);
    expectTrue("and the unstitched cascade opens most of that sagitta as a step",
               worstCrack > 0.5 * largestSag);
    expectTrue("so the stitched seam is three orders of magnitude tighter",
               worstStep * 1000.0 < worstCrack);
    // The sampler is not blind to a step: on the unstitched cascade the same two
    // probes, a fifth of a millimetre apart, report a third of a metre.
    expectTrue(label("the unstitched step is far outside what continuity allows, m", worstCrack),
               worstCrack > 100.0 * continuous);

    // And the sampler answers *on* the boundaries as well as either side of them.
    // A level chosen with < rather than <= sends a point standing exactly on a
    // seam to the level outside it, and a point on the outer edge to no level at
    // all; both are silent, and the epsilon probes above step right over them.
    long long onSeam = 0, missed = 0;
    double worstBoundary = 0.0;
    for (int level = 0; level < cascade.levels; ++level) {
        const double edge = cascade.halfExtent(level);
        for (int step = -12; step <= 12; ++step) {
            const double x = step * 0.077 * edge;
            double here = 0.0, inside = 0.0;
            if (!surface.sampleElevation(x, edge, here)) { ++missed; continue; }
            if (!surface.sampleElevation(x, edge - epsilon, inside)) { ++missed; continue; }
            worstBoundary = std::max(worstBoundary, std::abs(here - inside));
            ++onSeam;
        }
    }
    expectEqual("every point standing exactly on a level boundary is inside the cascade", missed,
                0);
    expectTrue("there were boundary points to sample", onSeam > 50);
    expectTrue(label("and the surface there is the one just inside it, worst", worstBoundary),
               worstBoundary <= continuous);
}

// A level carries a component only when it can reproduce it, and the geometry
// says so -- not just the bookkeeping.
void testCascadeDropsWhatItCannotResolve() {
    const double t = 2.0;
    const sim::WaveField field = spectralSea(3.0, 9.0, 48, 12, 0x5eaf00d);
    const double shortest = gpu::shortestWavelength(field);

    gpu::OceanCascade cascade;
    // Level 0 cells of 2 m, which is under a quarter of the 11.4 m shortest
    // component, so the finest level carries the whole spectrum and every
    // component that goes missing further out went missing for a reason.
    cascade.innerHalfExtent = 64.0;
    cascade.resolution = 64;
    cascade.levels = 8;  // out to 8.2 km, cells of 256 m
    cascade.minimumCellsPerWavelength = 4.0;
    expectTrue(label("level 0 resolves the shortest component, cell size", cascade.cellSize(0)),
               cascade.cellSize(0) * cascade.minimumCellsPerWavelength < shortest);
    gpu::OceanCascadeSurface surface;
    surface.build(field, cascade, t);

    std::size_t previous = field.components().size() + 1;
    for (int level = 0; level < cascade.levels; ++level) {
        const gpu::OceanCascadeSurface::Level& info = surface.levels()[static_cast<std::size_t>(level)];
        // Counted independently, from the field.
        std::size_t want = 0;
        for (const sim::WaveComponent& c : field.components())
            if (2.0 * sim::kPi / c.wavenumber >=
                cascade.minimumCellsPerWavelength * info.cellSize)
                ++want;
        expectEqual(label("the level carries every component it resolves and no others, at cell",
                          info.cellSize),
                    static_cast<long long>(info.components), static_cast<long long>(want));
        expectTrue("a coarser level never carries more than a finer one",
                   info.components <= previous);
        previous = info.components;
    }
    std::printf("     cascade at %g m cells: components per level", cascade.cellSize(0));
    for (const gpu::OceanCascadeSurface::Level& info : surface.levels())
        std::printf(" %zu", info.components);
    std::printf(" (field has %zu, shortest %.1f m)\n", field.components().size(), shortest);

    expectEqual("the finest level carries the whole spectrum",
                static_cast<long long>(surface.levels().front().components),
                static_cast<long long>(field.components().size()));
    expectEqual("the outermost level carries none of it",
                static_cast<long long>(surface.levels().back().components), 0);

    // A level that resolves nothing skips the recurrence entirely, which is what
    // makes the far rings nearly free -- and is a fast path with its own way of
    // being wrong. It shares the row accumulators with every other level, so a
    // version that forgot to write zeros would quietly emit the *previous* row's
    // heights and still render a plausible sea. Bit-exact, therefore, like the
    // calm sea's: not "near flat".
    long long flatVertices = 0;
    bool exactlyFlat = true, exactlyUpright = true;
    for (int level = 0; level < cascade.levels; ++level) {
        if (surface.levels()[static_cast<std::size_t>(level)].components != 0) continue;
        const int half = cascade.resolution / 2;
        for (int j = -half; j <= half; ++j)
            for (int i = -half; i <= half; ++i) {
                if (std::max(std::abs(i), std::abs(j)) <= cascade.resolution / 4) continue;
                const gpu::OceanVertex& v = surface.vertices()[surface.vertexIndexAt(level, i, j)];
                exactlyFlat = exactlyFlat && v.position[2] == 0.0f;
                exactlyUpright = exactlyUpright && v.normal[0] == 0.0f && v.normal[1] == 0.0f &&
                                 v.normal[2] == 1.0f;
                ++flatVertices;
            }
    }
    expectTrue(label("a level that resolves nothing is exactly flat, vertices checked",
                     static_cast<double>(flatVertices)),
               exactlyFlat);
    expectTrue("and its normals are exactly +z", exactlyUpright);
    expectTrue("there were flat levels to check", flatVertices > 1000);

    // The geometric consequence, which is the part that matters: a partial level's
    // vertices carry the band-limited sum and not the full one.
    // The *outermost* level that still carries something, because that is the one
    // that has dropped the most: the first partial level drops only a dozen of the
    // shortest components and their combined amplitude is a centimetre, which
    // would make the comparison below true without being able to fail.
    int partial = -1;
    for (int level = 0; level < cascade.levels; ++level) {
        const std::size_t count = surface.levels()[static_cast<std::size_t>(level)].components;
        if (count > 0 && count < field.components().size()) partial = level;
    }
    expectTrue("there is a level that drops some but not all of the spectrum", partial > 0);
    if (partial < 0) return;

    const double cell = cascade.cellSize(partial);
    const int lattice = cascade.resolution / 3;  // inside the ring, off the seam
    double worstBandLimited = 0.0, worstAgainstFull = 1e30;
    for (int j = -lattice; j <= lattice; j += 7)
        for (int i = lattice; i <= cascade.resolution / 2; i += 7) {
            const gpu::OceanVertex& v = surface.vertices()[surface.vertexIndexAt(partial, i, j)];
            const double x = i * cell, y = j * cell;
            double bandLimited = 0.0;
            for (const sim::WaveComponent& c : field.components())
                if (2.0 * sim::kPi / c.wavenumber >= cascade.minimumCellsPerWavelength * cell)
                    bandLimited += c.amplitude * std::cos(c.wavenumber * (c.dirX * x + c.dirY * y) -
                                                          c.omega * t + c.phase);
            worstBandLimited =
                std::max(worstBandLimited, std::abs(v.position[2] - bandLimited));
            worstAgainstFull =
                std::min(worstAgainstFull, std::abs(v.position[2] - field.elevation(x, y, t)));
        }
    expectTrue(label("a ring's vertices are the band-limited field, worst", worstBandLimited),
               worstBandLimited < 1e-6);
    // Guard: if the dropped components had no amplitude the check above would
    // pass against the full field too, and prove nothing.
    expectTrue(label("and are measurably not the full field, closest", worstAgainstFull),
               worstAgainstFull > 0.05);
}

// The lighting answers to the waves at every level, not just at the finest one.
//
// This exists because a mutation that replaced the cascade's normals with a flat
// +z was caught by nothing: `docs/03-renderer-audio.md` records that exact defect
// escaping once already on the uniform patch, and the cascade reintroduced the
// hole by having its own displacement loop.
void testCascadeNormalsAreTheBandLimitedSlope() {
    const double t = 6.5;
    const sim::WaveField field = spectralSea(3.0, 9.0, 24, 8, 0x510BE);
    gpu::OceanCascade cascade;
    cascade.innerHalfExtent = 64.0;
    cascade.resolution = 64;
    cascade.levels = 6;
    gpu::OceanCascadeSurface surface;
    surface.build(field, cascade, t);

    // Level 0 carries the whole spectrum, so its slope is the field's own,
    // checked against a central difference -- a different computation from the
    // sin() the builder accumulates -- with the tolerance from the difference's
    // truncation error rather than picked.
    const double step = 0.01;
    double thirdDerivative = 0.0;
    for (const sim::WaveComponent& c : field.components())
        thirdDerivative += c.amplitude * c.wavenumber * c.wavenumber * c.wavenumber;
    const double tolerance = step * step / 6.0 * thirdDerivative + 1e-6;

    const int half = cascade.resolution / 2;
    const double h0 = cascade.cellSize(0);
    double worstFine = 0.0, steepestFine = 0.0;
    for (int j = -half; j <= half; j += 5)
        for (int i = -half; i <= half; i += 5) {
            const gpu::OceanVertex& v = surface.vertices()[surface.vertexIndexAt(0, i, j)];
            const double x = i * h0, y = j * h0;
            const double gx = (field.elevation(x + step, y, t) - field.elevation(x - step, y, t)) /
                              (2.0 * step);
            const double gy = (field.elevation(x, y + step, t) - field.elevation(x, y - step, t)) /
                              (2.0 * step);
            const double inverse = 1.0 / std::sqrt(gx * gx + gy * gy + 1.0);
            const double want[3] = {-gx * inverse, -gy * inverse, inverse};
            for (int c = 0; c < 3; ++c)
                worstFine = std::max(worstFine, std::abs(static_cast<double>(v.normal[c]) - want[c]));
            steepestFine = std::max(steepestFine, std::hypot(gx, gy));
        }
    expectTrue(label("level 0's normals are the field's own slope, worst", worstFine),
               worstFine < tolerance);
    expectTrue(label("and they actually tilt, steepest slope", steepestFine), steepestFine > 0.05);

    // A ring's normal must be the slope of what the ring *carries* -- not of the
    // full field, which it does not draw, and not of nothing.
    int partial = -1;
    for (int level = 0; level < cascade.levels; ++level) {
        const std::size_t count = surface.levels()[static_cast<std::size_t>(level)].components;
        if (count > 0 && count < field.components().size()) partial = level;
    }
    expectTrue("there is a partly band-limited ring to check", partial > 0);
    if (partial < 0) return;

    const double cell = cascade.cellSize(partial);
    const double keep = cascade.minimumCellsPerWavelength * cell;
    double worstRing = 0.0, steepestRing = 0.0, againstFull = 0.0;
    for (int j = -half; j <= half; j += 5)
        for (int i = half / 2 + 1; i <= half; i += 5) {
            const gpu::OceanVertex& v = surface.vertices()[surface.vertexIndexAt(partial, i, j)];
            const double x = i * cell, y = j * cell;
            double sx = 0.0, sy = 0.0, fullX = 0.0, fullY = 0.0;
            for (const sim::WaveComponent& c : field.components()) {
                const double psi =
                    c.wavenumber * (c.dirX * x + c.dirY * y) - c.omega * t + c.phase;
                const double gain = -c.amplitude * c.wavenumber * std::sin(psi);
                fullX += gain * c.dirX;
                fullY += gain * c.dirY;
                if (2.0 * sim::kPi / c.wavenumber < keep) continue;
                sx += gain * c.dirX;
                sy += gain * c.dirY;
            }
            const double inverse = 1.0 / std::sqrt(sx * sx + sy * sy + 1.0);
            const double want[3] = {-sx * inverse, -sy * inverse, inverse};
            for (int c = 0; c < 3; ++c)
                worstRing = std::max(worstRing, std::abs(static_cast<double>(v.normal[c]) - want[c]));
            const double fullInverse = 1.0 / std::sqrt(fullX * fullX + fullY * fullY + 1.0);
            const double fullNormal[3] = {-fullX * fullInverse, -fullY * fullInverse, fullInverse};
            for (int c = 0; c < 3; ++c)
                againstFull = std::max(
                    againstFull, std::abs(static_cast<double>(v.normal[c]) - fullNormal[c]));
            steepestRing = std::max(steepestRing, std::hypot(sx, sy));
        }
    std::printf("     normals: level 0 worst %.2e (tolerance %.2e), ring %d worst %.2e,"
                " steepest ring slope %.3f\n",
                worstFine, tolerance, partial, worstRing, steepestRing);
    expectTrue(label("a ring's normals are the slope of what it carries, worst", worstRing),
               worstRing < 1e-6);
    // Guards: flat +z would satisfy nothing here only if the ring's slope is
    // real, and the band-limited answer has to be distinguishable from the full
    // one or the check is not about band limiting at all.
    expectTrue(label("a ring's surface really is sloped, steepest", steepestRing),
               steepestRing > 0.02);
    // Five orders of magnitude between the two, so the check above is about the
    // band limiting and not merely about the sea being a sea.
    expectTrue(label("and it is measurably not the full field's slope, worst", againstFull),
               againstFull > 0.01);
}

// How far is far enough, in terms of the camera and nothing else.
void testCascadeReachesTheHorizon() {
    const double eyeHeight = 40.0;
    const double fov = 50.0 * sim::kDegToRad;
    const int pixels = 720;
    const double reach = gpu::oceanHorizonReach(eyeHeight, fov, pixels);
    expectNear("the horizon reach is eyeHeight * pixels / (2 tan(fov / 2))", reach,
               eyeHeight * pixels / (2.0 * std::tan(0.5 * fov)), 1e-9 * reach);

    // Asserted through a real camera, not against its own formula: put the eye at
    // that height, look along +x at the horizon, and ask how many pixels separate
    // the patch edge from the vanishing point of the sea plane.
    const sim::Mat4 mvp =
        sim::perspective(fov, 1.0, 1.0, 4.0 * reach) *
        sim::lookAt({0, 0, eyeHeight}, {1, 0, eyeHeight}, {0, 0, 1});
    // The horizon is the vanishing point of a horizontal direction: the point at
    // infinity, w = 0.
    double horizonClip[4] = {0, 0, 0, 0};
    for (int r = 0; r < 4; ++r) horizonClip[r] = mvp(r, 0);
    const double horizonRow = (horizonClip[1] / horizonClip[3] * 0.5 + 0.5) * pixels;

    const auto rowOfEdge = [&](double distance) {
        double clip[4], px = 0, py = 0;
        mvp.transform({distance, 0, 0}, clip);
        if (!sim::clipToPixel(clip, pixels, pixels, px, py)) return 1e30;
        return py;
    };
    const double atReach = rowOfEdge(reach) - horizonRow;
    const double atAQuarter = rowOfEdge(0.25 * reach) - horizonRow;
    std::printf("     horizon: %.0f m reach for an eye %.0f m up at %d px puts its edge"
                " %.2f px below the horizon (%.2f px at a quarter of it)\n",
                reach, eyeHeight, pixels, atReach, atAQuarter);

    // Exactly one pixel, not "about one": with the eye level the offset works out
    // to pixelHeight * h / (2 D tan(fov / 2)), and the reach is defined as the
    // distance that makes that one. An identity, so it is asserted as one.
    expectNear("the sea's edge lands exactly one pixel below the horizon", atReach, 1.0, 1e-6);
    // Guard: the criterion is doing work. A quarter of the distance is four
    // pixels of sky between the sea and where the sea should be, and four pixels
    // of sky at the horizon is the defect this whole thing is about.
    expectNear("a quarter of that reach leaves four", atAQuarter, 4.0, 1e-5);

    // And the level count that gets there is the smallest one that does.
    const double inner = 200.0;
    const int levels = gpu::oceanCascadeLevelsFor(inner, reach);
    gpu::OceanCascade cascade;
    cascade.innerHalfExtent = inner;
    cascade.levels = levels;
    expectTrue(label("the cascade reaches at least as far as asked, m", cascade.reach()),
               cascade.reach() >= reach);
    gpu::OceanCascade shorter = cascade;
    shorter.levels = levels - 1;
    expectTrue(label("and one level fewer would not, m", shorter.reach()),
               shorter.reach() < reach);
    // Reach is exponential in the level count, which is the whole argument for
    // the shape: 63 times the distance for 6 more rings.
    expectEqual("each level doubles the reach",
                static_cast<long long>(std::llround(cascade.reach() / inner)),
                static_cast<long long>(1LL << (levels - 1)));
}

// The picture: with the cascade in place there is no background pixel anywhere
// below the horizon. Exact rather than approximate -- the elevation channel tags
// every sea fragment with blue 255 and the frame is cleared to blue 0, so "is
// this pixel sea" is an integer question.
void testCascadeLeavesNoSkyBelowTheHorizon() {
    gpu::Device device;
    gpu::OceanRenderer renderer;
    if (!setup(device, renderer)) return;

    const double t = 12.0;
    const double eyeHeight = 30.0;
    const double fov = 50.0 * sim::kDegToRad;
    const sim::WaveField field = spectralSea(4.0, 9.0, 16, 8, 0x5EA1EE);

    const double reach = gpu::oceanHorizonReach(eyeHeight, fov, static_cast<int>(kHeight));
    gpu::OceanCascade cascade;
    cascade.innerHalfExtent = 64.0;
    cascade.resolution = 64;
    cascade.levels = gpu::oceanCascadeLevelsFor(cascade.innerHalfExtent, reach);
    gpu::OceanCascadeSurface surface;
    surface.build(field, cascade, t);

    gpu::OceanView view;
    view.shading = gpu::OceanShading::Elevation;
    const float clear[4] = {0.0f, 0.0f, 0.0f, 1.0f};  // blue 0: not sea, unambiguously

    // **Four azimuths, not one.** A seam is a square, and one camera sees two of
    // its four sides. A mutation that dropped the transition cells along the two
    // sides behind the camera went straight through a single view -- it was the
    // edge census that caught it -- so the frame is asked about all four.
    const double pitch = 8.0 * sim::kDegToRad;
    const sim::Vec3 eye{0, 0, eyeHeight};
    const auto cameraFor = [&](int quadrant) {
        const double azimuth = quadrant * 0.5 * sim::kPi;
        const sim::Vec3 look{eye.x + std::cos(pitch) * std::cos(azimuth),
                             eye.y + std::cos(pitch) * std::sin(azimuth),
                             eye.z - std::sin(pitch)};
        return sim::perspective(fov, static_cast<double>(kWidth) / kHeight, 1.0,
                                3.0 * cascade.reach()) *
               sim::lookAt(eye, look, {0, 0, 1});
    };

    // Where the sea plane vanishes: the point at infinity along the view's own
    // horizontal direction, w = 0.
    const auto horizonRowOf = [&](const sim::Mat4& mvp, int quadrant) {
        const double azimuth = quadrant * 0.5 * sim::kPi;
        double clip[4] = {0, 0, 0, 0};
        for (int r = 0; r < 4; ++r)
            clip[r] = mvp(r, 0) * std::cos(azimuth) + mvp(r, 1) * std::sin(azimuth);
        return (clip[1] / clip[3] * 0.5 + 0.5) * kHeight;
    };

    const auto skyBelowHorizon = [&](const std::vector<gpu::OceanVertex>& vertices,
                                     const std::vector<std::uint32_t>& indices,
                                     const char* name, std::size_t& seaPixels,
                                     std::size_t& skyAbove) {
        long long sky = 0;
        seaPixels = 0;
        skyAbove = 0;
        for (int quadrant = 0; quadrant < 4; ++quadrant) {
            const sim::Mat4 mvp = cameraFor(quadrant);
            const double horizonRow = horizonRowOf(mvp, quadrant);
            float matrix[16];
            mvp.toFloats(matrix);
            Image rendered;
            if (!renderer.render(matrix, vertices, indices, view, clear, rendered)) return -1LL;
            if (quadrant == 0) core::writePng(outputDirectory() + "/" + name, rendered);
            for (std::uint32_t y = 0; y < rendered.height; ++y)
                for (std::uint32_t x = 0; x < rendered.width; ++x) {
                    double elevation = 0.0;
                    const bool sea =
                        gpu::decodeOceanElevation(view, rendered.pixel(x, y), elevation);
                    // Two rows of slack: the horizon is a line, the rasteriser
                    // fills half-open pixels, and the reach criterion is a
                    // one-pixel one.
                    if (static_cast<double>(y) > horizonRow + 2.0) {
                        if (sea) ++seaPixels; else ++sky;
                    } else if (!sea) {
                        ++skyAbove;
                    }
                }
        }
        return sky;
    };

    std::size_t seaPixels = 0, skyAbove = 0;
    const long long sky =
        skyBelowHorizon(surface.vertices(), surface.indices(), "ocean_cascade.png", seaPixels,
                        skyAbove);
    expectTrue("the cascade renders", sky >= 0);
    if (sky < 0) return;

    std::printf("     cascade: %d levels reaching %.0f m, horizon at row %.1f over 4 azimuths,"
                " %zu sea px below it, %lld sky px below it\n",
                cascade.levels, cascade.reach(), horizonRowOf(cameraFor(0), 0), seaPixels, sky);

    expectEqual("no pixel below the horizon is background", sky, 0);
    // Guards. The frame has to contain a horizon at all -- a camera pointed at
    // the water would pass the assertion above with nothing to say.
    expectTrue(label("most of the frame is sea, px", static_cast<double>(seaPixels)),
               static_cast<double>(seaPixels) > 4 * 0.5 * kWidth * kHeight);
    expectTrue(label("and there is sky above the horizon to bound it, px",
                     static_cast<double>(skyAbove)),
               static_cast<double>(skyAbove) > 4 * 0.1 * kWidth * kHeight);

    // Negative control 1: the uniform patch this replaces, at the same near-field
    // cell size and with *more* vertices than the whole cascade. It cannot reach,
    // and the frame says so.
    gpu::OceanGrid uniform{0.0, 0.0, cascade.innerHalfExtent * 8.0, 1};
    uniform.resolution =
        static_cast<int>(std::lround(2.0 * uniform.halfExtent / cascade.cellSize(0)));
    gpu::OceanSurface patch;
    patch.build(field, uniform, t);
    std::size_t patchSea = 0, patchSkyAbove = 0;
    const long long patchSky = skyBelowHorizon(patch.vertices(), patch.indices(),
                                               "ocean_cascade_uniform.png", patchSea,
                                               patchSkyAbove);
    std::printf("     the uniform patch it replaces: %d cells over %.0f m, %zu vertices,"
                " %lld sky px below the horizon\n",
                uniform.resolution, 2.0 * uniform.halfExtent, patch.vertices().size(), patchSky);
    expectTrue("the uniform patch at the same cell size has more vertices than the cascade",
               patch.vertices().size() > surface.vertices().size());
    expectTrue(label("and stops inside the frame, sky px", static_cast<double>(patchSky)),
               patchSky > 1000);

    // Negative control 2: the same cascade with its seams unstitched. Every
    // T-junction that opens downward is a hole straight through to the background.
    gpu::OceanCascade cracked = cascade;
    cracked.stitchSeams = false;
    gpu::OceanCascadeSurface unstitched;
    unstitched.build(field, cracked, t);
    std::size_t crackedSea = 0, crackedSkyAbove = 0;
    const long long crackedSky = skyBelowHorizon(unstitched.vertices(), unstitched.indices(),
                                                 "ocean_cascade_cracked.png", crackedSea,
                                                 crackedSkyAbove);
    std::printf("     the same cascade unstitched: %lld sky px below the horizon\n", crackedSky);
    expectTrue(label("an unstitched cascade shows background through its seams, px",
                     static_cast<double>(crackedSky)),
               crackedSky > 0);
}

// The near field must still be the sea the ship is in: the same comparison
// against `WaveField::elevation()` the uniform patch passes, over the cascade's
// level 0, to the same derived tolerance.
void testCascadeAgreesWithTheWaveFieldInItsNearField() {
    gpu::Device device;
    gpu::OceanRenderer renderer;
    if (!setup(device, renderer)) return;

    const double t = 21.5;
    const sim::WaveField field = spectralSea(3.0, 9.0, 32, 8, 0xC0FFEE);
    // The same 0.625 m cells over the same 140 m the uniform check uses, so the
    // tolerance and the result are comparable with the recorded ones.
    gpu::OceanCascade cascade;
    cascade.innerHalfExtent = 70.0;
    cascade.resolution = 224;
    cascade.levels = 5;
    gpu::OceanCascadeSurface surface;
    surface.build(field, cascade, t);

    gpu::OceanView view;
    view.shading = gpu::OceanShading::Elevation;
    view.elevationMin = -8.0f;
    view.elevationSpan = 16.0f;

    const sim::Mat4 mvp = overheadCamera(280.0, 30.0);
    const float clear[4] = {0.0f, 0.0f, 0.0f, 1.0f};

    float matrix[16];
    mvp.toFloats(matrix);
    Image rendered;
    expectTrue("the cascade's elevation channel renders",
               renderer.render(matrix, surface.vertices(), surface.indices(), view, clear,
                               rendered));
    if (!rendered.valid() || rendered.width != kWidth) return;
    core::writePng(outputDirectory() + "/ocean_cascade_elevation.png", rendered);

    const double h = cascade.cellSize(0);
    const double e = cascade.innerHalfExtent;
    const auto meshFromField = [&](double x, double y) {
        const double fx = std::clamp((x + e) / h, 0.0, cascade.resolution - 1e-9);
        const double fy = std::clamp((y + e) / h, 0.0, cascade.resolution - 1e-9);
        const int i = static_cast<int>(std::floor(fx));
        const int j = static_cast<int>(std::floor(fy));
        const double u = fx - i, v = fy - j;
        const double cornerX = -e + i * h, cornerY = -e + j * h;
        const double z00 = field.elevation(cornerX, cornerY, t);
        const double z11 = field.elevation(cornerX + h, cornerY + h, t);
        if (v <= u) {
            const double z10 = field.elevation(cornerX + h, cornerY, t);
            return z00 + u * (z10 - z00) + v * (z11 - z10);
        }
        const double z01 = field.elevation(cornerX, cornerY + h, t);
        return z00 + v * (z01 - z00) + u * (z11 - z01);
    };

    std::vector<sim::Vec3> probes;
    const int half = cascade.resolution / 2;
    for (int j = -half + 6; j < half; j += 11)
        for (int i = -half + 6; i < half; i += 11) {
            const gpu::OceanVertex& v = surface.vertices()[surface.vertexIndexAt(0, i, j)];
            probes.push_back({v.position[0], v.position[1], v.position[2]});
        }

    const ElevationAgreement found =
        compareRenderedElevation(rendered, mvp, view, field, t, probes, meshFromField);
    expectTrue("the cascade check had samples to work with", found.samples > 200);
    if (found.samples == 0) return;

    std::printf("     cascade near field vs WaveField::elevation: %lld samples, mean %+.5f m,"
                " rms %.5f m, worst %.5f m against a worst bound of %.5f m\n",
                found.samples, found.mean, found.rms, found.worstError, found.worstBound);

    expectEqual("every rendered cascade point agrees with WaveField::elevation inside its bound",
                found.exceeded, 0);
    expectNear("the cascade's near field carries no systematic offset", found.mean, 0.0, 0.002);
    // The recorded figures for the uniform patch under this camera are 0.6 mm rms
    // and 2.2 mm worst; the cascade's level 0 is the same grid and must not be
    // worse.
    expectTrue(label("the cascade holds the recorded near-field agreement, rms mm",
                     found.rms * 1000.0),
               found.rms < 0.001 && found.worstError < 0.004);
    expectTrue(label("the sampled surface varies by", found.highest - found.lowest),
               found.highest - found.lowest > 2.0);
    expectTrue(label("a one-second-stale field is rejected, fraction",
                     static_cast<double>(found.staleRejected) / static_cast<double>(found.samples)),
               found.staleRejected > found.samples * 8 / 10);
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

// The cascade against the uniform patch it replaces, at the same near-field cell
// size and the same reach. Measured on the shipped path, both of them.
void testCascadeCost() {
    // seaway_view's sea: Hs 4 m, Tp 9 s, 16 x 8.
    const sim::WaveField field = spectralSea(4.0, 9.0, 16, 8, 0x5eaf00d);
    const auto components = static_cast<double>(field.components().size());
    const double shortest = gpu::shortestWavelength(field);

    // The patch seaway_view drew: 3 Lpp across, eight cells across the shortest
    // component.
    const double reach = 525.0;
    gpu::OceanGrid uniform{0.0, 0.0, reach, gpu::oceanResolutionFor(reach, shortest, 8.0)};
    const double cell = uniform.cellSize();

    const auto timeIt = [](auto& surface, auto build) {
        double best = 1e30;
        for (int repeat = 0; repeat < 2; ++repeat) {
            build(surface, 1.0 + repeat);
            best = std::min(best, surface.buildSeconds());
        }
        return best;
    };

    gpu::OceanSurface patch;
    const double uniformSeconds = timeIt(patch, [&](gpu::OceanSurface& s, double t) {
        s.build(field, uniform, t);
    });

    // The same cell size at the centre and the same reach, as a cascade: level 0
    // is half the patch's extent, one ring takes it the rest of the way.
    gpu::OceanCascade cascade;
    cascade.innerHalfExtent = 0.5 * reach;
    cascade.resolution = static_cast<int>(std::lround(2.0 * cascade.innerHalfExtent / cell));
    cascade.levels = 2;
    gpu::OceanCascadeSurface matched;
    const double matchedSeconds = timeIt(matched, [&](gpu::OceanCascadeSurface& s, double t) {
        s.build(field, cascade, t);
    });

    // And the same cascade taken to the horizon for a 90 m eye at 720p, which is
    // where seaway_view's camera sits.
    gpu::OceanCascade far = cascade;
    far.levels = gpu::oceanCascadeLevelsFor(far.innerHalfExtent,
                                            gpu::oceanHorizonReach(90.0, 50.0 * sim::kDegToRad, 720));
    gpu::OceanCascadeSurface reaching;
    const double farSeconds = timeIt(reaching, [&](gpu::OceanCascadeSurface& s, double t) {
        s.build(field, far, t);
    });

    std::printf("     equal near-field cells (%.2f m) and equal reach (%.0f m), %g components:\n",
                cell, reach, components);
    std::printf("       uniform  %d cells   %7zu vertices  %7.2f ms\n", uniform.resolution,
                patch.vertices().size(), uniformSeconds * 1e3);
    std::printf("       cascade  %d levels  %7zu vertices  %7.2f ms\n", cascade.levels,
                matched.vertices().size(), matchedSeconds * 1e3);
    std::printf("       and to the horizon: %d levels, %.0f m, %7zu vertices  %7.2f ms\n",
                far.levels, far.reach(), reaching.vertices().size(), farSeconds * 1e3);

    // The cascade must actually be cheaper at equal reach, or the shape buys
    // nothing. It is cheaper twice over: three quarters of the outer area at a
    // quarter of the vertex density, and a ring that drops what it cannot carry.
    expectTrue(label("the cascade reaches as far as the patch, m", matched.reach()),
               matched.reach() >= reach - 1e-9);
    expectNear("at the same near-field cell size", matched.cascade().cellSize(0), cell,
               0.02 * cell);
    expectTrue(label("and costs less than the uniform patch, ratio",
                     matchedSeconds / uniformSeconds),
               matchedSeconds < uniformSeconds);
    // The claim the whole shape rests on: reach is exponential in the level count,
    // so the extra levels that take the sea from half a kilometre to the horizon
    // cost a bounded amount rather than the square of the distance. The uniform
    // patch reaching that far would need (reach / cell)^2 vertices.
    const double uniformAtHorizon = std::pow(2.0 * far.reach() / cell, 2.0);
    std::printf("       a uniform patch reaching %.0f m at %.2f m cells would be %.3g vertices,"
                " %.0f times the cascade\n",
                far.reach(), cell, uniformAtHorizon,
                uniformAtHorizon / static_cast<double>(reaching.vertices().size()));
    expectTrue("128 times the reach costs less than four times the cascade's time",
               farSeconds < 4.0 * matchedSeconds);
    expectTrue("and the uniform grid that reached that far would be four orders larger",
               uniformAtHorizon > 1e4 * static_cast<double>(reaching.vertices().size()));
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
    testCascadeIsOneWatertightSheet();
    testEveryResolutionRoundsUpAndStaysWatertight();
    testCascadeSeamVerticesAreShared();
    testRebuildingCarriesNothingStale();
    testCascadeElevationAgreesAcrossALevelBoundary();
    testCascadeDropsWhatItCannotResolve();
    testCascadeNormalsAreTheBandLimitedSlope();
    testCascadeReachesTheHorizon();
    testCascadeLeavesNoSkyBelowTheHorizon();
    testCascadeAgreesWithTheWaveFieldInItsNearField();
    testDisplacementCost();
    testCascadeCost();
}
