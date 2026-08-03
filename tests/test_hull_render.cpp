// SPDX-License-Identifier: MIT
//
// Validation of hull rendering and the material model.
//
// A lit hull looks like a ship whatever it is doing -- more so than the sea does,
// because a ship silhouette carries the recognition on its own and the shading can
// be arbitrarily wrong underneath it. So nothing here is eyeballed. The BRDF of
// `engine/gpu/material.hpp` is written out a second time in this file, from the
// formula rather than from the shader, and every shaded assertion is a pixel value
// predicted from it before the render:
//
//   * a hull face pointing exactly at the sun, and one exactly ninety degrees off
//     it, which is ambient and nothing else -- both exact;
//   * Lambert's cosine law over nineteen angles, with the specular lobe **cancelled
//     analytically** rather than modelled: the lobe is achromatic for a dielectric,
//     so the difference between two colour channels is pure diffuse, and
//     (red - blue) / (base_r - base_b) - sky is exactly (sun / pi) cos(theta);
//   * two materials on the same geometry under the same light, each predicted;
//   * the material a known hull point is painted with, found through
//     `sim::clipToPixel` and confirmed against the depth channel so an occluded
//     sample is skipped rather than silently compared against the wrong triangle.
//
// The instrument that is not a pixel comparison is the **energy balance**: the
// specular BRDF is integrated over the hemisphere and must not return more light
// than arrives. CLAUDE.md records that every subsystem here shipped green on its
// functional tests and was caught by a different instrument; a white-furnace
// integral is that instrument for a BRDF, and it is what would catch a wrong
// normalisation in D or a wrong Smith visibility while every pixel comparison
// stayed green against a shader making the same mistake.
//
// Skipped rather than failed when there is no GPU.
#include "engine/core/geometry.hpp"
#include "engine/core/math.hpp"
#include "engine/core/png.hpp"
#include "engine/gpu/hull.hpp"
#include "engine/gpu/material.hpp"
#include "engine/gpu/ocean.hpp"
#include "engine/sim/ship.hpp"
#include "engine/sim/waves.hpp"
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

bool announced = false;

std::string label(const char* what, double value) {
    char buffer[192];
    std::snprintf(buffer, sizeof buffer, "%s %.5g", what, value);
    return buffer;
}

bool setup(gpu::Device& device, gpu::HullRenderer& renderer) {
    std::string error;
    if (!device.create(error)) {
        if (!announced) {
            std::printf("     no usable GPU (%s) - hull render checks skipped\n", error.c_str());
            announced = true;
        }
        return false;
    }
    if (!renderer.create(device, kWidth, kHeight, SHIPSIM_SHADER_DIR, error)) {
        if (!announced) {
            std::printf("     hull renderer unavailable (%s) - render checks skipped\n",
                        error.c_str());
            announced = true;
        }
        return false;
    }
    return true;
}

// --- The reference BRDF -------------------------------------------------------
//
// Written from the formula in engine/gpu/material.hpp, not from shaders/hull.frag.
// The point of having two is that they are two: a test that asks the code under
// test what the answer should be cannot catch a wrong answer.

struct Reference {
    double base[3];
    double roughness;
    double metalness;
};

Reference referenceOf(const gpu::Material& material) {
    Reference out{};
    for (int c = 0; c < 3; ++c) out.base[c] = material.baseColour[c];
    out.roughness = material.roughness;
    out.metalness = material.metalness;
    return out;
}

// The specular lobe alone, per channel. Separated out because the energy integral
// needs it on its own.
void referenceSpecular(const Reference& m, const sim::Vec3& n, const sim::Vec3& l,
                       const sim::Vec3& v, double out[3]) {
    const sim::Vec3 h = normalize(l + v);
    const double ndl = std::max(dot(n, l), 0.0);
    const double ndv = std::max(dot(n, v), 1e-4);
    const double ndh = std::max(dot(n, h), 0.0);
    const double vdh = std::max(dot(v, h), 0.0);

    const double alpha = std::clamp(m.roughness, gpu::kMinRoughness, 1.0) *
                         std::clamp(m.roughness, gpu::kMinRoughness, 1.0);
    const double a2 = alpha * alpha;

    const double denominator = ndh * ndh * (a2 - 1.0) + 1.0;
    const double d = a2 / (sim::kPi * denominator * denominator);

    const double gv = ndl * std::sqrt(ndv * ndv * (1.0 - a2) + a2);
    const double gl = ndv * std::sqrt(ndl * ndl * (1.0 - a2) + a2);
    const double vis = 0.5 / std::max(gv + gl, 1e-9);

    for (int c = 0; c < 3; ++c) {
        const double f0 = 0.04 + (m.base[c] - 0.04) * m.metalness;
        const double fresnel = f0 + (1.0 - f0) * std::pow(1.0 - vdh, 5.0);
        out[c] = d * vis * fresnel;
    }
}

// Linear radiance before exposure and the clamp.
void referenceRadiance(const Reference& m, const sim::Vec3& n, const sim::Vec3& l,
                       const sim::Vec3& v, const gpu::SceneView& view, double out[3]) {
    double specular[3];
    referenceSpecular(m, n, l, v, specular);
    const double ndl = std::max(dot(n, l), 0.0);
    const double sky = 0.5 + 0.5 * n.z;
    for (int c = 0; c < 3; ++c) {
        const double diffuse = m.base[c] * (1.0 - m.metalness) / sim::kPi;
        out[c] = m.base[c] * static_cast<double>(view.skyColour[c]) * sky +
                 (diffuse + specular[c]) * ndl * static_cast<double>(view.sunColour[c]);
    }
}

// The UNORM8 code the renderer stores.
void referenceCode(const Reference& m, const sim::Vec3& n, const sim::Vec3& l, const sim::Vec3& v,
                   const gpu::SceneView& view, int out[3]) {
    double radiance[3];
    referenceRadiance(m, n, l, v, view, radiance);
    for (int c = 0; c < 3; ++c)
        out[c] = static_cast<int>(std::lround(
            std::clamp(radiance[c] * static_cast<double>(view.exposure), 0.0, 1.0) * 255.0));
}

sim::Vec3 sunOf(const gpu::SceneView& view) {
    return sim::normalize({view.sunDirection[0], view.sunDirection[1], view.sunDirection[2]});
}

// --- Scene helpers ------------------------------------------------------------

// A material set built the way a mod builds one: the shipped file first, then a
// second source that adds surfaces and overrides one. Nothing is compiled in.
bool testLibrary(gpu::MaterialLibrary& library, std::string& error) {
    if (!library.load(std::string(SHIPSIM_MATERIAL_DIR) + "/marine.materials", error)) return false;
    // Deliberately not round numbers, so a parser that dropped a digit shows up.
    static const char* kExtra =
        "# surfaces this suite needs, added the way a mod would add them\n"
        "material test_warm\n"
        "    base_colour 0.74 0.41 0.17\n"
        "    roughness   0.55\n"
        "    metalness   0.0\n"
        "material test_cool\n"
        "    base_colour 0.19 0.36 0.66\n"
        "    roughness   0.28\n"
        "    metalness   0.0\n"
        "material test_grey\n"
        "    base_colour 0.62 0.62 0.62\n"
        "    roughness   0.47\n"
        "    metalness   0.0\n";
    return library.parse(kExtra, "<test surfaces>", error);
}

// A 100 m hull: port/starboard symmetric by construction, and deliberately **not**
// fore-and-aft symmetric. The asymmetry is load bearing -- the mirror check views
// the ship abeam, so a fore-aft symmetric hull would make the image its own mirror
// and the comparison could not fail. The check guards against exactly that, and
// this is the hull that makes the guard pass honestly.
sim::Ship makeTestShip() {
    const std::vector<double> waterlines{0.0, 1.2, 2.4, 3.6, 5.0, 6.5, 8.5, 11.0};
    std::vector<sim::Station> stations;
    for (int i = 0; i <= 24; ++i) {
        sim::Station station;
        station.x = -50.0 + 100.0 * i / 24.0;
        const double u = station.x / 50.0;
        // A fine bow and a full stern, as a ship has.
        const double longitudinal =
            u > 0.35 ? 1.0 - 0.90 * std::pow((u - 0.35) / 0.65, 2.0)
                     : (u < -0.35 ? 1.0 - 0.45 * std::pow((-u - 0.35) / 0.65, 2.0) : 1.0);
        for (double z : waterlines) {
            const double vertical = z >= 3.6 ? 1.0 : 0.32 + 0.68 * std::pow(z / 3.6, 0.75);
            station.halfBeam.push_back(8.5 * longitudinal * vertical);
        }
        stations.push_back(station);
    }

    sim::Ship ship;
    ship.hull = sim::makeHullFromStations(stations, waterlines);
    ship.deckEdgeZ = 8.5;
    // Light displacement defined as whatever floats her at 5 m, so the hull form
    // and the loading stay consistent if the offsets above are edited.
    const double designVolume =
        sim::integrateBelowPlane(ship.hull, sim::Vec3{0, 0, 1}, 5.0).volume;
    ship.lightshipMass = designVolume * sim::kRhoSeawater;
    ship.lightshipCog = {-1.0, 0.0, 6.2};
    ship.gyradii = {6.0, 24.0, 25.0};
    return ship;
}

// A prism whose side quads are split into four triangles about their own centre.
// That construction is what makes the smooth normal at a mid-height ring vertex an
// exact closed form: every side quad is congruent and every corner takes two
// quarter-area triangles from each quad it touches, so the area-weighted average is
// exactly the bisector of the two adjacent facet normals -- which for a regular
// polygon is exactly the vertex's own radial direction.
sim::TriMesh makePrism(int columns, int rows, double radius, double height) {
    sim::TriMesh mesh;
    const auto ring = [&](int i, int j) {
        return static_cast<std::uint32_t>(j * columns + (i % columns));
    };
    for (int j = 0; j <= rows; ++j)
        for (int i = 0; i < columns; ++i) {
            const double angle = 2.0 * sim::kPi * i / columns;
            mesh.verts.push_back({radius * std::cos(angle), radius * std::sin(angle),
                                  height * j / rows});
        }
    for (int j = 0; j < rows; ++j)
        for (int i = 0; i < columns; ++i) {
            const std::uint32_t v00 = ring(i, j), v10 = ring(i + 1, j);
            const std::uint32_t v01 = ring(i, j + 1), v11 = ring(i + 1, j + 1);
            const auto centre = static_cast<std::uint32_t>(mesh.verts.size());
            mesh.verts.push_back((mesh.verts[v00] + mesh.verts[v10] + mesh.verts[v11] +
                                  mesh.verts[v01]) *
                                 0.25);
            mesh.tris.push_back({v00, v10, centre});
            mesh.tris.push_back({v10, v11, centre});
            mesh.tris.push_back({v11, v01, centre});
            mesh.tris.push_back({v01, v00, centre});
        }
    const auto top = static_cast<std::uint32_t>(mesh.verts.size());
    mesh.verts.push_back({0, 0, height});
    const auto bottom = static_cast<std::uint32_t>(mesh.verts.size());
    mesh.verts.push_back({0, 0, 0});
    for (int i = 0; i < columns; ++i) {
        mesh.tris.push_back({top, ring(i, rows), ring(i + 1, rows)});
        mesh.tris.push_back({bottom, ring(i + 1, 0), ring(i, 0)});
    }
    return mesh;
}

// Two strips meeting along a ridge, one of them `ratio` times the area of the
// other, with the same number of triangles on each. Area-weighted averaging and
// plain averaging then give measurably different normals at the ridge, and only
// one of them is right -- a fan of slivers must not outvote the large face it
// sits against. Quads are split about their own centres for the same reason as
// the prism: every corner then takes exactly half of its quad's area.
sim::TriMesh makeRidge(double wide, double narrow, double tiltRadians) {
    sim::TriMesh mesh;
    const sim::Vec3 ridge0{0, 0, 0}, ridge1{6, 0, 0};
    const sim::Vec3 farA{0, -wide, 0};
    const sim::Vec3 farB{0, narrow * std::cos(tiltRadians), narrow * std::sin(tiltRadians)};

    const auto strip = [&](const sim::Vec3& a, const sim::Vec3& b, const sim::Vec3& offset) {
        const auto base = static_cast<std::uint32_t>(mesh.verts.size());
        mesh.verts.push_back(a);
        mesh.verts.push_back(b);
        mesh.verts.push_back(b + offset);
        mesh.verts.push_back(a + offset);
        const auto centre = static_cast<std::uint32_t>(mesh.verts.size());
        mesh.verts.push_back((mesh.verts[base] + mesh.verts[base + 1] + mesh.verts[base + 2] +
                              mesh.verts[base + 3]) *
                             0.25);
        for (std::uint32_t k = 0; k < 4; ++k)
            mesh.tris.push_back({base + k, base + (k + 1) % 4, centre});
    };
    strip(ridge0, ridge1, farA - ridge0);
    strip(ridge1, ridge0, farB - ridge0);
    return mesh;
}

sim::Mat4 perspectiveFor(double distance) {
    return sim::perspective(45.0 * sim::kDegToRad, static_cast<double>(kWidth) / kHeight, 1.0,
                            8.0 * distance);
}

// Camera on a bearing, aimed at `target`. Returns the mvp and fills the eye the
// SceneView needs -- one function, so the two can never disagree.
sim::Mat4 camera(const sim::Vec3& eye, const sim::Vec3& target, double distance,
                 gpu::SceneView& view) {
    view.eye[0] = static_cast<float>(eye.x);
    view.eye[1] = static_cast<float>(eye.y);
    view.eye[2] = static_cast<float>(eye.z);
    // World up, except looking straight down it, where `lookAt`'s cross product
    // degenerates and every subsequent frame is empty. That cost this file two
    // whole tests on its first run -- the same shape of mistake as the bow-on
    // camera CLAUDE.md records, and just as silent: the renderer was fine and the
    // camera was nonsense.
    const sim::Vec3 forward = sim::normalize(target - eye);
    const sim::Vec3 up =
        std::abs(forward.z) > 0.999 ? sim::Vec3{1, 0, 0} : sim::Vec3{0, 0, 1};
    return perspectiveFor(distance) * sim::lookAt(eye, target, up);
}

bool pixelOf(const sim::Mat4& mvp, const sim::Vec3& world, std::uint32_t& column,
             std::uint32_t& row) {
    double clip[4], x = 0, y = 0;
    mvp.transform(world, clip);
    if (!sim::clipToPixel(clip, kWidth, kHeight, x, y)) return false;
    if (x < 0.5 || y < 0.5 || x > kWidth - 1.5 || y > kHeight - 1.5) return false;
    column = static_cast<std::uint32_t>(x);
    row = static_cast<std::uint32_t>(y);
    return true;
}

// Sub-pixel distance from a triangle's projected centroid to its nearest projected
// edge, or -1 if any corner falls behind the eye. A pixel is an area sample, so a
// comparison against the triangle one predicted is only meaningful where that pixel
// is unambiguously covered by it -- on a sliver, the centroid's pixel belongs as
// much to the neighbour, and at a paint-band boundary the neighbour is a different
// colour at almost exactly the same depth.
double centroidClearanceInPixels(const sim::Mat4& mvp, const sim::Vec3& a, const sim::Vec3& b,
                                 const sim::Vec3& c) {
    double corner[3][2];
    const sim::Vec3 world[3] = {a, b, c};
    for (int i = 0; i < 3; ++i) {
        double clip[4];
        mvp.transform(world[i], clip);
        if (!sim::clipToPixel(clip, kWidth, kHeight, corner[i][0], corner[i][1])) return -1.0;
    }
    double centre[2];
    double clip[4];
    mvp.transform((a + b + c) / 3.0, clip);
    if (!sim::clipToPixel(clip, kWidth, kHeight, centre[0], centre[1])) return -1.0;

    double clearance = 1e30;
    for (int i = 0; i < 3; ++i) {
        const int j = (i + 1) % 3;
        const double ex = corner[j][0] - corner[i][0], ey = corner[j][1] - corner[i][1];
        const double edge = std::hypot(ex, ey);
        if (edge < 1e-9) return 0.0;
        const double cross = std::abs(ex * (centre[1] - corner[i][1]) -
                                      ey * (centre[0] - corner[i][0]));
        clearance = std::min(clearance, cross / edge);
    }
    return clearance;
}

double clipDepthOf(const sim::Mat4& mvp, const sim::Vec3& world) {
    double clip[4];
    mvp.transform(world, clip);
    return clip[3] > 1e-12 ? clip[2] / clip[3] : -1.0;
}

bool renderThroughPng(gpu::HullRenderer& renderer, const sim::Mat4& mvp,
                      const gpu::SceneView& view, const gpu::SceneMesh& mesh,
                      const gpu::MaterialLibrary& library, const float clear[4],
                      const std::string& name, Image& out) {
    float matrix[16];
    mvp.toFloats(matrix);
    Image rendered;
    if (!renderer.render(matrix, view, mesh, library, clear, rendered)) return false;
    const std::string path = testing::scratchDir() + name;
    if (!core::writePng(path, rendered)) return false;
    return core::readPng(path, out);
}

bool sameColour(const std::uint8_t* a, const std::uint8_t* b, int tolerance) {
    return std::abs(int{a[0]} - int{b[0]}) <= tolerance &&
           std::abs(int{a[1]} - int{b[1]}) <= tolerance &&
           std::abs(int{a[2]} - int{b[2]}) <= tolerance;
}

// A single quad, four vertices and two triangles, in the plane its normal implies.
sim::TriMesh makeQuad(const sim::Vec3& centre, const sim::Vec3& across, const sim::Vec3& up) {
    sim::TriMesh mesh;
    mesh.verts = {centre - across - up, centre + across - up, centre + across + up,
                  centre - across + up};
    mesh.tris = {{0, 1, 2}, {0, 2, 3}};
    return mesh;
}

// --- Materials are data -------------------------------------------------------

void testMaterialSetIsData() {
    gpu::MaterialLibrary library;
    std::string error;
    expectTrue("the shipped marine material set loads: " + error,
               library.load(std::string(SHIPSIM_MATERIAL_DIR) + "/marine.materials", error));
    if (library.empty()) return;

    const char* wanted[] = {"painted_steel_topside", "boot_topping", "antifouling",
                            "bare_steel",            "rusted_steel", "timber_deck",
                            "glass",                 "sea_water"};
    for (const char* name : wanted)
        expectTrue(std::string("the shipped set has '") + name + "'", library.find(name) >= 0);

    // Properties of the *set*, which is what the values are for. Paint over steel
    // is a dielectric however metallic the plate underneath is, so exactly one
    // surface here may be a conductor; rust is an oxide and must not be one.
    int conductors = 0;
    double smallestRoughness = 2.0;
    std::string smoothest;
    for (const gpu::Material& material : library.materials()) {
        if (material.metalness > 0.5) ++conductors;
        if (material.roughness < smallestRoughness) {
            smallestRoughness = material.roughness;
            smoothest = material.name;
        }
    }
    expectEqual("exactly one shipped surface is a conductor", conductors, 1);
    expectTrue("and it is the bare steel",
               library[static_cast<std::size_t>(library.find("bare_steel"))].metalness == 1.0);
    expectTrue("rust is an oxide, not a conductor",
               library[static_cast<std::size_t>(library.find("rusted_steel"))].metalness == 0.0);
    expectTrue("glass is the smoothest thing on the ship, not " + smoothest, smoothest == "glass");
    expectTrue("antifouling is matt, as an ablating coating has to be",
               library[static_cast<std::size_t>(library.find("antifouling"))].roughness > 0.6);

    // The packed table is what the shader indexes, so a swap between roughness and
    // metalness on the way there would be invisible in `materials()` and wrong in
    // every pixel.
    bool packedAgrees = library.packed().size() == library.size();
    for (std::size_t i = 0; packedAgrees && i < library.size(); ++i) {
        const gpu::Material& source = library[i];
        const gpu::GpuMaterial& packed = library.packed()[i];
        for (int c = 0; c < 3; ++c)
            packedAgrees = packedAgrees &&
                           packed.baseColour[c] == static_cast<float>(source.baseColour[c]);
        packedAgrees = packedAgrees && packed.baseColour[3] == static_cast<float>(source.opacity) &&
                       packed.params[0] == static_cast<float>(source.roughness) &&
                       packed.params[1] == static_cast<float>(source.metalness);
    }
    expectTrue("the packed table agrees with the authored materials field for field",
               packedAgrees);

    // What a mod does: load a second file, add a surface, restate one. An index
    // resolved before the mod loaded must still point at the same surface.
    const int deckBefore = library.find("timber_deck");
    const std::size_t sizeBefore = library.size();
    const std::uint64_t revisionBefore = library.revision();
    const char* mod =
        "material scorched_deck\n"
        "    base_colour 0.06 0.05 0.045\n"
        "    roughness   0.95\n"
        "    metalness   0.0\n"
        "material timber_deck        # a mod restating a shipped surface\n"
        "    base_colour 0.44 0.33 0.19\n"
        "    roughness   0.71\n"
        "    metalness   0.0\n";
    expectTrue("a mod loads on top: " + error, library.parse(mod, "<mod>", error));
    expectEqual("the mod added exactly one surface", static_cast<long long>(library.size()),
                static_cast<long long>(sizeBefore) + 1);
    expectEqual("an index resolved before the mod still points at the same surface",
                library.find("timber_deck"), deckBefore);
    expectNear("and that surface now carries the mod's roughness",
               library[static_cast<std::size_t>(deckBefore)].roughness, 0.71, 0.0);
    expectTrue("the revision moved, so a renderer knows to re-upload",
               library.revision() != revisionBefore);
    expectTrue("the new surface is findable by name", library.find("scorched_deck") >= 0);
}

// A loader that fails open leaves a material with a name and default everything
// else, which renders, and looks almost right. CLAUDE.md records that exact
// failure in `World::load`; this is the same instrument pointed at this loader.
void testMaterialLoaderFailsClosed() {
    const std::string valid =
        "material alpha\n"
        "    base_colour 0.11 0.22 0.33\n"
        "    roughness   0.44\n"
        "    metalness   0.55\n"
        "material beta\n"
        "    base_colour 0.66 0.77 0.88\n"
        "    roughness   0.99\n"
        "    metalness   0.10\n"
        "    opacity     0.25\n";

    {
        gpu::MaterialLibrary library;
        std::string error;
        expectTrue("the reference file parses: " + error, library.parse(valid, "<ref>", error));
        expectEqual("it holds two materials", static_cast<long long>(library.size()), 2);
        expectNear("opacity defaults to one when it is not stated",
                   library[0].opacity, 1.0, 0.0);
        expectNear("and is read when it is", library[1].opacity, 0.25, 0.0);
        expectNear("metalness is not roughness", library[0].metalness, 0.55, 0.0);
        expectNear("roughness is not metalness", library[0].roughness, 0.44, 0.0);
    }

    struct Bad { const char* what; const char* text; };
    const Bad rejected[] = {
        {"an unknown key", "material a\n base_colour 1 1 1\n roughness 0.5\n metalness 0\n"
                           " shininess 4\n"},
        {"a key outside any block", "roughness 0.5\n"},
        {"a block with no roughness", "material a\n base_colour 1 1 1\n metalness 0\n"},
        {"a block with no metalness", "material a\n base_colour 1 1 1\n roughness 0.5\n"},
        {"a block with no base colour", "material a\n roughness 0.5\n metalness 0\n"},
        {"a colour with two components", "material a\n base_colour 1 1\n roughness 0.5\n"
                                         " metalness 0\n"},
        {"a colour above one", "material a\n base_colour 1.4 1 1\n roughness 0.5\n metalness 0\n"},
        {"a negative colour", "material a\n base_colour -0.1 1 1\n roughness 0.5\n metalness 0\n"},
        {"a roughness of zero", "material a\n base_colour 1 1 1\n roughness 0\n metalness 0\n"},
        {"a roughness above one", "material a\n base_colour 1 1 1\n roughness 1.5\n metalness 0\n"},
        {"a metalness above one", "material a\n base_colour 1 1 1\n roughness 0.5\n metalness 2\n"},
        {"a number with a trailing letter",
         "material a\n base_colour 1 1 1\n roughness 0.5x\n metalness 0\n"},
        {"a nameless material", "material\n base_colour 1 1 1\n roughness 0.5\n metalness 0\n"},
        {"a material named twice in one file",
         "material a\n base_colour 1 1 1\n roughness 0.5\n metalness 0\n"
         "material a\n base_colour 0 0 0\n roughness 0.5\n metalness 0\n"},
    };

    for (const Bad& bad : rejected) {
        // Loaded on top of a good library, because the failure that matters is not
        // "returns false" but "leaves something behind".
        gpu::MaterialLibrary library;
        std::string error;
        library.parse(valid, "<ref>", error);
        const std::uint64_t revision = library.revision();
        expectTrue(std::string("rejected: ") + bad.what,
                   !library.parse(bad.text, "<bad>", error));
        expectEqual(std::string("and left the library untouched after ") + bad.what,
                    static_cast<long long>(library.size()), 2);
        expectTrue(std::string("and did not bump the revision after ") + bad.what,
                   library.revision() == revision);
        expectTrue(std::string("and said why for ") + bad.what,
                   error.find("<bad>:") != std::string::npos);
    }

    // Every truncation of a valid file, at a line boundary. Whether a prefix is
    // legal is predictable without the loader: it is legal exactly when it ends
    // after a complete block, so the expectation is derived rather than observed.
    std::vector<std::size_t> lineEnds;
    for (std::size_t i = 0; i < valid.size(); ++i)
        if (valid[i] == '\n') lineEnds.push_back(i + 1);

    long long acceptedPrefixes = 0, rejectedPrefixes = 0;
    bool acceptanceAsPredicted = true, contentsComplete = true;
    for (std::size_t cut : lineEnds) {
        const std::string prefix = valid.substr(0, cut);
        // Count the required keys of the last block declared in this prefix.
        int blocks = 0, keys = 0;
        std::size_t position = 0;
        while (position < prefix.size()) {
            const std::size_t end = prefix.find('\n', position);
            const std::string line =
                prefix.substr(position, end == std::string::npos ? std::string::npos
                                                                 : end - position);
            position = end == std::string::npos ? prefix.size() : end + 1;
            if (line.find("material ") != std::string::npos) { ++blocks; keys = 0; }
            else if (line.find("base_colour") != std::string::npos ||
                     line.find("roughness") != std::string::npos ||
                     line.find("metalness") != std::string::npos)
                ++keys;
        }
        const bool shouldLoad = blocks == 0 || keys == 3;

        gpu::MaterialLibrary library;
        std::string error;
        const bool loaded = library.parse(prefix, "<truncated>", error);
        acceptanceAsPredicted = acceptanceAsPredicted && loaded == shouldLoad;
        if (loaded) {
            ++acceptedPrefixes;
            // Nothing half built: whatever survived must be the complete block it
            // came from, values and all.
            contentsComplete = contentsComplete && library.size() == static_cast<std::size_t>(blocks);
            if (library.size() >= 1)
                contentsComplete = contentsComplete && library[0].name == "alpha" &&
                                   library[0].roughness == 0.44 && library[0].metalness == 0.55;
            // Not opacity: it is optional, so `beta` is already complete one line
            // before its opacity is stated and a prefix cut there legitimately
            // carries the default.
            if (library.size() >= 2)
                contentsComplete = contentsComplete && library[1].name == "beta" &&
                                   library[1].roughness == 0.99 && library[1].metalness == 0.10;
        } else {
            ++rejectedPrefixes;
            expectEqual("a rejected prefix leaves nothing behind",
                        static_cast<long long>(library.size()), 0);
        }
    }
    expectTrue("every line truncation loads exactly when it ends after a complete block",
               acceptanceAsPredicted);
    expectTrue("no truncation produced a half-built material", contentsComplete);
    // Guards: a loader that accepted everything, or nothing, would satisfy one of
    // the two above on its own.
    expectTrue(label("some truncations were accepted", static_cast<double>(acceptedPrefixes)),
               acceptedPrefixes >= 3);
    expectTrue(label("and some were rejected", static_cast<double>(rejectedPrefixes)),
               rejectedPrefixes >= 4);

    // And every byte truncation, which can cut a number in half. The invariant is
    // weaker -- "0.4" is a legal truncation of "0.44" -- but the important half is
    // not: a failure must leave nothing behind, and a success must never invent a
    // material that was not declared.
    long long byteAccepted = 0;
    bool bytesSane = true;
    for (std::size_t cut = 0; cut <= valid.size(); ++cut) {
        gpu::MaterialLibrary library;
        std::string error;
        if (library.parse(valid.substr(0, cut), "<bytes>", error)) {
            ++byteAccepted;
            bytesSane = bytesSane && library.size() <= 2;
            for (std::size_t i = 0; i < library.size(); ++i)
                bytesSane = bytesSane &&
                            library[i].name == (i == 0 ? "alpha" : "beta") &&
                            library[i].roughness >= gpu::kMinRoughness;
        } else {
            bytesSane = bytesSane && library.empty();
        }
    }
    expectTrue("no byte truncation invents a material or survives a failure", bytesSane);
    expectTrue(label("byte truncations accepted", static_cast<double>(byteAccepted)),
               byteAccepted > 0 && byteAccepted < static_cast<long long>(valid.size()));
}

// --- The instrument that is not a pixel ---------------------------------------

// A BRDF that reflects more light than arrives is unphysical, and no pixel
// comparison can see it when the shader and the prediction make the same mistake.
// Integrating the specular lobe over the hemisphere is the check that can: with
// F0 driven to one -- a white furnace -- the directional albedo is a pure number
// bounded above by one, and a wrong normalisation in D or a wrong Smith visibility
// moves it immediately.
void testSpecularBrdfConservesEnergy() {
    // Sampled uniformly in cos(theta) so the cosine factor is absorbed into the
    // measure and the quadrature does not have to resolve a peak at grazing.
    constexpr int kCosineSamples = 512;
    constexpr int kAzimuthSamples = 512;

    const auto directionalAlbedo = [](double roughness, double ndv) {
        Reference furnace{{1.0, 1.0, 1.0}, roughness, 1.0};  // metal, base 1 -> F0 = 1
        const sim::Vec3 n{0, 0, 1};
        const sim::Vec3 v{std::sqrt(std::max(0.0, 1.0 - ndv * ndv)), 0.0, ndv};
        double total = 0.0;
        for (int i = 0; i < kCosineSamples; ++i) {
            const double cosTheta = (i + 0.5) / kCosineSamples;
            const double sinTheta = std::sqrt(std::max(0.0, 1.0 - cosTheta * cosTheta));
            for (int k = 0; k < kAzimuthSamples; ++k) {
                const double phi = 2.0 * sim::kPi * (k + 0.5) / kAzimuthSamples;
                const sim::Vec3 l{sinTheta * std::cos(phi), sinTheta * std::sin(phi), cosTheta};
                double specular[3];
                referenceSpecular(furnace, n, l, v, specular);
                // integral f cos(theta) dOmega, with dOmega = dcos dphi.
                total += specular[0] * cosTheta;
            }
        }
        return total * 2.0 * sim::kPi / (kCosineSamples * kAzimuthSamples);
    };

    const double roughnesses[] = {0.20, 0.35, 0.50, 0.70, 1.00};
    const double viewCosines[] = {0.25, 0.60, 0.95};
    double worst = 0.0, smallest = 1e30;
    std::printf("     white furnace, directional albedo of the specular lobe:\n");
    for (double roughness : roughnesses) {
        std::printf("       roughness %.2f: ", roughness);
        for (double ndv : viewCosines) {
            const double albedo = directionalAlbedo(roughness, ndv);
            std::printf(" n.v %.2f -> %.4f ", ndv, albedo);
            worst = std::max(worst, albedo);
            smallest = std::min(smallest, albedo);
            expectTrue(label("a white furnace reflects no more than it receives, at", albedo),
                       albedo <= 1.0 + 5e-3);
            // A lobe that had lost its normalisation entirely would come back near
            // zero and still satisfy the bound above, so the floor matters as much.
            // It is 0.25 and not more: single-scattering GGX genuinely loses about
            // two thirds of the energy at alpha = 1, which is the multi-scatter
            // deficit the literature tabulates. The first version of this check
            // demanded 0.45 and failed at 0.32 -- the expectation was wrong, not
            // the model, and the measurement below is the record of that.
            expectTrue(label("and it does reflect, at", albedo), albedo > 0.25);
        }
        std::printf("\n");
    }
    std::printf("     worst %.4f, smallest %.4f -- the shortfall is the single-scattering\n"
                "     energy loss GGX is known for, and is the multi-scatter term's job\n",
                worst, smallest);
    // The single-scattering deficit grows with roughness. That it does is a
    // property of the model rather than an accident, and it is worth asserting
    // because a visibility term with the wrong argument order would not have it.
    expectTrue("a smooth white furnace loses almost nothing",
               directionalAlbedo(0.20, 0.60) > 0.93);
    expectTrue("a rough one loses a great deal more",
               directionalAlbedo(1.00, 0.60) < directionalAlbedo(0.20, 0.60) - 0.15);

    // Helmholtz reciprocity: swapping the light and the eye must not change the
    // BRDF. An asymmetric visibility term is a real and quiet mistake.
    double worstReciprocity = 0.0;
    Reference surface{{0.6, 0.4, 0.2}, 0.4, 0.3};
    for (int i = 0; i < 24; ++i) {
        const double a = 0.13 * i, b = 0.37 * i + 0.4;
        const sim::Vec3 n{0, 0, 1};
        const sim::Vec3 l = sim::normalize({std::cos(a), std::sin(a), 0.35 + 0.02 * i});
        const sim::Vec3 v = sim::normalize({std::cos(b), std::sin(b), 0.25 + 0.03 * i});
        double forward[3], backward[3];
        referenceSpecular(surface, n, l, v, forward);
        referenceSpecular(surface, n, v, l, backward);
        worstReciprocity = std::max(worstReciprocity, std::abs(forward[0] - backward[0]));
    }
    expectTrue(label("the BRDF is reciprocal, worst", worstReciprocity), worstReciprocity < 1e-12);
}

// --- Normals come from the mesh -----------------------------------------------

void testFlatNormalsAreTheGeometry() {
    gpu::MaterialLibrary library;
    std::string error;
    if (!testLibrary(library, error)) {
        expectTrue("the material library loads: " + error, false);
        return;
    }
    const auto grey = static_cast<std::uint32_t>(library.find("test_grey"));

    const sim::Vec3 lo{-3, -5, -7}, hi{4, 6, 9};
    const sim::TriMesh box = sim::makeBox(lo, hi);
    expectTrue("the box is a closed, consistently wound manifold", sim::isClosedManifold(box));

    gpu::SceneMesh mesh;
    gpu::HullShading shading;
    shading.normals = gpu::HullNormals::Flat;
    mesh.appendMesh(box, sim::Mat3::identity(), {0, 0, 0}, grey, shading);

    expectEqual("flat shading emits three vertices per triangle",
                static_cast<long long>(mesh.vertices().size()),
                static_cast<long long>(box.tris.size()) * 3);
    expectEqual("and one index per corner", static_cast<long long>(mesh.indices().size()),
                static_cast<long long>(box.tris.size()) * 3);

    const sim::Vec3 centre = (lo + hi) * 0.5;
    int perAxis[6] = {0, 0, 0, 0, 0, 0};
    bool axisAligned = true, outward = true, materialKept = true;
    for (const gpu::HullVertex& vertex : mesh.vertices()) {
        const sim::Vec3 n{vertex.normal[0], vertex.normal[1], vertex.normal[2]};
        const sim::Vec3 p{vertex.position[0], vertex.position[1], vertex.position[2]};
        materialKept = materialKept && vertex.material == grey;
        int matched = -1;
        for (int axis = 0; axis < 3; ++axis)
            for (int sign = 0; sign < 2; ++sign) {
                sim::Vec3 want{};
                want[axis] = sign == 0 ? 1.0 : -1.0;
                if (std::abs(dot(n, want) - 1.0) < 1e-6) matched = axis * 2 + sign;
            }
        if (matched < 0) { axisAligned = false; continue; }
        ++perAxis[matched];
        outward = outward && dot(n, p - centre) > 0.0;
    }
    expectTrue("every flat normal on a box is exactly an axis direction", axisAligned);
    expectTrue("every flat normal points out of the box, so the winding is outward", outward);
    expectTrue("the material index survives the build", materialKept);
    bool sixPerFace = true;
    for (int face = 0; face < 6; ++face) sixPerFace = sixPerFace && perAxis[face] == 6;
    expectTrue("each of the six faces contributes two triangles' worth of corners", sixPerFace);
}

void testSmoothNormalsRespectCreases() {
    gpu::MaterialLibrary library;
    std::string error;
    if (!testLibrary(library, error)) return;
    const auto grey = static_cast<std::uint32_t>(library.find("test_grey"));

    const sim::TriMesh box = sim::makeBox({-2, -2, -2}, {2, 2, 2});
    gpu::SceneMesh flat, creased, blended;
    gpu::HullShading flatShading;
    flatShading.normals = gpu::HullNormals::Flat;
    gpu::HullShading tightCrease;
    tightCrease.normals = gpu::HullNormals::Smooth;
    // Thirty degrees. At a box corner three perpendicular faces meet, and the
    // area-weighted average is at least 35.3 degrees from each of them however the
    // triangulation splits the faces -- so every edge of a box creases, and the
    // smooth result must be bit-identical to the flat one.
    tightCrease.creaseAngle = 30.0 * sim::kDegToRad;
    gpu::HullShading looseCrease = tightCrease;
    looseCrease.creaseAngle = 70.0 * sim::kDegToRad;

    flat.appendMesh(box, sim::Mat3::identity(), {0, 0, 0}, grey, flatShading);
    creased.appendMesh(box, sim::Mat3::identity(), {0, 0, 0}, grey, tightCrease);
    blended.appendMesh(box, sim::Mat3::identity(), {0, 0, 0}, grey, looseCrease);

    bool identical = flat.vertices().size() == creased.vertices().size();
    for (std::size_t i = 0; identical && i < flat.vertices().size(); ++i)
        identical = std::memcmp(&flat.vertices()[i], &creased.vertices()[i],
                                sizeof(gpu::HullVertex)) == 0;
    expectTrue("a thirty-degree crease leaves every edge of a box hard", identical);

    // Guard: identical results prove nothing if the crease angle does nothing.
    std::size_t rounded = 0;
    for (std::size_t i = 0; i < flat.vertices().size() && i < blended.vertices().size(); ++i)
        if (std::memcmp(&flat.vertices()[i], &blended.vertices()[i], sizeof(gpu::HullVertex)) != 0)
            ++rounded;
    expectTrue(label("a seventy-degree crease rounds the box's corners instead, vertices",
                     static_cast<double>(rounded)),
               rounded > flat.vertices().size() / 2);

    // The closed form. On this prism the smooth normal at a mid-height ring vertex
    // is exactly the bisector of the two facet normals meeting there, which for a
    // regular polygon is exactly that vertex's own radial direction.
    const int columns = 32;
    const double radius = 5.0, height = 10.0;
    const sim::TriMesh prism = makePrism(columns, 2, radius, height);
    expectTrue("the prism is a closed, consistently wound manifold",
               sim::isClosedManifold(prism));

    gpu::SceneMesh smoothPrism, flatPrism;
    gpu::HullShading prismShading;
    prismShading.normals = gpu::HullNormals::Smooth;
    smoothPrism.appendMesh(prism, sim::Mat3::identity(), {0, 0, 0}, grey, prismShading);
    flatPrism.appendMesh(prism, sim::Mat3::identity(), {0, 0, 0}, grey, flatShading);

    long long ringSamples = 0;
    double worstRadial = 0.0, worstVertical = 0.0, worstFlatGap = 0.0;
    for (std::size_t i = 0; i < smoothPrism.vertices().size(); ++i) {
        const gpu::HullVertex& vertex = smoothPrism.vertices()[i];
        const double x = vertex.position[0], y = vertex.position[1], z = vertex.position[2];
        // Mid-height ring only: a rim vertex also touches the end cap, which is a
        // genuine crease and is a different question.
        if (std::abs(z - 0.5 * height) > 1e-4) continue;
        if (std::abs(std::hypot(x, y) - radius) > 1e-3) continue;
        ++ringSamples;
        const double inverse = 1.0 / std::hypot(x, y);
        worstRadial = std::max({worstRadial, std::abs(vertex.normal[0] - x * inverse),
                                std::abs(vertex.normal[1] - y * inverse)});
        worstVertical = std::max(worstVertical, std::abs(static_cast<double>(vertex.normal[2])));
        const gpu::HullVertex& flatVertex = flatPrism.vertices()[i];
        worstFlatGap = std::max(worstFlatGap,
                                std::abs(static_cast<double>(vertex.normal[0]) -
                                         static_cast<double>(flatVertex.normal[0])));
    }
    expectTrue(label("there were ring vertices to check", static_cast<double>(ringSamples)),
               ringSamples >= 4 * columns);
    expectTrue(label("the smooth normal on a prism's side is exactly radial, worst", worstRadial),
               worstRadial < 2e-6);
    expectTrue(label("and has no vertical component, worst", worstVertical), worstVertical < 2e-6);
    // Guard: on a fine enough prism the flat and smooth normals nearly agree, so
    // "exactly radial" would be nearly satisfied by doing nothing. The half-facet
    // angle here is 5.6 degrees, which is a 0.005 gap in the x component -- small,
    // but it has to be there.
    expectTrue(label("and it is not merely the face normal, worst gap", worstFlatGap),
               worstFlatGap > 1e-3);

    // The weighting, which nothing above can see: on the prism every triangle at a
    // shared vertex has the same area, so area-weighted and plain averaging agree
    // exactly. A mutation replacing the weighted sum with a plain one passed the
    // whole suite, which is what this is here for. The ridge has one side eight
    // times the area of the other with the same triangle count on each, so the two
    // answers are twenty degrees apart.
    const double tilt = 50.0 * sim::kDegToRad;
    const sim::TriMesh ridge = makeRidge(8.0, 1.0, tilt);
    gpu::HullShading ridgeShading;
    ridgeShading.normals = gpu::HullNormals::Smooth;
    // Wide enough that the fifty-degree ridge smooths rather than creasing, which
    // is the case this is about.
    ridgeShading.creaseAngle = 80.0 * sim::kDegToRad;
    gpu::SceneMesh ridgeMesh;
    ridgeMesh.appendMesh(ridge, sim::Mat3::identity(), {0, 0, 0}, grey, ridgeShading);

    // The two facet normals, taken from the mesh rather than assumed.
    const sim::Vec3 wideNormal = sim::normalize(
        cross(ridge.verts[ridge.tris[0].b] - ridge.verts[ridge.tris[0].a],
              ridge.verts[ridge.tris[0].c] - ridge.verts[ridge.tris[0].a]));
    const sim::Vec3 narrowNormal = sim::normalize(
        cross(ridge.verts[ridge.tris[4].b] - ridge.verts[ridge.tris[4].a],
              ridge.verts[ridge.tris[4].c] - ridge.verts[ridge.tris[4].a]));
    const sim::Vec3 weighted = sim::normalize(wideNormal * 8.0 + narrowNormal * 1.0);
    const sim::Vec3 unweighted = sim::normalize(wideNormal + narrowNormal);
    expectTrue(label("area weighting and plain averaging differ here, dot",
                     dot(weighted, unweighted)),
               dot(weighted, unweighted) < 0.96);

    long long ridgeSamples = 0;
    double worstRidge = 0.0, bestUnweighted = 1e30;
    for (const gpu::HullVertex& vertex : ridgeMesh.vertices()) {
        // The ridge itself: y = 0 and z = 0, shared by both strips.
        if (std::abs(vertex.position[1]) > 1e-6 || std::abs(vertex.position[2]) > 1e-6) continue;
        ++ridgeSamples;
        const sim::Vec3 n{vertex.normal[0], vertex.normal[1], vertex.normal[2]};
        worstRidge = std::max(worstRidge, length(n - weighted));
        bestUnweighted = std::min(bestUnweighted, length(n - unweighted));
    }
    expectTrue(label("there were ridge vertices to check", static_cast<double>(ridgeSamples)),
               ridgeSamples >= 4);
    expectTrue(label("the smooth normal at the ridge is the area-weighted average, worst",
                     worstRidge),
               worstRidge < 1e-6);
    expectTrue(label("and is not the unweighted one, nearest", bestUnweighted),
               bestUnweighted > 0.1);
}

// --- The two closed forms the light has ---------------------------------------

// A face pointing exactly at the sun, and one exactly ninety degrees off it. Both
// are stated before the render; the second is the stronger, because a surface with
// no light on it at all is ambient and nothing else, and that is one multiplication.
void testFaceAtTheSunIsFullyLitAndAtNinetyDegreesIsAmbientOnly() {
    gpu::Device device;
    gpu::HullRenderer renderer;
    if (!setup(device, renderer)) return;

    gpu::MaterialLibrary library;
    std::string error;
    if (!testLibrary(library, error)) return;

    // A box hull, so every face normal is exactly an axis and the two angles are
    // exact rather than nearly right. Driven through the ship path -- paint scheme
    // and rigid-body state included -- rather than through appendMesh, because that
    // is the path a ship actually takes.
    sim::Ship ship;
    ship.hull = sim::makeBox({-30, -9, 0}, {30, 9, 12});
    ship.state.position = {0, 0, 0};
    ship.state.orientation = {};

    gpu::HullPaint paint;
    paint.underwater = paint.bootTopping = paint.topside = paint.deck = "test_grey";
    gpu::HullShading shading;
    shading.normals = gpu::HullNormals::Flat;

    gpu::SceneMesh mesh;
    expectTrue("the box hull builds: " + error,
               mesh.appendShip(ship, paint, library, shading, error));

    gpu::SceneView view;
    view.mode = gpu::HullShadingMode::Shaded;
    // Straight up: the deck's normal is exactly the sun direction, and every side
    // of the box is exactly ninety degrees off it.
    view.sunDirection[0] = 0.0f;
    view.sunDirection[1] = 0.0f;
    view.sunDirection[2] = 1.0f;
    view.sunColour[0] = view.sunColour[1] = view.sunColour[2] = 2.10f;
    view.skyColour[0] = view.skyColour[1] = view.skyColour[2] = 0.30f;
    view.exposure = 1.0f;

    // Oblique, so the deck and one side are both visible at once and the two
    // predictions come out of the same frame under the same light.
    const sim::Vec3 eye{62.0, 58.0, 46.0};
    const sim::Mat4 mvp = camera(eye, {0, 0, 6}, 100.0, view);

    float matrix[16];
    mvp.toFloats(matrix);
    expectTrue("the eye the specular uses is the camera the matrix was built from",
               gpu::eyeAgreesWithCamera(matrix, view.eye, 100.0));

    const float clear[4] = {0.02f, 0.03f, 0.05f, 1.0f};
    Image image;
    expectTrue("the box hull renders",
               renderThroughPng(renderer, mvp, view, mesh, library, clear, "hull_box_lit.png",
                                image));
    if (!image.valid() || image.width != kWidth) return;

    const Reference grey =
        referenceOf(library[static_cast<std::size_t>(library.find("test_grey"))]);
    const sim::Vec3 sun = sunOf(view);

    struct Probe {
        const char* what;
        sim::Vec3 point;   // a point on the face, world frame
        sim::Vec3 normal;  // that face's normal
    };
    const Probe probes[] = {
        {"the deck, whose normal is exactly the sun direction", {6.0, 3.0, 12.0}, {0, 0, 1}},
        {"the port side, exactly ninety degrees off the sun", {6.0, 9.0, 6.0}, {0, 1, 0}},
        {"the bow face, also exactly ninety degrees off", {30.0, 3.0, 6.0}, {1, 0, 0}},
    };

    for (const Probe& probe : probes) {
        std::uint32_t column = 0, row = 0;
        expectTrue(std::string("the probe on ") + probe.what + " projects into the frame",
                   pixelOf(mvp, probe.point, column, row));
        if (!pixelOf(mvp, probe.point, column, row)) continue;

        const sim::Vec3 toEye = normalize(eye - probe.point);
        expectNear(std::string("the probe's normal really is ninety degrees off, on ") + probe.what,
                   dot(probe.normal, sun), probe.normal.z, 0.0);

        int predicted[3];
        referenceCode(grey, probe.normal, sun, toEye, view, predicted);
        const std::uint8_t* pixel = image.pixel(column, row);
        for (int c = 0; c < 3; ++c)
            expectTrue(label((std::string("predicted pixel on ") + probe.what +
                              ", channel error")
                                 .c_str(),
                             std::abs(int{pixel[c]} - predicted[c])),
                       std::abs(int{pixel[c]} - predicted[c]) <= 1);
    }

    // Stated as its own closed form rather than left inside the loop, because it
    // is the strongest single line here: a face with n.l = 0 has no light on it,
    // so its colour is base * sky * (0.5 + 0.5 n.z) and nothing else. On a vertical
    // side n.z is zero too, which halves it exactly.
    const double sideLinear = grey.base[0] * static_cast<double>(view.skyColour[0]) * 0.5;
    const auto sideCode = static_cast<int>(std::lround(std::clamp(sideLinear, 0.0, 1.0) * 255.0));
    std::uint32_t column = 0, row = 0;
    if (pixelOf(mvp, probes[1].point, column, row)) {
        const std::uint8_t* pixel = image.pixel(column, row);
        expectEqual("a face ninety degrees from the sun is exactly half the ambient",
                    static_cast<long long>(pixel[0]), sideCode);
    }

    // And a guard: the two would agree if the light were not doing anything.
    std::uint32_t deckColumn = 0, deckRow = 0;
    if (pixelOf(mvp, probes[0].point, deckColumn, deckRow) &&
        pixelOf(mvp, probes[1].point, column, row)) {
        const int deck = image.pixel(deckColumn, deckRow)[0];
        const int side = image.pixel(column, row)[0];
        std::printf("     box hull: deck %d, side %d (ambient-only closed form %d)\n", deck, side,
                    sideCode);
        expectTrue(label("the lit face is far brighter than the unlit one, by codes",
                         static_cast<double>(deck - side)),
                   deck - side > 60);
    }
}

// Lambert's cosine law, swept. The specular lobe is achromatic for a dielectric,
// so the difference between two colour channels removes it exactly -- and what is
// left is (base_r - base_b) * (sky + (sun / pi) cos(theta)), a straight line in
// cos(theta) whose slope and intercept are both known in advance.
void testLambertCosineLawOverASweep() {
    gpu::Device device;
    gpu::HullRenderer renderer;
    if (!setup(device, renderer)) return;

    gpu::MaterialLibrary library;
    std::string error;
    if (!testLibrary(library, error)) return;
    const int warmIndex = library.find("test_warm");
    const gpu::Material& warmMaterial = library[static_cast<std::size_t>(warmIndex)];
    const Reference warm = referenceOf(warmMaterial);

    // A plate at z = 0 with its normal exactly +z, seen from straight above so the
    // eye vector at the plate's centre is exactly +z as well.
    gpu::SceneMesh mesh;
    gpu::HullShading shading;
    shading.normals = gpu::HullNormals::Flat;
    mesh.appendMesh(makeQuad({0, 0, 0}, {14, 0, 0}, {0, 14, 0}), sim::Mat3::identity(), {0, 0, 0},
                    static_cast<std::uint32_t>(warmIndex), shading);

    gpu::SceneView view;
    // Equal in every channel, which is what makes the channel difference cancel the
    // sun's own colour along with the specular.
    view.sunColour[0] = view.sunColour[1] = view.sunColour[2] = 2.40f;
    view.skyColour[0] = view.skyColour[1] = view.skyColour[2] = 0.25f;
    view.exposure = 1.0f;

    const sim::Vec3 eye{0, 0, 70};
    const sim::Mat4 mvp = camera(eye, {0, 0, 0}, 70.0, view);
    std::uint32_t column = 0, row = 0;
    expectTrue("the plate's centre projects into the frame", pixelOf(mvp, {0, 0, 0}, column, row));
    if (!pixelOf(mvp, {0, 0, 0}, column, row)) return;

    const sim::Vec3 normal{0, 0, 1};
    const sim::Vec3 toEye{0, 0, 1};
    const double baseGap = warm.base[0] - warm.base[2];
    expectTrue(label("the sweep's material has channels far enough apart to divide by", baseGap),
               baseGap > 0.4);

    const float clear[4] = {0.0f, 0.0f, 0.0f, 1.0f};
    long long samples = 0, exceeded = 0, cosineExceeded = 0;
    long long squaredRejected = 0, linearRejected = 0;
    double worstCosineError = 0.0, worstChannelError = 0.0;

    std::printf("     Lambert sweep (theta, predicted rgb, measured rgb, cos from the"
                " channel difference):\n");
    for (int step = 0; step <= 19; ++step) {
        const double theta = 5.0 * step * sim::kDegToRad;  // 0 .. 95 degrees
        view.sunDirection[0] = static_cast<float>(std::sin(theta));
        view.sunDirection[1] = 0.0f;
        view.sunDirection[2] = static_cast<float>(std::cos(theta));
        const sim::Vec3 sun = sunOf(view);

        Image image;
        char name[64];
        std::snprintf(name, sizeof name, "hull_lambert_%02d.png", step);
        if (!renderThroughPng(renderer, mvp, view, mesh, library, clear, name, image)) {
            expectTrue("the sweep renders", false);
            return;
        }
        const std::uint8_t* pixel = image.pixel(column, row);

        int predicted[3];
        referenceCode(warm, normal, sun, toEye, view, predicted);
        ++samples;
        for (int c = 0; c < 3; ++c) {
            const double gap = std::abs(int{pixel[c]} - predicted[c]);
            worstChannelError = std::max(worstChannelError, gap);
            if (gap > 1.0) ++exceeded;
        }

        // The specular lobe cancels here rather than being modelled: F0 is 0.04 in
        // every channel for a dielectric and D and Vis carry no colour at all, so
        // the red-minus-blue difference is purely diffuse plus sky.
        const double sky = static_cast<double>(view.skyColour[0]) * (0.5 + 0.5 * normal.z);
        const double slope = static_cast<double>(view.sunColour[0]) / sim::kPi;
        const double measured =
            (static_cast<double>(pixel[0]) - static_cast<double>(pixel[2])) / 255.0;
        const double recovered = (measured / baseGap - sky) / slope;
        const double wanted = std::max(std::cos(theta), 0.0);
        // Two codes of rounding across two channels, divided through by the same
        // constants the recovery divides by.
        const double tolerance = 2.0 / 255.0 / baseGap / slope;
        worstCosineError = std::max(worstCosineError, std::abs(recovered - wanted));
        if (std::abs(recovered - wanted) > tolerance) ++cosineExceeded;

        // Negative controls, because "matches a curve" is worth nothing unless
        // some other curve does not. Both of these are darker as theta grows, so
        // an assertion that only said "darker" would accept either.
        if (std::abs(recovered - wanted * wanted) > tolerance) ++squaredRejected;
        if (std::abs(recovered - std::max(0.0, 1.0 - theta / (0.5 * sim::kPi))) > tolerance)
            ++linearRejected;

        if (step % 4 == 0)
            std::printf("       %5.1f deg  predicted %3d %3d %3d  measured %3d %3d %3d "
                        " cos %.4f (exact %.4f)\n",
                        theta * sim::kRadToDeg, predicted[0], predicted[1], predicted[2],
                        pixel[0], pixel[1], pixel[2], recovered, wanted);
    }

    expectEqual("the sweep ran over nineteen angles plus zero", samples, 20);
    expectEqual("every sampled pixel matches the BRDF written out from the formula", exceeded, 0);
    expectEqual("and the diffuse term follows cos(theta) at every one of them", cosineExceeded, 0);
    expectTrue(label("worst cosine discrepancy over the sweep", worstCosineError),
               worstCosineError < 0.02);
    expectTrue(label("worst channel discrepancy, codes", worstChannelError),
               worstChannelError <= 1.0);
    // A cos^2 falloff and a linear one are both monotone and both wrong. If either
    // survived the sweep, the sweep is not measuring the cosine law.
    expectTrue(label("a cos^2 falloff is rejected at, samples",
                     static_cast<double>(squaredRejected)),
               squaredRejected >= 14);
    expectTrue(label("a linear falloff is rejected at, samples",
                     static_cast<double>(linearRejected)),
               linearRejected >= 14);
    std::printf("     cos law: worst %.4f over 20 angles; cos^2 rejected by %lld,"
                " linear by %lld\n",
                worstCosineError, squaredRejected, linearRejected);
}

// Two materials, one geometry, one light, one frame. The assertion is the
// predicted value of each, not that they differ.
void testTwoMaterialsDifferPredictably() {
    gpu::Device device;
    gpu::HullRenderer renderer;
    if (!setup(device, renderer)) return;

    gpu::MaterialLibrary library;
    std::string error;
    if (!testLibrary(library, error)) return;
    const int warmIndex = library.find("test_warm");
    const int antifoulingIndex = library.find("antifouling");
    const int steelIndex = library.find("bare_steel");

    gpu::HullShading shading;
    shading.normals = gpu::HullNormals::Flat;
    gpu::SceneMesh mesh;
    // Three coplanar plates side by side: same normal, same light, same eye
    // direction to within the frustum, different rows of the material table.
    const sim::Vec3 centres[3] = {{-18, 0, 0}, {0, 0, 0}, {18, 0, 0}};
    const int indices[3] = {warmIndex, antifoulingIndex, steelIndex};
    for (int i = 0; i < 3; ++i)
        mesh.appendMesh(makeQuad(centres[i], {7, 0, 0}, {0, 7, 0}), sim::Mat3::identity(),
                        {0, 0, 0}, static_cast<std::uint32_t>(indices[i]), shading);

    gpu::SceneView view;
    // Well off the mirror direction. The first version put the sun near it, and the
    // bare steel -- a conductor with F0 = 0.56 and a narrow lobe -- came back at
    // 255 in every channel. Physically right, and a useless comparison: two clamped
    // values agree whatever either of them was.
    view.sunDirection[0] = 0.62f;
    view.sunDirection[1] = 0.10f;
    view.sunDirection[2] = 0.78f;
    view.sunColour[0] = 2.30f;
    view.sunColour[1] = 2.24f;
    view.sunColour[2] = 2.10f;
    view.skyColour[0] = 0.22f;
    view.skyColour[1] = 0.26f;
    view.skyColour[2] = 0.33f;

    const sim::Vec3 eye{0, -6, 62};
    const sim::Mat4 mvp = camera(eye, {0, 0, 0}, 62.0, view);
    const float clear[4] = {0.01f, 0.01f, 0.02f, 1.0f};

    Image image;
    expectTrue("the three plates render",
               renderThroughPng(renderer, mvp, view, mesh, library, clear, "hull_materials.png",
                                image));
    if (!image.valid() || image.width != kWidth) return;

    const sim::Vec3 sun = sunOf(view);
    int codes[3][3]{};
    bool sampled[3] = {false, false, false};
    for (int i = 0; i < 3; ++i) {
        std::uint32_t column = 0, row = 0;
        if (!pixelOf(mvp, centres[i], column, row)) continue;
        sampled[i] = true;
        const Reference material = referenceOf(library[static_cast<std::size_t>(indices[i])]);
        const sim::Vec3 toEye = normalize(eye - centres[i]);
        int predicted[3];
        referenceCode(material, {0, 0, 1}, sun, toEye, view, predicted);
        const std::uint8_t* pixel = image.pixel(column, row);
        for (int c = 0; c < 3; ++c) {
            codes[i][c] = pixel[c];
            expectTrue(label((library[static_cast<std::size_t>(indices[i])].name +
                              " matches its predicted pixel, channel error")
                                 .c_str(),
                             std::abs(int{pixel[c]} - predicted[c])),
                       std::abs(int{pixel[c]} - predicted[c]) <= 1);
            // A clamped channel agrees with anything that also clamped, so the
            // agreement above is only evidence while nothing is at a rail.
            expectTrue(library[static_cast<std::size_t>(indices[i])].name +
                           " is not clipped, so the agreement means something",
                       pixel[c] > 0 && pixel[c] < 255);
        }
    }
    expectTrue("all three plates were sampled", sampled[0] && sampled[1] && sampled[2]);
    std::printf("     same geometry and light: warm %3d %3d %3d, antifouling %3d %3d %3d,"
                " bare steel %3d %3d %3d\n",
                codes[0][0], codes[0][1], codes[0][2], codes[1][0], codes[1][1], codes[1][2],
                codes[2][0], codes[2][1], codes[2][2]);

    // Guard: matching three predictions proves nothing if the three predictions are
    // the same number. They are not, and they are not for reasons that are physical
    // -- antifouling is dark and matt, bare steel is a conductor and has no diffuse
    // lobe at all.
    expectTrue("antifouling is far darker than the warm dielectric",
               codes[0][0] - codes[1][0] > 40);
    expectTrue("the conductor's red and blue are much closer together than the warm one's",
               std::abs(codes[2][0] - codes[2][2]) < std::abs(codes[0][0] - codes[0][2]) / 2);
}

// Nothing is culled, because a hull is seen from inside once it is cut away or
// flooded. A back face is a real surface there, and it must be lit as the face it
// actually presents rather than as the one it was wound as. Nothing else in this
// file looks at a hull from inside, so a mutation removing the flip passed the
// whole suite -- this is that gap closed.
void testTheFarSideOfASurfaceIsLitAsTheFaceItPresents() {
    gpu::Device device;
    gpu::HullRenderer renderer;
    if (!setup(device, renderer)) return;

    gpu::MaterialLibrary library;
    std::string error;
    if (!testLibrary(library, error)) return;
    const int warmIndex = library.find("test_warm");
    const Reference warm = referenceOf(library[static_cast<std::size_t>(warmIndex)]);

    gpu::HullShading shading;
    shading.normals = gpu::HullNormals::Flat;
    gpu::SceneMesh mesh;
    // Wound so its geometric normal is +z; the camera is underneath it.
    mesh.appendMesh(makeQuad({0, 0, 0}, {12, 0, 0}, {0, 12, 0}), sim::Mat3::identity(),
                    {0, 0, 0}, static_cast<std::uint32_t>(warmIndex), shading);
    expectNear("the quad is wound normal-up", mesh.vertices()[0].normal[2], 1.0, 1e-6);

    gpu::SceneView view;
    view.sunDirection[0] = 0.24f;
    view.sunDirection[1] = 0.10f;
    // Below the plate, so the face the camera sees has the light on it.
    view.sunDirection[2] = -0.96f;
    view.sunColour[0] = view.sunColour[1] = view.sunColour[2] = 2.20f;
    view.skyColour[0] = view.skyColour[1] = view.skyColour[2] = 0.28f;

    const sim::Vec3 eye{4.0, 5.0, -40.0};
    const sim::Mat4 mvp = camera(eye, {0, 0, 0}, 40.0, view);
    const float clear[4] = {0.0f, 0.0f, 0.0f, 1.0f};
    Image image;
    expectTrue("the far side of the plate renders",
               renderThroughPng(renderer, mvp, view, mesh, library, clear, "hull_backface.png",
                                image));
    if (!image.valid() || image.width != kWidth) return;

    std::uint32_t column = 0, row = 0;
    if (!pixelOf(mvp, {0, 0, 0}, column, row)) {
        expectTrue("the plate's centre projects into the frame", false);
        return;
    }
    const sim::Vec3 sun = sunOf(view);
    const sim::Vec3 toEye = sim::normalize(eye);

    int presented[3], wound[3];
    referenceCode(warm, {0, 0, -1}, sun, toEye, view, presented);  // the face it presents
    referenceCode(warm, {0, 0, 1}, sun, toEye, view, wound);       // the face it was wound as
    const std::uint8_t* pixel = image.pixel(column, row);
    std::printf("     back face: measured %3d %3d %3d, presented-normal prediction %3d %3d %3d,"
                " wound-normal prediction %3d %3d %3d\n",
                pixel[0], pixel[1], pixel[2], presented[0], presented[1], presented[2], wound[0],
                wound[1], wound[2]);

    for (int c = 0; c < 3; ++c)
        expectTrue(label("the far side is lit by the normal it presents, channel error",
                         std::abs(int{pixel[c]} - presented[c])),
                   std::abs(int{pixel[c]} - presented[c]) <= 1);
    // Guard: the two predictions must actually differ, or the check is empty.
    expectTrue(label("the two predictions differ, codes",
                     static_cast<double>(presented[0] - wound[0])),
               std::abs(presented[0] - wound[0]) > 40);
}

// At a grazing view Schlick's term stops being negligible: it is what makes a wet
// deck flare along the horizon. Every other probe here sits near the surface
// normal, where (1 - v.h)^5 and (1 - v.h)^4 are both a rounding error apart -- a
// mutation changing the exponent passed the whole suite. This is the geometry
// where it does not.
void testGrazingViewFollowsSchlicksFresnel() {
    gpu::Device device;
    gpu::HullRenderer renderer;
    if (!setup(device, renderer)) return;

    gpu::MaterialLibrary library;
    std::string error;
    if (!testLibrary(library, error)) return;
    const int warmIndex = library.find("test_warm");
    const Reference warm = referenceOf(library[static_cast<std::size_t>(warmIndex)]);

    gpu::HullShading shading;
    shading.normals = gpu::HullNormals::Flat;
    gpu::SceneMesh mesh;
    mesh.appendMesh(makeQuad({0, 0, 0}, {40, 0, 0}, {0, 40, 0}), sim::Mat3::identity(), {0, 0, 0},
                    static_cast<std::uint32_t>(warmIndex), shading);

    gpu::SceneView view;
    // The sun on the far side from the camera, which is what pushes v.h down.
    view.sunDirection[0] = 0.0f;
    view.sunDirection[1] = -0.60f;
    view.sunDirection[2] = 0.80f;
    view.sunColour[0] = view.sunColour[1] = view.sunColour[2] = 2.40f;
    view.skyColour[0] = view.skyColour[1] = view.skyColour[2] = 0.24f;

    // Eight metres up at sixty out: n.v is 0.13, about eight degrees above the
    // plate.
    const sim::Vec3 eye{0, 60, 8};
    const sim::Mat4 mvp = camera(eye, {0, 0, 0}, 60.0, view);
    const sim::Vec3 sun = sunOf(view);
    const sim::Vec3 toEye = sim::normalize(eye);
    const sim::Vec3 half = sim::normalize(sun + toEye);
    std::printf("     grazing: n.v %.4f, v.h %.4f -- Schlick's (1 - v.h)^5 is %.4f against"
                " ^4 of %.4f\n",
                toEye.z, dot(toEye, half), std::pow(1.0 - dot(toEye, half), 5.0),
                std::pow(1.0 - dot(toEye, half), 4.0));
    expectTrue(label("the view really is grazing, n.v", toEye.z), toEye.z < 0.2);
    expectTrue(label("and v.h is far enough from one for the exponent to matter",
                     dot(toEye, half)),
               dot(toEye, half) < 0.7);

    const float clear[4] = {0.0f, 0.0f, 0.0f, 1.0f};
    Image image;
    expectTrue("the grazing view renders",
               renderThroughPng(renderer, mvp, view, mesh, library, clear, "hull_grazing.png",
                                image));
    if (!image.valid() || image.width != kWidth) return;

    std::uint32_t column = 0, row = 0;
    if (!pixelOf(mvp, {0, 0, 0}, column, row)) {
        expectTrue("the sampled point projects into the frame", false);
        return;
    }
    int predicted[3];
    referenceCode(warm, {0, 0, 1}, sun, toEye, view, predicted);
    const std::uint8_t* pixel = image.pixel(column, row);
    for (int c = 0; c < 3; ++c) {
        expectTrue(label("the grazing pixel matches the BRDF, channel error",
                         std::abs(int{pixel[c]} - predicted[c])),
                   std::abs(int{pixel[c]} - predicted[c]) <= 1);
        expectTrue("the grazing pixel is not clipped", pixel[c] > 0 && pixel[c] < 255);
    }

    // The instrument: recompute the same pixel with Schlick's exponent set to four
    // and require it to be visibly different. If it were not, this geometry would
    // not be testing the exponent at all.
    const double vdh = dot(toEye, half);
    const double fresnelFive = 0.04 + 0.96 * std::pow(1.0 - vdh, 5.0);
    const double fresnelFour = 0.04 + 0.96 * std::pow(1.0 - vdh, 4.0);
    double specular[3];
    referenceSpecular(warm, {0, 0, 1}, sun, toEye, specular);
    const double lobe = specular[0] / fresnelFive;  // D * Vis alone
    const double codesApart = 255.0 * lobe * (fresnelFour - fresnelFive) * std::max(sun.z, 0.0) *
                              static_cast<double>(view.sunColour[0]);
    std::printf("     grazing: a fourth-power Schlick term would move this pixel by"
                " %.1f codes\n", codesApart);
    expectTrue(label("this geometry can tell Schlick's exponent apart, codes", codesApart),
               codesApart > 4.0);
}

// --- Depth ---------------------------------------------------------------------

// Drawing near-then-far gives the right picture even with no depth buffer at all,
// so the far surface is appended *second* and must still lose. The material channel
// turns "which one won" into an integer, and the depth channel says by how much.
void testTheNearerSurfaceWinsWhicheverIsDrawnSecond() {
    gpu::Device device;
    gpu::HullRenderer renderer;
    if (!setup(device, renderer)) return;

    gpu::MaterialLibrary library;
    std::string error;
    if (!testLibrary(library, error)) return;
    const auto warm = static_cast<std::uint32_t>(library.find("test_warm"));
    const auto cool = static_cast<std::uint32_t>(library.find("test_cool"));

    gpu::HullShading shading;
    shading.normals = gpu::HullNormals::Flat;
    gpu::SceneMesh mesh;
    // Both plates face the camera and both fill the frame. Near first, far second.
    mesh.appendMesh(makeQuad({0, 0, 10}, {40, 0, 0}, {0, 40, 0}), sim::Mat3::identity(),
                    {0, 0, 0}, warm, shading);
    mesh.appendMesh(makeQuad({0, 0, -10}, {40, 0, 0}, {0, 40, 0}), sim::Mat3::identity(),
                    {0, 0, 0}, cool, shading);

    gpu::SceneView view;
    const sim::Vec3 eye{0, 0, 90};
    const sim::Mat4 mvp = camera(eye, {0, 0, 0}, 90.0, view);
    const float clear[4] = {0.0f, 0.0f, 0.0f, 1.0f};

    view.mode = gpu::HullShadingMode::MaterialId;
    Image identifiers;
    expectTrue("the depth scene renders its material channel",
               renderThroughPng(renderer, mvp, view, mesh, library, clear, "hull_depth_id.png",
                                identifiers));
    if (!identifiers.valid() || identifiers.width != kWidth) return;

    std::size_t nearPixels = 0, farPixels = 0, background = 0;
    for (std::uint32_t y = 0; y < identifiers.height; ++y)
        for (std::uint32_t x = 0; x < identifiers.width; ++x) {
            std::uint32_t material = 0;
            if (!gpu::decodeMaterialId(identifiers.pixel(x, y), material)) { ++background; continue; }
            if (material == warm) ++nearPixels;
            else if (material == cool) ++farPixels;
        }
    expectTrue(label("the near plate covers a large part of the frame",
                     static_cast<double>(nearPixels)),
               static_cast<double>(nearPixels) > 0.5 * kWidth * kHeight);
    expectEqual("the farther surface, drawn second, is completely hidden",
                static_cast<long long>(farPixels), 0);
    std::printf("     depth: near %zu px, far %zu px, background %zu px\n", nearPixels, farPixels,
                background);

    // And the depth actually written is the near plate's, to the resolution of the
    // channel. Predicted from the camera matrix, not read off the image.
    view.mode = gpu::HullShadingMode::Depth;
    Image depths;
    expectTrue("the depth channel renders",
               renderThroughPng(renderer, mvp, view, mesh, library, clear, "hull_depth.png",
                                depths));
    if (!depths.valid() || depths.width != kWidth) return;

    const double nearClip = clipDepthOf(mvp, {0, 0, 10});
    const double farClip = clipDepthOf(mvp, {0, 0, -10});
    // Stated in channel codes rather than in clip units, because that is the
    // resolution the question is asked at -- and because twenty metres of world
    // separation is only 0.0025 of the [0, 1] clip range here. A far plane at
    // eight times the eye distance puts the entire scene in the last one per cent
    // of the depth buffer, which is a real property of a perspective projection
    // and is precisely why this channel carries sixteen bits rather than eight.
    const double separationCodes = (farClip - nearClip) * 65535.0;
    std::printf("     depth: near %.6f, far %.6f -- %.0f codes apart of 65535\n", nearClip,
                farClip, separationCodes);
    expectTrue(label("the two plates are far enough apart to resolve, channel codes",
                     separationCodes),
               separationCodes > 100.0);

    std::uint32_t column = 0, row = 0;
    if (pixelOf(mvp, {0, 0, 10}, column, row)) {
        double decoded = 0.0;
        expectTrue("the centre pixel carries a depth",
                   gpu::decodeSceneDepth(depths.pixel(column, row), decoded));
        expectNear("the depth written is the near plate's, from the camera matrix", decoded,
                   nearClip, 2.0 / 65535.0);
        expectTrue(label("and is not the far plate's, away by codes",
                         (farClip - decoded) * 65535.0),
                   (farClip - decoded) * 65535.0 > 50.0);
    }
}

// --- Composition: the ship is in the sea, not beside it -----------------------

void testShipOccludesTheSeaBehindItAndNotInFrontOfIt() {
    gpu::Device device;
    gpu::HullRenderer renderer;
    if (!setup(device, renderer)) return;

    gpu::MaterialLibrary library;
    std::string error;
    if (!testLibrary(library, error)) return;
    const auto seaMaterial = static_cast<std::uint32_t>(library.find("sea_water"));

    sim::Ship ship = makeTestShip();
    expectTrue("the test hull is a closed, consistently wound manifold",
               sim::isClosedManifold(ship.hull));
    ship.initialise(sim::Sea{0.0});

    // Dead calm, so the sea is exactly the plane z = 0 and every prediction below is
    // exact rather than nearly right.
    sim::SeaState calm;
    calm.significantHeight = 0.0;
    const sim::WaveField field(calm);
    gpu::OceanSurface surface;
    surface.build(field, gpu::OceanGrid{0.0, 0.0, 260.0, 96}, 0.0);

    gpu::HullPaint paint;
    paint.waterlineZ = 5.0;
    paint.deckZ = 10.5;
    gpu::HullShading shading;
    shading.normals = gpu::HullNormals::Smooth;

    gpu::SceneMesh mesh;
    expectTrue("the ship builds: " + error, mesh.appendShip(ship, paint, library, shading, error));
    const std::size_t hullVertices = mesh.vertices().size();
    // The sea is appended **after** the hull, so a broken depth test paints the
    // water straight over the ship rather than hiding the mistake behind the draw
    // order.
    mesh.appendOcean(surface, seaMaterial);
    expectTrue("the sea was appended after the hull", mesh.vertices().size() > hullVertices);

    gpu::SceneView view;
    view.mode = gpu::HullShadingMode::MaterialId;
    const sim::Vec3 eye{0, 165, 42};
    const sim::Mat4 mvp = camera(eye, {0, 0, 2}, 170.0, view);
    const float clear[4] = {0.0f, 0.0f, 0.0f, 1.0f};

    Image identifiers;
    expectTrue("the composed scene renders its material channel",
               renderThroughPng(renderer, mvp, view, mesh, library, clear, "hull_sea_id.png",
                                identifiers));
    if (!identifiers.valid() || identifiers.width != kWidth) return;

    view.mode = gpu::HullShadingMode::Depth;
    Image depths;
    expectTrue("and its depth channel",
               renderThroughPng(renderer, mvp, view, mesh, library, clear, "hull_sea_depth.png",
                                depths));
    if (!depths.valid() || depths.width != kWidth) return;

    std::size_t seaPixels = 0, hullPixels = 0;
    for (std::uint32_t y = 0; y < identifiers.height; ++y)
        for (std::uint32_t x = 0; x < identifiers.width; ++x) {
            std::uint32_t material = 0;
            if (!gpu::decodeMaterialId(identifiers.pixel(x, y), material)) continue;
            if (material == seaMaterial) ++seaPixels;
            else ++hullPixels;
        }
    expectTrue(label("there is a sea in the frame, pixels", static_cast<double>(seaPixels)),
               static_cast<double>(seaPixels) > 0.25 * kWidth * kHeight);
    expectTrue(label("and a ship in it, pixels", static_cast<double>(hullPixels)),
               static_cast<double>(hullPixels) > 0.03 * kWidth * kHeight);

    // A point on the sea surface directly beyond the ship. The camera is at y =
    // +165, so this one is on the far side and the ray to it passes through the
    // hull -- which is a statement about the geometry, so the test checks it rather
    // than assuming it.
    const sim::Vec3 behind{0.0, -34.0, 0.0};
    const sim::Vec3 ahead{0.0, 62.0, 0.0};
    expectTrue("the far sea point really is beyond the ship's own extent",
               behind.y < ship.hullLo.y + ship.state.position.y);
    expectTrue("and the near one is short of it",
               ahead.y > ship.hullHi.y + ship.state.position.y);

    std::uint32_t column = 0, row = 0;
    if (pixelOf(mvp, behind, column, row)) {
        std::uint32_t material = 0;
        double drawn = 0.0;
        const double seaDepth = clipDepthOf(mvp, behind);
        expectTrue("the occluded sea point's pixel carries geometry",
                   gpu::decodeMaterialId(identifiers.pixel(column, row), material));
        expectTrue("the ship occludes the sea behind it", material != seaMaterial);
        expectTrue("that pixel carries a depth",
                   gpu::decodeSceneDepth(depths.pixel(column, row), drawn));
        // The stronger form, and the one that does not depend on knowing which
        // material the hull is painted there: whatever was drawn is *nearer* than
        // the water would have been.
        expectTrue(label("and what is drawn there is nearer than the water, by", seaDepth - drawn),
                   drawn < seaDepth - 4.0 / 65535.0);
    } else {
        expectTrue("the occluded sea point projects into the frame", false);
    }

    if (pixelOf(mvp, ahead, column, row)) {
        std::uint32_t material = 0;
        double drawn = 0.0;
        expectTrue("the near sea point's pixel carries geometry",
                   gpu::decodeMaterialId(identifiers.pixel(column, row), material));
        expectEqual("the ship does not occlude the sea in front of it",
                    static_cast<long long>(material), static_cast<long long>(seaMaterial));
        expectTrue("that pixel carries a depth",
                   gpu::decodeSceneDepth(depths.pixel(column, row), drawn));
        expectNear("and it is the sea's own depth, from the camera matrix", drawn,
                   clipDepthOf(mvp, ahead), 4.0 / 65535.0);
    } else {
        expectTrue("the near sea point projects into the frame", false);
    }

    // Guard against the whole thing passing on a frame with no sea in it: without
    // the ship, the far point must be water.
    gpu::SceneMesh seaOnly;
    seaOnly.appendOcean(surface, seaMaterial);
    view.mode = gpu::HullShadingMode::MaterialId;
    Image bare;
    if (renderThroughPng(renderer, mvp, view, seaOnly, library, clear, "hull_sea_only.png", bare) &&
        pixelOf(mvp, behind, column, row)) {
        std::uint32_t material = 0;
        expectTrue("with the ship removed, that same pixel is water",
                   gpu::decodeMaterialId(bare.pixel(column, row), material) &&
                       material == seaMaterial);
    }
}

// --- Paint is data, and lands where the scheme says ---------------------------

// A known hull point, projected with sim::clipToPixel, must be drawn in the
// material its z band was painted with. The depth channel confirms the triangle
// predicted is the one actually drawn there, so an occluded sample is skipped
// rather than quietly compared against something else.
void testPaintBandsLandWhereTheSchemeSaysTheyDo() {
    gpu::Device device;
    gpu::HullRenderer renderer;
    if (!setup(device, renderer)) return;

    gpu::MaterialLibrary library;
    std::string error;
    if (!testLibrary(library, error)) return;

    sim::Ship ship = makeTestShip();
    ship.state.position = {0, 0, 0};
    ship.state.orientation = {};

    gpu::HullPaint paint;
    paint.waterlineZ = 5.0;
    paint.bootTopDepth = 0.4;
    paint.bootTopHeight = 1.0;
    paint.deckZ = 10.5;
    const auto underwater = static_cast<std::uint32_t>(library.find(paint.underwater));
    const auto bootTopping = static_cast<std::uint32_t>(library.find(paint.bootTopping));
    const auto topside = static_cast<std::uint32_t>(library.find(paint.topside));
    const auto deck = static_cast<std::uint32_t>(library.find(paint.deck));

    gpu::HullShading shading;
    shading.normals = gpu::HullNormals::Flat;
    gpu::SceneMesh mesh;
    expectTrue("the painted ship builds: " + error,
               mesh.appendShip(ship, paint, library, shading, error));

    gpu::SceneView view;
    const sim::Vec3 eye{28.0, 96.0, 34.0};
    const sim::Mat4 mvp = camera(eye, {0, 0, 4}, 105.0, view);
    const float clear[4] = {0.0f, 0.0f, 0.0f, 1.0f};

    view.mode = gpu::HullShadingMode::MaterialId;
    Image identifiers;
    expectTrue("the painted hull renders its material channel",
               renderThroughPng(renderer, mvp, view, mesh, library, clear, "hull_paint_id.png",
                                identifiers));
    view.mode = gpu::HullShadingMode::Depth;
    Image depths;
    expectTrue("and its depth channel",
               renderThroughPng(renderer, mvp, view, mesh, library, clear, "hull_paint_depth.png",
                                depths));
    if (!identifiers.valid() || !depths.valid() || identifiers.width != kWidth) return;

    long long checked = 0, wrong = 0, skipped = 0, tooSmall = 0;
    long long seen[4] = {0, 0, 0, 0};
    for (const sim::Tri& tri : ship.hull.tris) {
        const sim::Vec3& a = ship.hull.verts[tri.a];
        const sim::Vec3& b = ship.hull.verts[tri.b];
        const sim::Vec3& c = ship.hull.verts[tri.c];
        const sim::Vec3 centroid = (a + b + c) / 3.0;
        const sim::Vec3 normal = normalize(cross(b - a, c - a));

        // The scheme, applied here independently of hull.cpp.
        std::uint32_t wanted = topside;
        int band = 2;
        if (centroid.z >= paint.deckZ && normal.z >= paint.deckNormalZ) { wanted = deck; band = 3; }
        else if (centroid.z < paint.waterlineZ - paint.bootTopDepth) { wanted = underwater; band = 0; }
        else if (centroid.z < paint.waterlineZ + paint.bootTopHeight) {
            wanted = bootTopping;
            band = 1;
        }

        std::uint32_t column = 0, row = 0;
        if (!pixelOf(mvp, centroid, column, row)) continue;
        // Two pixels of clearance. Below that the centroid's pixel is shared with a
        // neighbour, and at a band boundary the neighbour is a different material at
        // a depth a handful of codes away -- which is a limitation of sampling a
        // rendered image at a predicted point, not of the renderer. The first
        // version of this check omitted it and produced one mismatch in 389, on a
        // sliver nine pixels in from the frame edge.
        if (centroidClearanceInPixels(mvp, a, b, c) < 2.0) { ++tooSmall; continue; }
        double drawnDepth = 0.0;
        if (!gpu::decodeSceneDepth(depths.pixel(column, row), drawnDepth)) continue;
        // Only compare where the triangle we predicted is the surface the depth
        // buffer kept. Anything else is occluded, and comparing there would be
        // comparing against a different triangle.
        const double centroidDepth = clipDepthOf(mvp, centroid);
        if (std::abs(drawnDepth - centroidDepth) > 6.0 / 65535.0) { ++skipped; continue; }

        std::uint32_t material = 0;
        if (!gpu::decodeMaterialId(identifiers.pixel(column, row), material)) { ++skipped; continue; }
        ++checked;
        ++seen[band];
        if (material != wanted) {
            ++wrong;
            std::printf("       MISMATCH band %d z %.3f nz %.3f depth gap %.2e px %u,%u"
                        " got %u want %u\n",
                        band, centroid.z, normal.z, drawnDepth - centroidDepth, column, row,
                        material, wanted);
        }
    }

    std::printf("     paint: %lld triangles checked (%lld occluded, %lld too small to sample);"
                " %lld underwater, %lld boot-topping, %lld topside, %lld deck\n",
                checked, skipped, tooSmall, seen[0], seen[1], seen[2], seen[3]);
    expectTrue(label("enough hull triangles were visible to check", static_cast<double>(checked)),
               checked > 120);
    expectEqual("every visible hull triangle is drawn in the material its band was painted",
                wrong, 0);
    // Guards: a scheme that painted the whole ship one colour would satisfy the
    // check above perfectly.
    expectTrue("the underwater band is on the ship", seen[0] > 10);
    expectTrue("the boot-topping band is on the ship", seen[1] > 5);
    expectTrue("the topside band is on the ship", seen[2] > 10);

    // Paint is on the hull, not on the world: heel the ship over and the bands must
    // stay exactly where they were painted. Deciding them from world z instead would
    // slide the waterline stripe around the shell as she rolls, which looks almost
    // right and is entirely wrong.
    sim::Ship heeled = ship;
    heeled.state.orientation = sim::Quat::fromAxisAngle({1, 0, 0}, -25.0 * sim::kDegToRad);
    gpu::SceneMesh heeledMesh;
    expectTrue("the heeled ship builds: " + error,
               heeledMesh.appendShip(heeled, paint, library, shading, error));
    bool sameMaterials = heeledMesh.vertices().size() == mesh.vertices().size();
    for (std::size_t i = 0; sameMaterials && i < mesh.vertices().size(); ++i)
        sameMaterials = heeledMesh.vertices()[i].material == mesh.vertices()[i].material;
    expectTrue("heeling the ship does not repaint her", sameMaterials);
    bool moved = false;
    for (std::size_t i = 0; i < mesh.vertices().size() && !moved; ++i)
        moved = std::abs(heeledMesh.vertices()[i].position[1] -
                         mesh.vertices()[i].position[1]) > 0.5f;
    expectTrue("and the heel did move the hull, so that was not a comparison of two"
               " identical builds", moved);

    // A paint scheme naming a material nobody loaded is a broken ship definition,
    // not a grey default.
    gpu::HullPaint missing = paint;
    missing.bootTopping = "no_such_paint";
    gpu::SceneMesh rejected;
    expectTrue("a paint scheme naming an unknown material is refused",
               !rejected.appendShip(ship, missing, library, shading, error));
    expectTrue("and says which one: " + error, error.find("no_such_paint") != std::string::npos);
}

// --- Symmetry ------------------------------------------------------------------

// A handedness mistake in the camera, the winding or the normals looks perfectly
// fine from one viewpoint. It does not survive being looked at from both sides.
void testPortAndStarboardViewsMirror() {
    gpu::Device device;
    gpu::HullRenderer renderer;
    if (!setup(device, renderer)) return;

    gpu::MaterialLibrary library;
    std::string error;
    if (!testLibrary(library, error)) return;
    const auto seaMaterial = static_cast<std::uint32_t>(library.find("sea_water"));

    sim::Ship ship = makeTestShip();
    ship.initialise(sim::Sea{0.0});
    // Upright and on the centreline. The hull is symmetric by construction, so the
    // only thing that could break the symmetry of the scene is the ship's own
    // attitude -- and a fraction of a degree of heel would show up here as a
    // mismatch that had nothing to do with the renderer.
    const sim::Diagnostics diagnostics = ship.diagnostics(sim::Sea{0.0});
    expectTrue(label("the ship floats upright, heel deg", diagnostics.heelDeg),
               std::abs(diagnostics.heelDeg) < 1e-6);
    ship.state.orientation = {};
    ship.state.position.y = 0.0;

    sim::SeaState calm;
    calm.significantHeight = 0.0;
    const sim::WaveField field(calm);
    gpu::OceanSurface surface;
    surface.build(field, gpu::OceanGrid{0.0, 0.0, 220.0, 64}, 0.0);

    gpu::HullPaint paint;
    paint.waterlineZ = 5.0;
    paint.deckZ = 10.5;
    gpu::HullShading shading;
    shading.normals = gpu::HullNormals::Smooth;
    gpu::SceneMesh mesh;
    if (!mesh.appendShip(ship, paint, library, shading, error)) {
        expectTrue("the ship builds: " + error, false);
        return;
    }
    mesh.appendOcean(surface, seaMaterial);

    gpu::SceneView portView, starboardView;
    // The sun must lie in the ship's centreplane, or the two sides are genuinely
    // lit differently and the images have no business mirroring. This is a
    // constraint on the test, not a limitation of the renderer.
    for (gpu::SceneView* view : {&portView, &starboardView}) {
        view->sunDirection[0] = 0.36f;
        view->sunDirection[1] = 0.0f;
        view->sunDirection[2] = 0.93f;
    }
    const sim::Vec3 portEye{-22.0, 120.0, 26.0};
    const sim::Vec3 starboardEye{-22.0, -120.0, 26.0};
    const sim::Mat4 portMvp = camera(portEye, {0, 0, 3}, 130.0, portView);
    const sim::Mat4 starboardMvp = camera(starboardEye, {0, 0, 3}, 130.0, starboardView);

    const float clear[4] = {0.03f, 0.05f, 0.08f, 1.0f};
    Image port, starboard;
    expectTrue("the port view renders",
               renderThroughPng(renderer, portMvp, portView, mesh, library, clear,
                                "hull_port.png", port));
    expectTrue("the starboard view renders",
               renderThroughPng(renderer, starboardMvp, starboardView, mesh, library, clear,
                                "hull_starboard.png", starboard));
    if (!port.valid() || !starboard.valid() || port.width != kWidth) return;

    std::size_t drawn = 0, mismatched = 0, selfMismatched = 0;
    for (std::uint32_t y = 0; y < kHeight; ++y)
        for (std::uint32_t x = 0; x < kWidth; ++x) {
            const std::uint8_t* a = port.pixel(x, y);
            const std::uint8_t* mirrored = starboard.pixel(kWidth - 1 - x, y);
            const std::uint8_t* self = port.pixel(kWidth - 1 - x, y);
            if (!sameColour(a, port.pixel(0, 0), 3)) ++drawn;
            if (!sameColour(a, mirrored, 6)) ++mismatched;
            if (!sameColour(a, self, 6)) ++selfMismatched;
        }
    const double total = static_cast<double>(kWidth) * kHeight;
    std::printf("     mirror: %.1f%% of the frame is drawn, %.2f%% differs from the mirrored"
                " starboard view, %.1f%% differs from its own mirror\n",
                100.0 * static_cast<double>(drawn) / total,
                100.0 * static_cast<double>(mismatched) / total,
                100.0 * static_cast<double>(selfMismatched) / total);

    // Guard one: two blank frames mirror perfectly. CLAUDE.md records a mirror
    // check that would have passed on exactly that.
    expectTrue(label("the mirror comparison had something to compare, fraction drawn",
                     static_cast<double>(drawn) / total),
               static_cast<double>(drawn) > 0.25 * total);
    // Guard two, which the ferry's version does not have: a frame that is its own
    // mirror satisfies the check without the starboard view being involved at all.
    // The test hull is deliberately not fore-and-aft symmetric so that this cannot
    // happen quietly.
    expectTrue(label("and the frame is not its own mirror, fraction differing",
                     static_cast<double>(selfMismatched) / total),
               static_cast<double>(selfMismatched) > 0.05 * total);
    expectTrue(label("port and starboard views of a symmetric hull mirror, fraction differing",
                     static_cast<double>(mismatched) / total),
               static_cast<double>(mismatched) < 0.01 * total);
}

// --- Determinism, and the frame a person would look at ------------------------

void testRenderingIsRepeatableAndTheSceneComposes() {
    gpu::Device device;
    gpu::HullRenderer renderer;
    if (!setup(device, renderer)) return;

    gpu::MaterialLibrary library;
    std::string error;
    if (!testLibrary(library, error)) return;
    const auto seaMaterial = static_cast<std::uint32_t>(library.find("sea_water"));

    // A real seaway, and the ship floating in the same wave field the sea is drawn
    // from -- there is no second spectrum anywhere in this frame.
    sim::SeaState state;
    state.significantHeight = 2.6;
    state.peakPeriod = 8.0;
    state.frequencyCount = 24;
    state.directionCount = 6;
    state.seed = 0x5EA1;
    const sim::WaveField field(state);
    const double time = 12.0;

    sim::Sea sea;
    sea.level = 0.0;
    sea.waves = &field;
    sea.time = time;

    sim::Ship ship = makeTestShip();
    ship.initialise(sea);
    for (int step = 0; step < 60; ++step) ship.step(0.05, sea);

    gpu::OceanSurface surface;
    // Eight cells across the *shortest* component, which is the criterion
    // docs/03-renderer-audio.md arrived at -- resolving the dominant wavelength is
    // not enough, and a sea that invents a quarter of its own wave height under a
    // ship is not a frame worth looking at.
    const double shortest = gpu::shortestWavelength(field);
    const double halfExtent = 150.0;
    const int resolution = gpu::oceanResolutionFor(halfExtent, shortest, 8.0);
    surface.build(field, gpu::OceanGrid{0.0, 0.0, halfExtent, resolution}, sea.time);
    std::printf("     seaway: Hs %.1f m, Tp %.1f s, shortest component %.1f m -> %d cells"
                " over %.0f m, %zu vertices in %.1f ms\n",
                state.significantHeight, state.peakPeriod, shortest, resolution,
                2.0 * halfExtent, surface.vertices().size(), surface.buildSeconds() * 1e3);

    gpu::HullPaint paint;
    // Band edges put on the hull's own waterlines, at 3.6 m and 6.5 m. Assignment
    // is per triangle, so a boundary that falls between two waterlines is only as
    // straight as the tessellation and comes out as a sawtooth at exactly the line
    // a person looks at hardest. Fairing the bands onto the offsets table is how a
    // real hull is painted anyway; the general fix is the texture-space paint layer
    // docs/03-renderer-audio.md already plans.
    paint.waterlineZ = 5.0;
    paint.bootTopDepth = 1.4;
    paint.bootTopHeight = 1.5;
    paint.deckZ = 10.5;
    gpu::HullShading shading;
    shading.normals = gpu::HullNormals::Smooth;

    gpu::SceneMesh mesh;
    if (!mesh.appendShip(ship, paint, library, shading, error)) {
        expectTrue("the ship builds: " + error, false);
        return;
    }
    const std::size_t hullVertices = mesh.vertices().size();
    const std::size_t hullIndices = mesh.indices().size();
    mesh.appendOcean(surface, seaMaterial);

    // Two structural facts about the append that no rendered frame here can see,
    // and both were found by mutating the code until the suite stopped noticing.
    //
    // First, the sea's indices must be rebased onto the vertices it just appended.
    // Un-rebased they still address *valid* vertices -- the hull's -- so they draw
    // a plausible surface, and against a flat calm sea they draw the identical
    // surface, because shuffling the corners of a plane gives back the plane.
    std::uint32_t lowestSeaIndex = 0xffffffffu, highestSeaIndex = 0;
    for (std::size_t i = hullIndices; i < mesh.indices().size(); ++i) {
        lowestSeaIndex = std::min(lowestSeaIndex, mesh.indices()[i]);
        highestSeaIndex = std::max(highestSeaIndex, mesh.indices()[i]);
    }
    expectEqual("the sea's triangles address only the vertices the sea appended",
                static_cast<long long>(lowestSeaIndex), static_cast<long long>(hullVertices));
    expectTrue("and none of them addresses past the end",
               highestSeaIndex < mesh.vertices().size());

    // Second, the sea keeps the spectrum's own analytic slope. Re-deriving it from
    // the triangles would replace an exact answer with a finite difference of a
    // mesh that is deliberately too coarse for the shortest components, and the
    // picture would still look like the sea.
    bool normalsCarriedOver = true, positionsCarriedOver = true;
    double normalSpread = 0.0;
    for (std::size_t i = 0; i < surface.vertices().size(); ++i) {
        const gpu::OceanVertex& source = surface.vertices()[i];
        const gpu::HullVertex& copied = mesh.vertices()[hullVertices + i];
        for (int c = 0; c < 3; ++c) {
            normalsCarriedOver = normalsCarriedOver && copied.normal[c] == source.normal[c];
            positionsCarriedOver = positionsCarriedOver && copied.position[c] == source.position[c];
        }
        normalSpread = std::max(normalSpread, std::hypot(static_cast<double>(source.normal[0]),
                                                         static_cast<double>(source.normal[1])));
    }
    expectTrue("the sea's analytic normals reach the vertex buffer unchanged", normalsCarriedOver);
    expectTrue("and so do its positions", positionsCarriedOver);
    // Guard: on a calm sea every normal is +z and the check above would hold for a
    // renderer that simply wrote +z everywhere.
    expectTrue(label("this sea's normals actually tilt, worst horizontal component",
                     normalSpread),
               normalSpread > 0.05);

    gpu::SceneView view;
    view.sunDirection[0] = 0.44f;
    view.sunDirection[1] = 0.28f;
    view.sunDirection[2] = 0.85f;
    // Low and close: a sea seen from high up is a texture, and the point of the
    // frame is that the ship is *in* it.
    const sim::Vec3 eye{72.0, 82.0, 19.0};
    const sim::Mat4 mvp = camera(eye, {-4, 0, 3.0}, 110.0, view);
    const float clear[4] = {0.55f, 0.68f, 0.84f, 1.0f};

    Image first, second;
    expectTrue("the ship renders in the seaway",
               renderThroughPng(renderer, mvp, view, mesh, library, clear, "hull_ship_in_sea.png",
                                first));
    expectTrue("and again",
               renderThroughPng(renderer, mvp, view, mesh, library, clear, "hull_repeat.png",
                                second));
    expectTrue("rendering the same scene twice is byte-identical", first.rgba == second.rgba);

    // Guard: two blank frames are byte-identical too.
    std::size_t drawn = 0;
    for (std::uint32_t y = 0; y < first.height; ++y)
        for (std::uint32_t x = 0; x < first.width; ++x)
            if (!sameColour(first.pixel(x, y), first.pixel(0, 0), 3)) ++drawn;
    expectTrue(label("the repeated frames had a scene in them, fraction drawn",
                     static_cast<double>(drawn) / (double{kWidth} * kHeight)),
               static_cast<double>(drawn) > 0.25 * kWidth * kHeight);

    // Both are present, which is what "in the sea" means at the level of a frame.
    gpu::SceneView identifierView = view;
    identifierView.mode = gpu::HullShadingMode::MaterialId;
    Image identifiers;
    if (renderThroughPng(renderer, mvp, identifierView, mesh, library, clear, "hull_scene_id.png",
                         identifiers)) {
        std::size_t seaPixels = 0, hullPixels = 0;
        for (std::uint32_t y = 0; y < identifiers.height; ++y)
            for (std::uint32_t x = 0; x < identifiers.width; ++x) {
                std::uint32_t material = 0;
                if (!gpu::decodeMaterialId(identifiers.pixel(x, y), material)) continue;
                if (material == seaMaterial) ++seaPixels;
                else ++hullPixels;
            }
        std::printf("     seaway frame: %zu sea px, %zu hull px, ship at z %.2f m,"
                    " heel %.2f deg\n",
                    seaPixels, hullPixels, ship.state.position.z,
                    ship.diagnostics(sea).heelDeg);
        expectTrue("the seaway frame has both the sea and the ship in it",
                   seaPixels > 40000 && hullPixels > 8000);
    }

    // An empty scene afterwards must clear completely: no geometry from the
    // previous frame may survive in the buffers.
    gpu::SceneMesh empty;
    Image cleared;
    expectTrue("an empty scene renders",
               renderThroughPng(renderer, mvp, view, empty, library, clear, "hull_empty.png",
                                cleared));
    std::size_t background = 0;
    for (std::uint32_t y = 0; y < cleared.height; ++y)
        for (std::uint32_t x = 0; x < cleared.width; ++x)
            if (sameColour(cleared.pixel(x, y), cleared.pixel(0, 0), 0)) ++background;
    expectEqual("the previous frame does not survive the clear",
                static_cast<long long>(background), static_cast<long long>(kWidth) * kHeight);

    // Structural, and the reason it is asserted rather than assumed: the hull is
    // three vertices per triangle by construction, and the sea is one per grid
    // point. A regression to per-corner evaluation on the sea side would multiply
    // its vertex count by six for exactly the same picture -- the redundancy
    // CLAUDE.md records the physics tick paying.
    expectEqual("the hull is three vertices per triangle", static_cast<long long>(hullVertices),
                static_cast<long long>(ship.hull.tris.size()) * 3);
    expectEqual("the sea contributes one vertex per grid point",
                static_cast<long long>(mesh.vertices().size() - hullVertices),
                static_cast<long long>(surface.vertices().size()));
    expectEqual("and six indices per sea cell",
                static_cast<long long>(mesh.indices().size() - hullIndices),
                static_cast<long long>(surface.indices().size()));

    std::printf("     representative frame: %shull_ship_in_sea.png\n",
                testing::scratchDir().c_str());
}

// --- Cost ----------------------------------------------------------------------

// Measured on the shipped path, not extrapolated from a model of it. The numbers
// land in docs/03-renderer-audio.md.
void testFrameCost() {
    gpu::Device device;
    gpu::HullRenderer renderer;
    if (!setup(device, renderer)) return;

    gpu::MaterialLibrary library;
    std::string error;
    if (!testLibrary(library, error)) return;
    const auto seaMaterial = static_cast<std::uint32_t>(library.find("sea_water"));

    sim::SeaState state;
    state.significantHeight = 3.0;
    state.peakPeriod = 9.0;
    state.frequencyCount = 24;
    state.directionCount = 6;
    const sim::WaveField field(state);
    sim::Ship ship = makeTestShip();
    ship.initialise(sim::Sea{0.0});

    gpu::HullPaint paint;
    paint.waterlineZ = 5.0;
    gpu::HullShading shading;

    gpu::SceneView view;
    const sim::Vec3 eye{80.0, 95.0, 30.0};
    const sim::Mat4 mvp = camera(eye, {0, 0, 3}, 130.0, view);
    const float clear[4] = {0.55f, 0.68f, 0.84f, 1.0f};

    std::printf("     per-frame cost, %ux%u, hull %zu triangles:\n", kWidth, kHeight,
                ship.hull.tris.size());
    for (int cells : {64, 128, 192}) {
        gpu::OceanSurface surface;
        surface.build(field, gpu::OceanGrid{0.0, 0.0, 200.0, cells}, 3.0);

        gpu::SceneMesh mesh;
        double buildBest = 1e30, totalBest = 1e30;
        gpu::FrameCost best;
        for (int repeat = 0; repeat < 3; ++repeat) {
            mesh.clear();
            if (!mesh.appendShip(ship, paint, library, shading, error)) return;
            mesh.appendOcean(surface, seaMaterial);
            buildBest = std::min(buildBest, mesh.buildSeconds());

            Image image;
            float matrix[16];
            mvp.toFloats(matrix);
            if (!renderer.render(matrix, view, mesh, library, clear, image)) return;
            if (renderer.lastFrame().totalSeconds < totalBest) {
                totalBest = renderer.lastFrame().totalSeconds;
                best = renderer.lastFrame();
            }
        }
        std::printf("       sea %3d cells  %7zu vertices %8zu tris | build %6.2f ms"
                    "  upload %5.2f ms  gpu %5.2f ms  submit %5.2f ms  readback %5.2f ms"
                    "  total %6.2f ms\n",
                    cells, best.vertices, best.triangles, buildBest * 1e3,
                    best.uploadSeconds * 1e3, best.gpuSeconds * 1e3, best.submitSeconds * 1e3,
                    best.readbackSeconds * 1e3, best.totalSeconds * 1e3);

        expectEqual(label("the frame drew every triangle it was given, at cells", cells),
                    static_cast<long long>(best.triangles),
                    static_cast<long long>(mesh.indices().size() / 3));
        expectTrue(label("the GPU time was measured, ms at cells", cells),
                   best.gpuSeconds > 0.0 && best.gpuSeconds < best.totalSeconds);
    }

    // And at a resolution anyone would actually run at. This pass is fill bound
    // rather than vertex bound -- nothing is culled, so every hull triangle is
    // shaded twice -- so the small viewport above says almost nothing about the
    // number that matters.
    gpu::HullRenderer large;
    if (!large.create(device, 1920, 1080, SHIPSIM_SHADER_DIR, error)) {
        std::printf("     1920x1080 target unavailable (%s)\n", error.c_str());
        return;
    }
    gpu::OceanSurface surface;
    surface.build(field, gpu::OceanGrid{0.0, 0.0, 200.0, 192}, 3.0);
    gpu::SceneMesh mesh;
    if (!mesh.appendShip(ship, paint, library, shading, error)) return;
    mesh.appendOcean(surface, seaMaterial);

    gpu::FrameCost best;
    double bestGpu = 1e30;
    for (int repeat = 0; repeat < 4; ++repeat) {
        Image image;
        float matrix[16];
        mvp.toFloats(matrix);
        if (!large.render(matrix, view, mesh, library, clear, image)) return;
        if (large.lastFrame().gpuSeconds < bestGpu) {
            bestGpu = large.lastFrame().gpuSeconds;
            best = large.lastFrame();
        }
    }
    std::printf("       1920x1080     %7zu vertices %8zu tris |"
                "  upload %5.2f ms  gpu %5.2f ms  submit %5.2f ms  readback %5.2f ms"
                "  total %6.2f ms\n",
                best.vertices, best.triangles, best.uploadSeconds * 1e3, best.gpuSeconds * 1e3,
                best.submitSeconds * 1e3, best.readbackSeconds * 1e3, best.totalSeconds * 1e3);
    expectTrue(label("the 1080p GPU time was measured, ms", best.gpuSeconds * 1e3),
               best.gpuSeconds > 0.0 && best.gpuSeconds < best.totalSeconds);
    // Not a performance bar -- it is the statement that this pass is nowhere near
    // being the frame's problem, which is what decides whether any of it needs
    // optimising before the rest of Phase 3 lands on top of it.
    expectTrue(label("a whole ship and sea at 1080p costs well under a 60 Hz frame, ms",
                     best.gpuSeconds * 1e3),
               best.gpuSeconds < 0.016);
}

}  // namespace

void runHullRenderTests() {
    std::printf("\n--- hull rendering and materials ---\n");
    testMaterialSetIsData();
    testMaterialLoaderFailsClosed();
    testSpecularBrdfConservesEnergy();
    testFlatNormalsAreTheGeometry();
    testSmoothNormalsRespectCreases();
    testFaceAtTheSunIsFullyLitAndAtNinetyDegreesIsAmbientOnly();
    testLambertCosineLawOverASweep();
    testTwoMaterialsDifferPredictably();
    testTheFarSideOfASurfaceIsLitAsTheFaceItPresents();
    testGrazingViewFollowsSchlicksFresnel();
    testTheNearerSurfaceWinsWhicheverIsDrawnSecond();
    testShipOccludesTheSeaBehindItAndNotInFrontOfIt();
    testPaintBandsLandWhereTheSchemeSaysTheyDo();
    testPortAndStarboardViewsMirror();
    testRenderingIsRepeatableAndTheSceneComposes();
    testFrameCost();
}
