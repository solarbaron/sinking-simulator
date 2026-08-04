// SPDX-License-Identifier: MIT
//
// The lit solid: a `sim::Ship`'s hull, shaded by the material model, in the same
// frame and the same depth buffer as the sea.
//
// **The hull and the sea go through one pipeline.** That is the whole design.
// `SceneMesh` is the geometry container for everything the lit pass draws, and
// `appendOcean` converts an `OceanSurface` into it with a material index of its
// own -- so the ship is *in* the water rather than composited next to it, and the
// occlusion between them is the depth test doing its job rather than an ordering
// convention. It also means the sea is a material like any other, which is what
// makes "add a material without recompiling" true of the water too.
//
// This supersedes the flat-shaded path `tools/ferry_view` drives through
// `OffscreenRenderer`: `MeshVertex` carries a baked-in colour and no normal, so
// there is nothing there for a light to answer to. That tool is left alone; this
// is what replaces it.
//
// **Normals come from the mesh.** Flat (per face, the geometric normal) or smooth
// (area-weighted average over the faces meeting at a welded position, with a
// crease angle so the deck edge stays an edge). Flat is the default because it is
// what the geometry actually says and because a flat face is then exactly one
// colour, which is a closed form a test can state. Smooth exists because a hull
// tabulated at twenty-four stations is a smooth surface that has been sampled,
// not a faceted one.
//
// **Paint is data.** `HullPaint` names materials and the z bands they occupy --
// antifouling, boot-topping at the waterline, topside above it, and an optional
// deck material for upward-facing faces. Names are resolved against a
// `MaterialLibrary` at build time and a name that is not there is an error, not a
// silent default. A mod repaints a ship by naming different materials and adds a
// new one by loading another `.materials` file; see `material.hpp`.
//
// Pascal target: vertex/fragment over indexed triangles, one storage buffer for
// the material table. No mesh shaders, no ray tracing.
#pragma once

#include "../core/geometry.hpp"
#include "../core/math.hpp"
#include "../core/png.hpp"
#include "../sim/ship.hpp"
#include "damage.hpp"
#include "device.hpp"
#include "material.hpp"
#include "ocean.hpp"

#include <cstdint>
#include <string>
#include <vector>

struct VkRenderPass_T;
struct VkFramebuffer_T;
struct VkPipeline_T;
struct VkPipelineLayout_T;
struct VkDescriptorSetLayout_T;
struct VkDescriptorPool_T;
struct VkDescriptorSet_T;
struct VkQueryPool_T;

namespace gpu {

// World position, world normal, and the row of the material table to shade with.
// 28 bytes, every field naturally aligned, no padding.
struct HullVertex {
    float position[3];
    float normal[3];
    std::uint32_t material;
};
static_assert(sizeof(HullVertex) == 28, "the vertex must match shaders/hull.vert");

enum class HullNormals {
    // The geometric normal of each triangle, on all three of its corners. A flat
    // face is then exactly one colour.
    Flat,
    // Area-weighted average of the faces meeting at a welded position, except
    // across an edge sharper than the crease angle, where the face keeps its own
    // normal. Area weighting rather than plain averaging because a fan of slivers
    // must not outvote the one large face they sit against.
    Smooth,
};

struct HullShading {
    HullNormals normals = HullNormals::Flat;
    // Faces meeting at more than this angle keep their own normals. A ship has
    // both -- a smooth turn of the bilge and a hard deck edge -- so a single
    // choice of flat or smooth is wrong for one of them.
    double creaseAngle = 40.0 * sim::kDegToRad;
    // Positions closer than this are treated as the same point when smoothing.
    // Matches the default weld epsilon of `sim::isClosedManifold`, so a mesh that
    // passes the manifold check welds the same way here.
    double weldEpsilon = 1e-6;
};

// A paint scheme, in the ship's own body frame. Bands are half-open in z:
//
//     z <  waterlineZ - bootTopDepth        underwater
//     ...  within the boot-topping band     bootTopping
//     z >= waterlineZ + bootTopHeight       topside
//
// and any face above `deckZ` whose normal points up by at least `deckNormalZ`
// takes `deck` instead. Assignment is per triangle, from its centroid and its
// geometric normal, so it does not depend on how the hull happens to be
// tessellated.
struct HullPaint {
    std::string underwater = "antifouling";
    std::string bootTopping = "boot_topping";
    std::string topside = "painted_steel_topside";
    std::string deck = "timber_deck";
    // The plating around a hole: steel that was inside the plate a moment ago and
    // has no paint on it. A name like every other band, so a mod restates it or
    // replaces it in a `.materials` file rather than in a shader -- and it is
    // resolved **only on the damaged path**, so a library that has never heard of
    // it still paints an intact ship.
    std::string tornEdge = "torn_plate_edge";

    double waterlineZ = 5.5;      // m, body frame: the design waterline
    double bootTopDepth = 0.35;   // m below it
    double bootTopHeight = 0.85;  // m above it

    double deckZ = 1e30;          // m; the default puts the deck material nowhere
    double deckNormalZ = 0.85;    // how upward-facing an upward-facing face is
};

// Geometry for the lit pass. Named for the scene rather than for the hull because
// the sea goes through it too, which is the point.
class SceneMesh {
public:
    void clear();

    // The ship's hull, transformed into world space by its own rigid-body state
    // and painted by `paint`. False with a reason when a material name is not in
    // the library -- a missing material is a broken ship definition, not a grey
    // default.
    bool appendShip(const sim::Ship& ship, const HullPaint& paint,
                    const MaterialLibrary& library, const HullShading& shading,
                    std::string& error);

    // The same ship after `damage.hpp` has refined, displaced and cut her plating.
    // The geometry drawn is `damaged.deformed`; the paint bands are decided on
    // `damaged.rest`, because paint is on the hull and a dented plate keeps the
    // paint it was painted with -- deciding a band from a displaced z would slide
    // the boot-topping around as the plating moved, which is the world-frame
    // mistake one level down.
    //
    // Costs nothing when there is nothing wrong with her: a `DamagedHull` built
    // from an empty `HullDamage` is the ship's own mesh, and this produces byte
    // for byte what the call above does.
    bool appendShip(const sim::Ship& ship, const DamagedHull& damaged, const HullPaint& paint,
                    const MaterialLibrary& library, const HullShading& shading,
                    std::string& error);

    // Any mesh, one material, transformed by R * v + t.
    void appendMesh(const sim::TriMesh& mesh, const sim::Mat3& rotation,
                    const sim::Vec3& translation, std::uint32_t material,
                    const HullShading& shading);

    // The sea. Positions and normals come straight across -- `OceanSurface`
    // already carries the spectrum's own analytic slope, and re-deriving a normal
    // from the mesh here would throw that away.
    void appendOcean(const OceanSurface& surface, std::uint32_t material);
    // The cascade, which is the sea that reaches the horizon. Same conversion:
    // its levels already share their seam vertices, so there is nothing to weld
    // and nothing here that could unweld them.
    void appendOcean(const OceanCascadeSurface& surface, std::uint32_t material);

    const std::vector<HullVertex>& vertices() const { return vertices_; }
    const std::vector<std::uint32_t>& indices() const { return indices_; }
    // Wall-clock seconds spent in the append calls since the last clear().
    // Reported rather than asserted; docs/03-renderer-audio.md carries the
    // numbers.
    double buildSeconds() const { return buildSeconds_; }

private:
    // The one place normals and materials are decided. `perFaceMaterial`, when
    // present, is one index per triangle and overrides `singleMaterial`.
    void appendTriangles(const sim::TriMesh& mesh, const sim::Mat3& rotation,
                         const sim::Vec3& translation, const HullShading& shading,
                         const std::uint32_t* perFaceMaterial, std::uint32_t singleMaterial);
    void appendOceanGeometry(const std::vector<OceanVertex>& vertices,
                             const std::vector<std::uint32_t>& indices, std::uint32_t material);

    std::vector<HullVertex> vertices_;
    std::vector<std::uint32_t> indices_;
    double buildSeconds_ = 0.0;
    // Scratch, kept across calls so a per-frame rebuild does not reallocate.
    std::vector<sim::Vec3> world_, faceNormal_, vertexNormal_;
    std::vector<std::uint32_t> weld_;
    std::vector<std::uint32_t> perFaceMaterial_;
};

enum class HullShadingMode {
    // The metallic-roughness BRDF of material.hpp.
    Shaded,
    // The material index of whatever won the depth test, as a 16-bit code: high
    // byte in red, low byte in green, blue 255 as a tag. A measurement channel,
    // not a debug ramp -- it makes occlusion an exact integer question.
    MaterialId,
    // Clip-space depth as a 16-bit code, same encoding. Exact to 1 / 65535 of the
    // [0, 1] range, which is what lets a depth assertion be a number rather than
    // a colour comparison.
    Depth,
};

struct SceneView {
    HullShadingMode mode = HullShadingMode::Shaded;
    // Unit vector *toward* the sun, world frame. Normalised on submission, so a
    // caller who hands over a direction of length 1.001 does not quietly get
    // lighting that is 0.1% too bright.
    float sunDirection[3] = {0.42f, 0.30f, 0.86f};
    // Radiance, not a 0-1 colour: the diffuse lobe carries a 1 / pi, so a sun
    // bright enough to light a 0.7-albedo hull to 0.7 is about pi.
    float sunColour[3] = {3.4f, 3.3f, 3.1f};
    // Hemispheric sky irradiance, multiplied by (0.5 + 0.5 n.z).
    float skyColour[3] = {0.26f, 0.31f, 0.38f};
    float exposure = 1.0f;
    // World-space camera position. **Must be the eye the `mvp` was built from** or
    // the specular lobe answers to the wrong direction while everything else looks
    // right. `eyeAgreesWithCamera` checks it; the tests assert that check rather
    // than trusting the pairing.
    float eye[3] = {0.0f, 0.0f, 0.0f};
};

// True when `eye` is the centre of projection of `mvp`. For a perspective matrix
// the eye is exactly the point the projection sends to w = 0, so this is an
// identity rather than a heuristic; `scale` is the world size of the scene, used
// only to turn the tolerance into a relative one. Always false for an orthographic
// projection, which has no centre of projection.
bool eyeAgreesWithCamera(const float mvp[16], const float eye[3], double scale);

// Inverse of the MaterialId encoding. False when the pixel is background. Sixteen
// bits, so a library past 65535 surfaces renders correctly but cannot be read back
// through this channel -- the table is never truncated to suit it.
bool decodeMaterialId(const std::uint8_t* pixel, std::uint32_t& material);
// Inverse of the Depth encoding. False when the pixel is background.
bool decodeSceneDepth(const std::uint8_t* pixel, double& depth);

// Where a frame went. Measured with a GPU timestamp pair around the render pass
// and a steady clock around everything else, because the interesting question is
// not the total but which part of it is the total.
struct FrameCost {
    double uploadSeconds = 0;    // vertex, index and material staging
    double gpuSeconds = 0;       // timestamps around the render pass itself
    double submitSeconds = 0;    // command buffer build, submit and queue wait
    double readbackSeconds = 0;  // colour attachment to core::Image
    double totalSeconds = 0;
    std::size_t vertices = 0;
    std::size_t triangles = 0;
};

// Offscreen colour + depth target and the lit pipeline. A sibling of
// OceanRenderer rather than a user of it: this one needs a material index on the
// vertex, a storage buffer to look it up in, and an eye position for the
// view-dependent term, and none of those exist there.
class HullRenderer {
public:
    ~HullRenderer();

    bool create(Device& device, std::uint32_t width, std::uint32_t height,
                const std::string& shaderDirectory, std::string& error);
    void destroy();
    bool valid() const { return pipeline_ != nullptr; }

    std::uint32_t width() const { return width_; }
    std::uint32_t height() const { return height_; }

    // Clears to `clearColour`, draws `mesh` under `mvp` with depth testing, and
    // reads the colour attachment back into `out`. Triangles are drawn in index
    // order, so appending the sea after the hull makes the depth test carry the
    // occlusion rather than the draw order -- which is the arrangement the tests
    // use deliberately.
    //
    // An empty mesh is legal and yields the clear colour.
    bool render(const float mvp[16], const SceneView& view, const SceneMesh& mesh,
                const MaterialLibrary& library, const float clearColour[4], core::Image& out);

    const FrameCost& lastFrame() const { return lastFrame_; }

private:
    bool ensureGeometryCapacity(std::size_t vertexBytes, std::size_t indexBytes);
    bool ensureMaterials(const MaterialLibrary& library);

    Device*             device_ = nullptr;
    std::uint32_t       width_ = 0;
    std::uint32_t       height_ = 0;

    Device::Image2D     colour_;
    Device::Image2D     depth_;
    VkRenderPass_T*     renderPass_ = nullptr;
    VkFramebuffer_T*    framebuffer_ = nullptr;
    VkDescriptorSetLayout_T* setLayout_ = nullptr;
    VkDescriptorPool_T* descriptorPool_ = nullptr;
    VkDescriptorSet_T*  descriptorSet_ = nullptr;
    VkPipelineLayout_T* pipelineLayout_ = nullptr;
    VkPipeline_T*       pipeline_ = nullptr;
    VkQueryPool_T*      queryPool_ = nullptr;

    Device::Buffer      vertexBuffer_;
    Device::Buffer      indexBuffer_;
    Device::Buffer      materialBuffer_;
    Device::Buffer      readback_;

    std::uint64_t       materialRevision_ = 0;
    std::size_t         materialCapacity_ = 0;
    FrameCost           lastFrame_;
};

}  // namespace gpu
