// SPDX-License-Identifier: MIT
#include "hull.hpp"

#include <vulkan/vulkan.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstring>
#include <map>

namespace gpu {
namespace {

constexpr VkFormat kColourFormat = VK_FORMAT_R8G8B8A8_UNORM;
constexpr VkFormat kDepthFormat = VK_FORMAT_D32_SFLOAT;

// Mirrors the `Push` block in shaders/hull.vert and shaders/hull.frag. 128 bytes
// is exactly the push-constant size every Vulkan implementation guarantees, so
// this fits everywhere and has no room to spare -- anything further belongs in a
// uniform buffer rather than in here.
struct Push {
    float modelViewProjection[16];
    float sun[4];        // xyz unit vector toward the sun, w mode
    float sunColour[4];  // rgb radiance, w exposure
    float sky[4];        // rgb hemispheric ambient, w unused
    float eye[4];        // xyz world eye position, w unused
};
static_assert(sizeof(Push) == 128, "the push block must match the shaders");

std::array<long long, 3> weldKey(const sim::Vec3& p, double epsilon) {
    const double inverse = 1.0 / epsilon;
    return {std::llround(p.x * inverse), std::llround(p.y * inverse),
            std::llround(p.z * inverse)};
}

// The four bands, resolved once. A ship whose paint names a material nobody loaded
// is a broken ship definition; rendering it grey would hide that until someone
// looked at it, which is the failure mode docs/05 is written against.
bool resolveBands(const HullPaint& paint, const MaterialLibrary& library,
                  std::uint32_t resolved[4], std::string& error) {
    struct Named { const std::string& name; const char* role; };
    const Named wanted[4] = {{paint.underwater, "underwater"},
                             {paint.bootTopping, "boot topping"},
                             {paint.topside, "topside"},
                             {paint.deck, "deck"}};
    for (int i = 0; i < 4; ++i) {
        const int index = library.find(wanted[i].name);
        if (index < 0) {
            error = std::string("the ") + wanted[i].role + " material '" + wanted[i].name +
                    "' is not in the library";
            return false;
        }
        resolved[i] = static_cast<std::uint32_t>(index);
    }
    return true;
}

// Which band one triangle of the hull falls in, from its centroid and its
// geometric normal. **Decided in the body frame, not the world frame** -- paint is
// on the hull, so a ship heeled thirty degrees still has her boot-topping where it
// was painted, and deciding the bands from world z would slide the waterline
// stripe around the hull as she rolls.
std::uint32_t bandFor(const sim::Vec3& a, const sim::Vec3& b, const sim::Vec3& c,
                      const HullPaint& paint, const std::uint32_t resolved[4]) {
    const sim::Vec3 centroid = (a + b + c) / 3.0;
    const sim::Vec3 normal = normalize(cross(b - a, c - a));
    if (centroid.z >= paint.deckZ && normal.z >= paint.deckNormalZ) return resolved[3];
    if (centroid.z < paint.waterlineZ - paint.bootTopDepth) return resolved[0];
    if (centroid.z < paint.waterlineZ + paint.bootTopHeight) return resolved[1];
    return resolved[2];
}

}  // namespace

// --- Geometry -----------------------------------------------------------------

void SceneMesh::clear() {
    vertices_.clear();
    indices_.clear();
    buildSeconds_ = 0.0;
}

// Shared by the single-material and the painted paths. `perFaceMaterial`, when
// present, is one index per triangle and overrides `singleMaterial`.
void SceneMesh::appendTriangles(const sim::TriMesh& mesh, const sim::Mat3& rotation,
                                const sim::Vec3& translation, const HullShading& shading,
                                const std::uint32_t* perFaceMaterial,
                                std::uint32_t singleMaterial) {
    const std::size_t triangleCount = mesh.tris.size();
    if (triangleCount == 0) return;

    world_.resize(mesh.verts.size());
    for (std::size_t i = 0; i < mesh.verts.size(); ++i)
        world_[i] = rotation * mesh.verts[i] + translation;

    // Unnormalised, so its length is twice the triangle's area -- which is exactly
    // the weight the smooth path wants. A fan of slivers must not outvote the one
    // large face they sit against.
    faceNormal_.resize(triangleCount);
    for (std::size_t t = 0; t < triangleCount; ++t) {
        const sim::Tri& tri = mesh.tris[t];
        faceNormal_[t] = cross(world_[tri.b] - world_[tri.a], world_[tri.c] - world_[tri.a]);
    }

    const bool smooth = shading.normals == HullNormals::Smooth;
    if (smooth) {
        // Weld by position: `sim::TriMesh` shares vertices by index, but a mesh
        // that has been mirrored or clipped carries coincident-but-distinct
        // vertices along the seam, and smoothing across the seam is what makes a
        // mirrored hull look like one surface rather than two.
        std::map<std::array<long long, 3>, std::uint32_t> lookup;
        weld_.resize(mesh.verts.size());
        std::uint32_t representatives = 0;
        for (std::size_t i = 0; i < mesh.verts.size(); ++i) {
            const auto inserted =
                lookup.emplace(weldKey(world_[i], shading.weldEpsilon), representatives);
            weld_[i] = inserted.first->second;
            if (inserted.second) ++representatives;
        }
        vertexNormal_.assign(representatives, sim::Vec3{});
        for (std::size_t t = 0; t < triangleCount; ++t) {
            const sim::Tri& tri = mesh.tris[t];
            const sim::Vec3& weighted = faceNormal_[t];
            vertexNormal_[weld_[tri.a]] += weighted;
            vertexNormal_[weld_[tri.b]] += weighted;
            vertexNormal_[weld_[tri.c]] += weighted;
        }
    }

    const double creaseCos = std::cos(std::clamp(shading.creaseAngle, 0.0, sim::kPi));
    const auto base = static_cast<std::uint32_t>(vertices_.size());
    vertices_.reserve(vertices_.size() + triangleCount * 3);
    indices_.reserve(indices_.size() + triangleCount * 3);

    for (std::size_t t = 0; t < triangleCount; ++t) {
        const sim::Tri& tri = mesh.tris[t];
        const std::uint32_t material =
            perFaceMaterial != nullptr ? perFaceMaterial[t] : singleMaterial;
        const sim::Vec3 face = normalize(faceNormal_[t]);
        const std::uint32_t corner[3] = {tri.a, tri.b, tri.c};

        for (int k = 0; k < 3; ++k) {
            sim::Vec3 n = face;
            if (smooth) {
                const sim::Vec3 averaged = normalize(vertexNormal_[weld_[corner[k]]]);
                // A face that disagrees with the average by more than the crease
                // angle keeps its own normal, which is what leaves the deck edge an
                // edge while the turn of the bilge stays smooth.
                if (dot(averaged, face) >= creaseCos && length2(averaged) > 0.5) n = averaged;
            }
            const sim::Vec3& p = world_[corner[k]];
            HullVertex vertex{};
            vertex.position[0] = static_cast<float>(p.x);
            vertex.position[1] = static_cast<float>(p.y);
            vertex.position[2] = static_cast<float>(p.z);
            vertex.normal[0] = static_cast<float>(n.x);
            vertex.normal[1] = static_cast<float>(n.y);
            vertex.normal[2] = static_cast<float>(n.z);
            vertex.material = material;
            vertices_.push_back(vertex);
        }
        // Every triangle gets its own three vertices. That is what flat shading
        // needs, and it is also what lets two triangles sharing a position carry
        // different materials -- which a painted hull requires at every band
        // boundary.
        indices_.push_back(base + static_cast<std::uint32_t>(t) * 3 + 0);
        indices_.push_back(base + static_cast<std::uint32_t>(t) * 3 + 1);
        indices_.push_back(base + static_cast<std::uint32_t>(t) * 3 + 2);
    }
}

void SceneMesh::appendMesh(const sim::TriMesh& mesh, const sim::Mat3& rotation,
                           const sim::Vec3& translation, std::uint32_t material,
                           const HullShading& shading) {
    const auto started = std::chrono::steady_clock::now();
    appendTriangles(mesh, rotation, translation, shading, nullptr, material);
    buildSeconds_ += std::chrono::duration<double>(std::chrono::steady_clock::now() - started)
                         .count();
}

bool SceneMesh::appendShip(const sim::Ship& ship, const HullPaint& paint,
                           const MaterialLibrary& library, const HullShading& shading,
                           std::string& error) {
    const auto started = std::chrono::steady_clock::now();

    std::uint32_t resolved[4]{};
    if (!resolveBands(paint, library, resolved, error)) return false;

    const sim::Mat3 rotation = ship.state.orientation.toMat3();
    const sim::Vec3 translation = ship.state.position;

    perFaceMaterial_.resize(ship.hull.tris.size());
    for (std::size_t t = 0; t < ship.hull.tris.size(); ++t) {
        const sim::Tri& tri = ship.hull.tris[t];
        perFaceMaterial_[t] = bandFor(ship.hull.verts[tri.a], ship.hull.verts[tri.b],
                                      ship.hull.verts[tri.c], paint, resolved);
    }

    appendTriangles(ship.hull, rotation, translation, shading, perFaceMaterial_.data(), 0);
    buildSeconds_ += std::chrono::duration<double>(std::chrono::steady_clock::now() - started)
                         .count();
    return true;
}

bool SceneMesh::appendShip(const sim::Ship& ship, const DamagedHull& damaged,
                           const HullPaint& paint, const MaterialLibrary& library,
                           const HullShading& shading, std::string& error) {
    const auto started = std::chrono::steady_clock::now();

    if (damaged.rest.tris.size() != damaged.deformed.tris.size() ||
        damaged.exposed.size() != damaged.deformed.tris.size()) {
        error = "the damaged hull's rest and deformed meshes do not have the same topology";
        return false;
    }

    std::uint32_t resolved[4]{};
    if (!resolveBands(paint, library, resolved, error)) return false;
    // Resolved on this path only, so a library that predates the torn-edge material
    // still paints an intact ship. Missing is an error rather than a fallback, for
    // the reason the four bands are.
    const int tornEdge = library.find(paint.tornEdge);
    if (tornEdge < 0) {
        error = "the torn-edge material '" + paint.tornEdge + "' is not in the library";
        return false;
    }

    const sim::Mat3 rotation = ship.state.orientation.toMat3();
    const sim::Vec3 translation = ship.state.position;

    // Bands from the **undeformed** hull, for the same reason they come from the
    // body frame rather than the world one: a dented plate keeps its paint. Only
    // the exposed band is a property of the damage.
    perFaceMaterial_.resize(damaged.rest.tris.size());
    for (std::size_t t = 0; t < damaged.rest.tris.size(); ++t) {
        const sim::Tri& tri = damaged.rest.tris[t];
        perFaceMaterial_[t] =
            damaged.exposed[t] != 0
                ? static_cast<std::uint32_t>(tornEdge)
                : bandFor(damaged.rest.verts[tri.a], damaged.rest.verts[tri.b],
                          damaged.rest.verts[tri.c], paint, resolved);
    }

    appendTriangles(damaged.deformed, rotation, translation, shading, perFaceMaterial_.data(), 0);
    buildSeconds_ += std::chrono::duration<double>(std::chrono::steady_clock::now() - started)
                         .count();
    return true;
}

void SceneMesh::appendOcean(const OceanSurface& surface, std::uint32_t material) {
    appendOceanGeometry(surface.vertices(), surface.indices(), material);
}

void SceneMesh::appendOcean(const OceanCascadeSurface& surface, std::uint32_t material) {
    appendOceanGeometry(surface.vertices(), surface.indices(), material);
}

void SceneMesh::appendOceanGeometry(const std::vector<OceanVertex>& source,
                                    const std::vector<std::uint32_t>& sourceIndices,
                                    std::uint32_t material) {
    const auto started = std::chrono::steady_clock::now();
    const auto base = static_cast<std::uint32_t>(vertices_.size());
    vertices_.reserve(vertices_.size() + source.size());
    for (const OceanVertex& v : source) {
        HullVertex vertex{};
        // Straight across. `OceanSurface` already carries the spectrum's own
        // analytic slope as the normal, and re-deriving one from the triangles
        // here would replace an exact answer with a finite difference of a mesh
        // that is deliberately too coarse for the shortest components.
        std::memcpy(vertex.position, v.position, sizeof vertex.position);
        std::memcpy(vertex.normal, v.normal, sizeof vertex.normal);
        vertex.material = material;
        vertices_.push_back(vertex);
    }
    indices_.reserve(indices_.size() + sourceIndices.size());
    for (std::uint32_t index : sourceIndices) indices_.push_back(base + index);

    buildSeconds_ += std::chrono::duration<double>(std::chrono::steady_clock::now() - started)
                         .count();
}

// --- The readback channels ----------------------------------------------------

bool decodeMaterialId(const std::uint8_t* pixel, std::uint32_t& material) {
    if (pixel == nullptr || pixel[2] < 128) return false;
    material = static_cast<std::uint32_t>(pixel[0]) * 256u + static_cast<std::uint32_t>(pixel[1]);
    return true;
}

bool decodeSceneDepth(const std::uint8_t* pixel, double& depth) {
    if (pixel == nullptr || pixel[2] < 128) return false;
    depth = (static_cast<double>(pixel[0]) * 256.0 + static_cast<double>(pixel[1])) / 65535.0;
    return true;
}

bool eyeAgreesWithCamera(const float mvp[16], const float eye[3], double scale) {
    if (mvp == nullptr || eye == nullptr) return false;
    // Column-major, m[column * 4 + row], matching sim::Mat4 and GLSL.
    double w = mvp[15];
    for (int c = 0; c < 3; ++c) w += static_cast<double>(mvp[c * 4 + 3]) * eye[c];
    // A perspective projection sends its own centre of projection to w = 0
    // exactly; an orthographic one has row 3 equal to (0, 0, 0, 1) and gives 1
    // for every point, so it never agrees and never claims to.
    return std::abs(w) <= 1e-5 * std::max(scale, 1.0);
}

// --- The pipeline -------------------------------------------------------------

HullRenderer::~HullRenderer() { destroy(); }

bool HullRenderer::create(Device& device, std::uint32_t width, std::uint32_t height,
                          const std::string& shaderDirectory, std::string& error) {
    destroy();
    if (!device.valid()) {
        error = "device is not valid";
        return false;
    }
    if (!device.supportsGraphics()) {
        error = "queue family does not support graphics";
        return false;
    }
    device_ = &device;
    width_ = width;
    height_ = height;

    colour_ = device.createImage2D(width, height, kColourFormat,
                                   VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT |
                                       VK_IMAGE_USAGE_TRANSFER_SRC_BIT,
                                   VK_IMAGE_ASPECT_COLOR_BIT);
    depth_ = device.createImage2D(width, height, kDepthFormat,
                                  VK_IMAGE_USAGE_DEPTH_STENCIL_ATTACHMENT_BIT,
                                  VK_IMAGE_ASPECT_DEPTH_BIT);
    if (!colour_.valid() || !depth_.valid()) {
        error = "could not allocate the render target";
        destroy();
        return false;
    }

    // Same layout chain as OceanRenderer: the pass leaves the colour attachment
    // in TRANSFER_SRC_OPTIMAL so the readback copy needs no barrier.
    VkAttachmentDescription attachments[2]{};
    attachments[0].format = kColourFormat;
    attachments[0].samples = VK_SAMPLE_COUNT_1_BIT;
    attachments[0].loadOp = VK_ATTACHMENT_LOAD_OP_CLEAR;
    attachments[0].storeOp = VK_ATTACHMENT_STORE_OP_STORE;
    attachments[0].stencilLoadOp = VK_ATTACHMENT_LOAD_OP_DONT_CARE;
    attachments[0].stencilStoreOp = VK_ATTACHMENT_STORE_OP_DONT_CARE;
    attachments[0].initialLayout = VK_IMAGE_LAYOUT_UNDEFINED;
    attachments[0].finalLayout = VK_IMAGE_LAYOUT_TRANSFER_SRC_OPTIMAL;

    attachments[1].format = kDepthFormat;
    attachments[1].samples = VK_SAMPLE_COUNT_1_BIT;
    attachments[1].loadOp = VK_ATTACHMENT_LOAD_OP_CLEAR;
    attachments[1].storeOp = VK_ATTACHMENT_STORE_OP_DONT_CARE;
    attachments[1].stencilLoadOp = VK_ATTACHMENT_LOAD_OP_DONT_CARE;
    attachments[1].stencilStoreOp = VK_ATTACHMENT_STORE_OP_DONT_CARE;
    attachments[1].initialLayout = VK_IMAGE_LAYOUT_UNDEFINED;
    attachments[1].finalLayout = VK_IMAGE_LAYOUT_DEPTH_STENCIL_ATTACHMENT_OPTIMAL;

    VkAttachmentReference colourRef{0, VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL};
    VkAttachmentReference depthRef{1, VK_IMAGE_LAYOUT_DEPTH_STENCIL_ATTACHMENT_OPTIMAL};

    VkSubpassDescription subpass{};
    subpass.pipelineBindPoint = VK_PIPELINE_BIND_POINT_GRAPHICS;
    subpass.colorAttachmentCount = 1;
    subpass.pColorAttachments = &colourRef;
    subpass.pDepthStencilAttachment = &depthRef;

    VkRenderPassCreateInfo renderPassInfo{VK_STRUCTURE_TYPE_RENDER_PASS_CREATE_INFO};
    renderPassInfo.attachmentCount = 2;
    renderPassInfo.pAttachments = attachments;
    renderPassInfo.subpassCount = 1;
    renderPassInfo.pSubpasses = &subpass;
    if (vkCreateRenderPass(device.handle(), &renderPassInfo, nullptr, &renderPass_) != VK_SUCCESS) {
        error = "vkCreateRenderPass failed";
        destroy();
        return false;
    }

    VkImageView views[2] = {colour_.view, depth_.view};
    VkFramebufferCreateInfo framebufferInfo{VK_STRUCTURE_TYPE_FRAMEBUFFER_CREATE_INFO};
    framebufferInfo.renderPass = renderPass_;
    framebufferInfo.attachmentCount = 2;
    framebufferInfo.pAttachments = views;
    framebufferInfo.width = width;
    framebufferInfo.height = height;
    framebufferInfo.layers = 1;
    if (vkCreateFramebuffer(device.handle(), &framebufferInfo, nullptr, &framebuffer_) !=
        VK_SUCCESS) {
        error = "vkCreateFramebuffer failed";
        destroy();
        return false;
    }

    // One storage buffer holding the whole material table. This is the small end
    // of the bindless arrangement docs/03-renderer-audio.md plans: the draw
    // carries indices, the shader looks them up, and adding a material is adding
    // a row rather than editing a shader.
    VkDescriptorSetLayoutBinding binding{};
    binding.binding = 0;
    binding.descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
    binding.descriptorCount = 1;
    binding.stageFlags = VK_SHADER_STAGE_FRAGMENT_BIT;
    VkDescriptorSetLayoutCreateInfo setLayoutInfo{
        VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO};
    setLayoutInfo.bindingCount = 1;
    setLayoutInfo.pBindings = &binding;
    if (vkCreateDescriptorSetLayout(device.handle(), &setLayoutInfo, nullptr, &setLayout_) !=
        VK_SUCCESS) {
        error = "vkCreateDescriptorSetLayout failed";
        destroy();
        return false;
    }

    VkDescriptorPoolSize poolSize{VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 1};
    VkDescriptorPoolCreateInfo poolInfo{VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO};
    poolInfo.maxSets = 1;
    poolInfo.poolSizeCount = 1;
    poolInfo.pPoolSizes = &poolSize;
    if (vkCreateDescriptorPool(device.handle(), &poolInfo, nullptr, &descriptorPool_) !=
        VK_SUCCESS) {
        error = "vkCreateDescriptorPool failed";
        destroy();
        return false;
    }

    VkDescriptorSetAllocateInfo setAllocation{VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO};
    setAllocation.descriptorPool = descriptorPool_;
    setAllocation.descriptorSetCount = 1;
    setAllocation.pSetLayouts = &setLayout_;
    if (vkAllocateDescriptorSets(device.handle(), &setAllocation, &descriptorSet_) != VK_SUCCESS) {
        error = "vkAllocateDescriptorSets failed";
        destroy();
        return false;
    }

    VkShaderModule vertexModule = device.loadShader(shaderDirectory + "/hull.vert.spv", error);
    if (vertexModule == nullptr) { destroy(); return false; }
    VkShaderModule fragmentModule = device.loadShader(shaderDirectory + "/hull.frag.spv", error);
    if (fragmentModule == nullptr) {
        device.destroyShader(vertexModule);
        destroy();
        return false;
    }

    VkPushConstantRange push{VK_SHADER_STAGE_VERTEX_BIT | VK_SHADER_STAGE_FRAGMENT_BIT, 0,
                             sizeof(Push)};
    VkPipelineLayoutCreateInfo layoutInfo{VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO};
    layoutInfo.setLayoutCount = 1;
    layoutInfo.pSetLayouts = &setLayout_;
    layoutInfo.pushConstantRangeCount = 1;
    layoutInfo.pPushConstantRanges = &push;
    if (vkCreatePipelineLayout(device.handle(), &layoutInfo, nullptr, &pipelineLayout_) !=
        VK_SUCCESS) {
        error = "vkCreatePipelineLayout failed";
        device.destroyShader(vertexModule);
        device.destroyShader(fragmentModule);
        destroy();
        return false;
    }

    VkPipelineShaderStageCreateInfo stages[2]{};
    stages[0].sType = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO;
    stages[0].stage = VK_SHADER_STAGE_VERTEX_BIT;
    stages[0].module = vertexModule;
    stages[0].pName = "main";
    stages[1].sType = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO;
    stages[1].stage = VK_SHADER_STAGE_FRAGMENT_BIT;
    stages[1].module = fragmentModule;
    stages[1].pName = "main";

    VkVertexInputBindingDescription vertexBinding{0, sizeof(HullVertex),
                                                  VK_VERTEX_INPUT_RATE_VERTEX};
    VkVertexInputAttributeDescription vertexAttributes[3]{};
    vertexAttributes[0] = {0, 0, VK_FORMAT_R32G32B32_SFLOAT, offsetof(HullVertex, position)};
    vertexAttributes[1] = {1, 0, VK_FORMAT_R32G32B32_SFLOAT, offsetof(HullVertex, normal)};
    vertexAttributes[2] = {2, 0, VK_FORMAT_R32_UINT, offsetof(HullVertex, material)};

    VkPipelineVertexInputStateCreateInfo vertexInput{
        VK_STRUCTURE_TYPE_PIPELINE_VERTEX_INPUT_STATE_CREATE_INFO};
    vertexInput.vertexBindingDescriptionCount = 1;
    vertexInput.pVertexBindingDescriptions = &vertexBinding;
    vertexInput.vertexAttributeDescriptionCount = 3;
    vertexInput.pVertexAttributeDescriptions = vertexAttributes;

    VkPipelineInputAssemblyStateCreateInfo assembly{
        VK_STRUCTURE_TYPE_PIPELINE_INPUT_ASSEMBLY_STATE_CREATE_INFO};
    assembly.topology = VK_PRIMITIVE_TOPOLOGY_TRIANGLE_LIST;

    VkViewport viewport{0.0f, 0.0f, static_cast<float>(width), static_cast<float>(height),
                        0.0f, 1.0f};
    VkRect2D scissor{{0, 0}, {width, height}};
    VkPipelineViewportStateCreateInfo viewportState{
        VK_STRUCTURE_TYPE_PIPELINE_VIEWPORT_STATE_CREATE_INFO};
    viewportState.viewportCount = 1;
    viewportState.pViewports = &viewport;
    viewportState.scissorCount = 1;
    viewportState.pScissors = &scissor;

    VkPipelineRasterizationStateCreateInfo raster{
        VK_STRUCTURE_TYPE_PIPELINE_RASTERIZATION_STATE_CREATE_INFO};
    raster.polygonMode = VK_POLYGON_MODE_FILL;
    // Nothing is culled, for the reason the sea is not: a hull is seen from the
    // inside once it is cut away or flooded. The fragment shader flips the normal
    // for the far side instead. On a closed hull the front faces win the depth
    // test anyway, so this costs fill rate and no correctness.
    raster.cullMode = VK_CULL_MODE_NONE;
    raster.frontFace = VK_FRONT_FACE_COUNTER_CLOCKWISE;
    raster.lineWidth = 1.0f;

    VkPipelineMultisampleStateCreateInfo multisample{
        VK_STRUCTURE_TYPE_PIPELINE_MULTISAMPLE_STATE_CREATE_INFO};
    multisample.rasterizationSamples = VK_SAMPLE_COUNT_1_BIT;

    VkPipelineDepthStencilStateCreateInfo depthStencil{
        VK_STRUCTURE_TYPE_PIPELINE_DEPTH_STENCIL_STATE_CREATE_INFO};
    depthStencil.depthTestEnable = VK_TRUE;
    depthStencil.depthWriteEnable = VK_TRUE;
    depthStencil.depthCompareOp = VK_COMPARE_OP_LESS;

    VkPipelineColorBlendAttachmentState blendAttachment{};
    blendAttachment.colorWriteMask = VK_COLOR_COMPONENT_R_BIT | VK_COLOR_COMPONENT_G_BIT |
                                     VK_COLOR_COMPONENT_B_BIT | VK_COLOR_COMPONENT_A_BIT;
    VkPipelineColorBlendStateCreateInfo blend{
        VK_STRUCTURE_TYPE_PIPELINE_COLOR_BLEND_STATE_CREATE_INFO};
    blend.attachmentCount = 1;
    blend.pAttachments = &blendAttachment;

    VkGraphicsPipelineCreateInfo pipelineInfo{VK_STRUCTURE_TYPE_GRAPHICS_PIPELINE_CREATE_INFO};
    pipelineInfo.stageCount = 2;
    pipelineInfo.pStages = stages;
    pipelineInfo.pVertexInputState = &vertexInput;
    pipelineInfo.pInputAssemblyState = &assembly;
    pipelineInfo.pViewportState = &viewportState;
    pipelineInfo.pRasterizationState = &raster;
    pipelineInfo.pMultisampleState = &multisample;
    pipelineInfo.pDepthStencilState = &depthStencil;
    pipelineInfo.pColorBlendState = &blend;
    pipelineInfo.layout = pipelineLayout_;
    pipelineInfo.renderPass = renderPass_;

    const VkResult created = vkCreateGraphicsPipelines(device.handle(), VK_NULL_HANDLE, 1,
                                                       &pipelineInfo, nullptr, &pipeline_);
    device.destroyShader(vertexModule);
    device.destroyShader(fragmentModule);
    if (created != VK_SUCCESS) {
        error = "vkCreateGraphicsPipelines failed";
        destroy();
        return false;
    }

    readback_ = device.createBuffer(std::uint64_t{width} * height * 4,
                                    VK_BUFFER_USAGE_TRANSFER_DST_BIT,
                                    VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT |
                                        VK_MEMORY_PROPERTY_HOST_COHERENT_BIT);
    if (!readback_.valid()) {
        error = "could not allocate the readback buffer";
        destroy();
        return false;
    }

    // Timestamps around the render pass. A failure here is not fatal: the frame
    // still renders, `FrameCost::gpuSeconds` just stays zero rather than
    // reporting a number nobody measured.
    VkQueryPoolCreateInfo queryInfo{VK_STRUCTURE_TYPE_QUERY_POOL_CREATE_INFO};
    queryInfo.queryType = VK_QUERY_TYPE_TIMESTAMP;
    queryInfo.queryCount = 2;
    if (vkCreateQueryPool(device.handle(), &queryInfo, nullptr, &queryPool_) != VK_SUCCESS)
        queryPool_ = nullptr;
    return true;
}

void HullRenderer::destroy() {
    if (device_ == nullptr) return;
    VkDevice handle = device_->handle();
    if (handle != nullptr) {
        vkDeviceWaitIdle(handle);
        if (queryPool_ != nullptr) vkDestroyQueryPool(handle, queryPool_, nullptr);
        if (pipeline_ != nullptr) vkDestroyPipeline(handle, pipeline_, nullptr);
        if (pipelineLayout_ != nullptr) vkDestroyPipelineLayout(handle, pipelineLayout_, nullptr);
        if (descriptorPool_ != nullptr) vkDestroyDescriptorPool(handle, descriptorPool_, nullptr);
        if (setLayout_ != nullptr) vkDestroyDescriptorSetLayout(handle, setLayout_, nullptr);
        if (framebuffer_ != nullptr) vkDestroyFramebuffer(handle, framebuffer_, nullptr);
        if (renderPass_ != nullptr) vkDestroyRenderPass(handle, renderPass_, nullptr);
        device_->destroyBuffer(vertexBuffer_);
        device_->destroyBuffer(indexBuffer_);
        device_->destroyBuffer(materialBuffer_);
        device_->destroyBuffer(readback_);
        device_->destroyImage(colour_);
        device_->destroyImage(depth_);
    }
    queryPool_ = nullptr;
    pipeline_ = nullptr;
    pipelineLayout_ = nullptr;
    descriptorSet_ = nullptr;
    descriptorPool_ = nullptr;
    setLayout_ = nullptr;
    framebuffer_ = nullptr;
    renderPass_ = nullptr;
    device_ = nullptr;
    width_ = height_ = 0;
    materialRevision_ = 0;
    materialCapacity_ = 0;
}

bool HullRenderer::ensureGeometryCapacity(std::size_t vertexBytes, std::size_t indexBytes) {
    if (vertexBytes > vertexBuffer_.size) {
        device_->destroyBuffer(vertexBuffer_);
        vertexBuffer_ = device_->createBuffer(
            vertexBytes, VK_BUFFER_USAGE_VERTEX_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT,
            VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT);
        if (!vertexBuffer_.valid()) return false;
    }
    if (indexBytes > indexBuffer_.size) {
        device_->destroyBuffer(indexBuffer_);
        indexBuffer_ = device_->createBuffer(
            indexBytes, VK_BUFFER_USAGE_INDEX_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT,
            VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT);
        if (!indexBuffer_.valid()) return false;
    }
    return true;
}

bool HullRenderer::ensureMaterials(const MaterialLibrary& library) {
    // The whole table, however large. The MaterialId readback channel is 16 bits
    // and cannot report an index past 65535, but truncating the table here to suit
    // it would break *shading* to keep a measurement channel readable, which is
    // the wrong way round.
    const std::size_t count = library.packed().size();
    // The shader indexes this buffer unconditionally, so it must exist even for an
    // empty library. One default row is cheaper than a branch in the fragment
    // shader and cannot be reached by a mesh, which carries real indices.
    const std::size_t rows = std::max<std::size_t>(count, 1);
    if (materialBuffer_.valid() && rows <= materialCapacity_ &&
        materialRevision_ == library.revision())
        return true;

    if (!materialBuffer_.valid() || rows > materialCapacity_) {
        vkDeviceWaitIdle(device_->handle());
        device_->destroyBuffer(materialBuffer_);
        materialBuffer_ = device_->createBuffer(
            rows * sizeof(GpuMaterial),
            VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT,
            VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT);
        if (!materialBuffer_.valid()) return false;
        materialCapacity_ = rows;

        VkDescriptorBufferInfo bufferInfo{materialBuffer_.handle, 0, VK_WHOLE_SIZE};
        VkWriteDescriptorSet write{VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET};
        write.dstSet = descriptorSet_;
        write.dstBinding = 0;
        write.descriptorCount = 1;
        write.descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
        write.pBufferInfo = &bufferInfo;
        vkUpdateDescriptorSets(device_->handle(), 1, &write, 0, nullptr);
    }

    std::vector<GpuMaterial> rowsToUpload(rows, GpuMaterial{});
    std::copy_n(library.packed().begin(), count, rowsToUpload.begin());
    if (!device_->upload(materialBuffer_, rowsToUpload.data(), rows * sizeof(GpuMaterial)))
        return false;
    materialRevision_ = library.revision();
    return true;
}

bool HullRenderer::render(const float mvp[16], const SceneView& view, const SceneMesh& mesh,
                          const MaterialLibrary& library, const float clearColour[4],
                          core::Image& out) {
    if (!valid()) return false;
    const auto frameStarted = std::chrono::steady_clock::now();
    lastFrame_ = FrameCost{};

    const std::size_t vertexCount = mesh.vertices().size();
    const std::size_t indexCount = mesh.indices().size();
    const bool hasGeometry = vertexCount > 0 && indexCount > 0;
    lastFrame_.vertices = vertexCount;
    lastFrame_.triangles = indexCount / 3;

    const auto uploadStarted = std::chrono::steady_clock::now();
    if (!ensureMaterials(library)) return false;
    if (hasGeometry) {
        const std::size_t vertexBytes = vertexCount * sizeof(HullVertex);
        const std::size_t indexBytes = indexCount * sizeof(std::uint32_t);
        if (!ensureGeometryCapacity(vertexBytes, indexBytes)) return false;
        if (!device_->upload(vertexBuffer_, mesh.vertices().data(), vertexBytes)) return false;
        if (!device_->upload(indexBuffer_, mesh.indices().data(), indexBytes)) return false;
    }
    lastFrame_.uploadSeconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - uploadStarted).count();

    Push push{};
    std::memcpy(push.modelViewProjection, mvp, sizeof(push.modelViewProjection));
    // Normalised here rather than trusted, for the reason OceanRenderer gives: a
    // direction of length 1.001 is lighting that is quietly 0.1% too bright.
    const double length = std::sqrt(
        static_cast<double>(view.sunDirection[0]) * view.sunDirection[0] +
        static_cast<double>(view.sunDirection[1]) * view.sunDirection[1] +
        static_cast<double>(view.sunDirection[2]) * view.sunDirection[2]);
    const double inverse = length > 1e-12 ? 1.0 / length : 0.0;
    for (int i = 0; i < 3; ++i)
        push.sun[i] = static_cast<float>(view.sunDirection[i] * inverse);
    push.sun[3] = view.mode == HullShadingMode::MaterialId ? 1.0f
                  : view.mode == HullShadingMode::Depth    ? 2.0f
                                                           : 0.0f;
    for (int i = 0; i < 3; ++i) push.sunColour[i] = view.sunColour[i];
    push.sunColour[3] = view.exposure;
    for (int i = 0; i < 3; ++i) push.sky[i] = view.skyColour[i];
    push.sky[3] = 0.0f;
    for (int i = 0; i < 3; ++i) push.eye[i] = view.eye[i];
    push.eye[3] = 0.0f;

    VkClearValue clears[2]{};
    for (int i = 0; i < 4; ++i) clears[0].color.float32[i] = clearColour[i];
    clears[1].depthStencil.depth = 1.0f;

    const auto submitStarted = std::chrono::steady_clock::now();
    VkCommandBuffer commands = device_->beginOneShot();

    if (queryPool_ != nullptr) {
        vkCmdResetQueryPool(commands, queryPool_, 0, 2);
        vkCmdWriteTimestamp(commands, VK_PIPELINE_STAGE_TOP_OF_PIPE_BIT, queryPool_, 0);
    }

    VkRenderPassBeginInfo begin{VK_STRUCTURE_TYPE_RENDER_PASS_BEGIN_INFO};
    begin.renderPass = renderPass_;
    begin.framebuffer = framebuffer_;
    begin.renderArea = {{0, 0}, {width_, height_}};
    begin.clearValueCount = 2;
    begin.pClearValues = clears;
    vkCmdBeginRenderPass(commands, &begin, VK_SUBPASS_CONTENTS_INLINE);

    if (hasGeometry) {
        vkCmdBindPipeline(commands, VK_PIPELINE_BIND_POINT_GRAPHICS, pipeline_);
        vkCmdBindDescriptorSets(commands, VK_PIPELINE_BIND_POINT_GRAPHICS, pipelineLayout_, 0, 1,
                                &descriptorSet_, 0, nullptr);
        vkCmdPushConstants(commands, pipelineLayout_,
                           VK_SHADER_STAGE_VERTEX_BIT | VK_SHADER_STAGE_FRAGMENT_BIT, 0,
                           sizeof(Push), &push);
        const VkDeviceSize offset = 0;
        vkCmdBindVertexBuffers(commands, 0, 1, &vertexBuffer_.handle, &offset);
        vkCmdBindIndexBuffer(commands, indexBuffer_.handle, 0, VK_INDEX_TYPE_UINT32);
        vkCmdDrawIndexed(commands, static_cast<std::uint32_t>(indexCount), 1, 0, 0, 0);
    }
    vkCmdEndRenderPass(commands);

    if (queryPool_ != nullptr)
        vkCmdWriteTimestamp(commands, VK_PIPELINE_STAGE_BOTTOM_OF_PIPE_BIT, queryPool_, 1);

    VkBufferImageCopy copy{};
    copy.imageSubresource.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
    copy.imageSubresource.layerCount = 1;
    copy.imageExtent = {width_, height_, 1};
    vkCmdCopyImageToBuffer(commands, colour_.handle, VK_IMAGE_LAYOUT_TRANSFER_SRC_OPTIMAL,
                           readback_.handle, 1, &copy);

    device_->endOneShot(commands);
    lastFrame_.submitSeconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - submitStarted).count();

    if (queryPool_ != nullptr) {
        std::uint64_t stamps[2] = {0, 0};
        if (vkGetQueryPoolResults(device_->handle(), queryPool_, 0, 2, sizeof(stamps), stamps,
                                  sizeof(std::uint64_t),
                                  VK_QUERY_RESULT_64_BIT | VK_QUERY_RESULT_WAIT_BIT) ==
                VK_SUCCESS &&
            stamps[1] >= stamps[0])
            lastFrame_.gpuSeconds = static_cast<double>(stamps[1] - stamps[0]) *
                                    static_cast<double>(device_->timestampPeriodNanos()) * 1e-9;
    }

    const auto readbackStarted = std::chrono::steady_clock::now();
    out = core::Image(width_, height_);
    void* mapped = nullptr;
    if (vkMapMemory(device_->handle(), readback_.memory, 0, out.rgba.size(), 0, &mapped) !=
        VK_SUCCESS)
        return false;
    std::memcpy(out.rgba.data(), mapped, out.rgba.size());
    vkUnmapMemory(device_->handle(), readback_.memory);
    lastFrame_.readbackSeconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - readbackStarted).count();
    lastFrame_.totalSeconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - frameStarted).count();
    return true;
}

}  // namespace gpu
