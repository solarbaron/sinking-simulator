// SPDX-License-Identifier: MIT
#include "ocean.hpp"

#include <vulkan/vulkan.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <vector>

namespace gpu {
namespace {

constexpr VkFormat kColourFormat = VK_FORMAT_R8G8B8A8_UNORM;
constexpr VkFormat kDepthFormat = VK_FORMAT_D32_SFLOAT;

// A patch is a square grid, so the vertex count is the square of this. 2048 is
// 4.2 M vertices, well past anything the CPU path can displace in a frame; it
// exists to turn a nonsense argument into a clamp rather than an allocation
// failure.
constexpr int kMaxResolution = 2048;

// Mirrors the `Push` block in shaders/ocean.vert and shaders/ocean.frag. 112
// bytes, inside the 128 every Vulkan implementation guarantees, so the ocean
// needs no descriptor set at all.
struct Push {
    float modelViewProjection[16];
    float sun[4];     // xyz unit vector toward the sun, w strength
    float water[4];   // rgb water colour, w ambient
    float encode[4];  // x mode (0 shaded, 1 elevation), y zMin, z 1 / span, w unused
};
static_assert(sizeof(Push) == 112, "the push block must match the shaders");

}  // namespace

// --- Resolution ---------------------------------------------------------------

int oceanResolutionFor(double halfExtent, double wavelength, double cellsPerWavelength) {
    if (!(halfExtent > 0.0) || !(wavelength > 0.0) || !(cellsPerWavelength > 0.0)) return 1;
    const double cells = 2.0 * halfExtent * cellsPerWavelength / wavelength;
    if (!(cells > 1.0)) return 1;
    if (cells >= kMaxResolution) return kMaxResolution;
    return static_cast<int>(std::ceil(cells));
}

double dominantWavelength(const sim::WaveField& field) {
    // Equal-energy bins, so the peak of the density is the narrowest bin. Same
    // route tests/test_waves.cpp uses to find Tp in the discretisation.
    double bestDensity = -1.0;
    double bestOmega = 0.0;
    for (const sim::FrequencyBin& bin : field.frequencyBins()) {
        const double width = bin.omegaHigh - bin.omegaLow;
        if (!std::isfinite(width) || !(width > 0.0)) continue;  // the open-topped bin
        const double density = bin.energy / width;
        if (density > bestDensity) {
            bestDensity = density;
            bestOmega = bin.omega;
        }
    }
    if (bestDensity < 0.0)
        // Every bin is open-topped, which happens only when a sea state has a
        // single frequency: that bin runs from zero to infinity and so has no
        // density at all. Its centroid is the only frequency there is.
        for (const sim::FrequencyBin& bin : field.frequencyBins())
            if (bin.omega > 0.0) {
                bestOmega = bin.omega;
                break;
            }
    const double k = sim::deepWaterWavenumber(bestOmega);
    return k > 0.0 ? 2.0 * sim::kPi / k : 0.0;
}

double shortestWavelength(const sim::WaveField& field) {
    double maxWavenumber = 0.0;
    for (const sim::WaveComponent& component : field.components())
        maxWavenumber = std::max(maxWavenumber, component.wavenumber);
    return maxWavenumber > 0.0 ? 2.0 * sim::kPi / maxWavenumber : 0.0;
}

// --- The displaced grid -------------------------------------------------------

std::size_t OceanSurface::vertexIndex(int i, int j) const {
    return static_cast<std::size_t>(j) * static_cast<std::size_t>(grid_.resolution + 1) +
           static_cast<std::size_t>(i);
}

void OceanSurface::build(const sim::WaveField& field, const OceanGrid& grid, double time) {
    const auto started = std::chrono::steady_clock::now();

    grid_ = grid;
    grid_.resolution = std::clamp(grid.resolution, 1, kMaxResolution);
    time_ = time;

    const int cells = grid_.resolution;
    const int side = cells + 1;
    const double h = grid_.cellSize();
    const double x0 = grid_.centreX - grid_.halfExtent;
    const double y0 = grid_.centreY - grid_.halfExtent;

    vertices_.resize(static_cast<std::size_t>(side) * static_cast<std::size_t>(side));

    if (indexResolution_ != cells) {
        indices_.clear();
        indices_.reserve(static_cast<std::size_t>(cells) * static_cast<std::size_t>(cells) * 6);
        for (int j = 0; j < cells; ++j)
            for (int i = 0; i < cells; ++i) {
                const auto v00 = static_cast<std::uint32_t>(vertexIndex(i, j));
                const auto v10 = static_cast<std::uint32_t>(v00 + 1);
                const auto v01 = static_cast<std::uint32_t>(v00 + static_cast<std::uint32_t>(side));
                const auto v11 = static_cast<std::uint32_t>(v01 + 1);
                // Split on the v00-v11 diagonal, counter-clockwise seen from
                // above. The choice is load-bearing for the tests: the cell
                // centre then lies on the shared edge, where the interpolated
                // height is the mean of the two diagonal corners -- which is what
                // makes the linear-interpolation error a closed form rather than
                // a bound.
                indices_.push_back(v00);
                indices_.push_back(v10);
                indices_.push_back(v11);
                indices_.push_back(v00);
                indices_.push_back(v11);
                indices_.push_back(v01);
            }
        indexResolution_ = cells;
    }

    rowElevation_.resize(static_cast<std::size_t>(side));
    rowSlopeX_.resize(static_cast<std::size_t>(side));
    rowSlopeY_.resize(static_cast<std::size_t>(side));

    const std::vector<sim::WaveComponent>& components = field.components();

    // One evaluation per *unique grid vertex*. Six triangle corners share each
    // interior vertex, so evaluating per index instead would multiply the whole
    // cost by six -- which is exactly the redundancy the physics tick was just
    // caught paying (CLAUDE.md, "Sea surface queried 6x more than necessary").
    for (int j = 0; j < side; ++j) {
        const double y = y0 + h * j;
        std::fill(rowElevation_.begin(), rowElevation_.end(), 0.0);
        std::fill(rowSlopeX_.begin(), rowSlopeX_.end(), 0.0);
        std::fill(rowSlopeY_.begin(), rowSlopeY_.end(), 0.0);

        for (const sim::WaveComponent& c : components) {
            const double kx = c.wavenumber * c.dirX;
            const double ky = c.wavenumber * c.dirY;
            // eta = a cos(psi), so d eta/dx = -a kx sin(psi) and likewise in y.
            // The recurrence produces sin(psi) as a by-product, so the slope --
            // and therefore the normal -- is free.
            const double a = c.amplitude;
            const double slopeGainX = -a * kx;
            const double slopeGainY = -a * ky;

            const double psi = kx * x0 + ky * y - c.omega * time + c.phase;
            double cosPsi = std::cos(psi);
            double sinPsi = std::sin(psi);
            // Along a row only x moves, so psi advances by a constant and
            // (cos, sin) can be stepped by a fixed rotation. Restarted from a
            // real sincos on every row, so drift is bounded by the row length
            // rather than by the whole grid; measured worst error against direct
            // evaluation over a 257-wide row is 2e-14 m.
            const double stepCos = std::cos(kx * h);
            const double stepSin = std::sin(kx * h);

            for (int i = 0; i < side; ++i) {
                const auto index = static_cast<std::size_t>(i);
                rowElevation_[index] += a * cosPsi;
                rowSlopeX_[index] += slopeGainX * sinPsi;
                rowSlopeY_[index] += slopeGainY * sinPsi;
                const double nextCos = cosPsi * stepCos - sinPsi * stepSin;
                sinPsi = sinPsi * stepCos + cosPsi * stepSin;
                cosPsi = nextCos;
            }
        }

        for (int i = 0; i < side; ++i) {
            const auto index = static_cast<std::size_t>(i);
            // n = normalize(-d eta/dx, -d eta/dy, 1).
            const double nx = -rowSlopeX_[index];
            const double ny = -rowSlopeY_[index];
            const double inverse = 1.0 / std::sqrt(nx * nx + ny * ny + 1.0);
            OceanVertex& vertex = vertices_[vertexIndex(i, j)];
            vertex.position[0] = static_cast<float>(x0 + h * i);
            vertex.position[1] = static_cast<float>(y);
            vertex.position[2] = static_cast<float>(rowElevation_[index]);
            vertex.normal[0] = static_cast<float>(nx * inverse);
            vertex.normal[1] = static_cast<float>(ny * inverse);
            vertex.normal[2] = static_cast<float>(inverse);
        }
    }

    buildSeconds_ = std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
}

bool OceanSurface::sampleElevation(double x, double y, double& elevation) const {
    const int cells = grid_.resolution;
    if (vertices_.empty() || cells < 1) return false;
    const double h = grid_.cellSize();
    const double fx = (x - (grid_.centreX - grid_.halfExtent)) / h;
    const double fy = (y - (grid_.centreY - grid_.halfExtent)) / h;
    // Written so a NaN falls out of the patch rather than into an index.
    if (!(fx >= 0.0) || !(fy >= 0.0) || !(fx <= cells) || !(fy <= cells)) return false;

    int i = static_cast<int>(std::floor(fx));
    int j = static_cast<int>(std::floor(fy));
    if (i >= cells) i = cells - 1;
    if (j >= cells) j = cells - 1;
    const double u = fx - i;
    const double v = fy - j;

    const auto height = [&](int a, int b) {
        return static_cast<double>(vertices_[vertexIndex(a, b)].position[2]);
    };
    const double z00 = height(i, j);
    const double z11 = height(i + 1, j + 1);
    // The v00-v11 diagonal splits the cell along u == v; below it the triangle is
    // (v00, v10, v11) and above it (v00, v11, v01). Both agree on the diagonal,
    // where the value is (z00 + z11) / 2.
    elevation = v <= u ? z00 + u * (height(i + 1, j) - z00) + v * (z11 - height(i + 1, j))
                       : z00 + v * (height(i, j + 1) - z00) + u * (z11 - height(i, j + 1));
    return true;
}

// --- The elevation channel ----------------------------------------------------

bool decodeOceanElevation(const OceanView& view, const std::uint8_t* pixel, double& elevation) {
    // Blue is the surface tag. The clear colour has to stay away from it, which
    // is cheap to arrange and much safer than trying to tell a legitimate
    // elevation code apart from the background.
    if (pixel == nullptr || pixel[2] < 128) return false;
    const double code = static_cast<double>(pixel[0]) * 256.0 + static_cast<double>(pixel[1]);
    elevation = static_cast<double>(view.elevationMin) +
                code / 65535.0 * static_cast<double>(view.elevationSpan);
    return true;
}

// --- The pipeline -------------------------------------------------------------

OceanRenderer::~OceanRenderer() { destroy(); }

bool OceanRenderer::create(Device& device, std::uint32_t width, std::uint32_t height,
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

    // Same layout chain as OffscreenRenderer: the pass leaves the colour
    // attachment in TRANSFER_SRC_OPTIMAL so the readback copy needs no barrier.
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

    VkShaderModule vertexModule = device.loadShader(shaderDirectory + "/ocean.vert.spv", error);
    if (vertexModule == nullptr) { destroy(); return false; }
    VkShaderModule fragmentModule = device.loadShader(shaderDirectory + "/ocean.frag.spv", error);
    if (fragmentModule == nullptr) {
        device.destroyShader(vertexModule);
        destroy();
        return false;
    }

    // One range spanning both stages: the vertex shader reads the matrix, the
    // fragment shader the lighting and the encoding, and a single range keeps the
    // block declaration identical in the two shaders.
    VkPushConstantRange push{VK_SHADER_STAGE_VERTEX_BIT | VK_SHADER_STAGE_FRAGMENT_BIT, 0,
                             sizeof(Push)};
    VkPipelineLayoutCreateInfo layoutInfo{VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO};
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

    VkVertexInputBindingDescription binding{0, sizeof(OceanVertex), VK_VERTEX_INPUT_RATE_VERTEX};
    VkVertexInputAttributeDescription vertexAttributes[2]{};
    vertexAttributes[0] = {0, 0, VK_FORMAT_R32G32B32_SFLOAT, offsetof(OceanVertex, position)};
    vertexAttributes[1] = {1, 0, VK_FORMAT_R32G32B32_SFLOAT, offsetof(OceanVertex, normal)};

    VkPipelineVertexInputStateCreateInfo vertexInput{
        VK_STRUCTURE_TYPE_PIPELINE_VERTEX_INPUT_STATE_CREATE_INFO};
    vertexInput.vertexBindingDescriptionCount = 1;
    vertexInput.pVertexBindingDescriptions = &binding;
    vertexInput.vertexAttributeDescriptionCount = 2;
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
    // The sea is seen from below as often as from above -- from a flooded
    // compartment, from underwater -- so neither face is a back face. The
    // fragment shader flips the normal for the far side instead.
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
    return true;
}

void OceanRenderer::destroy() {
    if (device_ == nullptr) return;
    VkDevice handle = device_->handle();
    if (handle != nullptr) {
        vkDeviceWaitIdle(handle);
        if (pipeline_ != nullptr) vkDestroyPipeline(handle, pipeline_, nullptr);
        if (pipelineLayout_ != nullptr) vkDestroyPipelineLayout(handle, pipelineLayout_, nullptr);
        if (framebuffer_ != nullptr) vkDestroyFramebuffer(handle, framebuffer_, nullptr);
        if (renderPass_ != nullptr) vkDestroyRenderPass(handle, renderPass_, nullptr);
        device_->destroyBuffer(vertexBuffer_);
        device_->destroyBuffer(indexBuffer_);
        device_->destroyBuffer(readback_);
        device_->destroyImage(colour_);
        device_->destroyImage(depth_);
    }
    pipeline_ = nullptr;
    pipelineLayout_ = nullptr;
    framebuffer_ = nullptr;
    renderPass_ = nullptr;
    device_ = nullptr;
    width_ = height_ = 0;
}

bool OceanRenderer::ensureGeometryCapacity(std::size_t vertexBytes, std::size_t indexBytes) {
    if (vertexBytes > vertexBuffer_.size) {
        device_->destroyBuffer(vertexBuffer_);
        vertexBuffer_ = device_->createBuffer(
            vertexBytes,
            VK_BUFFER_USAGE_VERTEX_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT,
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

bool OceanRenderer::render(const float mvp[16], const OceanSurface& surface, const OceanView& view,
                           const float clearColour[4], core::Image& out) {
    if (!valid()) return false;

    const std::size_t vertexCount = surface.vertices().size();
    const std::size_t indexCount = surface.indices().size();
    const bool hasGeometry = vertexCount > 0 && indexCount > 0;
    if (hasGeometry) {
        const std::size_t vertexBytes = vertexCount * sizeof(OceanVertex);
        const std::size_t indexBytes = indexCount * sizeof(std::uint32_t);
        if (!ensureGeometryCapacity(vertexBytes, indexBytes)) return false;
        if (!device_->upload(vertexBuffer_, surface.vertices().data(), vertexBytes)) return false;
        if (!device_->upload(indexBuffer_, surface.indices().data(), indexBytes)) return false;
    }

    Push push{};
    std::memcpy(push.modelViewProjection, mvp, sizeof(push.modelViewProjection));
    // Normalised here rather than trusted: a caller who hands over a direction of
    // length 1.001 would otherwise get lighting that is subtly too bright, which
    // is the sort of thing nobody ever notices.
    const double length = std::sqrt(
        static_cast<double>(view.sunDirection[0]) * view.sunDirection[0] +
        static_cast<double>(view.sunDirection[1]) * view.sunDirection[1] +
        static_cast<double>(view.sunDirection[2]) * view.sunDirection[2]);
    const double inverse = length > 1e-12 ? 1.0 / length : 0.0;
    for (int i = 0; i < 3; ++i)
        push.sun[i] = static_cast<float>(view.sunDirection[i] * inverse);
    push.sun[3] = view.sunStrength;
    for (int i = 0; i < 3; ++i) push.water[i] = view.waterColour[i];
    push.water[3] = view.ambient;
    push.encode[0] = view.shading == OceanShading::Elevation ? 1.0f : 0.0f;
    push.encode[1] = view.elevationMin;
    push.encode[2] = view.elevationSpan != 0.0f ? 1.0f / view.elevationSpan : 0.0f;
    push.encode[3] = 0.0f;

    VkClearValue clears[2]{};
    for (int i = 0; i < 4; ++i) clears[0].color.float32[i] = clearColour[i];
    clears[1].depthStencil.depth = 1.0f;

    VkCommandBuffer commands = device_->beginOneShot();

    VkRenderPassBeginInfo begin{VK_STRUCTURE_TYPE_RENDER_PASS_BEGIN_INFO};
    begin.renderPass = renderPass_;
    begin.framebuffer = framebuffer_;
    begin.renderArea = {{0, 0}, {width_, height_}};
    begin.clearValueCount = 2;
    begin.pClearValues = clears;
    vkCmdBeginRenderPass(commands, &begin, VK_SUBPASS_CONTENTS_INLINE);

    if (hasGeometry) {
        vkCmdBindPipeline(commands, VK_PIPELINE_BIND_POINT_GRAPHICS, pipeline_);
        vkCmdPushConstants(commands, pipelineLayout_,
                           VK_SHADER_STAGE_VERTEX_BIT | VK_SHADER_STAGE_FRAGMENT_BIT, 0,
                           sizeof(Push), &push);
        const VkDeviceSize offset = 0;
        vkCmdBindVertexBuffers(commands, 0, 1, &vertexBuffer_.handle, &offset);
        vkCmdBindIndexBuffer(commands, indexBuffer_.handle, 0, VK_INDEX_TYPE_UINT32);
        vkCmdDrawIndexed(commands, static_cast<std::uint32_t>(indexCount), 1, 0, 0, 0);
    }
    vkCmdEndRenderPass(commands);

    VkBufferImageCopy copy{};
    copy.imageSubresource.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
    copy.imageSubresource.layerCount = 1;
    copy.imageExtent = {width_, height_, 1};
    vkCmdCopyImageToBuffer(commands, colour_.handle, VK_IMAGE_LAYOUT_TRANSFER_SRC_OPTIMAL,
                           readback_.handle, 1, &copy);

    device_->endOneShot(commands);

    out = core::Image(width_, height_);
    void* mapped = nullptr;
    if (vkMapMemory(device_->handle(), readback_.memory, 0, out.rgba.size(), 0, &mapped) !=
        VK_SUCCESS)
        return false;
    std::memcpy(out.rgba.data(), mapped, out.rgba.size());
    vkUnmapMemory(device_->handle(), readback_.memory);
    return true;
}

}  // namespace gpu
