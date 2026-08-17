// SPDX-License-Identifier: MIT
//
// The two-pass offscreen target: the lit solid, then the participating medium
// over it.
//
// Split from `smoke.cpp` for the reason `hull.cpp` is split from `damage.cpp` --
// everything in the other file is arithmetic that a machine with no GPU must
// still be able to run and check.
#include "smoke.hpp"

#include <vulkan/vulkan.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstring>

namespace gpu {
namespace {

constexpr VkFormat kColourFormat = VK_FORMAT_R8G8B8A8_UNORM;
constexpr VkFormat kDepthFormat = VK_FORMAT_D32_SFLOAT;

// Mirrors the `Push` block in shaders/hull.vert and shaders/hull.frag. The lit
// pass here is `HullRenderer`'s pass -- the same SPIR-V, the same push block --
// because an empty volume list has to produce the same bytes it does.
struct OpaquePush {
    float modelViewProjection[16];
    float sun[4];
    float sunColour[4];
    float sky[4];
    float eye[4];
};
static_assert(sizeof(OpaquePush) == 128, "the push block must match shaders/hull.*");

// Mirrors the `Push` block in shaders/smoke.vert and shaders/smoke.frag.
struct SmokePush {
    float modelViewProjection[16];
    float eyeA[4];      // xyz world eye, w the depth constant a
    float forwardB[4];  // xyz unit view direction, w the depth constant b
};
static_assert(sizeof(SmokePush) == 96, "the push block must match shaders/smoke.*");

// Mirrors `struct Volume` in shaders/smoke.{vert,frag}: seven vec4s, which std430
// packs tightly.
struct GpuSmokeVolume {
    float loInterface[4];
    float hi[4];
    float row0[4];
    float row1[4];
    float row2[4];
    float upper[4];
    float lower[4];
};
static_assert(sizeof(GpuSmokeVolume) == 112, "the volume row must match shaders/smoke.*");

GpuSmokeVolume pack(const SmokeVolume& v) {
    GpuSmokeVolume out{};
    out.loInterface[0] = static_cast<float>(v.lo.x);
    out.loInterface[1] = static_cast<float>(v.lo.y);
    out.loInterface[2] = static_cast<float>(v.lo.z);
    out.loInterface[3] = static_cast<float>(v.interfaceZ);
    out.hi[0] = static_cast<float>(v.hi.x);
    out.hi[1] = static_cast<float>(v.hi.y);
    out.hi[2] = static_cast<float>(v.hi.z);
    out.hi[3] = 0.0f;
    float* rows[3] = {out.row0, out.row1, out.row2};
    for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 3; ++c) rows[r][c] = static_cast<float>(v.rotation(r, c));
        rows[r][3] = static_cast<float>(v.translation[r]);
    }
    for (int c = 0; c < 3; ++c) {
        out.upper[c] = static_cast<float>(v.upper.emission[c]);
        out.lower[c] = static_cast<float>(v.lower.emission[c]);
    }
    out.upper[3] = static_cast<float>(v.upper.extinction);
    out.lower[3] = static_cast<float>(v.lower.extinction);
    return out;
}

}  // namespace

SmokeRenderer::~SmokeRenderer() { destroy(); }

bool SmokeRenderer::create(Device& device, std::uint32_t width, std::uint32_t height,
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

    colour_ = device.createImage2D(
        width, height, kColourFormat,
        VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT | VK_IMAGE_USAGE_TRANSFER_SRC_BIT,
        VK_IMAGE_ASPECT_COLOR_BIT);
    // Sampled as well as written: the volumetric pass reads it to find where each
    // ray stops, which is a clamp on an integration limit and cannot be a depth
    // test.
    depth_ = device.createImage2D(
        width, height, kDepthFormat,
        VK_IMAGE_USAGE_DEPTH_STENCIL_ATTACHMENT_BIT | VK_IMAGE_USAGE_SAMPLED_BIT,
        VK_IMAGE_ASPECT_DEPTH_BIT);
    if (!colour_.valid() || !depth_.valid()) {
        error = "could not allocate the render target";
        destroy();
        return false;
    }

    // --- pass one: the lit solid ---------------------------------------------
    //
    // Depth is STOREd, which is the one thing that differs from HullRenderer's
    // pass: it discards depth because nothing downstream of it reads one.
    VkAttachmentDescription opaqueAttachments[2]{};
    opaqueAttachments[0].format = kColourFormat;
    opaqueAttachments[0].samples = VK_SAMPLE_COUNT_1_BIT;
    opaqueAttachments[0].loadOp = VK_ATTACHMENT_LOAD_OP_CLEAR;
    opaqueAttachments[0].storeOp = VK_ATTACHMENT_STORE_OP_STORE;
    opaqueAttachments[0].stencilLoadOp = VK_ATTACHMENT_LOAD_OP_DONT_CARE;
    opaqueAttachments[0].stencilStoreOp = VK_ATTACHMENT_STORE_OP_DONT_CARE;
    opaqueAttachments[0].initialLayout = VK_IMAGE_LAYOUT_UNDEFINED;
    opaqueAttachments[0].finalLayout = VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL;

    opaqueAttachments[1].format = kDepthFormat;
    opaqueAttachments[1].samples = VK_SAMPLE_COUNT_1_BIT;
    opaqueAttachments[1].loadOp = VK_ATTACHMENT_LOAD_OP_CLEAR;
    opaqueAttachments[1].storeOp = VK_ATTACHMENT_STORE_OP_STORE;
    opaqueAttachments[1].stencilLoadOp = VK_ATTACHMENT_LOAD_OP_DONT_CARE;
    opaqueAttachments[1].stencilStoreOp = VK_ATTACHMENT_STORE_OP_DONT_CARE;
    opaqueAttachments[1].initialLayout = VK_IMAGE_LAYOUT_UNDEFINED;
    opaqueAttachments[1].finalLayout = VK_IMAGE_LAYOUT_DEPTH_STENCIL_ATTACHMENT_OPTIMAL;

    VkAttachmentReference colourRef{0, VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL};
    VkAttachmentReference depthRef{1, VK_IMAGE_LAYOUT_DEPTH_STENCIL_ATTACHMENT_OPTIMAL};

    VkSubpassDescription opaqueSubpass{};
    opaqueSubpass.pipelineBindPoint = VK_PIPELINE_BIND_POINT_GRAPHICS;
    opaqueSubpass.colorAttachmentCount = 1;
    opaqueSubpass.pColorAttachments = &colourRef;
    opaqueSubpass.pDepthStencilAttachment = &depthRef;

    VkRenderPassCreateInfo opaqueInfo{VK_STRUCTURE_TYPE_RENDER_PASS_CREATE_INFO};
    opaqueInfo.attachmentCount = 2;
    opaqueInfo.pAttachments = opaqueAttachments;
    opaqueInfo.subpassCount = 1;
    opaqueInfo.pSubpasses = &opaqueSubpass;
    if (vkCreateRenderPass(device.handle(), &opaqueInfo, nullptr, &opaquePass_) != VK_SUCCESS) {
        error = "vkCreateRenderPass failed for the lit pass";
        destroy();
        return false;
    }

    // --- pass two: the medium, composited over what pass one left -------------
    VkAttachmentDescription smokeAttachment{};
    smokeAttachment.format = kColourFormat;
    smokeAttachment.samples = VK_SAMPLE_COUNT_1_BIT;
    smokeAttachment.loadOp = VK_ATTACHMENT_LOAD_OP_LOAD;
    smokeAttachment.storeOp = VK_ATTACHMENT_STORE_OP_STORE;
    smokeAttachment.stencilLoadOp = VK_ATTACHMENT_LOAD_OP_DONT_CARE;
    smokeAttachment.stencilStoreOp = VK_ATTACHMENT_STORE_OP_DONT_CARE;
    smokeAttachment.initialLayout = VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL;
    smokeAttachment.finalLayout = VK_IMAGE_LAYOUT_TRANSFER_SRC_OPTIMAL;

    VkAttachmentReference smokeColourRef{0, VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL};
    VkSubpassDescription smokeSubpass{};
    smokeSubpass.pipelineBindPoint = VK_PIPELINE_BIND_POINT_GRAPHICS;
    smokeSubpass.colorAttachmentCount = 1;
    smokeSubpass.pColorAttachments = &smokeColourRef;

    // **The one dependency that is not implicit.** Two render pass instances in the
    // same command buffer are ordered but not synchronised: the lit pass writes the
    // colour attachment and this one blends against it, which is a read-after-write
    // across the boundary. The depth image gets an explicit image barrier below
    // because it also changes layout; the colour image does not change layout, so
    // it would be silently unsynchronised without this.
    VkSubpassDependency smokeDependency{};
    smokeDependency.srcSubpass = VK_SUBPASS_EXTERNAL;
    smokeDependency.dstSubpass = 0;
    smokeDependency.srcStageMask = VK_PIPELINE_STAGE_COLOR_ATTACHMENT_OUTPUT_BIT;
    smokeDependency.dstStageMask = VK_PIPELINE_STAGE_COLOR_ATTACHMENT_OUTPUT_BIT;
    smokeDependency.srcAccessMask = VK_ACCESS_COLOR_ATTACHMENT_WRITE_BIT;
    smokeDependency.dstAccessMask =
        VK_ACCESS_COLOR_ATTACHMENT_READ_BIT | VK_ACCESS_COLOR_ATTACHMENT_WRITE_BIT;

    VkRenderPassCreateInfo smokeInfo{VK_STRUCTURE_TYPE_RENDER_PASS_CREATE_INFO};
    smokeInfo.attachmentCount = 1;
    smokeInfo.pAttachments = &smokeAttachment;
    smokeInfo.subpassCount = 1;
    smokeInfo.pSubpasses = &smokeSubpass;
    smokeInfo.dependencyCount = 1;
    smokeInfo.pDependencies = &smokeDependency;
    if (vkCreateRenderPass(device.handle(), &smokeInfo, nullptr, &smokePass_) != VK_SUCCESS) {
        error = "vkCreateRenderPass failed for the volumetric pass";
        destroy();
        return false;
    }

    VkImageView opaqueViews[2] = {colour_.view, depth_.view};
    VkFramebufferCreateInfo framebufferInfo{VK_STRUCTURE_TYPE_FRAMEBUFFER_CREATE_INFO};
    framebufferInfo.renderPass = opaquePass_;
    framebufferInfo.attachmentCount = 2;
    framebufferInfo.pAttachments = opaqueViews;
    framebufferInfo.width = width;
    framebufferInfo.height = height;
    framebufferInfo.layers = 1;
    if (vkCreateFramebuffer(device.handle(), &framebufferInfo, nullptr, &opaqueFramebuffer_) !=
        VK_SUCCESS) {
        error = "vkCreateFramebuffer failed for the lit pass";
        destroy();
        return false;
    }
    framebufferInfo.renderPass = smokePass_;
    framebufferInfo.attachmentCount = 1;
    framebufferInfo.pAttachments = &colour_.view;
    if (vkCreateFramebuffer(device.handle(), &framebufferInfo, nullptr, &smokeFramebuffer_) !=
        VK_SUCCESS) {
        error = "vkCreateFramebuffer failed for the volumetric pass";
        destroy();
        return false;
    }

    // --- descriptors ----------------------------------------------------------
    VkDescriptorSetLayoutBinding materialBinding{};
    materialBinding.binding = 0;
    materialBinding.descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
    materialBinding.descriptorCount = 1;
    materialBinding.stageFlags = VK_SHADER_STAGE_FRAGMENT_BIT;
    VkDescriptorSetLayoutCreateInfo materialLayoutInfo{
        VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO};
    materialLayoutInfo.bindingCount = 1;
    materialLayoutInfo.pBindings = &materialBinding;
    if (vkCreateDescriptorSetLayout(device.handle(), &materialLayoutInfo, nullptr,
                                    &materialLayout_) != VK_SUCCESS) {
        error = "vkCreateDescriptorSetLayout failed for the materials";
        destroy();
        return false;
    }

    VkDescriptorSetLayoutBinding volumeBindings[2]{};
    volumeBindings[0].binding = 0;
    volumeBindings[0].descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
    volumeBindings[0].descriptorCount = 1;
    volumeBindings[0].stageFlags = VK_SHADER_STAGE_VERTEX_BIT | VK_SHADER_STAGE_FRAGMENT_BIT;
    volumeBindings[1].binding = 1;
    volumeBindings[1].descriptorType = VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER;
    volumeBindings[1].descriptorCount = 1;
    volumeBindings[1].stageFlags = VK_SHADER_STAGE_FRAGMENT_BIT;
    VkDescriptorSetLayoutCreateInfo volumeLayoutInfo{
        VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO};
    volumeLayoutInfo.bindingCount = 2;
    volumeLayoutInfo.pBindings = volumeBindings;
    if (vkCreateDescriptorSetLayout(device.handle(), &volumeLayoutInfo, nullptr, &volumeLayout_) !=
        VK_SUCCESS) {
        error = "vkCreateDescriptorSetLayout failed for the volumes";
        destroy();
        return false;
    }

    VkDescriptorPoolSize poolSizes[2] = {{VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 2},
                                         {VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1}};
    VkDescriptorPoolCreateInfo poolInfo{VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO};
    poolInfo.maxSets = 2;
    poolInfo.poolSizeCount = 2;
    poolInfo.pPoolSizes = poolSizes;
    if (vkCreateDescriptorPool(device.handle(), &poolInfo, nullptr, &descriptorPool_) !=
        VK_SUCCESS) {
        error = "vkCreateDescriptorPool failed";
        destroy();
        return false;
    }

    VkDescriptorSetLayout wanted[2] = {materialLayout_, volumeLayout_};
    VkDescriptorSet allocated[2] = {VK_NULL_HANDLE, VK_NULL_HANDLE};
    VkDescriptorSetAllocateInfo setAllocation{VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO};
    setAllocation.descriptorPool = descriptorPool_;
    setAllocation.descriptorSetCount = 2;
    setAllocation.pSetLayouts = wanted;
    if (vkAllocateDescriptorSets(device.handle(), &setAllocation, allocated) != VK_SUCCESS) {
        error = "vkAllocateDescriptorSets failed";
        destroy();
        return false;
    }
    materialSet_ = allocated[0];
    volumeSet_ = allocated[1];

    // Nearest and unnormalised-in-spirit: the shader uses texelFetch, so no
    // filtering is ever applied. A linear filter here would silently average two
    // depths across a silhouette and put smoke a metre in front of a bulkhead.
    VkSamplerCreateInfo samplerInfo{VK_STRUCTURE_TYPE_SAMPLER_CREATE_INFO};
    samplerInfo.magFilter = VK_FILTER_NEAREST;
    samplerInfo.minFilter = VK_FILTER_NEAREST;
    samplerInfo.mipmapMode = VK_SAMPLER_MIPMAP_MODE_NEAREST;
    samplerInfo.addressModeU = VK_SAMPLER_ADDRESS_MODE_CLAMP_TO_EDGE;
    samplerInfo.addressModeV = VK_SAMPLER_ADDRESS_MODE_CLAMP_TO_EDGE;
    samplerInfo.addressModeW = VK_SAMPLER_ADDRESS_MODE_CLAMP_TO_EDGE;
    samplerInfo.borderColor = VK_BORDER_COLOR_FLOAT_OPAQUE_WHITE;
    if (vkCreateSampler(device.handle(), &samplerInfo, nullptr, &depthSampler_) != VK_SUCCESS) {
        error = "vkCreateSampler failed";
        destroy();
        return false;
    }

    VkDescriptorImageInfo depthInfo{depthSampler_, depth_.view,
                                    VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL};
    VkWriteDescriptorSet writeDepth{VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET};
    writeDepth.dstSet = volumeSet_;
    writeDepth.dstBinding = 1;
    writeDepth.descriptorCount = 1;
    writeDepth.descriptorType = VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER;
    writeDepth.pImageInfo = &depthInfo;
    vkUpdateDescriptorSets(device.handle(), 1, &writeDepth, 0, nullptr);

    // --- pipelines ------------------------------------------------------------
    VkViewport viewport{0.0f, 0.0f, static_cast<float>(width), static_cast<float>(height), 0.0f,
                        1.0f};
    VkRect2D scissor{{0, 0}, {width, height}};
    VkPipelineViewportStateCreateInfo viewportState{
        VK_STRUCTURE_TYPE_PIPELINE_VIEWPORT_STATE_CREATE_INFO};
    viewportState.viewportCount = 1;
    viewportState.pViewports = &viewport;
    viewportState.scissorCount = 1;
    viewportState.pScissors = &scissor;

    VkPipelineMultisampleStateCreateInfo multisample{
        VK_STRUCTURE_TYPE_PIPELINE_MULTISAMPLE_STATE_CREATE_INFO};
    multisample.rasterizationSamples = VK_SAMPLE_COUNT_1_BIT;

    VkPipelineInputAssemblyStateCreateInfo assembly{
        VK_STRUCTURE_TYPE_PIPELINE_INPUT_ASSEMBLY_STATE_CREATE_INFO};
    assembly.topology = VK_PRIMITIVE_TOPOLOGY_TRIANGLE_LIST;

    // The lit pipeline, state for state as HullRenderer builds it. Anything that
    // differs here shows up as a frame that is *nearly* the hull renderer's, which
    // is the failure the bit-identity check exists to catch.
    {
        VkShaderModule vertexModule = device.loadShader(shaderDirectory + "/hull.vert.spv", error);
        if (vertexModule == nullptr) { destroy(); return false; }
        VkShaderModule fragmentModule =
            device.loadShader(shaderDirectory + "/hull.frag.spv", error);
        if (fragmentModule == nullptr) {
            device.destroyShader(vertexModule);
            destroy();
            return false;
        }

        VkPushConstantRange push{VK_SHADER_STAGE_VERTEX_BIT | VK_SHADER_STAGE_FRAGMENT_BIT, 0,
                                 sizeof(OpaquePush)};
        VkPipelineLayoutCreateInfo layoutInfo{VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO};
        layoutInfo.setLayoutCount = 1;
        layoutInfo.pSetLayouts = &materialLayout_;
        layoutInfo.pushConstantRangeCount = 1;
        layoutInfo.pPushConstantRanges = &push;
        if (vkCreatePipelineLayout(device.handle(), &layoutInfo, nullptr, &opaqueLayout_) !=
            VK_SUCCESS) {
            error = "vkCreatePipelineLayout failed for the lit pass";
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

        VkPipelineRasterizationStateCreateInfo raster{
            VK_STRUCTURE_TYPE_PIPELINE_RASTERIZATION_STATE_CREATE_INFO};
        raster.polygonMode = VK_POLYGON_MODE_FILL;
        raster.cullMode = VK_CULL_MODE_NONE;
        raster.frontFace = VK_FRONT_FACE_COUNTER_CLOCKWISE;
        raster.lineWidth = 1.0f;

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

        VkGraphicsPipelineCreateInfo pipelineInfo{
            VK_STRUCTURE_TYPE_GRAPHICS_PIPELINE_CREATE_INFO};
        pipelineInfo.stageCount = 2;
        pipelineInfo.pStages = stages;
        pipelineInfo.pVertexInputState = &vertexInput;
        pipelineInfo.pInputAssemblyState = &assembly;
        pipelineInfo.pViewportState = &viewportState;
        pipelineInfo.pRasterizationState = &raster;
        pipelineInfo.pMultisampleState = &multisample;
        pipelineInfo.pDepthStencilState = &depthStencil;
        pipelineInfo.pColorBlendState = &blend;
        pipelineInfo.layout = opaqueLayout_;
        pipelineInfo.renderPass = opaquePass_;

        const VkResult created = vkCreateGraphicsPipelines(device.handle(), VK_NULL_HANDLE, 1,
                                                           &pipelineInfo, nullptr,
                                                           &opaquePipeline_);
        device.destroyShader(vertexModule);
        device.destroyShader(fragmentModule);
        if (created != VK_SUCCESS) {
            error = "vkCreateGraphicsPipelines failed for the lit pass";
            destroy();
            return false;
        }
    }

    // The volumetric pipeline.
    {
        VkShaderModule vertexModule = device.loadShader(shaderDirectory + "/smoke.vert.spv", error);
        if (vertexModule == nullptr) { destroy(); return false; }
        VkShaderModule fragmentModule =
            device.loadShader(shaderDirectory + "/smoke.frag.spv", error);
        if (fragmentModule == nullptr) {
            device.destroyShader(vertexModule);
            destroy();
            return false;
        }

        VkPushConstantRange push{VK_SHADER_STAGE_VERTEX_BIT | VK_SHADER_STAGE_FRAGMENT_BIT, 0,
                                 sizeof(SmokePush)};
        VkPipelineLayoutCreateInfo layoutInfo{VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO};
        layoutInfo.setLayoutCount = 1;
        layoutInfo.pSetLayouts = &volumeLayout_;
        layoutInfo.pushConstantRangeCount = 1;
        layoutInfo.pPushConstantRanges = &push;
        if (vkCreatePipelineLayout(device.handle(), &layoutInfo, nullptr, &smokeLayout_) !=
            VK_SUCCESS) {
            error = "vkCreatePipelineLayout failed for the volumetric pass";
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

        // No vertex buffer at all: the box comes out of gl_VertexIndex.
        VkPipelineVertexInputStateCreateInfo vertexInput{
            VK_STRUCTURE_TYPE_PIPELINE_VERTEX_INPUT_STATE_CREATE_INFO};

        VkPipelineRasterizationStateCreateInfo raster{
            VK_STRUCTURE_TYPE_PIPELINE_RASTERIZATION_STATE_CREATE_INFO};
        raster.polygonMode = VK_POLYGON_MODE_FILL;
        // **Keep the far faces, drop the near ones**, so each covered pixel gets
        // exactly one fragment and the medium is applied once.
        //
        // The sense of this was settled by a test and not by reasoning about the
        // y flip, because **from outside the box both senses look identical**: the
        // shading depends only on the ray, so keeping the near faces instead gives
        // the same one fragment per pixel and the same answer. The two differ in
        // exactly one place -- with the camera *inside* the medium, where the near
        // faces are behind the eye and keeping them draws nothing at all. The
        // first version of this pipeline had it the other way round, every
        // closed-form check from outside passed at one least-significant bit, and
        // `testTheCameraCanBeInsideTheMedium` was off by 152.
        raster.cullMode = VK_CULL_MODE_FRONT_BIT;
        raster.frontFace = VK_FRONT_FACE_COUNTER_CLOCKWISE;
        raster.lineWidth = 1.0f;

        // No depth attachment in this pass, and no test: where the hull is nearer
        // than the box, the ray has to *stop* at it, not vanish.
        VkPipelineDepthStencilStateCreateInfo depthStencil{
            VK_STRUCTURE_TYPE_PIPELINE_DEPTH_STENCIL_STATE_CREATE_INFO};

        // Premultiplied alpha: dst is multiplied by exactly the transmittance the
        // shader computed, so Beer-Lambert is the blender's arithmetic.
        VkPipelineColorBlendAttachmentState blendAttachment{};
        blendAttachment.blendEnable = VK_TRUE;
        blendAttachment.srcColorBlendFactor = VK_BLEND_FACTOR_ONE;
        blendAttachment.dstColorBlendFactor = VK_BLEND_FACTOR_ONE_MINUS_SRC_ALPHA;
        blendAttachment.colorBlendOp = VK_BLEND_OP_ADD;
        blendAttachment.srcAlphaBlendFactor = VK_BLEND_FACTOR_ONE;
        blendAttachment.dstAlphaBlendFactor = VK_BLEND_FACTOR_ONE_MINUS_SRC_ALPHA;
        blendAttachment.alphaBlendOp = VK_BLEND_OP_ADD;
        blendAttachment.colorWriteMask = VK_COLOR_COMPONENT_R_BIT | VK_COLOR_COMPONENT_G_BIT |
                                         VK_COLOR_COMPONENT_B_BIT | VK_COLOR_COMPONENT_A_BIT;
        VkPipelineColorBlendStateCreateInfo blend{
            VK_STRUCTURE_TYPE_PIPELINE_COLOR_BLEND_STATE_CREATE_INFO};
        blend.attachmentCount = 1;
        blend.pAttachments = &blendAttachment;

        VkGraphicsPipelineCreateInfo pipelineInfo{
            VK_STRUCTURE_TYPE_GRAPHICS_PIPELINE_CREATE_INFO};
        pipelineInfo.stageCount = 2;
        pipelineInfo.pStages = stages;
        pipelineInfo.pVertexInputState = &vertexInput;
        pipelineInfo.pInputAssemblyState = &assembly;
        pipelineInfo.pViewportState = &viewportState;
        pipelineInfo.pRasterizationState = &raster;
        pipelineInfo.pMultisampleState = &multisample;
        pipelineInfo.pDepthStencilState = &depthStencil;
        pipelineInfo.pColorBlendState = &blend;
        pipelineInfo.layout = smokeLayout_;
        pipelineInfo.renderPass = smokePass_;

        const VkResult created = vkCreateGraphicsPipelines(device.handle(), VK_NULL_HANDLE, 1,
                                                           &pipelineInfo, nullptr,
                                                           &smokePipeline_);
        device.destroyShader(vertexModule);
        device.destroyShader(fragmentModule);
        if (created != VK_SUCCESS) {
            error = "vkCreateGraphicsPipelines failed for the volumetric pass";
            destroy();
            return false;
        }
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

    // Four timestamps, so the volumetric pass can be priced separately from the
    // solid it sits on. A failure here is not fatal; the figures stay zero rather
    // than becoming numbers nobody measured.
    VkQueryPoolCreateInfo queryInfo{VK_STRUCTURE_TYPE_QUERY_POOL_CREATE_INFO};
    queryInfo.queryType = VK_QUERY_TYPE_TIMESTAMP;
    queryInfo.queryCount = 4;
    if (vkCreateQueryPool(device.handle(), &queryInfo, nullptr, &queryPool_) != VK_SUCCESS)
        queryPool_ = nullptr;
    return true;
}

void SmokeRenderer::destroy() {
    if (device_ == nullptr) return;
    VkDevice handle = device_->handle();
    if (handle != nullptr) {
        vkDeviceWaitIdle(handle);
        if (queryPool_ != nullptr) vkDestroyQueryPool(handle, queryPool_, nullptr);
        if (smokePipeline_ != nullptr) vkDestroyPipeline(handle, smokePipeline_, nullptr);
        if (smokeLayout_ != nullptr) vkDestroyPipelineLayout(handle, smokeLayout_, nullptr);
        if (opaquePipeline_ != nullptr) vkDestroyPipeline(handle, opaquePipeline_, nullptr);
        if (opaqueLayout_ != nullptr) vkDestroyPipelineLayout(handle, opaqueLayout_, nullptr);
        if (depthSampler_ != nullptr) vkDestroySampler(handle, depthSampler_, nullptr);
        if (descriptorPool_ != nullptr) vkDestroyDescriptorPool(handle, descriptorPool_, nullptr);
        if (volumeLayout_ != nullptr) vkDestroyDescriptorSetLayout(handle, volumeLayout_, nullptr);
        if (materialLayout_ != nullptr)
            vkDestroyDescriptorSetLayout(handle, materialLayout_, nullptr);
        if (smokeFramebuffer_ != nullptr) vkDestroyFramebuffer(handle, smokeFramebuffer_, nullptr);
        if (opaqueFramebuffer_ != nullptr)
            vkDestroyFramebuffer(handle, opaqueFramebuffer_, nullptr);
        if (smokePass_ != nullptr) vkDestroyRenderPass(handle, smokePass_, nullptr);
        if (opaquePass_ != nullptr) vkDestroyRenderPass(handle, opaquePass_, nullptr);
        device_->destroyBuffer(vertexBuffer_);
        device_->destroyBuffer(indexBuffer_);
        device_->destroyBuffer(materialBuffer_);
        device_->destroyBuffer(volumeBuffer_);
        device_->destroyBuffer(readback_);
        device_->destroyImage(colour_);
        device_->destroyImage(depth_);
    }
    queryPool_ = nullptr;
    smokePipeline_ = nullptr;
    smokeLayout_ = nullptr;
    opaquePipeline_ = nullptr;
    opaqueLayout_ = nullptr;
    depthSampler_ = nullptr;
    volumeSet_ = nullptr;
    materialSet_ = nullptr;
    descriptorPool_ = nullptr;
    volumeLayout_ = nullptr;
    materialLayout_ = nullptr;
    smokeFramebuffer_ = nullptr;
    opaqueFramebuffer_ = nullptr;
    smokePass_ = nullptr;
    opaquePass_ = nullptr;
    device_ = nullptr;
    width_ = height_ = 0;
    materialRevision_ = 0;
    materialCapacity_ = 0;
    volumeCapacity_ = 0;
}

bool SmokeRenderer::ensureGeometryCapacity(std::size_t vertexBytes, std::size_t indexBytes) {
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

bool SmokeRenderer::ensureMaterials(const MaterialLibrary& library) {
    const std::size_t count = library.packed().size();
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
        write.dstSet = materialSet_;
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

bool SmokeRenderer::ensureVolumes(const std::vector<SmokeVolume>& sorted) {
    // One row even for an empty list: the vertex stage indexes this buffer
    // unconditionally, and a descriptor that points at nothing is a validation
    // error whether or not anything is drawn.
    const std::size_t rows = std::max<std::size_t>(sorted.size(), 1);
    if (!volumeBuffer_.valid() || rows > volumeCapacity_) {
        vkDeviceWaitIdle(device_->handle());
        device_->destroyBuffer(volumeBuffer_);
        volumeBuffer_ = device_->createBuffer(
            rows * sizeof(GpuSmokeVolume),
            VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT,
            VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT);
        if (!volumeBuffer_.valid()) return false;
        volumeCapacity_ = rows;

        VkDescriptorBufferInfo bufferInfo{volumeBuffer_.handle, 0, VK_WHOLE_SIZE};
        VkWriteDescriptorSet write{VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET};
        write.dstSet = volumeSet_;
        write.dstBinding = 0;
        write.descriptorCount = 1;
        write.descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
        write.pBufferInfo = &bufferInfo;
        vkUpdateDescriptorSets(device_->handle(), 1, &write, 0, nullptr);
    }

    std::vector<GpuSmokeVolume> packed(rows, GpuSmokeVolume{});
    for (std::size_t i = 0; i < sorted.size(); ++i) packed[i] = pack(sorted[i]);
    return device_->upload(volumeBuffer_, packed.data(), rows * sizeof(GpuSmokeVolume));
}

bool SmokeRenderer::render(const float mvp[16], const SceneView& view, const SceneMesh& mesh,
                           const MaterialLibrary& library,
                           const std::vector<SmokeVolume>& volumes, const float clearColour[4],
                           core::Image& out) {
    if (!valid()) return false;
    DepthBasis basis;
    // Refused rather than approximated: without a centre of projection there is no
    // way to turn a depth-buffer value back into a distance, and a silently wrong
    // distance would put the medium in front of or behind everything.
    if (!depthBasisFrom(mvp, basis)) return false;

    const auto frameStarted = std::chrono::steady_clock::now();
    lastFrame_ = SmokeFrameCost{};

    const std::size_t vertexCount = mesh.vertices().size();
    const std::size_t indexCount = mesh.indices().size();
    const bool hasGeometry = vertexCount > 0 && indexCount > 0;
    lastFrame_.vertices = vertexCount;
    lastFrame_.triangles = indexCount / 3;
    lastFrame_.volumes = volumes.size();

    // Back to front from the eye. Emission makes the composite order-dependent
    // even for media that never overlap, so this is arithmetic rather than a
    // transparency convention.
    const sim::Vec3 eye{view.eye[0], view.eye[1], view.eye[2]};
    sorted_ = volumes;
    std::stable_sort(sorted_.begin(), sorted_.end(),
                     [&eye](const SmokeVolume& a, const SmokeVolume& b) {
                         return sim::length2(a.centreWorld() - eye) >
                                sim::length2(b.centreWorld() - eye);
                     });

    const auto uploadStarted = std::chrono::steady_clock::now();
    if (!ensureMaterials(library)) return false;
    if (!ensureVolumes(sorted_)) return false;
    if (hasGeometry) {
        const std::size_t vertexBytes = vertexCount * sizeof(HullVertex);
        const std::size_t indexBytes = indexCount * sizeof(std::uint32_t);
        if (!ensureGeometryCapacity(vertexBytes, indexBytes)) return false;
        if (!device_->upload(vertexBuffer_, mesh.vertices().data(), vertexBytes)) return false;
        if (!device_->upload(indexBuffer_, mesh.indices().data(), indexBytes)) return false;
    }
    lastFrame_.uploadSeconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - uploadStarted).count();

    OpaquePush opaquePush{};
    std::memcpy(opaquePush.modelViewProjection, mvp, sizeof(opaquePush.modelViewProjection));
    const double length = std::sqrt(
        static_cast<double>(view.sunDirection[0]) * view.sunDirection[0] +
        static_cast<double>(view.sunDirection[1]) * view.sunDirection[1] +
        static_cast<double>(view.sunDirection[2]) * view.sunDirection[2]);
    const double inverse = length > 1e-12 ? 1.0 / length : 0.0;
    for (int i = 0; i < 3; ++i)
        opaquePush.sun[i] = static_cast<float>(view.sunDirection[i] * inverse);
    opaquePush.sun[3] = view.mode == HullShadingMode::MaterialId ? 1.0f
                        : view.mode == HullShadingMode::Depth    ? 2.0f
                                                                 : 0.0f;
    for (int i = 0; i < 3; ++i) opaquePush.sunColour[i] = view.sunColour[i];
    opaquePush.sunColour[3] = view.exposure;
    for (int i = 0; i < 3; ++i) opaquePush.sky[i] = view.skyColour[i];
    opaquePush.sky[3] = 0.0f;
    for (int i = 0; i < 3; ++i) opaquePush.eye[i] = view.eye[i];
    opaquePush.eye[3] = 0.0f;

    SmokePush smokePush{};
    std::memcpy(smokePush.modelViewProjection, mvp, sizeof(smokePush.modelViewProjection));
    for (int i = 0; i < 3; ++i) smokePush.eyeA[i] = view.eye[i];
    smokePush.eyeA[3] = static_cast<float>(basis.a);
    for (int i = 0; i < 3; ++i) smokePush.forwardB[i] = static_cast<float>(basis.forward[i]);
    smokePush.forwardB[3] = static_cast<float>(basis.b);

    VkClearValue clears[2]{};
    for (int i = 0; i < 4; ++i) clears[0].color.float32[i] = clearColour[i];
    clears[1].depthStencil.depth = 1.0f;

    const auto submitStarted = std::chrono::steady_clock::now();
    VkCommandBuffer commands = device_->beginOneShot();

    if (queryPool_ != nullptr) {
        vkCmdResetQueryPool(commands, queryPool_, 0, 4);
        vkCmdWriteTimestamp(commands, VK_PIPELINE_STAGE_TOP_OF_PIPE_BIT, queryPool_, 0);
    }

    VkRenderPassBeginInfo begin{VK_STRUCTURE_TYPE_RENDER_PASS_BEGIN_INFO};
    begin.renderPass = opaquePass_;
    begin.framebuffer = opaqueFramebuffer_;
    begin.renderArea = {{0, 0}, {width_, height_}};
    begin.clearValueCount = 2;
    begin.pClearValues = clears;
    vkCmdBeginRenderPass(commands, &begin, VK_SUBPASS_CONTENTS_INLINE);
    if (hasGeometry) {
        vkCmdBindPipeline(commands, VK_PIPELINE_BIND_POINT_GRAPHICS, opaquePipeline_);
        vkCmdBindDescriptorSets(commands, VK_PIPELINE_BIND_POINT_GRAPHICS, opaqueLayout_, 0, 1,
                                &materialSet_, 0, nullptr);
        vkCmdPushConstants(commands, opaqueLayout_,
                           VK_SHADER_STAGE_VERTEX_BIT | VK_SHADER_STAGE_FRAGMENT_BIT, 0,
                           sizeof(OpaquePush), &opaquePush);
        const VkDeviceSize offset = 0;
        vkCmdBindVertexBuffers(commands, 0, 1, &vertexBuffer_.handle, &offset);
        vkCmdBindIndexBuffer(commands, indexBuffer_.handle, 0, VK_INDEX_TYPE_UINT32);
        vkCmdDrawIndexed(commands, static_cast<std::uint32_t>(indexCount), 1, 0, 0, 0);
    }
    vkCmdEndRenderPass(commands);

    if (queryPool_ != nullptr)
        vkCmdWriteTimestamp(commands, VK_PIPELINE_STAGE_BOTTOM_OF_PIPE_BIT, queryPool_, 1);

    // Depth from written to read. The render pass left it in
    // DEPTH_STENCIL_ATTACHMENT_OPTIMAL; the fragment stage of the next pass wants
    // it as a texture.
    VkImageMemoryBarrier toRead{VK_STRUCTURE_TYPE_IMAGE_MEMORY_BARRIER};
    toRead.srcAccessMask = VK_ACCESS_DEPTH_STENCIL_ATTACHMENT_WRITE_BIT;
    toRead.dstAccessMask = VK_ACCESS_SHADER_READ_BIT;
    toRead.oldLayout = VK_IMAGE_LAYOUT_DEPTH_STENCIL_ATTACHMENT_OPTIMAL;
    toRead.newLayout = VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL;
    toRead.srcQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
    toRead.dstQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
    toRead.image = depth_.handle;
    toRead.subresourceRange = {VK_IMAGE_ASPECT_DEPTH_BIT, 0, 1, 0, 1};
    vkCmdPipelineBarrier(commands, VK_PIPELINE_STAGE_LATE_FRAGMENT_TESTS_BIT,
                         VK_PIPELINE_STAGE_FRAGMENT_SHADER_BIT, 0, 0, nullptr, 0, nullptr, 1,
                         &toRead);

    if (queryPool_ != nullptr)
        vkCmdWriteTimestamp(commands, VK_PIPELINE_STAGE_TOP_OF_PIPE_BIT, queryPool_, 2);

    begin.renderPass = smokePass_;
    begin.framebuffer = smokeFramebuffer_;
    begin.clearValueCount = 0;
    begin.pClearValues = nullptr;
    vkCmdBeginRenderPass(commands, &begin, VK_SUBPASS_CONTENTS_INLINE);
    if (!sorted_.empty()) {
        vkCmdBindPipeline(commands, VK_PIPELINE_BIND_POINT_GRAPHICS, smokePipeline_);
        vkCmdBindDescriptorSets(commands, VK_PIPELINE_BIND_POINT_GRAPHICS, smokeLayout_, 0, 1,
                                &volumeSet_, 0, nullptr);
        vkCmdPushConstants(commands, smokeLayout_,
                           VK_SHADER_STAGE_VERTEX_BIT | VK_SHADER_STAGE_FRAGMENT_BIT, 0,
                           sizeof(SmokePush), &smokePush);
        // One instance per volume, in the sorted order. Instances rasterize in
        // order, which is what makes the blend order the sort order.
        vkCmdDraw(commands, 36, static_cast<std::uint32_t>(sorted_.size()), 0, 0);
    }
    vkCmdEndRenderPass(commands);

    if (queryPool_ != nullptr)
        vkCmdWriteTimestamp(commands, VK_PIPELINE_STAGE_BOTTOM_OF_PIPE_BIT, queryPool_, 3);

    VkBufferImageCopy copy{};
    copy.imageSubresource.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
    copy.imageSubresource.layerCount = 1;
    copy.imageExtent = {width_, height_, 1};
    vkCmdCopyImageToBuffer(commands, colour_.handle, VK_IMAGE_LAYOUT_TRANSFER_SRC_OPTIMAL,
                           readback_.handle, 1, &copy);

    const bool submitted = device_->endOneShot(commands);
    lastFrame_.submitSeconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - submitStarted).count();
    // Timed first so the figure is recorded either way, then refused: everything
    // below reads back an image the device may never have drawn.
    if (!submitted) return false;

    if (queryPool_ != nullptr) {
        std::uint64_t stamps[4] = {0, 0, 0, 0};
        if (vkGetQueryPoolResults(device_->handle(), queryPool_, 0, 4, sizeof(stamps), stamps,
                                  sizeof(std::uint64_t),
                                  VK_QUERY_RESULT_64_BIT | VK_QUERY_RESULT_WAIT_BIT) ==
            VK_SUCCESS) {
            const double nanos = static_cast<double>(device_->timestampPeriodNanos()) * 1e-9;
            if (stamps[1] >= stamps[0])
                lastFrame_.opaqueGpuSeconds = static_cast<double>(stamps[1] - stamps[0]) * nanos;
            if (stamps[3] >= stamps[2])
                lastFrame_.smokeGpuSeconds = static_cast<double>(stamps[3] - stamps[2]) * nanos;
        }
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
