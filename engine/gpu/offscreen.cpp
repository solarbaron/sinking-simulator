// SPDX-License-Identifier: MIT
#include "offscreen.hpp"

#include <vulkan/vulkan.h>

#include <cstring>
#include <vector>

namespace gpu {
namespace {

constexpr VkFormat kColourFormat = VK_FORMAT_R8G8B8A8_UNORM;
constexpr VkFormat kDepthFormat = VK_FORMAT_D32_SFLOAT;

}  // namespace

OffscreenRenderer::~OffscreenRenderer() { destroy(); }

bool OffscreenRenderer::create(Device& device, std::uint32_t width, std::uint32_t height,
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

    // The colour attachment ends the pass already in TRANSFER_SRC_OPTIMAL, so
    // the readback copy needs no extra barrier and the layout chain stays
    // visible in one place.
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

    VkShaderModule vertexModule = device.loadShader(shaderDirectory + "/mesh.vert.spv", error);
    if (vertexModule == nullptr) { destroy(); return false; }
    VkShaderModule fragmentModule = device.loadShader(shaderDirectory + "/mesh.frag.spv", error);
    if (fragmentModule == nullptr) {
        device.destroyShader(vertexModule);
        destroy();
        return false;
    }

    VkPushConstantRange push{VK_SHADER_STAGE_VERTEX_BIT, 0, 16 * sizeof(float)};
    VkPipelineLayoutCreateInfo layoutInfo{VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO};
    layoutInfo.pushConstantRangeCount = 1;
    layoutInfo.pPushConstantRanges = &push;
    vkCreatePipelineLayout(device.handle(), &layoutInfo, nullptr, &pipelineLayout_);

    VkPipelineShaderStageCreateInfo stages[2]{};
    stages[0].sType = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO;
    stages[0].stage = VK_SHADER_STAGE_VERTEX_BIT;
    stages[0].module = vertexModule;
    stages[0].pName = "main";
    stages[1].sType = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO;
    stages[1].stage = VK_SHADER_STAGE_FRAGMENT_BIT;
    stages[1].module = fragmentModule;
    stages[1].pName = "main";

    VkVertexInputBindingDescription binding{0, sizeof(MeshVertex), VK_VERTEX_INPUT_RATE_VERTEX};
    VkVertexInputAttributeDescription vertexAttributes[2]{};
    vertexAttributes[0] = {0, 0, VK_FORMAT_R32G32B32_SFLOAT, offsetof(MeshVertex, position)};
    vertexAttributes[1] = {1, 0, VK_FORMAT_R32G32B32_SFLOAT, offsetof(MeshVertex, colour)};

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
    // No culling: a debug view of a hull's interior compartments needs to show
    // back faces, and it keeps the tests independent of triangle winding.
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

void OffscreenRenderer::destroy() {
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

bool OffscreenRenderer::ensureGeometryCapacity(std::size_t vertexBytes, std::size_t indexBytes) {
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

bool OffscreenRenderer::render(const float mvp[16], const Draw& draw, const float clearColour[4],
                               core::Image& out) {
    if (!valid()) return false;

    const bool hasGeometry = draw.vertexCount > 0 && draw.indexCount > 0 &&
                             draw.vertices != nullptr && draw.indices != nullptr;
    if (hasGeometry) {
        const std::size_t vertexBytes = draw.vertexCount * sizeof(MeshVertex);
        const std::size_t indexBytes = draw.indexCount * sizeof(std::uint32_t);
        if (!ensureGeometryCapacity(vertexBytes, indexBytes)) return false;
        if (!device_->upload(vertexBuffer_, draw.vertices, vertexBytes)) return false;
        if (!device_->upload(indexBuffer_, draw.indices, indexBytes)) return false;
    }

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
        vkCmdPushConstants(commands, pipelineLayout_, VK_SHADER_STAGE_VERTEX_BIT, 0,
                           16 * sizeof(float), mvp);
        const VkDeviceSize offset = 0;
        vkCmdBindVertexBuffers(commands, 0, 1, &vertexBuffer_.handle, &offset);
        vkCmdBindIndexBuffer(commands, indexBuffer_.handle, 0, VK_INDEX_TYPE_UINT32);
        vkCmdDrawIndexed(commands, static_cast<std::uint32_t>(draw.indexCount), 1, 0, 0, 0);
    }
    vkCmdEndRenderPass(commands);

    // The render pass left the colour attachment in TRANSFER_SRC_OPTIMAL.
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
