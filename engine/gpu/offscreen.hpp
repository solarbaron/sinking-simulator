// SPDX-License-Identifier: MIT
//
// Offscreen colour + depth render target and a minimal triangle pipeline.
//
// Rendering is verified without a window: draw here, read the colour attachment
// back into a core::Image, and assert on pixels. That keeps the visual checks
// runnable in CI, where a swapchain is not available and "it looked right" is
// not a test.
//
// Built on VkRenderPass and VkFramebuffer rather than dynamic rendering. Both
// work on the Pascal target, and docs/03-renderer-audio.md already notes that
// render passes still matter there; the explicit attachment layout transitions
// also make the readback path obvious rather than implicit.
#pragma once

#include "../core/png.hpp"
#include "device.hpp"

#include <cstdint>
#include <string>

struct VkRenderPass_T;
struct VkFramebuffer_T;
struct VkPipeline_T;
struct VkPipelineLayout_T;

namespace gpu {

// Position plus a flat colour. Deliberately the smallest vertex that can show a
// hull and its compartments apart; lighting and normals belong to Phase 2.
struct MeshVertex {
    float position[3];
    float colour[3];
};

class OffscreenRenderer {
public:
    ~OffscreenRenderer();

    bool create(Device& device, std::uint32_t width, std::uint32_t height,
                const std::string& shaderDirectory, std::string& error);
    void destroy();
    bool valid() const { return pipeline_ != nullptr; }

    std::uint32_t width() const { return width_; }
    std::uint32_t height() const { return height_; }

    struct Draw {
        const MeshVertex* vertices = nullptr;
        std::size_t vertexCount = 0;
        const std::uint32_t* indices = nullptr;
        std::size_t indexCount = 0;
    };

    // Clears to `clearColour`, draws the indexed triangles under `mvp` with depth
    // testing, and reads the colour attachment back into `out`. An empty draw is
    // legal and yields the clear colour.
    bool render(const float mvp[16], const Draw& draw, const float clearColour[4],
                core::Image& out);

private:
    bool ensureGeometryCapacity(std::size_t vertexBytes, std::size_t indexBytes);

    Device*             device_ = nullptr;
    std::uint32_t       width_ = 0;
    std::uint32_t       height_ = 0;

    Device::Image2D     colour_;
    Device::Image2D     depth_;
    VkRenderPass_T*     renderPass_ = nullptr;
    VkFramebuffer_T*    framebuffer_ = nullptr;
    VkPipelineLayout_T* pipelineLayout_ = nullptr;
    VkPipeline_T*       pipeline_ = nullptr;

    Device::Buffer      vertexBuffer_;
    Device::Buffer      indexBuffer_;
    Device::Buffer      readback_;
};

}  // namespace gpu
