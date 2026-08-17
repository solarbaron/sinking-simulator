// SPDX-License-Identifier: MIT
//
// Shared Vulkan context.
//
// This is the instance/device/queue/allocation code that was written for the FEM
// spike, generalised now that a second consumer exists. It graduated from
// `engine/gpu/fem_gpu.cpp` rather than being written fresh, because that code had
// already been run against a real GPU and had its timing verified -- see
// `docs/07-fem-spike-findings.md`.
//
// Deliberately small. It is not an abstraction layer over Vulkan and does not try
// to be: it owns the objects whose lifetime spans everything (instance, device,
// queue, command pool, staging buffer) and hands out raw handles for the rest.
// The render graph and bindless descriptor work belong a layer up.
//
// **Every entry point degrades rather than aborts when there is no GPU.** A
// headless CI box must still be able to run the rest of the test suite, so
// `create()` returns a failure reason instead of throwing or asserting.
#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

// Forward declarations, so callers that only pass a Device around do not pull in
// the whole Vulkan header.
struct VkInstance_T;
struct VkPhysicalDevice_T;
struct VkDevice_T;
struct VkQueue_T;
struct VkCommandPool_T;
struct VkBuffer_T;
struct VkDeviceMemory_T;
struct VkImage_T;
struct VkImageView_T;
struct VkCommandBuffer_T;
struct VkShaderModule_T;

namespace gpu {

class Device {
public:
    Device() = default;
    ~Device();
    Device(const Device&) = delete;
    Device& operator=(const Device&) = delete;

    // Returns false with a human-readable reason when no usable device exists.
    // Callers are expected to skip their work, not fail.
    bool create(std::string& error);
    void destroy();
    bool valid() const { return device_ != nullptr; }

    const std::string& deviceName() const { return deviceName_; }
    // Nanoseconds per timestamp tick, for turning query results into durations.
    float timestampPeriodNanos() const { return timestampPeriod_; }
    // True when the queue family supports graphics as well as compute. A
    // compute-only queue is enough for the FEM solver but not for rendering.
    bool supportsGraphics() const { return graphicsCapable_; }

    struct Buffer {
        VkBuffer_T* handle = nullptr;
        VkDeviceMemory_T* memory = nullptr;
        std::uint64_t size = 0;
        bool valid() const { return handle != nullptr; }
    };

    struct Image2D {
        VkImage_T* handle = nullptr;
        VkDeviceMemory_T* memory = nullptr;
        VkImageView_T* view = nullptr;
        std::uint32_t width = 0;
        std::uint32_t height = 0;
        std::uint32_t format = 0;  // VkFormat
        bool valid() const { return handle != nullptr; }
    };

    // `usage` and `memoryFlags` are VkBufferUsageFlags / VkMemoryPropertyFlags.
    // Passing them as plain integers keeps Vulkan out of this header.
    Buffer createBuffer(std::uint64_t bytes, std::uint32_t usage, std::uint32_t memoryFlags);
    void destroyBuffer(Buffer& buffer);
    // Both go through a shared staging buffer, so they are synchronous and are
    // meant for setup and readback, not for per-frame traffic.
    bool upload(Buffer& target, const void* data, std::uint64_t bytes);
    bool download(const Buffer& source, void* data, std::uint64_t bytes);

    Image2D createImage2D(std::uint32_t width, std::uint32_t height, std::uint32_t format,
                          std::uint32_t usage, std::uint32_t aspect);
    void destroyImage(Image2D& image);

    VkShaderModule_T* loadShader(const std::string& path, std::string& error);
    void destroyShader(VkShaderModule_T* module);

    // **Null on failure, and the submit's result comes back.** These used to
    // discard both: `beginOneShot` ignored `vkAllocateCommandBuffers` and then
    // recorded into whatever it got, and `endOneShot` dropped `vkQueueSubmit`, so
    // `upload` returned `true` for a transfer that never ran. A caller must check
    // the handle before recording and the result before believing the bytes moved.
    VkCommandBuffer_T* beginOneShot();
    [[nodiscard]] bool endOneShot(VkCommandBuffer_T* commands);

    // Raw handles, for the layers that do need Vulkan directly.
    VkInstance_T* instance() const { return instance_; }
    VkPhysicalDevice_T* physicalDevice() const { return physical_; }
    VkDevice_T* handle() const { return device_; }
    VkQueue_T* queue() const { return queue_; }
    VkCommandPool_T* commandPool() const { return commandPool_; }
    std::uint32_t queueFamily() const { return queueFamily_; }

    // Index of a memory type satisfying `typeBits` and `want`, or UINT32_MAX.
    std::uint32_t findMemoryType(std::uint32_t typeBits, std::uint32_t want) const;

private:
    bool ensureStaging(std::uint64_t bytes);

    VkInstance_T*       instance_ = nullptr;
    VkPhysicalDevice_T* physical_ = nullptr;
    VkDevice_T*         device_ = nullptr;
    VkQueue_T*          queue_ = nullptr;
    VkCommandPool_T*    commandPool_ = nullptr;
    std::uint32_t       queueFamily_ = 0;
    std::string         deviceName_;
    float               timestampPeriod_ = 1.0f;
    bool                graphicsCapable_ = false;

    Buffer staging_;
};

}  // namespace gpu
