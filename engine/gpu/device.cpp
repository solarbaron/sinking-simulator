// SPDX-License-Identifier: MIT
#include "device.hpp"

#include <vulkan/vulkan.h>

#include <algorithm>
#include <cstring>
#include <fstream>

namespace gpu {
namespace {

std::vector<std::uint32_t> readSpirv(const std::string& path) {
    std::ifstream file(path, std::ios::binary | std::ios::ate);
    if (!file) return {};
    const auto bytes = static_cast<std::size_t>(file.tellg());
    if (bytes == 0 || bytes % 4 != 0) return {};
    file.seekg(0);
    std::vector<std::uint32_t> words(bytes / 4);
    file.read(reinterpret_cast<char*>(words.data()), static_cast<std::streamsize>(bytes));
    return words;
}

}  // namespace

Device::~Device() { destroy(); }

bool Device::create(std::string& error) {
    VkApplicationInfo app{VK_STRUCTURE_TYPE_APPLICATION_INFO};
    app.pApplicationName = "shipsim";
    app.apiVersion = VK_API_VERSION_1_3;

    VkInstanceCreateInfo instanceInfo{VK_STRUCTURE_TYPE_INSTANCE_CREATE_INFO};
    instanceInfo.pApplicationInfo = &app;
    if (vkCreateInstance(&instanceInfo, nullptr, &instance_) != VK_SUCCESS) {
        error = "vkCreateInstance failed (no Vulkan loader or driver?)";
        instance_ = nullptr;
        return false;
    }

    std::uint32_t count = 0;
    vkEnumeratePhysicalDevices(instance_, &count, nullptr);
    if (count == 0) {
        error = "no Vulkan physical devices";
        destroy();
        return false;
    }
    std::vector<VkPhysicalDevice> devices(count);
    vkEnumeratePhysicalDevices(instance_, &count, devices.data());

    physical_ = devices[0];
    for (VkPhysicalDevice candidate : devices) {
        VkPhysicalDeviceProperties properties;
        vkGetPhysicalDeviceProperties(candidate, &properties);
        if (properties.deviceType == VK_PHYSICAL_DEVICE_TYPE_DISCRETE_GPU) {
            physical_ = candidate;
            break;
        }
    }

    VkPhysicalDeviceProperties properties;
    vkGetPhysicalDeviceProperties(physical_, &properties);
    deviceName_ = properties.deviceName;
    timestampPeriod_ = properties.limits.timestampPeriod;

    std::uint32_t families = 0;
    vkGetPhysicalDeviceQueueFamilyProperties(physical_, &families, nullptr);
    std::vector<VkQueueFamilyProperties> familyProperties(families);
    vkGetPhysicalDeviceQueueFamilyProperties(physical_, &families, familyProperties.data());

    // Prefer a family that can do graphics as well as compute and transfer, so
    // one queue serves both the solver and the renderer. Fall back to
    // compute+transfer, which is enough for everything except drawing.
    bool found = false;
    constexpr VkQueueFlags kComputeTransfer = VK_QUEUE_COMPUTE_BIT | VK_QUEUE_TRANSFER_BIT;
    for (std::uint32_t pass = 0; pass < 2 && !found; ++pass) {
        const VkQueueFlags wanted =
            pass == 0 ? (kComputeTransfer | VK_QUEUE_GRAPHICS_BIT) : kComputeTransfer;
        for (std::uint32_t i = 0; i < families; ++i)
            if ((familyProperties[i].queueFlags & wanted) == wanted) {
                queueFamily_ = i;
                graphicsCapable_ = pass == 0;
                found = true;
                break;
            }
    }
    if (!found) {
        error = "no queue family supports compute and transfer";
        destroy();
        return false;
    }

    const float priority = 1.0f;
    VkDeviceQueueCreateInfo queueInfo{VK_STRUCTURE_TYPE_DEVICE_QUEUE_CREATE_INFO};
    queueInfo.queueFamilyIndex = queueFamily_;
    queueInfo.queueCount = 1;
    queueInfo.pQueuePriorities = &priority;

    VkDeviceCreateInfo deviceInfo{VK_STRUCTURE_TYPE_DEVICE_CREATE_INFO};
    deviceInfo.queueCreateInfoCount = 1;
    deviceInfo.pQueueCreateInfos = &queueInfo;
    if (vkCreateDevice(physical_, &deviceInfo, nullptr, &device_) != VK_SUCCESS) {
        error = "vkCreateDevice failed";
        device_ = nullptr;
        destroy();
        return false;
    }
    vkGetDeviceQueue(device_, queueFamily_, 0, &queue_);

    VkCommandPoolCreateInfo poolInfo{VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO};
    poolInfo.queueFamilyIndex = queueFamily_;
    poolInfo.flags = VK_COMMAND_POOL_CREATE_RESET_COMMAND_BUFFER_BIT;
    if (vkCreateCommandPool(device_, &poolInfo, nullptr, &commandPool_) != VK_SUCCESS) {
        error = "vkCreateCommandPool failed";
        destroy();
        return false;
    }
    return true;
}

void Device::destroy() {
    if (device_ != nullptr) {
        vkDeviceWaitIdle(device_);
        destroyBuffer(staging_);
        if (commandPool_ != nullptr) vkDestroyCommandPool(device_, commandPool_, nullptr);
        vkDestroyDevice(device_, nullptr);
    }
    if (instance_ != nullptr) vkDestroyInstance(instance_, nullptr);
    instance_ = nullptr;
    physical_ = nullptr;
    device_ = nullptr;
    queue_ = nullptr;
    commandPool_ = nullptr;
    graphicsCapable_ = false;
    deviceName_.clear();
}

std::uint32_t Device::findMemoryType(std::uint32_t typeBits, std::uint32_t want) const {
    VkPhysicalDeviceMemoryProperties properties;
    vkGetPhysicalDeviceMemoryProperties(physical_, &properties);
    for (std::uint32_t i = 0; i < properties.memoryTypeCount; ++i)
        if ((typeBits & (1u << i)) && (properties.memoryTypes[i].propertyFlags & want) == want)
            return i;
    return UINT32_MAX;
}

Device::Buffer Device::createBuffer(std::uint64_t bytes, std::uint32_t usage,
                                    std::uint32_t memoryFlags) {
    Buffer buffer;
    if (device_ == nullptr) return buffer;
    bytes = std::max<std::uint64_t>(bytes, 4);

    VkBufferCreateInfo info{VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO};
    info.size = bytes;
    info.usage = usage;
    info.sharingMode = VK_SHARING_MODE_EXCLUSIVE;
    if (vkCreateBuffer(device_, &info, nullptr, &buffer.handle) != VK_SUCCESS) return {};

    VkMemoryRequirements requirements;
    vkGetBufferMemoryRequirements(device_, buffer.handle, &requirements);
    const std::uint32_t type = findMemoryType(requirements.memoryTypeBits, memoryFlags);
    if (type == UINT32_MAX) {
        vkDestroyBuffer(device_, buffer.handle, nullptr);
        return {};
    }

    VkMemoryAllocateInfo allocation{VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO};
    allocation.allocationSize = requirements.size;
    allocation.memoryTypeIndex = type;
    if (vkAllocateMemory(device_, &allocation, nullptr, &buffer.memory) != VK_SUCCESS) {
        vkDestroyBuffer(device_, buffer.handle, nullptr);
        return {};
    }
    vkBindBufferMemory(device_, buffer.handle, buffer.memory, 0);
    buffer.size = bytes;
    return buffer;
}

void Device::destroyBuffer(Buffer& buffer) {
    if (device_ == nullptr) return;
    if (buffer.handle != nullptr) vkDestroyBuffer(device_, buffer.handle, nullptr);
    if (buffer.memory != nullptr) vkFreeMemory(device_, buffer.memory, nullptr);
    buffer = {};
}

bool Device::ensureStaging(std::uint64_t bytes) {
    if (staging_.valid() && staging_.size >= bytes) return true;
    destroyBuffer(staging_);
    staging_ = createBuffer(bytes,
                            VK_BUFFER_USAGE_TRANSFER_SRC_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT,
                            VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT |
                                VK_MEMORY_PROPERTY_HOST_COHERENT_BIT);
    return staging_.valid();
}

bool Device::upload(Buffer& target, const void* data, std::uint64_t bytes) {
    if (device_ == nullptr || !target.valid() || bytes == 0) return false;
    if (bytes > target.size || !ensureStaging(bytes)) return false;

    void* mapped = nullptr;
    if (vkMapMemory(device_, staging_.memory, 0, bytes, 0, &mapped) != VK_SUCCESS) return false;
    std::memcpy(mapped, data, bytes);
    vkUnmapMemory(device_, staging_.memory);

    VkCommandBuffer commands = beginOneShot();
    VkBufferCopy copy{0, 0, bytes};
    vkCmdCopyBuffer(commands, staging_.handle, target.handle, 1, &copy);
    endOneShot(commands);
    return true;
}

bool Device::download(const Buffer& source, void* data, std::uint64_t bytes) {
    if (device_ == nullptr || !source.valid() || bytes == 0) return false;
    if (bytes > source.size || !ensureStaging(bytes)) return false;

    VkCommandBuffer commands = beginOneShot();
    VkBufferCopy copy{0, 0, bytes};
    vkCmdCopyBuffer(commands, source.handle, staging_.handle, 1, &copy);
    endOneShot(commands);

    void* mapped = nullptr;
    if (vkMapMemory(device_, staging_.memory, 0, bytes, 0, &mapped) != VK_SUCCESS) return false;
    std::memcpy(data, mapped, bytes);
    vkUnmapMemory(device_, staging_.memory);
    return true;
}

Device::Image2D Device::createImage2D(std::uint32_t width, std::uint32_t height,
                                      std::uint32_t format, std::uint32_t usage,
                                      std::uint32_t aspect) {
    Image2D image;
    if (device_ == nullptr || width == 0 || height == 0) return image;

    VkImageCreateInfo info{VK_STRUCTURE_TYPE_IMAGE_CREATE_INFO};
    info.imageType = VK_IMAGE_TYPE_2D;
    info.format = static_cast<VkFormat>(format);
    info.extent = {width, height, 1};
    info.mipLevels = 1;
    info.arrayLayers = 1;
    info.samples = VK_SAMPLE_COUNT_1_BIT;
    info.tiling = VK_IMAGE_TILING_OPTIMAL;
    info.usage = usage;
    info.sharingMode = VK_SHARING_MODE_EXCLUSIVE;
    info.initialLayout = VK_IMAGE_LAYOUT_UNDEFINED;
    if (vkCreateImage(device_, &info, nullptr, &image.handle) != VK_SUCCESS) return {};

    VkMemoryRequirements requirements;
    vkGetImageMemoryRequirements(device_, image.handle, &requirements);
    const std::uint32_t type =
        findMemoryType(requirements.memoryTypeBits, VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT);
    if (type == UINT32_MAX) {
        vkDestroyImage(device_, image.handle, nullptr);
        return {};
    }

    VkMemoryAllocateInfo allocation{VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO};
    allocation.allocationSize = requirements.size;
    allocation.memoryTypeIndex = type;
    if (vkAllocateMemory(device_, &allocation, nullptr, &image.memory) != VK_SUCCESS) {
        vkDestroyImage(device_, image.handle, nullptr);
        return {};
    }
    vkBindImageMemory(device_, image.handle, image.memory, 0);

    VkImageViewCreateInfo viewInfo{VK_STRUCTURE_TYPE_IMAGE_VIEW_CREATE_INFO};
    viewInfo.image = image.handle;
    viewInfo.viewType = VK_IMAGE_VIEW_TYPE_2D;
    viewInfo.format = static_cast<VkFormat>(format);
    viewInfo.subresourceRange.aspectMask = aspect;
    viewInfo.subresourceRange.levelCount = 1;
    viewInfo.subresourceRange.layerCount = 1;
    if (vkCreateImageView(device_, &viewInfo, nullptr, &image.view) != VK_SUCCESS) {
        vkFreeMemory(device_, image.memory, nullptr);
        vkDestroyImage(device_, image.handle, nullptr);
        return {};
    }

    image.width = width;
    image.height = height;
    image.format = format;
    return image;
}

void Device::destroyImage(Image2D& image) {
    if (device_ == nullptr) return;
    if (image.view != nullptr) vkDestroyImageView(device_, image.view, nullptr);
    if (image.handle != nullptr) vkDestroyImage(device_, image.handle, nullptr);
    if (image.memory != nullptr) vkFreeMemory(device_, image.memory, nullptr);
    image = {};
}

VkShaderModule_T* Device::loadShader(const std::string& path, std::string& error) {
    const std::vector<std::uint32_t> spirv = readSpirv(path);
    if (spirv.empty()) {
        error = "could not read SPIR-V from " + path;
        return nullptr;
    }
    VkShaderModuleCreateInfo info{VK_STRUCTURE_TYPE_SHADER_MODULE_CREATE_INFO};
    info.codeSize = spirv.size() * 4;
    info.pCode = spirv.data();
    VkShaderModule module = VK_NULL_HANDLE;
    if (vkCreateShaderModule(device_, &info, nullptr, &module) != VK_SUCCESS) {
        error = "vkCreateShaderModule failed for " + path;
        return nullptr;
    }
    return module;
}

void Device::destroyShader(VkShaderModule_T* module) {
    if (device_ != nullptr && module != nullptr) vkDestroyShaderModule(device_, module, nullptr);
}

VkCommandBuffer_T* Device::beginOneShot() {
    VkCommandBufferAllocateInfo allocation{VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO};
    allocation.commandPool = commandPool_;
    allocation.level = VK_COMMAND_BUFFER_LEVEL_PRIMARY;
    allocation.commandBufferCount = 1;
    VkCommandBuffer commands = VK_NULL_HANDLE;
    vkAllocateCommandBuffers(device_, &allocation, &commands);

    VkCommandBufferBeginInfo begin{VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO};
    begin.flags = VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT;
    vkBeginCommandBuffer(commands, &begin);
    return commands;
}

void Device::endOneShot(VkCommandBuffer_T* commands) {
    vkEndCommandBuffer(commands);
    VkSubmitInfo submit{VK_STRUCTURE_TYPE_SUBMIT_INFO};
    submit.commandBufferCount = 1;
    submit.pCommandBuffers = &commands;
    vkQueueSubmit(queue_, 1, &submit, VK_NULL_HANDLE);
    vkQueueWaitIdle(queue_);
    vkFreeCommandBuffers(device_, commandPool_, 1, &commands);
}

}  // namespace gpu
