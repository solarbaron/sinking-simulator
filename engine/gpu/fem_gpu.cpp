// SPDX-License-Identifier: MIT
#include "fem_gpu.hpp"

#include <vulkan/vulkan.h>

#include <cstring>
#include <fstream>
#include <vector>

namespace gpu {
namespace {

constexpr int kBufferCount = 10;
constexpr uint32_t kWorkgroupSize = 64;

struct PushConstants {
    float dt;
    float lambda;
    float mu;
    float gravity;
    float damping;
    uint32_t nodeCount;
    uint32_t tetCount;
};

std::vector<uint32_t> readSpirv(const std::string& path) {
    std::ifstream file(path, std::ios::binary | std::ios::ate);
    if (!file) return {};
    const auto bytes = static_cast<std::size_t>(file.tellg());
    file.seekg(0);
    std::vector<uint32_t> words(bytes / 4);
    file.read(reinterpret_cast<char*>(words.data()), static_cast<std::streamsize>(bytes));
    return words;
}

uint32_t groupsFor(std::size_t items) {
    return static_cast<uint32_t>((items + kWorkgroupSize - 1) / kWorkgroupSize);
}

}  // namespace

struct FemGpuSolver::Impl {
    VkInstance       instance = VK_NULL_HANDLE;
    VkPhysicalDevice physical = VK_NULL_HANDLE;
    VkDevice         device = VK_NULL_HANDLE;
    VkQueue          queue = VK_NULL_HANDLE;
    uint32_t         queueFamily = 0;
    VkCommandPool    commandPool = VK_NULL_HANDLE;
    VkQueryPool      queryPool = VK_NULL_HANDLE;
    float            timestampPeriod = 1.0f;

    VkDescriptorSetLayout setLayout = VK_NULL_HANDLE;
    VkDescriptorPool      descriptorPool = VK_NULL_HANDLE;
    VkDescriptorSet       descriptorSet = VK_NULL_HANDLE;
    VkPipelineLayout      pipelineLayout = VK_NULL_HANDLE;
    VkPipeline            forcePipeline = VK_NULL_HANDLE;
    VkPipeline            integratePipeline = VK_NULL_HANDLE;
    VkShaderModule        forceModule = VK_NULL_HANDLE;
    VkShaderModule        integrateModule = VK_NULL_HANDLE;

    struct Buffer {
        VkBuffer       handle = VK_NULL_HANDLE;
        VkDeviceMemory memory = VK_NULL_HANDLE;
        VkDeviceSize   size = 0;
    };
    Buffer buffers[kBufferCount];
    Buffer staging;

    std::size_t nodeCount = 0;
    std::size_t tetCount = 0;

    uint32_t findMemoryType(uint32_t typeBits, VkMemoryPropertyFlags want) const {
        VkPhysicalDeviceMemoryProperties props;
        vkGetPhysicalDeviceMemoryProperties(physical, &props);
        for (uint32_t i = 0; i < props.memoryTypeCount; ++i)
            if ((typeBits & (1u << i)) && (props.memoryTypes[i].propertyFlags & want) == want)
                return i;
        return UINT32_MAX;
    }

    bool createBuffer(Buffer& out, VkDeviceSize size, VkBufferUsageFlags usage,
                      VkMemoryPropertyFlags memoryFlags) {
        size = std::max<VkDeviceSize>(size, 4);
        VkBufferCreateInfo info{VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO};
        info.size = size;
        info.usage = usage;
        info.sharingMode = VK_SHARING_MODE_EXCLUSIVE;
        if (vkCreateBuffer(device, &info, nullptr, &out.handle) != VK_SUCCESS) return false;

        VkMemoryRequirements req;
        vkGetBufferMemoryRequirements(device, out.handle, &req);
        const uint32_t type = findMemoryType(req.memoryTypeBits, memoryFlags);
        // Unwound on both post-create paths, as `Device::createBuffer` does. These
        // two copies of it dropped the buffer on the floor: `out` is a member, so
        // the destructor reclaims it -- but it leaves a half-built `Buffer` with a
        // live handle, null memory and zero size, which `valid()` reads as usable.
        if (type == UINT32_MAX) {
            vkDestroyBuffer(device, out.handle, nullptr);
            out.handle = VK_NULL_HANDLE;
            return false;
        }

        VkMemoryAllocateInfo alloc{VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO};
        alloc.allocationSize = req.size;
        alloc.memoryTypeIndex = type;
        if (vkAllocateMemory(device, &alloc, nullptr, &out.memory) != VK_SUCCESS) {
            vkDestroyBuffer(device, out.handle, nullptr);
            out.handle = VK_NULL_HANDLE;
            return false;
        }
        if (vkBindBufferMemory(device, out.handle, out.memory, 0) != VK_SUCCESS) {
            vkFreeMemory(device, out.memory, nullptr);
            vkDestroyBuffer(device, out.handle, nullptr);
            out.memory = VK_NULL_HANDLE;
            out.handle = VK_NULL_HANDLE;
            return false;
        }
        out.size = size;
        return true;
    }

    VkCommandBuffer beginOneShot() {
        VkCommandBufferAllocateInfo alloc{VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO};
        alloc.commandPool = commandPool;
        alloc.level = VK_COMMAND_BUFFER_LEVEL_PRIMARY;
        alloc.commandBufferCount = 1;
        VkCommandBuffer cmd;
        vkAllocateCommandBuffers(device, &alloc, &cmd);
        VkCommandBufferBeginInfo begin{VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO};
        begin.flags = VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT;
        vkBeginCommandBuffer(cmd, &begin);
        return cmd;
    }

    void endOneShot(VkCommandBuffer cmd) {
        vkEndCommandBuffer(cmd);
        VkSubmitInfo submit{VK_STRUCTURE_TYPE_SUBMIT_INFO};
        submit.commandBufferCount = 1;
        submit.pCommandBuffers = &cmd;
        vkQueueSubmit(queue, 1, &submit, VK_NULL_HANDLE);
        vkQueueWaitIdle(queue);
        vkFreeCommandBuffers(device, commandPool, 1, &cmd);
    }

    // Host-visible staging is used for both directions; these transfers happen
    // once at setup and once at readback, so their cost is irrelevant.
    void uploadTo(Buffer& target, const void* data, VkDeviceSize bytes) {
        if (bytes == 0) return;
        void* mapped = nullptr;
        vkMapMemory(device, staging.memory, 0, bytes, 0, &mapped);
        std::memcpy(mapped, data, bytes);
        vkUnmapMemory(device, staging.memory);

        VkCommandBuffer cmd = beginOneShot();
        VkBufferCopy copy{0, 0, bytes};
        vkCmdCopyBuffer(cmd, staging.handle, target.handle, 1, &copy);
        endOneShot(cmd);
    }

    void downloadFrom(const Buffer& source, void* data, VkDeviceSize bytes) {
        if (bytes == 0) return;
        VkCommandBuffer cmd = beginOneShot();
        VkBufferCopy copy{0, 0, bytes};
        vkCmdCopyBuffer(cmd, source.handle, staging.handle, 1, &copy);
        endOneShot(cmd);

        void* mapped = nullptr;
        vkMapMemory(device, staging.memory, 0, bytes, 0, &mapped);
        std::memcpy(data, mapped, bytes);
        vkUnmapMemory(device, staging.memory);
    }
};

FemGpuSolver::~FemGpuSolver() {
    if (!impl_) return;
    Impl& d = *impl_;
    if (d.device) {
        vkDeviceWaitIdle(d.device);
        for (auto& b : d.buffers) {
            if (b.handle) vkDestroyBuffer(d.device, b.handle, nullptr);
            if (b.memory) vkFreeMemory(d.device, b.memory, nullptr);
        }
        if (d.staging.handle) vkDestroyBuffer(d.device, d.staging.handle, nullptr);
        if (d.staging.memory) vkFreeMemory(d.device, d.staging.memory, nullptr);
        if (d.forcePipeline) vkDestroyPipeline(d.device, d.forcePipeline, nullptr);
        if (d.integratePipeline) vkDestroyPipeline(d.device, d.integratePipeline, nullptr);
        if (d.forceModule) vkDestroyShaderModule(d.device, d.forceModule, nullptr);
        if (d.integrateModule) vkDestroyShaderModule(d.device, d.integrateModule, nullptr);
        if (d.pipelineLayout) vkDestroyPipelineLayout(d.device, d.pipelineLayout, nullptr);
        if (d.descriptorPool) vkDestroyDescriptorPool(d.device, d.descriptorPool, nullptr);
        if (d.setLayout) vkDestroyDescriptorSetLayout(d.device, d.setLayout, nullptr);
        if (d.queryPool) vkDestroyQueryPool(d.device, d.queryPool, nullptr);
        if (d.commandPool) vkDestroyCommandPool(d.device, d.commandPool, nullptr);
        vkDestroyDevice(d.device, nullptr);
    }
    if (d.instance) vkDestroyInstance(d.instance, nullptr);
    delete impl_;
    impl_ = nullptr;
}

bool FemGpuSolver::initialise(const sim::fem::TetMesh& mesh, const std::string& shaderDirectory,
                              std::string& error) {
    impl_ = new Impl();
    Impl& d = *impl_;
    d.nodeCount = mesh.nodeCount();
    d.tetCount = mesh.tetCount();

    VkApplicationInfo app{VK_STRUCTURE_TYPE_APPLICATION_INFO};
    app.pApplicationName = "shipsim-fem-spike";
    app.apiVersion = VK_API_VERSION_1_3;
    VkInstanceCreateInfo instanceInfo{VK_STRUCTURE_TYPE_INSTANCE_CREATE_INFO};
    instanceInfo.pApplicationInfo = &app;
    if (vkCreateInstance(&instanceInfo, nullptr, &d.instance) != VK_SUCCESS) {
        error = "vkCreateInstance failed (no Vulkan loader or driver?)";
        return false;
    }

    uint32_t count = 0;
    vkEnumeratePhysicalDevices(d.instance, &count, nullptr);
    if (count == 0) { error = "no Vulkan physical devices"; return false; }
    std::vector<VkPhysicalDevice> devices(count);
    vkEnumeratePhysicalDevices(d.instance, &count, devices.data());

    // Prefer a discrete GPU; fall back to whatever exists.
    d.physical = devices[0];
    for (VkPhysicalDevice candidate : devices) {
        VkPhysicalDeviceProperties props;
        vkGetPhysicalDeviceProperties(candidate, &props);
        if (props.deviceType == VK_PHYSICAL_DEVICE_TYPE_DISCRETE_GPU) { d.physical = candidate; break; }
    }
    VkPhysicalDeviceProperties props;
    vkGetPhysicalDeviceProperties(d.physical, &props);
    deviceName_ = props.deviceName;
    d.timestampPeriod = props.limits.timestampPeriod;

    uint32_t families = 0;
    vkGetPhysicalDeviceQueueFamilyProperties(d.physical, &families, nullptr);
    std::vector<VkQueueFamilyProperties> familyProps(families);
    vkGetPhysicalDeviceQueueFamilyProperties(d.physical, &families, familyProps.data());
    bool found = false;
    for (uint32_t i = 0; i < families; ++i)
        if ((familyProps[i].queueFlags & VK_QUEUE_COMPUTE_BIT) &&
            (familyProps[i].queueFlags & VK_QUEUE_TRANSFER_BIT)) {
            d.queueFamily = i;
            found = true;
            break;
        }
    if (!found) { error = "no compute-capable queue family"; return false; }

    const float priority = 1.0f;
    VkDeviceQueueCreateInfo queueInfo{VK_STRUCTURE_TYPE_DEVICE_QUEUE_CREATE_INFO};
    queueInfo.queueFamilyIndex = d.queueFamily;
    queueInfo.queueCount = 1;
    queueInfo.pQueuePriorities = &priority;
    VkDeviceCreateInfo deviceInfo{VK_STRUCTURE_TYPE_DEVICE_CREATE_INFO};
    deviceInfo.queueCreateInfoCount = 1;
    deviceInfo.pQueueCreateInfos = &queueInfo;
    if (vkCreateDevice(d.physical, &deviceInfo, nullptr, &d.device) != VK_SUCCESS) {
        error = "vkCreateDevice failed";
        return false;
    }
    device_ = d.device;
    vkGetDeviceQueue(d.device, d.queueFamily, 0, &d.queue);

    VkCommandPoolCreateInfo poolInfo{VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO};
    poolInfo.queueFamilyIndex = d.queueFamily;
    poolInfo.flags = VK_COMMAND_POOL_CREATE_RESET_COMMAND_BUFFER_BIT;
    vkCreateCommandPool(d.device, &poolInfo, nullptr, &d.commandPool);

    VkQueryPoolCreateInfo queryInfo{VK_STRUCTURE_TYPE_QUERY_POOL_CREATE_INFO};
    queryInfo.queryType = VK_QUERY_TYPE_TIMESTAMP;
    queryInfo.queryCount = 2;
    vkCreateQueryPool(d.device, &queryInfo, nullptr, &d.queryPool);

    // --- Buffers, in binding order ---
    const VkDeviceSize sizes[kBufferCount] = {
        mesh.position.size() * sizeof(float),          // 0 positions
        mesh.velocity.size() * sizeof(float),          // 1 velocities
        mesh.mass.size() * sizeof(float),              // 2 mass
        mesh.fixed.size() * sizeof(uint32_t),          // 3 fixed
        mesh.index.size() * sizeof(uint32_t),          // 4 indices
        mesh.restInv.size() * sizeof(float),           // 5 rest inverse
        mesh.restVolume.size() * sizeof(float),        // 6 rest volume
        mesh.tetForce.size() * sizeof(float),          // 7 element forces
        mesh.adjacencyOffset.size() * sizeof(uint32_t),// 8 CSR offsets
        mesh.adjacencyEntry.size() * sizeof(uint32_t), // 9 CSR entries
    };
    const void* data[kBufferCount] = {
        mesh.position.data(),  mesh.velocity.data(),   mesh.mass.data(),
        mesh.fixed.data(),     mesh.index.data(),      mesh.restInv.data(),
        mesh.restVolume.data(), mesh.tetForce.data(),  mesh.adjacencyOffset.data(),
        mesh.adjacencyEntry.data(),
    };

    VkDeviceSize largest = 0;
    for (VkDeviceSize s : sizes) largest = std::max(largest, s);
    if (!d.createBuffer(d.staging, largest,
                        VK_BUFFER_USAGE_TRANSFER_SRC_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT,
                        VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT)) {
        error = "staging buffer allocation failed";
        return false;
    }
    for (int i = 0; i < kBufferCount; ++i) {
        if (!d.createBuffer(d.buffers[i], sizes[i],
                            VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_SRC_BIT |
                                VK_BUFFER_USAGE_TRANSFER_DST_BIT,
                            VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT)) {
            error = "device buffer allocation failed";
            return false;
        }
        d.uploadTo(d.buffers[i], data[i], sizes[i]);
    }

    // --- Descriptors ---
    VkDescriptorSetLayoutBinding bindings[kBufferCount]{};
    for (int i = 0; i < kBufferCount; ++i) {
        bindings[i].binding = static_cast<uint32_t>(i);
        bindings[i].descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
        bindings[i].descriptorCount = 1;
        bindings[i].stageFlags = VK_SHADER_STAGE_COMPUTE_BIT;
    }
    VkDescriptorSetLayoutCreateInfo layoutInfo{VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO};
    layoutInfo.bindingCount = kBufferCount;
    layoutInfo.pBindings = bindings;
    vkCreateDescriptorSetLayout(d.device, &layoutInfo, nullptr, &d.setLayout);

    VkDescriptorPoolSize poolSize{VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, kBufferCount};
    VkDescriptorPoolCreateInfo descriptorPoolInfo{VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO};
    descriptorPoolInfo.maxSets = 1;
    descriptorPoolInfo.poolSizeCount = 1;
    descriptorPoolInfo.pPoolSizes = &poolSize;
    vkCreateDescriptorPool(d.device, &descriptorPoolInfo, nullptr, &d.descriptorPool);

    VkDescriptorSetAllocateInfo setAlloc{VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO};
    setAlloc.descriptorPool = d.descriptorPool;
    setAlloc.descriptorSetCount = 1;
    setAlloc.pSetLayouts = &d.setLayout;
    vkAllocateDescriptorSets(d.device, &setAlloc, &d.descriptorSet);

    VkDescriptorBufferInfo bufferInfos[kBufferCount]{};
    VkWriteDescriptorSet writes[kBufferCount]{};
    for (int i = 0; i < kBufferCount; ++i) {
        bufferInfos[i] = {d.buffers[i].handle, 0, VK_WHOLE_SIZE};
        writes[i].sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
        writes[i].dstSet = d.descriptorSet;
        writes[i].dstBinding = static_cast<uint32_t>(i);
        writes[i].descriptorCount = 1;
        writes[i].descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
        writes[i].pBufferInfo = &bufferInfos[i];
    }
    vkUpdateDescriptorSets(d.device, kBufferCount, writes, 0, nullptr);

    // --- Pipelines ---
    VkPushConstantRange pushRange{VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(PushConstants)};
    VkPipelineLayoutCreateInfo pipelineLayoutInfo{VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO};
    pipelineLayoutInfo.setLayoutCount = 1;
    pipelineLayoutInfo.pSetLayouts = &d.setLayout;
    pipelineLayoutInfo.pushConstantRangeCount = 1;
    pipelineLayoutInfo.pPushConstantRanges = &pushRange;
    vkCreatePipelineLayout(d.device, &pipelineLayoutInfo, nullptr, &d.pipelineLayout);

    auto makePipeline = [&](const std::string& file, VkShaderModule& module,
                            VkPipeline& pipeline) -> bool {
        const std::vector<uint32_t> spirv = readSpirv(shaderDirectory + "/" + file);
        if (spirv.empty()) { error = "could not read shader " + file; return false; }
        VkShaderModuleCreateInfo moduleInfo{VK_STRUCTURE_TYPE_SHADER_MODULE_CREATE_INFO};
        moduleInfo.codeSize = spirv.size() * 4;
        moduleInfo.pCode = spirv.data();
        if (vkCreateShaderModule(d.device, &moduleInfo, nullptr, &module) != VK_SUCCESS) {
            error = "vkCreateShaderModule failed for " + file;
            return false;
        }
        VkComputePipelineCreateInfo info{VK_STRUCTURE_TYPE_COMPUTE_PIPELINE_CREATE_INFO};
        info.stage.sType = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO;
        info.stage.stage = VK_SHADER_STAGE_COMPUTE_BIT;
        info.stage.module = module;
        info.stage.pName = "main";
        info.layout = d.pipelineLayout;
        if (vkCreateComputePipelines(d.device, VK_NULL_HANDLE, 1, &info, nullptr, &pipeline) !=
            VK_SUCCESS) {
            error = "vkCreateComputePipelines failed for " + file;
            return false;
        }
        return true;
    };

    if (!makePipeline("tet_forces.comp.spv", d.forceModule, d.forcePipeline)) return false;
    if (!makePipeline("node_integrate.comp.spv", d.integrateModule, d.integratePipeline)) return false;
    return true;
}

double FemGpuSolver::run(int substeps, const sim::fem::Material& material, float dt, float gravity,
                         float damping) {
    Impl& d = *impl_;

    PushConstants push{};
    push.dt = dt;
    push.lambda = material.lameLambda();
    push.mu = material.lameMu();
    push.gravity = gravity;
    push.damping = damping;
    push.nodeCount = static_cast<uint32_t>(d.nodeCount);
    push.tetCount = static_cast<uint32_t>(d.tetCount);

    VkCommandBuffer cmd = d.beginOneShot();
    vkCmdResetQueryPool(cmd, d.queryPool, 0, 2);
    vkCmdWriteTimestamp(cmd, VK_PIPELINE_STAGE_TOP_OF_PIPE_BIT, d.queryPool, 0);

    vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, d.pipelineLayout, 0, 1,
                            &d.descriptorSet, 0, nullptr);
    vkCmdPushConstants(cmd, d.pipelineLayout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(push), &push);

    // The two kernels alternate and each depends on the other's writes, so every
    // dispatch is separated by a full shader-write to shader-read barrier.
    VkMemoryBarrier barrier{VK_STRUCTURE_TYPE_MEMORY_BARRIER};
    barrier.srcAccessMask = VK_ACCESS_SHADER_WRITE_BIT;
    barrier.dstAccessMask = VK_ACCESS_SHADER_READ_BIT | VK_ACCESS_SHADER_WRITE_BIT;

    for (int step = 0; step < substeps; ++step) {
        vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, d.forcePipeline);
        vkCmdDispatch(cmd, groupsFor(d.tetCount), 1, 1);
        vkCmdPipelineBarrier(cmd, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
                             VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &barrier, 0, nullptr, 0,
                             nullptr);

        vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, d.integratePipeline);
        vkCmdDispatch(cmd, groupsFor(d.nodeCount), 1, 1);
        vkCmdPipelineBarrier(cmd, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
                             VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &barrier, 0, nullptr, 0,
                             nullptr);
    }

    vkCmdWriteTimestamp(cmd, VK_PIPELINE_STAGE_BOTTOM_OF_PIPE_BIT, d.queryPool, 1);
    d.endOneShot(cmd);

    uint64_t stamps[2] = {0, 0};
    vkGetQueryPoolResults(d.device, d.queryPool, 0, 2, sizeof(stamps), stamps, sizeof(uint64_t),
                          VK_QUERY_RESULT_64_BIT | VK_QUERY_RESULT_WAIT_BIT);
    return static_cast<double>(stamps[1] - stamps[0]) * d.timestampPeriod * 1e-6;
}

void FemGpuSolver::readback(sim::fem::TetMesh& mesh) {
    Impl& d = *impl_;
    d.downloadFrom(d.buffers[0], mesh.position.data(), mesh.position.size() * sizeof(float));
    d.downloadFrom(d.buffers[1], mesh.velocity.data(), mesh.velocity.size() * sizeof(float));
}

}  // namespace gpu
