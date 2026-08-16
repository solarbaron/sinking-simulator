// SPDX-License-Identifier: MIT
#include "zone_gpu.hpp"

#include <vulkan/vulkan.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>

namespace gpu {
namespace {

using sim::solidshell::kDof;
using sim::solidshell::kEas;
using sim::solidshell::kGauss;
using sim::solidshell::kNodes;

// 14 and 15 are the fp64 enhanced-strain block's: a double copy of G and the Gauss
// weights, and alpha's persistent per-element state. They are in the descriptor set
// layout unconditionally, because one layout serves every shader; they are **allocated
// with real contents only for the kernels that declare them**, and otherwise get a
// four-byte stub, so the float path pays nothing for them.
constexpr int kBufferCount = 16;

// Layout of one element's `RestForms` in the buffer, in floats. Must match the
// constants at the top of solidshell_forces.comp; the shader's copy is checked
// against this one by `tests/test_zone_gpu.cpp` rather than being trusted.
constexpr std::size_t kFormB = 0;
constexpr std::size_t kFormG = kGauss * 6 * kDof;                 // 1152
constexpr std::size_t kFormW = kFormG + kGauss * 6 * kEas;        // 1488
constexpr std::size_t kFormJ = kFormW + kGauss;                   // 1496
constexpr std::size_t kFormStride = kFormJ + 9;                   // 1505

constexpr std::size_t kPointStride = 15;                          // per Gauss point
constexpr std::size_t kStateEnhanced = kGauss * kPointStride;     // 120
constexpr std::size_t kStateFailure = kStateEnhanced + kEas;      // 127
constexpr std::size_t kStateTorn = kStateFailure + 1;             // 128
constexpr std::size_t kStateStride = kStateTorn + 1;              // 129

// The fp64 enhanced forms, in doubles per element. Must match the shader's
// kEasFormG / kEasFormW / kEasFormStride; `tests/test_zone_gpu.cpp` checks them.
constexpr std::size_t kEasFormG = 0;
constexpr std::size_t kEasFormW = kGauss * 6 * kEas;               // 336
constexpr std::size_t kEasFormStride = kEasFormW + kGauss;         // 344


// Which of the five compiled variants of the workgroup kernel a precision selects, and
// which fp64 buffers it needs. All five come from one GLSL source
// (`solidshell_forces_wg.comp`, `-DSHIPSIM_EAS_FP64=`), so the difference between any
// two of them is the precision and nothing else -- which is what makes the ladder a
// measurement rather than five kernels that might also differ in other ways.
const char* workgroupShaderFor(EasPrecision eas) {
    switch (eas) {
        case EasPrecision::FloatTight: return "solidshell_forces_wg_tight.comp.spv";
        case EasPrecision::Solve:      return "solidshell_forces_wg_f64solve.comp.spv";
        case EasPrecision::Condense:   return "solidshell_forces_wg_f64condense.comp.spv";
        case EasPrecision::Newton:     return "solidshell_forces_wg_f64newton.comp.spv";
        case EasPrecision::Float:      break;
    }
    return "solidshell_forces_wg.comp.spv";
}
bool needsDoubleForms(EasPrecision eas) {
    return eas == EasPrecision::Condense || eas == EasPrecision::Newton;
}
bool needsDoubleAlpha(EasPrecision eas) { return eas == EasPrecision::Newton; }
bool needsFloat64(EasPrecision eas) {
    return eas != EasPrecision::Float && eas != EasPrecision::FloatTight;
}

struct PushConstants {
    float dt;
    float mu;
    float kappa;
    float yieldStrength;
    float strengthCoefficient;
    float referenceStrain;
    float hardeningExponent;
    float kinematicModulus;
    float triaxialitySensitivity;
    float referenceTriaxiality;
    float cutoffTriaxiality;
    float punchSpeed;
    float axisX, axisY, axisZ;
    uint32_t nodeCount;
    uint32_t elementCount;
    uint32_t drivenCount;
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

uint32_t groupsFor(std::size_t items, uint32_t size) {
    return static_cast<uint32_t>((items + size - 1) / size);
}

}  // namespace

std::string describeElementPipelines(const std::string& shaderDirectory) {
    std::string out;
    VkInstance instance = VK_NULL_HANDLE;
    VkDevice device = VK_NULL_HANDLE;
    VkApplicationInfo app{VK_STRUCTURE_TYPE_APPLICATION_INFO};
    app.pApplicationName = "shipsim-zone-gpu-stats";
    app.apiVersion = VK_API_VERSION_1_3;
    VkInstanceCreateInfo instanceInfo{VK_STRUCTURE_TYPE_INSTANCE_CREATE_INFO};
    instanceInfo.pApplicationInfo = &app;
    if (vkCreateInstance(&instanceInfo, nullptr, &instance) != VK_SUCCESS)
        return "  (no Vulkan instance)\n";

    uint32_t count = 0;
    vkEnumeratePhysicalDevices(instance, &count, nullptr);
    if (count == 0) { vkDestroyInstance(instance, nullptr); return "  (no Vulkan device)\n"; }
    std::vector<VkPhysicalDevice> devices(count);
    vkEnumeratePhysicalDevices(instance, &count, devices.data());
    VkPhysicalDevice physical = devices[0];
    for (VkPhysicalDevice candidate : devices) {
        VkPhysicalDeviceProperties props;
        vkGetPhysicalDeviceProperties(candidate, &props);
        if (props.deviceType == VK_PHYSICAL_DEVICE_TYPE_DISCRETE_GPU) { physical = candidate; break; }
    }

    uint32_t extensionCount = 0;
    vkEnumerateDeviceExtensionProperties(physical, nullptr, &extensionCount, nullptr);
    std::vector<VkExtensionProperties> extensions(extensionCount);
    vkEnumerateDeviceExtensionProperties(physical, nullptr, &extensionCount, extensions.data());
    bool haveStats = false;
    for (const VkExtensionProperties& ext : extensions)
        if (std::strcmp(ext.extensionName, VK_KHR_PIPELINE_EXECUTABLE_PROPERTIES_EXTENSION_NAME) == 0)
            haveStats = true;
    if (!haveStats) {
        vkDestroyInstance(instance, nullptr);
        return "  (driver has no VK_KHR_pipeline_executable_properties)\n";
    }

    uint32_t families = 0;
    vkGetPhysicalDeviceQueueFamilyProperties(physical, &families, nullptr);
    std::vector<VkQueueFamilyProperties> familyProps(families);
    vkGetPhysicalDeviceQueueFamilyProperties(physical, &families, familyProps.data());
    uint32_t family = 0;
    for (uint32_t i = 0; i < families; ++i)
        if (familyProps[i].queueFlags & VK_QUEUE_COMPUTE_BIT) { family = i; break; }

    const float priority = 1.0f;
    VkDeviceQueueCreateInfo queueInfo{VK_STRUCTURE_TYPE_DEVICE_QUEUE_CREATE_INFO};
    queueInfo.queueFamilyIndex = family;
    queueInfo.queueCount = 1;
    queueInfo.pQueuePriorities = &priority;
    VkPhysicalDevicePipelineExecutablePropertiesFeaturesKHR feature{
        VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_PIPELINE_EXECUTABLE_PROPERTIES_FEATURES_KHR};
    feature.pipelineExecutableInfo = VK_TRUE;
    const char* wanted[] = {VK_KHR_PIPELINE_EXECUTABLE_PROPERTIES_EXTENSION_NAME};
    VkDeviceCreateInfo deviceInfo{VK_STRUCTURE_TYPE_DEVICE_CREATE_INFO};
    deviceInfo.pNext = &feature;
    deviceInfo.queueCreateInfoCount = 1;
    deviceInfo.pQueueCreateInfos = &queueInfo;
    deviceInfo.enabledExtensionCount = 1;
    deviceInfo.ppEnabledExtensionNames = wanted;
    // fp64 as well, or the enhanced-block variants below cannot be created and the one
    // number this report exists to give -- what the double block costs in registers and
    // spill -- would be missing exactly where it matters.
    VkPhysicalDeviceFeatures statsFeatures{};
    VkPhysicalDeviceFeatures availableFeatures{};
    vkGetPhysicalDeviceFeatures(physical, &availableFeatures);
    statsFeatures.shaderFloat64 = availableFeatures.shaderFloat64;
    deviceInfo.pEnabledFeatures = &statsFeatures;
    if (vkCreateDevice(physical, &deviceInfo, nullptr, &device) != VK_SUCCESS) {
        vkDestroyInstance(instance, nullptr);
        return "  (vkCreateDevice failed with the statistics extension)\n";
    }

    auto getStatistics = reinterpret_cast<PFN_vkGetPipelineExecutableStatisticsKHR>(
        vkGetDeviceProcAddr(device, "vkGetPipelineExecutableStatisticsKHR"));
    auto getExecutables = reinterpret_cast<PFN_vkGetPipelineExecutablePropertiesKHR>(
        vkGetDeviceProcAddr(device, "vkGetPipelineExecutablePropertiesKHR"));

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
    VkDescriptorSetLayout setLayout = VK_NULL_HANDLE;
    vkCreateDescriptorSetLayout(device, &layoutInfo, nullptr, &setLayout);
    VkPushConstantRange pushRange{VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(PushConstants)};
    VkPipelineLayoutCreateInfo pipelineLayoutInfo{VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO};
    pipelineLayoutInfo.setLayoutCount = 1;
    pipelineLayoutInfo.pSetLayouts = &setLayout;
    pipelineLayoutInfo.pushConstantRangeCount = 1;
    pipelineLayoutInfo.pPushConstantRanges = &pushRange;
    VkPipelineLayout pipelineLayout = VK_NULL_HANDLE;
    vkCreatePipelineLayout(device, &pipelineLayoutInfo, nullptr, &pipelineLayout);

    // **The subgroup size, because it decides whether `solidshell_forces_wg.comp`'s
    // barriers do anything on this machine.** That kernel declares
    // `local_size_x = 32`. If the device's subgroup is also 32 then a workgroup is
    // exactly one subgroup, executing in lockstep, and every `barrier()` in it is
    // synchronising threads the hardware has already synchronised -- so removing one
    // changes nothing here and would change everything on a device with independent
    // thread scheduling. Mutation testing found exactly that: three of the four
    // barrier mutants survive. This number is what turns that from a guess into a
    // mechanism, so it is printed next to the statistics that motivated looking.
    VkPhysicalDeviceSubgroupProperties subgroup{
        VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_SUBGROUP_PROPERTIES};
    VkPhysicalDeviceProperties2 props2{VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_PROPERTIES_2};
    props2.pNext = &subgroup;
    vkGetPhysicalDeviceProperties2(physical, &props2);
    // Built by concatenation rather than into a fixed buffer: `deviceName` is 256
    // bytes and any buffer short enough to be comfortable is short enough to
    // truncate it, which `-Wformat-truncation` says out loud.
    out += "  device ";
    out += props2.properties.deviceName;
    out += ": subgroup size " + std::to_string(subgroup.subgroupSize) +
           ", workgroup of the remapped kernel 32 -- ";
    out += subgroup.subgroupSize >= 32
               ? "so one workgroup is one subgroup and its barriers are no-ops here\n"
               : "so a workgroup spans subgroups and its barriers are load-bearing\n";

    // **`node_integrate.comp` is the calibration, and it is here rather than in a
    // constant on purpose.** This driver reports `Local Memory Size` with a large
    // fixed offset -- 2^36 on the 1070 Ti -- so the raw number says nothing on its
    // own. Rather than hardcode that, the first shader queried is one known to spill
    // nothing (it is a handful of loads, a multiply-add and a store), and every
    // spill figure below is quoted against whatever *it* reports. If a driver ever
    // stops adding the offset the baseline goes to zero and the deltas still read
    // correctly. Confirmed against `tet_forces.comp`, which reports the identical
    // baseline -- which is also §8's claim that a linear tet's state fits in
    // registers, checked rather than repeated.
    const char* files[6] = {"node_integrate.comp.spv",
                            "solidshell_forces.comp.spv",
                            "solidshell_forces_wg.comp.spv",
                            "solidshell_forces_wg_f64solve.comp.spv",
                            "solidshell_forces_wg_f64condense.comp.spv",
                            "solidshell_forces_wg_f64newton.comp.spv"};
    const char* labels[6] = {"node_integrate.comp (spill-free calibration)",
                             "one invocation per element",
                             "one workgroup per element",
                             "  + the 7x7 linear solve in fp64",
                             "  + Kaa and the residual condensed in fp64",
                             "  + alpha and the Newton loop in fp64"};
    long long localBaseline = -1;
    for (int which = 0; which < 6; ++which) {
        out += "  ";
        out += labels[which];
        out += " (";
        out += files[which];
        out += ")\n";
        const std::vector<uint32_t> spirv = readSpirv(shaderDirectory + "/" + files[which]);
        if (spirv.empty()) { out += "    (shader not found)\n"; continue; }
        VkShaderModuleCreateInfo moduleInfo{VK_STRUCTURE_TYPE_SHADER_MODULE_CREATE_INFO};
        moduleInfo.codeSize = spirv.size() * 4;
        moduleInfo.pCode = spirv.data();
        VkShaderModule module = VK_NULL_HANDLE;
        if (vkCreateShaderModule(device, &moduleInfo, nullptr, &module) != VK_SUCCESS) continue;
        VkComputePipelineCreateInfo info{VK_STRUCTURE_TYPE_COMPUTE_PIPELINE_CREATE_INFO};
        info.flags = VK_PIPELINE_CREATE_CAPTURE_STATISTICS_BIT_KHR;
        info.stage.sType = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO;
        info.stage.stage = VK_SHADER_STAGE_COMPUTE_BIT;
        info.stage.module = module;
        info.stage.pName = "main";
        info.layout = pipelineLayout;
        VkPipeline pipeline = VK_NULL_HANDLE;
        if (vkCreateComputePipelines(device, VK_NULL_HANDLE, 1, &info, nullptr, &pipeline) !=
            VK_SUCCESS) {
            out += "    (pipeline creation failed)\n";
            vkDestroyShaderModule(device, module, nullptr);
            continue;
        }
        VkPipelineInfoKHR pipelineInfo{VK_STRUCTURE_TYPE_PIPELINE_INFO_KHR};
        pipelineInfo.pipeline = pipeline;
        uint32_t executables = 0;
        getExecutables(device, &pipelineInfo, &executables, nullptr);
        std::vector<VkPipelineExecutablePropertiesKHR> props(
            executables, {VK_STRUCTURE_TYPE_PIPELINE_EXECUTABLE_PROPERTIES_KHR});
        getExecutables(device, &pipelineInfo, &executables, props.data());
        for (uint32_t i = 0; i < executables; ++i) {
            VkPipelineExecutableInfoKHR executableInfo{VK_STRUCTURE_TYPE_PIPELINE_EXECUTABLE_INFO_KHR};
            executableInfo.pipeline = pipeline;
            executableInfo.executableIndex = i;
            uint32_t statistics = 0;
            getStatistics(device, &executableInfo, &statistics, nullptr);
            std::vector<VkPipelineExecutableStatisticKHR> stats(
                statistics, {VK_STRUCTURE_TYPE_PIPELINE_EXECUTABLE_STATISTIC_KHR});
            getStatistics(device, &executableInfo, &statistics, stats.data());
            // **Built by concatenation, not into a fixed buffer.** `s.name` is a
            // 256-byte array the driver fills, so any `snprintf` into a buffer small
            // enough to be convenient can truncate it -- which
            // `-Wformat-truncation` says, but only at the optimisation level the
            // sanitizer build uses. `CLAUDE.md`: "a build cannot see a warning its
            // optimisation level does not produce", and this file proved it again.
            const auto padded = [](const char* name) {
                std::string s(name);
                if (s.size() < 28) s.append(28 - s.size(), ' ');
                return "    " + s + " ";
            };
            for (const VkPipelineExecutableStatisticKHR& s : stats) {
                out += padded(s.name);
                switch (s.format) {
                    case VK_PIPELINE_EXECUTABLE_STATISTIC_FORMAT_INT64_KHR:
                        out += std::to_string(s.value.i64);
                        break;
                    case VK_PIPELINE_EXECUTABLE_STATISTIC_FORMAT_UINT64_KHR:
                        out += std::to_string(s.value.u64);
                        break;
                    case VK_PIPELINE_EXECUTABLE_STATISTIC_FORMAT_FLOAT64_KHR:
                        out += std::to_string(s.value.f64);
                        break;
                    default:
                        out += s.value.b32 ? "true" : "false";
                        break;
                }
                out += "\n";
                if (std::strcmp(s.name, "Local Memory Size") == 0) {
                    const long long raw = static_cast<long long>(s.value.u64);
                    if (which == 0) {
                        localBaseline = raw;
                    } else if (localBaseline >= 0) {
                        out += padded("  spill over calibration") +
                               std::to_string(raw - localBaseline) + " bytes/thread (" +
                               std::to_string((raw - localBaseline) / 4) + " floats)\n";
                    }
                }
            }
        }
        vkDestroyPipeline(device, pipeline, nullptr);
        vkDestroyShaderModule(device, module, nullptr);
    }

    vkDestroyPipelineLayout(device, pipelineLayout, nullptr);
    vkDestroyDescriptorSetLayout(device, setLayout, nullptr);
    vkDestroyDevice(device, nullptr);
    vkDestroyInstance(instance, nullptr);
    return out;
}

struct ZoneGpuSolver::Impl {
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
    std::size_t elementCount = 0;
    std::size_t drivenCount = 0;
    PushConstants push{};
    Mapping mapping = Mapping::Workgroup;
    sim::zone::Indenter indenter{};
    double timestep = 0;
    double time = 0;
    double penetration = 0;
    double origin[3] = {0, 0, 0};  // the patch centroid every device coordinate is relative to

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
        if (type == UINT32_MAX) return false;

        VkMemoryAllocateInfo alloc{VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO};
        alloc.allocationSize = req.size;
        alloc.memoryTypeIndex = type;
        if (vkAllocateMemory(device, &alloc, nullptr, &out.memory) != VK_SUCCESS) return false;
        vkBindBufferMemory(device, out.handle, out.memory, 0);
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

ZoneGpuSolver::~ZoneGpuSolver() {
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

bool ZoneGpuSolver::initialise(const sim::zone::Patch& patch, const sim::zone::Solver& cpu,
                               const sim::plasticity::Material& material,
                               const sim::zone::SolveParams& params,
                               const std::string& shaderDirectory, std::string& error,
                               Mapping mapping, EasPrecision eas) {
    if (cpu.restForms().size() != patch.elementCount()) {
        error = "the CPU solver was built without SolveParams::cacheRestForms, so there is "
                "nothing to upload";
        return false;
    }
    if (!params.plastic) {
        error = "the GPU path implements the elastoplastic element only -- the elastic one is "
                "0.27 us against 7.3 and is not what costs anything";
        return false;
    }
    // The fp64 enhanced block exists for the workgroup mapping only. Saying so is
    // better than quietly running the float kernel, which is a comparison that would
    // pass every tolerance while measuring nothing -- this repo's most common shape of
    // vacuous test.
    if (eas != EasPrecision::Float && mapping != Mapping::Workgroup) {
        error = "the fp64 enhanced-strain variants are compiled for the workgroup mapping "
                "only; the one-invocation kernel is float throughout";
        return false;
    }

    impl_ = new Impl();
    Impl& d = *impl_;
    d.mapping = mapping;
    d.nodeCount = patch.nodeCount();
    d.elementCount = patch.elementCount();
    d.indenter = params.indenter;
    d.timestep = params.timestep > 0 ? params.timestep
                                     : patch.criticalTimestep * (params.timestepSafety / 0.9);

    VkApplicationInfo app{VK_STRUCTURE_TYPE_APPLICATION_INFO};
    app.pApplicationName = "shipsim-zone-gpu";
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
    d.physical = devices[0];
    for (VkPhysicalDevice candidate : devices) {
        VkPhysicalDeviceProperties props;
        vkGetPhysicalDeviceProperties(candidate, &props);
        if (props.deviceType == VK_PHYSICAL_DEVICE_TYPE_DISCRETE_GPU) {
            d.physical = candidate;
            break;
        }
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
    // **fp64 is a feature and has to be asked for.** A shader carrying the Float64
    // capability against a device created without `shaderFloat64` is undefined
    // behaviour that this driver happens to run, so it is requested explicitly and its
    // absence **skips rather than fails**, like a missing device. It is only requested
    // when a variant needs it, so the float path's device is unchanged.
    VkPhysicalDeviceFeatures wantedFeatures{};
    if (needsFloat64(eas)) {
        VkPhysicalDeviceFeatures available{};
        vkGetPhysicalDeviceFeatures(d.physical, &available);
        if (!available.shaderFloat64) {
            error = "this device has no shaderFloat64, so the fp64 enhanced-strain "
                    "variants cannot run on it";
            return false;
        }
        wantedFeatures.shaderFloat64 = VK_TRUE;
        deviceInfo.pEnabledFeatures = &wantedFeatures;
    }
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

    // --- The host-side mirrors of every buffer -----------------------------------
    //
    // Every one of these comes from the CPU solver's own setup rather than being
    // rebuilt here. The narrowing to float is the only transformation.
    const std::size_t n = d.nodeCount, e = d.elementCount;

    // **Everything is solved about the patch's own centroid, and that is not a
    // tidiness measure -- it is what makes float viable at all.**
    //
    // The element works on `u = R^T x - X`, a difference of two node positions.
    // On this ferry a side patch sits at (0, -9.9, 8), so |x| is about 13 m, where
    // float's 24-bit mantissa resolves 1e-6 m -- and the displacements the first
    // thousand steps produce are 1e-5 m. The difference is then computed from two
    // numbers that agree in six of their seven digits, and the answer is mostly
    // rounding. Measured before this shift: node positions 4e-7 m apart after a
    // **single** step, and by 5 500 steps the GPU had torn 58 elements where the
    // CPU tore none.
    //
    // Shifting to the centroid puts every coordinate inside the zone radius, so the
    // same absolute resolution now sits against a coordinate of ~2 m rather than
    // ~13. It costs nothing: the shift is common to `rest` and `position`, so `u`,
    // the Jacobians, the forms and the punch axis are all unchanged in exact
    // arithmetic. `readback` adds it back.
    d.origin[0] = patch.centre.x;
    d.origin[1] = patch.centre.y;
    d.origin[2] = patch.centre.z;

    std::vector<float> position(n * 3), velocity(n * 3, 0.0f), mass(n), restPosition(n * 3);
    for (std::size_t i = 0; i < n * 3; ++i) {
        const double shift = d.origin[i % 3];
        position[i] = static_cast<float>(cpu.position()[i] - shift);
        restPosition[i] = static_cast<float>(cpu.rest()[i] - shift);
        velocity[i] = static_cast<float>(cpu.velocity()[i]);
    }
    for (std::size_t i = 0; i < n; ++i) mass[i] = static_cast<float>(cpu.nodalMass()[i]);

    std::vector<uint32_t> pinned(n * 3, 0);
    for (std::size_t i = 0; i < n * 3 && i < cpu.pinnedDof().size(); ++i)
        pinned[i] = cpu.pinnedDof()[i] ? 1u : 0u;

    std::vector<float> forms(e * kFormStride, 0.0f);
    for (std::size_t el = 0; el < e; ++el) {
        const sim::solidshell::RestForms& f = cpu.restForms()[el];
        float* out = forms.data() + el * kFormStride;
        for (int gp = 0; gp < kGauss; ++gp)
            for (int i = 0; i < 6; ++i)
                for (int j = 0; j < kDof; ++j)
                    out[kFormB + static_cast<std::size_t>(gp) * 6 * kDof +
                        static_cast<std::size_t>(i) * kDof + static_cast<std::size_t>(j)] =
                        static_cast<float>(f.b[gp][i][j]);
        for (int gp = 0; gp < kGauss; ++gp)
            for (int i = 0; i < 6; ++i)
                for (int k = 0; k < kEas; ++k)
                    out[kFormG + static_cast<std::size_t>(gp) * 6 * kEas +
                        static_cast<std::size_t>(i) * kEas + static_cast<std::size_t>(k)] =
                        static_cast<float>(f.g[gp][i][k]);
        for (int gp = 0; gp < kGauss; ++gp)
            out[kFormW + static_cast<std::size_t>(gp)] = static_cast<float>(f.weight[gp]);
        for (int i = 0; i < 9; ++i)
            out[kFormJ + static_cast<std::size_t>(i)] = static_cast<float>(f.restJacobianInverse[i]);
    }

    std::vector<float> state(e * kStateStride, 0.0f);
    for (std::size_t el = 0; el < e; ++el) {
        const sim::solidshell::ElementPlasticState& s = cpu.elementState()[el];
        float* out = state.data() + el * kStateStride;
        for (int gp = 0; gp < kGauss; ++gp) {
            float* p = out + static_cast<std::size_t>(gp) * kPointStride;
            for (int i = 0; i < 6; ++i) {
                p[i] = static_cast<float>(s.point[gp].plasticStrain[i]);
                p[6 + i] = static_cast<float>(s.point[gp].backStress[i]);
            }
            p[12] = static_cast<float>(s.point[gp].equivalentPlasticStrain);
            p[13] = static_cast<float>(s.point[gp].damage);
            p[14] = s.point[gp].failed ? 1.0f : 0.0f;
        }
        for (int k = 0; k < kEas; ++k)
            out[kStateEnhanced + static_cast<std::size_t>(k)] = static_cast<float>(s.enhanced[k]);
        out[kStateFailure] = static_cast<float>(s.failureStrain);
        out[kStateTorn] = s.torn ? 1.0f : 0.0f;
    }

    // **The enhanced-strain operator in double, and only when a kernel reads it.** A
    // double accumulation of `sum w G^T C G` over a G that was rounded to float is an
    // operator known to seven digits however it is summed, so a level that widens the
    // arithmetic without widening G would measure nothing. B is deliberately *not*
    // copied: it is 1 152 of the 1 505 floats and the enhanced block never touches it.
    std::vector<double> easForms(needsDoubleForms(eas) ? e * kEasFormStride : 1, 0.0);
    if (needsDoubleForms(eas))
        for (std::size_t el = 0; el < e; ++el) {
            const sim::solidshell::RestForms& f = cpu.restForms()[el];
            double* out = easForms.data() + el * kEasFormStride;
            for (int gp = 0; gp < kGauss; ++gp)
                for (int i = 0; i < 6; ++i)
                    for (int k = 0; k < kEas; ++k)
                        out[kEasFormG + static_cast<std::size_t>(gp) * 6 * kEas +
                            static_cast<std::size_t>(i) * kEas + static_cast<std::size_t>(k)] =
                            f.g[gp][i][k];
            for (int gp = 0; gp < kGauss; ++gp)
                out[kEasFormW + static_cast<std::size_t>(gp)] = f.weight[gp];
        }
    // alpha's persistent state in double. The float copy in `plastic` is still written
    // by the kernel and is what `readback` reads; this is what the *next* step warm
    // starts from, which is the half a float round trip would quietly undo.
    std::vector<double> easAlpha(needsDoubleAlpha(eas) ? e * kEas : 1, 0.0);
    if (needsDoubleAlpha(eas))
        for (std::size_t el = 0; el < e; ++el)
            for (int k = 0; k < kEas; ++k)
                easAlpha[el * kEas + static_cast<std::size_t>(k)] =
                    cpu.elementState()[el].enhanced[k];

    const std::vector<uint32_t> index(patch.mesh.index.begin(), patch.mesh.index.end());
    const std::vector<uint32_t> driven(cpu.drivenNodes().begin(), cpu.drivenNodes().end());
    d.drivenCount = driven.size();
    const std::vector<float> elementForce(e * kDof, 0.0f);
    const std::vector<float> elementOut(e, 0.0f);
    const std::vector<int32_t> accumulator(2, 0);

    const VkDeviceSize sizes[kBufferCount] = {
        position.size() * sizeof(float),                 // 0
        velocity.size() * sizeof(float),                 // 1
        mass.size() * sizeof(float),                     // 2
        pinned.size() * sizeof(uint32_t),                // 3
        index.size() * sizeof(uint32_t),                 // 4
        forms.size() * sizeof(float),                    // 5
        restPosition.size() * sizeof(float),             // 6
        elementForce.size() * sizeof(float),             // 7
        cpu.adjacencyOffset().size() * sizeof(uint32_t), // 8
        cpu.adjacencyEntry().size() * sizeof(uint32_t),  // 9
        state.size() * sizeof(float),                    // 10
        elementOut.size() * sizeof(float),               // 11
        std::max<VkDeviceSize>(driven.size() * sizeof(uint32_t), 4),  // 12
        accumulator.size() * sizeof(int32_t),            // 13
        easForms.size() * sizeof(double),                // 14
        easAlpha.size() * sizeof(double),                // 15
    };
    const void* data[kBufferCount] = {
        position.data(),     velocity.data(),  mass.data(),         pinned.data(),
        index.data(),        forms.data(),     restPosition.data(), elementForce.data(),
        cpu.adjacencyOffset().data(), cpu.adjacencyEntry().data(),
        state.data(),        elementOut.data(), driven.data(),      accumulator.data(),
        easForms.data(),     easAlpha.data(),
    };

    VkDeviceSize largest = 0;
    for (VkDeviceSize s : sizes) largest = std::max(largest, s);
    if (!d.createBuffer(d.staging, largest,
                        VK_BUFFER_USAGE_TRANSFER_SRC_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT,
                        VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT |
                            VK_MEMORY_PROPERTY_HOST_COHERENT_BIT)) {
        error = "staging buffer allocation failed";
        return false;
    }
    for (int i = 0; i < kBufferCount; ++i) {
        if (!d.createBuffer(d.buffers[i], sizes[i],
                            VK_BUFFER_USAGE_STORAGE_BUFFER_BIT |
                                VK_BUFFER_USAGE_TRANSFER_SRC_BIT |
                                VK_BUFFER_USAGE_TRANSFER_DST_BIT,
                            VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT)) {
            error = "device buffer allocation failed";
            return false;
        }
        if (sizes[i] > 0 && data[i] != nullptr) d.uploadTo(d.buffers[i], data[i], sizes[i]);
    }

    VkDescriptorSetLayoutBinding bindings[kBufferCount]{};
    for (int i = 0; i < kBufferCount; ++i) {
        bindings[i].binding = static_cast<uint32_t>(i);
        bindings[i].descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
        bindings[i].descriptorCount = 1;
        bindings[i].stageFlags = VK_SHADER_STAGE_COMPUTE_BIT;
    }
    VkDescriptorSetLayoutCreateInfo layoutInfo{
        VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO};
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
    const char* forceShader = mapping == Mapping::Workgroup ? workgroupShaderFor(eas)
                                                            : "solidshell_forces.comp.spv";
    elementShader_ = forceShader;
    if (!makePipeline(forceShader, d.forceModule, d.forcePipeline)) return false;
    if (!makePipeline("solidshell_integrate.comp.spv", d.integrateModule, d.integratePipeline))
        return false;

    d.push.dt = static_cast<float>(d.timestep);
    d.push.mu = static_cast<float>(material.shearModulus());
    d.push.kappa = static_cast<float>(material.bulkModulus());
    d.push.yieldStrength = static_cast<float>(material.flow.yieldStrength);
    d.push.strengthCoefficient = static_cast<float>(material.flow.strengthCoefficient);
    d.push.referenceStrain = static_cast<float>(material.flow.referenceStrain);
    d.push.hardeningExponent = static_cast<float>(material.flow.hardeningExponent);
    d.push.kinematicModulus = static_cast<float>(material.flow.kinematicModulus);
    d.push.triaxialitySensitivity = static_cast<float>(material.failure.triaxialitySensitivity);
    d.push.referenceTriaxiality = static_cast<float>(material.failure.referenceTriaxiality);
    d.push.cutoffTriaxiality = static_cast<float>(material.failure.cutoffTriaxiality);
    d.push.axisX = static_cast<float>(patch.axis.x);
    d.push.axisY = static_cast<float>(patch.axis.y);
    d.push.axisZ = static_cast<float>(patch.axis.z);
    d.push.nodeCount = static_cast<uint32_t>(n);
    d.push.elementCount = static_cast<uint32_t>(e);
    d.push.drivenCount = static_cast<uint32_t>(d.drivenCount);
    return true;
}

double ZoneGpuSolver::run(int substeps) {
    Impl& d = *impl_;
    VkCommandBuffer cmd = d.beginOneShot();
    vkCmdResetQueryPool(cmd, d.queryPool, 0, 2);
    vkCmdWriteTimestamp(cmd, VK_PIPELINE_STAGE_TOP_OF_PIPE_BIT, d.queryPool, 0);
    vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, d.pipelineLayout, 0, 1,
                            &d.descriptorSet, 0, nullptr);

    VkMemoryBarrier barrier{VK_STRUCTURE_TYPE_MEMORY_BARRIER};
    barrier.srcAccessMask = VK_ACCESS_SHADER_WRITE_BIT;
    barrier.dstAccessMask = VK_ACCESS_SHADER_READ_BIT | VK_ACCESS_SHADER_WRITE_BIT;

    for (int step = 0; step < substeps; ++step) {
        // The ramp is a per-step quantity and push constants are a *command*, so it
        // varies inside one recorded buffer and the CPU still never sees a step.
        double speed = d.indenter.speed;
        if (d.indenter.rampTime > 0 && d.time < d.indenter.rampTime)
            speed = d.indenter.speed * (d.time / d.indenter.rampTime);
        d.push.punchSpeed = static_cast<float>(speed);
        vkCmdPushConstants(cmd, d.pipelineLayout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(d.push),
                           &d.push);

        vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, d.forcePipeline);
        // One workgroup per element, or one *invocation* per element packed 32 to a
        // workgroup. The two shaders declare the same `local_size_x = 32`, so the only
        // difference the host sees is how many groups an element costs.
        vkCmdDispatch(cmd,
                      d.mapping == Mapping::Workgroup
                          ? static_cast<uint32_t>(d.elementCount)
                          : groupsFor(d.elementCount, kElementGroup),
                      1, 1);
        vkCmdPipelineBarrier(cmd, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
                             VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &barrier, 0, nullptr, 0,
                             nullptr);

        vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, d.integratePipeline);
        vkCmdDispatch(cmd, groupsFor(d.nodeCount, kNodeGroup), 1, 1);
        vkCmdPipelineBarrier(cmd, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
                             VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &barrier, 0, nullptr, 0,
                             nullptr);

        d.penetration += speed * d.timestep;
        d.time += d.timestep;
        ++steps_;
    }

    vkCmdWriteTimestamp(cmd, VK_PIPELINE_STAGE_BOTTOM_OF_PIPE_BIT, d.queryPool, 1);
    d.endOneShot(cmd);

    uint64_t stamps[2] = {0, 0};
    vkGetQueryPoolResults(d.device, d.queryPool, 0, 2, sizeof(stamps), stamps, sizeof(uint64_t),
                          VK_QUERY_RESULT_64_BIT | VK_QUERY_RESULT_WAIT_BIT);
    return static_cast<double>(stamps[1] - stamps[0]) * d.timestampPeriod * 1e-6;
}

ZoneGpuState ZoneGpuSolver::readback() const {
    Impl& d = *impl_;
    ZoneGpuState out;
    const std::size_t n = d.nodeCount, e = d.elementCount;

    std::vector<float> position(n * 3), velocity(n * 3), state(e * kStateStride),
        elementOut(e);
    int32_t accumulator[2] = {0, 0};
    d.downloadFrom(d.buffers[0], position.data(), position.size() * sizeof(float));
    d.downloadFrom(d.buffers[1], velocity.data(), velocity.size() * sizeof(float));
    d.downloadFrom(d.buffers[10], state.data(), state.size() * sizeof(float));
    d.downloadFrom(d.buffers[11], elementOut.data(), elementOut.size() * sizeof(float));
    d.downloadFrom(d.buffers[13], accumulator, sizeof(accumulator));

    out.position.resize(n * 3);
    out.velocity.resize(n * 3);
    for (std::size_t i = 0; i < n * 3; ++i) {
        out.position[i] = static_cast<double>(position[i]) + d.origin[i % 3];
        out.velocity[i] = velocity[i];
    }
    out.plastic.resize(e);
    for (std::size_t el = 0; el < e; ++el) {
        const float* in = state.data() + el * kStateStride;
        sim::solidshell::ElementPlasticState& s = out.plastic[el];
        for (int gp = 0; gp < kGauss; ++gp) {
            const float* p = in + static_cast<std::size_t>(gp) * kPointStride;
            for (int i = 0; i < 6; ++i) {
                s.point[gp].plasticStrain[i] = p[i];
                s.point[gp].backStress[i] = p[6 + i];
            }
            s.point[gp].equivalentPlasticStrain = p[12];
            s.point[gp].damage = p[13];
            s.point[gp].failed = p[14] != 0.0f;
        }
        for (int k = 0; k < kEas; ++k) s.enhanced[k] = in[kStateEnhanced + static_cast<std::size_t>(k)];
        s.failureStrain = in[kStateFailure];
        s.torn = in[kStateTorn] != 0.0f;
        if (s.torn) ++out.tornElements;
        // Each element kept its own running total, so summing them here is a fixed
        // order over a fixed set -- no atomic on the device and no dependence on
        // which invocation finished first. The CPU sums the other way round (over
        // elements every step), so the two totals differ in rounding by
        // construction; that is a different summation order and not a different
        // model, and the comparison tool reports the gap rather than hiding it.
        out.dissipation += elementOut[el];
    }
    out.work = static_cast<double>(accumulator[0]) / kWorkScale;
    out.steps = steps_;
    out.time = d.time;
    out.penetration = d.penetration;
    return out;
}

}  // namespace gpu
