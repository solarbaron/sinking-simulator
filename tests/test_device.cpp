// SPDX-License-Identifier: MIT
//
// Validation of the shared Vulkan context.
//
// Every check here is skipped, not failed, when there is no usable GPU: the rest
// of the suite must stay runnable on a headless box. That "skip" is itself
// asserted — a Device that cannot be created must say why rather than crash, and
// must leave nothing behind.
//
// The substantive checks are round trips through real device memory, because the
// failure mode that matters is silent: a buffer copy that transfers the wrong
// size, or a staging path that reuses a too-small buffer, returns success and
// hands back plausible-looking garbage.
#include "engine/gpu/device.hpp"
#include "harness.hpp"

#include <cstdint>
#include <cstdio>
#include <numeric>
#include <string>
#include <vector>

// Only the values this test needs, rather than pulling the Vulkan header into a
// translation unit that otherwise does not want it.
namespace {
constexpr std::uint32_t kUsageStorage = 0x00000020u;   // VK_BUFFER_USAGE_STORAGE_BUFFER_BIT
constexpr std::uint32_t kUsageTransferSrc = 0x00000001u;
constexpr std::uint32_t kUsageTransferDst = 0x00000002u;
constexpr std::uint32_t kMemoryDeviceLocal = 0x00000001u;
constexpr std::uint32_t kFormatRgba8Unorm = 37u;       // VK_FORMAT_R8G8B8A8_UNORM
constexpr std::uint32_t kUsageColourAttachment = 0x00000010u;
constexpr std::uint32_t kAspectColour = 0x00000001u;
}  // namespace

using testing::expectEqual;
using testing::expectTrue;

namespace {

bool gpuReported = false;

// Creating a device is the gate for every other check in this file.
bool tryCreate(gpu::Device& device) {
    std::string error;
    if (device.create(error)) {
        if (!gpuReported) {
            std::printf("     GPU: %s%s\n", device.deviceName().c_str(),
                        device.supportsGraphics() ? " (graphics capable)" : " (compute only)");
            gpuReported = true;
        }
        return true;
    }
    if (!gpuReported) {
        std::printf("     no usable GPU (%s) - GPU checks skipped\n", error.c_str());
        gpuReported = true;
    }
    return false;
}

// A machine with no Vulkan must produce a reason and a clean, reusable object,
// not a crash and not a half-initialised device.
void testAbsentDeviceDegradesCleanly() {
    gpu::Device device;
    std::string error;
    const bool created = device.create(error);
    if (!created) {
        expectTrue("a failed create explains itself", !error.empty());
        expectTrue("a failed create leaves the device invalid", !device.valid());
        // Destroying an object that never initialised must be safe.
        device.destroy();
        expectTrue("destroying an uninitialised device is safe", !device.valid());
        return;
    }
    expectTrue("a successful create reports valid", device.valid());
    expectTrue("a successful create names the device", !device.deviceName().empty());
    expectTrue("timestamp period is positive", device.timestampPeriodNanos() > 0.0f);

    // Destroy then recreate: leaked handles or a stale staging buffer show up here.
    device.destroy();
    expectTrue("destroy leaves the device invalid", !device.valid());
    std::string second;
    expectTrue("a device can be created again after destruction", device.create(second));
}

// Upload then download must return exactly what went in. A copy that transfers
// the wrong byte count still succeeds and still returns data.
void testBufferRoundTrip() {
    gpu::Device device;
    if (!tryCreate(device)) return;

    constexpr std::size_t kCount = 4096;
    std::vector<std::uint32_t> source(kCount);
    std::iota(source.begin(), source.end(), 0x1234u);

    gpu::Device::Buffer buffer = device.createBuffer(
        source.size() * sizeof(std::uint32_t),
        kUsageStorage | kUsageTransferSrc | kUsageTransferDst, kMemoryDeviceLocal);
    expectTrue("a device-local buffer allocates", buffer.valid());
    if (!buffer.valid()) return;

    expectTrue("upload succeeds",
               device.upload(buffer, source.data(), source.size() * sizeof(std::uint32_t)));

    std::vector<std::uint32_t> readback(kCount, 0);
    expectTrue("download succeeds",
               device.download(buffer, readback.data(), readback.size() * sizeof(std::uint32_t)));
    expectTrue("device memory round trip is exact", readback == source);

    // Staging is shared and grows on demand; a second, larger transfer must not
    // reuse a buffer that is too small.
    std::vector<std::uint32_t> bigger(kCount * 4, 0xABCDEFu);
    gpu::Device::Buffer large = device.createBuffer(
        bigger.size() * sizeof(std::uint32_t),
        kUsageStorage | kUsageTransferSrc | kUsageTransferDst, kMemoryDeviceLocal);
    expectTrue("a larger buffer allocates", large.valid());
    expectTrue("a larger transfer grows the staging buffer",
               device.upload(large, bigger.data(), bigger.size() * sizeof(std::uint32_t)));
    std::vector<std::uint32_t> largeBack(bigger.size(), 0);
    expectTrue("the larger buffer reads back",
               device.download(large, largeBack.data(), largeBack.size() * sizeof(std::uint32_t)));
    expectTrue("the larger round trip is exact", largeBack == bigger);

    // And the first buffer must be undisturbed by the staging growth.
    std::vector<std::uint32_t> again(kCount, 0);
    device.download(buffer, again.data(), again.size() * sizeof(std::uint32_t));
    expectTrue("the earlier buffer survived the staging resize", again == source);

    device.destroyBuffer(large);
    device.destroyBuffer(buffer);
    expectTrue("destroying a buffer clears its handle", !buffer.valid());
}

// Transfers that cannot be satisfied must be refused rather than truncated.
void testTransfersRefuseImpossibleSizes() {
    gpu::Device device;
    if (!tryCreate(device)) return;

    gpu::Device::Buffer small =
        device.createBuffer(64, kUsageTransferSrc | kUsageTransferDst, kMemoryDeviceLocal);
    expectTrue("a small buffer allocates", small.valid());

    std::vector<std::uint8_t> tooMuch(4096, 0x5Au);
    expectTrue("uploading more than the buffer holds is refused",
               !device.upload(small, tooMuch.data(), tooMuch.size()));
    expectTrue("downloading more than the buffer holds is refused",
               !device.download(small, tooMuch.data(), tooMuch.size()));
    expectTrue("a zero-byte transfer is refused", !device.upload(small, tooMuch.data(), 0));

    gpu::Device::Buffer invalid;
    expectTrue("transferring to an unallocated buffer is refused",
               !device.upload(invalid, tooMuch.data(), 16));

    device.destroyBuffer(small);
}

void testImageAllocation() {
    gpu::Device device;
    if (!tryCreate(device)) return;

    gpu::Device::Image2D colour =
        device.createImage2D(320, 240, kFormatRgba8Unorm,
                             kUsageColourAttachment | kUsageTransferSrc, kAspectColour);
    expectTrue("a colour attachment allocates", colour.valid());
    expectEqual("the image records its width", static_cast<long long>(colour.width), 320);
    expectEqual("the image records its height", static_cast<long long>(colour.height), 240);
    expectTrue("the image has a view", colour.view != nullptr);

    device.destroyImage(colour);
    expectTrue("destroying an image clears its handle", !colour.valid());

    gpu::Device::Image2D degenerate =
        device.createImage2D(0, 240, kFormatRgba8Unorm, kUsageColourAttachment, kAspectColour);
    expectTrue("a zero-sized image is refused", !degenerate.valid());
}

void testShaderLoading() {
    gpu::Device device;
    if (!tryCreate(device)) return;

    std::string error;
    expectTrue("a missing shader file is reported rather than crashing",
               device.loadShader("/nonexistent/path/to/shader.spv", error) == nullptr);
    expectTrue("the failure names the file", error.find("shader.spv") != std::string::npos);

    // The FEM spike's compute shader is built alongside the tests, so it is a
    // real SPIR-V module to load.
    std::string ok;
    auto* module = device.loadShader(std::string(SHIPSIM_SHADER_DIR) + "/tet_forces.comp.spv", ok);
    if (module != nullptr) {
        expectTrue("a real SPIR-V module loads", true);
        device.destroyShader(module);
    } else {
        std::printf("     (shader dir not populated: %s)\n", ok.c_str());
    }
}

}  // namespace

void runDeviceTests() {
    std::printf("\n--- vulkan device ---\n");
    testAbsentDeviceDegradesCleanly();
    testBufferRoundTrip();
    testTransfersRefuseImpossibleSizes();
    testImageAllocation();
    testShaderLoading();
}
