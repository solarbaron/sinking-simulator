// SPDX-License-Identifier: MIT
#include "harness.hpp"

#include <cstdio>

// The Vulkan loader drags in the NVIDIA driver stack, which drags in libdbus,
// and that leaks on exit. Every frame of every reported leak is inside those
// libraries -- none is ours -- so they are suppressed by library name rather
// than by silencing the leak checker, which would hide our own leaks too.
// Defined here rather than in an external file so it applies however the binary
// is invoked.
#if defined(__SANITIZE_ADDRESS__)
extern "C" const char* __lsan_default_suppressions() {
    return "leak:libdbus-1\n"
           "leak:libnvidia\n"
           "leak:libGLX\n"
           "leak:libvulkan\n";
}
#endif

int main() {
    std::printf("shipsim validation\n");
    runCameraTests();
    runCoreTests();
    runJobTests();
    runArenaTests();
    runEcsTests();
    // ThreadSanitizer aborts inside the uninstrumented GPU driver ("nested bug
    // in the same thread"), which says nothing about our code and would mask the
    // reports that do. The device tests exercise no concurrency of ours, so they
    // sit out the TSan build; jobs, arenas and the scheduler still run under it.
#if defined(SHIPSIM_HAS_VULKAN) && !defined(__SANITIZE_THREAD__)
    runDeviceTests();
    runOffscreenTests();
    runOceanTests();
    runHullRenderTests();
#elif defined(SHIPSIM_HAS_VULKAN)
    std::printf("\n--- vulkan device ---\n     skipped under ThreadSanitizer"
                " (driver is not instrumented)\n");
#endif
    runPngTests();
    runPropulsionTests();
    runRadiationTests();
    runRollDampingTests();
    runSchedulerTests();
    runSerialiseTests();
    runShipFileTests();
    runWorldIoTests();
    runWaveTests();
    runGirderTests();
    runHullFormTests();
    runRaoTests();
    std::printf("%d checks, %d failures\n", testing::checkCount(), testing::failureCount());
    return testing::failureCount() == 0 ? 0 : 1;
}
