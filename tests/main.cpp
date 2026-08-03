// SPDX-License-Identifier: MIT
#include "harness.hpp"

#include <cstdio>

int main() {
    std::printf("shipsim validation\n");
    runCoreTests();
    runJobTests();
    runArenaTests();
    runEcsTests();
    runPngTests();
    runSchedulerTests();
    runSerialiseTests();
    runWorldIoTests();
    std::printf("%d checks, %d failures\n", testing::checkCount(), testing::failureCount());
    return testing::failureCount() == 0 ? 0 : 1;
}
