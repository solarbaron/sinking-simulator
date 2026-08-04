// SPDX-License-Identifier: MIT
//
// Shared assertion harness. Deliberately tiny: the value in this project's tests
// is that they compare against closed-form or independently-derived answers, not
// that the framework is clever.
#pragma once

#include <string>

namespace testing {

void expectNear(const std::string& what, double got, double want, double tolerance);
void expectTrue(const std::string& what, bool condition);
void expectEqual(const std::string& what, long long got, long long want);

int checkCount();
int failureCount();

// A writable directory for test output, with a trailing separator.
//
// Tests that write PNGs used to hard-code an absolute scratch path belonging to
// one machine, which meant they could only ever pass there -- and `full` runs
// the GPU suites, so the committed CI workflow would have failed on its first
// run. Honours $SHIPSIM_TEST_TMPDIR, then $TMPDIR, then falls back to /tmp.
const std::string& scratchDir();

}  // namespace testing

// Each suite lives in its own translation unit and is called from tests/main.cpp.
void runCameraTests();
void runCoreTests();
void runJobTests();
void runArenaTests();
void runBreachTests();
void runEcsTests();
void runDeviceTests();
void runOffscreenTests();
void runOceanTests();
void runHullRenderTests();
void runPngTests();
void runPropulsionTests();
void runRadiationTests();
void runRollDampingTests();
void runScantlingTests();
void runSchedulerTests();
void runSolidShellTests();
void runSerialiseTests();
void runShipFileTests();
void runWorldIoTests();
void runWaveTests();
void runBucklingTests();
void runCollisionTests();
void runGirderTests();
void runHullFormTests();
void runRaoTests();
