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

}  // namespace testing

// Each suite lives in its own translation unit and is called from tests/main.cpp.
void runCoreTests();
void runJobTests();
void runArenaTests();
void runEcsTests();
