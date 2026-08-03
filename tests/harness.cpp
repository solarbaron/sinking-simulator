// SPDX-License-Identifier: MIT
#include "harness.hpp"

#include <cmath>
#include <cstdio>

namespace testing {
namespace {
int checks = 0;
int failures = 0;
}  // namespace

void expectNear(const std::string& what, double got, double want, double tolerance) {
    ++checks;
    if (std::abs(got - want) > tolerance) {
        std::printf("  FAIL %-52s got %+.9g  want %+.9g  (tol %.2g)\n", what.c_str(), got, want,
                    tolerance);
        ++failures;
    }
}

void expectTrue(const std::string& what, bool condition) {
    ++checks;
    if (!condition) {
        std::printf("  FAIL %s\n", what.c_str());
        ++failures;
    }
}

void expectEqual(const std::string& what, long long got, long long want) {
    ++checks;
    if (got != want) {
        std::printf("  FAIL %-52s got %lld  want %lld\n", what.c_str(), got, want);
        ++failures;
    }
}

int checkCount() { return checks; }
int failureCount() { return failures; }

}  // namespace testing
