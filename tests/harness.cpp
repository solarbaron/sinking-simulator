// SPDX-License-Identifier: MIT
#include "harness.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>

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

const std::string& scratchDir() {
    static const std::string dir = [] {
        // First writable candidate wins. A test that cannot write its output
        // should say so rather than silently assert against a file it never
        // created, so the directory is created here and checked.
        for (const char* variable : {"SHIPSIM_TEST_TMPDIR", "TMPDIR"}) {
            const char* value = std::getenv(variable);
            if (value == nullptr || *value == '\0') continue;
            std::string candidate = value;
            if (candidate.back() != '/') candidate += '/';
            std::error_code ec;
            std::filesystem::create_directories(candidate, ec);
            if (std::filesystem::is_directory(candidate, ec)) return candidate;
        }
        return std::string("/tmp/");
    }();
    return dir;
}

int checkCount() { return checks; }
int failureCount() { return failures; }

}  // namespace testing
