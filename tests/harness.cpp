// SPDX-License-Identifier: MIT
#include "harness.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <string>
#include <vector>

#include <unistd.h>

namespace testing {
namespace {
int checks = 0;
int failures = 0;

// The scratch this process created, in a plain buffer rather than a std::string:
// the exit handler below runs *after* the string inside scratchDir() has been
// destroyed, and a handler that reads a dead object to decide what to delete is
// worse than no handler at all.
char ownScratch[4096] = {};

// What happens to a run's scratch when the run ends.
//
// A run that wrote nothing leaves nothing -- one directory per run is litter the
// fixed names never made, since they were overwritten in place. Everything else
// is *kept*, and announced so it can be found: these files are deliberate
// artefacts. `test_breach.cpp` writes the damage map because "re-deriving it by
// hand from a failing assertion is miserable", the GPU suites write frames and
// print the path so someone will open them, and docs/03-renderer-audio.md sends
// the reader to `hull_damaged_ship.png` by name. Deleting all that on the way
// out would turn every one of those printed paths into a lie, which is a worse
// trade than a directory in /tmp. Sweep them with
// `rm -rf ${TMPDIR:-/tmp}/shipsim-tests-*` -- no live suite is ever harmed by
// that, because a live suite is the only thing that writes into its own.
void keepOrRemoveScratchOnExit() {
    if (ownScratch[0] == '\0') return;
    std::error_code ec;
    const std::filesystem::path directory(ownScratch);
    if (std::filesystem::is_empty(directory, ec) && !ec) {
        std::filesystem::remove(directory, ec);
        return;
    }
    std::printf("test output kept in %s\n", ownScratch);
}

// A directory under `base` that belongs to this process and to nothing else.
//
// `mkdtemp` rather than a name assembled from the pid: it asks the kernel to
// *create* a directory that did not exist and fails if it cannot, so uniqueness
// is a property of the call rather than of a naming convention two processes
// might independently satisfy. A pid alone would not do it -- it is unique only
// within one pid namespace, so two containers sharing a /tmp mount both have a
// process 47, and pids are reused after wraparound. The pid is in the name
// anyway, for whoever is reading `ls /tmp` afterwards and wants to know whose
// scratch this was.
std::string privateDirectoryUnder(const std::string& base) {
    std::error_code ec;
    std::filesystem::create_directories(base, ec);
    if (!std::filesystem::is_directory(base, ec)) return {};
    const std::string pattern = base + "shipsim-tests-" + std::to_string(::getpid()) + "-XXXXXX";
    std::vector<char> buffer(pattern.begin(), pattern.end());
    buffer.push_back('\0');
    if (::mkdtemp(buffer.data()) == nullptr) return {};
    const std::string created = buffer.data();
    // Only claim it for cleanup if the whole path fits: a truncated one names a
    // directory we did not create.
    if (created.size() < sizeof ownScratch) {
        std::snprintf(ownScratch, sizeof ownScratch, "%s", created.c_str());
        std::atexit(keepOrRemoveScratchOnExit);
    }
    return created + '/';
}
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
        //
        // **And the directory is this process's own, not the one it was pointed
        // at.** Tests write *fixed* names into it -- `barge.ship`,
        // `ferry_damage_map.txt`, `collision_ram.csv`, half a dozen PNGs -- and
        // read some of them back, so two suites sharing one directory are two
        // suites deciding each other's verdicts. That has already happened: a
        // mutation sweep running `--workers 4` against one /tmp scored a mutant
        // KILLED on `FAIL a ship written to disk loads back: /tmp/barge.ship:
        // empty`, which is a verdict about a *neighbour's* file and not about the
        // mutation. It deserved killing on its own merits so the measured rate
        // did not move, and that is exactly why it matters -- the next one to
        // lose the race could be a control, and a control that dies reads as the
        // instrument being broken.
        //
        // A caller can avoid this by handing every suite its own $TMPDIR, and the
        // mutation harness does. That is a *declaration and not a property*: the
        // next caller who forgets, or the next test written against some other
        // variable, re-opens the race silently and nothing notices. So the
        // uniqueness is taken here, where no caller can fail to ask for it.
        for (const char* variable : {"SHIPSIM_TEST_TMPDIR", "TMPDIR"}) {
            const char* value = std::getenv(variable);
            if (value == nullptr || *value == '\0') continue;
            std::string candidate = value;
            if (candidate.back() != '/') candidate += '/';
            if (std::string own = privateDirectoryUnder(candidate); !own.empty()) return own;
        }
        if (std::string own = privateDirectoryUnder("/tmp/"); !own.empty()) return own;
        // Nowhere would give us one. Carry on in the shared directory rather than
        // refusing to run -- every test that writes here checks its own write --
        // but say so, because from here on a second suite can collide with this
        // one and the reader needs to know that before believing a failure.
        std::printf("     no private scratch directory could be created:"
                    " falling back to /tmp, where another suite can collide"
                    " with this one\n");
        return std::string("/tmp/");
    }();
    return dir;
}

int checkCount() { return checks; }
int failureCount() { return failures; }

}  // namespace testing
