// SPDX-License-Identifier: MIT
//
// `tests/test_flip.cpp` alone, as its own binary.
//
// The suite it normally lives in runs everything, which is right for the gate and
// wrong for a mutation sweep: nothing outside `flip.{hpp,cpp}` includes it, so no
// mutation of that file can reach another suite, and paying for all of them per
// mutant multiplies a fifty-one-substitution sweep by five. `tools/flip_probe/mutate.py`
// builds and runs this instead.
//
// It is the *same* translation unit and the same harness as the gate runs, not a
// reduced copy of it -- a mutation harness that runs a different set of assertions
// from the gate is reporting a kill rate for a suite nobody ships.
#include "harness.hpp"

#include <cstdio>

int main() {
    std::printf("shipsim validation (flip only)\n");
    runFlipTests();
    std::printf("%d checks, %d failures\n", testing::checkCount(), testing::failureCount());
    return testing::failureCount() == 0 ? 0 : 1;
}
