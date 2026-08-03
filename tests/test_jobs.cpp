// SPDX-License-Identifier: MIT
//
// Validation of the work-stealing job system.
//
// Concurrency bugs do not announce themselves, so every check here is an
// invariant that must hold on every run rather than a timing observation:
// exactly-once execution, complete range coverage, no deadlock under nesting,
// and bit-identical reductions regardless of worker count.
#include "engine/core/jobs.hpp"
#include "harness.hpp"

#include <atomic>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

using core::JobSystem;
using testing::expectEqual;
using testing::expectTrue;

namespace {

// Every submitted job must run exactly once -- not "at least once", which a
// double-steal bug would still satisfy, and not "roughly N times", which is what
// a racy counter looks like.
void testEveryJobRunsExactlyOnce() {
    for (unsigned workers : {0u, 1u, 4u, 15u}) {
        JobSystem jobs(workers);
        constexpr int kJobs = 20000;
        std::vector<std::atomic<int>> visits(kJobs);
        for (auto& v : visits) v.store(0);

        JobSystem::Counter counter;
        for (int i = 0; i < kJobs; ++i) {
            std::atomic<int>* slot = &visits[i];
            jobs.run(counter, [slot] { slot->fetch_add(1, std::memory_order_relaxed); });
        }
        jobs.wait(counter);

        int wrong = 0;
        for (auto& v : visits)
            if (v.load() != 1) ++wrong;
        expectEqual("every job ran exactly once, workers=" + std::to_string(workers), wrong, 0);
        expectTrue("counter drained, workers=" + std::to_string(workers), counter.complete());
    }
}

// parallelFor must tile its range: no element skipped, none visited twice, and
// no chunk straying outside the requested bounds.
void testParallelForCoversRangeExactlyOnce() {
    for (unsigned workers : {0u, 3u, 15u}) {
        JobSystem jobs(workers);
        constexpr std::size_t kBegin = 137, kEnd = 250137;
        std::vector<std::atomic<int>> visits(kEnd);
        for (auto& v : visits) v.store(0);
        std::atomic<int> outOfRange{0};

        jobs.parallelFor(kBegin, kEnd, 97, [&](std::size_t b, std::size_t e) {
            if (b < kBegin || e > kEnd || b >= e) outOfRange.fetch_add(1);
            for (std::size_t i = b; i < e; ++i) visits[i].fetch_add(1, std::memory_order_relaxed);
        });

        int wrong = 0;
        for (std::size_t i = 0; i < kEnd; ++i) {
            const int expected = i >= kBegin ? 1 : 0;
            if (visits[i].load() != expected) ++wrong;
        }
        expectEqual("parallelFor tiles its range, workers=" + std::to_string(workers), wrong, 0);
        expectEqual("no chunk left the requested bounds, workers=" + std::to_string(workers),
                    outOfRange.load(), 0);
    }
}

// A job that itself waits must not stall the pool. With a blocking wait this
// deadlocks as soon as the nesting depth reaches the worker count; with helping
// it completes. Run it with deliberately few workers so the distinction bites.
void testNestedParallelismDoesNotDeadlock() {
    for (unsigned workers : {0u, 1u, 2u, 8u}) {
        JobSystem jobs(workers);
        std::atomic<long long> total{0};

        jobs.parallelFor(0, 64, 1, [&](std::size_t, std::size_t) {
            jobs.parallelFor(0, 64, 1, [&](std::size_t, std::size_t) {
                jobs.parallelFor(0, 16, 1, [&](std::size_t b, std::size_t e) {
                    total.fetch_add(static_cast<long long>(e - b), std::memory_order_relaxed);
                });
            });
        });

        expectEqual("three-deep nesting completes, workers=" + std::to_string(workers),
                    total.load(), 64LL * 64LL * 16LL);
    }
}

// The deque is bounded. Submitting far past its capacity must still execute
// everything -- the overflow path runs jobs inline rather than dropping them.
void testDequeOverflowRunsEverything() {
    JobSystem jobs(4);
    constexpr int kJobs = 40000;  // ~10x the per-lane deque capacity
    std::atomic<int> ran{0};

    JobSystem::Counter counter;
    std::atomic<int>* target = &ran;
    for (int i = 0; i < kJobs; ++i)
        jobs.run(counter, [target] { target->fetch_add(1, std::memory_order_relaxed); });
    jobs.wait(counter);

    expectEqual("every job ran despite deque overflow", ran.load(), kJobs);
}

// The reduction has to be bit-identical across worker counts, or replays and
// multiplayer cannot be reproducible. Values are chosen so that summation order
// genuinely matters -- see the vacuity check below.
void testDeterministicReduction() {
    constexpr std::size_t kCount = 400000;
    constexpr std::size_t kGrain = 1000;
    constexpr std::size_t kChunks = kCount / kGrain;

    // Chunk c is 1000 copies of 2^-(c mod 64), so the partials span 19 orders of
    // magnitude. Added to a large accumulator the smallest ones fall below its
    // ULP and vanish; accumulated among themselves first they survive. That makes
    // the fold order genuinely observable, which is what gives this test teeth --
    // an earlier version used values so dominated by a few huge terms that every
    // order produced the same answer, and it would have passed against a
    // completely non-deterministic reduction.
    auto value = [](std::size_t i) {
        return std::ldexp(1.0, -static_cast<int>((i / kGrain) % 64));
    };
    auto sumRange = [&](std::size_t b, std::size_t e) {
        double s = 0.0;
        for (std::size_t i = b; i < e; ++i) s += value(i);
        return s;
    };
    auto add = [](double a, double b) { return a + b; };

    std::vector<double> partials(kChunks);
    for (std::size_t c = 0; c < kChunks; ++c)
        partials[c] = sumRange(c * kGrain, (c + 1) * kGrain);

    // The reference is the ordered fold -- the definition parallelReduce must
    // match exactly, not approximate.
    double forward = 0.0;
    for (std::size_t c = 0; c < kChunks; ++c) forward = add(forward, partials[c]);
    double reverse = 0.0;
    for (std::size_t c = kChunks; c-- > 0;) reverse = add(reverse, partials[c]);

    expectTrue("the reduction test is not vacuous: fold order changes the result",
               forward != reverse);

    for (unsigned workers : {0u, 1u, 2u, 7u, 15u, 23u}) {
        JobSystem jobs(workers);
        const double got = jobs.parallelReduce(0, kCount, kGrain, 0.0, sumRange, add);
        // Bit-identical, not approximately equal.
        expectTrue("reduction is bit-identical to the ordered fold, workers=" +
                       std::to_string(workers),
                   got == forward);
    }
}

// Work must actually reach the workers. This is the one property that cannot be
// asserted as an invariant -- a heavily loaded machine could legitimately run
// everything on the submitting lane -- so it is deliberately lenient and only
// asks that stealing happened at all.
void testWorkReachesWorkers() {
    constexpr unsigned kWorkers = 8;
    JobSystem jobs(kWorkers);
    std::vector<std::atomic<int>> perLane(jobs.laneCount());
    for (auto& c : perLane) c.store(0);

    jobs.parallelFor(0, 200000, 64, [&](std::size_t b, std::size_t e) {
        // A little arithmetic so a chunk is not so short that the submitting
        // thread drains the queue before any worker wakes.
        volatile double sink = 0;
        for (std::size_t i = b; i < e; ++i) sink += static_cast<double>(i) * 1.000001;
        (void)sink;
        perLane[jobs.currentLane()].fetch_add(1, std::memory_order_relaxed);
    });

    int lanesUsed = 0;
    for (auto& c : perLane)
        if (c.load() > 0) ++lanesUsed;
    expectTrue("work was distributed across more than one lane", lanesUsed >= 2);
}

// Zero workers must still be fully functional: it is the debugging and
// determinism-checking configuration, and every test above runs it.
void testSingleThreadedModeWorks() {
    JobSystem jobs(0);
    expectEqual("zero-worker system reports no workers", jobs.workerCount(), 0);

    std::atomic<int> ran{0};
    std::atomic<int>* target = &ran;
    JobSystem::Counter counter;
    jobs.run(counter, [target] { target->fetch_add(1); });
    jobs.wait(counter);
    expectEqual("single-threaded job still runs", ran.load(), 1);
}

// --- grain auto-tuning -------------------------------------------------------

// The tuner must still tile the range exactly, probe included. A probe that
// forgets to hand its elements to the body, or hands them over twice, is the
// obvious way to get this wrong.
void testAutoGrainCoversRangeExactlyOnce() {
    for (unsigned workers : {0u, 4u, 15u}) {
        JobSystem jobs(workers);
        constexpr std::size_t kBegin = 61, kEnd = 400061;
        std::vector<std::atomic<int>> visits(kEnd);
        for (auto& v : visits) v.store(0);

        jobs.parallelForAuto(kBegin, kEnd, [&](std::size_t b, std::size_t e) {
            for (std::size_t i = b; i < e; ++i) visits[i].fetch_add(1, std::memory_order_relaxed);
        });

        int wrong = 0;
        for (std::size_t i = 0; i < kEnd; ++i)
            if (visits[i].load() != (i >= kBegin ? 1 : 0)) ++wrong;
        expectEqual("auto grain tiles its range, workers=" + std::to_string(workers), wrong, 0);
    }
}

// Burn a controlled amount of time per element so the tuner has something real
// to measure. Volatile stops it being optimised into nothing.
inline void burn(int iterations) {
    volatile double x = 1.0;
    for (int i = 0; i < iterations; ++i) x = x * 1.0000001 + 0.5;
}

// The point of a tuner is that it responds to cost. Expensive elements must get
// a smaller grain than cheap ones -- if they do not, it is not tuning, it is
// just picking a constant.
void testAutoGrainRespondsToElementCost() {
    JobSystem jobs(8);
    constexpr std::size_t kCount = 200000;

    const auto cheap = jobs.parallelForAuto(0, kCount, [](std::size_t b, std::size_t e) {
        for (std::size_t i = b; i < e; ++i) burn(1);
    });
    const auto expensive = jobs.parallelForAuto(0, kCount / 20, [](std::size_t b, std::size_t e) {
        for (std::size_t i = b; i < e; ++i) burn(200);
    });

    expectTrue("cheap elements were measured as cheaper than expensive ones",
               cheap.nsPerElement < expensive.nsPerElement);
    expectTrue("expensive elements get a smaller grain than cheap ones",
               expensive.grain < cheap.grain);
    expectTrue("both runs dispatched rather than falling back to serial",
               !cheap.ranSerially && !expensive.ranSerially);
}

// Chunk count must land inside the band the benchmark justified: enough chunks
// per lane to balance the tail, few enough that dispatch stays negligible.
void testAutoGrainStaysInsideItsBand() {
    for (unsigned workers : {2u, 8u, 15u}) {
        JobSystem jobs(workers);
        const auto lanes = static_cast<std::size_t>(jobs.laneCount());
        const auto result = jobs.parallelForAuto(0, 400000, [](std::size_t b, std::size_t e) {
            for (std::size_t i = b; i < e; ++i) burn(4);
        });

        const std::string suffix = ", workers=" + std::to_string(workers);
        expectTrue("auto grain dispatched" + suffix, !result.ranSerially && result.chunks > 0);
        expectTrue("at least kMinChunksPerLane chunks per lane" + suffix,
                   result.chunks >= lanes * JobSystem::kMinChunksPerLane);
        expectTrue("at most kMaxChunksPerLane chunks per lane" + suffix,
                   result.chunks <= lanes * JobSystem::kMaxChunksPerLane);
        expectTrue("chunks times grain covers the remaining range" + suffix,
                   result.chunks * result.grain >= 400000 - result.probeElements);
    }
}

// Work too small to be worth dispatching must be run on the spot. Paying a
// dispatch round trip to parallelise a microsecond of work is a net loss.
void testAutoGrainRunsTinyWorkSerially() {
    JobSystem jobs(8);
    std::atomic<int> visited{0};
    const auto result = jobs.parallelForAuto(0, 32, [&](std::size_t b, std::size_t e) {
        visited.fetch_add(static_cast<int>(e - b), std::memory_order_relaxed);
    });
    expectTrue("a trivial range is run serially", result.ranSerially);
    expectEqual("a serial run still visits every element", visited.load(), 32);
    expectEqual("an empty range does nothing",
                static_cast<long long>(jobs.parallelForAuto(5, 5, [](std::size_t, std::size_t) {
                                            }).probeElements),
                0);
}

}  // namespace

void runJobTests() {
    std::printf("\n--- job system ---\n");
    testSingleThreadedModeWorks();
    testEveryJobRunsExactlyOnce();
    testParallelForCoversRangeExactlyOnce();
    testNestedParallelismDoesNotDeadlock();
    testDequeOverflowRunsEverything();
    testDeterministicReduction();
    testWorkReachesWorkers();
    testAutoGrainCoversRangeExactlyOnce();
    testAutoGrainRespondsToElementCost();
    testAutoGrainStaysInsideItsBand();
    testAutoGrainRunsTinyWorkSerially();
}
