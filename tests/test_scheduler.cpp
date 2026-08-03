// SPDX-License-Identifier: MIT
//
// Validation of the multi-rate scheduler.
//
// A scheduler that is subtly wrong does not crash; it runs a solver slightly too
// often or slightly too seldom and the simulation is quietly off. So the checks
// here are closed-form: a system at R Hz must run exactly R times per simulated
// second, a 10:1 rate ratio must produce exactly a 10:1 update ratio, and a
// given sequence of advances must always produce the same schedule.
#include "engine/core/jobs.hpp"
#include "engine/core/scheduler.hpp"
#include "harness.hpp"

#include <atomic>
#include <cstdio>
#include <utility>
#include <string>
#include <vector>

using core::JobSystem;
using core::Scheduler;
using core::SystemConfig;
using testing::expectEqual;
using testing::expectTrue;

namespace {

SystemConfig config(std::string name, double rateHz) {
    SystemConfig c;
    c.name = std::move(name);
    c.rateHz = rateHz;
    c.maxCatchUpSteps = 1000000;  // catch-up limits are exercised separately
    return c;
}

// A system at R Hz must run exactly R times per simulated second. Not "about" --
// exactly, because the count comes from accumulated simulation time and nothing
// else.
void testRateIsExact() {
    Scheduler scheduler;
    std::uint64_t hundred = 0, ten = 0, one = 0;
    const auto a = scheduler.add(config("hydrostatics", 100.0), [&](double) { ++hundred; });
    const auto b = scheduler.add(config("thermal", 10.0), [&](double) { ++ten; });
    const auto c = scheduler.add(config("crew", 1.0), [&](double) { ++one; });

    // Ten seconds of wall time at 1x, in 1/60 s slices.
    for (int i = 0; i < 600; ++i) scheduler.advance(1.0 / 60.0);

    // Closed form, in the scheduler's own integer units: a system runs exactly
    // (elapsed nanoseconds / its period) times. Stating it in seconds would
    // reintroduce the floating-point rounding the integer clock exists to avoid.
    const std::int64_t elapsed = scheduler.simulationNanos();
    expectEqual("100 Hz count equals elapsed nanos over its period",
                static_cast<long long>(hundred),
                static_cast<long long>(elapsed / scheduler.periodNanos(a)));
    expectEqual("10 Hz count equals elapsed nanos over its period",
                static_cast<long long>(ten),
                static_cast<long long>(elapsed / scheduler.periodNanos(b)));
    expectEqual("1 Hz count equals elapsed nanos over its period",
                static_cast<long long>(one),
                static_cast<long long>(elapsed / scheduler.periodNanos(c)));
    expectEqual("the scheduler agrees with the counter",
                static_cast<long long>(scheduler.updateCount(a)),
                static_cast<long long>(hundred));
    expectEqual("100 Hz over ten seconds is a thousand updates",
                static_cast<long long>(hundred), 1000);
    expectTrue("simulation time tracks wall time at 1x",
               std::abs(scheduler.simulationTime() - 10.0) < 1e-6);
}

// The whole point of multi-rate: the ratio between two systems' update counts
// must be the ratio of their rates, independently of the advance slicing.
void testRateRatiosHoldUnderIrregularSlicing() {
    Scheduler scheduler;
    std::uint64_t fast = 0, slow = 0;
    const auto fastId = scheduler.add(config("fast", 200.0), [&](double) { ++fast; });
    const auto slowId = scheduler.add(config("slow", 20.0), [&](double) { ++slow; });

    // Deliberately ragged frame times, as a real frame loop produces.
    const double slices[] = {0.011, 0.007, 0.033, 0.004, 0.021, 0.016};
    for (int repeat = 0; repeat < 100; ++repeat)
        for (double slice : slices) scheduler.advance(slice);

    const std::int64_t elapsed = scheduler.simulationNanos();
    expectEqual("fast count equals elapsed nanos over its period",
                static_cast<long long>(fast),
                static_cast<long long>(elapsed / scheduler.periodNanos(fastId)));
    expectEqual("slow count equals elapsed nanos over its period",
                static_cast<long long>(slow),
                static_cast<long long>(elapsed / scheduler.periodNanos(slowId)));
    // Ragged slices must not cost updates: 9.2 s of simulation at 200 Hz is
    // 1840 updates however it was sliced.
    expectEqual("ragged slicing loses no updates", static_cast<long long>(fast), 1840);
    expectEqual("the 10:1 rate ratio holds exactly",
                static_cast<long long>(fast), static_cast<long long>(slow) * 10);
}

// Dilation scales simulation time against wall time. Doubling it must double the
// updates per wall second; halving it must halve them.
void testTimeDilationScalesUpdates() {
    for (double dilation : {0.25, 0.5, 1.0, 4.0}) {
        Scheduler scheduler;
        std::uint64_t updates = 0;
        scheduler.add(config("sim", 100.0), [&](double) { ++updates; });
        scheduler.setTimeDilation(dilation);

        for (int i = 0; i < 100; ++i) scheduler.advance(0.01);  // 1 s of wall time

        const auto expected = static_cast<long long>(100.0 * dilation);
        expectEqual("dilation " + std::to_string(dilation) + " scales update count",
                    static_cast<long long>(updates), expected);
        expectTrue("dilation " + std::to_string(dilation) + " scales simulation time",
                   std::abs(scheduler.simulationTime() - dilation) < 1e-6);
    }
}

// A system that cannot honour its rate at a given dilation must be skipped, not
// run late. This is how a fast-forward drops the expensive solvers.
void testDilationBandGatesSystems() {
    Scheduler scheduler;
    std::uint64_t always = 0, onlyRealTime = 0, onlySlowMotion = 0;

    scheduler.add(config("hydrostatics", 100.0), [&](double) { ++always; });

    SystemConfig fem = config("fem_zone", 100.0);
    fem.maxDilation = 2.0;  // cannot keep up when time runs faster than 2x
    scheduler.add(fem, [&](double) { ++onlyRealTime; });

    SystemConfig fracture = config("fracture", 100.0);
    fracture.maxDilation = 0.1;  // only affordable in heavy slow motion
    scheduler.add(fracture, [&](double) { ++onlySlowMotion; });

    scheduler.setTimeDilation(1.0);
    for (int i = 0; i < 100; ++i) scheduler.advance(0.01);
    expectTrue("at 1x the always-on system runs", always > 0);
    expectTrue("at 1x the real-time-only system runs", onlyRealTime > 0);
    expectEqual("at 1x the slow-motion-only system is skipped",
                static_cast<long long>(onlySlowMotion), 0);

    const std::uint64_t femBefore = onlyRealTime;
    scheduler.setTimeDilation(60.0);
    for (int i = 0; i < 100; ++i) scheduler.advance(0.01);
    expectTrue("fast-forward keeps the cheap system running", always > 100);
    expectEqual("fast-forward drops the system that cannot honour its rate",
                static_cast<long long>(onlyRealTime), static_cast<long long>(femBefore));

    scheduler.setTimeDilation(0.02);
    for (int i = 0; i < 100; ++i) scheduler.advance(0.01);
    expectTrue("slow motion enables the expensive system", onlySlowMotion > 0);
}

// One long frame must not be able to demand unbounded catch-up, or the next
// frame is longer still and the loop never recovers.
void testCatchUpCeilingPreventsSpiral() {
    Scheduler scheduler;
    std::uint64_t updates = 0;
    SystemConfig c = config("sim", 1000.0);
    c.maxCatchUpSteps = 8;
    scheduler.add(c, [&](double) { ++updates; });

    // A ten second hitch would otherwise ask for 10 000 updates at once.
    const Scheduler::Report report = scheduler.advance(10.0);
    expectEqual("catch-up is capped at the configured ceiling",
                static_cast<long long>(updates), 8);
    expectTrue("the drop is reported rather than hidden", report.droppedCatchUp > 9000);

    // And the surplus must not be carried into the next advance.
    updates = 0;
    scheduler.advance(0.001);
    expectTrue("the discarded surplus does not resurface later", updates <= 8);
}

// Dependencies must be honoured, and a cycle must be refused rather than
// producing an arbitrary order.
void testDependencyOrdering() {
    Scheduler scheduler;
    std::vector<int> order;

    const auto flooding = scheduler.add(config("flooding", 100.0), [&](double) { order.push_back(0); });
    const auto rigidBody = scheduler.add(config("rigid_body", 100.0), [&](double) { order.push_back(1); });
    const auto render = scheduler.add(config("render_state", 100.0), [&](double) { order.push_back(2); });

    expectTrue("rigid body depends on flooding", scheduler.dependsOn(rigidBody, flooding));
    expectTrue("render state depends on rigid body", scheduler.dependsOn(render, rigidBody));
    expectTrue("a self-edge is refused", !scheduler.dependsOn(flooding, flooding));
    expectTrue("a cycle is refused", !scheduler.dependsOn(flooding, render));

    scheduler.advance(0.01);
    expectEqual("every system ran once", static_cast<long long>(order.size()), 3);
    expectTrue("execution followed the dependency order",
               order.size() == 3 && order[0] == 0 && order[1] == 1 && order[2] == 2);

    const auto& levels = scheduler.levels();
    expectEqual("a chain of three produces three levels",
                static_cast<long long>(levels.size()), 3);
}

// Independent systems share a level and may run in parallel; the results must
// not depend on whether they did.
void testParallelExecutionMatchesSerial() {
    auto run = [](JobSystem* jobs) {
        Scheduler scheduler(jobs);
        std::vector<std::atomic<std::uint64_t>> counters(16);
        for (auto& counter : counters) counter.store(0);
        for (int i = 0; i < 16; ++i)
            scheduler.add(config("system" + std::to_string(i), 50.0 + i),
                          [&counters, i](double) {
                              counters[static_cast<std::size_t>(i)].fetch_add(
                                  1, std::memory_order_relaxed);
                          });
        for (int i = 0; i < 200; ++i) scheduler.advance(0.01);

        std::vector<std::uint64_t> result;
        for (auto& counter : counters) result.push_back(counter.load());
        // Closed form alongside the counts, computed from the scheduler's own
        // integer clock so the comparison is exact.
        std::vector<std::uint64_t> expected;
        for (int i = 0; i < 16; ++i)
            expected.push_back(static_cast<std::uint64_t>(
                scheduler.simulationNanos() / scheduler.periodNanos(static_cast<unsigned>(i))));
        return std::pair{result, expected};
    };

    JobSystem jobs(8);
    const auto serial = run(nullptr);
    const auto parallel = run(&jobs);
    expectTrue("parallel execution produces identical update counts",
               serial.first == parallel.first);
    expectTrue("each system ran elapsed nanos over its period, serially",
               serial.first == serial.second);
    expectTrue("each system ran elapsed nanos over its period, in parallel",
               parallel.first == parallel.second);
}

// The same sequence of advances must always produce the same schedule, or
// replays diverge.
void testScheduleIsReproducible() {
    auto run = [] {
        Scheduler scheduler;
        std::vector<int> trace;
        const auto a = scheduler.add(config("a", 30.0), [&](double) { trace.push_back(0); });
        const auto b = scheduler.add(config("b", 70.0), [&](double) { trace.push_back(1); });
        scheduler.dependsOn(b, a);
        scheduler.setTimeDilation(1.5);

        const double slices[] = {0.013, 0.009, 0.027, 0.004};
        for (int repeat = 0; repeat < 50; ++repeat)
            for (double slice : slices) scheduler.advance(slice);
        return trace;
    };

    const auto first = run();
    const auto second = run();
    expectTrue("the schedule is non-trivial", first.size() > 100);
    expectTrue("identical inputs produce an identical schedule", first == second);
}

// A zero or negative advance must be a no-op rather than running anything or
// moving simulation time backwards.
void testDegenerateAdvances() {
    Scheduler scheduler;
    std::uint64_t updates = 0;
    scheduler.add(config("sim", 100.0), [&](double) { ++updates; });

    scheduler.advance(0.0);
    scheduler.advance(-1.0);
    expectEqual("a zero or negative advance runs nothing", static_cast<long long>(updates), 0);
    expectTrue("simulation time did not move", scheduler.simulationTime() == 0.0);

    // A paused simulation is dilation zero: wall time passes, sim time does not.
    scheduler.setTimeDilation(0.0);
    for (int i = 0; i < 100; ++i) scheduler.advance(0.01);
    expectEqual("dilation zero pauses the simulation", static_cast<long long>(updates), 0);
    expectTrue("paused simulation time stays put", scheduler.simulationTime() == 0.0);
}

}  // namespace

void runSchedulerTests() {
    std::printf("\n--- multi-rate scheduler ---\n");
    testRateIsExact();
    testRateRatiosHoldUnderIrregularSlicing();
    testTimeDilationScalesUpdates();
    testDilationBandGatesSystems();
    testCatchUpCeilingPreventsSpiral();
    testDependencyOrdering();
    testParallelExecutionMatchesSerial();
    testScheduleIsReproducible();
    testDegenerateAdvances();
}
