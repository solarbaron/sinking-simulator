// SPDX-License-Identifier: MIT
//
// Job system throughput and grain scaling.
//
// This exists to settle a question that was previously folklore: is the queue
// structure worth optimising, or is chunk granularity the thing that actually
// decides parallel efficiency? docs/06-roadmap.md makes the Chase-Lev revisit
// explicitly conditional on what this measures.
//
// The number that matters is **per-job overhead**. Everything else follows: if
// dispatching a job costs X, then a chunk must take much longer than X for the
// parallelism to be worth having, and the recommended grain is whatever makes
// that true.
#include "engine/core/jobs.hpp"

#include <algorithm>
#include <limits>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <string>
#include <thread>
#include <vector>

using core::JobSystem;

namespace {

int failures = 0;

void check(const char* what, bool ok, const std::string& detail = {}) {
    std::printf("  [%s] %-56s %s\n", ok ? "PASS" : "FAIL", what, detail.c_str());
    if (!ok) ++failures;
}

using Clock = std::chrono::steady_clock;

double millisSince(Clock::time_point start) {
    return std::chrono::duration<double, std::milli>(Clock::now() - start).count();
}

// Deterministic, memory-free arithmetic so a chunk's cost is predictable and
// nothing can be optimised away or turned into a cache-miss benchmark by
// accident. The dependent chain stops the compiler vectorising it into nothing.
inline double elementWork(std::size_t i) {
    double x = static_cast<double>(i) * 1e-6;
    for (int k = 0; k < 8; ++k) x = x * 1.0000001 + 0.5;
    return x;
}

double sumRange(std::size_t begin, std::size_t end) {
    double sum = 0.0;
    for (std::size_t i = begin; i < end; ++i) sum += elementWork(i);
    return sum;
}

// Every measured result must reach an observable location or the optimiser
// deletes the loop and the benchmark silently measures dispatch overhead alone --
// which is exactly what happened on the first run of this file.
volatile double g_sink = 0.0;

// One accumulator per lane, each on its own cache line. Accumulating per chunk
// into shared memory would make the fine-grain rows measure false sharing rather
// than the queue.
struct alignas(64) LaneAccumulator {
    double value = 0.0;
    char padding[56]{};
};

// Minimum and maximum across several runs. The development machine is not quiet
// -- see the load average printed below -- and the minimum is the best estimate
// of the cost without interference, where the mean would mostly measure the
// neighbours. The maximum is reported alongside so the reader can see how much
// of the number to believe: on the first version of this benchmark the
// single-threaded baseline varied between 56 ms and 99 ms across runs, which
// made every derived speedup figure meaningless while looking perfectly precise.
struct Timing {
    double best = 1e300;
    double worst = 0.0;
    double spread() const { return worst / best; }
};

template <typename F>
Timing timeRepeats(int repeats, const F& run) {
    Timing t;
    for (int r = 0; r < repeats; ++r) {
        const auto start = Clock::now();
        run();
        const double ms = millisSince(start);
        t.best = std::min(t.best, ms);
        t.worst = std::max(t.worst, ms);
    }
    return t;
}

template <typename F>
double bestOf(int repeats, const F& run) {
    return timeRepeats(repeats, run).best;
}

double loadAverage() {
    double one = 0;
    if (std::FILE* f = std::fopen("/proc/loadavg", "r")) {
        if (std::fscanf(f, "%lf", &one) != 1) one = 0;
        std::fclose(f);
    }
    return one;
}

// --- 1. Per-job dispatch overhead -------------------------------------------

double measureDispatchOverhead() {
    std::printf("\n1. Per-job dispatch overhead (empty jobs)\n");
    std::printf("     %8s %14s %16s\n", "workers", "ns/job", "M jobs/s");

    constexpr int kJobs = 200000;
    double singleLaneOverheadNs = 0;

    for (unsigned workers : {0u, 1u, 2u, 4u, 8u, 16u, 23u}) {
        JobSystem jobs(workers);
        std::atomic<int> sink{0};
        std::atomic<int>* target = &sink;

        const double ms = bestOf(3, [&] {
            JobSystem::Counter counter;
            for (int i = 0; i < kJobs; ++i)
                jobs.run(counter, [target] { target->fetch_add(1, std::memory_order_relaxed); });
            jobs.wait(counter);
        });

        const double nsPerJob = ms * 1e6 / kJobs;
        std::printf("     %8u %14.1f %16.2f\n", workers, nsPerJob, kJobs / (ms * 1e3));
        if (workers == 0) singleLaneOverheadNs = nsPerJob;
    }
    return singleLaneOverheadNs;
}

// --- 2. Grain sweep ----------------------------------------------------------

// Returns the span between the worst and best grain at this worker count -- the
// figure `docs/01-architecture.md` §2 publishes a bound on. It used to be printed
// in the verdict as a hard-coded "roughly 30x", which was neither measured nor the
// value the document carries: the doc retired the single-figure form outright
// ("a single figure hides that: ~40x was right for 23 workers and about double
// the truth at 8") and publishes a bound of at least 15x with medians of 21x and
// 44x. A literal in a printf is a figure nobody re-derives, in the one program
// that is in a position to.
double measureGrainSweep(unsigned workers, std::size_t elements, double baselineMs) {
    std::printf("\n   %u workers:\n", workers);
    std::printf("     %10s %12s %12s %10s %12s\n", "grain", "chunks", "chunk us", "ms", "speedup");

    JobSystem jobs(workers);
    std::vector<LaneAccumulator> lanes(jobs.laneCount());
    double best = 0.0, worst = std::numeric_limits<double>::infinity();

    for (std::size_t grain : {std::size_t{16}, std::size_t{64}, std::size_t{256}, std::size_t{1024},
                              std::size_t{4096}, std::size_t{16384}, std::size_t{65536},
                              std::size_t{262144}}) {
        const double ms = bestOf(3, [&] {
            jobs.parallelFor(0, elements, grain, [&](std::size_t b, std::size_t e) {
                lanes[jobs.currentLane()].value += sumRange(b, e);
            });
        });
        g_sink = lanes[0].value;

        const std::size_t chunks = (elements + grain - 1) / grain;
        const double chunkMicros = baselineMs * 1e3 / static_cast<double>(chunks);
        const double speedup = baselineMs / ms;
        best = std::max(best, speedup);
        worst = std::min(worst, speedup);
        std::printf("     %10zu %12zu %12.2f %10.2f %11.2fx\n", grain, chunks, chunkMicros, ms,
                    speedup);
    }
    // A ratio between two timings from the same sweep, which is what survives a
    // busy box: both ends move together. The absolute ms columns do not.
    return best / worst;
}

// --- 3. Worker scaling at a sensible grain -----------------------------------

void measureWorkerScaling(std::size_t elements, std::size_t grain, double baselineMs) {
    std::printf("\n3. Worker scaling at grain %zu\n", grain);
    std::printf("     %8s %10s %10s %12s\n", "workers", "ms", "speedup", "efficiency");

    for (unsigned workers : {1u, 2u, 4u, 8u, 16u, 23u}) {
        JobSystem jobs(workers);
        std::vector<LaneAccumulator> lanes(jobs.laneCount());
        const double ms = bestOf(3, [&] {
            jobs.parallelFor(0, elements, grain, [&](std::size_t b, std::size_t e) {
                lanes[jobs.currentLane()].value += sumRange(b, e);
            });
        });
        g_sink = lanes[0].value;
        // The submitting thread works too, so the pool is workers + 1 wide.
        const double poolWidth = workers + 1.0;
        const double speedup = baselineMs / ms;
        std::printf("     %8u %10.2f %9.2fx %11.0f%%\n", workers, ms, speedup,
                    100.0 * speedup / poolWidth);
    }
}

// --- 3b. Does the auto-tuner land in the plateau? ----------------------------

void measureAutoGrain(std::size_t elements, double baselineMs) {
    std::printf("\n3b. Auto-tuned grain vs the sweep\n");
    std::printf("     %8s %12s %12s %10s %12s\n", "workers", "grain", "chunks", "ms", "speedup");

    for (unsigned workers : {8u, 23u}) {
        JobSystem jobs(workers);
        std::vector<LaneAccumulator> lanes(jobs.laneCount());
        JobSystem::AutoGrain chosen;
        const double ms = bestOf(3, [&] {
            chosen = jobs.parallelForAuto(0, elements, [&](std::size_t b, std::size_t e) {
                lanes[jobs.currentLane()].value += sumRange(b, e);
            });
        });
        g_sink = lanes[0].value;
        std::printf("     %8u %12zu %12zu %10.2f %11.2fx\n", workers, chosen.grain, chosen.chunks,
                    ms, baselineMs / ms);
    }
    std::printf("     Compare against the sweep above: the tuner should land on the flat\n");
    std::printf("     part of the curve without having been told the element cost.\n");
}

// --- 4. Correctness under all of the above -----------------------------------

void checkResultsAreStable(std::size_t elements, std::size_t grain) {
    std::printf("\n4. Correctness across worker counts\n");
    auto add = [](double a, double b) { return a + b; };

    double reference = 0.0;
    for (std::size_t b = 0; b < elements; b += grain)
        reference = add(reference, sumRange(b, std::min(b + grain, elements)));

    bool identical = true;
    for (unsigned workers : {0u, 1u, 4u, 23u}) {
        JobSystem jobs(workers);
        const double got = jobs.parallelReduce(0, elements, grain, 0.0, sumRange, add);
        if (got != reference) identical = false;
    }
    check("parallelReduce is bit-identical across 0, 1, 4 and 23 workers", identical);
}

}  // namespace

int main() {
    const unsigned hardware = std::thread::hardware_concurrency();
    std::printf("shipsim - job system throughput and grain scaling\n");
    std::printf("  %u hardware threads, load average %.2f\n", hardware, loadAverage());
    std::printf("  (a loaded machine depresses scaling; timings are the best of 3 runs)\n");

    const double overheadNs = measureDispatchOverhead();

    constexpr std::size_t kElements = 20'000'000;
    std::printf("\n2. parallelFor grain sweep over %zu elements\n", kElements);
    // More repeats here than anywhere else: every speedup column is divided by
    // this one number, so its noise propagates into the entire table.
    const Timing baseline = timeRepeats(9, [&] { g_sink = sumRange(0, kElements); });
    const double baselineMs = baseline.best;
    std::printf("   single-threaded baseline: %.2f ms best, %.2f ms worst (%.2fx spread),"
                " %.1f ns/element\n",
                baseline.best, baseline.worst, baseline.spread(), baselineMs * 1e6 / kElements);
    if (baseline.spread() > 1.25)
        std::printf("   NOTE: baseline spread exceeds 25%%; treat the absolute ms columns as the\n"
                    "         reliable measurement and the speedup columns as indicative only.\n");

    const double span8 = measureGrainSweep(8, kElements, baselineMs);
    const double span23 = measureGrainSweep(23, kElements, baselineMs);

    // Grain that puts a chunk near the 50 us target from docs/01-architecture.md.
    const double nsPerElement = baselineMs * 1e6 / kElements;
    const auto targetGrain = static_cast<std::size_t>(50000.0 / nsPerElement);
    std::printf("\n   grain for a 50 us chunk at %.1f ns/element: %zu elements\n", nsPerElement,
                targetGrain);

    measureWorkerScaling(kElements, targetGrain, baselineMs);
    measureAutoGrain(kElements, baselineMs);
    checkResultsAreStable(kElements, targetGrain);

    // What this all means for the Chase-Lev question.
    const double chunkMicros = 50.0;
    const double overheadFraction = (overheadNs * 1e-3) / chunkMicros;
    std::printf("\n5. Verdict\n");
    std::printf("     uncontended dispatch     %.0f ns/job\n", overheadNs);
    std::printf("     at a %.0f us chunk that is %.3f%% of the work\n", chunkMicros,
                100.0 * overheadFraction);
    std::printf("     Chase-Lev could remove only part of that, so it cannot matter here.\n");
    std::printf("     Grain is what matters. Worst grain against best, this run:\n");
    // One line, both numbers, fixed columns: `check-figures.sh` re-derives the same
    // span off the printed sweep table and compares the two, so this has to be
    // parseable without a multi-line match. A figure the gate cannot read is a
    // figure back in the class this whole exercise is about.
    std::printf("     grain span: %.2fx at 8 workers, %.2fx at 23\n", span8, span23);
    std::printf("     Both move run to run -- the grain-16 row alone varies 3x at 23\n");
    std::printf("     workers -- so the published claim is the bound over many runs and\n");
    std::printf("     not either figure here.\n");

    // A conservative floor: dispatch must not be so slow that even coarse chunks
    // are dominated by it. This would fail on a genuinely broken queue.
    check("dispatch overhead is below 5% of a 50 us chunk", overheadFraction < 0.05,
          std::to_string(static_cast<int>(overheadNs)) + " ns/job");

    // The bound the document publishes, asserted where it is measured. 15x is not
    // a round number chosen for comfort: it was 17x on sixteen runs, survived a
    // further twenty-five at a minimum of 17.80x, and was contradicted once at
    // 16.51x on a box at load 0.47. It sits below that contradiction deliberately.
    //
    // It is also safe in the right direction. The load dependence has a *sign*:
    // grain-16 is dispatch-bound and the plateau is work-bound, so an idler box
    // makes the numerator cheaper faster than the denominator and the ratio
    // **falls**. Ten spinners moved it up to 20.6-28.1x. So contention cannot
    // falsely fail this check -- only an idle box can approach it, and this is the
    // one bound in the file for which that is the safe direction.
    const double span = std::min(span8, span23);
    check("the grain penalty is at least the 15x docs/01-architecture.md publishes",
          span >= 15.0, std::to_string(span) + "x");

    std::printf("\n%s\n", failures == 0 ? "all checks passed" : "CHECKS FAILED");
    return failures == 0 ? 0 : 1;
}
