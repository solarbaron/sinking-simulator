// SPDX-License-Identifier: MIT
//
// Work-stealing job system.
//
// The frame is a DAG of small jobs (target 50-500 us) rather than a set of
// long-lived per-system threads, because the simulation's load is wildly uneven:
// a ship with one active damage zone and fourteen quiet compartments cannot be
// balanced by static partitioning.
//
// Two properties matter more here than raw throughput:
//
//   * **Waiting never blocks a worker.** A thread that waits on a counter
//     executes other jobs while it waits. Without this, a nested parallelFor
//     inside a job would deadlock the pool.
//   * **Reductions are deterministic.** parallelReduce() bins per *chunk*, not
//     per lane, and folds the bins in chunk order. Floating-point addition is not
//     associative, so a reduction that depends on which worker finished first
//     would make replays and multiplayer irreproducible.
#pragma once

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <new>
#include <type_traits>
#include <vector>

namespace core {

class JobSystem {
public:
    // Tracks outstanding work for one batch of jobs. Not copyable, and it must
    // outlive every job submitted against it.
    class Counter {
    public:
        Counter() = default;
        Counter(const Counter&) = delete;
        Counter& operator=(const Counter&) = delete;
        bool complete() const { return outstanding_.load(std::memory_order_acquire) == 0; }

    private:
        friend class JobSystem;
        std::atomic<int> outstanding_{0};
    };

    // A job's payload is stored inline. Anything larger has to be referenced
    // through a pointer the caller keeps alive, which is the normal pattern for
    // parallelFor anyway.
    static constexpr std::size_t kPayloadBytes = 48;

    // Number of worker threads created when the caller does not specify:
    // hardware concurrency minus one, because the submitting thread participates.
    static unsigned defaultWorkerCount();

    explicit JobSystem(unsigned workerCount = defaultWorkerCount());
    ~JobSystem();
    JobSystem(const JobSystem&) = delete;
    JobSystem& operator=(const JobSystem&) = delete;

    unsigned workerCount() const { return workerCount_; }
    // Worker lanes plus a few spare lanes for external submitting threads.
    unsigned laneCount() const { return static_cast<unsigned>(lanes_.size()); }
    // Lane index of the calling thread, claiming one on first use.
    unsigned currentLane() const;

    // Submit one job. `body` is invoked as body(). It is copied into the job's
    // inline storage, so it must be trivially copyable and fit in kPayloadBytes.
    template <typename F>
    void run(Counter& counter, F&& body);

    // Submit [begin, end) split into chunks of at most `grain`. `body` is invoked
    // as body(chunkBegin, chunkEnd) and is referenced by pointer, so it must stay
    // alive until wait(counter) returns.
    template <typename F>
    void dispatch(Counter& counter, std::size_t begin, std::size_t end, std::size_t grain,
                  const F& body);

    // dispatch() followed by wait(). The common, safe form.
    template <typename F>
    void parallelFor(std::size_t begin, std::size_t end, std::size_t grain, const F& body);

    // What the auto-tuner decided. Returned rather than hidden so callers can
    // see it, and so tests can assert on it instead of inferring it.
    struct AutoGrain {
        std::size_t grain = 0;           // chosen chunk size, 0 if run serially
        std::size_t chunks = 0;          // dispatched chunks, excluding the probe
        std::size_t probeElements = 0;   // elements consumed measuring the cost
        double      nsPerElement = 0;    // measured cost
        bool        ranSerially = false; // too little work to be worth dispatching
    };

    // parallelFor with the grain chosen from a measurement rather than from the
    // caller's guess. `tools/job_bench` shows at least a 17x swing between a bad grain
    // and a good one, against ~0.2% from dispatch cost -- this is where the
    // leverage is.
    //
    // NOT for anything whose result depends on chunk boundaries. The grain here
    // comes from a wall-clock probe, so it varies run to run; a reduction chunked
    // this way would fold its partials in a different order each time. That is
    // why parallelReduce still requires an explicit grain and always will.
    template <typename F>
    AutoGrain parallelForAuto(std::size_t begin, std::size_t end, const F& body);

    // Deterministic parallel reduction. `body(chunkBegin, chunkEnd)` returns a
    // partial result per chunk; `combine(a, b)` folds two partials.
    //
    // Partials are stored per chunk and folded in chunk index order, never in
    // completion order, so the result is bit-identical regardless of worker count
    // or scheduling. That is the whole point -- see the class comment.
    template <typename T, typename Body, typename Combine>
    T parallelReduce(std::size_t begin, std::size_t end, std::size_t grain, const T& identity,
                     const Body& body, const Combine& combine);

    // Execute jobs until `counter` reaches zero. Helps rather than blocks.
    void wait(Counter& counter);

    // Tuning constants, all derived from tools/job_bench measurements.
    // Efficiency plateaus around 2 us chunks, so 10 us is comfortably inside the
    // flat region with margin for the probe being wrong.
    //
    // **On a large loop this target is never reached, and `kMaxChunksPerLane` is
    // why.** `parallelForAuto` clamps the cost-derived grain into
    // `[remaining/maxChunks, remaining/minChunks]`, and on `job_bench`'s 20 M
    // element sweep the lower bound binds every time: 12 lanes x 64 gives 768
    // chunks and 27 lanes x 64 gives 1728, so the grain comes out 26 040 and
    // 11 574 whatever the probe measured. At the measured 2.3 ns an element those
    // are **60 us and 27 us chunks against a 10 us target**. Measured across eight
    // runs the chosen grain was bit-identical every time while `nsPerElement`
    // moved, which is the signature of a binding clamp rather than a measurement.
    //
    // That is not obviously wrong -- both land inside the plateau, which is what
    // the target exists to guarantee -- but it means this constant does not decide
    // grain for any loop big enough to matter, and a change to it would move
    // nothing on that workload. The chunk-count ceiling is the real control.
    static constexpr double kTargetChunkNanos = 10'000.0;
    // Below this there is less work than a dispatch round trip is worth.
    static constexpr double kSerialThresholdNanos = 20'000.0;
    // Probe until the measurement is this long, so it is well clear of clock
    // resolution rather than measuring timer noise.
    static constexpr double kMinProbeNanos = 2'000.0;
    // At least this many chunks per lane, or the tail of the loop is one lane
    // working while the rest idle. The sweep shows balance degrading once chunks
    // get very large for exactly this reason.
    static constexpr std::size_t kMinChunksPerLane = 2;
    // And at most this many, so a cheap body does not drown in dispatch.
    static constexpr std::size_t kMaxChunksPerLane = 64;

private:
    struct Job {
        using Entry = void (*)(void* payload, std::size_t begin, std::size_t end);
        Entry entry = nullptr;
        Counter* counter = nullptr;
        std::size_t begin = 0;
        std::size_t end = 0;
        alignas(std::max_align_t) std::byte payload[kPayloadBytes]{};
    };

    // Bounded MPMC ring (Vyukov). Each cell carries a sequence number that acts
    // as its ownership token: a producer may write a cell only when it has won
    // the tail CAS *and* the cell's sequence equals its own ticket, and a
    // consumer may read it only when the sequence says a producer has published
    // it. That makes every cell access exclusive by construction.
    //
    // This replaced a Chase-Lev work-stealing deque, and the reason is worth
    // recording. Chase-Lev is faster -- the owner pushes and pops LIFO with no
    // atomic on the fast path -- but it is only safe when a slot is an atomic
    // that can be read speculatively and discarded. A job record is 80 bytes, so
    // the speculative read is a plain struct copy racing with the owner's write:
    // undefined behaviour, and ThreadSanitizer finds it in seconds. Claiming the
    // slot with the CAS first does not fix it either, because advancing `top`
    // does not reserve the slot -- once later thieves push `top` past it, the
    // owner may reuse it while the first thief has still not copied it out.
    //
    // Getting Chase-Lev back means storing atomic pointers plus an epoch or
    // hazard-pointer reclamation scheme for the records. That is worth doing when
    // job dispatch shows up in a profile; at 50-500 us per job, one CAS each way
    // does not.
    class Queue {
    public:
        static constexpr std::size_t kCapacity = 2048;  // power of two
        Queue();
        bool push(const Job& job);  // any thread
        bool pop(Job& out);         // any thread

    private:
        struct Cell {
            std::atomic<std::size_t> sequence;
            Job job;
        };
        static constexpr std::size_t kMask = kCapacity - 1;
        Cell cells_[kCapacity];
        alignas(64) std::atomic<std::size_t> head_{0};
        alignas(64) std::atomic<std::size_t> tail_{0};
    };

    struct alignas(64) Lane {
        Queue queue;
        std::uint64_t rng = 0x9E3779B97F4A7C15ull;
    };

    void submit(unsigned lane, const Job& job);
    // Finds one job: own deque first, then a randomised sweep of the others.
    bool acquireJob(unsigned lane, Job& out);
    void execute(Job& job);
    void workerMain(unsigned lane);
    void wakeWorkers(int count);

    unsigned workerCount_ = 0;
    std::vector<Lane> lanes_;
    struct Impl;
    Impl* impl_ = nullptr;
};

// --- inline template definitions --------------------------------------------

template <typename F>
void JobSystem::run(Counter& counter, F&& body) {
    using Body = std::decay_t<F>;
    static_assert(sizeof(Body) <= kPayloadBytes,
                  "job body too large for inline storage; capture by pointer instead");
    static_assert(std::is_trivially_copyable_v<Body>,
                  "job body must be trivially copyable; it is memcpy'd into the job record");
    static_assert(std::is_trivially_destructible_v<Body>,
                  "job body must be trivially destructible; jobs are never destructed");

    Job job;
    ::new (static_cast<void*>(job.payload)) Body(std::forward<F>(body));
    job.entry = [](void* payload, std::size_t, std::size_t) {
        (*static_cast<Body*>(payload))();
    };
    job.counter = &counter;
    counter.outstanding_.fetch_add(1, std::memory_order_relaxed);
    submit(currentLane(), job);
}

template <typename F>
void JobSystem::dispatch(Counter& counter, std::size_t begin, std::size_t end, std::size_t grain,
                         const F& body) {
    if (begin >= end) return;
    if (grain == 0) grain = 1;

    const unsigned lane = currentLane();
    const std::size_t chunks = (end - begin + grain - 1) / grain;
    counter.outstanding_.fetch_add(static_cast<int>(chunks), std::memory_order_relaxed);

    for (std::size_t c = 0; c < chunks; ++c) {
        Job job;
        const F* pointer = &body;
        std::memcpy(job.payload, &pointer, sizeof(pointer));
        job.entry = [](void* payload, std::size_t b, std::size_t e) {
            const F* fn = nullptr;
            std::memcpy(&fn, payload, sizeof(fn));
            (*fn)(b, e);
        };
        job.counter = &counter;
        job.begin = begin + c * grain;
        job.end = std::min(job.begin + grain, end);
        submit(lane, job);
    }
}

template <typename F>
void JobSystem::parallelFor(std::size_t begin, std::size_t end, std::size_t grain,
                            const F& body) {
    Counter counter;
    dispatch(counter, begin, end, grain, body);
    wait(counter);
}

template <typename F>
JobSystem::AutoGrain JobSystem::parallelForAuto(std::size_t begin, std::size_t end,
                                                const F& body) {
    AutoGrain result;
    if (begin >= end) return result;
    const std::size_t total = end - begin;

    // Probe by running a growing prefix on this thread and timing it. Growing
    // geometrically matters: a fixed probe size is either too small to measure
    // for cheap elements, or a large serial stall for expensive ones. This way
    // the probe costs a bounded few microseconds either way -- and it is real
    // work, not a throwaway sample.
    using Clock = std::chrono::steady_clock;
    std::size_t probed = 0;
    double elapsedNanos = 0.0;
    std::size_t step = 1;
    while (probed < total && elapsedNanos < kMinProbeNanos) {
        const std::size_t take = std::min(step, total - probed);
        const auto start = Clock::now();
        body(begin + probed, begin + probed + take);
        elapsedNanos +=
            std::chrono::duration<double, std::nano>(Clock::now() - start).count();
        probed += take;
        step *= 4;
    }

    result.probeElements = probed;
    result.nsPerElement = probed > 0 ? elapsedNanos / static_cast<double>(probed) : 0.0;

    const std::size_t remaining = total - probed;
    if (remaining == 0) {
        result.ranSerially = true;
        return result;
    }

    const double remainingNanos = result.nsPerElement * static_cast<double>(remaining);
    if (remainingNanos < kSerialThresholdNanos) {
        body(begin + probed, end);
        result.ranSerially = true;
        return result;
    }

    // Clamp so the chunk count stays inside the band that keeps both load
    // balance and dispatch cost reasonable.
    //
    // The bounds are computed in integers, and the lower one rounds *up*. Doing
    // this in doubles and truncating at the end lets the chosen grain fall just
    // below the lower bound, which puts the chunk count one over the ceiling.
    // The normal build never reached the clamp; a ThreadSanitizer build did,
    // because its instrumentation makes the body slow enough to drive the
    // tuner into that regime.
    const std::size_t lanes = laneCount();
    const std::size_t maxChunks = std::max<std::size_t>(1, lanes * kMaxChunksPerLane);
    const std::size_t minChunks = std::max<std::size_t>(1, lanes * kMinChunksPerLane);
    const std::size_t lowerGrain = std::max<std::size_t>(1, (remaining + maxChunks - 1) / maxChunks);
    const std::size_t upperGrain =
        std::max(lowerGrain, std::max<std::size_t>(1, remaining / minChunks));

    const double targetGrain = kTargetChunkNanos / std::max(result.nsPerElement, 1e-9);
    const auto fromCost = targetGrain >= static_cast<double>(upperGrain)
                              ? upperGrain
                              : static_cast<std::size_t>(targetGrain);
    result.grain = std::clamp(fromCost, lowerGrain, upperGrain);
    result.chunks = (remaining + result.grain - 1) / result.grain;
    parallelFor(begin + probed, end, result.grain, body);
    return result;
}

template <typename T, typename Body, typename Combine>
T JobSystem::parallelReduce(std::size_t begin, std::size_t end, std::size_t grain,
                            const T& identity, const Body& body, const Combine& combine) {
    if (begin >= end) return identity;
    if (grain == 0) grain = 1;

    const std::size_t chunks = (end - begin + grain - 1) / grain;
    std::vector<T> partials(chunks, identity);

    // One bin per chunk, indexed by chunk number rather than by whichever lane
    // happened to run it. Scheduling cannot influence the fold order.
    auto worker = [&](std::size_t chunkBegin, std::size_t chunkEnd) {
        partials[(chunkBegin - begin) / grain] = body(chunkBegin, chunkEnd);
    };
    parallelFor(begin, end, grain, worker);

    T result = identity;
    for (std::size_t c = 0; c < chunks; ++c) result = combine(result, partials[c]);
    return result;
}

}  // namespace core
