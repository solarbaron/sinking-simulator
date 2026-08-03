// SPDX-License-Identifier: MIT
#include "jobs.hpp"

#include <algorithm>
#include <condition_variable>
#include <mutex>
#include <thread>

namespace core {
namespace {

// How many spare lanes exist for threads that are not workers. One is enough for
// the normal engine model (a single main thread submits); the rest are slack for
// tools and tests.
constexpr unsigned kExternalLanes = 4;

// Failed acquire attempts before a worker stops spinning and sleeps. Spinning
// briefly is worth it because job durations are short and a sleep/wake round trip
// costs more than the job itself.
constexpr int kSpinsBeforeSleep = 64;

struct LaneBinding {
    const JobSystem* owner = nullptr;
    unsigned lane = 0;
};
thread_local LaneBinding t_binding;

std::uint64_t nextRandom(std::uint64_t& state) {
    // xorshift64*, enough for choosing a steal victim.
    state ^= state >> 12;
    state ^= state << 25;
    state ^= state >> 27;
    return state * 0x2545F4914F6CDD1Dull;
}

}  // namespace

struct JobSystem::Impl {
    std::vector<std::thread> workers;
    std::atomic<bool> stopping{false};
    std::atomic<unsigned> nextExternalLane{0};

    std::mutex sleepMutex;
    std::condition_variable sleepSignal;
    std::atomic<int> sleepers{0};
};

// --- bounded MPMC ring ------------------------------------------------------

JobSystem::Queue::Queue() {
    // Cell i starts holding ticket i, so the first producer to claim position i
    // finds a match and may write it.
    for (std::size_t i = 0; i < kCapacity; ++i)
        cells_[i].sequence.store(i, std::memory_order_relaxed);
}

bool JobSystem::Queue::push(const Job& job) {
    std::size_t pos = tail_.load(std::memory_order_relaxed);
    for (;;) {
        Cell& cell = cells_[pos & kMask];
        const std::size_t sequence = cell.sequence.load(std::memory_order_acquire);
        const auto delta = static_cast<std::intptr_t>(sequence) - static_cast<std::intptr_t>(pos);
        if (delta == 0) {
            if (tail_.compare_exchange_weak(pos, pos + 1, std::memory_order_relaxed)) {
                cell.job = job;
                // Publishes the write above to whichever consumer claims this cell.
                cell.sequence.store(pos + 1, std::memory_order_release);
                return true;
            }
        } else if (delta < 0) {
            return false;  // the ring is full
        } else {
            pos = tail_.load(std::memory_order_relaxed);
        }
    }
}

bool JobSystem::Queue::pop(Job& out) {
    std::size_t pos = head_.load(std::memory_order_relaxed);
    for (;;) {
        Cell& cell = cells_[pos & kMask];
        const std::size_t sequence = cell.sequence.load(std::memory_order_acquire);
        const auto delta =
            static_cast<std::intptr_t>(sequence) - static_cast<std::intptr_t>(pos + 1);
        if (delta == 0) {
            if (head_.compare_exchange_weak(pos, pos + 1, std::memory_order_relaxed)) {
                out = cell.job;
                // Hands the cell to the producer one lap ahead.
                cell.sequence.store(pos + kCapacity, std::memory_order_release);
                return true;
            }
        } else if (delta < 0) {
            return false;  // empty
        } else {
            pos = head_.load(std::memory_order_relaxed);
        }
    }
}

// --- lifecycle ---------------------------------------------------------------

unsigned JobSystem::defaultWorkerCount() {
    const unsigned hardware = std::thread::hardware_concurrency();
    return hardware > 1 ? hardware - 1 : 0;
}

JobSystem::JobSystem(unsigned workerCount) : workerCount_(workerCount) {
    impl_ = new Impl();
    lanes_ = std::vector<Lane>(workerCount_ + kExternalLanes);
    for (std::size_t i = 0; i < lanes_.size(); ++i)
        lanes_[i].rng = 0x9E3779B97F4A7C15ull + i * 0x0F1E2D3C4B5A6978ull;

    impl_->workers.reserve(workerCount_);
    for (unsigned w = 0; w < workerCount_; ++w)
        impl_->workers.emplace_back([this, w] { workerMain(w); });
}

JobSystem::~JobSystem() {
    impl_->stopping.store(true, std::memory_order_release);
    {
        std::lock_guard<std::mutex> lock(impl_->sleepMutex);
        impl_->sleepSignal.notify_all();
    }
    for (std::thread& worker : impl_->workers)
        if (worker.joinable()) worker.join();

    // Release any external lane this thread claimed, so a later JobSystem in the
    // same process starts from a clean binding.
    if (t_binding.owner == this) t_binding = {};
    delete impl_;
}

unsigned JobSystem::currentLane() const {
    if (t_binding.owner == this) return t_binding.lane;
    // Not a worker and not yet bound: claim one of the external lanes.
    const unsigned slot = impl_->nextExternalLane.fetch_add(1, std::memory_order_relaxed);
    t_binding.owner = this;
    t_binding.lane = workerCount_ + (slot % kExternalLanes);
    return t_binding.lane;
}

// --- submission and execution ------------------------------------------------

void JobSystem::submit(unsigned lane, const Job& job) {
    if (!lanes_[lane].queue.push(job)) {
        // Queue full. Running it inline is the correct back-pressure valve: it
        // cannot deadlock, and it degrades to sequential execution rather than
        // dropping work or growing without bound.
        Job local = job;
        execute(local);
        return;
    }
    if (impl_->sleepers.load(std::memory_order_acquire) > 0) wakeWorkers(1);
}

void JobSystem::execute(Job& job) {
    job.entry(job.payload, job.begin, job.end);
    job.counter->outstanding_.fetch_sub(1, std::memory_order_release);
}

bool JobSystem::acquireJob(unsigned lane, Job& out) {
    if (lanes_[lane].queue.pop(out)) return true;

    const unsigned count = static_cast<unsigned>(lanes_.size());
    if (count < 2) return false;

    // Randomised victim selection, then a full sweep. Random start avoids the
    // convoy effect where every idle lane hammers lane 0.
    const unsigned start = static_cast<unsigned>(nextRandom(lanes_[lane].rng) % count);
    for (unsigned i = 0; i < count; ++i) {
        const unsigned victim = (start + i) % count;
        if (victim == lane) continue;
        if (lanes_[victim].queue.pop(out)) return true;
    }
    return false;
}

void JobSystem::wakeWorkers(int count) {
    std::lock_guard<std::mutex> lock(impl_->sleepMutex);
    if (count == 1) impl_->sleepSignal.notify_one();
    else impl_->sleepSignal.notify_all();
}

void JobSystem::wait(Counter& counter) {
    const unsigned lane = currentLane();
    Job job;
    int idle = 0;
    while (!counter.complete()) {
        if (acquireJob(lane, job)) {
            execute(job);
            idle = 0;
            continue;
        }
        // Nothing to help with. Never sleep here: the counter may be finished by
        // a job already running on another lane, and there would be no submission
        // left to wake us.
        if (++idle > kSpinsBeforeSleep) {
            std::this_thread::yield();
            idle = 0;
        }
    }
}

void JobSystem::workerMain(unsigned lane) {
    t_binding.owner = this;
    t_binding.lane = lane;

    Job job;
    int idle = 0;
    while (!impl_->stopping.load(std::memory_order_acquire)) {
        if (acquireJob(lane, job)) {
            execute(job);
            idle = 0;
            continue;
        }
        if (++idle <= kSpinsBeforeSleep) {
            std::this_thread::yield();
            continue;
        }
        // Sleep with a timeout rather than relying purely on the wake signal.
        // A submitter only wakes workers when it sees sleepers > 0, and that
        // check races with a worker deciding to sleep; the timeout closes it.
        impl_->sleepers.fetch_add(1, std::memory_order_release);
        {
            std::unique_lock<std::mutex> lock(impl_->sleepMutex);
            impl_->sleepSignal.wait_for(lock, std::chrono::microseconds(200));
        }
        impl_->sleepers.fetch_sub(1, std::memory_order_release);
        idle = 0;
    }
}

}  // namespace core
