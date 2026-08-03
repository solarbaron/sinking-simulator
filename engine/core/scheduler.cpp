// SPDX-License-Identifier: MIT
#include "scheduler.hpp"

#include "jobs.hpp"

#include <algorithm>
#include <cmath>

namespace core {

Scheduler::SystemId Scheduler::add(const SystemConfig& config, UpdateFn update) {
    const auto id = static_cast<SystemId>(systems_.size());
    System system;
    system.config = config;
    if (system.config.rateHz <= 0.0) system.config.rateHz = 1.0;
    if (system.config.maxCatchUpSteps == 0) system.config.maxCatchUpSteps = 1;
    // Round the period once, here, and never touch the rate in floating point
    // again. Every subsequent step-count decision is integer division.
    system.periodNanos = std::max<std::int64_t>(
        1, static_cast<std::int64_t>(std::llround(1.0e9 / system.config.rateHz)));
    system.update = std::move(update);
    systems_.push_back(std::move(system));
    levelsDirty_ = true;
    return id;
}

bool Scheduler::wouldCycle(SystemId system, SystemId prerequisite) const {
    // A cycle appears exactly when `system` is already reachable from
    // `prerequisite` by following prerequisite edges.
    std::vector<SystemId> stack{prerequisite};
    std::vector<bool> seen(systems_.size(), false);
    while (!stack.empty()) {
        const SystemId current = stack.back();
        stack.pop_back();
        if (current == system) return true;
        if (seen[current]) continue;
        seen[current] = true;
        for (SystemId next : systems_[current].prerequisites) stack.push_back(next);
    }
    return false;
}

bool Scheduler::dependsOn(SystemId system, SystemId prerequisite) {
    if (system >= systems_.size() || prerequisite >= systems_.size()) return false;
    if (system == prerequisite) return false;
    if (wouldCycle(system, prerequisite)) return false;

    auto& list = systems_[system].prerequisites;
    if (std::find(list.begin(), list.end(), prerequisite) == list.end())
        list.push_back(prerequisite);
    levelsDirty_ = true;
    return true;
}

void Scheduler::setTimeDilation(double simulationSecondsPerWallSecond) {
    dilation_ = std::max(0.0, simulationSecondsPerWallSecond);
}

bool Scheduler::withinDilationBand(SystemId id) const {
    const SystemConfig& config = systems_[id].config;
    return dilation_ >= config.minDilation && dilation_ <= config.maxDilation;
}

void Scheduler::rebuildLevels() const {
    // Kahn's algorithm, but keeping the levels rather than flattening them: a
    // level is a set of systems with no dependency path between any two, which
    // is exactly the set that may run in parallel.
    const std::size_t count = systems_.size();
    std::vector<std::uint32_t> remaining(count, 0);
    std::vector<std::vector<SystemId>> dependents(count);

    for (std::size_t s = 0; s < count; ++s)
        for (SystemId prerequisite : systems_[s].prerequisites) {
            ++remaining[s];
            dependents[prerequisite].push_back(static_cast<SystemId>(s));
        }

    levels_.clear();
    std::vector<SystemId> current;
    for (std::size_t s = 0; s < count; ++s)
        if (remaining[s] == 0) current.push_back(static_cast<SystemId>(s));

    std::size_t placed = 0;
    while (!current.empty()) {
        // Ascending id within a level, so the grouping is reproducible rather
        // than dependent on the order edges happened to be added.
        std::sort(current.begin(), current.end());
        placed += current.size();
        std::vector<SystemId> next;
        for (SystemId id : current)
            for (SystemId dependent : dependents[id])
                if (--remaining[dependent] == 0) next.push_back(dependent);
        levels_.push_back(std::move(current));
        current = std::move(next);
    }

    // dependsOn() refuses to create cycles, so this cannot trigger; if it ever
    // does, running the stragglers last beats dropping them silently.
    if (placed < count) {
        std::vector<SystemId> stragglers;
        for (std::size_t s = 0; s < count; ++s)
            if (remaining[s] > 0) stragglers.push_back(static_cast<SystemId>(s));
        levels_.push_back(std::move(stragglers));
    }
    levelsDirty_ = false;
}

const std::vector<std::vector<Scheduler::SystemId>>& Scheduler::levels() const {
    if (levelsDirty_) rebuildLevels();
    return levels_;
}

Scheduler::Report Scheduler::advance(double wallSeconds) {
    Report report;
    if (wallSeconds <= 0.0 || systems_.empty()) return report;

    // The one floating-point step: convert this frame's wall delta into whole
    // nanoseconds of simulation time. Everything downstream is exact.
    const auto simulationDelta =
        static_cast<std::int64_t>(std::llround(wallSeconds * dilation_ * 1.0e9));
    if (simulationDelta <= 0) return report;
    simulationNanos_ += simulationDelta;

    // Decide every system's step count up front, from simulation time alone.
    // Doing this before any system runs is what keeps the schedule a pure
    // function of the inputs: a system that took longer than expected cannot
    // change how many times its neighbours run.
    std::vector<std::uint32_t> steps(systems_.size(), 0);
    for (std::size_t s = 0; s < systems_.size(); ++s) {
        System& system = systems_[s];
        if (!withinDilationBand(static_cast<SystemId>(s))) {
            ++report.systemsSkipped;
            continue;
        }
        system.accumulatorNanos += simulationDelta;

        auto wanted = static_cast<std::uint64_t>(system.accumulatorNanos / system.periodNanos);
        // Exact: the remainder carries forward with no drift, however many
        // advances have gone before.
        system.accumulatorNanos %= system.periodNanos;
        if (wanted > system.config.maxCatchUpSteps) {
            report.droppedCatchUp += wanted - system.config.maxCatchUpSteps;
            // The surplus is discarded, not carried: keeping it would only
            // guarantee the same overrun next time.
            wanted = system.config.maxCatchUpSteps;
        }
        steps[s] = static_cast<std::uint32_t>(wanted);
        report.updates += wanted;
    }

    for (const auto& level : levels()) {
        auto runSystem = [&](SystemId id) {
            System& system = systems_[id];
            const double dt = static_cast<double>(system.periodNanos) * 1e-9;
            for (std::uint32_t i = 0; i < steps[id]; ++i) {
                if (system.update) system.update(dt);
                ++system.updates;
            }
        };

        // Everything in a level is independent by declaration, so it is safe to
        // spread across the pool. A system's own repeats stay ordered.
        if (jobs_ != nullptr && level.size() > 1) {
            jobs_->parallelFor(0, level.size(), 1, [&](std::size_t begin, std::size_t end) {
                for (std::size_t i = begin; i < end; ++i) runSystem(level[i]);
            });
        } else {
            for (SystemId id : level) runSystem(id);
        }
    }
    return report;
}

}  // namespace core
