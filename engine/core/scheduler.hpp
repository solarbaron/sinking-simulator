// SPDX-License-Identifier: MIT
//
// Multi-rate simulation scheduler with time dilation.
//
// Nothing in this engine runs at one rate. Hydrostatics is happy at 100 Hz,
// thermal conduction at 10 Hz, an active FEM damage zone wants kilohertz, and a
// flooding fight worth watching lasts forty minutes while a hull rupture worth
// watching lasts fifty milliseconds. So the scheduler owns two things: a
// per-system rate in *simulation* time, and the ratio between simulation and
// wall time.
//
// Time dilation is deliberately a first-class input rather than a debug toggle.
// It is what makes the expensive solvers affordable -- a zone that cannot run at
// 1x can still run at 1/50x -- and it is simultaneously the player-facing feature
// of slowing down for a rupture and fast-forwarding through a long flood. A
// system declares the dilation band it can honour and is simply skipped outside
// it, which is how a fast-forward drops the FEM to its reduced model.
//
// Simulation time is accumulated in **integer nanoseconds**, not seconds in a
// double. Floating point is used exactly once per advance, converting the wall
// delta; everything after that is exact integer arithmetic. This matters because
// the accumulation is unbounded -- a flooding casualty runs for forty minutes of
// simulated time and a campaign far longer -- and a double accumulator both
// drifts and makes step counts depend on the order the slices arrived in. With
// integers, a given sequence of advance() calls produces exactly the same
// schedule, on any platform, however long it runs.
#pragma once

#include <cstdint>
#include <functional>
#include <limits>
#include <string>
#include <vector>

namespace core {

class JobSystem;

struct SystemConfig {
    std::string name;
    // Updates per second of *simulation* time.
    double rateHz = 100.0;
    // Dilation band over which this system runs at all. Outside it the system is
    // skipped entirely and its accumulator is held, which is how a 60x
    // fast-forward drops solvers that cannot honour their rate.
    double minDilation = 0.0;
    double maxDilation = std::numeric_limits<double>::infinity();
    // Ceiling on catch-up updates in a single advance(). Without this, one long
    // frame asks for more updates than the next frame has time to run, which
    // asks for more still -- the spiral of death.
    std::uint32_t maxCatchUpSteps = 8;
};

class Scheduler {
public:
    using SystemId = std::uint32_t;
    // Invoked with the system's own fixed simulation timestep, never a variable
    // one: a solver tuned for 1/rate seconds should never be handed anything else.
    using UpdateFn = std::function<void(double simulationDt)>;

    // `jobs` may be null, in which case everything runs on the calling thread.
    // Systems with no dependency path between them are run in parallel when it
    // is not.
    explicit Scheduler(JobSystem* jobs = nullptr) : jobs_(jobs) {}

    SystemId add(const SystemConfig& config, UpdateFn update);
    // Declares that `system` must not run until `prerequisite` has. Returns
    // false if the edge would close a cycle, leaving the graph untouched.
    bool dependsOn(SystemId system, SystemId prerequisite);

    void setTimeDilation(double simulationSecondsPerWallSecond);
    double timeDilation() const { return dilation_; }

    struct Report {
        std::uint64_t updates = 0;        // updates actually run
        std::uint64_t droppedCatchUp = 0; // updates skipped by the catch-up ceiling
        std::uint32_t systemsSkipped = 0; // systems outside their dilation band
    };
    Report advance(double wallSeconds);

    std::int64_t simulationNanos() const { return simulationNanos_; }
    double simulationTime() const { return static_cast<double>(simulationNanos_) * 1e-9; }
    // A system's period in nanoseconds, as the scheduler rounded it. Exposed so
    // callers and tests can predict step counts exactly rather than in floating
    // point.
    std::int64_t periodNanos(SystemId id) const { return systems_[id].periodNanos; }
    std::uint64_t updateCount(SystemId id) const { return systems_[id].updates; }
    bool withinDilationBand(SystemId id) const;
    std::size_t systemCount() const { return systems_.size(); }
    // Systems grouped into dependency levels; everything in one level is
    // independent of everything else in it. Exposed for tests and for the
    // profiler's frame view.
    const std::vector<std::vector<SystemId>>& levels() const;

private:
    struct System {
        SystemConfig config;
        UpdateFn update;
        std::int64_t periodNanos = 0;
        std::int64_t accumulatorNanos = 0;
        std::uint64_t updates = 0;
        std::vector<SystemId> prerequisites;
    };

    void rebuildLevels() const;
    bool wouldCycle(SystemId system, SystemId prerequisite) const;

    JobSystem* jobs_ = nullptr;
    std::vector<System> systems_;
    double dilation_ = 1.0;
    std::int64_t simulationNanos_ = 0;

    mutable std::vector<std::vector<SystemId>> levels_;
    mutable bool levelsDirty_ = true;
};

}  // namespace core
