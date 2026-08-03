// SPDX-License-Identifier: MIT
//
// Validation of the bump allocators.
//
// A bump allocator is a handful of lines and almost all of its failure modes are
// silent: misalignment that only bites on over-aligned types, an off-by-one that
// hands out overlapping blocks, an exhaustion path that wraps instead of failing.
// Every check here is an invariant, and the per-lane test runs real jobs on real
// threads so the "one arena per lane, no atomics" claim is exercised rather than
// asserted.
#include "engine/core/arena.hpp"
#include "engine/core/jobs.hpp"
#include "harness.hpp"

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

using core::Arena;
using core::FrameArenas;
using core::JobSystem;
using testing::expectEqual;
using testing::expectTrue;

namespace {

bool isAligned(const void* pointer, std::size_t alignment) {
    return (reinterpret_cast<std::uintptr_t>(pointer) & (alignment - 1)) == 0;
}

// Over-aligned types are the case that silently works until it does not: a
// naively written bump allocator returns 8-byte-aligned memory and nothing
// complains until an aligned SIMD store faults.
void testAlignmentIsHonoured() {
    Arena arena(1 << 20);
    for (std::size_t alignment : {1u, 2u, 4u, 8u, 16u, 32u, 64u, 128u, 256u}) {
        for (int repeat = 0; repeat < 8; ++repeat) {
            // An odd size between requests so the next offset is never already
            // aligned by luck.
            void* filler = arena.allocate(3, 1);
            void* pointer = arena.allocate(17, alignment);
            expectTrue("allocation honours alignment " + std::to_string(alignment),
                       pointer != nullptr && isAligned(pointer, alignment));
            expectTrue("filler allocation succeeded", filler != nullptr);
        }
    }
}

// Overlapping allocations are the failure that functional code hides best: the
// second write simply corrupts the first caller's data.
void testAllocationsDoNotOverlap() {
    Arena arena(1 << 20);
    struct Block {
        std::byte* pointer;
        std::size_t size;
        std::byte pattern;
    };
    std::vector<Block> blocks;

    for (int i = 0; i < 500; ++i) {
        const std::size_t size = 1 + static_cast<std::size_t>((i * 37) % 200);
        const std::size_t alignment = std::size_t{1} << (i % 7);
        auto* pointer = static_cast<std::byte*>(arena.allocate(size, alignment));
        if (pointer == nullptr) break;
        const auto pattern = static_cast<std::byte>(i & 0xFF);
        std::memset(pointer, static_cast<int>(pattern), size);
        blocks.push_back({pointer, size, pattern});
    }
    expectTrue("the arena served a useful number of blocks", blocks.size() > 400);

    // Every block must still read back as its own pattern.
    int corrupted = 0;
    for (const Block& block : blocks)
        for (std::size_t i = 0; i < block.size; ++i)
            if (block.pointer[i] != block.pattern) {
                ++corrupted;
                break;
            }
    expectEqual("no block was overwritten by another", corrupted, 0);

    // And no two blocks may share a byte.
    std::vector<Block> sorted = blocks;
    std::sort(sorted.begin(), sorted.end(),
              [](const Block& a, const Block& b) { return a.pointer < b.pointer; });
    int overlaps = 0;
    for (std::size_t i = 1; i < sorted.size(); ++i)
        if (sorted[i - 1].pointer + sorted[i - 1].size > sorted[i].pointer) ++overlaps;
    expectEqual("no two allocations overlap", overlaps, 0);
}

// Exhaustion must fail cleanly. The dangerous version wraps the size computation
// and returns a pointer near the end of the block with far less room than asked.
void testExhaustionFailsCleanly() {
    Arena arena(1024);
    void* big = arena.allocate(900, 8);
    expectTrue("first large allocation fits", big != nullptr);

    const std::size_t usedBefore = arena.used();
    void* tooBig = arena.allocate(500, 8);
    expectTrue("allocation past capacity returns nullptr", tooBig == nullptr);
    expectEqual("a failed allocation does not consume space",
                static_cast<long long>(arena.used()), static_cast<long long>(usedBefore));

    // The arena must still work afterwards.
    void* small = arena.allocate(16, 8);
    expectTrue("arena still usable after a failed allocation", small != nullptr);

    // Size arithmetic must not wrap.
    void* absurd = arena.allocate(SIZE_MAX - 8, 8);
    expectTrue("an absurd size fails rather than wrapping", absurd == nullptr);
    auto huge = arena.allocateArray<double>(SIZE_MAX / 4);
    expectTrue("an overflowing array count fails rather than wrapping", huge.empty());

    Arena empty(0);
    expectTrue("a zero-capacity arena allocates nothing", empty.allocate(1, 1) == nullptr);
}

void testResetAndHighWaterMark() {
    Arena arena(4096);
    expectEqual("a fresh arena is empty", static_cast<long long>(arena.used()), 0);

    void* first = arena.allocate(1000, 8);
    expectTrue("first allocation succeeded", first != nullptr);
    const std::size_t peak = arena.used();

    arena.reset();
    expectEqual("reset returns the arena to empty", static_cast<long long>(arena.used()), 0);
    expectEqual("reset preserves the high-water mark",
                static_cast<long long>(arena.highWaterMark()), static_cast<long long>(peak));

    // Memory must be genuinely reusable, and hand back the same region.
    void* second = arena.allocate(1000, 8);
    expectTrue("reset makes the same memory available again", second == first);

    // A larger peak must be recorded; a smaller one must not lower it.
    arena.allocate(2000, 8);
    const std::size_t higher = arena.used();
    expectTrue("high-water mark rises with a bigger peak", arena.highWaterMark() == higher);
    arena.reset();
    arena.allocate(8, 8);
    expectEqual("high-water mark does not fall back",
                static_cast<long long>(arena.highWaterMark()), static_cast<long long>(higher));
}

void testScopeRewinds() {
    Arena arena(4096);
    arena.allocate(64, 8);
    const std::size_t outer = arena.used();

    void* inner = nullptr;
    {
        Arena::Scope scope(arena);
        inner = arena.allocate(128, 8);
        expectTrue("scoped allocation succeeded", inner != nullptr);
        expectTrue("scoped allocation advances the arena", arena.used() > outer);

        {
            Arena::Scope nested(arena);
            arena.allocate(256, 8);
        }
        expectTrue("inner scope rewound", arena.used() < outer + 256);
    }
    expectEqual("outer scope restored the mark", static_cast<long long>(arena.used()),
                static_cast<long long>(outer));

    // After rewinding, the same address is handed out again.
    void* reused = arena.allocate(128, 8);
    expectTrue("rewound memory is reissued", reused == inner);
}

void testCreateAndArray() {
    Arena arena(1 << 16);

    struct Particle {
        double x, y, z;
        int id;
    };
    Particle* particle = arena.create<Particle>(1.0, 2.0, 3.0, 42);
    expectTrue("create<T> returned an object", particle != nullptr);
    expectTrue("create<T> ran the constructor arguments",
               particle->x == 1.0 && particle->y == 2.0 && particle->z == 3.0 &&
                   particle->id == 42);
    expectTrue("create<T> honours the type's alignment", isAligned(particle, alignof(Particle)));

    auto values = arena.allocateArray<std::uint64_t>(1000);
    expectEqual("allocateArray returned the requested count",
                static_cast<long long>(values.size()), 1000);
    expectTrue("array is correctly aligned", isAligned(values.data(), alignof(std::uint64_t)));
    // Writable across its whole extent -- under ASan this also proves the region
    // was unpoisoned exactly, since a short unpoison would fault here.
    for (std::size_t i = 0; i < values.size(); ++i) values[i] = i * 2654435761u;
    bool intact = true;
    for (std::size_t i = 0; i < values.size(); ++i)
        if (values[i] != i * 2654435761u) intact = false;
    expectTrue("array contents survive", intact);

    expectTrue("a zero-length array request is empty, not null-dereferencing",
               arena.allocateArray<int>(0).empty());
}

void testMoveTransfersOwnership() {
    Arena source(4096);
    void* pointer = source.allocate(100, 8);
    const std::size_t used = source.used();

    Arena target(std::move(source));
    expectEqual("moved-to arena keeps the offset", static_cast<long long>(target.used()),
                static_cast<long long>(used));
    expectEqual("moved-from arena has no capacity", static_cast<long long>(source.capacity()), 0);
    expectTrue("moved-from arena allocates nothing", source.allocate(1, 1) == nullptr);

    // The original pointer must still belong to the moved-to arena.
    void* next = target.allocate(8, 8);
    expectTrue("moved-to arena keeps serving the same block",
               next != nullptr && next > pointer);
}

// The load-bearing claim about FrameArenas is that one arena per lane means
// allocation needs no synchronisation. Run real jobs on real threads and check
// that nothing any lane handed out overlaps anything another lane handed out.
void testPerLaneArenasAreIsolated() {
    JobSystem jobs(8);
    FrameArenas arenas(jobs.laneCount(), 1 << 20);

    struct Block {
        std::byte* pointer;
        std::size_t size;
        unsigned lane;
    };
    // One record vector per lane, touched only by that lane's thread.
    std::vector<std::vector<Block>> perLane(jobs.laneCount());
    for (auto& v : perLane) v.reserve(4096);

    jobs.parallelFor(0, 20000, 16, [&](std::size_t begin, std::size_t end) {
        const unsigned lane = jobs.currentLane();
        Arena& arena = arenas.lane(lane);
        for (std::size_t i = begin; i < end; ++i) {
            const std::size_t size = 8 + (i % 64);
            auto* pointer = static_cast<std::byte*>(arena.allocate(size, 8));
            if (pointer == nullptr) continue;  // lane arena full; still valid
            std::memset(pointer, static_cast<int>(lane & 0xFF), size);
            perLane[lane].push_back({pointer, size, lane});
        }
    });

    std::vector<Block> all;
    for (const auto& lane : perLane) all.insert(all.end(), lane.begin(), lane.end());
    expectTrue("lanes allocated a useful number of blocks", all.size() > 10000);

    // Each block must still hold its own lane's byte: proof that no other lane
    // was handed the same memory.
    int corrupted = 0;
    for (const Block& block : all) {
        const auto expected = static_cast<std::byte>(block.lane & 0xFF);
        for (std::size_t i = 0; i < block.size; ++i)
            if (block.pointer[i] != expected) {
                ++corrupted;
                break;
            }
    }
    expectEqual("no lane overwrote another lane's memory", corrupted, 0);

    std::sort(all.begin(), all.end(),
              [](const Block& a, const Block& b) { return a.pointer < b.pointer; });
    int overlaps = 0;
    for (std::size_t i = 1; i < all.size(); ++i)
        if (all[i - 1].pointer + all[i - 1].size > all[i].pointer) ++overlaps;
    expectEqual("no allocation from any lane overlaps any other", overlaps, 0);

    const std::size_t peak = arenas.totalHighWaterMark();
    expectTrue("frame arenas report a plausible high-water mark",
               peak > 0 && peak <= arenas.totalCapacity());

    arenas.resetAll();
    std::size_t usedAfterReset = 0;
    for (unsigned lane = 0; lane < arenas.laneCount(); ++lane)
        usedAfterReset += arenas.lane(lane).used();
    expectEqual("resetAll empties every lane", static_cast<long long>(usedAfterReset), 0);
    expectEqual("resetAll preserves the high-water mark",
                static_cast<long long>(arenas.totalHighWaterMark()),
                static_cast<long long>(peak));
}

}  // namespace

void runArenaTests() {
    std::printf("\n--- arena allocators ---\n");
    testAlignmentIsHonoured();
    testAllocationsDoNotOverlap();
    testExhaustionFailsCleanly();
    testResetAndHighWaterMark();
    testScopeRewinds();
    testCreateAndArray();
    testMoveTransfersOwnership();
    testPerLaneArenasAreIsolated();
}
