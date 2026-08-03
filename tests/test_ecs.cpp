// SPDX-License-Identifier: MIT
//
// Validation of the archetype ECS.
//
// The dangerous failures here are all silent. A swap-remove that forgets to
// patch the moved entity's record does not crash -- two handles simply start
// aliasing one row. An archetype transition that drops a component leaves stale
// bytes that read as plausible values. A chunk capacity computed loosely
// overruns into the next array. Every check below targets one of those.
#include "engine/core/ecs.hpp"
#include "harness.hpp"

#include <algorithm>
#include <cstdio>
#include <set>
#include <string>
#include <vector>

using core::Entity;
using core::World;
using testing::expectEqual;
using testing::expectTrue;

namespace {

struct Position {
    double x = 0, y = 0, z = 0;
};
struct Velocity {
    double x = 0, y = 0, z = 0;
};
struct Mass {
    double kg = 0;
};
struct Tag {
    std::uint32_t value = 0;
};
// Over-aligned, to prove chunk offsets respect a component's own alignment.
struct alignas(64) Wide {
    double values[8]{};
};

void testEntityLifetime() {
    World world;
    expectEqual("a new world is empty", static_cast<long long>(world.entityCount()), 0);

    const Entity a = world.create();
    const Entity b = world.create();
    expectTrue("created entities are alive", world.alive(a) && world.alive(b));
    expectTrue("entities are distinct", !(a == b));
    expectEqual("entity count tracks creation", static_cast<long long>(world.entityCount()), 2);

    world.destroy(a);
    expectTrue("a destroyed entity is not alive", !world.alive(a));
    expectTrue("its neighbour is unaffected", world.alive(b));
    expectEqual("entity count tracks destruction", static_cast<long long>(world.entityCount()), 1);

    // Index reuse must not resurrect the old handle. Without a generation
    // counter the stale handle would silently address the new entity.
    const Entity c = world.create();
    expectEqual("a freed index is reused", static_cast<long long>(c.index),
                static_cast<long long>(a.index));
    expectTrue("the reused handle has a new generation", c.generation != a.generation);
    expectTrue("the stale handle is still dead", !world.alive(a));
    expectTrue("the new handle is alive", world.alive(c));

    world.destroy(a);  // must be a no-op, not a double free
    expectTrue("destroying a stale handle does not affect the live one", world.alive(c));
}

void testComponentAccess() {
    World world;
    const Entity e = world.create(Position{1, 2, 3}, Mass{500.0});

    expectTrue("created components are present", world.has<Position>(e) && world.has<Mass>(e));
    expectTrue("absent components report absent", !world.has<Velocity>(e));
    expectTrue("absent components return nullptr", world.get<Velocity>(e) == nullptr);

    const Position* position = world.get<Position>(e);
    expectTrue("component values survive creation",
               position != nullptr && position->x == 1 && position->y == 2 && position->z == 3);
    expectTrue("second component survives too", world.get<Mass>(e)->kg == 500.0);

    world.get<Position>(e)->y = 99.0;
    expectTrue("components are mutable in place", world.get<Position>(e)->y == 99.0);

    expectTrue("a dead entity yields no components",
               (world.destroy(e), world.get<Position>(e)) == nullptr);
}

// Adding or removing a component relocates the row to another archetype. The
// components that survive the transition must arrive intact.
void testArchetypeTransitionsPreserveData() {
    World world;
    const Entity e = world.create(Position{7, 8, 9}, Mass{1234.0});

    world.add(e, Velocity{0.5, 1.5, 2.5});
    expectTrue("added component is readable",
               world.has<Velocity>(e) && world.get<Velocity>(e)->y == 1.5);
    expectTrue("position survived the move to a wider archetype",
               world.get<Position>(e)->x == 7 && world.get<Position>(e)->z == 9);
    expectTrue("mass survived the move", world.get<Mass>(e)->kg == 1234.0);

    world.remove<Mass>(e);
    expectTrue("removed component is gone", !world.has<Mass>(e));
    expectTrue("position survived the move to a narrower archetype",
               world.get<Position>(e)->x == 7);
    expectTrue("velocity survived too", world.get<Velocity>(e)->z == 2.5);

    // Re-adding must not resurrect the old bytes.
    world.add(e, Mass{42.0});
    expectTrue("re-added component takes the new value", world.get<Mass>(e)->kg == 42.0);

    // Adding something already present is a value update, not a relocation.
    const std::size_t archetypes = world.archetypeCount();
    world.add(e, Position{-1, -2, -3});
    expectEqual("re-adding an existing component creates no archetype",
                static_cast<long long>(world.archetypeCount()),
                static_cast<long long>(archetypes));
    expectTrue("re-adding an existing component overwrites it", world.get<Position>(e)->x == -1);
}

// The classic archetype ECS bug: removing a row swaps the last row into the gap,
// and if the moved entity's record is not patched, two handles alias one row.
// Destroying from the front is what exposes it.
void testSwapRemovePatchesMovedEntity() {
    World world;
    constexpr int kCount = 500;
    std::vector<Entity> entities;
    for (int i = 0; i < kCount; ++i)
        entities.push_back(world.create(Tag{static_cast<std::uint32_t>(i)}));

    // Destroy from the front, so almost every removal triggers a swap.
    for (int i = 0; i < kCount; i += 2) world.destroy(entities[static_cast<std::size_t>(i)]);

    int wrong = 0;
    for (int i = 1; i < kCount; i += 2) {
        const Entity e = entities[static_cast<std::size_t>(i)];
        const Tag* tag = world.get<Tag>(e);
        if (tag == nullptr || tag->value != static_cast<std::uint32_t>(i)) ++wrong;
    }
    expectEqual("every survivor still reads its own tag after swap-removes", wrong, 0);
    expectEqual("entity count is right after the cull",
                static_cast<long long>(world.entityCount()), kCount / 2);

    // And every survivor must occupy a distinct row: if two aliased, a write
    // through one would be visible through the other.
    for (int i = 1; i < kCount; i += 2)
        world.get<Tag>(entities[static_cast<std::size_t>(i)])->value = 0xF0000000u + static_cast<std::uint32_t>(i);
    int aliased = 0;
    for (int i = 1; i < kCount; i += 2)
        if (world.get<Tag>(entities[static_cast<std::size_t>(i)])->value !=
            0xF0000000u + static_cast<std::uint32_t>(i))
            ++aliased;
    expectEqual("no two survivors share a row", aliased, 0);
}

// Iteration must visit exactly the entities holding every requested component,
// once each -- and hand back contiguous arrays.
void testEachVisitsMatchingEntitiesOnce() {
    World world;
    std::vector<Entity> withBoth, withPositionOnly;

    for (int i = 0; i < 300; ++i)
        withBoth.push_back(world.create(Position{static_cast<double>(i), 0, 0},
                                        Velocity{1, 0, 0}, Tag{static_cast<std::uint32_t>(i)}));
    for (int i = 0; i < 200; ++i)
        withPositionOnly.push_back(world.create(Position{-1, 0, 0}));

    std::size_t visited = 0;
    world.each<Position, Velocity>([&](std::size_t count, Position* p, Velocity* v) {
        for (std::size_t i = 0; i < count; ++i) {
            p[i].x += v[i].x;  // arrays must be parallel and contiguous
            ++visited;
        }
    });
    expectEqual("each visits exactly the matching entities",
                static_cast<long long>(visited), 300);

    int wrong = 0;
    for (std::size_t i = 0; i < withBoth.size(); ++i)
        if (world.get<Position>(withBoth[i])->x != static_cast<double>(i) + 1.0) ++wrong;
    expectEqual("each wrote through to the right entities", wrong, 0);

    int touched = 0;
    for (const Entity e : withPositionOnly)
        if (world.get<Position>(e)->x != -1.0) ++touched;
    expectEqual("non-matching entities were not touched", touched, 0);

    // A component nobody has must match nothing rather than everything.
    std::size_t massVisits = 0;
    world.each<Mass>([&](std::size_t count, Mass*) { massVisits += count; });
    expectEqual("querying an unused component visits nothing",
                static_cast<long long>(massVisits), 0);
}

// Chunks are fixed size, so a large population must span several of them. A
// capacity computed loosely would overrun one array into the next.
void testMultipleChunks() {
    World world;
    const std::size_t capacity = world.chunkCapacityFor<Position, Velocity>();
    expectTrue("chunk capacity is plausible for 48 bytes of components",
               capacity > 50 && capacity < core::Chunk::kBytes);

    const std::size_t count = capacity * 3 + 7;  // forces four chunks
    std::vector<Entity> entities;
    for (std::size_t i = 0; i < count; ++i)
        entities.push_back(world.create(Position{static_cast<double>(i), 0, 0},
                                        Velocity{static_cast<double>(i) * 2, 0, 0}));

    std::size_t visited = 0;
    world.each<Position, Velocity>([&](std::size_t n, Position*, Velocity*) { visited += n; });
    expectEqual("iteration spans every chunk", static_cast<long long>(visited),
                static_cast<long long>(count));

    int wrong = 0;
    for (std::size_t i = 0; i < count; ++i) {
        const Position* p = world.get<Position>(entities[i]);
        const Velocity* v = world.get<Velocity>(entities[i]);
        if (p == nullptr || v == nullptr) { ++wrong; continue; }
        if (p->x != static_cast<double>(i) || v->x != static_cast<double>(i) * 2) ++wrong;
    }
    expectEqual("no value was corrupted across chunk boundaries", wrong, 0);
}

// A component declaring alignas(64) must actually get 64-byte-aligned storage,
// or an aligned SIMD load over the array faults.
void testComponentAlignmentInChunks() {
    World world;
    for (int i = 0; i < 200; ++i) world.create(Tag{static_cast<std::uint32_t>(i)}, Wide{});

    int misaligned = 0;
    world.each<Wide>([&](std::size_t, Wide* wide) {
        if ((reinterpret_cast<std::uintptr_t>(wide) & 63u) != 0) ++misaligned;
    });
    expectEqual("over-aligned component arrays are correctly aligned in the chunk", misaligned, 0);
}

// Iteration order must be a function of the operations performed, not of hash
// ordering or allocator addresses, or replays diverge.
void testIterationOrderIsDeterministic() {
    auto collect = [] {
        World world;
        for (int i = 0; i < 400; ++i) {
            const Entity e = world.create(Tag{static_cast<std::uint32_t>(i)});
            if (i % 3 == 0) world.add(e, Position{static_cast<double>(i), 0, 0});
            if (i % 5 == 0) world.add(e, Mass{static_cast<double>(i)});
        }
        std::vector<std::uint32_t> order;
        world.each<Tag>([&](std::size_t count, Tag* tags) {
            for (std::size_t i = 0; i < count; ++i) order.push_back(tags[i].value);
        });
        return order;
    };

    const auto first = collect();
    const auto second = collect();
    expectEqual("iteration visits every entity", static_cast<long long>(first.size()), 400);
    expectTrue("iteration order is identical across identical runs", first == second);

    // And every entity appears exactly once, whatever the order.
    std::set<std::uint32_t> unique(first.begin(), first.end());
    expectEqual("iteration yields no duplicates", static_cast<long long>(unique.size()), 400);
}

}  // namespace

void runEcsTests() {
    std::printf("\n--- archetype ECS ---\n");
    testEntityLifetime();
    testComponentAccess();
    testArchetypeTransitionsPreserveData();
    testSwapRemovePatchesMovedEntity();
    testEachVisitsMatchingEntitiesOnce();
    testMultipleChunks();
    testComponentAlignmentInChunks();
    testIterationOrderIsDeterministic();
}
