// SPDX-License-Identifier: MIT
//
// Validation of World save/load.
//
// A save format that loses data does not announce it. The specific quiet
// failures here: entity handles that stop resolving after a reload, generation
// counters that reset so a destroyed entity comes back to life, a component this
// build does not know taking the rest of the entity down with it, and a
// truncated file loading "successfully" into a half-populated world.
#include "engine/core/ecs.hpp"
#include "engine/core/reflect.hpp"
#include "engine/core/serialise.hpp"
#include "harness.hpp"

#include <cstdio>
#include <string>
#include <vector>

using core::ByteReader;
using core::ByteWriter;
using core::Entity;
using core::World;
using testing::expectEqual;
using testing::expectTrue;

// Reflected, so they serialise field-wise and survive schema change.
struct Hydrostatics {
    double displacement = 0;
    double draft = 0;
    std::int32_t compartments = 0;
};
struct Attitude {
    double heel = 0;
    double trim = 0;
};
// Deliberately unreflected: exercises the opaque-blob path.
struct OpaqueState {
    std::uint64_t bits[3]{};
};

SHIPSIM_REFLECT_BEGIN(Hydrostatics)
    SHIPSIM_FIELD(displacement)
    SHIPSIM_FIELD(draft)
    SHIPSIM_FIELD(compartments)
SHIPSIM_REFLECT_END(Hydrostatics)

SHIPSIM_REFLECT_BEGIN(Attitude)
    SHIPSIM_FIELD(heel)
    SHIPSIM_FIELD(trim)
SHIPSIM_REFLECT_END(Attitude)

namespace {

// Builds a world with a mixed population and a hole punched in the entity index
// space, so the free list and generation table both have something to say.
std::vector<Entity> populate(World& world, int count) {
    std::vector<Entity> entities;
    for (int i = 0; i < count; ++i) {
        Entity e = world.create(Hydrostatics{8984.0 + i, 5.5, 16});
        if (i % 2 == 0) world.add(e, Attitude{0.1 * i, -0.05 * i});
        if (i % 5 == 0) world.add(e, OpaqueState{{static_cast<std::uint64_t>(i), 7, 9}});
        entities.push_back(e);
    }
    // Destroy a scattered third of them, so the reload has a free list and stale
    // generations to reproduce.
    for (int i = 0; i < count; i += 3) world.destroy(entities[static_cast<std::size_t>(i)]);
    return entities;
}

void testRoundTripPreservesEverything() {
    World source;
    const auto entities = populate(source, 200);
    const std::size_t liveBefore = source.entityCount();

    ByteWriter writer;
    source.save(writer);
    expectTrue("a populated world produces a save", writer.size() > 0);

    World target;
    ByteReader reader(writer.span());
    expectTrue("the save loads", target.load(reader));
    expectTrue("the whole buffer was consumed", reader.remaining() == 0);
    expectEqual("live entity count survives", static_cast<long long>(target.entityCount()),
                static_cast<long long>(liveBefore));

    // Handles captured before the save must still resolve, to the same values.
    int wrong = 0, aliveMismatch = 0;
    for (std::size_t i = 0; i < entities.size(); ++i) {
        const Entity e = entities[i];
        if (source.alive(e) != target.alive(e)) { ++aliveMismatch; continue; }
        if (!target.alive(e)) continue;

        const auto* before = source.get<Hydrostatics>(e);
        const auto* after = target.get<Hydrostatics>(e);
        if (before == nullptr || after == nullptr ||
            before->displacement != after->displacement || before->draft != after->draft ||
            before->compartments != after->compartments)
            ++wrong;

        if (source.has<Attitude>(e) != target.has<Attitude>(e)) ++wrong;
        if (source.has<Attitude>(e) &&
            source.get<Attitude>(e)->heel != target.get<Attitude>(e)->heel)
            ++wrong;
        // Unreflected component, carried as an opaque blob.
        if (source.has<OpaqueState>(e) != target.has<OpaqueState>(e)) ++wrong;
        if (source.has<OpaqueState>(e) &&
            source.get<OpaqueState>(e)->bits[0] != target.get<OpaqueState>(e)->bits[0])
            ++wrong;
    }
    expectEqual("liveness matches for every handle", aliveMismatch, 0);
    expectEqual("every component value survives the round trip", wrong, 0);
}

// Generation counters are part of the world state: a handle that was stale
// before the save must still be stale after it, or destroyed entities come back.
void testHandleStalenessSurvives() {
    World source;
    const Entity kept = source.create(Hydrostatics{1, 2, 3});
    const Entity doomed = source.create(Hydrostatics{4, 5, 6});
    source.destroy(doomed);
    // Reusing the index bumps the generation; the old handle must stay dead.
    const Entity reused = source.create(Hydrostatics{7, 8, 9});
    expectEqual("the index was reused", static_cast<long long>(reused.index),
                static_cast<long long>(doomed.index));

    ByteWriter writer;
    source.save(writer);
    World target;
    ByteReader reader(writer.span());
    expectTrue("save loads", target.load(reader));

    expectTrue("a live handle stays live", target.alive(kept));
    expectTrue("a stale handle stays stale after reload", !target.alive(doomed));
    expectTrue("the handle that reused the index stays live", target.alive(reused));
    expectTrue("the reused index still holds its own data",
               target.get<Hydrostatics>(reused)->displacement == 7);

    // And the free list must come back, so the next create reuses correctly and
    // does not hand out an index that is already occupied.
    World fresh;
    ByteReader again(writer.span());
    fresh.load(again);
    const std::size_t before = fresh.entityCount();
    const Entity next = fresh.create();
    expectEqual("creating after load grows the population",
                static_cast<long long>(fresh.entityCount()),
                static_cast<long long>(before) + 1);
    expectTrue("the new entity does not collide with a live one",
               fresh.alive(next) && fresh.alive(kept) && fresh.alive(reused));
}

// A save containing a component this build has never heard of must load, and the
// entities carrying it must arrive with everything else intact.
void testUnknownComponentIsSkipped() {
    World source;
    const Entity e = source.create(Hydrostatics{100, 6, 8}, Attitude{0.25, -0.5});

    ByteWriter writer;
    source.save(writer);
    std::vector<std::byte> bytes(writer.bytes());

    // Repoint one component's stable id at a type nothing registers, which is
    // what a save from a build with an extra component looks like from here.
    const std::uint64_t attitudeId = core::typeInfo<Attitude>().stableId;
    const std::uint64_t unknownId = core::stableHash("AComponentThisBuildHasNever");
    bool patched = false;
    for (std::size_t i = 0; i + 8 <= bytes.size() && !patched; ++i) {
        std::uint64_t candidate = 0;
        for (int b = 0; b < 8; ++b)
            candidate |= static_cast<std::uint64_t>(std::to_integer<std::uint8_t>(bytes[i + b]))
                         << (8 * b);
        if (candidate != attitudeId) continue;
        for (int b = 0; b < 8; ++b)
            bytes[i + static_cast<std::size_t>(b)] =
                static_cast<std::byte>((unknownId >> (8 * b)) & 0xFFu);
        patched = true;
    }
    expectTrue("the component id was found and repointed", patched);

    World target;
    ByteReader reader(bytes);
    expectTrue("a save with an unknown component still loads", target.load(reader));
    expectEqual("the entity carrying it still arrives",
                static_cast<long long>(target.entityCount()), 1);
    expectTrue("the entity handle still resolves", target.alive(e));
    expectTrue("the known component is intact",
               target.get<Hydrostatics>(e) != nullptr &&
                   target.get<Hydrostatics>(e)->displacement == 100);
    expectTrue("the unknown component is simply absent", !target.has<Attitude>(e));
}

// Every truncation must be refused, and must leave the world empty rather than
// half-populated. ASan proves the reader stayed in bounds while refusing.
void testTruncatedSaveIsRejected() {
    World source;
    populate(source, 40);
    ByteWriter writer;
    source.save(writer);
    const std::vector<std::byte> full(writer.bytes());

    int accepted = 0, leftPopulated = 0;
    for (std::size_t length = 0; length < full.size(); length += 7) {
        World target;
        ByteReader reader(std::span<const std::byte>(full.data(), length));
        if (target.load(reader)) {
            ++accepted;
        } else if (target.entityCount() != 0) {
            // A refused load must not leave a partly built world behind.
            ++leftPopulated;
        }
    }
    expectEqual("no truncation of a valid save is accepted", accepted, 0);
    expectEqual("a refused load leaves no half-built world", leftPopulated, 0);

    // Wrong magic must be refused outright rather than misread.
    std::vector<std::byte> wrongMagic(full);
    wrongMagic[0] = std::byte{0x00};
    World target;
    ByteReader reader(wrongMagic);
    expectTrue("a foreign file is refused", !target.load(reader));
}

// Saving is a pure function of world state, and loading must not perturb it:
// save, load, save again must produce identical bytes. This is the property that
// makes saves diffable and makes a save usable as a network baseline.
void testSaveIsStableUnderRoundTrip() {
    World source;
    populate(source, 120);

    ByteWriter first;
    source.save(first);
    ByteWriter firstAgain;
    source.save(firstAgain);
    expectTrue("saving the same world twice is byte-identical",
               first.bytes() == firstAgain.bytes());

    World target;
    ByteReader reader(first.span());
    expectTrue("round trip loads", target.load(reader));

    ByteWriter second;
    target.save(second);
    expectTrue("save, load, save produces identical bytes", first.bytes() == second.bytes());

    // A third pass, to catch anything that only settles after one reload.
    World third;
    ByteReader secondReader(second.span());
    third.load(secondReader);
    ByteWriter thirdWrite;
    third.save(thirdWrite);
    expectTrue("the format reaches a fixed point immediately",
               second.bytes() == thirdWrite.bytes());
}

// Loading into a world that already holds data must replace it, not merge.
void testLoadReplacesExistingContents() {
    World source;
    source.create(Hydrostatics{1, 1, 1});

    ByteWriter writer;
    source.save(writer);

    World target;
    populate(target, 50);
    expectTrue("the target starts populated", target.entityCount() > 1);

    ByteReader reader(writer.span());
    expectTrue("load succeeds over existing contents", target.load(reader));
    expectEqual("load replaced rather than merged",
                static_cast<long long>(target.entityCount()), 1);
}

}  // namespace

void runWorldIoTests() {
    std::printf("\n--- world save/load ---\n");
    testRoundTripPreservesEverything();
    testHandleStalenessSurvives();
    testUnknownComponentIsSkipped();
    testTruncatedSaveIsRejected();
    testSaveIsStableUnderRoundTrip();
    testLoadReplacesExistingContents();
}
