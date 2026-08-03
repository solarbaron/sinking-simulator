// SPDX-License-Identifier: MIT
#include "ecs.hpp"

#include <algorithm>

namespace core {
namespace {

std::size_t alignUp(std::size_t value, std::size_t alignment) {
    return (value + alignment - 1) & ~(alignment - 1);
}

}  // namespace

ComponentRegistry& ComponentRegistry::instance() {
    static ComponentRegistry registry;
    return registry;
}

ComponentId ComponentRegistry::registerType(std::size_t size, std::size_t alignment,
                                            const char* name) {
    const auto id = static_cast<ComponentId>(infos_.size());
    infos_.push_back({size, alignment, name});
    return id;
}

World::~World() {
    for (auto& archetype : archetypes_)
        for (Chunk& chunk : archetype->chunks)
            ::operator delete(chunk.data, std::align_val_t{64});
}

int World::slotOf(const Archetype& archetype, ComponentId id) {
    const auto it = std::lower_bound(archetype.ids.begin(), archetype.ids.end(), id);
    if (it == archetype.ids.end() || *it != id) return -1;
    return static_cast<int>(it - archetype.ids.begin());
}

std::uint32_t World::findOrCreateArchetype(const std::vector<ComponentId>& sortedIds) {
    const auto existing = archetypeLookup_.find(sortedIds);
    if (existing != archetypeLookup_.end()) return existing->second;

    auto archetype = std::make_unique<Archetype>();
    archetype->ids = sortedIds;

    const ComponentRegistry& registry = ComponentRegistry::instance();

    // Solve for the largest entity count whose arrays fit in one chunk. Solving
    // it by trial rather than algebra keeps the alignment padding exact instead
    // of over-reserving for a worst case that never happens.
    std::size_t perEntity = sizeof(Entity);
    for (ComponentId id : sortedIds) perEntity += registry.info(id).size;

    std::size_t capacity = perEntity > 0 ? Chunk::kBytes / perEntity : 1;
    while (capacity > 0) {
        std::size_t cursor = alignUp(0, alignof(Entity));
        cursor += capacity * sizeof(Entity);
        bool fits = true;
        for (ComponentId id : sortedIds) {
            const ComponentInfo& info = registry.info(id);
            cursor = alignUp(cursor, info.alignment);
            cursor += capacity * info.size;
            if (cursor > Chunk::kBytes) { fits = false; break; }
        }
        if (fits && cursor <= Chunk::kBytes) break;
        --capacity;
    }
    if (capacity == 0) capacity = 1;  // a single enormous entity still gets a chunk

    // Now record the offsets for the capacity we settled on.
    archetype->entityOffset = 0;
    std::size_t cursor = capacity * sizeof(Entity);
    archetype->offsets.reserve(sortedIds.size());
    for (ComponentId id : sortedIds) {
        const ComponentInfo& info = registry.info(id);
        cursor = alignUp(cursor, info.alignment);
        archetype->offsets.push_back(cursor);
        cursor += capacity * info.size;
    }
    archetype->capacity = capacity;

    const auto index = static_cast<std::uint32_t>(archetypes_.size());
    archetypes_.push_back(std::move(archetype));
    archetypeLookup_.emplace(sortedIds, index);
    return index;
}

std::pair<std::uint32_t, std::uint32_t> World::appendRow(std::uint32_t archetypeIndex,
                                                         Entity entity) {
    Archetype& archetype = *archetypes_[archetypeIndex];

    std::uint32_t chunkIndex = 0;
    bool found = false;
    for (std::uint32_t c = 0; c < archetype.chunks.size(); ++c)
        if (archetype.chunks[c].count < archetype.capacity) {
            chunkIndex = c;
            found = true;
            break;
        }
    if (!found) {
        Chunk chunk;
        chunk.data = static_cast<std::byte*>(
            ::operator new(Chunk::kBytes, std::align_val_t{64}));
        chunk.count = 0;
        archetype.chunks.push_back(chunk);
        chunkIndex = static_cast<std::uint32_t>(archetype.chunks.size() - 1);
    }

    Chunk& chunk = archetype.chunks[chunkIndex];
    const std::uint32_t row = chunk.count++;
    auto* entities = reinterpret_cast<Entity*>(chunk.data + archetype.entityOffset);
    entities[row] = entity;
    return {chunkIndex, row};
}

void World::removeRow(std::uint32_t archetypeIndex, std::uint32_t chunkIndex, std::uint32_t row) {
    Archetype& archetype = *archetypes_[archetypeIndex];
    Chunk& chunk = archetype.chunks[chunkIndex];
    const std::uint32_t last = chunk.count - 1;

    if (row != last) {
        // Swap the last row down. Every component array and the entity array
        // must move together, and the entity that moved needs its record
        // patched -- forgetting that is the classic archetype ECS bug, and it
        // produces silent aliasing rather than a crash.
        auto* entities = reinterpret_cast<Entity*>(chunk.data + archetype.entityOffset);
        const Entity moved = entities[last];
        entities[row] = moved;

        const ComponentRegistry& registry = ComponentRegistry::instance();
        for (std::size_t slot = 0; slot < archetype.ids.size(); ++slot) {
            const ComponentInfo& info = registry.info(archetype.ids[slot]);
            std::byte* base = chunk.data + archetype.offsets[slot];
            std::memcpy(base + row * info.size, base + last * info.size, info.size);
        }

        Record& record = records_[moved.index];
        record.chunk = chunkIndex;
        record.row = row;
    }
    --chunk.count;
}

void World::moveEntity(Entity entity, std::uint32_t targetArchetype) {
    Record& record = records_[entity.index];
    const std::uint32_t sourceArchetype = record.archetype;
    if (sourceArchetype == targetArchetype) return;

    const std::uint32_t sourceChunk = record.chunk;
    const std::uint32_t sourceRow = record.row;

    const auto [chunkIndex, row] = appendRow(targetArchetype, entity);

    // Carry over everything the two archetypes have in common.
    const Archetype& source = *archetypes_[sourceArchetype];
    const Archetype& target = *archetypes_[targetArchetype];
    const ComponentRegistry& registry = ComponentRegistry::instance();
    for (std::size_t slot = 0; slot < source.ids.size(); ++slot) {
        const int targetSlot = slotOf(target, source.ids[slot]);
        if (targetSlot < 0) continue;  // dropped by this transition
        const ComponentInfo& info = registry.info(source.ids[slot]);
        const std::byte* from = archetypes_[sourceArchetype]->chunks[sourceChunk].data +
                                source.offsets[slot] + sourceRow * info.size;
        std::byte* to = archetypes_[targetArchetype]->chunks[chunkIndex].data +
                        target.offsets[static_cast<std::size_t>(targetSlot)] + row * info.size;
        std::memcpy(to, from, info.size);
    }

    // Update the record before removing the old row: removeRow may swap another
    // entity into that slot and patch *its* record, and the two must not fight.
    records_[entity.index].archetype = targetArchetype;
    records_[entity.index].chunk = chunkIndex;
    records_[entity.index].row = row;

    removeRow(sourceArchetype, sourceChunk, sourceRow);
}

void* World::componentPointer(std::uint32_t archetypeIndex, std::uint32_t chunkIndex,
                              std::uint32_t row, ComponentId id) {
    Archetype& archetype = *archetypes_[archetypeIndex];
    const int slot = slotOf(archetype, id);
    if (slot < 0) return nullptr;
    const ComponentInfo& info = ComponentRegistry::instance().info(id);
    return archetype.chunks[chunkIndex].data +
           archetype.offsets[static_cast<std::size_t>(slot)] + row * info.size;
}

const void* World::componentPointer(std::uint32_t archetypeIndex, std::uint32_t chunkIndex,
                                    std::uint32_t row, ComponentId id) const {
    const Archetype& archetype = *archetypes_[archetypeIndex];
    const int slot = slotOf(archetype, id);
    if (slot < 0) return nullptr;
    const ComponentInfo& info = ComponentRegistry::instance().info(id);
    return archetype.chunks[chunkIndex].data +
           archetype.offsets[static_cast<std::size_t>(slot)] + row * info.size;
}

Entity World::create() {
    const std::uint32_t emptyArchetype = findOrCreateArchetype({});

    std::uint32_t index;
    if (!freeIndices_.empty()) {
        index = freeIndices_.back();
        freeIndices_.pop_back();
    } else {
        index = static_cast<std::uint32_t>(records_.size());
        records_.push_back({});
    }

    Entity entity{index, records_[index].generation};
    const auto [chunkIndex, row] = appendRow(emptyArchetype, entity);

    Record& record = records_[index];
    record.archetype = emptyArchetype;
    record.chunk = chunkIndex;
    record.row = row;
    record.alive = true;
    ++liveEntities_;
    return entity;
}

void World::destroy(Entity entity) {
    if (!alive(entity)) return;
    Record& record = records_[entity.index];
    removeRow(record.archetype, record.chunk, record.row);
    record.alive = false;
    // Bump the generation so any surviving handle is now detectably stale.
    ++record.generation;
    freeIndices_.push_back(entity.index);
    --liveEntities_;
}

bool World::alive(Entity entity) const {
    if (entity.index >= records_.size()) return false;
    const Record& record = records_[entity.index];
    return record.alive && record.generation == entity.generation;
}

}  // namespace core
