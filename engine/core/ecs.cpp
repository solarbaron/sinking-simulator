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
                                            const char* name, std::uint64_t stableId,
                                            const TypeInfo* type) {
    const auto id = static_cast<ComponentId>(infos_.size());
    infos_.push_back({size, alignment, name, stableId, type});
    return id;
}

ComponentId ComponentRegistry::findByStableId(std::uint64_t stableId) const {
    for (std::size_t i = 0; i < infos_.size(); ++i)
        if (infos_[i].stableId == stableId) return static_cast<ComponentId>(i);
    return kUnknownComponent;
}

World::~World() {
    for (auto& archetype : archetypes_)
        for (Chunk& chunk : archetype->chunks)
            ::operator delete(chunk.data, std::align_val_t{64});
}

int World::slotOf(const Archetype& archetype, ComponentId id) {
    // Linear rather than binary: the list is sorted by stable id, not by the
    // value being searched for, and archetypes hold few enough components that a
    // scan is faster than carrying the comparator around anyway.
    for (std::size_t i = 0; i < archetype.ids.size(); ++i)
        if (archetype.ids[i] == id) return static_cast<int>(i);
    return -1;
}

const std::vector<ComponentId>& World::archetypeIds(std::size_t index) const {
    return archetypes_[index]->ids;
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

// --- persistence -------------------------------------------------------------
//
// Format, all little-endian via ByteWriter:
//
//   magic u32, version u32
//   recordCount u32, then that many generation u32
//   freeCount u32, then that many index u32
//   entityCount u32, then per entity:
//     index u32, generation u32, componentCount u32
//     per component: stableId u64, reflected u8, byteLength u32, payload
//
// The generation table comes first because handle staleness is part of the world
// state: an Entity saved and then destroyed must still read as dead after a
// reload, which means the generation counters of *unoccupied* indices matter too.

namespace {
constexpr std::uint32_t kWorldMagic = 0x57504853u;  // "SHPW"
constexpr std::uint32_t kWorldVersion = 1u;
}  // namespace

void World::save(ByteWriter& writer) const {
    const ComponentRegistry& registry = ComponentRegistry::instance();

    writer.writeU32(kWorldMagic);
    writer.writeU32(kWorldVersion);

    writer.writeU32(static_cast<std::uint32_t>(records_.size()));
    for (const Record& record : records_) writer.writeU32(record.generation);

    writer.writeU32(static_cast<std::uint32_t>(freeIndices_.size()));
    for (std::uint32_t index : freeIndices_) writer.writeU32(index);

    writer.writeU32(static_cast<std::uint32_t>(liveEntities_));

    // Archetype then chunk then row: the same order queries use, so a saved
    // world reloads with its iteration order intact.
    for (const auto& archetype : archetypes_) {
        for (const Chunk& chunk : archetype->chunks) {
            const auto* entities =
                reinterpret_cast<const Entity*>(chunk.data + archetype->entityOffset);
            for (std::uint32_t row = 0; row < chunk.count; ++row) {
                writer.writeU32(entities[row].index);
                writer.writeU32(entities[row].generation);
                writer.writeU32(static_cast<std::uint32_t>(archetype->ids.size()));

                for (std::size_t slot = 0; slot < archetype->ids.size(); ++slot) {
                    const ComponentInfo& info = registry.info(archetype->ids[slot]);
                    const std::byte* source =
                        chunk.data + archetype->offsets[slot] + row * info.size;

                    writer.writeU64(info.stableId);
                    writer.writeU8(info.type != nullptr ? 1u : 0u);
                    const std::size_t lengthOffset = writer.size();
                    writer.writeU32(0);
                    const std::size_t payloadStart = writer.size();

                    if (info.type != nullptr) {
                        serialiseObject(*info.type, source, writer);
                    } else {
                        // Opaque: only this build can interpret it, so record the
                        // size and let the reader refuse a mismatch.
                        writer.writeU32(static_cast<std::uint32_t>(info.size));
                        writer.writeRaw(source, info.size);
                    }
                    writer.patchU32(lengthOffset,
                                    static_cast<std::uint32_t>(writer.size() - payloadStart));
                }
            }
        }
    }
}

void World::clear() {
    for (auto& archetype : archetypes_)
        for (Chunk& chunk : archetype->chunks)
            ::operator delete(chunk.data, std::align_val_t{64});
    archetypes_.clear();
    archetypeLookup_.clear();
    records_.clear();
    freeIndices_.clear();
    liveEntities_ = 0;
}

bool World::load(ByteReader& reader) {
    clear();
    if (loadImpl(reader)) return true;
    // A partly decoded world is worse than no world: a caller that ignores the
    // return value would then be simulating a ship missing an arbitrary suffix
    // of its entities. Fail closed.
    clear();
    return false;
}

bool World::loadImpl(ByteReader& reader) {
    const ComponentRegistry& registry = ComponentRegistry::instance();

    std::uint32_t magic = 0, version = 0;
    if (!reader.readU32(magic) || !reader.readU32(version)) return false;
    if (magic != kWorldMagic || version != kWorldVersion) return false;

    std::uint32_t recordCount = 0;
    if (!reader.readU32(recordCount)) return false;
    // A count is only believable if the bytes to back it exist. Trusting it as a
    // resize argument is how a corrupt file becomes an allocation failure.
    if (recordCount > reader.remaining() / sizeof(std::uint32_t)) return false;
    records_.resize(recordCount);
    for (std::uint32_t i = 0; i < recordCount; ++i) {
        if (!reader.readU32(records_[i].generation)) return false;
        records_[i].alive = false;
    }

    // **The free list's *values* are bounds-checked, not just its length.**
    //
    // This validated `freeCount` against the bytes remaining and then read the
    // indices without looking at them, while the entity loop twenty lines below
    // has always checked `index >= recordCount`, rejected a generation mismatch
    // and rejected a repeated index. The asymmetry was the bug: a save declaring
    // four records and one free index of 0xFFFFFFF0 **loaded successfully**, and
    // the next `World::create()` took that index off the back of the list and
    // used it to subscript `records_` — a read at `records_[0xFFFFFFF0]` and then
    // a write to `.alive`.
    //
    // It does not crash, which is what makes it worse than a crash: the index
    // lands inside the vector's own heap block, so **ASan does not trap** and
    // `create()` returns normally onto a corrupted entity table. That is the
    // "failed open" shape this loader was once fixed for, reached through a door
    // the fix did not cover — the hardening was against *truncation*, and this
    // takes a well-formed file with a wrong number in it.
    //
    // A duplicate is refused for the same reason a repeated entity index is: two
    // entities would be handed one `Record`, and the second `create()` would
    // silently alias the first.
    std::uint32_t freeCount = 0;
    if (!reader.readU32(freeCount)) return false;
    if (freeCount > reader.remaining() / sizeof(std::uint32_t)) return false;
    if (freeCount > recordCount) return false;  // cannot free more than exist
    freeIndices_.resize(freeCount);
    std::vector<bool> freed(recordCount, false);
    for (std::uint32_t i = 0; i < freeCount; ++i) {
        if (!reader.readU32(freeIndices_[i])) return false;
        const std::uint32_t index = freeIndices_[i];
        if (index >= recordCount) return false;  // would subscript out of range
        if (freed[index]) return false;          // the same index twice
        freed[index] = true;
    }

    std::uint32_t entityCount = 0;
    if (!reader.readU32(entityCount)) return false;
    if (entityCount > recordCount) return false;  // cannot have more live than records

    std::vector<ComponentId> ids;
    for (std::uint32_t e = 0; e < entityCount; ++e) {
        std::uint32_t index = 0, generation = 0, componentCount = 0;
        if (!reader.readU32(index) || !reader.readU32(generation) ||
            !reader.readU32(componentCount))
            return false;
        if (index >= recordCount) return false;
        if (records_[index].generation != generation) return false;
        if (records_[index].alive) return false;  // the same index twice
        // ...and not one the free list has already claimed. Without this a
        // well-formed-looking save can put an index in both lists: the entity
        // loads live, and the next `create()` pops the same index off the free
        // list and hands a second entity the same `Record`.
        if (freed[index]) return false;
        // Each component costs at least its 13-byte header.
        if (componentCount > reader.remaining() / 13 + 1) return false;

        // First pass: read the headers, note where each payload sits, and learn
        // which components this build can actually place. The archetype has to
        // be known in full before a row exists to decode into.
        struct SavedComponent {
            std::uint64_t stableId = 0;
            bool reflected = false;
            std::size_t offset = 0;
            std::size_t length = 0;
        };
        std::vector<SavedComponent> saved;
        saved.reserve(componentCount);
        ids.clear();

        for (std::uint32_t c = 0; c < componentCount; ++c) {
            SavedComponent component;
            std::uint8_t reflected = 0;
            std::uint32_t byteLength = 0;
            if (!reader.readU64(component.stableId) || !reader.readU8(reflected) ||
                !reader.readU32(byteLength))
                return false;
            component.reflected = reflected != 0;
            component.offset = reader.cursor();
            component.length = byteLength;
            if (!reader.skip(byteLength)) return false;
            saved.push_back(component);

            const ComponentId id = registry.findByStableId(component.stableId);
            if (id == ComponentRegistry::kUnknownComponent) continue;  // not in this build
            const ComponentInfo& info = registry.info(id);
            // A component that gained or lost reflection since the save cannot be
            // decoded either way round; skip it rather than guess.
            if (component.reflected != (info.type != nullptr)) continue;
            ids.push_back(id);
        }

        const std::uint32_t archetype = findOrCreateArchetype(detail::sortedIds(ids));
        const Entity entity{index, generation};
        const auto [chunkIndex, row] = appendRow(archetype, entity);
        records_[index].archetype = archetype;
        records_[index].chunk = chunkIndex;
        records_[index].row = row;
        records_[index].alive = true;
        ++liveEntities_;

        // Second pass: decode the payloads that have somewhere to land.
        for (const SavedComponent& component : saved) {
            const ComponentId id = registry.findByStableId(component.stableId);
            if (id == ComponentRegistry::kUnknownComponent) continue;
            const ComponentInfo& info = registry.info(id);
            if (component.reflected != (info.type != nullptr)) continue;

            void* target = componentPointer(archetype, chunkIndex, row, id);
            if (target == nullptr) continue;

            ByteReader payload = reader.sliceAt(component.offset, component.length);
            if (info.type != nullptr) {
                if (!deserialiseObject(*info.type, target, payload)) return false;
            } else {
                std::uint32_t blobSize = 0;
                if (!payload.readU32(blobSize)) return false;
                // Refuse a blob whose size no longer matches the component:
                // reading it in would be a silent reinterpretation.
                if (blobSize != info.size) continue;
                if (!payload.readRaw(target, blobSize)) return false;
            }
        }
    }

    return !reader.failed();
}

bool World::alive(Entity entity) const {
    if (entity.index >= records_.size()) return false;
    const Record& record = records_[entity.index];
    return record.alive && record.generation == entity.generation;
}

}  // namespace core
