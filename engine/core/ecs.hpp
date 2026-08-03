// SPDX-License-Identifier: MIT
//
// Archetype ECS with structure-of-arrays chunk storage.
//
// Chosen over sparse sets because the hot loops in this engine are numerical
// sweeps over homogeneous collections -- integrate every tetrahedron, evaluate
// every orifice, advance every fire cell. Archetype chunks give those loops
// contiguous arrays with no indirection, which is the whole point; the price is
// that adding or removing a component moves the entity's row to another
// archetype, so structural changes are comparatively expensive and belong at
// frame boundaries rather than in inner loops.
//
// Components are plain data. They are moved between chunks with memcpy, so they
// must be trivially copyable, and they are never destructed.
#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <map>
#include <memory>
#include <string>
#include <type_traits>
#include <typeinfo>
#include <vector>

namespace core {

// Index plus generation, so a handle to a destroyed entity is detectably stale
// rather than silently aliasing whatever was created in its place.
struct Entity {
    std::uint32_t index = 0;
    std::uint32_t generation = 0;

    friend bool operator==(const Entity&, const Entity&) = default;
};

inline constexpr Entity kNoEntity{0xFFFFFFFFu, 0u};

using ComponentId = std::uint16_t;

struct ComponentInfo {
    std::size_t size = 0;
    std::size_t alignment = 0;
    std::string name;
};

// Process-wide type registry. Ids are handed out on first use, so they are
// stable within a run but their *values* depend on which type is touched first.
// Nothing may depend on the numeric value; archetypes sort by id purely to get a
// canonical key. Stable cross-run ids are a job for the reflection pass.
class ComponentRegistry {
public:
    static ComponentRegistry& instance();

    ComponentId registerType(std::size_t size, std::size_t alignment, const char* name);
    const ComponentInfo& info(ComponentId id) const { return infos_[id]; }
    std::size_t count() const { return infos_.size(); }

private:
    std::vector<ComponentInfo> infos_;
};

template <typename T>
ComponentId componentId() {
    static_assert(std::is_trivially_copyable_v<T>,
                  "components are relocated between chunks with memcpy");
    static_assert(std::is_trivially_destructible_v<T>,
                  "components are never destructed");
    static const ComponentId id =
        ComponentRegistry::instance().registerType(sizeof(T), alignof(T), typeid(T).name());
    return id;
}

// A fixed-size block holding the SoA arrays for one archetype. 16 KB is a
// compromise: large enough that per-chunk overhead disappears, small enough that
// a chunk plus the arrays a loop touches stays inside L2.
struct Chunk {
    static constexpr std::size_t kBytes = 16 * 1024;
    std::byte* data = nullptr;
    std::uint32_t count = 0;
};

class World {
public:
    World() = default;
    ~World();
    World(const World&) = delete;
    World& operator=(const World&) = delete;

    // --- entity lifetime ---
    Entity create();
    template <typename... Cs>
    Entity create(const Cs&... components);
    void destroy(Entity entity);
    bool alive(Entity entity) const;

    // --- components ---
    template <typename T>
    bool has(Entity entity) const;
    template <typename T>
    T* get(Entity entity);
    template <typename T>
    const T* get(Entity entity) const;
    template <typename T>
    void add(Entity entity, const T& value);
    template <typename T>
    void remove(Entity entity);

    // Visit every chunk holding all of Cs. `fn` is called as
    // fn(count, Cs*...) with contiguous arrays -- the shape a numerical sweep
    // wants. Iteration order is archetype creation order then chunk order, which
    // is deterministic for a given sequence of operations.
    template <typename... Cs, typename F>
    void each(F&& fn);

    std::size_t entityCount() const { return liveEntities_; }
    std::size_t archetypeCount() const { return archetypes_.size(); }
    // Entities per chunk for the archetype holding exactly Cs. Exposed so tests
    // can drive the multi-chunk paths rather than guessing at the threshold.
    template <typename... Cs>
    std::size_t chunkCapacityFor();

private:
    struct Archetype {
        std::vector<ComponentId> ids;    // sorted, the canonical key
        std::vector<std::size_t> offsets;  // byte offset of each array in a chunk
        std::size_t entityOffset = 0;      // byte offset of the Entity array
        std::size_t capacity = 0;          // entities per chunk
        std::vector<Chunk> chunks;
    };

    struct Record {
        std::uint32_t archetype = 0;
        std::uint32_t chunk = 0;
        std::uint32_t row = 0;
        std::uint32_t generation = 0;
        bool alive = false;
    };

    std::uint32_t findOrCreateArchetype(const std::vector<ComponentId>& sortedIds);
    // Appends a row and returns (chunk, row). Allocates a chunk if needed.
    std::pair<std::uint32_t, std::uint32_t> appendRow(std::uint32_t archetypeIndex, Entity entity);
    // Removes a row by swapping the last row into its place, patching the record
    // of whichever entity moved.
    void removeRow(std::uint32_t archetypeIndex, std::uint32_t chunkIndex, std::uint32_t row);
    // Relocates an entity to a different archetype, carrying over every component
    // the two have in common.
    void moveEntity(Entity entity, std::uint32_t targetArchetype);

    void* componentPointer(std::uint32_t archetypeIndex, std::uint32_t chunkIndex,
                           std::uint32_t row, ComponentId id);
    const void* componentPointer(std::uint32_t archetypeIndex, std::uint32_t chunkIndex,
                                 std::uint32_t row, ComponentId id) const;
    // Index of `id` within the archetype's sorted id list, or -1.
    static int slotOf(const Archetype& archetype, ComponentId id);

    std::vector<Record> records_;
    std::vector<std::uint32_t> freeIndices_;
    std::vector<std::unique_ptr<Archetype>> archetypes_;
    // Ordered, not unordered: archetype creation order must not depend on hash
    // iteration if the engine is to stay reproducible.
    std::map<std::vector<ComponentId>, std::uint32_t> archetypeLookup_;
    std::size_t liveEntities_ = 0;
};

// --- inline template definitions --------------------------------------------

namespace detail {

inline std::vector<ComponentId> sortedIds(std::vector<ComponentId> ids) {
    std::sort(ids.begin(), ids.end());
    ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
    return ids;
}

}  // namespace detail

template <typename... Cs>
Entity World::create(const Cs&... components) {
    const auto ids = detail::sortedIds({componentId<Cs>()...});
    const std::uint32_t archetype = findOrCreateArchetype(ids);

    Entity entity = create();  // starts in the empty archetype
    moveEntity(entity, archetype);
    // Fold over the pack, writing each component into its array.
    (..., (void)std::memcpy(componentPointer(records_[entity.index].archetype,
                                             records_[entity.index].chunk,
                                             records_[entity.index].row, componentId<Cs>()),
                            &components, sizeof(Cs)));
    return entity;
}

template <typename T>
bool World::has(Entity entity) const {
    if (!alive(entity)) return false;
    const Record& record = records_[entity.index];
    return slotOf(*archetypes_[record.archetype], componentId<T>()) >= 0;
}

template <typename T>
T* World::get(Entity entity) {
    if (!alive(entity)) return nullptr;
    const Record& record = records_[entity.index];
    return static_cast<T*>(
        componentPointer(record.archetype, record.chunk, record.row, componentId<T>()));
}

template <typename T>
const T* World::get(Entity entity) const {
    if (!alive(entity)) return nullptr;
    const Record& record = records_[entity.index];
    return static_cast<const T*>(
        componentPointer(record.archetype, record.chunk, record.row, componentId<T>()));
}

template <typename T>
void World::add(Entity entity, const T& value) {
    if (!alive(entity)) return;
    const ComponentId id = componentId<T>();
    Record& record = records_[entity.index];
    if (slotOf(*archetypes_[record.archetype], id) < 0) {
        std::vector<ComponentId> ids = archetypes_[record.archetype]->ids;
        ids.push_back(id);
        moveEntity(entity, findOrCreateArchetype(detail::sortedIds(std::move(ids))));
    }
    const Record& moved = records_[entity.index];
    std::memcpy(componentPointer(moved.archetype, moved.chunk, moved.row, id), &value, sizeof(T));
}

template <typename T>
void World::remove(Entity entity) {
    if (!alive(entity)) return;
    const ComponentId id = componentId<T>();
    const Record& record = records_[entity.index];
    if (slotOf(*archetypes_[record.archetype], id) < 0) return;

    std::vector<ComponentId> ids = archetypes_[record.archetype]->ids;
    ids.erase(std::remove(ids.begin(), ids.end(), id), ids.end());
    moveEntity(entity, findOrCreateArchetype(ids));
}

template <typename... Cs, typename F>
void World::each(F&& fn) {
    const ComponentId wanted[] = {componentId<Cs>()...};
    for (const auto& archetype : archetypes_) {
        bool matches = true;
        for (ComponentId id : wanted)
            if (slotOf(*archetype, id) < 0) matches = false;
        if (!matches) continue;

        for (std::uint32_t c = 0; c < archetype->chunks.size(); ++c) {
            Chunk& chunk = archetype->chunks[c];
            if (chunk.count == 0) continue;
            fn(static_cast<std::size_t>(chunk.count),
               reinterpret_cast<Cs*>(
                   chunk.data + archetype->offsets[static_cast<std::size_t>(
                                    slotOf(*archetype, componentId<Cs>()))])...);
        }
    }
}

template <typename... Cs>
std::size_t World::chunkCapacityFor() {
    const auto ids = detail::sortedIds({componentId<Cs>()...});
    return archetypes_[findOrCreateArchetype(ids)]->capacity;
}

}  // namespace core
