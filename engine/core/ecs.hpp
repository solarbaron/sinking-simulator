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

#include "reflect.hpp"
#include "serialise.hpp"

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
    // Hash of the type's name. Unlike the runtime ComponentId this does not
    // depend on which type happened to be touched first.
    std::uint64_t stableId = 0;
    // Reflection for the component, when it has any. Null means the component
    // can still be stored and queried, but a save can only carry it as an opaque
    // blob that no other build can be trusted to interpret.
    const TypeInfo* type = nullptr;
};

// Process-wide type registry. A ComponentId is only an array index: it is handed
// out on first use, so its numeric value depends on which type was touched first
// and is meaningless outside the run that produced it. Everything that must be
// reproducible -- archetype ordering, and therefore iteration order -- keys off
// ComponentInfo::stableId instead, which is a hash of the type's name.
class ComponentRegistry {
public:
    static ComponentRegistry& instance();

    ComponentId registerType(std::size_t size, std::size_t alignment, const char* name,
                             std::uint64_t stableId, const TypeInfo* type);
    // Runtime id for a saved stable id, or kUnknownComponent when this build has
    // no such component. Loading must cope with the second case rather than
    // refusing the whole save.
    static constexpr ComponentId kUnknownComponent = 0xFFFF;
    ComponentId findByStableId(std::uint64_t stableId) const;
    // Orders two components by stable id, giving archetypes a canonical key that
    // is identical across builds and link orders.
    bool lessByStableId(ComponentId a, ComponentId b) const {
        return infos_[a].stableId < infos_[b].stableId;
    }
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
    // A reflected component names itself; anything else falls back to the
    // compiler's mangled name, which is stable within a toolchain but not across
    // them. Reflecting a component is what makes its data portable.
    static const ComponentId id = [] {
        if constexpr (Reflected<T>) {
            const TypeInfo& type = Reflect<T>::info();
            return ComponentRegistry::instance().registerType(sizeof(T), alignof(T), type.name,
                                                              type.stableId, &type);
        } else {
            return ComponentRegistry::instance().registerType(
                sizeof(T), alignof(T), typeid(T).name(), stableHash(typeid(T).name()), nullptr);
        }
    }();
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
    // Component ids of an archetype, in canonical (stable id) order. Exposed for
    // the editor and for tests that need to check the ordering invariant.
    const std::vector<ComponentId>& archetypeIds(std::size_t index) const;

    // --- persistence ---
    // Writes every entity, its components and the handle bookkeeping needed for
    // saved Entity values to keep resolving. Reflected components are written
    // field by field and survive schema changes; unreflected ones go out as
    // size-tagged opaque blobs, which is build-local and refused on load if the
    // size no longer matches.
    void save(ByteWriter& writer) const;
    // Replaces the world's contents. Returns false on a malformed stream,
    // leaving the world empty rather than half-loaded. Components this build
    // does not know are skipped; the entities carrying them still load.
    bool load(ByteReader& reader);
    void clear();
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

    // The decoding half of load(). Kept separate so load() can guarantee that a
    // failure leaves nothing behind, whatever point it failed at.
    bool loadImpl(ByteReader& reader);
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

// Canonical component order for an archetype. Sorting by stable id rather than
// by the runtime index is what keeps archetype identity -- and therefore query
// iteration order -- the same from one build to the next.
inline std::vector<ComponentId> sortedIds(std::vector<ComponentId> ids) {
    const ComponentRegistry& registry = ComponentRegistry::instance();
    std::sort(ids.begin(), ids.end(),
              [&](ComponentId a, ComponentId b) { return registry.lessByStableId(a, b); });
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
