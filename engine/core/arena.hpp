// SPDX-License-Identifier: MIT
//
// Bump allocators.
//
// The global allocator is off-limits during a tick (docs/01-architecture.md §5):
// it takes a lock, it is unpredictable, and a 100 Hz tick that occasionally
// stalls in malloc is a 100 Hz tick that occasionally misses. Everything
// transient allocates from an arena instead -- a pointer increment -- and the
// whole arena is released in one operation at the frame boundary.
//
// Three consequences worth stating up front, because they are the price:
//
//   * **Nothing is freed individually.** Only reset() and Scope release memory.
//   * **Destructors never run.** create<T>() therefore refuses anything that is
//     not trivially destructible, rather than silently leaking the destructor.
//   * **An Arena is not thread-safe.** That is the design, not an omission: one
//     arena per job-system lane means allocation needs no atomics at all. See
//     FrameArenas.
//
// Under AddressSanitizer the unused region is kept poisoned, so overruns and
// use-after-reset are caught inside the arena rather than being invisible
// scribbles within one big valid heap block.
#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <new>
#include <span>
#include <type_traits>
#include <vector>

namespace core {

class Arena {
public:
    // Default cache-line alignment for the underlying block, so per-lane arenas
    // never share a line and allocations start out favourably aligned.
    static constexpr std::size_t kBlockAlignment = 64;

    Arena() = default;
    explicit Arena(std::size_t capacityBytes);
    ~Arena();

    Arena(const Arena&) = delete;
    Arena& operator=(const Arena&) = delete;
    Arena(Arena&& other) noexcept;
    Arena& operator=(Arena&& other) noexcept;

    // Returns nullptr when the arena is exhausted rather than throwing or
    // growing, and also when `alignment` is not a power of two. Callers that
    // cannot handle that should check remaining() first; sizing is meant to come
    // from measured highWaterMark(), not from guesswork.
    //
    // Alignments larger than kBlockAlignment are honoured: the address is
    // aligned, not merely the offset within the block.
    void* allocate(std::size_t bytes, std::size_t alignment = alignof(std::max_align_t));

    template <typename T, typename... Args>
    T* create(Args&&... args);

    template <typename T>
    std::span<T> allocateArray(std::size_t count);

    // Releases everything at once. Pointers into the arena dangle afterwards.
    void reset();

    std::size_t capacity() const { return capacity_; }
    std::size_t used() const { return offset_; }
    std::size_t remaining() const { return capacity_ - offset_; }
    // Largest used() ever reached, across resets. This is the number that should
    // drive how big the arena is declared to be.
    std::size_t highWaterMark() const { return highWaterMark_; }

    // Restores the arena to its previous position on destruction. Nested scopes
    // must be destroyed in reverse order of creation, which RAII gives for free.
    class Scope {
    public:
        explicit Scope(Arena& arena) : arena_(&arena), mark_(arena.offset_) {}
        ~Scope();
        Scope(const Scope&) = delete;
        Scope& operator=(const Scope&) = delete;

    private:
        Arena* arena_;
        std::size_t mark_;
    };

private:
    friend class Scope;
    void rewindTo(std::size_t mark);
    void poisonFrom(std::size_t offset);

    std::byte*  block_ = nullptr;
    std::size_t capacity_ = 0;
    std::size_t offset_ = 0;
    std::size_t highWaterMark_ = 0;
};

// One arena per job-system lane, so job code allocates without contention and
// without atomics. Reset together at the frame boundary.
class FrameArenas {
public:
    FrameArenas() = default;
    FrameArenas(unsigned laneCount, std::size_t bytesPerLane);

    Arena& lane(unsigned index) { return arenas_[index]; }
    const Arena& lane(unsigned index) const { return arenas_[index]; }
    unsigned laneCount() const { return static_cast<unsigned>(arenas_.size()); }

    void resetAll();
    std::size_t totalCapacity() const;
    std::size_t totalHighWaterMark() const;

private:
    std::vector<Arena> arenas_;
};

// --- inline template definitions --------------------------------------------

template <typename T, typename... Args>
T* Arena::create(Args&&... args) {
    static_assert(std::is_trivially_destructible_v<T>,
                  "arena allocations are never destructed; use a trivially "
                  "destructible type, or manage the lifetime yourself");
    void* memory = allocate(sizeof(T), alignof(T));
    if (memory == nullptr) return nullptr;
    return ::new (memory) T(std::forward<Args>(args)...);
}

// The returned elements are **uninitialized**, like malloc and unlike vector.
// Array placement-new is deliberately avoided: the standard permits it to demand
// an unspecified amount of extra space for a cookie, which would silently
// overrun a bump allocator sized to sizeof(T) * count.
template <typename T>
std::span<T> Arena::allocateArray(std::size_t count) {
    static_assert(std::is_trivially_destructible_v<T>,
                  "arena allocations are never destructed; use a trivially "
                  "destructible type, or manage the lifetime yourself");
    static_assert(std::is_trivially_default_constructible_v<T>,
                  "arena arrays are handed back uninitialized; a type needing "
                  "construction must be built element by element with create<T>()");
    if (count == 0) return {};
    // Overflow guard: a caller computing count from simulation state should not
    // be able to wrap the byte count into a small successful allocation.
    if (count > SIZE_MAX / sizeof(T)) return {};
    void* memory = allocate(sizeof(T) * count, alignof(T));
    if (memory == nullptr) return {};
    return std::span<T>(static_cast<T*>(memory), count);
}

}  // namespace core
