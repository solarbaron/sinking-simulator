// SPDX-License-Identifier: MIT
#include "arena.hpp"

#include <algorithm>

// AddressSanitizer manual poisoning. A bump allocator is one large valid heap
// block as far as ASan is concerned, so without this an overrun past an
// allocation, or a read after reset(), is an invisible scribble inside memory
// ASan considers legitimately ours. Poisoning the unused tail turns both into
// ordinary reported errors.
#if defined(__SANITIZE_ADDRESS__)
#define SHIPSIM_HAS_ASAN 1
#elif defined(__has_feature)
#if __has_feature(address_sanitizer)
#define SHIPSIM_HAS_ASAN 1
#endif
#endif

#if defined(SHIPSIM_HAS_ASAN)
extern "C" void __asan_poison_memory_region(void const volatile* address, std::size_t size);
extern "C" void __asan_unpoison_memory_region(void const volatile* address, std::size_t size);
#define SHIPSIM_POISON(addr, size) __asan_poison_memory_region((addr), (size))
#define SHIPSIM_UNPOISON(addr, size) __asan_unpoison_memory_region((addr), (size))
#else
#define SHIPSIM_POISON(addr, size) ((void)(addr), (void)(size))
#define SHIPSIM_UNPOISON(addr, size) ((void)(addr), (void)(size))
#endif

namespace core {
namespace {

// Round `value` up to the next multiple of `alignment`, which must be a power of
// two. Returns SIZE_MAX on overflow so the caller can fail the allocation rather
// than wrap around into a small successful one.
std::size_t alignUp(std::size_t value, std::size_t alignment) {
    const std::size_t mask = alignment - 1;
    if (value > SIZE_MAX - mask) return SIZE_MAX;
    return (value + mask) & ~mask;
}

bool isPowerOfTwo(std::size_t value) { return value != 0 && (value & (value - 1)) == 0; }

}  // namespace

Arena::Arena(std::size_t capacityBytes) {
    if (capacityBytes == 0) return;
    const std::size_t rounded = alignUp(capacityBytes, kBlockAlignment);
    block_ = static_cast<std::byte*>(
        ::operator new(rounded, std::align_val_t{kBlockAlignment}));
    capacity_ = rounded;
    SHIPSIM_POISON(block_, capacity_);
}

Arena::~Arena() {
    if (block_ == nullptr) return;
    // Hand the block back unpoisoned; the allocator will want to write its own
    // bookkeeping into it.
    SHIPSIM_UNPOISON(block_, capacity_);
    ::operator delete(block_, std::align_val_t{kBlockAlignment});
}

Arena::Arena(Arena&& other) noexcept
    : block_(other.block_),
      capacity_(other.capacity_),
      offset_(other.offset_),
      highWaterMark_(other.highWaterMark_) {
    other.block_ = nullptr;
    other.capacity_ = 0;
    other.offset_ = 0;
    other.highWaterMark_ = 0;
}

Arena& Arena::operator=(Arena&& other) noexcept {
    if (this == &other) return *this;
    if (block_ != nullptr) {
        SHIPSIM_UNPOISON(block_, capacity_);
        ::operator delete(block_, std::align_val_t{kBlockAlignment});
    }
    block_ = other.block_;
    capacity_ = other.capacity_;
    offset_ = other.offset_;
    highWaterMark_ = other.highWaterMark_;
    other.block_ = nullptr;
    other.capacity_ = 0;
    other.offset_ = 0;
    other.highWaterMark_ = 0;
    return *this;
}

void* Arena::allocate(std::size_t bytes, std::size_t alignment) {
    if (block_ == nullptr || bytes == 0 || !isPowerOfTwo(alignment)) return nullptr;

    // Align the *address*, not the offset. Aligning the offset only works when
    // the block itself is at least as aligned as the request, so it silently
    // returns under-aligned memory for anything over kBlockAlignment -- which
    // then does not fault until someone does an aligned SIMD store to it.
    const auto base = reinterpret_cast<std::uintptr_t>(block_);
    const std::uintptr_t current = base + offset_;
    const std::uintptr_t mask = alignment - 1;
    if (current > UINTPTR_MAX - mask) return nullptr;
    const std::size_t aligned = static_cast<std::size_t>(((current + mask) & ~mask) - base);
    if (aligned > capacity_) return nullptr;
    if (bytes > capacity_ - aligned) return nullptr;  // exhausted, no wraparound

    std::byte* result = block_ + aligned;
    offset_ = aligned + bytes;
    highWaterMark_ = std::max(highWaterMark_, offset_);

    // Only the bytes actually handed out become writable. The alignment padding
    // stays poisoned, so an overrun off the end of a previous allocation is
    // caught rather than landing harmlessly in the gap.
    SHIPSIM_UNPOISON(result, bytes);
    return result;
}

void Arena::reset() {
    offset_ = 0;
    poisonFrom(0);
}

void Arena::rewindTo(std::size_t mark) {
    if (mark >= offset_) return;
    offset_ = mark;
    poisonFrom(mark);
}

void Arena::poisonFrom(std::size_t offset) {
    if (block_ == nullptr || offset >= capacity_) return;
    SHIPSIM_POISON(block_ + offset, capacity_ - offset);
}

Arena::Scope::~Scope() { arena_->rewindTo(mark_); }

// --- FrameArenas -------------------------------------------------------------

FrameArenas::FrameArenas(unsigned laneCount, std::size_t bytesPerLane) {
    arenas_.reserve(laneCount);
    for (unsigned i = 0; i < laneCount; ++i) arenas_.emplace_back(bytesPerLane);
}

void FrameArenas::resetAll() {
    for (Arena& arena : arenas_) arena.reset();
}

std::size_t FrameArenas::totalCapacity() const {
    std::size_t total = 0;
    for (const Arena& arena : arenas_) total += arena.capacity();
    return total;
}

std::size_t FrameArenas::totalHighWaterMark() const {
    std::size_t total = 0;
    for (const Arena& arena : arenas_) total += arena.highWaterMark();
    return total;
}

}  // namespace core
