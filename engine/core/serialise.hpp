// SPDX-License-Identifier: MIT
//
// Reflection-driven serialisation.
//
// The format is self-describing per field: name hash, kind, count and byte
// length. That costs bytes, and buys the property that actually matters for a
// project meant to run for years — **a save written by one build can be read by
// another**. Fields the reader does not know are skipped by their recorded
// length; fields the writer did not have keep whatever the reader's object
// already holds. Neither side has to agree on field order.
//
// Scalars are written little-endian explicitly rather than memcpy'd, so a save
// file is portable rather than accidentally x86-shaped.
//
// The reader never trusts the stream: every read is bounds-checked and a
// truncated or corrupt buffer fails cleanly instead of walking off the end.
#pragma once

#include "reflect.hpp"

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace core {

class ByteWriter {
public:
    void writeRaw(const void* data, std::size_t bytes);
    void writeU8(std::uint8_t value);
    void writeU32(std::uint32_t value);
    void writeU64(std::uint64_t value);
    // Patches a previously written u32 in place, for length prefixes whose value
    // is not known until the payload has been written.
    void patchU32(std::size_t offset, std::uint32_t value);

    std::size_t size() const { return bytes_.size(); }
    const std::vector<std::byte>& bytes() const { return bytes_; }
    std::span<const std::byte> span() const { return bytes_; }

private:
    std::vector<std::byte> bytes_;
};

class ByteReader {
public:
    explicit ByteReader(std::span<const std::byte> bytes) : bytes_(bytes) {}

    bool readRaw(void* out, std::size_t count);
    bool readU8(std::uint8_t& out);
    bool readU32(std::uint32_t& out);
    bool readU64(std::uint64_t& out);
    bool skip(std::size_t count);

    bool failed() const { return failed_; }
    std::size_t remaining() const { return bytes_.size() - cursor_; }
    std::size_t cursor() const { return cursor_; }

private:
    bool require(std::size_t count);

    std::span<const std::byte> bytes_;
    std::size_t cursor_ = 0;
    bool failed_ = false;
};

// Type-erased forms; the templates below are thin wrappers.
void serialiseObject(const TypeInfo& type, const void* object, ByteWriter& writer);
bool deserialiseObject(const TypeInfo& type, void* object, ByteReader& reader);

template <Reflected T>
void serialise(const T& value, ByteWriter& writer) {
    serialiseObject(typeInfo<T>(), &value, writer);
}

template <Reflected T>
bool deserialise(T& value, ByteReader& reader) {
    return deserialiseObject(typeInfo<T>(), &value, reader);
}

}  // namespace core
