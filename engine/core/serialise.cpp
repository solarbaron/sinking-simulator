// SPDX-License-Identifier: MIT
#include "serialise.hpp"

#include <cstring>

namespace core {
namespace {

// Explicit little-endian scalar codecs. memcpy would be faster and would make
// the format silently host-endian; a save file is meant to outlive the machine
// that wrote it.
void encodeLittleEndian(std::byte* out, std::uint64_t value, std::uint32_t width) {
    for (std::uint32_t i = 0; i < width; ++i)
        out[i] = static_cast<std::byte>((value >> (8 * i)) & 0xFFu);
}

std::uint64_t decodeLittleEndian(const std::byte* in, std::uint32_t width) {
    std::uint64_t value = 0;
    for (std::uint32_t i = 0; i < width; ++i)
        value |= static_cast<std::uint64_t>(std::to_integer<std::uint8_t>(in[i])) << (8 * i);
    return value;
}

// Reinterpret a scalar as its bit pattern. Floats go through their IEEE-754
// representation so they survive the byte-order round trip exactly, including
// negative zero and NaN payloads.
std::uint64_t scalarToBits(const void* source, TypeKind kind) {
    std::uint64_t bits = 0;
    switch (kind) {
        case TypeKind::Float32: {
            float value;
            std::memcpy(&value, source, sizeof(value));
            std::uint32_t raw;
            std::memcpy(&raw, &value, sizeof(raw));
            bits = raw;
            break;
        }
        case TypeKind::Float64: {
            double value;
            std::memcpy(&value, source, sizeof(value));
            std::memcpy(&bits, &value, sizeof(bits));
            break;
        }
        default:
            std::memcpy(&bits, source, kindByteWidth(kind));
            break;
    }
    return bits;
}

void bitsToScalar(void* target, TypeKind kind, std::uint64_t bits) {
    switch (kind) {
        case TypeKind::Float32: {
            const auto raw = static_cast<std::uint32_t>(bits);
            float value;
            std::memcpy(&value, &raw, sizeof(value));
            std::memcpy(target, &value, sizeof(value));
            break;
        }
        case TypeKind::Float64: {
            double value;
            std::memcpy(&value, &bits, sizeof(value));
            std::memcpy(target, &value, sizeof(value));
            break;
        }
        default:
            std::memcpy(target, &bits, kindByteWidth(kind));
            break;
    }
}

}  // namespace

// --- ByteWriter --------------------------------------------------------------

void ByteWriter::writeRaw(const void* data, std::size_t count) {
    const auto* source = static_cast<const std::byte*>(data);
    bytes_.insert(bytes_.end(), source, source + count);
}

void ByteWriter::writeU8(std::uint8_t value) { bytes_.push_back(static_cast<std::byte>(value)); }

void ByteWriter::writeU32(std::uint32_t value) {
    std::byte encoded[4];
    encodeLittleEndian(encoded, value, 4);
    writeRaw(encoded, 4);
}

void ByteWriter::writeU64(std::uint64_t value) {
    std::byte encoded[8];
    encodeLittleEndian(encoded, value, 8);
    writeRaw(encoded, 8);
}

void ByteWriter::patchU32(std::size_t offset, std::uint32_t value) {
    if (offset + 4 > bytes_.size()) return;
    encodeLittleEndian(bytes_.data() + offset, value, 4);
}

// --- ByteReader --------------------------------------------------------------

bool ByteReader::require(std::size_t count) {
    if (failed_) return false;
    if (count > remaining()) {
        failed_ = true;
        return false;
    }
    return true;
}

bool ByteReader::readRaw(void* out, std::size_t count) {
    if (!require(count)) return false;
    std::memcpy(out, bytes_.data() + cursor_, count);
    cursor_ += count;
    return true;
}

bool ByteReader::readU8(std::uint8_t& out) {
    if (!require(1)) return false;
    out = std::to_integer<std::uint8_t>(bytes_[cursor_]);
    ++cursor_;
    return true;
}

bool ByteReader::readU32(std::uint32_t& out) {
    if (!require(4)) return false;
    out = static_cast<std::uint32_t>(decodeLittleEndian(bytes_.data() + cursor_, 4));
    cursor_ += 4;
    return true;
}

bool ByteReader::readU64(std::uint64_t& out) {
    if (!require(8)) return false;
    out = decodeLittleEndian(bytes_.data() + cursor_, 8);
    cursor_ += 8;
    return true;
}

bool ByteReader::skip(std::size_t count) {
    if (!require(count)) return false;
    cursor_ += count;
    return true;
}

// --- object codec ------------------------------------------------------------

void serialiseObject(const TypeInfo& type, const void* object, ByteWriter& writer) {
    const auto* base = static_cast<const std::byte*>(object);
    writer.writeU64(type.stableId);
    writer.writeU32(type.fieldCount);

    for (std::uint32_t f = 0; f < type.fieldCount; ++f) {
        const FieldInfo& field = type.fields[f];
        writer.writeU64(field.nameHash);
        writer.writeU8(static_cast<std::uint8_t>(field.kind));
        writer.writeU32(field.count);

        // Length prefix, patched once the payload length is known. This is what
        // lets a future reader skip a field it has never heard of.
        const std::size_t lengthOffset = writer.size();
        writer.writeU32(0);
        const std::size_t payloadStart = writer.size();

        for (std::uint32_t i = 0; i < field.count; ++i) {
            const std::byte* element = base + field.offset + i * field.elementSize;
            if (field.kind == TypeKind::Struct) {
                if (field.nested != nullptr) serialiseObject(*field.nested, element, writer);
            } else {
                const std::uint32_t width = kindByteWidth(field.kind);
                std::byte encoded[8];
                encodeLittleEndian(encoded, scalarToBits(element, field.kind), width);
                writer.writeRaw(encoded, width);
            }
        }
        writer.patchU32(lengthOffset, static_cast<std::uint32_t>(writer.size() - payloadStart));
    }
}

bool deserialiseObject(const TypeInfo& type, void* object, ByteReader& reader) {
    auto* base = static_cast<std::byte*>(object);

    std::uint64_t streamTypeId = 0;
    std::uint32_t streamFieldCount = 0;
    if (!reader.readU64(streamTypeId) || !reader.readU32(streamFieldCount)) return false;
    if (streamTypeId != type.stableId) return false;  // wrong type entirely

    // A corrupt count must not be trusted as a loop bound: cap it at what the
    // remaining bytes could possibly describe. Each field costs at least its
    // 17-byte header.
    if (streamFieldCount > reader.remaining() / 17 + 1) return false;

    for (std::uint32_t f = 0; f < streamFieldCount; ++f) {
        std::uint64_t nameHash = 0;
        std::uint8_t rawKind = 0;
        std::uint32_t count = 0;
        std::uint32_t byteLength = 0;
        if (!reader.readU64(nameHash) || !reader.readU8(rawKind) || !reader.readU32(count) ||
            !reader.readU32(byteLength))
            return false;

        const auto kind = static_cast<TypeKind>(rawKind);
        const FieldInfo* field = type.findField(nameHash);

        // Unknown to this build, or changed shape since it was written: skip the
        // payload wholesale and leave the object's own value alone. This is the
        // whole point of the length prefix.
        const bool compatible = field != nullptr && field->kind == kind && field->count == count;
        if (!compatible) {
            if (!reader.skip(byteLength)) return false;
            continue;
        }

        const std::size_t payloadStart = reader.cursor();
        for (std::uint32_t i = 0; i < count; ++i) {
            std::byte* element = base + field->offset + i * field->elementSize;
            if (kind == TypeKind::Struct) {
                if (field->nested == nullptr) break;
                if (!deserialiseObject(*field->nested, element, reader)) return false;
            } else {
                const std::uint32_t width = kindByteWidth(kind);
                std::byte encoded[8];
                if (!reader.readRaw(encoded, width)) return false;
                bitsToScalar(element, kind, decodeLittleEndian(encoded, width));
            }
        }

        // Resynchronise on the recorded length even if the nested type has since
        // grown or shrunk, so one changed struct cannot desynchronise the rest of
        // the stream.
        const std::size_t consumed = reader.cursor() - payloadStart;
        if (consumed < byteLength) {
            if (!reader.skip(byteLength - consumed)) return false;
        } else if (consumed > byteLength) {
            return false;  // read past the field's own payload; the stream is corrupt
        }
    }
    return !reader.failed();
}

}  // namespace core
