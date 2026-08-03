// SPDX-License-Identifier: MIT
//
// Validation of reflection and serialisation.
//
// The failures that matter here are the quiet ones. A save format that silently
// drops a field, or reads a renamed one into the wrong slot, or walks off the
// end of a truncated buffer, produces a world that looks loaded. So: schema
// evolution is tested in both directions, and every malformed-input case is
// asserted to fail cleanly rather than merely "not crash on my machine" —
// AddressSanitizer is what proves the second half of that.
#include "engine/core/reflect.hpp"
#include "engine/core/serialise.hpp"
#include "harness.hpp"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

using core::ByteReader;
using core::ByteWriter;
using testing::expectEqual;
using testing::expectNear;
using testing::expectTrue;

namespace {

struct Vec3f {
    float x = 0, y = 0, z = 0;
};

struct Hull {
    double displacement = 0;
    Vec3f centre;
    std::int32_t compartments = 0;
    bool afloat = false;
    double offsets[4]{};
};

// Two shapes of the "same" record, standing in for the same type across two
// builds: `Evolved` drops `beta`, adds `delta`, and reorders what is left.
struct Original {
    std::int32_t alpha = 0;
    double beta = 0;
    std::uint64_t gamma = 0;
};

struct Evolved {
    std::uint64_t gamma = 0;
    std::int32_t alpha = 0;
    float delta = 0;
};

}  // namespace

SHIPSIM_REFLECT_BEGIN(Vec3f)
    SHIPSIM_FIELD(x)
    SHIPSIM_FIELD(y)
    SHIPSIM_FIELD(z)
SHIPSIM_REFLECT_END(Vec3f)

SHIPSIM_REFLECT_BEGIN(Hull)
    SHIPSIM_FIELD(displacement)
    SHIPSIM_FIELD(centre)
    SHIPSIM_FIELD(compartments)
    SHIPSIM_FIELD(afloat)
    SHIPSIM_FIELD(offsets)
SHIPSIM_REFLECT_END(Hull)

// Both spellings deliberately reflect under the *same* type name, because that
// is what "the same type in two builds" means to the format.
SHIPSIM_REFLECT_BEGIN(Original)
    SHIPSIM_FIELD(alpha)
    SHIPSIM_FIELD(beta)
    SHIPSIM_FIELD(gamma)
SHIPSIM_REFLECT_END(Original)

namespace core {
template <>
struct Reflect<Evolved> {
    using Self = Evolved;
    static constexpr bool kReflected = true;
    static const ::core::TypeInfo& info() {
        static const ::core::FieldInfo fields[] = {
            ::core::makeFieldInfo<decltype(Self::gamma)>("gamma", offsetof(Self, gamma)),
            ::core::makeFieldInfo<decltype(Self::alpha)>("alpha", offsetof(Self, alpha)),
            ::core::makeFieldInfo<decltype(Self::delta)>("delta", offsetof(Self, delta)),
        };
        // Same stable name as Original: this is the same logical type, later.
        static const ::core::TypeInfo type =
            ::core::makeTypeInfo("Original", sizeof(Self), alignof(Self), fields, 3);
        return type;
    }
};
}  // namespace core

namespace {

// The id must be a function of the name and nothing else, or it is not stable
// across builds. Checking it against the hash of the literal is the closed form.
void testStableIdsDependOnlyOnTheName() {
    expectTrue("type id is the hash of the type name",
               core::typeInfo<Hull>().stableId == core::stableHash("Hull"));
    expectTrue("field id is the hash of the field name",
               core::typeInfo<Hull>().findField(core::stableHash("displacement")) != nullptr);
    expectTrue("distinct type names give distinct ids",
               core::typeInfo<Hull>().stableId != core::typeInfo<Vec3f>().stableId);
    // The evolved shape must be recognised as the same logical type.
    expectTrue("a type keeps its id as its fields change",
               core::typeInfo<Original>().stableId == core::typeInfo<Evolved>().stableId);

    // And the hash must actually be constant-evaluated, not computed per call.
    static_assert(core::stableHash("Hull") != core::stableHash("Vec3f"));
    static_assert(core::stableHash("") == 1469598103934665603ull);
}

void testReflectionDescribesLayout() {
    const core::TypeInfo& type = core::typeInfo<Hull>();
    expectEqual("Hull reflects every field", static_cast<long long>(type.fieldCount), 5);
    expectEqual("reflected size matches the type", static_cast<long long>(type.size),
                static_cast<long long>(sizeof(Hull)));

    const core::FieldInfo* offsets = type.findField(core::stableHash("offsets"));
    expectTrue("array field found", offsets != nullptr);
    expectEqual("array extent is recorded", static_cast<long long>(offsets->count), 4);
    expectEqual("array element size is recorded", static_cast<long long>(offsets->elementSize),
                static_cast<long long>(sizeof(double)));
    expectEqual("array offset matches offsetof", static_cast<long long>(offsets->offset),
                static_cast<long long>(offsetof(Hull, offsets)));

    const core::FieldInfo* centre = type.findField(core::stableHash("centre"));
    expectTrue("nested struct field is marked as a struct",
               centre != nullptr && centre->kind == core::TypeKind::Struct);
    expectTrue("nested struct links to its own type info",
               centre != nullptr && centre->nested == &core::typeInfo<Vec3f>());
}

void testRoundTrip() {
    Hull written;
    written.displacement = 8984.25;
    written.centre = {1.5f, -2.5f, 3.25f};
    written.compartments = 16;
    written.afloat = true;
    written.offsets[0] = 0.125;
    written.offsets[1] = -17.0;
    written.offsets[2] = 1e300;
    written.offsets[3] = -0.0;

    ByteWriter writer;
    core::serialise(written, writer);
    expectTrue("something was written", writer.size() > 0);

    Hull read;
    ByteReader reader(writer.span());
    expectTrue("round trip deserialises", core::deserialise(read, reader));
    expectTrue("the whole buffer was consumed", reader.remaining() == 0);

    expectTrue("double survives exactly", read.displacement == written.displacement);
    expectTrue("nested float survives exactly",
               read.centre.x == 1.5f && read.centre.y == -2.5f && read.centre.z == 3.25f);
    expectEqual("int survives", read.compartments, 16);
    expectTrue("bool survives", read.afloat);
    expectTrue("array survives element by element",
               read.offsets[0] == 0.125 && read.offsets[1] == -17.0 && read.offsets[2] == 1e300);
    // Negative zero is the classic case a lazy codec loses.
    expectTrue("negative zero keeps its sign bit",
               read.offsets[3] == 0.0 && std::signbit(read.offsets[3]));
}

// The property the format exists for: data written by one build must load into
// a different build's version of the same type.
void testSchemaEvolution() {
    Original oldRecord;
    oldRecord.alpha = 42;
    oldRecord.beta = 3.5;
    oldRecord.gamma = 0xDEADBEEFCAFEull;

    ByteWriter writer;
    core::serialise(oldRecord, writer);

    // Forward: old data into a newer type that dropped beta, added delta and
    // reordered the rest.
    Evolved newRecord;
    newRecord.delta = 7.5f;  // pre-existing value must be left alone
    ByteReader reader(writer.span());
    expectTrue("old data loads into the evolved type", core::deserialise(newRecord, reader));
    expectEqual("shared field transfers despite reordering",
                static_cast<long long>(newRecord.alpha), 42);
    expectTrue("shared 64-bit field transfers", newRecord.gamma == 0xDEADBEEFCAFEull);
    expectTrue("a field the writer never had keeps its existing value",
               newRecord.delta == 7.5f);
    expectTrue("a field the reader no longer has is skipped cleanly",
               reader.remaining() == 0 && !reader.failed());

    // Backward: new data into the older type.
    ByteWriter forward;
    core::serialise(newRecord, forward);
    Original backRecord;
    backRecord.beta = -1.0;
    ByteReader backReader(forward.span());
    expectTrue("new data loads into the old type", core::deserialise(backRecord, backReader));
    expectEqual("shared field transfers backwards", static_cast<long long>(backRecord.alpha), 42);
    expectTrue("a field the new writer dropped keeps the reader's value",
               backRecord.beta == -1.0);
}

// Reading a foreign type must be refused, not silently reinterpreted.
void testWrongTypeIsRejected() {
    Vec3f source{1, 2, 3};
    ByteWriter writer;
    core::serialise(source, writer);

    Hull target;
    target.compartments = 99;
    ByteReader reader(writer.span());
    expectTrue("deserialising the wrong type fails", !core::deserialise(target, reader));
    expectEqual("the target is left alone on refusal", target.compartments, 99);
}

// Truncated and corrupt input must fail cleanly. ASan proves the reader did not
// step outside the buffer while doing so.
void testMalformedInputFailsCleanly() {
    Hull source;
    source.displacement = 1.0;
    ByteWriter writer;
    core::serialise(source, writer);
    const std::vector<std::byte> full(writer.bytes());

    int survived = 0;
    for (std::size_t length = 0; length < full.size(); ++length) {
        Hull target;
        ByteReader reader(std::span<const std::byte>(full.data(), length));
        // Every truncation must return false; none may read out of bounds.
        if (core::deserialise(target, reader)) ++survived;
    }
    expectEqual("every truncation of a valid buffer is rejected", survived, 0);

    // A corrupt field-count must not be believed as a loop bound.
    std::vector<std::byte> corrupt(full);
    for (int i = 0; i < 4; ++i) corrupt[8 + static_cast<std::size_t>(i)] = std::byte{0xFF};
    Hull target;
    ByteReader corruptReader(corrupt);
    expectTrue("an absurd field count is rejected rather than looped over",
               !core::deserialise(target, corruptReader));

    // A corrupt length prefix must not let the reader run past the buffer.
    std::vector<std::byte> badLength(full);
    if (badLength.size() > 30)
        for (int i = 0; i < 4; ++i) badLength[25 + static_cast<std::size_t>(i)] = std::byte{0x7F};
    Hull lengthTarget;
    ByteReader lengthReader(badLength);
    core::deserialise(lengthTarget, lengthReader);  // may fail; must not over-read
    expectTrue("a corrupt length prefix leaves the reader in bounds",
               lengthReader.cursor() <= badLength.size());
}

// The wire format must be byte-for-byte reproducible, or replays and network
// deltas built on it are not either.
void testEncodingIsDeterministicAndLittleEndian() {
    Vec3f value{1.0f, 0.0f, 0.0f};
    ByteWriter a, b;
    core::serialise(value, a);
    core::serialise(value, b);
    expectTrue("the same value encodes to the same bytes", a.bytes() == b.bytes());

    // 1.0f is 0x3F800000; little-endian that is 00 00 80 3F. Find and confirm it,
    // so the endianness claim is checked rather than assumed.
    const std::vector<std::byte>& bytes = a.bytes();
    bool found = false;
    for (std::size_t i = 0; i + 4 <= bytes.size(); ++i)
        if (bytes[i] == std::byte{0x00} && bytes[i + 1] == std::byte{0x00} &&
            bytes[i + 2] == std::byte{0x80} && bytes[i + 3] == std::byte{0x3F})
            found = true;
    expectTrue("floats are encoded little-endian, not host-endian", found);
}

}  // namespace

void runSerialiseTests() {
    std::printf("\n--- reflection and serialisation ---\n");
    testStableIdsDependOnlyOnTheName();
    testReflectionDescribesLayout();
    testRoundTrip();
    testSchemaEvolution();
    testWrongTypeIsRejected();
    testMalformedInputFailsCleanly();
    testEncodingIsDeterministicAndLittleEndian();
}
