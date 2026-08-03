// SPDX-License-Identifier: MIT
//
// Validation of the PNG codec that the renderer's offscreen harness depends on.
//
// If this is quietly wrong, every rendering test built on it is worthless: a
// decoder that drops the last scanline, or an encoder whose CRCs are bad, would
// still produce a file and still "pass" a check that only asserts the file
// exists. So the round trip is pixel-exact, the file structure is parsed and
// checked by hand, and corruption is required to be detected rather than
// tolerated.
#include "engine/core/png.hpp"
#include "harness.hpp"

#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

using core::Image;
using core::PngFilter;
using testing::expectEqual;
using testing::expectTrue;

namespace {

// A pattern with structure in both axes and in every channel, so a transposed,
// shifted or channel-swapped decode cannot accidentally match.
Image makePattern(std::uint32_t width, std::uint32_t height) {
    Image image(width, height);
    for (std::uint32_t y = 0; y < height; ++y)
        for (std::uint32_t x = 0; x < width; ++x)
            image.setPixel(x, y, static_cast<std::uint8_t>(x * 7 + y),
                           static_cast<std::uint8_t>(y * 13 + 5),
                           static_cast<std::uint8_t>((x * x + y * 3) & 0xFF),
                           static_cast<std::uint8_t>(200 + ((x + y) & 0x3F)));
    return image;
}

bool identical(const Image& a, const Image& b) {
    return a.width == b.width && a.height == b.height && a.rgba == b.rgba;
}

// Every filter must round-trip exactly. Filters are where an off-by-one in the
// left/up/up-left neighbour selection hides: the image still decodes, just wrong.
void testEveryFilterRoundTripsExactly() {
    const Image source = makePattern(61, 37);  // deliberately not a round number
    for (PngFilter filter : {PngFilter::None, PngFilter::Sub, PngFilter::Up,
                             PngFilter::Average, PngFilter::Paeth}) {
        const auto encoded = core::encodePng(source, filter);
        expectTrue("filter " + std::to_string(static_cast<int>(filter)) + " encodes",
                   !encoded.empty());

        Image decoded;
        expectTrue("filter " + std::to_string(static_cast<int>(filter)) + " decodes",
                   core::decodePng(encoded, decoded));
        expectTrue("filter " + std::to_string(static_cast<int>(filter)) + " is pixel-exact",
                   identical(source, decoded));
    }
}

// The bytes must actually be a PNG, not merely something this decoder accepts.
// Parsed by hand here rather than through the decoder, so an encoder and decoder
// that agree on a private format still fail.
void testFileStructureIsAValidPng() {
    const Image source = makePattern(16, 9);
    const auto encoded = core::encodePng(source, PngFilter::Paeth);
    const auto* bytes = reinterpret_cast<const std::uint8_t*>(encoded.data());

    const std::uint8_t signature[8] = {0x89, 'P', 'N', 'G', 0x0D, 0x0A, 0x1A, 0x0A};
    expectTrue("the file starts with the PNG signature",
               encoded.size() > 8 && std::memcmp(bytes, signature, 8) == 0);

    // First chunk must be IHDR, and must describe the image we handed in.
    const std::uint32_t ihdrLength =
        (bytes[8] << 24) | (bytes[9] << 16) | (bytes[10] << 8) | bytes[11];
    expectEqual("IHDR is 13 bytes", static_cast<long long>(ihdrLength), 13);
    expectTrue("the first chunk is IHDR", std::memcmp(bytes + 12, "IHDR", 4) == 0);

    const std::uint32_t width =
        (bytes[16] << 24) | (bytes[17] << 16) | (bytes[18] << 8) | bytes[19];
    const std::uint32_t height =
        (bytes[20] << 24) | (bytes[21] << 16) | (bytes[22] << 8) | bytes[23];
    expectEqual("IHDR width matches", static_cast<long long>(width), 16);
    expectEqual("IHDR height matches", static_cast<long long>(height), 9);
    expectEqual("bit depth is 8", static_cast<long long>(bytes[24]), 8);
    expectEqual("colour type is RGBA", static_cast<long long>(bytes[25]), 6);

    // And the file must end with a well-formed IEND.
    expectTrue("the file ends with IEND",
               encoded.size() >= 12 &&
                   std::memcmp(bytes + encoded.size() - 8, "IEND", 4) == 0);
}

// Corruption must be caught. A decoder that ignores CRCs will happily return a
// scrambled image, which is worse than failing.
void testCorruptionIsDetected() {
    const Image source = makePattern(24, 24);
    const auto encoded = core::encodePng(source, PngFilter::Paeth);

    int accepted = 0;
    // Flip a bit in each of a spread of positions past the signature.
    for (std::size_t i = 8; i < encoded.size(); i += 13) {
        std::vector<std::byte> damaged(encoded);
        damaged[i] = static_cast<std::byte>(std::to_integer<std::uint8_t>(damaged[i]) ^ 0x40u);
        Image decoded;
        if (core::decodePng(damaged, decoded)) ++accepted;
    }
    expectEqual("no single-bit corruption decodes successfully", accepted, 0);

    // Truncation likewise, and without reading past the buffer -- ASan checks
    // the second half of that claim.
    int truncationsAccepted = 0;
    for (std::size_t length = 0; length < encoded.size(); length += 5) {
        Image decoded;
        if (core::decodePng(std::span<const std::byte>(encoded.data(), length), decoded))
            ++truncationsAccepted;
    }
    expectEqual("no truncation decodes successfully", truncationsAccepted, 0);

    Image decoded;
    const std::byte garbage[32]{};
    expectTrue("garbage is rejected", !core::decodePng(garbage, decoded));
}

// Stored deflate blocks carry a 16-bit length, so an image whose scanline data
// exceeds 65535 bytes must span several blocks. Getting the block chaining wrong
// truncates the image at exactly 64 KB, which a small test would never notice.
void testImageSpanningMultipleDeflateBlocks() {
    // 200 x 200 RGBA with filter bytes is ~160 KB: three blocks.
    const Image source = makePattern(200, 200);
    expectTrue("the test image really does exceed one stored block",
               source.rgba.size() + source.height > 65535);

    const auto encoded = core::encodePng(source, PngFilter::Up);
    Image decoded;
    expectTrue("a multi-block image decodes", core::decodePng(encoded, decoded));
    expectTrue("a multi-block image is pixel-exact", identical(source, decoded));

    // Check the far corner explicitly: a block-chaining bug leaves it zeroed
    // while the first 64 KB looks perfect.
    const std::uint8_t* corner = decoded.pixel(199, 199);
    const std::uint8_t* expected = source.pixel(199, 199);
    expectTrue("the last pixel survived",
               std::memcmp(corner, expected, 4) == 0);
}

void testFileRoundTrip() {
    const Image source = makePattern(33, 17);
    const std::string path = testing::scratchDir() + "png_roundtrip.png";
    expectTrue("the image writes to disk", core::writePng(path, source));

    Image loaded;
    expectTrue("the image reads back from disk", core::readPng(path, loaded));
    expectTrue("the file round trip is pixel-exact", identical(source, loaded));

    Image missing;
    expectTrue("a missing file is reported rather than crashing",
               !core::readPng(path + ".does-not-exist", missing));
}

void testDegenerateImages() {
    Image empty;
    expectTrue("an empty image encodes to nothing", core::encodePng(empty).empty());

    const Image single(1, 1);
    const auto encoded = core::encodePng(single, PngFilter::Sub);
    Image decoded;
    expectTrue("a 1x1 image round trips",
               core::decodePng(encoded, decoded) && identical(single, decoded));

    // A single row and a single column each exercise one neighbour being absent.
    const Image row = makePattern(64, 1);
    const Image column = makePattern(1, 64);
    Image decodedRow, decodedColumn;
    expectTrue("a single-row image round trips",
               core::decodePng(core::encodePng(row, PngFilter::Paeth), decodedRow) &&
                   identical(row, decodedRow));
    expectTrue("a single-column image round trips",
               core::decodePng(core::encodePng(column, PngFilter::Paeth), decodedColumn) &&
                   identical(column, decodedColumn));
}

}  // namespace

void runPngTests() {
    std::printf("\n--- png codec ---\n");
    testEveryFilterRoundTripsExactly();
    testFileStructureIsAValidPng();
    testCorruptionIsDetected();
    testImageSpanningMultipleDeflateBlocks();
    testFileRoundTrip();
    testDegenerateImages();
}
