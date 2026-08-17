// SPDX-License-Identifier: MIT
//
// Minimal PNG encoder and decoder.
//
// This exists to be the renderer's verification harness. Rendering is verified
// offscreen -- draw to an image, write it out, read it back, and assert on pixel
// content -- because a windowed check cannot run in CI and "it looked right" is
// not a test.
//
// Written rather than pulled in because the dependency budget for this phase is
// a windowing library and nothing else. PNG needs a zlib stream, but deflate
// permits *stored* (uncompressed) blocks, so a valid file needs only CRC-32,
// Adler-32 and a few headers. The files are larger than a real encoder's; they
// are read by tests, not shipped.
//
// All five PNG scanline filters are implemented in both directions, so the
// decoder can also read images this encoder did not produce.
#pragma once

#include <cstddef>
#include <cstdint>
#include <span>
#include <string>
#include <vector>

namespace core {

struct Image {
    std::uint32_t width = 0;
    std::uint32_t height = 0;
    std::vector<std::uint8_t> rgba;  // width * height * 4, row-major, top-down

    Image() = default;
    Image(std::uint32_t w, std::uint32_t h) : width(w), height(h), rgba(std::size_t{w} * h * 4, 0) {}

    // **A size check must not be the constructor's own arithmetic.** This read
    // `rgba.size() == std::size_t{width} * height * 4` -- character for character
    // the expression the constructor allocates from -- so when that product wraps,
    // `valid()` recomputes the identical wrapped value and agrees with the buffer
    // built from it. `Image(1u << 31, 1u << 31)` wants 2^64 bytes, which is 0 in a
    // `size_t`: `rgba` comes out **empty**, both dimensions stay 2^31, and this
    // returned true. It is the one failure the function is positioned to catch and
    // the one it could not, because a guard that recomputes the bug agrees with it.
    //
    // Checked by division, which cannot wrap. It matters beyond tidiness because
    // `encodePng` uses `valid()` as its only size guard before indexing.
    bool valid() const {
        if (width == 0 || height == 0) return rgba.empty();
        const std::size_t pixels = std::size_t{width} * height;
        if (pixels / width != height) return false;                 // the product wrapped
        if (pixels > std::size_t(-1) / 4) return false;             // and so would the * 4
        return rgba.size() == pixels * 4;
    }

    std::uint8_t* pixel(std::uint32_t x, std::uint32_t y) {
        return rgba.data() + (std::size_t{y} * width + x) * 4;
    }
    const std::uint8_t* pixel(std::uint32_t x, std::uint32_t y) const {
        return rgba.data() + (std::size_t{y} * width + x) * 4;
    }
    void setPixel(std::uint32_t x, std::uint32_t y, std::uint8_t r, std::uint8_t g,
                  std::uint8_t b, std::uint8_t a = 255) {
        std::uint8_t* p = pixel(x, y);
        p[0] = r;
        p[1] = g;
        p[2] = b;
        p[3] = a;
    }
};

enum class PngFilter : std::uint8_t {
    None = 0,
    Sub = 1,
    Up = 2,
    Average = 3,
    Paeth = 4,
};

std::vector<std::byte> encodePng(const Image& image, PngFilter filter = PngFilter::Paeth);
// Returns false on anything malformed: bad signature, bad chunk CRC, truncation,
// an unsupported colour type, or a size that disagrees with the pixel data.
bool decodePng(std::span<const std::byte> bytes, Image& out);

bool writePng(const std::string& path, const Image& image,
              PngFilter filter = PngFilter::Paeth);
bool readPng(const std::string& path, Image& out);

}  // namespace core
