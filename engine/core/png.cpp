// SPDX-License-Identifier: MIT
#include "png.hpp"

#include <cstdio>
#include <cstring>

namespace core {
namespace {

constexpr std::uint8_t kSignature[8] = {0x89, 'P', 'N', 'G', 0x0D, 0x0A, 0x1A, 0x0A};
constexpr std::size_t kBytesPerPixel = 4;  // RGBA8
// Deflate stored blocks carry a 16-bit length, so this is the hard ceiling.
constexpr std::size_t kMaxStoredBlock = 65535;

std::uint32_t crc32(const std::uint8_t* data, std::size_t count, std::uint32_t crc = 0xFFFFFFFFu) {
    static std::uint32_t table[256];
    static bool ready = false;
    if (!ready) {
        for (std::uint32_t i = 0; i < 256; ++i) {
            std::uint32_t c = i;
            for (int k = 0; k < 8; ++k) c = (c & 1u) ? (0xEDB88320u ^ (c >> 1)) : (c >> 1);
            table[i] = c;
        }
        ready = true;
    }
    for (std::size_t i = 0; i < count; ++i)
        crc = table[(crc ^ data[i]) & 0xFFu] ^ (crc >> 8);
    return crc;
}

std::uint32_t adler32(const std::uint8_t* data, std::size_t count) {
    std::uint32_t a = 1, b = 0;
    for (std::size_t i = 0; i < count; ++i) {
        a = (a + data[i]) % 65521u;
        b = (b + a) % 65521u;
    }
    return (b << 16) | a;
}

void pushBigEndian32(std::vector<std::uint8_t>& out, std::uint32_t value) {
    out.push_back(static_cast<std::uint8_t>(value >> 24));
    out.push_back(static_cast<std::uint8_t>(value >> 16));
    out.push_back(static_cast<std::uint8_t>(value >> 8));
    out.push_back(static_cast<std::uint8_t>(value));
}

std::uint32_t readBigEndian32(const std::uint8_t* data) {
    return (static_cast<std::uint32_t>(data[0]) << 24) |
           (static_cast<std::uint32_t>(data[1]) << 16) |
           (static_cast<std::uint32_t>(data[2]) << 8) | static_cast<std::uint32_t>(data[3]);
}

void writeChunk(std::vector<std::uint8_t>& out, const char type[4],
                const std::vector<std::uint8_t>& data) {
    pushBigEndian32(out, static_cast<std::uint32_t>(data.size()));
    const std::size_t crcStart = out.size();
    out.insert(out.end(), type, type + 4);
    out.insert(out.end(), data.begin(), data.end());
    const std::uint32_t crc = crc32(out.data() + crcStart, out.size() - crcStart) ^ 0xFFFFFFFFu;
    pushBigEndian32(out, crc);
}

std::uint8_t paethPredictor(int a, int b, int c) {
    const int p = a + b - c;
    const int pa = std::abs(p - a), pb = std::abs(p - b), pc = std::abs(p - c);
    if (pa <= pb && pa <= pc) return static_cast<std::uint8_t>(a);
    if (pb <= pc) return static_cast<std::uint8_t>(b);
    return static_cast<std::uint8_t>(c);
}

}  // namespace

std::vector<std::byte> encodePng(const Image& image, PngFilter filter) {
    if (!image.valid() || image.width == 0 || image.height == 0) return {};

    const std::size_t stride = std::size_t{image.width} * kBytesPerPixel;

    // Filtered scanlines, each prefixed with its filter byte. This is the buffer
    // the zlib stream carries.
    std::vector<std::uint8_t> raw;
    raw.reserve((stride + 1) * image.height);
    std::vector<std::uint8_t> previous(stride, 0);

    for (std::uint32_t y = 0; y < image.height; ++y) {
        const std::uint8_t* line = image.rgba.data() + std::size_t{y} * stride;
        raw.push_back(static_cast<std::uint8_t>(filter));
        for (std::size_t x = 0; x < stride; ++x) {
            const int current = line[x];
            const int left = x >= kBytesPerPixel ? line[x - kBytesPerPixel] : 0;
            const int up = previous[x];
            const int upLeft = x >= kBytesPerPixel ? previous[x - kBytesPerPixel] : 0;
            int value = current;
            switch (filter) {
                case PngFilter::None: break;
                case PngFilter::Sub: value = current - left; break;
                case PngFilter::Up: value = current - up; break;
                case PngFilter::Average: value = current - ((left + up) >> 1); break;
                case PngFilter::Paeth: value = current - paethPredictor(left, up, upLeft); break;
            }
            raw.push_back(static_cast<std::uint8_t>(value & 0xFF));
        }
        std::memcpy(previous.data(), line, stride);
    }

    // zlib wrapper around stored deflate blocks.
    std::vector<std::uint8_t> zlib;
    zlib.push_back(0x78);  // CMF: deflate, 32K window
    zlib.push_back(0x01);  // FLG chosen so (CMF << 8 | FLG) % 31 == 0
    for (std::size_t offset = 0; offset < raw.size();) {
        const std::size_t count = std::min(kMaxStoredBlock, raw.size() - offset);
        const bool last = offset + count >= raw.size();
        zlib.push_back(last ? 1 : 0);  // BFINAL, BTYPE = 00 (stored)
        zlib.push_back(static_cast<std::uint8_t>(count & 0xFF));
        zlib.push_back(static_cast<std::uint8_t>(count >> 8));
        const std::uint16_t inverse = static_cast<std::uint16_t>(~count);
        zlib.push_back(static_cast<std::uint8_t>(inverse & 0xFF));
        zlib.push_back(static_cast<std::uint8_t>(inverse >> 8));
        zlib.insert(zlib.end(), raw.begin() + static_cast<long>(offset),
                    raw.begin() + static_cast<long>(offset + count));
        offset += count;
    }
    pushBigEndian32(zlib, adler32(raw.data(), raw.size()));

    std::vector<std::uint8_t> out(kSignature, kSignature + 8);

    std::vector<std::uint8_t> ihdr;
    pushBigEndian32(ihdr, image.width);
    pushBigEndian32(ihdr, image.height);
    ihdr.push_back(8);  // bit depth
    ihdr.push_back(6);  // colour type: RGBA
    ihdr.push_back(0);  // compression: deflate
    ihdr.push_back(0);  // filter method: adaptive
    ihdr.push_back(0);  // interlace: none
    writeChunk(out, "IHDR", ihdr);
    writeChunk(out, "IDAT", zlib);
    writeChunk(out, "IEND", {});

    std::vector<std::byte> bytes(out.size());
    std::memcpy(bytes.data(), out.data(), out.size());
    return bytes;
}

bool decodePng(std::span<const std::byte> bytes, Image& out) {
    const auto* data = reinterpret_cast<const std::uint8_t*>(bytes.data());
    const std::size_t size = bytes.size();
    if (size < 8 || std::memcmp(data, kSignature, 8) != 0) return false;

    std::uint32_t width = 0, height = 0;
    bool haveHeader = false, haveEnd = false;
    std::vector<std::uint8_t> compressed;

    std::size_t cursor = 8;
    while (cursor + 12 <= size) {
        const std::uint32_t length = readBigEndian32(data + cursor);
        // Bound the length against the buffer before trusting it for anything.
        if (length > size - cursor - 12) return false;
        const std::uint8_t* type = data + cursor + 4;
        const std::uint8_t* payload = type + 4;

        const std::uint32_t stored = readBigEndian32(payload + length);
        const std::uint32_t actual = crc32(type, length + 4u) ^ 0xFFFFFFFFu;
        if (stored != actual) return false;  // corrupt chunk

        if (std::memcmp(type, "IHDR", 4) == 0) {
            if (length != 13) return false;
            width = readBigEndian32(payload);
            height = readBigEndian32(payload + 4);
            if (width == 0 || height == 0) return false;
            if (payload[8] != 8 || payload[9] != 6) return false;  // only RGBA8
            if (payload[10] != 0 || payload[11] != 0 || payload[12] != 0) return false;
            haveHeader = true;
        } else if (std::memcmp(type, "IDAT", 4) == 0) {
            compressed.insert(compressed.end(), payload, payload + length);
        } else if (std::memcmp(type, "IEND", 4) == 0) {
            haveEnd = true;
        }
        cursor += 12 + length;
    }
    if (!haveHeader || !haveEnd || compressed.size() < 6) return false;

    // zlib header, then stored deflate blocks.
    if ((compressed[0] & 0x0F) != 8) return false;                       // must be deflate
    if (((compressed[0] << 8) | compressed[1]) % 31 != 0) return false;  // header check
    if ((compressed[1] & 0x20) != 0) return false;                       // no preset dictionary

    std::vector<std::uint8_t> raw;
    std::size_t at = 2;
    bool finished = false;
    while (!finished) {
        if (at + 5 > compressed.size()) return false;
        const std::uint8_t header = compressed[at];
        if (((header >> 1) & 0x03) != 0) return false;  // only stored blocks are supported
        finished = (header & 1u) != 0;
        const std::size_t count =
            static_cast<std::size_t>(compressed[at + 1]) | (static_cast<std::size_t>(compressed[at + 2]) << 8);
        const std::size_t inverse =
            static_cast<std::size_t>(compressed[at + 3]) | (static_cast<std::size_t>(compressed[at + 4]) << 8);
        if ((count ^ 0xFFFFu) != inverse) return false;  // LEN/NLEN disagree
        at += 5;
        if (at + count > compressed.size()) return false;
        raw.insert(raw.end(), compressed.begin() + static_cast<long>(at),
                   compressed.begin() + static_cast<long>(at + count));
        at += count;
    }
    if (at + 4 > compressed.size()) return false;
    const std::uint32_t expectedAdler = readBigEndian32(compressed.data() + at);
    if (adler32(raw.data(), raw.size()) != expectedAdler) return false;

    const std::size_t stride = std::size_t{width} * kBytesPerPixel;
    if (raw.size() != (stride + 1) * height) return false;

    out = Image(width, height);
    std::vector<std::uint8_t> previous(stride, 0);
    for (std::uint32_t y = 0; y < height; ++y) {
        const std::uint8_t filter = raw[(stride + 1) * y];
        if (filter > 4) return false;
        const std::uint8_t* source = raw.data() + (stride + 1) * y + 1;
        std::uint8_t* line = out.rgba.data() + std::size_t{y} * stride;

        for (std::size_t x = 0; x < stride; ++x) {
            const int left = x >= kBytesPerPixel ? line[x - kBytesPerPixel] : 0;
            const int up = previous[x];
            const int upLeft = x >= kBytesPerPixel ? previous[x - kBytesPerPixel] : 0;
            int value = source[x];
            switch (filter) {
                case 0: break;
                case 1: value += left; break;
                case 2: value += up; break;
                case 3: value += (left + up) >> 1; break;
                default: value += paethPredictor(left, up, upLeft); break;
            }
            line[x] = static_cast<std::uint8_t>(value & 0xFF);
        }
        std::memcpy(previous.data(), line, stride);
    }
    return true;
}

bool writePng(const std::string& path, const Image& image, PngFilter filter) {
    const std::vector<std::byte> bytes = encodePng(image, filter);
    if (bytes.empty()) return false;
    std::FILE* file = std::fopen(path.c_str(), "wb");
    if (file == nullptr) return false;
    const std::size_t written = std::fwrite(bytes.data(), 1, bytes.size(), file);
    std::fclose(file);
    return written == bytes.size();
}

bool readPng(const std::string& path, Image& out) {
    std::FILE* file = std::fopen(path.c_str(), "rb");
    if (file == nullptr) return false;
    std::fseek(file, 0, SEEK_END);
    const long size = std::ftell(file);
    std::fseek(file, 0, SEEK_SET);
    if (size <= 0) {
        std::fclose(file);
        return false;
    }
    std::vector<std::byte> bytes(static_cast<std::size_t>(size));
    const std::size_t read = std::fread(bytes.data(), 1, bytes.size(), file);
    std::fclose(file);
    if (read != bytes.size()) return false;
    return decodePng(bytes, out);
}

}  // namespace core
