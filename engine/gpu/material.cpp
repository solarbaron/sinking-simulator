// SPDX-License-Identifier: MIT
#include "material.hpp"

#include <algorithm>
#include <charconv>
#include <fstream>
#include <sstream>
#include <system_error>

namespace gpu {
namespace {

// Whitespace-separated tokens of one line, with everything from '#' dropped.
void tokenise(std::string_view line, std::vector<std::string_view>& out) {
    out.clear();
    const std::size_t comment = line.find('#');
    if (comment != std::string_view::npos) line = line.substr(0, comment);
    std::size_t i = 0;
    while (i < line.size()) {
        while (i < line.size() && (line[i] == ' ' || line[i] == '\t' || line[i] == '\r')) ++i;
        if (i >= line.size()) break;
        const std::size_t start = i;
        while (i < line.size() && line[i] != ' ' && line[i] != '\t' && line[i] != '\r') ++i;
        out.push_back(line.substr(start, i - start));
    }
}

// std::from_chars rather than strtod: the decimal point must not depend on the
// host's locale, or the same ship file means different things on two machines.
// The whole token has to be consumed, so "0.4x" is rejected rather than read as
// 0.4 -- a loader that accepts a typo is a loader that fails open.
bool parseNumber(std::string_view token, double& out) {
    const char* begin = token.data();
    const char* end = begin + token.size();
    const std::from_chars_result result = std::from_chars(begin, end, out);
    return result.ec == std::errc{} && result.ptr == end;
}

std::string at(const std::string& origin, int line, const std::string& what) {
    std::ostringstream message;
    message << origin << ":" << line << ": " << what;
    return message.str();
}

GpuMaterial pack(const Material& material) {
    GpuMaterial packed{};
    for (int c = 0; c < 3; ++c) packed.baseColour[c] = static_cast<float>(material.baseColour[c]);
    packed.baseColour[3] = static_cast<float>(material.opacity);
    packed.params[0] = static_cast<float>(material.roughness);
    packed.params[1] = static_cast<float>(material.metalness);
    packed.params[2] = 0.0f;
    packed.params[3] = 0.0f;
    return packed;
}

}  // namespace

void MaterialLibrary::clear() {
    materials_.clear();
    packed_.clear();
    ++revision_;
}

int MaterialLibrary::find(std::string_view name) const {
    for (std::size_t i = 0; i < materials_.size(); ++i)
        if (materials_[i].name == name) return static_cast<int>(i);
    return -1;
}

bool MaterialLibrary::load(const std::string& path, std::string& error) {
    std::ifstream file(path, std::ios::binary);
    if (!file) {
        error = path + ": cannot open";
        return false;
    }
    std::ostringstream contents;
    contents << file.rdbuf();
    if (file.bad()) {
        error = path + ": read failed";
        return false;
    }
    return parse(contents.str(), path, error);
}

bool MaterialLibrary::parse(std::string_view text, const std::string& origin, std::string& error) {
    // Everything lands here first. `materials_` is only touched once the whole
    // file has parsed, so a failure half way through cannot leave a material with
    // a name and default everything else -- which would render, and look almost
    // right, which is the failure mode worth engineering against.
    struct Pending {
        Material material;
        bool hasBaseColour = false;
        bool hasRoughness = false;
        bool hasMetalness = false;
        int declaredAt = 0;
    };
    std::vector<Pending> pending;
    std::vector<std::string_view> tokens;

    const auto complete = [&](const Pending& block, std::string& why) {
        if (!block.hasBaseColour) { why = "material '" + block.material.name +
                                         "' has no base_colour"; return false; }
        if (!block.hasRoughness)  { why = "material '" + block.material.name +
                                         "' has no roughness";  return false; }
        if (!block.hasMetalness)  { why = "material '" + block.material.name +
                                         "' has no metalness";  return false; }
        return true;
    };

    int lineNumber = 0;
    std::size_t cursor = 0;
    while (cursor <= text.size()) {
        const std::size_t newline = text.find('\n', cursor);
        const std::string_view line =
            text.substr(cursor, newline == std::string_view::npos ? std::string_view::npos
                                                                  : newline - cursor);
        cursor = newline == std::string_view::npos ? text.size() + 1 : newline + 1;
        ++lineNumber;

        tokenise(line, tokens);
        if (tokens.empty()) continue;

        if (tokens[0] == "material") {
            if (tokens.size() != 2) {
                error = at(origin, lineNumber, "'material' takes exactly one name");
                return false;
            }
            if (!pending.empty()) {
                std::string why;
                if (!complete(pending.back(), why)) {
                    error = at(origin, pending.back().declaredAt, why);
                    return false;
                }
            }
            const std::string name(tokens[1]);
            for (const Pending& earlier : pending)
                if (earlier.material.name == name) {
                    error = at(origin, lineNumber, "material '" + name +
                                                       "' is defined twice in this file");
                    return false;
                }
            Pending block;
            block.material.name = name;
            block.declaredAt = lineNumber;
            pending.push_back(block);
            continue;
        }

        if (pending.empty()) {
            error = at(origin, lineNumber,
                       "'" + std::string(tokens[0]) + "' outside any material block");
            return false;
        }
        Pending& block = pending.back();

        const auto scalar = [&](double lo, double hi, double& out, bool& flag) {
            if (tokens.size() != 2) {
                error = at(origin, lineNumber,
                           "'" + std::string(tokens[0]) + "' takes one number");
                return false;
            }
            double value = 0.0;
            if (!parseNumber(tokens[1], value)) {
                error = at(origin, lineNumber, "'" + std::string(tokens[1]) + "' is not a number");
                return false;
            }
            if (!(value >= lo) || !(value <= hi)) {
                std::ostringstream why;
                why << "'" << tokens[0] << "' must be within [" << lo << ", " << hi << "], got "
                    << value;
                error = at(origin, lineNumber, why.str());
                return false;
            }
            out = value;
            flag = true;
            return true;
        };

        if (tokens[0] == "base_colour") {
            if (tokens.size() != 4) {
                error = at(origin, lineNumber, "'base_colour' takes three numbers");
                return false;
            }
            for (int c = 0; c < 3; ++c) {
                double value = 0.0;
                if (!parseNumber(tokens[static_cast<std::size_t>(c) + 1], value)) {
                    error = at(origin, lineNumber,
                               "'" + std::string(tokens[static_cast<std::size_t>(c) + 1]) +
                                   "' is not a number");
                    return false;
                }
                if (!(value >= 0.0) || !(value <= 1.0)) {
                    error = at(origin, lineNumber, "base_colour components are reflectances and"
                                                   " must be within [0, 1]");
                    return false;
                }
                block.material.baseColour[c] = value;
            }
            block.hasBaseColour = true;
        } else if (tokens[0] == "roughness") {
            if (!scalar(kMinRoughness, 1.0, block.material.roughness, block.hasRoughness))
                return false;
        } else if (tokens[0] == "metalness") {
            if (!scalar(0.0, 1.0, block.material.metalness, block.hasMetalness)) return false;
        } else if (tokens[0] == "opacity") {
            bool ignored = false;
            if (!scalar(0.0, 1.0, block.material.opacity, ignored)) return false;
        } else {
            error = at(origin, lineNumber, "unknown key '" + std::string(tokens[0]) + "'");
            return false;
        }
    }

    if (!pending.empty()) {
        std::string why;
        if (!complete(pending.back(), why)) {
            error = at(origin, pending.back().declaredAt, why);
            return false;
        }
    }

    // Committed. A name already present is replaced in place so indices resolved
    // before the mod loaded stay pointing at the same surface.
    for (const Pending& block : pending) {
        const int existing = find(block.material.name);
        if (existing >= 0) {
            materials_[static_cast<std::size_t>(existing)] = block.material;
            packed_[static_cast<std::size_t>(existing)] = pack(block.material);
        } else {
            materials_.push_back(block.material);
            packed_.push_back(pack(block.material));
        }
    }
    ++revision_;
    return true;
}

}  // namespace gpu
