// SPDX-License-Identifier: MIT
//
// Surface materials, and the analytic BRDF the hull is shaded with.
//
// **The material set is data, not code.** `docs/05-data-modding-validation.md`
// commits to a ship being a directory rather than a compiled artefact, and a
// paint scheme is part of a ship. So there is no enumeration of materials in the
// shader and none in this header either: the shipped set lives in
// `engine/gpu/materials/marine.materials`, `MaterialLibrary::load` merges any
// number of such files, and the fragment shader indexes a storage buffer built
// from whatever was loaded. A mod adds a material by adding a block to a text
// file and naming it; nothing recompiles and no shader changes.
//
// **The load path fails closed.** CLAUDE.md records a loader that failed open and
// left a half-built world, so `load` parses into a scratch copy and commits only
// on success -- an unknown key, a number out of range, a duplicate name inside one
// file, or a block that stops before it is complete all leave the library exactly
// as it was, with a line-numbered reason. That is checked against every truncation
// of a valid file in `tests/test_hull_render.cpp`.
//
// **The BRDF.** Metallic-roughness, which is the parameterisation glTF, Filament
// and every DCC tool already speak, so an imported material means something here.
// Per channel c, for a surface normal n, unit vector l toward the light, unit
// vector v toward the eye, and h = normalize(l + v):
//
//     L_c = base_c * sky_c * (0.5 + 0.5 n.z)                          sky
//         + (base_c (1 - metal) / pi + D * Vis * F_c) * max(n.l, 0) * sun_c
//
//     alpha = roughness^2,   a2 = alpha^2
//     D     = a2 / (pi * ((n.h)^2 (a2 - 1) + 1)^2)        GGX / Trowbridge-Reitz
//     Vis   = 0.5 / ( n.l sqrt((n.v)^2 (1 - a2) + a2)
//                   + n.v sqrt((n.l)^2 (1 - a2) + a2) )   Smith, height-correlated,
//                                                         with 1 / (4 n.l n.v) folded in
//     F0_c  = mix(0.04, base_c, metal)
//     F_c   = F0_c + (1 - F0_c) (1 - v.h)^5               Schlick
//
// and the stored pixel is `round(255 * clamp(exposure * L_c, 0, 1))`.
//
// Three deliberate choices, each of which buys a closed form:
//
//   * **The diffuse lobe is exactly Lambert** -- no (1 - F) factor, which glTF's
//     reference BRDF carries and Filament drops. Keeping it would couple the
//     diffuse term to the view direction and cost the exact cosine law that
//     `tests/test_hull_render.cpp` sweeps over nineteen angles. It is worth a few
//     per cent at grazing incidence, and the sweep is worth more.
//   * **No tone mapping and no sRGB curve.** The pixel is a clamped linear
//     radiance, so every assertion is an equality rather than an inversion of a
//     curve. The transfer function belongs with the post chain, and arrives with
//     it -- the same trade `docs/03-renderer-audio.md` records for the sea's
//     specular.
//   * **The hemispheric sky term is not scaled by (1 - metal).** A metal has no
//     diffuse lobe, but it does mirror the sky, and its mirror reflectance *is*
//     its base colour; base * sky is a cheap stand-in for that until image-based
//     lighting exists.
//
// Roughness below 0.03 is rejected rather than clamped. At alpha = 9e-4 the GGX
// peak is already 4e5 and the lobe is a fraction of a pixel wide, so a smaller
// value is not a smoother surface, it is a sampling problem -- and silently
// clamping a mod's number is how data stops meaning what it says.
#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

namespace gpu {

// One surface. Doubles because this is authoring data; it is packed to float on
// the way to the GPU and nowhere else.
struct Material {
    std::string name;
    double baseColour[3] = {0.5, 0.5, 0.5};  // linear reflectance, each in [0, 1]
    double roughness = 0.5;                  // perceptual; alpha = roughness^2
    double metalness = 0.0;                  // 0 dielectric, 1 conductor
    // Carried because it is a property of the material and the data must not lie
    // about the ship. **The current pass ignores it** -- there is no transparency
    // sort and no refraction, so glass renders as the smooth dark dielectric it is
    // from outside. See docs/03-renderer-audio.md.
    double opacity = 1.0;
};

// The std430 row the fragment shader indexes. Mirrors `Material` in
// shaders/hull.frag; both are two vec4s, which std430 packs tightly.
struct GpuMaterial {
    float baseColour[4];  // rgb linear reflectance, w opacity
    float params[4];      // x roughness, y metalness, zw reserved
};
static_assert(sizeof(GpuMaterial) == 32, "the material row must match shaders/hull.frag");

// Smallest roughness a material may declare. See the header comment.
inline constexpr double kMinRoughness = 0.03;

class MaterialLibrary {
public:
    // Parse `path` and merge it in. A name already present is **replaced**, which
    // is how a mod overrides a shipped material; a name repeated inside one file
    // is an error, because that is a typo rather than an intent.
    //
    // On any failure the library is left exactly as it was and `error` carries
    // the origin, the line number and the reason.
    bool load(const std::string& path, std::string& error);
    // Same, from memory. `origin` only labels the error messages.
    bool parse(std::string_view text, const std::string& origin, std::string& error);

    // Index of `name`, or -1. Callers resolve names once and pass indices to the
    // GPU; the shader never sees a string.
    int find(std::string_view name) const;

    std::size_t size() const { return materials_.size(); }
    bool empty() const { return materials_.empty(); }
    const std::vector<Material>& materials() const { return materials_; }
    const Material& operator[](std::size_t index) const { return materials_[index]; }

    // The packed table, kept in step with `materials()` so a renderer can upload
    // it directly. `revision()` changes whenever either does, which is how the
    // renderer knows to re-upload without comparing contents.
    const std::vector<GpuMaterial>& packed() const { return packed_; }
    std::uint64_t revision() const { return revision_; }

    void clear();

private:
    std::vector<Material> materials_;
    std::vector<GpuMaterial> packed_;
    std::uint64_t revision_ = 0;
};

}  // namespace gpu
