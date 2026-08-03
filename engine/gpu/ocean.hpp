// SPDX-License-Identifier: MIT
//
// A renderable sea surface, driven by the wave field the physics reads.
//
// **The grid is displaced by `sim::WaveField` itself, not by a visual copy of
// it.** `docs/02-simulation.md` and `docs/03-renderer-audio.md` both state this
// and it is the one constraint here that is not negotiable: if the wave the ship
// responds to is not the wave drawn under the bow, the premise collapses. So
// there is no second spectrum, no separate amplitudes, and no re-seeded phases --
// `OceanSurface::build` walks the same `WaveComponent` array
// `WaveField::elevation()` walks, and `tests/test_ocean.cpp` asserts the two
// agree at the vertices to 1e-6 m, which is float32 storage and nothing else.
//
// **Displacement is on the CPU, for now, and it is the whole cost.** Measured on
// the development machine (one core, -O2), a 576-component sea (Hs 3 m, Tp 9 s,
// N = 48, M = 12) over a 129 x 129 grid:
//
//     obvious loop, sincos per vertex per component   190 ms    19.8 ns/vertex-component
//     rotation recurrence along each grid row          19 ms     2.0 ns/vertex-component
//     four recurrences interleaved (not shipped)        8.7 ms   0.9 ns/vertex-component
//
// The recurrence is what is implemented. Along a grid row the phase advances by a
// constant, so `(cos psi, sin psi)` can be stepped by a fixed rotation instead of
// re-evaluated: one `sincos` per row per component instead of one per vertex, and
// the sine that falls out is exactly what the slope needs. It is an algebraic
// rearrangement of the same sum, not an approximation -- measured against direct
// evaluation the worst disagreement over a 257 x 257 grid is 2e-14 m.
//
// Even so this does not scale to a full-resolution sea (see the note on
// resolution below), which is why `docs/03-renderer-audio.md` plans a Tessendorf
// FFT cascade on the compute queue. Evaluating the spectrum in a vertex shader
// would be the same sum with the same components -- correct, but it re-does the
// per-vertex sincos the recurrence just removed, so it is a step sideways rather
// than forward. The FFT is the step forward, and it still reads this spectrum.
//
// **Resolution is a correctness constraint, not a quality knob.**
// `engine/core/geometry.hpp` records the same hazard on the hull side: a panel
// spanning several wavelengths reports whatever its three corners happen to say,
// and the error is a systematic phase-dependent bias rather than noise. Measured
// here on Hs 3 m, Tp 9 s, N = 48 (dominant wavelength 126 m, *shortest component
// 11.4 m*), worst cell-centre error against the analytic field:
//
//     cells per dominant wavelength    16      32      64     128
//     grid spacing                   7.90 m  3.95 m  1.98 m  0.99 m
//     worst error                    24%     8.2%    2.2%    0.56%   of Hs
//
// Resolving the dominant wavelength is *not* enough. The criterion is the
// shortest component the discretisation carries, and eight cells across it is
// about where the surface stops inventing height -- which for that sea state is a
// 1.4 m grid, 80 000 vertices over a 400 m patch, and 92 ms per frame on one
// core. `shortestWavelength()` exists so callers can compute that rather than
// guess it.
//
// Pascal target, so: no mesh shaders and no hardware ray tracing. This is a
// plain vertex/fragment pipeline over an indexed triangle grid.
#pragma once

#include "../core/png.hpp"
#include "../sim/waves.hpp"
#include "device.hpp"

#include <cstdint>
#include <string>
#include <vector>

struct VkRenderPass_T;
struct VkFramebuffer_T;
struct VkPipeline_T;
struct VkPipelineLayout_T;

namespace gpu {

// World position and the analytic surface normal. The normal comes from the
// spectrum's own slope, not from a finite difference of the mesh, so it stays
// right even where the grid is too coarse for the geometry -- and the lighting
// then responds to the actual waves rather than to the tessellation.
struct OceanVertex {
    float position[3];
    float normal[3];
};

// A square patch of sea, in world metres. Vertices are (resolution + 1)^2.
struct OceanGrid {
    double centreX = 0.0;
    double centreY = 0.0;
    double halfExtent = 100.0;  // the patch spans 2 * halfExtent in x and in y
    int resolution = 128;       // cells per side

    double cellSize() const { return 2.0 * halfExtent / (resolution > 0 ? resolution : 1); }
};

// Cells per side that puts `cellsPerWavelength` cells across `wavelength`.
// Rounded up, never below one.
int oceanResolutionFor(double halfExtent, double wavelength, double cellsPerWavelength);

// Wavelength at the peak of the *discretised* spectrum: the frequency bin with
// the highest energy per unit frequency, which with equal-energy bins is the
// narrowest one. Zero for a still field.
double dominantWavelength(const sim::WaveField& field);
// Wavelength of the shortest component the field carries. This, and not the
// dominant wavelength, is what a grid has to resolve to reproduce the field.
double shortestWavelength(const sim::WaveField& field);

// The displaced grid. Separated from the renderer so the geometry can be
// asserted against closed forms with no GPU in the picture.
class OceanSurface {
public:
    // Evaluates `field` over `grid` at `time`. Reuses its own storage across
    // calls; the index buffer is rebuilt only when the resolution changes.
    void build(const sim::WaveField& field, const OceanGrid& grid, double time);

    const OceanGrid& grid() const { return grid_; }
    double time() const { return time_; }
    const std::vector<OceanVertex>& vertices() const { return vertices_; }
    const std::vector<std::uint32_t>& indices() const { return indices_; }
    // Wall-clock seconds the last build() took. Reported rather than asserted:
    // it is a measurement, and docs/03-renderer-audio.md carries the numbers.
    double buildSeconds() const { return buildSeconds_; }

    // Elevation of the *mesh* at a world point -- the piecewise-linear surface
    // the rasteriser interpolates, which is not the analytic elevation and
    // differs from it by exactly the amount the grid is too coarse by. False
    // outside the patch.
    bool sampleElevation(double x, double y, double& elevation) const;

    // Index of the vertex at grid coordinates (i, j), i and j in [0, resolution].
    std::size_t vertexIndex(int i, int j) const;

private:
    OceanGrid grid_;
    double time_ = 0.0;
    double buildSeconds_ = 0.0;
    int indexResolution_ = -1;
    std::vector<OceanVertex> vertices_;
    std::vector<std::uint32_t> indices_;
    // Per-row accumulators, so a whole row's worth of components stays in cache.
    std::vector<double> rowElevation_, rowSlopeX_, rowSlopeY_;
};

enum class OceanShading {
    // Lambert against a directional sun plus a hemispheric sky term.
    // **Deliberately view-independent** -- no specular, no Fresnel -- so a flat
    // sea has exactly one colour over its whole surface and that colour is a
    // closed form a test can state:
    //
    //     colour = waterColour * (ambient * (0.5 + 0.5 * n.z) + sunStrength * max(dot(n, sun), 0))
    //
    // A specular lobe is the obvious next thing and it costs that assertion, so
    // it arrives with the reflection work rather than before it.
    Shaded,
    // World elevation as a 16-bit code across [elevationMin, +elevationSpan],
    // high byte in red, low byte in green, blue set to 255 as a surface tag.
    // This is how a test reads geometry back out of a rendered frame and holds
    // it against `WaveField::elevation()`; it is a measurement channel, not a
    // debug colour ramp, and it is exact to span / 65535.
    Elevation,
};

struct OceanView {
    OceanShading shading = OceanShading::Shaded;
    // Unit vector *toward* the sun, world frame. Normalised on submission.
    float sunDirection[3] = {0.48f, 0.36f, 0.80f};
    float sunStrength = 0.55f;
    float waterColour[3] = {0.10f, 0.32f, 0.48f};
    float ambient = 0.45f;
    // Range of the Elevation encoding. Elevations outside it saturate.
    float elevationMin = -8.0f;
    float elevationSpan = 16.0f;
};

// Inverse of the Elevation encoding. False when the pixel is not sea.
bool decodeOceanElevation(const OceanView& view, const std::uint8_t* pixel, double& elevation);

// Offscreen colour + depth target and the grid pipeline. A sibling of
// OffscreenRenderer rather than a user of it: the sea needs a normal on the
// vertex and lighting in the fragment shader, and MeshVertex has neither.
class OceanRenderer {
public:
    ~OceanRenderer();

    bool create(Device& device, std::uint32_t width, std::uint32_t height,
                const std::string& shaderDirectory, std::string& error);
    void destroy();
    bool valid() const { return pipeline_ != nullptr; }

    std::uint32_t width() const { return width_; }
    std::uint32_t height() const { return height_; }

    // Clears to `clearColour`, draws `surface` under `mvp` with depth testing,
    // and reads the colour attachment back into `out`. An unbuilt surface is
    // legal and yields the clear colour.
    bool render(const float mvp[16], const OceanSurface& surface, const OceanView& view,
                const float clearColour[4], core::Image& out);

private:
    bool ensureGeometryCapacity(std::size_t vertexBytes, std::size_t indexBytes);

    Device*             device_ = nullptr;
    std::uint32_t       width_ = 0;
    std::uint32_t       height_ = 0;

    Device::Image2D     colour_;
    Device::Image2D     depth_;
    VkRenderPass_T*     renderPass_ = nullptr;
    VkFramebuffer_T*    framebuffer_ = nullptr;
    VkPipelineLayout_T* pipelineLayout_ = nullptr;
    VkPipeline_T*       pipeline_ = nullptr;

    Device::Buffer      vertexBuffer_;
    Device::Buffer      indexBuffer_;
    Device::Buffer      readback_;
};

}  // namespace gpu
