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

// --- The cascade --------------------------------------------------------------

// Concentric square rings, each level twice the cell size and twice the extent of
// the one inside it.
//
// **A uniform patch cannot reach the horizon.** The resolution criterion above is
// about the shortest component the *field* carries, and a uniform grid applies it
// to the water twenty kilometres away as readily as to the water under the bow --
// where it is nonsense, because an 11 m wave at 20 km is a thousandth of a pixel.
// Cost then goes as the square of the reach and the sea stops inside the frame
// instead.
//
// Halving the resolution every time the extent doubles inverts that: **reach is
// exponential in the number of levels and cost is linear in it.** Level 0 is a
// full grid of `resolution` cells; level L is a square annulus of the same cell
// count minus the quarter its predecessor already covers, so each ring costs
// three quarters of level 0 and covers four times the area. Ten levels reach 512
// times as far as one.
//
// Two things have to be true for that to be a sea rather than a pile of patches:
//
//   * **The levels must not crack.** A ring's inner boundary vertices *are* the
//     previous level's outer boundary vertices -- the same entries in the same
//     array, not a second evaluation at the same place -- and every coarse cell
//     that meets the seam is split into three triangles so it uses the fine
//     level's midpoint too. The seam is then a set of shared edges with shared
//     endpoints: no T-junction, no step, and no dependence on two evaluations
//     agreeing to the last bit. `tests/test_ocean.cpp` asserts it combinatorially,
//     as an edge-manifold check over the whole cascade.
//   * **A level must drop what it cannot resolve.** Carrying a component the grid
//     samples below its Nyquist does not render the wave, it renders an artefact
//     with the wave's amplitude. `minimumCellsPerWavelength` is where that line
//     is drawn, and the far rings end up carrying nothing at all, which is what
//     makes them nearly free.
//
// The same `sim::WaveField` drives every level. Dropping a component is not a
// second spectrum; it is this one, band-limited to what the level can carry.
struct OceanCascade {
    double centreX = 0.0;
    double centreY = 0.0;
    // Half-extent of level 0, the finest patch. The cascade reaches
    // `innerHalfExtent * 2^(levels - 1)`.
    double innerHalfExtent = 128.0;
    // Cells per side, the same at every level. Rounded up to a multiple of four:
    // the ring's hole is a quarter of its extent and its boundary vertices have
    // to fall on even lattice indices, and both need that.
    int resolution = 128;
    int levels = 8;
    // A level carries a component only when its wavelength spans at least this
    // many of the level's cells.
    //
    // Two is the floor rather than the answer. The mesh's own interpolation error
    // at a cell centre is `a (1 - cos(pi / n))` for n cells per wavelength -- the
    // identity `tests/test_ocean.cpp` asserts -- which at n = 2 is exactly `a`:
    // the surface carries none of the wave and all of its amplitude as artefact,
    // so keeping the component is precisely as wrong as dropping it, and worse to
    // look at because the artefact has structure. Four leaves 29% and is where a
    // dropped component starts costing more than it returns.
    double minimumCellsPerWavelength = 4.0;
    // Off leaves the T-junctions at the level boundaries in place, which is the
    // defect the stitching exists to prevent. **It is here as a negative control**
    // -- a crack test that has never seen a crack proves nothing -- and nothing
    // that renders should ever set it.
    bool stitchSeams = true;

    // The resolution actually used: `resolution` rounded up to a multiple of four
    // and clamped to something a machine can allocate.
    int cellsPerSide() const;
    double cellSize(int level) const;
    double halfExtent(int level) const;
    double reach() const { return halfExtent(levels - 1); }
};

// Levels needed for a cascade of `innerHalfExtent` to reach at least `reach`.
int oceanCascadeLevelsFor(double innerHalfExtent, double reach);

// How far a flat sea has to reach for its edge to be indistinguishable from the
// horizon, given the camera rather than a taste in numbers.
//
// A flat sea's horizon *is* the eye's own horizontal plane, so a surface point at
// horizontal distance D sits `atan(eyeHeight / D)` below it and a patch that stops
// at D leaves exactly that much sky between its edge and the water that should
// still be there. One pixel at the centre of the frame subtends
// `2 tan(fov / 2) / pixelHeight`, so the edge disappears into the horizon once
//
//     D >= eyeHeight * pixelHeight / (2 tan(verticalFov / 2)).
//
// Height and field of view, not a distance: doubling the eye height doubles the
// reach required, and so does doubling the vertical resolution. Nothing about the
// sea state enters it -- a wave crest is below the eye whatever it is doing, so it
// is the plane and not the waves that has to reach.
double oceanHorizonReach(double eyeHeight, double verticalFov, int pixelHeight);

// The built cascade: one vertex array and one index array for the whole thing,
// because the levels share vertices along their seams and a per-level buffer
// could not.
class OceanCascadeSurface {
public:
    struct Level {
        double cellSize = 0;
        double halfExtent = 0;
        std::size_t components = 0;  // of the field's, the ones this level resolves
        std::size_t vertices = 0;    // this level's own; seam vertices belong inward
        std::size_t triangles = 0;
    };

    // Evaluates `field` over `cascade` at `time`. Reuses its own storage; the
    // index buffer is rebuilt only when the topology changes.
    void build(const sim::WaveField& field, const OceanCascade& cascade, double time);

    const OceanCascade& cascade() const { return cascade_; }
    double time() const { return time_; }
    const std::vector<OceanVertex>& vertices() const { return vertices_; }
    const std::vector<std::uint32_t>& indices() const { return indices_; }
    const std::vector<Level>& levels() const { return levels_; }
    double buildSeconds() const { return buildSeconds_; }
    double reach() const { return cascade_.reach(); }

    // Index of the vertex at lattice coordinates (i, j) of `level`, both in
    // [-cells/2, +cells/2]. **A point on or inside the level's hole resolves to
    // the finer level that owns it**, which is what makes the seam shared rather
    // than duplicated, and is why this is the only way to address a vertex.
    std::size_t vertexIndexAt(int level, int i, int j) const;

    // Elevation of the *mesh* at a world point: the piecewise-linear surface the
    // rasteriser interpolates, transition triangles included. False outside the
    // cascade.
    bool sampleElevation(double x, double y, double& elevation) const;
    // Which level owns a world point; -1 outside the cascade.
    int levelAt(double x, double y) const;

private:
    void layOutVertices();
    void buildIndices();

    OceanCascade cascade_;
    int cells_ = 0;
    double time_ = 0.0;
    double buildSeconds_ = 0.0;
    // What the cached index buffer was built for.
    int indexCells_ = -1, indexLevels_ = -1;
    bool indexStitched_ = true;

    std::vector<Level> levels_;
    // Counted while the index buffer was written, and therefore surviving the
    // rebuilds the index buffer is cached through.
    std::vector<std::size_t> levelTriangles_;
    std::vector<OceanVertex> vertices_;
    std::vector<std::uint32_t> indices_;
    // Absolute offset of each stored row, level-major. Rows a level does not own
    // in full are stored as two segments; see the .cpp.
    std::vector<std::size_t> rowStart_;
    // The field's components ordered by wavenumber, so "the ones a level can
    // resolve" is a prefix and the level sets are nested by construction.
    std::vector<sim::WaveComponent> byWavenumber_;
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
    bool render(const float mvp[16], const OceanCascadeSurface& surface, const OceanView& view,
                const float clearColour[4], core::Image& out);
    // The same, over geometry from wherever. The two above are this one.
    bool render(const float mvp[16], const std::vector<OceanVertex>& vertices,
                const std::vector<std::uint32_t>& indices, const OceanView& view,
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
