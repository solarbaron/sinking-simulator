// SPDX-License-Identifier: MIT
#include "ocean.hpp"

#include <vulkan/vulkan.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <vector>

namespace gpu {
namespace {

constexpr VkFormat kColourFormat = VK_FORMAT_R8G8B8A8_UNORM;
constexpr VkFormat kDepthFormat = VK_FORMAT_D32_SFLOAT;

// A patch is a square grid, so the vertex count is the square of this. 2048 is
// 4.2 M vertices, well past anything the CPU path can displace in a frame; it
// exists to turn a nonsense argument into a clamp rather than an allocation
// failure.
constexpr int kMaxResolution = 2048;

// Mirrors the `Push` block in shaders/ocean.vert and shaders/ocean.frag. 112
// bytes, inside the 128 every Vulkan implementation guarantees, so the ocean
// needs no descriptor set at all.
struct Push {
    float modelViewProjection[16];
    float sun[4];     // xyz unit vector toward the sun, w strength
    float water[4];   // rgb water colour, w ambient
    float encode[4];  // x mode (0 shaded, 1 elevation), y zMin, z 1 / span, w unused
};
static_assert(sizeof(Push) == 112, "the push block must match the shaders");

// One grid row's elevation and slope, accumulated over `count` components.
//
// **The row recurrence lives here and nowhere else.** Along a row only x moves, so
// the phase advances by a constant and (cos psi, sin psi) can be stepped by a
// fixed rotation instead of re-evaluated: one sincos per row per component rather
// than one per vertex, and the sine that falls out is exactly what the slope
// needs, so the normal is free. It is an algebraic rearrangement of the same sum,
// not an approximation -- worst measured disagreement against direct evaluation
// over a 257-wide row is 2e-14 m -- and it is an eighth of the cost. Restarted
// from a real sincos on every row and every segment, so drift is bounded by the
// run length rather than by the whole grid.
void accumulateRow(const sim::WaveComponent* components, std::size_t count, double x0, double y,
                   double h, int points, double time, double* elevation, double* slopeX,
                   double* slopeY) {
    std::fill(elevation, elevation + points, 0.0);
    std::fill(slopeX, slopeX + points, 0.0);
    std::fill(slopeY, slopeY + points, 0.0);

    for (std::size_t c = 0; c < count; ++c) {
        const sim::WaveComponent& component = components[c];
        const double kx = component.wavenumber * component.dirX;
        const double ky = component.wavenumber * component.dirY;
        // eta = a cos(psi), so d eta/dx = -a kx sin(psi) and likewise in y.
        const double a = component.amplitude;
        const double slopeGainX = -a * kx;
        const double slopeGainY = -a * ky;

        const double psi = kx * x0 + ky * y - component.omega * time + component.phase;
        double cosPsi = std::cos(psi);
        double sinPsi = std::sin(psi);
        const double stepCos = std::cos(kx * h);
        const double stepSin = std::sin(kx * h);

        for (int i = 0; i < points; ++i) {
            const auto index = static_cast<std::size_t>(i);
            elevation[index] += a * cosPsi;
            slopeX[index] += slopeGainX * sinPsi;
            slopeY[index] += slopeGainY * sinPsi;
            const double nextCos = cosPsi * stepCos - sinPsi * stepSin;
            sinPsi = sinPsi * stepCos + cosPsi * stepSin;
            cosPsi = nextCos;
        }
    }
}

// n = normalize(-d eta/dx, -d eta/dy, 1), written into a vertex at (x, y).
void writeVertex(OceanVertex& vertex, double x, double y, double z, double slopeX, double slopeY) {
    const double nx = -slopeX;
    const double ny = -slopeY;
    const double inverse = 1.0 / std::sqrt(nx * nx + ny * ny + 1.0);
    vertex.position[0] = static_cast<float>(x);
    vertex.position[1] = static_cast<float>(y);
    vertex.position[2] = static_cast<float>(z);
    vertex.normal[0] = static_cast<float>(nx * inverse);
    vertex.normal[1] = static_cast<float>(ny * inverse);
    vertex.normal[2] = static_cast<float>(inverse);
}

// A cascade level's ring has a hole a quarter of its extent across, and its
// boundary vertices have to land on even lattice indices of the finer level
// inside it. Both need the cell count to be a multiple of four.
constexpr int kMaxCascadeLevels = 24;

}  // namespace

// --- Resolution ---------------------------------------------------------------

int oceanResolutionFor(double halfExtent, double wavelength, double cellsPerWavelength) {
    if (!(halfExtent > 0.0) || !(wavelength > 0.0) || !(cellsPerWavelength > 0.0)) return 1;
    const double cells = 2.0 * halfExtent * cellsPerWavelength / wavelength;
    if (!(cells > 1.0)) return 1;
    if (cells >= kMaxResolution) return kMaxResolution;
    return static_cast<int>(std::ceil(cells));
}

double dominantWavelength(const sim::WaveField& field) {
    // Equal-energy bins, so the peak of the density is the narrowest bin. Same
    // route tests/test_waves.cpp uses to find Tp in the discretisation.
    double bestDensity = -1.0;
    double bestOmega = 0.0;
    for (const sim::FrequencyBin& bin : field.frequencyBins()) {
        const double width = bin.omegaHigh - bin.omegaLow;
        if (!std::isfinite(width) || !(width > 0.0)) continue;  // the open-topped bin
        const double density = bin.energy / width;
        if (density > bestDensity) {
            bestDensity = density;
            bestOmega = bin.omega;
        }
    }
    if (bestDensity < 0.0)
        // Every bin is open-topped, which happens only when a sea state has a
        // single frequency: that bin runs from zero to infinity and so has no
        // density at all. Its centroid is the only frequency there is.
        for (const sim::FrequencyBin& bin : field.frequencyBins())
            if (bin.omega > 0.0) {
                bestOmega = bin.omega;
                break;
            }
    const double k = sim::deepWaterWavenumber(bestOmega);
    return k > 0.0 ? 2.0 * sim::kPi / k : 0.0;
}

double shortestWavelength(const sim::WaveField& field) {
    double maxWavenumber = 0.0;
    for (const sim::WaveComponent& component : field.components())
        maxWavenumber = std::max(maxWavenumber, component.wavenumber);
    return maxWavenumber > 0.0 ? 2.0 * sim::kPi / maxWavenumber : 0.0;
}

// --- The displaced grid -------------------------------------------------------

std::size_t OceanSurface::vertexIndex(int i, int j) const {
    return static_cast<std::size_t>(j) * static_cast<std::size_t>(grid_.resolution + 1) +
           static_cast<std::size_t>(i);
}

void OceanSurface::build(const sim::WaveField& field, const OceanGrid& grid, double time) {
    const auto started = std::chrono::steady_clock::now();

    grid_ = grid;
    grid_.resolution = std::clamp(grid.resolution, 1, kMaxResolution);
    time_ = time;

    const int cells = grid_.resolution;
    const int side = cells + 1;
    const double h = grid_.cellSize();
    const double x0 = grid_.centreX - grid_.halfExtent;
    const double y0 = grid_.centreY - grid_.halfExtent;

    vertices_.resize(static_cast<std::size_t>(side) * static_cast<std::size_t>(side));

    if (indexResolution_ != cells) {
        indices_.clear();
        indices_.reserve(static_cast<std::size_t>(cells) * static_cast<std::size_t>(cells) * 6);
        for (int j = 0; j < cells; ++j)
            for (int i = 0; i < cells; ++i) {
                const auto v00 = static_cast<std::uint32_t>(vertexIndex(i, j));
                const auto v10 = static_cast<std::uint32_t>(v00 + 1);
                const auto v01 = static_cast<std::uint32_t>(v00 + static_cast<std::uint32_t>(side));
                const auto v11 = static_cast<std::uint32_t>(v01 + 1);
                // Split on the v00-v11 diagonal, counter-clockwise seen from
                // above. The choice is load-bearing for the tests: the cell
                // centre then lies on the shared edge, where the interpolated
                // height is the mean of the two diagonal corners -- which is what
                // makes the linear-interpolation error a closed form rather than
                // a bound.
                indices_.push_back(v00);
                indices_.push_back(v10);
                indices_.push_back(v11);
                indices_.push_back(v00);
                indices_.push_back(v11);
                indices_.push_back(v01);
            }
        indexResolution_ = cells;
    }

    rowElevation_.resize(static_cast<std::size_t>(side));
    rowSlopeX_.resize(static_cast<std::size_t>(side));
    rowSlopeY_.resize(static_cast<std::size_t>(side));

    const std::vector<sim::WaveComponent>& components = field.components();

    // One evaluation per *unique grid vertex*. Six triangle corners share each
    // interior vertex, so evaluating per index instead would multiply the whole
    // cost by six -- which is exactly the redundancy the physics tick was just
    // caught paying (CLAUDE.md, "Sea surface queried 6x more than necessary").
    for (int j = 0; j < side; ++j) {
        const double y = y0 + h * j;
        accumulateRow(components.data(), components.size(), x0, y, h, side, time,
                      rowElevation_.data(), rowSlopeX_.data(), rowSlopeY_.data());
        for (int i = 0; i < side; ++i) {
            const auto index = static_cast<std::size_t>(i);
            writeVertex(vertices_[vertexIndex(i, j)], x0 + h * i, y, rowElevation_[index],
                        rowSlopeX_[index], rowSlopeY_[index]);
        }
    }

    buildSeconds_ = std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
}

bool OceanSurface::sampleElevation(double x, double y, double& elevation) const {
    const int cells = grid_.resolution;
    if (vertices_.empty() || cells < 1) return false;
    const double h = grid_.cellSize();
    const double fx = (x - (grid_.centreX - grid_.halfExtent)) / h;
    const double fy = (y - (grid_.centreY - grid_.halfExtent)) / h;
    // Written so a NaN falls out of the patch rather than into an index.
    if (!(fx >= 0.0) || !(fy >= 0.0) || !(fx <= cells) || !(fy <= cells)) return false;

    int i = static_cast<int>(std::floor(fx));
    int j = static_cast<int>(std::floor(fy));
    if (i >= cells) i = cells - 1;
    if (j >= cells) j = cells - 1;
    const double u = fx - i;
    const double v = fy - j;

    const auto height = [&](int a, int b) {
        return static_cast<double>(vertices_[vertexIndex(a, b)].position[2]);
    };
    const double z00 = height(i, j);
    const double z11 = height(i + 1, j + 1);
    // The v00-v11 diagonal splits the cell along u == v; below it the triangle is
    // (v00, v10, v11) and above it (v00, v11, v01). Both agree on the diagonal,
    // where the value is (z00 + z11) / 2.
    elevation = v <= u ? z00 + u * (height(i + 1, j) - z00) + v * (z11 - height(i + 1, j))
                       : z00 + v * (height(i, j + 1) - z00) + u * (z11 - height(i, j + 1));
    return true;
}

// --- The cascade --------------------------------------------------------------

namespace {

// Which edge of cell (i, j) of `level` lies on the seam with the finer level
// inside it. **At most one can**, and that is what makes the whole construction
// tractable: the four cells at the corners of the hole touch it at a single point
// rather than along an edge, so there is no corner case to get wrong.
enum class SeamEdge { None, Bottom, Top, Left, Right };

SeamEdge seamEdgeOf(int cells, int level, int i, int j, bool stitch) {
    if (level < 1 || !stitch) return SeamEdge::None;
    const int hole = cells / 4;
    if (j == hole && i >= -hole && i <= hole - 1) return SeamEdge::Bottom;
    if (j == -hole - 1 && i >= -hole && i <= hole - 1) return SeamEdge::Top;
    if (i == hole && j >= -hole && j <= hole - 1) return SeamEdge::Left;
    if (i == -hole - 1 && j >= -hole && j <= hole - 1) return SeamEdge::Right;
    return SeamEdge::None;
}

// True when cell (i, j) of `level` is inside the hole the finer level fills.
bool cellIsInHole(int cells, int level, int i, int j) {
    if (level < 1) return false;
    const int hole = cells / 4;
    return i >= -hole && i <= hole - 1 && j >= -hole && j <= hole - 1;
}

}  // namespace

int OceanCascade::cellsPerSide() const {
    const int atLeast = std::max(resolution, 4);
    return std::min((atLeast + 3) / 4 * 4, kMaxResolution);
}

double OceanCascade::halfExtent(int level) const {
    const int clamped = std::clamp(level, 0, kMaxCascadeLevels - 1);
    return innerHalfExtent * static_cast<double>(1u << static_cast<unsigned>(clamped));
}

double OceanCascade::cellSize(int level) const {
    return 2.0 * halfExtent(level) / cellsPerSide();
}

int oceanCascadeLevelsFor(double innerHalfExtent, double reach) {
    if (!(innerHalfExtent > 0.0) || !(reach > innerHalfExtent)) return 1;
    // The tolerance is there so an exact power of two asks for the level that
    // reaches it, not for one more.
    const int doublings = static_cast<int>(std::ceil(std::log2(reach / innerHalfExtent) - 1e-9));
    return std::clamp(1 + doublings, 1, kMaxCascadeLevels);
}

double oceanHorizonReach(double eyeHeight, double verticalFov, int pixelHeight) {
    if (!(eyeHeight > 0.0) || !(verticalFov > 0.0) || !(verticalFov < sim::kPi) || pixelHeight < 1)
        return 0.0;
    return eyeHeight * pixelHeight / (2.0 * std::tan(0.5 * verticalFov));
}

std::size_t OceanCascadeSurface::vertexIndexAt(int level, int i, int j) const {
    const int half = cells_ / 2, hole = cells_ / 4;
    // A lattice point on or inside the hole belongs to the finer level, at twice
    // the index. This is the whole of the anti-cracking scheme: there is no second
    // copy of a seam vertex to disagree with the first.
    while (level > 0 && std::max(std::abs(i), std::abs(j)) <= hole) {
        i *= 2;
        j *= 2;
        --level;
    }
    const std::size_t row =
        rowStart_[static_cast<std::size_t>(level) * static_cast<std::size_t>(cells_ + 1) +
                  static_cast<std::size_t>(j + half)];
    if (level == 0 || std::abs(j) > hole) return row + static_cast<std::size_t>(i + half);
    // A row crossing the hole is stored as two segments, left then right.
    return i < 0 ? row + static_cast<std::size_t>(i + half)
                 : row + static_cast<std::size_t>(cells_ / 4) +
                       static_cast<std::size_t>(i - (hole + 1));
}

void OceanCascadeSurface::layOutVertices() {
    const int half = cells_ / 2, hole = cells_ / 4;
    rowStart_.assign(static_cast<std::size_t>(cascade_.levels) *
                         static_cast<std::size_t>(cells_ + 1),
                     0);
    std::size_t offset = 0;
    for (int level = 0; level < cascade_.levels; ++level) {
        const std::size_t levelBase = offset;
        for (int r = 0; r <= cells_; ++r) {
            rowStart_[static_cast<std::size_t>(level) * static_cast<std::size_t>(cells_ + 1) +
                      static_cast<std::size_t>(r)] = offset;
            const int j = r - half;
            const bool full = level == 0 || std::abs(j) > hole;
            offset += full ? static_cast<std::size_t>(cells_ + 1)
                           : static_cast<std::size_t>(cells_ / 2);
        }
        levels_[static_cast<std::size_t>(level)].vertices = offset - levelBase;
    }
    vertices_.resize(offset);
}

void OceanCascadeSurface::buildIndices() {
    if (indexCells_ == cells_ && indexLevels_ == cascade_.levels &&
        indexStitched_ == cascade_.stitchSeams)
        return;

    const int half = cells_ / 2;
    indices_.clear();
    levelTriangles_.assign(static_cast<std::size_t>(cascade_.levels), 0);
    for (int level = 0; level < cascade_.levels; ++level) {
        const std::size_t before = indices_.size();
        const auto index = [&](int i, int j) {
            return static_cast<std::uint32_t>(vertexIndexAt(level, i, j));
        };
        const auto finer = [&](int i, int j) {
            return static_cast<std::uint32_t>(vertexIndexAt(level - 1, i, j));
        };
        const auto triangle = [&](std::uint32_t a, std::uint32_t b, std::uint32_t c) {
            indices_.push_back(a);
            indices_.push_back(b);
            indices_.push_back(c);
        };

        for (int j = -half; j < half; ++j)
            for (int i = -half; i < half; ++i) {
                if (cellIsInHole(cells_, level, i, j)) continue;
                const std::uint32_t v00 = index(i, j);
                const std::uint32_t v10 = index(i + 1, j);
                const std::uint32_t v11 = index(i + 1, j + 1);
                const std::uint32_t v01 = index(i, j + 1);
                // Counter-clockwise seen from above in every case, and the
                // no-seam split is the v00-v11 diagonal OceanSurface uses, so a
                // cell centre still lands on the shared edge.
                switch (seamEdgeOf(cells_, level, i, j, cascade_.stitchSeams)) {
                    case SeamEdge::None:
                        triangle(v00, v10, v11);
                        triangle(v00, v11, v01);
                        break;
                    case SeamEdge::Bottom: {
                        const std::uint32_t m = finer(2 * i + 1, 2 * j);
                        triangle(v00, m, v01);
                        triangle(m, v11, v01);
                        triangle(m, v10, v11);
                        break;
                    }
                    case SeamEdge::Top: {
                        const std::uint32_t m = finer(2 * i + 1, 2 * (j + 1));
                        triangle(v00, v10, m);
                        triangle(v00, m, v01);
                        triangle(v10, v11, m);
                        break;
                    }
                    case SeamEdge::Left: {
                        const std::uint32_t m = finer(2 * i, 2 * j + 1);
                        triangle(v00, v10, m);
                        triangle(m, v10, v11);
                        triangle(m, v11, v01);
                        break;
                    }
                    case SeamEdge::Right: {
                        const std::uint32_t m = finer(2 * (i + 1), 2 * j + 1);
                        triangle(v00, v10, m);
                        triangle(v00, m, v11);
                        triangle(v00, v11, v01);
                        break;
                    }
                }
            }
        levelTriangles_[static_cast<std::size_t>(level)] = (indices_.size() - before) / 3;
    }
    indexCells_ = cells_;
    indexLevels_ = cascade_.levels;
    indexStitched_ = cascade_.stitchSeams;
}

void OceanCascadeSurface::build(const sim::WaveField& field, const OceanCascade& cascade,
                                double time) {
    const auto started = std::chrono::steady_clock::now();

    cascade_ = cascade;
    cascade_.levels = std::clamp(cascade.levels, 1, kMaxCascadeLevels);
    cells_ = cascade_.cellsPerSide();
    cascade_.resolution = cells_;  // report what was used, not what was asked for
    time_ = time;

    const int half = cells_ / 2, hole = cells_ / 4;

    // Sorted by wavenumber, so "the components this level resolves" is a prefix
    // and the level sets nest: a component a ring carries is carried by every
    // finer level too, which is what stops a seam from being a step in the
    // spectrum as well as in the mesh.
    byWavenumber_ = field.components();
    std::sort(byWavenumber_.begin(), byWavenumber_.end(),
              [](const sim::WaveComponent& a, const sim::WaveComponent& b) {
                  return a.wavenumber < b.wavenumber;
              });

    levels_.assign(static_cast<std::size_t>(cascade_.levels), Level{});
    const double cellsPerWavelength = std::max(cascade_.minimumCellsPerWavelength, 1e-9);
    for (int level = 0; level < cascade_.levels; ++level) {
        Level& info = levels_[static_cast<std::size_t>(level)];
        info.cellSize = cascade_.cellSize(level);
        info.halfExtent = cascade_.halfExtent(level);
        // wavelength >= cellsPerWavelength * h  <=>  k <= 2 pi / (cellsPerWavelength * h)
        const double maxWavenumber = 2.0 * sim::kPi / (cellsPerWavelength * info.cellSize);
        const auto end = std::upper_bound(byWavenumber_.begin(), byWavenumber_.end(), maxWavenumber,
                                          [](double bound, const sim::WaveComponent& c) {
                                              return bound < c.wavenumber;
                                          });
        info.components = static_cast<std::size_t>(end - byWavenumber_.begin());
    }

    layOutVertices();
    buildIndices();
    // Counted while the indices were written, not predicted from the shape: a
    // reported figure that comes from a second derivation of the same thing is a
    // figure that can disagree with the geometry. The closed forms live in
    // tests/test_ocean.cpp, where disagreeing with them is the point.
    for (int level = 0; level < cascade_.levels; ++level)
        levels_[static_cast<std::size_t>(level)].triangles =
            levelTriangles_[static_cast<std::size_t>(level)];

    rowElevation_.resize(static_cast<std::size_t>(cells_ + 1));
    rowSlopeX_.resize(static_cast<std::size_t>(cells_ + 1));
    rowSlopeY_.resize(static_cast<std::size_t>(cells_ + 1));

    // One segment of one row: the recurrence, then the vertices it produced.
    const auto segment = [&](std::size_t out, int i0, int points, double h, double y,
                             std::size_t componentCount) {
        const double x0 = cascade_.centreX + static_cast<double>(i0) * h;
        if (componentCount == 0) {
            // A level that resolves nothing is dead flat, and skipping the
            // recurrence entirely is what makes the far rings almost free.
            for (int p = 0; p < points; ++p)
                writeVertex(vertices_[out + static_cast<std::size_t>(p)], x0 + h * p, y, 0.0, 0.0,
                            0.0);
            return;
        }
        accumulateRow(byWavenumber_.data(), componentCount, x0, y, h, points, time,
                      rowElevation_.data(), rowSlopeX_.data(), rowSlopeY_.data());
        for (int p = 0; p < points; ++p) {
            const auto index = static_cast<std::size_t>(p);
            writeVertex(vertices_[out + index], x0 + h * p, y, rowElevation_[index],
                        rowSlopeX_[index], rowSlopeY_[index]);
        }
    };

    for (int level = 0; level < cascade_.levels; ++level) {
        const double h = levels_[static_cast<std::size_t>(level)].cellSize;
        const std::size_t componentCount = levels_[static_cast<std::size_t>(level)].components;
        for (int r = 0; r <= cells_; ++r) {
            const int j = r - half;
            const double y = cascade_.centreY + static_cast<double>(j) * h;
            const std::size_t rowBase =
                rowStart_[static_cast<std::size_t>(level) * static_cast<std::size_t>(cells_ + 1) +
                          static_cast<std::size_t>(r)];
            if (level == 0 || std::abs(j) > hole) {
                segment(rowBase, -half, cells_ + 1, h, y, componentCount);
            } else {
                // The hole and its boundary belong to the finer level, so this row
                // is two runs with the middle missing.
                segment(rowBase, -half, cells_ / 4, h, y, componentCount);
                segment(rowBase + static_cast<std::size_t>(cells_ / 4), hole + 1, cells_ / 4, h, y,
                        componentCount);
            }
        }
    }

    buildSeconds_ = std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
}

int OceanCascadeSurface::levelAt(double x, double y) const {
    if (levels_.empty()) return -1;
    const double radius = std::max(std::abs(x - cascade_.centreX), std::abs(y - cascade_.centreY));
    // Written so a NaN falls out of the cascade rather than into a level.
    if (!(radius >= 0.0)) return -1;
    for (int level = 0; level < cascade_.levels; ++level)
        if (radius <= levels_[static_cast<std::size_t>(level)].halfExtent) return level;
    return -1;
}

bool OceanCascadeSurface::sampleElevation(double x, double y, double& elevation) const {
    const int level = levelAt(x, y);
    if (level < 0 || vertices_.empty()) return false;

    const int half = cells_ / 2;
    const double h = levels_[static_cast<std::size_t>(level)].cellSize;
    const double fi = (x - cascade_.centreX) / h;
    const double fj = (y - cascade_.centreY) / h;
    const int i = std::clamp(static_cast<int>(std::floor(fi)), -half, half - 1);
    const int j = std::clamp(static_cast<int>(std::floor(fj)), -half, half - 1);
    const double u = fi - i, v = fj - j;

    const auto height = [&](int a, int b) {
        return static_cast<double>(vertices_[vertexIndexAt(level, a, b)].position[2]);
    };
    const double z00 = height(i, j), z10 = height(i + 1, j);
    const double z11 = height(i + 1, j + 1), z01 = height(i, j + 1);

    const SeamEdge seam = seamEdgeOf(cells_, level, i, j, cascade_.stitchSeams);
    if (seam == SeamEdge::None) {
        elevation = v <= u ? z00 + u * (z10 - z00) + v * (z11 - z10)
                           : z00 + v * (z01 - z00) + u * (z11 - z01);
        return true;
    }

    // A transition cell is three triangles, so the containing one has to be found
    // rather than chosen by a diagonal test. The midpoint comes from the finer
    // level, which is the point of the whole arrangement.
    const auto midpoint = [&](int a, int b) {
        return static_cast<double>(vertices_[vertexIndexAt(level - 1, a, b)].position[2]);
    };
    double soup[9][3];  // three triangles, nine corners, each (u, v, z)
    const auto put = [&](int slot, double cu, double cv, double cz) {
        soup[slot][0] = cu;
        soup[slot][1] = cv;
        soup[slot][2] = cz;
    };
    switch (seam) {
        case SeamEdge::Bottom: {
            const double zm = midpoint(2 * i + 1, 2 * j);
            put(0, 0, 0, z00); put(1, 0.5, 0, zm);  put(2, 0, 1, z01);
            put(3, 0.5, 0, zm); put(4, 1, 1, z11);  put(5, 0, 1, z01);
            put(6, 0.5, 0, zm); put(7, 1, 0, z10);  put(8, 1, 1, z11);
            break;
        }
        case SeamEdge::Top: {
            const double zm = midpoint(2 * i + 1, 2 * (j + 1));
            put(0, 0, 0, z00); put(1, 1, 0, z10);   put(2, 0.5, 1, zm);
            put(3, 0, 0, z00); put(4, 0.5, 1, zm);  put(5, 0, 1, z01);
            put(6, 1, 0, z10); put(7, 1, 1, z11);   put(8, 0.5, 1, zm);
            break;
        }
        case SeamEdge::Left: {
            const double zm = midpoint(2 * i, 2 * j + 1);
            put(0, 0, 0, z00); put(1, 1, 0, z10);   put(2, 0, 0.5, zm);
            put(3, 0, 0.5, zm); put(4, 1, 0, z10);  put(5, 1, 1, z11);
            put(6, 0, 0.5, zm); put(7, 1, 1, z11);  put(8, 0, 1, z01);
            break;
        }
        case SeamEdge::Right: {
            const double zm = midpoint(2 * (i + 1), 2 * j + 1);
            put(0, 0, 0, z00); put(1, 1, 0, z10);   put(2, 1, 0.5, zm);
            put(3, 0, 0, z00); put(4, 1, 0.5, zm);  put(5, 1, 1, z11);
            put(6, 0, 0, z00); put(7, 1, 1, z11);   put(8, 0, 1, z01);
            break;
        }
        case SeamEdge::None:
            break;
    }

    for (int t = 0; t < 3; ++t) {
        const double* a = soup[3 * t];
        const double* b = soup[3 * t + 1];
        const double* c = soup[3 * t + 2];
        const double area = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
        if (std::abs(area) < 1e-12) continue;
        const double wa = ((b[0] - u) * (c[1] - v) - (b[1] - v) * (c[0] - u)) / area;
        const double wb = ((c[0] - u) * (a[1] - v) - (c[1] - v) * (a[0] - u)) / area;
        const double wc = 1.0 - wa - wb;
        if (wa < -1e-9 || wb < -1e-9 || wc < -1e-9) continue;
        elevation = wa * a[2] + wb * b[2] + wc * c[2];
        return true;
    }
    return false;
}

// --- The elevation channel ----------------------------------------------------

bool decodeOceanElevation(const OceanView& view, const std::uint8_t* pixel, double& elevation) {
    // Blue is the surface tag. The clear colour has to stay away from it, which
    // is cheap to arrange and much safer than trying to tell a legitimate
    // elevation code apart from the background.
    if (pixel == nullptr || pixel[2] < 128) return false;
    const double code = static_cast<double>(pixel[0]) * 256.0 + static_cast<double>(pixel[1]);
    elevation = static_cast<double>(view.elevationMin) +
                code / 65535.0 * static_cast<double>(view.elevationSpan);
    return true;
}

// --- The pipeline -------------------------------------------------------------

OceanRenderer::~OceanRenderer() { destroy(); }

bool OceanRenderer::create(Device& device, std::uint32_t width, std::uint32_t height,
                           const std::string& shaderDirectory, std::string& error) {
    destroy();
    if (!device.valid()) {
        error = "device is not valid";
        return false;
    }
    if (!device.supportsGraphics()) {
        error = "queue family does not support graphics";
        return false;
    }
    device_ = &device;
    width_ = width;
    height_ = height;

    colour_ = device.createImage2D(width, height, kColourFormat,
                                   VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT |
                                       VK_IMAGE_USAGE_TRANSFER_SRC_BIT,
                                   VK_IMAGE_ASPECT_COLOR_BIT);
    depth_ = device.createImage2D(width, height, kDepthFormat,
                                  VK_IMAGE_USAGE_DEPTH_STENCIL_ATTACHMENT_BIT,
                                  VK_IMAGE_ASPECT_DEPTH_BIT);
    if (!colour_.valid() || !depth_.valid()) {
        error = "could not allocate the render target";
        destroy();
        return false;
    }

    // Same layout chain as OffscreenRenderer: the pass leaves the colour
    // attachment in TRANSFER_SRC_OPTIMAL so the readback copy needs no barrier.
    VkAttachmentDescription attachments[2]{};
    attachments[0].format = kColourFormat;
    attachments[0].samples = VK_SAMPLE_COUNT_1_BIT;
    attachments[0].loadOp = VK_ATTACHMENT_LOAD_OP_CLEAR;
    attachments[0].storeOp = VK_ATTACHMENT_STORE_OP_STORE;
    attachments[0].stencilLoadOp = VK_ATTACHMENT_LOAD_OP_DONT_CARE;
    attachments[0].stencilStoreOp = VK_ATTACHMENT_STORE_OP_DONT_CARE;
    attachments[0].initialLayout = VK_IMAGE_LAYOUT_UNDEFINED;
    attachments[0].finalLayout = VK_IMAGE_LAYOUT_TRANSFER_SRC_OPTIMAL;

    attachments[1].format = kDepthFormat;
    attachments[1].samples = VK_SAMPLE_COUNT_1_BIT;
    attachments[1].loadOp = VK_ATTACHMENT_LOAD_OP_CLEAR;
    attachments[1].storeOp = VK_ATTACHMENT_STORE_OP_DONT_CARE;
    attachments[1].stencilLoadOp = VK_ATTACHMENT_LOAD_OP_DONT_CARE;
    attachments[1].stencilStoreOp = VK_ATTACHMENT_STORE_OP_DONT_CARE;
    attachments[1].initialLayout = VK_IMAGE_LAYOUT_UNDEFINED;
    attachments[1].finalLayout = VK_IMAGE_LAYOUT_DEPTH_STENCIL_ATTACHMENT_OPTIMAL;

    VkAttachmentReference colourRef{0, VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL};
    VkAttachmentReference depthRef{1, VK_IMAGE_LAYOUT_DEPTH_STENCIL_ATTACHMENT_OPTIMAL};

    VkSubpassDescription subpass{};
    subpass.pipelineBindPoint = VK_PIPELINE_BIND_POINT_GRAPHICS;
    subpass.colorAttachmentCount = 1;
    subpass.pColorAttachments = &colourRef;
    subpass.pDepthStencilAttachment = &depthRef;

    VkRenderPassCreateInfo renderPassInfo{VK_STRUCTURE_TYPE_RENDER_PASS_CREATE_INFO};
    renderPassInfo.attachmentCount = 2;
    renderPassInfo.pAttachments = attachments;
    renderPassInfo.subpassCount = 1;
    renderPassInfo.pSubpasses = &subpass;
    if (vkCreateRenderPass(device.handle(), &renderPassInfo, nullptr, &renderPass_) != VK_SUCCESS) {
        error = "vkCreateRenderPass failed";
        destroy();
        return false;
    }

    VkImageView views[2] = {colour_.view, depth_.view};
    VkFramebufferCreateInfo framebufferInfo{VK_STRUCTURE_TYPE_FRAMEBUFFER_CREATE_INFO};
    framebufferInfo.renderPass = renderPass_;
    framebufferInfo.attachmentCount = 2;
    framebufferInfo.pAttachments = views;
    framebufferInfo.width = width;
    framebufferInfo.height = height;
    framebufferInfo.layers = 1;
    if (vkCreateFramebuffer(device.handle(), &framebufferInfo, nullptr, &framebuffer_) !=
        VK_SUCCESS) {
        error = "vkCreateFramebuffer failed";
        destroy();
        return false;
    }

    VkShaderModule vertexModule = device.loadShader(shaderDirectory + "/ocean.vert.spv", error);
    if (vertexModule == nullptr) { destroy(); return false; }
    VkShaderModule fragmentModule = device.loadShader(shaderDirectory + "/ocean.frag.spv", error);
    if (fragmentModule == nullptr) {
        device.destroyShader(vertexModule);
        destroy();
        return false;
    }

    // One range spanning both stages: the vertex shader reads the matrix, the
    // fragment shader the lighting and the encoding, and a single range keeps the
    // block declaration identical in the two shaders.
    VkPushConstantRange push{VK_SHADER_STAGE_VERTEX_BIT | VK_SHADER_STAGE_FRAGMENT_BIT, 0,
                             sizeof(Push)};
    VkPipelineLayoutCreateInfo layoutInfo{VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO};
    layoutInfo.pushConstantRangeCount = 1;
    layoutInfo.pPushConstantRanges = &push;
    if (vkCreatePipelineLayout(device.handle(), &layoutInfo, nullptr, &pipelineLayout_) !=
        VK_SUCCESS) {
        error = "vkCreatePipelineLayout failed";
        device.destroyShader(vertexModule);
        device.destroyShader(fragmentModule);
        destroy();
        return false;
    }

    VkPipelineShaderStageCreateInfo stages[2]{};
    stages[0].sType = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO;
    stages[0].stage = VK_SHADER_STAGE_VERTEX_BIT;
    stages[0].module = vertexModule;
    stages[0].pName = "main";
    stages[1].sType = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO;
    stages[1].stage = VK_SHADER_STAGE_FRAGMENT_BIT;
    stages[1].module = fragmentModule;
    stages[1].pName = "main";

    VkVertexInputBindingDescription binding{0, sizeof(OceanVertex), VK_VERTEX_INPUT_RATE_VERTEX};
    VkVertexInputAttributeDescription vertexAttributes[2]{};
    vertexAttributes[0] = {0, 0, VK_FORMAT_R32G32B32_SFLOAT, offsetof(OceanVertex, position)};
    vertexAttributes[1] = {1, 0, VK_FORMAT_R32G32B32_SFLOAT, offsetof(OceanVertex, normal)};

    VkPipelineVertexInputStateCreateInfo vertexInput{
        VK_STRUCTURE_TYPE_PIPELINE_VERTEX_INPUT_STATE_CREATE_INFO};
    vertexInput.vertexBindingDescriptionCount = 1;
    vertexInput.pVertexBindingDescriptions = &binding;
    vertexInput.vertexAttributeDescriptionCount = 2;
    vertexInput.pVertexAttributeDescriptions = vertexAttributes;

    VkPipelineInputAssemblyStateCreateInfo assembly{
        VK_STRUCTURE_TYPE_PIPELINE_INPUT_ASSEMBLY_STATE_CREATE_INFO};
    assembly.topology = VK_PRIMITIVE_TOPOLOGY_TRIANGLE_LIST;

    VkViewport viewport{0.0f, 0.0f, static_cast<float>(width), static_cast<float>(height),
                        0.0f, 1.0f};
    VkRect2D scissor{{0, 0}, {width, height}};
    VkPipelineViewportStateCreateInfo viewportState{
        VK_STRUCTURE_TYPE_PIPELINE_VIEWPORT_STATE_CREATE_INFO};
    viewportState.viewportCount = 1;
    viewportState.pViewports = &viewport;
    viewportState.scissorCount = 1;
    viewportState.pScissors = &scissor;

    VkPipelineRasterizationStateCreateInfo raster{
        VK_STRUCTURE_TYPE_PIPELINE_RASTERIZATION_STATE_CREATE_INFO};
    raster.polygonMode = VK_POLYGON_MODE_FILL;
    // The sea is seen from below as often as from above -- from a flooded
    // compartment, from underwater -- so neither face is a back face. The
    // fragment shader flips the normal for the far side instead.
    raster.cullMode = VK_CULL_MODE_NONE;
    raster.frontFace = VK_FRONT_FACE_COUNTER_CLOCKWISE;
    raster.lineWidth = 1.0f;

    VkPipelineMultisampleStateCreateInfo multisample{
        VK_STRUCTURE_TYPE_PIPELINE_MULTISAMPLE_STATE_CREATE_INFO};
    multisample.rasterizationSamples = VK_SAMPLE_COUNT_1_BIT;

    VkPipelineDepthStencilStateCreateInfo depthStencil{
        VK_STRUCTURE_TYPE_PIPELINE_DEPTH_STENCIL_STATE_CREATE_INFO};
    depthStencil.depthTestEnable = VK_TRUE;
    depthStencil.depthWriteEnable = VK_TRUE;
    depthStencil.depthCompareOp = VK_COMPARE_OP_LESS;

    VkPipelineColorBlendAttachmentState blendAttachment{};
    blendAttachment.colorWriteMask = VK_COLOR_COMPONENT_R_BIT | VK_COLOR_COMPONENT_G_BIT |
                                     VK_COLOR_COMPONENT_B_BIT | VK_COLOR_COMPONENT_A_BIT;
    VkPipelineColorBlendStateCreateInfo blend{
        VK_STRUCTURE_TYPE_PIPELINE_COLOR_BLEND_STATE_CREATE_INFO};
    blend.attachmentCount = 1;
    blend.pAttachments = &blendAttachment;

    VkGraphicsPipelineCreateInfo pipelineInfo{VK_STRUCTURE_TYPE_GRAPHICS_PIPELINE_CREATE_INFO};
    pipelineInfo.stageCount = 2;
    pipelineInfo.pStages = stages;
    pipelineInfo.pVertexInputState = &vertexInput;
    pipelineInfo.pInputAssemblyState = &assembly;
    pipelineInfo.pViewportState = &viewportState;
    pipelineInfo.pRasterizationState = &raster;
    pipelineInfo.pMultisampleState = &multisample;
    pipelineInfo.pDepthStencilState = &depthStencil;
    pipelineInfo.pColorBlendState = &blend;
    pipelineInfo.layout = pipelineLayout_;
    pipelineInfo.renderPass = renderPass_;

    const VkResult created = vkCreateGraphicsPipelines(device.handle(), VK_NULL_HANDLE, 1,
                                                       &pipelineInfo, nullptr, &pipeline_);
    device.destroyShader(vertexModule);
    device.destroyShader(fragmentModule);
    if (created != VK_SUCCESS) {
        error = "vkCreateGraphicsPipelines failed";
        destroy();
        return false;
    }

    readback_ = device.createBuffer(std::uint64_t{width} * height * 4,
                                    VK_BUFFER_USAGE_TRANSFER_DST_BIT,
                                    VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT |
                                        VK_MEMORY_PROPERTY_HOST_COHERENT_BIT);
    if (!readback_.valid()) {
        error = "could not allocate the readback buffer";
        destroy();
        return false;
    }
    return true;
}

void OceanRenderer::destroy() {
    if (device_ == nullptr) return;
    VkDevice handle = device_->handle();
    if (handle != nullptr) {
        vkDeviceWaitIdle(handle);
        if (pipeline_ != nullptr) vkDestroyPipeline(handle, pipeline_, nullptr);
        if (pipelineLayout_ != nullptr) vkDestroyPipelineLayout(handle, pipelineLayout_, nullptr);
        if (framebuffer_ != nullptr) vkDestroyFramebuffer(handle, framebuffer_, nullptr);
        if (renderPass_ != nullptr) vkDestroyRenderPass(handle, renderPass_, nullptr);
        device_->destroyBuffer(vertexBuffer_);
        device_->destroyBuffer(indexBuffer_);
        device_->destroyBuffer(readback_);
        device_->destroyImage(colour_);
        device_->destroyImage(depth_);
    }
    pipeline_ = nullptr;
    pipelineLayout_ = nullptr;
    framebuffer_ = nullptr;
    renderPass_ = nullptr;
    device_ = nullptr;
    width_ = height_ = 0;
}

bool OceanRenderer::ensureGeometryCapacity(std::size_t vertexBytes, std::size_t indexBytes) {
    if (vertexBytes > vertexBuffer_.size) {
        device_->destroyBuffer(vertexBuffer_);
        vertexBuffer_ = device_->createBuffer(
            vertexBytes,
            VK_BUFFER_USAGE_VERTEX_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT,
            VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT);
        if (!vertexBuffer_.valid()) return false;
    }
    if (indexBytes > indexBuffer_.size) {
        device_->destroyBuffer(indexBuffer_);
        indexBuffer_ = device_->createBuffer(
            indexBytes, VK_BUFFER_USAGE_INDEX_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT,
            VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT);
        if (!indexBuffer_.valid()) return false;
    }
    return true;
}

bool OceanRenderer::render(const float mvp[16], const OceanSurface& surface, const OceanView& view,
                           const float clearColour[4], core::Image& out) {
    return render(mvp, surface.vertices(), surface.indices(), view, clearColour, out);
}

bool OceanRenderer::render(const float mvp[16], const OceanCascadeSurface& surface,
                           const OceanView& view, const float clearColour[4], core::Image& out) {
    return render(mvp, surface.vertices(), surface.indices(), view, clearColour, out);
}

bool OceanRenderer::render(const float mvp[16], const std::vector<OceanVertex>& surfaceVertices,
                           const std::vector<std::uint32_t>& surfaceIndices, const OceanView& view,
                           const float clearColour[4], core::Image& out) {
    if (!valid()) return false;

    const std::size_t vertexCount = surfaceVertices.size();
    const std::size_t indexCount = surfaceIndices.size();
    const bool hasGeometry = vertexCount > 0 && indexCount > 0;
    if (hasGeometry) {
        const std::size_t vertexBytes = vertexCount * sizeof(OceanVertex);
        const std::size_t indexBytes = indexCount * sizeof(std::uint32_t);
        if (!ensureGeometryCapacity(vertexBytes, indexBytes)) return false;
        if (!device_->upload(vertexBuffer_, surfaceVertices.data(), vertexBytes)) return false;
        if (!device_->upload(indexBuffer_, surfaceIndices.data(), indexBytes)) return false;
    }

    Push push{};
    std::memcpy(push.modelViewProjection, mvp, sizeof(push.modelViewProjection));
    // Normalised here rather than trusted: a caller who hands over a direction of
    // length 1.001 would otherwise get lighting that is subtly too bright, which
    // is the sort of thing nobody ever notices.
    const double length = std::sqrt(
        static_cast<double>(view.sunDirection[0]) * view.sunDirection[0] +
        static_cast<double>(view.sunDirection[1]) * view.sunDirection[1] +
        static_cast<double>(view.sunDirection[2]) * view.sunDirection[2]);
    const double inverse = length > 1e-12 ? 1.0 / length : 0.0;
    for (int i = 0; i < 3; ++i)
        push.sun[i] = static_cast<float>(view.sunDirection[i] * inverse);
    push.sun[3] = view.sunStrength;
    for (int i = 0; i < 3; ++i) push.water[i] = view.waterColour[i];
    push.water[3] = view.ambient;
    push.encode[0] = view.shading == OceanShading::Elevation ? 1.0f : 0.0f;
    push.encode[1] = view.elevationMin;
    push.encode[2] = view.elevationSpan != 0.0f ? 1.0f / view.elevationSpan : 0.0f;
    push.encode[3] = 0.0f;

    VkClearValue clears[2]{};
    for (int i = 0; i < 4; ++i) clears[0].color.float32[i] = clearColour[i];
    clears[1].depthStencil.depth = 1.0f;

    VkCommandBuffer commands = device_->beginOneShot();

    VkRenderPassBeginInfo begin{VK_STRUCTURE_TYPE_RENDER_PASS_BEGIN_INFO};
    begin.renderPass = renderPass_;
    begin.framebuffer = framebuffer_;
    begin.renderArea = {{0, 0}, {width_, height_}};
    begin.clearValueCount = 2;
    begin.pClearValues = clears;
    vkCmdBeginRenderPass(commands, &begin, VK_SUBPASS_CONTENTS_INLINE);

    if (hasGeometry) {
        vkCmdBindPipeline(commands, VK_PIPELINE_BIND_POINT_GRAPHICS, pipeline_);
        vkCmdPushConstants(commands, pipelineLayout_,
                           VK_SHADER_STAGE_VERTEX_BIT | VK_SHADER_STAGE_FRAGMENT_BIT, 0,
                           sizeof(Push), &push);
        const VkDeviceSize offset = 0;
        vkCmdBindVertexBuffers(commands, 0, 1, &vertexBuffer_.handle, &offset);
        vkCmdBindIndexBuffer(commands, indexBuffer_.handle, 0, VK_INDEX_TYPE_UINT32);
        vkCmdDrawIndexed(commands, static_cast<std::uint32_t>(indexCount), 1, 0, 0, 0);
    }
    vkCmdEndRenderPass(commands);

    VkBufferImageCopy copy{};
    copy.imageSubresource.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
    copy.imageSubresource.layerCount = 1;
    copy.imageExtent = {width_, height_, 1};
    vkCmdCopyImageToBuffer(commands, colour_.handle, VK_IMAGE_LAYOUT_TRANSFER_SRC_OPTIMAL,
                           readback_.handle, 1, &copy);

    if (!device_->endOneShot(commands)) return false;

    out = core::Image(width_, height_);
    void* mapped = nullptr;
    if (vkMapMemory(device_->handle(), readback_.memory, 0, out.rgba.size(), 0, &mapped) !=
        VK_SUCCESS)
        return false;
    std::memcpy(out.rgba.data(), mapped, out.rgba.size());
    vkUnmapMemory(device_->handle(), readback_.memory);
    return true;
}

}  // namespace gpu
