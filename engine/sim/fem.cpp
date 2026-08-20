// SPDX-License-Identifier: MIT
#include "fem.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace sim::fem {
namespace {

// Column-major 3x3 helpers, matching GLSL's mat3 so the CPU reference and the
// compute kernel can be compared term by term.
struct M3 {
    float m[9];  // m[col * 3 + row]
    float& operator()(int r, int c) { return m[c * 3 + r]; }
    float  operator()(int r, int c) const { return m[c * 3 + r]; }
};

M3 multiply(const M3& a, const M3& b) {
    M3 r{};
    for (int c = 0; c < 3; ++c)
        for (int i = 0; i < 3; ++i) {
            float s = 0;
            for (int k = 0; k < 3; ++k) s += a(i, k) * b(k, c);
            r(i, c) = s;
        }
    return r;
}

M3 transpose(const M3& a) {
    M3 r{};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) r(i, j) = a(j, i);
    return r;
}

float determinant(const M3& a) {
    return a(0, 0) * (a(1, 1) * a(2, 2) - a(1, 2) * a(2, 1)) -
           a(0, 1) * (a(1, 0) * a(2, 2) - a(1, 2) * a(2, 0)) +
           a(0, 2) * (a(1, 0) * a(2, 1) - a(1, 1) * a(2, 0));
}

// Inverse, or the zero matrix when the argument is too near singular to invert.
//
// **The threshold is a shape ratio, not a size, because the two callers below
// differ by the units of a volume.** `polarRotation` hands this a near-rotation:
// nondimensional, det ~ 1. `computeRestState` hands it `dm`, whose columns are
// node position *differences* in metres, so det(dm) = 6 * restVolume and the
// units are m^3. That is the defect in the 1e-30f that used to stand here -- not
// that it is denormal-ish for float, which is true but incidental, but that it
// is an absolute number being asked to bound a quantity that scales as h^3.
//
// For an element of edge h, det ~ h^3; at the h = 0.1 m of a plate mesh that is
// 1e-3 m^3, twenty-seven decades above the old floor. A *fully flat or inverted*
// tet -- precisely what this guard exists to catch -- leaves |det| ~ eps_float *
// h^3 ~ 1.2e-10 m^3, still twenty decades clear of it. It passed: `restInv` came
// back with entries ~1e10 m^-1 where healthy is ~1/h ~ 10, F = Ds * restInv was
// amplified ~1e9, and restVolume ~ 0 meant the element had no mass to resist it.
//
// `scale` is the product of the three column lengths, so it carries the same m^3
// the determinant does and the ratio is nondimensional. Hadamard's inequality
// makes |det| <= scale for every matrix, with equality exactly when the columns
// are orthogonal, so the ratio is a shape quality in [0, 1]: 1 for a perfect
// corner, 0 for a flat one. That is what makes one form serve both callers --
// polarRotation's near-unit columns give scale ~ 1 against det ~ 1.
//
// The L1 column norm is used rather than the Euclidean one to keep three square
// roots out of `tetForces`, which runs this four times per element per step.
// L1 >= L2 >= L2-product, so the ratio stays bounded by 1 and the guard stays
// conservative; the penalty is at worst 3^(3/2) = 5.2, and the constant has six
// decades of room for it. Measured: a general rotation gives 0.25, the Kuhn
// corner tet `makeBoxTetMesh` actually emits gives 0.083, and the near-flat tet
// gives 1.7e-9 -- so 1e-6f sits five decades below every healthy case and nearly
// three above the degenerate one.
//
// **The alternative -- `solid_shell.cpp`'s `invert3`, which tests `det == 0.0`
// exactly and leaves the magnitude question to callers that assert `det > 0.0`
// -- is right there and is not right here.** That works because inversion of an
// element Jacobian is a sign question, and solid_shell's callers ask the sign
// question separately. `computeRestState` asks nobody: it takes the inverse and
// moves on. An exact-zero test would let the 1.2e-10 m^3 flat tet through
// untouched, which is the whole failure. Written `>` so a NaN refuses too.
M3 inverse(const M3& a) {
    const float det = determinant(a);
    float scale = 1.0f;
    for (int c = 0; c < 3; ++c)
        scale *= std::abs(a(0, c)) + std::abs(a(1, c)) + std::abs(a(2, c));
    const float id = std::abs(det) > 1e-6f * scale ? 1.0f / det : 0.0f;
    M3 r{};
    r(0, 0) = (a(1, 1) * a(2, 2) - a(1, 2) * a(2, 1)) * id;
    r(0, 1) = (a(0, 2) * a(2, 1) - a(0, 1) * a(2, 2)) * id;
    r(0, 2) = (a(0, 1) * a(1, 2) - a(0, 2) * a(1, 1)) * id;
    r(1, 0) = (a(1, 2) * a(2, 0) - a(1, 0) * a(2, 2)) * id;
    r(1, 1) = (a(0, 0) * a(2, 2) - a(0, 2) * a(2, 0)) * id;
    r(1, 2) = (a(0, 2) * a(1, 0) - a(0, 0) * a(1, 2)) * id;
    r(2, 0) = (a(1, 0) * a(2, 1) - a(1, 1) * a(2, 0)) * id;
    r(2, 1) = (a(0, 1) * a(2, 0) - a(0, 0) * a(2, 1)) * id;
    r(2, 2) = (a(0, 0) * a(1, 1) - a(0, 1) * a(1, 0)) * id;
    return r;
}

// Rotation part of F, by Higham's Newton iteration R <- (R + R^-T)/2.
// Quadratic convergence, so four iterations reaches float precision. The kernel
// uses the same fixed count -- an explicit solver amplifies any difference, so
// the two implementations have to agree step for step, not just eventually.
M3 polarRotation(const M3& f, int iterations = 4) {
    M3 r = f;
    for (int i = 0; i < iterations; ++i) {
        const M3 rit = transpose(inverse(r));
        for (int k = 0; k < 9; ++k) r.m[k] = 0.5f * (r.m[k] + rit.m[k]);
    }
    return r;
}

}  // namespace

void tetForces(const TetMesh& mesh, std::size_t t, float lambda, float mu, float out[12]) {
    const uint32_t* idx = &mesh.index[t * 4];
    const float* p = mesh.position.data();

    auto node = [&](int corner, int axis) { return p[idx[corner] * 3 + axis]; };

    // Current shape matrix: edges from corner 0, as columns.
    M3 ds{};
    for (int c = 0; c < 3; ++c)
        for (int a = 0; a < 3; ++a) ds(a, c) = node(c + 1, a) - node(0, a);

    M3 dmInv{};
    for (int k = 0; k < 9; ++k) dmInv.m[k] = mesh.restInv[t * 9 + k];

    const M3 f = multiply(ds, dmInv);
    const M3 r = polarRotation(f);

    // Co-rotational small strain: unrotate, then use the linear strain tensor.
    // This is what lets a linear material survive large rotations without the
    // spurious volume growth that plain linear elasticity produces.
    const M3 s = multiply(transpose(r), f);
    M3 strain{};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) strain(i, j) = 0.5f * (s(i, j) + s(j, i));
    for (int i = 0; i < 3; ++i) strain(i, i) -= 1.0f;

    const float trace = strain(0, 0) + strain(1, 1) + strain(2, 2);
    M3 stress{};
    for (int k = 0; k < 9; ++k) stress.m[k] = 2.0f * mu * strain.m[k];
    for (int i = 0; i < 3; ++i) stress(i, i) += lambda * trace;

    // First Piola-Kirchhoff, rotated back into the current frame.
    const M3 piola = multiply(r, stress);
    const M3 h = multiply(piola, transpose(dmInv));

    const float scale = -mesh.restVolume[t];
    float f0[3] = {0, 0, 0};
    for (int c = 0; c < 3; ++c)
        for (int a = 0; a < 3; ++a) {
            const float v = scale * h(a, c);
            out[(c + 1) * 3 + a] = v;
            f0[a] -= v;
        }
    for (int a = 0; a < 3; ++a) out[a] = f0[a];
}

void TetMesh::computeRestState(const Material& material) {
    const std::size_t nodes = position.size() / 3;
    const std::size_t tets = index.size() / 4;

    velocity.assign(nodes * 3, 0.0f);
    mass.assign(nodes, 0.0f);
    if (fixed.size() != nodes) fixed.assign(nodes, 0u);
    restInv.resize(tets * 9);
    restVolume.resize(tets);
    tetForce.assign(tets * 12, 0.0f);

    for (std::size_t t = 0; t < tets; ++t) {
        const uint32_t* idx = &index[t * 4];
        M3 dm{};
        for (int c = 0; c < 3; ++c)
            for (int a = 0; a < 3; ++a)
                dm(a, c) = position[idx[c + 1] * 3 + a] - position[idx[0] * 3 + a];

        const float volume = determinant(dm) / 6.0f;
        restVolume[t] = std::abs(volume);
        const M3 inv = inverse(dm);
        for (int k = 0; k < 9; ++k) restInv[t * 9 + k] = inv.m[k];

        // Lumped mass: each tet gives a quarter of its mass to each corner.
        const float quarter = material.density * restVolume[t] * 0.25f;
        for (int c = 0; c < 4; ++c) mass[idx[c]] += quarter;
    }

    // CSR adjacency, node -> incident corners.
    adjacencyOffset.assign(nodes + 1, 0u);
    for (std::size_t t = 0; t < tets; ++t)
        for (int c = 0; c < 4; ++c) ++adjacencyOffset[index[t * 4 + c] + 1];
    for (std::size_t n = 0; n < nodes; ++n) adjacencyOffset[n + 1] += adjacencyOffset[n];

    adjacencyEntry.resize(adjacencyOffset.back());
    std::vector<uint32_t> cursor(adjacencyOffset.begin(), adjacencyOffset.end() - 1);
    for (std::size_t t = 0; t < tets; ++t)
        for (uint32_t c = 0; c < 4; ++c)
            adjacencyEntry[cursor[index[t * 4 + c]]++] = static_cast<uint32_t>(t) * 4 + c;
}

// **A minimum fold cannot carry a NaN, in either argument order.** This mattered
// because a NaN node used to vanish here and `criticalTimestep` then returned an
// entirely plausible stable step for geometry that cannot be integrated at all --
// the shortest edge of the mesh's *finite* tets, with nothing anywhere saying the
// rest had been dropped.
//
// The obvious repair is wrong, and it was tried first. `std::min(a, b)` is
// `b < a ? b : a`, so `std::min(shortest, x)` returns `shortest` for a NaN `x`
// (dropped) while `std::min(x, shortest)` returns `x` (carried) -- which looks like
// the argument order fixes it. It does not: once the accumulator *is* a NaN, the
// very next finite edge evaluates `std::min(finite, NaN)` as
// `NaN < finite ? NaN : finite` and overwrites it. A NaN survives only if it
// happens to be the last edge examined, which is a lottery rather than a guard.
//
// So the non-finite case is tested explicitly and kept out of the fold entirely.
float TetMesh::shortestEdge() const {
    float shortest = std::numeric_limits<float>::infinity();
    bool nonFinite = false;
    for (std::size_t t = 0; t < tetCount(); ++t) {
        const uint32_t* idx = &index[t * 4];
        for (int a = 0; a < 4; ++a)
            for (int b = a + 1; b < 4; ++b) {
                float d2 = 0;
                for (int k = 0; k < 3; ++k) {
                    const float d = position[idx[a] * 3 + k] - position[idx[b] * 3 + k];
                    d2 += d * d;
                }
                const float edge = std::sqrt(d2);
                if (!std::isfinite(edge))
                    nonFinite = true;
                else
                    shortest = std::min(shortest, edge);
            }
    }
    // An empty mesh keeps the `+infinity` it started with, which is what
    // `criticalTimestep` turns into its zero sentinel. A mesh with a non-finite
    // node is a different thing and says so.
    return nonFinite ? std::numeric_limits<float>::quiet_NaN() : shortest;
}

float TetMesh::totalMass() const {
    float m = 0;
    for (float v : mass) m += v;
    return m;
}

float criticalTimestep(const TetMesh& mesh, const Material& material, float safety) {
    // **Zero when there is no answer, which is what the hex solver's identical
    // function has always done.** `solidshell::criticalTimestep` ends
    // `return std::isfinite(smallest) ? smallest : 0.0;`; this one ended with a
    // bare division, and there are three ways to reach it without an answer:
    //
    //   * An empty mesh never enters `shortestEdge`'s loop, so it returns
    //     `+infinity` and this returned `+infinity`. A caller writing
    //     `dt = std::min(dt, criticalTimestep(...))` keeps its own `dt` unchanged
    //     and never learns the mesh was empty -- the step it runs at is a step
    //     nothing chose.
    //   * A NaN node coordinate used to be **swallowed** by `shortestEdge`'s
    //     minimum fold, so the shortest edge came back as the shortest of the
    //     *finite* tets and this returned an entirely plausible stable timestep for
    //     a mesh with NaN geometry. That is the bad one: the number looks right,
    //     `stepCpu` then produces NaN everywhere, and nothing points back here.
    //     `shortestEdge` now reports a NaN for that mesh -- see the note there for
    //     why no argument order of `std::min` was enough -- and the guard below
    //     turns it into the sentinel.
    //   * `waveSpeed()` is `sqrt(E / rho)` with no guard on either. Zero density
    //     gives `+infinity` and a step of 0; zero modulus gives 0 and a step of
    //     `+infinity`; both zero gives NaN.
    //
    // Zero is the sentinel because it is the one every consumer already tests:
    // `!(dt > 0)` catches it, and catches NaN with it, while no comparison against
    // an infinity is false.
    const float dt = safety * mesh.shortestEdge() / material.waveSpeed();
    return std::isfinite(dt) && dt > 0.0f ? dt : 0.0f;
}

void stepCpu(TetMesh& mesh, const Material& material, float dt, float gravity, float damping) {
    const float lambda = material.lameLambda();
    const float mu = material.lameMu();
    const std::size_t tets = mesh.tetCount();
    const std::size_t nodes = mesh.nodeCount();

    for (std::size_t t = 0; t < tets; ++t)
        tetForces(mesh, t, lambda, mu, &mesh.tetForce[t * 12]);

    for (std::size_t n = 0; n < nodes; ++n) {
        if (mesh.fixed[n]) {
            for (int k = 0; k < 3; ++k) mesh.velocity[n * 3 + k] = 0.0f;
            continue;
        }
        float force[3] = {0.0f, 0.0f, -gravity * mesh.mass[n]};
        for (uint32_t e = mesh.adjacencyOffset[n]; e < mesh.adjacencyOffset[n + 1]; ++e) {
            const uint32_t entry = mesh.adjacencyEntry[e];
            const float* f = &mesh.tetForce[(entry >> 2) * 12 + (entry & 3u) * 3];
            for (int k = 0; k < 3; ++k) force[k] += f[k];
        }
        const float invMass = 1.0f / mesh.mass[n];
        for (int k = 0; k < 3; ++k) {
            float v = mesh.velocity[n * 3 + k] + force[k] * invMass * dt;
            v *= damping;
            mesh.velocity[n * 3 + k] = v;
            mesh.position[n * 3 + k] += v * dt;
        }
    }
}

TetMesh makeBoxTetMesh(float lx, float ly, float lz, int nx, int ny, int nz) {
    TetMesh mesh;
    const int sx = nx + 1, sy = ny + 1, sz = nz + 1;
    const auto nodeAt = [&](int i, int j, int k) {
        return static_cast<uint32_t>((k * sy + j) * sx + i);
    };

    mesh.position.resize(static_cast<std::size_t>(sx) * sy * sz * 3);
    for (int k = 0; k < sz; ++k)
        for (int j = 0; j < sy; ++j)
            for (int i = 0; i < sx; ++i) {
                float* p = &mesh.position[nodeAt(i, j, k) * 3];
                p[0] = lx * static_cast<float>(i) / static_cast<float>(nx);
                p[1] = ly * static_cast<float>(j) / static_cast<float>(ny);
                p[2] = lz * static_cast<float>(k) / static_cast<float>(nz);
            }

    // Kuhn subdivision: six tets per cell, all sharing the cell's main diagonal
    // 0-7. Adjacent cells produce matching faces, so the mesh is conforming.
    static const int kCellTets[6][4] = {{0, 1, 3, 7}, {0, 1, 7, 5}, {0, 5, 7, 4},
                                        {0, 3, 2, 7}, {0, 2, 6, 7}, {0, 4, 6, 7}};

    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i) {
                uint32_t corner[8];
                for (int c = 0; c < 8; ++c)
                    corner[c] = nodeAt(i + (c & 1), j + ((c >> 1) & 1), k + ((c >> 2) & 1));

                for (const auto& tet : kCellTets) {
                    uint32_t v[4];
                    for (int c = 0; c < 4; ++c) v[c] = corner[tet[c]];

                    // Enforce positive orientation so rest volumes and the
                    // deformation gradient come out with a consistent sign.
                    M3 dm{};
                    for (int c = 0; c < 3; ++c)
                        for (int a = 0; a < 3; ++a)
                            dm(a, c) = mesh.position[v[c + 1] * 3 + a] - mesh.position[v[0] * 3 + a];
                    if (determinant(dm) < 0) std::swap(v[1], v[2]);

                    mesh.index.insert(mesh.index.end(), {v[0], v[1], v[2], v[3]});
                }
            }

    return mesh;
}

}  // namespace sim::fem
