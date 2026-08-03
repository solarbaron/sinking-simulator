// SPDX-License-Identifier: MIT
#include "fem.hpp"

#include <algorithm>
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

M3 inverse(const M3& a) {
    const float det = determinant(a);
    const float id = std::abs(det) > 1e-30f ? 1.0f / det : 0.0f;
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

float TetMesh::shortestEdge() const {
    float shortest = std::numeric_limits<float>::infinity();
    for (std::size_t t = 0; t < tetCount(); ++t) {
        const uint32_t* idx = &index[t * 4];
        for (int a = 0; a < 4; ++a)
            for (int b = a + 1; b < 4; ++b) {
                float d2 = 0;
                for (int k = 0; k < 3; ++k) {
                    const float d = position[idx[a] * 3 + k] - position[idx[b] * 3 + k];
                    d2 += d * d;
                }
                shortest = std::min(shortest, std::sqrt(d2));
            }
    }
    return shortest;
}

float TetMesh::totalMass() const {
    float m = 0;
    for (float v : mass) m += v;
    return m;
}

float criticalTimestep(const TetMesh& mesh, const Material& material, float safety) {
    return safety * mesh.shortestEdge() / material.waveSpeed();
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
