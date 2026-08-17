// SPDX-License-Identifier: MIT
//
// Explicit co-rotational tetrahedral FEM -- the Tier-2 solver from
// docs/02-simulation.md section 3, in its simplest honest form.
//
// Linear tets, lumped mass, symplectic-Euler integration, co-rotational
// linear elasticity. No plasticity or fracture yet: the question this answers
// is the throughput question, because that is what decides whether adaptive
// full-3D FEM on a ship is a plan or a wish.
//
// Storage is deliberately flat float arrays rather than vec3/mat3 structs. They
// upload to std430 storage buffers with no alignment translation, so the CPU
// reference and the GPU kernels operate on bit-identical layouts and can be
// compared directly.
#pragma once

#include <array>
#include <cmath>
#include <cstdint>
#include <vector>

namespace sim::fem {

struct Material {
    float youngsModulus = 210.0e9f;  // Pa, mild steel
    float poissonRatio  = 0.30f;
    float density       = 7850.0f;   // kg/m^3

    float lameLambda() const {
        return youngsModulus * poissonRatio /
               ((1.0f + poissonRatio) * (1.0f - 2.0f * poissonRatio));
    }
    float lameMu() const { return youngsModulus / (2.0f * (1.0f + poissonRatio)); }
    // Dilatational wave speed; sets the explicit stability limit.
    float waveSpeed() const { return std::sqrt(youngsModulus / density); }
};

struct TetMesh {
    // Per node, 3 floats each.
    std::vector<float>    position;
    std::vector<float>    velocity;
    std::vector<float>    mass;
    std::vector<uint32_t> fixed;    // non-zero pins the node

    // Per tet.
    std::vector<uint32_t> index;      // 4 node indices
    std::vector<float>    restInv;    // 9 floats: inverse of the rest shape matrix
    std::vector<float>    restVolume; // 1 float

    // Scratch: 12 floats per tet, the four nodal force vectors it contributes.
    std::vector<float> tetForce;

    // Node -> incident corners, as CSR. Gathering rather than scattering keeps
    // the accumulation order fixed, which is what makes the solver deterministic
    // and avoids needing float atomics at all.
    std::vector<uint32_t> adjacencyOffset;  // nodeCount + 1
    std::vector<uint32_t> adjacencyEntry;   // tetIndex * 4 + corner

    std::size_t nodeCount() const { return mass.size(); }
    std::size_t tetCount()  const { return restVolume.size(); }

    // Fills restInv, restVolume, lumped nodal mass and the adjacency, from the
    // current positions taken as the rest configuration.
    void computeRestState(const Material& material);

    float shortestEdge() const;
    float totalMass() const;
};

// A box of lx by ly by lz split into a regular grid of cells, each cell into six
// tetrahedra (the Kuhn subdivision, which tiles space consistently so shared
// faces always match).
TetMesh makeBoxTetMesh(float lx, float ly, float lz, int nx, int ny, int nz);

// Largest stable explicit step: the time for a dilatational wave to cross the
// smallest element, with the usual safety factor.
//
// **Zero means there is no answer** -- an empty mesh, a non-finite node, or a
// material with no wave speed -- and is the same sentinel
// `solidshell::criticalTimestep` uses. Test it with `!(dt > 0)`, which catches a
// NaN too; a comparison against an infinity is never false, which is what this
// used to return for an empty mesh.
float criticalTimestep(const TetMesh& mesh, const Material& material, float safety = 0.5f);

// One explicit step. `damping` multiplies velocity each step (1.0 = undamped);
// it stands in for the mass-proportional damping a real solver would use, and is
// how the cantilever test settles to its static deflection.
void stepCpu(TetMesh& mesh, const Material& material, float dt, float gravity, float damping);

// Exposed for the GPU kernel to be checked against: the four nodal forces a
// single tet contributes, given the current node positions.
void tetForces(const TetMesh& mesh, std::size_t tet, float lambda, float mu, float out[12]);

}  // namespace sim::fem
