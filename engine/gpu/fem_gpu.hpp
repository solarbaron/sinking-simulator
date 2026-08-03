// SPDX-License-Identifier: MIT
//
// Vulkan compute back-end for the explicit tet FEM.
//
// Scoped to the spike: it owns its own instance and device rather than sitting on
// a general renderer, because the question being answered ("how many element
// updates per second does Pascal give us?") does not need one, and building the
// general thing first is how the answer gets postponed.
#pragma once

#include "../sim/fem.hpp"

#include <string>

namespace gpu {

class FemGpuSolver {
public:
    ~FemGpuSolver();

    // Uploads the mesh and builds both pipelines. Returns false with a reason on
    // any failure, including "no Vulkan device" -- this must degrade, not abort.
    bool initialise(const sim::fem::TetMesh& mesh, const std::string& shaderDirectory,
                    std::string& error);

    // Records and submits `substeps` integration steps as one command buffer, and
    // returns the GPU-measured wall time in milliseconds. Batching matters: a
    // single explicit step at these element sizes is far shorter than a
    // submission round trip, so per-submit overhead would otherwise dominate the
    // measurement and hide the actual kernel cost.
    double run(int substeps, const sim::fem::Material& material, float dt, float gravity,
               float damping);

    // Copies node positions and velocities back into the mesh.
    void readback(sim::fem::TetMesh& mesh);

    const std::string& deviceName() const { return deviceName_; }
    bool valid() const { return device_ != nullptr; }

private:
    struct Impl;
    Impl* impl_ = nullptr;
    void* device_ = nullptr;  // mirrors impl_->device, for valid()
    std::string deviceName_;
};

}  // namespace gpu
