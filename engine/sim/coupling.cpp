// SPDX-License-Identifier: MIT
#include "coupling.hpp"

#include <algorithm>
#include <cmath>

namespace sim::coupling {

using solidshell::kNodes;

// --- 1. A boundary condition on a reduced model ---------------------------------

bool prescribedStaticSolve(const reduction::Assembly& assembly, const std::vector<double>& load,
                           const std::vector<std::uint32_t>& held,
                           const std::vector<Prescribed>& prescribed, std::vector<double>& state,
                           std::string* problem) {
    return reduction::assembledStaticSolve(assembly, load, held, prescribed, state, problem);
}

// --- 2. Matching a zone to the structure round it --------------------------------

Coupling couple(const reduction::Substructure& surroundings,
                const reduction::Reduction& reducedSurroundings,
                const reduction::Substructure& zone, const reduction::Reduction& reducedZone,
                double tolerance) {
    Coupling out;
    out.map = reduction::matchBoundaries(surroundings, zone, tolerance);
    out.worstGap = out.map.worstGap;
    out.sharedDof = out.map.shared;

    // The map is expressed in the *substructures'* boundary DOF and the assembly in
    // the *reductions'*. They agree by construction, and a reduction that fell back
    // to an empty one -- which `craigBampton` does rather than throwing -- is the
    // case where they do not.
    if (reducedSurroundings.boundary != static_cast<int>(surroundings.boundaryCount()) ||
        reducedZone.boundary != static_cast<int>(zone.boundaryCount())) {
        out.problems.push_back(
            "a reduction does not carry its substructure's boundary, so the interface map does "
            "not describe it");
        return out;
    }

    out.assembly = reduction::assemble(reducedSurroundings, reducedZone, out.map);
    for (const std::string& p : out.assembly.problems) out.problems.push_back(p);
    if (out.assembly.empty()) return out;

    out.surroundDof.assign(static_cast<std::size_t>(reducedSurroundings.boundary), -1);
    for (int i = 0; i < reducedSurroundings.boundary; ++i)
        out.surroundDof[static_cast<std::size_t>(i)] =
            out.assembly.fromA()[static_cast<std::size_t>(i)];
    out.zoneDof.assign(static_cast<std::size_t>(reducedZone.boundary), -1);
    for (int j = 0; j < reducedZone.boundary; ++j)
        out.zoneDof[static_cast<std::size_t>(j)] = out.assembly.fromB()[static_cast<std::size_t>(j)];

    out.zoneShared.assign(static_cast<std::size_t>(reducedZone.boundary), 0u);
    for (int j : out.map.aToB)
        if (j >= 0 && j < reducedZone.boundary) out.zoneShared[static_cast<std::size_t>(j)] = 1u;
    out.zoneUnshared = 0;
    for (std::uint8_t shared : out.zoneShared)
        if (!shared) ++out.zoneUnshared;
    return out;
}

// --- 3. Driving the zone's edge ---------------------------------------------------

EdgeDrive edgeDrive(const Coupling& coupling, const reduction::Substructure& zone,
                    const std::vector<double>& assembledState) {
    EdgeDrive out;
    out.displacement.assign(zone.dofCount(), 0.0);
    out.driven.assign(zone.dofCount(), 0u);
    if (!coupling.ready()) {
        out.problems.push_back("the coupling is not ready, so there is no interface to drive from");
        return out;
    }
    if (assembledState.size() != static_cast<std::size_t>(coupling.assembly.size())) {
        // The same failure `componentState` guards: every index is in range, so a
        // state of the wrong length comes back as a plausible field rather than as
        // an error.
        out.problems.push_back("the state is not this assembly's, so its interface DOF are not "
                               "the ones being read");
        return out;
    }
    const std::vector<std::uint32_t>& boundary = zone.boundaryDof();
    if (boundary.size() != coupling.zoneDof.size()) {
        out.problems.push_back("the coupling was not built against this zone");
        return out;
    }

    for (std::size_t j = 0; j < boundary.size(); ++j) {
        if (!coupling.zoneShared[j]) continue;
        const int a = coupling.zoneDof[j];
        if (a < 0) continue;
        const std::size_t g = boundary[j];
        if (g >= out.displacement.size()) continue;
        out.displacement[g] = assembledState[static_cast<std::size_t>(a)];
        out.driven[g] = 1u;
        ++out.count;
        out.largest = std::max(out.largest, std::fabs(out.displacement[g]));
    }
    if (out.count == 0)
        out.problems.push_back("nothing was driven: the zone shares no boundary DOF");
    return out;
}

std::size_t applyEdgeDrive(const EdgeDrive& drive, solidshell::HexMesh& mesh) {
    const std::size_t dofs = mesh.nodeCount() * 3;
    if (mesh.fixed.size() < dofs) mesh.fixed.resize(dofs, 0u);
    if (mesh.prescribed.size() < dofs) mesh.prescribed.resize(dofs, 0.0);

    std::size_t written = 0;
    const std::size_t limit = std::min(dofs, drive.driven.size());
    for (std::size_t d = 0; d < limit; ++d) {
        if (!drive.driven[d]) continue;
        mesh.fixed[d] = 1u;
        mesh.prescribed[d] = drive.displacement[d];
        ++written;
    }
    return written;
}

// --- 4. What a torn zone hands back ------------------------------------------------

DamagedMesh withoutElements(const solidshell::HexMesh& mesh,
                            const std::vector<std::uint8_t>& removed) {
    DamagedMesh out;
    const std::size_t elements = mesh.elementCount();
    const std::size_t nodes = mesh.nodeCount();

    std::vector<std::uint8_t> used(nodes, 0u);
    for (std::size_t e = 0; e < elements; ++e) {
        if (e < removed.size() && removed[e]) {
            ++out.removedElements;
            continue;
        }
        out.element.push_back(static_cast<std::uint32_t>(e));
        for (int a = 0; a < kNodes; ++a) {
            const std::uint32_t n = mesh.index[e * kNodes + static_cast<std::size_t>(a)];
            if (n < nodes) used[n] = 1u;
        }
    }

    // Ascending in the original numbering, which keeps the node ordering the mesh
    // was built with. `makePlateMesh`'s ordering is what holds the assembled
    // bandwidth down (`solid_shell.hpp`), and a compaction that reordered would
    // throw that away for nothing.
    std::vector<std::uint32_t> slot(nodes, 0u);
    for (std::size_t n = 0; n < nodes; ++n) {
        if (!used[n]) {
            ++out.orphanedNodes;
            continue;
        }
        slot[n] = static_cast<std::uint32_t>(out.node.size());
        out.node.push_back(static_cast<std::uint32_t>(n));
    }

    out.mesh.position.reserve(out.node.size() * 3);
    out.mesh.fixed.reserve(out.node.size() * 3);
    out.mesh.prescribed.reserve(out.node.size() * 3);
    for (std::uint32_t n : out.node)
        for (int axis = 0; axis < 3; ++axis) {
            const std::size_t d = static_cast<std::size_t>(n) * 3 + static_cast<std::size_t>(axis);
            out.mesh.position.push_back(mesh.position[d]);
            out.mesh.fixed.push_back(d < mesh.fixed.size() ? mesh.fixed[d] : std::uint8_t{0});
            out.mesh.prescribed.push_back(d < mesh.prescribed.size() ? mesh.prescribed[d] : 0.0);
        }

    // Bounds-checked on the same rule as the pass that filled `used` above: a
    // malformed index would otherwise be caught there and read past `slot` here,
    // which is a worse failure than the one being guarded against.
    out.mesh.index.reserve(out.element.size() * kNodes);
    for (std::uint32_t e : out.element)
        for (int a = 0; a < kNodes; ++a) {
            const std::uint32_t n =
                mesh.index[static_cast<std::size_t>(e) * kNodes + static_cast<std::size_t>(a)];
            out.mesh.index.push_back(n < nodes ? slot[n] : 0u);
        }

    if (out.element.empty()) out.problems.push_back("every element was removed: nothing is left");
    if (out.orphanedNodes > 0)
        out.problems.push_back(std::to_string(out.orphanedNodes) +
                               " node(s) were left with no element and have been dropped");
    return out;
}

DamagedMesh withoutTornElements(const zone::Patch& patch, const zone::Solver& solver) {
    const std::vector<solidshell::ElementPlasticState>& state = solver.elementState();
    std::vector<std::uint8_t> removed(patch.elementCount(), 0u);
    if (state.empty()) {
        // The elastic path keeps no plastic state at all, so nothing can have torn.
        // Saying so is the difference between "no damage" and "no damage model",
        // and a caller reading an empty reduction back into Tier 1 would not
        // otherwise be able to tell.
        DamagedMesh out = withoutElements(patch.mesh, removed);
        out.problems.push_back("the zone was solved elastically, so no element can have torn");
        return out;
    }
    for (std::size_t e = 0; e < removed.size() && e < state.size(); ++e)
        removed[e] = state[e].torn ? std::uint8_t{1} : std::uint8_t{0};
    return withoutElements(patch.mesh, removed);
}

std::vector<std::uint32_t> carriedInterface(const DamagedMesh& damaged,
                                            const std::vector<std::uint32_t>& interfaceNodes) {
    std::vector<std::uint32_t> out;
    if (damaged.node.empty()) return out;
    // `damaged.node` is ascending, so the lookup is a binary search rather than a
    // reverse table sized by the original mesh.
    for (std::uint32_t n : interfaceNodes) {
        const auto found = std::lower_bound(damaged.node.begin(), damaged.node.end(), n);
        if (found == damaged.node.end() || *found != n) continue;
        out.push_back(static_cast<std::uint32_t>(found - damaged.node.begin()));
    }
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
}

}  // namespace sim::coupling
