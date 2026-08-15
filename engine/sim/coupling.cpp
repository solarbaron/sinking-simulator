// SPDX-License-Identifier: MIT
#include "coupling.hpp"

#include "plasticity.hpp"

#include <algorithm>
#include <cmath>
#include <string>

namespace sim::coupling {

using solidshell::kDof;
using solidshell::kGauss;
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
    // **No path through `couple()` is known to reach this, and it is kept anyway.**
    // `ready()` above already requires `sharedDof > 0`, and `sharedDof` counts the
    // matched entries that set `zoneShared[j] = 1` twenty lines up -- so a ready
    // coupling has at least one shared DOF, and the `a < 0` and `g >= size()`
    // continues above are guarded by invariants `assemble` and `Substructure`
    // maintain. The physically meaningful case this message describes -- a zone
    // sharing no boundary DOF -- is caught earlier by `ready()`, whose message is
    // the more accurate one.
    //
    // It stays because `Coupling` is an aggregate with public fields and no
    // invariant of its own: anything that builds one by hand, or a future `couple`
    // that fills the arrays partially, lands here rather than driving nothing in
    // silence. What it is *not* is a tested guard, and it should not be mistaken for
    // one -- the suite reaches the three above it and not this. Its negative control
    // is incidental: the happy-path tests assert `problems.empty()`, so this firing
    // would turn them red.
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

// --- 5. What a yielded zone hands back ---------------------------------------------
namespace {

// The shear modulus one element has left, as a fraction of the elastic one.
//
// **The average is taken over the plastic *compliance*, not over the modulus.**
// The closed form in `plasticity.hpp` is additive in compliance -- 1/G_s = 1/G +
// 3 eps_p / sigma_y -- so averaging the added term is averaging the thing that is
// linear in the damage, and it is also the softer of the two averages, which is
// the direction §5 exists to move in. Averaging the modulus instead would report a
// patchily yielded element stiffer than a uniformly yielded one carrying the same
// total flow.
//
// **An element none of whose points has flowed returns exactly 1.0, by the
// short-circuit rather than by the arithmetic.** An elastic point adds no
// compliance under either modulus, so it is skipped before either closed form is
// called; the alternative -- letting it add `1/G - 1/G` -- is not the same thing,
// because the tangent of a point that is not flowing is `G` and is not the limit of
// `tangentShearModulus` at zero. Exactness here is what makes the unyielded case
// bit-identical downstream rather than merely close, and mutation testing kills
// both the removal of this line and the removal of its partner in
// `plasticity::secantShearModulus`.
double shearRatio(const plasticity::Material& material,
                  const solidshell::ElementPlasticState& state, const double volume[kGauss],
                  Modulus modulus, double* peakPlasticStrain) {
    const double shear = material.shearModulus();
    double totalVolume = 0.0, plasticCompliance = 0.0;
    bool dead = false;
    for (int g = 0; g < kGauss; ++g) {
        const double plastic = state.point[g].equivalentPlasticStrain;
        if (peakPlasticStrain) *peakPlasticStrain = std::max(*peakPlasticStrain, plastic);
        totalVolume += volume[g];
        if (plastic <= 0.0) continue;  // elastic: no added compliance, either modulus
        const double left = modulus == Modulus::Secant
                                ? plasticity::secantShearModulus(material, plastic)
                                : plasticity::tangentShearModulus(material, plastic);
        if (!(left > 0.0)) {
            dead = true;  // a perfectly plastic tangent: infinite added compliance
            continue;
        }
        plasticCompliance += volume[g] * (1.0 / left - 1.0 / shear);
    }
    if (dead) return 0.0;
    if (!(totalVolume > 0.0) || plasticCompliance <= 0.0) return 1.0;
    return 1.0 / (1.0 + shear * plasticCompliance / totalVolume);
}

}  // namespace

Softening softening(const solidshell::HexMesh& mesh, const StructuralMaterial& elastic,
                    const plasticity::Material& material,
                    const std::vector<solidshell::ElementPlasticState>& state, Modulus modulus,
                    solidshell::Formulation form) {
    Softening out;
    const std::size_t elements = mesh.elementCount();
    if (state.size() != elements) {
        // Refused rather than read short. A state array of the wrong length is a
        // zone whose history belongs to a different mesh, and the blocks it would
        // produce would soften whichever elements happened to line up -- a
        // plausible field and the wrong one, which is the failure mode §5 is here
        // to remove rather than to introduce.
        out.problems.push_back("the plastic state is " + std::to_string(state.size()) +
                               " long for a mesh of " + std::to_string(elements) +
                               " elements, so it is not this mesh's history");
        return out;
    }
    if (elements == 0) return out;

    const double bulk = elastic.youngsModulus / (3.0 * (1.0 - 2.0 * elastic.poissonRatio));
    const double shear = elastic.youngsModulus / (2.0 * (1.0 + elastic.poissonRatio));

    double softVolume = 0.0, totalVolume = 0.0;
    for (std::size_t e = 0; e < elements; ++e) {
        double nodes[kDof], volume[kGauss];
        mesh.gather(e, mesh.position, nodes);
        solidshell::gaussVolumes(nodes, volume);
        double elementVolume = 0.0;
        for (int g = 0; g < kGauss; ++g) elementVolume += volume[g];
        totalVolume += elementVolume;

        if (state[e].torn) {
            ++out.torn;
            for (int g = 0; g < kGauss; ++g)
                out.peakPlasticStrain =
                    std::max(out.peakPlasticStrain, state[e].point[g].equivalentPlasticStrain);
            continue;  // §5: a tear is `withoutTornElements`, not a knockdown
        }

        const double ratio = shearRatio(material, state[e], volume, modulus, &out.peakPlasticStrain);
        softVolume += elementVolume * ratio;
        if (ratio >= 1.0) continue;  // never flowed: no block at all, so no rounding
        out.worstRatio = std::min(out.worstRatio, ratio);
        if (!(ratio > 0.0))
            out.problems.push_back(
                "element " + std::to_string(e) +
                " has no shear stiffness left, so the reduced zone has a hole in it that the "
                "interior factorisation will not reliably catch");

        double softened[kDof * kDof], intact[kDof * kDof];
        StructuralMaterial weakened = elastic;
        plasticity::isotropicFromBulkShear(bulk, ratio * shear, &weakened.youngsModulus,
                                           &weakened.poissonRatio);
        solidshell::elementStiffness(nodes, weakened, form, softened);
        solidshell::elementStiffness(nodes, elastic, form, intact);

        solidshell::DofBlock block;
        block.dof.resize(kDof);
        for (int a = 0; a < kNodes; ++a) {
            const std::uint32_t n = mesh.index[e * kNodes + static_cast<std::size_t>(a)];
            for (int k = 0; k < 3; ++k)
                block.dof[static_cast<std::size_t>(a * 3 + k)] = n * 3 + static_cast<std::uint32_t>(k);
        }
        // Symmetrised on the way out. `elementStiffness` is symmetric to rounding
        // rather than to the last bit, and `Substructure` scatters both triangles of
        // a block, so an asymmetry here would survive into a matrix the reduction
        // assumes is symmetric -- unlike the element assembly beside it, which never
        // sees the two halves apart.
        block.stiffness.resize(static_cast<std::size_t>(kDof) * kDof);
        for (int i = 0; i < kDof; ++i)
            for (int j = 0; j <= i; ++j) {
                const double d = 0.5 * ((softened[i * kDof + j] - intact[i * kDof + j]) +
                                        (softened[j * kDof + i] - intact[j * kDof + i]));
                block.stiffness[static_cast<std::size_t>(i) * kDof + static_cast<std::size_t>(j)] = d;
                block.stiffness[static_cast<std::size_t>(j) * kDof + static_cast<std::size_t>(i)] = d;
            }
        out.attachment.stiffness.push_back(std::move(block));
        out.attachment.stress.emplace_back(static_cast<std::size_t>(kDof), 0.0);
        out.element.push_back(static_cast<std::uint32_t>(e));
        ++out.softened;
    }

    if (totalVolume > 0.0) out.meanRatio = softVolume / totalVolume;
    if (!out.attachment.stiffness.empty()) out.attachment.mass.assign(mesh.nodeCount(), 0.0);
    if (out.torn > 0)
        out.problems.push_back(
            std::to_string(out.torn) +
            " torn element(s) were not softened: a tear is element deletion, not a knockdown -- "
            "see `withoutTornElements`");
    return out;
}

Softening softening(const zone::Patch& patch, const zone::Solver& solver,
                    const plasticity::Material& material, Modulus modulus) {
    if (solver.elementState().empty()) {
        // The same distinction `withoutTornElements` draws: "no damage" and "no
        // damage model" are different answers and a caller cannot tell them apart
        // from an empty result.
        Softening out;
        out.problems.push_back("the zone was solved elastically, so no element can have yielded");
        return out;
    }
    return softening(patch.mesh, patch.material, material, solver.elementState(), modulus);
}

}  // namespace sim::coupling
