// SPDX-License-Identifier: MIT
//
// Validation of Tier-1 to Tier-2 coupling: a zone's edge driven by the reduced
// model round it, and a torn zone handed back to that model.
//
// **The reference is the monolithic plate**, meshed in one piece by
// `makePlateMesh` and solved in one piece by `solidshell::solveStatic`. It uses
// none of the reduction, none of the assembly and none of this file's splitting,
// which is the property that makes it a reference at all -- the same argument
// `test_reduction.cpp` makes for its free-free spectrum.
//
// The claim being tested is stronger than "close", and saying so is the whole
// design of the assertions here. `reduction.hpp` §1 property 1 is that static
// condensation is exact at the interface for *any* load at *any* mode count, and
// assembly is scatter-add, so a coupled zone does not approach the monolithic
// answer -- it **is** the monolithic answer, to the conditioning of two
// independent solves of the same system. Everything below is therefore asserted
// against a floor that was measured rather than against a tolerance that would
// have passed on a model which had lost the property and was merely well
// converged.
//
// Three guards, each because the obvious version of the test proves nothing:
//
//   * **The bracket has to be wide.** `zone.hpp` §4 says a clamped edge is too
//     stiff and a free one too soft. If the load died away inside the patch the
//     three answers would agree and an exactly-clamped zone would pass a test of
//     coupling. So the punch sits one element from the zone's edge and the
//     clamped-to-free spread is asserted to be large before anything is compared
//     to it.
//   * **The drive has to be non-trivial.** A coupling that computed a beautiful
//     interface displacement of zero is a clamped zone wearing a better name, so
//     the drive's own magnitude is checked against the field it came out of, and
//     the DOF it touches are counted.
//   * **The negative control has to fail.** The clamped zone is compared to the
//     monolithic answer with the same instrument, and has to *miss* it by orders
//     of magnitude more than the coupled one. Without that, the exactness
//     assertion could be measuring nothing but a well-converged mesh.
#include "engine/sim/constraint.hpp"
#include "engine/sim/coupling.hpp"
#include "engine/sim/plasticity.hpp"
#include "engine/sim/reduction.hpp"
#include "engine/sim/scantlings.hpp"
#include "engine/sim/solid_shell.hpp"
#include "engine/sim/zone.hpp"
#include "harness.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

using namespace sim;
using solidshell::HexMesh;
using solidshell::kDof;
using solidshell::kNodes;
using testing::expectEqual;
using testing::expectNear;
using testing::expectTrue;

namespace {

void expectEqualCount(const std::string& what, std::size_t got, std::size_t want) {
    testing::expectEqual(what, static_cast<long long>(got), static_cast<long long>(want));
}

// --- The structure ---------------------------------------------------------------
//
// A square plate, meshed three ways that have to agree: whole, as a Tier-2 patch
// over its middle, and as everything else. Square and flat because then the
// element extrusion is exact (`zone.hpp` §2) and nothing in the comparison is
// competing with a geometry error.

constexpr double kSide = 1.2;        // m
constexpr double kThickness = 0.012; // m
constexpr int kPanels = 4;           // panels a side in the StructuralMesh
constexpr int kSubdivision = 2;      // elements per panel edge -> 8 x 8 elements
constexpr int kElements = kPanels * kSubdivision;
constexpr double kPunchTravel = 2.0e-4;  // m the punch is pushed, along -z

const StructuralMaterial& steel() {
    static const StructuralMaterial material = ah36Steel();
    return material;
}

// The plate as a `StructuralMesh`, which is what `zone::buildPatch` meshes.
StructuralMesh flatPlate() {
    StructuralMesh mesh;
    mesh.materials.push_back(steel());
    mesh.frameSpacing = kSide;
    for (int i = 0; i < kPanels; ++i)
        for (int j = 0; j < kPanels; ++j) {
            PlatePanel p;
            const double x0 = -0.5 * kSide + kSide * i / kPanels;
            const double x1 = -0.5 * kSide + kSide * (i + 1) / kPanels;
            const double y0 = -0.5 * kSide + kSide * j / kPanels;
            const double y1 = -0.5 * kSide + kSide * (j + 1) / kPanels;
            p.corner[0] = {x0, y0, 0};
            p.corner[1] = {x1, y0, 0};
            p.corner[2] = {x1, y1, 0};
            p.corner[3] = {x0, y1, 0};
            p.thickness = kThickness;
            p.material = 0;
            p.role = PanelRole::Shell;
            mesh.panels.push_back(p);
        }
    return mesh;
}

// The same plate as one hex mesh, centred on the origin. This is the reference and
// it is built by the element library, not by anything here.
HexMesh wholePlate(double thickness = kThickness) {
    HexMesh mesh = solidshell::makePlateMesh(kSide, kSide, thickness, kElements, kElements, 1);
    for (std::size_t n = 0; n < mesh.nodeCount(); ++n) {
        mesh.position[3 * n] -= 0.5 * kSide;
        mesh.position[3 * n + 1] -= 0.5 * kSide;
    }
    return mesh;
}

zone::MeshParams patchParams() {
    zone::MeshParams params;
    // 0.3 m reaches the four middle panels, whose centroids are 0.212 m out, and
    // stops short of the next ring at 0.474 m. So the patch is exactly the middle
    // quarter of the plate and its perimeter is entirely interior -- which is what
    // makes every perimeter DOF something the surroundings can drive.
    params.radius = 0.3;
    params.subdivision = kSubdivision;
    params.outward = {0, 0, 1};
    return params;
}

std::size_t nodeAt(const HexMesh& mesh, double x, double y, double z, double tolerance = 1e-9) {
    for (std::size_t n = 0; n < mesh.nodeCount(); ++n)
        if (std::fabs(mesh.position[3 * n] - x) < tolerance &&
            std::fabs(mesh.position[3 * n + 1] - y) < tolerance &&
            std::fabs(mesh.position[3 * n + 2] - z) < tolerance)
            return n;
    return mesh.nodeCount();
}

void elementCentroid(const HexMesh& mesh, std::size_t element, double out[3]) {
    out[0] = out[1] = out[2] = 0.0;
    for (int a = 0; a < kNodes; ++a) {
        const std::size_t n = mesh.index[element * kNodes + static_cast<std::size_t>(a)];
        for (int k = 0; k < 3; ++k) out[k] += mesh.position[3 * n + static_cast<std::size_t>(k)];
    }
    for (int k = 0; k < 3; ++k) out[k] /= kNodes;
}

// --- Reading a reaction out of a solved field --------------------------------------
//
// f = K u, assembled from `solidshell::elementStiffness` the same way `solveStatic`
// assembled the system it solved. At a free DOF this is the residual, which has to
// be zero; at a prescribed DOF it is the reaction, which is what the structure
// pushes back with. Doing it here rather than asking the solver keeps the reaction
// an independent reading of the displacement field, and it also gives the
// equilibrium residual for free -- which is the check that says the reference
// solution really is one.
std::vector<double> nodalForce(const HexMesh& mesh, const StructuralMaterial& material,
                               const std::vector<double>& u) {
    std::vector<double> f(mesh.nodeCount() * 3, 0.0);
    for (std::size_t e = 0; e < mesh.elementCount(); ++e) {
        double nodePosition[kDof], displacement[kDof], ke[kDof * kDof];
        mesh.gather(e, mesh.position, nodePosition);
        mesh.gather(e, u, displacement);
        solidshell::elementStiffness(nodePosition, material, solidshell::Formulation::SolidShell,
                                     ke);
        for (int i = 0; i < kDof; ++i) {
            double sum = 0;
            for (int j = 0; j < kDof; ++j) sum += ke[i * kDof + j] * displacement[j];
            const std::size_t node =
                mesh.index[e * kNodes + static_cast<std::size_t>(i / 3)];
            f[node * 3 + static_cast<std::size_t>(i % 3)] += sum;
        }
    }
    return f;
}

struct Reaction {
    double force = 0;     // N resisting the punch, positive
    double energy = 0;    // J stored, from the same field
    double residual = 0;  // N left unbalanced on a free DOF
};

// `punchNodes` are indices into `mesh`; the punch is prescribed at (0, 0, -travel).
Reaction punchReaction(const HexMesh& mesh, const StructuralMaterial& material,
                       const std::vector<double>& u,
                       const std::vector<std::size_t>& punchNodes) {
    const std::vector<double> f = nodalForce(mesh, material, u);
    Reaction out;
    for (std::size_t n : punchNodes) out.force -= f[3 * n + 2];
    for (std::size_t d = 0; d < f.size(); ++d) {
        out.energy += 0.5 * u[d] * f[d];
        if (d < mesh.fixed.size() && mesh.fixed[d]) continue;
        out.residual = std::max(out.residual, std::fabs(f[d]));
    }
    return out;
}

// --- The fixture ------------------------------------------------------------------

struct Fixture {
    StructuralMesh structure;
    HexMesh whole;
    zone::Patch patch;              // the Tier-2 zone
    HexMesh surround;               // the plate minus the zone's elements

    std::vector<std::uint8_t> inZone;      // per whole-plate element
    std::vector<std::uint32_t> zoneInterface;    // patch nodes: perimeter + punch
    std::vector<std::uint32_t> zonePerimeter;    // patch nodes: perimeter alone
    std::vector<std::uint32_t> surroundInterface;  // surroundings nodes
    std::vector<std::uint32_t> surroundOuter;      // the clamped rim, in surroundings numbering

    std::vector<std::size_t> punchWhole, punchPatch;  // the same nodes, two numberings
    std::vector<std::string> problems;
};

// The punch: a two-by-two block of mid-surface points, one element in from the
// zone's own edge and deliberately **off centre**, so that a symmetry nobody
// intended cannot flatter the answer.
bool isPunchPoint(double x, double y) {
    const double h = kSide / kElements;  // 0.15 m
    return (std::fabs(x + h) < 1e-9 || std::fabs(x) < 1e-9) &&
           (std::fabs(y) < 1e-9 || std::fabs(y - h) < 1e-9);
}

bool onRim(double x, double y) {
    const double half = 0.5 * kSide;
    return std::fabs(std::fabs(x) - half) < 1e-9 || std::fabs(std::fabs(y) - half) < 1e-9;
}

Fixture buildFixture() {
    Fixture f;
    f.structure = flatPlate();
    f.whole = wholePlate();
    f.patch = zone::buildPatch(f.structure, {0, 0, 0}, patchParams());

    // Which whole-plate elements the patch covers, by centroid. The two meshes are
    // built by different code and share no numbering, which is exactly the
    // situation `matchBoundaries` exists for and the reason nothing here assumes an
    // ordering.
    f.inZone.assign(f.whole.elementCount(), 0u);
    for (std::size_t p = 0; p < f.patch.elementCount(); ++p) {
        double centre[3];
        elementCentroid(f.patch.mesh, p, centre);
        std::size_t best = f.whole.elementCount();
        for (std::size_t e = 0; e < f.whole.elementCount(); ++e) {
            double other[3];
            elementCentroid(f.whole, e, other);
            if (std::fabs(centre[0] - other[0]) < 1e-9 && std::fabs(centre[1] - other[1]) < 1e-9 &&
                std::fabs(centre[2] - other[2]) < 1e-9)
                best = e;
        }
        if (best == f.whole.elementCount())
            f.problems.push_back("a patch element has no counterpart in the whole plate");
        else
            f.inZone[best] = 1u;
    }

    const coupling::DamagedMesh cut = coupling::withoutElements(f.whole, f.inZone);
    f.surround = cut.mesh;

    // The zone's interface: its perimeter -- which `buildPatch` has already pinned,
    // so `nodesPinned` reads it straight off -- plus the punch, because a DOF that
    // something outside the reduced model drives has to *be* in the reduced model
    // (`coupling.hpp` §2).
    f.zonePerimeter = reduction::nodesPinned(f.patch.mesh);
    f.zoneInterface = f.zonePerimeter;
    for (std::size_t n = 0; n < f.patch.nodeCount(); ++n)
        if (isPunchPoint(f.patch.mesh.position[3 * n], f.patch.mesh.position[3 * n + 1])) {
            f.punchPatch.push_back(n);
            f.zoneInterface.push_back(static_cast<std::uint32_t>(n));
        }
    std::sort(f.zoneInterface.begin(), f.zoneInterface.end());
    f.zoneInterface.erase(std::unique(f.zoneInterface.begin(), f.zoneInterface.end()),
                          f.zoneInterface.end());

    // The surroundings' interface: everything it shares with the zone, plus the
    // plate's rim, because the rim is what holds the whole model up and a reduced
    // model can only be held at a boundary DOF.
    for (std::size_t n = 0; n < f.surround.nodeCount(); ++n) {
        const double x = f.surround.position[3 * n], y = f.surround.position[3 * n + 1];
        const double z = f.surround.position[3 * n + 2];
        if (onRim(x, y)) {
            f.surroundInterface.push_back(static_cast<std::uint32_t>(n));
            f.surroundOuter.push_back(static_cast<std::uint32_t>(n));
            continue;
        }
        if (nodeAt(f.patch.mesh, x, y, z) < f.patch.nodeCount())
            f.surroundInterface.push_back(static_cast<std::uint32_t>(n));
    }
    std::sort(f.surroundInterface.begin(), f.surroundInterface.end());

    for (std::size_t n : f.punchPatch) {
        const std::size_t w = nodeAt(f.whole, f.patch.mesh.position[3 * n],
                                     f.patch.mesh.position[3 * n + 1],
                                     f.patch.mesh.position[3 * n + 2]);
        if (w >= f.whole.nodeCount())
            f.problems.push_back("a punch node has no counterpart in the whole plate");
        else
            f.punchWhole.push_back(w);
    }
    return f;
}

// Clamp the rim and prescribe the punch, on any of the three meshes.
void clampRim(HexMesh& mesh) {
    mesh.fixed.assign(mesh.nodeCount() * 3, 0u);
    mesh.prescribed.assign(mesh.nodeCount() * 3, 0.0);
    for (std::size_t n = 0; n < mesh.nodeCount(); ++n)
        if (onRim(mesh.position[3 * n], mesh.position[3 * n + 1]))
            for (int k = 0; k < 3; ++k) mesh.fixed[3 * n + static_cast<std::size_t>(k)] = 1u;
}

void prescribePunch(HexMesh& mesh, const std::vector<std::size_t>& nodes) {
    for (std::size_t n : nodes)
        for (int k = 0; k < 3; ++k) {
            mesh.fixed[3 * n + static_cast<std::size_t>(k)] = 1u;
            mesh.prescribed[3 * n + static_cast<std::size_t>(k)] = k == 2 ? -kPunchTravel : 0.0;
        }
}

double peakOf(const std::vector<double>& u) {
    double peak = 0;
    for (double v : u) peak = std::max(peak, std::fabs(v));
    return peak;
}

// --- 1. The whole plate, three ways ------------------------------------------------

struct Coupled {
    double force = 0;
    std::vector<double> zoneField;      // the zone's own nodal displacement
    std::vector<double> surroundField;  // the surroundings', recovered
    double driveLargest = 0;
    std::size_t driveCount = 0;
    int assemblySize = 0;  // so a mode sweep that kept no modes cannot pass as one
    std::vector<std::string> problems;
};

// One coupled solve: reduce both components at `modes`, assemble, prescribe the
// punch on the assembly, drive the zone's perimeter with what comes out, and solve
// the zone.
Coupled solveCoupled(const Fixture& f, const HexMesh& zoneMesh,
                     const std::vector<std::uint32_t>& zoneInterface,
                     const std::vector<std::size_t>& punchZoneNodes, int modes) {
    Coupled out;
    reduction::Substructure surroundings(f.surround, steel(), f.surroundInterface);
    reduction::Substructure zoneSub(zoneMesh, steel(), zoneInterface);
    if (!surroundings.ready() || !zoneSub.ready()) {
        out.problems.push_back("a substructure would not factor");
        return out;
    }
    reduction::ReduceParams params;
    params.modes = modes;
    params.cutoffFrequency = 0;
    const reduction::Reduction rs = reduction::craigBampton(surroundings, params);
    const reduction::Reduction rz = reduction::craigBampton(zoneSub, params);

    const coupling::Coupling link = coupling::couple(surroundings, rs, zoneSub, rz);
    out.assemblySize = link.assembly.size();
    if (!link.ready()) {
        out.problems.push_back("the coupling is not ready");
        for (const std::string& p : link.problems) out.problems.push_back(p);
        return out;
    }

    // Hold the rim at zero and prescribe the punch, both through boundary DOF of
    // their own component -- which is the only place a reduced model has anything
    // to hold.
    std::vector<std::uint32_t> held;
    for (std::size_t b = 0; b < surroundings.boundaryDof().size(); ++b) {
        const std::uint32_t global = surroundings.boundaryDof()[b];
        const std::size_t node = global / 3;
        if (!onRim(f.surround.position[3 * node], f.surround.position[3 * node + 1])) continue;
        if (link.surroundDof[b] >= 0)
            held.push_back(static_cast<std::uint32_t>(link.surroundDof[b]));
    }
    std::vector<coupling::Prescribed> punch;
    for (std::size_t b = 0; b < zoneSub.boundaryDof().size(); ++b) {
        const std::uint32_t global = zoneSub.boundaryDof()[b];
        const std::size_t node = global / 3, axis = global % 3;
        if (std::find(punchZoneNodes.begin(), punchZoneNodes.end(), node) == punchZoneNodes.end())
            continue;
        if (link.zoneDof[b] < 0) continue;
        punch.push_back({static_cast<std::uint32_t>(link.zoneDof[b]),
                         axis == 2 ? -kPunchTravel : 0.0});
    }
    if (held.empty() || punch.empty()) {
        out.problems.push_back("the assembly has nothing held or nothing prescribed");
        return out;
    }

    std::vector<double> state;
    std::string problem;
    const std::vector<double> load(static_cast<std::size_t>(link.assembly.size()), 0.0);
    if (!coupling::prescribedStaticSolve(link.assembly, load, held, punch, state, &problem)) {
        out.problems.push_back("the assembled solve failed: " + problem);
        return out;
    }

    const coupling::EdgeDrive drive = coupling::edgeDrive(link, zoneSub, state);
    out.driveLargest = drive.largest;
    out.driveCount = drive.count;
    for (const std::string& p : drive.problems) out.problems.push_back(p);

    // The zone, driven rather than clamped. The punch is prescribed here as well:
    // it is a load on the zone, not on the surroundings, and the assembly carried
    // it only so the surroundings could feel it.
    HexMesh driven = zoneMesh;
    driven.fixed.assign(driven.nodeCount() * 3, 0u);
    driven.prescribed.assign(driven.nodeCount() * 3, 0.0);
    if (coupling::applyEdgeDrive(drive, driven) != drive.count)
        out.problems.push_back("the drive named a DOF the zone mesh does not have");
    prescribePunch(driven, punchZoneNodes);

    std::vector<double> u;
    if (!solidshell::solveStatic(driven, steel(), solidshell::Formulation::SolidShell, {}, u,
                                 &problem)) {
        out.problems.push_back("the driven zone would not solve: " + problem);
        return out;
    }
    out.zoneField = u;
    out.force = punchReaction(driven, steel(), u, punchZoneNodes).force;

    // What the surroundings believe, recovered through the same reduction. This is
    // the half of the coupling that only exists if the zone's stiffness is really
    // in the assembled model.
    out.surroundField = reduction::recover(
        surroundings, rs, reduction::componentState(link.assembly, link.assembly.fromA, state));
    return out;
}

void testCoupledZoneIsTheWholePlate() {
    std::printf("\n--- coupling: a driven zone against the monolithic plate ---\n");
    const Fixture f = buildFixture();
    expectTrue("the fixture built", f.problems.empty());
    expectTrue("the patch meshed", !f.patch.empty());
    // A stiffener line pinning interior nodes would make the zone a different
    // structure from the piece of plate it stands for, and every comparison below
    // meaningless. The plate carries no members, so this is zero rather than small.
    expectEqual("no interior node is pinned by a stiffener",
                static_cast<long long>(f.patch.stiffenerNodes), 0);
    expectEqualCount("the patch is the middle quarter of the plate", f.patch.elementCount(),
                     static_cast<std::size_t>(kElements / 2) *
                         static_cast<std::size_t>(kElements / 2));
    expectEqualCount("the surroundings are the rest of it", f.surround.elementCount(),
                     f.whole.elementCount() - f.patch.elementCount());
    expectEqualCount("four punch points, two nodes each", f.punchPatch.size(), 8u);
    expectEqualCount("and the whole plate has the same ones", f.punchWhole.size(), 8u);

    // Splitting by element must not create or destroy steel. A node split would
    // double the interface mass and every answer below would be quietly soft.
    {
        reduction::Substructure surroundings(f.surround, steel(), f.surroundInterface);
        reduction::Substructure zoneSub(f.patch.mesh, steel(), f.zoneInterface);
        expectTrue("both substructures factor", surroundings.ready() && zoneSub.ready());
        const double exact = steel().density * kSide * kSide * kThickness;
        expectNear("the two pieces weigh what the plate weighs",
                   surroundings.totalMass() + zoneSub.totalMass(), exact, 1e-9 * exact);
    }

    // --- The reference ------------------------------------------------------------
    HexMesh mono = f.whole;
    clampRim(mono);
    prescribePunch(mono, f.punchWhole);
    std::vector<double> uMono;
    std::string problem;
    expectTrue("the monolithic plate solves",
               solidshell::solveStatic(mono, steel(), solidshell::Formulation::SolidShell, {},
                                       uMono, &problem));
    const Reaction reference = punchReaction(mono, steel(), uMono, f.punchWhole);
    const double peak = peakOf(uMono);
    expectTrue("the punch actually moves the plate", peak > 0.5 * kPunchTravel);
    // The reaction is read out of the field rather than out of the solver, so the
    // residual at every free DOF says whether that field is in equilibrium at all.
    expectTrue("the reference field is in equilibrium",
               reference.residual < 1e-9 * reference.force);
    // Two independent routes to the same force: the sum of the reactions, and the
    // strain energy, which for a single prescribed travel is `F * travel / 2`.
    expectNear("and its strain energy agrees with its reaction", reference.energy,
               0.5 * reference.force * kPunchTravel, 1e-9 * reference.energy);
    std::printf("     monolithic: %.6f MN at %.1f um, peak %.3e m\n", reference.force * 1e-6,
                kPunchTravel * 1e6, peak);

    // --- The two bounds -----------------------------------------------------------
    HexMesh clamped = f.patch.mesh;
    clamped.fixed.assign(clamped.nodeCount() * 3, 0u);
    clamped.prescribed.assign(clamped.nodeCount() * 3, 0.0);
    for (std::uint32_t n : f.zonePerimeter)
        for (int k = 0; k < 3; ++k) clamped.fixed[3 * n + static_cast<std::size_t>(k)] = 1u;
    prescribePunch(clamped, f.punchPatch);
    std::vector<double> uClamped;
    expectTrue("the clamped zone solves",
               solidshell::solveStatic(clamped, steel(), solidshell::Formulation::SolidShell, {},
                                       uClamped, &problem));
    const double clampedForce = punchReaction(clamped, steel(), uClamped, f.punchPatch).force;

    HexMesh loose = f.patch.mesh;
    loose.fixed.assign(loose.nodeCount() * 3, 0u);
    loose.prescribed.assign(loose.nodeCount() * 3, 0.0);
    prescribePunch(loose, f.punchPatch);
    std::vector<double> uFree;
    // A free-edged patch is not a mechanism here: the punch grips eight nodes in
    // all three directions, which removes every rigid body mode by itself. That is
    // what makes "free" a *bound* rather than a singular system, and it is why the
    // load is a prescribed travel rather than a force.
    expectTrue("the free-edged zone solves",
               solidshell::solveStatic(loose, steel(), solidshell::Formulation::SolidShell, {},
                                       uFree, &problem));
    const double freeForce = punchReaction(loose, steel(), uFree, f.punchPatch).force;

    // **The free end of the force bracket is exactly zero, and that is a fact about
    // the problem rather than a slack tolerance.** A patch with nothing holding its
    // edge, given the same displacement at every gripped punch node, answers with a
    // rigid translation: no strain anywhere, no reaction at all. So "the coupled
    // force lies between the free one and the clamped one" is a *vacuous* lower
    // bound however the numbers come out, and asserting it would prove nothing.
    // It is asserted here as the identity it is, and the bracket that carries the
    // content is the edge-following ratio below.
    std::printf("     clamped edge %.6f MN, free edge %.3e N\n", clampedForce * 1e-6, freeForce);
    expectTrue("a clamped edge is stiffer than the real one", clampedForce > reference.force);
    expectTrue("a free edge resists a rigid punch with exactly nothing",
               std::fabs(freeForce) < 1e-9 * reference.force);
    // The vacuity guard on the stiff end. If the punch's answer did not reach the
    // zone's edge the two would agree and an uncoupled zone would pass every test
    // below. Measured at 433x, asserted at 10.
    expectTrue("and clamping the edge changes the answer by a factor, not a rounding",
               clampedForce > 10.0 * reference.force);

    // --- The bracket that is not vacuous: how far the edge follows ------------------
    //
    // The quantity at issue is not a force but the perimeter's own displacement,
    // and there both bounds are exact closed forms rather than measurements: a
    // clamped edge follows nothing (0) and a free one follows the punch exactly (1,
    // the rigid translation above). The truth is strictly inside, and a coupled
    // answer that landed on either end would be one of the two models this file
    // exists to replace.
    double edgeFollow = 0;
    for (std::uint32_t n : f.zonePerimeter) {
        const std::size_t w =
            nodeAt(f.whole, f.patch.mesh.position[3 * n], f.patch.mesh.position[3 * n + 1],
                   f.patch.mesh.position[3 * n + 2]);
        if (w >= f.whole.nodeCount()) continue;
        edgeFollow = std::max(edgeFollow, std::fabs(uMono[3 * w + 2]) / kPunchTravel);
    }
    std::printf("     the edge follows %.3f of the punch: clamped says 0, free says 1\n",
                edgeFollow);
    expectTrue("the real edge is nowhere near clamped", edgeFollow > 0.1);
    expectTrue("and nowhere near free", edgeFollow < 0.9);

    // --- The coupled answer, at three mode counts ----------------------------------
    //
    // Mode count is swept not because the answer should converge but because it
    // should not move at all: `reduction.hpp` §1 property 1 says the interface is
    // exact at zero modes, and a coupling that quietly lost that property would
    // still improve with modes and still look reasonable.
    double firstForce = 0;
    int firstSize = 0;
    std::vector<double> firstZone;
    for (int modes : {0, 4, 12}) {
        const Coupled coupled = solveCoupled(f, f.patch.mesh, f.zoneInterface, f.punchPatch, modes);
        expectTrue("the coupled solve raises nothing", coupled.problems.empty());
        if (!coupled.problems.empty()) {
            for (const std::string& p : coupled.problems) std::printf("     ! %s\n", p.c_str());
            return;
        }

        // The drive has to be a real displacement field, not a well-formed zero.
        expectEqualCount("every perimeter DOF is driven", coupled.driveCount,
                         3 * f.zonePerimeter.size());
        expectTrue("and the drive is a real fraction of the plate's own movement",
                   coupled.driveLargest > 0.05 * peak);

        double worst = 0;
        for (std::size_t n = 0; n < f.patch.nodeCount(); ++n) {
            const std::size_t w =
                nodeAt(f.whole, f.patch.mesh.position[3 * n], f.patch.mesh.position[3 * n + 1],
                       f.patch.mesh.position[3 * n + 2]);
            if (w >= f.whole.nodeCount()) continue;
            for (int k = 0; k < 3; ++k)
                worst = std::max(worst,
                                 std::fabs(coupled.zoneField[3 * n + static_cast<std::size_t>(k)] -
                                           uMono[3 * w + static_cast<std::size_t>(k)]));
        }
        std::printf("     %2d modes: coupled %.6f MN (monolithic %.6f), field to %.2e m of %.3e\n",
                    modes, coupled.force * 1e-6, reference.force * 1e-6, worst, peak);

        // Inside the bracket at the end where the bracket means something, and on
        // the ratio where both ends are exact.
        expectTrue("the coupled zone is softer than a clamped one",
                   coupled.force > freeForce && coupled.force < clampedForce);
        double coupledFollow = 0;
        for (std::uint32_t n : f.zonePerimeter)
            coupledFollow = std::max(
                coupledFollow, std::fabs(coupled.zoneField[3 * n + 2]) / kPunchTravel);
        expectNear("and its edge follows the plate by exactly as much as the plate's does",
                   coupledFollow, edgeFollow, 1e-9 * edgeFollow);
        // Exact, not close. 1e-9 relative is three orders past what "well
        // converged" would deliver and is what was measured; a 1e-3 assertion here
        // would pass on a coupling that had lost the property entirely.
        expectNear("and reproduces the monolithic reaction", coupled.force, reference.force,
                   1e-9 * reference.force);
        expectTrue("and the monolithic displacement field inside the zone", worst < 1e-9 * peak);

        if (modes == 0) {
            firstForce = coupled.force;
            firstSize = coupled.assemblySize;
            firstZone = coupled.zoneField;
        } else {
            // Non-vacuity: the modes have to actually be *in* the assembled model,
            // or "adding modes changes nothing" is a statement about a sweep that
            // never added any.
            expectEqual("the modes really were kept",
                        static_cast<long long>(coupled.assemblySize),
                        static_cast<long long>(firstSize + 2 * modes));
            expectTrue("adding modes does not move a static coupled answer",
                       std::fabs(coupled.force - firstForce) < 1e-9 * firstForce);
            double drift = 0;
            for (std::size_t d = 0; d < firstZone.size(); ++d)
                drift = std::max(drift, std::fabs(coupled.zoneField[d] - firstZone[d]));
            expectTrue("nor its displacement field", drift < 1e-9 * peak);
        }
    }

    // --- The negative control -------------------------------------------------------
    //
    // The same instrument, pointed at the zone the engine had before this file
    // existed. If it did not miss by orders of magnitude, the agreement above would
    // be measuring a well-converged mesh rather than a coupling.
    double clampedWorst = 0;
    for (std::size_t n = 0; n < f.patch.nodeCount(); ++n) {
        const std::size_t w =
            nodeAt(f.whole, f.patch.mesh.position[3 * n], f.patch.mesh.position[3 * n + 1],
                   f.patch.mesh.position[3 * n + 2]);
        if (w >= f.whole.nodeCount()) continue;
        for (int k = 0; k < 3; ++k)
            clampedWorst =
                std::max(clampedWorst, std::fabs(uClamped[3 * n + static_cast<std::size_t>(k)] -
                                                 uMono[3 * w + static_cast<std::size_t>(k)]));
    }
    std::printf("     control: the clamped zone misses the monolithic field by %.2e m (%.0f%%)\n",
                clampedWorst, 100.0 * clampedWorst / peak);
    expectTrue("a clamped zone is not the monolithic answer", clampedWorst > 0.1 * peak);
    expectTrue("nor is its reaction", clampedForce > 1.05 * reference.force);
}

// --- 2. The zone reports back ------------------------------------------------------

void testDamageIsFeltOutsideTheZone() {
    std::printf("\n--- coupling: a torn zone hands its softness back ---\n");
    const Fixture f = buildFixture();
    expectTrue("the fixture built", f.problems.empty());

    // A slit across the zone: the row of elements whose centroids lie in the strip
    // 0 < y < h. It sits directly under the punch, so the softening it causes
    // cannot be mistaken for a rounding, and it stops short of the perimeter, so no
    // interface node is lost -- the case where a tear reaches the edge is a
    // different one and is exercised below.
    const double h = kSide / kElements;
    std::vector<std::uint8_t> tornPatch(f.patch.elementCount(), 0u);
    std::vector<std::uint8_t> tornWhole(f.whole.elementCount(), 0u);
    std::size_t cutCount = 0;
    for (std::size_t e = 0; e < f.patch.elementCount(); ++e) {
        double centre[3];
        elementCentroid(f.patch.mesh, e, centre);
        if (!(centre[1] > 0.0 && centre[1] < h)) continue;
        tornPatch[e] = 1u;
        ++cutCount;
        for (std::size_t w = 0; w < f.whole.elementCount(); ++w) {
            double other[3];
            elementCentroid(f.whole, w, other);
            if (std::fabs(centre[0] - other[0]) < 1e-9 && std::fabs(centre[1] - other[1]) < 1e-9)
                tornWhole[w] = 1u;
        }
    }
    expectEqualCount("the slit is a full row of the zone", cutCount,
                     static_cast<std::size_t>(kElements / 2));

    const coupling::DamagedMesh damaged = coupling::withoutElements(f.patch.mesh, tornPatch);
    expectEqualCount("the damaged zone lost exactly the torn elements", damaged.removedElements,
                     cutCount);
    expectEqualCount("and no node with it", damaged.orphanedNodes, 0u);
    expectEqualCount("so it keeps every node it had", damaged.mesh.nodeCount(),
                     f.patch.nodeCount());

    const std::vector<std::uint32_t> damagedInterface =
        coupling::carriedInterface(damaged, f.zoneInterface);
    expectEqualCount("and every interface node", damagedInterface.size(), f.zoneInterface.size());
    std::vector<std::size_t> punchDamaged;
    for (std::size_t n : f.punchPatch) {
        const auto found = std::lower_bound(damaged.node.begin(), damaged.node.end(),
                                            static_cast<std::uint32_t>(n));
        expectTrue("the punch survives the slit",
                   found != damaged.node.end() && *found == static_cast<std::uint32_t>(n));
        if (found != damaged.node.end() && *found == static_cast<std::uint32_t>(n))
            punchDamaged.push_back(static_cast<std::size_t>(found - damaged.node.begin()));
    }

    // The references: the whole plate, intact and slit, solved by the element code.
    const auto monolithic = [&](const std::vector<std::uint8_t>& removed) {
        const coupling::DamagedMesh cut = coupling::withoutElements(f.whole, removed);
        HexMesh mesh = cut.mesh;
        clampRim(mesh);
        std::vector<std::size_t> punch;
        for (std::size_t n : f.punchWhole) {
            const auto found = std::lower_bound(cut.node.begin(), cut.node.end(),
                                                static_cast<std::uint32_t>(n));
            if (found != cut.node.end() && *found == static_cast<std::uint32_t>(n))
                punch.push_back(static_cast<std::size_t>(found - cut.node.begin()));
        }
        prescribePunch(mesh, punch);
        std::vector<double> u;
        std::string problem;
        const bool ok = solidshell::solveStatic(mesh, steel(),
                                                solidshell::Formulation::SolidShell, {}, u,
                                                &problem);
        expectTrue("the monolithic reference solves", ok);
        struct Answer {
            HexMesh mesh;
            std::vector<std::uint32_t> node;
            std::vector<double> u;
            double force;
        };
        return Answer{mesh, cut.node, u, punchReaction(mesh, steel(), u, punch).force};
    };
    const auto intactRef = monolithic(std::vector<std::uint8_t>(f.whole.elementCount(), 0u));
    const auto slitRef = monolithic(tornWhole);
    std::printf("     monolithic: intact %.6f MN, slit %.6f MN (%.1f%% softer)\n",
                intactRef.force * 1e-6, slitRef.force * 1e-6,
                100.0 * (1.0 - slitRef.force / intactRef.force));
    // Vacuity guard: a slit that barely changed the plate would make everything
    // below agree for the wrong reason.
    expectTrue("the slit really softens the plate", slitRef.force < 0.9 * intactRef.force);

    const Coupled intact = solveCoupled(f, f.patch.mesh, f.zoneInterface, f.punchPatch, 0);
    const Coupled slit = solveCoupled(f, damaged.mesh, damagedInterface, punchDamaged, 0);
    expectTrue("both coupled solves raise nothing",
               intact.problems.empty() && slit.problems.empty());
    if (!intact.problems.empty() || !slit.problems.empty()) {
        for (const std::string& p : intact.problems) std::printf("     ! %s\n", p.c_str());
        for (const std::string& p : slit.problems) std::printf("     ! %s\n", p.c_str());
        return;
    }
    expectNear("the coupled slit zone reproduces the monolithic slit plate", slit.force,
               slitRef.force, 1e-9 * slitRef.force);

    // Direction 2, and it is a statement about the *surroundings*, not the zone.
    // The zone's own answer would change if it were merely driven; what says the
    // damage was reported back is that the plating outside it moved.
    double moved = 0, worstOutside = 0;
    const double peakOutside = peakOf(slitRef.u);
    for (std::size_t n = 0; n < f.surround.nodeCount(); ++n) {
        const std::size_t r = nodeAt(slitRef.mesh, f.surround.position[3 * n],
                                     f.surround.position[3 * n + 1],
                                     f.surround.position[3 * n + 2]);
        if (r >= slitRef.mesh.nodeCount()) continue;
        for (int k = 0; k < 3; ++k) {
            const std::size_t d = static_cast<std::size_t>(k);
            moved = std::max(moved, std::fabs(slit.surroundField[3 * n + d] -
                                              intact.surroundField[3 * n + d]));
            worstOutside = std::max(worstOutside,
                                    std::fabs(slit.surroundField[3 * n + d] - slitRef.u[3 * r + d]));
        }
    }
    std::printf("     surroundings: moved %.2e m when the zone tore, and match the slit plate to "
                "%.2e m of %.3e\n",
                moved, worstOutside, peakOutside);
    // The measurement that says direction 2 exists at all: without the zone's own
    // stiffness in the assembled model, tearing it could not move anything outside.
    expectTrue("the surroundings feel the tear", moved > 0.02 * peakOutside);
    // And they feel the *right* tear: the recovered field is the slit plate's own.
    expectTrue("and what they feel is the slit plate's own field",
               worstOutside < 1e-9 * peakOutside);
    // The negative control for direction 2. The intact surroundings are what a
    // one-way coupling would still be reporting, and they are wrong by far more
    // than the tolerance above -- so that tolerance is measuring something.
    double stale = 0;
    for (std::size_t n = 0; n < f.surround.nodeCount(); ++n) {
        const std::size_t r = nodeAt(slitRef.mesh, f.surround.position[3 * n],
                                     f.surround.position[3 * n + 1],
                                     f.surround.position[3 * n + 2]);
        if (r >= slitRef.mesh.nodeCount()) continue;
        for (int k = 0; k < 3; ++k) {
            const std::size_t d = static_cast<std::size_t>(k);
            stale = std::max(stale,
                             std::fabs(intact.surroundField[3 * n + d] - slitRef.u[3 * r + d]));
        }
    }
    expectTrue("where a zone that never reported back would be wrong", stale > 1e3 * worstOutside);
}

// --- 3. Why the Tier-0 shape does not transfer --------------------------------------
//
// `promotion.hpp` §5 reports damage to Tier 0 as an effective thickness, and the
// question this file has to answer is whether that is the right shape for Tier 1.
// It is not, and the reason is geometric rather than a matter of taste: a
// solid-shell carries its thickness in the positions of its nodes, so thinning a
// zone moves the very nodes it shares with the structure round it.
void testAThinnerZoneNoLongerFits() {
    std::printf("\n--- coupling: a thinner zone is a different component ---\n");
    const Fixture f = buildFixture();
    expectTrue("the fixture built", f.problems.empty());

    const double knockdown = 0.8;
    StructuralMesh thinStructure = flatPlate();
    for (PlatePanel& p : thinStructure.panels) p.thickness = kThickness * knockdown;
    const zone::Patch thin = zone::buildPatch(thinStructure, {0, 0, 0}, patchParams());
    expectEqualCount("the thinned zone meshes the same elements", thin.elementCount(),
                     f.patch.elementCount());

    reduction::Substructure surroundings(f.surround, steel(), f.surroundInterface);
    reduction::Substructure thinZone(thin.mesh, steel(), reduction::nodesPinned(thin.mesh));
    expectTrue("both factor", surroundings.ready() && thinZone.ready());

    const reduction::InterfaceMap match = reduction::matchBoundaries(surroundings, thinZone, 1e-9);
    expectEqualCount("a thinned zone shares no interface DOF at all with its surroundings",
                     match.shared, 0u);
    expectTrue("and says so rather than assembling two models side by side",
               !match.problems.empty());

    // The offset is exactly half the thickness lost, on both faces, which is the
    // closed form and not an observation: the mid-surface stays put and the faces
    // move to +/- t'/2.
    const double offset = 0.5 * kThickness * (1.0 - knockdown);
    const reduction::InterfaceMap loose =
        reduction::matchBoundaries(surroundings, thinZone, 2.0 * offset);
    expectNear("the interface nodes are exactly (t - t')/2 out of place", loose.worstGap, offset,
               1e-12);
    std::printf("     a %.0f%% thickness knockdown moves every interface node %.3f mm\n",
                100.0 * (1.0 - knockdown), offset * 1e3);
    // Which is the whole point: at a tolerance loose enough to match, the interface
    // is no longer coincident, so the assembled model would be tying together nodes
    // that are a millimetre apart and calling it one displacement.
    expectTrue("and a tolerance loose enough to match them is a millimetre of geometry error",
               loose.shared > 0 && offset > 1e-3 * kThickness);

    // Element deletion is the shape that *does* transfer, and here is why in one
    // number: it leaves every surviving interface node exactly where it was.
    std::vector<std::uint8_t> torn(f.patch.elementCount(), 0u);
    torn[0] = 1u;
    const coupling::DamagedMesh damaged = coupling::withoutElements(f.patch.mesh, torn);
    reduction::Substructure damagedZone(
        damaged.mesh, steel(), coupling::carriedInterface(damaged, f.zonePerimeter));
    expectTrue("the damaged zone factors", damagedZone.ready());
    const reduction::InterfaceMap kept =
        reduction::matchBoundaries(surroundings, damagedZone, 1e-9);
    expectTrue("a zone that lost an element still matches at the coincidence tolerance",
               kept.shared > 0);
    expectTrue("with no gap at all", kept.worstGap < 1e-12);
    std::printf("     deleting an element instead: %zu DOF still shared, worst gap %.1e m\n",
                kept.shared, kept.worstGap);
}

// --- 4. The explicit Tier-2 solver follows the drive ----------------------------------

void testTheExplicitSolverFollowsItsBoundary() {
    std::printf("\n--- coupling: the Tier-2 explicit solver on a driven edge ---\n");
    const StructuralMesh structure = flatPlate();
    zone::Patch patch = zone::buildPatch(structure, {0, 0, 0}, patchParams());
    expectTrue("the patch meshed", !patch.empty());
    const std::vector<std::uint32_t> perimeter = reduction::nodesPinned(patch.mesh);
    expectTrue("it has a perimeter to drive", !perimeter.empty());

    // A rigid translation of the whole perimeter. The static answer is a rigid
    // translation of the whole patch, exactly, carrying no strain energy -- which
    // is a statement the co-rotational element makes and a deflection test cannot.
    const Vec3 travel{7.0e-4, -3.0e-4, 5.0e-4};
    for (std::uint32_t n : perimeter)
        for (int k = 0; k < 3; ++k) patch.mesh.prescribed[3 * n + static_cast<std::size_t>(k)] =
            travel[k];

    zone::SolveParams solve;
    solve.plastic = false;
    solve.indenter.halfLength = 0.0;  // no punch: the boundary is the only agent
    solve.edgeRamp = 4000.0 * patch.criticalTimestep;
    solve.duration = 30000.0 * patch.criticalTimestep;
    solve.damping = 0.9985;
    zone::Solver solver(patch, plasticity::shipSteel(), solve);
    expectEqual("every perimeter DOF is driven", static_cast<long long>(solver.result().drivenEdgeDof),
                static_cast<long long>(3 * perimeter.size()));
    solver.run();

    // The identity first: a driven DOF ends exactly where it was told to be. Not
    // "close" -- the target is computed from the meshed position and the step is
    // chosen to land on it, so anything but machine precision is a defect in the
    // drive rather than in the physics.
    double worstBoundary = 0;
    for (std::uint32_t n : perimeter)
        for (int k = 0; k < 3; ++k) {
            const std::size_t d = 3 * n + static_cast<std::size_t>(k);
            worstBoundary = std::max(
                worstBoundary,
                std::fabs(solver.position()[d] - (patch.mesh.position[d] + travel[k])));
        }
    const double magnitude = std::sqrt(dot(travel, travel));
    expectTrue("the drive asked for a real displacement", magnitude > 1e-5);
    expectTrue("and every driven DOF is exactly where the drive put it",
               worstBoundary < 1e-15 * magnitude);

    // And the physics: the interior followed, and stored nothing doing it.
    double worstInterior = 0;
    for (std::size_t n = 0; n < patch.nodeCount(); ++n)
        for (int k = 0; k < 3; ++k) {
            const std::size_t d = 3 * n + static_cast<std::size_t>(k);
            worstInterior =
                std::max(worstInterior,
                         std::fabs(solver.position()[d] - (patch.mesh.position[d] + travel[k])));
        }
    // The scale a "no strain energy" claim is a zero against: what the patch would
    // store if it were strained by the same distance across one element.
    const double reference = 0.5 * steel().youngsModulus * patch.area * kThickness *
                             (magnitude / (kSide / kElements)) * (magnitude / (kSide / kElements));
    std::printf("     rigid drive: interior to %.2e m of %.2e, strain energy %.2e J of %.2e\n",
                worstInterior, magnitude, solver.result().strainEnergy, reference);
    // Measured at 3.4e-14 m against a 9.1e-4 m travel -- machine precision, because
    // a rigid translation is an *exact* equilibrium of the discrete system and the
    // damped run settles onto it rather than approaching it. Asserted at 1e-9
    // relative, which is still five orders inside what a 1e-3 tolerance would have
    // allowed a patch that was merely being dragged along late.
    expectTrue("the whole patch followed its boundary", worstInterior < 1e-9 * magnitude);
    expectTrue("and stored no strain energy doing it",
               std::fabs(solver.result().strainEnergy) < 1e-10 * reference);
    expectTrue("the reference energy is a real number", reference > 1.0);

    // The energy account still closes with a second agent on the patch. Nothing
    // else in the result would notice a boundary that was quietly doing work.
    const zone::SolveResult& result = solver.result();
    const double throughput = std::fabs(result.boundaryWork) + result.strainEnergy +
                              result.kinetic + result.dampingLoss;
    std::printf("     boundary work %.3e J, residual %+.3f%% of %.3e J\n", result.boundaryWork,
                100.0 * result.energyResidual() / throughput, throughput);
    // One per cent, the same tolerance `test_zone.cpp` holds the indenter's account
    // to and for the same reason: what is left is the integrator's own order, not a
    // missing term.
    expectTrue("and the energy account closes",
               std::fabs(result.energyResidual()) < 0.01 * throughput);
    // Which it would not, by two orders of magnitude, if the boundary's work were
    // left out of the account -- so the closure above is measuring something.
    expectTrue("with the boundary genuinely doing the work being accounted for",
               std::fabs(result.boundaryWork) > 50.0 * std::fabs(result.energyResidual()));
}

// The whole chain: a Tier-1 reduced model of the plate round the patch drives the
// patch's perimeter, the Tier-2 **explicit** solver runs on it, and the answer is
// compared to the monolithic plate. Every link is exercised on the path a caller
// would actually take.
void testTheWholeChainAgainstTheWholePlate() {
    std::printf("\n--- coupling: Tier 1 driving the Tier-2 explicit solver ---\n");
    const Fixture f = buildFixture();
    expectTrue("the fixture built", f.problems.empty());

    HexMesh mono = f.whole;
    clampRim(mono);
    prescribePunch(mono, f.punchWhole);
    std::vector<double> uMono;
    std::string problem;
    expectTrue("the monolithic plate solves",
               solidshell::solveStatic(mono, steel(), solidshell::Formulation::SolidShell, {},
                                       uMono, &problem));
    const double peak = peakOf(uMono);

    const Coupled coupled = solveCoupled(f, f.patch.mesh, f.zoneInterface, f.punchPatch, 0);
    expectTrue("the coupled static solve raises nothing", coupled.problems.empty());
    if (!coupled.problems.empty()) return;

    // The same boundary condition, handed to the explicit solver instead of the
    // direct one. The punch is prescribed rather than driven by `Indenter` because
    // the comparison is against a *static* plate: a kinematic punch travelling at a
    // speed would be answering a different question.
    zone::Patch patch = f.patch;
    patch.mesh.fixed.assign(patch.nodeCount() * 3, 0u);
    patch.mesh.prescribed.assign(patch.nodeCount() * 3, 0.0);
    for (std::size_t n = 0; n < patch.nodeCount(); ++n)
        for (int k = 0; k < 3; ++k) {
            const std::size_t d = 3 * n + static_cast<std::size_t>(k);
            const double target = coupled.zoneField[d];
            const bool perimeter = f.patch.mesh.fixed[d] != 0u;
            const bool punch = std::find(f.punchPatch.begin(), f.punchPatch.end(), n) !=
                               f.punchPatch.end();
            if (!perimeter && !punch) continue;
            patch.mesh.fixed[d] = 1u;
            patch.mesh.prescribed[d] = perimeter ? target : (k == 2 ? -kPunchTravel : 0.0);
        }

    zone::SolveParams solve;
    solve.plastic = false;
    solve.indenter.halfLength = 0.0;
    solve.edgeRamp = 4000.0 * patch.criticalTimestep;
    solve.duration = 30000.0 * patch.criticalTimestep;
    solve.damping = 0.9985;
    zone::Solver solver(patch, plasticity::shipSteel(), solve);
    expectTrue("the drive reaches both the perimeter and the punch",
               solver.result().drivenEdgeDof > static_cast<int>(3 * f.punchPatch.size()));
    solver.run();

    double worst = 0;
    for (std::size_t n = 0; n < patch.nodeCount(); ++n) {
        const std::size_t w =
            nodeAt(f.whole, patch.mesh.position[3 * n], patch.mesh.position[3 * n + 1],
                   patch.mesh.position[3 * n + 2]);
        if (w >= f.whole.nodeCount()) continue;
        for (int k = 0; k < 3; ++k) {
            const std::size_t d = static_cast<std::size_t>(k);
            worst = std::max(worst, std::fabs((solver.position()[3 * n + d] -
                                               patch.mesh.position[3 * n + d]) -
                                              uMono[3 * w + d]));
        }
    }
    std::printf("     explicit zone settles to the monolithic field within %.2e m of %.3e "
                "(%.2f%%), %d steps\n",
                worst, peak, 100.0 * worst / peak, solver.result().steps);
    expectTrue("the explicit solver ran", solver.result().steps > 10000);
    // An explicit solve settling under damping is not a direct solve and does not
    // pretend to be: what is asserted is that it reaches the monolithic answer to
    // the accuracy the settling delivers, measured at 5.5e-5 of the peak and
    // asserted at 2e-4.
    expectTrue("and settles onto the monolithic displacement field", worst < 2e-4 * peak);
    // Against the bound it would have had with the edge clamped, so this is not a
    // tolerance wide enough to swallow the difference it is testing.
    double clampedGap = 0;
    for (std::size_t n = 0; n < patch.nodeCount(); ++n) {
        const std::size_t w =
            nodeAt(f.whole, patch.mesh.position[3 * n], patch.mesh.position[3 * n + 1],
                   patch.mesh.position[3 * n + 2]);
        if (w >= f.whole.nodeCount()) continue;
        if (f.patch.mesh.fixed[3 * n] == 0u) continue;
        for (int k = 0; k < 3; ++k)
            clampedGap = std::max(clampedGap, std::fabs(uMono[3 * w + static_cast<std::size_t>(k)]));
    }
    std::printf("     a clamped edge would have been %.2e m out at the perimeter alone\n",
                clampedGap);
    expectTrue("and the clamped edge it replaces is much further out than that tolerance",
               clampedGap > 1000.0 * worst);
}

// --- 4b. The ramp, and the tear a real solve reports ----------------------------------

void testTheRampAndTheTearReadBack() {
    std::printf("\n--- coupling: the edge ramp, and reading a tear off a solver ---\n");
    const StructuralMesh structure = flatPlate();
    zone::Patch patch = zone::buildPatch(structure, {0, 0, 0}, patchParams());
    const std::vector<std::uint32_t> perimeter = reduction::nodesPinned(patch.mesh);
    const Vec3 travel{0.0, 0.0, 4.0e-4};
    for (std::uint32_t n : perimeter)
        for (int k = 0; k < 3; ++k) patch.mesh.prescribed[3 * n + static_cast<std::size_t>(k)] =
            travel[k];

    // Half way through the ramp the drive must be at smoothstep(0.5), which is
    // exactly one half -- the one point of a smoothstep with a closed form, and the
    // one a linear ramp would agree with, so the *quarter* point is checked too
    // where they differ: smoothstep(0.25) = 0.15625 against a linear 0.25.
    zone::SolveParams solve;
    solve.plastic = false;
    solve.indenter.halfLength = 0.0;
    solve.damping = 1.0;
    const double dt = patch.criticalTimestep;
    solve.edgeRamp = 4000.0 * dt;
    solve.duration = 1e9;  // stepped by hand
    zone::Solver solver(patch, plasticity::shipSteel(), solve);
    const std::size_t probe = 3 * perimeter[0] + 2;
    const double start = patch.mesh.position[probe];
    for (int i = 0; i < 1000; ++i) solver.step();
    expectNear("a quarter through the ramp the edge is at smoothstep(1/4), not at a quarter",
               (solver.position()[probe] - start) / travel.z, 0.15625, 1e-12);
    for (int i = 0; i < 1000; ++i) solver.step();
    expectNear("and half way through it is at exactly half",
               (solver.position()[probe] - start) / travel.z, 0.5, 1e-12);
    for (int i = 0; i < 2000; ++i) solver.step();
    expectNear("and at the end of the ramp it is all the way there",
               (solver.position()[probe] - start) / travel.z, 1.0, 1e-12);
    for (int i = 0; i < 500; ++i) solver.step();
    expectNear("and stays there", (solver.position()[probe] - start) / travel.z, 1.0, 1e-12);

    // **A pre-load moves the rest configuration out from under the mesh**, and the
    // drive is a statement about where the surrounding *ship* is, not about where
    // the stress-free configuration went. So the target is the meshed position plus
    // the drive, and a drive measured from `rest()` would carry the pre-load's own
    // displacement field into the boundary condition. The two differ by `e0 * x`,
    // which at 84 MPa over this patch is 8e-5 m -- a fifth of the travel below, so
    // the assertion is not a rounding.
    {
        zone::Patch preloaded = zone::buildPatch(structure, {0, 0, 0}, patchParams());
        const double travelX = 4.0e-4;
        for (std::uint32_t n : perimeter) preloaded.mesh.prescribed[3 * n] = travelX;
        zone::SolveParams pre;
        pre.plastic = false;
        pre.indenter.halfLength = 0.0;
        pre.preload.stress = 84.0e6;
        pre.duration = 40.0 * preloaded.criticalTimestep;
        zone::Solver solver2(preloaded, plasticity::shipSteel(), pre);
        solver2.run();
        double worstDrive = 0, restShift = 0;
        for (std::uint32_t n : perimeter) {
            const std::size_t d = 3 * n;
            worstDrive = std::max(
                worstDrive,
                std::fabs(solver2.position()[d] - (preloaded.mesh.position[d] + travelX)));
            restShift = std::max(restShift,
                                 std::fabs(solver2.rest()[d] - preloaded.mesh.position[d]));
        }
        std::printf("     pre-loaded: the rest configuration moved %.2e m, the drive landed to "
                    "%.2e m\n",
                    restShift, worstDrive);
        expectTrue("a pre-load really did move the rest configuration", restShift > 1e-5);
        expectTrue("and the drive is still measured from the meshed position",
                   worstDrive < 1e-12 * travelX);
    }

    // `withoutTornElements` is on the caller's own path and nothing above reaches
    // it. Driving a real zone to a real tear costs core-minutes, so the tear is
    // handed to the solver through `adopt` -- the same route `gpu::ZoneGpuSolver`
    // uses -- which exercises the reading rather than the tearing. What tears is
    // `test_zone.cpp`'s problem; that the reading matches is this one's.
    zone::Patch plain = zone::buildPatch(structure, {0, 0, 0}, patchParams());
    zone::SolveParams plastic;
    plastic.plastic = true;
    plastic.indenter.halfLength = 0.0;
    plastic.duration = 2.0 * plain.criticalTimestep;
    zone::Solver torn(plain, plasticity::shipSteel(), plastic);
    std::vector<solidshell::ElementPlasticState> state = torn.elementState();
    expectEqualCount("the solver keeps a state per element", state.size(), plain.elementCount());
    state[0].torn = true;
    state[3].torn = true;
    torn.adopt(torn.position(), std::vector<double>(3 * plain.nodeCount(), 0.0), state, 0, 0, 0, 0,
               0);
    const coupling::DamagedMesh read = coupling::withoutTornElements(plain, torn);
    expectEqualCount("what the zone reports torn is what comes out", read.removedElements, 2u);
    expectEqualCount("and the survivors are everything else", read.mesh.elementCount(),
                     plain.elementCount() - 2);
    expectTrue("the elements kept are the ones that did not tear",
               std::find(read.element.begin(), read.element.end(), 0u) == read.element.end() &&
                   std::find(read.element.begin(), read.element.end(), 3u) == read.element.end() &&
                   std::find(read.element.begin(), read.element.end(), 1u) != read.element.end());

    // An elastically solved zone has no plastic state at all, so "nothing torn" and
    // "no damage model" are the same array. They are not the same answer, and the
    // difference is reported rather than left to be assumed.
    zone::SolveParams elastic;
    elastic.plastic = false;
    elastic.indenter.halfLength = 0.0;
    elastic.duration = 2.0 * plain.criticalTimestep;
    zone::Solver elasticSolver(plain, plasticity::shipSteel(), elastic);
    const coupling::DamagedMesh none = coupling::withoutTornElements(plain, elasticSolver);
    expectEqualCount("an elastic zone reports no tear", none.removedElements, 0u);
    expectTrue("and says why rather than looking undamaged", !none.problems.empty());
}

// --- 5. The pieces on their own ------------------------------------------------------

void testPrescribedSolveAndMeshSurgery() {
    std::printf("\n--- coupling: prescribing on a reduced model, and mesh surgery ---\n");

    // `prescribedStaticSolve` against `assembledStaticSolve` on a problem both can
    // answer: prescribing zero has to be exactly holding. Anything else means the
    // right-hand side elimination has a term it should not.
    const Fixture f = buildFixture();
    reduction::Substructure surroundings(f.surround, steel(), f.surroundInterface);
    reduction::Substructure zoneSub(f.patch.mesh, steel(), f.zoneInterface);
    reduction::ReduceParams params;
    params.modes = 0;
    params.cutoffFrequency = 0;
    const reduction::Reduction rs = reduction::craigBampton(surroundings, params);
    const reduction::Reduction rz = reduction::craigBampton(zoneSub, params);
    const coupling::Coupling link = coupling::couple(surroundings, rs, zoneSub, rz);
    expectTrue("the coupling is ready", link.ready());
    expectEqualCount("every perimeter DOF of the zone is shared with the surroundings",
                     link.sharedDof, 3 * f.zonePerimeter.size());
    expectEqualCount("and the punch is the only part of the zone's interface that is not",
                     link.zoneUnshared, 3 * f.punchPatch.size());
    expectTrue("the matched nodes are coincident", link.worstGap < 1e-12);

    std::vector<std::uint32_t> held;
    for (std::size_t b = 0; b < surroundings.boundaryDof().size(); ++b) {
        const std::size_t node = surroundings.boundaryDof()[b] / 3;
        if (onRim(f.surround.position[3 * node], f.surround.position[3 * node + 1]) &&
            link.surroundDof[b] >= 0)
            held.push_back(static_cast<std::uint32_t>(link.surroundDof[b]));
    }
    std::vector<double> load(static_cast<std::size_t>(link.assembly.size()), 0.0);
    // A load on a boundary DOF of the zone, so both components have to carry it.
    load[static_cast<std::size_t>(link.zoneDof[2])] = 1.0e4;

    std::vector<double> byHolding, byPrescribing;
    std::string problem;
    std::vector<std::uint32_t> alsoHeld = held;
    alsoHeld.push_back(static_cast<std::uint32_t>(link.zoneDof[5]));
    expectTrue("the held solve runs", reduction::assembledStaticSolve(link.assembly, load,
                                                                      alsoHeld, byHolding,
                                                                      &problem));
    expectTrue("the prescribed solve runs",
               coupling::prescribedStaticSolve(
                   link.assembly, load, held,
                   {{static_cast<std::uint32_t>(link.zoneDof[5]), 0.0}}, byPrescribing, &problem));
    double difference = 0, magnitude = 0;
    for (std::size_t i = 0; i < byHolding.size(); ++i) {
        difference = std::max(difference, std::fabs(byHolding[i] - byPrescribing[i]));
        magnitude = std::max(magnitude, std::fabs(byHolding[i]));
    }
    expectTrue("the load moves the assembly", magnitude > 1e-12);
    expectTrue("prescribing zero is exactly holding", difference == 0.0);

    // And prescribing a value puts it there, which is the property the whole
    // coupling stands on. Superposition gives the closed form: doubling the
    // prescribed value with no load doubles every displacement.
    const std::uint32_t target = static_cast<std::uint32_t>(link.zoneDof[5]);
    const std::vector<double> noLoad(static_cast<std::size_t>(link.assembly.size()), 0.0);
    std::vector<double> once, twice;
    expectTrue("a prescribed value solves",
               coupling::prescribedStaticSolve(link.assembly, noLoad, held, {{target, 1.0e-3}},
                                               once, &problem));
    expectTrue("and twice the value solves",
               coupling::prescribedStaticSolve(link.assembly, noLoad, held, {{target, 2.0e-3}},
                                               twice, &problem));
    expectNear("the prescribed DOF takes the value it was given",
               once[target], 1.0e-3, 1e-18);
    double worstScale = 0, scale = 0;
    for (std::size_t i = 0; i < once.size(); ++i) {
        worstScale = std::max(worstScale, std::fabs(twice[i] - 2.0 * once[i]));
        scale = std::max(scale, std::fabs(once[i]));
    }
    expectTrue("and it moves the rest of the model", scale > 1e-6);
    expectTrue("linearly, which is what a right-hand side elimination has to be",
               worstScale < 1e-12 * scale);

    // A DOF named twice takes the **last** value rather than the sum. A caller
    // building a prescribed list from two sources -- a punch and a hull girder
    // pulling on the same cut -- will name one twice, and a right-hand side that
    // accumulated would put the DOF somewhere neither source asked for while
    // leaving every other equation looking perfectly reasonable.
    std::vector<double> namedTwice;
    expectTrue("a repeated prescribed DOF solves",
               coupling::prescribedStaticSolve(link.assembly, noLoad, held,
                                               {{target, 1.0e-3}, {target, 2.0e-3}}, namedTwice,
                                               &problem));
    expectNear("and takes the last value, not the sum", namedTwice[target], 2.0e-3, 0.0);
    double worstTwice = 0;
    for (std::size_t i = 0; i < namedTwice.size(); ++i)
        worstTwice = std::max(worstTwice, std::fabs(namedTwice[i] - twice[i]));
    expectTrue("and the whole field is the one that value alone would have given",
               worstTwice < 1e-12 * scale);
    // A prescribed DOF past the end is refused rather than ignored: silently
    // dropping it would leave the model unloaded and plausible.
    std::vector<double> nowhere;
    expectTrue("a prescribed DOF past the end of the assembly is refused",
               !coupling::prescribedStaticSolve(
                   link.assembly, noLoad, held,
                   {{static_cast<std::uint32_t>(link.assembly.size()), 1.0}}, nowhere, &problem));

    // --- The guards, which mutation testing found nothing was reaching ------------
    //
    // Each of these three survived a mutant that removed it, which means the suite
    // was carrying the guard's cost and none of its value.

    // A reduction that is not this substructure's. `assemble` catches the case
    // where component *A*'s map is the wrong size, and does not catch component B
    // being over-large -- every map entry is still in range, so it assembles a
    // model whose zone boundary is longer than the zone has -- and the failure then
    // surfaces two calls later, or not at all.
    {
        const coupling::Coupling wrong = coupling::couple(surroundings, rs, zoneSub, rs);
        expectTrue("a reduction that is not the zone's is refused", !wrong.ready());
        expectTrue("and says which precondition failed", !wrong.problems.empty());
    }

    // A state of the wrong length. Every index in the coupling is in range for the
    // assembly it was built from, so a short state is not a bounds error waiting to
    // happen -- it is a read past the end of the caller's vector, and a plausible
    // drive on the other side of it.
    {
        std::vector<double> state;
        std::vector<coupling::Prescribed> nothing;
        expectTrue("the reference solve runs", coupling::prescribedStaticSolve(
                                                   link.assembly, load, held, nothing, state,
                                                   &problem));
        const coupling::EdgeDrive good = coupling::edgeDrive(link, zoneSub, state);
        expectTrue("a full state drives the edge", good.count > 0);
        std::vector<double> truncated(state.begin(), state.end() - 1);
        const coupling::EdgeDrive bad = coupling::edgeDrive(link, zoneSub, truncated);
        expectEqualCount("a state that is not this assembly's drives nothing", bad.count, 0u);
        expectTrue("and says so", !bad.problems.empty());
    }

    // Mesh surgery, against counts that are known rather than observed.
    const std::size_t elements = f.patch.elementCount();
    std::vector<std::uint8_t> all(elements, 1u);
    const coupling::DamagedMesh nothing = coupling::withoutElements(f.patch.mesh, all);
    expectEqualCount("removing every element leaves no element", nothing.mesh.elementCount(), 0u);
    expectEqualCount("and no node", nothing.mesh.nodeCount(), 0u);
    expectTrue("and says so", !nothing.problems.empty());

    std::vector<std::uint8_t> none(elements, 0u);
    const coupling::DamagedMesh whole = coupling::withoutElements(f.patch.mesh, none);
    expectEqualCount("removing nothing keeps every element", whole.mesh.elementCount(), elements);
    expectEqualCount("and every node", whole.mesh.nodeCount(), f.patch.nodeCount());
    expectTrue("and raises nothing", whole.problems.empty());
    double drift = 0;
    for (std::size_t d = 0; d < whole.mesh.position.size(); ++d)
        drift = std::max(drift, std::fabs(whole.mesh.position[d] - f.patch.mesh.position[d]));
    expectTrue("with every node exactly where it was", drift == 0.0);
    bool sameIndex = whole.mesh.index == f.patch.mesh.index;
    expectTrue("and the same connectivity", sameIndex);

    // A corner element: removing it orphans exactly the nodes only it used. On a
    // patch quad grid that is the one mid-surface corner point, two nodes.
    std::vector<std::uint8_t> corner(elements, 0u);
    std::size_t cornerElement = elements;
    for (std::size_t e = 0; e < elements; ++e) {
        double centre[3];
        elementCentroid(f.patch.mesh, e, centre);
        if (centre[0] < -0.2 && centre[1] < -0.2) cornerElement = e;
    }
    expectTrue("there is a corner element", cornerElement < elements);
    corner[cornerElement] = 1u;
    const coupling::DamagedMesh cut = coupling::withoutElements(f.patch.mesh, corner);
    expectEqualCount("removing a corner element orphans its own corner and nothing else",
                     cut.orphanedNodes, 2u);
    expectEqualCount("the node map is one entry per surviving node", cut.node.size(),
                     f.patch.nodeCount() - 2);
    expectTrue("and every surviving index is in range",
               *std::max_element(cut.mesh.index.begin(), cut.mesh.index.end()) <
                   cut.mesh.nodeCount());
    // An orphan left in would be a zero stiffness row, and `reduction.hpp` §3
    // records that the factorisation does *not* reliably catch one.
    reduction::Substructure trimmed(cut.mesh, steel(),
                                    coupling::carriedInterface(cut, f.zonePerimeter));
    expectTrue("so what is left factors", trimmed.ready());

    // The boundary conditions have to survive the surgery, and the check that says
    // so is the one a caller would actually make: a damaged patch's perimeter is
    // still `nodesPinned` of its own mesh, and that has to be the same set as the
    // perimeter carried across by index. Nothing reached this until a mutant that
    // dropped `fixed` entirely passed the whole suite.
    const std::vector<std::uint32_t> byPin = reduction::nodesPinned(cut.mesh);
    const std::vector<std::uint32_t> byIndex = coupling::carriedInterface(cut, f.zonePerimeter);
    expectTrue("a damaged patch still knows its own perimeter", !byPin.empty());
    expectTrue("and it is the perimeter it had", byPin == byIndex);

    // The prescribed values have to survive the surgery for the same reason as the
    // pins: a damaged patch that was being driven is still being driven, and a
    // drive that quietly became a clamp is the defect this whole file exists to
    // remove.
    std::vector<std::uint32_t> orphans, survivors;
    for (std::uint32_t n : f.zonePerimeter) {
        const auto found = std::lower_bound(cut.node.begin(), cut.node.end(), n);
        if (found != cut.node.end() && *found == n)
            survivors.push_back(n);
        else
            orphans.push_back(n);
    }
    expectEqualCount("the corner element took two perimeter nodes with it", orphans.size(), 2u);
    expectTrue("and left the rest", !survivors.empty());
    {
        HexMesh marked = f.patch.mesh;
        marked.prescribed[3 * survivors.front() + 2] = 1.5e-3;
        const coupling::DamagedMesh carried = coupling::withoutElements(marked, corner);
        const auto where =
            std::lower_bound(carried.node.begin(), carried.node.end(), survivors.front());
        expectTrue("the driven node survived", where != carried.node.end() &&
                                                   *where == survivors.front());
        expectNear("and is still driven to the same place",
                   carried.mesh.prescribed[3 * static_cast<std::size_t>(
                                                   where - carried.node.begin()) +
                                           2],
                   1.5e-3, 0.0);
    }

    // A node that is gone must be **dropped**, not mapped to whichever node now
    // sits at its place in the ordering. Asked with the orphan on its own, because
    // asked as part of the whole perimeter the wrong answer lands on a node that is
    // in the list anyway and the duplicate is collapsed -- which is exactly why a
    // mutant that dropped the equality test survived the first version of this.
    const std::vector<std::uint32_t> gone =
        coupling::carriedInterface(cut, {orphans.front()});
    expectTrue("an interface node the tear removed is dropped, not mapped to its neighbour",
               gone.empty());

    // `applyEdgeDrive` writes what the drive names and **never un-pins**: a
    // perimeter DOF the surroundings do not share stays clamped rather than going
    // quietly free, which is the difference between a partly matched interface and
    // a hole in the boundary condition.
    {
        HexMesh partly = f.patch.mesh;  // buildPatch pinned the whole perimeter
        const std::size_t before = static_cast<std::size_t>(
            std::count(partly.fixed.begin(), partly.fixed.end(), std::uint8_t{1}));
        expectTrue("the patch arrives with its perimeter pinned", before > 3u);
        coupling::EdgeDrive one;
        one.displacement.assign(partly.nodeCount() * 3, 0.0);
        one.driven.assign(partly.nodeCount() * 3, 0u);
        one.driven[3 * survivors.front()] = 1u;
        one.displacement[3 * survivors.front()] = 1.0e-4;
        one.count = 1;
        expectEqualCount("one DOF is written", coupling::applyEdgeDrive(one, partly), 1u);
        expectEqualCount("and nothing that was pinned is now free",
                         static_cast<std::size_t>(std::count(partly.fixed.begin(),
                                                             partly.fixed.end(),
                                                             std::uint8_t{1})),
                         before);
    }
}

// --- 6. What Tier 1 cannot see: the stiffeners -----------------------------------
//
// A `reduction::Substructure` is built from a `HexMesh` and a material. A
// `zone::Patch` under `Stiffeners::Modelled` carries its members as
// `constraint::Stiffening` -- fibres condensed onto the plating's own DOF, with no
// nodes of their own and no elements. **They are therefore invisible to Tier 1**,
// and a stiffened zone is handed to the surrounding model as bare plating.
//
// This is measured rather than asserted in prose, because nothing tests a comment
// and because the size of the error is the whole question: if the fibres were
// worth a per cent it would be a footnote.
void testTierOneCannotSeeTheStiffeners() {
    std::printf("\n--- coupling: the stiffeners a reduced model does not carry ---\n");
    StructuralMesh stiffened = flatPlate();
    StructuralMember member;
    member.a = {-0.5 * kSide, 0.0, 0.0};
    member.b = {0.5 * kSide, 0.0, 0.0};
    member.rise = {0, 0, 1};
    member.profile = flatBar(0.200, 0.010);
    member.attachedPlateThickness = kThickness;
    member.role = MemberRole::Longitudinal;
    stiffened.members.push_back(member);

    zone::MeshParams params = patchParams();
    params.stiffeners = zone::Stiffeners::Modelled;
    const zone::Patch patch = zone::buildPatch(stiffened, {0, 0, 0}, params);
    expectTrue("the patch carries the member as fibres", patch.stiffening.fiberCount() > 0);

    // The same mesh and the same load, solved with and without the fibres'
    // condensed stiffness. `solveStatic` takes exactly the `DofBlock` list
    // `constraint::stiffnessBlocks` produces -- and `reduction::Substructure` takes
    // no such list, which is the gap.
    const constraint::RestFibers forms =
        constraint::restFibers(patch.stiffening, patch.mesh.position);
    expectTrue("the fibres have rest lengths", forms.ok);
    const std::vector<solidshell::DofBlock> blocks = constraint::stiffnessBlocks(
        patch.stiffening, patch.mesh.position, forms, steel().youngsModulus);
    expectEqualCount("one condensed block per fibre", blocks.size(),
                     patch.stiffening.fiberCount());

    std::vector<std::size_t> punch;
    for (std::size_t n = 0; n < patch.nodeCount(); ++n)
        if (isPunchPoint(patch.mesh.position[3 * n], patch.mesh.position[3 * n + 1]))
            punch.push_back(n);
    expectTrue("the punch landed on the patch", !punch.empty());

    HexMesh mesh = patch.mesh;
    prescribePunch(mesh, punch);
    std::vector<double> withFibres, withoutFibres;
    std::string problem;
    expectTrue("the stiffened patch solves",
               solidshell::solveStatic(mesh, steel(), solidshell::Formulation::SolidShell, blocks,
                                       {}, withFibres, &problem));
    expectTrue("and the bare one",
               solidshell::solveStatic(mesh, steel(), solidshell::Formulation::SolidShell, {},
                                       withoutFibres, &problem));
    const double stiff = punchReaction(mesh, steel(), withFibres, punch).force;
    const double bare = punchReaction(mesh, steel(), withoutFibres, punch).force;
    // `punchReaction` assembles the elements alone, so for the stiffened field it
    // reports only the plating's share; the fibres' share is the difference in the
    // *fields*, which is what the ratio below is taken on instead.
    double worst = 0, peak = 0;
    for (std::size_t d = 0; d < withFibres.size(); ++d) {
        worst = std::max(worst, std::fabs(withFibres[d] - withoutFibres[d]));
        peak = std::max(peak, std::fabs(withoutFibres[d]));
    }
    std::printf("     a 200x10 bar changes the patch's field by %.1f%% (plating reaction "
                "%.4g N against %.4g N)\n",
                100.0 * worst / peak, stiff, bare);
    expectTrue("the fibres are worth a great deal more than a rounding", worst > 0.05 * peak);

    // And the substructure a coupling would build from the same patch carries none
    // of it: its operator is the element assembly and nothing else, to the last
    // bit. That is the measurement, not the argument.
    reduction::Substructure zoneSub(patch.mesh, steel(), reduction::nodesPinned(patch.mesh));
    expectTrue("the substructure factors", zoneSub.ready());
    std::vector<double> probe(zoneSub.dofCount(), 0.0), got;
    for (std::size_t d = 0; d < probe.size(); ++d)
        probe[d] = 1e-6 * static_cast<double>((d * 37u) % 19u) - 9e-6;
    zoneSub.stiffnessTimes(probe, got);
    const std::vector<double> elementsOnly = nodalForce(patch.mesh, steel(), probe);
    double gap = 0, scale = 0;
    for (std::size_t d = 0; d < got.size(); ++d) {
        gap = std::max(gap, std::fabs(got[d] - elementsOnly[d]));
        scale = std::max(scale, std::fabs(elementsOnly[d]));
    }
    expectTrue("the probe loads the operator", scale > 1.0);
    // 1e-12, because the measurement is 3e-16 relative: the two are the *same*
    // assembly, not two well-agreeing ones, and a loose tolerance here would pass
    // on a substructure that carried a little of the stiffener.
    expectTrue("and a Tier-1 substructure of a stiffened patch is the bare plating",
               gap < 1e-12 * scale);
    std::printf("     and the substructure a coupling builds from it is the plating alone, "
                "to %.1e of %.3e N\n",
                gap, scale);
}

}  // namespace

void runCouplingTests() {
    std::printf("\n=== Tier-1 to Tier-2 coupling ===\n");
    testCoupledZoneIsTheWholePlate();
    testDamageIsFeltOutsideTheZone();
    testAThinnerZoneNoLongerFits();
    testTheExplicitSolverFollowsItsBoundary();
    testTheRampAndTheTearReadBack();
    testTheWholeChainAgainstTheWholePlate();
    testPrescribedSolveAndMeshSurgery();
    testTierOneCannotSeeTheStiffeners();
}
