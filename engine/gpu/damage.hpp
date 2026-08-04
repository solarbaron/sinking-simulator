// SPDX-License-Identifier: MIT
//
// **Structural damage, made into geometry.** `engine/sim/zone.{hpp,cpp}` solves
// solid-shell elements over a real patch of plating and reports displaced node
// positions and which panels tore; `engine/sim/indentation.{hpp,cpp}` reports a
// depth and a torn set for microseconds instead of core-minutes. Neither of them
// reached the screen. This turns either into a `sim::TriMesh` the lit pass already
// knows how to draw: a hull that has been pushed in, with holes in it.
//
// It contains no Vulkan, deliberately -- it is mesh arithmetic, so it builds and
// is testable on a machine with no GPU at all, the same argument that keeps
// `material.cpp` out of the Vulkan-gated sources. `gpu::SceneMesh::appendShip`
// takes the result.
//
// --- 1. Joining a deformed patch to an undeformed hull -------------------------
//
// The zone deforms a few square metres of a 120 m ship. Everything about making
// that look like one hull rather than a patch stuck onto one comes down to two
// separate hazards, and they have separate answers:
//
//   * **A crack** -- two triangles that shared an edge no longer agreeing where
//     that edge is. Prevented by construction: the displacement is a **pure
//     function of the undeformed position**, `displacementAt()`. Two triangles
//     sharing a position -- by index, or as two coincident-but-distinct vertex
//     records on a mirrored seam -- are handed the same position and get back the
//     same answer, so nothing can pull them apart. Nothing is displaced per
//     triangle and no vertex is added or removed by the deformation.
//
//   * **A seam** -- a step in the surface where the deformed region stops.
//     Prevented by the **boundary condition the solve already has**:
//     `zone::Edge::Clamped` pins the patch perimeter, so the outermost nodes carry
//     displacement exactly zero and the field reaches zero *inside* its own
//     support rather than being cut off at it. `HullDamage::boundaryDisplacement()`
//     reports the largest displacement on the patch's boundary loop, so this is
//     measured rather than assumed -- a patch solved with a free edge says so.
//
// **The field is the FEM's own, not a fit to it.** `displacementAt` is piecewise
// linear over the patch's own elements, each split into two triangles, evaluated
// by barycentric weights. That makes it *interpolating*: at a mid-surface node the
// weights are exactly (1, 0, 0) and the answer is exactly that node's displacement,
// which is the identity `tests/test_hull_render.cpp` asserts to 1e-9 m over every
// node in the patch. Inverse-distance (Shepard) weighting was written first and
// rejected: it is only interpolating in the limit, it puts a flat spot at every
// data point, and it turns "does the picture agree with the solver" into a
// tolerance instead of an equality.
//
// Splitting a quad into triangles and interpolating linearly is C0 across every
// seam it can produce -- along a shared element edge both sides are linear in the
// same edge parameter through the same two nodal values, and the diagonal is
// shared by the two triangles of one element -- so the field is continuous
// everywhere and identically zero off the patch.
//
// --- 2. Resolution: a dent needs somewhere to happen ---------------------------
//
// The reference ferry's hull is tabulated at 25 stations and 12 waterlines, so its
// triangles are about 5 m by 1.5 m. A 0.2 m dent over a 2 m punch drawn on that
// mesh moves nothing: there is no vertex inside the dent to move. The plating has
// to be refined where the damage is, and refining part of a mesh is exactly the
// problem the ocean cascade solves at a level boundary -- **a T-junction is a
// vertex halfway along a neighbour's edge that the neighbour interpolates straight
// past**, and it leaves a wedge of background showing through.
//
// The answer here is the same shape as the cascade's: make the two sides of every
// edge decide identically.
//
//   * Whether an edge splits, `needsSplit`, is a **pure function of its two
//     endpoints** -- their distance against a target size that depends only on
//     position. The two triangles sharing an edge therefore cannot disagree.
//   * The midpoint is `(a + b) * 0.5` on the *welded* endpoint positions, which is
//     commutative in IEEE arithmetic and so is bit-identical from both sides, and
//     it is interned so both sides address the same vertex rather than two vertices
//     that happen to coincide.
//   * A triangle with one or two split edges is divided into two or three -- the
//     "green" cases -- so it **uses the midpoint too**. That is the cascade's
//     transition cell, in a triangle mesh.
//
// The recursion is uniform: every child is examined at depth + 1 whatever template
// produced it, so a depth cap cuts both sides of a shared edge at the same depth.
// `HullDamageParams::stitch` exists only to switch the green cases off, which
// reintroduces the T-junctions on purpose -- it is the negative control the crack
// tests need, and the edge census counts the difference.
//
// **Cost is zero when there is no damage.** The target size away from the damage
// defaults to "never split", so an undamaged hull is walked once and emitted
// unchanged -- and unchanged means *bit-identical*, which is asserted rather than
// hoped for. The one place that needs care is that `-0.0 + 0.0` is `+0.0`, so a
// zero displacement must not go through the addition at all.
//
// --- 3. A tear is a hole ------------------------------------------------------
//
// A torn panel is removed from the drawn mesh. It is not recoloured, not made
// transparent, and not moved: the triangles are gone, so the pixel behind them
// shows whatever is behind them -- the far side of the hull from the inside, the
// compartment structure if any is drawn, or the sea. The hull's rasterisation
// already permits this: `cullMode` is `NONE` and the fragment shader flips the
// normal on a back face, both of which exist because "a hull is seen from inside
// once it is cut away", and this is that case arriving.
//
// What it costs: the drawn hull stops being a closed surface, so the interior is
// shaded and filled where it used to be free. Nothing else changes -- no second
// pass, no sort, no stencil, no discard in the fragment shader, and no extra
// vertex attribute. `sim::Ship::hull` is untouched, so buoyancy still integrates a
// closed mesh; the hole exists in the render mesh only, which is the same
// separation `breach.hpp` keeps between a hole and the flooding network.
//
// The hole's edge is only as sharp as the refinement: a triangle is dropped when
// its centroid falls inside a torn panel, so the boundary is accurate to about
// `fineSize`. That is why a torn panel is also a refinement feature.
//
// --- 4. Exposed metal is a material, not a shader branch -----------------------
//
// `docs/03-renderer-audio.md` asks for exposed-metal blending along a torn edge.
// The material system is data (`engine/gpu/materials/marine.materials`), so this
// is a *material name* on the triangles within `exposedWidth` of a hole, resolved
// against the library like every other paint band. A mod renames it or redefines
// it without recompiling anything.
//
// **It is a band and not a blend, and that is a real limit.** `HullVertex` carries
// one material index, so two surfaces cannot be mixed per pixel; a plastic-strain
// driven blend needs a second index and a weight on the vertex and a change to
// `shaders/hull.frag`. Recorded in `docs/03-renderer-audio.md` rather than
// pretended away.
#pragma once

#include "../core/geometry.hpp"
#include "../core/math.hpp"
#include "../sim/scantlings.hpp"
#include "../sim/zone.hpp"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace gpu {

struct HullDamageParams {
    // Target edge length inside the damaged region, m. Everything about whether a
    // dent is visible is this number against the size of the dent: the reference
    // ferry's own plating is drawn at about 5 m by 1.5 m, which cannot show a
    // 0.2 m dent because there is no vertex inside it.
    double fineSize = 0.20;
    // Target edge length away from it. The default is larger than any ship, which
    // means "never split", which is what makes an undamaged hull cost nothing and
    // a damaged one cost only where it is damaged.
    double coarseSize = 1e30;
    // How fast the target grows with distance from the damage. Below one the mesh
    // grades out instead of stepping, so no single edge halves more than once per
    // level and the green cases stay small.
    double grading = 0.6;
    // Recursion guard. Reaching it is reported (`DamagedHull::depthLimited`)
    // rather than silently accepted, and it is applied at the same depth on both
    // sides of every edge, so hitting it cannot itself open a crack.
    int maxDepth = 8;

    // **The through-thickness rule, and it is one rule for all three sources.**
    // How far off a torn panel's plane -- or a patch's, or a tent's -- a hull
    // triangle may be and still be part of it, m.
    //
    // It is not optional and it is not a tolerance. A hull *wraps*: a point on her
    // port shell projects into the same in-plane cell as the point opposite it on
    // the starboard shell, so without this, ramming one side dents, cuts and
    // exposes the other identically. Panel corners lie on the plate mid-surface,
    // which *is* the hull mesh surface (`scantlings.hpp`), so on the struck side
    // the two disagree only by the sagitta of the hull's curvature across one
    // panel -- centimetres. Half a metre is generous for that and far less than a
    // ship's breadth.
    //
    // Raise it for a field whose own axis is not the plating's normal: the guard
    // is measured along that axis, so a displacement directed *along* the shell
    // needs a slab deep enough to contain the plating it is meant to move.
    double panelReach = 0.5;
    // Width of the exposed-metal band outside a hole, m.
    double exposedWidth = 0.25;

    // **A negative control, not a setting.** False emits a triangle whole unless
    // all three of its edges want splitting, which leaves a T-junction against
    // every neighbour that did split -- exactly the crack the green cases exist to
    // remove. The tests build both and count the difference; nothing in the engine
    // sets it false.
    bool stitch = true;
};

// A dent and a set of holes, in the ship's **body frame** -- the frame paint bands
// are decided in, the frame `sim::StructuralMesh` is built in, and the frame the
// zone meshes in. The rigid-body transform is applied once, later, by
// `SceneMesh::appendShip`, so this survives the ship moving.
class HullDamage {
public:
    HullDamageParams params;

    // The zone's answer: the patch, and the solver's node positions
    // (`zone::Solver::position()`, 3 per node in the patch's own numbering).
    // False with a reason when the two do not describe the same patch.
    //
    // A patch contributes a piecewise-linear displacement over its own elements
    // and a refinement feature covering its extent.
    bool addZone(const sim::zone::Patch& patch, const std::vector<double>& deformed,
                 std::string& error);

    // The membrane model's answer. `indentation.hpp` reports a penetration and a
    // torn set and **not a shape**, so the shape is stated here rather than
    // implied: it is that model's own kinematics -- a plate pinned at `halfSpan`
    // and pushed to `depth` at the centre, linear between, which is the tent the
    // closed forms in `indentation.hpp` are integrals of. `outward` is the
    // direction the plating is pushed *from*; the dent goes the other way.
    //
    // Its support is a **ball** of radius `halfSpan` about `centre`, so the field
    // reaches zero at the edge of it and is continuous everywhere, with Lipschitz
    // constant exactly `depth / halfSpan` -- which is what a test can hold every
    // edge of the drawn mesh to. Keep `halfSpan` under the ship's half breadth, or
    // a dent on one side reaches the plating on the other; there is no separate
    // guard, because the support already is one.
    void addTent(const sim::Vec3& centre, const sim::Vec3& outward, double halfSpan,
                 double depth);

    // Panels that tore, as indices into `structure.panels`. Out-of-range indices
    // are ignored rather than fatal, matching `breachesFromFailedPanels`.
    void addTornPanels(const sim::StructuralMesh& structure, const std::vector<int>& panels);

    bool empty() const { return zones_.empty() && tents_.empty() && torn_.empty(); }
    std::size_t tornPanelCount() const { return torn_.size(); }
    // m^2 of torn panel, as the structural mesh measures it. The area actually
    // removed from the drawn hull is `DamagedHull::holeArea`, and the two agree to
    // the refinement.
    double tornPanelArea() const { return tornArea_; }

    // The displacement at an undeformed body-frame position. **A pure function of
    // its argument**, which is what makes the join crack-free; exactly zero off
    // every patch and every tent.
    sim::Vec3 displacementAt(const sim::Vec3& rest) const;

    // Largest displacement over every mid-surface node of every patch, and largest
    // over the nodes on a patch's **boundary loop**. The second is zero for a
    // clamped edge, which is the statement that the field reaches zero before its
    // support ends and therefore joins the undeformed hull without a step. Both are
    // reported so a test can assert the first is non-trivial and the second is not.
    double largestNodeDisplacement() const;
    double boundaryDisplacement() const;

    // Target edge length at a position: `fineSize` inside the damage, growing at
    // `grading` outside it, capped at `coarseSize`.
    double targetEdgeSize(const sim::Vec3& at) const;

    // True when `at` is inside a torn panel -- the test that removes a triangle.
    bool insideTear(const sim::Vec3& at) const;
    // Distance from `at` to the nearest torn panel's outline, m, or a large number
    // when there is none within reach of the panel's own plane.
    double distanceToTearEdge(const sim::Vec3& at) const;

private:
    // One meshed patch, flattened into its own in-plane frame so a query is a
    // point-in-triangle test rather than a search.
    struct Zone {
        sim::Vec3 centre{}, right{}, up{}, axis{};
        std::vector<sim::Vec3> rest;          // per mid-surface node
        std::vector<sim::Vec3> displacement;  // per mid-surface node
        std::vector<double> u, v;             // in-plane coordinates of `rest`
        std::vector<std::uint32_t> tri;       // 3 per triangle, into the above
        // Uniform bucket grid over the (u, v) bounding box.
        double uLo = 0, vLo = 0, cell = 1;
        int nu = 1, nv = 1;
        std::vector<std::uint32_t> bucketStart, bucketItem;
        double radius = 0;     // in-plane radius about `centre`, for refinement
        double largest = 0, boundary = 0;
    };
    struct Tent {
        sim::Vec3 centre{}, outward{};
        double halfSpan = 0, depth = 0;
    };
    struct Quad {
        sim::Vec3 corner[4]{};
        sim::Vec3 normal{}, centroid{}, right{}, up{};
        double radius = 0;
        double u[4]{}, v[4]{};
    };

    std::vector<Zone> zones_;
    std::vector<Tent> tents_;
    std::vector<Quad> torn_;
    double tornArea_ = 0;
};

// The hull as it is drawn once it has been damaged: refined where the damage is,
// displaced by the field, and with the torn panels removed.
//
// `rest` and `deformed` have **identical topology** and differ only in vertex
// position. Both are kept because they answer different questions: paint bands are
// decided on the undeformed hull -- a dented plate keeps the paint it was painted
// with, and deciding the boot-topping from a displaced z would slide the stripe
// around as the plating moved -- while the geometry drawn and the normals lit are
// the deformed one.
struct DamagedHull {
    sim::TriMesh rest;
    sim::TriMesh deformed;
    // Per triangle of either mesh: 1 where it is within `exposedWidth` of a hole
    // and takes the exposed-metal material instead of its paint band.
    std::vector<std::uint8_t> exposed;
    // The triangles that were cut out, indexing the same vertex arrays. Kept
    // because the hole's boundary is only computable from them: the census of what
    // is left cannot tell an edge the tear exposed from an edge a bad refinement
    // opened, and separating those two is the whole question. A caller with debris
    // or a cutaway view has them for free.
    std::vector<sim::Tri> removed;

    std::size_t sourceTriangles = 0;   // triangles of the hull it came from
    std::size_t droppedTriangles = 0;  // removed because they fell inside a tear
    int deepestSplit = 0;
    bool depthLimited = false;         // `maxDepth` was reached somewhere
    double holeArea = 0;               // m^2 of undeformed surface removed
    double largestDisplacement = 0;    // m
    double buildSeconds = 0;

    std::size_t triangleCount() const { return deformed.tris.size(); }
};

// Refine, displace and cut. A rebuild, not a per-frame step: the result is in the
// body frame, so a ship that has moved needs the transform reapplied (which
// `SceneMesh::appendShip` does every frame anyway) and not this.
//
// With an empty `damage` the result is the input, **bit-identical** -- same vertex
// array, same triangle array, no triangle dropped and none exposed.
DamagedHull buildDamagedHull(const sim::TriMesh& hull, const HullDamage& damage);

}  // namespace gpu
