// SPDX-License-Identifier: MIT
//
// The Phase 3 milestone, in one act: **ram the ferry.**
//
// `docs/06-roadmap.md` states it as "the hull deforms, tears where the stress
// says it should, and the resulting hole floods at a rate the hole's own area
// determines." Every piece of that exists and each was built and validated
// separately. Nothing had ever run them as one event, which is where the
// interesting failures live.
//
// The chain, and not a step of it is reimplemented here:
//
//   collision   two closed hulls overlap; the penetration volume gives a force,
//               a patch, and the energy that went into the meeting
//   indentation that energy spent outward panel by panel; membrane stretching
//               says how deep it got and which bays tore
//   breach      torn panels become openings in the flooding network, merged,
//               placed at their own centroid, connected by geometry
//   ship        the openings flood at Cd A sqrt(2 dp / rho), the free surfaces
//               re-level, the trapped air compresses, and she lists or does not
//
// The striking ship is rigid here. She is not: a real bow crushes and takes a
// large share of the energy, and `indentation.hpp` records what leaving that out
// costs. This over-states the damage, knowingly, and the run prints enough to see
// by how much.
//
// **No damage control is applied.** Nobody closes a door, starts a pump or
// counterfloods, so this is Phase 0's `none` scenario with the breach computed
// rather than authored -- and `none` loses her too. She is lost at every speed and
// every aiming point tried, which is consistent rather than suspicious: the
// scenarios that let her live are the ones where somebody acts. What does respond
// is *how* she is lost. Struck amidships she takes 2100 t and lolls 6 degrees;
// struck at the quarter she takes 8100 t and goes over 26 the other way.
//
// It also shows where the outcome stops caring about the strike. From 1.5 to
// 6 m/s the hole grows from 3.4 to 113 m2 and the floodwater barely moves,
// because a 3.4 m2 breach fills the compartment behind it inside 900 s just as a
// 113 m2 one does. Beyond a threshold, damage stability is decided by *which*
// compartments are opened and not by how big the hole is -- which is exactly what
// the subdivision rules are written around, and it falls out here rather than
// being put in.
//
// **It draws, now.** `--frames=N --out=DIR` writes the flooding sequence with the
// damage in it: the struck plating dished in by the penetration the membrane model
// reports, the torn bays cut out as holes you see the compartment through, and
// exposed metal round their edges. The dent's *shape* is the membrane model's own
// tent kinematics -- `indentation.hpp` reports a depth and a torn set and not a
// surface -- and `zone.hpp` is what replaces the assumption when a run can afford
// it; `gpu::HullDamage::addZone` takes its displaced nodes instead.
//
// Without `--frames` nothing GPU happens at all, which is why the gate can run
// this on a machine with no device.
//
//   ./ram_view [--speed=M_PER_S] [--aim=X_METRES] [--reach=METRES]
//              [--duration=SECONDS] [--out=DIR] [--frames=N]
#include "engine/core/math.hpp"
#include "engine/core/png.hpp"
#include "engine/gpu/damage.hpp"
#include "engine/gpu/hull.hpp"
#include "engine/gpu/ocean.hpp"
#include "engine/sim/breach.hpp"
#include "engine/sim/collision.hpp"
#include "engine/sim/hullform.hpp"
#include "engine/sim/indentation.hpp"
#include "engine/sim/scantlings.hpp"
#include "engine/sim/waves.hpp"
#include "game/prototype/ferry.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

namespace {

struct Options {
    double speed = 6.0;       // m/s closing
    double aim = 0.0;         // m along the ferry, where the bow lands
    double reach = 12.0;      // m, how far damage may propagate from the impact
    double duration = 900.0;  // s of flooding after the strike
    std::string out;
    int frames = 0;
};

bool parse(int argc, char** argv, Options& o) {
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        const auto value = [&](const char* key) -> const char* {
            const std::string prefix = std::string("--") + key + "=";
            return a.rfind(prefix, 0) == 0 ? a.c_str() + prefix.size() : nullptr;
        };
        if (const char* v = value("speed")) o.speed = std::atof(v);
        else if (const char* v = value("aim")) o.aim = std::atof(v);
        else if (const char* v = value("reach")) o.reach = std::atof(v);
        else if (const char* v = value("duration")) o.duration = std::atof(v);
        else if (const char* v = value("out")) o.out = v;
        else if (const char* v = value("frames")) o.frames = std::atoi(v);
        else {
            std::printf("unknown argument: %s\n", a.c_str());
            return false;
        }
    }
    return true;
}

// A striking ship: a smaller hull, built from particulars like any other.
sim::Ship striker() {
    sim::HullParticulars p;
    p.lengthPp = 90.0;
    p.beam = 15.0;
    p.draft = 5.5;
    p.depth = 11.0;
    p.blockCoefficient = 0.68;
    p.midshipCoefficient = 0.97;
    p.parallelMiddleBodyFraction = 0.25;
    p.stationCount = 41;

    sim::Ship s;
    s.hull = sim::makeHullFromParticulars(p);
    s.deckEdgeZ = p.depth;
    s.lightshipMass = p.blockCoefficient * p.lengthPp * p.beam * p.draft * sim::kRhoSeawater;
    s.lightshipCog = {0.0, 0.0, 0.55 * p.depth};
    s.gyradii = {0.35 * p.beam, 0.25 * p.lengthPp, 0.25 * p.lengthPp};
    return s;
}

}  // namespace

int main(int argc, char** argv) {
    Options options;
    if (!parse(argc, argv, options)) return 2;

    sim::Ship ferry = game::buildFerry();
    ferry.initialise(0.0);
    const sim::Scantlings scantlings = sim::ferryScantlings();
    const sim::StructuralMesh structure = sim::makeStructuralMesh(ferry.hull, scantlings);

    sim::Ship bow = striker();
    bow.initialise(0.0);

    const sim::Diagnostics before = ferry.diagnostics(0.0);
    std::printf("=== ram the ferry ===\n");
    std::printf("struck : 120 m ro-pax, %.0f t, GM %.2f m, %zu structural panels\n",
                before.displacementMass / 1000.0, before.gmTransverse, structure.panels.size());
    std::printf("striker: 90 m hull, %.0f t, closing at %.1f m/s, aimed at x = %+.0f m\n",
                bow.diagnostics(0.0).displacementMass / 1000.0, options.speed, options.aim);

    // Place the striker abeam, on the ferry's starboard side, closing.
    const double standoff = 0.5 * (ferry.hullHi.y - ferry.hullLo.y) + 45.0;
    bow.state.position = {options.aim, -standoff, bow.state.position.z};
    bow.state.orientation = sim::Quat::fromAxisAngle({0, 0, 1}, sim::kPi / 2.0);
    bow.state.velocity = {0.0, options.speed, 0.0};

    // --- 1. The strike -------------------------------------------------------
    const sim::Sea sea(0.0);
    const double dt = 0.01;
    sim::ContactMaterial material;
    sim::ContactHistory history;

    bool touched = false;
    double peakForce = 0;
    for (int i = 0; i < 4000 && (!touched || history.duration > 0); ++i) {
        const sim::HullContact contact = sim::applyContact(ferry, bow, material, dt, &history);
        if (contact.touching) {
            touched = true;
            const sim::ContactBody a = sim::contactBodyOf(ferry), b = sim::contactBodyOf(bow);
            peakForce = std::max(peakForce, length(sim::contactLoad(contact, a, b, material).force));
        } else if (touched) {
            break;   // through and clear
        }
        ferry.step(dt, sea);
        bow.step(dt, sea);
    }

    if (!touched) {
        std::printf("\nno contact -- she was missed\n");
        return 0;
    }
    std::printf("\ncontact: %.2f s, peak %.0f MN, %.0f MJ taken out of the two ships\n",
                history.duration, peakForce / 1e6, history.work / 1e6);

    // --- 2. What that energy did to the plating ------------------------------
    const sim::Vec3 impact = history.loadAtPeak.point - ferry.state.position;
    const sim::ImpactDamage damage =
        sim::impactDamage(structure, impact, options.reach, history.work, scantlings);
    std::printf("damage : %.3f m into her over %zu bays, %zu torn, %.1f m2 of hole\n",
                damage.penetration, damage.panels.size(), damage.torn.size(), damage.tornArea);
    if (damage.energyUnspent > 0)
        std::printf("         %.0f MJ unspent -- reach bounded the answer, the hull did not\n",
                    damage.energyUnspent / 1e6);

    // --- 3. The hole -------------------------------------------------------
    const sim::BreachSet breaches = sim::breachesFromFailedPanels(ferry, structure, damage.torn);
    double openArea = 0;
    for (const sim::Breach& b : breaches.breaches) openArea += b.opening.area;
    std::printf("breach : %zu opening(s), %.2f m2 reaching a compartment\n",
                breaches.breaches.size(), openArea);
    for (const std::string& problem : breaches.problems)
        std::printf("         ! %s\n", problem.c_str());
    sim::applyBreaches(ferry, breaches);

    // --- 4. What that looks like ---------------------------------------------
    //
    // The damage as geometry, in the ferry's own body frame, so it survives her
    // heeling over. Two things go into it and they come from different places:
    // the torn bays are exactly the panels `impactDamage` reported, cut out; the
    // dent is a *shape*, which that model does not report, so it is stated as the
    // tent its own closed forms are integrals of, over the equivalent radius of
    // the bays it actually struck.
    gpu::HullDamage visible;
    // 0.35 m rather than the 0.20 m default: at 4 m/s she opens 60-odd bays over
    // a hundred square metres, and refining all of that to 0.20 m is four times
    // the triangles for a hole whose edge is a panel boundary anyway.
    visible.params.fineSize = 0.35;
    visible.addTornPanels(structure, damage.torn);
    double struckArea = 0;
    for (int index : damage.panels)
        if (index >= 0 && static_cast<std::size_t>(index) < structure.panels.size())
            struckArea += structure.panels[static_cast<std::size_t>(index)].area();
    const double dentRadius = std::max(std::sqrt(struckArea / sim::kPi), 1.0);
    // She is struck on the starboard side, which is -y in the body frame.
    visible.addTent(impact, {0.0, -1.0, 0.0}, dentRadius, damage.penetration);
    const gpu::DamagedHull damagedHull = gpu::buildDamagedHull(ferry.hull, visible);
    std::printf("drawn  : %.2f m dent over a %.1f m radius; her %zu hull triangles refine to"
                " %zu, of which %zu are cut out (%.1f m2 of hole), in %.1f ms\n",
                damagedHull.largestDisplacement, dentRadius, damagedHull.sourceTriangles,
                damagedHull.deformed.tris.size() + damagedHull.droppedTriangles,
                damagedHull.droppedTriangles, damagedHull.holeArea,
                damagedHull.buildSeconds * 1e3);

    gpu::Device device;
    gpu::MaterialLibrary library;
    gpu::HullRenderer renderer;
    constexpr std::uint32_t kWidth = 1280, kHeight = 720;
    std::string error;
    bool drawing = options.frames > 0 && !options.out.empty();
    if (drawing && !device.create(error)) {
        std::printf("       no usable GPU (%s) -- not drawing\n", error.c_str());
        drawing = false;
    }
    if (drawing &&
        !library.load(std::string(SHIPSIM_MATERIAL_DIR) + "/marine.materials", error)) {
        std::printf("materials: %s\n", error.c_str());
        return 1;
    }
    if (drawing && !renderer.create(device, kWidth, kHeight, SHIPSIM_SHADER_DIR, error)) {
        std::printf("renderer: %s\n", error.c_str());
        return 1;
    }

    sim::SeaState calm;
    calm.significantHeight = 0.0;
    const sim::WaveField field(calm);
    gpu::OceanSurface water;
    gpu::SceneMesh scene;
    gpu::HullPaint paint;
    paint.waterlineZ = 5.5;
    paint.deckZ = 15.0;
    gpu::SceneView view;
    // The sun over the struck side. The default is on her port bow, which is the
    // side the camera is not on, and a hole lit only by the sky term is a black
    // patch in a black patch.
    view.sunDirection[0] = -0.32f;
    view.sunDirection[1] = -0.58f;
    view.sunDirection[2] = 0.75f;
    gpu::HullShading shading;
    shading.normals = gpu::HullNormals::Smooth;

    // The compartments behind the shell, so a hole is a hole into *something*. They
    // are the same meshes the flooding solves against -- `Compartment::mesh` -- not
    // a set of boxes drawn to look like an interior, and they are drawn flat rather
    // than smoothed because a compartment is a box and its corners are corners.
    //
    // They are **inset** before being drawn, and the reason is worth recording: a
    // compartment is `clipToBox(hull, ...)`, so its outboard face is not merely
    // near the shell, it *is* the shell -- the same surface, to the last bit. Drawn
    // as they are, the two z-fight over the whole side and the ship comes out
    // speckled. Three per cent about each compartment's own centroid is a
    // rendering-only inset; nothing the sim does uses these copies.
    sim::TriMesh interior;
    for (const sim::Compartment& compartment : ferry.compartments) {
        if (compartment.mesh.verts.empty()) continue;
        sim::TriMesh inset = compartment.mesh;
        sim::Vec3 centroid{0, 0, 0};
        for (const sim::Vec3& v : inset.verts) centroid += v;
        centroid = centroid / static_cast<double>(inset.verts.size());
        for (sim::Vec3& v : inset.verts) v = centroid + (v - centroid) * 0.97;
        interior.append(inset);
    }

    // --- 5. Whether she lives ------------------------------------------------
    std::printf("\n%8s %10s %8s %8s %8s\n", "t (s)", "flood t", "heel", "trim", "GM");
    const auto steps = static_cast<int>(options.duration / dt);
    const int framesEvery = drawing ? std::max(1, steps / options.frames) : 0;
    int frame = 0;
    for (int i = 1; i <= steps; ++i) {
        ferry.step(dt, sea);
        if (i % (steps / 6) == 0 || i == steps) {
            const sim::Diagnostics d = ferry.diagnostics(sea);
            std::printf("%8.0f %10.0f %7.2f° %7.2f° %8.2f\n", i * dt,
                        d.floodwaterMass / 1000.0, d.heelDeg, d.trimDeg, d.gmTransverse);
        }
        if (!drawing || frame >= options.frames || i % framesEvery != 0) continue;

        // On the struck bow quarter and well above her, in the ship's own frame, so
        // the hole stays in shot as she lists away from it.
        const sim::Mat3 R = ferry.state.orientation.toMat3();
        const sim::Vec3 eye =
            ferry.state.position + R * sim::Vec3{impact.x + 33.0, -46.0, 19.0};
        const sim::Vec3 at = ferry.state.position + R * sim::Vec3{impact.x, -8.0, 6.0};
        water.build(field, gpu::OceanGrid{ferry.state.position.x, ferry.state.position.y, 600.0,
                                          96},
                    0.0);
        scene.clear();
        // The interior first, the shell over it, the sea last: the depth test then
        // carries every one of those occlusions rather than the draw order, which is
        // the arrangement the render tests use deliberately.
        scene.appendMesh(interior, R, ferry.state.position,
                         static_cast<std::uint32_t>(library.find("rusted_steel")),
                         gpu::HullShading{});
        if (!scene.appendShip(ferry, damagedHull, paint, library, shading, error)) {
            std::printf("scene: %s\n", error.c_str());
            return 1;
        }
        scene.appendOcean(water, static_cast<std::uint32_t>(library.find("sea_water")));

        const sim::Mat4 mvp =
            sim::perspective(45.0 * sim::kDegToRad, double(kWidth) / kHeight, 1.0, 900.0) *
            sim::lookAt(eye, at, {0, 0, 1});
        float matrix[16];
        mvp.toFloats(matrix);
        view.eye[0] = float(eye.x);
        view.eye[1] = float(eye.y);
        view.eye[2] = float(eye.z);
        const float clear[4] = {0.55f, 0.68f, 0.82f, 1.0f};

        core::Image image;
        if (!renderer.render(matrix, view, scene, library, clear, image)) {
            std::printf("render failed on frame %d\n", frame);
            return 1;
        }
        char path[512];
        std::snprintf(path, sizeof(path), "%s/ram_%03d.png", options.out.c_str(), frame);
        if (!core::writePng(path, image)) {
            std::printf("could not write %s\n", path);
            return 1;
        }
        if (frame == 0)
            std::printf("       %zu vertices, %zu triangles, gpu %.3f ms -> %s\n",
                        renderer.lastFrame().vertices, renderer.lastFrame().triangles,
                        renderer.lastFrame().gpuSeconds * 1e3, path);
        ++frame;
    }
    if (drawing) std::printf("       %d frame(s) written to %s\n", frame, options.out.c_str());

    const sim::Diagnostics after = ferry.diagnostics(sea);
    std::printf("\noutcome: %s -- %.0f t of water, heel %.1f deg, GM %.2f m\n",
                (!after.afloat || after.gmTransverse < 0) ? "LOST" : "SURVIVED",
                after.floodwaterMass / 1000.0, after.heelDeg, after.gmTransverse);
    std::printf("ok\n");
    return 0;
}
