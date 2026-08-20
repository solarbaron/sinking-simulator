// SPDX-License-Identifier: MIT
//
// **The Phase 4 milestone, as one run.**
//
// `docs/06-roadmap.md`: *"an engine room fire that heats a bulkhead until it fails
// under the head of water behind it, and the flooding spreads. Three subsystems,
// none of which know about each other, producing one consequence."*
//
// Every piece existed and nothing joined them. `fire.{hpp,cpp}` computes an upper
// layer temperature and had a boundary loss term with no wall behind it;
// `thermal.{hpp,cpp}` conducts through steel and reduces its strength;
// `breach.{hpp,cpp}` turns failed panels into `Opening`s. Not one file mentioned
// another. This tool is the four links joined and driven on the reference ferry:
//
//   1. **fire -> steel.** `fire::wallExchange` splits the bulkhead's boundary faces
//      at the smoke layer interface and hands `thermal::Problem` two films -- hot
//      gas over the upper part, cool air over the lower -- with the radiative half
//      carried exactly by the factored film coefficient. The steel's own surface
//      temperature goes back the other way, into `GasCompartment::wallTemperature`,
//      so the loss the fire was already computing is now the heat the steel is
//      actually taking.
//   2. **steel -> load.** Each bulkhead stiffener with its attached plating is a
//      `thermal::HeatedMember`: restrained thermal expansion as an axial
//      compression, the head of water as a lateral moment integrated off the ship's
//      *own* free surface, and the two joined by the exact beam-column magnifier.
//   3. **failure -> hole.** The panel at the hinge of every failed member goes to
//      `breachesFromFailedPanels`, unchanged, and the openings it returns into
//      `Ship::openings`.
//   4. **hole -> flooding.** `Ship::step` from there, with nothing added.
//
// --- What it is for beyond running the chain -------------------------------------
//
// **The milestone's own sentence is an acceptance test and this tool is the thing
// that runs it.** "Fails *under the head of water behind it*" is a claim that
// neither cause is sufficient alone, so the run is done three times -- fire with a
// dry hold behind the bulkhead, water with no fire, and both -- and `ok` is refused
// unless the two controls survive and the pair does not.
//
// It also measures the *width* of that window, which is the part that could not be
// guessed. The restraint on a heated member is set by the stiffness of the structure
// at its ends -- a bulkhead deck and a tank top -- and that is in none of the three
// subsystems. So `--restraint` is swept, and what comes out is the band of it over
// which the milestone's sentence is true at all.
//
//   ./bulkhead_probe [--duration=S] [--power=W] [--steady=S] [--decay=S]
//                    [--fill=FRACTION] [--restraint=R] [--coupling=S] [--ship-step=S]
//                    [--element=M] [--water-film=W] [--air-film=W] [--cable-transit]
//                    [--sweep] [--quiet]
#include "engine/sim/breach.hpp"
#include "engine/sim/buckling.hpp"
#include "engine/sim/fire.hpp"
#include "engine/sim/scantlings.hpp"
#include "engine/sim/solid_shell.hpp"
#include "engine/sim/thermal.hpp"
#include "game/prototype/ferry.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <string>
#include <vector>

namespace {

using namespace sim;

// The bulkhead this is about: the ferry's watertight bulkhead at x = -8, between
// each engine room and the aft hold behind it. 9.5 mm plating on 180 x 10 flat bars
// at 700 mm, spanning the tank top at 1.8 m to the bulkhead deck at 7.0 m.
constexpr double kStation = -8.0;
constexpr double kTankTop = 1.80;
constexpr double kBulkheadDeck = 7.00;
constexpr double kSpacing = 0.70;
// The meshed band: whole panel rows from 1.4 m up, and whole panel columns out to
// |y| = 7.7 m, which is inside the compartments' own 8 m boundary on both sides.
constexpr double kMeshZLo = 1.40, kMeshZHi = 7.00;
constexpr double kMeshYHalf = 7.70;

struct Options {
    double duration = 3600.0;    // s of run
    double power = 4.0e6;        // W at the plateau
    // **The fire has to end, or there is no control to run.** Steel asymptotes to
    // the gas temperature it stands in, and the restrained-buckling limit is a fixed
    // temperature -- so a fire that burns for ever fells the bulkhead on its own,
    // always, and the only question is when. Measured on the first version of this
    // tool, which held 4 MW indefinitely: the fire alone took the bulkhead at 2760 s
    // and the milestone's sentence was true only over a window in *time* that the
    // tool itself had chosen. A design fire that grows, burns and decays gives the
    // steel a **peak**, and "the fire alone does not fell it" becomes a statement
    // about the whole fire rather than about how long anyone watched.
    double steady = 1800.0;      // s at the plateau
    double decay = 900.0;        // s of linear decay to nothing
    double fill = 0.45;          // fraction of each aft hold's floodable volume
    // The measured centre of the window the tool itself reports: the milestone's own
    // sentence is true for 0.2194 <= r < 0.2505 on this ship under this fire, and this
    // is the geometric middle of it. **It is an input the model cannot derive** -- see
    // `thermal::HeatedMember::restraint` -- so the window is reported on every run and
    // the default is not allowed to be the only evidence.
    double restraint = 0.234;
    double coupling = 5.0;       // s between gas/steel/structure exchanges
    double shipStep = 0.05;      // s
    double element = 0.35;       // m, nominal conduction element size
    double waterFilm = 500.0;    // W/(m^2 K), still water against a steel plate
    double airFilm = 8.0;        // W/(m^2 K), still air in an unheated space
    // The ferry carries an unsealed 0.04 m^2 cable transit through this very
    // bulkhead, authored long before any of this. Left open it is a second, entirely
    // independent leak path, and over a long enough run it floods the machinery
    // spaces, downfloods the vehicle deck and lolls her -- which fells the bulkhead
    // in the *water-only control*, at 3645 s, with nothing burning anywhere. That is
    // a real property of this ship and it is reported rather than hidden, but it is
    // a different experiment, so the milestone run seals it.
    bool   cableTransit = false;
    bool   sweep = false;
    bool   quiet = false;
};

bool parse(int argc, char** argv, Options& o) {
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        const auto value = [&](const char* key) -> const char* {
            const std::string prefix = std::string("--") + key + "=";
            return a.rfind(prefix, 0) == 0 ? a.c_str() + prefix.size() : nullptr;
        };
        if (const char* v = value("duration")) o.duration = std::atof(v);
        else if (const char* v = value("power")) o.power = std::atof(v);
        else if (const char* v = value("fill")) o.fill = std::atof(v);
        else if (const char* v = value("restraint")) o.restraint = std::atof(v);
        else if (const char* v = value("coupling")) o.coupling = std::atof(v);
        else if (const char* v = value("ship-step")) o.shipStep = std::atof(v);
        else if (const char* v = value("element")) o.element = std::atof(v);
        else if (const char* v = value("water-film")) o.waterFilm = std::atof(v);
        else if (const char* v = value("air-film")) o.airFilm = std::atof(v);
        else if (const char* v = value("steady")) o.steady = std::atof(v);
        else if (const char* v = value("decay")) o.decay = std::atof(v);
        else if (a == "--cable-transit") o.cableTransit = true;
        else if (a == "--sweep") o.sweep = true;
        else if (a == "--quiet") o.quiet = true;
        else {
            std::printf("unknown option %s\n", a.c_str());
            return false;
        }
    }
    return true;
}

// --- The bulkhead, as three views of the same steel ------------------------------
//
// A column of plating panels, the stiffener down its edge, and the strip of
// conduction elements over it. One index space, built once, so nothing downstream
// has to search for the panel a member owns.
struct Column {
    double y = 0;                  // m, the centre of the strip
    std::vector<int> panel;        // ascending in z, indices into StructuralMesh::panels
    std::vector<double> panelZ;    // m, each panel's centroid height
    std::vector<double> panelArea; // m^2
};

// The conduction mesh, in the body frame, and the map from a (y, z) point to the
// element over it.
struct Slab {
    solidshell::HexMesh mesh;
    int nx = 0, ny = 0;            // elements across y and up z
    double y0 = 0, z0 = 0, dy = 0, dz = 0;

    int elementAt(double y, double z) const {
        const int i = std::clamp(static_cast<int>((y - y0) / dy), 0, nx - 1);
        const int j = std::clamp(static_cast<int>((z - z0) / dz), 0, ny - 1);
        return i * ny + j;   // makePlateMesh runs the second index fastest
    }
};

// `makePlateMesh` builds in its own frame: the plate in x-y with the thickness
// along z and the mid-surface at zero. The bulkhead is at a station, so the two are
// related by the cyclic permutation (x, y, z) -> (z, x, y), which is a rotation and
// therefore leaves every element Jacobian positive.
Slab buildSlab(double element) {
    Slab s;
    const double ly = kMeshYHalf * 2.0, lz = kMeshZHi - kMeshZLo;
    s.nx = std::max(2, static_cast<int>(std::lround(ly / element)));
    s.ny = std::max(2, static_cast<int>(std::lround(lz / element)));
    s.mesh = solidshell::makePlateMesh(ly, lz, 0.0095, s.nx, s.ny, 1);
    s.y0 = -kMeshYHalf;
    s.z0 = kMeshZLo;
    s.dy = ly / s.nx;
    s.dz = lz / s.ny;
    for (std::size_t n = 0; n < s.mesh.nodeCount(); ++n) {
        double* p = &s.mesh.position[n * 3];
        const double a = p[0], b = p[1], c = p[2];
        p[0] = kStation + c;   // through the thickness, mid-surface on the station
        p[1] = s.y0 + a;
        p[2] = s.z0 + b;
    }
    return s;
}

// --- The head of water, and the moment it puts in a member -----------------------

// Differential pressure across the bulkhead at a body-frame point, Pa, positive
// pushing from the flooded side into the compartment on fire.
//
// The same comparison `Ship::sideStateAt` makes and for the reason
// `fire::Scupper` records: a *position* against the compartment's own
// `surfaceWorldZ`, not a height against `surfaceOffset`, so it stays exact when the
// ship takes up an angle. The gas buoyancy head `sideStateAt` also carries is the
// one term not reproduced -- it is tens of pascals against tens of kilopascals here
// and it is not public.
double waterHead(const Ship& ship, int wet, const Vec3& bodyPoint) {
    const Mat3 R = ship.state.orientation.toMat3();
    const Vec3 world = ship.state.position + R * bodyPoint;
    const Compartment& w = ship.compartments[static_cast<std::size_t>(wet)];
    if (w.waterVolume <= 1e-9) return 0.0;
    return std::max(0.0, w.surfaceWorldZ - world.z);
}

double differentialPressure(const Ship& ship, int wet, int dry, const Vec3& bodyPoint) {
    const Compartment& w = ship.compartments[static_cast<std::size_t>(wet)];
    const Compartment& d = ship.compartments[static_cast<std::size_t>(dry)];
    return w.airPressure + ship.seaDensity * kGravity * waterHead(ship, wet, bodyPoint) -
           d.airPressure;
}

// Largest bending moment in a pin-ended member spanning `zLo` to `zHi` under a
// lateral load `w(z)`, by direct integration of the beam equations.
//
// It is not `q L^2 / 8`: a bulkhead's load is hydrostatic over the part of its span
// that is under water and zero above it, so both the magnitude and the *place* of
// the peak move as the compartment behind fills. The place is what decides which
// panel the hinge is in, so it is returned rather than assumed to be mid-span.
struct BeamMoment {
    double moment = 0;   // N m
    double at = 0;       // m, the height of the peak
    double load = 0;     // N, the total lateral load on the member
};

template <typename Load>
BeamMoment peakMoment(double zLo, double zHi, const Load& w, int samples = 512) {
    BeamMoment out;
    const double span = zHi - zLo;
    if (!(span > 0)) return out;
    const double h = span / samples;
    // Reaction at the lower support of a simply supported span: the load's moment
    // about the upper one, divided by the span.
    double total = 0, aboutTop = 0;
    std::vector<double> q(static_cast<std::size_t>(samples));
    for (int i = 0; i < samples; ++i) {
        const double z = zLo + (i + 0.5) * h;
        q[static_cast<std::size_t>(i)] = w(z) * h;
        total += q[static_cast<std::size_t>(i)];
        aboutTop += q[static_cast<std::size_t>(i)] * (zHi - z);
    }
    out.load = total;
    const double lower = span > 0 ? aboutTop / span : 0.0;
    // M(z) = R_lower (z - zLo) - sum over the load below z of q (z - z_q).
    double carried = 0, firstMoment = 0;
    for (int i = 0; i < samples; ++i) {
        const double z = zLo + (i + 1) * h;
        const double zq = zLo + (i + 0.5) * h;
        carried += q[static_cast<std::size_t>(i)];
        firstMoment += q[static_cast<std::size_t>(i)] * zq;
        const double m = lower * (z - zLo) - (carried * z - firstMoment);
        if (std::abs(m) > std::abs(out.moment)) {
            out.moment = m;
            out.at = z;
        }
    }
    return out;
}

// --- One run of the chain ---------------------------------------------------------

struct Outcome {
    double failureTime = -1;        // s, first member to reach utilisation 1
    double failureTemperature = 0;  // K, that member's equivalent uniform temperature
    double failureUtilisation = 0;
    double failureAdditive = 0;     // the same sum with the magnifier forced to 1
    double failureMagnifier = 0;
    double failureMoment = 0;       // N m
    double failureAxial = 0;        // Pa
    double failureY = 0, failureZ = 0;
    thermal::MemberLimit failureLimit = thermal::MemberLimit::None;

    int    failedMembers = 0;
    int    failedPanels = 0;
    double openingArea = 0;         // m^2 reaching a compartment
    int    openings = 0;

    double peakSteel = 0;           // K, hottest element on the bulkhead
    // The hottest *member*: its equivalent uniform temperature, which is what the
    // failure check is a function of and which is well below the hottest element,
    // because a member spans from the deckhead the smoke is against to the foot the
    // water behind is cooling.
    double peakMember = 0;          // K
    double peakGas = 0;             // K, hottest upper layer
    double peakUtilisation = 0;
    double peakHead = 0;            // m of water on the bulkhead
    double peakMoment = 0;          // N m in any member
    double plateBucklingTime = -1;  // s the plating first went, which is not a hole

    double filmHeat = 0;            // J into the steel
    double enthalpyGain = 0;        // J the steel actually holds
    double gasWallLoss = 0;         // J the fire booked as boundary loss
    double worstLinearisation = 0;  // W, |heat - exactHeat| at any exchange
    double exchange = 0;            // W, the largest exchange either half ever carried

    // Where the water went, at the end.
    double waterBehind = 0, waterFireSide = 0, waterElsewhere = 0;   // m^3
    int    compartmentsWet = 0;
    double heel = 0, draft = 0, gm = 0;
    bool   afloat = true;

    int    couplingSteps = 0, thermalFactorisations = 0;
    double wallSeconds = 0;

    // **The three solves are stepped by one `o.coupling` each and printed under one
    // `t`, and only one of them could fall behind without saying so.** `gas.step()`
    // may return having advanced less than it was asked for -- its substep budget
    // bounds trials, and a rejected trial spends a slot without committing time --
    // so the gas temperature, the interface height and the steel temperature on a
    // single row of the trace above would be readings at three different model
    // times, and the row would look exactly like every other one. Counted here and
    // required to be zero below, which is the discipline `flip_probe` already keeps
    // over `flip::StepResult::incomplete`.
    int    gasShortSteps = 0;
    // The gas pressure solve failed to bracket a root. `fire.hpp` publishes this
    // "rather than swallowed so that 'should never' can be asserted", and until now
    // nothing outside `tests/test_fire.cpp` asserted it.
    int    gasPressureCapped = 0;
    // The worst disagreement between the gas's own clock and the `t` this loop
    // prints. Zero to rounding on a healthy run; it is the direct measurement of
    // the failure the two counters above only imply.
    double gasClockGap = 0;          // s

    // Every (equivalent uniform temperature, lateral moment) this run ever put into
    // a member, so that the restraint window can be found from the run itself rather
    // than from a summary of it. The two peaks do not happen at the same moment --
    // the head is greatest before the compartment behind starts emptying into the
    // one on fire, and the steel is hottest well after -- so a window computed from
    // `peakMember` and `peakMoment` together would be a state the run never reached.
    std::vector<std::pair<double, double>> history;
};

// The restraint at which a run's own history first reaches utilisation 1 somewhere:
// bisected on the run's recorded states, so it is the answer for the states that
// actually occurred and not for an idealisation of them. Zero if no restraint in
// (0, 4] would have done it.
double restraintWindowBound(const std::vector<std::pair<double, double>>& history,
                            const StructuralMaterial& steel, double modulus,
                            double eulerStress) {
    const auto worst = [&](double r) {
        double u = 0;
        thermal::HeatedMember m;
        m.modulus = modulus;
        m.eulerStress = eulerStress;
        m.restraint = r;
        for (const std::pair<double, double>& s : history) {
            m.lateralMoment = s.second;
            u = std::max(u, thermal::memberState(m, steel, s.first).utilisation);
        }
        return u;
    };
    if (worst(4.0) < 1.0) return 0.0;
    double lo = 0.0, hi = 4.0;
    for (int i = 0; i < 60; ++i) {
        const double m = 0.5 * (lo + hi);
        if (worst(m) >= 1.0) hi = m;
        else lo = m;
    }
    return 0.5 * (lo + hi);
}

struct Chain {
    const Options* o = nullptr;
    Scantlings scantlings;
    StructuralMesh structure;
    Slab slab;
    std::vector<Column> column;          // ascending in y
    std::vector<thermal::BoundaryFace> fireFace, wetFace;
    std::vector<int> fireFaceColumn;   // which column each fire-side face stands over
    StiffenedSection section;
    double eulerStress = 0;
    double plateElastic = 0;

    // The bulkhead's stiffener with its attached plating, as the failure check wants
    // it. One place, so the tool's own reporting and its member loop cannot disagree
    // about which member they are talking about.
    thermal::HeatedMember member(double restraint) const {
        thermal::HeatedMember m;
        m.modulus = section.modulusStiffener;
        m.eulerStress = eulerStress;
        m.restraint = restraint;
        return m;
    }
};

// Every exterior face of the slab that looks along +x (into the engine rooms) or
// -x (into the aft holds), with the column it stands over. The four edge faces
// carry no film: they are the cut this idealisation makes and heating them would be
// inventing a boundary.
void classifyFaces(Chain& c) {
    const std::vector<thermal::BoundaryFace> faces = thermal::boundaryFaces(c.slab.mesh);
    for (const thermal::BoundaryFace& f : faces) {
        if (!(f.area > 0)) continue;
        int col = -1;
        for (std::size_t k = 0; k < c.column.size(); ++k)
            if (std::abs(f.centroid.y - c.column[k].y) <= 0.5 * kSpacing + 1e-9) {
                col = static_cast<int>(k);
                break;
            }
        if (f.normal.x > 0.5) {
            c.fireFace.push_back(f);
            c.fireFaceColumn.push_back(col);
        } else if (f.normal.x < -0.5) {
            c.wetFace.push_back(f);
        }
    }
}

Chain buildChain(const Ship& ferry, const Options& o) {
    Chain c;
    c.o = &o;
    c.scantlings = ferryScantlings();
    c.structure = makeStructuralMesh(ferry.hull, c.scantlings);
    c.slab = buildSlab(o.element);

    // Panel columns of the bulkhead, from the structural mesh's own panels.
    std::vector<double> centres;
    for (int k = -11; k <= 11; ++k) {
        const double y = (k >= 0 ? 0.5 : -0.5) * kSpacing + k * kSpacing;
        if (std::abs(y) < kMeshYHalf) centres.push_back(y);
    }
    std::sort(centres.begin(), centres.end());
    for (double y : centres) {
        Column col;
        col.y = y;
        for (std::size_t i = 0; i < c.structure.panels.size(); ++i) {
            const PlatePanel& p = c.structure.panels[i];
            if (p.role != PanelRole::Bulkhead) continue;
            const Vec3 g = p.centroid();
            if (std::abs(g.x - kStation) > 0.1) continue;
            if (std::abs(g.y - y) > 1e-6) continue;
            if (g.z < kMeshZLo || g.z > kMeshZHi) continue;
            col.panel.push_back(static_cast<int>(i));
            col.panelZ.push_back(g.z);
            col.panelArea.push_back(p.area());
        }
        if (!col.panel.empty()) c.column.push_back(col);
    }
    classifyFaces(c);

    const StructuralMaterial steel = c.scantlings.materials[0];
    c.section = stiffenedSection(flatBar(0.180, 0.010), 0.0095, kSpacing);
    c.eulerStress =
        columnBuckling(c.section, kBulkheadDeck - kTankTop, 1.0, steel).elasticStress;
    c.plateElastic = plateBuckling(0.0095, kSpacing, kSpacing, 1.0, steel).elasticStress;
    return c;
}

// Run the chain once. `fire` lights the design fire; `water` floods the aft holds.
// `applyDamage` false evaluates every member exactly as usual and then does *not*
// open the hole. It exists because the failure relieves what caused it -- the
// compartment behind empties into the one on fire and the head falls -- so a run
// that fails cannot say how far past its limit it *would* have gone, and the
// restraint window's lower bound computed from a relieved run comes back equal to
// whatever restraint that run was given. Self-referential, and it read as a
// suspiciously exact answer twice before it was noticed.
Outcome run(const Chain& chain, const Options& o, bool fire, bool water, bool verbose,
            bool applyDamage = true) {
    Outcome out;
    const double started = static_cast<double>(std::clock()) / CLOCKS_PER_SEC;
    const StructuralMaterial steel = chain.scantlings.materials[0];
    const Sea sea{0.0};

    Ship ship = game::buildFerry();
    // This casualty is a fire, not a collision: the ferry's authored side breach is
    // shut, so the only way water reaches the machinery space is the bulkhead.
    for (Opening& op : ship.openings) {
        if (op.name == "breach_er_s") op.open = false;
        if (op.name == "cable_transit" && !o.cableTransit) op.open = false;
    }

    const int aftP = ship.findCompartment("aft_hold_p");
    const int aftS = ship.findCompartment("aft_hold_s");
    const int erP = ship.findCompartment("engine_room_p");
    const int erS = ship.findCompartment("engine_room_s");
    // **The water goes in before `initialise`, and that is not a tidiness point.**
    // `initialise` sets each compartment's air mass from the gas volume it actually
    // has and then solves the flooded equilibrium draft. Filling a hold afterwards
    // leaves the air mass that was there when the hold was empty, so the trapped gas
    // reads 184 kPa -- 8.3 m of head that is not water -- and the bulkhead fails in
    // the first step of every run, controls included. Measured, on the first version
    // of this tool.
    ship.initialise(sea);   // caches the gross volumes the fill is a fraction of
    if (water) {
        // A grounding aft, both holds together, which is what keeps her upright: the
        // asymmetric case is a stability problem and this is a structural one.
        // **Checked, because `findCompartment` returns -1 and this is a store.**
        // `static_cast<std::size_t>(-1)` is `SIZE_MAX`, so a renamed or missing
        // compartment writes to a wild address rather than reading a wrong one.
        // The ferry has been short of an authored compartment before -- her mid
        // wing tanks were never written -- and `water_probe` guards the same call.
        for (int c : {aftP, aftS}) {
            if (c < 0) {
                std::printf("no aft hold to flood: expected aft_hold_p and aft_hold_s\n");
                return out;
            }
            ship.compartments[static_cast<std::size_t>(c)].waterVolume =
                o.fill * ship.compartments[static_cast<std::size_t>(c)].floodableVolume();
        }
        ship.initialise(sea);
    }

    fire::Model gas;
    gas.attach(ship, {erP, erS});
    const int burning = gas.gasIndexOf(erP);
    if (fire) {
        fire::DesignFire d;
        d.name = "machinery";
        d.compartment = burning;
        d.baseZ = 2.5;
        d.diameter = 2.5;
        d.growthCoefficient = fire::kGrowthFast;
        d.peakHeatRelease = o.power;
        d.steadyDuration = o.steady;
        d.decayDuration = o.decay;
        gas.fires.push_back(d);
    }

    // --- The conduction problem, and why every film's membership is fixed --------
    //
    // Films are banded by height at the element size: one per element row per half
    // on the fire side, and one per row on the wet side. **Nothing about that
    // membership depends on the state**, so the whole run is one `prepare` and every
    // step is `setFilm`, which is the difference between an energy account that
    // spans the run and one that restarts whenever the smoke layer moves. The first
    // version of this tool re-prepared on every interface crossing and reported
    // 0.0020 GJ through the films against 0.2015 GJ held by the steel -- a hundredfold
    // "residual" that was nothing but an account being zeroed twelve times.
    //
    // Banding is also what makes the radiative coefficient honest: `fire.hpp`'s note
    // at `WallExchange` measures 8% of the exchange lost to one coefficient over a
    // whole layer, on this very bulkhead, because the foot of it is held near
    // ambient by the water behind while its head is at 500 K.
    const auto half = [&](int c) { return chain.column[static_cast<std::size_t>(c)].y > 0 ? 0 : 1; };
    std::vector<thermal::BoundaryFace> fireFaceHalf[2];
    std::vector<double> fireTempHalf[2];
    for (std::size_t i = 0; i < chain.fireFace.size(); ++i) {
        const int col = chain.fireFaceColumn[i];
        if (col < 0) continue;
        fireFaceHalf[half(col)].push_back(chain.fireFace[i]);
        fireTempHalf[half(col)].push_back(kTAmbient);
    }

    thermal::Problem problem;
    problem.mesh = &chain.slab.mesh;
    problem.material = steel;
    problem.temperatureDependent = true;

    // A cold-start exchange, only to learn how many bands each half has and in what
    // order -- `wallExchange` bands by the face's own centroid, so the answer does
    // not depend on the state it is asked in.
    std::size_t filmBase[2] = {0, 0};
    std::size_t bandCount[2] = {0, 0};
    for (int h = 0; h < 2; ++h) {
        const int gi = gas.gasIndexOf(h == 0 ? erP : erS);
        const fire::WallExchange x = fire::wallExchange(gas.gas[static_cast<std::size_t>(gi)],
                                                       fireFaceHalf[h], fireTempHalf[h], {},
                                                       chain.slab.dz);
        filmBase[h] = problem.film.size();
        bandCount[h] = x.film.size();
        for (const thermal::Film& f : x.film) problem.film.push_back(f);
    }
    // The wet side: the same banding, one film per element row, so that the film
    // holding a row can be switched between water and air as the level behind the
    // bulkhead rises without any face changing hands.
    const std::size_t wetBase = problem.film.size();
    std::vector<std::vector<thermal::BoundaryFace>> wetBand(static_cast<std::size_t>(chain.slab.ny));
    for (const thermal::BoundaryFace& f : chain.wetFace) {
        const auto b = static_cast<std::size_t>((f.centroid.z - kMeshZLo) / chain.slab.dz);
        wetBand[std::min(b, wetBand.size() - 1)].push_back(f);
    }
    for (const std::vector<thermal::BoundaryFace>& b : wetBand) {
        thermal::Film f;
        f.face = b;
        f.ambient = kTAmbient;
        problem.film.push_back(f);
    }

    thermal::Solver solver;
    std::string why;
    if (!solver.prepare(problem, kTAmbient, &why)) {
        std::printf("thermal solve refused: %s\n", why.c_str());
        return out;
    }

    std::vector<double> elementT;
    std::vector<char> failed(chain.column.size(), 0);
    std::vector<int> failedPanel;
    const double startEnthalpy = solver.account().enthalpy;

    // **Four of these columns are not the value at `t`, and the header now says so.**
    // `peakSteel`, `peakHead` and `peakUtilisation` are running maxima over the whole
    // run -- monotone, never the instantaneous reading -- and `failureAxial` is
    // latched at the moment a member goes, so it prints 0.00 until then. Only
    // `gas`, `iface` and `M` are recomputed each report. Read as a time series the
    // peak columns erase the very thing §the coupling says matters: the head is
    // greatest *before* the compartment behind starts emptying, and the steel is
    // hottest well *after*, and a monotone trace can show neither.
    if (verbose)
        std::printf("\n%8s %9s %9s %9s %9s %9s %9s %9s %7s\n", "t (s)", "gas (K)",
                    "pkSteel(K)", "iface (m)", "pkHead(m)", "M (kN m)", "sigN@fail", "pkUtil",
                    "failed");

    const int steps = static_cast<int>(std::lround(o.duration / o.coupling));
    for (int step = 0; step < steps; ++step) {
        // --- 1. the gas -------------------------------------------------------
        // Read rather than discarded, for the reason `Outcome::gasShortSteps` gives:
        // this loop's time axis is arithmetic and the gas's is not.
        const fire::StepResult gasStep = gas.step(o.coupling, ship, sea);
        if (gasStep.incomplete) ++out.gasShortSteps;
        if (gasStep.pressureSolveCapped) ++out.gasPressureCapped;
        gas.applyTo(ship);
        for (const fire::GasCompartment& g : gas.gas)
            out.peakGas = std::max(out.peakGas, g.upper.temperature());

        // --- 2. the boundary, both ways ---------------------------------------
        const double interfaceZ = gas.gas[static_cast<std::size_t>(burning)].interfaceZ();
        const std::vector<double>& nodal = solver.temperature();
        for (int h = 0; h < 2; ++h) {
            // Surface temperatures, by averaging the solve's own nodal field over
            // each face's four corners.
            for (std::size_t i = 0; i < fireFaceHalf[h].size(); ++i) {
                const thermal::BoundaryFace& f = fireFaceHalf[h][i];
                double t = 0;
                for (int k = 0; k < 4; ++k) t += nodal[f.node[k]];
                fireTempHalf[h][i] = 0.25 * t;
            }
            const int gi = gas.gasIndexOf(h == 0 ? erP : erS);
            const fire::WallExchange x = fire::wallExchange(gas.gas[static_cast<std::size_t>(gi)],
                                                            fireFaceHalf[h], fireTempHalf[h], {},
                                                            chain.slab.dz);
            for (std::size_t b = 0; b < x.film.size() && b < bandCount[h]; ++b)
                solver.setFilm(filmBase[h] + b, 0.0, x.film[b].coefficient, x.film[b].ambient);
            // The other direction: the steel's own surface is now what the fire is
            // losing heat to, and the gas-side film is what it is losing it through.
            fire::GasCompartment& g = gas.gas[static_cast<std::size_t>(gi)];
            g.wallTemperature = x.wallTemperature;
            g.wallConductance = x.wallConductance;
            out.worstLinearisation =
                std::max(out.worstLinearisation, std::abs(x.linearisationError));
            out.exchange = std::max(out.exchange, std::abs(x.exactHeat));
        }
        // The wet side, row by row. Water below the level behind the bulkhead, still
        // air above it -- and the water is the reason the foot of the bulkhead, which
        // is where the head is greatest, is also the part that never gets hot.
        const double wetSurface =
            water && aftP >= 0
                ? ship.compartments[static_cast<std::size_t>(aftP)].surfaceWorldZ -
                        ship.state.position.z
                  : kMeshZLo - 1.0;
        for (std::size_t b = 0; b < wetBand.size(); ++b) {
            const double z = kMeshZLo + (static_cast<double>(b) + 0.5) * chain.slab.dz;
            solver.setFilm(wetBase + b, 0.0, z <= wetSurface ? o.waterFilm : o.airFilm, kTAmbient);
        }

        // --- 3. the steel -----------------------------------------------------
        if (!solver.step(o.coupling, &why)) {
            std::printf("thermal step refused: %s\n", why.c_str());
            return out;
        }
        // Checked, like `solver.step` six lines up. On refusal `elementTemperatures`
        // *clears* its output, so the loop below runs zero times and `peakSteel`
        // silently keeps the value it already had -- its initialiser on the first
        // step, a stale maximum on any later one. `peakSteel` is a published
        // milestone figure this file's own `require` only bounds from below, so a
        // partial failure would have published a wrong peak rather than no peak.
        if (!thermal::elementTemperatures(chain.slab.mesh, solver.temperature(), elementT)) {
            std::printf("element temperatures refused: the field is not one value per node\n");
            return out;
        }
        for (double t : elementT) out.peakSteel = std::max(out.peakSteel, t);

        // --- 4. every member of the bulkhead ----------------------------------
        const double now = (step + 1) * o.coupling;
        // The printed axis against the model that is meant to be on it.
        out.gasClockGap = std::max(out.gasClockGap, std::abs(gasStep.time - now));
        for (std::size_t k = 0; k < chain.column.size(); ++k) {
            const Column& col = chain.column[k];
            const bool port = col.y > 0;
            const int wet = port ? aftP : aftS, dry = port ? erP : erS;

            // The member's equivalent uniform temperature: the mean of its own
            // elongations, inverted. Not the mean of its temperatures -- the
            // elongation curve is quadratic and a member is restrained by the length
            // it wanted to grow, which is an integral.
            double elong = 0;
            int taken = 0;
            for (double z = kTankTop + 0.5 * chain.slab.dz; z < kBulkheadDeck;
                 z += chain.slab.dz) {
                const double t = elementT[static_cast<std::size_t>(
                    chain.slab.elementAt(col.y, std::min(z, kMeshZHi - 1e-6)))];
                elong += thermal::carbonSteelElongation(t);
                ++taken;
            }
            const double equivalent =
                taken > 0 ? thermal::temperatureForElongation(elong / taken) : kTAmbient;

            const BeamMoment bending = peakMoment(kTankTop, kBulkheadDeck, [&](double z) {
                return kSpacing *
                       std::max(0.0, differentialPressure(ship, wet, dry,
                                                          Vec3{kStation, col.y, z}));
            });

            thermal::HeatedMember member = chain.member(o.restraint);
            member.lateralMoment = bending.moment;
            const thermal::MemberState state =
                thermal::memberState(member, steel, equivalent);

            out.peakUtilisation = std::max(out.peakUtilisation, state.utilisation);
            out.peakMoment = std::max(out.peakMoment, std::abs(bending.moment));
            out.peakMember = std::max(out.peakMember, equivalent);
            out.history.emplace_back(equivalent, bending.moment);
            if (failed[k] || state.utilisation < 1.0) continue;

            failed[k] = 1;
            ++out.failedMembers;
            if (out.failureTime < 0) {
                out.failureTime = now;
                out.failureTemperature = equivalent;
                out.failureUtilisation = state.utilisation;
                out.failureAdditive = state.additiveUtilisation;
                out.failureMagnifier = state.magnifier;
                out.failureMoment = bending.moment;
                out.failureAxial = state.axialStress;
                out.failureY = col.y;
                out.failureZ = bending.at;
                out.failureLimit = state.limit;
            }
            // The panel the hinge is in. A plastic hinge is a fold, not a collapse
            // of the whole strip, so one panel goes and `breach.hpp` merges whatever
            // its neighbours do -- which is the mechanism that makes several failing
            // members one hole rather than several.
            int best = -1;
            double bestGap = 1e30;
            for (std::size_t p = 0; p < col.panel.size(); ++p)
                if (std::abs(col.panelZ[p] - bending.at) < bestGap) {
                    bestGap = std::abs(col.panelZ[p] - bending.at);
                    best = col.panel[p];
                }
            if (best >= 0) failedPanel.push_back(best);
        }

        // The plating between the stiffeners, under the horizontal compression the
        // hot band takes from the cold band beside it -- `twoStripStress` at the hot
        // area fraction the fire's own interface height gives. Recorded, and
        // deliberately **not** a hole: a buckled plate panel goes out of plane and
        // sheds its in-plane load, and it is still watertight.
        //
        // Taken at the *hottest element on the bulkhead* rather than at a mean, so
        // it is the earliest moment any panel could have gone rather than a
        // statement about the panel that did. A diagnostic wants the bound.
        if (out.plateBucklingTime < 0) {
            const double hotFraction =
                std::clamp((kBulkheadDeck - interfaceZ) / (kBulkheadDeck - kTankTop), 0.0, 1.0);
            const double hottest = out.peakSteel;
            if (std::abs(thermal::twoStripStress(steel, hottest, hotFraction)) >=
                johnsonOstenfeld(thermal::carbonSteelModulusFactor(hottest) * chain.plateElastic,
                                 thermal::carbonSteelYieldFactor(hottest) * steel.yieldStrength))
                out.plateBucklingTime = now;
        }

        // --- 5. the hole, and the flooding ------------------------------------
        if (!failedPanel.empty() && applyDamage) {
            const BreachSet set = breachesFromFailedPanels(ship, chain.structure, failedPanel);
            // Rebuild the ship's own opening list from scratch each time rather than
            // appending: `breachesFromFailedPanels` merges edge-adjacent panels, so
            // the *set* of openings changes when a member joins one, and appending
            // would count the merged area twice.
            ship.openings.erase(std::remove_if(ship.openings.begin(), ship.openings.end(),
                                               [](const Opening& op) {
                                                   return op.name.rfind("breach_", 0) == 0 &&
                                                          op.name != "breach_er_s";
                                               }),
                                ship.openings.end());
            applyBreaches(ship, set);
            out.failedPanels = static_cast<int>(failedPanel.size());
            out.openings = static_cast<int>(set.breaches.size());
            out.openingArea = set.totalArea();
        }

        const int shipSteps = std::max(1, static_cast<int>(std::lround(o.coupling / o.shipStep)));
        for (int i = 0; i < shipSteps; ++i) ship.step(o.coupling / shipSteps, sea);

        for (std::size_t k = 0; k < chain.column.size(); ++k) {
            const bool port = chain.column[k].y > 0;
            out.peakHead = std::max(out.peakHead,
                                    waterHead(ship, port ? aftP : aftS,
                                              Vec3{kStation, chain.column[k].y, kTankTop}));
        }

        ++out.couplingSteps;
        if (verbose && (step % std::max(1, steps / 20) == 0 || step == steps - 1)) {
            // `size() - 6` wraps to near `SIZE_MAX` on a bulkhead with fewer than
            // six columns -- a narrower one, a different ship, or a wider spacing.
            if (chain.column.size() < 6) continue;
            const Column& probe = chain.column[chain.column.size() - 6];
            const BeamMoment b = peakMoment(kTankTop, kBulkheadDeck, [&](double z) {
                return kSpacing * std::max(0.0, differentialPressure(ship, aftP, erP,
                                                                     Vec3{kStation, probe.y, z}));
            });
            std::printf("%8.0f %9.2f %9.2f %9.3f %9.3f %9.2f %9.2f %9.4f %7d\n", now,
                        gas.gas[static_cast<std::size_t>(burning)].upper.temperature(),
                        out.peakSteel, interfaceZ, out.peakHead, b.moment / 1e3,
                        out.failureAxial / 1e6, out.peakUtilisation, out.failedMembers);
        }
    }

    out.filmHeat = solver.account().filmHeat;
    out.enthalpyGain = solver.account().enthalpy - startEnthalpy;
    out.gasWallLoss = gas.account.wallLoss;
    out.thermalFactorisations = solver.factorisations();

    for (std::size_t i = 0; i < ship.compartments.size(); ++i) {
        const Compartment& c = ship.compartments[i];
        if (c.waterVolume <= 1.0) continue;
        ++out.compartmentsWet;
        const int idx = static_cast<int>(i);
        if (idx == aftP || idx == aftS) out.waterBehind += c.waterVolume;
        else if (idx == erP || idx == erS) out.waterFireSide += c.waterVolume;
        else out.waterElsewhere += c.waterVolume;
    }
    const Diagnostics d = ship.diagnostics(sea);
    out.heel = d.heelDeg;
    out.draft = d.draftMidship;
    out.gm = d.gmTransverse;
    out.afloat = d.afloat;
    out.wallSeconds = static_cast<double>(std::clock()) / CLOCKS_PER_SEC - started;
    return out;
}

const char* limitName(thermal::MemberLimit l) {
    switch (l) {
        case thermal::MemberLimit::Column: return "column";
        case thermal::MemberLimit::Interaction: return "interaction";
        default: return "none";
    }
}

}  // namespace

int main(int argc, char** argv) {
    Options o;
    if (!parse(argc, argv, o)) return 2;

    const Sea sea{0.0};
    Ship reference = game::buildFerry();
    reference.initialise(sea);
    const Chain chain = buildChain(reference, o);

    std::printf("bulkhead: x = %.1f, %zu stiffener strips at %.2f m, %.1f mm plating,"
                " span %.2f m\n",
                kStation, chain.column.size(), kSpacing, 9.5, kBulkheadDeck - kTankTop);
    std::printf("          section A %.4e m2, I %.4e m4, Z %.4e m3; Euler %.2f MPa,"
                " plate %.2f MPa\n",
                chain.section.area, chain.section.secondMoment, chain.section.modulusStiffener,
                chain.eulerStress / 1e6, chain.plateElastic / 1e6);
    std::printf("conduct : %zu elements, %zu nodes, %d x %d over %.2f x %.2f m\n",
                chain.slab.mesh.elementCount(), chain.slab.mesh.nodeCount(), chain.slab.nx,
                chain.slab.ny, 2.0 * kMeshYHalf, kMeshZHi - kMeshZLo);
    std::printf("films   : %zu faces on the fire side, %zu on the wet side\n",
                chain.fireFace.size(), chain.wetFace.size());
    std::printf("run     : %.0f s at a %.1f s coupling step, %.0f MW, aft holds %.0f%% full,"
                " restraint %.3f\n",
                o.duration, o.coupling, o.power / 1e6, 100.0 * o.fill, o.restraint);

    // --- The two controls, and the pair ------------------------------------------
    std::printf("\n=== fire alone (the hold behind the bulkhead is dry) ===\n");
    const Outcome fireOnly = run(chain, o, true, false, !o.quiet);
    std::printf("\n=== the head alone (nothing is burning) ===\n");
    const Outcome waterOnly = run(chain, o, false, true, !o.quiet);
    std::printf("\n=== both ===\n");
    const Outcome both = run(chain, o, true, true, !o.quiet);
    // The same case with the damage not applied, so the restraint window's lower
    // bound is a measurement rather than a restatement of the run's own input.
    const Outcome unrelieved = run(chain, o, true, true, false, /*applyDamage=*/false);

    const auto line = [](const char* label, const Outcome& r) {
        std::printf("%-14s %7s %9.0f %9.1f %9.1f %9.4f %8d %9.2f %10.1f %5d\n", label,
                    r.failureTime >= 0 ? "yes" : "no",
                    r.failureTime >= 0 ? r.failureTime : 0.0, r.peakSteel - thermal::kCelsius,
                    r.peakMember - thermal::kCelsius, r.peakUtilisation, r.failedMembers,
                    r.openingArea, r.waterFireSide, r.compartmentsWet);
    };
    std::printf("\n%-14s %7s %9s %9s %9s %9s %8s %9s %10s %5s\n", "", "failed", "at (s)",
                "steel C", "member C", "worst u", "members", "hole m2", "into ER m3", "wet");
    line("fire alone", fireOnly);
    line("head alone", waterOnly);
    line("both", both);

    if (both.failureTime >= 0)
        std::printf("\nthe bulkhead went at t = %.0f s, y = %+.2f m, hinge at z = %.2f m,"
                    " %s limit\n",
                    both.failureTime, both.failureY, both.failureZ, limitName(both.failureLimit));
    std::printf("  member at %.1f C: axial %.2f MPa, moment %.2f kN m, magnifier %.3f\n",
                both.failureTemperature - thermal::kCelsius, both.failureAxial / 1e6,
                both.failureMoment / 1e3, both.failureMagnifier);
    std::printf("  utilisation %.4f, of which the magnification is %.4f -- a purely additive"
                " check reads %.4f and has not failed\n",
                both.failureUtilisation, both.failureUtilisation - both.failureAdditive,
                both.failureAdditive);
    if (both.plateBucklingTime >= 0)
        std::printf("  the plating between the stiffeners buckled at t = %.0f s, which is not"
                    " a hole: a buckled panel sheds its in-plane load and stays watertight\n",
                    both.plateBucklingTime);
    else
        std::printf("  the plating between the stiffeners never buckled: the two-strip restraint"
                    " a hot band takes from the cold one below it never reached %.1f MPa\n",
                    chain.plateElastic / 1e6);

    std::printf("\nflooding: %.1f m3 behind the bulkhead, %.1f m3 into the machinery spaces,"
                " %.1f m3 beyond, %d compartments wet\n",
                both.waterBehind, both.waterFireSide, both.waterElsewhere, both.compartmentsWet);
    std::printf("          draft %.3f m, heel %.2f deg, GM %.3f m, afloat %d\n", both.draft,
                both.heel, both.gm, static_cast<int>(both.afloat));
    std::printf("          %d panel(s) failed -> %d opening(s), %.3f m2 reaching a compartment\n",
                both.failedPanels, both.openings, both.openingArea);

    std::printf("\ngas     : the upper layer peaked at %.1f C and the steel it stood against"
                " at %.1f C, %.1f K behind it\n",
                both.peakGas - thermal::kCelsius, both.peakSteel - thermal::kCelsius,
                both.peakGas - both.peakSteel);
    std::printf("\naccount : %.4f GJ through the films, %.4f GJ held by the steel"
                " (residual %+.3e of scale)\n",
                both.filmHeat / 1e9, both.enthalpyGain / 1e9,
                std::abs(both.filmHeat) > 0
                    ? (both.filmHeat - both.enthalpyGain) / std::abs(both.filmHeat)
                    : 0.0);
    std::printf("          the fire booked %.4f GJ of boundary loss over the whole enclosure,"
                " against a peak of %.1f kW through the meshed bulkhead\n",
                both.gasWallLoss / 1e9, both.exchange / 1e3);
    std::printf("          worst film linearisation %.3e W, %.4f%% of that exchange -- one film"
                " per element row rather than one per layer\n",
                both.worstLinearisation,
                both.exchange > 0 ? 100.0 * both.worstLinearisation / both.exchange : 0.0);
    std::printf("cost    : %.2f s wall for the pair, %d coupling steps, %d factorisations\n",
                both.wallSeconds, both.couplingSteps, both.thermalFactorisations);
    std::printf("clocks  : the gas fell short of its %.1f s tick %d time(s) and capped its"
                " pressure solve %d time(s);\n          worst gap between the gas's clock and"
                " the printed time axis %.3e s\n",
                o.coupling, both.gasShortSteps, both.gasPressureCapped, both.gasClockGap);

    // --- The same state, decomposed ------------------------------------------------
    //
    // The three runs above are the control the milestone's sentence asks for, and
    // they differ in *two* things rather than one: a dry hold does not load the
    // bulkhead and it does not cool it either, so the fire-alone control runs
    // measurably hotter than the case it is a control for. This is the sharper
    // statement, because it removes one cause at a time from the state the run
    // actually reached.
    const StructuralMaterial steel = chain.scantlings.materials[0];
    if (both.failureTime >= 0) {
        thermal::HeatedMember at = chain.member(o.restraint);
        at.lateralMoment = both.failureMoment;
        const thermal::MemberState s = thermal::memberState(at, steel, both.failureTemperature);
        const double axialTerm = s.axialStress / s.columnCapacity;
        std::printf("\nat the state the bulkhead actually failed in -- %.1f C and %.2f kN m --"
                    " each cause alone:\n",
                    both.failureTemperature - thermal::kCelsius, both.failureMoment / 1e3);
        std::printf("  restrained expansion alone   %.4f of the member's capacity\n", axialTerm);
        std::printf("  the head alone, unmagnified  %.4f\n", s.additiveUtilisation - axialTerm);
        std::printf("  the head alone, magnified    %.4f   (x%.3f, and the axial load is what"
                    " magnifies it)\n",
                    s.utilisation - axialTerm, s.magnifier);
        std::printf("  together                     %.4f\n", s.utilisation);
    }

    // --- The window, measured off the runs themselves --------------------------------
    //
    // The restraint is the one number in the chain that none of the three subsystems
    // owns -- it is set by the stiffness of the bulkhead deck and the tank top, which
    // are in neither the fire model nor the conduction one -- so the run above is a
    // point on a curve. These are the two ends of the curve over which the
    // milestone's sentence is true, bisected against each run's own recorded states
    // rather than against a summary of them.
    const double upper = restraintWindowBound(fireOnly.history, steel,
                                              chain.section.modulusStiffener, chain.eulerStress);
    const double lower = restraintWindowBound(unrelieved.history, steel,
                                              chain.section.modulusStiffener, chain.eulerStress);
    std::printf("\nrestraint window: the milestone's sentence is true for %.4f <= r < %.4f,"
                " a factor of %.3f; this run used %.4f\n",
                lower, upper, upper > 0 && lower > 0 ? upper / lower : 0.0, o.restraint);
    std::printf("  below %.4f nothing fells the bulkhead at all; at or above %.4f the fire"
                " alone does and the head is not what fails it\n", lower, upper);
    std::printf("  the window is narrow because the two controls are not one change apart:"
                " a dry hold leaves the member %.1f K hotter (%.1f C against %.1f C), because"
                " the water that loads the bulkhead also cools it\n",
                fireOnly.peakMember - both.peakMember, fireOnly.peakMember - thermal::kCelsius,
                both.peakMember - thermal::kCelsius);

    thermal::HeatedMember axial = chain.member(1.0);
    thermal::HeatedMember loaded = axial;
    loaded.lateralMoment = both.failureMoment != 0.0 ? both.failureMoment : both.peakMoment;
    std::printf("\nthe temperatures behind it, at the moment the run developed (%.2f kN m):\n",
                loaded.lateralMoment / 1e3);
    std::printf("%10s %14s %14s %14s %10s\n", "restraint", "head alone", "fire alone",
                "both", "additive");
    for (double r : {1.0, 0.6, 0.4, 0.3, 0.25, 0.22, 0.2, 0.15, 0.1, 0.05}) {
        axial.restraint = r;
        loaded.restraint = r;
        const double tf = thermal::memberFailureTemperature(axial, steel);
        const double tc = thermal::memberFailureTemperature(loaded, steel);
        double ta = 0;
        for (double k = thermal::kCelsius + 20.0; k < thermal::kCelsius + 1200.0; k += 0.1)
            if (thermal::memberState(loaded, steel, k).additiveUtilisation >= 1.0) {
                ta = k;
                break;
            }
        std::printf("%10.3f %14s %14.2f %14.2f %10.2f\n", r, "never",
                    tf > 0 ? tf - thermal::kCelsius : -1.0, tc > 0 ? tc - thermal::kCelsius : -1.0,
                    ta > 0 ? ta - thermal::kCelsius : -1.0);
    }
    std::printf("the head alone never fells it at any temperature: the member is elastic"
                " under %.2f kN m and %.4f of yield\n",
                loaded.lateralMoment / 1e3,
                std::abs(loaded.lateralMoment) / chain.section.modulusStiffener /
                    steel.yieldStrength);

    if (o.sweep) {
        std::printf("\nrestraint sweep, the chain run end to end at each:\n");
        std::printf("%10s %10s %10s %10s %10s\n", "restraint", "fire only", "both", "members",
                    "hole (m2)");
        // Chosen to bracket both ends of the window the run above measures, so the
        // sweep is a check on that measurement rather than a picture beside it: the
        // fire alone has to start failing between 0.24 and 0.26, and the pair has to
        // stop failing between 0.20 and 0.22.
        Options s = o;
        for (double r : {0.40, 0.30, 0.26, 0.24, 0.22, 0.20, 0.15}) {
            s.restraint = r;
            const Outcome f = run(chain, s, true, false, false);
            const Outcome b = run(chain, s, true, true, false);
            std::printf("%10.3f %10s %10s %10d %10.2f\n", r,
                        f.failureTime >= 0 ? "fails" : "holds",
                        b.failureTime >= 0 ? "fails" : "holds", b.failedMembers, b.openingArea);
        }
    }

    // --- The verdict --------------------------------------------------------------
    //
    // Every one of these is a way the run could report success while proving
    // nothing, and each has been the actual failure of some earlier version of it.
    int bad = 0;
    const auto require = [&](bool ok, const char* what) {
        if (!ok) {
            std::printf("  ! %s\n", what);
            ++bad;
        }
    };
    require(both.peakSteel > kTAmbient + 100.0, "the bulkhead never got hot");
    require(both.peakHead > 1.0, "there was never a metre of water behind the bulkhead");
    require(both.peakMoment > 1.0e3, "the head never put a moment into a member");
    require(both.failureTime >= 0, "the bulkhead did not fail with both causes present");
    require(waterOnly.failureTime < 0, "the head alone felled the bulkhead");
    require(fireOnly.failureTime < 0, "the fire alone felled the bulkhead");
    require(both.openingArea > 0.1, "the failure opened no hole into a compartment");
    require(both.waterFireSide > 10.0, "no water reached the compartment that was on fire");
    require(both.compartmentsWet >= 3, "the flooding did not spread beyond the first compartment");
    require(both.failureAdditive < 1.0,
            "a purely additive check would have failed too: the coupling bought nothing");
    require(std::abs(both.filmHeat - both.enthalpyGain) <= 1e-6 * std::abs(both.filmHeat),
            "the steel's energy account does not close");
    // **The one way this tool can publish a wrong number while every check above
    // passes.** The gas, the steel and the ship are each stepped by one `o.coupling`
    // and printed under a `t` this file computes as `(step + 1) * o.coupling`; only
    // the gas can advance by less than it was asked for, and if it does, the three
    // columns of a row are readings at three different model times and the trace
    // stays perfectly smooth. `flip_probe` keeps the same guard over
    // `flip::StepResult::incomplete` and states it the same way.
    //
    // All three runs, not just `both`: a control that under-advanced is not a
    // control, and "the fire alone does not fell it" would then be a statement about
    // a fire that ran for less time than it was given.
    require(both.gasShortSteps == 0 && fireOnly.gasShortSteps == 0 &&
                waterOnly.gasShortSteps == 0,
            "a gas step was short of the time it was asked for: the three solves are no"
            " longer on one clock");
    // 1e-6 s against a 5 s coupling. The gap is the difference between an arithmetic
    // axis and a sum of substeps, so on a healthy run it is rounding and nothing
    // else -- measured 0 s over 720 couplings -- and the bound is six orders of
    // magnitude below the smallest shortfall a single dropped substep could produce.
    require(both.gasClockGap < 1e-6 && fireOnly.gasClockGap < 1e-6 &&
                waterOnly.gasClockGap < 1e-6,
            "the gas's clock and the printed time axis have come apart");
    // Published by `fire.hpp` "so that 'should never' can be asserted", and until now
    // asserted only by `tests/test_fire.cpp` and by none of its three tools.
    require(both.gasPressureCapped == 0 && fireOnly.gasPressureCapped == 0 &&
                waterOnly.gasPressureCapped == 0,
            "the gas pressure solve failed to bracket its root");
    require(fireOnly.peakSteel > kTAmbient + 100.0,
            "the fire-alone control was not a fire (vacuous survival)");
    require(waterOnly.peakMoment > 1.0e3,
            "the head-alone control had no head in it (vacuous survival)");

    if (bad != 0) {
        std::printf("\n%d check(s) failed\n", bad);
        return 1;
    }
    std::printf("\nok\n");
    return 0;
}
