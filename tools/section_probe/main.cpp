// SPDX-License-Identifier: MIT
//
// **Cut the reference ferry into a hold and reduce it.**
//
// `tests/test_section.cpp` checks the mesher against closed forms at unit scale and
// against `hullGirderSection` on two frame bays, and it has to stay inside the unit
// gate. This is the thing that runs at ship scale: eleven bays of the ferry, the
// resolution sweep that says what a Tier-1 section actually needs, and the
// Craig-Bampton reduction of the result.
//
// **Every *ship-scale* figure `docs/02-simulation.md` §3 publishes about the section
// mesher comes out of this program**, so re-running it is how a claim there is
// checked rather than quoted -- which is the failure mode this repo has hit three
// times, most recently with a factor of ten that three documents repeated from each
// other while the code beneath them had been right for weeks. The box-girder figures
// there -- the cost of welding a corner, the Bredt ratios, what the thickness taper
// costs -- come from `tests/test_section.cpp`, where they are assertions rather than
// output.
//
//   ./section_probe [--from=X] [--to=X] [--sub=N] [--sweep=N] [--reduce] [--modes=N]
//                   [--chain=N] [--match=M] [--scan=BAYS] [--invariance=BAYS]
//                   [--no-interface-ties] [--whole=N] [--whole-solve] [--profile=BAYS]
//                   [--wave[=AMPLITUDE]] [--wave-sections=N] [--no-frames]
//
// The last four modes are the whole-ship half and each returns on its own: `--whole`
// builds the entire 120 m as one piece and as two partitions and times both,
// `--profile` compares Tier 1 against Tier 0 along the length, and `--wave` drives the
// Tier-1 model with `girder.hpp`'s own wave loading. `--no-frames` strips every
// athwartships member, which is the control that separates the two tiers' answers --
// `hullGirderSection` drops those members itself, so Tier 0 cannot tell the difference
// and Tier 1 can.
#include "engine/sim/girder.hpp"
#include "engine/sim/reduction.hpp"
#include "engine/sim/scantlings.hpp"
#include "engine/sim/section.hpp"
#include "engine/sim/ship.hpp"
#include "engine/sim/waves.hpp"
#include "game/prototype/ferry.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <numbers>
#include <string>
#include <vector>

namespace {

double now() {
    using namespace std::chrono;
    return duration<double>(steady_clock::now().time_since_epoch()).count();
}

struct Options {
    // Eleven frame bays about midship. Frame stations are multiples of 2.4 m; the
    // watertight bulkheads are at -44, -38, -8, 20, 44 and none of them is, so
    // cutting "at a bulkhead" would put the plane through 188 panels. See
    // `section.hpp` and the `--from=-8 --to=20` run, which reports exactly that.
    double from = -7.2, to = 19.2;
    int subdivision = 1;
    int sweep = 3;       // refine to this subdivision; 0 skips the sweep
    bool reduce = true;  // build the substructure and reduce it
    int modes = 0;       // -1 takes ReduceParams' frequency cutoff
    // Cut [from, to] into this many sections, reduce each and assemble them, against
    // the same length in one piece. 0 skips it. It is off by default because the
    // assembled model is dense at (N+1) x 1 170 boundary DOF on this ship and the
    // Cholesky at the end of it is the most expensive thing this program does.
    int chain = 0;
    // Build the chain's interior cut planes with the in-plane line ties of
    // `section.hpp` §9 off. The negative control every figure §9 publishes is measured
    // against, and what this tool reported before they existed.
    bool interfaceTies = true;
    // How far apart two boundary DOF on a shared cut plane may be and still be the
    // same unknown. `matchBoundaries`' 1e-9 default used to be right amidships and to
    // leave a third of the interface unmatched at the ends -- see `section.hpp` §6
    // note 1 and the `--from=36 --to=48` run, which is what measured it. Since the
    // halo of §8 the two sections put that plane's nodes at the same double, so it is
    // not a tolerance any more and the default matches everything. It stays an option
    // because `--invariance` runs the mesher with the halo off, and that is the
    // configuration it was ever needed for.
    double match = 1e-9;
    // **The reach.** Slide a window this many frame bays long along the whole ship
    // and report, station by station, whether the mesher delivers something that
    // reduces and solves. 0 skips it. This is the measurement that says how much of
    // the ship a Tier-1 model can be built on, and before the collapsed-element
    // work (`section.hpp` §7) the answer was 21 windows of 49 -- 62.4 m of 120 in
    // five islands, whose longest unbroken run was 26.4 m.
    int scan = 0;
    // **The invariance.** At every frame station of the ship, mesh the window aft of
    // it, the window forward of it and the window spanning it, and ask whether the
    // three agree about the station they share. They did not: a node is the
    // mid-surface offset by t/2 along the nodal normal, both of which `buildSection`
    // averaged over the sub-quads *inside* the window, so the same node came out in a
    // different place from either side -- up to 1.0 mm at the bow -- and the
    // stiffener steel of a window depended on where it was cut, by 25.2% at the bow
    // shoulder. `SectionParams::halo` is the fix; this is the sweep that says so
    // along the whole hull rather than at the two stations the unit suite can afford.
    // 0 skips it. Meshing only: no solve, no reduction.
    int invariance = 0;
    // **The whole ship, both ways.** Mesh, reduce and solve the entire 120 m as one
    // section and as a chain of N, and time every stage of each. This is the model
    // `docs/02-simulation.md` says the tier exists to produce, and until it was built
    // once, end to end, nothing said what it costs. 0 skips it.
    int whole = 0;
    // Solve the monolith as well as meshing it. Off by default: the band is 5 384 over
    // 56 340 degrees of freedom and a banded factorisation of that is ninety times the
    // eleven-bay hold's, which is an hour rather than five seconds. It is the cost the
    // chain exists to avoid, so it is measured once and not gated.
    bool wholeSolve = false;
    // **The comparison against Tier 0 along the length.** Tile the hull with windows
    // this many frame bays long and report `A`, the neutral axis and `I` from both
    // tiers at every one of them. They agree amidships and they do not agree at the
    // ends, and the point of the mode is that the disagreement comes out as a
    // measurement with a mechanism rather than as a tolerance. 0 skips it.
    int profile = 0;
    // Drop every member with **no extent along x** -- the frames, the deck beams, the
    // bulkhead stiffeners. `hullGirderSection` already drops them, because an
    // athwartships member carries no longitudinal stress, so Tier 0's answer is
    // unchanged by this and Tier 1's is not: what is left is the transverse restraint
    // a ring of frames puts on the section's own Poisson contraction. `--profile`
    // runs both and reports the difference; this flag forces the reduced structure.
    bool frames = true;
    // **The hull-girder response.** Poise the ferry on a wave of her own length with
    // this amplitude in metres, take Tier 0's bending moment and its `M/Z` stress, and
    // ask the Tier-1 model the same question. 0 skips it; `--wave` alone means 3.0 m,
    // which is what `tests/test_promotion.cpp` uses and where the 84 MPa deck stress
    // `docs/02-simulation.md` publishes comes from.
    double wave = 0;
    // How many sections the `--wave` chain is cut into. Ten gives eleven cut planes --
    // enough to draw a deflection curve, and few enough that the dense assembled solve
    // is seconds rather than minutes.
    int waveSections = 10;
};

bool parse(int argc, char** argv, Options& o) {
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        const auto value = [&](const char* key) -> const char* {
            const std::string prefix = std::string("--") + key + "=";
            return a.rfind(prefix, 0) == 0 ? a.c_str() + prefix.size() : nullptr;
        };
        if (const char* v = value("from")) o.from = std::atof(v);
        else if (const char* v = value("to")) o.to = std::atof(v);
        else if (const char* v = value("sub")) o.subdivision = std::atoi(v);
        else if (const char* v = value("sweep")) o.sweep = std::atoi(v);
        else if (const char* v = value("modes")) o.modes = std::atoi(v);
        else if (const char* v = value("chain")) o.chain = std::atoi(v);
        else if (const char* v = value("match")) o.match = std::atof(v);
        else if (const char* v = value("scan")) o.scan = std::atoi(v);
        else if (const char* v = value("invariance")) o.invariance = std::atoi(v);
        else if (const char* v = value("whole")) o.whole = std::atoi(v);
        else if (const char* v = value("profile")) o.profile = std::atoi(v);
        else if (const char* v = value("wave")) o.wave = std::atof(v);
        else if (const char* v = value("wave-sections")) o.waveSections = std::atoi(v);
        else if (a == "--wave") o.wave = 3.0;
        else if (a == "--no-frames") o.frames = false;
        else if (a == "--whole-solve") o.wholeSolve = true;
        else if (a == "--no-interface-ties") o.interfaceTies = false;
        else if (a == "--no-reduce") o.reduce = false;
        else if (a == "--reduce") o.reduce = true;
        else {
            std::printf("unknown option %s\n", a.c_str());
            return false;
        }
    }
    return true;
}

// The hull's own extent, taken off the plating rather than off the ship, because
// this program only ever has the structure and the two need not agree at the stem.
void hullExtent(const sim::StructuralMesh& structure, double& lo, double& hi) {
    lo = 1e300;
    hi = -1e300;
    for (const sim::PlatePanel& p : structure.panels)
        for (int c = 0; c < 4; ++c) {
            lo = std::min(lo, p.corner[c].x);
            hi = std::max(hi, p.corner[c].x);
        }
}

// Every member with an extent along x, and nothing else.
//
// `hullGirderSection` already drops the rest -- an athwartships member carries no
// longitudinal stress, so a frame, a deck beam or a bulkhead stiffener is worth
// nothing to a beam idealisation and is *excluded geometrically* rather than by
// label. **Tier 0 is therefore identical on this structure and Tier 1 is not**, and
// the difference is one number: the transverse restraint a ring of frames puts on
// the section's own Poisson contraction, which raises `sigma_xx` above `E eps` on
// every strip the ring holds. It is a real stiffness that a beam has no way to
// represent, and this is the control that measures it instead of asserting it.
sim::StructuralMesh withoutTransverseMembers(const sim::StructuralMesh& source) {
    sim::StructuralMesh out = source;
    out.members.clear();
    for (const sim::StructuralMember& m : source.members)
        if (std::abs(m.b.x - m.a.x) > 1e-9) out.members.push_back(m);
    return out;
}

// Tier 0 over a **window** rather than at a point, and it is two different objects.
//
// A Tier-1 section spans frame bays and reports one `EA` for the whole length; the
// beam model varies continuously along it, so a comparison has to say which average
// it means -- and the two averages differ by a factor of two at the ends of this
// ship, which is most of what a naive sweep reads as disagreement.
//
//   * **The mean section**: `A`, `A z` and `A z^2 + I` averaged along the window and
//     recombined by the parallel axis theorem. It is what anyone writes down first
//     and it is *not* what a finite length of ship delivers.
//   * **The compliance mean**: the inverse of the mean of the section *compliance*.
//     A section with a plane-sections field prescribed on both cut planes carries a
//     constant `N` and a constant `M` along its length -- there is no load in
//     between -- so the generalised strains vary as `C(x) [N; M]` and what the two
//     planes prescribe is their **mean**. The stiffness that comes back is therefore
//     `(mean C)^-1`, a harmonic mean, and it is dominated by the softest station in
//     the window rather than by the average one. On the parallel middle body the two
//     agree to nothing at all; at the stem, where the section runs out over one bay,
//     the harmonic mean is less than half the arithmetic one.
//
// The 2x2 form is `[[A, A z], [A z, I + A z^2]]` about the **baseline**, so every
// window is combined in one frame and shifted afterwards; averaging neutral axes
// directly would be averaging ratios, which is not a section.
//
// Sampled at fractions of a bay that are never a frame seam, for the reason `--scan`
// gives: `hullGirderSection` asked exactly on a seam loses every panel either side of
// it -- 0.429 m^2 against 1.801 amidships -- on a floating-point knife edge in
// `sectionElements`' half-open `straddles` test. That is a Tier-0 defect this program
// declines to walk into and does not pretend to have fixed.
struct TierZeroWindow {
    bool ok = false;
    double x = 0;
    double area = 0, neutralAxis = 0, secondMoment = 0;              // the mean section
    double areaSeries = 0, neutralAxisSeries = 0, secondMomentSeries = 0;  // the compliance mean
    double zDeck = 0, zKeel = 0;
    double modulusDeck = 0, modulusKeel = 0;
    int samples = 0, severed = 0;  // stations where Tier 0 finds no section at all
};

TierZeroWindow tierZeroOverWindow(const sim::StructuralMesh& structure, double from, double bay,
                                  int bays, int perBay = 8) {
    TierZeroWindow out;
    const int n = std::max(bays, 0) * std::max(perBay, 1);
    if (n == 0) return out;
    out.x = from + 0.5 * bays * bay;
    double k00 = 0, k01 = 0, k11 = 0;    // the mean stiffness
    double c00 = 0, c01 = 0, c11 = 0;    // the mean compliance
    out.zKeel = 1e300;
    out.zDeck = -1e300;
    for (int i = 0; i < n; ++i) {
        const double x = from + bays * bay * (i + 0.5) / n;
        const sim::HullGirderSection one = sim::hullGirderSection(structure, x);
        ++out.samples;
        if (!(one.area > 0) || !(one.secondMoment > 0)) {
            ++out.severed;
            continue;
        }
        const double a = one.area, s = one.area * one.neutralAxis;
        const double j = one.secondMoment + one.area * one.neutralAxis * one.neutralAxis;
        k00 += a;
        k01 += s;
        k11 += j;
        const double det = a * j - s * s;
        if (!(det > 0)) {
            ++out.severed;
            continue;
        }
        c00 += j / det;
        c01 += -s / det;
        c11 += a / det;
        out.zKeel = std::min(out.zKeel, one.zKeel);
        out.zDeck = std::max(out.zDeck, one.zDeck);
    }
    const int good = out.samples - out.severed;
    if (good <= 0 || !(k00 > 0)) return out;
    k00 /= good;
    k01 /= good;
    k11 /= good;
    out.area = k00;
    out.neutralAxis = k01 / k00;
    out.secondMoment = k11 - k01 * k01 / k00;
    if (out.zDeck > out.neutralAxis) out.modulusDeck = out.secondMoment / (out.zDeck - out.neutralAxis);
    if (out.neutralAxis > out.zKeel) out.modulusKeel = out.secondMoment / (out.neutralAxis - out.zKeel);
    out.ok = true;

    // **A severed station makes the compliance mean infinite, and that is the answer
    // rather than a failure to compute one.** Tier 0 finding no section at all
    // somewhere inside the window says the beam is cut through there, and a cut beam
    // carries no axial force however stout it is on either side. It is left at zero
    // and counted, because a window like that is exactly where the two tiers have
    // something to say to each other.
    if (out.severed > 0) return out;
    c00 /= good;
    c01 /= good;
    c11 /= good;
    const double det = c00 * c11 - c01 * c01;
    if (!(det > 0)) return out;
    const double s00 = c11 / det, s01 = -c01 / det, s11 = c00 / det;
    out.areaSeries = s00;
    out.neutralAxisSeries = s01 / s00;
    out.secondMomentSeries = s11 - s01 * s01 / s00;
    return out;
}

// The second moment of the compliance-mean section about an axis that is **not** its
// own -- which is what `applyBeamLoad` reports, because the field it prescribes turns
// about the `reference` it was given rather than about whatever axis the section
// settles on. Shifting Tier 0 to the same axis is what makes the two comparable.
double secondMomentAbout(double area, double neutralAxis, double second, double about) {
    const double arm = neutralAxis - about;
    return second + area * arm * arm;
}

// What one Tier-1 section says about itself, in the three quantities Tier 0 offers an
// independent answer to. Three solves: a unit strain, a unit curvature and a unit
// twist. The neutral axis is *found* rather than assumed -- with a pure axial strain
// the reaction's first moment about `reference` is zero only if `reference` is the
// true axis -- so it is a third number and not a restatement of the first two.
struct TierOne {
    bool ok = false;
    double area = 0;         // m^2, EA / E
    double neutralAxis = 0;  // m
    double secondMoment = 0; // m^4, EI / E
    double torsion = 0;      // N m^2, GJ
    double residual = 0, restraint = 0;
    std::string problem;
};

TierOne measure(const sim::section::Section& piece, const sim::StructuralMaterial& material,
                double reference, bool wantTorsion = true) {
    TierOne out;
    if (piece.empty()) {
        out.problem = "empty section";
        return out;
    }
    sim::section::BeamLoad axial;
    axial.strain = 1e-6;
    axial.reference = reference;
    const sim::section::BeamResponse stretched =
        sim::section::applyBeamLoad(piece, material, axial);
    sim::section::BeamLoad bending;
    bending.curvature = 1e-6;
    bending.reference = reference;
    const sim::section::BeamResponse bent = sim::section::applyBeamLoad(piece, material, bending);
    if (!stretched.ok || !bent.ok) {
        out.problem = stretched.ok ? bent.problem : stretched.problem;
        return out;
    }
    out.area = stretched.axialStiffness / material.youngsModulus;
    out.neutralAxis = stretched.bendingMoment / stretched.axialForce + reference;
    out.secondMoment = bent.bendingStiffness / material.youngsModulus;
    out.residual = std::max(stretched.residual, bent.residual);
    out.restraint = std::max(stretched.restraintReaction, bent.restraintReaction);
    if (wantTorsion) {
        const sim::section::TorsionResponse twisted =
            sim::section::applyTwist(piece, material, 1e-6, reference);
        if (twisted.ok) out.torsion = twisted.torsionalStiffness;
    }
    out.ok = true;
    return out;
}

// Six assembled rows that take out the rigid body motions of a ship floating free.
//
// **Six and not three, which is the whole difference from `applyBeamLoad`.** A beam
// load prescribes both end planes, so only the three motions those planes leave are
// unheld; a ship on a wave is held by nothing at all. They are a *statically
// determinate* set -- three axes at one point kill the translations, two more at the
// most distant point kill the rotations about y and z, and one at a point off that
// line kills the rotation about x -- so on a balanced load their reaction is exactly
// zero, and reading it is the end-to-end check that the load balanced.
//
// A row an interior plane's line tie eliminated cannot restrain anything: it is
// isolated in the assembly, so holding it removes no motion and the factorisation
// reports a singular stiffness with no hint of where it came from. That is the same
// trap `section.cpp`'s `chainRestraints` records falling into.
//
// The three axes of one node come out of one `boundaryIdentity` entry, so they carry
// the same three doubles bit for bit and grouping them by exact position is right
// rather than lazy.
bool sixRestraints(const sim::section::Chain& chain, std::vector<std::uint32_t>& held) {
    held.clear();
    const std::vector<sim::reduction::BoundaryDof>& point = chain.assembly.boundaryPoint;
    if (point.size() != static_cast<std::size_t>(chain.assembly.boundary)) return false;
    const auto eliminated = [&](int b) {
        return std::binary_search(chain.planeTieDof.begin(), chain.planeTieDof.end(),
                                  static_cast<std::uint32_t>(b));
    };
    std::map<std::array<double, 3>, std::array<int, 3>> node;
    for (int b = 0; b < chain.assembly.boundary; ++b) {
        const sim::reduction::BoundaryDof& d = point[static_cast<std::size_t>(b)];
        if (d.axis > 2 || eliminated(b)) continue;
        std::array<int, 3>& slot =
            node.try_emplace({d.position.x, d.position.y, d.position.z},
                             std::array<int, 3>{-1, -1, -1}).first->second;
        slot[d.axis] = b;
    }
    const std::array<double, 3>* anchor = nullptr;
    const std::array<int, 3>* anchorRows = nullptr;
    for (const auto& [where, rows] : node)
        if (rows[0] >= 0 && rows[1] >= 0 && rows[2] >= 0) {
            anchor = &where;
            anchorRows = &rows;
            break;
        }
    if (anchor == nullptr) return false;

    const std::array<int, 3>* far = nullptr;
    double reach = 0;
    for (const auto& [where, rows] : node) {
        if (rows[1] < 0 || rows[2] < 0) continue;
        const double along = std::abs(where[0] - (*anchor)[0]);
        if (along > reach) {
            reach = along;
            far = &rows;
        }
    }
    if (far == nullptr || !(reach > 0)) return false;

    int lever = -1;
    double best = 0;
    for (const auto& [where, rows] : node) {
        const double dy = where[1] - (*anchor)[1], dz = where[2] - (*anchor)[2];
        const double arm = std::sqrt(dy * dy + dz * dz);
        // Hold whichever transverse axis the arm is **perpendicular** to: a deck node's
        // arm is nearly all in y, so holding its `u_y` would carry the rotation about x
        // through a lever of nothing.
        const int axis = std::abs(dy) >= std::abs(dz) ? 2 : 1;
        if (rows[axis] < 0 || !(arm > best) || rows[axis] == (*far)[1] || rows[axis] == (*far)[2] ||
            rows[axis] == (*anchorRows)[1] || rows[axis] == (*anchorRows)[2])
            continue;
        best = arm;
        lever = rows[axis];
    }
    if (lever < 0 || !(best > 0)) return false;

    held = {static_cast<std::uint32_t>((*anchorRows)[0]),
            static_cast<std::uint32_t>((*anchorRows)[1]),
            static_cast<std::uint32_t>((*anchorRows)[2]),
            static_cast<std::uint32_t>((*far)[1]), static_cast<std::uint32_t>((*far)[2]),
            static_cast<std::uint32_t>(lever)};
    std::sort(held.begin(), held.end());
    return std::unique(held.begin(), held.end()) == held.end();
}

}  // namespace

int main(int argc, char** argv) {
    Options options;
    if (!parse(argc, argv, options)) return 2;

    sim::Ship ferry = game::buildFerry();
    std::vector<std::string> problems;
    // `--no-frames` strips every athwartships member, and it applies to **every** mode
    // rather than only to `--profile`, because it is the one control that separates a
    // Tier-1 answer from a Tier-0 one: `hullGirderSection` drops those members itself,
    // so Tier 0 scores this structure identically and any move belongs to the mesh.
    const sim::StructuralMesh structure = [&] {
        sim::StructuralMesh built =
            sim::makeStructuralMesh(ferry.hull, sim::ferryScantlings(), &problems);
        return options.frames ? built : withoutTransverseMembers(built);
    }();
    const sim::StructuralMaterial material = sim::ah36Steel();
    const double youngs = material.youngsModulus;
    const double middle = 0.5 * (options.from + options.to);
    const sim::HullGirderSection girder = sim::hullGirderSection(structure, middle);

    std::printf("ferry structure: %zu panels, %zu members, frame spacing %.3f m\n",
                structure.panels.size(), structure.members.size(), structure.frameSpacing);
    std::printf("Tier 0 at x = %.2f m: A = %.5f m^2, neutral axis %.5f m, I = %.5f m^4\n",
                middle, girder.area, girder.neutralAxis, girder.secondMoment);
    std::printf("             so EA = %.5e N, EI = %.5e N m^2\n\n", youngs * girder.area,
                youngs * girder.secondMoment);

    // --- The invariance: a section is worth what it is worth wherever it was cut ----
    //
    // Slide along the hull a station at a time and, at each, mesh the window aft of
    // it, the window forward of it and the window spanning both. Two questions, and
    // they fail in different places:
    //
    //   * **Do the two halves put the shared plane's nodes in the same place?** They
    //     did not, wherever the hull turns: the nodal normal was averaged over the
    //     sub-quads inside the window, so it leaned aft for one section and forward
    //     for the other. Reported as the furthest any node on the plane is from its
    //     partner across it -- the same question `matchBoundaries` asks, asked
    //     without a tolerance.
    //   * **Does the spanning window carry the steel its two halves carry?** It did
    //     not, wherever the plating steps: the nodal thickness was averaged the same
    //     way, so a run of stiffener broke in the spanning window and not in either
    //     half. Reported as the relative disagreement in `memberMass`.
    //
    // `SectionParams::halo` averages both over one bay beyond each plane and meshes
    // only what is inside. The node gap then reads **exactly zero**, and "exactly" is
    // the claim rather than a rounding: the same node is reached by the same panels in
    // the same order whichever window contains it, so its normal and its thickness are
    // formed from the same terms and come out at the same double. The steel is the
    // same fibres summed in a different order and agrees to 3.9e-14 of itself, which
    // is why that one is counted against 1e-9 and the node against nothing at all.
    if (options.invariance > 0) {
        double lo = 1e300, hi = -1e300;
        for (const sim::PlatePanel& p : structure.panels)
            for (int c = 0; c < 4; ++c) {
                lo = std::min(lo, p.corner[c].x);
                hi = std::max(hi, p.corner[c].x);
            }
        const double bay = structure.frameSpacing;
        const double span = options.invariance * bay;
        std::printf("=== invariance: %g m windows either side of every station, halo on and off"
                    " ===\n", span);
        std::printf("  %8s %6s %7s %13s %13s %11s %11s\n", "station", "halo", "nodes", "worst gap m",
                    "steel halves", "one piece", "relative");
        const auto cut = [&](double from, double to, bool halo) {
            sim::section::SectionParams p;
            p.xFrom = from;
            p.xTo = to;
            p.subdivision = options.subdivision;
            p.junctions = false;  // the tie moves no node and this is about nodes
            p.halo = halo;
            return sim::section::buildSection(structure, p);
        };
        int stations = 0, moved[2] = {0, 0}, lost[2] = {0, 0}, withHalo = 0;
        double worstGap[2] = {0, 0}, worstSteel[2] = {0, 0};
        for (double x = std::floor(lo / bay + 0.5) * bay + span; x + span <= hi + 1e-9; x += bay) {
            ++stations;
            for (int halo = 1; halo >= 0; --halo) {
                const sim::section::Section aft = cut(x - span, x, halo != 0);
                const sim::section::Section forward = cut(x, x + span, halo != 0);
                const sim::section::Section whole = cut(x - span, x + span, halo != 0);
                if (aft.empty() || forward.empty() || whole.empty()) continue;
                if (halo == 1 && whole.haloPanels > 0) ++withHalo;
                // Every node of the aft section's forward plane, against the nearest
                // node of the forward section's aft plane. Nearest and not "the same
                // index": the two meshes are numbered for bandwidth and share no
                // numbering at all.
                const auto place = [](const sim::section::Section& s, std::uint32_t node) {
                    return sim::Vec3{s.mesh.position[node * 3], s.mesh.position[node * 3 + 1],
                                     s.mesh.position[node * 3 + 2]};
                };
                double gap = 0;
                for (std::uint32_t a : aft.forwardNodes) {
                    double nearest = 1e300;
                    for (std::uint32_t b : forward.aftNodes)
                        nearest = std::min(nearest, sim::length(place(aft, a) - place(forward, b)));
                    gap = std::max(gap, nearest);
                }
                const double halves = aft.memberMass + forward.memberMass;
                const double steel = halves > 0 ? std::abs(whole.memberMass / halves - 1.0) : 0.0;
                const auto slot = static_cast<std::size_t>(halo);
                worstGap[slot] = std::max(worstGap[slot], gap);
                worstSteel[slot] = std::max(worstSteel[slot], steel);
                if (gap > 1e-9) ++moved[slot];
                if (steel > 1e-9) ++lost[slot];
                if (halo == 0 && (gap > 1e-9 || steel > 1e-9))
                    std::printf("  %8.1f %6d %7zu %13.3e %13.1f %11.1f %+11.3e\n", x, halo,
                                aft.forwardNodes.size(), gap, halves, whole.memberMass,
                                halves > 0 ? whole.memberMass / halves - 1.0 : 0.0);
            }
            std::fflush(stdout);
        }
        std::printf("\n  %d stations swept at subdivision %d\n", stations, options.subdivision);
        for (int halo = 1; halo >= 0; --halo) {
            const auto slot = static_cast<std::size_t>(halo);
            std::printf("  halo %s: %3d stations move a node (worst %.3e m), %3d lose stiffener"
                        " steel (worst %.3e)\n",
                        halo ? "on " : "off", moved[slot], worstGap[slot], lost[slot],
                        worstSteel[slot]);
        }

        // --- The success contract -------------------------------------------------
        //
        // Zero, and two guards against zero being free. The first is that the halo
        // was made of something; the second is that the *same sweep without it* is
        // not zero, because a hull with no shoulder and no strake step would report a
        // clean sweep from a mesher that had never been fixed.
        if (moved[1] != 0 || lost[1] != 0) {
            std::printf("       ! %d stations put a node in a different place depending on the"
                        " window (worst %.3e m) and %d carry different stiffener steel (worst"
                        " %.3e): a section is not a property of the ship\n",
                        moved[1], worstGap[1], lost[1], worstSteel[1]);
            return 1;
        }
        if (withHalo != stations) {
            std::printf("       ! only %d of %d stations had any plating beyond the cut to"
                        " average, so the sweep is a statement about nothing\n",
                        withHalo, stations);
            return 1;
        }
        if (moved[0] == 0 || lost[0] == 0) {
            std::printf("       ! without the halo the same sweep moves %d nodes and loses steel"
                        " at %d stations. Both have to be non-zero or this ship cannot tell a"
                        " fixed mesher from an unfixed one\n",
                        moved[0], lost[0]);
            return 1;
        }
        if (stations < 8) {
            std::printf("       ! %d stations is not a ship\n", stations);
            return 1;
        }
        std::printf("\nok\n");
        return 0;
    }

    // --- The reach: how much of the ship can be meshed at all ----------------------
    //
    // Everything below this measures one section. This measures how many sections
    // there are: a window of `--scan` frame bays slid along the whole hull, asked
    // for a mesh, a reduction and two solves at every station. A window counts as
    // reached only when all four succeed *and* the mesh is one connected piece --
    // `EA` alone would score a mesh that joined nothing (see `section.hpp` §2), so
    // it is reported but is not the criterion.
    if (options.scan > 0) {
        double lo = 1e300, hi = -1e300;
        for (const sim::PlatePanel& p : structure.panels)
            for (int c = 0; c < 4; ++c) {
                lo = std::min(lo, p.corner[c].x);
                hi = std::max(hi, p.corner[c].x);
            }
        const double bay = structure.frameSpacing;
        const double span = options.scan * bay;
        std::printf("=== reach: a %g m window (%d bays) along a hull of %.1f m, subdivision %d ===\n",
                    span, options.scan, hi - lo, options.subdivision);
        std::printf("  %8s %8s %8s %7s %7s %6s %11s %9s %9s %10s %9s %8s\n", "xFrom", "xTo",
                    "elems", "collps", "invert", "comps", "minGaussJ", "A_eff m2", "Tier0 m2",
                    "GJ", "first Hz", "verdict");
        // Which bays are inside at least one working window. Consecutive windows
        // overlap by `span - bay`, so summing their lengths would count the middle
        // body several times over; a union is the only honest reach.
        const int bays = static_cast<int>(std::lround((hi - lo) / bay));
        std::vector<bool> covered(static_cast<std::size_t>(std::max(bays, 0)), false);
        int windows = 0, good = 0, usable = 0, collapsedSeen = 0;
        double firstGood = 1e300, lastGood = -1e300, bestTorsion = 0;
        for (double x = std::floor(lo / bay + 0.5) * bay; x + span <= hi + 1e-9; x += bay) {
            sim::section::SectionParams window;
            window.xFrom = x;
            window.xTo = x + span;
            window.subdivision = options.subdivision;
            const sim::section::Section piece = sim::section::buildSection(structure, window);
            ++windows;
            // **Tier 0 at the centre of a bay, not at the window's midpoint.** An
            // even number of bays puts that midpoint on a frame station, and
            // `hullGirderSection` sampled exactly on a panel seam loses every panel
            // either side of it: 0.42932 m^2 against 1.80133 amidships, a 76%
            // shortfall that is a knife-edge on the half-open `straddles` test in
            // `sectionElements` rather than anything about the ship. It is not this
            // file's defect to fix -- and it is why the column below reports both
            // areas rather than their ratio.
            const sim::HullGirderSection tier0 =
                sim::hullGirderSection(structure, x + 0.5 * (span - bay));
            double effective = 0, gj = 0, hz = 0;
            bool solved = false, reduces = false;
            if (!piece.empty()) {
                sim::section::BeamLoad axial;
                axial.strain = 1e-6;
                axial.reference = tier0.neutralAxis;
                const sim::section::BeamResponse stretched =
                    sim::section::applyBeamLoad(piece, material, axial);
                const sim::section::TorsionResponse twisted =
                    sim::section::applyTwist(piece, material, 1e-6, tier0.neutralAxis);
                solved = stretched.ok && twisted.ok;
                if (stretched.ok) effective = stretched.axialStiffness / youngs;
                if (twisted.ok) gj = twisted.torsionalStiffness;
                const sim::reduction::Substructure substructure(
                    piece.mesh, piece.material, piece.interfaceNodes, piece.attachment);
                reduces = substructure.ready();
                if (reduces) {
                    const sim::reduction::Eigenpairs modes = substructure.fixedInterfaceModes(1);
                    if (!modes.value.empty())
                        hz = std::sqrt(std::max(0.0, modes.value[0])) / (2.0 * std::numbers::pi);
                }
            }
            // Two verdicts, because they answer different questions. **Usable** is
            // what a Tier-1 model needs: it meshes, it reduces, it solves, and no
            // piece of it floats free of the interface (a floating component is a
            // mechanism `reduction::Substructure` does not catch -- `section.hpp`).
            // **One piece** is stricter and is what a chain wants; a section in two
            // spanning pieces still reduces and still carries load, it just has a
            // junction the mesher declined to close.
            const bool works = !piece.empty() && piece.invertedElements == 0 && solved && reduces &&
                               piece.floatingComponents == 0;
            const bool ok = works && piece.components == 1;
            if (works) ++usable;
            if (piece.collapsedElements > 0) ++collapsedSeen;
            bestTorsion = std::max(bestTorsion, gj);
            if (ok) {
                ++good;
                firstGood = std::min(firstGood, window.xFrom);
                lastGood = std::max(lastGood, window.xTo);
                for (int b = 0; b < options.scan; ++b) {
                    const auto index = static_cast<std::size_t>(std::lround((x - lo) / bay) + b);
                    if (index < covered.size()) covered[index] = true;
                }
            }
            std::printf("  %8.2f %8.2f %8zu %7d %7d %6d %11.3e %9.5f %9.5f %10.3e %9.4f %8s\n",
                        window.xFrom, window.xTo, piece.elementCount(), piece.collapsedElements,
                        piece.invertedElements, piece.components, piece.worstJacobian, effective,
                        tier0.area, gj, hz,
                        ok ? "ok" : (piece.empty() ? "empty" : (works ? "npieces" : "REFUSED")));
            std::fflush(stdout);
        }
        double reached = 0;
        for (bool b : covered)
            if (b) reached += bay;
        std::printf("\n  %d of %d windows mesh, reduce and solve with nothing floating (%.1f%%)\n",
                    usable, windows, 100.0 * usable / windows);
        std::printf("  %d of %d are also a single connected piece (%.1f%%)\n", good, windows,
                    100.0 * good / windows);
        std::printf("  outermost working cut planes x = %.1f .. %.1f m\n",
                    good > 0 ? firstGood : 0.0, good > 0 ? lastGood : 0.0);
        std::printf("  reach: %.1f m of a %.1f m ship (%.1f%%)\n", reached, hi - lo,
                    100.0 * reached / (hi - lo));

        // --- The success contract -------------------------------------------------
        //
        // The reach is the claim this mode exists to make, so it is the thing the
        // gate checks. **With two guards against it being vacuous**, because "every
        // window worked" is also what an empty ship, or a ship with no degenerate
        // panels on it, would report -- and the second of those is exactly the state
        // this code was in before the collapsed element was understood.
        if (usable != windows) {
            std::printf("       ! %d of %d windows do not mesh, reduce and solve\n",
                        windows - usable, windows);
            return 1;
        }
        if (reached < hi - lo - 0.5 * bay) {
            std::printf("       ! the mesher reaches %.1f m of a %.1f m ship\n", reached, hi - lo);
            return 1;
        }
        if (collapsedSeen == 0) {
            std::printf("       ! not one window contained a collapsed element, so the reach"
                        " measured nothing: the ship this ran on has no degenerate plating and"
                        " the whole comparison is vacuous\n");
            return 1;
        }
        if (windows < 8 || !(bestTorsion > 1e11)) {
            std::printf("       ! %d windows and a best GJ of %.3e: the sections carry no real"
                        " stiffness, so 'every window worked' is a statement about nothing\n",
                        windows, bestTorsion);
            return 1;
        }
        std::printf("\nok\n");
        return 0;
    }

    // --- The whole ship, as one piece and as a chain ---------------------------------
    //
    // Every part of this existed before today and nobody had put them together: the
    // mesher reaches all 120 m (§7), the halo makes a station's nodes a property of
    // the ship rather than of the cut (§8), and the in-plane line tie stops an
    // interior plane costing torsion (§9). What was never built is the model those
    // three exist to make.
    //
    // **Two costs, and they are not the same shape.** One piece is one mesh, one
    // banded factorisation and one Guyan condensation over an interface of two
    // planes. A chain of N is N of each, every one of them small, and then a *dense*
    // assembled model of `(N+1) x n_b` -- so cutting finer makes every piece cheaper
    // to reduce and the assembly quadratically dearer to hold and cubically dearer to
    // solve. Both halves are timed here rather than argued, because the trade is what
    // a caller is choosing whether or not anyone has stated it.
    if (options.whole > 0) {
        double lo = 0, hi = 0;
        hullExtent(structure, lo, hi);
        const double bay = structure.frameSpacing;
        const double from = std::floor(lo / bay + 0.5) * bay;
        const double to = std::floor(hi / bay + 0.5) * bay;
        const sim::HullGirderSection amidships = sim::hullGirderSection(structure, 0.5 * bay);
        std::printf("=== the whole ship: %.1f .. %.1f m, %d bays, subdivision %d ===\n", from, to,
                    static_cast<int>(std::lround((to - from) / bay)), options.subdivision);
        std::fflush(stdout);

        sim::section::SectionParams shipParams;
        shipParams.xFrom = from;
        shipParams.xTo = to;
        shipParams.subdivision = options.subdivision;

        // --- one piece ---------------------------------------------------------------
        double at = now();
        const sim::section::Section ship = sim::section::buildSection(structure, shipParams);
        const double shipMesh = now() - at;
        if (ship.empty()) {
            std::printf("  ! the whole ship did not mesh\n");
            for (const std::string& problem : ship.problems) std::printf("      ! %s\n", problem.c_str());
            return 1;
        }
        std::printf("  one piece   mesh %7.2f s  %zu elements, %zu nodes, band %zu, %d components"
                    " (%d spanning, %d floating)\n",
                    shipMesh, ship.elementCount(), ship.nodeCount(), ship.halfBandwidth,
                    ship.components, ship.spanningComponents, ship.floatingComponents);
        std::printf("  %-11s %.0f kg of plate + %.0f kg of member; %.1f m of junction edge, %.1f m"
                    " tied\n", "", ship.plateMass, ship.memberMass, ship.junctionEdges,
                    ship.tiedEdges);
        std::fflush(stdout);

        // **The monolith meshes in half a second and does not solve in an hour, and the
        // number that says why is the band.** `solidshell::solveStatic` factors a
        // banded system at `n * b^2`: the eleven-bay hold is 146 untied and 1 520 tied,
        // and the whole ship tied is **5 384** over 56 340 degrees of freedom, which is
        // 1.6e12 flops per factorisation against the hold's 1.8e10 -- a factor of
        // ninety. The band grows because a tie joins nodes no element edge joins and the
        // ordering has to carry the whole cross-section plus every junction across 50
        // bays; it is not a rounding and it is not a property of this ordering, it is
        // what a ship-length wavefront costs.
        //
        // So it is behind a flag. `--whole-solve` runs it, and the run it was measured
        // on is reported below rather than repeated in a gate.
        TierOne shipSays;
        double shipSolve = 0;
        if (options.wholeSolve) {
            at = now();
            shipSays = measure(ship, material, amidships.neutralAxis);
            shipSolve = now() - at;
            if (!shipSays.ok) {
                std::printf("  ! the whole ship did not solve: %s\n", shipSays.problem.c_str());
                return 1;
            }
            std::printf("  %-11s solve %6.2f s (three banded solves at band %zu) EA %.5e"
                        "  z_na %.5f  EI %.5e  GJ %.4e\n", "", shipSolve, ship.halfBandwidth,
                        youngs * shipSays.area, shipSays.neutralAxis,
                        youngs * shipSays.secondMoment, shipSays.torsion);
            at = now();
            const sim::reduction::Substructure whole(ship.mesh, ship.material, ship.interfaceNodes,
                                                     ship.attachment);
            const double shipAssemble = now() - at;
            if (whole.ready()) {
                sim::reduction::ReduceParams guyan;
                guyan.modes = 0;
                guyan.cutoffFrequency = 0;
                at = now();
                const sim::reduction::Reduction reduced =
                    sim::reduction::craigBampton(whole, guyan);
                std::printf("  %-11s reduce: %.2f s to assemble + %.2f s Guyan, %zu boundary +"
                            " %zu interior DOF, interior band %zu, %d reduced DOF\n", "",
                            shipAssemble, now() - at, whole.boundaryCount(),
                            whole.interiorCount(), whole.halfBandwidth(), reduced.size());
            } else {
                std::printf("  %-11s reduce: the substructure refused after %.2f s\n", "",
                            shipAssemble);
                for (const std::string& problem : whole.problems())
                    std::printf("      ! %s\n", problem.c_str());
            }
        } else {
            std::printf("  %-11s solve: skipped. A banded factorisation of this mesh is"
                        " %.2e flops (n = %zu, b = %zu) against %.2e for the eleven-bay hold;"
                        " --whole-solve runs it anyway\n", "",
                        3.0 * ship.nodeCount() * static_cast<double>(ship.halfBandwidth) *
                            static_cast<double>(ship.halfBandwidth),
                        3 * ship.nodeCount(), ship.halfBandwidth, 8000.0 * 1520.0 * 1520.0 * 3);
        }
        std::fflush(stdout);

        // --- and as a chain, at two partitions ----------------------------------------
        //
        // **Two partitions and not one, because the monolith is not affordable to solve
        // and a comparison needs a second opinion.** Consecutive sections share every
        // degree of freedom on the plane they meet at, static condensation is exact
        // there, and the in-plane line ties of §9 close the ring of junctions an interior
        // plane opens -- so a chain of N and a chain of 2N are the *same model* cut two
        // different ways, and they have to agree to the conditioning of two dense
        // Choleskys rather than to a truncation. Where they would not agree is exactly
        // where a cut plane costs something, which is what the whole of §9 is about.
        //
        // It is also where the trade lives. Every extra cut makes each piece cheaper to
        // reduce -- the interior shrinks and the band with it -- and makes the assembled
        // model quadratically larger to hold and cubically dearer to factor. Both halves
        // are timed.
        struct Built {
            int pieces = 0;
            double build = 0, mesh = 0, reduce = 0, assemble = 0, solve = 0;
            int size = 0, boundary = 0, components = 0;
            std::size_t unmatched = 0;
            double area = 0, neutralAxis = 0, secondMoment = 0, torsion = 0;
            double residual = 0, restraint = 0, tied = 0, edges = 0, gap = 0;
            int planeTieNodes = 0, disagreeing = 0;
            bool ok = false;
        };
        const auto buildAndSolve = [&](int pieces) {
            Built out;
            out.pieces = pieces;
            sim::section::ChainParams params;
            params.section = shipParams;
            params.reduce.modes = 0;
            params.reduce.cutoffFrequency = 0;
            params.matchTolerance = options.match;
            for (int i = 0; i <= pieces; ++i)
                params.station.push_back(from + (to - from) * i / pieces);
            const double start = now();
            const sim::section::Chain chain = sim::section::buildChain(structure, params);
            out.build = now() - start;
            out.mesh = chain.meshSeconds;
            out.reduce = chain.reduceSeconds;
            out.assemble = chain.assembleSeconds;
            out.size = chain.assembly.size();
            out.boundary = chain.assembly.boundary;
            out.components = chain.components;
            out.tied = chain.tiedEdges;
            out.edges = chain.junctionEdges;
            out.gap = chain.worstGap;
            out.planeTieNodes = chain.planeTieNodes;
            out.disagreeing = chain.planeTiesDisagreeing;
            for (std::size_t u : chain.unmatched) out.unmatched += u;
            std::printf("  chain of %-2d built %6.2f s (mesh %.2f, reduce %.2f, assemble %.2f):"
                        " %d pieces, %d boundary DOF, size %d (%.0f MB dense), %zu unmatched,"
                        " worst gap %.2e m\n", pieces, out.build, out.mesh, out.reduce,
                        out.assemble, out.components, out.boundary, out.size,
                        2.0 * out.size * out.size * 8.0 / 1048576.0, out.unmatched, out.gap);
            for (const std::string& problem : chain.problems)
                std::printf("      ! %s\n", problem.c_str());
            std::fflush(stdout);
            if (!chain.ready()) return out;

            sim::section::BeamLoad axial;
            axial.strain = 1e-6;
            axial.reference = amidships.neutralAxis;
            sim::section::BeamLoad bending;
            bending.curvature = 1e-6;
            bending.reference = amidships.neutralAxis;
            const double solving = now();
            const sim::section::BeamResponse stretched = sim::section::applyBeamLoad(chain, axial);
            out.solve = now() - solving;
            const sim::section::BeamResponse bent = sim::section::applyBeamLoad(chain, bending);
            const sim::section::TorsionResponse twisted =
                sim::section::applyTwist(chain, 1e-6, amidships.neutralAxis);
            if (!stretched.ok || !bent.ok || !twisted.ok) {
                std::printf("      ! %s %s %s\n", stretched.problem.c_str(), bent.problem.c_str(),
                            twisted.problem.c_str());
                return out;
            }
            out.area = stretched.axialStiffness / youngs;
            out.neutralAxis = stretched.bendingMoment / stretched.axialForce + amidships.neutralAxis;
            out.secondMoment = bent.bendingStiffness / youngs;
            out.torsion = twisted.torsionalStiffness;
            out.residual = stretched.residual;
            out.restraint = stretched.restraintReaction;
            out.ok = true;
            std::printf("  %-11s solve %6.2f s (dense) EA %.5e  z_na %.5f  EI %.5e  GJ %.4e;"
                        " residual %.2e N, restraint %.2e N; ties %.1f m of %.1f m, %d line-tie"
                        " nodes, %d planes disagreeing\n", "", out.solve, youngs * out.area,
                        out.neutralAxis, youngs * out.secondMoment, out.torsion, out.residual,
                        out.restraint, out.tied, out.edges, out.planeTieNodes, out.disagreeing);
            std::fflush(stdout);
            return out;
        };
        const Built coarse = buildAndSolve(options.whole);
        const Built fine = buildAndSolve(2 * options.whole);
        if (!coarse.ok || !fine.ok) {
            std::printf("       ! one of the two partitions did not solve\n");
            return 1;
        }
        std::printf("  the two partitions against each other: EA %+.3e  z_na %+.3e  EI %+.3e"
                    "  GJ %+.3e\n", fine.area / coarse.area - 1.0,
                    fine.neutralAxis / coarse.neutralAxis - 1.0,
                    fine.secondMoment / coarse.secondMoment - 1.0,
                    fine.torsion / coarse.torsion - 1.0);
        if (options.wholeSolve && shipSays.ok)
            std::printf("  and against the same ship in one piece: EA %+.3e  z_na %+.3e"
                        "  EI %+.3e  GJ %+.3e\n", coarse.area / shipSays.area - 1.0,
                        coarse.neutralAxis / shipSays.neutralAxis - 1.0,
                        coarse.secondMoment / shipSays.secondMoment - 1.0,
                        coarse.torsion / shipSays.torsion - 1.0);
        std::fflush(stdout);

        // --- The success contract -----------------------------------------------------
        //
        // **`EA` is checked loosely and the joins are checked tightly**, which is the
        // rule `section.hpp` §2 sets and which this mode is the largest possible instance
        // of: prescribing plane sections on the two ends of a 120 m ship makes every
        // continuous strip carry `sigma = E eps` whatever it is attached to, so a model
        // that joined nothing would score best on `EA` here. What sees the joins is `GJ`,
        // the component count, the unmatched interface and the tied edge.
        if (ship.components != 1) {
            std::printf("       ! the whole ship meshes as %d components, not one\n",
                        ship.components);
            return 1;
        }
        if (coarse.components != 1 || fine.components != 1 || coarse.unmatched != 0 ||
            fine.unmatched != 0) {
            std::printf("       ! the chains are %d and %d pieces with %zu and %zu unmatched\n",
                        coarse.components, fine.components, coarse.unmatched, fine.unmatched);
            return 1;
        }
        if (coarse.disagreeing != 0 || fine.disagreeing != 0) {
            std::printf("       ! %d and %d interior planes had their two sections derive"
                        " different ties\n", coarse.disagreeing, fine.disagreeing);
            return 1;
        }
        if (!(coarse.torsion > 1e11) || !(fine.torsion > 1e11)) {
            std::printf("       ! GJ came out %.3e and %.3e: a ship that carries no torsion has"
                        " not been joined, and EA would not have said so\n", coarse.torsion,
                        fine.torsion);
            return 1;
        }
        // **Measured at 2.4e-4 between five pieces and ten, and the first expectation
        // written here was 1e-6, which was wrong for a reason `section.hpp` §9 already
        // states.** On the box a cut plane's line tie and the monolith's face tie are the
        // *same constraint* -- the slave lands exactly on the master face's edge, where
        // the bilinear shape functions collapse onto the two nodes the line uses -- so
        // cutting the box anywhere costs 1e-11. On a real hull the slave lands a little
        // off that edge, so the two are slightly different constraints on the same
        // junction and every extra cut plane trades one for the other. §9 measures that
        // at 0.095% for a two-section chain of the hold; five extra planes in a 120 m
        // ship come to 0.024%, the same effect an order smaller and in the same direction
        // -- the finer chain is the softer one.
        //
        // So the bound is 1e-3: four orders below the 3.3% a *single open* cut plane used
        // to cost, and above the residual that closing it honestly leaves.
        const double torsionGap = std::abs(fine.torsion / coarse.torsion - 1.0);
        if (!(torsionGap < 1e-3)) {
            std::printf("       ! the two partitions differ by %+.3e in GJ. An interior cut plane"
                        " is *shared*, not prescribed, so these are the same model twice and it"
                        " has to close\n", torsionGap);
            return 1;
        }
        // And the tied edge has to reach the monolith's, or the line ties are not closing
        // what a cut opens -- which is the one thing about a chain that used to be worse
        // than one piece rather than equal to it.
        if (!(coarse.tied > 0.99 * ship.tiedEdges) || !(fine.tied > 0.99 * ship.tiedEdges)) {
            std::printf("       ! the chains tie %.1f m and %.1f m of junction edge against the"
                        " monolith's %.1f m. A cut plane that stays open is a ship carrying less"
                        " shear than the same ship in one piece\n", coarse.tied, fine.tied,
                        ship.tiedEdges);
            return 1;
        }
        std::printf("\nok\n");
        return 0;
    }

    // --- Along the length: Tier 1 against Tier 0 -------------------------------------
    //
    // **This is the comparison the tier was built to be checked by, and it had never
    // been run.** `hullGirderSection` sums `A`, `A z` and `A z^2` over a transverse
    // cut; `applyBeamLoad` prescribes a plane-sections field on a 3D mesh of the same
    // scantlings and reads the reaction. The two share no line of code, so where they
    // agree it means something -- and where they disagree the *disagreement* is the
    // finding, not an error to be driven to zero. A beam idealisation carries no shear
    // lag, no local buckling and no warping, and a reduced 3D model carries all three;
    // a difference tuned away would mean the finer model had stopped being finer.
    //
    // **Three things have to be got right before a difference means anything, and two
    // of them are accounting.**
    //
    //   1. The section is **short by the members it could not attach**. On this ship
    //      those are exactly the girders, which sit off the longitudinal spacing and so
    //      pass through no node of a mesh built from panel seams. `missedMemberArea`
    //      and `missedMemberSecondMoment` are what to subtract -- *both*, because the
    //      girders are 4.4% of her area and 5.3% of her second moment. Subtracting only
    //      the area made the two tiers look 5% apart amidships when they are 0.4%.
    //   2. A finite length of ship does **not** deliver the average section. Both cut
    //      planes are prescribed and nothing is loaded in between, so `N` and `M` are
    //      constant and the generalised strains vary as `C(x) [N; M]`: what comes back
    //      is `(mean C)^-1`, a harmonic mean of the stiffness, dominated by the softest
    //      station in the window. Amidships that is the same number; at the stem it is
    //      less than half of it. Comparing against the arithmetic mean reads as an 80%
    //      disagreement and is a comparison of two different quantities.
    //   3. What is left over is real, and the control that names it is
    //      `--no-frames`. `hullGirderSection` drops every member with no extent along x
    //      -- an athwartships member carries no longitudinal stress -- so Tier 0 scores
    //      that structure identically and Tier 1 does not. The difference is the
    //      transverse restraint a ring of frames puts on the section's own Poisson
    //      contraction, which raises `sigma_xx` above `E eps` on every strip it holds.
    //      **The FEM is right and the beam is wrong about it**, and it is worth a few
    //      tenths of a per cent.
    if (options.profile > 0) {
        double lo = 0, hi = 0;
        hullExtent(structure, lo, hi);
        const double bay = structure.frameSpacing;
        const double span = options.profile * bay;
        const sim::StructuralMesh bare = withoutTransverseMembers(structure);
        std::printf("=== along the length: %g m windows (%d bays) tiling a %.1f m hull ===\n", span,
                    options.profile, hi - lo);
        std::printf("  A0/I0 are `hullGirderSection` averaged over the window; A0s/I0s are the"
                    " same sections' *compliance* mean, which is the section-level upper bound on"
                    " what a finite length can carry.\n");
        std::printf("  d is against A0s less the members the mesh could not attach. I is taken"
                    " about Tier 0's own mean neutral axis, on both sides.\n");
        std::printf("  %7s %7s %6s %3s %8s %8s %8s %8s %7s %7s %7s %7s %8s %8s %8s %7s %9s\n",
                    "xFrom", "xTo", "elems", "cmp", "A1", "A0", "A0s", "A0s-ms", "dA %", "dAnofr",
                    "z_na1", "z_na0s", "I1", "I0", "I0s-ms", "dI %", "GJ");
        int windows = 0, agreeing = 0, severed = 0;
        double worstArea = 0, worstAreaX = 0, worstSecond = 0, worstSecondX = 0;
        double worstFrame = 0, meanFrame = 0, meanArea = 0, meanSecond = 0;
        double worstNaive = 0;
        double amidshipsArea = 0, amidshipsSecond = 0, endArea = 0, endSecond = 0;
        for (double x = std::floor(lo / bay + 0.5) * bay; x + span <= hi + 1e-9; x += span) {
            sim::section::SectionParams window;
            window.xFrom = x;
            window.xTo = x + span;
            window.subdivision = options.subdivision;
            const sim::section::Section piece = sim::section::buildSection(structure, window);
            const TierZeroWindow tier0 = tierZeroOverWindow(structure, x, bay, options.profile);
            const TierOne tier1 = measure(piece, material, tier0.neutralAxis);
            if (!tier1.ok || !tier0.ok) {
                std::printf("  %7.2f %7.2f %6zu %3d   -- %s\n", x, x + span, piece.elementCount(),
                            piece.components,
                            tier1.ok ? "Tier 0 has no section here" : tier1.problem.c_str());
                std::fflush(stdout);
                continue;
            }
            ++windows;
            // The same window on a structure with no athwartships members. Tier 0's
            // answer is identical on it by construction, so the whole of the move is
            // Tier 1's own.
            const sim::section::Section stripped = sim::section::buildSection(bare, window);
            const TierOne noFrames =
                measure(stripped, material, tier0.neutralAxis, /*wantTorsion=*/false);

            // Tier 0, on the same structure the mesh actually built, about the same
            // axis. The missed members come off the compliance mean as though they
            // were missing everywhere in the window, which is what they are.
            const double areaReference = tier0.areaSeries - piece.missedMemberArea;
            const double secondReference =
                tier0.severed > 0
                    ? 0.0
                    : secondMomentAbout(tier0.areaSeries, tier0.neutralAxisSeries,
                                        tier0.secondMomentSeries, tier0.neutralAxis) -
                          (piece.missedMemberSecondMoment -
                           2.0 * tier0.neutralAxis * piece.missedMemberFirstMoment +
                           tier0.neutralAxis * tier0.neutralAxis * piece.missedMemberArea);
            const bool comparable = tier0.severed == 0 && areaReference > 0 && secondReference > 0;
            if (!comparable) ++severed;
            const double dArea = comparable ? tier1.area / areaReference - 1.0 : 0.0;
            const double dSecond = comparable ? tier1.secondMoment / secondReference - 1.0 : 0.0;
            const double dNaive =
                tier0.area > 0 ? tier1.area / (tier0.area - piece.missedMemberArea) - 1.0 : 0.0;
            const double dFrames =
                comparable && noFrames.ok
                    ? noFrames.area / (tier0.areaSeries - stripped.missedMemberArea) - 1.0
                    : dArea;
            if (comparable) {
                if (std::abs(dArea) < 0.01 && std::abs(dSecond) < 0.02) ++agreeing;
                if (std::abs(dArea) > std::abs(worstArea)) { worstArea = dArea; worstAreaX = x; }
                if (std::abs(dSecond) > std::abs(worstSecond)) {
                    worstSecond = dSecond;
                    worstSecondX = x;
                }
                worstNaive = std::max(worstNaive, std::abs(dNaive));
                worstFrame = std::max(worstFrame, std::abs(dArea - dFrames));
                meanFrame += dArea - dFrames;
                meanArea += dArea;
                meanSecond += dSecond;
                if (std::abs(x + 0.5 * span) < 12.0) {
                    amidshipsArea = dArea;
                    amidshipsSecond = dSecond;
                }
                if (std::abs(x + 0.5 * span) > 0.4 * (hi - lo)) {
                    endArea = std::max(endArea, std::abs(dArea));
                    endSecond = std::max(endSecond, std::abs(dSecond));
                }
            }
            std::printf("  %7.2f %7.2f %6zu %3d %8.5f %8.5f %8.5f %8.5f %+7.3f %+7.3f %7.4f %7.4f"
                        " %8.4f %8.4f %8.4f %+7.3f %9.3e\n",
                        x, x + span, piece.elementCount(), piece.components, tier1.area, tier0.area,
                        tier0.areaSeries, areaReference, 100.0 * dArea, 100.0 * dFrames,
                        tier1.neutralAxis, tier0.neutralAxisSeries, tier1.secondMoment,
                        tier0.secondMoment, secondReference, 100.0 * dSecond, tier1.torsion);
            std::fflush(stdout);
        }
        const int compared = windows - severed;
        if (compared <= 0) {
            std::printf("       ! nothing was measured\n");
            return 1;
        }
        std::printf("\n  %d windows meshed and solved; %d comparable, %d contain a station where"
                    " Tier 0 finds no section at all\n", windows, compared, severed);
        std::printf("  %d agree inside 1%% in A and 2%% in I\n", agreeing);
        std::printf("  mean dA %+.3f%%, mean dI %+.3f%%; worst dA %+.3f%% at x = %.1f, worst dI"
                    " %+.3f%% at x = %.1f\n", 100.0 * meanArea / compared,
                    100.0 * meanSecond / compared, 100.0 * worstArea, worstAreaX,
                    100.0 * worstSecond, worstSecondX);
        std::printf("  against the *arithmetic* mean section instead, the worst window reads"
                    " %.1f%%. The two references differ by under a per cent everywhere, which is"
                    " the measured answer to the first hypothesis and it is no: the shortfall at"
                    " the ends is not section-level averaging\n", 100.0 * worstNaive);
        std::printf("  the transverse members are worth %+.3f%% of A on average and %.3f%% at"
                    " worst: the Poisson restraint a beam cannot carry\n",
                    100.0 * meanFrame / compared, 100.0 * worstFrame);
        std::printf("  amidships dA %+.3f%% dI %+.3f%%; over the outer fifths, |dA| up to %.3f%%"
                    " and |dI| up to %.3f%%\n", 100.0 * amidshipsArea, 100.0 * amidshipsSecond,
                    100.0 * endArea, 100.0 * endSecond);


        // --- Why the ends fall short, which is a measurement rather than an argument ----
        //
        // Two hypotheses, and one experiment separates them.
        //
        //   * **Section-level averaging**: a length is as stiff as the harmonic mean of
        //     its own stations' section properties. The `A0s` column has already answered
        //     this and the answer is no -- it is within a per cent of `A0` everywhere,
        //     because the section *total* barely changes over two bays even at the stem.
        //   * **Plane sections is the assumption, and it is worth more the further apart
        //     the planes are.** Tier 0 asserts plane sections at *every* station. A
        //     Tier-1 window asserts it at two, and lets the plating do what it likes in
        //     between -- so the longer the window, the more freedom it has and the softer
        //     it comes out. Where the hull is prismatic there is nothing to gain by that
        //     freedom and the two agree; where a girth band is closing to nothing and a
        //     strake is turning to meet the stem, there is a great deal.
        //
        // **The experiment is nested, so window content is held fixed and only the
        // spacing of the prescribed planes changes.** One eight-bay window, and against
        // it the same eight bays cut into k pieces whose stiffnesses are combined in
        // series. At k = 8 the planes are 2.4 m apart, which is nearly what Tier 0
        // assumes; at k = 1 they are 19.2 m apart. Series is exact for the comparison
        // because each piece carries the same constant force, and amidships is the
        // control: the parallel middle body has to show none of the fall.
        //
        // **It is not monotone in between, and that was predicted wrongly.** Prescribing
        // more planes can only stiffen the *same* partition -- but eight bays cut into
        // four pieces and into two are cut at different frames, so which station a plane
        // lands on matters as well as how far apart the planes are. At the bow the four-
        // and two-piece readings sit within 0.02 of each other and *below* the one-piece
        // reading. The claim is carried by the two ends of the sweep, 2.4 m against
        // 19.2 m, and not by the shape in between.
        std::printf("\n  the same eight bays with the prescribed planes at different spacings"
                    " (A_eff over Tier 0's own answer for the length):\n");
        std::printf("  %8s %10s %12s %12s %12s\n", "pieces", "spacing m", "stern", "amidships",
                    "bow");
        const double reach = 8 * bay;
        const double place[3] = {std::floor(lo / bay + 0.5) * bay,
                                 std::round(-0.5 * reach / bay) * bay,
                                 std::floor(hi / bay + 0.5) * bay - reach};
        double reference[3] = {0, 0, 0};
        for (int p = 0; p < 3; ++p) {
            const TierZeroWindow zero = tierZeroOverWindow(structure, place[p], bay, 8);
            sim::section::SectionParams window;
            window.xFrom = place[p];
            window.xTo = place[p] + reach;
            window.subdivision = options.subdivision;
            reference[p] =
                zero.areaSeries - sim::section::buildSection(structure, window).missedMemberArea;
        }
        double nested[4][3] = {};
        const int pieces[4] = {8, 4, 2, 1};
        for (int k = 0; k < 4; ++k) {
            for (int p = 0; p < 3; ++p) {
                const int count = pieces[k];
                double compliance = 0;
                bool ok = true;
                for (int i = 0; i < count && ok; ++i) {
                    sim::section::SectionParams window;
                    window.xFrom = place[p] + reach * i / count;
                    window.xTo = place[p] + reach * (i + 1) / count;
                    window.subdivision = options.subdivision;
                    const sim::section::Section part =
                        sim::section::buildSection(structure, window);
                    const TierZeroWindow zero =
                        tierZeroOverWindow(structure, window.xFrom, bay, 8 / count);
                    const TierOne one = measure(part, material, zero.neutralAxis, false);
                    if (!one.ok || !(one.area > 0)) ok = false;
                    else compliance += 1.0 / one.area;
                }
                if (ok && compliance > 0 && reference[p] > 0)
                    nested[k][p] = (count / compliance) / reference[p];
            }
            std::printf("  %8d %10.1f %12.4f %12.4f %12.4f\n", pieces[k], reach / pieces[k],
                        nested[k][0], nested[k][1], nested[k][2]);
            std::fflush(stdout);
        }
        const double sternFall = nested[0][0] - nested[3][0];
        const double middleFall = std::abs(nested[0][1] - nested[3][1]);
        const double bowFall = nested[0][2] - nested[3][2];
        std::printf("  moving the prescribed planes from 2.4 m apart to 19.2 m costs the stern"
                    " %.3f and the bow %.3f of Tier 0's own answer, and the middle body %.4f."
                    " Plane sections is the assumption, and the ends are where it is worth"
                    " something\n", sternFall, bowFall, middleFall);
        // --- The success contract -----------------------------------------------------
        //
        // **What is asserted is the accounting and the mechanism, not the agreement.**
        // The two tiers are allowed to differ and the mode exists to say by how much,
        // so what has to hold is: amidships, where the hull is prismatic and a beam is
        // exactly right, the difference is small; the control that explains what is
        // left is not itself zero; and the ends are reached and are *not* like the
        // middle, because a sweep that came back uniform would be a sweep that never
        // left the parallel middle body.
        if (windows < 8) {
            std::printf("       ! %d windows is not a ship\n", windows);
            return 1;
        }
        if (!(std::abs(amidshipsArea) < 0.01) || !(std::abs(amidshipsSecond) < 0.02)) {
            std::printf("       ! amidships the two tiers differ by %+.3f%% in A and %+.3f%% in I."
                        " That is the parallel middle body, where a beam idealisation is exactly"
                        " right and an unexplained difference is a defect\n",
                        100.0 * amidshipsArea, 100.0 * amidshipsSecond);
            return 1;
        }
        if (!(worstFrame > 1e-4)) {
            std::printf("       ! removing every transverse member moved A by %.2e. The Poisson"
                        " restraint is the mechanism this mode names for the residual difference,"
                        " and a control worth nothing explains nothing\n", worstFrame);
            return 1;
        }
        if (!(endArea > 2.0 * std::abs(amidshipsArea)) &&
            !(endSecond > 2.0 * std::abs(amidshipsSecond))) {
            std::printf("       ! the ends agree with the beam as well as the middle body does"
                        " (%.3f%% against %.3f%% in A). Either the hull has no ends or this sweep"
                        " is not reaching them\n", 100.0 * endArea, 100.0 * amidshipsArea);
            return 1;
        }
        // **The sweep is the mechanism and the middle body is what stops it being a
        // story.** A shortfall that grows with window length is longitudinal continuity;
        // one that does not is section-level averaging, which the `A0s` column has
        // already ruled out. If the parallel middle body drifted under the same sweep
        // the trend would be the mesher's and not the ship's, so both halves are checked.
        if (!(sternFall > 0.1) || !(bowFall > 0.1)) {
            std::printf("       ! moving the prescribed planes from 2.4 m apart to 19.2 m costs"
                        " the stern %.3f and the bow %.3f. A shortfall that does not depend on"
                        " where plane sections is asserted is section-level averaging, and that"
                        " hypothesis is already measured and rejected\n", sternFall, bowFall);
            return 1;
        }
        if (!(middleFall < 0.02)) {
            std::printf("       ! the parallel middle body moved %.4f under the same nesting. The"
                        " trend at the ends is then a property of the mesher rather than of the"
                        " ship, and the mechanism above is not established\n", middleFall);
            return 1;
        }
        std::printf("\nok\n");
        return 0;
    }

    // --- The hull-girder response: the two tiers asked the same question --------------
    //
    // `girder.hpp` poises the ship on a wave, integrates weight minus buoyancy into a
    // shear and a bending moment, and divides by a section modulus. **Everything about
    // that is a beam**: plane sections remain plane, the load reaches the section
    // through nothing in particular, and `sigma` is linear in `z` and constant across
    // the breadth. This drives the Tier-1 model with the *same* load and asks what it
    // says instead. It is the first time the two tiers have been asked one question.
    //
    // **The load is applied, not the answer.** Handing the Tier-1 model Tier 0's own
    // bending moment would be assuming most of what is being compared, so what is
    // handed over is the distribution: the net force per station, weight minus
    // buoyancy, spread over the elements in that station's slab in proportion to their
    // volume. The resultant per station is then exact -- which is the whole of what
    // sets `V(x)` and `M(x)` -- and the local distribution is not, because where a
    // buoyancy pressure lands on the shell is a Tier-2 question and this is not one.
    //
    // **Guyan is exact here and it is worth knowing why, because it looks as though it
    // should not be.** `reduction.hpp` property 1: with zero modes the *boundary*
    // response is the same solve as eliminating the interior from the full system, for
    // **any** load, including one applied inside. So the deflection of every cut plane
    // is exact. What zero modes does not buy is the interior recovery under an interior
    // load -- `u_i = Psi u_b` misses `K_ii^-1 f_i` -- so the stress *inside* a section
    // is not read off the reduced model. It is recovered by prescribing the chain's own
    // (exact) interface displacement on that section's mesh and solving it directly,
    // which is exact and costs one banded factorisation.
    //
    // **Three restraints per piece is not enough here and that is the difference from
    // `applyBeamLoad`.** A beam load prescribes both end planes, so only three motions
    // are left; a floating ship is held by nothing at all and has six. They are
    // statically determinate, so on a balanced load their reaction is zero -- and that
    // reading is the end-to-end check that the load really was in equilibrium.
    if (options.wave > 0) {
        sim::Ship afloat = game::buildFerry();
        afloat.initialise(0.0);
        const double shipLength = afloat.hullHi.x - afloat.hullLo.x;
        // A wave of the ship's own length with the crest amidships: the classical
        // standard hogging condition, and the same construction `test_girder.cpp` and
        // `test_promotion.cpp` use. The crest is *found* by scanning a period rather
        // than assumed to be at t = 0.
        const double omega = std::sqrt(sim::kGravity * 2.0 * std::numbers::pi / shipLength);
        const sim::WaveField field = sim::WaveField::regular(options.wave, omega, 0.0);
        sim::Sea sea;
        sea.waves = &field;
        const double period = 2.0 * std::numbers::pi / omega;
        double crest = -1e30;
        for (int i = 0; i < 720; ++i) {
            const double t = period * i / 720.0;
            const double eta = field.elevation(afloat.state.position.x, 0.0, t);
            if (eta > crest) {
                crest = eta;
                sea.time = t;
            }
        }
        if (!sim::balanceOnWave(afloat, sea)) {
            std::printf("  ! she will not balance on this wave\n");
            return 1;
        }

        // Fifty stations and not fifty-one. The hull is 120 m on a 2.4 m frame
        // spacing, so fifty-one stations would put *every* station on a panel seam,
        // where `hullGirderSection` falls through the knife edge in `sectionElements`'
        // half-open `straddles` test and reports 0.429 m^2 against 1.801. At fifty the
        // spacing is 120/49 and the only stations on a seam are the two
        // perpendiculars, where there is no structure to lose.
        const int stationCount = 50;
        const std::vector<double> sx = sim::girderStations(afloat, stationCount);
        std::vector<double> weight = sim::weightDistribution(afloat, sx);
        const std::vector<double> buoyancy = sim::buoyancyDistribution(afloat, sea, sx);
        std::vector<double> width(sx.size(), 0.0);
        for (std::size_t i = 0; i < sx.size(); ++i) {
            const double left = i == 0 ? sx[0] : 0.5 * (sx[i - 1] + sx[i]);
            const double right = i + 1 == sx.size() ? sx.back() : 0.5 * (sx[i] + sx[i + 1]);
            width[i] = right - left;
        }
        const sim::HullGirder raw = sim::integrateGirder(sx, weight, buoyancy);
        std::printf("=== the hull girder on a %.1f m wave of her own length, crest amidships ===\n",
                    options.wave);
        std::printf("  Tier 0: %.0f t displacement, %s, peak M %.4e N m at x = %.1f, peak V %.3e N"
                    " at x = %.1f\n", raw.totalWeight / sim::kGravity / 1000.0,
                    raw.hogging() ? "hogging" : "sagging", raw.maxMoment, raw.maxMomentX,
                    raw.maxShear, raw.maxShearX);
        std::printf("  closure at the forward perpendicular: shear %+.3e, moment %+.3e of the"
                    " reference\n", raw.shearClosure, raw.momentClosure);
        for (const std::string& problem : sim::validateGirder(raw))
            std::printf("      ! %s\n", problem.c_str());
        std::fflush(stdout);

        // --- The Tier-1 model ---------------------------------------------------------
        double lo = 0, hi = 0;
        hullExtent(structure, lo, hi);
        const double bay = structure.frameSpacing;
        const double from = std::floor(lo / bay + 0.5) * bay;
        const double to = std::floor(hi / bay + 0.5) * bay;
        sim::section::ChainParams chainParams;
        chainParams.section.subdivision = options.subdivision;
        chainParams.reduce.modes = 0;
        chainParams.reduce.cutoffFrequency = 0;
        chainParams.matchTolerance = options.match;
        for (int i = 0; i <= options.waveSections; ++i)
            chainParams.station.push_back(from + (to - from) * i / options.waveSections);
        const double built = now();
        const sim::section::Chain chain = sim::section::buildChain(structure, chainParams);
        std::printf("  Tier 1: a chain of %d over %.1f .. %.1f m built in %.2f s -- %d pieces,"
                    " %d assembled DOF, %zu unmatched\n", options.waveSections, from, to,
                    now() - built, chain.components, chain.assembly.size(),
                    [&] { std::size_t u = 0; for (std::size_t n : chain.unmatched) u += n; return u; }());
        for (const std::string& problem : chain.problems) std::printf("      ! %s\n", problem.c_str());
        if (!chain.ready() || chain.components != 1) {
            std::printf("       ! the chain is not one ready piece\n");
            return 1;
        }
        std::fflush(stdout);

        // --- Where the load lands -----------------------------------------------------
        struct Bit {
            std::size_t section = 0, element = 0;
            int slab = 0;
            double volume = 0;
        };
        std::vector<Bit> bits;
        std::vector<double> slabVolume(sx.size(), 0.0);
        const double dx = sx.size() > 1 ? sx[1] - sx[0] : 1.0;
        for (std::size_t s = 0; s < chain.section.size(); ++s) {
            const sim::solidshell::HexMesh& mesh = chain.section[s].mesh;
            double nodes[sim::solidshell::kDof], volume[sim::solidshell::kGauss];
            for (std::size_t e = 0; e < mesh.elementCount(); ++e) {
                mesh.gather(e, mesh.position, nodes);
                sim::solidshell::gaussVolumes(nodes, volume);
                double total = 0, centre = 0;
                for (int g = 0; g < sim::solidshell::kGauss; ++g) total += volume[g];
                for (int n = 0; n < sim::solidshell::kNodes; ++n) centre += nodes[n * 3];
                centre /= sim::solidshell::kNodes;
                const int slab = std::clamp(static_cast<int>(std::lround((centre - sx.front()) / dx)),
                                            0, static_cast<int>(sx.size()) - 1);
                bits.push_back({s, e, slab, total});
                slabVolume[static_cast<std::size_t>(slab)] += total;
            }
        }

        // The net upward force each station's slab carries. Weight minus buoyancy is
        // the *downward* load `girder.hpp` integrates, so the force on the structure is
        // its negative.
        std::vector<double> force(sx.size(), 0.0);
        for (std::size_t i = 0; i < sx.size(); ++i)
            force[i] = -(weight[i] - buoyancy[i]) * width[i];

        // **Balanced before it is applied, and how much had to be removed is reported
        // rather than absorbed.** She floats and she floats level, so the exact load
        // sums to zero force and zero moment; what a quadrature over fifty slabs leaves
        // is a residual, and on a *free-free* structure a residual is not a small error
        // -- it is carried by whatever six rows hold the ship and reappears as a
        // bending moment that is not physical. The correction is spread in proportion
        // to steel volume, which is the same rule the load itself uses.
        double f0 = 0, f1 = 0, v0 = 0, v1 = 0, v2 = 0, scale = 0;
        for (std::size_t i = 0; i < sx.size(); ++i) {
            f0 += force[i];
            f1 += sx[i] * force[i];
            v0 += slabVolume[i];
            v1 += sx[i] * slabVolume[i];
            v2 += sx[i] * sx[i] * slabVolume[i];
            scale += std::abs(force[i]);
        }
        const double det = v0 * v2 - v1 * v1;
        if (!(det > 0) || !(scale > 0)) {
            std::printf("       ! the structure has no volume to balance against\n");
            return 1;
        }
        const double alpha = (f0 * v2 - f1 * v1) / det, beta = (v0 * f1 - v1 * f0) / det;
        double removed = 0;
        for (std::size_t i = 0; i < sx.size(); ++i) {
            const double correction = (alpha + beta * sx[i]) * slabVolume[i];
            force[i] -= correction;
            removed += std::abs(correction);
            weight[i] += correction / width[i];  // so Tier 0 integrates the same load
        }
        const sim::HullGirder tier0 = sim::integrateGirder(sx, weight, buoyancy);
        std::printf("  the discrete load was out of balance by %.3e N and %.3e N m; rebalancing"
                    " moved %.3f%% of it, and Tier 0 is re-integrated from the balanced load"
                    " (peak M %.4e, closure %+.2e)\n", f0, f1, 100.0 * removed / scale,
                    tier0.maxMoment, tier0.momentClosure);
        std::fflush(stdout);

        // --- Reduce it and solve ------------------------------------------------------
        std::vector<std::vector<double>> meshLoad(chain.section.size());
        for (std::size_t s = 0; s < chain.section.size(); ++s)
            meshLoad[s].assign(chain.section[s].mesh.nodeCount() * 3, 0.0);
        {
            double nodes[sim::solidshell::kDof], share[sim::solidshell::kNodes];
            for (const Bit& bit : bits) {
                if (!(slabVolume[static_cast<std::size_t>(bit.slab)] > 0) || !(bit.volume > 0))
                    continue;
                const sim::solidshell::HexMesh& mesh = chain.section[bit.section].mesh;
                mesh.gather(bit.element, mesh.position, nodes);
                sim::solidshell::elementMass(nodes, 1.0, share);
                const double f = force[static_cast<std::size_t>(bit.slab)] * bit.volume /
                                 slabVolume[static_cast<std::size_t>(bit.slab)];
                double sum = 0;
                for (int n = 0; n < sim::solidshell::kNodes; ++n) sum += share[n];
                if (!(sum > 0)) continue;
                for (int n = 0; n < sim::solidshell::kNodes; ++n) {
                    const std::uint32_t node =
                        mesh.index[bit.element * sim::solidshell::kNodes + static_cast<std::size_t>(n)];
                    meshLoad[bit.section][static_cast<std::size_t>(node) * 3 + 2] +=
                        f * share[n] / sum;
                }
            }
        }
        // **A junction tie's slave has no row, so a load left on it is a load thrown
        // away.** `reduction::reduceLoad` reads the boundary and interior partitions and
        // an eliminated degree of freedom is in neither, so it would silently vanish --
        // it belongs to the masters by the transpose of the constraint, exactly as the
        // reaction does in `section.cpp`'s `stiffnessTimes`. No master is ever a slave,
        // so one pass moves everything.
        double folded = 0;
        for (std::size_t s = 0; s < chain.section.size(); ++s)
            for (const sim::solidshell::Mpc& mpc : chain.section[s].attachment.constrained) {
                const double carried = meshLoad[s][mpc.slave];
                if (carried == 0.0) continue;
                folded += std::abs(carried);
                for (std::size_t a = 0; a < mpc.master.size(); ++a)
                    meshLoad[s][mpc.master[a]] += mpc.weight[a] * carried;
                meshLoad[s][mpc.slave] = 0.0;
            }

        std::vector<double> load(static_cast<std::size_t>(chain.assembly.size()), 0.0);
        for (std::size_t s = 0; s < chain.section.size(); ++s) {
            const std::vector<double> reduced =
                sim::reduction::reduceLoad(chain.substructure[s], chain.reduced[s], meshLoad[s]);
            const std::vector<int>& to_ = chain.assembly.from[s];
            for (std::size_t j = 0; j < reduced.size() && j < to_.size(); ++j)
                if (to_[j] >= 0) load[static_cast<std::size_t>(to_[j])] += reduced[j];
        }
        // And the same fold over the interior planes' line ties, whose rows the
        // assembly has already eliminated.
        for (const sim::solidshell::Mpc& mpc : chain.planeTies) {
            const double carried = load[mpc.slave];
            if (carried == 0.0) continue;
            for (std::size_t a = 0; a < mpc.master.size(); ++a)
                load[mpc.master[a]] += mpc.weight[a] * carried;
            load[mpc.slave] = 0.0;
        }

        std::vector<std::uint32_t> held;
        if (!sixRestraints(chain, held)) {
            std::printf("       ! no determinate set of six restraints could be found\n");
            return 1;
        }
        std::vector<std::uint32_t> solved = held;
        solved.insert(solved.end(), chain.planeTieDof.begin(), chain.planeTieDof.end());
        std::sort(solved.begin(), solved.end());
        solved.erase(std::unique(solved.begin(), solved.end()), solved.end());

        std::vector<double> state;
        std::string problem;
        const double solving = now();
        if (!sim::reduction::assembledStaticSolve(chain.assembly, load, solved, {}, state,
                                                  &problem)) {
            std::printf("       ! the loaded ship did not solve: %s\n", problem.c_str());
            return 1;
        }
        const double solveSeconds = now() - solving;
        // The reaction at the six restraints. Zero is the right answer and not a
        // tolerance: they select one of a family of zero-energy motions rather than
        // resisting anything, so anything here is load that did not balance.
        double restraint = 0, biggest = 0;
        {
            const std::size_t n = static_cast<std::size_t>(chain.assembly.size());
            for (std::uint32_t d : held) {
                double sum = 0;
                for (std::size_t j = 0; j < n; ++j)
                    sum += chain.assembly.stiffness[static_cast<std::size_t>(d) * n + j] * state[j];
                restraint = std::max(restraint, std::abs(sum - load[d]));
            }
            for (double f : load) biggest = std::max(biggest, std::abs(f));
        }
        // The interior planes' eliminated rows, filled in from their masters, so
        // everything below reads a displacement and not the zero the elimination left.
        for (const sim::solidshell::Mpc& mpc : chain.planeTies) {
            double sum = 0;
            for (std::size_t a = 0; a < mpc.master.size(); ++a)
                sum += mpc.weight[a] * state[mpc.master[a]];
            state[mpc.slave] = sum;
        }
        std::printf("  the reduced ship solved in %.2f s; %.1f kN of junction-tie load was folded"
                    " onto masters; the six restraints carry %.3e N against a largest applied"
                    " %.3e N\n", solveSeconds, folded / 1000.0, restraint, biggest);
        std::fflush(stdout);

        // --- What she does: the deflection, which is exact at every cut plane ----------
        //
        // Against the beam's own `w'' = -M / EI`, integrated over the same stations with
        // `I` from `hullGirderSection`. Both curves are defined only up to a rigid heave
        // and trim -- the ship is floating -- so the best-fit line is removed from each
        // before they are compared, and what is left is the *bending shape* and nothing
        // else. Removing it from only one of them would compare a shape with a shape
        // plus a rotation.
        // **And against a second beam, which is what separates the two causes.** A
        // difference here can only be `EI(x)` or shear: those are the two things a
        // Bernoulli integration of `M / EI` leaves out, and estimating them instead of
        // separating them is the mistake CLAUDE.md records costing a factor of ten. So
        // the same integration is run twice -- once with Tier 0's `I` and once with the
        // *section's own* `EI`, measured by `applyBeamLoad` on each piece of the chain --
        // and what survives the second is shear deflection and nothing else.
        std::vector<double> tierOneSecond(chain.section.size(), 0.0);
        for (std::size_t i = 0; i < chain.section.size(); ++i) {
            const double middle = 0.5 * (chain.section[i].xFrom + chain.section[i].xTo);
            const sim::HullGirderSection here = sim::hullGirderSection(structure, middle);
            const TierOne says = measure(chain.section[i], material, here.neutralAxis, false);
            if (says.ok) tierOneSecond[i] = says.secondMoment;
        }
        const auto integrate = [&](const std::vector<double>& kappa) {
            std::vector<double> slope(sx.size(), 0.0), deflection(sx.size(), 0.0);
            for (std::size_t i = 1; i < sx.size(); ++i) {
                const double h = sx[i] - sx[i - 1];
                slope[i] = slope[i - 1] + 0.5 * h * (kappa[i] + kappa[i - 1]);
                deflection[i] = deflection[i - 1] + 0.5 * h * (slope[i] + slope[i - 1]);
            }
            return deflection;
        };
        std::vector<double> curvature(sx.size(), 0.0), ownCurvature(sx.size(), 0.0);
        int noSection = 0;
        for (std::size_t i = 0; i < sx.size(); ++i) {
            const sim::HullGirderSection s = sim::hullGirderSection(structure, sx[i]);
            if (s.secondMoment > 0)
                curvature[i] = -tier0.stations[i].moment / (youngs * s.secondMoment);
            else
                ++noSection;
            std::size_t which = 0;
            for (std::size_t j = 0; j < chain.section.size(); ++j)
                if (chain.section[j].xFrom <= sx[i] && sx[i] < chain.section[j].xTo) which = j;
            if (tierOneSecond[which] > 0)
                ownCurvature[i] = -tier0.stations[i].moment / (youngs * tierOneSecond[which]);
        }
        const std::vector<double> beamDeflection = integrate(curvature);
        const std::vector<double> ownDeflection = integrate(ownCurvature);

        // **The Tier-1 deflection of each cut plane, from `u_x` and not from `u_z`, and
        // the difference between the two is most of what this measurement is worth.**
        //
        // A plane's mean `u_z` is the obvious reading and it is the wrong one here. The
        // load is applied to every node in a slab in proportion to its steel, which puts
        // the right resultant on each station and a deliberately unphysical distribution
        // inside it -- buoyancy on the deck as well as on the shell -- so the plating
        // deflects locally under it and the mean carries that as well as the girder's
        // bending. Measured: the second difference of the mean `u_z` is 3.8 times the
        // curvature the *stress* field carries at the same station, which is not
        // something bending can do.
        //
        // The bending is in `u_x`. A beam's axial displacement is `-z dw/dx`, so a
        // volume-blind least squares of `u_x` against `z` over the plane's own rows
        // gives the slope directly, and local plate motion contributes nothing to it --
        // a panel bulging between frames moves `u_z`, not the plane's `du_x/dz`. `w` is
        // its integral. The mean `u_z` is reported alongside, because the gap between
        // them *is* the local response and hiding it would be worse than publishing it.
        std::vector<double> planeX, planeW, planeHeave;
        {
            std::vector<double> slope;
            for (int i = 0; i <= options.waveSections; ++i) {
                const double x = from + (to - from) * i / options.waveSections;
                double n = 0, sz = 0, su = 0, szz = 0, szu = 0, heave = 0, heaveCount = 0;
                for (int b = 0; b < chain.assembly.boundary; ++b) {
                    const sim::reduction::BoundaryDof& d =
                        chain.assembly.boundaryPoint[static_cast<std::size_t>(b)];
                    if (std::abs(d.position.x - x) > 1e-6) continue;
                    const double u = state[static_cast<std::size_t>(b)];
                    if (d.axis == 2) {
                        heave += u;
                        heaveCount += 1.0;
                    }
                    if (d.axis != 0) continue;
                    n += 1.0;
                    sz += d.position.z;
                    su += u;
                    szz += d.position.z * d.position.z;
                    szu += d.position.z * u;
                }
                const double det = n * szz - sz * sz;
                if (!(n > 2) || !(det > 0) || !(heaveCount > 0)) continue;
                planeX.push_back(x);
                slope.push_back(-(n * szu - sz * su) / det);
                planeHeave.push_back(heave / heaveCount);
            }
            planeW.assign(planeX.size(), 0.0);
            for (std::size_t i = 1; i < planeX.size(); ++i)
                planeW[i] = planeW[i - 1] +
                            0.5 * (planeX[i] - planeX[i - 1]) * (slope[i] + slope[i - 1]);
        }
        // Both beams' deflection at those same planes, linearly interpolated.
        std::vector<double> beamAt, ownAt;
        for (double x : planeX) {
            std::size_t i = 1;
            while (i + 1 < sx.size() && sx[i] < x) ++i;
            const double t = (x - sx[i - 1]) / (sx[i] - sx[i - 1]);
            beamAt.push_back(beamDeflection[i - 1] + t * (beamDeflection[i] - beamDeflection[i - 1]));
            ownAt.push_back(ownDeflection[i - 1] + t * (ownDeflection[i] - ownDeflection[i - 1]));
        }
        const auto deTrim = [](const std::vector<double>& x, std::vector<double>& w) {
            const double n = static_cast<double>(x.size());
            if (n < 3) return;
            double sx0 = 0, sy = 0, sxx = 0, sxy = 0;
            for (std::size_t i = 0; i < x.size(); ++i) {
                sx0 += x[i];
                sy += w[i];
                sxx += x[i] * x[i];
                sxy += x[i] * w[i];
            }
            const double d = n * sxx - sx0 * sx0;
            if (!(std::abs(d) > 0)) return;
            const double b = (n * sxy - sx0 * sy) / d, a = (sy - b * sx0) / n;
            for (std::size_t i = 0; i < x.size(); ++i) w[i] -= a + b * x[i];
        };
        deTrim(planeX, planeW);
        deTrim(planeX, planeHeave);
        deTrim(planeX, beamAt);
        deTrim(planeX, ownAt);

        std::printf("\n  the deflected shape, heave and trim removed from every column:\n");
        std::printf("  %8s %12s %12s %12s %10s %10s %12s %12s\n", "x m", "Tier 1 mm", "beam, I0",
                    "beam, I1", "vs I0", "vs I1", "mean u_z", "M N m");
        double peak = 0, worstGap = 0, worstOwn = 0, worstHeave = 0, rms = 0, worstMiddle = 0;
        for (std::size_t i = 0; i < planeX.size(); ++i) {
            std::size_t j = 1;
            while (j + 1 < sx.size() && sx[j] < planeX[i]) ++j;
            std::printf("  %8.1f %12.4f %12.4f %12.4f %10.4f %10.4f %12.4f %12.4e\n", planeX[i],
                        1000.0 * planeW[i], 1000.0 * beamAt[i], 1000.0 * ownAt[i],
                        1000.0 * (planeW[i] - beamAt[i]), 1000.0 * (planeW[i] - ownAt[i]),
                        1000.0 * planeHeave[i], tier0.stations[j].moment);
            peak = std::max(peak, std::abs(planeW[i]));
            worstGap = std::max(worstGap, std::abs(planeW[i] - beamAt[i]));
            worstOwn = std::max(worstOwn, std::abs(planeW[i] - ownAt[i]));
            worstHeave = std::max(worstHeave, std::abs(planeHeave[i] - planeW[i]));
            if (std::abs(planeX[i]) < 0.3 * shipLength)
                worstMiddle = std::max(worstMiddle, std::abs(planeW[i] - beamAt[i]));
            rms += (planeW[i] - beamAt[i]) * (planeW[i] - beamAt[i]);
        }
        rms = std::sqrt(rms / std::max<std::size_t>(planeX.size(), 1));
        std::printf("  peak Tier-1 bending deflection %.3f mm; against Tier 0's own EI it is out"
                    " by %.3f mm over the middle two thirds and %.3f mm at worst, and putting the"
                    " *sections'* own EI into the same integration takes the latter to %.3f mm --"
                    " so about a quarter of the difference is EI(x) and the rest is shear and the"
                    " ends, which are not separated here\n", 1000.0 * peak, 1000.0 * worstMiddle,
                    1000.0 * worstGap, 1000.0 * worstOwn);
        std::printf("  the mean u_z of the same planes differs from the bending shape by up to"
                    " %.3f mm, which is the local response to a load spread over steel rather"
                    " than applied as a pressure -- it is not girder deflection and it is why"
                    " the column above is taken from u_x\n", 1000.0 * worstHeave);
        std::printf("  Tier 0 had no section at %d of %d stations; rms difference %.3f mm\n",
                    noSection, static_cast<int>(sx.size()), 1000.0 * rms);
        std::fflush(stdout);

        // --- And the stress, which is the half a beam cannot produce --------------------
        //
        // Section by section along the length. The chain's interface displacement is
        // exact -- `reduction.hpp` property 1 holds for a load applied inside, which is
        // the whole reason this is affordable -- so prescribing it on that section's own
        // mesh and solving with that section's own share of the load gives the **exact**
        // stress field of the whole-ship model there. Shear included, which is what
        // matters: shear lag is driven by `dM/dx`, and a section handed a constant
        // moment has none of it to show.
        //
        // Four things are compared and only the first is one a beam can also produce:
        //
        //   * `dsigma/dz` against `M / I`. This is the *local* form of the deflection
        //     comparison above and it is far less delicate, because it is not a double
        //     integral of anything.
        //   * the deck fibre, `M / Z`, against the worst stress really found on the deck.
        //   * that worst against the deck's own **mean**, which is shear lag and which
        //     Tier 0 cannot report at all: a beam says every point at one `z` carries
        //     one stress.
        //   * the residual of the best straight line through the field, which is
        //     identically zero for a beam and is everything the plane-sections
        //     assumption threw away.
        std::printf("\n  the stress along the length, from the whole-ship model, at each section's"
                    " mid-length:\n");
        std::printf("  %8s %11s %9s %9s %8s %9s %9s %7s %8s %8s\n", "x m", "M N m", "sig0 MPa",
                    "sig1 MPa", "mean", "worst/mn", "M/I MPa/m", "Tier 1", "kappa1/0", "resid %");
        double worstLag = 0, worstLagX = 0, worstResidual = 0;
        double gradientAmidships = 0, gradientEnds = 0;
        sim::section::BeamFit peakFit;
        sim::section::FibreStress peakDeck, peakKeel;
        double peakBeamDeck = 0, peakBeamKeel = 0, peakStation = 0, peakMoment = 0;
        double peakNeutralAxis = 0;
        int stressed = 0;
        for (std::size_t s = 0; s < chain.section.size(); ++s) {
            const sim::section::Section& piece = chain.section[s];
            const std::vector<double> interfaceField =
                sim::section::sectionDisplacement(chain, s, state);
            if (interfaceField.empty()) {
                std::printf("       ! the chain's own state does not expand onto section %zu\n", s);
                return 1;
            }
            sim::solidshell::HexMesh loaded = piece.mesh;
            for (std::uint32_t node : piece.interfaceNodes)
                for (int k = 0; k < 3; ++k)
                    loaded.pin(node, k, interfaceField[static_cast<std::size_t>(node) * 3 +
                                                       static_cast<std::size_t>(k)]);
            std::vector<double> exact;
            if (!sim::solidshell::solveStatic(loaded, material,
                                              sim::solidshell::Formulation::SolidShell,
                                              piece.attachment.stiffness,
                                              piece.attachment.constrained, meshLoad[s], exact,
                                              &problem)) {
                std::printf("       ! section %zu did not solve: %s\n", s, problem.c_str());
                return 1;
            }
            const double station = 0.5 * (piece.xFrom + piece.xTo);
            const std::vector<sim::section::StressSample> samples =
                sim::section::axialStress(piece, material, exact, station);
            const sim::HullGirderSection cut = sim::hullGirderSection(structure, station);
            if (samples.empty() || !(cut.secondMoment > 0)) {
                std::printf("  %8.2f   -- no cut to compare\n", station);
                continue;
            }
            const sim::section::BeamFit fit = sim::section::fitBeam(samples, cut.neutralAxis);
            const sim::section::FibreStress deck = sim::section::fibreStress(samples, true);
            const sim::section::FibreStress keel = sim::section::fibreStress(samples, false);
            std::size_t index = 0;
            for (std::size_t i = 0; i < sx.size(); ++i)
                if (std::abs(sx[i] - station) < std::abs(sx[index] - station)) index = i;
            const double moment = tier0.stations[index].moment;
            const double beamDeck = cut.modulusDeck > 0 ? moment / cut.modulusDeck : 0.0;
            const double beamGradient = moment / cut.secondMoment;
            const double lag = deck.mean != 0 ? deck.worst / deck.mean : 0.0;
            ++stressed;
            std::printf("  %8.2f %11.4e %9.3f %9.3f %8.3f %9.3f %9.3f %7.3f %8.3f %8.2f\n", station,
                        moment, beamDeck / 1e6, deck.worst / 1e6, deck.mean / 1e6, lag,
                        beamGradient / 1e6, fit.gradient / 1e6,
                        beamGradient != 0 ? fit.gradient / beamGradient : 0.0,
                        fit.peak != 0 ? 100.0 * fit.residualRms / std::abs(fit.peak) : 0.0);
            std::fflush(stdout);
            // **Only where the beam has an answer worth dividing by.** Towards the ends
            // the deck's mean stress passes through zero, so `worst / mean` runs off to
            // seventy and says nothing about shear lag; and the highest fibre there is a
            // stem or a forecastle rather than the strength deck, which is a place the
            // beam's own notion of a deck modulus has already stopped meaning anything.
            // The statistic is therefore taken over the sections carrying at least
            // **half** the peak moment, and the ends are reported rather than averaged
            // in. (This said "a quarter" against a `0.5 *` two lines below, which is
            // the wrong window to picture when reading the number it produces.)
            if (std::abs(moment) > 0.5 * std::abs(tier0.maxMoment)) {
                if (std::abs(lag) > std::abs(worstLag)) {
                    worstLag = lag;
                    worstLagX = station;
                }
                worstResidual = std::max(worstResidual,
                                         fit.peak != 0 ? fit.residualRms / std::abs(fit.peak) : 0.0);
            }
            if (std::abs(station) < 0.2 * shipLength && beamGradient != 0)
                gradientAmidships = fit.gradient / beamGradient;
            if (std::abs(station) > 0.3 * shipLength && beamGradient != 0)
                gradientEnds = std::max(gradientEnds, std::abs(fit.gradient / beamGradient));
            if (piece.xFrom <= tier0.maxMomentX && tier0.maxMomentX < piece.xTo) {
                peakFit = fit;
                peakDeck = deck;
                peakKeel = keel;
                peakBeamDeck = beamDeck;
                peakBeamKeel = cut.modulusKeel > 0 ? -moment / cut.modulusKeel : 0.0;
                peakStation = station;
                peakNeutralAxis = cut.neutralAxis;
                peakMoment = moment;
            }
        }
        if (!peakFit.ok) {
            std::printf("       ! the section carrying the peak moment produced no stress field\n");
            return 1;
        }
        std::printf("\n  at x = %.2f m, where the moment peaks at %.4e N m:\n", peakStation,
                    peakMoment);
        std::printf("  %-28s %12s %12s\n", "", "Tier 0", "Tier 1");
        std::printf("  %-28s %12.3f %12.3f   MPa\n", "deck fibre", peakBeamDeck / 1e6,
                    peakDeck.worst / 1e6);
        std::printf("  %-28s %12s %12.3f   MPa\n", "  its mean over the deck", "(one number)",
                    peakDeck.mean / 1e6);
        std::printf("  %-28s %12.3f %12.3f   MPa\n", "keel fibre", peakBeamKeel / 1e6,
                    peakKeel.worst / 1e6);
        std::printf("  %-28s %12.4f %12.4f   m\n", "neutral axis", peakNeutralAxis,
                    peakFit.neutralAxis);
        std::printf("  %-28s %12s %12.3f   MPa\n", "what a beam cannot carry", "0.000",
                    peakFit.residualRms / 1e6);
        std::printf("  %-28s %12s %12.3f   MPa at (%.2f, %.2f, %.2f)\n", "worst |sigma_xx|", "-",
                    peakFit.peak / 1e6, peakFit.peakAt.x, peakFit.peakAt.y, peakFit.peakAt.z);
        std::printf("  worst shear lag %.3f at x = %.1f; worst residual %.1f%% of the peak"
                    " stress; the local curvature is %.3f of the beam's amidships and up to %.3f"
                    " towards the ends\n", worstLag, worstLagX, 100.0 * worstResidual,
                    gradientAmidships, gradientEnds);
        std::fflush(stdout);
        // --- The success contract -------------------------------------------------------
        //
        // What is asserted is that the two tiers were asked *the same* question and that
        // the answer is not vacuous: she really hogs, the load really balanced, the
        // restraints really carry nothing, and the Tier-1 field really is not a beam. The
        // *size* of the disagreement is reported and not asserted, because it is the
        // finding.
        if (!tier0.hogging() || !(std::abs(tier0.maxMoment) > 1e7)) {
            std::printf("       ! peak moment %.3e N m: a crest amidships has to hog and it has to"
                        " hog by something, or nothing below is being driven\n", tier0.maxMoment);
            return 1;
        }
        if (!(restraint < 1e-3 * biggest)) {
            std::printf("       ! the six restraints carry %.3e N against %.3e applied. They are"
                        " statically determinate, so anything they carry is load that did not"
                        " balance and every moment below is wrong by it\n", restraint, biggest);
            return 1;
        }
        if (!(peak > 1e-4)) {
            std::printf("       ! the ship deflects %.3e m under this wave, which is not a"
                        " response\n", peak);
            return 1;
        }
        if (stressed < 3) {
            std::printf("       ! only %d sections produced a stress field to compare\n", stressed);
            return 1;
        }
        if (!(std::abs(peakDeck.worst) > 1e7) || !(std::abs(peakBeamDeck) > 1e7)) {
            std::printf("       ! deck stress %.3e Tier 0 and %.3e Tier 1: one of the two is not"
                        " carrying the moment at all\n", peakBeamDeck, peakDeck.worst);
            return 1;
        }
        if (!(peakDeck.worst * peakBeamDeck > 0) || !(peakKeel.worst * peakBeamKeel > 0)) {
            std::printf("       ! the two tiers disagree about the *sign* at the deck (%.3e /"
                        " %.3e) or the keel (%.3e / %.3e). Hogging tensions the deck and"
                        " compresses the keel, and a magnitude alone cannot tell those apart\n",
                        peakBeamDeck, peakDeck.worst, peakBeamKeel, peakKeel.worst);
            return 1;
        }
        // **The deck's mean is the beam's own answer and its worst is not**, and both
        // halves have to hold or the comparison is measuring one model twice: a mean
        // that did not agree would mean the moment is not being carried, and a worst
        // that agreed would mean the 3D field had come back a beam.
        if (!(std::abs(peakDeck.mean / peakBeamDeck - 1.0) < 0.15)) {
            // "differs from M/Z by", not "is % of M/Z": the value is `100*(r - 1)`,
            // and the guard only fires above 15, so "+18.3% of M/Z" read literally
            // says the moment is nearly absent when it is agreeing to 18%. On a path
            // that returns 1, that sends the reader after a missing load rather than
            // a distribution error.
            std::printf("       ! the deck's *mean* stress differs from M/Z by %+.1f%%. The mean"
                        " over a fibre is the one thing a beam does get right, so a disagreement"
                        " here is the moment not arriving rather than a distribution\n",
                        100.0 * (peakDeck.mean / peakBeamDeck - 1.0));
            return 1;
        }
        if (!(worstLag > 1.05)) {
            std::printf("       ! the worst deck stress is %.3f of the deck's mean. A reduced 3D"
                        " model that reproduces a beam exactly has not been asked a 3D question --"
                        " and shear lag needs dM/dx, so this is also what a constant-moment drive"
                        " would report\n", worstLag);
            return 1;
        }
        if (!(peakFit.residualRms > 0) || !(peakFit.samples > 100)) {
            std::printf("       ! the Tier-1 stress field is a beam to the last bit over %zu"
                        " samples. Either the recovery is not reaching the interior or this"
                        " comparison is measuring one model twice\n", peakFit.samples);
            return 1;
        }
        std::printf("\nok\n");
        return 0;
    }

    // --- What the mesher built -----------------------------------------------------

    sim::section::SectionParams params;
    params.xFrom = options.from;
    params.xTo = options.to;
    params.subdivision = options.subdivision;

    const double built = now();
    const sim::section::Section hold = sim::section::buildSection(structure, params);
    const double meshSeconds = now() - built;
    std::printf("=== section [%.2f, %.2f] m, subdivision %d, meshed in %.3f s ===\n", options.from,
                options.to, options.subdivision, meshSeconds);
    if (hold.empty()) {
        for (const std::string& problem : hold.problems) std::printf("  ! %s\n", problem.c_str());
        return 1;
    }
    std::printf("  %zu elements, %zu nodes, DOF half-bandwidth %zu\n", hold.elementCount(),
                hold.nodeCount(), hold.halfBandwidth);
    std::printf("  %zu interface nodes (%zu boundary DOF), %d surfaces, %d components,"
                " %d floating, %d spanning\n",
                hold.interfaceNodes.size(), 3 * hold.interfaceNodes.size(), hold.surfaces,
                hold.components, hold.floatingComponents, hold.spanningComponents);
    std::printf("  %.1f m^2 of mid-surface, %.0f kg of plate + %.0f kg of member = %.0f kg\n",
                hold.area, hold.plateMass, hold.memberMass, hold.mass());
    std::printf("  worst Jacobian %.3e, worst aspect %.2f, normal spread %.4f rad on %d elements"
                " (%.0f%% excess plate bending)\n",
                hold.worstJacobian, hold.worstAspect, hold.worstNormalSpread,
                hold.distortedElements, 100.0 * hold.spuriousStiffness);
    std::printf("  %d tapered elements, worst dt/t %.4f (%.0f%% excess plate bending)\n",
                hold.taperedElements, hold.worstTaper, 100.0 * hold.taperStiffness);
    std::printf("  free edge %.1f m, of which %.1f m sits on plating it is not welded to"
                " (worst gap %.4f m)\n",
                hold.freeEdgeLength, hold.junctionEdges, hold.worstJunctionGap);
    std::printf("  junction ties: %d nodes tied joining %.1f m of that edge; %d left on a cut"
                " plane, %d refused as a chain, %d outside a master face\n",
                hold.junctionTies, hold.tiedEdges, hold.junctionsOnInterface,
                hold.junctionsChained, hold.junctionsOutsideFace);
    std::printf("  worst tie: overshoot %.4f of a face, through-thickness weight %.4f\n",
                hold.worstJunctionOvershoot, hold.worstJunctionWeight);
    std::printf("  in-plane line ties: %d nodes on the cut planes, joining %.1f m aft + %.1f m"
                " forward + %.1f m needing both when a chain applies them; %d unreached,"
                " %d off the end of a line, %d over weight, %d chained\n",
                hold.planeTieNodes, hold.planeTiedEdgesAft, hold.planeTiedEdgesForward,
                hold.planeTiedEdgesBoth, hold.planeTiesUnreached, hold.planeTiesOutsideLine,
                hold.planeTiesThroughThickness, hold.planeTiesChained);
    std::printf("  worst line tie: overshoot %.4f of a segment, weight %.4f, slip %.3e m\n",
                hold.worstPlaneTieOvershoot, hold.worstPlaneTieWeight, hold.worstPlaneTieSlip);
    std::printf("  panels straddling a cut plane: %d; halo panels averaged into the nodal"
                " normals and thicknesses and then dropped: %d\n",
                hold.straddlingPanels, hold.haloPanels);
    std::printf("  members: %d attached, %d refused, %d missed; effective area attached %.5f"
                " + missed %.5f = %.5f m^2\n",
                hold.membersAttached, hold.membersRefused, hold.membersMissed,
                hold.attachedMemberArea, hold.missedMemberArea,
                hold.attachedMemberArea + hold.missedMemberArea);
    for (const std::string& problem : hold.problems) std::printf("  ! %s\n", problem.c_str());

    // --- Against Tier 0 --------------------------------------------------------------

    const auto report = [&](const sim::section::Section& section, const char* label) {
        sim::section::BeamLoad axial;
        axial.strain = 1e-6;
        axial.reference = girder.neutralAxis;
        const sim::section::BeamResponse stretched =
            sim::section::applyBeamLoad(section, material, axial);
        sim::section::BeamLoad bending;
        bending.curvature = 1e-6;
        bending.reference = girder.neutralAxis;
        const sim::section::BeamResponse bent =
            sim::section::applyBeamLoad(section, material, bending);
        if (!stretched.ok || !bent.ok) {
            std::printf("  %-16s refused: %s %s\n", label, stretched.problem.c_str(),
                        bent.problem.c_str());
            return;
        }
        const double area = stretched.axialStiffness / youngs;
        const double neutralAxis =
            stretched.bendingMoment / stretched.axialForce + girder.neutralAxis;
        const double second = bent.bendingStiffness / youngs;
        const double predicted = girder.area - section.missedMemberArea;
        std::printf("  %-16s A_eff %.5f m^2 (%+.3f%% of Tier 0, %+.3f%% of Tier 0 less the"
                    " members it could not attach)\n",
                    label, area, 100.0 * (area / girder.area - 1.0),
                    100.0 * (area / predicted - 1.0));
        std::printf("  %-16s z_na  %.5f m (Tier 0 %.5f)   I_eff %.5f m^4 (%+.3f%% of Tier 0)\n", "",
                    neutralAxis, girder.neutralAxis, second,
                    100.0 * (second / girder.secondMoment - 1.0));
        std::printf("  %-16s residual %.2e N, rigid-body restraint reaction %.2e N\n", "",
                    stretched.residual, stretched.restraintReaction);
    };

    std::printf("\n=== against Tier 0 ===\n");
    report(hold, "with members");
    sim::section::SectionParams bare = params;
    bare.members = false;
    const sim::section::Section plating = sim::section::buildSection(structure, bare);
    report(plating, "bare plating");
    double stiffenerArea = 0, plateArea = 0;
    for (const sim::SectionElement& element : sim::sectionElements(structure, middle))
        (element.stiffener ? stiffenerArea : plateArea) += element.area;
    std::printf("  Tier 0 says the stiffeners are %.5f m^2 of %.5f, %.1f%% -- which is what"
                " bare plating must be short by\n",
                stiffenerArea, girder.area, 100.0 * stiffenerArea / girder.area);

    // --- Resolution ------------------------------------------------------------------

    // --- What the junction tie costs and buys, at one element per panel -------------
    //
    // The two quantities a prescribed plane-sections field cannot see. `EA` and `EI`
    // are above and they move 0.19% and 0.12%; these move by factors.
    // What the junction tie is worth, kept for the success contract at the end of
    // main: `EA` cannot see whether the plating is joined (§2 of `section.hpp`) and
    // torsion and the lowest fixed-interface frequency can, so those two are what
    // the gate is allowed to assert on.
    double torsionCut = 0, torsionTied = 0, hertzCut = 0, hertzTied = 0;
    int componentsCut = 0, componentsTied = 0;
    {
        std::printf("\n=== the junction tie: cut against tied ===\n");
        std::printf("  %-6s %7s %6s %9s %9s %8s %12s %10s %10s\n", "", "band", "comps", "tied m",
                    "A_eff", "solve s", "GJ", "first Hz", "reduce s");
        for (int tie = 0; tie < 2; ++tie) {
            sim::section::SectionParams part = params;
            part.junctions = tie != 0;
            const sim::section::Section piece = sim::section::buildSection(structure, part);
            const double solving = now();
            sim::section::BeamLoad axial;
            axial.strain = 1e-6;
            const sim::section::BeamResponse stretched =
                sim::section::applyBeamLoad(piece, material, axial);
            const sim::section::TorsionResponse twisted =
                sim::section::applyTwist(piece, material, 1e-6, girder.neutralAxis);
            const double solveSeconds = now() - solving;
            const double reducing = now();
            const sim::reduction::Substructure substructure(piece.mesh, piece.material,
                                                            piece.interfaceNodes, piece.attachment);
            double hz = 0;
            if (substructure.ready()) {
                const sim::reduction::Eigenpairs modes = substructure.fixedInterfaceModes(1);
                if (!modes.value.empty())
                    hz = std::sqrt(std::max(0.0, modes.value[0])) / (2.0 * std::numbers::pi);
            }
            // The two orderings that lost, which is what "the node numbering is a
            // hundredfold, not a rounding" is a claim *about*. The chooser kept only
            // the winner, so the figure the documents quote for the alternative --
            // 1 382 on the untied hold -- could be produced by no invocation at all,
            // and so could be neither checked nor refuted.
            // **Not led by "cut"/"tied".** `check-figures.sh`'s `tierow` selects on
            // the first field and takes the first match, so a second line starting
            // with the same word would have been read as the data row and quietly
            // broken eight gated figures.
            std::printf("  orderings-%-4s x-fastest %zu, y-fastest %zu, RCM %zu;"
                        " best without RCM %zu\n",
                        tie ? "tied" : "cut", piece.candidateBandwidth[0],
                        piece.candidateBandwidth[1], piece.candidateBandwidth[2],
                        std::min(piece.candidateBandwidth[0], piece.candidateBandwidth[1]));
            std::printf("  %-6s %7zu %6d %9.1f %9.5f %8.2f %12.4e %10.4f %10.2f\n",
                        tie ? "tied" : "cut", piece.halfBandwidth, piece.components,
                        piece.tiedEdges, stretched.axialStiffness / youngs, solveSeconds,
                        twisted.torsionalStiffness, hz, now() - reducing);
            std::fflush(stdout);
            (tie ? torsionTied : torsionCut) = twisted.torsionalStiffness;
            (tie ? hertzTied : hertzCut) = hz;
            (tie ? componentsTied : componentsCut) = piece.components;
        }
    }

    // The resolution sweep runs **untied**. It is a study of the mesher's own
    // convergence, and the tie's cost is a band: 146 against 1 520 at subdivision 1,
    // which is 0.16 s of banded factorisation against 5.34, and 278 against 3 188 at
    // subdivision 2, which is 1.1 s against 149. Sweeping tied would take an hour and
    // measure the same four columns to within 0.2%.
    if (options.sweep >= 1) {
        std::printf("\n=== resolution: refining the reduced answer (junctions untied) ===\n");
        std::printf("  %4s %9s %7s %11s %10s %12s %13s %9s\n", "sub", "elements", "band", "A_eff",
                    "z_na", "I_eff", "GJ", "solve s");
        for (int subdivision = 1; subdivision <= options.sweep; ++subdivision) {
            sim::section::SectionParams refined = params;
            refined.subdivision = subdivision;
            refined.junctions = false;
            const sim::section::Section section = sim::section::buildSection(structure, refined);
            const double start = now();
            sim::section::BeamLoad axial;
            axial.strain = 1e-6;
            axial.reference = girder.neutralAxis;
            const sim::section::BeamResponse stretched =
                sim::section::applyBeamLoad(section, material, axial);
            sim::section::BeamLoad bending;
            bending.curvature = 1e-6;
            bending.reference = girder.neutralAxis;
            const sim::section::BeamResponse bent =
                sim::section::applyBeamLoad(section, material, bending);
            const sim::section::TorsionResponse twisted =
                sim::section::applyTwist(section, material, 1e-6, girder.neutralAxis);
            std::printf("  %4d %9zu %7zu %11.5f %10.5f %12.5f %13.4e %9.2f\n", subdivision,
                        section.elementCount(), section.halfBandwidth,
                        stretched.axialStiffness / youngs,
                        stretched.bendingMoment / stretched.axialForce + girder.neutralAxis,
                        bent.bendingStiffness / youngs, twisted.torsionalStiffness, now() - start);
            std::fflush(stdout);
        }
    }

    // --- A ship: the same length as a chain of sections ------------------------------
    //
    // The end-to-end claim of the tier, at ship scale: cut a length into N pieces,
    // reduce each once, assemble, and get the model the same length in one piece
    // gives. The reference here owes nothing to the assembly -- it is
    // `applyBeamLoad` and `applyTwist` on the monolithic section, through
    // `solidshell::solveStatic`.
    //
    // Everything is untied unless asked otherwise, because tied at ship scale is a
    // band of 1 520 and a 5.3 s banded factorisation per section, and the two
    // questions -- does the assembly reproduce the monolith, and what do the cut
    // planes cost the ties -- are answered separately below.
    if (options.chain > 0) {
        std::printf("\n=== a chain of %d sections against the same length in one piece ===\n",
                    options.chain);
        for (int tie = 0; tie < 2; ++tie) {
            sim::section::ChainParams chainParams;
            chainParams.section = params;
            chainParams.section.junctions = tie != 0;
            chainParams.section.interfaceTies = options.interfaceTies;
            chainParams.reduce.modes = 0;
            chainParams.reduce.cutoffFrequency = 0;
            chainParams.matchTolerance = options.match;
            for (int i = 0; i <= options.chain; ++i)
                chainParams.station.push_back(options.from + (options.to - options.from) * i /
                                                                 options.chain);

            const double start = now();
            const sim::section::Chain chain = sim::section::buildChain(structure, chainParams);
            std::printf("  %-5s built in %.2f s (mesh %.2f, reduce %.2f, assemble %.2f):"
                        " %d pieces, %d assembled boundary DOF, size %d\n",
                        tie ? "tied" : "cut", now() - start, chain.meshSeconds,
                        chain.reduceSeconds, chain.assembleSeconds, chain.components,
                        chain.assembly.boundary, chain.assembly.size());
            std::size_t unmatched = 0;
            for (std::size_t u : chain.unmatched) unmatched += u;
            std::printf("  %-5s interior planes: %zu shared DOF each, %zu unmatched, worst gap"
                        " %.3e m; ties %.1f m of %.1f m of junction edge\n",
                        "", chain.shared.empty() ? 0 : chain.shared.front(), unmatched,
                        chain.worstGap, chain.tiedEdges, chain.junctionEdges);
            std::printf("  %-5s in-plane line ties: %d nodes over %d interior planes, %d planes the"
                        " two sides disagreed about (worst %.3e)\n",
                        "", chain.planeTieNodes, options.chain - 1, chain.planeTiesDisagreeing,
                        chain.worstPlaneTieDisagreement);
            for (const std::string& problem : chain.problems)
                std::printf("      ! %s\n", problem.c_str());
            std::fflush(stdout);
            if (!chain.ready()) continue;

            sim::section::SectionParams monoParams = params;
            monoParams.junctions = tie != 0;
            const sim::section::Section mono = sim::section::buildSection(structure, monoParams);
            std::printf("  %-5s one piece: %d components, ties %.1f m of %.1f m\n", "",
                        mono.components, mono.tiedEdges, mono.junctionEdges);

            sim::section::BeamLoad axial;
            axial.strain = 1e-6;
            axial.reference = girder.neutralAxis;
            sim::section::BeamLoad bending;
            bending.curvature = 1e-6;
            bending.reference = girder.neutralAxis;

            double at = now();
            const sim::section::BeamResponse chainAxial =
                sim::section::applyBeamLoad(chain, axial);
            const double chainSolve = now() - at;
            at = now();
            const sim::section::BeamResponse monoAxial =
                sim::section::applyBeamLoad(mono, material, axial);
            const double monoSolve = now() - at;
            const sim::section::BeamResponse chainBend = sim::section::applyBeamLoad(chain, bending);
            const sim::section::BeamResponse monoBend =
                sim::section::applyBeamLoad(mono, material, bending);
            const sim::section::TorsionResponse chainTwist =
                sim::section::applyTwist(chain, 1e-6, girder.neutralAxis);
            const sim::section::TorsionResponse monoTwist =
                sim::section::applyTwist(mono, material, 1e-6, girder.neutralAxis);
            if (!chainAxial.ok || !monoAxial.ok) {
                std::printf("      ! %s / %s\n", chainAxial.problem.c_str(),
                            monoAxial.problem.c_str());
                continue;
            }
            std::printf("  %-5s %-8s %14s %14s %10s\n", "", "", "chain", "one piece", "relative");
            const auto line = [&](const char* what, double a, double b) {
                std::printf("  %-5s %-8s %14.6e %14.6e %+10.3e\n", "", what, a, b,
                            b != 0 ? a / b - 1.0 : 0.0);
            };
            line("EA", chainAxial.axialStiffness, monoAxial.axialStiffness);
            line("z_na", chainAxial.bendingMoment / chainAxial.axialForce + girder.neutralAxis,
                 monoAxial.bendingMoment / monoAxial.axialForce + girder.neutralAxis);
            line("EI", chainBend.bendingStiffness, monoBend.bendingStiffness);
            line("GJ", chainTwist.torsionalStiffness, monoTwist.torsionalStiffness);
            std::printf("  %-5s dense assembled solve %.2f s against a banded monolithic %.2f s;"
                        " chain restraint reaction %.2e N, residual %.2e N\n",
                        "", chainSolve, monoSolve, chainAxial.restraintReaction,
                        chainAxial.residual);
            // One anchored line per variant, because `scripts/check-figures.sh` reads
            // these and a checker that misparses is worse than no checker: every figure
            // above is spread over a block whose rows carry no label of their own, and
            // the first version of that script picked the intact ship's GM out of a
            // banner because it matched on a word rather than on a line.
            std::printf("chain summary: N=%d junctions=%d lines=%d tiedEdges=%.1f"
                        " onePieceTiedEdges=%.1f planeTieNodes=%d disagree=%d GJrel=%+.4e"
                        " EArel=%+.4e\n",
                        options.chain, tie, options.interfaceTies ? 1 : 0, chain.tiedEdges,
                        mono.tiedEdges, chain.planeTieNodes, chain.planeTiesDisagreeing,
                        monoTwist.torsionalStiffness != 0
                            ? chainTwist.torsionalStiffness / monoTwist.torsionalStiffness - 1.0
                            : 0.0,
                        monoAxial.axialStiffness != 0
                            ? chainAxial.axialStiffness / monoAxial.axialStiffness - 1.0
                            : 0.0);
            std::fflush(stdout);
        }
    }

    // --- The reduction ----------------------------------------------------------------

    if (options.reduce) {
        std::printf("\n=== Craig-Bampton ===\n");
        // Which piece of the section owns the softest fixed-interface mode. Untied,
        // the answer is the decks and it is the *same number* with and without the
        // shell they should be welded to -- the junctions' cost stated as a
        // frequency. Tied, the whole section is stiffer than either piece, which is
        // what a joined structure does. See `section.hpp` §2 and §5.
        struct Case {
            const char* label;
            bool shell, deck, bulkhead, junctions;
        };
        const Case cases[] = {{"shell only", true, false, false, false},
                              {"decks only", false, true, false, false},
                              {"bulkheads only", false, false, true, false},
                              {"whole, untied", true, true, true, false},
                              {"whole, tied", true, true, true, true}};
        for (const Case& one : cases) {
            sim::section::SectionParams part = params;
            part.shell = one.shell;
            part.deck = one.deck;
            part.bulkhead = one.bulkhead;
            part.junctions = one.junctions;
            const sim::section::Section piece = sim::section::buildSection(structure, part);
            if (piece.empty()) continue;
            const double assembled = now();
            const sim::reduction::Substructure substructure(piece.mesh, piece.material,
                                                            piece.interfaceNodes, piece.attachment);
            const double assembleSeconds = now() - assembled;
            if (!substructure.ready()) {
                std::printf("  %-15s substructure refused\n", one.label);
                for (const std::string& problem : substructure.problems())
                    std::printf("      ! %s\n", problem.c_str());
                continue;
            }
            sim::reduction::ReduceParams reduceParams;
            reduceParams.modes = one.shell && one.deck && one.bulkhead ? options.modes : 0;
            reduceParams.verifyModes = reduceParams.modes != 0;
            const double reducing = now();
            const sim::reduction::Reduction reduced =
                sim::reduction::craigBampton(substructure, reduceParams);
            const double reduceSeconds = now() - reducing;
            std::printf("  %-15s %zu boundary + %zu interior DOF, band %zu, mass %.0f kg"
                        " (attached %.0f), assemble %.2f s\n",
                        one.label, substructure.boundaryCount(), substructure.interiorCount(),
                        substructure.halfBandwidth(), substructure.totalMass(),
                        substructure.attachedMass(), assembleSeconds);
            std::printf("  %-15s first fixed-interface mode %.4f rad/s = %.4f Hz;"
                        " %d modes kept in %.2f s\n",
                        "", reduced.firstFixedFrequency,
                        reduced.firstFixedFrequency / (2.0 * std::numbers::pi), reduced.modes,
                        reduceSeconds);
            // **The subspace iteration says it did not converge, so the frequency
            // above is not evidence on its own.** `eigenvaluesBelow` counts rather
            // than converges -- the inertia of an LDL^T factorisation, Sylvester's
            // law -- so bracketing the reported value with it is a different
            // instrument answering the same question, which is the only kind of
            // agreement worth having here.
            const double omega = reduced.firstFixedFrequency;
            bool exactBelow = false, exactAbove = false;
            const int below = substructure.eigenvaluesBelow(0.99 * 0.99 * omega * omega, &exactBelow);
            const int above = substructure.eigenvaluesBelow(1.01 * 1.01 * omega * omega, &exactAbove);
            std::printf("  %-15s inertia count brackets it: %d modes below %.4f rad/s, %d below"
                        " %.4f (exact %d/%d)\n",
                        "", below, 0.99 * omega, above, 1.01 * omega, exactBelow, exactAbove);
            for (const std::string& problem : reduced.problems)
                std::printf("      ! %s\n", problem.c_str());
            std::fflush(stdout);
        }
    }

    // --- The success contract ------------------------------------------------------
    //
    // **This program had none until it was put in the gate**, and it published a
    // ship-scale frequency that nothing ever re-ran. What it asserts is chosen the
    // same way `section.hpp` §2 chooses what a test may assert on: `EA` and `EI` are
    // exact on a section whose plating is joined to nothing, so they are checked
    // only loosely and against Tier 0, while everything about *joining* is checked
    // on torsion, on the lowest fixed-interface frequency and on the component
    // count -- the three quantities a prescribed plane-sections field cannot see.
    if (hold.invertedElements != 0 || !(hold.worstJacobian > 0)) {
        std::printf("       ! the section has %d inverted elements and a worst Gauss determinant"
                    " of %.3e\n", hold.invertedElements, hold.worstJacobian);
        return 1;
    }
    if (hold.elementCount() < 1000 || !(hold.area > 100.0)) {
        std::printf("       ! %zu elements over %.1f m2 is not a ship section, so nothing below"
                    " means anything\n", hold.elementCount(), hold.area);
        return 1;
    }
    if (!(torsionCut > 0) || !(torsionTied > 1.3 * torsionCut)) {
        std::printf("       ! tying the junctions moved GJ from %.4e to %.4e. A tie that does not"
                    " close the cell has not joined the plating, and EA would not have said so\n",
                    torsionCut, torsionTied);
        return 1;
    }
    if (!(hertzCut > 0) || !(hertzTied > 1.3 * hertzCut)) {
        std::printf("       ! the first fixed-interface mode went %.4f -> %.4f Hz. Untied it is"
                    " the decks' own frequency; a tie that does not raise it has joined nothing\n",
                    hertzCut, hertzTied);
        return 1;
    }
    if (componentsTied != 1 || componentsCut <= 1) {
        std::printf("       ! components went %d untied -> %d tied, and the pair has to be"
                    " many -> one or the tie is being credited with a mesh that was already"
                    " joined\n", componentsCut, componentsTied);
        return 1;
    }
    std::printf("\nok\n");
    return 0;
}
