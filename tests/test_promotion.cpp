// SPDX-License-Identifier: MIT
//
// Adaptive zone promotion: `engine/sim/promotion.{hpp,cpp}`.
//
// Three separate things are under test and they fail in different ways, so they
// are tested separately:
//
//   1. **The criterion.** Its interesting failure is not missing a hot spot -- it
//      is promoting *everything*, because the cost is linear in the number of
//      zones and a criterion that always fires is unaffordable rather than merely
//      wrong. So the tests here are mostly negative: a ship at rest promotes
//      nothing, a uniformly loaded ship promotes nothing, a broad peak promotes
//      one zone and not eight. Each of those carries its own vacuity guard,
//      because "nothing promoted" is also what a criterion that never fires
//      reports.
//   2. **The anti-chatter mechanism.** Tested against its own negative control:
//      the same signal with the mechanism switched off has to chatter, or the
//      test is measuring a signal that was never going to chatter anyway.
//   3. **The pre-load and the reaction back.** Both against closed forms -- the
//      stress a patch is handed, the energy it stores, the capacity it spends --
//      rather than against the previous run's output.
#include "harness.hpp"

#include "engine/sim/promotion.hpp"
#include "engine/sim/solid_shell.hpp"
#include "game/prototype/ferry.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <string>
#include <vector>

using namespace sim;
using testing::expectNear;
using testing::expectTrue;

namespace {

// --- Shared fixtures -----------------------------------------------------------

// The reference ferry's structure, built once: 8900 panels at 3.7 ms a time, and
// six tests want it.
const StructuralMesh& ferryStructure() {
    static const StructuralMesh mesh = [] {
        Ship ferry = game::buildFerry();
        ferry.initialise(0.0);
        return makeStructuralMesh(ferry.hull, ferryScantlings());
    }();
    return mesh;
}

Ship ferryAfloat() {
    Ship ferry = game::buildFerry();
    ferry.initialise(0.0);
    return ferry;
}

// A wave of the ship's length with the crest amidships, found by scanning a full
// period rather than assumed -- the same construction `test_girder.cpp` uses.
struct Crest {
    WaveField field;
    Sea sea;
};

Crest crestAmidships(const Ship& ship, double amplitude) {
    const double length = ship.hullHi.x - ship.hullLo.x;
    const double omega = std::sqrt(kGravity * 2.0 * kPi / length);
    Crest out{WaveField::regular(amplitude, omega, 0.0), Sea{}};
    out.sea.waves = &out.field;
    const double period = 2.0 * kPi / omega;
    double best = -1e30;
    for (int i = 0; i < 720; ++i) {
        const double t = period * i / 720.0;
        const double eta = out.field.elevation(ship.state.position.x, 0.0, t);
        if (eta > best) {
            best = eta;
            out.sea.time = t;
        }
    }
    return out;
}

// A `TierZero` with a chosen buckling profile and nothing else, so the criterion
// can be handed a shape rather than a ship. Everything else is zero, which means
// the yield and collapse triggers cannot fire and the buckling one is on its own.
promotion::TierZero syntheticBuckling(const std::vector<double>& x,
                                      const std::vector<double>& utilisation) {
    promotion::TierZero tier;
    tier.yieldStrength = ah36Steel().yieldStrength;
    for (std::size_t i = 0; i < x.size(); ++i) {
        GirderBuckling b;
        b.x = x[i];
        b.utilisation = utilisation[i];
        b.deckInCompression = false;  // the keel, which is unambiguous on the centreline
        b.compressiveStress = utilisation[i] * 1e8;
        tier.buckling.push_back(b);
        if (utilisation[i] > tier.buckleUtilisation) {
            tier.buckleUtilisation = utilisation[i];
            tier.buckleX = x[i];
        }
    }
    return tier;
}

std::vector<double> stationsAlong(double from, double to, int count) {
    std::vector<double> out;
    for (int i = 0; i < count; ++i)
        out.push_back(from + (to - from) * i / (count - 1));
    return out;
}

// A flat strip as a `StructuralMesh`, the same fixture `test_zone.cpp` uses --
// built directly rather than through `makeStructuralMesh` because a rectangle's
// area is `lengthX * spanY` and nothing else.
StructuralMesh flatStrip(double lengthX, double spanY, double thickness, int nx, int ny,
                         const StiffenerProfile* profile = nullptr) {
    StructuralMesh mesh;
    mesh.materials.push_back(ah36Steel());
    mesh.frameSpacing = lengthX;
    for (int i = 0; i < nx; ++i)
        for (int j = 0; j < ny; ++j) {
            PlatePanel p;
            const double x0 = -0.5 * lengthX + lengthX * i / nx;
            const double x1 = -0.5 * lengthX + lengthX * (i + 1) / nx;
            const double y0 = -0.5 * spanY + spanY * j / ny;
            const double y1 = -0.5 * spanY + spanY * (j + 1) / ny;
            p.corner[0] = {x0, y0, 0};
            p.corner[1] = {x1, y0, 0};
            p.corner[2] = {x1, y1, 0};
            p.corner[3] = {x0, y1, 0};
            p.thickness = thickness;
            p.material = 0;
            p.role = PanelRole::Shell;
            mesh.panels.push_back(p);
        }
    if (profile != nullptr) {
        // One longitudinal down the seam at y = 0, its web rising *inward* against
        // the patch's +z outward normal -- the sign a shell longitudinal has, and
        // the one every tie in `constraint.cpp` has to survive.
        StructuralMember member;
        member.a = {-0.5 * lengthX, 0.0, 0.0};
        member.b = {0.5 * lengthX, 0.0, 0.0};
        member.rise = {0, 0, -1};
        member.profile = *profile;
        member.attachedPlateThickness = thickness;
        member.role = MemberRole::Longitudinal;
        mesh.members.push_back(member);
    }
    return mesh;
}

zone::MeshParams flatParams(int subdivision) {
    zone::MeshParams params;
    params.radius = 1e3;
    params.subdivision = subdivision;
    params.outward = {0, 0, 1};
    return params;
}

bool anyPointYielded(const zone::Solver& solver) {
    for (const solidshell::ElementPlasticState& e : solver.elementState())
        for (int gp = 0; gp < solidshell::kGauss; ++gp)
            if (e.point[gp].equivalentPlasticStrain > 0) return true;
    return false;
}

// --- 1. Nothing promotes on a ship at rest ------------------------------------
//
// The failure this guards against is a criterion that always fires, which is
// unaffordable rather than merely noisy. The failure it *could* hide is a
// criterion that never fires, so both vacuity guards are here: the ship is really
// being evaluated (her still-water utilisations are non-zero), and the very same
// criterion on the very same still water does promote when a contact patch
// arrives.

void testNothingPromotesOnAShipAtRestInStillWater() {
    std::printf("\n   a ship at rest in still water\n");
    const Ship ferry = ferryAfloat();
    const StructuralMesh& structure = ferryStructure();
    const Scantlings scantlings = ferryScantlings();
    const promotion::TierZero still =
        promotion::tierZero(ferry, Sea(0.0), structure, scantlings);

    std::printf("     still water: yield %.4f, buckling %.4f, collapse %.4f;"
                " Tier-0 cost %.0f ms\n", still.yieldUtilisation, still.buckleUtilisation,
                still.collapseUtilisation, still.seconds * 1e3);

    // Vacuity guard one: she is being evaluated. A ferry floating on her own
    // weight distribution carries a real still-water bending moment, and a Tier-0
    // answer of exactly zero would mean the calculation did not happen.
    expectTrue("the still-water calculation produced stations",
               still.stress.size() > 20 && still.buckling.size() > 20 &&
                   still.strength.size() > 20);
    expectTrue("and she carries a real still-water bending moment",
               still.yieldUtilisation > 0.01 && still.buckleUtilisation > 0.01);
    expectTrue("but nowhere near anything", still.buckleUtilisation < 0.5);

    promotion::Criterion criterion;
    criterion.mesh.radius = 3.0;
    promotion::Promoter promoter(criterion);
    for (int i = 0; i < 6; ++i) {
        const promotion::Review review = promoter.review(structure, still);
        expectTrue("nothing is even a candidate in still water", review.considered.empty());
    }
    testing::expectEqual("and nothing has been promoted after six reviews",
                         promoter.promotions(), 0);
    testing::expectEqual("with no zones active", static_cast<long long>(promoter.active().size()),
                         0);

    // Vacuity guard two, and the important one: the same criterion, the same
    // still water, and a bow against her side. If this did not promote, the test
    // above would be measuring a criterion that cannot fire at all.
    promotion::ContactPatch bow;
    bow.centre = {0.0, -9.9, 8.0};
    bow.radius = 1.5;
    bow.force = 20.0e6;
    promotion::Promoter struck(criterion);
    for (int i = 0; i < 4; ++i) struck.review(structure, still, {bow});
    std::printf("     the same still water with a 20 MN bow against her side promotes %d zone(s)\n",
                static_cast<int>(struck.active().size()));
    testing::expectEqual("a contact in still water does promote", struck.promotions(), 1);

    // And a gentle one does not. 200 kN over the same footprint is a tug pushing,
    // and the plating hinges at 402 kPa -- which is the closed form the criterion
    // uses, not a tuned number.
    promotion::ContactPatch tug = bow;
    tug.force = 0.2e6;
    promotion::Promoter berthing(criterion);
    for (int i = 0; i < 4; ++i) berthing.review(structure, still, {tug});
    testing::expectEqual("a tug pushing alongside does not", berthing.promotions(), 0);

    // The threshold is the bay's own plastic collapse pressure, so it is
    // *bracketed* rather than sampled once: a force 20% under it must not promote
    // and one 20% over must. A single far-away "no" passes just as well against a
    // criterion using the wrong power of the contact radius.
    int nearest = 0;
    double closest = 1e30;
    for (std::size_t i = 0; i < structure.panels.size(); ++i) {
        const double d = length(structure.panels[i].centroid() - bow.centre);
        if (d < closest) {
            closest = d;
            nearest = static_cast<int>(i);
        }
    }
    const PlatePanel& bay = structure.panels[static_cast<std::size_t>(nearest)];
    const double hinge = promotion::platingCollapsePressure(
        ah36Steel().yieldStrength, bay.thickness, promotion::panelSpan(bay));
    const double footprint = kPi * bow.radius * bow.radius;
    std::printf("     the struck bay is %.0f mm on a %.3f m span and hinges at %.0f kPa over a"
                " %.2f m2 footprint\n", bay.thickness * 1e3, promotion::panelSpan(bay),
                hinge / 1e3, footprint);
    for (double factor : {0.8, 1.2}) {
        promotion::ContactPatch load = bow;
        load.force = factor * hinge * footprint;
        promotion::Promoter promoter(criterion);
        for (int i = 0; i < 4; ++i) promoter.review(structure, still, {load});
        std::printf("     %.0f%% of the hinging force (%.2f MN): %d promotion(s)\n", 100.0 * factor,
                    load.force / 1e6, promoter.promotions());
        testing::expectEqual(factor < 1.0 ? "just under the hinging pressure nothing promotes"
                                          : "just over it one zone does",
                             promoter.promotions(), factor < 1.0 ? 0 : 1);
    }
}

// --- 2. A uniformly loaded ship must not promote everywhere --------------------

void testAUniformlyLoadedShipDoesNotPromoteEverywhere() {
    std::printf("\n   0.9 of capacity everywhere is a well designed ship, not forty zones\n");
    const StructuralMesh& structure = ferryStructure();
    const std::vector<double> x = stationsAlong(-45.0, 45.0, 31);

    std::vector<double> flat(x.size(), 0.90);
    const promotion::TierZero uniform = syntheticBuckling(x, flat);

    promotion::Criterion criterion;
    criterion.mesh.radius = 3.0;
    expectTrue("the flat profile really is past the absolute threshold",
               0.90 >= criterion.bucklePromote);

    const std::vector<promotion::Candidate> none =
        promotion::candidates(structure, uniform, {}, criterion);
    testing::expectEqual("a flat profile at 0.9 produces no candidates at all",
                         static_cast<long long>(none.size()), 0);

    // The negative control. With the local-excess test switched off, every one of
    // those stations qualifies -- which is what proves the result above is the
    // guard working rather than the threshold never being reached.
    promotion::Criterion noGuard = criterion;
    noGuard.localExcess = 0.0;
    const std::vector<promotion::Candidate> all =
        promotion::candidates(structure, uniform, {}, noGuard);
    std::printf("     flat at 0.90: %zu candidates with the median guard, %zu without\n",
                none.size(), all.size());
    expectTrue("without the guard the same profile floods the criterion", all.size() >= 10);

    // A real hull girder's profile is a peak, not a plateau, and a peak has to be
    // found. Same threshold, same guard.
    std::vector<double> peaked;
    for (double station : x)
        peaked.push_back(0.95 * std::exp(-(station * station) / (2.0 * 18.0 * 18.0)));
    const promotion::TierZero hump = syntheticBuckling(x, peaked);
    const std::vector<promotion::Candidate> found =
        promotion::candidates(structure, hump, {}, criterion);
    expectTrue("a peak at 0.95 over a 0.3 background is found", !found.empty());
    std::printf("     a Gaussian peak at 0.95: %zu stations qualify\n", found.size());

    // But the peak is broad -- a hull girder's always is -- so several adjacent
    // stations qualify, and two zones a station apart would pay twice for the same
    // plating. The expectation is *not* that this collapses to one zone: a peak
    // that clears the threshold over twenty metres of ship needs covering over
    // twenty metres of ship. The expectation is that no two zones overlap, and the
    // count that follows from that is a derived bound rather than an eyeballed
    // number.
    expectTrue("and the peak is broad enough for this to be a real question",
               found.size() >= 3);
    promotion::Promoter promoter(criterion);
    for (int i = 0; i < 4; ++i) promoter.review(structure, hump);
    double qualifyingSpan = 0;
    for (const promotion::Candidate& a : found)
        for (const promotion::Candidate& b : found)
            qualifyingSpan = std::max(qualifyingSpan, std::abs(a.x - b.x));
    const double pitch = criterion.separation * criterion.mesh.radius;
    const std::size_t bound = static_cast<std::size_t>(qualifyingSpan / pitch) + 1;
    std::printf("     %zu qualifying stations over %.0f m become %zu zone(s) at a %.0f m"
                " minimum pitch; the geometry allows at most %zu\n", found.size(), qualifyingSpan,
                promoter.active().size(), pitch, bound);
    expectTrue("separation thins the qualifying stations", promoter.active().size() < found.size());
    expectTrue("to no more than the stretch of ship divided by the pitch",
               promoter.active().size() <= bound);
    // And no two of them overlap, which is the property the count is derived from.
    for (std::size_t i = 0; i < promoter.active().size(); ++i)
        for (std::size_t j = i + 1; j < promoter.active().size(); ++j)
            expectTrue("no two zones are within the separation distance",
                       length(promoter.active()[i].impact - promoter.active()[j].impact) >= pitch);
    // The strongest claim on a zone wins, so the first one built is the peak.
    expectTrue("and the first zone built is the peak, not the first station that cleared",
               std::abs(promoter.active().front().x) < 6.0);
}

// --- 2b. The background statistic, the threshold, and the fibre -----------------
//
// Every check here was added because a mutant survived without it. A flat profile
// cannot tell a median from a mean from a minimum -- they are all the same number
// on it -- so the guard needs a *skewed* profile before it says anything, and the
// site a zone lands on needs a case where the two fibres are in different places.

void testTheBackgroundIsAMedianAndTheThresholdIsWhereItSays() {
    std::printf("\n   the background statistic, and where the threshold is\n");
    const StructuralMesh& structure = ferryStructure();
    const std::vector<double> x = stationsAlong(-45.0, 45.0, 31);
    promotion::Criterion criterion;
    criterion.mesh.radius = 3.0;

    // A ship uniformly worked to 0.88 with fine, lightly loaded ends. The three
    // statistics differ by a lot here, and only the median says "flat":
    // median 0.88, mean 0.72, minimum 0.05.
    std::vector<double> worked;
    for (std::size_t i = 0; i < x.size(); ++i) worked.push_back(i < 6 ? 0.05 : 0.88);
    double mean = 0;
    for (double u : worked) mean += u / static_cast<double>(worked.size());
    std::printf("     a hull worked to 0.88 amidships with light ends: median 0.88, mean %.3f,"
                " minimum %.3f\n", mean, *std::min_element(worked.begin(), worked.end()));
    expectTrue("the three statistics genuinely disagree on this profile", mean < 0.88 - 0.10);
    testing::expectEqual(
        "and a median background still promotes nothing",
        static_cast<long long>(
            promotion::candidates(structure, syntheticBuckling(x, worked), {}, criterion).size()),
        0);

    // A ship *past* capacity everywhere, with one station a little worse. The
    // excess is a difference and not a ratio, and here the two disagree: 2.15
    // against a 2.00 background is +0.15 of utilisation and only +7.5%.
    std::vector<double> overloaded(x.size(), 2.00);
    overloaded[x.size() / 2] = 2.15;
    const std::vector<promotion::Candidate> overloadedCandidates =
        promotion::candidates(structure, syntheticBuckling(x, overloaded), {}, criterion);
    std::printf("     2.15 against a 2.00 background is +0.15 of utilisation and +7.5%%:"
                " %zu candidate(s)\n", overloadedCandidates.size());
    testing::expectEqual("a difference of 0.15 clears a localExcess of 0.10",
                         static_cast<long long>(overloadedCandidates.size()), 1);

    // The threshold is `>=`, so a station sitting exactly on it qualifies. A ship
    // whose worst panel is exactly at its critical stress is a real case and the
    // interesting side of the boundary.
    std::vector<double> onTheLine(x.size(), 0.20);
    onTheLine[x.size() / 2] = criterion.bucklePromote;
    testing::expectEqual(
        "a station exactly at the promote threshold qualifies",
        static_cast<long long>(
            promotion::candidates(structure, syntheticBuckling(x, onTheLine), {}, criterion).size()),
        1);
    std::vector<double> justUnder = onTheLine;
    justUnder[x.size() / 2] = criterion.bucklePromote * (1.0 - 1e-9);
    testing::expectEqual(
        "and one a billionth below it does not",
        static_cast<long long>(
            promotion::candidates(structure, syntheticBuckling(x, justUnder), {}, criterion).size()),
        0);

    // **Candidates are one per panel, whatever the trigger, and the strongest claim
    // on a panel is the one that survives.** Two triggers are aimed at the same
    // station and the same fibre -- yield at 0.99 of a 0.90 threshold, buckling at
    // 0.83 of a 0.80 one -- so they resolve to the same panel and have to be
    // reduced to one zone. Without both triggers pointing at one panel nothing here
    // is a test: a single trigger is unique by construction.
    promotion::TierZero doubled;
    doubled.yieldStrength = ah36Steel().yieldStrength;
    for (double station : x) {
        const bool hot = std::abs(station) < 1e-9;
        GirderStress s;
        s.x = station;
        s.stressDeck = 1.0e8;
        s.stressKeel = -3.4e8;  // the keel is the worse fibre, so yield sites there
        s.utilisation = hot ? 0.99 : 0.20;
        doubled.stress.push_back(s);
        GirderBuckling b;
        b.x = station;
        b.utilisation = hot ? 0.83 : 0.20;
        b.deckInCompression = false;  // and so does buckling
        doubled.buckling.push_back(b);
    }
    const std::vector<promotion::Candidate> unique =
        promotion::candidates(structure, doubled, {}, criterion);
    std::printf("     yield at %.3f and buckling at %.3f on the same station and the same fibre"
                " become %zu candidate(s), trigger %s at score %.3f\n", 0.99, 0.83, unique.size(),
                unique.empty() ? "-" : promotion::name(unique.front().trigger),
                unique.empty() ? 0.0 : unique.front().score);
    testing::expectEqual("two triggers on one panel are one candidate",
                         static_cast<long long>(unique.size()), 1);
    expectTrue("and the strongest claim on it is the one that survives",
               unique.front().trigger == promotion::Trigger::Yield);
    expectNear("carrying its own score", unique.front().score, 0.99 / criterion.yieldPromote, 1e-12);
    // Both really did qualify, or "one candidate" is a statement about the other
    // one never firing.
    promotion::Criterion yieldOnly = criterion;
    yieldOnly.bucklePromote = 1e9;
    promotion::Criterion buckleOnly = criterion;
    buckleOnly.yieldPromote = 1e9;
    testing::expectEqual(
        "the yield trigger fires on its own",
        static_cast<long long>(promotion::candidates(structure, doubled, {}, buckleOnly).size()), 1);
    testing::expectEqual(
        "and so does the buckling one",
        static_cast<long long>(promotion::candidates(structure, doubled, {}, yieldOnly).size()), 1);
    expectTrue("and the two would have been ranked the other way round by raw utilisation",
               0.99 / criterion.yieldPromote > 0.83 / criterion.bucklePromote);
}

void testTheZoneLandsOnTheFibreThatIsInTrouble() {
    std::printf("\n   which fibre a girder trigger points at\n");
    const StructuralMesh& structure = ferryStructure();
    const HullGirderSection midship = hullGirderSection(structure, 0.0);
    promotion::Criterion criterion;
    criterion.mesh.radius = 3.0;
    std::printf("     midship: keel at %.2f m, deck at %.2f m, neutral axis %.2f m\n",
                midship.zKeel, midship.zDeck, midship.neutralAxis);
    expectTrue("the two fibres are far apart, so this is a real question",
               midship.zDeck - midship.zKeel > 10.0);

    // Buckling names its own compressed fibre.
    for (bool deck : {false, true}) {
        promotion::TierZero tier;
        tier.yieldStrength = ah36Steel().yieldStrength;
        for (double station : stationsAlong(-45.0, 45.0, 31)) {
            GirderBuckling b;
            b.x = station;
            b.utilisation = station == 0.0 ? 0.95 : 0.20;
            b.deckInCompression = deck;
            tier.buckling.push_back(b);
        }
        const std::vector<promotion::Candidate> found =
            promotion::candidates(structure, tier, {}, criterion);
        expectTrue("the buckling trigger found its station", found.size() == 1);
        const double wanted = deck ? midship.zDeck : midship.zKeel;
        std::printf("     %s in compression puts the zone at z = %.2f m (fibre at %.2f)\n",
                    deck ? "deck" : "keel", found.front().impact.z, wanted);
        expectTrue("and the zone lands on the compressed fibre",
                   std::abs(found.front().impact.z - wanted) < 3.0);
    }

    // Collapse takes it from the sign of the moment: hogging arches the hull and
    // compresses the keel. Getting this backwards puts every zone at the wrong end
    // of the section and nothing else would say so.
    for (double moment : {+1.0e9, -1.0e9}) {
        promotion::TierZero tier;
        tier.yieldStrength = ah36Steel().yieldStrength;
        for (double station : stationsAlong(-45.0, 45.0, 31)) {
            StrengthStation s;
            s.x = station;
            s.appliedMoment = moment;
            s.ultimateMoment = moment * (station == 0.0 ? 1.05 : 5.0);
            s.margin = std::abs(s.ultimateMoment / s.appliedMoment);
            tier.strength.push_back(s);
        }
        const std::vector<promotion::Candidate> found =
            promotion::candidates(structure, tier, {}, criterion);
        expectTrue("the collapse trigger found its station", found.size() == 1);
        const double wanted = moment > 0 ? midship.zKeel : midship.zDeck;
        std::printf("     %s puts the zone at z = %.2f m (the compressed fibre is at %.2f)\n",
                    moment > 0 ? "hogging" : "sagging", found.front().impact.z, wanted);
        expectTrue("hogging compresses the keel and sagging the deck",
                   std::abs(found.front().impact.z - wanted) < 3.0);
    }

    // Yield takes it from whichever fibre carries the larger stress, which on a
    // section whose two moduli differ is not the same as the compressed one.
    for (int worse = 0; worse < 2; ++worse) {
        promotion::TierZero tier;
        tier.yieldStrength = ah36Steel().yieldStrength;
        for (double station : stationsAlong(-45.0, 45.0, 31)) {
            GirderStress s;
            s.x = station;
            s.stressDeck = worse == 0 ? 3.4e8 : 1.0e8;
            s.stressKeel = worse == 0 ? -1.0e8 : -3.4e8;
            s.utilisation = station == 0.0 ? 0.96 : 0.20;
            tier.stress.push_back(s);
        }
        const std::vector<promotion::Candidate> found =
            promotion::candidates(structure, tier, {}, criterion);
        expectTrue("the yield trigger found its station", found.size() == 1);
        const double wanted = worse == 0 ? midship.zDeck : midship.zKeel;
        expectTrue("and lands on the fibre carrying the larger stress",
                   std::abs(found.front().impact.z - wanted) < 3.0);
    }
}

void testCandidatesAreRankedByHowFarPastTheirOwnThresholdTheyAre() {
    std::printf("\n   ranking: past its own threshold, not raw utilisation\n");
    const StructuralMesh& structure = ferryStructure();
    const std::vector<double> x = stationsAlong(-45.0, 45.0, 31);
    promotion::Criterion criterion;
    criterion.mesh.radius = 3.0;

    // Yield at 0.92 against a 0.90 threshold is 1.022 of its capacity; buckling at
    // 0.85 against 0.80 is 1.063 of its. Raw utilisation ranks the yield station
    // first and the criterion has to rank the buckling one first -- which is the
    // whole reason the score is normalised, since a buckling utilisation and a
    // yield utilisation are not the same currency.
    promotion::TierZero tier;
    tier.yieldStrength = ah36Steel().yieldStrength;
    for (double station : x) {
        GirderStress s;
        s.x = station;
        s.stressDeck = 3.4e8;
        s.stressKeel = -1.0e8;
        s.utilisation = std::abs(station + 30.0) < 1e-9 ? 0.92 : 0.20;
        tier.stress.push_back(s);
        GirderBuckling b;
        b.x = station;
        b.utilisation = std::abs(station - 30.0) < 1e-9 ? 0.85 : 0.20;
        b.deckInCompression = false;
        tier.buckling.push_back(b);
    }
    const std::vector<promotion::Candidate> found =
        promotion::candidates(structure, tier, {}, criterion);
    expectTrue("both stations qualified", found.size() == 2);
    std::printf("     yield %.2f/%.2f = %.3f at x = %.0f; buckling %.2f/%.2f = %.3f at x = %.0f;"
                " ranked %s first\n", 0.92, criterion.yieldPromote, found.back().score,
                found.back().x, 0.85, criterion.bucklePromote, found.front().score,
                found.front().x, promotion::name(found.front().trigger));
    expectTrue("the higher raw utilisation is the yield one",
               0.92 > 0.85);
    expectTrue("but the buckling station is further past its own threshold and ranks first",
               found.front().trigger == promotion::Trigger::Buckling);
    expectTrue("and the scores say why", found.front().score > found.back().score);
}

// --- 3. Idempotence and determinism -------------------------------------------

void testPromotionIsIdempotentAndDeterministic() {
    std::printf("\n   the same loads promote the same patches\n");
    const StructuralMesh& structure = ferryStructure();
    const std::vector<double> x = stationsAlong(-45.0, 45.0, 31);
    std::vector<double> peaked;
    for (double station : x)
        peaked.push_back(0.95 * std::exp(-((station - 12.0) * (station - 12.0)) / (2.0 * 12.0 * 12.0)));
    const promotion::TierZero tier = syntheticBuckling(x, peaked);

    promotion::Criterion criterion;
    criterion.mesh.radius = 3.0;
    promotion::Promoter first(criterion), second(criterion);
    for (int i = 0; i < 12; ++i) {
        first.review(structure, tier);
        second.review(structure, tier);
    }
    expectTrue("something was promoted, so this is not a comparison of two empty sets",
               !first.active().empty());
    testing::expectEqual("two promoters given the same reviews agree on the count",
                         static_cast<long long>(first.active().size()),
                         static_cast<long long>(second.active().size()));
    for (std::size_t i = 0; i < first.active().size(); ++i)
        testing::expectEqual("and on which panel", first.active()[i].panel,
                             second.active()[i].panel);

    // Idempotent: after the dwell, an unchanging load changes nothing. The check
    // is on the *cumulative* counters, so a promote-demote-promote cycle that left
    // the same set behind would still fail it.
    const int promotionsAfterSettling = first.promotions();
    for (int i = 0; i < 20; ++i) first.review(structure, tier);
    testing::expectEqual("twenty more reviews of the same load promote nothing further",
                         first.promotions(), promotionsAfterSettling);
    testing::expectEqual("and demote nothing", first.demotions(), 0);
    std::printf("     settled at %zu zone(s) after %d promotions and %d demotions in 32 reviews\n",
                first.active().size(), first.promotions(), first.demotions());
}

// --- 4. Chatter, against its own negative control -------------------------------
//
// Two oscillations, because the two mechanisms catch different things. Each is
// run with the mechanism on and off, and the off case has to chatter -- a test
// that only ran the "on" case would pass just as well against a signal that was
// never going to cross the threshold.

void testHysteresisAndDwellPreventChatter() {
    std::printf("\n   chatter, and the two mechanisms against it\n");
    const StructuralMesh& structure = ferryStructure();
    const std::vector<double> x = stationsAlong(-45.0, 45.0, 31);

    // A profile with one hot station, whose level is driven by `peak`.
    const auto profile = [&](double peak) {
        std::vector<double> u;
        for (double station : x) u.push_back(station == 0.0 ? peak : 0.20);
        return syntheticBuckling(x, u);
    };

    promotion::Criterion criterion;
    criterion.mesh.radius = 3.0;

    // (a) A small oscillation *inside* the hysteresis band: 0.82 +- 0.03 about a
    // promote threshold of 0.80 and a hold of 0.70. Dwell alone cannot help --
    // the signal stays above the threshold for many consecutive reviews at a
    // time -- so this is the hysteresis on its own.
    const auto runSlow = [&](const promotion::Criterion& c) {
        promotion::Promoter promoter(c);
        for (int i = 0; i < 40; ++i) {
            const double level = 0.82 + 0.03 * std::sin(2.0 * kPi * i / 8.0);
            promoter.review(structure, profile(level));
        }
        return promoter;
    };
    promotion::Criterion noBand = criterion;
    noBand.buckleHold = noBand.bucklePromote;  // no hysteresis at all
    const promotion::Promoter withBand = runSlow(criterion);
    const promotion::Promoter without = runSlow(noBand);
    std::printf("     0.82 +- 0.03 across a 0.80 threshold: %d promotions with the band,"
                " %d without\n", withBand.promotions(), without.promotions());
    expectTrue("without a hysteresis band the zone chatters", without.promotions() >= 4);
    testing::expectEqual("with one it is promoted once and stays", withBand.promotions(), 1);
    testing::expectEqual("and is never demoted", withBand.demotions(), 0);

    // (b) A fast oscillation *wider* than the band: 0.75 +- 0.20, so the signal
    // genuinely leaves the hold level every other review. Hysteresis alone cannot
    // help; the dwell can, because the excursion is shorter than it.
    const auto runFast = [&](const promotion::Criterion& c) {
        promotion::Promoter promoter(c);
        for (int i = 0; i < 40; ++i)
            promoter.review(structure, profile(i % 2 == 0 ? 0.95 : 0.55));
        return promoter;
    };
    promotion::Criterion noDwell = criterion;
    noDwell.dwell = 1;
    noDwell.hold = 1;
    const promotion::Promoter dwelt = runFast(criterion);
    const promotion::Promoter undwelt = runFast(noDwell);
    std::printf("     0.95/0.55 alternating: %d promotions and %d demotions with a dwell of %d,"
                " %d and %d without\n", dwelt.promotions(), dwelt.demotions(), criterion.dwell,
                undwelt.promotions(), undwelt.demotions());
    expectTrue("with no dwell the zone is rebuilt every other review",
               undwelt.promotions() >= 10 && undwelt.demotions() >= 10);
    testing::expectEqual("a dwell of two never lets it promote at all", dwelt.promotions(), 0);

    // Guard against the dwell simply switching the criterion off: the same dwell
    // against a load that *stays* up promotes, and promptly.
    promotion::Promoter sustained(criterion);
    for (int i = 0; i < 6; ++i) sustained.review(structure, profile(0.95));
    testing::expectEqual("but a load that stays up is promoted", sustained.promotions(), 1);
    testing::expectEqual("on the review the dwell is satisfied",
                         sustained.active().front().promotedAtReview, criterion.dwell);

    // (c) And it is dropped after exactly `hold` reviews of not qualifying -- not
    // on the first, which would throw away a zone's plastic history for one quiet
    // tick, and not never.
    for (int i = 0; i < criterion.hold - 1; ++i) {
        sustained.review(structure, profile(0.05));
        testing::expectEqual("a zone survives a quiet review", sustained.demotions(), 0);
    }
    sustained.review(structure, profile(0.05));
    testing::expectEqual("and is dropped on the hold-th", sustained.demotions(), 1);
    testing::expectEqual("leaving nothing running",
                         static_cast<long long>(sustained.active().size()), 0);
    std::printf("     a zone whose load has gone survives %d quiet reviews and is dropped on"
                " the %dth\n", criterion.hold - 1, criterion.hold);
}

// --- 5. Cost: linear in the number of zones, and bounded ------------------------

void testCostIsLinearInTheNumberOfZones() {
    std::printf("\n   cost: n patches cost n times one patch\n");
    // The prediction first, because it is what the budget is spent against.
    const StructuralMesh strip = flatStrip(1.2, 0.8, 0.012, 3, 2);
    const zone::Patch one = zone::buildPatch(strip, {0, 0, 0}, flatParams(2));
    const double cost = zone::estimatedCost(one);
    expectTrue("a patch costs something", cost > 0);
    for (int n : {1, 2, 5, 17}) {
        double total = 0;
        for (int i = 0; i < n; ++i) total += zone::estimatedCost(one);
        expectNear("the predicted cost of n identical zones is n times one",
                   total, n * cost, 1e-9 * n * cost);
    }

    // And measured, because a prediction that is linear by construction proves
    // nothing about the solver. Two patches solved must take about twice one --
    // the assertion is two orders loose per the rule `test_plasticity.cpp` wrote
    // down, but it would still catch anything with a global assembly in it.
    zone::SolveParams solve;
    solve.indenter.halfLength = 0.06;
    solve.indenter.halfWidth = 1e3;
    solve.indenter.speed = 20.0;
    solve.indenter.stopAt = 0.01;
    // **The minimum of several runs, not one sample.** A wall clock measures the
    // machine as much as the code, and this assertion used to take one timing of
    // each and compare them: a single scheduling hiccup on the *one-solve* sample
    // pushed the ratio under 1.4 and failed the run. It did exactly that to three
    // agents working in parallel on this box -- one spent eighteen reproduction
    // runs and a second full gate on it -- and it sits in the gate step whose whole
    // purpose is detecting flakiness, which is the worst possible place for a test
    // that fails at random. A test that cries wolf under load teaches everyone to
    // discount gate failures, and that is the habit that lets a real one through.
    //
    // The minimum is the right statistic and not a widened tolerance: contention
    // can only ever make a run *slower*, so the fastest of several is the closest
    // estimate of the unloaded cost, and it converges from above rather than
    // wandering. The band stays where it was, because it was never too tight --
    // the sampling was too thin. Timing is still what is measured, deliberately:
    // the point is to catch a global assembly, and a step count would be identical
    // either way and so could not see one.
    const auto timeSolves = [&](int count) {
        double best = std::numeric_limits<double>::infinity();
        for (int attempt = 0; attempt < 3; ++attempt) {
            double seconds = 0;
            for (int i = 0; i < count; ++i) {
                zone::Solver solver(one, plasticity::shipSteel(), solve);
                seconds += solver.run().wallSeconds;
            }
            best = std::min(best, seconds);
        }
        return best;
    };
    const double single = timeSolves(1);
    const double pair = timeSolves(2);
    std::printf("     %zu elements: one solve %.3f s, two %.3f s, ratio %.2f (best of 3 each)\n",
                one.elementCount(), single, pair, single > 0 ? pair / single : 0.0);
    expectTrue("a solve takes a measurable time at all", single > 0);
    expectTrue("two zones cost about two zones", pair > 1.4 * single && pair < 3.0 * single);

    // The budget is in elements, and it binds. A criterion that promoted past it
    // would be unaffordable rather than wrong.
    const StructuralMesh& structure = ferryStructure();
    const std::vector<double> x = stationsAlong(-45.0, 45.0, 31);
    std::vector<double> hot;
    for (double station : x) hot.push_back(0.20 + 0.75 * std::exp(-std::abs(station) / 4.0));
    const promotion::TierZero tier = syntheticBuckling(x, hot);

    promotion::Criterion tight;
    tight.mesh.radius = 3.0;
    tight.elementBudget = 10;  // less than one zone
    promotion::Promoter starved(tight);
    promotion::Review review;
    for (int i = 0; i < 4; ++i) review = starved.review(structure, tier);
    expectTrue("a qualifying zone was found", !review.considered.empty());
    testing::expectEqual("but the budget refused it", starved.promotions(), 0);
    expectTrue("and said so rather than silently dropping it", !review.problems.empty());

    promotion::Criterion generous = tight;
    generous.elementBudget = 100000;
    promotion::Promoter funded(generous);
    for (int i = 0; i < 4; ++i) funded.review(structure, tier);
    expectTrue("the same load with a budget promotes", funded.promotions() > 0);
    expectTrue("and the active cost is the sum of the zones'",
               funded.activeCost() > 0 && funded.activeElements() > 0);

    // **The budget is spent against what is already running**, not against this
    // review's promotions alone. Two hot spots far apart, and room for one.
    std::vector<double> twoSpots;
    for (double station : x)
        twoSpots.push_back(0.20 + 0.75 * std::exp(-std::abs(std::abs(station) - 30.0) / 2.0));
    const promotion::TierZero spots = syntheticBuckling(x, twoSpots);
    promotion::Criterion oneZone = generous;
    promotion::Promoter unbounded(oneZone);
    for (int i = 0; i < 4; ++i) unbounded.review(structure, spots);
    expectTrue("both hot spots are promotable when there is room",
               unbounded.active().size() == 2);
    oneZone.elementBudget = unbounded.active().front().elements + 1;
    promotion::Promoter bounded(oneZone);
    for (int i = 0; i < 6; ++i) bounded.review(structure, spots);
    std::printf("     two hot spots %d elements each: %zu promoted with room for both,"
                " %zu with a budget of %d\n", unbounded.active().front().elements,
                unbounded.active().size(), bounded.active().size(), oneZone.elementBudget);
    testing::expectEqual("a budget with room for one promotes one and stops",
                         static_cast<long long>(bounded.active().size()), 1);

    // And a zone is promoted once, however many reviews agree with it. Separation
    // would hide this -- a second zone on the same panel is zero metres from the
    // first -- so it is switched off, which is the only way to ask the question.
    promotion::Criterion overlapping = generous;
    overlapping.separation = 0.0;
    promotion::Promoter once(overlapping);
    for (int i = 0; i < 10; ++i) once.review(structure, tier);
    for (std::size_t i = 0; i < once.active().size(); ++i)
        for (std::size_t j = i + 1; j < once.active().size(); ++j)
            expectTrue("no panel carries two zones even with separation switched off",
                       once.active()[i].panel != once.active()[j].panel);
}

// --- 6. What the criterion costs to run ----------------------------------------
//
// Two numbers, and they are three orders apart. The **decision** is cheap enough
// to run every tick; the **Tier-0 answer it reads** is not, and that is the
// honest form of "it is explicitly not run every tick".

void testTheDecisionIsCheapAndTheTierZeroAnswerIsNot() {
    std::printf("\n   what promotion costs to decide, and what it costs to know\n");
    const Ship ferry = ferryAfloat();
    const StructuralMesh& structure = ferryStructure();
    const Scantlings scantlings = ferryScantlings();
    const Crest crest = crestAmidships(ferry, 3.0);

    const promotion::TierZero tier =
        promotion::tierZero(ferry, crest.sea, structure, scantlings);
    promotion::TierZeroParams cheap;
    cheap.collapse = false;
    const promotion::TierZero withoutCollapse =
        promotion::tierZero(ferry, crest.sea, structure, scantlings, cheap);

    promotion::Criterion criterion;
    criterion.mesh.radius = 3.0;
    promotion::Promoter promoter(criterion);
    double worst = 0;
    for (int i = 0; i < 20; ++i)
        worst = std::max(worst, promoter.review(structure, tier).microseconds);

    std::printf("     Tier-0 %.0f ms, of which the Smith sweep is %.0f ms;"
                " the decision itself %.0f us over %zu panels\n",
                tier.seconds * 1e3, (tier.seconds - withoutCollapse.seconds) * 1e3, worst,
                structure.panels.size());
    // Two orders loose, deliberately: this is a statement about which side of a
    // tick each of them is on, not a benchmark.
    expectTrue("the decision is tick-cheap", worst < 20000.0);
    expectTrue("the Tier-0 answer it reads is not", tier.seconds > 0.01);
    expectTrue("and the progressive-collapse sweep is most of it",
               tier.seconds > 2.0 * withoutCollapse.seconds);
    expectTrue("turning the sweep off turns the collapse trigger off with it",
               withoutCollapse.collapseUtilisation == 0.0 && tier.collapseUtilisation > 0);
}

// --- 7. The pre-load a patch is handed, against closed forms --------------------

void testAPreLoadedPatchCarriesTheStressItWasHanded() {
    std::printf("\n   the pre-load: the stress asked for, the energy it stores, and no motion\n");
    const double lengthX = 1.6, spanY = 0.8, thickness = 0.020, stress = 84.0e6;
    const StructuralMesh strip = flatStrip(lengthX, spanY, thickness, 2, 1);
    const zone::Patch patch = zone::buildPatch(strip, {0, 0, 0}, flatParams(4));
    expectTrue("the patch meshed", !patch.empty());

    for (bool plastic : {false, true}) {
        zone::SolveParams solve;
        solve.plastic = plastic;
        solve.indenter.halfLength = -1;  // no punch at all
        solve.indenter.halfWidth = -1;
        solve.preload.stress = stress;
        zone::Solver solver(patch, plasticity::shipSteel(), solve);

        double mean[6];
        solver.meanStress(mean);
        std::printf("     %-8s mean stress xx %.3f MPa, yy %.3f, zz %.3f, shears %.4f MPa\n",
                    plastic ? "plastic" : "elastic", mean[0] / 1e6, mean[1] / 1e6, mean[2] / 1e6,
                    std::max({std::abs(mean[3]), std::abs(mean[4]), std::abs(mean[5])}) / 1e6);
        expectNear("the patch carries exactly the longitudinal stress it was handed", mean[0],
                   stress, 0.005 * stress);
        // Uniaxial is the claim, so the other five components are the test. A
        // pre-load applied as a *strain* without the Poisson contraction would put
        // 0.43 of it into yy and zz and this is what would say so.
        for (int i = 1; i < 6; ++i)
            expectTrue("and nothing else: the state is uniaxial",
                       std::abs(mean[i]) < 0.002 * stress);

        // U = sigma^2 V / 2E, exactly, and it is in the account as an initial
        // condition rather than as work the indenter did.
        const double volume = patch.area * patch.thickness;
        expectNear("its stored energy is sigma^2 V / 2E",
                   solver.result().initialStrainEnergy,
                   stress * stress * volume / (2.0 * ah36Steel().youngsModulus),
                   0.01 * stress * stress * volume / (2.0 * ah36Steel().youngsModulus));
    }

    // **And it does not move.** The strain field is compatible and equilibrated,
    // so a clamped patch handed it and nothing else sits still. Run for 20 ms,
    // which is several fundamental periods of a 0.8 m span of 20 mm plate --
    // "it did not move" over a tenth of a period is a statement about not having
    // waited.
    zone::SolveParams settle;
    settle.plastic = false;
    settle.indenter.halfLength = -1;
    settle.indenter.halfWidth = -1;
    settle.preload.stress = stress;
    settle.duration = 0.02;
    zone::Solver solver(patch, plasticity::shipSteel(), settle);
    const std::vector<double> start = solver.position();
    solver.run();
    double moved = 0;
    for (std::size_t n = 0; n < patch.nodeCount(); ++n) {
        double d = 0;
        for (int k = 0; k < 3; ++k) {
            const double e = solver.position()[n * 3 + static_cast<std::size_t>(k)] -
                             start[n * 3 + static_cast<std::size_t>(k)];
            d += e * e;
        }
        moved = std::max(moved, std::sqrt(d));
    }
    std::printf("     after %d steps (%.0f ms) the worst node has moved %.2e m and the patch holds"
                " %.4g J of kinetic energy against %.4g J stored\n",
                solver.result().steps, solver.result().time * 1e3, moved,
                solver.result().kinetic, solver.result().strainEnergy);
    expectTrue("a pre-loaded patch under no other load does not move", moved < 1e-6);
    expectTrue("and carries no kinetic energy worth the name",
               solver.result().kinetic < 1e-6 * solver.result().strainEnergy);

    // The gradient form -- which is the one a hull girder actually produces, since
    // the stress varies linearly from the neutral axis. A patch standing on its
    // edge, so the gradient runs *across* it and the field is genuinely a bending
    // one rather than a uniform stress in disguise.
    StructuralMesh upright;
    upright.materials.push_back(ah36Steel());
    for (int i = 0; i < 2; ++i) {
        PlatePanel p;
        const double z0 = 4.0 + 0.8 * i, z1 = z0 + 0.8;
        p.corner[0] = {-0.8, 0.0, z0};
        p.corner[1] = {0.8, 0.0, z0};
        p.corner[2] = {0.8, 0.0, z1};
        p.corner[3] = {-0.8, 0.0, z1};
        p.thickness = thickness;
        p.role = PanelRole::Shell;
        upright.panels.push_back(p);
    }
    zone::MeshParams uprightMesh;
    uprightMesh.radius = 1e3;
    uprightMesh.subdivision = 4;
    uprightMesh.outward = {0, -1, 0};
    const zone::Patch wall = zone::buildPatch(upright, {0.0, 0.0, 4.8}, uprightMesh);
    expectTrue("the upright patch meshed", !wall.empty());

    zone::SolveParams graded;
    graded.plastic = false;
    graded.indenter.halfLength = -1;
    graded.indenter.halfWidth = -1;
    graded.preload.gradient = 1.0e7;   // Pa/m
    graded.preload.reference = -1.0;   // m; the patch centre at z = 4.8 sees 58 MPa
    graded.duration = 0.02;
    zone::Solver gradedSolver(wall, plasticity::shipSteel(), graded);
    double mean[6];
    gradedSolver.meanStress(mean);
    std::printf("     a gradient of %.0f MPa/m about z = %.0f puts %.2f MPa through a patch"
                " centred at z = %.2f, varying %.1f MPa across it\n",
                graded.preload.gradient / 1e6, graded.preload.reference, mean[0] / 1e6,
                wall.centre.z, graded.preload.gradient * 1.6 / 1e6);
    expectNear("a gradient pre-load puts gradient x height through the patch", mean[0],
               graded.preload.at(wall.centre.z), 0.005 * graded.preload.at(wall.centre.z));
    expectTrue("and the stress really does vary across it, so this is not a uniform one",
               graded.preload.gradient * 1.6 > 0.2 * mean[0]);

    // And a *bending* pre-load is equilibrated too. It is a different statement
    // from the uniform case: the displacement field carries a curvature term in
    // x^2, and a field missing it is not compatible with the strain it claims, so
    // the patch would spring the moment it was let go.
    const std::vector<double> gradedStart = gradedSolver.position();
    gradedSolver.run();
    double gradedMoved = 0;
    for (std::size_t n = 0; n < wall.nodeCount(); ++n) {
        double d = 0;
        for (int k = 0; k < 3; ++k) {
            const double e = gradedSolver.position()[n * 3 + static_cast<std::size_t>(k)] -
                             gradedStart[n * 3 + static_cast<std::size_t>(k)];
            d += e * e;
        }
        gradedMoved = std::max(gradedMoved, std::sqrt(d));
    }
    std::printf("     after %.0f ms under the gradient alone the worst node has moved %.2e m\n",
                gradedSolver.result().time * 1e3, gradedMoved);
    expectTrue("a patch under a bending pre-load alone does not move either", gradedMoved < 1e-6);
}

// --- 8. What a pre-load spends, exactly ----------------------------------------

void testAPreLoadSpendsYieldCapacityExactly() {
    std::printf("\n   a pre-load spends yield capacity, and the amount is closed form\n");
    const StructuralMesh strip = flatStrip(1.6, 0.8, 0.020, 2, 1);
    const zone::Patch patch = zone::buildPatch(strip, {0, 0, 0}, flatParams(2));
    const double yieldStrength = 355.0e6;  // `plasticity::shipSteel()`'s own initial yield

    const auto yieldsAt = [&](double preload) {
        zone::SolveParams solve;
        solve.indenter.halfLength = -1;
        solve.indenter.halfWidth = -1;
        solve.preload.stress = preload;
        zone::Solver solver(patch, plasticity::shipSteel(), solve);
        return anyPointYielded(solver);
    };

    expectTrue("a patch pre-loaded to 0.99 of yield has not yielded",
               !yieldsAt(0.99 * yieldStrength));
    expectTrue("one pre-loaded to 1.01 of yield has", yieldsAt(1.01 * yieldStrength));
    expectTrue("and the same in compression", !yieldsAt(-0.99 * yieldStrength) &&
                                                  yieldsAt(-1.01 * yieldStrength));
    std::printf("     the pre-load alone yields the patch between 0.99 and 1.01 of sigma_y"
                " = %.0f MPa\n", yieldStrength / 1e6);

    // So the capacity a zone starts with is `sigma_y - sigma_0` and not `sigma_y`.
    // On the reference ferry at 84 MPa of deck stress that is 76% of it, which is
    // the whole of the argument for handing a zone the girder's stress.
    std::printf("     a zone handed the ferry's 84 MPa therefore starts with %.1f%% of the"
                " uniaxial capacity an unloaded one claims\n",
                100.0 * (yieldStrength - 84.0e6) / yieldStrength);
}

void testAPreLoadedZoneYieldsEarlierUnderAPunch() {
    std::printf("\n   a pre-loaded zone under a punch\n");
    // Six millimetre plating on a 0.3 m span, the fixture `test_zone.cpp` tears.
    // Subdivision four rather than two, because at two the answer is one element
    // dying and the tearing depth moves several per cent with the mesh.
    const StructuralMesh strip = flatStrip(0.6, 0.3, 0.006, 3, 1);
    const zone::Patch patch = zone::buildPatch(strip, {0.05, 0, 0}, flatParams(4));

    struct Outcome {
        double yieldAt = 0, tearAt = 0;
    };
    const auto drive = [&](double preload, double stopAt) {
        zone::SolveParams solve;
        solve.indenter.halfLength = 0.03;
        solve.indenter.halfWidth = 0.08;
        solve.indenter.speed = 20.0;
        solve.indenter.rampTime = 5.0e-4;
        solve.indenter.stopAt = stopAt;
        solve.preload.stress = preload;
        zone::Solver solver(patch, plasticity::shipSteel(), solve);
        Outcome out;
        while (solver.step()) {
            if (out.yieldAt == 0 && anyPointYielded(solver))
                out.yieldAt = solver.result().penetration;
            if (out.tearAt == 0 && solver.result().tornElements > 0)
                out.tearAt = solver.result().penetration;
        }
        return out;
    };

    // Yield onset only: the punch is stopped as soon as the unloaded case has
    // yielded, which is what keeps this affordable.
    const Outcome none = drive(0.0, 0.008);
    const Outcome loaded = drive(0.90 * 355.0e6, 0.008);
    std::printf("     first yield: unloaded at %.5f m, pre-loaded to 0.9 sigma_y at %.5f m"
                " (%.3f of it)\n", none.yieldAt, loaded.yieldAt,
                none.yieldAt > 0 ? loaded.yieldAt / none.yieldAt : 0.0);
    expectTrue("the unloaded patch does reach first yield", none.yieldAt > 0);
    expectTrue("a patch already at 0.9 of yield yields far sooner",
               loaded.yieldAt > 0 && loaded.yieldAt < 0.25 * none.yieldAt);

    // **And the correction this measurement makes.** The pre-strain is *elastic*
    // -- 0.9 sigma_y is 1.6e-3 -- against a failure strain around 0.15, so it is
    // two orders below what tears the plate and it barely moves the tear. "A
    // pre-loaded zone fails earlier" is true of yielding and false of tearing.
    const Outcome tornNone = drive(0.0, 0.06);
    const Outcome tornLoaded = drive(0.90 * 355.0e6, 0.06);
    const double ratio = tornNone.tearAt > 0 ? tornLoaded.tearAt / tornNone.tearAt : 0.0;
    std::printf("     first tear:  unloaded at %.5f m, pre-loaded at %.5f m (%.4f of it) --"
                " the pre-strain is 1.6e-3 against a failure strain of 0.15\n",
                tornNone.tearAt, tornLoaded.tearAt, ratio);
    expectTrue("both cases tore, so this compares two tears", tornNone.tearAt > 0 &&
                                                                  tornLoaded.tearAt > 0);
    expectTrue("and the pre-load moves the tear by a few per cent, not by a factor",
               ratio > 0.9 && ratio < 1.1);
}

// --- 9. The pre-load comes from the girder, not from a parameter ----------------

void testThePreLoadIsTheGirdersOwnStress() {
    std::printf("\n   where the pre-load comes from\n");
    const Ship ferry = ferryAfloat();
    const StructuralMesh& structure = ferryStructure();
    const Crest crest = crestAmidships(ferry, 3.0);
    const HullGirder girder = hullGirder(ferry, crest.sea, 41);
    expectTrue("she is hogging on the crest", girder.hogging());

    zone::MeshParams mesh;
    mesh.radius = 2.0;
    mesh.subdivision = 2;
    const zone::Patch side = zone::buildPatch(structure, {0.0, -9.9, 8.0}, mesh);
    expectTrue("a patch meshed on her side", !side.empty());
    const promotion::PreloadCheck check = promotion::preloadFor(girder, structure, side);

    // Against `girderStress()`, which computes the same thing a completely
    // different way: M / Z at the deck, where Z = I / (z_deck - z_na).
    const std::vector<GirderStress> stresses =
        girderStress(girder, structure, ah36Steel().yieldStrength);
    const HullGirderSection section = hullGirderSection(structure, side.centre.x);
    double deckStress = 0;
    for (const GirderStress& s : stresses)
        if (std::abs(s.x - side.centre.x) < 1e-9) deckStress = s.stressDeck;
    std::printf("     M %.4e N m, I %.3f m4, neutral axis %.3f m; the pre-load reads %.2f MPa at"
                " the deck where M/Z says %.2f MPa\n", check.moment, section.secondMoment,
                check.neutralAxis, check.preload.at(section.zDeck) / 1e6, deckStress / 1e6);
    expectTrue("the deck stress is the 84 MPa the milestone publishes",
               deckStress > 80.0e6 && deckStress < 90.0e6);
    expectNear("the pre-load at the deck fibre is M/Z", check.preload.at(section.zDeck),
               deckStress, 1e-6 * std::abs(deckStress));
    expectNear("and it is zero at the neutral axis", check.preload.at(section.neutralAxis), 0.0,
               1.0);
    expectTrue("it is applied on plating that faces athwartships", check.applied);
    expectNear("whose normal lies in the athwartships plane", check.obliquity, 0.0, 1e-6);
    expectTrue("so no traction is left unbalanced on it", check.tractionError < 1.0);

    // Hogging tensions the deck and compresses the keel, and the pre-load has to
    // carry that sign or a zone below the neutral axis would be handed the wrong
    // one entirely.
    expectTrue("hogging tensions the deck", check.preload.at(section.zDeck) > 0);
    expectTrue("and compresses the keel", check.preload.at(section.zKeel) < 0);

    // A patch whose normal leans along the ship carries no hull girder stress
    // across it, and is refused rather than given one. A synthetic panel is used
    // rather than a real bulkhead, because the point is the geometry and not
    // whether this ferry happens to have a bulkhead where a zone would fit.
    StructuralMesh transverse;
    transverse.materials.push_back(ah36Steel());
    for (int i = 0; i < 2; ++i) {
        PlatePanel p;
        const double y0 = -1.0 + i, y1 = y0 + 1.0;
        p.corner[0] = {0.0, y0, 5.0};
        p.corner[1] = {0.0, y1, 5.0};
        p.corner[2] = {0.0, y1, 7.0};
        p.corner[3] = {0.0, y0, 7.0};
        p.thickness = 0.009;
        p.role = PanelRole::Bulkhead;
        transverse.panels.push_back(p);
    }
    zone::MeshParams bulkheadMesh;
    bulkheadMesh.radius = 1e3;
    bulkheadMesh.subdivision = 2;
    bulkheadMesh.role = PanelRole::Bulkhead;
    bulkheadMesh.outward = {1, 0, 0};
    const zone::Patch bulkhead = zone::buildPatch(transverse, {0.0, 0.0, 6.0}, bulkheadMesh);
    expectTrue("the bulkhead patch meshed", !bulkhead.empty());
    const promotion::PreloadCheck refused =
        promotion::preloadFor(girder, structure, bulkhead);
    std::printf("     a patch whose normal is along the ship: obliquity %.4f rad, applied %d\n",
                refused.obliquity, static_cast<int>(refused.applied));
    expectNear("its normal is a right angle out of the athwartships plane", refused.obliquity,
               0.5 * kPi, 1e-6);
    expectTrue("so it is refused a pre-load", !refused.applied && !refused.preload.active());
    expectTrue("and told why", !refused.problems.empty());

    // Between the two: a panel leaning `phi` out of the athwartships plane keeps
    // `sigma sin^2(phi)` of unbalanced traction on its own face, because the
    // traction is `sigma n_x` along x and its component on the face is another
    // factor of `n_x`. The square is the whole point -- linear in the normal would
    // over-state a raked bow's error by a factor of three at 20 degrees.
    for (double lean : {0.10, 0.19}) {
        StructuralMesh raked;
        raked.materials.push_back(ah36Steel());
        // A wall of plating standing on the athwartships direction rotated by
        // `lean` about z, so its outward normal is (-sin, -cos, 0) and leans
        // exactly `lean` out of the athwartships plane.
        const double c = std::cos(lean), s = std::sin(lean);
        for (int i = 0; i < 2; ++i) {
            PlatePanel p;
            const double z0 = 5.0 + i, z1 = z0 + 1.0;
            p.corner[0] = {-c, s, z0};
            p.corner[1] = {c, -s, z0};
            p.corner[2] = {c, -s, z1};
            p.corner[3] = {-c, s, z1};
            p.thickness = 0.012;
            p.role = PanelRole::Shell;
            raked.panels.push_back(p);
        }
        zone::MeshParams rakedMesh;
        rakedMesh.radius = 1e3;
        rakedMesh.subdivision = 2;
        rakedMesh.outward = {-s, -c, 0.0};
        const zone::Patch leaning = zone::buildPatch(raked, {0.0, 0.0, 6.0}, rakedMesh);
        expectTrue("the leaning patch meshed", !leaning.empty());
        const promotion::PreloadCheck oblique =
            promotion::preloadFor(girder, structure, leaning, 0.15);
        const double sine = std::sin(oblique.obliquity);
        std::printf("     a panel leaning %.2f rad: obliquity %.3f, %.1f kPa unbalanced of"
                    " %.1f MPa (sin^2 = %.4f), applied %d\n", lean, oblique.obliquity,
                    oblique.tractionError / 1e3, std::abs(oblique.surfaceStress) / 1e6, sine * sine,
                    static_cast<int>(oblique.applied));
        expectNear("the unbalanced traction is sigma sin^2 of the lean", oblique.tractionError,
                   std::abs(oblique.surfaceStress) * sine * sine,
                   1e-6 * std::abs(oblique.surfaceStress));
        expectTrue("and sin^2 is not sin, so the square is doing work",
                   sine * sine < 0.6 * sine);
        testing::expectEqual(lean < 0.15 ? "inside the limit it is pre-loaded"
                                         : "outside it, it is not",
                             static_cast<long long>(oblique.applied), lean < 0.15 ? 1 : 0);
    }
}

// --- 9b. What leaving the longitudinals at full strength was worth ---------------

// **The headline of the fibre-damage work: how wrong the section reduction was.**
//
// `06-roadmap.md` recorded the gap as "a collision that opens fourteen bays leaves
// their longitudinals at full strength, which is the un-conservative direction", and
// said nothing about how far un-conservative. This measures it, on the reference
// ferry, and the answer is not small: about a fifth of the strength a ram removes
// from the hull girder was invisible, and in sagging about a third.
//
// Driven synthetically -- the panels and the members are set to nothing rather than
// solved to nothing -- for the same reason `testTheReductionMakesTierZeroReportLessStrength`
// below is: what is under test is the *coupling* into Tier 0, and a solve would only
// make it slower and less severe. `testASolvedZoneNamesTheLongitudinalItTore` is the
// end-to-end half.
void testTheLongitudinalsARamDestroysAreWorthAFifthOfTheGirder() {
    std::printf("\n   what leaving a ram's longitudinals at full strength was worth\n");
    const StructuralMesh& structure = ferryStructure();
    const Scantlings scantlings = ferryScantlings();

    // A ram amidships on the port side at the sheer: every shell panel within 8 m,
    // and the shell longitudinals actually running through them.
    const Vec3 centre{0.0, -9.0, 13.0};
    promotion::SectionReduction plateOnly, both;
    double panelArea = 0, memberLength = 0;
    for (std::size_t i = 0; i < structure.panels.size(); ++i) {
        const PlatePanel& p = structure.panels[i];
        if (p.role != PanelRole::Shell) continue;
        if (length(p.centroid() - centre) > 8.0) continue;
        promotion::PanelDamage damage;
        damage.panel = static_cast<int>(i);
        damage.effectiveness = 0.0;
        plateOnly.panels.push_back(damage);
        both.panels.push_back(damage);
        panelArea += p.area();
    }
    // **The member selection is where this test could go vacuous, and it is checked
    // rather than trusted.** Picking members by distance from the centre alone
    // catches deck beams and bulkhead stiffeners metres inboard that the ram never
    // touched, and over-states the loss by construction -- measured, that selection
    // reported a 2.5x error where the honest one reports 1.2x. So a member counts
    // only if it is a shell longitudinal whose midpoint lands on a panel that was
    // opened, and the guard below is that the length picked up divides the opened
    // area by the ferry's own longitudinal spacing.
    for (std::size_t i = 0; i < structure.members.size(); ++i) {
        const StructuralMember& member = structure.members[i];
        if (member.role != MemberRole::Longitudinal) continue;
        const Vec3 mid = (member.a + member.b) * 0.5;
        bool onOpenedPlating = false;
        for (const promotion::PanelDamage& damage : plateOnly.panels)
            if (length(mid - structure.panels[static_cast<std::size_t>(damage.panel)].centroid()) <
                0.75) {
                onOpenedPlating = true;
                break;
            }
        if (!onOpenedPlating || length(mid - centre) > 8.0) continue;
        promotion::MemberDamage damage;
        damage.member = static_cast<int>(i);
        damage.effectiveness = 0.0;
        both.members.push_back(damage);
        memberLength += member.length();
    }

    const double spacing = panelArea / memberLength;
    std::printf("     the ram opens %.1f m2 of shell over %zu panels and %.1f m of longitudinal"
                " over %zu members -- %.3f m apart, against a design spacing of %.3f m\n",
                panelArea, plateOnly.panels.size(), memberLength, both.members.size(), spacing,
                scantlings.longitudinalSpacing);
    expectTrue("the sweep covers a real patch of her side", plateOnly.panels.size() > 40);
    expectTrue("and the longitudinals it picks up are the ones on that plating, at her own"
               " spacing", spacing > 0.5 * scantlings.longitudinalSpacing &&
                               spacing < 1.6 * scantlings.longitudinalSpacing);

    struct Row { double area, second, modulusDeck, hog, sag; };
    const auto measure = [&](const promotion::SectionReduction& reduction) {
        const StructuralMesh damaged = promotion::reduce(structure, reduction);
        const HullGirderSection section = hullGirderSection(damaged, 0.0);
        const std::vector<CollapseElement> elements = collapseElementsAt(damaged, scantlings, 0.0);
        return Row{section.area, section.secondMoment, section.modulusDeck,
                   collapseCurve(elements, 1.0).ultimateMoment,
                   std::abs(collapseCurve(elements, -1.0).ultimateMoment)};
    };
    const Row intact = measure(promotion::SectionReduction{});
    const Row plating = measure(plateOnly);
    const Row all = measure(both);

    std::printf("     %-14s %10s %10s %10s %13s %13s\n", "", "area (m2)", "I (m4)", "Zdeck (m3)",
                "M_hog (N m)", "M_sag (N m)");
    const auto say = [&](const char* label, const Row& row) {
        std::printf("     %-14s %10.5f %10.4f %10.5f %13.5e %13.5e\n", label, row.area, row.second,
                    row.modulusDeck, row.hog, row.sag);
    };
    say("intact", intact);
    say("plating only", plating);
    say("+ members", all);

    const auto lost = [](double damaged, double whole) { return 1.0 - damaged / whole; };
    struct Gap { const char* name; double plate, both_; };
    const Gap gaps[] = {
        {"section area", lost(plating.area, intact.area), lost(all.area, intact.area)},
        {"second moment", lost(plating.second, intact.second), lost(all.second, intact.second)},
        {"hogging ultimate moment", lost(plating.hog, intact.hog), lost(all.hog, intact.hog)},
        {"sagging ultimate moment", lost(plating.sag, intact.sag), lost(all.sag, intact.sag)},
    };
    double worstRatio = 1.0;
    for (const Gap& gap : gaps) {
        std::printf("     %-26s plating alone %.3f%%, with the longitudinals %.3f%%  (%.2fx)\n",
                    gap.name, 100.0 * gap.plate, 100.0 * gap.both_, gap.both_ / gap.plate);
        expectTrue("losing the longitudinals too costs strictly more, on every measure",
                   gap.both_ > gap.plate);
        worstRatio = std::max(worstRatio, gap.both_ / gap.plate);
    }
    // The direction is the whole point: the model that left the longitudinals alone
    // reported *more* strength left than there is, never less.
    expectTrue("and every quantity falls when the members go", all.area < plating.area &&
                                                                   all.second < plating.second &&
                                                                   all.hog < plating.hog &&
                                                                   all.sag < plating.sag);
    // Non-vacuous in the other direction: the error has to be worth reporting. A
    // ratio near 1.0 would mean the longitudinals never mattered and this whole task
    // bought nothing.
    std::printf("     so the un-conservative error was 1.2x on the section properties and %.2fx"
                " at worst -- about a fifth of what the ram takes out of her, and about a third"
                " of it in sagging, was invisible\n", worstRatio);
    expectTrue("the error is a fifth of the damage at least, not a rounding one",
               gaps[0].both_ > 1.15 * gaps[0].plate);
    expectTrue("and the sagging ultimate moment is where it bites hardest, because a panel that"
               " has lost its stiffener buckles far earlier than one that has merely thinned",
               gaps[3].both_ / gaps[3].plate > 1.35);
    // The reduction is not allowed to invent material: a member with no damage in the
    // list is left exactly as it was, bit for bit.
    const StructuralMesh damaged = promotion::reduce(structure, both);
    bool untouchedElsewhere = damaged.members.size() == structure.members.size();
    std::size_t changed = 0;
    for (std::size_t i = 0; untouchedElsewhere && i < structure.members.size(); ++i) {
        const bool named = std::any_of(both.members.begin(), both.members.end(),
                                       [&](const promotion::MemberDamage& d) {
                                           return static_cast<std::size_t>(d.member) == i;
                                       });
        const bool same =
            damaged.members[i].profile.webThickness == structure.members[i].profile.webThickness &&
            damaged.members[i].profile.flangeThickness ==
                structure.members[i].profile.flangeThickness;
        if (named && !same) ++changed;
        if (!named && !same) untouchedElsewhere = false;
    }
    expectTrue("a member nobody named comes back bit-identical", untouchedElsewhere);
    testing::expectEqual("and every member that was named is thinner",
                         static_cast<long long>(changed),
                         static_cast<long long>(both.members.size()));
    // A member reduced to nothing has to leave the section outright, not linger as
    // an epsilon of area -- the same trap the plating's `1e-6` snap exists for.
    expectTrue("and one reduced to nothing has no section left at all",
               profileSection(damaged.members[static_cast<std::size_t>(both.members[0].member)]
                                  .profile)
                       .area == 0.0);
}

// The end-to-end half: a real solve, a real tear, and the member named back to the
// structure. Plus the control -- `SolveParams::fiberFailure` false is the model that
// shipped before the fibres carried damage -- so the size of the zone-level error is
// a measurement rather than an inference from the section one.
void testASolvedZoneNamesTheLongitudinalItTore() {
    std::printf("\n   a solved zone names the longitudinal it tore\n");
    // **The strip is three times the zone on purpose.** A longitudinal runs on past
    // the damage, and the part of it nobody looked at has to count as intact -- the
    // same rule the plating has, and the defect `indentation.hpp` records when it is
    // missing, where the zone *radius* would decide how much of the ship is broken. A
    // fixture whose member the zone meshes end to end cannot tell the two apart:
    // mutation testing found the rule deleted and surviving on exactly that fixture.
    const double lengthX = 4.8, spanY = 0.8, thickness = 0.012, depth = 0.42;
    const StiffenerProfile bar = flatBar(0.200, 0.010);
    const StructuralMesh strip = flatStrip(lengthX, spanY, thickness, 12, 4, &bar);
    const plasticity::Material steel = plasticity::shipSteel();
    const double memberLength = strip.members[0].length();

    struct Run {
        double work = 0, fibreEnergy = 0, bundleForce = 0, peakPlastic = 0, peakDamage = 0;
        int tornElements = 0, tornFibers = 0;
        double tornVolume = 0, lostVolume = 0, meshedFraction = 0;
        double memberEffectiveness = 1, memberArea = 0;
        double firstPlateTear = 0, plateFailureStrain = 0;
        double lostSteelMass = 0, lostPlateMass = 0;
        std::size_t members = 0, panels = 0;
        std::size_t anonymousMembers = 0, anonymousProblems = 0;
        bool anonymousDiagnosed = false;
        std::size_t strangerMembers = 0;
        bool strangerDiagnosed = false;
    };
    const auto go = [&](bool allowFailure) {
        zone::MeshParams params = flatParams(2);
        params.radius = 0.9;
        params.stiffeners = zone::Stiffeners::Modelled;
        params.edge = zone::Edge::Clamped;
        const zone::Patch patch = zone::buildPatch(strip, {0, 0, 0}, params);
        expectTrue("the stiffened patch meshed with fibres on it",
                   !patch.empty() && !patch.stiffening.empty());

        zone::SolveParams solve;
        solve.indenter.halfLength = 0.08;
        solve.indenter.halfWidth = 1e3;
        solve.indenter.speed = 20.0;
        solve.indenter.rampTime = 1.0e-3;
        solve.indenter.stopAt = depth;
        solve.fiberFailure = allowFailure;
        zone::Solver solver(patch, steel, solve);
        const zone::SolveResult& result = solver.run();
        // **`run()` on a finished solve is idempotent**, and the fibre counters have
        // to be too. It steps nothing and re-runs `collectTorn`, so a count that
        // accumulated instead of being reset would double here -- which is exactly
        // what mutation testing found, on a counter that is otherwise only ever read
        // once. The plating's `tornElements` has the same property and is checked
        // alongside so the two cannot drift apart.
        const int tornOnce = result.tornFibers, elementsOnce = result.tornElements;
        const double volumeOnce = result.tornFiberVolume;
        solver.run();
        testing::expectEqual("running a finished solve again does not double the torn fibre count",
                    static_cast<long long>(result.tornFibers), static_cast<long long>(tornOnce));
        testing::expectEqual("nor the torn element count",
                    static_cast<long long>(result.tornElements),
                    static_cast<long long>(elementsOnce));
        expectNear("nor the volume of stiffener steel it says has gone",
                   result.tornFiberVolume, volumeOnce, 0.0);

        Run out;
        out.work = result.work;
        out.tornElements = result.tornElements;
        out.tornFibers = result.tornFibers;
        out.tornVolume = result.tornFiberVolume;
        out.fibreEnergy = solver.fiberStrainEnergy();
        for (const constraint::FiberState& state : solver.fiberState()) {
            out.peakPlastic = std::max(out.peakPlastic, state.equivalentPlasticStrain);
            out.peakDamage = std::max(out.peakDamage, state.damage);
        }
        // The plating's own failure strain, off its own geometry, and the strain the
        // *first* plate point actually tore at. The two together are the measurement
        // behind "the failure strains are not what separates the plating from the
        // stiffener" -- see below.
        {
            double nodes[solidshell::kDof], inPlane = 0, thick = 0;
            patch.mesh.gather(0, patch.mesh.position, nodes);
            solidshell::elementSize(nodes, &inPlane, &thick);
            out.plateFailureStrain =
                plasticity::regularisedFailureStrain(steel.failure, inPlane, thick);
        }
        out.firstPlateTear = std::numeric_limits<double>::infinity();
        for (const solidshell::ElementPlasticState& element : solver.elementState())
            for (int gp = 0; gp < solidshell::kGauss; ++gp)
                if (element.point[gp].failed)
                    out.firstPlateTear =
                        std::min(out.firstPlateTear, element.point[gp].equivalentPlasticStrain);
        // What the bundle of bars is putting into the plating right now, as the norm
        // of the nodal force it contributes. It is the honest scalar for "how much
        // load is this longitudinal still taking", and unlike the punch reaction it
        // is not an instantaneous sample of a ringing explicit solve.
        std::vector<double> force(patch.nodeCount() * 3, 0.0);
        std::vector<constraint::FiberState> copy = solver.fiberState();
        constraint::fiberForces(patch.stiffening, solver.restFibers(), solver.position(), steel,
                                &copy, force, allowFailure);
        double sum = 0;
        for (double f : force) sum += f * f;
        out.bundleForce = std::sqrt(sum);

        const promotion::SectionReduction reduction = promotion::reactionOf(strip, patch, solver);
        out.lostSteelMass = reduction.lostSteelMass;
        // The plating's share of that mass, computed here from the panel damages the
        // same reduction reports. What is left over must be the stiffener's, and the
        // stiffener's must be its lost volume times its density -- so the split is
        // checked rather than the total, and a member share with the wrong sense of
        // `effectiveness` cannot hide inside a plating figure ten times its size.
        for (const promotion::PanelDamage& damage : reduction.panels) {
            const PlatePanel& panel = strip.panels[static_cast<std::size_t>(damage.panel)];
            out.lostPlateMass += (1.0 - damage.effectiveness) * panel.area() * panel.thickness *
                                 strip.materials[0].density;
        }

        // **Fibres with no member index cannot be reduced, and must not be.** The
        // same solve read through a patch whose fibres have been made anonymous: the
        // fibre list, its order and the state all line up, only the label is gone.
        // `reactionOf` has to name no member and say so, because a zone that built
        // fibres without member indices would otherwise be indistinguishable from a
        // zone whose stiffeners are fine.
        zone::Patch anonymous = patch;
        for (constraint::Fiber& fibre : anonymous.stiffening.fiber) fibre.member = -1;
        const promotion::SectionReduction unnamed =
            promotion::reactionOf(strip, anonymous, solver);
        // And the other end of the same guard: a member index the structure does not
        // have. `buildPatch` cannot produce one, but a patch and a structure that
        // have drifted apart can, and the alternative to checking is indexing a
        // vector out of bounds. The panel path already carries the same check for
        // the same reason.
        zone::Patch stranger = patch;
        for (constraint::Fiber& fibre : stranger.stiffening.fiber) fibre.member = 9999;
        const promotion::SectionReduction lost =
            promotion::reactionOf(strip, stranger, solver);
        out.strangerMembers = lost.members.size();
        for (const std::string& problem : lost.problems)
            if (problem.find("does not have") != std::string::npos)
                out.strangerDiagnosed = true;
        out.anonymousMembers = unnamed.members.size();
        out.anonymousProblems = unnamed.problems.size();
        // **The diagnosis, not just the count.** Deleting the `member < 0` branch
        // leaves `-1` to fall through to the range check, where the cast to
        // `size_t` makes it enormous and the reduction reports "the patch names a
        // member the structure does not have" -- one problem either way, and a
        // caller told the wrong thing about its own mesh. Mutation testing found
        // that edit surviving a count-only assertion, so the text is what is
        // checked.
        for (const std::string& problem : unnamed.problems)
            if (problem.find("no member index") != std::string::npos)
                out.anonymousDiagnosed = true;

        out.panels = reduction.panels.size();
        out.members = reduction.members.size();
        out.memberEffectiveness = reduction.worstMemberEffectiveness;
        if (!reduction.members.empty()) {
            out.lostVolume = reduction.members[0].lostVolume;
            out.meshedFraction = reduction.members[0].meshedFraction;
        }
        out.memberArea =
            profileSection(promotion::reduce(strip, reduction).members[0].profile).area;
        return out;
    };

    const Run live = go(true);
    const Run immortal = go(false);
    const double failureStrain = constraint::fiberFailureStrain(steel.failure, 0.200, 0.010);
    const double whole = profileSection(bar).area;

    std::printf("     %-24s %12s %12s\n", "", "fibres tear", "control");
    const auto say = [&](const char* label, double a, double b, const char* fmt) {
        std::printf("     %-24s ", label);
        std::printf(fmt, a);
        std::printf(" ");
        std::printf(fmt, b);
        std::printf("\n");
    };
    say("punch work (J)", live.work, immortal.work, "%12.5e");
    say("torn plate elements", live.tornElements, immortal.tornElements, "%12.0f");
    say("torn fibres", live.tornFibers, immortal.tornFibers, "%12.0f");
    say("peak fibre eps_p", live.peakPlastic, immortal.peakPlastic, "%12.4f");
    say("fibre strain energy (J)", live.fibreEnergy, immortal.fibreEnergy, "%12.5e");
    say("fibre nodal force (N)", live.bundleForce, immortal.bundleForce, "%12.5e");
    say("member section left (m2)", live.memberArea, immortal.memberArea, "%12.5e");

    // **Vacuity, both halves.** The fibres have to be carrying real load, and the
    // case has to actually fail them. A patch where the plating governed would pass
    // on any implementation at all.
    expectTrue("the plating really tore", live.tornElements > 0);
    expectTrue("and the fibres really tore too", live.tornFibers > 0);
    expectTrue("and they were carrying a real share of the load, not a token one",
               live.bundleForce > 1.0e6);
    expectTrue("the two runs tore the same plating, so the difference is the stiffeners alone",
               live.tornElements == immortal.tornElements);

    // **What separates the plating from the stiffener is not their failure strains.**
    // They are within 12% of each other at these sizes. It is Rice-Tracey: a bar's
    // triaxiality is exactly the reference so its multiplier is exactly one, while a
    // plate point under a punch is somewhere else entirely and tears at a small
    // fraction of its own regularised strain. Asserting it here rather than only
    // writing it in `02-simulation.md` -- nothing tests a comment.
    std::printf("     the plate's own eps_f is %.6f and its first point tore at %.6f (%.1f%% of"
                " it); the fibre's is %.6f and it tore at %.6f (%.1f%%)\n",
                live.plateFailureStrain, live.firstPlateTear,
                100.0 * live.firstPlateTear / live.plateFailureStrain, failureStrain,
                live.peakPlastic, 100.0 * live.peakPlastic / failureStrain);
    expectTrue("the two failure strains are close, so they are not the difference",
               live.plateFailureStrain < 1.2 * failureStrain &&
                   failureStrain < 1.2 * live.plateFailureStrain);
    expectTrue("but the plating tears at a small fraction of its own, where the fibre tears at"
               " all of its -- which is Rice-Tracey and nothing else",
               live.firstPlateTear < 0.2 * live.plateFailureStrain);

    // **The un-conservative direction, at the zone.** The control drives the
    // longitudinal to a plastic strain multiples of its own failure strain and still
    // reports it whole.
    std::printf("     the control leaves the longitudinal at eps_p = %.4f, %.2fx its own"
                " regularised failure strain of %.4f, and still reports it at full section\n",
                immortal.peakPlastic, immortal.peakPlastic / failureStrain, failureStrain);
    expectTrue("the control carries a longitudinal well past the strain that fails it",
               immortal.peakPlastic > 2.0 * failureStrain);
    // **It caps it to one increment, not to the bit, and that is the criterion
    // rather than a defect.** Damage is accumulated and tested at the end of a step,
    // so the step that crosses one leaves `eps_p` past `eps_f` by whatever that step
    // added -- exactly as `plasticity::update` does for an element. What *is* exact
    // is the damage variable, which is why the assertion below is on that.
    std::printf("     the criterion stops it at eps_p = %.6f against eps_f = %.6f: one step's"
                " increment past, %.3f%% -- damage is accumulated and tested at the end of a"
                " step, as it is for an element\n", live.peakPlastic, failureStrain,
                100.0 * (live.peakPlastic / failureStrain - 1.0));
    expectTrue("the criterion stops it within one increment of its failure strain",
               live.peakPlastic >= failureStrain && live.peakPlastic < 1.01 * failureStrain);
    expectTrue("and the damage it accumulated reached one exactly", live.peakDamage == 1.0);
    testing::expectEqual("while the control tore nothing", static_cast<long long>(immortal.tornFibers), 0LL);

    std::printf("     so the un-conservative model claims %.2fx the fibre force, %.2fx the stored"
                " energy, and %.1f%% more punch work to open the same hole\n",
                immortal.bundleForce / live.bundleForce, immortal.fibreEnergy / live.fibreEnergy,
                100.0 * (immortal.work / live.work - 1.0));
    expectTrue("the control has the longitudinal carrying substantially more force",
               immortal.bundleForce > 1.5 * live.bundleForce);
    expectTrue("and storing substantially more energy",
               immortal.fibreEnergy > 3.0 * live.fibreEnergy);
    expectTrue("and it costs the ram more energy to reach the same penetration, which is the"
               " un-conservative direction on the collision as well as on the section",
               immortal.work > 1.1 * live.work);

    // **And it reaches `reduce()`, which is the thing the roadmap said was missing.**
    expectTrue("the reduction names the longitudinal", live.members == 1);
    testing::expectEqual("and the control names none", static_cast<long long>(immortal.members), 0LL);
    expectTrue("its section is reduced", live.memberArea < whole);
    expectTrue("and the control's is untouched, bit for bit", immortal.memberArea == whole);
    // **The closed form for what is left, and it goes through the unmeshed part.**
    // The member is `whole = A * L_member` of steel; the zone tore `tornFiberVolume`
    // of it and looked at the rest of what it meshed; everything outside the zone is
    // intact because nothing looked at it. So the section left is
    // `1 - tornVolume / whole`, and the zone's radius does not appear.
    const double wholeVolume = whole * memberLength;
    const double predicted = 1.0 - live.tornVolume / wholeVolume;
    std::printf("     the zone meshed %.1f%% of a %.2f m longitudinal and tore %.5e m3 of its"
                " %.5e m3; that leaves %.6f of its section, and the reduction says %.6f\n",
                100.0 * live.meshedFraction, memberLength, live.tornVolume, wholeVolume,
                predicted, live.memberArea / whole);
    // Vacuity for the rule this fixture exists to test: the zone must *not* have
    // meshed the whole member, or "the part nobody looked at is intact" is untested.
    expectTrue("the zone looked at part of the longitudinal and not all of it",
               live.meshedFraction > 0.1 && live.meshedFraction < 0.9);
    expectNear("what is left is the untorn share of the whole member, unmeshed part included",
               live.memberArea / whole, predicted, 1e-12);
    // Two independent routes to the volume that went: the solver's own running total
    // over the fibres, and the reduction's roll-up per member. They are computed in
    // different files from the same state and must agree exactly.
    expectNear("the solver and the reduction agree on how much steel went", live.tornVolume,
               live.lostVolume, 0.0);
    expectTrue("and it is a real volume, not a count", live.tornVolume > 0);
    expectTrue("the plating went too, so this is the case the roadmap describes",
               live.panels > 0);

    // **The lost mass splits, and the stiffener's half is its lost steel.** `1 - eff`
    // is `lost / whole` by construction, so `(1 - eff) * whole * density` is exactly
    // the torn volume times the density -- a closed form for a number that would
    // otherwise be a plausible-looking total.
    const double memberMass = live.lostSteelMass - live.lostPlateMass;
    const double wantMemberMass = live.tornVolume * ah36Steel().density;
    std::printf("     the reduction says %.4f kg of steel is gone, of which the plating is %.4f"
                " and the longitudinal %.4f (its torn volume x density is %.4f)\n",
                live.lostSteelMass, live.lostPlateMass, memberMass, wantMemberMass);
    expectNear("the stiffener's share of the lost mass is its torn volume times its density",
               memberMass, wantMemberMass, 1e-9 * wantMemberMass);
    expectTrue("and it is a real share, not a rounding residue", memberMass > 0.5);
    // The control has no member damage at all, so its total is the plating's alone.
    expectNear("the control's lost mass is the plating's and nothing else",
               immortal.lostSteelMass, immortal.lostPlateMass, 1e-9 * immortal.lostSteelMass);

    // Anonymous fibres: named nothing, reported loudly.
    std::printf("     read through a patch whose fibres carry no member index the same solve"
                " names %zu members and files %zu problems\n", live.anonymousMembers,
                live.anonymousProblems);
    testing::expectEqual("fibres with no member index reduce nothing",
                         static_cast<long long>(live.anonymousMembers), 0LL);
    expectTrue("but the reduction says so rather than reporting healthy stiffeners",
               live.anonymousProblems > 0);
    expectTrue("and it names the right fault -- fibres with no member index, not a member the"
               " structure does not have", live.anonymousDiagnosed);
    // Vacuity: the named read of the same solve *did* find a damaged member, so the
    // zero above is the missing label and not a zone that tore nothing.
    testing::expectEqual("while the same solve read with its labels names one",
                         static_cast<long long>(live.members), 1LL);
    testing::expectEqual("a member index the structure does not have reduces nothing either",
                         static_cast<long long>(live.strangerMembers), 0LL);
    expectTrue("and is diagnosed as that rather than indexed out of bounds",
               live.strangerDiagnosed);
}

// --- 10. The reaction back: conservative, in the right direction ----------------

void testTheReductionMakesTierZeroReportLessStrength() {
    std::printf("\n   a damaged zone makes Tier-0 report less strength\n");
    const StructuralMesh& structure = ferryStructure();
    const Scantlings scantlings = ferryScantlings();

    // Every side panel within 6 m of the sheer at midship, swept from intact to
    // gone. Driven synthetically rather than by a solve, because the coupling is
    // what is under test and a solve would only make it slower and less severe.
    promotion::SectionReduction reduction;
    for (std::size_t i = 0; i < structure.panels.size(); ++i) {
        const PlatePanel& p = structure.panels[i];
        if (p.role != PanelRole::Shell) continue;
        if (length(p.centroid() - Vec3{0.0, -9.0, 13.0}) > 6.0) continue;
        promotion::PanelDamage damage;
        damage.panel = static_cast<int>(i);
        reduction.panels.push_back(damage);
    }
    expectTrue("the sweep covers a real patch of her side", reduction.panels.size() > 20);

    const auto measure = [&](double effectiveness) {
        for (promotion::PanelDamage& d : reduction.panels) d.effectiveness = effectiveness;
        const StructuralMesh damaged = promotion::reduce(structure, reduction);
        const HullGirderSection section = hullGirderSection(damaged, 0.0);
        const std::vector<CollapseElement> elements = collapseElementsAt(damaged, scantlings, 0.0);
        struct Row {
            double area, second, modulusDeck, modulusKeel, hog, sag;
        };
        return Row{section.area, section.secondMoment, section.modulusDeck, section.modulusKeel,
                   collapseCurve(elements, 1.0).ultimateMoment,
                   collapseCurve(elements, -1.0).ultimateMoment};
    };

    const auto intact = measure(1.0);
    std::printf("     %8s %10s %10s %10s %10s %12s %12s\n", "eff", "area", "I", "Zdeck", "Zkeel",
                "M_ult hog", "M_ult sag");
    auto previous = intact;
    double worstSagRise = 0;
    for (double effectiveness : {1.0, 0.875, 0.75, 0.625, 0.5, 0.375, 0.25, 0.125, 0.0}) {
        const auto row = measure(effectiveness);
        std::printf("     %8.2f %10.5f %10.4f %10.5f %10.5f %12.5e %12.5e\n", effectiveness,
                    row.area, row.second, row.modulusDeck, row.modulusKeel, row.hog, row.sag);
        if (effectiveness < 1.0) {
            // The exact section quantities, which are linear in the thickness and
            // so must be strictly monotone.
            expectTrue("section area falls with damage", row.area < previous.area);
            expectTrue("second moment falls with damage", row.second < previous.second);
            expectTrue("and the modulus at the damaged fibre with it",
                       row.modulusDeck < previous.modulusDeck);
            // Strength, in both directions, against the *intact* ship. That is the
            // claim -- never more than intact -- and it is a stronger statement
            // than step-to-step monotonicity, which the ultimate moment does not
            // quite have (see below).
            expectTrue("the hogging ultimate moment never exceeds the intact one",
                       row.hog < intact.hog);
            expectTrue("nor the sagging one", std::abs(row.sag) < std::abs(intact.sag));
            expectTrue("and hogging falls monotonically", row.hog < previous.hog);
            worstSagRise = std::max(worstSagRise,
                                    std::abs(row.sag) / std::abs(previous.sag) - 1.0);
        }
        previous = row;
    }
    std::printf("     losing that plating costs %.1f%% of the hogging ultimate moment and"
                " %.1f%% of the sagging one\n", 100.0 * (1.0 - previous.hog / intact.hog),
                100.0 * (1.0 - std::abs(previous.sag / intact.sag)));
    expectTrue("and the sweep is not a rounding error", previous.hog < 0.97 * intact.hog);

    // **Two things that are not monotone, measured rather than assumed.**
    //
    // The section modulus at the *undamaged* fibre can rise: taking material away
    // above the neutral axis pulls the axis down, and the keel's lever arm shrinks
    // faster than the second moment does. So `modulusKeel` is not a conservative
    // reading of damage above it, and the ultimate moment is.
    promotion::SectionReduction nearAxis;
    for (std::size_t i = 0; i < structure.panels.size(); ++i) {
        const PlatePanel& p = structure.panels[i];
        if (p.role != PanelRole::Shell) continue;
        if (length(p.centroid() - Vec3{0.0, -10.0, 8.0}) > 2.0) continue;
        promotion::PanelDamage damage;
        damage.panel = static_cast<int>(i);
        damage.effectiveness = 0.0;
        nearAxis.panels.push_back(damage);
    }
    const HullGirderSection before = hullGirderSection(structure, 0.0);
    const HullGirderSection after = hullGirderSection(promotion::reduce(structure, nearAxis), 0.0);
    std::printf("     damage just above the neutral axis moves it %.4f -> %.4f m and takes the"
                " keel modulus %.5f -> %.5f m3 -- *upwards*\n", before.neutralAxis,
                after.neutralAxis, before.modulusKeel, after.modulusKeel);
    expectTrue("the neutral axis moves towards the damage-free side",
               after.neutralAxis < before.neutralAxis);
    expectTrue("and the far fibre's section modulus rises, so it is not the conservative reading",
               after.modulusKeel > before.modulusKeel);
    expectTrue("while the near one falls", after.modulusDeck < before.modulusDeck);

    // The other is the sagging ultimate moment, which does not fall at every step:
    // a thinner panel buckles *earlier* and therefore sheds earlier, and Smith's
    // method lets the neutral axis migrate in response, so a little more damage can
    // leave a little more moment. It is in the model rather than in the quadrature
    // -- eightfold the curvature steps leaves it where it is -- and it is why the
    // conservatism claim above is against the *intact* section rather than
    // step-to-step.
    std::printf("     the sagging branch's worst step-to-step rise with damage is %+.3f%%:"
                " load shedding is not monotone in thickness\n", 100.0 * worstSagRise);
    expectTrue("and where the sagging branch rises it rises by well under a per cent",
               worstSagRise < 0.01);
}

// --- 11. What a solved zone reports ---------------------------------------------

void testASolvedZoneReportsTheThicknessItHasLeft() {
    std::printf("\n   what a solved zone hands back\n");
    // Subdivision two, the fixture `test_zone.cpp` uses for tearing, and driven
    // deep enough that a bay goes through: a reduction that never reaches zero
    // would never exercise the case where the two answers -- what `breach.hpp`
    // opens and what the section loses -- have to agree.
    const StructuralMesh strip = flatStrip(0.6, 0.3, 0.006, 3, 1);
    const zone::Patch patch = zone::buildPatch(strip, {0.05, 0, 0}, flatParams(2));

    zone::SolveParams solve;
    solve.indenter.halfLength = 0.03;
    solve.indenter.halfWidth = 0.08;
    solve.indenter.speed = 20.0;
    solve.indenter.rampTime = 5.0e-4;
    solve.indenter.stopAt = 0.09;
    zone::Solver solver(patch, plasticity::shipSteel(), solve);
    const zone::SolveResult& result = solver.run();
    expectTrue("something tore", result.tornElements > 0);

    const promotion::SectionReduction reduction =
        promotion::reactionOf(strip, patch, solver);
    std::printf("     %d of %zu elements torn; %zu panel(s) reduced, worst effectiveness %.4f,"
                " worst dent %.4f m, %.1f kg of steel gone\n", result.tornElements,
                patch.elementCount(), reduction.panels.size(), reduction.worstEffectiveness,
                reduction.worstOutOfPlane, reduction.lostSteelMass);
    for (const promotion::PanelDamage& damage : reduction.panels)
        std::printf("       panel %d: effectiveness %.4f, %.4f m2 torn, %.4f m of dent,"
                    " %.5f m thinner\n", damage.panel, damage.effectiveness, damage.tornArea,
                    damage.outOfPlane, damage.thinning);

    expectTrue("the zone reports damage", !reduction.panels.empty());
    expectTrue("every effectiveness is a fraction",
               std::all_of(reduction.panels.begin(), reduction.panels.end(),
                           [](const promotion::PanelDamage& d) {
                               return d.effectiveness >= 0.0 && d.effectiveness < 1.0;
                           }));

    // The panels `breach.hpp` will open are the ones the reduction has written
    // down to nothing or nearly so -- the two answers come from the same element
    // states and must not disagree about which bay is gone.
    for (int torn : result.tornPanels) {
        double effectiveness = 1.0;
        for (const promotion::PanelDamage& damage : reduction.panels)
            if (damage.panel == torn) effectiveness = damage.effectiveness;
        expectTrue("a panel breach.hpp will open has lost at least half its section",
                   effectiveness <= 0.5);
    }
    // Guard: the run has to have torn *some* panels through, or the loop above
    // asserted nothing at all.
    expectTrue("and there were panels to check", !result.tornPanels.empty());

    // A panel nothing happened to is absent rather than present at 1.0, so a
    // caller folding the reduction into a mesh touches only what was damaged.
    expectTrue("only damaged panels are reported",
               reduction.panels.size() <= patch.panels.size());

    // **A pre-load is not damage, and neither is numerical dust.** A patch under a
    // real pre-load and no punch at all, stepped forward far enough that its nodes
    // have actually moved, must report nothing: the elements shift by nanometres
    // and an unfloored reduction names every panel it can see. Measured on the
    // ferry that doubled the damaged-panel count under her own hogging stress,
    // every extra one intact to six figures.
    zone::SolveParams stressed;
    stressed.indenter.halfLength = -1;
    stressed.indenter.halfWidth = -1;
    stressed.preload.stress = 0.5 * 355.0e6;
    stressed.duration = 2.0e-3;
    zone::Solver unpunched(patch, plasticity::shipSteel(), stressed);
    unpunched.run();
    const promotion::SectionReduction quiet = promotion::reactionOf(strip, patch, unpunched);
    std::printf("     the same patch pre-loaded to %.0f MPa and stepped %d times with no punch:"
                " %zu panel(s) reduced, worst node moved %.2e m\n",
                stressed.preload.stress / 1e6, unpunched.result().steps, quiet.panels.size(),
                unpunched.largestDisplacement());
    expectTrue("the pre-loaded patch really did carry a strain", 
               unpunched.largestDisplacement() > 1e-9);
    testing::expectEqual("but a pre-loaded patch nothing has struck has lost nothing",
                         static_cast<long long>(quiet.panels.size()), 0);
    expectNear("and its section is untouched", quiet.worstEffectiveness, 1.0, 0.0);

    // **And nor is curvature.** On a curved patch the elements are extruded along
    // nodal normals that disagree, so volume over mid-surface area is not quite the
    // thickness the strake was authored with. Measuring against the nominal would
    // report that geometry as section lost the moment the zone was promoted, before
    // anything had touched it. A cylinder, unsolved, must report nothing at all.
    StructuralMesh curved;
    curved.materials.push_back(ah36Steel());
    const double radius = 3.0, sweep = 1.0;
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 4; ++j) {
            const double x0 = -0.3 + 0.3 * i, x1 = x0 + 0.3;
            const double a0 = -0.5 * sweep + sweep * j / 4.0, a1 = a0 + sweep / 4.0;
            PlatePanel p;
            p.corner[0] = {x0, radius * std::sin(a0), radius * std::cos(a0)};
            p.corner[1] = {x1, radius * std::sin(a0), radius * std::cos(a0)};
            p.corner[2] = {x1, radius * std::sin(a1), radius * std::cos(a1)};
            p.corner[3] = {x0, radius * std::sin(a1), radius * std::cos(a1)};
            p.thickness = 0.020;
            p.role = PanelRole::Shell;
            curved.panels.push_back(p);
        }
    zone::MeshParams curvedMesh;
    curvedMesh.radius = 1e3;
    curvedMesh.subdivision = 2;
    curvedMesh.outward = {0, 0, 1};
    const zone::Patch bilge = zone::buildPatch(curved, {0.0, 0.0, radius}, curvedMesh);
    expectTrue("the curved patch meshed and is genuinely curved",
               !bilge.empty() && bilge.worstNormalSpread > 0.01);
    zone::SolveParams untouched;
    untouched.indenter.halfLength = -1;
    untouched.indenter.halfWidth = -1;
    zone::Solver never(bilge, plasticity::shipSteel(), untouched);
    const promotion::SectionReduction pristine = promotion::reactionOf(curved, bilge, never);
    std::printf("     a curved patch (normals spread %.4f rad), unsolved: %zu panel(s) reduced\n",
                bilge.worstNormalSpread, pristine.panels.size());
    testing::expectEqual("a patch nothing has happened to has lost nothing, curved or not",
                         static_cast<long long>(pristine.panels.size()), 0);

    // **The invariant the "part nobody looked at" term rests on.** The mesher takes
    // whole panels or none, so every panel in a patch is meshed edge to edge and
    // the unmeshed share is zero. That is why a fully torn panel comes out at
    // exactly zero rather than at an epsilon -- and it is an invariant of the
    // *mesher*, not of this file, so it is asserted here rather than assumed.
    for (const promotion::PanelDamage& damage : reduction.panels)
        expectNear("every panel in a patch is meshed whole", damage.meshedFraction, 1.0, 1e-9);

    // **The epsilon that was not tidiness.** A completely torn panel has to come
    // out at exactly zero: a panel meshed edge to edge has its meshed area equal
    // to its own area only to rounding, and an effectiveness of 1e-16 leaves a
    // plate 1e-18 m thick, which is not zero, so `collapseElementsAt` keeps it and
    // a plate that thin has a critical stress of 1e-16 Pa.
    for (const promotion::PanelDamage& damage : reduction.panels) {
        double meshedArea = 0, tornArea = 0;
        for (std::size_t e = 0; e < patch.elementCount(); ++e)
            if (patch.panelOf[e] == damage.panel) {
                meshedArea += patch.elementArea[e];
                if (solver.elementState()[e].torn) tornArea += patch.elementArea[e];
            }
        if (tornArea < 0.999 * meshedArea) continue;
        expectTrue("a panel torn through reports exactly zero, not an epsilon",
                   damage.effectiveness == 0.0);
    }

    // **A zone that has *not* torn still has something to hand back**, and it is a
    // different pair of quantities: the plating is dented out of plane and it is
    // thinner where it stretched. Both are measured from the deformed geometry
    // rather than inferred from a plastic strain and an assumed Poisson ratio.
    //
    // A separate, finer, shallower run, because the dent is measured on what
    // *survives* -- a deleted element's node positions mean nothing -- so a zone
    // torn through reports almost no dent, which the coarse run above shows.
    const zone::Patch fine = zone::buildPatch(strip, {0.05, 0, 0}, flatParams(4));
    zone::SolveParams dent = solve;
    dent.indenter.stopAt = 0.020;
    zone::Solver denting(fine, plasticity::shipSteel(), dent);
    const zone::SolveResult& dented = denting.run();
    const promotion::SectionReduction bent = promotion::reactionOf(strip, fine, denting);
    std::printf("     driven only to %.3f m: %d torn, worst effectiveness %.4f, dent %.4f m,"
                " thinned by %.5f m of %.4f\n", dented.penetration, dented.tornElements,
                bent.worstEffectiveness, bent.worstOutOfPlane,
                bent.panels.empty() ? 0.0 : bent.panels.front().thinning, fine.thickness);
    testing::expectEqual("nothing tore in the shallow run", dented.tornElements, 0);
    // The punch grips the nodes it touches, so their travel *along the patch
    // normal* is exactly the penetration -- and the dent is the component along
    // that normal, never the length of a displacement that also has the membrane
    // pull-in in it. Bracketed both ways, since one side alone would pass against
    // an unprojected magnitude.
    expectNear("the dent is the punch's own travel along the patch normal",
               bent.worstOutOfPlane, dented.penetration, 1e-6);
    expectTrue("and there is in-plane motion for the projection to have removed",
               denting.largestDisplacement() > dented.penetration * (1.0 + 1e-4));
    expectTrue("and it is thinner for having stretched", !bent.panels.empty() &&
                                                             bent.worstEffectiveness < 1.0);
    expectTrue("by a fraction of its thickness, not by all of it",
               bent.worstEffectiveness > 0.9);
}

// --- 12. The Tier-0 defect this coupling found ----------------------------------
//
// Not a promotion test so much as the regression for what promotion exposed:
// `longitudinalStrength` sized its collapse sweep from `firstYieldCurvature`,
// which is set by the *weakest* element in the section. A damaged bay buckles at a
// couple of MPa, which puts first yield five orders below the curvature at which
// the section actually collapses, and the sweep then never reached the peak.

void testACollapseSweepReachesThePeakOnADamagedSection() {
    std::printf("\n   a damaged section still has to reach its own peak\n");
    const StructuralMesh& structure = ferryStructure();
    const Scantlings scantlings = ferryScantlings();

    promotion::SectionReduction reduction;
    for (std::size_t i = 0; i < structure.panels.size(); ++i) {
        const PlatePanel& p = structure.panels[i];
        if (p.role != PanelRole::Shell) continue;
        if (length(p.centroid() - Vec3{0.0, -9.0, 13.0}) > 6.0) continue;
        promotion::PanelDamage damage;
        damage.panel = static_cast<int>(i);
        damage.effectiveness = 0.125;   // an eighth of its thickness left
        reduction.panels.push_back(damage);
    }
    const StructuralMesh eighth = promotion::reduce(structure, reduction);
    // The same panels taken away *entirely*, which is unambiguously worse damage.
    for (promotion::PanelDamage& damage : reduction.panels) damage.effectiveness = 0.0;
    const StructuralMesh gone = promotion::reduce(structure, reduction);

    const std::vector<CollapseElement> intact = collapseElementsAt(structure, scantlings, 0.0);
    const std::vector<CollapseElement> hurt = collapseElementsAt(eighth, scantlings, 0.0);
    const std::vector<CollapseElement> removed = collapseElementsAt(gone, scantlings, 0.0);

    const double intactYield = firstYieldCurvature(intact);
    const double hurtYield = firstYieldCurvature(hurt);
    const double robust = extremeFibreYieldCurvature(hurt);
    std::printf("     first yield curvature: intact %.4e, at an eighth thickness %.4e"
                " (%.0fx smaller); the extreme fibre's own is %.4e\n", intactYield, hurtYield,
                intactYield / hurtYield, robust);
    expectTrue("thinning one strake drops the section's first yield by more than an order",
               hurtYield < 0.1 * intactYield);
    expectTrue("while the extreme fibre's yield curvature barely moves",
               robust > 0.5 * extremeFibreYieldCurvature(intact));

    // The old sizing, reproduced: a sweep to six times first yield, which is where
    // the peak of an *undamaged* section lives.
    const auto oldSizing = [](const std::vector<CollapseElement>& e) {
        return progressiveCollapse(e, 6.0 * firstYieldCurvature(e), 150).ultimateMoment;
    };
    const CollapseCurve full = collapseCurve(hurt, 1.0);
    std::printf("     ultimate moment at an eighth thickness: %.4e swept from first yield,"
                " %.4e swept until the peak is inside; intact is %.4e\n", oldSizing(hurt),
                full.ultimateMoment, collapseCurve(intact, 1.0).ultimateMoment);
    expectTrue("the sweep sized from first yield never reaches the peak",
               oldSizing(hurt) < 0.2 * full.ultimateMoment);
    expectTrue("the peak of the corrected sweep is strictly inside it",
               std::abs(full.ultimateCurvature) < 0.999 * std::abs(full.points.back().curvature));

    // **The signature, and the reason this is a defect rather than a tolerance.**
    // Under the old sizing, taking the plating away *entirely* reports more
    // strength than leaving an eighth of it -- a hull girder that gets stronger
    // when material is removed. Under the corrected one it reports less, which is
    // the only orderable answer there is.
    std::printf("     removing that plating entirely: %.4e under the old sizing against %.4e with"
                " an eighth of it left, and %.4e once the sweep reaches the peak\n",
                oldSizing(removed), oldSizing(hurt), collapseCurve(removed, 1.0).ultimateMoment);
    expectTrue("the old sizing reported a section as stronger with less material in it",
               oldSizing(removed) > 2.0 * oldSizing(hurt));
    expectTrue("and the corrected one puts them back in order",
               collapseCurve(removed, 1.0).ultimateMoment < full.ultimateMoment);

    // The fix must not be a way of getting the strength back.
    expectTrue("a damaged section is still weaker than an intact one",
               full.ultimateMoment < collapseCurve(intact, 1.0).ultimateMoment);
    expectTrue("but not by an order of magnitude",
               full.ultimateMoment > 0.5 * collapseCurve(intact, 1.0).ultimateMoment);

    // An intact section's sweep is unchanged to the last bit, which is what makes
    // this a fix rather than a re-tuning.
    const CollapseCurve before = progressiveCollapse(intact, 6.0 * intactYield, 150);
    expectTrue("and an intact section's answer is bit-identical to the old sizing",
               collapseCurve(intact, 1.0).ultimateMoment == before.ultimateMoment);
}

// --- 13. The small closed forms the criterion is built out of -------------------

void testTheContactTriggerIsAClosedForm() {
    std::printf("\n   the contact trigger\n");
    // The span of a bay is its *shorter* side, and taking it from the panel itself
    // is the one definition that cannot be got the wrong way round -- which is the
    // defect `indentation.hpp` records, where the plating was given the frame
    // spacing to span when it spans the longitudinals.
    PlatePanel bay;
    bay.corner[0] = {-1.2, 0.0, 0.0};
    bay.corner[1] = {1.2, 0.0, 0.0};
    bay.corner[2] = {1.2, 0.7, 0.0};
    bay.corner[3] = {-1.2, 0.7, 0.0};
    expectNear("a 2.4 by 0.7 bay spans 0.7", promotion::panelSpan(bay), 0.70, 1e-12);

    // A tapered panel reports the mean of the pair, not the narrow end: half of a
    // wedge is not a span.
    PlatePanel tapered = bay;
    tapered.corner[2] = {1.2, 0.5, 0.0};
    expectNear("and a tapered one its mean side", promotion::panelSpan(tapered), 0.60, 1e-12);

    // A wedge, where the *averaging* is what decides. One pair runs 0.3 and 2.4 --
    // mean 1.35 -- against 0.7 and 2.214, mean 1.457. So the span is 1.35, and a
    // routine that took one edge of each pair rather than the mean would answer
    // 0.3: a quarter of it, and a bay sixteen times too strong.
    PlatePanel wedge;
    wedge.corner[0] = {-1.2, 0.0, 0.0};
    wedge.corner[1] = {-0.9, 0.0, 0.0};
    wedge.corner[2] = {1.2, 0.7, 0.0};
    wedge.corner[3] = {-1.2, 0.7, 0.0};
    const double narrowPair = 0.5 * (0.3 + 2.4);
    const double widePair = 0.5 * (0.7 + std::sqrt(2.1 * 2.1 + 0.7 * 0.7));
    expectTrue("the wedge's two pairs are the way round the test needs", narrowPair < widePair);
    expectNear("a wedge spans the mean of its narrower pair", promotion::panelSpan(wedge),
               narrowPair, 1e-12);
    expectTrue("and the mean is nowhere near either of the edges it averages",
               narrowPair > 4.0 * 0.3 && narrowPair < 0.6 * 2.4);

    // 4 sigma_y (t/L)^2 is the three-hinge mechanism of a clamped strip:
    // w L^2 / 16 = M_p with M_p = sigma_y t^2 / 4.
    const double yieldStrength = 355.0e6, thickness = 0.012, span = 0.70;
    const double plastic = yieldStrength * thickness * thickness / 4.0;
    expectNear("the plating's collapse pressure is 16 M_p / L^2",
               promotion::platingCollapsePressure(yieldStrength, thickness, span),
               16.0 * plastic / (span * span), 1e-6 * 16.0 * plastic / (span * span));
    std::printf("     12 mm plating on a 0.70 m span hinges at %.0f kPa\n",
                promotion::platingCollapsePressure(yieldStrength, thickness, span) / 1e3);

    // Twice the thickness is four times the pressure, which is the (t/L)^2 the
    // whole scantling trade lives on.
    expectNear("and it goes as the square of the thickness",
               promotion::platingCollapsePressure(yieldStrength, 2.0 * thickness, span) /
                   promotion::platingCollapsePressure(yieldStrength, thickness, span),
               4.0, 1e-9);
}

void testTheDentedCapacityKnockdownIsContinuousAndMonotone() {
    std::printf("\n   what a dent costs a panel in compression\n");
    const double yieldStrength = 355.0e6, critical = 200.0e6, thickness = 0.012;
    expectNear("an undented panel keeps min(yield, critical)",
               promotion::dentedCompressiveCapacity(yieldStrength, critical, 0.0, thickness),
               critical, 1e-9 * critical);
    expectNear("and the other way round when yield comes first",
               promotion::dentedCompressiveCapacity(critical, yieldStrength, 0.0, thickness),
               critical, 1e-9 * critical);

    double previous = promotion::dentedCompressiveCapacity(yieldStrength, critical, 0.0, thickness);
    for (double deviation : {0.001, 0.003, 0.010, 0.030, 0.100}) {
        const double capacity =
            promotion::dentedCompressiveCapacity(yieldStrength, critical, deviation, thickness);
        std::printf("     dent %.3f m (%.1f t): %.1f MPa of %.1f\n", deviation,
                    deviation / thickness, capacity / 1e6, critical / 1e6);
        expectTrue("a deeper dent is always weaker", capacity < previous);
        previous = capacity;
    }
    // Independently: the root of `sigma (1 + eta/(1 - sigma/sigma_cr)) = sigma_y`,
    // found by bisection rather than by the quadratic the implementation solves.
    const double deviation = 0.010;
    const double eta = 6.0 * deviation / thickness;
    double lo = 0.0, hi = std::min(yieldStrength, critical);
    for (int i = 0; i < 200; ++i) {
        const double mid = 0.5 * (lo + hi);
        const double demand = mid * (1.0 + eta / (1.0 - mid / critical));
        (demand < yieldStrength ? lo : hi) = mid;
    }
    expectNear("and it is the root of the Perry-Robertson equation",
               promotion::dentedCompressiveCapacity(yieldStrength, critical, deviation, thickness),
               0.5 * (lo + hi), 1e-6 * critical);
}

void testTheElementEstimateOverStatesRatherThanUnder() {
    std::printf("\n   the element estimate the budget is spent against\n");
    // On flat plating the flood fill reaches everything, so the estimate is exact.
    const StructuralMesh strip = flatStrip(1.2, 0.8, 0.012, 3, 2);
    zone::MeshParams mesh = flatParams(3);
    const zone::Patch flat = zone::buildPatch(strip, {0, 0, 0}, mesh);
    testing::expectEqual("on flat plating the estimate is exact",
                         promotion::estimateElements(strip, {0, 0, 0}, mesh),
                         static_cast<long long>(flat.elementCount()));

    // It counts the role the zone will be built from and no other. A deck landing
    // on top of a shell strip is the ordinary case at a deck edge, and counting
    // both would charge a shell zone for plating it will never mesh.
    StructuralMesh mixed = strip;
    for (const PlatePanel& panel : strip.panels) {
        PlatePanel deck = panel;
        deck.role = PanelRole::Deck;
        mixed.panels.push_back(deck);
    }
    zone::MeshParams shellOnly = mesh;
    shellOnly.role = PanelRole::Shell;
    std::printf("     %zu shell panels under %zu deck panels: the estimate counts %d\n",
                strip.panels.size(), strip.panels.size(),
                promotion::estimateElements(mixed, {0, 0, 0}, shellOnly));
    testing::expectEqual("only the role the zone will be meshed from is counted",
                         promotion::estimateElements(mixed, {0, 0, 0}, shellOnly),
                         promotion::estimateElements(strip, {0, 0, 0}, shellOnly));

    // On a real ship a fold or a thickness seam truncates the patch, so the
    // estimate over-states -- which is the safe direction for a budget, and means
    // a promoted zone can be smaller than it was charged for.
    const StructuralMesh& structure = ferryStructure();
    zone::MeshParams shipMesh;
    shipMesh.radius = 3.0;
    shipMesh.subdivision = 4;
    for (const Vec3& impact : {Vec3{0.0, -9.9, 8.0}, Vec3{0.0, 0.0, 0.0}, Vec3{30.0, -9.0, 5.0}}) {
        const zone::Patch patch = zone::buildPatch(structure, impact, shipMesh);
        const int estimate = promotion::estimateElements(structure, impact, shipMesh);
        std::printf("     at (%.0f, %.0f, %.0f): estimated %d elements, meshed %zu\n", impact.x,
                    impact.y, impact.z, estimate, patch.elementCount());
        expectTrue("the estimate is never below what is meshed",
                   estimate >= static_cast<int>(patch.elementCount()));
        expectTrue("and never wildly above it",
                   estimate <= 4 * static_cast<int>(patch.elementCount()) + 16);
    }
}

}  // namespace

void runPromotionTests() {
    std::printf("\n--- adaptive zone promotion ---\n");
    testNothingPromotesOnAShipAtRestInStillWater();
    testAUniformlyLoadedShipDoesNotPromoteEverywhere();
    testTheBackgroundIsAMedianAndTheThresholdIsWhereItSays();
    testTheZoneLandsOnTheFibreThatIsInTrouble();
    testCandidatesAreRankedByHowFarPastTheirOwnThresholdTheyAre();
    testPromotionIsIdempotentAndDeterministic();
    testHysteresisAndDwellPreventChatter();
    testCostIsLinearInTheNumberOfZones();
    testTheDecisionIsCheapAndTheTierZeroAnswerIsNot();
    testAPreLoadedPatchCarriesTheStressItWasHanded();
    testAPreLoadSpendsYieldCapacityExactly();
    testAPreLoadedZoneYieldsEarlierUnderAPunch();
    testThePreLoadIsTheGirdersOwnStress();
    testTheLongitudinalsARamDestroysAreWorthAFifthOfTheGirder();
    testASolvedZoneNamesTheLongitudinalItTore();
    testTheReductionMakesTierZeroReportLessStrength();
    testASolvedZoneReportsTheThicknessItHasLeft();
    testACollapseSweepReachesThePeakOnADamagedSection();
    testTheContactTriggerIsAClosedForm();
    testTheDentedCapacityKnockdownIsContinuousAndMonotone();
    testTheElementEstimateOverStatesRatherThanUnder();
}
