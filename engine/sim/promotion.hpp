// SPDX-License-Identifier: MIT
//
// Adaptive zone promotion: the piece that makes the structural tiers one system
// rather than three models that do not speak.
//
// `girder.hpp`, `buckling.hpp` and `collapse.hpp` are Tier-0 -- the ship as a
// beam, everywhere, always, for microseconds. `zone.hpp` is Tier-2 -- solid-shell
// elements with plasticity over a few dozen square metres of real plating, for
// core-minutes. Until this file existed nothing decided **when** to build a zone,
// nothing told the beam what a torn zone had cost it, and nothing told the zone
// what the beam was already putting through it. Those are three separate holes
// and this file is three separate answers:
//
//   1. `Promoter` -- which patches deserve Tier-2, under a cost budget, without
//      chattering.
//   2. `reactionOf` / `reduce` -- what a solved zone hands back, as a structure
//      Tier-0 already knows how to read.
//   3. `preloadFor` -- the girder stress the zone starts from, as
//      `zone::Preload`.
//
// Where the three-tier plan in `02-simulation.md` §3 says a zone couples to Tier-1
// through retained interface DOF, this couples it to Tier-0 through a section:
// cruder, and the honest thing available. `reduction.{hpp,cpp}` now builds the
// Craig-Bampton reduced models themselves -- boundary DOF kept exactly, so they are
// couplable by construction -- but nothing yet drives a zone's edge from one, and
// **this file is unchanged by their arrival**. Everything below was written so that
// inserting Tier-1 replaces the coupling and not the criterion, and that is still
// the work outstanding.
//
// --- 1. Cost is the whole design ----------------------------------------------
//
// A Tier-2 element costs `4.0 x elementCount` core-seconds per simulated second
// (`zone.hpp` §1), and that is **linear in the number of elements and therefore
// linear in the number of zones**. There is no economy of scale to hide behind:
// promoting ten patches costs ten times one patch, exactly, and
// `tests/test_promotion.cpp` measures that rather than assuming it.
//
// The consequence shapes the criterion. A ship at 0.9 of her buckling capacity
// *everywhere* must not promote everywhere -- that is a ship-shaped Tier-2 model
// and it is unaffordable by three orders of magnitude. But a uniformly utilised
// ship is not a badly designed one; it is a **well** designed one, since making
// utilisation uniform is what scantling design is for. So high utilisation on its
// own cannot be the trigger.
//
// --- 2. The criterion, and what it misses -------------------------------------
//
// What Tier-2 adds over Tier-0 is an answer where the response is **local and
// non-linear**. Where the load is smooth and the structure is intact, a beam is
// right and a zone tells you nothing you did not already have for free. So a
// station qualifies only when it is *both*:
//
//   * past an absolute threshold for its trigger -- yield, plate/column buckling,
//     or the progressive-collapse margin; and
//   * standing at least `Criterion::localExcess` above the **median** over the
//     ship's own stations for that same trigger.
//
// The median is the flat-ship guard, and it is deliberately robust rather than a
// mean: a hull girder's utilisation curve is roughly parabolic, so the median sits
// near half the peak and a real peak clears it easily, while a *flat* profile at
// 0.9 has median 0.9 and excess zero and promotes nothing at all. That last case
// is not a failure to notice a problem. A ship uniformly past her buckling
// capacity has a **girder** answer, and `collapse.hpp` already gives it; twenty-four
// square metres of solid-shell elements would not add to it.
//
// A **contact patch** is the other trigger and it is not treated the same way. A
// beam cannot represent a local load at all -- it smears it over a station -- so
// there is no background to stand above, and the test is absolute: the mean
// contact pressure against the struck bay's own plastic collapse pressure,
//
//     p_L = 4 sigma_y (t / span)^2
//
// the closed form for a clamped strip hinging at its ends and its middle. `span`
// is taken as the **shorter side of the struck panel itself**, which is the one
// definition that cannot be got the wrong way round -- and getting it the wrong
// way round is precisely the defect `indentation.hpp` records in its own history,
// where the plating was given the frame spacing to span when it spans the
// longitudinals.
//
// Candidates are then ranked by `utilisation / promoteThreshold`, thinned so that
// no two zones sit closer than `Criterion::separation` zone radii, and accepted in
// order until the element budget is spent. The separation rule is what stops a
// broad moment peak from promoting eight overlapping zones down the length of a
// ship, and the budget is what makes the cost bounded rather than hoped for.
//
// **What it misses**, and none of these is small:
//
//  1. **Everything Tier-0 cannot see.** No torsion, no shear lag, no racking, no
//     transverse strength, no local pressure. A load that does not change `M(x)`
//     and is not handed in as a contact patch is invisible here, which includes
//     slamming, sloshing, a dropped weight, and a fire softening a bulkhead.
//  2. **It is one-dimensional.** A station says *where along the ship*, never
//     where around the girth. The site is taken as the panel nearest the extreme
//     fibre on the centreline, and a panel that is locally the weakest somewhere
//     else on the girth is not found. Nor can a beam say which **side** is in
//     trouble; only a contact patch can.
//  3. **A uniform overload promotes nothing**, by construction -- see above. That
//     is a decision, not an oversight, and it is the one most worth revisiting.
//  4. **It is a snapshot at a cadence.** `tierZero()` is 170 ms on the reference
//     ferry, of which 138 ms is the Smith's-method sweep, so it is emphatically
//     *not* run every tick -- see §3. A load that peaks and passes between two
//     reviews is missed entirely, and whipping and springing both live down there.
//     `Review::microseconds` measures the decision alone, which *is* tick-cheap;
//     what is not is the Tier-0 answer it reads.
//  5. **It has no memory.** Corrosion, fatigue cracking and previous damage enter
//     only once something has already changed the structural mesh -- which
//     `reduce()` below does, so damage is at least visible to the *next* review.
//  6. **The element count is estimated, not meshed.** `buildPatch` is 7 ms and the
//     criterion runs over every station, so the budget is spent against a panel
//     count times `subdivision^2`. That over-states wherever a fold or a thickness
//     seam truncates the patch, which is the safe direction for a budget and means
//     a promoted zone can be smaller than it was charged for.
//  7. **Nothing is evicted.** A zone holds its budget until its trigger clears. A
//     collision arriving while the budget is full is refused and reported, not
//     traded against a lower-scoring zone that is already running -- because
//     stopping a zone mid-solve throws away its plastic history.
//
// --- 3. Chatter, and the two mechanisms against it ----------------------------
//
// A criterion evaluated repeatedly on a load hovering at its threshold will
// promote and demote and promote again, and every promotion pays the meshing cost
// and throws away the previous zone's plastic state. Two mechanisms, because they
// catch different things:
//
//   * **Hysteresis.** Promote at `promote`, keep until below `hold`, with
//     `hold < promote`. This kills any oscillation whose amplitude fits inside the
//     band, whatever its frequency.
//   * **Dwell.** A candidate must qualify on `dwell` consecutive reviews before it
//     is promoted, and an active zone must fail to qualify on `hold`
//     consecutive reviews before it is dropped. This kills an oscillation *wider*
//     than the band, provided it is faster than the dwell.
//
// Neither subsumes the other and the tests exercise both against their own
// negative control -- the same signal with the mechanism switched off, which has
// to chatter, or the test proves nothing.
//
// --- 4. The pre-load ----------------------------------------------------------
//
// `preloadFor()` reads the bending moment and the section at the patch's station
// and returns `sigma_xx(z) = M (z - z_na) / I` as a `zone::Preload`. On the
// reference ferry poised on a 3 m crest that is 84 MPa at the deck against a
// 355 MPa yield, so a zone that starts unstressed starts with capacity it does not
// have.
//
// **How much capacity, and in which direction, is not the 24% the ratio
// suggests**, and measuring it is what showed that. Three separate measurements,
// because they say three different things:
//
//   * **The capacity it spends is exactly its magnitude.** With no punch at all a
//     patch yields when the pre-load reaches sigma_y and not before -- measured
//     between 0.99 and 1.01 of it, both signs -- so a zone handed the ferry's
//     84 MPa really does start with 76% of what an unloaded one claims. That half
//     of the argument holds, and `tests/test_promotion.cpp` asserts it against the
//     material's own yield stress.
//   * **Under a punch it moves yield onset, and only where the criterion would
//     have promoted anyway.** A patch pre-loaded to 0.9 sigma_y first yields at
//     0.059 of the penetration an unloaded one needs. At the ferry's 84 MPa the
//     shift is a few per cent -- and at 84 MPa nothing would have promoted a
//     girder-triggered zone, since her worst buckling utilisation on that crest is
//     0.26.
//   * **It barely moves tearing at all**, which is the correction. The pre-strain
//     is elastic -- 0.9 sigma_y is 1.6e-3 -- against a regularised failure strain
//     around 0.15, two orders larger; the first element lets go at 0.9942 of the
//     unloaded penetration. "A pre-loaded zone fails earlier" is true of yielding
//     and false of tearing.
//
// **And at ship scale it makes the zone stronger, not weaker**, which no argument
// from the ratio would have predicted. `tools/zone_probe --no-preload` exists so
// this is a measurement: a 2 m punch into the ferry's own side resists at
// **18.90 MN** at 0.078 m told it starts unstressed, and **20.25 MN** carrying her
// own 13.1 MPa of hogging tension.
//
// The reason is geometry. The girder's stress runs along x; the membrane stress a
// punch raises in longitudinally framed side plating runs across the bay, which is
// *vertical*, because the plating spans between longitudinals. The two are
// perpendicular, and von Mises does not add perpendicular stresses, it subtracts
// their product:
//
//     sigma_vm^2 = sigma_x^2 - sigma_x sigma_z + sigma_z^2
//
// so a **tensile** girder stress of 84 MPa raises the transverse stress at which
// the plating yields from 355 to 389 MPa and a **compressive** one lowers it to
// 305. Ignoring the pre-load therefore *under*-states the ferry's side above the
// neutral axis in hogging and over-states it below, and in sagging the other way
// round. The sign is a property of where the patch is and which way she is
// bending, and it is not available from the magnitude.
//
// One thing the closed form does **not** predict: first yield under a punch is a
// *bending* event at the clamped edge, which presents both signs of surface stress
// at once, so both signs of pre-load bring first yield forward. A membrane
// argument about the sign of that would be wrong and is not made.
//
// --- 5. The reaction back -----------------------------------------------------
//
// A zone that has torn is no longer carrying what the beam thinks it is.
// `reactionOf()` measures, per panel the zone meshed, the fraction of its
// thickness still working:
//
//   * a torn element carries nothing, so it counts zero;
//   * an intact element counts the thickness it **actually has now**, taken as
//     volume over mid-surface area from its deformed geometry -- measured, not
//     modelled from a plastic strain and an assumed Poisson ratio, and measured
//     against the thickness it had *when the zone was promoted*, because an
//     element's volume over its area is not the nominal plate thickness on a
//     curved patch and that difference is not damage;
//   * the part of a panel the zone did not mesh counts as intact, because nothing
//     looked at it.
//
// `reduce()` then folds that into the `StructuralMesh` as a thinner ship, and
// **nothing in Tier-0 changes at all**: `hullGirderSection`, `girderStress`,
// `girderBuckling`, `collapseElementsAt` and `longitudinalStrength` all already
// read a thickness. Section area falls, second moment falls, the plate's critical
// buckling stress falls as `t^2`, the ultimate moment falls and the applied
// utilisation rises -- none of it a new model.
//
// **Two of those are not monotone, and both are measured rather than assumed:**
//
//   * The section modulus at the **undamaged** fibre can *rise*. Damage above the
//     neutral axis pulls the axis down, and the far fibre's lever shrinks faster
//     than the second moment does: measured on the ferry, a zone at z = 8 m moves
//     the axis 6.713 -> 6.681 m and takes `modulusKeel` from 6.875 to 6.890 m^3.
//     So a far-fibre section modulus is not a conservative reading of damage. The
//     **ultimate moment** is, and it is what the conservatism claim is made
//     against.
//   * The **sagging** ultimate moment rises by up to 0.35% over part of a damage
//     sweep, because a thinner panel buckles earlier, therefore sheds earlier, and
//     Smith's method lets the neutral axis migrate in response. It is in the model
//     rather than in the quadrature -- eightfold the curvature steps leaves it
//     where it is -- which is why the claim asserted is "never more than intact"
//     and not "monotone".
//
// **What the reduction does not carry**, in the un-conservative direction and
// worth being loud about: the zone meshes plating only (`zone.hpp` §3), so the
// stiffeners running through a torn panel are left at full strength. A collision
// that opens fourteen bays of side shell has certainly destroyed the longitudinals
// in them, and `reduce()` will not say so. It needs the multi-point constraint
// that would let the zone mesh a web in the first place.
//
// The other thing it does not carry is that a **dented** panel is far weaker in
// compression than a merely thinner one, because the dent is an initial
// imperfection. `dentedCompressiveCapacity()` computes that from the measured
// deviation by Perry-Robertson, and it is deliberately *not* folded in: it is a
// compression-only knockdown and an effective thickness is not, so applying it
// would quietly weaken the panel in tension too. Folding it in needs a
// `collapse.hpp` that takes a per-element imperfection, which is a change to that
// file rather than to this one.
//
// SI units, body frame per CLAUDE.md.
#pragma once

#include "buckling.hpp"
#include "collapse.hpp"
#include "girder.hpp"
#include "scantlings.hpp"
#include "zone.hpp"

#include <string>
#include <vector>

namespace sim::promotion {

// --- What Tier-0 knows ---------------------------------------------------------

struct TierZeroParams {
    int stations = 41;
    double shedExponent = 0.45;
    int curvatureSteps = 150;
    // The Smith's-method sweep, which is 138 ms of the 170 ms this call costs on
    // the reference ferry. False skips it; `collapseUtilisation` is then zero and
    // the collapse trigger cannot fire, which is a caller's choice to make
    // knowingly rather than a cheaper answer to the same question.
    bool collapse = true;
    // Zero takes the plating's own. Explicit because a mixed-material ship has no
    // single yield strength and the utilisation is a ratio to one.
    double yieldStrength = 0;
};

// Everything the beam knows about a ship on a sea, in one call, plus the three
// utilisations the criterion reads. This is the expensive half of promotion and
// its cost is reported rather than assumed -- see §2 item 4.
struct TierZero {
    HullGirder girder;
    std::vector<GirderStress> stress;
    std::vector<GirderBuckling> buckling;
    std::vector<StrengthStation> strength;

    double yieldUtilisation = 0, yieldX = 0;        // worst |sigma| / sigma_y
    double buckleUtilisation = 0, buckleX = 0;      // worst applied / critical
    double collapseUtilisation = 0, collapseX = 0;  // worst |applied| / |ultimate|

    double yieldStrength = 0;   // Pa, what the utilisations are ratios to
    double seconds = 0;         // what this cost to compute
    std::vector<std::string> problems;
};

TierZero tierZero(const Ship& ship, const Sea& sea, const StructuralMesh& structure,
                  const Scantlings& scantlings, const TierZeroParams& params = {});

// --- The criterion -------------------------------------------------------------

enum class Trigger { None, Yield, Buckling, Collapse, Contact };
const char* name(Trigger trigger);

// A load applied to a place rather than to a section: a contact patch from
// `collision.hpp`, a grounding, a berthing. `force` is the resultant normal force
// and `radius` the footprint it is spread over; a touch carrying no force is not a
// promotion.
struct ContactPatch {
    Vec3 centre{};
    double radius = 0;  // m
    double force = 0;   // N
};

struct Criterion {
    // Absolute thresholds, per trigger, and the level each is held to once it has
    // promoted. `hold < promote` is the hysteresis band and the constructor of
    // `Promoter` says so if it is not.
    //
    // Yield sits highest because Tier-0 is *trustworthy* there -- a section
    // modulus is not an approximation to something better. Buckling sits lower
    // because that is where the beam stops being able to answer: past the critical
    // stress the problem is non-linear and Tier-0 only reports a ratio.
    double yieldPromote = 0.90, yieldHold = 0.80;
    double bucklePromote = 0.80, buckleHold = 0.70;
    double collapsePromote = 0.80, collapseHold = 0.70;

    // How far above the ship's own median a station has to stand, in utilisation
    // points. Zero switches the flat-ship guard off, which is what the negative
    // control in the tests does.
    double localExcess = 0.10;

    // A contact promotes when its mean pressure reaches this multiple of the
    // struck bay's plastic collapse pressure. One is the plating hinging.
    double contactPressure = 1.0, contactHold = 0.70;

    // Core-seconds of simulated second per element, for reporting only: 4.0 is the
    // measured figure for 12 mm plating (`zone.hpp` §1), and it scales as
    // `costThickness / t` because the stable step does. The **budget is in
    // elements**, which is what the zone is actually bounded by; this is the
    // conversion a caller wants printed and `zone::estimatedCost` is the exact
    // answer once a patch exists.
    double coreSecondsPerElement = 4.0;
    double costThickness = 0.012;  // m

    // Consecutive reviews a candidate must qualify for before it is promoted, and
    // consecutive reviews an active zone must fail to qualify for before it is
    // dropped. One and one is no dwell at all.
    int dwell = 2;
    int hold = 3;

    // Elements across every zone at once. 4000 is 16 000 core-seconds per
    // simulated second, so it is a statement about how much wall time a caller
    // will spend, not about how much memory it has.
    int elementBudget = 4000;

    // Minimum separation between two zones, in multiples of the zone radius. Two
    // is the smallest value at which they do not overlap and pay twice for the
    // same plating, which is what it is set to.
    //
    // It does **not** reduce a broad peak to one zone and is not meant to. A hull
    // girder's utilisation curve is broad -- on the reference ferry the buckling
    // check runs 0.31, 0.44, 0.33 over forty metres -- so a peak that clears the
    // threshold clears it over a real stretch of ship, and covering that stretch
    // takes several zones. What bounds the answer is `elementBudget`; what this
    // bounds is paying for the same plating twice.
    double separation = 2.0;

    // What a promoted zone gets. `role` is overridden per candidate from the panel
    // that was chosen, since a deck and a side shell are different roles.
    zone::MeshParams mesh;
};

struct Candidate {
    int panel = -1;      // index into StructuralMesh::panels -- the zone's identity
    Vec3 impact{};       // where the zone would be centred
    double x = 0;        // m, the station it came from
    PanelRole role = PanelRole::Shell;

    Trigger trigger = Trigger::None;
    double utilisation = 0;  // what the trigger measured
    double background = 0;   // the ship's own median for that trigger
    double score = 0;        // utilisation / the trigger's promote threshold; >= 1 qualifies

    int elements = 0;   // what a zone here would cost, estimated
    double cost = 0;    // core-seconds per simulated second
    std::string why;
};

// Every station and every contact that qualifies, ranked, deduplicated by panel,
// **before** separation and budget are applied. Exposed on its own because a
// caller -- or a test -- is entitled to see what was considered and rejected, and
// because it is the pure function inside the state machine.
std::vector<Candidate> candidates(const StructuralMesh& structure, const TierZero& tier,
                                  const std::vector<ContactPatch>& contacts,
                                  const Criterion& criterion);

// The plastic collapse pressure of a plate strip clamped at both ends:
// `4 sigma_y (t / span)^2`, from the three-hinge mechanism `w L^2 / 16 = M_p`.
// A real bay carries two-way action and is stronger, so this promotes slightly
// eagerly, which is the right direction for a criterion.
double platingCollapsePressure(double yieldStrength, double thickness, double span);

// The shorter of a panel's two side lengths -- the span the plating between
// stiffeners actually has. Averaged over each opposite pair, so a tapered panel
// gives its mean rather than the narrow end.
double panelSpan(const PlatePanel& panel);

// Elements a zone at `impact` would build, without building it: panels of the
// requested role within the radius, times `subdivision^2`. Over-states wherever
// the flood fill is truncated, which is the safe direction for a budget.
int estimateElements(const StructuralMesh& structure, const Vec3& impact,
                     const zone::MeshParams& mesh);

// --- The state machine ---------------------------------------------------------

struct Active {
    int panel = -1;
    Vec3 impact{};
    double x = 0;
    PanelRole role = PanelRole::Shell;
    Trigger trigger = Trigger::None;
    double score = 0;
    int elements = 0;
    double cost = 0;
    int promotedAtReview = 0;
    int idleReviews = 0;  // consecutive reviews it has not qualified for
};

struct Review {
    std::vector<Candidate> considered;  // ranked, before separation and budget
    std::vector<Active> promoted;       // newly promoted this review
    std::vector<Active> demoted;        // dropped this review
    int elementsActive = 0;
    double costActive = 0;  // core-seconds per simulated second, over all zones
    double microseconds = 0;  // what the decision itself cost -- not `TierZero`
    std::vector<std::string> problems;
};

// The decision, and nothing else: it does not mesh, solve or own a zone. That is
// what makes it cheap enough to run often and testable without a solver.
class Promoter {
public:
    explicit Promoter(Criterion criterion = {});

    Review review(const StructuralMesh& structure, const TierZero& tier,
                  const std::vector<ContactPatch>& contacts = {});

    const std::vector<Active>& active() const { return active_; }
    const Criterion& criterion() const { return criterion_; }
    int reviews() const { return reviews_; }
    // Cumulative, so a chatter test counts transitions rather than inferring them.
    int promotions() const { return promotions_; }
    int demotions() const { return demotions_; }
    int activeElements() const;
    double activeCost() const;
    void clear();

    // The mesh parameters a candidate's zone gets: the criterion's, with the role
    // and the outward direction taken from the panel that was chosen.
    zone::MeshParams meshFor(const Candidate& candidate) const;

private:
    Criterion criterion_;
    std::vector<Active> active_;
    std::vector<std::pair<int, int>> qualifying_;  // panel -> consecutive reviews, ascending
    int reviews_ = 0, promotions_ = 0, demotions_ = 0;
};

// --- The pre-load --------------------------------------------------------------

struct PreloadCheck {
    zone::Preload preload;
    bool applied = false;
    double obliquity = 0;      // rad the patch normal leans out of the athwartships plane
    double tractionError = 0;  // Pa left unbalanced on the patch's own faces: sigma sin^2(phi)
    double surfaceStress = 0;  // Pa the girder puts through the patch centre
    double neutralAxis = 0;    // m above the baseline
    double moment = 0;         // N m at the patch's station
    std::vector<std::string> problems;
};

// The girder stress through a patch, as a `zone::Preload`.
//
// Declines -- `applied` false, preload zero -- when the patch's normal leans
// further than `obliquityLimit` out of the athwartships plane, because the
// uniaxial state is then not traction-free on the patch's own faces and the zone
// would ring rather than sit. A transverse bulkhead is the extreme case and it is
// also the right answer: a transverse bulkhead carries no hull girder stress,
// which is exactly why `hullGirderSection` leaves it out of the section.
PreloadCheck preloadFor(const HullGirder& girder, const StructuralMesh& structure,
                        const zone::Patch& patch, double obliquityLimit = 0.20);

// --- The reaction back ---------------------------------------------------------

struct PanelDamage {
    int panel = -1;
    double effectiveness = 1;   // share of its thickness still carrying, in [0, 1]
    double meshedFraction = 0;  // how much of the panel the zone looked at
    double outOfPlane = 0;      // m, worst residual deviation along the patch axis
    double thinning = 0;        // m of thickness lost where it is still intact
    double tornArea = 0;        // m^2 of it deleted
};

struct SectionReduction {
    std::vector<PanelDamage> panels;  // ascending by panel index
    double xLo = 0, xHi = 0;          // m, the longitudinal reach of the damage
    double lostPlateArea = 0;         // m^2 of mid-surface no longer carrying
    double lostSteelMass = 0;         // kg
    double worstEffectiveness = 1;
    double worstOutOfPlane = 0;       // m
    std::vector<std::string> problems;

    bool empty() const { return panels.empty(); }
};

// What a solved zone says about the section it sits in. The solver is needed, not
// just its `SolveResult`: the residual geometry is where the thinning and the dent
// are measured, and a result carries neither. The structure is needed because a
// panel at the edge of a patch is only partly meshed, and the part nobody looked
// at has to count as intact -- otherwise the zone *radius* decides how much of the
// ship is damaged, which is the defect `indentation.hpp` records in its history.
SectionReduction reactionOf(const StructuralMesh& structure, const zone::Patch& patch,
                            const zone::Solver& solver);

// The ship Tier-0 should now believe in. Panel indices are preserved, so
// `breach.hpp`'s view of the same structure still lines up; a fully lost panel
// keeps its place with zero thickness rather than being erased.
StructuralMesh reduce(const StructuralMesh& structure, const SectionReduction& reduction);

// Compressive capacity left in a panel dented out of plane, by Perry-Robertson on
// the measured deviation as an initial imperfection:
//
//     sigma (1 + eta / (1 - sigma / sigma_cr)) = sigma_y,   eta = 6 * deviation / t
//
// `eta` is `w0 c / r^2` for a strip bending about its own mid-plane, where
// `c = t/2` and `r^2 = t^2/12`. Zero deviation returns `min(sigma_y, sigma_cr)`
// exactly, so it is a knockdown with no free coefficient in it.
//
// **Not folded into the section reduction**, and §5 says why.
double dentedCompressiveCapacity(double yieldStrength, double criticalStress, double deviation,
                                 double thickness);

}  // namespace sim::promotion
