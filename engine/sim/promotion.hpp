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
// A fourth answer now sits beside them and it is not structural at all:
// `GasPromoter` decides which *burning compartment* deserves `les.hpp`'s resolved
// flow instead of `fire.hpp`'s two zones. It is here rather than in `les.hpp` for
// the reason `Promoter` is here rather than in `zone.hpp` -- a promotion criterion
// is a decision about cost and evidence, it needs no solver, and the two are the
// same state machine over different evidence. §6 has the criterion.
//
// Where the three-tier plan in `02-simulation.md` §3 says a zone couples to Tier-1
// through retained interface DOF, this couples it to Tier-0 through a section.
// **That coupling now exists -- `coupling.{hpp,cpp}` -- and this file is still
// unchanged by it, which is what it was written for.** The two are not the same
// coupling wearing different clothes and neither replaces the other:
//
//   * Tier 0 reads a **thickness** and nothing else, so damage goes back to it as
//     a thinner ship (§5) and `hullGirderSection` needs no new model.
//   * Tier 1 reads a **mesh**, so damage goes back to it as the elements that
//     tore, deleted. A thickness knockdown is not even expressible there: a
//     solid-shell carries its thickness in the positions of its nodes, so thinning
//     a zone moves the very interface nodes it is coupled through, by (t - t')/2
//     -- 1.2 mm for a 20% knockdown on 12 mm plating, measured, against a
//     coincidence tolerance of 1e-9 m. `coupling.hpp` §3 has the argument in full.
//
// So the section reduction below is Tier 0's answer, permanently, and not a
// placeholder for the interface one. What is still outstanding is the *criterion*
// side: nothing here yet consults a Tier-1 model when deciding what to promote,
// and `reduction::checkValidity` is the trigger that would.
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
//     slamming, sloshing, a dropped weight, and a fire against a bulkhead.
//
//     **And "softening" was the wrong word for the last of those**, which the
//     Phase 4 milestone settled: the ferry's engine room bulkhead fails at a
//     member temperature of 151.6 C, where `k_y` is *exactly 1* and the steel has
//     lost none of its strength. What a fire does to a restrained member is put it
//     into compression -- `thermal::HeatedMember` -- and what fells it is that
//     compression magnifying the bending the head of water behind it is already
//     applying, by 2.086. So a promotion criterion that looked for lost strength
//     would find nothing to promote right up to the moment the bulkhead went. The
//     trigger this file would need is a *temperature field*, not a knock-down.
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
// **The reduction now carries the stiffeners too, and this is where the
// un-conservative direction used to be.** It used to read: the members running
// through a torn panel are left at full strength, so a collision that opens
// fourteen bays of side shell leaves their longitudinals intact and `reduce()` will
// not say so. The multi-point constraint had arrived and `zone::Stiffeners::Modelled`
// was meshing the member; what was missing was on the other side -- a fibre had no
// damage variable and was never deleted, so there was no "this longitudinal is
// gone" to read. `constraint.hpp` §2b closed that, and this consumes it.
//
// A member is reduced **exactly as a panel is**: `MemberDamage::effectiveness` is
// the share of the steel the zone meshed for it that is still carrying, the part
// nobody looked at counts as intact, and `reduce()` folds it in by scaling the
// profile's **thicknesses** -- web and flange together. That is the exact analogue
// of scaling a plate's thickness, and for the same reason: area, the profile's own
// second moment and its Steiner term `A d^2` are all *linear* in the rectangle
// thicknesses at fixed height, so one factor moves all three consistently and the
// profile's centroid does not move at all. Scaling the web *height* instead would
// move the centroid and is a different damage.
//
// **Both thicknesses, and snapped to zero together.** `profileSection` reads a tee's
// flange whenever `flangeThickness > 0` regardless of the web, so scaling the web
// alone to nothing leaves a flange floating on no web with area and a Steiner term
// intact. `sectionElements` already drops a member whose `ProfileSection::area` is
// zero, which is what makes a fully torn longitudinal leave the hull girder
// altogether rather than linger as an epsilon -- the same trap `reactionOf`'s
// `1e-6` snap exists for on the plating side, and the same fix.
//
// **What it still does not carry.** The knockdown is uniform over the profile, and
// the fibre that tears first is the *outer* one -- the one furthest from the plate,
// carrying the most of the member's own second moment. So a partly torn member is
// reported slightly strong about its own axis. About the hull girder's axis, which
// is the one Tier 0 asks about, the error is the ratio of `I_own` to `A d^2` and for
// a 200 mm bar metres from the neutral axis it is small. And a fibre fails on axial
// damage alone: a longitudinal whose plating has been *deleted* from under it is
// unsupported, not merely undamaged, and that is `coupling::withoutTornElements`'
// question rather than this one.
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
#include "les.hpp"
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

    // Core-seconds of simulated second per element, for reporting only: 1.7 is
    // the measured figure for 12 mm plating (`zone.hpp` §1), and it scales as
    // `costThickness / t` because the stable step does. The **budget is in
    // elements**, which is what the zone is actually bounded by; this is the
    // conversion a caller wants printed and `zone::estimatedCost` is the exact
    // answer once a patch exists.
    //
    // **This was 4.0, and 4.0 is the figure `zone.hpp` §1 warns about by name.**
    // An elastoplastic element cost 7.3 µs per step until `cacheRestForms`
    // stopped rebuilding the step-invariant element forms every step; at 5.5e5
    // steps per simulated second that is 4.0 core-seconds per element, and it is
    // what this field held while citing as its source the very paragraph that
    // says "every figure below that predates it is 2.4x pessimistic per element".
    // The live cost is `kPlasticMicroseconds = 3.1` µs (`zone.cpp:27`), so
    // 3.1e-6 x 5.5e5 = 1.7.
    //
    // Nothing caught it because nothing reads this field outside its own
    // arithmetic: `grep coreSecondsPerElement tests/ tools/` is empty, and the
    // one cost assertion in the suite checks `activeCost() > 0`. The tier has a
    // second, independent cost estimator in `zone::estimatedCost` (`zone.cpp:644`)
    // which computes off `kPlasticMicroseconds` and the patch's own critical
    // timestep -- so the two disagreed by 2.35x, `zone_probe` printed both on the
    // same run, and no one had read them side by side.
    double coreSecondsPerElement = 1.7;
    double costThickness = 0.012;  // m

    // Consecutive reviews a candidate must qualify for before it is promoted, and
    // consecutive reviews an active zone must fail to qualify for before it is
    // dropped. One and one is no dwell at all.
    int dwell = 2;
    int hold = 3;

    // Elements across every zone at once. 4000 is ~6 800 core-seconds per
    // simulated second, so it is a statement about how much wall time a caller
    // will spend, not about how much memory it has. (It read 16 000 while
    // `coreSecondsPerElement` was the stale 4.0; the budget itself has not moved.)
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

// The same statement about a stiffener. `effectiveness` is the share of the
// member's *section* still carrying, measured on the steel volume the zone's fibres
// stand for -- see the header, and `constraint::memberDamage` for why volume rather
// than any other weighting.
struct MemberDamage {
    int member = -1;
    double effectiveness = 1;   // share of its section still carrying, in [0, 1]
    double meshedFraction = 0;  // how much of the member's steel the zone looked at
    double lostVolume = 0;      // m^3 of stiffener steel deleted
    int fibers = 0, tornFibers = 0;
};

struct SectionReduction {
    std::vector<PanelDamage> panels;  // ascending by panel index
    std::vector<MemberDamage> members;  // ascending by member index
    double xLo = 0, xHi = 0;          // m, the longitudinal reach of the damage
    double lostPlateArea = 0;         // m^2 of mid-surface no longer carrying
    double lostSteelMass = 0;         // kg, plating and stiffeners together
    double worstEffectiveness = 1;
    double worstMemberEffectiveness = 1;
    double worstOutOfPlane = 0;       // m
    std::vector<std::string> problems;

    // Still "are there damaged panels": the plating is what `breach.hpp` and the
    // flooding model read, and a reduction with members in it and no panels is a
    // zone that bent its longitudinals without opening anything. `nothing()` is the
    // stricter question, for a caller that wants to know whether `reduce()` would
    // change the ship at all.
    bool empty() const { return panels.empty(); }
    bool nothing() const { return panels.empty() && members.empty(); }
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

// --- 6. The gas side: when a two-zone compartment deserves a resolved model -----
//
// Everything above is structural, and everything below is not; what they share is
// the *shape of the decision*, which is why they live in one file. A fidelity tier
// is worth paying for when the cheap model is answering a question it cannot
// answer, the decision has to be cheap enough to take often, and it has to not
// chatter. `Promoter` and `GasPromoter` are the same state machine over different
// evidence -- absolute threshold, hysteresis band, dwell, budget -- and neither
// owns a solver.
//
// **What the cheap model asserts here is that the upper layer is one temperature.**
// `fire.hpp` has two well-mixed control volumes per compartment, so a bulkhead
// forty metres from a fire sees exactly the same gas as the deckhead above it. That
// is right in a small room and wrong in a long one, and *how* wrong is a published
// closed form rather than a judgement: Alpert's ceiling-jet correlation has two
// branches either side of `r/H = 0.18`, and their ratio is
//
//     dT(0) / dT(r) = (16.9 / 5.38) (r / H)^(2/3)
//
// with the heat release **cancelling exactly**, and passing through one at the
// crossover to 0.1427% -- which is the closest thing to a validation two independent
// branches of one correlation can give. So the geometric trigger is a
// property of the compartment and not of the fire in it -- a bigger fire does not
// switch it on -- and it needs an activity trigger beside it, which is the layer
// temperature rise. Neither implies the other and both are required: a cold ro-ro
// deck has the shape and nothing to resolve, a fierce fire in a tank has plenty to
// resolve and a shape a zone model gets right.
//
// `radius` is the **half-diagonal of the plan rectangle** `les::planRectangle`
// derives from the compartment's own area and perimeter: the distance from a fire
// at the centre to the farthest corner, which is how far the ceiling jet has to
// travel before it stops being one temperature. A `fire::DesignFire` carries no
// horizontal position at all, so the fire is at the centre -- which is the *least*
// favourable reading for promotion, since a fire in a corner is further from the far
// bulkhead than this says, and understating the case is the right direction for a
// criterion that spends core-seconds.
//
// **The threshold is anchored on the published guidance rather than picked.** Zone
// models are held to be sound for compartments up to about `L/H = 3`; for a square
// room the half-diagonal is `L/sqrt(2)`, so `r/H = 2.121` there and the spread is
// `3.1413 * 2.121^(2/3) = 5.16`. `spreadPromote = 5.0` is that boundary, and the
// consequence is that a 6 x 5 x 4 m tank (3.09) and a 12 x 10 x 5 m machinery space
// (4.24) stay two-zone while a 20 x 8 x 5 m hold (5.27) and the ferry's own
// 100 x 19 m vehicle deck (14.8) do not.
//
// **Cost here is not linear in the resolution, and that is the difference from §1.**
// A Tier-2 zone costs `4.0 x elements` core-seconds per simulated second, exactly,
// because the step is fixed. Halving this model's cell size multiplies the cells by
// eight *and* halves the advective step *and* lengthens the pressure solve, so the
// measured cost per cell per simulated second **rises** with refinement: 1.5e-6 at
// 168 cells, 8.9e-6 at 600, 1.2e-5 at 1456, on one core. `cellBudget` therefore
// bounds the arithmetic honestly and the wall time only loosely, and
// `coreSecondsPerCell` is a figure for the resolution the default parameters give
// rather than a constant of the model.
//
// The budget is in **cells**, on the same terms the structural budget is in
// elements, and `les::estimateCells` is what it is spent against.

struct GasCriterion {
    // Alpert's near/far ratio at which the two-zone model's single upper-layer
    // temperature is asserting away a factor this large. Five is the `L/H = 3` a
    // zone model is conventionally trusted to -- see above; it is a threshold with a
    // published anchor rather than a round number.
    double spreadPromote = 5.0, spreadHold = 4.0;

    // And the fire has to be doing something: the upper layer standing this far
    // above the lower one. A compartment with no stratification has no layer whose
    // uniformity could be wrong.
    double risePromote = 20.0, riseHold = 10.0;   // K

    // Consecutive reviews a candidate must qualify for before it is promoted, and
    // consecutive reviews a resolved compartment must fail to qualify for before it
    // is dropped. One and one is no dwell at all.
    int dwell = 2;
    int hold = 3;

    // Cells across every resolved compartment at once. One compartment at 0.5 m on
    // a machinery space is a few thousand, so this is a statement about how many
    // compartments may be resolved rather than about how fine any of them is.
    int cellBudget = 12000;

    // Core-seconds of wall time per simulated second per cell, for reporting only,
    // exactly as `Criterion::coreSecondsPerElement` is. Measured in
    // `tests/test_promotion.cpp` rather than assumed.
    double coreSecondsPerCell = 1.0e-5;

    les::Params grid;
};

struct GasCandidate {
    int compartment = -1;    // index into fire::Model::gas -- the identity
    std::string name;
    double spread = 0;       // Alpert's near/far ratio for this compartment
    double rise = 0;         // K, upper layer over lower
    // How far past its own threshold the **weaker** of the two triggers is, so a
    // score of one means both have just cleared. Ranking on the weaker one is what
    // stops a compartment with an enormous fire and a cubical shape from outranking
    // one the zone model is actually wrong about.
    double score = 0;
    int    cells = 0;        // what a resolved model here would cost, estimated
    double cost = 0;         // core-seconds per simulated second
    std::string why;
};

// Every tracked compartment that qualifies, ranked, before the budget is applied.
// Exposed for the same reason `candidates()` is: a caller is entitled to see what
// was considered and rejected.
std::vector<GasCandidate> gasCandidates(const fire::Model& model, const GasCriterion& criterion);

// Half the diagonal of the plan rectangle a compartment's area and perimeter imply
// -- the distance from a fire at its centre to the farthest corner. Falls back to
// `sqrt(A/2)`, the square's own half-diagonal, when no rectangle has both.
double compartmentReach(double floorArea, double perimeter);

struct GasActive {
    int compartment = -1;
    std::string name;
    double spread = 0, rise = 0, score = 0;
    int cells = 0;
    double cost = 0;
    int promotedAtReview = 0;
    int idleReviews = 0;
};

struct GasReview {
    std::vector<GasCandidate> considered;  // ranked, before the budget
    std::vector<GasActive> promoted;       // newly promoted this review
    std::vector<GasActive> demoted;        // dropped this review
    int cellsActive = 0;
    double costActive = 0;
    double microseconds = 0;
    std::vector<std::string> problems;
};

// The decision, and nothing else: it does not build a grid, step one, or own one.
class GasPromoter {
public:
    explicit GasPromoter(GasCriterion criterion = {});

    GasReview review(const fire::Model& model);

    const std::vector<GasActive>& active() const { return active_; }
    const GasCriterion& criterion() const { return criterion_; }
    int reviews() const { return reviews_; }
    int promotions() const { return promotions_; }
    int demotions() const { return demotions_; }
    int activeCells() const;
    double activeCost() const;
    void clear();

private:
    GasCriterion criterion_;
    std::vector<GasActive> active_;
    std::vector<std::pair<int, int>> qualifying_;  // compartment -> consecutive reviews
    int reviews_ = 0, promotions_ = 0, demotions_ = 0;
};

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
