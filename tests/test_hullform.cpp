// SPDX-License-Identifier: MIT
//
// Validation of hull generation from principal particulars.
//
// A hull generator is easy to write and hard to check, because any smooth blob
// looks like a ship. So nothing here compares against a picture or against the
// numbers that were requested -- that would be a tautology, since the request is
// the input. Everything is checked by *measuring the mesh back*, with the same
// integrator the hydrostatics use, and comparing against either a closed form or
// an independently computed value:
//
//   * the bilge radius and the midship coefficient are exact inverses;
//   * the area curve's Cp and LCB are analytic, and are checked against
//     quadrature of the curve itself -- two independent routes to one number;
//   * Cm = Cb = 1 makes a rectangular box, whose displacement is L*B*T exactly;
//   * everything else has to *converge* under refinement, not merely be small,
//     because a small constant error and a discretisation error look identical
//     at a single resolution.
#include "engine/sim/hullform.hpp"
#include "harness.hpp"

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

using namespace sim;
using testing::expectEqual;
using testing::expectNear;
using testing::expectTrue;

namespace {

// --- Closed forms ------------------------------------------------------------

void testBilgeRadiusAndMidshipCoefficientAreInverses() {
    const double beam = 25.4, draft = 9.5;
    double worst = 0;
    for (double cm : {0.999, 0.99, 0.95, 0.90, 0.88}) {
        const double radius = bilgeRadiusForMidshipCoefficient(beam, draft, cm);
        const double back = midshipCoefficientForBilgeRadius(beam, draft, radius);
        worst = std::max(worst, std::abs(back - cm));
        expectTrue("the bilge radius fits inside the section",
                   radius <= std::min(0.5 * beam, draft) + 1e-12);
    }
    expectTrue("radius and midship coefficient invert exactly", worst < 1e-12);

    // A full rectangle has no bilge at all.
    expectNear("Cm of 1 needs no bilge radius",
               bilgeRadiusForMidshipCoefficient(beam, draft, 1.0), 0.0, 1e-12);

    // And the section area integrates to what the coefficient claims. This is the
    // link between the formula and the geometry it is supposed to describe.
    const double cm = 0.92;
    const double radius = bilgeRadiusForMidshipCoefficient(beam, draft, cm);
    const int steps = 200000;
    double area = 0;
    for (int i = 0; i < steps; ++i) {
        const double z = draft * (i + 0.5) / steps;
        area += 2.0 * midshipHalfBreadth(beam, draft, radius, z) * (draft / steps);
    }
    expectNear("the section shape really has that midship coefficient", area / (beam * draft),
               cm, 1e-5);
}

// The area curve reports Cp and LCB from closed forms. Quadrature of the curve
// is a completely separate route to the same two numbers, so agreement means the
// algebra is right rather than merely self-consistent.
void testAreaCurveIntegralsMatchQuadrature() {
    for (double cp : {0.55, 0.65, 0.80}) {
        for (double lcb : {-0.03, 0.0, 0.025}) {
            for (double parallel : {0.0, 0.35}) {
                const double transom = parallel > 0.0 ? 0.15 : 0.0;
                const AreaCurve curve = solveAreaCurve(cp, lcb, transom, 0.01, parallel);

                const int steps = 400000;
                double area = 0, moment = 0;
                for (int i = 0; i < steps; ++i) {
                    const double u = -1.0 + 2.0 * (i + 0.5) / steps;
                    const double f = curve(u);
                    area += f * (2.0 / steps);
                    moment += u * f * (2.0 / steps);
                }
                expectNear("analytic Cp matches quadrature", curve.prismaticCoefficient(),
                           area / 2.0, 1e-6);
                expectNear("analytic LCB matches quadrature", curve.lcbFraction(),
                           moment / area * 0.5, 1e-6);

                // And the solve actually hit what it was asked for.
                expectNear("the solve reached the requested Cp", curve.prismaticCoefficient(), cp,
                           1e-6);
                expectNear("the solve reached the requested LCB", curve.lcbFraction(), lcb, 1e-6);
            }
        }
    }
}

// A parallel middle body is a span of *constant* section, so the area curve must
// be flat over it and the hull's stations identical there. Nothing about Cb, Cp
// or LCB can see the difference -- which is why its absence went unnoticed until
// a rendered frame showed a canoe.
void testParallelMiddleBodyIsActuallyParallel() {
    const AreaCurve curve = solveAreaCurve(0.70, 0.0, 0.12, 0.01, 0.40);
    expectNear("the solve still reaches Cp with a parallel middle body",
               curve.prismaticCoefficient(), 0.70, 1e-6);
    expectNear("and still reaches LCB", curve.lcbFraction(), 0.0, 1e-6);

    for (double u : {-0.39, -0.2, 0.0, 0.2, 0.39})
        expectNear("the section is constant over the parallel middle body", curve(u), 1.0, 1e-12);
    expectTrue("and tapers outside it", curve(0.7) < 0.999 && curve(-0.7) < 0.999);

    // Without one, the same request has no flat span at all.
    const AreaCurve tapered = solveAreaCurve(0.70, 0.0, 0.12, 0.01, 0.0);
    expectTrue("a curve without a middle body is nowhere flat", tapered(0.2) < 0.999);

    // And the hull built from it really does repeat its stations. This is the
    // check that ties the curve to the geometry rather than to itself.
    HullParticulars p = s175Particulars();
    p.parallelMiddleBodyFraction = 0.40;
    p.stationCount = 41;
    const TriMesh hull = makeHullFromParticulars(p);
    expectTrue("the hull is still a closed manifold", isClosedManifold(hull));

    const HullCoefficients c = measureHull(hull, p.draft, p.lengthPp, p.beam);
    expectNear("and still measures as the ship that was asked for", c.blockCoefficient,
               p.blockCoefficient, 0.005 * p.blockCoefficient);
    expectNear("with LCB where it was asked", c.lcbFraction, p.lcbFraction, 1e-3);

    // Two slabs inside the parallel span must have equal sectional area; one
    // outside must not.
    const double halfLength = 0.5 * p.lengthPp;
    Vec3 lo = hull.verts[0], hi = hull.verts[0];
    for (const Vec3& v : hull.verts) {
        lo = {std::min(lo.x, v.x), std::min(lo.y, v.y), std::min(lo.z, v.z)};
        hi = {std::max(hi.x, v.x), std::max(hi.y, v.y), std::max(hi.z, v.z)};
    }
    const double slab = 0.004 * p.lengthPp;
    const auto sectionArea = [&](double x) {
        const TriMesh cut = clipToBox(hull, {x - 0.5 * slab, lo.y - 1e3, lo.z - 1e3},
                                      {x + 0.5 * slab, hi.y + 1e3, lo.z + p.draft});
        return integrate(cut).volume / slab;
    };
    const double inner = sectionArea(0.0);
    const double alsoInner = sectionArea(0.30 * halfLength);
    const double outer = sectionArea(0.80 * halfLength);
    expectTrue("the section really is non-trivial", inner > 1.0);
    expectNear("two stations inside the parallel body have the same area", alsoInner, inner,
               1e-6 * inner);
    expectTrue("and one outside it does not", outer < 0.9 * inner);
}

// --- The mesh ----------------------------------------------------------------

// Cm = Cb = 1 is a rectangular box, and a box of L x B x T displaces exactly
// L*B*T. No tolerance argument, no tessellation error: flat faces are exact.
void testAFullBlockIsExactlyABox() {
    HullParticulars p;
    p.lengthPp = 80.0;
    p.beam = 20.0;
    p.draft = 5.0;
    p.depth = 10.0;
    p.blockCoefficient = 0.999;
    p.midshipCoefficient = 1.0;
    p.lcbFraction = 0.0;
    p.stemFraction = 1.0;      // no taper at all
    p.transomFraction = 1.0;

    const TriMesh hull = makeHullFromParticulars(p);
    expectTrue("the generated hull is a closed manifold", isClosedManifold(hull));

    const HullCoefficients c = measureHull(hull, p.draft, p.lengthPp, p.beam);
    expectNear("a full block displaces exactly L*B*T", c.displacedVolume,
               p.lengthPp * p.beam * p.draft, 1e-6 * p.lengthPp * p.beam * p.draft);
    expectNear("its block coefficient is 1", c.blockCoefficient, 1.0, 1e-9);
    expectNear("its midship coefficient is 1", c.midshipCoefficient, 1.0, 1e-6);
    expectNear("its waterplane coefficient is 1", c.waterplaneCoefficient, 1.0, 1e-6);
    expectNear("and its centre of buoyancy is amidships", c.lcbFraction, 0.0, 1e-9);
}

// The generated hull must measure as the ship that was asked for -- and the
// residual must *shrink* under refinement, which is what distinguishes a
// discretisation error from a mistake.
void testGeneratedHullMeasuresAsRequestedAndConverges() {
    const HullParticulars base = s175Particulars();

    double previous = 1e30;
    int improved = 0;
    for (int waterlines : {11, 21, 41}) {
        HullParticulars p = base;
        p.waterlineCount = waterlines;
        const TriMesh hull = makeHullFromParticulars(p);
        expectTrue("every refinement is still a closed manifold", isClosedManifold(hull));

        const HullCoefficients c = measureHull(hull, p.draft, p.lengthPp, p.beam);
        const double error = std::abs(c.blockCoefficient - p.blockCoefficient);
        if (error < previous) ++improved;
        previous = error;

        // LCB comes from the area curve analytically and does not depend on the
        // bilge tessellation, so it is right at every resolution.
        expectNear("LCB lands where it was asked to", c.lcbFraction, p.lcbFraction, 5e-4);
    }
    expectEqual("refining the waterlines improves the block coefficient every time", improved, 3);
    expectTrue("and the finest is within a quarter of a percent", previous < 0.0025 * base.blockCoefficient);
}

// A ship's coefficients are not independent: Cp = Cb / Cm is a definition, so a
// generated hull has to satisfy it as measured, not merely as requested.
void testMeasuredCoefficientsAreSelfConsistent() {
    HullParticulars p = kvlcc2Particulars();
    const TriMesh hull = makeHullFromParticulars(p);
    const HullCoefficients c = measureHull(hull, p.draft, p.lengthPp, p.beam);

    expectNear("Cp equals Cb over Cm as measured", c.prismaticCoefficient,
               c.blockCoefficient / c.midshipCoefficient, 1e-9);
    expectNear("displacement equals Cb L B T", c.displacedVolume,
               c.blockCoefficient * p.lengthPp * p.beam * p.draft,
               1e-9 * c.displacedVolume);
    expectNear("the tanker comes out at its published block coefficient",
               c.blockCoefficient, p.blockCoefficient, 0.005 * p.blockCoefficient);
}

// LCB has a sign, and getting it backwards produces a ship that trims the wrong
// way while looking entirely normal.
void testLcbMovesTheRightWay() {
    HullParticulars aft = s175Particulars(), forward = s175Particulars();
    aft.lcbFraction = -0.03;
    forward.lcbFraction = 0.03;

    const HullCoefficients a =
        measureHull(makeHullFromParticulars(aft), aft.draft, aft.lengthPp, aft.beam);
    const HullCoefficients f =
        measureHull(makeHullFromParticulars(forward), forward.draft, forward.lengthPp,
                    forward.beam);

    expectTrue("a forward LCB really is forward of an aft one", f.lcbFraction > a.lcbFraction);
    expectNear("and both land where they were asked", a.lcbFraction, -0.03, 1e-3);
    expectNear("forward too", f.lcbFraction, 0.03, 1e-3);
    // Guard: the two must differ by roughly what was asked, or "moves the right
    // way" would pass on two nearly identical hulls.
    expectNear("the separation is the one requested", f.lcbFraction - a.lcbFraction, 0.06, 2e-3);
}

// --- Refusals ----------------------------------------------------------------

void testUnreasonableParticularsAreReported() {
    HullParticulars p;
    p.blockCoefficient = 0.85;
    p.midshipCoefficient = 0.80;   // Cb > Cm means Cp > 1
    const std::vector<std::string> problems = validateParticulars(p);
    expectTrue("Cb above Cm is refused", !problems.empty());

    HullParticulars ordinary = s175Particulars();
    expectTrue("an ordinary ship raises nothing", validateParticulars(ordinary).empty());

    // A bilge radius cannot exceed the draft, so a low enough Cm on a shallow
    // hull is simply unreachable -- and must be said, not silently rounded up.
    HullParticulars shallow;
    shallow.lengthPp = 142.0;
    shallow.beam = 19.1;
    shallow.draft = 6.15;
    shallow.depth = 12.5;
    shallow.blockCoefficient = 0.507;
    shallow.midshipCoefficient = 0.83;
    std::vector<std::string> reported;
    const TriMesh hull = makeHullFromParticulars(shallow, &reported);
    expectTrue("an unreachable midship coefficient is reported", !reported.empty());
    expectTrue("and the hull is still built and usable", isClosedManifold(hull));

    // Having said so, the hull it does build must be consistent: Cb is honoured
    // against the section it actually has.
    const HullCoefficients c =
        measureHull(hull, shallow.draft, shallow.lengthPp, shallow.beam);
    expectNear("the achievable hull still hits the requested block coefficient",
               c.blockCoefficient, shallow.blockCoefficient, 0.01 * shallow.blockCoefficient);
}

}  // namespace

void runHullFormTests() {
    std::printf("\n--- hull generation from particulars ---\n");
    testBilgeRadiusAndMidshipCoefficientAreInverses();
    testAreaCurveIntegralsMatchQuadrature();
    testParallelMiddleBodyIsActuallyParallel();
    testAFullBlockIsExactlyABox();
    testGeneratedHullMeasuresAsRequestedAndConverges();
    testMeasuredCoefficientsAreSelfConsistent();
    testLcbMovesTheRightWay();
    testUnreasonableParticularsAreReported();
}
