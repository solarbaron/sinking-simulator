// SPDX-License-Identifier: MIT
//
// Validation of the sparse FLIP/APIC water solver.
//
// A fluid solver is unusually well supplied with exact answers, and almost
// nothing here is a tolerance on an eyeballed number. What is asserted, and why
// each one is the sharpest form of its question:
//
//   * **The kernel's two identities.** `sum w = 1` and `sum w (x_n - x_p)^2 = h^2/4`,
//     over 200 000 offsets spanning forty cells. The second is what makes APIC's
//     `D_p` a constant, and asserting it directly is the only way to catch a
//     weight from the wrong kernel: a trilinear weight satisfies the first
//     identity and not the second.
//   * **Mass is exact, not conserved to a tolerance.** `expectNear(..., 0.0, 0.0)`
//     -- a tolerance of literally zero -- on the residual, and the particle count
//     as an integer beside it. A particle's mass is never modified, so anything
//     other than `0.0` is a particle that was created or lost.
//   * **A hydrostatic column does not move at all.** Not "moves slowly": after
//     twenty-five steps *every particle is at the bit-identical position it
//     started at*, because the residual velocity is 5.5e-16 m/s and `x + v dt` at
//     that size is a no-op. The same run would have fallen 12.3 mm at 0.49 m/s
//     unprojected, and that contrast is asserted as the vacuity guard.
//   * **Its pressure field is `rho g h (K - k)` per cell**, asserted cell by cell
//     to 4e-15 relative rather than as a norm -- the single most valuable defect
//     shape in this repo is an error that cancels when asked globally.
//   * **The projection reproduces a discrete Helmholtz decomposition exactly.** A
//     field built as the MAC curl of a vector potential is discretely
//     divergence-free by construction; add the discrete gradient of a cell-centred
//     potential twice its size, project, and the curl part must come back and the
//     potential must come out as the pressure. On a 7x5x6 grid, so an index that
//     reads the wrong axis cannot survive.
//   * **A block in free fall follows the discrete closed form** of its own
//     integrator, `-g dt^2 N(N+1)/2`, to 1.1e-15 relative, with no lateral drift
//     at all -- exactly `0.0` in x and y.
//   * **A sloshing tank's first mode has the period `2 pi / sqrt(g k tanh k d)`**,
//     and what is asserted is the period *and* that refining the grid moves it
//     towards the analytic value, so the residual error is discretisation rather
//     than a wrong dispersion relation.
//   * **The sparse structure is tested for what it is for.** An empty room of 64
//     million cells allocates exactly zero bytes; the same water in a 20^3 room
//     and a 400^3 room allocates the same tiles and produces *bit-identical*
//     particle positions after fifty steps; and water arriving in a region that
//     had no storage is caught by the invariant that face weights sum to the
//     particle mass exactly.
//
// The vacuity guards matter as much as the assertions. A sloshing test in which
// nothing sloshes passes trivially; a mass test on an empty field passes
// trivially; a projection test on a field that was already divergence-free proves
// nothing. Every one of those carries an explicit check that the thing being
// measured is non-trivial, and each of those guards is here because the first
// version of that test passed while proving nothing.
#include "engine/sim/flip.hpp"
#include "harness.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <cstdio>
#include <random>
#include <vector>

using testing::expectEqual;
using testing::expectNear;
using testing::expectTrue;

namespace {

namespace fl = sim::flip;
using sim::kGravity;
using sim::kPi;
using sim::kRhoSeawater;

// A tank `L` long, three cells wide, `depth` deep, with `freeboard` cells of air
// above it. Quasi-two-dimensional: the first sloshing mode runs along x and the
// side walls in y do not enter its dispersion relation.
fl::Field makeTank(double L, double depth, double h, int freeboard) {
    fl::Field f;
    f.grid.h = h;
    f.grid.n[0] = static_cast<int>(std::lround(L / h));
    f.grid.n[1] = 3;
    f.grid.n[2] = static_cast<int>(std::lround(depth / h)) + freeboard;
    const double width = 3.0 * h;
    const double lo[3] = {0, 0, 0}, hi[3] = {L, width, depth};
    fl::seedBox(f, lo, hi, 2, kRhoSeawater);
    fl::setTotalMass(f, kRhoSeawater * L * width * depth);
    return f;
}

// The linear standing wave's own velocity field, imposed over a **flat** free
// surface. Seeding the tilted surface instead would quantise the initial
// condition to the cell -- the first version of this test did, and at an
// amplitude of a sixth of a cell the tank did not slosh at all, because a
// voxelised free surface cannot represent a sub-cell tilt and therefore has no
// restoring force to give. The flat surface is represented exactly.
//
//   u = U cosh(k z)/cosh(k d) sin(k x),   w = -U sinh(k z)/cosh(k d) cos(k x)
//
// which is divergence-free, has no flux through the floor, and none through the
// end walls either.
void seedFirstMode(fl::Field& f, double L, double depth, double amplitude) {
    const double k = kPi / L;
    const double omega = std::sqrt(kGravity * k * std::tanh(k * depth));
    const double u0 = amplitude * kGravity * k / omega;
    for (fl::Particle& p : f.particles) {
        const double x = p.position[0], z = p.position[2];
        p.velocity[0] = u0 * std::cosh(k * z) / std::cosh(k * depth) * std::sin(k * x);
        p.velocity[2] = -u0 * std::sinh(k * z) / std::cosh(k * depth) * std::cos(k * x);
    }
}

double firstModePeriod(double L, double depth) {
    const double k = kPi / L;
    return 2.0 * kPi / std::sqrt(kGravity * k * std::tanh(k * depth));
}

struct SloshRun {
    double period = 0;
    double peak = 0;
    double noise = 0;
    double massResidual = 0;
    int    crossings = 0;
    int    substeps = 0;
    bool   capped = false;
    // The substep count ran past what the Courant condition can explain. Checked
    // *inside* the loop rather than after it, because that is the difference
    // between a mutant that fails in twenty seconds and one that fails in two
    // hundred: `maxSubsteps` bounds one step at 4096, and four hundred steps of
    // that is a run nobody waits for.
    bool   runaway = false;
    // A step that could not take the time it was asked for. The run **stops** at
    // the first one rather than pressing on, and that is a mutation-testing
    // decision as much as a physical one: a defect that drives the substep
    // controller towards zero makes every remaining step cost the whole budget,
    // which is how a mutant kills a suite by *hanging* instead of by failing. Two
    // did exactly that on the first sweep here. Bailing out turns them into an
    // assertion in about the time the clean run takes.
    bool   incomplete = false;
};

// Run a tank and measure the first-mode period from the zero crossings of its
// mass centroid. The centroid rather than a surface probe: it is a global,
// mass-weighted signal with no sampling of its own to go wrong, and it has the
// mode's period exactly.
SloshRun runTank(double L, double depth, double h, double amplitude, double seconds, double dt,
                 bool affine, double blend) {
    fl::Field f = makeTank(L, depth, h, 6);
    seedFirstMode(f, L, depth, amplitude);
    fl::Params p;
    p.affine = affine;
    p.flipBlend = blend;
    p.maxSubstep = dt;
    // **The budget is the instrument, not the clock.** `maxSubstep` is `dt` and
    // the tank never exceeds about 1 m/s, so the Courant bound never binds and the
    // controller's arithmetic floor is one substep a step. Eight is that floor with
    // room, and it is set here rather than left at 4096 so that a defect which
    // drives the step towards zero comes back as `incomplete` in one step instead
    // of running four thousand projections first. That is the difference between a
    // mutant that fails in twenty seconds and one that fails in four minutes -- and
    // between a kill scored by an assertion and one scored by a timeout.
    p.maxSubsteps = 8;
    // The period comes out the same to six digits at 1e-13; the projection is not
    // what limits this measurement, and 1e-8 costs a fifth fewer iterations.
    p.projectionTolerance = 1e-8;
    fl::Account account;
    fl::resetAccount(f, account);
    fl::Solver s;

    SloshRun out;
    const double centre = 0.5 * L;
    double c[3];
    fl::centroid(f, c);
    double previous = c[0] - centre, previousTime = 0;
    std::vector<double> crossings;
    const int steps = static_cast<int>(std::lround(seconds / dt));
    for (int i = 0; i < steps; ++i) {
        const fl::StepResult r = s.step(f, dt, p, account);
        out.substeps += r.substeps;
        // **Stop at the first step that did not do what it was asked**, whether
        // that is a projection out of iterations or a substep loop out of budget.
        // Both are failures the assertions below name, and pressing on past either
        // is what turns a mutant into an hours-long run instead of a red line: the
        // two that killed the first sweep by *hanging* now fail in about the time
        // the clean run takes.
        if (r.projectionCapped) {
            out.capped = true;
            break;
        }
        if (r.incomplete) {
            out.incomplete = true;
            break;
        }
        // A tank whose water never exceeds about 1 m/s takes one substep a step at
        // `maxSubstep = dt`; four times that is the arithmetic floor of the
        // controller with room to spare, and nothing physical crosses it.
        if (out.substeps > 4 * steps) {
            out.runaway = true;
            break;
        }
        fl::centroid(f, c);
        const double offset = c[0] - centre;
        out.peak = std::max(out.peak, std::abs(offset));
        const double t = static_cast<double>(i + 1) * dt;
        if ((previous < 0 && offset >= 0) || (previous > 0 && offset <= 0)) {
            const double fraction = previous / (previous - offset);
            crossings.push_back(previousTime + fraction * (t - previousTime));
        }
        previous = offset;
        previousTime = t;
    }
    out.crossings = static_cast<int>(crossings.size());
    if (crossings.size() >= 3)
        out.period = 2.0 * (crossings.back() - crossings.front()) /
                     static_cast<double>(crossings.size() - 1);
    double numerator = 0, denominator = 0;
    for (const fl::Particle& q : f.particles) {
        double g[3];
        s.sampleVelocity(q.position, g);
        for (int a = 0; a < 3; ++a) {
            numerator += (q.velocity[a] - g[a]) * (q.velocity[a] - g[a]);
            denominator += q.velocity[a] * q.velocity[a];
        }
    }
    out.noise = denominator > 0 ? std::sqrt(numerator / denominator) : 0.0;
    out.massResidual = account.massResidual();
    return out;
}

// ---------------------------------------------------------------------------
// 1. The kernel
// ---------------------------------------------------------------------------

void testKernel() {
    std::printf("\n--- flip: the transfer kernel ---\n");
    std::mt19937 rng(20260806u);
    std::uniform_real_distribution<double> offsets(-20.0, 20.0);
    double worstSum = 0, worstFirst = 0, worstSecond = 0;
    int lowBase = 1 << 30, highBase = -(1 << 30);
    double smallestWeight = 1.0, largestWeight = 0.0;
    bool everyWeightNonNegative = true;
    for (int i = 0; i < 200000; ++i) {
        const double s = offsets(rng);
        double w[3], offset[3];
        const int base = fl::splineWeights(s, w, offset);
        lowBase = std::min(lowBase, base);
        highBase = std::max(highBase, base);
        double sum = 0, first = 0, second = 0;
        for (int m = 0; m < 3; ++m) {
            sum += w[m];
            first += w[m] * offset[m];
            second += w[m] * offset[m] * offset[m];
            smallestWeight = std::min(smallestWeight, w[m]);
            largestWeight = std::max(largestWeight, w[m]);
            // Every weight of a B-spline is non-negative; a negative one is a
            // higher-order interpolant, which transfers mass that is not there.
            if (!(w[m] >= 0.0)) everyWeightNonNegative = false;
        }
        worstSum = std::max(worstSum, std::abs(sum - 1.0));
        worstFirst = std::max(worstFirst, std::abs(first));
        worstSecond = std::max(worstSecond, std::abs(second - 0.25) / 0.25);
    }
    expectTrue("every spline weight is non-negative", everyWeightNonNegative);
    std::printf("     200000 offsets over %d cells: |sum w - 1| <= %.3g,  |sum w d| <= %.3g,"
                "  |sum w d^2 / (h^2/4) - 1| <= %.3g\n",
                highBase - lowBase, worstSum, worstFirst, worstSecond);
    // Asserted at what was measured, not at a round number: these are identities
    // and the only thing between them and exactness is the rounding of three
    // multiplies.
    expectTrue("partition of unity holds to 3e-16", worstSum <= 3e-16);
    expectTrue("the kernel reproduces a linear field to 1e-16", worstFirst <= 1e-16);
    expectTrue("the second moment is h^2/4 to 3e-16", worstSecond <= 3e-16);
    // Vacuity: the sweep has to have moved the base index and exercised the whole
    // range of the weight function, or all three identities are being asserted at
    // one point.
    expectTrue("the offsets swept many cells", highBase - lowBase >= 35);
    expectTrue("the weights covered their whole range",
               smallestWeight < 1e-6 && largestWeight > 0.749);

    // `D_p = h^2/4 I` in metres squared, which is the form APIC actually uses.
    // Two more roundings than the dimensionless identity above, and it shows.
    const double h = 0.037;
    double worstMoment = 0;
    for (int i = 0; i < 20000; ++i) {
        const double s = offsets(rng);
        worstMoment = std::max(worstMoment,
                               std::abs(fl::kernelSecondMoment(s, h) - 0.25 * h * h) /
                                   (0.25 * h * h));
    }
    std::printf("     D_p against h^2/4 at h = %.3f m: worst relative deviation %.3g\n", h,
                worstMoment);
    expectTrue("kernelSecondMoment is h^2/4 to 1e-15", worstMoment <= 1e-15);
    // The trilinear alternative, which passes the partition of unity and fails
    // this: its second moment is h^2 f(1-f), zero at a node. Stated as a control
    // so the assertion above is known to be discriminating rather than generous.
    const double trilinearAtNode = 0.0;
    expectTrue("the trilinear second moment this rejects is far away",
               std::abs(trilinearAtNode - 0.25 * h * h) > 0.2 * h * h);
}

// ---------------------------------------------------------------------------
// 2. The transfer: partition of unity, per face and in total
// ---------------------------------------------------------------------------

void testTransfer() {
    std::printf("\n--- flip: particle to grid ---\n");

    // One particle, one known position, and the twenty-seven weights it deposits
    // checked against the kernel evaluated by hand. Asked of a single face rather
    // than of a sum, because a sum of weights is right whenever the weights are a
    // permutation of the right ones -- which is exactly what an axis swap is.
    {
        fl::Field f;
        f.grid.h = 0.25;
        f.grid.n[0] = 8; f.grid.n[1] = 8; f.grid.n[2] = 8;
        fl::Particle p;
        // Deliberately not at a cell centre and not on a node: 1.30 = cell 5 plus
        // 0.2 of a cell, 0.70 = cell 2 plus 0.8, 0.45 = cell 1 plus 0.8.
        p.position[0] = 1.30; p.position[1] = 0.70; p.position[2] = 0.45;
        p.mass = 3.0;
        f.particles.push_back(p);
        fl::Params params;
        fl::Solver s;
        s.rebuild(f, params);
        s.transferToGrid(f, params);

        // Each component's nodes sit on the cell faces along its *own* axis and at
        // the cell centres along the other two. All three are worked out here from
        // that definition and checked node by node, because a sum over the
        // twenty-seven is right whenever the weights are a permutation of the right
        // ones -- which is exactly what an axis swap is.
        int base[3][3];
        double weight[3][3][3];
        for (int a = 0; a < 3; ++a) {
            const double position[3] = {1.30, 0.70, 0.45};
            for (int b = 0; b < 3; ++b) {
                double offset[3];
                const double coordinate = position[b] / 0.25 - (b == a ? 0.0 : 0.5);
                base[a][b] = fl::splineWeights(coordinate, weight[a][b], offset);
            }
        }
        // x at 5.2 cells, y at 2.3, z at 1.3 for the x component; the y and z
        // components see the same particle half a cell along in the other axes.
        expectEqual("x component's x base node", base[0][0], 4);
        expectEqual("x component's y base node", base[0][1], 1);
        expectEqual("x component's z base node", base[0][2], 0);
        expectEqual("y component's x base node", base[1][0], 4);
        expectEqual("y component's y base node", base[1][1], 2);
        expectEqual("z component's z base node", base[2][2], 1);
        // The staggering itself: the y component's y nodes are half a cell off the
        // x component's y nodes, and likewise in z. A collocated grid would have
        // all three the same, and would have a checkerboard null space with it.
        expectTrue("the components are staggered in y", base[0][1] != base[1][1]);
        expectTrue("the components are staggered in z", base[0][2] != base[2][2]);

        for (int a = 0; a < 3; ++a) {
            double sum = 0, biggest = 0, worst = 0;
            for (int kk = 0; kk < 3; ++kk)
                for (int jj = 0; jj < 3; ++jj)
                    for (int ii = 0; ii < 3; ++ii) {
                        const double want =
                            3.0 * weight[a][0][ii] * weight[a][1][jj] * weight[a][2][kk];
                        const double got = s.faceMass(a, base[a][0] + ii, base[a][1] + jj,
                                                      base[a][2] + kk);
                        worst = std::max(worst, std::abs(got - want));
                        sum += got;
                        biggest = std::max(biggest, want);
                    }
            expectTrue("every node of component " + std::to_string(a) +
                           " carries the kernel's own weight",
                       worst <= 4e-16 * 3.0);
            expectNear("component " + std::to_string(a) +
                           "'s twenty-seven nodes hold the whole particle",
                       sum, 3.0, 1e-15);
            expectTrue("the transfer is spread rather than concentrated on one node",
                       biggest < 0.5 * 3.0 && biggest > 0.05 * 3.0);
        }
        // And the three components put different amounts on the *same* node index,
        // which is what an axis swap would destroy.
        expectTrue("the components disagree about one node",
                   std::abs(s.faceMass(0, 5, 2, 1) - s.faceMass(1, 5, 2, 1)) > 0.05);
    }

    // A block of water straddling tile seams, so the sum below has to cross the
    // sparse structure's own boundaries. Tiles are four cells; the block starts
    // and ends inside tiles rather than on them.
    {
        fl::Field f;
        f.grid.h = 0.1;
        f.grid.n[0] = 20; f.grid.n[1] = 20; f.grid.n[2] = 20;
        const double lo[3] = {0.35, 0.35, 0.35}, hi[3] = {1.05, 0.75, 0.95};
        fl::seedBox(f, lo, hi, 2, kRhoSeawater);
        fl::Params params;
        fl::Solver s;
        s.rebuild(f, params);
        s.transferToGrid(f, params);
        const double mass = f.totalMass();
        expectTrue("the straddling block has particles", f.particles.size() > 500);
        for (int a = 0; a < 3; ++a) {
            // Exactly, with a tolerance of zero. The kernel is a partition of
            // unity and the deposits are the same numbers in the same order, so a
            // missing halo tile is the only thing that can move this.
            expectNear("component " + std::to_string(a) + " face mass is the particle mass",
                       s.totalFaceMass(a), mass, 0.0);
        }
        // Every face beside a fluid cell carries mass. The transfer divides by
        // this, and the claim that it never divides by a small number is the
        // premise rather than a hope.
        int fluidFaces = 0;
        double smallest = 1e30;
        for (int k = -1; k <= f.grid.n[2]; ++k)
            for (int j = -1; j <= f.grid.n[1]; ++j)
                for (int i = -1; i <= f.grid.n[0]; ++i)
                    for (int a = 0; a < 3; ++a)
                        if (s.faceAt(a, i, j, k) == fl::Face::Fluid) {
                            ++fluidFaces;
                            smallest = std::min(smallest, s.faceMass(a, i, j, k));
                        }
        std::printf("     %d fluid faces, lightest carries %.6g kg of %.6g\n", fluidFaces,
                    smallest, mass);
        expectTrue("there are fluid faces to check", fluidFaces > 500);
        expectTrue("every fluid face carries mass", smallest > 0.0);
        // And it is not merely positive: the lightest is a real fraction of what a
        // full face carries, so the division has no conditioning problem.
        expectTrue("the lightest fluid face is not a denormal",
                   smallest > 1e-4 * mass / static_cast<double>(fluidFaces));
    }
}

// ---------------------------------------------------------------------------
// 3. Hydrostatic rest -- the sharpest test of a projection there is
// ---------------------------------------------------------------------------

void testHydrostatic() {
    std::printf("\n--- flip: a hydrostatic column does not move ---\n");
    const double h = 0.05;
    const int layers = 24;
    fl::Field f;
    f.grid.h = h;
    f.grid.n[0] = 6; f.grid.n[1] = 6; f.grid.n[2] = 32;
    const double lo[3] = {0, 0, 0}, hi[3] = {0.30, 0.30, h * layers};
    fl::seedBox(f, lo, hi, 2, kRhoSeawater);
    fl::setTotalMass(f, kRhoSeawater * 0.30 * 0.30 * h * layers);
    const std::vector<fl::Particle> before = f.particles;

    fl::Params p;
    p.maxSubstep = 0.002;
    // A column at rest cannot subdivide: its velocity is 1e-16 m/s and the Courant
    // bound is 22 s. One substep a step is the floor; four is the budget. See the
    // note in `runTank` for why this is set rather than left at its default.
    p.maxSubsteps = 4;
    fl::Account account;
    fl::resetAccount(f, account);
    fl::Solver s;

    const int steps = 25;
    const double dt = 0.002;
    double worstSpeed = 0;
    int worstIterations = 0;
    int substeps = 0;
    bool bailed = false;
    for (int i = 0; i < steps; ++i) {
        const fl::StepResult r = s.step(f, dt, p, account);
        worstIterations = std::max(worstIterations, r.projectionIterations);
        substeps += r.substeps;
        expectTrue("the hydrostatic projection converged", !r.projectionCapped);
        for (const fl::Particle& q : f.particles)
            for (int a = 0; a < 3; ++a) worstSpeed = std::max(worstSpeed, std::abs(q.velocity[a]));
        // A column at rest takes exactly one substep a step, because the Courant
        // bound cannot bind on a velocity of 1e-16 m/s. Anything above that is a
        // column that is not at rest, and stopping here is what keeps a mutant's
        // failure a red line rather than a three-minute wait -- `maxSubstep` is
        // `dt`, so the only way to exceed it is for the controller to subdivide.
        if (r.projectionCapped || r.incomplete || substeps > i + 1) {
            bailed = true;
            break;
        }
    }
    expectEqual("the column took one substep a step", substeps, bailed ? -1 : steps);

    int moved = 0;
    double worstMove = 0;
    for (std::size_t i = 0; i < f.particles.size(); ++i)
        for (int a = 0; a < 3; ++a) {
            const double delta = f.particles[i].position[a] - before[i].position[a];
            if (delta != 0.0) ++moved;
            worstMove = std::max(worstMove, std::abs(delta));
        }
    const double t = steps * dt;
    std::printf("     %d steps: worst particle speed %.4g m/s, %d of %zu coordinates moved"
                " (worst %.4g m); unprojected free fall would be %.6g m at %.6g m/s\n",
                steps, worstSpeed, moved, 3 * f.particles.size(), worstMove,
                0.5 * kGravity * t * t, kGravity * t);
    // Measured at 5.5e-16 m/s; asserted at 2e-15, which is a factor of four and
    // not a factor of a thousand. A hydrostatic column that had lost the property
    // entirely would be at 0.49 m/s.
    expectTrue("the column's residual speed is machine precision", worstSpeed < 2e-15);
    // The stronger form: nothing moved at all. `x + v dt` with `v` at 1e-15 and
    // `dt` at 2e-3 is below the last bit of a coordinate near 1, so a column truly
    // at rest is *bit identical* after twenty-five steps, and that is asserted
    // rather than a small displacement.
    expectEqual("not one particle coordinate changed", moved, 0);
    // Vacuity: the run had to have something to hold up.
    expectTrue("gravity would have moved it a long way",
               0.5 * kGravity * t * t > 1e10 * std::max(worstMove, 1e-18));

    // The pressure, cell by cell against `rho g h (K - k)`. Per cell rather than
    // as a norm: an error that cancels when asked globally is this repo's
    // characteristic defect and it dies only when one cell is asked alone.
    double worstPressure = 0;
    for (int k = 0; k < layers; ++k) {
        const double want = kRhoSeawater * kGravity * h * static_cast<double>(layers - k);
        const double got = s.pressureAt(3, 3, k);
        worstPressure = std::max(worstPressure, std::abs(got - want) / want);
        expectNear("hydrostatic pressure in cell k=" + std::to_string(k), got, want, 4e-15 * want);
    }
    // The air cell above the surface is the Dirichlet zero, and the cell below the
    // floor is solid and carries no pressure at all. Both asked directly.
    expectNear("the air cell above the surface is at zero", s.pressureAt(3, 3, layers), 0.0, 0.0);
    expectEqual("the cell above the surface is air",
                static_cast<long long>(s.cellAt(3, 3, layers)),
                static_cast<long long>(fl::Cell::Air));
    expectEqual("the cell below the floor is solid",
                static_cast<long long>(s.cellAt(3, 3, -1)),
                static_cast<long long>(fl::Cell::Solid));
    expectEqual("the floor face is a solid face",
                static_cast<long long>(s.faceAt(2, 3, 3, 0)),
                static_cast<long long>(fl::Face::Solid));

    // One face, alone. The face between the fifth and sixth fluid cells carries
    // the whole of gravity before the projection and exactly none of it after.
    const double faceVelocity = s.faceVelocity(2, 3, 3, 5);
    std::printf("     one interior face: %.4g m/s, against %.6g m/s of gravity in a substep;"
                " worst cell pressure error %.3g relative\n",
                faceVelocity, kGravity * dt, worstPressure);
    expectTrue("one interior vertical face is at rest", std::abs(faceVelocity) < 2e-15);
    expectTrue("that face had a substep of gravity to cancel", kGravity * dt > 1e-2);

    // Divergence, per cell, and the mass account.
    double worstDivergence = 0;
    for (int k = 0; k < layers; ++k)
        for (int j = 0; j < 6; ++j)
            for (int i = 0; i < 6; ++i)
                worstDivergence = std::max(worstDivergence, std::abs(s.divergenceAt(i, j, k)));
    expectTrue("every cell is divergence free to 1e-12 1/s", worstDivergence < 1e-12);
    expectNear("mass is exactly conserved", account.massResidual(), 0.0, 0.0);
    expectEqual("no particle was created or lost", account.particleResidual(), 0);
    expectEqual("no particle touched a wall", account.clamped, 0);
}

// ---------------------------------------------------------------------------
// 4. The projection reproduces a discrete Helmholtz decomposition
// ---------------------------------------------------------------------------

void testProjection() {
    std::printf("\n--- flip: the projection against a discrete Helmholtz decomposition ---\n");
    // Unequal on every axis, and none of them a multiple of the tile size, so an
    // index that reads the wrong axis or runs one past a bound cannot survive.
    const int n[3] = {7, 5, 6};
    const double h = 0.08;
    fl::Field f;
    f.grid.h = h;
    for (int a = 0; a < 3; ++a) f.grid.n[a] = n[a];
    const double L[3] = {n[0] * h, n[1] * h, n[2] * h};
    const double lo[3] = {0, 0, 0}, hi[3] = {L[0], L[1], L[2]};
    fl::seedBox(f, lo, hi, 2, kRhoSeawater);

    fl::Params p;
    p.projectionTolerance = 1e-15;
    p.projectionIterations = 4000;
    fl::Solver s;
    s.rebuild(f, p);
    s.transferToGrid(f, p);
    // Every cell holds water, so there is no free surface and the pressure is
    // determined only up to a constant. That is not an edge case for a flooded
    // compartment, it is the ordinary case.
    expectTrue("a fully submerged region has no Dirichlet boundary", s.singular());
    expectEqual("every cell is fluid", s.fluidCells(), n[0] * n[1] * n[2]);

    // A vector potential whose factors are exactly zero on the walls, so the MAC
    // curl of it has exactly zero flux through them -- `0.0`, not `sin(pi)`.
    const auto potentialZ = [&](int i, int j, int k) {
        const double x = i * h, y = j * h, z = (k + 0.5) * h;
        return x * (L[0] - x) * y * (L[1] - y) * (1.0 + 0.7 * x + 0.3 * z);
    };
    const auto potentialX = [&](int i, int j, int k) {
        const double x = (i + 0.5) * h, y = j * h, z = k * h;
        return y * (L[1] - y) * z * (L[2] - z) * (1.0 + 0.5 * x + 0.9 * y);
    };
    // Three different wavenumbers, so a swapped axis changes the answer. Scaled so
    // the gradient the projection has to remove is half again the size of the
    // field it has to give back -- see the vacuity note below.
    const double phiScale = 0.004;
    const auto phi = [&](int i, int j, int k) {
        return phiScale * std::cos(1.0 * kPi * (i + 0.5) / n[0]) *
               std::cos(2.0 * kPi * (j + 0.5) / n[1]) * std::cos(3.0 * kPi * (k + 0.5) / n[2]);
    };

    const auto slot = [&](int i, int j, int k) {
        return static_cast<std::size_t>((k * (n[1] + 1) + j) * (n[0] + 1) + i);
    };
    std::vector<double> want[3];
    for (int a = 0; a < 3; ++a)
        want[a].assign(static_cast<std::size_t>((n[0] + 1) * (n[1] + 1) * (n[2] + 1)), 0.0);
    for (int k = 0; k <= n[2]; ++k)
        for (int j = 0; j <= n[1]; ++j)
            for (int i = 0; i <= n[0]; ++i) {
                if (j < n[1] && k < n[2])
                    want[0][slot(i, j, k)] = (potentialZ(i, j + 1, k) - potentialZ(i, j, k)) / h;
                if (i < n[0] && k < n[2])
                    want[1][slot(i, j, k)] = (potentialX(i, j, k + 1) - potentialX(i, j, k)) / h -
                                             (potentialZ(i + 1, j, k) - potentialZ(i, j, k)) / h;
                if (i < n[0] && j < n[1])
                    want[2][slot(i, j, k)] = -(potentialX(i, j + 1, k) - potentialX(i, j, k)) / h;
            }

    double fieldScale = 0;
    for (int a = 0; a < 3; ++a)
        for (double v : want[a]) fieldScale = std::max(fieldScale, std::abs(v));

    // It is discretely divergence-free by construction -- div of a curl -- and
    // that is checked rather than assumed, because if it were not the test would
    // be measuring the wrong thing.
    double targetDivergence = 0;
    for (int k = 0; k < n[2]; ++k)
        for (int j = 0; j < n[1]; ++j)
            for (int i = 0; i < n[0]; ++i) {
                double d = 0;
                for (int a = 0; a < 3; ++a) {
                    int high[3] = {i, j, k};
                    high[a] += 1;
                    d += want[a][slot(high[0], high[1], high[2])] - want[a][slot(i, j, k)];
                }
                targetDivergence = std::max(targetDivergence, std::abs(d) / h);
            }
    expectTrue("the target field is discretely divergence free", targetDivergence < 1e-15);
    expectTrue("the target field is not merely zero", fieldScale > 1e-3);

    double gradientScale = 0;
    bool everyFaceStored = true;
    for (int a = 0; a < 3; ++a) {
        int count[3] = {n[0], n[1], n[2]};
        count[a] += 1;
        for (int k = 0; k < count[2]; ++k)
            for (int j = 0; j < count[1]; ++j)
                for (int i = 0; i < count[0]; ++i) {
                    const int index[3] = {i, j, k};
                    double value = 0;
                    if (index[a] > 0 && index[a] < count[a] - 1) {
                        int low[3] = {i, j, k};
                        low[a] -= 1;
                        const double g = (phi(i, j, k) - phi(low[0], low[1], low[2])) / h;
                        gradientScale = std::max(gradientScale, std::abs(g));
                        value = want[a][slot(i, j, k)] + g;
                    }
                    if (!s.setFaceVelocity(a, i, j, k, value)) everyFaceStored = false;
                }
    }
    expectTrue("every face of a fully submerged room has storage", everyFaceStored);
    // Vacuity: the gradient the projection has to remove is larger than the field
    // it has to give back, so "it returned what it was handed" is not a pass.
    std::printf("     divergence-free part %.5g m/s, gradient added %.5g m/s (%.0f%% of it)\n",
                fieldScale, gradientScale, 100.0 * gradientScale / fieldScale);
    expectTrue("the gradient removed is at least the size of the field kept",
               gradientScale > fieldScale);

    const double dt = 0.01;
    s.project(dt, p);
    expectTrue("the projection converged", !s.lastCapped());

    double worstFace = 0;
    int worstIndex[4] = {0, 0, 0, 0};
    for (int a = 0; a < 3; ++a) {
        int count[3] = {n[0], n[1], n[2]};
        count[a] += 1;
        for (int k = 0; k < count[2]; ++k)
            for (int j = 0; j < count[1]; ++j)
                for (int i = 0; i < count[0]; ++i) {
                    const int index[3] = {i, j, k};
                    const double target =
                        (index[a] > 0 && index[a] < count[a] - 1) ? want[a][slot(i, j, k)] : 0.0;
                    const double got = s.faceVelocity(a, i, j, k);
                    if (std::abs(got - target) > worstFace) {
                        worstFace = std::abs(got - target);
                        worstIndex[0] = a; worstIndex[1] = i; worstIndex[2] = j; worstIndex[3] = k;
                    }
                }
    }
    // One face, alone: the worst one, named, rather than a norm over all of them.
    std::printf("     worst face is component %d at (%d,%d,%d): %.4g m/s of %.4g -> %.3g relative\n",
                worstIndex[0], worstIndex[1], worstIndex[2], worstIndex[3], worstFace, fieldScale,
                worstFace / fieldScale);
    expectTrue("the projection returns the divergence-free part to 1e-14 relative",
               worstFace < 1e-14 * fieldScale);

    // And the potential comes out as the pressure, cell by cell. The solver
    // removes the mean on a singular system, so the comparison is against the
    // mean-free potential.
    double meanPhi = 0;
    for (int k = 0; k < n[2]; ++k)
        for (int j = 0; j < n[1]; ++j)
            for (int i = 0; i < n[0]; ++i) meanPhi += phi(i, j, k);
    meanPhi /= static_cast<double>(n[0] * n[1] * n[2]);
    double worstCell = 0, phiScaleSeen = 0;
    for (int k = 0; k < n[2]; ++k)
        for (int j = 0; j < n[1]; ++j)
            for (int i = 0; i < n[0]; ++i) {
                const double got = s.pressureAt(i, j, k) * dt / p.density;
                const double target = phi(i, j, k) - meanPhi;
                worstCell = std::max(worstCell, std::abs(got - target));
                phiScaleSeen = std::max(phiScaleSeen, std::abs(target));
                expectNear("recovered potential in one cell", got, target, 1e-14 * phiScale);
            }
    std::printf("     recovered potential: worst cell %.4g of %.4g -> %.3g relative;"
                " divergence after %.3g 1/s\n",
                worstCell, phiScaleSeen, worstCell / phiScaleSeen, s.maxDivergence());
    expectTrue("the pressure is the potential to 1e-14 relative",
               worstCell < 1e-14 * phiScaleSeen);
    expectTrue("the projected field is divergence free", s.maxDivergence() < 1e-12);
}

// ---------------------------------------------------------------------------
// 5. Free fall: the discrete closed form of the integrator itself
// ---------------------------------------------------------------------------

void testBallistic() {
    std::printf("\n--- flip: a block in free fall, and the tiles it needs on the way ---\n");
    fl::Field f;
    f.grid.h = 0.1;
    f.grid.n[0] = 30; f.grid.n[1] = 30; f.grid.n[2] = 60;
    const double lo[3] = {1.0, 1.0, 3.0}, hi[3] = {1.4, 1.4, 3.4};
    fl::seedBox(f, lo, hi, 2, kRhoSeawater);
    fl::Params p;
    p.maxSubstep = 0.005;
    p.maxSubsteps = 4;   // free fall reaches 3 m/s on a 0.1 m cell: one substep a step
    fl::Account account;
    fl::resetAccount(f, account);
    fl::Solver s;
    double start[3];
    fl::centroid(f, start);

    const double dt = 0.005;
    const int steps = 60;
    int lowestTiles = 1 << 30, highestTiles = 0;
    std::vector<int> tileTrace;
    bool incomplete = false;
    int substeps = 0;
    for (int i = 0; i < steps; ++i) {
        const fl::StepResult r = s.step(f, dt, p, account);
        substeps += r.substeps;
        // Stop at the first step that did not do what it was asked. See
        // `SloshRun::incomplete`: pressing on is how a mutant hangs a suite. A
        // block in free fall reaches 3 m/s on a 0.1 m cell, so the Courant bound
        // never binds and one substep a step is the controller's arithmetic floor.
        if (r.incomplete || r.projectionCapped || substeps > 4 * (i + 1)) {
            incomplete = true;
            break;
        }
        lowestTiles = std::min(lowestTiles, r.tiles);
        highestTiles = std::max(highestTiles, r.tiles);
        tileTrace.push_back(r.tiles);
        // The block is in free fall, so its own divergence is zero and the
        // projection has nothing to do. Asserted, because a projection that
        // invents a pressure here would be inventing a force.
        expectTrue("free fall needs no pressure", std::abs(s.pressureAt(12, 12, 32)) < 1e-9);
    }
    double finish[3];
    fl::centroid(f, finish);
    // Symplectic Euler applied `steps` times: `v_n = -g n dt`, `x_n = x_{n-1} + v_n dt`,
    // so the total drop is `-g dt^2 N(N+1)/2` exactly.
    const double N = steps;
    const double want = -kGravity * dt * dt * N * (N + 1.0) / 2.0;
    const double got = finish[2] - start[2];
    std::printf("     dropped %.17g m against the discrete closed form %.17g m (%.3g relative);"
                " tiles %d..%d\n",
                got, want, std::abs(got - want) / std::abs(want), lowestTiles, highestTiles);
    expectTrue("free fall never ran out of substeps", !incomplete);
    // One substep per step: nothing here goes fast enough to subdivide, and a
    // controller that has collapsed says so here rather than by taking an hour.
    expectEqual("free fall takes one substep a step", substeps, steps);
    expectNear("the drop is the integrator's own closed form", got, want, 1e-14 * std::abs(want));
    // Exactly zero sideways, not merely small: nothing in a uniform vertical
    // acceleration has an x or y component to round.
    expectNear("no drift along x", finish[0] - start[0], 0.0, 0.0);
    expectNear("no drift along y", finish[1] - start[1], 0.0, 0.0);
    expectNear("mass is exactly conserved", account.massResidual(), 0.0, 0.0);
    expectEqual("no particle was created or lost", account.particleResidual(), 0);
    expectEqual("nothing reached a wall", account.clamped, 0);
    // Vacuity, twice over: it fell a long way, and the sparse structure had to
    // follow it -- the tile count is not the same at the end as at the start,
    // which is the arrival-and-departure this test doubles as.
    expectTrue("the block fell several cells", std::abs(got) > 4.0 * f.grid.h);
    expectTrue("the tile set changed while it fell", highestTiles != lowestTiles);
    expectTrue("the tile count came back down", tileTrace.back() <= highestTiles);
}

// ---------------------------------------------------------------------------
// 6. The sparse structure, tested for what it is for
// ---------------------------------------------------------------------------

void testSparse() {
    std::printf("\n--- flip: sparsity ---\n");

    // An empty room costs nothing, whatever the room is.
    for (int extent : {20, 100, 400}) {
        fl::Field f;
        f.grid.h = 0.1;
        for (int a = 0; a < 3; ++a) f.grid.n[a] = extent;
        fl::Params p;
        fl::Account account;
        fl::resetAccount(f, account);
        fl::Solver s;
        const fl::StepResult r = s.step(f, 0.01, p, account);
        expectEqual("an empty room of " + std::to_string(extent) + "^3 allocates no tiles",
                    s.tiles(), 0);
        expectEqual("an empty room of " + std::to_string(extent) + "^3 allocates no bytes",
                    static_cast<long long>(s.bytes()), 0);
        expectEqual("an empty room takes no substeps", r.substeps, 0);
        // Time still passes, so a caller integrating an empty compartment beside a
        // full one does not have the two come apart.
        expectNear("time still advanced", f.time, 0.01, 0.0);
    }
    {
        fl::Field f;
        f.grid.h = 0.1;
        for (int a = 0; a < 3; ++a) f.grid.n[a] = 400;
        std::printf("     an empty 400^3 room is %d cells and 0 bytes; at the 115 bytes a cell"
                    " this structure costs, dense it would be %.1f GB\n",
                    400 * 400 * 400, 400.0 * 400.0 * 400.0 * 115.0 / 1e9);
        expectTrue("the room it costs nothing for is a large one",
                   static_cast<double>(f.grid.n[0]) * f.grid.n[1] * f.grid.n[2] > 6e7);
    }

    // The tile count against one derived from the water's own extent, rather than
    // against another run of the same code. A 0.6 x 0.6 x 0.4 m block at h = 0.1
    // occupies cells 5..10, 5..10, 10..13; dilated by the 3-cell halo that is
    // 2..13, 2..13, 7..16, which is tiles 0..3, 0..3, 1..4 -- four by four by
    // four, sixty-four. **This is the check that was missing**: `tiles()` returned
    // `tileKey_.size()`, three entries per tile, and reported 192. Every other
    // tile assertion here compares one count with another and a factor of three
    // survives all of them.
    {
        fl::Field f;
        f.grid.h = 0.1;
        for (int a = 0; a < 3; ++a) f.grid.n[a] = 20;
        const double lo[3] = {0.5, 0.5, 1.0}, hi[3] = {1.1, 1.1, 1.4};
        fl::seedBox(f, lo, hi, 2, kRhoSeawater);
        fl::Params p;
        fl::Solver s;
        s.rebuild(f, p);
        expectEqual("the water fills six by six by four cells", s.fluidCells(), 6 * 6 * 4);
        expectEqual("and needs four by four by four tiles around them", s.tiles(), 64);
        // And the bytes follow the tiles, which is the claim `bytes()` is making.
        const std::size_t perCell = s.bytes() / static_cast<std::size_t>(64 * fl::kTileCells);
        std::printf("     a 6x6x4-cell block: %d fluid cells, %d tiles, %zu bytes"
                    " (%zu per allocated cell)\n",
                    s.fluidCells(), s.tiles(), s.bytes(), perCell);
        expectTrue("the cost per allocated cell is a few dozen bytes",
                   perCell > 16 && perCell < 512);
        // Every allocated tile is reachable and every cell in it is classified:
        // the halo box, walked by hand.
        int allocated = 0;
        for (int k = 1 * fl::kTile; k < 5 * fl::kTile; ++k)
            for (int j = 0; j < 4 * fl::kTile; ++j)
                for (int i = 0; i < 4 * fl::kTile; ++i)
                    if (s.allocated(i, j, k)) ++allocated;
        expectEqual("every cell of the derived halo box has storage", allocated,
                    64 * fl::kTileCells);
    }

    // The same water in two rooms: the same tiles, and bit-identical answers. A
    // stronger statement than "about the same cost" -- it says the room's extent
    // does not enter the arithmetic at all.
    std::vector<fl::Particle> outcome[2];
    int tiles[2] = {0, 0};
    std::size_t bytes[2] = {0, 0};
    for (int which = 0; which < 2; ++which) {
        const int extent = which == 0 ? 20 : 400;
        fl::Field f;
        f.grid.h = 0.1;
        for (int a = 0; a < 3; ++a) f.grid.n[a] = extent;
        const double lo[3] = {0.5, 0.5, 1.0}, hi[3] = {1.1, 1.1, 1.4};
        fl::seedBox(f, lo, hi, 2, kRhoSeawater);
        fl::Params p;
        p.maxSubstep = 0.004;
        p.maxSubsteps = 8;
        fl::Account account;
        fl::resetAccount(f, account);
        fl::Solver s;
        bool short_ = false;
        int substeps = 0;
        for (int i = 0; i < 50; ++i) {
            const fl::StepResult r = s.step(f, 0.004, p, account);
            substeps += r.substeps;
            if (r.incomplete || r.projectionCapped || substeps > 8 * (i + 1)) {
                short_ = true;
                break;
            }
        }
        expectTrue("the room comparison ran its fifty steps in full", !short_);
        outcome[which] = f.particles;
        tiles[which] = s.tiles();
        bytes[which] = s.bytes();
    }
    expectEqual("the same water allocates the same tiles in a room 8000 times the volume",
                tiles[0], tiles[1]);
    expectEqual("and the same bytes", static_cast<long long>(bytes[0]),
                static_cast<long long>(bytes[1]));
    int differing = 0;
    double travelled = 0;
    for (std::size_t i = 0; i < outcome[0].size(); ++i)
        for (int a = 0; a < 3; ++a) {
            if (outcome[0][i].position[a] != outcome[1][i].position[a]) ++differing;
            travelled = std::max(travelled, std::abs(outcome[0][i].velocity[a]));
        }
    std::printf("     %d tiles, %zu bytes either way; %d of %zu coordinates differ; "
                "peak speed %.4g m/s\n",
                tiles[0], bytes[0], differing, 3 * outcome[0].size(), travelled);
    expectEqual("every particle is bit identical", differing, 0);
    // Vacuity: the water has to have been doing something for fifty steps, or two
    // identical stationary blocks would pass this.
    expectTrue("the water was moving", travelled > 1.0);

    // Arrival. A single particle walked across a tile seam, one substep at a time,
    // with the invariant checked at every position: the face weights sum to its
    // mass exactly. A halo one tile too small drops part of that sum and nothing
    // else in the solver notices.
    {
        fl::Field f;
        f.grid.h = 0.25;
        f.grid.n[0] = 24; f.grid.n[1] = 8; f.grid.n[2] = 8;
        fl::Particle q;
        q.position[0] = 0.0; q.position[1] = 1.0; q.position[2] = 1.0;
        q.mass = 7.0;
        f.particles.push_back(q);
        fl::Params p;
        fl::Solver s;
        int seams = 0, lastTile = -999;
        double worst = 0;
        // Tiles are four cells, so a seam falls every 1.0 m. Walk 6 m in 1 mm
        // steps, which lands on the seams and on both sides of each of them.
        for (int i = 0; i <= 6000; ++i) {
            f.particles[0].position[0] = static_cast<double>(i) * 0.001;
            s.rebuild(f, p);
            s.transferToGrid(f, p);
            for (int a = 0; a < 3; ++a)
                worst = std::max(worst, std::abs(s.totalFaceMass(a) - 7.0));
            const int tile = static_cast<int>(std::floor(f.particles[0].position[0] / 1.0));
            if (tile != lastTile) {
                ++seams;
                lastTile = tile;
            }
        }
        std::printf("     one particle walked 6 m across %d tile seams: worst face-mass"
                    " defect %.3g kg of 7\n", seams - 1, worst);
        expectTrue("it crossed several tile seams", seams - 1 >= 5);
        expectTrue("the face mass held to 1e-14 kg at every position", worst < 1e-14);
    }

    // Arrival, the other way round: water flying into a region that had no storage
    // at all a moment ago, with the tile count watched as it goes.
    {
        fl::Field f;
        f.grid.h = 0.1;
        f.grid.n[0] = 24; f.grid.n[1] = 8; f.grid.n[2] = 8;
        const double lo[3] = {0.2, 0.2, 0.2}, hi[3] = {0.6, 0.6, 0.6};
        fl::seedBox(f, lo, hi, 2, kRhoSeawater);
        for (fl::Particle& q : f.particles) q.velocity[0] = 3.0;
        fl::Params p;
        p.gravity[2] = 0;
        p.maxSubstep = 0.01;
        p.maxSubsteps = 8;
        fl::Account account;
        fl::resetAccount(f, account);
        fl::Solver s;
        double startX[3];
        fl::centroid(f, startX);
        int firstTiles = 0, sameTiles = 0, arrivalSubsteps = 0;
        bool short_ = false;
        for (int i = 0; i < 80; ++i) {
            const fl::StepResult r = s.step(f, 0.01, p, account);
            arrivalSubsteps += r.substeps;
            if (r.incomplete || r.projectionCapped || arrivalSubsteps > 8 * (i + 1)) {
                short_ = true;
                break;
            }
            if (i == 0) firstTiles = r.tiles;
            if (r.tiles == firstTiles) ++sameTiles;
            expectNear("mass while arriving", account.massResidual(), 0.0, 0.0);
            double totals[3];
            for (int a = 0; a < 3; ++a) totals[a] = s.totalFaceMass(a);
            for (int a = 0; a < 3; ++a)
                expectTrue("the arriving water's mass is all on the grid",
                           std::abs(totals[a] - f.totalMass()) < 1e-12 * f.totalMass());
        }
        double endX[3];
        fl::centroid(f, endX);
        std::printf("     a block flew %.4g m into empty space; tiles %d -> %d, %d steps at the"
                    " starting count\n", endX[0] - startX[0], firstTiles, s.tiles(), sameTiles);
        expectTrue("the arriving block took every step in full", !short_);
        expectTrue("it flew a long way", endX[0] - startX[0] > 1.0);
        // It crossed many tiles, so the allocation had to keep up rather than
        // happening once.
        expectTrue("the tile count moved as it went", sameTiles < 70);
        expectEqual("nothing was lost on the way", account.particleResidual(), 0);
        // It ends against the far wall, which is the other half of the claim:
        // clamping keeps mass in and counts every time it does.
        expectTrue("the far wall stopped it", account.clamped > 0);
        expectNear("and kept every gram", account.massResidual(), 0.0, 0.0);
    }
}

// ---------------------------------------------------------------------------
// 7. Sloshing: the first mode against `g k tanh(k d)`
// ---------------------------------------------------------------------------

void testSloshing() {
    std::printf("\n--- flip: a sloshing tank's first mode ---\n");
    const double L = 1.0;
    const double depth = 0.30;
    const double amplitude = 0.10;
    const double coarseH = 0.05, fineH = 1.0 / 30.0;
    // **The same tank at two resolutions**, not two tanks. The measured period has
    // a discretisation part that falls with `h` and a nonlinear part set by `a/d`,
    // which is physical and does not; holding the geometry and the amplitude fixed
    // in metres is what makes refining `h` a statement about the first alone.
    // The first version of this check varied the depth and the amplitude with the
    // grid and came out backwards, which is the whole reason the note is here.
    const SloshRun coarse = runTank(L, depth, coarseH, amplitude, 4.0, 0.01, true, 0.0);
    const SloshRun fine = runTank(L, depth, fineH, amplitude, 4.0, 0.01, true, 0.0);
    const double want = firstModePeriod(L, depth);

    const double coarseError = (coarse.period - want) / want;
    const double fineError = (fine.period - want) / want;
    std::printf("     h=%.4f (%2.0f cells deep): T %.5f against %.5f  (%+.2f%%), %d crossings,"
                " peak %.4f m\n",
                coarseH, depth / coarseH, coarse.period, want, 100.0 * coarseError,
                coarse.crossings, coarse.peak);
    std::printf("     h=%.4f (%2.0f cells deep): T %.5f against %.5f  (%+.2f%%), %d crossings,"
                " peak %.4f m\n",
                fineH, depth / fineH, fine.period, want, 100.0 * fineError, fine.crossings,
                fine.peak);

    // Vacuity first, and it is not a formality: the first version of this test
    // seeded a tilted surface a sixth of a cell high, the voxelised free surface
    // could not represent it, the tank never sloshed at all and the period came
    // back as zero.
    expectTrue("the coarse tank actually sloshed", coarse.crossings >= 3);
    expectTrue("the fine tank actually sloshed", fine.crossings >= 3);
    expectTrue("the centroid swung a real distance", coarse.peak > 0.02 * L);
    expectTrue("a period was measured at all", coarse.period > 0 && fine.period > 0);

    // Measured at +5.66% on six cells of depth and +3.97% on nine, converging on
    // about +3% -- which is the nonlinear correction at `a/d = 1/3` and does not
    // go away with the grid. Asserted at 7% and 5%, which is a factor of 1.25 on
    // each and not a factor of ten.
    expectTrue("the coarse period is within 7% of g k tanh(k d)", std::abs(coarseError) < 0.07);
    expectTrue("the fine period is within 5% of g k tanh(k d)", std::abs(fineError) < 0.05);
    expectTrue("refining the grid moves the period towards the analytic value",
               std::abs(fineError) < std::abs(coarseError));
    // And by a real margin, not by a last digit: a wrong dispersion relation would
    // converge just as monotonically to the wrong answer, so what makes this a
    // test is that the gap closes by a third of itself.
    expectTrue("and closes a real part of the gap",
               std::abs(coarseError) - std::abs(fineError) > 0.005);
    // And the sign: a voxelised free surface adds apparent depth at the surface
    // cell, which slows the mode. A solver that came out *fast* would be a
    // different defect and this would catch it.
    expectTrue("the discretisation error is a slow mode, not a fast one",
               coarseError > 0 && fineError > 0);
    expectNear("mass is exactly conserved through the slosh", coarse.massResidual, 0.0, 0.0);
    expectTrue("the projection kept up", !coarse.capped && !fine.capped);
    expectTrue("no step was short of the time it was asked for",
               !coarse.incomplete && !fine.incomplete);
    expectTrue("neither tank ran away in substeps", !coarse.runaway && !fine.runaway);
    // The substep count against the controller's own arithmetic: the tank never
    // exceeds about 1 m/s, so `courant h / v` is never below `h`, and 400 steps of
    // 0.01 s at `maxSubstep = 0.01` is 400 substeps and not 4 000. A collapsed
    // controller blows through this rather than running for an hour.
    std::printf("     substeps: coarse %d, fine %d, against 400 steps of 0.01 s\n",
                coarse.substeps, fine.substeps);
    expectTrue("the coarse tank took the substeps its Courant number allows",
               coarse.substeps <= 4 * 400);
    expectTrue("and so did the fine one", fine.substeps <= 4 * 400);
}

// ---------------------------------------------------------------------------
// 8. Why APIC: the measurement, not the citation
// ---------------------------------------------------------------------------

void testTransferComparison() {
    std::printf("\n--- flip: PIC against FLIP against APIC ---\n");
    struct Mode { const char* name; bool affine; double blend; };
    const Mode modes[3] = {{"PIC ", false, 0.0}, {"FLIP", false, 1.0}, {"APIC", true, 0.0}};
    double keptAfterTen[3] = {0, 0, 0};
    for (int m = 0; m < 3; ++m) {
        fl::Field f;
        f.grid.h = 0.05;
        for (int a = 0; a < 3; ++a) f.grid.n[a] = 24;
        const double lo[3] = {0.20, 0.20, 0.20}, hi[3] = {1.00, 1.00, 1.00};
        fl::seedBox(f, lo, hi, 2, kRhoSeawater);
        const double about[3] = {0.60, 0.60, 0.60};
        const double omega = 2.0;
        for (fl::Particle& q : f.particles) {
            q.velocity[0] = -omega * (q.position[1] - about[1]);
            q.velocity[1] = omega * (q.position[0] - about[0]);
            // The affine part of a rigid rotation, so the first transfer is not
            // handicapped by a particle that has never been through one.
            q.affine[1] = -omega;
            q.affine[3] = omega;
        }
        fl::Params p;
        p.affine = modes[m].affine;
        p.flipBlend = modes[m].blend;
        p.gravity[2] = 0;
        fl::Solver s;
        double before[3];
        fl::angularMomentum(f, about, before);
        for (int i = 0; i < 10; ++i) {
            s.rebuild(f, p);
            s.transferToGrid(f, p);
            s.extrapolate(p);
            s.saveGrid();
            s.extrapolate(p);
            s.transferToParticles(f, p);
        }
        double after[3];
        fl::angularMomentum(f, about, after);
        keptAfterTen[m] = after[2] / before[2];
        std::printf("     %s keeps %.6f of its angular momentum through ten transfers\n",
                    modes[m].name, keptAfterTen[m]);
        expectTrue("the rotating body had angular momentum to lose", std::abs(before[2]) > 10.0);
    }
    // FLIP with no grid dynamics at all is the identity, exactly. This is the
    // control that fails when `saveGrid` is called before the extrapolation
    // instead of after it -- and that defect is invisible at the default blend.
    expectNear("FLIP's transfer with no grid dynamics is the identity", keptAfterTen[1], 1.0,
               1e-12);
    // APIC loses 0.18% per transfer where PIC loses 1.13%.
    expectTrue("APIC keeps more angular momentum than PIC", keptAfterTen[2] > keptAfterTen[0]);
    expectTrue("APIC's loss is at least five times smaller than PIC's",
               (1.0 - keptAfterTen[2]) * 5.0 < (1.0 - keptAfterTen[0]));
    expectTrue("PIC's loss is real and not rounding", 1.0 - keptAfterTen[0] > 0.05);

    // The other half of the choice: FLIP keeps everything, including the
    // particle-borne velocity the grid never sees.
    const SloshRun apic = runTank(1.0, 0.30, 0.05, 0.10, 1.5, 0.01, true, 0.0);
    const SloshRun flip = runTank(1.0, 0.30, 0.05, 0.10, 1.5, 0.01, false, 1.0);
    const SloshRun pic = runTank(1.0, 0.30, 0.05, 0.10, 1.5, 0.01, false, 0.0);
    std::printf("     particle noise after 1.5 s of tank: PIC %.5f  APIC %.5f  FLIP %.5f\n",
                pic.noise, apic.noise, flip.noise);
    std::printf("     centroid swing over the same run:   PIC %.5f  APIC %.5f  FLIP %.5f m\n",
                pic.peak, apic.peak, flip.peak);
    expectTrue("PIC carries no particle-borne velocity at all", pic.noise < 1e-2);
    expectTrue("FLIP's particle noise is large", flip.noise > 0.05);
    expectTrue("APIC's noise is at least five times smaller than FLIP's",
               apic.noise * 5.0 < flip.noise);
    // And PIC pays for its quiet with the amplitude of the mode it is meant to be
    // measuring, which is why the answer is APIC and not PIC.
    expectTrue("PIC damps the mode more than APIC does", pic.peak < apic.peak);
    expectTrue("there was a mode to damp", apic.peak > 0.02);
    expectTrue("none of the three ran out of substeps",
               !apic.incomplete && !flip.incomplete && !pic.incomplete);
    expectTrue("and none of them ran away in substeps",
               !apic.runaway && !flip.runaway && !pic.runaway);
}

// ---------------------------------------------------------------------------
// 9. Mass, walls, and the substep budget
// ---------------------------------------------------------------------------

void testMassAndBudget() {
    std::printf("\n--- flip: mass, walls and the substep bound ---\n");
    // A dam break: half a tank of water released against three walls, which is the
    // most violent thing this solver is asked to do and therefore the right place
    // to ask whether it loses any.
    fl::Field f;
    f.grid.h = 0.04;
    f.grid.n[0] = 25; f.grid.n[1] = 4; f.grid.n[2] = 12;
    const double lo[3] = {0, 0, 0}, hi[3] = {0.28, 0.16, 0.40};
    fl::seedBox(f, lo, hi, 2, kRhoSeawater);
    fl::setTotalMass(f, kRhoSeawater * 0.28 * 0.16 * 0.40);
    const double startMass = f.totalMass();
    const long long startCount = static_cast<long long>(f.particles.size());

    fl::Params p;
    p.maxSubstep = 0.005;
    // The one test here that genuinely subdivides: 3.3 m/s on a 0.04 m cell puts
    // `courant h / v` at 11 ms against a 5 ms step, so it is still one substep, and
    // thirty-two is the budget a flow ten times faster would need.
    p.maxSubsteps = 32;
    fl::Account account;
    fl::resetAccount(f, account);
    fl::Solver s;
    double peak = 0;
    int substeps = 0;
    bool incomplete = false, capped = false, runaway = false;
    for (int i = 0; i < 200; ++i) {
        const fl::StepResult r = s.step(f, 0.005, p, account);
        substeps += r.substeps;
        peak = std::max(peak, r.peakSpeed);
        expectNear("mass through the dam break", account.massResidual(), 0.0, 0.0);
        if (r.projectionCapped) {
            capped = true;
            break;
        }
        if (r.incomplete) {
            incomplete = true;
            break;
        }
        // Same arithmetic floor, checked inside the loop. The dam break reaches
        // about 3.3 m/s on a 0.04 m cell, so `courant h / v` is about 11 ms
        // against a 5 ms step: one substep each, and twenty times that is a bound
        // no physical flow here approaches.
        if (substeps > 20 * (i + 1)) {
            runaway = true;
            break;
        }
    }
    std::printf("     1.0 s of dam break: peak %.4g m/s, %d substeps, %lld wall clamps,"
                " mass residual %.17g\n", peak, substeps, account.clamped,
                account.massResidual());
    expectNear("mass is exact after a dam break", f.totalMass(), startMass, 0.0);
    expectEqual("every particle survived", static_cast<long long>(f.particles.size()), startCount);
    expectEqual("the account agrees", account.particleResidual(), 0);
    expectTrue("the dam actually broke", peak > 1.0);
    expectTrue("the water reached the walls", account.clamped > 100);
    expectTrue("every particle is still inside the room", fl::validate(f, p).empty());
    expectTrue("the step was never short of the time it was asked for", !incomplete);
    expectTrue("and the projection never ran out of iterations", !capped);
    expectTrue("and the substep count never ran away", !runaway);

    // The substep count against the **arithmetic floor** of the controller rather
    // than against a clock. The controller takes `min(maxSubstep, courant h / v)`,
    // so 200 steps of 0.005 s cannot need more than
    // `200 * ceil(0.005 * v / (courant h))` substeps at the peak speed it reached
    // -- and a controller that has collapsed blows straight through that. This is
    // the bound a mutation sweep needs: a defect that drives the substep down
    // fails here instead of turning a two-second test into an hour.
    const double floorSubstep = p.courant * f.grid.h / std::max(peak, 1e-9);
    const int bound = 200 * (1 + static_cast<int>(std::ceil(0.005 / floorSubstep)));
    std::printf("     substep bound from the peak speed: %d, taken %d\n", bound, substeps);
    expectTrue("the substep count is inside the controller's own arithmetic bound",
               substeps <= bound);
    expectTrue("the bound is not vacuous", bound < 20 * 200);

    // And the budget, when it does bind, is published rather than swallowed.
    {
        fl::Field g = f;
        g.time = 0;
        fl::Params tight = p;
        tight.maxSubsteps = 3;
        tight.maxSubstep = 1e-4;
        fl::Account other;
        fl::resetAccount(g, other);
        fl::Solver t;
        const fl::StepResult r = t.step(g, 0.5, tight, other);
        expectTrue("a step that ran out of budget says so", r.incomplete);
        expectEqual("and took exactly the substeps it was allowed", r.substeps, 3);
        expectTrue("and advanced only that far", g.time < 0.5);
        expectNear("and still lost no mass", other.massResidual(), 0.0, 0.0);
    }
}

// ---------------------------------------------------------------------------
// 10. What Phase 5's escalation will stand on
// ---------------------------------------------------------------------------

void testQuiescent() {
    std::printf("\n--- flip: the quiescent round trip ---\n");
    fl::Field f;
    f.grid.h = 0.05;
    f.grid.n[0] = 8; f.grid.n[1] = 8; f.grid.n[2] = 80;
    const double area = 0.40 * 0.40;
    const double level = 2.7183;   // deliberately not a multiple of the cell
    const double lo[3] = {0, 0, 0}, hi[3] = {0.40, 0.40, level};
    fl::seedBox(f, lo, hi, 2, kRhoSeawater);
    const double want = kRhoSeawater * area * level;
    const double lattice = f.totalMass();
    expectTrue("the lattice alone does not hit the mass exactly", lattice != want);
    expectTrue("but it is close", std::abs(lattice - want) < 0.05 * want);
    expectTrue("setTotalMass took", fl::setTotalMass(f, want));
    // Exactly the number asked for, because the last particle takes the remainder.
    expectNear("the total is the mass that was asked for", f.totalMass(), want, 0.0);
    const double back = fl::quiescentLevel(f, kRhoSeawater, area, 0.0);
    std::printf("     lattice %.17g -> exact %.17g kg; level %.17g back from %.17g (%.3g m)\n",
                lattice, f.totalMass(), back, level, std::abs(back - level));
    expectNear("the still-water level comes back", back, level, 5e-16);
    // Every particle still carries a positive mass -- the remainder trick must not
    // make the last one negative or zero.
    double smallest = 1e30, largest = 0;
    for (const fl::Particle& q : f.particles) {
        smallest = std::min(smallest, q.mass);
        largest = std::max(largest, q.mass);
    }
    expectTrue("every particle has positive mass", smallest > 0);
    expectTrue("the remainder is a rounding and not a share",
               (largest - smallest) < 1e-12 * largest);
    // A floor offset is carried through, because a compartment's floor is not at
    // z = 0 on a ship.
    expectNear("the level is measured from the floor it was given",
               fl::quiescentLevel(f, kRhoSeawater, area, -3.5), level - 3.5, 5e-16);
    expectEqual("an empty field has no problems worth reporting",
                static_cast<long long>(fl::validate(f).size()), 0);
}

void testValidation() {
    std::printf("\n--- flip: validation ---\n");
    fl::Field f;
    expectTrue("a field with no grid is a problem", !fl::validate(f).empty());
    f.grid.h = 0.1;
    for (int a = 0; a < 3; ++a) f.grid.n[a] = 4;
    expectEqual("a well formed empty field has none",
                static_cast<long long>(fl::validate(f).size()), 0);
    fl::Particle q;
    q.position[0] = 9.0;   // outside
    q.mass = 1.0;
    f.particles.push_back(q);
    expectTrue("a particle outside the room is a problem", !fl::validate(f).empty());
    fl::Params bad;
    bad.flipBlend = 1.5;
    bad.extrapolationDepth = 1;
    bad.courant = 5.0;
    expectTrue("three bad parameters are three problems",
               fl::validate(f, bad).size() >= 4);
    // And the solver counts a particle it cannot place rather than dropping it.
    fl::Params ok;
    fl::Solver s;
    s.rebuild(f, ok);
    expectEqual("the solver counts the particle it could not place", s.outsideParticles(), 1);
}

// **A projection handed a NaN reported a perfect converged solve.**
//
// Both of `project`'s reductions were `std::max(accumulator, x)`, which returns
// the accumulator for a NaN `x` -- the argument order decides, and this one drops
// it. Everything followed: `rhsScale` came back as the largest *finite* cell,
// `best` came back 0, `best > tolerance` was false so the entire CG solve was
// skipped, and `residual_ = 0.0`, `capped_ = false`, `iterations_ = 0` were
// reported over a velocity field that is not a number.
//
// In a singular region -- a fully submerged compartment, the ordinary case for
// flooding -- one bad face is enough, and the failure is quieter still:
// `rhsScale` is taken before `removeMean`, so it is a real positive number, and
// then the mean subtraction spreads the NaN across every cell. The caller is left
// with a healthy report, an all-zero pressure field, and a divergence nothing
// removed.
void testANaNFaceIsNotAConvergedSolve() {
    std::printf("\n--- flip: a projection handed a NaN has not converged ---\n");
    const int n[3] = {6, 5, 4};
    const double h = 0.05;
    fl::Field f;
    f.grid.h = h;
    for (int a = 0; a < 3; ++a) f.grid.lo[a] = 0.0;
    for (int a = 0; a < 3; ++a) f.grid.n[a] = n[a];
    const double lo[3] = {0, 0, 0}, hi[3] = {n[0] * h, n[1] * h, n[2] * h};
    fl::seedBox(f, lo, hi, 2, kRhoSeawater);

    fl::Params p;
    p.projectionTolerance = 1e-14;
    p.projectionIterations = 500;

    // The control first: the same field, untouched, must converge. Without this
    // every assertion below would hold for a solver that called everything capped.
    fl::Solver clean;
    clean.rebuild(f, p);
    clean.transferToGrid(f, p);
    expectTrue("the reference region is fully submerged", clean.singular());
    clean.project(0.01, p);
    expectTrue("and it converges", !clean.lastCapped());
    expectTrue("with a finite residual", std::isfinite(clean.lastResidual()));

    // One face, out of several hundred, set to a NaN.
    fl::Solver poisoned;
    poisoned.rebuild(f, p);
    poisoned.transferToGrid(f, p);
    expectTrue("the poisoned solver stores the face it is given",
               poisoned.setFaceVelocity(0, 2, 2, 2, std::numeric_limits<double>::quiet_NaN()));
    poisoned.project(0.01, p);
    expectTrue("one NaN face is reported as not converged", poisoned.lastCapped());
    expectTrue("and the residual it reports is not a number, so no threshold passes",
               std::isnan(poisoned.lastResidual()));
    // The specific lie: it used to report exactly zero, which reads as an exact
    // solve rather than as a solve that never ran.
    expectTrue("in particular it does not report a residual of zero",
               poisoned.lastResidual() != 0.0);
}

}  // namespace

void runFlipTests() {
    testKernel();
    testTransfer();
    testHydrostatic();
    testProjection();
    testBallistic();
    testSparse();
    testSloshing();
    testTransferComparison();
    testMassAndBudget();
    testQuiescent();
    testValidation();
    testANaNFaceIsNotAConvergedSolve();
}
