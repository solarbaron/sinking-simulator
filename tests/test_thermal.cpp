// SPDX-License-Identifier: MIT
//
// Validation of the implicit conduction solve.
//
// Conduction is unusually well supplied with exact answers, so almost nothing
// here is a tolerance on an eyeballed number. Seven things are asserted as
// identities or as closed forms:
//
//   * **A linear temperature field is reproduced exactly.** The conduction patch
//     test: prescribe `T = a + b.x` on the boundary of a randomly distorted block
//     and every interior node comes back to rounding. It is run on an *oblique*
//     gradient and on a mesh containing a collapsed element, because a gradient
//     down an axis is reproduced by an operator with two of its three Jacobian
//     columns wrong.
//   * **`T^T K T = k |grad T|^2 V` for that linear field**, on a distorted
//     element, which is the closed form that ties the Gauss weights to the
//     element's actual volume. A wrong weight cancels out of the patch test and
//     does not cancel out of this.
//   * **`K 1 = 0` and `sum_a grad N_a = 0`.** Partition of unity. These are what
//     make the energy account close, so they are asserted directly rather than
//     inferred from it.
//   * **The semi-infinite solid**, `T = Ts + (T0 - Ts) erf(x / 2 sqrt(alpha t))`.
//     Exact, and the only test here that exercises the time integration against
//     something outside the code.
//   * **The energy account closes to machine precision**, not to the
//     integrator's order -- see `thermal.hpp`'s `Account`. Asserted at 1e-12 of
//     the enthalpy moved, with each of the four heat channels separately checked
//     to be carrying a real share. With temperature-dependent properties it
//     closes to the *Picard tolerance* instead, and that is asserted as a
//     scaling -- tighten the tolerance a thousandfold and the residual follows,
//     which a wrong quadrature would not.
//   * **The orders are asserted, not the tolerances.** Backward Euler is first
//     order in time and the element second in space, and a scheme that quietly
//     came out at the wrong order is a defect a loose tolerance hides. Both are
//     measured by refinement and the *slope* is what is asserted.
//   * **The explicit stability limit is `h_min^2 / 2 alpha` exactly**, and the
//     Kirchhoff transform gives the *nonlinear* steady operator a closed form of
//     its own. Two more places where an exact answer was available and taken.
//
// Two of these caught real defects while being written, both recorded where they
// happened: a power iteration whose start vector was orthogonal by symmetry to
// the eigenvalue it was hunting, and a Simpson rule in this file that sampled a
// material curve exactly at one of its own discontinuities.
//
// The vacuity guards matter as much as the assertions, because a conduction test
// in which nothing conducts passes trivially. Every transient here checks that
// the front is inside the domain and has not reached the far end; every steady
// one checks the gradient is not zero; the nonlinear ones check that the
// temperature dependence actually moved the answer.
#include "engine/sim/thermal.hpp"
#include "engine/sim/buckling.hpp"
#include "engine/sim/collapse.hpp"
#include "engine/sim/constraint.hpp"
#include "engine/sim/plasticity.hpp"
#include "engine/sim/promotion.hpp"
#include "engine/sim/scantlings.hpp"
#include "engine/sim/solid_shell.hpp"
#include "game/prototype/ferry.hpp"
#include "harness.hpp"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <numeric>
#include <tuple>
#include <random>
#include <string>
#include <vector>

using testing::expectEqual;
using testing::expectNear;
using testing::expectTrue;

namespace {

using sim::StructuralMaterial;
using sim::Vec3;
using sim::ah36Steel;
namespace ss = sim::solidshell;
namespace th = sim::thermal;

constexpr double kC = th::kCelsius;

// A material with round numbers, so a closed form can be checked against
// arithmetic done by hand rather than against another floating-point pipeline.
StructuralMaterial roundSteel() {
    StructuralMaterial m = ah36Steel();
    m.density = 8000.0;
    m.conductivity = 50.0;
    m.specificHeat = 500.0;  // alpha = 50 / (8000 * 500) = 1.25e-5 m^2/s exactly
    return m;
}
constexpr double kRoundAlpha = 1.25e-5;

// A bar `length` long by `w` by `w`, `n` elements along x and one across each
// other direction. `grade` stretches the spacing geometrically, so the elements
// are all different sizes and no error can cancel between a pair of them.
ss::HexMesh makeBar(double length, double w, int n, double grade = 1.0) {
    ss::HexMesh mesh;
    std::vector<double> x(static_cast<std::size_t>(n) + 1, 0.0);
    double step = 1.0, total = 0.0;
    for (int i = 0; i < n; ++i) {
        total += step;
        step *= grade;
    }
    step = 1.0;
    for (int i = 0; i < n; ++i) {
        x[static_cast<std::size_t>(i) + 1] = x[static_cast<std::size_t>(i)] + length * step / total;
        step *= grade;
    }
    const auto at = [&](int i, int j, int k) {
        return static_cast<std::uint32_t>((i * 2 + j) * 2 + k);
    };
    mesh.position.resize((static_cast<std::size_t>(n) + 1) * 4 * 3);
    for (int i = 0; i <= n; ++i)
        for (int j = 0; j < 2; ++j)
            for (int k = 0; k < 2; ++k) {
                double* p = &mesh.position[at(i, j, k) * 3];
                p[0] = x[static_cast<std::size_t>(i)];
                p[1] = j * w;
                p[2] = k * w;
            }
    for (int i = 0; i < n; ++i) {
        const std::uint32_t v[8] = {at(i, 0, 0),     at(i + 1, 0, 0), at(i + 1, 1, 0),
                                    at(i, 1, 0),     at(i, 0, 1),     at(i + 1, 0, 1),
                                    at(i + 1, 1, 1), at(i, 1, 1)};
        mesh.index.insert(mesh.index.end(), v, v + 8);
    }
    mesh.fixed.assign(mesh.nodeCount() * 3, 0u);
    mesh.prescribed.assign(mesh.nodeCount() * 3, 0.0);
    return mesh;
}

// A structured block, nx by ny by nz elements over `size`. Interior nodes are
// displaced by up to `jitter` of the local spacing, which makes every element a
// different distorted hexahedron -- the asymmetric fixture a symmetric error
// hides in.
ss::HexMesh makeBlock(Vec3 size, int nx, int ny, int nz, double jitter, unsigned seed = 7u) {
    ss::HexMesh mesh;
    const int sx = nx + 1, sy = ny + 1, sz = nz + 1;
    const auto at = [&](int i, int j, int k) {
        return static_cast<std::uint32_t>((i * sy + j) * sz + k);
    };
    mesh.position.resize(static_cast<std::size_t>(sx) * sy * sz * 3);
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> u(-1.0, 1.0);
    for (int i = 0; i < sx; ++i)
        for (int j = 0; j < sy; ++j)
            for (int k = 0; k < sz; ++k) {
                double* p = &mesh.position[at(i, j, k) * 3];
                p[0] = size.x * i / nx;
                p[1] = size.y * j / ny;
                p[2] = size.z * k / nz;
                const bool interior = i > 0 && i < nx && j > 0 && j < ny && k > 0 && k < nz;
                if (!interior || jitter == 0.0) continue;
                p[0] += jitter * u(rng) * size.x / nx;
                p[1] += jitter * u(rng) * size.y / ny;
                p[2] += jitter * u(rng) * size.z / nz;
            }
    for (int i = 0; i < nx; ++i)
        for (int j = 0; j < ny; ++j)
            for (int k = 0; k < nz; ++k) {
                const std::uint32_t v[8] = {
                    at(i, j, k),         at(i + 1, j, k),         at(i + 1, j + 1, k),
                    at(i, j + 1, k),     at(i, j, k + 1),         at(i + 1, j, k + 1),
                    at(i + 1, j + 1, k + 1), at(i, j + 1, k + 1)};
                mesh.index.insert(mesh.index.end(), v, v + 8);
            }
    mesh.fixed.assign(mesh.nodeCount() * 3, 0u);
    mesh.prescribed.assign(mesh.nodeCount() * 3, 0.0);
    return mesh;
}

Vec3 nodeAt(const ss::HexMesh& mesh, std::size_t n) {
    return {mesh.position[n * 3], mesh.position[n * 3 + 1], mesh.position[n * 3 + 2]};
}

// Everything a `Problem` needs to hold a set of nodes at a temperature.
void hold(th::Problem& problem, const ss::HexMesh& mesh, std::size_t node, double kelvin) {
    if (problem.prescribed.empty()) {
        problem.prescribed.assign(mesh.nodeCount(), 0u);
        problem.prescribedValue.assign(mesh.nodeCount(), 0.0);
    }
    problem.prescribed[node] = 1u;
    problem.prescribedValue[node] = kelvin;
}

double elementVolume(const ss::HexMesh& mesh, std::size_t e) {
    double nodes[ss::kDof];
    mesh.gather(e, mesh.position, nodes);
    double w[ss::kGauss];
    ss::gaussVolumes(nodes, w);
    return std::accumulate(w, w + ss::kGauss, 0.0);
}

double meshVolume(const ss::HexMesh& mesh) {
    double v = 0.0;
    for (std::size_t e = 0; e < mesh.elementCount(); ++e) v += elementVolume(mesh, e);
    return v;
}

// The slope of log(error) against log(refinement), by least squares over the
// whole sequence rather than from the two finest points -- a two-point order is
// one rounding away from being an artefact.
double orderOf(const std::vector<double>& error, double ratio) {
    const std::size_t n = error.size();
    double sx = 0, sy = 0, sxx = 0, sxy = 0;
    for (std::size_t i = 0; i < n; ++i) {
        const double x = static_cast<double>(i) * std::log(ratio);
        const double y = std::log(error[i]);
        sx += x;
        sy += y;
        sxx += x * x;
        sxy += x * y;
    }
    const double d = static_cast<double>(n) * sxx - sx * sx;
    return -(static_cast<double>(n) * sxy - sx * sy) / d;
}

// --- 1. The material curves ------------------------------------------------------

void testMaterialCurves() {
    const StructuralMaterial steel = ah36Steel();

    // The constants on `StructuralMaterial` and the curves in `thermal.cpp` are two
    // statements of the same thing, so they are tied together rather than each
    // asserted against a transcribed number. Drifting apart is exactly the failure
    // `CLAUDE.md` records for the three documents that quoted each other.
    expectNear("k(20 C) ties to StructuralMaterial", th::carbonSteelConductivity(kC + 20.0),
               steel.conductivity, 1e-12);
    expectNear("c(20 C) ties to StructuralMaterial", th::carbonSteelSpecificHeat(kC + 20.0),
               steel.specificHeat, 1e-9);

    // EN 1993-1-2 3.4.1.3, evaluated by hand: 54 - 0.0333*600.
    expectNear("k(600 C)", th::carbonSteelConductivity(kC + 600.0), 34.02, 1e-10);
    // **The 800 C branch point is not reachable, and that is worth asserting
    // rather than discovering twice.** `kelvin - kCelsius` never lands on exactly
    // 800.0 for any double `kelvin`: the nearest values either side are
    // 799.9999999999999 and 800.0000000000001, and the subtraction steps straight
    // over. So `t >= 800` and `t > 800` are the same function, and no test can
    // tell them apart -- which is what a mutation run reported, correctly, as an
    // undetected change.
    //
    // 600 C and 735 C *are* representable, so their comparisons do decide
    // something and are tested on both sides below. Recording which boundaries are
    // reachable is the same lesson `CLAUDE.md` draws from `sectionElements`: two
    // pieces of code that decide one thing must agree about what is on the edge.
    bool reachable = false;
    for (double k = std::nextafter(kC + 800.0, 0.0), i = 0; i < 8; k = std::nextafter(k, 2000.0), ++i)
        if (k - kC == 800.0) reachable = true;
    expectTrue("no double sits exactly on the 800 C branch point", !reachable);
    expectTrue("600 C does, and 735 C does",
               (kC + 600.0) - kC == 600.0 && (kC + 735.0) - kC == 735.0);
    expectNear("k just above 800 C is the flat branch",
               th::carbonSteelConductivity(kC + 800.0), 27.3, 1e-12);
    expectNear("k just below it is the linear one",
               th::carbonSteelConductivity(std::nextafter(kC + 800.0, 0.0)), 27.36, 1e-9);
    // The claim `thermal.hpp` makes -- and a guard against a curve that has gone
    // flat, which would make every temperature-dependent test below vacuous.
    const double fall = 1.0 - th::carbonSteelConductivity(kC + 600.0) / steel.conductivity;
    expectTrue("conductivity falls by more than 30% at 600 C", fall > 0.30 && fall < 0.40);

    // 3.4.1.2's two hyperbolae both reach 5000 at 735 C. That they meet is the
    // standard's, not this file's, and it is asserted because the transition is the
    // whole reason a secant capacity is used.
    expectNear("c(735 C) from below", th::carbonSteelSpecificHeat(kC + 735.0 - 1e-7), 5000.0,
               1e-3);
    expectNear("c(735 C) from above", th::carbonSteelSpecificHeat(kC + 735.0), 5000.0, 1e-9);
    expectTrue("c(735 C) is more than ten times c(20 C)",
               th::carbonSteelSpecificHeat(kC + 735.0) > 10.0 * steel.specificHeat);

    // **`h` is the integral of `c`.** Against composite Simpson over each smooth
    // piece, not against a finite difference of `h` -- a central difference near
    // the 735 C transition carries a truncation error of its own that is *larger*
    // than the defect it would be looking for, so it can only ever assert a loose
    // tolerance. Quadrature is the operation being checked, so quadrature is what
    // it is checked against, and to nine figures.
    //
    // The piece boundaries are integrated across explicitly, which is also the
    // check that `h`'s constants of integration are the ones that make it
    // continuous: an integral assembled from four independent antiderivatives
    // would agree with this on each piece and disagree cumulatively.
    expectNear("h(20 C) is the datum", th::carbonSteelEnthalpy(kC + 20.0), 0.0, 1e-9);
    // The endpoints are sampled a nanokelvin *inside* the interval. `c` is
    // discontinuous at 600, 735 and 900 C -- the standard's own pieces step by a
    // few tenths -- so sampling the endpoint itself takes the wrong branch and puts
    // a `(jump) * d/3` error into the quadrature. That is a defect in the check and
    // not in `h`: it showed up as a *constant* 2.875e-4 J/kg offset appearing at
    // exactly 600 C and never growing, which is the signature of one bad sample and
    // not of a wrong antiderivative.
    const auto simpson = [](double a, double b) {
        const int panels = 20000;  // even
        const double d = (b - a) / panels;
        const double inside = 1e-9;
        double sum = th::carbonSteelSpecificHeat(kC + a + inside) +
                     th::carbonSteelSpecificHeat(kC + b - inside);
        for (int i = 1; i < panels; ++i)
            sum += ((i & 1) ? 4.0 : 2.0) *
                   th::carbonSteelSpecificHeat(kC + a + static_cast<double>(i) * d);
        return sum * d / 3.0;
    };
    const double kink[3] = {600.0, 735.0, 900.0};
    double worst = 0.0, worstAt = 0.0;
    for (double target : {150.0, 400.0, 599.9, 650.0, 700.0, 734.9, 760.0, 850.0, 899.9, 1000.0,
                          1200.0}) {
        double integral = 0.0, from = 20.0;
        for (double k : kink)
            if (target > k) {
                integral += simpson(from, k);
                from = k;
            }
        integral += simpson(from, target);
        const double h = th::carbonSteelEnthalpy(kC + target);
        const double relative = std::abs(h - integral) / std::abs(integral);
        if (relative > worst) {
            worst = relative;
            worstAt = target;
        }
    }
    std::printf("     h against composite Simpson of c: worst relative error %.3e at %.0f C\n",
                worst, worstAt);
    expectTrue("h is the integral of c", worst < 1e-9);

    // Monotone, so the secant capacity is never negative -- which would turn an
    // unconditionally stable scheme into one that is not.
    double last = th::carbonSteelEnthalpy(0.0);
    bool monotone = true;
    for (double k = 5.0; k <= 1600.0; k += 1.0) {
        const double h = th::carbonSteelEnthalpy(k);
        if (!(h > last)) monotone = false;
        last = h;
    }
    expectTrue("h is strictly increasing over 0..1600 K", monotone);

    // Continuity at the three junctions: the standard's own pieces step by a few
    // tenths in `c`, and `h` must not step at all.
    for (double t : {600.0, 735.0, 900.0}) {
        const double jump = th::carbonSteelEnthalpy(kC + t + 1e-7) -
                            th::carbonSteelEnthalpy(kC + t - 1e-7);
        expectTrue("h is continuous at " + std::to_string(static_cast<int>(t)) + " C",
                   std::abs(jump) < 1e-3);
    }

    // The point of the secant, quantified: a 30 K step across the transition
    // carries far more enthalpy than either endpoint's `c` would charge for.
    const double secant =
        (th::carbonSteelEnthalpy(kC + 750.0) - th::carbonSteelEnthalpy(kC + 720.0)) / 30.0;
    const double endpoints =
        0.5 * (th::carbonSteelSpecificHeat(kC + 720.0) + th::carbonSteelSpecificHeat(kC + 750.0));
    std::printf("     secant c over 720..750 C = %.0f J/(kg K), endpoint mean %.0f\n", secant,
                endpoints);
    expectTrue("the secant across the spike exceeds the endpoint mean", secant > 1.5 * endpoints);
}

// --- 2. The element -----------------------------------------------------------------

void testElement() {
    // A thoroughly distorted single hexahedron. Nothing about it is symmetric, so
    // nothing cancels.
    const double nodes[ss::kDof] = {0.03, -0.02, 0.01,  1.11, 0.05,  -0.04, 1.02, 0.93,
                                    0.07, 0.08,  1.05,  0.02, 0.02,  0.06,  0.61, 0.97,
                                    0.03, 0.55,  1.06,  1.02, 0.66,  0.04,  0.95, 0.58};
    th::Forms forms;
    expectTrue("a distorted element integrates", th::computeForms(nodes, forms));

    // Partition of unity, both of it. `sum_a grad N_a = 0` is what makes `K 1 = 0`
    // and therefore what makes the energy account close.
    double worstShape = 0.0, worstGradient = 0.0, scale = 0.0;
    for (int gp = 0; gp < ss::kGauss; ++gp) {
        double s = 0.0, g[3] = {0, 0, 0};
        for (int a = 0; a < ss::kNodes; ++a) {
            s += forms.shape[gp][a];
            for (int i = 0; i < 3; ++i) {
                g[i] += forms.gradient[gp][a][i];
                scale = std::max(scale, std::abs(forms.gradient[gp][a][i]));
            }
        }
        worstShape = std::max(worstShape, std::abs(s - 1.0));
        for (int i = 0; i < 3; ++i) worstGradient = std::max(worstGradient, std::abs(g[i]));
    }
    expectTrue("sum N_a = 1", worstShape < 1e-15);
    expectTrue("sum grad N_a = 0", worstGradient < 1e-14 * scale);

    // The Gauss weights are the same object `solidshell::gaussVolumes` publishes,
    // and it is validated. Comparing them is what says this file's Jacobian is the
    // element's and not merely something plausible.
    double volumes[ss::kGauss];
    ss::gaussVolumes(nodes, volumes);
    double volume = 0.0, worstWeight = 0.0;
    for (int gp = 0; gp < ss::kGauss; ++gp) {
        volume += volumes[gp];
        worstWeight = std::max(worstWeight, std::abs(forms.weight[gp] - volumes[gp]));
    }
    expectTrue("Gauss weights match solidshell::gaussVolumes", worstWeight < 1e-16 * volume);
    expectTrue("the element has a volume", volume > 0.3 && volume < 3.0);

    double k[ss::kGauss], rc[ss::kGauss];
    for (int gp = 0; gp < ss::kGauss; ++gp) {
        k[gp] = 50.0;
        rc[gp] = 4.0e6;
    }
    double ke[64], ce[64], cl[64];
    th::conductance(forms, k, ke);
    th::capacity(forms, rc, false, ce);
    th::capacity(forms, rc, true, cl);

    double asymmetry = 0.0, magnitude = 0.0;
    for (int a = 0; a < 8; ++a)
        for (int b = 0; b < 8; ++b) {
            asymmetry = std::max(asymmetry, std::abs(ke[a * 8 + b] - ke[b * 8 + a]));
            magnitude = std::max(magnitude, std::abs(ke[a * 8 + b]));
        }
    expectTrue("K is symmetric", asymmetry < 1e-14 * magnitude);

    double rowMax = 0.0;
    for (int a = 0; a < 8; ++a) {
        double s = 0.0;
        for (int b = 0; b < 8; ++b) s += ke[a * 8 + b];
        rowMax = std::max(rowMax, std::abs(s));
    }
    expectTrue("K 1 = 0", rowMax < 1e-13 * magnitude);

    // `1^T C 1 = rho c V`, exactly. This is where a wrong Gauss weight shows: it
    // cancels out of every field the element reproduces exactly and does not cancel
    // out of the total heat capacity.
    double totalC = 0.0, totalL = 0.0;
    for (int a = 0; a < 8; ++a)
        for (int b = 0; b < 8; ++b) {
            totalC += ce[a * 8 + b];
            totalL += cl[a * 8 + b];
        }
    expectNear("1^T C 1 = rho c V", totalC, 4.0e6 * volume, 1e-9 * 4.0e6 * volume);
    expectNear("lumping preserves the total capacity", totalL, totalC, 1e-9 * totalC);

    // **Row-sum lumping leaves `1^T C dT` exactly unchanged** -- the property the
    // energy account depends on, asserted against an asymmetric `dT` so that a
    // uniform one cannot hide it.
    const double dT[8] = {13.0, -41.0, 7.5, 202.0, -3.25, 88.0, 0.5, -160.0};
    double consistent = 0.0, lumped = 0.0;
    for (int a = 0; a < 8; ++a)
        for (int b = 0; b < 8; ++b) {
            consistent += ce[a * 8 + b] * dT[b];
            lumped += cl[a * 8 + b] * dT[b];
        }
    expectTrue("1^T C dT is not trivially zero", std::abs(consistent) > 1e5);
    expectNear("lumped and consistent move the same enthalpy", lumped, consistent,
               1e-11 * std::abs(consistent));

    // `T^T K T = k |grad T|^2 V` for a linear field, on this distorted element.
    // Exact, because a trilinear hexahedron reproduces a linear field exactly, and
    // it is the closed form that ties the quadrature to the element's volume.
    const Vec3 grad{137.0, -83.0, 211.0};
    double field[8];
    for (int a = 0; a < 8; ++a)
        field[a] = 300.0 + grad.x * nodes[a * 3] + grad.y * nodes[a * 3 + 1] +
                   grad.z * nodes[a * 3 + 2];
    double quadratic = 0.0;
    for (int a = 0; a < 8; ++a)
        for (int b = 0; b < 8; ++b) quadratic += field[a] * ke[a * 8 + b] * field[b];
    const double exact = 50.0 * dot(grad, grad) * volume;
    expectNear("T^T K T = k |grad T|^2 V", quadratic, exact, 1e-10 * exact);

    // The capacity's row sums are `integral rho c N_a dV`, which is exactly what
    // `solidshell::elementMass` computes for a density. Comparing them puts this
    // file's shape functions against validated ones.
    double mass[ss::kNodes];
    ss::elementMass(nodes, 4.0e6, mass);
    double worstRow = 0.0;
    for (int a = 0; a < 8; ++a) worstRow = std::max(worstRow, std::abs(cl[a * 8 + a] - mass[a]));
    expectTrue("capacity row sums match solidshell::elementMass",
               worstRow < 1e-12 * 4.0e6 * volume);

    // An inverted element is refused rather than integrated with a negative weight.
    double inverted[ss::kDof];
    for (int i = 0; i < ss::kDof; ++i) inverted[i] = nodes[i];
    for (int a = 0; a < 4; ++a) std::swap(inverted[a * 3 + 2], inverted[(a + 4) * 3 + 2]);
    th::Forms bad;
    expectTrue("an inverted element is refused", !th::computeForms(inverted, bad) && !bad.ok);
}

// --- 3. The conduction patch test ----------------------------------------------------

void testSteadyLinearField() {
    // An oblique gradient on a randomly distorted block. A gradient down one axis
    // is reproduced by an operator with two of its three Jacobian columns wrong;
    // this one is not.
    const Vec3 gradient{137.0, -83.0, 211.0};
    const double base = 400.0;
    ss::HexMesh mesh = makeBlock({1.3, 0.9, 0.7}, 4, 3, 3, 0.28);

    th::Problem problem;
    problem.mesh = &mesh;
    problem.material = roundSteel();
    const auto exact = [&](std::size_t n) {
        return base + dot(gradient, nodeAt(mesh, n));
    };

    std::size_t held = 0, interior = 0;
    double lo = 1e30, hi = -1e30;
    for (std::size_t n = 0; n < mesh.nodeCount(); ++n) {
        const Vec3 p = nodeAt(mesh, n);
        const bool onBoundary = p.x < 1e-9 || p.x > 1.3 - 1e-9 || p.y < 1e-9 || p.y > 0.9 - 1e-9 ||
                                p.z < 1e-9 || p.z > 0.7 - 1e-9;
        lo = std::min(lo, exact(n));
        hi = std::max(hi, exact(n));
        if (onBoundary) {
            hold(problem, mesh, n, exact(n));
            ++held;
        } else {
            ++interior;
        }
    }
    expectTrue("the patch has an interior to get wrong", interior >= 12);
    expectTrue("the prescribed field spans hundreds of kelvin", hi - lo > 200.0);

    th::Solver solver;
    std::string why;
    expectTrue("the patch prepares: " + why, solver.prepare(problem, base, &why));
    expectTrue("the patch solves: " + why, solver.solveSteady(&why));

    double worst = 0.0;
    for (std::size_t n = 0; n < mesh.nodeCount(); ++n)
        worst = std::max(worst, std::abs(solver.temperature()[n] - exact(n)));
    std::printf("     patch test: %zu held, %zu free, worst %.3e K over a %.0f K span\n", held,
                interior, worst, hi - lo);
    expectTrue("a linear field is reproduced exactly", worst < 1e-9);
    expectTrue("the assembly dropped nothing", solver.account().equilibriumResidual < 1e-6);

    // The same on a mesh whose last element is **collapsed** -- the wedge a
    // degenerate plate panel extrudes to, of which the reference ferry has 166. Its
    // `det J` is exactly zero at the coincident corners and positive everywhere the
    // element is integrated, so it must still reproduce a linear field.
    ss::HexMesh wedge = makeBlock({1.0, 1.0, 1.0}, 2, 2, 2, 0.0);
    const std::size_t collapseFrom = wedge.index[7 * 8 + 1], collapseTo = wedge.index[7 * 8 + 2];
    for (int i = 0; i < 3; ++i)
        wedge.position[collapseFrom * 3 + static_cast<std::size_t>(i)] =
            wedge.position[collapseTo * 3 + static_cast<std::size_t>(i)];

    th::Problem wp;
    wp.mesh = &wedge;
    wp.material = roundSteel();
    std::size_t wedgeInterior = 0;
    for (std::size_t n = 0; n < wedge.nodeCount(); ++n) {
        const Vec3 p = nodeAt(wedge, n);
        const bool onBoundary = p.x < 1e-9 || p.x > 1 - 1e-9 || p.y < 1e-9 || p.y > 1 - 1e-9 ||
                                p.z < 1e-9 || p.z > 1 - 1e-9;
        if (onBoundary)
            hold(wp, wedge, n, base + dot(gradient, p));
        else
            ++wedgeInterior;
    }
    expectTrue("the collapsed mesh has a free node", wedgeInterior >= 1);
    th::Solver ws;
    expectTrue("a collapsed element prepares: " + why, ws.prepare(wp, base, &why));
    expectTrue("a collapsed element solves: " + why, ws.solveSteady(&why));
    double wedgeWorst = 0.0;
    for (std::size_t n = 0; n < wedge.nodeCount(); ++n)
        wedgeWorst = std::max(wedgeWorst,
                              std::abs(ws.temperature()[n] - (base + dot(gradient, nodeAt(wedge, n)))));
    std::printf("     collapsed element: worst %.3e K\n", wedgeWorst);
    expectTrue("a collapsed element reproduces a linear field", wedgeWorst < 1e-9);
}

// --- 4. Steady conduction through a plate, and its flux -------------------------------

void testSteadyPlate() {
    // Graded elements, so the twelve of them are all different sizes: an error that
    // cancels between two equal elements does not cancel here.
    const double length = 0.6, width = 0.05;
    ss::HexMesh mesh = makeBar(length, width, 12, 1.25);
    const double hot = kC + 500.0, cold = kC + 20.0;

    th::Problem problem;
    problem.mesh = &mesh;
    problem.material = roundSteel();
    for (std::size_t n = 0; n < mesh.nodeCount(); ++n) {
        const double x = mesh.position[n * 3];
        if (x < 1e-12) hold(problem, mesh, n, hot);
        if (x > length - 1e-12) hold(problem, mesh, n, cold);
    }

    th::Solver solver;
    std::string why;
    expectTrue("the plate prepares: " + why, solver.prepare(problem, cold, &why));
    const double startEnthalpy = solver.account().enthalpy;
    expectTrue("the plate solves: " + why, solver.solveSteady(&why));

    double worst = 0.0;
    for (std::size_t n = 0; n < mesh.nodeCount(); ++n) {
        const double x = mesh.position[n * 3];
        worst = std::max(worst, std::abs(solver.temperature()[n] - (hot + (cold - hot) * x / length)));
    }
    expectTrue("steady conduction is linear in x", worst < 1e-9);
    expectTrue("the profile is not flat", std::abs(hot - cold) > 400.0);
    expectEqual("a constant-property steady solve is one pass", solver.iterations(), 1);

    // The enthalpy of a linear field over a uniform-section bar is
    // `rho c V (Tmean - 20 C)` with `Tmean` the mean of the two ends, exactly --
    // the Gauss rule integrates a linear field without error.
    const double capacity =
        problem.material.density * problem.material.specificHeat * meshVolume(mesh);
    const double expected = capacity * (0.5 * (hot + cold) - (kC + 20.0));
    expectNear("the steady field's enthalpy", solver.account().enthalpy, expected,
               1e-9 * std::abs(expected));
    // The *change* is against the state `prepare` actually took, which is not the
    // uniform field it was handed: `prepare` applies the prescribed values, so the
    // hot end already stands at 500 C before a single step. Asserting against
    // `capacity * (cold - datum)` was wrong by 1.8% for exactly that reason -- the
    // test's expectation, not the solver's arithmetic.
    expectNear("the enthalpy change is the difference the account carries",
               solver.account().enthalpyChange, solver.account().enthalpy - startEnthalpy,
               1e-9 * std::abs(expected));
    expectTrue("the enthalpy change is not zero", std::abs(solver.account().enthalpyChange) > 1e6);

    // One step of a transient at a huge timestep reaches the same steady state, and
    // the reaction it reports is the analytic conduction rate `k A dT / L`. That is
    // the closed form for the *magnitude* of the boundary flux, which the patch
    // test above says nothing about.
    th::Solver transient;
    expectTrue("the transient prepares", transient.prepare(problem, cold, &why));
    for (int i = 0; i < 40; ++i) expectTrue("a long step runs: " + why, transient.step(1e5, &why));
    const double conduction = problem.material.conductivity * width * width * (hot - cold) / length;
    // The account's `prescribedHeat` is a time integral, so the *rate* is what the
    // steady answer predicts; take it from the last step alone.
    th::Solver rate;
    rate.prepare(problem, cold, &why);
    for (int i = 0; i < 39; ++i) rate.step(1e5, &why);
    const double before = rate.account().prescribedHeat;
    rate.step(1e5, &why);
    const double power = (rate.account().prescribedHeat - before) / 1e5;
    std::printf("     steady flux: net %.6e W into a body carrying %.6e W through\n", power,
                conduction);
    expectTrue("at steady state the net heat in is nothing", std::abs(power) < 1e-6 * conduction);

    // The hot end alone, which is the conduction rate itself. Summed from the
    // residual over the nodes at x = 0.
    // And the heat the hot end has to supply, computed from the solved gradient
    // rather than from the boundary values: `-k grad T . A`.
    const double solvedGradient =
        (solver.temperature()[4] - solver.temperature()[0]) / mesh.position[4 * 3];
    const double hotSide = -problem.material.conductivity * solvedGradient * width * width;
    expectNear("the solved gradient carries the analytic conduction rate", hotSide, conduction,
               1e-9 * conduction);
}

// --- 5. The semi-infinite solid ---------------------------------------------------------

// One run of the step-surface problem; returns the field and the mesh's x values.
struct StepRun {
    std::vector<double> x, temperature;
    double time = 0;
};

StepRun runStep(int elements, double timestep, double duration, bool lumped = false) {
    const double length = 0.30, width = 0.01;
    ss::HexMesh mesh = makeBar(length, width, elements);
    const double t0 = kC + 20.0, ts = kC + 900.0;

    th::Problem problem;
    problem.mesh = &mesh;
    problem.material = roundSteel();
    problem.lumpedCapacity = lumped;
    for (std::size_t n = 0; n < mesh.nodeCount(); ++n)
        if (mesh.position[n * 3] < 1e-12) hold(problem, mesh, n, ts);

    th::Solver solver;
    std::string why;
    if (!solver.prepare(problem, t0, &why)) return {};
    const int steps = static_cast<int>(std::llround(duration / timestep));
    for (int i = 0; i < steps; ++i)
        if (!solver.step(timestep, &why)) return {};

    StepRun run;
    run.time = solver.time();
    // The field on the x axis: nodes are laid out with the two transverse indices
    // fastest, so node 4i is the (i, 0, 0) corner.
    for (int i = 0; i <= elements; ++i) {
        run.x.push_back(mesh.position[static_cast<std::size_t>(i) * 4 * 3]);
        run.temperature.push_back(solver.temperature()[static_cast<std::size_t>(i) * 4]);
    }
    return run;
}

double erfSolution(double x, double t) {
    const double t0 = kC + 20.0, ts = kC + 900.0;
    return ts + (t0 - ts) * std::erf(x / (2.0 * std::sqrt(kRoundAlpha * t)));
}

void testSemiInfinite() {
    const double duration = 60.0;
    const StepRun run = runStep(240, 0.05, duration);
    expectTrue("the step run completed", run.temperature.size() == 241);
    if (run.temperature.size() != 241) return;

    // **Vacuity guards.** The domain must be semi-infinite over the run -- the far
    // end still at its initial value -- and the front must actually be inside it,
    // or the comparison is between two constants.
    const double t0 = kC + 20.0, ts = kC + 900.0;
    expectTrue("the far end has not felt the step",
               std::abs(run.temperature.back() - t0) < 0.05);
    const double mid = run.temperature[run.temperature.size() / 8];
    const double fraction = (mid - t0) / (ts - t0);
    expectTrue("the front is inside the domain", fraction > 0.15 && fraction < 0.85);

    double worst = 0.0, worstAt = 0.0;
    for (std::size_t i = 1; i < run.x.size(); ++i) {
        const double e = std::abs(run.temperature[i] - erfSolution(run.x[i], run.time));
        if (e > worst) {
            worst = e;
            worstAt = run.x[i];
        }
    }
    std::printf("     erf: 240 elements, dt = 0.05 s, t = %.0f s -- worst %.4f K at x = %.3f m"
                " (%.3f%% of the step)\n",
                run.time, worst, worstAt, 100.0 * worst / (ts - t0));
    expectTrue("the transient matches the semi-infinite closed form", worst < 1.5);

    // The same problem at **eighty times the explicit stability limit**, which is
    // the case for solving this implicitly and is stated as the measurement it is.
    // It has to remain stable, and it has to stay within a fraction of a kelvin on
    // an 880 K step.
    ss::HexMesh probe = makeBar(0.30, 0.01, 240);
    const double limit = th::explicitLimit(probe, roundSteel());
    const StepRun coarse = runStep(240, 80.0 * limit, duration);
    expectTrue("the coarse run completed", coarse.temperature.size() == 241);
    if (coarse.temperature.size() != 241) return;
    double coarseWorst = 0.0;
    for (std::size_t i = 1; i < coarse.x.size(); ++i)
        coarseWorst = std::max(coarseWorst,
                               std::abs(coarse.temperature[i] - erfSolution(coarse.x[i], coarse.time)));
    std::printf("     explicit limit %.4f s for a 1.25 mm element; backward Euler at 80x it"
                " (%.1f s steps, %.0f of them) -- worst %.2f K, %.2f%% of the step\n",
                limit, 80.0 * limit, duration / (80.0 * limit), coarseWorst,
                100.0 * coarseWorst / (ts - t0));
    // Stable at any step -- that is what implicit buys -- and first order, so
    // twelve steps over a *discontinuous* boundary condition costs about one
    // percent of the step it applied. Asserted at the measured value rather than
    // at a round one, and bracketed below as well: an error far *smaller* than
    // first order predicts would mean the coarse run had not actually taken
    // coarse steps.
    expectTrue("backward Euler at 80x the explicit limit is stable and first-order accurate",
               coarseWorst > 3.0 && coarseWorst < 12.0);
}

// --- 5b. The explicit limit itself, as a closed form -----------------------------------

void testExplicitLimit() {
    // For a rectangular hexahedron the highest mode is the jump through the
    // thinnest direction, so `2/lambda_max` is exactly `h_min^2 / (2 alpha)`. That
    // is an identity and is asserted as one -- it is what caught a power iteration
    // whose checkerboard start vector was orthogonal, by symmetry, to precisely the
    // mode it was supposed to find, and which therefore reported a stable step
    // three times too long.
    const StructuralMaterial steel = roundSteel();
    const double alpha = steel.conductivity / (steel.density * steel.specificHeat);
    expectNear("the round material's diffusivity", alpha, kRoundAlpha, 1e-20);

    struct Case { double lx, ly, t; int nz; };
    const Case cases[] = {{0.700, 2.400, 0.012, 1},  // a ferry shell panel as meshed
                          {0.050, 0.050, 0.012, 1},  // a crush-zone element
                          {0.050, 0.050, 0.012, 4},  // the same, resolved through
                          {0.0015, 0.05, 0.05, 1}};  // a thin surface layer
    double worst = 0.0;
    for (const Case& c : cases) {
        ss::HexMesh mesh = ss::makePlateMesh(c.lx, c.ly, c.t, 1, 1, c.nz);
        const double got = th::explicitLimit(mesh, steel);
        const double h = std::min({c.lx, c.ly, c.t / c.nz});
        worst = std::max(worst, std::abs(got - h * h / (2.0 * alpha)) / got);
    }
    expectTrue("the explicit limit is h_min^2 / 2 alpha exactly", worst < 1e-12);

    // The two scaling laws, which no transcribed constant can fake.
    ss::HexMesh small = ss::makePlateMesh(0.05, 0.05, 0.012, 1, 1, 1);
    ss::HexMesh large = ss::makePlateMesh(0.10, 0.10, 0.024, 1, 1, 1);
    StructuralMaterial hotter = steel;
    hotter.conductivity *= 4.0;
    expectNear("the limit scales as h^2",
               th::explicitLimit(large, steel) / th::explicitLimit(small, steel), 4.0, 1e-12);
    expectNear("the limit scales as 1/alpha",
               th::explicitLimit(small, hotter) / th::explicitLimit(small, steel), 0.25, 1e-12);

    // **The correction this file exists to record.** The usual claim is that
    // explicit conduction on ship plating is limited to milliseconds. For AH36 at
    // 12 mm it is 4.66 *seconds*, a thousand times more, because steel's
    // diffusivity is 1.5e-5 m^2/s and 12 mm is not a small length.
    ss::HexMesh plating = ss::makePlateMesh(0.700, 2.400, 0.012, 1, 1, 1);
    const double ferry = th::explicitLimit(plating, ah36Steel());
    ss::HexMesh resolved = ss::makePlateMesh(0.700, 2.400, 0.012, 1, 1, 4);
    const double refined = th::explicitLimit(resolved, ah36Steel());
    std::printf("     explicit limit: 12 mm AH36 plating %.3f s (one element through),"
                " %.3f s at four -- SECONDS, not milliseconds\n",
                ferry, refined);
    expectTrue("12 mm plating is seconds, not milliseconds", ferry > 4.0 && ferry < 5.0);
    expectNear("and four elements through it costs a factor of sixteen", ferry / refined, 16.0,
               1e-9);
}

// --- 5c. What the lumped capacity is for -------------------------------------------------

double undershootOf(bool lumped, double timestep, int steps) {
    // A step change in surface temperature, and how far *below* its initial value
    // the interior dips. A consistent capacity gives a mass matrix with negative
    // off-diagonal influence on the step, so backward Euler produces a physically
    // impossible dip ahead of the front unless the step is long enough to smear it.
    const double length = 0.30, width = 0.01;
    const int n = 60;  // h = 5 mm
    ss::HexMesh mesh = makeBar(length, width, n);
    const double t0 = kC + 20.0, ts = kC + 900.0;
    th::Problem problem;
    problem.mesh = &mesh;
    problem.material = roundSteel();
    problem.lumpedCapacity = lumped;
    for (std::size_t node = 0; node < mesh.nodeCount(); ++node)
        if (mesh.position[node * 3] < 1e-12) hold(problem, mesh, node, ts);
    th::Solver solver;
    std::string why;
    if (!solver.prepare(problem, t0, &why)) return 0.0;
    double dip = 0.0;
    for (int i = 0; i < steps; ++i) {
        if (!solver.step(timestep, &why)) return 0.0;
        for (double v : solver.temperature()) dip = std::max(dip, t0 - v);
    }
    return dip;
}

void testUndershoot() {
    // The criterion is `dt >= rho c h^2 / (6 k)`, which for these 5 mm elements is
    // `h^2 / (6 alpha)` = 0.333 s. Measured on both sides of it, and with the
    // lumped capacity as the control -- lumped backward Euler is unconditionally
    // monotone and must not dip at all.
    const double h = 0.30 / 60.0;
    const double criterion = h * h / (6.0 * kRoundAlpha);
    const double shortStep = 0.1 * criterion, longStep = 3.0 * criterion;

    const double dipConsistentShort = undershootOf(false, shortStep, 40);
    const double dipLumpedShort = undershootOf(true, shortStep, 40);
    const double dipConsistentLong = undershootOf(false, longStep, 8);
    std::printf("     undershoot below the initial 20 C, h = 5 mm, criterion dt = %.3f s:\n"
                "       consistent at %.4f s (0.1x): %.3f K;  lumped, same step: %.2e K;"
                "  consistent at %.3f s (3x): %.2e K\n",
                criterion, shortStep, dipConsistentShort, dipLumpedShort, longStep,
                dipConsistentLong);

    expectTrue("a consistent capacity undershoots below the criterion",
               dipConsistentShort > 1.0);
    expectTrue("a lumped capacity never undershoots", dipLumpedShort < 1e-9);
    expectTrue("a consistent capacity above the criterion does not undershoot",
               dipConsistentLong < 1e-9);
    // **Vacuity guard.** The two steps have to bracket the criterion, or this is
    // two runs of the same thing.
    expectTrue("the two steps bracket the criterion",
               shortStep < criterion && longStep > criterion);
}

// --- 6. Convergence orders --------------------------------------------------------------

void testOrders() {
    // **Time.** By Richardson on a fixed mesh, so the spatial error cancels
    // identically and what is measured is the time integration alone. Backward
    // Euler is first order, so successive halvings differ by a factor of two.
    const double duration = 20.0;
    std::vector<double> difference;
    std::vector<double> previous;
    for (int level = 0; level < 5; ++level) {
        const double dt = 0.8 / (1 << level);
        const StepRun run = runStep(160, dt, duration);
        if (run.temperature.empty()) {
            expectTrue("the time refinement ran", false);
            return;
        }
        if (!previous.empty()) {
            double norm = 0.0;
            for (std::size_t i = 0; i < run.temperature.size(); ++i) {
                const double d = run.temperature[i] - previous[i];
                norm += d * d;
            }
            difference.push_back(std::sqrt(norm / static_cast<double>(run.temperature.size())));
        }
        previous = run.temperature;
    }
    std::printf("     time refinement (dt 0.8 s halved): ");
    for (double d : difference) std::printf("%.4e ", d);
    const double timeOrder = orderOf(difference, 2.0);
    std::printf("-> order %.3f\n", timeOrder);
    expectTrue("the time refinement is not measuring rounding", difference.front() > 1e-3);
    expectTrue("backward Euler is first order in time", timeOrder > 0.90 && timeOrder < 1.10);

    // **Space.** Against the erf itself, with the step small enough that the time
    // error is well below the spatial one at the finest mesh -- checked, not
    // assumed: the time error at dt = 0.01 s is bounded by the Richardson
    // difference above, which is far below the coarsest spatial error here.
    std::vector<double> spatial;
    for (int level = 0; level < 4; ++level) {
        const int n = 30 << level;
        const StepRun run = runStep(n, 0.01, duration);
        if (run.temperature.empty()) {
            expectTrue("the space refinement ran", false);
            return;
        }
        // L2 over the bar, by the trapezium rule on a uniform grid.
        double norm = 0.0;
        for (std::size_t i = 0; i < run.x.size(); ++i) {
            const double e = run.temperature[i] - erfSolution(run.x[i], run.time);
            const double w = (i == 0 || i + 1 == run.x.size()) ? 0.5 : 1.0;
            norm += w * e * e;
        }
        spatial.push_back(std::sqrt(norm / static_cast<double>(run.x.size())));
    }
    std::printf("     space refinement (30 elements doubled): ");
    for (double d : spatial) std::printf("%.4e ", d);
    const double spaceOrder = orderOf(spatial, 2.0);
    std::printf("-> order %.3f\n", spaceOrder);
    expectTrue("the space refinement is not measuring rounding", spatial.front() > 1e-2);
    expectTrue("the element is second order in space", spaceOrder > 1.80 && spaceOrder < 2.20);
}

// --- 7. The energy account -----------------------------------------------------------------

void testEnergyAccount(bool lumped) {
    const std::string tag = lumped ? " (lumped)" : " (consistent)";
    // Four channels at once, all asymmetric: a flux on one face, convection on
    // another, a held temperature on a third, and a volumetric source in a wedge of
    // the elements. Nothing about this fixture is symmetric, so nothing cancels.
    ss::HexMesh mesh = makeBlock({0.8, 0.5, 0.012}, 8, 5, 1, 0.15);
    th::Problem problem;
    problem.mesh = &mesh;
    problem.material = roundSteel();
    problem.lumpedCapacity = lumped;

    for (std::size_t n = 0; n < mesh.nodeCount(); ++n)
        if (mesh.position[n * 3] < 1e-12) hold(problem, mesh, n, kC + 20.0);

    // **Element zero carries one.** A source on every third element starting at the
    // second leaves element 0 empty, and a power sum that skipped element 0 then
    // reported exactly the right total -- which mutation testing found and the
    // account could not.
    problem.volumetricSource.assign(mesh.elementCount(), 0.0);
    for (std::size_t e = 0; e < mesh.elementCount(); ++e)
        if (e % 3 == 0) problem.volumetricSource[e] = 4.0e6;

    const std::vector<th::BoundaryFace> boundary = th::boundaryFaces(mesh);
    th::Film fire, cool;
    fire.flux = 22000.0;  // W/m^2, roughly a post-flashover compartment on steel
    cool.coefficient = 35.0;
    cool.ambient = kC + 20.0;
    double fireArea = 0.0, coolArea = 0.0, totalArea = 0.0;
    for (const th::BoundaryFace& f : boundary) {
        totalArea += f.area;
        if (f.normal.z > 0.5) {
            fire.face.push_back(f);
            fireArea += f.area;
        } else if (f.normal.z < -0.5) {
            cool.face.push_back(f);
            coolArea += f.area;
        }
    }
    problem.film = {fire, cool};

    // The boundary the mesher found must be the boundary of the block, or the two
    // films are applied to some other surface.
    const double exactArea = 2 * (0.8 * 0.5 + 0.8 * 0.012 + 0.5 * 0.012);
    expectNear("the boundary faces close the block" + tag, totalArea, exactArea, 1e-9 * exactArea);
    expectNear("the fire covers one whole face" + tag, fireArea, 0.4, 1e-9);
    expectNear("the coolant covers the other" + tag, coolArea, 0.4, 1e-9);

    th::Solver solver;
    std::string why;
    expectTrue("the account problem prepares: " + why, solver.prepare(problem, kC + 20.0, &why));

    // Varying steps, so a scheme that only balanced at a fixed step is caught.
    // A constant-property step is one solve and no iteration at all.
    double worstEquilibrium = 0.0;
    for (int i = 0; i < 60; ++i) {
        const double dt = (i % 4 == 0) ? 2.0 : ((i % 4 == 1) ? 0.5 : 5.0);
        expectTrue("an accounted step runs: " + why, solver.step(dt, &why));
        expectEqual("a linear step is one solve" + tag, solver.iterations(), 1);
        worstEquilibrium = std::max(worstEquilibrium, solver.account().equilibriumResidual);
    }
    const th::Account& a = solver.account();
    const double largest =
        std::max({std::abs(a.enthalpyChange), std::abs(a.prescribedHeat), std::abs(a.filmHeat),
                  std::abs(a.sourceHeat)});
    std::printf("     account%s: dH %.6e = held %.4e + film %.4e + source %.4e, residual %.3e J\n",
                tag.c_str(), a.enthalpyChange, a.prescribedHeat, a.filmHeat, a.sourceHeat,
                a.residual());

    // **Vacuity guards.** Every channel has to be carrying a real share, or the
    // balance is between three zeros.
    expectTrue("the block actually heated up" + tag, a.enthalpyChange > 1e5);
    expectTrue("the film carries a real share" + tag, std::abs(a.filmHeat) > 0.05 * largest);
    expectTrue("the source carries a real share" + tag, std::abs(a.sourceHeat) > 0.05 * largest);
    expectTrue("the held boundary carries a real share" + tag,
               std::abs(a.prescribedHeat) > 0.01 * largest);

    expectTrue("the energy account closes to machine precision" + tag,
               std::abs(a.residual()) < 1e-12 * largest);
    // A free row of `A T1 - b` that is not zero means the band silently dropped a
    // term -- the exact shape of the defect `CLAUDE.md` records for the node
    // ordering. Scaled off the heat actually flowing, not a constant.
    const double heatRate = std::abs(a.filmHeat + a.sourceHeat) / solver.time();
    expectTrue("the band held every term" + tag, worstEquilibrium < 1e-9 * heatRate);

    // Only the *last* factorisation is kept -- caching one per distinct step would
    // be an unbounded amount of memory decided by a caller's step schedule. So the
    // 2, 0.5, 5, 5 cycle refactors three times in every four steps, and the 5-to-5
    // repeat is the one that does not. `testCost` is the case that matters: a
    // fixed step over 201 steps factors once.
    std::printf("     %d factorisations over 60 steps, half-bandwidth %zu of %zu free nodes\n",
                solver.factorisations(), solver.halfBandwidth(), solver.freeNodes());
    expectEqual("only a changed step refactors" + tag, solver.factorisations(), 45);
}

// --- 8. Newton cooling ----------------------------------------------------------------------

void testNewtonCooling() {
    // A thin plate cooling by convection on its two large faces. At small Biot
    // number the plate is isothermal and
    // `T = Tinf + (T0 - Tinf) exp(-(h1 A1 + h2 A2) t / (rho c V))`.
    //
    // The edges are deliberately left adiabatic. With convection on all six faces
    // the closed form does not apply: heat leaving an edge has to travel 0.2 m in
    // plane, a Biot number of 0.05 rather than 0.001, and the plate develops a
    // real 15 K spread that is physics and not error. Choosing the fixture the
    // closed form is the answer to is the point; widening the tolerance until the
    // wrong fixture passed would have been the mistake.
    //
    // The two faces get *different* coefficients, so the fixture is asymmetric and
    // a solver that used one film for both would be caught.
    ss::HexMesh mesh = makeBlock({0.4, 0.3, 0.010}, 4, 3, 1, 0.0);
    th::Problem problem;
    problem.mesh = &mesh;
    problem.material = roundSteel();

    const double hTop = 12.0, hBottom = 7.0, ambient = kC + 20.0, start = kC + 620.0;
    th::Film top, bottom;
    top.coefficient = hTop;
    top.ambient = ambient;
    bottom.coefficient = hBottom;
    bottom.ambient = ambient;
    double topArea = 0.0, bottomArea = 0.0;
    for (const th::BoundaryFace& f : th::boundaryFaces(mesh)) {
        if (f.normal.z > 0.5) {
            top.face.push_back(f);
            topArea += f.area;
        } else if (f.normal.z < -0.5) {
            bottom.face.push_back(f);
            bottomArea += f.area;
        }
    }
    problem.film = {top, bottom};
    expectNear("the cooled faces are the plate's own", topArea, 0.12, 1e-12);
    expectNear("both of them", bottomArea, 0.12, 1e-12);

    const double volume = meshVolume(mesh);
    expectNear("the plate's volume", volume, 0.4 * 0.3 * 0.010, 1e-15);
    // Biot on the half-thickness, which is the conduction length for a plate
    // cooled from both sides. It has to be small or the lumped closed form is not
    // the answer to this problem at all.
    const double biot = hTop * 0.005 / problem.material.conductivity;
    expectTrue("the Biot number is small enough for a lumped answer", biot < 0.01);

    const double capacity =
        problem.material.density * problem.material.specificHeat * volume;
    const double lambda = (hTop * topArea + hBottom * bottomArea) / capacity;
    th::Solver solver;
    std::string why;
    expectTrue("the cooling plate prepares: " + why, solver.prepare(problem, start, &why));
    const double dt = 4.0;
    const int steps = 1580;  // just over three time constants
    for (int i = 0; i < steps; ++i) expectTrue("a cooling step runs: " + why, solver.step(dt, &why));

    const double t = solver.time();
    const double analytic = ambient + (start - ambient) * std::exp(-lambda * t);
    double lo = 1e30, hi = -1e30;
    for (double v : solver.temperature()) {
        lo = std::min(lo, v);
        hi = std::max(hi, v);
    }
    const double mean = 0.5 * (lo + hi);
    std::printf("     Newton cooling: Bi = %.4f, tau = %.0f s, after %.0f s (%.2f tau) got"
                " %.4f K, exact %.4f K, spread %.3e K\n",
                biot, 1.0 / lambda, t, lambda * t, mean, analytic, hi - lo);

    // **Vacuity guards.** A plate that has barely cooled agrees with an exponential
    // it has not yet distinguished from a constant.
    expectTrue("the plate cooled through nine tenths of its excess",
               (start - hi) > 0.9 * (start - ambient));
    // Through-thickness only: with unequal coefficients the profile is linear
    // across the plate, so the extremes bracket the mean and the spread is the
    // Biot number's worth of it.
    expectTrue("the plate is isothermal to its Biot number",
               hi - lo < 0.02 * (mean - ambient));

    // Backward Euler's own decay rate is `ln(1 + lambda dt)/dt`, so the exact
    // discrete answer is the continuum exponential with that rate. The allowance is
    // the first-order error that follows -- `lambda^2 dt t / 2` in the exponent --
    // computed rather than chosen, and then checked not to be a licence.
    const double allowance = (mean - ambient) * (lambda * lambda * dt * t * 0.5) * 1.3;
    expectNear("Newton cooling matches the lumped closed form", mean, analytic,
               allowance + 1e-9);
    expectTrue("the allowance is not a licence", allowance < 0.01 * (analytic - ambient));

    // And the account still closes: convection is the only channel here.
    const th::Account& a = solver.account();
    expectTrue("cooling removes energy", a.enthalpyChange < -1e6);
    expectTrue("the cooling account closes",
               std::abs(a.residual()) < 1e-12 * std::abs(a.enthalpyChange));
}

// --- 9. Temperature-dependent properties: the Kirchhoff transform ----------------------------

// The exact steady profile with `k = k0 - k1 theta`. The Kirchhoff potential
// `psi(T) = integral k dT` is linear in x whatever `k(T)` is, so inverting a
// quadratic gives the temperature.
double kirchhoff(double celsiusTemperature) {
    const double t = celsiusTemperature;
    return 54.0 * t - 0.5 * 3.33e-2 * t * t;
}
double kirchhoffInverse(double psi) {
    // 54 t - 0.01665 t^2 = psi, the root below the turning point at t = 1621 C.
    const double a = 0.5 * 3.33e-2;
    return (54.0 - std::sqrt(54.0 * 54.0 - 4.0 * a * psi)) / (2.0 * a);
}

void testKirchhoffSteady() {
    // Both ends below 800 C, where EN 1993-1-2's conductivity is exactly linear and
    // the closed form is exact.
    //
    // **The discrete answer is exact at the nodes, and that is a stronger result
    // than the second-order convergence this test was first written to assert.**
    // The reason is worth recording because it is a property of the quadrature and
    // not luck: the 1D element conductance is `(1/h) integral k(T(xi)) dxi`, so the
    // discrete flux through an element is the two-point Gauss rule applied to
    // `integral k dT` with `T` linear in `xi`. For `k` linear in `T` that
    // integrand is linear in `xi`, and a two-point Gauss rule is exact on a cubic.
    // So every element carries exactly the Kirchhoff potential difference the
    // continuum does, the nodal balance is the exact one, and there is no
    // discretisation error left to converge. Asserting an order here would have
    // been asserting a defect.
    const double length = 0.4, width = 0.02;
    const double hotC = 700.0, coldC = 100.0;

    std::vector<double> error;
    for (int level = 0; level < 4; ++level) {
        const int n = 4 << level;
        ss::HexMesh mesh = makeBar(length, width, n);
        th::Problem problem;
        problem.mesh = &mesh;
        problem.material = ah36Steel();
        problem.temperatureDependent = true;
        for (std::size_t node = 0; node < mesh.nodeCount(); ++node) {
            const double x = mesh.position[node * 3];
            if (x < 1e-12) hold(problem, mesh, node, kC + hotC);
            if (x > length - 1e-12) hold(problem, mesh, node, kC + coldC);
        }
        th::Solver solver;
        std::string why;
        expectTrue("the Kirchhoff bar prepares: " + why, solver.prepare(problem, kC + coldC, &why));
        expectTrue("the Kirchhoff bar solves: " + why, solver.solveSteady(&why));

        double worst = 0.0, worstLinear = 0.0;
        for (int i = 0; i <= n; ++i) {
            const double x = mesh.position[static_cast<std::size_t>(i) * 4 * 3];
            const double s = x / length;
            const double psi = kirchhoff(hotC) + s * (kirchhoff(coldC) - kirchhoff(hotC));
            const double exact = kC + kirchhoffInverse(psi);
            const double linear = kC + hotC + s * (coldC - hotC);
            worst = std::max(worst,
                             std::abs(solver.temperature()[static_cast<std::size_t>(i) * 4] - exact));
            worstLinear = std::max(worstLinear, std::abs(exact - linear));
        }
        if (level == 0) {
            // **Vacuity guard.** If the nonlinear profile were the linear one, this
            // test would be `testSteadyPlate` again and would pass on a solver that
            // ignored `k(T)` entirely.
            std::printf("     Kirchhoff: the nonlinear profile departs from the linear one by"
                        " %.2f K, %d Picard iterations\n",
                        worstLinear, solver.iterations());
            expectTrue("k(T) actually bends the profile", worstLinear > 8.0);
            expectTrue("the nonlinear solve iterated", solver.iterations() > 1);
        }
        error.push_back(worst);
    }
    std::printf("     Kirchhoff (4 elements doubled, worst nodal error): ");
    for (double e : error) std::printf("%.3e ", e);
    std::printf("K\n");
    // At the Picard tolerance, which is what limits it -- not at the mesh size,
    // which does not enter.
    for (double e : error)
        expectTrue("the nonlinear steady profile is exact at the nodes", e < 1e-6);
    (void)0;
    expectTrue("and does not degrade with refinement", error.back() < 10.0 * error.front() + 1e-9);

    // A nonlinear *steady* solve is allowed to give up too, and must say so rather
    // than return whatever the last iterate happened to be as if it had converged.
    ss::HexMesh mesh = makeBar(length, width, 8);
    th::Problem problem;
    problem.mesh = &mesh;
    problem.material = ah36Steel();
    problem.temperatureDependent = true;
    for (std::size_t node = 0; node < mesh.nodeCount(); ++node) {
        const double x = mesh.position[node * 3];
        if (x < 1e-12) hold(problem, mesh, node, kC + hotC);
        if (x > length - 1e-12) hold(problem, mesh, node, kC + coldC);
    }
    th::Solver stubborn;
    std::string why;
    stubborn.setPicard(1e-30, 2);
    expectTrue("the stubborn steady solver prepares", stubborn.prepare(problem, kC + coldC, &why));
    why.clear();
    expectTrue("an unconverged steady solve is reported",
               !stubborn.solveSteady(&why) && !why.empty());
    expectEqual("it used its whole budget", stubborn.iterations(), 2);
}

// --- 10. Temperature-dependent transient across the phase change --------------------------------

void testSpikeAccount() {
    // Heat 12 mm plating with a post-flashover flux until it is through the 735 C
    // transition, with the EN 1993-1-2 curves live. The point is the energy account:
    // a tangent capacity leaves a residual of order `dt` per step here, and the
    // secant leaves none.
    ss::HexMesh mesh = makeBlock({0.6, 0.4, 0.012}, 6, 4, 1, 0.12);
    th::Problem problem;
    problem.mesh = &mesh;
    problem.material = ah36Steel();
    problem.temperatureDependent = true;

    th::Film fire;
    fire.flux = 45000.0;  // W/m^2
    for (const th::BoundaryFace& f : th::boundaryFaces(mesh))
        if (f.normal.z > 0.5) fire.face.push_back(f);
    problem.film = {fire};

    th::Solver solver;
    std::string why;
    expectTrue("the fire problem prepares: " + why, solver.prepare(problem, kC + 20.0, &why));

    double crossedFrom = 1e30, crossedTo = -1e30;
    int totalIterations = 0;
    bool sawSpike = false;
    double spikeSeconds = 0.0;
    const int fireSteps = 600;
    double worstEquilibrium = 0.0;
    for (int i = 0; i < fireSteps; ++i) {
        expectTrue("a fire step runs: " + why, solver.step(2.0, &why));
        totalIterations += solver.iterations();
        worstEquilibrium = std::max(worstEquilibrium, solver.account().equilibriumResidual);
        double lo = 1e30, hi = -1e30;
        for (double v : solver.temperature()) {
            lo = std::min(lo, v);
            hi = std::max(hi, v);
        }
        crossedFrom = std::min(crossedFrom, lo);
        crossedTo = std::max(crossedTo, hi);
        if (lo > kC + 700.0 && hi < kC + 770.0) {
            sawSpike = true;
            spikeSeconds += 2.0;
        }
    }
    const th::Account& a = solver.account();
    std::printf("     fire: %.0f s of 45 kW/m^2 took the plate %.0f C -> %.0f C in %d Picard"
                " iterations (%.1f per step); %.0f s inside the 735 C transition\n",
                solver.time(), crossedFrom - kC, crossedTo - kC, totalIterations,
                totalIterations / 400.0, spikeSeconds);

    // **Vacuity guards.** The run has to have crossed the transition, or the secant
    // is never tested on the thing it exists for.
    expectTrue("the plate crossed the 735 C transition",
               crossedFrom < kC + 700.0 && crossedTo > kC + 800.0);
    expectTrue("the transition took real time", sawSpike && spikeSeconds > 20.0);

    const double largest = std::max(std::abs(a.enthalpyChange), std::abs(a.filmHeat));
    std::printf("     fire account: dH %.6e J, film %.6e J, residual %.3e J (%.2e relative)\n",
                a.enthalpyChange, a.filmHeat, a.residual(), a.residual() / largest);

    // **With temperature-dependent properties the account closes to the Picard
    // tolerance, not to machine precision**, and the difference is worth being
    // exact about. The secant capacity makes `1^T C dT` identically the enthalpy
    // change *of the state the system was solved for*; a Picard iteration stopped
    // at 1e-8 K has not quite reached that state, and the leftover free-node
    // residual is the whole of what is missing.
    //
    // So it is asserted as a *scaling* rather than as a constant: tightening the
    // tolerance a thousandfold has to tighten the residual with it. A residual
    // coming from a wrong quadrature would not move at all.
    expectTrue("a nonlinear account closes to the Picard tolerance",
               std::abs(a.residual()) < 1e-10 * largest);

    const auto residualAt = [&](double tolerance) {
        th::Solver s;
        std::string w;
        s.setPicard(tolerance, 200);
        s.prepare(problem, kC + 20.0, &w);
        for (int i = 0; i < 120; ++i) s.step(2.0, &w);
        const th::Account& acc = s.account();
        return std::abs(acc.residual()) / std::abs(acc.filmHeat);
    };
    const double loose = residualAt(1e-6), tight = residualAt(1e-11);
    std::printf("     nonlinear account residual: %.3e relative at a 1e-6 K Picard tolerance,"
                " %.3e at 1e-11\n",
                loose, tight);
    expectTrue("the loose tolerance leaves a residual to shrink", loose > 1e-13);
    expectTrue("the nonlinear residual is the Picard tolerance and nothing else",
               tight < 0.02 * loose);

    // A nonlinear step re-factors once per Picard iteration, which is the cost of
    // the temperature dependence and is worth being explicit about.
    expectEqual("a nonlinear step factors once per iteration", solver.factorisations(),
                totalIterations);

    // The properties are re-formed at the *accepted* state after the iteration
    // stops, so the system the account is taken from is the one the answer solves.
    // Without that the free rows of `A T1 - b` are left at the Picard tolerance
    // rather than at rounding -- which no channel of the account above can see when
    // nothing is prescribed, because film and enthalpy are both computed
    // independently of `C`.
    std::printf("     fire equilibrium residual: %.3e W against %.3e W of heat flowing\n",
                worstEquilibrium, std::abs(a.filmHeat) / solver.time());
    expectTrue("the nonlinear system is in equilibrium at the accepted state",
               worstEquilibrium < 1e-6 * std::abs(a.filmHeat) / solver.time());

    // The same, with a **held** boundary, which is the case where the reaction term
    // does depend on `C` and the account can therefore see the difference.
    th::Problem held = problem;
    held.prescribed.assign(mesh.nodeCount(), 0u);
    held.prescribedValue.assign(mesh.nodeCount(), kC + 20.0);
    for (std::size_t n = 0; n < mesh.nodeCount(); ++n)
        if (mesh.position[n * 3] < 1e-12) held.prescribed[n] = 1u;
    th::Solver heldSolver;
    expectTrue("the held fire prepares: " + why, heldSolver.prepare(held, kC + 20.0, &why));
    for (int i = 0; i < 200; ++i) expectTrue("a held fire step runs: " + why,
                                             heldSolver.step(2.0, &why));
    const th::Account& hb = heldSolver.account();
    const double biggest = std::max({std::abs(hb.enthalpyChange), std::abs(hb.filmHeat),
                                     std::abs(hb.prescribedHeat)});
    std::printf("     held fire account: dH %.6e = held %.4e + film %.4e, residual %.3e J\n",
                hb.enthalpyChange, hb.prescribedHeat, hb.filmHeat, hb.residual());
    expectTrue("the held boundary carries a real share", std::abs(hb.prescribedHeat) > 0.01 * biggest);
    expectTrue("a nonlinear account with a reaction still closes",
               std::abs(hb.residual()) < 1e-10 * biggest);

    // **Picard is allowed to give up, and has to say so.** Two iterations at an
    // impossible tolerance cannot converge; the step must fail with a reason and
    // the state must still have advanced, which is the documented choice.
    th::Solver stubborn;
    stubborn.setPicard(1e-30, 2);
    expectTrue("the stubborn solver prepares", stubborn.prepare(problem, kC + 20.0, &why));
    const double was = stubborn.temperature()[0];
    why.clear();
    expectTrue("an unconverged step is reported", !stubborn.step(50.0, &why) && !why.empty());
    expectTrue("and the state advanced anyway", stubborn.temperature()[0] > was + 1.0);
    expectEqual("it used its whole budget", stubborn.iterations(), 2);
}

// --- 11. Invariance, and what is refused ----------------------------------------------------

void testRotationInvariance() {
    // Rotating the mesh and the boundary data together must leave the temperature
    // field untouched: it is the one check that the Cartesian gradient transform is
    // a transform and not a relabelling of axes.
    ss::HexMesh mesh = makeBlock({0.9, 0.6, 0.4}, 3, 3, 2, 0.25);
    ss::HexMesh turned = mesh;
    const double ca = std::cos(0.7), sa = std::sin(0.7);
    const double cb = std::cos(-1.1), sb = std::sin(-1.1);
    for (std::size_t n = 0; n < mesh.nodeCount(); ++n) {
        const Vec3 p = nodeAt(mesh, n);
        const Vec3 q{p.x * ca - p.y * sa, p.x * sa + p.y * ca, p.z};
        const Vec3 r{q.x, q.y * cb - q.z * sb, q.y * sb + q.z * cb};
        turned.position[n * 3] = r.x;
        turned.position[n * 3 + 1] = r.y;
        turned.position[n * 3 + 2] = r.z;
    }

    const auto build = [&](const ss::HexMesh& m) {
        th::Problem p;
        p.mesh = &m;
        p.material = roundSteel();
        for (std::size_t n = 0; n < mesh.nodeCount(); ++n) {
            const Vec3 q = nodeAt(mesh, n);  // classified in the *original* frame
            if (q.x < 1e-9) hold(p, m, n, kC + 700.0);
            if (q.x > 0.9 - 1e-9) hold(p, m, n, kC + 20.0);
        }
        return p;
    };
    th::Problem straight = build(mesh), rotated = build(turned);

    th::Solver a, b;
    std::string why;
    expectTrue("the straight mesh solves", a.prepare(straight, kC + 20.0, &why) &&
                                               a.solveSteady(&why));
    expectTrue("the rotated mesh solves", b.prepare(rotated, kC + 20.0, &why) &&
                                              b.solveSteady(&why));
    double worst = 0.0, spread = 0.0;
    for (std::size_t n = 0; n < mesh.nodeCount(); ++n) {
        worst = std::max(worst, std::abs(a.temperature()[n] - b.temperature()[n]));
        spread = std::max(spread, std::abs(a.temperature()[n] - (kC + 20.0)));
    }
    std::printf("     rotation invariance: worst %.3e K over a %.0f K field\n", worst, spread);
    expectTrue("the field is not trivial", spread > 300.0);
    expectTrue("rotating the mesh does not change the answer", worst < 1e-9);
}

void testRefusals() {
    ss::HexMesh mesh = makeBlock({1.0, 1.0, 1.0}, 2, 2, 2, 0.0);
    std::string why;

    th::Problem none;
    th::Solver s0;
    expectTrue("no mesh is refused", !s0.prepare(none, 300.0, &why) && !why.empty());

    // **A zero conductivity is a solve in which nothing conducts**, and every test
    // of it passes trivially. It has to be refused rather than return a field of
    // zeros that reads as a converged answer.
    th::Problem dead;
    dead.mesh = &mesh;
    dead.material = roundSteel();
    dead.material.conductivity = 0.0;
    th::Solver s1;
    expectTrue("a zero conductivity is refused", !s1.prepare(dead, 300.0, &why));
    dead.material.conductivity = 50.0;
    dead.material.specificHeat = 0.0;
    th::Solver s2;
    expectTrue("a zero specific heat is refused", !s2.prepare(dead, 300.0, &why));

    // Nothing holds the temperature: `K` alone is singular on the constant field.
    // A transient is fine -- the capacity regularises it -- and the steady solve is
    // not, and the difference has to be reported rather than returning whatever the
    // factorisation stumbled into.
    th::Problem floating;
    floating.mesh = &mesh;
    floating.material = roundSteel();
    th::Solver s3;
    expectTrue("a floating problem prepares", s3.prepare(floating, 300.0, &why));
    expectTrue("a floating transient step is legal", s3.step(1.0, &why));
    expectTrue("a floating steady solve is refused", !s3.solveSteady(&why));

    // A step that is not positive.
    th::Solver s4;
    expectTrue("a sane problem prepares", s4.prepare(floating, 300.0, &why));
    expectTrue("a zero step is refused", !s4.step(0.0, &why));
    expectTrue("a negative step is refused", !s4.step(-1.0, &why));

    // An inverted element.
    ss::HexMesh folded = mesh;
    for (int a = 0; a < 4; ++a) {
        const std::uint32_t lo = folded.index[static_cast<std::size_t>(a)];
        const std::uint32_t hi = folded.index[static_cast<std::size_t>(a) + 4];
        std::swap(folded.position[lo * 3 + 2], folded.position[hi * 3 + 2]);
    }
    th::Problem bent;
    bent.mesh = &folded;
    bent.material = roundSteel();
    th::Solver s5;
    expectTrue("an inverted element is refused", !s5.prepare(bent, 300.0, &why));

    // Mis-sized boundary data.
    th::Problem ragged;
    ragged.mesh = &mesh;
    ragged.material = roundSteel();
    ragged.prescribed.assign(3, 1u);
    ragged.prescribedValue.assign(3, 300.0);
    th::Solver s6;
    expectTrue("a short prescribed vector is refused", !s6.prepare(ragged, 300.0, &why));

    th::Problem badSource;
    badSource.mesh = &mesh;
    badSource.material = roundSteel();
    badSource.volumetricSource.assign(3, 1.0);
    th::Solver s7;
    expectTrue("a short source vector is refused", !s7.prepare(badSource, 300.0, &why));

    // A solver that was never prepared.
    th::Solver s8;
    expectTrue("an unprepared step is refused", !s8.step(1.0, &why));
    expectTrue("an unprepared steady solve is refused", !s8.solveSteady(&why));

    // **A solver whose `prepare` failed must not be usable.** Every refusal above
    // returns false, and a caller that ignores the return value must still not get
    // an answer computed on a half-built solver -- the `World::load` failure mode
    // `CLAUDE.md` records, arriving somewhere else.
    th::Solver s8b;
    expectTrue("this prepare fails", !s8b.prepare(dead, 300.0, &why));
    expectTrue("a step after a failed prepare is refused", !s8b.step(1.0, &why));
    expectTrue("a steady solve after a failed prepare is refused", !s8b.solveSteady(&why));

    ss::HexMesh empty;
    empty.position.assign(24, 0.0);
    th::Problem noElements;
    noElements.mesh = &empty;
    noElements.material = roundSteel();
    th::Solver s8c;
    expectTrue("a mesh with nodes but no elements is refused",
               !s8c.prepare(noElements, 300.0, &why));

    th::Problem sane;
    sane.mesh = &mesh;
    sane.material = roundSteel();
    th::Solver s8d;
    expectTrue("an initial field of the wrong length is refused",
               !s8d.prepare(sane, std::vector<double>(mesh.nodeCount() - 1, 300.0), &why));

    th::Problem weightless = sane;
    weightless.material.density = 0.0;
    th::Solver s8e;
    expectTrue("a zero density is refused", !s8e.prepare(weightless, 300.0, &why));

    // Only *one* of the two prescribed vectors the right length -- the pair has to
    // be checked as a pair.
    th::Problem halfPrescribed = sane;
    halfPrescribed.prescribed.assign(mesh.nodeCount(), 0u);
    halfPrescribed.prescribedValue.assign(mesh.nodeCount() - 2, 300.0);
    th::Solver s8f;
    expectTrue("a prescribed value vector of the wrong length is refused",
               !s8f.prepare(halfPrescribed, 300.0, &why));

    // **Every node prescribed** is a legal problem with nothing to solve, and it
    // has to run rather than be refused or silently do something else.
    th::Problem allHeld = sane;
    for (std::size_t n = 0; n < mesh.nodeCount(); ++n)
        hold(allHeld, mesh, n, kC + 400.0 + 10.0 * static_cast<double>(n));
    th::Solver s8g;
    expectTrue("a fully prescribed problem prepares: " + why, s8g.prepare(allHeld, 300.0, &why));
    expectEqual("it has no free nodes", static_cast<long long>(s8g.freeNodes()), 0LL);
    expectTrue("and a step on it runs: " + why, s8g.step(1.0, &why));
    double worstHeld = 0.0;
    for (std::size_t n = 0; n < mesh.nodeCount(); ++n)
        worstHeld = std::max(worstHeld,
                             std::abs(s8g.temperature()[n] -
                                      (kC + 400.0 + 10.0 * static_cast<double>(n))));
    expectTrue("every node stays where it was held", worstHeld == 0.0);

    // A single element, which every mesh above has more than one of.
    ss::HexMesh one = ss::makePlateMesh(0.05, 0.05, 0.012, 1, 1, 1);
    th::Problem alone;
    alone.mesh = &one;
    alone.material = roundSteel();
    for (std::size_t n = 0; n < one.nodeCount(); ++n)
        if (one.position[n * 3] < 1e-12) hold(alone, one, n, kC + 500.0);
    th::Solver s8h;
    expectTrue("a one-element problem prepares: " + why, s8h.prepare(alone, kC + 20.0, &why));
    expectTrue("and solves steady: " + why, s8h.solveSteady(&why));
    for (std::size_t n = 0; n < one.nodeCount(); ++n)
        expectNear("one element with one face held is isothermal", s8h.temperature()[n],
                   kC + 500.0, 1e-9);
}

// --- 11a. Four holes mutation testing found ------------------------------------------
//
// Every test in this section exists because a mutant survived the ones above it.
// They are grouped rather than scattered because what they have in common is
// worth stating: each is a quantity the earlier tests measured only in *total*,
// where the defect was in the *distribution*.

void testBoundaryOrientation() {
    // **The divergence theorem.** `integral x . n dA = 3 V` over a closed surface,
    // so the boundary's areas, its centroids and the *sign* of its normals are all
    // tied to a volume this file computes a completely different way, through
    // `solidshell::gaussVolumes`.
    //
    // This exists because flipping every outward normal survived the suite: the
    // tests above select faces by `normal.z > 0.5`, and reversing all of them
    // merely swaps which face gets the fire on a mesh that is symmetric about the
    // swap. A test that reads the normals only to sort faces cannot see their sign.
    for (int variant = 0; variant < 2; ++variant) {
        ss::HexMesh mesh = variant == 0 ? makeBlock({1.3, 0.9, 0.7}, 3, 2, 2, 0.0)
                                        : makeBar(0.6, 0.05, 5, 1.3);
        double flux = 0.0, area = 0.0;
        for (const th::BoundaryFace& f : th::boundaryFaces(mesh)) {
            flux += dot(f.centroid, f.normal) * f.area;
            area += f.area;
        }
        const double volume = meshVolume(mesh);
        expectTrue("the boundary has area", area > 0.1);
        expectNear("the outward normals satisfy the divergence theorem", flux, 3.0 * volume,
                   1e-9 * volume);
        // **Vacuity guard.** The identity has to be able to tell the two signs
        // apart, which it cannot if the volume is near zero.
        expectTrue("and 3V is far from -3V", volume > 1e-3);
    }
}

void testNodeOrdering() {
    // **Cuthill-McKee has to actually win somewhere.** Nothing above builds a mesh
    // whose natural numbering is bad, so every mutant that disabled the reordering
    // -- or that mis-measured one of the two candidate bandwidths -- survived by
    // falling back on a numbering that was already good.
    //
    // The same mesh is solved twice, once with its nodes shuffled. The answer must
    // be identical to rounding and the bandwidth must be far smaller than the
    // shuffle left it.
    ss::HexMesh straight = makeBlock({1.2, 0.6, 0.4}, 6, 4, 3, 0.0);
    ss::HexMesh shuffled = straight;
    std::vector<std::uint32_t> to(straight.nodeCount());
    for (std::uint32_t n = 0; n < to.size(); ++n) to[n] = n;
    std::mt19937 rng(20260807u);
    std::shuffle(to.begin(), to.end(), rng);
    for (std::size_t n = 0; n < straight.nodeCount(); ++n)
        for (int i = 0; i < 3; ++i)
            shuffled.position[to[n] * 3 + static_cast<std::size_t>(i)] =
                straight.position[n * 3 + static_cast<std::size_t>(i)];
    for (std::size_t i = 0; i < straight.index.size(); ++i) shuffled.index[i] = to[straight.index[i]];

    const auto solve = [&](const ss::HexMesh& mesh, const std::vector<std::uint32_t>* map) {
        th::Problem problem;
        problem.mesh = &mesh;
        problem.material = roundSteel();
        for (std::size_t n = 0; n < straight.nodeCount(); ++n) {
            const double x = straight.position[n * 3];
            const std::size_t target = map ? (*map)[n] : n;
            if (x < 1e-12) hold(problem, mesh, target, kC + 700.0);
            if (x > 1.2 - 1e-12) hold(problem, mesh, target, kC + 20.0);
        }
        th::Solver solver;
        std::string why;
        expectTrue("the ordering problem solves: " + why,
                   solver.prepare(problem, kC + 20.0, &why) && solver.solveSteady(&why));
        std::vector<double> field(straight.nodeCount());
        for (std::size_t n = 0; n < straight.nodeCount(); ++n)
            field[n] = solver.temperature()[map ? (*map)[n] : n];
        return std::make_tuple(field, solver.halfBandwidth(), solver.freeNodes(),
                               solver.account().equilibriumResidual);
    };
    const auto [fieldA, bandA, freeA, residualA] = solve(straight, nullptr);
    const auto [fieldB, bandB, freeB, residualB] = solve(shuffled, &to);

    // What the shuffled mesh's *own* numbering would have delivered, derived here
    // rather than remembered: free nodes take slots in ascending node index, and
    // the band is the widest element's spread over them.
    const auto identityBand = [&](const ss::HexMesh& mesh, const std::vector<std::uint32_t>* map) {
        std::vector<std::ptrdiff_t> slot(mesh.nodeCount(), -1);
        std::vector<bool> free(mesh.nodeCount(), true);
        for (std::size_t n = 0; n < straight.nodeCount(); ++n) {
            const double x = straight.position[n * 3];
            if (x < 1e-12 || x > 1.2 - 1e-12) free[map ? (*map)[n] : n] = false;
        }
        std::ptrdiff_t next = 0;
        for (std::size_t n = 0; n < mesh.nodeCount(); ++n)
            if (free[n]) slot[n] = next++;
        std::size_t width = 0;
        for (std::size_t e = 0; e < mesh.elementCount(); ++e) {
            std::ptrdiff_t lo = -1, hi = -1;
            for (int a = 0; a < 8; ++a) {
                const std::ptrdiff_t d = slot[mesh.index[e * 8 + static_cast<std::size_t>(a)]];
                if (d < 0) continue;
                if (lo < 0 || d < lo) lo = d;
                if (d > hi) hi = d;
            }
            if (lo >= 0) width = std::max(width, static_cast<std::size_t>(hi - lo));
        }
        return width;
    };
    const std::size_t rawA = identityBand(straight, nullptr);
    const std::size_t rawB = identityBand(shuffled, &to);

    std::printf("     node ordering: half-bandwidth %zu in the mesh's own order (its raw band is"
                " %zu), %zu after a random shuffle (raw %zu), over %zu free nodes\n",
                bandA, rawA, bandB, rawB, freeA);
    expectEqual("the same problem has the same free count", static_cast<long long>(freeB),
                static_cast<long long>(freeA));
    // The solver has to keep the *narrower* of the two candidates it builds, which
    // means never doing worse than the mesh arrived with, and much better when the
    // mesh arrived badly numbered.
    expectTrue("the shuffle really did wreck the natural order", rawB > 3 * rawA);
    expectTrue("Cuthill-McKee rescues it", bandB * 2 < rawB);
    expectTrue("and never does worse than the numbering it was given",
               bandA <= rawA && bandB <= rawB);

    double worst = 0.0, spread = 0.0;
    for (std::size_t n = 0; n < fieldA.size(); ++n) {
        worst = std::max(worst, std::abs(fieldA[n] - fieldB[n]));
        spread = std::max(spread, std::abs(fieldA[n] - (kC + 20.0)));
    }
    expectTrue("the field is not trivial", spread > 300.0);
    expectTrue("the node ordering does not change the answer", worst < 1e-9);
    // A band one element too narrow drops terms silently; both orderings must show
    // the same, machine-level, free-node residual.
    expectTrue("neither ordering dropped a term", residualA < 1e-6 && residualB < 1e-6);
}

void testVolumetricSourceProfile() {
    // **A uniform source between two held ends gives a parabola**,
    // `T = T0 + Q x (L - x) / 2k`, and in one dimension with linear elements the
    // *nodal* values of the Galerkin solution are exact.
    //
    // This exists because `testEnergyAccount` only ever checked the source's
    // *total*: every mutant that scattered the consistent source load onto the
    // wrong nodes still deposited the right number of joules, so the account
    // closed and the profile was never looked at.
    const double length = 0.4, width = 0.02, q = 5.0e6;
    const double hold0 = kC + 20.0;
    ss::HexMesh mesh = makeBar(length, width, 16);
    th::Problem problem;
    problem.mesh = &mesh;
    problem.material = roundSteel();
    problem.volumetricSource.assign(mesh.elementCount(), q);
    for (std::size_t n = 0; n < mesh.nodeCount(); ++n) {
        const double x = mesh.position[n * 3];
        if (x < 1e-12 || x > length - 1e-12) hold(problem, mesh, n, hold0);
    }

    th::Solver solver;
    std::string why;
    expectTrue("the source bar prepares: " + why, solver.prepare(problem, hold0, &why));
    expectTrue("the source bar solves: " + why, solver.solveSteady(&why));

    double worst = 0.0, peak = 0.0;
    for (int i = 0; i <= 16; ++i) {
        const double x = mesh.position[static_cast<std::size_t>(i) * 4 * 3];
        const double exact = hold0 + q * x * (length - x) / (2.0 * problem.material.conductivity);
        peak = std::max(peak, exact - hold0);
        worst = std::max(worst,
                         std::abs(solver.temperature()[static_cast<std::size_t>(i) * 4] - exact));
    }
    std::printf("     volumetric source: a %.0f K parabola reproduced to %.3e K at the nodes\n",
                peak, worst);
    expectTrue("the parabola has a real peak", peak > 300.0);
    expectTrue("a uniform source gives the exact parabola at the nodes", worst < 1e-9);

    // And the account: the source is the only channel in, the held ends the only
    // one out, and at steady state they cancel exactly.
    const double total = q * length * width * width;
    expectNear("the source's total power", total, q * meshVolume(mesh), 1e-9 * total);
}

void testFin() {
    // **A fin.** A bar held hot at one end and losing heat by convection along its
    // sides settles at `theta = theta0 cosh(m(L-x)) / cosh(mL)` with
    // `m = sqrt(hP/kA)` -- valid where the Biot number across the section is small
    // and the tip is adiabatic.
    //
    // This exists because nothing above made a *convective surface* carry a
    // varying temperature. `testNewtonCooling` and the mixed-film check both
    // depend on the film matrix only through its total `h A`, and the row sums of
    // `integral N_a N_b dA` are independent of where the Gauss points sit -- so
    // every mutation of the face quadrature's sample locations survived. A fin's
    // sides are cooled by a temperature that varies along them, which is exactly
    // the term those mutations break.
    const double length = 0.30, side = 0.010, h = 50.0;
    const double base = kC + 700.0, ambient = kC + 20.0;
    ss::HexMesh mesh = makeBar(length, side, 60);

    th::Problem problem;
    problem.mesh = &mesh;
    problem.material = roundSteel();
    th::Film sides;
    sides.coefficient = h;
    sides.ambient = ambient;
    double sideArea = 0.0;
    for (const th::BoundaryFace& f : th::boundaryFaces(mesh)) {
        if (std::abs(f.normal.x) > 0.5) continue;  // base and tip stay out of it
        sides.face.push_back(f);
        sideArea += f.area;
    }
    problem.film = {sides};
    for (std::size_t n = 0; n < mesh.nodeCount(); ++n)
        if (mesh.position[n * 3] < 1e-12) hold(problem, mesh, n, base);

    const double perimeter = 4.0 * side, section = side * side;
    const double m = std::sqrt(h * perimeter / (problem.material.conductivity * section));
    const double biot = h * 0.5 * side / problem.material.conductivity;
    expectNear("the cooled sides are the bar's four flanks", sideArea, perimeter * length, 1e-9);
    expectTrue("the section's Biot number is small", biot < 0.02);
    expectTrue("the fin is long enough to decay", m * length > 4.0);

    th::Solver solver;
    std::string why;
    expectTrue("the fin prepares: " + why, solver.prepare(problem, ambient, &why));
    expectTrue("the fin solves: " + why, solver.solveSteady(&why));

    double worst = 0.0;
    for (int i = 0; i <= 60; ++i) {
        const double x = mesh.position[static_cast<std::size_t>(i) * 4 * 3];
        const double exact = ambient + (base - ambient) * std::cosh(m * (length - x)) /
                                           std::cosh(m * length);
        worst = std::max(worst,
                         std::abs(solver.temperature()[static_cast<std::size_t>(i) * 4] - exact));
    }
    const double heatOut = problem.material.conductivity * section * m * (base - ambient) *
                           std::tanh(m * length);
    std::printf("     fin: mL = %.2f, Bi = %.4f -- worst %.3f K against cosh, %.1f W leaving a"
                " %.0f K base\n",
                m * length, biot, worst, heatOut, base - ambient);
    // The fin solution is a *continuum* answer and the mesh discretises it, so this
    // is a convergence tolerance rather than an identity: 60 elements at
    // `mL = 6` is `m h = 0.1`, and the second-order element gives about
    // `(m h)^2 / 12` of the 680 K base, which is 0.6 K. Asserted at twice that,
    // and bracketed below so a mutant that made the fin *too* good is caught too.
    expectTrue("the fin matches its closed form", worst < 1.5);
    expectTrue("and the fin actually decayed",
               solver.temperature()[60 * 4] - ambient < 0.02 * (base - ambient));

    // A second closed form off the same fixture: the heat entering the base is
    // `k A m theta0 tanh(mL)`. It is taken from the *account* rather than from a
    // one-sided gradient at the base -- a difference over the first element misses
    // that element's own convective loss, which is `h P (dx/2) theta0` and here is
    // 5% of the answer, so it would have needed a tolerance five times looser than
    // the quantity it was checking.
    //
    // The field is already steady, so a long transient step changes nothing and
    // the account reports the steady power. That the two channels are equal and
    // opposite is exact -- it is the `1^T K = 0` identity -- and that either
    // equals the fin's closed form is the physics.
    const double before = solver.account().prescribedHeat, beforeFilm = solver.account().filmHeat;
    expectTrue("a step off the steady state runs: " + why, solver.step(1e4, &why));
    const double basePower = (solver.account().prescribedHeat - before) / 1e4;
    const double filmPower = (solver.account().filmHeat - beforeFilm) / 1e4;
    std::printf("     fin heat: %.4f W in at the base, %.4f W out through the sides,"
                " %.4f W from k A m theta tanh(mL)\n",
                basePower, -filmPower, heatOut);
    expectNear("what enters the base leaves through the sides", basePower, -filmPower,
               1e-9 * heatOut);
    // Tightened to what was measured, with the reason: the element is second order
    // and `m h` is 0.1 here, so `(m h)^2 / 24` is 4.2e-4 -- which is the 0.04% the
    // solve delivers, to one figure.
    const double allowance = heatOut * (m * mesh.position[4 * 3]) * (m * mesh.position[4 * 3]) /
                             24.0 * 1.5;
    expectNear("and it is the fin's closed-form heat rate", basePower, heatOut, allowance);
    expectTrue("the fin allowance is not a licence", allowance < 1e-3 * heatOut);
}

void testWarpedFilmFace() {
    // A film face that is **not planar**. Every boundary face above is flat, and a
    // flat face's quadrature is insensitive to where its sample points sit -- so
    // the 2x2 rule's own `1/sqrt(3)` survived every mutation. On a warped face it
    // does not: the area element varies over the face and the rule has to be the
    // right rule.
    //
    // The area is checked against a 64x64 midpoint sum of the same bilinear
    // surface, computed here, which is an independent derivation and not a
    // remembered number.
    ss::HexMesh mesh = ss::makePlateMesh(0.4, 0.3, 0.012, 1, 1, 1);
    // Lift one corner of the +zeta face well out of plane.
    const std::uint32_t lifted = mesh.index[6];  // node 6 of element 0
    mesh.position[lifted * 3 + 2] += 0.09;

    const std::vector<th::BoundaryFace> faces = th::boundaryFaces(mesh);
    const th::BoundaryFace* top = nullptr;
    for (const th::BoundaryFace& f : faces)
        if (f.face == 1) top = &f;
    expectTrue("the warped face is on the boundary", top != nullptr);
    if (top == nullptr) return;

    Vec3 corner[4];
    for (int i = 0; i < 4; ++i) corner[i] = nodeAt(mesh, top->node[i]);
    double fine = 0.0;
    const int n = 400;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            const double u = (i + 0.5) / n * 2.0 - 1.0, v = (j + 0.5) / n * 2.0 - 1.0;
            Vec3 du{}, dv{};
            static constexpr double c[4][2] = {{-1, -1}, {+1, -1}, {+1, +1}, {-1, +1}};
            for (int k = 0; k < 4; ++k) {
                du += corner[k] * (0.25 * c[k][0] * (1.0 + v * c[k][1]));
                dv += corner[k] * (0.25 * c[k][1] * (1.0 + u * c[k][0]));
            }
            fine += length(cross(du, dv)) * (2.0 / n) * (2.0 / n);
        }
    std::printf("     warped face: 2x2 Gauss gives %.9f m^2, a 400x400 midpoint sum %.9f m^2"
                " (flat would be %.9f)\n",
                top->area, fine, 0.4 * 0.3);
    // **Vacuity guard.** A face that is barely warped is a flat face and proves
    // nothing about the quadrature.
    expectTrue("the face really is warped", top->area > 1.02 * 0.4 * 0.3);
    // The 2x2 rule is not exact on the square root of a quartic, so this is a
    // quadrature comparison and not an identity; a wrong sample location moves it
    // by percent, not by a part in a thousand.
    expectNear("the warped face's area", top->area, fine, 3e-5 * fine);
}

// --- 11b. The paths the tests above do not reach ---------------------------------------
//
// Written by going through the diff and asking what in it is *not* exercised,
// rather than by asking whether the tests pass -- the habit `CLAUDE.md` records
// against a commit whose headline feature was well tested and whose two new
// functions on the caller's own path had no test at all.

void testUncovered() {
    ss::HexMesh mesh = makeBlock({0.5, 0.4, 0.012}, 4, 3, 1, 0.1);
    std::string why;

    // **The non-uniform initial field.** Every test above starts the solver from a
    // single temperature, so the vector overload of `prepare` was untested -- and
    // it is the one a restart or a coupled solve uses.
    //
    // Started from a linear field with no boundary condition at all, conduction
    // must level it to the *volume-averaged* temperature and conserve every joule
    // doing so. Both are closed forms: the mean is the enthalpy divided by the
    // heat capacity, and the account must show no heat crossing any boundary.
    th::Problem problem;
    problem.mesh = &mesh;
    problem.material = roundSteel();

    std::vector<double> initial(mesh.nodeCount());
    for (std::size_t n = 0; n < mesh.nodeCount(); ++n)
        initial[n] = kC + 20.0 + 900.0 * mesh.position[n * 3] / 0.5;

    th::Solver solver;
    expectTrue("a non-uniform start prepares: " + why, solver.prepare(problem, initial, &why));
    const double startEnthalpy = solver.account().enthalpy;
    // The slowest mode of an insulated bar decays as `exp(-alpha pi^2 t / L^2)`, a
    // time constant of 2027 s here, so the run has to be tens of those before
    // "level" means anything -- 8000 s left 12 K of the original 900 K standing,
    // which is the first version of this test passing nothing.
    const double slowest = 0.5 * 0.5 / (kRoundAlpha * 3.14159265358979 * 3.14159265358979);
    for (int i = 0; i < 300; ++i) expectTrue("a levelling step runs: " + why, solver.step(200.0, &why));
    expectTrue("the run is long against the slowest mode", solver.time() > 25.0 * slowest);

    double lo = 1e30, hi = -1e30;
    for (double v : solver.temperature()) {
        lo = std::min(lo, v);
        hi = std::max(hi, v);
    }
    const double capacity =
        problem.material.density * problem.material.specificHeat * meshVolume(mesh);
    const double mean = kC + 20.0 + startEnthalpy / capacity;
    std::printf("     insulated levelling: %.1f K spread became %.3e K, settling at %.4f K"
                " against a mean of %.4f K; the account moved %.3e J\n",
                900.0, hi - lo, 0.5 * (lo + hi), mean, solver.account().enthalpyChange);
    expectTrue("it started far from level", startEnthalpy > 1e6);
    expectTrue("an insulated body levels", hi - lo < 1e-6);
    expectNear("and levels at the volume mean", 0.5 * (lo + hi), mean, 1e-9);
    // Nothing crosses the boundary, so the enthalpy is a conserved quantity and not
    // merely a balanced one.
    expectTrue("an insulated body conserves its enthalpy exactly",
               std::abs(solver.account().enthalpyChange) < 1e-9 * startEnthalpy);
    expectEqual("no heat channel opened", static_cast<long long>(solver.account().filmHeat != 0.0 ||
                                                                solver.account().prescribedHeat != 0.0 ||
                                                                solver.account().sourceHeat != 0.0),
                0LL);

    // **The clamped ends of the material curves.** Below 20 C and above 1200 C the
    // standard says nothing, so `c` holds its end value and `h` continues linearly.
    // Both branches are outside every solve above.
    expectNear("c is clamped below 20 C", th::carbonSteelSpecificHeat(kC - 50.0),
               th::carbonSteelSpecificHeat(kC + 20.0), 1e-12);
    expectNear("c is clamped above 1200 C", th::carbonSteelSpecificHeat(kC + 1500.0), 650.0, 1e-12);
    expectNear("k is clamped below 20 C", th::carbonSteelConductivity(kC - 50.0),
               th::carbonSteelConductivity(kC + 20.0), 1e-12);
    expectNear("h continues linearly below 20 C", th::carbonSteelEnthalpy(kC - 30.0),
               -50.0 * th::carbonSteelSpecificHeat(kC + 20.0), 1e-9);
    expectNear("h continues linearly above 1200 C",
               th::carbonSteelEnthalpy(kC + 1400.0) - th::carbonSteelEnthalpy(kC + 1200.0),
               200.0 * 650.0, 1e-9);

    // **`BoundaryFace::centroid` and the single-element `explicitLimit`**, neither
    // of which any test above reads.
    const std::vector<th::BoundaryFace> faces = th::boundaryFaces(mesh);
    double area = 0.0;
    Vec3 moment{};
    for (const th::BoundaryFace& f : faces) {
        area += f.area;
        moment += f.centroid * f.area;
    }
    moment *= 1.0 / area;
    // The block's surface is symmetric about its own middle, so the area-weighted
    // centroid of the boundary is the block's centre. A face whose centroid was
    // taken from the wrong four nodes would move it.
    expectNear("the boundary's area centroid is the block's centre, x", moment.x, 0.25, 1e-12);
    expectNear("the boundary's area centroid is the block's centre, y", moment.y, 0.20, 1e-12);
    expectNear("the boundary's area centroid is the block's centre, z", moment.z, 0.006, 1e-12);

    double one[ss::kDof];
    mesh.gather(0, mesh.position, one);
    expectTrue("the single-element explicit limit agrees with the mesh's",
               std::abs(th::explicitLimit(one, roundSteel()) -
                        th::explicitLimit(mesh, roundSteel())) <
                   1e-9 * th::explicitLimit(mesh, roundSteel()) + 1e-6);

    // **A film face that does not belong to the mesh it was applied to.** Silently
    // integrating it would put heat into whichever nodes the stale indices name.
    th::Problem stale;
    stale.mesh = &mesh;
    stale.material = roundSteel();
    th::Film wrong;
    wrong.flux = 1000.0;
    th::BoundaryFace forged = faces.front();
    forged.node[2] += 1;
    wrong.face.push_back(forged);
    stale.film = {wrong};
    th::Solver s9;
    expectTrue("a film face that does not match the mesh is refused",
               !s9.prepare(stale, 300.0, &why));

    th::Film absent = wrong;
    absent.face.front() = faces.front();
    absent.face.front().element = static_cast<std::uint32_t>(mesh.elementCount());
    stale.film = {absent};
    th::Solver s10;
    why.clear();
    expectTrue("a film face naming an element that is not there is refused",
               !s10.prepare(stale, 300.0, &why));
    // For the *stated* reason. An off-by-one in the bound lets the index through
    // and the node-match check refuses it instead -- same answer, different code
    // path, and an out-of-bounds read on the way.
    expectTrue("and says so: " + why, why.find("names element") != std::string::npos);

    th::Film sixth = wrong;
    sixth.face.front() = faces.front();
    sixth.face.front().face = 6;  // a hexahedron has six, numbered 0..5
    stale.film = {sixth};
    th::Solver s10b;
    expectTrue("a film face naming a seventh face is refused",
               !s10b.prepare(stale, 300.0, &why));

    // An element naming a node the mesh does not have reads out of bounds in every
    // loop in the file, so it is refused rather than left to a sanitizer.
    ss::HexMesh broken = mesh;
    broken.index[5] = static_cast<std::uint32_t>(broken.nodeCount());
    th::Problem outside;
    outside.mesh = &broken;
    outside.material = roundSteel();
    th::Solver s10c;
    why.clear();
    expectTrue("an element naming a node past the end is refused",
               !s10c.prepare(outside, 300.0, &why));
    expectTrue("and says so rather than calling the element degenerate: " + why,
               why.find("names node") != std::string::npos);

    // **A film carrying flux and convection at once**, which is what a real fire
    // boundary is -- radiation delivered as a flux plus a convective term. The two
    // are additive by construction, and the check is that the steady state is the
    // one the algebra gives: `T = ambient + flux / coefficient`, everywhere.
    th::Problem both;
    both.mesh = &mesh;
    both.material = roundSteel();
    th::Film mixed;
    mixed.flux = 9000.0;
    mixed.coefficient = 30.0;
    mixed.ambient = kC + 20.0;
    for (const th::BoundaryFace& f : faces) mixed.face.push_back(f);
    both.film = {mixed};
    th::Solver s11;
    expectTrue("a mixed film prepares: " + why, s11.prepare(both, kC + 20.0, &why));
    expectTrue("a mixed film has a steady state: " + why, s11.solveSteady(&why));
    double worst = 0.0;
    for (double v : s11.temperature())
        worst = std::max(worst, std::abs(v - (mixed.ambient + mixed.flux / mixed.coefficient)));
    expectTrue("flux and convection on one surface balance where the algebra says",
               worst < 1e-9);
    expectTrue("and that is not the ambient temperature", mixed.flux / mixed.coefficient > 100.0);
}

// --- 12. What it costs on a real patch --------------------------------------------------------

void testCost() {
    // One frame bay of the ferry's side shell: 2.4 m by 2.4 m of 12 mm plating at
    // the 50 mm element size `zone.hpp` meshes a crush zone with. Not the ferry's
    // own mesh -- this file has no business reaching into the zone builder -- but
    // the same element count and the same aspect ratio.
    ss::HexMesh mesh = ss::makePlateMesh(2.4, 2.4, 0.012, 48, 48, 1);
    th::Problem problem;
    problem.mesh = &mesh;
    problem.material = ah36Steel();

    th::Film fire;
    fire.flux = 45000.0;
    for (const th::BoundaryFace& f : th::boundaryFaces(mesh))
        if (f.normal.z > 0.5) fire.face.push_back(f);
    problem.film = {fire};
    for (std::size_t n = 0; n < mesh.nodeCount(); ++n)
        if (mesh.position[n * 3] < 1e-12) hold(problem, mesh, n, kC + 20.0);

    using Clock = std::chrono::steady_clock;
    th::Solver solver;
    std::string why;
    const auto t0 = Clock::now();
    expectTrue("the patch prepares: " + why, solver.prepare(problem, kC + 20.0, &why));
    const auto t1 = Clock::now();
    expectTrue("the first step runs: " + why, solver.step(1.0, &why));
    const auto t2 = Clock::now();
    const int steps = 200;
    for (int i = 0; i < steps; ++i) expectTrue("a patch step runs: " + why, solver.step(1.0, &why));
    const auto t3 = Clock::now();

    const auto ms = [](auto a, auto b) {
        return std::chrono::duration<double, std::milli>(b - a).count();
    };
    std::printf("     cost: %zu elements, %zu nodes, half-bandwidth %zu -- prepare %.1f ms,"
                " first step %.1f ms, %d further steps %.2f ms each\n",
                mesh.elementCount(), mesh.nodeCount(), solver.halfBandwidth(), ms(t0, t1),
                ms(t1, t2), steps, ms(t2, t3) / steps);
    // The honest comparison, at this mesh: explicit is *viable* here, and the case
    // for implicit is what happens under the refinement the coupling will want.
    const double limit = th::explicitLimit(mesh, problem.material);
    ss::HexMesh resolved = ss::makePlateMesh(2.4, 2.4, 0.012, 48, 48, 4);
    const double refined = th::explicitLimit(resolved, problem.material);
    std::printf("     explicit limit here %.3f s, so %d s of fire is %.0f explicit steps against"
                " %d implicit -- explicit is viable on this mesh. Four elements through the"
                " plate takes it to %.3f s and %.0f explicit steps.\n",
                limit, steps + 1, (steps + 1) / limit, steps + 1, refined,
                (steps + 1) / refined);

    // The whole point of a reused factorisation: 201 steps, one factorisation.
    //
    // **Asserted on the count and not on the clock.** The first version of this
    // compared the two wall times and it failed once, under load, on a machine
    // sharing three other jobs -- a ratio between two timings is a statement about
    // the box and not about the code, and a test that is sometimes red is worse
    // than no test. The times are printed because they are the answer to "what
    // does it cost"; the count is what carries the claim.
    expectEqual("a constant-property run at a fixed step factors once", solver.factorisations(), 1);

    expectTrue("the patch heated", solver.temperature()[mesh.nodeCount() - 1] > kC + 100.0);
    const th::Account& a = solver.account();
    expectTrue("the patch account closes",
               std::abs(a.residual()) < 1e-12 * std::abs(a.enthalpyChange));
}

// ================================================================================
// Hot steel: EN 1993-1-2 §3.2, and what it does to a structure
// ================================================================================
//
// The conduction solve above says what the steel's temperature is. Everything
// below says what the steel then carries, and the tests divide into three kinds:
// the standard reproduced exactly, the coupling shown to be exact at 20 C, and a
// structural consequence with hand-computed ends.

// --- 1. The standard's own table ------------------------------------------------

// **This is a second, independent transcription of EN 1993-1-2:2005 Table 3.1.**
// `thermal.cpp` has its own and does not export it, deliberately: a test that read
// the implementation's array would move with it, and a mutant that changed one row
// would change both sides and survive. Two transcriptions is the whole point, and
// it is the reason the numbers below are written out rather than generated.
struct TabulatedReduction {
    double celsius, yield, proportional, modulus;
};
constexpr TabulatedReduction kTable32[] = {
    {  20.0, 1.000, 1.0000, 1.0000},
    { 100.0, 1.000, 1.0000, 1.0000},
    { 200.0, 1.000, 0.8070, 0.9000},
    { 300.0, 1.000, 0.6130, 0.8000},
    { 400.0, 1.000, 0.4200, 0.7000},
    { 500.0, 0.780, 0.3600, 0.6000},
    { 600.0, 0.470, 0.1800, 0.3100},
    { 700.0, 0.230, 0.0750, 0.1300},
    { 800.0, 0.110, 0.0500, 0.0900},
    { 900.0, 0.060, 0.0375, 0.0675},
    {1000.0, 0.040, 0.0250, 0.0450},
    {1100.0, 0.020, 0.0125, 0.0225},
    {1200.0, 0.000, 0.0000, 0.0000},
};
constexpr int kTableRows = static_cast<int>(sizeof(kTable32) / sizeof(kTable32[0]));

void testReductionTable() {
    // **Every tabulated point exactly -- and the three where "exactly" is not
    // available, with the reason.** The interpolant takes a station branch rather
    // than trusting `a + (b - a) * 1.0` to land on `b`, so given the exact Celsius
    // station it returns the standard's literal. What it is given is
    // `celsius(kelvin)`, and **the Celsius-to-Kelvin round trip is not exact at
    // 800, 900 and 1000 C**: the nearest double to 1073.15 differs from
    // `800 + 273.15` by 1.1e-13 K, and no double at all satisfies
    // `k - 273.15 == 800` because the spacing at 1073 is twice the spacing at 800.
    // That is a property of the Kelvin interface `thermal.hpp` chose, not of the
    // table, and the resulting error is bounded by the steepest slope in Table 3.1
    // -- 0.0031 per K, between 500 and 600 C -- times 1.1e-13 K, so 3.5e-16, which
    // is three units in the last place of the factor.
    //
    // So the test asserts bit-equality wherever the round trip is exact, counts
    // those so the check cannot quietly become vacuous, and holds the remaining
    // three to that derived bound rather than loosening all thirteen to it.
    int exactStations = 0;
    for (const TabulatedReduction& row : kTable32) {
        const double k = row.celsius + kC;
        const th::SteelReduction r = th::carbonSteelReduction(k);
        const bool roundTrips = (k - kC) == row.celsius;
        const double tolerance = roundTrips ? 0.0 : 4e-16;
        if (roundTrips) ++exactStations;
        expectNear("k_y at a tabulated temperature is the tabulated value", r.effectiveYield,
                   row.yield, tolerance);
        expectNear("k_p at a tabulated temperature is the tabulated value",
                   r.proportionalLimit, row.proportional, tolerance);
        expectNear("k_E at a tabulated temperature is the tabulated value", r.youngsModulus,
                   row.modulus, tolerance);
        // The three scalar accessors are the same numbers, bit for bit. They are
        // separate entry points and a caller may use either; nothing but this
        // stops them drifting.
        expectTrue("carbonSteelYieldFactor agrees with the struct bit for bit",
                   th::carbonSteelYieldFactor(k) == r.effectiveYield);
        expectTrue("carbonSteelProportionalFactor agrees with the struct bit for bit",
                   th::carbonSteelProportionalFactor(k) == r.proportionalLimit);
        expectTrue("carbonSteelModulusFactor agrees with the struct bit for bit",
                   th::carbonSteelModulusFactor(k) == r.youngsModulus);
    }
    expectEqual("ten of the thirteen stations are reached exactly through kelvin",
                exactStations, 10);

    // **Linear between them, and monotone.** Linearity is asserted at the midpoint
    // of every interval against the mean of its ends -- which is what linear
    // interpolation *means*, and is a different statement from "it is between
    // them". A quadratic through the same points passes a betweenness check and
    // fails this one.
    for (int i = 1; i < kTableRows; ++i) {
        const double lo = kTable32[i - 1].celsius, hi = kTable32[i].celsius;
        const double mid = 0.5 * (lo + hi);
        const th::SteelReduction r = th::carbonSteelReduction(mid + kC);
        expectNear("k_y at the midpoint is the mean of the ends", r.effectiveYield,
                   0.5 * (kTable32[i - 1].yield + kTable32[i].yield), 1e-15);
        expectNear("k_p at the midpoint is the mean of the ends", r.proportionalLimit,
                   0.5 * (kTable32[i - 1].proportional + kTable32[i].proportional), 1e-15);
        expectNear("k_E at the midpoint is the mean of the ends", r.youngsModulus,
                   0.5 * (kTable32[i - 1].modulus + kTable32[i].modulus), 1e-15);

        // And at a quarter and three quarters, so an implementation that split the
        // interval in half and recursed would still have to be linear.
        for (double s : {0.25, 0.75}) {
            const double t = lo + s * (hi - lo);
            const th::SteelReduction q = th::carbonSteelReduction(t + kC);
            expectNear("k_y interpolates linearly across the interval", q.effectiveYield,
                       kTable32[i - 1].yield + s * (kTable32[i].yield - kTable32[i - 1].yield),
                       1e-15);
        }
    }

    // **Why the interpolant's `t == station` branch is currently unobservable**,
    // asserted rather than assumed. `s` is exactly 1 at the upper station, so the
    // general expression reduces to `a + (b - a)` -- and for all thirty-six
    // adjacent pairs of Table 3.1 that is bit-exactly `b`. Mutation testing found
    // the branch dead and found the comment that justified it inventing a unit in
    // the last place between 0.0675 and 0.045 that does not exist. The branch is
    // kept because this is a property of these particular numbers; this check is
    // what fires the day a row is added that does not have it.
    for (int i = 1; i < kTableRows; ++i) {
        expectTrue("k_y interpolates to its upper station exactly without the branch",
                   kTable32[i - 1].yield + (kTable32[i].yield - kTable32[i - 1].yield) ==
                       kTable32[i].yield);
        expectTrue("k_p interpolates to its upper station exactly without the branch",
                   kTable32[i - 1].proportional +
                       (kTable32[i].proportional - kTable32[i - 1].proportional) ==
                       kTable32[i].proportional);
        expectTrue("k_E interpolates to its upper station exactly without the branch",
                   kTable32[i - 1].modulus + (kTable32[i].modulus - kTable32[i - 1].modulus) ==
                       kTable32[i].modulus);
    }

    // Monotone non-increasing over the whole range, sampled at 1 K. Steel does not
    // get stronger as it heats, and a transposed pair of rows would show here.
    double previousY = 2.0, previousP = 2.0, previousE = 2.0;
    for (int c = 20; c <= 1200; ++c) {
        const th::SteelReduction r = th::carbonSteelReduction(c + kC);
        expectTrue("k_y never rises with temperature", r.effectiveYield <= previousY);
        expectTrue("k_p never rises with temperature", r.proportionalLimit <= previousP);
        expectTrue("k_E never rises with temperature", r.youngsModulus <= previousE);
        expectTrue("the proportional limit never exceeds the effective yield",
                   r.proportionalLimit <= r.effectiveYield);
        previousY = r.effectiveYield;
        previousP = r.proportionalLimit;
        previousE = r.youngsModulus;
    }

    // Clamped outside the standard's range, the same way `carbonSteelConductivity`
    // clamps. Below 20 C the factors are 1 -- cold steel is stronger and the
    // standard says nothing, so the room-temperature value is the only defensible
    // answer and it is also the conservative one.
    for (double c : {-273.0, -100.0, 0.0, 19.99}) {
        expectTrue("below 20 C every factor is exactly 1",
                   th::carbonSteelYieldFactor(c + kC) == 1.0 &&
                       th::carbonSteelProportionalFactor(c + kC) == 1.0 &&
                       th::carbonSteelModulusFactor(c + kC) == 1.0);
    }
    for (double c : {1200.01, 1500.0, 5000.0}) {
        expectTrue("above 1200 C every factor is exactly 0",
                   th::carbonSteelYieldFactor(c + kC) == 0.0 &&
                       th::carbonSteelProportionalFactor(c + kC) == 0.0 &&
                       th::carbonSteelModulusFactor(c + kC) == 0.0);
    }

    // **Not vacuous.** The whole task is that these numbers move; a table of ones
    // would satisfy every assertion above about ordering and clamping.
    expectNear("k_y has lost 53% of the strength by 600 C",
               th::carbonSteelYieldFactor(600.0 + kC), 0.470, 0.0);
    expectTrue("and k_E has lost more than k_y has, from 500 C up",
               th::carbonSteelModulusFactor(600.0 + kC) < th::carbonSteelYieldFactor(600.0 + kC));
    expectTrue("while at 400 C the yield has not moved at all and the modulus has",
               th::carbonSteelYieldFactor(400.0 + kC) == 1.0 &&
                   th::carbonSteelModulusFactor(400.0 + kC) == 0.7);
}

// --- 2. The standard's stress-strain curve --------------------------------------

void testStressStrainCurve() {
    const double fy = 355.0e6, ea = 210.0e9;

    for (const TabulatedReduction& row : kTable32) {
        const double k = row.celsius + kC;
        if (row.modulus == 0.0) {
            expectTrue("at 1200 C the curve is identically zero",
                       th::carbonSteelStress(0.001, k, fy, ea) == 0.0 &&
                           th::carbonSteelStress(0.05, k, fy, ea) == 0.0);
            continue;
        }
        const double fyT = row.yield * fy, fpT = row.proportional * fy, eaT = row.modulus * ea;
        const double epsP = fpT / eaT;

        // The three landmarks the four branches have to meet at. Each is a value
        // the standard states, so each is asserted against arithmetic rather than
        // against the neighbouring branch.
        expectNear("the elastic branch is exactly E_a,theta times the strain",
                   th::carbonSteelStress(0.5 * epsP, k, fy, ea), 0.5 * epsP * eaT,
                   1e-9 * fyT);
        expectNear("the curve reaches the proportional limit at eps_p",
                   th::carbonSteelStress(epsP, k, fy, ea), fpT, 1e-9 * fy);
        expectNear("and the effective yield at exactly 2% strain",
                   th::carbonSteelStress(0.02, k, fy, ea), fyT, 1e-9 * fy);
        expectNear("the plateau holds it to 15%", th::carbonSteelStress(0.15, k, fy, ea), fyT,
                   1e-9 * fy);
        expectNear("and it falls linearly to half of it at 17.5%",
                   th::carbonSteelStress(0.175, k, fy, ea), 0.5 * fyT, 1e-9 * fy);
        expectTrue("and to exactly zero at 20%", th::carbonSteelStress(0.20, k, fy, ea) == 0.0);
        // **And stays there.** The standard's curve ends at eps_u; past it the
        // material is gone. Without the `eps >= eps_u` branch the falling line
        // would simply keep going and hand back a *negative* stress at 25% strain,
        // which is a mutant that survived until this line existed.
        for (double e : {0.2001, 0.25, 0.5, 2.0}) {
            expectTrue("past 20% strain there is nothing left, in either direction",
                       th::carbonSteelStress(e, k, fy, ea) == 0.0 &&
                           th::carbonSteelStress(-e, k, fy, ea) == 0.0);
        }

        // **The ellipse is tangent at both ends, and that is what makes it an
        // ellipse rather than any curve through the same two points.** At eps_p
        // its slope is E_a,theta and at 2% it is zero, so a **centred** difference
        // straddling the junction has to come back with the branch slope: a kink
        // of size `delta` there would return `E - delta/2` and be caught at first
        // order, where a one-sided difference on the smooth side would see
        // nothing. `h = 1e-9` because the ellipse's curvature is highest at its
        // steep end -- at 200 C a forward difference at 1e-6 is 6% low purely from
        // truncation, which is what a first attempt at this read as a failure.
        //
        // **At 20 C and 100 C there is no ellipse and there should be a kink.**
        // `k_p == k_y == 1` collapses the transition to nothing and the standard's
        // own curve is elastic-perfectly-plastic, corner included. Asserting
        // tangency there would be asserting against a curve the standard does not
        // draw.
        const double h = 1e-9;
        const double slopeAtY = (th::carbonSteelStress(0.02, k, fy, ea) -
                                 th::carbonSteelStress(0.02 - 1e-8, k, fy, ea)) / 1e-8;
        expectTrue("the curve arrives at 2% strain horizontally",
                   std::abs(slopeAtY) < 1e-6 * eaT);
        if (row.proportional < row.yield) {
            const double slopeAtP = (th::carbonSteelStress(epsP + h, k, fy, ea) -
                                     th::carbonSteelStress(epsP - h, k, fy, ea)) / (2.0 * h);
            expectTrue("the ellipse leaves eps_p along the elastic line",
                       std::abs(slopeAtP - eaT) < 1e-3 * eaT);
        } else {
            // The degenerate case, asserted rather than skipped: the curve is
            // exactly the effective yield everywhere past eps_p.
            expectTrue("with no ellipse the curve is flat immediately past eps_p",
                       th::carbonSteelStress(epsP + h, k, fy, ea) == fyT &&
                           th::carbonSteelStress(0.01, k, fy, ea) == fyT);
        }

        // Monotone up to 2% and odd in the strain, both over a fine sweep.
        double previous = 0.0;
        for (int i = 0; i <= 4000; ++i) {
            const double e = 0.02 * i / 4000.0;
            const double s = th::carbonSteelStress(e, k, fy, ea);
            expectTrue("the curve never falls before 2% strain", s >= previous - 1e-6);
            expectTrue("the curve is odd in the strain",
                       th::carbonSteelStress(-e, k, fy, ea) == -s);
            expectTrue("and never exceeds the effective yield", s <= fyT * (1.0 + 1e-12));
            previous = s;
        }

        // **The plateau, swept rather than sampled at three points.** Removing the
        // `eps >= eps_y` branch does not break the plateau everywhere -- past
        // about 3.9% strain the ellipse's argument goes negative and its guard
        // returns `f_y` anyway, so 5%, 10% and 15% all still look right. What sags
        // is the narrow band just past 2%, by 3e-4 of `f_y` at 2.1%, and three
        // spot checks walked straight past it. This sweep is what caught it.
        for (int i = 0; i <= 1300; ++i) {
            const double e = 0.02 + (0.15 - 0.02) * i / 1300.0;
            expectTrue("the plateau holds the effective yield across the whole of it",
                       std::abs(th::carbonSteelStress(e, k, fy, ea) - fyT) <= 1e-9 * fy);
        }
    }

    // **What the J2 model this file couples to gives up, measured.**
    // `plasticity::Material` carries one yield strength; scaled by `k_y` it is an
    // elastic-perfectly-plastic curve at `k_E E` and `k_y f_y`. Where it differs
    // from the standard is between `k_p f_y` and `k_y f_y`, and the difference is
    // bounded by `(k_y - k_p) f_y` by construction. The measurement is the point:
    // the worst temperature is **400 C**, not the hottest one.
    double worstGap = 0, worstAt = 0, worstTemperature = 0;
    for (const TabulatedReduction& row : kTable32) {
        if (row.modulus == 0.0) continue;
        const double k = row.celsius + kC;
        double gap = 0, at = 0;
        for (int i = 0; i <= 4000; ++i) {
            const double e = 0.02 * i / 4000.0;
            const double j2 = std::min(row.modulus * ea * e, row.yield * fy);
            const double d = std::abs(j2 - th::carbonSteelStress(e, k, fy, ea));
            if (d > gap) { gap = d; at = e; }
        }
        expectTrue("the J2 gap never exceeds (k_y - k_p) f_y",
                   gap <= (row.yield - row.proportional) * fy + 1e-6);
        if (gap > worstGap) { worstGap = gap; worstAt = at; worstTemperature = row.celsius; }
        // Past 2% strain the two are the same curve, identically.
        for (double e : {0.02, 0.05, 0.10, 0.15}) {
            expectNear("past 2% strain the J2 curve and the standard's agree",
                       th::carbonSteelStress(e, k, fy, ea), row.yield * fy, 1e-9 * fy);
        }
    }
    expectNear("the J2 model is exact at 20 C, where k_p == k_y",
               th::carbonSteelStress(0.001, 20.0 + kC, fy, ea), 0.001 * ea, 1e-9 * fy);
    expectNear("the worst J2 gap is 0.388 f_y", worstGap / fy, 0.388, 0.002);
    expectNear("and it is at 400 C, where the yield has not moved and the limit has",
               worstTemperature, 400.0, 0.0);
    std::printf("     EN 1993-1-2 §3.2.1 vs a k_y-scaled J2 curve: worst gap %.3f f_y at %.0f C,"
                " %.3f%% strain; identical at and past 2%%\n",
                worstGap / fy, worstTemperature, 100.0 * worstAt);
}

// --- 3. The coupling, and the exact control at 20 C -------------------------------

void testMaterialsAtTemperature() {
    const StructuralMaterial cold = ah36Steel();

    // **The exact control.** At 20 C -- and anywhere below it, where the factors
    // clamp to 1 -- a reduced material is the same material, bit for bit, on every
    // field. Not "close": `x *= 1.0` is exact for every finite double, and a
    // reduction that came back merely nearly equal would put a difference into
    // every downstream result on an unheated ship.
    for (double c : {20.0, 0.0, -50.0}) {
        const StructuralMaterial same = th::atTemperature(cold, c + kC);
        expectTrue("an unheated StructuralMaterial is bit-identical",
                   same.youngsModulus == cold.youngsModulus &&
                       same.yieldStrength == cold.yieldStrength &&
                       same.poissonRatio == cold.poissonRatio && same.density == cold.density &&
                       same.conductivity == cold.conductivity &&
                       same.specificHeat == cold.specificHeat && same.name == cold.name);
    }

    // The two that move, and the four that do not.
    //
    // **Against the published accessor, not against the table literal.** The claim
    // being made is that `atTemperature` applies exactly the factor
    // `carbonSteelModulusFactor` reports -- one multiply and no rounding of its
    // own. Comparing against `row.modulus` instead would fold in the Kelvin
    // round-trip error `testReductionTable` already isolates, and would make this
    // check fail for a reason that has nothing to do with the coupling.
    for (const TabulatedReduction& row : kTable32) {
        const double k = row.celsius + kC;
        const StructuralMaterial hot = th::atTemperature(cold, k);
        expectTrue("E is scaled by exactly k_E",
                   hot.youngsModulus == cold.youngsModulus * th::carbonSteelModulusFactor(k));
        expectTrue("the yield strength is scaled by exactly k_y",
                   hot.yieldStrength == cold.yieldStrength * th::carbonSteelYieldFactor(k));
        expectTrue("density, Poisson ratio and the 20 C thermal properties are untouched",
                   hot.density == cold.density && hot.poissonRatio == cold.poissonRatio &&
                       hot.conductivity == cold.conductivity &&
                       hot.specificHeat == cold.specificHeat);
    }

    // --- the flow model ---
    //
    // Two curves, because they scale through different fields: Linear carries a
    // hardening modulus, Swift a strength coefficient, and a reduction that
    // touched only `yieldStrength` would pass every Swift test while leaving a
    // Linear curve hardening back to its cold strength.
    sim::plasticity::Material swift = sim::plasticity::shipSteel();
    sim::plasticity::Material linear = swift;
    linear.flow = sim::plasticity::linearHardening(355.0e6, 1.2e9);
    // **Kinematic hardening is off by default and that made its scaling
    // untestable.** A mutant that deleted `kinematicModulus *= k_y` survived the
    // whole suite because every material in it carried zero there. `plasticity.hpp`
    // keeps the path live and tested precisely so switching it on is a parameter
    // rather than a rewrite, so the reduction has to reach it too.
    linear.flow.kinematicModulus = 8.0e8;
    swift.flow.kinematicModulus = 6.0e8;

    for (const sim::plasticity::Material& base : {swift, linear}) {
        for (double c : {20.0, 0.0}) {
            const sim::plasticity::Material same = th::atTemperature(base, c + kC);
            expectTrue("an unheated plasticity::Material is bit-identical",
                       same.youngsModulus == base.youngsModulus &&
                           same.poissonRatio == base.poissonRatio &&
                           same.flow.yieldStrength == base.flow.yieldStrength &&
                           same.flow.hardeningModulus == base.flow.hardeningModulus &&
                           same.flow.strengthCoefficient == base.flow.strengthCoefficient &&
                           same.flow.referenceStrain == base.flow.referenceStrain &&
                           same.flow.hardeningExponent == base.flow.hardeningExponent &&
                           same.flow.kinematicModulus == base.flow.kinematicModulus &&
                           same.failure.uniformStrain == base.failure.uniformStrain &&
                           same.failure.fractureStrain == base.failure.fractureStrain);
        }

        for (const TabulatedReduction& row : kTable32) {
            const double k = row.celsius + kC;
            const sim::plasticity::Material hot = th::atTemperature(base, k);
            expectTrue("E is scaled by exactly k_E",
                       hot.youngsModulus == base.youngsModulus * th::carbonSteelModulusFactor(k));

            // **The whole curve scales, at every plastic strain.** Sampling only
            // at eps_p = 0 would be satisfied by scaling the intercept alone.
            for (double p : {0.0, 0.001, 0.01, 0.05, 0.2, 0.5}) {
                const double coldStress = sim::plasticity::flowStress(base.flow, p);
                expectNear("the flow stress scales by k_y at every plastic strain",
                           sim::plasticity::flowStress(hot.flow, p), row.yield * coldStress,
                           1e-9 * coldStress + 1e-9);
                const double coldSlope = sim::plasticity::flowSlope(base.flow, p);
                expectNear("and so does the hardening slope",
                           sim::plasticity::flowSlope(hot.flow, p), row.yield * coldSlope,
                           1e-9 * std::abs(coldSlope) + 1e-9);
            }

            // **Considere's necking strain does not move.** `dsigma/deps = sigma`
            // has both sides multiplied by the same positive number, so its root
            // cannot move -- and that is why `Failure` is left alone rather than
            // adjusted. At 1200 C the curve has no strength at all and the root is
            // undefined, which is the one case excluded.
            //
            // For **Swift** it is bit-exact, and structurally so: `n - eps_0` is
            // built from two fields `atTemperature` does not touch. For **Linear**
            // it is `1 - sigma_y0/H`, and `(k sigma_y0)/(k H)` is not `sigma_y0/H`
            // for every double -- the invariance is exact in the reals and one
            // unit in the last place in floating point. Claiming bit-exactness for
            // both was this test's first mistake; the two are asserted apart
            // rather than to the looser of the two, which would have hidden a
            // genuine drift on the Swift path.
            if (row.yield > 0.0) {
                const double got = sim::plasticity::uniformElongation(hot.flow);
                const double want = sim::plasticity::uniformElongation(base.flow);
                if (base.flow.kind == sim::plasticity::Hardening::Swift)
                    expectTrue("a scaled Swift curve necks at exactly the same strain",
                               got == want);
                else
                    expectNear("a scaled linear curve necks at the same strain to rounding",
                               got, want, 4e-16 * std::abs(want));
            }
            expectTrue("and the failure constants are untouched",
                       hot.failure.uniformStrain == base.failure.uniformStrain &&
                           hot.failure.fractureStrain == base.failure.fractureStrain &&
                           hot.failure.triaxialitySensitivity ==
                               base.failure.triaxialitySensitivity);

            // Every field of the flow curve that carries a stress, including the
            // kinematic modulus a default material leaves at zero.
            expectTrue("the kinematic modulus is scaled by exactly k_y too",
                       hot.flow.kinematicModulus ==
                           base.flow.kinematicModulus * th::carbonSteelYieldFactor(k));
            expectTrue("the two strain-valued Swift parameters are not scaled",
                       hot.flow.referenceStrain == base.flow.referenceStrain &&
                           hot.flow.hardeningExponent == base.flow.hardeningExponent);
            expectTrue("and the base curve really has a kinematic modulus to scale",
                       base.flow.kinematicModulus > 0.0);
        }
    }

    // A dead material at 1200 C has to come back finite rather than NaN: the
    // return map divides by the shear modulus, and the standard's table ends at
    // zero for all three factors.
    const sim::plasticity::Material dead = th::atTemperature(swift, 1200.0 + kC);
    double strain[6] = {0.01, 0.0, 0.0, 0.0, 0.0, 0.0}, stress[6] = {};
    sim::plasticity::State state;
    sim::plasticity::update(dead, sim::plasticity::kNeverFails, strain, state, stress);
    expectTrue("a material with no strength and no stiffness returns zero, not NaN",
               stress[0] == 0.0 && stress[1] == 0.0 && stress[2] == 0.0 &&
                   std::isfinite(state.equivalentPlasticStrain));
}

// --- 4. Does the return map need a term it does not have? ------------------------

void testReturnMapUnderAMovingYieldSurface() {
    using namespace sim::plasticity;
    const Material cold = shipSteel();
    const double e = 0.006;  // deviatoric, well past yield at every temperature here

    // **The map's single strongest property survives the reduction.** At a fixed
    // temperature the reduced curve is still monotonically hardening -- `k_y` is a
    // positive constant multiplying both `sigma_y` and its slope -- so the scalar
    // consistency equation still has exactly one root and radial return still
    // lands on it whatever the step count. That is `testStepIndependence` from
    // `test_plasticity.cpp` re-run at 600 C, and it is what says the return map
    // needs no new term: the softening is *between* steps, not inside one.
    double reference[6] = {};
    double referenceStrain = 0;
    for (int n : {1, 7, 100}) {
        State s;
        double stress[6] = {};
        for (int i = 1; i <= n; ++i) {
            const double f = 0.01 * i / n;
            double strain[6] = {f, -0.5 * f, -0.5 * f, 0, 0, 0};
            update(th::atTemperature(cold, 600.0 + kC), kNeverFails, strain, s, stress);
        }
        if (n == 1) {
            for (int i = 0; i < 6; ++i) reference[i] = stress[i];
            referenceStrain = s.equivalentPlasticStrain;
        } else {
            expectNear("the return map at 600 C is step-independent to rounding",
                       vonMises(stress), vonMises(reference), 1e-9 * vonMises(reference));
            expectNear("and so is the plastic strain it accumulated",
                       s.equivalentPlasticStrain, referenceStrain, 1e-12 * referenceStrain);
        }
    }
    expectTrue("and it actually yielded, so there is something to be independent of",
               referenceStrain > 1e-3);

    // **The shrinking surface.** 400 C to 500 C at constant total strain: `k_y`
    // falls from 1.00 to 0.78 while `k_E` falls only from 0.70 to 0.60, so the
    // elastic demand falls *slower* than the capacity and the point has to flow
    // again. This is the case a hardening-only return map is supposed not to
    // handle, and it needs nothing: the stored stress simply starts outside the
    // new surface and the map returns it, which is stress relaxation and is what a
    // heated restrained member physically does.
    {
        const Material at400 = th::atTemperature(cold, 400.0 + kC);
        const Material at500 = th::atTemperature(cold, 500.0 + kC);
        double strain[6] = {e, -0.5 * e, -0.5 * e, 0, 0, 0}, stress[6] = {};
        State s;
        update(at400, kNeverFails, strain, s, stress);
        const double before = s.equivalentPlasticStrain;
        expectNear("at 400 C the point sits exactly on its yield surface", vonMises(stress),
                   flowStress(at400.flow, before), 1e-9 * flowStress(at400.flow, before));

        const Increment inc = update(at500, kNeverFails, strain, s, stress);
        expectTrue("heating at constant strain makes it flow again", inc.yielded);
        expectTrue("and the plastic strain grows by a real amount, not a rounding",
                   s.equivalentPlasticStrain > before * 1.02);
        expectNear("landing exactly on the smaller surface", vonMises(stress),
                   flowStress(at500.flow, s.equivalentPlasticStrain),
                   1e-9 * flowStress(at500.flow, s.equivalentPlasticStrain));
        expectTrue("and the increment dissipated energy rather than absorbing it",
                   inc.dissipation >= 0.0);
    }

    // **The other direction, and it is the one worth naming.** 200 C to 400 C:
    // `k_y` is 1 at both ends while `k_E` falls 0.90 to 0.70, so the elastic
    // demand falls *faster* than the capacity and the point unloads. The plastic
    // strain must then be **bit-identical** -- not nearly -- and the stress must be
    // exactly the elastic one at the frozen plastic strain.
    {
        const Material at200 = th::atTemperature(cold, 200.0 + kC);
        const Material at400 = th::atTemperature(cold, 400.0 + kC);
        double strain[6] = {e, -0.5 * e, -0.5 * e, 0, 0, 0}, stress[6] = {};
        State s;
        update(at200, kNeverFails, strain, s, stress);
        const double frozen = s.equivalentPlasticStrain;
        expectTrue("it yielded at 200 C first", frozen > 1e-4);

        const Increment inc = update(at400, kNeverFails, strain, s, stress);
        expectTrue("heating past it unloads elastically instead", !inc.yielded);
        expectTrue("so the plastic strain is bit-identical",
                   s.equivalentPlasticStrain == frozen);
        // Uniaxial deviatoric: the von Mises stress of an elastic point at this
        // strain is 3 G (e - e_p), and e_p here is the axial component of the
        // stored plastic strain, which for this path is the equivalent one.
        expectNear("and the stress is exactly the elastic stress of the hot material",
                   vonMises(stress), 3.0 * at400.shearModulus() * (e - frozen),
                   1e-9 * vonMises(stress));
    }

    // **What this means for a caller, stated as a measurement rather than a
    // caveat.** Because `k_E` falls faster than `k_y` from 100 C to 400 C, a point
    // held at constant strain accumulates *no* further plastic strain over that
    // range -- the whole of its plastic history is set by the coldest, stiffest
    // state it passed through. That is right for this model and wrong for real
    // steel, and the missing term is thermal elongation, which is not here. The
    // number is in `thermal.hpp`: 7.08e-3 of free expansion strain over a 500 K
    // rise against a 1.72e-3 yield strain, a factor of 4.1, so for restrained
    // structure the missing term arrives first and is four times the size.
    {
        const Material at100 = th::atTemperature(cold, 100.0 + kC);
        double strain[6] = {e, -0.5 * e, -0.5 * e, 0, 0, 0}, stress[6] = {};
        State s;
        update(at100, kNeverFails, strain, s, stress);
        const double at = s.equivalentPlasticStrain;
        for (double c : {200.0, 300.0, 400.0})
            update(th::atTemperature(cold, c + kC), kNeverFails, strain, s, stress);
        expectTrue("no plastic strain is accumulated from 100 C to 400 C at constant strain",
                   s.equivalentPlasticStrain == at);
        expectTrue("and there was some to begin with", at > 1e-4);
    }
}

// --- 5. Does the same material path reach a stiffener fibre? ----------------------

void testFibresSoften() {
    using namespace sim::constraint;
    using sim::plasticity::Material;

    // The smallest fibre set that is still a fibre set: one bar between two node
    // pairs, tied at the mid-surface. `fiberForces` reads its stiffness, its flow
    // curve and its failure constants from the `plasticity::Material` argument and
    // from nothing else -- `Stiffening::material` is used for mass -- so the
    // reduced material reaches the fibres through exactly the same door as the
    // solid elements, with no interface change at all. That is what this asserts.
    Stiffening stiffening;
    stiffening.material = ah36Steel();
    Fiber fiber;
    fiber.area = 2.0e-3;
    fiber.neckWidth = 0.010;
    fiber.end[0].bottom = 0; fiber.end[0].top = 1; fiber.end[0].weight = 0.5;
    fiber.end[1].bottom = 2; fiber.end[1].top = 3; fiber.end[1].weight = 0.5;
    stiffening.fiber.push_back(fiber);

    const std::vector<double> rest = {0, 0, 0,  0, 0, 0.012,  1, 0, 0,  1, 0, 0.012};
    const RestFibers forms = restFibers(stiffening, rest);
    expectTrue("the fibre set formed", forms.ok && forms.length.size() == 1);
    expectNear("with the rest length it was built at", forms.length[0], 1.0, 1e-12);

    const double stretch = 0.004;  // 0.4% strain: past yield at every temperature here
    std::vector<double> current = rest;
    current[6] += stretch;
    current[9] += stretch;

    double coldForce = 0, hotForce = 0;
    for (int hot = 0; hot < 2; ++hot) {
        const Material material =
            th::atTemperature(sim::plasticity::shipSteel(), (hot ? 600.0 : 20.0) + kC);
        std::vector<FiberState> state(1);
        std::vector<double> force(rest.size(), 0.0);
        fiberForces(stiffening, forms, current, material, &state, force);
        // The axial force appears on the far pair, split by the tie weight.
        const double axial = std::abs(force[6]) + std::abs(force[9]);
        (hot ? hotForce : coldForce) = axial;
        expectTrue("the fibre yielded", state[0].equivalentPlasticStrain > 0.0);
        // On the yield surface: N = sigma_y(eps_p) A, in closed form.
        const double want =
            sim::plasticity::flowStress(material.flow, state[0].equivalentPlasticStrain) *
            fiber.area;
        expectNear("a yielded fibre carries sigma_y(eps_p) times its area", axial, want,
                   1e-9 * want);
    }

    // **The consequence, and the guard against vacuity.** The same stretch through
    // the same fibre carries roughly half the force at 600 C. If `atTemperature`
    // ignored the reduction factors the two would be identical, which this cannot
    // pass.
    expectTrue("a hot fibre carries substantially less", hotForce < 0.6 * coldForce);
    expectTrue("and it is not zero either", hotForce > 0.2 * coldForce);
    std::printf("     one stiffener fibre at 0.4%% strain: %.1f kN cold, %.1f kN at 600 C"
                " (%.0f%%)\n",
                coldForce * 1e-3, hotForce * 1e-3, 100.0 * hotForce / coldForce);

    // And an unheated fibre is bit-identical to one solved with the unreduced
    // material -- the exact 20 C control, on this path too.
    {
        const Material base = sim::plasticity::shipSteel();
        std::vector<double> a(rest.size(), 0.0), b(rest.size(), 0.0);
        std::vector<FiberState> sa(1), sb(1);
        fiberForces(stiffening, forms, current, base, &sa, a);
        fiberForces(stiffening, forms, current, th::atTemperature(base, 20.0 + kC), &sb, b);
        bool identical = true;
        for (std::size_t i = 0; i < a.size(); ++i) identical = identical && a[i] == b[i];
        expectTrue("an unheated fibre is bit-identical to the unreduced one",
                   identical && sa[0].equivalentPlasticStrain == sb[0].equivalentPlasticStrain);
    }
}

// --- 6. The structural consequence, with hand-computed ends ----------------------

void testHotStructureIsWeaker() {
    // (a) **A plate strip's collapse pressure**, which is a closed form: the
    // three-hinge mechanism of a strip clamped at both ends, `4 f_y (t/s)^2`.
    // Linear in the yield strength and in nothing else, so it falls by exactly
    // `k_y` -- the one place in this test where the yield factor is the whole
    // answer, and it is here to be contrasted with (b) and (c) where it is not.
    {
        const double t = 0.012, s = 0.70;
        const double cold = sim::promotion::platingCollapsePressure(
            th::atTemperature(ah36Steel(), 20.0 + kC).yieldStrength, t, s);
        const double hot = sim::promotion::platingCollapsePressure(
            th::atTemperature(ah36Steel(), 600.0 + kC).yieldStrength, t, s);
        // By hand: 4 * 355e6 * (0.012/0.70)^2 = 417 306.12... Pa.
        expectNear("the cold plate collapses at 4 f_y (t/s)^2", cold, 417306.122448979, 1e-6);
        expectNear("and the hot one at exactly 0.47 of that", hot, 0.47 * 417306.122448979,
                   1e-6);
    }

    // (b) **The compressive capacity of the plating between two stiffeners**,
    // which is *not* linear in the yield strength, and this is the finding. The
    // elastic buckling stress
    // goes with `E` and the Johnson-Ostenfeld cap with `f_y`, and above 500 C
    // `k_E < k_y`. On this panel it changes the regime outright: cold, it is
    // squash-governed and the correction bites; at 600 C the elastic stress has
    // fallen below half the hot yield and the panel is a pure elastic buckling
    // problem again.
    {
        const StructuralMaterial cold = th::atTemperature(ah36Steel(), 20.0 + kC);
        const StructuralMaterial hot = th::atTemperature(ah36Steel(), 600.0 + kC);
        const sim::BucklingCheck c = sim::plateBuckling(0.012, 2.4, 0.70, 0.0, cold);
        const sim::BucklingCheck h = sim::plateBuckling(0.012, 2.4, 0.70, 0.0, hot);

        // Hand-computed. alpha = 2.4/0.7 = 3.4286, so m = 3 half-waves and
        // k = (3/alpha + alpha/3)^2 = 4.071747...
        //   sigma_e = k pi^2 E / (12 (1 - nu^2)) (t/b)^2
        //           = 4.071747 * 9.8696044 * 206e9 / 10.92 * (0.012/0.70)^2
        //           = 2.22788e8 Pa
        // and sigma_e > f_y/2 = 1.775e8, so Johnson-Ostenfeld applies:
        //   sigma_c = 355e6 (1 - 355e6 / (4 * 2.22788e8)) = 2.13582e8 Pa.
        expectNear("the cold panel's buckling coefficient is the m = 3 minimum", c.coefficient,
                   4.0717474490, 1e-8);
        expectNear("its elastic buckling stress is the hand value", c.elasticStress, 2.227879e8,
                   1e3);
        expectNear("and Johnson-Ostenfeld caps it at the hand value", c.criticalStress,
                   2.135819e8, 1e3);
        expectTrue("because the cold panel is squash-governed",
                   c.elasticStress > 0.5 * cold.yieldStrength);

        // At 600 C: sigma_e scales by exactly k_E = 0.31 to 6.90642e7, and
        // f_y/2 = 0.5 * 0.47 * 355e6 = 8.3425e7, which is now *above* it -- so the
        // correction switches off and the capacity is the raw eigenvalue.
        expectNear("the hot elastic buckling stress is exactly 0.31 of the cold one",
                   h.elasticStress, 0.31 * c.elasticStress, 1e-9 * c.elasticStress);
        expectTrue("and the hot panel is elastic-buckling-governed instead",
                   h.elasticStress <= 0.5 * hot.yieldStrength);
        expectTrue("so no correction is applied at all",
                   h.criticalStress == h.elasticStress);
        expectNear("the hot capacity is the hand value", h.criticalStress, 6.906424e7, 1e2);

        // **The point.** The capacity falls to 0.323, well below k_y = 0.47: a
        // model that reduced only the yield strength would over-predict this panel
        // by 45%.
        const double ratio = h.criticalStress / c.criticalStress;
        expectNear("the hot panel keeps 0.323 of its capacity", ratio, 0.3233, 5e-4);
        expectTrue("which is well below the yield reduction factor",
                   ratio < 0.85 * th::carbonSteelYieldFactor(600.0 + kC));
        std::printf("     0.7 x 2.4 m x 12 mm panel: capacity %.1f MPa cold, %.1f MPa at 600 C"
                    " (%.3f, against k_y = %.2f) -- squash-governed becomes"
                    " buckling-governed\n",
                    c.criticalStress * 1e-6, h.criticalStress * 1e-6, ratio,
                    th::carbonSteelYieldFactor(600.0 + kC));
    }

    // (c) **The ferry's midship section**, all the way through Smith's method. A
    // hot ship is the same ship with a reduced material: nothing in `collapse.cpp`
    // or `scantlings.cpp` changes, because every element already looks its
    // material up by index in `StructuralMesh::materials`.
    {
        const sim::StructuralMesh mesh =
            sim::makeStructuralMesh(game::buildFerry().hull, sim::ferryScantlings());
        const sim::Scantlings scantlings = sim::ferryScantlings();
        expectTrue("the reference section was built", !mesh.panels.empty());

        const auto at = [&](double c) {
            sim::StructuralMesh hot = mesh;
            for (StructuralMaterial& m : hot.materials) m = th::atTemperature(m, c + kC);
            return sim::collapseElementsAt(hot, scantlings, 0.0);
        };

        const std::vector<sim::CollapseElement> cold = at(20.0);
        const std::vector<sim::CollapseElement> warm = at(400.0);
        const std::vector<sim::CollapseElement> hot = at(600.0);

        // **The exact control**, on the whole chain rather than on one function:
        // the section built through `atTemperature` at 20 C is bit-identical to the
        // one built without it.
        const std::vector<sim::CollapseElement> control =
            sim::collapseElementsAt(mesh, scantlings, 0.0);
        expectEqual("the section has the same elements at 20 C",
                    static_cast<long long>(cold.size()),
                    static_cast<long long>(control.size()));
        bool identical = !cold.empty();
        for (std::size_t i = 0; i < cold.size() && i < control.size(); ++i)
            identical = identical && cold[i].area == control[i].area &&
                        cold[i].height == control[i].height &&
                        cold[i].curve.yieldStrength == control[i].curve.yieldStrength &&
                        cold[i].curve.youngsModulus == control[i].curve.youngsModulus &&
                        cold[i].curve.bucklingStress == control[i].curve.bucklingStress;
        expectTrue("and every one of them is bit-identical to the unheated section",
                   identical);

        // The fully plastic moment is linear in every element's yield strength, so
        // it scales by exactly `k_y`. That is the one closed form Smith's method
        // has, and it is the anchor.
        const double mpCold = sim::fullyPlasticMoment(cold);
        expectTrue("the section has a fully plastic moment", mpCold > 1e9);
        expectNear("the fully plastic moment scales by exactly k_y at 400 C",
                   sim::fullyPlasticMoment(warm), mpCold, 1e-9 * mpCold);
        expectNear("and by exactly k_y at 600 C", sim::fullyPlasticMoment(hot), 0.47 * mpCold,
                   1e-9 * mpCold);

        // The ultimate moment does not, and that is the consequence worth having.
        const double uCold = std::abs(sim::collapseCurve(cold, -1.0).ultimateMoment);
        const double uWarm = std::abs(sim::collapseCurve(warm, -1.0).ultimateMoment);
        const double uHot = std::abs(sim::collapseCurve(hot, -1.0).ultimateMoment);
        expectTrue("the cold section collapses well below its fully plastic moment",
                   uCold > 0.3 * mpCold && uCold < 0.8 * mpCold);

        // **The test that fails if the modulus factor is ignored.** At 400 C
        // `k_y` is exactly 1, so a yield-only reduction leaves the section
        // completely unchanged -- and it is not: it has lost 17.6% of its ultimate
        // strength through `k_E` alone, because the plate panels' buckling stress
        // is proportional to E.
        expectTrue("at 400 C the yield strength has not moved at all",
                   th::carbonSteelYieldFactor(400.0 + kC) == 1.0);
        expectNear("yet the section has lost 17.6% of its ultimate strength", uWarm / uCold,
                   0.824, 0.01);

        // And at 600 C it loses more than `k_y` says, for the same reason.
        expectNear("at 600 C it keeps 0.376 of its ultimate strength", uHot / uCold, 0.376,
                   0.01);
        expectTrue("which is below the yield reduction factor of 0.47",
                   uHot / uCold < 0.9 * th::carbonSteelYieldFactor(600.0 + kC));
        std::printf("     ferry midship, sagging: Mu = %.3f GN m cold, %.3f at 400 C (%.3f),"
                    " %.3f at 600 C (%.3f vs k_y = 0.47); Mp scales by k_y exactly\n",
                    uCold * 1e-9, uWarm * 1e-9, uWarm / uCold, uHot * 1e-9, uHot / uCold);
    }
}

// --- 7. From a nodal field to an element's temperature ---------------------------

void testElementTemperatureField() {
    // A box, so the volume average of a linear field has a closed form: the field
    // at the box centre. On a general trilinear hex the interpolant of a linear
    // field is not that field, which is why the closed form is claimed here and
    // not there.
    const double h = 0.4, w = 0.3, d = 0.012;
    double nodes[ss::kDof];
    double nodal[ss::kNodes];
    const double x[8] = {0, h, h, 0, 0, h, h, 0};
    const double y[8] = {0, 0, w, w, 0, 0, w, w};
    const double z[8] = {0, 0, 0, 0, d, d, d, d};
    for (int a = 0; a < ss::kNodes; ++a) {
        nodes[a * 3] = x[a];
        nodes[a * 3 + 1] = y[a];
        nodes[a * 3 + 2] = z[a];
        // An oblique linear field, so a routine that got one axis right and two
        // wrong does not pass.
        nodal[a] = kC + 300.0 + 500.0 * x[a] / h + 200.0 * y[a] / w + 700.0 * z[a] / d;
    }
    th::Forms forms;
    expectTrue("the box forms", th::computeForms(nodes, forms));

    expectNear("the element mean of a linear field is the field at the centre",
               th::elementTemperature(forms, nodal), kC + 300.0 + 250.0 + 100.0 + 350.0, 1e-9);

    double gauss[ss::kGauss];
    th::gaussTemperature(nodal, gauss);
    const double q = 1.0 / std::sqrt(3.0);
    for (int gp = 0; gp < ss::kGauss; ++gp) {
        // **The Gauss points are at +-1/sqrt(3) in each direction, in the bit
        // order `solid_shell.cpp` uses.** Written as a closed form in the field's
        // own values rather than by reading either file's constants, so a routine
        // that sampled the corners, the centre, or the right points in the wrong
        // order fails.
        const double want = kC + 300.0 + 250.0 * (1.0 + ((gp & 1) ? q : -q)) +
                            100.0 * (1.0 + ((gp & 2) ? q : -q)) +
                            350.0 * (1.0 + ((gp & 4) ? q : -q));
        expectNear("a nodal field interpolates to the 2x2x2 Gauss points", gauss[gp], want,
                   1e-9);
        // And it is exactly what the conduction element's own shape functions say,
        // which is the tie to the operator this field will be read beside.
        double byShape = 0.0;
        for (int a = 0; a < ss::kNodes; ++a) byShape += forms.shape[gp][a] * nodal[a];
        expectNear("and agrees with Forms::shape to rounding", gauss[gp], byShape, 1e-9);
    }

    // The mesh-level bridge, on a plate with four elements through the thickness:
    // a field linear in z gives each element the temperature of its own mid-plane,
    // in closed form, and the four are distinct.
    //
    // **The expectation is derived from the mesh's own coordinates**, because
    // `makePlateMesh` centres the plate on `z = 0` rather than resting it on it --
    // this test first asserted the wrong four numbers, all of them 500 K out, for
    // exactly that reason. Reading the extent back means the closed form cannot be
    // wrong about where the plate is while still being a closed form about what
    // the average of a linear field over an element is.
    ss::HexMesh mesh = ss::makePlateMesh(0.2, 0.2, 0.012, 1, 1, 4);
    double zLo = 1e30, zHi = -1e30;
    for (std::size_t n = 0; n < mesh.nodeCount(); ++n) {
        zLo = std::min(zLo, mesh.position[n * 3 + 2]);
        zHi = std::max(zHi, mesh.position[n * 3 + 2]);
    }
    expectNear("the plate is 12 mm thick however it is placed", zHi - zLo, 0.012, 1e-15);
    std::vector<double> field(mesh.nodeCount());
    for (std::size_t n = 0; n < mesh.nodeCount(); ++n)
        field[n] = kC + 20.0 + 1000.0 * (mesh.position[n * 3 + 2] - zLo) / (zHi - zLo);
    std::vector<double> perElement;
    expectTrue("the mesh-level bridge works", th::elementTemperatures(mesh, field, perElement));
    expectEqual("one temperature per element", static_cast<long long>(perElement.size()),
                static_cast<long long>(mesh.elementCount()));
    std::vector<double> sorted = perElement;
    std::sort(sorted.begin(), sorted.end());
    for (int i = 0; i < 4; ++i)
        expectNear("each element sees its own layer's mid-plane temperature", sorted[i],
                   kC + 20.0 + 1000.0 * (i + 0.5) / 4.0, 1e-9);
    expectTrue("and the four layers are genuinely different temperatures",
               sorted[3] - sorted[0] > 700.0);
    expectTrue("a field of the wrong length is refused",
               !th::elementTemperatures(mesh, std::vector<double>(3, kC), perElement) &&
                   perElement.empty());

    // **The refusal paths, which nothing reached until mutation testing said so.**
    // Three mutants survived here: deleting `if (!forms.ok) return 0`, deleting the
    // `ok = false` on a bad element, and deleting the `continue` that skips it. All
    // three are on the same road and nothing in this file had driven down it.
    {
        // **A `Forms` with good weights and `ok` false**, not a default-constructed
        // one. Only `ok` has an initialiser, so a default `Forms` has indeterminate
        // weights -- which in practice come out zero, so `volume > 0` is false and
        // the answer is zero whether or not the `ok` guard is there. That made
        // deleting the guard an accidental survivor. Copying a good set and
        // clearing the flag is what makes the flag the only thing being tested.
        th::Forms broken = forms;
        broken.ok = false;
        expectTrue("an element whose forms did not build has no temperature",
                   th::elementTemperature(broken, nodal) == 0.0);
        expectTrue("and a good one does, so the check is not passing by accident",
                   th::elementTemperature(forms, nodal) != 0.0);
    }
    {
        // Invert one element by turning its through-thickness direction inside
        // out: `computeForms` refuses a non-positive Jacobian, which is the same
        // condition `solidshell::computeRestForms` refuses on.
        ss::HexMesh bad = ss::makePlateMesh(0.2, 0.2, 0.012, 2, 1, 1);
        for (int a = 0; a < 4; ++a)
            std::swap(bad.index[static_cast<std::size_t>(a)],
                      bad.index[static_cast<std::size_t>(a) + 4]);
        std::vector<double> hot(bad.nodeCount());
        for (std::size_t n = 0; n < bad.nodeCount(); ++n)
            hot[n] = kC + 500.0 + 100.0 * bad.position[n * 3];
        std::vector<double> got;
        expectTrue("a mesh with an inverted element is reported, not silently averaged",
                   !th::elementTemperatures(bad, hot, got));
        expectEqual("but the good elements are still filled in",
                    static_cast<long long>(got.size()),
                    static_cast<long long>(bad.elementCount()));
        expectTrue("the inverted element is left at the fill value and the sound one is not",
                   got[0] == 0.0 && got[1] > kC + 100.0);
    }
}

// --- 8. What the gradient actually is, and what per-element averaging costs -------

void testGradientThroughPlating() {
    // A 12 mm plate with a post-flashover compartment on one face and nothing on
    // the other. This is the measurement that decides whether a per-element
    // temperature is enough, and it is run rather than argued.
    ss::HexMesh mesh = ss::makePlateMesh(0.10, 0.10, 0.012, 1, 1, 4);
    th::Problem problem;
    problem.mesh = &mesh;
    problem.material = ah36Steel();
    problem.temperatureDependent = true;

    th::Film fire;
    fire.coefficient = 200.0;  // EN's 25 W/(m^2 K) convective plus ~175 radiative
    fire.ambient = kC + 900.0;
    for (const th::BoundaryFace& f : th::boundaryFaces(mesh))
        if (f.normal.z < -0.9) fire.face.push_back(f);
    expectTrue("the fire lands on the one face it should", fire.face.size() == 1);
    problem.film.push_back(fire);

    th::Solver solver;
    std::string why;
    expectTrue("the fire problem prepares", solver.prepare(problem, kC + 20.0, &why));

    double worstSpread = 0, worstFactorSpread = 0, hottest = 0;
    for (int step = 0; step < 1800; ++step) {
        expectTrue("the fire steps", solver.step(1.0, &why));
        double lo = 1e30, hi = -1e30;
        for (double v : solver.temperature()) {
            lo = std::min(lo, v);
            hi = std::max(hi, v);
        }
        worstSpread = std::max(worstSpread, hi - lo);
        worstFactorSpread = std::max(worstFactorSpread, std::abs(th::carbonSteelYieldFactor(lo) -
                                                                 th::carbonSteelYieldFactor(hi)));
        hottest = std::max(hottest, hi);
    }
    expectTrue("the plate actually got hot enough to lose strength",
               th::carbonSteelYieldFactor(hottest) < 0.2);
    // Biot number 0.07: steel plating is thermally thin and has no
    // through-thickness gradient worth resolving. The bound is on the *whole
    // plate*, so it bounds any element in it whatever the layer count.
    expectTrue("the spread across the whole plate never exceeds 20 K", worstSpread < 20.0);
    expectTrue("and the spread in k_y across it never exceeds 0.04", worstFactorSpread < 0.04);
    expectTrue("but there was a gradient at all, so this is not vacuous", worstSpread > 1.0);
    std::printf("     12 mm plate, 900 C compartment, h = 200 W/(m^2 K): worst spread across the"
                " whole plate %.1f K, worst k_y spread %.3f, Bi = %.3f\n",
                worstSpread, worstFactorSpread,
                200.0 * 0.012 / th::carbonSteelConductivity(kC + 600.0));

    // **The same plate under the ISO 834 standard fire**, which is the ramp a
    // furnace test actually follows: `T_g = 20 + 345 log10(8 t + 1)` with `t` in
    // minutes. A gentler boundary condition gives a gentler gradient, and the
    // heating *rate* it produces is the number that decides whether creep needs
    // modelling at all -- EN 1993-1-2 §3.2.1's curves carry creep implicitly for
    // 2 to 50 K/min. The solver copies its `Problem`, so a moving ambient means
    // re-preparing; on twenty nodes that is cheap and it is what makes the claim a
    // measurement rather than a quotation.
    {
        th::Problem ramp = problem;
        th::Solver iso;
        std::vector<double> state(mesh.nodeCount(), kC + 20.0);
        double atFourteen = 0, worstRamp = 0;
        for (int step = 0; step < 1200; ++step) {
            const double minutes = (step + 1) / 60.0;
            ramp.film[0].ambient = kC + 20.0 + 345.0 * std::log10(8.0 * minutes + 1.0);
            expectTrue("the ramped fire prepares", iso.prepare(ramp, state, &why));
            expectTrue("and steps", iso.step(1.0, &why));
            state = iso.temperature();
            double lo = 1e30, hi = -1e30;
            for (double v : state) {
                lo = std::min(lo, v);
                hi = std::max(hi, v);
            }
            worstRamp = std::max(worstRamp, hi - lo);
            if (step == 839) atFourteen = hi;  // fourteen minutes
        }
        const double rate = (atFourteen - (kC + 20.0)) / 14.0;
        expectTrue("the ISO 834 curve is gentler than a step change", worstRamp < 8.0);
        expectTrue("and it heats the plate at a rate inside EN's 2-50 K/min band, so the"
                   " curves' implicit creep applies",
                   rate > 2.0 && rate < 50.0);
        expectTrue("having reached the strength-losing range by then",
                   atFourteen > kC + 550.0);
        std::printf("     ISO 834 on the same plate: %.0f C at 14 min, %.0f K/min, worst spread"
                    " %.1f K -- inside EN 1993-1-2 §3.2.1's 2-50 K/min band\n",
                    atFourteen - kC, rate, worstRamp);
    }

    // The quadrature error the per-element choice actually costs: `k_y` of the
    // element mean against the mean of `k_y` at the two Gauss levels of a linear
    // profile. Worst at the 400 C kink, which is the sharpest in Table 3.1.
    for (double spread : {20.0, 100.0}) {
        double worst = 0, at = 0;
        for (double centre = 20.0; centre <= 1200.0; centre += 0.5) {
            const double q = 1.0 / std::sqrt(3.0);
            const double mean = 0.5 * (th::carbonSteelYieldFactor(centre - 0.5 * spread * q + kC) +
                                       th::carbonSteelYieldFactor(centre + 0.5 * spread * q + kC));
            const double d = std::abs(mean - th::carbonSteelYieldFactor(centre + kC));
            if (d > worst) { worst = d; at = centre; }
        }
        if (spread == 20.0) {
            expectNear("a 20 K spread across an element costs 0.006 in k_y", worst, 0.0064,
                       5e-4);
            expectNear("and the worst place for it is the 400 C kink", at, 400.0, 1.0);
        } else {
            expectNear("a 100 K spread costs 0.032", worst, 0.0318, 5e-4);
        }
    }
}

}  // namespace

void runThermalTests() {
    std::printf("\n=== implicit thermal conduction ===\n");
    testMaterialCurves();
    testElement();
    testSteadyLinearField();
    testSteadyPlate();
    testSemiInfinite();
    testExplicitLimit();
    testUndershoot();
    testOrders();
    testEnergyAccount(false);
    testEnergyAccount(true);
    testNewtonCooling();
    testKirchhoffSteady();
    testSpikeAccount();
    testRotationInvariance();
    testBoundaryOrientation();
    testNodeOrdering();
    testVolumetricSourceProfile();
    testFin();
    testWarpedFilmFace();
    testRefusals();
    testUncovered();
    testCost();

    std::printf("\n=== hot steel: EN 1993-1-2 §3.2 strength at temperature ===\n");
    testReductionTable();
    testStressStrainCurve();
    testMaterialsAtTemperature();
    testReturnMapUnderAMovingYieldSurface();
    testFibresSoften();
    testHotStructureIsWeaker();
    testElementTemperatureField();
    testGradientThroughPlating();
}
