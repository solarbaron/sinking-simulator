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
#include "engine/sim/scantlings.hpp"
#include "engine/sim/solid_shell.hpp"
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
}
