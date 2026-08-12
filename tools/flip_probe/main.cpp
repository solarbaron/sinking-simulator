// SPDX-License-Identifier: MIT
//
// The sparse FLIP/APIC water solver at sizes a unit test cannot afford, and the
// place every figure `engine/sim/flip.hpp` publishes about itself comes from.
//
// `tests/test_flip.cpp` owns the closed forms; this owns the *studies*: the
// resolution sweep behind the sloshing period, the transfer comparison that chose
// APIC, the linear-solver comparison that chose conjugate gradients over the
// red-black SOR `les.cpp` uses, and the sparsity table. Each of those is several
// runs of the same thing at different settings, which is core-seconds rather than
// milliseconds and is why it is a tool.
//
//   ./flip_probe                 everything, ~2 minutes
//   ./flip_probe --transfer      PIC / FLIP / APIC
//   ./flip_probe --slosh         the first-mode period against g k tanh(k d)
//   ./flip_probe --sparse        what an empty room costs, and what water costs
//   ./flip_probe --solver        conjugate gradients against SOR, and the
//                                hydrostatic column against the tolerance
//   ./flip_probe --dam           a dam break, with the mass account
//   ./flip_probe --quick         the cheap half
#include "engine/sim/flip.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

namespace fl = sim::flip;
using sim::kGravity;
using sim::kPi;
using sim::kRhoSeawater;

namespace {

int failures = 0;

// Every parameter that produced the figures below, **including the ones nobody
// passed**. Two published tables in this repository lost the parameters that made
// them because only the overridden ones were printed, and a default that later
// moves silently re-derives the table; `tools/zone_probe` prints the same line for
// the same reason. So the whole of `flip::Params` goes out as it is constructed,
// and each study names what it overrides on top of it.
void printRun(const char* what, const char* overrides = "none") {
    const fl::Params d;
    std::printf("run: %s\n"
                "     defaults: density=%g gravity=(%g,%g,%g) affine=%d flipBlend=%g"
                " courant=%g maxSubstep=%g\n"
                "               maxSubsteps=%d projectionIterations=%d projectionTolerance=%g"
                " extrapolationDepth=%d\n"
                "               wallMargin=%g tile=%d halo=%d particlesPerCellAxis=2\n"
                "     overrides: %s\n",
                what, d.density, d.gravity[0], d.gravity[1], d.gravity[2], d.affine ? 1 : 0,
                d.flipBlend, d.courant, d.maxSubstep, d.maxSubsteps, d.projectionIterations,
                d.projectionTolerance, d.extrapolationDepth, d.wallMargin, fl::kTile, fl::kHalo,
                overrides);
}

void require(const char* what, bool condition) {
    if (!condition) {
        std::printf("  FAIL %s\n", what);
        ++failures;
    }
}

double seconds(const std::chrono::steady_clock::time_point& from) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - from).count();
}

// ---------------------------------------------------------------------------
// A quasi-two-dimensional tank and the linear standing wave in it
// ---------------------------------------------------------------------------

struct Tank {
    double period = 0, analytic = 0, peak = 0, noise = 0, decay = 0, massResidual = 0;
    int crossings = 0, fluidCells = 0, worstIterations = 0, particles = 0;
    double wall = 0;
    bool capped = false;
    // A step short of the time it was asked for. The run stops there rather than
    // publishing a period measured over a clock that stopped advancing.
    bool incomplete = false;
};

Tank runTank(double length, double depth, double h, double amplitude, double duration, double dt,
             bool affine, double blend, double tolerance = 1e-8) {
    fl::Field f;
    f.grid.h = h;
    f.grid.n[0] = static_cast<int>(std::lround(length / h));
    f.grid.n[1] = 3;
    f.grid.n[2] = static_cast<int>(std::lround(depth / h)) + 6;
    const double width = 3.0 * h;
    const double lo[3] = {0, 0, 0}, hi[3] = {length, width, depth};
    fl::seedBox(f, lo, hi, 2, kRhoSeawater);
    fl::setTotalMass(f, kRhoSeawater * length * width * depth);

    const double k = kPi / length;
    const double omega = std::sqrt(kGravity * k * std::tanh(k * depth));
    const double u0 = amplitude * kGravity * k / omega;
    for (fl::Particle& p : f.particles) {
        const double x = p.position[0], z = p.position[2];
        p.velocity[0] = u0 * std::cosh(k * z) / std::cosh(k * depth) * std::sin(k * x);
        p.velocity[2] = -u0 * std::sinh(k * z) / std::cosh(k * depth) * std::cos(k * x);
    }

    fl::Params p;
    p.affine = affine;
    p.flipBlend = blend;
    p.maxSubstep = dt;
    p.projectionTolerance = tolerance;
    fl::Account account;
    fl::resetAccount(f, account);
    fl::Solver s;

    Tank out;
    out.analytic = 2.0 * kPi / omega;
    out.particles = static_cast<int>(f.particles.size());
    const double centre = 0.5 * length;
    double c[3];
    fl::centroid(f, c);
    double previous = c[0] - centre, previousTime = 0, running = 0;
    std::vector<double> crossings, extrema;
    const auto started = std::chrono::steady_clock::now();
    for (int i = 0; i < static_cast<int>(std::lround(duration / dt)); ++i) {
        const fl::StepResult r = s.step(f, dt, p, account);
        out.worstIterations = std::max(out.worstIterations, r.projectionIterations);
        out.fluidCells = r.fluidCells;
        if (r.projectionCapped) {
            out.capped = true;
            break;
        }
        if (r.incomplete) {
            out.incomplete = true;
            break;
        }
        fl::centroid(f, c);
        const double offset = c[0] - centre;
        out.peak = std::max(out.peak, std::abs(offset));
        running = std::max(running, std::abs(offset));
        const double t = static_cast<double>(i + 1) * dt;
        if ((previous < 0 && offset >= 0) || (previous > 0 && offset <= 0)) {
            const double fraction = previous / (previous - offset);
            crossings.push_back(previousTime + fraction * (t - previousTime));
            extrema.push_back(running);
            running = 0;
        }
        previous = offset;
        previousTime = t;
    }
    out.wall = seconds(started);
    out.crossings = static_cast<int>(crossings.size());
    if (crossings.size() >= 3)
        out.period = 2.0 * (crossings.back() - crossings.front()) /
                     static_cast<double>(crossings.size() - 1);
    if (extrema.size() >= 3) out.decay = extrema.back() / extrema[1];
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
// 1. PIC, FLIP, APIC
// ---------------------------------------------------------------------------

void transferStudy(bool quick) {
    std::printf("\n=== the transfer: PIC, FLIP, APIC ===\n\n");
    printRun("the transfer study",
             "rotation: gravity=0, transfer only, no projection, h=0.05, 24^3 room;\n                affine/flipBlend per row. tank: h=0.04 d=0.32 a=0.08 L=1.0 dt=0.01\n                projectionTolerance=1e-8 duration=5 s");
    std::printf("A 0.8 m cube of water in solid-body rotation at 2 rad/s on a 0.05 m grid,\n"
                "carried through N particle->grid->particle transfers and nothing else.\n\n");
    std::printf("  %-6s %14s %14s %14s\n", "", "after 1", "after 10", "kinetic energy");
    struct Mode { const char* name; bool affine; double blend; };
    const Mode modes[3] = {{"PIC", false, 0.0}, {"FLIP", false, 1.0}, {"APIC", true, 0.0}};
    double kept[3] = {0, 0, 0};
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
        const double energy = fl::kineticEnergy(f);
        double after1 = 0;
        for (int i = 0; i < 10; ++i) {
            s.rebuild(f, p);
            s.transferToGrid(f, p);
            s.extrapolate(p);
            s.saveGrid();
            s.extrapolate(p);
            s.transferToParticles(f, p);
            if (i == 0) {
                double now[3];
                fl::angularMomentum(f, about, now);
                after1 = now[2] / before[2];
            }
        }
        double after[3];
        fl::angularMomentum(f, about, after);
        kept[m] = after[2] / before[2];
        std::printf("  %-6s %14.9f %14.9f %14.6f\n", modes[m].name, after1, kept[m],
                    fl::kineticEnergy(f) / energy);
    }
    std::printf("\n  PIC loses %.2f%% of the angular momentum per transfer, APIC %.2f%% --"
                " %.1f times less.\n",
                100.0 * (1.0 - kept[0]) / 10.0, 100.0 * (1.0 - kept[2]) / 10.0,
                (1.0 - kept[0]) / std::max(1.0 - kept[2], 1e-30));
    std::printf("  FLIP's transfer with no grid dynamics is the identity, and that is the\n"
                "  control: it is 1.0 only because `saveGrid` runs after the first\n"
                "  extrapolation. Before that fix it read %s.\n", "2.456");
    require("APIC keeps more angular momentum than PIC", kept[2] > kept[0]);
    require("FLIP's dynamics-free transfer is the identity", std::abs(kept[1] - 1.0) < 1e-12);

    if (quick) return;
    std::printf("\n  And in a tank, where the flow is not a rigid rotation:\n\n");
    std::printf("  %-6s %9s %9s %9s %8s %8s\n", "", "period", "err", "noise", "peak", "cross");
    const double L = 1.0, depth = 0.32, h = 0.04, amplitude = 0.08;
    for (int m = 0; m < 3; ++m) {
        const Tank t = runTank(L, depth, h, amplitude, 5.0, 0.01, modes[m].affine,
                               modes[m].blend);
        std::printf("  %-6s %9.5f %+8.2f%% %9.4f %8.4f %8d\n", modes[m].name, t.period,
                    100.0 * (t.period - t.analytic) / t.analytic, t.noise, t.peak, t.crossings);
        if (m == 2) require("APIC's period is inside 8%",
                            std::abs(t.period - t.analytic) < 0.08 * t.analytic);
    }
    std::printf("\n  analytic first-mode period 2 pi / sqrt(g k tanh k d) = %.5f s\n",
                2.0 * kPi / std::sqrt(kGravity * (kPi / L) * std::tanh(kPi * depth / L)));
}

// ---------------------------------------------------------------------------
// 2. The sloshing period against the dispersion relation
// ---------------------------------------------------------------------------

void sloshStudy(bool quick) {
    std::printf("\n=== the first sloshing mode against g k tanh(k d) ===\n\n");
    printRun("the sloshing study",
             "L=1.0 m, 3 cells wide, 6 cells freeboard, dt=0.01 s, duration=5 s,\n                projectionTolerance=1e-8; h, depth and amplitude per row");
    std::printf("A 1 m tank, three cells wide, seeded with the linear mode's own velocity\n"
                "field over a flat surface -- a tilted surface would be quantised to the\n"
                "cell, and below about one cell of amplitude a voxelised free surface has\n"
                "no restoring force to give at all.\n\n");
    std::printf("  %6s %6s %6s %10s %10s %9s %7s %8s %7s\n", "h", "depth", "a/h", "measured",
                "analytic", "error", "cross", "fluid", "wall");
    struct Case { double h, depth, amplitude; };
    std::vector<Case> cases = {{0.05, 0.30, 0.10}, {0.04, 0.32, 0.08}, {0.04, 0.32, 0.12},
                               {0.03125, 0.3125, 0.09375}};
    if (!quick) cases.push_back({0.025, 0.30, 0.075});
    double previousError = 1e9;
    bool monotone = true;
    for (const Case& c : cases) {
        const Tank t = runTank(1.0, c.depth, c.h, c.amplitude, 5.0, 0.01, true, 0.0);
        const double error = (t.period - t.analytic) / t.analytic;
        std::printf("  %6.4f %6.3f %6.1f %10.5f %10.5f %+8.3f%% %7d %8d %6.1fs\n", c.h, c.depth,
                    c.amplitude / c.h, t.period, t.analytic, 100.0 * error, t.crossings,
                    t.fluidCells, t.wall);
        require("the tank sloshed", t.crossings >= 3);
        require("no step was short of the time it was asked for", !t.incomplete);
        require("mass is exact through the slosh", t.massResidual == 0.0);
        // The two-cell cases are the convergence sequence; the three-cell one at
        // the same h is the amplitude control beside it.
        if (std::abs(c.amplitude / c.h - 2.0) < 0.6) {
            if (std::abs(error) > previousError) monotone = false;
            previousError = std::abs(error);
        }
    }
    std::printf("\n  At a fixed two cells of amplitude the error falls with h; raising the\n"
                "  amplitude to three cells at the same h raises it again, which is the\n"
                "  nonlinear correction and not the discretisation.\n");
    require("refinement improves the period at fixed amplitude in cells", monotone);
}

// ---------------------------------------------------------------------------
// 3. Sparsity
// ---------------------------------------------------------------------------

void sparseStudy() {
    std::printf("\n=== sparsity: what an empty room costs, and what water costs ===\n\n");
    printRun("the sparsity study", "water 0.6 x 0.6 x 0.4 m at h=0.1, dt=0.004 s, 50 steps");
    std::printf("  %10s %12s %12s %10s\n", "room", "cells", "tiles", "bytes");
    for (int extent : {20, 100, 400}) {
        fl::Field f;
        f.grid.h = 0.1;
        for (int a = 0; a < 3; ++a) f.grid.n[a] = extent;
        fl::Params p;
        fl::Solver s;
        s.rebuild(f, p);
        std::printf("  %8d^3 %12.0f %12d %10zu\n", extent,
                    static_cast<double>(extent) * extent * extent, s.tiles(), s.bytes());
        require("an empty room allocates nothing", s.tiles() == 0 && s.bytes() == 0);
    }
    std::printf("\n  The same 0.6 x 0.6 x 0.4 m of water, in rooms of three sizes:\n\n");
    std::printf("  %10s %10s %10s %12s %10s\n", "room", "tiles", "fluid", "bytes", "wall/step");
    std::vector<fl::Particle> reference;
    double perCell = 0;
    for (int extent : {20, 100, 400}) {
        fl::Field f;
        f.grid.h = 0.1;
        for (int a = 0; a < 3; ++a) f.grid.n[a] = extent;
        const double lo[3] = {0.5, 0.5, 1.0}, hi[3] = {1.1, 1.1, 1.4};
        fl::seedBox(f, lo, hi, 2, kRhoSeawater);
        fl::Params p;
        p.maxSubstep = 0.004;
        fl::Account account;
        fl::resetAccount(f, account);
        fl::Solver s;
        const auto started = std::chrono::steady_clock::now();
        fl::StepResult r{};
        for (int i = 0; i < 50; ++i) r = s.step(f, 0.004, p, account);
        const double wall = seconds(started) / 50.0;
        std::printf("  %8d^3 %10d %10d %12zu %8.3f ms\n", extent, s.tiles(), r.fluidCells,
                    s.bytes(), 1e3 * wall);
        perCell = static_cast<double>(s.bytes()) /
                  std::max(1.0, static_cast<double>(s.tiles()) * fl::kTileCells);
        if (reference.empty()) {
            reference = f.particles;
        } else {
            int differing = 0;
            for (std::size_t i = 0; i < reference.size(); ++i)
                for (int a = 0; a < 3; ++a)
                    if (reference[i].position[a] != f.particles[i].position[a]) ++differing;
            require("the answer is bit identical across room sizes", differing == 0);
        }
    }
    std::printf("\n  Bit identical in all three, so the room's extent does not enter the\n"
                "  arithmetic. The sparse structure costs %.1f bytes per allocated cell, so\n"
                "  a dense grid over the largest room would be %.1f GB before a single\n"
                "  particle existed.\n",
                perCell, 400.0 * 400.0 * 400.0 * perCell / 1e9);
}

// ---------------------------------------------------------------------------
// 4. Conjugate gradients against SOR, and the hydrostatic column
// ---------------------------------------------------------------------------

void solverStudy() {
    std::printf("\n=== the pressure solve ===\n\n");
    printRun("the pressure-solve study",
             "column 6x6x24 cells at h=0.05, dt=0.002 s, 25 steps;\n                projectionTolerance per row. The SOR comparison is a standalone\n                1-D operator at relaxation 1.7, not this solver");
    std::printf("The hydrostatic column's own Poisson problem, written out standalone: K\n"
                "cells, Neumann at the floor, a Dirichlet zero in the air above the surface,\n"
                "the whole divergence in the bottom cell. `les.cpp`'s red-black SOR at its\n"
                "own relaxation of 1.7, against the Jacobi-preconditioned CG this file uses.\n\n");
    std::printf("  %6s %10s %14s %14s\n", "K", "residual", "SOR sweeps", "CG iterations");
    for (int K : {24, 64}) {
        std::vector<double> b(static_cast<std::size_t>(K), 0.0);
        b[0] = 1.0;
        const auto diagonal = [&](int i) { return i == 0 ? 1.0 : 2.0; };
        const auto residual = [&](const std::vector<double>& x) {
            double worst = 0;
            for (int i = 0; i < K; ++i) {
                double v = diagonal(i) * x[static_cast<std::size_t>(i)] -
                           b[static_cast<std::size_t>(i)];
                if (i > 0) v -= x[static_cast<std::size_t>(i - 1)];
                if (i + 1 < K) v -= x[static_cast<std::size_t>(i + 1)];
                worst = std::max(worst, std::abs(v));
            }
            return worst;
        };
        for (double target : {1e-13, 1e-15}) {
            std::vector<double> x(static_cast<std::size_t>(K), 0.0);
            const int cap = 200000;
            int sweeps = 0;
            for (; sweeps < cap; ++sweeps) {
                for (int colour = 0; colour < 2; ++colour)
                    for (int i = colour; i < K; i += 2) {
                        double total = b[static_cast<std::size_t>(i)];
                        if (i > 0) total += x[static_cast<std::size_t>(i - 1)];
                        if (i + 1 < K) total += x[static_cast<std::size_t>(i + 1)];
                        x[static_cast<std::size_t>(i)] +=
                            1.7 * (total / diagonal(i) - x[static_cast<std::size_t>(i)]);
                    }
                if (residual(x) <= target) break;
            }
            std::vector<double> p(static_cast<std::size_t>(K), 0.0), r = b,
                                z(static_cast<std::size_t>(K), 0.0),
                                d(static_cast<std::size_t>(K), 0.0),
                                q(static_cast<std::size_t>(K), 0.0);
            for (int i = 0; i < K; ++i)
                z[static_cast<std::size_t>(i)] = r[static_cast<std::size_t>(i)] / diagonal(i);
            d = z;
            double rho = 0;
            for (int i = 0; i < K; ++i)
                rho += r[static_cast<std::size_t>(i)] * z[static_cast<std::size_t>(i)];
            int iterations = 0;
            for (; iterations < 100000; ++iterations) {
                for (int i = 0; i < K; ++i) {
                    double v = diagonal(i) * d[static_cast<std::size_t>(i)];
                    if (i > 0) v -= d[static_cast<std::size_t>(i - 1)];
                    if (i + 1 < K) v -= d[static_cast<std::size_t>(i + 1)];
                    q[static_cast<std::size_t>(i)] = v;
                }
                double denominator = 0;
                for (int i = 0; i < K; ++i)
                    denominator += d[static_cast<std::size_t>(i)] * q[static_cast<std::size_t>(i)];
                if (!(std::abs(denominator) > 0)) break;
                const double alpha = rho / denominator;
                for (int i = 0; i < K; ++i) {
                    p[static_cast<std::size_t>(i)] += alpha * d[static_cast<std::size_t>(i)];
                    r[static_cast<std::size_t>(i)] -= alpha * q[static_cast<std::size_t>(i)];
                }
                if (residual(p) <= target) { ++iterations; break; }
                for (int i = 0; i < K; ++i)
                    z[static_cast<std::size_t>(i)] = r[static_cast<std::size_t>(i)] / diagonal(i);
                double next = 0;
                for (int i = 0; i < K; ++i)
                    next += r[static_cast<std::size_t>(i)] * z[static_cast<std::size_t>(i)];
                const double beta = next / rho;
                rho = next;
                for (int i = 0; i < K; ++i)
                    d[static_cast<std::size_t>(i)] =
                        z[static_cast<std::size_t>(i)] + beta * d[static_cast<std::size_t>(i)];
            }
            if (sweeps >= cap)
                std::printf("  %6d %10.0e %14s %14d\n", K, target, "not in 200000", iterations);
            else
                std::printf("  %6d %10.0e %14d %14d\n", K, target, sweeps, iterations);
            require("CG beats SOR on this problem", iterations <= sweeps);
        }
    }

    std::printf("\n  And what that buys, on the real solver: a 24-cell hydrostatic column\n"
                "  stepped 25 times at 2 ms, against the projection tolerance.\n\n");
    std::printf("  %10s %8s %14s %16s %14s\n", "tolerance", "maxit", "|v| (m/s)", "drift (m)",
                "p error");
    for (double tolerance : {1e-8, 1e-10, 1e-12, 1e-13, 1e-15}) {
        fl::Field f;
        f.grid.h = 0.05;
        f.grid.n[0] = 6; f.grid.n[1] = 6; f.grid.n[2] = 32;
        const double lo[3] = {0, 0, 0}, hi[3] = {0.30, 0.30, 1.20};
        fl::seedBox(f, lo, hi, 2, kRhoSeawater);
        fl::setTotalMass(f, kRhoSeawater * 0.30 * 0.30 * 1.20);
        fl::Params p;
        p.maxSubstep = 0.002;
        p.projectionTolerance = tolerance;
        fl::Account account;
        fl::resetAccount(f, account);
        fl::Solver s;
        double start[3];
        fl::centroid(f, start);
        double worst = 0;
        int iterations = 0;
        for (int i = 0; i < 25; ++i) {
            const fl::StepResult r = s.step(f, 0.002, p, account);
            iterations = std::max(iterations, r.projectionIterations);
            for (const fl::Particle& q : f.particles)
                for (int a = 0; a < 3; ++a) worst = std::max(worst, std::abs(q.velocity[a]));
        }
        double finish[3];
        fl::centroid(f, finish);
        double worstPressure = 0;
        for (int k = 0; k < 24; ++k) {
            const double want = kRhoSeawater * kGravity * 0.05 * (24 - k);
            worstPressure = std::max(worstPressure, std::abs(s.pressureAt(3, 3, k) - want) / want);
        }
        std::printf("  %10.0e %8d %14.4g %16.4g %14.4g\n", tolerance, iterations, worst,
                    std::abs(finish[2] - start[2]), worstPressure);
    }
    const double t = 25 * 0.002;
    std::printf("\n  Unprojected, the same column falls %.6g m and reaches %.6g m/s.\n",
                0.5 * kGravity * t * t, kGravity * t);
}

// ---------------------------------------------------------------------------
// 5. A dam break, with the mass account
// ---------------------------------------------------------------------------

void damStudy() {
    std::printf("\n=== a dam break ===\n\n");
    printRun("the dam break",
             "column 0.30 x 0.10 x 0.50 m into 1.20 x 0.10 x 0.60 m at h=0.02,\n                dt=0.002 s, 250 steps");
    fl::Field f;
    f.grid.h = 0.02;
    f.grid.n[0] = 60; f.grid.n[1] = 5; f.grid.n[2] = 30;
    const double lo[3] = {0, 0, 0}, hi[3] = {0.30, 0.10, 0.50};
    fl::seedBox(f, lo, hi, 2, kRhoSeawater);
    fl::setTotalMass(f, kRhoSeawater * 0.30 * 0.10 * 0.50);
    const double mass = f.totalMass();
    fl::Params p;
    p.maxSubstep = 0.002;
    fl::Account account;
    fl::resetAccount(f, account);
    fl::Solver s;
    std::printf("  column 0.30 x 0.50 m released into a 1.20 m tank, %zu particles,"
                " %.4f kg\n\n", f.particles.size(), mass);
    std::printf("  %8s %10s %10s %10s %10s %12s %10s\n", "t (s)", "front (m)", "peak (m/s)",
                "tiles", "fluid", "mass res", "clamps");
    const double celerity = 2.0 * std::sqrt(kGravity * 0.50);
    const auto started = std::chrono::steady_clock::now();
    for (int i = 0; i < 250; ++i) {
        const fl::StepResult r = s.step(f, 0.002, p, account);
        if ((i + 1) % 50 == 0) {
            double front = 0;
            for (const fl::Particle& q : f.particles) front = std::max(front, q.position[0]);
            std::printf("  %8.3f %10.4f %10.4f %10d %10d %12.3g %10lld\n", f.time, front,
                        r.peakSpeed, r.tiles, r.fluidCells, account.massResidual(),
                        account.clamped);
            // Ritter's dry-bed front celerity is the fastest anything in an
            // inviscid dam break can go; a front past it is a solver adding energy
            // rather than a solver resolving a wave.
            require("the front stays inside Ritter's celerity",
                    front <= 0.30 + celerity * f.time + 0.05);
        }
        require("mass is exact through the dam break", account.massResidual() == 0.0);
    }
    std::printf("\n  Ritter's dry-bed front celerity 2 sqrt(g h0) = %.4f m/s; the front never\n"
                "  passes it, which is the only closed form a three-dimensional dam break\n"
                "  has to offer -- the shallow-water solution is not the same problem.\n",
                celerity);
    std::printf("  %.1f s of wall clock for 0.5 s of water.\n", seconds(started));
    require("no particle was lost", account.particleResidual() == 0);
}

}  // namespace

int main(int argc, char** argv) {
    bool all = argc == 1;
    bool quick = false;
    bool transfer = false, slosh = false, sparse = false, solver = false, dam = false;
    for (int i = 1; i < argc; ++i) {
        const char* a = argv[i];
        if (std::strcmp(a, "--transfer") == 0) transfer = true;
        else if (std::strcmp(a, "--slosh") == 0) slosh = true;
        else if (std::strcmp(a, "--sparse") == 0) sparse = true;
        else if (std::strcmp(a, "--solver") == 0) solver = true;
        else if (std::strcmp(a, "--dam") == 0) dam = true;
        else if (std::strcmp(a, "--quick") == 0) { quick = true; all = true; }
        else {
            std::printf("usage: flip_probe [--transfer] [--slosh] [--sparse] [--solver]"
                        " [--dam] [--quick]\n");
            return 2;
        }
    }
    std::printf("shipsim -- sparse FLIP/APIC water solver probe\n");
    if (all || sparse) sparseStudy();
    if (all || solver) solverStudy();
    if (all || transfer) transferStudy(quick);
    if ((all && !quick) || slosh) sloshStudy(quick);
    if ((all && !quick) || dam) damStudy();

    if (failures == 0) {
        std::printf("\nall checks passed\n");
        return 0;
    }
    std::printf("\n%d check(s) failed\n", failures);
    return 1;
}
