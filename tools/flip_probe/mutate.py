#!/usr/bin/env python3
"""Mutation testing for the sparse FLIP/APIC water solver, `engine/sim/flip.{hpp,cpp}`.

Same harness discipline as `tools/smoke_view/mutate.py` and
`scripts/mutate-stability.py`: the source tree is copied **outside the repository**
before anything is edited, every pattern is validated against clean source before
any of them is applied, a crash counts as a kill, a negative control on clean
source runs first so an already-red suite cannot score a perfect pass, and the
originals are restored after every mutant.

Three things are different here, and each is a lesson `CLAUDE.md` records.

**1. The characteristic kill in this repository is a hang, and a wall clock is the
wrong instrument for it.** Five of 196 mutants in `fire.cpp`/`thermal.cpp` killed
the suite with *zero failing assertions*, by collapsing a substep controller and
turning a nine-second run into an hours-long one -- and whether those score as
kills or as survivors depends entirely on where the timeout lands. The answer here
is structural rather than procedural:

  * `flip.cpp` **cannot** run unboundedly. `Params::maxSubsteps` bounds the substep
    loop, `Params::projectionIterations` bounds the conjugate gradients and
    `Params::extrapolationDepth` bounds the extrapolation. A mutant that collapses
    the substep controller therefore *returns*, with `StepResult::incomplete` set.
  * `tests/test_flip.cpp` asserts the substep count against **the controller's own
    arithmetic floor** -- `courant * h / v_peak` at the peak speed the run actually
    reached -- rather than against a clock, and asserts `!incomplete` beside it. So
    mutant 30 below, which is the substep-collapse shape, fails an assertion in
    about the time the clean run takes.
  * The timeout that remains is a *backstop*, not the measuring instrument, and it
    is set from the measured baseline with a wide margin (`--timeout-factor`,
    default 20, floor 300 s) precisely so that a busy box does not manufacture
    false kills. A false kill inflates the rate and hides a real gap, which is the
    one direction of error a kill rate cannot afford.

**2. A sweep that is interrupted leaves a tree that is green because the mutant
survived.** `--scan` re-derives every substitution in this file and greps a source
tree for each, reporting any `new` text that is present or any `old` text that has
gone missing. Run it before trusting a green tree, and run it against the *repo*
as well as against the scratch copy -- an agent killed mid-sweep came back to a
clean-looking worktree with substitution 43 still in `les.cpp`.

**3. The mutants are aimed at the shapes this solver can get wrong quietly.** A
staggered grid has three axes and three components and every index in it can read
the wrong one; a sparse grid has a halo that is right until it is one cell too
small; a projection has four terms and dropping any of them still produces a
pressure field. So: sign flips, a neighbour index reading the wrong axis, an
off-by-one on the base node, a weight from the wrong kernel, a term dropped from
the projection, a halo one cell short, and a substep controller collapsed.

`expect` is what the mutant *should* do. A mutant marked `survive` is a deliberate
control: it is a genuinely equivalent edit, or one whose effect is unreachable from
this suite, and if the harness kills it then the harness is measuring something
other than the mutation. Each control says which of those it is.

    ./mutate.py [--list] [--only N,M] [--keep] [--scan PATH] [--timeout-factor F]
"""
import argparse
import pathlib
import shutil
import subprocess
import sys
import tempfile
import time

REPO = pathlib.Path(__file__).resolve().parents[2]

CPP = "engine/sim/flip.cpp"
HPP = "engine/sim/flip.hpp"

# (label, file, old, new, expect)
MUTANTS = [
    # --- the kernel, and the two identities APIC stands on ----------------------
    ("the middle spline weight breaks the partition of unity", CPP,
     "    w[1] = 0.75 - f * f;",
     "    w[1] = 0.5 - f * f;",
     "kill"),
    ("the base node is the cell rather than half a cell below it", CPP,
     "    const int base = static_cast<int>(std::floor(s - 0.5));",
     "    const int base = static_cast<int>(std::floor(s));",
     "kill"),
    ("the weights are centred on the base node instead of the middle one", CPP,
     "    const double f = s - (static_cast<double>(base) + 1.0);",
     "    const double f = s - static_cast<double>(base);",
     "kill"),
    ("the outer two weights are swapped, mirroring the kernel", CPP,
     "    w[0] = 0.5 * (0.5 - f) * (0.5 - f);\n"
     "    w[1] = 0.75 - f * f;\n"
     "    w[2] = 0.5 * (0.5 + f) * (0.5 + f);",
     "    w[0] = 0.5 * (0.5 + f) * (0.5 + f);\n"
     "    w[1] = 0.75 - f * f;\n"
     "    w[2] = 0.5 * (0.5 - f) * (0.5 - f);",
     "kill"),
    ("the node offsets lose their sign", CPP,
     "    offset[0] = -1.0 - f;\n"
     "    offset[1] = -f;\n"
     "    offset[2] = 1.0 - f;",
     "    offset[0] = 1.0 + f;\n"
     "    offset[1] = f;\n"
     "    offset[2] = -1.0 + f;",
     "kill"),
    ("the second moment is taken about the node rather than the particle", CPP,
     "    for (int m = 0; m < 3; ++m) total += w[m] * offset[m] * offset[m];",
     "    for (int m = 0; m < 3; ++m) total += w[m] * offset[m];",
     "kill"),

    # --- the staggering ---------------------------------------------------------
    ("the transfer puts a component's own axis at the cell centre", CPP,
     "                    const double s = (p.position[b] - grid_.lo[b]) / h - (b == a ? 0.0 : 0.5);\n"
     "                    base[b] = splineWeights(s, w[b], offset[b]);\n"
     "                }\n"
     "                for (int kk = 0; kk < 3; ++kk)\n"
     "                    for (int jj = 0; jj < 3; ++jj)\n"
     "                        for (int ii = 0; ii < 3; ++ii) {\n"
     "                            const double weight = w[0][ii] * w[1][jj] * w[2][kk];\n"
     "                            if (weight == 0.0) continue;\n"
     "                            const int flat = haloIndex(&halo, base[0] + ii, base[1] + jj,\n"
     "                                                       base[2] + kk);",
     "                    const double s = (p.position[b] - grid_.lo[b]) / h - (b == a ? 0.5 : 0.0);\n"
     "                    base[b] = splineWeights(s, w[b], offset[b]);\n"
     "                }\n"
     "                for (int kk = 0; kk < 3; ++kk)\n"
     "                    for (int jj = 0; jj < 3; ++jj)\n"
     "                        for (int ii = 0; ii < 3; ++ii) {\n"
     "                            const double weight = w[0][ii] * w[1][jj] * w[2][kk];\n"
     "                            if (weight == 0.0) continue;\n"
     "                            const int flat = haloIndex(&halo, base[0] + ii, base[1] + jj,\n"
     "                                                       base[2] + kk);",
     "kill"),
    ("the grid-to-particle transfer staggers the wrong axis", CPP,
     "                    const double s = (p.position[b] - grid_.lo[b]) / h - (b == a ? 0.0 : 0.5);\n"
     "                    base[b] = splineWeights(s, w[b], offset[b]);\n"
     "                }\n"
     "                for (int kk = 0; kk < 3; ++kk)\n"
     "                    for (int jj = 0; jj < 3; ++jj)\n"
     "                        for (int ii = 0; ii < 3; ++ii) {\n"
     "                            const double weight = w[0][ii] * w[1][jj] * w[2][kk];\n"
     "                            if (weight == 0.0) continue;\n"
     "                            const int flat =\n"
     "                                haloIndex(&halo, base[0] + ii, base[1] + jj, base[2] + kk);",
     "                    const double s = (p.position[b] - grid_.lo[b]) / h - (b == a ? 0.5 : 0.0);\n"
     "                    base[b] = splineWeights(s, w[b], offset[b]);\n"
     "                }\n"
     "                for (int kk = 0; kk < 3; ++kk)\n"
     "                    for (int jj = 0; jj < 3; ++jj)\n"
     "                        for (int ii = 0; ii < 3; ++ii) {\n"
     "                            const double weight = w[0][ii] * w[1][jj] * w[2][kk];\n"
     "                            if (weight == 0.0) continue;\n"
     "                            const int flat =\n"
     "                                haloIndex(&halo, base[0] + ii, base[1] + jj, base[2] + kk);",
     "kill"),
    ("the sampler staggers the wrong axis", CPP,
     "            const double s = (x[b] - grid_.lo[b]) / h - (b == a ? 0.0 : 0.5);",
     "            const double s = (x[b] - grid_.lo[b]) / h - (b == a ? 0.5 : 0.0);",
     "kill"),

    # --- the sparse structure ---------------------------------------------------
    ("the tile of a negative cell index truncates towards zero", CPP,
     "int tileOf(int c) { return c >= 0 ? c / kTile : -((-c + kTile - 1) / kTile); }",
     "int tileOf(int c) { return c / kTile; }",
     "kill"),
    ("a cell's index within its tile reads x and z the wrong way round", CPP,
     "int localIndex(int lx, int ly, int lz) { return (lz * kTile + ly) * kTile + lx; }",
     "int localIndex(int lx, int ly, int lz) { return (lx * kTile + ly) * kTile + lz; }",
     "kill"),
    ("the halo is one cell short of what the transfer reaches", HPP,
     "inline constexpr int kHalo = 3;",
     "inline constexpr int kHalo = 2;",
     "kill"),
    ("the halo is one tile-corner short", HPP,
     "inline constexpr int kHalo = 3;",
     "inline constexpr int kHalo = 1;",
     "kill"),
    ("two tile coordinates alias onto one key", CPP,
     "    return (a << 42) | (b << 21) | c;",
     "    return (a << 40) | (b << 20) | c;",
     "kill"),
    ("the halo box is taken from the tile rather than from the particles in it", CPP,
     "            from[a] = tileOf(base + lowLocal[at] - kHalo);\n"
     "            to[a] = tileOf(base + highLocal[at] + kHalo);",
     "            from[a] = tileOf(base + lowLocal[at]);\n"
     "            to[a] = tileOf(base + highLocal[at]);",
     "kill"),

    # --- neighbours and face classification -------------------------------------
    ("a cell's plus and minus neighbours are swapped", CPP,
     "                            int n[3] = {lx, ly, lz};\n"
     "                            n[a] += side == 0 ? -1 : 1;",
     "                            int n[3] = {lx, ly, lz};\n"
     "                            n[a] += side == 0 ? 1 : -1;",
     "kill"),
    ("a neighbour crossing a tile boundary lands at the wrong end of it", CPP,
     "                            if (n[a] < 0) {\n"
     "                                host = adjacent[a * 2];\n"
     "                                n[a] = kTile - 1;",
     "                            if (n[a] < 0) {\n"
     "                                host = adjacent[a * 2];\n"
     "                                n[a] = 0;",
     "kill"),
    ("a face is classified from the cell on the wrong side of it", CPP,
     "            const std::int32_t minus = nb_[flat * 6 + static_cast<std::size_t>(a * 2)];\n"
     "            const bool minusSolid =",
     "            const std::int32_t minus = nb_[flat * 6 + static_cast<std::size_t>(a * 2 + 1)];\n"
     "            const bool minusSolid =",
     "kill"),
    ("a face beside a solid cell is treated as an ordinary fluid face", CPP,
     "            face_[a][flat] = static_cast<std::uint8_t>(minusSolid ? Face::Solid : Face::Fluid);",
     "            face_[a][flat] = static_cast<std::uint8_t>(Face::Fluid);",
     "kill"),
    # **This entered the file as a real mutant and survived, and it is the harness
    # working rather than a gap.** `fnb_` of -1 (air, a Dirichlet zero) and -2
    # (solid, the term drops) are distinguished nowhere: `apply` subtracts only for
    # `other >= 0`, and the diagonal is counted from `nonSolid` above this branch,
    # before either value is written. With `singular_` still cleared the edit is
    # exactly equivalent, so it is a control now and the real defect it was meant to
    # stand for is the mutant below it.
    ("an air neighbour is recorded as solid, with the surface still noticed", CPP,
     "                fnb_[s * 6 + static_cast<std::size_t>(d)] = -1;\n"
     "                singular_ = false;",
     "                fnb_[s * 6 + static_cast<std::size_t>(d)] = -2;\n"
     "                singular_ = false;",
     "survive"),
    ("a free surface is never noticed, so its pressure is fixed up to a constant", CPP,
     "                fnb_[s * 6 + static_cast<std::size_t>(d)] = -1;\n"
     "                singular_ = false;\n"
     "            }",
     "                fnb_[s * 6 + static_cast<std::size_t>(d)] = -1;\n"
     "            }",
     "kill"),

    # --- the projection ---------------------------------------------------------
    ("the divergence is taken with the wrong sign", CPP,
     "            const double hi = plus < 0 ? 0.0 : vel_[a][static_cast<std::size_t>(plus)];\n"
     "            divergence += hi - vel_[a][flat];\n"
     "        }\n"
     "        divergence /= h;",
     "            const double hi = plus < 0 ? 0.0 : vel_[a][static_cast<std::size_t>(plus)];\n"
     "            divergence += vel_[a][flat] - hi;\n"
     "        }\n"
     "        divergence /= h;",
     "kill"),
    ("the right-hand side loses a factor of the cell size", CPP,
     "    const double scale = params.density * h * h / dt;",
     "    const double scale = params.density * h / dt;",
     "kill"),
    ("the pressure gradient is added instead of subtracted", CPP,
     "            vel_[a][flat] -= gradient * (pressure_[flat] - low);",
     "            vel_[a][flat] += gradient * (pressure_[flat] - low);",
     "kill"),
    ("the pressure gradient is taken across the face backwards", CPP,
     "            vel_[a][flat] -= gradient * (pressure_[flat] - low);",
     "            vel_[a][flat] -= gradient * (low - pressure_[flat]);",
     "kill"),
    ("the gradient scale carries an extra cell size", CPP,
     "    const double gradient = dt / (params.density * h);",
     "    const double gradient = dt / (params.density * h * h);",
     "kill"),
    ("the Laplacian's diagonal counts only fluid neighbours", CPP,
     "            ++nonSolid;\n"
     "            if (cell_[static_cast<std::size_t>(other)] ==",
     "            if (cell_[static_cast<std::size_t>(other)] == static_cast<std::uint8_t>(Cell::Fluid))\n"
     "                ++nonSolid;\n"
     "            if (cell_[static_cast<std::size_t>(other)] ==",
     "kill"),
    ("the off-diagonal term of the Laplacian changes sign", CPP,
     "                if (other >= 0) total -= x[static_cast<std::size_t>(other)];",
     "                if (other >= 0) total += x[static_cast<std::size_t>(other)];",
     "kill"),
    ("a fully submerged region's compatibility condition is dropped", CPP,
     "        if (!singular_ || n == 0) return;",
     "        if (true || n == 0) return;",
     "kill"),
    ("the wall velocity is left at whatever the transfer put on it", CPP,
     "            if (face_[a][flat] == static_cast<std::uint8_t>(Face::Solid)) {\n"
     "                vel_[a][flat] = 0.0;\n"
     "                continue;\n"
     "            }\n"
     "            if (face_[a][flat] != static_cast<std::uint8_t>(Face::Fluid)) continue;\n"
     "            const std::int32_t minus",
     "            if (face_[a][flat] == static_cast<std::uint8_t>(Face::Solid)) {\n"
     "                continue;\n"
     "            }\n"
     "            if (face_[a][flat] != static_cast<std::uint8_t>(Face::Fluid)) continue;\n"
     "            const std::int32_t minus",
     "kill"),

    # --- body force and the substep controller ----------------------------------
    ("gravity is integrated twice over the step", CPP,
     "        const double delta = params.gravity[a] * dt;",
     "        const double delta = params.gravity[a] * dt * dt;",
     "kill"),
    ("gravity reaches the air side of the surface, which the projection cannot answer", CPP,
     "            if (face_[a][flat] == static_cast<std::uint8_t>(Face::Fluid))\n"
     "                vel_[a][flat] += delta;",
     "            if (face_[a][flat] != static_cast<std::uint8_t>(Face::Solid))\n"
     "                vel_[a][flat] += delta;",
     "kill"),
    # The substep-collapse shape: the Courant bound becomes a Courant *multiplier*,
    # so a fast flow drives the step towards zero instead of away from it. In
    # `fire.cpp` this shape killed the suite by hanging; here the substep budget
    # bounds it and `tests/test_flip.cpp` asserts the count against the
    # controller's own arithmetic floor, so it fails instead.
    ("the substep controller multiplies by the speed instead of dividing", CPP,
     "        if (fastest > 0) h = std::min(h, params.courant * field.grid.h / fastest);",
     "        if (fastest > 0) h = std::min(h, params.courant * field.grid.h * fastest);",
     "kill"),
    ("a substep advances the field further than it charges the clock", CPP,
     "        field.time += h;\n"
     "        remaining -= h;",
     "        field.time += h;\n"
     "        remaining -= 0.5 * h;",
     "kill"),

    # --- transfer arithmetic ----------------------------------------------------
    ("the face velocity is a momentum rather than a velocity", CPP,
     "            vel_[a][flat] = mass > 0 ? vel_[a][flat] / mass : 0.0;",
     "            vel_[a][flat] = mass > 0 ? vel_[a][flat] * mass : 0.0;",
     "kill"),
    ("the affine term is applied in cells rather than in metres", CPP,
     "                                value += (p.affine[a * 3 + 0] * offset[0][ii] +\n"
     "                                          p.affine[a * 3 + 1] * offset[1][jj] +\n"
     "                                          p.affine[a * 3 + 2] * offset[2][kk]) * h;",
     "                                value += (p.affine[a * 3 + 0] * offset[0][ii] +\n"
     "                                          p.affine[a * 3 + 1] * offset[1][jj] +\n"
     "                                          p.affine[a * 3 + 2] * offset[2][kk]);",
     "kill"),
    ("the affine matrix is read transposed", CPP,
     "                                value += (p.affine[a * 3 + 0] * offset[0][ii] +\n"
     "                                          p.affine[a * 3 + 1] * offset[1][jj] +\n"
     "                                          p.affine[a * 3 + 2] * offset[2][kk]) * h;",
     "                                value += (p.affine[0 * 3 + a] * offset[0][ii] +\n"
     "                                          p.affine[1 * 3 + a] * offset[1][jj] +\n"
     "                                          p.affine[2 * 3 + a] * offset[2][kk]) * h;",
     "kill"),
    ("D inverse is halved, so the affine field is half strength", CPP,
     "                for (int e = 0; e < 9; ++e) p.affine[e] = affine[e] * 4.0 / h;",
     "                for (int e = 0; e < 9; ++e) p.affine[e] = affine[e] * 2.0 / h;",
     "kill"),
    ("the affine gather reads one axis's offset for all three", CPP,
     "                                affine[a * 3 + 0] += share * offset[0][ii];\n"
     "                                affine[a * 3 + 1] += share * offset[1][jj];\n"
     "                                affine[a * 3 + 2] += share * offset[2][kk];",
     "                                affine[a * 3 + 0] += share * offset[0][ii];\n"
     "                                affine[a * 3 + 1] += share * offset[0][ii];\n"
     "                                affine[a * 3 + 2] += share * offset[0][ii];",
     "kill"),
    ("the FLIP blend is applied the wrong way round", CPP,
     "                p.velocity[a] = blend * (p.velocity[a] + delta[a]) + (1.0 - blend) * picked[a];",
     "                p.velocity[a] = (1.0 - blend) * (p.velocity[a] + delta[a]) + blend * picked[a];",
     "kill"),
    ("the FLIP difference is saved before the surface is extrapolated", CPP,
     "        extrapolate(params);\n"
     "        saveGrid();\n"
     "        addBodyForce(h, params);",
     "        saveGrid();\n"
     "        extrapolate(params);\n"
     "        addBodyForce(h, params);",
     "kill"),

    # --- extrapolation and advection --------------------------------------------
    ("the air side of the surface is never filled", CPP,
     "    const int depth = std::max(params.extrapolationDepth, 0);",
     "    const int depth = 0 * std::max(params.extrapolationDepth, 0);",
     "kill"),
    ("the extrapolation seeds itself from the air rather than from the fluid", CPP,
     "            if (face_[a][flat] == static_cast<std::uint8_t>(Face::Fluid)) valid_[flat] = 1u;",
     "            if (face_[a][flat] == static_cast<std::uint8_t>(Face::Air)) valid_[flat] = 1u;",
     "kill"),
    ("the RK2 midpoint is a whole step", CPP,
     "        for (int a = 0; a < 3; ++a) mid[a] = p.position[a] + 0.5 * dt * first[a];",
     "        for (int a = 0; a < 3; ++a) mid[a] = p.position[a] + 1.0 * dt * first[a];",
     "kill"),
    ("a particle pushed off one wall is put against the other", CPP,
     "            if (p.position[a] < low) {\n"
     "                p.position[a] = low;",
     "            if (p.position[a] < low) {\n"
     "                p.position[a] = high;",
     "kill"),
    # The replacement carries a marker rather than simply deleting the line: `--scan`
    # decides a substitution is present by the *absence* of the original, and a `new`
    # that is a substring of its own `old` makes both halves of that test ambiguous.
    ("a particle keeps the velocity that drove it into a wall", CPP,
     "                if (p.velocity[a] < 0) p.velocity[a] = 0;\n"
     "                clamped = true;",
     "                clamped = true;  /* mutant: the velocity is kept */",
     "kill"),

    # --- the exact statements the escalation will stand on ----------------------
    ("the last particle takes a share rather than the remainder", CPP,
     "    field.particles[n - 1].mass = mass - sum.total();",
     "    field.particles[n - 1].mass = share;",
     "kill"),
    ("the still-water level forgets the floor it is measured from", CPP,
     "    return floorZ + field.totalMass() / (density * planArea);",
     "    return field.totalMass() / (density * planArea);",
     "kill"),

    # --- controls: these are supposed to survive --------------------------------
    # Genuinely equivalent: the mass total is the same numbers in the same order at
    # every step, so a naive sum is as constant as a compensated one and the
    # residual is 0.0 either way. If this dies, the mass tests are measuring the
    # summation and not the conservation.
    ("the mass total drops its compensated summation", CPP,
     "double Field::totalMass() const {\n"
     "    Accumulator total;\n"
     "    for (const Particle& p : particles) total.add(p.mass);\n"
     "    return total.total();\n"
     "}",
     "double Field::totalMass() const {\n"
     "    double total = 0;\n"
     "    for (const Particle& p : particles) total += p.mass;\n"
     "    return total;\n"
     "}",
     "survive"),
    # Unreachable: no test comes within three orders of magnitude of the substep
    # budget except the one that sets it to three explicitly.
    ("the substep budget is doubled", HPP,
     "    int    maxSubsteps = 4096;",
     "    int    maxSubsteps = 8192;",
     "survive"),
    # Unreachable from the suite: nothing here has a fully submerged region *and* a
    # free surface at once, so the flag is set on the first air neighbour either
    # way. Kept as a control because it looks like a real edit.
    ("the singular flag is cleared before the loop rather than inside it", CPP,
     "    singular_ = true;\n"
     "    for (std::size_t s = 0; s < n; ++s) {",
     "    singular_ = (n > 0);\n"
     "    for (std::size_t s = 0; s < n; ++s) {",
     "survive"),
    # Equivalent within the tests' reach: the margin only decides how far inside a
    # wall a clamped particle is put, 4e-5 m at the dam break's cell size, and no
    # assertion is that sharp about a clamped particle's position.
    ("the wall margin is doubled", HPP,
     "    double wallMargin = 1e-3;",
     "    double wallMargin = 2e-3;",
     "survive"),
]

# `flip_tests` is `tests/test_flip.cpp` and `tests/harness.cpp` in their own binary
# -- the same translation unit the gate compiles, not a reduced copy. Nothing
# outside `flip.{hpp,cpp}` includes it, so no substitution below can reach another
# suite, and running all forty of them per mutant multiplies fifty-one
# substitutions by five for no coverage at all.
TARGETS = ["flip_tests"]
SUITE = "flip_tests"


def build(build_dir, timeout):
    r = subprocess.run(["ninja", "-C", str(build_dir)] + TARGETS,
                       capture_output=True, text=True, timeout=timeout)
    return r.returncode == 0, (r.stdout + r.stderr)[-2000:]


def run_suite(build_dir, timeout):
    """Returns (killed, reason, seconds)."""
    began = time.time()
    try:
        r = subprocess.run([str(build_dir / SUITE)], capture_output=True, text=True,
                           timeout=timeout)
    except subprocess.TimeoutExpired:
        return True, f"TIMED OUT after {timeout:.0f}s (a hang, not a failure)", time.time() - began
    took = time.time() - began
    out = r.stdout + r.stderr
    summary = [l for l in out.splitlines() if "checks," in l and "failures" in l]
    if summary:
        failures = int(summary[-1].split("checks,")[1].split("failures")[0].strip())
        if failures > 0:
            names = [l.strip() for l in out.splitlines() if l.strip().startswith("FAIL")]
            first = names[0][5:95] if names else "?"
            return True, f"{failures} failing check(s), first: {first}", took
        if r.returncode != 0:
            return True, f"suite green but exited {r.returncode}", took
        return False, "suite green", took
    # A crash prints no FAIL line. Counting only FAIL lines is how an earlier
    # harness in this repository scored eight false survivors.
    return True, f"crashed before the summary (exit {r.returncode})", took


def scan(root):
    """Re-derive every substitution and look for leftovers in `root`.

    The whole premise of mutation testing is that a surviving mutant is
    indistinguishable from correct code by looking at test output, and that
    applies to one's own leftovers. Never trust a green tree after an interrupted
    sweep without running this.
    """
    root = pathlib.Path(root)
    problems = 0
    texts = {}
    for _label, name, _old, _new, _e in MUTANTS:
        if name not in texts:
            texts[name] = (root / name).read_text()
    for i, (label, name, old, new, _e) in enumerate(MUTANTS):
        text = texts[name]
        # The presence of the *original* is the test, not the presence of the
        # replacement: several replacements are the original with a line removed,
        # and "is the new text there" cannot tell those apart from clean source.
        count = text.count(old)
        if count == 1:
            continue
        if count == 0 and new in text:
            print(f"  LEFTOVER  {i:3d} {label}: substitution {i} is still applied in {name}")
        elif count == 0:
            print(f"  CHANGED   {i:3d} {label}: the original text has gone from {name}, and the"
                  f" replacement is not there either")
        else:
            print(f"  AMBIGUOUS {i:3d} {label}: the original text appears {count} times in {name}")
        problems += 1
    if problems == 0:
        print(f"clean: none of the {len(MUTANTS)} substitutions is present in {root}")
    return problems


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--only", default="")
    ap.add_argument("--keep", action="store_true", help="leave the scratch copy behind")
    ap.add_argument("--scan", default="", help="grep a source tree for leftover mutations")
    ap.add_argument("--timeout-factor", type=float, default=20.0,
                    help="backstop timeout as a multiple of the measured baseline")
    args = ap.parse_args()

    if args.list:
        for i, m in enumerate(MUTANTS):
            print(f"{i:3d} [{m[4]:7s}] {m[0]}")
        return 0
    if args.scan:
        return 2 if scan(args.scan) else 0

    chosen = list(range(len(MUTANTS)))
    if args.only:
        chosen = [int(x) for x in args.only.split(",")]

    scratch = pathlib.Path(tempfile.mkdtemp(prefix="shipsim-flip-mutate-"))
    work = scratch / "src"
    build_dir = scratch / "build"
    print(f"copying the source tree to {work}")
    shutil.copytree(REPO, work,
                    ignore=shutil.ignore_patterns(".git", "build*", ".claude", "*.o", "*.spv"))

    try:
        files = {name: work / name for name in {m[1] for m in MUTANTS}}
        originals = {path: path.read_text() for path in files.values()}
        # Every pattern is checked before any is applied, so a stale one is a loud
        # error at the start rather than a silent survivor an hour in.
        bad = []
        for i, (label, name, old, _new, _e) in enumerate(MUTANTS):
            count = originals[files[name]].count(old)
            if count != 1:
                bad.append(f"{i}: {label!r} matches {count} times in {name}")
        if bad:
            print("MUTANT PATTERNS DO NOT APPLY:")
            for b in bad:
                print("  " + b)
            return 2

        r = subprocess.run(["cmake", "-S", str(work), "-B", str(build_dir), "-G", "Ninja",
                            "-DCMAKE_BUILD_TYPE=RelWithDebInfo"],
                           capture_output=True, text=True, timeout=1800)
        if r.returncode != 0:
            print("cmake configure failed in the copy\n" + r.stdout[-2000:] + r.stderr[-2000:])
            return 2
        ok, log = build(build_dir, 1800)
        if not ok:
            print("baseline build failed\n" + log)
            return 2
        # The negative control on clean source, and the measurement the backstop
        # timeout is derived from.
        killed_baseline, reason, baseline = run_suite(build_dir, 3600)
        if killed_baseline:
            print(f"BASELINE IS NOT GREEN ({reason}) -- every mutant would score a kill")
            return 2
        timeout = max(300.0, args.timeout_factor * baseline)
        print(f"baseline: suite green in {baseline:.0f}s; backstop timeout {timeout:.0f}s "
              f"({args.timeout_factor:g}x)\n")
        print("The backstop is not the instrument. `flip.cpp` cannot loop unboundedly --\n"
              "maxSubsteps, projectionIterations and extrapolationDepth are hard bounds --\n"
              "and the suite asserts the substep count against the controller's own\n"
              "arithmetic floor, so a collapsed controller fails rather than hangs.\n")

        results = []
        for i in chosen:
            label, name, old, new, expect = MUTANTS[i]
            path = files[name]
            began = time.time()
            path.write_text(originals[path].replace(old, new, 1))
            after = path.read_text()
            assert after != originals[path], f"mutant {i} changed nothing"
            assert after.count(old) == originals[path].count(old) - 1, \
                f"mutant {i} did not remove the text it targets"
            built, log = build(build_dir, timeout)
            if not built:
                killed, reason, took = True, "did not compile", time.time() - began
            else:
                killed, reason, took = run_suite(build_dir, timeout)
            path.write_text(originals[path])
            verdict = "KILLED " if killed else "SURVIVED"
            agrees = (killed and expect == "kill") or (not killed and expect == "survive")
            flag = "" if agrees else "   <-- UNEXPECTED"
            hang = "  [HANG]" if "TIMED OUT" in reason else ""
            print(f"{i:3d} {verdict} ({time.time()-began:6.1f}s, suite {took:5.1f}s vs "
                  f"{baseline:.0f}s){hang} {label}: {reason}{flag}")
            sys.stdout.flush()
            results.append((i, label, killed, expect, reason))

        # Restored and verified, every time, not only on the happy path.
        leftovers = scan(work)
        if leftovers:
            print("THE SCRATCH COPY STILL CARRIES MUTATIONS -- do not trust anything above")

        kills = sum(1 for r in results if r[2] and r[3] == "kill")
        want = sum(1 for r in results if r[3] == "kill")
        controls = [r for r in results if r[3] == "survive"]
        hangs = sum(1 for r in results if "TIMED OUT" in r[4])
        print(f"\n{kills} of {want} real mutants killed "
              f"({100.0 * kills / max(want, 1):.1f}%); "
              f"{sum(1 for c in controls if not c[2])} of {len(controls)} controls survived "
              f"as they should; {hangs} killed by the backstop rather than by an assertion")
        for i, label, killed, expect, reason in results:
            if killed and expect == "survive":
                print(f"  CONTROL WRONGLY KILLED {i}: {label} -- {reason}")
            if not killed and expect == "kill":
                print(f"  SURVIVOR {i}: {label}")
        return 0 if kills == want else 1
    finally:
        if args.keep:
            print(f"scratch copy left at {scratch}")
        else:
            shutil.rmtree(scratch, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
