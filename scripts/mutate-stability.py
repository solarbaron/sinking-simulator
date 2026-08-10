#!/usr/bin/env python3
"""Mutation testing for the stability block of `engine/sim/ship.cpp`.

Same harness discipline as `tools/zone_gpu_probe/mutate.py` and
`tools/smoke_view/mutate.py`, for the same reasons: the source tree is copied
**outside the repository** before anything is edited, a crash counts as a kill,
every mutant runs under a timeout, and a negative control on clean source runs
first so that an already-red suite cannot score a perfect pass.

What is different about this target:

  * **The defect being guarded against is a silent wrong number, not a crash.**
    `Diagnostics::gmTransverse` was finite-differenced at a fixed +/-0.03 rad and
    read **+0.59 m on the ferry where her initial GM is -3.77 m**, because a 2.9 cm
    layer on a 19 m deck pockets at 0.0031 rad. Every functional test passed. So
    the mutants below are aimed at edits that leave a plausible metacentric height
    behind: a sampling angle put back, a convergence tolerance loosened, a
    refinement that reports the angle it did not use.
  * **Half of them target `rightingArmAtHeel` rather than the sampling.** The
    closed forms in `test_core.cpp` are assertions about the free-surface *moment*
    -- `rho mu b^3 l / 12` and the wedge collapse `3/tau - 2/tau^1.5` -- and those
    have to be caught by algebra rather than by the sampling agreeing with itself.
  * **The bit-identity control is a mutant too.** An intact ship must come out
    exactly as it did before any of this existed, so a mutant that makes the
    refinement run on a dry ship has to die, and it is the only thing standing
    between the gated published figures and a silent re-derivation.

**What it found on its first pass**, which is the reason a control that is expected
to *survive* is worth the two minutes it costs. Nineteen of nineteen real mutants
died and three of four controls behaved -- and the fourth, "raise the halving bound
from fifteen to twenty-four", was killed. The justification written next to that
bound said it was a safety rail nothing reaches. It was not: at fifteen halvings
the refinement bottoms out at 9.2e-7 rad, and between 10 and 100 kg of water on the
ferry's vehicle deck that is not deep enough to find an answer that is perfectly
well posed -- 100 kg came out right and was flagged unreliable, and 10 kg came out
at -3.25 m against a true -3.78. The bound is now the round-off limit (1.79e-9 rad,
where the noise costs 1.4e-7 m of slope against a 1e-6 m tolerance), the old value
is a mutant that must die, and this slot holds a control that is actually equivalent.

`expect` is what the mutant *should* do. A mutant marked `survive` is a deliberate
control: it is a genuinely equivalent edit, or one whose effect is unreachable
from this suite, and if the harness kills it then the harness is measuring
something other than the mutation. Each control says which of those it is.

    ./scripts/mutate-stability.py [--list] [--only N,M] [--keep]
"""
import argparse
import pathlib
import shutil
import subprocess
import sys
import tempfile
import time

REPO = pathlib.Path(__file__).resolve().parents[1]
# The clean suite is ~110 s. No mutant here can make it much worse -- the only
# loop any of them touches is bounded at twenty-four halvings -- but a mutant that
# turns the refinement into a non-terminating search is exactly what this is for.
TIMEOUT = 1200

SHIP = "engine/sim/ship.cpp"

# (label, file, old, new, expect)
MUTANTS = [
    # --- the sampling angle itself ---------------------------------------------
    ("the refinement is skipped entirely: back to a fixed 0.03 rad", SHIP,
     "    if (anyFreeSurface) {",
     "    if (false && anyFreeSurface) {",
     "kill"),
    ("the refinement runs on every ship, dry ones included", SHIP,
     "    if (anyFreeSurface) {",
     "    if (true || anyFreeSurface) {",
     "kill"),
    ("a compartment with water in it does not count as a free surface", SHIP,
     "        if (c.waterVolume * seaDensity > 0) { anyFreeSurface = true; break; }",
     "        if (c.waterVolume * seaDensity > 1e30) { anyFreeSurface = true; break; }",
     "kill"),
    ("the angle is quartered rather than halved each pass", SHIP,
     "            const double halved = 0.5 * eps;",
     "            const double halved = 0.25 * eps;",
     "kill"),
    ("the angle barely moves, so the loop runs out before it converges", SHIP,
     "            const double halved = 0.5 * eps;",
     "            const double halved = 0.9 * eps;",
     "kill"),
    ("two halvings instead of twenty-four", SHIP,
     "        for (int i = 0; i < 24; ++i) {",
     "        for (int i = 0; i < 2; ++i) {",
     "kill"),
    # **This entered the file as a control and was killed.** It is the design
    # error the harness found: 15 halvings bottoms out at 9.2e-7 rad, which is
    # not deep enough for a layer between 10 and 100 kg on the ferry's vehicle
    # deck -- a range where the answer is perfectly well posed. It stays as a
    # mutant so the bound cannot quietly go back.
    ("the halving bound goes back to fifteen", SHIP,
     "        for (int i = 0; i < 24; ++i) {",
     "        for (int i = 0; i < 15; ++i) {",
     "kill"),
    ("the starting angle is smaller, so a dry ship moves too", SHIP,
     "    double eps = 0.03;",
     "    double eps = 0.01;",
     "kill"),

    # --- the convergence test ---------------------------------------------------
    ("the tolerance is loosened to a centimetre of GM", SHIP,
     "            if (std::abs(gmHalved - gm) <= std::max(1e-6, 1e-4 * std::abs(gm))) {",
     "            if (std::abs(gmHalved - gm) <= std::max(1e-2, 1e-4 * std::abs(gm))) {",
     "kill"),
    ("the relative tolerance is loosened to 10%", SHIP,
     "            if (std::abs(gmHalved - gm) <= std::max(1e-6, 1e-4 * std::abs(gm))) {",
     "            if (std::abs(gmHalved - gm) <= std::max(1e-6, 1e-1 * std::abs(gm))) {",
     "kill"),
    ("the convergence test is signed, so an increasing slope never stops it", SHIP,
     "            if (std::abs(gmHalved - gm) <= std::max(1e-6, 1e-4 * std::abs(gm))) {",
     "            if (gmHalved - gm <= std::max(1e-6, 1e-4 * std::abs(gm))) {",
     "kill"),

    # --- what is reported -------------------------------------------------------
    ("the published angle is the finer of the two that agreed", SHIP,
     "    d.gmSampledAtRad = eps;",
     "    d.gmSampledAtRad = 0.5 * eps;",
     "kill"),
    ("the published angle is the one the refinement started from", SHIP,
     "    d.gmSampledAtRad = eps;",
     "    d.gmSampledAtRad = 0.03;",
     "kill"),
    ("the ship always claims its slope converged", SHIP,
     "        d.gmSlopeConverged = false;\n"
     "        for (int i = 0; i < 24; ++i) {",
     "        d.gmSlopeConverged = true;\n"
     "        for (int i = 0; i < 24; ++i) {",
     "kill"),

    # --- the difference the whole thing rests on --------------------------------
    ("the central difference is halved", SHIP,
     "        return (rightingArmAtHeel(e, sea) - rightingArmAtHeel(-e, sea)) / (2 * e);",
     "        return (rightingArmAtHeel(e, sea) - rightingArmAtHeel(-e, sea)) / (4 * e);",
     "kill"),
    ("the difference is one-sided about a non-zero arm", SHIP,
     "        return (rightingArmAtHeel(e, sea) - rightingArmAtHeel(-e, sea)) / (2 * e);",
     "        return (rightingArmAtHeel(e, sea) - rightingArmAtHeel(0.0, sea)) / e;",
     "kill"),

    # --- the free surface itself, which the closed forms are about ---------------
    ("permeability leaves the water body's geometry, so the layer is thinner", SHIP,
     "        const double region = c.waterVolume / std::max(c.permeability, 1e-6);",
     "        const double region = c.waterVolume;",
     "kill"),
    # Anchored to the following line, because `m = waterVolume * seaDensity` is
    # written three times in this file and only the one inside rightingArmAtHeel()
    # is being asked about. The harness's own pattern check found that before a
    # single mutant ran, which is what it is for.
    ("the floodwater is weighed at the geometric region rather than its own volume", SHIP,
     "        const double m = c.waterVolume * seaDensity;\n"
     "        if (m <= 0) continue;\n"
     "        const double region = c.waterVolume / std::max(c.permeability, 1e-6);",
     "        const double m = c.waterVolume * seaDensity / std::max(c.permeability, 1e-6);\n"
     "        if (m <= 0) continue;\n"
     "        const double region = c.waterVolume / std::max(c.permeability, 1e-6);",
     "kill"),
    ("the floodwater is re-levelled against the body frame, not gravity", SHIP,
     "        const double off = solvePlaneOffsetForVolume(c.mesh, up, region, -kInf, kInf);\n"
     "        const VolumeIntegral vi = integrateBelowPlane(c.mesh, up, off);",
     "        const double off = solvePlaneOffsetForVolume(c.mesh, Vec3{0, 0, 1}, region,"
     " -kInf, kInf);\n"
     "        const VolumeIntegral vi = integrateBelowPlane(c.mesh, Vec3{0, 0, 1}, off);",
     "kill"),
    ("the hull sinks to a draft that ignores the floodwater", SHIP,
     "    const double offset = solvePlaneOffsetForVolume(hull, up, mass / seaDensity,"
     " -kInf, kInf);",
     "    const double offset = solvePlaneOffsetForVolume(hull, up, lightshipMass / seaDensity,"
     " -kInf, kInf);",
     "kill"),

    # --- controls ---------------------------------------------------------------
    # seaDensity is a positive constant of the ship, so multiplying by it cannot
    # change which side of zero a volume is on. It is written that way to match the
    # test rightingArmAtHeel() itself uses, so the two can never disagree about
    # whether a compartment contributes. If this dies, they already do.
    ("CONTROL: the free-surface test drops the density it multiplies by", SHIP,
     "        if (c.waterVolume * seaDensity > 0) { anyFreeSurface = true; break; }",
     "        if (c.waterVolume > 0) { anyFreeSurface = true; break; }",
     "survive"),
    # The loop breaks out on the first wet compartment; without the break it
    # visits the rest and sets the same flag to the same value. Genuinely
    # equivalent, and it is here because the *previous* occupant of this slot --
    # "raise the halving bound" -- turned out not to be equivalent at all.
    ("CONTROL: the free-surface scan does not stop at the first wet compartment", SHIP,
     "        if (c.waterVolume * seaDensity > 0) { anyFreeSurface = true; break; }",
     "        if (c.waterVolume * seaDensity > 0) { anyFreeSurface = true; }",
     "survive"),
    # `max(1e-6, 1e-4 |GM|)` takes the relative branch whenever |GM| > 0.01 m, and
    # the smallest |GM| any ship with a free surface reaches in this suite is about
    # 0.11 m -- the box barge, whose 1.52 m of free-surface loss very nearly
    # cancels its 1.41 m of GM. So the absolute floor is unreachable here, by a
    # factor of ten and not by much more than that. It is not an equivalent edit:
    # it exists so the loop cannot chase round-off on a ship whose GM is near zero,
    # and that ship is one nothing here builds. If this ever dies, something has
    # been added that does.
    ("CONTROL: the absolute tolerance floor drops a decade", SHIP,
     "            if (std::abs(gmHalved - gm) <= std::max(1e-6, 1e-4 * std::abs(gm))) {",
     "            if (std::abs(gmHalved - gm) <= std::max(1e-7, 1e-4 * std::abs(gm))) {",
     "survive"),
    # Two writes to different fields of the same struct, neither read in between.
    ("CONTROL: two independent publications swapped", SHIP,
     "    d.gmTransverse = gm;\n    d.gmSampledAtRad = eps;",
     "    d.gmSampledAtRad = eps;\n    d.gmTransverse = gm;",
     "survive"),
]


def build(build_dir):
    args = ["ninja", "-C", str(build_dir), "shipsim_tests"]
    try:
        r = subprocess.run(args, capture_output=True, text=True, timeout=TIMEOUT)
    except subprocess.TimeoutExpired:
        return False, "build timed out"
    return r.returncode == 0, (r.stdout + r.stderr)[-2000:]


def run_suite(build_dir):
    """Returns (killed, reason)."""
    try:
        r = subprocess.run([str(build_dir / "shipsim_tests")], capture_output=True, text=True,
                           timeout=TIMEOUT)
    except subprocess.TimeoutExpired:
        return True, "timed out"
    out = r.stdout + r.stderr
    summary = [l for l in out.splitlines() if "checks," in l and "failures" in l]
    if summary:
        failures = int(summary[-1].split("checks,")[1].split("failures")[0].strip())
        if failures > 0:
            names = [l.strip() for l in out.splitlines() if l.strip().startswith("FAIL")]
            return True, f"{failures} failing check(s), first: {names[0][5:95] if names else '?'}"
        if r.returncode != 0:
            return True, f"suite green but exited {r.returncode}"
        return False, "suite green"
    # A crash prints no FAIL line. Counting only FAIL lines is how an earlier
    # harness in this repository scored eight false survivors.
    return True, f"crashed before the summary (exit {r.returncode})"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--only", default="")
    ap.add_argument("--keep", action="store_true", help="leave the scratch copy behind")
    args = ap.parse_args()

    if args.list:
        for i, m in enumerate(MUTANTS):
            print(f"{i:3d} [{m[4]:7s}] {m[0]}")
        return 0

    chosen = list(range(len(MUTANTS)))
    if args.only:
        chosen = [int(x) for x in args.only.split(",")]

    # The copy, outside the repository. Nothing below ever writes inside REPO.
    scratch = pathlib.Path(tempfile.mkdtemp(prefix="shipsim-stability-mutate-"))
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
        for i, (labelled, name, old, _new, _e) in enumerate(MUTANTS):
            text = originals[files[name]]
            if text.count(old) != 1:
                bad.append(f"{i}: {labelled!r} matches {text.count(old)} times")
        if bad:
            print("MUTANT PATTERNS DO NOT APPLY:")
            for b in bad:
                print("  " + b)
            return 2

        r = subprocess.run(["cmake", "-S", str(work), "-B", str(build_dir), "-G", "Ninja",
                            "-DCMAKE_BUILD_TYPE=RelWithDebInfo"],
                           capture_output=True, text=True, timeout=TIMEOUT)
        if r.returncode != 0:
            print("cmake configure failed in the copy\n" + r.stdout[-2000:] + r.stderr[-2000:])
            return 2
        ok, log = build(build_dir)
        if not ok:
            print("baseline build failed\n" + log)
            return 2
        # The negative control on clean source. Without it, a suite that was already
        # red would score a kill on every mutant and the pass would look perfect.
        killed_baseline, reason = run_suite(build_dir)
        if killed_baseline:
            print(f"BASELINE IS NOT GREEN ({reason}) -- every mutant would score a kill")
            return 2
        print("baseline: suite green\n")

        results = []
        for i in chosen:
            labelled, name, old, new, expect = MUTANTS[i]
            path = files[name]
            began = time.time()
            path.write_text(originals[path].replace(old, new, 1))
            after = path.read_text()
            assert after != originals[path], f"mutant {i} changed nothing"
            assert after.count(old) == originals[path].count(old) - 1, \
                f"mutant {i} did not remove the text it targets"
            built, log = build(build_dir)
            if not built:
                killed, reason = True, "did not compile"
            else:
                killed, reason = run_suite(build_dir)
            path.write_text(originals[path])
            verdict = "KILLED " if killed else "SURVIVED"
            agrees = (killed and expect == "kill") or (not killed and expect == "survive")
            flag = "" if agrees else "   <-- UNEXPECTED"
            print(f"{i:3d} {verdict} ({time.time()-began:5.1f}s) {labelled}: {reason}{flag}")
            sys.stdout.flush()
            results.append((i, labelled, killed, expect, reason))

        kills = sum(1 for r in results if r[2] and r[3] == "kill")
        want = sum(1 for r in results if r[3] == "kill")
        controls = [r for r in results if r[3] == "survive"]
        print(f"\n{kills} of {want} real mutants killed; "
              f"{sum(1 for c in controls if not c[2])} of {len(controls)} controls survived "
              f"as they should")
        for i, labelled, killed, expect, reason in results:
            if killed and expect == "survive":
                print(f"  CONTROL WRONGLY KILLED {i}: {labelled} -- {reason}")
            if not killed and expect == "kill":
                print(f"  SURVIVOR {i}: {labelled}")
        return 0
    finally:
        if args.keep:
            print(f"\nscratch copy kept at {scratch}")
        else:
            shutil.rmtree(scratch, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
