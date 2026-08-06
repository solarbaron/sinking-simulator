#!/usr/bin/env python3
"""Mutation testing for the remapped solid-shell GPU kernel.

Each mutant is one plausible single edit. The harness applies it, rebuilds only
what the edit touched, runs the suite under a per-mutant timeout, and classifies
the outcome.

Two things this harness does that earlier ones in this repo did not, both because
getting them wrong has cost real findings here:

  * **A crash counts as a kill.** An earlier harness counted `FAIL` lines only,
    and a crash prints none -- so it scored eight surviving mutants that had in
    fact aborted the binary. Here anything that is not "exit 0 with a summary
    line reporting zero failures" is a kill, and the reason is recorded.
  * **Every mutant has a timeout.** One mutant elsewhere in this repo dropped a
    square root in a timestep and took a one-minute suite to twenty. A GPU kernel
    can do worse than that.

`expect` is what the mutant *should* do. A mutant marked `survive` is a
deliberate control: it is a genuinely equivalent edit, and if the harness kills
it the harness is measuring something other than the mutation.

    ./mutate.py [--list] [--only N,M]
"""
import argparse
import pathlib
import subprocess
import sys
import time

ROOT = pathlib.Path(__file__).resolve().parents[2]
WG = ROOT / "engine/gpu/shaders/solidshell_forces_wg.comp"
HOST = ROOT / "engine/gpu/zone_gpu.cpp"
BUILD = ROOT / "build"
TIMEOUT = 420  # seconds; the clean suite is ~55

# (label, file, old, new, targets-to-rebuild, expect)
MUTANTS = [
    # --- synchronisation: the failure mode this mapping has and the other does not
    ("barrier after the return map removed", WG,
     "                returnMap(int(t), strain);\n            }\n            barrier();",
     "                returnMap(int(t), strain);\n            }",
     "shaders", "kill"),
    ("barrier after the Kaa phase removed", WG,
     "                sKaa[entry] = acc;\n            }\n            barrier();",
     "                sKaa[entry] = acc;\n            }",
     "shaders", "kill"),
    ("barrier after the 7x7 solve removed", WG,
     "            }\n            barrier();\n            if (sStop != 0) break;\n\n            if (t < uint(sEasCount)) sAlpha[t] += sDelta[t];",
     "            }\n            if (sStop != 0) break;\n\n            if (t < uint(sEasCount)) sAlpha[t] += sDelta[t];",
     "shaders", "kill"),
    ("barrier before the internal force removed", WG,
     "    if (t < uint(kEas)) plastic[sb + kStateEnhanced + t] = sAlpha[t];\n    barrier();",
     "    if (t < uint(kEas)) plastic[sb + kStateEnhanced + t] = sAlpha[t];",
     "shaders", "kill"),

    # --- the work split: entries or threads that the strided loops miss
    ("Kaa covers only the first 32 of its 49 entries", WG,
     "for (uint entry = t; entry < uint(kEas * kEas); entry += gl_WorkGroupSize.x) {",
     "for (uint entry = t; entry < 32u; entry += gl_WorkGroupSize.x) {",
     "shaders", "kill"),
    ("the eighth Gauss point is never updated", WG,
     "            if (t < uint(kGauss)) {\n                float strain[6];",
     "            if (t < uint(kGauss) - 1u) {\n                float strain[6];",
     "shaders", "kill"),
    ("the strain phase misses its last task", WG,
     "for (uint task = t; task < uint(kGauss * 6); task += gl_WorkGroupSize.x) {",
     "for (uint task = t; task < uint(kGauss * 6) - 1u; task += gl_WorkGroupSize.x) {",
     "shaders", "kill"),
    ("the last degree of freedom is never loaded", WG,
     "    for (uint i = t; i < uint(kDof); i += gl_WorkGroupSize.x) {\n        const uint a = i / 3u, k = i % 3u;",
     "    for (uint i = t; i < uint(kDof) - 1u; i += gl_WorkGroupSize.x) {\n        const uint a = i / 3u, k = i % 3u;",
     "shaders", "kill"),

    # --- the element's own arithmetic
    ("kRoot23 wrong by 0.4%", WG,
     "const float kRoot23 = 0.816496580927726f;",
     "const float kRoot23 = 0.819761567650437f;",
     "shaders", "kill"),
    ("the engineering-shear factor of two dropped", WG,
     "           2.0f * (v[3] * v[3] + v[4] * v[4] + v[5] * v[5]);",
     "           1.0f * (v[3] * v[3] + v[4] * v[4] + v[5] * v[5]);",
     "shaders", "kill"),
    ("Kaa loses its Gauss weight", WG,
     "                    acc += w * s;\n                }\n                sKaa[entry] = acc;",
     "                    acc += s;\n                }\n                sKaa[entry] = acc;",
     "shaders", "kill"),
    ("the dissipation loses its Gauss weight", WG,
     "            dissipation += forms[fb + kFormW + uint(gp)] * sDissipation[gp];",
     "            dissipation += sDissipation[gp];",
     "shaders", "kill"),
    ("the Newton no longer restarts from the step's own state", WG,
     "    gTrial = gStart;\n    for (int i = 0; i < 6; ++i) sStress[gp * 6 + i] = 0.0f;",
     "    for (int i = 0; i < 6; ++i) sStress[gp * 6 + i] = 0.0f;",
     "shaders", "kill"),
    ("the enhanced modes are dropped from the strain", WG,
     "                for (int k = 0; k < sEasCount; ++k)\n                    s += forms[gb + i * uint(kEas) + uint(k)] * sAlpha[k];",
     "                for (int k = 0; k < 0; ++k)\n                    s += forms[gb + i * uint(kEas) + uint(k)] * sAlpha[k];",
     "shaders", "kill"),
    ("Kaa is factored without equilibration", WG,
     "                    sScaleD[i] = 1.0f / sqrt(diagonal);",
     "                    sScaleD[i] = 1.0f;",
     "shaders", "kill"),
    ("the residual reads G transposed", WG,
     "                        inner += forms[gb + uint(i) * uint(kEas) + t] * sStress[gp * 6 + i];",
     "                        inner += forms[gb + t * 6u + uint(i)] * sStress[gp * 6 + i];",
     "shaders", "kill"),
    ("the internal force reads B transposed", WG,
     "            for (int i = 0; i < 6; ++i) s += forms[bb + uint(i) * uint(kDof) + j] * sStress[gp * 6 + i];",
     "            for (int i = 0; i < 6; ++i) s += forms[bb + j * 6u + uint(i)] * sStress[gp * 6 + i];",
     "shaders", "kill"),
    ("the force is rotated by R rather than R^T", WG,
     "        for (int k = 0; k < 3; ++k) s += sR[uint(k) * 3u + i] * sU[a * 3u + uint(k)];\n        elementForce[e * uint(kDof) + d] = -s;",
     "        for (int k = 0; k < 3; ++k) s += sR[i * 3u + uint(k)] * sU[a * 3u + uint(k)];\n        elementForce[e * uint(kDof) + d] = -s;",
     "shaders", "kill"),
    ("the polar decomposition stops after one iteration", WG,
     "    for (int i = 0; i < 6; ++i) r = 0.5 * (r + transpose(inverse(r)));",
     "    for (int i = 0; i < 1; ++i) r = 0.5 * (r + transpose(inverse(r)));",
     "shaders", "kill"),
    ("damage is overwritten rather than accumulated", WG,
     "        if (critical > 0.0f) gTrial.damage += dEps / critical;",
     "        if (critical > 0.0f) gTrial.damage = dEps / critical;",
     "shaders", "kill"),
    ("the enhanced tolerance loses the yield stress from its scale", WG,
     "        sYieldEnergy = params.yieldStrength * volume;",
     "        sYieldEnergy = volume;",
     "shaders", "kill"),
    ("a degraded element keeps its stale enhanced parameters", WG,
     "    if (t < uint(kEas)) plastic[sb + kStateEnhanced + t] = sAlpha[t];",
     "    if (t < uint(sEasCount)) plastic[sb + kStateEnhanced + t] = sAlpha[t];",
     "shaders", "kill"),
    ("a torn element stops being reported as torn", WG,
     "        plastic[sb + kStateTorn] = torn ? 1.0f : 0.0f;",
     "        plastic[sb + kStateTorn] = 0.0f;",
     "shaders", "kill"),
    ("the rest position is read where the current one is wanted", WG,
     "        sCurrent[i] = position[n * 3u + k];",
     "        sCurrent[i] = restPosition[n * 3u + k];",
     "shaders", "kill"),

    # --- the host
    ("the workgroup mapping dispatches one group per 32 elements", HOST,
     "                      d.mapping == Mapping::Workgroup\n                          ? static_cast<uint32_t>(d.elementCount)\n                          : groupsFor(d.elementCount, kElementGroup),",
     "                      groupsFor(d.elementCount, kElementGroup),",
     "all", "kill"),
    ("the two mappings load the same shader", HOST,
     '    const char* forceShader = mapping == Mapping::Workgroup ? "solidshell_forces_wg.comp.spv"\n                                                            : "solidshell_forces.comp.spv";',
     '    const char* forceShader = "solidshell_forces.comp.spv";',
     "all", "kill"),

    # --- deliberate controls: equivalent edits that MUST survive
    ("CONTROL: gl_WorkGroupSize.x written as the literal 32", WG,
     "    for (uint d = t; d < uint(kDof); d += gl_WorkGroupSize.x) {\n        const uint a = d / 3u, i = d % 3u;\n        float s = 0.0f;\n        for (int k = 0; k < 3; ++k) s += sR[i * 3u + uint(k)] * sCurrent[a * 3u + uint(k)];",
     "    for (uint d = t; d < uint(kDof); d += 32u) {\n        const uint a = d / 3u, i = d % 3u;\n        float s = 0.0f;\n        for (int k = 0; k < 3; ++k) s += sR[i * 3u + uint(k)] * sCurrent[a * 3u + uint(k)];",
     "shaders", "survive"),
    ("CONTROL: two independent shared loads swapped", WG,
     "        sRest[i] = restPosition[n * 3u + k];\n        sCurrent[i] = position[n * 3u + k];",
     "        sCurrent[i] = position[n * 3u + k];\n        sRest[i] = restPosition[n * 3u + k];",
     "shaders", "survive"),
]


def build(targets):
    """**Always builds both targets, and the `targets` argument is now only a label.**

    It used to build `shipsim_shaders` alone for a shader-only mutant, which is
    correct in isolation and wrong in sequence: after a mutant that edited the
    *host*, restoring `zone_gpu.cpp` on disk does not relink `shipsim_tests` unless
    something asks ninja to. So the next shader mutant ran against the previous
    mutant's host object. It was caught by a deliberate control being killed for a
    reason that had nothing to do with it -- which is the entire purpose of having
    controls, and the second time in this harness that a bug in the *harness* was
    the finding rather than a bug in the code under test.

    Rebuilding both costs about two seconds when nothing changed, against a suite
    run of fifty-five.
    """
    del targets
    args = ["ninja", "-C", str(BUILD), "shipsim_shaders", "shipsim_tests"]
    r = subprocess.run(args, capture_output=True, text=True, timeout=TIMEOUT)
    return r.returncode == 0, (r.stdout + r.stderr)[-2000:]


def run_suite():
    """Returns (killed, reason)."""
    try:
        r = subprocess.run([str(BUILD / "shipsim_tests")], capture_output=True, text=True,
                           timeout=TIMEOUT)
    except subprocess.TimeoutExpired:
        return True, "timed out"
    out = r.stdout + r.stderr
    summary = [l for l in out.splitlines() if "checks," in l and "failures" in l]
    # **The summary is parsed before the exit code is judged, and the order matters
    # for the report rather than for the verdict.** This suite exits 1 on a failing
    # check as well as on a crash, so keying off the exit code alone would label
    # every clean failure "exit 1" and lose the distinction between "an assertion
    # caught it" and "it took the process down". Both are kills; only one of them
    # tells you the suite is doing its job.
    if summary:
        failures = int(summary[-1].split("checks,")[1].split("failures")[0].strip())
        if failures > 0:
            names = [l.strip() for l in out.splitlines() if l.strip().startswith("FAIL")]
            return True, f"{failures} failing check(s), first: {names[0][5:75] if names else '?'}"
        if r.returncode != 0:
            return True, f"suite green but exited {r.returncode}"
        return False, "suite green"
    # No summary line at all: it never reached the end. A crash prints no FAIL line,
    # and counting only FAIL lines is exactly how an earlier harness in this repo
    # scored eight false survivors.
    return True, f"crashed before the summary (exit {r.returncode})"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--only", default="")
    args = ap.parse_args()

    if args.list:
        for i, m in enumerate(MUTANTS):
            print(f"{i:3d} [{m[5]:7s}] {m[0]}")
        return 0

    chosen = range(len(MUTANTS))
    if args.only:
        chosen = [int(x) for x in args.only.split(",")]

    originals = {WG: WG.read_text(), HOST: HOST.read_text()}
    # Applicability is checked for *every* mutant before any is applied, so a
    # stale pattern is a loud error at the start rather than a silent survivor
    # two hours in. CLAUDE.md: after a scripted edit, grep to confirm it landed.
    bad = []
    for i, (label, path, old, new, _t, _e) in enumerate(MUTANTS):
        if originals[path].count(old) != 1:
            bad.append(f"{i}: {label!r} matches {originals[path].count(old)} times")
    if bad:
        print("MUTANT PATTERNS DO NOT APPLY:")
        for b in bad:
            print("  " + b)
        return 2

    ok, log = build("all")
    if not ok:
        print("baseline build failed\n" + log)
        return 2
    killed_baseline, reason = run_suite()
    if killed_baseline:
        print(f"BASELINE IS NOT GREEN ({reason}) -- every mutant would score a kill")
        return 2
    print("baseline: suite green\n")

    results = []
    # **The restore is in a `finally`, and that is not defensive style.** The first
    # version of this harness restored the file on the line after the suite ran, so
    # an exception between the two left the mutant applied in the working tree --
    # which is precisely the accident `CLAUDE.md` records ("the mutant was left
    # applied in the tree by a harness accident for the better part of an hour and
    # nothing noticed"). It happened here too, on the very first run.
    try:
        for i in chosen:
            label, path, old, new, targets, expect = MUTANTS[i]
            began = time.time()
            path.write_text(originals[path].replace(old, new, 1))
            # Checked by the *disappearance of the original*, not by counting the
            # replacement. Counting the replacement is wrong whenever the new text is
            # a substring of something else in the file -- one mutant here removes a
            # line and leaves a statement that also appears, more deeply indented,
            # inside a branch further down, so `count(new)` was 2 and the harness
            # aborted mid-mutant. A non-matching edit must be loud, not fatal.
            after = path.read_text()
            assert after != originals[path], f"mutant {i} changed nothing"
            assert after.count(old) == originals[path].count(old) - 1, \
                f"mutant {i} did not remove the text it targets"
            built, log = build(targets)
            if not built:
                killed, reason = True, "did not compile"
            else:
                killed, reason = run_suite()
            path.write_text(originals[path])
            verdict = "KILLED " if killed else "SURVIVED"
            agrees = (killed and expect == "kill") or (not killed and expect == "survive")
            flag = "" if agrees else "   <-- UNEXPECTED"
            print(f"{i:3d} {verdict} ({time.time()-began:5.1f}s) {label}: {reason}{flag}")
            sys.stdout.flush()
            results.append((i, label, killed, expect, reason))
    finally:
        for path, text in originals.items():
            path.write_text(text)
        build("all")

    kills = sum(1 for r in results if r[2] and r[3] == "kill")
    want = sum(1 for r in results if r[3] == "kill")
    controls = [r for r in results if r[3] == "survive"]
    print(f"\n{kills} of {want} real mutants killed; "
          f"{sum(1 for c in controls if not c[2])} of {len(controls)} controls survived as they should")
    for i, label, killed, expect, reason in results:
        if killed and expect == "survive":
            print(f"  CONTROL WRONGLY KILLED {i}: {label} -- {reason}")
        if not killed and expect == "kill":
            print(f"  SURVIVOR {i}: {label}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
