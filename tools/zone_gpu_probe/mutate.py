#!/usr/bin/env python3
"""Mutation testing for the remapped solid-shell GPU kernel and its fp64 variants.

Each mutant is one plausible single edit. The harness applies it, rebuilds, runs the
suite under a per-mutant timeout, and classifies the outcome.

Three things this harness does that earlier ones in this repo did not, each because
getting it wrong has cost a real finding here:

  * **It never edits the repository.** The source tree is copied to a scratch
    directory outside it and every mutant is applied, compiled and run there. Two
    harnesses in this repo have now left a mutant applied in the working tree after
    being killed before their cleanup ran -- a `finally` does not help when the
    process is killed, and the second time it happened the mutant sat in the tree for
    the better part of an hour with nothing noticing. A copy cannot have that failure
    mode at all, and it is also what makes the harness safe to run while something
    else is building the real tree.
  * **A crash counts as a kill.** An earlier harness counted `FAIL` lines only, and a
    crash prints none -- so it scored eight surviving mutants that had in fact aborted
    the binary. Here anything that is not "exit 0 with a summary line reporting zero
    failures" is a kill, and the reason is recorded.
  * **Every mutant has a timeout.** One mutant elsewhere in this repo dropped a square
    root in a timestep and took a one-minute suite to twenty. A GPU kernel can do
    worse than that, and the fp64 variants are five to ten times slower before any
    mutation.

`expect` is what the mutant *should* do. A mutant marked `survive` is a deliberate
control: it is a genuinely equivalent edit, and if the harness kills it the harness is
measuring something other than the mutation.

    ./mutate.py [--list] [--only N,M] [--keep]
"""
import argparse
import pathlib
import shutil
import subprocess
import sys
import tempfile
import time

REPO = pathlib.Path(__file__).resolve().parents[2]
TIMEOUT = 900  # seconds; the clean suite is ~75 and an fp64 mutant can be much worse

# Paths are relative to the *copy*, which is why they are strings here and only become
# paths once the copy exists.
WG = "engine/gpu/shaders/solidshell_forces_wg.comp"
HOST = "engine/gpu/zone_gpu.cpp"

# (label, file, old, new, expect)
MUTANTS = [
    # --- synchronisation: the failure mode this mapping has and the other does not
    ("barrier after the return map removed", WG,
     "                returnMap(int(t), strain);\n            }\n            barrier();",
     "                returnMap(int(t), strain);\n            }",
     "kill"),
    ("barrier after the Kaa phase removed", WG,
     "                sKaa[entry] = acc;\n            }\n#endif\n            barrier();",
     "                sKaa[entry] = acc;\n            }\n#endif",
     "kill"),
    ("barrier after the 7x7 solve removed", WG,
     "            }\n            barrier();\n            if (sStop != 0) break;\n\n            if (t < uint(sEasCount)) sAlpha[t] += alphareal(sDelta[t]);",
     "            }\n            if (sStop != 0) break;\n\n            if (t < uint(sEasCount)) sAlpha[t] += alphareal(sDelta[t]);",
     "kill"),
    ("barrier before the internal force removed", WG,
     "    if (t < uint(kEas)) plastic[sb + kStateEnhanced + t] = sAlpha[t];\n#endif\n    barrier();",
     "    if (t < uint(kEas)) plastic[sb + kStateEnhanced + t] = sAlpha[t];\n#endif",
     "kill"),

    # --- the work split: entries or threads that the strided loops miss
    ("Kaa covers only the first 32 of its 49 entries", WG,
     "#else\n            for (uint entry = t; entry < uint(kEas * kEas); entry += gl_WorkGroupSize.x) {",
     "#else\n            for (uint entry = t; entry < 32u; entry += gl_WorkGroupSize.x) {",
     "kill"),
    ("the eighth Gauss point is never updated", WG,
     "            if (t < uint(kGauss)) {\n                float strain[6];",
     "            if (t < uint(kGauss) - 1u) {\n                float strain[6];",
     "kill"),
    ("the strain phase misses its last task", WG,
     "#else\n            for (uint task = t; task < uint(kGauss * 6); task += gl_WorkGroupSize.x) {",
     "#else\n            for (uint task = t; task < uint(kGauss * 6) - 1u; task += gl_WorkGroupSize.x) {",
     "kill"),
    ("the last degree of freedom is never loaded", WG,
     "    for (uint i = t; i < uint(kDof); i += gl_WorkGroupSize.x) {\n        const uint a = i / 3u, k = i % 3u;",
     "    for (uint i = t; i < uint(kDof) - 1u; i += gl_WorkGroupSize.x) {\n        const uint a = i / 3u, k = i % 3u;",
     "kill"),

    # --- the element's own arithmetic
    ("kRoot23 wrong by 0.4%", WG,
     "const float kRoot23 = 0.816496580927726f;",
     "const float kRoot23 = 0.819761567650437f;",
     "kill"),
    ("the engineering-shear factor of two dropped", WG,
     "           2.0f * (v[3] * v[3] + v[4] * v[4] + v[5] * v[5]);",
     "           1.0f * (v[3] * v[3] + v[4] * v[4] + v[5] * v[5]);",
     "kill"),
    ("Kaa loses its Gauss weight", WG,
     "                    acc += w * s;\n                }\n                sKaa[entry] = acc;\n            }\n#endif",
     "                    acc += s;\n                }\n                sKaa[entry] = acc;\n            }\n#endif",
     "kill"),
    ("the dissipation loses its Gauss weight", WG,
     "            dissipation += forms[fb + kFormW + uint(gp)] * sDissipation[gp];",
     "            dissipation += sDissipation[gp];",
     "kill"),
    ("the Newton no longer restarts from the step's own state", WG,
     "    gTrial = gStart;\n    for (int i = 0; i < 6; ++i) sStress[gp * 6 + i] = 0.0f;",
     "    for (int i = 0; i < 6; ++i) sStress[gp * 6 + i] = 0.0f;",
     "kill"),
    ("the enhanced modes are dropped from the strain", WG,
     "                for (int k = 0; k < sEasCount; ++k)\n                    s += forms[gb + i * uint(kEas) + uint(k)] * sAlpha[k];",
     "                for (int k = 0; k < 0; ++k)\n                    s += forms[gb + i * uint(kEas) + uint(k)] * sAlpha[k];",
     "kill"),
    ("Kaa is factored without equilibration", WG,
     "                    sScaleD[i] = easreal(1.0) / sqrt(diagonal);",
     "                    sScaleD[i] = easreal(1.0);",
     "kill"),
    ("the residual reads G transposed", WG,
     "                        inner += forms[gb + uint(i) * uint(kEas) + t] * sStress[gp * 6 + i];",
     "                        inner += forms[gb + t * 6u + uint(i)] * sStress[gp * 6 + i];",
     "kill"),
    ("the internal force reads B transposed", WG,
     "            for (int i = 0; i < 6; ++i) s += forms[bb + uint(i) * uint(kDof) + j] * sStress[gp * 6 + i];",
     "            for (int i = 0; i < 6; ++i) s += forms[bb + j * 6u + uint(i)] * sStress[gp * 6 + i];",
     "kill"),
    ("the force is rotated by R rather than R^T", WG,
     "        for (int k = 0; k < 3; ++k) s += sR[uint(k) * 3u + i] * sU[a * 3u + uint(k)];\n        elementForce[e * uint(kDof) + d] = -s;",
     "        for (int k = 0; k < 3; ++k) s += sR[i * 3u + uint(k)] * sU[a * 3u + uint(k)];\n        elementForce[e * uint(kDof) + d] = -s;",
     "kill"),
    ("the polar decomposition stops after one iteration", WG,
     "    for (int i = 0; i < 6; ++i) r = 0.5 * (r + transpose(inverse(r)));",
     "    for (int i = 0; i < 1; ++i) r = 0.5 * (r + transpose(inverse(r)));",
     "kill"),
    ("damage is overwritten rather than accumulated", WG,
     "        if (critical > 0.0f) gTrial.damage += dEps / critical;",
     "        if (critical > 0.0f) gTrial.damage = dEps / critical;",
     "kill"),
    ("the enhanced tolerance loses the yield stress from its scale", WG,
     "        sYieldEnergy = params.yieldStrength * volume;",
     "        sYieldEnergy = volume;",
     "kill"),
    ("a degraded element keeps its stale enhanced parameters", WG,
     "    if (t < uint(kEas)) plastic[sb + kStateEnhanced + t] = sAlpha[t];\n#endif",
     "    if (t < uint(sEasCount)) plastic[sb + kStateEnhanced + t] = sAlpha[t];\n#endif",
     "kill"),
    ("a torn element stops being reported as torn", WG,
     "        plastic[sb + kStateTorn] = torn ? 1.0f : 0.0f;",
     "        plastic[sb + kStateTorn] = 0.0f;",
     "kill"),
    ("the rest position is read where the current one is wanted", WG,
     "        sCurrent[i] = position[n * 3u + k];",
     "        sCurrent[i] = restPosition[n * 3u + k];",
     "kill"),

    # --- the fp64 enhanced block, which is what this pass exists for -----------------
    #
    # These are the edits that would make the precision ladder *look* like a
    # measurement while measuring nothing. Every one of them leaves a kernel that runs,
    # produces plausible numbers, and answers a different question from the one asked.
    ("the tight stopping rule is not tight", WG,
     "const easreal kEnhancedWorkTol = easreal(1.0e-16);",
     "const easreal kEnhancedWorkTol = easreal(1.0e-9);",
     "kill"),
    ("the tight variant keeps the float kernel's 12-iteration cap", WG,
     "#if SHIPSIM_EAS_TIGHT\nconst int   kMaxIterations = 40;",
     "#if SHIPSIM_EAS_TIGHT\nconst int   kMaxIterations = 12;",
     "kill"),
    ("alpha's fp64 state is never warm-started from the previous step", WG,
     "    if (t < uint(kEas)) sAlpha[t] = easAlpha[e * uint(kEas) + t];",
     "    if (t < uint(kEas)) sAlpha[t] = alphareal(0.0);",
     "kill"),
    ("alpha's fp64 state is never written back", WG,
     "        easAlpha[e * uint(kEas) + t] = sAlpha[t];\n        plastic[sb + kStateEnhanced + t] = float(sAlpha[t]);",
     "        plastic[sb + kStateEnhanced + t] = float(sAlpha[t]);",
     "kill"),
    ("the fp64 condensation reads its Gauss weights from the wrong offset", WG,
     "                    const double w = easForms[eb + kEasFormW + uint(gp)];\n"
     "                    double s = 0.0;",
     "                    const double w = easForms[eb + kEasFormW];\n"
     "                    double s = 0.0;",
     "kill"),
    ("the fp64 enhanced forms are indexed with the float stride", WG,
     "    const uint eb = e * kEasFormStride;",
     "    const uint eb = e * kFormStride;",
     "kill"),
    ("the fp64 strain drops the enhanced part it was widened for", WG,
     "                sStrain[task] = float(double(s) + enhanced);",
     "                sStrain[task] = s;",
     "kill"),

    # --- the host
    ("the workgroup mapping dispatches one group per 32 elements", HOST,
     "                      d.mapping == Mapping::Workgroup\n                          ? static_cast<uint32_t>(d.elementCount)\n                          : groupsFor(d.elementCount, kElementGroup),",
     "                      groupsFor(d.elementCount, kElementGroup),",
     "kill"),
    ("the two mappings load the same shader", HOST,
     '    const char* forceShader = mapping == Mapping::Workgroup ? workgroupShaderFor(eas)\n                                                            : "solidshell_forces.comp.spv";',
     '    const char* forceShader = "solidshell_forces.comp.spv";',
     "kill"),
    ("every enhanced-block precision loads the float kernel", HOST,
     '        case EasPrecision::FloatTight: return "solidshell_forces_wg_tight.comp.spv";',
     '        case EasPrecision::FloatTight: return "solidshell_forces_wg.comp.spv";',
     "kill"),
    ("the fp64 kernels are handed an unfilled enhanced-forms buffer", HOST,
     "bool needsDoubleForms(EasPrecision eas) {\n    return eas == EasPrecision::Condense || eas == EasPrecision::Newton;\n}",
     "bool needsDoubleForms(EasPrecision eas) {\n    return false && eas == EasPrecision::Condense;\n}",
     "kill"),
    # In bounds on purpose: `k * 6 + i` stays inside the 42-entry block, so this is a
    # transposed operator and not an overrun that ASan would find for free.
    ("the host uploads the enhanced-strain operator transposed", HOST,
     "                        out[kEasFormG + static_cast<std::size_t>(gp) * 6 * kEas +\n"
     "                            static_cast<std::size_t>(i) * kEas + static_cast<std::size_t>(k)] =",
     "                        out[kEasFormG + static_cast<std::size_t>(gp) * 6 * kEas +\n"
     "                            static_cast<std::size_t>(k) * 6 + static_cast<std::size_t>(i)] =",
     "kill"),
    # **Is the ladder a ladder?** If `Newton` quietly loads the float kernel that only
    # tightened the tolerance, every figure it produces is still plausible and the
    # comparison has become float against float -- §8's most expensive class of mistake,
    # and the exact mutant that survived the last pass on the two mappings.
    ("the fp64 Newton kernel is really the float one with a tighter tolerance", HOST,
     '        case EasPrecision::Newton:     return "solidshell_forces_wg_f64newton.comp.spv";',
     '        case EasPrecision::Newton:     return "solidshell_forces_wg_tight.comp.spv";',
     "kill"),
    ("shaderFloat64 is never requested of the device", HOST,
     "        wantedFeatures.shaderFloat64 = VK_TRUE;",
     "        wantedFeatures.shaderFloat64 = VK_FALSE;",
     "kill"),

    # --- deliberate controls: equivalent edits that MUST survive
    # **It is the *stride* that becomes the literal, not the bound.** Rewriting this
    # file for the fp64 pass mistyped it as `d < 32u`, which runs the loop eight times
    # past `kDof` and writes off the end of a shared array -- not an equivalent edit at
    # all, and the suite killed it. That is the third time in this repo that a control
    # behaving unexpectedly has been a defect in the *harness* rather than in the code
    # under test, and the only reason it was visible is that the control carries an
    # expectation.
    ("CONTROL: gl_WorkGroupSize.x written as the literal 32", WG,
     "    for (uint d = t; d < uint(kDof); d += gl_WorkGroupSize.x) {\n        const uint a = d / 3u, i = d % 3u;\n        float s = 0.0f;\n        for (int k = 0; k < 3; ++k) s += sR[i * 3u + uint(k)] * sCurrent[a * 3u + uint(k)];",
     "    for (uint d = t; d < uint(kDof); d += 32u) {\n        const uint a = d / 3u, i = d % 3u;\n        float s = 0.0f;\n        for (int k = 0; k < 3; ++k) s += sR[i * 3u + uint(k)] * sCurrent[a * 3u + uint(k)];",
     "survive"),
    ("CONTROL: two independent shared loads swapped", WG,
     "        sRest[i] = restPosition[n * 3u + k];\n        sCurrent[i] = position[n * 3u + k];",
     "        sCurrent[i] = position[n * 3u + k];\n        sRest[i] = restPosition[n * 3u + k];",
     "survive"),
    # alpha is zero for every element at t = 0 -- `ElementPlasticState` value-initialises
    # it and nothing has run -- so seeding the fp64 buffer from the CPU's state and
    # seeding it from nothing are the same buffer. If this one dies, the seed is being
    # reached with something in it and the claim above is wrong.
    ("CONTROL: the fp64 alpha buffer is seeded with zero rather than the CPU's state", HOST,
     "                    cpu.elementState()[el].enhanced[k];",
     "                    0.0;",
     "survive"),
]


def build(work, build_dir):
    args = ["ninja", "-C", str(build_dir), "shipsim_shaders", "shipsim_tests"]
    r = subprocess.run(args, capture_output=True, text=True, timeout=TIMEOUT)
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
    ap.add_argument("--keep", action="store_true", help="leave the scratch copy behind")
    args = ap.parse_args()

    if args.list:
        for i, m in enumerate(MUTANTS):
            print(f"{i:3d} [{m[4]:7s}] {m[0]}")
        return 0

    chosen = range(len(MUTANTS))
    if args.only:
        chosen = [int(x) for x in args.only.split(",")]

    # **The copy, outside the repository.** Nothing below ever writes inside REPO.
    scratch = pathlib.Path(tempfile.mkdtemp(prefix="shipsim-mutate-"))
    work = scratch / "src"
    build_dir = scratch / "build"
    print(f"copying the source tree to {work}")
    shutil.copytree(REPO, work,
                    ignore=shutil.ignore_patterns(".git", "build*", ".claude", "*.o", "*.spv"))

    try:
        files = {name: work / name for name in {m[1] for m in MUTANTS}}
        originals = {path: path.read_text() for path in files.values()}
        # Applicability is checked for *every* mutant before any is applied, so a stale
        # pattern is a loud error at the start rather than a silent survivor two hours
        # in. CLAUDE.md: after a scripted edit, grep to confirm it landed.
        bad = []
        for i, (label, name, old, new, _e) in enumerate(MUTANTS):
            text = originals[files[name]]
            if text.count(old) != 1:
                bad.append(f"{i}: {label!r} matches {text.count(old)} times")
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
        ok, log = build(work, build_dir)
        if not ok:
            print("baseline build failed\n" + log)
            return 2
        killed_baseline, reason = run_suite(build_dir)
        if killed_baseline:
            print(f"BASELINE IS NOT GREEN ({reason}) -- every mutant would score a kill")
            return 2
        print("baseline: suite green\n")

        results = []
        for i in chosen:
            label, name, old, new, expect = MUTANTS[i]
            path = files[name]
            began = time.time()
            path.write_text(originals[path].replace(old, new, 1))
            # Checked by the *disappearance of the original*, not by counting the
            # replacement. Counting the replacement is wrong whenever the new text is a
            # substring of something else in the file -- one mutant here removes a line
            # and leaves a statement that also appears, more deeply indented, inside a
            # branch further down, so `count(new)` was 2 and the harness aborted
            # mid-mutant. A non-matching edit must be loud, not fatal.
            after = path.read_text()
            assert after != originals[path], f"mutant {i} changed nothing"
            assert after.count(old) == originals[path].count(old) - 1, \
                f"mutant {i} did not remove the text it targets"
            built, log = build(work, build_dir)
            if not built:
                killed, reason = True, "did not compile"
            else:
                killed, reason = run_suite(build_dir)
            path.write_text(originals[path])
            verdict = "KILLED " if killed else "SURVIVED"
            agrees = (killed and expect == "kill") or (not killed and expect == "survive")
            flag = "" if agrees else "   <-- UNEXPECTED"
            print(f"{i:3d} {verdict} ({time.time()-began:5.1f}s) {label}: {reason}{flag}")
            sys.stdout.flush()
            results.append((i, label, killed, expect, reason))

        kills = sum(1 for r in results if r[2] and r[3] == "kill")
        want = sum(1 for r in results if r[3] == "kill")
        controls = [r for r in results if r[3] == "survive"]
        print(f"\n{kills} of {want} real mutants killed; "
              f"{sum(1 for c in controls if not c[2])} of {len(controls)} controls survived "
              f"as they should")
        for i, label, killed, expect, reason in results:
            if killed and expect == "survive":
                print(f"  CONTROL WRONGLY KILLED {i}: {label} -- {reason}")
            if not killed and expect == "kill":
                print(f"  SURVIVOR {i}: {label}")
        return 0
    finally:
        # There is nothing to restore in the repository, because nothing in it was ever
        # written. The scratch copy goes unless it is wanted for a post-mortem.
        if args.keep:
            print(f"\nscratch copy kept at {scratch}")
        else:
            shutil.rmtree(scratch, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
