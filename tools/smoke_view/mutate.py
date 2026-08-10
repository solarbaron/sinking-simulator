#!/usr/bin/env python3
"""Mutation testing for the volumetric fire and smoke pass.

Same harness discipline as `tools/zone_gpu_probe/mutate.py`, for the same reasons:
the source tree is copied **outside the repository** before anything is edited, a
crash counts as a kill, and every mutant runs under a timeout. Read that file's
docstring for why each of those exists; each cost a real finding here.

What is different about this pass, and why it needed its own set of mutants:

  * **A rendering test can be right about the wrong thing.** A previous GPU agent's
    device suite "only asked whether the plate moved" and five shader mutants
    survived it. So the mutants below are aimed at the edits that leave a picture
    which still looks like smoke: a Beer-Lambert linearised, an emission not
    weighted by its own emissivity, the two layers composited in the wrong order,
    a rotation used where its transpose was wanted.
  * **Two of the mutants target the *model*, not the shader.** `engine/gpu/smoke.cpp`
    states the same transfer integral in double and the pixel tests predict from
    it, so an edit there has to be caught by a closed form rather than by the GPU
    agreeing with itself.
  * **Five of these found real holes in the suite** and the checks that kill them
    were written afterwards. The face-culling sense (19) was found first, by a
    camera placed inside the medium; the first pass of this harness then reported
    four survivors, every one of which was a question no test asked:

      - 5, the layer *order*: every camera in the render suite happened to produce
        ascending rays once the heel was taken out of the body frame, so forcing
        `upperFirst` false left the whole suite green. The per-pixel sweep now runs
        two cameras and counts both orders.
      - 28, `expm1` against `exp(x) - 1`: the comment justifying it claimed the
        whole-spectrum integral would fail, and that was simply false -- the two
        differ by 1e-16 there. What does see it is the Rayleigh-Jeans *limit*,
        which is now asserted as a convergence.
      - 31, the prism's aspect ratio: a square footprint of the same area satisfies
        the volume, the layer split and the interface, and draws a 24 m engine room
        as a 15 m square.
      - 32, the cool layer's own volume: the ferry's lower layer never carries any
        tracer, so reading the whole gas volume for it produces the same zero. The
        state is now set by hand to ask the question the ship cannot.

`expect` is what the mutant *should* do. A mutant marked `survive` is a deliberate
control: it is a genuinely equivalent edit, or one whose effect is unreachable, and
if the harness kills it then the harness is measuring something other than the
mutation. Each control says which of those it is.

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
# The clean suite is ~100 s and no mutant here can make it much worse -- there is no
# iteration in any of this code -- but a mutant that makes the driver hang is exactly
# what a timeout is for.
TIMEOUT = 900

FRAG = "engine/gpu/shaders/smoke.frag"
VERT = "engine/gpu/shaders/smoke.vert"
HOST = "engine/gpu/smoke_gpu.cpp"
MODEL = "engine/gpu/smoke.cpp"

# (label, file, old, new, expect)
MUTANTS = [
    # --- Beer-Lambert itself ---------------------------------------------------
    ("the optical depth drops the path length", FRAG,
     "    const float tUpper = exp(-v.upper.a * upperPath);",
     "    const float tUpper = exp(-v.upper.a);",
     "kill"),
    ("Beer-Lambert is linearised", FRAG,
     "    const float tUpper = exp(-v.upper.a * upperPath);\n"
     "    const float tLower = exp(-v.lower.a * lowerPath);",
     "    const float tUpper = max(1.0 - v.upper.a * upperPath, 0.0);\n"
     "    const float tLower = max(1.0 - v.lower.a * lowerPath, 0.0);",
     "kill"),
    ("the alpha forgets the cool layer", FRAG,
     "    outColour = vec4(source, 1.0 - tUpper * tLower);",
     "    outColour = vec4(source, 1.0 - tUpper);",
     "kill"),

    # --- Kirchhoff: a transparent gas cannot glow -------------------------------
    ("the emission is not weighted by its own emissivity", FRAG,
     "    const vec3 eUpper = v.upper.rgb * (1.0 - tUpper);",
     "    const vec3 eUpper = v.upper.rgb;",
     "kill"),

    # --- the layer split, which is the whole of what a two-zone model has -------
    ("the two layers are composited in the wrong order", FRAG,
     "    const vec3 source = upperFirst ? eUpper + tUpper * eLower : eLower + tLower * eUpper;",
     "    const vec3 source = upperFirst ? eLower + tLower * eUpper : eUpper + tUpper * eLower;",
     "kill"),
    ("which layer the ray meets first is decided the wrong way round", FRAG,
     "        if (ray.z < 0.0) {\n            upperFirst = true;",
     "        if (ray.z < 0.0) {\n            upperFirst = false;",
     "kill"),
    ("the interface crossing is clamped to the wrong end", FRAG,
     "            upperPath = clamp(min(tExit, crossing), tEnter, tExit) - tEnter;",
     "            upperPath = clamp(max(tExit, crossing), tEnter, tExit) - tEnter;",
     "kill"),
    ("the ascending branch measures the descending one's path", FRAG,
     "            upperPath = tExit - clamp(max(tEnter, crossing), tEnter, tExit);",
     "            upperPath = clamp(max(tEnter, crossing), tEnter, tExit) - tEnter;",
     "kill"),
    ("the cool layer gets whatever is left of the ray rather than of the segment", FRAG,
     "    const float lowerPath = span - upperPath;",
     "    const float lowerPath = tExit - upperPath;",
     "kill"),

    # --- where the ray starts and stops -----------------------------------------
    ("the ray no longer stops at whatever is solid", FRAG,
     "    const float tExit = min(min(exit.x, exit.y), min(exit.z, maxDistance));",
     "    const float tExit = min(min(exit.x, exit.y), exit.z);",
     "kill"),
    ("the ray is allowed to start behind the eye", FRAG,
     "    const float tEnter = max(max(entry.x, entry.y), max(entry.z, 0.0));",
     "    const float tEnter = max(max(entry.x, entry.y), entry.z);",
     "kill"),
    ("the depth is read as a distance rather than inverted", FRAG,
     "    const float axial = push.forwardB.w / (depth + push.eyeA.w);",
     "    const float axial = push.forwardB.w / depth;",
     "kill"),
    ("the axial depth is used without correcting for the ray's obliquity", FRAG,
     "    const float along = max(dot(direction, forward), 1e-6);\n"
     "    const float maxDistance = axial / along;",
     "    const float maxDistance = axial;",
     "kill"),
    ("the depth texture is read from the wrong channel", FRAG,
     "    const float depth = texelFetch(opaqueDepth, ivec2(gl_FragCoord.xy), 0).x;",
     "    const float depth = texelFetch(opaqueDepth, ivec2(gl_FragCoord.xy), 0).y;",
     "kill"),

    # --- the body frame, which is what makes a heeled ship's layer horizontal ----
    ("the ray is rotated by R rather than by R transpose", FRAG,
     "    const mat3 bodyFromWorld = mat3(v.row0.xyz, v.row1.xyz, v.row2.xyz);",
     "    const mat3 bodyFromWorld = transpose(mat3(v.row0.xyz, v.row1.xyz, v.row2.xyz));",
     "kill"),
    ("the volume's translation is dropped from the ray's origin", FRAG,
     "    const vec3 origin = bodyFromWorld * (eye - translation);",
     "    const vec3 origin = bodyFromWorld * eye;",
     "kill"),

    # --- the box the medium lives in ---------------------------------------------
    ("one face of the box is wound into the next", VERT,
     "    2, 6, 7,  2, 7, 3,      // +y",
     "    2, 6, 7,  2, 7, 6,      // +y",
     "kill"),
    ("the box corners take the low bound on z whatever the index says", VERT,
     "                     (c & 4) != 0 ? v.hi.z : v.loInterface.z);",
     "                     v.loInterface.z);",
     "kill"),
    ("the box is placed by the transpose of its own rotation", VERT,
     "    vec3 world = vec3(dot(v.row0.xyz, body), dot(v.row1.xyz, body), dot(v.row2.xyz, body)) +",
     "    vec3 world = mat3(v.row0.xyz, v.row1.xyz, v.row2.xyz) * body +",
     "kill"),

    # --- the pipeline ---------------------------------------------------------
    # This is the defect the first version of this pass shipped with. From outside
    # the box both senses give exactly one fragment per pixel and exactly the same
    # answer; only a camera inside the medium can tell them apart.
    ("the near faces are kept instead of the far ones", HOST,
     "        raster.frontFace = VK_FRONT_FACE_COUNTER_CLOCKWISE;\n"
     "        raster.lineWidth = 1.0f;\n"
     "\n"
     "        // No depth attachment in this pass",
     "        raster.frontFace = VK_FRONT_FACE_CLOCKWISE;\n"
     "        raster.lineWidth = 1.0f;\n"
     "\n"
     "        // No depth attachment in this pass",
     "kill"),
    ("nothing is culled, so every volume is applied twice", HOST,
     "        raster.cullMode = VK_CULL_MODE_FRONT_BIT;\n"
     "        raster.frontFace = VK_FRONT_FACE_COUNTER_CLOCKWISE;",
     "        raster.cullMode = VK_CULL_MODE_NONE;\n"
     "        raster.frontFace = VK_FRONT_FACE_COUNTER_CLOCKWISE;",
     "kill"),
    ("the destination is not attenuated at all", HOST,
     "        blendAttachment.dstColorBlendFactor = VK_BLEND_FACTOR_ONE_MINUS_SRC_ALPHA;",
     "        blendAttachment.dstColorBlendFactor = VK_BLEND_FACTOR_ONE;",
     "kill"),
    ("the source is treated as unpremultiplied", HOST,
     "        blendAttachment.srcColorBlendFactor = VK_BLEND_FACTOR_ONE;",
     "        blendAttachment.srcColorBlendFactor = VK_BLEND_FACTOR_SRC_ALPHA;",
     "kill"),
    ("the volumes are sorted front to back", HOST,
     "                         return sim::length2(a.centreWorld() - eye) >\n"
     "                                sim::length2(b.centreWorld() - eye);",
     "                         return sim::length2(a.centreWorld() - eye) <\n"
     "                                sim::length2(b.centreWorld() - eye);",
     "kill"),
    ("the depth constant a is never pushed", HOST,
     "    smokePush.eyeA[3] = static_cast<float>(basis.a);",
     "    smokePush.eyeA[3] = 0.0f;",
     "kill"),
    ("both layers are packed with the hot one's extinction", HOST,
     "    out.lower[3] = static_cast<float>(v.lower.extinction);",
     "    out.lower[3] = static_cast<float>(v.upper.extinction);",
     "kill"),
    ("the interface height is packed as the deckhead", HOST,
     "    out.loInterface[3] = static_cast<float>(v.interfaceZ);",
     "    out.loInterface[3] = static_cast<float>(v.hi.z);",
     "kill"),
    ("the lit pass throws its depth away, as HullRenderer's does", HOST,
     "    opaqueAttachments[1].storeOp = VK_ATTACHMENT_STORE_OP_STORE;",
     "    opaqueAttachments[1].storeOp = VK_ATTACHMENT_STORE_OP_DONT_CARE;",
     "kill"),

    # --- the model, which the pixel tests predict from --------------------------
    ("Planck loses the Rayleigh-Jeans tail to cancellation", MODEL,
     "    const double denominator = std::expm1(x);",
     "    const double denominator = std::exp(x) - 1.0;",
     "kill"),
    ("Stefan-Boltzmann is derived with the wrong denominator", MODEL,
     "    return 2.0 * pi5 * k4 / (15.0 * kLightSpeed * kLightSpeed * h3);",
     "    return 2.0 * pi5 * k4 / (16.0 * kLightSpeed * kLightSpeed * h3);",
     "kill"),
    ("the display bands are indexed blue-first", MODEL,
     "        out[c] = blackbodyBandRadiance(kBandEdges[2 - c], kBandEdges[3 - c], temperature,\n"
     "                                       kDisplayBandIntervals);",
     "        out[c] = blackbodyBandRadiance(kBandEdges[c], kBandEdges[c + 1], temperature,\n"
     "                                       kDisplayBandIntervals);",
     "kill"),
    ("the prism forgets the compartment's aspect ratio", MODEL,
     "        const double sy = std::sqrt(g.floorArea / aspect);\n"
     "        const double sx = aspect * sy;",
     "        const double sy = std::sqrt(g.floorArea);\n"
     "        const double sx = sy;",
     "kill"),
    ("the cool layer's concentration is taken over the whole gas volume", MODEL,
     "        v.lower = layerFrom(g.lower, std::max(g.gasVolume - upperVolume, 0.0), shading);",
     "        v.lower = layerFrom(g.lower, g.gasVolume, shading);",
     "kill"),
    ("a ray pointing away from the prism still hits it", MODEL,
     "    double t0 = 0.0, t1 = maxDistance;",
     "    double t0 = -1e30, t1 = maxDistance;",
     "kill"),
    ("the depth basis takes the sign of a from the wrong row", MODEL,
     "    out.a = -z[best] / w[best];",
     "    out.a = z[best] / w[best];",
     "kill"),

    # --- deliberate controls ---------------------------------------------------
    # **Unreachable, not equivalent.** A ray whose body-frame z component is under
    # 1e-30 is exactly horizontal, and no pixel of a rendered frame is: the set has
    # measure zero. `engine/gpu/smoke.cpp`'s copy of this branch *is* reached, by
    # `testThePathSplitsAtTheInterface`, which is why the same edit there is a kill
    # and this one is not. If this dies, some frame is aiming a ray at the horizon
    # exactly and the claim above is wrong.
    ("CONTROL: the exactly-horizontal ray picks the wrong layer", FRAG,
     "        upperPath = origin.z >= v.loInterface.w ? span : 0.0;",
     "        upperPath = origin.z >= v.loInterface.w ? 0.0 : span;",
     "survive"),
    # texelFetch takes no filtering, so the sampler's filter is never consulted.
    # If this dies, something is sampling the depth with `texture()` instead.
    ("CONTROL: the depth sampler filters linearly", HOST,
     "    samplerInfo.magFilter = VK_FILTER_NEAREST;\n"
     "    samplerInfo.minFilter = VK_FILTER_NEAREST;",
     "    samplerInfo.magFilter = VK_FILTER_LINEAR;\n"
     "    samplerInfo.minFilter = VK_FILTER_LINEAR;",
     "survive"),
    # A species loading is a mass and masses are not negative, so the floor never
    # binds. It is there because a NaN would propagate into an extinction and out
    # into every pixel, not because the state is expected to go negative.
    ("CONTROL: the species floor is removed", MODEL,
     "        const double concentration = std::max(layer.products, 0.0) / volume;",
     "        const double concentration = layer.products / volume;",
     "survive"),
    # Two writes to different fields of the same row.
    ("CONTROL: two independent packs swapped", HOST,
     "    out.upper[3] = static_cast<float>(v.upper.extinction);\n"
     "    out.lower[3] = static_cast<float>(v.lower.extinction);",
     "    out.lower[3] = static_cast<float>(v.lower.extinction);\n"
     "    out.upper[3] = static_cast<float>(v.upper.extinction);",
     "survive"),
]


def build(build_dir):
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
    if summary:
        failures = int(summary[-1].split("checks,")[1].split("failures")[0].strip())
        if failures > 0:
            names = [l.strip() for l in out.splitlines() if l.strip().startswith("FAIL")]
            return True, f"{failures} failing check(s), first: {names[0][5:85] if names else '?'}"
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
    scratch = pathlib.Path(tempfile.mkdtemp(prefix="shipsim-smoke-mutate-"))
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
