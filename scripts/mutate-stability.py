#!/usr/bin/env python3
"""Mutation testing for the stability block of `engine/sim/ship.cpp`, and for the
verdict the two tools draw from it.

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
  * **The consumer is mutated as well as the measurement, and it needs the tools
    to be run.** `gmSlopeConverged` spent a while with no reader at all: both
    tools judged survival on `gmTransverse < 0`, which scores an unusable GM
    exactly as confidently as a usable one. The policy now lives in
    `sim::judgeStability` and the two verdict strings in `game::floodingOutcome`
    / `game::stabilityWord`, all three of which the suite compiles -- but **the
    call sites in the two `main()`s do not link into the suite at all**, so
    `run_probes` runs both tools on every mutant. See PROBES for which runs and
    why those.

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
FERRY = "game/prototype/ferry.cpp"
RAM = "tools/ram_view/main.cpp"
SHIPSIM = "game/prototype/main.cpp"

# The suite, and the two tools whose verdicts it cannot reach. `ram_view` exists
# only where Vulkan was configured, so it is asked for rather than assumed --
# `ninja` fails the whole build on an unknown target, which would score a kill on
# every mutant on a machine with no device.
TARGETS = ["shipsim_tests", "shipsim"]


def add_optional_targets(build_dir):
    for target in ("ram_view",):
        r = subprocess.run(["ninja", "-C", str(build_dir), "-t", "query", target],
                           capture_output=True, text=True)
        if r.returncode == 0:
            TARGETS.append(target)
        else:
            print(f"  {target} is not configured here -- its call site will not be probed")

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

    # --- the consumer: what a verdict may be drawn from -------------------------
    #
    # The flag had no consumer at all for a while: both tools judged survival on
    # `gmTransverse < 0` and scored a GM the ship had already flagged as unusable
    # exactly as confidently as one it had not. These are aimed at the wiring that
    # fixed that, and half of them are only reachable through the tool probes.
    ("an unresolved GM is judged by its sign after all", SHIP,
     "    if (!d.gmSlopeConverged) return StabilityJudgement::Unresolved;",
     "    if (false) return StabilityJudgement::Unresolved;",
     "kill"),
    ("the convergence flag is read inverted", SHIP,
     "    if (!d.gmSlopeConverged) return StabilityJudgement::Unresolved;",
     "    if (d.gmSlopeConverged) return StabilityJudgement::Unresolved;",
     "kill"),
    ("a GM of exactly zero counts as negative", SHIP,
     "    return d.gmTransverse < 0 ? StabilityJudgement::Negative"
     " : StabilityJudgement::Positive;",
     "    return d.gmTransverse <= 0 ? StabilityJudgement::Negative"
     " : StabilityJudgement::Positive;",
     "kill"),
    ("the sign is read the wrong way round", SHIP,
     "    return d.gmTransverse < 0 ? StabilityJudgement::Negative"
     " : StabilityJudgement::Positive;",
     "    return d.gmTransverse < 0 ? StabilityJudgement::Positive"
     " : StabilityJudgement::Negative;",
     "kill"),
    ("every resolved GM is judged positive", SHIP,
     "    return d.gmTransverse < 0 ? StabilityJudgement::Negative"
     " : StabilityJudgement::Positive;",
     "    return StabilityJudgement::Positive;",
     "kill"),

    # --- the verdict the two tools print ----------------------------------------
    ("an unresolved GM is reported as lost", FERRY,
     "            return \"UNDETERMINED - the righting arm is not linear at any angle it can \"\n"
     "                   \"be sampled at, so there is no GM to judge her by\";",
     "            return \"LOST - negative GM, loll imminent\";",
     "kill"),
    ("an unresolved GM is reported as survival", FERRY,
     "            return \"UNDETERMINED - the righting arm is not linear at any angle it can \"\n"
     "                   \"be sampled at, so there is no GM to judge her by\";",
     "            return \"SURVIVED - positive GM, deck edge dry\";",
     "kill"),
    ("the sentence drops its foundered guard", FERRY,
     "const char* floodingOutcome(const sim::Diagnostics& d) {\n"
     "    if (!d.afloat) return \"FOUNDERED\";",
     "const char* floodingOutcome(const sim::Diagnostics& d) {",
     "kill"),
    ("the loll threshold goes to sixty degrees", FERRY,
     "            return std::abs(d.heelDeg) > 20.0",
     "            return std::abs(d.heelDeg) > 60.0",
     "kill"),
    ("the loll test is blind to a list to port", FERRY,
     "            return std::abs(d.heelDeg) > 20.0",
     "            return d.heelDeg > 20.0",
     "kill"),
    ("the deck edge test is inverted", FERRY,
     "    return d.freeboardMin < 0 ? \"SURVIVED but the deck edge is under; no margin left\"",
     "    return d.freeboardMin > 0 ? \"SURVIVED but the deck edge is under; no margin left\"",
     "kill"),
    ("the one-word verdict calls an unresolved GM SURVIVED", FERRY,
     "        case sim::StabilityJudgement::Unresolved: return \"UNDETERMINED\";",
     "        case sim::StabilityJudgement::Unresolved: return \"SURVIVED\";",
     "kill"),
    # The conservative reading, and the one that is *wrong on the measurement*:
    # the flag fires where the layer is microns deep, so defaulting to LOST cries
    # wolf exactly where the water is negligible. It is a mutant and not a control
    # because the difference is observable -- `ram_view --duration=97` is a ship in
    # no more danger than she was one step earlier.
    ("the one-word verdict calls an unresolved GM LOST", FERRY,
     "        case sim::StabilityJudgement::Unresolved: return \"UNDETERMINED\";",
     "        case sim::StabilityJudgement::Unresolved: return \"LOST\";",
     "kill"),
    ("the one-word verdict drops its afloat guard", FERRY,
     "const char* stabilityWord(const sim::Diagnostics& d) {\n"
     "    if (!d.afloat) return \"LOST\";",
     "const char* stabilityWord(const sim::Diagnostics& d) {",
     "kill"),
    ("the one-word verdict swaps lost for survived", FERRY,
     "        case sim::StabilityJudgement::Negative:   return \"LOST\";",
     "        case sim::StabilityJudgement::Negative:   return \"SURVIVED\";",
     "kill"),

    # --- the call sites, which no unit test links ------------------------------
    ("ram_view goes back to reading the sign of GM", RAM,
     "                game::stabilityWord(after), after.floodwaterMass / 1000.0, after.heelDeg,",
     "                (!after.afloat || after.gmTransverse < 0) ? \"LOST\" : \"SURVIVED\","
     " after.floodwaterMass / 1000.0, after.heelDeg,",
     "kill"),
    ("ram_view stops saying why the GM it printed is not one", RAM,
     "    if (!after.gmSlopeConverged)",
     "    if (false && !after.gmSlopeConverged)",
     "kill"),
    ("shipsim stops asking for a verdict at all", SHIPSIM,
     "    if (std::string_view(outcome).starts_with(\"still\")) outcome ="
     " game::floodingOutcome(fin);",
     "    if (false) outcome = game::floodingOutcome(fin);",
     "kill"),

    # --- what `--gm-detail` publishes -------------------------------------------
    #
    # Three README figures rest on these: the GM a fixed +/-0.03 rad sample reads,
    # the angle the layer pockets at, and the fraction of the free surface that
    # sample sees. Until the flag existed they were re-derivable only by writing
    # C++ against the library, which is the least checkable position a load-bearing
    # number can be in.
    ("the halving count stops being kept", SHIP,
     "            ++d.gmHalvings;",
     "            ;",
     "kill"),
    ("the halving count includes the pass that stopped the search", SHIP,
     "                gm = gmHalved;\n"
     "                d.gmSlopeConverged = true;",
     "                gm = gmHalved;\n"
     "                ++d.gmHalvings;\n"
     "                d.gmSlopeConverged = true;",
     "kill"),
    ("the free surface reported is the first wet one rather than the largest", SHIP,
     "        if (area <= best.planArea || area <= 0 || length <= 0) continue;",
     "        if (best.compartment >= 0 || area <= 0 || length <= 0) continue;",
     "kill"),
    ("the layer depth forgets the permeability the water actually occupies", SHIP,
     "        best.depth = c.waterVolume / std::max(c.permeability, 1e-6) / area;",
     "        best.depth = c.waterVolume / area;",
     "kill"),
    ("the mean breadth becomes the widest part of the deck", SHIP,
     "        best.breadth = area / length;",
     "        best.breadth = c.bboxHi.y - c.bboxLo.y;",
     "kill"),
    ("the pocketing angle loses its factor of two", SHIP,
     "        best.pocketingRad = std::atan(2.0 * best.depth / best.breadth);",
     "        best.pocketingRad = std::atan(best.depth / best.breadth);",
     "kill"),
    ("the floor area is measured through the whole compartment", SHIP,
     "                                {c.bboxHi.x, c.bboxHi.y, c.bboxLo.z + kSlab})).volume"
     " / kSlab;",
     "                                {c.bboxHi.x, c.bboxHi.y, c.bboxHi.z})).volume / kSlab;",
     "kill"),
    ("the wedge is applied below the pocketing angle as well as above it", SHIP,
     "    if (tau <= 1.0) return 1.0;",
     "    if (tau <= 0.0) return 1.0;",
     "kill"),
    ("the pocketed moment keeps only its leading term", SHIP,
     "    return 3.0 / tau - 2.0 / (tau * std::sqrt(tau));",
     "    return 3.0 / tau;",
     "kill"),
    ("tau is taken on the angle rather than on its tangent", SHIP,
     "    const double tau = std::abs(std::tan(heelRad)) * breadth / (2.0 * depth);",
     "    const double tau = std::abs(heelRad) * breadth / (2.0 * depth);",
     "kill"),
    ("the control is taken at a tenth of the angle it names", SHIPSIM,
     "    constexpr double kFixed = 0.03;",
     "    constexpr double kFixed = 0.003;",
     "kill"),
    ("the control's difference is one-sided", SHIPSIM,
     "        (s.rightingArmAtHeel(kFixed, sea) - s.rightingArmAtHeel(-kFixed, sea))"
     " / (2 * kFixed);",
     "        (s.rightingArmAtHeel(kFixed, sea) - s.rightingArmAtHeel(0.0, sea)) / kFixed;",
     "kill"),

    # --- controls ---------------------------------------------------------------
    # Heel is positive starboard-down and the wedge is symmetric, so `tau` is the
    # same on either side; nothing in this repository asks for the moment fraction
    # at a negative angle, and the two callers that ask at all -- `--gm-detail` and
    # the closed-form test -- both ask at +0.03 rad and above. The abs is there so
    # that a caller sweeping a GZ curve to port cannot get a negative tau and a
    # nonsensical fraction, which is a case that does not exist yet.
    ("CONTROL: the moment fraction drops the abs on the heel it is asked at", SHIP,
     "    const double tau = std::abs(std::tan(heelRad)) * breadth / (2.0 * depth);",
     "    const double tau = std::tan(heelRad) * breadth / (2.0 * depth);",
     "survive"),
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
    # **The gap, stated as a control rather than left implicit.** `shipsim`'s
    # refusal branch is not reachable by any invocation this harness can afford:
    # the only states with an unresolved GM in a 900 s scenario run last about
    # 0.1 s -- t+679.4 in `none`, while the first litres reach the vehicle deck --
    # and a run has to *end* on one of those steps for the branch to print. The
    # probe therefore covers the call site (a `shipsim` that stopped asking for a
    # verdict dies) but not the wording of the refusal, which is covered instead by
    # `game::floodingOutcome`'s own tests. If this is ever killed, something has
    # been added that reaches it and the note above is out of date.
    ("CONTROL: shipsim's unresolved-GM line is reworded", SHIPSIM,
     "        std::printf(\"  effective GM unresolved: %.2f m is the best slope available,"
     " and it is\\n\"",
     "        std::printf(\"  effective GM unknown: %.2f m is the best slope available,"
     " and it is\\n\"",
     "survive"),
    # Two writes to different fields of the same struct, neither read in between.
    ("CONTROL: two independent publications swapped", SHIP,
     "    d.gmTransverse = gm;\n    d.gmSampledAtRad = eps;",
     "    d.gmSampledAtRad = eps;\n    d.gmTransverse = gm;",
     "survive"),
]


def build(build_dir, targets):
    args = ["ninja", "-C", str(build_dir), *targets]
    try:
        r = subprocess.run(args, capture_output=True, text=True, timeout=TIMEOUT)
    except subprocess.TimeoutExpired:
        return False, "build timed out"
    return r.returncode == 0, (r.stdout + r.stderr)[-2000:]


# --- The two tools, which the unit suite cannot reach --------------------------
#
# `shipsim` and `ram_view` each draw a verdict from one `Diagnostics`, and each
# used to do it with its own copy of `gmTransverse < 0`. The policy now lives in
# `game::floodingOutcome` / `game::stabilityWord`, which `shipsim_tests` compiles
# and tests -- but **the call sites do not link into the suite at all**, so a
# mutant that puts the old rule back in `main()` would survive every unit test
# there is. That is the same shape as the defect this whole exercise is about: the
# thing nothing reaches.
#
# So both tools are run on every mutant, and the runs are chosen to be the cheap
# ones that carry the verdict:
#
#   ram_view --speed=4.0 --duration=97   3.0 s, and it ends *inside* the 0.1 s
#       window where the first litres reach the vehicle deck and the GM is
#       unresolved. It is the only invocation of anything this repository ships
#       that prints UNDETERMINED, and it was found by instrumenting the run rather
#       than by guessing at it.
#   shipsim --scenario=none --duration=60   0.85 s, converged, deck edge dry: the
#       control. It cannot exercise the refusal -- no shipsim run whose *final*
#       state is unresolved is reachable in under 679 s -- but it does catch a call
#       site that stopped calling.
PROBES = [
    ("ram_view", ["--speed=4.0", "--duration=97"],
     ["outcome: UNDETERMINED -- ", "not a metacentric height"]),
    ("shipsim", ["--scenario=none", "--duration=60"],
     ["=== Outcome at t+60s: SURVIVED - positive GM, deck edge dry ==="]),
    #   shipsim --scenario=doors --duration=200 --gm-detail   3.7 s, and the only
    #       run that exercises `--gm-detail` on a wet ship: an 18 cm layer on the
    #       vehicle deck, pocketing at 1.9e-2 rad, one halving, and a fixed +/-0.03
    #       rad sample that sees 90% of the free surface and reads 0.55 m high.
    #       Every quantity is checked as a literal string, because the whole purpose
    #       of the flag is that these numbers can be re-derived.
    ("shipsim", ["--scenario=doors", "--duration=200", "--gm-detail"],
     ["gm-detail: gm_converged_m -3.324897",
      "gm-detail: gm_sampled_at_rad 1.500000e-02",
      "gm-detail: gm_halvings 1",
      "gm-detail: gm_converged yes",
      "gm-detail: gm_at_fixed_0.03rad_m -2.778259",
      "gm-detail: layer_compartment vehicle_deck",
      "gm-detail: layer_plan_area_m2 1868.365",
      "gm-detail: layer_depth_m 1.811410e-01",
      "gm-detail: layer_breadth_m 18.684",
      "gm-detail: pockets_at_rad 1.938789e-02",
      "gm-detail: fixed_sample_sees_frac 0.899655"]),
]


def run_probes(build_dir):
    """Returns (killed, reason). A tool that was not built is skipped, loudly."""
    for tool, argv, expected in PROBES:
        exe = build_dir / tool
        if not exe.exists():
            print(f"      (probe {tool} not built -- its call site is unchecked)")
            continue
        try:
            r = subprocess.run([str(exe), *argv], capture_output=True, text=True, timeout=600)
        except subprocess.TimeoutExpired:
            return True, f"{tool} timed out"
        out = r.stdout + r.stderr
        if r.returncode != 0:
            return True, f"{tool} exited {r.returncode}"
        for want in expected:
            if want not in out:
                return True, f"{tool} no longer prints {want!r}"
    return False, "probes agree"


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
        add_optional_targets(build_dir)
        ok, log = build(build_dir, TARGETS)
        if not ok:
            print("baseline build failed\n" + log)
            return 2
        # The negative control on clean source. Without it, a suite that was already
        # red would score a kill on every mutant and the pass would look perfect.
        # The probes get the same treatment for the same reason: a probe whose
        # expectation was already wrong would score a kill on every mutant too.
        killed_baseline, reason = run_suite(build_dir)
        if killed_baseline:
            print(f"BASELINE IS NOT GREEN ({reason}) -- every mutant would score a kill")
            return 2
        killed_baseline, reason = run_probes(build_dir)
        if killed_baseline:
            print(f"BASELINE PROBES DISAGREE ({reason}) -- every mutant would score a kill")
            return 2
        print("baseline: suite green, tool verdicts as published\n")

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
            built, log = build(build_dir, TARGETS)
            if not built:
                killed, reason = True, "did not compile"
            else:
                killed, reason = run_suite(build_dir)
                if not killed:
                    killed, reason = run_probes(build_dir)
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
