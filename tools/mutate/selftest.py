#!/usr/bin/env python3
"""The mutation harness, measured against substitutions whose answers are known.

A mutation harness is a measuring instrument and its output -- a kill rate -- is a
number nobody can check by looking at it. That is the whole premise of mutation
testing turned on the tool itself: you cannot tell a surviving mutant from correct
code by reading test output, and you cannot tell a harness that is measuring from
one that is miscounting by reading its kill rate. So this builds a toy project
whose every substitution has an answer decided in advance, and asserts the harness
returns that answer.

The toy is not arbitrary. It is a miniature of the thing that makes this
repository's mutants hard: an explicit integrator with a substep controller, where
a wrong sign does not produce a wrong number, it collapses the step and the run
stops finishing.

The known answers, one per failure mode the harness has to separate:

  * a substitution that fails an assertion            -> KILLED
  * one that aborts, printing no failing check at all -> KILLED (an earlier harness
                                                         here scored eight such as
                                                         survivors)
  * one that does not compile                         -> KILLED
  * one that compiles with a *warning* and is
    otherwise behaviour-identical                     -> KILLED, because warnings
                                                         are failures in this build
                                                         (without that rule it would
                                                         read SURVIVED)
  * one that is genuinely equivalent                  -> SURVIVED
  * one that is equivalent but 2.5x slower            -> SURVIVED, and a wall-clock
                                                         bound would have killed it
  * one that spins forever                            -> HANG (cpu bound)
  * one that blocks forever, burning no CPU           -> HANG (wall backstop) -- the
                                                         case a CPU bound alone
                                                         cannot see
  * one killed by a check that disagrees with itself  -> FLAKY, and out of the kill
                                                         count; with the re-run
                                                         switched off the same
                                                         fixture reads as a dead
                                                         control, which is the
                                                         number a harness without
                                                         it publishes
  * a catalogue pattern that no longer matches        -> refused before anything runs
  * a substitution that replaces text with nothing    -> refused at load: it could
                                                         not be told later from a
                                                         pattern that had drifted
  * a suite that is red before any mutation           -> refused, since every mutant
                                                         would otherwise score a kill

and four properties rather than outcomes:

  * SIGKILL the harness mid-sweep -> the tree is byte-identical afterwards
  * `scan` finds a mutation that *is* applied, and says a clean tree is clean
  * `manifest` notices a single changed byte
  * the *repository* is byte-identical after driving every mode over every
    catalogue in it, and a mutant's own relative write lands in the copy -- the
    docstring's central claim, which was false in two places until it was checked
  * the cheap paths -- stopping at the first failing check, and running mutants in
    parallel copies -- reach the same verdicts as the expensive ones
  * two workers cannot write the same file: each gets its own $TMPDIR, which the
    first real sweep needed after one worker's `barge.ship` scored a kill on
    another's mutant

Run it with no arguments. It writes nothing outside /tmp.
"""
import hashlib
import json
import os
import pathlib
import re
import shutil
import signal
import subprocess
import sys
import tempfile
import time

# Before `import mutate` below: importing it writes `tools/mutate/__pycache__` into
# the repository, which is the very claim the case at the bottom of this file
# measures. A checker that violates what it checks would report the tree dirty and
# be right.
sys.dont_write_bytecode = True

HERE = pathlib.Path(__file__).resolve().parent
MUTATE = HERE / "mutate.py"

checks = 0
failures = 0


def check(what, condition, detail=""):
    global checks, failures
    checks += 1
    if not condition:
        failures += 1
        print(f"  FAIL {what}" + (f" -- {detail}" if detail else ""))
    return condition


def check_equal(what, got, want):
    return check(f"{what}: got {got!r}, want {want!r}", got == want)


# --- the toy project ---------------------------------------------------------------

SIM_HPP = """\
#pragma once
namespace sim {
double relax(double start, double target, double rate, double capacity, double dt);
double mean(const double* values, int count);
double twice(double x);
double checked(double x);
}  // namespace sim
"""

SIM_CPP = """\
#include "sim.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <thread>

namespace sim {

// A relaxation towards `target` integrated explicitly, with the substep controller
// that makes this repository's mutants hard: the step is set by how fast the state
// is moving, so a wrong sign does not give a wrong answer -- the state runs away,
// the admissible step collapses, and the loop stops finishing.
double relax(double start, double target, double rate, double capacity, double dt) {
    double value = start, elapsed = 0.0;
    while (elapsed < dt) {
        const double flux = rate * (target - value);
        const double step = std::min(dt - elapsed, 0.25 * capacity * capacity /
                                                       (std::abs(flux) + 1e-300));
        value += step * flux / capacity;
        elapsed += step;
    }
    return value;
}

double mean(const double* values, int count) {
    double sum = 0.0;
    for (int i = 0; i < count; ++i) sum += values[i];
    return sum / count;
}

double twice(double x) { return 2.0 * x; }

// A guard on a state the caller is not allowed to reach, of the kind this codebase
// puts in front of a square root. Tripping it aborts, and an abort prints no
// failing check.
double checked(double x) {
    if (!(x >= 0.0)) std::abort();
    return std::sqrt(x);
}

}  // namespace sim
"""

# The toy's suite speaks the same protocol as tests/harness.cpp -- "  FAIL ..." per
# failing check and "N checks, M failures" at the end -- because that protocol is
# exactly what the harness parses, and a fixture that spoke a different one would
# leave the parsing untested.
TEST_CPP = """\
#include "sim.hpp"

#include <cmath>
#include <cstdio>

static int checks = 0, failures = 0;

static void near(const char* what, double got, double want, double tolerance) {
    ++checks;
    if (!(std::abs(got - want) <= tolerance)) {
        std::printf("  FAIL %-40s got %+.9g want %+.9g\\n", what, got, want);
        ++failures;
    }
}

int main() {
    std::printf("toy validation\\n");

    // Closed form, not eyeballed output: the relaxation is exp(-rate dt / capacity).
    const double rate = 3.0, capacity = 0.02, dt = 0.5;
    const double want = 1.0 + (0.0 - 1.0) * std::exp(-rate * dt / capacity);
    near("relax reaches the analytic value", sim::relax(0.0, 1.0, rate, capacity, dt), want, 2e-3);

    const double values[4] = {1.0, 2.0, 3.0, 6.0};
    near("mean of four values", sim::mean(values, 4), 3.0, 1e-12);
    near("twice", sim::twice(1.5), 3.0, 1e-12);
    near("checked square root", sim::checked(0.5), std::sqrt(0.5), 1e-12);

    // Enough work that the run has a CPU cost worth bounding, and enough of it in
    // the integrator that a change to the substep controller shows up as one.
    const int repeats = 400000;
    double total = 0.0;
    for (int i = 0; i < repeats; ++i) total += sim::relax(0.0, 1.0, rate, capacity, dt);
    near("the integrator repeats itself", total / repeats, want, 2e-3);

#ifdef FLAKY_PATH
    // A check that disagrees with itself from the third run onwards. This is what a
    // load-sensitive assertion looks like from outside -- `test_jobs.cpp` has one
    // that its own comment says a busy machine may legitimately fail -- and the
    // point of the fixture is that the two baseline runs pass, so the harness only
    // meets it once it is scoring mutants.
    {
        int count = 0;
        if (FILE* handle = std::fopen(FLAKY_PATH, "r")) {
            if (std::fscanf(handle, "%d", &count) != 1) count = 0;
            std::fclose(handle);
        }
        ++count;
        if (FILE* handle = std::fopen(FLAKY_PATH, "w")) {
            std::fprintf(handle, "%d", count);
            std::fclose(handle);
        }
        near("a check that does not always agree with itself",
             (count > 2 && count % 2 == 1) ? 1.0 : 0.0, 0.0, 1e-12);
    }
#endif

    std::printf("%d checks, %d failures\\n", checks, failures);
    return failures == 0 ? 0 : 1;
}
"""

CATALOGUE = '''\
"""Toy catalogue: every substitution here has an answer decided in advance."""
NAME = "selftest"
TREE = {tree!r}
CONFIGURE = []
BUILD = ["g++", "-std=c++20", "-Wall", "-Wextra", "-O1", "-o", "{{build}}/toy",
         "-I", "{{src}}", "{{src}}/sim.cpp", "{{src}}/test.cpp"]
TEST = ["{{build}}/toy"]
COPY_EXCLUDES = ["build*", "__pycache__"]
# Small, because the toy suite is small; the real catalogues take the defaults. The
# wall backstop is the one that costs time -- it is what a *blocked* mutant sits
# against -- so it is 45 s here rather than the ten minutes a real suite wants.
CPU_FLOOR = 6.0
WALL_FLOOR = 45.0
WALL_FACTOR = 6.0

SIM = "sim.cpp"

MUTANTS = [
    # 0 -- an assertion catches it
    ("mean divides by one fewer than it averaged", SIM,
     "    return sum / count;", "    return sum / (count - 1);", "kill"),
    # 1 -- an abort prints no failing check at all
    ("the guard admits a value it exists to reject", SIM,
     "    if (!(x >= 0.0)) std::abort();", "    if (!(x >= 1.0)) std::abort();", "kill"),
    # 2 -- does not compile
    ("the return loses its semicolon", SIM,
     "double twice(double x) {{ return 2.0 * x; }}",
     "double twice(double x) {{ return 2.0 * x }}", "kill"),
    # 3 -- compiles, warns, and is otherwise behaviour-identical. Without "a warning
    #      is a failed build" this reads SURVIVED, which is why it is here.
    ("an unused local, warned about and otherwise harmless", SIM,
     "double twice(double x) {{ return 2.0 * x; }}",
     "double twice(double x) {{ int spare = 3; return 2.0 * x; }}", "kill"),
    # 4 -- spins forever: the sign flip the substep controller turns into a hang
    ("the relaxation drives away from its target", SIM,
     "        const double flux = rate * (target - value);",
     "        const double flux = rate * (value - target);", "hang"),
    # 5 -- blocks forever, burning no CPU. A CPU bound alone cannot see this one.
    ("the integrator waits for something that never arrives", SIM,
     "    double value = start, elapsed = 0.0;",
     "    double value = start, elapsed = 0.0;\\n"
     "    for (;;) std::this_thread::sleep_for(std::chrono::seconds(1));", "hang"),
    # 6 -- CONTROL: genuinely equivalent
    ("CONTROL: twice as a sum rather than a product", SIM,
     "double twice(double x) {{ return 2.0 * x; }}",
     "double twice(double x) {{ return x + x; }}", "survive"),
    # 7 -- CONTROL: equivalent answer, 2.5x the CPU. A wall-clock bound on a busy box
    #      kills this one, and killing a control is the error a kill rate cannot afford.
    ("CONTROL: a finer substep, same answer, 2.5x the work", SIM,
     "        const double step = std::min(dt - elapsed, 0.25 * capacity * capacity /",
     "        const double step = std::min(dt - elapsed, 0.10 * capacity * capacity /",
     "survive"),
    # 8 -- CONTROL: the suite does the same work and spends most of it descheduled,
    #      which is what a busy box does to it. Sleeping and being descheduled are
    #      the same thing seen from here: neither accumulates CPU time.
    ("CONTROL: the same work, mostly not running", SIM,
     "    double value = start, elapsed = 0.0;",
     "    double value = start, elapsed = 0.0;\\n"
     "    if (dt > 1e9) std::this_thread::sleep_for(std::chrono::seconds(20));",
     "survive"),
]
'''

# Mutant 8 above only sleeps when the suite asks it to. In the known-answer sweep
# nothing does, so it is an ordinary equivalent control; for the wall-policy case
# below the toy's own suite is given this extra call, which makes the mutated build
# spend twenty seconds not running while computing exactly the same answers.
SLEEPY_CALL = ("    near(\"the integrator repeats itself\", total / repeats, want, 2e-3);",
               "    near(\"the integrator repeats itself\", total / repeats, want, 2e-3);\n"
               "    (void)sim::relax(0.0, 1.0, rate, capacity, 2e9);")

# A suite that writes a file to its *current directory*. Nothing in the real suite
# does -- all three of its writers go through `testing::scratchDir()` -- but that is
# an audit of today's tests rather than a property of the harness, and this fixture
# is what turns it into one: wherever this file lands is where the mutant binary was
# run from, and it must not be the repository.
STRAY_WRITE = ('    std::printf("toy validation\\n");',
               '    std::printf("toy validation\\n");\n'
               '    if (FILE* stray = std::fopen("stray-from-cwd.txt", "w")) {\n'
               '        std::fprintf(stray, "a relative path, written by the suite\\n");\n'
               '        std::fclose(stray);\n'
               '    }')

def with_tmpdir_marker(test):
    """A suite that writes a *fixed* name under $TMPDIR, as three real ones do.

    `test_shipfile.cpp` writes `barge.ship` and reads it back, `test_breach.cpp`
    writes `ferry_damage_map.txt` and reads it back, `test_collision.cpp` writes
    `collision_ram.csv`. Under one shared /tmp, N workers are N suites over the
    same three files.
    """
    return test.replace(
        "#include <cstdio>", "#include <cstdio>\n#include <cstdlib>", 1).replace(
        '    std::printf("toy validation\\n");',
        '    std::printf("toy validation\\n");\n'
        '    {\n'
        '        const char* dir = std::getenv("TMPDIR");\n'
        '        char marker[512];\n'
        '        std::snprintf(marker, sizeof marker, "%s/shipsim-mutate-worker-marker.txt",\n'
        '                      dir != nullptr ? dir : "/tmp");\n'
        '        if (FILE* handle = std::fopen(marker, "w")) {\n'
        '            std::fprintf(handle, "one worker\\n");\n'
        '            std::fclose(handle);\n'
        '        }\n'
        '    }', 1)


EXPECTED = {
    0: ("KILLED", "failing check"),
    1: ("KILLED", "crash"),
    2: ("KILLED", "compile"),
    3: ("KILLED", "warning"),
    4: ("HANG", "cpu"),
    5: ("HANG", "wall"),
    6: ("SURVIVED", ""),
    7: ("SURVIVED", ""),
    8: ("SURVIVED", ""),
}


def write_toy(root, sim=SIM_CPP, test=TEST_CPP):
    tree = root / "toy"
    tree.mkdir(parents=True, exist_ok=True)
    (tree / "sim.hpp").write_text(SIM_HPP)
    (tree / "sim.cpp").write_text(sim)
    (tree / "test.cpp").write_text(test)
    return tree


def write_catalogue(root, tree, name="toy_catalogue.py", body=None):
    path = root / name
    path.write_text(body if body is not None else CATALOGUE.format(tree=str(tree)))
    return path


def harness(*arguments, timeout=900, cwd=None):
    started = time.monotonic()
    done = subprocess.run([sys.executable, str(MUTATE), *arguments],
                          capture_output=True, text=True, timeout=timeout, cwd=cwd)
    return done, time.monotonic() - started


# --- the cases ---------------------------------------------------------------------

def case_known_answers(root):
    """Every outcome the taxonomy has, against a substitution whose answer is known."""
    tree = write_toy(root)
    catalogue = write_catalogue(root, tree)
    result = root / "known.json"
    done, elapsed = harness("run", str(catalogue), "--json", str(result),
                            "--scratch", str(root))
    print(done.stdout)
    if not check("the sweep ran", result.exists(), done.stderr[-800:]):
        return None
    report = json.loads(result.read_text())
    by_index = {m["index"]: m for m in report["mutants"]}
    for index, (outcome, hint) in EXPECTED.items():
        record = by_index.get(index, {})
        got = record.get("outcome")
        check_equal(f"mutant {index} ({record.get('label', '?')})", got, outcome)
        if hint and got == outcome:
            check(f"mutant {index} is {outcome} for the right reason: {record.get('reason')}",
                  hint in record.get("reason", "").lower(),
                  f"expected the reason to mention {hint!r}")
    # The taxonomy is only useful if the three outcomes are actually distinct: the
    # two hangs must not be counted as kills, and the controls must not be either.
    summary = report["summary"]
    check_equal("real mutants counted", summary["real"], 6)
    check_equal("kills counted", summary["killed"], 4)
    check_equal("hangs counted apart from kills", summary["hung"], 2)
    check_equal("controls counted", summary["controls"], 3)
    check_equal("no control was killed", summary["controls_wrongly_killed"], 0)
    check_equal("the kill rate is the fraction it says it is",
                summary["kill_rate_text"], "4/6")
    check("the tree is byte-identical after the sweep", summary["tree_byte_identical_after"])
    check_equal("nothing was left applied", summary["mutations_applied_after"], [])
    print(f"    (known-answer sweep: {elapsed:.0f}s)")
    return report


def case_hang_is_not_a_kill_or_a_survival(report):
    """The point of a third outcome: the same mutant must not be scoreable either way."""
    if report is None:
        return
    hangs = [m for m in report["mutants"] if m["outcome"] == "HANG"]
    check_equal("both hangs are reported", len(hangs), 2)
    for record in hangs:
        check(f"hang {record['index']} is out of the kill count",
              record["label"] not in report["summary"]["survivor_labels"])
        check(f"hang {record['index']} names what it cost",
              re.search(r"\d+s CPU", record["reason"]) is not None, record["reason"])
    spinner = next((m for m in hangs if "away from its target" in m["label"]), None)
    blocked = next((m for m in hangs if "never arrives" in m["label"]), None)
    if check("the spinning hang was caught by the CPU bound", spinner is not None):
        check("the spinning hang burned CPU", spinner["cpu"] > 1.0, str(spinner["cpu"]))
        check("the CPU bound is what stopped it", "cpu bound" in spinner["reason"],
              spinner["reason"])
    if check("the blocked hang was caught at all", blocked is not None):
        # This is the case the CPU bound cannot see, and the reason a wall backstop
        # exists at all: a deadlocked process accumulates no CPU time, ever.
        check("the blocked hang burned almost no CPU", blocked["cpu"] < 2.0, str(blocked["cpu"]))
        check("the wall backstop is what stopped it", "wall bound" in blocked["reason"],
              blocked["reason"])


def case_wall_policy_would_have_lied(root):
    """The measured claim: a wall-clock bound kills a control that the CPU bound does not.

    Same mutant, same machine, two policies. The mutant is equivalent -- the suite
    does identical work and gets identical answers -- and it spends most of its time
    not running, which is what a shared box does to a process. A wall-clock harness
    calls that a kill and reports strength it does not have.
    """
    tree = write_toy(root, test=TEST_CPP.replace(*SLEEPY_CALL))
    check("the sleepy call landed in the fixture",
          "2e9" in (tree / "test.cpp").read_text())
    catalogue = write_catalogue(root, tree, "sleepy_catalogue.py")
    result = root / "sleepy.json"
    # Mutant 8 only: the control that sleeps.
    done, _ = harness("run", str(catalogue), "--only", "8", "--json", str(result),
                      "--scratch", str(root))
    print(done.stdout)
    if not check("the sleepy sweep ran", result.exists(), done.stderr[-800:]):
        return
    record = json.loads(result.read_text())["mutants"][0]
    check_equal("the CPU bound lets the sleeping control live", record["outcome"], "SURVIVED")
    check("it really did spend its time not running",
          record["wall"] > 2.0 * max(record["cpu"], 0.01),
          f"wall {record['wall']}s vs cpu {record['cpu']}s")
    check("the harness records that a wall-clock policy would have hung it",
          record["wall_policy_would_have_hung"], json.dumps(record))

    # And now actually run that policy, by moving the wall backstop down to where a
    # wall-clock harness would have put it. Same mutant, same machine, opposite
    # verdict -- which is the whole argument for bounding on CPU time.
    other = root / "sleepy-wall.json"
    done, _ = harness("run", str(catalogue), "--only", "8", "--json", str(other),
                      "--scratch", str(root), "--wall-factor", "2", "--cpu-factor", "4",
                      "--wall-floor", "8")
    if check("the wall-policy sweep ran", other.exists(), done.stderr[-800:]):
        wall_record = json.loads(other.read_text())["mutants"][0]
        check_equal("a wall-clock bound at the same factor wrongly hangs the control",
                    wall_record["outcome"], "HANG")
        check("...and reports it as a control that should not have died",
              json.loads(other.read_text())["summary"]["controls_wrongly_killed"] == 1)


def case_load_does_not_manufacture_kills(root, idle):
    """The same control again, under real contention rather than a sleep.

    The contention is made deterministic rather than hoped for: the sweep and five
    busy loops are confined to **one** CPU, which is a busy box in miniature and
    leaves the other twenty-three alone. `idle` is the baseline wall time measured
    earlier on a quiet machine -- what a wall-clock harness would have calibrated
    against before the box got busy, which is exactly the situation that falsely
    killed two controls during the GM sweep.
    """
    tree = write_toy(root)
    catalogue = write_catalogue(root, tree, "loaded_catalogue.py")
    result = root / "loaded.json"
    core = (os.cpu_count() or 2) - 1
    hogs = []
    spin = "x = 0\nwhile True:\n    x += 1\n"
    original = os.sched_getaffinity(0)
    try:
        os.sched_setaffinity(0, {core})  # inherited by everything spawned below
        for _ in range(5):
            hogs.append(subprocess.Popen([sys.executable, "-c", spin],
                                         stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL))
        # Mutant 7: equivalent, 2.5x the CPU. Under load its *wall* time is what a
        # wall-clock harness sees, and that is the number that moves.
        done, _ = harness("run", str(catalogue), "--only", "7", "--json", str(result),
                          "--scratch", str(root))
    finally:
        for hog in hogs:
            hog.kill()
        for hog in hogs:
            hog.wait()
        os.sched_setaffinity(0, original)
    print(done.stdout)
    if not check("the loaded sweep ran", result.exists(), done.stderr[-800:]):
        return
    report = json.loads(result.read_text())
    record = report["mutants"][0]
    check_equal("an equivalent-but-slower control survives under load", record["outcome"],
                "SURVIVED")
    ratio = record["wall"] / max(record["cpu"], 0.01)
    print(f"    (six processes on one core: {record['cpu']:.1f}s CPU, {record['wall']:.1f}s wall, "
          f"wall/cpu {ratio:.2f}, CPU against an idle baseline {record['cpu_ratio']}x, "
          f"wall against an idle baseline {record['wall'] / idle:.2f}x)")
    check("contention really did stretch the wall clock", ratio > 2.0, f"{ratio:.2f}")
    check("the CPU cost of an equivalent mutant stays inside the bound",
          record["cpu_ratio"] is not None and record["cpu_ratio"] < 4.0,
          str(record["cpu_ratio"]))
    check("a wall-clock bound calibrated on the idle machine would have killed this control",
          record["wall"] > 4.0 * idle, f"{record['wall']:.1f}s wall vs 4 x {idle:.1f}s idle")


def case_parallel_and_whole_suite_agree(root):
    """Two options that change *how* a verdict is reached must not change the verdict.

    Running the suite to the end after a failing check, and running mutants in
    parallel copies, are both there for cost rather than for meaning. If either
    moved an outcome, the cheap path would be measuring something else -- which is
    the shape of the mistake `verify.sh` made when it summarised a failing step by
    grepping for the text it expected.
    """
    tree = write_toy(root)
    catalogue = write_catalogue(root, tree, "parallel_catalogue.py")
    full = root / "full.json"
    done, _ = harness("run", str(catalogue), "--only", "0,6", "--json", str(full),
                      "--scratch", str(root), "--no-early-kill", timeout=600)
    if check("the whole-suite sweep ran", full.exists(), done.stderr[-800:]):
        records = {m["index"]: m for m in json.loads(full.read_text())["mutants"]}
        check_equal("a failing check is a kill with the suite run to the end",
                    records[0]["outcome"], "KILLED")
        check("...and the count of failing checks is what says so",
              "failing check(s)" in records[0]["reason"], records[0]["reason"])
        check_equal("the control still survives", records[6]["outcome"], "SURVIVED")

    parallel = root / "parallel.json"
    done, _ = harness("run", str(catalogue), "--only", "0,2,6,7", "--json", str(parallel),
                      "--scratch", str(root), "--workers", "2", timeout=600)
    if check("the parallel sweep ran", parallel.exists(), done.stderr[-800:]):
        report = json.loads(parallel.read_text())
        records = {m["index"]: m["outcome"] for m in report["mutants"]}
        check_equal("two workers reach the same four verdicts", records,
                    {0: "KILLED", 2: "KILLED", 6: "SURVIVED", 7: "SURVIVED"})
        check("the tree is still byte-identical after a parallel sweep",
              report["summary"]["tree_byte_identical_after"])


def case_a_verdict_that_does_not_reproduce(root):
    """A kill that does not happen twice is not a kill.

    This is the GM sweep's failure with the cause put back in: a check in the suite
    that is not a function of the source. The fixture's flaky check passes on both
    baseline runs and then disagrees with itself, so the harness meets it exactly
    where it hurts -- scoring a mutant. The mutant chosen is a **control**, because
    a control dying is the loudest thing this harness can say, and the difference
    between "your equivalence argument is wrong" and "your suite is flaky" is one
    re-run.
    """
    tree = write_toy(root)
    marker = root / "flaky-count"
    body = CATALOGUE.format(tree=str(tree)).replace(
        '"-I", "{src}", "{src}/sim.cpp", "{src}/test.cpp"]',
        f'"-I", "{{src}}", "-DFLAKY_PATH=\\"{marker}\\"",'
        ' "{src}/sim.cpp", "{src}/test.cpp"]')
    check("the flaky fixture was configured", "FLAKY_PATH" in body)
    catalogue = write_catalogue(root, tree, "flaky_catalogue.py", body)
    result = root / "flaky.json"
    done, _ = harness("run", str(catalogue), "--only", "6", "--json", str(result),
                      "--scratch", str(root), "--confirm", "1", timeout=600)
    print(done.stdout)
    if not check("the flaky sweep ran", result.exists(), done.stderr[-800:]):
        return
    report = json.loads(result.read_text())
    record = report["mutants"][0]
    check_equal("a verdict that does not reproduce is FLAKY, not a kill",
                record["outcome"], "FLAKY")
    check_equal("both verdicts are kept", record.get("confirmations"), ["KILLED", "SURVIVED"])
    check_equal("a flaky control is not counted as a control that died",
                report["summary"]["controls_wrongly_killed"], 0)
    check_equal("but the sweep is not clean either", report["summary"]["clean"], False)
    check("...and the flaky mutant is named", report["summary"]["flaky_labels"] != [])

    # And the negative control for the re-run itself: with --confirm 0 the same
    # fixture reports a kill, which is the number a harness without this pass
    # publishes.
    marker.unlink(missing_ok=True)
    naive = root / "naive.json"
    done, _ = harness("run", str(catalogue), "--only", "6", "--json", str(naive),
                      "--scratch", str(root), "--confirm", "0", timeout=600)
    if check("the unconfirmed sweep ran", naive.exists(), done.stderr[-800:]):
        naive_report = json.loads(naive.read_text())
        check_equal("without the re-run the same flaky check reads as a dead control",
                    naive_report["summary"]["controls_wrongly_killed"], 1)


def case_a_sweep_of_nothing_has_no_rate(root):
    """Zero mutants must not come back as a perfect score.

    The figure gate here printed a success line and exited 0 having checked nothing
    at all, because every tool it reads was missing. Zero checks is a failure there
    now, and zero mutants is a failure here: an empty numerator over an empty
    denominator is the most flattering number a harness can print.
    """
    tree = write_toy(root)
    catalogue = write_catalogue(root, tree, "empty_catalogue.py")

    done, _ = harness("run", str(catalogue), "--only", "", "--scratch", str(root), timeout=300)
    check_equal("a selection of no mutants is refused", done.returncode, 2)
    check("...and says why", "no kill rate" in done.stdout, done.stdout[-400:])
    check("...before building anything", "baseline:" not in done.stdout)

    done, _ = harness("run", str(catalogue), "--only", "999", "--scratch", str(root), timeout=300)
    check_equal("an index the catalogue does not have is refused", done.returncode, 2)
    check("...and names it", "999" in done.stdout, done.stdout[-400:])

    # And the subtler one: a sweep of controls *only* runs something, finishes, and
    # has no real mutant in it. It must report no rate rather than 0/0 or 100%.
    result = root / "controls-only.json"
    done, _ = harness("run", str(catalogue), "--only", "6,7", "--json", str(result),
                      "--scratch", str(root), timeout=900)
    print(done.stdout)
    if check("the controls-only sweep ran", result.exists(), done.stderr[-800:]):
        summary = json.loads(result.read_text())["summary"]
        check_equal("no real mutants means no rate", summary["kill_rate"], None)
        check("...and the text says so, rather than a fraction",
              "no real mutants" in summary["kill_rate_text"], summary["kill_rate_text"])
        check_equal("...and the sweep is not clean", summary["clean"], False)
        check("...even though both controls did survive",
              summary["controls"] == 2 and summary["controls_wrongly_killed"] == 0)


def case_the_result_says_how_to_rerun_itself(root):
    """Everything needed to run it again, including the defaults it ran with.

    `zone_probe` printed the parameters it was given and not the ones it used, and
    two published tables lost the defaults they were taken at -- one of them turned
    out to describe a different experiment wearing the same invocation. A kill rate
    has the same exposure.
    """
    tree = write_toy(root)
    catalogue = write_catalogue(root, tree, "rerun_catalogue.py")
    result = root / "rerun.json"
    done, _ = harness("run", str(catalogue), "--only", "0", "--json", str(result),
                      "--scratch", str(root), timeout=900)
    if not check("the sweep ran", result.exists(), done.stderr[-800:]):
        return
    report = json.loads(result.read_text())
    record = report.get("invocation", {})
    for field in ("rerun", "options", "catalogue_file", "catalogue_sha256", "tree",
                  "commands", "mutants_selected", "host"):
        check(f"the result records its {field}", field in record, json.dumps(record)[:400])
    options = record.get("options", {})
    # The defaults are the point: none of these was named on the command line.
    for field, value in (("workers", 1), ("cpu_factor", 4.0), ("wall_factor", 6.0),
                         ("baseline_runs", 2), ("confirm", 1), ("early_kill", True)):
        check_equal(f"the default it actually used for {field} is recorded",
                    options.get(field), value)
    check("the rerun line names the options explicitly rather than relying on defaults",
          all(flag in record.get("rerun", "")
              for flag in ("--workers", "--cpu-factor", "--wall-factor", "--cpu-floor",
                           "--wall-floor", "--baseline-runs", "--confirm", "--only")),
          record.get("rerun", ""))
    check("the catalogue is pinned by content, not just by name",
          record.get("catalogue_sha256") == hashlib.sha256(
              pathlib.Path(catalogue).read_bytes()).hexdigest())
    check_equal("the build and test commands are recorded",
                record.get("commands", {}).get("test"), ["{build}/toy"])

    # The rerun line has to actually run. That is the difference between recording a
    # command and recording one that works.
    again = root / "again.json"
    command = record["rerun"].split() + ["--json", str(again), "--scratch", str(root)]
    repeat = subprocess.run(command, capture_output=True, text=True, timeout=900)
    if check("the recorded rerun line runs", again.exists(), repeat.stderr[-500:]):
        check_equal("...and reaches the same verdict",
                    json.loads(again.read_text())["mutants"][0]["outcome"],
                    report["mutants"][0]["outcome"])


def case_refuses_a_stale_pattern(root):
    """A substitution that no longer matches must be loud at the start, not silent."""
    tree = write_toy(root)
    body = CATALOGUE.format(tree=str(tree)).replace(
        '("mean divides by one fewer than it averaged", SIM,\n'
        '     "    return sum / count;", "    return sum / (count - 1);", "kill"),',
        '("a pattern that no longer occurs", SIM,\n'
        '     "    return sum / countt;", "    return 0.0;", "kill"),')
    check("the stale-pattern fixture was built", "countt" in body)
    catalogue = write_catalogue(root, tree, "stale_catalogue.py", body)
    done, _ = harness("run", str(catalogue), "--scratch", str(root), timeout=300)
    check_equal("a stale pattern is refused", done.returncode, 2)
    check("...and says which one", "a pattern that no longer occurs" in done.stdout,
          done.stdout[-600:])
    check("...before building anything", "baseline:" not in done.stdout)


def case_refuses_a_red_baseline(root):
    """If the suite is red before any mutation, every mutant would score a kill."""
    tree = write_toy(root, test=TEST_CPP.replace("near(\"twice\", sim::twice(1.5), 3.0, 1e-12);",
                                                 "near(\"twice\", sim::twice(1.5), 4.0, 1e-12);"))
    catalogue = write_catalogue(root, tree, "red_catalogue.py")
    done, _ = harness("run", str(catalogue), "--only", "6", "--scratch", str(root), timeout=300)
    check("a red baseline is refused", "BASELINE IS NOT GREEN" in done.stdout + done.stderr,
          (done.stdout + done.stderr)[-600:])
    check_equal("...with a non-zero status", done.returncode, 1)


def case_scan_sees_a_leftover(root):
    """The scanner needs both controls: it must fire, and it must not fire falsely."""
    tree = write_toy(root, )
    catalogue = write_catalogue(root, tree, "scan_catalogue.py")
    done, _ = harness("scan", str(catalogue), "--tree", str(tree), timeout=120)
    check_equal("a pristine tree scans clean", done.returncode, 0)
    check("...and says so", "clean" in done.stdout, done.stdout)

    # Now leave one applied, exactly as a sweep killed mid-iteration would.
    source = tree / "sim.cpp"
    original = source.read_text()
    source.write_text(original.replace("    return sum / count;", "    return sum / (count - 1);"))
    done, _ = harness("scan", str(catalogue), "--tree", str(tree), timeout=120)
    check_equal("an applied mutation is found", done.returncode, 1)
    check("...and named", "APPLIED" in done.stdout and "mean divides" in done.stdout,
          done.stdout)
    done, _ = harness("scan", str(catalogue), "--tree", str(tree), "--revert", timeout=120)
    check("the revert puts the source back byte for byte", source.read_text() == original)
    done, _ = harness("scan", str(catalogue), "--tree", str(tree), timeout=120)
    check_equal("and the tree scans clean again", done.returncode, 0)


def case_interruption_leaves_the_tree_alone(root):
    """SIGKILL mid-sweep. This is the failure that put a mutant in master's history."""
    sys.path.insert(0, str(HERE))
    import mutate as engine

    tree = write_toy(root)
    catalogue_path = write_catalogue(root, tree, "kill_catalogue.py")
    loaded = engine.load_catalogue(str(catalogue_path))
    before = engine.manifest(tree, loaded.copy_excludes)

    process = subprocess.Popen([sys.executable, str(MUTATE), "run", str(catalogue_path),
                                "--scratch", str(root), "--json", str(root / "killed.json")],
                               stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
                               start_new_session=True)
    # Wait until it is demonstrably past the baseline and into applying mutants, so
    # that "nothing changed" is not just "nothing happened yet".
    deadline = time.monotonic() + 300
    saw_a_mutant = False
    while time.monotonic() < deadline:
        line = process.stdout.readline()
        if not line:
            break
        if re.match(r"\s*\d+ (KILLED|SURVIVED|HANG)", line):
            saw_a_mutant = True
            break
    check("the sweep reached its mutants before being killed", saw_a_mutant)
    time.sleep(0.2)
    os.killpg(os.getpgid(process.pid), signal.SIGKILL)
    process.wait()

    after = engine.manifest(tree, loaded.copy_excludes)
    check("SIGKILL mid-sweep leaves the tree byte-identical", before == after)
    applied, stale = engine.scan(loaded, tree)
    check_equal("...and no substitution applied to it", [m.label for m in applied], [])
    check_equal("...and no pattern disturbed", [m.label for m, _ in stale], [])


def case_manifest_notices_a_byte(root):
    """The reporter needs a negative control: a hash that never differs proves nothing."""
    sys.path.insert(0, str(HERE))
    import mutate as engine

    tree = write_toy(root / "manifest")
    before = engine.manifest(tree, [])
    source = tree / "sim.cpp"
    text = source.read_text()
    source.write_text(text.replace("0.25 * capacity", "0.25* capacity", 1))
    after = engine.manifest(tree, [])
    check("one changed byte changes the manifest", after != before)
    # And it must say *which* file. "Something moved" is a failure path that reports
    # only that a failure happened, which is the hole `verify.sh` had twice.
    check_equal("...and names the file that moved", engine.moved(before, after), ["sim.cpp"])
    source.write_text(text)
    check("and putting it back restores the manifest", engine.manifest(tree, []) == before)
    # A file that appears must move it too -- a manifest over content alone would
    # miss a leftover file, which is how 176 junk files once reached master.
    (tree / "stray.txt").write_text("")
    check_equal("a new file is named too", engine.moved(before, engine.manifest(tree, [])),
                ["stray.txt"])
    (tree / "stray.txt").unlink()


def fingerprint(directory):
    """sha256 per file, skipping nothing at all.

    `mutate.manifest` deliberately skips `__pycache__` and `build*` at the top of
    the tree it measures, which is right for a source tree and exactly wrong here:
    the bytecode *is* the thing under test.
    """
    return {str(path.relative_to(directory)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in sorted(directory.rglob("*")) if path.is_file() and not path.is_symlink()}


def case_the_repository_is_never_written_to(root):
    """The docstring's central claim, measured rather than believed.

    It was false in two places when this was first run against the tree, and both
    are the same shape: a write nobody had asked whether it happened.

    A catalogue is a Python module in the repository, so importing one made CPython
    drop a `.pyc` beside it -- on `list`, which builds nothing and runs nothing. The
    sweep's own manifest could not report it, because the catalogue is loaded in
    `main` and the manifest is taken in `sweep`: by the time the tree is first
    hashed the file is already part of it. `.gitignore` names that debris and
    ignores it, which is not the same as it not being there.

    And every build and test inherited the working directory the harness was
    launched from, which is the repository. Today no test writes a relative path --
    all three writers in the suite go through `testing::scratchDir()` -- so nothing
    landed, but that is an audit of the suite and it expires the next time someone
    writes a test. Both halves are checked here, the second with a fixture that
    writes a relative path on purpose.
    """
    # First, the instrument's own negative control: a fingerprint that can never
    # differ proves nothing, and this one is about to be used to prove a negative.
    # In /tmp, because a checker that wrote into the repository to test whether the
    # repository is written to would be reporting on itself.
    probe = root / "probe"
    probe.mkdir(parents=True, exist_ok=True)
    (probe / "a.py").write_text("x = 1\n")
    reference = fingerprint(probe)
    (probe / "__pycache__").mkdir(exist_ok=True)
    (probe / "__pycache__" / "a.pyc").write_bytes(b"\x00")
    check("a stray .pyc changes the fingerprint", fingerprint(probe) != reference)
    check_equal("...and is named", sorted(set(fingerprint(probe)) - set(reference)),
                ["__pycache__/a.pyc"])

    # The half that needs no build: every mode that imports a catalogue, over every
    # catalogue in the tree. Set difference rather than "no .pyc exists anywhere",
    # so debris left by an older harness reads as what it is -- already there --
    # rather than as this run's doing.
    before = fingerprint(HERE)
    catalogues = sorted((HERE / "catalogues").glob("*.py"))
    check("there are catalogues in the tree to import", catalogues != [])
    for path in catalogues:
        harness("list", str(path), timeout=180)
        harness("scan", str(path), timeout=600)
    after = fingerprint(HERE)
    check_equal("importing every catalogue writes no new file into the repository",
                sorted(set(after) - set(before)), [])
    check_equal("...and changes no file already there",
                sorted(n for n in set(before) & set(after) if before[n] != after[n]), [])
    check_equal("...and removes nothing", sorted(set(before) - set(after)), [])

    # The half that needs a sweep: where does a mutant binary's relative write go?
    tree = write_toy(root, test=TEST_CPP.replace(*STRAY_WRITE))
    check("the stray write landed in the fixture",
          "stray-from-cwd" in (tree / "test.cpp").read_text())
    catalogue = write_catalogue(root, tree, "cwd_catalogue.py")
    launch = root / "launched-from"
    scratch = root / "scratch"
    launch.mkdir(parents=True, exist_ok=True)
    scratch.mkdir(parents=True, exist_ok=True)
    done, _ = harness("run", str(catalogue), "--only", "6", "--json", str(root / "cwd.json"),
                      "--scratch", str(scratch), "--keep", timeout=900, cwd=str(launch))
    check_equal("nothing is written to the directory the harness was launched from",
                sorted(q.name for q in launch.iterdir()), [])
    check("the suite's relative write landed inside the copy instead",
          sorted(scratch.rglob("stray-from-cwd.txt")) != [], done.stdout[-500:])
    check("...and the sweep still reported the tree unmoved",
          json.loads((root / "cwd.json").read_text())["summary"]["tree_byte_identical_after"])


def case_workers_do_not_share_a_scratch_directory(root):
    """Parallelism must not be able to decide a verdict.

    On the first `indentation` sweep it did: with four workers on one `/tmp`, four
    copies of the suite wrote `/tmp/barge.ship` and read it back, and mutant 4 was
    scored KILLED on `FAIL a ship written to disk loads back: /tmp/barge.ship:
    empty` -- a failing check about a *neighbour*. That mutant deserved killing
    anyway so the rate did not move, which is precisely why it is worth a check:
    the next one to lose that race could be a control, and a falsely killed control
    inflates the rate while hiding a real gap.

    Deterministic rather than hopeful. Rather than running two workers and hoping
    they collide, this asserts the mechanism that makes a collision impossible:
    each workspace is handed its own directory, the suite's writes land inside the
    sweep's scratch, and the shared name is never created at all.
    """
    tree = write_toy(root, test=with_tmpdir_marker(TEST_CPP))
    check("the marker fixture was configured", "TMPDIR" in (tree / "test.cpp").read_text())
    catalogue = write_catalogue(root, tree, "workers_catalogue.py")
    scratch = root / "scratch"
    scratch.mkdir(parents=True, exist_ok=True)
    shared = pathlib.Path("/tmp/shipsim-mutate-worker-marker.txt")
    shared.unlink(missing_ok=True)

    result = root / "workers.json"
    done, _ = harness("run", str(catalogue), "--only", "0,2,6,7", "--workers", "2",
                      "--json", str(result), "--scratch", str(scratch), "--keep", timeout=900)
    if not check("the two-worker sweep ran", result.exists(), done.stderr[-800:]):
        return
    markers = sorted(scratch.rglob("shipsim-mutate-worker-marker.txt"))
    check("the suite's $TMPDIR write happened at all", markers != [], done.stdout[-500:])
    check("every one of them landed inside the sweep's own scratch",
          all(str(marker).startswith(str(scratch)) for marker in markers),
          str(markers))
    check("and none of them under the shared /tmp the workers would have raced over",
          not shared.exists())
    private = sorted(p for p in scratch.glob("*/tmp*") if p.is_dir())
    check_equal("each of the two workers was given a scratch of its own",
                len(private), 2)
    check_equal("...and the markers sit in those and nowhere else",
                sorted({marker.parent for marker in markers}) == private, True)
    check_equal("the verdicts are the ones the serial sweep reaches",
                {m["index"]: m["outcome"] for m in json.loads(result.read_text())["mutants"]},
                {0: "KILLED", 2: "KILLED", 6: "SURVIVED", 7: "SURVIVED"})


def main():
    started = time.time()
    root = pathlib.Path(tempfile.mkdtemp(prefix="shipsim-mutate-selftest-"))
    print(f"mutation harness self-test, in {root}\n")
    try:
        print("--- every outcome the taxonomy has, against known answers ---")
        report = case_known_answers(root / "known")
        case_hang_is_not_a_kill_or_a_survival(report)
        print("\n--- a wall-clock bound would have killed a control ---")
        case_wall_policy_would_have_lied(root / "sleepy")
        print("\n--- and real contention does not manufacture a kill either ---")
        idle = report["bound"]["baseline_wall"] if report else 1.0
        case_load_does_not_manufacture_kills(root / "loaded", idle)
        print("\n--- a verdict that does not reproduce ---")
        case_a_verdict_that_does_not_reproduce(root / "flaky")
        print("\n--- the cheap paths reach the same verdicts as the expensive ones ---")
        case_parallel_and_whole_suite_agree(root / "parallel")
        print("\n--- a sweep of nothing, and a result that can rerun itself ---")
        case_a_sweep_of_nothing_has_no_rate(root / "empty")
        case_the_result_says_how_to_rerun_itself(root / "rerun")
        print("\n--- what the harness refuses to do ---")
        case_refuses_a_stale_pattern(root / "stale")
        case_refuses_a_red_baseline(root / "red")
        print("\n--- is a mutation applied right now? ---")
        case_scan_sees_a_leftover(root / "scan")
        print("\n--- interruption, and the tree afterwards ---")
        case_interruption_leaves_the_tree_alone(root / "interrupt")
        case_manifest_notices_a_byte(root / "hash")
        print("\n--- and the repository itself, which the docstring claims is never written ---")
        case_the_repository_is_never_written_to(root / "unwritten")
        print("\n--- parallelism must not be able to decide a verdict ---")
        case_workers_do_not_share_a_scratch_directory(root / "workers")
    finally:
        shutil.rmtree(root, ignore_errors=True)
    print(f"\n{checks} checks, {failures} failures ({time.time() - started:.0f}s)")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
