#!/usr/bin/env python3
"""The mutation harness this repository kept rewriting and throwing away.

Around thirty mutation sweeps have run here and roughly seven hundred defects and
test gaps have come out of them, and until now exactly one harness was checked in
(`tools/zone_gpu_probe/mutate.py`, one GPU shader path). Every other sweep was a
throwaway script that was deleted after it ran, which is why an audit of
`docs/07-fem-spike-findings.md` could not re-run a single one of the 52 + 68
mutants that document publishes kill rates for. A kill rate nobody can re-derive
is a number, not a measurement.

    ./tools/mutate/mutate.py list   plasticity
    ./tools/mutate/mutate.py scan   plasticity          # is a mutation applied *now*?
    ./tools/mutate/mutate.py run    plasticity [--workers 4]
    ./tools/mutate/selftest.py                          # the harness against known answers

A catalogue is a small Python module (see `catalogues/`) naming the substitutions,
the build and the test command. Each substitution carries an `expect`: `kill` for a
real mutant, `survive` for a **control** -- an edit argued to be genuinely
equivalent. A harness that only counts kills cannot tell a strong suite from a
harness that is miscounting, and controls are how the difference becomes visible:
a control that dies is a bug in the harness or a flaky test, and it is reported as
loudly as a survivor.

Nothing here is decoration. Each of the following is a mistake this repository has
already paid for, most of them recorded in CLAUDE.md:

**The repository is never written to.** The tree is copied to /tmp and every
mutant is applied, built and run in the copy. Two harnesses here have left a
mutant applied in the working tree after being killed before their cleanup ran --
`les.cpp` substitution 43 sat in a clean-*looking* worktree, and the suite passed
*because the mutant survived*. A `finally` does not survive SIGKILL; a copy has no
such failure mode at all. `run` additionally hashes every file of the real tree
before and after the sweep and refuses to report if a byte moved.

That claim was checked rather than believed, and it was **false in two places**
when this was first run against the tree. Importing a catalogue made CPython write
a `.pyc` beside it -- in the repository, on a `list` that builds nothing -- and the
manifest could not see it because the import happens before the manifest is taken;
`sys.dont_write_bytecode` below is the fix. And every build and test ran with the
repository as its working directory, so one test writing a relative path would have
landed a file there; they now run in the copy. The self-test asserts both, because
a docstring is exactly the kind of claim this repository has learned nothing tests.

**A hang is a distinct outcome from a kill and from a survival.** The
characteristic mutation kill in this codebase is not a failing assertion; it is a
nine-second run becoming an hours-long one, because these are explicit integrators
with substep controllers and a wrong sign collapses the step rather than
perturbing a number. Six such mutants are on record (one in `les.cpp`, five of 196
in `fire.cpp`/`thermal.cpp`), all with zero failing assertions. A harness that only
distinguishes pass from fail scores those six as kills or as survivors depending
on where its timeout happens to land, which is to say it reports a number it did
not measure. HANG is its own outcome here and is excluded from the kill count.

**The bound is CPU time, not wall time.** A wall-clock assertion on a shared box
falsely killed two *controls* during the GM sweep, and a false kill inflates the
rate while hiding a real gap -- the one direction of error a kill rate cannot
afford. The same gate here has gone from 180 s to 358 s purely from sharing the
machine. A descheduled process accumulates no CPU time, so a CPU bound calibrated
on a measured baseline survives a busy box; the wall bound is kept only as a
backstop for a mutant that *blocks* (a deadlock burns no CPU and would otherwise
sit there forever), and is set loose enough that contention cannot reach it. Every
run records both, and the JSON records whether a wall-clock policy at the same
factor *would* have called it -- so the margin is measured rather than assumed.

**A build warning is a failed build.** The build is `-Wall -Wextra -Wpedantic` and
has been warning-clean since Phase 0. A mutant that compiles with a warning is not
a mutant the suite survived.

**A sweep of nothing has no kill rate.** The figure gate here once printed a success
line and exited 0 having checked nothing at all, because every tool it reads was
missing. Zero mutants is refused, an unknown index is refused by name, `--only ''`
selects nothing rather than everything, and a run containing only controls reports
that there is no rate rather than dividing by zero.

**A kill can be right for the wrong reason, and the confirmation pass was blind to
exactly that.** `--confirm` re-runs the verdicts that did *not* match their recorded
expectation, on the argument that a verdict which matched is not the one at risk.
That argument is false in one direction and it has already cost us: on the first
`indentation` sweep, mutant 4 was scored KILLED on `FAIL a ship written to disk
loads back: /tmp/barge.ship: empty, or not a ship file` -- a failing check about a
*neighbouring worker's* file, four suites having shared one `/tmp`. The mutant
deserved killing on its own merits, so the rate did not move and nothing was
re-run, because the outcome matched expectation. A kill for the wrong reason is
invisible to a pass that only looks at surprises. The race itself is closed
(`testing::scratchDir()` now mkdtemps per process), the blindness was not, and the
next cause will be a different one. Three things now look at kills that surprised
nobody:

  * the strings that *decided* a verdict -- the failing check that stopped the
    suite, the warning that failed the build -- are scanned for absolute paths, and
    a kill whose evidence names a file outside the workspace that produced it is
    reported and re-run. `/tmp/barge.ship` is not under any workspace; a
    neighbour's `tmp1/` is not under worker 0's either. This is narrow, and the
    narrowness is honest: it catches a verdict decided by a file the mutant did not
    own, and nothing else.
  * a **sample** of the expected kills is re-run too (`--confirm-kills`, a fraction,
    default 0.1, rounded up so the pass is never entirely blind). A kill that does
    not reproduce is FLAKY and leaves the kill count, exactly as an unexpected one
    does. The sample is chosen by hashing `seed:index`, so it is the same sample on
    every re-run of the same invocation and the `rerun` line carries the seed.
  * a re-run that reaches the same verdict is also asked whether it reached it the
    *same way*. The `/tmp/barge.ship` kill would have reproduced as a kill --
    mutant 4 was a real mutant -- while its reason changed completely. Outcome
    equality is a weaker check than it looks, and comparing reasons costs nothing
    once the re-run has happened.

What none of this can see is a kill that is wrong *and* deterministic *and* names
no path: an assertion that fails for a reason unrelated to the mutation, every
time. That needs the third option -- proving the failing check is one the mutated
code can reach -- and the catalogue does not record which *symbol* it edited, only
a file and a fragment of text, so there is nothing to start a call graph from.

**The result says how to run itself again, defaults included.** `zone_probe` printed
the parameters it was *given* rather than the ones it *used*, and two published
tables lost the defaults they were taken at -- one of them turned out to describe a
different experiment wearing the same invocation. The JSON therefore carries a
`rerun` line naming every option explicitly, the catalogue's sha256, the tree and
its HEAD, and the commands that built and ran it.
"""
import argparse
import concurrent.futures
import dataclasses
import hashlib
import importlib.util
import json
import os
import pathlib
import re
import shutil
import signal
import subprocess
import sys
import tempfile
import threading
import time

# **The one writer that got past the claim above.** A catalogue is a Python module
# loaded from the tree, and importing it makes CPython drop a `__pycache__` beside
# it -- so `mutate.py list plasticity`, which builds nothing and runs nothing,
# left `tools/mutate/catalogues/__pycache__/plasticity.cpython-314.pyc` in the
# repository. The manifest cannot see it: `load_catalogue` runs in `main`, before
# `sweep` takes the "before" hash, so the file is already there when the tree is
# first measured and the sweep reports byte-identical. `.gitignore` catches it at
# the commit and says in as many words that this is the same class of debris as
# the 176 files a `git add -A` once swept in -- which makes it ignored, not
# absent. This is the fix, and the self-test drives `list` and `scan` over every
# catalogue afterwards and asserts nothing appeared.
sys.dont_write_bytecode = True

HERE = pathlib.Path(__file__).resolve().parent
REPO = HERE.parents[1]
CLOCK_TICKS = os.sysconf("SC_CLK_TCK")

KILLED, SURVIVED, HANG = "KILLED", "SURVIVED", "HANG"
ERROR, FLAKY = "ERROR", "FLAKY"
EXPECTED_OUTCOME = {"kill": KILLED, "survive": SURVIVED, "hang": HANG}

# Defaults a catalogue may override. `{src}` and `{build}` are substituted.
DEFAULTS = {
    "CONFIGURE": ["cmake", "-S", "{src}", "-B", "{build}", "-G", "Ninja",
                  "-DCMAKE_BUILD_TYPE=RelWithDebInfo"],
    "BUILD": ["ninja", "-C", "{build}", "shipsim_tests"],
    "TEST": ["{build}/shipsim_tests"],
    "COPY_EXCLUDES": [".git", "build", "build*", ".claude", "*.o", "*.spv", "__pycache__"],
    "WARNING": r"warning:",
    # The line the suite prints for a failing check, and the line it ends with.
    "FAIL_LINE": r"^\s*FAIL\b",
    "SUMMARY": r"(\d+) checks, (\d+) failures",
    # Multipliers on the *measured* baseline. See the module docstring: the CPU one
    # is the bound, the wall one is a backstop against a mutant that blocks.
    "CPU_FACTOR": 4.0,
    "WALL_FACTOR": 20.0,
    "CPU_FLOOR": 30.0,
    "WALL_FLOOR": 300.0,
    # **Every worker gets its own scratch, and these are the variables pointed at
    # it.** `tests/harness.cpp` resolves `SHIPSIM_TEST_TMPDIR`, then `TMPDIR`, then
    # `/tmp`, and three tests write *fixed* names under it -- `barge.ship`,
    # `ferry_damage_map.txt`, `collision_ram.csv` -- and read two of them back. Four
    # workers sharing one `/tmp` therefore run four suites over the same three
    # files, and on the first `indentation` sweep one of them lost: mutant 4 was
    # scored KILLED on `FAIL a ship written to disk loads back: /tmp/barge.ship:
    # empty, or not a ship file`, which is a verdict about a *neighbour* and not
    # about the mutation. That one deserved killing on its own merits, so the rate
    # did not move -- but the same race landing on a control is exactly the GM
    # sweep's two false kills, and a false kill inflates the rate while hiding a
    # real gap. Parallelism must not be able to decide a verdict.
    "SCRATCH_ENV": ["TMPDIR", "TMP", "TEMP", "SHIPSIM_TEST_TMPDIR"],
    # Absolute paths a verdict may legitimately name without owning them: the
    # read-only roots a toolchain diagnostic points into. `/usr/include/c++/14/bits/
    # stl_algo.h:1234: warning:` is a real warning about the mutant and must not be
    # mistaken for a neighbour's scratch file. Everything else absolute and outside
    # the workspace is foreign, `/tmp` emphatically included -- that is where the
    # incident happened.
    "FOREIGN_PATH_EXEMPT": ["/usr/", "/lib/", "/lib64/", "/opt/", "/bin/", "/sbin/",
                            "/etc/", "/proc/", "/sys/", "/dev/", "/nix/store/",
                            "/gnu/store/", "/run/"],
}


# A path with at least two components, not preceded by a word character, another
# slash or a colon -- so `https://host/thing` does not read as `//host/thing`, and a
# bare `/tmp` mentioned in prose does not read as a file. Trailing punctuation is
# left out of the character class on purpose: the incident's line was
# `... loads back: /tmp/barge.ship: empty`, and the path there ends at the colon.
PATH_IN_TEXT = re.compile(r"(?<![\w/:.~-])(/[\w.+@%-]+(?:/[\w.+@%-]+)+)")


# --- running something under a bound -------------------------------------------

@dataclasses.dataclass
class Run:
    """One execution of a build or a test, with what it cost and how it ended."""
    code: int = 0            # negative means killed by that signal
    wall: float = 0.0
    cpu: float = 0.0
    bound: str = ""          # "", "cpu" or "wall": which bound stopped it
    stopped_early: str = ""  # the line that made us stop it (a failing check)
    lines: list = dataclasses.field(default_factory=list)  # tail, for the report
    kept: list = dataclasses.field(default_factory=list)   # every line matching `keep`

    @property
    def hung(self):
        return self.bound != ""


def process_cpu(pid):
    """User+system CPU seconds of a process *and its reaped children*.

    /proc is the only place this is available without a third-party package, and it
    is the whole reason the bound can be CPU rather than wall: a process that is
    not running does not accumulate any. `cutime`/`cstime` matter because the build
    step is `ninja` spawning compilers -- without them a build's CPU reads as zero.
    """
    try:
        with open(f"/proc/{pid}/stat", "rb") as handle:
            text = handle.read().decode("utf-8", "replace")
    except OSError:
        return None  # exited between the poll and the read
    # The command name is in parentheses and may itself contain spaces, so fields
    # are counted from the last ')'. utime stime cutime cstime are proc(5) fields
    # 14-17, and field 3 is the first after the name.
    try:
        fields = text[text.rindex(")") + 2:].split()
        return sum(int(fields[i]) for i in (11, 12, 13, 14)) / CLOCK_TICKS
    except (ValueError, IndexError):
        return None


def supervise(command, cpu_limit, wall_limit, stop_pattern=None, keep_pattern=None, tail=40,
              env=None, cwd=None):
    """Run `command`, bounded by CPU seconds first and wall seconds as a backstop.

    `cwd` is what makes "the repository is never written to" structural rather than
    a property of every test in the suite. Without it the mutant binary inherits
    the directory `mutate.py` was launched from, which is the repository, and one
    test writing a *relative* path would put a file there. Three tests here write
    files and all three go through `testing::scratchDir()`, so today nothing lands
    -- but that is an audit of the suite rather than a property of the harness, and
    an audit is only true until the next test is written.

    `stop_pattern`, when given, stops the run at the first matching output line. A
    suite that has already printed a failing check is already a kill and the
    remaining hundred seconds buy nothing -- this is what makes a sweep affordable.
    It cannot change a verdict: the line it stops on *is* the verdict.

    `keep_pattern` lines are kept whatever else scrolls past. Only the *tail* of the
    output is retained otherwise, and a compiler warning in the middle of a
    hundred-line build would scroll off it -- which would quietly turn "warnings are
    failures" into "warnings are failures if they happen last".
    """
    started = time.monotonic()
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                               text=True, bufsize=1, start_new_session=True, env=env, cwd=cwd)
    run = Run()
    matcher = re.compile(stop_pattern) if stop_pattern else None
    keeper = re.compile(keep_pattern) if keep_pattern else None
    seen = []
    hit = {"line": ""}

    def reader():
        for line in process.stdout:
            line = line.rstrip("\n")
            seen.append(line)
            del seen[:-tail]
            run.lines = list(seen)
            if keeper is not None and len(run.kept) < 50 and keeper.search(line):
                run.kept.append(line)
            if matcher and not hit["line"] and matcher.search(line):
                hit["line"] = line

    pump = threading.Thread(target=reader, daemon=True)
    pump.start()

    stopped = False

    def stop(reason):
        nonlocal stopped
        stopped = True
        run.bound = reason
        try:
            os.killpg(os.getpgid(process.pid), signal.SIGKILL)
        except (ProcessLookupError, PermissionError):
            pass

    status, usage = None, None
    while True:
        done, raw, rusage = os.wait4(process.pid, os.WNOHANG)
        if done == process.pid:
            status, usage = raw, rusage
            break
        cpu = process_cpu(process.pid)
        if cpu is not None:
            run.cpu = cpu
        run.wall = time.monotonic() - started
        # Once a decision has been taken the bounds stop being consulted, so a slow
        # death cannot relabel a kill as a hang.
        if not stopped:
            if hit["line"]:
                run.stopped_early = hit["line"]
                stop("")  # not a bound: a verdict already printed by the suite
            elif run.cpu > cpu_limit:
                stop("cpu")
            elif run.wall > wall_limit:
                stop("wall")
        time.sleep(0.05)

    run.wall = time.monotonic() - started
    # rusage is reported for a killed child too, and is exact where the /proc poll
    # is only as fine as its interval.
    if usage is not None:
        run.cpu = max(run.cpu, usage.ru_utime + usage.ru_stime)
    try:
        process.returncode = os.waitstatus_to_exitcode(status)
    except ValueError:
        process.returncode = -1
    run.code = process.returncode
    pump.join(timeout=5)
    try:
        process.stdout.close()
    except OSError:
        pass
    if hit["line"] and not run.stopped_early:
        run.stopped_early = hit["line"]
    return run


# --- catalogues ------------------------------------------------------------------

@dataclasses.dataclass
class Mutant:
    index: int
    label: str
    path: str
    old: str
    new: str
    expect: str
    note: str = ""


class Catalogue:
    def __init__(self, module, source):
        self.module = module
        self.source = source
        self.name = getattr(module, "NAME", pathlib.Path(source).stem)
        for key, value in DEFAULTS.items():
            setattr(self, key.lower(), getattr(module, key, value))
        self.tree = pathlib.Path(getattr(module, "TREE", REPO))
        self.mutants = []
        for index, entry in enumerate(module.MUTANTS):
            label, path, old, new, expect = entry[:5]
            # `hang` is an expectation and not only an outcome, because six mutants
            # in this repo are on record as killing the suite by hanging. Writing
            # that down is what makes the next sweep's verdict comparable to this
            # one's instead of depending on where a timeout landed.
            if expect not in ("kill", "survive", "hang"):
                raise SystemExit(f"mutant {index}: expect must be kill, survive or hang,"
                                 f" not {expect!r}")
            # **Every substitution must leave something to find.** A mutant whose
            # replacement is empty is indistinguishable, afterwards, from a
            # catalogue pattern that has drifted out of the source -- both look like
            # "the original text is not there". That ambiguity would land in exactly
            # the check this harness exists to make trustworthy: is a mutation
            # applied right now? Write the deletion as an edit that leaves a mark:
            # `if (x == 0.0)` -> `if (x < 0.0)`, or a loop bound of 0.
            if new == "" or new == old:
                raise SystemExit(f"mutant {index} ({label}): a substitution must replace the"
                                 f" text with something detectable, not with nothing")
            note = entry[5] if len(entry) > 5 else ""
            self.mutants.append(Mutant(index, label, path, old, new, expect, note))
        if not self.mutants:
            raise SystemExit(f"{source}: no mutants")

    def command(self, key, src, build):
        return [part.format(src=src, build=build) for part in getattr(self, key)]


def load_catalogue(name):
    path = pathlib.Path(name)
    if not path.exists():
        path = HERE / "catalogues" / (name if name.endswith(".py") else name + ".py")
    if not path.exists():
        available = sorted(p.stem for p in (HERE / "catalogues").glob("*.py"))
        raise SystemExit(f"no catalogue {name!r}; have: {', '.join(available)}")
    spec = importlib.util.spec_from_file_location("catalogue_" + path.stem, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return Catalogue(module, str(path))


# --- "is a mutation applied right now?" -------------------------------------------

def scan(catalogue, tree, revert=False):
    """First-class answer to the question a killed sweep leaves behind.

    An interrupted sweep leaves a mutation in the source, and the suite then passes
    *because the mutant survived* -- the one situation where a green suite is
    evidence against you. It has happened here. This re-derives every substitution
    in the catalogue and looks for it in the tree, which is the check CLAUDE.md
    says to keep next to the harness rather than in your head.

    Returns (applied, stale). `applied` is the loud one.
    """
    applied, stale = [], []
    cache = {}
    for mutant in catalogue.mutants:
        path = tree / mutant.path
        if path not in cache:
            try:
                cache[path] = path.read_text()
            except OSError as error:
                stale.append((mutant, f"unreadable: {error}"))
                cache[path] = ""
        text = cache[path]
        old, new = text.count(mutant.old), text.count(mutant.new)
        # **Presence of the substituted text, not absence of the original.** Some
        # substitutions *insert* -- their replacement contains the line they
        # replace -- so "the original is gone" scores those as clean however
        # thoroughly they are applied. The self-test carries two such mutants
        # because this harness got it wrong the same way its predecessors did.
        inserts = mutant.old in mutant.new
        if new >= 1 and (inserts or old == 0):
            applied.append(mutant)
            if revert:
                text = text.replace(mutant.new, mutant.old, 1)
                path.write_text(text)
                cache[path] = text
        elif new >= 1:
            # The substituted text already occurs in a tree that has the original in
            # it too, so this substitution can never be detected. That is a defect in
            # the catalogue and it has to be said, not averaged away.
            stale.append((mutant, f"the substituted text already occurs {new} times in the"
                                  f" pristine source: this mutant cannot be detected"))
        elif old != 1:
            # Not necessarily an applied mutant -- more often a catalogue that has
            # drifted from the code it describes. Either way it means this mutant
            # cannot be trusted to answer the question, which must be said out loud
            # rather than counted as "clean".
            stale.append((mutant, f"the original text occurs {old} times, the mutated one {new}"))
    return applied, stale


def manifest(tree, excludes):
    """sha256 per file under `tree`, so "the tree did not move" is checkable.

    Per file rather than one digest over everything: "the tree moved" on its own is
    a failure path that reports only that something happened, which is the shape of
    reporting `verify.sh` had to have knocked out of it twice. Naming the files is
    the difference between "do not trust this result" and "the files under mutation
    are untouched and the operator edited a document".
    """
    skip = set(excludes) | {".git", ".claude"}
    entries = {}
    for path in sorted(tree.rglob("*")):
        top = path.relative_to(tree).parts[0]
        if top in skip or top.startswith("build") or top == "__pycache__":
            continue
        if not path.is_file() or path.is_symlink():
            continue
        entries[str(path.relative_to(tree))] = hashlib.sha256(path.read_bytes()).hexdigest()
    return entries


def moved(before, after):
    """Which files differ between two manifests, as a sorted list of names."""
    return sorted(set(before) ^ set(after)
                  | {name for name in set(before) & set(after) if before[name] != after[name]})


# --- was this kill about the mutation, or about something else? -------------------

def verdict_evidence(catalogue, reason, build, test):
    """The untruncated strings a verdict was actually read from.

    Not the whole output. A path scan over a hundred lines of build chatter is a
    false-positive generator, and `tests/harness.cpp` alone prints `test output kept
    in /tmp/...` on every run without that line deciding anything. Only lines that
    *are* a verdict qualify: the failing check that stopped the suite, the warnings
    that failed the build, and -- when the suite was run to the end and the verdict
    is a count of failures -- the failing checks themselves. The reason string is
    included too, but it is truncated to 90 characters and a path can be cut in
    half by that, which is why the lines it was built from are here as well.
    """
    items = [reason]
    if build is not None:
        items += list(build.kept)
    if test is not None:
        if test.stopped_early:
            items.append(test.stopped_early)
        failing = re.compile(catalogue.fail_line)
        items += [line for line in test.lines if failing.search(line)]
    return items


def paths_outside(items, owned, exempt):
    """Absolute paths in `items` that no directory in `owned` contains.

    `owned` is *this workspace's* three directories rather than the sweep's scratch
    root, deliberately: worker 0 being killed by a file under worker 1's `tmp1/` is
    the original incident with the shared `/tmp` swapped for a shared root, and it
    would read as "inside the sweep" to a check that only looked at the root.
    """
    owned = [str(root).rstrip("/") for root in owned]
    found = []
    for text in items:
        for match in PATH_IN_TEXT.finditer(text or ""):
            path = match.group(1)
            # Containment in both directions. The forward test is the question being
            # asked. The backward one is about truncation: the reason string is cut
            # at 90 characters, so a workspace path can arrive as
            # `/tmp/shipsim-mutate-ab/tm` -- which no longer *starts with* any owned
            # directory and would otherwise be accused of being somebody else's. A
            # string that could be a prefix of a directory this workspace owns is not
            # evidence that it is not, and a false accusation here costs more than a
            # missed one: it is the same asymmetry as a false kill against a survivor.
            if any(path == root or path.startswith(root + "/") or root.startswith(path)
                   for root in owned):
                continue
            if any(path.startswith(prefix) for prefix in exempt):
                continue
            if path not in found:
                found.append(path)
    # The reason string is truncated to 90 characters and a path can be cut in half
    # by that, so the same file can arrive twice -- once whole from the line the
    # verdict was read from, once as its own prefix. Report the whole one. A
    # directory named alongside a file inside it collapses the same way, which is
    # what we want in a report of what a verdict touched.
    return [path for path in found
            if not any(other != path and other.startswith(path) for other in found)]


# Timings and counts move between two runs of the same verdict for reasons that are
# not the verdict: `no failing check and no exit after 31s CPU / 44s wall` is the
# same reason twice. Digits are the only part that legitimately drifts, so they are
# the only part normalised away -- a *different failing check* still reads as
# different, which is the whole point.
DIGITS = re.compile(r"\d+(?:\.\d+)?")


def same_reason(first, second):
    return DIGITS.sub("#", first) == DIGITS.sub("#", second)


def confirmation_sample(indices, fraction, seed):
    """Which of the kills that surprised nobody get re-run anyway.

    Deterministic, and deterministic in a way that does not depend on the Python
    build: `random.sample` is not specified to give the same answer across versions,
    and a sample nobody can reproduce would put the `rerun` line back in the state
    `zone_probe` left its published tables in. Ordering by sha256 of `seed:index` is
    fully specified by this line.

    Rounded up rather than down, so a fraction above zero always samples at least
    one kill. Rounding down would make a small catalogue silently fall back to the
    behaviour this exists to replace.
    """
    if not indices or fraction <= 0:
        return []
    count = min(len(indices), -(-int(round(fraction * len(indices) * 1000)) // 1000))
    count = max(1, count)
    order = sorted(indices,
                   key=lambda i: hashlib.sha256(f"{seed}:{i}".encode()).hexdigest())
    return sorted(order[:count])


# --- the sweep ---------------------------------------------------------------------

@dataclasses.dataclass
class Bound:
    cpu: float
    wall: float
    build_cpu: float
    build_wall: float
    baseline_cpu: float
    baseline_wall: float
    baseline_spread: float


class Workspace:
    """A copy of the tree with its own build directory. Nothing here is the repo."""

    def __init__(self, catalogue, root, tag):
        self.catalogue = catalogue
        self.src = root / f"src{tag}"
        self.build = root / f"build{tag}"
        self.tmp = root / f"tmp{tag}"
        shutil.copytree(catalogue.tree, self.src,
                        ignore=shutil.ignore_patterns(*catalogue.copy_excludes))
        self.build.mkdir(parents=True, exist_ok=True)
        self.tmp.mkdir(parents=True, exist_ok=True)
        self.originals = {}

    def environment(self):
        """This workspace's own scratch, so two workers cannot write one file.

        See `SCRATCH_ENV`. Started from the real environment rather than replacing
        it, because a build needs `PATH` and a Vulkan test needs its ICD.
        """
        return os.environ | {name: str(self.tmp) for name in self.catalogue.scratch_env}

    def configure(self, bound):
        command = self.catalogue.command("configure", self.src, self.build)
        if not command:
            return Run()
        return supervise(command, bound.cpu, bound.wall, cwd=self.src,
                         env=self.environment())

    def build_once(self, bound):
        return supervise(self.catalogue.command("build", self.src, self.build),
                         bound.build_cpu, bound.build_wall,
                         keep_pattern=self.catalogue.warning, cwd=self.src,
                         env=self.environment())

    def test_once(self, bound, early=True):
        return supervise(self.catalogue.command("test", self.src, self.build),
                         bound.cpu, bound.wall,
                         stop_pattern=self.catalogue.fail_line if early else None,
                         keep_pattern=self.catalogue.summary, cwd=self.src,
                         env=self.environment())

    def apply(self, mutant):
        path = self.src / mutant.path
        if path not in self.originals:
            self.originals[path] = path.read_text()
        text = self.originals[path]
        if text.count(mutant.old) != 1:
            raise RuntimeError(f"mutant {mutant.index}: pattern occurs {text.count(mutant.old)}x")
        path.write_text(text.replace(mutant.old, mutant.new, 1))
        # Read it back and confirm. CLAUDE.md: after a scripted edit, grep to confirm
        # it landed -- a non-matching search-and-replace once silently dropped a test
        # assertion here and the failure looked like a bug in working code.
        after = path.read_text()
        if after == text or mutant.new not in after:
            raise RuntimeError(f"mutant {mutant.index}: the edit did not land")

    def restore(self):
        for path, text in self.originals.items():
            if path.read_text() != text:
                path.write_text(text)


def judge(catalogue, build, test):
    """Turn a build and a test run into one of KILLED / SURVIVED / HANG."""
    if build.hung:
        return HANG, f"the build hit its {build.bound} bound", build
    if build.code != 0:
        return KILLED, f"did not compile (exit {build.code})", build
    if build.kept:
        # Warnings are failures in this build, so a mutant that only warns has not
        # been survived by anything. `kept` rather than the tail: see `supervise`.
        return KILLED, f"built with a warning: {build.kept[0].strip()[:90]}", build
    if test is None:
        return KILLED, "not run", build
    if test.hung:
        return HANG, (f"no failing check and no exit after {test.cpu:.0f}s CPU / "
                      f"{test.wall:.0f}s wall ({test.bound} bound)"), test
    if test.stopped_early:
        return KILLED, "failing check: " + test.stopped_early.strip()[:90], test
    summary = None
    for line in test.kept or test.lines:
        found = re.search(catalogue.summary, line)
        if found:
            summary = found
    if summary is None:
        # A crash prints no failing check at all. An earlier harness here counted
        # FAIL lines only and scored eight aborted binaries as survivors.
        signalled = f" (signal {-test.code})" if test.code < 0 else ""
        return KILLED, f"no summary line: crashed at exit {test.code}{signalled}", test
    failures = int(summary.group(2))
    if failures > 0:
        return KILLED, f"{failures} failing check(s)", test
    if test.code != 0:
        return KILLED, f"suite green but exited {test.code}", test
    return SURVIVED, f"{summary.group(1)} checks, none failing", test


def measure_baseline(catalogue, workspace, args):
    """Calibrate the bound off this machine, now, rather than off a constant.

    The suite is run twice. The second run is not redundant: the difference between
    the two is this box's own variance, and a bound has to clear that before it can
    claim to be measuring a mutant. It is also the check that the baseline is green
    -- if it is not, every mutant would score a kill and the sweep is worthless.
    """
    generous = Bound(cpu=3600, wall=7200, build_cpu=3600, build_wall=7200,
                     baseline_cpu=0, baseline_wall=0, baseline_spread=0)
    configure = workspace.configure(generous)
    if configure.code != 0:
        raise SystemExit("configure failed in the copy:\n  " + "\n  ".join(configure.lines[-15:]))
    build = workspace.build_once(generous)
    if build.code != 0:
        raise SystemExit("baseline build failed in the copy:\n  " + "\n  ".join(build.lines[-15:]))
    # `kept`, not the tail: a warning in the middle of a hundred-line build scrolls
    # off the tail, and the baseline is the one build where a warning means the whole
    # sweep is measuring a tree that was already degraded.
    if build.kept:
        raise SystemExit("the baseline build warns, and warnings are failures here:\n  "
                         + build.kept[0])
    runs = []
    for _ in range(args.baseline_runs):
        run = workspace.test_once(generous, early=False)
        verdict, reason, _ = judge(catalogue, build, run)
        if verdict != SURVIVED:
            raise SystemExit(f"BASELINE IS NOT GREEN ({reason}) -- every mutant would score a kill")
        runs.append(run)
    cpu = max(r.cpu for r in runs)
    wall = max(r.wall for r in runs)
    spread = (cpu - min(r.cpu for r in runs)) / cpu if cpu > 0 else 0.0
    return Bound(cpu=max(catalogue.cpu_floor, args.cpu_factor * cpu),
                 wall=max(catalogue.wall_floor, args.wall_factor * wall),
                 build_cpu=max(catalogue.cpu_floor, args.cpu_factor * max(build.cpu, 1.0)),
                 build_wall=max(catalogue.wall_floor, args.wall_factor * max(build.wall, 1.0)),
                 baseline_cpu=cpu, baseline_wall=wall, baseline_spread=spread)


def invocation(catalogue, args, chosen):
    """Everything needed to run this sweep again, defaults included.

    `zone_probe` learned this the expensive way: it printed the parameters it was
    *given* and not the ones it *used*, so two published tables lost the defaults
    they were taken at, and one of them turned out to describe a different
    experiment wearing the same invocation. A kill rate has the same exposure --
    which tree, which catalogue, which bound, how many workers, whether the suite
    was stopped at its first failing check -- so the JSON carries the lot, and a
    `rerun` line that states every option explicitly rather than relying on what
    today's defaults happen to be.
    """
    options = {"only": "(all)" if args.only is None else args.only, "workers": args.workers,
               "cpu_factor": args.cpu_factor, "wall_factor": args.wall_factor,
               "cpu_floor": catalogue.cpu_floor, "wall_floor": catalogue.wall_floor,
               "baseline_runs": args.baseline_runs, "confirm": args.confirm,
               # The sample is only reproducible if the seed travels with it, and an
               # unreproducible sample is the state `zone_probe` left its tables in.
               "confirm_kills": args.confirm_kills, "confirm_seed": args.confirm_seed,
               "early_kill": not args.no_early_kill,
               "scratch_env": catalogue.scratch_env,
               "foreign_path_exempt": catalogue.foreign_path_exempt}
    rerun = [sys.executable, str(pathlib.Path(__file__).resolve()), "run", args.catalogue,
             "--workers", str(args.workers), "--cpu-factor", repr(args.cpu_factor),
             "--wall-factor", repr(args.wall_factor), "--cpu-floor", repr(catalogue.cpu_floor),
             "--wall-floor", repr(catalogue.wall_floor),
             "--baseline-runs", str(args.baseline_runs), "--confirm", str(args.confirm),
             "--confirm-kills", repr(args.confirm_kills),
             "--confirm-seed", str(args.confirm_seed)]
    if args.only is not None:
        rerun += ["--only", args.only]
    if args.no_early_kill:
        rerun += ["--no-early-kill"]
    head, dirty = "", None
    try:
        head = subprocess.run(["git", "-C", str(catalogue.tree), "rev-parse", "HEAD"],
                              capture_output=True, text=True, timeout=30).stdout.strip()
        status = subprocess.run(["git", "-C", str(catalogue.tree), "status", "--porcelain"],
                                capture_output=True, text=True, timeout=60)
        dirty = bool(status.stdout.strip())
    except (OSError, subprocess.SubprocessError):
        pass
    return {
        "argv": sys.argv,
        "rerun": " ".join(rerun),
        "options": options,
        "catalogue_file": catalogue.source,
        "catalogue_sha256": hashlib.sha256(
            pathlib.Path(catalogue.source).read_bytes()).hexdigest(),
        "mutants_selected": [m.index for m in chosen],
        "mutant_count": len(chosen),
        "tree": str(catalogue.tree),
        "tree_head": head,
        "tree_had_uncommitted_changes": dirty,
        "commands": {"configure": catalogue.configure, "build": catalogue.build,
                     "test": catalogue.test},
        "host": {"cpus": os.cpu_count(), "loadavg": [round(x, 2) for x in os.getloadavg()]},
    }


def sweep(catalogue, args):
    chosen = catalogue.mutants
    # `--only ''` is a selection of nothing, not a selection of everything. A script
    # that computes an index list -- "re-run the survivors" -- and computes an empty
    # one would otherwise sweep the whole catalogue while believing it had run the
    # three mutants it named.
    if args.only is not None:
        wanted = {int(x) for x in args.only.replace(" ", "").split(",") if x != ""}
        missing = sorted(wanted - {m.index for m in catalogue.mutants})
        if missing:
            print(f"--only names {missing}, which this catalogue does not have")
            return 2
        chosen = [m for m in catalogue.mutants if m.index in wanted]
    # **A sweep that runs nothing must not report a rate.** The figure gate here
    # once printed a success line and exited 0 having checked nothing at all,
    # because every tool it reads was missing; zero checks is now a failure there
    # and zero mutants is a failure here. An empty numerator over an empty
    # denominator is the most flattering number a harness can print.
    if not chosen:
        print("REFUSING TO RUN: no mutants selected, and a sweep of nothing has no kill rate")
        return 2

    tree = catalogue.tree
    before = manifest(tree, catalogue.copy_excludes)
    applied, stale = scan(catalogue, tree)
    if applied:
        print("REFUSING TO START: a mutation from this catalogue is already applied to " + str(tree))
        for mutant in applied:
            print(f"  {mutant.index:3d} {mutant.label}")
        return 3
    if stale:
        # Every pattern is checked before any is applied, so a catalogue that has
        # drifted from the code is a loud failure at second one rather than a silent
        # survivor two hours in.
        print("catalogue patterns that do not apply cleanly to the tree:")
        for mutant, why in stale:
            print(f"  {mutant.index:3d} {mutant.label}: {why}")
        print("every substitution must match the pristine source exactly once")
        return 2

    root = pathlib.Path(tempfile.mkdtemp(prefix="shipsim-mutate-", dir=args.scratch))
    print(f"catalogue {catalogue.name}: {len(chosen)} of {len(catalogue.mutants)} mutants")
    print(f"working in {root} (the repository is never written to)")
    started = time.time()
    results = []
    try:
        workspaces = []
        for worker in range(args.workers):
            workspaces.append(Workspace(catalogue, root, f"{worker}"))
        bound = measure_baseline(catalogue, workspaces[0], args)
        print(f"baseline: {bound.baseline_cpu:.0f}s CPU, {bound.baseline_wall:.0f}s wall, "
              f"{bound.baseline_spread * 100:.0f}% spread over {args.baseline_runs} runs")
        print(f"bound:    {bound.cpu:.0f}s CPU per mutant "
              f"(x{args.cpu_factor:g}), wall backstop {bound.wall:.0f}s\n")
        for workspace in workspaces[1:]:
            if workspace.configure(bound).code != 0 or workspace.build_once(bound).code != 0:
                raise SystemExit("a worker's copy failed to build")

        free = list(workspaces)
        lock = threading.Lock()
        printed = threading.Lock()

        def one(mutant):
            with lock:
                workspace = free.pop()
            build = test = None
            try:
                began = time.time()
                workspace.apply(mutant)
                build = workspace.build_once(bound)
                if build.code == 0 and not build.hung:
                    test = workspace.test_once(bound, early=not args.no_early_kill)
                verdict, reason, run = judge(catalogue, build, test)
            except Exception as error:  # noqa: BLE001 - one bad mutant must not cost the sweep
                # An hour of results is not worth losing to a single unusable entry,
                # and swallowing it silently would be worse: ERROR is neither a kill
                # nor a survival, and it keeps the sweep from being reported clean.
                verdict, reason, run = "ERROR", f"{type(error).__name__}: {error}", Run()
                build = build or Run()
            finally:
                workspace.restore()
                with lock:
                    free.append(workspace)
            record = {
                "index": mutant.index, "label": mutant.label, "file": mutant.path,
                "expect": mutant.expect, "outcome": verdict, "reason": reason,
                "note": mutant.note,
                "phase": "build" if run is build else "test",
                "cpu": round(run.cpu, 1), "wall": round(run.wall, 1),
                "build_cpu": round(build.cpu, 1),
                # Against the baseline *suite*, so it only means anything for a run
                # that got as far as the suite.
                "cpu_ratio": (round(test.cpu / bound.baseline_cpu, 2)
                              if test is not None and bound.baseline_cpu else None),
                "elapsed": round(time.time() - began, 1),
                # What a wall-clock bound at the same factor would have decided. The
                # GM sweep's two falsely killed controls are what this column exists
                # to make visible, and it costs nothing to record.
                "wall_policy_would_have_hung":
                    bool(test is not None and not test.hung
                         and test.wall > args.cpu_factor * bound.baseline_wall),
                # Absolute paths in the evidence that this workspace does not own.
                # Recorded for every verdict, acted on for kills: a survivor whose
                # output mentions a neighbour's file is a curiosity, a *kill*
                # decided by one is a verdict about somebody else's mutation.
                "foreign_paths": paths_outside(
                    verdict_evidence(catalogue, reason, build, test),
                    [workspace.src, workspace.build, workspace.tmp],
                    catalogue.foreign_path_exempt),
            }
            wanted = EXPECTED_OUTCOME[mutant.expect]
            flag = ""
            if verdict != wanted:
                flag = "   <-- UNEXPECTED" + (" CONTROL DEATH" if mutant.expect == "survive" else "")
            elif verdict == HANG:
                flag = "   <-- hung, as recorded; neither killed nor survived"
            with printed:
                print(f"{mutant.index:3d} {verdict:8s} ({record['elapsed']:6.1f}s, "
                      f"{record['cpu']:6.1f}s CPU) {mutant.label}: {reason}{flag}")
                sys.stdout.flush()
            return record

        if args.workers == 1:
            results = [one(m) for m in chosen]
        else:
            with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
                results = list(pool.map(one, chosen))

        # **Run the surprises again.** A suite is not a deterministic function of the
        # source: this one contains at least one check that is documented as
        # load-sensitive (`test_jobs.cpp` only asks that stealing happened at all,
        # and says in its own comment that a heavily loaded machine may legitimately
        # run everything on one lane). Under a three-worker sweep it fired, and it
        # killed a *control* -- an edit argued equivalent. Counted naively that is a
        # false kill, which inflates the rate and hides a real gap.
        #
        # Only the surprises are re-run, because they are few and because a verdict
        # that matches its recorded expectation is not the one at risk of being an
        # accident. A verdict that does not reproduce is reported FLAKY and kept out
        # of the kill count entirely. The re-runs are serial even when the sweep was
        # not: the thing being re-examined may well be load-sensitive, so this is the
        # wrong moment to add load.
        #
        # **And the surprises are not the only verdicts at risk.** The argument above
        # is sound in one direction only. A kill that *matched* its expectation is
        # never re-examined, so a kill that happened for the wrong reason cannot be
        # seen at all -- which is not hypothetical: `indentation` mutant 4 was scored
        # KILLED on a failing check about a neighbouring worker's `/tmp/barge.ship`.
        # It deserved killing anyway, the rate did not move, and nothing looked. Two
        # more populations are therefore re-run: the kills whose evidence names a
        # file outside the workspace that produced them, and a sample of the rest.
        for record in results:
            record["kill_names_a_foreign_path"] = bool(record["outcome"] == KILLED
                                                       and record["foreign_paths"])
        # Nothing is *marked* for re-checking unless it is going to be re-checked.
        # `--confirm 0` switches the whole pass off, and a report that named a sample
        # it never ran would be claiming an examination that did not happen -- the
        # same shape as the figure gate printing a success line having checked
        # nothing.
        settled, sample, recheck = [], [], []
        if args.confirm > 0:
            for record in results:
                if record["outcome"] != EXPECTED_OUTCOME[record["expect"]]:
                    record["rechecked_because"] = "the verdict did not match its expectation"
                elif record["kill_names_a_foreign_path"]:
                    record["rechecked_because"] = ("the failing check names "
                                                   + ", ".join(record["foreign_paths"][:3])
                                                   + ", which this workspace does not own")
            settled = [r["index"] for r in results
                       if r["expect"] == "kill" and r["outcome"] == KILLED
                       and "rechecked_because" not in r]
            sample = confirmation_sample(settled, args.confirm_kills, args.confirm_seed)
            for record in results:
                if record["index"] in sample:
                    record["rechecked_because"] = "sampled from the kills that surprised nobody"
            recheck = [r for r in results if "rechecked_because" in r]
        surprises = [r for r in results if r["outcome"] != EXPECTED_OUTCOME[r["expect"]]]
        if recheck:
            print(f"\nre-running {len(recheck)} verdict(s) {args.confirm}x to see whether they"
                  f" reproduce: {len(surprises)} unexpected,"
                  f" {sum(1 for r in recheck if r['kill_names_a_foreign_path'])} naming a path"
                  f" outside their workspace, {len(sample)} sampled from {len(settled)}"
                  f" kill(s) that matched expectation")
            index = {m.index: m for m in chosen}
            for record in recheck:
                seen = [record["outcome"]]
                reasons = [record["reason"]]
                for _ in range(args.confirm):
                    again = one(index[record["index"]])
                    seen.append(again["outcome"])
                    reasons.append(again["reason"])
                record["confirmations"] = seen
                record["confirmation_reasons"] = reasons
                if len(set(seen)) > 1:
                    record["outcome"] = FLAKY
                    record["reason"] = ("does not reproduce: " + ", ".join(seen)
                                        + " -- " + record["reason"])
                    print(f"{record['index']:3d} FLAKY    {record['label']}: {record['reason']}")
                    continue
                # **The same verdict is not the same measurement.** Mutant 4 above
                # would have reproduced as a kill -- it was a real mutant, and the
                # second run would have killed it on its own merits while the
                # neighbour's file was gone. Outcome equality would have called that
                # confirmed. Only the digits are normalised away, because a hang's
                # reason carries its seconds; a *different failing check* still reads
                # as different, which is the entire point.
                record["reason_reproduced"] = all(same_reason(reasons[0], other)
                                                  for other in reasons[1:])
                if not record["reason_reproduced"]:
                    print(f"{record['index']:3d} SAME VERDICT, ANOTHER REASON  {record['label']}:")
                    for line in reasons:
                        print(f"      {line}")
    finally:
        if args.keep:
            print(f"\nscratch kept at {root}")
        else:
            shutil.rmtree(root, ignore_errors=True)

    changed = moved(before, manifest(tree, catalogue.copy_excludes))
    applied_after, _ = scan(catalogue, tree)
    report = summarise(catalogue, results, bound, changed, applied_after,
                       time.time() - started, args)
    report["invocation"] = invocation(catalogue, args, chosen)
    path = pathlib.Path(args.json or (pathlib.Path(tempfile.gettempdir()) /
                                      f"shipsim-mutate-{catalogue.name}.json"))
    path.write_text(json.dumps(report, indent=2) + "\n")
    print(f"\nmachine-readable result: {path}")
    return 0 if report["summary"]["clean"] else 1


def summarise(catalogue, results, bound, changed, applied_after, elapsed, args):
    real = [r for r in results if r["expect"] in ("kill", "hang")]
    controls = [r for r in results if r["expect"] == "survive"]
    killed = [r for r in real if r["outcome"] == KILLED]
    survivors = [r for r in real if r["outcome"] == SURVIVED]
    hangs = [r for r in results if r["outcome"] == HANG]
    errors = [r for r in results if r["outcome"] == ERROR]
    flaky = [r for r in results if r["outcome"] == FLAKY]
    dead_controls = [r for r in controls if r["outcome"] not in (SURVIVED, FLAKY)]
    # These three are the confirmation pass's new eyes. `kill_names_a_foreign_path`
    # is recorded before any re-run can change an outcome, so a foreign-path kill
    # that then turned FLAKY is still counted here as what it was.
    foreign = [r for r in results if r.get("kill_names_a_foreign_path")]
    lying = [r for r in results if r.get("reason_reproduced") is False]
    sampled = [r for r in results
               if r.get("rechecked_because", "").startswith("sampled")]
    expected_kills = [r for r in real if r["expect"] == "kill"
                      and r["outcome"] in (KILLED, FLAKY)]
    rechecked_kills = [r for r in expected_kills if "confirmations" in r]
    unreproduced = [r for r in rechecked_kills if r["outcome"] == FLAKY]
    # Hangs are reported apart from the kill rate rather than folded into it, and
    # *both* readings are printed. A sweep that hides the distinction is reporting a
    # number that depends on where its timeout landed; a sweep that publishes only
    # the flattering one is rounding in its author's favour.
    # **No real mutants means no rate, not a perfect one.** `killed/real` with an
    # empty denominator is the one arithmetic a kill rate must never perform: a
    # sweep that ran only controls, or that was filtered down to nothing, would
    # otherwise publish 0/0 as though it had measured something.
    rate = len(killed) / len(real) if real else None
    with_hangs = (len(killed) + len(hangs)) / len(real) if real else None
    print()
    if not real:
        print("NO REAL MUTANTS RAN -- there is no kill rate to report, only "
              f"{len(controls)} control(s)")
    else:
        print(f"{len(killed)} of {len(real)} real mutants killed ({rate * 100:.1f}%), "
              f"{len(hangs)} hung, {len(controls) - len(dead_controls)} of {len(controls)} "
              f"controls survived as they should")
    if hangs and real:
        print(f"  counting a hang as a kill, which is a choice and not a measurement, "
              f"that is {len(killed) + len(hangs)}/{len(real)} ({with_hangs * 100:.1f}%)")
    for record in survivors:
        print(f"  SURVIVOR {record['index']}: {record['label']}")
    for record in hangs:
        print(f"  HANG     {record['index']}: {record['label']} -- {record['reason']}")
    for record in dead_controls:
        print(f"  CONTROL WRONGLY {record['outcome']} {record['index']}: {record['label']}"
              f" -- {record['reason']}")
    for record in errors:
        print(f"  ERROR    {record['index']}: {record['label']} -- {record['reason']}")
    for record in flaky:
        print(f"  FLAKY    {record['index']}: {record['label']} -- {record['reason']}")
    for record in foreign:
        print(f"  KILLED BY A PATH IT DOES NOT OWN {record['index']}: {record['label']}"
              f" -- {', '.join(record['foreign_paths'])}")
    if foreign:
        print("  a verdict decided by a file outside the workspace is a verdict about"
              " something other than the mutation")
    for record in lying:
        print(f"  SAME VERDICT, ANOTHER REASON {record['index']}: {record['label']}"
              f" -- {' | '.join(record['confirmation_reasons'])}")
    if rechecked_kills:
        print(f"  {len(rechecked_kills)} of {len(expected_kills)} expected kill(s) were re-run"
              f" ({len(sampled)} of them sampled); {len(unreproduced)} did not reproduce."
              f" The rate above is not corrected for the"
              f" {len(expected_kills) - len(rechecked_kills)} that were not re-run")
    elif expected_kills:
        print("  no expected kill was re-run: --confirm-kills 0 leaves a kill for the wrong"
              " reason invisible")
    if changed:
        mutated = {m["file"] for m in results}
        print(f"  THE TREE MOVED DURING THE SWEEP: {len(changed)} file(s) differ from how they"
              f" started")
        for name in changed[:10]:
            print(f"    {name}{'   <-- A FILE UNDER MUTATION' if name in mutated else ''}")
        if any(name in mutated for name in changed):
            print("  a file under mutation moved: this result is not trustworthy")
    for mutant in applied_after:
        print(f"  A MUTATION IS APPLIED TO THE TREE: {mutant.index} {mutant.label}")
    ratios = [r["cpu_ratio"] for r in results if r["outcome"] != HANG and r["cpu_ratio"]]
    near = max(ratios) if ratios else 0.0
    false_kills = [r for r in results if r["wall_policy_would_have_hung"]]
    print(f"  the slowest mutant that did not hang cost {near:.2f}x the baseline's CPU, "
          f"against a bound of {args.cpu_factor:g}x")
    if false_kills:
        print(f"  a wall-clock bound at the same factor would have wrongly hung "
              f"{len(false_kills)}: " + ", ".join(str(r["index"]) for r in false_kills))
    print(f"  {elapsed / 60:.1f} minutes, {args.workers} worker(s)")
    return {
        "catalogue": catalogue.name,
        "source": catalogue.source,
        "tree": str(catalogue.tree),
        "when": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "bound": dataclasses.asdict(bound) | {"cpu_factor": args.cpu_factor,
                                              "wall_factor": args.wall_factor},
        "summary": {
            "real": len(real), "killed": len(killed), "survived": len(survivors),
            "hung": len(hangs), "controls": len(controls),
            "controls_wrongly_killed": len(dead_controls),
            "kill_rate": round(rate, 4) if rate is not None else None,
            "kill_rate_text": (f"{len(killed)}/{len(real)}" if real
                               else "no real mutants ran: there is no rate"),
            "kill_rate_counting_hangs": round(with_hangs, 4) if with_hangs is not None else None,
            "kill_rate_counting_hangs_text": (f"{len(killed) + len(hangs)}/{len(real)}" if real
                                              else "no real mutants ran: there is no rate"),
            "errors": len(errors),
            "flaky": len(flaky),
            "flaky_labels": [r["label"] for r in flaky],
            "unexpected": [r["index"] for r in results
                           if r["outcome"] != EXPECTED_OUTCOME[r["expect"]]],
            "survivor_labels": [r["label"] for r in survivors],
            "hang_labels": [r["label"] for r in hangs],
            "control_death_labels": [r["label"] for r in dead_controls],
            "worst_cpu_ratio_of_a_finishing_mutant": near,
            "wall_policy_false_hangs": [r["index"] for r in false_kills],
            "kills_naming_a_path_outside_the_workspace": [r["index"] for r in foreign],
            "foreign_paths_named": sorted({q for r in foreign for q in r["foreign_paths"]}),
            "kills_whose_reason_did_not_reproduce": [r["index"] for r in lying],
            "expected_kills": len(expected_kills),
            "expected_kills_rechecked": [r["index"] for r in rechecked_kills],
            "expected_kills_sampled": [r["index"] for r in sampled],
            "expected_kills_that_did_not_reproduce": [r["index"] for r in unreproduced],
            "tree_byte_identical_after": not changed,
            "tree_files_that_moved": changed,
            "a_mutated_file_moved": any(m["file"] in set(changed) for m in results),
            "mutations_applied_after": [m.label for m in applied_after],
            "minutes": round(elapsed / 60, 2),
            "workers": args.workers,
            # A sweep is only clean if nothing surprised us: no dead control, no
            # mutant the harness could not run, no mutation left behind, no tree
            # movement. Survivors are a finding, not a harness fault, so they do not
            # clear this flag.
            # A sweep with nothing real in it is never clean, whatever else went
            # right: there is no measurement to be clean about.
            # A kill decided by a file the workspace does not own, and a kill that
            # reproduced as a kill by a different route, are both results this sweep
            # cannot stand behind -- whatever the rate came out at.
            "clean": (bool(real) and not dead_controls and not errors and not flaky
                      and not applied_after and not changed and not foreign and not lying),
        },
        "mutants": results,
    }


def main(argv=None):
    # A sweep is an hour long and is watched through a redirect, where Python's
    # stdout is block-buffered: without this the log stays empty until the run ends,
    # which is indistinguishable from a sweep that has hung -- in a tool whose whole
    # job is telling a hang from a result.
    sys.stdout.reconfigure(line_buffering=True)
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("mode", choices=["run", "list", "scan"])
    parser.add_argument("catalogue")
    parser.add_argument("--only", default=None,
                        help="comma-separated mutant indices; an empty list selects nothing")
    parser.add_argument("--workers", type=int, default=1,
                        help="parallel copies; sound only because the bound is CPU time")
    parser.add_argument("--json", default="", help="where to write the result (default /tmp)")
    parser.add_argument("--scratch", default=None, help="where the copies go (default /tmp)")
    parser.add_argument("--keep", action="store_true", help="leave the copies behind")
    parser.add_argument("--cpu-factor", type=float, default=None)
    parser.add_argument("--wall-factor", type=float, default=None)
    parser.add_argument("--cpu-floor", type=float, default=None)
    parser.add_argument("--wall-floor", type=float, default=None)
    parser.add_argument("--baseline-runs", type=int, default=2)
    parser.add_argument("--confirm", type=int, default=1,
                        help="re-run each unexpected verdict this many times; a verdict that"
                             " does not reproduce is reported FLAKY rather than counted")
    parser.add_argument("--confirm-kills", type=float, default=0.1,
                        help="also re-run this fraction of the kills that matched their"
                             " expectation, rounded up; 0 restores the behaviour that"
                             " cannot see a kill for the wrong reason")
    parser.add_argument("--confirm-seed", type=int, default=None,
                        help="which sample; default is derived from the catalogue's sha256"
                             " and is recorded in the result")
    parser.add_argument("--no-early-kill", action="store_true",
                        help="run every suite to the end even after a failing check")
    parser.add_argument("--tree", default="", help="which tree to scan (default: the repo)")
    parser.add_argument("--revert", action="store_true", help="scan: undo what it finds")
    args = parser.parse_args(argv)

    catalogue = load_catalogue(args.catalogue)
    if args.cpu_factor is None:
        args.cpu_factor = catalogue.cpu_factor
    if args.wall_factor is None:
        args.wall_factor = catalogue.wall_factor
    if args.cpu_floor is not None:
        catalogue.cpu_floor = args.cpu_floor
    if args.wall_floor is not None:
        catalogue.wall_floor = args.wall_floor
    if not 0.0 <= args.confirm_kills <= 1.0:
        raise SystemExit(f"--confirm-kills is a fraction of the expected kills, "
                         f"between 0 and 1, not {args.confirm_kills!r}")
    if args.confirm_seed is None:
        # Derived rather than random, so two runs of the same invocation re-run the
        # same mutants; different catalogues still get different samples.
        # `--confirm-seed` is how you ask for a *different* tenth of one catalogue.
        args.confirm_seed = int(hashlib.sha256(
            pathlib.Path(catalogue.source).read_bytes()).hexdigest()[:8], 16)

    if args.mode == "list":
        for mutant in catalogue.mutants:
            print(f"{mutant.index:3d} [{mutant.expect:7s}] {mutant.path}: {mutant.label}")
        return 0

    if args.mode == "scan":
        tree = pathlib.Path(args.tree) if args.tree else catalogue.tree
        applied, stale = scan(catalogue, tree, revert=args.revert)
        for mutant, why in stale:
            print(f"STALE   {mutant.index:3d} {mutant.label}: {why}")
        for mutant in applied:
            print(f"APPLIED {mutant.index:3d} {mutant.label}"
                  + ("  (reverted)" if args.revert else ""))
        if not applied and not stale:
            print(f"clean: none of the {len(catalogue.mutants)} substitutions is applied to {tree}")
        return 1 if applied and not args.revert else (2 if stale else 0)

    return sweep(catalogue, args)


if __name__ == "__main__":
    sys.exit(main())
