#!/usr/bin/env bash
# Every number the documents publish about a tool, checked against a re-run of that
# tool. It started on docs/06-roadmap.md's Phase 3 milestone and now covers the FEM
# findings, the simulation and renderer docs, and README.md.
#
# **Why this exists.** The repo's rule is that docs must not drift from the code,
# and nothing enforced it. Two drifts had accumulated by the time anyone re-ran
# the tool:
#
#   - `indentation.hpp` said its own span fix was "left for a session that can
#     re-validate them" while `indentation.cpp` had applied it for weeks, and two
#     other documents repeated the claim because each quoted the previous one
#     rather than the measurement.
#   - the milestone's flooding figures were written before a wing tank authored
#     *inside* a hold was fixed. That removed 217 m3 which had been flooding
#     twice, and moved every published tonnage -- 8100 t at the quarter became
#     7102 t -- with nobody re-running them.
#
# A green test suite says nothing about either, because **nothing tests a
# comment**. This does, and it is why it belongs in `verify.sh full` rather than
# in a notebook somewhere.
#
# **What it learned about itself.** The first version checked figures alone, and
# the amidships strike turned out to sit about 4% above a capsize threshold: at
# 5.5 m/s she lolls, at 5.75 she goes over, and the tonnage jumps 7415 -> 16425 t
# across that gap. A tonnage on its own is nearly meaningless there, so the two
# speeds are checked as a *pair* and what is asserted is which side of 90 degrees
# she ends on. A change that moved the threshold outside the bracket would leave
# both tonnages looking individually plausible and still have altered what the
# ship does.
#
# **And what it counted about itself.** The count this script prints is not a count
# of figures. At 155 it was 77 figures and 78 checks of its own pointers -- useful,
# and not the same thing -- and neither number says *where* the coverage points.
# Counted by document, which is the only way the README hole was ever going to be
# visible, and against the crudest honest denominator there is, a line of the
# document carrying a digit:
#
#   document                            lines  with a digit   figures gated
#   README.md                             321            73        35
#   docs/01-architecture.md               540            98         0
#   docs/02-simulation.md                6151          1686         9
#   docs/03-renderer-audio.md            1179           217         6  ->  56
#   docs/04-multiplayer.md                130            18         0
#   docs/05-data-modding-validation.md    423            72         0
#   docs/06-roadmap.md                   1392           451        22
#   docs/07-fem-spike-findings.md        1661           585         5
#
# 77 figures before, 127 after, inside a printed count that went 155 -> 229. Read
# the two columns apart: a rising total is not coverage arriving anywhere in
# particular, which is the whole of what the README hole taught.
#
# Two things fall out of reading it. **`06-roadmap.md` is the default sixth argument
# to `check`, so counting explicit document arguments scores it zero** -- 22 of its
# figures are gated by calls that name no document at all, and a coverage audit that
# counts arguments rather than destinations gets the roadmap exactly backwards.
#
# And **`03-renderer-audio.md` was the worst covered of the documents that carry
# tool output at all**: 217 lines with a digit against six figures, and the six were
# the headline numbers of three findings rather than any of the measurements the
# findings rest on. That was worth more than the raw ratio suggests, because it is
# the only document describing a subsystem whose figures cannot be re-derived by
# reading the source -- they are pixels and vertex counts. Two of them had drifted:
# `ram_view`'s scene triangle count, and a 10 MW counterfactual that no longer holds
# either in its temperature or in its conclusion. The three blocks below now cover
# it, and the tools they run -- `seaway_view` most of all -- had published figures
# that nothing in the repository re-ran.
#
# **`01` no longer stays at zero, and the reason it used to is a lesson about
# document-level reasoning.** The rule -- `01` and `04` are design documents whose
# digits are dates, section numbers and budgets nothing computes -- is true of most
# of `01` and false of the table at §104-109 and the paragraphs under it, which are
# `job_bench` output and nothing else. Judging a document as a whole hid a table of
# live measurements inside a correct generalisation about the document containing
# it; CI has been running `job_bench` since the beginning under a comment saying its
# figures are what the docs quote, checking only the exit code. Four of those
# figures had drifted by the time anyone re-derived them. `04` and `05` do stay at
# zero: `04`'s digits really are budgets and dates, and `05`'s are the fields of a
# file format `test_shipfile.cpp` already parses. A figure is worth gating here when
# a *tool* produces it -- which is a question about the line, not about the file.
#
# Cost: the damage figures need no flooding at all and run with `--duration=1` in
# about a second each; three runs pay for the full 900 s. Around 105 s for the
# `ram_view` half, and another ~230 s for the README block below it -- the front
# page quotes three 1800 s scenario runs and the size of the test suite, and none of
# those can be had without running the thing. The renderer blocks add ~18 s: ~5 s of
# `seaway_view`, ~12 s of `smoke_view` and ~1 s of a sixth `ram_view` run.
set -u

RAM=${RAM:-./build/ram_view}
WATER=${WATER:-./build/water_probe}
GAS=${GAS:-./build/gas_probe}
GASDOC=engine/sim/promotion.hpp
ZONE=${ZONE:-./build/zone_probe}
ZONE_DOC=docs/02-simulation.md
JOBS=${JOBS:-./build/job_bench}
ARCH_DOC=docs/01-architecture.md
FLIP=${FLIP:-./build/flip_probe}
FLIP_DOC=engine/sim/flip.hpp
ESC_DOC=docs/flip-escalation-design.md
SECTION=${SECTION:-./build/section_probe}
SECTION_DOC=docs/02-simulation.md
# The same file under the name the roll-damping block reads it by: those figures
# are §2 seakeeping and have nothing to do with the section mesher.
SIM_DOC=docs/02-simulation.md
RENDER_DOC=docs/03-renderer-audio.md
DOC=docs/06-roadmap.md
fails=0
checks=0

# **The binaries this reads must not change while it reads them**, and until a
# figure check reported three section-mesher figures as drifted, nothing said so.
# What had happened was a `ninja` finishing mid-run: the check started at
# 02:16:06, `section_probe` was relinked at 02:17:44 from a tree carrying an
# experimental edit, and the section block sits ninety-odd checks in. The figures
# were a correct measurement of a *different program*, and the obvious reading was
# a defect in a mesher nobody had touched for sixty-five commits.
#
# One `stat` at each end turns that into a line of output. It is deliberately not
# a hash: the point is to catch a rebuild, and mtime is what a rebuild changes.
# The tool paths are resolved further down, so this stats `build/` rather than
# naming them: the first version listed `$SHIPSIM` and friends here and died on
# `unbound variable` under `set -u`, three hundred lines before those are set.
# Globbing the directory is also stricter -- it notices a tool that appears or
# disappears mid-run, which naming a fixed list cannot.
toolStamps() {
  local out="" t
  for t in ./build/*; do
    [ -x "$t" ] && [ -f "$t" ] || continue
    out="$out $t:$(stat -c '%Y:%s' "$t" 2>/dev/null || echo '?')"
  done
  printf '%s' "$out"
}
stampsBefore=$(toolStamps)
# **Every block that does not run is named in the summary.** See the note above the
# summary itself: a gate that skipped everything used to print `ok` and a count of
# zero, which `verify.sh` accepts, and the count is the only thing that would have
# said so.
skipped=""

red=$'\033[0;31m'; green=$'\033[0;32m'; dim=$'\033[2m'; off=$'\033[0m'

# check <label> <expected> <tolerance> <actual> <doc quote to find>
#
# **The tolerance is absolute, always** -- `awk` on |actual - expected|. This
# comment used to say it was "relative for magnitudes", which it has never been:
# a reader sizing a tolerance for a figure of 356736 from that sentence would
# have written 0.01 meaning 1% and got a gate a hundred-thousand times tighter
# than intended. It is deliberately tight: these
# runs are deterministic, so a figure that moves at all is either a real change in
# the physics or a real change in the ship, and both want a human to look.
check() {
  local label="$1" expect="$2" tol="$3" actual="$4" quote="$5" where="${6:-$DOC}"
  checks=$((checks + 1))
  # Validate the pointer on every call, deduplicated. See the note above `hint`.
  #
  # **The key is the document as well as the quote**, and it was the quote alone.
  # Nothing had collided yet -- all 57 pointers were distinct strings -- but two
  # documents quoting the same table row is the ordinary case here, not an exotic
  # one, and the second of them would have been recorded as validated on the
  # strength of the first. A dedup key that is narrower than the thing it stands
  # for is a checker that reports coverage it does not have.
  case "$checkedPointers" in
    *"|$where|$quote|"*) ;;
    *) checkedPointers="$checkedPointers|$where|$quote|"
       checks=$((checks + 1))
       if ! grep -qF -- "$quote" "$where"; then
         printf '  %s✗%s the pointer %s"%s"%s no longer occurs in %s\n' \
                "$red" "$off" "$dim" "$quote" "$off" "$where"
         fails=$((fails + 1))
       fi ;;
  esac
  if [ -z "$actual" ]; then
    printf '  %s✗%s %s — could not parse the figure out of the tool output\n' "$red" "$off" "$label"
    fails=$((fails + 1)); return 1
  fi
  # **An empty parse was rejected and a garbled one was not.** `awk` coerces any
  # non-numeric string to 0, so every check whose expected value *is* zero -- and
  # there are six of them -- passed on a tool whose output format had moved out
  # from under the `sed`. `n/a`, `-`, a stray unit, an error message: all zero,
  # all green. That is the same shape as the reporters `verify.sh selftest`
  # exists for, one level down: a failure path that only recognises the failure
  # it was written for.
  if ! awk -v a="$actual" \
       'BEGIN { exit !(a ~ /^[+-]?([0-9]+\.?[0-9]*|\.[0-9]+)([eE][+-]?[0-9]+)?$/) }'; then
    printf '  %s✗%s %s — the figure parsed out of the tool is not a number: %s%s%s\n' \
           "$red" "$off" "$label" "$dim" "$actual" "$off"
    fails=$((fails + 1)); return 1
  fi
  if awk -v a="$actual" -v e="$expect" -v t="$tol" \
        'BEGIN { d = a - e; if (d < 0) d = -d; exit !(d <= t) }'; then
    printf '  %s✓%s %s = %s\n' "$green" "$off" "$label" "$actual"
    return 0
  fi
  printf '  %s✗%s %s: %s publishes %s, the tool now says %s\n' \
         "$red" "$off" "$label" "$where" "$expect" "$actual"
  printf '      %sthe doc line to update contains: %s%s\n' "$dim" "$quote" "$off"
  fails=$((fails + 1))
  return 1
}

# **Every pointer is checked, including the ones only a failure would print.**
# `hint` below validates the pointers it is given, and that was only half the job:
# `check`'s fifth argument is a pointer too, and it is *only* read when a figure has
# already drifted -- so it rots unobserved and misleads exactly when it is needed.
# One had: `"46 of 49"` occurred nowhere in the document it named, in the script
# written to prevent that. Found by an agent doing something else. Now every `check`
# validates its pointer up front, so a rotted one fails on a green run rather than
# waiting for a red one.
checkedPointers=""

# A pointer into the doc is itself a claim about the doc, and it rots the same
# way. The first draft here pointed at "the hole grows from 3.4 to", which never
# matched anything because the paragraph wraps between "grows" and "from" -- so
# every failure message would have sent the reader to a line that does not exist.
# Each hint is checked against the doc before any figure is.
hint() {
  local where="${2:-$DOC}"
  checks=$((checks + 1))
  if grep -qF -- "$1" "$where"; then return 0; fi
  printf '  %s✗%s the pointer %s"%s"%s no longer occurs in %s\n' \
         "$red" "$off" "$dim" "$1" "$off" "$where"
  fails=$((fails + 1))
  return 1
}

# A hook for `verify.sh selftest`, which drives every reporter in this repository
# with a failure it was **not** written for. `check` is a reporter like any other
# and had no control on it at all: it rejected an empty parse and accepted a
# garbled one, because `awk` turns any non-numeric string into 0 and six checks
# here expect exactly 0.
#
# Driven through the real `check` rather than a mock, so the control breaks if the
# predicate is ever loosened. Placed after the definitions and before any tool
# runs; it costs nothing on a normal invocation and exits without touching a
# binary.
if [ "${1:-}" = --selftest-parse ]; then
  check "control: the tool printed something that is not a figure" 0 0 \
        "${CTL_FIGURE:-n/a}" "shipsim" README.md
  if [ "$fails" -gt 0 ]; then exit 1; fi
  echo "control did not go red: a non-numeric parse was accepted as zero"
  exit 0
fi

# --- a control on this file's own inputs -------------------------------------
#
# Every figure below is produced by running a tool with flags. That makes the flag
# a load-bearing part of the measurement, and two of the eleven tools here that
# parse arguments had **no `else` on the parse loop**: `gas_probe` and
# `water_probe` dropped an unrecognised option silently and ran their defaults.
#
# The quiet outcome is the dangerous one. A typo in `water_probe --wave=2.0` does
# not error -- it measures still water, hands the number back under the name of the
# wave study, and this file compares it against a published figure. Sometimes that
# goes red and sends the reader hunting for a defect in the model. Sometimes the
# flag barely moved the answer and it goes green, and the gate has certified an
# experiment nobody ran.
#
# So the tools are asked, before anything is measured, whether they can tell a flag
# they know from one they do not. This is `verify.sh selftest`'s idea one level
# down: a negative control on the instrument rather than on the reading. All twelve
# exit 2 in about 2 ms, so the whole block costs a fortieth of a second.
refuses() {
  local tool="$1"
  [ -x "$tool" ] || return 0
  checks=$((checks + 1))
  if "$tool" --not-a-real-flag >/dev/null 2>&1; then
    printf '  %s✗%s %s accepts an unknown option and runs its defaults instead\n' \
           "$red" "$off" "$tool"
    fails=$((fails + 1))
    return 1
  fi
  return 0
}
for badflag in ./build/gas_probe ./build/water_probe ./build/section_probe \
               ./build/zone_probe ./build/zone_gpu_probe ./build/bulkhead_probe \
               ./build/flip_probe ./build/ram_view ./build/seaway_view \
               ./build/smoke_view ./build/ferry_view ./build/shipsim; do
  refuses "$badflag"
done

# --- the GPU zone solver's published torn counts ---------------------------------
#
# `verify.sh` already runs `zone_gpu_probe`, and its comment there says why: every
# figure in `07-fem-spike-findings.md` §8 comes out of this tool and nothing ran it.
# But that gate only asserts the tool still prints `ok`. **It proves the tool runs;
# it does not prove §8's numbers are still the tool's numbers**, and that is the half
# that kept failing -- three times now in one section:
#
#   - "float tears 60 elements where double tears none", carried by three documents,
#     each quoting the previous rather than the tool. Withdrawn; it never reproduced.
#   - five cells marked as not reproducing, which did reproduce. The re-run that
#     marked them had let the probe *derive* its step count from the punch depth
#     while the table was taken at a fixed 5 505 -- 5 513 steps at 768 and 5 545 at
#     3 072, because the critical timestep falls with element size. Eight extra
#     steps move the torn count by one and the dissipation by 0.07%.
#   - and the second of those propagated into `02-simulation.md` before anyone
#     re-ran it.
#
# So `--steps` is passed explicitly here and not derived. That is the whole defect
# in the second case: a parameter of the experiment that was not written next to the
# number. These are deterministic -- identical across six repeats, and at 1, 4, 12
# and 23 workers -- so the tolerance is zero and any movement at all wants a human.
#
# ~20 s of the `full` gate. The 3 072 invocation cell is the expensive one at ~15 s
# and it is the one that drifted 13%, which is exactly why it is here.
GPU=${GPU:-./build/zone_gpu_probe}
FEM_DOC=docs/07-fem-spike-findings.md
if [ -x "$GPU" ]; then
  gpu768=$("$GPU" --radius=2.5 --sub=8 --steps=5505 --mapping=workgroup --jitter=2e-7 2>&1)
  if printf '%s\n' "$gpu768" | grep -q '^skipped: '; then
    echo "  - no Vulkan device, skipping the §8 torn counts"
    skipped="$skipped zone_gpu_probe"
  else
    for h in "GPU float, torn (workgroup) | **40**" "CPU dissipation | **1.5194**" \
             "GPU float, torn (invocation) | **44**"; do
      hint "$h" "$FEM_DOC"
    done
    # Anchored to its own line for the reason the damage figures above are: the
    # word "torn" appears in the control block as well, and a checker that
    # misparses is worse than no checker.
    # `torn elements` prints two columns and nothing else; `plastic dissipation
    # (MJ)` prints three, and its CPU value is the *fourth* field rather than
    # `NF-1`. Counting back from the end picked the GPU column here once.
    torn_cpu() { printf '%s\n' "$1" | grep '^torn elements' | awk '{ print $3 }'; }
    torn_gpu() { printf '%s\n' "$1" | grep '^torn elements' | awk '{ print $4 }'; }
    diss_cpu() { printf '%s\n' "$1" | grep '^plastic dissipation' | awk '{ print $4 }'; }
    ctl() { printf '%s\n' "$1" | sed -n '/^control:/,$p' | grep "$2" | awk '{ print $NF }'; }
    check "§8 768 el, CPU double torn"  32 0 \
          "$(torn_cpu "$gpu768")" "CPU double, torn | **32**" "$FEM_DOC"
    check "§8 768 el, GPU float torn"   40 0 \
          "$(torn_gpu "$gpu768")" "GPU float, torn (workgroup) | **40**" "$FEM_DOC"
    check "§8 768 el, CPU dissipation (MJ)" 1.5194 0.0001 \
          "$(diss_cpu "$gpu768")" "CPU dissipation | **1.5194**" "$FEM_DOC"
    # **The negative control, and it is the load-bearing cell of the argument.** A
    # geometric jitter the size of float's own representation error, applied to the
    # double solver, must not move the torn count at all -- otherwise "float tears
    # 25% too many" says nothing about float and everything about the problem.
    check "§8 768 el, jittered control torn" 32 0 \
          "$(ctl "$gpu768" 'torn elements')" \
          "control: double, mesh jittered 2 × 10⁻⁷ m, torn" "$FEM_DOC"

    gpu3072=$("$GPU" --radius=2.5 --sub=16 --steps=5505 --mapping=invocation 2>&1)
    check "§8 3 072 el, GPU float torn, invocation mapping" 213 0 \
          "$(torn_gpu "$gpu3072")" "GPU float, torn (invocation) | **44**" \
          "$FEM_DOC"

    # **The 3 072-element workgroup pair, which is the other half of §8's headline
    # and was ungated.** The block above checks 768 and the *invocation* mapping at
    # 3 072; the "52% over" that the section's conclusion rests on is the workgroup
    # kernel's 247 against the double reference's 162, and its jittered control.
    # The control is the load-bearing cell: without it, 247-against-162 is a
    # statement about the problem rather than about float.
    wg3072=$("$GPU" --radius=2.5 --sub=16 --steps=5505 --mapping=workgroup --jitter=2e-7 2>&1)
    check "§8 3 072 el, CPU double torn" 162 0 \
          "$(torn_cpu "$wg3072")" "| CPU double, torn | **32** | **162** |" "$FEM_DOC"
    check "§8 3 072 el, GPU float torn (workgroup)" 247 0 \
          "$(torn_gpu "$wg3072")" "| GPU float, torn (workgroup) | **40** | **247** |" "$FEM_DOC"
    check "§8 3 072 el, jittered control torn" 162 0 \
          "$(ctl "$wg3072" 'torn elements')" \
          "control: double, mesh jittered 2 × 10⁻⁷ m, torn" "$FEM_DOC"

    # **The sharpest assertion available in this document, and it costs half a
    # second.** Below 50 steps the float kernel's enhanced modes never switch on at
    # all: alpha comes back *bit-zero* on every element where the CPU has 7.671e-9.
    # A bit-zero is not a tolerance, it is an identity, and it is the mechanism
    # behind §8's closing negative -- the divergence cannot be the enhanced block,
    # because at this step count there is no enhanced block. `tight` is the control
    # that says so: the same arithmetic with the CPU's stopping rule does *not*
    # stall, and returns 2.010e-8.
    alpha() { printf '%s\n' "$1" | grep '^peak |alpha|' | awk "{ print \$$2 }"; }
    a50f=$("$GPU" --radius=2.5 --sub=4 --steps=50 --mapping=workgroup --eas=float 2>&1)
    a50t=$("$GPU" --radius=2.5 --sub=4 --steps=50 --mapping=workgroup --eas=tight 2>&1)
    check "§8 peak |alpha| at 50 steps, CPU" 7.671e-9 1e-12 \
          "$(alpha "$a50f" 3)" "| 7.671e-9 |" "$FEM_DOC"
    check "§8 peak |alpha| at 50 steps, GPU float (bit-zero)" 0 0 \
          "$(alpha "$a50f" 4)" "| 7.671e-9 |" "$FEM_DOC"
    check "§8 peak |alpha| at 50 steps, GPU tight" 2.010e-8 1e-11 \
          "$(alpha "$a50t" 4)" "2.010e-8" "$FEM_DOC"

    # --- the compiler's own account of each kernel -------------------------------
    #
    # `--stats` was the last output in this tree that nothing ran at all. It is the
    # sole caller of `gpu::describeElementPipelines`, it has no test, and no script
    # passed the flag -- so the two tables it produces, §8's register/spill table
    # and §12's fp64 cost table, were carried on the authority of the single run
    # that first took them. That is the shape this file exists for, and this is a
    # cheap instance of it: one shader compile, no dispatch, nothing timed.
    #
    # The document says of these rows that they "come from the driver's compiler and
    # not from a clock, so they are the one part of this section that a busy box
    # cannot move". That claim is what earns the zero tolerance, and it held: every
    # one of the fourteen figures came back on the digit.
    #
    # `tet_forces.comp` is deliberately absent. §8's footnote already records that
    # `--stats` does not query it and that its cell therefore has a different
    # provenance from its neighbours; gating the rest does not change that, and
    # inventing a check for it here would hide it again.
    stats=$("$GPU" --stats 2>&1)
    # Each field label appears in all six blocks, so every reader below is bounded
    # to one kernel's block and exits at the next `.spv)` header. An unbounded range
    # would silently take whichever block came last -- which is how a check in this
    # file once passed by luck on a multi-line value.
    field() {  # $1 = .spv basename, $2 = first word of the row, $3 = which column
      printf '%s\n' "$stats" | awk -v s="$1" -v w="$2" -v c="$3" '
        !k && index($0, "(" s ")") { k = 1; next }
        k && index($0, ".spv)")    { exit }
        k && $1 == w               { print $c; exit }'
    }
    reg()   { field "$1" Register  3; }
    shmem() { field "$1" Shared    4; }
    spill() { field "$1" spill     4; }

    check "§8 calibration registers/thread" 32 0 \
          "$(reg node_integrate.comp.spv)" \
          '| `node_integrate.comp` (calibration) | 32 | 0 B |' "$FEM_DOC"
    check "§8 one invocation per element, registers" 128 0 \
          "$(reg solidshell_forces.comp.spv)" \
          '| `solidshell_forces.comp` — one invocation | 128 |' "$FEM_DOC"
    check "§8 one invocation per element, spill (B)" 1936 0 \
          "$(spill solidshell_forces.comp.spv)" \
          '**1 936 B = 484 floats**' "$FEM_DOC"
    check "§8 one workgroup per element, registers" 64 0 \
          "$(reg solidshell_forces_wg.comp.spv)" \
          '| `solidshell_forces_wg.comp` — one workgroup | 64 |' "$FEM_DOC"
    check "§8 one workgroup per element, spill (B)" 96 0 \
          "$(spill solidshell_forces_wg.comp.spv)" \
          '**96 B = 24 floats**' "$FEM_DOC"

    # §12's fp64 ladder: three variants against the float baseline, all three
    # columns each. Shared memory is in the table and was ungated with the rest.
    check "§12 fp64 baseline, shared (B)" 2832 0 \
          "$(shmem solidshell_forces_wg.comp.spv)" \
          '| `solidshell_forces_wg.comp` | 64 | 96 B | 2 832 B |' "$FEM_DOC"
    check "§12 fp64 7x7 solve, registers" 72 0 \
          "$(reg solidshell_forces_wg_f64solve.comp.spv)" \
          '| + fp64 7×7 solve | 72 | 128 B | 2 928 B |' "$FEM_DOC"
    check "§12 fp64 7x7 solve, spill (B)" 128 0 \
          "$(spill solidshell_forces_wg_f64solve.comp.spv)" \
          '| + fp64 7×7 solve | 72 | 128 B | 2 928 B |' "$FEM_DOC"
    check "§12 fp64 7x7 solve, shared (B)" 2928 0 \
          "$(shmem solidshell_forces_wg_f64solve.comp.spv)" \
          '| + fp64 7×7 solve | 72 | 128 B | 2 928 B |' "$FEM_DOC"
    check "§12 fp64 condensation, registers" 80 0 \
          "$(reg solidshell_forces_wg_f64condense.comp.spv)" \
          '| + fp64 condensation | 80 | 128 B | 2 928 B |' "$FEM_DOC"
    check "§12 fp64 condensation, spill (B)" 128 0 \
          "$(spill solidshell_forces_wg_f64condense.comp.spv)" \
          '| + fp64 condensation | 80 | 128 B | 2 928 B |' "$FEM_DOC"
    check "§12 fp64 condensation, shared (B)" 2928 0 \
          "$(shmem solidshell_forces_wg_f64condense.comp.spv)" \
          '| + fp64 condensation | 80 | 128 B | 2 928 B |' "$FEM_DOC"
    check "§12 fp64 alpha and Newton, registers" 80 0 \
          "$(reg solidshell_forces_wg_f64newton.comp.spv)" \
          '| + fp64 alpha and Newton | 80 | 128 B | 2 960 B |' "$FEM_DOC"
    check "§12 fp64 alpha and Newton, spill (B)" 128 0 \
          "$(spill solidshell_forces_wg_f64newton.comp.spv)" \
          '| + fp64 alpha and Newton | 80 | 128 B | 2 960 B |' "$FEM_DOC"
    check "§12 fp64 alpha and Newton, shared (B)" 2960 0 \
          "$(shmem solidshell_forces_wg_f64newton.comp.spv)" \
          '| + fp64 alpha and Newton | 80 | 128 B | 2 960 B |' "$FEM_DOC"

    # **And the three claims the section actually argues from, which are derived
    # rather than printed.** Gating the register counts alone would leave the
    # sentences that use them free to drift, which is the failure this whole file
    # was written against -- so the arithmetic is done here, off the measured
    # numbers, rather than trusted. Pascal's register file is 65 536 per SM and a
    # warp is 32 threads, so warps per SM is 65536/(32*registers), floored.
    occupancy() { awk -v r="$1" 'BEGIN { printf "%d\n", int(65536 / (32 * r)) }'; }
    check "§8 warp occupancy, one invocation per element" 16 0 \
          "$(occupancy "$(reg solidshell_forces.comp.spv)")" \
          'occupancy goes from 16 of 64 warps per SM to 32' "$FEM_DOC"
    check "§8 warp occupancy, one workgroup per element" 32 0 \
          "$(occupancy "$(reg solidshell_forces_wg.comp.spv)")" \
          'occupancy goes from 16 of 64 warps per SM to 32' "$FEM_DOC"
    check "§12 warp occupancy under fp64" 25 0 \
          "$(occupancy "$(reg solidshell_forces_wg_f64newton.comp.spv)")" \
          'from 32 warps per SM to 25' "$FEM_DOC"
    # "The remap takes the spill down 20x and halves the register count" -- both
    # halves, off the same two blocks the rows above read.
    ratio2() { awk -v a="$1" -v b="$2" 'BEGIN { printf "%.4f\n", a / b }'; }
    check "§8 the remap's spill reduction" 20.1667 0.001 \
          "$(ratio2 "$(spill solidshell_forces.comp.spv)" \
                    "$(spill solidshell_forces_wg.comp.spv)")" \
          'takes the spill down **20×** and halves the register' "$FEM_DOC"
    check "§8 the remap's register reduction" 2.0 0.001 \
          "$(ratio2 "$(reg solidshell_forces.comp.spv)" \
                    "$(reg solidshell_forces_wg.comp.spv)")" \
          'takes the spill down **20×** and halves the register' "$FEM_DOC"
  fi
else
  echo "  - zone_gpu_probe not built, skipping the §8 torn counts"
  skipped="$skipped zone_gpu_probe"
fi

# --- §1: the number the whole element choice rests on ----------------------------
#
# **`fem_spike` was not run by this gate at all**, which left §1 ungated -- and §1 is
# where "linear tetrahedra are the wrong element for ship plating" comes from. That
# conclusion is one number: a two-elements-through-thickness tet mesh gets **63.6%**
# of the way wrong against a closed form, and the whole solid-shell programme in §4
# and §6 exists because of it. It was the single most load-bearing figure in the
# document and nothing re-derived it.
#
# ~90 s, nearly all of it GPU wait (6% CPU), which is why it sits behind the same
# device check as everything else here -- and **the deflections need the device
# too**, which this said the opposite of. `fem_spike/main.cpp:104-108` builds a
# `FemGpuSolver` inside the refinement loop and returns on failure, so the whole §1
# table is GPU-solved; only the closed-form `beam theory:` line at `:79` is printed
# before any device call. The tool has no way to ask for that one line, so the
# device check stands, but the reason it stands is not the one written here.
SPIKE=${SPIKE:-./build/fem_spike}
if [ -x "$SPIKE" ]; then
  spike=$("$SPIKE" 2>&1)
  if printf '%s\n' "$spike" | grep -qiE 'no usable GPU|no Vulkan'; then
    echo "  - no Vulkan device, skipping §1's deflection table"
    skipped="$skipped fem_spike"
  else
    # Anchored to the row's own tet count, because "0.0801" and "63.6" both appear
    # in prose elsewhere in the document and a checker that misparses is worse than
    # no checker.
    row() { printf '%s\n' "$1" | awk -v t="$2" '$2 == t { print $4 }'; }
    err() { printf '%s\n' "$1" | awk -v t="$2" '$2 == t { print $5 }' | tr -d '%'; }
    check "§1 beam theory tip deflection (mm)" 0.2199 0 \
          "$(printf '%s\n' "$spike" | sed -n 's/.*beam theory: \([0-9.]*\) mm.*/\1/p')" \
          "Theory: 0.2199 mm." "$FEM_DOC"
    check "§1 2 elements through thickness, deflection (mm)" 0.0801 0 \
          "$(row "$spike" 480)"  "| 2 | 480 | 0.0801 mm | 63.6% |" "$FEM_DOC"
    check "§1 2 through, error vs theory (%)"  63.6 0 \
          "$(err "$spike" 480)"  "| 2 | 480 | 0.0801 mm | 63.6% |" "$FEM_DOC"
    check "§1 4 through, error vs theory (%)"  32.2 0 \
          "$(err "$spike" 3840)" "| 4 | 3 840 | 0.1491 mm | 32.2% |" "$FEM_DOC"
    check "§1 8 through, error vs theory (%)"  10.8 0 \
          "$(err "$spike" 30720)" "| 8 | 30 720 | 0.1962 mm | 10.8% |" "$FEM_DOC"

    # --- §2 and §3, out of the same run ------------------------------------
    #
    # This block parsed §1 and threw the rest of `$spike` away, while §2 and §3
    # sat in it costing nothing. The document says at :72-74 that every other cell
    # of both tables "comes back on `./build/fem_spike` exactly as printed, to
    # every digit", which is the strongest possible claim to leave unchecked.
    #
    # **Range-restricted, and the `exit` matters.** `awk '$1 == 1'` also matches
    # the section *title* line, because `"1."` is a numeric string to awk and
    # compares equal to 1 -- and §4's line carries `$1 = 46080`, which collides
    # with §3's third row. Anchoring on the section heading and stopping at the
    # first hit is what makes each of these one number instead of several.
    sec2() {
      printf '%s\n' "$spike" |
        awk -v steps="$1" -v col="$2" \
            '/^2\. GPU kernel/ { f = 1 }
             /^3\. Throughput/ { f = 0 }
             f && NF == 3 && $1 == steps { print $col; exit }'
    }
    sec3() {
      printf '%s\n' "$spike" |
        awk -v tets="$1" -v col="$2" \
            '/^3\. Throughput/ { f = 1 }
             /^4\. Single-threaded/ { f = 0 }
             f && NF == 5 && $1 == tets { print $col; exit }'
    }

    # §2 is a CPU-against-GPU difference and contains no wall clock at all, so
    # every cell of it is gateable. Tolerances are the rounding the document
    # published at -- three significant figures on the velocity column, two on
    # the position column, where the tool prints three.
    check "§2 1 step, RMS dv / RMS v" 1.96e-5 5e-8 "$(sec2 1 2)" \
          "| 1 | 1.96 × 10⁻⁵ | 1.5 × 10⁻⁸ m |" "$FEM_DOC"
    check "§2 1 step, max |dx| (m)" 1.5e-8 5e-10 "$(sec2 1 3)" \
          "| 1 | 1.96 × 10⁻⁵ | 1.5 × 10⁻⁸ m |" "$FEM_DOC"
    check "§2 10 steps, RMS dv / RMS v" 1.45e-4 5e-7 "$(sec2 10 2)" \
          "| 10 | 1.45 × 10⁻⁴ | 3.0 × 10⁻⁸ m |" "$FEM_DOC"
    check "§2 10 steps, max |dx| (m)" 3.0e-8 5e-10 "$(sec2 10 3)" \
          "| 10 | 1.45 × 10⁻⁴ | 3.0 × 10⁻⁸ m |" "$FEM_DOC"
    check "§2 100 steps, RMS dv / RMS v" 1.21e-3 5e-6 "$(sec2 100 2)" \
          "| 100 | 1.21 × 10⁻³ | 2.7 × 10⁻⁷ m |" "$FEM_DOC"
    check "§2 100 steps, max |dx| (m)" 2.7e-7 5e-9 "$(sec2 100 3)" \
          "| 100 | 1.21 × 10⁻³ | 2.7 × 10⁻⁷ m |" "$FEM_DOC"
    check "§2 1000 steps, RMS dv / RMS v" 4.61e-3 5e-6 "$(sec2 1000 2)" \
          "| 1 000 | 4.61 × 10⁻³ | 6.0 × 10⁻⁷ m |" "$FEM_DOC"
    check "§2 1000 steps, max |dx| (m)" 6.0e-7 5e-9 "$(sec2 1000 3)" \
          "| 1 000 | 4.61 × 10⁻³ | 6.0 × 10⁻⁷ m |" "$FEM_DOC"

    # **§3's two count columns and not one of its three clock columns.** The
    # document says so itself at :100-102 -- "the two count columns re-derive; the
    # three clock columns never have" -- and a run here confirms it: `ms/step`
    # came back 0.0079 against a published 0.0083 and `Melem-upd/s` 121.6 against
    # 116, on the same binary and the same mesh. The counts are pure functions of
    # the literal `sizes[]` in `fem_spike/main.cpp`, so they carry tolerance 0.
    for pair in "960 315" "11520 2745" "46080 10285" "230400 45225" "491520 95337"; do
      set -- $pair
      check "§3 $1 tets, node count" "$2" 0 "$(sec3 "$1" 2)" \
            "| $(printf '%s' "$1" | sed -E 's/([0-9])([0-9]{3})$/\1 \2/') |" "$FEM_DOC"
    done
  fi
else
  echo "  - fem_spike not built, skipping §1's deflection table"
  skipped="$skipped fem_spike"
fi

# --- §6's locking table, which the suite already computes -------------------------
#
# Free: `shipsim_tests` prints this on every run and the gate already runs it. Only
# the document-side comparison was missing. It is the justification for the element
# change -- at the slenderness of real plating a plain hex is 1 400x too stiff and a
# linear tet 3 800x -- and the two stiff columns are quoted in the document to one
# significant figure, so the expectation is the *published* value at the rounding it
# was published at, exactly as the front page's Cb is.
#
# `TESTS` is set for the README block further down, which is *after* this one, and
# `set -u` turns reading it early into an abort. Defaulted here rather than moved,
# because the assignment below must keep working when this block is skipped. The
# abort was caught the honest way, by running it: the gate exited 1 and printed no
# success line, which is the behaviour the zero-checks fix above exists to guarantee.
: "${TESTS:=./build/shipsim_tests}"
if [ -x "$TESTS" ]; then
  lock=$("$TESTS" 2>&1 | awk '/locking, demonstrated/,/explicit stability limit/' | awk '$1 == 500')
  check "§6 locking at L/t=500, solid-shell"   0.996  5e-4 "$(echo "$lock" | awk '{print $2}')" \
        "| 500 | 0.996 | 0.813 | **0.0007** | **0.0003** |" "$FEM_DOC"
  check "§6 locking at L/t=500, ANS no EAS"    0.813  5e-4 "$(echo "$lock" | awk '{print $3}')" \
        "| 500 | 0.996 | 0.813 | **0.0007** | **0.0003** |" "$FEM_DOC"
  check "§6 locking at L/t=500, plain hex"     0.0007 5e-5 "$(echo "$lock" | awk '{print $4}')" \
        "a plain hex is 1 400× too stiff" "$FEM_DOC"
  check "§6 locking at L/t=500, linear tet"    0.0003 5e-5 "$(echo "$lock" | awk '{print $5}')" \
        "tet 3 800× too stiff" "$FEM_DOC"
else
  echo "  - shipsim_tests not built, skipping §6's locking table"
  skipped="$skipped shipsim_tests-locking"
fi

# --- the Phase 4 milestone -------------------------------------------------------
#
# `docs/06-roadmap.md`'s Phase 4 milestone publishes a table and a restraint window,
# and both come out of `bulkhead_probe`. **The window is the figure most worth
# gating**: it is the band of a parameter the model cannot derive over which the
# milestone's own sentence is true at all, and a change anywhere in the fire, the
# conduction solve or the strength model moves it without moving anything that fails
# a test. `verify.sh` already runs the tool for its `ok` contract; this is the other
# half, and it is the same failure the smoke and section figures above were added
# for -- a number nobody re-ran.
#
# ~45 s, and it is a second full run rather than a parse of the first because
# `expect_ok` throws its log away. The alternative -- one run feeding both -- would
# make the acceptance test and the figures the same invocation, and the acceptance
# test is the one that must be able to fail on its own.
MILESTONE=${MILESTONE:-./build/bulkhead_probe}
if [ -x "$MILESTONE" ]; then
  bulk=$("$MILESTONE" --quiet 2>&1)
  for h in "member peak | 252.2 °C" "**0.2194 ≤ r < 0.2505**" "**77.4 K hotter**" \
           "**2.086**" "0.747 and has not failed"; do
    hint "$h"
  done
  # Each reader anchored to its own line and the **label stripped before any field
  # is counted**, because two of the three labels are two words and one is one word,
  # so `$2` is "alone" on two rows and "yes" on the third. Counting back from the end
  # would work here and is the thing that picked the wrong column in the §8 block
  # above, so the label goes instead and the fields are counted forward: 1 failed,
  # 2 at, 3 steel C, 4 member C, 5 worst u, 6 members, 7 hole, 8 into ER, 9 wet.
  row() { printf '%s\n' "$bulk" | grep "^$1  *" | sed "s/^$1  *//"; }

  # **The second leak path, which is the milestone's other control and was ungated.**
  # The ferry carries an unsealed 0.04 m2 cable transit through this bulkhead. The
  # milestone seals it; `--cable-transit` reopens it, and then the water alone fells
  # the bulkhead at 3645 s with nothing burning -- the figure `06-roadmap.md`
  # publishes and nothing re-derived.
  #
  # **This run exits non-zero on purpose, and that is asserted rather than ignored.**
  # `bulkhead_probe`'s own `require` says the head alone must not fell the bulkhead,
  # which is true of the milestone configuration and is exactly what reopening the
  # transit inverts. So the check is that it fails *once*, for *that* reason: an
  # ignored exit status would be a hole, and a run that started failing for some
  # other reason would slip through a bare `|| true`.
  #
  # `--duration=3700` rather than the default 3600, because the event is at 3645 and
  # a run that stops before it reports nothing. 99 s.
  transit=$("$MILESTONE" --quiet --cable-transit --duration=3700 2>&1)
  trow() { printf '%s\n' "$transit" | grep "^$1  *" | sed "s/^$1  *//"; }
  check "cable transit open: the water alone fells the bulkhead at (s)" 3645 0 \
        "$(trow 'head alone' | awk '{ print $2 }')" \
        "the **water-only control at 3645 s with nothing burning anywhere**"
  checks=$((checks + 1))
  if printf '%s\n' "$transit" | grep -q '^1 check(s) failed$' &&
     printf '%s\n' "$transit" | grep -q '^  ! the head alone felled the bulkhead$'; then
    printf '  %s✓%s reopening the cable transit inverts exactly one milestone control\n' \
           "$green" "$off"
  else
    printf '  %s✗%s the cable-transit run did not fail in the one way it is supposed to\n' \
           "$red" "$off"
    printf '      %sexpected "1 check(s) failed" and "the head alone felled the bulkhead"%s\n' \
           "$dim" "$off"
    fails=$((fails + 1))
  fi
  window_lo=$(printf '%s\n' "$bulk" | sed -n 's/^restraint window: .* true for \([0-9.]*\) <= r.*/\1/p')
  window_hi=$(printf '%s\n' "$bulk" | sed -n 's/^restraint window: .*<= r < \([0-9.]*\),.*/\1/p')
  check "the restraint window's lower bound" 0.2194 0.0005 "$window_lo" \
        "**0.2194 ≤ r < 0.2505**"
  check "the restraint window's upper bound" 0.2505 0.0005 "$window_hi" \
        "**0.2194 ≤ r < 0.2505**"
  check "the fire-alone control's peak member temperature (C)" 252.2 0.2 \
        "$(row 'fire alone' | awk '{ print $4 }')" "member peak | 252.2 °C"
  check "the fire-alone control's worst utilisation" 0.934 0.002 \
        "$(row 'fire alone' | awk '{ print $5 }')" "| 0.934 |"
  check "the head-alone control's worst utilisation" 0.274 0.002 \
        "$(row 'head alone' | awk '{ print $5 }')" "| 0.274 |"
  check "the coupled case's worst utilisation" 1.001 0.002 \
        "$(row 'both' | awk '{ print $5 }')" "| **1.001** |"
  check "the failure time (s)" 1935 5 "$(row 'both' | awk '{ print $2 }')" "at 1935 s"
  check "water into the machinery spaces (m3)" 307.7 1.0 \
        "$(row 'both' | awk '{ print $8 }')" "307.7 m³"
  check "the hole the failure opened (m2)" 0.490 0.001 \
        "$(row 'both' | awk '{ print $7 }')" "0.490 m²"
  magnifier=$(printf '%s\n' "$bulk" | sed -n 's/.*magnified    [0-9.]*   (x\([0-9.]*\),.*/\1/p')
  additive=$(printf '%s\n' "$bulk" | sed -n 's/^  utilisation [0-9.]*, .* additive check reads \([0-9.]*\) and.*/\1/p')
  check "the beam-column magnification at failure" 2.086 0.005 "$magnifier" "**2.086**"
  check "what a purely additive check would have read" 0.747 0.002 "$additive" \
        "0.747 and has not failed"
  cooler=$(printf '%s\n' "$bulk" | sed -n 's/.*leaves the member \([0-9.]*\) K hotter.*/\1/p')
  check "how much hotter a dry hold leaves the member (K)" 77.4 0.3 "$cooler" \
        "**77.4 K hotter**"
  # **And the finding itself, not only its numbers.** Each figure above could match
  # while the milestone's own sentence had stopped reproducing, so the relation is
  # asserted as a relation.
  checks=$((checks + 1))
  if printf '%s\n' "$bulk" | grep -q '^ok$' &&
     [ "$(row 'fire alone' | awk '{ print $1 }')" = "no" ] &&
     [ "$(row 'head alone' | awk '{ print $1 }')" = "no" ] &&
     [ "$(row 'both' | awk '{ print $1 }')" = "yes" ]; then
    printf '  %s✓%s and the milestone still reproduces its own sentence: neither cause alone,'\
' both together\n' "$green" "$off"
  else
    printf '  %s✗%s the milestone no longer reproduces its own sentence\n' "$red" "$off"
    printf '      %s%s publishes the three-run table%s\n' "$dim" "$DOC" "$off"
    fails=$((fails + 1))
  fi
else
  echo "  - bulkhead_probe not built, skipping the Phase 4 milestone figures"
  skipped="$skipped bulkhead_probe"
fi

# --- the front page ---------------------------------------------------------------
#
# **98 gated figures and not one of them on `README.md`.** This script grew a
# document at a time -- the roadmap, then the FEM findings, then the simulation and
# renderer docs -- and the one page every reader opens first was never in it. It had
# drifted about as far as that suggests:
#
#   - `16 compartments` for 18. The mid wing tanks were added when a ram amidships
#     turned out to tear 26 bays open onto nothing, and the front page still counted
#     the ship without them.
#   - a compartment table from *before* the wing tank authored inside a hold was
#     pulled back out of it: `fwd_hold_s` published at 1074 m3, which is 965 plus the
#     108 m3 tank that was double-counted inside it. The same drift the header of
#     this file was written about, one document further out.
#   - two mutually contradictory counts of the same test suite, `116` and `1309`,
#     against a measured 198 869.
#
# **The durations are the interesting part and they are now written into the README.**
# `full` is still flooding when the run ends, so its figures are a function of how
# long anyone watched: -3.30 m, 7.2 deg and 1442 t at t+900, -3.23 m, 9.3 deg and
# 1556 t at t+1800. `06-roadmap.md` publishes the first and the README the second,
# which reads as a contradiction until the duration is written next to each -- and
# a reader reconciling them by editing one to match the other would be destroying a
# correct figure. Both are re-run here, each at the length its own document claims
# it at.
#
# **And the verdicts, not only the figures.** `full` went SURVIVED -> LOST one commit
# before this block existed, on a change to how GM is sampled that moved no tonnage,
# heel or draft at all. Every figure in that row would have matched while the row's
# own verdict had reversed, so each scenario's verdict is asserted as well.
#
# ~230 s, of which 112 s is the test suite: the front page publishes its check count
# and nothing else here re-derives it. Adding a test is *expected* to fail this, and
# updating the README line is part of adding the test.
SHIPSIM=${SHIPSIM:-./build/shipsim}
TESTS=${TESTS:-./build/shipsim_tests}
FRONT=README.md
SHIPFILE=ships/ferry.ship

if [ -x "$SHIPSIM" ]; then
  # The intact banner, which needs no flooding: `--duration=0` runs no steps.
  intact=$("$SHIPSIM" --scenario=none --duration=0 2>&1)
  disp=$(printf '%s\n' "$intact" | sed -n 's/^  displacement  *\([0-9.]*\) t/\1/p')
  cb=$(printf '%s\n'   "$intact" | sed -n 's/^  block coeff Cb  *\([0-9.]*\)/\1/p')
  # The trailing `.*` is not decoration: `s/.../\1/` replaces only what it matched
  # and leaves the rest of the line, so without it this reads back "2.00 m" -- which
  # `awk` then converts to 2.00 and passes, having parsed a unit as a number.
  gm0=$(printf '%s\n'  "$intact" | sed -n 's/^  GM (transverse)  *\(-\{0,1\}[0-9.]*\).*/\1/p')
  # Cb and GM are quoted to two figures on the front page where the tool prints
  # three, so the expectation is the *published* value and the tolerance is the
  # rounding it published at. Anything looser would accept a different ship.
  check "the ferry's intact displacement (t)" 8984 1 "$disp" \
        "8984 t, Cb 0.66" "$FRONT"
  check "the ferry's block coefficient" 0.66 0.005 "$cb" \
        "Cb 0.66, intact GM 2.00 m" "$FRONT"
  check "the ferry's intact GM (m)" 2.00 0.005 "$gm0" \
        "intact GM 2.00 m, 18 compartments" "$FRONT"
  # Counted off `ships/ferry.ship` rather than the compiled ferry, because no tool
  # prints a compartment count -- and `verify.sh` already requires the two to reach
  # identical outcome strings over 900 s in all three scenarios, so the file's count
  # is the ship's count or that gate is already red.
  check "compartments in the ferry" 18 0 "$(grep -c '^compartment ' "$SHIPFILE")" \
        "18 compartments carved" "$FRONT"

  # Each row of the scenario table, at the duration the table now names. The
  # readers are anchored to their own lines for the reason the damage figures are:
  # `heel` opens the outcome summary *and* the two righting-arm curves, and `GM`
  # appears in the banner as well as the verdict.
  fin_heel()  { printf '%s\n' "$1" | sed -n 's/^  heel \(-\{0,1\}[0-9.]*\) deg.*/\1/p'; }
  # **`\([0-9]*\) t of floodwater` on its own reads back empty**, and silently: the
  # greedy `.*` in front of it eats the digits, the group matches zero characters and
  # is happy, and `check` then reports "could not parse" for a figure that is right
  # there. It needs a delimiter the `.*` cannot cross and at least one digit.
  fin_water() { printf '%s\n' "$1" | sed -n 's/.*, \([0-9][0-9]*\) t of floodwater.*/\1/p'; }
  fin_gm()    { printf '%s\n' "$1" | sed -n 's/^  effective GM \(-\{0,1\}[0-9.]*\) m.*/\1/p'; }
  fin_t()     { printf '%s\n' "$1" | sed -n 's/^=== Outcome at t+\([0-9]*\)s.*/\1/p'; }
  verdict()   { printf '%s\n' "$1" | sed -n 's/^=== Outcome at t+[0-9]*s: \(.*\) ===$/\1/p'; }
  # A compartment row is `name gross fill% kPa water`, and every name here is one
  # token, so the fields count forward from it.
  space()     { printf '%s\n' "$1" | grep "^  $2  *" | awk -v c="$3" '{ print $c }'; }

  none=$("$SHIPSIM" --scenario=none --duration=1800 2>&1)
  # The first report at which GM has gone negative. The run table is the only place
  # this is visible and it is printed every 15 s, which is the tolerance. `NF == 8`
  # alone would also match a five-word event line, so the first field has to be a
  # time as well.
  gmneg=$(printf '%s\n' "$none" |
          awk 'NF == 8 && $1 ~ /^[0-9]+$/ && $5 + 0 < 0 { print $1; exit }')
  check "'none': when GM goes negative (s)" 690 15 "$gmneg" \
        "GM negative at t+690 s" "$FRONT"
  check "'none': the angle she lolls to (deg)" 53 0.5 "$(fin_heel "$none")" \
        "lolls to **53° by t+1800 s**" "$FRONT"
  check "'none': floodwater at t+1800 (t)" 6311 1 "$(fin_water "$none")" \
        "6311 t aboard at the end" "$FRONT"

  doors=$("$SHIPSIM" --scenario=doors --duration=1800 2>&1)
  check "'doors': when she capsizes (s)" 930 15 "$(fin_t "$doors")" \
        "**capsizes at t+930 s**" "$FRONT"
  # The pointer carries the deck's dimensions as well as the tonnage, because the
  # 19 m is a figure too and no tool prints it: the vehicle deck's plan area comes
  # out at 1868 m² over the 100 m box, 18.7 m of mean breadth.
  check "'doors': water on the vehicle deck at capsize (t)" 3950 1 \
        "$(space "$doors" vehicle_deck 5)" "100 × 19 m undivided deck — 3950 t of it by the" \
        "$FRONT"
  # The trapped-air table, every cell of it. This is the block that had gone stale,
  # and it went stale in the *gross volume* column -- a number that cannot move
  # unless the ship itself changes -- so all four columns are checked and not only
  # the ones the flooding drives.
  #
  # One pointer per row and it is the row's whole numeric tail, verbatim as the
  # README prints it. That is what makes an edit *to the README* fail here: `check`
  # compares the tool against a constant in this file, so a doc that drifts on its
  # own is caught by the pointer or it is not caught at all.
  air_p="965     20.6    127.7       194"
  air_w="108     32.5    150.1        35"
  air_a="1232     25.9    123.2       311"
  check "'doors': fwd_hold_s gross volume (m3)" 965   1   "$(space "$doors" fwd_hold_s 2)" "$air_p" "$FRONT"
  check "'doors': fwd_hold_s fill (%)"          20.6  0.1 "$(space "$doors" fwd_hold_s 3)" "$air_p" "$FRONT"
  check "'doors': fwd_hold_s air pressure (kPa)" 127.7 0.2 "$(space "$doors" fwd_hold_s 4)" "$air_p" "$FRONT"
  check "'doors': fwd_hold_s water (t)"         194   1   "$(space "$doors" fwd_hold_s 5)" "$air_p" "$FRONT"
  check "'doors': wing_tank_fwd_s gross volume (m3)" 108 1 "$(space "$doors" wing_tank_fwd_s 2)" "$air_w" "$FRONT"
  check "'doors': wing_tank_fwd_s fill (%)"     32.5  0.1 "$(space "$doors" wing_tank_fwd_s 3)" "$air_w" "$FRONT"
  check "'doors': wing_tank_fwd_s air pressure (kPa)" 150.1 0.2 "$(space "$doors" wing_tank_fwd_s 4)" "$air_w" "$FRONT"
  check "'doors': wing_tank_fwd_s water (t)"    35    1   "$(space "$doors" wing_tank_fwd_s 5)" "$air_w" "$FRONT"
  check "'doors': aft_hold_s gross volume (m3)" 1232  1   "$(space "$doors" aft_hold_s 2)" "$air_a" "$FRONT"
  check "'doors': aft_hold_s fill (%)"          25.9  0.1 "$(space "$doors" aft_hold_s 3)" "$air_a" "$FRONT"
  check "'doors': aft_hold_s air pressure (kPa)" 123.2 0.2 "$(space "$doors" aft_hold_s 4)" "$air_a" "$FRONT"
  check "'doors': aft_hold_s water (t)"         311   1   "$(space "$doors" aft_hold_s 5)" "$air_a" "$FRONT"
  # No separate check that these stand at "1.2 to 1.5 atmospheres": three pressures
  # gated to 0.2 kPa already say it, and a relation that cannot fail while the
  # figures above hold is coverage that is not there.

  full=$("$SHIPSIM" --scenario=full --duration=1800 2>&1)
  check "'full': GM at t+1800 (m)" -3.23 0.005 "$(fin_gm "$full")" \
        "GM −3.23 m under 6 mm of water" "$FRONT"
  check "'full': list at t+1800 (deg)" 9.3 0.05 "$(fin_heel "$full")" \
        "9.3° list, 1556 t" "$FRONT"
  check "'full': floodwater at t+1800 (t)" 1556 1 "$(fin_water "$full")" \
        "9.3° list, 1556 t" "$FRONT"
  check "'full': water on the vehicle deck (t)" 11 1 "$(space "$full" vehicle_deck 5)" \
        "under 6 mm of water on the vehicle deck" "$FRONT"
  # **And the same run at t+900, which is the figure the README and the roadmap read
  # as disagreeing about.** Gating both is the point: with only one of them checked,
  # the obvious way to reconcile two documents is to overwrite whichever one is not
  # gated, and it is a correct measurement of a different instant.
  #
  # Free rather than a second run -- the row the 1800 s run prints at t = 900 is the
  # state a 900 s run ends on, the trajectory to there being the same one. Note the
  # table prints heel to 2 dp (7.19) where the outcome line rounds to 1 (7.2).
  row900=$(printf '%s\n' "$full" | awk 'NF == 8 && $1 == 900')
  check "'full' at t+900: GM (m)" -3.30 0.005 \
        "$(printf '%s\n' "$row900" | awk '{ print $5 }')" \
        "at −3.30 m, 7.2° and 1442 t" "$FRONT"
  check "'full' at t+900: heel (deg)" 7.2 0.05 \
        "$(printf '%s\n' "$row900" | awk '{ print $3 }')" \
        "at −3.30 m, 7.2° and 1442 t" "$FRONT"
  check "'full' at t+900: floodwater (t)" 1442 1 \
        "$(printf '%s\n' "$row900" | awk '{ print $7 }')" \
        "at −3.30 m, 7.2° and 1442 t" "$FRONT"

  # **And the three verdicts.** The figures above are what the ship is doing; these
  # are what the front page says has happened to her, and the two moved apart once
  # already -- `full` went SURVIVED -> LOST with every tonnage, heel and draft in
  # that row unchanged to the digit.
  checks=$((checks + 1))
  if [ "$(verdict "$none")"  = "LOST - lolled over with negative GM, flooding continuing" ] &&
     [ "$(verdict "$doors")" = "CAPSIZED" ] &&
     [ "$(verdict "$full")"  = "LOST - negative GM, loll imminent" ]; then
    printf '  %s✓%s and the three verdicts still read lolled over / capsized / lost\n' \
           "$green" "$off"
  else
    printf '  %s✗%s a scenario verdict has changed: none=%s doors=%s full=%s\n' \
           "$red" "$off" "$(verdict "$none")" "$(verdict "$doors")" "$(verdict "$full")"
    printf '      %s%s publishes these as nothing / capsizes / lost%s\n' \
           "$dim" "$FRONT" "$off"
    fails=$((fails + 1))
  fi
else
  echo "  - shipsim not built, skipping the README's scenario figures"
  skipped="$skipped shipsim-scenarios"
fi

# The front page's own count of the validation suite, which is the one figure here
# that a *test* changes rather than the physics. 112 s, and it is the only thing in
# the repository that re-derives it: `verify.sh` prints the count and asserts only
# that the failures are zero.
#
# **The pointer is derived from the figure rather than written beside it**, and it
# is the only call here that needs to be. Every other figure is a physical result
# that moves when the physics moves; this one moves whenever *anyone adds a test*,
# which in practice means it goes stale on the same commit that changes it. Held
# as two independent copies of one number it then failed twice -- once as the
# figure and once as the pointer string that no longer occurred in the README --
# and that happened three times in this session alone, each time to whoever had
# just added a test. One `expected` variable cannot disagree with itself.
if [ -x "$TESTS" ]; then
  # **One run, read many times.** This piped the suite straight into `sed` and
  # threw the rest away. The suite is 112 s and it already prints every figure
  # `02-simulation.md` publishes about Ikeda's roll damping, so capturing stdout
  # instead of streaming it costs nothing and closes the largest coverage hole
  # this script had: **not one figure in it was computed with Ikeda attached.**
  # `attachRollDamping` is reachable only through `shipsim --bilge-keels=`, a flag
  # this script never passes, so no run it performed ever built a
  # `RollDampingHull` at all -- a whole validated subsystem, fifteen published
  # numbers, and the only thing re-deriving them was whoever last read the test's
  # own stdout. The tell was not a wrong figure; it was asking which *flags* the
  # gate never passes.
  suiteout=$("$TESTS" 2>&1)
  suite=$(printf '%s\n' "$suiteout" | sed -n 's/^\([0-9]*\) checks, [0-9]* failures$/\1/p' | tail -1)
  expected=201029
  check "closed-form validation checks in the suite" "$expected" 0 "$suite" \
        "$expected validation checks" "$FRONT"
  # **The roadmap publishes the same count in a different format**, and it was
  # 198 869 while the README said 200 336 -- stale by two commits and invisible to
  # every grep aimed at the README's spelling, because this one is written with a
  # thin space. A figure gated in one document and loose in another is the drift
  # the whole script exists to stop; the count is derived from the same
  # `$expected` so the two cannot part company.
  hint "$(printf '%s %s closed-form validation checks' \
          "${expected%???}" "${expected#${expected%???}}")" "$DOC"

  # --- what a ram takes out of the hull girder, out of the same run ----------------
  #
  # `02-simulation.md`'s stiffener-loss table and `06-roadmap.md`'s restatement of
  # it. These were *printed* all along -- `test_promotion.cpp` emits every one on
  # every run -- and four of the twelve had still drifted, which is the sharper
  # version of this file's usual finding: printed is not the same as read.
  #
  # What hid it was the assertion around them. The test guards these with one-sided
  # bounds -- `gap.both_ > gap.plate`, `> 1.15x`, `> 1.35x` -- against actuals of
  # 1.23x, 1.25x and 1.44x. All still pass. So when `4b35bc7` (the plate's `b`),
  # `60c96ad` (the column-buckling mode) and `3519e9d` (an explicit yield strength)
  # each correctly moved `collapseCurve`, the hogging and sagging rows moved with
  # them and nothing anywhere went red. The section-area and second-moment rows did
  # *not* move, correctly: they come from `hullGirderSection` and never touch
  # buckling, which is the cross-check that says the code is right and the table was
  # three fixes stale.
  #
  # The bounds are left alone deliberately. They assert a claim about the physics --
  # that losing the longitudinals costs strictly more, by a factor and not a
  # rounding -- and that claim should not be re-tuned every time the collapse model
  # improves. Pinning the published digits is this file's job, not theirs.
  #
  # Free: the suite is already captured above.
  gap() {  # $1 = row label, $2 = 1 plating alone, 2 with members, 3 the ratio
    printf '%s\n' "$suiteout" |
      sed -n "s/^ *$1 *plating alone \([0-9.]*\)%, with the longitudinals \([0-9.]*\)%  (\([0-9.]*\)x)$/\\$2/p" |
      head -1
  }
  check "ram damage: section area lost, plating alone" 6.861 0 "$(gap 'section area' 1)" \
        '| section area lost | 6.861% | **8.455%** | 1.23× |' "$SIM_DOC"
  check "ram damage: section area lost, with the longitudinals" 8.455 0 \
        "$(gap 'section area' 2)" '| section area lost | 6.861% | **8.455%** | 1.23× |' "$SIM_DOC"
  check "ram damage: section area ratio" 1.23 0 "$(gap 'section area' 3)" \
        '| section area lost | 6.861% | **8.455%** | 1.23× |' "$SIM_DOC"
  check "ram damage: second moment lost, plating alone" 5.429 0 "$(gap 'second moment' 1)" \
        '| second moment lost | 5.429% | **6.590%** | 1.21× |' "$SIM_DOC"
  check "ram damage: second moment lost, with the longitudinals" 6.590 0 \
        "$(gap 'second moment' 2)" '| second moment lost | 5.429% | **6.590%** | 1.21× |' "$SIM_DOC"
  check "ram damage: second moment ratio" 1.21 0 "$(gap 'second moment' 3)" \
        '| second moment lost | 5.429% | **6.590%** | 1.21× |' "$SIM_DOC"
  check "ram damage: hogging ultimate moment lost, plating alone" 5.388 0 \
        "$(gap 'hogging ultimate moment' 1)" \
        '| hogging ultimate moment lost | 5.388% | **6.747%** | 1.25× |' "$SIM_DOC"
  check "ram damage: hogging ultimate moment lost, with the longitudinals" 6.747 0 \
        "$(gap 'hogging ultimate moment' 2)" \
        '| hogging ultimate moment lost | 5.388% | **6.747%** | 1.25× |' "$SIM_DOC"
  check "ram damage: hogging ratio" 1.25 0 "$(gap 'hogging ultimate moment' 3)" \
        '| hogging ultimate moment lost | 5.388% | **6.747%** | 1.25× |' "$SIM_DOC"
  check "ram damage: sagging ultimate moment lost, plating alone" 12.369 0 \
        "$(gap 'sagging ultimate moment' 1)" \
        '| sagging ultimate moment lost | 12.369% | **17.854%** | **1.44×** |' "$SIM_DOC"
  check "ram damage: sagging ultimate moment lost, with the longitudinals" 17.854 0 \
        "$(gap 'sagging ultimate moment' 2)" \
        '| sagging ultimate moment lost | 12.369% | **17.854%** | **1.44×** |' "$SIM_DOC"
  check "ram damage: sagging ratio" 1.44 0 "$(gap 'sagging ultimate moment' 3)" \
        '| sagging ultimate moment lost | 12.369% | **17.854%** | **1.44×** |' "$SIM_DOC"

  # The extent of the damage the table is *about*. Both documents quote it, and the
  # test bounds it only through a derived spacing window wide enough that the area
  # and the length could move a long way together and still pass.
  ram() {  # $1 = 1 area m2, 2 panels, 3 length m, 4 members
    printf '%s\n' "$suiteout" |
      sed -n "s/^ *the ram opens \([0-9.]*\) m2 of shell over \([0-9]*\) panels and \([0-9.]*\) m of longitudinal over \([0-9]*\) members.*/\\$1/p" |
      head -1
  }
  check "ram damage: shell opened (m2)" 125.6 0.05 "$(ram 1)" \
        'a ram opening 125.6 m² of side shell' "$SIM_DOC"
  check "ram damage: panels opened" 72 0 "$(ram 2)" \
        'over 72 panels and the 163 m of longitudinal running through them' "$SIM_DOC"
  check "ram damage: longitudinal opened (m)" 163.2 0.05 "$(ram 3)" \
        'over 72 panels and the 163 m of longitudinal running through them' "$SIM_DOC"

  # And the zone control: the shipped-before model with `fiberFailure` off. Its
  # four figures are the argument for the criterion existing at all.
  ctlfibre() {  # $1 = 1 eps_p, 2 the multiple, 3 the failure strain
    printf '%s\n' "$suiteout" |
      sed -n "s/^ *the control leaves the longitudinal at eps_p = \([0-9.]*\), \([0-9.]*\)x its own regularised failure strain of \([0-9.]*\).*/\\$1/p" |
      head -1
  }
  unconservative() {  # $1 = 1 force, 2 stored energy, 3 punch work %
    printf '%s\n' "$suiteout" |
      sed -n "s/^ *so the un-conservative model claims \([0-9.]*\)x the fibre force, \([0-9.]*\)x the stored energy, and \([0-9.]*\)% more punch work.*/\\$1/p" |
      head -1
  }
  check "fiberFailure off: plastic strain left on the longitudinal" 0.5128 0.0001 \
        "$(ctlfibre 1)" '`ε_p = 0.513`, **2.83× its own failure strain**' "$SIM_DOC"
  check "fiberFailure off: multiple of its own failure strain" 2.83 0.005 "$(ctlfibre 2)" \
        '`ε_p = 0.513`, **2.83× its own failure strain**' "$SIM_DOC"
  check "fiberFailure off: multiple of the nodal force" 2.23 0.005 "$(unconservative 1)" \
        'carrying 2.23× the nodal force and 5.74× the stored energy' "$SIM_DOC"
  check "fiberFailure off: multiple of the stored energy" 5.74 0.005 "$(unconservative 2)" \
        'carrying 2.23× the nodal force and 5.74× the stored energy' "$SIM_DOC"
  check "fiberFailure off: extra punch work (%)" 25.3 0.05 "$(unconservative 3)" \
        '**25.3% more energy to open the same hole**' "$SIM_DOC"

  # **And the roadmap's restatement of the same table**, pinned so the two copies
  # cannot part company. That is not hypothetical here: the suite count was gated on
  # the README and loose on the roadmap once already, and went stale by two commits
  # in a spelling no grep aimed at the README would find. Every figure above appears
  # again in `06-roadmap.md`'s Phase 3 prose, wrapped differently, so each pointer is
  # a single line of that wrapping and nothing wider.
  for h in 'ram opening 125.6 m² of her side and the 163 m of longitudinal in it: the section' \
           'area lost goes 6.861% → **8.455%**, the second moment 5.429% → **6.590%**, the' \
           'hogging ultimate moment 5.388% → **6.747%** and the sagging one 12.369% →' \
           '`ε_p = 0.513` — **2.83× its own failure strain** — still at full section, carrying' \
           '2.23× the force and costing the ram **25.3% more energy to open the same hole**.'; do
    hint "$h" "$DOC"
  done

  # --- the clamped-zone control, out of the same run -------------------------------
  #
  # "74% out and 433× too stiff" is published on `README.md`, in the roadmap, and in
  # `02-simulation.md`'s table with the 625 671 N behind it. The 74% and the 1.47e-4
  # were printed all along. **The 433× was not**, and the only thing holding it was
  # the vacuity guard in `test_coupling.cpp` -- `clampedForce > 10.0 * reference`,
  # a bound forty-three times looser than the figure it stands next to, with its own
  # comment saying "Measured at 433x, asserted at 10."
  #
  # The guard is left where it is: it is a *vacuity* check, and its job is to fail
  # when an uncoupled zone would pass everything below it, not to pin a published
  # digit. Pinning the digit is this file's job, and the test now prints the ratio
  # so there is something to pin.
  clampedmn=$(printf '%s\n' "$suiteout" |
              sed -n 's/^ *clamped edge \([0-9.]*\) MN (\([0-9]*\)x the coupled answer).*/\1/p' | head -1)
  clampedx=$(printf '%s\n' "$suiteout" |
             sed -n 's/^ *clamped edge \([0-9.]*\) MN (\([0-9]*\)x the coupled answer).*/\2/p' | head -1)
  missm=$(printf '%s\n' "$suiteout" |
          sed -n 's/^ *control: the clamped zone misses the monolithic field by \([0-9.e+-]*\) m (\([0-9]*\)%)$/\1/p' | head -1)
  misspct=$(printf '%s\n' "$suiteout" |
            sed -n 's/^ *control: the clamped zone misses the monolithic field by \([0-9.e+-]*\) m (\([0-9]*\)%)$/\2/p' | head -1)
  check "clamped zone: punch reaction (N)" 625671 1 \
        "$(awk -v m="$clampedmn" 'BEGIN { printf "%.0f", m * 1e6 }')" \
        '| zone clamped, as before | 625 671 N (**433×**) | 1.47e-4 m out, 74% of peak |' \
        "$SIM_DOC"
  check "clamped zone: multiple of the coupled answer" 433 0 "$clampedx" \
        'clamped zone it replaces is 74% out and 433× too stiff. That it does not merely' \
        "$FRONT"
  check "clamped zone: field error at the perimeter (m)" 1.47e-4 1e-6 "$missm" \
        '| zone clamped, as before | 625 671 N (**433×**) | 1.47e-4 m out, 74% of peak |' \
        "$SIM_DOC"
  check "clamped zone: field error as a fraction of peak (%)" 74 0 "$misspct" \
        '| zone clamped, as before | 625 671 N (**433×**) | 1.47e-4 m out, 74% of peak |' \
        "$SIM_DOC"
  # The roadmap's third copy of the same pair.
  hint 'out and 433× too stiff. The *mesher* that was missing now exists —' "$DOC"

  # --- static condensation, out of the same run ------------------------------------
  #
  # **The assertion in the test cannot pin these, and that is the point.** They are
  # conditioning-limited residuals, so the in-test bound has to leave room for the
  # `-O3` engine build `verify.sh full` also compiles; it was `1e-6` of the peak
  # against a measurement 132x smaller, and even tightened to `1e-7` it spans the
  # whole decade the published figure moved across. Four documents said 2e-10 m
  # while the suite printed 2.34e-9 m, bit-stable, on every run.
  #
  # So the property is asserted in the test with headroom, and the *digit* is pinned
  # here, on the one build this file drives. That split is the honest one: a bound
  # loose enough to survive two optimisation levels cannot also be a bound tight
  # enough to notice a tenfold drift.
  cond() {  # $1 = 1 boundary, 2 interior, 3 peak
    printf '%s\n' "$suiteout" |
      sed -n "s/^ *boundary \([0-9.e+-]*\) m, interior \([0-9.e+-]*\) m, of a peak \([0-9.e+-]*\) m$/\\$1/p" |
      head -1
  }
  check "static condensation: interface error at zero modes (m)" 2.34e-9 5e-11 \
        "$(cond 1)" 'a 0.31 m** deflection' "$SIM_DOC"
  check "static condensation: interior error for an interface load (m)" 1.90e-9 5e-11 \
        "$(cond 2)" 'a 0.31 m** deflection' "$SIM_DOC"
  check "static condensation: the peak it is a fraction of (m)" 3.085e-1 1e-4 \
        "$(cond 3)" 'a 0.31 m** deflection' "$SIM_DOC"
  check "static condensation: interface error on the ferry's own plating (m)" 2.47e-9 5e-11 \
        "$(printf '%s\n' "$suiteout" |
            sed -n 's/^ *interface response \([0-9.e+-]*\) m of a peak .*/\1/p' | head -1)" \
        '**2.5 × 10⁻⁹ m of 0.027 m** on a real patch of the' "$SIM_DOC"
  # The three README copies of the same pair, so they cannot part company again.
  hint 'static condensation is exact at the interface for *any* load (2.3 × 10⁻⁹ m of a' "$FRONT"
  hint 'the interface for any load** (2.3 × 10⁻⁹ m of a 0.31 m deflection against an' "$DOC"

  # --- the coupled zone's own figures, same run ------------------------------------
  #
  # `0.737` and `1.06e-15` are published to three figures each; both are now pinned
  # in the test as well, so these are the cross-check that the *document* still says
  # what the suite measures.
  check "coupling: the monolithic edge follows the punch by" 0.737 0.0005 \
        "$(printf '%s\n' "$suiteout" |
            sed -n 's/^ *the edge follows \([0-9.]*\) of the punch.*/\1/p' | head -1)" \
        "the plate's own answer is **0.737**" "$SIM_DOC"
  check "coupling: the coupled field against the monolithic one (m)" 1.06e-15 1e-17 \
        "$(printf '%s\n' "$suiteout" |
            sed -n 's/^ *0 modes: .* field to \([0-9.e+-]*\) m of .*/\1/p' | head -1)" \
        '**1.06e-15 m** of a 2.0e-4 m peak' "$SIM_DOC"

  # --- the shoulder's spurious bending stiffness, same run -------------------------
  #
  # The front page's "319% too stiff in bending across this hull's shoulder". It was
  # printed and asserted nowhere: the only guard was on `worstNormalSpread`, a
  # different quantity, and the flat patch's exact 0.0000 turned that guard into
  # `> 5e-4` against a measured 0.1884.
  check "zone over the shoulder: offset/t" 0.1884 0.0002 \
        "$(printf '%s\n' "$suiteout" |
            sed -n 's/^ *flat of side offset.t .*over the shoulder \([0-9.]*\) (.*/\1/p' | head -1)" \
        '| across the shoulder at z ≈ 4.2 m | 0.188 |' "$SIM_DOC"
  check "zone over the shoulder: spurious bending stiffness (%)" 319 0.5 \
        "$(printf '%s\n' "$suiteout" |
            sed -n 's/^ *flat of side offset.t .*over the shoulder [0-9.]* (+\([0-9]*\)%).*/\1/p' | head -1)" \
        "319% too stiff in bending across this hull's shoulder" "$FRONT"

  # --- where Kawahara's eddy fit dies, out of the same run --------------------------
  #
  # A third of the declared validity box has no eddy damping at all, and the root is
  # inside the published Cb bound. Gated because it is a *limit* rather than a
  # result: if a future edit to the regression moves it, the paragraph in §2 that
  # states where the method stops working has to move with it.
  eddyroot() {  # $1 = 1 centre, 2 at B/d 2.5, 3 at B/d 4.5
    printf '%s\n' "$suiteout" |
      sed -n "s/^ *the eddy coefficient goes non-positive at Cb = \([0-9.]*\) (B.d [0-9.]*), \([0-9.]*\) at B.d [0-9.]* and \([0-9.]*\) at .*/\\$1/p" |
      head -1
  }
  check "the eddy fit's Cb root at B/d 3.5" 0.8438 0.0005 "$(eddyroot 1)" \
        'root at **0.8426 to 0.8461** depending on `B/d`' "$SIM_DOC"
  check "the eddy fit's Cb root at B/d 2.5" 0.8426 0.0005 "$(eddyroot 2)" \
        'root at **0.8426 to 0.8461** depending on `B/d`' "$SIM_DOC"
  check "the eddy fit's Cb root at B/d 4.5" 0.8461 0.0005 "$(eddyroot 3)" \
        'root at **0.8426 to 0.8461** depending on `B/d`' "$SIM_DOC"

  # --- what the tile budget is actually denominated in, out of the same run --------
  #
  # The estimator bills tiles per unit *volume*; the solver allocates them over the
  # compartment's *footprint*. On the ferry's forepeak those differ by 32.8x at 3 m3
  # and 8.2x at 12 m3 -- and the allocated count is the same 12 300 both times,
  # because both fills are thinner than one 0.2 m tile. Gated because the numbers
  # are the argument: a single figure would read as "the estimate is a bit low",
  # and the pair is what shows it is a function of the wrong variable.
  tiletbl() {  # $1 = 1 billed small, 2 alloc small, 3 billed large, 4 alloc large
    printf '%s\n' "$suiteout" |
      sed -n "s/^ *forepeak tiles, billed against allocated: 3 m3 \([0-9]*\) vs \([0-9]*\) (.*), 12 m3 \([0-9]*\) vs \([0-9]*\) (.*/\\$1/p" |
      head -1
  }
  check "forepeak at 3 m3: tiles the estimator bills" 375 0 "$(tiletbl 1)" \
        '3 m3   billed   375, allocated 12 300' engine/sim/water_promotion.hpp
  check "forepeak at 3 m3: tiles the solver allocates" 12300 0 "$(tiletbl 2)" \
        '3 m3   billed   375, allocated 12 300' engine/sim/water_promotion.hpp
  check "forepeak at 12 m3: tiles the estimator bills" 1500 0 "$(tiletbl 3)" \
        '12 m3  billed  1 500, allocated 12 300' engine/sim/water_promotion.hpp
  check "forepeak at 12 m3: tiles the solver allocates, unmoved" 12300 0 "$(tiletbl 4)" \
        '12 m3  billed  1 500, allocated 12 300' engine/sim/water_promotion.hpp

  # --- the spectral wave field, out of the same run --------------------------------
  #
  # Section 2's wave-field figures. `test_waves.cpp` printed **nothing** but its own
  # section header, so every one of these was published and produced by no run at
  # all -- the same condition the damping-ratio table was in. Printing them found
  # the document wrong immediately: see the rate correction below.
  #
  # **Anchored on words unique to these lines.** `test_ocean.cpp:1047` prints
  # "Hs 3 m Tp 9 s:" -- the very fixture the wave-field paragraph names -- and
  # `runOceanTests` runs *before* `runWaveTests`, so a grep on the fixture
  # description would take the ocean line. Worse, that block is inside
  # `#if defined(SHIPSIM_HAS_VULKAN) && !defined(__SANITIZE_THREAD__)`, so the
  # collision exists on a GPU box and vanishes on a headless one, and the gate would
  # read a different number on different machines without ever saying why.
  check "PM quantile: tabulated inverse against the closed form (relative)" 6e-15 1e-15 \
        "$(printf '%s\n' "$suiteout" | sed -n 's/.*PM quantile:.*worst \([0-9.e+-]*\) relative.*/\1/p')" \
        "reproduces it to 6 × 10⁻¹⁵ relative" "$SIM_DOC"
  # gamma is on the line because the edge depends on it: PM gives 2.985 and the
  # JONSWAP 3.3 the document names gives 2.683. It does *not* depend on Tp or Hs,
  # which is why this Tp = 13 s fixture answers for the Tp = 9 s one published.
  check "the open-topped bin starts at this multiple of omega_p" 2.7 0.05 \
        "$(printf '%s\n' "$suiteout" | sed -n 's/.*open bin at gamma 3.3 starts at \([0-9.]*\) omega_p.*/\1/p')" \
        "starts at 2.7ω_p and runs to infinity" "$SIM_DOC"
  check "PM mean period ratio T1/Tp" 0.77177 5e-6 \
        "$(printf '%s\n' "$suiteout" | sed -n 's/.*PM T1 = \([0-9.]*\) Tp.*/\1/p')" \
        "T1 = 0.77177 Tp" "$SIM_DOC"
  check "and how exactly the discretisation reproduces it (relative)" 2e-12 5e-13 \
        "$(printf '%s\n' "$suiteout" | sed -n 's/.*Tp, exact to \([0-9.e+-]*\) relative.*/\1/p')" \
        "2 × 10⁻¹² relative" "$SIM_DOC"

  # **The zero-crossing bias, and the rate the document had wrong.** It said the
  # bias halves as N doubles; its own three data points say otherwise, and printing
  # the sweep settled it -- 12 to 48 is a factor of four in N for a factor of two in
  # bias. All four fit 6.47/sqrt(N). Printed unsigned, because the document writes
  # them signed and no extractor in this file admits a leading `+`; the sign is
  # asserted in the test instead.
  t2bias() {
    printf '%s\n' "$suiteout" |
      sed -n "s/.*PM T2 bias by component count:.*N $1 \([0-9.]*\)%.*/\1/p"
  }
  check "PM zero-crossing bias at N = 12 (%)" 1.92 5e-3 "$(t2bias 12)" \
        "+1.92% at N = 12" "$SIM_DOC"
  check "at N = 48 (%)" 0.93 5e-3 "$(t2bias 48)" "+0.93% at" "$SIM_DOC"
  check "at N = 96 (%)" 0.66 5e-3 "$(t2bias 96)" "+0.66% at N = 96" "$SIM_DOC"
  check "at N = 384, a quarter the bias for thirty-two times the components (%)" 0.33 5e-3 \
        "$(t2bias 384)" "+0.33% at N = 384" "$SIM_DOC"

  # --- hull-to-hull contact, out of the same run ----------------------------------
  #
  # All eight rows of section 2's "What a ram looks like in numbers", printed by
  # `test_collision.cpp` on every run and gated by nothing. Two S-175-like hulls,
  # one struck on her port side 30 m forward of midship at 6 m/s, both floating
  # free with hydrostatics, damping and drag.
  #
  # The patch location is the row that matters most and reads least like a figure:
  # it is what the structural model consumes, and the paragraph under the table
  # claims the patch "lands within 5 cm of the station aimed at, on the port side,
  # inside the hull's depth". Three coordinates, gated at the centimetre.
  ram1=$(printf '%s\n' "$suiteout" | grep -m1 'ram: closing')
  ram2=$(printf '%s\n' "$suiteout" | grep -m1 'm penetration over')
  ram3=$(printf '%s\n' "$suiteout" | grep -m1 'patch at x =')

  check "ram: closing speed (m/s)" 6.0 0.05 \
        "$(printf '%s\n' "$ram1" | sed -n 's/.*closing \([0-9.]*\) m\/s.*/\1/p')" \
        "| closing speed | 6.0 m/s |" "$SIM_DOC"
  check "steps in contact" 162 0 \
        "$(printf '%s\n' "$ram1" | sed -n 's/.*m\/s, \([0-9]*\) steps.*/\1/p')" \
        "| 1.62 s (162 steps at dt = 10 ms) |" "$SIM_DOC"
  check "and the duration that is (s)" 1.62 5e-3 \
        "$(printf '%s\n' "$ram1" | sed -n 's/.*contact (\([0-9.]*\) s).*/\1/p')" \
        "| 1.62 s (162 steps at dt = 10 ms) |" "$SIM_DOC"
  check "peak normal force (MN)" 323 0.5 \
        "$(printf '%s\n' "$ram1" | sed -n 's/.*peak \([0-9.]*\) MN.*/\1/p')" \
        "| peak normal force | 323 MN |" "$SIM_DOC"

  check "penetration at peak (m)" 0.398 5e-4 \
        "$(printf '%s\n' "$ram2" | sed -n 's/ *\([0-9.]*\) m penetration.*/\1/p')" \
        "| penetration at peak | 0.398 m |" "$SIM_DOC"
  check "projected patch at peak (m2)" 11.4 0.05 \
        "$(printf '%s\n' "$ram2" | sed -n 's/.*over \([0-9.]*\) m2.*/\1/p')" \
        "| projected patch at peak | 11.4 m² |" "$SIM_DOC"
  check "mean contact pressure at peak (MPa)" 28.3 0.05 \
        "$(printf '%s\n' "$ram2" | sed -n 's/.*at \([0-9.]*\) MPa.*/\1/p')" \
        "| mean contact pressure at peak | 28.3 MPa |" "$SIM_DOC"
  check "energy absorbed (MJ)" 233 0.5 \
        "$(printf '%s\n' "$ram2" | sed -n 's/.*MPa, \([0-9.]*\) MJ.*/\1/p')" \
        "| energy absorbed | 233 MJ |" "$SIM_DOC"

  patch="| x = +30.0 m, y = +9.3 m, z = +8.3 m |"
  check "patch x in the struck ship's frame (m)" 30.0 0.05 \
        "$(printf '%s\n' "$ram3" | sed -n 's/.*x = +\([0-9.]*\) m.*/\1/p')" "$patch" "$SIM_DOC"
  check "patch y, which puts it on the port side (m)" 9.3 0.05 \
        "$(printf '%s\n' "$ram3" | sed -n 's/.*y = +\([0-9.]*\) m.*/\1/p')" "$patch" "$SIM_DOC"
  check "patch z, which puts it inside the hull's depth (m)" 8.3 0.05 \
        "$(printf '%s\n' "$ram3" | sed -n 's/.*z = +\([0-9.]*\) m.*/\1/p')" "$patch" "$SIM_DOC"

  # --- propulsion and manoeuvring, out of the same run ----------------------------
  #
  # All eight rows of section 7's "Measured behaviour" table. Six were printed by
  # `test_propulsion.cpp` already; the two bollard coefficients were computed and
  # then only *bracketed against the published B4-70 reference* -- 0.35 +- 0.03 and
  # 0.49 +- 0.04 -- so the model's own 0.347 and 0.493, which are what the table
  # publishes, appeared nowhere a reader or a gate could see them. A bracket around
  # a reference is not a check that the model still produces what the document says.
  #
  # The section's own caveat is left standing and is worth repeating: the MMG
  # coefficient set "was transcribed from the MMG standard-method literature and
  # has not been checked against a primary source in this worktree". These checks
  # say the model has not drifted from its own published figures. They say nothing
  # about whether those figures are right.
  boll=$(printf '%s\n' "$suiteout" | grep -m1 'bollard K_T')
  eta=$(printf '%s\n' "$suiteout" | grep -m1 'open-water peak eta')
  turn=$(printf '%s\n' "$suiteout" | grep -m1 '35 deg turn:')
  turn2=$(printf '%s\n' "$suiteout" | grep -m1 '20 deg: R =')

  check "bollard K_T at P/D 1.0, A_E/A_0 0.70" 0.347 5e-4 \
        "$(printf '%s\n' "$boll" | sed -n 's/.*bollard K_T \([0-9.]*\),.*/\1/p')" \
        "| 0.347 |" "$SIM_DOC"
  check "bollard 10 K_Q, same" 0.493 5e-4 \
        "$(printf '%s\n' "$boll" | sed -n 's/.*10 K_Q \([0-9.]*\),.*/\1/p')" \
        "| 0.493 |" "$SIM_DOC"
  check "zero-thrust advance ratio" 0.849 5e-4 \
        "$(printf '%s\n' "$boll" | sed -n 's/.*zero thrust at J \([0-9.]*\).*/\1/p')" \
        "| 0.849 |" "$SIM_DOC"
  check "peak open-water efficiency" 0.672 5e-4 \
        "$(printf '%s\n' "$eta" | sed -n 's/.*peak eta = \([0-9.]*\) at.*/\1/p')" \
        "| 0.672 at J = 0.710 |" "$SIM_DOC"
  check "and the advance ratio it peaks at" 0.710 5e-4 \
        "$(printf '%s\n' "$eta" | sed -n 's/.*at J = \([0-9.]*\) (zero.*/\1/p')" \
        "| 0.672 at J = 0.710 |" "$SIM_DOC"

  # The turning circle, which is the whole point of having a manoeuvring model.
  check "steady turning radius at 35 deg rudder (m)" 360 0.5 \
        "$(printf '%s\n' "$turn" | sed -n 's/.*R=\([0-9.]*\) m =.*/\1/p')" \
        "| 1.13 L (360 m) |" "$SIM_DOC"
  check "the same as a multiple of Lpp" 1.13 5e-3 \
        "$(printf '%s\n' "$turn" | sed -n 's/.*m = \([0-9.]*\) L.*/\1/p')" \
        "| 1.13 L (360 m) |" "$SIM_DOC"
  check "drift angle at 35 deg (deg)" 19.4 0.05 \
        "$(printf '%s\n' "$turn" | sed -n 's/.*drift=-\([0-9.]*\) deg.*/\1/p')" \
        "| 19.4° |" "$SIM_DOC"
  # Not printed as a percentage anywhere, so it is derived here from the two speeds
  # on the same line rather than from a second run. **The model's figure is 40.6 and
  # the table publishes 41**, so this check prints a different number from the one
  # in the document by design: the expectation is the published value and the
  # tolerance is the whole percent it was published at, which is the rule the rest
  # of this file follows.
  check "speed retained at 35 deg (% of approach)" 41 0.5 \
        "$(printf '%s\n' "$turn" |
           awk '{ for (i = 1; i <= NF; ++i) { if ($i == "approach") a = $(i+1)
                                              if ($i ~ /^U=/) { u = $i; sub(/^U=/, "", u) } }
                  if (a > 0) printf "%.1f", 100 * u / a }')" \
        "| 41 % |" "$SIM_DOC"
  check "steady turning radius at 20 deg (L)" 1.76 5e-3 \
        "$(printf '%s\n' "$turn2" | sed -n 's/.*20 deg: R = \([0-9.]*\) L.*/\1/p')" \
        "| 1.76 L / 2.73 L |" "$SIM_DOC"
  check "and at 10 deg (L)" 2.73 5e-3 \
        "$(printf '%s\n' "$turn2" | sed -n 's/.*10 deg: R = \([0-9.]*\) L.*/\1/p')" \
        "| 1.76 L / 2.73 L |" "$SIM_DOC"

  # --- strip-theory radiation, out of the same run --------------------------------
  #
  # Twelve of the thirteen rows of section 2's "What was measured" table, every one
  # of which `test_radiation.cpp` already printed and none of which was gated. They
  # all reproduce, so nothing here is a correction -- which is worth saying, because
  # the two tables gated just before this one did not.
  #
  # The runtime row is excluded: 0.53 microseconds per tick is a wall clock.
  radline() { printf '%s\n' "$suiteout" | grep -m1 -- "$1"; }

  semi=$(radline 'semicircle: Ca(inf)')
  check "heaving semicircle, Ca at infinite frequency" 1.0022 5e-5 \
        "$(printf '%s\n' "$semi" | sed -n 's/.*Ca(inf) \([0-9.]*\) .*/\1/p')" \
        "1.0022 × ρπa²/2" "$SIM_DOC"
  check "and its added-mass minimum" 0.5983 5e-5 \
        "$(printf '%s\n' "$semi" | sed -n 's/.*min Ca \([0-9.]*\) .*/\1/p')" \
        "**0.5983** at" "$SIM_DOC"
  check "at this nondimensional frequency" 0.894 5e-4 \
        "$(printf '%s\n' "$semi" | sed -n 's/.*omega sqrt(a\/g) \([0-9.]*\);.*/\1/p')" \
        "= 0.894\`, against Ursell" "$SIM_DOC"
  # The near-field/far-field energy balance -- the instrument that found the
  # radiation solver returning *negative* damping at four frequencies.
  check "energy residual at 40 panels" 4.5e-3 5e-5 \
        "$(printf '%s\n' "$semi" | sed -n 's/.*residual \([0-9.e+-]*\) at 40 panels.*/\1/p')" \
        "4.5 × 10⁻³ at 40 panels" "$SIM_DOC"
  check "and at 80" 2.2e-3 5e-5 \
        "$(printf '%s\n' "$semi" | sed -n 's/.*at 40 panels, \([0-9.e+-]*\) at 80.*/\1/p')" \
        "2.2 × 10⁻³ at 80" "$SIM_DOC"

  # Published as 3.9e-4 and measured 3.85e-4, which is exactly half a unit in the
  # last published place -- the policy tolerance would sit on the boundary, so this
  # one is a whole unit and says why.
  check "reciprocity A24 against A42 at 80 panels" 3.9e-4 1e-5 \
        "$(radline 'reciprocity (A24 vs A42)' | sed -n 's/.*A42) \([0-9.e+-]*\) at.*/\1/p')" \
        "3.9 × 10⁻⁴ at 80 panels" "$SIM_DOC"

  ogil=$(radline 'Ogilvie: B->K->B')
  check "Ogilvie round trip B to K to B, worst of peak" 6.3e-3 5e-5 \
        "$(printf '%s\n' "$ogil" | sed -n 's/.*worst \([0-9.e+-]*\) of peak.*/\1/p')" \
        "6.3 × 10⁻³ of peak" "$SIM_DOC"
  check "A_inf, rigid lid against Ogilvie (%)" 0.57 5e-3 \
        "$(printf '%s\n' "$ogil" | sed -n 's/.*lid), \([0-9.]*\)% apart.*/\1/p')" \
        "0.57% apart" "$SIM_DOC"
  check "and how far the Ogilvie value varies across omega (%)" 2.3 5e-2 \
        "$(printf '%s\n' "$ogil" | sed -n 's/.*apart, \([0-9.]*\)% spread.*/\1/p')" \
        "varies 2.3% across" "$SIM_DOC"

  trans=$(radline 'closed-form transforms:')
  check "closed-form transform of B" 2.8e-7 5e-9 \
        "$(printf '%s\n' "$trans" | sed -n 's/.*error \([0-9.e+-]*\) (B),.*/\1/p')" \
        "2.8 × 10⁻⁷ (B)" "$SIM_DOC"
  check "and of A" 4.4e-6 5e-8 \
        "$(printf '%s\n' "$trans" | sed -n 's/.*, \([0-9.e+-]*\) (A).*/\1/p')" \
        "4.4 × 10⁻⁶ (A)" "$SIM_DOC"

  mem=$(radline 'memory: K falls to')
  check "memory: K falls to 1% of peak at (s)" 20.3 5e-2 \
        "$(printf '%s\n' "$mem" | sed -n 's/.*peak at \([0-9.]*\) s.*/\1/p')" \
        "**20.3 s**" "$SIM_DOC"
  check "and to 0.1% at (s)" 56.6 5e-2 \
        "$(printf '%s\n' "$mem" | sed -n 's/.*and 0.1% at \([0-9.]*\) s.*/\1/p')" \
        "0.1% at 56.6 s" "$SIM_DOC"

  ss=$(radline "state space on the ferry's K33")
  check "state space, relative rms" 7.6e-2 5e-4 \
        "$(printf '%s\n' "$ss" | sed -n 's/.*relative rms \([0-9.e+-]*\),.*/\1/p')" \
        "7.6% relative RMS" "$SIM_DOC"
  check "and its peak error as a percentage of K(0)" 2.8 5e-2 \
        "$(printf '%s\n' "$ss" | sed -n 's/.*peak error \([0-9.]*\)% of K(0).*/\1/p')" \
        "**2.8% of K(0)**" "$SIM_DOC"
  check "Prony on a planted 4-pole signal, relative rms" 3.0e-9 5e-11 \
        "$(radline 'Prony on a planted' | sed -n 's/.*relative rms \([0-9.e+-]*\).*/\1/p')" \
        "3.0 × 10⁻⁹ relative RMS" "$SIM_DOC"

  scale=$(radline 'geometric scaling by 2.5:')
  check "geometric scaling, added mass exact to" 7.5e-15 5e-17 \
        "$(printf '%s\n' "$scale" | sed -n 's/.*added mass to \([0-9.e+-]*\),.*/\1/p')" \
        "7.5 × 10⁻¹⁵" "$SIM_DOC"
  check "and damping to" 1.9e-12 5e-14 \
        "$(printf '%s\n' "$scale" | sed -n 's/.*damping to \([0-9.e+-]*\) relative.*/\1/p')" \
        "1.9 × 10⁻¹²" "$SIM_DOC"

  # --- the barge RAO sweep, out of the same run -----------------------------------
  #
  # **Section 2 says an RAO is "the one place this simulator can be checked against
  # the outside world", and then published a table nothing produced.** The
  # assertions around it checked that heave collapses relative to its neighbours
  # and that the asymptotes hold -- never what the cells say. Re-derived, row
  # three's `lambda/L` turned out to be upside down (0.97 is 60/61.62) and row
  # two's phase read +9 against a measured -38.4.
  #
  # The wavelength column is not a measurement -- `lambda = 2 pi g / omega^2` -- but
  # it is gated because the notch claim two paragraphs down is stated in metres and
  # was otherwise two divisions away from anything printed.
  raorow() {
    printf '%s\n' "$suiteout" |
      awk -v w="$1" -v c="$2" '/barge sweep, 60 x 16 m box/ { f = 1 }
                               /^ *$/ { if (f && seen) f = 0 }
                               f && NF == 8 && $1 == w { seen = 1; print $c; exit }'
  }
  for row in "0.25 985.87 16.431 1.034 1.073 -9.5" \
             "0.70 125.75 2.096 0.899 1.352 -38.4" \
             "1.00 61.62 1.027 0.036 0.517 -74.7" \
             "1.45 29.31 0.488 0.019 0.076 -129.4" \
             "2.20 12.73 0.212 0.015 0.004 -152.9"; do
    set -- $row
    w=$1
    ptr="| $w | $2 | $3 |"
    check "barge RAO omega $w: wavelength (m)" "$2" 5e-3 "$(raorow "$w" 2)" "$ptr" "$SIM_DOC"
    check "barge RAO omega $w: heave" "$4" 5e-4 "$(raorow "$w" 4)" "$ptr" "$SIM_DOC"
    check "barge RAO omega $w: pitch" "$5" 5e-4 "$(raorow "$w" 5)" "$ptr" "$SIM_DOC"
    check "barge RAO omega $w: heave phase (deg)" "$6" 0.05 "$(raorow "$w" 6)" "$ptr" "$SIM_DOC"
  done
  # The two claims the section rests on, which are closed forms rather than table
  # cells: the Froude-Krylov sinc zeros land at a wavelength of L and L/2, and
  # pitch follows the wave *slope* so it lags the elevation by a quarter period.
  check "the first sinc notch, at a wavelength of the ship's length (m)" 61.62 5e-3 \
        "$(raorow 1.00 2)" "at 61.62 m and" "$SIM_DOC"
  check "and the second, at half of it (m)" 29.31 5e-3 "$(raorow 1.45 2)" \
        "29.31 m against a 60 m ship" "$SIM_DOC"
  check "pitch lags the surface elevation (deg)" -96.3 0.05 "$(raorow 0.25 7)" \
        "lags the surface elevation by 96.3°" "$SIM_DOC"

  # --- the hull-form convergence table, out of the same run ----------------------
  #
  # **`docs/05` had zero gated figures and this script's own header explains why,
  # wrongly.** It reasons at :78-79 that `05`'s digits "are the fields of a file
  # format `test_shipfile.cpp` already parses" -- true of §1 and false of §2, which
  # is nine measured convergence percentages. That is the same document-level
  # generalisation the header already confesses to about `01`: "A figure is worth
  # gating here when a *tool* produces it -- which is a question about the line,
  # not about the file."
  #
  # Nothing produced these until now. `test_hullform.cpp` swept the same waterline
  # counts, computed the middle row as an absolute residual, asserted only that
  # refining improves it, and printed none of the numbers. Re-derived, they came
  # back two to five times smaller than what had been published for as long as the
  # table existed.
  hullrow() {
    printf '%s\n' "$suiteout" |
      awk -v st="$1" -v col="$2" \
          '/S-175 block-coefficient error/ { f = 1 }
           /worst \|LCB - asked\|/ { f = 0 }
           f && NF == 4 && $1 == st { print $col; exit }'
  }
  MOD_DOC=docs/05-data-modding-validation.md
  # Tolerance is the rounding the table publishes at, three decimals of a percent.
  for row in "21 0.195 0.104 0.180" "41 0.351 0.053 0.022" "161 0.392 0.094 0.019"; do
    set -- $row
    st=$1
    ptr="| $st | $2% | $3% |"
    check "S-175 Cb error, $st stations, 11 waterlines (%)" "$2" 5e-4 "$(hullrow "$st" 2)" \
          "$ptr" "$MOD_DOC"
    check "S-175 Cb error, $st stations, 21 waterlines (%)" "$3" 5e-4 "$(hullrow "$st" 3)" \
          "$ptr" "$MOD_DOC"
    check "S-175 Cb error, $st stations, 41 waterlines (%)" "$4" 5e-4 "$(hullrow "$st" 4)" \
          "$ptr" "$MOD_DOC"
  done
  # The two LCB figures are reported apart because they are different claims: the
  # area curve is refined by stations and not by waterlines, and the document used
  # to publish the waterline-sweep number as holding "at every resolution".
  lcb=$(printf '%s\n' "$suiteout" |
        sed -n 's/.*worst |LCB - asked|: \([0-9.e+-]*\) of Lpp over the waterline sweep, \([0-9.e+-]*\) .*/\1 \2/p' |
        tail -1)
  check "worst |LCB - asked| over the waterline sweep (of Lpp)" 4.0e-5 5e-7 \
        "$(printf '%s\n' "$lcb" | cut -d' ' -f1)" "worst 4.0 × 10⁻⁵ of Lpp over that sweep" "$MOD_DOC"
  check "and over all nine meshes, where stations move it" 2.3e-4 5e-6 \
        "$(printf '%s\n' "$lcb" | cut -d' ' -f2)" "2.3 × 10⁻⁴ over all nine meshes" "$MOD_DOC"

  # --- Ikeda's viscous roll damping, out of the same run -------------------------
  #
  # Two hulls and two tables: the ro-pax figures come from `test_roll_damping.cpp`
  # and the damping-ratio table from `test_rao.cpp`, both printing into the stdout
  # captured above.
  #
  # **The damping-ratio table had no producer at all until it was given one.**
  # Grepping its own published figures -- 0.0259, 0.0323, 0.0470, 0.0810, 0.00062,
  # 0.00193 -- across `tests/`, `tools/` and `engine/` returned nothing for any of
  # them. It was a hand-made table, and re-deriving it meant writing a driver
  # against the library. Printed at last, it came back 1.5% from what had been
  # published, on a hull where the `B0` correction in the same commit moves it by
  # 0.2%. So it had drifted from the code by something that was not that change,
  # and nothing in the repository could have said so.
  #
  # Anchored on `$2 == "m/s"` and on the numeric amplitude rather than on line
  # offsets: both tables sit in the middle of the suite's output, and a checker
  # that misparses is worse than no checker.
  #
  # **`exit` on the first match, and it is load-bearing.** Without it the flag `f`
  # opens a range that never closes, so every *later* line in 200 000 lines of
  # suite output whose first field is numerically 5 or 10 prints as well -- and
  # `check` then hands awk a multi-line `$actual`, whose numeric coercion takes the
  # leading token and passes. The figure was right by luck and the parse was
  # wrong; the tell was stray fragments (`m3:`, `x`) interleaved with the ticks.
  # This is the collision the file already warns about at :335 for `fem_spike`,
  # where `awk '$1 == 1'` matches the *section title* `1.` because "1." is a
  # numeric string.
  ropax() {
    printf '%s\n' "$suiteout" |
      awk -v deg="$1" -v col="$2" \
          '/ro-pax 170x25x6.5/ { f = 1 }
           f && $2 == "m/s" && $1 == 0.0 && $3 == deg { print $col; exit }' | tr -d '%'
  }
  zetarow() {
    printf '%s\n' "$suiteout" |
      awk -v deg="$1" -v col="$2" \
          '/roll damping ratio at omega_n/ { f = 1 }
           f && $1 == deg { print $col; exit }' | tr -d '%'
  }
  hulls=$(printf '%s\n' "$suiteout" |
          sed -n 's/.* bare \([0-9.]*\)  keeled \([0-9.]*\).*/\1 \2/p' | tail -1)

  check "the reference ro-pax's bare-hull viscous B44hat" 0.0081 5e-5 \
        "$(printf '%s\n' "$hulls" | cut -d' ' -f1)" "0.0081 bare" "$SIM_DOC"
  check "and the same hull with bilge keels, at omegahat 1.1 and 20 deg" 0.0429 5e-5 \
        "$(printf '%s\n' "$hulls" | cut -d' ' -f2)" "0.0429 with keels" "$SIM_DOC"
  # The figure `Ship::zetaRoll = 0.08` was anchored to, written down in three
  # places: here, again at the stand-in paragraph, and in `ship.cpp`.
  check "fraction of critical at the natural frequency and 10 deg" 0.061 5e-4 \
        "$(ropax 10.0 8)" "6.1% of critical" "$SIM_DOC"
  check "the bilge keels' share of the ro-pax total there (%)" 80.2 0.5 \
        "$(ropax 10.0 16)" "Bilge keels are 80% of the total" "$SIM_DOC"

  check "ferry damping ratio with keels at 2.5 deg" 0.0263 5e-5 "$(zetarow 2.5 2)" \
        "| 2.5° | 0.0263 |" "$SIM_DOC"
  check "at 5 deg" 0.0327 5e-5 "$(zetarow 5.0 2)" "| 5° | 0.0327 |" "$SIM_DOC"
  check "at 10 deg" 0.0476 5e-5 "$(zetarow 10.0 2)" "| 10° | 0.0476 |" "$SIM_DOC"
  check "at 20 deg" 0.0820 5e-5 "$(zetarow 20.0 2)" "| 20° | 0.0820 |" "$SIM_DOC"
  # The bare column is the control: it never reaches `B0` at all, so a change to
  # the bilge-keel pressure model must leave it exactly alone. It is also what
  # showed the published table had been derived some other way, because it had
  # moved too.
  check "and bare at 10 deg, which no bilge-keel term can touch" 0.00196 5e-6 \
        "$(zetarow 10.0 3)" "| 10° | 0.0476 | 0.00196 |" "$SIM_DOC"
else
  echo "  - shipsim_tests not built, skipping the README's check count"
  skipped="$skipped shipsim_tests"
fi

# The three numbers the SURVIVED -> LOST verdict actually rests on. Until `--gm-detail`
# existed these were the least checkable figures on the front page: nothing in the
# repository printed the angle GM was sampled at, so re-deriving them meant writing a
# C++ driver against `libshipsim_engine`. They are gated here because an argument that
# only its author can reproduce is not evidence.
#
# **The angle is the one that had drifted, and it looked derived.** The README published
# 6.8e-4 rad, taken off the deck's *nominal* 100 x 19 m. The deck the ship carries is
# 1868.4 m2 with a mean breadth of 18.684 m, and the run measures 7.07e-4 -- 4% out. A
# bounding box would say 20.00 m, the deck at its widest, which is the one width nothing
# in this calculation wants; `test_ship.cpp` asserts the mean and not the box for exactly
# that reason.
#
# The fixed-angle figure is deliberately a *control*: it is the current code asked for the
# answer it used to give, not a number kept from before the change. That is what makes it
# reproducible rather than historical.
if [ -x "$SHIPSIM" ]; then
  gmd=$("$SHIPSIM" --scenario=full --duration=1800 --gm-detail 2>/dev/null)
  detail() { printf '%s\n' "$gmd" | sed -n "s/^gm-detail: $1 \(.*\)$/\1/p"; }

  check "GM at a fixed +-0.03 rad, the answer that used to be published" \
        1.3799 5e-4 "$(detail gm_at_fixed_0.03rad_m)" "reports +1.38 m" "$FRONT"
  check "the angle the layer pockets at (rad)" \
        7.0655e-4 5e-8 "$(detail pockets_at_rad)" "7.07e-4 rad" "$FRONT"
  check "fraction of the deck's free surface a fixed sample sees" \
        0.06341 5e-5 "$(detail fixed_sample_sees_frac)" "sees **6%**" "$FRONT"
  check "the layer's depth (mm)" \
        6.6005 5e-4 "$(awk -v d="$(detail layer_depth_m)" 'BEGIN{print d*1000}')" \
        "a **6.6 mm** layer" "$FRONT"
  check "the vehicle deck's mean breadth (m)" \
        18.684 5e-4 "$(detail layer_breadth_m)" "mean** breadth of 18.684 m" "$FRONT"

  # **`06-roadmap.md` publishes the same measurement for all three scenarios and
  # none of it was gated.** The `full` row was already known to be wrong: this
  # block has checked the README's copy since `--gm-detail` existed, and
  # `README.md` names 6.8e-4 rad in its own prose as the superseded value, while
  # the roadmap still carried 6.8e-4 along with the two GM figures that moved with
  # it. One document held the correction and the other held the figure it
  # corrected, for as long as the gate was pointed at one file.
  #
  # `doors` had drifted 2.8%. **`none` had drifted by a factor of 34** -- 119 t to
  # 4025 t, a 6.9 cm layer to 2.34 m -- which is not drift so much as a different
  # ship: it predates authoring the mid wing tanks. The verdict column never
  # changed, so nothing about the row looked wrong.
  #
  # Two extra runs, ~150 s. The alternative was adding `--gm-detail` to the three
  # scenario runs this script already makes, which would have been free; it is not
  # done because those runs feed four positional parsers apiece and a misparse
  # bought with saved seconds is the worst trade on offer here.
  for scen in none doors; do
    sd=$("$SHIPSIM" --scenario=$scen --duration=1800 --gm-detail 2>/dev/null)
    sdetail() { printf '%s\n' "$sd" | sed -n "s/^gm-detail: $1 \(.*\)$/\1/p"; }
    case "$scen" in
      none)  want_w=4025; want_l=2.34; want_p=0.245; want_f=-2.234; want_c=-2.234 ;;
      doors) want_w=3950; want_l=2.29; want_p=0.241; want_f=-2.475; want_c=-2.475 ;;
    esac
    # **The pointer carries the figures, not just the row label.** `| \`none\`` on
    # its own also matches the superseded table in the correction block right
    # below, so it would go on passing while the live row was edited -- a pointer
    # narrower than the thing it stands for, which this file already warns about
    # for the dedup key at :154.
    case "$scen" in
      none)  row='| 4025 t | 2.34 m | 0.245 rad  | −2.234 m |' ;;
      doors) row='| 3950 t | 2.29 m | 0.241 rad  | −2.475 m |' ;;
    esac
    check "'$scen': water on the vehicle deck (t)" "$want_w" 1 "$(sdetail layer_water_t)" \
          "$row" "$DOC"
    check "'$scen': the layer's depth (m)" "$want_l" 5e-3 "$(sdetail layer_depth_m)" \
          "$row" "$DOC"
    check "'$scen': the angle it pockets at (rad)" "$want_p" 5e-4 \
          "$(sdetail pockets_at_rad)" "$row" "$DOC"
    # These two are the point of the row: with the deck this deeply flooded there
    # is no free surface left to pocket, so the fixed-angle sample and the
    # converged value agree to under a millimetre. On `full` they do not agree at
    # all, which is the finding the table exists to carry.
    check "'$scen': GM at a fixed +-0.03 rad (m)" "$want_f" 5e-4 \
          "$(sdetail gm_at_fixed_0.03rad_m)" "$row" "$DOC"
    check "'$scen': and the converged GM (m)" "$want_c" 5e-4 "$(sdetail gm_converged_m)" \
          "$row" "$DOC"
  done
  # `full`'s roadmap copy, from the run already made above.
  fullrow='| 6.6 mm | 7.07e-4 rad | **+1.38 m** | **−3.23 m** |'
  check "'full': the roadmap's copy of the layer depth (mm)" 6.6005 5e-4 \
        "$(awk -v d="$(detail layer_depth_m)" 'BEGIN{print d*1000}')" "$fullrow" "$DOC"
  check "'full': the roadmap's copy of the pocketing angle (rad)" 7.0655e-4 5e-8 \
        "$(detail pockets_at_rad)" "$fullrow" "$DOC"
  check "'full': the roadmap's copy of the converged GM (m)" -3.2350 5e-4 \
        "$(detail gm_converged_m)" "$fullrow" "$DOC"
else
  echo "  - shipsim not built, skipping the README's GM-detail figures"
  skipped="$skipped shipsim-gmdetail"
fi

# **`ram_view` is a Vulkan target, so a machine with no device never builds it --
# and this used to `exit 0` at that point.** That took the section and the smoke
# blocks below down with it, neither of which needs a GPU, and it threw away any
# failure already counted: a red gate exiting 0. The README figures above are what
# made it matter, because they are pure CPU and a machine with no device is exactly
# where nobody is watching the output.
if [ -x "$RAM" ]; then
  for h in "opens 3.4 m" "63 bays torn, 107.7" "amidships she goes over" \
           "at the quarter she takes" "she lolls to 63"; do
    hint "$h"
  done

  # Every figure below comes off one of two lines, and *only* off that line. The
  # first version of this script matched `GM` anywhere in the output and picked up
  # the intact ship's GM from the banner as well as the flooded one from the
  # outcome, reporting "2.00" and "-1.62" together. A checker that misparses is
  # worse than no checker, so each reader is anchored to its own line.
  damage_line() { printf '%s\n' "$1" | grep '^damage'; }
  outcome_line() { printf '%s\n' "$1" | grep '^outcome'; }

  # --- damage, which needs no flooding -------------------------------------------
  quick6=$("$RAM" --speed=6 --duration=1 2>&1)
  quick15=$("$RAM" --speed=1.5 --duration=1 2>&1)

  area6=$(damage_line "$quick6" | sed -n 's/.*torn, \([0-9.]*\) m2 of hole.*/\1/p')
  area15=$(damage_line "$quick15" | sed -n 's/.*torn, \([0-9.]*\) m2 of hole.*/\1/p')
  torn=$(damage_line "$quick6" | sed -n 's/.*over [0-9]* bays, \([0-9]*\) torn.*/\1/p')
  check "hole at 1.5 m/s (m2)" 3.4 0.3 "$area15" "opens 3.4 m"
  check "hole at 6 m/s (m2)" 107.7 0.6 "$area6" "107.7 m² of hole"
  check "bays torn at 6 m/s" 63 1 "$torn" "63 bays torn, 107.7"

  # --- the headline outcome, which needs the full 900 s ---------------------------
  #
  # Three full runs, which is the same budget the first version of this script
  # spent, but aimed at the thing that turned out to matter: the amidships strike
  # sits just above a capsize threshold, so a figure alone says very little. 5.5 m/s
  # and 6 m/s straddle it, and asserting both is what makes the pair a finding
  # rather than two numbers.
  out=$("$RAM" --speed=6 2>&1)
  water=$(outcome_line "$out" | sed -n 's/.*-- \([0-9]*\) t of water.*/\1/p')
  heel=$(outcome_line "$out" | sed -n 's/.*heel \(-\{0,1\}[0-9.]*\) deg.*/\1/p')
  check "floodwater amidships at 6 m/s (t)" 16613 400 "$water" "amidships she goes over"
  checks=$((checks + 1))
  if [ -n "$heel" ] && awk -v h="$heel" 'BEGIN { if (h < 0) h = -h; exit !(h > 90) }'; then
    printf '  %s✓%s and she goes right over: heel %s deg, past 90\n' "$green" "$off" "$heel"
  else
    printf '  %s✗%s she no longer capsizes at 6 m/s: heel %s deg\n' "$red" "$off" "$heel"
    printf '      %s%s describes this as a capsize under "amidships she goes over"%s\n' \
           "$dim" "$DOC" "$off"
    fails=$((fails + 1))
  fi

  # Below the threshold, where she lolls instead. The pair is the point: a change
  # that moved the threshold past 6 m/s, or below 5.5, would leave both tonnages
  # looking reasonable on their own and still have altered what the ship does.
  lo=$("$RAM" --speed=5.5 2>&1)
  lowater=$(outcome_line "$lo" | sed -n 's/.*-- \([0-9]*\) t of water.*/\1/p')
  loheel=$(outcome_line "$lo" | sed -n 's/.*heel \(-\{0,1\}[0-9.]*\) deg.*/\1/p')
  check "floodwater amidships at 5.5 m/s (t)" 7415 300 "$lowater" "she lolls to 63"
  checks=$((checks + 1))
  if [ -n "$loheel" ] && awk -v h="$loheel" 'BEGIN { if (h < 0) h = -h; exit !(h < 90) }'; then
    printf '  %s✓%s and at 5.5 m/s she lolls instead: heel %s deg, short of 90\n' \
           "$green" "$off" "$loheel"
  else
    printf '  %s✗%s the capsize threshold has moved below 5.5 m/s: heel %s deg\n' \
           "$red" "$off" "$loheel"
    printf '      %sthe threshold between 5.5 and 5.75 m/s is a finding in %s%s\n' \
           "$dim" "$DOC" "$off"
    fails=$((fails + 1))
  fi

  # --- and that being struck at the quarter is worse, and the other way -----------
  out=$("$RAM" --speed=6 --aim=-30 2>&1)
  qwater=$(outcome_line "$out" | sed -n 's/.*-- \([0-9]*\) t of water.*/\1/p')
  qheel=$(outcome_line "$out" | sed -n 's/.*heel \(-\{0,1\}[0-9.]*\) deg.*/\1/p')
  check "floodwater at the quarter (t)" 7323 150 "$qwater" "at the quarter she takes"
  check "heel at the quarter (deg)"    -47.8 1.0 "$qheel" "at the quarter she takes"

  # --- and what the *renderer* makes of that damage, which is 03's figure not 06's
  #
  # `docs/03-renderer-audio.md`'s damage-cost paragraph is the only place the
  # renderer is measured on a real ship, and none of it was gated. It had drifted,
  # in the one column that cannot move unless the ship does: 27 631 scene triangles
  # for 28 019. The interior drawn behind the shell *is* her compartment set, so
  # authoring the mid wing tanks -- the same gap that had 41% of a ram amidships
  # tearing open onto nothing -- added 388 triangles here, and the front page's
  # `16 compartments` was the identical drift one document in.
  #
  # ~1 s: `--duration=1` runs almost no flooding, and every count below is fixed
  # before the first step -- the refinement is of the undeformed plating and the
  # scene is rebuilt identically each frame. Both were checked against
  # `--duration=120` and against `--frames=3`, and neither moves.
  #
  # **The two millisecond figures in that paragraph are deliberately not here.**
  # `3.7 ms` of rebuild and `0.15 ms` of GPU are wall clocks; on a box with other
  # work on it the rebuild reads 3.5, which is contention and not a change in the
  # mesh, and a gate that goes red for that teaches people to ignore it.
  drawn=$("$RAM" --speed=4.0 --duration=1 --frames=1 --out=/tmp 2>&1)
  drawn_line() { printf '%s\n' "$drawn" | grep '^drawn'; }
  src=$(drawn_line | sed -n 's/.*her \([0-9]*\) hull triangles.*/\1/p')
  refined=$(drawn_line | sed -n 's/.*refine to \([0-9]*\),.*/\1/p')
  cutout=$(drawn_line | sed -n 's/.*of which \([0-9]*\) are cut out.*/\1/p')
  holem2=$(drawn_line | sed -n 's/.*cut out (\([0-9.]*\) m2 of hole).*/\1/p')
  check "the ferry's hull triangles before refinement" 1196 0 "$src" \
        "Her 1 196 hull triangles" "$RENDER_DOC"
  check "what they refine to at 4 m/s" 7568 0 "$refined" \
        "refine to 7 568, of which 3 345 are cut out" "$RENDER_DOC"
  check "how many of those are cut out" 3345 0 "$cutout" \
        "refine to 7 568, of which 3 345 are cut out" "$RENDER_DOC"
  check "the drawn hole at 4 m/s (m2)" 44.8 0.05 "$holem2" \
        "44.8 m² of hole at 4 m/s" "$RENDER_DOC"
  # The scene total needs a device, because it is counted off a frame that was
  # actually submitted. The four counts above do not and are checked either way --
  # `damage.cpp` contains no Vulkan on purpose, and this keeps that purpose.
  if printf '%s\n' "$drawn" | grep -q 'no usable GPU'; then
    echo "  - no Vulkan device, skipping ram_view's drawn scene total"
    skipped="$skipped ram_view-scene"
  else
    scene=$(printf '%s\n' "$drawn" | sed -n 's/^ *[0-9]* vertices, \([0-9]*\) triangles, gpu.*/\1/p')
    # **`06-roadmap.md` publishes the same count and nothing re-derived it.** It
    # read 27 631 -- the value `03-renderer-audio.md` names, in prose, as the one
    # superseded when the mid wing tanks were authored and added 388 triangles to
    # this frame. One document carried the correction and the other carried the
    # figure it corrected, for as long as the gate was pointed at one spelling in
    # one file. Derived from the same literal so the two cannot part company, the
    # way the suite count at :660 is.
    sceneWant=28019
    check "triangles in the drawn scene: hull, interior and sea" "$sceneWant" 0 "$scene" \
          "28 019 triangles" "$RENDER_DOC"
    hint "$(printf '%s %s triangles at 1280' \
            "${sceneWant%???}" "${sceneWant#${sceneWant%???}}")" "$DOC"
  fi
else
  echo "  - ram_view not built, skipping the collision-milestone figures"
  skipped="$skipped ram_view"
fi

# --- what the water promoter sees, and the control that says it can see anything ---
#
# The roadmap's Phase 5 entry rests on a *pair* of runs, and the pair is the only
# form in which either half means anything. "The flooding scenario promotes
# nothing" is the claim that activating the promoter is safe for every figure
# above, and on its own it is equally consistent with a criterion that is simply
# dead -- which is the more likely of the two failures and the one that would
# never announce itself. So the beam-sea control is gated beside it and has to
# keep promoting.
#
# The same shape as the capsize-threshold pair and the `--no-interface-ties`
# control: a figure that can only fail in the safe direction is not gated at all.
#
# ~25 s: 900 s of flooding at dt=0.01 and 300 s of beam sea.
if [ -x "$WATER" ]; then
  for h in "0.0042 rad/s of roll rate" "0.1525 rad/s, 3× over the threshold"; do
    hint "$h" "$DOC"
  done
  # The flooding half: she must promote nothing, and the roll rate she gets
  # nowhere near is the reason. Gating the promotion count alone would pass on a
  # criterion that had been switched off.
  wflood=$("$WATER" --duration=900 2>&1)
  wroll=$(printf '%s\n' "$wflood" | sed -n 's/^ *peak roll rate seen *\([0-9.]*\).*/\1/p')
  wprom=$(printf '%s\n' "$wflood" | sed -n 's/^ *promotions *\([0-9]*\).*/\1/p')
  check "the flooding ferry's peak roll rate (rad/s)" 0.0042 0.0001 "$wroll" \
        "0.0042 rad/s of roll rate" "$DOC"
  check "and what it promotes" 0 0 "$wprom" "0.0042 rad/s of roll rate" "$DOC"

  # The control: the same ship in a beam sea at her own roll period. If this ever
  # reads zero promotions the block above has stopped being evidence of anything.
  wsea=$("$WATER" --wave=2.0 --duration=300 2>&1)
  searoll=$(printf '%s\n' "$wsea" | sed -n 's/^ *peak roll rate seen *\([0-9.]*\).*/\1/p')
  seaprom=$(printf '%s\n' "$wsea" | sed -n 's/^ *promotions *\([0-9]*\).*/\1/p')
  check "the beam-sea control's peak roll rate (rad/s)" 0.1525 0.002 "$searoll" \
        "0.1525 rad/s, 3× over the threshold" "$DOC"
  check "and that it does promote" 7 0 "$seaprom" \
        "0.1525 rad/s, 3× over the threshold" "$DOC"

  # **The cost, which is the figure that decides whether this tier can be turned
  # on at all.** It replaced a `5.0` that had carried "estimate, will be measured"
  # since it was written, and a tolerance here is not a rounding allowance: the
  # claim is that a promoted compartment costs three orders of magnitude more than
  # realtime, so what has to stay true is the order of magnitude. Gated at the 1 m³
  # end because that is the cheapest case and therefore the strongest form of the
  # claim -- if even 1 m³ is 28x realtime, nothing larger is affordable either.
  #
  # ~90 s: five sizes up to 100 m³, six steps each.
  hint "| 1 m³ | 65 650 | 525 | 279 | **27.9** | 5.6× |" "$DOC"
  wcost=$("$WATER" --cost 2>&1)
  onem3=$(printf '%s\n' "$wcost" | sed -n 's/^ *1 m3 *[0-9]* *[0-9]* *[0-9.]* *\([0-9.]*\).*/\1/p')
  check "one cubic metre, core-seconds per simulated second" 27.9 4.0 "$onem3" \
        "| 1 m³ | 65 650 | 525 | 279 | **27.9** | 5.6× |" "$DOC"
else
  echo "  - water_probe not built, skipping the promoter's figures"
  skipped="$skipped water_probe"
fi

# --- the gas tier's budget, which admits exactly one compartment -------------------
#
# `gas_probe` is the gas tier's `water_probe`, and it exists because `GasPromoter`
# had been exercised by one test file and nothing else -- so every compartment it
# had ever judged was built by a test to have the property under test. Run against
# the real ferry it reproduced the water tier's defect in a different unit: the
# machinery spaces are 7 670 cells each against a 12 000-cell budget, so the first
# to promote takes 64% and the second is refused however hard it qualifies.
#
# **The refusal count is the figure, not the promotion count.** One promotion on
# its own is consistent with a criterion that is simply quiet; 35 refusals in the
# same run is what says the tier is budget-bound rather than idle, and it is the
# number that would move if either limit were retuned. Both are gated, as a pair.
#
# ~20 s: 600 s of fire at 2 s steps.
if [ -x "$GAS" ]; then
  for h in "7 670 cells each and the" "**35 times**"; do
    hint "$h" "$GASDOC"
  done
  gasrun=$("$GAS" --duration=600 2>&1)
  gasrc=$?
  # `gas_probe` gained a failure path when `fire::StepResult::incomplete` landed:
  # it exits 1 rather than publishing figures indexed by a clock that stopped
  # keeping up with the model. Capturing its output without its status would parse
  # those figures as though nothing had happened -- the tool would print its
  # refusal and this gate would read straight past it, which is the same defect
  # the exit was added to close, one level up.
  checks=$((checks + 1))
  if [ "$gasrc" -eq 0 ]; then
    printf '  %s✓%s gas_probe published its figures off a clock that kept up\n' \
           "$green" "$off"
  else
    printf '  %s✗%s gas_probe exited %s — the figures below are indexed by a clock that fell behind\n' \
           "$red" "$off" "$gasrc"
    fails=$((fails + 1))
  fi
  gasprom=$(printf '%s\n' "$gasrun" | sed -n 's/^ *promotions *\([0-9]*\).*/\1/p')
  gasref=$(printf '%s\n' "$gasrun" | sed -n 's/^ *reviews that refused on budget *\([0-9]*\).*/\1/p')
  gascells=$(printf '%s\n' "$gasrun" | sed -n 's/^ *engine_room_s *\([0-9]*\) *[0-9.]*.*/\1/p' | tail -1)
  check "the gas tier promotes one compartment" 1 0 "$gasprom" \
        "7 670 cells each and the" "$GASDOC"
  check "and refuses the other on budget, repeatedly" 35 2 "$gasref" \
        "**35 times**" "$GASDOC"
  check "a machinery space's resolved cell count" 7670 0 "$gascells" \
        "7 670 cells each and the" "$GASDOC"
else
  echo "  - gas_probe not built, skipping the gas tier's budget figures"
  skipped="$skipped gas_probe"
fi

# --- the section mesher's reach ---------------------------------------------------
#
# `section_probe` joined the gate earlier today, after it published a 3.46 Hz shell
# frequency that nothing had ever re-run. **Gating a tool stops the tool rotting; it
# does not stop the documents quoting it.** That half stayed open, and it showed:
# `section.hpp` §4's resolution table had already drifted from the program it names
# -- 1.73075 where it published 1.73122 -- with nothing to catch it. So the figures
# the reach scan publishes are checked here the same way `ram_view`'s are.
#
# The reach is the right thing to check rather than the section properties: it is
# what the whole-ship claim rests on, it is quoted in two documents and a header, and
# it is an integer count rather than a tolerance argument.
if [ -x "$SECTION" ]; then
  for h in "49 of 49 two-bay windows" "120.0 m of 120.0 m"; do hint "$h" "$SECTION_DOC"; done
  scan=$("$SECTION" --scan=2 2>&1)
  meshed=$(printf '%s\n' "$scan" | sed -n 's/^ *\([0-9]*\) of \([0-9]*\) windows mesh.*/\1/p')
  windows=$(printf '%s\n' "$scan" | sed -n 's/^ *\([0-9]*\) of \([0-9]*\) windows mesh.*/\2/p')
  onepiece=$(printf '%s\n' "$scan" | sed -n 's/^ *\([0-9]*\) of [0-9]* are also a single connected piece.*/\1/p')
  reach=$(printf '%s\n' "$scan" | sed -n 's/^ *reach: \([0-9.]*\) m of.*/\1/p')
  check "windows that mesh and solve" 49 0 "$meshed" "49 of 49 two-bay windows" "$SECTION_DOC"
  check "windows in the scan"         49 0 "$windows" "49 of 49 two-bay windows" "$SECTION_DOC"
  check "windows in one piece"        46 0 "$onepiece" "46 in one piece" "$SECTION_DOC"
  check "reach along the hull (m)"    120.0 0.05 "$reach" "120.0 m of 120.0 m" "$SECTION_DOC"

  # --- and what an interior cut plane costs the junctions, which is torsion or nothing
  #
  # Four figures over two runs, and the pair is the point exactly as the capsize
  # threshold above is: the tied row on its own would pass on a chain that had never
  # lost anything, so the negative control -- `--no-interface-ties`, the same tool with
  # the line ties off -- is checked alongside it and has to still lose 3.3% of `GJ`.
  #
  # `GJ` and not `EA`: `EA` is exact to 1e-9 with the cut planes open *and* closed, so
  # a figure taken off it would score a do-nothing implementation best. That is the
  # trap the section it points at is written around.
  for h in "72.0 m** (= one piece) | **−0.095%" "9.6 m of 134.4 | −3.334%"; do
    hint "$h" "$SECTION_DOC"
  done
  summary() { printf '%s\n' "$1" | grep '^chain summary: ' | grep "junctions=1 lines=$2"; }
  field() { printf '%s\n' "$1" | sed -n "s/.* $2=\\(-\\{0,1\\}[0-9.e+-]*\\).*/\\1/p"; }
  tied=$(summary "$("$SECTION" --chain=2 --from=-7.2 --to=2.4 --sweep=0 --no-reduce 2>&1)" 1)
  open=$(summary "$("$SECTION" --chain=2 --from=-7.2 --to=2.4 --sweep=0 --no-reduce \
                    --no-interface-ties 2>&1)" 0)
  check "ferry chain of 2, junction edge joined (m)" 72.0 0.05 "$(field "$tied" tiedEdges)" \
        "72.0 m** (= one piece)" "$SECTION_DOC"
  check "ferry chain of 2, GJ against one piece" -0.00095 0.0004 "$(field "$tied" GJrel)" \
        "72.0 m** (= one piece) | **−0.095%" "$SECTION_DOC"
  check "ferry chain of 2, planes the two sides disagreed about" 0 0 \
        "$(field "$tied" disagree)" "worstPlaneTieDisagreement\` is **0.0" "$SECTION_DOC"
  check "ferry chain of 2 with the line ties off, edge joined (m)" 9.6 0.05 \
        "$(field "$open" tiedEdges)" "9.6 m of 134.4" "$SECTION_DOC"
  check "ferry chain of 2 with the line ties off, GJ" -0.03334 0.0005 \
        "$(field "$open" GJrel)" "9.6 m of 134.4 | −3.334%" "$SECTION_DOC"

  # --- the modes, the tie and the refinement, which the flags above switch off ----
  #
  # **This gate ran `section_probe` and passed the flags that disable the stages
  # producing the figures.** `--sweep=0` skips the refinement sweep and
  # `--no-reduce` skips the substructure, so four document sections that the tool
  # names itself as the source of -- the Tier-1 mesher, the junction tie, the
  # Craig-Bampton reduction and the resolution table, 400-odd digit-carrying lines
  # between them -- had between zero and one gated figure each. The tool's own
  # header says "Every *ship-scale* figure `docs/02-simulation.md` section 3
  # publishes about the section mesher comes out of this program".
  #
  # It cost a figure while nobody was looking. `02-simulation.md` carries a
  # **Correction** recording that the shell-alone frequency read 3.4600 Hz until the
  # tool was re-run, ending "this is the fourth time that has cost this repository a
  # published number". The value that correction installed, 1.6999 Hz, was written
  # into seven places and gated in none, and it has since moved to 1.6980 -- a fifth
  # time, to the figure the correction itself put there.
  #
  # The exact inertia count the document cites as what makes these frequencies
  # trustworthy cannot catch that: it brackets the shell mode between 10.5619 and
  # 10.7753 rad/s, a window 2% wide, and both values sit inside it. A bracket that
  # proves a mode exists is not a check that it has not moved.
  #
  # 103 s, and every physics figure in it is bit-identical across runs -- two runs
  # diff only in `meshed in`, `solve s` and `reduce s` -- so these carry tolerance 0
  # except where the document rounds. The three wall-clock columns are not gated,
  # for the reason the rest of this file gives.
  sweep=$("$SECTION" --sweep=2 2>&1)
  cbmode() {
    printf '%s\n' "$sweep" |
      awk -v p="$1" 'index($0, "  " p) == 1 { getline
                     if (match($0, /= [0-9.]+ Hz/)) { print substr($0, RSTART + 2, RLENGTH - 5); exit } }'
  }
  resrow() {
    printf '%s\n' "$sweep" |
      awk -v s="$1" -v c="$2" '/=== resolution:/ { f = 1 } /^ok$/ { f = 0 }
                               f && NF == 8 && $1 == s { print $c; exit }'
  }
  tierow() {
    printf '%s\n' "$sweep" |
      awk -v r="$1" -v c="$2" '/=== the junction tie:/ { f = 1 } /=== resolution:/ { f = 0 }
                               f && $1 == r { print $c; exit }'
  }

  # The five fixed-interface frequencies, and the two that matter most are the ones
  # the paragraph turns on: the decks alone are the softest thing in the section,
  # and tying the junctions lifts the whole above *both* pieces.
  check "shell alone, first fixed-interface mode (Hz)" 1.6980 5e-5 "$(cbmode 'shell only')" \
        "the shell on its own is **1.6980 Hz**" "$SECTION_DOC"
  check "decks alone, first fixed-interface mode (Hz)" 0.7785 5e-5 "$(cbmode 'decks only')" \
        "decks *on their own* are" "$SECTION_DOC"
  check "the untied whole, which the decks set" 0.7785 5e-5 "$(cbmode 'whole, untied')" \
        "first fixed-interface mode is **0.7785 Hz**" "$SECTION_DOC"
  check "and tied, above both pieces" 2.3026 5e-5 "$(cbmode 'whole, tied')" \
        "the same section is **2.3026 Hz**" "$SECTION_DOC"

  # --- the hull girder on a wave, which nothing ran -------------------------------
  #
  # `section_probe --wave` is the only producer of §2's Tier-0-against-Tier-1 table,
  # and no script passed the flag. Ten published figures, including the two the
  # section argues from: the deck's *mean* agreeing with the beam to 1.5% says the
  # moment is arriving, and the worst standing 42% above it says the field is not a
  # beam. Either alone is a different claim.
  #
  # 162 s, of which nearly all is a chain of ten built and solved over the whole
  # hull. Its two siblings are measured and declined: `--profile=2` is 20 s but
  # produces the along-the-length window table rather than the planes-apart grid the
  # document publishes, so gating it would pin numbers §2 does not quote; and
  # `--whole=5` is **487 s** for the chain-against-monolith table, which is the same
  # cost bracket as the 1 044 s monolithic solve the document already declines by
  # name. Both re-derive correctly by hand -- `--whole=5` reproduces EA 1.47255e+11,
  # z_na 4.97283, EI 2.38679e+12 and all four deltas -- so what is missing there is
  # budget, not evidence.
  wave=$("$SECTION" --wave --sweep=0 --no-reduce 2>&1)
  wrow() {  # $1 = the row's leading words, $2 = Tier 0 (1) or Tier 1 (2)
    printf '%s\n' "$wave" | grep "^  $1  *" | sed "s/^  $1  *//" | awk -v c="$2" '{ print $c }'
  }
  check "hull girder at x=6: Tier 0 deck fibre (MPa)" 82.057 0.005 \
        "$(wrow 'deck fibre' 1)" '| deck fibre | 82.06 MPa | **118.42 MPa** |' "$SECTION_DOC"
  check "hull girder at x=6: Tier 1 deck fibre (MPa)" 118.415 0.005 \
        "$(wrow 'deck fibre' 2)" '| deck fibre | 82.06 MPa | **118.42 MPa** |' "$SECTION_DOC"
  check "hull girder at x=6: the deck's own mean (MPa)" 83.255 0.005 \
        "$(printf '%s\n' "$wave" | sed -n 's/^ *its mean over the deck  *(one number)  *\([0-9.]*\).*/\1/p')" \
        "| the deck's own mean | (one number) | 83.26 MPa |" "$SECTION_DOC"
  check "hull girder at x=6: Tier 0 keel fibre (MPa)" -66.516 0.005 \
        "$(wrow 'keel fibre' 1)" '| keel fibre | −66.52 MPa | −88.49 MPa |' "$SECTION_DOC"
  check "hull girder at x=6: Tier 1 keel fibre (MPa)" -88.489 0.005 \
        "$(wrow 'keel fibre' 2)" '| keel fibre | −66.52 MPa | −88.49 MPa |' "$SECTION_DOC"
  check "hull girder at x=6: Tier 0 neutral axis (m)" 6.7132 0.0005 \
        "$(wrow 'neutral axis' 1)" '| neutral axis | 6.7132 m | 6.7961 m |' "$SECTION_DOC"
  check "hull girder at x=6: Tier 1 neutral axis (m)" 6.7961 0.0005 \
        "$(wrow 'neutral axis' 2)" '| neutral axis | 6.7132 m | 6.7961 m |' "$SECTION_DOC"
  # The two halves of the argument. Tier 0 must be exactly zero -- a beam cannot
  # carry a field a beam does not have -- and Tier 1's 8.06 is the residual that
  # says the field is not one.
  check "hull girder at x=6: what a beam cannot carry, Tier 0" 0 0 \
        "$(wrow 'what a beam cannot carry' 1)" \
        '| what a beam cannot carry | 0 | **8.06 MPa rms** |' "$SECTION_DOC"
  check "hull girder at x=6: what a beam cannot carry, Tier 1 (MPa)" 8.057 0.005 \
        "$(wrow 'what a beam cannot carry' 2)" \
        '| what a beam cannot carry | 0 | **8.06 MPa rms** |' "$SECTION_DOC"
  # And the six rigid restraints, which are determinate and so must carry nothing.
  check "the six restraints carry (N)" 3847 5 \
        "$(printf '%s\n' "$wave" |
           sed -n 's/.*six restraints carry \([0-9.e+]*\) N against.*/\1/p' | head -1)" \
        'They are determinate, so on a balanced load they' "$SECTION_DOC"

  # **The orderings that lost.** `1 382` is quoted in three places as what the node
  # numbering is worth, and until the chooser was made to keep its losing scores no
  # invocation could produce it: it kept the narrowest of three candidates and threw
  # the other two away. It could not be checked and it could not be refuted -- a
  # *tied* candidate at 1 382 is arithmetically impossible, since the tied winner is
  # gated at 1 520 and the winner is the smallest. Measured, it is the untied
  # y-fastest ordering exactly.
  #
  # Read on its own prefix rather than on "cut"/"tied", because `tierow` selects on
  # the first field and takes the first match -- a second line starting with the same
  # word would have been read as the data row and quietly broken eight figures above.
  ordering() {  # $1 = cut|tied, $2 = which column
    printf '%s\n' "$sweep" |
      awk -v r="orderings-$1" -v c="$2" '$1 == r { print $c; exit }' | tr -d ',;'
  }
  check "untied hold: x-fastest ordering" 2528 0 "$(ordering cut 3)" \
        'half-bandwidth from 146 to **1 382**' "$SECTION_DOC"
  check "untied hold: y-fastest ordering, the 1382 three documents quote" 1382 0 \
        "$(ordering cut 5)" 'half-bandwidth from 146 to **1 382**' "$SECTION_DOC"
  check "untied hold: RCM, which is what it ships with" 146 0 "$(ordering cut 7)" \
        'half-bandwidth from 146 to **1 382**' "$SECTION_DOC"
  check "tied hold: best ordering without RCM" 2540 0 "$(ordering tied 11)" \
        '(2 528 x-fastest, 1 382 y-fastest,' "$SECTION_DOC"

  # The junction tie, cut against tied. `GJ` is the figure the tie exists for.
  check "cut section: half-bandwidth" 146 0 "$(tierow cut 2)" "| cut" "$SECTION_DOC"
  check "cut section: components" 7 0 "$(tierow cut 3)" "| cut" "$SECTION_DOC"
  check "cut section: GJ (N m^2)" 3.6164e12 5e8 "$(tierow cut 7)" "3.6164e12" "$SECTION_DOC"
  check "tied section: half-bandwidth" 1520 0 "$(tierow tied 2)" "| tied" "$SECTION_DOC"
  check "tied section: components" 1 0 "$(tierow tied 3)" "| tied" "$SECTION_DOC"
  check "tied section: edge joined (m)" 309.6 0.05 "$(tierow tied 4)" "309.6" "$SECTION_DOC"
  check "tied section: A_eff (m^2)" 1.73266 5e-6 "$(tierow tied 5)" "**1.73266**" "$SECTION_DOC"
  check "tied section: GJ (N m^2)" 5.2387e12 5e8 "$(tierow tied 7)" "5.2387e12" "$SECTION_DOC"

  # The resolution table, whose own Correction block records that it had drifted
  # from the program it names and that "no gate re-runs" it. Subdivisions 3 and 4
  # are 13 s and 38 s of solve on top and are left out; a drifted mesher shows in
  # the first row.
  check "resolution sub 1: elements" 2068 0 "$(resrow 1 2)" "| 1 | 2 068 |" "$SECTION_DOC"
  check "resolution sub 1: A_eff (m^2)" 1.72945 5e-6 "$(resrow 1 4)" "| 1 | 2 068 |" "$SECTION_DOC"
  check "resolution sub 1: z_na (m)" 6.86238 5e-6 "$(resrow 1 5)" "| 1 | 2 068 |" "$SECTION_DOC"
  check "resolution sub 1: I_eff (m^4)" 43.8169 5e-5 "$(resrow 1 6)" "| 1 | 2 068 |" "$SECTION_DOC"
  check "resolution sub 2: elements" 8272 0 "$(resrow 2 2)" "| 2 | 8 272 |" "$SECTION_DOC"
  check "resolution sub 2: A_eff (m^2)" 1.73081 5e-6 "$(resrow 2 4)" "| 2 | 8 272 |" "$SECTION_DOC"
  check "resolution sub 2: z_na (m)" 6.85933 5e-6 "$(resrow 2 5)" "| 2 | 8 272 |" "$SECTION_DOC"
  check "resolution sub 2: I_eff (m^4)" 43.8598 5e-5 "$(resrow 2 6)" "| 2 | 8 272 |" "$SECTION_DOC"
else
  echo "  - section_probe not built, skipping the reach figures"
  skipped="$skipped section_probe"
fi

# --- the ocean cascade, which is the largest ungated block in 03 -------------------
#
# `docs/03-renderer-audio.md` publishes the cascade's shape as a table -- nine
# levels, what each carries, what each costs in vertices -- and until now not one
# cell of it was checked. It is the section that most wants checking, because the
# whole argument is that **reach is exponential in the number of levels while cost
# is linear in it**, and the evidence for that is a set of integers: the level
# count, the ring vertex count that is the same at every level, and the component
# count falling to zero on the rings that resolve nothing.
#
# **Every parameter of the experiment is passed explicitly, including the ones that
# are already the defaults.** That is the §8 lesson kept: a probe that derives part
# of its own experiment is a probe whose figures cannot be compared against the
# table they came from. `--hs` sets the shortest component, the shortest component
# sets the cell size, the cell size sets the cell count -- so a default that moved
# would move every figure below while the document went on naming Hs 4 m.
#
# **The wall clocks in that table are not gated and cannot be.** `46 ms` of
# displacement comes back between 41 and 52 on this box depending on what else is
# running -- a 25% spread with the mesh bit-identical -- and `2.5 ms` of upload
# spreads further. Vertices and reach survive contention; milliseconds do not.
#
# ~5 s at `--frames=1`: the counts below are set before the loop and do not depend
# on it, but the level table and the vertex total are printed from a frame that was
# actually built, so one frame is the floor. Verified identical at `--frames=2`.
SEAWAY=${SEAWAY:-./build/seaway_view}
if [ -x "$SEAWAY" ]; then
  sea=$("$SEAWAY" --out=/tmp --frames=1 --ship=s175 --hs=4 --tp=9 --heading=180 --revs=4.6 2>&1)
  casc() { printf '%s\n' "$sea" | grep '^ *ocean cascade:'; }
  levels=$(casc | sed -n 's/^ *ocean cascade: \([0-9]*\) levels.*/\1/p')
  cells=$(casc | sed -n 's/^ *ocean cascade: [0-9]* levels, \([0-9]*\) cells.*/\1/p')
  eyeup=$(casc | sed -n 's/.*for an eye \([0-9]*\) m up.*/\1/p')
  reachkm=$(casc | sed -n 's/.*reaching \([0-9]*\) km .*/\1/p')
  # The cascade line is printed before any device is touched, so these four hold on
  # a machine with no Vulkan as well as on this one -- which is the half of the
  # section that is pure arithmetic over the camera and the spectrum.
  check "seaway_view's cascade levels"        9   0 "$levels"  "Nine levels of 212 cells" "$RENDER_DOC"
  check "seaway_view's cells per level"       212 0 "$cells"   "Nine levels of 212 cells" "$RENDER_DOC"
  check "the eye height the reach is set by (m)" 69 0.5 "$eyeup" \
        "camera 69 m up at 1280" "$RENDER_DOC"
  check "the cascade's reach (km)"            67  0.5 "$reachkm" "| 67 km | 316 729 |" "$RENDER_DOC"

  if printf '%s\n' "$sea" | grep -q 'no usable GPU'; then
    echo "  - no Vulkan device, skipping the per-level cascade table"
    skipped="$skipped seaway_view-levels"
  else
    # `[cell components vertices]` per level, so the n-th bracket is field 2n when
    # the line is split on brackets. Counted forward from the level rather than back
    # from the end, for the reason the §8 columns are: the tail of this line is nine
    # identical-looking rings and the interesting one is the first.
    level() { printf '%s\n' "$1" | grep '^ *cascade levels' |
              awk -F'[][]' -v n="$2" '{ print $(2 * n) }' | awk -v c="$3" '{ print $c }'; }
    comp_row="| components carried | 128 | 128 | 120 | 96 | 16 | **0** |"
    vert_row="| vertices | 45 369 | 33 920 | 33 920 | 33 920 | 33 920 | 33 920 each |"
    # **The components column is the finding and the vertices column is the guard.**
    # Levels 5 and 6 dropping to 16 and then 0 is what makes 60 km of sea nearly
    # free; the vertices column is flat by construction, so a ring that quietly
    # stopped being a ring shows up there and nowhere else.
    check "level 1: components carried" 128 0 "$(level "$sea" 1 2)" "$comp_row" "$RENDER_DOC"
    check "level 2: components carried" 128 0 "$(level "$sea" 2 2)" "$comp_row" "$RENDER_DOC"
    check "level 3: components carried" 120 0 "$(level "$sea" 3 2)" "$comp_row" "$RENDER_DOC"
    check "level 4: components carried"  96 0 "$(level "$sea" 4 2)" "$comp_row" "$RENDER_DOC"
    check "level 5: components carried"  16 0 "$(level "$sea" 5 2)" "$comp_row" "$RENDER_DOC"
    check "level 6: components carried"   0 0 "$(level "$sea" 6 2)" "$comp_row" "$RENDER_DOC"
    check "level 1: vertices" 45369 0 "$(level "$sea" 1 3)" "$vert_row" "$RENDER_DOC"
    check "level 2: vertices" 33920 0 "$(level "$sea" 2 3)" "$vert_row" "$RENDER_DOC"
    check "level 3: vertices" 33920 0 "$(level "$sea" 3 3)" "$vert_row" "$RENDER_DOC"
    check "level 4: vertices" 33920 0 "$(level "$sea" 4 3)" "$vert_row" "$RENDER_DOC"
    check "level 5: vertices" 33920 0 "$(level "$sea" 5 3)" "$vert_row" "$RENDER_DOC"
    check "level 6: vertices" 33920 0 "$(level "$sea" 6 3)" "$vert_row" "$RENDER_DOC"
    verts=$(printf '%s\n' "$sea" | sed -n 's/^ *best frame: \([0-9]*\) ocean vertices.*/\1/p')
    check "ocean vertices in the drawn cascade" 316729 0 "$verts" \
          "| 67 km | 316 729 |" "$RENDER_DOC"
  fi
else
  echo "  - seaway_view not built, skipping the ocean cascade figures"
  skipped="$skipped seaway_view"
fi

# --- the volumetric fire and smoke figures ---------------------------------------
#
# `docs/03-renderer-audio.md`'s fire section publishes a table of what a two-zone
# fire looks like and three findings taken off it, and the findings are the fragile
# part: "this fire has no fire in it" is a *negative* result, and a negative result
# is exactly what a well-meaning change to an exposure constant would quietly
# reverse. So the glow threshold and the peak layer temperature are checked as a
# pair, on either side of each other, rather than as two numbers.
#
# **The table those findings rest on is now checked cell by cell, and so is the
# counterfactual beside them.** Six figures of this document were gated when 155
# were gated in total, against 217 of its lines carrying a digit -- and the six were
# the findings' headline numbers rather than any of the measurements underneath
# them. Two things had drifted in that gap, neither of which any test could see:
# `ram_view`'s scene triangle count below, and this section's own "10 MW reaches
# 862 K and the layer does then glow", which is 876 K and does not.
#
# ~12 s: the 4 MW run is the one `verify.sh` already makes, and the two `--power`
# runs cost under a second each. **They are read for their output and not for their
# exit status, deliberately.**
#
# **That limitation has since been fixed in the tool, and this note is kept because
# the reason it was true is the interesting part.** `smoke_view` used to fail its own
# checks above its design fire: it bisected one threshold -- the temperature at which
# emission first puts a single byte of red on the screen, 834 K -- and then asserted
# against a completely different one, `red > blue + 8`, red *dominance* against a
# sky-blue background. Between them sat a real state the tool had no branch for: a
# layer that is emitting and is not red-dominant. A 10 MW fire lands exactly there, so
# the tool asserted a glow it had never claimed it would draw and exited non-zero on a
# correct picture. At 12 MW it failed the opposite way -- "the smoke reached the frame"
# counted only pixels the layer *blacked out*, and a glowing layer is not black, so a
# frame with 16 309 lit pixels was reported as never reached. The tool now bisects the
# dominance threshold too (935 K), carries the middle branch, and counts both of the
# things its own comment says a layer can do to a frame.
#
# The runs are still parsed rather than trusted, because an acceptance contract and a
# figure are different claims and this block only wants the figure.
SMOKE=${SMOKE:-./build/smoke_view}
if [ -x "$SMOKE" ]; then
  smoke=$("$SMOKE" --out=/tmp --frames=8 --duration=600 2>&1)
  if printf '%s\n' "$smoke" | grep -q 'no usable GPU'; then
    echo "  - no Vulkan device, skipping the fire and smoke figures"
    skipped="$skipped smoke_view"
  else
    # The last two are quoted by the *relation* blocks below rather than by any
    # `check`, and that is exactly the case the header of this file is about: a
    # pointer only a failure prints is a pointer nobody reads until it is wrong.
    # `"This fire has no fire in it"` had been printed by the glow relation since it
    # was written, against nothing.
    for h in "at **834 K**" "peaks at **531 K**" "reaches 511" \
             "84.1% of the engine" "This fire has no fire in it" \
             "passes 10 before t = 100 s"; do
      hint "$h" "$RENDER_DOC"
    done
    # Each reader anchored to its own line, for the reason the damage figures are:
    # `531` also appears in the table two rows above the one being read.
    #
    # A cell of the casualty table, addressed by the frame's own file name and by
    # column number counted forward: 1 frame, 2 t, 3 Q, 4 T_u, 5 z_i, 6 k, 7 tau,
    # 8 visibility, 9 emission. The file name is what makes a row addressable at
    # all -- every other field repeats down the table, and `t` is the very thing
    # being asserted on two of these rows.
    cell() { printf '%s\n' "$1" | grep "$2" | awk '{ print $'"$3"' }'; }
    glow=$(printf '%s\n' "$smoke" | sed -n 's/.*one byte of red on the screen at \([0-9]*\) K.*/\1/p')
    peak=$(printf '%s\n' "$smoke" | sed -n 's/.*this one peaked at \([0-9]*\) K.*/\1/p')
    check "the glow threshold (K)" 834 1 "$glow" "at **834 K**" "$RENDER_DOC"
    check "the ferry fire's peak layer temperature (K)" 531 2 "$peak" \
          "peaks at **531 K**" "$RENDER_DOC"
    checks=$((checks + 1))
    if [ -n "$glow" ] && [ -n "$peak" ] && [ "$peak" -lt "$glow" ]; then
      printf '  %s✓%s and the finding still holds: %s K of gas against a %s K threshold,'\
' so there is no glow to draw\n' "$green" "$off" "$peak" "$glow"
    else
      printf '  %s✗%s the ferry fire now reaches the glow threshold (%s K vs %s K)\n' \
             "$red" "$off" "$peak" "$glow"
      printf '      %s%s says "This fire has no fire in it"%s\n' "$dim" "$RENDER_DOC" "$off"
      fails=$((fails + 1))
    fi
    check "optical depth across the room at t = 600 s" 511 3 "$(cell "$smoke" smoke_07.png 7)" \
          "reaches 511" "$RENDER_DOC"
    check "the layer's extinction at t = 600 s (1/m)" 19.90 0.1 "$(cell "$smoke" smoke_07.png 6)" \
          "| 19.90 | 511" "$RENDER_DOC"
    plan=$(printf '%s\n' "$smoke" | sed -n 's/.*\([0-9][0-9]\.[0-9][0-9]\)% of the bounding box.*/\1/p')
    check "the drawn prism, % of the bounding box in plan" 84.11 0.05 "$plan" \
          "84.1% of the engine" "$RENDER_DOC"

    # --- and the casualty table itself, cell by cell --------------------------
    #
    # **One pointer per row, and it is the row verbatim**, the arrangement the
    # README's trapped-air table above uses and for the same reason: `check`
    # compares the tool against a constant in this file, so a row edited in the
    # *document* alone is caught by its pointer or it is not caught at all. Every
    # number in a row is inside that row's pointer, so there is no cell here that a
    # doc-only edit can move quietly.
    #
    # The t column is gated too, on the two interior rows, because it is not
    # decoration: the sample instants are `duration * frame / (frames - 1)`, and a
    # figure read at a different instant is the mistake the §8 block above exists
    # to prevent. Gating 171 and 343 pins the schedule the rest of the row is read
    # on, in the document, next to the numbers.
    row0="| 0   | 0.00 | 288 | 7.00 | 0.00  | 0    | clear |"
    row171="| 171 | 1.38 | 343 | 3.43 | 2.34  | 60   | 1.28  |"
    row343="| 343 | 4.00 | 523 | 3.05 | 11.22 | 288  | 0.267 |"
    row600="| 600 | 4.00 | 531 | 3.13 | 19.90 | 511  | 0.151 |"
    # t = 0 is the ambient control: the room before anything happens, at the
    # compartment's own height and with no extinction at all. It is the row that
    # cannot move unless the ship or the atmosphere does.
    check "t=0: the layer temperature (K)"      288  0.5   "$(cell "$smoke" smoke_00.png 4)" "$row0"   "$RENDER_DOC"
    check "t=0: the interface height (m)"       7.00 0.005 "$(cell "$smoke" smoke_00.png 5)" "$row0"   "$RENDER_DOC"
    check "t=171: the sample instant (s)"       171  0.5   "$(cell "$smoke" smoke_02.png 2)" "$row171" "$RENDER_DOC"
    check "t=171: heat release (MW)"            1.38 0.005 "$(cell "$smoke" smoke_02.png 3)" "$row171" "$RENDER_DOC"
    check "t=171: the layer temperature (K)"    343  0.5   "$(cell "$smoke" smoke_02.png 4)" "$row171" "$RENDER_DOC"
    check "t=171: the interface height (m)"     3.43 0.005 "$(cell "$smoke" smoke_02.png 5)" "$row171" "$RENDER_DOC"
    check "t=171: extinction (1/m)"             2.34 0.005 "$(cell "$smoke" smoke_02.png 6)" "$row171" "$RENDER_DOC"
    check "t=171: optical depth across the room" 60  0.5   "$(cell "$smoke" smoke_02.png 7)" "$row171" "$RENDER_DOC"
    check "t=171: visibility 3/k (m)"           1.28 0.005 "$(cell "$smoke" smoke_02.png 8)" "$row171" "$RENDER_DOC"
    check "t=343: the sample instant (s)"       343  0.5   "$(cell "$smoke" smoke_04.png 2)" "$row343" "$RENDER_DOC"
    check "t=343: heat release (MW)"            4.00 0.005 "$(cell "$smoke" smoke_04.png 3)" "$row343" "$RENDER_DOC"
    check "t=343: the layer temperature (K)"    523  0.5   "$(cell "$smoke" smoke_04.png 4)" "$row343" "$RENDER_DOC"
    check "t=343: the interface height (m)"     3.05 0.005 "$(cell "$smoke" smoke_04.png 5)" "$row343" "$RENDER_DOC"
    check "t=343: extinction (1/m)"             11.22 0.005 "$(cell "$smoke" smoke_04.png 6)" "$row343" "$RENDER_DOC"
    check "t=343: optical depth across the room" 288 0.5   "$(cell "$smoke" smoke_04.png 7)" "$row343" "$RENDER_DOC"
    check "t=343: visibility 3/k (m)"           0.267 0.0005 "$(cell "$smoke" smoke_04.png 8)" "$row343" "$RENDER_DOC"
    check "t=600: the interface height (m)"     3.13 0.005 "$(cell "$smoke" smoke_07.png 5)" "$row600" "$RENDER_DOC"
    check "t=600: visibility 3/k (m)"           0.151 0.0005 "$(cell "$smoke" smoke_07.png 8)" "$row600" "$RENDER_DOC"

    # --- the three findings' own numbers, which are not in the table -----------
    #
    # The descent and the recovery are the whole of "the layer does not descend
    # monotonically", and neither was gated. That finding has been wrong once
    # already in the other direction -- the first version of the test asserted a
    # monotone descent -- so it is the one on this page most worth holding to a
    # measurement. Gating both numbers *is* gating the relation: the recovery is
    # above the minimum by construction of what is printed.
    lowest=$(printf '%s\n' "$smoke" | sed -n 's/.*layer descended  *[0-9.]* m to \([0-9.]*\) m.*/\1/p')
    back=$(printf '%s\n' "$smoke" | sed -n 's/.*steady state back to \([0-9.]*\) m.*/\1/p')
    emit=$(printf '%s\n' "$smoke" | sed -n 's/.*(emission \([0-9.e+-]*\)).*/\1/p')
    check "the lowest the interface reaches (m)" 2.96 0.005 "$lowest" \
          "reaches 2.96 m at about t = 300" "$RENDER_DOC"
    check "what it recovers to by t = 600 (m)"   3.13 0.005 "$back" \
          "*recovers* to 3.13 m" "$RENDER_DOC"
    # 9.6e-10 of full scale is the *quantitative* half of "no fire in it": the
    # temperature pair above says the layer is under the threshold, this says by
    # how far. A tolerance of 5e-12 is the digit the document publishes to.
    check "the peak emission, fraction of full scale" 9.6e-10 5e-12 "$emit" \
          "emits 9.6e-10 of full scale" "$RENDER_DOC"

    # **"τ passes 10 before t = 100 s" is a relation and nothing above implies it.**
    # Every gated cell could match while the first sample to go optically thick had
    # moved past 100 s, because no row of the table is at that instant. So the
    # printed table is scanned for the first row whose optical depth exceeds 10 and
    # its own timestamp is asserted -- which is also the only thing here that would
    # notice the frame schedule stretching underneath the figures.
    checks=$((checks + 1))
    thick=$(printf '%s\n' "$smoke" |
            awk 'NF == 10 && $1 ~ /^[0-9]+$/ && $7 + 0 > 10 { print $2; exit }')
    if [ -n "$thick" ] && awk -v t="$thick" 'BEGIN { exit !(t > 0 && t < 100) }'; then
      printf '  %s✓%s and the layer is optically thick by t = %s s, inside the minute claimed\n' \
             "$green" "$off" "$thick"
    else
      printf '  %s✗%s the layer no longer passes tau = 10 before t = 100 s (first at "%s")\n' \
             "$red" "$off" "$thick"
      printf '      %s%s says "passes 10 before t = 100 s"%s\n' "$dim" "$RENDER_DOC" "$off"
      fails=$((fails + 1))
    fi

    # --- the counterfactual, which is the control for the negative finding -----
    #
    # "This fire has no fire in it" is only a statement about *this fire* if the
    # renderer can draw a glow when there is one, and the document's evidence for
    # that was a sentence nothing re-ran. Two runs, ~0.6 s each: at 10 MW the layer
    # is past the one-byte threshold and still puts no red-dominant pixel on the
    # screen; at 12 MW it is drawn as a light source. The pair is the point exactly
    # as the capsize bracket above is -- either number alone would look plausible
    # while the boundary between smoke and fire had moved.
    #
    # The pixel counts are rasterised on the device, so they carry a tolerance where
    # the temperatures do not: 100 pixels of a 16 000-pixel layer is far tighter
    # than the difference between a glow and none, and does not pretend a different
    # card would round every edge identically.
    # The emissions are quoted to two figures in the document where the tool prints
    # three, so -- as the front page's Cb and GM are -- the expectation is the
    # *published* value and the tolerance is the rounding it was published at.
    peakof() { printf '%s\n' "$1" | sed -n 's/.*this one peaked at \([0-9]*\) K.*/\1/p'; }
    emitof() { printf '%s\n' "$1" | sed -n 's/.*(emission \([0-9.e+-]*\)).*/\1/p'; }
    litof()  { printf '%s\n' "$1" | sed -n 's/.*  *\([0-9][0-9]*\) red-dominant pixels.*/\1/p'; }
    glow10=$("$SMOKE" --out=/tmp --frames=8 --duration=600 --power=10e6 2>&1)
    glow12=$("$SMOKE" --out=/tmp --frames=8 --duration=600 --power=12e6 2>&1)
    # **Zero tolerance on both temperatures**, which took a perturbation to settle:
    # `876 ± 1` cannot fail on a one-kelvin move, and one kelvin is the whole digit
    # the tool prints. Both come back exact over six repeats and over every frame
    # count from 8 to 48, so the tolerance is what was measured -- the §8 argument
    # for gating a deterministic integer at 0.
    check "at 10 MW, the peak layer temperature (K)" 876 0 "$(peakof "$glow10")" \
          "**876 K**" "$RENDER_DOC"
    check "at 10 MW, the peak emission, fraction of full scale" 6.7e-3 5e-5 \
          "$(emitof "$glow10")" "6.7e-3 of full scale" "$RENDER_DOC"
    check "at 10 MW, red-dominant pixels in the last frame" 0 0 "$(litof "$glow10")" \
          "**0** pixels of the frame come back red-dominant" "$RENDER_DOC"
    check "at 12 MW, the peak layer temperature (K)" 990 0 "$(peakof "$glow12")" \
          "**990 K**" "$RENDER_DOC"
    check "at 12 MW, the peak emission, fraction of full scale" 1.1e-1 5e-3 \
          "$(emitof "$glow12")" "1.1e-1 of full scale" "$RENDER_DOC"
    check "at 12 MW, red-dominant pixels in the last frame" 16309 100 "$(litof "$glow12")" \
          "**16 309** pixels of it are" "$RENDER_DOC"
  fi
else
  echo "  - smoke_view not built, skipping the fire and smoke figures"
  skipped="$skipped smoke_view"
fi

# **A gate that ran nothing used to report success.** With every tool absent this
# printed `ok — 0 published figures still match the tool` and exited 0, and
# `verify.sh` greps for `^ok — ` and passed it. That is the same shape as the
# `exit 0` on an unbuilt `ram_view` recorded above -- a green result standing for
# work that did not happen -- except that this one survives *every* block being
# skipped rather than the ones after a single point, so it is strictly the larger
# hole. The count was the only thing that would have shown it, and a count nobody
# diffs is not a signal.
#
# Two things close it. **A gate that checked nothing has not passed**, so zero
# checks is a failure in its own right. And a skip is named on the success line
# rather than only where it happened, because the summary is the line that gets
# read and pasted -- the per-block notes scroll off the top of a 900 s run. GPU
# blocks must still *skip* rather than fail on a machine with no device, which is
# why a named skip is not itself red: what was missing is that it be visible.
# --- the job system's grain figures, which are the first gate on `01` ----------
#
# `01-architecture.md` was one of the documents the note at the top of this file
# recorded as gating zero figures, on a reason -- design documents whose digits are
# section numbers and budgets nothing computes -- that is true of most of it and
# false of the table at §104-109. CI runs `job_bench`
# (`.github/workflows/ci.yml:49`) under a comment saying its figures are what the
# docs quote, and asserts only the exit code: the `zone_gpu_probe` hole above, in
# a second place.
#
# **This is the only tool here whose figures are load-dependent**, and that splits
# the block rather than making it unwritable.
#
# *Machine-independent, always checked.* `parallelForAuto`'s chunk count is not a
# measurement. `jobs.hpp:125-133` says why: on a 20 M element loop the lower clamp
# binds, so 12 lanes x 64 gives 768 chunks and 27 lanes x 64 gives 1728 whatever
# the probe measured. Confirmed bit-identical over thirty-one runs spanning load
# 1.7 to 10. The *grain* carries a tolerance of 2 where the chunk count carries 0,
# and they are not the same claim: grain is `ceil(remaining / 768)` and
# `remaining` moves by the length of the geometric probe prefix, which steps from
# 1365 elements to 341 on a machine slower than 5.9 ns/element. The chunk count is
# invariant under that by construction, since `ceil(r / ceil(r/768)) == 768`.
#
# *Load-dependent, checked only on a box wide enough to mean it.* `verify.sh full`
# runs on a 2-core GitHub runner where `JobSystem(23)` is 23 threads on 2 cores
# and every millisecond column is meaningless. A gate that goes red there teaches
# people to ignore red, which is worse than the hole it closes -- so the sweep is
# skipped by name rather than fudged with a tolerance wide enough to survive it.
#
# ~3 s.
if [ -x "$JOBS" ]; then
  bench=$("$JOBS" 2>&1)
  # Keyed on `$2 == "workers:"` and not on line numbers: §2's per-worker tables
  # and §3b's rows are both five numeric columns, and only §2's sit under a
  # `   N workers:` heading. The `grain for a` line closes the last table, which
  # is what stops 3b being read as a continuation of it.
  sweep() {
    printf '%s\n' "$bench" |
      awk -v w="$1" '$2 == "workers:" { inblock = ($1 == w); next }
                     /grain for a/ { inblock = 0 }
                     inblock && NF == 5 && $1 ~ /^[0-9]+$/ { print $1, $3, $4 }'
  }
  auto() {
    printf '%s\n' "$bench" |
      awk -v w="$1" -v c="$2" '/^3b\./ { in3b = 1; next } /^4\./ { in3b = 0 }
                               in3b && NF == 5 && $1 == w { print $c }'
  }

  check "auto-grain chunk count, 8 workers (12 lanes x 64)"  768 0 "$(auto 8 3)" \
        '12 lanes x 64 gives 768' engine/core/jobs.hpp
  check "auto-grain chunk count, 23 workers (27 lanes x 64)" 1728 0 "$(auto 23 3)" \
        '27 lanes x 64 gives 1728' engine/core/jobs.hpp
  check "the grain that falls out of it, 8 workers"  26040 2 "$(auto 8 2)" \
        'the grain comes out 26 040 and' engine/core/jobs.hpp
  check "the grain that falls out of it, 23 workers" 11574 2 "$(auto 23 2)" \
        'the grain comes out 26 040 and' engine/core/jobs.hpp

  if [ "$(nproc 2>/dev/null || echo 0)" -ge 24 ]; then
    uncont=$(printf '%s\n' "$bench" | awk '/^1\. Per-job/{i=1;next} /^2\./{i=0}
                                           i && NF==3 && $1=="0" {print $2}')
    # **The plateau is asserted as a *position*, not as a millisecond.** The
    # milliseconds move 1.2-1.5x between runs; which grain first comes within 15%
    # of the best does not. That is also the figure `01` actually publishes: the
    # chunk *size* where efficiency plateaus, whose microsecond column is a pure
    # function of the baseline and the chunk count.
    pgrain=$(sweep 8 | awk '{ g[NR]=$1; t[NR]=$3; if (b=="" || $3<b) b=$3 }
                            END { for (i=1;i<=NR;i++) if (t[i] <= 1.15*b) { print g[i]; exit } }')
    # A tolerance of 8 on 17 is not a rounding allowance, for the reason the
    # `water_probe` block above gives about 27.9 +/- 4.0: what has to stay true is
    # that dispatching a job costs tens of nanoseconds and not hundreds, because
    # that is the whole of the cancelled-Chase-Lev argument at §138-143.
    # Twenty-three runs here gave 14.8-16.1; a further fifteen elsewhere gave up
    # to 19.2.
    check "uncontended dispatch (ns/job)" 17 8 "$uncont" \
          '| Dispatch, uncontended (0 workers) | **15–19 ns/job** |' "$ARCH_DOC"
    check "grain at which 8 workers reach the plateau" 1024 0 "$pgrain" \
          '| Chunk size at which efficiency plateaus | **≈ 2 µs** |' "$ARCH_DOC"

    # **The grain penalty is gated as a bound and can only be gated as a bound.**
    # The grain-16 row moved 1.5x and 2.5x across runs while the plateau it is
    # divided by held to 1.05x, so the ratio spans 17.8-28.1x and 28.0-73.7x here.
    # The floor gated is *below* the published 15x, deliberately: the published
    # bound is the claim, and a gate set at the claim with 5% of headroom is a
    # tripwire that goes red on a correct build. See §112-129 for why the bound
    # moves down as the box gets quieter -- CI is not the quiet end of that curve,
    # but this check runs on developer boxes that might be.
    checks=$((checks + 1))
    pen8=$(sweep 8   | awk '$1 == 16 { w = $3 } { if (b=="" || $3<b) b=$3 } END { printf "%.1f", w/b }')
    pen23=$(sweep 23 | awk '$1 == 16 { w = $3 } { if (b=="" || $3<b) b=$3 } END { printf "%.1f", w/b }')
    if awk -v a="$pen8" -v b="$pen23" 'BEGIN { exit !(a >= 13 && b >= 20) }'; then
      printf '  %s✓%s the grain-16 penalty is at least 13x at 8 workers (%sx) and 20x at 23 (%sx)\n' \
             "$green" "$off" "$pen8" "$pen23"
    else
      printf '  %s✗%s the grain penalty fell through its floor: %sx at 8 workers, %sx at 23\n' \
             "$red" "$off" "$pen8" "$pen23"
      printf '      %s%s publishes "**at least 15×**" — below 13x the claim is gone%s\n' \
             "$dim" "$ARCH_DOC" "$off"
      fails=$((fails + 1))
    fi

    # **The tool now computes that span too, so the two are compared rather than
    # both printed.** `job_bench`'s verdict used to end on a hard-coded "the sweep
    # above spans roughly 30x" -- a literal in a `printf`, in the one program in a
    # position to measure it, at a third value the document does not carry and had
    # explicitly retired. It now tracks best and worst speedup as it sweeps and
    # prints them. That creates a second route to a number this block already
    # derives, which is the thing this file distrusts most, so it is closed the way
    # the repo's own rule says: assert they agree.
    #
    # **The relation is one-sided, and the tolerance is a rounding budget rather
    # than a fudge.** In exact arithmetic the tool's span cannot be smaller: the
    # block above takes grain 16 as the worst grain, the tool takes the worst of
    # all eight, so `vspan >= pen` always. But the two do not read the same
    # numbers. The tool tracks full-precision speedups; this block divides the
    # **printed** ms column, which carries two decimals. The denominator is the
    # small one -- 3.64 ms against 144.14 at 23 workers -- so its rounding is worth
    # 0.005/3.64 = 0.14% of the ratio, and `pen` is then itself rounded to one
    # decimal. Both push `pen` above the truth, and the first version of this check
    # allowed a flat 0.05 and went red on 36.84 against 36.9, which is the whole
    # error and none of the physics.
    #
    # So the budget is derived from the printed precision at the values actually
    # seen, not chosen: half a unit in the last place of each operand, propagated
    # through the quotient, plus half a unit of `pen`'s own rounding.
    checks=$((checks + 1))
    vspan8=$(printf '%s\n' "$bench"  | awk '/^ *grain span:/ { print $3 + 0; exit }')
    vspan23=$(printf '%s\n' "$bench" | awk '/^ *grain span:/ { print $7 + 0; exit }')
    slack() {  # $1 = ms at grain 16, $2 = smallest ms in the sweep, $3 = the span
      # Three terms, one per rounding in the chain and nothing else: the two ms
      # operands at half a unit in their second decimal, propagated through the
      # quotient; `pen`'s own `%.1f`; and the tool's `%.2f` on the span it prints.
      awk -v w="$1" -v b="$2" -v s="$3" \
          'BEGIN { printf "%.4f\n", s * (0.005 / w + 0.005 / b) + 0.05 + 0.005 }'
    }
    worst8=$(sweep 8   | awk '$1 == 16 { print $3; exit }')
    best8=$(sweep 8    | awk '{ if (b == "" || $3 < b) b = $3 } END { print b }')
    worst23=$(sweep 23 | awk '$1 == 16 { print $3; exit }')
    best23=$(sweep 23  | awk '{ if (b == "" || $3 < b) b = $3 } END { print b }')
    slack8=$(slack "$worst8" "$best8" "$pen8")
    slack23=$(slack "$worst23" "$best23" "$pen23")
    if [ -z "$vspan8" ] || [ -z "$vspan23" ]; then
      printf '  %s✗%s job_bench printed no parseable "grain span:" line\n' "$red" "$off"
      fails=$((fails + 1))
    elif awk -v a="$vspan8" -v b="$pen8" -v s="$slack8" \
             -v c="$vspan23" -v d="$pen23" -v t="$slack23" \
           'BEGIN { exit !(a >= b - s && c >= d - t) }'; then
      printf '  %s✓%s job_bench'"'"'s own span agrees with the swept table (%sx/%sx against %sx/%sx)\n' \
             "$green" "$off" "$vspan8" "$vspan23" "$pen8" "$pen23"
    else
      printf '  %s✗%s job_bench'"'"'s verdict span is below the table it printed: %sx/%sx against %sx/%sx\n' \
             "$red" "$off" "$vspan8" "$vspan23" "$pen8" "$pen23"
      printf '      %sthe worst of eight grains cannot be smaller than grain 16 alone;\n' "$dim"
      printf '      allowed %s and %s for the printed table'"'"'s two decimals%s\n' \
             "$slack8" "$slack23" "$off"
      fails=$((fails + 1))
    fi

    # **And that the tuner lands in the plateau**, which is the claim §123-129
    # makes and is a relation rather than a millisecond.
    #
    # **1.40, and 1.15 was inside the noise.** This is not a ratio between two
    # timings of the same thing, which is the form this repo trusts. It is a
    # *single* timed run of the tuner's grain over the **minimum of eight** swept
    # timings, and a minimum over eight noisy samples is biased low -- so the ratio
    # sits above 1 even when the tuner is perfect, and the bias grows with the
    # noise. Measured, fifteen runs on an idle box:
    #
    #     ratio8   0.969 .. 1.125   (four of fifteen above 1.07)
    #     ratio23  0.940 .. 1.028
    #
    # against a bound of 1.15. It went red at **1.282** inside a gate run, which is
    # the only sample taken under the load a gate actually runs at. The check three
    # lines above this one warns in those words that "a gate set at the claim with
    # 5% of headroom is a tripwire that goes red on a correct build"; this one was
    # set at 2.6% of headroom from a quiet-box sample.
    #
    # **There is real separation to aim at.** The tuner is deterministic -- grain
    # 26040 and 11574, 768 and 1728 chunks, on every one of twenty runs, and those
    # counts are gated exactly a few lines up, which is the structural half of this
    # claim. For the timing half, the nearest genuinely wrong landing is the
    # neighbouring swept row: 256 chunks costs 7.96 ms against the plateau's 5.31,
    # a ratio of **1.50**. Everything further off is 5.6x or worse. So 1.40 sits
    # above the loaded noise and below the cheapest real failure, and says so.
    checks=$((checks + 1))
    ratio8=$(sweep 8 | awk -v a="$(auto 8 4)" '{ if (b=="" || $3<b) b=$3 } END { printf "%.3f", a/b }')
    ratio23=$(sweep 23 | awk -v a="$(auto 23 4)" '{ if (b=="" || $3<b) b=$3 } END { printf "%.3f", a/b }')
    if awk -v a="$ratio8" -v b="$ratio23" 'BEGIN { exit !(a <= 1.40 && b <= 1.40) }'; then
      printf '  %s✓%s the auto-tuner lands in the plateau (%s and %s of the best swept grain)\n' \
             "$green" "$off" "$ratio8" "$ratio23"
    else
      printf '  %s✗%s the auto-tuner no longer lands in the plateau: %s and %s of the best grain\n' \
             "$red" "$off" "$ratio8" "$ratio23"
      fails=$((fails + 1))
    fi
  else
    echo "  - fewer than 24 hardware threads, skipping job_bench's load-dependent figures"
    skipped="$skipped job_bench(sweep)"
  fi
else
  echo "  - job_bench not built, skipping the job system's grain figures"
  skipped="$skipped job_bench"
fi

# --- the FLIP solver's studies, which nothing in the repository re-ran ------------
#
# `tools/flip_probe` is the purest case of the shape this file exists for:
# `06-roadmap.md`'s Phase 5 item, `engine/sim/flip.hpp` §1-§3 and
# `flip-escalation-design.md` all quote it, and nothing has ever run the probe.
#
# **This comment used to say `verify.sh` runs `flip_tests`, and it does not.**
# `grep flip scripts/verify.sh` returns nothing at all; `flip_tests`
# (`CMakeLists.txt:131`) is built by every gate and executed by no script, no CI
# step and nothing else in the tree. The coverage is not missing -- `test_flip.cpp`
# is compiled into `shipsim_tests` as well (`CMakeLists.txt:176`), which is what
# actually runs -- so the target is a developer convenience for running the FLIP
# tests alone, and is kept as one. What was wrong was the sentence, which named a
# thing that runs to explain why something else does not, and nothing tests a
# comment.
#
# **Unlike `job_bench` above, every figure here is gateable.** Two 70 s runs diff
# to nothing but wall-clock columns; filtered of those, three runs are
# byte-identical. So these carry the digit the document rounded to and nothing
# more, and several carry a tolerance of zero.
#
# **What runs is four of the five studies.** `--slosh` is the other 39 s -- 56% of
# the tool's runtime for two figures -- and is left out. Its `h=0.04`, two-cell row
# is the same configuration the transfer study's tank already reports as APIC
# `+4.45%`, so it is covered here for free, and its convergence claim is asserted
# inside the tool by `require("refinement improves the period at fixed amplitude
# in cells", monotone)` rather than published as a figure. `--quick` would be 3.4 s
# but drops the tank table, which is five of the eight transfer figures and the
# whole "the ordering reverses" argument. ~31 s.
if [ -x "$FLIP" ]; then
  flip=$("$FLIP" --sparse --solver --transfer --dam 2>&1)
  # Four columns is the empty-room table, six is the same water in three rooms;
  # both keyed on the room, neither on a line number.
  empty()   { printf '%s\n' "$flip" | awk -v c="$1" 'NF == 4 && $1 == "400^3" { print $c; exit }'; }
  wet()     { printf '%s\n' "$flip" | awk -v c="$1" 'NF == 6 && $1 == "400^3" { print $c }'; }
  poisson() { printf '%s\n' "$flip" | awk -v k="$1" -v c="$2" '$1 == k && $2 == "1e-13" && NF == 4 { print $c }'; }
  rot()     { printf '%s\n' "$flip" | awk -v m="$1" -v c="$2" '$1 == m && NF == 4 { print $c }'; }
  tank()    { printf '%s\n' "$flip" | awk -v m="$1" -v c="$2" '$1 == m && NF == 6 { gsub(/[+%]/, "", $c); print $c }'; }

  check "an empty 400^3 room allocates tiles" 0 0 "$(empty 3)" \
        '**zero tiles and zero bytes**' "$DOC"
  check "an empty 400^3 room allocates bytes" 0 0 "$(empty 4)" \
        '**zero tiles and zero bytes**' "$DOC"
  check "the same water in a 400^3 room, tiles" 48 0 "$(wet 2)" \
        'allocates the same 48 tiles and produces **bit-identical** particle positions' "$DOC"
  check "bytes per allocated cell" 116.1 0.05 \
        "$(printf '%s\n' "$flip" | sed -n 's/.*costs \([0-9.]*\) bytes per allocated cell.*/\1/p')" \
        '7.4 GB dense at the 116 bytes a cell this structure costs' "$DOC"

  # **"the same 48 tiles" is a claim about three rooms, and gating one of them
  # does not make it.** A structure that had started counting the domain would
  # give 48 in the 20^3 room and something larger in the 400^3 -- and the row
  # above, taken from the 400^3, would be the one that moved, while the claim
  # under test is that the three are *equal*. Same argument the `section_probe`
  # block makes for checking the negative control beside the tied one.
  checks=$((checks + 1))
  distinct=$(printf '%s\n' "$flip" | awk 'NF == 6 && $1 ~ /\^3$/ { print $2, $3, $4 }' | sort -u | wc -l)
  if [ "$distinct" = 1 ]; then
    printf '  %s✓%s the 20^3, 100^3 and 400^3 rooms allocate one identical row of tiles/cells/bytes\n' \
           "$green" "$off"
  else
    printf '  %s✗%s the room extent has entered the arithmetic — %s distinct rows across three rooms\n' \
           "$red" "$off" "$distinct"
    fails=$((fails + 1))
  fi

  check "SOR sweeps to a 1e-13 residual on 24 cells" 1170 0 "$(poisson 24 3)" \
        'relaxation of 1.7 needs **1 170 sweeps** to reach a 1e-13 residual on 24' "$FLIP_DOC"
  check "Jacobi-CG iterations for the same"            24 0 "$(poisson 24 4)" \
        '- CG solve dominates (24 iterations for machine precision)' "$ESC_DOC"
  check "SOR sweeps to a 1e-13 residual on 64 cells"  8221 0 "$(poisson 64 3)" \
        'At 64 cells it is 8 221 sweeps' "$FLIP_DOC"
  check "Jacobi-CG iterations for the same"             64 0 "$(poisson 64 4)" \
        'At 64 cells it is 8 221 sweeps' "$FLIP_DOC"

  check "PIC angular momentum kept over ten transfers"  0.892 0.0005 "$(rot PIC 3)" \
        'keeps **89.2%** of its angular momentum through ten particle→grid→particle' "$DOC"
  check "APIC angular momentum kept over ten transfers" 0.982 0.0005 "$(rot APIC 3)" \
        'transfers and APIC **98.2%** — 1.08% against 0.18% per transfer, **6.1× less**.' "$DOC"
  # FLIP is the identity here by construction, and it is the control that says the
  # other two rows are measuring a transfer at all rather than a decaying number.
  check "FLIP is the identity, which is the control"      1.0 0.0 "$(rot FLIP 3)" \
        'transfers and APIC **98.2%** — 1.08% against 0.18% per transfer, **6.1× less**.' "$DOC"
  check "how much less APIC loses per transfer"           6.1 0.05 \
        "$(printf '%s\n' "$flip" | sed -n 's/.*-- \([0-9.]*\) times less\..*/\1/p')" \
        'transfers and APIC **98.2%** — 1.08% against 0.18% per transfer, **6.1× less**.' "$DOC"

  # The tank pair is the point exactly as the capsize bracket above is: the whole
  # argument for APIC is that the ordering *reverses* between the rotation study
  # and this one, and either figure alone is consistent with a solver that had
  # lost one of them. The crossings go with the periods -- a period read off a
  # signal that crossed the centreline three times is not the same measurement as
  # one read off eight.
  check "FLIP particle noise in the tank"  0.876 0.0005 "$(tank FLIP 4)" \
        'In a sloshing tank the ordering reverses: FLIP carries **87.6%** particle-borne' "$DOC"
  check "APIC particle noise in the tank"  0.013 0.0006 "$(tank APIC 4)" \
        'velocity the grid never sees, APIC **1.3%**' "$DOC"
  check "APIC first-mode period error (%)"  4.45 0.005 "$(tank APIC 3)" \
        'comes out **+4.45%**, PIC **+9.16%** and FLIP **+9.81%** -- and PIC crosses' "$FLIP_DOC"
  check "PIC first-mode period error (%)"   9.16 0.005 "$(tank PIC 3)" \
        'comes out **+4.45%**, PIC **+9.16%** and FLIP **+9.81%** -- and PIC crosses' "$FLIP_DOC"
  check "FLIP first-mode period error (%)"  9.81 0.005 "$(tank FLIP 3)" \
        'comes out **+4.45%**, PIC **+9.16%** and FLIP **+9.81%** -- and PIC crosses' "$FLIP_DOC"
  check "PIC centreline crossings in five seconds"  3 0 "$(tank PIC 6)" \
        'the centreline three times in five seconds where APIC crosses eight, because' "$FLIP_DOC"
  check "APIC centreline crossings in five seconds" 8 0 "$(tank APIC 6)" \
        'the centreline three times in five seconds where APIC crosses eight, because' "$FLIP_DOC"

  # **The dam break has no figure worth gating and one relation that is.** Its
  # front position, peak speed and clamp count are outputs of a scheme with no
  # closed form to check them against -- the tool says so itself. What it does
  # have is Ritter's dry-bed celerity as a ceiling the front must never pass,
  # which is the one thing here a wrong advection would break and a tolerance
  # would hide.
  checks=$((checks + 1))
  ritter=$(printf '%s\n' "$flip" | sed -n "s/.*2 sqrt(g h0) = \([0-9.]*\) m\/s.*/\1/p")
  peak=$(printf '%s\n' "$flip" | awk 'NF == 7 && $1 ~ /^0\./ { if ($3 > m) m = $3 } END { print m + 0 }')
  if awk -v p="$peak" -v r="$ritter" 'BEGIN { exit !(p > 0 && p < r) }'; then
    printf '  %s✓%s the dam-break front stays under Ritter (%s m/s against %s)\n' \
           "$green" "$off" "$peak" "$ritter"
  else
    printf '  %s✗%s the dam-break front reached %s m/s against Ritter %s\n' \
           "$red" "$off" "$peak" "$ritter"
    fails=$((fails + 1))
  fi
else
  echo "  - flip_probe not built, skipping the FLIP solver's studies"
  skipped="$skipped flip_probe"
fi

# --- the Tier-2 zone at ship scale, which verify.sh runs for its exit code -------
#
# `zone_probe` is the third and last tool in the shape this file exists for.
# `verify.sh:671` runs it -- under a comment saying "it is why the zone's cost has
# to stay visible" -- and asserts `^ok$`. Nothing checked a figure. It publishes
# into `02-simulation.md`'s cost table, `06-roadmap.md`'s Phase 4 item,
# `07-fem-spike-findings.md` §6, and the headers of `zone.hpp` and
# `solid_shell.hpp` -- **the two engine headers this gate reads nothing from**,
# and between them the densest concentration of ungated tool output in the repo.
# One of its figures was already wrong when this block was written: the cost table
# predicted 900 core-seconds against its own measurement of 385, on the row below.
#
# **It is bit-deterministic in everything but wall clock.** Two runs diff to the
# `ms to know`, `decided in us`, `s wall` and profile columns and nothing else, so
# every physics figure here carries tolerance 0. The seconds are not gated at all:
# they are this box on 23 threads, and the same solve is 4.8 s with
# `--forms-cache=never`, which is the configuration the table used to publish
# without saying so.
#
# ~8 s for the whole thing, which is the cheapest block in this file per figure.
if [ -x "$ZONE" ]; then
  zone=$("$ZONE" --radius=3.0 --depth=0.22 2>&1)
  zfield() { printf '%s\n' "$zone" | sed -n "s/$1/\1/p" | head -1; }

  check "the zone the criterion promotes, elements" 224 0 \
        "$(zfield '^zone *: [0-9]* panels -> \([0-9]*\) elements.*')" \
        '| Zone | 14 panels → **224 elements**, 522 nodes, 24.0 m² of 12 mm plating |' "$ZONE_DOC"
  check "and nodes" 522 0 "$(zfield '^zone *: .* \([0-9]*\) nodes.*')" \
        '| Zone | 14 panels → **224 elements**, 522 nodes, 24.0 m² of 12 mm plating |' "$ZONE_DOC"
  check "and plating area (m2)" 24.0 0.05 "$(zfield '^zone *: .*, \([0-9.]*\) m2 of 12 mm.*')" \
        '| Zone | 14 panels → **224 elements**, 522 nodes, 24.0 m² of 12 mm plating |' "$ZONE_DOC"
  check "the critical timestep (us)" 1.816 0.0005 \
        "$(zfield '^cost *: dt = \([0-9.]*\) us.*')" \
        '| 2 (0.25 × 0.175 m) | 1.816 µs | 1.816 µs | 1.000 |' "$ZONE_DOC"
  check "steps the solve actually took" 21290 0 \
        "$(zfield '^solve *: \([0-9]*\) steps.*')" \
        '| Delivered, 23 workers | 21 290 steps, **2.6 s of wall time**, 0.55 µs/element/step |' "$ZONE_DOC"
  check "the work the punch cost (MJ)" 2.755 0.0005 \
        "$(zfield '^energy *: in \([0-9.]*\) MJ.*')" \
        'the 2.755 MJ a 0.22 m punch costs' "$DOC"
  check "elements torn" 80 0 "$(zfield '^damage *: \([0-9]*\) of 224 elements deleted.*')" \
        '| Delivered, 23 workers | 21 290 steps, **2.6 s of wall time**, 0.55 µs/element/step |' "$ZONE_DOC"
  check "and what the energy-driven drive tears instead" 42 0 \
        "$(zfield '^ *\([0-9]*\) element(s) torn against 80.*')" \
        '45 970 steps against 21 290 — because a striker that has nearly stopped crawls,' "$ZONE_DOC"
  check "the steps that drive takes" 45970 0 \
        "$(zfield '.*against 6; \([0-9]*\) steps against 21290.*')" \
        '45 970 steps against 21 290 — because a striker that has nearly stopped crawls,' "$ZONE_DOC"
  check "the pre-load the girder puts through the patch (MPa)" 13.1 0.05 \
        "$(zfield '^preload: .* -> \([0-9.]*\) MPa through the patch.*')" \
        "| zone handed the girder's 13.1 MPa | **20.25 MN**, +7.1% |" "$ZONE_DOC"
  check "and the force at 0.078 m under it (MN)" 20.25 0.005 \
        "$(printf '%s\n' "$zone" | awk '$1=="solid-shell" && $2=="FEM" && NF==5 { print $4 }')" \
        "| zone handed the girder's 13.1 MPa | **20.25 MN**, +7.1% |" "$ZONE_DOC"

  # **The control, which was ungated while its result was gated.** §2's table is two
  # rows and the finding is the difference between them: 18.90 MN with the patch
  # told it starts unstressed, 20.25 MN carrying the girder's own 13.1 MPa, +7.1%.
  # The loaded row and the 13.1 MPa were both gated; the unstressed row was not, so
  # the *published* half of a two-row comparison stood on a single measurement. 8 s.
  #
  # The `preload:` line is read too, because until this commit it announced
  # `applied 1` on this very run -- `preloadFor`'s verdict that a pre-stress exists,
  # printed where a reader would take it for what the solve received.
  bare=$("$ZONE" --radius=3.0 --depth=0.22 --no-preload 2>&1)
  bareforce=$(printf '%s\n' "$bare" |
              awk '$1=="solid-shell" && $2=="FEM" && NF==5 { print $4 }')
  check "the force at 0.078 m with no pre-stress (MN)" 18.90 0.005 "$bareforce" \
        '| zone told it starts unstressed | **18.90 MN** — the figure this file published |' \
        "$ZONE_DOC"
  check "and the solve is told it got none" 0 0 \
        "$(printf '%s\n' "$bare" | sed -n 's/^preload: .* given to the solve \([0-9]*\)$/\1/p')" \
        '| zone told it starts unstressed | **18.90 MN** — the figure this file published |' \
        "$ZONE_DOC"
  # The +7.1% the table publishes, derived here rather than transcribed: it is the
  # whole content of the comparison and neither row alone carries it.
  # 0.06, because both rows are published to two decimals: 18.90 and 20.25 each
  # carry +/-0.005, which is +/-0.04 on the ratio. Measured 7.14 against a doc that
  # rounds to 7.1, so the bound is the rounding and nothing more.
  check "the pre-stress is worth +7.1% at 0.078 m" 7.1 0.06 \
        "$(awk -v a="$bareforce" -v b="20.25" 'BEGIN { printf "%.2f", 100.0 * (b / a - 1.0) }')" \
        "| zone handed the girder's 13.1 MPa | **20.25 MN**, +7.1% |" "$ZONE_DOC"

  # **The two cost estimators, checked against each other and not only against the
  # document.** `zone::estimatedCost` works off `kPlasticMicroseconds` and the
  # patch's own critical timestep; `promotion.hpp`'s `coreSecondsPerElement` is a
  # per-element constant. The tool prints both on the same run -- 382 on the `cost`
  # line and 381 on the `promote` line -- and they were **2.35x apart** until
  # someone read them side by side. That is this repo's most-repeated defect: two
  # computations of one quantity that nothing compares. Gating either alone would
  # not have caught it, so the ratio is gated too.
  predicted=$(zfield '^ *predicted \([0-9]*\) core-seconds.*')
  promoted=$(printf '%s\n' "$zone" | sed -n 's/^promote: .*224 element(s) active, \([0-9]*\) core-s\/s.*/\1/p' | head -1)
  check "the cost the zone predicts for itself" 382 2 "$predicted" \
        '| Predicted | **382 core-seconds** per simulated second — `zone::estimatedCost`, which `zone_probe` prints |' "$ZONE_DOC"
  check "and the cost the promoter charges for the same zone" 381 2 "$promoted" \
        'them side by side. 0.3% apart is a corroboration; 2.35x apart was two numbers for' "$ZONE_DOC"
  checks=$((checks + 1))
  if awk -v a="$predicted" -v b="$promoted" 'BEGIN { exit !(a > 0 && b > 0 && a/b < 1.05 && b/a < 1.05) }'; then
    printf '  %s✓%s the promoter and the zone agree on what the zone costs (%s against %s)\n' \
           "$green" "$off" "$promoted" "$predicted"
  else
    printf '  %s✗%s the two cost estimators have parted again: promoter %s, zone %s\n' \
           "$red" "$off" "$promoted" "$predicted"
    printf '      %sthey were 2.35x apart once and both were printed on this same run%s\n' "$dim" "$off"
    fails=$((fails + 1))
  fi

  # **The predicted elastic/plastic ratio against the measured one**, which is the
  # check that would have caught `kElasticMicroseconds` at 0.273 -- a value that
  # made the elastic path 11.4x cheaper than the plastic one where it is 2.4x, and
  # that `estimatedCost` used to refuse zones at promotion. The unit test beside it
  # asserted `3.1 / 0.273` against `estimatedCost(p,true)/estimatedCost(p,false)`,
  # which is those same two constants with everything else cancelled: it could not
  # fail, and it held for years at a figure wrong by nine.
  #
  # **Read at `--depth=0.05`, and the depth is the point.** At the standard 0.22 the
  # plastic run tears 80 of its 224 elements and stops integrating them while the
  # elastic run never tears and carries all 224, so the ratio there -- 1.98 -- is
  # between two different amounts of work. At 0.05 and 0.08 nothing tears in either
  # and the ratio is 2.414 and 2.395, agreeing to 0.4%. The `deleted 0` check below
  # is that precondition asserted rather than assumed.
  #
  # `--threads=1` for the same reason a ratio is used at all: at 23 workers the
  # elastic solve is barrier-bound -- less work per element to amortise the sync --
  # and comes out *slower* per element-step than the plastic one, 0.76 us against
  # 0.55. Wall time over element-steps is not a per-element cost unless the thread
  # count is held down. ~10 s.
  zplastic=$("$ZONE" --radius=3.0 --depth=0.05 --threads=1 2>&1)
  zelastic=$("$ZONE" --radius=3.0 --depth=0.05 --threads=1 --elastic 2>&1)
  wallOf() { printf '%s\n' "$1" | sed -n 's/^solve *: [0-9]* steps in \([0-9.]*\) s wall.*/\1/p'; }
  costOf() { printf '%s\n' "$1" | sed -n 's/^ *predicted \([0-9]*\) core-seconds.*/\1/p'; }
  tornOf() { printf '%s\n' "$1" | sed -n 's/^damage *: \([0-9]*\) of [0-9]* elements deleted.*/\1/p'; }

  check "the shallow punch tears nothing, plastic"  0 0 "$(tornOf "$zplastic")" \
        'depth 0.05   2.68 s / 1.11 s = 2.414     0 of 224 elements deleted' engine/sim/zone.cpp
  check "nor elastic, so the two runs integrate the same elements" 0 0 "$(tornOf "$zelastic")" \
        'depth 0.05   2.68 s / 1.11 s = 2.414     0 of 224 elements deleted' engine/sim/zone.cpp

  checks=$((checks + 1))
  measuredRatio=$(awk -v p="$(wallOf "$zplastic")" -v e="$(wallOf "$zelastic")" \
                      'BEGIN { if (e > 0) printf "%.3f", p / e; else print "0" }')
  predictedRatio=$(awk -v p="$(costOf "$zplastic")" -v e="$(costOf "$zelastic")" \
                       'BEGIN { if (e > 0) printf "%.3f", p / e; else print "0" }')
  if awk -v m="$measuredRatio" -v p="$predictedRatio" \
        'BEGIN { exit !(m > 0 && p > 0 && m / p < 1.15 && p / m < 1.15) }'; then
    printf '  %s✓%s the elastic/plastic cost the constants predict matches the solver (%s against %s)\n' \
           "$green" "$off" "$predictedRatio" "$measuredRatio"
  else
    printf '  %s✗%s the cost constants have drifted from the solver: predicted %sx, measured %sx\n' \
           "$red" "$off" "$predictedRatio" "$measuredRatio"
    printf '      %skPlasticMicroseconds / kElasticMicroseconds in engine/sim/zone.cpp%s\n' "$dim" "$off"
    fails=$((fails + 1))
  fi
else
  echo "  - zone_probe not built, skipping the Tier-2 zone's figures"
  skipped="$skipped zone_probe"
fi

# Did anything we read get rebuilt while we were reading it? Reported before the
# figure verdict, because if it did then the verdict is about two programs and
# every other line below is unsafe to act on.
stampsAfter=$(toolStamps)
if [ "$stampsBefore" != "$stampsAfter" ]; then
  printf '%s✗%s a tool changed on disk while the gate was reading it — this run\n' "$red" "$off"
  printf '      compared figures against more than one build, and its verdict\n'
  printf '      means nothing. Re-run on a tree nobody is building.\n'
  printf '      %sbefore:%s%s\n' "$dim" "$stampsBefore" "$off"
  printf '      %safter: %s%s\n' "$dim" "$stampsAfter" "$off"
  exit 1
fi

if [ "$checks" -eq 0 ]; then
  printf '%s✗%s the gate checked nothing at all — no tool it reads was built\n' "$red" "$off"
  echo "0 of 0 published figures could be checked;$skipped did not run"
  exit 1
fi
if [ "$fails" -eq 0 ]; then
  if [ -n "$skipped" ]; then
    echo "ok — $checks published figures still match the tool, and these did not run:$skipped"
  else
    echo "ok — $checks published figures still match the tool"
  fi
  exit 0
fi
echo "$fails of $checks published figures have drifted"
[ -n "$skipped" ] && echo "   and these blocks did not run:$skipped"
exit 1
