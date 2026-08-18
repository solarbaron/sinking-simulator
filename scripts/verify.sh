#!/usr/bin/env bash
#
# The verification gate.
#
# Everything this project claims about itself is checked here. Each level is a
# superset of the one before:
#
#   quick     configure, build warning-clean, run the test suite      (~1 min)
#   full      + fresh rebuild, GPU tools, flooding scenarios, repeats (~40 min)
#   sanitize  + AddressSanitizer and ThreadSanitizer builds           (~70 min)
#   all       everything above
#   selftest  only the gate's own negative controls                   (~2 s)
#
# Measured, on a box with five other builds running: `full` took 2317 s, of which
# the figures check is about a third and the six repeat runs another third. The
# older ~400 s figure predates a suite that now runs 198 869 checks. Re-measure
# rather than trusting either number.
#
# Warnings are failures, not notes. The build is -Wall -Wextra -Wpedantic and has
# been warning-clean since Phase 0; the moment that stops being true it stops
# being a useful signal, so the gate treats a single warning as a red build.
#
# Incremental builds only recompile what changed, so they cannot see a warning in
# a file nobody touched. `full` and above therefore configure a throwaway build
# directory and compile everything from scratch.
#
# **Every check in here has been made to go red on purpose.** Three holes of the
# same shape were found in one week -- a check that could not report the failure
# it existed to catch -- and reading the code is what missed them all three
# times. `selftest` is the standing answer: it drives each reporter with a
# failure it was *not* written for and asserts on what came back. It runs at the
# head of every level, because a gate whose reporters are broken is worse than no
# gate: it is a green light with nothing behind it.
set -uo pipefail
# Note on `-e`, which is deliberately *not* set: the gate must keep going after a
# failure so one red step does not hide the other twelve, and it counts failures
# rather than propagating exit status. The cost is that every failure has to be
# noticed explicitly, so: no `fail` may be called inside a `$( )` or a pipeline
# (the FAILURES increment would happen in a subshell and be discarded), and no
# step may end in `|| true` unless its status has already been examined. Both
# rules have been broken here before; `selftest` now checks the first and the
# scenario block below is where the second one cost a real diagnosis.

readonly ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT" || { echo "cannot cd to $ROOT" >&2; exit 2; }

LEVEL="${1:-quick}"
case "$LEVEL" in
  quick|full|sanitize|all|selftest) ;;
  -h|--help|help)
    # Print the header comment, however long it grows. A fixed line range meant
    # `--help` ended in `set -uo pipefail` the first time someone added a
    # paragraph.
    awk 'NR == 1 { next } /^#/ { sub(/^# ?/, ""); print; next } { exit }' "${BASH_SOURCE[0]}"
    exit 0 ;;
  *) echo "usage: $0 [quick|full|sanitize|all|selftest]" >&2; exit 2 ;;
esac

readonly dim=$'\033[2m' off=$'\033[0m'
readonly BOLD=$'\033[1m' RED=$'\033[0;31m' GREEN=$'\033[0;32m' YELLOW=$'\033[0;33m' OFF=$'\033[0m'
FAILURES=0
SKIPPED=0
SKIPPED_WHAT=""
STARTED=$SECONDS

section() { printf '\n%s== %s ==%s\n' "$BOLD" "$1" "$OFF"; }
pass()    { printf '  %s✓%s %s\n' "$GREEN" "$OFF" "$1"; }
fail()    { printf '  %s✗ %s%s\n' "$RED" "$1" "$OFF"; FAILURES=$((FAILURES + 1)); }
# A skip is counted and *named on the result line*, not only where it happened.
# The result line is what gets read and pasted; the per-step notes are three
# hundred lines up in a run that takes an hour, and "passed" over six silently
# skipped checks reads exactly like "passed".
skip()    { printf '  %s—%s %s\n' "$YELLOW" "$OFF" "$1"; SKIPPED=$((SKIPPED + 1))
            SKIPPED_WHAT="$SKIPPED_WHAT ${1%% *}"; }
# The verdict line every level ends on. Never just "passed".
verdict_line() {
  local what="$1"
  if [ "$FAILURES" -ne 0 ]; then
    fail "$FAILURES failure(s) in $((SECONDS - STARTED))s"
  elif [ "$SKIPPED" -ne 0 ]; then
    pass "$what passed in $((SECONDS - STARTED))s — but $SKIPPED check(s) did not run:$SKIPPED_WHAT"
  else
    pass "$what passed in $((SECONDS - STARTED))s"
  fi
}
note()    { printf '    %s%s%s\n' "$dim" "$1" "$off"; }

# Nothing in this gate runs unbounded. The characteristic failure in this repo is
# a **hang**, not an assertion -- five of 196 mutations in fire.cpp/thermal.cpp
# killed the suite with zero failing checks by turning a nine-second run into an
# hours-long one -- and a gate that hangs reports nothing at all, which is worse
# than a gate that reports the wrong thing. These are bounds on "this cannot
# possibly still be working", not expectations: the box is usually shared.
BUILD_TIMEOUT=3600      # a clean rebuild of the whole tree on a busy box
STEP_TIMEOUT=1800       # any one tool run or test-suite run
SCENARIO_TIMEOUT=900    # one 900 s flooding scenario, which runs far faster than real time

# ------------------------------------------------------------- reporting ------
# One place that turns an exit status into a sentence, because there were five
# copies of this and they did not agree. `timeout` returns 124 when it fires, and
# the gate used to print a bare `exited 124` and leave the reader to know that --
# a hang and an assertion failure want completely different next steps.
why() {
  local status="$1" limit="${2:-}"
  if [ "$status" -eq 124 ] && [ -n "$limit" ]; then
    printf 'HUNG: no result after %s s, killed' "$limit"
  elif [ "$status" -eq 127 ]; then
    printf 'exited 127 (command not found)'
  elif [ "$status" -gt 128 ] && [ "$status" -lt 160 ]; then
    printf 'died on signal %d' "$((status - 128))"
  else
    printf 'exited %d' "$status"
  fi
}

# Everything a failing step is allowed to tell us, in one function, because the
# three copies of this drifted apart and two of them summarised a failure by
# grepping the log for the text a *foreseen* failure would contain: a compiler
# killed by a full disk never says `error:`, so the gate printed a bare red cross
# and nothing else. Named lines first, then the tail **regardless of whether
# anything matched**, and the status either way.
#
# The named lines are not a nicety. An AddressSanitizer report puts its `ERROR:`
# thirty lines above the shadow-byte legend, so `tail` alone diagnoses the
# legend; that was measured with a negative control, not guessed.
#
# The exit status is not printed here: every caller puts `why` in its own `fail`
# line, so the headline already carries it. Keep it that way -- a status only in
# the detail lines is a status nobody reads.
diagnose() {
  local log="$1" status="$2" limit="${3:-}"
  if [ ! -r "$log" ]; then
    printf '      %s(no log: nowhere to write one, or it was removed under us)%s\n' "$dim" "$off"
    return
  fi
  if [ ! -s "$log" ]; then
    printf '      %s(it printed nothing at all)%s\n' "$dim" "$off"
    return
  fi
  # -a because a crashed process can leave NUL bytes in the log, and grep then
  # says "binary file matches" and prints not one line of the diagnosis.
  # The marker list, and then the thing the marker list cannot know about.
  #
  # **`check-figures.sh` reports a drift with a line beginning `  ✗` and the word
  # "publishes", and matched not one pattern here** -- so a drifted figure printed
  # its diagnosis into the log and the reader got `--- last lines ---` and eight
  # lines of the passing tail, twice in one session. That is `expect_ok`'s own
  # recorded defect ("a failure path that only reports the failures it already
  # anticipated is not a failure path") in the function written to fix it.
  #
  # `✗` is added to the markers, and the tail is no longer the only fallback: when
  # nothing matches at all, the *first* lines are shown too, because a report that
  # fails early and then prints two hundred successes hides its own verdict behind
  # `tail`.
  local matched
  matched=$(grep -aE '✗|error:|ERROR: .*Sanitizer|SUMMARY: .*Sanitizer|runtime error:|^  FAIL|Killed|No space left|Cannot allocate|cannot open|fatal error|ninja: build stopped|Assertion .* failed' \
                 "$log" | head -10)
  if [ -n "$matched" ]; then
    printf '%s\n' "$matched"
  else
    printf '      %s(nothing matched the failure markers; showing the head)%s\n' "$dim" "$off"
    head -8 "$log"
  fi
  printf '      %s--- last lines ---%s\n' "$dim" "$off"
  tail -8 "$log"
}

# Run one command with a bound, capturing everything it says. Sets STEP_STATUS
# and STEP_LOG; returns non-zero only when the step could not be *started*.
#
# `mktemp` failing is the full-disk case arriving one function early: `log` comes
# back empty, every redirection into it fails, and the step gets reported as
# whatever ran first rather than as "there is nowhere to write". Note the
# `|| { ...; }` rather than a helper: `fail` inside `$( )` increments the counter
# in a subshell, and the increment is then thrown away.
STEP_LOG=""
STEP_STATUS=0
run_step() {
  local label="$1" limit="$2"; shift 2
  STEP_LOG=""
  STEP_STATUS=0
  STEP_LOG="$(mktemp -t shipsim-gate.XXXXXXXX)" \
    || { fail "$label: cannot create a log file — is ${TMPDIR:-/tmp} full or unwritable?"; return 1; }
  timeout "$limit" "$@" >"$STEP_LOG" 2>&1
  STEP_STATUS=$?
  return 0
}

# The verdict on a completed step: non-zero status, or zero status without the
# line that means success. Both print the same diagnosis, because "it exited 0
# and never said it passed" is every bit as much an unexplained failure.
verdict() {
  local label="$1" expect="$2" limit="$3"
  if [ "$STEP_STATUS" -ne 0 ]; then
    fail "$label — $(why "$STEP_STATUS" "$limit")"
    diagnose "$STEP_LOG" "$STEP_STATUS" "$limit"
    return 1
  fi
  if [ -n "$expect" ] && ! grep -qaE "$expect" "$STEP_LOG"; then
    fail "$label exited 0 but never reported '$expect'"
    diagnose "$STEP_LOG" "$STEP_STATUS" "$limit"
    return 1
  fi
  # `${x:-ok}` rather than `... || echo ok`: the fallback in a pipeline fires on
  # the *pipeline's* status, which depends on `pipefail` still being set three
  # hundred lines away.
  local summary; summary="$(grep -oaE '[0-9]+ checks, [0-9]+ failures' "$STEP_LOG" | tail -1)"
  pass "$label — ${summary:-ok}"
  return 0
}

# Run a command, requiring success and (optionally) a pattern in its output.
expect_ok() {
  local label="$1" expect="$2"; shift 2
  run_step "$label" "$STEP_TIMEOUT" "$@" || return 1
  verdict "$label" "$expect" "$STEP_TIMEOUT"
  local rc=$?
  rm -f "$STEP_LOG"
  return $rc
}

# A tool that needs no GPU is built on every machine, so a missing one is a build
# failure and never a skip. `check-figures.sh` was gated on `ram_view` -- a
# Vulkan-only binary -- which meant every CPU-only figure in it went unchecked,
# with no line printed at all, on any machine without a device.
require_built() {
  [ -x "$1" ] && return 0
  fail "$1 is missing or not executable — it needs no GPU, so this is a broken build, not a machine without a device"
  return 1
}

# ------------------------------------------------------------------ build -----
# Build into `dir`, failing on any compiler warning. `--target t` builds one
# target; extra cmake args follow.
#
# There is exactly one of these now. The -O3 step below had its own copy, and
# that copy still had the hole this one was fixed for: it greped `error:` and
# printed nothing when the build died without saying it, and its configure step
# sent cmake's output to /dev/null so there was nothing to print at all.
#
# `--fresh` says the directory was just removed, so the build has to *compile*
# something. Without it the two steps whose whole reason for existing is that
# they recompile everything -- the clean rebuild and the -O3 pass -- report
# `built warning-clean` off a `ninja: no work to do.`, which is what a `rm -rf`
# that did not take leaves behind. A build that compiled nothing cannot have seen
# a warning either, and "an incremental build cannot see a warning it did not
# recompile" is the sentence that step is named after.
build_into() {
  local dir="$1"; shift
  local target="" fresh=0
  while :; do
    case "${1:-}" in
      --target) target="$2"; shift 2 ;;
      --fresh)  fresh=1; shift ;;
      *) break ;;
    esac
  done
  local what="$dir${target:+ ($target)}"
  local log
  log="$(mktemp -t shipsim-gate.XXXXXXXX)" \
    || { fail "cannot create a log file for $what — is ${TMPDIR:-/tmp} full or unwritable?"; return 1; }

  timeout "$BUILD_TIMEOUT" cmake -S . -B "$dir" -G Ninja "$@" >"$log" 2>&1
  local status=$?
  if [ "$status" -ne 0 ]; then
    fail "cmake configure failed for $what — $(why "$status" "$BUILD_TIMEOUT")"
    diagnose "$log" "$status" "$BUILD_TIMEOUT"
    rm -f "$log"; return 1
  fi

  local goal=()
  [ -n "$target" ] && goal=("$target")
  timeout "$BUILD_TIMEOUT" ninja -C "$dir" ${goal[@]+"${goal[@]}"} >"$log" 2>&1
  status=$?
  if [ "$status" -ne 0 ]; then
    fail "build failed in $what — $(why "$status" "$BUILD_TIMEOUT")"
    diagnose "$log" "$status" "$BUILD_TIMEOUT"
    rm -f "$log"; return 1
  fi

  # A log we cannot read is not a log with no warnings in it. Left to itself,
  # `warnings` would come back empty, `[ "" -ne 0 ]` would error, the `if` would
  # take the *else*, and an unreadable build log would report warning-clean.
  if [ ! -r "$log" ]; then
    fail "the build log for $what disappeared — cannot say whether it was warning-clean"
    return 1
  fi
  local warnings; warnings="$(grep -acE 'warning:' "$log")"
  case "$warnings" in ''|*[!0-9]*) warnings=0 ;; esac
  if [ "$warnings" -ne 0 ]; then
    fail "$warnings compiler warning(s) in $what"
    grep -aE 'warning:' "$log" | head -8
    rm -f "$log"; return 1
  fi
  # ninja prints one `[n/m] ...` line per edge it actually runs, and only
  # `ninja: no work to do.` when it runs none.
  if [ "$fresh" -eq 1 ] && ! grep -qaE '^\[[0-9]+/[0-9]+\]' "$log"; then
    fail "$what was supposed to be a fresh build and compiled nothing — the directory was not removed, so this proved nothing"
    tail -3 "$log"
    rm -f "$log"; return 1
  fi
  pass "$what built warning-clean$([ "$fresh" -eq 1 ] && printf ', %s edges' "$(grep -caE '^\[[0-9]+/[0-9]+\]' "$log")")"
  rm -f "$log"; return 0
}

# -------------------------------------------------------------------- GPU -----
# GPU work must **skip, not fail**, when there is no Vulkan device -- the suite
# has to stay runnable in CI, which is where it is worth the most. The trap is
# that "skip" and "broken" then look identical, and for four tools they were:
# `no usable GPU` matched the expected pattern and came out as a green tick, on a
# machine with a working card. A driver that stopped enumerating devices would
# have taken every GPU check with it and the gate would have passed.
#
# So the machine is asked directly, and the answer decides which contract is in
# force. `VK_DRIVER_FILES` / `VK_ICD_FILENAMES` come first because the Vulkan
# loader honours them exclusively when set: pointing one at a nonexistent file
# genuinely deprives the process tree of Vulkan, which is how the deviceless
# branch below is tested on a machine that has a card.
GPU=absent
GPU_WHY="not yet determined"
detect_gpu() {
  local icd="" node="" list="" f d
  if   [ -n "${VK_DRIVER_FILES:-}" ];   then list="$VK_DRIVER_FILES"
  elif [ -n "${VK_ICD_FILENAMES:-}" ]; then list="$VK_ICD_FILENAMES"
  fi
  if [ -n "$list" ]; then
    local IFS=:
    for f in $list; do [ -r "$f" ] && { icd="$f"; break; }; done
  else
    for f in /usr/share/vulkan/icd.d/*.json /etc/vulkan/icd.d/*.json \
             "${HOME:-/nonexistent}"/.local/share/vulkan/icd.d/*.json; do
      [ -r "$f" ] && { icd="$f"; break; }
    done
  fi
  for d in /dev/dri/renderD* /dev/nvidia[0-9]*; do [ -e "$d" ] && { node="$d"; break; }; done

  if [ -n "$icd" ] && [ -n "$node" ]; then
    GPU=present
    GPU_WHY="this machine has a Vulkan driver ($icd) and a device node ($node)"
  elif [ -z "$icd" ]; then
    GPU=absent
    GPU_WHY="no readable Vulkan ICD manifest${list:+ among $list}"
  else
    GPU=absent
    GPU_WHY="a Vulkan driver ($icd) is installed but there is no device node (/dev/dri/renderD*, /dev/nvidia*)"
  fi
}

# What a tool says when it found no device. Deliberately covers the raw loader
# error as well as each tool's own wording, because `fem_spike` has no skip path
# at all: it runs its checks against a device that is not there and prints
# `CHECKS FAILED`. On a deviceless machine that is a skip; on this one it is a
# failure, and the difference is the whole point of this block.
NO_DEVICE_RE='no usable GPU|^skipped: |vkCreateInstance failed|VK_ERROR_INCOMPATIBLE_DRIVER|no Vulkan (loader|driver|device)'

gpu_ok() {   # <label> <success pattern> <binary> [args...]
  local label="$1" expect="$2"; shift 2
  local bin="$1" glslc="glslc is installed"
  command -v glslc >/dev/null 2>&1 || glslc="glslc is not installed"
  if [ ! -x "$bin" ]; then
    if [ "$GPU" = present ] && [ "$glslc" = "glslc is installed" ]; then
      fail "$label — $bin was not built, yet $GPU_WHY and $glslc: that is a build failure, not a machine without a device"
      return 1
    fi
    skip "$label — $bin not built ($GPU_WHY, $glslc)"
    return 0
  fi
  run_step "$label" "$STEP_TIMEOUT" "$@" || return 1
  if grep -qaE "$NO_DEVICE_RE" "$STEP_LOG"; then
    if [ "$GPU" = present ]; then
      fail "$label — reported no usable GPU ($(why "$STEP_STATUS" "$STEP_TIMEOUT")), but $GPU_WHY"
      diagnose "$STEP_LOG" "$STEP_STATUS" "$STEP_TIMEOUT"
      rm -f "$STEP_LOG"; return 1
    fi
    skip "$label — no Vulkan device on this machine ($GPU_WHY)"
    rm -f "$STEP_LOG"; return 0
  fi
  verdict "$label" "$expect" "$STEP_TIMEOUT"
  local rc=$?
  rm -f "$STEP_LOG"
  return $rc
}

# --------------------------------------------------------------- selftest -----
# The negative controls. Each drives one reporter with a failure it was **not**
# written for and asserts on what came back; a reporter that stays quiet, or that
# reports the wrong thing, fails the gate here rather than in six months when it
# is standing between a real defect and a commit.
#
# The one that has already earned its keep twice is the first: a step that exits
# **3** must be reported as 3. `if ! cmd; then ... $?` makes `$?` the status of
# the *negation* -- 0, or 1 -- and that was written into two functions here,
# hours apart, by someone who had already fixed it once. Reading the code did not
# catch it either time.
CONTROLS=0
control() {   # <what it proves> <regex the report must contain> <command...>
  local what="$1" want="$2"; shift 2
  CONTROLS=$((CONTROLS + 1))
  local out rc
  out="$("$@" 2>&1)"; rc=$?          # a subshell: the reporter's FAILURES++ is discarded
  if [ "$rc" -eq 0 ]; then
    fail "control [$what]: the reporter returned SUCCESS on a step that failed"
    printf '%s\n' "$out" | head -5
    return 1
  fi
  if ! printf '%s\n' "$out" | grep -qE "$want"; then
    fail "control [$what]: the report never said /$want/ — an unforeseen failure came out as silence"
    printf '%s\n' "$out" | head -8
    return 1
  fi
  return 0
}
control_pass() {   # a reporter that cannot go green is as broken as one that cannot go red
  local what="$1"; shift
  CONTROLS=$((CONTROLS + 1))
  local out rc
  out="$("$@" 2>&1)"; rc=$?
  if [ "$rc" -ne 0 ] || ! printf '%s\n' "$out" | grep -q '✓'; then
    fail "control [$what]: the reporter went RED on a step that succeeded"
    printf '%s\n' "$out" | head -5
    return 1
  fi
  return 0
}

selftest() {
  local sandbox
  sandbox="$(mktemp -d -t shipsim-gate-selftest.XXXXXX)" \
    || { fail "selftest: cannot create a sandbox in ${TMPDIR:-/tmp}"; return 1; }
  # Every log the controls make lands in the sandbox and leaves with it.
  local saved_tmpdir="${TMPDIR:-}"
  export TMPDIR="$sandbox"

  # --- expect_ok, driven by failures it was not written for -----------------
  control "a step that fails with no output at all" 'exited 3' \
          expect_ok ctl '' sh -c 'exit 3'
  control "a step killed by a signal" 'signal 11' \
          expect_ok ctl '' sh -c 'kill -SEGV $$'
  # `check-figures.sh` reports a drift as `  ✗ <label>: <doc> publishes X, the
  # tool now says Y` and matched none of `diagnose`'s markers, so a real drift
  # came out as `--- last lines ---` and eight lines of the passing tail. Twice.
  # This drives the same shape: a failure whose only evidence is a ✗ line, buried
  # under enough successes that `tail -8` cannot reach it.
  control "a failure whose only marker is a drift line, buried under successes" 'publishes' \
          expect_ok ctl '^ok' sh -c \
          'printf "  ✗ a figure: README.md publishes 41, the tool now says 42\n"; \
           for i in $(seq 20); do printf "  ✓ something else = %s\n" "$i"; done; exit 1'
  local saved_step=$STEP_TIMEOUT
  STEP_TIMEOUT=1
  control "a step that hangs" 'HUNG' expect_ok ctl '' sleep 30
  STEP_TIMEOUT=$saved_step
  control "a step that exits 0 without saying it passed" 'never reported' \
          expect_ok ctl 'all checks passed' sh -c 'echo working; echo done'
  # The pattern is a substring match unless it is anchored, and `0 failures`
  # is a substring of `10 failures`. This control is the reason the call sites
  # below all say `checks, 0 failures`.
  control "a suite that exits 0 while reporting 10 failures" 'never reported' \
          expect_ok ctl '[1-9][0-9]* checks, 0 failures' sh -c 'echo "4210 checks, 10 failures"'
  # **A gate that checked nothing has not passed.** `0 checks, 0 failures` is
  # what a suite whose runners were all compiled out prints, and it satisfied
  # `checks, 0 failures` exactly -- a degenerate run emitting the success string,
  # which is the fourth instance of that shape found in this repo in a week.
  control "a suite that exits 0 having run no checks at all" 'never reported' \
          expect_ok ctl '[1-9][0-9]* checks, 0 failures' sh -c 'echo "0 checks, 0 failures"'
  control_pass "a step that succeeds" \
          expect_ok ctl '[1-9][0-9]* checks, 0 failures' sh -c 'echo "4210 checks, 0 failures"'

  # **A figure gate that skipped most of its blocks has not passed either.**
  # `^ok — [1-9]` accepted `ok — 10 published figures ...` from a tree where only
  # `shipsim_tests` was built and twelve tool blocks never ran -- a tenth of the
  # coverage, reported in the same words as all of it. This is the `0 checks, 0
  # failures` shape one level up: not a degenerate run emitting the success
  # string, but a *partial* one. The pattern below is the one used against
  # `check-figures.sh` on a machine with a device, where nothing may be skipped.
  control "a figure gate that ran a tenth of its blocks" 'never reported' \
          expect_ok ctl '^ok — [1-9][0-9]* published figures still match the tool$' \
          sh -c 'echo "ok — 10 published figures still match the tool, and these did not run: ram_view zone_probe"'
  control_pass "a figure gate that skipped nothing" \
          expect_ok ctl '^ok — [1-9][0-9]* published figures still match the tool$' \
          sh -c 'echo "ok — 381 published figures still match the tool"'

  # **And `check()` rejected an empty parse while accepting a garbled one.**
  # `awk` coerces any non-numeric string to 0, so a check whose expected value is
  # zero -- there are six -- passed on `n/a`, on a stray unit, on an error
  # message. Driven here through the real script so the control breaks if that
  # predicate is ever loosened: a tool whose output the `sed` cannot parse into a
  # number must go red, not green.
  control "a figure gate handed a non-numeric parse" 'not a number' \
          sh -c 'cd '"$PWD"' && CTL_FIGURE=n/a ./scripts/check-figures.sh --selftest-parse'

  # --- build_into, against a compiler that dies without saying `error:` -----
  # The real case was a full disk. This is the same class without one: the
  # child is killed outright, ninja says nothing, and the reporter has only the
  # exit status to work with.
  printf '#!/bin/sh\ncase "${GATE_CONTROL_CMAKE:-ok}" in fail) echo "CMake Error at CMakeLists.txt:1"; exit 1 ;; esac\nexit 0\n' > "$sandbox/cmake"
  # `eatlog` unlinks the log it is writing to and then *succeeds*. `grep -c` on a
  # file that is no longer there prints nothing, `[ "" -ne 0 ]` is an error, and
  # an `if` reads that error as false -- so the gate used to report a build whose
  # warning count it could not determine as warning-clean. fd 2, not fd 1: inside
  # a command substitution fd 1 is the substitution's own pipe.
  # `noop` is the rebuild that was not clean: `ninja: no work to do.`, exit 0, no
  # warnings and not one file compiled.
  printf '#!/bin/sh\ncase "${GATE_CONTROL_NINJA:-ok}" in\n  kill) kill -9 $$ ;;\n  warn) echo "engine/sim/hull.cpp:1:1: warning: control [-Wall]"; exit 0 ;;\n  hang) sleep 30 ;;\n  quiet) exit 1 ;;\n  eatlog) rm -f "$(readlink /proc/self/fd/2)"; exit 0 ;;\n  noop) echo "ninja: no work to do."; exit 0 ;;\nesac\necho "[1/1] Building CXX object control.cpp.o"\nexit 0\n' > "$sandbox/ninja"
  chmod +x "$sandbox/cmake" "$sandbox/ninja"

  # PATH rather than `env`, because `env` cannot run a shell function and these
  # controls have to drive the real `build_into`, not a copy of it.
  # PATH rather than `env`, and an explicit export rather than an assignment
  # prefix, because `env` cannot run a shell function and a prefix on a function
  # call is not reliably exported to that function's own children.
  local saved_path="$PATH"
  PATH="$sandbox:$PATH"
  export GATE_CONTROL_CMAKE=ok GATE_CONTROL_NINJA=ok

  GATE_CONTROL_NINJA=kill
  control "a build whose compiler is killed outright" 'signal 9' build_into "$sandbox/b1"
  GATE_CONTROL_NINJA=quiet
  control "a build that fails without the word error:" 'exited 1' build_into "$sandbox/b2"
  GATE_CONTROL_NINJA=warn
  control "a compiler warning" 'warning' build_into "$sandbox/b3"
  GATE_CONTROL_NINJA=eatlog
  control "a build whose log vanished before it could be read" 'disappeared' \
          build_into "$sandbox/b3a"
  GATE_CONTROL_NINJA=ok GATE_CONTROL_CMAKE=fail
  control "a configure that fails" 'configure failed' build_into "$sandbox/b4"
  GATE_CONTROL_CMAKE=ok GATE_CONTROL_NINJA=hang
  local saved_build=$BUILD_TIMEOUT
  BUILD_TIMEOUT=1
  control "a build that hangs" 'HUNG' build_into "$sandbox/b5"
  BUILD_TIMEOUT=$saved_build
  GATE_CONTROL_NINJA=noop
  control "a clean rebuild that compiled nothing" 'compiled nothing' \
          build_into "$sandbox/b7" --fresh
  # ...and the same no-op build without `--fresh`, which is an *incremental*
  # build with nothing to do and must stay green. A guard that fires on both is
  # not distinguishing anything.
  control_pass "an incremental build with nothing to do" build_into "$sandbox/b8"
  GATE_CONTROL_NINJA=ok
  control_pass "a build that succeeds" build_into "$sandbox/b6"
  control_pass "a fresh build that did compile something" build_into "$sandbox/b9" --fresh

  unset GATE_CONTROL_CMAKE GATE_CONTROL_NINJA
  PATH="$saved_path"

  # --- the GPU contract: a broken driver must not look like a missing one ---
  printf '#!/bin/sh\necho "  no usable GPU (vkCreateInstance failed)"\nexit 0\n' > "$sandbox/nogpu"
  chmod +x "$sandbox/nogpu"
  local saved_gpu=$GPU saved_why=$GPU_WHY
  GPU=present; GPU_WHY="a device node exists"
  control "a GPU tool reporting no device on a machine that has one" 'reported no usable GPU' \
          gpu_ok ctl '^ok$' "$sandbox/nogpu"
  control "a GPU tool missing on a machine that has a device" 'build failure' \
          gpu_ok ctl '^ok$' "$sandbox/never-built"
  # ...and the other direction, which must stay a skip and must stay green: a
  # machine with no device is not a red gate, or nobody runs the gate in CI.
  GPU=absent; GPU_WHY="no ICD manifest"
  CONTROLS=$((CONTROLS + 1))
  local out rc
  out="$(gpu_ok ctl '^ok$' "$sandbox/nogpu" 2>&1)"; rc=$?
  if [ "$rc" -ne 0 ] || ! printf '%s\n' "$out" | grep -q 'no Vulkan device'; then
    fail "control [a machine with no device]: must skip and stay green, not fail"
    printf '%s\n' "$out" | head -5
  fi
  GPU=$saved_gpu; GPU_WHY=$saved_why

  # --- and the counter itself, since every verdict above is built on it ------
  # `fail` inside `$( )` or a pipeline increments a *subshell's* copy and the
  # increment is thrown away. That is how a counted failure goes missing, so both
  # halves are asserted: the subshell must not change the count, and the plain
  # call must.
  CONTROLS=$((CONTROLS + 1))
  local before=$FAILURES
  ( fail "control (expected: discarded)" ) >/dev/null
  if [ "$FAILURES" -ne "$before" ]; then
    printf '  %s✗ control [subshell]: a fail() in a subshell changed the parent count%s\n' "$RED" "$OFF"
    FAILURES=$((before + 1)); before=$FAILURES
  fi
  fail "control (expected: unwound)" >/dev/null
  if [ "$FAILURES" -ne $((before + 1)) ]; then
    printf '  %s✗ control [counter]: fail() did not increment FAILURES%s\n' "$RED" "$OFF"
    FAILURES=$((before + 1))
  else
    FAILURES=$before      # unwind the deliberate one
  fi

  rm -rf "$sandbox"
  if [ -n "$saved_tmpdir" ]; then export TMPDIR="$saved_tmpdir"; else unset TMPDIR; fi
  return 0
}

section "the gate's own failure paths (negative controls)"
before_controls=$FAILURES
selftest
if [ "$CONTROLS" -eq 0 ]; then
  fail "the negative controls did not run at all — the reporters below are unverified"
elif [ "$FAILURES" -eq "$before_controls" ]; then
  pass "$CONTROLS controls: every reporter went red on a failure it was not written for"
else
  fail "$((FAILURES - before_controls)) of $CONTROLS negative controls failed — nothing below this line can be trusted"
fi

if [ "$LEVEL" = "selftest" ]; then
  section "result"
  verdict_line "selftest"
  exit $(( FAILURES > 0 ))
fi

# ---------------------------------------------------------------- quick -------
section "build and unit tests"
if build_into build; then
  # `[1-9][0-9]*` and not `[0-9]*`. `0 checks, 0 failures` is what a suite that
  # ran nothing prints -- a `#if` that took out every runner, a `main` that
  # returned early -- and it satisfied `checks, 0 failures` exactly. **A gate
  # that checked nothing has not passed**: the same hole, in the same week, as
  # `check-figures.sh` reporting `ok — 0 published figures`.
  expect_ok "shipsim_tests" '[1-9][0-9]* checks, 0 failures' ./build/shipsim_tests
else
  # Not silent. The `build && test` chain used to leave no line at all for the
  # test suite, so a reader had to notice which check was *absent*.
  skip "shipsim_tests not run: the build step did not report success"
fi

if [ "$LEVEL" = "quick" ]; then
  section "result"
  verdict_line "quick gate"
  exit $(( FAILURES > 0 ))
fi

# ----------------------------------------------------------------- full -------
section "clean rebuild (an incremental build cannot see a warning it did not recompile)"
rm -rf build-verify
build_into build-verify --fresh
rm -rf build-verify

# ...and a build cannot see a warning its optimisation level does not produce.
#
# Everything above compiles `RelWithDebInfo`. GCC's `-Waggressive-loop-optimizations`
# needs `-O3`, and it was sitting on `reduction.cpp`'s triangular solves for as long
# as they existed: with a signed loop counter it cannot bound the trip count, so it
# reasons `n` might be `INT_MAX` and the indexing would then be undefined. Warnings
# are failures here, so a warning no gate ever compiles is a blind spot rather than
# a curiosity.
#
# Engine only, and that is a deliberate trade: it is where the numerics live and
# where the optimiser has something to reason about, and it costs ~6 s where the
# whole tree would cost far more. Tools and tests are not covered at -O3.
section "optimiser warnings (-O3 sees what RelWithDebInfo does not)"
rm -rf build-o3
build_into build-o3 --target shipsim_engine --fresh \
           -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-O3"
rm -rf build-o3

section "GPU tools"
detect_gpu
note "Vulkan: $GPU — $GPU_WHY"
if [ "$GPU" = present ]; then
  note "so a tool that reports no device, or was never built, is a failure here and not a skip"
fi
gpu_ok "fem_spike" 'all checks passed' ./build/fem_spike
gpu_ok "ferry_view" 'all checks passed' \
       ./build/ferry_view --out=/tmp --frames=3 --duration=600
# The integration case: hull from particulars, radiation, propulsion, spectral
# sea, ocean and hull rendering, all in one run. Unit tests cover each of those
# alone; this is the only thing that runs them against each other.
# The Phase 3 milestone as one act: collision, indentation, breach and flooding
# in a single run. Every one of those is unit-tested alone; this is the only thing
# that runs them against each other.
gpu_ok "ram_view" '^ok$' ./build/ram_view --speed=4.0 --duration=120
# And that the numbers docs/06-roadmap.md publishes about this milestone are
# still the numbers the tool produces. Docs are the source of truth here and
# nothing enforced it: the previous set went stale when a wing tank authored
# inside a hold stopped flooding twice, and a claim in indentation.hpp outlived
# the code that contradicted it by weeks. ~330 s, and it names the line to fix.
#
# **Not gated on ram_view.** It was, and ram_view is a Vulkan-only binary, so on
# any machine without a device every figure this checks -- including the ones
# that come from CPU-only tools -- went unchecked with no line printed at all.
# That is the same shape as the `exit 0` this script's sibling had: a check that
# could not report the failure it existed to catch. What it can and cannot
# re-derive on a deviceless machine is `check-figures.sh`'s own business to say.
# The `[1-9]` is not decoration. `check-figures.sh` guards every tool it drives
# with `[ -x ]` and skips the block when one is missing, so with nothing built at
# all it prints `ok — 0 published figures still match the tool` and exits 0. That
# is a green verdict on an empty set, and `^ok — ` accepted it. Measured, not
# supposed: with the GPU tools hidden it still checks 78 of them, so any number
# the gate should accept starts with a digit other than zero.
#
# **And `[1-9]` closes `ok — 0` and nothing above it.** A tree with only
# `shipsim_tests` built runs two blocks, skips the other twelve, prints
# `ok — 10 published figures ...` and passes: a green verdict on a tenth of the
# coverage, which reads exactly like a green verdict on all of it. A floor on the
# count would be a guess that goes stale every time a block is added. What
# actually needs asserting is the **skipped set**, which `check-figures.sh`
# already prints and which it deliberately does not judge -- whether a skip is
# legitimate is a question about *this machine*, and this script is the one that
# knows. Same division of labour as `gpu_ok` above.
#
# With a device present nothing may be skipped at all, which the `$` anchor says:
# the skipped form of that line continues ", and these did not run: ...".
if [ "$GPU" = present ]; then
  figures_ok='^ok — [1-9][0-9]* published figures still match the tool$'
else
  figures_ok='^ok — [1-9]'
fi
require_built ./scripts/check-figures.sh &&
  expect_ok "published figures" "$figures_ok" ./scripts/check-figures.sh
# The Tier-2 zone at ship scale: solid-shell elements over the ferry's own
# plating, driven to tearing, and the torn panels fed to breach. The unit suite
# tests the pieces at unit scale because a real solve is core-minutes; this is the
# only thing that runs the whole chain on a real ship, and it is why the zone's
# cost has to stay visible. No GPU: a missing binary is a broken build.
require_built ./build/zone_probe &&
  expect_ok "zone_probe" '^ok$' ./build/zone_probe --radius=3.0 --depth=0.22
# The same zone on the GPU against the CPU double reference. **Every figure in
# `07-fem-spike-findings.md` §8 comes out of this tool and nothing ran it**, so it
# could rot silently between the measurements that cite it — which is the shape of
# the wing-tank warning that printed on every run for months with nobody reading it.
# 1 500 steps, and the number was measured rather than picked. The tool refuses to
# report `ok` unless the patch both moved and yielded — its own guard against
# comparing two solvers that did nothing — and 600 steps does not yield it: peak
# equivalent plastic strain 0.0000, and the guard fires exactly as it should. 1 200
# is the first count that yields; 1 500 gives margin and still costs under a second,
# against the ~30 s a run to tearing would.
gpu_ok "zone_gpu_probe" '^ok$' \
       ./build/zone_gpu_probe --radius=2.5 --sub=4 --steps=1500
# The Phase 4 milestone as one act: a two-zone fire against the ferry's own engine
# room bulkhead, an implicit conduction solve on that bulkhead, EN 1993-1-2 steel
# under restrained expansion *and* the head of water behind it, and the panels that
# go handed to `breach` and on into the flooding network. Every piece is unit-tested
# alone; this is the only thing that runs the chain on a real ship.
#
# **And it is the only thing that runs the milestone's own sentence as an acceptance
# test.** "Fails under the head of water behind it" is a claim that neither cause is
# sufficient, so the tool runs the case three times -- fire with a dry hold, water
# with no fire, both -- and refuses `ok` unless the two controls survive and the pair
# does not. A fourth pass runs the pair with the damage evaluated but not applied,
# because the failure relieves what caused it and a window measured off a relieved
# run comes back equal to whatever restraint that run was given. ~45 s.
require_built ./build/bulkhead_probe &&
  expect_ok "bulkhead_probe" '^ok$' ./build/bulkhead_probe --quiet
# The Tier-1 section at ship scale. `tests/test_section.cpp` checks the mesher
# against closed forms at unit scale and on two frame bays; this is the only thing
# that asks the whole hull, and **every ship-scale figure docs/02-simulation.md §3
# publishes about the section mesher comes out of it**.
#
# It was not in the gate until it published a 3.46 Hz shell frequency that nothing
# had ever re-run — the same failure `zone_gpu_probe` was added for, and the third
# time this session that a document quoted a tool nobody executed. It had no
# success contract at all either: it exited 0 whatever it measured.
#
# Two runs, 45 s together. `--scan=2` is the **reach**: a two-bay window slid along
# the whole ship, every one of which has to mesh, reduce and solve. The default run
# is what the reach is *for* — that the junction tie closes the cell — and it is
# checked on torsion, on the lowest fixed-interface frequency and on the component
# count, because `section.hpp` §2 shows EA and EI are exact on a section whose
# plating is joined to nothing. `--no-reduce` skips the Craig-Bampton block, which
# is 70 s and which the unit suite already covers at unit scale.
#
# The third run is the **invariance** and it is 2 s: at every frame station of the
# ship, the window aft of it, the window forward of it and the window spanning both,
# asked whether they agree about the station they share. They did not — up to 1.0 mm
# of node displacement and 25.2% of the stiffener steel depended on where the section
# was cut — and the unit suite can only afford two stations. This sweeps 47, with the
# `halo = false` control alongside so that "every station agreed" cannot be reported
# by a mesher that was never fixed.
if require_built ./build/section_probe; then
  expect_ok "section_probe reach"      '^ok$' ./build/section_probe --scan=2
  expect_ok "section_probe invariance" '^ok$' ./build/section_probe --invariance=2
  expect_ok "section_probe tie"        '^ok$' ./build/section_probe --sweep=0 --no-reduce
fi
# The volumetric fire and smoke pass on the ferry. Its own closed-form checks live
# in the unit suite; what this adds is the whole chain on a real casualty — a
# two-zone fire model driving a renderer that must not invent anything it does not
# have — and every figure `docs/03-renderer-audio.md`'s fire section publishes.
# It asserts the *absence* of a glow at 4 MW, which is the finding that would go
# stale first if someone reached for the exposure knob.
gpu_ok "smoke_view" '^ok$' ./build/smoke_view --out=/tmp --frames=8 --duration=600
gpu_ok "seaway_view" '^ok$' ./build/seaway_view --out=/tmp --frames=2

section "flooding scenarios still behave"
# Run twice: once from the compiled reference ferry, once from ships/ferry.ship.
# The outcomes must be identical strings, which is the check that the ship
# definition format captured what actually matters rather than what was easy to
# serialise. The unit suite compares the two ships in lockstep too, but only for
# the first 150 s -- these are the 900 s runs the outcomes are quoted from, and
# 'none' does not reach its verdict until t+880.
#
# **What this block cannot see, by construction**: a change in the physics that
# moves both runs the same way. It compares the two paths against each other and
# never against a published value, so a ferry that stopped capsizing altogether
# would pass here in identical strings. The absolute outcomes are `check-figures.sh`'s
# job — which is the other reason it must not be gated on a GPU binary.
#
# **The run's exit status is part of the outcome.** This block used to be
# `... | grep '^=== Outcome' || true`, which threw away the status of a pipeline
# and of the run inside it: a scenario that printed its verdict and then died on
# the way out passed, twice, in identical strings. It threw away stderr too, on a
# ship whose missing wing tanks were announced by a warning nobody read for months.
outcome_of() {   # <label> <cmd...>; sets OUTCOME
  OUTCOME=""
  run_step "$1" "$SCENARIO_TIMEOUT" "${@:2}" || return 1
  if [ "$STEP_STATUS" -ne 0 ]; then
    fail "$1 — $(why "$STEP_STATUS" "$SCENARIO_TIMEOUT")"
    diagnose "$STEP_LOG" "$STEP_STATUS" "$SCENARIO_TIMEOUT"
    rm -f "$STEP_LOG"; return 1
  fi
  OUTCOME="$(grep -aE '^=== Outcome' "$STEP_LOG" | tail -1)"
  if [ -z "$OUTCOME" ]; then
    fail "$1 exited 0 but never printed an outcome line"
    diagnose "$STEP_LOG" "$STEP_STATUS" "$SCENARIO_TIMEOUT"
    rm -f "$STEP_LOG"; return 1
  fi
  # A marker with no verdict behind it compares equal to the same marker from the
  # other run, and the pair passes on two empty strings. The comparison below is
  # only worth anything if there is something in it to compare.
  local verdict="${OUTCOME#=== Outcome}"
  if [ -z "${verdict//[[:space:]=]/}" ]; then
    fail "$1 printed an outcome marker with no verdict behind it: '$OUTCOME'"
    diagnose "$STEP_LOG" "$STEP_STATUS" "$SCENARIO_TIMEOUT"
    rm -f "$STEP_LOG"; return 1
  fi
  rm -f "$STEP_LOG"; return 0
}
for scenario in none doors full; do
  outcome_of "$scenario (compiled ferry)" \
             ./build/shipsim --scenario="$scenario" --duration=900 || continue
  compiled="$OUTCOME"
  outcome_of "$scenario (ships/ferry.ship)" \
             ./build/shipsim --ship=ships/ferry.ship --scenario="$scenario" --duration=900 || continue
  if [ "$compiled" != "$OUTCOME" ]; then
    fail "$scenario: ships/ferry.ship diverges from the compiled ferry"
    printf '      compiled: %s\n      file:     %s\n' "$compiled" "$OUTCOME"
  else
    pass "$scenario (compiled and from file): ${compiled#=== Outcome }"
  fi
done

section "repeat runs (concurrency must not be flaky)"
# **Keep the evidence.** This step used to discard every run's output, so a failure
# here said "1 of 6 runs failed" and nothing else -- which is the one thing a flake
# report must not do, because a flake that cannot be named cannot be told from a
# regression and the next run is the only chance to see it. The exit status
# separates a crash from a failed assertion from a **hang**, and the failing run's
# own `FAIL` lines say which check went. The bound matters as much as the report:
# without it a suite that stopped rather than failed took the whole gate with it,
# silently, for as long as anyone was willing to wait.
flaky=0
for run in 1 2 3 4 5 6; do
  if ! run_step "repeat $run" "$STEP_TIMEOUT" ./build/shipsim_tests; then
    flaky=$((flaky + 1)); continue
  fi
  # Status *and* content. A run that exits 0 having checked nothing is the worst
  # flake there is, and this loop used to ask only whether the process came back.
  if [ "$STEP_STATUS" -eq 0 ] && grep -qaE '[1-9][0-9]* checks, 0 failures' "$STEP_LOG"; then
    rm -f "$STEP_LOG"; continue
  fi
  flaky=$((flaky + 1))
  if [ "$STEP_STATUS" -eq 0 ]; then
    printf '      run %d exited 0 but reported no completed checks\n' "$run"
  else
    printf '      run %d %s\n' "$run" "$(why "$STEP_STATUS" "$STEP_TIMEOUT")"
  fi
  grep -aE '^  FAIL' "$STEP_LOG" | head -5
  tail -1 "$STEP_LOG"   # the check/failure count, or the last line before a crash
  rm -f "$STEP_LOG"
done
[ "$flaky" -eq 0 ] && pass "6/6 stable" || fail "$flaky of 6 runs failed"

if [ "$LEVEL" = "full" ]; then
  section "result"
  verdict_line "full gate"
  exit $(( FAILURES > 0 ))
fi

# ------------------------------------------------------------- sanitize -------
# Both sanitizers have found real bugs the functional tests passed straight
# through: a slot-reuse data race in the job queue, and -- not a race at all --
# an off-by-one that only appeared because TSan's instrumentation made the code
# slow enough to reach a clamp the fast build never hit.
section "AddressSanitizer + UBSan"
if build_into build-asan -DSHIPSIM_SANITIZE=address; then
  expect_ok "shipsim_tests (asan)" '[1-9][0-9]* checks, 0 failures' ./build-asan/shipsim_tests
else
  skip "shipsim_tests (asan) not run: the build-asan build failed"
fi

section "ThreadSanitizer"
if build_into build-tsan -DSHIPSIM_SANITIZE=thread; then
  # halt_on_error=0 so one race does not hide the next; the count below is then
  # the number of distinct reports rather than the first one.
  export TSAN_OPTIONS="halt_on_error=0"
  run_step "shipsim_tests (tsan)" "$STEP_TIMEOUT" ./build-tsan/shipsim_tests
  unset TSAN_OPTIONS
  races="$(grep -acE '^SUMMARY: ThreadSanitizer' "$STEP_LOG")"
  warned="$(grep -acE '^WARNING: ThreadSanitizer' "$STEP_LOG")"
  case "$races"  in ''|*[!0-9]*) races=0 ;; esac
  case "$warned" in ''|*[!0-9]*) warned=0 ;; esac
  [ "$warned" -gt "$races" ] && races="$warned"   # a report whose SUMMARY never printed
  if [ "$races" -ne 0 ]; then
    fail "$races ThreadSanitizer report(s)"
    grep -aE '^(WARNING|SUMMARY): ThreadSanitizer' "$STEP_LOG" | head -5
  # **The exit status is part of the answer.** It was not looked at here at all,
  # so a TSan run that printed its summary and then died on a signal -- or that
  # hung and was killed after printing it -- came out as a green tick.
  elif [ "$STEP_STATUS" -ne 0 ]; then
    fail "shipsim_tests (tsan) — $(why "$STEP_STATUS" "$STEP_TIMEOUT")"
    diagnose "$STEP_LOG" "$STEP_STATUS" "$STEP_TIMEOUT"
  elif ! grep -qaE '[1-9][0-9]* checks, 0 failures' "$STEP_LOG"; then
    fail "shipsim_tests (tsan) reported failures, or ran no checks at all"
    diagnose "$STEP_LOG" "$STEP_STATUS" "$STEP_TIMEOUT"
  else
    pass "shipsim_tests (tsan) — $(grep -oaE '[0-9]+ checks, [0-9]+ failures' "$STEP_LOG" | tail -1), 0 races"
  fi
  rm -f "$STEP_LOG"
else
  skip "shipsim_tests (tsan) not run: the build-tsan build failed"
fi

section "result"
verdict_line "$LEVEL gate"
exit $(( FAILURES > 0 ))
