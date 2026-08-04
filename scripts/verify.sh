#!/usr/bin/env bash
#
# The verification gate.
#
# Everything this project claims about itself is checked here. Each level is a
# superset of the one before:
#
#   quick     configure, build warning-clean, run the test suite      (~1 min)
#   full      + fresh rebuild, GPU tools, flooding scenarios, repeats (~5 min)
#   sanitize  + AddressSanitizer and ThreadSanitizer builds           (~10 min)
#   all       everything above
#
# Warnings are failures, not notes. The build is -Wall -Wextra -Wpedantic and has
# been warning-clean since Phase 0; the moment that stops being true it stops
# being a useful signal, so the gate treats a single warning as a red build.
#
# Incremental builds only recompile what changed, so they cannot see a warning in
# a file nobody touched. `full` and above therefore configure a throwaway build
# directory and compile everything from scratch.
set -uo pipefail

readonly ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

LEVEL="${1:-quick}"
case "$LEVEL" in
  quick|full|sanitize|all) ;;
  -h|--help|help)
    sed -n '2,20p' "${BASH_SOURCE[0]}" | sed 's/^# \?//'
    exit 0 ;;
  *) echo "usage: $0 [quick|full|sanitize|all]" >&2; exit 2 ;;
esac

readonly BOLD=$'\033[1m' RED=$'\033[0;31m' GREEN=$'\033[0;32m' YELLOW=$'\033[0;33m' OFF=$'\033[0m'
FAILURES=0
STARTED=$SECONDS

section() { printf '\n%s== %s ==%s\n' "$BOLD" "$1" "$OFF"; }
pass()    { printf '  %s✓%s %s\n' "$GREEN" "$OFF" "$1"; }
fail()    { printf '  %s✗ %s%s\n' "$RED" "$1" "$OFF"; FAILURES=$((FAILURES + 1)); }
skip()    { printf '  %s—%s %s\n' "$YELLOW" "$OFF" "$1"; }

# Build into `dir`, failing on any compiler warning. Extra cmake args follow.
build_into() {
  local dir="$1"; shift
  local log; log="$(mktemp)"
  if ! cmake -S . -B "$dir" -G Ninja "$@" >"$log" 2>&1; then
    fail "cmake configure failed for $dir"
    tail -20 "$log"; rm -f "$log"; return 1
  fi
  if ! ninja -C "$dir" >"$log" 2>&1; then
    fail "build failed in $dir"
    grep -E 'error:' "$log" | head -10; rm -f "$log"; return 1
  fi
  local warnings; warnings="$(grep -cE 'warning:' "$log" || true)"
  if [ "$warnings" -ne 0 ]; then
    fail "$warnings compiler warning(s) in $dir"
    grep -E 'warning:' "$log" | head -8
    rm -f "$log"; return 1
  fi
  pass "$dir built warning-clean"
  rm -f "$log"; return 0
}

# Run a command, requiring success and (optionally) a pattern in its output.
expect_ok() {
  local label="$1" expect="$2"; shift 2
  local log; log="$(mktemp)"
  if ! timeout 1800 "$@" >"$log" 2>&1; then
    fail "$label exited non-zero"
    tail -15 "$log"; rm -f "$log"; return 1
  fi
  if [ -n "$expect" ] && ! grep -qE "$expect" "$log"; then
    fail "$label did not report '$expect'"
    tail -15 "$log"; rm -f "$log"; return 1
  fi
  pass "$label — $(grep -oE '[0-9]+ checks, [0-9]+ failures' "$log" | tail -1 || echo ok)"
  rm -f "$log"; return 0
}

# ---------------------------------------------------------------- quick -------
section "build and unit tests"
build_into build && expect_ok "shipsim_tests" '0 failures' ./build/shipsim_tests

if [ "$LEVEL" = "quick" ]; then
  section "result"
  [ "$FAILURES" -eq 0 ] && pass "quick gate passed in $((SECONDS - STARTED))s" \
                        || fail "$FAILURES failure(s)"
  exit $(( FAILURES > 0 ))
fi

# ----------------------------------------------------------------- full -------
section "clean rebuild (an incremental build cannot see a warning it did not recompile)"
rm -rf build-verify
build_into build-verify
rm -rf build-verify

section "GPU tools"
# These must skip, not fail, on a machine without Vulkan — the suite has to stay
# runnable in CI, which is where it is worth the most.
if [ -x ./build/fem_spike ]; then
  expect_ok "fem_spike" 'all checks passed' ./build/fem_spike
else
  skip "fem_spike not built (no Vulkan or no glslc)"
fi
if [ -x ./build/ferry_view ]; then
  expect_ok "ferry_view" 'all checks passed|no usable GPU' \
            ./build/ferry_view --out=/tmp --frames=3 --duration=600
else
  skip "ferry_view not built"
fi
# The integration case: hull from particulars, radiation, propulsion, spectral
# sea, ocean and hull rendering, all in one run. Unit tests cover each of those
# alone; this is the only thing that runs them against each other.
if [ -x ./build/seaway_view ]; then
  expect_ok "seaway_view" '^ok$|no usable GPU' \
            ./build/seaway_view --out=/tmp --frames=2
else
  skip "seaway_view not built"
fi

section "flooding scenarios still behave"
# Run twice: once from the compiled reference ferry, once from ships/ferry.ship.
# The outcomes must be identical strings, which is the check that the ship
# definition format captured what actually matters rather than what was easy to
# serialise. The unit suite compares the two ships in lockstep too, but only for
# the first 150 s -- these are the 900 s runs the outcomes are quoted from, and
# 'none' does not reach its verdict until t+880.
for scenario in none doors full; do
  out="$(timeout 900 ./build/shipsim --scenario="$scenario" --duration=900 2>/dev/null \
         | grep -E '^=== Outcome' || true)"
  from_file="$(timeout 900 ./build/shipsim --ship=ships/ferry.ship --scenario="$scenario" \
               --duration=900 2>/dev/null | grep -E '^=== Outcome' || true)"
  if [ -z "$out" ]; then fail "$scenario produced no outcome"
  elif [ "$out" != "$from_file" ]; then
    fail "$scenario: ships/ferry.ship diverges from the compiled ferry"
    printf '      compiled: %s\n      file:     %s\n' "$out" "$from_file"
  else pass "$scenario (compiled and from file): ${out#=== Outcome }"; fi
done

section "repeat runs (concurrency must not be flaky)"
flaky=0
for _ in $(seq 1 6); do ./build/shipsim_tests >/dev/null 2>&1 || flaky=$((flaky + 1)); done
[ "$flaky" -eq 0 ] && pass "6/6 stable" || fail "$flaky of 6 runs failed"

if [ "$LEVEL" = "full" ]; then
  section "result"
  [ "$FAILURES" -eq 0 ] && pass "full gate passed in $((SECONDS - STARTED))s" \
                        || fail "$FAILURES failure(s)"
  exit $(( FAILURES > 0 ))
fi

# ------------------------------------------------------------- sanitize -------
# Both sanitizers have found real bugs the functional tests passed straight
# through: a slot-reuse data race in the job queue, and -- not a race at all --
# an off-by-one that only appeared because TSan's instrumentation made the code
# slow enough to reach a clamp the fast build never hit.
section "AddressSanitizer + UBSan"
if build_into build-asan -DSHIPSIM_SANITIZE=address; then
  expect_ok "shipsim_tests (asan)" '0 failures' ./build-asan/shipsim_tests
fi

section "ThreadSanitizer"
if build_into build-tsan -DSHIPSIM_SANITIZE=thread; then
  log="$(mktemp)"
  TSAN_OPTIONS="halt_on_error=0" timeout 1800 ./build-tsan/shipsim_tests >"$log" 2>&1
  races="$(grep -cE '^SUMMARY: ThreadSanitizer' "$log" || true)"
  if [ "$races" -ne 0 ]; then
    fail "$races ThreadSanitizer report(s)"
    grep -E '^SUMMARY' "$log" | head -5
  elif ! grep -qE '0 failures' "$log"; then
    fail "shipsim_tests (tsan) reported failures"
    grep -E 'FAIL' "$log" | head -8
  else
    pass "shipsim_tests (tsan) — $(grep -oE '[0-9]+ checks, [0-9]+ failures' "$log" | tail -1), 0 races"
  fi
  rm -f "$log"
fi

section "result"
if [ "$FAILURES" -eq 0 ]; then
  pass "$LEVEL gate passed in $((SECONDS - STARTED))s"
else
  fail "$FAILURES failure(s) in $((SECONDS - STARTED))s"
fi
exit $(( FAILURES > 0 ))
