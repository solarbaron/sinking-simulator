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
if cmake -S . -B build-o3 -G Ninja -DCMAKE_BUILD_TYPE=Release \
         -DCMAKE_CXX_FLAGS="-O3" >/dev/null 2>&1; then
  o3log="$(mktemp)"
  if ninja -C build-o3 shipsim_engine >"$o3log" 2>&1; then
    o3warnings="$(grep -cE 'warning:' "$o3log" || true)"
    [ "$o3warnings" -eq 0 ] && pass "engine built warning-clean at -O3" \
                            || { fail "$o3warnings warning(s) at -O3"; grep -E 'warning:' "$o3log" | head -5; }
  else
    fail "the -O3 engine build failed"
    grep -E 'error:' "$o3log" | head -5
  fi
  rm -f "$o3log"
else
  fail "cmake configure failed for build-o3"
fi
rm -rf build-o3

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
# The Phase 3 milestone as one act: collision, indentation, breach and flooding
# in a single run. Every one of those is unit-tested alone; this is the only thing
# that runs them against each other.
if [ -x ./build/ram_view ]; then
  expect_ok "ram_view" '^ok$' ./build/ram_view --speed=4.0 --duration=120
  # And that the numbers docs/06-roadmap.md publishes about this milestone are
  # still the numbers the tool produces. Docs are the source of truth here and
  # nothing enforced it: the previous set went stale when a wing tank authored
  # inside a hold stopped flooding twice, and a claim in indentation.hpp outlived
  # the code that contradicted it by weeks. ~100 s, and it names the line to fix.
  expect_ok "published figures" '^ok — ' ./scripts/check-figures.sh
else
  skip "ram_view not built"
fi
# The Tier-2 zone at ship scale: solid-shell elements over the ferry's own
# plating, driven to tearing, and the torn panels fed to breach. The unit suite
# tests the pieces at unit scale because a real solve is core-minutes; this is the
# only thing that runs the whole chain on a real ship, and it is why the zone's
# cost has to stay visible.
if [ -x ./build/zone_probe ]; then
  expect_ok "zone_probe" '^ok$' ./build/zone_probe --radius=3.0 --depth=0.22
else
  skip "zone_probe not built"
fi
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
if [ -x ./build/zone_gpu_probe ]; then
  expect_ok "zone_gpu_probe" '^ok$|^skipped: ' \
            ./build/zone_gpu_probe --radius=2.5 --sub=4 --steps=1500
else
  skip "zone_gpu_probe not built"
fi
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
if [ -x ./build/section_probe ]; then
  expect_ok "section_probe reach" '^ok$' ./build/section_probe --scan=2
  expect_ok "section_probe tie"   '^ok$' ./build/section_probe --sweep=0 --no-reduce
else
  skip "section_probe not built"
fi
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
