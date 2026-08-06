#!/usr/bin/env bash
# Every number docs/06-roadmap.md publishes about the Phase 3 milestone, checked
# against the tool that produces them.
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
# Cost: the damage figures need no flooding at all and run with `--duration=1` in
# about a second each; three runs pay for the full 900 s. Around 105 s against the
# ~580 s the `full` gate already costs, and the same budget the first version
# spent on figures that said less.
set -u

RAM=${RAM:-./build/ram_view}
SECTION=${SECTION:-./build/section_probe}
SECTION_DOC=docs/02-simulation.md
DOC=docs/06-roadmap.md
fails=0
checks=0

red=$'\033[0;31m'; green=$'\033[0;32m'; dim=$'\033[2m'; off=$'\033[0m'

# check <label> <expected> <tolerance> <actual> <doc quote to find>
#
# The tolerance is relative for magnitudes and the caller passes an absolute one
# for angles by giving a tolerance larger than 1. It is deliberately tight: these
# runs are deterministic, so a figure that moves at all is either a real change in
# the physics or a real change in the ship, and both want a human to look.
check() {
  local label="$1" expect="$2" tol="$3" actual="$4" quote="$5" where="${6:-$DOC}"
  checks=$((checks + 1))
  if [ -z "$actual" ]; then
    printf '  %s✗%s %s — could not parse the figure out of the tool output\n' "$red" "$off" "$label"
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

if [ ! -x "$RAM" ]; then
  echo "  - ram_view not built, skipping published-figure check"
  exit 0
fi

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
  check "windows in one piece"        46 0 "$onepiece" "46 of 49" "$SECTION_DOC"
  check "reach along the hull (m)"    120.0 0.05 "$reach" "120.0 m of 120.0 m" "$SECTION_DOC"
else
  echo "  - section_probe not built, skipping the reach figures"
fi

if [ "$fails" -eq 0 ]; then
  echo "ok — $checks published figures still match the tool"
  exit 0
fi
echo "$fails of $checks published figures have drifted"
exit 1
