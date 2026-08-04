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
# Kept cheap deliberately: the damage figures need no flooding at all, so they run
# with `--duration=1` in about a second each, and only the one headline outcome
# pays for the full 900 s of flooding. About 40 s in total against the ~400 s the
# `full` gate already costs.
set -u

RAM=${RAM:-./build/ram_view}
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
  local label="$1" expect="$2" tol="$3" actual="$4" quote="$5"
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
         "$red" "$off" "$label" "$DOC" "$expect" "$actual"
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
  checks=$((checks + 1))
  if grep -qF -- "$1" "$DOC"; then return 0; fi
  printf '  %s✗%s the pointer %s"%s"%s no longer occurs in %s\n' \
         "$red" "$off" "$dim" "$1" "$off" "$DOC"
  fails=$((fails + 1))
  return 1
}

if [ ! -x "$RAM" ]; then
  echo "  - ram_view not built, skipping published-figure check"
  exit 0
fi

for h in "from 3.4 to" "63 bays torn, 107.7" "amidships she takes" \
         "at the quarter she takes" "GM −1.62 m"; do
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
check "hole at 1.5 m/s (m2)" 3.4 0.3 "$area15" "from 3.4 to"
check "hole at 6 m/s (m2)" 107.7 0.6 "$area6" "from 3.4 to"
check "bays torn at 6 m/s" 63 1 "$torn" "63 bays torn, 107.7"

# --- the headline outcome, which needs the full 900 s ---------------------------
out=$("$RAM" --speed=6 2>&1)
water=$(outcome_line "$out" | sed -n 's/.*-- \([0-9]*\) t of water.*/\1/p')
heel=$(outcome_line "$out" | sed -n 's/.*heel \(-\{0,1\}[0-9.]*\) deg.*/\1/p')
gm=$(outcome_line "$out" | sed -n 's/.*GM \(-\{0,1\}[0-9.]*\) m.*/\1/p')
check "floodwater amidships (t)" 2210 45 "$water" "amidships she takes"
check "heel amidships (deg)"     8.1  0.5 "$heel" "amidships she takes"
check "GM at loss (m)"          -1.62 0.05 "$gm"  "GM −1.62 m"

# --- and that being struck at the quarter is worse, and the other way -----------
out=$("$RAM" --speed=6 --aim=-30 2>&1)
qwater=$(outcome_line "$out" | sed -n 's/.*-- \([0-9]*\) t of water.*/\1/p')
qheel=$(outcome_line "$out" | sed -n 's/.*heel \(-\{0,1\}[0-9.]*\) deg.*/\1/p')
check "floodwater at the quarter (t)" 7102 150 "$qwater" "at the quarter she takes"
check "heel at the quarter (deg)"    -45.0 1.0 "$qheel" "at the quarter she takes"

# The claim that survives all of it, and the only one that is a *finding* rather
# than a figure: a 32x range of hole size barely moves the floodwater, because a
# small breach fills the compartment behind it inside 900 s just as a large one
# does. Asserted as a spread rather than as values, so it keeps meaning something
# when the ship is re-tuned.
lo=$(outcome_line "$("$RAM" --speed=1.5 2>&1)" | sed -n 's/.*-- \([0-9]*\) t of water.*/\1/p')
if [ -n "$lo" ] && [ -n "$water" ]; then
  checks=$((checks + 1))
  if awk -v a="$lo" -v b="$water" \
        'BEGIN { r = a / b; if (r < 1) r = 1 / r; exit !(r < 1.15) }'; then
    printf '  %s✓%s a 32x range of hole size moves the floodwater under 15%% (%s t vs %s t)\n' \
           "$green" "$off" "$lo" "$water"
  else
    printf '  %s✗%s the threshold finding no longer holds: %s t at 1.5 m/s against %s t at 6 m/s\n' \
           "$red" "$off" "$lo" "$water"
    printf '      %sthis is a finding, not a figure -- %s explains it under "Beyond a threshold"%s\n' \
           "$dim" "$DOC" "$off"
    fails=$((fails + 1))
  fi
fi

if [ "$fails" -eq 0 ]; then
  echo "ok — $checks published figures still match the tool"
  exit 0
fi
echo "$fails of $checks published figures have drifted"
exit 1
