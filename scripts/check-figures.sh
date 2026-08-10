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
  # Validate the pointer on every call, deduplicated. See the note above `hint`.
  case "$checkedPointers" in
    *"|$quote|"*) ;;
    *) checkedPointers="$checkedPointers|$quote|"
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
  fi
else
  echo "  - zone_gpu_probe not built, skipping the §8 torn counts"
fi

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
else
  echo "  - section_probe not built, skipping the reach figures"
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
# ~10 s. `smoke_view` asserts these itself and exits non-zero if they move; what
# this adds is that the *document* still says what the tool says.
SMOKE=${SMOKE:-./build/smoke_view}
RENDER_DOC=docs/03-renderer-audio.md
if [ -x "$SMOKE" ]; then
  smoke=$("$SMOKE" --out=/tmp --frames=8 --duration=600 2>&1)
  if printf '%s\n' "$smoke" | grep -q 'no usable GPU'; then
    echo "  - no Vulkan device, skipping the fire and smoke figures"
  else
    for h in "at **834 K**" "peaks at **531 K**" "reaches 511" \
             "84.1% of the engine"; do
      hint "$h" "$RENDER_DOC"
    done
    # Each reader anchored to its own line, for the reason the damage figures are:
    # `531` also appears in the table two rows above the one being read.
    last_frame() { printf '%s\n' "$1" | grep 'smoke_07.png' | awk '{ print $'"$2"' }'; }
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
    check "optical depth across the room at t = 600 s" 511 3 "$(last_frame "$smoke" 7)" \
          "reaches 511" "$RENDER_DOC"
    check "the layer's extinction at t = 600 s (1/m)" 19.90 0.1 "$(last_frame "$smoke" 6)" \
          "| 19.90 | 511" "$RENDER_DOC"
    plan=$(printf '%s\n' "$smoke" | sed -n 's/.*\([0-9][0-9]\.[0-9][0-9]\)% of the bounding box.*/\1/p')
    check "the drawn prism, % of the bounding box in plan" 84.11 0.05 "$plan" \
          "84.1% of the engine" "$RENDER_DOC"
  fi
else
  echo "  - smoke_view not built, skipping the fire and smoke figures"
fi

if [ "$fails" -eq 0 ]; then
  echo "ok — $checks published figures still match the tool"
  exit 0
fi
echo "$fails of $checks published figures have drifted"
exit 1
