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
# `01`, `04` and `05` stay at zero on purpose: `01` and `04` are design documents
# whose digits are dates, section numbers and budgets nothing computes, and `05`'s
# are the fields of a file format that `test_shipfile.cpp` already parses. A figure
# is worth gating here when a *tool* produces it.
#
# Cost: the damage figures need no flooding at all and run with `--duration=1` in
# about a second each; three runs pay for the full 900 s. Around 105 s for the
# `ram_view` half, and another ~230 s for the README block below it -- the front
# page quotes three 1800 s scenario runs and the size of the test suite, and none of
# those can be had without running the thing. The renderer blocks add ~18 s: ~5 s of
# `seaway_view`, ~12 s of `smoke_view` and ~1 s of a sixth `ram_view` run.
set -u

RAM=${RAM:-./build/ram_view}
SECTION=${SECTION:-./build/section_probe}
SECTION_DOC=docs/02-simulation.md
RENDER_DOC=docs/03-renderer-audio.md
DOC=docs/06-roadmap.md
fails=0
checks=0
# **Every block that does not run is named in the summary.** See the note above the
# summary itself: a gate that skipped everything used to print `ok` and a count of
# zero, which `verify.sh` accepts, and the count is the only thing that would have
# said so.
skipped=""

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
# device check as everything else here. The deflections are a CPU closed-form
# comparison and would run without a device, but the tool has no way to ask for §1
# alone -- recorded as the reason rather than worked around.
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
if [ -x "$TESTS" ]; then
  suite=$("$TESTS" 2>&1 | sed -n 's/^\([0-9]*\) checks, [0-9]* failures$/\1/p' | tail -1)
  check "closed-form validation checks in the suite" 198869 0 "$suite" \
        "198869 validation checks" "$FRONT"
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
    check "triangles in the drawn scene: hull, interior and sea" 28019 0 "$scene" \
          "28 019 triangles" "$RENDER_DOC"
  fi
else
  echo "  - ram_view not built, skipping the collision-milestone figures"
  skipped="$skipped ram_view"
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
