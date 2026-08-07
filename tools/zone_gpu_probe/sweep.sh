#!/usr/bin/env bash
# **Every cell of one table, from one run of one tree, at one set of parameters.**
#
# `07-fem-spike-findings.md` §8 has published tables assembled cell by cell from
# separate invocations, and it has twice concluded that a cell "did not come from
# the run the rest of it did" -- once correctly (a torn count of 60 that no
# configuration reproduces, withdrawn) and once **wrongly**: five cells were
# annotated as not reproducing when the tool reproduces all five exactly. What had
# changed was not the tree and not the run but the *step count*, which the probe
# derives from the punch depth unless `--steps` says otherwise, and which differs
# by size because the critical timestep does -- 5 505, 5 513, 5 545 at 192, 768 and
# 3 072 elements. Eight extra steps move a torn count by one element and the
# dissipation by 0.07%.
#
# **So the failure mode is not a stale cell. It is a parameter of the experiment
# that was not written down next to the number**, and that is worse, because every
# control anyone thought to run varied something else: the re-run that produced the
# wrong correction proved its figures were not noise, not the reduction order, and
# not the tree -- it built the old commit from a `git archive` -- and all three were
# true and none could see it.
#
# This runs the whole sweep in one invocation and emits one CSV with a run id, a
# commit, a timestamp **and the step count** on every row. A table transcribed from
# one such file is a table from one run at one set of parameters by construction,
# and the file is the thing to re-run rather than the document.
#
#   tools/zone_gpu_probe/sweep.sh [--repeats=N] [--out=FILE] [--probe=PATH]
#                                 [--sizes=LIST] [--mappings=LIST] [--eas=LIST]
#                                 [--steps-scale=F] [--max-load=X] [--no-idle-check]
#
# **Timing figures need an idle box, and this checks rather than assumes.** The
# fp64 sweep §8 published was taken with an unrelated consumer holding the GPU at
# 100%: repeat passes of the same configuration differed by up to 35% and the CPU
# reference moved 47%, so nothing finer than one significant figure survived it --
# and the re-run showed that pass had *understated* the fp64 cost by 70% at one
# size, because contention inflates a short kernel more than a long one and so
# biases a ratio rather than merely blurring it. Three guards, because the
# resources fail differently and at different times:
#
#   - **Before each run**, while nothing of ours is running, the CPU's busy
#     fraction is sampled from `/proc/stat` over one second and the GPU's from the
#     driver, and the run waits for both to fall below `--max-load` per cent. A
#     build finishing next door is enough -- measured here, the first probe run
#     after a 19 s `ninja` read **3.5x** the settled CPU wall, and the settled
#     value repeats to 0.3%.
#
#     It is deliberately **not** the load average, which is the obvious instrument
#     and the wrong one: a 1-minute average still carries this sweep's own
#     previous run, so it would flag every row after the first and be waiting for
#     a quiet box that it is itself the reason is not quiet. A `/proc/stat` delta
#     taken in the gap between runs sees only other people's work, which is the
#     question being asked.
#   - **During each run**, `other_cpu_pct`: everything the machine burned minus
#     everything the probe burned. The before-check cannot see a neighbour that
#     starts two seconds later, and on this box neighbours do -- four rows of a
#     sixty-row sweep came back 15-83% slow on a before-check that read under 3%.
#   - **Repeats are interleaved**, not blocked: the whole sweep is run `--repeats`
#     times rather than each configuration being repeated in place. A box that
#     goes busy halfway then shows up as a spread on every configuration instead
#     of as a bias on the ones that happened to run during it.
#
# The spread across repeats is the useful output and the reason the CSV is per-run
# rather than pre-averaged: it is also the *sharpest* idleness check available,
# because it measures the thing that actually matters rather than a proxy for it.
# A GPU column that repeats to 2% could not have been taken against a competing
# load. **Reduce with the minimum, not the mean**: contention can only add time, so
# an inflated run is a run with someone else's work in it and averaging it in
# publishes a number that is partly the neighbour's.
set -u

PROBE=${PROBE:-./build/zone_gpu_probe}
OUT=
REPEATS=6
STEPS_SCALE=1
MAX_LOAD=
IDLE_CHECK=1
# elements:radius:sub:steps. The step counts are §8's and are not derived from the
# depth, because a run's cost has to be the same across the sweep for the columns
# to be comparable: the critical timestep falls with the element size, so a fixed
# punch depth would silently give the larger meshes more steps.
SIZES="192:2.5:4:5505 768:2.5:8:5505 3072:2.5:16:5505 8192:4.0:16:1500 16384:6.0:16:1000"
MAPPINGS="workgroup invocation"
EAS="float"

for arg in "$@"; do
  case "$arg" in
    --repeats=*)     REPEATS="${arg#*=}" ;;
    --out=*)         OUT="${arg#*=}" ;;
    --probe=*)       PROBE="${arg#*=}" ;;
    --sizes=*)       SIZES="$(echo "${arg#*=}" | tr ',' ' ')" ;;
    --mappings=*)    MAPPINGS="$(echo "${arg#*=}" | tr ',' ' ')" ;;
    --eas=*)         EAS="$(echo "${arg#*=}" | tr ',' ' ')" ;;
    --steps-scale=*) STEPS_SCALE="${arg#*=}" ;;
    --max-load=*)    MAX_LOAD="${arg#*=}" ;;
    --no-idle-check) IDLE_CHECK=0 ;;
    *) echo "unknown option $arg" >&2; exit 2 ;;
  esac
done

if [ ! -x "$PROBE" ]; then
  echo "sweep: no probe at $PROBE -- build zone_gpu_probe first" >&2
  exit 2
fi

cores=$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 8)
# Per cent busy, on both CPU and GPU, not a load average. 20% leaves room for a
# desktop -- a compositor, a browser and an editor idle at a few per cent here --
# and does not leave room for a build or a competing benchmark.
[ -n "$MAX_LOAD" ] || MAX_LOAD=20
run_id="$(date -u +%Y%m%dT%H%M%SZ)"
commit="$(git rev-parse --short HEAD 2>/dev/null || echo unknown)"
dirty=""
git diff --quiet 2>/dev/null || dirty="+dirty"
[ -n "$OUT" ] || OUT="sweep-${run_id}.csv"

# Whatever the driver will say about the device, once, at the top -- so a file
# that turns out to have been taken on a busy box says so in its own header
# rather than in someone's memory of the session.
gpu_line() {
  nvidia-smi --query-gpu=utilization.gpu,memory.used,clocks.sm,power.draw \
             --format=csv,noheader,nounits 2>/dev/null | head -1 | tr -d ' '
}
gpu_others() {
  # Who else has a context on the device. Read it as an inventory, not as a
  # verdict: on this driver an Electron app's GPU process appears here (VSCodium
  # does, at 56 MiB and 0.2% of a core) and it is not the hazard -- desktop
  # clients idle the device between frames. What the utilisation and clock
  # columns are for is telling an idle context from a busy one.
  nvidia-smi --query-compute-apps=pid,used_memory --format=csv,noheader 2>/dev/null \
    | grep -v '^$' | tr '\n' ';' | tr -d ' '
}

{
  echo "# zone_gpu_probe sweep ${run_id}"
  echo "# commit ${commit}${dirty}"
  echo "# probe ${PROBE}"
  echo "# host $(uname -srm), ${cores} logical cores, load at start $(cut -d' ' -f1-3 /proc/loadavg 2>/dev/null)"
  echo "# gpu at start (util%,mem,clock,power): $(gpu_line)"
  echo "# gpu compute clients at start: $(gpu_others)"
  echo "# repeats ${REPEATS}, interleaved; max load ${MAX_LOAD}"
  echo "run_id,commit,repeat,elements,radius,sub,steps,mapping,eas,seq,epoch,load_before,other_cpu_pct,gpu_util_before,gpu_clock_before,cpu_wall_s,gpu_wall_s,cpu_elem_ms,kernel_ms,cpu_torn,gpu_torn,cpu_diss_MJ,gpu_diss_MJ,diss_rel,cpu_peak_damage,gpu_peak_damage,status"
} > "$OUT"

# System-wide CPU busy fraction over one second, as a percentage. Called only in
# the gap between runs, so what it sees is other people's work.
cpu_busy() {
  awk '/^cpu /{ i1=$5; t1=0; for (f=2; f<=NF; f++) t1+=$f }
       END { print i1, t1 }' /proc/stat > /tmp/.sweep_stat1.$$
  sleep 1
  awk '/^cpu /{ i2=$5; t2=0; for (f=2; f<=NF; f++) t2+=$f }
       END { print i2, t2 }' /proc/stat > /tmp/.sweep_stat2.$$
  awk 'NR==1 { i1=$1; t1=$2 } NR==2 { i2=$1; t2=$2 }
       END { d = t2 - t1; if (d <= 0) { print "0.0"; exit }
             printf "%.1f", 100.0 * (1.0 - (i2 - i1) / d) }' \
      /tmp/.sweep_stat1.$$ /tmp/.sweep_stat2.$$
  rm -f /tmp/.sweep_stat1.$$ /tmp/.sweep_stat2.$$
}

# Wait for the box to go quiet on both resources. Returns the CPU busy fraction
# it settled at, or the one it gave up at after 120 s -- never blocks forever,
# because a permanently busy box is a result to record, not a reason to hang.
settle() {
  local waited=0 busy util
  while [ "$waited" -lt 120 ]; do
    busy=$(cpu_busy)
    util=$(nvidia-smi --query-gpu=utilization.gpu --format=csv,noheader,nounits 2>/dev/null | head -1 | tr -d ' ')
    if awk -v b="$busy" -v u="${util:-0}" -v m="$MAX_LOAD" \
           'BEGIN { exit !(b <= m && u <= m) }'; then
      echo "$busy"; return 0
    fi
    sleep 4; waited=$((waited + 5))
  done
  echo "$busy"
  return 1
}

# One probe invocation, parsed into one CSV row.
one() {
  local repeat="$1" elements="$2" radius="$3" sub="$4" steps="$5" mapping="$6" eas="$7" seq="$8"
  local load=0 log status prefix=
  if [ "$IDLE_CHECK" -eq 1 ]; then
    if ! load=$(settle); then
      printf '  ! %s%% of the box still busy (limit %s%%) after 120 s;'\
' running anyway and flagging the row\n' "$load" "$MAX_LOAD" >&2
      # Flagged rather than dropped, and flagged in the row itself: a sweep that
      # silently discarded contended runs would look like a clean sweep, which is
      # the failure this file exists to prevent.
      prefix="busy-"
    fi
  fi
  local before; before=$(gpu_line)
  local util clock
  util=$(echo "$before" | cut -d, -f1); clock=$(echo "$before" | cut -d, -f3)
  log=$(mktemp)
  local epoch; epoch=$(date +%s)
  # **How much of the box somebody else used while this row was being measured.**
  # The check above only sees the instant before the run; a sibling's test suite
  # that starts two seconds later lands inside the window and the row looks clean.
  # Measured as (all CPU time the machine burned) minus (the CPU time this probe
  # burned), over the run's own wall clock -- so our own 23 workers cancel and what
  # is left is other people's work. This found a sibling agent's `shipsim_tests`
  # inflating four CPU columns by 15-83% on rows whose before-check read under 3%.
  #
  # Read it as a change, not as an absolute: a desktop floor of about 10% of this
  # box is a browser, a compositor and a process monitor, and it is there whether
  # or not anything is competing for the solver's cores. What a competing build
  # looks like is this number in the fifties and above.
  # Busy jiffies, which is every column of the `cpu` line *except* idle (5) and
  # iowait (6). Summing the whole line instead makes `total` the wall clock times
  # the core count, so "everything the machine did" becomes "everything the machine
  # could have done" and an idle box reports 80% contention. It did, once.
  local busy1 busy2 wall1 wall2 childcpu
  busy_jiffies() {
    awk '/^cpu /{ t=0; for (f=2; f<=NF; f++) if (f != 5 && f != 6) t+=$f; print t; exit }' \
        /proc/stat
  }
  busy1=$(busy_jiffies)
  wall1=$(date +%s.%N)
  TIMEFORMAT='%3U %3S'
  # The status of a command substitution is the status of its last command, so
  # `$?` here is the probe's own exit code and not the shell's.
  childcpu=$( { time "$PROBE" --radius="$radius" --sub="$sub" --steps="$steps" \
                     --mapping="$mapping" --eas="$eas" > "$log" 2>&1; } 2>&1 )
  [ $? -eq 0 ] && status=ok || status=failed
  wall2=$(date +%s.%N)
  busy2=$(busy_jiffies)
  local other
  other=$(awk -v b1="$busy1" -v b2="$busy2" -v w1="$wall1" -v w2="$wall2" \
              -v c="$childcpu" -v n="$cores" -v hz="$(getconf CLK_TCK)" \
    'BEGIN { split(c, p, " "); mine = p[1] + p[2];
             total = (b2 - b1) / hz; wall = w2 - w1;
             if (wall <= 0 || n <= 0) { print "0.0"; exit }
             o = total - mine; if (o < 0) o = 0;
             printf "%.1f", 100.0 * o / (wall * n) }')
  grep -q '^skipped: ' "$log" && status=skipped
  # Guard against the row that looks like a measurement and is not: the probe's
  # own vacuity check refuses `ok` when the run neither deformed nor yielded the
  # patch, and a sweep that recorded those as timings would publish the cost of a
  # solver doing nothing.
  grep -q '^ok$' "$log" || [ "$status" = skipped ] || status=vacuous
  status="$prefix$status"

  local cpu_wall gpu_wall cpu_elem kernel cpu_torn gpu_torn cpu_diss gpu_diss diss_rel cpu_dam gpu_dam
  cpu_wall=$(awk '/^wall seconds/ { print $3; exit }' "$log")
  gpu_wall=$(awk '/^wall seconds/ { print $4; exit }' "$log")
  cpu_elem=$(awk '/of which kernel/ { print $5; exit }' "$log")
  kernel=$(awk  '/of which kernel/ { print $6; exit }' "$log")
  cpu_torn=$(awk '/^torn elements/ { print $3; exit }' "$log")
  gpu_torn=$(awk '/^torn elements/ { print $4; exit }' "$log")
  cpu_diss=$(awk '/^plastic dissipation/ { print $4; exit }' "$log")
  gpu_diss=$(awk '/^plastic dissipation/ { print $5; exit }' "$log")
  diss_rel=$(awk '/^plastic dissipation/ { print $6; exit }' "$log")
  cpu_dam=$(awk  '/^peak damage/ { print $5; exit }' "$log")
  gpu_dam=$(awk  '/^peak damage/ { print $6; exit }' "$log")

  printf '%s,%s,%d,%d,%s,%d,%d,%s,%s,%d,%d,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
    "$run_id" "$commit$dirty" "$repeat" "$elements" "$radius" "$sub" "$steps" \
    "$mapping" "$eas" "$seq" "$epoch" "$load" "${other:-}" "$util" "$clock" \
    "${cpu_wall:-}" "${gpu_wall:-}" "${cpu_elem:-}" "${kernel:-}" \
    "${cpu_torn:-}" "${gpu_torn:-}" "${cpu_diss:-}" "${gpu_diss:-}" "${diss_rel:-}" \
    "${cpu_dam:-}" "${gpu_dam:-}" "$status" >> "$OUT"
  printf '  r%d %6s %-10s %-8s cpu %8ss  gpu %8ss  kernel %10s ms  torn %4s/%-4s  else %4s%%  %s\n' \
    "$repeat" "$elements" "$mapping" "$eas" "${cpu_wall:--}" "${gpu_wall:--}" \
    "${kernel:--}" "${cpu_torn:--}" "${gpu_torn:--}" "${other:--}" "$status"
  rm -f "$log"
}

# A warm-up that is thrown away. The first run after anything else on the box is
# not a measurement -- see the 3.5x above -- and one discarded run is cheaper than
# an outlier that has to be explained afterwards.
printf 'warm-up (discarded)\n'
"$PROBE" --radius=2.5 --sub=4 --steps=500 >/dev/null 2>&1 || true

seq=0
r=1
while [ "$r" -le "$REPEATS" ]; do
  printf 'repeat %d of %d\n' "$r" "$REPEATS"
  for size in $SIZES; do
    elements=$(echo "$size" | cut -d: -f1)
    radius=$(echo "$size" | cut -d: -f2)
    sub=$(echo "$size" | cut -d: -f3)
    steps=$(echo "$size" | cut -d: -f4)
    steps=$(awk -v s="$steps" -v f="$STEPS_SCALE" 'BEGIN { printf "%d", s * f }')
    for mapping in $MAPPINGS; do
      for eas in $EAS; do
        seq=$((seq + 1))
        one "$r" "$elements" "$radius" "$sub" "$steps" "$mapping" "$eas" "$seq"
      done
    done
  done
  r=$((r + 1))
done

{
  echo "# gpu at end (util%,mem,clock,power): $(gpu_line)"
  echo "# gpu compute clients at end: $(gpu_others)"
  echo "# load at end $(cut -d' ' -f1-3 /proc/loadavg 2>/dev/null)"
} >> "$OUT"

printf '\nok -- %d rows in %s\n' "$seq" "$OUT"
