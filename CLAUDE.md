# shipsim — working notes

A ship simulator where the ship is actually simulated: hull with material
properties, compartments that flood through real orifices at real rates, air that
compresses and finds its way out, buoyancy integrated over the instantaneous
wetted surface. `README.md` is the overview; `docs/` is the plan and the source
of truth.

## The one command

```sh
./scripts/verify.sh selftest   # the gate's own negative controls      ~2 s
./scripts/verify.sh            # quick    build + tests              ~120 s
./scripts/verify.sh full       # + clean rebuild, GPU, scenarios     ~2300 s
./scripts/verify.sh sanitize   # + ASan and TSan                     ~3200 s
```

**Run `selftest` when you touch the gate.** It drives each reporter with a failure
it was *not* written for and asserts it goes red — 21 controls in two seconds, at
the head of every level. It exists because an audit found **fourteen checks in
this gate that could not go red at all**, including three that accepted
`0 checks, 0 failures` and two that called a rebuild compiling *nothing*
"warning-clean".

The figures above were re-measured after that audit; the ones they replace were
stale by roughly five times, having been written when the suite was a fraction of
its present size. Re-measure rather than trusting them, and budget generously on
a shared box — the same gate took 180 s and 358 s at an earlier size depending
only on what else was running.

Run `quick` constantly and `sanitize` **before every commit** — the extra
coverage over `full` is about a third more wall time, which is cheap enough that
there is no reason to reserve it for concurrency or raw-memory work.

**Warnings are failures** — the build is `-Wall -Wextra -Wpedantic` and has been
warning-clean since Phase 0, so a single warning means the signal is degrading.

An incremental build cannot see a warning in a file it did not recompile, which
is why `full` configures a throwaway build directory and compiles everything.

**And a build cannot see a warning its optimisation level does not produce.** Every
gate compiles `RelWithDebInfo`; GCC's `-Waggressive-loop-optimizations` needs `-O3`,
and it sat on `reduction.cpp`'s triangular solves for as long as they existed —
a signed loop counter it could not bound, so it reasoned `n` might be `INT_MAX` and
the indexing would then be undefined. `full` now builds the engine at `-O3` as well,
which costs 6 s. Tools and tests are still not covered there.

## Conventions

- **SI units everywhere**, no exceptions. Angles in radians except at display
  boundaries.
- **Body frame**: +x forward (bow), +y to port, +z up, origin at midship on the
  baseline. Heel positive starboard-down, trim positive bow-down.
- **Vulkan clip space** for the renderer: y points *down*, z in [0, 1]. Both
  differ from OpenGL and both produce a plausible-looking wrong image if muddled.
- Components in the ECS are plain data, trivially copyable, never destructed.
- Match the surrounding code's comment density and idiom. Comments explain *why*,
  especially when the obvious implementation is wrong.

## How testing works here

Assert against **closed-form or independently-derived answers**, never against
eyeballed output. Where a rendered or simulated result has an analytic value —
a triangle covering half a viewport, `Hs = 4√m₀`, a log decrement of
`2πζ/√(1−ζ²)` — assert the analytic value.

Aim tests at **failure modes that stay silent rather than crash**. Prefer an
assertion that would fail against a plausible-but-wrong implementation:

- Draw the *far* surface second, so the depth test has to actually work.
- Destroy entities from the *front*, so swap-remove has to patch the moved record.
- Feed the loader *every truncation* of a valid file.
- Sweep alignments 1–256, not the three anyone writes by hand.
- Use `cos` not `sin` when a symmetric range would make the wave integrate to zero.

**Guard against vacuous tests.** Several tests here carry an explicit check that
the thing being measured is non-trivial — that fold order changes the result,
that the analytic answer differs from the still-water value, that the mirror
comparison had something to compare. Each of those guards exists because the
first version of that test passed while proving nothing.

**When a test fails, decide whether the test's expectation or the design is
wrong before loosening any assertion.** That habit has caught, in this repo: a
floating-point accumulator that needed an integer clock; a load path that failed
open; a test that mis-computed a projected area; a camera pointed down the ship's
axis; and a wave whose contribution cancelled by symmetry.

## What has actually gone wrong here

Every subsystem shipped green on its functional tests while still containing a
real defect. In each case a *different instrument* found it. This is the single
most useful thing to know about this codebase:

| Defect | Found by |
|---|---|
| Job records recycled while still queued | running with **zero workers** |
| Chase-Lev slot-reuse data race | ThreadSanitizer |
| Auto-grain clamp off by one | TSan as a *slow operating point*, not as a race detector |
| Arena aligning offsets instead of addresses | sweeping alignments 1–256 |
| Arena overruns invisible to ASan | manual poisoning + a deliberate negative control |
| Hull wound inconsistently, displacement 40% high | a manifold check, added while doing something else |
| Scheduler drifting over long runs | a 999-vs-1000 count a tolerance would have hidden |
| `World::load` leaving a half-built world | every truncation of a valid save |
| Under-tessellated hull inventing ±6% displacement | a short-wave test that should have cancelled |
| Sea surface queried 6× more than necessary | timing the real tick instead of extrapolating |
| Radiation solver returning **negative** damping at four frequencies in the seakeeping band | a near-field/far-field **energy balance**, not any coefficient test |
| Retardation fitted over 1.3 s of a 20 s decay, turning a damper into an integrator | a free-decay test, after the ship reached NaN in five steps |
| Radiation damping added *on top of* the modal damping standing in for it | comparing an RAO sweep before and after, not any single run |
| Reading `state.velocity.x` as "speed" — it is a *world* vector | a steady turn that looked like chaos until surge was taken along the bow |
| `makeHullFromStations` padding a short station with **zeros** | writing a *second* way to build the same ship and comparing them |
| A wing tank authored **inside** a hold, 217 m³ flooding twice | a pairwise overlap check; the total-volume one saw 89% and was happy |
| `Mat3{}` is the **identity**, so an accumulator started at I | a second moment measured two ways: every volume right, every diagonal one too large |
| Two coincident faces counted twice, overlap 5/3 too large | an overlap coming out *larger than one of the solids* |
| Roll stiffness finite-differenced about the **body origin** while the moment was taken about the cog, so `zetaRoll = 0.08` delivered 0.144 | timing a free decay's log decrement against the ζ the label claims |
| A published damage table whose energy column was ½·v²·the **struck** ship's mass — the wrong ship, and as if all of it reached the plating | re-running the tool that generated the table instead of quoting the table |
| Three documents claiming a fix was still unapplied, above code that applied it | the same re-run; nothing tests a comment |
| "The hole is 10× too big" inferred from a **force** ratio, when the hole is set by energy per unit area and moves 5% | measuring the two effects separately instead of propagating one number |
| The ferry's **mid wing tanks were never authored**, so 41% of a ram amidships tore open onto no compartment and could not flood | reading a warning `ram_view` had printed on every run for months |
| A published "hole size stops mattering" finding that was an **artefact of that gap** — the extra hole opened onto nothing by construction | re-deriving the finding after the defect was fixed, not just the figures |
| A rigid-body count that found **three of six** modes: translations come out at exactly zero, rotations at 8e-4…9e-3, and a fixed 1e-3 rad/s cutoff landed between them | scaling the threshold off an independently-derived frequency instead of a constant |
| Two functions on the caller's own path shipped **with no test at all**, in a commit whose headline feature was well tested | asking what in the diff was *not* exercised, rather than whether the tests passed |
| A section cut **admitted by a tolerant membership test and then dropped by an exact geometry one**, losing all 188 plate panels on 11 of 51 stations while all 181 stiffeners survived | a sweep along the length with the two element populations counted *apart* |
| The bow and stern refusing to mesh — diagnosed for months as *inverted* elements, when **not one Jacobian on the ship is negative**: 166 input panels are triangles wearing four corners | sampling the nine points the element is *integrated* at, rather than the eight corners |
| A node ordering chosen on the element graph while the solver assembled a *constrained* one, so a reported half-bandwidth of 146 was really 10 769 and the suite stopped finishing | a test that **rebuilds the assembled band from the mesh and the constraints** instead of asserting the number the mesher reports |
| `reduction::recover` leaving every MPC-eliminated DOF at **zero** — a hole in the displacement field that read as 850 788 MPa against a true 427, on every tied section the mesher builds | adding a *second* reader of the same recovered field, and asking what else already read it |
| A compartment filled *after* `Ship::initialise` keeping the air mass it had while empty: 184 kPa of trapped gas, **8.3 m of head that is not water**, and a bulkhead that failed in the first step of every run | running the *controls*, which failed too — a scenario that fails is a result, a control that fails is a bug |
| A fire with no decay phase felling that bulkhead on its own at 2760 s, because steel asymptotes to the gas temperature and the restrained-buckling limit is a fixed one — so "the fire alone does not do it" was a statement about how long anyone watched | running the same control **three times longer** |
| One convective film per gas layer misbooking **8%** of the boundary exchange: the radiative coefficient goes as `T_s²`, and the foot of a bulkhead is held near ambient by the water behind it while its head is at 500 K | publishing the linearisation error next to the exact integral instead of assuming a mean was good enough |
| A restraint window that came back equal to whatever restraint the run was given, on every run — because the failure floods the compartment and **relieves the head that caused it** | asking the same question at two different inputs and getting the input back both times |
| The gate's own build step failing and printing **nothing**: it greps the log for `error:`, and a compiler killed by a full disk never says that word | a negative control that fails *without* the word the reporter looks for |
| 98 gated figures and **not one of them on `README.md`** — the front page, where the verdict a reader sees first had just changed | counting the gate's own coverage by document, rather than by figure |
| A **mutation left applied** in the source by a sweep that was killed mid-iteration — and it did not fail loudly, it made the suite *hang* | a script that re-derives every injected substitution and greps the source for each, run *before* trusting a green tree |
| The gate reporting `built warning-clean` off a rebuild that **compiled nothing** | a `ninja` stub that says `no work to do.` — the fix asserts an edge count |
| `0 checks, 0 failures` satisfying a grep for `checks, 0 failures` in **three** places — unit, ASan and TSan — in the same week as `ok — 0 published figures` | driving each reporter with a *degenerate* run rather than a failing one |
| A scenario that printed its verdict and then **segfaulted**, and a TSan run that crashed after its summary — both green, because the exit status was never read | asking what the status was, having already parsed the output |
| A tool reporting `no usable GPU` **with a card present** passing as a skip, on four tools | classifying from outside the tool: an ICD and a device node exist, so a skip is a failure |
| A boundary relaxation using `c_p` where `c_v` belongs, invisible to every existing test because they all stepped short enough that `C(1−e^{−rΔt/C}) → rΔt` **whatever `C` is** | asking the question at a step long enough for the exponential to bend |
| A promotion criterion comparing `length(state.velocity)` — a **speed** — against a threshold in m/s², and taking world `angularVelocity.x` as roll on a ship at 58° of loll. Every test agreed, because the fixture set a velocity and *called* it an acceleration | running the criterion on the real ferry next to a beam-sea control, rather than on the synthetic ship the tests drive it with |
| Two budgets for one cost model, set independently: 2000 tiles admitted 16 m³ where 100 000 particles admitted 100, so **no compartment over 16 m³ could ever be promoted** and the tier could not reach the vehicle deck it exists for. The 154 refusals were silent — nothing read `problems` — and read exactly like hysteresis working | seeding the compartment the tier is *for* and asking why it never promoted, then printing the refusal the review had been returning all along |
| A cost estimate of `5.0 core-s/sim-s per compartment` that was low by **5.6× to 606×**, and was the wrong *shape*: a per-compartment constant where the real cost scales with the water. It had carried `// estimate, will be measured` for as long as it existed | running it — the comment naming the measurement as outstanding is not the measurement, and the whole tier's affordability rested on it |
| A correction that was itself the error: "tiles are **92%** of the footprint" computed at the estimator's 1000 particles/m³, in the same comment block that had just established byte claims must use the solver's real 64 000/m³. At the real density tiles are **10%** and the line it "corrected" was right | an audit asking which claims in the subsystem nothing measures — the contradiction was two paragraphs apart in one entry, written in one sitting |
| `flip::Particle` documented as "~80 bytes" in a comment that **itemises the sixteen doubles adding to 128** — every megabyte figure in the tier was built on it | `sizeof`, once someone asked. The parenthesis had contradicted itself since it was written |
| `coreSecondsPerElement = 4.0` in `promotion.hpp`, citing as its source the `zone.hpp` paragraph that names 4.0 as the pre-`cacheRestForms` figure and says "every figure below that predates it is **2.4x pessimistic**" | asking whether the tier's *two* cost estimators agree. They stood 2.35× apart, `zone_probe` printed both on the same run, and nobody had read them side by side |
| Ikeda's `B0` carrying `m2^2` where every other term in it and in `A0` is degree 2 — a girth times a lever, both over `d` — so the one term setting the bilge keels' hull-pressure damping did not scale like the rest | sweeping `OG/d`, the *only* input the two readings differ in. The check that was there swept frequency and amplitude, nine points, and held the hull at the one `OG/d` where the wrong reading sits within 0.05 of the right one: 3.6% agreement there, **15.6%** at the edge of the same validated range |

**Two estimators for one quantity will disagree, and printing both is not
comparing them.** Three separate instances here: the water tier's particle and
tile budgets (6.25× apart, making the tier unreachable), the structural tier's
`coreSecondsPerElement` against `zone::estimatedCost` (2.35×), and
`estimateFlipCost`'s truncation against its own exact form (0.08%, from
`10.0/(64·0.05³)` landing on 1249.9999999999998). Each was visible in output
nobody read as a pair. **When a second way to compute something exists, assert
they agree** — and drive the real one rather than re-deriving its arithmetic, because
a test that recomputes the model by hand agrees with a broken model and disagrees
with nothing. The water tier's budget test claimed in its own comment to assert
"against the arithmetic rather than the constants" while doing exactly that.

A green functional test is evidence the code does what you thought of, not that
it is correct. **Nothing at all tests a comment** — three documents here repeated a
wrong factor of ten because each quoted the previous one rather than the
measurement, and the code beneath them had been right for weeks.

**A loose assertion is nearly a vacuous one.** An interface displacement that is
exact to 8e-12 was first asserted at 1e-3, which would have passed on a model that
had lost the property entirely and was merely well converged. When a measurement
comes back far better than the tolerance, tighten the tolerance to what was
measured and say why it is allowed to be that tight.

**A tolerant test followed by an exact one is a trapdoor.** `sectionElements` admitted
a panel with a `1e-9` membership tolerance and then looked for a sign change with
exact arithmetic; a query one ULP outside the bay passed the first and failed the
second, and the panel vanished without a word. The stiffener branch survived the
same query only because it happened to `clamp` its parameter. When two tests in
sequence decide one thing, they need the same notion of "on the boundary".

**Writing a lesson down does not make it learned.** The rigid-body threshold above
is the second time that mistake was made in this repo — the first was recorded, in
this table, by the work that immediately preceded it.

**A grid that sweeps only the axes a defect is constant along is a one-point
grid.** `testBilgeKeelAgreesWithTheSimplifiedRegression` compared Ikeda's sectional
bilge-keel model against Kawahara's independent regression at nine points — three
frequencies by three amplitudes — and passed at 3.6% against a 6% band for as long
as it existed. The term it was pointed at is a function of `OG/d` alone. Neither
swept axis appears in it, so the nine points were nine copies of one point, and the
single hull they were taken on sits within 0.05 of where the wrong reading and the
right one cross. Swept over the *validated range of that one parameter* and nothing
else, the same defect reads 15.6%. **Before trusting a comparison, ask which input
the two things being compared could disagree about, and check that the grid moves
it** — a wide sweep of the wrong parameter buys nothing, and it buys it while
looking thorough.

**A figure a program prints is checked by everyone who looks; a figure a person
copied is checked by nobody.** Seven published tables were re-derived from their
own tools in one session. The four that had drifted — the roll damping-ratio table,
the roadmap's GM-detail rows, the hull-form convergence table, and the barge RAO
sweep — were every one of them transcribed once by hand from a run nobody repeated,
and two were badly out: a `none` row off by a factor of 34, and nine convergence
percentages two to five times too large. The three that had not moved a digit —
strip-theory radiation, propulsion and manoeuvring, hull-to-hull contact — were
every one of them printed by a test on every run.

The distinction is not the subsystem, the age of the code, or how load-bearing the
claim is. It is whether anyone sees the number between the day it was written and
the day it is questioned. **So the repair for a stale figure is not to correct it;
it is to make a tool print it and then gate it** — correcting it alone puts the
next reader in exactly the position of the last one. Four of the fixes above
therefore added a `printf` before adding a `check`, and two of those found the
document wrong the moment the number was printed for the first time.

**And nothing gated any of it.** The two published figures that comparison was
protecting — `B44hat = 0.0439` and `6.3% of critical`, each written down in three
places — are not among the 346 the figure gate re-derives, because Ikeda is only
ever reached through `shipsim --bilge-keels=`, a flag `check-figures.sh` never
passes. A whole validated subsystem, fifteen published numbers, and not one run in
the gate that constructs a `RollDampingHull` at all. This is the README hole again
in a different place: coverage counted as a total rises without arriving anywhere
in particular, and **the way to see it is to ask which flags the gate never passes**,
not how many checks it has.

**"Restore *and* rebuild" is not enough if the restore moves the mtime
backwards.** The harness that produced the kills above did call rebuild after each
mutant, and the tree still came back **not green**. The restore was
`shutil.move(f + '.bak', f)` — a rename, which carries the backup's *own* mtime,
taken before the mutant was compiled. `ninja` compared that against an object file
built from the mutant, called it up to date, and did nothing; the source read clean
and the binary kept the mutant. The damage was not to the mutant that caused it but
to the *next* one in a different file, which then ran with two mutations live and
scored a kill that was partly somebody else's. `touch` the file after restoring, or
copy the contents back instead of renaming, and assert the tree is green after each
restore rather than only at the end. This is the same trapdoor the paragraph above
describes, reached through the clock instead of through a missing call — and it was
written by someone who had just read that paragraph.

**`cmd | tail -20` reports the exit status of `tail`.** The sanitize gate above was
launched as `./scripts/verify.sh sanitize 2>&1 | tail -20` and run in the
background. It came back **exit 0**, and its last line was
`✗ 1 failure(s) in 4322s`. A pipeline's status is its *last* stage, so the harness
recorded a passing run of `tail` and the gate's own verdict reached nobody — and
because `tail` keeps only the end, every line identifying *which* check failed was
discarded at the same time. Two losses from one pipe: the status, and the evidence.
`verify.sh`'s own step logs are `mktemp`'d and cleaned up, so there was nothing to
go back to. Redirect to a file and read the file; if a pipeline is unavoidable, set
`pipefail`. This is the fourth shape of the same defect in one session — the
reporter that only sees expected failures, the `grep FAIL` that scored a segfault
as a survivor, the waiter below that matched itself, and now a status that was
never the gate's.

**A waiter whose predicate matches its own command line never fires.** Seven
background waiters had been armed across a session as
`until ! pgrep -f "verify.sh sanitize"; do sleep 20; done`. `pgrep -f` matches the
whole command line, and each waiter's own command line contains that literal
string, so every one of them matched itself, the condition was permanently false,
and none could ever fire — one had been spinning for over a day. The tell was
`pgrep -c -f 'verify.sh sanitize'` returning **8** with exactly one gate running.
Wait on the pid, which cannot self-match: `while kill -0 "$pid"; do sleep 30; done`.
This is the same shape as the reporters above — a detector that includes itself in
what it detects — and it survived seven arming attempts because its failure mode is
*silence*, which is exactly what a waiter looks like when the thing it waits for has
not happened yet.

**A gate reads the tree and the build directory for its whole duration, so
editing either while one runs produces a failure that is yours and looks like the
code's.** Three times in one session: twice `verify.sh` caught a check count that
had moved since the run began, and once a figure check reported three
section-mesher figures as drifted when `ninja` had relinked `section_probe`
underneath it — the run started at 02:16:06 and the binary was rebuilt at
02:17:45, with the section block at line 93 of about a hundred. That third one
cost the most, because the figures were *correctly* reporting an experimental
change to the strake fractions and the obvious reading was a defect in a mesher
that had not been touched for sixty-five commits. All three reproduced their
published values the moment the tree was clean. Start the gate, then stop typing.

**A mutation run that is interrupted leaves a tree that is green because the
mutant survived.** This is the one place where a passing suite is actively
evidence *against* you, and it has happened here: an agent killed mid-sweep by a
session limit came back to a clean-looking worktree with substitution 43 still in
`les.cpp`. Never resume a sweep, and never commit from one, without first
re-deriving every substitution you injected and grepping the source for each.
Keep that scanner next to the harness rather than in your head — the whole
premise of mutation testing is that you cannot tell a surviving mutant from
correct code by looking at test output, and that applies to *your own* leftovers.
The same run is also why the harness must be able to detect a **hang**: mutant 43
did not fail, it drove a rejection test into a corner and stopped.

**A mutation harness that only distinguishes pass from fail is under-powered here,
because the characteristic kill in this codebase is a hang.** That is now measured
rather than anecdotal: one leftover mutant in `les.cpp`, and *five of 196* in
`fire.cpp`/`thermal.cpp` — a zone swap, an interface height taken from the floor, a
band donor swap, a direction flip and a sign on the wall relaxation — killed the
suite with **zero failing assertions** by turning a nine-second run into an
hours-long one. Whether those score as kills or survivors depends entirely on where
the timeout lands, so a sweep without a per-mutant time bound is reporting a number
it did not measure. Assert against the arithmetic floor of the substep controller,
not against a wall clock, so the bound survives a busy box.

**A ratio is only safe when both sides are the same kind of measurement.** The
note below says what survives instrumentation is a ratio between two timings,
because it scales on both sides. That is true and it is not sufficient: the figure
gate compared the auto-tuner's *single* timed run against the **minimum of eight**
swept timings, and a minimum over eight noisy samples is biased low, so the ratio
sits above 1 even when the tuner is perfect and the bias grows with the noise. The
bound was 1.15. Fifteen runs on an idle box spread 0.969–1.125, four of them above
1.07; inside a gate run it reached **1.282** and went red on code nobody had
touched. Both sides must be the same *statistic*, not merely the same units — and
where they cannot be, set the bound from the loaded distribution and record the
distance to the nearest genuine failure. Here that distance was there to use: the
tuner is deterministic across twenty runs, its chunk counts are gated exactly, and
the cheapest wrong landing costs 1.50× against noise reaching 1.28.

**A timed quantity cannot be bracketed against a constant, and a minimum over
attempts does not fix it.** The gas tier's `coreSecondsPerCell` was documented as
"measured in the tests" while the test only asserted it was positive, so a bracket
was added: the constant against the coarsest and finest per-cell timings, taking a
minimum over three attempts on the argument that contention can only slow a run.
It holds on an ordinary build and fails under **both** sanitizers. Instrumentation
multiplies a timed solve by about ten and leaves the constant exactly where it is,
so both bounds rise past it together -- a minimum over attempts defends against a
busy neighbour and does nothing whatever against a slower binary. This is the same
false-kill hazard as above, reached from the opposite direction: not a wall clock
asserted as a bound, but a constant asserted against a wall clock. **What survives
instrumentation is a ratio between two timings, because it scales on both sides.**

**And a mutation harness on a busy box reports strength it does not have.** A
wall-clock assertion falsely killed two *controls* during the GM sweep — controls
being mutants that are supposed to survive. A false kill inflates the rate and
hides a real gap, which is the one direction of error a kill rate cannot afford.

**A failure path that only reports the failures it already anticipated is not a
failure path.** `verify.sh` had this twice in adjacent functions: `expect_ok` and
`build_into` each summarised a failing step by grepping its log for the text a
*foreseen* failure would contain, so an unforeseen one — a signal, a full disk —
came out as silence. Both now print the exit status, separate a signal from an
assertion, and show the tail regardless of whether anything matched. **Test a
reporter with a negative control that fails in a way it was not written for**; the
control is also what caught `$?` being the status of the *negation* after
`if ! cmd`, which was written into both functions, hours apart, by someone who had
already fixed it once.

**A mutation harness that restores the source and not the binary is a fourth
door into the same room.** The sweep script above copied the file back and stopped
there, so the next `./build/shipsim_tests` -- run to confirm the tree was clean --
executed the *mutant's* binary against restored source and reported one failure.
It reads exactly like a real regression, and the instinct it trains is to go
looking for a bug in the code that was just written. Restore and **rebuild**, then
check; a clean source tree is not a clean tree.

**This has now happened a third time, in the ad-hoc detector rather than in the
code.** A mutation sweep over `edgeDrive` piped each run through `grep "FAIL"` and
reported one mutant as surviving. It had not survived: removing that guard is a
heap overread, and the binary died at **exit 139** before printing a summary, so
there was no `FAIL` line to match and the loudest possible kill scored as a
escape. The detector recognised only the failure shape that was expected —
exactly the defect the paragraph above describes, written by the same hands that
had just written the paragraph. **Judge a mutation run on the exit code and the
summary line, never on the presence of a string**; a mutant that crashes, hangs or
dies before it reports is the normal case here, not the exotic one. The same run
also shows why it matters that a check be pointed at the right claim: `grep -c`
returning 0 was, on a previous occasion, proof a mutation *had not applied* — the
identical observation meaning the opposite thing.

## Settled decisions — do not reopen

Recorded in `docs/01-architecture.md` §1 with full reasoning:

- **Task helping, not fibers.** Fibers are *rejected*, not deferred — `thread_local`
  stops being stable across a yield, which is a broad class of subtle bugs for a
  problem this workload does not have. The real answer is continuation-style
  scheduling so nothing blocks.
- **MPMC ring, not Chase-Lev,** as the job queue. Chase-Lev is only safe when a
  slot is an atomic that can be read speculatively; an 80-byte job record is not.
  Recoverable later via 32-bit handles into the frame arena — *not* via epoch
  reclamation.
- **Chase-Lev revisit cancelled on evidence**: dispatch is 0.2% of a 10 µs chunk,
  and the sweep shows no dispatch-limited regime. Grain matters at least 17× and
  worse with more workers (21× at 8, 41× at 23); the queue
  does not.

The renderer targets a **GTX 1070 Ti (Pascal)**: no mesh shaders, no hardware ray
tracing.

## Docs are the source of truth

`docs/` must not drift from the code. When a measurement or an obstacle changes a
decision, record *what was learned*, not just what was done — several docs carry
corrections where a claim was measured and turned out wrong (convergence order,
job-size targets, the Ikeda nondimensionalisation).

Prefer measuring over asserting: if a cheap experiment settles a design question,
run the experiment. **An extrapolated figure can be numerically right and still
point at the wrong fix** — the wave-cost estimate predicted 23 ms/tick correctly
and prescribed a vectorised sincos, while the actual first fix was a 6× query
redundancy that no per-component figure could reveal. Extrapolate to decide
whether to measure, not to decide what to do.

## Notes

- GPU work must **skip, not fail**, when there is no Vulkan device, and stays out
  of the TSan build because the driver is uninstrumented.
- If a sanitizer reports something, read the stacks and establish whether it is
  *ours* before suppressing. Suppress by library name, never by disabling the
  check. The libdbus leaks come from the driver stack, not from us.
- Commit with a message file (`git commit -F`), not an inline `-m` — embedded
  quotes have broken the shell here before.
- After a scripted/`sed` edit, grep to confirm it landed. A non-matching
  search-and-replace once silently dropped a test assertion, and the resulting
  failure looked like a bug in working code.
- `git add -A` once swept agent worktrees into a commit. `.gitignore` covers it
  now, and `scripts/install-git-hooks.sh` installs a pre-commit check.
