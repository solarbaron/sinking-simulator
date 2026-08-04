# shipsim — working notes

A ship simulator where the ship is actually simulated: hull with material
properties, compartments that flood through real orifices at real rates, air that
compresses and finds its way out, buoyancy integrated over the instantaneous
wetted surface. `README.md` is the overview; `docs/` is the plan and the source
of truth.

## The one command

```sh
./scripts/verify.sh            # quick    build + tests            ~6 s
./scripts/verify.sh full       # + clean rebuild, GPU, scenarios   ~150 s
./scripts/verify.sh sanitize   # + ASan and TSan                   ~180 s
```

Run `quick` constantly and `sanitize` **before every commit** — on an idle
machine it is 180 s against `full`'s 150 s, so the extra coverage is nearly free
and there is no reason to reserve it for concurrency or raw-memory work. (Those
are idle figures: sharing the box with one busy benchmark took the same gate to
358 s. Budget accordingly before assuming it has hung.) **Warnings are
failures** — the build is
`-Wall -Wextra -Wpedantic` and has been warning-clean since Phase 0, so a single
warning means the signal is degrading.

An incremental build cannot see a warning in a file it did not recompile, which
is why `full` configures a throwaway build directory and compiles everything.

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
| Roll stiffness finite-differenced about the **body origin** while the moment was taken about the cog, so `zetaRoll = 0.08` delivered 0.144 | timing a free decay's log decrement against the ζ the label claims |

A green functional test is evidence the code does what you thought of, not that
it is correct.

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
  and the sweep shows no dispatch-limited regime. Grain matters ~40×; the queue
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
