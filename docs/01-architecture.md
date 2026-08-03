# 01 — Engine Architecture

The engine exists to serve one unusual constraint: **the simulation is the
expensive part, not the rendering.** A conventional game engine budgets ~14 ms to
the GPU and a couple of milliseconds to gameplay. Here it is closer to inverted —
a damaged 200 m hull with an active FEM zone, a compartment flooding network, a
fire in a machinery space and a wave field is tens of milliseconds of CPU work per
frame if written naively. Every architectural decision below follows from that.

---

## 1. Threading model

**Work-stealing job system, no per-system threads.** Implemented in
`engine/core/jobs.cpp`.

- One worker thread per hardware thread minus one. On the 24-thread development
  machine, 23 workers. Zero workers is a supported configuration and every
  concurrency test runs it, because it is the debugging and determinism baseline.
- Jobs are small and expressed as parallel-for over ranges, not as long-lived
  tasks. The original target here was 50–500 µs; **measurement says that is about
  25× more conservative than it needs to be** — see "What the benchmark said"
  below. Chunks reach full efficiency from roughly 2 µs.
- The frame is a DAG of jobs built fresh each tick. No system ever calls another
  system directly; they declare read/write sets over component types and the
  scheduler derives ordering.
- Deterministic reduction order everywhere results are summed across threads.
  Floating-point addition is not associative and multiplayer will not forgive it.

### Two decisions the implementation reached, with reasons

**Task helping, not fibers. Fibers are rejected, not deferred.** A thread that
waits on a counter executes other jobs while it waits, rather than parking. This
gets the property fibers were wanted for -- a nested `parallelFor` inside a job
must not deadlock the pool -- without a fiber runtime. (The earlier text here
named both Naughty Dog and enkiTS as the model; those are two different designs,
and enkiTS -- which helps rather than switching fibers -- is the one being
followed.)

Fibers do buy three things helping does not: bounded stack growth under deep
nesting, wait latency bounded by the actual dependency rather than by the longest
unrelated job a helper picks up, and freedom from the reentrancy hazard where a
helping thread runs a job conflicting with the one it is suspended inside.

They are still the wrong trade here. The cost is per-fiber stacks (hundreds of
64-512 KB allocations), degraded debugger and profiler support, sanitizer
annotations at every switch, and -- worst -- `thread_local` ceasing to be stable,
because a fiber may resume on a different thread. `currentLane()` is thread_local
today, and that pattern is everywhere in an engine. That is a broad class of
subtle bugs adopted across the whole codebase to fix a problem this workload does
not have: wide parallel-fors over homogeneous data, shallow nesting, comparable
job durations.

**The real answer is to make waiting rare rather than cheap.** Continuation-style
scheduling -- a job declares "when this counter drains, schedule that job" -- means
nothing blocks and there is nothing to help with. This is already what the DAG
below implies: if the scheduler derives ordering from declared read/write sets,
systems should not be calling `wait()` in the steady state at all. Helping then
survives only at the frame boundary and in tools and tests, which is exactly where
it is safe. See the multi-rate scheduler item in `06-roadmap.md` Phase 1.

**Bounded MPMC ring per lane, not a Chase-Lev deque.** Chase-Lev is the faster
structure and was implemented first. It is only safe when a queue slot is an
atomic that can be read speculatively and thrown away; a job record here is 80
bytes, so that speculative read is a plain struct copy racing with the owner's
write. ThreadSanitizer found it immediately. Claiming the slot with the CAS
before reading does not fix it either -- advancing `top` does not reserve the
slot, so once later thieves push `top` past it the owner may reuse it while the
first thief has still not copied it out. The Vyukov ring gives every cell a
sequence number that makes each access exclusive by construction.

Chase-Lev is worth recovering, and the route does **not** require epoch or
hazard-pointer reclamation. Store a **32-bit handle** (lane + index) in the deque
and keep payloads in the per-frame arena described in §5. The speculative read
then becomes a well-defined relaxed load of four bytes; a stale handle is
harmless because it is discarded on CAS failure; and a valid handle always points
at live storage, because nothing in the arena is recycled until the frame
quiesces. The allocator that makes this safe is already specified — the queue was
simply built before it existed. The MPMC ring is the correct thing to have in the
meantime, not a permanent answer.

**Grain matters more than either.** A `parallelFor` over 250 k elements at grain
97 produces 2578 chunks of roughly 1-2 µs each, against a queue overhead near
1 µs — that is the only regime where queue choice is visible at all, and the fix
is not a faster queue but auto-tuning grain toward the 50 µs target. Both claims
are currently folklore: there is no job throughput benchmark yet, and until there
is, the Chase-Lev question stays open on purpose rather than being settled by
assertion.

**Reductions bin per chunk, not per lane.** `parallelReduce` stores a partial per
chunk and folds them in chunk index order. Binning per lane would still be
non-deterministic, because a lane accumulates whichever chunks it happened to
win, in whatever order it won them. Tested by requiring bit-identical results
from 0, 1, 2, 7, 15 and 23 workers against the ordered single-threaded fold.

### What the benchmark said

`tools/job_bench` (run it; it prints the machine's load average, because these
numbers depend on it). Measured on the 24-thread development box under a load
average of ~11, so treat the absolute millisecond columns as the reliable
measurement and derived speedups as indicative — the single-threaded baseline
alone varied by 2.1× between runs, which the tool now reports rather than hides.

| Quantity | Measured |
|---|---|
| Dispatch, uncontended (0 workers) | **17–23 ns/job** |
| Dispatch, 16–23 workers, all jobs from one lane | ~250 ns/job |
| Chunk size at which efficiency plateaus | **≈ 2 µs** |
| Penalty for getting grain wrong | **~40×** (grain 16 vs the plateau) |

Three conclusions:

**Grain dominates everything else.** At a fixed worker count, the span between
the worst and best grain is around 40×. Nothing else in the job system comes
close to mattering that much, which is why grain auto-tuning was built and the
queue rewrite was not.

`parallelForAuto()` acts on this: it probes with a geometrically growing prefix
on the calling thread, derives cost per element, and picks a grain targeting
10 µs chunks — clamped so the chunk count stays between 2 and 64 per lane, for
tail balance at one end and dispatch cost at the other. The probe is real work,
not a discarded sample. Measured against the fixed-grain sweep on the same
workload, it matches or beats the best hand-picked grain (5.60 ms vs 6.07 ms at
23 workers) without being told anything about the body.

**Auto-grain is for `parallelFor` only, never `parallelReduce`.** The grain comes
from a wall-clock probe, so it varies between runs; a reduction chunked that way
would fold its partials in a different order each time and lose the bit-identical
property that the whole determinism argument rests on. `parallelReduce` therefore
requires an explicit grain and always will. This is the one place where the
performance work and the determinism requirement genuinely conflict, and
determinism wins.

**The Chase-Lev revisit is cancelled, on evidence.** Uncontended dispatch is
0.2% of a 10 µs chunk; even the worst-case contended figure is ~2.5%, and that
case is empty jobs all submitted from one lane, which maximises steal pressure
and does not occur with real work. At plateau grains the sweep shows no
dispatch-limited regime at all, so there is nothing for a faster queue to
recover. Revisit only if a profile shows queue contention.

**Efficiency falls off well before 23 workers.** 8 workers reach ~5.7×, 16 reach
~6.7×, 23 reach ~6.9×. How much of that is the job system and how much is a
machine already running at load 11 cannot be separated here; it needs a quiet
box before any conclusion is drawn about worker-count tuning.

### Verification

Concurrency correctness is not something functional tests establish on their own.
The job system is checked three ways, all in CI:

- 29 invariant assertions (`tests/test_jobs.cpp`): exactly-once execution, range
  tiling with no gaps or overlaps, three-deep nesting without deadlock at 0, 1, 2
  and 8 workers, ring overflow losing no work, and bit-identical reductions.
- ThreadSanitizer, via `cmake -DSHIPSIM_SANITIZE=thread`. This is what caught the
  Chase-Lev race; the functional tests passed the whole time it was there.
- AddressSanitizer and UBSan, via `cmake -DSHIPSIM_SANITIZE=address`.

## 2. Entity model

**Archetype ECS, structure-of-arrays storage.**

Entities are ships, compartments, openings, structural elements, fire cells, crew,
cargo items, particles. An archetype is a unique set of component types; all
entities sharing one live in contiguous SoA chunks of 16 KB.

Why archetypes rather than sparse sets: the hot loops here are numerical sweeps
over homogeneous collections (integrate 400 000 tetrahedra, evaluate 3 000 orifices,
advance 200 000 fire cells). Archetype chunks give those loops linear reads with
no indirection and vectorise without help. The price is that adding or removing a
component relocates the entity's row to a different archetype, so structural
changes belong at frame boundaries rather than in inner loops.

Implemented in `engine/core/ecs.cpp`. Entities are an index plus a generation, so
a handle to a destroyed entity is detectably stale rather than silently aliasing
whatever was created in its place. `each<Cs...>()` hands back contiguous parallel
arrays per chunk, which is the shape a numerical sweep wants.

Two things the tests are pointed squarely at, because both fail silently rather
than loudly: the swap-remove that fills a deleted row must patch the *moved*
entity's record (otherwise two handles alias one row, and nothing crashes), and
an archetype transition must carry over every component the two archetypes share
without resurrecting stale bytes for the ones they do not.

A `ComponentId` is only an array index: it is handed out on first use, so its
numeric value depends on which type was touched first and means nothing outside
the run that produced it. Everything that must be reproducible — archetype
identity, and therefore query iteration order — keys off a **stable id**, the
FNV-1a hash of the type's name. Archetype component lists are sorted by stable
id, so two builds agree on iteration order regardless of link order or which
type was registered first. Reflecting a component is what gives it a real name
to hash; anything unreflected falls back to the compiler's mangled name, which
is stable within a toolchain but not across them.

### Reflection and serialisation

`engine/core/reflect.hpp` describes types once — field names, kinds, offsets,
array extents, nested types — and `engine/core/serialise.cpp` consumes that for
the save format. Replication and editor panels are meant to consume the same
tables.

It is macro-based rather than parser-generated. A build-time parser over
annotated headers is still the right end state, but it needs libclang, a build
step, and a sync story, none of which earn their cost before anything consumes
the tables. The macros emit exactly the structures a generator would, so
swapping the front end later changes no consumer.

The wire format is self-describing per field — name hash, kind, count, byte
length — which costs bytes and buys the property that actually matters over a
project's lifetime: **a save written by one build loads into another**. Unknown
fields are skipped by their recorded length; fields the writer lacked keep the
reader's existing value; neither side has to agree on field order. Scalars are
written explicitly little-endian rather than memcpy'd, so a save file is portable
rather than accidentally x86-shaped, and floats go through their IEEE-754 bit
pattern so negative zero and NaN payloads survive exactly.

The reader trusts nothing. Every read is bounds-checked, the field count is
capped at what the remaining bytes could describe rather than believed as a loop
bound, and each field resynchronises on its recorded length so one changed nested
struct cannot desynchronise the rest of the stream. Tested by feeding it every
truncation of a valid buffer and requiring all of them to be rejected — with
AddressSanitizer proving it did not step outside the buffer while rejecting.

### World save/load

`World::save` / `World::load` put a whole ECS world through that codec.

**Entity handles survive.** The generation table is written for every index, not
just live ones, because handle staleness is world state: an `Entity` that was
destroyed before the save must still read as dead afterwards, and the free list
must come back so the next `create()` reuses indices the same way. Without this
a reload silently resurrects destroyed entities.

**Unknown components do not take the entity with them.** A save may contain a
component this build has no code for — an older or newer build, or a disabled
module. The component is skipped and the entity loads with everything else
intact. Reflected components are written field-wise and get the schema tolerance
above; unreflected ones go out as size-tagged opaque blobs, which is build-local
and refused on load if the size no longer matches, because reading it in anyway
would be a silent reinterpretation.

**Load fails closed.** The first implementation returned `false` on a malformed
stream but left the entities it had already decoded in place — 554 of the
truncation cases produced a half-built world. A caller that ignores the return
value would then be simulating a ship missing an arbitrary suffix of itself.
Loading now clears on any failure path.

**Save is a fixed point.** Save, load, save again is byte-identical, and so is
the third pass. That is what makes saves diffable and usable as a network
baseline, and it is a sharp test: anything that reorders, drops or normalises
during a round trip shows up immediately.

Components are plain data with no virtual functions. Behaviour lives in systems.
Reflection is generated by a small build-time parser over annotated headers and
drives serialisation, network replication, the editor's property panels, and the
save format — one description, four consumers.

## 3. The tick

**Fixed 100 Hz simulation tick, decoupled variable-rate rendering, state
interpolation for presentation.**

Rendering never observes a partially updated simulation. The sim writes into a
double-buffered snapshot; the renderer reads the previous one and interpolates
toward it. This costs memory and buys the freedom to run the sim on all cores
while the render thread is doing something else entirely.

### Multi-rate solvers

Nothing runs at one rate. Each solver declares its own step and the scheduler
subdivides:

| Solver | Rate | Notes |
|---|---|---|
| Rigid body, hydrostatics | 100 Hz | cheap, already implemented |
| Flooding / air orifice network | 100 Hz | sub-steps adaptively when a large breach opens |
| Structural FEM, quiescent | 100 Hz | reduced modal model, ~200 DOF |
| Structural FEM, active damage zone | 2–20 kHz | explicit central difference, CFL-limited by element size |
| Thermal conduction | 10 Hz | implicit, unconditionally stable, large steps fine |
| Fire / gas transport | 20–50 Hz | operator-split advection + reaction |
| Particle fluid (sloshing, spray) | 100–200 Hz | substepped XPBD or FLIP |
| Wave field | 30 Hz | FFT regeneration; evaluation is analytic at any time |
| Crew / AI | 10 Hz | staggered across entities |
| Network | 20 Hz | snapshot + event stream |

The active-FEM rate is the hard one: an explicit solver on 10 mm elements in steel
has a stable step around 1.5 µs, which is ~700 substeps per 1 ms of sim time. This
is why the damage zone is spatially adaptive and time-boxed rather than always-on
(see `02-simulation.md` §3).

Implemented in `engine/core/scheduler.cpp`. Each system declares a rate in
simulation time, a dilation band it can honour, and a catch-up ceiling; systems
declare prerequisites and the scheduler derives dependency levels, running
everything within a level in parallel on the job system.

**Simulation time is accumulated in integer nanoseconds.** Floating point is used
exactly once per advance, converting the wall delta; every step-count decision
after that is integer division. This is not fastidiousness — the first version
accumulated seconds in a `double`, and 600 advances of 1/60 s summed to fractionally
under 10 s, so a 100 Hz system ran 999 times instead of 1000. Over a forty-minute
flooding casualty, let alone a campaign, that drift compounds and makes the
schedule depend on how the frame times happened to be sliced. With an integer
clock the remainder carries exactly and a given sequence of advances produces the
same schedule on any platform, however long it runs.

Two related properties fall out and are tested: step counts for the whole frame
are decided *before* any system runs, so a system that overran cannot change how
often its neighbours update; and the catch-up ceiling discards its surplus rather
than carrying it, since carrying it only guarantees the same overrun next frame.

### Time dilation as a first-class feature

A structural failure worth watching happens over 50 ms. Flooding worth watching
happens over 40 minutes. The engine treats the sim/wall clock ratio as a
controllable input, and solvers advertise which rates they can honour at a given
dilation. Accelerating to 60× drops the FEM to its reduced model and disables
particle fluids; slowing to 0.02× enables everything. The player-facing feature —
slow-motion on a hull rupture, fast-forward through a long flooding fight — is the
same mechanism the developers use to make the expensive solvers affordable.

## 4. Determinism

Required for the multiplayer model in `04-multiplayer.md` and, more immediately,
for reproducing bugs in a system with 10⁶ interacting state variables.

- `-ffp-contract=off`, no `-ffast-math`, no x87. SSE2/AVX only, which is
  IEEE-754-exact.
- No `std::unordered_map` iteration in anything that affects state. Ordered or
  index-based containers only.
- Deterministic parallel reduction as described above.
- The RNG is an explicit, seeded, per-system counter-based generator
  (Philox/Threefry), never a global.
- A replay is a seed plus the input stream. CI re-runs a corpus of recorded
  scenarios every night and diffs the final state bit-for-bit.

Transcendental functions are the residual risk (`libm` differs across platforms
and versions). The sim uses an in-tree implementation of `sin`, `cos`, `exp`,
`log`, `pow`, `atan2` and `sqrt`-adjacent helpers for anything on the state path.

**The boundary of this guarantee, measured rather than assumed.** Determinism
covers the CPU simulation path. It does **not** extend to GPU field solvers. The
FEM spike (`07-fem-spike-findings.md` §2) shows the GPU and CPU kernels agreeing
to 2 × 10⁻⁵ on a single step — the kernel is correct — and then diverging to
5 × 10⁻³ over a thousand steps, because the GPU contracts multiply-add pairs
differently and an explicit scheme at the CFL limit amplifies that in its
highest-frequency modes. This is a property of explicit dynamics, not a bug to be
fixed. Anything requiring reproducibility must therefore consume *events* from
the field solvers (a tear, its geometry, the resulting orifice area) rather than
their raw state.

## 5. Memory

- **Arena allocators per frame and per lane** — implemented in
  `engine/core/arena.cpp`. The global allocator is off-limits during a tick: it
  takes a lock and it is unpredictable, and a 100 Hz tick that occasionally stalls
  in `malloc` is a 100 Hz tick that occasionally misses. `FrameArenas` gives one
  arena per job-system lane, so allocation needs no atomics at all — the arena
  itself is deliberately not thread-safe.
- Pool allocators for the churn-heavy things (particles, fire cells, FEM zones).
- Chunked SoA for entity data as above.
- Asset streaming is virtual-texture style for the exterior world; ship interiors
  are resident because the player can be anywhere in them and a hitch while
  flooding is unacceptable.

Budget target for a 300 m vessel: 6 GB resident, 2 GB of that simulation state.
Arena sizing is meant to come from measured `highWaterMark()`, not from guesses.

### Two things the arena implementation settled

**Align the address, not the offset.** The obvious implementation rounds the
running offset up to the requested alignment. That is correct only while the
requested alignment does not exceed the block's own alignment; past that it
silently returns under-aligned memory, and nothing faults until someone does an
aligned SIMD store. Caught by testing alignments from 1 to 256 rather than the
handful anyone writes by hand.

**Arenas are poisoned under AddressSanitizer.** A bump allocator is one large
valid heap block as far as ASan is concerned, so an overrun past an allocation,
or a read after `reset()`, is an invisible scribble inside memory ASan considers
legitimately ours — precisely the class of bug arenas make easier to write. The
unused tail is kept poisoned and only the bytes actually handed out are
unpoisoned, including the alignment padding, so overruns land in poisoned space.
Verified with a negative control: a 64-byte write into a 16-byte arena allocation
is reported as use-after-poison.

**No destructors, ever.** `create<T>()` refuses anything not trivially
destructible rather than silently leaking it, and `allocateArray<T>()` avoids
array placement-new entirely — the standard permits it to demand an unspecified
cookie, which would overrun an allocator sized to `sizeof(T) * count`.

## 6. Module boundaries

```
engine/core/      math, containers, allocators, jobs, reflection, serialisation
engine/sim/       hydrostatics, flooding, seakeeping, FEM, thermal, fire, fluids,
                  aero, machinery — each a separate library with no sibling deps
engine/render/    Vulkan backend, render graph, passes
engine/audio/     DSP graph, propagation
engine/net/       replication, prediction, transport
engine/editor/    ship editor, scenario editor, debug visualisers
game/             ship definitions, scenarios, UI, campaign
```

Simulation libraries depend only on `core`. They communicate through explicitly
declared coupling interfaces registered with the scheduler — the FEM solver never
includes a flooding header. Coupling is the part that will rot if it is allowed to
be implicit, so it is the part with the most ceremony.

## 7. Coupling

The physically interesting behaviour is almost entirely in the couplings:

| From | To | Mechanism |
|---|---|---|
| FEM | Flooding | a torn element becomes an orifice with the area of the tear |
| Flooding | FEM | hydrostatic pressure on wetted internal surfaces becomes a nodal load |
| Flooding | Rigid body | floodwater mass and centroid |
| Rigid body | Flooding | attitude sets every free-surface plane |
| Fire | FEM | steel temperature scales yield strength and Young's modulus |
| FEM | Fire | a breached bulkhead is a new gas flow path |
| Waves | FEM | wave-induced hull girder bending, hogging and sagging |
| Flooding | Fire | water suppresses; steam displaces oxygen |
| Aero | Rigid body | wind heeling moment on the superstructure |

Couplings are evaluated **explicitly** (Gauss–Seidel over subsystems within a
tick, in a fixed declared order) rather than monolithically. This is the standard
partitioned approach in multiphysics and it is stable here because the coupling
stiffnesses are modest relative to the individual solvers' internal stiffness.
The one exception is floodwater ↔ rigid body, which is stiff enough that the
current implementation already re-levels free surfaces twice per tick (before and
after the flow solve) — a cheap predictor-corrector.
