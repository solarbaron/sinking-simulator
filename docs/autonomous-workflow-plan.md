# Autonomous Workflow Plan

## Current Status (Session Start)

### Completed
- ✅ CO2 and inert-gas total flooding (merged, pushed to GitHub)
- ✅ FLIP/APIC solver (merged, pushed to GitHub)
- ✅ Horizontal openings with buoyant exchange (debugged, merged, pushed to GitHub)
- ✅ All 200,249 checks passing, warning-clean build
- ✅ Repository at github.com:solarbaron/sinking-simulator.git

### Active Work
- 🔄 3 background agents analyzing codebase:
  - Agent analyzing promotion patterns and Ship integration points
  - Agent studying architecture docs (memory, performance, constraints)
  - Agent running comprehensive test analysis

### Phase 5 Priority: FLIP Escalation
**Goal:** Enable quiescent ↔ dynamic water transitions to reach milestone "stand in a flooding compartment as it fills around you, in a rolling ship, and have the water behave."

## Autonomous Work Strategy

### 1. Multi-Agent Parallelism

**Pattern:**
- Spawn 3-5 agents simultaneously on independent subtasks
- Design/analysis agents: read-only, no conflicts
- Implementation agents: separate files or worktree isolation
- Synthesis agent: consolidates findings after completion

**Current deployment:**
```
Main: Planning & coordination
├── Agent 1: Phase 5 analysis (COMPLETED)
├── Agent 2: Code scanning (IN PROGRESS)
├── Agent 3: Test analysis (IN PROGRESS)
└── Agent 4: Architecture study (IN PROGRESS)
```

**Next wave:**
```
Main: Synthesis & implementation planning
├── Agent 5: WaterPromoter class implementation
├── Agent 6: Ship integration wiring
├── Agent 7: State transfer functions
├── Agent 8: Test suite for escalation
└── Agent 9: Documentation updates
```

### 2. Task Tracking

Using TaskCreate/TaskUpdate/TaskList to maintain work backlog:
- Each major feature → top-level task
- Implementation steps → subtasks via description
- Status transitions: pending → in_progress → completed
- Blocked tasks track dependencies

**Current tasks:**
1. ✅ Analyze Phase 5 priorities
2. 🔄 Code improvements scan
3. 🔄 Test analysis
4. 🔄 Architecture study
5. ⏳ Design escalation criterion
6. ⏳ Plan Ship-FLIP integration
7. ⏳ Synthesize findings
8. ⏳ Autonomous workflow setup

### 3. Git Workflow

**Branch strategy:**
- `main` — stable, tested, pushed to GitHub
- `feature/flip-escalation` — work-in-progress
- Worktrees for parallel agent development (when needed)

**Commit cadence:**
- Commit after each logical unit (function, test, file)
- Use `git commit -F` with descriptive messages
- Push to GitHub after tests pass
- Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>

### 4. Testing Discipline

**Before any commit:**
```bash
ninja -C build && ./build/shipsim_tests
```

**For escalation work:**
```bash
# Round-trip conservation
./build/shipsim_tests 2>&1 | grep "flip.*escalation\|water.*promote"

# Performance regression
time ./build/shipsim_tests > /dev/null
```

**Integration validation:**
```bash
./scripts/verify.sh full  # if it exists
```

### 5. Autonomous Loop

**Using cron-style scheduling:**

```javascript
// Every 30 minutes: check progress, spawn new agents
CronCreate({
  cron: "27 * * * *",  // :27 past each hour (off-peak)
  prompt: `Check autonomous workflow progress:
1. List completed agents and extract findings
2. Identify next actionable tasks
3. Spawn new agents for unblocked work
4. Report status summary
5. Commit completed work to git`,
  recurring: true
})
```

**Fallback: Manual check-in prompts**

If cron isn't suitable, use explicit instructions:
- "Continue autonomous work"
- "Status update and next steps"
- "What's blocking progress?"

### 6. Decision Points & Escalation

**Auto-proceed:**
- Implementation following established patterns (GasPromoter → WaterPromoter)
- Test additions for new features
- Documentation updates
- Performance measurements
- Bug fixes with clear root cause

**Require user confirmation:**
- New architectural patterns
- Breaking API changes
- Performance tradeoffs >10%
- External dependencies
- Deviations from roadmap

**Block and report:**
- Test failures (any failure blocks merge)
- Build errors
- Conflicting agent findings
- Missing design decisions

### 7. Progress Metrics

**Track:**
- Lines added/changed per session
- Tests added/passing
- Features completed vs. roadmap
- Agent utilization (spawned vs. completed)
- Commit frequency

**Report periodically:**
```markdown
## Session Progress

- **Duration:** 2h 15m
- **Agents spawned:** 9 (7 completed, 2 active)
- **Code changes:** +847 lines (3 files)
- **Tests:** +12 new, 200,261 total passing
- **Commits:** 4 to feature branch
- **Status:** WaterPromoter implemented, wiring in progress
```

### 8. Context Management

**Preserve across compaction:**
- Task list (TaskList)
- Git state (current branch, uncommitted work)
- Active agent IDs and topics
- Open design questions

**Reset after compaction:**
- Re-check test status
- Re-verify build passes
- Re-read task list
- Don't assume agent completion

## Implementation Plan: FLIP Escalation

### Phase 1: Foundation (1-2 days)

**Files to create:**
- `engine/sim/water_promotion.hpp` — WaterCriterion, WaterPromoter class
- `engine/sim/water_promotion.cpp` — implementation
- `tests/test_water_promotion.cpp` — unit tests

**Dependencies:**
- Read: promotion.hpp, flip.hpp, ship.hpp
- Pattern: GasPromoter (§6 of promotion.hpp)

**Agents:**
- Agent A: Implement WaterCriterion + WaterCandidate structs
- Agent B: Implement WaterPromoter state machine
- Agent C: Write unit tests for criterion and hysteresis

### Phase 2: State Transfer (1-2 days)

**Files to modify:**
- `engine/sim/ship.hpp` — add WaterPromoter member, activeWater_ map
- `engine/sim/ship.cpp` — escalation logic in step()
- `engine/sim/flip.hpp` — add API for state transfer if needed
- `engine/sim/flip.cpp` — implement state transfer helpers

**Functions to add:**
```cpp
// In ship.cpp
flip::Solver* promoteWater(Compartment& c, const WaterCriterion& crit);
void demoteWater(Compartment& c, flip::Solver* solver);
void applyWaterReview(const WaterReview& review);
```

**Agents:**
- Agent D: Implement promoteWater (Compartment → FLIP)
- Agent E: Implement demoteWater (FLIP → Compartment)
- Agent F: Wire into Ship::step() with review cadence

### Phase 3: Testing & Validation (2-3 days)

**Files to create:**
- `tests/test_ship_flip.cpp` — integration tests
- `tools/water_escalation_probe.cpp` — scenario tool (optional)

**Test scenarios:**
1. Round-trip mass conservation (exact)
2. Escalation on ship roll
3. Demotion when motion stops
4. Budget enforcement
5. No chatter across threshold
6. Multiple compartments

**Agents:**
- Agent G: Write integration tests
- Agent H: Run performance measurements
- Agent I: Update documentation

### Phase 4: Milestone Validation (1 day)

**Verify:**
- Flooding compartment escalates to FLIP during ship motion
- Water "behaves" (sloshes, responds to ship roll)
- Demotes cleanly when ship settles
- Performance acceptable (<10 core-s/sim-s per compartment)

**Deliverable:**
- `tools/flooding_scenario` demonstrating milestone
- Video/screenshots (if renderer available)
- Updated roadmap marking escalation complete

## Estimated Timeline

- **Phase 1-2:** 2-4 days (foundation + wiring)
- **Phase 3:** 2-3 days (testing)
- **Phase 4:** 1 day (validation)
- **Total:** 5-8 days with 3-5 agents working in parallel

## Success Criteria

1. ✅ WaterPromoter follows GasPromoter pattern
2. ✅ Round-trip mass conservation exact (0.0 tolerance)
3. ✅ No test failures, warning-clean build
4. ✅ All checks pass (>200,000)
5. ✅ Performance within budget (<10 core-s/sim-s)
6. ✅ Milestone scenario demonstrates water behavior
7. ✅ Committed to git, pushed to GitHub
8. ✅ Documentation updated (roadmap, architecture docs)

## Blockers & Risks

**Known:**
- FLIP boundary conditions for actual TriMesh (vs AABB) — may need Phase 2
- Ship motion in FLIP frame — coordinate transformation TBD
- Threading for multiple FLIP solvers — serial initially, optimize later

**Monitoring:**
- Agent completion rate (stalls indicate design gap)
- Test failure rate (indicates integration issues)
- Build time (indicates scope creep)

## Next Actions (Immediate)

1. ⏳ Wait for 3 active agents to complete
2. 📋 Synthesize findings into flip-escalation-design.md
3. 🚀 Spawn Phase 1 agents (WaterPromoter implementation)
4. 🔄 Set up autonomous loop (30-min check-ins)
5. 📊 Report progress after first implementation wave

---

**Last updated:** Start of autonomous work session
**Active agents:** 3 (code scan, tests, architecture)
**Next review:** After agent completion notifications
