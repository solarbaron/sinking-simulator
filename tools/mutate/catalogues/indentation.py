"""The rigid-plastic indentation model, one plausible edit at a time.

`engine/sim/indentation.cpp` is 247 lines and every quantity in it has a closed
form, three of them inverses of each other. That is the strongest structure a
mutation catalogue can be handed: a substitution either moves a number an inverse
disagrees with, or it does not move anything at all.

It is also the file whose four closed forms were rewritten today (`cce7a16`)
because each was evaluated as a difference of two nearly equal lengths and lost
its digits for a shallow dent -- `sqrt(1 + r^2) - 1` is out by 2.3e-2 at a
0.12 um dent, which is no digits at all. The conjugate forms that replaced them
are the *same numbers*, not approximations, so **the conjugate rewrite is
itself a mutation with a known answer**: putting any of the four back is an edit
that changes no published figure, passes every test this file had before today,
and must now be killed by the ulp-bounded round trips that arrived with the fix.
Mutants 2, 3, 4 and 5 are exactly those four reversions. If the suite ever stops
killing them the repair has quietly come undone.

Three substitutions are **controls** (`expect: survive`) and each carries its
equivalence argument at the site. A control that dies is reported as loudly as a
survivor: it means either the argument is wrong or the harness is measuring
something other than the mutation.

    ./tools/mutate/mutate.py list indentation
    ./tools/mutate/mutate.py scan indentation
    ./tools/mutate/mutate.py run  indentation --workers 4
"""
NAME = "indentation"

I = "engine/sim/indentation.cpp"

# No bound overrides. The defaults are `max(floor, factor x measured baseline)`,
# and the floor is there for a suite so fast that four times it is a fraction of a
# second -- it is a minimum, not a cap. This suite measures 137.6 s CPU and 141.9 s
# wall for 200 937 checks on an idle box, so 4x puts the CPU bound at about 550 s
# and the wall backstop at about 2840 s, both far above the 30 s and 300 s floors.
# Setting the floors higher would only *weaken* hang detection, and setting the
# factors from a constant would be asserting a wall clock, which is the mistake the
# gas tier's `coreSecondsPerCell` made.

MUTANTS = [
    # --- the half-span everything else is written about --------------------------
    ("the half-span is the whole span", I,
     "double half(double span) { return 0.5 * std::max(span, 1e-12); }",
     "double half(double span) { return std::max(span, 1e-12); }", "kill"),

    # CONTROL 1. `x * 0.5` and `x / 2.0` are the *same IEEE-754 operation* on this
    # input: 2 is a power of two, so x/2 is exact whenever it does not underflow
    # into the subnormal range, and both forms therefore round the identical exact
    # value. `std::max(span, 1e-12)` is at least 1e-12, so the result is at least
    # 5e-13 -- 295 orders above the subnormal boundary at 2.2e-308 -- and division
    # by 2 cannot overflow at any magnitude. This is bit-identical for every input
    # the function can be given, which is a stronger claim than "equivalent within
    # tolerance" and is the one a control should make where it can.
    ("CONTROL: the half-span is written as a division rather than a halving", I,
     "double half(double span) { return 0.5 * std::max(span, 1e-12); }",
     "double half(double span) { return std::max(span, 1e-12) / 2.0; }", "survive"),

    # --- the four conjugate forms, put back the way they cancelled ---------------
    # Each of these is *algebraically exact* and was the code until today. None of
    # them moves a published figure. They are killable only by an assertion that
    # measures the arithmetic rather than the algebra, which is what arrived with
    # the fix: two round trips bounded at 8 ulps and two small-dent series checks
    # bounded at 8 ulps, all of them geometric down to a micron.
    ("the membrane strain cancels as sqrt(1 + r^2) - 1 again", I,
     "    return ratio * ratio / (1.0 + std::sqrt(1.0 + ratio * ratio));",
     "    return std::sqrt(1.0 + ratio * ratio) - 1.0;", "kill",
     "algebraically exact; 2.2e-5 relative at a micron dent against an 8 ulp bound"),
    ("the strain inverse cancels as (1 + eps)^2 - 1 again", I,
     "    return half(span) * std::sqrt(strain) * std::sqrt(strain + 2.0);",
     "    return half(span) * std::sqrt((1.0 + strain) * (1.0 + strain) - 1.0);", "kill",
     "algebraically exact; loses the strain below 2.2e-16 absolute"),
    ("the energy cancels as sqrt(h^2 + d^2) - h again", I,
     "    return 2.0 * p.yieldStrength * p.thickness * p.contactWidth *\n"
     "           (penetration * penetration / (std::sqrt(h * h + penetration * penetration) + h));",
     "    return 2.0 * p.yieldStrength * p.thickness * p.contactWidth *\n"
     "           (std::sqrt(h * h + penetration * penetration) - h);", "kill",
     "algebraically exact; 3.3e-5 relative at a micron dent"),
    ("the energy inverse cancels as (u + h)^2 - h^2 again", I,
     "    const double penetration = std::sqrt(u) * std::sqrt(u + 2.0 * h);",
     "    const double penetration = std::sqrt((u + h) * (u + h) - h * h);", "kill",
     "algebraically exact; differences two O(h^2) numbers whose difference is 2 h u"),

    # CONTROL 2. The source argues this equivalence itself: "Split as two roots
    # rather than sqrt(u * (u + 2h)): both factors are positive ... it is the same
    # value, and neither factor can overflow when the other is large." The split
    # is an overflow guard, not an accuracy one -- `u * (u + 2h)` only overflows
    # for u above about 1.3e154 m, a penetration no ship has -- and neither form
    # cancels, because `u` and `u + 2h` are both formed without a subtraction.
    # What separates them is one rounding: sqrt is correctly rounded, so the two
    # forms differ by the rounding of a single multiplication.
    #
    # **Measured, not argued.** Re-running the test's own sweep -- 601 geometric
    # depths from 2.4 um to the tearing penetration -- outside the suite gives the
    # shipped two-root form a worst round trip of 1.070 ulp, which is the figure
    # `test_indentation.cpp` publishes, and the product form 0.986 ulp: the control
    # is if anything a hair *better* here, and both sit against a bound of 8. The
    # small-dent series check reads 0.000 ulp shipped and -1.000 ulp mutated. So the
    # claim is "the same value to within one rounding at every depth this model is
    # used at", with 7 ulps of headroom, and it is a measurement.
    #
    # The contrast that makes the control worth having: the *cancelling* form of the
    # same inverse, mutant 5 above, reads 1.33e11 ulp on that identical sweep. One
    # substitution to this line must die and the other must live, and the file cannot
    # tell them apart by looking -- both are exact algebra.
    ("CONTROL: the energy inverse takes one root of the product, not two roots", I,
     "    const double penetration = std::sqrt(u) * std::sqrt(u + 2.0 * h);",
     "    const double penetration = std::sqrt(u * (u + 2.0 * h));", "survive",
     "the source argues this equivalence; the split exists for overflow, not accuracy"),

    # --- the force, and the tent it resolves ------------------------------------
    ("only one leg of the tent carries the membrane tension", I,
     "    return 2.0 * p.yieldStrength * p.thickness * p.contactWidth * penetration / slant;",
     "    return 1.0 * p.yieldStrength * p.thickness * p.contactWidth * penetration / slant;",
     "kill"),
    ("the slant is a difference of squares rather than a sum", I,
     "    const double slant = std::sqrt(h * h + penetration * penetration);",
     "    const double slant = std::sqrt(h * h - penetration * penetration);", "kill"),

    # --- tearing ----------------------------------------------------------------
    ("a plate tears below its failure strain and heals above it", I,
     "    s.torn = p.failureStrain > 0 && s.strain >= p.failureStrain;",
     "    s.torn = p.failureStrain > 0 && s.strain <= p.failureStrain;", "kill"),
    ("the energy to tear is measured at half the failure strain", I,
     "    return indentationEnergy(p, penetrationForStrain(p.span, p.failureStrain));",
     "    return indentationEnergy(p, penetrationForStrain(p.span, 0.5 * p.failureStrain));",
     "kill"),
    ("a torn plate keeps resisting, so more energy drives the dent deeper", I,
     "    return tearing > 0 ? std::min(penetration, tearing) : penetration;",
     "    return tearing > 0 ? std::max(penetration, tearing) : penetration;", "kill"),

    # --- the strike against a real hull ------------------------------------------
    ("a strike with no energy at all is admitted", I,
     "    if (!(radius > 0) || !(energy > 0)) return damage;",
     "    if (!(radius > 0) || !(energy >= 0)) return damage;", "kill"),
    ("the radius stops bounding how far the damage reaches", I,
     "        if (distance > radius) continue;",
     "        if (distance > 2.0 * radius) continue;", "kill"),
    ("the bow spends its energy on the farthest panel first", I,
     "    std::sort(reachable.begin(), reachable.end());",
     "    std::sort(reachable.rbegin(), reachable.rend());", "kill"),

    # CONTROL 3. `reachable` holds `std::pair<double, int>` and pairs compare
    # lexicographically, so two entries compare equivalent only if they agree in
    # *both* members. The second member is the panel index and each panel is
    # pushed exactly once, so no two entries can be equivalent: the order is a
    # strict total order, the sorted sequence is unique, and a stable sort and an
    # unstable one produce the identical vector. This is the same question today's
    # `water_promotion` tie-break commit asked one file away, and here it has the
    # opposite answer -- which is worth having a control say rather than assuming.
    ("CONTROL: the panel order is sorted stably", I,
     "    std::sort(reachable.begin(), reachable.end());",
     "    std::stable_sort(reachable.begin(), reachable.end());", "survive",
     "pairs carrying a distinct index admit no ties, so the sorted order is unique"),

    ("the plate spans the long way between its supports", I,
     "        const double span = std::min(frameSpacing, stiffenerSpacing);",
     "        const double span = std::max(frameSpacing, stiffenerSpacing);", "kill",
     "the defect ce02678's neighbour fixed: authored as the frame spacing"),
    ("the struck width multiplies by the span instead of dividing", I,
     "        model.contactWidth = std::max(panel.area() / std::max(span, 1e-6), 1e-6);",
     "        model.contactWidth = std::max(panel.area() * std::max(span, 1e-6), 1e-6);", "kill"),
    ("the failure strain is regularised with its two lengths swapped", I,
     "            plasticity::regularisedFailureStrain(material.failure, span, panel.thickness);",
     "            plasticity::regularisedFailureStrain(material.failure, panel.thickness, span);",
     "kill"),
    ("a torn bay costs the bow half the energy it took to tear it", I,
     "            remaining -= toTear;",
     "            remaining -= 0.5 * toTear;", "kill"),
    ("the deepest dent is reported as the shallowest", I,
     "            damage.penetration = std::max(damage.penetration, reached);",
     "            damage.penetration = std::min(damage.penetration, reached);", "kill"),
    ("every strike reports all its energy as unspent", I,
     "    damage.energyUnspent = remaining;",
     "    damage.energyUnspent = energy;", "kill"),

    # --- the model's own refusals ------------------------------------------------
    ("a span only a few thicknesses long is no longer refused", I,
     "    if (p.thickness > 0 && p.span / p.thickness < 20.0)",
     "    if (p.thickness > 0 && p.span / p.thickness < 0.0)", "kill",
     "the bending-resistance guard: span/t below 20 is a block, not a membrane"),
]
