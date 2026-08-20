"""§7 of docs/07-fem-spike-findings.md, re-derived: 68 mutants of the plasticity work.

That section publishes "68 mutants, each a single plausible edit... the suite went
from 58 to 65 kills", and the sweep that produced it was a throwaway script that
was deleted. An audit could not re-run one of them. This is the catalogue that
should have been checked in beside the figure.

**It is not the original list, and it cannot be.** The original substitutions were
never written down anywhere, so what is reproducible from the document is the
*population* -- single plausible edits to `plasticity.cpp` and to the solid-shell
element's plasticity hooks -- and, exactly, the nine mutants §7 describes in
words. Those nine are marked `DOC` below and are the sharp part of the
re-derivation: six are holes the section says its new tests now close, and three
are the equivalents it argues must survive. A rate agreeing to a mutant would be a
coincidence; the nine agreeing one by one would not be.

The three `DOC` equivalents are carried as **controls** (`survive`). They are the
harness's own negative controls: if the suite kills one, either the equivalence
argument in §7 is wrong or the harness is measuring something other than the
mutation, and both need saying out loud.
"""
NAME = "plasticity"

P = "engine/sim/plasticity.cpp"
S = "engine/sim/solid_shell.cpp"

MUTANTS = [
    # --- the named constants and the small algebra everything else rests on -------
    ("kRoot23 is sqrt(3/2) rather than sqrt(2/3)", P,
     "const double kRoot23 = std::sqrt(2.0 / 3.0);",
     "const double kRoot23 = std::sqrt(3.0 / 2.0);", "kill"),
    ("the tensor norm loses the factor two on its shear components", P,
     "           2.0 * (v[3] * v[3] + v[4] * v[4] + v[5] * v[5]);",
     "           1.0 * (v[3] * v[3] + v[4] * v[4] + v[5] * v[5]);", "kill"),
    # DOC: "The failure plane's normal was only ever asked for on an axis-aligned
    # pull, where the stress tensor is already diagonal, the Jacobi sweep has nothing
    # to do and the eigenvectors come back as the identity however badly they are
    # accumulated. Freezing them entirely passed. The same tear, rotated 0.9 rad
    # about an arbitrary axis, kills it."
    ("DOC: the eigenvectors are frozen at the identity", P,
     "        if (off <= 0.0) break;",
     "        if (off >= 0.0) break;", "kill"),
    ("the Swift strength coefficient multiplies where it should divide", P,
     "    curve.strengthCoefficient = yieldStrength / std::pow(referenceStrain, exponent);",
     "    curve.strengthCoefficient = yieldStrength * std::pow(referenceStrain, exponent);",
     "kill"),
    ("the strength ratio is inverted inside the logarithm", P,
     "        return n * std::log(n / eps0);",
     "        return n * std::log(eps0 / n);", "kill"),
    ("the bisection moves the wrong bracket", P,
     "        if (ratio(mid) > target)\n            lo = mid;\n        else\n            hi = mid;",
     "        if (ratio(mid) > target)\n            hi = mid;\n        else\n            lo = mid;",
     "kill"),
    # DOC: "The geometric midpoint in the Swift fit is conditioning rather than
    # correctness: 200 arithmetic bisections reach the same root."
    ("DOC CONTROL: the Swift fit bisects arithmetically rather than geometrically", P,
     "        const double mid = std::sqrt(lo * hi);",
     "        const double mid = 0.5 * (lo + hi);", "survive"),
    ("the Swift flow stress forgets it starts at the yield point", P,
     "    return curve.strengthCoefficient *\n"
     "           std::pow(curve.referenceStrain + p, curve.hardeningExponent);",
     "    return curve.strengthCoefficient *\n"
     "           std::pow(p, curve.hardeningExponent);", "kill"),
    ("the Swift slope is not the derivative of the Swift stress", P,
     "    return curve.strengthCoefficient * curve.hardeningExponent *\n"
     "           std::pow(curve.referenceStrain + p, curve.hardeningExponent - 1.0);",
     "    return curve.strengthCoefficient * curve.hardeningExponent *\n"
     "           std::pow(curve.referenceStrain + p, curve.hardeningExponent);", "kill"),
    ("linear hardening reports no incremental stiffness", P,
     "        case Hardening::Linear: return curve.hardeningModulus;",
     "        case Hardening::Linear: return 0.0;", "kill"),
    ("the Swift uniform elongation adds the reference strain back", P,
     "    return std::max(0.0, curve.hardeningExponent - curve.referenceStrain);",
     "    return std::max(0.0, curve.hardeningExponent + curve.referenceStrain);", "kill"),

    # --- failure: regularisation and triaxiality ------------------------------------
    ("the mesh regularisation is upside down", P,
     "        share = std::min(1.0, thickness / elementLength);",
     "        share = std::min(1.0, elementLength / thickness);", "kill"),
    ("the regularisation interpolates from the wrong end", P,
     "    return failure.uniformStrain + (failure.fractureStrain - failure.uniformStrain) * share;",
     "    return failure.fractureStrain + (failure.uniformStrain - failure.fractureStrain) * share;",
     "kill"),
    ("the triaxiality cutoff admits the states it exists to exclude", P,
     "    if (eta <= failure.cutoffTriaxiality) return std::numeric_limits<double>::infinity();",
     "    if (eta >= failure.cutoffTriaxiality) return std::numeric_limits<double>::infinity();",
     "kill"),
    ("triaxiality makes steel tougher rather than more brittle", P,
     "    return std::exp(-failure.triaxialitySensitivity * (eta - failure.referenceTriaxiality));",
     "    return std::exp(failure.triaxialitySensitivity * (eta - failure.referenceTriaxiality));",
     "kill"),
    ("ship steel is fitted to engineering stress rather than true stress", P,
     "    material.flow = swiftFromTensile(355.0e6,\n"
     "                                     engineeringUltimate * (1.0 + engineeringUniform),",
     "    material.flow = swiftFromTensile(355.0e6,\n"
     "                                     engineeringUltimate,", "kill"),
    # The comment above this line in the source records that deleting it *survived*
    # the first round of mutation testing, and that the modulus sweep now carries a
    # case it fails. That is a claim the source makes about the suite, and this is
    # the mutant that checks it. `plastic` is `max(0, .)`, so the guard below can
    # never fire: this is the deletion, written so that it leaves a mark.
    ("the exact early return in the secant shear modulus is disabled", P,
     "    if (plastic == 0.0) return shear;",
     "    if (plastic < 0.0) return shear;", "kill"),
    ("the tangent modulus divides by three instead of into it", P,
     "    return 1.0 / (1.0 / shear + 3.0 / slope);",
     "    return 1.0 / (1.0 / shear + slope / 3.0);", "kill"),
    ("E from bulk and shear is out by a factor of three", P,
     "        *youngsModulus = denominator != 0.0 ? 9.0 * bulk * shearModulus / denominator : 0.0;",
     "        *youngsModulus = denominator != 0.0 ? 3.0 * bulk * shearModulus / denominator : 0.0;",
     "kill"),
    ("Lame's lambda is written with the plane-stress denominator", P,
     "    const double lambda = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));",
     "    const double lambda = e * nu / ((1.0 + nu) * (1.0 - nu));", "kill"),
    ("the elastic matrix's diagonal gets mu rather than 2mu", P,
     "        c[i * kVoigt + i] += 2.0 * mu;",
     "        c[i * kVoigt + i] += mu;", "kill"),
    # DOC: "elasticStress was reachable by no test at all. Dropping its deviatoric
    # split -- returning K tr(eps) + 2 mu eps_ii instead of K tr(eps) + 2 mu(eps_ii
    # - tr/3) -- passed everything, on a public entry point. It is now tied both to
    # elasticModuli computed by a different route and, bit for bit, to the elastic
    # branch of the return map."
    ("DOC: elasticStress drops its deviatoric split", P,
     "        stress[i] = pressure + 2.0 * mu * (strain[i] - volumetric / 3.0);",
     "        stress[i] = pressure + 2.0 * mu * strain[i];", "kill"),
    ("elasticStress doubles its shear stresses", P,
     "        stress[3 + i] = mu * strain[3 + i];",
     "        stress[3 + i] = 2.0 * mu * strain[3 + i];", "kill"),
    ("the deviator takes the mean off the shear components too", P,
     "        out[3 + i] = stress[3 + i];\n    }\n}",
     "        out[3 + i] = stress[3 + i] - mean;\n    }\n}", "kill"),
    ("von Mises uses sqrt(2/3) where sqrt(3/2) belongs", P,
     "double vonMises(const double stress[kVoigt]) { return std::sqrt(1.5) * deviatorNorm(stress); }",
     "double vonMises(const double stress[kVoigt]) { return kRoot23 * deviatorNorm(stress); }",
     "kill"),
    ("triaxiality comes back with the wrong sign", P,
     "    return q > 0.0 ? meanStress(stress) / q : 0.0;",
     "    return q > 0.0 ? -meanStress(stress) / q : 0.0;", "kill"),
    ("the principal direction is the most compressive one", P,
     "        if (values[i] > values[best]) best = i;",
     "        if (values[i] < values[best]) best = i;", "kill"),
    ("the principal direction reads the eigenvector matrix transposed", P,
     "    for (int i = 0; i < 3; ++i) out[i] = vectors[i][best];",
     "    for (int i = 0; i < 3; ++i) out[i] = vectors[best][i];", "kill"),
    ("the eigenvector's arbitrary sign is fixed the other way", P,
     "    if (out[largest] < 0.0) norm = -norm;",
     "    if (out[largest] > 0.0) norm = -norm;", "kill"),

    # --- the return map ------------------------------------------------------------------
    # `increment` is default-constructed a few lines above and `failedNow` is false
    # there, so this is the early-out deleted -- a torn point re-enters the return
    # map and starts carrying stress again.
    ("a torn point starts carrying stress again", P,
     "    if (state.failed) {\n"
     "        for (int i = 0; i < kVoigt; ++i) stress[i] = 0.0;",
     "    if (state.failed && increment.failedNow) {\n"
     "        for (int i = 0; i < kVoigt; ++i) stress[i] = 0.0;", "kill"),
    ("the elastic strain adds the plastic strain instead of removing it", P,
     "    for (int i = 0; i < kVoigt; ++i) elastic[i] = strain[i] - state.plasticStrain[i];",
     "    for (int i = 0; i < kVoigt; ++i) elastic[i] = strain[i] + state.plasticStrain[i];",
     "kill"),
    ("the pressure is taken from the total strain rather than the elastic part", P,
     "    const double volumetric = elastic[0] + elastic[1] + elastic[2];\n"
     "    const double pressure = kappa * volumetric;",
     "    const double volumetric = strain[0] + strain[1] + strain[2];\n"
     "    const double pressure = kappa * volumetric;", "kill"),
    ("the yield surface is offset the wrong way by the back stress", P,
     "    for (int i = 0; i < kVoigt; ++i) relative[i] = trial[i] - state.backStress[i];",
     "    for (int i = 0; i < kVoigt; ++i) relative[i] = trial[i] + state.backStress[i];",
     "kill"),
    ("the yield radius is the flow stress rather than sqrt(2/3) of it", P,
     "    const double yieldRadius = kRoot23 * flowStress(material.flow, state.equivalentPlasticStrain);",
     "    const double yieldRadius = flowStress(material.flow, state.equivalentPlasticStrain);",
     "kill"),
    ("the yield test is inverted", P,
     "    const double excess = relativeNorm - yieldRadius;",
     "    const double excess = yieldRadius - relativeNorm;", "kill"),
    ("the elastic branch reports the invariant of the relative stress", P,
     "        increment.vonMises = vonMises(stress);\n"
     "        increment.triaxiality = triaxiality(stress);\n"
     "        return increment;",
     "        increment.vonMises = vonMises(relative);\n"
     "        increment.triaxiality = triaxiality(stress);\n"
     "        return increment;", "kill"),
    ("the consistency solve takes one iteration and commits to it", P,
     "    for (; iterations < 50; ++iterations) {",
     "    for (; iterations < 1; ++iterations) {", "kill"),
    ("the hardening is evaluated at the multiplier rather than the plastic strain", P,
     "        const double accumulated = state.equivalentPlasticStrain + kRoot23 * gamma;",
     "        const double accumulated = state.equivalentPlasticStrain + gamma;", "kill"),
    ("the Newton slope has the hardening term with the wrong sign", P,
     "        const double slope = -linearPart - (2.0 / 3.0) * flowSlope(material.flow, accumulated);",
     "        const double slope = -linearPart + (2.0 / 3.0) * flowSlope(material.flow, accumulated);",
     "kill"),
    ("DOC CONTROL: the Newton step tolerance is loosened from 1e-15 to 1e-6", P,
     "        if (std::abs(next - gamma) <= 1e-15 * gamma) {",
     "        if (std::abs(next - gamma) <= 1e-6 * gamma) {", "survive"),
    ("the flow direction is taken from the stress rather than the relative stress", P,
     "    for (int i = 0; i < kVoigt; ++i) direction[i] = relative[i] / relativeNorm;",
     "    for (int i = 0; i < kVoigt; ++i) direction[i] = trial[i] / relativeNorm;", "kill"),
    ("the plastic shear increment loses its engineering factor of two", P,
     "        plasticIncrement[3 + i] = 2.0 * gamma * direction[3 + i];",
     "        plasticIncrement[3 + i] = gamma * direction[3 + i];", "kill"),
    ("the return map moves away from the yield surface", P,
     "        stress[i] = pressure + trial[i] - 2.0 * mu * gamma * direction[i];",
     "        stress[i] = pressure + trial[i] + 2.0 * mu * gamma * direction[i];", "kill"),
    ("the back stress evolves without its two thirds", P,
     "        back[i] = state.backStress[i] + (2.0 / 3.0) * kinematic * gamma * direction[i];",
     "        back[i] = state.backStress[i] + kinematic * gamma * direction[i];", "kill"),
    ("plastic work is measured against the effective stress, like the dissipation", P,
     "    increment.plasticWork = dot(stress, plasticIncrement);",
     "    increment.plasticWork = dot(effective, plasticIncrement);", "kill"),
    ("the plastic strain is overwritten rather than accumulated", P,
     "        state.plasticStrain[i] += plasticIncrement[i];",
     "        state.plasticStrain[i] = plasticIncrement[i];", "kill"),
    ("damage is overwritten rather than accumulated", P,
     "        state.damage += increment.equivalentPlasticStrainIncrement / critical;",
     "        state.damage = increment.equivalentPlasticStrainIncrement / critical;", "kill"),
    ("triaxiality divides the failure strain rather than scaling it", P,
     "    const double critical = failureStrain * triaxialityFactor(material.failure,\n"
     "                                                              increment.triaxiality);",
     "    const double critical = failureStrain / triaxialityFactor(material.failure,\n"
     "                                                              increment.triaxiality);",
     "kill"),
    ("a point has to be twice as damaged before it tears", P,
     "    if (state.damage >= 1.0) {",
     "    if (state.damage >= 2.0) {", "kill"),
    ("the failure plane is computed after the stress has been zeroed", P,
     "        maxPrincipalDirection(stress, increment.failureNormal);\n"
     "        for (int i = 0; i < kVoigt; ++i) stress[i] = 0.0;",
     "        for (int i = 0; i < kVoigt; ++i) stress[i] = 0.0;\n"
     "        maxPrincipalDirection(stress, increment.failureNormal);", "kill"),
    ("the algorithmic tangent's theta hardens instead of softening", P,
     "    const double theta = 1.0 - 2.0 * mu * gamma / relativeNorm;",
     "    const double theta = 1.0 + 2.0 * mu * gamma / relativeNorm;", "kill"),
    ("the tangent's rank-one term stiffens the flow direction", P,
     "            tangent[i * kVoigt + j] -= 2.0 * mu * thetaBar * direction[i] * direction[j];",
     "            tangent[i * kVoigt + j] += 2.0 * mu * thetaBar * direction[i] * direction[j];",
     "kill"),
    ("the plastic modulus ignores the kinematic part", P,
     "        flowSlope(material.flow, state.equivalentPlasticStrain) + kinematic;",
     "        flowSlope(material.flow, state.equivalentPlasticStrain);", "kill"),

    # --- the element's plasticity hooks ---------------------------------------------------
    # DOC: "The element's size could be taken off its bottom face, because every
    # element under test was prismatic."
    ("DOC: the element's size is taken off its bottom face", S,
     "            mid[a][i] = 0.5 * (nodes[a * 3 + i] + nodes[(a + 4) * 3 + i]);",
     "            mid[a][i] = nodes[a * 3 + i];", "kill"),
    ("the in-plane size is the area rather than its square root", S,
     "    if (inPlane != nullptr) *inPlane = std::sqrt(area);",
     "    if (inPlane != nullptr) *inPlane = area;", "kill"),
    ("DOC: initialisePlasticState does not clear the history it initialises", S,
     "    state = ElementPlasticState{};",
     "    state.torn = false;", "kill"),
    ("the regularised failure strain is given its lengths the wrong way round", S,
     "        plasticity::regularisedFailureStrain(material.failure, inPlane, thickness);",
     "        plasticity::regularisedFailureStrain(material.failure, thickness, inPlane);", "kill"),
    # DOC: "An element could call itself torn on its first dead integration point.
    # Every tearing test until then strained the element uniformly."
    ("DOC: an element calls itself torn on its first dead integration point", S,
     "    state.torn = result.failedPoints == kGauss;",
     "    state.torn = result.failedPoints > 0;", "kill"),
    # DOC: "the pre-loop check that skips the enhanced modes for an already-degraded
    # element is output-identical to letting the in-step retry catch it -- it saves a
    # wasted Newton pass per step for the rest of a torn element's life, which is a
    # cost difference and not a behaviour one."
    # `easCount` is `forms.easCount` two lines above, so this is the pre-loop check
    # doing nothing -- the degraded element now reaches the retry instead.
    ("DOC CONTROL: the degraded element is caught by the retry rather than up front", S,
     "    if (degraded) easCount = 0;",
     "    if (degraded) easCount = forms.easCount;", "survive"),
    ("the element does not re-run the step in which a point died", S,
     "        if (!(diedThisStep && easCount > 0)) break;",
     "        if (!(diedThisStep && easCount < 0)) break;", "kill"),
    ("a degraded element keeps the enhanced parameters it had before it tore", S,
     "    for (int k = 0; k < kEas; ++k) state.enhanced[k] = alpha[k];",
     "    for (int k = 0; k < easCount; ++k) state.enhanced[k] = alpha[k];", "kill"),
    # **This one survives, and it is the catalogue's fault rather than the suite's.**
    # `state.point[gp]` is not written until after the Newton has finished, so it *is*
    # the start-of-step state and this substitution is exactly equivalent -- it does
    # not say what it means. The mutation intended here is
    # `if (iteration == 0) trial[gp] = start[gp];`, which restarts only once and lets
    # the history carry over afterwards. Left as it was run, because the published
    # sweep was run against it and a catalogue edited after the fact is a catalogue
    # that no longer matches the figure.
    ("the enhanced Newton updates the history in place", S,
     "            trial[gp] = start[gp];",
     "            trial[gp] = state.point[gp];", "kill"),
    ("the internal force comes back with the wrong sign", S,
     "            force[a * 3 + i] = -s;",
     "            force[a * 3 + i] = s;", "kill"),
    ("the co-rotated displacement uses R rather than R transpose", S,
     "            for (int k = 0; k < 3; ++k) s += r[i * 3 + k] * current[a * 3 + k];\n"
     "            u[a * 3 + i] = s - rest[a * 3 + i];",
     "            for (int k = 0; k < 3; ++k) s += r[k * 3 + i] * current[a * 3 + k];\n"
     "            u[a * 3 + i] = s - rest[a * 3 + i];", "kill"),
    ("the enhanced strain is dropped from the strain the law integrates", S,
     "                    for (int k = 0; k < easCount; ++k) s += forms.g[gp][i][k] * alpha[k];",
     "                    for (int k = 0; k < 0; ++k) s += forms.g[gp][i][k] * alpha[k];", "kill"),
    ("the enhanced tolerance loses the yield strength from its scale", S,
     "    const double scale = material.flow.yieldStrength * volume;",
     "    const double scale = volume;", "kill"),
    ("the eigenstrain is added to the mechanical strain rather than removed", S,
     "                    strain[i] = eigenstrain != nullptr ? s - eigenstrain[i] : s;",
     "                    strain[i] = eigenstrain != nullptr ? s + eigenstrain[i] : s;", "kill"),
    ("the dissipation is summed over the Gauss points without their weights", S,
     "        result.dissipation += forms.weight[gp] * increments[gp].dissipation;",
     "        result.dissipation += increments[gp].dissipation;", "kill"),
]
