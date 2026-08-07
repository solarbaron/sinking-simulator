// SPDX-License-Identifier: MIT
#include "constraint.hpp"

#include <algorithm>
#include <cmath>

namespace sim::constraint {
namespace {

// Two-point Gauss through a rectangle of the profile. The stations are the roots
// of the second Legendre polynomial mapped onto [lo, hi], which is what makes the
// area, the first moment *and* the second moment of the rectangle exact -- see the
// header. One station would give the first two and lose the third, and the third
// is `I_own` -- 27% of a 200x10 flat bar's stiffened panel about its own neutral
// axis, and a 22.8% error in the panel measured about the plate's mid-surface,
// which is exactly what the mutant that sets this to zero delivers.
constexpr double kGaussOffset = 0.28867513459481288225;  // 1 / (2 sqrt 3)

Vec3 nodeAt(const std::vector<double>& nodal, std::uint32_t node) {
    const std::size_t base = static_cast<std::size_t>(node) * 3;
    return Vec3{nodal[base], nodal[base + 1], nodal[base + 2]};
}

// The uniaxial return map, backward Euler on the flow rule and solved rather than
// stepped, exactly as `plasticity::update` is for the three-dimensional case. It
// is written here rather than reused from there because that law is strain-driven
// in three dimensions: getting a uniaxial *stress* state out of it means iterating
// on the two transverse strains, and the scalar problem below is the closed form
// that iteration would converge to.
//
// Isotropic hardening only. A fibre has no back stress, so a reversal reyields at
// `-sigma_y(eps_p)` rather than at `sigma - 2 sigma_y0`; that is the same choice
// `plasticity.hpp` defaults to and for the same reason -- there is no measurement
// that sets a kinematic modulus.
//
// `failureStrain` is the regularised value from `fiberFailureStrain`, passed rather
// than looked up because it belongs to the fibre and not to the material -- exactly
// as `plasticity::update` takes it. `plasticity::kNeverFails` exercises the flow rule
// alone.
double uniaxialReturn(const plasticity::Material& material, double failureStrain, double strain,
                      FiberState& state, double* dissipationPerVolume, bool* failedNow) {
    if (dissipationPerVolume) *dissipationPerVolume = 0.0;
    if (failedNow) *failedNow = false;
    // A torn fibre carries nothing and learns nothing. Checked before anything else
    // so that "unloading does not heal it" is structural rather than a consequence
    // of the arithmetic below -- the same shape, and for the same reason, as
    // `plasticity::update`'s first branch.
    if (state.failed) return 0.0;

    const double youngs = material.youngsModulus;
    const double trial = youngs * (strain - state.plasticStrain);
    const double yield = plasticity::flowStress(material.flow, state.equivalentPlasticStrain);
    if (std::abs(trial) <= yield) return trial;

    // First estimate is exact for a linear curve, so the loop below confirms it and
    // stops; for Swift it is the Newton start.
    const double slope0 = plasticity::flowSlope(material.flow, state.equivalentPlasticStrain);
    double increment = (std::abs(trial) - yield) / (youngs + slope0);
    for (int iteration = 0; iteration < 32; ++iteration) {
        const double residual =
            std::abs(trial) - youngs * increment -
            plasticity::flowStress(material.flow, state.equivalentPlasticStrain + increment);
        if (std::abs(residual) <= 1e-13 * (yield + std::abs(trial))) break;
        const double slope =
            youngs + plasticity::flowSlope(material.flow, state.equivalentPlasticStrain + increment);
        if (!(slope > 0)) break;
        increment += residual / slope;
        if (increment < 0) increment = 0;
    }

    const double direction = trial >= 0 ? 1.0 : -1.0;
    const double stress = trial - direction * youngs * increment;
    state.plasticStrain += direction * increment;
    state.equivalentPlasticStrain += increment;
    if (dissipationPerVolume) *dissipationPerVolume = std::abs(stress) * increment;

    // --- Damage ------------------------------------------------------------------
    //
    // Accumulated rather than compared, so a path that yields in tension, unloads and
    // yields again spends the right fraction of its life in each -- and so that
    // compression, which contributes no damage at all, does not un-do what tension
    // did. The same two lines as `plasticity::update`'s, with the bar's own closed
    // form for the triaxiality: exactly `referenceTriaxiality` in tension and exactly
    // `cutoffTriaxiality` in compression. See the header §2b.
    const double critical =
        failureStrain * plasticity::triaxialityFactor(material.failure, fiberTriaxiality(stress));
    if (critical > 0.0 && std::isfinite(critical)) state.damage += increment / critical;

    if (state.damage >= 1.0) {
        state.damage = 1.0;
        state.failed = true;
        if (failedNow) *failedNow = true;
        // The stress drops discontinuously, which an explicit scheme feels as a
        // small shock -- element deletion, on one bar.
        return 0.0;
    }
    return stress;
}

}  // namespace

double fiberTriaxiality(double axialStress) {
    if (axialStress == 0.0) return 0.0;
    // `+/-(1.0/3.0)` to the bit: `1.0/3.0` here rounds identically to the `1.0/3.0`
    // and `-1.0/3.0` the `Failure` defaults are built from, and negation is exact. So
    // tension lands exactly on `referenceTriaxiality` -- multiplier exactly 1 -- and
    // compression exactly on `cutoffTriaxiality`, whose branch is a `<=`.
    return (axialStress > 0.0 ? 1.0 : -1.0) / 3.0;
}

double fiberFailureStrain(const plasticity::Failure& failure, double restLength,
                          double neckWidth) {
    return plasticity::regularisedFailureStrain(failure, restLength, neckWidth);
}

// --- 1. The tie ----------------------------------------------------------------

double tieWeight(double offset, double thickness) {
    if (!(thickness > 0)) return 0.5;
    return (offset + 0.5 * thickness) / thickness;
}

Vec3 tiedPoint(const Tie& tie, const std::vector<double>& nodal) {
    const Vec3 bottom = nodeAt(nodal, tie.bottom);
    const Vec3 top = nodeAt(nodal, tie.top);
    return bottom * (1.0 - tie.weight) + top * tie.weight;
}

void applyTiedForce(const Tie& tie, const Vec3& force, std::vector<double>& nodal) {
    const std::size_t bottom = static_cast<std::size_t>(tie.bottom) * 3;
    const std::size_t top = static_cast<std::size_t>(tie.top) * 3;
    for (int k = 0; k < 3; ++k) {
        nodal[bottom + static_cast<std::size_t>(k)] += (1.0 - tie.weight) * force[k];
        nodal[top + static_cast<std::size_t>(k)] += tie.weight * force[k];
    }
}

// --- 2. The eccentric stiffener -------------------------------------------------

ProfileFibers profileFibers(const StiffenerProfile& profile, double plateThickness, double sign) {
    ProfileFibers fibers;
    const double webHeight = std::max(0.0, profile.webHeight);
    const double webThickness = std::max(0.0, profile.webThickness);
    const bool flanged = profile.kind != ProfileKind::FlatBar && profile.flangeWidth > 0 &&
                         profile.flangeThickness > 0;
    const double flangeWidth = flanged ? profile.flangeWidth : 0.0;
    const double flangeThickness = flanged ? profile.flangeThickness : 0.0;

    // Each rectangle of the profile, as (from, to, area), measured from the face
    // of the plate the web is welded to. An angle's flange is offset sideways from
    // the web and a tee's straddles it: that difference lives entirely in the
    // weak-axis second moment, which fibres on the stiffener line do not carry, so
    // the two decompose identically here. Recorded because it is a real omission,
    // not because it is free.
    //
    // `width` is the rectangle's own thickness -- the dimension a neck localises
    // across, and the fibre's half of the plating's `elementSize` pair. It is taken
    // from the profile here rather than derived from `area` and the station spacing
    // downstream, because a derivation would divide by a height that a degenerate
    // profile can make zero.
    struct Rect { double from, to, area, width; };
    Rect part[2];
    int count = 0;
    if (webHeight > 0 && webThickness > 0)
        part[count++] = Rect{0.0, webHeight, webHeight * webThickness, webThickness};
    if (flanged)
        part[count++] = Rect{webHeight, webHeight + flangeThickness, flangeWidth * flangeThickness,
                             flangeThickness};

    const double faceOffset = 0.5 * plateThickness;
    for (int i = 0; i < count; ++i) {
        const double centre = 0.5 * (part[i].from + part[i].to);
        const double extent = part[i].to - part[i].from;
        for (int side = -1; side <= 1; side += 2) {
            const double height = centre + side * kGaussOffset * extent;
            fibers.offset[fibers.count] = sign * (faceOffset + height);
            fibers.area[fibers.count] = 0.5 * part[i].area;
            fibers.width[fibers.count] = part[i].width;
            ++fibers.count;
        }
    }
    return fibers;
}

RestFibers restFibers(const Stiffening& stiffening, const std::vector<double>& rest) {
    RestFibers forms;
    forms.length.assign(stiffening.fiber.size(), 0.0);
    forms.ok = true;
    for (std::size_t i = 0; i < stiffening.fiber.size(); ++i) {
        const Fiber& fiber = stiffening.fiber[i];
        const double span = length(tiedPoint(fiber.end[1], rest) - tiedPoint(fiber.end[0], rest));
        forms.length[i] = span;
        if (!(span > 0)) forms.ok = false;
    }
    return forms;
}

FiberForces fiberForces(const Stiffening& stiffening, const RestFibers& forms,
                        const std::vector<double>& current, const plasticity::Material& material,
                        std::vector<FiberState>* state, std::vector<double>& force,
                        bool allowFailure) {
    FiberForces out;
    const double youngs = material.youngsModulus;
    if (state != nullptr && state->size() < stiffening.fiber.size())
        state->resize(stiffening.fiber.size());

    for (std::size_t i = 0; i < stiffening.fiber.size(); ++i) {
        const Fiber& fiber = stiffening.fiber[i];
        const double restLength = i < forms.length.size() ? forms.length[i] : 0.0;
        if (!(restLength > 0) || !(fiber.area > 0)) continue;

        const Vec3 first = tiedPoint(fiber.end[0], current);
        const Vec3 second = tiedPoint(fiber.end[1], current);
        const Vec3 span = second - first;
        const double span_ = length(span);
        if (!(span_ > 0)) continue;
        const Vec3 axis = span / span_;
        // Nominal strain from the change in length. Exactly zero for any rigid body
        // motion, finite rotations included, which is the same statement the
        // co-rotational element makes and is why a fibre on a swinging patch does
        // not invent force.
        const double strain = (span_ - restLength) / restLength;

        double stress = 0;
        if (state != nullptr) {
            double dissipationPerVolume = 0;
            bool failedNow = false;
            const double before = (*state)[i].equivalentPlasticStrain;
            // The failure strain is the fibre's own -- its rest length as the
            // averaging length, its rectangle's thickness as the neck width. Both
            // are properties of *this* fibre, so it is computed here rather than
            // once for the set: a run over a curved patch has fibres of different
            // lengths, and a member with a flange has two neck widths in it.
            const double failureStrain =
                allowFailure ? fiberFailureStrain(material.failure, restLength, fiber.neckWidth)
                             : plasticity::kNeverFails;
            stress = uniaxialReturn(material, failureStrain, strain, (*state)[i],
                                    &dissipationPerVolume, &failedNow);
            out.dissipation += dissipationPerVolume * fiber.area * restLength;
            if ((*state)[i].equivalentPlasticStrain > before) ++out.yielded;
            if (failedNow) ++out.tornNow;
            if ((*state)[i].failed) {
                ++out.torn;
                out.tornVolume += fiber.area * restLength;
            }
        } else {
            stress = youngs * strain;
        }

        const double axial = stress * fiber.area;
        // The force the fibre applies to its own ends: a fibre in tension pulls
        // them together. Same sign convention as `solidshell::internalForce`.
        applyTiedForce(fiber.end[0], axis * axial, force);
        applyTiedForce(fiber.end[1], axis * -axial, force);
        out.strainEnergy += 0.5 * stress * stress / youngs * fiber.area * restLength;
    }
    return out;
}

void lumpFiberMass(const Stiffening& stiffening, const RestFibers& forms, double density,
                   std::vector<double>& nodalMass) {
    for (std::size_t i = 0; i < stiffening.fiber.size(); ++i) {
        const Fiber& fiber = stiffening.fiber[i];
        const double restLength = i < forms.length.size() ? forms.length[i] : 0.0;
        if (!(restLength > 0) || !(fiber.area > 0)) continue;
        // Half the fibre to each end, then that half split equally over the pair.
        // Not the consistent condensation: see the header -- `T^T M T` puts a
        // negative mass on one of the two nodes for any eccentric tie, and an
        // explicit scheme cannot integrate that.
        const double quarter = 0.25 * density * fiber.area * restLength;
        for (const Tie& tie : fiber.end) {
            nodalMass[tie.bottom] += quarter;
            nodalMass[tie.top] += quarter;
        }
    }
}

// --- 3. Stiffness ---------------------------------------------------------------

FiberStiffness fiberStiffness(const Fiber& fiber, const std::vector<double>& rest,
                              double restLength, double youngsModulus) {
    FiberStiffness out;
    if (!(restLength > 0) || !(fiber.area > 0)) return out;
    const Vec3 span = tiedPoint(fiber.end[1], rest) - tiedPoint(fiber.end[0], rest);
    const double span_ = length(span);
    if (!(span_ > 0)) return out;
    const Vec3 axis = span / span_;

    out.scale = youngsModulus * fiber.area / restLength;
    int slot = 0;
    for (int end = 0; end < 2; ++end) {
        const Tie& tie = fiber.end[end];
        const double sign = end == 0 ? -1.0 : 1.0;
        const std::uint32_t node[2] = {tie.bottom, tie.top};
        const double share[2] = {1.0 - tie.weight, tie.weight};
        for (int which = 0; which < 2; ++which)
            for (int k = 0; k < 3; ++k) {
                out.dof[slot] = node[which] * 3u + static_cast<std::uint32_t>(k);
                out.vector[slot] = sign * share[which] * axis[k];
                ++slot;
            }
    }
    return out;
}

double fiberFrequencySquared(const FiberStiffness& stiffness,
                             const std::vector<double>& nodalMass) {
    if (!(stiffness.scale > 0)) return 0.0;
    double sum = 0;
    for (int i = 0; i < 12; ++i) {
        const std::size_t node = stiffness.dof[i] / 3u;
        if (node >= nodalMass.size() || !(nodalMass[node] > 0)) continue;
        sum += stiffness.vector[i] * stiffness.vector[i] / nodalMass[node];
    }
    return stiffness.scale * sum;
}

AttachedForms attachedForms(const Stiffening& stiffening, const std::vector<double>& rest,
                            const RestFibers& forms, double youngsModulus,
                            const std::vector<FiberState>* state) {
    AttachedForms out;
    out.stiffness.reserve(stiffening.fiber.size());
    out.stress.reserve(stiffening.fiber.size());
    for (std::size_t i = 0; i < stiffening.fiber.size(); ++i) {
        const double restLength = i < forms.length.size() ? forms.length[i] : 0.0;
        const bool gone = state != nullptr && i < state->size() && (*state)[i].failed;
        const FiberStiffness fiber =
            fiberStiffness(stiffening.fiber[i], rest, restLength, youngsModulus);
        // One skip test for both halves. Splitting this loop in two would let the
        // two lists disagree about which fibres they dropped, and then every
        // stress after the first disagreement would belong to a different member
        // from the block beside it -- a plausible number for the wrong fibre.
        // A torn fibre is one more clause on that same test, never a second one.
        if (gone || !(fiber.scale > 0)) continue;
        solidshell::DofBlock block;
        block.dof.assign(fiber.dof, fiber.dof + 12);
        block.stiffness.assign(144, 0.0);
        for (int p = 0; p < 12; ++p)
            for (int q = 0; q < 12; ++q)
                block.stiffness[static_cast<std::size_t>(p) * 12 + static_cast<std::size_t>(q)] =
                    fiber.scale * fiber.vector[p] * fiber.vector[q];
        out.stiffness.push_back(std::move(block));

        // sigma = E * (v . u) / L. `fiber.scale` is EA/L and the area cancels out
        // of the stress, so E/L is taken directly rather than as scale/area --
        // one fewer place a zero area could divide.
        std::vector<double> gradient(12, 0.0);
        const double perMetre = youngsModulus / restLength;
        for (int p = 0; p < 12; ++p) gradient[static_cast<std::size_t>(p)] =
            perMetre * fiber.vector[p];
        out.stress.push_back(std::move(gradient));
    }
    return out;
}

std::vector<solidshell::DofBlock> stiffnessBlocks(const Stiffening& stiffening,
                                                  const std::vector<double>& rest,
                                                  const RestFibers& forms, double youngsModulus,
                                                  const std::vector<FiberState>* state) {
    return attachedForms(stiffening, rest, forms, youngsModulus, state).stiffness;
}

std::vector<MemberFibers> memberDamage(const Stiffening& stiffening, const RestFibers& forms,
                                       const std::vector<FiberState>& state) {
    std::vector<MemberFibers> out;
    const auto slot = [&](int member) -> MemberFibers& {
        for (MemberFibers& entry : out)
            if (entry.member == member) return entry;
        out.push_back(MemberFibers{});
        out.back().member = member;
        return out.back();
    };

    for (std::size_t i = 0; i < stiffening.fiber.size(); ++i) {
        const Fiber& fiber = stiffening.fiber[i];
        const double restLength = i < forms.length.size() ? forms.length[i] : 0.0;
        // The same skip test `fiberForces` and `attachedForms` use, so a fibre that
        // carries nothing by geometry is not counted as a member that tore.
        if (!(restLength > 0) || !(fiber.area > 0)) continue;
        MemberFibers& entry = slot(fiber.member);
        const double volume = fiber.area * restLength;
        entry.volume += volume;
        ++entry.fibers;
        // A short state is "nothing failed" -- an elastic solve keeps no state at
        // all, and reading past its end as failure would report a torn ship.
        if (i < state.size() && state[i].failed) {
            ++entry.torn;
            entry.lost += volume;
        } else {
            entry.carrying += volume;
        }
    }

    for (MemberFibers& entry : out)
        entry.effectiveness = entry.volume > 0 ? entry.carrying / entry.volume : 1.0;
    std::sort(out.begin(), out.end(),
              [](const MemberFibers& a, const MemberFibers& b) { return a.member < b.member; });
    return out;
}

double criticalTimestep(const Stiffening& stiffening, const RestFibers& forms,
                        const std::vector<double>& rest, const std::vector<double>& nodalMass,
                        double youngsModulus, double safety) {
    double worst = 0;
    for (std::size_t i = 0; i < stiffening.fiber.size(); ++i) {
        const double restLength = i < forms.length.size() ? forms.length[i] : 0.0;
        const FiberStiffness stiffness =
            fiberStiffness(stiffening.fiber[i], rest, restLength, youngsModulus);
        worst = std::max(worst, fiberFrequencySquared(stiffness, nodalMass));
    }
    if (!(worst > 0)) return 0.0;
    return safety * 2.0 / std::sqrt(worst);
}

// --- 4. Building ----------------------------------------------------------------

std::size_t addStiffener(const SeamRun& run, const StiffenerProfile& profile,
                         double plateThickness, const std::vector<double>& rest,
                         Stiffening& out) {
    const std::size_t stations = std::min(run.bottom.size(), run.top.size());
    if (stations < 2) return 0;
    const ProfileFibers profileSet = profileFibers(profile, plateThickness, run.sign);
    if (profileSet.count == 0) return 0;

    // The mid-surface of the pair, which is the stiffener's own line. Taking the
    // segment off the bottom nodes instead would be a fraction of a thickness out
    // on flat plating and the plate's turning angle out on a curved one.
    const auto midpoint = [&](std::size_t station) {
        return (nodeAt(rest, run.bottom[station]) + nodeAt(rest, run.top[station])) * 0.5;
    };

    std::size_t added = 0;
    for (std::size_t s = 0; s + 1 < stations; ++s) {
        const double segment = length(midpoint(s + 1) - midpoint(s));
        if (!(segment > 0)) continue;
        out.length += segment;
        for (int j = 0; j < profileSet.count; ++j) {
            Fiber fiber;
            for (int end = 0; end < 2; ++end) {
                const std::size_t station = s + static_cast<std::size_t>(end);
                fiber.end[end].bottom = run.bottom[station];
                fiber.end[end].top = run.top[station];
                fiber.end[end].weight = tieWeight(profileSet.offset[j], plateThickness);
            }
            fiber.area = profileSet.area[j];
            fiber.offset = profileSet.offset[j];
            fiber.neckWidth = profileSet.width[j];
            fiber.member = run.member;
            // The fibre's own length, which on a curved patch is not the seam's:
            // the tie extrapolates away from the mid-surface, so a fibre on the
            // inside of a bend is shorter than the plating it is welded to. That
            // is the eccentricity doing its job and it is what carries the moment.
            const double span =
                length(tiedPoint(fiber.end[1], rest) - tiedPoint(fiber.end[0], rest));
            out.mass += out.material.density * fiber.area * span;
            out.fiber.push_back(fiber);
            ++added;
        }
    }
    return added;
}

}  // namespace sim::constraint
