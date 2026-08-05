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
double uniaxialReturn(const plasticity::Material& material, double strain, FiberState& state,
                      double* dissipationPerVolume) {
    const double youngs = material.youngsModulus;
    const double trial = youngs * (strain - state.plasticStrain);
    const double yield = plasticity::flowStress(material.flow, state.equivalentPlasticStrain);
    if (dissipationPerVolume) *dissipationPerVolume = 0.0;
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
    return stress;
}

}  // namespace

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
    struct Rect { double from, to, area; };
    Rect part[2];
    int count = 0;
    if (webHeight > 0 && webThickness > 0)
        part[count++] = Rect{0.0, webHeight, webHeight * webThickness};
    if (flanged)
        part[count++] = Rect{webHeight, webHeight + flangeThickness, flangeWidth * flangeThickness};

    const double faceOffset = 0.5 * plateThickness;
    for (int i = 0; i < count; ++i) {
        const double centre = 0.5 * (part[i].from + part[i].to);
        const double extent = part[i].to - part[i].from;
        for (int side = -1; side <= 1; side += 2) {
            const double height = centre + side * kGaussOffset * extent;
            fibers.offset[fibers.count] = sign * (faceOffset + height);
            fibers.area[fibers.count] = 0.5 * part[i].area;
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
                        std::vector<FiberState>* state, std::vector<double>& force) {
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
            const double before = (*state)[i].equivalentPlasticStrain;
            stress = uniaxialReturn(material, strain, (*state)[i], &dissipationPerVolume);
            out.dissipation += dissipationPerVolume * fiber.area * restLength;
            if ((*state)[i].equivalentPlasticStrain > before) ++out.yielded;
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

std::vector<solidshell::DofBlock> stiffnessBlocks(const Stiffening& stiffening,
                                                  const std::vector<double>& rest,
                                                  const RestFibers& forms, double youngsModulus) {
    std::vector<solidshell::DofBlock> blocks;
    blocks.reserve(stiffening.fiber.size());
    for (std::size_t i = 0; i < stiffening.fiber.size(); ++i) {
        const double restLength = i < forms.length.size() ? forms.length[i] : 0.0;
        const FiberStiffness fiber =
            fiberStiffness(stiffening.fiber[i], rest, restLength, youngsModulus);
        if (!(fiber.scale > 0)) continue;
        solidshell::DofBlock block;
        block.dof.assign(fiber.dof, fiber.dof + 12);
        block.stiffness.assign(144, 0.0);
        for (int p = 0; p < 12; ++p)
            for (int q = 0; q < 12; ++q)
                block.stiffness[static_cast<std::size_t>(p) * 12 + static_cast<std::size_t>(q)] =
                    fiber.scale * fiber.vector[p] * fiber.vector[q];
        blocks.push_back(std::move(block));
    }
    return blocks;
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
