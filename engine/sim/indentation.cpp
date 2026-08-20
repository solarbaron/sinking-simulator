// SPDX-License-Identifier: MIT
#include "indentation.hpp"

#include <algorithm>
#include <cmath>

namespace sim {
namespace {

// Half the span, which every expression below is written about.
double half(double span) { return 0.5 * std::max(span, 1e-12); }

}  // namespace

// A dent is small compared with the bay it is in, and every closed form in this
// file is a difference of two nearly equal lengths when it is. Written the way the
// algebra reads, each one cancels: `sqrt(1 + r^2) - 1` at r = 1e-5 subtracts two
// numbers that agree to eleven digits and keeps the noise. Multiplying above and
// below by the conjugate telescopes the numerator and leaves the identical value
// with the subtraction gone -- `(sqrt(A) - B)(sqrt(A) + B) = A - B^2`, so
// `sqrt(A) - B = (A - B^2) / (sqrt(A) + B)` whenever `sqrt(A) + B > 0`, and here
// `A - B^2` is a bare square that no cancellation touches. These are the same
// numbers, not approximations of them; the measured cost of the direct forms is
// quoted at each one and swept in `tests/test_indentation.cpp`.

double membraneStrain(double span, double penetration) {
    if (!(span > 0)) return 0.0;
    const double ratio = penetration / half(span);
    // eps = sqrt(1 + r^2) - 1 = ((1 + r^2) - 1) / (sqrt(1 + r^2) + 1) = r^2 / (1 + sqrt(1 + r^2)).
    //
    // The direct form loses the answer into the 1: sqrt(1 + r^2) rounds to within
    // half an ulp of 1 while the answer is r^2/2, so the relative error is about
    // 2.2e-16 / r^2. On the 2.40 m reference span that is 2.0e-12 at a 12 mm dent,
    // 3.6e-9 at 0.12 mm, 8.3e-8 at 12 um and 2.3e-2 at 0.12 um -- no digits left at
    // all a little below a micron, where the true strain is still a clean 1e-14.
    return ratio * ratio / (1.0 + std::sqrt(1.0 + ratio * ratio));
}

double penetrationForStrain(double span, double strain) {
    if (!(span > 0) || strain <= 0) return 0.0;
    // Inverting eps = sqrt(1 + (d/h)^2) - 1 gives d = h sqrt((1 + eps)^2 - 1),
    // exactly. No search, and no tolerance to argue about.
    //
    // (1 + eps)^2 - 1 = eps^2 + 2 eps = eps (eps + 2), identically -- the 1 that
    // carries no information is never formed. It is the same cancellation as
    // above, in the other direction and just as expensive: forming `1.0 + strain`
    // rounds away everything below eps * 2.2e-16 / eps = 2.2e-16 absolute, so the
    // strain arrives with a relative error of 2.2e-16 / eps and the square root
    // halves it into the answer -- 3.0e-9 of the returned penetration at eps = 5e-9
    // (a 0.12 mm dent on the reference span), 1.1e-2 at eps = 5e-15.
    //
    // `strain > 0` above is what makes both roots real, so no clamp is needed.
    return half(span) * std::sqrt(strain) * std::sqrt(strain + 2.0);
}

double indentationForce(const IndentedPanel& p, double penetration) {
    if (!(p.span > 0) || !(p.thickness > 0) || penetration <= 0) return 0.0;
    const double h = half(p.span);
    const double slant = std::sqrt(h * h + penetration * penetration);
    // Two legs of the tent, each pulling at the fully plastic membrane tension
    // sigma_y * t * w, resolved into the direction the indenter moves.
    return 2.0 * p.yieldStrength * p.thickness * p.contactWidth * penetration / slant;
}

double indentationEnergy(const IndentedPanel& p, double penetration) {
    if (!(p.span > 0) || !(p.thickness > 0) || penetration <= 0) return 0.0;
    const double h = half(p.span);
    // The force integrates in closed form: d/dd of sqrt(h^2 + d^2) is exactly
    // d / sqrt(h^2 + d^2), which is the shape of the force. So the energy is
    // scale * (sqrt(h^2 + d^2) - h), and `penetrationForEnergy` inverts that by
    // solving rather than by searching.
    //
    // This comment used to end "that is what makes the inverse below exact too",
    // which is a statement about the algebra wearing the clothes of one about the
    // arithmetic. The closed form is exact; *evaluating* it as that subtraction is
    // not, and for a shallow dent it is not even close, because sqrt(h^2 + d^2) and
    // h agree to their first 2 log10(h/d) digits and the answer is what is left.
    // Measured on the reference bay (h = 1.20 m): the energy carries a relative
    // error of 2.0e-12 at a 12 mm dent, 4.3e-11 at 1.2 mm, 1.1e-8 at 0.12 mm and
    // 8.3e-8 at 12 um -- about nine of its sixteen digits gone at 12 um and all but
    // three at 0.12 um. Nothing about the closed form prevents that; only writing
    // it without the subtraction does.
    //
    // sqrt(h^2 + d^2) - h = ((h^2 + d^2) - h^2) / (sqrt(h^2 + d^2) + h)
    //                     = d^2 / (sqrt(h^2 + d^2) + h),
    // identically, and h > 0 so the denominator cannot vanish.
    return 2.0 * p.yieldStrength * p.thickness * p.contactWidth *
           (penetration * penetration / (std::sqrt(h * h + penetration * penetration) + h));
}

double penetrationForEnergy(const IndentedPanel& p, double energy) {
    if (!(p.span > 0) || !(p.thickness > 0) || energy <= 0) return 0.0;
    const double h = half(p.span);
    const double scale = 2.0 * p.yieldStrength * p.thickness * p.contactWidth;
    if (scale <= 0) return 0.0;
    // E = scale (sqrt(h^2 + d^2) - h)  =>  d = sqrt((E/scale + h)^2 - h^2).
    //
    // With u = E/scale the bracket is (u + h)^2 - h^2 = u^2 + 2 h u = u (u + 2h),
    // identically -- and that matters here more than anywhere else in the file,
    // because the direct form differences two O(h^2) numbers whose true difference
    // is only 2 h u. A 12 um dent on the reference bay has u = 6e-11 m against
    // h = 1.20 m, so (u + h)^2 and h^2 agree to ten digits and the answer is the
    // rounding; the round trip through `indentationEnergy` comes back 2.7e-7 out
    // there, and 8.8e-7 at the worst point in that decade, against about one ulp
    // for the form below.
    //
    // Split as two roots rather than sqrt(u * (u + 2h)): both factors are positive
    // (energy > 0 and scale > 0 give u > 0, and h > 0), it is the same value, and
    // neither factor can overflow when the other is large.
    const double u = energy / scale;
    const double penetration = std::sqrt(u) * std::sqrt(u + 2.0 * h);
    // Past tearing the model has nothing further to say: the membrane is gone and
    // what happens next is a different mechanism. Report the tearing penetration
    // rather than extrapolating a torn plate's resistance.
    const double tearing = penetrationForStrain(p.span, p.failureStrain);
    return tearing > 0 ? std::min(penetration, tearing) : penetration;
}

IndentationState indentAt(const IndentedPanel& p, double penetration) {
    IndentationState s;
    s.penetration = std::max(0.0, penetration);
    s.strain = membraneStrain(p.span, s.penetration);
    s.force = indentationForce(p, s.penetration);
    s.energy = indentationEnergy(p, s.penetration);
    s.torn = p.failureStrain > 0 && s.strain >= p.failureStrain;
    return s;
}

double energyToTear(const IndentedPanel& p) {
    return indentationEnergy(p, penetrationForStrain(p.span, p.failureStrain));
}

ImpactDamage impactDamage(const StructuralMesh& structure, const Vec3& impact, double radius,
                          double energy, const Scantlings& scantlings,
                          const plasticity::Material& material) {
    ImpactDamage damage;
    if (!(radius > 0) || !(energy > 0)) return damage;

    const double frameSpacing =
        structure.frameSpacing > 0 ? structure.frameSpacing : scantlings.frameSpacing;
    const double stiffenerSpacing =
        scantlings.longitudinalSpacing > 0 ? scantlings.longitudinalSpacing : frameSpacing;

    // Panels the strike can reach, nearest first.
    //
    // The first version shared the energy over a fixed patch by area, and
    // measurement showed that was the wrong mechanism: because each panel is
    // capped at its own tearing penetration, once the patch had torn there was
    // nowhere left for energy to go. A 2 m/s strike and an 8 m/s strike tore the
    // same sixteen panels and opened the same 27.3 m2 -- the hole was a property
    // of the contact radius and not of the collision at all.
    //
    // A bow that has punched through does not stop; it keeps going into the next
    // bay. So the energy is spent outward, panel by panel, and the hole grows
    // with the strike. `radius` is now a bound on how far damage may reach rather
    // than the thing that decides its size.
    std::vector<std::pair<double, int>> reachable;
    for (std::size_t i = 0; i < structure.panels.size(); ++i) {
        const PlatePanel& p = structure.panels[i];
        if (p.role != PanelRole::Shell) continue;
        const double distance = length(p.centroid() - impact);
        if (distance > radius) continue;
        reachable.push_back({distance, static_cast<int>(i)});
    }
    std::sort(reachable.begin(), reachable.end());
    if (reachable.empty()) return damage;

    double remaining = energy;
    for (const auto& [distance, index] : reachable) {
        (void)distance;
        const PlatePanel& panel = structure.panels[static_cast<std::size_t>(index)];

        // A plate spans the *short* way between its supports, so the span is the
        // smaller of the two spacings and not the frame spacing by convention.
        // This was authored as the frame spacing and it is wrong: on the ferry the
        // plating is bounded by longitudinals at 0.70 m, not by frames at 2.40 m.
        //
        // Found by the zone FEM, which has no span in it at all -- only plating
        // and where it is held -- so it settles the question from outside the
        // model. What it costs is worth being exact about, because the headline
        // reading of that finding is too broad: the energy to tear a bay moves
        // only 5%, because the failure strain is nearly flat over this range of
        // element size, so the *number* of bays torn barely changes. What moves is
        // penetration, 0.686 m to 0.205 m, and resisting force, 2.96 MN to
        // 10.35 MN -- both by a factor of 3.4, and both in the direction of the
        // hull being far stiffer and far less deeply dented than reported.
        const double span = std::min(frameSpacing, stiffenerSpacing);
        IndentedPanel model;
        model.span = span;
        model.thickness = panel.thickness;
        // The struck width is the panel's extent along the *other* direction,
        // which for a bay of known area is its area over the span.
        model.contactWidth = std::max(panel.area() / std::max(span, 1e-6), 1e-6);
        model.yieldStrength =
            panel.material >= 0 && panel.material < static_cast<int>(structure.materials.size())
                ? structure.materials[static_cast<std::size_t>(panel.material)].yieldStrength
                : ah36Steel().yieldStrength;
        // Regularised on this panel's own geometry, not taken as one number: the
        // failure strain is a property of the length the strain is smeared over.
        // Regularised on the span, which is the length the membrane strain is
        // smeared over -- not on the frame spacing, which is now the other side.
        model.failureStrain =
            plasticity::regularisedFailureStrain(material.failure, span, panel.thickness);

        const double toTear = energyToTear(model);
        damage.panels.push_back(index);

        if (remaining >= toTear && toTear > 0) {
            // Through it, and on to the next.
            remaining -= toTear;
            damage.energyAbsorbed += toTear;
            damage.torn.push_back(index);
            damage.tornArea += panel.area();
            damage.penetration = std::max(damage.penetration,
                                          penetrationForStrain(model.span, model.failureStrain));
        } else {
            // Stopped here: this panel takes what is left and holds.
            const double reached = penetrationForEnergy(model, remaining);
            damage.energyAbsorbed += indentationEnergy(model, reached);
            damage.penetration = std::max(damage.penetration, reached);
            remaining = 0.0;
            break;
        }
    }
    // Whatever is left ran out of hull to spend itself on, not out of energy.
    damage.energyUnspent = remaining;
    return damage;
}

std::vector<std::string> validateIndentation(const IndentedPanel& p) {
    std::vector<std::string> problems;
    if (!(p.span > 0)) problems.push_back("span is not positive");
    if (!(p.thickness > 0)) problems.push_back("thickness is not positive");
    if (!(p.contactWidth > 0)) problems.push_back("contact width is not positive");
    if (!(p.yieldStrength > 0)) problems.push_back("yield strength is not positive");
    if (p.failureStrain <= 0 || p.failureStrain > 1.0)
        problems.push_back("failure strain outside (0, 1]");
    if (p.thickness > 0 && p.span / p.thickness < 20.0)
        problems.push_back("span/thickness below 20 is a block, not a membrane -- bending "
                           "resistance is no longer negligible and this model omits it");
    if (p.contactWidth > 4.0 * p.span)
        problems.push_back("contact much wider than the span: the surrounding structure would "
                           "share the load and this model gives it no way to");
    return problems;
}

}  // namespace sim
