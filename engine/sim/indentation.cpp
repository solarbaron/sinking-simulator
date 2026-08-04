// SPDX-License-Identifier: MIT
#include "indentation.hpp"

#include <algorithm>
#include <cmath>

namespace sim {
namespace {

// Half the span, which every expression below is written about.
double half(double span) { return 0.5 * std::max(span, 1e-12); }

}  // namespace

double membraneStrain(double span, double penetration) {
    if (!(span > 0)) return 0.0;
    const double ratio = penetration / half(span);
    return std::sqrt(1.0 + ratio * ratio) - 1.0;
}

double penetrationForStrain(double span, double strain) {
    if (!(span > 0) || strain <= 0) return 0.0;
    // Inverting eps = sqrt(1 + (d/h)^2) - 1 gives d = h sqrt((1+eps)^2 - 1),
    // exactly. No search, and no tolerance to argue about.
    const double stretched = 1.0 + strain;
    return half(span) * std::sqrt(std::max(0.0, stretched * stretched - 1.0));
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
    // d / sqrt(h^2 + d^2), which is the shape of the force. That is what makes
    // the inverse below exact too.
    return 2.0 * p.yieldStrength * p.thickness * p.contactWidth *
           (std::sqrt(h * h + penetration * penetration) - h);
}

double penetrationForEnergy(const IndentedPanel& p, double energy) {
    if (!(p.span > 0) || !(p.thickness > 0) || energy <= 0) return 0.0;
    const double h = half(p.span);
    const double scale = 2.0 * p.yieldStrength * p.thickness * p.contactWidth;
    if (scale <= 0) return 0.0;
    // E = scale (sqrt(h^2 + d^2) - h)  =>  d = sqrt((E/scale + h)^2 - h^2).
    const double slant = energy / scale + h;
    const double penetration = std::sqrt(std::max(0.0, slant * slant - h * h));
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
