// SPDX-License-Identifier: MIT
// What an impact does to a panel: how far it dents, how hard it resists, how much
// energy that costs, and when it tears.
//
// `collision` (when it lands) says two hulls met and how much energy went into
// the meeting. `plasticity` says a point of steel has torn. `breach` turns torn
// panels into flooding openings. **Nothing joins the first to the third**, and
// this is the cheapest honest thing that does.
//
// The model is rigid-plastic membrane stretching, which is the standard first
// answer in ship-collision work and is exact enough to be tested rather than
// calibrated. A plate spanning `L` between frames, struck over a width `w`, is
// pushed into a tent of depth `delta`. Bending resistance is neglected -- for
// ship plating at penetrations of interest it is one to two orders below the
// membrane term -- so the plate carries the load purely by stretching, at the
// fully plastic membrane tension `sigma_y t` per unit width.
//
// Every quantity then has a closed form and they are all exact:
//
//     strain      eps(delta) = sqrt(1 + (2 delta / L)^2) - 1
//     force       F(delta)   = 2 sigma_y t w delta / sqrt((L/2)^2 + delta^2)
//     energy      E(delta)   = 2 sigma_y t w (sqrt((L/2)^2 + delta^2) - L/2)
//     tears at    delta_f    = (L/2) sqrt((1 + eps_f)^2 - 1)
//
// The energy integral is the force integrated over the penetration and comes out
// in closed form, which is what makes the inverse -- *given* the energy an impact
// delivered, how far did it get in and did it tear -- exact rather than a search.
//
// **What this is not.** It is one panel, struck normally, by something rigid,
// with the plating already fully plastic and the boundaries held. It has no
// bending, no stiffeners resisting in their own right, no membrane relief from
// neighbouring bays pulling in, no strain-rate effect, and no oblique impact. Of
// those the missing stiffeners matter most: on a real side structure the
// longitudinals carry a comparable share and this therefore *under*-predicts the
// resistance, which means it over-predicts penetration for a given energy. That
// is the conservative direction for a damage-stability question and the wrong one
// for a survivability question, and either way it is the number to replace when
// the zone FEM can be driven over a real patch.
//
// **How badly, measured.** `tools/ram_view` drives a 5175 t hull into the 8984 t
// ferry abeam and reports the energy the contact model actually delivered, which
// is well short of the striker's kinetic energy because both ships keep moving:
//
//     2 m/s,   7 MJ  ->    5 bays,    8.6 m2
//     4 m/s,  31 MJ  ->   26 bays,   44.4 m2
//     6 m/s,  74 MJ  ->   63 bays,  107.7 m2
//     8 m/s, 135 MJ  ->  114 bays,  196.5 m2
//
// An earlier table here read 27 / 107 / 234 / 415 m2 and **both** of its columns
// were wrong. The energy was `0.5 m v^2` over the *struck* ship's 8984 t -- the
// wrong ship, and as though all of it reached the plating -- overstating what is
// delivered by 2.1 to 2.6 times. Nearly all of the difference is that.
//
// **The span, and what correcting it actually cost.** `impactDamage()` once took
// `span = frameSpacing`, 2.4 m, where a longitudinally framed side spans between
// *longitudinals* at 0.70 m. `zone.{hpp,cpp}` settled it from outside, having no
// span in it at all -- only plating and where it is held: a 2 m punch into the
// ferry's side resists at **18.9 MN at 0.078 m**, against 10.6 MN for this model
// on the short span and **1.10 MN** on the long one.
//
// It is fixed, and the natural reading of that finding is far too broad, so it is
// worth being exact. The energy to tear a bay is `sigma_y t A eps_f`, in which
// **the span cancels** -- it reaches the answer only through the failure-strain
// regularisation, which is nearly flat here. Per unit struck area a bay absorbs
// 0.646 MJ/m2 on the long span against 0.680 on the short, so the hole moves 5%:
// the old energies through the corrected span give 25.6 / 102.0 / 222.4 / 392.8 m2
// against the old 27 / 107 / 234 / 415. Not the order of magnitude an earlier
// draft of this comment claimed -- that reasoned from force to area, and they do
// not scale together. What *was* wrong by 3.4x is force and penetration: 0.686 m
// of denting against 0.205, in the direction of reporting the hull far softer than
// it is.
//
// The low end is credible. **The high end is not** -- 196 m2 is a seventh of her
// side, and a real 135 MJ collision does not open that, because the striking bow crushes
// too, the longitudinals resist in their own right, and friction takes a share.
// None of those are here, so every joule goes into tearing the struck plating and
// the hole grows linearly with energy when it should saturate. Treat the answer as
// an upper bound that is tight at low energy and loose at high, and do not quote a
// hole from it above a few tens of megajoules.
#pragma once

#include "plasticity.hpp"
#include "scantlings.hpp"

#include <string>
#include <vector>

namespace sim {

// The panel being struck. `span` is between supports -- *longitudinal* spacing for
// a longitudinally framed side, frame spacing for a transversely framed one -- and
// `contactWidth` is how much of it the striking body actually touches. Plating
// spans the short way, between whichever stiffeners are closer together; getting
// these two round the wrong way is the defect recorded above.
struct IndentedPanel {
    double span = 2.40;          // m, L
    double thickness = 0.012;    // m, t
    double contactWidth = 2.00;  // m, w
    double yieldStrength = 355.0e6;   // Pa
    // Strain at which the membrane tears. Comes from
    // `plasticity::regularisedFailureStrain()`, because it is a property of the
    // length over which the strain is smeared and not of the steel alone.
    double failureStrain = 0.20;
};

struct IndentationState {
    double penetration = 0;   // m
    double strain = 0;        // membrane strain, dimensionless
    double force = 0;         // N
    double energy = 0;        // J absorbed to reach this penetration
    bool torn = false;
};

// Membrane strain at a given penetration. Exact.
double membraneStrain(double span, double penetration);

// The penetration at which a given membrane strain is reached -- the inverse of
// the above, in closed form rather than by search.
double penetrationForStrain(double span, double strain);

// Resisting force and absorbed energy at a penetration. Both exact.
double indentationForce(const IndentedPanel& panel, double penetration);
double indentationEnergy(const IndentedPanel& panel, double penetration);

// The penetration that absorbs a given energy: the inverse of the energy
// expression, and closed form for the same reason. Returns the tearing
// penetration if the energy exceeds what the panel can absorb intact.
double penetrationForEnergy(const IndentedPanel& panel, double energy);

// Everything at once, including whether it has torn.
IndentationState indentAt(const IndentedPanel& panel, double penetration);

// Energy the panel can absorb before it tears. The quantity a collision has to
// exceed for anything to flood at all.
double energyToTear(const IndentedPanel& panel);

// --- Against a real ship --------------------------------------------------------

// Panels within `radius` of `impact` on the struck hull, and the penetration the
// given energy drives into them.
//
// Energy is spent outward, panel by panel: the nearest bay resists until it
// tears, then the next, until the strike runs out of energy and the last panel
// takes what is left without letting go. `radius` bounds how far damage may
// reach; it does not decide how big the hole is.
//
// The first version shared the energy over a fixed patch by area instead, and
// that was measurably the wrong mechanism. Each panel is capped at its own
// tearing penetration, so once the patch had torn there was nowhere for further
// energy to go: a 2 m/s strike and an 8 m/s strike tore the same sixteen panels
// and opened the same 27.3 m2. The hole was a property of the contact radius and
// not of the collision. A bow that has punched through keeps going.
struct ImpactDamage {
    std::vector<int> panels;       // indices into StructuralMesh::panels
    std::vector<int> torn;         // the subset that tore
    double penetration = 0;        // m, of the worst-struck panel
    double energyAbsorbed = 0;     // J, total over the patch
    double tornArea = 0;           // m^2, the hole that results
    // Energy the strike still had when it ran out of reachable panels. Non-zero
    // means `radius` truncated the answer, not the structure -- the hole would
    // have been bigger and the caller is looking at a bound rather than a result.
    // Worth checking, because a tight radius silently reinstates exactly the
    // defect the outward march was written to remove.
    double energyUnspent = 0;      // J
};

// `material` supplies the failure strain, which is regularised on each panel's
// own span and thickness rather than taken as one number.
ImpactDamage impactDamage(const StructuralMesh& structure, const Vec3& impact, double radius,
                          double energy, const Scantlings& scantlings,
                          const plasticity::Material& material = plasticity::shipSteel());

// Every way this model is being used outside what it can honestly claim.
std::vector<std::string> validateIndentation(const IndentedPanel& panel);

}  // namespace sim
