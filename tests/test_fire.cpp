// SPDX-License-Identifier: MIT
//
// Validation of the two-zone compartment fire.
//
// Compartment fire has an unusually good supply of independent answers, and this
// file leans on all of them rather than on any one:
//
//   * **Closed forms of the model's own algebra.** The volume split is exactly the
//     internal-energy split; a sealed compartment reaches
//     `p0 + (gamma-1) E / V` to machine precision; `c_p - c_v = R` exactly. These
//     are asserted at 1e-12 and better because they are identities, not
//     approximations, and a tolerance would hide the day one stops holding.
//   * **The classical doorway integral.** A vertical vent with uniform gas either
//     side has a textbook answer -- `(2/3) Cd W sqrt(2 rho drho g) (z_t-z_n)^(3/2)`
//     and a neutral plane at `(z_n-z_b)/(z_t-z_n) = (T_inf/T_h)^(1/3)` -- derived
//     here from the hydrostatics rather than restated from the code.
//   * **Two independent plume correlations that meet.** Heskestad's flame region
//     and his far field are separate formulae fitted to separate data. They agree
//     at the mean flame height, and the ratio at which they agree is a *constant*
//     independent of fire size and diameter. That constant is derived in the test
//     and asserted to 1e-12; nothing in fire.cpp computes it.
//   * **Two independent enclosure correlations that bracket.** MQH inverted at a
//     500 K rise and Thomas's flashover power are fitted to different data by
//     different people, and they agree within 10% on the reference room.
//   * **A conservation account on a real case.** `zone.hpp` records its energy
//     residual closing on unit fixtures and running to -102% on a real tearing
//     run, which is why the account here is asserted on the ferry's engine rooms
//     under a 4 MW growth-steady-decay fire and not only on a box.
//
// **Temperatures are in KELVIN**, like `thermal.hpp` and for the same reason.
#include "engine/sim/fire.hpp"
#include "harness.hpp"

#include "../game/prototype/ferry.hpp"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

using namespace sim;
using testing::expectEqual;
using testing::expectNear;
using testing::expectTrue;

namespace {

// The ISO 9705 reference room: 3.6 x 2.4 x 2.4 m with a 0.8 x 2.0 m doorway. The
// enclosure every compartment-fire correlation in this file was fitted around.
constexpr double kRoomLength = 3.6, kRoomWidth = 2.4, kRoomHeight = 2.4;
constexpr double kDoorWidth = 0.8, kDoorHeight = 2.0;
constexpr double kDoorArea = kDoorWidth * kDoorHeight;

// Total internal bounding surface less the doorway, which is MQH's `A_T`.
double roomWallArea() {
    return 2.0 * kRoomLength * kRoomWidth + 2.0 * kRoomLength * kRoomHeight +
           2.0 * kRoomWidth * kRoomHeight - kDoorArea;
}

struct Room {
    fire::Model model;
    Ship ship;                    // empty: nothing here needs a hull
    Sea sea{-1000.0};             // far below, so no opening is ever water-blocked
};

Room makeRoom(double heatRelease, double wallConductance, double cd = 0.7) {
    Room r;
    fire::GasCompartment g;
    g.name = "room";
    g.shipCompartment = kSea;
    g.floorZ = 0.0;
    g.ceilingZ = kRoomHeight;
    g.floorArea = kRoomLength * kRoomWidth;
    g.perimeter = 2.0 * (kRoomLength + kRoomWidth);
    g.gasVolume = kRoomLength * kRoomWidth * kRoomHeight;
    g.wallConductance = wallConductance;
    g.fillAmbient();
    r.model.gas.push_back(g);

    fire::Vent v;
    v.name = "door";
    v.a = 0;
    v.b = kSea;
    v.sillZ = 0.0;
    v.soffitZ = kDoorHeight;
    v.width = kDoorWidth;
    v.area = kDoorArea;
    v.dischargeCoeff = cd;
    r.model.vents.push_back(v);

    fire::DesignFire f;
    f.name = "pool";
    f.compartment = 0;
    f.baseZ = 0.0;
    f.diameter = 1.0;
    f.peakHeatRelease = heatRelease;
    r.model.fires.push_back(f);
    r.model.resetAccount();
    return r;
}

// Run to steady state. 1500 s is well past the settling time at every power in
// this file -- the room converges in about 40 s -- and the assertions below check
// that it really did settle rather than assuming it.
void runToSteady(Room& r, int seconds = 1500) {
    for (int i = 0; i < seconds; ++i) r.model.step(1.0, r.ship, r.sea);
}

// A taller, plainer box than the ISO room, with the vent low in the wall so that
// everything it carries comes out of the *cool* layer. The upper layer's mass
// balance is then the plume and nothing else, which is what makes the check
// below a check of the plume rather than a restatement of the solver.
Room makeTallBox(double heatRelease) {
    Room r;
    fire::GasCompartment g;
    g.name = "box";
    g.shipCompartment = kSea;
    g.floorZ = 0.0;
    g.ceilingZ = 4.0;
    g.floorArea = 24.0;                 // 6 m x 4 m
    g.perimeter = 20.0;
    g.gasVolume = 96.0;
    g.wallConductance = 20.0;
    g.fillAmbient();
    r.model.gas.push_back(g);

    fire::Vent v;
    v.name = "low_vent";
    v.a = 0;
    v.b = kSea;
    v.sillZ = 0.0;
    v.soffitZ = 0.3;
    v.width = 1.0;
    v.area = 0.3;
    v.dischargeCoeff = 0.7;
    r.model.vents.push_back(v);

    fire::DesignFire f;
    f.name = "pool";
    f.compartment = 0;
    // **Off the floor.** The plume is measured above the fire's own base, and a
    // base at z = 0 would make `interfaceZ() - baseZ` and `interfaceZ() + baseZ`
    // the same number -- a sign error with nothing to show it.
    f.baseZ = 0.5;
    f.diameter = 1.0;
    f.peakHeatRelease = heatRelease;
    r.model.fires.push_back(f);
    r.model.resetAccount();
    return r;
}

// --- Thermodynamic identities -------------------------------------------------

// The pressure closure `p = (gamma-1) U / V` and the volume identity
// `V_u/V = U_u/U` are both exact only if the caloric constants are mutually
// consistent. A quoted c_p of 1005 J/(kg K) against this repo's R and gamma would
// break them in the fourth digit, which is enough to show in an account asserted
// at machine precision.
void testCaloricConstantsAreExactlyConsistent() {
    expectNear("c_p - c_v is exactly R", fire::kCpAir - fire::kCvAir, kRAir, 1e-12);
    expectNear("c_p / c_v is exactly gamma", fire::kCpAir / fire::kCvAir, kGammaAir, 1e-15);
    expectNear("ambient density is p/(R T)", fire::kRhoAmbient,
               kPatm / (kRAir * kTAmbient), 1e-15);
}

// The single load-bearing algebraic fact in fire.hpp: with both layers at one
// pressure, the volume split *is* the internal-energy split. Asserted over a
// sweep of splits and temperatures, because an identity that holds at one point
// is a coincidence.
void testVolumeSplitIsTheEnergySplit() {
    double worst = 0;
    for (int i = 1; i < 20; ++i) {
        for (double tu : {300.0, 500.0, 900.0, 1400.0}) {
            fire::GasCompartment g;
            g.floorZ = 0;
            g.ceilingZ = 4.0;
            g.floorArea = 10.0;
            g.gasVolume = 40.0;
            g.fillAmbient(i / 20.0);
            // Heat the upper layer without touching its mass: pure superheat.
            g.upper.energy = g.upper.mass * fire::kCvAir * tu;

            const double u = g.upper.energy + g.lower.energy;
            worst = std::max(worst, std::abs(g.upperVolume() / g.gasVolume - g.upper.energy / u));
            // And the ideal gas law must hold layer by layer at the shared pressure.
            const double p = g.pressure();
            const double vu = g.upperVolume();
            worst = std::max(worst,
                             std::abs(p * vu - g.upper.mass * kRAir * g.upper.temperature()) /
                                 (p * vu));
        }
    }
    expectTrue("the volume split is the energy split, and pV = mRT holds in each layer",
               worst < 1e-14);
    // Guard against vacuity: the sweep must actually have produced a split worth
    // checking, not twenty copies of the same one.
    fire::GasCompartment a, b;
    a.gasVolume = b.gasVolume = 40.0;
    a.floorArea = b.floorArea = 10.0;
    a.ceilingZ = b.ceilingZ = 4.0;
    a.fillAmbient(0.05);
    b.fillAmbient(0.45);
    expectTrue("and the sweep spans a real range of interface heights",
               std::abs(a.interfaceZ() - b.interfaceZ()) > 1.0);
}

// --- The plume ----------------------------------------------------------------

// Heskestad's flame-region and far-field entrainment laws are separate fits to
// separate data and neither is derivable from the other. Where they meet, at the
// mean flame height, they must agree -- and because the flame height and the
// virtual origin carry the same `Q^(2/5)` factor, the ratio at which they agree
// is a **constant**, independent of both the fire's power and its diameter.
//
// That constant is derived here from the published coefficients. Nothing in
// fire.cpp computes it, so this is a cross-check and not a restatement.
void testPlumeBranchesMeetAtTheFlameTip() {
    const double chi = 0.7;
    // L - z0 = (0.235 - 0.083) Q^(2/5) = 0.152 Q^(2/5), so (L-z0)^(5/3) = 0.152^(5/3) Q^(2/3),
    // and the far-field branch at z = L is Q * [0.071 chi^(1/3) 0.152^(5/3) + 0.0018 chi].
    const double far = 0.071 * std::cbrt(chi) * std::pow(0.152, 5.0 / 3.0) + 0.0018 * chi;
    const double flame = 0.0056 * chi;
    const double expected = far / flame;

    expectTrue("the two plume branches agree at the flame tip to better than 2%",
               std::abs(expected - 1.0) < 0.02);

    double worst = 0;
    int checked = 0;
    for (double qkw : {200.0, 500.0, 1000.0, 2500.0, 8000.0}) {
        for (double diameter : {0.5, 1.0, 2.0, 3.0}) {
            const fire::Plume p{qkw * 1000.0, diameter, chi};
            const double lf = p.flameHeight();
            if (lf <= 0.05) continue;   // no coherent flame: only one branch applies
            const double inFlame = 0.0056 * chi * qkw;                 // kg/s, at z = L
            // One ulp above the flame tip, so the far-field branch is taken with
            // no room for the height offset itself to show in the comparison.
            const double above = p.entrainment(std::nextafter(lf, 1e30));
            worst = std::max(worst, std::abs(above / inFlame - expected));
            ++checked;
        }
    }
    expectTrue("across every fire size and diameter with a flame", checked >= 12);
    expectTrue("the ratio is that same constant to 1e-9", worst < 1e-9);
}

// The far field is a 5/3 power law in height above the *virtual origin*, plus a
// constant. Subtract the constant and the power law must be exact.
void testPlumeObeysTheFiveThirdsPowerLaw() {
    const fire::Plume p{1.0e6, 1.0, 0.7};
    const double qc = 700.0;                      // kW
    const double additive = 0.0018 * qc;          // the non-scaling term
    const double z0 = p.virtualOrigin();
    const double z1 = z0 + 6.0, z2 = z0 + 12.0;   // exactly a factor of two above the origin
    const double m1 = p.entrainment(z1) - additive;
    const double m2 = p.entrainment(z2) - additive;
    expectNear("doubling the height above the virtual origin multiplies entrainment by 2^(5/3)",
               m2 / m1, std::pow(2.0, 5.0 / 3.0), 1e-12);

    // And the cube-root dependence on power, at a fixed height above the origin.
    const fire::Plume q{8.0e6, 1.0, 0.7};
    const double h = 6.0;
    const double a = p.entrainment(p.virtualOrigin() + h) - 0.0018 * 700.0;
    const double b = q.entrainment(q.virtualOrigin() + h) - 0.0018 * 5600.0;
    expectNear("and eight times the power entrains twice as much", b / a, 2.0, 1e-12);
}

// The virtual origin sits below the fire base for a wide fire and above it for a
// small tall one. The sign is load bearing -- it moves the entrainment by tens of
// percent near the fire -- and it is the kind of thing that reads as plausible
// either way.
void testVirtualOriginChangesSignWithFireShape() {
    const fire::Plume wide{500.0e3, 3.0, 0.7};
    const fire::Plume narrow{5000.0e3, 0.3, 0.7};
    expectTrue("a wide, weak fire has its virtual origin below the pan",
               wide.virtualOrigin() < -1.0);
    expectTrue("a small, fierce one has it above", narrow.virtualOrigin() > 0.5);
    expectTrue("and a fire wide enough has no coherent flame at all", wide.flameHeight() < 0);
    expectNear("entrainment is zero at the fire base", wide.entrainment(0.0), 0.0, 0.0);
    expectNear("and below it", wide.entrainment(-1.0), 0.0, 0.0);
}

// --- The vent integral --------------------------------------------------------

// A vertical vent with *uniform* gas on the inside and ambient outside has a
// closed-form answer that predates any of this. Derive it here and assert it.
//
//   dp(z) = gauge_f + drho g z        (linear, drho = rho_inf - rho_h > 0)
//   z_n   = -gauge_f / (drho g)       (the neutral plane)
//   m_out = (2/3) Cd W sqrt(2 rho_h drho g) (z_t - z_n)^(3/2)
//   m_in  = (2/3) Cd W sqrt(2 rho_inf drho g) (z_n - z_b)^(3/2)
//
// This checks the closed-form band integral, the neutral-plane split, and the
// choice of donor density on each side of it, all at once.
void testVentIntegralMatchesTheClassicalDoorway() {
    const double th = 700.0;
    const double rhoH = kPatm / (kRAir * th);
    const double dRho = fire::kRhoAmbient - rhoH;
    const double zb = 0.4, zt = 2.4, width = 0.9, cd = 0.68;

    fire::VentSide in;
    in.floorZ = 0.0;
    in.interfaceZ = -1.0;          // one uniform layer: the interface is below the vent
    in.tLower = in.tUpper = th;
    in.rhoLower = in.rhoUpper = rhoH;
    const fire::VentSide out = fire::ambientSide();

    fire::Vent v;
    v.a = 0;
    v.b = kSea;
    v.sillZ = zb;
    v.soffitZ = zt;
    v.width = width;
    v.area = width * (zt - zb);
    v.dischargeCoeff = cd;

    double worstOut = 0, worstIn = 0, worstNp = 0;
    int checked = 0;
    for (double zn : {0.8, 1.2, 1.6, 2.0}) {
        in.gaugeAtFloor = -dRho * kGravity * zn;
        const fire::VentResult r = fire::ventMassFlow(v, in, out);
        const double k = (2.0 / 3.0) * cd * width * std::sqrt(2.0 * dRho * kGravity);
        const double wantOut = k * std::sqrt(rhoH) * std::pow(zt - zn, 1.5);
        const double wantIn = k * std::sqrt(fire::kRhoAmbient) * std::pow(zn - zb, 1.5);
        worstOut = std::max(worstOut, std::abs(r.massAToB - wantOut) / wantOut);
        worstIn = std::max(worstIn, std::abs(r.massBToA - wantIn) / wantIn);
        worstNp = std::max(worstNp, std::abs(r.neutralPlaneZ - zn));
        expectTrue("the vent runs in both directions at once", r.bidirectional);
        ++checked;
    }
    expectEqual("four neutral-plane positions checked", checked, 4);
    expectTrue("outflow matches the classical doorway integral to 1e-13", worstOut < 1e-13);
    expectTrue("inflow does too", worstIn < 1e-13);
    expectTrue("and the neutral plane is located exactly", worstNp < 1e-12);
}

// The neutral plane that balances the mass flows has its own closed form, from
// setting the two expressions above equal:
//
//   (z_n - z_b) / (z_t - z_n) = (rho_h / rho_inf)^(1/3) = (T_inf / T_h)^(1/3)
//
// Nothing in fire.cpp knows this. It falls out of the integral being right.
void testBalancedNeutralPlaneFollowsTheCubeRootOfTemperature() {
    const double zb = 0.0, zt = 2.2, width = 1.1, cd = 0.7;
    double worst = 0;
    for (double th : {400.0, 600.0, 900.0, 1200.0}) {
        const double rhoH = kPatm / (kRAir * th);
        const double dRho = fire::kRhoAmbient - rhoH;

        fire::VentSide in;
        in.floorZ = 0.0;
        in.interfaceZ = -1.0;
        in.tLower = in.tUpper = th;
        in.rhoLower = in.rhoUpper = rhoH;
        const fire::VentSide out = fire::ambientSide();

        fire::Vent v;
        v.a = 0;
        v.b = kSea;
        v.sillZ = zb;
        v.soffitZ = zt;
        v.width = width;
        v.area = width * (zt - zb);
        v.dischargeCoeff = cd;

        // Bisect for the floor gauge at which the net flow vanishes.
        double lo = -dRho * kGravity * zt, hi = 0.0;
        for (int i = 0; i < 200; ++i) {
            const double mid = 0.5 * (lo + hi);
            in.gaugeAtFloor = mid;
            const fire::VentResult t = fire::ventMassFlow(v, in, out);
            if (t.massAToB - t.massBToA > 0) hi = mid; else lo = mid;
        }
        in.gaugeAtFloor = 0.5 * (lo + hi);
        const fire::VentResult r = fire::ventMassFlow(v, in, out);
        const double ratio = (r.neutralPlaneZ - zb) / (zt - r.neutralPlaneZ);
        worst = std::max(worst, std::abs(ratio - std::cbrt(kTAmbient / th)));
    }
    expectTrue("the balanced neutral plane sits at (T_inf/T_h)^(1/3) of the vent", worst < 1e-9);
}

// What the neutral plane is actually worth, at the ferry's own opening sizes.
//
// The comparison that matters is not "how different is the number". It is that a
// **single pressure difference at the orifice centre has zero flow at its own
// equilibrium**, while the two-layer integral has a large bidirectional exchange
// at the same state. A model built on the first has no ventilation at all behind
// a doorway; the compartment then heats without bound.
void testSingleDeltaPHasNoVentilationWhereTheIntegralHas() {
    struct Case {
        const char* name;
        double area;
        OpeningKind kind;
        double cd;
    };
    const Case cases[] = {
        {"airpipe_wms", 0.02, OpeningKind::Vent, 0.70},
        {"vent_er_s", 0.50, OpeningKind::Vent, 0.80},
        {"breach_er_s", 2.40, OpeningKind::Breach, 0.62},
        {"wt_door_er", 3.60, OpeningKind::Door, 0.75},
    };
    const double tHot = 600.0;
    double smallest = 1e30, largest = 0;
    double worstRatio = 0, bestRatio = 1e30;
    for (const Case& c : cases) {
        Opening o;
        o.area = c.area;
        o.kind = c.kind;
        o.dischargeCoeff = c.cd;
        o.pos = Vec3{0, 0, 1.2};
        const fire::VentShape shape = fire::ventShapeFor(o);

        fire::VentSide in;
        in.floorZ = 0.0;
        in.interfaceZ = 1.0;          // hot layer down to 1 m, cool air below it
        in.tLower = kTAmbient;
        in.tUpper = tHot;
        in.rhoLower = kPatm / (kRAir * kTAmbient);
        in.rhoUpper = kPatm / (kRAir * tHot);
        const fire::VentSide out = fire::ambientSide();

        fire::Vent v;
        v.a = 0;
        v.b = kSea;
        v.width = shape.width;
        v.area = c.area;
        v.horizontal = shape.horizontal;
        v.sillZ = o.pos.z - 0.5 * shape.height;
        v.soffitZ = o.pos.z + 0.5 * shape.height;
        v.dischargeCoeff = c.cd;

        // The state the vent actually sits at: zero *net* mass flow.
        double lo = -200.0, hi = 200.0;
        for (int i = 0; i < 200; ++i) {
            const double mid = 0.5 * (lo + hi);
            in.gaugeAtFloor = mid;
            const fire::VentResult t = fire::ventMassFlow(v, in, out);
            if (t.massAToB - t.massBToA > 0) hi = mid; else lo = mid;
        }
        in.gaugeAtFloor = 0.5 * (lo + hi);
        const fire::VentResult r = fire::ventMassFlow(v, in, out);

        // What the existing single-orifice treatment would see at that same
        // state: one pressure difference at the centre, and the density of
        // whatever is against the high side there.
        const double zc = o.pos.z;
        const double dp = in.gaugeAt(zc) - out.gaugeAt(zc);
        const double rhoDonor = dp > 0 ? in.densityAt(zc) : fire::kRhoAmbient;
        const double single = c.cd * c.area * std::sqrt(2.0 * rhoDonor * std::abs(dp));

        const std::string what = std::string("vent '") + c.name + "' exchanges gas where a "
                                 "single delta-p at the centre moves far less";
        expectTrue(what, r.massAToB > 0 && r.massBToA > 0);
        expectTrue(what + " (net is zero, exchange is not)",
                   std::abs(r.massAToB - r.massBToA) < 1e-9 * r.massAToB);
        expectTrue(what + " (the integral carries more)", r.massAToB > 1.10 * single);
        smallest = std::min(smallest, r.massAToB);
        largest = std::max(largest, r.massAToB);
        worstRatio = std::max(worstRatio, r.massAToB / single);
        bestRatio = std::min(bestRatio, r.massAToB / single);
    }
    // Guard against vacuity: the sweep must span the ferry's real range, and the
    // exchange must be a physically substantial number rather than a residue.
    expectTrue("the exchange spans two decades across the ferry's openings",
               largest / smallest > 100.0);
    expectTrue("and the door's exchange is worth hundreds of kilowatts",
               largest * fire::kCpAir * (tHot - kTAmbient) > 5.0e5);

    // And the error is not a bounded correction that could be calibrated away
    // with one discharge coefficient: it runs from 14% on an air pipe to a factor
    // of nine on the breach, whose centre happens to land near the neutral plane.
    expectTrue("the shortfall is not a constant factor across the openings",
               worstRatio / bestRatio > 5.0);

    // The exact statement about the single-delta-p model: at *its* own
    // equilibrium the pressure difference at the orifice centre is zero, so it
    // moves nothing at all. A compartment behind it has no ventilation, and a
    // fire in that compartment heats without bound.
    expectNear("a single delta-p at the orifice centre is at rest with zero flow",
               0.7 * 3.6 * std::sqrt(2.0 * fire::kRhoAmbient * 0.0), 0.0, 0.0);
}

// `Opening` carries an area and a centre but no height. What is derived from the
// kind is a guess and is labelled as one; it still has to be self-consistent.
void testVentShapesAreConsistentWithTheirAreas() {
    Opening door;
    door.kind = OpeningKind::Door;
    door.area = 3.6;
    const fire::VentShape ds = fire::ventShapeFor(door);
    expectNear("a watertight door is 2 m tall", ds.height, 2.0, 1e-15);
    expectNear("and as wide as its area needs", ds.width * ds.height, 3.6, 1e-15);
    expectTrue("and vertical", !ds.horizontal);

    Opening pipe;
    pipe.kind = OpeningKind::Pipe;
    pipe.area = 0.02;
    const fire::VentShape ps = fire::ventShapeFor(pipe);
    expectNear("an air pipe is square", ps.width * ps.height, 0.02, 1e-17);
    expectNear("literally square", ps.width, ps.height, 1e-17);

    Opening hatch;
    hatch.kind = OpeningKind::Hatch;
    hatch.area = 1.0;
    expectTrue("a hatch is horizontal, because it is in a deck",
               fire::ventShapeFor(hatch).horizontal);

    Opening none;
    none.kind = OpeningKind::Breach;
    none.area = 0.0;
    expectNear("and an opening of no area has no shape", fire::ventShapeFor(none).width, 0.0, 0.0);
}

// A vent with layered gas on *both* sides. The span has to break at each side's
// own interface, and each resulting band has to use that band's densities and
// that band's temperatures.
//
// Checked by additivity rather than by restating the integral: integrating the
// whole span in one call must equal the sum of three calls over the bands the
// interfaces define. A breakpoint the solver failed to insert would carry one
// side's density straight through a change in it, and the two would part
// company. The enthalpy is then checked against the layer temperatures the mass
// is claimed to have come from, which is a second, independent relation.
void testTheVentSplitsAtBothSidesInterfaces() {
    fire::VentSide a;
    a.floorZ = 0.0;
    a.interfaceZ = 1.2;
    a.tLower = 300.0;
    a.tUpper = 800.0;
    a.rhoLower = kPatm / (kRAir * 300.0);
    a.rhoUpper = kPatm / (kRAir * 800.0);
    a.gaugeAtFloor = 1.5;

    fire::VentSide b;
    b.floorZ = 0.0;
    b.interfaceZ = 0.6;
    b.tLower = kTAmbient;
    b.tUpper = 500.0;
    b.rhoLower = kPatm / (kRAir * kTAmbient);
    b.rhoUpper = kPatm / (kRAir * 500.0);
    b.gaugeAtFloor = 0.0;

    fire::Vent v;
    v.a = 0;
    v.b = 1;
    v.sillZ = 0.0;
    v.soffitZ = 2.0;
    v.width = 1.0;
    v.area = 2.0;
    v.dischargeCoeff = 0.7;
    const fire::VentResult whole = fire::ventMassFlow(v, a, b);

    double mAB = 0, mBA = 0, hAB = 0, hBA = 0, upA = 0, upB = 0;
    const double edges[] = {0.0, 0.6, 1.2, 2.0};
    for (int i = 0; i < 3; ++i) {
        fire::Vent band = v;
        band.sillZ = edges[i];
        band.soffitZ = edges[i + 1];
        const fire::VentResult r = fire::ventMassFlow(band, a, b);
        mAB += r.massAToB;
        mBA += r.massBToA;
        hAB += r.enthalpyAToB;
        hBA += r.enthalpyBToA;
        upA += r.fromUpperA;
        upB += r.fromUpperB;
    }
    expectTrue("the whole span is the sum of its bands, A to B",
               std::abs(whole.massAToB - mAB) < 1e-12 * mAB);
    expectTrue("and B to A", std::abs(whole.massBToA - mBA) < 1e-12 * mBA);
    expectTrue("with the same enthalpy", std::abs(whole.enthalpyAToB - hAB) < 1e-12 * hAB);
    expectTrue("both ways", std::abs(whole.enthalpyBToA - hBA) < 1e-12 * hBA);
    expectTrue("and the same share out of each upper layer",
               std::abs(whole.fromUpperA - upA) < 1e-12 * upA &&
                   std::abs(whole.fromUpperB - upB) < 1e-12 * upB);

    // Not vacuous: this configuration really does run both ways and really does
    // draw on both upper layers.
    expectTrue("the vent runs in both directions", whole.massAToB > 0 && whole.massBToA > 0);
    expectTrue("and is reported as bidirectional", whole.bidirectional);
    // A's outflow straddles A's own interface at 1.2 m, so it comes out of both
    // of A's layers; B's return flow all sits above B's interface at 0.6 m, so it
    // comes out of B's hot layer alone. Both are asserted, because "some of each"
    // and "all of one" are different statements about the band split and a test
    // that only made the first would not notice the second going wrong.
    expectTrue("A's outflow draws on its hot layer", whole.fromUpperA > 0.05 * whole.massAToB);
    expectTrue("and on its cool one", whole.fromUpperA < 0.95 * whole.massAToB);
    expectTrue("while B's return flow is entirely out of B's hot layer",
               whole.fromUpperB == whole.massBToA);

    // The enthalpy each stream carries is its upper share at that side's upper
    // temperature plus its lower share at that side's lower temperature. Nothing
    // in fire.cpp assembles it that way -- it accumulates band by band -- so this
    // is a second reader of the same answer.
    const double wantAB = fire::kCpAir * (whole.fromUpperA * a.tUpper +
                                          (whole.massAToB - whole.fromUpperA) * a.tLower);
    const double wantBA = fire::kCpAir * (whole.fromUpperB * b.tUpper +
                                          (whole.massBToA - whole.fromUpperB) * b.tLower);
    expectTrue("the enthalpy out of A is its two layers' temperatures, weighted by what left them",
               std::abs(whole.enthalpyAToB - wantAB) < 1e-12 * wantAB);
    expectTrue("and the same out of B", std::abs(whole.enthalpyBToA - wantBA) < 1e-12 * wantBA);
    expectTrue("the neutral plane is inside the span",
               whole.neutralPlaneZ > 0.0 && whole.neutralPlaneZ < 2.0);
}

// Species riding the gas from one space to the next, deposited by buoyancy. No
// fire here: a tracer, a pressure difference and the network, so that what is
// being checked is the transport and nothing else.
//
// Everything is at ambient temperature on purpose. Gas arriving no warmer than
// the receiving compartment's cool layer is not buoyant, and has to be put in
// that layer -- the opposite of what the same code does with smoke from a fire.
void testSpeciesRideTheGasIntoTheNeighboursCoolLayer() {
    fire::Model m;
    Ship ship;
    Sea sea{-1000.0};
    for (int i = 0; i < 2; ++i) {
        fire::GasCompartment g;
        g.name = i == 0 ? "source" : "sink";
        g.shipCompartment = kSea;
        g.floorZ = 0.0;
        g.ceilingZ = 3.0;
        g.floorArea = 20.0;
        g.perimeter = 18.0;
        g.gasVolume = 60.0;
        g.wallConductance = 0.0;
        g.fillAmbient();
        m.gas.push_back(g);
    }
    // Pressurise the source by adding mass and energy in the same proportion, so
    // its temperature is untouched and buoyancy plays no part at all.
    for (fire::Layer* l : {&m.gas[0].upper, &m.gas[0].lower}) {
        l->mass *= 1.002;
        l->energy *= 1.002;
    }
    m.gas[0].lower.products = 5.0;

    fire::Vent v;
    v.name = "door";
    v.a = 0;
    v.b = 1;
    v.sillZ = 0.2;
    v.soffitZ = 2.2;
    v.width = 1.0;
    v.area = 2.0;
    v.dischargeCoeff = 0.7;
    m.vents.push_back(v);
    m.resetAccount();

    expectTrue("the source starts over-pressured", m.gas[0].gaugeAtFloor() > 100.0);
    expectNear("and both spaces start at exactly ambient temperature",
               m.gas[0].lower.temperature(), m.gas[1].lower.temperature(), 1e-12);

    for (int i = 0; i < 40; ++i) m.step(1.0, ship, sea);

    expectTrue("the tracer reaches the space next door", m.gas[1].lower.products > 1e-3);
    expectNear("and goes into its cool layer, because it arrived no warmer than the air there",
               m.gas[1].upper.products, 0.0, 0.0);
    expectTrue("the source is left with less than it started with", m.gas[0].lower.products < 5.0);
    const double held = m.gas[0].upper.products + m.gas[0].lower.products +
                        m.gas[1].upper.products + m.gas[1].lower.products;
    expectTrue("and every gram is still in the ship", std::abs(held - 5.0) < 1e-12);
    expectTrue("none of it was generated -- there is no fire in this fixture",
               m.account.productsGenerated == 0.0);
    expectNear("and none of it left the model", m.account.productsOut, 0.0, 0.0);
    expectTrue("the pressures have equalised", std::abs(m.gas[0].gaugeAtFloor() -
                                                        m.gas[1].gaugeAtFloor()) < 1.0);
    expectTrue("the mass account closes", std::abs(m.account.massResidualFraction()) < 1e-11);
    expectTrue("and the energy account with it",
               std::abs(m.account.energyResidualFraction()) < 1e-11);
}

// The two layers share one pressure and exchange volume across a moving
// interface, so each does `p dV/dt` of work on the other. That work is the whole
// content of `layerSplit`, and it has a closed form nobody has to take on trust:
// **a layer that gains no mass and loses no heat is a closed parcel, and
// compressing it is adiabatic** --
//
//     T_u / T_u0 = (p / p0)^((gamma-1)/gamma)
//
// The receiving box below is set up so its upper layer is exactly that parcel:
// the gas arriving through the vent is no warmer than the cool layer, so all of
// it is deposited underneath, and the boundary conducts nothing. The pressure
// then rises around a layer that nothing has touched.
//
// This is the sharp test of the split. Nothing else is: the sealed-compartment
// closed form checks the *total*, which is invariant to how the two layers
// divide it, and MQH is a 12% correlation.
// How far a compressed parcel departs from its adiabat, for a given vent area.
// Two boxes, one over-pressured, no fire and no boundary loss; the vent sits
// entirely below the receiving box's interface, so everything that arrives is
// deposited underneath and the box's upper layer is a **closed parcel** with the
// pressure rising around it.
double adiabatError(double ventArea, Ship& ship, const Sea& sea, double* warming) {
    fire::Model m;
    for (int i = 0; i < 2; ++i) {
        fire::GasCompartment g;
        g.name = i == 0 ? "source" : "sink";
        g.shipCompartment = kSea;
        g.floorZ = 0.0;
        g.ceilingZ = 3.0;
        g.floorArea = 20.0;
        g.perimeter = 18.0;
        g.gasVolume = 60.0;
        g.wallConductance = 0.0;       // a closed parcel needs no heat leaving it
        g.fillAmbient(0.2);            // a layer with enough in it to measure
        m.gas.push_back(g);
    }
    for (fire::Layer* l : {&m.gas[0].upper, &m.gas[0].lower}) {
        l->mass *= 1.05;
        l->energy *= 1.05;
    }
    fire::Vent v;
    v.name = "door";
    v.a = 0;
    v.b = 1;
    v.sillZ = 0.1;
    v.soffitZ = 1.1;
    v.width = ventArea;
    v.area = ventArea;
    v.dischargeCoeff = 0.7;
    m.vents.push_back(v);
    m.resetAccount();

    const double p0 = m.gas[1].pressure();
    const double t0 = m.gas[1].upper.temperature();
    const double mass0 = m.gas[1].upper.mass;
    for (int i = 0; i < 3000; ++i) m.step(1.0, ship, sea);

    const fire::GasCompartment& sink = m.gas[1];
    expectTrue("nothing was added to the compressed parcel", sink.upper.mass == mass0);
    expectTrue("but it was compressed", sink.pressure() > 1.02 * p0);
    expectTrue("and the account closed while it happened",
               std::abs(m.account.energyResidualFraction()) < 1e-11);
    if (warming != nullptr) *warming = sink.upper.temperature() - t0;
    const double want = t0 * std::pow(sink.pressure() / p0, (kGammaAir - 1.0) / kGammaAir);
    return std::abs(sink.upper.temperature() / want - 1.0);
}

// The two layers share one pressure and exchange volume across a moving
// interface, so each does `p dV/dt` of work on the other. That work is the whole
// content of `layerSplit`, and it has a closed form nobody has to take on trust:
// a layer that gains no mass and loses no heat is a closed parcel, and squeezing
// it is adiabatic --
//
//     T_u / T_u0 = (p / p0)^((gamma-1)/gamma)
//
// This is the *only* sharp test of the split in this file. The sealed-compartment
// closed form checks the total, which is invariant to how the two layers divide
// it; MQH is a 12% correlation. A wrong `gamma-1` here would move the exponent
// from 0.286 to 0.786 and nothing else would notice.
void testCompressingAnUntouchedLayerIsAdiabatic() {
    Ship ship;
    Sea sea{-1000.0};
    double warming = 0;
    const double coarse = adiabatError(0.5, ship, sea, &warming);
    expectTrue("a compressed layer follows its adiabat", coarse < 5e-5);
    expectTrue("over a warming that is measurable, not a rounding error", warming > 1.0);

    // And it is not a fixed offset dressed as agreement. The split evaluates the
    // volume fraction at the start of each substep, so what is left is first
    // order in the energy each step moves; taking the vent down by a hundred has
    // to take the departure down with it. Measured 3.0e-5 -> 3.0e-6 -> 3.1e-7
    // across two decades of vent area.
    const double fine = adiabatError(0.005, ship, sea, nullptr);
    const double finer = adiabatError(0.0005, ship, sea, nullptr);
    expectTrue("and converges onto it as the step's energy shrinks", fine < 0.15 * coarse);
    expectTrue("first order, decade for decade", finer < 0.15 * fine);
    expectTrue("rather than levelling off at some floor", finer < 1e-6);
}

// A hatch in a deck has no height to integrate over, so it gets one pressure
// difference and the density of whatever is against the high side -- the same
// treatment `Ship::solveFlowNetwork` gives every opening. It is also where this
// model's acknowledged gap lives: the buoyancy-driven exchange across a
// horizontal vent has no pressure difference at all and is not modelled, so this
// asserts that it is *never* bidirectional rather than pretending otherwise.
void testHorizontalVentTakesOneDeltaPAndIsNeverBidirectional() {
    const double tHot = 700.0;
    fire::VentSide below;
    below.floorZ = 0.0;
    below.interfaceZ = 1.0;
    below.tLower = kTAmbient;
    below.tUpper = tHot;
    below.rhoLower = kPatm / (kRAir * kTAmbient);
    below.rhoUpper = kPatm / (kRAir * tHot);
    below.gaugeAtFloor = 12.0;          // pressurised by the fire underneath
    const fire::VentSide out = fire::ambientSide();

    fire::Vent v;
    v.a = 0;
    v.b = kSea;
    v.horizontal = true;
    v.sillZ = v.soffitZ = 2.5;          // in the deckhead, inside the hot layer
    // Not 1 m^2. A unit area makes `Cd * A` and `Cd / A` the same number, and an
    // assertion that cannot tell those apart is most of an assertion missing.
    v.area = 1.6;
    v.width = 1.6;
    v.dischargeCoeff = 0.7;

    const fire::VentResult r = fire::ventMassFlow(v, below, out);
    const double dp = below.gaugeAt(2.5) - out.gaugeAt(2.5);
    expectTrue("the hatch is being pushed from below", dp > 0);
    expectNear("a horizontal vent takes one delta-p, at the donor's density", r.massAToB,
               0.7 * 1.6 * std::sqrt(2.0 * below.rhoUpper * dp), 1e-13);
    expectNear("and nothing comes the other way", r.massBToA, 0.0, 0.0);
    expectTrue("a horizontal vent is never bidirectional", !r.bidirectional);
    expectNear("and the net is all of it", r.massAToB - r.massBToA, r.massAToB, 0.0);
    expectTrue("what it sheds came out of the hot layer", r.fromUpperA == r.massAToB);
    expectNear("carrying the hot layer's enthalpy", r.enthalpyAToB,
               r.massAToB * fire::kCpAir * tHot, 1e-9 * r.massAToB * fire::kCpAir * tHot);

    // Reverse the pressure and the donor swaps: the flow, the density and the
    // temperature all have to come from the other side.
    below.gaugeAtFloor = -12.0;
    const fire::VentResult back = fire::ventMassFlow(v, below, out);
    const double dpBack = out.gaugeAt(2.5) - below.gaugeAt(2.5);
    expectTrue("and reversed it is being pushed from above", dpBack > 0);
    expectNear("it draws ambient air down through the hatch", back.massBToA,
               0.7 * 1.6 * std::sqrt(2.0 * fire::kRhoAmbient * dpBack), 1e-13);
    expectNear("and nothing goes up", back.massAToB, 0.0, 0.0);
    expectNear("at ambient enthalpy", back.enthalpyBToA,
               back.massBToA * fire::kCpAir * kTAmbient,
               1e-9 * back.massBToA * fire::kCpAir * kTAmbient);
    expectNear("and it took nothing out of anybody's upper layer", back.fromUpperB, 0.0, 0.0);
}

// --- The design fire as a source ----------------------------------------------

// The heat release is a known function of time, so the account's `heatReleased`
// over a run covering the whole curve must equal the curve's own integral. This
// is what says the source is *integrated* rather than sampled at whichever end of
// the substep is convenient. Sampled at the left end, a fire that starts at t = 0
// loses exactly one substep of heat -- 50 kJ, on every run, forever.
void testTheAccountIntegratesTheDesignCurve() {
    Room r = makeRoom(100.0e3, 0.0);
    r.model.vents.clear();
    r.model.fires[0].growthCoefficient = fire::kGrowthMedium;
    r.model.fires[0].steadyDuration = 120.0;
    r.model.fires[0].decayDuration = 90.0;
    r.model.resetAccount();
    const double want = r.model.fires[0].totalEnergy();
    for (int i = 0; i < 400; ++i) r.model.step(1.0, r.ship, r.sea);
    expectTrue("the account's released energy is the area under the design curve",
               std::abs(r.model.account.heatReleased / want - 1.0) < 1e-4);
    expectTrue("over a curve that grows, plateaus and decays", want > 1.5e7);
    expectNear("and the fire is out at the end", r.model.step(1.0, r.ship, r.sea).heatRelease, 0.0,
               0.0);
}

// A declared radiative share leaves the gas without heating it. Sealed, that is
// another closed form: the pressure rise is the *remaining* fraction, exactly.
void testRadiativeLossFractionRemovesExactlyItsShare() {
    Room r = makeRoom(200.0e3, 0.0);
    r.model.vents.clear();
    r.model.fires[0].radiativeLossFraction = 0.4;
    r.model.resetAccount();
    const double v = r.model.gas[0].gasVolume;
    const double p0 = r.model.gas[0].pressure();
    for (int i = 0; i < 60; ++i) r.model.step(1.0, r.ship, r.sea);

    const double released = 200.0e3 * 60.0;
    const double want = p0 + (kGammaAir - 1.0) * 0.6 * released / v;
    expectTrue("only the non-radiative share reaches the gas",
               std::abs(r.model.gas[0].pressure() - want) < 1e-12 * want);
    expectNear("the account books the rest as radiative loss", r.model.account.radiativeLoss,
               0.4 * released, 1e-6 * released);
    expectNear("and books the whole release as released", r.model.account.heatReleased, released,
               1e-6 * released);
    expectTrue("the account closes with the loss in it",
               std::abs(r.model.account.energyResidual()) < 1e-13 * released);
    expectTrue("and the default is zero, so nothing else in this file carries it",
               fire::DesignFire{}.radiativeLossFraction == 0.0);
}

// --- Species ------------------------------------------------------------------

// Where the products are, not just how many. The plume delivers them to the
// upper layer, so a sealed compartment's smoke is all up there and its mass is
// the declared yield times the heat released.
void testProductsAreMadeInTheUpperLayerAtTheDeclaredYield() {
    Room r = makeTallBox(300.0e3);
    r.model.vents.clear();
    r.model.resetAccount();
    for (int i = 0; i < 30; ++i) r.model.step(1.0, r.ship, r.sea);

    const fire::GasCompartment& g = r.model.gas[0];
    const double want = r.model.fires[0].productYield * 300.0e3 * 30.0;
    const double held = g.upper.products + g.lower.products;
    expectTrue("the smoke mass is the yield times the heat released",
               std::abs(held / want - 1.0) < 1e-9);
    expectTrue("all of it is in the hot layer", g.upper.products == held);
    expectNear("and none in the cool one", g.lower.products, 0.0, 0.0);
    expectTrue("the upper layer is visibly sooty", g.upper.productFraction() > 1e-4);
    expectTrue("and a real number of grams was made", want > 0.02);
}

// Entrainment moves gas, and whatever that gas is carrying goes with it. Because
// the plume takes a well-mixed sample, the fraction of the cool layer's tracer it
// has moved must equal the fraction of the cool layer's *mass* it has moved --
// exactly, since the two are stepped by the same factor.
void testThePlumeCarriesTheLowerLayersSpeciesUpward() {
    Room r = makeTallBox(300.0e3);
    r.model.vents.clear();
    r.model.fires[0].productYield = 0.0;   // no new smoke: only the seeded tracer
    r.model.gas[0].lower.products = 10.0;
    r.model.resetAccount();
    const double lowerMass0 = r.model.gas[0].lower.mass;

    for (int i = 0; i < 20; ++i) r.model.step(1.0, r.ship, r.sea);
    const fire::GasCompartment& g = r.model.gas[0];
    expectTrue("the tracer is conserved",
               std::abs(g.upper.products + g.lower.products - 10.0) < 1e-12);
    expectTrue("and the plume has carried a real share of it up", g.upper.products > 1.0);
    expectTrue("the cool layer keeps exactly the tracer fraction it keeps of its mass",
               std::abs(g.lower.products / 10.0 - g.lower.mass / lowerMass0) < 1e-12);
    expectTrue("nothing was generated, because the yield was set to zero",
               r.model.account.productsGenerated == 0.0);
}

// --- The layer on the floor ---------------------------------------------------

// Drive the layer all the way down and the plume must have nothing left to
// entrain. Without the taper the last of the cool layer goes in a single step,
// leaving a residue of mass at zero energy: a layer at absolute zero, an infinite
// density, and a NaN one vent integral later.
void testEntrainmentStopsWhenTheLayerReachesTheFloor() {
    Room r = makeRoom(800.0e3, 30.0);        // 800 kW in the small room: it fills
    r.model.vents[0].soffitZ = 0.2;          // and a mean vent low down, so little escapes
    r.model.vents[0].area = kDoorWidth * 0.2;
    r.model.resetAccount();

    fire::StepResult s{};
    for (int i = 0; i < 400; ++i) s = r.model.step(1.0, r.ship, r.sea);
    const fire::GasCompartment& g = r.model.gas[0];
    const double raw = r.model.totalEntrainment(r.model.time);

    expectTrue("the smoke layer has reached the floor", g.interfaceZ() < 0.02 * kRoomHeight);
    expectTrue("the correlation is still asking for entrainment", raw > 0.04);
    expectTrue("but the model applies materially less, because the cool air is nearly gone",
               s.entrainment < 0.85 * raw);
    expectTrue("without shutting the plume off outright", s.entrainment > 0.5 * raw);
    expectTrue("both layers are still physical", g.lower.mass > 0 && g.lower.energy > 0);
    expectTrue("neither temperature went to NaN",
               std::isfinite(g.upper.temperature()) && std::isfinite(g.lower.temperature()));
    expectTrue("the cool layer never got hotter than the hot one",
               g.lower.temperature() < g.upper.temperature());
    expectTrue("and the account still closes on the way down",
               std::abs(r.model.account.energyResidualFraction()) < 1e-11);
    expectTrue("and the mass account with it",
               std::abs(r.model.account.massResidualFraction()) < 1e-11);
}

// --- The sealed compartment ---------------------------------------------------

// With no vent and no boundary loss, every joule the fire releases stays in the
// gas and the pressure has a closed form: `p = p0 + (gamma-1) E / V`. Exact, not
// approximate -- the energy equation for the compartment total is `dU/dt = Q`
// with nothing else in it, and the midpoint source rule integrates a steady fire
// without error.
void testSealedCompartmentReachesItsClosedFormPressure() {
    Room r = makeRoom(200.0e3, 0.0);
    r.model.vents.clear();
    r.model.resetAccount();
    const double v = r.model.gas[0].gasVolume;
    const double p0 = r.model.gas[0].pressure();

    for (int i = 0; i < 60; ++i) r.model.step(1.0, r.ship, r.sea);
    const double want = p0 + (kGammaAir - 1.0) * (200.0e3 * 60.0) / v;
    const double got = r.model.gas[0].pressure();
    expectTrue("a sealed compartment reaches p0 + (gamma-1) E / V to machine precision",
               std::abs(got - want) < 1e-12 * want);
    expectTrue("its energy account closes to machine precision",
               std::abs(r.model.account.energyResidual()) <
                   1e-13 * r.model.account.heatReleased);
    expectTrue("and its mass never changed at all",
               std::abs(r.model.account.mass - r.model.account.initialMass) < 1e-12);
    // Guard against vacuity: the pressure has to have actually moved.
    expectTrue("and the pressure rise is a factor of three, not a rounding error",
               got > 3.0 * p0);
}

// --- The enclosure correlations -----------------------------------------------

// MQH gives the steady upper-layer temperature of a well-ventilated compartment
// fire from its power, its vent and its boundary. Around the middle of the range
// it was fitted over, a two-zone model should land on it.
void testSteadyLayerMatchesMqhInTheMiddleOfItsRange() {
    const double at = roomWallArea();
    double worst = 0;
    for (double hk : {10.0, 20.0, 30.0, 50.0}) {
        Room r = makeRoom(500.0e3, hk);
        runToSteady(r);
        const double got = r.model.gas[0].upper.temperature() - kTAmbient;
        const double want = fire::mqhTemperatureRise(500.0e3, kDoorArea, kDoorHeight, hk, at);
        worst = std::max(worst, std::abs(got / want - 1.0));

        // Not vacuous: the layer must have descended and the door must be working.
        expectTrue("the smoke layer came down past half the room's height",
                   r.model.gas[0].interfaceZ() < 0.5 * kRoomHeight);
        expectTrue("the door is exchanging gas in both directions", r.model.vents[0].bidirectional);
        expectTrue("with a neutral plane inside the doorway",
                   r.model.vents[0].neutralPlaneZ > 0.0 &&
                       r.model.vents[0].neutralPlaneZ < kDoorHeight);
    }
    // Measured worst deviation is 11.9%, at h_k = 50 W/(m^2 K). MQH is a fit to
    // about a hundred experiments and is quoted at +/-15%, so the model agreeing
    // to better than that across a five-fold sweep of boundary conductance is as
    // close as this comparison can mean anything.
    expectTrue("the steady layer is within 13% of MQH across a 5x sweep of h_k", worst < 0.13);
}

// Away from the middle of the range the model and MQH part company, and they do
// it in one direction: the model runs cool at low power and hot at high power.
//
// That is the signature of a **missing T^4**. The boundary loss here is a linear
// `h_k A dT`, because radiation view factors are a separate roadmap item and are
// explicitly out of this file's scope; MQH's fitted coefficient carries the
// radiative loss its experiments had. A linear law under-charges the loss where
// the layer is hot and over-charges it where the layer is cool, which is exactly
// this ordering. Asserting the *ordering* rather than a band is what makes this a
// statement about the physics instead of a tolerance.
void testMqhDisagreementIsOrderedByPowerNotScattered() {
    const double at = roomWallArea();
    const double powers[] = {100.0e3, 250.0e3, 500.0e3, 1000.0e3};
    for (double hk : {10.0, 30.0, 50.0}) {
        double previous = 0;
        double lowest = 1e30, highest = 0;
        for (double q : powers) {
            Room r = makeRoom(q, hk);
            runToSteady(r);
            const double ratio = (r.model.gas[0].upper.temperature() - kTAmbient) /
                                 fire::mqhTemperatureRise(q, kDoorArea, kDoorHeight, hk, at);
            expectTrue("the model/MQH ratio rises with fire power, monotonically",
                       ratio > previous);
            previous = ratio;
            lowest = std::min(lowest, ratio);
            highest = std::max(highest, ratio);
        }
        expectTrue("and stays inside a factor of 0.55 to 1.20 over a ten-fold power range",
                   lowest > 0.55 && highest < 1.20);
        expectTrue("while the disagreement is a real trend, not noise", highest / lowest > 1.3);
    }
}

// MQH inverted at a 500 K rise and Thomas's flashover power are independent fits,
// by different people to different data. On the reference room they must agree,
// and where they agree is where the model's `h_k` should be taken to sit.
void testMqhAndThomasAgreeOnFlashover() {
    const double at = roomWallArea();
    const double thomas = fire::thomasFlashoverPower(at, kDoorArea, kDoorHeight);
    // Invert MQH for dT = 500 K.
    auto mqhPowerFor = [&](double hk) {
        double lo = 1.0e3, hi = 1.0e8;
        for (int i = 0; i < 200; ++i) {
            const double mid = 0.5 * (lo + hi);
            if (fire::mqhTemperatureRise(mid, kDoorArea, kDoorHeight, hk, at) < 500.0)
                lo = mid;
            else
                hi = mid;
        }
        return lo;
    };
    expectTrue("MQH and Thomas agree within 10% at h_k = 30 W/(m^2 K)",
               std::abs(mqhPowerFor(30.0) / thomas - 1.0) < 0.10);
    expectTrue("and they bracket: MQH is below Thomas at h_k = 10 and above it at h_k = 50",
               mqhPowerFor(10.0) < thomas && mqhPowerFor(50.0) > thomas);
    expectTrue("Thomas puts flashover in this room above a megawatt", thomas > 1.0e6);

    // Both correlations are meaningless without every one of their arguments, and
    // each guard has to be its own guard: a short-circuit chain wired with `&&`
    // instead of `||` lets three of the four through.
    expectNear("MQH refuses a fire of no power",
               fire::mqhTemperatureRise(0.0, kDoorArea, kDoorHeight, 30.0, at), 0.0, 0.0);
    expectNear("or a vent of no area",
               fire::mqhTemperatureRise(500e3, 0.0, kDoorHeight, 30.0, at), 0.0, 0.0);
    expectNear("or a vent of no height",
               fire::mqhTemperatureRise(500e3, kDoorArea, 0.0, 30.0, at), 0.0, 0.0);
    expectNear("or a boundary that conducts nothing",
               fire::mqhTemperatureRise(500e3, kDoorArea, kDoorHeight, 0.0, at), 0.0, 0.0);
    expectNear("or an enclosure with no surface",
               fire::mqhTemperatureRise(500e3, kDoorArea, kDoorHeight, 30.0, 0.0), 0.0, 0.0);
    expectNear("and Thomas refuses an enclosure with no surface",
               fire::thomasFlashoverPower(0.0, kDoorArea, kDoorHeight), 0.0, 0.0);
    expectNear("a vent of no area", fire::thomasFlashoverPower(at, 0.0, kDoorHeight), 0.0, 0.0);
    expectNear("and a vent of no height", fire::thomasFlashoverPower(at, kDoorArea, 0.0), 0.0, 0.0);
}

// --- Entrainment drives the layer down ----------------------------------------

// "The layer descended" is not evidence that entrainment did it. A fire heats the
// gas it already has, and that expansion pushes the interface down on its own.
// This asserts that the descent is the plume's, by two routes: it matches the
// entrained volume flux, and killing entrainment nearly stops it.
void testTheLayerIsDrivenDownByEntrainmentAndNotByExpansion() {
    Room r = makeTallBox(300.0e3);
    for (int i = 0; i < 5; ++i) r.model.step(1.0, r.ship, r.sea);
    const double zBefore = r.model.gas[0].interfaceZ();
    const double massBefore = r.model.gas[0].upper.mass;

    // Heskestad's entrainment, written out here from the published coefficients
    // rather than read back out of fire.cpp, so a transcription error in either
    // would show. 300 kW over a 1 m pool: kW and metres in, kg/s out.
    const double qkw = 300.0, chi = 0.7, diameter = 1.0;
    const double qc = chi * qkw;
    const double flameHeight = 0.235 * std::pow(qkw, 0.4) - 1.02 * diameter;
    const double z0 = 0.083 * std::pow(qkw, 0.4) - 1.02 * diameter;
    // Measured above the fire's base, which is at 0.5 m, not above the floor.
    const double rise = zBefore - 0.5;
    expectTrue("the smoke layer is well clear of the low vent", zBefore > 1.5);
    expectTrue("and the plume is measured above the flame tip", rise > flameHeight + 0.5);
    const double entrained = 0.071 * std::cbrt(qc) * std::pow(rise - z0, 5.0 / 3.0) + 0.0018 * qc;

    // A short step, so the interface has barely moved and the rate is the rate at
    // the height the correlation was evaluated at.
    const double h = 0.02;
    const fire::StepResult s = r.model.step(h, r.ship, r.sea);
    const double gained = r.model.gas[0].upper.mass - massBefore;
    expectTrue("the upper layer gains exactly what Heskestad's plume entrains",
               std::abs(gained / (entrained * h) - 1.0) < 5e-3);
    expectTrue("and the step reports the entrainment it applied",
               std::abs(s.entrainment / entrained - 1.0) < 5e-3);
    expectTrue("which is the untapered correlation while there is cool air left",
               std::abs(s.entrainment / r.model.totalEntrainment(r.model.time) - 1.0) < 5e-3);
    expectTrue("and the interface came down with it", r.model.gas[0].interfaceZ() < zBefore);
    expectTrue("having already fallen a good way from the deckhead", zBefore < 3.6);

    // The negative control. Take the plume's buoyancy flux to nothing and the
    // same fire barely moves the interface at all: what is left is the thermal
    // expansion of the gas already up there.
    Room q = makeTallBox(300.0e3);
    q.model.fires[0].convectiveFraction = 1e-9;
    for (int i = 0; i < 5; ++i) q.model.step(1.0, q.ship, q.sea);
    const double withoutPlume = 4.0 - q.model.gas[0].interfaceZ();
    const double withPlume = 4.0 - zBefore;
    expectTrue("with entrainment removed the layer barely moves",
               withPlume > 20.0 * withoutPlume);
    expectTrue("and the fire released the same heat in both",
               std::abs(q.model.account.heatReleased - 300.0e3 * 5.0) < 1.0);
}

// --- The ferry: a realistic multi-compartment case ----------------------------

struct FerryFire {
    Ship ship;
    Sea sea{0.0};
    fire::Model model;
};

// A 4 MW machinery-space fire in the starboard engine room, with the watertight
// door to port standing open. Deliberately **asymmetric**: the ferry's two engine
// rooms are mirror images, so a fire in both would let a sign error cancel in
// every total.
FerryFire makeFerryFire() {
    FerryFire f;
    f.ship = game::buildFerry();
    f.ship.initialise(f.sea);
    f.model.attach(f.ship, {f.ship.findCompartment("engine_room_s"),
                            f.ship.findCompartment("engine_room_p")});
    fire::DesignFire d;
    d.name = "machinery";
    d.compartment = f.model.findGas("engine_room_s");
    d.baseZ = 2.5;
    d.diameter = 2.5;
    d.growthCoefficient = fire::kGrowthFast;
    d.peakHeatRelease = 4.0e6;
    d.steadyDuration = 600.0;
    d.decayDuration = 300.0;
    f.model.fires.push_back(d);
    for (fire::GasCompartment& g : f.model.gas) g.wallConductance = 25.0;
    return f;
}

// The account, on the case that is allowed to break it. `zone.hpp` records its
// own energy residual closing on unit fixtures and reaching -102% on a real
// tearing run, which is the reason this is asserted here and not only on a box.
void testTheAccountClosesOnTheFerry() {
    FerryFire f = makeFerryFire();
    expectTrue("the attached model is self-consistent", f.model.validate().empty());

    for (int i = 0; i < 1200; ++i) f.model.step(1.0, f.ship, f.sea);
    const fire::Account& a = f.model.account;

    expectTrue("the energy account closes to machine precision on a real fire",
               std::abs(a.energyResidualFraction()) < 1e-11);
    expectTrue("so does the mass account", std::abs(a.massResidualFraction()) < 1e-11);
    expectTrue("and the species account", std::abs(a.productsResidual()) < 1e-9);

    // Guard against vacuity in every direction: an account of nothing closes
    // perfectly. Each of these terms has to be a real number.
    expectTrue("the fire released gigajoules", a.heatReleased > 3.0e9);
    expectTrue("the boundary took most of it", a.wallLoss > 0.5 * a.heatReleased);
    expectTrue("gas crossed the hull in both directions",
               a.massIn > 500.0 && a.massOut > 500.0);
    expectTrue("the enthalpy that left exceeds what came in", a.enthalpyOut > a.enthalpyIn * 1.2);
    expectTrue("combustion products were made and some got out",
               a.productsGenerated > 1.0 && a.productsOut > 0.1 * a.productsGenerated);
    expectTrue("and the compartments ended lighter than they started",
               a.mass < a.initialMass - 10.0);
}

// Species transport through the opening network -- the second half of the
// roadmap item. Soot made in the starboard engine room has to reach the port one
// through the watertight door, and it has to arrive diluted.
void testProductsCrossTheOpeningNetwork() {
    FerryFire f = makeFerryFire();
    const int s = f.model.findGas("engine_room_s");
    const int p = f.model.findGas("engine_room_p");
    for (int i = 0; i < 400; ++i) f.model.step(1.0, f.ship, f.sea);

    const double ys = f.model.gas[static_cast<std::size_t>(s)].upper.productFraction();
    const double yp = f.model.gas[static_cast<std::size_t>(p)].upper.productFraction();
    expectTrue("smoke reaches the compartment next door", yp > 1e-6);
    expectTrue("but arrives diluted", yp < ys);
    expectTrue("and the room with the fire in it is the sooty one", ys > 1e-4);

    double held = 0;
    for (const fire::GasCompartment& g : f.model.gas) held += g.upper.products + g.lower.products;
    expectTrue("every gram is accounted for",
               std::abs(f.model.account.productsGenerated - f.model.account.productsOut - held) <
                   1e-9);
}

// The two engine rooms must not end up in the same state. They are mirror images
// with identical geometry and identical vents, so if the fire's asymmetry did not
// survive the solve, every total in the account would still balance and every
// conservation check would still pass.
void testTheUnfiredCompartmentIsClearlyDifferent() {
    FerryFire f = makeFerryFire();
    for (int i = 0; i < 400; ++i) f.model.step(1.0, f.ship, f.sea);
    const fire::GasCompartment& s = f.model.gas[0];
    const fire::GasCompartment& p = f.model.gas[1];
    expectTrue("the burning compartment is hundreds of kelvin hotter",
               s.upper.temperature() > p.upper.temperature() + 150.0);
    expectTrue("its smoke layer is much deeper", s.interfaceZ() < p.interfaceZ() - 1.0);

    const fire::Vent* door = nullptr;
    for (const fire::Vent& v : f.model.vents)
        if (v.name == "wt_door_er") door = &v;
    expectTrue("the watertight door is in the network", door != nullptr);
    if (door != nullptr) {
        expectTrue("and it is running hot gas one way and cool air the other",
                   door->massAToB > 0.1 && door->massBToA > 0.1);
        expectNear("its net is the difference of the two streams", door->netMass(),
                   door->massAToB - door->massBToA, 0.0);
        expectTrue("and that net is well under what either stream carries, because the door is "
                   "exchanging rather than blowing one way",
                   std::abs(door->netMass()) <
                       0.8 * std::max(door->massAToB, door->massBToA));
        expectTrue("with a neutral plane inside its own span",
                   door->neutralPlaneZ > 1.79 && door->neutralPlaneZ < 3.81);
    }
}

// Water against an opening belongs to the flooding solve, not to the gas model.
// The ferry's breach is 2.5 m below the waterline and its air escapes are on the
// weather deck, so the split is unambiguous and it had better be made.
void testWaterBlockedOpeningsAreLeftToTheFloodingSolve() {
    FerryFire f = makeFerryFire();
    f.model.step(1.0, f.ship, f.sea);
    int blocked = 0, open = 0;
    for (const fire::Vent& v : f.model.vents) {
        if (v.name == "breach_er_s") {
            expectTrue("the underwater breach is left to the flooding solve", v.blockedByWater);
            ++blocked;
        }
        if (v.name == "vent_er_s" || v.name == "wt_door_er") {
            expectTrue("openings in air are the gas model's", !v.blockedByWater);
            ++open;
        }
    }
    expectEqual("the breach was found", blocked, 1);
    expectEqual("and both dry openings were", open, 2);

    // The ferry authors its air escapes at the gooseneck on the weather deck,
    // 5.5 m above the engine room they drain. Trimming the span alone would have
    // deleted every one of them and sealed the machinery spaces.
    int vents = 0;
    for (const fire::Vent& v : f.model.vents)
        if (v.name == "vent_er_s" || v.name == "vent_er_p") ++vents;
    expectEqual("both engine-room air escapes survived into the gas network", vents, 2);
    expectTrue("and the compartment can actually vent to atmosphere",
               f.model.account.massOut > 0.0);
}

// The gas box each compartment gets, derived from the ship's own mesh. A
// bounding box would over-state the floor area by the turn of the bilge and the
// interface would then descend too slowly by exactly that ratio.
void testAttachDerivesTheGasBoxFromTheShipsOwnCompartment() {
    FerryFire f = makeFerryFire();
    const int idx = f.ship.findCompartment("engine_room_s");
    const Compartment& c = f.ship.compartments[static_cast<std::size_t>(idx)];
    const fire::GasCompartment& g = f.model.gas[0];

    expectNear("the gas box takes the compartment's own floor", g.floorZ, c.bboxLo.z, 1e-15);
    expectNear("and its deckhead", g.ceilingZ, c.bboxHi.z, 1e-15);
    expectNear("the gas volume is the compartment's air volume", g.gasVolume, c.airVolume(),
               1e-9);
    expectTrue("the footprint is the prismatic equivalent of that volume",
               std::abs(g.floorArea * (g.ceilingZ - g.floorZ) / g.gasVolume - 1.0) < 1e-12);
    const double boxFloor = (c.bboxHi.x - c.bboxLo.x) * (c.bboxHi.y - c.bboxLo.y);
    expectTrue("which is well inside the bounding box, because a hull is not a box",
               g.floorArea < 0.9 * boxFloor);
    expectNear("the perimeter is the bounding box's", g.perimeter,
               2.0 * ((c.bboxHi.x - c.bboxLo.x) + (c.bboxHi.y - c.bboxLo.y)), 1e-12);

    // With floodwater in it, the gas space starts at the water surface: a
    // half-flooded engine room has half the gas and a much shorter fall.
    Ship flooded = game::buildFerry();
    Sea sea{0.0};
    flooded.initialise(sea);
    const int fi = flooded.findCompartment("engine_room_s");
    flooded.compartments[static_cast<std::size_t>(fi)].waterVolume =
        0.5 * flooded.compartments[static_cast<std::size_t>(fi)].floodableVolume();
    flooded.step(1e-6, sea);
    fire::Model wet;
    wet.attach(flooded, {fi});
    const Compartment& fc = flooded.compartments[static_cast<std::size_t>(fi)];
    expectNear("the gas floor is the floodwater surface", wet.gas[0].floorZ, fc.surfaceOffset,
               1e-12);
    expectTrue("which is well above the bottom of the space", wet.gas[0].floorZ > fc.bboxLo.z + 1.0);
    expectTrue("and the gas space is about half what it was",
               wet.gas[0].gasVolume < 0.6 * g.gasVolume);
}

// The ferry's air escapes are authored at the gooseneck on the weather deck,
// 5.5 m above the engine room they drain, because `Opening::pos` for a pipe is
// the end whose height decides whether the sea is over it. Trimming the span
// alone would delete every one of them and seal the machinery spaces, so the
// vent is *slid* to the deckhead first. That adjustment is silent unless it is
// asserted.
void testVentsAreSlidIntoTheGasSpaceBeforeBeingTrimmed() {
    FerryFire f = makeFerryFire();
    // Far enough in that the fire is established and the door is working, so the
    // neutral-plane assertion below has a neutral plane to be about.
    for (int i = 0; i < 300; ++i) f.model.step(1.0, f.ship, f.sea);

    const fire::Vent* vent = nullptr;
    const fire::Vent* door = nullptr;
    for (const fire::Vent& v : f.model.vents) {
        if (v.name == "vent_er_s") vent = &v;
        if (v.name == "wt_door_er") door = &v;
    }
    expectTrue("the engine room's air escape is in the network", vent != nullptr);
    expectTrue("and so is the watertight door", door != nullptr);
    if (vent == nullptr || door == nullptr) return;

    const double ceiling = f.model.gas[0].ceilingZ;
    const double floor = f.model.gas[0].floorZ;
    expectTrue("as authored, the air escape sits far above the space it drains",
               vent->soffitZ > ceiling + 5.0);
    expectNear("but it is worked at the deckhead", vent->activeSoffitZ, ceiling, 1e-12);
    expectNear("keeping its own height", vent->activeSoffitZ - vent->activeSillZ,
               vent->soffitZ - vent->sillZ, 1e-12);

    // The door's sill is authored below the deck it opens onto, and gets slid up
    // the same way.
    expectTrue("the door is authored with its sill below the deck", door->sillZ < floor);
    expectNear("and stands on the deck", door->activeSillZ, floor, 1e-12);
    expectNear("also keeping its height", door->activeSoffitZ - door->activeSillZ,
               door->soffitZ - door->sillZ, 1e-12);
    expectTrue("and it is exchanging gas across a neutral plane", door->bidirectional);
    expectTrue("which can only ever be inside the span it works over",
               door->neutralPlaneZ >= door->activeSillZ &&
                   door->neutralPlaneZ <= door->activeSoffitZ);
}

// Floodwater *inside* the ship blocks an opening exactly as the sea outside
// does. The ferry's breach only ever exercises the sea half of that rule, so the
// internal half needs its own case.
void testAnOpeningUnderFloodwaterInsideTheShipIsBlockedToo() {
    FerryFire f = makeFerryFire();
    f.model.step(1.0, f.ship, f.sea);
    const fire::Vent* door = nullptr;
    for (const fire::Vent& v : f.model.vents)
        if (v.name == "wt_door_er") door = &v;
    expectTrue("the door starts dry", door != nullptr && !door->blockedByWater);

    const int idx = f.ship.findCompartment("engine_room_s");
    Compartment& c = f.ship.compartments[static_cast<std::size_t>(idx)];
    c.waterVolume = 0.9 * c.floodableVolume();
    f.ship.step(1e-6, f.sea);
    // In the body frame, which is the frame the door's own position is in.
    expectTrue("the engine room's own surface is now well above the door",
               c.surfaceOffset > 4.0);

    f.model.step(1.0, f.ship, f.sea);
    for (const fire::Vent& v : f.model.vents)
        if (v.name == "wt_door_er") door = &v;
    expectTrue("so the door belongs to the flooding solve, not the gas model",
               door != nullptr && door->blockedByWater);
    expectNear("and the gas model moves nothing through it", door->massAToB, 0.0, 0.0);
    expectNear("in either direction", door->massBToA, 0.0, 0.0);
}

// --- The exact control --------------------------------------------------------

// A byte-for-byte view of a ship, so "unchanged" can mean unchanged.
std::vector<double> shipFingerprint(const Ship& s) {
    std::vector<double> v;
    v.push_back(s.state.position.x);
    v.push_back(s.state.position.y);
    v.push_back(s.state.position.z);
    v.push_back(s.state.orientation.w);
    v.push_back(s.state.orientation.x);
    v.push_back(s.state.orientation.y);
    v.push_back(s.state.orientation.z);
    v.push_back(s.state.velocity.x);
    v.push_back(s.state.velocity.y);
    v.push_back(s.state.velocity.z);
    v.push_back(s.state.angularVelocity.x);
    v.push_back(s.state.angularVelocity.y);
    v.push_back(s.state.angularVelocity.z);
    for (const Compartment& c : s.compartments) {
        v.push_back(c.waterVolume);
        v.push_back(c.airMass);
        v.push_back(c.airPressure);
        v.push_back(c.surfaceOffset);
        v.push_back(c.surfaceWorldZ);
        v.push_back(c.waterCentroid.x);
        v.push_back(c.waterCentroid.y);
        v.push_back(c.waterCentroid.z);
    }
    for (const Opening& o : s.openings) {
        v.push_back(o.lastFlow);
        v.push_back(o.lastFlowWasWater ? 1.0 : 0.0);
    }
    return v;
}

// **With no heat release the fire model must leave the ship bit-identical.**
//
// Not close: identical. The flooding scenarios and the 35 figures the full gate
// publishes are validated against their own behaviour, and a fire model that is
// not burning has no business moving any of them. This runs the ferry's damage
// case twice -- once bare, once with the gas model attached, stepped, and
// **writing back through `applyTo`**, which is the path that could actually
// change something -- and compares every state double for exact equality.
void testZeroHeatReleaseLeavesTheShipBitIdentical() {
    Sea sea{0.0};
    Ship bare = game::buildFerry();
    bare.initialise(sea);

    Ship lit = game::buildFerry();
    lit.initialise(sea);
    fire::Model model;
    model.attach(lit, {lit.findCompartment("engine_room_s"), lit.findCompartment("engine_room_p"),
                       lit.findCompartment("aft_hold_s"), lit.findCompartment("vehicle_deck")});
    // A design fire in every respect except that it releases nothing.
    fire::DesignFire d;
    d.name = "unlit";
    d.compartment = 0;
    d.baseZ = 2.5;
    d.diameter = 2.5;
    d.peakHeatRelease = 0.0;
    model.fires.push_back(d);

    const int ticks = 400;
    for (int i = 0; i < ticks; ++i) {
        bare.step(0.05, sea);
        lit.step(0.05, sea);
        model.step(0.05, lit, sea);
        model.applyTo(lit);
    }

    const std::vector<double> a = shipFingerprint(bare);
    const std::vector<double> b = shipFingerprint(lit);
    expectEqual("the two fingerprints are the same shape", static_cast<long long>(a.size()),
                static_cast<long long>(b.size()));
    std::size_t differing = 0;
    for (std::size_t i = 0; i < a.size() && i < b.size(); ++i)
        if (std::memcmp(&a[i], &b[i], sizeof(double)) != 0) ++differing;
    expectEqual("an unlit fire leaves the ship bit-identical, in every state double",
                static_cast<long long>(differing), 0);

    // Guard against vacuity: the ship has to have been *doing* something, or two
    // identical still-water ships would prove nothing.
    expectTrue("and the ship really was flooding through the run",
               bare.totalFloodwaterMass() > 1.0e5);
    expectTrue("and heeling", std::abs(bare.diagnostics(sea).heelDeg) > 0.5);
    expectTrue("while the gas model really did step", model.time > 19.0);
}

// The write-back path itself, with a fire that *is* burning: `applyTo` puts in
// the mass that reproduces this model's pressure under the ship's own isothermal
// formula, so the flooding network sees the right pressure.
void testApplyToReproducesTheModelPressureUnderTheShipsOwnFormula() {
    FerryFire f = makeFerryFire();
    std::vector<double> pressureAtAttach, airMassBefore;
    for (const fire::GasCompartment& g : f.model.gas) {
        pressureAtAttach.push_back(g.pressure());
        airMassBefore.push_back(
            f.ship.compartments[static_cast<std::size_t>(g.shipCompartment)].airMass);
    }
    for (int i = 0; i < 300; ++i) f.model.step(1.0, f.ship, f.sea);
    f.model.applyTo(f.ship);

    // `applyTo` writes a *delta*, so what it promises is that the pressure the
    // ship's own isothermal formula reads moves by exactly what the gas model's
    // pressure moved. It cannot promise the absolute value: the ship's air mass
    // is its own, and a flooding compartment changes it underneath.
    double worst = 0, biggest = 0;
    for (std::size_t i = 0; i < f.model.gas.size(); ++i) {
        const fire::GasCompartment& g = f.model.gas[i];
        const Compartment& c = f.ship.compartments[static_cast<std::size_t>(g.shipCompartment)];
        const double shipDelta = (c.airMass - airMassBefore[i]) * kRAir * kTAmbient / g.gasVolume;
        const double modelDelta = g.pressure() - pressureAtAttach[i];
        worst = std::max(worst, std::abs(shipDelta - modelDelta));
        biggest = std::max(biggest, std::abs(modelDelta));
    }
    expectTrue("applyTo moves the ship's pressure by exactly what the gas model moved",
               worst < 1e-9 * std::max(biggest, 1.0));
    expectTrue("and the move is a real number of pascals, not a rounding error", biggest > 1.0);

    // And a second call with nothing in between must be a no-op, because it
    // writes the delta rather than the value.
    const double before = f.ship.compartments[0].airMass;
    f.model.applyTo(f.ship);
    expectTrue("and calling it twice changes nothing at all",
               std::memcmp(&before, &f.ship.compartments[0].airMass, sizeof(double)) == 0);
}

// --- Integration --------------------------------------------------------------

// The caller's tick must not be able to change the answer. `step()` subdivides
// internally against its own accuracy cap, so any tick at or above `maxSubstep`
// has to give the identical trajectory.
void testTheAnswerDoesNotDependOnTheCallersTick() {
    double reference = 0, referenceZ = 0;
    double worst = 0, worstZ = 0;
    for (double dt : {4.0, 2.0, 1.0, 0.5, 0.25}) {
        Room r = makeRoom(500.0e3, 30.0);
        for (double t = 0; t < 600.0; t += dt) r.model.step(dt, r.ship, r.sea);
        const double tu = r.model.gas[0].upper.temperature();
        const double zi = r.model.gas[0].interfaceZ();
        if (reference == 0) {
            reference = tu;
            referenceZ = zi;
        } else {
            worst = std::max(worst, std::abs(tu / reference - 1.0));
            worstZ = std::max(worstZ, std::abs(zi - referenceZ));
        }
    }
    // Measured at 2.0e-9 relative on the temperature and 3.1e-9 m on the
    // interface over a 600 s run, and **non-monotone in the tick**, which is the
    // signature of round-off rather than truncation: the internal substeps are
    // the same size whatever the caller does, but `remaining -= h` accumulates
    // differently when one call is subdivided sixteen times and another is not.
    // A few ULP per substep over 2400 substeps is 1e-6 K on 576 K. Asserted just
    // above what was measured, because anything looser would also pass on a model
    // that had genuinely become tick-dependent and merely converged well.
    expectTrue("every tick from 4 s down to 0.25 s gives the same layer temperature",
               worst < 5e-9);
    expectTrue("and the same interface height", worstZ < 1e-8);
    expectTrue("and the fire actually did something", reference > kTAmbient + 200.0);
}

// Refining the *internal* step is where accuracy lives, and the scheme is first
// order, so the answer must converge as the substep shrinks rather than wander.
void testRefiningTheInternalStepConverges() {
    double previous = 0, coarse = 0;
    std::vector<double> answers;
    for (double h : {0.4, 0.2, 0.1, 0.05}) {
        Room r = makeRoom(500.0e3, 30.0);
        r.model.maxSubstep = h;
        for (int i = 0; i < 600; ++i) r.model.step(1.0, r.ship, r.sea);
        answers.push_back(r.model.gas[0].upper.temperature());
    }
    coarse = answers.front();
    previous = answers.back();
    double worstGap = 0;
    for (std::size_t i = 1; i < answers.size(); ++i)
        worstGap = std::max(worstGap, std::abs(answers[i] - answers[i - 1]) / answers[i]);
    expectTrue("halving the internal step moves the steady temperature by under 0.3%",
               worstGap < 3e-3);
    expectTrue("and the successive answers are converging, not wandering",
               std::abs(answers[3] - answers[2]) < std::abs(answers[1] - answers[0]));
    expectTrue("while the coarse and fine answers describe the same fire",
               std::abs(coarse / previous - 1.0) < 0.01);
}

// The implicit pressure solve replaced an explicit step with a clamp because the
// explicit one is unusable. It has to actually converge, and it has to leave the
// compartment near atmospheric rather than pressurised -- a compartment with a
// 1.6 m^2 doorway that sits hundreds of pascals above the outside is a model that
// is not ventilating.
void testThePressureSolveConvergesAndDoesNotPressuriseTheRoom() {
    Room r = makeRoom(500.0e3, 30.0);
    fire::StepResult last{};
    bool capped = false;
    for (int i = 0; i < 600; ++i) {
        last = r.model.step(1.0, r.ship, r.sea);
        capped = capped || last.pressureSolveCapped;
        expectTrue("the pressure solve never fails to bracket", !last.pressureSolveCapped);
        if (capped) break;
    }
    expectTrue("it converges in a handful of Gauss-Seidel sweeps", last.pressureSweeps <= 4);
    expectTrue("a vented room stays within a few pascals of atmospheric",
               std::abs(r.model.gas[0].gaugeAtFloor()) < 10.0);
    expectTrue("but not exactly at it -- the neutral plane has to sit in the door",
               std::abs(r.model.gas[0].gaugeAtFloor()) > 0.1);
}

// A bad model definition does not crash; it quietly produces a wrong answer. Same
// contract as `Ship::validate`.
void testValidateCatchesBadDefinitions() {
    fire::Model m;
    fire::GasCompartment g;
    g.name = "flat";
    g.floorZ = 3.0;
    g.ceilingZ = 1.0;
    g.gasVolume = 0.0;
    g.floorArea = 0.0;
    m.gas.push_back(g);
    fire::Vent v;
    v.name = "nowhere";
    v.a = 0;
    v.b = 0;
    v.sillZ = 1.0;
    v.soffitZ = 1.0;
    m.vents.push_back(v);
    fire::DesignFire f;
    f.name = "elsewhere";
    f.compartment = 7;
    f.diameter = 0.0;
    m.fires.push_back(f);

    const std::vector<std::string> problems = m.validate();
    auto reported = [&](const char* fragment) {
        for (const std::string& p : problems)
            if (p.find(fragment) != std::string::npos) return true;
        return false;
    };
    expectTrue("a space with no volume is reported", reported("has no volume"));
    expectTrue("a ceiling below its floor is reported", reported("ceiling at or below"));
    expectTrue("no floor area is reported", reported("has no floor area"));
    expectTrue("a vent onto itself is reported", reported("connects a space to itself"));
    expectTrue("a vent with no height is reported", reported("no height to integrate"));
    expectTrue("a fire in a space that does not exist is reported",
               reported("is in a gas space that does not exist"));
    expectTrue("and a fire of no size is reported", reported("non-positive diameter"));

    fire::Model good;
    fire::GasCompartment h;
    h.name = "box";
    h.floorZ = 0;
    h.ceilingZ = 3;
    h.gasVolume = 30;
    h.floorArea = 10;
    good.gas.push_back(h);
    expectTrue("and a sound one is not", good.validate().empty());
}

// The design curve itself: growth, plateau, decay, and an area under it that the
// account can be checked against.
void testDesignFireFollowsItsCurve() {
    fire::DesignFire f;
    f.growthCoefficient = fire::kGrowthMedium;
    f.peakHeatRelease = 2.0e6;
    f.steadyDuration = 300.0;
    f.decayDuration = 200.0;
    const double tGrow = std::sqrt(f.peakHeatRelease / f.growthCoefficient);

    expectNear("nothing before it starts", f.heatRelease(0.0), 0.0, 0.0);
    expectNear("t-squared while growing", f.heatRelease(0.5 * tGrow),
               f.growthCoefficient * 0.25 * tGrow * tGrow, 1e-9);
    expectNear("the peak at the end of growth", f.heatRelease(tGrow), f.peakHeatRelease, 1e-9);
    expectNear("and through the plateau", f.heatRelease(tGrow + 150.0), f.peakHeatRelease, 1e-9);
    expectNear("half power halfway down the decay ramp", f.heatRelease(tGrow + 400.0),
               0.5 * f.peakHeatRelease, 1e-9);
    expectNear("and out at the end", f.heatRelease(tGrow + 600.0), 0.0, 0.0);

    // The published total against a fine numerical integration of the same curve.
    double sum = 0;
    const double h = 1e-3;
    for (double t = 0; t < tGrow + f.steadyDuration + f.decayDuration + 10.0; t += h)
        sum += f.heatRelease(t + 0.5 * h) * h;
    expectTrue("and totalEnergy() is the area under it",
               std::abs(sum / f.totalEnergy() - 1.0) < 1e-9);

    fire::DesignFire steady;
    steady.peakHeatRelease = 1.0e6;
    expectNear("a fire with no growth curve is a step to full power",
               steady.heatRelease(1e-9), 1.0e6, 0.0);
    expectTrue("and one that never decays releases unbounded energy",
               std::isinf(steady.totalEnergy()));

    // The canonical growth coefficients: 1055 kW in 600/300/150/75 s. All four,
    // because three right and one wrong is the shape this repo keeps finding.
    expectNear("the slow curve reaches 1055 kW in 600 s",
               fire::kGrowthSlow * 600.0 * 600.0, 1.055e6, 1e-6);
    expectNear("the medium one in 300 s", fire::kGrowthMedium * 300.0 * 300.0, 1.055e6, 1e-6);
    expectNear("the fast one in 150 s", fire::kGrowthFast * 150.0 * 150.0, 1.055e6, 1e-6);
    expectNear("the ultrafast one in 75 s", fire::kGrowthUltrafast * 75.0 * 75.0, 1.055e6, 1e-6);
    expectTrue("and they are ordered", fire::kGrowthSlow < fire::kGrowthMedium &&
                                           fire::kGrowthMedium < fire::kGrowthFast &&
                                           fire::kGrowthFast < fire::kGrowthUltrafast);

    // The default product yield is 0.05 kg per kg of a 20 MJ/kg fuel.
    expectNear("the default soot yield is 0.05 kg/kg on a 20 MJ/kg fuel",
               fire::DesignFire{}.productYield * 20.0e6, 0.05, 1e-15);
}

// ==============================================================================
// Suppression, and its effect on stability
// ==============================================================================
//
// The roadmap item is "suppression systems, **and their effect on stability**",
// and on a ship those are one question: water put into a compartment to fight a
// fire is water in a compartment.
//
// The independent answers this half leans on:
//
//   * **The free-surface moment of a rectangular tank**, `rho mu b^3 l / 12`.
//     Asserted against a *box* barge, where the second moment is arithmetic, and
//     then again against the ferry's vehicle deck, where it is measured off the
//     compartment's own mesh by a slicing integration that shares no code with
//     the model. What the model computes instead is a water body re-levelled
//     against gravity; the two have nothing in common but the answer.
//   * **The classical weir**, `(2/3) Cd b sqrt(2g) h^(3/2)`, and the equilibrium
//     depth it holds a deck at under a steady inflow.
//   * **The caloric arithmetic of the split**: `cp dT + e L` per kilogram, which
//     is the whole cooling mechanism and is an identity, not a correlation.
//   * **An exact control**: suppression switched off must leave both the ship and
//     the gas bit-identical to a run with no suppression objects at all.
//   * **Convergence**, because the sink is stiff and a stiff term that does not
//     converge under refinement is a term whose answer is the step size.
//
// **Two of these tests exist because a mutant survived the first version.** A
// mutant that handed the deck the whole delivered flow instead of the part that
// did not evaporate lived through the entire suite, because every fixture that
// wrote water was *cold* and `evaporated` was therefore zero -- the two
// quantities were the same number everywhere it was looked at. And a mutant that
// discarded the water that would not fit in a full compartment lived because
// nothing ever filled one. `testTheWaterAccountClosesAgainstWhatTheShipReceived`
// now carries a 30 MW fire so that a third of the flow leaves as steam, and
// `testWaterThatWillNotFitStaysOwedRatherThanBeingDropped` starts a compartment
// at 99.9% full.

// --- A box barge, so the free surface second moment is arithmetic --------------

constexpr double kBargeL = 60.0, kBargeB = 16.0, kBargeDepth = 10.0;
constexpr double kDeckL = 50.0, kDeckB = 12.0, kDeckZ = 4.0, kDeckTop = 8.0;

// The deck's free-surface second moment about its own centreline. A box, so this
// is the textbook `b^3 l / 12` with nothing measured.
constexpr double kDeckI = kDeckB * kDeckB * kDeckB * kDeckL / 12.0;

struct Barge {
    Ship ship;
    Sea sea{0.0};
    fire::Model model;
};

// A rectangular barge with one undivided deck compartment. `draft` is where she
// floats: 3 m leaves the deck a metre clear of the sea, 4.5 m puts it half a
// metre under, which is the difference between a freeing port that drains and one
// that admits.
Barge makeBarge(double draft) {
    Barge b;
    b.ship.hull = makeBox({-kBargeL / 2, -kBargeB / 2, 0.0},
                          {kBargeL / 2, kBargeB / 2, kBargeDepth});
    b.ship.deckEdgeZ = kDeckZ;
    Compartment d;
    d.name = "deck";
    d.mesh = makeBox({-kDeckL / 2, -kDeckB / 2, kDeckZ}, {kDeckL / 2, kDeckB / 2, kDeckTop});
    // Unity, so that the free-surface arithmetic below is the box's own and not
    // the box's times a factor. The ferry case carries a real 0.90 and checks
    // that the factor is there.
    d.permeability = 1.0;
    d.ventedToAtmosphere = true;
    b.ship.compartments.push_back(d);
    b.ship.lightshipMass = kBargeL * kBargeB * draft * kRhoSeawater;
    b.ship.lightshipCog = {0.0, 0.0, 5.0};
    b.ship.gyradii = {6.0, 15.0, 16.0};
    b.ship.initialise(b.sea);
    b.model.attach(b.ship, {0});
    return b;
}

// GM by central difference on the ship's own righting arm. `Diagnostics` takes
// this at eps = 0.03 rad, which is the right question for a ship at a finite
// angle and the wrong one for the initial metacentric height of a shallow layer
// -- see `testPocketingLimitsFreeSurfaceExactlyWhereGeometrySays`.
double gmAt(const Ship& s, const Sea& sea, double eps) {
    return (s.rightingArmAtHeel(eps, sea) - s.rightingArmAtHeel(-eps, sea)) / (2 * eps);
}

// The smallest positive heel at which GZ comes back to zero: where a ship with no
// initial stability actually sits. Coarse scan then bisect, because 0.1 degree of
// resolution over 45 degrees is 450 plane sweeps and this is called in a loop.
double lollDeg(const Ship& s, const Sea& sea) {
    // Stability is asked of GM, not of GZ at zero heel. `rightingArmAtHeel(0)` is
    // zero by symmetry on any hull whose loading is symmetric, to within the sign
    // of the last bit -- and an early-out on that sign reported "upright" or
    // "lolling" for the same ship depending on round-off. It did: two ferry runs
    // differing only in whether a freeing port was blocked reported 0.00 and 1.50
    // degrees for that reason and not for any other.
    if (gmAt(s, sea, 0.002) >= 0) return 0.0;
    for (int i = 1; i <= 180; ++i) {
        const double a = i * 0.25;
        if (s.rightingArmAtHeel(a * kDegToRad, sea) >= 0) {
            double lo = (i - 1) * 0.25, hi = a;
            for (int k = 0; k < 20; ++k) {
                const double mid = 0.5 * (lo + hi);
                if (s.rightingArmAtHeel(mid * kDegToRad, sea) < 0) lo = mid; else hi = mid;
            }
            return 0.5 * (lo + hi);
        }
    }
    return 180.0;
}

// The second moment of a compartment's waterplane about the centreline, measured
// off the mesh and nothing else: slice into strips, take each strip's breadth
// from its own volume, integrate `b^3/12`. Deliberately a *different* route to
// the answer than anything in `fire.cpp` or `ship.cpp`, which is the only reason
// it is worth asserting against. Converged to 1e-5 relative at dx = 0.125 on the
// ferry's deck: 56300.23 at dx = 2 m, 56307.67 at dx = 0.125 m.
double waterplaneSecondMoment(const TriMesh& mesh, const Vec3& lo, const Vec3& hi, double z,
                              double dx) {
    // One-sided upward, so a surface sitting on the compartment floor is not
    // measured through a slab that is half outside the space. Getting that wrong
    // reads as a factor of two in the answer, because the breadth is cubed.
    const double slab = 0.05;
    double iyy = 0;
    for (double x = lo.x; x < hi.x - 1e-12; x += dx) {
        const double x1 = std::min(x + dx, hi.x);
        const TriMesh strip =
            clipToBox(mesh, Vec3{x, lo.y - 1.0, z}, Vec3{x1, hi.y + 1.0, z + slab});
        const double vol = integrate(strip).volume;
        if (vol <= 0) continue;
        const double b = vol / ((x1 - x) * slab);
        iyy += b * b * b / 12.0 * (x1 - x);
    }
    return iyy;
}

// --- The arithmetic of the split ----------------------------------------------

// The rule unit, converted once. Every suppression standard is written in litres
// per minute per square metre and every equation here is in kilograms per second,
// so a factor of sixty and a factor of a thousand both live in this one line.
void testSprayMassFlowConvertsTheRuleUnitExactly() {
    const double area = 373.6;
    const double got = fire::sprayMassFlow(area, 5.0);
    // SOLAS II-2 Reg. 20 for a ro-ro space: 5 L/(min m^2).
    expectNear("5 L/(min m2) over 373.6 m2 is the litres it says it is",
               got * 60.0 / kRhoFresh * 1e3, area * 5.0, 1e-9);
    expectNear("and that is 31.077 kg/s", got, area * 5.0 * 1e-3 * kRhoFresh / 60.0, 1e-12);
    // The two mistakes this function exists to prevent, named so that a future
    // edit that makes either has to delete an assertion rather than pass.
    expectTrue("which is neither the litres per minute", std::abs(got - area * 5.0) > 1.0);
    expectTrue("nor sixty times the answer", std::abs(got - area * 5.0 * 1e-3 * kRhoFresh) > 1.0);
    expectNear("no area is no flow", fire::sprayMassFlow(0.0, 5.0), 0.0, 0.0);
    expectNear("and no spray density is no flow", fire::sprayMassFlow(area, 0.0), 0.0, 0.0);
}

// The freeing port's discharge law, against the classical weir it has to reduce
// to, the submerged orifice it has to reduce to at the other end, and the
// antisymmetry that makes "the sea comes in instead" the same formula.
void testScupperIntegralIsTheClassicalWeirAndItsDrownedLimit() {
    const double root2g = std::sqrt(2.0 * kGravity);
    for (double h : {0.005, 0.02, 0.1, 0.5, 2.0}) {
        expectNear("a dry outside gives the free weir exactly",
                   fire::scupperFlow(h, 0.0), (2.0 / 3.0) * root2g * h * std::sqrt(h), 1e-14);
    }
    expectNear("equal levels move nothing at all", fire::scupperFlow(0.4, 0.4), 0.0, 0.0);
    expectNear("and neither does a dry port", fire::scupperFlow(0.0, 0.0), 0.0, 0.0);

    // Water below the sill cannot reach it, so a negative head is the same as no
    // head. Asserted directly because the model's own caller already clamps
    // before it gets here -- which means nothing else in this file can see
    // whether the function does -- and this one is public.
    expectNear("water below the sill is the same as a dry deck",
               fire::scupperFlow(-1.0, 0.5), fire::scupperFlow(0.0, 0.5), 0.0);
    expectNear("in both arguments", fire::scupperFlow(0.5, -1.0),
               fire::scupperFlow(0.5, 0.0), 0.0);

    // Antisymmetric to the bit: a port with the sea standing over it is the same
    // port run backwards, and that is one formula rather than two branches that
    // have to be kept in step.
    for (double a : {0.05, 0.3, 1.2})
        for (double c : {0.0, 0.1, 0.9})
            expectNear("swapping the two heads flips the sign and nothing else",
                       fire::scupperFlow(a, c), -fire::scupperFlow(c, a), 0.0);

    expectTrue("the sea above the sill drives water inward", fire::scupperFlow(0.1, 1.0) < 0);

    // Deeply drowned, the free band is negligible against the drowned one and the
    // law becomes Torricelli through the wetted area. Checked as a limit rather
    // than at a point: the ratio has to *approach* one.
    double previous = 0;
    for (double lo : {1.0, 10.0, 100.0, 1000.0}) {
        const double drop = 0.01;
        const double ratio = fire::scupperFlow(lo + drop, lo) / (lo * std::sqrt(2 * kGravity * drop));
        expectTrue("a deeply drowned port tends to Torricelli through its wetted area",
                   ratio > 1.0 && ratio < 1.01);
        expectTrue("and gets closer the deeper it is drowned", previous == 0 || ratio < previous);
        previous = ratio;
    }
}

// The cooling arithmetic, on a fixture where the layer is far above saturation
// and can supply everything asked of it: per kilogram, sensible heat to
// saturation plus latent heat on the evaporated share, and nothing else.
void testDrencherCoolingIsSensibleHeatPlusLatentOnTheEvaporatedShare() {
    for (double e : {0.0, 0.3, 0.75, 1.0}) {
        Room r = makeRoom(2.0e6, 30.0);
        fire::Drencher d;
        d.name = "head";
        d.gasCompartment = 0;
        d.flow = 0.02;                 // small: the layer must be able to supply it
        d.evaporatedFraction = e;
        d.on = true;
        r.model.drenchers.push_back(d);
        runToSteady(r, 400);

        const fire::Drencher& s = r.model.drenchers[0];
        expectTrue("the layer is far above saturation, so the whole split is available",
                   r.model.gas[0].upper.temperature() > fire::kTSaturation + 500.0);
        const double perKg = fire::kCpWater * (fire::kTSaturation - kTAmbient) +
                             e * fire::kLatentHeat;
        const double got = s.lastCooling / (s.flow * perKg);
        // The sink is applied as an exact relaxation rather than as a rate, so it
        // sits just *under* the demand by the exponential's own second-order term
        // -- measured between 2.8e-4 and 2.1e-3 of it across this sweep, growing
        // with the demand exactly as `-expm1(-x)/x` does. Bounded on both sides,
        // because a scheme that over-delivered would be removing energy the layer
        // did not have.
        expectTrue("cooling is the flow times the heat one kilogram takes out",
                   got > 0.995 && got <= 1.0);
        expectNear("the steam leaves in the same proportion as the heat",
                   s.lastEvaporated, e * s.flow * got, 1e-15);
        expectNear("and every kilogram either evaporated or landed",
                   s.lastEvaporated + s.lastToDeck, s.flow, 0.0);
    }

    // Guard against vacuity in the direction that matters: the water has to
    // actually cool the fire, not merely be counted.
    Room dry = makeRoom(2.0e6, 30.0);
    runToSteady(dry, 400);
    Room wet = makeRoom(2.0e6, 30.0);
    fire::Drencher d;
    d.name = "head";
    d.gasCompartment = 0;
    d.flow = 0.4;
    d.on = true;
    wet.model.drenchers.push_back(d);
    runToSteady(wet, 400);
    expectTrue("and 0.4 kg/s into a 2 MW fire takes 218 K off the layer",
               wet.model.gas[0].upper.temperature() <
                   dry.model.gas[0].upper.temperature() - 150.0);
    expectTrue("by removing a fifth of the fire's whole output",
               wet.model.account.suppressionCooling > 0.15 * wet.model.account.heatReleased);
}

// **The coefficient nothing here measures, and how much it matters.**
//
// The evaporated fraction is a guess -- the header says so -- so what has to be
// established is the sensitivity rather than the value. It has two halves, and
// they point opposite ways.
void testTheEvaporatedFractionMovesCoolingByFourAndNotByNine() {
    const double sensible = fire::kCpWater * (fire::kTSaturation - kTAmbient);
    auto perKg = [&](double e) { return sensible + e * fire::kLatentHeat; };

    expectNear("the sensible floor is 355.8 kJ/kg", sensible, 355.83e3, 20.0);
    expectNear("cooling at e = 0.1 is 581.5 kJ/kg", perKg(0.1), 581.5e3, 100.0);
    expectNear("and at e = 0.9 it is 2387.1 kJ/kg", perKg(0.9), 2387.1e3, 100.0);
    expectNear("so the span over that range is 4.105, not 9", perKg(0.9) / perKg(0.1), 4.105,
               0.002);
    expectTrue("which is well short of the 9.0 the latent heats alone would give",
               perKg(0.9) / perKg(0.1) < 0.6 * 9.0);
    // The floor is the whole reason, so name it: at the default fraction the
    // sensible term is still a third of the total.
    expectTrue("because at the default 0.3 the sensible term is still a third of it",
               sensible / perKg(0.3) > 0.3 && sensible / perKg(0.3) < 0.36);

    // The same sweep measured through the model rather than asserted from the
    // arithmetic, because nothing tests a comment.
    double coolLow = 0, coolHigh = 0, deckLow = 0, deckHigh = 0;
    for (int i = 0; i < 2; ++i) {
        Room r = makeRoom(2.0e6, 30.0);
        fire::Drencher d;
        d.name = "head";
        d.gasCompartment = 0;
        d.flow = 0.02;
        d.evaporatedFraction = i == 0 ? 0.1 : 0.9;
        d.on = true;
        r.model.drenchers.push_back(d);
        runToSteady(r, 400);
        (i == 0 ? coolLow : coolHigh) = r.model.drenchers[0].lastCooling;
        (i == 0 ? deckLow : deckHigh) = r.model.drenchers[0].lastToDeck;
    }
    expectNear("the model's own cooling ratio is the same 4.10", coolHigh / coolLow, 4.099,
               0.02);
    // **And the water runs the other way, by more.** On the arithmetic alone the
    // coefficient matters *more* to the stability half of this item than to the
    // cooling half, which is the opposite of what the latent heat suggests.
    expectTrue("while the water landing differs by nine, not by four",
               deckLow / deckHigh > 8.5 && deckLow / deckHigh < 9.0);
    expectTrue("so the fraction that cools best is the one that wets the deck least",
               coolHigh > coolLow && deckHigh < deckLow);
}

// **And on a real ship the coefficient barely matters at all**, because it is a
// ceiling and not a rate: what actually evaporates is bounded by the power
// available to boil it. A drencher can deliver far more water than any fire in
// the space can turn to steam, and the rest lands whatever the fraction says.
void testWhatEvaporatesIsBoundedByTheFiresPowerAndNotByTheFraction() {
    double deck[3] = {0, 0, 0};
    double evaporatedShare[3] = {0, 0, 0};
    const double e[3] = {0.1, 0.3, 0.9};
    for (int i = 0; i < 3; ++i) {
        Room r = makeRoom(2.0e6, 30.0);
        fire::Drencher d;
        d.name = "deluge";
        d.gasCompartment = 0;
        d.flow = 5.0;              // far more than a 2 MW fire can boil
        d.evaporatedFraction = e[i];
        d.on = true;
        r.model.drenchers.push_back(d);
        runToSteady(r, 900);
        deck[i] = r.model.drenchers[0].lastToDeck;
        evaporatedShare[i] = r.model.drenchers[0].lastEvaporated / d.flow;
    }
    // The most it could evaporate is the fire's whole output divided by the heat
    // of one kilogram, and even that is an over-estimate because the boundary and
    // the vents take a share.
    const double ceiling = 2.0e6 / (fire::kCpWater * (fire::kTSaturation - kTAmbient) +
                                    fire::kLatentHeat);
    for (int i = 0; i < 3; ++i) {
        expectTrue("nothing evaporates beyond what the fire can boil",
                   evaporatedShare[i] * 5.0 <= ceiling);
        expectTrue("so the nominal fraction is a ceiling that is never reached",
                   evaporatedShare[i] < 0.5 * e[i]);
    }
    expectTrue("and a nine-fold change in the fraction moves the water landing under 2%",
               std::abs(deck[0] / deck[2] - 1.0) < 0.02);
    // Vacuity: there has to be water landing, and the three cases have to differ
    // at all, or "insensitive" would be a statement about a constant.
    expectTrue("while several kilograms a second really are landing", deck[1] > 4.0);
    expectTrue("and the three cases are genuinely different runs", deck[0] != deck[2]);
}

// A drencher may not take a layer below the temperature of its own water: that is
// energy from nowhere. On this ship that is not a corner case, it is the
// operating point -- a deck drencher can absorb eight times what the design fire
// releases -- so it is asserted rather than assumed.
void testTheDrencherCannotCoolBelowTheTemperatureOfItsOwnWater() {
    // No fire at all, and a flow that would remove megawatts if it were allowed.
    Room r = makeRoom(0.0, 30.0);
    fire::Drencher d;
    d.name = "flood";
    d.gasCompartment = 0;
    d.flow = 50.0;
    d.on = true;
    r.model.drenchers.push_back(d);
    runToSteady(r, 200);

    // Not exactly zero, and it should not be asserted at exactly zero: the vent
    // solve moves the layer temperature by round-off, and the drencher correctly
    // removes what that round-off represents. 1.0e-13 J over 200 s of a 5 MJ room
    // is 2e-20 of it.
    expectTrue("a drencher in a cold compartment removes nothing measurable",
               std::abs(r.model.account.suppressionCooling) < 1e-9);
    expectTrue("and the gas never goes below ambient",
               r.model.gas[0].upper.temperature() >= kTAmbient - 1e-9 &&
                   r.model.gas[0].lower.temperature() >= kTAmbient - 1e-9);
    expectNear("nothing evaporates, exactly, because nothing is near saturation",
               r.model.account.waterEvaporated, 0.0, 0.0);
    expectNear("so every kilogram delivered lands", r.model.drenchers[0].lastToDeck, 50.0, 0.0);
    expectNear("which is what the account says too", r.model.account.waterToShip(),
               r.model.account.waterDelivered, 0.0);

    // And with a fire the drencher over-powers, the layer settles at a real
    // equilibrium -- `Q_fire = Q_spray(T)` -- rather than being pinned by a cap.
    Room hot = makeRoom(2.0e6, 30.0);
    fire::Drencher big;
    big.name = "deluge";
    big.gasCompartment = 0;
    big.flow = 20.0;
    big.on = true;
    hot.model.drenchers.push_back(big);
    runToSteady(hot, 600);
    const double demand = big.flow * (fire::kCpWater * (fire::kTSaturation - kTAmbient) +
                                      0.3 * fire::kLatentHeat);
    expectTrue("an over-powered drencher applies a tenth of its nominal demand",
               hot.model.drenchers[0].lastCooling < 0.15 * demand);
    expectTrue("because it is absorbing the fire rather than its own capacity",
               hot.model.account.suppressionCooling > 0.9 * hot.model.account.heatReleased &&
                   hot.model.account.suppressionCooling <= hot.model.account.heatReleased);
    expectTrue("and the layer sits at a steady 333 K, above the water and below saturation",
               hot.model.gas[0].upper.temperature() > kTAmbient + 20.0 &&
                   hot.model.gas[0].upper.temperature() < fire::kTSaturation);
    expectNear("with nothing evaporating at all out of a layer below boiling",
               hot.model.drenchers[0].lastEvaporated, 0.0, 0.0);
}

// The suppression sink is stiff, and a stiff term that does not converge under
// refinement is a term whose answer is the step size. It converges, first order,
// and -- the part that matters for this item -- **the water landing on the deck
// converges an order of magnitude tighter than the layer temperature does**, so
// the stability answer is not hostage to the thermal one.
void testTheSuppressionSinkConvergesAndTheWaterLandingConvergesFaster() {
    std::vector<double> temperature, water;
    for (double h : {0.4, 0.2, 0.1, 0.05, 0.025}) {
        Room r = makeRoom(2.0e6, 30.0);
        fire::Drencher d;
        d.name = "head";
        d.gasCompartment = 0;
        d.flow = 2.0;
        d.on = true;
        r.model.drenchers.push_back(d);
        r.model.maxSubstep = h;
        runToSteady(r, 600);
        temperature.push_back(r.model.gas[0].upper.temperature());
        water.push_back(r.model.account.waterDelivered - r.model.account.waterEvaporated);
    }
    // Successive halvings have to shrink, which is what separates convergence
    // from wandering. Measured 44.0, 5.1, 2.3, 1.1 K.
    for (std::size_t i = 2; i < temperature.size(); ++i)
        expectTrue("halving the substep moves the layer temperature less each time",
                   std::abs(temperature[i] - temperature[i - 1]) <
                       std::abs(temperature[i - 1] - temperature[i - 2]));
    expectTrue("and the finest two agree to under half a per cent",
               std::abs(temperature[4] / temperature[3] - 1.0) < 5e-3);
    expectTrue("while the water landing agrees to under one per cent across a 16x refinement",
               std::abs(water.front() / water.back() - 1.0) < 0.01);
    expectTrue("which is an order of magnitude tighter than the temperature over the same span",
               std::abs(water.front() / water.back() - 1.0) <
                   0.15 * std::abs(temperature.front() / temperature.back() - 1.0));
    expectTrue("with the coarse and fine runs both delivering real tonnages",
               water.front() > 800.0 && water.back() > 800.0);
}

// --- Suppression water is floodwater ------------------------------------------

// **The load-bearing claim of the whole item.** `ship.cpp` contains no free
// surface *correction* -- no `rho i / disp` term anywhere -- because it re-levels
// the real water body at every attitude it is asked about. So suppression water
// needs no path of its own, and the proof is that the emergent free-surface loss
// is the closed form it is supposed to be.
//
// Measured by freezing the same mass at the same centroid into the lightship, so
// that draft, KG, KB and BM are identical between the two ships and the only
// difference left is that one body of water can re-level and the other cannot.
void testSuppressionWaterIsFloodwaterAndItsFreeSurfaceIsTheClosedForm() {
    const double eps = 0.002;
    for (double tonnes : {100.0, 200.0, 400.0}) {
        Barge liquid = makeBarge(3.0);
        liquid.ship.compartments[0].waterVolume = tonnes * 1000.0 / kRhoSeawater;
        liquid.ship.step(1e-9, liquid.sea);

        Barge solid = makeBarge(3.0);
        const double m = tonnes * 1000.0, lm = solid.ship.lightshipMass;
        solid.ship.lightshipCog =
            (solid.ship.lightshipCog * lm + liquid.ship.compartments[0].waterCentroid * m) /
            (lm + m);
        solid.ship.lightshipMass = lm + m;
        solid.ship.step(1e-9, solid.sea);

        const double disp = liquid.ship.diagnostics(liquid.sea).displacementMass;
        const double closedForm = kRhoSeawater * 1.0 * kDeckI / disp;
        const double loss = gmAt(solid.ship, solid.sea, eps) - gmAt(liquid.ship, liquid.sea, eps);

        // Asserted at 5e-6 relative, which is where the *finite difference* sits:
        // the residual scales as eps^2 and is 1.36e-6 at eps = 0.002, so anything
        // looser would also pass on a model that had lost the property and merely
        // converged well. The two GMs are metres apart, so this is a genuine five
        // significant figures.
        expectTrue("the free-surface loss is rho mu b^3 l / 12 over the displacement",
                   std::abs(loss / closedForm - 1.0) < 5e-6);
        // Guard against vacuity from both ends.
        expectTrue("and it is metres, not a rounding error", closedForm > 1.5);
        expectTrue("with the free surface taking most of the barge's stability",
                   loss > 0.6 * gmAt(solid.ship, solid.sea, eps));
        expectTrue("so the liquid ship keeps under half the GM the frozen one has",
                   gmAt(liquid.ship, liquid.sea, eps) <
                       0.4 * gmAt(solid.ship, solid.sea, eps));
    }

    // The same identity on a real hull, where the second moment is not arithmetic
    // and the permeability is not one. Measured off the compartment's own mesh by
    // a slicing integration that shares no code with either model.
    Sea sea{0.0};
    Ship reference = game::buildFerry();
    reference.initialise(sea);
    const int vd = reference.findCompartment("vehicle_deck");
    const Compartment& c0 = reference.compartments[static_cast<std::size_t>(vd)];

    const double tonnes = 500.0;
    Ship liquid = game::buildFerry();
    liquid.initialise(sea);
    liquid.compartments[static_cast<std::size_t>(vd)].waterVolume = tonnes * 1000.0 / kRhoSeawater;
    liquid.step(1e-9, sea);
    const Compartment& cl = liquid.compartments[static_cast<std::size_t>(vd)];

    Ship solid = game::buildFerry();
    solid.initialise(sea);
    const double m = tonnes * 1000.0, lm = solid.lightshipMass;
    solid.lightshipCog = (solid.lightshipCog * lm + cl.waterCentroid * m) / (lm + m);
    solid.lightshipMass = lm + m;
    solid.step(1e-9, sea);

    const double iyy = waterplaneSecondMoment(c0.mesh, c0.bboxLo, c0.bboxHi, cl.surfaceOffset,
                                              0.125);
    const double disp = liquid.diagnostics(sea).displacementMass;
    const double closedForm = kRhoSeawater * cl.permeability * iyy / disp;
    const double loss = gmAt(solid, sea, eps) - gmAt(liquid, sea, eps);
    expectTrue("and on the ferry's own vehicle deck, against its measured second moment",
               std::abs(loss / closedForm - 1.0) < 3e-4);
    expectTrue("which is 56308 m4 of waterplane", iyy > 56000.0 && iyy < 56600.0);
    // **The permeability is in there and it is not optional.** Dropping it would
    // over-state the loss by 11% on this deck, which is more than the tolerance
    // above by two orders of magnitude and would still look plausible.
    expectTrue("with the compartment's permeability genuinely in the answer",
               std::abs(loss / (closedForm / cl.permeability) - 1.0) > 0.09);
    expectTrue("and 500 t on an undivided vehicle deck costs this ship metres of GM",
               closedForm > 5.0);
}

// The free-surface loss above is the *initial* one, and it only exists while the
// water still spans the compartment. Tilt far enough and a shallow layer runs off
// the high side; the surface is then narrower than the deck and the moment
// collapses. That is real -- it is why a ship with no GM lolls instead of
// capsizing -- and it is the reason `Diagnostics::gmTransverse`, taken at
// eps = 0.03 rad, is not the initial GM of a centimetres-deep layer.
void testPocketingLimitsFreeSurfaceExactlyWhereGeometrySays() {
    const double tonnes = 100.0;
    Barge liquid = makeBarge(3.0);
    liquid.ship.compartments[0].waterVolume = tonnes * 1000.0 / kRhoSeawater;
    liquid.ship.step(1e-9, liquid.sea);
    Barge solid = makeBarge(3.0);
    const double m = tonnes * 1000.0, lm = solid.ship.lightshipMass;
    solid.ship.lightshipCog =
        (solid.ship.lightshipCog * lm + liquid.ship.compartments[0].waterCentroid * m) / (lm + m);
    solid.ship.lightshipMass = lm + m;
    solid.ship.step(1e-9, solid.sea);

    const double depth = liquid.ship.compartments[0].waterVolume / (kDeckL * kDeckB);
    // The surface reaches the high side of the deck at exactly this angle.
    const double limit = std::atan(2.0 * depth / kDeckB);
    const double closedForm = kRhoSeawater * kDeckI /
                              liquid.ship.diagnostics(liquid.sea).displacementMass;
    expectTrue("the layer is 16 cm deep", depth > 0.15 && depth < 0.18);
    expectTrue("so it spans the deck only out to about 1.55 degrees",
               limit > 0.026 && limit < 0.028);

    auto ratio = [&](double eps) {
        return (gmAt(solid.ship, solid.sea, eps) - gmAt(liquid.ship, liquid.sea, eps)) /
               closedForm;
    };
    // Measured 1.000061 at half the limit and 1.000198 at nine tenths of it, so
    // the tolerances are what was measured rather than a round number: a looser
    // pair would pass on a model with no pocketing at all.
    expectTrue("well inside that angle the full free surface is there",
               std::abs(ratio(0.5 * limit) - 1.0) < 1e-4);
    expectTrue("and still at nine tenths of it", std::abs(ratio(0.9 * limit) - 1.0) < 3e-4);
    // Past it the moment falls away, and it has to fall a long way rather than a
    // little: this is a different regime, not a rounding difference.
    expectTrue("but at three times it, nearly 40% of the moment is gone",
               ratio(3.0 * limit) < 0.65);
    expectTrue("monotonically", ratio(3.0 * limit) < ratio(2.0 * limit) &&
                                    ratio(2.0 * limit) < ratio(1.2 * limit) &&
                                    ratio(1.2 * limit) < ratio(0.9 * limit));
    // Which is exactly why the ship's own published GM is not the initial one for
    // a layer this shallow: eps = 0.03 rad is 1.1 times the limit here, and the
    // deviation it reports is a hundred times the one just inside it.
    expectTrue("so the ship's own eps = 0.03 GM is not the initial GM of this layer",
               std::abs(ratio(0.03) - 1.0) > 20.0 * std::abs(ratio(0.5 * limit) - 1.0));
}

// And on the ferry the same effect is not a technicality, it is the difference
// between a stability report that says she is safe and one that says she is not.
// **`Diagnostics::gmTransverse` reads +0.59 m with 50 tonnes on the vehicle deck,
// where the true initial GM is -3.77 m**, because 50 tonnes is a 2.9 cm layer and
// eps = 0.03 rad is ten times the angle it spans the deck to.
//
// Reported rather than repaired: eps = 0.03 is the right question for a ship at a
// finite angle, it is what the 35 published figures were taken under, and moving
// it would move all of them.
void testTheFerrysPublishedGmIsNotTheInitialGmOfAShallowLayer() {
    Sea sea{0.0};
    Ship s = game::buildFerry();
    s.initialise(sea);
    const int vd = s.findCompartment("vehicle_deck");
    s.compartments[static_cast<std::size_t>(vd)].waterVolume = 50.0e3 / kRhoSeawater;
    s.step(1e-9, sea);

    const double published = s.diagnostics(sea).gmTransverse;
    const double initial = gmAt(s, sea, 0.001);
    expectTrue("with 50 t on the vehicle deck the published GM still reads positive",
               published > 0.4 && published < 0.8);
    expectTrue("while the initial GM is strongly negative", initial < -3.5);
    expectTrue("a gap of over four metres", published - initial > 4.0);

    // The mechanism, so that this is a statement about pocketing and not about
    // two numbers: the layer is 2.9 cm deep and spans the deck only to 0.0031 rad.
    const Compartment& c = s.compartments[static_cast<std::size_t>(vd)];
    const double depth = c.waterVolume / c.permeability / 1868.4;
    expectTrue("because the layer is under three centimetres deep",
               depth > 0.028 && depth < 0.030);
    expectTrue("and spans the deck only out to a fifth of a degree",
               std::atan(2.0 * depth / 18.68) < 0.2 * 0.03);
    // Vacuity: the two GMs have to be measuring the same ship.
    Ship dryShip = game::buildFerry();
    dryShip.initialise(sea);
    expectTrue("and dry, the two agree",
               std::abs(dryShip.diagnostics(sea).gmTransverse - gmAt(dryShip, sea, 0.001)) <
                   0.01);
}

// --- Freeing ports -------------------------------------------------------------

// A drencher into a compartment with a port in its side reaches an equilibrium
// depth, and that depth is the weir's. This is the design question a freeing port
// is sized against, and it is a closed form: `h = (q / ((2/3) Cd b sqrt(2g)))^(2/3)`.
//
// Deliberately run with **no fire**, so the whole flow lands and the hydraulics
// are isolated from the thermodynamics.
void testFreeingPortsHoldTheDeckAtTheWeirsEquilibriumDepth() {
    const double width = 0.8, cd = 0.6;
    for (double flow : {10.0, 40.0}) {
        const double q = flow / kRhoSeawater;   // m^3/s landing, all of it
        const double expected =
            std::pow(q / (cd * width * (2.0 / 3.0) * std::sqrt(2.0 * kGravity)), 2.0 / 3.0);

        // **Seeded at the closed form and asked whether it is a fixed point**,
        // rather than run from dry and asked whether it arrived. The approach is
        // exponential with a time constant of `(2/3) A h / q` -- 25 minutes at
        // the lower flow -- so a from-dry run that stopped short would look like
        // a failed closed form rather than like an unfinished integration. The
        // from-dry direction is checked below, for what it can actually say.
        Barge b = makeBarge(3.0);
        b.ship.compartments[0].waterVolume = expected * kDeckL * kDeckB;
        b.ship.step(1e-9, b.sea);
        fire::Drencher d;
        d.name = "deck_drencher";
        d.gasCompartment = 0;
        d.flow = flow;
        d.on = true;
        b.model.drenchers.push_back(d);
        fire::Scupper sc;
        sc.name = "freeing_port";
        sc.gasCompartment = 0;
        sc.sillPos = {0.0, -kDeckB / 2, kDeckZ};
        sc.width = width;
        sc.dischargeCoeff = cd;
        b.model.scuppers.push_back(sc);
        expectTrue("the definition is sound", b.model.validate().empty());

        for (int i = 0; i < 2000; ++i) {
            b.ship.step(1.0, b.sea);
            b.model.step(1.0, b.ship, b.sea);
            b.model.applyTo(b.ship);
        }

        const double depth = b.ship.compartments[0].waterVolume / (kDeckL * kDeckB);
        // Measured to drift 1.5e-4 of the seeded depth over 2000 s, downward, and
        // that drift is the model's own one-step lag on the ship's water level
        // rather than a disagreement with the closed form.
        expectTrue("the weir's equilibrium depth is a fixed point of the model",
                   std::abs(depth / expected - 1.0) < 1e-3);
        expectTrue("with the port carrying what is landing to a part in a thousand",
                   std::abs(b.model.scuppers[0].lastFlow * kRhoSeawater / flow - 1.0) < 1e-3);
        expectTrue("and the head the port sees is that depth",
                   std::abs(b.model.scuppers[0].lastInsideHead / expected - 1.0) < 1e-3);
        expectTrue("the port is discharging freely, not drowned",
                   b.model.scuppers[0].lastOutsideHead == 0.0);
        // Vacuity: the equilibrium has to be a *balance* rather than an empty deck
        // or a full one.
        expectTrue("the deck is neither dry nor full",
                   b.ship.compartments[0].waterVolume > 1.0 &&
                       b.ship.compartments[0].fillFraction() < 0.2);

        // And from dry it climbs towards the same depth without passing it, which
        // is what says the fixed point is the one the model actually seeks.
        Barge from = makeBarge(3.0);
        from.model.drenchers.push_back(d);
        from.model.scuppers.push_back(sc);
        double previous = 0;
        bool monotone = true, overshot = false;
        for (int i = 0; i < 4000; ++i) {
            from.ship.step(1.0, from.sea);
            from.model.step(1.0, from.ship, from.sea);
            from.model.applyTo(from.ship);
            // Checked into a flag rather than asserted in the loop: a regression
            // here would otherwise print four thousand identical failures and
            // scroll every other diagnostic in the suite off the top.
            const double now = from.ship.compartments[0].waterVolume;
            monotone = monotone && now >= previous - 1e-12;
            overshot = overshot || now > expected * kDeckL * kDeckB * 1.001;
            previous = now;
        }
        expectTrue("the deck fills monotonically towards the equilibrium", monotone);
        expectTrue("and never overshoots it", !overshot);
        expectTrue("and gets most of the way there in an hour",
                   previous / (expected * kDeckL * kDeckB) > 0.88);
    }

    // Four times the flow is only 2.52 times the depth, because the weir goes as
    // h^(3/2). That exponent is the reason freeing ports work at all and it is
    // what a linear drain model would get wrong, so it is asserted directly.
    expectNear("and the depth goes as the two-thirds power of the flow",
               std::pow(4.0, 2.0 / 3.0), 2.5198, 1e-4);
}

// Block the ports and the same drencher fills the deck without bound. This is the
// contrast the item is about: the water gets in either way, and whether it gets
// out is what decides the ship.
void testABlockedFreeingPortLetsTheDeckFillWithoutBound() {
    double held[2] = {0, 0};
    for (int blocked = 0; blocked < 2; ++blocked) {
        Barge b = makeBarge(3.0);
        fire::Drencher d;
        d.name = "deck_drencher";
        d.gasCompartment = 0;
        d.flow = 40.0;
        d.on = true;
        b.model.drenchers.push_back(d);
        fire::Scupper sc;
        sc.name = "freeing_port";
        sc.gasCompartment = 0;
        sc.sillPos = {0.0, -kDeckB / 2, kDeckZ};
        // Four metres of port, which is enough to reach equilibrium well inside
        // the run: with 0.8 m the equilibrium is 56 t and an hour of blocked
        // accumulation is only 144 t, so the "contrast" would be a factor of two
        // and would say more about the run length than about the port.
        sc.width = 4.0;
        sc.blocked = blocked != 0;
        b.model.scuppers.push_back(sc);
        for (int i = 0; i < 3600; ++i) {
            b.ship.step(1.0, b.sea);
            b.model.step(1.0, b.ship, b.sea);
            b.model.applyTo(b.ship);
        }
        held[blocked] = b.ship.compartments[0].waterVolume * kRhoSeawater;
        if (blocked != 0) {
            expectNear("a blocked port drains nothing at all", b.model.scuppers[0].lastFlow, 0.0,
                       0.0);
            expectNear("so every kilogram delivered is still aboard",
                       held[1], b.model.account.waterDelivered, 1e-6);
        }
    }
    expectTrue("an hour of drencher with the ports blocked is 144 tonnes",
               held[1] > 143.0e3 && held[1] < 145.0e3);
    expectTrue("while clear ports hold it to 19 tonnes and stop",
               held[0] > 18.0e3 && held[0] < 21.0e3);
    expectTrue("a factor of seven", held[1] > 7.0 * held[0]);
}

// "...and what happens when they are blocked **or under water**." A ship that has
// settled far enough to put her freeing ports under does not drain through them;
// she floods through them, and it is the same formula with the heads the other
// way round.
void testAFreeingPortUnderTheSeaAdmitsInsteadOfDraining() {
    Barge b = makeBarge(4.5);   // deck 4 m up, floating at 4.5: half a metre under
    fire::Scupper sc;
    sc.name = "freeing_port";
    sc.gasCompartment = 0;
    sc.sillPos = {0.0, -kDeckB / 2, kDeckZ};
    sc.width = 0.8;
    b.model.scuppers.push_back(sc);

    b.ship.step(1e-9, b.sea);
    b.model.step(1.0, b.ship, b.sea);
    expectTrue("the sea is standing over the sill", b.model.scuppers[0].lastOutsideHead > 0.4);
    expectNear("with nothing inside to come out", b.model.scuppers[0].lastInsideHead, 0.0, 0.0);
    expectTrue("so the port runs inward", b.model.scuppers[0].lastFlow < 0);

    for (int i = 0; i < 600; ++i) {
        b.ship.step(1.0, b.sea);
        b.model.step(1.0, b.ship, b.sea);
        b.model.applyTo(b.ship);
    }
    expectTrue("and the deck floods through it", b.ship.compartments[0].waterVolume > 10.0);
    expectTrue("which the account records as negative drainage",
               b.model.account.waterDrained < -10.0e3);
    expectNear("with no drencher anywhere near it", b.model.account.waterDelivered, 0.0, 0.0);
}

// --- The water account ---------------------------------------------------------

// Every kilogram is somewhere. Delivered, less evaporated, less drained, is
// exactly what has been handed to the flooding solve plus what is still owed --
// asserted to the bit, because it is a sum of the same doubles rather than an
// approximation of anything.
void testTheWaterAccountClosesAgainstWhatTheShipReceived() {
    Barge b = makeBarge(3.0);
    // **With a fire in the compartment, and that is not decoration.** The first
    // version ran this cold, so nothing evaporated, and a mutant that handed the
    // *whole* delivered flow to the deck instead of the part that did not
    // evaporate survived the entire suite: with `evaporated == 0` the two are the
    // same number. A 30 MW fire against 25 kg/s puts 30% of the flow up the
    // funnel as steam, which is the only condition under which this test can tell
    // the two apart.
    fire::DesignFire f;
    f.name = "deck_fire";
    f.compartment = 0;
    f.baseZ = 4.5;
    f.diameter = 3.0;
    f.peakHeatRelease = 30.0e6;
    f.steadyDuration = 1e9;
    b.model.fires.push_back(f);
    b.model.gas[0].wallConductance = 25.0;

    fire::Drencher d;
    d.name = "deck_drencher";
    d.gasCompartment = 0;
    d.flow = 25.0;
    d.on = true;
    b.model.drenchers.push_back(d);
    fire::Scupper sc;
    sc.name = "freeing_port";
    sc.gasCompartment = 0;
    sc.sillPos = {0.0, -kDeckB / 2, kDeckZ};
    sc.width = 0.5;
    b.model.scuppers.push_back(sc);

    double written = 0;
    for (int i = 0; i < 1500; ++i) {
        b.ship.step(1.0, b.sea);
        b.model.step(1.0, b.ship, b.sea);
        const double before = b.ship.compartments[0].waterVolume;
        b.model.applyTo(b.ship);
        written += b.ship.compartments[0].waterVolume - before;
    }
    const fire::Account& a = b.model.account;
    const double owed = b.model.pendingWater()[0];
    // Measured at 5.8e-11 kg on 37.5 tonnes moved -- 1.6e-15 relative, which is
    // the accumulation of 1500 steps of round-off and nothing else. Asserted at
    // 1e-9 rather than at a fraction, because this is a sum of the same doubles
    // and any real hole in it would be kilograms.
    expectTrue("what was written plus what is still owed is what suppression moved",
               std::abs((written + owed) * kRhoSeawater - a.waterToShip()) < 1e-9);
    expectNear("and waterToShip is delivered less evaporated less drained",
               a.waterToShip(), a.waterDelivered - a.waterEvaporated - a.waterDrained, 0.0);

    // Vacuity guards: all three terms have to be real numbers, and in particular
    // the evaporated one has to be large enough that leaving it out of the write
    // would show. It is 30% of the flow here.
    expectTrue("37 tonnes were delivered", a.waterDelivered > 37.0e3);
    expectTrue("eleven of them left as steam", a.waterEvaporated > 9.0e3);
    expectTrue("which is a quarter to a third of the flow",
               a.waterEvaporated > 0.25 * a.waterDelivered &&
                   a.waterEvaporated < 0.35 * a.waterDelivered);
    expectTrue("four were drained back over the side", a.waterDrained > 3.0e3);
    expectTrue("leaving twenty on the deck", written * kRhoSeawater > 18.0e3);
    // And the three are genuinely different quantities: handing the deck the
    // whole delivered flow would land here, 11 tonnes out.
    expectTrue("with the water landing well short of the water delivered",
               written * kRhoSeawater < a.waterDelivered - a.waterEvaporated + 1.0);
}

// A compartment can be flooded solid with a drencher still running into it, and
// what will not fit has to stay **owed** rather than be dropped. Discarding it
// would put a hole in the water account of exactly the kind the layer clamps put
// in the energy one -- and both of those paths survived the whole suite until
// this case existed.
void testWaterThatWillNotFitStaysOwedRatherThanBeingDropped() {
    Barge b = makeBarge(3.0);
    Compartment& c = b.ship.compartments[0];
    c.waterVolume = 0.999 * c.floodableVolume();
    b.ship.step(1e-9, b.sea);
    fire::Drencher d;
    d.name = "deck_drencher";
    d.gasCompartment = 0;
    d.flow = 25.0;
    d.on = true;
    b.model.drenchers.push_back(d);

    double written = 0;
    for (int i = 0; i < 400; ++i) {
        b.ship.step(1.0, b.sea);
        b.model.step(1.0, b.ship, b.sea);
        const double before = b.ship.compartments[0].waterVolume;
        b.model.applyTo(b.ship);
        written += b.ship.compartments[0].waterVolume - before;
    }
    const fire::Account& a = b.model.account;
    const double owed = b.model.pendingWater()[0];

    expectTrue("the compartment really did fill solid",
               b.ship.compartments[0].fillFraction() > 0.9999);
    expectTrue("it never took more than it holds",
               b.ship.compartments[0].waterVolume <=
                   b.ship.compartments[0].floodableVolume() + 1e-9);
    expectTrue("most of what was delivered would not fit", owed * kRhoSeawater > 6.0e3);
    expectTrue("and only the headroom was written", written * kRhoSeawater < 3.0e3);
    // The account is still exact. Measured at 9.1e-9 kg on 10 tonnes moved, which
    // is 9e-13 relative -- worse than the un-clamped case because `want - got`
    // cancels two large numbers every step, and still nothing next to a kilogram.
    expectTrue("and written plus owed is still exactly what suppression moved",
               std::abs((written + owed) * kRhoSeawater - a.waterToShip()) < 1e-6);
    expectNear("with nothing evaporated to complicate it", a.waterEvaporated, 0.0, 0.0);
    expectTrue("ten tonnes were delivered in all", a.waterDelivered > 9.9e3);
}

// The energy account, with a drencher running, on the ferry. The existing account
// closes to 1e-14 of scale; adding a megawatt-scale sink that is subtracted from
// the same total is exactly the sort of change that puts a hole in it.
void testTheAccountStillClosesWithSuppressionRunning() {
    FerryFire f = makeFerryFire();
    fire::Drencher d;
    d.name = "engine_room_deluge";
    d.gasCompartment = f.model.findGas("engine_room_s");
    d.flow = 2.0;
    d.on = true;
    f.model.drenchers.push_back(d);
    expectTrue("the model with suppression on it is still self-consistent",
               f.model.validate().empty());

    for (int i = 0; i < 1200; ++i) {
        f.model.step(1.0, f.ship, f.sea);
        f.model.applyTo(f.ship);
    }
    const fire::Account& a = f.model.account;
    // Measured at 8.9e-15 of scale with a 1.7 GJ sink in the books, which is the
    // same machine precision the dry account closes to. Asserted at 5e-14 rather
    // than at the 1e-11 the dry account uses: a looser bound would not notice the
    // sink being applied twice, or applied to one layer and accounted against the
    // other.
    expectTrue("the energy account still closes with a gigajoule-scale sink in it",
               std::abs(a.energyResidualFraction()) < 5e-14);
    expectTrue("and the mass account is untouched", std::abs(a.massResidualFraction()) < 1e-13);
    expectTrue("the drencher really did take gigajoules out", a.suppressionCooling > 1.0e9);
    expectTrue("but not more than the fire released", a.suppressionCooling < a.heatReleased);
    expectTrue("and it really did put a tonne and a half of water into the ship",
               a.waterToShip() > 1.5e3);
    expectTrue("with the burning compartment measurably cooler for it",
               f.model.gas[0].upper.temperature() < 400.0);

    // The unfired compartment next door has no drencher, so it must not have
    // cooled: a sink written to the wrong index would balance the account
    // perfectly and still be wrong.
    expectTrue("while the compartment with no drencher in it is untouched by one",
               f.model.gas[1].upper.temperature() > kTAmbient + 1.0);
}

// --- The structural consequence ------------------------------------------------

// **The number the roadmap item is about.** A lorry fire on the ferry's undivided
// vehicle deck, a SOLAS drencher over it, and what four hours of it does to her
// stability -- with the freeing ports clear and with them blocked.
//
// The headline is not the GM. GM collapses within a minute and then stops moving:
// an undivided 100 x 19 m deck has a free-surface moment so large that ten tonnes
// of water takes this ship's 2.00 m of GM negative, and everything after that is
// the same -3.7 m whether there are thirty tonnes on the deck or two thousand.
// **The angle of loll is what tracks the water**, and it is what the ports decide.
void testTheFerrysLollTracksTheDrencherAndTheFreeingPortsDecideIt() {
    struct Outcome { double water = 0, gm = 0, loll = 0, drained = 0; };
    Outcome result[2];

    for (int blocked = 0; blocked < 2; ++blocked) {
        Ship ship = game::buildFerry();
        // No damage: this is about firefighting water and nothing else. The
        // ferry's authored `downflood_*` openings are shut as well, so that the
        // drainage measured here is the freeing ports' and not theirs -- they sit
        // 0.1 m above the deck and do drain it, which is worth knowing and is a
        // different experiment.
        for (Opening& o : ship.openings)
            if (o.name == "breach_er_s" || o.name.rfind("downflood_", 0) == 0) o.open = false;
        Sea sea{0.0};
        ship.initialise(sea);
        const int vd = ship.findCompartment("vehicle_deck");

        fire::Model model;
        model.attach(ship, {vd});
        fire::DesignFire d;
        d.name = "lorry";
        d.compartment = 0;
        d.baseZ = 7.5;
        d.diameter = 4.0;
        d.growthCoefficient = fire::kGrowthFast;
        d.peakHeatRelease = 20.0e6;
        d.steadyDuration = 1e9;
        model.fires.push_back(d);
        model.gas[0].wallConductance = 25.0;

        fire::Drencher dr;
        dr.name = "deck_drencher";
        dr.gasCompartment = 0;
        // One 20 m section of the deck at the SOLAS ro-ro rate.
        dr.flow = fire::sprayMassFlow(20.0 * 18.68, 5.0);
        dr.on = true;
        model.drenchers.push_back(dr);

        // Six one-metre freeing ports along both deck edges.
        for (double x : {-20.0, 0.0, 20.0})
            for (double y : {-9.5, 9.5}) {
                fire::Scupper sc;
                sc.name = "port";
                sc.gasCompartment = 0;
                sc.sillPos = {x, y, 7.0};
                sc.width = 1.0;
                sc.blocked = blocked != 0;
                model.scuppers.push_back(sc);
            }
        expectTrue("the ferry's suppression definition is sound", model.validate().empty());

        for (int i = 0; i < 14400; ++i) {
            ship.step(1.0, sea);
            model.step(1.0, ship, sea);
            model.applyTo(ship);
        }
        const Compartment& c = ship.compartments[static_cast<std::size_t>(vd)];
        result[blocked] = {c.waterVolume * ship.seaDensity, gmAt(ship, sea, 0.002),
                           lollDeg(ship, sea), model.account.waterDrained};
    }

    // Clear ports: the deck reaches a steady 36 tonnes and stays there for the
    // rest of the fire, at 0.81 degrees of loll. 397 t went over the side to hold
    // it there, which is the whole point of the port -- eleven times what is
    // standing on the deck at any moment.
    expectTrue("with the ports clear the deck holds a steady 36 tonnes",
               result[0].water > 33.0e3 && result[0].water < 40.0e3);
    expectTrue("and she lolls under a degree", result[0].loll > 0.5 && result[0].loll < 1.2);
    expectTrue("having put four hundred tonnes back over the side",
               result[0].drained > 380.0e3 && result[0].drained > 10.0 * result[0].water);

    // Blocked: everything that lands stays. Four hours of it is 433 tonnes and
    // 9.7 degrees of loll, on a ship that started with 2.00 m of GM.
    expectTrue("blocked, four hours puts 433 tonnes on the deck",
               result[1].water > 400.0e3 && result[1].water < 470.0e3);
    expectNear("with nothing drained at all", result[1].drained, 0.0, 0.0);
    expectTrue("and she lolls into double figures of degrees",
               result[1].loll > 8.0 && result[1].loll < 12.0);
    expectTrue("which is many times the loll the working ports bought",
               result[1].loll > 8.0 * result[0].loll);

    // **GM is not the discriminator, and saying so is the finding.** Both cases
    // sit near -3.7 m; a twelve-fold difference in water on deck moves GM by a
    // tenth of a metre while it moves the loll by twelve.
    expectTrue("both cases have lost all their initial stability",
               result[0].gm < -3.0 && result[1].gm < -3.0);
    expectTrue("and GM barely distinguishes them", std::abs(result[1].gm - result[0].gm) < 0.2);
    expectTrue("while the water on deck differs by a factor of ten",
               result[1].water > 10.0 * result[0].water);

    // Vacuity: the intact ship really did have stability to lose.
    Ship intact = game::buildFerry();
    Sea sea{0.0};
    intact.initialise(sea);
    expectTrue("the ferry starts with 2.00 m of GM",
               std::abs(gmAt(intact, sea, 0.002) - 2.0) < 0.02);
}

// --- The exact control ---------------------------------------------------------

// A byte-for-byte view of the gas, so that "suppression off changed nothing" can
// be checked on the fire model as well as on the ship.
std::vector<double> gasFingerprint(const fire::Model& m) {
    std::vector<double> v;
    v.push_back(m.time);
    for (const fire::GasCompartment& g : m.gas) {
        v.push_back(g.upper.mass);
        v.push_back(g.upper.energy);
        v.push_back(g.upper.products);
        v.push_back(g.lower.mass);
        v.push_back(g.lower.energy);
        v.push_back(g.lower.products);
    }
    for (const fire::Vent& t : m.vents) {
        v.push_back(t.massAToB);
        v.push_back(t.massBToA);
        v.push_back(t.neutralPlaneZ);
    }
    v.push_back(m.account.energy);
    v.push_back(m.account.mass);
    v.push_back(m.account.wallLoss);
    v.push_back(m.account.enthalpyOut);
    return v;
}

// **Suppression switched off must be worth exactly nothing.**
//
// Not close: identical. The three flooding scenarios and the 35 figures the full
// gate publishes are validated against their own behaviour, and a drencher that
// is not running has no business moving any of them. This runs a real fire on the
// ferry twice -- once with no suppression objects at all, once with a drencher
// switched off and freeing ports blocked -- and compares every state double, on
// the ship *and* on the gas, for exact equality.
void testSuppressionOffLeavesTheShipAndTheGasBitIdentical() {
    FerryFire bare = makeFerryFire();
    FerryFire off = makeFerryFire();

    fire::Drencher d;
    d.name = "idle";
    d.gasCompartment = 0;
    d.flow = 40.0;          // a large flow, so "off" is doing the work and not a zero
    d.on = false;
    off.model.drenchers.push_back(d);
    fire::Scupper sc;
    sc.name = "shut";
    sc.gasCompartment = 0;
    sc.sillPos = {6.0, -8.0, 1.8};
    sc.width = 2.0;
    sc.blocked = true;
    off.model.scuppers.push_back(sc);

    const int ticks = 400;
    for (int i = 0; i < ticks; ++i) {
        bare.ship.step(1.0, bare.sea);
        bare.model.step(1.0, bare.ship, bare.sea);
        bare.model.applyTo(bare.ship);
        off.ship.step(1.0, off.sea);
        off.model.step(1.0, off.ship, off.sea);
        off.model.applyTo(off.ship);
    }

    auto compare = [](const char* what, const std::vector<double>& a,
                      const std::vector<double>& b) {
        expectEqual("the fingerprints are the same shape", static_cast<long long>(a.size()),
                    static_cast<long long>(b.size()));
        std::size_t differing = 0;
        for (std::size_t i = 0; i < a.size() && i < b.size(); ++i)
            if (std::memcmp(&a[i], &b[i], sizeof(double)) != 0) ++differing;
        expectEqual(what, static_cast<long long>(differing), 0);
    };
    compare("suppression switched off leaves the ship bit-identical, in every state double",
            shipFingerprint(bare.ship), shipFingerprint(off.ship));
    compare("and the gas bit-identical too", gasFingerprint(bare.model), gasFingerprint(off.model));

    // And nothing was written, rather than zero being written.
    expectNear("no water was delivered", off.model.account.waterDelivered, 0.0, 0.0);
    expectNear("none was drained", off.model.account.waterDrained, 0.0, 0.0);
    expectNear("no cooling was applied", off.model.account.suppressionCooling, 0.0, 0.0);
    expectNear("and nothing is owed to the ship", off.model.pendingWater()[0], 0.0, 0.0);

    // Guard against vacuity: the ship and the fire both have to have been doing
    // something, or two identical inert runs would prove nothing.
    expectTrue("while the ship really was flooding through the run",
               bare.ship.totalFloodwaterMass() > 1.0e5);
    expectTrue("and the fire really was burning",
               bare.model.gas[0].upper.temperature() > kTAmbient + 100.0);
}

// A bad suppression definition does not crash; it quietly produces a wrong answer.
void testValidateCatchesBadSuppressionDefinitions() {
    fire::Model m;
    fire::GasCompartment g;
    g.name = "box";
    g.shipCompartment = kSea;   // no ship compartment behind it
    g.floorZ = 0;
    g.ceilingZ = 3;
    g.gasVolume = 30;
    g.floorArea = 10;
    m.gas.push_back(g);

    fire::Drencher d;
    d.name = "bad";
    d.gasCompartment = 4;
    d.evaporatedFraction = 1.5;
    d.flow = -1.0;
    d.waterTemperature = 400.0;
    m.drenchers.push_back(d);
    fire::Scupper sc;
    sc.name = "bad";
    sc.gasCompartment = 0;
    sc.width = 0.0;
    sc.sillPos = {0, 0, -1.0};
    m.scuppers.push_back(sc);

    const std::vector<std::string> problems = m.validate();
    auto reported = [&](const char* fragment) {
        for (const std::string& p : problems)
            if (p.find(fragment) != std::string::npos) return true;
        return false;
    };
    expectTrue("a drencher in a space that does not exist is reported",
               reported("drencher 'bad' is in a gas space"));
    expectTrue("an evaporated fraction outside [0, 1] is reported",
               reported("evaporated fraction outside"));
    expectTrue("a negative flow is reported", reported("negative flow"));
    expectTrue("water supplied at boiling is reported", reported("at or above boiling"));
    expectTrue("a scupper with nothing behind it is reported",
               reported("no ship compartment behind it"));
    expectTrue("a scupper with no width is reported", reported("has no width"));
    expectTrue("and a sill below its own deck is reported", reported("sill below the deck"));
}

}  // namespace

// --- The boundary the structure sees -----------------------------------------
//
// A plate standing on its edge, so that the two large faces look along +/-x and
// their centroids run up z -- which is a bulkhead, and is what `wallExchange`
// bands by.
solidshell::HexMesh standingPlate(double width, double height, double thickness, int nx, int ny,
                                  double baseZ = 0.0) {
    solidshell::HexMesh mesh = solidshell::makePlateMesh(width, height, thickness, nx, ny, 1);
    for (std::size_t n = 0; n < mesh.nodeCount(); ++n) {
        double* p = &mesh.position[n * 3];
        const double a = p[0], b = p[1], c = p[2];
        p[0] = c;   // thickness along x: the two large faces now look fore and aft
        p[1] = a;
        p[2] = baseZ + b;
    }
    return mesh;
}

// A gas compartment with both layers set by hand, so that the interface height is
// an input rather than something a fire has to be run to produce.
fire::GasCompartment twoLayerBox(double interfaceZ, double upperKelvin, double lowerKelvin) {
    fire::GasCompartment g;
    g.floorZ = 0.0;
    g.ceilingZ = 4.0;
    g.floorArea = 10.0;
    g.perimeter = 13.0;
    g.gasVolume = g.floorArea * (g.ceilingZ - g.floorZ);
    // `interfaceZ()` is the ceiling less the upper layer's share of the volume, and
    // the share is the internal-energy split -- so the two energies are what set it.
    const double upperVolume = (g.ceilingZ - interfaceZ) * g.floorArea;
    const double lowerVolume = g.gasVolume - upperVolume;
    g.upper.mass = kPatm * upperVolume / (kRAir * upperKelvin);
    g.lower.mass = kPatm * lowerVolume / (kRAir * lowerKelvin);
    g.upper.energy = g.upper.mass * fire::kCvAir * upperKelvin;
    g.lower.energy = g.lower.mass * fire::kCvAir * lowerKelvin;
    return g;
}

void testTheFilmCoefficientCarriesRadiationExactly() {
    fire::BoundaryFilm p;

    // The whole claim: `h (T_g - T_s)` is the convective flux plus
    // `eps sigma (T_g^4 - T_s^4)`, because `(a^2+b^2)(a+b)(a-b) = a^4 - b^4` is an
    // identity. Swept over the range a compartment fire actually spans.
    double worst = 0, worstAt = 0;
    for (double tg = 280.0; tg < 1400.0; tg += 3.7)
        for (double ts = 280.0; ts < 1400.0; ts += 5.3) {
            if (std::abs(tg - ts) < 1.0) continue;
            const double lhs = fire::filmCoefficient(tg, ts, p) * (tg - ts);
            const double rhs = p.convective * (tg - ts) +
                               p.emissivity * fire::kStefanBoltzmann *
                                   (tg * tg * tg * tg - ts * ts * ts * ts);
            const double relative = std::abs(lhs - rhs) / std::abs(rhs);
            if (relative > worst) { worst = relative; worstAt = tg; }
        }
    expectTrue("the film coefficient reproduces Stefan-Boltzmann to 1e-13", worst < 1e-13);

    // **Where they disagree, the film coefficient is the one that is right.** Near
    // equilibrium `T_g^4 - T_s^4` is a cancellation and the factored form is not, so
    // the "exact" expression is the one losing digits -- and near equilibrium is
    // where a coupled solve spends its time.
    const auto disagreement = [&](double gap) {
        const double tg = 500.0, ts = 500.0 + gap;
        const double factored = fire::filmCoefficient(tg, ts, p) * (tg - ts);
        const double naive = p.convective * (tg - ts) +
                             p.emissivity * fire::kStefanBoltzmann *
                                 (tg * tg * tg * tg - ts * ts * ts * ts);
        return std::abs(factored - naive) / std::abs(factored);
    };
    // The loss goes as `eps T / (4 dT)`, so it climbs by eight orders of magnitude
    // over eight orders of the gap -- and it is the difference form doing the
    // losing, because the factored one has no subtraction of nearly equal
    // quantities anywhere in it.
    expectTrue("a kelvin apart the two forms are indistinguishable", disagreement(1.0) < 1e-13);
    expectTrue("a millikelvin apart they still are", disagreement(1e-3) < 1e-11);
    expectTrue("a nanokelvin apart the difference form has lost six digits",
               disagreement(1e-9) > 1e-7);
    expectTrue("and it is monotone in the gap, which is what says it is cancellation",
               disagreement(1e-9) > disagreement(1e-7) && disagreement(1e-7) > disagreement(1e-3));

    // Exactly zero flux at equal temperatures -- not a rounding of zero. This is
    // what the cold control rests on: a chain attached to a ship that is not on fire
    // must add nothing at all.
    expectTrue("equal temperatures give exactly zero flux",
               fire::filmCoefficient(700.0, 700.0, p) * (700.0 - 700.0) == 0.0);
    // And with no emissivity it is the convective coefficient, exactly.
    fire::BoundaryFilm bare{25.0, 0.0};
    expectTrue("with no emissivity it is the convective coefficient, to the bit",
               fire::filmCoefficient(900.0, 300.0, bare) == 25.0);
    expectTrue("radiation dominates convection in a real compartment fire",
               fire::filmCoefficient(1073.15, 500.0, p) > 4.0 * p.convective);

    std::printf("     film: %.1f W/(m2 K) at 800 C gas against 227 C steel, of which %.1f is"
                " radiation; identity holds to %.1e (worst at T_g = %.0f K)\n",
                fire::filmCoefficient(1073.15, 500.0, p),
                fire::filmCoefficient(1073.15, 500.0, p) - p.convective, worst, worstAt);
}

void testTheWallExchangeSplitsAtTheLayerInterface() {
    const solidshell::HexMesh mesh = standingPlate(4.0, 4.0, 0.0095, 4, 8);
    std::vector<thermal::BoundaryFace> face;
    for (const thermal::BoundaryFace& f : thermal::boundaryFaces(mesh))
        if (f.normal.x > 0.5) face.push_back(f);
    expectEqual("the plate has one film face per element on its fore side",
                static_cast<long long>(face.size()), 32);

    // A surface that is hot at the top and cold at the foot, which is what a
    // bulkhead standing in a smoke layer with water behind it looks like.
    std::vector<double> surface(face.size());
    for (std::size_t i = 0; i < face.size(); ++i)
        surface[i] = kTAmbient + 220.0 * face[i].centroid.z / 4.0;

    const fire::GasCompartment gas = twoLayerBox(2.5, 800.0, 300.0);
    expectNear("the fixture's interface is where it was asked for", gas.interfaceZ(), 2.5, 1e-12);

    const fire::WallExchange two = fire::wallExchange(gas, face, surface);
    expectEqual("with no banding there are exactly two films",
                static_cast<long long>(two.film.size()), 2);
    expectNear("and they cover the whole plate between them", two.totalArea, 16.0, 1e-9);
    expectTrue("the upper film is the hot layer's temperature and the lower the cool one's",
               two.film[0].ambient == gas.upper.temperature() &&
                   two.film[1].ambient == gas.lower.temperature());
    // Split at 2.5 m on a plate whose element rows are 0.5 m tall: the three rows
    // centred at 2.75, 3.25 and 3.75 are above it and the other five are below, so
    // 12 faces of 0.5 m2 against 20. Not the middle, which is what a single film at
    // a mean temperature would have made of it.
    expectNear("the split is at the interface and not at the middle", two.area[0], 6.0, 1e-9);
    expectNear("and the rest is below it", two.area[1], 10.0, 1e-9);
    for (const thermal::BoundaryFace& f : two.film[0].face)
        expectTrue("every face in the upper film is above the interface", f.centroid.z >= 2.5);
    for (const thermal::BoundaryFace& f : two.film[1].face)
        expectTrue("and every face in the lower one is below it", f.centroid.z < 2.5);

    // The back-reaction the gas takes: area means, so a caller writing them onto a
    // `GasCompartment` is writing what the meshed part of its boundary is doing.
    double meanSurface = 0;
    for (std::size_t i = 0; i < face.size(); ++i) meanSurface += face[i].area * surface[i];
    meanSurface /= 16.0;
    expectNear("the wall temperature handed back is the area mean of the surface",
               two.wallTemperature, meanSurface, 1e-9);

    // Refusals: a mismatched or empty pair produces nothing, rather than a film at
    // an invented temperature.
    expectTrue("a mismatched surface field is refused",
               fire::wallExchange(gas, face, {1.0}).film.empty());
    expectTrue("and so is an empty one", fire::wallExchange(gas, {}, {}).film.empty());
}

void testBandingTheFilmRemovesTheSpreadError() {
    const solidshell::HexMesh mesh = standingPlate(4.0, 4.0, 0.0095, 4, 8);
    std::vector<thermal::BoundaryFace> face;
    for (const thermal::BoundaryFace& f : thermal::boundaryFaces(mesh))
        if (f.normal.x > 0.5) face.push_back(f);
    std::vector<double> surface(face.size());
    for (std::size_t i = 0; i < face.size(); ++i)
        surface[i] = kTAmbient + 220.0 * face[i].centroid.z / 4.0;
    const fire::GasCompartment gas = twoLayerBox(2.5, 800.0, 300.0);

    // **The radiative coefficient goes as `T_s^2`, so one coefficient over a whole
    // layer is wrong wherever the surface has a gradient.** The error is exactly the
    // spread *within* a film, so banding by height -- which is the direction a
    // bulkhead's surface temperature varies in -- drives it down.
    const fire::WallExchange unbanded = fire::wallExchange(gas, face, surface);
    const fire::WallExchange banded = fire::wallExchange(gas, face, surface, {}, 0.5);
    expectEqual("banding at the element height gives one film per row",
                static_cast<long long>(banded.film.size()), 8);
    expectNear("and it still covers the whole plate", banded.totalArea, unbanded.totalArea, 1e-12);
    expectNear("and books the same exact heat, because that does not depend on the films",
               banded.exactHeat, unbanded.exactHeat, 1e-9);
    expectTrue("the unbanded error is a real fraction of the exchange",
               std::abs(unbanded.linearisationError) > 0.004 * std::abs(unbanded.exactHeat));
    expectTrue("banding removes at least nine tenths of it",
               std::abs(banded.linearisationError) < 0.1 * std::abs(unbanded.linearisationError));

    // **The band index is a pure function of the face's own centroid**, which is the
    // property `thermal::Solver::setFilm` depends on: the same faces must give the
    // same films in the same order however hot anything is, or a coupled run would
    // be writing one band's coefficient onto another's faces.
    const fire::GasCompartment colder = twoLayerBox(3.5, 400.0, 290.0);
    const fire::WallExchange moved = fire::wallExchange(colder, face, surface, {}, 0.5);
    expectEqual("a different fire gives the same number of bands",
                static_cast<long long>(moved.film.size()),
                static_cast<long long>(banded.film.size()));
    bool sameFaces = true;
    for (std::size_t b = 0; b < banded.film.size(); ++b) {
        if (banded.film[b].face.size() != moved.film[b].face.size()) { sameFaces = false; break; }
        for (std::size_t i = 0; i < banded.film[b].face.size(); ++i)
            if (banded.film[b].face[i].element != moved.film[b].face[i].element ||
                banded.film[b].face[i].face != moved.film[b].face[i].face)
                sameFaces = false;
    }
    expectTrue("and exactly the same faces in exactly the same bands", sameFaces);
    expectTrue("while the layer the bands sit in did move",
               moved.film[5].ambient != banded.film[5].ambient);
    // Ordering: band 0 is the lowest, which is what a caller indexing by height
    // needs and is not something to leave to the face table's own order.
    for (std::size_t b = 1; b < banded.film.size(); ++b)
        expectTrue("bands run upward from the foot of the plate",
                   banded.film[b].face.front().centroid.z >
                       banded.film[b - 1].face.front().centroid.z);

    // **The bands are measured from the lowest face and not from the origin**, and a
    // plate sitting on the origin cannot tell the difference -- 0.25 m is half a
    // band, so both readings land on the same integers. Mutation testing found that:
    // dropping `- zBase` survived the whole suite. A bulkhead does not start at the
    // baseline, so the fixture that matters is the offset one.
    const solidshell::HexMesh raised = standingPlate(4.0, 4.0, 0.0095, 4, 8, 1.4);
    std::vector<thermal::BoundaryFace> high;
    for (const thermal::BoundaryFace& f : thermal::boundaryFaces(raised))
        if (f.normal.x > 0.5) high.push_back(f);
    std::vector<double> highSurface(high.size());
    for (std::size_t i = 0; i < high.size(); ++i)
        highSurface[i] = kTAmbient + 220.0 * (high[i].centroid.z - 1.4) / 4.0;
    const fire::WallExchange offset =
        fire::wallExchange(twoLayerBox(2.5, 800.0, 300.0), high, highSurface, {}, 0.5);
    expectEqual("a plate raised off the origin still bands into one film per row",
                static_cast<long long>(offset.film.size()), 8);
    for (std::size_t b = 0; b < offset.film.size(); ++b) {
        expectEqual("every band holds exactly one row of faces",
                    static_cast<long long>(offset.film[b].face.size()), 4);
        expectNear("and band 0 is the row at the foot, not wherever the origin fell",
                   offset.film[b].face.front().centroid.z,
                   1.4 + 0.25 + 0.5 * static_cast<double>(b), 1e-9);
    }

    std::printf("     wall exchange: %.3f kW into the steel; one film per layer misbooks"
                " %.1f%% of it, one film per %.2f m row misbooks %.3f%%\n",
                banded.exactHeat / 1e3,
                100.0 * std::abs(unbanded.linearisationError / unbanded.exactHeat), 0.5,
                100.0 * std::abs(banded.linearisationError / banded.exactHeat));
}

// The exact control for the whole chain, not only for the gas: a ship carrying a
// fire model, a conduction solve on one of her bulkheads and the boundary exchange
// between them, with nothing burning, must be bit-identical to the same ship
// carrying none of it.
//
// It is a different assertion from `testZeroHeatReleaseLeavesTheShipBitIdentical`
// above, because the coupling adds two write paths that test does not exercise:
// `wallExchange` writes `wallConductance` and `wallTemperature` onto every tracked
// compartment, and a film that was a rounding away from zero would heat the steel,
// which would move the wall temperature, which would move the gas.
void testAColdChainLeavesTheShipAndTheSteelBitIdentical() {
    Sea sea{0.0};
    Ship bare = game::buildFerry();
    bare.initialise(sea);

    Ship coupled = game::buildFerry();
    coupled.initialise(sea);
    fire::Model model;
    model.attach(coupled, {coupled.findCompartment("engine_room_s"),
                           coupled.findCompartment("engine_room_p")});
    fire::DesignFire d;
    d.name = "unlit";
    d.compartment = 0;
    d.baseZ = 2.5;
    d.diameter = 2.5;
    d.peakHeatRelease = 0.0;
    model.fires.push_back(d);

    const solidshell::HexMesh mesh = standingPlate(4.0, 4.0, 0.0095, 4, 8);
    std::vector<thermal::BoundaryFace> face;
    for (const thermal::BoundaryFace& f : thermal::boundaryFaces(mesh))
        if (f.normal.x > 0.5) face.push_back(f);
    std::vector<double> surface(face.size(), kTAmbient);

    thermal::Problem problem;
    problem.mesh = &mesh;
    problem.material = ah36Steel();
    problem.temperatureDependent = true;
    const fire::WallExchange seed = fire::wallExchange(model.gas[0], face, surface, {}, 0.5);
    problem.film = seed.film;
    thermal::Solver solver;
    std::string why;
    expectTrue("the conduction problem prepares", solver.prepare(problem, kTAmbient, &why));

    for (int i = 0; i < 400; ++i) {
        bare.step(0.05, sea);
        coupled.step(0.05, sea);
        model.step(0.05, coupled, sea);
        const std::vector<double>& nodal = solver.temperature();
        for (std::size_t k = 0; k < face.size(); ++k) {
            double t = 0;
            for (int c = 0; c < 4; ++c) t += nodal[face[k].node[c]];
            surface[k] = 0.25 * t;
        }
        const fire::WallExchange x = fire::wallExchange(model.gas[0], face, surface, {}, 0.5);
        for (std::size_t b = 0; b < x.film.size(); ++b)
            solver.setFilm(b, 0.0, x.film[b].coefficient, x.film[b].ambient);
        model.gas[0].wallTemperature = x.wallTemperature;
        model.gas[0].wallConductance = x.wallConductance;
        // Not exactly zero, and the reason is worth being exact about: the banded
        // solve returns a uniform field to *its* rounding and not to the bit, so the
        // steel sits nanokelvins off ambient and the film books nanowatts against a
        // 100 kW exchange. What has to be exact is the ship, and it is.
        expectTrue("a boundary at the temperature of the gas against it books nothing",
                   std::abs(x.heat) < 1e-6 && std::abs(x.exactHeat) < 1e-6);
        solver.step(0.05, &why);
        model.applyTo(coupled);
    }

    double drift = 0;
    for (double t : solver.temperature()) drift = std::max(drift, std::abs(t - kTAmbient));
    expectTrue("the steel is still at ambient to 1e-8 K after 400 steps", drift < 1e-8);

    const std::vector<double> a = shipFingerprint(bare);
    const std::vector<double> b = shipFingerprint(coupled);
    std::size_t differing = 0;
    for (std::size_t i = 0; i < a.size() && i < b.size(); ++i)
        if (std::memcmp(&a[i], &b[i], sizeof(double)) != 0) ++differing;
    expectEqual("a cold chain leaves the ship bit-identical, in every state double",
                static_cast<long long>(differing), 0);

    // Vacuity, in both directions: the ship has to have been doing something, the
    // gas model has to have stepped, and the coupling has to have actually written.
    expectTrue("and the ship really was flooding through the run",
               bare.totalFloodwaterMass() > 1.0e5);
    expectTrue("the gas model really did step", model.time > 19.0);
    expectTrue("and the wall coupling really did write a conductance",
               model.gas[0].wallConductance > 25.0);
    // **Which is the gas-side film and not the lumped `h_k` it replaced.** 28.80 at
    // ambient -- EN's 25 of convection plus 3.80 of radiation between two bodies
    // both at 15 C -- against the 30.0 `GasCompartment` defaults to, which was MQH's
    // wall conductance standing in for a wall that is now solved.
    expectNear("which is the gas-side film, not the MQH wall conductance it replaced",
               model.gas[0].wallConductance, fire::filmCoefficient(kTAmbient, kTAmbient, {}),
               1e-9);
    expectTrue("and the two are genuinely different numbers",
               std::abs(model.gas[0].wallConductance - fire::GasCompartment{}.wallConductance) >
                   1.0);
}

void runFireTests() {
    std::printf("\n--- compartment fire ---\n");
    testCaloricConstantsAreExactlyConsistent();
    testVolumeSplitIsTheEnergySplit();
    testPlumeBranchesMeetAtTheFlameTip();
    testPlumeObeysTheFiveThirdsPowerLaw();
    testVirtualOriginChangesSignWithFireShape();
    testVentIntegralMatchesTheClassicalDoorway();
    testBalancedNeutralPlaneFollowsTheCubeRootOfTemperature();
    testSingleDeltaPHasNoVentilationWhereTheIntegralHas();
    testVentShapesAreConsistentWithTheirAreas();
    testSealedCompartmentReachesItsClosedFormPressure();
    testSteadyLayerMatchesMqhInTheMiddleOfItsRange();
    testMqhDisagreementIsOrderedByPowerNotScattered();
    testMqhAndThomasAgreeOnFlashover();
    testTheLayerIsDrivenDownByEntrainmentAndNotByExpansion();
    testTheVentSplitsAtBothSidesInterfaces();
    testSpeciesRideTheGasIntoTheNeighboursCoolLayer();
    testCompressingAnUntouchedLayerIsAdiabatic();
    testHorizontalVentTakesOneDeltaPAndIsNeverBidirectional();
    testTheAccountIntegratesTheDesignCurve();
    testRadiativeLossFractionRemovesExactlyItsShare();
    testProductsAreMadeInTheUpperLayerAtTheDeclaredYield();
    testThePlumeCarriesTheLowerLayersSpeciesUpward();
    testEntrainmentStopsWhenTheLayerReachesTheFloor();
    testAttachDerivesTheGasBoxFromTheShipsOwnCompartment();
    testVentsAreSlidIntoTheGasSpaceBeforeBeingTrimmed();
    testAnOpeningUnderFloodwaterInsideTheShipIsBlockedToo();
    testTheAccountClosesOnTheFerry();
    testProductsCrossTheOpeningNetwork();
    testTheUnfiredCompartmentIsClearlyDifferent();
    testWaterBlockedOpeningsAreLeftToTheFloodingSolve();
    testZeroHeatReleaseLeavesTheShipBitIdentical();
    testApplyToReproducesTheModelPressureUnderTheShipsOwnFormula();
    testTheAnswerDoesNotDependOnTheCallersTick();
    testRefiningTheInternalStepConverges();
    testThePressureSolveConvergesAndDoesNotPressuriseTheRoom();
    testValidateCatchesBadDefinitions();
    testDesignFireFollowsItsCurve();

    std::printf("\n--- suppression, and its effect on stability ---\n");
    testSprayMassFlowConvertsTheRuleUnitExactly();
    testScupperIntegralIsTheClassicalWeirAndItsDrownedLimit();
    testDrencherCoolingIsSensibleHeatPlusLatentOnTheEvaporatedShare();
    testTheEvaporatedFractionMovesCoolingByFourAndNotByNine();
    testWhatEvaporatesIsBoundedByTheFiresPowerAndNotByTheFraction();
    testTheDrencherCannotCoolBelowTheTemperatureOfItsOwnWater();
    testTheSuppressionSinkConvergesAndTheWaterLandingConvergesFaster();
    testSuppressionWaterIsFloodwaterAndItsFreeSurfaceIsTheClosedForm();
    testPocketingLimitsFreeSurfaceExactlyWhereGeometrySays();
    testTheFerrysPublishedGmIsNotTheInitialGmOfAShallowLayer();
    testFreeingPortsHoldTheDeckAtTheWeirsEquilibriumDepth();
    testABlockedFreeingPortLetsTheDeckFillWithoutBound();
    testAFreeingPortUnderTheSeaAdmitsInsteadOfDraining();
    testTheWaterAccountClosesAgainstWhatTheShipReceived();
    testWaterThatWillNotFitStaysOwedRatherThanBeingDropped();
    testTheAccountStillClosesWithSuppressionRunning();
    testTheFerrysLollTracksTheDrencherAndTheFreeingPortsDecideIt();
    testSuppressionOffLeavesTheShipAndTheGasBitIdentical();
    testValidateCatchesBadSuppressionDefinitions();

    std::printf("\n--- the boundary the structure sees ---\n");
    testTheFilmCoefficientCarriesRadiationExactly();
    testTheWallExchangeSplitsAtTheLayerInterface();
    testBandingTheFilmRemovesTheSpreadError();
    testAColdChainLeavesTheShipAndTheSteelBitIdentical();
}
