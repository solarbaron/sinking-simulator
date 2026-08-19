// SPDX-License-Identifier: MIT
#include "les.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <cstddef>

namespace sim::les {
namespace {

// Below this a cell is treated as empty rather than as a very thin gas, on exactly
// the terms `fire.cpp`'s `kMassFloor` is: `U/(m c_v)` with `m` at machine epsilon
// is a number, not a temperature. The substep controller rejects any step that
// would approach it, so this is the second line of defence.
constexpr double kCellMassFloor = 1e-12;   // kg

constexpr double kWattsPerKilowatt = 1000.0;

double clamp01(double x) { return x < 0.0 ? 0.0 : (x > 1.0 ? 1.0 : x); }

// Compensated (Kahan) summation, and it is not a nicety here.
//
// A naive sum over a few thousand cells carries a relative error around
// `sqrt(n) eps`, and on this model that noise is **larger than the quantity it is
// being used to measure**: the round-trip interface came back to 8.6e-14 m on a
// 400-cell box with a naive sum and to a few times 1e-16 m with this, because the
// two-layer reduction divides `M/A - rho_u H` -- a difference of two numbers three
// times its own size -- by a density gap. The conservation residuals fall by the
// same factor. A test cannot be tighter than its own arithmetic, and every
// assertion this file is under is a conservation assertion.
class Accumulator {
public:
    void add(double v) {
        const double y = v - carry_;
        const double t = total_ + y;
        carry_ = (t - total_) - y;
        total_ = t;
    }
    double total() const { return total_; }

private:
    double total_ = 0, carry_ = 0;
};

}  // namespace

// ---------------------------------------------------------------------------
// Field accessors
// ---------------------------------------------------------------------------

double Field::pressure() const {
    const double v = grid.volume();
    if (!(v > 0)) return kPatm;
    return (kGammaAir - 1.0) * energy / v;
}

double Field::cellEnergy() const {
    const int n = grid.cells();
    return n > 0 ? energy / static_cast<double>(n) : 0.0;
}

double Field::totalMass() const {
    Accumulator total;
    for (double m : mass) total.add(m);
    return total.total();
}

double Field::totalProducts() const {
    Accumulator total;
    for (double s : products) total.add(s);
    return total.total();
}

double Field::density(int cell) const {
    const double v = grid.cellVolume();
    if (!(v > 0) || cell < 0 || static_cast<std::size_t>(cell) >= mass.size()) return 0.0;
    return mass[static_cast<std::size_t>(cell)] / v;
}

double Field::temperature(int cell) const {
    if (cell < 0 || static_cast<std::size_t>(cell) >= mass.size()) return kTAmbient;
    const double m = mass[static_cast<std::size_t>(cell)];
    const double u = cellEnergy();
    if (m <= kCellMassFloor || u <= 0) return kTAmbient;
    return u / (m * fire::kCvAir);
}

double Field::meanTemperature() const {
    const double m = totalMass();
    if (m <= kCellMassFloor || energy <= 0) return kTAmbient;
    return energy / (m * fire::kCvAir);
}

double Field::rowDensity(int k) const {
    if (k < 0 || k >= grid.n[2]) return 0.0;
    Accumulator total;
    for (int j = 0; j < grid.n[1]; ++j)
        for (int i = 0; i < grid.n[0]; ++i)
            total.add(mass[static_cast<std::size_t>(grid.cell(i, j, k))]);
    const double rowVolume =
        grid.cellVolume() * static_cast<double>(grid.n[0]) * static_cast<double>(grid.n[1]);
    return rowVolume > 0 ? total.total() / rowVolume : 0.0;
}

// ---------------------------------------------------------------------------
// Account
// ---------------------------------------------------------------------------

double Account::energyResidualFraction() const {
    const double scale = std::max({std::abs(heatReleased), std::abs(wallLoss),
                                   std::abs(energy - initialEnergy), std::abs(energy), 1.0});
    return energyResidual() / scale;
}

double Account::massResidualFraction() const {
    const double scale = std::max({std::abs(mass - initialMass), std::abs(mass), 1e-9});
    return massResidual() / scale;
}

// ---------------------------------------------------------------------------
// Geometry
// ---------------------------------------------------------------------------

bool planRectangle(double floorArea, double perimeter, double* lengthX, double* lengthY) {
    const double side = floorArea > 0 ? std::sqrt(floorArea) : 0.0;
    *lengthX = side;
    *lengthY = side;
    if (!(floorArea > 0)) return false;
    const double half = 0.5 * perimeter;
    // `L + W = P/2`, `L W = A`: the sides are the roots of `t^2 - (P/2) t + A = 0`.
    // A negative discriminant means no rectangle has both, which is the ordinary
    // case on a ship -- the perimeter is the bounding box's and the area is the
    // prismatic equivalent. The square keeps the *area*, because the area is what
    // the interface height and the wall exchange are computed from.
    const double discriminant = half * half - 4.0 * floorArea;
    if (!(half > 0) || !(discriminant >= 0)) return false;
    const double root = std::sqrt(discriminant);
    const double a = 0.5 * (half + root);
    const double b = 0.5 * (half - root);
    if (!(b > 0)) return false;
    *lengthX = a;
    *lengthY = b;
    return true;
}

Grid gridFor(const fire::GasCompartment& gas, const Params& params) {
    Grid g;
    const double height = gas.ceilingZ - gas.floorZ;
    if (!(height > 0) || !(gas.floorArea > 0)) return g;

    double lx = 0, ly = 0;
    g.squareFallback = !planRectangle(gas.floorArea, gas.perimeter, &lx, &ly);
    if (!(lx > 0) || !(ly > 0)) return g;

    const int minPer = std::max(1, params.minPerAxis);
    double size = params.cellSize > 0 ? params.cellSize : 0.5;
    int n[3] = {minPer, minPer, minPer};
    const double extent[3] = {lx, ly, height};
    for (int attempt = 0; attempt < 64; ++attempt) {
        long long total = 1;
        for (int a = 0; a < 3; ++a) {
            const long long want = std::lround(extent[a] / size);
            n[a] = static_cast<int>(std::clamp<long long>(want, minPer, 512));
            total *= n[a];
        }
        if (params.maxCells <= 0 || total <= params.maxCells) break;
        // Grow the cell towards the budget rather than stepping blindly: the count
        // goes as the cube of the size, so one cube root is nearly the answer and
        // the 1.02 is what stops rounding from leaving it one cell over.
        size *= 1.02 * std::cbrt(static_cast<double>(total) /
                                 static_cast<double>(std::max(params.maxCells, 1)));
    }

    for (int a = 0; a < 3; ++a) {
        g.n[a] = n[a];
        g.h[a] = extent[a] / n[a];
    }
    // Centred on the body-frame origin in plan: a `fire::GasCompartment` carries no
    // horizontal position, so there is none to honour. `HeatSource` is therefore
    // written as an offset from the plan centre rather than as a body-frame point.
    g.lo[0] = -0.5 * lx;
    g.lo[1] = -0.5 * ly;
    g.lo[2] = gas.floorZ;
    return g;
}

int estimateCells(const fire::GasCompartment& gas, const Params& params) {
    return gridFor(gas, params).cells();
}

// ---------------------------------------------------------------------------
// Across the boundary
// ---------------------------------------------------------------------------

namespace {

// The fraction of cell row `k` that lies above `z`.
double aboveFraction(const Grid& g, int k, double z) {
    const double lo = g.lo[2] + static_cast<double>(k) * g.h[2];
    return clamp01((lo + g.h[2] - z) / g.h[2]);
}

}  // namespace

Field promote(const fire::GasCompartment& gas, const Params& params) {
    Field f;
    f.grid = gridFor(gas, params);
    f.floorZ = gas.floorZ;
    f.ceilingZ = gas.ceilingZ;
    f.wallConductance = gas.wallConductance;
    f.wallTemperature = gas.wallTemperature;
    f.energy = gas.totalEnergy();
    f.interfaceZ = gas.interfaceZ();
    if (f.grid.empty()) {
        f.problems.push_back("the compartment has no box to resolve: no floor area or no height");
        return f;
    }
    double lx = 0, ly = 0;
    if (!planRectangle(gas.floorArea, gas.perimeter, &lx, &ly))
        f.problems.push_back("no rectangle has this compartment's area and perimeter together;"
                             " the grid is a square of the same area");

    const Grid& g = f.grid;
    const int cells = g.cells();
    f.mass.assign(static_cast<std::size_t>(cells), 0.0);
    f.products.assign(static_cast<std::size_t>(cells), 0.0);
    for (int a = 0; a < 3; ++a)
        f.velocity[a].assign(static_cast<std::size_t>(g.faces(a)), 0.0);

    const double cellVolume = g.cellVolume();
    // The upper region's volume as the **grid** measures it, not as
    // `floorArea * (ceiling - z_i)` measures it. Using the grid's own is what makes
    // the shares sum to one and therefore what makes the mapping exact.
    double upperVolume = 0;
    for (int k = 0; k < g.n[2]; ++k)
        upperVolume += aboveFraction(g, k, f.interfaceZ) * cellVolume *
                       static_cast<double>(g.n[0]) * static_cast<double>(g.n[1]);
    const double lowerVolume = g.volume() - upperVolume;

    // A layer with volume but no mass is ordinary (the seed upper layer); a layer
    // with mass and no volume is not, and merging it into the other one is the only
    // answer that keeps the mass -- an infinitely dense layer is not a state.
    const bool haveUpper = upperVolume > 0;
    const bool haveLower = lowerVolume > 0;
    const double upperMass = haveUpper ? gas.upper.mass : 0.0;
    const double lowerMass = haveLower ? gas.lower.mass : 0.0;
    const double spillMass = gas.upper.mass + gas.lower.mass - upperMass - lowerMass;
    const double upperProducts = haveUpper ? gas.upper.products : 0.0;
    const double lowerProducts = haveLower ? gas.lower.products : 0.0;
    const double spillProducts =
        gas.upper.products + gas.lower.products - upperProducts - lowerProducts;

    for (int k = 0; k < g.n[2]; ++k) {
        const double above = aboveFraction(g, k, f.interfaceZ);
        const double below = 1.0 - above;
        double wUpper = haveUpper ? above * cellVolume / upperVolume : 0.0;
        double wLower = haveLower ? below * cellVolume / lowerVolume : 0.0;
        // Whichever region exists takes the homeless mass, spread by volume.
        if (haveUpper && !haveLower) wUpper = cellVolume / upperVolume;
        if (haveLower && !haveUpper) wLower = cellVolume / lowerVolume;
        const double m = upperMass * wUpper + lowerMass * wLower +
                         spillMass * (haveUpper ? wUpper : wLower);
        const double s = upperProducts * wUpper + lowerProducts * wLower +
                         spillProducts * (haveUpper ? wUpper : wLower);
        for (int j = 0; j < g.n[1]; ++j)
            for (int i = 0; i < g.n[0]; ++i) {
                const std::size_t c = static_cast<std::size_t>(g.cell(i, j, k));
                f.mass[c] = m;
                f.products[c] = s;
            }
    }
    return f;
}

namespace {

// The two-layer profile equivalent to a resolved field: the one that carries the
// field's own ceiling and floor densities and the field's own total mass.
//
// **Reading it off the profile rather than integrating the cells is the whole
// point, and it was found by measuring.** The obvious demotion sums the cell masses
// above the interface, and it comes back **9% wrong** on a field it has just
// promoted -- because promotion has to *smear* the interface across the cell row it
// falls inside, and no weighted sum of a smeared row can pull the two layers back
// apart. The reconstruction can: the smear leaves the top and bottom rows untouched
// and leaves the total untouched, and those three numbers are exactly what the
// equivalent profile is made of. Mass comes back to 1e-16 rather than to 1e-1.
struct Layers {
    double interfaceZ = 0;
    double upperDensity = 0, lowerDensity = 0;
    double upperFraction = 0;     // kg of products per kg of gas, in the top row
    double upperVolume = 0;       // m^3
    bool   determined = false;    // false when the box is uniform and there is none
};

double rowMass(const Field& f, int k) {
    Accumulator total;
    for (int j = 0; j < f.grid.n[1]; ++j)
        for (int i = 0; i < f.grid.n[0]; ++i)
            total.add(f.mass[static_cast<std::size_t>(f.grid.cell(i, j, k))]);
    return total.total();
}

double rowProducts(const Field& f, int k) {
    Accumulator total;
    for (int j = 0; j < f.grid.n[1]; ++j)
        for (int i = 0; i < f.grid.n[0]; ++i)
            total.add(f.products[static_cast<std::size_t>(f.grid.cell(i, j, k))]);
    return total.total();
}

Layers reduceToLayers(const Field& f, double held) {
    Layers out;
    out.interfaceZ = held;
    const Grid& g = f.grid;
    const double area = g.planArea();
    if (g.empty() || !(area > 0)) return out;

    const double rowVolume = area * g.h[2];
    const double topMass = rowMass(f, g.n[2] - 1);
    out.upperDensity = topMass / rowVolume;
    out.lowerDensity = rowMass(f, 0) / rowVolume;
    out.upperFraction = topMass > kCellMassFloor ? rowProducts(f, g.n[2] - 1) / topMass : 0.0;

    const double zFloor = g.lo[2];
    const double zCeiling = zFloor + static_cast<double>(g.n[2]) * g.h[2];
    // `rho_u (H - z_i) + rho_l (z_i - z_f) = M / A`, one linear equation in `z_i`.
    // A uniform box makes the coefficient zero and determines nothing at all, which
    // is the physics rather than a numerical difficulty -- every split of a uniform
    // gas is the same gas -- so the last determined interface is held. The floor is
    // relative, because "uniform" for a 1.2 kg/m^3 gas is not the same number as
    // for a 0.3 kg/m^3 one.
    const double gap = out.lowerDensity - out.upperDensity;
    const double scale = std::max(std::abs(out.lowerDensity), std::abs(out.upperDensity));
    if (gap > 1e-12 * std::max(scale, 1e-12)) {
        const double z =
            (f.totalMass() / area - out.upperDensity * zCeiling + out.lowerDensity * zFloor) / gap;
        out.interfaceZ = std::clamp(z, zFloor, zCeiling);
        out.determined = true;
    }
    out.upperVolume = area * (zCeiling - std::clamp(out.interfaceZ, zFloor, zCeiling));
    return out;
}

}  // namespace

double equivalentInterface(const Field& field, double held) {
    return reduceToLayers(field, held).interfaceZ;
}

void demote(const Field& field, fire::GasCompartment& gas) {
    if (field.empty()) return;
    const Layers layers = reduceToLayers(field, field.interfaceZ);
    const double volume = field.grid.volume();
    const double totalMass = field.totalMass();
    const double totalProducts = field.totalProducts();

    // The upper layer is *reconstructed* and the lower is **the remainder**, so the
    // two add to the total in floating point and not merely in algebra. Same trick,
    // same reason, as `fire.cpp`'s `layerSplit`.
    const double upperMass =
        std::clamp(layers.upperDensity * layers.upperVolume, 0.0, std::max(totalMass, 0.0));
    gas.upper.mass = upperMass;
    gas.lower.mass = totalMass - upperMass;

    const double upperProducts =
        std::clamp(layers.upperFraction * upperMass, 0.0, std::max(totalProducts, 0.0));
    gas.upper.products = upperProducts;
    gas.lower.products = totalProducts - upperProducts;

    // The energy split *is* the volume split, because the internal energy density is
    // uniform at a uniform pressure -- which is the identity `fire.hpp` states as
    // `V_u / V = U_u / U`. So this is not an approximation of the two-zone closure,
    // it is the closure.
    gas.upper.energy = volume > 0 ? field.energy * (layers.upperVolume / volume) : 0.0;
    gas.lower.energy = field.energy - gas.upper.energy;
}

std::vector<HeatSource> sourcesFor(const fire::Model& model, int gasCompartment, double atTime) {
    std::vector<HeatSource> out;
    for (const fire::DesignFire& fire : model.fires) {
        if (fire.compartment != gasCompartment) continue;
        const double q = fire.heatRelease(atTime);
        if (!(q > 0)) continue;
        HeatSource s;
        s.baseZ = fire.baseZ;
        s.diameter = fire.diameter;
        // Over the flame, not over the pan. `fire::Plume::flameHeight` is negative
        // for a wide, weak fire -- there is then no coherent flame -- and the clamp
        // to zero leaves the heat in the row the pan is in, which is the right
        // reading of that case.
        const fire::Plume plume{q, fire.diameter, fire.convectiveFraction};
        s.flameHeight = std::max(plume.flameHeight(), 0.0);
        s.power = q * (1.0 - fire.radiativeLossFraction);
        s.productRate = fire.productYield * q;
        out.push_back(s);
    }
    return out;
}

void resetAccount(const Field& field, Account& account) {
    account = Account{};
    account.initialEnergy = field.energy;
    account.initialMass = field.totalMass();
    account.energy = account.initialEnergy;
    account.mass = account.initialMass;
    account.products = field.totalProducts();
}

std::vector<std::string> validate(const Field& field) {
    std::vector<std::string> problems = field.problems;
    if (field.grid.empty()) {
        problems.push_back("the field has no cells");
        return problems;
    }
    if (!(field.ceilingZ > field.floorZ))
        problems.push_back("the field has a ceiling at or below its floor");
    if (!(field.energy > 0)) problems.push_back("the field carries no internal energy");
    for (std::size_t c = 0; c < field.mass.size(); ++c)
        if (!(field.mass[c] > kCellMassFloor)) {
            problems.push_back("a cell holds no gas, so its temperature is not defined");
            break;
        }
    for (int a = 0; a < 3; ++a)
        if (static_cast<int>(field.velocity[a].size()) != field.grid.faces(a))
            problems.push_back("the velocity component does not match the grid's face count");
    return problems;
}

// ---------------------------------------------------------------------------
// The correlations the criterion reads
// ---------------------------------------------------------------------------

double alpertCeilingJetRise(double heatRelease, double radius, double height) {
    if (!(heatRelease > 0) || !(height > 0)) return 0.0;
    const double qkw = heatRelease / kWattsPerKilowatt;
    const double r = std::max(radius, 0.0);
    if (r <= kAlpertCrossover * height)
        return kAlpertNear * std::cbrt(qkw * qkw) / (height * std::cbrt(height * height));
    return kAlpertFar * std::cbrt(qkw * qkw / (r * r)) / height;
}

double ceilingJetSpread(double radius, double height) {
    if (!(height > 0) || !(radius > 0)) return 1.0;
    const double ratio = radius / height;
    if (ratio <= kAlpertCrossover) return 1.0;
    return (kAlpertNear / kAlpertFar) * std::cbrt(ratio * ratio);
}

// ---------------------------------------------------------------------------
// The solve
// ---------------------------------------------------------------------------

namespace {

// Scratch, allocated once per `step()` and reused by every substep. Allocating per
// substep would dominate: a run that rejects its way down to millisecond steps
// takes tens of thousands of them.
struct Work {
    std::vector<double> heat, source, rhs, potential, accepted, eddy;
    std::vector<double> cellVelocity[3];
    std::vector<double> trial[3], flux[3];
    std::vector<double> deltaMass, deltaProducts, productSource;

    void size(const Grid& g) {
        const std::size_t n = static_cast<std::size_t>(g.cells());
        heat.assign(n, 0.0);
        source.assign(n, 0.0);
        rhs.assign(n, 0.0);
        potential.assign(n, 0.0);
        accepted.assign(n, 0.0);
        eddy.assign(n, 0.0);
        deltaMass.assign(n, 0.0);
        deltaProducts.assign(n, 0.0);
        productSource.assign(n, 0.0);
        for (int a = 0; a < 3; ++a) {
            cellVelocity[a].assign(n, 0.0);
            const std::size_t faces = static_cast<std::size_t>(g.faces(a));
            trial[a].assign(faces, 0.0);
            flux[a].assign(faces, 0.0);
        }
    }
};

// The exterior area of one cell: the box faces it owns. Summed over every cell this
// is `2(LxLy + LyLz + LxLz)`, which is the same wetted enclosure the two-zone model
// charges as `2 A + P H` -- so the two models pay the same boundary, cell by cell
// against layer by layer, and a comparison between them is not a comparison of two
// different rooms.
double exteriorArea(const Grid& g, int i, int j, int k) {
    double area = 0;
    if (i == 0) area += g.faceArea(0);
    if (i == g.n[0] - 1) area += g.faceArea(0);
    if (j == 0) area += g.faceArea(1);
    if (j == g.n[1] - 1) area += g.faceArea(1);
    if (k == 0) area += g.faceArea(2);
    if (k == g.n[2] - 1) area += g.faceArea(2);
    return area;
}

bool substep(Field& f, double dt, const std::vector<HeatSource>& sources, const Params& params,
             Work& w, Account& account, StepResult& out) {
    const Grid& g = f.grid;
    const int cells = g.cells();
    const double cellVolume = g.cellVolume();
    const double volume = g.volume();
    const double p0 = f.pressure();
    if (!(p0 > 0) || !(cellVolume > 0)) return true;

    std::fill(w.heat.begin(), w.heat.end(), 0.0);
    std::fill(w.productSource.begin(), w.productSource.end(), 0.0);

    // --- The fires ----------------------------------------------------------
    //
    // Spread over the **flame**, not over the pan. A megawatt released into one
    // 0.125 m^3 cell holding 0.15 kg is 6600 K/s and the substep controller would
    // chase it down to microseconds; released over the flame volume Heskestad's own
    // height gives, it is a hundredth of that and the step is set by the flow.
    double released = 0, generated = 0;
    for (const HeatSource& s : sources) {
        if (!(s.power > 0) && !(s.productRate > 0)) continue;
        const double radius = 0.5 * std::max(s.diameter, 0.0);
        const double top = s.baseZ + std::max(s.flameHeight, 0.0);
        int count = 0;
        const auto visit = [&](auto&& fn) {
            for (int k = 0; k < g.n[2]; ++k) {
                const double zLo = g.lo[2] + static_cast<double>(k) * g.h[2];
                if (zLo + g.h[2] <= s.baseZ || zLo > top) continue;
                for (int j = 0; j < g.n[1]; ++j)
                    for (int i = 0; i < g.n[0]; ++i) {
                        const double dx = g.centre(0, i) - s.offsetX;
                        const double dy = g.centre(1, j) - s.offsetY;
                        if (dx * dx + dy * dy > radius * radius) continue;
                        fn(g.cell(i, j, k));
                    }
            }
        };
        visit([&](int) { ++count; });
        if (count == 0) {
            // The fire falls between cell centres, or outside the box: the single
            // nearest cell takes it. Dropping it instead would be a fire that
            // releases nothing, which is indistinguishable from no fire at all.
            int best = -1;
            double nearest = 0;
            for (int k = 0; k < g.n[2]; ++k)
                for (int j = 0; j < g.n[1]; ++j)
                    for (int i = 0; i < g.n[0]; ++i) {
                        const double dx = g.centre(0, i) - s.offsetX;
                        const double dy = g.centre(1, j) - s.offsetY;
                        const double dz = g.centre(2, k) - 0.5 * (s.baseZ + top);
                        const double d = dx * dx + dy * dy + dz * dz;
                        if (best < 0 || d < nearest) {
                            nearest = d;
                            best = g.cell(i, j, k);
                        }
                    }
            if (best < 0) continue;
            w.heat[static_cast<std::size_t>(best)] += s.power;
            w.productSource[static_cast<std::size_t>(best)] += s.productRate;
        } else {
            const double share = 1.0 / static_cast<double>(count);
            visit([&](int c) {
                w.heat[static_cast<std::size_t>(c)] += s.power * share;
                w.productSource[static_cast<std::size_t>(c)] += s.productRate * share;
            });
        }
        released += s.power;
        generated += s.productRate;
    }

    // --- The boundary -------------------------------------------------------
    //
    // Exactly, not explicitly, and for the reason `fire.cpp` gives at the same term:
    // in isolation this is a linear relaxation towards the wall temperature whose
    // solution over a step is an exponential, so taking the average rate from that
    // solution costs one `expm1` and is unconditionally stable. The capacity is
    // `m c_p` and not `m c_v` because a **cell** is at constant pressure -- the
    // thermodynamic pressure is a property of the whole box -- while the box as a
    // whole is at constant volume and loses `q` from its internal energy. Both are
    // the first law; they are asking about different systems.
    double wallLoss = 0;
    if (f.wallConductance > 0) {
        for (int k = 0; k < g.n[2]; ++k)
            for (int j = 0; j < g.n[1]; ++j)
                for (int i = 0; i < g.n[0]; ++i) {
                    const double area = exteriorArea(g, i, j, k);
                    if (!(area > 0)) continue;
                    const int c = g.cell(i, j, k);
                    const double capacity = f.mass[static_cast<std::size_t>(c)] * fire::kCpAir;
                    const double rate = f.wallConductance * area;
                    if (!(capacity > 0) || !(rate > 0)) continue;
                    const double removed = capacity * (f.temperature(c) - f.wallTemperature) *
                                           -std::expm1(-rate * dt / capacity);
                    const double q = removed / dt;
                    w.heat[static_cast<std::size_t>(c)] -= q;
                    wallLoss += q;
                }
    }

    // The net, taken as the **sum of the per-cell sources** rather than as
    // `released - wallLoss`. The divergence source below is the deviation of the
    // heating from its own mean, and it is only a compatible right-hand side for an
    // all-Neumann projection if the mean is taken over the same numbers.
    Accumulator sum;
    for (int c = 0; c < cells; ++c) sum.add(w.heat[static_cast<std::size_t>(c)]);
    const double net = sum.total();

    // --- The divergence the heating demands ---------------------------------
    //
    //   div u = (gamma-1) q''' / (gamma p)  -  dp/dt / (gamma p),
    //   dp/dt = (gamma-1) Q / V   (the box is sealed, so the net divergence is zero)
    //
    // which leaves the deviation of the local heating from the box mean. That
    // `dp/dt` is `fire.hpp`'s own sealed closed form, recovered rather than imposed.
    const double expansion = (kGammaAir - 1.0) / (kGammaAir * p0);
    for (int c = 0; c < cells; ++c)
        w.source[static_cast<std::size_t>(c)] =
            expansion * (w.heat[static_cast<std::size_t>(c)] / cellVolume - net / volume);

    // --- Cell-centred velocities and the subgrid viscosity ------------------
    for (int a = 0; a < 3; ++a) {
        for (int k = 0; k < g.n[2]; ++k)
            for (int j = 0; j < g.n[1]; ++j)
                for (int i = 0; i < g.n[0]; ++i) {
                    int lo[3] = {i, j, k}, hi[3] = {i, j, k};
                    hi[a] += 1;
                    const double a0 = f.velocity[a][static_cast<std::size_t>(
                        g.face(a, lo[0], lo[1], lo[2]))];
                    const double a1 = f.velocity[a][static_cast<std::size_t>(
                        g.face(a, hi[0], hi[1], hi[2]))];
                    w.cellVelocity[a][static_cast<std::size_t>(g.cell(i, j, k))] = 0.5 * (a0 + a1);
                }
    }

    const double delta = std::cbrt(cellVolume);
    const double smagorinsky = params.smagorinsky * delta;
    const double smagorinsky2 = smagorinsky * smagorinsky;
    for (int k = 0; k < g.n[2]; ++k)
        for (int j = 0; j < g.n[1]; ++j)
            for (int i = 0; i < g.n[0]; ++i) {
                const int c = g.cell(i, j, k);
                const int index[3] = {i, j, k};
                double gradient[3][3];
                for (int a = 0; a < 3; ++a)
                    for (int b = 0; b < 3; ++b) {
                        int lo[3] = {i, j, k}, hi[3] = {i, j, k};
                        lo[b] = std::max(index[b] - 1, 0);
                        hi[b] = std::min(index[b] + 1, g.n[b] - 1);
                        const int span = hi[b] - lo[b];
                        if (span <= 0) {
                            gradient[a][b] = 0.0;
                            continue;
                        }
                        const double vLo =
                            w.cellVelocity[a][static_cast<std::size_t>(g.cell(lo[0], lo[1], lo[2]))];
                        const double vHi =
                            w.cellVelocity[a][static_cast<std::size_t>(g.cell(hi[0], hi[1], hi[2]))];
                        gradient[a][b] = (vHi - vLo) / (static_cast<double>(span) * g.h[b]);
                    }
                double norm = 0;
                for (int a = 0; a < 3; ++a)
                    for (int b = 0; b < 3; ++b) {
                        const double s = 0.5 * (gradient[a][b] + gradient[b][a]);
                        norm += 2.0 * s * s;
                    }
                w.eddy[static_cast<std::size_t>(c)] = smagorinsky2 * std::sqrt(norm);
            }

    // --- The buoyancy reference ---------------------------------------------
    //
    // Arbitrary up to a constant: a spatially uniform body force in a closed box is
    // exactly balanced by the pressure, so only the *horizontal* gradient of the
    // density survives the projection. The box's own mid-density is chosen because
    // it is exactly the cell density when the box is uniform, which makes the
    // buoyant acceleration exactly 0.0 on every face of a quiescent compartment
    // rather than a rounding of it.
    double reference = params.buoyancyReference;
    if (!(reference > 0)) {
        double lowest = f.density(0), highest = f.density(0);
        for (int c = 1; c < cells; ++c) {
            const double rho = f.density(c);
            lowest = std::min(lowest, rho);
            highest = std::max(highest, rho);
        }
        reference = 0.5 * (lowest + highest);
    }
    if (!(reference > 0)) return true;

    // --- Momentum -----------------------------------------------------------
    for (int a = 0; a < 3; ++a) {
        const int count[3] = {g.faceCount(a, 0), g.faceCount(a, 1), g.faceCount(a, 2)};
        for (int k = 0; k < count[2]; ++k)
            for (int j = 0; j < count[1]; ++j)
                for (int i = 0; i < count[0]; ++i) {
                    const int index[3] = {i, j, k};
                    const std::size_t face = static_cast<std::size_t>(g.face(a, i, j, k));
                    // The box is sealed: the normal velocity on every wall is zero,
                    // and stays zero through the projection as well.
                    if (index[a] == 0 || index[a] == count[a] - 1) {
                        w.trial[a][face] = 0.0;
                        continue;
                    }
                    int left[3] = {i, j, k}, right[3] = {i, j, k};
                    left[a] -= 1;
                    const int cl = g.cell(left[0], left[1], left[2]);
                    const int cr = g.cell(right[0], right[1], right[2]);
                    const double phi = f.velocity[a][face];

                    double advecting[3];
                    for (int b = 0; b < 3; ++b)
                        advecting[b] = b == a ? phi
                                              : 0.5 * (w.cellVelocity[b][static_cast<std::size_t>(cl)] +
                                                       w.cellVelocity[b][static_cast<std::size_t>(cr)]);

                    // Upwind advection and a Laplacian, accumulated **axis pair by
                    // axis pair**. Floating-point addition is commutative and not
                    // associative, so pairing the axes is what makes a mirror of the
                    // grid map the sum onto itself term for term -- and therefore
                    // what lets the symmetry test be an exact statement rather than
                    // a nearly-exact one.
                    double advection = 0, laplacian = 0;
                    for (int b = 0; b < 3; ++b) {
                        int lo[3] = {i, j, k}, hi[3] = {i, j, k};
                        lo[b] = std::max(index[b] - 1, 0);
                        hi[b] = std::min(index[b] + 1, count[b] - 1);
                        const double vLo =
                            f.velocity[a][static_cast<std::size_t>(g.face(a, lo[0], lo[1], lo[2]))];
                        const double vHi =
                            f.velocity[a][static_cast<std::size_t>(g.face(a, hi[0], hi[1], hi[2]))];
                        advection += advecting[b] *
                                     (advecting[b] > 0 ? (phi - vLo) : (vHi - phi)) / g.h[b];
                        laplacian += ((vLo - phi) + (vHi - phi)) / (g.h[b] * g.h[b]);
                    }

                    const double viscosity =
                        0.5 * (w.eddy[static_cast<std::size_t>(cl)] +
                               w.eddy[static_cast<std::size_t>(cr)]) + params.viscosity;
                    double buoyancy = 0;
                    if (a == 2) {
                        const double rho = 0.5 * (f.density(cl) + f.density(cr));
                        buoyancy = kGravity * (reference - rho) / reference;
                    }
                    w.trial[a][face] =
                        phi + dt * ((viscosity * laplacian - advection) + buoyancy);
                }
    }

    // --- The projection -----------------------------------------------------
    const double inverse[3] = {1.0 / (g.h[0] * g.h[0]), 1.0 / (g.h[1] * g.h[1]),
                               1.0 / (g.h[2] * g.h[2])};
    double rhsScale = 0;
    // Set by any non-finite residual below; one is already fatal to the answer.
    bool nonFinite = false;
    for (int k = 0; k < g.n[2]; ++k)
        for (int j = 0; j < g.n[1]; ++j)
            for (int i = 0; i < g.n[0]; ++i) {
                const int c = g.cell(i, j, k);
                const double dx = (w.trial[0][static_cast<std::size_t>(g.face(0, i + 1, j, k))] -
                                   w.trial[0][static_cast<std::size_t>(g.face(0, i, j, k))]) / g.h[0];
                const double dy = (w.trial[1][static_cast<std::size_t>(g.face(1, i, j + 1, k))] -
                                   w.trial[1][static_cast<std::size_t>(g.face(1, i, j, k))]) / g.h[1];
                const double dz = (w.trial[2][static_cast<std::size_t>(g.face(2, i, j, k + 1))] -
                                   w.trial[2][static_cast<std::size_t>(g.face(2, i, j, k))]) / g.h[2];
                const double value =
                    ((dx + dy) + dz) - w.source[static_cast<std::size_t>(c)];
                w.rhs[static_cast<std::size_t>(c)] = value;
                rhsScale = std::max(rhsScale, std::abs(value));
            }

    w.potential = w.accepted;
    const int strideY = g.n[0];
    const int strideZ = g.n[0] * g.n[1];
    const auto neighbourSum = [&](int i, int j, int k, int c, double* diagonal) {
        const int index[3] = {i, j, k};
        const int stride[3] = {1, strideY, strideZ};
        double total = 0;
        *diagonal = 0;
        for (int b = 0; b < 3; ++b) {
            double pair = 0;
            if (index[b] > 0) {
                pair += w.potential[static_cast<std::size_t>(c - stride[b])];
                *diagonal += inverse[b];
            }
            if (index[b] + 1 < g.n[b]) {
                pair += w.potential[static_cast<std::size_t>(c + stride[b])];
                *diagonal += inverse[b];
            }
            total += pair * inverse[b];
        }
        return total;
    };

    double residual = 0;
    int sweep = 0;
    // Relative to what the right-hand side actually is, with an absolute floor. The
    // floor is not tidiness: a divergence of 1e-14 per second is a velocity of
    // 1e-14 m/s across a metre, and without it a compartment doing very nearly
    // nothing would spend every one of its sweeps chasing round-off.
    const double tolerance = std::max(params.projectionTolerance * rhsScale, 1e-14);
    for (; sweep < std::max(params.projectionSweeps, 1); ++sweep) {
        // Red-black: within a colour every cell reads only the other colour, so a
        // half-sweep is a pure function of that colour's values and the answer does
        // not depend on the order the cells are visited in. That is what keeps the
        // solve deterministic under any traversal and what makes it commute with a
        // reflection of the grid.
        for (int colour = 0; colour < 2; ++colour)
            for (int k = 0; k < g.n[2]; ++k)
                for (int j = 0; j < g.n[1]; ++j)
                    for (int i = 0; i < g.n[0]; ++i) {
                        if (((i + j + k) & 1) != colour) continue;
                        const int c = g.cell(i, j, k);
                        double diagonal = 0;
                        const double total = neighbourSum(i, j, k, c, &diagonal);
                        if (!(diagonal > 0)) continue;
                        const double target =
                            (total - w.rhs[static_cast<std::size_t>(c)]) / diagonal;
                        w.potential[static_cast<std::size_t>(c)] +=
                            params.relaxation * (target - w.potential[static_cast<std::size_t>(c)]);
                    }
        residual = 0;
        for (int k = 0; k < g.n[2]; ++k)
            for (int j = 0; j < g.n[1]; ++j)
                for (int i = 0; i < g.n[0]; ++i) {
                    const int c = g.cell(i, j, k);
                    double diagonal = 0;
                    const double total = neighbourSum(i, j, k, c, &diagonal);
                    const double value = total - diagonal * w.potential[static_cast<std::size_t>(c)] -
                                         w.rhs[static_cast<std::size_t>(c)];
                    // **A non-finite residual must not vanish into the maximum.**
                    // `std::max(residual, x)` is `x < residual ? x : residual`, so
                    // `NaN < residual` is false and the NaN is dropped -- and every
                    // consequence follows: `residual` comes back 0, `rhsScale`
                    // above came back 0 the same way so `tolerance` is its 1e-14
                    // floor, `residual <= tolerance` is true, and the relaxation
                    // **breaks out after a single sweep** reporting `sweeps = 1`
                    // and `capped = false` over a field that is not a number.
                    //
                    // Same defect and same repair as `flip.cpp`'s pressure
                    // projection, which is the other Poisson solve in this tree
                    // and was written from this one. Only this fold is guarded:
                    // `w.rhs` feeds `value` directly, so a non-finite right-hand
                    // side reaches here too, and a second test beside `rhsScale`
                    // would be one no test could tell from its own absence --
                    // which is what a mutation sweep scored it there.
                    if (std::isfinite(value))
                        residual = std::max(residual, std::abs(value));
                    else
                        nonFinite = true;
                }
        // Breaking here on a non-finite field is right and not a shortcut: the
        // residual cannot fall, so the remaining sweeps are work that changes
        // nothing. The report below is what carries the failure.
        if (residual <= tolerance) break;
    }
    const int sweeps = sweep + 1;
    // A relaxation handed a NaN has not converged, and said it had. `capped` rather
    // than a new flag because every consumer already asserts on it --
    // `test_promotion.cpp` at two sites -- and "it ran out of sweeps" and "it never
    // started" are both *the answer is not trustworthy*. The residual is reported
    // as NaN so that a threshold test on it fails too.
    const bool capped = nonFinite || residual > tolerance;
    if (nonFinite) residual = std::numeric_limits<double>::quiet_NaN();

    for (int a = 0; a < 3; ++a) {
        const int count[3] = {g.faceCount(a, 0), g.faceCount(a, 1), g.faceCount(a, 2)};
        for (int k = 0; k < count[2]; ++k)
            for (int j = 0; j < count[1]; ++j)
                for (int i = 0; i < count[0]; ++i) {
                    const int index[3] = {i, j, k};
                    if (index[a] == 0 || index[a] == count[a] - 1) continue;
                    int left[3] = {i, j, k}, right[3] = {i, j, k};
                    left[a] -= 1;
                    const std::size_t face = static_cast<std::size_t>(g.face(a, i, j, k));
                    const double lo =
                        w.potential[static_cast<std::size_t>(g.cell(left[0], left[1], left[2]))];
                    const double hi =
                        w.potential[static_cast<std::size_t>(g.cell(right[0], right[1], right[2]))];
                    w.trial[a][face] -= (hi - lo) / g.h[a];
                }
    }

    // --- Transport ----------------------------------------------------------
    double courant = 0, peak = 0;
    for (int a = 0; a < 3; ++a) {
        const int count[3] = {g.faceCount(a, 0), g.faceCount(a, 1), g.faceCount(a, 2)};
        for (int k = 0; k < count[2]; ++k)
            for (int j = 0; j < count[1]; ++j)
                for (int i = 0; i < count[0]; ++i) {
                    const int index[3] = {i, j, k};
                    const std::size_t face = static_cast<std::size_t>(g.face(a, i, j, k));
                    if (index[a] == 0 || index[a] == count[a] - 1) {
                        w.flux[a][face] = 0.0;
                        continue;
                    }
                    int left[3] = {i, j, k}, right[3] = {i, j, k};
                    left[a] -= 1;
                    const double speed = w.trial[a][face];
                    peak = std::max(peak, std::abs(speed));
                    courant = std::max(courant, std::abs(speed) * dt / g.h[a]);
                    const int donor = speed > 0 ? g.cell(left[0], left[1], left[2])
                                                : g.cell(right[0], right[1], right[2]);
                    w.flux[a][face] =
                        speed * g.faceArea(a) * f.mass[static_cast<std::size_t>(donor)] / cellVolume;
                }
    }

    for (int k = 0; k < g.n[2]; ++k)
        for (int j = 0; j < g.n[1]; ++j)
            for (int i = 0; i < g.n[0]; ++i) {
                const std::size_t c = static_cast<std::size_t>(g.cell(i, j, k));
                // Gathered per cell from the two faces of each axis, rather than
                // scattered per face into two cells. The gather is what makes the
                // update a *pure function of the cell*, so a reflected grid produces
                // a reflected answer to the last bit; the scatter would depend on
                // the order the faces happened to be visited in.
                const std::size_t xm = static_cast<std::size_t>(g.face(0, i, j, k));
                const std::size_t xp = static_cast<std::size_t>(g.face(0, i + 1, j, k));
                const std::size_t ym = static_cast<std::size_t>(g.face(1, i, j, k));
                const std::size_t yp = static_cast<std::size_t>(g.face(1, i, j + 1, k));
                const std::size_t zm = static_cast<std::size_t>(g.face(2, i, j, k));
                const std::size_t zp = static_cast<std::size_t>(g.face(2, i, j, k + 1));
                const double dx = w.flux[0][xm] - w.flux[0][xp];
                const double dy = w.flux[1][ym] - w.flux[1][yp];
                const double dz = w.flux[2][zm] - w.flux[2][zp];
                w.deltaMass[c] = (dx + dy) + dz;

                // The donor of a face is the cell **behind** it, which is not the
                // same cell for the two faces of an axis: the low face's donor at
                // positive flux is the neighbour, the high face's is this cell. The
                // clamped indices are only ever reached on a wall face, where the
                // flux is exactly zero and the lambda has already returned.
                const int here = g.cell(i, j, k);
                const auto carried = [&](std::size_t face, int a, int minus, int plus) {
                    const double flux = w.flux[a][face];
                    if (flux == 0.0) return 0.0;
                    const std::size_t d = static_cast<std::size_t>(flux > 0 ? minus : plus);
                    return f.mass[d] > kCellMassFloor ? flux * (f.products[d] / f.mass[d]) : 0.0;
                };
                const double sx = carried(xm, 0, g.cell(std::max(i - 1, 0), j, k), here) -
                                  carried(xp, 0, here, g.cell(std::min(i + 1, g.n[0] - 1), j, k));
                const double sy = carried(ym, 1, g.cell(i, std::max(j - 1, 0), k), here) -
                                  carried(yp, 1, here, g.cell(i, std::min(j + 1, g.n[1] - 1), k));
                const double sz = carried(zm, 2, g.cell(i, j, std::max(k - 1, 0)), here) -
                                  carried(zp, 2, here, g.cell(i, j, std::min(k + 1, g.n[2] - 1)));
                w.deltaProducts[c] = ((sx + sy) + sz) + w.productSource[c];
            }

    // --- Accept or reject ---------------------------------------------------
    //
    // Nothing above has touched the field or the account, so a rejection is free and
    // leaves no trace -- including the projection's warm start, which is restored
    // from the last *accepted* substep at the top of every attempt.
    if (dt > 1e-9) {
        if (courant > params.courant) return false;
        for (int c = 0; c < cells; ++c) {
            const std::size_t index = static_cast<std::size_t>(c);
            const double m = f.mass[index];
            if (std::abs(w.deltaMass[index]) * dt > params.maxRelativeChange * m) return false;
            if (m + w.deltaMass[index] * dt <= kCellMassFloor) return false;
            if (f.products[index] + w.deltaProducts[index] * dt < 0) return false;
        }
        if (f.energy + net * dt <= 0) return false;
    }

    // --- Commit -------------------------------------------------------------
    for (int c = 0; c < cells; ++c) {
        const std::size_t index = static_cast<std::size_t>(c);
        f.mass[index] += w.deltaMass[index] * dt;
        f.products[index] += w.deltaProducts[index] * dt;
    }
    for (int a = 0; a < 3; ++a) f.velocity[a] = w.trial[a];
    w.accepted = w.potential;
    f.energy += net * dt;
    f.time += dt;
    f.interfaceZ = equivalentInterface(f, f.interfaceZ);

    account.heatReleased += released * dt;
    account.wallLoss += wallLoss * dt;
    account.productsGenerated += generated * dt;
    account.energy = f.energy;
    account.mass = f.totalMass();
    account.products = f.totalProducts();

    out.courant = std::max(out.courant, courant);
    out.peakSpeed = std::max(out.peakSpeed, peak);
    out.projectionSweeps = sweeps;
    out.projectionResidual = residual;
    if (capped) out.projectionCapped = true;
    return true;
}

}  // namespace

StepResult step(Field& field, double dt, const std::vector<HeatSource>& sources,
                const Params& params, Account& account) {
    StepResult out;
    out.time = field.time;
    for (const HeatSource& s : sources) out.heatRelease += s.power;
    if (dt <= 0 || field.empty()) return out;

    Work work;
    work.size(field.grid);

    double remaining = dt;
    double h = std::min(dt, std::max(params.maxSubstep, 1e-9));
    int budget = std::max(params.maxSubsteps, 1);
    while (remaining > 1e-12 * dt && budget-- > 0) {
        h = std::min(h, remaining);
        if (substep(field, h, sources, params, work, account, out)) {
            remaining -= h;
            ++out.substeps;
            // Creep back up, so one transient does not pin the step for the whole
            // run. The same 1.5 `fire.cpp` settled on, and for the same reason:
            // doubling oscillates against the rejection test.
            h = std::min(h * 1.5, params.maxSubstep);
        } else {
            ++out.rejections;
            h *= 0.5;
        }
    }
    // **Say so when the budget ran out**, rather than returning a field that is
    // short of the time it was asked for and looks like any other. `fire.cpp`'s
    // substep loop has the same shape and the same exposure; mutation testing found
    // it here, because a defect that drives the rejection test into a corner does
    // not fail -- it silently under-advances, and a caller integrating a fire
    // against a conduction solve would never know its two halves had come apart.
    out.incomplete = remaining > 1e-12 * dt;
    out.time = field.time;
    return out;
}

}  // namespace sim::les
