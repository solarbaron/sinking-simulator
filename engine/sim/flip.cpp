// SPDX-License-Identifier: MIT
#include "flip.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>

namespace sim::flip {
namespace {

// Compensated (Kahan) summation, on exactly the terms `les.cpp`'s `Accumulator`
// states them: every assertion this file is under is a conservation assertion, and
// a test cannot be tighter than its own arithmetic. The mass total is the one place
// it is *not* needed -- the same numbers are summed in the same order at every
// step, so a naive sum would be equally constant -- and it is used there anyway,
// because "exactly 0.0" is a claim about the sum and not about the summand.
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

// Floor division onto the tile lattice. `c / kTile` truncates towards zero, which
// puts cells -1..-3 in tile 0 alongside cells 0..3 -- eight cells in a four-cell
// tile, and every one of them aliased onto the wrong storage. Negative cell
// indices are ordinary here: the halo reaches outside the domain by construction.
int tileOf(int c) { return c >= 0 ? c / kTile : -((-c + kTile - 1) / kTile); }

// Tile coordinates are packed into 21 bits each with a bias, which covers
// +/-1048576 tiles -- +/-4.2 million cells, 210 km at a 5 cm cell. `validate`
// reports a field that would overflow it rather than letting two tiles alias.
inline constexpr std::int32_t kTileBias = 1 << 20;
inline constexpr std::int32_t kTileLimit = kTileBias - 2;

std::uint64_t packTile(int ti, int tj, int tk) {
    const std::uint64_t a = static_cast<std::uint64_t>(ti + kTileBias);
    const std::uint64_t b = static_cast<std::uint64_t>(tj + kTileBias);
    const std::uint64_t c = static_cast<std::uint64_t>(tk + kTileBias);
    return (a << 42) | (b << 21) | c;
}

int localIndex(int lx, int ly, int lz) { return (lz * kTile + ly) * kTile + lx; }

}  // namespace

// ---------------------------------------------------------------------------
// Grid, Field, Account
// ---------------------------------------------------------------------------

int Grid::cellOf(int axis, double x) const {
    if (!(h > 0)) return 0;
    return static_cast<int>(std::floor((x - lo[axis]) / h));
}

double Field::totalMass() const {
    Accumulator total;
    for (const Particle& p : particles) total.add(p.mass);
    return total.total();
}

double Field::volumeAt(double density) const {
    return density > 0 ? totalMass() / density : 0.0;
}

void resetAccount(const Field& field, Account& account) {
    account = Account{};
    account.initialMass = field.totalMass();
    account.mass = account.initialMass;
    account.initialParticles = static_cast<long long>(field.particles.size());
    account.particles = account.initialParticles;
}

double Account::massResidualFraction() const {
    const double scale = std::max({std::abs(mass), std::abs(initialMass), 1e-12});
    return massResidual() / scale;
}

// ---------------------------------------------------------------------------
// The kernel
// ---------------------------------------------------------------------------

int splineWeights(double s, double w[3], double offset[3]) {
    const int base = static_cast<int>(std::floor(s - 0.5));
    // `f` is the particle's offset from the *middle* node, in cells, and lies in
    // [-0.5, 0.5]. Writing the weights around the middle node rather than around
    // the base is what makes the three of them sum to one identically instead of
    // to one plus a rounding that depends on where the particle is.
    const double f = s - (static_cast<double>(base) + 1.0);
    w[0] = 0.5 * (0.5 - f) * (0.5 - f);
    w[1] = 0.75 - f * f;
    w[2] = 0.5 * (0.5 + f) * (0.5 + f);
    offset[0] = -1.0 - f;
    offset[1] = -f;
    offset[2] = 1.0 - f;
    return base;
}

double kernelSecondMoment(double s, double h) {
    double w[3], offset[3];
    splineWeights(s, w, offset);
    double total = 0;
    for (int m = 0; m < 3; ++m) total += w[m] * offset[m] * offset[m];
    return total * h * h;
}

// ---------------------------------------------------------------------------
// Seeding and the quiescent closed forms
// ---------------------------------------------------------------------------

void seedBox(Field& field, const double lo[3], const double hi[3], int perAxis,
             double density) {
    const Grid& g = field.grid;
    if (g.empty() || perAxis < 1 || !(density > 0)) return;

    int cellLo[3], cellHi[3];
    for (int a = 0; a < 3; ++a) {
        cellLo[a] = std::max(g.cellOf(a, lo[a]), 0);
        // The cell the *upper* face falls in, exclusive. A box whose top lands
        // exactly on a cell boundary must not seed the cell above it, so the
        // ceiling is taken on the fractional coordinate rather than on the cell
        // index -- the same "on the boundary" question `section.cpp` got wrong by
        // asking it two different ways in two places.
        const double top = (hi[a] - g.lo[a]) / g.h;
        cellHi[a] = std::min(static_cast<int>(std::ceil(top - 1e-12)), g.n[a]);
    }

    const double share = g.cellVolume() * density /
                         (static_cast<double>(perAxis) * perAxis * perAxis);
    const double step = 1.0 / static_cast<double>(perAxis);
    for (int k = cellLo[2]; k < cellHi[2]; ++k)
        for (int j = cellLo[1]; j < cellHi[1]; ++j)
            for (int i = cellLo[0]; i < cellHi[0]; ++i)
                for (int c = 0; c < perAxis; ++c)
                    for (int b = 0; b < perAxis; ++b)
                        for (int a = 0; a < perAxis; ++a) {
                            const int index[3] = {i, j, k};
                            const int sub[3] = {a, b, c};
                            Particle p;
                            bool in = true;
                            for (int d = 0; d < 3; ++d) {
                                p.position[d] = g.lo[d] +
                                                (static_cast<double>(index[d]) +
                                                 (static_cast<double>(sub[d]) + 0.5) * step) * g.h;
                                if (p.position[d] < lo[d] || p.position[d] >= hi[d]) in = false;
                            }
                            if (!in) continue;
                            p.mass = share;
                            field.particles.push_back(p);
                        }
}

bool setTotalMass(Field& field, double mass) {
    const std::size_t n = field.particles.size();
    if (n == 0 || !(mass > 0)) return false;
    const double share = mass / static_cast<double>(n);
    Accumulator sum;
    for (std::size_t i = 0; i + 1 < n; ++i) {
        field.particles[i].mass = share;
        sum.add(share);
    }
    // The remainder, not the share. `n` equal shares do not add to `mass` in
    // floating point; the last particle absorbs the difference so that the total
    // is the number that was asked for. Same discipline as `fire.cpp`'s
    // `layerSplit` and `les.cpp`'s two-layer reconstruction.
    field.particles[n - 1].mass = mass - sum.total();
    return true;
}

double quiescentLevel(const Field& field, double density, double planArea, double floorZ) {
    if (!(density > 0) || !(planArea > 0)) return floorZ;
    return floorZ + field.totalMass() / (density * planArea);
}

double kineticEnergy(const Field& field) {
    Accumulator total;
    for (const Particle& p : field.particles) {
        const double v2 = p.velocity[0] * p.velocity[0] + p.velocity[1] * p.velocity[1] +
                          p.velocity[2] * p.velocity[2];
        total.add(0.5 * p.mass * v2);
    }
    return total.total();
}

void linearMomentum(const Field& field, double out[3]) {
    Accumulator total[3];
    for (const Particle& p : field.particles)
        for (int a = 0; a < 3; ++a) total[a].add(p.mass * p.velocity[a]);
    for (int a = 0; a < 3; ++a) out[a] = total[a].total();
}

void angularMomentum(const Field& field, const double about[3], double out[3]) {
    Accumulator total[3];
    for (const Particle& p : field.particles) {
        const double r[3] = {p.position[0] - about[0], p.position[1] - about[1],
                             p.position[2] - about[2]};
        total[0].add(p.mass * (r[1] * p.velocity[2] - r[2] * p.velocity[1]));
        total[1].add(p.mass * (r[2] * p.velocity[0] - r[0] * p.velocity[2]));
        total[2].add(p.mass * (r[0] * p.velocity[1] - r[1] * p.velocity[0]));
    }
    for (int a = 0; a < 3; ++a) out[a] = total[a].total();
}

void centroid(const Field& field, double out[3]) {
    Accumulator total[3];
    Accumulator mass;
    for (const Particle& p : field.particles) {
        mass.add(p.mass);
        for (int a = 0; a < 3; ++a) total[a].add(p.mass * p.position[a]);
    }
    const double m = mass.total();
    for (int a = 0; a < 3; ++a) out[a] = m > 0 ? total[a].total() / m : 0.0;
}

std::vector<std::string> validate(const Field& field, const Params& params) {
    std::vector<std::string> problems;
    const Grid& g = field.grid;
    if (!(g.h > 0)) problems.push_back("cell size is not positive");
    for (int a = 0; a < 3; ++a)
        if (g.n[a] <= 0) problems.push_back("domain has no cells along axis " +
                                            std::to_string(a));
    if (!(params.density > 0)) problems.push_back("density is not positive");
    if (params.flipBlend < 0 || params.flipBlend > 1)
        problems.push_back("flipBlend is outside [0, 1]");
    if (params.extrapolationDepth < 2)
        problems.push_back("extrapolationDepth below 2 leaves the transfer reading "
                           "faces nothing filled");
    if (params.courant <= 0 || params.courant >= static_cast<double>(kHalo))
        problems.push_back("courant must be positive and below the halo depth");

    int outside = 0, massless = 0;
    for (const Particle& p : field.particles) {
        bool in = true;
        for (int a = 0; a < 3; ++a)
            if (p.position[a] < g.lo[a] || p.position[a] > g.hi(a)) in = false;
        if (!in) ++outside;
        if (!(p.mass > 0)) ++massless;
    }
    if (outside > 0)
        problems.push_back(std::to_string(outside) + " particles are outside the domain");
    if (massless > 0)
        problems.push_back(std::to_string(massless) + " particles carry no mass");

    if (!g.empty() && g.h > 0) {
        for (int a = 0; a < 3; ++a) {
            const int high = tileOf(g.n[a] + kHalo);
            const int low = tileOf(-kHalo - 1);
            if (high > kTileLimit || low < -kTileLimit)
                problems.push_back("the domain needs more tiles than the key packs");
        }
    }
    return problems;
}

// ---------------------------------------------------------------------------
// Solver -- the sparse structure
// ---------------------------------------------------------------------------

int Solver::tileIndex(int ti, int tj, int tk) const {
    const auto it = map_.find(packTile(ti, tj, tk));
    return it == map_.end() ? -1 : it->second;
}

int Solver::index(int i, int j, int k) const {
    const int t = tileIndex(tileOf(i), tileOf(j), tileOf(k));
    if (t < 0) return -1;
    const int lx = i - tileOf(i) * kTile;
    const int ly = j - tileOf(j) * kTile;
    const int lz = k - tileOf(k) * kTile;
    return t * kTileCells + localIndex(lx, ly, lz);
}

std::size_t Solver::bytes() const {
    std::size_t total = 0;
    total += tileKey_.size() * sizeof(std::int32_t);
    total += cell_.size() * sizeof(std::uint8_t);
    total += nb_.size() * sizeof(std::int32_t);
    for (int a = 0; a < 3; ++a) {
        total += vel_[a].size() * sizeof(double);
        total += old_[a].size() * sizeof(double);
        total += wgt_[a].size() * sizeof(double);
        total += face_[a].size() * sizeof(std::uint8_t);
    }
    total += pressure_.size() * sizeof(double);
    total += fluid_.size() * sizeof(std::int32_t);
    total += fluidOf_.size() * sizeof(std::int32_t);
    total += fnb_.size() * sizeof(std::int32_t);
    total += (diag_.size() + rhs_.size() + sol_.size() + res_.size() + dir_.size() +
              mat_.size() + pre_.size()) * sizeof(double);
    return total;
}

void Solver::rebuild(const Field& field, const Params& params) {
    (void)params;
    grid_ = field.grid;
    map_.clear();
    tileKey_.clear();
    cell_.clear();
    nb_.clear();
    for (int a = 0; a < 3; ++a) {
        vel_[a].clear();
        old_[a].clear();
        wgt_[a].clear();
        face_[a].clear();
    }
    pressure_.clear();
    fluid_.clear();
    fluidOf_.clear();
    fnb_.clear();
    diag_.clear();
    rhs_.clear();
    sol_.clear();
    res_.clear();
    dir_.clear();
    mat_.clear();
    pre_.clear();
    order_.clear();
    tileStart_.clear();
    particleTile_.clear();
    outside_ = 0;
    if (grid_.empty() || field.particles.empty()) return;

    // --- the tiles the particles themselves land in --------------------------
    //
    // Insertion order is the order the particles are visited in, which is the
    // caller's array order. Everything downstream iterates the tile *vector* and
    // never the hash map, so the whole solve is a pure function of the particle
    // array and not of the map's bucket layout -- which is what lets the same
    // water in two differently sized rooms come out bit for bit identical.
    const std::size_t count = field.particles.size();
    particleTile_.resize(count);
    std::vector<std::int32_t> particleLocal(count);
    for (std::size_t p = 0; p < count; ++p) {
        const Particle& q = field.particles[p];
        const int c[3] = {grid_.cellOf(0, q.position[0]), grid_.cellOf(1, q.position[1]),
                          grid_.cellOf(2, q.position[2])};
        const int t[3] = {tileOf(c[0]), tileOf(c[1]), tileOf(c[2])};
        const std::uint64_t key = packTile(t[0], t[1], t[2]);
        auto it = map_.find(key);
        int slot;
        if (it == map_.end()) {
            slot = static_cast<int>(tileKey_.size() / 3);
            map_.emplace(key, slot);
            tileKey_.push_back(t[0]);
            tileKey_.push_back(t[1]);
            tileKey_.push_back(t[2]);
        } else {
            slot = it->second;
        }
        particleTile_[p] = slot;
        particleLocal[p] = static_cast<std::int32_t>(
            localIndex(c[0] - t[0] * kTile, c[1] - t[1] * kTile, c[2] - t[2] * kTile));
    }

    const int fluidTiles = static_cast<int>(tileKey_.size() / 3);

    // --- which cells hold particles, and how far each tile has to reach -------
    mark_.assign(static_cast<std::size_t>(fluidTiles) * kTileCells, 0u);
    // Per fluid tile, the extreme local coordinates its particles occupy. The
    // halo is then taken from those rather than from the whole tile, so a tile
    // holding one particle in a corner does not pull in the neighbours on the
    // other three sides.
    std::vector<std::int32_t> lowLocal(static_cast<std::size_t>(fluidTiles) * 3, kTile);
    std::vector<std::int32_t> highLocal(static_cast<std::size_t>(fluidTiles) * 3, -1);
    for (std::size_t p = 0; p < count; ++p) {
        const int slot = particleTile_[p];
        const int local = particleLocal[p];
        mark_[static_cast<std::size_t>(slot) * kTileCells + static_cast<std::size_t>(local)] = 1u;
        const int lx = local % kTile;
        const int ly = (local / kTile) % kTile;
        const int lz = local / (kTile * kTile);
        const int l[3] = {lx, ly, lz};
        for (int a = 0; a < 3; ++a) {
            const std::size_t at = static_cast<std::size_t>(slot) * 3 + static_cast<std::size_t>(a);
            lowLocal[at] = std::min(lowLocal[at], l[a]);
            highLocal[at] = std::max(highLocal[at], l[a]);
        }
    }

    // --- the halo ------------------------------------------------------------
    for (int t = 0; t < fluidTiles; ++t) {
        int from[3], to[3];
        for (int a = 0; a < 3; ++a) {
            const std::size_t at = static_cast<std::size_t>(t) * 3 + static_cast<std::size_t>(a);
            const int base = tileKey_[at] * kTile;
            from[a] = tileOf(base + lowLocal[at] - kHalo);
            to[a] = tileOf(base + highLocal[at] + kHalo);
        }
        for (int tk = from[2]; tk <= to[2]; ++tk)
            for (int tj = from[1]; tj <= to[1]; ++tj)
                for (int ti = from[0]; ti <= to[0]; ++ti) {
                    const std::uint64_t key = packTile(ti, tj, tk);
                    if (map_.find(key) != map_.end()) continue;
                    map_.emplace(key, static_cast<int>(tileKey_.size() / 3));
                    tileKey_.push_back(ti);
                    tileKey_.push_back(tj);
                    tileKey_.push_back(tk);
                }
    }

    // --- storage -------------------------------------------------------------
    const int tileCount = static_cast<int>(tileKey_.size() / 3);
    const std::size_t cells = static_cast<std::size_t>(tileCount) * kTileCells;
    cell_.assign(cells, static_cast<std::uint8_t>(Cell::Air));
    pressure_.assign(cells, 0.0);
    nb_.assign(cells * 6, -1);
    for (int a = 0; a < 3; ++a) {
        vel_[a].assign(cells, 0.0);
        old_[a].assign(cells, 0.0);
        wgt_[a].assign(cells, 0.0);
        face_[a].assign(cells, static_cast<std::uint8_t>(Face::Air));
    }

    for (int t = 0; t < tileCount; ++t) {
        const int base[3] = {tileKey_[static_cast<std::size_t>(t) * 3] * kTile,
                             tileKey_[static_cast<std::size_t>(t) * 3 + 1] * kTile,
                             tileKey_[static_cast<std::size_t>(t) * 3 + 2] * kTile};
        for (int lz = 0; lz < kTile; ++lz)
            for (int ly = 0; ly < kTile; ++ly)
                for (int lx = 0; lx < kTile; ++lx) {
                    const std::size_t flat = static_cast<std::size_t>(t) * kTileCells +
                                             static_cast<std::size_t>(localIndex(lx, ly, lz));
                    if (!grid_.inside(base[0] + lx, base[1] + ly, base[2] + lz))
                        cell_[flat] = static_cast<std::uint8_t>(Cell::Solid);
                }
    }
    for (std::size_t flat = 0; flat < mark_.size(); ++flat) {
        if (!mark_[flat]) continue;
        if (cell_[flat] == static_cast<std::uint8_t>(Cell::Solid)) {
            // A particle outside the solid box. Nothing here creates one -- the
            // advection clamps -- so this is a seeded field that was never
            // validated, and it is counted rather than silently treated as water.
            ++outside_;
            continue;
        }
        cell_[flat] = static_cast<std::uint8_t>(Cell::Fluid);
    }

    buildNeighbours();
    classifyFaces();
    buildFluidList();

    // --- particles bucketed by tile ------------------------------------------
    //
    // The transfer's inner loop needs the twenty-seven tiles around a particle's
    // own. Looked up per particle that is eighty-one hash probes each; looked up
    // per *tile*, with the particles gathered behind it, it is twenty-seven per
    // tile and arithmetic thereafter. A counting sort, so the order within a tile
    // stays the caller's order and the whole thing stays deterministic.
    tileStart_.assign(static_cast<std::size_t>(tileCount) + 1, 0);
    for (std::size_t p = 0; p < count; ++p)
        ++tileStart_[static_cast<std::size_t>(particleTile_[p]) + 1];
    for (std::size_t t = 0; t < static_cast<std::size_t>(tileCount); ++t)
        tileStart_[t + 1] += tileStart_[t];
    order_.resize(count);
    std::vector<std::int32_t> cursor(tileStart_.begin(), tileStart_.end() - 1);
    for (std::size_t p = 0; p < count; ++p) {
        const std::size_t slot = static_cast<std::size_t>(particleTile_[p]);
        order_[static_cast<std::size_t>(cursor[slot]++)] = static_cast<std::int32_t>(p);
    }
}

void Solver::buildNeighbours() {
    const int tileCount = static_cast<int>(tileKey_.size() / 3);
    for (int t = 0; t < tileCount; ++t) {
        const int tc[3] = {tileKey_[static_cast<std::size_t>(t) * 3],
                           tileKey_[static_cast<std::size_t>(t) * 3 + 1],
                           tileKey_[static_cast<std::size_t>(t) * 3 + 2]};
        int adjacent[6];
        for (int a = 0; a < 3; ++a)
            for (int side = 0; side < 2; ++side) {
                int n[3] = {tc[0], tc[1], tc[2]};
                n[a] += side == 0 ? -1 : 1;
                adjacent[a * 2 + side] = tileIndex(n[0], n[1], n[2]);
            }
        for (int lz = 0; lz < kTile; ++lz)
            for (int ly = 0; ly < kTile; ++ly)
                for (int lx = 0; lx < kTile; ++lx) {
                    const std::size_t flat = static_cast<std::size_t>(t) * kTileCells +
                                             static_cast<std::size_t>(localIndex(lx, ly, lz));
                    for (int a = 0; a < 3; ++a)
                        for (int side = 0; side < 2; ++side) {
                            int n[3] = {lx, ly, lz};
                            n[a] += side == 0 ? -1 : 1;
                            int host = t;
                            if (n[a] < 0) {
                                host = adjacent[a * 2];
                                n[a] = kTile - 1;
                            } else if (n[a] >= kTile) {
                                host = adjacent[a * 2 + 1];
                                n[a] = 0;
                            }
                            nb_[flat * 6 + static_cast<std::size_t>(a * 2 + side)] =
                                host < 0 ? -1
                                         : static_cast<std::int32_t>(
                                               static_cast<std::size_t>(host) * kTileCells +
                                               static_cast<std::size_t>(
                                                   localIndex(n[0], n[1], n[2])));
                        }
                }
    }
}

void Solver::classifyFaces() {
    for (std::size_t flat = 0; flat < cell_.size(); ++flat) {
        if (cell_[flat] != static_cast<std::uint8_t>(Cell::Fluid)) continue;
        for (int a = 0; a < 3; ++a) {
            // The face stored *on* this cell is its low face, and it separates
            // this cell from the neighbour on the minus side. The face stored on
            // the plus neighbour is this cell's high face. Getting that the wrong
            // way round makes every wall one cell thick in the wrong direction and
            // a hydrostatic column drift; it is the single easiest sign error in a
            // staggered grid and it is why the two are written out separately.
            const std::int32_t minus = nb_[flat * 6 + static_cast<std::size_t>(a * 2)];
            const bool minusSolid =
                minus < 0 || cell_[static_cast<std::size_t>(minus)] ==
                                 static_cast<std::uint8_t>(Cell::Solid);
            face_[a][flat] = static_cast<std::uint8_t>(minusSolid ? Face::Solid : Face::Fluid);

            const std::int32_t plus = nb_[flat * 6 + static_cast<std::size_t>(a * 2 + 1)];
            if (plus >= 0) {
                const bool plusSolid = cell_[static_cast<std::size_t>(plus)] ==
                                       static_cast<std::uint8_t>(Cell::Solid);
                face_[a][static_cast<std::size_t>(plus)] =
                    static_cast<std::uint8_t>(plusSolid ? Face::Solid : Face::Fluid);
            }
        }
    }
}

void Solver::buildFluidList() {
    fluidOf_.assign(cell_.size(), -1);
    for (std::size_t flat = 0; flat < cell_.size(); ++flat)
        if (cell_[flat] == static_cast<std::uint8_t>(Cell::Fluid)) {
            fluidOf_[flat] = static_cast<std::int32_t>(fluid_.size());
            fluid_.push_back(static_cast<std::int32_t>(flat));
        }

    const std::size_t n = fluid_.size();
    fnb_.assign(n * 6, -2);
    diag_.assign(n, 0.0);
    rhs_.assign(n, 0.0);
    sol_.assign(n, 0.0);
    res_.assign(n, 0.0);
    dir_.assign(n, 0.0);
    mat_.assign(n, 0.0);
    pre_.assign(n, 0.0);
    singular_ = true;
    for (std::size_t s = 0; s < n; ++s) {
        const std::size_t flat = static_cast<std::size_t>(fluid_[s]);
        int nonSolid = 0;
        for (int d = 0; d < 6; ++d) {
            const std::int32_t other = nb_[flat * 6 + static_cast<std::size_t>(d)];
            if (other < 0 || cell_[static_cast<std::size_t>(other)] ==
                                 static_cast<std::uint8_t>(Cell::Solid)) {
                fnb_[s * 6 + static_cast<std::size_t>(d)] = -2;   // solid: the term drops
                continue;
            }
            ++nonSolid;
            if (cell_[static_cast<std::size_t>(other)] ==
                static_cast<std::uint8_t>(Cell::Fluid)) {
                fnb_[s * 6 + static_cast<std::size_t>(d)] = fluidOf_[static_cast<std::size_t>(other)];
            } else {
                // Air: a Dirichlet zero. It carries the diagonal and no
                // off-diagonal, and its presence is what makes the system
                // non-singular -- a body of water with a free surface anywhere
                // determines its own pressure, one with none does not.
                fnb_[s * 6 + static_cast<std::size_t>(d)] = -1;
                singular_ = false;
            }
        }
        diag_[s] = static_cast<double>(nonSolid);
    }
    if (n == 0) singular_ = false;
}

// ---------------------------------------------------------------------------
// Queries
// ---------------------------------------------------------------------------

Cell Solver::cellAt(int i, int j, int k) const {
    const int flat = index(i, j, k);
    if (flat < 0) return grid_.inside(i, j, k) ? Cell::Air : Cell::Solid;
    return static_cast<Cell>(cell_[static_cast<std::size_t>(flat)]);
}

Face Solver::faceAt(int axis, int i, int j, int k) const {
    const int flat = index(i, j, k);
    if (flat < 0) return Face::Air;
    return static_cast<Face>(face_[axis][static_cast<std::size_t>(flat)]);
}

double Solver::faceVelocity(int axis, int i, int j, int k) const {
    const int flat = index(i, j, k);
    return flat < 0 ? 0.0 : vel_[axis][static_cast<std::size_t>(flat)];
}

double Solver::faceMass(int axis, int i, int j, int k) const {
    const int flat = index(i, j, k);
    return flat < 0 ? 0.0 : wgt_[axis][static_cast<std::size_t>(flat)];
}

double Solver::totalFaceMass(int axis) const {
    Accumulator total;
    for (double w : wgt_[axis]) total.add(w);
    return total.total();
}

double Solver::pressureAt(int i, int j, int k) const {
    const int flat = index(i, j, k);
    return flat < 0 ? 0.0 : pressure_[static_cast<std::size_t>(flat)];
}

double Solver::divergenceAt(int i, int j, int k) const {
    const int flat = index(i, j, k);
    if (flat < 0 || !(grid_.h > 0)) return 0.0;
    double total = 0;
    for (int a = 0; a < 3; ++a) {
        const std::int32_t plus =
            nb_[static_cast<std::size_t>(flat) * 6 + static_cast<std::size_t>(a * 2 + 1)];
        const double hi = plus < 0 ? 0.0 : vel_[a][static_cast<std::size_t>(plus)];
        total += hi - vel_[a][static_cast<std::size_t>(flat)];
    }
    return total / grid_.h;
}

double Solver::maxDivergence() const {
    double worst = 0;
    for (std::size_t s = 0; s < fluid_.size(); ++s) {
        const std::size_t flat = static_cast<std::size_t>(fluid_[s]);
        double total = 0;
        for (int a = 0; a < 3; ++a) {
            const std::int32_t plus = nb_[flat * 6 + static_cast<std::size_t>(a * 2 + 1)];
            const double hi = plus < 0 ? 0.0 : vel_[a][static_cast<std::size_t>(plus)];
            total += hi - vel_[a][flat];
        }
        worst = std::max(worst, std::abs(total) / grid_.h);
    }
    return worst;
}

// ---------------------------------------------------------------------------
// Particle -> grid
// ---------------------------------------------------------------------------

namespace {

// The twenty-seven tiles around tile `t`, as flat tile indices or -1. Everything
// a particle in `t` can reach: the kernel spans +/-2 cells and the advection
// midpoint one more, and `kTile` is 4, so no node is ever two tiles away.
struct TileHalo {
    int index[27];
    int base[3];
};

}  // namespace

void Solver::gatherHalo(int tile, void* out) const {
    TileHalo& halo = *static_cast<TileHalo*>(out);
    const int tc[3] = {tileKey_[static_cast<std::size_t>(tile) * 3],
                       tileKey_[static_cast<std::size_t>(tile) * 3 + 1],
                       tileKey_[static_cast<std::size_t>(tile) * 3 + 2]};
    for (int a = 0; a < 3; ++a) halo.base[a] = tc[a] * kTile;
    for (int dz = -1; dz <= 1; ++dz)
        for (int dy = -1; dy <= 1; ++dy)
            for (int dx = -1; dx <= 1; ++dx)
                halo.index[((dz + 1) * 3 + (dy + 1)) * 3 + (dx + 1)] =
                    tileIndex(tc[0] + dx, tc[1] + dy, tc[2] + dz);
}

// The flat cell index of cell (i, j, k) given a tile halo it is known to fall in,
// or -1. Pure arithmetic: no hash probe.
int Solver::haloIndex(const void* halo, int i, int j, int k) const {
    const TileHalo& h = *static_cast<const TileHalo*>(halo);
    const int c[3] = {i, j, k};
    int offset[3], local[3];
    for (int a = 0; a < 3; ++a) {
        const int rel = c[a] - h.base[a];
        // `rel` is in [-kTile, 2 kTile) for anything the kernel can reach, so the
        // tile offset is one of the three the halo carries. Anything else is a
        // caller error and comes back absent rather than aliased.
        if (rel < -kTile || rel >= 2 * kTile) return -1;
        offset[a] = rel < 0 ? -1 : (rel < kTile ? 0 : 1);
        local[a] = rel - offset[a] * kTile;
    }
    const int tile = h.index[((offset[2] + 1) * 3 + (offset[1] + 1)) * 3 + (offset[0] + 1)];
    if (tile < 0) return -1;
    return tile * kTileCells + localIndex(local[0], local[1], local[2]);
}

void Solver::transferToGrid(const Field& field, const Params& params) {
    for (int a = 0; a < 3; ++a) {
        std::fill(vel_[a].begin(), vel_[a].end(), 0.0);
        std::fill(wgt_[a].begin(), wgt_[a].end(), 0.0);
    }
    if (tileKey_.empty()) return;
    const double h = grid_.h;
    const int tileCount = static_cast<int>(tileKey_.size() / 3);

    TileHalo halo;
    for (int t = 0; t < tileCount; ++t) {
        const std::size_t from = static_cast<std::size_t>(tileStart_[static_cast<std::size_t>(t)]);
        const std::size_t to = static_cast<std::size_t>(tileStart_[static_cast<std::size_t>(t) + 1]);
        if (from == to) continue;
        gatherHalo(t, &halo);
        for (std::size_t o = from; o < to; ++o) {
            const Particle& p = field.particles[static_cast<std::size_t>(order_[o])];
            for (int a = 0; a < 3; ++a) {
                double w[3][3], offset[3][3];
                int base[3];
                for (int b = 0; b < 3; ++b) {
                    // The component's own axis carries its nodes on the cell
                    // faces; the other two carry them at the cell centres. Half a
                    // cell, three times, and reading it off the wrong axis is the
                    // classic staggered-grid transposition.
                    const double s = (p.position[b] - grid_.lo[b]) / h - (b == a ? 0.0 : 0.5);
                    base[b] = splineWeights(s, w[b], offset[b]);
                }
                for (int kk = 0; kk < 3; ++kk)
                    for (int jj = 0; jj < 3; ++jj)
                        for (int ii = 0; ii < 3; ++ii) {
                            const double weight = w[0][ii] * w[1][jj] * w[2][kk];
                            if (weight == 0.0) continue;
                            const int flat = haloIndex(&halo, base[0] + ii, base[1] + jj,
                                                       base[2] + kk);
                            if (flat < 0) continue;
                            double value = p.velocity[a];
                            if (params.affine)
                                value += (p.affine[a * 3 + 0] * offset[0][ii] +
                                          p.affine[a * 3 + 1] * offset[1][jj] +
                                          p.affine[a * 3 + 2] * offset[2][kk]) * h;
                            const double share = weight * p.mass;
                            wgt_[a][static_cast<std::size_t>(flat)] += share;
                            vel_[a][static_cast<std::size_t>(flat)] += share * value;
                        }
            }
        }
    }

    for (int a = 0; a < 3; ++a)
        for (std::size_t flat = 0; flat < vel_[a].size(); ++flat) {
            if (face_[a][flat] != static_cast<std::uint8_t>(Face::Fluid)) {
                vel_[a][flat] = 0.0;
                continue;
            }
            // Every face beside a fluid cell has strictly positive weight -- a
            // particle in that cell is within one cell of the face and the kernel
            // reaches one and a half. So this division has no small-denominator
            // branch, and `tests/test_flip.cpp` asserts the premise rather than
            // trusting it.
            const double mass = wgt_[a][flat];
            vel_[a][flat] = mass > 0 ? vel_[a][flat] / mass : 0.0;
        }
}

void Solver::saveGrid() {
    for (int a = 0; a < 3; ++a) old_[a] = vel_[a];
}

bool Solver::setFaceVelocity(int axis, int i, int j, int k, double value) {
    const int flat = index(i, j, k);
    if (flat < 0) return false;
    vel_[axis][static_cast<std::size_t>(flat)] = value;
    return true;
}

// ---------------------------------------------------------------------------
// Body force, projection, extrapolation
// ---------------------------------------------------------------------------

void Solver::addBodyForce(double dt, const Params& params) {
    for (int a = 0; a < 3; ++a) {
        const double delta = params.gravity[a] * dt;
        if (delta == 0.0) continue;
        for (std::size_t flat = 0; flat < vel_[a].size(); ++flat)
            if (face_[a][flat] == static_cast<std::uint8_t>(Face::Fluid))
                vel_[a][flat] += delta;
    }
}

double Solver::dotFluid(const std::vector<double>& a, const std::vector<double>& b) const {
    Accumulator total;
    for (std::size_t s = 0; s < a.size(); ++s) total.add(a[s] * b[s]);
    return total.total();
}

void Solver::project(double dt, const Params& params) {
    iterations_ = 0;
    residual_ = 0;
    capped_ = false;
    std::fill(pressure_.begin(), pressure_.end(), 0.0);
    const std::size_t n = fluid_.size();
    if (n == 0 || !(dt > 0) || !(params.density > 0) || !(grid_.h > 0)) return;

    const double h = grid_.h;
    const double scale = params.density * h * h / dt;
    double rhsScale = 0;
    // Set by any non-finite right-hand side or residual seen below. It is a flag
    // and not a count because one is already fatal to the answer.
    bool nonFinite = false;
    for (std::size_t s = 0; s < n; ++s) {
        const std::size_t flat = static_cast<std::size_t>(fluid_[s]);
        double divergence = 0;
        for (int a = 0; a < 3; ++a) {
            const std::int32_t plus = nb_[flat * 6 + static_cast<std::size_t>(a * 2 + 1)];
            const double hi = plus < 0 ? 0.0 : vel_[a][static_cast<std::size_t>(plus)];
            divergence += hi - vel_[a][flat];
        }
        divergence /= h;
        rhs_[s] = -divergence * scale;
        sol_[s] = 0.0;
        // `std::max(rhsScale, x)` is `x < rhsScale ? x : rhsScale`, so a NaN here
        // is dropped and `rhsScale` comes back as the largest of the *finite*
        // cells. That much is wanted, and this is deliberately **not** where the
        // non-finite flag is raised: `res_` is a copy of `rhs_`, so
        // `worstResidual` below sees every entry this loop could, and a second
        // test here would be a guard no test can tell from its own absence.
        // Written that way first and removed on the evidence -- a mutation sweep
        // scored it a survivor, because taking it out changes nothing observable.
        rhsScale = std::max(rhsScale, std::abs(rhs_[s]));
    }

    // A region with no free surface anywhere determines its pressure only up to a
    // constant, and accepts only a right-hand side that sums to zero. The sum is
    // zero here by telescoping -- every interior face appears twice with opposite
    // signs and every wall face is zero -- but only in exact arithmetic, so the
    // mean is removed explicitly. Same compatibility discipline as `les.cpp`'s
    // all-Neumann heating source, arrived at from the other direction.
    const auto removeMean = [&](std::vector<double>& v) {
        if (!singular_ || n == 0) return;
        Accumulator sum;
        for (std::size_t s = 0; s < n; ++s) sum.add(v[s]);
        const double mean = sum.total() / static_cast<double>(n);
        for (std::size_t s = 0; s < n; ++s) v[s] -= mean;
    };
    const auto apply = [&](const std::vector<double>& x, std::vector<double>& y) {
        for (std::size_t s = 0; s < n; ++s) {
            double total = diag_[s] * x[s];
            for (int d = 0; d < 6; ++d) {
                const std::int32_t other = fnb_[s * 6 + static_cast<std::size_t>(d)];
                if (other >= 0) total -= x[static_cast<std::size_t>(other)];
            }
            y[s] = total;
        }
        removeMean(y);
    };

    removeMean(rhs_);
    res_ = rhs_;
    // The largest |residual|, with the same NaN discipline as `rhsScale` above and
    // for the same reason: this number decides whether the solve runs at all.
    const auto worstResidual = [&] {
        double worst = 0;
        for (std::size_t s = 0; s < n; ++s) {
            if (std::isfinite(res_[s]))
                worst = std::max(worst, std::abs(res_[s]));
            else
                nonFinite = true;
        }
        return worst;
    };
    double best = worstResidual();
    const double tolerance = params.projectionTolerance * rhsScale;
    if (best > tolerance) {
        for (std::size_t s = 0; s < n; ++s)
            pre_[s] = diag_[s] > 0 ? res_[s] / diag_[s] : 0.0;
        removeMean(pre_);
        dir_ = pre_;
        double rho = dotFluid(res_, pre_);
        const int limit = std::max(params.projectionIterations, 1);
        for (int it = 0; it < limit; ++it) {
            apply(dir_, mat_);
            const double denominator = dotFluid(dir_, mat_);
            if (!(std::abs(denominator) > 0)) break;
            const double alpha = rho / denominator;
            for (std::size_t s = 0; s < n; ++s) {
                sol_[s] += alpha * dir_[s];
                res_[s] -= alpha * mat_[s];
            }
            ++iterations_;
            best = worstResidual();
            if (best <= tolerance) break;
            for (std::size_t s = 0; s < n; ++s)
                pre_[s] = diag_[s] > 0 ? res_[s] / diag_[s] : 0.0;
            removeMean(pre_);
            const double next = dotFluid(res_, pre_);
            if (!(std::abs(rho) > 0)) break;
            const double beta = next / rho;
            rho = next;
            for (std::size_t s = 0; s < n; ++s) dir_[s] = pre_[s] + beta * dir_[s];
        }
    }
    removeMean(sol_);
    // **A projection that was handed a NaN has not converged, and used to say it
    // had.** Both reductions above are `std::max(accumulator, x)`, which drops a
    // NaN rather than carrying it, and every consequence followed from that:
    // `rhsScale` came back 0 on an all-NaN field, `best` came back 0, `tolerance`
    // came back 0, `best > tolerance` was `0 > 0` and **the whole CG solve was
    // skipped**, and then `residual_` reported `0.0` and `capped_` reported
    // `false` over zero iterations. A perfect converged solve, on a velocity field
    // that is not a number.
    //
    // One NaN cell is enough, and in a singular region it is worse: `rhsScale` is
    // taken before `removeMean`, so it is a real positive number, and then the
    // mean subtraction spreads the NaN to every cell. The solve is skipped with a
    // *plausible* tolerance in hand and the pressures come back all-zero and
    // finite -- so the caller sees a healthy report, a clean field, and a
    // divergence that was never removed.
    //
    // `capped_` rather than a new flag, because every consumer already asserts on
    // it -- `test_flip.cpp` at six sites, `test_promotion.cpp` at two,
    // `flip_probe` at one -- and "it ran out of iterations" and "it never started"
    // are both `the answer is not trustworthy`, which is what those call sites are
    // asking. `residual_` is NaN so that any threshold test on it fails too.
    if (nonFinite) {
        residual_ = std::numeric_limits<double>::quiet_NaN();
        capped_ = true;
    } else {
        residual_ = rhsScale > 0 ? best / rhsScale : 0.0;
        capped_ = best > tolerance;
    }

    for (std::size_t s = 0; s < n; ++s)
        pressure_[static_cast<std::size_t>(fluid_[s])] = sol_[s];

    // One pass over every face, so a face shared by two fluid cells is corrected
    // once. Gathering per cell and correcting both of its faces would correct an
    // interior face twice, which halves nothing and doubles the gradient -- and it
    // would still leave a hydrostatic column looking nearly right.
    const double gradient = dt / (params.density * h);
    for (int a = 0; a < 3; ++a)
        for (std::size_t flat = 0; flat < vel_[a].size(); ++flat) {
            if (face_[a][flat] == static_cast<std::uint8_t>(Face::Solid)) {
                vel_[a][flat] = 0.0;
                continue;
            }
            if (face_[a][flat] != static_cast<std::uint8_t>(Face::Fluid)) continue;
            const std::int32_t minus = nb_[flat * 6 + static_cast<std::size_t>(a * 2)];
            const double low = minus < 0 ? 0.0 : pressure_[static_cast<std::size_t>(minus)];
            vel_[a][flat] -= gradient * (pressure_[flat] - low);
        }
}

void Solver::extrapolate(const Params& params) {
    const int depth = std::max(params.extrapolationDepth, 0);
    if (depth == 0) return;
    for (int a = 0; a < 3; ++a) {
        valid_.assign(vel_[a].size(), 0u);
        for (std::size_t flat = 0; flat < vel_[a].size(); ++flat)
            if (face_[a][flat] == static_cast<std::uint8_t>(Face::Fluid)) valid_[flat] = 1u;
        for (int layer = 0; layer < depth; ++layer) {
            scratch_.assign(vel_[a].size(), 0.0);
            filled_.assign(vel_[a].size(), 0u);
            bool any = false;
            for (std::size_t flat = 0; flat < vel_[a].size(); ++flat) {
                if (valid_[flat]) continue;
                if (face_[a][flat] == static_cast<std::uint8_t>(Face::Solid)) continue;
                double total = 0;
                int count = 0;
                for (int d = 0; d < 6; ++d) {
                    const std::int32_t other = nb_[flat * 6 + static_cast<std::size_t>(d)];
                    if (other < 0 || !valid_[static_cast<std::size_t>(other)]) continue;
                    total += vel_[a][static_cast<std::size_t>(other)];
                    ++count;
                }
                if (count == 0) continue;
                scratch_[flat] = total / static_cast<double>(count);
                filled_[flat] = 1u;
                any = true;
            }
            if (!any) break;
            // Two phase, so a layer is a pure function of the layer before it and
            // not of the order the faces happen to be visited in. A single-phase
            // sweep would make the answer depend on the tile insertion order, and
            // §2's bit-identity claim would stop being true.
            for (std::size_t flat = 0; flat < vel_[a].size(); ++flat)
                if (filled_[flat]) {
                    vel_[a][flat] = scratch_[flat];
                    valid_[flat] = 1u;
                }
        }
    }
}

// ---------------------------------------------------------------------------
// Grid -> particle, and advection
// ---------------------------------------------------------------------------

void Solver::sampleVelocity(const double x[3], double out[3]) const {
    for (int a = 0; a < 3; ++a) out[a] = 0.0;
    if (tileKey_.empty() || !(grid_.h > 0)) return;
    const double h = grid_.h;
    const int c[3] = {grid_.cellOf(0, x[0]), grid_.cellOf(1, x[1]), grid_.cellOf(2, x[2])};
    const int t = tileIndex(tileOf(c[0]), tileOf(c[1]), tileOf(c[2]));
    if (t < 0) return;
    TileHalo halo;
    gatherHalo(t, &halo);
    for (int a = 0; a < 3; ++a) {
        double w[3][3], offset[3][3];
        int base[3];
        for (int b = 0; b < 3; ++b) {
            const double s = (x[b] - grid_.lo[b]) / h - (b == a ? 0.0 : 0.5);
            base[b] = splineWeights(s, w[b], offset[b]);
        }
        double total = 0;
        for (int kk = 0; kk < 3; ++kk)
            for (int jj = 0; jj < 3; ++jj)
                for (int ii = 0; ii < 3; ++ii) {
                    const double weight = w[0][ii] * w[1][jj] * w[2][kk];
                    if (weight == 0.0) continue;
                    const int flat = haloIndex(&halo, base[0] + ii, base[1] + jj, base[2] + kk);
                    if (flat < 0) continue;
                    total += weight * vel_[a][static_cast<std::size_t>(flat)];
                }
        out[a] = total;
    }
}

void Solver::transferToParticles(Field& field, const Params& params) const {
    if (tileKey_.empty() || !(grid_.h > 0)) return;
    const double h = grid_.h;
    const double blend = std::min(std::max(params.flipBlend, 0.0), 1.0);
    const int tileCount = static_cast<int>(tileKey_.size() / 3);
    TileHalo halo;
    for (int t = 0; t < tileCount; ++t) {
        const std::size_t from = static_cast<std::size_t>(tileStart_[static_cast<std::size_t>(t)]);
        const std::size_t to = static_cast<std::size_t>(tileStart_[static_cast<std::size_t>(t) + 1]);
        if (from == to) continue;
        gatherHalo(t, &halo);
        for (std::size_t o = from; o < to; ++o) {
            Particle& p = field.particles[static_cast<std::size_t>(order_[o])];
            double picked[3] = {0, 0, 0}, delta[3] = {0, 0, 0}, affine[9] = {0, 0, 0, 0, 0,
                                                                             0, 0, 0, 0};
            for (int a = 0; a < 3; ++a) {
                double w[3][3], offset[3][3];
                int base[3];
                for (int b = 0; b < 3; ++b) {
                    const double s = (p.position[b] - grid_.lo[b]) / h - (b == a ? 0.0 : 0.5);
                    base[b] = splineWeights(s, w[b], offset[b]);
                }
                for (int kk = 0; kk < 3; ++kk)
                    for (int jj = 0; jj < 3; ++jj)
                        for (int ii = 0; ii < 3; ++ii) {
                            const double weight = w[0][ii] * w[1][jj] * w[2][kk];
                            if (weight == 0.0) continue;
                            const int flat =
                                haloIndex(&halo, base[0] + ii, base[1] + jj, base[2] + kk);
                            if (flat < 0) continue;
                            const double now = vel_[a][static_cast<std::size_t>(flat)];
                            picked[a] += weight * now;
                            delta[a] += weight * (now - old_[a][static_cast<std::size_t>(flat)]);
                            if (params.affine) {
                                const double share = weight * now;
                                affine[a * 3 + 0] += share * offset[0][ii];
                                affine[a * 3 + 1] += share * offset[1][jj];
                                affine[a * 3 + 2] += share * offset[2][kk];
                            }
                        }
            }
            for (int a = 0; a < 3; ++a) {
                // `blend` 0 takes the grid velocity outright and 1 takes the
                // particle's own plus the change: PIC and FLIP, and the whole
                // difference between them, on one line.
                p.velocity[a] = blend * (p.velocity[a] + delta[a]) + (1.0 - blend) * picked[a];
            }
            if (params.affine) {
                // `C = B D^-1` with `D = h^2/4 I`, so the inverse is a multiply by
                // `4/h^2` -- and the offsets above are in cells, which absorbs one
                // factor of `h`. See flip.hpp §1 for why only a quadratic kernel
                // has a constant `D` at all.
                for (int e = 0; e < 9; ++e) p.affine[e] = affine[e] * 4.0 / h;
            } else {
                for (int e = 0; e < 9; ++e) p.affine[e] = 0.0;
            }
        }
    }
}

void Solver::advect(Field& field, double dt, const Params& params, Account& account) const {
    if (!(dt > 0)) return;
    const Grid& g = field.grid;
    const double margin = params.wallMargin * g.h;
    for (Particle& p : field.particles) {
        double first[3], mid[3], second[3];
        sampleVelocity(p.position, first);
        for (int a = 0; a < 3; ++a) mid[a] = p.position[a] + 0.5 * dt * first[a];
        sampleVelocity(mid, second);
        for (int a = 0; a < 3; ++a) p.position[a] += dt * second[a];

        // The wall. A particle is put back inside and the component that drove it
        // out is zeroed -- not reflected, because a reflection puts energy back
        // into a body of water that has just been stopped by a bulkhead. Every
        // clamp is counted, so "no mass left the domain" is a claim with a witness.
        bool clamped = false;
        for (int a = 0; a < 3; ++a) {
            const double low = g.lo[a] + margin;
            const double high = g.hi(a) - margin;
            if (p.position[a] < low) {
                p.position[a] = low;
                if (p.velocity[a] < 0) p.velocity[a] = 0;
                clamped = true;
            } else if (p.position[a] > high) {
                p.position[a] = high;
                if (p.velocity[a] > 0) p.velocity[a] = 0;
                clamped = true;
            }
        }
        if (clamped) ++account.clamped;
    }
}

// ---------------------------------------------------------------------------
// The step
// ---------------------------------------------------------------------------

StepResult Solver::step(Field& field, double dt, const Params& params, Account& account) {
    StepResult out;
    out.time = field.time;
    if (!(dt > 0) || field.grid.empty()) {
        out.tiles = tiles();
        out.fluidCells = fluidCells();
        return out;
    }
    if (field.particles.empty()) {
        // No water: no tiles, no bytes, no work. Time still passes, because a
        // caller integrating an empty compartment alongside a full one must not
        // have the two come apart.
        rebuild(field, params);
        field.time += dt;
        out.time = field.time;
        account.mass = field.totalMass();
        account.particles = static_cast<long long>(field.particles.size());
        return out;
    }

    double remaining = dt;
    int budget = std::max(params.maxSubsteps, 1);
    while (remaining > 1e-12 * dt && budget-- > 0) {
        double fastest = 0;
        for (const Particle& p : field.particles)
            for (int a = 0; a < 3; ++a) fastest = std::max(fastest, std::abs(p.velocity[a]));
        double h = params.maxSubstep;
        if (fastest > 0) h = std::min(h, params.courant * field.grid.h / fastest);
        h = std::min(h, remaining);
        if (!(h > 0)) break;

        rebuild(field, params);
        transferToGrid(field, params);
        // Extrapolate *before* the difference is saved as well as after it -- see
        // `saveGrid`. The first pass costs a sweep at `flipBlend = 0`, where the
        // difference is never read; it is taken anyway rather than made conditional,
        // because a branch on a parameter is a code path the default configuration
        // does not test.
        extrapolate(params);
        saveGrid();
        addBodyForce(h, params);
        project(h, params);
        out.projectionIterations = iterations_;
        out.projectionResidual = residual_;
        if (capped_) out.projectionCapped = true;
        out.maxDivergence = maxDivergence();
        extrapolate(params);
        transferToParticles(field, params);
        advect(field, h, params, account);

        field.time += h;
        remaining -= h;
        ++out.substeps;
        out.courant = std::max(out.courant, fastest * h / field.grid.h);
        out.peakSpeed = std::max(out.peakSpeed, fastest);
    }
    // Say so when the budget ran out rather than returning a field that is short of
    // the time it was asked for and looks like any other -- `les::StepResult`'s own
    // reason, and the reason the bound exists at all.
    out.incomplete = remaining > 1e-12 * dt;
    out.time = field.time;
    out.tiles = tiles();
    out.fluidCells = fluidCells();
    out.clamped = account.clamped;

    account.mass = field.totalMass();
    account.particles = static_cast<long long>(field.particles.size());
    return out;
}

}  // namespace sim::flip
