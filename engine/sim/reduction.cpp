// SPDX-License-Identifier: MIT
#include "reduction.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <queue>

namespace sim::reduction {
namespace {

using solidshell::kDof;
using solidshell::kGauss;
using solidshell::kNodes;

double now() {
    using clock = std::chrono::steady_clock;
    return std::chrono::duration<double>(clock::now().time_since_epoch()).count();
}

// --- Dense linear algebra -------------------------------------------------------
//
// Row-major throughout, and small enough that nothing here is blocked or tiled:
// the dense matrices are (boundary + modes) square, which is hundreds, while the
// sparse ones are thousands and are never formed densely.

// Cholesky of a symmetric positive definite `n` x `n` matrix, lower triangle
// written in place over the whole matrix (the upper triangle is zeroed). False on
// a non-positive pivot.
bool cholesky(std::vector<double>& a, int n) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= i; ++j) {
            double s = a[static_cast<std::size_t>(i) * static_cast<std::size_t>(n) +
                         static_cast<std::size_t>(j)];
            for (int k = 0; k < j; ++k)
                s -= a[static_cast<std::size_t>(i) * static_cast<std::size_t>(n) +
                       static_cast<std::size_t>(k)] *
                     a[static_cast<std::size_t>(j) * static_cast<std::size_t>(n) +
                       static_cast<std::size_t>(k)];
            if (i == j) {
                if (!(s > 0.0)) return false;
                a[static_cast<std::size_t>(i) * static_cast<std::size_t>(n) +
                  static_cast<std::size_t>(i)] = std::sqrt(s);
            } else {
                a[static_cast<std::size_t>(i) * static_cast<std::size_t>(n) +
                  static_cast<std::size_t>(j)] =
                    s / a[static_cast<std::size_t>(j) * static_cast<std::size_t>(n) +
                          static_cast<std::size_t>(j)];
            }
        }
        for (int j = i + 1; j < n; ++j)
            a[static_cast<std::size_t>(i) * static_cast<std::size_t>(n) +
              static_cast<std::size_t>(j)] = 0.0;
    }
    return true;
}

// Both triangular solves count in `std::size_t` rather than `int`, and that is not
// only tidiness: with a signed counter GCC cannot bound the trip count, so at `-O3`
// it reasons that `n` might be `INT_MAX`, that `b[i]` would then run past any
// allocation, and warns that "iteration 2147483647 invokes undefined behavior"
// -- `-Waggressive-loop-optimizations`, on the inlined call from
// `generalisedEigen`. An unsigned counter bounded by the same value used to index
// leaves it nothing to prove. It also deletes a wall of casts.
//
// The warning did not fire in any gate, because `verify.sh` builds
// `RelWithDebInfo` and this needs `-O3`. **Warnings are failures here, so a
// warning that only exists in a configuration no gate compiles is a blind spot,
// not a curiosity** -- see the note in `CLAUDE.md`.

// L y = b, in place, with `l` a lower-triangular factor from `cholesky`.
void forwardSolve(const std::vector<double>& l, int n, double* b) {
    if (n <= 0) return;
    const std::size_t size = static_cast<std::size_t>(n);
    for (std::size_t i = 0; i < size; ++i) {
        double s = b[i];
        for (std::size_t k = 0; k < i; ++k) s -= l[i * size + k] * b[k];
        b[i] = s / l[i * size + i];
    }
}

// L^T x = y, in place.
void backwardSolve(const std::vector<double>& l, int n, double* b) {
    if (n <= 0) return;
    const std::size_t size = static_cast<std::size_t>(n);
    for (std::size_t ii = size; ii > 0; --ii) {
        const std::size_t i = ii - 1;
        double s = b[i];
        for (std::size_t k = i + 1; k < size; ++k) s -= l[k * size + i] * b[k];
        b[i] = s / l[i * size + i];
    }
}

// --- The symmetric eigensolver ---------------------------------------------------
//
// Householder tridiagonalisation with the transform accumulated, then implicit QL
// with Wilkinson shifts. Written out rather than taken from a library for the same
// reason this repo has its own PNG codec: it is a bounded, well-specified piece of
// numerical linear algebra with exact properties to test against.

// `a` (n x n, row-major, lower triangle read) is overwritten with the accumulated
// orthogonal transform; `d` takes the diagonal and `e` the sub-diagonal.
void tridiagonalise(std::vector<double>& a, int n, std::vector<double>& d, std::vector<double>& e) {
    const auto A = [&](int i, int j) -> double& {
        return a[static_cast<std::size_t>(i) * static_cast<std::size_t>(n) +
                 static_cast<std::size_t>(j)];
    };
    d.assign(static_cast<std::size_t>(n), 0.0);
    e.assign(static_cast<std::size_t>(n), 0.0);

    for (int i = n - 1; i >= 1; --i) {
        const int l = i - 1;
        double h = 0.0, scale = 0.0;
        if (l > 0) {
            for (int k = 0; k <= l; ++k) scale += std::fabs(A(i, k));
            if (scale == 0.0) {
                e[static_cast<std::size_t>(i)] = A(i, l);
            } else {
                for (int k = 0; k <= l; ++k) {
                    A(i, k) /= scale;
                    h += A(i, k) * A(i, k);
                }
                double f = A(i, l);
                double g = (f >= 0.0) ? -std::sqrt(h) : std::sqrt(h);
                e[static_cast<std::size_t>(i)] = scale * g;
                h -= f * g;
                A(i, l) = f - g;
                f = 0.0;
                for (int j = 0; j <= l; ++j) {
                    A(j, i) = A(i, j) / h;
                    double gg = 0.0;
                    for (int k = 0; k <= j; ++k) gg += A(j, k) * A(i, k);
                    for (int k = j + 1; k <= l; ++k) gg += A(k, j) * A(i, k);
                    e[static_cast<std::size_t>(j)] = gg / h;
                    f += e[static_cast<std::size_t>(j)] * A(i, j);
                }
                const double hh = f / (h + h);
                for (int j = 0; j <= l; ++j) {
                    const double fj = A(i, j);
                    const double gj = e[static_cast<std::size_t>(j)] - hh * fj;
                    e[static_cast<std::size_t>(j)] = gj;
                    for (int k = 0; k <= j; ++k)
                        A(j, k) -= (fj * e[static_cast<std::size_t>(k)] + gj * A(i, k));
                }
            }
        } else {
            e[static_cast<std::size_t>(i)] = A(i, l);
        }
        d[static_cast<std::size_t>(i)] = h;
    }

    d[0] = 0.0;
    e[0] = 0.0;
    for (int i = 0; i < n; ++i) {
        const int l = i - 1;
        if (d[static_cast<std::size_t>(i)] != 0.0) {
            for (int j = 0; j <= l; ++j) {
                double g = 0.0;
                for (int k = 0; k <= l; ++k) g += A(i, k) * A(k, j);
                for (int k = 0; k <= l; ++k) A(k, j) -= g * A(k, i);
            }
        }
        d[static_cast<std::size_t>(i)] = A(i, i);
        A(i, i) = 1.0;
        for (int j = 0; j <= l; ++j) {
            A(j, i) = 0.0;
            A(i, j) = 0.0;
        }
    }
}

// Implicit QL with Wilkinson shifts on the tridiagonal (d, e), accumulating into
// `z` (n x n, row-major; eigenvector j is column j). Returns the sweep count, or
// -1 if any eigenvalue failed to converge in 50 sweeps.
int implicitQl(std::vector<double>& d, std::vector<double>& e, std::vector<double>& z, int n) {
    const auto Z = [&](int i, int j) -> double& {
        return z[static_cast<std::size_t>(i) * static_cast<std::size_t>(n) +
                 static_cast<std::size_t>(j)];
    };
    for (int i = 1; i < n; ++i) e[static_cast<std::size_t>(i - 1)] = e[static_cast<std::size_t>(i)];
    e[static_cast<std::size_t>(n - 1)] = 0.0;

    int sweeps = 0;
    for (int l = 0; l < n; ++l) {
        int iter = 0;
        int m = l;
        do {
            for (m = l; m < n - 1; ++m) {
                const double dd = std::fabs(d[static_cast<std::size_t>(m)]) +
                                  std::fabs(d[static_cast<std::size_t>(m + 1)]);
                if (std::fabs(e[static_cast<std::size_t>(m)]) <=
                    std::numeric_limits<double>::epsilon() * dd)
                    break;
            }
            if (m == l) break;
            if (iter++ == 50) return -1;
            ++sweeps;

            double g = (d[static_cast<std::size_t>(l + 1)] - d[static_cast<std::size_t>(l)]) /
                       (2.0 * e[static_cast<std::size_t>(l)]);
            double r = std::hypot(g, 1.0);
            g = d[static_cast<std::size_t>(m)] - d[static_cast<std::size_t>(l)] +
                e[static_cast<std::size_t>(l)] / (g + (g >= 0.0 ? std::fabs(r) : -std::fabs(r)));
            double s = 1.0, c = 1.0, p = 0.0;
            int i = m - 1;
            for (; i >= l; --i) {
                double f = s * e[static_cast<std::size_t>(i)];
                const double b = c * e[static_cast<std::size_t>(i)];
                r = std::hypot(f, g);
                e[static_cast<std::size_t>(i + 1)] = r;
                if (r == 0.0) {
                    d[static_cast<std::size_t>(i + 1)] -= p;
                    e[static_cast<std::size_t>(m)] = 0.0;
                    break;
                }
                s = f / r;
                c = g / r;
                g = d[static_cast<std::size_t>(i + 1)] - p;
                r = (d[static_cast<std::size_t>(i)] - g) * s + 2.0 * c * b;
                p = s * r;
                d[static_cast<std::size_t>(i + 1)] = g + p;
                g = c * r - b;
                for (int k = 0; k < n; ++k) {
                    f = Z(k, i + 1);
                    Z(k, i + 1) = s * Z(k, i) + c * f;
                    Z(k, i) = c * Z(k, i) - s * f;
                }
            }
            if (r == 0.0 && i >= l) continue;
            d[static_cast<std::size_t>(l)] -= p;
            e[static_cast<std::size_t>(l)] = g;
            e[static_cast<std::size_t>(m)] = 0.0;
        } while (m != l);
    }
    return sweeps;
}

// --- A banded symmetric matrix with an inertia count -----------------------------
//
// `solidshell::BandedSpd` covers the positive definite solve and is reused for it.
// What it cannot do is count eigenvalues, because a Sturm sequence needs the
// factorisation of an **indefinite** shifted matrix: LDL^T rather than Cholesky,
// and the sign of every pivot rather than the first non-positive one.
class BandedLdl {
public:
    BandedLdl(std::size_t n, std::size_t band) : n_(n), b_(band), a_(n * (band + 1), 0.0) {}

    void add(std::size_t row, std::size_t column, double value) {
        if (row < column) return;
        const std::size_t offset = row - column;
        if (offset > b_) return;
        a_[row * (b_ + 1) + offset] += value;
    }

    // Factor and return the number of negative pivots -- by Sylvester's law of
    // inertia, the number of eigenvalues below the shift already folded in. `exact`
    // is false if a pivot came out at zero and had to be perturbed, which makes the
    // count a very strong estimate rather than a proof.
    int negativePivots(bool* exact) {
        if (exact) *exact = true;
        double scale = 0.0;
        for (std::size_t i = 0; i < n_; ++i) scale = std::max(scale, std::fabs(a_[i * (b_ + 1)]));
        const double tiny = scale > 0.0 ? scale * 1e-300 : std::numeric_limits<double>::min();

        int negative = 0;
        for (std::size_t i = 0; i < n_; ++i) {
            const std::size_t lo = i >= b_ ? i - b_ : 0;
            for (std::size_t j = lo; j <= i; ++j) {
                double s = a_[i * (b_ + 1) + (i - j)];
                for (std::size_t k = lo; k < j; ++k)
                    s -= a_[i * (b_ + 1) + (i - k)] * a_[j * (b_ + 1) + (j - k)] *
                         a_[k * (b_ + 1)];
                if (j < i) {
                    a_[i * (b_ + 1) + (i - j)] = s / a_[j * (b_ + 1)];
                } else {
                    if (s == 0.0) {
                        s = tiny;
                        if (exact) *exact = false;
                    }
                    a_[i * (b_ + 1)] = s;
                    if (s < 0.0) ++negative;
                }
            }
        }
        return negative;
    }

private:
    std::size_t n_, b_;
    std::vector<double> a_;  // L below the diagonal, D on it
};

// --- Reverse Cuthill-McKee -------------------------------------------------------
//
// The interior block is solved banded, and banded storage is `n * (band + 1)`, so
// the numbering the caller's mesh happens to arrive in decides how much memory the
// reduction needs and how long its factorisation takes. `makePlateMesh` numbers
// nodes so the band is small; a flood-filled patch or a hold cut out of a ship has
// no such guarantee.
//
// **RCM is not unconditionally better and it was measured not to be.** With a
// minimum-degree start it came out 59 against 41 on the plate and 83 against 53 on
// one ferry patch -- *worse* than the numbering it replaced -- while beating it
// 89 to 173 and 137 to 341 on two others. Two things follow, and both are here:
// the start is the George-Liu **pseudo-peripheral** node rather than the first
// minimum-degree one, and the result is compared against leaving the numbering
// alone and the narrower of the two is kept. Comparing is free -- a bandwidth is a
// pass over the adjacency -- and it makes the ordering incapable of being a
// regression.
//
// The **reversal** is the other half of the name and it is measured to buy exactly
// nothing here: plain Cuthill-McKee gives the identical bandwidth on all five test
// meshes. What reversal reduces is the *profile*, and constant-band storage does
// not exploit a profile. It is kept because it costs one line and a skyline solver
// would want it, not because it was observed to help.

// Breadth-first level structure from `start`: the level of every reachable node,
// and the node in the last level with the smallest degree.
std::uint32_t sweepLevels(const std::vector<std::vector<std::uint32_t>>& adjacency,
                          std::uint32_t start, std::vector<std::int32_t>& level, int& depth) {
    level.assign(adjacency.size(), -1);
    std::queue<std::uint32_t> queue;
    queue.push(start);
    level[start] = 0;
    depth = 0;
    std::uint32_t last = start;
    while (!queue.empty()) {
        const std::uint32_t node = queue.front();
        queue.pop();
        if (level[node] > depth) {
            depth = level[node];
            last = node;
        } else if (level[node] == depth && adjacency[node].size() < adjacency[last].size()) {
            last = node;
        }
        for (std::uint32_t neighbour : adjacency[node])
            if (level[neighbour] < 0) {
                level[neighbour] = level[node] + 1;
                queue.push(neighbour);
            }
    }
    return last;
}

// A node at (or near) the end of a longest path in the component containing
// `seed`. George-Liu: sweep, restart from the last level, stop when the depth
// stops growing.
std::uint32_t pseudoPeripheral(const std::vector<std::vector<std::uint32_t>>& adjacency,
                               std::uint32_t seed) {
    std::vector<std::int32_t> level;
    int depth = 0;
    std::uint32_t node = seed;
    for (int pass = 0; pass < 8; ++pass) {
        int nextDepth = 0;
        const std::uint32_t candidate = sweepLevels(adjacency, node, level, nextDepth);
        if (nextDepth <= depth && pass > 0) break;
        depth = nextDepth;
        node = candidate;
    }
    return node;
}

std::vector<std::uint32_t> reverseCuthillMcKee(const std::vector<std::vector<std::uint32_t>>& adjacency) {
    const std::size_t n = adjacency.size();
    std::vector<std::uint32_t> order;
    order.reserve(n);
    std::vector<std::uint8_t> seen(n, 0);

    while (order.size() < n) {
        std::size_t seed = n;
        std::size_t bestDegree = std::numeric_limits<std::size_t>::max();
        for (std::size_t i = 0; i < n; ++i)
            if (!seen[i] && adjacency[i].size() < bestDegree) {
                bestDegree = adjacency[i].size();
                seed = i;
            }
        if (seed == n) break;
        const std::uint32_t start = pseudoPeripheral(adjacency, static_cast<std::uint32_t>(seed));

        std::queue<std::uint32_t> queue;
        queue.push(start);
        seen[start] = 1;
        while (!queue.empty()) {
            const std::uint32_t node = queue.front();
            queue.pop();
            order.push_back(node);
            std::vector<std::uint32_t> next;
            for (std::uint32_t neighbour : adjacency[node])
                if (!seen[neighbour]) {
                    seen[neighbour] = 1;
                    next.push_back(neighbour);
                }
            std::sort(next.begin(), next.end(), [&](std::uint32_t a, std::uint32_t b) {
                if (adjacency[a].size() != adjacency[b].size())
                    return adjacency[a].size() < adjacency[b].size();
                return a < b;
            });
            for (std::uint32_t neighbour : next) queue.push(neighbour);
        }
    }
    std::reverse(order.begin(), order.end());
    return order;
}

// The half-bandwidth an interior ordering delivers, over the assembled pattern.
// Taken from the pattern rather than estimated, because `BandedSpd::add` silently
// drops anything outside the band it is given.
std::size_t bandwidthOf(const std::vector<std::vector<std::uint32_t>>& adjacency,
                        const std::vector<std::uint32_t>& order) {
    std::vector<std::size_t> position(adjacency.size(), 0);
    for (std::size_t i = 0; i < order.size(); ++i) position[order[i]] = i;
    std::size_t band = 0;
    for (std::size_t i = 0; i < adjacency.size(); ++i)
        for (std::uint32_t neighbour : adjacency[i]) {
            const std::size_t a = position[i], b = position[neighbour];
            band = std::max(band, a > b ? a - b : b - a);
        }
    return band;
}

double vonMises(const double s[6]) {
    const double a = s[0] - s[1], b = s[1] - s[2], c = s[2] - s[0];
    return std::sqrt(0.5 * (a * a + b * b + c * c) +
                     3.0 * (s[3] * s[3] + s[4] * s[4] + s[5] * s[5]));
}

}  // namespace

std::vector<std::uint32_t> bandwidthReducingOrder(
    const std::vector<std::vector<std::uint32_t>>& adjacency) {
    return reverseCuthillMcKee(adjacency);
}

// --- Public dense eigensolvers ---------------------------------------------------

Eigenpairs symmetricEigen(const std::vector<double>& a, int n) {
    Eigenpairs out;
    out.size = n;
    if (n <= 0 || a.size() < static_cast<std::size_t>(n) * static_cast<std::size_t>(n)) {
        out.problem = "matrix is smaller than the size given";
        return out;
    }
    std::vector<double> z(a.begin(),
                          a.begin() + static_cast<std::ptrdiff_t>(static_cast<std::size_t>(n) *
                                                                  static_cast<std::size_t>(n)));
    std::vector<double> d, e;
    tridiagonalise(z, n, d, e);
    const int sweeps = implicitQl(d, e, z, n);
    out.iterations = sweeps < 0 ? 50 : sweeps;
    out.converged = sweeps >= 0;
    if (!out.converged) out.problem = "QL did not converge in 50 sweeps";

    std::vector<int> order(static_cast<std::size_t>(n));
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](int p, int q) {
        return d[static_cast<std::size_t>(p)] < d[static_cast<std::size_t>(q)];
    });

    out.count = n;
    out.value.resize(static_cast<std::size_t>(n));
    out.vector.resize(static_cast<std::size_t>(n) * static_cast<std::size_t>(n));
    for (int j = 0; j < n; ++j) {
        const int src = order[static_cast<std::size_t>(j)];
        out.value[static_cast<std::size_t>(j)] = d[static_cast<std::size_t>(src)];
        // Sign convention: the largest-magnitude component positive, so two runs
        // of the same problem produce the same vectors rather than a sign flip
        // that no assertion could distinguish from a defect.
        int largest = 0;
        for (int k = 1; k < n; ++k)
            if (std::fabs(z[static_cast<std::size_t>(k) * static_cast<std::size_t>(n) +
                            static_cast<std::size_t>(src)]) >
                std::fabs(z[static_cast<std::size_t>(largest) * static_cast<std::size_t>(n) +
                            static_cast<std::size_t>(src)]))
                largest = k;
        const double sign =
            z[static_cast<std::size_t>(largest) * static_cast<std::size_t>(n) +
              static_cast<std::size_t>(src)] < 0.0 ? -1.0 : 1.0;
        for (int k = 0; k < n; ++k)
            out.vector[static_cast<std::size_t>(j) * static_cast<std::size_t>(n) +
                       static_cast<std::size_t>(k)] =
                sign * z[static_cast<std::size_t>(k) * static_cast<std::size_t>(n) +
                         static_cast<std::size_t>(src)];
    }
    return out;
}

Eigenpairs generalisedEigen(const std::vector<double>& a, const std::vector<double>& b, int n) {
    Eigenpairs out;
    out.size = n;
    const std::size_t nn = static_cast<std::size_t>(n) * static_cast<std::size_t>(n);
    if (n <= 0 || a.size() < nn || b.size() < nn) {
        out.problem = "matrix is smaller than the size given";
        return out;
    }
    std::vector<double> l(b.begin(), b.begin() + static_cast<std::ptrdiff_t>(nn));
    if (!cholesky(l, n)) {
        out.problem = "the right-hand matrix is not positive definite";
        return out;
    }

    // C = L^-1 A L^-T, formed as two forward substitution passes: first
    // Y = L^-1 A by columns, then L^-1 Y^T, which is C^T = C because A is
    // symmetric. Nothing is transposed explicitly.
    std::vector<double> c(a.begin(), a.begin() + static_cast<std::ptrdiff_t>(nn));
    std::vector<double> column(static_cast<std::size_t>(n));
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i)
            column[static_cast<std::size_t>(i)] =
                c[static_cast<std::size_t>(i) * static_cast<std::size_t>(n) +
                  static_cast<std::size_t>(j)];
        forwardSolve(l, n, column.data());
        for (int i = 0; i < n; ++i)
            c[static_cast<std::size_t>(i) * static_cast<std::size_t>(n) +
              static_cast<std::size_t>(j)] = column[static_cast<std::size_t>(i)];
    }
    for (int i = 0; i < n; ++i)
        forwardSolve(l, n,
                     c.data() + static_cast<std::size_t>(i) * static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < i; ++j)
            c[static_cast<std::size_t>(j) * static_cast<std::size_t>(n) +
              static_cast<std::size_t>(i)] =
                c[static_cast<std::size_t>(i) * static_cast<std::size_t>(n) +
                  static_cast<std::size_t>(j)];

    out = symmetricEigen(c, n);
    if (!out.converged) return out;
    // x = L^-T y, which restores the B-orthonormality: y^T y = 1 becomes x^T B x = 1.
    for (int j = 0; j < n; ++j)
        backwardSolve(l, n,
                      out.vector.data() + static_cast<std::size_t>(j) * static_cast<std::size_t>(n));
    return out;
}

// --- Choosing an interface --------------------------------------------------------

std::vector<std::uint32_t> nodesNearPlanes(const solidshell::HexMesh& mesh,
                                           const std::vector<Plane>& planes, double tolerance) {
    std::vector<std::uint32_t> out;
    const std::size_t nodes = mesh.nodeCount();
    for (std::size_t n = 0; n < nodes; ++n) {
        const Vec3 p{mesh.position[n * 3], mesh.position[n * 3 + 1], mesh.position[n * 3 + 2]};
        for (const Plane& plane : planes) {
            const Vec3 unit = normalize(plane.normal);
            if (length2(unit) == 0.0) continue;
            if (std::fabs(dot(p - plane.point, unit)) <= tolerance) {
                out.push_back(static_cast<std::uint32_t>(n));
                break;
            }
        }
    }
    return out;
}

std::vector<std::uint32_t> nodesPinned(const solidshell::HexMesh& mesh) {
    std::vector<std::uint32_t> out;
    const std::size_t nodes = mesh.nodeCount();
    for (std::size_t n = 0; n < nodes; ++n)
        for (int axis = 0; axis < 3; ++axis) {
            const std::size_t d = n * 3 + static_cast<std::size_t>(axis);
            if (d < mesh.fixed.size() && mesh.fixed[d]) {
                out.push_back(static_cast<std::uint32_t>(n));
                break;
            }
        }
    return out;
}

// --- The substructure -------------------------------------------------------------

struct Substructure::Impl {
    const solidshell::HexMesh* mesh = nullptr;
    StructuralMaterial material;
    solidshell::Formulation form = solidshell::Formulation::SolidShell;

    std::size_t nodes = 0, dofs = 0;
    std::vector<double> massDiag;  // row-sum lumped, per global DOF

    // The whole stiffness, both triangles, in global DOF numbering. Blocks are
    // extracted from here rather than assembled separately, so there is one
    // assembly and one place a scatter could be wrong.
    std::vector<std::size_t> rowStart;
    std::vector<std::uint32_t> column;
    std::vector<double> value;

    std::vector<std::uint32_t> boundary, interior;  // global DOF, in reduced order
    std::vector<std::int32_t> interiorSlot;         // global DOF -> interior index, or -1

    // What each degree of freedom stands for. Identity unless the `Attachment`
    // constrained it; an eliminated one is in neither partition.
    solidshell::DofExpansion expansion;

    std::size_t band = 0;
    // What the mesher's own numbering would have cost. `band` is the narrower of
    // the two by construction, so this is never the smaller of the pair -- and it
    // is the only way to see whether the renumbering bought anything on a given
    // mesh without hand-instrumenting the ordering.
    // The two candidate orderings, **as node bandwidths over the interior
    // adjacency** -- which is not the same quantity as `band` above and must not be
    // compared with it. `band` is a DOF bandwidth taken from the assembled pattern,
    // roughly three times the node figure but not exactly, because eliminated DOF
    // and attached blocks both move it. Kept so a caller can see which ordering was
    // taken and by how much, without hand-instrumenting the ordering.
    std::size_t naturalNodeBand = 0;
    std::size_t renumberedNodeBand = 0;
    std::unique_ptr<solidshell::BandedSpd> factor;
    bool ready = false;
    double assemblySeconds = 0;
    std::size_t attachedBlocks = 0;
    double attachedMass = 0;
    std::vector<std::string> problems;

    // What each attached member is carrying, as a linear form on the displacement
    // field: `memberForm[m][k]` multiplies `u[memberDof[m][k]]`. The blocks
    // themselves are scattered into the CSR and thrown away -- a rank-one block
    // `s v v^T` has forgotten how to separate the stress from the volume it acts
    // in -- so this is kept rather than reconstructed. See §9.
    std::vector<std::vector<std::uint32_t>> memberDof;
    std::vector<std::vector<double>> memberForm;

    // The slot for (row, col), or `value.size()` when the pattern has no such
    // entry. **The miss is reported rather than assumed away**: `lower_bound` on a
    // column the row does not carry returns the slot of the next one, so an
    // unchecked scatter would add an attached stiffener's term to a *neighbouring*
    // degree of freedom -- a plausible field, silently wrong, which is worse than
    // dropping it -- and on the last row it returns `value.size()`, which is a
    // write past the end. Element scatters can never miss, which is why the
    // unchecked form was safe while elements were all there was; an arbitrary
    // `DofBlock` is why it is not any more. See §8 item 1.
    std::size_t csrSlot(std::uint32_t row, std::uint32_t col) const {
        const auto begin = column.begin() + static_cast<std::ptrdiff_t>(rowStart[row]);
        const auto end = column.begin() + static_cast<std::ptrdiff_t>(rowStart[row + 1]);
        const auto it = std::lower_bound(begin, end, col);
        if (it == end || *it != col) return value.size();
        return static_cast<std::size_t>(it - column.begin());
    }
};

Substructure::Substructure(const solidshell::HexMesh& mesh, const StructuralMaterial& material,
                           const std::vector<std::uint32_t>& interfaceNodes,
                           solidshell::Formulation form)
    : Substructure(mesh, material, interfaceNodes, Attachment{}, form) {}

Substructure::Substructure(const solidshell::HexMesh& mesh, const StructuralMaterial& material,
                           const std::vector<std::uint32_t>& interfaceNodes,
                           const Attachment& attached, solidshell::Formulation form)
    : impl_(std::make_unique<Impl>()) {
    Impl& s = *impl_;
    const double started = now();
    s.mesh = &mesh;
    s.material = material;
    s.form = form;
    s.nodes = mesh.nodeCount();
    s.dofs = s.nodes * 3;

    const std::size_t elements = mesh.elementCount();
    if (elements == 0 || s.nodes == 0) {
        s.problems.push_back("the mesh has no elements");
        return;
    }
    for (std::size_t e = 0; e < elements; ++e) {
        double nodePos[kDof];
        mesh.gather(e, mesh.position, nodePos);
        // `integrable`, not `smallestJacobian > 0`: a collapsed hexahedron -- the
        // wedge a degenerate plate panel extrudes to -- has a zero determinant at
        // the closed edge and a sound one everywhere it is integrated. See
        // `solidshell::ElementShape`.
        if (!solidshell::elementShape(nodePos).integrable) {
            s.problems.push_back("element " + std::to_string(e) + " is inverted or degenerate");
            return;
        }
    }
    // See §3: the mesh's own pins are not consumed. Saying so is the point --
    // a constraint quietly dropped is indistinguishable from one honoured.
    for (std::uint8_t f : mesh.fixed)
        if (f) {
            s.problems.push_back(
                "the mesh pins some degrees of freedom; a substructure is free, so they are "
                "ignored -- constrain the reduced model instead");
            break;
        }

    // --- The attachment, checked before anything is sized from it (§8) ---
    //
    // A block naming a degree of freedom this mesh does not have, or carrying too
    // small a stiffness array, is refused rather than skipped. `solveStatic` skips
    // -- it is a one-shot solve and the caller sees the answer -- but a
    // substructure is built once and asked thousands of times, and the whole point
    // of this section is that a stiffener which quietly does not arrive is
    // indistinguishable from bare plating.
    for (std::size_t b = 0; b < attached.stiffness.size(); ++b) {
        const solidshell::DofBlock& block = attached.stiffness[b];
        const std::size_t n = block.dof.size();
        if (block.stiffness.size() < n * n) {
            s.problems.push_back("attached block " + std::to_string(b) + " has " +
                                 std::to_string(block.stiffness.size()) + " stiffness entries for " +
                                 std::to_string(n) + " degrees of freedom, not " +
                                 std::to_string(n * n));
            return;
        }
        for (std::uint32_t d : block.dof)
            if (d >= s.dofs) {
                s.problems.push_back("attached block " + std::to_string(b) +
                                     " names degree of freedom " + std::to_string(d) +
                                     ", which is not in a mesh of " + std::to_string(s.dofs));
                return;
            }
    }
    if (!attached.mass.empty() && attached.mass.size() != s.nodes) {
        s.problems.push_back("the attached mass is " + std::to_string(attached.mass.size()) +
                             " long; it is one entry per node, so it must be " +
                             std::to_string(s.nodes) + " or empty");
        return;
    }
    if (!attached.stiffness.empty() && attached.mass.empty())
        s.problems.push_back(
            "attached stiffness with no attached mass: the reduced model is stiffer than the "
            "plating but no heavier, so its frequencies come out high -- "
            "`constraint::lumpFiberMass` is what fills this in");

    // The stress forms, checked the same way and for the same reason (§9). A form
    // paired with the wrong block reads the right degrees of freedom for a
    // *different* member and reports a number that is plausible and wrong, which
    // is worse than reporting none -- so a list that is not parallel to the blocks
    // is refused rather than truncated to the shorter of the two.
    if (!attached.stress.empty()) {
        if (attached.stress.size() != attached.stiffness.size()) {
            s.problems.push_back(
                "the attachment has " + std::to_string(attached.stress.size()) +
                " stress forms for " + std::to_string(attached.stiffness.size()) +
                " stiffness blocks; they are parallel, so it must be either count or empty");
            return;
        }
        for (std::size_t b = 0; b < attached.stress.size(); ++b)
            if (attached.stress[b].size() != attached.stiffness[b].dof.size()) {
                s.problems.push_back("attached stress form " + std::to_string(b) + " has " +
                                     std::to_string(attached.stress[b].size()) +
                                     " entries for a block naming " +
                                     std::to_string(attached.stiffness[b].dof.size()) +
                                     " degrees of freedom");
                return;
            }
    } else if (!attached.stiffness.empty()) {
        s.problems.push_back(
            "attached stiffness with no stress forms: `checkValidity` cannot see what these "
            "members are carrying, so a stiffened region is judged by its plating and its "
            "utilisation comes out low -- `constraint::attachedForms` is what fills this in");
    }

    s.attachedBlocks = attached.stiffness.size();
    for (double m : attached.mass) s.attachedMass += m;
    s.memberDof.reserve(attached.stress.size());
    s.memberForm.reserve(attached.stress.size());
    for (std::size_t b = 0; b < attached.stress.size(); ++b) {
        s.memberDof.push_back(attached.stiffness[b].dof);
        s.memberForm.push_back(attached.stress[b]);
    }

    // --- The constraints, and what a node stands for once they are applied ---
    s.expansion = solidshell::DofExpansion(s.dofs, attached.constrained);
    if (!s.expansion.ok()) {
        s.problems.push_back(s.expansion.problem());
        return;
    }
    // Node granularity for the adjacency below: the CSR gives every row of a node
    // the same length, so a node's reach is the union over its three axes.
    std::vector<std::vector<std::uint32_t>> nodeMasters(s.nodes);
    for (std::size_t n = 0; n < s.nodes; ++n) {
        for (int k = 0; k < 3; ++k) {
            const auto d = static_cast<std::uint32_t>(n * 3 + static_cast<std::size_t>(k));
            for (const solidshell::DofExpansion::Term* t = s.expansion.begin(d);
                 t != s.expansion.end(d); ++t)
                nodeMasters[n].push_back(t->dof / 3);
        }
        std::sort(nodeMasters[n].begin(), nodeMasters[n].end());
        nodeMasters[n].erase(std::unique(nodeMasters[n].begin(), nodeMasters[n].end()),
                             nodeMasters[n].end());
    }

    // --- Node adjacency, then the CSR pattern over its degrees of freedom ---
    //
    // The adjacency is over what each node *expands to*, not over the nodes the
    // element names, because that is what the scatter below will touch. An element
    // with a constrained corner couples its other seven nodes to that corner's
    // masters, which are eight nodes on another surface entirely -- so the pattern
    // has to carry a coupling no element edge exists for.
    std::vector<std::vector<std::uint32_t>> nodeAdjacency(s.nodes);
    for (std::size_t e = 0; e < elements; ++e)
        for (int a = 0; a < kNodes; ++a) {
            const std::uint32_t na = mesh.index[e * kNodes + static_cast<std::size_t>(a)];
            for (int b = 0; b < kNodes; ++b) {
                const std::uint32_t nb = mesh.index[e * kNodes + static_cast<std::size_t>(b)];
                for (std::uint32_t p : nodeMasters[na])
                    for (std::uint32_t q : nodeMasters[nb]) nodeAdjacency[p].push_back(q);
            }
        }
    // The attachment joins the adjacency *here*, before the pattern is sized, not
    // as an afterthought once it exists. Everything downstream is derived from
    // this one list -- the CSR, the interior renumbering, and the bandwidth taken
    // from the assembled pattern -- so a block that reaches outside the elements'
    // own reach is carried by all three or by none of them. §8 items 1 and 2.
    //
    // Node granularity, because the CSR gives every row of a node the same length:
    // a block naming one axis of a node gets the whole 3 x 3 sub-block, which is a
    // superset of what it needs and costs a handful of structural zeros.
    for (const solidshell::DofBlock& block : attached.stiffness)
        for (std::uint32_t p : block.dof)
            for (std::uint32_t q : block.dof)
                for (std::uint32_t pm : nodeMasters[p / 3])
                    for (std::uint32_t qm : nodeMasters[q / 3])
                        nodeAdjacency[pm].push_back(qm);
    for (auto& list : nodeAdjacency) {
        std::sort(list.begin(), list.end());
        list.erase(std::unique(list.begin(), list.end()), list.end());
    }

    s.rowStart.assign(s.dofs + 1, 0);
    for (std::size_t n = 0; n < s.nodes; ++n)
        for (int axis = 0; axis < 3; ++axis)
            s.rowStart[n * 3 + static_cast<std::size_t>(axis) + 1] = nodeAdjacency[n].size() * 3;
    for (std::size_t d = 0; d < s.dofs; ++d) s.rowStart[d + 1] += s.rowStart[d];
    s.column.resize(s.rowStart[s.dofs]);
    s.value.assign(s.rowStart[s.dofs], 0.0);
    for (std::size_t n = 0; n < s.nodes; ++n)
        for (int axis = 0; axis < 3; ++axis) {
            std::size_t at = s.rowStart[n * 3 + static_cast<std::size_t>(axis)];
            for (std::uint32_t neighbour : nodeAdjacency[n])
                for (int k = 0; k < 3; ++k)
                    s.column[at++] = neighbour * 3 + static_cast<std::uint32_t>(k);
        }

    // --- Element stiffness and lumped mass ---
    //
    // Every scatter below goes through `csrSlot`, which returns `value.size()` on a
    // pattern miss. The element scatters cannot miss -- the pattern is built from
    // their own adjacency -- and are checked anyway, because the counter is what
    // makes the attachment's "the pattern covers it by construction" a measurement
    // instead of a claim: delete the adjacency fold above and this fires.
    std::size_t patternMisses = 0;
    s.massDiag.assign(s.dofs, 0.0);
    // One scatter for elements and blocks alike, expanded through the constraints,
    // so the matrix assembled **is** `T^T K T` rather than `K` with a correction
    // bolted on beside it. With nothing constrained every expansion is one term of
    // weight one and this is the loop it always was.
    const auto scatter = [&](std::uint32_t row, std::uint32_t column, double value) {
        for (const solidshell::DofExpansion::Term* p = s.expansion.begin(row);
             p != s.expansion.end(row); ++p)
            for (const solidshell::DofExpansion::Term* q = s.expansion.begin(column);
                 q != s.expansion.end(column); ++q) {
                const std::size_t slot = s.csrSlot(p->dof, q->dof);
                if (slot >= s.value.size()) {
                    ++patternMisses;
                    continue;
                }
                s.value[slot] += p->weight * q->weight * value;
            }
    };
    // The mass of an eliminated degree of freedom goes to its masters by the same
    // weights. `T^T M T` is not diagonal, and its row sums for one slave are
    // `weight[a] * m`, which sum to `m` exactly when the weights do -- so total mass
    // is preserved for any constraint and the lumping is the row sum, as everywhere
    // else here. **Unlike `constraint::lumpFiberMass` there is usually nothing to
    // give up:** an eccentric fibre extrapolates *far* -- a weight of 8.83 for a
    // 200 mm bar on 12 mm plating -- so its consistent row sums are `-7.83 m` and
    // `+8.83 m` and the mass has to be split equally instead, giving up the first
    // moment. A junction tie interpolates a point in a face, so its weights are
    // near a convex combination and the row sums keep the first moment as well as
    // the total. **Near, not always**: a face overshoot of `d` puts `-d/2` on a
    // master, and the reference ferry's 9 mm junction gap makes a through-thickness
    // weight 1.69 rather than 0.5. So this is not an argument that the diagonal
    // stays positive -- the check below is, and it is why it is a check.
    const auto lump = [&](std::uint32_t d, double m) {
        for (const solidshell::DofExpansion::Term* t = s.expansion.begin(d);
             t != s.expansion.end(d); ++t)
            s.massDiag[t->dof] += t->weight * m;
    };
    std::vector<double> ke(static_cast<std::size_t>(kDof) * kDof);
    for (std::size_t e = 0; e < elements; ++e) {
        double nodePos[kDof];
        mesh.gather(e, mesh.position, nodePos);
        solidshell::elementStiffness(nodePos, material, form, ke.data());

        std::uint32_t dof[kDof];
        for (int a = 0; a < kNodes; ++a) {
            const std::uint32_t n = mesh.index[e * kNodes + static_cast<std::size_t>(a)];
            for (int k = 0; k < 3; ++k)
                dof[a * 3 + k] = n * 3 + static_cast<std::uint32_t>(k);
        }
        for (int p = 0; p < kDof; ++p)
            for (int q = 0; q < kDof; ++q)
                scatter(dof[p], dof[q],
                        ke[static_cast<std::size_t>(p) * kDof + static_cast<std::size_t>(q)]);

        double lumped[kNodes];
        solidshell::elementMass(nodePos, material.density, lumped);
        for (int a = 0; a < kNodes; ++a) {
            const std::uint32_t n = mesh.index[e * kNodes + static_cast<std::size_t>(a)];
            for (int k = 0; k < 3; ++k)
                lump(n * 3 + static_cast<std::uint32_t>(k), lumped[a]);
        }
    }

    // --- The attachment's stiffness, into the same CSR the elements went into ---
    for (const solidshell::DofBlock& block : attached.stiffness) {
        const std::size_t n = block.dof.size();
        for (std::size_t p = 0; p < n; ++p)
            for (std::size_t q = 0; q < n; ++q)
                scatter(block.dof[p], block.dof[q], block.stiffness[p * n + q]);
    }
    // Like the band check further down, this cannot fire while the adjacency fold
    // above is right, and mutation testing confirms no single edit reaches it. It
    // is the executable form of the invariant that fold exists to establish, and
    // the alternative to reporting a miss is a stiffener landing on a neighbouring
    // degree of freedom, which is a field nothing downstream could tell from an
    // answer.
    if (patternMisses != 0) {
        s.problems.push_back(
            std::to_string(patternMisses) +
            " stiffness entries had no slot in the sparsity pattern and were dropped; the "
            "pattern does not cover what was assembled into it");
        return;
    }

    // --- And its mass, onto the same lumped diagonal (§8 item 3) ---
    // One entry per node, spread over that node's three translations, exactly as
    // the element mass above is: an isotropic point mass, which is what a lumped
    // mass matrix means.
    for (std::size_t n = 0; n < attached.mass.size(); ++n)
        for (int k = 0; k < 3; ++k)
            lump(static_cast<std::uint32_t>(n * 3 + static_cast<std::size_t>(k)), attached.mass[n]);

    for (std::size_t d = 0; d < s.dofs; ++d) {
        // An eliminated degree of freedom carries no mass because it has no unknown;
        // its steel has gone to its masters above. Everything else must have some,
        // and a *negative* diagonal is the failure a constraint with a negative
        // weight would produce -- an extrapolating tie -- so the test is `> 0` on
        // the live degrees of freedom rather than `!= 0`.
        if (s.expansion.eliminated(static_cast<std::uint32_t>(d))) continue;
        if (!(s.massDiag[d] > 0.0)) {
            s.problems.push_back("node " + std::to_string(d / 3) +
                                 " carries no lumped mass: the element around it is degenerate");
            return;
        }
    }

    // --- Partition ---
    std::vector<std::uint8_t> isInterface(s.nodes, 0);
    for (std::uint32_t n : interfaceNodes) {
        if (n >= s.nodes) {
            s.problems.push_back("interface node " + std::to_string(n) + " is not in the mesh");
            return;
        }
        isInterface[n] = 1;
    }
    // The interface is what a reduction keeps exactly, so a boundary degree of
    // freedom that is secretly a weighted sum of others is not kept -- and a caller
    // who tied one would get a reduced model whose boundary silently has fewer
    // unknowns than it names. Refused rather than dropped.
    for (std::uint32_t n : interfaceNodes)
        for (int k = 0; k < 3; ++k)
            if (s.expansion.eliminated(n * 3 + static_cast<std::uint32_t>(k))) {
                s.problems.push_back("interface node " + std::to_string(n) +
                                     " has a degree of freedom eliminated by a constraint; the"
                                     " interface is kept exactly and cannot be a function of"
                                     " something else");
                return;
            }
    // A node every one of whose degrees of freedom is eliminated is in neither
    // partition: it is not an unknown at all. Leaving it in the interior would put
    // an empty row in K_ii and the factorisation would find a zero pivot.
    const auto live = [&](std::size_t n) {
        for (int k = 0; k < 3; ++k)
            if (!s.expansion.eliminated(static_cast<std::uint32_t>(n * 3 + static_cast<std::size_t>(k))))
                return true;
        return false;
    };
    std::vector<std::uint32_t> interiorNodes;
    for (std::size_t n = 0; n < s.nodes; ++n)
        if (!isInterface[n] && live(n)) interiorNodes.push_back(static_cast<std::uint32_t>(n));
    for (std::size_t n = 0; n < s.nodes; ++n)
        if (isInterface[n])
            for (int k = 0; k < 3; ++k)
                s.boundary.push_back(static_cast<std::uint32_t>(n * 3 + static_cast<std::size_t>(k)));
    if (s.boundary.empty()) {
        s.problems.push_back("the interface is empty: there is nothing to keep exactly");
        return;
    }
    if (interiorNodes.empty()) {
        s.problems.push_back("the interface is the whole mesh: there is nothing to reduce");
        return;
    }

    // K_ii is the stiffness with the interface held, so it is singular unless the
    // interface removes all six rigid body modes -- three non-collinear nodes.
    //
    // **The factorisation does not catch this and it was measured not to.** A
    // mechanism leaves an exactly zero pivot in exact arithmetic and a *tiny
    // positive* one in floating point, so `BandedSpd::factor` returns true and the
    // solve returns garbage that nothing downstream can tell from an answer. The
    // precondition is geometric, so it is checked geometrically, before anything
    // is factored. (It is necessary and, for a well-formed hex mesh, sufficient: a
    // mechanism *inside* the mesh would survive it, and a caller who suspects one
    // can count `eigenvaluesBelow` a small positive shift.)
    {
        // Collinearity is the whole test, and a node count is not a second one: two
        // nodes are collinear, one is, and none is. Testing the count separately
        // reads as belt and braces and is in fact dead -- which mutation testing
        // said, by finding that changing the threshold from three to two altered
        // nothing at all. The count survives only in the message, where it does tell
        // a caller something the word "collinear" does not.
        std::size_t interfaceNodeCount = 0;
        Vec3 first{}, axis{};
        bool haveFirst = false, haveAxis = false, collinear = true;
        double extent = 0;
        for (std::size_t n = 0; n < s.nodes; ++n) {
            if (!isInterface[n]) continue;
            ++interfaceNodeCount;
            const Vec3 p{mesh.position[n * 3], mesh.position[n * 3 + 1], mesh.position[n * 3 + 2]};
            if (!haveFirst) {
                first = p;
                haveFirst = true;
                continue;
            }
            extent = std::max(extent, length(p - first));
            if (!haveAxis) {
                if (length(p - first) > 0.0) {
                    axis = normalize(p - first);
                    haveAxis = true;
                }
                continue;
            }
            if (length(cross(p - first, axis)) > 1e-9 * std::max(extent, 1.0)) collinear = false;
        }
        if (collinear) {
            s.problems.push_back(
                "the interface is " +
                std::string(interfaceNodeCount < 3 ? "fewer than three nodes" : "collinear") +
                ", so it does not restrain all six rigid body modes and the interior stiffness is "
                "singular");
            return;
        }
    }

    // Interior numbering by RCM, so the band is a property of the geometry rather
    // than of the order the caller's mesher happened to emit nodes in.
    std::vector<std::int32_t> localOf(s.nodes, -1);
    for (std::size_t i = 0; i < interiorNodes.size(); ++i)
        localOf[interiorNodes[i]] = static_cast<std::int32_t>(i);
    std::vector<std::vector<std::uint32_t>> interiorAdjacency(interiorNodes.size());
    for (std::size_t i = 0; i < interiorNodes.size(); ++i)
        for (std::uint32_t neighbour : nodeAdjacency[interiorNodes[i]]) {
            const std::int32_t local = localOf[neighbour];
            if (local >= 0 && static_cast<std::size_t>(local) != i)
                interiorAdjacency[i].push_back(static_cast<std::uint32_t>(local));
        }
    // Whichever of the two is narrower -- see the note above `reverseCuthillMcKee`:
    // it is not unconditionally better than the numbering the mesher chose, and the
    // comparison costs one pass over the adjacency.
    std::vector<std::uint32_t> natural(interiorNodes.size());
    std::iota(natural.begin(), natural.end(), std::uint32_t{0});
    std::vector<std::uint32_t> renumbered = reverseCuthillMcKee(interiorAdjacency);
    // An ordering that does not cover every interior node is a defect, not a
    // choice, and taking the narrower of two would quietly paper over it -- an
    // interior in two disconnected pieces is exactly what a substructure cut at
    // midspan has, and a breadth-first sweep that stops after the first component
    // would deliver half a model with nothing to say so.
    if (renumbered.size() != interiorNodes.size()) {
        s.problems.push_back("the interior renumbering covered " +
                             std::to_string(renumbered.size()) + " of " +
                             std::to_string(interiorNodes.size()) +
                             " nodes; the mesher's own numbering is used instead");
        renumbered = natural;
    }
    const std::size_t naturalBand = bandwidthOf(interiorAdjacency, natural);
    const std::size_t renumberedBand = bandwidthOf(interiorAdjacency, renumbered);
    s.naturalNodeBand = naturalBand;
    s.renumberedNodeBand = renumberedBand;
    const std::vector<std::uint32_t>& order =
        renumberedBand <= naturalBand ? renumbered : natural;

    s.interiorSlot.assign(s.dofs, -1);
    s.interior.reserve(interiorNodes.size() * 3);
    for (std::uint32_t local : order) {
        const std::uint32_t node = interiorNodes[local];
        for (int k = 0; k < 3; ++k) {
            const std::uint32_t d = node * 3 + static_cast<std::uint32_t>(k);
            if (s.expansion.eliminated(d)) continue;
            s.interiorSlot[d] = static_cast<std::int32_t>(s.interior.size());
            s.interior.push_back(d);
        }
    }

    // --- K_ii, banded. The bandwidth is taken from the assembled pattern rather
    // than estimated from the elements: `BandedSpd::add` silently drops anything
    // outside the band it was given, so an estimate that came out one too small
    // would produce a plausible wrong answer with nothing to say so.
    const std::size_t ni = s.interior.size();
    for (std::size_t p = 0; p < ni; ++p) {
        const std::uint32_t row = s.interior[p];
        for (std::size_t k = s.rowStart[row]; k < s.rowStart[row + 1]; ++k) {
            const std::int32_t q = s.interiorSlot[s.column[k]];
            if (q < 0) continue;
            const std::size_t qq = static_cast<std::size_t>(q);
            s.band = std::max(s.band, p > qq ? p - qq : qq - p);
        }
    }
    s.factor = std::make_unique<solidshell::BandedSpd>(ni, s.band);
    // **`BandedSpd::add` drops an entry outside the band without a word**, and an
    // attached block can tie degrees of freedom no element shares -- which is
    // precisely how a band computed from the elements alone ends up one short and
    // the interior solve ends up quietly back at the bare plating. The band above
    // is taken from the assembled pattern, which carries the blocks, so it covers
    // them; this counts what was offered against what the band admits, so that
    // "covers them" is measured rather than asserted.
    //
    // **Which means it cannot fire while the band above is right, and mutation
    // testing says so rather than leaving it to be assumed**: no single edit to
    // this file reaches it. Two together do -- a band capped at what the elements
    // imply *and* this check deleted -- and the tests kill that combination as
    // well, so what the check buys is a stated reason in `problems()` in place of
    // what that mutant actually delivered, which was a segmentation fault three
    // tests later. It is kept as the executable form of the invariant the band is
    // chosen to satisfy, and it costs one comparison per stored entry.
    std::size_t outsideBand = 0;
    for (std::size_t p = 0; p < ni; ++p) {
        const std::uint32_t row = s.interior[p];
        for (std::size_t k = s.rowStart[row]; k < s.rowStart[row + 1]; ++k) {
            const std::int32_t q = s.interiorSlot[s.column[k]];
            if (q < 0) continue;
            const std::size_t qq = static_cast<std::size_t>(q);
            if ((p > qq ? p - qq : qq - p) > s.band) {
                ++outsideBand;
                continue;
            }
            s.factor->add(p, qq, s.value[k]);
        }
    }
    if (outsideBand != 0) {
        s.problems.push_back(std::to_string(outsideBand) +
                             " interior stiffness entries fell outside the half bandwidth of " +
                             std::to_string(s.band) +
                             " and were dropped; the band does not cover the assembly");
        s.factor.reset();
        return;
    }
    if (!s.factor->factor()) {
        s.problems.push_back(
            "the interior stiffness is not positive definite: the interface does not restrain "
            "all six rigid body modes (it needs at least three non-collinear nodes)");
        s.factor.reset();
        return;
    }

    s.ready = true;
    s.assemblySeconds = now() - started;
}

Substructure::~Substructure() = default;
Substructure::Substructure(Substructure&&) noexcept = default;
Substructure& Substructure::operator=(Substructure&&) noexcept = default;

bool Substructure::ready() const { return impl_->ready; }
const std::vector<std::string>& Substructure::problems() const { return impl_->problems; }
std::size_t Substructure::nodeCount() const { return impl_->nodes; }
std::size_t Substructure::dofCount() const { return impl_->dofs; }
std::size_t Substructure::boundaryCount() const { return impl_->boundary.size(); }
std::size_t Substructure::interiorCount() const { return impl_->interior.size(); }
std::size_t Substructure::halfBandwidth() const { return impl_->band; }
std::size_t Substructure::naturalNodeBandwidth() const { return impl_->naturalNodeBand; }
std::size_t Substructure::renumberedNodeBandwidth() const { return impl_->renumberedNodeBand; }
double Substructure::assemblySeconds() const { return impl_->assemblySeconds; }
std::size_t Substructure::attachedBlocks() const { return impl_->attachedBlocks; }
double Substructure::attachedMass() const { return impl_->attachedMass; }
const std::vector<std::uint32_t>& Substructure::boundaryDof() const { return impl_->boundary; }
const std::vector<std::uint32_t>& Substructure::interiorDof() const { return impl_->interior; }
const std::vector<double>& Substructure::mass() const { return impl_->massDiag; }
const solidshell::HexMesh& Substructure::mesh() const { return *impl_->mesh; }
const StructuralMaterial& Substructure::material() const { return impl_->material; }
solidshell::Formulation Substructure::formulation() const { return impl_->form; }
const solidshell::DofExpansion& Substructure::expansion() const { return impl_->expansion; }

std::size_t Substructure::attachedMembers() const { return impl_->memberForm.size(); }

double Substructure::memberStress(std::size_t member,
                                  const std::vector<double>& displacement) const {
    if (member >= impl_->memberForm.size()) return 0.0;
    const std::vector<double>& form = impl_->memberForm[member];
    const std::vector<std::uint32_t>& dof = impl_->memberDof[member];
    double stress = 0;
    for (std::size_t k = 0; k < form.size() && k < dof.size(); ++k)
        if (dof[k] < displacement.size()) stress += form[k] * displacement[dof[k]];
    return stress;
}

// **The two accessors below are bounded by what was actually built, not by
// `nodes`, and that is a fix rather than a nicety.** Every early `return` in the
// constructor -- an inverted element, a malformed `Attachment`, a constraint that
// does not compose -- leaves `nodes` and `dofs` set from the mesh while `massDiag`
// and the CSR are still empty, so a caller that reads a refused substructure
// instead of checking `ready()` used to run off the front of an empty vector.
// Mutation testing found it, as a **segmentation fault three tests after** the
// mutant it was chasing, and it was reachable long before any constraint existed.
// A caller who ignores `ready()` now gets a zero, which is wrong and says so, in
// place of a read of address zero.
double Substructure::totalMass() const {
    double total = 0.0;
    for (std::size_t n = 0; n * 3 < impl_->massDiag.size(); ++n) total += impl_->massDiag[n * 3];
    return total;
}

void Substructure::stiffnessTimes(const std::vector<double>& x, std::vector<double>& y) const {
    const Impl& s = *impl_;
    y.assign(s.dofs, 0.0);
    if (s.rowStart.size() != s.dofs + 1) return;
    for (std::size_t row = 0; row < s.dofs; ++row) {
        double sum = 0.0;
        for (std::size_t k = s.rowStart[row]; k < s.rowStart[row + 1]; ++k)
            sum += s.value[k] * x[s.column[k]];
        y[row] = sum;
    }
}

void Substructure::boundaryRow(std::size_t b, std::vector<double>& toBoundary,
                               std::vector<double>& toInterior) const {
    const Impl& s = *impl_;
    toBoundary.assign(s.boundary.size(), 0.0);
    toInterior.assign(s.interior.size(), 0.0);
    if (b >= s.boundary.size()) return;
    // Boundary DOF are ordered by global index, so the position of a boundary
    // column within the boundary partition is a lower-bound search on that list.
    const std::uint32_t row = s.boundary[b];
    for (std::size_t k = s.rowStart[row]; k < s.rowStart[row + 1]; ++k) {
        const std::uint32_t col = s.column[k];
        const std::int32_t interior = s.interiorSlot[col];
        if (interior >= 0) {
            toInterior[static_cast<std::size_t>(interior)] = s.value[k];
        } else {
            const auto it = std::lower_bound(s.boundary.begin(), s.boundary.end(), col);
            if (it != s.boundary.end() && *it == col)
                toBoundary[static_cast<std::size_t>(it - s.boundary.begin())] = s.value[k];
        }
    }
}

bool Substructure::interiorSolve(std::vector<double>& r) const {
    if (!impl_->factor) return false;
    impl_->factor->solve(r);
    return true;
}

void Substructure::interiorStiffnessTimes(const std::vector<double>& x,
                                          std::vector<double>& y) const {
    const Impl& s = *impl_;
    const std::size_t ni = s.interior.size();
    y.assign(ni, 0.0);
    for (std::size_t p = 0; p < ni; ++p) {
        const std::uint32_t row = s.interior[p];
        double sum = 0.0;
        for (std::size_t k = s.rowStart[row]; k < s.rowStart[row + 1]; ++k) {
            const std::int32_t q = s.interiorSlot[s.column[k]];
            if (q >= 0) sum += s.value[k] * x[static_cast<std::size_t>(q)];
        }
        y[p] = sum;
    }
}

int Substructure::eigenvaluesBelow(double shift, bool* exact) const {
    const Impl& s = *impl_;
    if (!s.ready) {
        if (exact) *exact = false;
        return 0;
    }
    const std::size_t ni = s.interior.size();
    BandedLdl shifted(ni, s.band);
    for (std::size_t p = 0; p < ni; ++p) {
        const std::uint32_t row = s.interior[p];
        for (std::size_t k = s.rowStart[row]; k < s.rowStart[row + 1]; ++k) {
            const std::int32_t q = s.interiorSlot[s.column[k]];
            if (q >= 0) shifted.add(p, static_cast<std::size_t>(q), s.value[k]);
        }
        shifted.add(p, p, -shift * s.massDiag[row]);
    }
    return shifted.negativePivots(exact);
}

Eigenpairs Substructure::fixedInterfaceModes(int count, double tolerance, int maxIterations) const {
    const Impl& s = *impl_;
    Eigenpairs out;
    const std::size_t ni = s.interior.size();
    out.size = static_cast<int>(ni);
    if (!s.ready || count <= 0 || ni == 0) {
        out.converged = count <= 0;
        return out;
    }
    if (static_cast<std::size_t>(count) > ni) {
        out.problem = "asked for more modes than the interior has degrees of freedom";
        count = static_cast<int>(ni);
    }

    // A block of `q` vectors rather than `count`, because subspace iteration
    // converges at the rate (lambda_j / lambda_q)^2 -- the extra vectors are what
    // separate the modes that are wanted from the ones just above them.
    const int q = static_cast<int>(std::min<std::size_t>(ni, static_cast<std::size_t>(
                                       std::max(count + 8, 2 * count))));
    const std::size_t nq = static_cast<std::size_t>(q);

    std::vector<double> mi(ni);
    for (std::size_t p = 0; p < ni; ++p) mi[p] = s.massDiag[s.interior[p]];
    std::vector<double> kd(ni);
    for (std::size_t p = 0; p < ni; ++p) {
        const std::uint32_t row = s.interior[p];
        kd[p] = 1.0;
        for (std::size_t k = s.rowStart[row]; k < s.rowStart[row + 1]; ++k)
            if (s.column[k] == row) kd[p] = s.value[k];
    }

    // Bathe's starting block: the mass itself, then unit vectors at the degrees of
    // freedom with the largest mass-to-stiffness ratio -- the ones the lowest modes
    // move most -- and one spread vector so the block cannot start rank deficient.
    //
    // The spread vector is a deterministic hash rather than alternating signs, and
    // rather than `waves.cpp`'s Threefry -- which is the repo's generator but is
    // private to that file, and what is wanted here is not a random number but a
    // vector no other column happens to lie in. Alternating signs is what this had
    // first, and it **fails**: with q close to the interior size the unit vectors
    // cover all but a couple of coordinates, and two vectors whose entries in the
    // remaining two happen to match leave the block singular, which shows up as a
    // projected mass matrix that will not factor. Deterministic on purpose -- a
    // random start would make the iteration count irreproducible, and this file's
    // whole claim is that its answers do not move.
    //
    // When the block is as wide as the interior there is nothing to iterate towards
    // and the identity is both exactly non-singular and exactly right: one pass
    // then solves the whole eigenproblem in the projected space.
    std::vector<double> x(ni * nq, 0.0);
    if (static_cast<std::size_t>(q) == ni) {
        for (std::size_t p = 0; p < ni; ++p) x[p * nq + p] = 1.0;
    } else {
        for (std::size_t p = 0; p < ni; ++p) x[p] = mi[p];
        std::vector<std::size_t> byRatio(ni);
        std::iota(byRatio.begin(), byRatio.end(), std::size_t{0});
        std::sort(byRatio.begin(), byRatio.end(), [&](std::size_t a, std::size_t b) {
            const double ra = mi[a] / kd[a], rb = mi[b] / kd[b];
            if (ra != rb) return ra > rb;
            return a < b;
        });
        for (int j = 1; j < q - 1; ++j)
            x[byRatio[static_cast<std::size_t>(j - 1) % ni] * nq + static_cast<std::size_t>(j)] =
                1.0;
        for (std::size_t p = 0; p < ni; ++p) {
            std::uint64_t h = p * 0x9e3779b97f4a7c15ull + 0x632be59bd9b4e019ull;
            h ^= h >> 30;
            h *= 0xbf58476d1ce4e5b9ull;
            h ^= h >> 27;
            x[p * nq + nq - 1] =
                2.0 * (static_cast<double>(h >> 11) / static_cast<double>(1ull << 53)) - 1.0;
        }
    }

    std::vector<double> previous(static_cast<std::size_t>(count),
                                 std::numeric_limits<double>::max());
    std::vector<double> rhs(ni * nq), bar(ni * nq), small(nq * nq);
    std::vector<double> work(std::max(ni, nq)), column(ni), product(ni), dot(nq);

    for (int iteration = 1; iteration <= maxIterations; ++iteration) {
        for (std::size_t p = 0; p < ni; ++p)
            for (std::size_t j = 0; j < nq; ++j) rhs[p * nq + j] = mi[p] * x[p * nq + j];

        for (std::size_t j = 0; j < nq; ++j) {
            for (std::size_t p = 0; p < ni; ++p) work[p] = rhs[p * nq + j];
            s.factor->solve(work);
            for (std::size_t p = 0; p < ni; ++p) bar[p * nq + j] = work[p];
        }

        // **The block is re-orthogonalised in the mass inner product before it is
        // projected, and that is not a refinement.** Without it, repeated
        // multiplication by K^-1 M drives every column towards the lowest modes and
        // the block goes numerically rank deficient -- measured, on a 210-DOF
        // interior with a block of 200: the projected mass matrix stopped being
        // positive definite, `generalisedEigen` refused it, and the whole reduction
        // fell back to Guyan with only a line in `problems` to say so. Two passes of
        // classical Gram-Schmidt, which is as stable as the modified form and reads
        // down the rows rather than across them; a column whose norm collapses is
        // replaced by a fresh deterministic vector rather than kept as noise.
        for (std::size_t j = 0; j < nq; ++j) {
            for (std::size_t p = 0; p < ni; ++p) column[p] = bar[p * nq + j];
            double before = 0.0;
            for (std::size_t p = 0; p < ni; ++p) before += mi[p] * column[p] * column[p];
            for (int attempt = 0; attempt < 2; ++attempt) {
                for (int pass = 0; pass < 2; ++pass) {
                    std::fill(dot.begin(), dot.begin() + static_cast<std::ptrdiff_t>(j), 0.0);
                    for (std::size_t p = 0; p < ni; ++p) {
                        const double mv = mi[p] * column[p];
                        const double* row = bar.data() + p * nq;
                        for (std::size_t k = 0; k < j; ++k) dot[k] += mv * row[k];
                    }
                    for (std::size_t p = 0; p < ni; ++p) {
                        const double* row = bar.data() + p * nq;
                        double sum = 0.0;
                        for (std::size_t k = 0; k < j; ++k) sum += dot[k] * row[k];
                        column[p] -= sum;
                    }
                }
                double norm = 0.0;
                for (std::size_t p = 0; p < ni; ++p) norm += mi[p] * column[p] * column[p];
                if (norm > 1e-24 * before || attempt == 1) {
                    const double scale = norm > 0.0 ? 1.0 / std::sqrt(norm) : 0.0;
                    for (std::size_t p = 0; p < ni; ++p) column[p] *= scale;
                    break;
                }
                for (std::size_t p = 0; p < ni; ++p) {
                    std::uint64_t h = (p * nq + j) * 0x9e3779b97f4a7c15ull + 0x94d049bb133111ebull;
                    h ^= h >> 29;
                    h *= 0xbf58476d1ce4e5b9ull;
                    h ^= h >> 32;
                    column[p] =
                        2.0 * (static_cast<double>(h >> 11) / static_cast<double>(1ull << 53)) - 1.0;
                }
                before = 0.0;
                for (std::size_t p = 0; p < ni; ++p) before += mi[p] * column[p] * column[p];
            }
            for (std::size_t p = 0; p < ni; ++p) bar[p * nq + j] = column[p];
        }

        // With the block M-orthonormal the projected mass matrix is the identity by
        // construction, so the projected problem is a standard symmetric one and
        // there is no Cholesky left to fail. K* = Xbar^T K Xbar has to be formed
        // against K rather than reused from `M X`, because the orthogonalisation
        // above broke the identity `K Xbar = M X` that shortcut relied on.
        std::fill(small.begin(), small.end(), 0.0);
        for (std::size_t j = 0; j < nq; ++j) {
            for (std::size_t p = 0; p < ni; ++p) column[p] = bar[p * nq + j];
            interiorStiffnessTimes(column, product);
            for (std::size_t i = 0; i <= j; ++i) {
                double sum = 0.0;
                for (std::size_t p = 0; p < ni; ++p) sum += bar[p * nq + i] * product[p];
                small[i * nq + j] = sum;
                small[j * nq + i] = sum;
            }
        }

        const Eigenpairs ritz = symmetricEigen(small, q);
        if (!ritz.converged) {
            out.problem = "the projected eigenproblem failed: " + ritz.problem;
            return out;
        }

        for (std::size_t p = 0; p < ni; ++p) {
            for (std::size_t j = 0; j < nq; ++j) {
                double sum = 0.0;
                for (std::size_t k = 0; k < nq; ++k)
                    sum += bar[p * nq + k] * ritz.vector[j * nq + k];
                work[j] = sum;
            }
            for (std::size_t j = 0; j < nq; ++j) x[p * nq + j] = work[j];
        }

        double worst = 0.0;
        for (int j = 0; j < count; ++j) {
            const double value = ritz.value[static_cast<std::size_t>(j)];
            const double reference = std::max(std::fabs(value), 1e-300);
            worst = std::max(worst,
                             std::fabs(value - previous[static_cast<std::size_t>(j)]) / reference);
            previous[static_cast<std::size_t>(j)] = value;
        }
        out.iterations = iteration;
        out.lastChange = worst;
        if (worst < tolerance) {
            out.converged = true;
            break;
        }
    }
    if (out.iterations == 0) {
        // No iteration ran, so `previous` still holds its sentinel and the
        // "eigenvalues" would be the largest representable double. Returning that
        // as an answer is worse than returning nothing: `craigBampton` would build
        // a modal block out of it. Nothing is what is returned.
        out.problem = "subspace iteration was given no iterations to run";
        return out;
    }
    if (!out.converged && out.problem.empty())
        // **The size of the change is the message.** An iteration that has hit the
        // floor `cond(K_ii)` puts under it -- the ordinary case on slender plating,
        // see `fixedInterfaceModes` in the header -- leaves this around 1e-8 while
        // the eigenvalues themselves are right to seven figures; a solve that is
        // genuinely lost leaves it at order one. Without the number the two read
        // identically, and a warning that cannot be triaged is a warning that gets
        // ignored, which is how this repo sat on a `ram_view` finding for months.
        out.problem = "subspace iteration did not converge in " + std::to_string(maxIterations) +
                      " iterations: the wanted eigenvalues still moved by " +
                      std::to_string(out.lastChange) + " relative, against a tolerance of " +
                      std::to_string(tolerance) +
                      " -- around 1e-8 is the accuracy floor a slender K_ii's conditioning puts"
                      " under the banded solve, and the values are still good to most of their"
                      " figures; order one is a lost solve";

    out.count = count;
    out.value.assign(previous.begin(), previous.end());
    out.vector.assign(ni * static_cast<std::size_t>(count), 0.0);
    for (int j = 0; j < count; ++j)
        for (std::size_t p = 0; p < ni; ++p)
            out.vector[static_cast<std::size_t>(j) * ni + p] = x[p * nq + static_cast<std::size_t>(j)];
    return out;
}

// --- The reduction -----------------------------------------------------------------

Reduction craigBampton(const Substructure& sub, const ReduceParams& params) {
    Reduction out;
    const double started = now();
    if (!sub.ready()) {
        out.problems.push_back("the substructure is not usable");
        for (const std::string& problem : sub.problems()) out.problems.push_back(problem);
        return out;
    }

    const std::size_t nb = sub.boundaryCount();
    const std::size_t ni = sub.interiorCount();
    out.boundary = static_cast<int>(nb);

    // --- How many modes. A frequency cutoff is answered exactly by an inertia
    // count, before a single eigenvector is computed -- see §2.
    int wanted = params.modes;
    if (wanted < 0) {
        bool exact = true;
        wanted = sub.eigenvaluesBelow(params.cutoffFrequency * params.cutoffFrequency, &exact);
        out.modesBelowCutoff = wanted;
        if (!exact)
            out.problems.push_back(
                "the inertia count hit a zero pivot, so the mode count below the cutoff is an "
                "estimate");
    }
    if (wanted > params.maxModes) {
        out.problems.push_back("clamped " + std::to_string(wanted) + " modes below the cutoff to " +
                               std::to_string(params.maxModes));
        wanted = params.maxModes;
    }
    if (static_cast<std::size_t>(wanted) > ni) wanted = static_cast<int>(ni);
    if (wanted < 0) wanted = 0;

    // --- Constraint modes: Psi = -K_ii^-1 K_ib, one banded solve per interface DOF.
    // K_ib's column b is K_bi's row b, by symmetry, and the row is what the CSR has
    // to hand.
    out.constraintModes.assign(ni * nb, 0.0);
    std::vector<double> rowBoundary, rowInterior;
    for (std::size_t b = 0; b < nb; ++b) {
        sub.boundaryRow(b, rowBoundary, rowInterior);
        for (double& v : rowInterior) v = -v;
        sub.interiorSolve(rowInterior);
        for (std::size_t p = 0; p < ni; ++p) out.constraintModes[p * nb + b] = rowInterior[p];
    }

    // --- Fixed-interface normal modes. One more than is kept when the Sturm check
    // is on, because the shift has to sit between the last kept and the first
    // discarded, and the first discarded is not otherwise known.
    const bool verify = params.verifyModes && wanted > 0 && static_cast<std::size_t>(wanted) < ni;
    const int asked = std::max(1, verify ? wanted + 1 : wanted);
    const Eigenpairs modes =
        sub.fixedInterfaceModes(asked, params.eigenTolerance, params.maxIterations);
    if (!modes.converged && !modes.problem.empty()) out.problems.push_back(modes.problem);
    if (modes.count > 0) out.firstFixedFrequency = std::sqrt(std::max(0.0, modes.value[0]));
    // A failed eigensolve returns no pairs at all, and the reduction must not go on
    // to read them: the first version of this walked off the end of an empty
    // spectrum, which is a crash where the right answer is a Guyan reduction and a
    // stated reason.
    if (modes.count < wanted) {
        out.problems.push_back("only " + std::to_string(modes.count) + " of " +
                               std::to_string(wanted) +
                               " fixed-interface modes were computed; the rest are dropped");
        wanted = modes.count;
    }

    if (verify && modes.count > wanted) {
        const double kept = modes.value[static_cast<std::size_t>(wanted) - 1];
        const double next = modes.value[static_cast<std::size_t>(wanted)];
        if (next > kept) {
            bool exact = true;
            const int below = sub.eigenvaluesBelow(0.5 * (kept + next), &exact);
            out.modesVerified = exact && below == wanted;
            if (!out.modesVerified)
                out.problems.push_back("the Sturm check counted " + std::to_string(below) +
                                       " eigenvalues below the last mode kept, not " +
                                       std::to_string(wanted) +
                                       ": subspace iteration skipped one");
        } else {
            out.problems.push_back(
                "the last mode kept and the first discarded are numerically equal, so the mode "
                "count cannot be verified by an inertia count");
        }
    }

    const int m = wanted;
    out.modes = m;
    out.normalModes.assign(ni * static_cast<std::size_t>(m), 0.0);
    out.frequency.assign(static_cast<std::size_t>(m), 0.0);
    for (int j = 0; j < m; ++j) {
        out.frequency[static_cast<std::size_t>(j)] =
            std::sqrt(std::max(0.0, modes.value[static_cast<std::size_t>(j)]));
        for (std::size_t p = 0; p < ni; ++p)
            out.normalModes[p * static_cast<std::size_t>(m) + static_cast<std::size_t>(j)] =
                modes.vector[static_cast<std::size_t>(j) * ni + p];
    }

    // --- The reduced pair.
    const std::size_t n = nb + static_cast<std::size_t>(m);
    out.stiffness.assign(n * n, 0.0);
    out.mass.assign(n * n, 0.0);

    // Khat = K_bb + K_bi Psi, upper triangle then mirrored. It is analytically
    // symmetric; mirroring makes it symmetric to the last bit, and the check that
    // it is the *right* matrix is `tests/test_reduction.cpp` forming T^T K T the
    // long way, not the symmetry.
    for (std::size_t b = 0; b < nb; ++b) {
        sub.boundaryRow(b, rowBoundary, rowInterior);
        for (std::size_t c = b; c < nb; ++c) out.stiffness[b * n + c] = rowBoundary[c];
        for (std::size_t p = 0; p < ni; ++p) {
            const double kbi = rowInterior[p];
            if (kbi == 0.0) continue;
            const double* psi = out.constraintModes.data() + p * nb;
            for (std::size_t c = b; c < nb; ++c) out.stiffness[b * n + c] += kbi * psi[c];
        }
    }
    // The modal stiffness block is Lambda exactly, and the coupling block is
    // identically zero: K_bi Phi + Psi^T K_ii Phi = 0. Both are stored as what they
    // are rather than computed, which is the whole reason the reduced eigenproblem
    // stays well conditioned.
    for (int j = 0; j < m; ++j) {
        const double omega = out.frequency[static_cast<std::size_t>(j)];
        out.stiffness[(nb + static_cast<std::size_t>(j)) * n + nb + static_cast<std::size_t>(j)] =
            omega * omega;
    }

    // M_bb + Psi^T M_ii Psi, and Psi^T M_ii Phi. M is diagonal, so M_bi vanishes
    // and there is no boundary-interior cross term to carry.
    const std::vector<double>& mass = sub.mass();
    for (std::size_t b = 0; b < nb; ++b)
        out.mass[b * n + b] += mass[sub.boundaryDof()[b]];
    for (std::size_t p = 0; p < ni; ++p) {
        const double mp = mass[sub.interiorDof()[p]];
        const double* psi = out.constraintModes.data() + p * nb;
        for (std::size_t b = 0; b < nb; ++b) {
            const double v = mp * psi[b];
            if (v == 0.0) continue;
            for (std::size_t c = b; c < nb; ++c) out.mass[b * n + c] += v * psi[c];
            for (int j = 0; j < m; ++j)
                out.mass[b * n + nb + static_cast<std::size_t>(j)] +=
                    v * out.normalModes[p * static_cast<std::size_t>(m) +
                                        static_cast<std::size_t>(j)];
        }
        // Phi^T M_ii Phi is the identity by the mass normalisation, but it is
        // computed rather than assumed: it is cheap, and asserting on it is what
        // says the eigenvectors really came back mass-normalised.
        for (int i = 0; i < m; ++i) {
            const double v = mp * out.normalModes[p * static_cast<std::size_t>(m) +
                                                  static_cast<std::size_t>(i)];
            for (int j = i; j < m; ++j)
                out.mass[(nb + static_cast<std::size_t>(i)) * n + nb + static_cast<std::size_t>(j)] +=
                    v * out.normalModes[p * static_cast<std::size_t>(m) +
                                        static_cast<std::size_t>(j)];
        }
    }

    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < i; ++j) {
            out.stiffness[i * n + j] = out.stiffness[j * n + i];
            out.mass[i * n + j] = out.mass[j * n + i];
        }

    out.reduceSeconds = now() - started;
    return out;
}

// --- Using it -----------------------------------------------------------------------

std::vector<double> reduceLoad(const Substructure& sub, const Reduction& reduced,
                               const std::vector<double>& load) {
    const std::size_t nb = static_cast<std::size_t>(reduced.boundary);
    const std::size_t m = static_cast<std::size_t>(reduced.modes);
    const std::size_t ni = sub.interiorCount();
    std::vector<double> out(nb + m, 0.0);
    if (load.size() < sub.dofCount()) return out;

    for (std::size_t b = 0; b < nb; ++b) out[b] = load[sub.boundaryDof()[b]];
    for (std::size_t p = 0; p < ni; ++p) {
        const double f = load[sub.interiorDof()[p]];
        if (f == 0.0) continue;
        for (std::size_t b = 0; b < nb; ++b) out[b] += f * reduced.constraintModes[p * nb + b];
        for (std::size_t j = 0; j < m; ++j)
            out[nb + j] += f * reduced.normalModes[p * m + j];
    }
    return out;
}

std::vector<double> recover(const Substructure& sub, const Reduction& reduced,
                            const std::vector<double>& state) {
    const std::size_t nb = static_cast<std::size_t>(reduced.boundary);
    const std::size_t m = static_cast<std::size_t>(reduced.modes);
    const std::size_t ni = sub.interiorCount();
    std::vector<double> out(sub.dofCount(), 0.0);
    if (state.size() < nb + m) return out;

    for (std::size_t b = 0; b < nb; ++b) out[sub.boundaryDof()[b]] = state[b];
    for (std::size_t p = 0; p < ni; ++p) {
        double sum = 0.0;
        for (std::size_t b = 0; b < nb; ++b) sum += reduced.constraintModes[p * nb + b] * state[b];
        for (std::size_t j = 0; j < m; ++j) sum += reduced.normalModes[p * m + j] * state[nb + j];
        out[sub.interiorDof()[p]] = sum;
    }
    // A degree of freedom an `Attachment::constrained` eliminated is in neither
    // partition, so nothing above has written it. Leaving it at zero is not a
    // small error: it is a hole in the middle of the displacement field, and every
    // element and every fibre touching that node then sees an artificial gradient
    // across one element length. Measured on a plate with one interior node tied
    // to two neighbours, `checkValidity` came back with 850 788 MPa against the
    // 427 MPa the same state really carries -- a factor of 1992 -- and it would have
    // done so on every tied section `section.cpp` builds. Same call `solveStatic`
    // already makes for the same reason.
    sub.expansion().recover(out);
    return out;
}

bool staticSolve(const Reduction& reduced, const std::vector<double>& load,
                 const std::vector<std::uint32_t>& held, std::vector<double>& state,
                 std::string* problem) {
    const std::size_t n = static_cast<std::size_t>(reduced.size());
    state.assign(n, 0.0);
    if (n == 0) {
        if (problem) *problem = "the reduction is empty";
        return false;
    }
    std::vector<std::uint8_t> fixed(n, 0);
    for (std::uint32_t d : held)
        if (d < n) fixed[d] = 1;

    std::vector<std::ptrdiff_t> map(n, -1);
    std::size_t free = 0;
    for (std::size_t d = 0; d < n; ++d)
        if (!fixed[d]) map[d] = static_cast<std::ptrdiff_t>(free++);
    if (free == 0) return true;

    std::vector<double> a(free * free, 0.0), rhs(free, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        if (map[i] < 0) continue;
        const std::size_t ii = static_cast<std::size_t>(map[i]);
        rhs[ii] = i < load.size() ? load[i] : 0.0;
        for (std::size_t j = 0; j < n; ++j) {
            if (map[j] < 0) continue;
            a[ii * free + static_cast<std::size_t>(map[j])] = reduced.stiffness[i * n + j];
        }
    }
    if (!cholesky(a, static_cast<int>(free))) {
        if (problem)
            *problem = "the reduced stiffness is singular with what is held: a rigid body mode is "
                       "free";
        return false;
    }
    forwardSolve(a, static_cast<int>(free), rhs.data());
    backwardSolve(a, static_cast<int>(free), rhs.data());
    for (std::size_t d = 0; d < n; ++d)
        if (map[d] >= 0) state[d] = rhs[static_cast<std::size_t>(map[d])];
    return true;
}

std::vector<double> reducedFrequencies(const Reduction& reduced,
                                       const std::vector<std::uint32_t>& held) {
    const std::size_t n = static_cast<std::size_t>(reduced.size());
    std::vector<std::uint8_t> fixed(n, 0);
    for (std::uint32_t d : held)
        if (d < n) fixed[d] = 1;
    std::vector<std::size_t> keep;
    for (std::size_t d = 0; d < n; ++d)
        if (!fixed[d]) keep.push_back(d);
    const std::size_t f = keep.size();
    if (f == 0) return {};

    std::vector<double> k(f * f), m(f * f);
    for (std::size_t i = 0; i < f; ++i)
        for (std::size_t j = 0; j < f; ++j) {
            k[i * f + j] = reduced.stiffness[keep[i] * n + keep[j]];
            m[i * f + j] = reduced.mass[keep[i] * n + keep[j]];
        }
    const Eigenpairs spectrum = generalisedEigen(k, m, static_cast<int>(f));
    std::vector<double> out;
    out.reserve(f);
    for (double value : spectrum.value) out.push_back(std::sqrt(std::max(0.0, value)));
    return out;
}

Validity checkValidity(const Substructure& sub, const Reduction& reduced,
                       const std::vector<double>& state) {
    Validity out;
    if (!sub.ready() || reduced.empty()) return out;
    const std::vector<double> u = recover(sub, reduced, state);
    const solidshell::HexMesh& mesh = sub.mesh();

    for (std::size_t d = 0; d + 2 < u.size(); d += 3) {
        const double magnitude =
            std::sqrt(u[d] * u[d] + u[d + 1] * u[d + 1] + u[d + 2] * u[d + 2]);
        out.peakDisplacement = std::max(out.peakDisplacement, magnitude);
    }

    const std::size_t elements = mesh.elementCount();
    for (std::size_t e = 0; e < elements; ++e) {
        double nodePos[kDof], displacement[kDof], stress[kGauss * 6];
        mesh.gather(e, mesh.position, nodePos);
        mesh.gather(e, u, displacement);
        solidshell::elementStress(nodePos, displacement, sub.material(), sub.formulation(), stress);
        for (int g = 0; g < kGauss; ++g) {
            const double mises = vonMises(stress + g * 6);
            if (mises > out.platingVonMises) {
                out.platingVonMises = mises;
                out.worstElement = static_cast<int>(e);
            }
        }
    }

    // The other half of the structure (§9). A fibre is an axial bar, so its stress
    // tensor has one non-zero entry and its von Mises is `|sigma|` identically --
    // the same equivalent stress the elements are measured by, not an analogue of
    // it. The magnitude is taken here and not in `memberStress`, which reports the
    // sign because which side of the neutral axis a fibre is on is information.
    const std::size_t members = sub.attachedMembers();
    for (std::size_t m = 0; m < members; ++m) {
        const double mises = std::fabs(sub.memberStress(m, u));
        if (mises > out.memberVonMises) {
            out.memberVonMises = mises;
            out.worstMember = static_cast<int>(m);
        }
    }

    // With no members this is `platingVonMises` to the last bit, which is what
    // makes the whole of §9 inert on a substructure that has none.
    out.peakVonMises = std::max(out.platingVonMises, out.memberVonMises);
    const double yield = sub.material().yieldStrength;
    out.utilisation = yield > 0.0 ? out.peakVonMises / yield : 0.0;
    out.linear = out.utilisation < 1.0;
    return out;
}

// --- Synthesis -------------------------------------------------------------------

namespace {

// Solve a dense symmetric system with some DOF held at zero. Shared by the one
// component and the assembled paths, which differ only in where the matrix lives.
bool solveHeld(const std::vector<double>& stiffness, std::size_t n,
               const std::vector<double>& load, const std::vector<std::uint32_t>& held,
               std::vector<double>& state, std::string* problem) {
    state.assign(n, 0.0);
    if (n == 0) {
        if (problem) *problem = "the model is empty";
        return false;
    }
    std::vector<std::uint8_t> fixed(n, 0);
    for (std::uint32_t d : held)
        if (d < n) fixed[d] = 1;
    std::vector<std::ptrdiff_t> map(n, -1);
    std::size_t free = 0;
    for (std::size_t d = 0; d < n; ++d)
        if (!fixed[d]) map[d] = static_cast<std::ptrdiff_t>(free++);
    if (free == 0) return true;

    std::vector<double> a(free * free, 0.0), rhs(free, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        if (map[i] < 0) continue;
        const std::size_t ii = static_cast<std::size_t>(map[i]);
        rhs[ii] = i < load.size() ? load[i] : 0.0;
        for (std::size_t j = 0; j < n; ++j)
            if (map[j] >= 0)
                a[ii * free + static_cast<std::size_t>(map[j])] = stiffness[i * n + j];
    }
    if (!cholesky(a, static_cast<int>(free))) {
        if (problem)
            *problem = "the stiffness is singular with what is held: a rigid body mode is free";
        return false;
    }
    forwardSolve(a, static_cast<int>(free), rhs.data());
    backwardSolve(a, static_cast<int>(free), rhs.data());
    for (std::size_t d = 0; d < n; ++d)
        if (map[d] >= 0) state[d] = rhs[static_cast<std::size_t>(map[d])];
    return true;
}

// Scatter-add one component's reduced pair into the assembly. Shared boundary DOF
// land on the same row and column, which is the whole of the coupling: components
// meeting at an interface have the same displacement there, so those are the same
// unknown. Every route into an `Assembly` goes through this, so "assembling is
// addition" is a property of the code rather than of three copies of it.
void scatterInto(Assembly& out, const Reduction& r, const std::vector<int>& from) {
    const std::size_t n = static_cast<std::size_t>(out.size());
    const std::size_t m = static_cast<std::size_t>(r.size());
    for (std::size_t i = 0; i < m; ++i) {
        const int ri = from[i];
        if (ri < 0) continue;
        for (std::size_t j = 0; j < m; ++j) {
            const int rj = from[j];
            if (rj < 0) continue;
            out.stiffness[static_cast<std::size_t>(ri) * n + static_cast<std::size_t>(rj)] +=
                r.stiffness[i * m + j];
            out.mass[static_cast<std::size_t>(ri) * n + static_cast<std::size_t>(rj)] +=
                r.mass[i * m + j];
        }
    }
}

// Union-find over the boundary DOF of every component at once, which is what makes
// three components possible where two needed only a map: a DOF shared by A and B and
// by B and C is *one* assembled row, and no pairwise map says so.
struct DofClasses {
    std::vector<int> parent;
    explicit DofClasses(std::size_t n) : parent(n) {
        for (std::size_t i = 0; i < n; ++i) parent[i] = static_cast<int>(i);
    }
    int find(int a) {
        while (parent[static_cast<std::size_t>(a)] != a) {
            parent[static_cast<std::size_t>(a)] =
                parent[static_cast<std::size_t>(parent[static_cast<std::size_t>(a)])];
            a = parent[static_cast<std::size_t>(a)];
        }
        return a;
    }
    void join(int a, int b) {
        a = find(a);
        b = find(b);
        if (a != b) parent[static_cast<std::size_t>(a)] = b;
    }
};

}  // namespace

std::vector<BoundaryDof> boundaryIdentity(const Substructure& substructure) {
    std::vector<BoundaryDof> out;
    if (!substructure.ready()) return out;
    const std::vector<std::uint32_t>& dof = substructure.boundaryDof();
    const std::vector<double>& position = substructure.mesh().position;
    out.reserve(dof.size());
    for (std::uint32_t d : dof) {
        const std::size_t node = d / 3u;
        // A boundary DOF past the end of the mesh cannot happen from a ready
        // substructure, and coming back short rather than reading past the array is
        // what makes that a statement rather than an assumption.
        if (3 * node + 2 >= position.size()) return {};
        out.push_back({Vec3{position[3 * node], position[3 * node + 1], position[3 * node + 2]},
                       d % 3u});
    }
    return out;
}

InterfaceMap matchBoundaries(const std::vector<BoundaryDof>& a, const std::vector<BoundaryDof>& b,
                             double tolerance) {
    InterfaceMap out;
    out.aToB.assign(a.size(), -1);

    const double tol2 = tolerance * tolerance;
    std::vector<std::uint8_t> takenB(b.size(), 0);
    for (std::size_t i = 0; i < a.size(); ++i) {
        for (std::size_t j = 0; j < b.size(); ++j) {
            // The axis must agree. Matching on position alone would let a
            // coincident node couple x to y, and the assembled model would be
            // wrong in a way that still produces plausible numbers.
            if (takenB[j] || b[j].axis != a[i].axis) continue;
            const double dx = a[i].position.x - b[j].position.x;
            const double dy = a[i].position.y - b[j].position.y;
            const double dz = a[i].position.z - b[j].position.z;
            const double d2 = dx * dx + dy * dy + dz * dz;
            if (d2 > tol2) continue;
            out.aToB[i] = static_cast<int>(j);
            takenB[j] = 1;
            out.shared++;
            out.worstGap = std::max(out.worstGap, std::sqrt(d2));
            break;
        }
    }
    if (out.shared == 0)
        out.problems.push_back("the two substructures share no boundary DOF at this tolerance, so "
                               "assembling them would produce two independent models side by side");
    return out;
}

InterfaceMap matchBoundaries(const Substructure& a, const Substructure& b, double tolerance) {
    InterfaceMap out;
    if (!a.ready() || !b.ready()) {
        out.problems.push_back("a substructure is not ready, so its boundary means nothing");
        return out;
    }
    return matchBoundaries(boundaryIdentity(a), boundaryIdentity(b), tolerance);
}

InterfaceMap matchBoundaries(const Assembly& a, const Substructure& b, double tolerance) {
    InterfaceMap out;
    if (a.boundaryPoint.size() != static_cast<std::size_t>(a.boundary)) {
        out.problems.push_back("the assembly carries no boundary DOF identity, so there is nothing "
                               "to match against: it was built from reductions alone");
        return out;
    }
    if (!b.ready()) {
        out.problems.push_back("a substructure is not ready, so its boundary means nothing");
        return out;
    }
    return matchBoundaries(a.boundaryPoint, boundaryIdentity(b), tolerance);
}

const std::vector<int>& Assembly::fromA() const {
    static const std::vector<int> none;
    return from.size() > 0 ? from[0] : none;
}

const std::vector<int>& Assembly::fromB() const {
    static const std::vector<int> none;
    return from.size() > 1 ? from[1] : none;
}

Assembly assemble(const Reduction& a, const Reduction& b, const InterfaceMap& map) {
    Assembly out;
    if (a.empty() || b.empty()) {
        out.problems.push_back("a component reduction is empty");
        return out;
    }
    if (map.aToB.size() != static_cast<std::size_t>(a.boundary)) {
        out.problems.push_back("the interface map does not describe component A's boundary");
        return out;
    }
    for (int i : map.aToB)
        if (i >= b.boundary) {
            out.problems.push_back("the interface map points past component B's boundary");
            return out;
        }

    // A's boundary first, then B's unshared boundary, then the modal blocks.
    std::vector<int> bBoundary(static_cast<std::size_t>(b.boundary), -1);
    for (std::size_t i = 0; i < map.aToB.size(); ++i)
        if (map.aToB[i] >= 0) bBoundary[static_cast<std::size_t>(map.aToB[i])] =
                static_cast<int>(i);
    int next = a.boundary;
    for (std::size_t j = 0; j < bBoundary.size(); ++j)
        if (bBoundary[j] < 0) bBoundary[j] = next++;

    out.boundary = next;
    out.modes = a.modes + b.modes;
    out.parts = 2;
    const std::size_t n = static_cast<std::size_t>(out.size());
    out.stiffness.assign(n * n, 0.0);
    out.mass.assign(n * n, 0.0);

    out.from.assign(2, {});
    std::vector<int>& fromA = out.from[0];
    std::vector<int>& fromB = out.from[1];
    fromA.assign(static_cast<std::size_t>(a.size()), -1);
    for (int i = 0; i < a.boundary; ++i) fromA[static_cast<std::size_t>(i)] = i;
    for (int j = 0; j < a.modes; ++j)
        fromA[static_cast<std::size_t>(a.boundary + j)] = out.boundary + j;

    fromB.assign(static_cast<std::size_t>(b.size()), -1);
    for (int i = 0; i < b.boundary; ++i)
        fromB[static_cast<std::size_t>(i)] = bBoundary[static_cast<std::size_t>(i)];
    for (int j = 0; j < b.modes; ++j)
        fromB[static_cast<std::size_t>(b.boundary + j)] = out.boundary + a.modes + j;

    scatterInto(out, a, fromA);
    scatterInto(out, b, fromB);

    for (const std::string& p : map.problems) out.problems.push_back(p);
    return out;
}

std::vector<double> assembledFrequencies(const Assembly& assembly,
                                         const std::vector<std::uint32_t>& held) {
    const std::size_t n = static_cast<std::size_t>(assembly.size());
    std::vector<std::uint8_t> fixed(n, 0);
    for (std::uint32_t d : held)
        if (d < n) fixed[d] = 1;
    std::vector<std::size_t> keep;
    for (std::size_t d = 0; d < n; ++d)
        if (!fixed[d]) keep.push_back(d);
    const std::size_t f = keep.size();
    if (f == 0) return {};

    std::vector<double> k(f * f), m(f * f);
    for (std::size_t i = 0; i < f; ++i)
        for (std::size_t j = 0; j < f; ++j) {
            k[i * f + j] = assembly.stiffness[keep[i] * n + keep[j]];
            m[i * f + j] = assembly.mass[keep[i] * n + keep[j]];
        }
    const Eigenpairs spectrum = generalisedEigen(k, m, static_cast<int>(f));
    std::vector<double> out;
    out.reserve(f);
    for (double value : spectrum.value) out.push_back(std::sqrt(std::max(0.0, value)));
    return out;
}

bool assembledStaticSolve(const Assembly& assembly, const std::vector<double>& load,
                          const std::vector<std::uint32_t>& held, std::vector<double>& state,
                          std::string* problem) {
    return solveHeld(assembly.stiffness, static_cast<std::size_t>(assembly.size()), load, held,
                     state, problem);
}

bool assembledStaticSolve(const Assembly& assembly, const std::vector<double>& load,
                          const std::vector<std::uint32_t>& held,
                          const std::vector<Prescribed>& prescribed, std::vector<double>& state,
                          std::string* problem) {
    const std::size_t n = static_cast<std::size_t>(assembly.size());
    state.assign(n, 0.0);
    if (n == 0 || assembly.stiffness.size() != n * n) {
        if (problem) *problem = "the assembly is empty";
        return false;
    }

    // The prescribed field, zero everywhere else. Written by DOF rather than pushed,
    // so a caller naming the same DOF twice gets the last value instead of two
    // contributions to the right-hand side.
    std::vector<double> xp(n, 0.0);
    std::vector<std::uint8_t> isPrescribed(n, 0u);
    for (const Prescribed& p : prescribed) {
        if (p.dof >= n) {
            if (problem) *problem = "a prescribed DOF is past the end of the assembly";
            return false;
        }
        xp[p.dof] = p.value;
        isPrescribed[p.dof] = 1u;
    }

    // f_f - (K x_p)_f. The rows belonging to held or prescribed DOF are discarded by
    // the solve below, so they are computed and thrown away rather than skipped; the
    // assembly is dense and the branch would cost more than the multiply.
    std::vector<double> rhs(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        double sum = i < load.size() ? load[i] : 0.0;
        const double* row = &assembly.stiffness[i * n];
        for (std::size_t j = 0; j < n; ++j)
            if (isPrescribed[j]) sum -= row[j] * xp[j];
        rhs[i] = sum;
    }

    // A DOF that is both held and prescribed is prescribed: held is the special case
    // value zero, and the other resolution would discard the value silently.
    std::vector<std::uint32_t> fixed = held;
    fixed.reserve(fixed.size() + prescribed.size());
    for (const Prescribed& p : prescribed) fixed.push_back(p.dof);

    if (!assembledStaticSolve(assembly, rhs, fixed, state, problem)) return false;
    for (std::size_t d = 0; d < n; ++d)
        if (isPrescribed[d]) state[d] = xp[d];
    return true;
}

std::vector<double> componentState(const Assembly& assembly, const std::vector<int>& from,
                                   const std::vector<double>& state) {
    // The assembly is here to be checked against, not for symmetry. A state of the
    // wrong length is the one mistake this call cannot otherwise notice: every
    // index in `from` is in range for the assembly, so a short state would silently
    // return zeros for the tail and the caller would recover a plausible field that
    // is missing its modal content.
    if (state.size() != static_cast<std::size_t>(assembly.size())) return {};
    std::vector<double> out(from.size(), 0.0);
    for (std::size_t i = 0; i < from.size(); ++i) {
        const int a = from[i];
        if (a >= 0 && static_cast<std::size_t>(a) < state.size())
            out[i] = state[static_cast<std::size_t>(a)];
    }
    return out;
}

// --- Assembling many -------------------------------------------------------------

namespace {

// The joints of a component list, restricted to the pairs `wanted` accepts. One
// implementation so that "neighbours only" is a filter rather than a second match.
std::vector<Joint> matchPairs(const std::vector<Component>& parts, double tolerance,
                              bool neighboursOnly) {
    std::vector<Joint> out;
    std::vector<std::vector<BoundaryDof>> identity(parts.size());
    for (std::size_t c = 0; c < parts.size(); ++c)
        if (parts[c].substructure) identity[c] = boundaryIdentity(*parts[c].substructure);

    for (std::size_t a = 0; a < parts.size(); ++a)
        for (std::size_t b = a + 1; b < parts.size(); ++b) {
            if (neighboursOnly && b != a + 1) continue;
            if (identity[a].empty() || identity[b].empty()) continue;
            InterfaceMap map = matchBoundaries(identity[a], identity[b], tolerance);
            // Sharing nothing is the expected answer for two components that are not
            // neighbours, so the complaint `matchBoundaries` raises about it is not a
            // problem here and the joint is simply not made.
            if (map.shared == 0) continue;
            map.problems.clear();
            out.push_back({static_cast<int>(a), static_cast<int>(b), std::move(map)});
        }
    return out;
}

}  // namespace

std::vector<Joint> matchComponents(const std::vector<Component>& parts, double tolerance) {
    return matchPairs(parts, tolerance, false);
}

std::vector<Joint> matchNeighbours(const std::vector<Component>& parts, double tolerance) {
    return matchPairs(parts, tolerance, true);
}

Assembly assemble(const std::vector<Component>& parts, const std::vector<Joint>& joints) {
    Assembly out;
    if (parts.empty()) {
        out.problems.push_back("no components to assemble");
        return out;
    }
    for (const Component& part : parts)
        if (!part.reduced || part.reduced->empty()) {
            out.problems.push_back("a component reduction is empty");
            return out;
        }
    // A reduction that fell back to an empty one -- which `craigBampton` does rather
    // than throwing -- is how a map and the matrices it indexes come to disagree, so
    // the substructure a joint was matched against has to be the one that was reduced.
    for (const Component& part : parts)
        if (part.substructure &&
            part.reduced->boundary != static_cast<int>(part.substructure->boundaryCount())) {
            out.problems.push_back("a component reduction does not carry its substructure's "
                                   "boundary, so an interface map does not describe it");
            return out;
        }

    const std::size_t k = parts.size();
    std::vector<std::size_t> offset(k + 1, 0);
    for (std::size_t c = 0; c < k; ++c)
        offset[c + 1] = offset[c] + static_cast<std::size_t>(parts[c].reduced->boundary);
    const std::size_t total = offset[k];

    DofClasses classes(total);
    for (const Joint& joint : joints) {
        if (joint.a < 0 || joint.b < 0 || static_cast<std::size_t>(joint.a) >= k ||
            static_cast<std::size_t>(joint.b) >= k) {
            out.problems.push_back("a joint names a component that is not in the list");
            return out;
        }
        // A component joined to itself would merge two of its own boundary DOF into
        // one unknown, which is a constraint and not an assembly.
        if (joint.a == joint.b) {
            out.problems.push_back("a joint joins a component to itself");
            return out;
        }
        const std::size_t ca = static_cast<std::size_t>(joint.a), cb = static_cast<std::size_t>(joint.b);
        if (joint.map.aToB.size() != static_cast<std::size_t>(parts[ca].reduced->boundary)) {
            out.problems.push_back("a joint's interface map does not describe its component's "
                                   "boundary");
            return out;
        }
        for (std::size_t i = 0; i < joint.map.aToB.size(); ++i) {
            const int j = joint.map.aToB[i];
            if (j < 0) continue;
            if (j >= parts[cb].reduced->boundary) {
                out.problems.push_back("a joint's interface map points past its component's "
                                       "boundary");
                return out;
            }
            classes.join(static_cast<int>(offset[ca] + i),
                         static_cast<int>(offset[cb] + static_cast<std::size_t>(j)));
        }
    }

    // Number the classes in component order, so two components reproduce exactly the
    // ordering the two-component `assemble` has always produced.
    std::vector<int> row(total, -1);
    int boundary = 0;
    for (std::size_t d = 0; d < total; ++d) {
        const int root = classes.find(static_cast<int>(d));
        if (row[static_cast<std::size_t>(root)] < 0)
            row[static_cast<std::size_t>(root)] = boundary++;
        row[d] = row[static_cast<std::size_t>(root)];
    }
    out.boundary = boundary;
    out.parts = static_cast<int>(k);
    for (const Component& part : parts) out.modes += part.reduced->modes;

    const std::size_t n = static_cast<std::size_t>(out.size());
    out.stiffness.assign(n * n, 0.0);
    out.mass.assign(n * n, 0.0);

    out.from.assign(k, {});
    int nextMode = out.boundary;
    for (std::size_t c = 0; c < k; ++c) {
        const Reduction& r = *parts[c].reduced;
        std::vector<int>& from = out.from[c];
        from.assign(static_cast<std::size_t>(r.size()), -1);
        for (int i = 0; i < r.boundary; ++i)
            from[static_cast<std::size_t>(i)] = row[offset[c] + static_cast<std::size_t>(i)];
        for (int j = 0; j < r.modes; ++j)
            from[static_cast<std::size_t>(r.boundary + j)] = nextMode++;
    }

    // The identity, and the check that the merges it describes were merges of the
    // same physical DOF. `matchBoundaries` guarantees that; a hand-built `Joint` does
    // not, and an assembly that coupled x to y would still solve.
    bool haveIdentity = true;
    for (const Component& part : parts)
        if (!part.substructure) haveIdentity = false;
    if (haveIdentity) {
        out.boundaryPoint.assign(static_cast<std::size_t>(out.boundary), BoundaryDof{});
        std::vector<std::uint8_t> filled(static_cast<std::size_t>(out.boundary), 0u);
        for (std::size_t c = 0; c < k; ++c) {
            const std::vector<BoundaryDof> identity = boundaryIdentity(*parts[c].substructure);
            if (identity.size() != static_cast<std::size_t>(parts[c].reduced->boundary)) {
                out.problems.push_back("a component's boundary identity is not its boundary");
                out.boundaryPoint.clear();
                haveIdentity = false;
                break;
            }
            for (std::size_t i = 0; i < identity.size(); ++i) {
                const auto r = static_cast<std::size_t>(row[offset[c] + i]);
                if (!filled[r]) {
                    out.boundaryPoint[r] = identity[i];
                    filled[r] = 1u;
                    continue;
                }
                const BoundaryDof& had = out.boundaryPoint[r];
                if (had.axis != identity[i].axis) ++out.axisDisagreements;
                out.worstMergedGap =
                    std::max(out.worstMergedGap, length(had.position - identity[i].position));
            }
        }
    }
    if (out.axisDisagreements > 0)
        out.problems.push_back(std::to_string(out.axisDisagreements) +
                               " assembled boundary rows merge degrees of freedom along different "
                               "axes: the model couples one direction to another and still solves");

    for (std::size_t c = 0; c < k; ++c) scatterInto(out, *parts[c].reduced, out.from[c]);
    for (const Joint& joint : joints)
        for (const std::string& p : joint.map.problems) out.problems.push_back(p);
    return out;
}

Assembly assemble(const Assembly& a, const Component& b, const InterfaceMap& map) {
    Assembly out;
    if (a.empty()) {
        out.problems.push_back("the assembly is empty");
        return out;
    }
    if (!b.reduced || b.reduced->empty()) {
        out.problems.push_back("a component reduction is empty");
        return out;
    }
    if (b.substructure &&
        b.reduced->boundary != static_cast<int>(b.substructure->boundaryCount())) {
        out.problems.push_back("a component reduction does not carry its substructure's boundary, "
                               "so the interface map does not describe it");
        return out;
    }
    if (map.aToB.size() != static_cast<std::size_t>(a.boundary)) {
        out.problems.push_back("the interface map does not describe the assembly's boundary");
        return out;
    }
    for (int j : map.aToB)
        if (j >= b.reduced->boundary) {
            out.problems.push_back("the interface map points past the component's boundary");
            return out;
        }

    // Every existing assembled row keeps its index: the assembly's boundary first,
    // then the component's unshared boundary, then the assembly's modal block, then
    // the component's. That is what makes the fold usable by a caller that already
    // holds indices into the assembly it started with.
    const Reduction& r = *b.reduced;
    std::vector<int> bBoundary(static_cast<std::size_t>(r.boundary), -1);
    for (std::size_t i = 0; i < map.aToB.size(); ++i)
        if (map.aToB[i] >= 0)
            bBoundary[static_cast<std::size_t>(map.aToB[i])] = static_cast<int>(i);
    int next = a.boundary;
    for (std::size_t j = 0; j < bBoundary.size(); ++j)
        if (bBoundary[j] < 0) bBoundary[j] = next++;

    out.boundary = next;
    out.modes = a.modes + r.modes;
    out.parts = a.parts + 1;
    out.worstMergedGap = a.worstMergedGap;
    out.axisDisagreements = a.axisDisagreements;
    const std::size_t n = static_cast<std::size_t>(out.size());
    out.stiffness.assign(n * n, 0.0);
    out.mass.assign(n * n, 0.0);

    // The assembly enters as one component whose reduced DOF are its own, remapped:
    // its boundary rows are unchanged and its modal block slides up by however many
    // boundary rows the new component brought.
    std::vector<int> fromAssembly(static_cast<std::size_t>(a.size()), -1);
    for (int i = 0; i < a.boundary; ++i) fromAssembly[static_cast<std::size_t>(i)] = i;
    for (int j = 0; j < a.modes; ++j)
        fromAssembly[static_cast<std::size_t>(a.boundary + j)] = out.boundary + j;

    out.from.assign(static_cast<std::size_t>(out.parts), {});
    for (int c = 0; c < a.parts && static_cast<std::size_t>(c) < a.from.size(); ++c) {
        std::vector<int> mapped = a.from[static_cast<std::size_t>(c)];
        for (int& d : mapped)
            if (d >= 0) d = fromAssembly[static_cast<std::size_t>(d)];
        out.from[static_cast<std::size_t>(c)] = std::move(mapped);
    }
    std::vector<int>& fromB = out.from[static_cast<std::size_t>(a.parts)];
    fromB.assign(static_cast<std::size_t>(r.size()), -1);
    for (int i = 0; i < r.boundary; ++i)
        fromB[static_cast<std::size_t>(i)] = bBoundary[static_cast<std::size_t>(i)];
    for (int j = 0; j < r.modes; ++j)
        fromB[static_cast<std::size_t>(r.boundary + j)] = out.boundary + a.modes + j;

    // The identity carries only when both sides have one -- an assembly built from
    // two bare `Reduction`s cannot acquire one by being folded.
    if (a.boundaryPoint.size() == static_cast<std::size_t>(a.boundary) && b.substructure) {
        const std::vector<BoundaryDof> identity = boundaryIdentity(*b.substructure);
        if (identity.size() == static_cast<std::size_t>(r.boundary)) {
            out.boundaryPoint.assign(static_cast<std::size_t>(out.boundary), BoundaryDof{});
            for (int i = 0; i < a.boundary; ++i)
                out.boundaryPoint[static_cast<std::size_t>(i)] =
                    a.boundaryPoint[static_cast<std::size_t>(i)];
            for (std::size_t i = 0; i < identity.size(); ++i) {
                const auto place = static_cast<std::size_t>(bBoundary[i]);
                if (static_cast<int>(place) >= a.boundary) {
                    out.boundaryPoint[place] = identity[i];
                    continue;
                }
                const BoundaryDof& had = out.boundaryPoint[place];
                if (had.axis != identity[i].axis) ++out.axisDisagreements;
                out.worstMergedGap =
                    std::max(out.worstMergedGap, length(had.position - identity[i].position));
            }
        }
    }

    // The assembly's own pair is scattered as though it were a reduction, which it is
    // -- boundary rows first, then modal, symmetric, dense. This is the copy the
    // N-way route does not make.
    Reduction asComponent;
    asComponent.boundary = a.boundary;
    asComponent.modes = a.modes;
    asComponent.stiffness = a.stiffness;
    asComponent.mass = a.mass;
    scatterInto(out, asComponent, fromAssembly);
    scatterInto(out, r, fromB);

    for (const std::string& p : a.problems) out.problems.push_back(p);
    for (const std::string& p : map.problems) out.problems.push_back(p);
    // The same complaint the N-way route makes, and it has to be made here too: an
    // assembly that counted a crossed axis and said nothing about it is an assembly
    // whose caller has to know to look.
    if (out.axisDisagreements > a.axisDisagreements)
        out.problems.push_back(std::to_string(out.axisDisagreements) +
                               " assembled boundary rows merge degrees of freedom along different "
                               "axes: the model couples one direction to another and still solves");
    return out;
}

int assembledComponents(const Assembly& assembly) {
    const std::size_t n = static_cast<std::size_t>(assembly.size());
    if (n == 0) return 0;
    DofClasses classes(n);
    for (const std::vector<int>& from : assembly.from) {
        int first = -1;
        for (int d : from) {
            if (d < 0 || static_cast<std::size_t>(d) >= n) continue;
            if (first < 0) {
                first = d;
                continue;
            }
            classes.join(first, d);
        }
    }
    std::vector<std::uint8_t> seen(n, 0u);
    int pieces = 0;
    for (std::size_t d = 0; d < n; ++d) {
        const auto root = static_cast<std::size_t>(classes.find(static_cast<int>(d)));
        if (!seen[root]) {
            seen[root] = 1u;
            ++pieces;
        }
    }
    return pieces;
}

}  // namespace sim::reduction
