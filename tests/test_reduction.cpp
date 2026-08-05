// SPDX-License-Identifier: MIT
//
// Validation of Tier 1: Craig-Bampton component mode synthesis.
//
// This is unusually well supplied with exact answers, so almost nothing here is a
// tolerance on an eyeballed number. Four things are *identities* and are asserted
// as such:
//
//   * **Zero modes is Guyan static condensation, and static condensation is exact
//     at the interface** -- for any load, at any mode count. What limits the
//     agreement is not the reduction but the conditioning of two independent
//     solves of the same system, and the tests measure that floor rather than
//     assuming one.
//   * **The rigid body modes survive.** A free-free substructure keeps exactly six
//     zero eigenvalues, `Psi` reproduces a rigid interface motion in the interior
//     exactly, and the reduced mass of a rigid translation is the substructure's
//     own mass.
//   * **The reduced pair is T^T K T and T^T M T.** Both are formed here the long
//     way, from a T built out of the constraint and normal modes, and compared to
//     what `craigBampton` assembled from closed forms. That is what makes the
//     symmetry check non-vacuous: the stored matrices are symmetric by
//     construction, so symmetry alone would prove nothing, and the content is that
//     they are the *right* symmetric matrices.
//   * **Reduced frequencies fall from above, monotonically.** The direction is
//     asserted, not the closeness: a reduction can only stiffen, so an
//     implementation whose error had the wrong sign would pass a "close enough"
//     test and fail this one.
//
// The vacuity guards matter as much as the assertions. A substructure with no
// interior proves nothing about a reduction; neither does one whose modes are
// never excited, nor a convergence study whose first point is already converged.
// Each of those is checked explicitly, with the negative control beside it.
#include "engine/sim/reduction.hpp"
#include "engine/sim/scantlings.hpp"
#include "engine/sim/solid_shell.hpp"
#include "engine/sim/zone.hpp"
#include "game/prototype/ferry.hpp"
#include "harness.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <numbers>
#include <string>
#include <vector>

using namespace sim;
using sim::reduction::Eigenpairs;
using sim::reduction::Plane;
using sim::reduction::ReduceParams;
using sim::reduction::Reduction;
using sim::reduction::Substructure;
using testing::expectEqual;
using testing::expectNear;
using testing::expectTrue;

namespace {

void expectEqualCount(const std::string& what, std::size_t got, std::size_t want) {
    expectEqual(what, static_cast<long long>(got), static_cast<long long>(want));
}

// --- The substructures under test -----------------------------------------------
//
// A plate is used rather than a ferry patch for everything with a closed form,
// because a flat rectangle's mass is `rho * lx * ly * t` and nothing else, and a
// rigid body mode of it is a rigid body mode exactly. The ferry's own plating gets
// the structural properties at the end, where the point is that they survive real
// geometry.

constexpr double kLx = 1.2, kLy = 0.6, kThickness = 0.012;
constexpr int kNx = 8, kNy = 4;

solidshell::HexMesh testPlate() {
    return solidshell::makePlateMesh(kLx, kLy, kThickness, kNx, kNy, 1);
}

// The same plate turned 45 degrees about z, for everything that looks at a
// **stress**. Nothing in the formulation is axis-aligned and the stress comes back
// in the global frame, so a stress running along this plate is equal parts
// sigma_xx, sigma_yy and sigma_xy -- which is what makes a von Mises formula
// testable. On the axis-aligned plate it is not.
const Vec3 kAxis{0.7071067811865476, 0.7071067811865475, 0.0};
constexpr double kEccentricLoad = 900.0;  // N per interior node, along the plate

solidshell::HexMesh rotatedPlate() {
    solidshell::HexMesh mesh = testPlate();
    for (std::size_t n = 0; n < mesh.nodeCount(); ++n) {
        const double x = mesh.position[n * 3], y = mesh.position[n * 3 + 1];
        mesh.position[n * 3] = x * kAxis.x - y * kAxis.y;
        mesh.position[n * 3 + 1] = x * kAxis.y + y * kAxis.x;
    }
    return mesh;
}

// The rotated plate's two ends, found with planes that are not axis-aligned
// either.
std::vector<std::uint32_t> rotatedEndInterface(const solidshell::HexMesh& mesh) {
    const std::vector<Plane> planes{{Vec3{0, 0, 0}, kAxis},
                                    {Vec3{kLx * kAxis.x, kLx * kAxis.y, 0}, kAxis}};
    return reduction::nodesNearPlanes(mesh, planes, 1e-9);
}

std::vector<std::uint32_t> heldAtRotatedOrigin(const Substructure& sub,
                                               const solidshell::HexMesh& mesh) {
    std::vector<std::uint32_t> held;
    for (std::size_t b = 0; b < sub.boundaryCount(); ++b) {
        const std::size_t node = sub.boundaryDof()[b] / 3;
        const double along = mesh.position[node * 3] * kAxis.x + mesh.position[node * 3 + 1] * kAxis.y;
        if (along < 0.5 * kLx) held.push_back(static_cast<std::uint32_t>(b));
    }
    return held;
}

// The two short edges: a strip of plating held at both ends, which is the shape a
// substructure cut out between two bulkheads has.
std::vector<std::uint32_t> endInterface(const solidshell::HexMesh& mesh) {
    const std::vector<Plane> planes{{Vec3{0, 0, 0}, Vec3{1, 0, 0}}, {Vec3{kLx, 0, 0}, Vec3{1, 0, 0}}};
    return reduction::nodesNearPlanes(mesh, planes, 1e-9);
}

// The reduced DOF that lie on the x = 0 edge, which is what a test holds when it
// wants a cantilever out of a two-ended substructure.
std::vector<std::uint32_t> heldAtOrigin(const Substructure& sub, const solidshell::HexMesh& mesh) {
    std::vector<std::uint32_t> held;
    for (std::size_t b = 0; b < sub.boundaryCount(); ++b)
        if (mesh.position[(sub.boundaryDof()[b] / 3) * 3] < 0.5 * kLx)
            held.push_back(static_cast<std::uint32_t>(b));
    return held;
}

// The full model's own natural frequencies, with `pinned` global DOF removed: a
// dense generalised eigensolve of the assembled pair, which is the reference every
// convergence claim here is made against. Independent of the reduction in
// everything except the assembly the two necessarily share.
std::vector<double> fullFrequencies(const Substructure& sub,
                                    const std::vector<std::uint8_t>& pinned) {
    std::vector<std::size_t> keep;
    for (std::size_t d = 0; d < sub.dofCount(); ++d)
        if (!pinned[d]) keep.push_back(d);
    const std::size_t f = keep.size();
    std::vector<double> k(f * f, 0.0), m(f * f, 0.0), x(sub.dofCount(), 0.0), y;
    for (std::size_t j = 0; j < f; ++j) {
        std::fill(x.begin(), x.end(), 0.0);
        x[keep[j]] = 1.0;
        sub.stiffnessTimes(x, y);
        for (std::size_t i = 0; i < f; ++i) k[i * f + j] = y[keep[i]];
    }
    for (std::size_t i = 0; i < f; ++i) m[i * f + i] = sub.mass()[keep[i]];
    const Eigenpairs spectrum = reduction::generalisedEigen(k, m, static_cast<int>(f));
    std::vector<double> out;
    out.reserve(spectrum.value.size());
    for (double value : spectrum.value) out.push_back(std::sqrt(std::max(0.0, value)));
    return out;
}

double misesOf(const double* s) {
    const double a = s[0] - s[1], b = s[1] - s[2], c = s[2] - s[0];
    return std::sqrt(0.5 * (a * a + b * b + c * c) +
                     3 * (s[3] * s[3] + s[4] * s[4] + s[5] * s[5]));
}

// The peak over every Gauss point of every element, and -- for the vacuity guards
// the mutation testing demanded -- two deliberately wrong versions of it: one that
// forgets the shear terms and one that looks only at the first Gauss point. A test
// that compares `checkValidity` against a stress state with no shear in it, or
// whose peak happens to sit at Gauss point zero, has not tested the formula.
struct StressPeak {
    double mises = 0;
    double withoutShear = 0;  // the same field, evaluated with the shear terms dropped
    double firstGaussOnly = 0;
};

StressPeak peakVonMises(const solidshell::HexMesh& mesh, const StructuralMaterial& material,
                        const std::vector<double>& displacement) {
    StressPeak peak;
    for (std::size_t e = 0; e < mesh.elementCount(); ++e) {
        double nodePos[24], disp[24], stress[48];
        mesh.gather(e, mesh.position, nodePos);
        mesh.gather(e, displacement, disp);
        solidshell::elementStress(nodePos, disp, material, solidshell::Formulation::SolidShell,
                                  stress);
        for (int g = 0; g < 8; ++g) {
            const double* s = stress + g * 6;
            const double m = misesOf(s);
            if (m > peak.mises) {
                peak.mises = m;
                const double a = s[0] - s[1], b = s[1] - s[2], c = s[2] - s[0];
                peak.withoutShear = std::sqrt(0.5 * (a * a + b * b + c * c));
            }
            if (g == 0) peak.firstGaussOnly = std::max(peak.firstGaussOnly, m);
        }
    }
    return peak;
}

// Volume and first moment of a trilinear hexahedron by 3 x 3 x 3 Gauss quadrature.
// Written out here rather than taken from `solid_shell.cpp` on purpose: the
// property being checked is that row-sum lumping puts the mass where the *geometry*
// says it is, and a check that shared the geometry routine would only be checking
// that a number equals itself.
void hexVolumeMoment(const double n[24], double& volume, double moment[3]) {
    static const double g[3] = {-0.7745966692414834, 0.0, 0.7745966692414834};
    static const double w[3] = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
    static const double xi[8] = {-1, 1, 1, -1, -1, 1, 1, -1};
    static const double eta[8] = {-1, -1, 1, 1, -1, -1, 1, 1};
    static const double zeta[8] = {-1, -1, -1, -1, 1, 1, 1, 1};
    volume = 0;
    moment[0] = moment[1] = moment[2] = 0;
    for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b)
            for (int c = 0; c < 3; ++c) {
                double shape[8], d[8][3], jac[3][3] = {}, point[3] = {};
                for (int i = 0; i < 8; ++i) {
                    const double sx = 1 + xi[i] * g[a], sy = 1 + eta[i] * g[b],
                                 sz = 1 + zeta[i] * g[c];
                    shape[i] = 0.125 * sx * sy * sz;
                    d[i][0] = 0.125 * xi[i] * sy * sz;
                    d[i][1] = 0.125 * eta[i] * sx * sz;
                    d[i][2] = 0.125 * zeta[i] * sx * sy;
                }
                for (int i = 0; i < 8; ++i)
                    for (int r = 0; r < 3; ++r) {
                        point[r] += shape[i] * n[i * 3 + r];
                        for (int q = 0; q < 3; ++q) jac[q][r] += d[i][q] * n[i * 3 + r];
                    }
                const double det = jac[0][0] * (jac[1][1] * jac[2][2] - jac[1][2] * jac[2][1]) -
                                   jac[0][1] * (jac[1][0] * jac[2][2] - jac[1][2] * jac[2][0]) +
                                   jac[0][2] * (jac[1][0] * jac[2][1] - jac[1][1] * jac[2][0]);
                const double weight = w[a] * w[b] * w[c] * det;
                volume += weight;
                for (int r = 0; r < 3; ++r) moment[r] += weight * point[r];
            }
}

// --- 1. The eigensolver, against closed forms -----------------------------------
//
// Everything downstream rests on this, and a symmetric eigensolver is one of the
// few things with an exactly known spectrum available: the second-difference
// operator's is `2 - 2 cos(k pi / (n + 1))` and its eigenvectors are sines.

void testDenseEigensolver() {
    std::printf("\n--- reduction: the eigensolver against closed forms ---\n");

    const int n = 24;
    std::vector<double> a(static_cast<std::size_t>(n) * n, 0.0);
    for (int i = 0; i < n; ++i) {
        a[static_cast<std::size_t>(i) * n + i] = 2.0;
        if (i + 1 < n) {
            a[static_cast<std::size_t>(i) * n + i + 1] = -1.0;
            a[static_cast<std::size_t>(i + 1) * n + i] = -1.0;
        }
    }
    const Eigenpairs spectrum = reduction::symmetricEigen(a, n);
    expectTrue("symmetricEigen converges", spectrum.converged);
    expectEqual("symmetricEigen returns every pair", spectrum.count, n);

    double worstValue = 0, worstVector = 0;
    for (int k = 0; k < n; ++k) {
        const double want = 2.0 - 2.0 * std::cos((k + 1) * std::numbers::pi / (n + 1));
        worstValue = std::max(worstValue, std::fabs(spectrum.value[static_cast<std::size_t>(k)] - want));
        // The eigenvector is sin(i (k+1) pi / (n+1)), normalised. Compared up to a
        // sign, because a sign convention is a choice and not a property.
        std::vector<double> want_v(static_cast<std::size_t>(n));
        double norm = 0;
        for (int i = 0; i < n; ++i) {
            want_v[static_cast<std::size_t>(i)] =
                std::sin((i + 1) * (k + 1) * std::numbers::pi / (n + 1));
            norm += want_v[static_cast<std::size_t>(i)] * want_v[static_cast<std::size_t>(i)];
        }
        norm = std::sqrt(norm);
        double dot = 0;
        for (int i = 0; i < n; ++i)
            dot += spectrum.mode(k)[i] * want_v[static_cast<std::size_t>(i)] / norm;
        worstVector = std::max(worstVector, std::fabs(std::fabs(dot) - 1.0));
    }
    expectNear("second-difference eigenvalues", worstValue, 0.0, 1e-13);
    // Two or three sweeps an eigenvalue is what a correctly shifted QL costs; an
    // unshifted or half-shifted one converges to the same answer and takes far
    // longer, which no assertion on the answer can see. The count is deterministic,
    // so it can be asserted where a wall clock could not.
    // 54 sweeps as written and 64 with the Wilkinson shift halved; dropping the
    // shift altogether does not converge at all. A *wrong* shift converges to the
    // same answer and only costs time, which no assertion on the answer can see --
    // but the count is deterministic, so it can be asserted where a wall clock
    // could not.
    expectTrue("and the Wilkinson shift is doing its job", spectrum.iterations <= 60);
    std::printf("     %d QL sweeps for %d eigenvalues\n", spectrum.iterations, n);
    expectNear("second-difference eigenvectors", worstVector, 0.0, 1e-12);

    // Orthonormality. A QL implementation that dropped the rotation accumulation
    // would still deliver the eigenvalues and nothing here else would notice.
    double worstOrth = 0;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            double d = 0;
            for (int k = 0; k < n; ++k) d += spectrum.mode(i)[k] * spectrum.mode(j)[k];
            worstOrth = std::max(worstOrth, std::fabs(d - (i == j ? 1.0 : 0.0)));
        }
    expectNear("eigenvectors are orthonormal", worstOrth, 0.0, 1e-13);

    // A **repeated** eigenvalue. A shifted QL that assumed distinct eigenvalues
    // would still get every test above right; a degenerate pair is where it stops.
    // Q diag(3, 3, 3, 7, 11) Q^T for a Householder Q, whose spectrum is known
    // exactly by construction.
    {
        const int d = 5;
        const double want[5] = {3, 3, 3, 7, 11};
        std::vector<double> v{1.0, 2.0, -1.0, 0.5, 3.0}, q(static_cast<std::size_t>(d) * d, 0.0);
        double vv = 0;
        for (double e : v) vv += e * e;
        for (int i = 0; i < d; ++i)
            for (int j = 0; j < d; ++j)
                q[static_cast<std::size_t>(i) * d + j] =
                    (i == j ? 1.0 : 0.0) - 2.0 * v[static_cast<std::size_t>(i)] *
                                               v[static_cast<std::size_t>(j)] / vv;
        std::vector<double> m(static_cast<std::size_t>(d) * d, 0.0);
        for (int i = 0; i < d; ++i)
            for (int j = 0; j < d; ++j) {
                double s = 0;
                for (int k = 0; k < d; ++k)
                    s += q[static_cast<std::size_t>(i) * d + k] * want[k] *
                         q[static_cast<std::size_t>(j) * d + k];
                m[static_cast<std::size_t>(i) * d + j] = s;
            }
        const Eigenpairs degenerate = reduction::symmetricEigen(m, d);
        double worst = 0;
        for (int k = 0; k < d; ++k)
            worst = std::max(worst, std::fabs(degenerate.value[static_cast<std::size_t>(k)] - want[k]));
        expectNear("a triply repeated eigenvalue is resolved", worst, 0.0, 1e-13);
        double orth = 0;
        for (int i = 0; i < d; ++i)
            for (int j = 0; j < d; ++j) {
                double dd = 0;
                for (int k = 0; k < d; ++k) dd += degenerate.mode(i)[k] * degenerate.mode(j)[k];
                orth = std::max(orth, std::fabs(dd - (i == j ? 1.0 : 0.0)));
            }
        expectNear("a degenerate eigenspace comes back orthonormal", orth, 0.0, 1e-13);
    }

    // The generalised solver: A x = lambda B x with B = 4I halves nothing and
    // quarters the spectrum, and the vectors come back B-orthonormal rather than
    // orthonormal -- which is the property the whole modal formulation rests on and
    // which an implementation that forgot the back substitution would fail.
    {
        std::vector<double> b(static_cast<std::size_t>(n) * n, 0.0);
        for (int i = 0; i < n; ++i) b[static_cast<std::size_t>(i) * n + i] = 4.0;
        const Eigenpairs g = reduction::generalisedEigen(a, b, n);
        expectTrue("generalisedEigen converges", g.converged);
        double worst = 0;
        for (int k = 0; k < n; ++k)
            worst = std::max(worst, std::fabs(g.value[static_cast<std::size_t>(k)] -
                                              0.25 * spectrum.value[static_cast<std::size_t>(k)]));
        expectNear("B = 4I quarters the spectrum", worst, 0.0, 1e-14);
        double worstB = 0;
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) {
                double d = 0;
                for (int k = 0; k < n; ++k) d += g.mode(i)[k] * 4.0 * g.mode(j)[k];
                worstB = std::max(worstB, std::fabs(d - (i == j ? 1.0 : 0.0)));
            }
        expectNear("generalised eigenvectors are B-orthonormal", worstB, 0.0, 1e-13);

        // Two negative controls, because they fail at different pivots: an
        // **indefinite** B gives a negative one, and a **singular** B gives an
        // exactly zero one. A Cholesky that accepts a zero pivot passes the first
        // and divides by zero on the second, which is a NaN with nothing to say
        // where it came from.
        std::vector<double> bad = b;
        bad[0] = -1.0;
        const Eigenpairs refused = reduction::generalisedEigen(a, bad, n);
        expectTrue("generalisedEigen refuses an indefinite right-hand matrix",
                   !refused.converged && !refused.problem.empty());
        // The zero goes on the **last** diagonal, and that is not arbitrary. Put it
        // anywhere else and the row below divides by it, produces a NaN, and the
        // very next pivot test rejects the matrix anyway -- so a Cholesky that
        // accepted a zero pivot would still look like one that refused it. On the
        // last diagonal there is no row below to catch it, and the failure escapes
        // into the solve.
        std::vector<double> singular = b;
        singular[static_cast<std::size_t>(n - 1) * n + (n - 1)] = 0.0;
        const Eigenpairs alsoRefused = reduction::generalisedEigen(a, singular, n);
        expectTrue("generalisedEigen refuses a singular right-hand matrix",
                   !alsoRefused.converged && !alsoRefused.problem.empty());
        // The *reason* matters and not merely the refusal: a Cholesky that accepted
        // a zero pivot would divide by it, fill the matrix with NaN, and be refused
        // a few lines later because the QL failed to converge -- which looks like
        // the same answer and is not.
        expectTrue("and names the pivot rather than the QL",
                   refused.problem.find("positive definite") != std::string::npos &&
                       alsoRefused.problem.find("positive definite") != std::string::npos);
    }
}

// --- 2. Assembly and partition ---------------------------------------------------

void testSubstructure() {
    std::printf("\n--- reduction: what the substructure assembled ---\n");
    solidshell::HexMesh mesh = testPlate();
    const StructuralMaterial steel = ah36Steel();
    Substructure sub(mesh, steel, endInterface(mesh));
    expectTrue("the plate substructure is usable", sub.ready());

    // The lumped mass is the plate's, exactly: row-sum lumping integrates the shape
    // functions, which sum to one, so the total is rho * volume however distorted
    // the elements are.
    expectNear("lumped mass is rho * lx * ly * t", sub.totalMass(),
               steel.density * kLx * kLy * kThickness, 1e-9);

    // Vacuity guard on every reduction claim that follows: there has to be an
    // interior worth reducing. Two nodes through the thickness at each of 9 x 5
    // stations is 270 DOF, of which the two end columns are the interface.
    expectEqualCount("plate DOF", sub.dofCount(), 270);
    // `makePlateMesh` numbers nodes so the band is already small, and reverse
    // Cuthill-McKee on it comes out *worse* -- 65 against 41. That is why the
    // ordering is the narrower of the two rather than the renumbering, and this is
    // the assertion that says so: a reduction that always renumbered would fail it.
    expectTrue("a mesh already well numbered keeps its numbering", sub.halfBandwidth() <= 45);
    expectEqualCount("interface DOF", sub.boundaryCount(), 60);
    expectEqualCount("interior DOF", sub.interiorCount(), 210);
    expectTrue("the interior is most of the substructure", sub.interiorCount() > sub.boundaryCount());

    // The assembled stiffness carries all six rigid body modes. This is a property
    // of the element, but it is *this* assembly that has to preserve it, and a
    // scatter that transposed a block or lost a term would not.
    {
        std::vector<double> u(sub.dofCount(), 0.0), f;
        double scale = 0;
        for (int axis = 0; axis < 3; ++axis) {
            std::fill(u.begin(), u.end(), 0.0);
            for (std::size_t d = static_cast<std::size_t>(axis); d < u.size(); d += 3) u[d] = 1.0;
            sub.stiffnessTimes(u, f);
            double worst = 0;
            for (double e : f) worst = std::max(worst, std::fabs(e));
            // The scale to compare against: the force a *unit strain* would raise,
            // which is what "zero" has to be small against.
            std::fill(u.begin(), u.end(), 0.0);
            for (std::size_t node = 0; node * 3 < u.size(); ++node)
                u[node * 3] = mesh.position[node * 3];
            std::vector<double> strained;
            sub.stiffnessTimes(u, strained);
            for (double e : strained) scale = std::max(scale, std::fabs(e));
            expectTrue("rigid translation carries no force (axis " + std::to_string(axis) + ")",
                       worst < 1e-12 * scale);
        }
        // And a rotation, which is the one a co-rotational formulation gets right
        // for a different reason and a linear one gets right only to first order.
        std::fill(u.begin(), u.end(), 0.0);
        for (std::size_t node = 0; node * 3 < u.size(); ++node) {
            const double x = mesh.position[node * 3], y = mesh.position[node * 3 + 1];
            u[node * 3] = -y;
            u[node * 3 + 1] = x;
        }
        sub.stiffnessTimes(u, f);
        double worst = 0;
        for (double e : f) worst = std::max(worst, std::fabs(e));
        expectTrue("infinitesimal rotation carries no force", worst < 1e-12 * scale);
    }

    // The mesh's own pins are ignored, and saying so is the contract. A constraint
    // silently consumed would make the free-free property below untestable.
    {
        solidshell::HexMesh pinned = testPlate();
        pinned.pin(0, 0, 0.0);
        Substructure withPins(pinned, steel, endInterface(pinned));
        bool reported = false;
        for (const std::string& problem : withPins.problems())
            if (problem.find("ignored") != std::string::npos) reported = true;
        expectTrue("a substructure says it ignored the mesh's pins", reported);
        expectTrue("and is still usable", withPins.ready());
    }

    // Negative control: an interface that does not restrain the six rigid body
    // modes leaves K_ii singular, and the failure has to be named rather than
    // discovered as a NaN twenty steps later. **The factorisation does not catch
    // this** -- a mechanism leaves a tiny *positive* pivot in floating point, so
    // `BandedSpd::factor` succeeds and the solve returns nonsense. The two ways an
    // interface can be degenerate are separated because they say different things
    // to whoever chose it.
    {
        Substructure tooFew(mesh, steel, {0u, 1u});
        expectTrue("an interface of one mid-surface point is refused", !tooFew.ready());
        bool named = false, counted = false;
        for (const std::string& problem : tooFew.problems()) {
            if (problem.find("rigid body") != std::string::npos) named = true;
            if (problem.find("fewer than three nodes") != std::string::npos) counted = true;
        }
        expectTrue("and the reason names the rigid body modes", named);
        expectTrue("and says it was a node count", counted);

        // Nodes are not enough if they lie on a line, however many there are: the
        // rotation about that line is still free. A whole row of them along one
        // edge of the lower face is exactly that, and it is the shape a caller who
        // counted nodes rather than looking at them would produce.
        std::vector<std::uint32_t> line;
        for (std::size_t n = 0; n < mesh.nodeCount(); ++n)
            if (std::fabs(mesh.position[n * 3 + 1]) < 1e-12 && mesh.position[n * 3 + 2] < 0)
                line.push_back(static_cast<std::uint32_t>(n));
        expectTrue("the collinear interface has plenty of nodes", line.size() >= 5);
        Substructure collinear(mesh, steel, line);
        expectTrue("a collinear interface is refused however many nodes it has",
                   !collinear.ready());
        bool saidCollinear = false;
        for (const std::string& problem : collinear.problems())
            if (problem.find("the interface is collinear") != std::string::npos) saidCollinear = true;
        expectTrue("and says so", saidCollinear);

        // **Nearly** collinear is the case that matters, because exactly collinear
        // is what a synthetic mesh produces and a real one never does. A node ten
        // femtometres off the line leaves the rotation about it restrained by
        // nothing worth having, and a tolerance compared against zero rather than
        // against the interface's own extent would accept it.
        solidshell::HexMesh nudged = testPlate();
        nudged.position[line[2] * 3 + 1] += 1e-14;
        Substructure nearly(nudged, steel, line);
        expectTrue("and so is an interface collinear to within a rounding of its own length",
                   !nearly.ready());
        // For the *stated* reason. A tolerance compared against zero would call it
        // non-collinear, hand it to the factorisation, and have that fail instead --
        // which is the same refusal with a reason that points at the wrong thing.
        bool stillCollinear = false;
        for (const std::string& problem : nearly.problems())
            if (problem.find("the interface is collinear") != std::string::npos) stillCollinear = true;
        expectTrue("for the reason that is actually true of it", stillCollinear);
    }
    {
        std::vector<std::uint32_t> everything(mesh.nodeCount());
        for (std::size_t n = 0; n < mesh.nodeCount(); ++n)
            everything[n] = static_cast<std::uint32_t>(n);
        Substructure whole(mesh, steel, everything);
        expectTrue("a substructure with no interior is refused", !whole.ready());
    }

    // The answer must not depend on the order the mesher happened to emit nodes in.
    // The interior is renumbered by reverse Cuthill-McKee before it is banded, and
    // a renumbering that lost a node or mismapped one would still factor and would
    // still be symmetric -- but it would not give the same frequencies.
    {
        solidshell::HexMesh shuffled = testPlate();
        const std::size_t nodes = shuffled.nodeCount();
        std::vector<std::uint32_t> newOf(nodes);
        for (std::size_t n = 0; n < nodes; ++n)
            newOf[n] = static_cast<std::uint32_t>((n * 37 + 11) % nodes);
        // (n * 37 + 11) mod 90 is a permutation because gcd(37, 90) = 1.
        std::vector<std::uint8_t> hit(nodes, 0);
        for (std::size_t n = 0; n < nodes; ++n) hit[newOf[n]] = 1;
        bool permutation = true;
        for (std::uint8_t h : hit) permutation = permutation && h;
        expectTrue("the shuffle really is a permutation", permutation);

        solidshell::HexMesh out = shuffled;
        for (std::size_t n = 0; n < nodes; ++n)
            for (int k = 0; k < 3; ++k)
                out.position[newOf[n] * 3 + static_cast<std::size_t>(k)] =
                    shuffled.position[n * 3 + static_cast<std::size_t>(k)];
        for (std::size_t i = 0; i < out.index.size(); ++i) out.index[i] = newOf[shuffled.index[i]];

        Substructure other(out, steel, endInterface(out));
        expectTrue("the shuffled substructure is usable", other.ready());
        ReduceParams params;
        params.modes = 6;
        const Reduction a = craigBampton(sub, params);
        const Reduction b = craigBampton(other, params);
        double worst = 0;
        for (int j = 0; j < a.modes; ++j)
            worst = std::max(worst,
                             std::fabs(a.frequency[static_cast<std::size_t>(j)] -
                                       b.frequency[static_cast<std::size_t>(j)]) /
                                 a.frequency[static_cast<std::size_t>(j)]);
        // 1e-7, not the 1e-10 this first asked for. The two assemblies sum the same
        // terms in different orders and the interiors are banded in different
        // renumberings, so the floor is the conditioning of a plate stiffness at
        // a/t = 100 -- around 1e-7 in an eigenvalue. The measured shift is printed
        // beside it, and it sits three orders under that floor.
        expectTrue("node numbering does not move the fixed-interface frequencies", worst < 1e-7);
        std::printf("     worst frequency shift under a node permutation: %.2e\n", worst);
    }
}

// --- 3. Zero modes is Guyan, and Guyan is exact at the interface -----------------

void testGuyanIsExact() {
    std::printf("\n--- reduction: static condensation is exact, and stays exact ---\n");
    solidshell::HexMesh mesh = testPlate();
    const StructuralMaterial steel = ah36Steel();
    Substructure sub(mesh, steel, endInterface(mesh));
    ReduceParams zero;
    zero.modes = 0;
    const Reduction guyan = craigBampton(sub, zero);
    expectEqual("zero modes keeps only the interface", guyan.size(), 60);

    // Psi applied to a rigid interface translation *is* the rigid interior field.
    // This is the sharpest statement about the sign of Psi there is: get it
    // backwards and the interior moves the other way.
    {
        std::vector<double> state(static_cast<std::size_t>(guyan.size()), 0.0);
        for (std::size_t b = 0; b < sub.boundaryCount(); ++b)
            if (sub.boundaryDof()[b] % 3 == 2) state[b] = 1.0;
        const std::vector<double> u = reduction::recover(sub, guyan, state);
        double worst = 0;
        for (std::size_t d = 0; d < u.size(); ++d)
            worst = std::max(worst, std::fabs(u[d] - (d % 3 == 2 ? 1.0 : 0.0)));
        expectTrue("Psi carries a rigid interface translation into the interior exactly",
                   worst < 1e-8);
        std::printf("     rigid interior recovery error: %.2e m of 1 m\n", worst);
    }

    // The static identity, against `solidshell::solveStatic` -- a different
    // assembly, a different numbering and a different factorisation of the same
    // problem.
    solidshell::HexMesh full = testPlate();
    std::vector<double> load(sub.dofCount(), 0.0);
    const std::vector<std::uint32_t> held = heldAtOrigin(sub, mesh);
    std::vector<std::uint8_t> pinned(sub.dofCount(), 0);
    for (std::uint32_t b : held) {
        const std::uint32_t d = sub.boundaryDof()[b];
        pinned[d] = 1;
        full.pin(d / 3, static_cast<int>(d % 3), 0.0);
    }
    for (std::size_t b = 0; b < sub.boundaryCount(); ++b) {
        const std::uint32_t d = sub.boundaryDof()[b];
        if (!pinned[d] && d % 3 == 2) load[d] = 1000.0;
    }
    std::vector<double> uFull;
    std::string problem;
    expectTrue("the full model solves", solidshell::solveStatic(full, steel,
                                                                solidshell::Formulation::SolidShell,
                                                                load, uFull, &problem));

    double peak = 0;
    for (double e : uFull) peak = std::max(peak, std::fabs(e));
    expectTrue("the load actually moves the plate", peak > 1e-3);

    std::vector<double> state;
    expectTrue("the reduced model solves",
               reduction::staticSolve(guyan, reduction::reduceLoad(sub, guyan, load), held, state,
                                      &problem));
    double worstBoundary = 0, worstInterior = 0;
    const std::vector<double> u = reduction::recover(sub, guyan, state);
    for (std::size_t b = 0; b < sub.boundaryCount(); ++b)
        worstBoundary = std::max(worstBoundary,
                                 std::fabs(u[sub.boundaryDof()[b]] - uFull[sub.boundaryDof()[b]]));
    for (std::size_t p = 0; p < sub.interiorCount(); ++p)
        worstInterior = std::max(worstInterior,
                                 std::fabs(u[sub.interiorDof()[p]] - uFull[sub.interiorDof()[p]]));
    // In exact arithmetic these are zero. In double they are limited by the
    // conditioning of a plate stiffness at a/t = 100, which is what the printed
    // figure measures; 1e-6 of the peak is two orders of margin over it and three
    // orders below any modelling error a reduction could make.
    expectTrue("zero modes reproduces the full interface response", worstBoundary < 1e-6 * peak);
    expectTrue("and the interior too, for a load applied at the interface",
               worstInterior < 1e-6 * peak);
    std::printf("     boundary %.2e m, interior %.2e m, of a peak %.3e m\n", worstBoundary,
                worstInterior, peak);

    // The correction to the usual statement of this property: with the load moved
    // *into* the interior, the interface response is **still** exact at zero modes
    // and does not improve as modes are added. What improves is the interior.
    {
        solidshell::HexMesh interiorLoaded = testPlate();
        for (std::uint32_t b : held) {
            const std::uint32_t d = sub.boundaryDof()[b];
            interiorLoaded.pin(d / 3, static_cast<int>(d % 3), 0.0);
        }
        std::vector<double> interiorLoad(sub.dofCount(), 0.0);
        for (std::size_t p = 0; p < sub.interiorCount(); ++p) {
            const std::uint32_t d = sub.interiorDof()[p];
            if (d % 3 == 2) interiorLoad[d] = 50.0;
        }
        std::vector<double> reference;
        solidshell::solveStatic(interiorLoaded, steel, solidshell::Formulation::SolidShell,
                                interiorLoad, reference, &problem);
        double peakB = 0, peakI = 0;
        for (std::size_t b = 0; b < sub.boundaryCount(); ++b)
            peakB = std::max(peakB, std::fabs(reference[sub.boundaryDof()[b]]));
        for (std::size_t p = 0; p < sub.interiorCount(); ++p)
            peakI = std::max(peakI, std::fabs(reference[sub.interiorDof()[p]]));
        expectTrue("the interior load moves the interface at all", peakB > 1e-4);

        double boundaryAtZero = 0, interiorAtZero = 0, boundaryAtMany = 0, interiorAtMany = 0;
        for (int m : {0, 32}) {
            ReduceParams params;
            params.modes = m;
            const Reduction r = craigBampton(sub, params);
            std::vector<double> s;
            reduction::staticSolve(r, reduction::reduceLoad(sub, r, interiorLoad), held, s, &problem);
            const std::vector<double> got = reduction::recover(sub, r, s);
            double eb = 0, ei = 0;
            for (std::size_t b = 0; b < sub.boundaryCount(); ++b)
                eb = std::max(eb, std::fabs(got[sub.boundaryDof()[b]] -
                                            reference[sub.boundaryDof()[b]]));
            for (std::size_t p = 0; p < sub.interiorCount(); ++p)
                ei = std::max(ei, std::fabs(got[sub.interiorDof()[p]] -
                                            reference[sub.interiorDof()[p]]));
            (m == 0 ? boundaryAtZero : boundaryAtMany) = eb;
            (m == 0 ? interiorAtZero : interiorAtMany) = ei;
        }
        expectTrue("an interior load leaves the interface exact at zero modes",
                   boundaryAtZero < 1e-6 * peakB);
        expectTrue("and thirty-two modes do not improve it",
                   boundaryAtMany < 1e-6 * peakB &&
                       boundaryAtMany > 0.1 * boundaryAtZero);
        // The negative control for the whole convergence story: at zero modes the
        // interior really is wrong, by a margin nothing could mistake for rounding.
        expectTrue("the interior at zero modes is wrong by percent, not by rounding",
                   interiorAtZero > 0.01 * peakI);
        expectTrue("and thirty-two modes fix it", interiorAtMany < 1e-4 * interiorAtZero);
        std::printf("     interior load: interface %.2e -> %.2e m (unchanged), interior %.2e -> "
                    "%.2e m of %.3e\n",
                    boundaryAtZero, boundaryAtMany, interiorAtZero, interiorAtMany, peakI);
    }
}

// --- 4. Rigid body modes ----------------------------------------------------------

void testRigidBody() {
    std::printf("\n--- reduction: a free component keeps its six rigid body modes ---\n");
    solidshell::HexMesh mesh = testPlate();
    const StructuralMaterial steel = ah36Steel();
    // The interface is a ring through the middle of the plate, so the substructure
    // is genuinely free-free: nothing anywhere holds it. That also leaves the
    // interior in *two* disconnected halves, which is the case an interior
    // renumbering can silently get wrong.
    const std::vector<Plane> mid{{Vec3{0.5 * kLx, 0, 0}, Vec3{1, 0, 0}}};
    Substructure sub(mesh, steel, reduction::nodesNearPlanes(mesh, mid, 1e-9));
    expectTrue("a free-free substructure cut at midspan is usable", sub.ready());
    expectTrue("its interior is in two pieces and both are reduced",
               sub.interiorCount() == sub.dofCount() - sub.boundaryCount());
    // And nothing complained on the way. A breadth-first renumbering that stopped
    // after the first component would leave half the interior unvisited, which is
    // reported rather than silently narrowing the band it is compared against.
    expectTrue("a two-piece interior renumbers without complaint", sub.problems().empty());

    ReduceParams params;
    params.modes = 6;
    const Reduction reduced = craigBampton(sub, params);
    const std::vector<double> omega = reduction::reducedFrequencies(reduced, {});
    expectTrue("the free-free reduced model has something to say",
               omega.size() >= 8);

    // Six zero eigenvalues, and the seventh is not. The comparison is on the
    // *eigenvalue*, not on the frequency: a square root turns 1e-11 into 3e-6 and
    // makes an exact zero look like a sloppy one.
    const double seventh = omega[6] * omega[6];
    for (int i = 0; i < 6; ++i)
        expectTrue("free-free eigenvalue " + std::to_string(i) + " is zero",
                   omega[static_cast<std::size_t>(i)] * omega[static_cast<std::size_t>(i)] <
                       1e-9 * seventh);
    expectTrue("and the seventh is not", seventh > 0.0);
    std::printf("     rigid eigenvalue / first elastic: %.2e, first elastic %.1f rad/s\n",
                omega[5] * omega[5] / seventh, omega[6]);

    // The mass went through the same transformation, so a rigid translation of the
    // reduced model weighs what the substructure weighs. This is the measurement
    // that decides the "reduce mass too" question: keeping M_bb alone would report
    // the interface's share and nothing else.
    const int n = reduced.size();
    double boundaryShare = 0;
    for (std::size_t b = 0; b < sub.boundaryCount(); ++b)
        if (sub.boundaryDof()[b] % 3 == 0) boundaryShare += sub.mass()[sub.boundaryDof()[b]];
    for (int axis = 0; axis < 3; ++axis) {
        std::vector<double> state(static_cast<std::size_t>(n), 0.0);
        for (std::size_t b = 0; b < sub.boundaryCount(); ++b)
            if (static_cast<int>(sub.boundaryDof()[b] % 3) == axis) state[b] = 1.0;
        double m = 0;
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                m += state[static_cast<std::size_t>(i)] *
                     reduced.mass[static_cast<std::size_t>(i) * n + j] *
                     state[static_cast<std::size_t>(j)];
        expectNear("rigid translation carries the whole mass (axis " + std::to_string(axis) + ")",
                   m, sub.totalMass(), 1e-6 * sub.totalMass());
        double k = 0;
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                k += state[static_cast<std::size_t>(i)] *
                     reduced.stiffness[static_cast<std::size_t>(i) * n + j] *
                     state[static_cast<std::size_t>(j)];
        expectTrue("and no strain energy (axis " + std::to_string(axis) + ")",
                   std::fabs(k) < 1e-12 * reduced.stiffness[0] * static_cast<double>(n));
    }
    expectTrue("keeping only the interface mass would have lost most of it",
               boundaryShare < 0.5 * sub.totalMass());
    std::printf("     the interface carries %.1f%% of the mass; the reduction carries 100%%\n",
                100.0 * boundaryShare / sub.totalMass());

    // A rigid *rotation* is the one lumping does not make trivial: the reduced
    // rotary inertia has to match the full lumped model's, which is a different
    // number from the mass and from the continuum's.
    {
        std::vector<double> state(static_cast<std::size_t>(n), 0.0);
        double want = 0;
        for (std::size_t node = 0; node < mesh.nodeCount(); ++node) {
            const double x = mesh.position[node * 3], y = mesh.position[node * 3 + 1];
            want += sub.mass()[node * 3] * (x * x + y * y);
        }
        for (std::size_t b = 0; b < sub.boundaryCount(); ++b) {
            const std::uint32_t d = sub.boundaryDof()[b];
            const std::size_t node = d / 3;
            if (d % 3 == 0) state[b] = -mesh.position[node * 3 + 1];
            if (d % 3 == 1) state[b] = mesh.position[node * 3];
        }
        double got = 0;
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                got += state[static_cast<std::size_t>(i)] *
                       reduced.mass[static_cast<std::size_t>(i) * n + j] *
                       state[static_cast<std::size_t>(j)];
        expectNear("rigid rotation carries the full model's rotary inertia", got, want,
                   1e-6 * want);
        expectTrue("and the rotary inertia is not the mass", std::fabs(want - sub.totalMass()) >
                                                                 0.1 * sub.totalMass());
    }
}

// --- 5. The reduced pair is T^T K T -----------------------------------------------

void testProjection() {
    std::printf("\n--- reduction: the reduced pair, formed the long way ---\n");
    solidshell::HexMesh mesh = testPlate();
    const StructuralMaterial steel = ah36Steel();
    Substructure sub(mesh, steel, endInterface(mesh));
    ReduceParams params;
    params.modes = 8;
    const Reduction reduced = craigBampton(sub, params);
    const int n = reduced.size();
    expectEqual("boundary DOF plus modes", n, 68);

    // T, column by column, out of the constraint and normal modes.
    std::vector<std::vector<double>> t(static_cast<std::size_t>(n));
    for (int c = 0; c < n; ++c) {
        std::vector<double> unit(static_cast<std::size_t>(n), 0.0);
        unit[static_cast<std::size_t>(c)] = 1.0;
        t[static_cast<std::size_t>(c)] = reduction::recover(sub, reduced, unit);
    }

    double worstK = 0, worstM = 0, scaleK = 0, scaleM = 0, worstCoupling = 0, worstModal = 0;
    std::vector<double> product;
    for (int c = 0; c < n; ++c) {
        sub.stiffnessTimes(t[static_cast<std::size_t>(c)], product);
        for (int d = 0; d < n; ++d) {
            double k = 0, m = 0;
            for (std::size_t i = 0; i < product.size(); ++i) {
                k += t[static_cast<std::size_t>(d)][i] * product[i];
                m += t[static_cast<std::size_t>(d)][i] * sub.mass()[i] *
                     t[static_cast<std::size_t>(c)][i];
            }
            const double storedK = reduced.stiffness[static_cast<std::size_t>(c) * n + d];
            const double storedM = reduced.mass[static_cast<std::size_t>(c) * n + d];
            worstK = std::max(worstK, std::fabs(k - storedK));
            worstM = std::max(worstM, std::fabs(m - storedM));
            scaleK = std::max(scaleK, std::fabs(k));
            scaleM = std::max(scaleM, std::fabs(m));
            // The two structural claims the assembly *asserts* rather than
            // computes: the boundary-modal stiffness block is identically zero and
            // the modal mass block is the identity.
            const bool boundaryC = c < reduced.boundary, boundaryD = d < reduced.boundary;
            if (boundaryC != boundaryD) worstCoupling = std::max(worstCoupling, std::fabs(k));
            if (!boundaryC && !boundaryD)
                worstModal = std::max(worstModal, std::fabs(m - (c == d ? 1.0 : 0.0)));
        }
    }
    expectTrue("the stored stiffness is T^T K T", worstK < 1e-9 * scaleK);
    expectTrue("the stored mass is T^T M T", worstM < 1e-9 * scaleM);
    expectTrue("the boundary-modal stiffness block really is zero",
               worstCoupling < 1e-9 * scaleK);
    expectTrue("the modal mass block really is the identity", worstModal < 1e-9);
    std::printf("     K: %.2e of %.2e, M: %.2e of %.2e, coupling %.2e, modal mass %.2e\n", worstK,
                scaleK, worstM, scaleM, worstCoupling, worstModal);

    // Symmetry, to the last bit. On its own this proves nothing -- the analytically
    // symmetric blocks are computed on one triangle and mirrored -- and it is here
    // because a caller factoring the reduced pair is entitled to rely on it, and
    // because the check above is what says they are the right matrices.
    long long asymmetric = 0;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            if (reduced.stiffness[static_cast<std::size_t>(i) * n + j] !=
                reduced.stiffness[static_cast<std::size_t>(j) * n + i])
                ++asymmetric;
            if (reduced.mass[static_cast<std::size_t>(i) * n + j] !=
                reduced.mass[static_cast<std::size_t>(j) * n + i])
                ++asymmetric;
        }
    expectEqual("the reduced matrices are symmetric bit for bit", asymmetric, 0);

    // The mass is positive definite and the stiffness -- with this interface, which
    // does restrain the rigid modes only when something is held -- is positive
    // definite once a cantilever's worth of it is held.
    {
        const Eigenpairs massSpectrum = reduction::symmetricEigen(reduced.mass, n);
        expectTrue("the reduced mass is positive definite", massSpectrum.value[0] > 0.0);
        std::printf("     reduced mass eigenvalues span %.3e to %.3e kg\n", massSpectrum.value[0],
                    massSpectrum.value.back());

        const std::vector<std::uint32_t> holding = heldAtOrigin(sub, mesh);
        const std::vector<double> held = [&] {
            std::vector<double> k;
            std::vector<std::size_t> keep;
            for (std::size_t d = 0; d < static_cast<std::size_t>(n); ++d)
                if (std::find(holding.begin(), holding.end(), static_cast<std::uint32_t>(d)) ==
                    holding.end())
                    keep.push_back(d);
            const std::size_t f = keep.size();
            k.assign(f * f, 0.0);
            for (std::size_t i = 0; i < f; ++i)
                for (std::size_t j = 0; j < f; ++j)
                    k[i * f + j] = reduced.stiffness[keep[i] * static_cast<std::size_t>(n) + keep[j]];
            const Eigenpairs s = reduction::symmetricEigen(k, static_cast<int>(f));
            return s.value;
        }();
        expectTrue("the held reduced stiffness is positive definite", held[0] > 0.0);
    }
}

// --- 6. Convergence from above, monotonically -------------------------------------

void testConvergence() {
    std::printf("\n--- reduction: frequencies come down from above ---\n");
    solidshell::HexMesh mesh = testPlate();
    const StructuralMaterial steel = ah36Steel();
    Substructure sub(mesh, steel, endInterface(mesh));
    const std::vector<std::uint32_t> held = heldAtOrigin(sub, mesh);
    std::vector<std::uint8_t> pinned(sub.dofCount(), 0);
    for (std::uint32_t b : held) pinned[sub.boundaryDof()[b]] = 1;
    const std::vector<double> exact = fullFrequencies(sub, pinned);
    expectTrue("the full model has frequencies to be compared against", exact.size() > 10);

    const int counts[] = {0, 1, 2, 4, 8, 16, 32};
    std::vector<std::vector<double>> got;
    for (int m : counts) {
        ReduceParams params;
        params.modes = m;
        got.push_back(reduction::reducedFrequencies(craigBampton(sub, params), held));
    }

    // Every reduced frequency is an upper bound on the full model's. The floor
    // below zero is the conditioning of the two eigensolves, not slack in the
    // claim: at 32 modes the two agree to a part in 1e9 and the sign of a
    // difference that size is noise.
    const int compare = 6;
    double worstBelow = 0;
    for (std::size_t c = 0; c < got.size(); ++c)
        for (int i = 0; i < compare; ++i) {
            const double rel = (got[c][static_cast<std::size_t>(i)] -
                                exact[static_cast<std::size_t>(i)]) /
                               exact[static_cast<std::size_t>(i)];
            worstBelow = std::min(worstBelow, rel);
        }
    expectTrue("no reduced frequency falls below the full model's", worstBelow > -1e-7);

    // And they fall monotonically as modes are added: the subspaces nest, so the
    // Rayleigh quotients can only improve.
    long long rises = 0;
    for (std::size_t c = 1; c < got.size(); ++c)
        for (int i = 0; i < compare; ++i)
            if (got[c][static_cast<std::size_t>(i)] >
                got[c - 1][static_cast<std::size_t>(i)] * (1.0 + 1e-9))
                ++rises;
    expectEqual("adding modes never raises a reduced frequency", rises, 0);

    // The vacuity guard: without it this passes on a reduction that was already
    // exact at zero modes and never had anything to converge.
    double firstWorst = 0, lastWorst = 0;
    for (int i = 0; i < compare; ++i) {
        firstWorst = std::max(firstWorst, (got.front()[static_cast<std::size_t>(i)] -
                                           exact[static_cast<std::size_t>(i)]) /
                                              exact[static_cast<std::size_t>(i)]);
        lastWorst = std::max(lastWorst, (got.back()[static_cast<std::size_t>(i)] -
                                         exact[static_cast<std::size_t>(i)]) /
                                            exact[static_cast<std::size_t>(i)]);
    }
    expectTrue("Guyan alone is wrong by tens of per cent, so there is something to converge",
               firstWorst > 0.3);
    expectTrue("and the modes converge it by five orders", lastWorst < 1e-4 * firstWorst);
    std::printf("     worst of the first %d: %.3f at 0 modes, %.3e at 32\n", compare, firstWorst,
                lastWorst);

    // Holding the whole interface leaves the reduced model's frequencies equal to
    // the fixed-interface frequencies themselves -- which are exact eigenvalues of
    // the full constrained problem, not approximations. A reduction that had
    // mis-normalised its modes would fail this and pass everything above.
    {
        ReduceParams params;
        params.modes = 5;
        const Reduction reduced = craigBampton(sub, params);
        std::vector<std::uint32_t> everything;
        for (int b = 0; b < reduced.boundary; ++b) everything.push_back(static_cast<std::uint32_t>(b));
        const std::vector<double> modal = reduction::reducedFrequencies(reduced, everything);
        std::vector<std::uint8_t> allPinned(sub.dofCount(), 0);
        for (std::size_t b = 0; b < sub.boundaryCount(); ++b) allPinned[sub.boundaryDof()[b]] = 1;
        const std::vector<double> fixedExact = fullFrequencies(sub, allPinned);
        double worst = 0;
        for (int j = 0; j < 5; ++j)
            worst = std::max(worst, std::fabs(modal[static_cast<std::size_t>(j)] -
                                              fixedExact[static_cast<std::size_t>(j)]) /
                                        fixedExact[static_cast<std::size_t>(j)]);
        expectTrue("with the interface held the reduced frequencies are exact", worst < 1e-9);
        std::printf("     fixed-interface frequencies reproduced to %.2e\n", worst);
    }
}

// --- 7. The mode count is verified, not hoped ------------------------------------

void testModeCount() {
    std::printf("\n--- reduction: the Sturm check, and choosing modes by frequency ---\n");
    solidshell::HexMesh mesh = testPlate();
    const StructuralMaterial steel = ah36Steel();
    Substructure sub(mesh, steel, endInterface(mesh));

    // The inertia count against the dense spectrum of the same interior problem.
    // Two entirely different computations of the same integer.
    std::vector<std::uint8_t> allPinned(sub.dofCount(), 0);
    for (std::size_t b = 0; b < sub.boundaryCount(); ++b) allPinned[sub.boundaryDof()[b]] = 1;
    const std::vector<double> interior = fullFrequencies(sub, allPinned);
    expectTrue("the fixed-interface problem has a spectrum", interior.size() > 20);

    long long disagreements = 0;
    for (int k : {1, 2, 3, 5, 8, 13, 21}) {
        // A shift half way between two consecutive eigenvalues, so the count is
        // unambiguous.
        const double lo = interior[static_cast<std::size_t>(k) - 1];
        const double hi = interior[static_cast<std::size_t>(k)];
        const double shift = 0.5 * (lo * lo + hi * hi);
        bool exactCount = false;
        if (sub.eigenvaluesBelow(shift, &exactCount) != k) ++disagreements;
        if (!exactCount) ++disagreements;
    }
    expectEqual("the inertia count agrees with the dense spectrum at every shift", disagreements, 0);
    expectEqual("nothing lies below zero", sub.eigenvaluesBelow(0.0, nullptr), 0);
    expectEqual("everything lies below infinity",
                sub.eigenvaluesBelow(1e30, nullptr),
                static_cast<long long>(sub.interiorCount()));

    // Subspace iteration finds the lowest modes, and says so -- at every mode count,
    // not at one. The shift the check uses sits half way between the last mode kept
    // and the first discarded; putting it *on* the last mode kept gives a count that
    // may come out either side by a rounding, so one sample would not settle it.
    for (int m : {2, 5, 12, 30}) {
        ReduceParams params;
        params.modes = m;
        const Reduction reduced = craigBampton(sub, params);
        expectTrue("the mode count is verified at " + std::to_string(m) + " modes",
                   reduced.modesVerified);
    }
    {
        ReduceParams params;
        params.modes = 12;
        const Reduction reduced = craigBampton(sub, params);
        expectTrue("the mode count is verified by the Sturm check", reduced.modesVerified);
        double worst = 0;
        for (int j = 0; j < 12; ++j)
            worst = std::max(worst, std::fabs(reduced.frequency[static_cast<std::size_t>(j)] -
                                              interior[static_cast<std::size_t>(j)]) /
                                        interior[static_cast<std::size_t>(j)]);
        expectTrue("and subspace iteration agrees with the dense solve", worst < 1e-9);
        std::printf("     twelve fixed-interface modes to %.2e of a dense solve\n", worst);
    }

    // A frequency cutoff is answered by the inertia count, so the number of modes
    // kept is exactly the number below it -- not one more, not one less.
    {
        const double cutoff = 0.5 * (interior[6] + interior[7]);
        ReduceParams params;
        params.modes = -1;
        params.cutoffFrequency = cutoff;
        const Reduction reduced = craigBampton(sub, params);
        expectEqual("a frequency cutoff keeps exactly the modes below it", reduced.modes, 7);
        expectTrue("every mode kept is below the cutoff",
                   reduced.frequency.back() < cutoff);
        expectTrue("and the first one discarded is above it",
                   interior[7] > cutoff);
    }

    // Zero modes still reports where the first one was, because "no modes were
    // needed" is only defensible as a number.
    {
        ReduceParams params;
        params.modes = 0;
        const Reduction reduced = craigBampton(sub, params);
        expectEqual("zero modes kept", reduced.modes, 0);
        expectNear("and the first fixed-interface frequency is still reported",
                   reduced.firstFixedFrequency, interior[0], 1e-6 * interior[0]);
    }

    // Asking for more modes than the interior has is a request that cannot be met,
    // and the answer is the whole interior rather than a walk off the end of the
    // spectrum. This is the path that crashed the first version of the file.
    {
        ReduceParams params;
        params.modes = static_cast<int>(sub.interiorCount()) + 10;
        // `maxModes` would otherwise bite first, and this is a test of the other clamp.
        params.maxModes = 10000;
        const Reduction reduced = craigBampton(sub, params);
        expectEqualCount("more modes than the interior has is clamped to the interior",
                         static_cast<std::size_t>(reduced.modes), sub.interiorCount());
        expectEqualCount("and the frequency list matches",
                         reduced.frequency.size(), sub.interiorCount());
        expectTrue("and the reduction is still a reduction",
                   reduced.frequency.front() > 0 &&
                       reduced.frequency.back() >= reduced.frequency.front());
        // Clamping a request that cannot be met is not a failure and must not be
        // reported as one. Without the clamp the eigensolver clamps instead, and the
        // reduction reports modes it never asked for as modes it could not compute.
        expectTrue("and nothing is reported as having gone wrong", reduced.problems.empty());
    }

    // An eigensolve that did not converge must not have its mode count reported as
    // verified. One iteration of subspace iteration is nowhere near the answer, and
    // the Sturm check is precisely the thing that is supposed to notice.
    {
        ReduceParams params;
        params.modes = 8;
        params.maxIterations = 1;
        const Reduction reduced = craigBampton(sub, params);
        expectTrue("an unconverged eigensolve is not reported as verified",
                   !reduced.modesVerified);
        expectTrue("and says why", !reduced.problems.empty());
    }

    // An eigensolve that ran no iterations at all has no modes to give, and the
    // reduction has to fall back to Guyan and say so rather than build a modal block
    // out of whatever the uninitialised values were. That was a crash once.
    {
        ReduceParams params;
        params.modes = 8;
        params.maxIterations = 0;
        const Reduction reduced = craigBampton(sub, params);
        expectEqual("an eigensolve given no iterations yields no modes", reduced.modes, 0);
        expectTrue("and the reduction says so rather than inventing them",
                   !reduced.problems.empty());
        expectTrue("and what is left is a usable Guyan condensation",
                   reduced.size() == reduced.boundary && reduced.stiffness.size() ==
                       static_cast<std::size_t>(reduced.boundary) *
                           static_cast<std::size_t>(reduced.boundary));
    }

    // The extra vectors in the subspace block are what make it converge: the rate
    // is `(lambda_j / lambda_q)^2`, so a block exactly as wide as the modes wanted
    // has a ratio of one for the last of them. The iteration count is deterministic,
    // so it can be asserted where a wall clock could not.
    {
        const Eigenpairs modes = sub.fixedInterfaceModes(12);
        expectTrue("subspace iteration converges", modes.converged);
        expectTrue("in a handful of iterations", modes.iterations <= 12);
        std::printf("     twelve modes in %d subspace iterations\n", modes.iterations);
    }
}

// --- 8. Where the linear model stops ---------------------------------------------

void testValidity() {
    std::printf("\n--- reduction: what a reduced region cannot do ---\n");
    // The **rotated** plate, and the load below is membrane plus bending. Both are
    // there because mutation testing killed neither the shear terms of the von
    // Mises formula nor the sweep over the eight Gauss points on the axis-aligned
    // bending case that was here first, and both survivors were the same mistake:
    // the peak of a bent, axis-aligned plate is very nearly a uniaxial stress at a
    // free surface, so `sqrt(0.5 sum (s_i - s_j)^2)` agrees with the real thing and
    // the corner Gauss point is as stressed as any other.
    //
    // Rotating the plate 45 degrees fixes the first: the stress is reported in the
    // **global** frame, so a stress along the plate is equal parts sigma_xx,
    // sigma_yy and sigma_xy and dropping the shear halves it. Adding a membrane
    // load of the same size as the bending fixes the second: the two add on one
    // face and cancel on the other, so the zeta = -1 Gauss points -- which is what
    // index 0 is -- carry almost nothing.
    solidshell::HexMesh mesh = rotatedPlate();
    const StructuralMaterial steel = ah36Steel();
    Substructure sub(mesh, steel, rotatedEndInterface(mesh));
    const std::vector<std::uint32_t> held = heldAtRotatedOrigin(sub, mesh);

    solidshell::HexMesh full = rotatedPlate();
    for (std::uint32_t b : held) {
        const std::uint32_t d = sub.boundaryDof()[b];
        full.pin(d / 3, static_cast<int>(d % 3), 0.0);
    }
    // The load is an **eccentric in-plane pull**: along the plate's own axis, on
    // the outer face only. A force at t/2 off the mid-surface puts 4 P/A on the
    // face it acts on and -2 P/A on the other, so the zeta = -1 integration points
    // -- which is what Gauss index 0 is -- carry less than half the peak, and the
    // sweep over the eight of them is under test rather than assumed. A transverse
    // load would not do it: bending is symmetric about the mid-surface and both
    // faces peak together.
    std::vector<double> load(sub.dofCount(), 0.0);
    for (std::size_t p = 0; p < sub.interiorCount(); ++p) {
        const std::uint32_t d = sub.interiorDof()[p];
        if (d % 3 == 2 || mesh.position[(d / 3) * 3 + 2] <= 0) continue;
        load[d] = kEccentricLoad * ((d % 3 == 0) ? kAxis.x : kAxis.y);
    }
    std::vector<double> reference;
    std::string problem;
    solidshell::solveStatic(full, steel, solidshell::Formulation::SolidShell, load, reference,
                            &problem);
    const StressPeak peak = peakVonMises(mesh, steel, reference);
    const double referencePeak = peak.mises;
    expectTrue("the reference load raises a real stress", referencePeak > 1e6);
    expectTrue("and the peak carries real shear, so the von Mises formula is under test",
               peak.withoutShear < 0.85 * referencePeak);
    expectTrue("and the peak is not at the first Gauss point, so the sweep is under test",
               peak.firstGaussOnly < 0.85 * referencePeak);
    std::printf("     eccentric pull: peak %.2f MPa, %.2f without the shear terms, %.2f at "
                "Gauss point 0\n",
                referencePeak * 1e-6, peak.withoutShear * 1e-6, peak.firstGaussOnly * 1e-6);

    // The recovered stress converges to the full model's, and it converges from
    // **below**: a truncated basis cannot represent a concentration, so the warning
    // a caller gets is late rather than early. That is the direction that matters
    // and it is asserted rather than described.
    double guyanPeak = 0, convergedPeak = 0;
    for (int m : {0, 64}) {
        ReduceParams params;
        params.modes = m;
        const Reduction reduced = craigBampton(sub, params);
        std::vector<double> state;
        reduction::staticSolve(reduced, reduction::reduceLoad(sub, reduced, load), held, state,
                               &problem);
        const reduction::Validity validity = reduction::checkValidity(sub, reduced, state);
        (m == 0 ? guyanPeak : convergedPeak) = validity.peakVonMises;
        expectTrue("the reduced model reports itself linear under a working load", validity.linear);
    }
    expectTrue("Guyan alone under-states the peak stress", guyanPeak < 0.95 * referencePeak);
    expectTrue("enough modes recover it", std::fabs(convergedPeak - referencePeak) <
                                              0.02 * referencePeak);
    std::printf("     peak von Mises: %.2f MPa at 0 modes, %.2f at 64, %.2f in the full model\n",
                guyanPeak * 1e-6, convergedPeak * 1e-6, referencePeak * 1e-6);

    // **It cannot yield**, and the assertion is that it is exactly linear: a
    // thousandfold load gives a thousandfold stress, straight through the yield
    // strength without a word. `linear` going false is the only thing that happens,
    // and it is the caller's cue to promote the region to Tier 2.
    {
        ReduceParams params;
        params.modes = 8;
        const Reduction reduced = craigBampton(sub, params);
        std::vector<double> small, large;
        std::vector<double> heavy = load;
        for (double& f : heavy) f *= 1000.0;
        reduction::staticSolve(reduced, reduction::reduceLoad(sub, reduced, load), held, small,
                               &problem);
        reduction::staticSolve(reduced, reduction::reduceLoad(sub, reduced, heavy), held, large,
                               &problem);
        const reduction::Validity light = reduction::checkValidity(sub, reduced, small);
        const reduction::Validity heavyResult = reduction::checkValidity(sub, reduced, large);
        expectNear("a thousandfold load gives a thousandfold stress",
                   heavyResult.peakVonMises / light.peakVonMises, 1000.0, 1e-6);
        expectTrue("past yield it reports itself invalid rather than yielding",
                   !heavyResult.linear && heavyResult.utilisation > 1.0);
        expectTrue("and the light case really was inside yield", light.linear);
        std::printf("     a load %.0fx past yield is reported at utilisation %.1f, still linear "
                    "in every number -- the region has to be promoted, not asked again\n",
                    heavyResult.utilisation, heavyResult.utilisation);
    }
}

// --- 9. A real patch of the ferry's own plating ----------------------------------

void testFerryPatch() {
    std::printf("\n--- reduction: the ferry's own side shell ---\n");
    const Ship ship = game::buildFerry();
    const StructuralMesh structure = makeStructuralMesh(ship.hull, ferryScantlings(), nullptr);
    zone::MeshParams params;
    params.radius = 2.5;
    params.stiffeners = zone::Stiffeners::Ignored;
    const zone::Patch patch = zone::buildPatch(structure, Vec3{10.0, 9.0, 4.0}, params);
    expectTrue("the ferry patch meshed", !patch.empty());

    // The patch arrives with its perimeter clamped, and that perimeter is exactly
    // the set of nodes the plating outside it holds -- so it is the interface.
    solidshell::HexMesh mesh = patch.mesh;
    const std::vector<std::uint32_t> interface = reduction::nodesPinned(mesh);
    mesh.fixed.assign(mesh.nodeCount() * 3, 0);
    Substructure sub(mesh, patch.material, interface);
    expectTrue("the ferry patch reduces", sub.ready());
    expectTrue("its perimeter is a genuine interface, not the whole patch",
               sub.interiorCount() > sub.boundaryCount());

    // The lumped mass is the mass of the *elements*, so it is checked against the
    // quadrature that measures their volume rather than against `Patch::mass`.
    // Those are two different quantities and it is worth saying which: `Patch::mass`
    // is mid-surface area times thickness, and on curved plating extruded along
    // averaged nodal normals the hexahedra do not have exactly that volume. Row-sum
    // lumping conserves whatever volume the elements actually have, which is the
    // property being asserted; the difference from the plate-area figure is
    // reported, not asserted away.
    {
        // Two properties, and the second is the one that says the lumping is a
        // *row sum* rather than volume over eight. Row-sum lumping is
        // `m_a = rho integral N_a dV`, and since `sum_a N_a x_a` is the isoparametric
        // map, `sum_a m_a x_a = rho integral x dV` **exactly**: the lumped masses
        // reproduce the meshed volume's own centre of gravity. Splitting a distorted
        // element's mass evenly does not, and on the ferry's curved, tapered plating
        // the difference is measurable -- the guard below insists it is, because on
        // a uniform brick mesh the two lumpings are identical and the check would be
        // vacuous.
        double volume = 0, moment[3] = {0, 0, 0}, evenMoment[3] = {0, 0, 0};
        for (std::size_t e = 0; e < mesh.elementCount(); ++e) {
            double nodePos[24], elementVolume = 0, elementMoment[3];
            mesh.gather(e, mesh.position, nodePos);
            hexVolumeMoment(nodePos, elementVolume, elementMoment);
            volume += elementVolume;
            for (int r = 0; r < 3; ++r) {
                moment[r] += elementMoment[r];
                double corners = 0;
                for (int a = 0; a < 8; ++a) corners += nodePos[a * 3 + r];
                evenMoment[r] += elementVolume * corners / 8.0;
            }
        }
        expectNear("lumped mass is the density times the meshed volume", sub.totalMass(),
                   patch.material.density * volume, 1e-8 * sub.totalMass());
        expectTrue("and within a per cent of area times thickness",
                   std::fabs(sub.totalMass() - patch.mass) < 0.01 * patch.mass);

        double centre[3] = {0, 0, 0};
        for (std::size_t n = 0; n < mesh.nodeCount(); ++n)
            for (int r = 0; r < 3; ++r)
                centre[r] += sub.mass()[n * 3] * mesh.position[n * 3 + static_cast<std::size_t>(r)];
        double worst = 0, evenDifference = 0, extent = 0;
        for (int r = 0; r < 3; ++r) {
            worst = std::max(worst, std::fabs(centre[r] / sub.totalMass() - moment[r] / volume));
            evenDifference =
                std::max(evenDifference, std::fabs(evenMoment[r] / volume - moment[r] / volume));
            extent = std::max(extent, std::fabs(moment[r] / volume));
        }
        expectTrue("the lumped masses sit at the meshed volume's centre of gravity",
                   worst < 1e-12 * std::max(extent, 1.0));
        expectTrue("and an even split of each element's mass would not",
                   evenDifference > 1e4 * worst);
        std::printf("     lumped %.2f kg against area x thickness %.2f kg (%.3f%%); centre of "
                    "gravity to %.2e m, where an even split is out by %.2e m\n",
                    sub.totalMass(), patch.mass, 100.0 * (sub.totalMass() - patch.mass) / patch.mass,
                    worst, evenDifference);
    }

    // The interior is renumbered before it is banded, and the band is what the
    // factorisation costs and stores. Measured on this patch: the flood fill's own
    // numbering gives 173, Cuthill-McKee started at a minimum-degree node gives 89,
    // and started at a pseudo-peripheral one gives 71. The assertion sits below the
    // middle of those, so dropping either half of the ordering is a failure rather
    // than a slowdown nobody sees.
    expectTrue("the interior renumbering earns its place", sub.halfBandwidth() <= 80);
    std::printf("     half-bandwidth %zu over %zu interior DOF (%.2f MB of banded store)\n",
                sub.halfBandwidth(), sub.interiorCount(),
                8.0 * static_cast<double>(sub.interiorCount()) *
                    static_cast<double>(sub.halfBandwidth() + 1) / 1e6);

    ReduceParams reduceParams;
    reduceParams.modes = 4;
    const Reduction reduced = craigBampton(sub, reduceParams);
    expectTrue("the mode count is verified on real geometry", reduced.modesVerified);
    expectTrue("its frequencies are real and ascending",
               reduced.frequency.front() > 0 && reduced.frequency.back() >= reduced.frequency.front());

    // The identity survives curved, unevenly sized plating: static condensation is
    // still exact at the interface.
    solidshell::HexMesh full = patch.mesh;
    std::vector<double> load(sub.dofCount(), 0.0);
    std::vector<std::uint32_t> held;
    double xMid = 0;
    for (std::size_t b = 0; b < sub.boundaryCount(); ++b)
        xMid += mesh.position[(sub.boundaryDof()[b] / 3) * 3];
    xMid /= static_cast<double>(sub.boundaryCount());
    for (std::size_t b = 0; b < sub.boundaryCount(); ++b) {
        const std::uint32_t d = sub.boundaryDof()[b];
        if (mesh.position[(d / 3) * 3] < xMid) {
            held.push_back(static_cast<std::uint32_t>(b));
        } else {
            full.fixed[d] = 0;
            if (d % 3 == 1) load[d] = 50.0;
        }
    }
    for (std::size_t d = 0; d < sub.dofCount(); ++d) full.fixed[d] = 0;
    for (std::uint32_t b : held) full.fixed[sub.boundaryDof()[b]] = 1;
    std::vector<double> uFull;
    std::string problem;
    expectTrue("the full patch solves",
               solidshell::solveStatic(full, patch.material, solidshell::Formulation::SolidShell,
                                       load, uFull, &problem));
    double peak = 0;
    for (double e : uFull) peak = std::max(peak, std::fabs(e));
    expectTrue("the load moves the patch", peak > 1e-5);

    ReduceParams guyanParams;
    guyanParams.modes = 0;
    const Reduction guyan = craigBampton(sub, guyanParams);
    std::vector<double> state;
    reduction::staticSolve(guyan, reduction::reduceLoad(sub, guyan, load), held, state, &problem);
    const std::vector<double> u = reduction::recover(sub, guyan, state);
    double worst = 0;
    for (std::size_t b = 0; b < sub.boundaryCount(); ++b)
        worst = std::max(worst, std::fabs(u[sub.boundaryDof()[b]] - uFull[sub.boundaryDof()[b]]));
    expectTrue("static condensation is exact on the ferry's own plating too", worst < 1e-6 * peak);
    std::printf("     %zu elements, %zu interface DOF, %zu interior, band %zu\n",
                patch.elementCount(), sub.boundaryCount(), sub.interiorCount(),
                sub.halfBandwidth());
    std::printf("     interface response %.2e m of a peak %.3e m; reduce %.3f s\n", worst, peak,
                reduced.reduceSeconds);

    // Cost, printed rather than asserted -- `test_plasticity.cpp` records what a
    // tight timing assertion costs on a shared machine. The comparison that matters
    // is against the same patch at Tier 2.
    std::printf("     Tier-2 on this patch: %.0f core-seconds per simulated second\n",
                zone::estimatedCost(patch));
}

// --- Synthesis: two components joined at an interface ----------------------------
//
// The load-bearing test of the whole tier. Reduction on its own is one component,
// and one component proves nothing about *synthesis*: the claim is that a
// structure can be cut into pieces, each reduced independently, and reassembled
// into a model of the original. So the reference is the original -- the same plate
// meshed in one piece, whose free-free spectrum is computed the long way and owes
// nothing to any of the machinery under test.
//
// The negative control matters as much. Two halves that are *not* coupled are
// still two perfectly good reduced models sitting in the same matrix, and they
// produce a full set of plausible frequencies. What gives them away is that an
// uncoupled pair has **twelve** rigid body modes rather than six, because each
// half floats free of the other. Without that check this test would pass on an
// assembly that never joined anything.
void testTwoComponentsAssembleIntoTheWhole() {
    std::printf("\n--- reduction: two components joined at an interface ---\n");

    const double L = 1.0, W = 0.2, T = 0.01;
    const int NX = 8, NY = 2, NZ = 1;
    const StructuralMaterial steel = ah36Steel();

    // The same plate three ways: whole, and split down the middle by ELEMENT. The
    // node spacing matches by construction -- L*i/NX on the whole, and
    // (L/2)*i/(NX/2) on each half is the same number -- so the halves have the
    // interface nodes in common and nothing is interpolated.
    solidshell::HexMesh whole = solidshell::makePlateMesh(L, W, T, NX, NY, NZ);
    solidshell::HexMesh left = solidshell::makePlateMesh(0.5 * L, W, T, NX / 2, NY, NZ);
    solidshell::HexMesh right = solidshell::makePlateMesh(0.5 * L, W, T, NX / 2, NY, NZ);
    for (std::size_t i = 0; i + 2 < right.position.size(); i += 3) right.position[i] += 0.5 * L;

    const std::vector<reduction::Plane> split{{{0.5 * L, 0, 0}, {1, 0, 0}}};
    const std::vector<std::uint32_t> ifLeft = reduction::nodesNearPlanes(left, split, 1e-9);
    const std::vector<std::uint32_t> ifRight = reduction::nodesNearPlanes(right, split, 1e-9);
    expectEqualCount("the split plane carries the interface nodes", ifLeft.size(),
                     static_cast<std::size_t>((NY + 1) * (NZ + 1)));
    expectEqualCount("and the other half has the same ones", ifRight.size(), ifLeft.size());

    reduction::Substructure sa(left, steel, ifLeft);
    reduction::Substructure sb(right, steel, ifRight);
    expectTrue("the left component is ready", sa.ready());
    expectTrue("the right component is ready", sb.ready());

    // Splitting by element must not create or destroy mass. A mesh split by NODE
    // would double the interface mass, and every frequency below would come out
    // low while still looking entirely reasonable.
    reduction::Substructure sWhole(whole, steel, reduction::nodesNearPlanes(whole, split, 1e-9));
    expectTrue("the whole plate is ready", sWhole.ready());
    expectNear("the two halves weigh what the whole plate weighs",
               sa.totalMass() + sb.totalMass(), sWhole.totalMass(), 1e-9 * sWhole.totalMass());

    const reduction::InterfaceMap map = reduction::matchBoundaries(sa, sb);
    expectTrue("matching raises nothing", map.problems.empty());
    expectEqualCount("every boundary DOF of one half is shared with the other", map.shared,
                     3 * ifLeft.size());
    expectTrue("and the matched nodes are actually coincident", map.worstGap < 1e-12);

    // The whole plate's own free-free spectrum, formed densely from the operator
    // and solved once. This is the reference and it uses none of the reduction.
    const std::size_t nFull = sWhole.dofCount();
    std::vector<double> kFull(nFull * nFull, 0.0), mFull(nFull * nFull, 0.0);
    {
        std::vector<double> e(nFull, 0.0), col(nFull, 0.0);
        for (std::size_t j = 0; j < nFull; ++j) {
            std::fill(e.begin(), e.end(), 0.0);
            e[j] = 1.0;
            sWhole.stiffnessTimes(e, col);
            for (std::size_t i = 0; i < nFull; ++i) kFull[i * nFull + j] = col[i];
        }
        const std::vector<double>& lumped = sWhole.mass();
        for (std::size_t i = 0; i < nFull; ++i) mFull[i * nFull + i] = lumped[i];
    }
    const reduction::Eigenpairs exact =
        reduction::generalisedEigen(kFull, mFull, static_cast<int>(nFull));
    expectTrue("the reference spectrum converged", exact.converged);

    // A free solid has exactly six rigid body modes, so the seventh eigenvalue is
    // the first elastic one. That count is known physics rather than something to
    // discover, so what gets asserted is the *separation* -- and that separation is
    // then the scale everything below is measured against, instead of a magic
    // number in rad/s that would silently stop meaning anything if the plate
    // changed.
    //
    // The first version of this test thresholded rigid modes at 1e-3 rad/s and
    // found three of the six. The three translations come out at exactly zero and
    // the three rotations at 8e-4 to 9e-3 -- near zero, but not to machine
    // precision, and a fixed cutoff landed in the middle of them. This is the
    // mistake the reduction work had already recorded once: compare eigenvalues,
    // or scale by something derived, because a square root flatters a zero.
    expectTrue("the reference spectrum has at least seven modes", exact.value.size() > 6);
    const double firstElastic = std::sqrt(exact.value[6]);
    expectTrue("the sixth mode is rigid and the seventh is not, by orders of magnitude",
               std::sqrt(exact.value[5]) < 1e-4 * firstElastic);
    const double rigidCutoff = 1e-3 * firstElastic;

    std::vector<double> exactElastic;
    for (std::size_t i = 6; i < exact.value.size(); ++i)
        exactElastic.push_back(std::sqrt(exact.value[i]));

    // Now the same structure as two reduced components, at rising mode counts. The
    // assembled answer must approach the reference from above and keep improving.
    double previous = 0;
    for (int modes : {0, 4, 12}) {
        reduction::ReduceParams params;
        params.modes = modes;
        params.cutoffFrequency = 0;
        const reduction::Reduction ra = reduction::craigBampton(sa, params);
        const reduction::Reduction rb = reduction::craigBampton(sb, params);
        const reduction::Assembly asm2 = reduction::assemble(ra, rb, map);
        expectTrue("the assembly raises nothing", asm2.problems.empty());
        expectEqualCount("the assembly counts the shared boundary once",
                         static_cast<std::size_t>(asm2.boundary),
                         static_cast<std::size_t>(ra.boundary + rb.boundary) - map.shared);

        const std::vector<double> omega = reduction::assembledFrequencies(asm2, {});
        int rigid = 0;
        for (double w : omega)
            if (w < rigidCutoff) ++rigid;
        expectEqual("a free assembly of two components has six rigid body modes, not twelve",
                    static_cast<long long>(rigid), 6);

        double worst = 0;
        const std::size_t compare = 4;
        for (std::size_t i = 0; i < compare && i < exactElastic.size(); ++i) {
            const double got = omega[static_cast<std::size_t>(rigid) + i];
            worst = std::max(worst, std::fabs(got - exactElastic[i]) / exactElastic[i]);
            expectTrue("an assembled frequency is an upper bound on the true one",
                       got > exactElastic[i] * (1.0 - 1e-6));
        }
        std::printf("     %2d modes per half: worst of the first %zu elastic modes %.3e\n", modes,
                    compare, worst);
        if (modes > 0)
            expectTrue("adding modes to each component improves the assembled spectrum",
                       worst < previous);
        previous = worst;
    }
    expectTrue("and twelve modes a side gets the assembled spectrum to four figures",
               previous < 1e-4);

    // The negative control. Assemble the same two components with nothing matched:
    // the result is still a well-formed reduced model of the right size, and it is
    // wrong in the one way that is hard to see -- the halves are no longer joined.
    reduction::ReduceParams params;
    params.modes = 12;
    params.cutoffFrequency = 0;
    const reduction::Reduction ra = reduction::craigBampton(sa, params);
    const reduction::Reduction rb = reduction::craigBampton(sb, params);
    reduction::InterfaceMap none;
    none.aToB.assign(static_cast<std::size_t>(ra.boundary), -1);
    const reduction::Assembly loose = reduction::assemble(ra, rb, none);
    int looseRigid = 0;
    for (double w : reduction::assembledFrequencies(loose, {}))
        if (w < rigidCutoff) ++looseRigid;
    expectEqual("two components that were never joined float free of each other, twelve modes",
                static_cast<long long>(looseRigid), 12);
    expectTrue("and an unmatched interface is reported rather than assembled silently",
               !reduction::matchBoundaries(sa, sWhole, 1e-12).problems.empty() ||
                   reduction::matchBoundaries(sa, sb, 1e-9).shared > 0);
}

}  // namespace

void runReductionTests() {
    std::printf("\n=== Tier-1 reduction (Craig-Bampton) ===\n");
    testDenseEigensolver();
    testSubstructure();
    testGuyanIsExact();
    testRigidBody();
    testProjection();
    testConvergence();
    testModeCount();
    testValidity();
    testFerryPatch();
    testTwoComponentsAssembleIntoTheWhole();
}
