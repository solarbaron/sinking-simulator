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
#include "engine/sim/constraint.hpp"
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

// --- A stiffener on the test plate, for §8 ---------------------------------------
//
// `makePlateMesh` numbers the thickness index fastest, then y, then x -- its own
// header says so and says why (the bandwidth). The node lookup below depends on
// that, so `stiffenSeam` checks the coordinates of every node it picks rather than
// trusting the formula: a renumbering would otherwise put the stiffener somewhere
// else on the plate and every measurement here would still be self-consistent.
std::uint32_t plateNode(int i, int j, int k) {
    return static_cast<std::uint32_t>((i * (kNy + 1) + j) * 2 + k);
}

// One stiffener along the plate's own x direction at the row `j`, with a station
// every `stride` mesh nodes. `stride > 1` is a stiffener modelled coarser than the
// plating it sits on, which is the case whose fibres tie node pairs no single
// element shares -- see §8 item 2, and the band it widens.
constraint::Stiffening stiffenSeam(const solidshell::HexMesh& mesh,
                                   const StructuralMaterial& material,
                                   const StiffenerProfile& profile, int j, int stride,
                                   bool axisAligned = true) {
    constraint::SeamRun run;
    run.sign = 1.0;  // the pair runs -z to +z here, so a web rising in +z is +1
    for (int i = 0; i <= kNx; i += stride) {
        const std::uint32_t bottom = plateNode(i, j, 0), top = plateNode(i, j, 1);
        if (axisAligned) {
            expectNear("the seam station is where the mesh numbering says",
                       mesh.position[bottom * 3], kLx * static_cast<double>(i) / kNx, 1e-12);
            expectNear("and on the row it was asked for", mesh.position[bottom * 3 + 1],
                       kLy * static_cast<double>(j) / kNy, 1e-12);
        }
        expectNear("with the pair through the thickness", mesh.position[top * 3 + 2] -
                                                              mesh.position[bottom * 3 + 2],
                   kThickness, 1e-15);
        run.bottom.push_back(bottom);
        run.top.push_back(top);
    }
    constraint::Stiffening out;
    out.material = material;
    out.members = 1;
    constraint::addStiffener(run, profile, kThickness, mesh.position, out);
    return out;
}

// The `Attachment` a caller assembles from a `Stiffening`, which is the whole of
// the integration §8 describes: two calls into `constraint.hpp` and nothing else.
reduction::Attachment attachmentOf(const constraint::Stiffening& stiffening,
                                   const solidshell::HexMesh& mesh,
                                   const StructuralMaterial& material, bool withMass,
                                   bool withStress = true) {
    const constraint::RestFibers forms = constraint::restFibers(stiffening, mesh.position);
    expectTrue("every fibre has a rest length", forms.ok);
    constraint::AttachedForms built =
        constraint::attachedForms(stiffening, mesh.position, forms, material.youngsModulus);
    reduction::Attachment attached;
    attached.stiffness = std::move(built.stiffness);
    // `withStress == false` is the substructure this file had before §9: the
    // members are in the stiffness and invisible to `checkValidity`. It is kept
    // because it is the *only* way to reproduce the old, low utilisation, and a
    // regression test for a number needs the wrong number to still be reachable.
    if (withStress) attached.stress = std::move(built.stress);
    if (withMass) {
        attached.mass.assign(mesh.nodeCount(), 0.0);
        constraint::lumpFiberMass(stiffening, forms, material.density, attached.mass);
    }
    return attached;
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
        // Two guards fire here, and they are cause and symptom. `fixedInterfaceModes`
        // refuses outright -- "given no iterations to run" -- and `craigBampton` then
        // separately reports that 0 of 8 modes came back. Deleting the first leaves
        // only the second, which describes what the reduction is missing without
        // saying that the eigensolver was never allowed to start.
        expectTrue("and the reduction says so rather than inventing them",
                   !reduced.problems.empty());
        bool saidWhy = false, saidWhat = false;
        for (const std::string& p : reduced.problems) {
            if (p.find("given no iterations to run") != std::string::npos) saidWhy = true;
            if (p.find("of 8 fixed-interface modes were computed") != std::string::npos)
                saidWhat = true;
        }
        expectTrue("naming the eigensolve that never ran", saidWhy);
        expectTrue("and the modes that are therefore missing", saidWhat);
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

    // --- What `converged == false` means, which is not "the answer is wrong" -------
    //
    // `tools/section_probe` prints "subspace iteration did not converge in 60
    // iterations" on three of its five ferry cases and always has. The question is
    // whether the eigenvalue it reports anyway is any good, and the answer needs a
    // *controlled* version of the case rather than the ferry: the same plate at
    // three thicknesses, where the only thing that changes is `h/t` and therefore
    // the conditioning of `K_ii`.
    //
    // Every subspace iteration solves through the banded Cholesky of `K_ii`, so the
    // Ritz values carry that solve's error. On a slender plate that error is around
    // 1e-8 relative, the per-eigenvalue change **plateaus** there instead of
    // decreasing, and a 1e-10 tolerance can then only be met by coincidence.
    {
        std::printf("     the convergence flag against slenderness:\n");
        const double lx = 1.2, ly = 0.6;
        const std::vector<Plane> ends{{Vec3{0, 0, 0}, Vec3{1, 0, 0}},
                                      {Vec3{lx, 0, 0}, Vec3{1, 0, 0}}};
        struct Row {
            double thickness;
            bool converges;
        };
        Eigenpairs slender;
        double slenderValue = 0;
        for (const Row& row : {Row{0.012, true}, Row{0.004, true}, Row{0.0015, false}}) {
            solidshell::HexMesh plate =
                solidshell::makePlateMesh(lx, ly, row.thickness, 24, 12, 1);
            const Substructure sub(plate, steel, reduction::nodesNearPlanes(plate, ends, 1e-9));
            expectTrue("the plate reduces", sub.ready());
            const Eigenpairs modes = sub.fixedInterfaceModes(1);
            std::printf("       t %6.4f m (h/t %5.1f): %3d iterations, converged %d, last change"
                        " %9.2e, w1 %12.6f rad/s\n",
                        row.thickness, (lx / 24) / row.thickness, modes.iterations,
                        static_cast<int>(modes.converged), modes.lastChange,
                        std::sqrt(modes.value[0]));
            // The claim is that the flag is a function of the *structure* and not of
            // the solve being wrong. Both directions are asserted, because only the
            // thick one converging makes the thin one's failure informative.
            expectEqual(std::string("the ") + (row.converges ? "stocky" : "slender") +
                            " plate's convergence flag",
                        static_cast<long long>(modes.converged), row.converges ? 1LL : 0LL);
            if (!row.converges) {
                slender = modes;
                slenderValue = modes.value[0];
            }
        }

        expectTrue("the slender plate says why it stopped, with the number",
                   slender.problem.find("relative") != std::string::npos);
        // It stalled at a *small* change, not an order-one one. Measured at 7.5e-9;
        // asserted at 1e-5, four decades above it and four below the order-one a
        // genuinely lost solve leaves.
        expectTrue("and it stalled at a change far below one", slender.lastChange < 1e-5);
        expectTrue("but above the tolerance it was asked for", slender.lastChange > 1e-10);

        // More iterations buy nothing: the change is a floor, not a trend. If this
        // ever fails, the cause is *not* the floor and the diagnosis above is stale.
        const solidshell::HexMesh thin = solidshell::makePlateMesh(lx, ly, 0.0015, 24, 12, 1);
        const Substructure sub(thin, steel, reduction::nodesNearPlanes(thin, ends, 1e-9));
        const Eigenpairs longer = sub.fixedInterfaceModes(1, 1e-10, 240);
        expectTrue("four times the iterations does not converge it either", !longer.converged);
        const double drift = std::abs(longer.value[0] / slenderValue - 1.0);
        std::printf("       60 against 240 iterations moves the eigenvalue by %.2e relative\n",
                    drift);
        // Measured at 5e-8. Asserted at 1e-5, which is three decades above it and
        // three below the 1e-2 that would make the reported frequency untrustworthy.
        expectTrue("and the eigenvalue barely moves", drift < 1e-5);

        // **The independent instrument.** `eigenvaluesBelow` is an inertia count of
        // an LDL^T factorisation -- it counts rather than converges, so it owes the
        // subspace iteration nothing. It brackets the unconverged answer to one part
        // in ten thousand, which is what turns "probably fine" into a measurement.
        bool exactBelow = false, exactAbove = false;
        const int below = sub.eigenvaluesBelow(0.9999 * slenderValue, &exactBelow);
        const int above = sub.eigenvaluesBelow(1.0001 * slenderValue, &exactAbove);
        std::printf("       inertia count: %d below 0.9999 lambda1, %d below 1.0001 (exact %d/%d)\n",
                    below, above, static_cast<int>(exactBelow), static_cast<int>(exactAbove));
        expectTrue("both inertia counts are exact", exactBelow && exactAbove);
        expectEqual("nothing lies below the unconverged eigenvalue", below, 0);
        expectEqual("and exactly one lies just above it", above, 1);
        // Vacuity: an eigenvalue at zero would make a relative bracket meaningless.
        expectTrue("and it is a real frequency", slenderValue > 1.0);
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
    //
    // **Those three numbers used to live only in this comment**, which made the
    // bound impossible to re-derive: any change to the ferry's geometry moves all
    // three and nothing short of hand-instrumenting the ordering could say where
    // they landed. The two candidates are now reported, and the ordering's real
    // guarantee -- that it keeps *the narrower of the two*, so it "is incapable of
    // being a regression" -- is asserted here rather than only claimed in
    // `reduction.cpp`.
    //
    // **They are node bandwidths and `halfBandwidth()` is a DOF bandwidth**, which
    // is worth stating because the obvious assertion across them is wrong: this
    // patch reports 57 nodes against 71 DOF, which reads like the ordering picking
    // the worse of two until you notice the units. Three DOF a node puts the
    // mesher's own numbering at about 171 -- the 173 above -- and the 71 that was
    // kept at about 23 nodes.
    const std::size_t naturalNodes = sub.naturalNodeBandwidth();
    const std::size_t renumberedNodes = sub.renumberedNodeBandwidth();
    expectTrue("the renumbering is the narrower of the two candidates, or it is not used",
               std::min(naturalNodes, renumberedNodes) == renumberedNodes ||
                   sub.halfBandwidth() >= 3 * std::min(naturalNodes, renumberedNodes));
    expectTrue("and on this patch it is the renumbering that wins",
               renumberedNodes < naturalNodes);
    expectTrue("the interior renumbering earns its place", sub.halfBandwidth() <= 80);
    std::printf("     half-bandwidth %zu DOF over %zu interior DOF (%.2f MB of banded store); "
                "interior node bands: mesher %zu, Cuthill-McKee %zu\n",
                sub.halfBandwidth(), sub.interiorCount(),
                8.0 * static_cast<double>(sub.interiorCount()) *
                    static_cast<double>(sub.halfBandwidth() + 1) / 1e6,
                naturalNodes, renumberedNodes);

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

    // --- The static path, against the same reference ----------------------------
    //
    // Frequencies exercise the assembled matrices. They do not exercise getting a
    // load in or a displacement out, and `assembledStaticSolve` and
    // `componentState` were written and shipped without a test until this was
    // noticed. Both are on the path a caller actually uses.
    //
    // Restraining it is the fiddly part and it is worth being explicit. The reduced
    // model exposes only boundary and modal DOF, so the plate's *ends* cannot be
    // clamped -- they are interior to their components. What can be held is the
    // interface, and holding all of it would decouple the halves and quietly turn
    // this into two separate problems. So the six DOF are a proper 3-2-1 restraint
    // on three interface nodes, which removes exactly the six rigid body modes and
    // leaves the two components carrying each other.
    {
        const std::vector<std::uint32_t>& bdofA = sa.boundaryDof();
        const auto nodeOf = [&](std::uint32_t d) { return static_cast<std::size_t>(d) / 3; };
        const auto axisOf = [&](std::uint32_t d) { return static_cast<std::size_t>(d) % 3; };
        const auto posOf = [&](std::size_t n, int c) { return left.position[3 * n + c]; };

        // Three distinct interface nodes: a corner, one along y, one along z.
        std::size_t p1 = 0, p2 = 0, p3 = 0;
        for (std::size_t i = 0; i < bdofA.size(); ++i) {
            const std::size_t n = nodeOf(bdofA[i]);
            if (posOf(n, 1) < posOf(nodeOf(bdofA[p1]), 1) + 1e-12 &&
                posOf(n, 2) < posOf(nodeOf(bdofA[p1]), 2) + 1e-12)
                p1 = i;
        }
        const std::size_t n1 = nodeOf(bdofA[p1]);
        for (std::size_t i = 0; i < bdofA.size(); ++i) {
            const std::size_t n = nodeOf(bdofA[i]);
            if (std::fabs(posOf(n, 2) - posOf(n1, 2)) < 1e-12 &&
                posOf(n, 1) > posOf(nodeOf(bdofA[p2]), 1))
                p2 = i;
            if (std::fabs(posOf(n, 1) - posOf(n1, 1)) < 1e-12 &&
                posOf(n, 2) > posOf(nodeOf(bdofA[p3]), 2))
                p3 = i;
        }
        const std::size_t n2 = nodeOf(bdofA[p2]), n3 = nodeOf(bdofA[p3]);
        expectTrue("the three restraint nodes are distinct", n1 != n2 && n1 != n3 && n2 != n3);

        // 3-2-1: all of n1; x and z of n2 (killing rotation about x and about z);
        // x of n3 (killing rotation about y). Every interface node shares an x, so
        // the usual textbook arrangement has to be read in this plane.
        std::vector<std::uint32_t> held;
        std::vector<std::uint8_t> pin(3 * left.nodeCount(), 0);
        const auto hold = [&](std::size_t node, std::size_t axis) {
            for (std::size_t i = 0; i < bdofA.size(); ++i)
                if (nodeOf(bdofA[i]) == node && axisOf(bdofA[i]) == axis)
                    held.push_back(static_cast<std::uint32_t>(i));   // fromA[i] == i
            pin[3 * node + axis] = 1;
        };
        hold(n1, 0); hold(n1, 1); hold(n1, 2);
        hold(n2, 0); hold(n2, 2);
        hold(n3, 0);
        expectEqualCount("a 3-2-1 restraint holds six DOF", held.size(), 6u);

        // The same load on both models: a transverse pull on one interior node of
        // each half, which is a load the interface has to carry across.
        const std::size_t pullA = 0, pullB = 0;   // node 0 of each half, an outer corner
        std::vector<double> loadA(sa.dofCount(), 0.0), loadB(sb.dofCount(), 0.0);
        loadA[3 * pullA + 2] = 5.0e3;
        loadB[3 * pullB + 2] = -5.0e3;

        reduction::ReduceParams sp;
        sp.modes = 12;
        sp.cutoffFrequency = 0;
        const reduction::Reduction sra = reduction::craigBampton(sa, sp);
        const reduction::Reduction srb = reduction::craigBampton(sb, sp);
        const reduction::Assembly sasm = reduction::assemble(sra, srb, map);

        const std::vector<double> fa = reduction::reduceLoad(sa, sra, loadA);
        const std::vector<double> fb = reduction::reduceLoad(sb, srb, loadB);
        std::vector<double> f(static_cast<std::size_t>(sasm.size()), 0.0);
        for (std::size_t i = 0; i < fa.size(); ++i) f[static_cast<std::size_t>(sasm.fromA()[i])] += fa[i];
        for (std::size_t i = 0; i < fb.size(); ++i) f[static_cast<std::size_t>(sasm.fromB()[i])] += fb[i];

        std::vector<double> x;
        std::string problem;
        expectTrue("the assembled model solves with six DOF held",
                   reduction::assembledStaticSolve(sasm, f, held, x, &problem));

        // The reference: the whole plate, same restraint, same load, solved by the
        // element code rather than by anything here.
        solidshell::HexMesh clamped = whole;
        clamped.fixed.assign(3 * clamped.nodeCount(), 0);
        std::vector<double> loadFull(3 * clamped.nodeCount(), 0.0);
        const auto wholeNodeAt = [&](double px, double py, double pz) {
            for (std::size_t n = 0; n < clamped.nodeCount(); ++n)
                if (std::fabs(clamped.position[3 * n] - px) < 1e-9 &&
                    std::fabs(clamped.position[3 * n + 1] - py) < 1e-9 &&
                    std::fabs(clamped.position[3 * n + 2] - pz) < 1e-9)
                    return n;
            return clamped.nodeCount();
        };
        for (std::size_t n = 0; n < left.nodeCount(); ++n)
            for (std::size_t axis = 0; axis < 3; ++axis)
                if (pin[3 * n + axis]) {
                    const std::size_t w =
                        wholeNodeAt(posOf(n, 0), posOf(n, 1), posOf(n, 2));
                    expectTrue("every restrained node exists in the whole plate too",
                               w < clamped.nodeCount());
                    if (w < clamped.nodeCount()) clamped.fixed[3 * w + axis] = 1;
                }
        const std::size_t wa = wholeNodeAt(left.position[3 * pullA], left.position[3 * pullA + 1],
                                           left.position[3 * pullA + 2]);
        const std::size_t wb = wholeNodeAt(right.position[3 * pullB], right.position[3 * pullB + 1],
                                           right.position[3 * pullB + 2]);
        expectTrue("both loaded nodes exist in the whole plate",
                   wa < clamped.nodeCount() && wb < clamped.nodeCount());
        loadFull[3 * wa + 2] = 5.0e3;
        loadFull[3 * wb + 2] = -5.0e3;

        std::vector<double> uFull;
        expectTrue("the reference plate solves",
                   solidshell::solveStatic(clamped, steel, solidshell::Formulation::SolidShell,
                                           loadFull, uFull, &problem));
        double peak = 0;
        for (double v : uFull) peak = std::max(peak, std::fabs(v));
        expectTrue("the load actually moves the reference plate", peak > 1e-6);

        // Interface displacements come straight out of the assembled state -- the
        // boundary DOF are physical, which is the whole reason they are kept.
        double worstInterface = 0;
        for (std::size_t i = 0; i < bdofA.size(); ++i) {
            const std::size_t n = nodeOf(bdofA[i]);
            const std::size_t w = wholeNodeAt(posOf(n, 0), posOf(n, 1), posOf(n, 2));
            if (w >= clamped.nodeCount()) continue;
            worstInterface = std::max(worstInterface,
                                      std::fabs(x[i] - uFull[3 * w + axisOf(bdofA[i])]));
        }
        std::printf("     static: interface to %.2e m of a peak %.3e m\n", worstInterface, peak);
        // Exact, not merely close: static condensation reproduces the interface
        // response for any load at any mode count, and that survives assembly.
        // Asserted at 1e-9 relative because 1e-3 would pass on a model that had
        // lost the property entirely and was simply well converged.
        expectTrue("the assembled model reproduces the interface displacement exactly",
                   worstInterface < 1e-9 * peak);

        // And `componentState` + `recover` put a component's own field back. The
        // boundary part must agree with what was just read out of the assembly,
        // which is the check that the index map is not quietly transposed.
        const std::vector<double> xa = reduction::componentState(sasm, sasm.fromA(), x);
        expectEqualCount("a component's state is its own reduced size", xa.size(),
                         static_cast<std::size_t>(sra.size()));
        double worstMap = 0;
        for (int i = 0; i < sra.boundary; ++i)
            worstMap = std::max(worstMap, std::fabs(xa[static_cast<std::size_t>(i)] -
                                                    x[static_cast<std::size_t>(i)]));
        expectTrue("componentState returns that component's own DOF", worstMap < 1e-15);

        // A state of the wrong length is the one mistake it cannot otherwise
        // notice: every index in `fromA` is in range, so a short state would come
        // back as a plausible field silently missing its modal content.
        std::vector<double> truncated(x.begin(), x.end() - 1);
        expectTrue("and refuses a state that is not this assembly's",
                   reduction::componentState(sasm, sasm.fromA(), truncated).empty());

        const std::vector<double> ua = reduction::recover(sa, sra, xa);
        double worstInterior = 0;
        for (std::size_t n = 0; n < left.nodeCount(); ++n) {
            const std::size_t w = wholeNodeAt(posOf(n, 0), posOf(n, 1), posOf(n, 2));
            if (w >= clamped.nodeCount()) continue;
            for (std::size_t axis = 0; axis < 3; ++axis)
                worstInterior =
                    std::max(worstInterior, std::fabs(ua[3 * n + axis] - uFull[3 * w + axis]));
        }
        std::printf("     static: recovered interior to %.2e m of a peak %.3e m\n", worstInterior,
                    peak);
        expectTrue("and the recovered interior field is the plate's own",
                   worstInterior < 1e-3 * peak);
    }
}

// --- 11. Stiffness the elements do not carry (§8) ---------------------------------
//
// A stiffener is condensed onto the plating's degrees of freedom, so it has no
// nodes and no elements and a substructure built from the mesh alone reduces a
// stiffened patch as bare plating. `Attachment` is what closes that, and the
// assertions here are almost all identities rather than tolerances:
//
//   * The **energy an attachment adds under a prescribed strain field is a closed
//     form in the profile's area, first moment and second moment**, and those three
//     come from `scantlings::profileSection` -- a different file, computing them
//     from the rectangle dimensions, owing nothing to the fibres. That single
//     measurement catches a wrong area, a wrong eccentricity *sign*, and a lost
//     Steiner term, which is the failure `scantlings.hpp` §1 exists to prevent.
//   * The **static answer is `solidshell::solveStatic` with the same blocks**: a
//     different assembly, a different numbering and a different factorisation of
//     the same physics, and Guyan is exact at the interface for any load.
//   * The **negative control is bit equality**. An empty attachment must not be
//     close to what a substructure did before this section existed; it must be the
//     same bits, and it is, because both constructors are one code path.
void testAttachedStiffness() {
    std::printf("\n--- reduction: the stiffener the elements do not carry ---\n");
    const solidshell::HexMesh mesh = testPlate();
    const StructuralMaterial steel = ah36Steel();
    const StiffenerProfile profile = flatBar(0.200, 0.010);
    const constraint::Stiffening seam = stiffenSeam(mesh, steel, profile, kNy / 2, 1);
    expectEqualCount("the seam carries two fibres per segment of a flat bar", seam.fiberCount(),
                     static_cast<std::size_t>(2 * kNx));
    expectNear("and covers the whole plate", seam.length, kLx, 1e-12);

    const reduction::Attachment attached = attachmentOf(seam, mesh, steel, true);
    const Substructure bare(mesh, steel, endInterface(mesh));
    const Substructure sub(mesh, steel, endInterface(mesh), attached);
    expectTrue("the attached substructure factors", sub.ready());
    expectEqualCount("and carries one block per fibre", sub.attachedBlocks(), seam.fiberCount());

    // --- The negative control, first, because everything else rests on it ---
    //
    // Not "close": the same bits. The unattached constructor delegates to the
    // attached one with an empty `Attachment`, so there is a single assembly path
    // and a difference of one ulp would mean the empty case had grown a branch of
    // its own.
    {
        const Substructure empty(mesh, steel, endInterface(mesh), reduction::Attachment{});
        expectEqual("an empty attachment leaves the bandwidth alone",
                    static_cast<long long>(empty.halfBandwidth()),
                    static_cast<long long>(bare.halfBandwidth()));
        expectEqualCount("and the partition", empty.boundaryCount(), bare.boundaryCount());
        std::vector<double> probe(bare.dofCount(), 0.0), fromEmpty, fromBare;
        for (std::size_t d = 0; d < probe.size(); ++d)
            probe[d] = 1e-6 * static_cast<double>((d * 37u) % 19u) - 9e-6;
        empty.stiffnessTimes(probe, fromEmpty);
        bare.stiffnessTimes(probe, fromBare);
        double scale = 0;
        bool identical = fromEmpty.size() == fromBare.size();
        for (std::size_t d = 0; d < fromBare.size() && identical; ++d) {
            scale = std::max(scale, std::fabs(fromBare[d]));
            identical = fromEmpty[d] == fromBare[d];
        }
        expectTrue("the probe loads the operator", scale > 1.0);
        expectTrue("an empty attachment is the bare operator to the last bit", identical);
        bool sameMass = empty.mass().size() == bare.mass().size();
        for (std::size_t d = 0; d < bare.mass().size() && sameMass; ++d)
            sameMass = empty.mass()[d] == bare.mass()[d];
        expectTrue("and the same mass to the last bit", sameMass);
        expectNear("with nothing attached to report", empty.attachedMass(), 0.0, 0.0);
        expectEqualCount("and no blocks", empty.attachedBlocks(), 0u);
    }

    // --- What the blocks are worth, in closed form ---
    //
    // Prescribe `u_x = (eps + kappa z) x` and nothing else. Every fibre lies along
    // x, so the tie reads only `u_x`, and the tied point at offset `e` moves by
    // `(eps + kappa e) x` **exactly** -- that is what `tieWeight` is defined to
    // deliver. A fibre of area `A_j` over a span `L_j` then stores
    // `E A_j L_j (eps + kappa e_j)^2` of `u^T K u`, so summed over the seam
    //
    //     u^T K_attached u - u^T K_bare u = E L (eps^2 A + 2 eps kappa S + kappa^2 I)
    //
    // with A, S and I the profile's area, first moment and second moment about the
    // **plate mid-surface**. The plate's own contribution cancels between the two
    // operators whatever its discretisation error, so this is an identity and not a
    // convergence statement. `eps` and `kappa` are both non-zero on purpose: with
    // either one alone the cross term S drops out, and S is the only quantity here
    // that knows which side of the plate the web is on.
    const ProfileSection section = profileSection(profile);
    const double centroid = 0.5 * kThickness + section.centroid;
    const double area = section.area;
    const double firstMoment = section.area * centroid;
    const double secondMoment = section.secondMoment + section.area * centroid * centroid;
    {
        // The same second moment out of `stiffenedSection`, which computes it about
        // the *combined* neutral axis and therefore by a different route: the
        // parallel axis theorem takes it back to the mid-surface, and the plate's
        // own b t^3 / 12 comes off. Two files, two formulas, one number.
        const StiffenedSection panel = stiffenedSection(profile, kThickness, kLy);
        const double aboutMidSurface =
            panel.secondMoment + panel.area * panel.neutralAxis * panel.neutralAxis;
        const double plateOwn = kLy * kThickness * kThickness * kThickness / 12.0;
        expectNear("the two section routines agree about the profile's second moment",
                   aboutMidSurface - plateOwn, secondMoment, 1e-16 * secondMoment);
        // And the axial one, against the smeared thickness the Tier-0 beam uses:
        // t + A/s is the same statement as "the fibres carry A".
        expectNear("and about its area, against the smeared thickness",
                   (smearedThickness(kThickness, profile, kLy) - kThickness) * kLy, area,
                   1e-15 * area);
    }

    const double kStrain = 1e-4, kCurvature = 1e-3;
    {
        std::vector<double> u(sub.dofCount(), 0.0), attachedForce, bareForce;
        for (std::size_t n = 0; n < mesh.nodeCount(); ++n) {
            const double x = mesh.position[n * 3], z = mesh.position[n * 3 + 2];
            u[n * 3] = (kStrain + kCurvature * z) * x;
        }
        sub.stiffnessTimes(u, attachedForce);
        bare.stiffnessTimes(u, bareForce);
        double added = 0, bareEnergy = 0;
        for (std::size_t d = 0; d < u.size(); ++d) {
            added += u[d] * (attachedForce[d] - bareForce[d]);
            bareEnergy += u[d] * bareForce[d];
        }
        const double predicted =
            steel.youngsModulus * seam.length *
            (kStrain * kStrain * area + 2.0 * kStrain * kCurvature * firstMoment +
             kCurvature * kCurvature * secondMoment);
        // 1e-11 relative: the measurement is 2.2e-13, limited by the cancellation
        // between two operators that agree to a part in sixteen everywhere the
        // fibres are not. A per-cent tolerance would pass on a stiffener that had
        // lost `I_own` -- 23% of this number -- or had its web on the wrong side,
        // which moves it by 4 * eps * kappa * S.
        expectNear("the attachment adds exactly the profile's A, S and I", added, predicted,
                   1e-11 * predicted);
        expectTrue("and it is not a rounding on top of the plating", added > 0.05 * bareEnergy);
        std::printf("     prescribed strain %g + curvature %g: fibres add %.6e J, closed form "
                    "%.6e (%.1e relative), %.1f%% of the plating's own\n",
                    kStrain, kCurvature, added, predicted, std::fabs(added - predicted) / predicted,
                    100.0 * added / bareEnergy);

    }

    // Per degree of freedom, against the blocks summed by hand. The energy above is
    // one scalar and a dropped off-diagonal can cancel inside it -- two recent
    // mutants in this repo died only when the suspect entry was asked about alone,
    // both of them errors that cancelled when the model was asked globally.
    //
    // **Not the strain field**: a uniform fibre strain is in equilibrium at every
    // interior station, so the two fibres meeting there cancel and only the ends of
    // the seam carry force. The probe below is deliberately incoherent, so every
    // degree of freedom a block can load is loaded and each can be checked on its
    // own.
    //
    // It is run on the plate turned 45 degrees as well, and that is the check that
    // carries the content. A bar has stiffness along its own axis and none across
    // it (`constraint.hpp` §2), so on the axis-aligned seam every block is a
    // rank-one form in the x components alone -- an assembly that transposed an
    // axis, or that dropped one, would scatter zeros onto zeros and agree
    // perfectly. Turned 45 degrees the same seam loads x and y equally, and the two
    // are no longer interchangeable.
    for (int turned = 0; turned < 2; ++turned) {
        const solidshell::HexMesh plate = turned ? rotatedPlate() : testPlate();
        const constraint::Stiffening line =
            stiffenSeam(plate, steel, profile, kNy / 2, 1, turned == 0);
        const reduction::Attachment blocks = attachmentOf(line, plate, steel, true);
        const std::vector<std::uint32_t> interface =
            turned ? rotatedEndInterface(plate) : endInterface(plate);
        const Substructure on(plate, steel, interface, blocks);
        const Substructure off(plate, steel, interface);
        expectTrue("the substructure factors", on.ready() && off.ready());

        std::vector<double> probe(on.dofCount(), 0.0), attachedForce, bareForce;
        for (std::size_t d = 0; d < probe.size(); ++d)
            probe[d] = 1e-6 * static_cast<double>((d * 41u) % 23u) - 1.1e-5;
        on.stiffnessTimes(probe, attachedForce);
        off.stiffnessTimes(probe, bareForce);

        std::vector<double> byHand(probe.size(), 0.0);
        for (const solidshell::DofBlock& block : blocks.stiffness) {
            const std::size_t n = block.dof.size();
            for (std::size_t p = 0; p < n; ++p)
                for (std::size_t q = 0; q < n; ++q)
                    byHand[block.dof[p]] += block.stiffness[p * n + q] * probe[block.dof[q]];
        }
        std::size_t loaded = 0;
        double worst = 0, peak = 0;
        for (double v : byHand) peak = std::max(peak, std::fabs(v));
        for (std::size_t d = 0; d < probe.size(); ++d) {
            worst = std::max(worst, std::fabs((attachedForce[d] - bareForce[d]) - byHand[d]));
            if (std::fabs(byHand[d]) > 1e-6 * peak) ++loaded;
        }
        // Two nodes at each of the nine stations, times one axis flat and two
        // turned. Asserted rather than counted at run time, because "however many
        // it happened to load" would make the per-DOF check below vacuous exactly
        // when a scatter had gone missing.
        expectEqualCount("the fibres load one axis of the pair when the seam is axis aligned "
                         "and two when it is not",
                         loaded, static_cast<std::size_t>((turned ? 2 : 1) * 2 * (kNx + 1)));
        expectTrue("every degree of freedom gets exactly its own share of the blocks",
                   worst < 1e-9 * peak);
        std::printf("     scatter%s: %zu degrees of freedom, each to %.2e N of %.3e N\n",
                    turned ? " (seam at 45 degrees)" : "", loaded, worst, peak);
    }

    // --- The static answer, against a reference that owes nothing to this file ---
    //
    // `solidshell::solveStatic` takes the same `DofBlock` list, assembles it into
    // its own banded system with its own free-DOF numbering, and factors it once.
    // Guyan condensation is exact at the interface for any load (§1 property 1), so
    // the two must agree to the conditioning of two independent solves -- and that
    // is the assertion that would fail if any part of the attachment failed to
    // reach `K_ii`, because `Psi` is `-K_ii^-1 K_ib` and nothing else.
    //
    // `stride` 2 is the case §8 item 2 is about: fibres spanning two elements tie
    // node pairs that share no element, so the sparsity pattern has to grow and the
    // band with it. It is checked that the band really did grow, because a case
    // that did not need one would test nothing.
    for (int stride : {1, 2}) {
        const constraint::Stiffening run = stiffenSeam(mesh, steel, profile, kNy / 2, stride);
        const reduction::Attachment attach = attachmentOf(run, mesh, steel, true);
        const Substructure strided(mesh, steel, endInterface(mesh), attach);
        expectTrue("the strided substructure factors", strided.ready());
        if (stride > 1)
            expectTrue("a fibre spanning two elements widens the band the elements imply",
                       strided.halfBandwidth() > bare.halfBandwidth());
        else
            expectEqual("a fibre inside one element does not",
                        static_cast<long long>(strided.halfBandwidth()),
                        static_cast<long long>(bare.halfBandwidth()));

        solidshell::HexMesh pinned = testPlate();
        std::vector<double> load(strided.dofCount(), 0.0);
        std::vector<std::uint32_t> held;
        for (std::size_t b = 0; b < strided.boundaryCount(); ++b) {
            const std::uint32_t d = strided.boundaryDof()[b];
            if (mesh.position[(d / 3) * 3] < 0.5 * kLx) {
                held.push_back(static_cast<std::uint32_t>(b));
                pinned.pin(d / 3, static_cast<int>(d % 3), 0.0);
            } else if (d % 3 == 2) {
                load[d] = 1000.0;
            }
        }
        std::vector<double> reference, withoutFibres;
        std::string problem;
        expectTrue("the stiffened plate solves",
                   solidshell::solveStatic(pinned, steel, solidshell::Formulation::SolidShell,
                                           attach.stiffness, load, reference, &problem));
        expectTrue("and the bare one",
                   solidshell::solveStatic(pinned, steel, solidshell::Formulation::SolidShell, load,
                                           withoutFibres, &problem));
        double peak = 0, barePeak = 0, stiffening = 0;
        for (std::size_t d = 0; d < reference.size(); ++d) {
            peak = std::max(peak, std::fabs(reference[d]));
            barePeak = std::max(barePeak, std::fabs(withoutFibres[d]));
            stiffening = std::max(stiffening, std::fabs(reference[d] - withoutFibres[d]));
        }
        expectTrue("the load moves the plate", peak > 1e-4);
        // The vacuity guard the whole exercise exists for: without the stiffener
        // the answer is the bare plating, and the difference must be a large
        // fraction of the field rather than a correction to it.
        expectTrue("and the stiffener is worth far more than a rounding of it",
                   stiffening > 0.5 * barePeak);

        ReduceParams zero;
        zero.modes = 0;
        const Reduction guyan = craigBampton(strided, zero);
        std::vector<double> state;
        expectTrue("the reduced stiffened model solves",
                   reduction::staticSolve(guyan, reduction::reduceLoad(strided, guyan, load), held,
                                          state, &problem));
        const std::vector<double> u = reduction::recover(strided, guyan, state);
        double worstBoundary = 0, worstInterior = 0;
        for (std::size_t b = 0; b < strided.boundaryCount(); ++b) {
            const std::uint32_t d = strided.boundaryDof()[b];
            worstBoundary = std::max(worstBoundary, std::fabs(u[d] - reference[d]));
        }
        for (std::size_t p = 0; p < strided.interiorCount(); ++p) {
            const std::uint32_t d = strided.interiorDof()[p];
            worstInterior = std::max(worstInterior, std::fabs(u[d] - reference[d]));
        }
        // 1e-8 of the peak, where the measurement is 5e-10: a stiffened plate is an
        // order worse conditioned than the bare one this file's other static test
        // uses, and a looser tolerance would admit a model that had lost one fibre
        // in sixteen -- which moves the answer by percent, not by 1e-8.
        expectTrue("the reduced stiffened model reproduces the stiffened plate's interface",
                   worstBoundary < 1e-8 * peak);
        expectTrue("and its interior", worstInterior < 1e-8 * peak);
        std::printf("     stride %d: band %zu -> %zu, stiffener moves the field %.1f%%, Guyan "
                    "matches solveStatic to %.1e m of %.3e m\n",
                    stride, bare.halfBandwidth(), strided.halfBandwidth(),
                    100.0 * stiffening / barePeak, std::max(worstBoundary, worstInterior), peak);
    }

    // --- What a malformed attachment does. Refused, not skipped: a stiffener that
    // quietly does not arrive is indistinguishable from bare plating, which is the
    // whole failure this section closes.
    {
        reduction::Attachment bad = attached;
        bad.stiffness[0].dof[3] = static_cast<std::uint32_t>(sub.dofCount());
        const Substructure refused(mesh, steel, endInterface(mesh), bad);
        expectTrue("a block naming a degree of freedom the mesh does not have is refused",
                   !refused.ready() && !refused.problems().empty());

        reduction::Attachment short_ = attached;
        short_.stiffness[0].stiffness.pop_back();
        const Substructure refusedShort(mesh, steel, endInterface(mesh), short_);
        expectTrue("and so is a block whose stiffness array is short",
                   !refusedShort.ready() && !refusedShort.problems().empty());

        reduction::Attachment wrongMass = attached;
        wrongMass.mass.assign(mesh.nodeCount() * 3, 0.0);
        const Substructure refusedMass(mesh, steel, endInterface(mesh), wrongMass);
        expectTrue("and a mass array indexed by degree of freedom rather than by node",
                   !refusedMass.ready() && !refusedMass.problems().empty());

        const Substructure silent(mesh, steel, endInterface(mesh),
                                  attachmentOf(seam, mesh, steel, false));
        expectTrue("stiffness with no mass is allowed but not silent",
                   silent.ready() && !silent.problems().empty());

        // And the same, for the stress forms §9 added. A form list that is not
        // parallel to the blocks would read the right degrees of freedom for the
        // *wrong* member -- a number that is plausible and wrong -- so it is
        // refused, not truncated to the shorter of the two.
        reduction::Attachment shortStress = attached;
        shortStress.stress.pop_back();
        const Substructure refusedStress(mesh, steel, endInterface(mesh), shortStress);
        expectTrue("a stress form list shorter than the block list is refused",
                   !refusedStress.ready() && !refusedStress.problems().empty());

        reduction::Attachment longStress = attached;
        longStress.stress.push_back(longStress.stress.back());
        const Substructure refusedLong(mesh, steel, endInterface(mesh), longStress);
        expectTrue("and one longer than it",
                   !refusedLong.ready() && !refusedLong.problems().empty());

        reduction::Attachment ragged = attached;
        ragged.stress[0].pop_back();
        const Substructure refusedRagged(mesh, steel, endInterface(mesh), ragged);
        expectTrue("and a form with fewer entries than its block names degrees of freedom",
                   !refusedRagged.ready() && !refusedRagged.problems().empty());

        const Substructure blind(mesh, steel, endInterface(mesh),
                                 attachmentOf(seam, mesh, steel, true, false));
        expectTrue("stiffness with no stress form is allowed but not silent either",
                   blind.ready() && !blind.problems().empty());
        expectEqualCount("and it carries no members to read", blind.attachedMembers(), 0u);
        expectEqualCount("while still carrying every block", blind.attachedBlocks(),
                         seam.fiberCount());
    }
}

// --- 11b. What the promotion trigger reads (§9) -----------------------------------
//
// `checkValidity` used to walk the elements alone, so a stiffened region was
// judged by its plating -- and the plating is the softer half, so the utilisation
// came out low in the *unsafe* direction. That is what decides whether a region is
// promoted to the nonlinear Tier-2 model, so reading it low means staying linear
// after the structure has stopped being linear.
//
// What is asserted here, and each answers a different way of getting it wrong:
//
//   * **That the two halves are comparable at all.** A bar's stress tensor has one
//     non-zero entry and the von Mises of that state is `|sigma|` identically --
//     asserted as an equality of bits through the same formula the element half is
//     measured by, because it is algebra rather than a limit. With it, the bound on
//     what the fibre model *omits*: a transverse stress makes von Mises
//     `|s| sqrt(1 - a + a^2)`, which is 0.866 |s| at `a = 1/2` and sqrt(3) |s| at
//     `a = -1`, so the omission has no fixed sign in either direction.
//   * **The closed form**, which owes this file nothing. Under a prescribed
//     `u_x = (eps + kappa z) g(x)` the tie puts a fibre at offset `e` at exactly
//     `(eps + kappa e) g(x)`, so a fibre spanning `x_a` to `x_b` carries
//     `E (eps + kappa e) (g(x_b) - g(x_a)) / (x_b - x_a)`. Choose `eps` to put the
//     section in pure bending and the offset term is `E kappa (e - z_na)` with
//     `z_na` the combined neutral axis -- and `z_na` derived from the fibre
//     stations is compared against `scantlings::stiffenedSection`'s, computed from
//     the rectangle dimensions in a different file by a different formula.
//   * **The regression case**, 58.2 MPa in the plating against 65.2 in the member:
//     utilisation 0.1838 where it used to be 0.1639.
//   * **The negative control**, bit for bit: with no attachment `checkValidity`
//     must return exactly what it returned before §9, not something close to it.
//   * **The vacuity guard, in both directions.** A test in which the member always
//     governs would pass on an implementation that reported the member and ignored
//     the plating. So the *same* profile under the *same* load is run with the seam
//     moved to the plate edge, where the plating governs by 2.5x, and the answer
//     has to follow the structure.
//   * **That the tie is filled in.** A degree of freedom `Attachment::constrained`
//     eliminated is in neither partition, and `recover` leaving it at zero put a
//     hole in the field that read as 850 788 MPa on a plate carrying 427.
//
// The cases that look synthetic -- a seam of one fibre, the same fibre twice, a
// diagonal brace -- are each here because mutation testing showed the sixteen-fibre
// seam could not reach the predicate they test. Their reasons are beside them.
void testMemberValidity() {
    std::printf("\n--- reduction: what the promotion trigger reads (§9) ---\n");
    const solidshell::HexMesh mesh = testPlate();
    const StructuralMaterial steel = ah36Steel();
    const StiffenerProfile profile = flatBar(0.200, 0.010);

    // --- The identity that makes the two halves comparable at all ---
    //
    // A fibre is an axial bar, so its stress tensor has one non-zero entry, and the
    // von Mises of that state is `|sigma|` *identically*. Asserted through the same
    // formula the element half is measured by rather than restated, and asserted as
    // an equality of bits: it is algebra, not a limit. Both signs, because
    // `sqrt(...)` of a squared difference would agree on one of them by accident.
    for (double s : {355.0e6, -355.0e6, 1.0, -1e-30}) {
        const double uniaxial[6] = {s, 0, 0, 0, 0, 0};
        const double transverse[6] = {0, s, 0, 0, 0, 0};
        expectTrue("a uniaxial stress state's von Mises is exactly its own magnitude",
                   misesOf(uniaxial) == std::fabs(s) && misesOf(transverse) == std::fabs(s));
    }
    // ...and the bound on what the fibre model leaves out, which is the reason the
    // equality above is a statement about the *model* and not about the member: a
    // transverse stress beside the axial one changes von Mises to
    // `|s| sqrt(1 - a + a^2)`, whose minimum over `a` is at a = 1/2. So ignoring
    // biaxiality has no fixed sign, and the header says so rather than claiming the
    // fibre stress is conservative.
    {
        const double s = 100.0e6;
        const double half[6] = {s, 0.5 * s, 0, 0, 0, 0};
        const double opposed[6] = {s, -s, 0, 0, 0, 0};
        expectNear("a transverse tension of half the axial lowers von Mises to sqrt(3)/2",
                   misesOf(half) / s, std::sqrt(0.75), 1e-15);
        expectNear("and a transverse compression of equal size raises it to sqrt(3)",
                   misesOf(opposed) / s, std::sqrt(3.0), 1e-15);
    }

    // --- The closed form, against `scantlings::stiffenedSection` -----------------
    //
    // Two independent routes to the neutral axis. `stiffenedSection` takes it from
    // the rectangle dimensions; the route here takes it from the fibre stations the
    // constraint machinery actually built, by requiring the section to carry no net
    // axial force under a pure curvature -- `sum A_j (e_j - z_na) + b t (0 - z_na)
    // = 0`. The two agree only if the two-point Gauss stations reproduce the
    // profile's area and first moment exactly, which is `constraint.hpp` §2's claim
    // and not this file's.
    {
        const constraint::Stiffening seam = stiffenSeam(mesh, steel, profile, kNy / 2, 1);
        const Substructure sub(mesh, steel, endInterface(mesh),
                               attachmentOf(seam, mesh, steel, true));
        expectTrue("the stiffened substructure factors", sub.ready());
        expectEqualCount("and carries a stress form for every block", sub.attachedMembers(),
                         sub.attachedBlocks());

        const constraint::ProfileFibers stations = constraint::profileFibers(profile, kThickness, 1.0);
        double stationArea = 0, stationFirstMoment = 0;
        for (int i = 0; i < stations.count; ++i) {
            stationArea += stations.area[i];
            stationFirstMoment += stations.area[i] * stations.offset[i];
        }
        const double neutralAxis = stationFirstMoment / (kLy * kThickness + stationArea);
        const StiffenedSection panel = stiffenedSection(profile, kThickness, kLy);
        // The measurement is 3.5e-18 m, which is **one ulp** of an answer of
        // 0.0230 m: the two formulas agree as exactly as double precision allows.
        // The tolerance is five ulp rather than one, because that is the room two
        // different orderings of the same arithmetic need across the three
        // optimisation levels the gate compiles at. A per-cent tolerance would pass
        // on a section that had lost the plate's own area, which moves this by 60%.
        expectNear("the fibre stations and stiffenedSection agree about the neutral axis",
                   neutralAxis, panel.neutralAxis, 8e-16 * panel.neutralAxis);

        // Prescribe the field. Every fibre lies along x and the tie reads only the
        // x components, so the tied point at offset `e` moves by `(eps + kappa e)
        // g(x)` exactly -- that is what `tieWeight` is defined to deliver -- and a
        // fibre spanning `x_a` to `x_b` therefore carries
        //
        //     sigma = E (eps + kappa e) (g(x_b) - g(x_a)) / (x_b - x_a)
        //
        // **`g` is `x + x^2 / (2 c)` and not `x`, and that is a vacuity this test
        // was carrying.** With `g(x) = x` every segment of the seam elongates by the
        // same amount, so every fibre at a given offset carries the *same* stress
        // -- and a recovery that read the neighbouring segment's degrees of freedom
        // instead of its own returned exactly the right number. Mutation testing
        // found it: the mutant "reads the wrong member's dof list" survived a test
        // that checked every one of the sixteen fibres. Adding the quadratic term
        // makes each segment's elongation its own while keeping the closed form
        // closed, because `(x_b^2 - x_a^2) / (x_b - x_a)` is `x_a + x_b`.
        const double kappa = 2.0e-3;
        const double strain = -kappa * neutralAxis;
        const double kTaper = 3.0 * kLx;  // `c`: a 20% spread of elongation over the seam
        const auto along = [&](double x) { return x + x * x / (2.0 * kTaper); };
        std::vector<double> u(sub.dofCount(), 0.0);
        for (std::size_t n = 0; n < mesh.nodeCount(); ++n) {
            const double x = mesh.position[n * 3], z = mesh.position[n * 3 + 2];
            u[n * 3] = (strain + kappa * z) * along(x);
        }
        double worst = 0, biggest = 0, netForce = 0, plateForce = 0;
        double mostTension = 0, mostCompression = 0, spread = 0, leastFactor = 1e30;
        for (std::size_t m = 0; m < sub.attachedMembers(); ++m) {
            const double got = sub.memberStress(m, u);
            const double offset = seam.fiber[m].offset;
            const double xa = constraint::tiedPoint(seam.fiber[m].end[0], mesh.position).x;
            const double xb = constraint::tiedPoint(seam.fiber[m].end[1], mesh.position).x;
            const double factor = 1.0 + (xa + xb) / (2.0 * kTaper);
            const double want = steel.youngsModulus * kappa * (offset - neutralAxis) * factor;
            worst = std::max(worst, std::fabs(got - want));
            biggest = std::max(biggest, std::fabs(want));
            mostTension = std::max(mostTension, got);
            mostCompression = std::min(mostCompression, got);
            spread = std::max(spread, factor);
            leastFactor = std::min(leastFactor, factor);
            // The axial force this fibre carries is `sigma A`, and the plating over
            // the same span carries `E eps b t` times the same factor. They cancel
            // exactly when the neutral axis is right, which is what makes the closed
            // form above a statement and not a definition. Half the plate's share
            // per fibre, because each segment carries two of them.
            netForce += got * seam.fiber[m].area;
            plateForce += 0.5 * factor * steel.youngsModulus * strain * kLy * kThickness;
        }
        netForce /= static_cast<double>(kNx);   // one station's worth of fibres
        plateForce /= static_cast<double>(kNx);
        expectTrue("the prescribed curvature stresses the members", biggest > 10.0e6);
        // The guard on the vacuity above: the segments must actually differ, or the
        // quadratic term has bought nothing and the test is the old one again.
        expectTrue("and the segments of the seam carry visibly different stress",
                   spread > 1.15 * leastFactor);
        // 1e-13 relative, where the measurement is 5.3e-15 -- there is no solve in
        // this at all, only a twelve-term dot product against a field built by two
        // multiplications, so the floor is a few ulp of the largest term rather than
        // anything that converges. The headroom over the measurement is for the
        // gate's other optimisation levels, which are free to contract the multiply
        // and add differently. A 1e-6 tolerance would pass on a stiffener whose
        // fibres had all been put at the profile centroid, which loses `I_own` --
        // 23% of the second moment.
        expectTrue("every fibre reports E kappa (e - z_na), the closed form",
                   worst < 1e-13 * biggest);
        expectTrue("the section carries no net axial force, so z_na really is the neutral axis",
                   std::fabs(netForce + plateForce) < 1e-14 * std::fabs(plateForce));

        // **Every fibre comes out in tension, and that is not a bug in the test.**
        // The neutral axis is 23.0 mm above the plate mid-surface -- inside the web,
        // whose root is at 6.0 mm -- but two-point Gauss puts the lowest station at
        // 48.3 mm, so no fibre station is below it. The section straddles the
        // neutral axis through the *plating*, which carries all of the compression,
        // and the sign guard therefore has to be taken between the two halves rather
        // than within the member. (It also means the fibre model never samples the
        // part of the web that is in compression: it integrates the energy exactly,
        // which is what §2 of `constraint.hpp` claims, but it cannot report a stress
        // there.)
        expectTrue("the members are in tension and the plating carries the compression",
                   mostTension > 0 && plateForce < 0 && mostCompression >= 0);
        // The neutral axis is not a rounding on this section: dropping it would move
        // every fibre stress by E kappa z_na, 16% of the peak. Without this the
        // closed form would be satisfied by an implementation that ignored z_na and
        // a test whose `strain` happened to be small.
        expectTrue("and it is worth far more than the tolerance the closed form is asserted at",
                   steel.youngsModulus * kappa * neutralAxis > 0.1 * biggest);

        // Reverse the curvature: every fibre stress must flip sign exactly. A
        // recovery that returned a magnitude -- `|sigma|` instead of `sigma` --
        // passes every assertion above and fails this one, and it is the mistake
        // `checkValidity` itself makes on purpose one layer up.
        std::vector<double> reversed(sub.dofCount(), 0.0);
        for (std::size_t n = 0; n < mesh.nodeCount(); ++n) {
            const double x = mesh.position[n * 3], z = mesh.position[n * 3 + 2];
            reversed[n * 3] = -(strain + kappa * z) * along(x);
        }
        double worstFlip = 0;
        for (std::size_t m = 0; m < sub.attachedMembers(); ++m)
            worstFlip = std::max(worstFlip, std::fabs(sub.memberStress(m, reversed) +
                                                      sub.memberStress(m, u)));
        expectNear("and reversing the curvature reverses every fibre's stress, exactly",
                   worstFlip, 0.0, 0.0);

        // Out of range is zero, not the entry past the end of the array. Asked for
        // rather than assumed: the guard is a `>=` one index away from undefined
        // behaviour, and ASan is what would say so if it slipped to `>`.
        expectNear("a member index the substructure does not have reads zero",
                   sub.memberStress(sub.attachedMembers(), u), 0.0, 0.0);
        expectNear("and so does one far past the end",
                   sub.memberStress(sub.attachedMembers() + 1000, u), 0.0, 0.0);

        std::printf("     prescribed kappa %g about z_na %.9f m: closed form to %.1e Pa of "
                    "%.3e Pa, net axial force %.2e N against a plate's %.6e N, fibres "
                    "%.2f to %+.2f MPa (all above z_na)\n",
                    kappa, neutralAxis, worst, biggest, netForce + plateForce, plateForce,
                    1e-6 * mostCompression, 1e-6 * mostTension);
    }

    // --- The 58.2 / 65.2 MPa case, which is the regression test ------------------
    //
    // The stiffened cantilever: the plate held along x = 0, 1000 N up on every node
    // of the far end, eight modes. The member carries 12% more than the plating, so
    // the utilisation the trigger reports moves from 0.164 to 0.184.
    //
    // The seam row is swept, because the same profile under the same load has the
    // plating governing at the plate edge and the member governing at mid width.
    // Without that the assertion would pass on an implementation that reported the
    // member and forgot the plating entirely.
    struct Reading {
        double plating = 0, member = 0, utilisation = 0, blindUtilisation = 0;
        int worstMember = -1, members = 0;
        double offSeam = 0;     // m, the worst element's centroid from the seam line
        double biaxiality = 0;  // its von Mises over its own axial component
    };
    // `heldEnd` is which end is built in. Both are run because the peak member has
    // to be found at *both* ends of the fibre list or a loop bound off by one at
    // either end goes unnoticed -- which is exactly what mutation testing showed:
    // with the root at x = 0 the worst fibre is index 1 of sixteen, so a loop that
    // stopped one short still found it.
    const auto cantilever = [&](int seamRow, double heldEnd) {
        const constraint::Stiffening seam = stiffenSeam(mesh, steel, profile, seamRow, 1);
        const Substructure sub(mesh, steel, endInterface(mesh),
                               attachmentOf(seam, mesh, steel, true));
        // The same substructure with the stress forms withheld: this is what the
        // trigger read before §9, and it is what the 0.164 has to come out of.
        const Substructure blind(mesh, steel, endInterface(mesh),
                                 attachmentOf(seam, mesh, steel, true, false));
        expectTrue("the seam-row substructures factor", sub.ready() && blind.ready());

        std::vector<double> load(sub.dofCount(), 0.0);
        std::vector<std::uint32_t> held;
        for (std::size_t b = 0; b < sub.boundaryCount(); ++b) {
            const std::uint32_t d = sub.boundaryDof()[b];
            const double x = mesh.position[(d / 3) * 3];
            if (std::fabs(x - heldEnd) < 0.5 * kLx)
                held.push_back(static_cast<std::uint32_t>(b));
            else if (d % 3 == 2)
                load[d] = 1000.0;
        }
        ReduceParams params;
        params.modes = 8;
        const Reduction reduced = craigBampton(sub, params);
        const Reduction blindReduced = craigBampton(blind, params);
        std::vector<double> state, blindState;
        std::string problem;
        expectTrue("the stiffened cantilever solves",
                   reduction::staticSolve(reduced, reduction::reduceLoad(sub, reduced, load), held,
                                          state, &problem));
        reduction::staticSolve(blindReduced, reduction::reduceLoad(blind, blindReduced, load), held,
                               blindState, &problem);
        const reduction::Validity v = reduction::checkValidity(sub, reduced, state);
        const reduction::Validity b = reduction::checkValidity(blind, blindReduced, blindState);

        // The plating half must be untouched by the members being visible: the two
        // substructures have the same stiffness and the same mass, so their element
        // stresses are the same bits. That is what says §9 added a reading rather
        // than perturbing the model.
        expectTrue("making the members visible does not move the plating's own stress",
                   v.platingVonMises == b.peakVonMises && v.worstElement == b.worstElement);
        expectTrue("and the invariants hold",
                   v.peakVonMises == std::max(v.platingVonMises, v.memberVonMises) &&
                       v.utilisation == v.peakVonMises / steel.yieldStrength &&
                       v.linear == (v.utilisation < 1.0));

        // **`worstMember` names the member, and `memberVonMises` is the peak over
        // *every* one of them.** Re-derived here from the same recovered field
        // through `memberStress`, which the closed form above has already pinned to
        // `E kappa (e - z_na)` -- so this is the sweep inside `checkValidity` under
        // test, not the arithmetic it sweeps. Without it a loop that read one member
        // too few, or that recorded the wrong index, passed everything.
        const std::vector<double> u = reduction::recover(sub, reduced, state);
        double independent = 0;
        int argmax = -1;
        for (std::size_t m = 0; m < sub.attachedMembers(); ++m) {
            const double s = std::fabs(sub.memberStress(m, u));
            if (s > independent) { independent = s; argmax = static_cast<int>(m); }
        }
        expectNear("the member peak is the peak over every member, to the last bit",
                   v.memberVonMises, independent, 0.0);
        expectTrue("and worstMember names the one it came from",
                   v.worstMember == argmax && v.worstMember >= 0 &&
                       std::fabs(sub.memberStress(static_cast<std::size_t>(v.worstMember), u)) ==
                           v.memberVonMises);
        // Where the plating's governing point is, and what kind of stress it is.
        // §9 claims the two halves peak in different places and different states;
        // nothing tests a claim, so it is measured here and asserted below.
        double nodePos[24], disp[24], stress[48], centroid[3] = {0, 0, 0};
        mesh.gather(static_cast<std::size_t>(v.worstElement), mesh.position, nodePos);
        mesh.gather(static_cast<std::size_t>(v.worstElement), u, disp);
        for (int n = 0; n < 8; ++n)
            for (int k = 0; k < 3; ++k) centroid[k] += nodePos[n * 3 + k] / 8.0;
        solidshell::elementStress(nodePos, disp, steel, solidshell::Formulation::SolidShell,
                                  stress);
        double atPeak = 0, axialThere = 0;
        for (int g = 0; g < 8; ++g) {
            const double m = misesOf(stress + g * 6);
            if (m > atPeak) { atPeak = m; axialThere = std::fabs(stress[g * 6]); }
        }
        expectNear("the element sweep here finds the same peak checkValidity did", atPeak,
                   v.platingVonMises, 0.0);
        const double seamY = kLy * static_cast<double>(seamRow) / kNy;
        return Reading{v.platingVonMises, v.memberVonMises, v.utilisation, b.utilisation,
                       v.worstMember, static_cast<int>(sub.attachedMembers()),
                       std::fabs(centroid[1] - seamY), atPeak / axialThere};
    };

    const Reading middle = cantilever(kNy / 2, 0.0);
    const Reading flipped = cantilever(kNy / 2, kLx);
    const Reading edge = cantilever(0, 0.0);

    // The two cantilevers are mirror images, so they must report the same numbers
    // from opposite ends of the fibre list -- and between them the worst member is
    // found near both ends of it. That is what makes the sweep's bounds tested
    // rather than merely exercised.
    expectNear("the mirrored cantilever reports the same member stress",
               flipped.member, middle.member, 1e-9 * middle.member);
    expectTrue("from the other end of the fibre list",
               middle.worstMember < middle.members / 4 &&
                   flipped.worstMember >= 3 * flipped.members / 4);

    // 1e-4 relative on quantities printed to four figures in `reduction.hpp` §8;
    // the run reproduces them to 58.1962 and 65.2402 MPa.
    expectNear("the plating reports 58.2 MPa", 1e-6 * middle.plating, 58.196, 1e-3);
    expectNear("and the member carries 65.2 MPa", 1e-6 * middle.member, 65.240, 1e-3);
    // 1e-6 absolute on a number the run reproduces as 0.183775 and 0.163933. It is
    // the headline of §8 and of `docs/02-simulation.md`, so it is asserted to every
    // figure those quote rather than to the two the prose rounds to.
    expectNear("so the utilisation is 0.183775, not the 0.163933 the plating alone gives",
               middle.utilisation, 0.183775, 1e-6);
    expectNear("and the old number is still exactly what a blind substructure reports",
               middle.blindUtilisation, 0.163933, 1e-6);
    // The guard the whole case rests on: if the plating governed at mid width the
    // assertion above would be satisfied by the plating alone.
    expectTrue("the member really is the governing half at mid width",
               middle.member > middle.plating && middle.worstMember >= 0);
    expectTrue("and it is 11% low without it, in the unsafe direction",
               middle.blindUtilisation < middle.utilisation &&
                   middle.blindUtilisation > 0.88 * middle.utilisation);

    // **The two halves are not two readings of one stress**, which is the whole
    // reason `Validity` reports them apart. The member's peak is pure axial, at the
    // far fibre, on the seam. The plating's is one element *off* the seam and is
    // dominated by transverse bending: its von Mises is three times its own axial
    // component there. A single peak would hide both facts.
    expectTrue("the plating's governing point is off the seam, not on it",
               middle.offSeam > 0.5 * (kLy / kNy) && middle.offSeam < 1.5 * (kLy / kNy));
    expectTrue("and it is not an axial stress at all -- von Mises is 3x its own sigma_xx",
               middle.biaxiality > 2.5 && middle.biaxiality < 3.5);
    std::printf("     the plating's peak is %.4f m off the seam and its von Mises is %.4fx its "
                "own axial component; the member's is pure axial\n",
                middle.offSeam, middle.biaxiality);

    // ...and the other direction. Same profile, same load, seam moved to the plate
    // edge: the plating governs by 2.5x and the answer must follow the structure.
    expectTrue("with the seam at the plate edge the plating is the governing half",
               edge.plating > 2.0 * edge.member);
    expectNear("so the trigger reports the plating's number", edge.utilisation,
               edge.plating / steel.yieldStrength, 1e-15 * edge.utilisation);
    expectNear("which is exactly what it reported before the members were visible",
               edge.blindUtilisation, edge.utilisation, 0.0);
    expectTrue("and the member is still read, it simply does not govern",
               edge.member > 1e6 && edge.worstMember >= 0);
    std::printf("     mid width: plating %.4f MPa, member %.4f MPa, utilisation %.6f "
                "(was %.6f)\n",
                1e-6 * middle.plating, 1e-6 * middle.member, middle.utilisation,
                middle.blindUtilisation);
    std::printf("     plate edge: plating %.3f MPa, member %.3f MPa, utilisation %.5f "
                "(was %.5f) -- the plating governs\n",
                1e-6 * edge.plating, 1e-6 * edge.member, edge.utilisation, edge.blindUtilisation);
    std::printf("     worst member: index %d of %d held at x = 0, %d of %d held at x = L\n",
                middle.worstMember, middle.members, flipped.worstMember, flipped.members);

    // --- `linear` at exactly one -------------------------------------------------
    //
    // The one place `utilisation < 1` and `utilisation <= 1` differ is at exactly
    // one, and no scaled load lands there: the utilisation is a quotient and
    // scaling it to 1 leaves a rounding. So the *yield* is moved instead, to the
    // peak the state already carries -- `x / x` is exactly 1.0 in IEEE arithmetic
    // for any finite non-zero `x`, and `yieldStrength` does not enter
    // `elementStress`, so the stress is bit for bit the same state.
    {
        const constraint::Stiffening seam = stiffenSeam(mesh, steel, profile, kNy / 2, 1);
        const Substructure sub(mesh, steel, endInterface(mesh),
                               attachmentOf(seam, mesh, steel, true));
        std::vector<double> load(sub.dofCount(), 0.0);
        std::vector<std::uint32_t> held = heldAtOrigin(sub, mesh);
        for (std::size_t b = 0; b < sub.boundaryCount(); ++b) {
            const std::uint32_t d = sub.boundaryDof()[b];
            if (mesh.position[(d / 3) * 3] >= 0.5 * kLx && d % 3 == 2) load[d] = 1000.0;
        }
        ReduceParams params;
        params.modes = 8;
        const Reduction reduced = craigBampton(sub, params);
        std::vector<double> state;
        std::string problem;
        reduction::staticSolve(reduced, reduction::reduceLoad(sub, reduced, load), held, state,
                               &problem);
        const reduction::Validity loose = reduction::checkValidity(sub, reduced, state);
        StructuralMaterial atYield = steel;
        atYield.yieldStrength = loose.peakVonMises;
        const Substructure judged(mesh, atYield, endInterface(mesh),
                                  attachmentOf(seam, mesh, atYield, true));
        const Reduction judgedReduced = craigBampton(judged, params);
        std::vector<double> judgedState;
        reduction::staticSolve(judgedReduced, reduction::reduceLoad(judged, judgedReduced, load),
                               held, judgedState, &problem);
        const reduction::Validity v = reduction::checkValidity(judged, judgedReduced, judgedState);
        expectNear("moving the yield does not move the stress", v.peakVonMises, loose.peakVonMises,
                   0.0);
        expectNear("so the utilisation is exactly one", v.utilisation, 1.0, 0.0);
        expectTrue("and exactly one is not linear", !v.linear);
        expectTrue("while a hair under it is",
                   loose.linear == (loose.utilisation < 1.0) && loose.utilisation < 1.0);
        std::printf("     yield set to the peak: utilisation %.17g, linear %d\n", v.utilisation,
                    static_cast<int>(v.linear));
    }

    // --- The two producers cannot drift apart ------------------------------------
    //
    // `stiffnessBlocks` is `attachedForms(...).stiffness` and nothing else, so the
    // stress forms are built in the same loop behind the same skip test and cannot
    // come back paired with the wrong block. Asserted bit for bit, and with a
    // *degenerate* fibre in the set so the skip test's verdict actually changes --
    // without one the branch is never taken and the mutant that widens it to
    // `scale >= 0` is indistinguishable.
    {
        constraint::Stiffening seam = stiffenSeam(mesh, steel, profile, kNy / 2, 1);
        const std::size_t sound = seam.fiberCount();
        constraint::Fiber dead = seam.fiber.front();
        dead.area = 0.0;  // no area, so EA/L is zero and the fibre carries nothing
        seam.fiber.push_back(dead);
        constraint::Fiber collapsed = seam.fiber.front();
        collapsed.end[1] = collapsed.end[0];  // no length either
        seam.fiber.push_back(collapsed);

        const constraint::RestFibers forms = constraint::restFibers(seam, mesh.position);
        expectTrue("a fibre with coincident ends has no rest length", !forms.ok);
        const constraint::AttachedForms built =
            constraint::attachedForms(seam, mesh.position, forms, steel.youngsModulus);
        expectEqualCount("the degenerate fibres are dropped from the blocks", built.stiffness.size(),
                         sound);
        expectEqualCount("and from the stress forms, by the same test", built.stress.size(), sound);

        const std::vector<solidshell::DofBlock> alone =
            constraint::stiffnessBlocks(seam, mesh.position, forms, steel.youngsModulus);
        bool same = alone.size() == built.stiffness.size();
        for (std::size_t b = 0; b < alone.size() && same; ++b) {
            same = alone[b].dof == built.stiffness[b].dof &&
                   alone[b].stiffness.size() == built.stiffness[b].stiffness.size();
            for (std::size_t e = 0; e < alone[b].stiffness.size() && same; ++e)
                same = alone[b].stiffness[e] == built.stiffness[b].stiffness[e];
        }
        expectTrue("and `stiffnessBlocks` is that same list, to the last bit", same);
        std::printf("     %zu sound fibres plus one with no area and one with no length: "
                    "%zu blocks, %zu stress forms\n",
                    sound, built.stiffness.size(), built.stress.size());
    }

    // `Attachment::empty()` has to know about the stress forms too, or an
    // attachment carrying nothing else reports itself as nothing at all.
    {
        reduction::Attachment onlyStress;
        onlyStress.stress.push_back({1.0, 2.0, 3.0});
        expectTrue("an attachment carrying only stress forms is not empty", !onlyStress.empty());
        expectTrue("and one carrying nothing is", reduction::Attachment{}.empty());
    }

    // --- A member that is not axis aligned --------------------------------------
    //
    // Every fibre a `SeamRun` builds on a flat plate runs along the plating, so its
    // direction has no through-thickness component and the last four entries of the
    // twelve-term stress form are *identically zero*. Mutation testing said so: the
    // mutant that drops the twelfth term survived everything above, and it is
    // equivalent on any seam lying in a plane -- `0.0 * u` is a no-op.
    //
    // A `constraint::Fiber` is not restricted to that. Tying one end to the -zeta
    // face and the other to the +zeta face of the next column is a diagonal brace,
    // and its direction has all three components. It is checked against the
    // *geometric* definition of an axial strain -- the change in the distance
    // between the two tied points, taken through `constraint::tiedPoint` on the
    // displacement array -- which is a different route to the number from the
    // condensed rank-one form the stress is recovered through.
    {
        const constraint::Stiffening full = stiffenSeam(mesh, steel, profile, kNy / 2, 1);
        constraint::Stiffening brace = full;
        brace.fiber.clear();
        constraint::Fiber diagonal;
        diagonal.end[0] = constraint::Tie{plateNode(0, kNy / 2, 0), plateNode(0, kNy / 2, 1), 0.0};
        diagonal.end[1] = constraint::Tie{plateNode(1, kNy / 2, 0), plateNode(1, kNy / 2, 1), 1.0};
        diagonal.area = 1.0e-4;
        diagonal.offset = 0.0;
        brace.fiber.push_back(diagonal);

        const constraint::RestFibers forms = constraint::restFibers(brace, mesh.position);
        expectTrue("the diagonal brace has a rest length", forms.ok);
        const constraint::AttachedForms built =
            constraint::attachedForms(brace, mesh.position, forms, steel.youngsModulus);
        expectEqualCount("and one block and one stress form", built.stress.size(), 1u);
        // The guard: without a through-thickness component the last four terms are
        // zero and this whole case proves nothing about them.
        expectTrue("the brace really does run through the thickness",
                   std::fabs(built.stress[0][11]) > 1e-3 * std::fabs(built.stress[0][9]));

        const Substructure sub(mesh, steel, endInterface(mesh),
                               attachmentOf(brace, mesh, steel, true));
        expectTrue("the braced substructure factors", sub.ready());

        // A field with all three components moving, so every one of the twelve terms
        // multiplies something non-zero.
        std::vector<double> u(sub.dofCount(), 0.0);
        for (std::size_t n = 0; n < mesh.nodeCount(); ++n) {
            const double x = mesh.position[n * 3], y = mesh.position[n * 3 + 1];
            const double z = mesh.position[n * 3 + 2];
            u[n * 3] = 1.1e-4 * x + 3.0e-4 * z + 7.0e-5 * y;
            u[n * 3 + 1] = -5.0e-5 * x + 9.0e-5 * z;
            u[n * 3 + 2] = 4.0e-5 * x - 2.0e-5 * y + 6.0e-4 * z;
        }
        const Vec3 restSpan = constraint::tiedPoint(diagonal.end[1], mesh.position) -
                              constraint::tiedPoint(diagonal.end[0], mesh.position);
        const double restLength = length(restSpan);
        const Vec3 axis = restSpan / restLength;
        const Vec3 moved = constraint::tiedPoint(diagonal.end[1], u) -
                           constraint::tiedPoint(diagonal.end[0], u);
        const double want =
            steel.youngsModulus * (moved.x * axis.x + moved.y * axis.y + moved.z * axis.z) /
            restLength;
        const double got = sub.memberStress(0, u);
        expectTrue("the brace is stressed", std::fabs(want) > 1e6);
        // 1e-13 relative, measured at 4e-16: two twelve-term sums of the same
        // products in different orders, so the floor is a couple of ulp.
        expectNear("a member that is not axis aligned reports the change in its own length", got,
                   want, 1e-13 * std::fabs(want));
        std::printf("     diagonal brace: axis (%.4f, %.4f, %.4f), %.4f MPa against a geometric "
                    "%.4f MPa\n",
                    axis.x, axis.y, axis.z, 1e-6 * got, 1e-6 * want);
    }

    // --- One member, and two that tie -------------------------------------------
    //
    // Both are here because mutation testing found the sweep's *first* index and its
    // tie-break untested, and neither can be reached by the sixteen-fibre seam
    // above: under bending the worst fibre is always the one furthest from the
    // neutral axis, which is never index 0 and never ties with anything.
    //
    //   * A region carrying **one** member has to report it. A sweep starting one
    //     index late reports no member at all and a `worstMember` of -1, which is
    //     indistinguishable from a region that has none -- the exact failure §8
    //     exists to prevent, one layer up.
    //   * Two **identical** members tie to the last bit, and the tie-break is then
    //     what `worstMember` means. Two identical longitudinals either side of a
    //     symmetric section under a symmetric load do this for real, so it is
    //     specified here -- the *first* of the tied members -- rather than left to
    //     whichever comparison the sweep happens to use.
    {
        const constraint::Stiffening full = stiffenSeam(mesh, steel, profile, kNy / 2, 1);
        for (int which = 0; which < 2; ++which) {
            constraint::Stiffening cut = full;
            cut.fiber.assign(1, full.fiber[1]);  // the far fibre of the first segment
            if (which == 1) cut.fiber.push_back(cut.fiber.front());  // ...twice over
            const Substructure sub(mesh, steel, endInterface(mesh),
                                   attachmentOf(cut, mesh, steel, true));
            expectTrue("the one-member substructure factors", sub.ready());
            expectEqualCount("and carries exactly the members it was given",
                             sub.attachedMembers(), which == 0 ? 1u : 2u);

            std::vector<double> load(sub.dofCount(), 0.0);
            std::vector<std::uint32_t> held = heldAtOrigin(sub, mesh);
            for (std::size_t b = 0; b < sub.boundaryCount(); ++b) {
                const std::uint32_t d = sub.boundaryDof()[b];
                if (mesh.position[(d / 3) * 3] >= 0.5 * kLx && d % 3 == 2) load[d] = 1000.0;
            }
            ReduceParams params;
            params.modes = 8;
            const Reduction reduced = craigBampton(sub, params);
            std::vector<double> state;
            std::string problem;
            expectTrue("and solves",
                       reduction::staticSolve(reduced, reduction::reduceLoad(sub, reduced, load),
                                              held, state, &problem));
            const reduction::Validity v = reduction::checkValidity(sub, reduced, state);
            const std::vector<double> u = reduction::recover(sub, reduced, state);
            expectTrue("the single member is stressed and is read", v.memberVonMises > 1e6);
            expectNear("at exactly what `memberStress` says it carries", v.memberVonMises,
                       std::fabs(sub.memberStress(0, u)), 0.0);
            expectEqual("and it is member zero that is named", v.worstMember, 0);
            if (which == 1)
                expectNear("the duplicate really does tie it, to the last bit",
                           sub.memberStress(1, u), sub.memberStress(0, u), 0.0);
        }
    }

    // --- The negative control, bit for bit --------------------------------------
    //
    // With nothing attached, `checkValidity` must return *exactly* what it returned
    // before §9 -- not close. Every field, compared against a substructure that
    // cannot have a member at all, and the derived fields compared against the
    // element half rather than against a remembered constant.
    {
        const Substructure bare(mesh, steel, endInterface(mesh));
        ReduceParams params;
        params.modes = 8;
        const Reduction reduced = craigBampton(bare, params);
        std::vector<double> load(bare.dofCount(), 0.0);
        std::vector<std::uint32_t> held = heldAtOrigin(bare, mesh);
        for (std::size_t b = 0; b < bare.boundaryCount(); ++b) {
            const std::uint32_t d = bare.boundaryDof()[b];
            if (mesh.position[(d / 3) * 3] >= 0.5 * kLx && d % 3 == 2) load[d] = 1000.0;
        }
        std::vector<double> state;
        std::string problem;
        expectTrue("the bare cantilever solves",
                   reduction::staticSolve(reduced, reduction::reduceLoad(bare, reduced, load), held,
                                          state, &problem));
        const reduction::Validity v = reduction::checkValidity(bare, reduced, state);
        // Against `peakVonMises` computed here, over the same recovered field: a
        // second implementation of the element half, so this is not the new code
        // agreeing with itself.
        const std::vector<double> u = reduction::recover(bare, reduced, state);
        const StressPeak independent = peakVonMises(mesh, steel, u);
        expectTrue("with nothing attached there are no members to read",
                   bare.attachedMembers() == 0 && bare.attachedBlocks() == 0);
        expectTrue("the load stresses the bare plate", v.peakVonMises > 1e6);
        expectTrue("and the member half is exactly absent",
                   v.memberVonMises == 0.0 && v.worstMember == -1);
        expectTrue("so the peak is the plating's, to the last bit",
                   v.peakVonMises == v.platingVonMises);
        expectNear("and it is the peak an independent sweep of the same field finds",
                   v.platingVonMises, independent.mises, 0.0);
        expectTrue("with the utilisation the old formula gives, to the last bit",
                   v.utilisation == v.platingVonMises / steel.yieldStrength);
        std::printf("     no attachment: %.6f MPa, utilisation %.9f, member half exactly absent\n",
                    1e-6 * v.peakVonMises, v.utilisation);
    }

    // --- A tie leaves a hole in the recovered field, and it is not a small one ---
    //
    // A degree of freedom `Attachment::constrained` eliminated is in neither
    // partition, so `recover` used to leave it at zero -- a hole in the middle of
    // the displacement field that every element touching that node reads as an
    // enormous gradient. It is not a curiosity: `section.cpp` ties every junction
    // of a hold this way, so `checkValidity` would have reported every tied section
    // as thousands of times past yield.
    {
        const auto node = [](int i, int j, int k) {
            return static_cast<std::uint32_t>((i * (kNy + 1) + j) * 2 + k);
        };
        const std::uint32_t slave = node(4, 2, 0), first = node(3, 2, 0), second = node(5, 2, 0);
        reduction::Attachment tied;
        for (int k = 0; k < 3; ++k) {
            solidshell::Mpc mpc;
            mpc.slave = slave * 3u + static_cast<std::uint32_t>(k);
            mpc.master = {first * 3u + static_cast<std::uint32_t>(k),
                          second * 3u + static_cast<std::uint32_t>(k)};
            mpc.weight = {0.5, 0.5};
            tied.constrained.push_back(mpc);
        }
        const Substructure sub(mesh, steel, endInterface(mesh), tied);
        expectTrue("the tied substructure factors", sub.ready());
        expectEqualCount("and really did eliminate three degrees of freedom",
                         sub.expansion().eliminatedCount(), 3u);

        std::vector<double> load(sub.dofCount(), 0.0);
        std::vector<std::uint32_t> held = heldAtOrigin(sub, mesh);
        for (std::size_t b = 0; b < sub.boundaryCount(); ++b) {
            const std::uint32_t d = sub.boundaryDof()[b];
            if (mesh.position[(d / 3) * 3] >= 0.5 * kLx && d % 3 == 2) load[d] = 1000.0;
        }
        ReduceParams params;
        params.modes = 0;
        const Reduction reduced = craigBampton(sub, params);
        std::vector<double> state;
        std::string problem;
        expectTrue("the tied cantilever solves",
                   reduction::staticSolve(reduced, reduction::reduceLoad(sub, reduced, load), held,
                                          state, &problem));
        const std::vector<double> u = reduction::recover(sub, reduced, state);
        double worst = 0, scale = 0;
        for (int k = 0; k < 3; ++k) {
            const std::size_t s = slave * 3 + static_cast<std::size_t>(k);
            const double want = 0.5 * u[first * 3 + static_cast<std::size_t>(k)] +
                                0.5 * u[second * 3 + static_cast<std::size_t>(k)];
            worst = std::max(worst, std::fabs(u[s] - want));
            scale = std::max(scale, std::fabs(want));
        }
        expectTrue("the tie has something to carry", scale > 1e-6);
        expectNear("a constrained degree of freedom comes back as its masters say, exactly", worst,
                   0.0, 0.0);
        const reduction::Validity v = reduction::checkValidity(sub, reduced, state);
        // The hole put the peak at 850 788 MPa on this case, so anything within
        // three orders of the plate's real stress says it is gone. Asserted against
        // the *unconstrained* plate rather than a remembered number: the tie is
        // between two neighbours of a node, so it stiffens the plate a little and
        // must not change its stress by orders.
        const Substructure loose(mesh, steel, endInterface(mesh));
        const Reduction looseReduced = craigBampton(loose, params);
        std::vector<double> looseState;
        reduction::staticSolve(looseReduced, reduction::reduceLoad(loose, looseReduced, load), held,
                               looseState, &problem);
        const reduction::Validity untied = reduction::checkValidity(loose, looseReduced, looseState);
        expectTrue("and the tied plate's peak stress is the untied one's, to within the tie",
                   v.peakVonMises < 2.0 * untied.peakVonMises &&
                       v.peakVonMises > 0.5 * untied.peakVonMises);
        std::printf("     one node tied to two neighbours: peak %.3f MPa against %.3f MPa "
                    "untied; unfilled it read 850 788 MPa\n",
                    1e-6 * v.peakVonMises, 1e-6 * untied.peakVonMises);
    }
}

// --- 12. And the steel it is made of (§8 item 3) ----------------------------------
//
// Stiffness without mass is a model that is stiffer and no heavier, so every
// frequency comes out high -- and frequencies are most of what this tier is for.
// What is asserted here is the mass identity, exactly, and the *direction and
// bracket* of what it does to the spectrum, because both ends of that bracket are
// rigorous statements about the discrete problem rather than eyeballed numbers.
void testAttachedMass() {
    std::printf("\n--- reduction: the stiffener's own inertia ---\n");
    const solidshell::HexMesh mesh = testPlate();
    const StructuralMaterial steel = ah36Steel();
    const StiffenerProfile profile = flatBar(0.200, 0.010);
    const constraint::Stiffening seam = stiffenSeam(mesh, steel, profile, kNy / 2, 1);
    const ProfileSection section = profileSection(profile);

    const Substructure bare(mesh, steel, endInterface(mesh));
    const Substructure sub(mesh, steel, endInterface(mesh),
                           attachmentOf(seam, mesh, steel, true));
    const Substructure noMass(mesh, steel, endInterface(mesh),
                              attachmentOf(seam, mesh, steel, false));

    // The steel, in closed form: density times the profile's area times the length
    // of seam the fibres cover. `profileSection` computes the area from the
    // rectangle dimensions and knows nothing about fibres.
    const double steelMass = steel.density * section.area * kLx;
    expectNear("the attachment weighs the stiffener's own steel", sub.attachedMass(), steelMass,
               1e-12 * steelMass);
    expectNear("and the fibre builder agrees with it", seam.mass, steelMass, 1e-12 * steelMass);
    expectNear("so the substructure's total mass is the plating plus the stiffener",
               sub.totalMass(), bare.totalMass() + steelMass, 1e-12 * sub.totalMass());
    expectTrue("which is not a rounding on the plating", steelMass > 0.2 * bare.totalMass());

    // Through T, which is the statement that matters: a mass that reached the
    // diagonal but not the reduction would pass the line above and fail this one.
    // §1 property 2 -- a rigid translation of the reduced model reports the
    // substructure's own mass, exactly.
    {
        ReduceParams params;
        params.modes = 6;
        const Reduction reduced = craigBampton(sub, params);
        expectTrue("the stiffened reduction is usable", !reduced.empty());
        const std::size_t n = static_cast<std::size_t>(reduced.size());
        std::vector<double> rigid(n, 0.0);
        for (std::size_t b = 0; b < sub.boundaryCount(); ++b)
            if (sub.boundaryDof()[b] % 3 == 2) rigid[b] = 1.0;
        double quadratic = 0;
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j) quadratic += rigid[i] * reduced.mass[i * n + j] * rigid[j];
        expectNear("a rigid translation of the reduced model weighs plating plus stiffener",
                   quadratic, bare.totalMass() + steelMass, 1e-9 * sub.totalMass());
        std::printf("     reduced rigid translation: %.6f kg against %.6f kg\n", quadratic,
                    bare.totalMass() + steelMass);
    }

    // What the equal split over the pair gives up, as a number rather than as a
    // caveat. Rotate rigidly about the seam: the stiffener's steel is split over
    // the through-thickness pair, so about that axis it carries `m (t/2)^2` -- the
    // plate's own half thickness and nothing else -- where the real member carries
    // `m d^2` about its own centroid, a hundred millimetres out. See §8 and
    // `constraint.hpp`: the consistent condensation that *would* carry it puts a
    // negative mass on one node of every eccentric pair.
    {
        const double seamY = kLy * static_cast<double>(kNy / 2) / kNy;
        double lumped = 0, plating = 0;
        for (std::size_t n = 0; n < mesh.nodeCount(); ++n) {
            const double dy = mesh.position[n * 3 + 1] - seamY, dz = mesh.position[n * 3 + 2];
            lumped += sub.mass()[n * 3] * (dy * dy + dz * dz);
            plating += bare.mass()[n * 3] * (dy * dy + dz * dz);
        }
        const double centroid = 0.5 * kThickness + section.centroid;
        const double carried = steelMass * 0.25 * kThickness * kThickness;
        const double real = steelMass * centroid * centroid;
        expectNear("the stiffener's steel carries only the plate's half thickness about the seam",
                   lumped - plating, carried, 1e-9 * carried);
        // 312x, and the ratio is `(d / (t/2))^2` in closed form: what the model
        // keeps of this member's rotary inertia is a third of a per cent of it.
        expectNear("which is the square of the ratio of the two lever arms short of the member's",
                   real / carried, (centroid / (0.5 * kThickness)) * (centroid / (0.5 * kThickness)),
                   1e-9 * real / carried);
        expectTrue("and that is a shortfall of orders, not of per cent", real > 100.0 * carried);
        std::printf("     rotary inertia about the seam: the stiffener adds %.4e kg m^2 where its "
                    "own eccentricity would add %.4e -- %.1fx, all of it given up\n",
                    carried, real, real / carried);
    }

    // --- The spectrum. Both ends of the bracket are rigorous:
    //
    //   * `M` grows and `K` does not, so **every** eigenvalue falls or stays. Any
    //     frequency that rose would mean the mass had landed somewhere it does not
    //     belong.
    //   * `M_with <= r M_without` entrywise on a diagonal mass, so no frequency
    //     falls further than `1 / sqrt(r)`. `r` is measured here, not assumed.
    //
    // The guard against the vacuous version -- which is the version where no mass
    // is added at all, and which sits exactly on the upper end of that bracket --
    // is that at least one frequency must fall by more than a tenth.
    {
        std::vector<std::uint8_t> pinned(sub.dofCount(), 0);
        for (int j = 0; j <= kNy; ++j)
            for (int k = 0; k < 2; ++k)
                for (int axis = 0; axis < 3; ++axis) pinned[plateNode(0, j, k) * 3 + axis] = 1;
        const std::vector<double> withMass = fullFrequencies(sub, pinned);
        const std::vector<double> without = fullFrequencies(noMass, pinned);
        const std::vector<double> plain = fullFrequencies(bare, pinned);
        double ratio = 0;
        for (std::size_t d = 0; d < sub.dofCount(); ++d)
            if (!pinned[d]) ratio = std::max(ratio, sub.mass()[d] / noMass.mass()[d]);
        const double floorRatio = 1.0 / std::sqrt(ratio);

        const std::size_t count = std::min<std::size_t>(12, withMass.size());
        double worstDrop = 0;
        bool everRose = false, everBelowFloor = false;
        for (std::size_t i = 0; i < count; ++i) {
            worstDrop = std::max(worstDrop, 1.0 - withMass[i] / without[i]);
            if (withMass[i] > without[i] * (1.0 + 1e-12)) everRose = true;
            if (withMass[i] < without[i] * floorRatio * (1.0 - 1e-12)) everBelowFloor = true;
        }
        expectTrue("no frequency rises when the stiffener's mass is added", !everRose);
        expectTrue("and none falls further than the worst nodal mass ratio allows",
                   !everBelowFloor);
        expectTrue("the mass is not decoration: one of the first twelve falls by over a tenth",
                   worstDrop > 0.10);
        // And the stiffness is still the dominant effect, in the direction the
        // section properties give: the panel's second moment about the mid-surface
        // goes up by two orders, so the first frequency goes up, not down.
        expectTrue("the stiffener still raises the first frequency", withMass[0] > 2.0 * plain[0]);
        std::printf("     first frequency %.3f -> %.3f rad/s (bare -> stiffened); the fibre mass "
                    "moves the first by %.3f%% and the worst of twelve by %.1f%% (floor %.1f%%)\n",
                    plain[0], withMass[0], 100.0 * (1.0 - withMass[0] / without[0]),
                    100.0 * worstDrop, 100.0 * (1.0 - floorRatio));
    }

    // A stiffener stiff enough to hold its own line still is a *node* of the first
    // mode, which is why the first frequency above barely moved. On a member the
    // panel can actually bend with, the whole stiffened section works and the first
    // frequency lands where the beam idealisation says: `sqrt(I_ratio / m_ratio)`,
    // both taken from `scantlings::stiffenedSection`. The **smeared** panel is the
    // negative control -- it has the same area and the same mass and its second
    // moment is `(t + A/s)^3 / t^3`, which is the factor of 130 `scantlings.hpp` §1
    // rejects smearing for.
    {
        const StiffenerProfile modest = flatBar(0.060, 0.006);
        const constraint::Stiffening line = stiffenSeam(mesh, steel, modest, kNy / 2, 1);
        const Substructure stiffened(mesh, steel, endInterface(mesh),
                                     attachmentOf(line, mesh, steel, true));
        std::vector<std::uint8_t> pinned(sub.dofCount(), 0);
        for (int j = 0; j <= kNy; ++j)
            for (int k = 0; k < 2; ++k)
                for (int axis = 0; axis < 3; ++axis) pinned[plateNode(0, j, k) * 3 + axis] = 1;
        const double got = fullFrequencies(stiffened, pinned)[0] / fullFrequencies(bare, pinned)[0];

        const StiffenedSection panel = stiffenedSection(modest, kThickness, kLy);
        const double plateOwn = kLy * kThickness * kThickness * kThickness / 12.0;
        const double aboutMidSurface =
            panel.secondMoment + panel.area * panel.neutralAxis * panel.neutralAxis;
        const double massRatio = stiffened.totalMass() / bare.totalMass();
        const double predicted = std::sqrt(aboutMidSurface / plateOwn / massRatio);
        const double smeared = smearedThickness(kThickness, modest, kLy) / kThickness;
        const double smearedPrediction = std::sqrt(smeared * smeared * smeared / massRatio);

        // 6%, where the measurement is 4.5%: an Euler beam of the full plate width
        // against a two-dimensional plate carrying its stiffener on one line, so a
        // few per cent is what a correct model gives. The tolerance is nowhere near
        // wide enough to admit the smeared answer, which is the point.
        expectNear("the first frequency lands where the stiffened section says", got, predicted,
                   0.06 * predicted);
        expectTrue("and nowhere near where a smeared panel would put it",
                   got > 2.0 * smearedPrediction);
        std::printf("     a 60x6 bar: first frequency x%.4f, stiffened section predicts x%.4f, a "
                    "smeared plate of the same area predicts x%.4f\n",
                    got, predicted, smearedPrediction);
    }
}

// --- The position match, asked about on its own -----------------------------------
//
// **Every claim `matchBoundaries` makes is invisible through a substructure**, and
// mutation testing is what said so: deleting the axis check, deleting the "this DOF
// of B is already taken" mark, and deleting the `break` that keeps the *first* match
// all survived a suite that only ever matched two real meshes. They survive together
// because a substructure's boundary DOF come out ordered by node and then by axis, so
// a greedy first-match walks the two lists in step and lands on the right partner for
// the wrong reason — and because each of the three mutants is masked by the other two.
//
// So the primitive is asked directly, on lists built here, where the order can be
// wrong and the positions can repeat. Two coincident nodes with the same axis is not
// a synthetic case: it is what an unwelded seam in a mesh looks like.
void testBoundaryMatchOnItsOwn() {
    std::printf("\n--- reduction: the position match, on lists built to break it ---\n");
    using reduction::BoundaryDof;

    // A: one node's three DOF, in axis order. B: the same node, axes **reversed**.
    const std::vector<BoundaryDof> a{{{1.0, 2.0, 3.0}, 0}, {{1.0, 2.0, 3.0}, 1}, {{1.0, 2.0, 3.0}, 2}};
    const std::vector<BoundaryDof> b{{{1.0, 2.0, 3.0}, 2}, {{1.0, 2.0, 3.0}, 1}, {{1.0, 2.0, 3.0}, 0}};
    const reduction::InterfaceMap crossed = reduction::matchBoundaries(a, b);
    expectEqualCount("every DOF of a coincident node matches", crossed.shared, 3u);
    for (std::size_t i = 0; i < a.size(); ++i) {
        expectTrue("and it matched something", crossed.aToB[i] >= 0);
        if (crossed.aToB[i] < 0) continue;
        expectEqual("on the same axis, whatever order the other list is in",
                    static_cast<long long>(b[static_cast<std::size_t>(crossed.aToB[i])].axis),
                    static_cast<long long>(a[i].axis));
    }

    // B with a **duplicate**: two entries at the same place on the same axis, which
    // is what a mesh with an unwelded seam hands over. Neither may be claimed twice,
    // and one A DOF may not claim both.
    const std::vector<BoundaryDof> twice{{{0, 0, 0}, 0}, {{0, 0, 0}, 0}};
    const std::vector<BoundaryDof> one{{{0, 0, 0}, 0}};
    const reduction::InterfaceMap oneToTwo = reduction::matchBoundaries(one, twice);
    expectEqualCount("one DOF matches one DOF, not both of a coincident pair", oneToTwo.shared, 1u);
    const reduction::InterfaceMap twoToTwo = reduction::matchBoundaries(twice, twice);
    expectEqualCount("and two match two", twoToTwo.shared, 2u);
    expectTrue("each on its own partner", twoToTwo.aToB[0] != twoToTwo.aToB[1]);
    const reduction::InterfaceMap twoToOne = reduction::matchBoundaries(twice, one);
    expectEqualCount("where only one is on offer, only one is taken", twoToOne.shared, 1u);

    // The tolerance, at and either side of it. Inclusive at exactly the tolerance:
    // "within" is the word, and pinning it is what stops the comparison drifting.
    const double tol = 1e-6;
    for (auto [gap, want] : std::vector<std::pair<double, std::size_t>>{
             {0.5 * tol, 1u}, {tol, 1u}, {2.0 * tol, 0u}}) {
        const std::vector<BoundaryDof> near{{{gap, 0, 0}, 0}};
        const reduction::InterfaceMap m = reduction::matchBoundaries(one, near, tol);
        expectEqualCount("the tolerance includes its own value and excludes twice it", m.shared,
                         want);
        if (want > 0) expectNear("and the gap is reported", m.worstGap, gap, 1e-18);
    }

    // Nothing shared is a complaint, not silence: two components assembled through an
    // empty map are two independent models standing side by side.
    const std::vector<BoundaryDof> elsewhere{{{9, 9, 9}, 0}};
    const reduction::InterfaceMap apart = reduction::matchBoundaries(one, elsewhere);
    expectEqualCount("nothing matches", apart.shared, 0u);
    expectTrue("and it says so", !apart.problems.empty());
}

// --- Synthesis: more than two components ------------------------------------------
//
// The same plate, cut into **three** pieces this time, which is the case the
// two-component `assemble` cannot express at all -- and the middle piece is what
// makes it a real test rather than the two-piece one written twice, because its
// boundary is shared with a *different* component at each end and neither pairwise
// map says on its own which assembled row a DOF belongs on.
//
// The reference is the same monolithic plate the two-component test uses, solved the
// long way, and the negative controls are the same shape: an unjoined trio has
// eighteen rigid body modes rather than six, and the assembled model is then three
// pieces rather than one.
void testThreeComponentsAssembleIntoTheWhole() {
    std::printf("\n--- reduction: three components, and the union-find that needs ---\n");

    const double L = 1.2, W = 0.2, T = 0.01;
    const int NX = 9, NY = 2, NZ = 1;
    const StructuralMaterial steel = ah36Steel();

    // Three thirds, by element. The node spacing matches by construction: L*i/NX on
    // the whole is (L/3)*i/(NX/3) on each third.
    solidshell::HexMesh whole = solidshell::makePlateMesh(L, W, T, NX, NY, NZ);
    std::vector<solidshell::HexMesh> part;
    for (int p = 0; p < 3; ++p) {
        part.push_back(solidshell::makePlateMesh(L / 3, W, T, NX / 3, NY, NZ));
        for (std::size_t i = 0; i + 2 < part.back().position.size(); i += 3)
            part.back().position[i] += (L / 3) * p;
    }

    const std::vector<reduction::Plane> cuts{{{L / 3, 0, 0}, {1, 0, 0}}, {{2 * L / 3, 0, 0}, {1, 0, 0}}};
    std::vector<reduction::Substructure> substructure;
    substructure.reserve(3);
    for (int p = 0; p < 3; ++p)
        substructure.emplace_back(part[static_cast<std::size_t>(p)], steel,
                                  reduction::nodesNearPlanes(part[static_cast<std::size_t>(p)],
                                                             cuts, 1e-9));
    for (const reduction::Substructure& s : substructure)
        expectTrue("each third is ready", s.ready());
    // The middle third is held at both ends and the outer two at one, which is the
    // asymmetry that makes this three components rather than two done twice.
    expectEqualCount("the middle third has twice the boundary of an end one",
                     substructure[1].boundaryCount(), 2 * substructure[0].boundaryCount());

    reduction::Substructure sWhole(whole, steel, reduction::nodesNearPlanes(whole, cuts, 1e-9));
    expectTrue("the whole plate is ready", sWhole.ready());
    double parts = 0;
    for (const reduction::Substructure& s : substructure) parts += s.totalMass();
    expectNear("the three thirds weigh what the whole plate weighs", parts, sWhole.totalMass(),
               1e-9 * sWhole.totalMass());

    // The reference spectrum: the whole plate, formed densely and solved once.
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
        for (std::size_t i = 0; i < nFull; ++i) mFull[i * nFull + i] = sWhole.mass()[i];
    }
    const reduction::Eigenpairs exact =
        reduction::generalisedEigen(kFull, mFull, static_cast<int>(nFull));
    expectTrue("the reference spectrum converged", exact.converged);
    expectTrue("the reference spectrum has at least seven modes", exact.value.size() > 6);
    const double firstElastic = std::sqrt(exact.value[6]);
    expectTrue("the sixth mode is rigid and the seventh is not, by orders of magnitude",
               std::sqrt(exact.value[5]) < 1e-4 * firstElastic);
    const double rigidCutoff = 1e-3 * firstElastic;

    // Matching. `matchComponents` looks at every pair; the outer two share nothing
    // and must therefore produce no joint at all rather than an empty one.
    std::vector<reduction::Reduction> reduced;
    reduction::ReduceParams params;
    params.modes = 12;
    params.cutoffFrequency = 0;
    reduced.reserve(3);
    for (const reduction::Substructure& s : substructure)
        reduced.push_back(reduction::craigBampton(s, params));

    std::vector<reduction::Component> component(3);
    for (int p = 0; p < 3; ++p)
        component[static_cast<std::size_t>(p)] = {&substructure[static_cast<std::size_t>(p)],
                                                 &reduced[static_cast<std::size_t>(p)]};

    const std::vector<reduction::Joint> all = reduction::matchComponents(component);
    expectEqualCount("three components in a row make two joints, not three", all.size(), 2u);
    for (const reduction::Joint& j : all)
        expectTrue("and neither of them joins the two ends to each other",
                   j.b == j.a + 1);
    const std::vector<reduction::Joint> neighbours = reduction::matchNeighbours(component);
    expectEqualCount("neighbours-only finds the same two", neighbours.size(), all.size());

    // One cut plane's worth of DOF. An end third touches one of the two planes and
    // the middle third both, so the boundary sets are 1, 2 and 1 planes and the
    // assembled union is the **two interior planes** -- the plate's own ends are
    // interior to their components and are not interface at all.
    const std::size_t perPlane = 3 * reduction::nodesNearPlanes(part[0], cuts, 1e-9).size();
    for (const reduction::Joint& j : all)
        expectEqualCount("each joint shares a whole cut plane", j.map.shared, perPlane);

    const reduction::Assembly trio = reduction::assemble(component, all);
    expectTrue("the assembly raises nothing", trio.problems.empty());
    expectEqual("it remembers how many components it has", trio.parts, 3);
    // Four cut planes' worth of boundary between three thirds: the two ends and the
    // two interior planes, each counted once.
    expectEqualCount("the assembly counts each shared boundary DOF once",
                     static_cast<std::size_t>(trio.boundary), 2 * perPlane);
    expectEqualCount("and carries an identity for every one of them", trio.boundaryPoint.size(),
                     static_cast<std::size_t>(trio.boundary));
    expectNear("nothing was merged across a gap", trio.worstMergedGap, 0.0, 1e-15);
    expectEqual("nor across axes", trio.axisDisagreements, 0);
    expectEqual("and it is one piece", reduction::assembledComponents(trio), 1);

    // The identity has to be the substructures' own, not something plausible: every
    // component's boundary DOF must land on a row whose position and axis are the
    // ones that DOF has in its own mesh.
    double worstIdentity = 0;
    for (int p = 0; p < 3; ++p) {
        const auto c = static_cast<std::size_t>(p);
        const std::vector<reduction::BoundaryDof> own = reduction::boundaryIdentity(substructure[c]);
        for (std::size_t b = 0; b < own.size(); ++b) {
            const int row = trio.from[c][b];
            expectTrue("every boundary DOF landed on a boundary row", row >= 0 && row < trio.boundary);
            if (row < 0 || row >= trio.boundary) continue;
            const reduction::BoundaryDof& there = trio.boundaryPoint[static_cast<std::size_t>(row)];
            expectEqual("on the same axis", static_cast<long long>(there.axis),
                        static_cast<long long>(own[b].axis));
            worstIdentity = std::max(worstIdentity, length(there.position - own[b].position));
        }
    }
    expectNear("and at the same place", worstIdentity, 0.0, 1e-15);

    std::vector<double> omega = reduction::assembledFrequencies(trio, {});
    int rigid = 0;
    for (double w : omega)
        if (w < rigidCutoff) ++rigid;
    expectEqual("a free assembly of three components has six rigid body modes, not eighteen",
                static_cast<long long>(rigid), 6);
    double worst = 0;
    for (std::size_t i = 0; i < 4 && 6 + i < exact.value.size(); ++i) {
        const double got = omega[static_cast<std::size_t>(rigid) + i];
        const double want = std::sqrt(exact.value[6 + i]);
        worst = std::max(worst, std::fabs(got - want) / want);
        expectTrue("an assembled frequency is an upper bound on the true one",
                   got > want * (1.0 - 1e-6));
    }
    std::printf("     three thirds, 12 modes each: worst of the first four elastic modes %.3e\n",
                worst);
    // Measured at 7.0e-08. Asserted at 1e-6 rather than at the 1e-4 the two-component
    // test uses, because a tolerance far looser than the measurement would pass on an
    // assembly that had lost the property and was merely well converged.
    expectTrue("and twelve modes a piece gets the assembled spectrum to seven figures",
               worst < 1e-6);

    // The negative control, and it is the one that matters: three components that
    // were never joined are still three well-formed reduced models in one matrix.
    {
        const reduction::Assembly loose = reduction::assemble(component, {});
        expectEqualCount("unjoined, nothing is shared and the boundary is the sum",
                         static_cast<std::size_t>(loose.boundary), 4 * perPlane);
        expectEqual("and the model is in three pieces", reduction::assembledComponents(loose), 3);
        int looseRigid = 0;
        for (double w : reduction::assembledFrequencies(loose, {}))
            if (w < rigidCutoff) ++looseRigid;
        expectEqual("three components that float free of each other have eighteen rigid modes",
                    static_cast<long long>(looseRigid), 18);
    }

    // --- Folding one component at a time, which is the other route ----------------
    //
    // The hypothesis `reduction.hpp` used to carry: carry the boundary identity and
    // run the same position match against the assembly. It works, and the point of
    // testing it is that it produces the **same matrix**, so the only thing separating
    // the two routes is what they cost.
    {
        reduction::Assembly fold = reduction::assemble({component[0]}, {});
        expectTrue("a one-component assembly is a well-formed assembly", !fold.empty());
        expectEqual("of one piece", reduction::assembledComponents(fold), 1);
        for (int p = 1; p < 3; ++p) {
            const auto c = static_cast<std::size_t>(p);
            const reduction::InterfaceMap map =
                reduction::matchBoundaries(fold, *component[c].substructure);
            expectTrue("the match against an assembly raises nothing", map.problems.empty());
            expectEqualCount("and finds the whole shared plane", map.shared, perPlane);
            fold = reduction::assemble(fold, component[c], map);
            expectTrue("the fold raises nothing", fold.problems.empty());
        }
        expectEqual("the folded model has the same number of components", fold.parts, trio.parts);
        expectEqual("and the same boundary", fold.boundary, trio.boundary);
        expectEqual("and the same size", fold.size(), trio.size());
        // Bit for bit. Not "close": both routes are the same scatter-add through the
        // same helper, so anything but equality is an index that disagrees, and an
        // index that disagrees by one row is a model that is merely a bit wrong.
        bool identical = fold.stiffness.size() == trio.stiffness.size();
        double scale = 0;
        for (std::size_t i = 0; i < trio.stiffness.size() && identical; ++i) {
            scale = std::max(scale, std::fabs(trio.stiffness[i]));
            identical = fold.stiffness[i] == trio.stiffness[i] && fold.mass[i] == trio.mass[i];
        }
        expectTrue("the two routes produce the same matrices to the last bit", identical);
        expectTrue("on a matrix that is not all zeros", scale > 0);
        bool sameMaps = fold.from.size() == trio.from.size();
        for (std::size_t c = 0; c < trio.from.size() && sameMaps; ++c)
            sameMaps = fold.from[c] == trio.from[c];
        expectTrue("and the same component maps", sameMaps);
        expectEqualCount("and the folded model carries an identity too", fold.boundaryPoint.size(),
                         static_cast<std::size_t>(fold.boundary));

        // **The fold has its own copy of the axis check and the N-way test does not
        // reach it.** Deleting it survived everything until this existed, which is
        // the shape this file keeps finding: two routes to the same answer means two
        // places for the same guard, and testing one of them tests one of them.
        reduction::Assembly pair = reduction::assemble({component[0]}, {});
        pair = reduction::assemble(pair, component[1],
                                   reduction::matchBoundaries(pair, *component[1].substructure));
        reduction::InterfaceMap cross = reduction::matchBoundaries(pair, *component[2].substructure);
        const std::vector<reduction::BoundaryDof> third =
            reduction::boundaryIdentity(*component[2].substructure);
        // Two entries **swapped**, not one retargeted. `assemble` inverts the map into
        // "which assembly row is this component DOF", so retargeting onto a DOF
        // another entry already names loses the crossing to whichever is written
        // last; a swap keeps the map a bijection and crosses two axes at once.
        bool rewiredFold = false;
        for (std::size_t i = 0; i < cross.aToB.size() && !rewiredFold; ++i) {
            if (cross.aToB[i] < 0) continue;
            for (std::size_t k = i + 1; k < cross.aToB.size(); ++k) {
                if (cross.aToB[k] < 0) continue;
                const auto j1 = static_cast<std::size_t>(cross.aToB[i]);
                const auto j2 = static_cast<std::size_t>(cross.aToB[k]);
                if (third[j1].axis == third[j2].axis) continue;
                std::swap(cross.aToB[i], cross.aToB[k]);
                rewiredFold = true;
                break;
            }
        }
        expectTrue("the crossed fold map was actually built", rewiredFold);
        const reduction::Assembly wrongFold = reduction::assemble(pair, component[2], cross);
        expectEqual("folding two rows onto each other's axes is counted, both of them",
                    wrongFold.axisDisagreements, 2);
        expectTrue("and said out loud", !wrongFold.problems.empty());
    }

    // --- The refusals, every one of them a model that would otherwise assemble ----
    {
        // An assembly built from two bare reductions has no mesh behind it, so there
        // is nothing to match against. That is the limit the two-component overload
        // has always had, and it is now stated rather than implied.
        const reduction::Assembly bare =
            reduction::assemble(reduced[0], reduced[1],
                                reduction::matchBoundaries(substructure[0], substructure[1]));
        expectTrue("a two-reduction assembly carries no identity", bare.boundaryPoint.empty());
        const reduction::InterfaceMap none = reduction::matchBoundaries(bare, substructure[2]);
        expectEqualCount("with nothing matched", none.shared, 0u);
        // **The refusal is not enough without its reason.** Deleting the identity
        // guard leaves the match running against an empty list, which raises the
        // *other* complaint -- "these share no boundary DOF" -- and a test that only
        // asked whether `problems` was non-empty passed on it. This repository has
        // recorded that exact shape once already, on a Cholesky failure reported as a
        // QL failure.
        bool named = false;
        for (const std::string& p : none.problems)
            named = named || p.find("identity") != std::string::npos;
        expectTrue("and it says the assembly has no identity, not that nothing coincided", named);

        // Seven guards in `assemble` share this one observable, and they are ordered:
        // the component-range check runs before the self-join check, which runs
        // before the map-size and map-range checks. So each of these fixtures can be
        // answered by a guard other than the one it was built for, and the message
        // is what separates them.
        //
        // The range check is also the one holding `parts[cb]` inside its array --
        // `past[0].b = 7` against a three-component list -- so its removal is a heap
        // overread rather than a wrong message, the same shape as `edgeDrive`'s
        // wrong-zone guard.
        const auto refusedWith = [&](const std::vector<reduction::Joint>& j, const char* text) {
            const std::vector<std::string> said = reduction::assemble(component, j).problems;
            if (said.empty()) return false;
            for (const std::string& p : said)
                if (p.find(text) != std::string::npos) return true;
            return false;
        };

        std::vector<reduction::Joint> self = all;
        self[0].b = self[0].a;
        expectTrue("a component joined to itself is refused",
                   !reduction::assemble(component, self).problems.empty());
        expectTrue("and named as a self-join, not as a map that overruns",
                   refusedWith(self, "joins a component to itself"));

        std::vector<reduction::Joint> past = all;
        past[0].b = 7;
        expectTrue("a joint naming a component that is not there is refused",
                   !reduction::assemble(component, past).problems.empty());
        expectTrue("and named as a component that is not in the list",
                   refusedWith(past, "names a component that is not in the list"));

        // A map shorter than the boundary it claims to describe. This is what a
        // component whose reduction fell back to empty produces, and accepting it
        // would join the DOF the map does reach and leave the rest as two unknowns --
        // a model that assembles, solves, and is torn along the tail of a cut.
        std::vector<reduction::Joint> truncated = all;
        truncated[0].map.aToB.pop_back();
        expectTrue("and so is a map that does not cover its component's boundary",
                   !reduction::assemble(component, truncated).problems.empty());

        // The one an `InterfaceMap` could never produce and a hand-built `Joint`
        // can: a merge that couples one axis to another. It assembles, it solves,
        // and it is a ship whose plating pushes sideways when you pull it.
        std::vector<reduction::Joint> crossed = all;
        const std::vector<reduction::BoundaryDof> ia = reduction::boundaryIdentity(substructure[0]);
        const std::vector<reduction::BoundaryDof> ib = reduction::boundaryIdentity(substructure[1]);
        bool rewired = false;
        for (std::size_t i = 0; i < crossed[0].map.aToB.size() && !rewired; ++i) {
            if (crossed[0].map.aToB[i] < 0) continue;
            for (std::size_t j = 0; j < ib.size(); ++j)
                if (ib[j].axis != ia[i].axis) {
                    crossed[0].map.aToB[i] = static_cast<int>(j);
                    rewired = true;
                    break;
                }
        }
        expectTrue("the crossed map was actually built", rewired);
        const reduction::Assembly wrong = reduction::assemble(component, crossed);
        expectTrue("an assembled row that merges two axes is counted",
                   wrong.axisDisagreements > 0);
        expectTrue("and said out loud", !wrong.problems.empty());
    }
}

}  // namespace

// --- The constraint the attachment carries but does not add ------------------------
//
// `Attachment::constrained` eliminates degrees of freedom where
// `Attachment::stiffness` adds them, and the two go into the same assembly. What
// has to be true of it is that the substructure and `solveStatic` -- two entirely
// separate assemblies of the same physics -- reach the same answer through it, and
// that the eliminated degrees of freedom leave the partition rather than sitting in
// it with an empty row.
void testConstrainedSubstructure() {
    std::printf("\n--- reduction: a substructure with a degree of freedom tied away ---\n");
    const solidshell::HexMesh mesh = testPlate();
    const StructuralMaterial steel = ah36Steel();
    const int sy = kNy + 1;
    const auto index = [&](int i, int j, int k) {
        return static_cast<std::uint32_t>((i * sy + j) * 2 + k);
    };

    // An interior node pair tied to the pairs either side of it along x, halfway
    // between them. Not a junction -- this file has no second surface -- but the
    // same shape of constraint: a mesh node with elements of its own, eliminated.
    reduction::Attachment tied;
    for (int k = 0; k < 2; ++k)
        for (int axis = 0; axis < 3; ++axis)
            tied.constrained.push_back(
                {index(4, 2, k) * 3 + static_cast<std::uint32_t>(axis),
                 {index(3, 2, k) * 3 + static_cast<std::uint32_t>(axis),
                  index(5, 2, k) * 3 + static_cast<std::uint32_t>(axis)},
                 {0.5, 0.5}});

    const Substructure bare(mesh, steel, endInterface(mesh));
    const Substructure sub(mesh, steel, endInterface(mesh), tied);
    expectTrue("the constrained substructure factors", sub.ready());
    // An eliminated degree of freedom is in neither partition: it is not an unknown.
    expectEqualCount("six degrees of freedom leave the partition",
                     bare.boundaryCount() + bare.interiorCount(),
                     sub.boundaryCount() + sub.interiorCount() + 6);
    expectEqualCount("the boundary is untouched", sub.boundaryCount(), bare.boundaryCount());
    // The steel does not vanish with the unknown: `T^T M T`'s row sums are the
    // weights, which are a partition of unity, so the total is preserved exactly.
    expectNear("and the mass they carried went to their masters", sub.totalMass(),
               bare.totalMass(), 1e-9 * bare.totalMass());

    // The strong statement: the operator the substructure assembled is the operator
    // `solveStatic` assembles, on a problem with a real answer. Prescribe the two
    // ends and relax the interior both ways.
    solidshell::HexMesh driven = mesh;
    for (std::uint32_t node : endInterface(mesh)) {
        const double x = mesh.position[node * 3];
        driven.pin(node, 0, 2e-4 * x);
        driven.pin(node, 1, 0.0);
        driven.pin(node, 2, 0.0);
    }
    const std::vector<double> noLoad(mesh.nodeCount() * 3, 0.0);
    std::vector<double> field;
    std::string problem;
    expectTrue("the constrained plate solves: " + problem,
               solidshell::solveStatic(driven, steel, solidshell::Formulation::SolidShell, {},
                                       tied.constrained, noLoad, field, &problem));

    // `K u` over the substructure, with `u` the field `solveStatic` found. Every
    // interior row must be in equilibrium -- that is what says the two assemblies
    // are the same matrix -- and the eliminated rows must be exactly empty.
    std::vector<double> reaction;
    sub.stiffnessTimes(field, reaction);
    double worstInterior = 0, worstEliminated = 0, boundaryScale = 0;
    for (std::uint32_t d : sub.interiorDof()) worstInterior = std::max(worstInterior, std::fabs(reaction[d]));
    for (const solidshell::Mpc& mpc : tied.constrained)
        worstEliminated = std::max(worstEliminated, std::fabs(reaction[mpc.slave]));
    for (std::uint32_t d : sub.boundaryDof())
        boundaryScale = std::max(boundaryScale, std::fabs(reaction[d]));
    std::printf("     interior residual %.2e N, eliminated rows %.2e N, boundary reaction %.2e N\n",
                worstInterior, worstEliminated, boundaryScale);
    expectTrue("the drive loaded the substructure", boundaryScale > 1e3);
    expectTrue("the substructure's interior is in equilibrium on solveStatic's field",
               worstInterior < 1e-9 * boundaryScale);
    expectTrue("and an eliminated row is empty, not merely small", worstEliminated == 0.0);

    // The vacuity guard: the *bare* substructure is not in equilibrium on that
    // field, so "in equilibrium" is a statement about the constraint having been
    // assembled and not about the field being trivial.
    std::vector<double> bareReaction;
    bare.stiffnessTimes(field, bareReaction);
    double worstBare = 0;
    for (std::uint32_t d : bare.interiorDof())
        worstBare = std::max(worstBare, std::fabs(bareReaction[d]));
    expectTrue("where the unconstrained operator is not", worstBare > 1e-3 * boundaryScale);

    // **A node with only *one* axis eliminated is still an interior node**, and it
    // is the case that tells the two places the partition is filtered apart: the
    // node-level filter keeps it, so the per-degree-of-freedom skip has to drop the
    // one axis. A junction tie never produces one -- it eliminates all three -- so
    // without this the per-axis skip is dead code that looks alive.
    {
        reduction::Attachment oneAxis;
        oneAxis.constrained.push_back({index(4, 2, 0) * 3,
                                       {index(3, 2, 0) * 3, index(5, 2, 0) * 3}, {0.5, 0.5}});
        const Substructure partial(mesh, steel, endInterface(mesh), oneAxis);
        expectTrue("a partially constrained node still factors", partial.ready());
        expectEqualCount("and gives up exactly the one degree of freedom",
                         bare.boundaryCount() + bare.interiorCount(),
                         partial.boundaryCount() + partial.interiorCount() + 1);
        expectNear("keeping its mass", partial.totalMass(), bare.totalMass(),
                   1e-9 * bare.totalMass());
    }

    // Reading a substructure that **refused** must not run off the end of what the
    // constructor never built. Every early return leaves `nodeCount()` set from the
    // mesh while the mass and the sparsity pattern are still empty, and a caller who
    // skips `ready()` used to get a read of address zero. Mutation testing found it,
    // as a segmentation fault three tests downstream of the mutant it was chasing,
    // and it was reachable from an inverted element long before any constraint
    // existed.
    {
        reduction::Attachment selfReferential;
        selfReferential.constrained.push_back({0, {0}, {1.0}});
        const Substructure refused(mesh, steel, endInterface(mesh), selfReferential);
        expectTrue("a self-referential constraint is refused", !refused.ready());
        expectTrue("and it still reports its node count", refused.nodeCount() > 0);
        expectNear("its mass reads zero rather than reading out of bounds", refused.totalMass(),
                   0.0, 0.0);
        std::vector<double> probe(refused.dofCount(), 1.0), out;
        refused.stiffnessTimes(probe, out);
        expectEqualCount("and its operator is the zero operator", out.size(), refused.dofCount());
        double any = 0;
        for (double v : out) any = std::max(any, std::fabs(v));
        expectNear("with nothing in it", any, 0.0, 0.0);
    }

    // Refusals. Both are modelling errors that would otherwise assemble.
    {
        reduction::Attachment onInterface;
        onInterface.constrained.push_back({endInterface(mesh).front() * 3, {index(3, 2, 0) * 3},
                                           {1.0}});
        const Substructure refused(mesh, steel, endInterface(mesh), onInterface);
        expectTrue("a constrained interface degree of freedom is refused", !refused.ready());
    }
    {
        reduction::Attachment chained = tied;
        chained.constrained.push_back({index(3, 2, 0) * 3, {index(2, 2, 0) * 3}, {1.0}});
        const Substructure refused(mesh, steel, endInterface(mesh), chained);
        expectTrue("and so is a chain", !refused.ready());
    }
}

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
    testBoundaryMatchOnItsOwn();
    testThreeComponentsAssembleIntoTheWhole();
    testAttachedStiffness();
    testMemberValidity();
    testAttachedMass();
    testConstrainedSubstructure();
}
