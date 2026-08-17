// SPDX-License-Identifier: MIT
//
// Validation of the solid-shell element.
//
// Elasticity hands out exact answers and this file spends them. Nothing here is
// asserted against a number that came out of this code:
//
//   * the **patch test** -- a patch of distorted elements under a linear
//     displacement field must reproduce constant stress *exactly*. An element that
//     fails it is wrong at every mesh size, and no amount of refinement rescues
//     it. Run on a regular grid, on an in-plane distorted grid (which is what a
//     real plate mesh is), and on a warped one, because a patch of rectangles
//     passes tests a real mesh does not;
//   * **rigid body motion**, including a finite rotation, must produce exactly
//     zero internal force -- the one thing the co-rotational formulation exists
//     to get right -- and the stronger form of it, **frame indifference**:
//     rotating an element that is already 25% stretched must rotate its internal
//     force and change nothing else. The weaker version is not enough on its own,
//     because rotating an *undeformed* element hands the polar decomposition a
//     deformation gradient that is already a rotation, and it comes back
//     untouched however badly the decomposition is implemented;
//   * **beam and plate closed forms**: PL^3/3EI, PL^3/3Db for a strip in
//     cylindrical bending, 5qL^4/384D and qL^4/384D for supported and clamped
//     strips, and the Navier double series for a simply supported plate, summed
//     here rather than quoted;
//   * **the rank of the element stiffness**: exactly six zero eigenvalues and
//     eighteen positive ones. That is what says the element has no spurious
//     zero-energy mode, which is the failure a formulation with assumed strains
//     can acquire and which no deflection test would notice.
//
// Two of those closed forms turned out to be the *approximate* party, and both are
// checked as such rather than being given a slack tolerance. Beam theory is 1.4%
// stiffer than three-dimensional elasticity for a strip eight thicknesses wide,
// because the clamped root suppresses anticlastic curvature; the residual is
// therefore swept against b/t, where it must vanish. Kirchhoff plate theory
// neglects transverse shear; that residual is swept against a/t, where it must
// fall as (t/a)^2. In both cases the sweep is the evidence that what is left over
// belongs to the theory and not to the element.
//
// Roughly a third of this file exists because **mutation testing** said the rest of
// it was not enough. Fifty-two single-edit mutants of solid_shell.cpp; the first
// pass killed forty-one and every one of the eleven survivors was a real hole --
// four of the seven enhanced parameters were never loaded in the plane where they
// act, the enhanced strains could be discarded on recovery without any deflection
// test noticing, and the mass lumping and pressure load were only ever asked about
// on rectangles, where every plausible wrong answer is also the right one.
// docs/07-fem-spike-findings.md section 6 records the full account.
#include "engine/sim/fem.hpp"
#include "engine/sim/solid_shell.hpp"
#include "harness.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <limits>
#include <map>
#include <numbers>
#include <random>
#include <string>
#include <utility>
#include <vector>

using namespace sim;
using namespace sim::solidshell;
using testing::expectEqual;
using testing::expectNear;
using testing::expectTrue;

namespace {

StructuralMaterial steelMaterial() { return ah36Steel(); }

// --- small helpers ------------------------------------------------------------

// Jacobi eigenvalues of a symmetric n x n matrix, ascending. Used to establish the
// rank of the element stiffness, which is the test for spurious modes.
std::vector<double> symmetricEigenvalues(const double* matrix, int n) {
    std::vector<double> a(matrix, matrix + n * n);
    for (int sweep = 0; sweep < 200; ++sweep) {
        double off = 0;
        for (int p = 0; p < n; ++p)
            for (int q = p + 1; q < n; ++q) off += a[p * n + q] * a[p * n + q];
        if (off < 1e-30) break;
        for (int p = 0; p < n; ++p)
            for (int q = p + 1; q < n; ++q) {
                const double apq = a[p * n + q];
                if (std::abs(apq) < 1e-300) continue;
                const double theta = (a[q * n + q] - a[p * n + p]) / (2.0 * apq);
                const double t = (theta >= 0 ? 1.0 : -1.0) /
                                 (std::abs(theta) + std::sqrt(theta * theta + 1.0));
                const double c = 1.0 / std::sqrt(t * t + 1.0), s = t * c;
                for (int k = 0; k < n; ++k) {
                    const double akp = a[k * n + p], akq = a[k * n + q];
                    a[k * n + p] = c * akp - s * akq;
                    a[k * n + q] = s * akp + c * akq;
                }
                for (int k = 0; k < n; ++k) {
                    const double apk = a[p * n + k], aqk = a[q * n + k];
                    a[p * n + k] = c * apk - s * aqk;
                    a[q * n + k] = s * apk + c * aqk;
                }
            }
    }
    std::vector<double> eigenvalues(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) eigenvalues[static_cast<std::size_t>(i)] = a[i * n + i];
    std::sort(eigenvalues.begin(), eigenvalues.end());
    return eigenvalues;
}

// Analytic stress of a constant strain state, for the patch test to compare with.
void constantStress(const StructuralMaterial& material, const double gradient[3][3],
                    double out[6]) {
    const double nu = material.poissonRatio, e = material.youngsModulus;
    const double lambda = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    const double mu = e / (2.0 * (1.0 + nu));
    const double strain[6] = {gradient[0][0],
                              gradient[1][1],
                              gradient[2][2],
                              gradient[0][1] + gradient[1][0],
                              gradient[1][2] + gradient[2][1],
                              gradient[2][0] + gradient[0][2]};
    const double trace = strain[0] + strain[1] + strain[2];
    for (int i = 0; i < 3; ++i) {
        out[i] = lambda * trace + 2.0 * mu * strain[i];
        out[3 + i] = mu * strain[3 + i];
    }
}

void clampFace(HexMesh& mesh, int axis, double value, bool onlyNormal = false) {
    for (std::size_t n = 0; n < mesh.nodeCount(); ++n)
        if (std::abs(mesh.position[n * 3 + static_cast<std::size_t>(axis)] - value) < 1e-12) {
            if (onlyNormal) {
                mesh.pin(n, axis, 0.0);
            } else {
                for (int a = 0; a < 3; ++a) mesh.pin(n, a, 0.0);
            }
        }
}

// Average vertical displacement over the nodes at x = value.
double sectionDeflection(const HexMesh& mesh, const std::vector<double>& displacement,
                         double value) {
    double sum = 0;
    int count = 0;
    for (std::size_t n = 0; n < mesh.nodeCount(); ++n)
        if (std::abs(mesh.position[n * 3] - value) < 1e-12) {
            sum += displacement[n * 3 + 2];
            ++count;
        }
    return count > 0 ? sum / count : 0.0;
}

std::vector<double> endLoad(const HexMesh& mesh, double x, double force) {
    std::vector<double> load(mesh.nodeCount() * 3, 0.0);
    int count = 0;
    for (std::size_t n = 0; n < mesh.nodeCount(); ++n)
        if (std::abs(mesh.position[n * 3] - x) < 1e-12) ++count;
    for (std::size_t n = 0; n < mesh.nodeCount(); ++n)
        if (std::abs(mesh.position[n * 3] - x) < 1e-12)
            load[n * 3 + 2] = -force / std::max(count, 1);
    return load;
}

// Cylindrical bending: no displacement across the width anywhere, so the strip is
// in exact plane strain in y and the flexural rigidity is E t^3 / 12(1 - nu^2)
// with no anticlastic ambiguity left for a tolerance to absorb.
void restrainWidth(HexMesh& mesh) {
    for (std::size_t n = 0; n < mesh.nodeCount(); ++n) mesh.pin(n, 1, 0.0);
}

double flexuralRigidity(const StructuralMaterial& material, double thickness) {
    const double nu = material.poissonRatio;
    return material.youngsModulus * thickness * thickness * thickness /
           (12.0 * (1.0 - nu * nu));
}

// --- 1. The element on its own: rank, mass, geometry, loads --------------------

// Symmetry, the six rigid body modes, and -- the assertion that a deflection test
// cannot make -- that there are exactly six of them. An assumed-strain element
// with too many enhanced modes goes rank deficient and develops a zero-energy
// mode that a cantilever would happily bend around without ever showing it.
void testElementStiffnessRank() {
    std::printf("\n--- solid-shell elements ---\n");
    const StructuralMaterial steel = steelMaterial();

    struct Case {
        const char* label;
        double nodes[kDof];
    };
    // A cube; a plate element at the aspect ratio the formulation is for; and one
    // distorted in every direction at once, because a rank deficiency can hide on
    // a regular shape.
    const Case cases[] = {
        {"cube",
         {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1}},
        {"plate element 50x50x10 mm",
         {0, 0, -0.005, 0.05, 0, -0.005, 0.05, 0.05, -0.005, 0, 0.05, -0.005, 0, 0, 0.005, 0.05, 0,
          0.005, 0.05, 0.05, 0.005, 0, 0.05, 0.005}},
        {"distorted in all three directions",
         {0.03, -0.02, 0.01, 1.07, 0.04, -0.03, 0.94, 1.11, 0.05, -0.05, 0.96, -0.02, 0.02, 0.06,
          1.05, 1.09, -0.03, 0.92, 1.03, 0.97, 1.12, 0.05, 1.04, 0.95}},
    };

    for (const Case& c : cases) {
        expectTrue(std::string("the ") + c.label + " element is not inverted",
                   smallestJacobian(c.nodes) > 0.0);
        double k[kDof * kDof];
        elementStiffness(c.nodes, steel, Formulation::SolidShell, k);

        double scale = 0, asymmetry = 0;
        for (int i = 0; i < kDof; ++i)
            for (int j = 0; j < kDof; ++j) {
                scale = std::max(scale, std::abs(k[i * kDof + j]));
                asymmetry = std::max(asymmetry, std::abs(k[i * kDof + j] - k[j * kDof + i]));
            }
        expectTrue(std::string("stiffness is symmetric (") + c.label + ")",
                   asymmetry < 1e-12 * scale);

        const std::vector<double> eigenvalues = symmetricEigenvalues(k, kDof);
        // Six rigid body modes and no more. The gap between the sixth and the
        // seventh is the whole assertion, so it is measured rather than assumed:
        // a spurious mode would land in between.
        const double smallest = eigenvalues[6];
        double largestZero = 0;
        for (int i = 0; i < 6; ++i) largestZero = std::max(largestZero, std::abs(eigenvalues[i]));
        expectTrue(std::string("exactly six zero-energy modes (") + c.label + ")",
                   largestZero < 1e-9 * scale && smallest > 1e-6 * scale);
        expectTrue(std::string("no negative eigenvalue (") + c.label + ")",
                   eigenvalues[0] > -1e-9 * scale);
        std::printf("     %-34s rigid modes to %8.1e, first real mode %8.1e (of %8.1e)\n", c.label,
                    largestZero / scale, smallest / scale, 1.0);

        // Rigid body modes explicitly, not just by counting: three translations and
        // three infinitesimal rotations must all be annihilated.
        for (int mode = 0; mode < 6; ++mode) {
            double u[kDof] = {};
            for (int a = 0; a < kNodes; ++a) {
                if (mode < 3) {
                    u[a * 3 + mode] = 1.0;
                } else {
                    const int axis = mode - 3;
                    const int i = (axis + 1) % 3, j = (axis + 2) % 3;
                    u[a * 3 + i] = -c.nodes[a * 3 + j];
                    u[a * 3 + j] = c.nodes[a * 3 + i];
                }
            }
            double worst = 0;
            for (int i = 0; i < kDof; ++i) {
                double s = 0;
                for (int j = 0; j < kDof; ++j) s += k[i * kDof + j] * u[j];
                worst = std::max(worst, std::abs(s));
            }
            expectTrue(std::string("rigid mode ") + std::to_string(mode) + " carries no force (" +
                           c.label + ")",
                       worst < 1e-10 * scale);
        }
    }
}

// A trilinear hexahedron's volume and first moment, by a four-point Gauss rule
// written out here rather than reusing the element's own quadrature. It is exact
// for both -- det J is quadratic in each coordinate and x det J cubic, and a
// four-point rule is exact to seventh order -- so it is a genuine second opinion
// on the element's lumped mass, not the same integral run twice.
void hexVolumeAndMoment(const double nodes[kDof], double& volume, double moment[3]) {
    static constexpr double kAbscissa[4] = {-0.8611363115940526, -0.3399810435848563,
                                            0.3399810435848563, 0.8611363115940526};
    static constexpr double kWeight[4] = {0.3478548451374538, 0.6521451548625461,
                                          0.6521451548625461, 0.3478548451374538};
    static constexpr double kXi[kNodes] = {-1, +1, +1, -1, -1, +1, +1, -1};
    static constexpr double kEta[kNodes] = {-1, -1, +1, +1, -1, -1, +1, +1};
    static constexpr double kZta[kNodes] = {-1, -1, -1, -1, +1, +1, +1, +1};
    volume = 0;
    for (int i = 0; i < 3; ++i) moment[i] = 0;
    for (int a = 0; a < 4; ++a)
        for (int b = 0; b < 4; ++b)
            for (int c = 0; c < 4; ++c) {
                const double xi = kAbscissa[a], eta = kAbscissa[b], zeta = kAbscissa[c];
                double shape[kNodes], derivative[kNodes][3];
                for (int n = 0; n < kNodes; ++n) {
                    const double x = 1.0 + xi * kXi[n], y = 1.0 + eta * kEta[n],
                                 z = 1.0 + zeta * kZta[n];
                    shape[n] = 0.125 * x * y * z;
                    derivative[n][0] = 0.125 * kXi[n] * y * z;
                    derivative[n][1] = 0.125 * kEta[n] * x * z;
                    derivative[n][2] = 0.125 * kZta[n] * x * y;
                }
                double jacobian[3][3] = {};
                double point[3] = {0, 0, 0};
                for (int r = 0; r < 3; ++r) {
                    for (int k = 0; k < 3; ++k)
                        for (int n = 0; n < kNodes; ++n)
                            jacobian[r][k] += derivative[n][k] * nodes[n * 3 + r];
                    for (int n = 0; n < kNodes; ++n) point[r] += shape[n] * nodes[n * 3 + r];
                }
                const double determinant =
                    jacobian[0][0] *
                        (jacobian[1][1] * jacobian[2][2] - jacobian[1][2] * jacobian[2][1]) -
                    jacobian[0][1] *
                        (jacobian[1][0] * jacobian[2][2] - jacobian[1][2] * jacobian[2][0]) +
                    jacobian[0][2] *
                        (jacobian[1][0] * jacobian[2][1] - jacobian[1][1] * jacobian[2][0]);
                const double w = kWeight[a] * kWeight[b] * kWeight[c] * determinant;
                volume += w;
                for (int r = 0; r < 3; ++r) moment[r] += w * point[r];
            }
}

// Lumped mass, on a *distorted* element. On a box every reasonable lumping gives
// volume/8 and any of them passes, which is why this is asserted on a shape where
// they differ: the row sum reproduces the element's centre of mass exactly and an
// even eighth does not.
void testLumpedMass() {
    const StructuralMaterial steel = steelMaterial();
    const double nodes[kDof] = {0.03, -0.02, 0.01,  1.07,  0.04, -0.03, 0.94, 1.11,
                                0.05, -0.05, 0.96,  -0.02, 0.02, 0.06,  1.05, 1.09,
                                -0.03, 0.92, 1.03,  0.97,  1.12, 0.05,  1.04, 0.95};
    double mass[kNodes];
    elementMass(nodes, steel.density, mass);

    double volume = 0, moment[3];
    hexVolumeAndMoment(nodes, volume, moment);
    double total = 0, centre[3] = {0, 0, 0};
    for (int a = 0; a < kNodes; ++a) {
        expectTrue("every lumped nodal mass is positive", mass[a] > 0.0);
        total += mass[a];
        for (int i = 0; i < 3; ++i) centre[i] += mass[a] * nodes[a * 3 + i];
    }
    expectNear("lumped mass sums to density times volume", total, steel.density * volume,
               1e-10 * steel.density * volume);
    for (int i = 0; i < 3; ++i)
        expectNear("and sits at the element's own centre of mass, axis " + std::to_string(i),
                   centre[i] / total, moment[i] / volume, 1e-10);

    // Non-vacuous: on this element an even eighth would land somewhere else.
    double evenCentre[3] = {0, 0, 0};
    for (int a = 0; a < kNodes; ++a)
        for (int i = 0; i < 3; ++i) evenCentre[i] += nodes[a * 3 + i] / kNodes;
    double gap = 0;
    for (int i = 0; i < 3; ++i)
        gap = std::max(gap, std::abs(evenCentre[i] - moment[i] / volume));
    expectTrue("the element is distorted enough for the lumping to matter", gap > 1e-3);
}

// The inverted-element guard. It has to look at every corner: a hexahedron can be
// perfectly well behaved where it is first sampled and folded through itself at
// the far one, and everything computed from it after that is meaningless.
void testInvertedElementIsCaught() {
    const StructuralMaterial steel = steelMaterial();
    double folded[kDof] = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0,
                           0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1};
    // Node 6 dragged back across the diagonal. The Jacobian at corner 0 depends
    // only on nodes 0, 1, 3 and 4, so it stays positive and a check that sampled
    // one corner would see nothing wrong.
    folded[6 * 3 + 0] = -0.4;
    folded[6 * 3 + 1] = -0.4;
    expectTrue("the folded element is caught", smallestJacobian(folded) <= 0.0);

    double good[kDof] = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0,
                         0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1};
    expectTrue("a sound element is not", smallestJacobian(good) > 0.0);

    // And it reaches the caller rather than producing a plausible wrong answer.
    // Note the indirection: the mesh's node numbering is not the element's, and
    // reaching into `position` by the local index would move a different corner.
    HexMesh mesh = makePlateMesh(1.0, 1.0, 0.2, 1, 1, 1);
    const std::size_t corner = mesh.index[6];
    mesh.position[corner * 3 + 0] = -0.4;
    mesh.position[corner * 3 + 1] = -0.4;
    clampFace(mesh, 2, -0.1);
    std::vector<double> load(mesh.nodeCount() * 3, 0.0), displacement;
    std::string problem;
    const bool solved =
        solveStatic(mesh, steel, Formulation::SolidShell, load, displacement, &problem);
    expectTrue("solveStatic refuses an inverted element", !solved && !problem.empty());
}

// --- The collapsed element, which is a wedge and not a fold ------------------------
//
// `smallestJacobian` samples the corners and so cannot tell a **collapsed**
// hexahedron -- the triangular prism a degenerate plate panel extrudes to -- from an
// inverted one. `elementShape` can, and everything below is a closed form on the
// simplest such element there is: the right triangle (0,0), (1,0), (0,1) written as a
// quad with its fourth corner on its first, extruded +/- t/2 in z.
//
// For that element the bilinear map collapses the whole `xi = -1` edge onto the
// origin, and
//
//     dX/dxi = ((1-eta)/4, (1+eta)/4, 0)
//     dX/deta = (-(1+xi)/4, (1+xi)/4, 0)
//     dX/dzta = (0, 0, t/2)
//     det J   = t (1 + xi) / 16
//
// -- exactly zero on the collapsed edge, positive everywhere else, and **linear in
// xi**, so the 2x2x2 rule integrates the wedge's volume exactly. Every number below
// is that expression rather than a reading.
void testCollapsedElementIsAWedgeAndNotAFold() {
    std::printf("\n   collapsed element: a triangular prism written in eight nodes\n");
    const double t = 0.012;
    // Corners 0-3 on zeta = -1 counter-clockwise seen from +zeta, 4-7 above them.
    const double wedge[kDof] = {0, 0, -t / 2, 1, 0, -t / 2, 0, 1, -t / 2, 0, 0, -t / 2,
                                0, 0, t / 2,  1, 0, t / 2,  0, 1, t / 2,  0, 0, t / 2};
    const ElementShape shape = elementShape(wedge);

    // Node 0 sits on node 3 and node 4 on node 7: four corners of eight.
    expectEqual("the wedge has four collapsed corners", shape.collapsedNodes, 4);
    expectTrue("and they are the ones that coincide",
               shape.collapsed[0] && shape.collapsed[3] && shape.collapsed[4] && shape.collapsed[7] &&
                   !shape.collapsed[1] && !shape.collapsed[2] && !shape.collapsed[5] &&
                   !shape.collapsed[6]);

    // The nodal determinant is *exactly* zero, not small: at xi = -1 the expression
    // above is 0. So `smallestJacobian` refuses this element, and that refusal is
    // the whole defect -- 64% of the reference ferry could not be meshed because of
    // it. Asserted at 0.0 tolerance because the arithmetic is exact.
    expectNear("the collapsed corners give a determinant of exactly zero",
               shape.smallestNodal, 0.0, 0.0);
    expectTrue("so the corner-sampling check refuses it", !(smallestJacobian(wedge) > 0.0));

    // And the quadrature does not see it. The smallest of the nine samples is the
    // Gauss point at xi = -1/sqrt(3).
    const double smallest = t * (1.0 - 1.0 / std::sqrt(3.0)) / 16.0;
    expectNear("the smallest quadrature determinant is t(1 - 1/sqrt 3)/16",
               shape.smallestGauss, smallest, 1e-15 * smallest);
    expectTrue("which is positive", shape.smallestGauss > 0.0);
    expectTrue("so the element is integrable", shape.integrable);

    // The volume, integrated by the rule the element actually uses. det J is linear
    // in xi and constant in eta and zeta, so 2x2x2 Gauss is exact and the answer is
    // the triangle's own area times the thickness -- to the last bit, which is why
    // this is asserted at 1e-15 rather than at a convergence tolerance.
    double volumes[kGauss];
    gaussVolumes(wedge, volumes);
    double volume = 0;
    for (double v : volumes) volume += v;
    expectNear("and the rule integrates the wedge's volume exactly", volume, 0.5 * t,
               1e-15 * 0.5 * t);

    // A sound hexahedron is unchanged by all of this: nothing collapsed, and the two
    // determinants agree because the element is a cuboid.
    const double cube[kDof] = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0,
                               0, 0, t, 1, 0, t, 1, 1, t, 0, 1, t};
    const ElementShape sound = elementShape(cube);
    expectEqual("a sound element has no collapsed corner", sound.collapsedNodes, 0);
    expectTrue("and is integrable", sound.integrable);
    expectNear("with the same determinant at its corners and its Gauss points",
               sound.smallestNodal, sound.smallestGauss, 1e-15 * sound.smallestNodal);

    // --- The three negative controls, because a check that accepts a wedge could
    // just as easily have been a check that accepts anything -------------------------

    // 1. The folded element from `testInvertedElementIsCaught`. Nothing coincides,
    //    so there is nothing to excuse the non-positive corner.
    double folded[kDof] = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0,
                           0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1};
    folded[6 * 3 + 0] = -0.4;
    folded[6 * 3 + 1] = -0.4;
    const ElementShape bad = elementShape(folded);
    expectEqual("a folded element has no collapsed corner", bad.collapsedNodes, 0);
    expectTrue("and is refused", !bad.integrable);

    // 2. **Collapsed *and* folded.** This is the one that matters: an element that
    //    has a genuine wedge edge and is also inside out must not be waved through
    //    on the strength of the wedge. Node 2 is dragged back past the collapsed
    //    edge, which turns the triangle over without un-collapsing it.
    double foldedWedge[kDof];
    for (int i = 0; i < kDof; ++i) foldedWedge[i] = wedge[i];
    foldedWedge[2 * 3 + 1] = -1.0;
    foldedWedge[6 * 3 + 1] = -1.0;
    const ElementShape both = elementShape(foldedWedge);
    expectEqual("a folded wedge is still collapsed", both.collapsedNodes, 4);
    expectTrue("and its quadrature is negative", both.smallestGauss < 0.0);
    expectTrue("so it is refused", !both.integrable);

    // 3. **The control that says `integrable` is not just `smallestGauss > 0`.** A
    //    wedge folded *through its thickness* at one sound corner: node 5 pushed
    //    0.6 t below node 1, so the thickness direction reverses there while the
    //    other three corners hold the quadrature positive. No in-plane move would
    //    do -- a collapsed quad's determinant is `(1 + xi) det(P1, P2) / 8`, whose
    //    sign does not vary over the element, so a triangle folds all at once or
    //    not at all and the thickness is the only direction that can fold locally.
    double corneredWedge[kDof];
    for (int i = 0; i < kDof; ++i) corneredWedge[i] = wedge[i];
    corneredWedge[5 * 3 + 2] = -t / 2 - 0.6 * t;
    const ElementShape cornered = elementShape(corneredWedge);
    expectEqual("the third control is still a wedge", cornered.collapsedNodes, 4);
    expectTrue("its quadrature alone would accept it", cornered.smallestGauss > 0.0);
    expectTrue("but a sound corner has folded", cornered.smallestNodal < 0.0);
    expectTrue("so it is refused", !cornered.integrable);

    // 4. **A sliver is not a wedge.** Coincidence is tested by exact equality, not by
    //    a tolerance, and this is what that buys: pull the collapsed pair a nanometre
    //    apart and the element stops being a wedge, so the same fold at the same
    //    corner is judged on its own determinant again. A tolerance of even 1e-12 m
    //    here would excuse a *nearly* collapsed element -- turning "this is a
    //    triangle" into "this is almost a triangle", which is the loosening the whole
    //    classification exists to avoid. Mutation testing found this: replacing the
    //    equality with a 1e-6 m distance survived everything else in the suite.
    double sliver[kDof];
    for (int i = 0; i < kDof; ++i) sliver[i] = corneredWedge[i];
    // **Into the element, not out of it.** Node 3 is the collapsed one, and pulling
    // it anywhere except towards the triangle's interior -- away from node 1, towards
    // node 2 -- turns the quad itself inside out, which would make the element
    // refused for the nudge rather than for the fold under test. Measured: of the
    // eight nanometre offsets in the plane, only this one leaves a positive corner.
    sliver[3 * 3 + 0] = -1e-9;
    sliver[3 * 3 + 1] = 1e-9;
    sliver[7 * 3 + 0] = -1e-9;
    sliver[7 * 3 + 1] = 1e-9;
    const ElementShape thin = elementShape(sliver);
    expectEqual("a pair a nanometre apart is not collapsed", thin.collapsedNodes, 0);
    expectTrue("its corner determinant is tiny but the sign is what is read",
               thin.smallestNodal < 0.0);
    expectTrue("and with nothing to excuse it, the element is refused", !thin.integrable);
    // The vacuity guard for *this* control: pulled apart the other way, with no fold,
    // a sliver is still a legal element and must not be refused -- otherwise the
    // assertion above would pass on a rule that rejected everything near-degenerate.
    double soundSliver[kDof];
    for (int i = 0; i < kDof; ++i) soundSliver[i] = wedge[i];
    soundSliver[3 * 3 + 0] = -1e-9;
    soundSliver[3 * 3 + 1] = 1e-9;
    soundSliver[7 * 3 + 0] = -1e-9;
    soundSliver[7 * 3 + 1] = 1e-9;
    const ElementShape ok = elementShape(soundSliver);
    expectEqual("a sound sliver is not collapsed either", ok.collapsedNodes, 0);
    expectTrue("but it is accepted, because nothing about it is inside out", ok.integrable);
}

// **Which nine points `elementShape` looks at, checked against a second evaluation of
// the same determinant.**
//
// Mutation testing put this here. Dropping the centre sample, dropping the last Gauss
// point, sampling one zeta level instead of two, and dropping the quadrature test from
// `integrable` altogether all survived a suite that only ever fed the predicate
// elements it was happy with. A rule is only tested by an input whose verdict it
// changes.
//
// `determinantAt` below is the trilinear map's Jacobian written out independently --
// three partial derivatives and a triple product -- so the comparison is two
// implementations of the same quantity rather than one restated.
double determinantAt(const double nodes[kDof], double xi, double eta, double zta) {
    static constexpr double kXiOf[kNodes] = {-1, +1, +1, -1, -1, +1, +1, -1};
    static constexpr double kEtaOf[kNodes] = {-1, -1, +1, +1, -1, -1, +1, +1};
    static constexpr double kZtaOf[kNodes] = {-1, -1, -1, -1, +1, +1, +1, +1};
    double d[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    for (int a = 0; a < kNodes; ++a) {
        const double x = 1.0 + xi * kXiOf[a], y = 1.0 + eta * kEtaOf[a], z = 1.0 + zta * kZtaOf[a];
        const double dn[3] = {0.125 * kXiOf[a] * y * z, 0.125 * kEtaOf[a] * x * z,
                              0.125 * kZtaOf[a] * x * y};
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) d[i][j] += dn[i] * nodes[a * 3 + j];
    }
    return d[0][0] * (d[1][1] * d[2][2] - d[1][2] * d[2][1]) -
           d[0][1] * (d[1][0] * d[2][2] - d[1][2] * d[2][0]) +
           d[0][2] * (d[1][0] * d[2][1] - d[1][1] * d[2][0]);
}

void testElementShapeSamplesTheQuadratureAndNotOnlyTheCorners() {
    std::printf("\n   element shape: the nine points it looks at\n");
    const double q = 1.0 / std::sqrt(3.0);

    // **The element that separates the two halves of `integrable`**: every one of its
    // eight corners is sound and its interior quadrature is not, so the corner rule
    // alone accepts it. Found by searching a quarter-unit lattice rather than derived
    // -- a trilinear hexahedron can be convex at all eight corners and still fold
    // through itself inside, and no closed form hands one over.
    const double concave[kDof] = {2.00,  1.00,  1.75, 0.00,  1.00, 1.50, 1.50,  2.75, 0.00,
                                  -0.75, 0.00,  -1.00, -0.75, -0.75, 1.50, 1.50,  -0.75, 2.00,
                                  0.25,  3.00,  2.75, -0.75, 2.75, 2.25};
    const ElementShape bent = elementShape(concave);
    std::printf("     corner-sound, quadrature-folded: nodal %+.6f, Gauss %+.6f\n",
                bent.smallestNodal, bent.smallestGauss);
    expectEqual("nothing about it is collapsed", bent.collapsedNodes, 0);
    expectTrue("every corner is sound, so the corner rule would accept it",
               bent.smallestNodal > 0.0);
    expectTrue("but the quadrature folds through itself", bent.smallestGauss < 0.0);
    expectTrue("so the element is refused", !bent.integrable);
    // And it reaches a caller: `computeForms` requires the same positivity, so an
    // element like this must not reach a solve either.
    {
        double stiffness[kDof * kDof];
        elementStiffness(concave, steelMaterial(), Formulation::SolidShell, stiffness);
        double sum = 0;
        for (int i = 0; i < kDof * kDof; ++i) sum += std::abs(stiffness[i]);
        expectNear("and the element refuses to produce a stiffness at all", sum, 0.0, 0.0);
    }

    // **The sample set itself, swept.** Two hundred deterministic pseudo-random
    // elements: for each, the two determinants `elementShape` reports must equal an
    // independent minimum over the eight corners and over the centre plus the 2x2x2
    // rule. A sweep rather than three hand-written cases, because which of the nine
    // points is the smallest depends on the element, and a rule that dropped one
    // would be caught only by an element whose minimum was there.
    std::uint64_t state = 0x9e3779b97f4a7c15ull;
    const auto next = [&] {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        return 2.0 * (static_cast<double>(state >> 11) / static_cast<double>(1ull << 53)) - 1.0;
    };
    const double base[kDof] = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0,
                               0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1};
    int centreWorst = 0, lowLevelWorst = 0, highLevelWorst = 0;
    double worstNodalGap = 0, worstGaussGap = 0;
    for (int trial = 0; trial < 200; ++trial) {
        double nodes[kDof];
        for (int i = 0; i < kDof; ++i) nodes[i] = base[i] + 0.8 * next();
        const ElementShape shape = elementShape(nodes);

        static constexpr double kXiOf[kNodes] = {-1, +1, +1, -1, -1, +1, +1, -1};
        static constexpr double kEtaOf[kNodes] = {-1, -1, +1, +1, -1, -1, +1, +1};
        static constexpr double kZtaOf[kNodes] = {-1, -1, -1, -1, +1, +1, +1, +1};
        double nodal = 1e300;
        for (int a = 0; a < kNodes; ++a)
            nodal = std::min(nodal, determinantAt(nodes, kXiOf[a], kEtaOf[a], kZtaOf[a]));
        double gauss = determinantAt(nodes, 0.0, 0.0, 0.0);
        int argmin = -1;  // -1 the centre, otherwise the Gauss index
        for (int gp = 0; gp < kGauss; ++gp) {
            const double value = determinantAt(nodes, (gp & 1) ? q : -q, (gp & 2) ? q : -q,
                                               (gp & 4) ? q : -q);
            if (value < gauss) {
                gauss = value;
                argmin = gp;
            }
        }
        const double scale = std::max(1.0, std::abs(gauss));
        worstNodalGap = std::max(worstNodalGap, std::abs(shape.smallestNodal - nodal) / scale);
        worstGaussGap = std::max(worstGaussGap, std::abs(shape.smallestGauss - gauss) / scale);
        if (argmin < 0) ++centreWorst;
        else if (argmin & 4) ++highLevelWorst;
        else ++lowLevelWorst;
    }
    std::printf("     over 200 elements: worst disagreement %.2e (corners) %.2e (quadrature);"
                " the smallest sample was the centre %d times, zeta=-q %d, zeta=+q %d\n",
                worstNodalGap, worstGaussGap, centreWorst, lowLevelWorst, highLevelWorst);
    expectTrue("the corner minimum is the eight corners", worstNodalGap < 1e-12);
    expectTrue("and the quadrature minimum is the centre and the 2x2x2 rule",
               worstGaussGap < 1e-12);
    // **The vacuity guards, and they are what make the sweep worth running.** Unless
    // the smallest sample lands on both zeta levels somewhere in the sweep, dropping
    // one of them would change nothing and the comparison above would be eight points
    // against eight points.
    expectTrue("the smallest sample is on the lower zeta level for some element",
               lowLevelWorst > 0);
    expectTrue("and on the upper for some other", highLevelWorst > 0);

    // The centre is not the smallest sample of a random element -- it did not happen
    // once in two hundred -- so it needs one chosen for it, or dropping it would cost
    // nothing here. This one's centre determinant is a quarter of its smallest Gauss
    // point. Found by the same lattice search.
    expectEqual("no random element had its minimum at the centre", centreWorst, 0);
    const double pinched[kDof] = {0.25, 1.00, 1.50,  0.00, -0.25, -0.25, 2.50, 2.75, 2.00,
                                  1.00, 1.75, 0.25,  1.50, 0.00,  2.00,  2.75, -1.00, 0.00,
                                  1.00, 3.00, 2.25,  -0.25, 2.50, 0.00};
    double pinchedGauss = 1e300;
    for (int gp = 0; gp < kGauss; ++gp)
        pinchedGauss = std::min(pinchedGauss, determinantAt(pinched, (gp & 1) ? q : -q,
                                                            (gp & 2) ? q : -q, (gp & 4) ? q : -q));
    const double pinchedCentre = determinantAt(pinched, 0.0, 0.0, 0.0);
    std::printf("     an element pinched at its centre: centre %.6f, smallest Gauss %.6f\n",
                pinchedCentre, pinchedGauss);
    expectTrue("the pinched element really is worst at its centre", pinchedCentre < pinchedGauss);
    expectNear("and that is the determinant `elementShape` reports",
               elementShape(pinched).smallestGauss, pinchedCentre, 1e-12);
}

// And the same statement at the level a caller sees: a mesh of collapsed elements
// solves, and solves *right*.
//
// The unit square as two triangles, extruded to a plate, under the linear field
// `u = (eps x, -nu eps y, -nu eps z)` prescribed on every boundary node with the one
// interior-most node free. That field is exact for uniaxial stress and it is in the
// trilinear span of both elements, so an element that passes the patch test
// reproduces it to round-off whatever its shape -- and an element whose Jacobian is
// wrong does not.
void testCollapsedElementSolvesAndPassesThePatchTest() {
    const StructuralMaterial steel = steelMaterial();
    const double t = 0.02, strain = 1e-4, nu = steel.poissonRatio;

    HexMesh mesh;
    // Five mid-surface points: the unit square's corners and its centre, so the
    // square is four triangles round a free interior node rather than two with
    // nothing free -- a patch with no free degree of freedom proves nothing.
    const double p[5][2] = {{0, 0}, {1, 0}, {1, 1}, {0, 1}, {0.37, 0.44}};
    for (const auto& q : p) {
        mesh.position.push_back(q[0]);
        mesh.position.push_back(q[1]);
        mesh.position.push_back(-t / 2);
        mesh.position.push_back(q[0]);
        mesh.position.push_back(q[1]);
        mesh.position.push_back(t / 2);
    }
    const auto triangle = [&](std::uint32_t a, std::uint32_t b) {
        // (a, b, centre, a): the fourth corner on the first, which is exactly what a
        // degenerate `PlatePanel` hands the section mesher.
        const std::uint32_t corner[4] = {a, b, 4u, a};
        for (int k = 0; k < 4; ++k) mesh.index.push_back(corner[k] * 2);
        for (int k = 0; k < 4; ++k) mesh.index.push_back(corner[k] * 2 + 1);
    };
    triangle(0, 1);
    triangle(1, 2);
    triangle(2, 3);
    triangle(3, 0);

    std::size_t collapsed = 0;
    for (std::size_t e = 0; e < mesh.elementCount(); ++e) {
        double nodes[kDof];
        mesh.gather(e, mesh.position, nodes);
        const ElementShape shape = elementShape(nodes);
        if (shape.collapsedNodes > 0) ++collapsed;
        expectTrue("every triangle of the patch is integrable", shape.integrable);
    }
    expectEqual("all four elements of the patch are collapsed",
                static_cast<long long>(collapsed), 4LL);

    // Prescribe the exact field on the four corners; the centre node pair is free.
    mesh.fixed.assign(mesh.nodeCount() * 3, 0u);
    mesh.prescribed.assign(mesh.nodeCount() * 3, 0.0);
    for (std::uint32_t n = 0; n < 8u; ++n) {
        for (int k = 0; k < 3; ++k) {
            mesh.fixed[n * 3 + static_cast<std::size_t>(k)] = 1u;
            const double x = mesh.position[n * 3 + static_cast<std::size_t>(k)];
            mesh.prescribed[n * 3 + static_cast<std::size_t>(k)] =
                (k == 0 ? strain : -nu * strain) * x;
        }
    }
    std::vector<double> load(mesh.nodeCount() * 3, 0.0), displacement;
    std::string problem;
    const bool solved =
        solveStatic(mesh, steel, Formulation::SolidShell, load, displacement, &problem);
    expectTrue("solveStatic accepts a mesh of collapsed elements: " + problem, solved);

    double worst = 0;
    for (std::size_t n = 8; n < mesh.nodeCount(); ++n)
        for (int k = 0; k < 3; ++k) {
            const double x = mesh.position[n * 3 + static_cast<std::size_t>(k)];
            worst = std::max(worst, std::abs(displacement[n * 3 + static_cast<std::size_t>(k)] -
                                             (k == 0 ? strain : -nu * strain) * x));
        }
    std::printf("     collapsed patch: free node off the exact field by %.2e m of %.2e\n", worst,
                strain);
    // Measured at 3.4e-19 m, which is 3.4e-15 of the field. Asserted a hundredfold
    // above that and still eleven decades below anything a mis-integrated wedge
    // could hide in -- a wrong Jacobian shows up in the *third* figure, not the
    // fifteenth.
    expectTrue("and reproduces the linear field exactly", worst < 1e-13 * strain);
    // Vacuity: the free node has to be somewhere the field is non-zero, or "exact"
    // would mean "zero equals zero".
    expectTrue("the free node is somewhere the field is not zero",
               std::abs(mesh.position[8 * 3]) > 0.1);

    // And the steel is all there: the four triangles are the unit square.
    double volume = 0;
    for (std::size_t e = 0; e < mesh.elementCount(); ++e) {
        double nodes[kDof], gauss[kGauss];
        mesh.gather(e, mesh.position, nodes);
        gaussVolumes(nodes, gauss);
        for (double v : gauss) volume += v;
    }
    expectNear("and the four wedges are exactly the plate they cover", volume, t, 1e-15 * t);
}

// The pressure load's *line of action*, not just its total. Both are closed forms
// -- the resultant is p times the face area vector and it acts through the face's
// area centroid -- and on a rectangular face an even quarter-each lumping gets both
// right, so the test is run on a face that is not rectangular.
void testPressureLoadResultant() {
    HexMesh mesh = makePlateMesh(1.0, 1.0, 0.2, 1, 1, 1);
    // Skew the top face into a general quadrilateral. Through `index`, because the
    // element's local node 5 is not the mesh's node 5.
    const std::size_t face[4] = {mesh.index[4], mesh.index[5], mesh.index[6], mesh.index[7]};
    mesh.position[face[1] * 3 + 0] = 1.4;
    mesh.position[face[2] * 3 + 1] = 1.6;
    mesh.position[face[3] * 3 + 0] = -0.2;
    const double pressure = 3.0e4;
    const std::vector<double> load = uniformPressureLoad(mesh, pressure);

    // Independent quadrature of the same face.
    static constexpr double kAbscissa[4] = {-0.8611363115940526, -0.3399810435848563,
                                            0.3399810435848563, 0.8611363115940526};
    static constexpr double kWeight[4] = {0.3478548451374538, 0.6521451548625461,
                                          0.6521451548625461, 0.3478548451374538};
    static constexpr double kCorner[4][2] = {{-1, -1}, {+1, -1}, {+1, +1}, {-1, +1}};
    double area = 0, firstMoment[2] = {0, 0};
    for (int a = 0; a < 4; ++a)
        for (int b = 0; b < 4; ++b) {
            const double s = kAbscissa[a], t = kAbscissa[b];
            double ds[4], dt[4], shape[4];
            for (int c = 0; c < 4; ++c) {
                shape[c] = 0.25 * (1 + s * kCorner[c][0]) * (1 + t * kCorner[c][1]);
                ds[c] = 0.25 * kCorner[c][0] * (1 + t * kCorner[c][1]);
                dt[c] = 0.25 * kCorner[c][1] * (1 + s * kCorner[c][0]);
            }
            double xs[3] = {0, 0, 0}, xt[3] = {0, 0, 0}, point[3] = {0, 0, 0};
            for (int c = 0; c < 4; ++c)
                for (int i = 0; i < 3; ++i) {
                    const double p =
                        mesh.position[face[c] * 3 + static_cast<std::size_t>(i)];
                    xs[i] += ds[c] * p;
                    xt[i] += dt[c] * p;
                    point[i] += shape[c] * p;
                }
            const double normalZ = xs[0] * xt[1] - xs[1] * xt[0];
            const double w = kWeight[a] * kWeight[b] * normalZ;
            area += w;
            for (int i = 0; i < 2; ++i) firstMoment[i] += w * point[i];
        }

    double totalZ = 0, centre[2] = {0, 0};
    for (std::size_t node = 0; node < mesh.nodeCount(); ++node) {
        const double fz = load[node * 3 + 2];
        totalZ += fz;
        for (int i = 0; i < 2; ++i)
            centre[i] += fz * mesh.position[node * 3 + static_cast<std::size_t>(i)];
    }
    expectNear("pressure resultant is p times the face area", totalZ, -pressure * area,
               1e-9 * pressure * std::abs(area));
    for (int i = 0; i < 2; ++i)
        expectNear("and acts through the face centroid, axis " + std::to_string(i),
                   centre[i] / totalZ, firstMoment[i] / area, 1e-9);

    // Non-vacuous: on this face the centroid is not the average of the corners.
    double cornerAverage = 0;
    for (int c = 0; c < 4; ++c) cornerAverage += mesh.position[face[c] * 3] / 4.0;
    expectTrue("the face is skewed enough for the distinction to exist",
               std::abs(cornerAverage - firstMoment[0] / area) > 1e-3);
}

// --- 1b. Eigenstrain ------------------------------------------------------------
//
// A strain that is not a stress. Everything here is an identity, because every one
// of these has an exact answer and a tolerance on any of them would hide the one
// defect the whole mechanism is exposed to: a subtraction that went the wrong way,
// or that happened in one of the two places it has to happen and not the other.

// A mesh whose elements are all distorted **in plane** -- the distortion a hull
// panel actually meshes to -- with every element still of constant thickness. That
// last part is load-bearing: a solid-shell whose thickness varies *within* the
// element cannot represent a free dilatation exactly, and
// `testTaperedElementCannotExpandFreely` below is where that is measured rather
// than tripped over.
//
// Interior nodes only, so the block stays a block and its boundary stays planar;
// the test needs a distorted *element*, not a distorted *body*.
HexMesh distortedBlock(int nx, int ny, int nz, double lx, double ly, double lz) {
    HexMesh mesh = makePlateMesh(lx, ly, lz, nx, ny, nz);
    std::mt19937 rng(20260809u);
    std::uniform_real_distribution<double> jitter(-0.22, 0.22);
    const double h[2] = {lx / nx, ly / ny};
    const double extent[2] = {lx, ly};
    // The in-plane position of a node, keyed by the (x, y) it started at, so that
    // every node in a through-thickness stack moves together and the columns stay
    // straight. Otherwise the elements taper and the test measures the element's
    // ANS sampling instead of the eigenstrain.
    std::map<std::pair<long long, long long>, std::pair<double, double>> moved;
    for (std::size_t n = 0; n < mesh.nodeCount(); ++n) {
        const double x0 = mesh.position[n * 3], y0 = mesh.position[n * 3 + 1];
        const std::pair<long long, long long> key{std::llround(x0 * 1e9),
                                                  std::llround(y0 * 1e9)};
        auto it = moved.find(key);
        if (it == moved.end()) {
            double x = x0, y = y0;
            if (x0 > 1e-12 && x0 < extent[0] - 1e-12) x += jitter(rng) * h[0];
            if (y0 > 1e-12 && y0 < extent[1] - 1e-12) y += jitter(rng) * h[1];
            it = moved.emplace(key, std::pair<double, double>{x, y}).first;
        }
        mesh.position[n * 3] = it->second.first;
        mesh.position[n * 3 + 1] = it->second.second;
    }
    return mesh;
}

// The largest |sigma| over every Gauss point of every element.
double worstStress(const HexMesh& mesh, const StructuralMaterial& material, Formulation form,
                   const std::vector<double>& displacement, const double* eigenstrain) {
    double worst = 0;
    for (std::size_t e = 0; e < mesh.elementCount(); ++e) {
        double nodes[kDof], u[kDof], stress[kGauss * 6];
        for (int a = 0; a < kNodes; ++a) {
            const std::size_t n = mesh.index[e * kNodes + static_cast<std::size_t>(a)];
            for (int i = 0; i < 3; ++i) {
                nodes[a * 3 + i] = mesh.position[n * 3 + static_cast<std::size_t>(i)];
                u[a * 3 + i] = displacement[n * 3 + static_cast<std::size_t>(i)];
            }
        }
        elementStress(nodes, u, material, form, stress, eigenstrain);
        for (int i = 0; i < kGauss * 6; ++i) worst = std::max(worst, std::abs(stress[i]));
    }
    return worst;
}

// Free expansion: the three constraints that remove rigid body motion and nothing
// else. Node at the origin fully, one along +x in y and z, one along +y in z --
// statically determinate, so the body is free to grow in every direction.
void supportFreely(HexMesh& mesh) {
    std::size_t origin = 0, alongX = 0, alongY = 0;
    double bestX = -1e30, bestY = -1e30, best = 1e30;
    for (std::size_t n = 0; n < mesh.nodeCount(); ++n) {
        const double x = mesh.position[n * 3], y = mesh.position[n * 3 + 1],
                     z = mesh.position[n * 3 + 2];
        if (x + y + z < best) { best = x + y + z; origin = n; }
    }
    const double oy = mesh.position[origin * 3 + 1], oz = mesh.position[origin * 3 + 2];
    for (std::size_t n = 0; n < mesh.nodeCount(); ++n) {
        const double x = mesh.position[n * 3], y = mesh.position[n * 3 + 1],
                     z = mesh.position[n * 3 + 2];
        if (std::abs(y - oy) < 1e-12 && std::abs(z - oz) < 1e-12 && x > bestX) {
            bestX = x;
            alongX = n;
        }
        if (std::abs(x - mesh.position[origin * 3]) < 1e-12 && std::abs(z - oz) < 1e-12 &&
            y > bestY) {
            bestY = y;
            alongY = n;
        }
    }
    mesh.pin(origin, 0);
    mesh.pin(origin, 1);
    mesh.pin(origin, 2);
    mesh.pin(alongX, 1);
    mesh.pin(alongX, 2);
    mesh.pin(alongY, 2);
}

void testFreeExpansionCarriesNoStress() {
    const StructuralMaterial steel = steelMaterial();
    const double e = 5.1984e-3;  // the EN elongation at 400 C, so the scale is a real one
    const double eigen[6] = {e, e, e, 0, 0, 0};

    // What the stress would be if the eigenstrain were simply *not* subtracted --
    // the whole point of asserting zero rather than something small.
    const double nu = steel.poissonRatio, ee = steel.youngsModulus;
    const double bulkStress = ee / (1.0 - 2.0 * nu) * e;  // (3 lambda + 2 mu) e

    for (int distorted = 0; distorted < 2; ++distorted) {
        HexMesh mesh = distorted ? distortedBlock(3, 2, 2, 0.9, 0.6, 0.048)
                                 : makePlateMesh(0.9, 0.6, 0.048, 3, 2, 2);
        supportFreely(mesh);
        const std::string tag = distorted ? " (distorted)" : " (rectangular)";

        std::vector<double> eigenstrain(mesh.elementCount() * 6);
        for (std::size_t el = 0; el < mesh.elementCount(); ++el)
            for (int i = 0; i < 6; ++i) eigenstrain[el * 6 + i] = eigen[i];

        const std::vector<double> load =
            thermalLoad(mesh, steel, Formulation::SolidShell, eigenstrain);
        expectEqual("the thermal load has one value per degree of freedom" + tag,
                    static_cast<long long>(load.size()),
                    static_cast<long long>(mesh.nodeCount() * 3));

        // The load is in equilibrium with itself: heating a free body applies no
        // net force and no net moment to it.
        double force[3] = {0, 0, 0}, moment[3] = {0, 0, 0}, biggest = 0;
        for (std::size_t n = 0; n < mesh.nodeCount(); ++n) {
            const double* f = &load[n * 3];
            const double* p = &mesh.position[n * 3];
            for (int i = 0; i < 3; ++i) {
                force[i] += f[i];
                biggest = std::max(biggest, std::abs(f[i]));
            }
            moment[0] += p[1] * f[2] - p[2] * f[1];
            moment[1] += p[2] * f[0] - p[0] * f[2];
            moment[2] += p[0] * f[1] - p[1] * f[0];
        }
        expectTrue("the thermal load is not itself trivial" + tag, biggest > 1e5);
        for (int i = 0; i < 3; ++i) {
            expectNear("the thermal load carries no net force, axis " + std::to_string(i) + tag,
                       force[i], 0.0, 1e-9 * biggest);
            expectNear("nor any net moment, axis " + std::to_string(i) + tag, moment[i], 0.0,
                       1e-9 * biggest * 0.9);
        }

        std::vector<double> displacement;
        std::string problem;
        expectTrue("the free-expansion solve succeeds" + tag,
                   solveStatic(mesh, steel, Formulation::SolidShell, load, displacement,
                               &problem));

        // The exact solution is a uniform dilatation about the pinned origin,
        // u = e (x - x0), and it is in the trilinear space, so the discrete answer
        // is that field to rounding.
        std::size_t origin = 0;
        double best = 1e30;
        for (std::size_t n = 0; n < mesh.nodeCount(); ++n) {
            const double s = mesh.position[n * 3] + mesh.position[n * 3 + 1] +
                             mesh.position[n * 3 + 2];
            if (s < best) { best = s; origin = n; }
        }
        double worstDisplacement = 0, span = 0;
        for (std::size_t n = 0; n < mesh.nodeCount(); ++n)
            for (int i = 0; i < 3; ++i) {
                const double want =
                    e * (mesh.position[n * 3 + static_cast<std::size_t>(i)] -
                         mesh.position[origin * 3 + static_cast<std::size_t>(i)]);
                worstDisplacement = std::max(
                    worstDisplacement,
                    std::abs(displacement[n * 3 + static_cast<std::size_t>(i)] - want));
                span = std::max(span, std::abs(want));
            }
        expectTrue("the body actually expanded" + tag, span > 1e-6);
        expectNear("free expansion is exactly the uniform dilatation" + tag, worstDisplacement,
                   0.0, 1e-11 * span);

        // **The assertion this whole mechanism exists for.**
        const double zero = worstStress(mesh, steel, Formulation::SolidShell, displacement,
                                        eigen);
        expectNear("a freely expanding body carries no stress at all" + tag, zero, 0.0,
                   1e-9 * bulkStress);

        // Non-vacuous, twice over. Reading the same displacement *without* the
        // eigenstrain gives the full bulk stress -- so the test fails if the
        // subtraction in `elementStress` is dropped -- and solving without the
        // thermal load gives no displacement at all, so it fails if the load is.
        const double ignored =
            worstStress(mesh, steel, Formulation::SolidShell, displacement, nullptr);
        expectNear("and it is 2.68 GPa if the eigenstrain is not subtracted" + tag, ignored,
                   bulkStress, 1e-6 * bulkStress);

        std::vector<double> none(mesh.nodeCount() * 3, 0.0), unloaded;
        expectTrue("the unloaded solve succeeds" + tag,
                   solveStatic(mesh, steel, Formulation::SolidShell, none, unloaded, &problem));
        const double withoutLoad =
            worstStress(mesh, steel, Formulation::SolidShell, unloaded, eigen);
        expectNear("and it is the same magnitude, negated, if the load is dropped" + tag,
                   withoutLoad, bulkStress, 1e-6 * bulkStress);
    }
}

// Six element geometries: what the mesher makes, and two that it does not.
struct NamedElement {
    const char* name;
    double node[kDof];
};

const NamedElement kGeometries[] = {
    {"rectangular",
     {0, 0, 0, 0.60, 0, 0, 0.60, 0.45, 0, 0, 0.45, 0, 0, 0, 0.012, 0.60, 0, 0.012, 0.60, 0.45,
      0.012, 0, 0.45, 0.012}},
    {"in-plane trapezoid",
     {0, 0, 0, 0.62, 0.03, 0, 0.55, 0.48, 0, 0.04, 0.44, 0, 0, 0, 0.012, 0.62, 0.03, 0.012, 0.55,
      0.48, 0.012, 0.04, 0.44, 0.012}},
    {"warped mid-surface, 1 mm in 12",
     {0, 0, 0.000, 0.60, 0, 0.001, 0.60, 0.45, 0.000, 0, 0.45, 0.001, 0, 0, 0.012, 0.60, 0, 0.013,
      0.60, 0.45, 0.012, 0, 0.45, 0.013}},
    {"uniform director tilt",
     {0, 0, 0, 0.60, 0, 0, 0.60, 0.45, 0, 0, 0.45, 0, 0.004, 0.003, 0.012, 0.604, 0.003, 0.012,
      0.604, 0.453, 0.012, 0.004, 0.453, 0.012}},
};

void testEigenstrainAndTheEnhancedModes() {
    // **`int G^T C eps* dV` is identically zero for a constant eigenstrain, on every
    // geometry**, because `int G^T dV = 0` is the Simo-Rifai condition that makes
    // EAS pass the patch test. That was not the expectation when the term was
    // written, so it is asserted rather than argued -- and asserted on the integral
    // itself, so that an enhanced basis which stopped satisfying it is caught here
    // and not in a caller.
    const StructuralMaterial steel = steelMaterial();
    const double e = 5.1984e-3;  // the EN elongation at 400 C
    const double eigen[6] = {e, e, e, 0, 0, 0};
    const double bulkStress = steel.youngsModulus / (1.0 - 2.0 * steel.poissonRatio) * e;

    for (const NamedElement& g : kGeometries) {
        const std::string tag = std::string(" -- ") + g.name;

        // `int G^T dV = 0`, read out of the element's own forms rather than
        // re-derived, one component at a time.
        RestForms forms;
        expectTrue("the forms build" + tag,
                   computeRestForms(g.node, Formulation::SolidShell, forms));
        double volume = 0;
        for (int gp = 0; gp < kGauss; ++gp) volume += forms.weight[gp];
        double worstIntegral = 0, worstMode = 0;
        for (int k = 0; k < forms.easCount; ++k)
            for (int i = 0; i < 6; ++i) {
                double integral = 0, magnitude = 0;
                for (int gp = 0; gp < kGauss; ++gp) {
                    integral += forms.weight[gp] * forms.g[gp][i][k];
                    magnitude += forms.weight[gp] * std::abs(forms.g[gp][i][k]);
                }
                worstIntegral = std::max(worstIntegral, std::abs(integral));
                worstMode = std::max(worstMode, magnitude);
            }
        expectTrue("the enhanced modes are not themselves zero" + tag, worstMode > 1e-6);
        expectNear("int G dV vanishes, which is why the enhanced thermal term does" + tag,
                   worstIntegral, 0.0, 1e-14 * worstMode);
        expectTrue("and the element has a volume to integrate over" + tag, volume > 1e-6);

        // Free expansion, exactly. The displacement is the closed-form dilatation,
        // so no solve stands between the element and the assertion.
        double u[kDof], stress[kGauss * 6];
        for (int i = 0; i < kDof; ++i) u[i] = e * g.node[i];
        elementStress(g.node, u, steel, Formulation::SolidShell, stress, eigen);
        double worst = 0;
        for (int i = 0; i < kGauss * 6; ++i) worst = std::max(worst, std::abs(stress[i]));
        expectNear("free expansion is exactly stress free" + tag, worst, 0.0,
                   1e-13 * bulkStress);
    }
}

void testTaperedElementCannotExpandFreely() {
    // **The one geometry that is not exact, and the eigenstrain is not why.** A
    // solid-shell whose thickness varies within the element samples
    // `E_zeta,zeta` at four in-plane corners and interpolates it bilinearly, which
    // is exact only while `|g_zeta|` is constant. The plain patch test -- a linear
    // displacement field, no eigenstrain anywhere -- already fails on it, and by
    // more than the free-expansion case does. Both are measured so the day the
    // mesher produces one, the size is written down rather than discovered.
    //
    // Nothing reaches it today: `zone.cpp` extrudes one patch at one thickness by
    // construction and warns when a thickness seam is crossed.
    const StructuralMaterial steel = steelMaterial();
    const double e = 5.1984e-3;
    const double eigen[6] = {e, e, e, 0, 0, 0};
    const double bulkStress = steel.youngsModulus / (1.0 - 2.0 * steel.poissonRatio) * e;
    // 10 mm at one corner to 15 mm at another: a +-20% taper, far past anything a
    // hull mesh would carry.
    const double tapered[kDof] = {0,    0,    0,     0.60, 0,    0,     0.60, 0.45,
                                  0,    0,    0.45,  0,    0,    0,     0.014, 0.60,
                                  0,    0.010, 0.60, 0.45, 0.015, 0,    0.45,  0.011};

    double u[kDof], stress[kGauss * 6];
    for (int i = 0; i < kDof; ++i) u[i] = e * tapered[i];
    elementStress(tapered, u, steel, Formulation::SolidShell, stress, eigen);
    double worst = 0;
    for (int i = 0; i < kGauss * 6; ++i) worst = std::max(worst, std::abs(stress[i]));
    expectNear("a +-20% tapered element leaves 1.75% of the bulk stress under free expansion",
               worst / bulkStress, 0.01748, 5e-4);

    // The same element, no eigenstrain at all, under a general linear displacement:
    // the stress ought to be constant over the eight points and is not.
    const double gradient[3][3] = {{1.1e-3, 3.0e-4, -2.0e-4},
                                   {5.0e-4, -7.0e-4, 1.0e-4},
                                   {2.0e-4, -4.0e-4, 9.0e-4}};
    double linear[kDof], plain[kGauss * 6];
    for (int a = 0; a < kNodes; ++a)
        for (int i = 0; i < 3; ++i)
            linear[a * 3 + i] = gradient[i][0] * tapered[a * 3] +
                                gradient[i][1] * tapered[a * 3 + 1] +
                                gradient[i][2] * tapered[a * 3 + 2];
    elementStress(tapered, linear, steel, Formulation::SolidShell, plain, nullptr);
    double spread = 0, magnitude = 0;
    for (int i = 0; i < 6; ++i) {
        double lo = plain[i], hi = plain[i];
        for (int gp = 1; gp < kGauss; ++gp) {
            lo = std::min(lo, plain[gp * 6 + i]);
            hi = std::max(hi, plain[gp * 6 + i]);
        }
        spread = std::max(spread, hi - lo);
        magnitude = std::max(magnitude, std::abs(hi));
    }
    expectNear("and the plain patch test spreads by 1.09% on the same element, with no"
               " eigenstrain in sight",
               spread / magnitude, 0.01089, 5e-4);

    // Non-vacuous: the same two checks on a constant-thickness element are at the
    // rounding floor, so the taper is what is being measured.
    const double even[kDof] = {0,    0,    0,     0.60, 0,    0,     0.60, 0.45,
                               0,    0,    0.45,  0,    0,    0,     0.012, 0.60,
                               0,    0.012, 0.60, 0.45, 0.012, 0,    0.45,  0.012};
    double ue[kDof], se[kGauss * 6];
    for (int i = 0; i < kDof; ++i) ue[i] = e * even[i];
    elementStress(even, ue, steel, Formulation::SolidShell, se, eigen);
    double evenWorst = 0;
    for (int i = 0; i < kGauss * 6; ++i) evenWorst = std::max(evenWorst, std::abs(se[i]));
    expectNear("the same element at constant thickness is exact", evenWorst / bulkStress, 0.0,
               1e-13);
}

void testRestraintAndNotTemperatureIsWhatGeneratesStress() {
    // A bar heated uniformly and held at both ends to a *prescribed* extension. The
    // degree of restraint R = 1 - delta / (eps* L) sweeps the whole spectrum from
    // free to fully restrained, and the closed form is exact at every point:
    //
    //     sigma_xx = -R E eps*
    //
    // -- because eps_xx = (1 - R) eps*, the lateral faces are free so sigma_yy and
    // sigma_zz vanish, and for that state the uniaxial modulus is exactly E. This
    // is the statement that a temperature on its own carries no stress and that
    // the whole of the effect is the boundary condition.
    const StructuralMaterial steel = steelMaterial();
    const double lx = 2.40, ly = 0.20, lz = 0.012;
    const double e = 1.842374e-3;  // the EN elongation at the 164.63 C anchor
    const double eigen[6] = {e, e, e, 0, 0, 0};

    for (double restraint : {1.0, 0.75, 0.5, 0.25, 0.0}) {
        HexMesh mesh = makePlateMesh(lx, ly, lz, 4, 2, 1);
        const double delta = (1.0 - restraint) * e * lx;
        // Both x faces held: the -x face at zero, the +x face at `delta`. Free in y
        // and z apart from three pins that remove the remaining rigid motions --
        // and the third of them is not optional. Pinning only the corner in y and z
        // leaves the bar free to **spin about its own axis**, which the x-face pins
        // do not touch, and `solveStatic` refuses the singular system rather than
        // returning a plausible answer. That refusal is what caught it.
        std::size_t corner = mesh.nodeCount(), offAxis = mesh.nodeCount();
        for (std::size_t n = 0; n < mesh.nodeCount(); ++n) {
            const double x = mesh.position[n * 3], y = mesh.position[n * 3 + 1],
                         z = mesh.position[n * 3 + 2];
            if (x < 1e-12) mesh.pin(n, 0, 0.0);
            if (x > lx - 1e-12) mesh.pin(n, 0, delta);
            if (x < 1e-12 && y < 1e-12 && z < 1e-12) corner = n;
            if (x < 1e-12 && y > ly - 1e-12 && z < 1e-12) offAxis = n;
        }
        expectTrue("both support nodes were found",
                   corner < mesh.nodeCount() && offAxis < mesh.nodeCount());
        mesh.pin(corner, 1, 0.0);
        mesh.pin(corner, 2, 0.0);
        mesh.pin(offAxis, 2, 0.0);  // kills rotation about x; both nodes sit at z = 0,
                                    // where the exact solution has u_z = 0 anyway

        std::vector<double> eigenstrain(mesh.elementCount() * 6);
        for (std::size_t el = 0; el < mesh.elementCount(); ++el)
            for (int i = 0; i < 6; ++i) eigenstrain[el * 6 + i] = eigen[i];
        const std::vector<double> load =
            thermalLoad(mesh, steel, Formulation::SolidShell, eigenstrain);

        std::vector<double> displacement;
        std::string problem;
        expectTrue("the restrained-bar solve succeeds",
                   solveStatic(mesh, steel, Formulation::SolidShell, load, displacement,
                               &problem));

        const double want = -restraint * steel.youngsModulus * e;
        double worstAxial = 0, worstLateral = 0;
        for (std::size_t el = 0; el < mesh.elementCount(); ++el) {
            double nodes[kDof], u[kDof], stress[kGauss * 6];
            for (int a = 0; a < kNodes; ++a) {
                const std::size_t n = mesh.index[el * kNodes + static_cast<std::size_t>(a)];
                for (int i = 0; i < 3; ++i) {
                    nodes[a * 3 + i] = mesh.position[n * 3 + static_cast<std::size_t>(i)];
                    u[a * 3 + i] = displacement[n * 3 + static_cast<std::size_t>(i)];
                }
            }
            elementStress(nodes, u, steel, Formulation::SolidShell, stress, eigen);
            for (int gp = 0; gp < kGauss; ++gp) {
                worstAxial = std::max(worstAxial, std::abs(stress[gp * 6] - want));
                for (int i = 1; i < 6; ++i)
                    worstLateral = std::max(worstLateral, std::abs(stress[gp * 6 + i]));
            }
        }
        const double scale = steel.youngsModulus * e;
        expectNear("a bar restrained to R = " + std::to_string(restraint) +
                       " carries exactly -R E eps*",
                   worstAxial, 0.0, 1e-9 * scale);
        expectNear("and nothing at all on the free faces, R = " + std::to_string(restraint),
                   worstLateral, 0.0, 1e-9 * scale);
    }

    // Non-vacuous: full restraint is 379 MPa of compression and free is zero, so
    // the sweep is not measuring the same thing five times.
    expectNear("full restraint at this eigenstrain is 379.5 MPa", steel.youngsModulus * e * 1e-6,
               379.53, 0.05);
}

void testNullEigenstrainIsTheSameArithmetic() {
    // The exact control. A zero eigenstrain and no eigenstrain must be **bit**
    // identical, on every output of both entry points, or an unheated ship is not
    // the ship master computes.
    const StructuralMaterial steel = steelMaterial();
    const double zeroEigen[6] = {0, 0, 0, 0, 0, 0};
    // Deliberately the *awkward* geometry here -- warped, tapered and sheared. Bit
    // identity has to hold on an element the arithmetic is hardest on, and this one
    // is not required to expand freely.
    const double nodes[kDof] = {0.00, 0.00, 0.000, 0.62, 0.03, 0.001, 0.55, 0.48, -0.002,
                                0.04, 0.44, 0.003, 0.02, 0.05, 0.014, 0.60, 0.01, 0.013,
                                0.58, 0.46, 0.015, 0.01, 0.41, 0.012};
    // A displacement with bending, shear and stretch in it, so every row of B and
    // every enhanced mode is exercised.
    double current[kDof];
    std::mt19937 rng(770213u);
    std::uniform_real_distribution<double> jitter(-2.0e-4, 2.0e-4);
    for (int i = 0; i < kDof; ++i) current[i] = nodes[i] + jitter(rng);

    double u[kDof];
    for (int i = 0; i < kDof; ++i) u[i] = current[i] - nodes[i];
    double withNull[kGauss * 6], withZero[kGauss * 6];
    elementStress(nodes, u, steel, Formulation::SolidShell, withNull, nullptr);
    elementStress(nodes, u, steel, Formulation::SolidShell, withZero, zeroEigen);
    bool identical = true, moved = false;
    for (int i = 0; i < kGauss * 6; ++i) {
        identical = identical && withNull[i] == withZero[i];
        moved = moved || std::abs(withNull[i]) > 1.0;
    }
    expectTrue("elementStress at zero eigenstrain is bit-identical to no eigenstrain",
               identical);
    expectTrue("and the element was carrying stress, so that is not vacuous", moved);

    const plasticity::Material flow = plasticity::shipSteel();
    ElementPlasticState stateNull, stateZero;
    initialisePlasticState(nodes, flow, stateNull);
    initialisePlasticState(nodes, flow, stateZero);
    double forceNull[kDof], forceZero[kDof], stressNull[kGauss * 6], stressZero[kGauss * 6];
    const PlasticUpdate a = elementPlasticUpdate(nodes, current, flow, Formulation::SolidShell,
                                                 stateNull, forceNull, stressNull, nullptr);
    const PlasticUpdate b = elementPlasticUpdate(nodes, current, flow, Formulation::SolidShell,
                                                 stateZero, forceZero, stressZero, zeroEigen);
    bool same = a.iterations == b.iterations && a.yieldedPoints == b.yieldedPoints &&
                a.dissipation == b.dissipation && a.enhancedWork == b.enhancedWork;
    for (int i = 0; i < kDof; ++i) same = same && forceNull[i] == forceZero[i];
    for (int i = 0; i < kGauss * 6; ++i) same = same && stressNull[i] == stressZero[i];
    for (int gp = 0; gp < kGauss; ++gp) {
        same = same && stateNull.point[gp].equivalentPlasticStrain ==
                           stateZero.point[gp].equivalentPlasticStrain;
        for (int i = 0; i < plasticity::kVoigt; ++i)
            same = same && stateNull.point[gp].plasticStrain[i] ==
                               stateZero.point[gp].plasticStrain[i];
    }
    for (int k = 0; k < kEas; ++k) same = same && stateNull.enhanced[k] == stateZero.enhanced[k];
    expectTrue("elementPlasticUpdate is bit-identical too, history included", same);
    expectTrue("and it did some work, so that is not vacuous either",
               a.enhancedWork > 0.0 || a.yieldedPoints > 0);
}

void testEigenstrainThroughThePlasticUpdate() {
    // The eigenstrain has to reach the *nonlinear* path as well, and there the
    // enhanced parameters are found by Newton on `r(alpha) = 0` rather than by a
    // closed form -- so if the subtraction is in the right place the enhanced modes
    // need no term of their own. Two checks, both exact.
    const plasticity::Material flow = plasticity::shipSteel();
    const StructuralMaterial steel = steelMaterial();
    // In-plane trapezoid at constant thickness: what a hull panel meshes to, and
    // the geometry on which free expansion is exact.
    const double nodes[kDof] = {0.00, 0.00, 0.000, 0.62, 0.03, 0.000, 0.55, 0.48, 0.000,
                                0.04, 0.44, 0.000, 0.00, 0.00, 0.012, 0.62, 0.03, 0.012,
                                0.55, 0.48, 0.012, 0.04, 0.44, 0.012};

    // (a) An element held at its rest position and given an eigenstrain is the same
    //     element displaced by -eps* x with no eigenstrain: the constitutive law
    //     sees the same mechanical strain either way, so it must answer the same.
    const double e = 4.0e-4;  // elastic, so the comparison is not about the return map
    const double eigen[6] = {e, e, e, 0, 0, 0};
    double shrunk[kDof];
    for (int i = 0; i < kDof; ++i) shrunk[i] = nodes[i] * (1.0 - e);

    ElementPlasticState heated, pulled;
    initialisePlasticState(nodes, flow, heated);
    initialisePlasticState(nodes, flow, pulled);
    double fHeated[kDof], fPulled[kDof], sHeated[kGauss * 6], sPulled[kGauss * 6];
    elementPlasticUpdate(nodes, nodes, flow, Formulation::SolidShell, heated, fHeated, sHeated,
                         eigen);
    elementPlasticUpdate(nodes, shrunk, flow, Formulation::SolidShell, pulled, fPulled, sPulled,
                         nullptr);
    double worst = 0, scale = 0;
    for (int i = 0; i < kGauss * 6; ++i) {
        worst = std::max(worst, std::abs(sHeated[i] - sPulled[i]));
        scale = std::max(scale, std::abs(sHeated[i]));
    }
    expectTrue("the heated element is actually stressed", scale > 1e8);
    // Not bit-identical: the two reach the same mechanical strain by different
    // arithmetic -- one subtracts a constant from B u, the other scales the
    // geometry -- and the co-rotation of a shrunk element is not exactly the
    // identity. Nine digits is what that costs and it is measured, not chosen.
    expectNear("heating an element is contracting it, to nine digits", worst, 0.0,
               2e-9 * scale);

    // (b) The stress of an element under free expansion, through the plastic path,
    //     is exactly zero -- the same identity as the elastic one, and the check
    //     that the Newton found alpha = 0 rather than something that happened to
    //     be small.
    //
    //     Both bounds are taken against the **same element fully restrained**,
    //     which is `heated` above: it is the only reference that gets the units and
    //     the magnitude right for a force as well as for a stress, and it makes the
    //     assertion "the free element is 1e-11 of the restrained one" rather than
    //     "small against a constant somebody chose". Measured: 3.1e-13 in both,
    //     which is the rounding of a 24-term dot product against a 206 MPa stress
    //     and 13.5 MN of force, so 1e-11 is thirty times the floor and not more.
    double grown[kDof];
    for (int i = 0; i < kDof; ++i) grown[i] = nodes[i] * (1.0 + e);
    ElementPlasticState free;
    initialisePlasticState(nodes, flow, free);
    double fFree[kDof], sFree[kGauss * 6];
    elementPlasticUpdate(nodes, grown, flow, Formulation::SolidShell, free, fFree, sFree,
                         eigen);
    double freeStress = 0, freeForce = 0, heldStress = 0, heldForce = 0;
    for (int i = 0; i < kGauss * 6; ++i) {
        freeStress = std::max(freeStress, std::abs(sFree[i]));
        heldStress = std::max(heldStress, std::abs(sHeated[i]));
    }
    for (int i = 0; i < kDof; ++i) {
        freeForce = std::max(freeForce, std::abs(fFree[i]));
        heldForce = std::max(heldForce, std::abs(fHeated[i]));
    }
    expectTrue("the restrained reference carries a real stress and a real force",
               heldStress > 1e8 && heldForce > 1e6);
    expectNear("free expansion through the plastic path carries no stress",
               freeStress / heldStress, 0.0, 1e-11);
    expectNear("and therefore no internal force", freeForce / heldForce, 0.0, 1e-11);
    // The same element is the reference for both, so a units slip in one would show
    // as a ratio that is not a pure number.
    expectNear("the two residuals are the same rounding, not two different ones",
               freeForce / heldForce, freeStress / heldStress, 1e-12);
}

// --- 2. Rigid body motion under the co-rotational formulation ------------------

void testFiniteRotationCarriesNoForce() {
    const StructuralMaterial steel = steelMaterial();
    // Distorted, so the result cannot come from a symmetry of the shape.
    double rest[kDof] = {0.01,  -0.02, -0.004, 0.052, 0.003, -0.006, 0.048, 0.047,
                         -0.005, 0.002, 0.051, -0.003, -0.001, 0.004, 0.006, 0.049,
                         -0.002, 0.005, 0.053, 0.052, 0.004, 0.003, 0.048, 0.007};
    double k[kDof * kDof];
    elementStiffness(rest, steel, Formulation::SolidShell, k);

    // The control: what a genuine strain of 1e-6 produces, so "zero" is measured
    // against a force this element is capable of and not against nothing.
    const double probeStrain = 1.0e-6;
    double strained[kDof];
    for (int i = 0; i < kDof; ++i)
        strained[i] = rest[i] * (i % 3 == 0 ? 1.0 + probeStrain : 1.0);
    double reference[kDof];
    internalForce(k, rest, strained, reference);
    double referenceNorm = 0;
    for (int i = 0; i < kDof; ++i) referenceNorm = std::max(referenceNorm, std::abs(reference[i]));
    expectTrue("a 1e-6 strain produces a measurable force", referenceNorm > 1.0);

    // Zero is not free. f = R K (R^T x - X) forms K times a vector whose entries
    // are the rigid translation, and the rows of K cancel to zero only to
    // rounding, so the floor is |K| |u| eps. Asserting a fixed small number
    // instead would silently be a statement about how far the element was moved.
    double stiffnessScale = 0;
    for (int i = 0; i < kDof * kDof; ++i)
        stiffnessScale = std::max(stiffnessScale, std::abs(k[i]));
    const double translation[3] = {0.30, -0.75, 1.125};
    double translationScale = 0;
    for (double component : translation)
        translationScale = std::max(translationScale, std::abs(component));
    const double floor =
        stiffnessScale * translationScale * std::numeric_limits<double>::epsilon() * kDof;

    for (double angle : {0.001, 0.3, 1.2, 2.5, 3.0}) {
        // Rotation about a skew axis, by Rodrigues.
        const double axis[3] = {0.5773502691896258, -0.5773502691896258, 0.5773502691896258};
        const double c = std::cos(angle), s = std::sin(angle);
        double r[3][3];
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) r[i][j] = (i == j ? c : 0.0) + (1.0 - c) * axis[i] * axis[j];
        r[0][1] += -s * axis[2];
        r[0][2] += s * axis[1];
        r[1][0] += s * axis[2];
        r[1][2] += -s * axis[0];
        r[2][0] += -s * axis[1];
        r[2][1] += s * axis[0];

        double moved[kDof];
        for (int a = 0; a < kNodes; ++a)
            for (int i = 0; i < 3; ++i) {
                double sum = translation[i];
                for (int j = 0; j < 3; ++j) sum += r[i][j] * rest[a * 3 + j];
                moved[a * 3 + i] = sum;
            }

        double force[kDof];
        internalForce(k, rest, moved, force);
        double worst = 0;
        for (int i = 0; i < kDof; ++i) worst = std::max(worst, std::abs(force[i]));
        // Two bounds, because they catch different things. Against the rounding
        // floor: a polar decomposition stopped short -- fem.cpp's fixed four
        // Higham iterations, say -- leaves a residue that shows here and nowhere
        // else. Against a real strain: the implied spurious strain must be
        // negligible against anything the element is asked to resolve.
        expectTrue("a finite rotation of " + std::to_string(angle) +
                       " rad carries no force (rounding floor)",
                   worst < 50.0 * floor);
        expectTrue("a finite rotation of " + std::to_string(angle) +
                       " rad implies under 1e-13 of spurious strain",
                   worst / referenceNorm * probeStrain < 1e-13);
        if (angle > 1.0 && angle < 1.5)
            std::printf("     finite rotation %.1f rad: force %.2e N against a rounding floor of"
                        " %.2e N, %.2e of a 1e-6 strain\n",
                        angle, worst, floor, worst / referenceNorm);

        // And the rotation it recovered is the one applied.
        double recovered[9];
        elementRotation(rest, moved, recovered);
        double worstRotation = 0;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                worstRotation = std::max(worstRotation, std::abs(recovered[j * 3 + i] - r[i][j]));
        expectTrue("the polar decomposition recovers the rotation applied", worstRotation < 1e-12);
    }
}

// Frame indifference, which is a stronger statement than "a rigid motion carries
// no force" and catches what that one cannot. Rotating an *undeformed* element
// hands the polar decomposition a deformation gradient that is already a rotation,
// and Higham's iteration returns it unchanged on the first step -- so a
// decomposition truncated to one iteration, or stopped at 1e-6, passes. Rotating a
// *deformed* element does not let it off: F = R U with U away from the identity,
// and the recovered force has to be the deformed one rotated, exactly.
void testFrameIndifference() {
    const StructuralMaterial steel = steelMaterial();
    double rest[kDof] = {0.01,   -0.02, -0.004, 0.052, 0.003,  -0.006, 0.048, 0.047,
                         -0.005, 0.002, 0.051,  -0.003, -0.001, 0.004, 0.006, 0.049,
                         -0.002, 0.005, 0.053,  0.052, 0.004,  0.003, 0.048, 0.007};
    double k[kDof * kDof];
    elementStiffness(rest, steel, Formulation::SolidShell, k);

    // A genuinely large deformation, not a small one: 25% stretch and 15% shear, so
    // the stretch factor U is a long way from the identity and Higham's iteration
    // has real work to do. A small strain would leave F close to a rotation
    // already, and then a decomposition that barely iterated would pass -- which is
    // exactly what happened to the first version of this test. The identity being
    // checked holds for any deformation, so there is no reason to pick a polite one.
    double deformed[kDof];
    for (int a = 0; a < kNodes; ++a) {
        const double x = rest[a * 3], y = rest[a * 3 + 1], z = rest[a * 3 + 2];
        deformed[a * 3 + 0] = x * 1.25 + 0.15 * y - 0.08 * z;
        deformed[a * 3 + 1] = y * 0.82 + 0.11 * z;
        deformed[a * 3 + 2] = z * 1.13 + 0.17 * x;
    }
    double reference[kDof];
    internalForce(k, rest, deformed, reference);
    double referenceNorm = 0;
    for (int i = 0; i < kDof; ++i) referenceNorm = std::max(referenceNorm, std::abs(reference[i]));
    expectTrue("the deformed element carries a real force", referenceNorm > 1.0);

    for (double angle : {0.4, 1.7, 2.9}) {
        const double axis[3] = {0.2672612419124244, 0.5345224838248488, 0.8017837257372732};
        const double c = std::cos(angle), s = std::sin(angle);
        double r[3][3];
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) r[i][j] = (i == j ? c : 0.0) + (1.0 - c) * axis[i] * axis[j];
        r[0][1] += -s * axis[2];
        r[0][2] += s * axis[1];
        r[1][0] += s * axis[2];
        r[1][2] += -s * axis[0];
        r[2][0] += -s * axis[1];
        r[2][1] += s * axis[0];

        double moved[kDof];
        for (int a = 0; a < kNodes; ++a)
            for (int i = 0; i < 3; ++i) {
                double sum = (i == 0 ? 2.0 : i == 1 ? -3.0 : 5.0);
                for (int j = 0; j < 3; ++j) sum += r[i][j] * deformed[a * 3 + j];
                moved[a * 3 + i] = sum;
            }
        double rotatedForce[kDof];
        internalForce(k, rest, moved, rotatedForce);

        double worst = 0;
        for (int a = 0; a < kNodes; ++a)
            for (int i = 0; i < 3; ++i) {
                double want = 0;
                for (int j = 0; j < 3; ++j) want += r[i][j] * reference[a * 3 + j];
                worst = std::max(worst, std::abs(rotatedForce[a * 3 + i] - want));
            }
        expectTrue("rotating a deformed element rotates its force, " + std::to_string(angle) +
                       " rad",
                   worst < 1e-12 * referenceNorm);
        if (angle > 1.5 && angle < 2.0)
            std::printf("     a 25%% stretch rotated %.1f rad: force differs from the rotated"
                        " force by %.2e of it\n",
                        angle, worst / referenceNorm);
    }
}

// The internal force's sign and its equilibrium. Neither shows up in any test that
// only asks for zero: a sign flip turns the solver into an anti-spring and every
// rigid-body check still passes.
void testInternalForceSignAndEquilibrium() {
    const StructuralMaterial steel = steelMaterial();
    double rest[kDof] = {0, 0, -0.005, 0.05, 0, -0.005, 0.05, 0.05, -0.005, 0, 0.05, -0.005,
                         0, 0, 0.005,  0.05, 0, 0.005,  0.05, 0.05, 0.005,  0, 0.05, 0.005};
    double k[kDof * kDof];
    elementStiffness(rest, steel, Formulation::SolidShell, k);

    // Stretched 1e-4 along x.
    double stretched[kDof];
    for (int i = 0; i < kDof; ++i) stretched[i] = rest[i] * (i % 3 == 0 ? 1.0 + 1.0e-4 : 1.0);
    double force[kDof];
    internalForce(k, rest, stretched, force);

    // A stretched elastic body pulls its ends back together.
    double pullBack = 0, total[3] = {0, 0, 0}, moment[3] = {0, 0, 0};
    for (int a = 0; a < kNodes; ++a) {
        if (rest[a * 3] > 0.025) pullBack += force[a * 3];
        for (int i = 0; i < 3; ++i) total[i] += force[a * 3 + i];
        const double* x = &rest[a * 3];
        const double* f = &force[a * 3];
        moment[0] += x[1] * f[2] - x[2] * f[1];
        moment[1] += x[2] * f[0] - x[0] * f[2];
        moment[2] += x[0] * f[1] - x[1] * f[0];
    }
    expectTrue("a stretched element pulls its far face back", pullBack < 0.0);
    // The magnitude has a closed form, and it is *not* E eps A. Scaling only the x
    // coordinates imposes a uniaxial **strain**, not a uniaxial stress: the
    // transverse strains are held at zero, so the modulus is the constrained one,
    // lambda + 2 mu. (This test first asserted E eps A and was 34.6% out, which is
    // exactly (1-nu)/((1+nu)(1-2nu)) -- the expectation was wrong, not the code.)
    const double nu = steel.poissonRatio, e = steel.youngsModulus;
    const double constrained = e * (1.0 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu));
    const double expected = constrained * 1.0e-4 * 0.05 * 0.01;
    expectNear("and by exactly (lambda + 2 mu) eps A", -pullBack, expected, 1e-8 * expected);

    double worstTotal = 0, worstMoment = 0;
    for (int i = 0; i < 3; ++i) {
        worstTotal = std::max(worstTotal, std::abs(total[i]));
        worstMoment = std::max(worstMoment, std::abs(moment[i]));
    }
    expectTrue("the internal force is self-equilibrated", worstTotal < 1e-9 * std::abs(pullBack));
    expectTrue("and carries no net moment",
               worstMoment < 1e-9 * std::abs(pullBack) * 0.05);
}

// --- 3. The patch test --------------------------------------------------------

struct PatchResult {
    double displacementError = 0;
    double stressError = 0;
    int freeNodes = 0;
    double distortion = 0;
};

// `warp` moves interior node layers through the thickness, which is the distortion
// the assumed natural strain interpolation is *not* exact for; `distort` moves
// nodes in plane, which it is.
PatchResult runPatch(Formulation form, int n, int nz, double distort, double warp) {
    const StructuralMaterial steel = steelMaterial();
    const double size = 1.0, thickness = 0.3;
    HexMesh mesh = makePlateMesh(size, size, thickness, n, n, nz);
    const int sy = n + 1, sz = nz + 1;

    std::mt19937 rng(20260803u);
    std::uniform_real_distribution<double> jitter(-1.0, 1.0);
    PatchResult result;
    for (int i = 0; i <= n; ++i)
        for (int j = 0; j <= n; ++j) {
            const double dx = (i == 0 || i == n) ? 0.0 : distort * jitter(rng) * size / n;
            const double dy = (j == 0 || j == n) ? 0.0 : distort * jitter(rng) * size / n;
            for (int k = 0; k <= nz; ++k) {
                const std::size_t node = static_cast<std::size_t>((i * sy + j) * sz + k);
                mesh.position[node * 3 + 0] += dx;
                mesh.position[node * 3 + 1] += dy;
                if (k != 0 && k != nz)
                    mesh.position[node * 3 + 2] += warp * thickness * std::sin(2.7 * i + 1.3 * j);
                result.distortion =
                    std::max(result.distortion, std::max(std::abs(dx), std::abs(dy)));
            }
        }

    // A general linear field: three normal strains, three shears, a rigid rotation
    // and a translation, all at once. A field with any of those missing would let a
    // whole family of errors through.
    const double gradient[3][3] = {{1.1e-4, 3.0e-5, -2.0e-5},
                                   {-4.0e-5, 0.7e-4, 1.5e-5},
                                   {2.5e-5, -1.0e-5, -0.9e-4}};
    const double offset[3] = {1.0e-3, -2.0e-3, 3.0e-3};
    const auto field = [&](const double* p, int i) {
        return gradient[i][0] * p[0] + gradient[i][1] * p[1] + gradient[i][2] * p[2] + offset[i];
    };

    for (std::size_t node = 0; node < mesh.nodeCount(); ++node) {
        const double* p = &mesh.position[node * 3];
        const std::size_t k = node % static_cast<std::size_t>(sz);
        const bool boundary = p[0] < 1e-12 || p[0] > size - 1e-12 || p[1] < 1e-12 ||
                              p[1] > size - 1e-12 || k == 0 || k == static_cast<std::size_t>(nz);
        if (!boundary) {
            ++result.freeNodes;
            continue;
        }
        for (int i = 0; i < 3; ++i) mesh.pin(node, i, field(p, i));
    }

    double want[6];
    constantStress(steel, gradient, want);
    double stressScale = 0;
    for (double s : want) stressScale = std::max(stressScale, std::abs(s));

    std::vector<double> load(mesh.nodeCount() * 3, 0.0), displacement;
    std::string problem;
    if (!solveStatic(mesh, steel, form, load, displacement, &problem)) {
        result.displacementError = result.stressError = 1.0;
        return result;
    }

    double displacementScale = 0;
    for (std::size_t node = 0; node < mesh.nodeCount(); ++node)
        for (int i = 0; i < 3; ++i) {
            const double exact = field(&mesh.position[node * 3], i);
            displacementScale = std::max(displacementScale, std::abs(exact));
            result.displacementError =
                std::max(result.displacementError, std::abs(displacement[node * 3 + i] - exact));
        }
    result.displacementError /= displacementScale;

    for (std::size_t e = 0; e < mesh.elementCount(); ++e) {
        double nodes[kDof], u[kDof], stress[kGauss * 6];
        mesh.gather(e, mesh.position, nodes);
        mesh.gather(e, displacement, u);
        elementStress(nodes, u, steel, form, stress);
        for (int gp = 0; gp < kGauss; ++gp)
            for (int i = 0; i < 6; ++i)
                result.stressError =
                    std::max(result.stressError, std::abs(stress[gp * 6 + i] - want[i]));
    }
    result.stressError /= stressScale;
    return result;
}

void testPatchTest() {
    std::printf("\n   patch test: linear field prescribed on the boundary, interior free\n");

    const PatchResult regular = runPatch(Formulation::SolidShell, 3, 2, 0.0, 0.0);
    const PatchResult distorted = runPatch(Formulation::SolidShell, 3, 2, 0.22, 0.0);
    const PatchResult control = runPatch(Formulation::Displacement, 3, 2, 0.22, 0.0);
    const PatchResult thick = runPatch(Formulation::SolidShell, 3, 3, 0.22, 0.0);

    // Guard against a vacuous patch: interior nodes have to exist, and the
    // distorted patch has to actually be distorted.
    expectTrue("the patch has free interior nodes", regular.freeNodes >= 4);
    expectTrue("the distorted patch is genuinely distorted", distorted.distortion > 0.01);
    expectTrue("the regular patch is not distorted", regular.distortion == 0.0);

    std::printf("     %-40s %12s %12s\n", "", "displacement", "stress");
    const struct {
        const char* label;
        const PatchResult& result;
    } rows[] = {{"regular grid, solid-shell", regular},
                {"in-plane distorted, solid-shell", distorted},
                {"in-plane distorted, plain hex (control)", control},
                {"in-plane distorted, 3 layers", thick}};
    for (const auto& row : rows)
        std::printf("     %-40s %12.2e %12.2e\n", row.label, row.result.displacementError,
                    row.result.stressError);

    expectTrue("regular patch reproduces the linear field", regular.displacementError < 1e-12);
    expectTrue("regular patch reproduces constant stress", regular.stressError < 1e-11);
    // The one that matters. In-plane distortion is exactly what a plate mesh on a
    // ship's shell has, and the assumed natural strain interpolation is exact for
    // it because the element is prismatic through its thickness.
    expectTrue("distorted patch reproduces the linear field", distorted.displacementError < 1e-12);
    expectTrue("distorted patch reproduces constant stress exactly",
               distorted.stressError < 1e-11);
    expectTrue("a plain hex passes the same patch (the patch itself is sound)",
               control.stressError < 1e-11);
    expectTrue("three layers through the thickness pass too", thick.stressError < 1e-11);

    // Warping the element -- moving nodes through the thickness so the mid-surface
    // is no longer flat -- is the case ANS is *not* exact for. The consistency
    // error is proportional to the warp, so on a mesh of a smooth surface it
    // vanishes as O(h): it is a bound on how coarsely a curved shell may be
    // meshed, not a defect that refinement cannot reach.
    std::printf("     warped patch (ANS is not exact for a non-prismatic element):\n");
    double previous = -1;
    bool monotone = true;
    for (double warp : {0.02, 0.05, 0.10, 0.20}) {
        const PatchResult warped = runPatch(Formulation::SolidShell, 3, 2, 0.22, warp);
        std::printf("       warp %.2f of the layer     stress error %10.3e   error/warp %6.2f\n",
                    warp, warped.stressError, warped.stressError / warp);
        if (previous >= 0 && !(warped.stressError > previous)) monotone = false;
        previous = warped.stressError;
    }
    expectTrue("warped-patch error grows with the warp", monotone);
    const PatchResult smallWarp = runPatch(Formulation::SolidShell, 3, 2, 0.22, 0.02);
    const PatchResult bigWarp = runPatch(Formulation::SolidShell, 3, 2, 0.22, 0.20);
    // Ten times the warp for roughly ten times the error: first order, not a
    // catastrophe and not a free lunch.
    expectTrue("warped-patch error is close to linear in the warp",
               bigWarp.stressError / smallWarp.stressError > 5.0 &&
                   bigWarp.stressError / smallWarp.stressError < 40.0);
    // And the plain hex passes it, which is what makes this the ANS interpolation's
    // limitation specifically rather than something about the warped mesh.
    const PatchResult warpedControl = runPatch(Formulation::Displacement, 3, 2, 0.22, 0.20);
    expectTrue("a plain hex still passes the warped patch", warpedControl.stressError < 1e-11);
    std::printf("       plain hex at warp 0.20     stress error %10.3e (the ANS interpolation is"
                " the term at fault)\n",
                warpedControl.stressError);
}

// --- 4. Bending against closed forms ------------------------------------------

// PL^3/3EI on a narrow strip, and the residual swept against b/t to establish
// that what is left over is beam theory's error and not the element's.
void testCantileverAgainstBeamTheory() {
    std::printf("\n   cantilever tip deflection against PL^3/3EI\n");
    const StructuralMaterial steel = steelMaterial();
    const double length = 1.0, thickness = 0.02, force = 100.0;

    std::printf("     %8s %10s %10s %10s\n", "b/t", "n=16", "n=32", "n=64");
    std::vector<double> plateau;
    for (double widthRatio : {0.5, 1.0, 2.0, 4.0, 8.0}) {
        const double width = widthRatio * thickness;
        const double inertia = width * thickness * thickness * thickness / 12.0;
        // Timoshenko shear correction; 0.03% here, included so the reference is the
        // better theory rather than the more convenient one.
        const double shear = 0.6 * (1.0 + steel.poissonRatio) * (thickness / length) *
                             (thickness / length);
        const double theory =
            force * length * length * length / (3.0 * steel.youngsModulus * inertia) *
            (1.0 + shear);
        std::printf("     %8.1f", widthRatio);
        double last = 0;
        for (int n : {16, 32, 64}) {
            HexMesh mesh = makePlateMesh(length, width, thickness, n, std::max(n / 8, 1), 1);
            clampFace(mesh, 0, 0.0);
            std::vector<double> load = endLoad(mesh, length, force), displacement;
            std::string problem;
            expectTrue("cantilever solves", solveStatic(mesh, steel, Formulation::SolidShell, load,
                                                        displacement, &problem));
            last = -sectionDeflection(mesh, displacement, length) / theory;
            std::printf(" %10.5f", last);
        }
        std::printf("\n");
        plateau.push_back(std::abs(1.0 - last));
    }

    // A narrow strip *is* a beam, and there the answer is beam theory's to within
    // a fifth of a percent. It gets steadily worse as the strip widens, because a
    // wide strip is a plate: the clamped root stops it curving anticlastically and
    // the rigidity climbs towards E/(1-nu^2). That is a statement about the
    // reference, so it is asserted as a trend rather than absorbed by a tolerance.
    // (Convergence here is first order, not the second order the cylindrical
    // bending case shows, because a free edge meeting a clamped root is a corner
    // singularity. That is a property of the problem, not of the element.)
    expectTrue("a narrow strip matches beam theory to 0.2%", plateau.front() < 2.0e-3);
    for (std::size_t i = 1; i < plateau.size(); ++i)
        expectTrue("the departure from beam theory grows with b/t",
                   plateau[i] > plateau[i - 1] * 0.99);
    expectTrue("and reaches about 1.4% at b/t = 8, not more",
               plateau.back() > 5.0e-3 && plateau.back() < 2.5e-2);
}

// The same cantilever with the width fully restrained, which is exact cylindrical
// bending: rigidity D = E t^3 / 12(1 - nu^2) with nothing left to argue about.
// This is also the instrument that isolates thickness locking, because a solid
// shell without the enhanced thickness strain is plane strain through its
// thickness as well, and the ratio of the two moduli is a closed form.
void testCylindricalBendingAndThicknessLocking() {
    std::printf("\n   strip in cylindrical bending, PL^3/3Db, and where the locking goes\n");
    const StructuralMaterial steel = steelMaterial();
    const double length = 1.0, width = 0.25, thickness = 0.01, force = 100.0;
    const double nu = steel.poissonRatio;
    const double theory =
        force * length * length * length / (3.0 * flexuralRigidity(steel, thickness) * width);

    // Without the enhanced thickness strain the element cannot let eps_zz vary
    // through the plate, so sigma_zz is not released and the modulus is the
    // oedometer one. Both are closed forms; their ratio is the prediction.
    const double plateModulus = steel.youngsModulus / (1.0 - nu * nu);
    const double lockedModulus =
        steel.youngsModulus * (1.0 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu));
    const double predicted = plateModulus / lockedModulus;
    std::printf("     theory %.8g m; predicted thickness-locked ratio %.4f\n", theory, predicted);

    std::vector<double> errors;
    double ansRatio = 0;
    for (int r : {1, 2, 4, 8}) {
        HexMesh mesh = makePlateMesh(length, width, thickness, 8 * r, 2 * r, 1);
        clampFace(mesh, 0, 0.0);
        restrainWidth(mesh);
        std::vector<double> load = endLoad(mesh, length, force), displacement;
        std::string problem;
        expectTrue("cylindrical bending solves",
                   solveStatic(mesh, steel, Formulation::SolidShell, load, displacement, &problem));
        const double w = -sectionDeflection(mesh, displacement, length);
        errors.push_back(std::abs(w / theory - 1.0));

        expectTrue("ANS-only solves", solveStatic(mesh, steel, Formulation::AssumedNaturalStrain,
                                                  load, displacement, &problem));
        ansRatio = -sectionDeflection(mesh, displacement, length) / theory;
        std::printf("     %4d elements: solid-shell %.6f of theory, ANS without EAS %.6f\n",
                    static_cast<int>(mesh.elementCount()), w / theory, ansRatio);
    }

    expectNear("the finest mesh reaches the closed form", 1.0 - errors.back(), 1.0, 1.0e-3);
    // Second order, measured across the halvings rather than asserted from one.
    // The last refinement is excluded from the *rate*: its error is 1e-5, which is
    // where the tip load's Saint-Venant end effect and the reference's own
    // slenderness assumption live, so the ratio there measures the noise floor and
    // not the element. Reporting an order of 4.03 from it would be reading a
    // convergence rate out of two numbers that are both essentially zero.
    int rates = 0;
    for (std::size_t i = 1; i < errors.size(); ++i) {
        const double order = std::log2(errors[i - 1] / errors[i]);
        const bool meaningful = errors[i] > 5.0e-5;
        std::printf("       refinement %zu -> %zu: error %.3e, observed order %.2f%s\n", i - 1, i,
                    errors[i], order, meaningful ? "" : "   (at the noise floor, not counted)");
        if (!meaningful) continue;
        ++rates;
        expectTrue("convergence is second order", order > 1.7 && order < 2.6);
    }
    expectTrue("the convergence rate was measured on more than one refinement", rates >= 2);
    // The whole reason the enhanced modes are carried: without them the element is
    // 22% too stiff, and the number is not a fudge, it is the ratio of two moduli.
    expectNear("thickness locking without EAS matches the closed-form penalty", ansRatio,
               predicted, 5.0e-3);
}

// --- 5. Stress recovery, and what geometric distortion costs -------------------

// The stress the element actually carries through its thickness. This is where the
// enhanced parameters have to be recovered and used: with one element through the
// plate, the *only* thing that lets sigma_zz relax to zero on the free faces is the
// enhanced thickness strain. A deflection test cannot see whether it was recovered
// with the right sign, or at all -- the condensed stiffness is the same either way.
void testThroughThicknessStress() {
    std::printf("\n   stress through a bent plate: sigma_zz must relax to zero\n");
    const StructuralMaterial steel = steelMaterial();
    const double length = 1.0, width = 0.25, thickness = 0.01, force = 100.0;
    const double nu = steel.poissonRatio;
    const double inertia = width * thickness * thickness * thickness / 12.0;
    const double q = 1.0 / std::sqrt(3.0);

    double axialError[2] = {};
    for (int refinement = 0; refinement < 2; ++refinement) {
        const int nx = 16 << refinement;
        HexMesh mesh = makePlateMesh(length, width, thickness, nx, 2, 1);
        clampFace(mesh, 0, 0.0);
        restrainWidth(mesh);
        std::vector<double> load = endLoad(mesh, length, force), displacement;
        std::string problem;
        expectTrue("bent plate solves", solveStatic(mesh, steel, Formulation::SolidShell, load,
                                                    displacement, &problem));

        double worstNormal = 0, worstAxial = 0, worstTransverse = 0, largestAxial = 0;
        for (std::size_t e = 0; e < mesh.elementCount(); ++e) {
            double nodes[kDof], u[kDof], stress[kGauss * 6];
            mesh.gather(e, mesh.position, nodes);
            double lo = nodes[0], hi = nodes[0];
            for (int a = 1; a < kNodes; ++a) {
                lo = std::min(lo, nodes[a * 3]);
                hi = std::max(hi, nodes[a * 3]);
            }
            const double centre = 0.5 * (lo + hi);
            // Stay clear of the clamped root and the loaded tip, where the real
            // solution is not beam-like and the reference would be the wrong one.
            if (centre < 0.25 * length || centre > 0.75 * length) continue;
            mesh.gather(e, displacement, u);
            elementStress(nodes, u, steel, Formulation::SolidShell, stress);

            for (int gp = 0; gp < kGauss; ++gp) {
                const double x = centre + ((gp & 1) ? q : -q) * 0.5 * (hi - lo);
                const double z = ((gp & 4) ? q : -q) * 0.5 * thickness;
                const double axial = force * (length - x) * z / inertia;
                largestAxial = std::max(largestAxial, std::abs(axial));
                worstAxial = std::max(worstAxial, std::abs(stress[gp * 6 + 0] - axial));
                // Plane strain across the width with the faces free: sigma_yy is
                // exactly nu sigma_xx, and sigma_zz is exactly zero.
                worstTransverse = std::max(worstTransverse,
                                           std::abs(stress[gp * 6 + 1] - nu * stress[gp * 6 + 0]));
                worstNormal = std::max(worstNormal, std::abs(stress[gp * 6 + 2]));
            }
        }
        expectTrue("bending stress is large enough to be a real test", largestAxial > 1.0e7);
        axialError[refinement] = worstAxial / largestAxial;
        std::printf("     nx=%2d  peak sigma_xx %.4g Pa;  |sigma_xx - Mz/I| %5.2f%%;"
                    "  |sigma_zz| %8.2e Pa (%.1e of sigma_xx)\n",
                    nx, largestAxial, 100.0 * axialError[refinement], worstNormal,
                    worstNormal / largestAxial);

        expectTrue("sigma_yy is nu sigma_xx across the restrained width",
                   worstTransverse < 1e-9 * largestAxial);
        // The assertion the enhanced parameters exist for, and it is not a
        // tolerance -- sigma_zz comes out at machine zero. Without the recovered
        // enhanced strain the element is in plane strain through its thickness and
        // sigma_zz would be nu/(1-nu) of sigma_xx, which is 43% for steel.
        expectTrue("sigma_zz relaxes to zero on the free faces",
                   worstNormal < 1e-9 * largestAxial);
    }

    // Stress converges one order below displacement, which is what a trilinear
    // element does: the recovered sigma_xx is piecewise constant along the span
    // while the real bending moment is linear, so the residual halves with the mesh
    // rather than quartering. Asserting the *rate* rather than a tolerance is what
    // distinguishes that from a formulation error, which would not move at all.
    const double order = std::log2(axialError[0] / axialError[1]);
    std::printf("     sigma_xx converges at order %.2f (first order: the stress is one order"
                " below the displacement)\n",
                order);
    expectTrue("bending stress converges at first order", order > 0.85 && order < 1.35);
    expectTrue("and is within 1.5% of M z / I on the finer mesh", axialError[1] < 0.015);
}

// In-plane bending on a mesh of trapezoids. Four of the seven enhanced parameters
// exist only for this: a non-parallelogram hex develops parasitic in-plane shear
// under membrane bending, and nothing else in this file loads the element in its
// own plane. It is the MacNeal-Harder trapezoidal cantilever, which is where
// distortion-sensitive elements are conventionally caught.
void testInPlaneBendingOnDistortedElements() {
    std::printf("\n   in-plane bending on trapezoidal elements\n");
    const StructuralMaterial steel = steelMaterial();
    const double length = 1.0, depth = 0.1, thickness = 0.02, force = 100.0;
    const double inertia = thickness * depth * depth * depth / 12.0;
    const double shear = 0.6 * (1.0 + steel.poissonRatio) * (depth / length) * (depth / length);
    const double theory =
        force * length * length * length / (3.0 * steel.youngsModulus * inertia) * (1.0 + shear);

    std::printf("     %-24s %12s %12s %12s\n", "", "regular", "parallelogram", "trapezoidal");
    double ratio[2][3] = {};
    int row = 0;
    for (Formulation form : {Formulation::SolidShell, Formulation::Displacement}) {
        int column = 0;
        for (int shape = 0; shape < 3; ++shape) {
            const int nx = 8, ny = 1;
            HexMesh mesh = makePlateMesh(length, depth, thickness, nx, ny, 1);
            const int sy = ny + 1, sz = 2;
            // Two distortions, because an enhanced-strain element treats them
            // differently: it is exact for a parallelogram (the Jacobian is
            // constant, so the Simo-Rifai scaling is too) and only approximate for
            // a trapezoid. Sliding the far edge uniformly makes parallelograms;
            // sliding it alternately makes trapezoids. The root and tip stations
            // are left alone so the clamp and the load stay where they were.
            const double skew = shape == 0 ? 0.0 : 0.35;
            for (int i = 1; i < nx; ++i)
                for (int k = 0; k < sz; ++k) {
                    const std::size_t node = static_cast<std::size_t>((i * sy + 1) * sz + k);
                    const double sign = shape == 2 ? ((i % 2 == 0) ? 1.0 : -1.0) : 1.0;
                    mesh.position[node * 3] += skew * (length / nx) * sign;
                }
            clampFace(mesh, 0, 0.0);
            std::vector<double> load(mesh.nodeCount() * 3, 0.0), displacement;
            int tip = 0;
            for (std::size_t n = 0; n < mesh.nodeCount(); ++n)
                if (std::abs(mesh.position[n * 3] - length) < 1e-12) ++tip;
            for (std::size_t n = 0; n < mesh.nodeCount(); ++n)
                if (std::abs(mesh.position[n * 3] - length) < 1e-12)
                    load[n * 3 + 1] = force / tip;
            std::string problem;
            expectTrue("in-plane cantilever solves",
                       solveStatic(mesh, steel, form, load, displacement, &problem));
            double sum = 0;
            int count = 0;
            for (std::size_t n = 0; n < mesh.nodeCount(); ++n)
                if (std::abs(mesh.position[n * 3] - length) < 1e-12) {
                    sum += displacement[n * 3 + 1];
                    ++count;
                }
            ratio[row][column++] = sum / count / theory;
        }
        std::printf("     %-24s %12.5f %12.5f %12.5f\n", name(form), ratio[row][0], ratio[row][1],
                    ratio[row][2]);
        ++row;
    }

    expectTrue("the solid-shell bends in its own plane on a regular mesh",
               std::abs(ratio[0][0] - 1.0) < 0.06);
    // Parallelogram distortion is nearly free -- the Jacobian is constant over the
    // element, so the enhanced modes are mapped exactly -- but not entirely, because
    // the element's axes no longer line up with the beam's. Measured at 5%, which is
    // recorded rather than assumed: the first version of this assertion said "almost
    // nothing" and had to be corrected by the measurement.
    expectTrue("parallelogram distortion costs the solid-shell a few percent",
               ratio[0][1] > 0.90 * ratio[0][0] && ratio[0][1] < ratio[0][0]);
    // Trapezoids are not free, for the same reason the warped patch is not: the
    // Jacobian varies and the assumed strains stop being exact.
    expectTrue("trapezoidal distortion costs it about a third, and is recorded rather"
               " than hidden",
               ratio[0][2] > 0.55 && ratio[0][2] < 0.75);
    // The enhanced in-plane modes are what stops it being far worse. On a regular
    // mesh a plain hex already loses 40% to in-plane shear locking; the four
    // Simo-Rifai modes are what recover it.
    expectTrue("the enhanced modes are worth having on a regular in-plane mesh",
               ratio[0][0] > 1.5 * ratio[1][0]);
    expectTrue("and on trapezoids", ratio[0][2] > 1.3 * ratio[1][2]);
}

// The element's headline limit, measured. The assumed natural strain interpolation
// is exact only for an element that is prismatic through its thickness; where the
// two faces are not parallel it introduces a parasitic strain, and the warped patch
// test above already showed the consistency error is first order in the
// non-parallelism. In bending that error is squared into a spurious stiffness. This
// pins the number, so that a change in the formulation cannot quietly move it.
void testTrapezoidalElementLimit() {
    std::printf("\n   the limit: elements whose faces are not parallel\n");
    const StructuralMaterial steel = steelMaterial();
    const double length = 1.0, width = 0.25, thickness = 0.01, force = 100.0;
    const double theory =
        force * length * length * length / (3.0 * flexuralRigidity(steel, thickness) * width);

    const auto solve = [&](double shift, int nx) {
        const int ny = 2, nz = 2;
        HexMesh mesh = makePlateMesh(length, width, thickness, nx, ny, nz);
        const int sy = ny + 1, sz = nz + 1;
        // Two elements through the thickness with the interior layer displaced up
        // and down alternately: every element is a trapezoid in section while the
        // strip itself is still a uniform box, so the closed form is unchanged and
        // whatever is left over is the element's.
        for (int i = 0; i <= nx; ++i)
            for (int j = 0; j <= ny; ++j) {
                const std::size_t node = static_cast<std::size_t>((i * sy + j) * sz + 1);
                mesh.position[node * 3 + 2] += shift * thickness * ((i % 2 == 0) ? 1.0 : -1.0);
            }
        clampFace(mesh, 0, 0.0);
        restrainWidth(mesh);
        std::vector<double> load = endLoad(mesh, length, force), displacement;
        std::string problem;
        if (!solveStatic(mesh, steel, Formulation::SolidShell, load, displacement, &problem))
            return 0.0;
        return -sectionDeflection(mesh, displacement, length) / theory;
    };

    const double prismatic = solve(0.0, 8);
    std::printf("     faces parallel:              %.5f of theory\n", prismatic);
    expectTrue("two elements through the thickness are still right when prismatic",
               std::abs(prismatic - 1.0) < 0.03);

    std::printf("     %14s %12s %20s\n", "face offset", "of theory", "excess / offset^2");
    double coefficients[3] = {};
    int index = 0;
    for (double shift : {0.05, 0.10, 0.20}) {
        const double ratio = solve(shift, 8);
        coefficients[index++] = (1.0 / ratio - 1.0) / (shift * shift);
        std::printf("     %13.2ft %12.5f %20.1f\n", shift, ratio, coefficients[index - 1]);
    }
    // Quadratic in the offset, which is the signature of a parasitic strain: the
    // strain is linear in the distortion and the energy is its square. A rate that
    // was not quadratic would mean something else was wrong.
    for (int i = 0; i < 3; ++i)
        expectTrue("the spurious stiffness is quadratic in the face offset",
                   coefficients[i] > 60.0 && coefficients[i] < 130.0);

    // It is a consistency error, not an inconsistency: refining while holding the
    // mesh's *shape* fixed converges. Refining while holding the offset fixed does
    // not, because that steepens the elements as it goes.
    std::printf("     refined at a fixed element shape:");
    double previous = 0;
    bool improving = true;
    for (int nx : {4, 8, 16, 32}) {
        const double ratio = solve(0.2 * 8.0 / nx, nx);
        std::printf("  nx=%2d %.4f", nx, ratio);
        if (nx > 4 && !(ratio > previous - 1e-9)) improving = false;
        previous = ratio;
    }
    std::printf("\n");
    expectTrue("refinement at a fixed element shape converges toward the answer", improving);
    expectTrue("and gets most of the way back", previous > 0.85);
}

// --- 6. Plates and strips against series solutions -----------------------------

// Navier's double series for a simply supported rectangular plate under uniform
// pressure, summed here. Quoting 0.00406 would be asserting against a table.
double navierCentreDeflection(double pressure, double side, double rigidity) {
    double sum = 0;
    for (int m = 1; m <= 299; m += 2)
        for (int n = 1; n <= 299; n += 2) {
            const double sm = ((m - 1) / 2) % 2 == 0 ? 1.0 : -1.0;
            const double sn = ((n - 1) / 2) % 2 == 0 ? 1.0 : -1.0;
            const double mm = m, nn = n;
            sum += sm * sn / (mm * nn * (mm * mm + nn * nn) * (mm * mm + nn * nn));
        }
    const double pi = std::numbers::pi;
    return 16.0 * pressure * std::pow(side, 4.0) / (std::pow(pi, 6.0) * rigidity) * sum;
}

// A quarter plate: symmetry planes at x = 0 and y = 0, support along the outer
// edges. Symmetry rather than a whole plate with a minimal constraint set, because
// u_x = 0 on both node layers of the symmetry plane is exactly dw/dx = 0 there,
// and it removes every rigid body mode without adding a single arbitrary pin.
double quarterPlateCentreDeflection(double side, double thickness, double pressure, int n,
                                    bool clamped, Formulation form) {
    const StructuralMaterial steel = steelMaterial();
    HexMesh mesh = makePlateMesh(0.5 * side, 0.5 * side, thickness, n, n, 1);
    for (std::size_t node = 0; node < mesh.nodeCount(); ++node) {
        const double x = mesh.position[node * 3], y = mesh.position[node * 3 + 1];
        if (x < 1e-12) mesh.pin(node, 0, 0.0);
        if (y < 1e-12) mesh.pin(node, 1, 0.0);
        if (!(x > 0.5 * side - 1e-12 || y > 0.5 * side - 1e-12)) continue;
        mesh.pin(node, 2, 0.0);
        if (clamped) {
            mesh.pin(node, 0, 0.0);
            mesh.pin(node, 1, 0.0);
        }
    }
    std::vector<double> load = uniformPressureLoad(mesh, pressure), displacement;
    std::string problem;
    if (!solveStatic(mesh, steel, form, load, displacement, &problem)) return 0.0;
    double sum = 0;
    int count = 0;
    for (std::size_t node = 0; node < mesh.nodeCount(); ++node)
        if (mesh.position[node * 3] < 1e-12 && mesh.position[node * 3 + 1] < 1e-12) {
            sum += displacement[node * 3 + 2];
            ++count;
        }
    return count > 0 ? -sum / count : 0.0;
}

void testPlateSeriesSolutions() {
    std::printf("\n   square plate under uniform pressure, against the series solutions\n");
    const StructuralMaterial steel = steelMaterial();
    const double side = 1.0, pressure = 1.0e5;

    // Guard: the total load must actually reach the plate. A pressure routine that
    // found no boundary face would silently solve an unloaded plate.
    {
        HexMesh mesh = makePlateMesh(0.5, 0.5, 0.01, 4, 4, 1);
        const std::vector<double> load = uniformPressureLoad(mesh, pressure);
        double total = 0;
        for (std::size_t i = 2; i < load.size(); i += 3) total += load[i];
        expectNear("uniform pressure integrates to pressure times area", total,
                   -pressure * 0.25, 1e-6 * pressure * 0.25);
    }

    const double thickness = side / 200.0;
    const double rigidity = flexuralRigidity(steel, thickness);
    const double supported = navierCentreDeflection(pressure, side, rigidity);
    std::printf("     Navier series: %.8g m, i.e. %.6f q a^4 / D (the table says 0.00406)\n",
                supported, supported * rigidity / (pressure * std::pow(side, 4.0)));
    expectNear("the summed Navier series is the tabulated coefficient",
               supported * rigidity / (pressure * std::pow(side, 4.0)), 0.00406, 5e-6);

    std::vector<double> errors;
    for (int n : {4, 8, 16}) {
        const double w =
            quarterPlateCentreDeflection(side, thickness, pressure, n, false,
                                         Formulation::SolidShell);
        errors.push_back(std::abs(w / supported - 1.0));
        std::printf("     simply supported, %2d x %2d per quarter: %.8g  (%+.3f%%)\n", n, n, w,
                    (w / supported - 1.0) * 100.0);
    }
    expectTrue("the plate converges toward the series solution",
               errors[1] < errors[0] && errors[2] < errors[1]);
    expectTrue("and reaches it within 0.2%", errors.back() < 2.0e-3);

    // The residual is Kirchhoff's, not the element's: plate theory drops transverse
    // shear, whose contribution is O((t/a)^2). Sweeping the thickness at a fixed
    // mesh is the instrument that says so -- if the residual belonged to the
    // element it would not care how thin the plate is.
    std::printf("     residual against Kirchhoff at a fixed 16 x 16 mesh, vs slenderness:\n");
    double previous = 1e30;
    for (double slenderness : {25.0, 50.0, 100.0, 200.0}) {
        const double t = side / slenderness;
        const double d = flexuralRigidity(steel, t);
        const double want = navierCentreDeflection(pressure, side, d);
        const double w = quarterPlateCentreDeflection(side, t, pressure, 16, false,
                                                      Formulation::SolidShell);
        const double residual = std::abs(w / want - 1.0);
        std::printf("       a/t = %5.0f   %+.3f%%\n", slenderness, (w / want - 1.0) * 100.0);
        expectTrue("the Kirchhoff residual falls as the plate thins", residual < previous);
        previous = residual;
    }

    // Clamped. There is no elementary closed form, so this is against Timoshenko's
    // tabulated 0.00126 -- a coefficient quoted to three figures, which is the
    // tolerance it is worth.
    const double clampedCoefficient =
        quarterPlateCentreDeflection(side, thickness, pressure, 16, true,
                                     Formulation::SolidShell) *
        rigidity / (pressure * std::pow(side, 4.0));
    std::printf("     clamped, 16 x 16 per quarter: coefficient %.6f (Timoshenko 0.00126)\n",
                clampedCoefficient);
    expectNear("the clamped plate matches the tabulated coefficient", clampedCoefficient, 0.00126,
               1.0e-5);
}

// Strips under uniform pressure, where both support conditions have elementary
// closed forms and the clamped one is the only check in this file on how the
// element behaves against a fully built-in edge.
void testStripsUnderPressure() {
    std::printf("\n   strips under uniform pressure: 5qL^4/384D and qL^4/384D\n");
    const StructuralMaterial steel = steelMaterial();
    const double span = 1.0, width = 0.1, thickness = 0.005, pressure = 2.0e4;
    const double rigidity = flexuralRigidity(steel, thickness);

    for (int mode = 0; mode < 2; ++mode) {
        const double theory = (mode == 0 ? 5.0 : 1.0) * pressure * std::pow(span, 4.0) /
                              (384.0 * rigidity);
        double last = 0;
        std::printf("     %-18s theory %.8g m:", mode == 0 ? "simply supported" : "clamped",
                    theory);
        for (int n : {4, 8, 16, 32}) {
            // Half the span, with x = 0 the midspan symmetry plane.
            HexMesh mesh = makePlateMesh(0.5 * span, width, thickness, n, 1, 1);
            restrainWidth(mesh);
            clampFace(mesh, 0, 0.0, /*onlyNormal=*/true);  // symmetry at midspan
            if (mode == 0) {
                clampFace(mesh, 0, 0.5 * span, /*onlyNormal=*/false);
                // A simple support restrains w only; undo the two the clamp added.
                for (std::size_t node = 0; node < mesh.nodeCount(); ++node)
                    if (std::abs(mesh.position[node * 3] - 0.5 * span) < 1e-12)
                        mesh.fixed[node * 3] = 0u;
            } else {
                clampFace(mesh, 0, 0.5 * span, /*onlyNormal=*/false);
            }
            std::vector<double> load = uniformPressureLoad(mesh, pressure), displacement;
            std::string problem;
            expectTrue("strip solves",
                       solveStatic(mesh, steel, Formulation::SolidShell, load, displacement,
                                   &problem));
            last = -sectionDeflection(mesh, displacement, 0.0) / theory;
            std::printf("  n=%2d %.5f", n, last);
        }
        std::printf("\n");
        expectNear(mode == 0 ? "simply supported strip is 5qL^4/384D"
                             : "clamped strip is qL^4/384D",
                   last, 1.0, 3.0e-3);
    }
}

// --- 7. Locking, demonstrated -------------------------------------------------

// A linear tetrahedron's stiffness, built here rather than borrowed from fem.cpp:
// that path is explicit, float, and lumped-mass, and running the two elements down
// two different solver paths would measure the paths. The connectivity *is*
// fem.cpp's, so the element being compared against is the one the spike rejected.
void linearTetStiffness(const double p[12], const StructuralMaterial& material, double out[144]) {
    double m[3][3];
    for (int c = 0; c < 3; ++c)
        for (int i = 0; i < 3; ++i) m[i][c] = p[(c + 1) * 3 + i] - p[i];
    const double det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
                       m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
                       m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    const double volume = std::abs(det) / 6.0;
    const double id = 1.0 / det;
    double gradient[4][3];
    gradient[1][0] = (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * id;
    gradient[1][1] = (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * id;
    gradient[1][2] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * id;
    gradient[2][0] = (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * id;
    gradient[2][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * id;
    gradient[2][2] = (m[0][2] * m[1][0] - m[0][0] * m[1][2]) * id;
    gradient[3][0] = (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * id;
    gradient[3][1] = (m[0][1] * m[2][0] - m[0][0] * m[2][1]) * id;
    gradient[3][2] = (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * id;
    for (int k = 0; k < 3; ++k)
        gradient[0][k] = -(gradient[1][k] + gradient[2][k] + gradient[3][k]);

    double b[6][12] = {};
    for (int a = 0; a < 4; ++a) {
        b[0][a * 3 + 0] = gradient[a][0];
        b[1][a * 3 + 1] = gradient[a][1];
        b[2][a * 3 + 2] = gradient[a][2];
        b[3][a * 3 + 0] = gradient[a][1];
        b[3][a * 3 + 1] = gradient[a][0];
        b[4][a * 3 + 1] = gradient[a][2];
        b[4][a * 3 + 2] = gradient[a][1];
        b[5][a * 3 + 2] = gradient[a][0];
        b[5][a * 3 + 0] = gradient[a][2];
    }
    const double e = material.youngsModulus, nu = material.poissonRatio;
    const double lambda = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu)), mu = e / (2.0 * (1.0 + nu));
    double c[6][6] = {};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) c[i][j] = lambda;
        c[i][i] += 2.0 * mu;
        c[3 + i][3 + i] = mu;
    }
    for (int i = 0; i < 12; ++i)
        for (int j = 0; j < 12; ++j) {
            double s = 0;
            for (int r = 0; r < 6; ++r)
                for (int k = 0; k < 6; ++k) s += b[r][i] * c[r][k] * b[k][j];
            out[i * 12 + j] = volume * s;
        }
}

// The same strip, the same load, the same solver, made of linear tets.
double tetStripDeflection(double length, double width, double thickness, double force, int nx,
                          int ny, int nz) {
    const StructuralMaterial steel = steelMaterial();
    fem::TetMesh tets =
        fem::makeBoxTetMesh(static_cast<float>(length), static_cast<float>(width),
                            static_cast<float>(thickness), nx, ny, nz);
    const std::size_t nodes = tets.position.size() / 3;
    std::vector<double> position(nodes * 3);
    for (std::size_t i = 0; i < nodes * 3; ++i) position[i] = tets.position[i];

    std::vector<std::ptrdiff_t> map(nodes * 3, -1);
    std::size_t free = 0;
    for (std::size_t n = 0; n < nodes; ++n)
        for (int i = 0; i < 3; ++i) {
            // Root clamped, and no motion across the width anywhere: the same
            // cylindrical bending condition the hex strip is under.
            if (position[n * 3] < 1e-9 || i == 1) continue;
            map[n * 3 + i] = static_cast<std::ptrdiff_t>(free++);
        }
    const std::size_t tetCount = tets.index.size() / 4;
    std::size_t band = 0;
    for (std::size_t e = 0; e < tetCount; ++e) {
        std::size_t lo = free, hi = 0;
        bool any = false;
        for (int a = 0; a < 4; ++a)
            for (int i = 0; i < 3; ++i) {
                const std::ptrdiff_t d = map[tets.index[e * 4 + static_cast<std::size_t>(a)] * 3 +
                                             static_cast<std::size_t>(i)];
                if (d < 0) continue;
                lo = std::min(lo, static_cast<std::size_t>(d));
                hi = std::max(hi, static_cast<std::size_t>(d));
                any = true;
            }
        if (any) band = std::max(band, hi - lo);
    }

    BandedSpd system(free, band);
    std::vector<double> rhs(free, 0.0);
    int tip = 0;
    for (std::size_t n = 0; n < nodes; ++n)
        if (position[n * 3] > length - 1e-9) ++tip;
    for (std::size_t n = 0; n < nodes; ++n)
        if (position[n * 3] > length - 1e-9 && map[n * 3 + 2] >= 0)
            rhs[static_cast<std::size_t>(map[n * 3 + 2])] -= force / tip;

    for (std::size_t e = 0; e < tetCount; ++e) {
        double p[12];
        for (int a = 0; a < 4; ++a)
            for (int i = 0; i < 3; ++i)
                p[a * 3 + i] = position[tets.index[e * 4 + static_cast<std::size_t>(a)] * 3 +
                                        static_cast<std::size_t>(i)];
        double ke[144];
        linearTetStiffness(p, steel, ke);
        for (int pi = 0; pi < 12; ++pi) {
            const std::ptrdiff_t r =
                map[tets.index[e * 4 + static_cast<std::size_t>(pi / 3)] * 3 +
                    static_cast<std::size_t>(pi % 3)];
            if (r < 0) continue;
            for (int qi = 0; qi < 12; ++qi) {
                const std::ptrdiff_t c =
                    map[tets.index[e * 4 + static_cast<std::size_t>(qi / 3)] * 3 +
                        static_cast<std::size_t>(qi % 3)];
                if (c < 0) continue;
                system.add(static_cast<std::size_t>(r), static_cast<std::size_t>(c),
                           ke[pi * 12 + qi]);
            }
        }
    }
    if (!system.factor()) return 0.0;
    system.solve(rhs);
    double sum = 0;
    int count = 0;
    for (std::size_t n = 0; n < nodes; ++n)
        if (position[n * 3] > length - 1e-9) {
            sum += rhs[static_cast<std::size_t>(map[n * 3 + 2])];
            ++count;
        }
    return count > 0 ? -sum / count : 0.0;
}

void testLockingDemonstration() {
    std::printf("\n   locking, demonstrated rather than assumed: the same 8 x 2 x 1 mesh\n");
    std::printf("     %8s %13s %13s %13s %13s\n", "L/t", "solid-shell", "ANS, no EAS",
                "plain hex", "linear tets");
    const StructuralMaterial steel = steelMaterial();
    const double length = 1.0, width = 0.25, force = 100.0;

    double worstSolidShell = 0, worstHex = 0, worstTet = 0;
    for (double slenderness : {10.0, 20.0, 50.0, 100.0, 200.0, 500.0}) {
        const double thickness = length / slenderness;
        const double theory = force * length * length * length /
                              (3.0 * flexuralRigidity(steel, thickness) * width);
        double ratio[3] = {};
        int column = 0;
        for (Formulation form : {Formulation::SolidShell, Formulation::AssumedNaturalStrain,
                                 Formulation::Displacement}) {
            HexMesh mesh = makePlateMesh(length, width, thickness, 8, 2, 1);
            clampFace(mesh, 0, 0.0);
            restrainWidth(mesh);
            std::vector<double> load = endLoad(mesh, length, force), displacement;
            std::string problem;
            expectTrue("locking sweep solves",
                       solveStatic(mesh, steel, form, load, displacement, &problem));
            ratio[column++] = -sectionDeflection(mesh, displacement, length) / theory;
        }
        const double tet = tetStripDeflection(length, width, thickness, force, 8, 2, 1) / theory;
        std::printf("     %8.0f %13.5f %13.5f %13.5f %13.5f\n", slenderness, ratio[0], ratio[1],
                    ratio[2], tet);
        worstSolidShell = std::max(worstSolidShell, std::abs(ratio[0] - 1.0));
        if (slenderness >= 100.0) {
            worstHex = std::max(worstHex, ratio[2]);
            worstTet = std::max(worstTet, tet);
        }
    }

    // The claim in docs/07 §4, as an assertion: the solid-shell is insensitive to
    // slenderness over a range where the alternatives collapse by three orders.
    expectTrue("the solid-shell holds the closed form across L/t = 10 to 500",
               worstSolidShell < 0.05);
    expectTrue("a plain hex on the same mesh locks below 2% of the answer", worstHex < 0.02);
    expectTrue("linear tets on the same mesh lock at least as hard", worstTet < worstHex);
    // Guard against a vacuous comparison: the alternatives must have been given a
    // fair chance, i.e. they are not simply broken at every slenderness.
    const double thick = tetStripDeflection(length, width, length / 4.0, force, 8, 2, 1);
    expectTrue("the tet comparison is not vacuous -- tets are fine on a stubby beam",
               thick > 0.0);
}

// --- 8. The explicit stability limit ------------------------------------------

// docs/07 §4 claimed that a solid-shell's timestep is "set by the in-plane element
// size rather than the plate thickness", worth 5-10x. That is true of a *shell*
// element with no thickness stretch. This element has one -- deliberately, because
// a crush zone needs it -- and its through-thickness dilatational mode is the
// highest frequency in the element however large it is in plane. The measurement
// is here, and the doc is corrected.
void testExplicitStabilityLimit() {
    std::printf("\n   explicit stability limit: what actually sets the timestep\n");
    const StructuralMaterial steel = steelMaterial();
    const double nu = steel.poissonRatio;
    const double constrained =
        steel.youngsModulus * (1.0 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu));
    const double waveSpeed = std::sqrt(constrained / steel.density);
    std::printf("     constrained (oedometer) wave speed %.0f m/s\n", waveSpeed);
    std::printf("     %10s %10s %14s %12s %12s\n", "thickness", "in-plane", "dt (s)",
                "dt c_p / t", "dt c_p / h");

    for (double thickness : {0.02, 0.01}) {
        double previous = 0;
        for (double inPlane : {0.02, 0.05, 0.10, 0.20, 0.50}) {
            double nodes[kDof];
            const double corner[4][2] = {{0, 0}, {inPlane, 0}, {inPlane, inPlane}, {0, inPlane}};
            for (int a = 0; a < kNodes; ++a) {
                nodes[a * 3 + 0] = corner[a % 4][0];
                nodes[a * 3 + 1] = corner[a % 4][1];
                nodes[a * 3 + 2] = (a < 4 ? -0.5 : 0.5) * thickness;
            }
            const double dt = criticalTimestep(nodes, steel, Formulation::SolidShell, 1.0);
            std::printf("     %10.3f %10.3f %14.6g %12.4f %12.4f\n", thickness, inPlane, dt,
                        dt * waveSpeed / thickness, dt * waveSpeed / inPlane);
            if (inPlane >= 5.0 * thickness) {
                // The closed form for a two-node bar of length t with lumped mass
                // is dt = t/c exactly; the plate element reaches it because the
                // thickness stretch mode is that bar.
                expectNear("dt is the thickness crossing time t / c_p", dt * waveSpeed / thickness,
                           1.0, 0.01);
                expectTrue("and stops depending on the in-plane size",
                           previous == 0.0 || std::abs(dt / previous - 1.0) < 0.01);
                previous = dt;
            }
        }
    }

    // Against the tet mesh that would be needed for the same accuracy: the spike
    // measured 11% error at eight linear tets through the thickness, so that is the
    // comparison. The win is real; it is just not the win the doc described.
    const double thickness = 0.02;
    const double tetSize = thickness / 8.0;
    const double tetStep = tetSize / std::sqrt(steel.youngsModulus / steel.density);
    double plateNodes[kDof];
    const double corner[4][2] = {{0, 0}, {0.05, 0}, {0.05, 0.05}, {0, 0.05}};
    for (int a = 0; a < kNodes; ++a) {
        plateNodes[a * 3 + 0] = corner[a % 4][0];
        plateNodes[a * 3 + 1] = corner[a % 4][1];
        plateNodes[a * 3 + 2] = (a < 4 ? -0.5 : 0.5) * thickness;
    }
    const double shellStep = criticalTimestep(plateNodes, steel, Formulation::SolidShell, 1.0);
    std::printf("     20 mm plate: solid-shell dt %.3g s; linear tets at 8 through thickness"
                " %.3g s -> %.1fx\n",
                shellStep, tetStep, shellStep / tetStep);
    expectTrue("the timestep win over a tet mesh that resolves bending is real but modest",
               shellStep / tetStep > 3.0 && shellStep / tetStep < 10.0);
}

// --- 9. What it costs ---------------------------------------------------------

void testCostPerElement() {
    std::printf("\n   measured cost, one core\n");
    const StructuralMaterial steel = steelMaterial();
    double nodes[kDof];
    const double corner[4][2] = {{0, 0}, {0.05, 0}, {0.051, 0.049}, {0.001, 0.05}};
    for (int a = 0; a < kNodes; ++a) {
        nodes[a * 3 + 0] = corner[a % 4][0];
        nodes[a * 3 + 1] = corner[a % 4][1];
        nodes[a * 3 + 2] = (a < 4 ? -0.01 : 0.01);
    }

    const auto time = [](const char* label, int repeats, auto&& body) {
        const auto begin = std::chrono::steady_clock::now();
        body(repeats);
        const auto end = std::chrono::steady_clock::now();
        const double nanoseconds =
            std::chrono::duration<double, std::nano>(end - begin).count() / repeats;
        std::printf("     %-46s %8.0f ns/element\n", label, nanoseconds);
        return nanoseconds;
    };

    double sink = 0;
    double stiffness[kDof * kDof];
    const double formation = time("solid-shell stiffness (once, at build)", 20000, [&](int n) {
        for (int i = 0; i < n; ++i) {
            elementStiffness(nodes, steel, Formulation::SolidShell, stiffness);
            sink += stiffness[0];
        }
    });
    double plain[kDof * kDof];
    const double plainFormation = time("plain hex stiffness, for comparison", 20000, [&](int n) {
        for (int i = 0; i < n; ++i) {
            elementStiffness(nodes, steel, Formulation::Displacement, plain);
            sink += plain[0];
        }
    });

    double current[kDof];
    for (int i = 0; i < kDof; ++i) current[i] = nodes[i] * 1.0001 + 1e-5;
    double force[kDof];
    const double perStep = time("solid-shell internal force (every step)", 200000, [&](int n) {
        for (int i = 0; i < n; ++i) {
            internalForce(stiffness, nodes, current, force);
            sink += force[0];
        }
    });

    // The tet the spike measured, through its own code path, so the comparison is
    // against the thing that exists rather than against a re-implementation.
    fem::Material femSteel;
    fem::TetMesh tets = fem::makeBoxTetMesh(0.4f, 0.1f, 0.02f, 16, 4, 2);
    tets.computeRestState(femSteel);
    const double tetStep =
        time("fem.cpp linear tet internal force (every step)", 200000, [&](int n) {
            float out[12];
            for (int i = 0; i < n; ++i) {
                fem::tetForces(tets, static_cast<std::size_t>(i) % tets.tetCount(),
                               femSteel.lameLambda(), femSteel.lameMu(), out);
                sink += out[0];
            }
        });

    expectTrue("the cost measurement did work that could not be optimised away",
               std::isfinite(sink));
    // Deliberately loose: a timing assertion tight enough to be interesting is a
    // flaky test on a shared machine. Three orders of headroom, so it only fires if
    // something has gone structurally wrong.
    expectTrue("stiffness formation is under a millisecond per element", formation < 1.0e6);
    expectTrue("an internal force evaluation is under 10 microseconds", perStep < 1.0e4);

    std::printf("     EAS and the assumed strains cost %.2fx a plain hex to form\n",
                formation / plainFormation);
    std::printf("     one solid-shell force step = %.1f tet force steps\n", perStep / tetStep);
    // What that buys, per square metre of 20 mm plate per second of simulated time.
    // The tet mesh has to resolve bending through the thickness, so it pays in
    // element count *and* in timestep; the solid-shell pays in neither.
    const double thickness = 0.02, inPlane = 0.05;
    const double tetEdge = thickness / 8.0;
    const double tetsPerSquareMetre = 6.0 / (tetEdge * tetEdge) * (thickness / tetEdge);
    const double shellsPerSquareMetre = 1.0 / (inPlane * inPlane);
    const double waveSpeed = std::sqrt(steel.youngsModulus / steel.density);
    const double tetSteps = 1.0 / (0.5 * tetEdge / waveSpeed);
    double plateNodes[kDof];
    for (int a = 0; a < kNodes; ++a) {
        plateNodes[a * 3 + 0] = (a % 4 == 1 || a % 4 == 2) ? inPlane : 0.0;
        plateNodes[a * 3 + 1] = (a % 4 == 2 || a % 4 == 3) ? inPlane : 0.0;
        plateNodes[a * 3 + 2] = (a < 4 ? -0.5 : 0.5) * thickness;
    }
    const double shellSteps =
        1.0 / criticalTimestep(plateNodes, steel, Formulation::SolidShell, 0.9);
    const double tetCost = tetsPerSquareMetre * tetSteps * tetStep * 1e-9;
    const double shellCost = shellsPerSquareMetre * shellSteps * perStep * 1e-9;
    std::printf("     one m^2 of 20 mm plate, one simulated second, one core:\n");
    std::printf("       linear tets at %.1f mm  %10.0f s     solid-shell at %.0f mm  %10.2f s"
                "   (%.0fx)\n",
                tetEdge * 1000.0, tetCost, inPlane * 1000.0, shellCost, tetCost / shellCost);
    expectTrue("the solid-shell is orders of magnitude cheaper per square metre of plate",
               tetCost / shellCost > 1.0e3);
}


// --- The cached rest forms ------------------------------------------------------
//
// `RestForms` is a pure optimisation: it holds the part of an element that depends
// only on its rest configuration, so that an explicit solver stops rebuilding it
// every step. The whole claim is that it is **the same arithmetic on the same
// numbers**, so the assertion is bit equality and not a tolerance. A tolerance
// here would pass a cache that had transposed the rest Jacobian, or that had been
// filled from the *deformed* geometry, which are the two ways this can be wrong.
//
// Three guards against a vacuous pass, each because the obvious version of this
// test proves nothing:
//
//   * the element is **distorted, rotated and past yield**, so every term the
//     cache carries is actually loaded. On an undeformed axis-aligned cube the
//     rest Jacobian is a multiple of the identity, where a transposed cache is
//     the same cache;
//   * the rotation is asserted to be **genuinely a rotation and not the identity**,
//     because a cached `restJacobianInverse` that came back as zeros would give
//     F = 0, and a polar decomposition of zero returns without touching its input;
//   * the plastic state is asserted to have **actually yielded and advanced**, so
//     the history comparison has something in it to compare.
void testCachedRestFormsAreBitIdentical() {
    const plasticity::Material material = plasticity::shipSteel();
    const StructuralMaterial steel = steelMaterial();

    // A plate element of the proportions a zone meshes -- 0.6 m x 0.175 m of 12 mm
    // plating -- distorted in all three directions so the Jacobian is not diagonal
    // and the enhanced modes are not decoupled from it.
    double rest[kDof];
    const double cx[kNodes] = {0, 1, 1, 0, 0, 1, 1, 0};
    const double cy[kNodes] = {0, 0, 1, 1, 0, 0, 1, 1};
    const double cz[kNodes] = {0, 0, 0, 0, 1, 1, 1, 1};
    for (int a = 0; a < kNodes; ++a) {
        rest[a * 3 + 0] = cx[a] * 0.60 + 0.07 * cy[a];   // in-plane trapezoid
        rest[a * 3 + 1] = cy[a] * 0.175;
        rest[a * 3 + 2] = cz[a] * 0.012 + 0.004 * cx[a];  // and a slope through it
    }

    // Deform it: a finite rotation on top of a bending-and-stretching field, so the
    // co-rotational frame is doing real work rather than returning the identity.
    const double angle = 0.7;
    const double c = std::cos(angle), sn = std::sin(angle);
    double current[kDof];
    for (int a = 0; a < kNodes; ++a) {
        const double x = rest[a * 3 + 0], y = rest[a * 3 + 1], z = rest[a * 3 + 2];
        const double dx = x * 1.02 + 0.9 * x * x;
        const double dy = y * 0.995;
        const double dz = z + 0.05 * x;
        current[a * 3 + 0] = c * dx - sn * dz;
        current[a * 3 + 1] = dy;
        current[a * 3 + 2] = sn * dx + c * dz;
    }

    RestForms forms;
    expectTrue("the rest forms build on a distorted element",
               computeRestForms(rest, Formulation::SolidShell, forms));
    expectTrue("and carry all seven enhanced modes", forms.easCount == kEas);

    // 1. The rotation.
    double direct[9], cached[9];
    elementRotation(rest, current, direct);
    elementRotation(forms, current, cached);
    bool sameRotation = true;
    for (int i = 0; i < 9; ++i) sameRotation = sameRotation && direct[i] == cached[i];
    expectTrue("the cached rotation is bit-identical to the recomputed one", sameRotation);
    // Guard: a cache full of zeros gives F = 0, whose polar factor comes back as
    // whatever `polarRotation` was handed. The rotation has to be a real one.
    double offDiagonal = 0, orthogonality = 0;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            double s = 0;
            for (int k = 0; k < 3; ++k) s += cached[k * 3 + i] * cached[k * 3 + j];
            orthogonality = std::max(orthogonality, std::abs(s - (i == j ? 1.0 : 0.0)));
            if (i != j) offDiagonal = std::max(offDiagonal, std::abs(cached[j * 3 + i]));
        }
    expectTrue("and it is orthogonal", orthogonality < 1e-12);
    expectTrue("and it is not the identity, so the comparison had a rotation in it",
               offDiagonal > 0.5);

    // 2. The elastic internal force.
    double stiffness[kDof * kDof];
    elementStiffness(rest, steel, Formulation::SolidShell, stiffness);
    double forceDirect[kDof], forceCached[kDof];
    internalForce(stiffness, rest, current, forceDirect);
    internalForce(forms, stiffness, rest, current, forceCached);
    bool sameForce = true;
    double biggest = 0;
    for (int i = 0; i < kDof; ++i) {
        sameForce = sameForce && forceDirect[i] == forceCached[i];
        biggest = std::max(biggest, std::abs(forceDirect[i]));
    }
    expectTrue("the cached elastic internal force is bit-identical", sameForce);
    expectTrue("and the element was carrying a force worth comparing", biggest > 1.0e5);

    // 3. The elastoplastic update, including the history it commits. Stepped
    // several times, because a stale cache is indistinguishable from a fresh one on
    // the first call -- the failure mode is a cache that stops matching once the
    // state has moved.
    ElementPlasticState stateDirect, stateCached;
    initialisePlasticState(rest, material, stateDirect);
    initialisePlasticState(rest, material, stateCached);
    bool sameUpdate = true, sameHistory = true;
    double peakPlasticStrain = 0;
    int totalYielded = 0;
    for (int step = 1; step <= 6; ++step) {
        double moving[kDof];
        const double fraction = static_cast<double>(step) / 6.0;
        for (int i = 0; i < kDof; ++i) moving[i] = rest[i] + fraction * (current[i] - rest[i]);

        double fA[kDof], fB[kDof], sA[kGauss * 6], sB[kGauss * 6];
        const PlasticUpdate a =
            elementPlasticUpdate(rest, moving, material, Formulation::SolidShell, stateDirect, fA, sA);
        const PlasticUpdate b =
            elementPlasticUpdate(forms, rest, moving, material, stateCached, fB, sB);
        for (int i = 0; i < kDof; ++i) sameUpdate = sameUpdate && fA[i] == fB[i];
        for (int i = 0; i < kGauss * 6; ++i) sameUpdate = sameUpdate && sA[i] == sB[i];
        sameUpdate = sameUpdate && a.dissipation == b.dissipation &&
                     a.iterations == b.iterations && a.yieldedPoints == b.yieldedPoints;
        for (int k = 0; k < kEas; ++k)
            sameHistory = sameHistory && stateDirect.enhanced[k] == stateCached.enhanced[k];
        for (int gp = 0; gp < kGauss; ++gp) {
            sameHistory = sameHistory && stateDirect.point[gp].equivalentPlasticStrain ==
                                             stateCached.point[gp].equivalentPlasticStrain;
            for (int i = 0; i < plasticity::kVoigt; ++i)
                sameHistory = sameHistory &&
                              stateDirect.point[gp].backStress[i] == stateCached.point[gp].backStress[i];
            peakPlasticStrain =
                std::max(peakPlasticStrain, stateCached.point[gp].equivalentPlasticStrain);
        }
        totalYielded += a.yieldedPoints;
    }
    expectTrue("the cached elastoplastic update is bit-identical, six steps in", sameUpdate);
    expectTrue("and so is the history it commits", sameHistory);
    // Guard: two elastic elements agree trivially. The comparison is only worth
    // anything if the return map, the enhanced-strain Newton and the history were
    // all exercised.
    std::printf("     cached rest forms: %d yielded point-steps, peak eps_p %.4f, rotation"
                " %.2f rad\n", totalYielded, peakPlasticStrain, angle);
    expectTrue("and the element genuinely yielded, so the history had something in it",
               totalYielded > 0 && peakPlasticStrain > 1.0e-3);

    // 4. An inverted element is refused, the same condition `elementStiffness`
    // declines to work on -- so a caller cannot cache nonsense and then solve on it.
    double folded[kDof];
    for (int i = 0; i < kDof; ++i) folded[i] = rest[i];
    std::swap(folded[0], folded[3]);   // node 0 and node 1 exchange x
    std::swap(folded[1], folded[4]);
    RestForms bad;
    expectTrue("an inverted element has no rest forms",
               !computeRestForms(folded, Formulation::SolidShell, bad) && !bad.ok);
    // And the rotation off refused forms is the **identity**, not whatever a polar
    // decomposition of a zero matrix happens to leave behind. Mutation testing put
    // this here: dropping the guard passed everything, and what it returns then is a
    // matrix of zeros -- which an explicit solver would happily use as a rotation.
    double refusedRotation[9];
    elementRotation(bad, current, refusedRotation);
    bool identity = true;
    for (int i = 0; i < 9; ++i)
        identity = identity && refusedRotation[i] == (i % 4 == 0 ? 1.0 : 0.0);
    expectTrue("a rotation off refused forms is exactly the identity", identity);
    // Guard: the good forms must give something that is *not* the identity, or the
    // check above is satisfied by a function that ignores its argument.
    double goodRotation[9];
    elementRotation(forms, current, goodRotation);
    bool differs = false;
    for (int i = 0; i < 9; ++i)
        differs = differs || std::abs(goodRotation[i] - (i % 4 == 0 ? 1.0 : 0.0)) > 0.1;
    expectTrue("and a rotation off good forms is not", differs);

    // And an update off refused forms is a no-op rather than nonsense.
    ElementPlasticState after;
    initialisePlasticState(rest, material, after);
    double refusedForce[kDof];
    const PlasticUpdate refused =
        elementPlasticUpdate(bad, rest, current, material, after, refusedForce);
    bool zeroed = !refused.converged;
    for (int i = 0; i < kDof; ++i) zeroed = zeroed && refusedForce[i] == 0.0;
    expectTrue("and an update off them returns no force rather than a wrong one", zeroed);
}


// The enhanced modes are a basis, and `computeRestForms` now normalises them so
// that Kaa is O(1) instead of O((h/t)^4). That is only allowed to be free if the
// element it produces is *identical*, so the identity is asserted rather than the
// conditioning alone.
//
// The check re-derives the condensation independently from the public `RestForms`
// and a constitutive matrix built here, then does it again with the enhanced basis
// rescaled by wildly different per-column factors. Scaling column j of G by s_j
// takes Kua -> Kua S and Kaa -> S Kaa S, so the condensation
// Kua S (S Kaa S)^-1 S Kua^T collapses back to Kua Kaa^-1 Kua^T exactly -- the
// scaling cancels, which is the whole reason the normalisation costs nothing.
void testEnhancedBasisScalingCannotChangeTheElement() {
    const StructuralMaterial steel = ah36Steel();
    // Slender, because that is where the conditioning is bad and where a mistake
    // in the cancellation would show.
    const double t = 0.010, h = 0.600;
    const double rest[24] = {
        0, 0, -t / 2,  h, 0, -t / 2,  h, h, -t / 2,  0, h, -t / 2,
        0, 0,  t / 2,  h, 0,  t / 2,  h, h,  t / 2,  0, h,  t / 2};

    solidshell::RestForms f;
    expectTrue("the slender element builds its rest forms",
               solidshell::computeRestForms(rest, solidshell::Formulation::SolidShell, f));
    expectTrue("and it carries enhanced modes", f.easCount > 0);

    // Isotropic C, written here so the comparison is against an independent
    // derivation rather than against the same code twice.
    double c[6][6] = {};
    {
        const double E = steel.youngsModulus, nu = steel.poissonRatio;
        const double a = E * (1.0 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu));
        const double b = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
        const double g = 0.5 * E / (1.0 + nu);
        c[0][0] = c[1][1] = c[2][2] = a;
        c[0][1] = c[0][2] = c[1][0] = c[1][2] = c[2][0] = c[2][1] = b;
        c[3][3] = c[4][4] = c[5][5] = g;
    }

    const int n = f.easCount;
    // Condense from the rest forms, with column j of G scaled by scale[j].
    const auto condense = [&](const std::vector<double>& columnScale, std::vector<double>& out,
                              double* conditionProxy) {
        std::vector<double> kuu(24 * 24, 0.0), kua(24 * static_cast<std::size_t>(n), 0.0),
            kaa(static_cast<std::size_t>(n) * n, 0.0);
        for (int gp = 0; gp < 8; ++gp) {
            double cb[6][24] = {}, cg[6][7] = {};
            for (int i = 0; i < 6; ++i) {
                for (int d = 0; d < 24; ++d)
                    for (int k = 0; k < 6; ++k) cb[i][d] += c[i][k] * f.b[gp][k][d];
                for (int j = 0; j < n; ++j)
                    for (int k = 0; k < 6; ++k)
                        cg[i][j] += c[i][k] * f.g[gp][k][j] * columnScale[static_cast<std::size_t>(j)];
            }
            for (int i = 0; i < 24; ++i) {
                for (int j = 0; j < 24; ++j) {
                    double s = 0;
                    for (int k = 0; k < 6; ++k) s += f.b[gp][k][i] * cb[k][j];
                    kuu[static_cast<std::size_t>(i) * 24 + j] += f.weight[gp] * s;
                }
                for (int j = 0; j < n; ++j) {
                    double s = 0;
                    for (int k = 0; k < 6; ++k) s += f.b[gp][k][i] * cg[k][j];
                    kua[static_cast<std::size_t>(i) * n + j] += f.weight[gp] * s;
                }
            }
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j) {
                    double s = 0;
                    for (int k = 0; k < 6; ++k)
                        s += f.g[gp][k][i] * columnScale[static_cast<std::size_t>(i)] * cg[k][j];
                    kaa[static_cast<std::size_t>(i) * n + j] += f.weight[gp] * s;
                }
        }
        if (conditionProxy) {
            double lo = kaa[0], hi = kaa[0];
            for (int i = 0; i < n; ++i) {
                lo = std::min(lo, kaa[static_cast<std::size_t>(i) * n + i]);
                hi = std::max(hi, kaa[static_cast<std::size_t>(i) * n + i]);
            }
            *conditionProxy = hi / lo;
        }
        // Kaa^-1 Kua^T by Gaussian elimination on the small system.
        std::vector<double> a(kaa), x(static_cast<std::size_t>(n) * 24);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < 24; ++j)
                x[static_cast<std::size_t>(i) * 24 + j] = kua[static_cast<std::size_t>(j) * n + i];
        for (int col = 0; col < n; ++col) {
            const double pivot = a[static_cast<std::size_t>(col) * n + col];
            for (int r = 0; r < n; ++r) {
                if (r == col) continue;
                const double factor = a[static_cast<std::size_t>(r) * n + col] / pivot;
                for (int k = 0; k < n; ++k)
                    a[static_cast<std::size_t>(r) * n + k] -= factor * a[static_cast<std::size_t>(col) * n + k];
                for (int k = 0; k < 24; ++k)
                    x[static_cast<std::size_t>(r) * 24 + k] -= factor * x[static_cast<std::size_t>(col) * 24 + k];
            }
        }
        for (int r = 0; r < n; ++r)
            for (int k = 0; k < 24; ++k)
                x[static_cast<std::size_t>(r) * 24 + k] /= a[static_cast<std::size_t>(r) * n + r];

        out.assign(24 * 24, 0.0);
        for (int i = 0; i < 24; ++i)
            for (int j = 0; j < 24; ++j) {
                double s = 0;
                for (int k = 0; k < n; ++k)
                    s += kua[static_cast<std::size_t>(i) * n + k] * x[static_cast<std::size_t>(k) * 24 + j];
                out[static_cast<std::size_t>(i) * 24 + j] =
                    kuu[static_cast<std::size_t>(i) * 24 + j] - s;
            }
    };

    std::vector<double> unit(static_cast<std::size_t>(n), 1.0), plain, skewed;
    double conditionPlain = 0, conditionSkewed = 0;
    condense(unit, plain, &conditionPlain);

    // Deliberately awful scaling, spanning twelve orders of magnitude.
    std::vector<double> wild(static_cast<std::size_t>(n));
    for (int j = 0; j < n; ++j) wild[static_cast<std::size_t>(j)] = std::pow(10.0, 2 * j - 6);
    condense(wild, skewed, &conditionSkewed);

    // Vacuity: the rescaling has to have actually wrecked Kaa, or invariance is
    // being asserted about nothing.
    expectTrue("the normalised basis leaves Kaa well conditioned", conditionPlain < 10.0);
    expectTrue("and the deliberate rescaling ruins it", conditionSkewed > 1e12);

    // Vacuity: the enhanced modes have to be doing something, or Kuu alone would
    // pass every comparison below.
    double norm = 0, difference = 0;
    for (std::size_t i = 0; i < plain.size(); ++i) norm = std::max(norm, std::fabs(plain[i]));
    expectTrue("the element has stiffness at all", norm > 0);

    // The independent condensation must reproduce the shipped element...
    double shipped[24 * 24];
    solidshell::elementStiffness(rest, steel, solidshell::Formulation::SolidShell, shipped);
    double worstShipped = 0;
    for (std::size_t i = 0; i < plain.size(); ++i)
        worstShipped = std::max(worstShipped, std::fabs(plain[i] - shipped[i]));
    expectTrue("an independent condensation reproduces the shipped element",
               worstShipped < 1e-9 * norm);

    // ...and the enhanced basis's scaling must not reach the answer.
    for (std::size_t i = 0; i < plain.size(); ++i)
        difference = std::max(difference, std::fabs(plain[i] - skewed[i]));
    std::printf("     enhanced basis: kappa %.2e normalised, %.2e wrecked; element moves %.2e of %.2e\n",
                conditionPlain, conditionSkewed, difference, norm);
    expectTrue("rescaling the enhanced basis cannot change the element", difference < 1e-8 * norm);
}

// --- Multi-point constraints ------------------------------------------------------
//
// `DofExpansion` is the whole of what an assembler has to understand about a
// constraint, so it is tested on its own before anything is solved through it.
// Every refusal here is a modelling error that would otherwise assemble into a
// plausible matrix.

void testDofExpansionRefusesWhatItCannotCompose() {
    const auto refuses = [&](const std::string& what, const std::vector<solidshell::Mpc>& set) {
        const solidshell::DofExpansion expansion(12, set);
        expectTrue(what + " is refused", !expansion.ok());
        expectTrue(what + " says why", !expansion.problem().empty());
    };
    refuses("a slave outside the system", {{12, {0}, {1.0}}});
    refuses("a master outside the system", {{0, {12}, {1.0}}});
    refuses("a weight count that does not match", {{0, {1, 2}, {1.0}}});
    refuses("a constraint with no masters", {{0, {}, {}}});
    refuses("a degree of freedom eliminated twice", {{0, {1}, {1.0}}, {0, {2}, {1.0}}});
    refuses("a constraint that names its own slave", {{0, {0}, {1.0}}});
    // Chained, written both ways round. The second order is the one a check made as
    // the constraints arrive would accept, because when 3 is declared the slave of 6
    // nothing yet knows 6 is a master.
    refuses("a chain", {{3, {6}, {1.0}}, {6, {9}, {1.0}}});
    refuses("the same chain in the other order", {{6, {9}, {1.0}}, {3, {6}, {1.0}}});

    const solidshell::DofExpansion nothing(12, {});
    expectTrue("no constraints is the identity", nothing.ok() && nothing.empty());
    expectEqual("and eliminates nothing", static_cast<long long>(nothing.eliminatedCount()), 0);
    for (std::uint32_t d = 0; d < 12; ++d) {
        expectTrue("every degree of freedom stands for itself", !nothing.eliminated(d));
        expectEqual("as one term", static_cast<long long>(nothing.end(d) - nothing.begin(d)), 1);
        expectTrue("with weight one", nothing.begin(d)->dof == d && nothing.begin(d)->weight == 1.0);
    }

    // `recover` is what turns a solution over the free degrees of freedom back into
    // one over the mesh's, and it is the half of the transformation a solver would
    // still look right without: the free part would be correct and the tied node
    // would sit at zero.
    const solidshell::DofExpansion mixed(12, {{6, {0, 3}, {0.25, 0.75}}});
    expectTrue("a well-formed set is accepted", mixed.ok());
    expectEqual("and eliminates one degree of freedom",
                static_cast<long long>(mixed.eliminatedCount()), 1);
    std::vector<double> values(12, 0.0);
    values[0] = 4.0;
    values[3] = 8.0;
    values[6] = -1.0;  // whatever was there is overwritten, not added to
    mixed.recover(values);
    expectNear("recover reads the slave off its masters", values[6], 0.25 * 4.0 + 0.75 * 8.0,
               1e-15);
}

void testConstrainedSolveIsExactWhereTheConstraintIsTrue() {
    // A patch test with a constraint through it. The plate's whole boundary is
    // prescribed to a linear field; one interior node pair is then eliminated in
    // favour of the two node pairs either side of it along x, with the weights that
    // put it where it actually is. A linear field satisfies that constraint
    // identically, so the constrained solve must return **exactly** the same answer
    // as the unconstrained one -- which is the strongest statement available about a
    // transformation, and it fails for any error in `T`, in `T^T`, or in the band.
    constexpr int n = 4, nz = 3;
    const int sy = n + 1, sz = nz + 1;
    solidshell::HexMesh mesh = solidshell::makePlateMesh(1.0, 1.0, 0.3, n, n, nz);
    const StructuralMaterial steel = ah36Steel();
    // A general linear field: three normal strains, three shears, a rigid rotation
    // and a translation, all at once -- the same field `runPatch` uses, and for the
    // same reason.
    const double gradient[3][3] = {{1.1e-4, 3.0e-5, -2.0e-5},
                                   {-4.0e-5, 0.7e-4, 1.5e-5},
                                   {2.5e-5, -1.0e-5, -0.9e-4}};
    const double offset[3] = {1.0e-3, -2.0e-3, 3.0e-3};
    const auto field = [&](std::size_t node, int axis) {
        const double* p = &mesh.position[node * 3];
        return gradient[axis][0] * p[0] + gradient[axis][1] * p[1] + gradient[axis][2] * p[2] +
               offset[axis];
    };
    // The whole boundary is prescribed, the two through-thickness faces included:
    // pinning only the in-plane perimeter would leave the top and bottom free, where
    // a constant stress state has traction, and the linear field would stop being the
    // exact solution.
    int freeNodes = 0;
    for (int i = 0; i <= n; ++i)
        for (int j = 0; j <= n; ++j)
            for (int k = 0; k <= nz; ++k) {
                const auto node = static_cast<std::size_t>((i * sy + j) * sz + k);
                if (i != 0 && i != n && j != 0 && j != n && k != 0 && k != nz) {
                    ++freeNodes;
                    continue;
                }
                for (int axis = 0; axis < 3; ++axis) mesh.pin(node, axis, field(node, axis));
            }
    expectTrue("the patch has free interior nodes to constrain", freeNodes >= 8);

    const std::vector<double> noLoad(mesh.nodeCount() * 3, 0.0);
    std::vector<double> free;
    std::string problem;
    expectTrue("the unconstrained patch solves: " + problem,
               solidshell::solveStatic(mesh, steel, solidshell::Formulation::SolidShell, noLoad,
                                       free, &problem));

    // The interior node at (i, j, k) = (2, 2, 1), eliminated in favour of (1, 2, 1)
    // and (3, 2, 1): it is exactly halfway between them on a regular grid, so the
    // weights are a half each and a linear field satisfies the constraint
    // identically.
    const auto index = [&](int i, int j, int k) {
        return static_cast<std::uint32_t>((i * sy + j) * sz + k);
    };
    std::vector<solidshell::Mpc> constrained;
    for (int k = 0; k < 3; ++k)
        constrained.push_back({index(2, 2, 1) * 3 + static_cast<std::uint32_t>(k),
                               {index(1, 2, 1) * 3 + static_cast<std::uint32_t>(k),
                                index(3, 2, 1) * 3 + static_cast<std::uint32_t>(k)},
                               {0.5, 0.5}});
    std::vector<double> tied;
    expectTrue("the constrained patch solves: " + problem,
               solidshell::solveStatic(mesh, steel, solidshell::Formulation::SolidShell, {},
                                       constrained, noLoad, tied, &problem));

    double worstFree = 0, worstTied = 0, scale = 0;
    for (std::size_t node = 0; node < mesh.nodeCount(); ++node)
        for (int k = 0; k < 3; ++k) {
            const double want = field(node, k);
            scale = std::max(scale, std::fabs(want));
            worstFree = std::max(worstFree, std::fabs(free[node * 3 + static_cast<std::size_t>(k)] - want));
            worstTied = std::max(worstTied, std::fabs(tied[node * 3 + static_cast<std::size_t>(k)] - want));
        }
    std::printf("     constrained patch: unconstrained %.2e m, constrained %.2e m, of %.2e m\n",
                worstFree, worstTied, scale);
    expectTrue("the unconstrained patch is exact", worstFree < 1e-12 * scale);
    // The eliminated node is filled from its masters, so this also says `recover`
    // ran: without it the tied node would sit at zero and this would be `scale`.
    expectTrue("and so is the constrained one, because a linear field satisfies the tie",
               worstTied < 1e-12 * scale);

    // The vacuity guard, and it is the whole reason the test means anything: a
    // constraint the exact solution does **not** satisfy has to change the answer.
    // Tie the same node to the two either side with the weights swapped in a way
    // that puts it a quarter of a bay off where it is.
    for (auto& mpc : constrained) mpc.weight = {0.25, 0.75};
    std::vector<double> wrong;
    expectTrue("the off-centre constraint also solves",
               solidshell::solveStatic(mesh, steel, solidshell::Formulation::SolidShell, {},
                                       constrained, noLoad, wrong, &problem));
    double moved = 0;
    for (std::size_t d = 0; d < wrong.size(); ++d)
        moved = std::max(moved, std::fabs(wrong[d] - tied[d]));
    expectTrue("a constraint the exact field does not satisfy moves the answer",
               moved > 1e-4 * scale);

    // **A load applied at an eliminated degree of freedom has to go to its masters
    // by the transpose of the constraint**, or the virtual work is not the same and
    // the tied point silently carries nothing. Nothing above would see it: the patch
    // is driven entirely by its boundary and its load vector is zero. So the load is
    // put where only that path can carry it, and compared against the same load
    // spread over the masters by hand -- which must give the identical field,
    // because `f^T u` is the same number written two ways.
    for (auto& mpc : constrained) mpc.weight = {0.5, 0.5};
    std::vector<double> atSlave(mesh.nodeCount() * 3, 0.0), atMasters = atSlave;
    for (int k = 0; k < 3; ++k) {
        const double f = 1.0e4 * (k + 1);
        atSlave[index(2, 2, 1) * 3 + static_cast<std::size_t>(k)] = f;
        atMasters[index(1, 2, 1) * 3 + static_cast<std::size_t>(k)] = 0.5 * f;
        atMasters[index(3, 2, 1) * 3 + static_cast<std::size_t>(k)] = 0.5 * f;
    }
    std::vector<double> pushedAtSlave, pushedAtMasters;
    expectTrue("a load at the tied node solves",
               solidshell::solveStatic(mesh, steel, solidshell::Formulation::SolidShell, {},
                                       constrained, atSlave, pushedAtSlave, &problem));
    expectTrue("and so does the same load at its masters",
               solidshell::solveStatic(mesh, steel, solidshell::Formulation::SolidShell, {},
                                       constrained, atMasters, pushedAtMasters, &problem));
    double loadGap = 0, loadScale = 0;
    for (std::size_t d = 0; d < pushedAtSlave.size(); ++d) {
        loadScale = std::max(loadScale, std::fabs(pushedAtSlave[d] - tied[d]));
        loadGap = std::max(loadGap, std::fabs(pushedAtSlave[d] - pushedAtMasters[d]));
    }
    std::printf("     load through the tie: two routes agree to %.2e m, having moved %.2e m\n",
                loadGap, loadScale);
    // The vacuity guard: the load has to have done something, or "the two agree" is
    // a statement about two copies of the boundary-driven answer.
    expectTrue("the load moved the plate", loadScale > 1e-8);
    expectTrue("a load at a tied node is its masters' load", loadGap < 1e-12 * loadScale);

    // And the two ways of saying "this degree of freedom is spoken for" may not both
    // be used at once: the pin and the tie would each claim it and there is no
    // reading of "both" that is not a guess.
    solidshell::HexMesh clash = mesh;
    clash.pin(index(2, 2, 1), 0, 0.0);
    for (auto& mpc : constrained) mpc.weight = {0.5, 0.5};
    std::vector<double> refused;
    expectTrue("a degree of freedom both pinned and tied is refused",
               !solidshell::solveStatic(clash, steel, solidshell::Formulation::SolidShell, {},
                                        constrained, noLoad, refused, &problem));
    expectTrue("and says so", problem.find("prescribed") != std::string::npos);
}

// **Two `criticalTimestep` functions, and until now two contracts.**
//
// `solidshell::criticalTimestep` has always ended
// `std::isfinite(smallest) ? smallest : 0.0`. `fem::criticalTimestep` computed the
// same kind of quantity for a tet mesh, ended with a bare division, and had no test
// anywhere in the tree -- so it answered `+infinity` for an empty mesh, and a caller
// writing `dt = std::min(dt, criticalTimestep(...))` kept its own `dt` and never
// learned. Zero is the sentinel because `!(dt > 0)` catches it and catches a NaN with
// it, while **no comparison against an infinity is false**.
//
// The one worth the test is the NaN, because it was silent rather than absurd.
// `TetMesh::shortestEdge` folded with `std::min(shortest, sqrt(d2))`, and
// `std::min(a, b)` is `b < a ? b : a` -- with the accumulator *first*, `NaN < shortest`
// is false and the NaN is dropped from the minimum. The shortest edge came back as the
// shortest of the finite tets, and the step looked entirely reasonable for geometry
// that cannot be integrated. This asserts the fold propagates it now, which is what
// gives the sentinel below anything to catch.
void testBothCriticalTimestepsShareTheirSentinel() {
    std::printf("\n   two criticalTimestep functions, one sentinel\n");
    const fem::Material tetSteel;
    // `makeBoxTetMesh` lays out nodes and connectivity; `tetCount()` reads
    // `restVolume.size()`, which stays empty until `computeRestState` runs. Without
    // it `shortestEdge` iterates nothing and returns `+infinity` -- which is how the
    // vacuity guard below caught this test's own setup before it caught anything else.
    fem::TetMesh box = fem::makeBoxTetMesh(2.0f, 0.5f, 0.05f, 20, 4, 2);
    box.computeRestState(tetSteel);
    expectEqual("the reference mesh is the 960-tet box fem_spike's first row uses",
                static_cast<int>(box.tetCount()), 960);

    // Vacuity first: a real mesh must give a real answer, or every zero below
    // would hold for a routine that always said zero.
    const float good = fem::criticalTimestep(box, tetSteel, 0.5f);
    expectTrue("a real tet mesh has a positive finite step", good > 0 && std::isfinite(good));
    expectTrue("and its shortest edge is a real length",
               box.shortestEdge() > 0 && std::isfinite(box.shortestEdge()));

    // Empty: nothing to cross. Both solvers now answer zero rather than infinity.
    expectNear("an empty tet mesh answers zero, not infinity",
               fem::criticalTimestep(fem::TetMesh{}, tetSteel, 0.5f), 0.0f, 0.0f);
    expectNear("which is what the hex solver has always answered",
               criticalTimestep(HexMesh{}, steelMaterial(), Formulation::SolidShell, 0.5), 0.0,
               0.0);

    // A single NaN coordinate. The fold must carry it, and the sentinel must catch it.
    fem::TetMesh poisoned = box;
    poisoned.position[0] = std::numeric_limits<float>::quiet_NaN();
    expectTrue("one NaN node makes the shortest edge a NaN rather than vanishing",
               std::isnan(poisoned.shortestEdge()));
    expectNear("so the step it feeds answers zero", fem::criticalTimestep(poisoned, tetSteel, 0.5f),
               0.0f, 0.0f);

    // A material with no wave speed divides by zero. `sqrt(E/rho)` is unguarded, so
    // zero modulus gives a wave speed of zero and a step of +infinity.
    fem::Material noModulus = tetSteel;
    noModulus.youngsModulus = 0.0f;
    expectNear("a material with no wave speed answers zero too",
               fem::criticalTimestep(box, noModulus, 0.5f), 0.0f, 0.0f);
}

}  // namespace

void runSolidShellTests() {
    testElementStiffnessRank();
    testLumpedMass();
    testInvertedElementIsCaught();
    testCollapsedElementIsAWedgeAndNotAFold();
    testElementShapeSamplesTheQuadratureAndNotOnlyTheCorners();
    testCollapsedElementSolvesAndPassesThePatchTest();
    testPressureLoadResultant();
    testFreeExpansionCarriesNoStress();
    testEigenstrainAndTheEnhancedModes();
    testTaperedElementCannotExpandFreely();
    testRestraintAndNotTemperatureIsWhatGeneratesStress();
    testNullEigenstrainIsTheSameArithmetic();
    testEigenstrainThroughThePlasticUpdate();
    testFiniteRotationCarriesNoForce();
    testFrameIndifference();
    testInternalForceSignAndEquilibrium();
    testPatchTest();
    testCantileverAgainstBeamTheory();
    testCylindricalBendingAndThicknessLocking();
    testThroughThicknessStress();
    testInPlaneBendingOnDistortedElements();
    testTrapezoidalElementLimit();
    testStripsUnderPressure();
    testPlateSeriesSolutions();
    testLockingDemonstration();
    testExplicitStabilityLimit();
    testCostPerElement();
    testCachedRestFormsAreBitIdentical();
    testEnhancedBasisScalingCannotChangeTheElement();
    testDofExpansionRefusesWhatItCannotCompose();
    testConstrainedSolveIsExactWhereTheConstraintIsTrue();
    testBothCriticalTimestepsShareTheirSentinel();
}
