// SPDX-License-Identifier: MIT
#include "solid_shell.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <limits>
#include <map>

namespace sim::solidshell {
namespace {

// Natural coordinates of the eight nodes, in the ordering the header fixes.
constexpr double kXi[kNodes]  = {-1, +1, +1, -1, -1, +1, +1, -1};
constexpr double kEta[kNodes] = {-1, -1, +1, +1, -1, -1, +1, +1};
constexpr double kZta[kNodes] = {-1, -1, -1, -1, +1, +1, +1, +1};

// Which cures each formulation switches on.
struct Cures {
    bool assumedShear;
    bool assumedThickness;
    bool enhanced;
};

Cures curesFor(Formulation form) {
    switch (form) {
        case Formulation::Displacement: return {false, false, false};
        case Formulation::AssumedNaturalStrain: return {true, true, false};
        case Formulation::SolidShell: break;
    }
    return {true, true, true};
}

// --- 3x3 helpers, row-major here (unlike fem.cpp's column-major M3, which exists
// to mirror GLSL's mat3; nothing here talks to a shader). --------------------

double determinant3(const double m[3][3]) {
    return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
           m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
           m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
}

// Returns the determinant; `inv` is left untouched when it is zero.
double invert3(const double m[3][3], double inv[3][3]) {
    const double det = determinant3(m);
    if (det == 0.0) return 0.0;
    const double d = 1.0 / det;
    inv[0][0] = (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * d;
    inv[0][1] = (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * d;
    inv[0][2] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * d;
    inv[1][0] = (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * d;
    inv[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * d;
    inv[1][2] = (m[0][2] * m[1][0] - m[0][0] * m[1][2]) * d;
    inv[2][0] = (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * d;
    inv[2][1] = (m[0][1] * m[2][0] - m[0][0] * m[2][1]) * d;
    inv[2][2] = (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * d;
    return det;
}

// --- Element geometry at a natural point --------------------------------------

struct Point {
    double dN[kNodes][3];  // dN_a / d(xi, eta, zeta)
    double jac[3][3];      // jac[i][k] = dx_i / dxi_k, so column k is the base vector g_k
};

void evaluate(const double nodes[kDof], double xi, double eta, double zta, Point& p) {
    for (int a = 0; a < kNodes; ++a) {
        const double x = 1.0 + xi * kXi[a];
        const double y = 1.0 + eta * kEta[a];
        const double z = 1.0 + zta * kZta[a];
        p.dN[a][0] = 0.125 * kXi[a] * y * z;
        p.dN[a][1] = 0.125 * kEta[a] * x * z;
        p.dN[a][2] = 0.125 * kZta[a] * x * y;
    }
    for (int i = 0; i < 3; ++i)
        for (int k = 0; k < 3; ++k) {
            double s = 0.0;
            for (int a = 0; a < kNodes; ++a) s += p.dN[a][k] * nodes[a * 3 + i];
            p.jac[i][k] = s;
        }
}

// Strain-displacement in *covariant* (natural) components, Voigt-ordered
// [E11, E22, E33, 2E12, 2E23, 2E31] with 1,2,3 = xi, eta, zeta.
//
//   E_kl = 1/2 (g_k . du/dxi_l + g_l . du/dxi_k)
//
// so dE_kl/du_a = 1/2 (N_a,l g_k + N_a,k g_l). The assumed-strain modifications
// all act on these rows, before the transformation to Cartesian, because that is
// where the parasitic terms are identifiable as whole components.
void naturalB(const Point& p, double b[6][kDof]) {
    for (int a = 0; a < kNodes; ++a)
        for (int i = 0; i < 3; ++i) {
            const int c = a * 3 + i;
            const double g0 = p.jac[i][0], g1 = p.jac[i][1], g2 = p.jac[i][2];
            const double n0 = p.dN[a][0], n1 = p.dN[a][1], n2 = p.dN[a][2];
            b[0][c] = n0 * g0;
            b[1][c] = n1 * g1;
            b[2][c] = n2 * g2;
            b[3][c] = n1 * g0 + n0 * g1;
            b[4][c] = n2 * g1 + n1 * g2;
            b[5][c] = n0 * g2 + n2 * g0;
        }
}

// Voigt transformation from natural to Cartesian strain: eps = A^T E A with
// A = J^-1. Built by pushing the six unit natural strains through the tensor
// transformation rather than transcribing a closed-form 6x6, which is where this
// kind of code usually acquires a transposition that only shows on a distorted
// element.
void voigtTransform(const double a[3][3], double t[6][6]) {
    static constexpr int kRow[6] = {0, 1, 2, 0, 1, 2};
    static constexpr int kCol[6] = {0, 1, 2, 1, 2, 0};
    for (int j = 0; j < 6; ++j) {
        double e[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
        // Voigt entries 3..5 are engineering shears, so the tensor component is
        // half of them.
        const double v = j < 3 ? 1.0 : 0.5;
        e[kRow[j]][kCol[j]] = v;
        e[kCol[j]][kRow[j]] = v;

        double ea[3][3];
        for (int r = 0; r < 3; ++r)
            for (int c = 0; c < 3; ++c) {
                double s = 0.0;
                for (int k = 0; k < 3; ++k) s += e[r][k] * a[k][c];
                ea[r][c] = s;
            }
        double m[3][3];
        for (int r = 0; r < 3; ++r)
            for (int c = 0; c < 3; ++c) {
                double s = 0.0;
                for (int k = 0; k < 3; ++k) s += a[k][r] * ea[k][c];
                m[r][c] = s;
            }
        t[0][j] = m[0][0];
        t[1][j] = m[1][1];
        t[2][j] = m[2][2];
        t[3][j] = 2.0 * m[0][1];
        t[4][j] = 2.0 * m[1][2];
        t[5][j] = 2.0 * m[2][0];
    }
}

void isotropic(const StructuralMaterial& material, double c[6][6]) {
    const double e = material.youngsModulus, nu = material.poissonRatio;
    const double lambda = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    const double mu = e / (2.0 * (1.0 + nu));
    for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j) c[i][j] = 0.0;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) c[i][j] = lambda;
        c[i][i] += 2.0 * mu;
        c[3 + i][3 + i] = mu;
    }
}

// The seven enhanced modes, in natural Voigt components. Every entry is odd in at
// least one natural coordinate, so each integrates to zero over the reference
// cube -- which, with the Simo-Rifai det(J0)/det(J) scaling below, is precisely
// the condition that makes EAS pass the patch test on *any* geometry.
//
//   alpha 0,1     xi on E11, eta on E22          in-plane bending
//   alpha 2,3,4   zeta, xi*zeta, eta*zeta on E33 thickness strain linear through
//                                                the plate, varying in plane
//   alpha 5,6     xi and eta on 2E12             in-plane shear of a distorted element
void enhancedM(double xi, double eta, double zta, double m[6][kEas]) {
    for (int i = 0; i < 6; ++i)
        for (int j = 0; j < kEas; ++j) m[i][j] = 0.0;
    m[0][0] = xi;
    m[1][1] = eta;
    m[2][2] = zta;
    m[2][3] = xi * zta;
    m[2][4] = eta * zta;
    m[3][5] = xi;
    m[3][6] = eta;
}

// `RestForms` is everything an element's stiffness, stress and internal force
// need, and it is public because it is also the whole of what an explicit solver
// should be caching. The shape functions are wanted only by `elementMass`, which
// runs once per element per solve, so they are an optional out-parameter rather
// than 512 bytes an element carried for the life of a zone.
using Shape = double[kGauss][kNodes];

// ANS sampling points. Transverse shear (Dvorkin-Bathe):
//   2E_23 (eta-zeta) is taken constant in eta and linear in xi, from (xi, eta) =
//   (-1, 0) and (+1, 0);
//   2E_31 (zeta-xi) is taken constant in xi and linear in eta, from (0, -1) and
//   (0, +1).
// Under cylindrical bending about y the parasitic shear is proportional to xi, so
// sampling 2E_31 at xi = 0 annihilates it exactly.
//
// Thickness strain (Betsch-Stein): E_33 is biquadratic in (xi, eta) for a warped
// element and the biquadratic part is curvature thickness locking, so it is
// sampled at the four in-plane corners and interpolated bilinearly. On a flat
// prismatic element E_33 is already bilinear, so this is exact and costs only the
// sampling.
//
// Every sample keeps the integration point's own zeta rather than being taken at
// the mid-surface: for the prismatic plate geometry this element is for, the
// covariant shear is zeta-independent anyway, so the two agree, and keeping zeta
// is the smaller assumption for the thick, warped elements a crush zone produces.
// A 2x2x2 rule has only two distinct zeta levels, so the samples are taken twice
// per element rather than once per integration point -- 8 evaluations instead of
// 64, which is most of the assumed-strain machinery's cost.
constexpr double kShearSample[4][2] = {{-1, 0}, {+1, 0}, {0, -1}, {0, +1}};
constexpr double kCornerSample[4][2] = {{-1, -1}, {+1, -1}, {+1, +1}, {-1, +1}};

struct AnsSamples {
    double shear[2][4][kDof];      // [zeta level][sample]: row 4 from 0,1 and row 5 from 2,3
    double thickness[2][4][kDof];  // [zeta level][in-plane corner]: row 2
};

void sampleAns(const double nodes[kDof], double zeta, bool shear, bool thickness,
               AnsSamples& out, int level) {
    for (int s = 0; s < 4; ++s) {
        if (shear) {
            Point p;
            evaluate(nodes, kShearSample[s][0], kShearSample[s][1], zeta, p);
            double b[6][kDof];
            naturalB(p, b);
            for (int c = 0; c < kDof; ++c) out.shear[level][s][c] = b[s < 2 ? 4 : 5][c];
        }
        if (thickness) {
            Point p;
            evaluate(nodes, kCornerSample[s][0], kCornerSample[s][1], zeta, p);
            double b[6][kDof];
            naturalB(p, b);
            for (int c = 0; c < kDof; ++c) out.thickness[level][s][c] = b[2][c];
        }
    }
}

void computeForms(const double nodes[kDof], Formulation form, RestForms& out,
                  Shape* shape = nullptr) {
    const Cures cures = curesFor(form);
    out.easCount = cures.enhanced ? kEas : 0;
    out.ok = true;

    Point centre;
    evaluate(nodes, 0.0, 0.0, 0.0, centre);
    double a0[3][3];
    const double det0 = invert3(centre.jac, a0);
    if (!(det0 > 0.0)) {
        out.ok = false;
        return;
    }
    // The polar decomposition wants exactly this matrix every step, and
    // `elementRotation` re-derives it every step. Keeping it here is what lets the
    // cached path be bit-identical rather than merely equivalent.
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) out.restJacobianInverse[i * 3 + j] = a0[i][j];
    double t0[6][6];
    voigtTransform(a0, t0);

    const double q = 1.0 / std::sqrt(3.0);
    AnsSamples samples{};
    if (cures.assumedShear || cures.assumedThickness)
        for (int level = 0; level < 2; ++level)
            sampleAns(nodes, level ? q : -q, cures.assumedShear, cures.assumedThickness, samples,
                      level);

    for (int gp = 0; gp < kGauss; ++gp) {
        const double xi = (gp & 1) ? q : -q;
        const double eta = (gp & 2) ? q : -q;
        const double zta = (gp & 4) ? q : -q;
        const int level = (gp & 4) ? 1 : 0;

        Point p;
        evaluate(nodes, xi, eta, zta, p);
        double a[3][3];
        const double det = invert3(p.jac, a);
        if (!(det > 0.0)) {
            out.ok = false;
            return;
        }
        out.weight[gp] = det;  // 2x2x2 Gauss weights are all 1

        if (shape != nullptr)
            for (int n = 0; n < kNodes; ++n)
                (*shape)[gp][n] = 0.125 * (1.0 + xi * kXi[n]) * (1.0 + eta * kEta[n]) *
                                  (1.0 + zta * kZta[n]);

        double bn[6][kDof];
        naturalB(p, bn);
        if (cures.assumedShear) {
            const double wa = 0.5 * (1.0 - xi), wc = 0.5 * (1.0 + xi);
            const double wb = 0.5 * (1.0 - eta), wd = 0.5 * (1.0 + eta);
            for (int c = 0; c < kDof; ++c) {
                bn[4][c] = wa * samples.shear[level][0][c] + wc * samples.shear[level][1][c];
                bn[5][c] = wb * samples.shear[level][2][c] + wd * samples.shear[level][3][c];
            }
        }
        if (cures.assumedThickness) {
            for (int c = 0; c < kDof; ++c) bn[2][c] = 0.0;
            for (int s = 0; s < 4; ++s) {
                const double w = 0.25 * (1.0 + xi * kCornerSample[s][0]) *
                                 (1.0 + eta * kCornerSample[s][1]);
                for (int c = 0; c < kDof; ++c) bn[2][c] += w * samples.thickness[level][s][c];
            }
        }

        double t[6][6];
        voigtTransform(a, t);
        for (int i = 0; i < 6; ++i)
            for (int c = 0; c < kDof; ++c) {
                double s = 0.0;
                for (int k = 0; k < 6; ++k) s += t[i][k] * bn[k][c];
                out.b[gp][i][c] = s;
            }

        if (out.easCount == 0) {
            for (int i = 0; i < 6; ++i)
                for (int j = 0; j < kEas; ++j) out.g[gp][i][j] = 0.0;
            continue;
        }
        double m[6][kEas];
        enhancedM(xi, eta, zta, m);
        const double scale = det0 / det;
        for (int i = 0; i < 6; ++i)
            for (int j = 0; j < kEas; ++j) {
                double s = 0.0;
                for (int k = 0; k < 6; ++k) s += t0[i][k] * m[k][j];
                out.g[gp][i][j] = scale * s;
            }
    }

    // --- Normalise the enhanced modes ------------------------------------------
    //
    // The enhanced modes are a *basis* and their scaling is a free choice: scaling
    // column j of G by s and alpha_j by 1/s leaves G alpha, and therefore the
    // element, exactly unchanged. Unnormalised they are the natural-coordinate
    // polynomials, and for a thin element the thickness modes pick up the
    // Jacobian's t/2 against the in-plane h/2, so Kaa's condition number goes as
    // **(h/t)^4** -- measured on this element at 2.19e3, 3.50e4, 2.84e6 and 4.54e7
    // for h/t of 5, 10, 30 and 60, which is that power to three figures.
    //
    // At the ferry's 300 mm elements on 10 mm plating that is 2.8e6, and it is why
    // the float GPU path had to equilibrate Kaa in the shader before factoring it
    // and still lost the enhanced parameters entirely (`07-fem-spike-findings.md`
    // §8). Fixing it here fixes it for every path.
    //
    // **Why this is safe, which is the part that is not obvious.** `alpha` is
    // persistent per-element state -- `ElementPlasticState::enhanced` carries it
    // between steps -- so a scale that is recomputed differently on any step would
    // silently reinterpret the history already stored. This one is a function of
    // the **rest** configuration alone, so it is constant for the element's life
    // and the stored alpha stays in one basis by construction. Deriving it from the
    // current configuration is the version that looks equivalent and is not.
    //
    // Nothing downstream needs to know: `Kua` and `Kaa` are both built from `G`, so
    // scaling `G` scales them consistently, and alpha comes back out of a system
    // built the same way. That is what makes this a basis change rather than an
    // approximation, and the tests assert it as an identity on the recovered
    // displacement rather than as a tolerance.
    //
    // The norm is the weighted L2 one over the Gauss points -- the material is not
    // available here, and does not need to be: Kaa is diagonal to within the
    // material's own anisotropy, so removing the geometric spread is what matters.
    if (out.easCount > 0) {
        for (int j = 0; j < kEas; ++j) {
            double norm = 0.0;
            for (int gp = 0; gp < kGauss; ++gp)
                for (int i = 0; i < 6; ++i)
                    norm += out.weight[gp] * out.g[gp][i][j] * out.g[gp][i][j];
            if (!(norm > 0.0)) continue;   // an unused mode stays as it is
            const double inverse = 1.0 / std::sqrt(norm);
            for (int gp = 0; gp < kGauss; ++gp)
                for (int i = 0; i < 6; ++i) out.g[gp][i][j] *= inverse;
        }
    }
}

struct Blocks {
    double kuu[kDof * kDof];
    double kua[kDof * kEas];
    double kaa[kEas * kEas];
    int easCount;
};

void computeBlocks(const RestForms& f, const double c[6][6], Blocks& out) {
    out.easCount = f.easCount;
    std::fill(std::begin(out.kuu), std::end(out.kuu), 0.0);
    std::fill(std::begin(out.kua), std::end(out.kua), 0.0);
    std::fill(std::begin(out.kaa), std::end(out.kaa), 0.0);

    for (int gp = 0; gp < kGauss; ++gp) {
        const double w = f.weight[gp];
        double cb[6][kDof];
        for (int i = 0; i < 6; ++i)
            for (int j = 0; j < kDof; ++j) {
                double s = 0.0;
                for (int k = 0; k < 6; ++k) s += c[i][k] * f.b[gp][k][j];
                cb[i][j] = s;
            }
        for (int i = 0; i < kDof; ++i)
            for (int j = 0; j < kDof; ++j) {
                double s = 0.0;
                for (int k = 0; k < 6; ++k) s += f.b[gp][k][i] * cb[k][j];
                out.kuu[i * kDof + j] += w * s;
            }
        if (f.easCount == 0) continue;

        double cg[6][kEas];
        for (int i = 0; i < 6; ++i)
            for (int j = 0; j < kEas; ++j) {
                double s = 0.0;
                for (int k = 0; k < 6; ++k) s += c[i][k] * f.g[gp][k][j];
                cg[i][j] = s;
            }
        for (int i = 0; i < kDof; ++i)
            for (int j = 0; j < kEas; ++j) {
                double s = 0.0;
                for (int k = 0; k < 6; ++k) s += f.b[gp][k][i] * cg[k][j];
                out.kua[i * kEas + j] += w * s;
            }
        for (int i = 0; i < kEas; ++i)
            for (int j = 0; j < kEas; ++j) {
                double s = 0.0;
                for (int k = 0; k < 6; ++k) s += f.g[gp][k][i] * cg[k][j];
                out.kaa[i * kEas + j] += w * s;
            }
    }
}

// Dense Cholesky solve of the small enhanced-parameter system, n <= kEas.
// False on a non-positive pivot, which would mean the enhanced modes are not
// linearly independent and the condensation is not defined.
bool solveSmall(const double a[kEas * kEas], int n, double* rhs, int columns) {
    double l[kEas * kEas] = {};
    for (int i = 0; i < n; ++i)
        for (int j = 0; j <= i; ++j) {
            double s = a[i * kEas + j];
            for (int k = 0; k < j; ++k) s -= l[i * kEas + k] * l[j * kEas + k];
            if (i == j) {
                if (!(s > 0.0)) return false;
                l[i * kEas + i] = std::sqrt(s);
            } else {
                l[i * kEas + j] = s / l[j * kEas + j];
            }
        }
    for (int c = 0; c < columns; ++c) {
        for (int i = 0; i < n; ++i) {
            double s = rhs[i * columns + c];
            for (int k = 0; k < i; ++k) s -= l[i * kEas + k] * rhs[k * columns + c];
            rhs[i * columns + c] = s / l[i * kEas + i];
        }
        for (int i = n - 1; i >= 0; --i) {
            double s = rhs[i * columns + c];
            for (int k = i + 1; k < n; ++k) s -= l[k * kEas + i] * rhs[k * columns + c];
            rhs[i * columns + c] = s / l[i * kEas + i];
        }
    }
    return true;
}

// Rotation part of F by Higham's Newton iteration R <- (R + R^-T)/2, as fem.cpp
// does. fem.cpp uses a fixed four iterations because its GPU kernel has to agree
// step for step; there is no kernel here, so this runs to convergence in double
// instead -- a finite rotation has to give *exactly* zero force, and a fixed
// count leaves a residue that a rigid-body test would then have to tolerate.
// Truncating it to one iteration is caught by the frame-indifference test.
//
// The 1e-16 is *not* load-bearing and is not worth tightening or defending:
// convergence is quadratic, so by the time the step falls below 1e-6 the answer is
// already at machine precision. Mutation testing confirmed it -- loosening this to
// 1e-6 moves the frame-indifference residual from 8.3e-13 to 5.9e-13, which is to
// say it moves it around inside the rounding floor. The strict bound simply buys
// one more iteration on an element that is nearly degenerate, and that is cheap.
void polarRotation(const double f[3][3], double r[3][3]) {
    std::memcpy(r, f, sizeof(double) * 9);
    for (int iteration = 0; iteration < 64; ++iteration) {
        double inv[3][3];
        if (invert3(r, inv) == 0.0) return;
        double worst = 0.0;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) {
                const double next = 0.5 * (r[i][j] + inv[j][i]);
                worst = std::max(worst, std::abs(next - r[i][j]));
                r[i][j] = next;
            }
        if (worst < 1e-16) return;
    }
}

}  // namespace

const char* name(Formulation form) {
    switch (form) {
        case Formulation::Displacement: return "displacement hex";
        case Formulation::AssumedNaturalStrain: return "ANS hex";
        case Formulation::SolidShell: break;
    }
    return "solid-shell";
}

void elementStiffness(const double nodes[kDof], const StructuralMaterial& material,
                      Formulation form, double out[kDof * kDof]) {
    std::fill(out, out + kDof * kDof, 0.0);
    RestForms forms;
    computeForms(nodes, form, forms);
    if (!forms.ok) return;

    double c[6][6];
    isotropic(material, c);
    Blocks blocks;
    computeBlocks(forms, c, blocks);

    std::memcpy(out, blocks.kuu, sizeof(double) * kDof * kDof);
    if (blocks.easCount == 0) return;

    // Static condensation: K = Kuu - Kua Kaa^-1 Kua^T.
    double x[kEas * kDof];
    for (int i = 0; i < blocks.easCount; ++i)
        for (int j = 0; j < kDof; ++j) x[i * kDof + j] = blocks.kua[j * kEas + i];
    if (!solveSmall(blocks.kaa, blocks.easCount, x, kDof)) return;
    for (int i = 0; i < kDof; ++i)
        for (int j = 0; j < kDof; ++j) {
            double s = 0.0;
            for (int k = 0; k < blocks.easCount; ++k) s += blocks.kua[i * kEas + k] * x[k * kDof + j];
            out[i * kDof + j] -= s;
        }
}

namespace {

// `int B^T C eps_star dV` and `int G^T C eps_star dV` -- the two halves of an
// eigenstrain's equivalent force, before condensation. One place, so that
// `elementStress`'s enhanced recovery and `elementThermalLoad`'s nodal force
// cannot be built from two different quadratures of the same integral.
//
// `fa` is **identically zero for an eigenstrain constant over the element**, on any
// geometry, because `int G^T dV = 0` is the condition that makes EAS pass the patch
// test -- see the note at `elementStress` in the header, and the test that asserts
// the integral rather than assuming it. It is computed rather than skipped because
// an eigenstrain that varies within the element is the next refinement and would
// make it real, and because a loop that is provably zero costs nothing next to the
// Cholesky it feeds.
void eigenLoads(const RestForms& f, const double c[6][6], const double eigen[6],
                double fu[kDof], double fa[kEas]) {
    std::fill(fu, fu + kDof, 0.0);
    std::fill(fa, fa + kEas, 0.0);
    double ce[6];
    for (int i = 0; i < 6; ++i) {
        double s = 0.0;
        for (int k = 0; k < 6; ++k) s += c[i][k] * eigen[k];
        ce[i] = s;
    }
    for (int gp = 0; gp < kGauss; ++gp) {
        const double w = f.weight[gp];
        for (int j = 0; j < kDof; ++j) {
            double s = 0.0;
            for (int i = 0; i < 6; ++i) s += f.b[gp][i][j] * ce[i];
            fu[j] += w * s;
        }
        for (int k = 0; k < f.easCount; ++k) {
            double s = 0.0;
            for (int i = 0; i < 6; ++i) s += f.g[gp][i][k] * ce[i];
            fa[k] += w * s;
        }
    }
}

}  // namespace

void elementStress(const double nodes[kDof], const double displacement[kDof],
                   const StructuralMaterial& material, Formulation form,
                   double out[kGauss * 6], const double eigenstrain[6]) {
    std::fill(out, out + kGauss * 6, 0.0);
    RestForms forms;
    computeForms(nodes, form, forms);
    if (!forms.ok) return;

    double c[6][6];
    isotropic(material, c);
    Blocks blocks;
    computeBlocks(forms, c, blocks);

    // Recover the enhanced parameters the element would have settled on:
    // alpha = Kaa^-1 (int G^T C eps_star dV - Kua^T u). Without an eigenstrain that
    // is -Kaa^-1 Kua^T u, and for a linear displacement field it is exactly zero,
    // which is what makes the patch test a statement about the element and not
    // about the condensation.
    //
    // The `-s` is written out rather than reached through `0.0 - s` so that a null
    // eigenstrain is the same arithmetic as before this parameter existed, signed
    // zeros included.
    double fu[kDof], fa[kEas];
    if (eigenstrain != nullptr) eigenLoads(forms, c, eigenstrain, fu, fa);
    double alpha[kEas] = {};
    if (blocks.easCount > 0) {
        for (int i = 0; i < blocks.easCount; ++i) {
            double s = 0.0;
            for (int j = 0; j < kDof; ++j) s += blocks.kua[j * kEas + i] * displacement[j];
            alpha[i] = -s;
            if (eigenstrain != nullptr) alpha[i] += fa[i];
        }
        if (!solveSmall(blocks.kaa, blocks.easCount, alpha, 1))
            std::fill(std::begin(alpha), std::end(alpha), 0.0);
    }

    for (int gp = 0; gp < kGauss; ++gp) {
        double strain[6];
        for (int i = 0; i < 6; ++i) {
            double s = 0.0;
            for (int j = 0; j < kDof; ++j) s += forms.b[gp][i][j] * displacement[j];
            for (int k = 0; k < blocks.easCount; ++k) s += forms.g[gp][i][k] * alpha[k];
            strain[i] = eigenstrain != nullptr ? s - eigenstrain[i] : s;
        }
        for (int i = 0; i < 6; ++i) {
            double s = 0.0;
            for (int k = 0; k < 6; ++k) s += c[i][k] * strain[k];
            out[gp * 6 + i] = s;
        }
    }
}

void elementThermalLoad(const RestForms& forms, const StructuralMaterial& material,
                        const double eigenstrain[6], double out[kDof]) {
    std::fill(out, out + kDof, 0.0);
    if (!forms.ok || eigenstrain == nullptr) return;

    double c[6][6];
    isotropic(material, c);
    double fu[kDof], fa[kEas];
    eigenLoads(forms, c, eigenstrain, fu, fa);
    std::memcpy(out, fu, sizeof(double) * kDof);
    if (forms.easCount == 0) return;

    // Condense: f = fu - Kua Kaa^-1 fa, the same elimination `elementStiffness`
    // performs on the stiffness and necessarily the same Kaa, or the pair would
    // not solve the system it claims to.
    Blocks blocks;
    computeBlocks(forms, c, blocks);
    double x[kEas];
    for (int k = 0; k < blocks.easCount; ++k) x[k] = fa[k];
    if (!solveSmall(blocks.kaa, blocks.easCount, x, 1)) return;
    for (int j = 0; j < kDof; ++j) {
        double s = 0.0;
        for (int k = 0; k < blocks.easCount; ++k) s += blocks.kua[j * kEas + k] * x[k];
        out[j] -= s;
    }
}

void elementThermalLoad(const double nodes[kDof], const StructuralMaterial& material,
                        Formulation form, const double eigenstrain[6], double out[kDof]) {
    RestForms forms;
    computeForms(nodes, form, forms);
    elementThermalLoad(forms, material, eigenstrain, out);
}

void elementMass(const double nodes[kDof], double density, double out[kNodes]) {
    std::fill(out, out + kNodes, 0.0);
    RestForms forms;
    Shape shape;
    computeForms(nodes, Formulation::Displacement, forms, &shape);
    if (!forms.ok) return;
    for (int gp = 0; gp < kGauss; ++gp)
        for (int a = 0; a < kNodes; ++a) out[a] += density * forms.weight[gp] * shape[gp][a];
}

void gaussVolumes(const double nodes[kDof], double out[kGauss]) {
    std::fill(out, out + kGauss, 0.0);
    RestForms forms;
    // The weights are geometry, so the formulation is irrelevant and the cheapest
    // one is used. `weight` already carries det J times the 2x2x2 Gauss weight of 1.
    computeForms(nodes, Formulation::Displacement, forms);
    if (!forms.ok) return;
    for (int gp = 0; gp < kGauss; ++gp) out[gp] = forms.weight[gp];
}

bool computeRestForms(const double rest[kDof], Formulation form, RestForms& out) {
    computeForms(rest, form, out);
    return out.ok;
}

namespace {
// The half of `elementRotation` that is not the rest element's Jacobian: shared by
// the cached and uncached entry points so there is one polar decomposition here
// and not two that could drift apart.
void rotationFrom(const double restInv[3][3], const double current[kDof], double out[9]) {
    Point c;
    evaluate(current, 0.0, 0.0, 0.0, c);
    double f[3][3];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            double s = 0.0;
            for (int k = 0; k < 3; ++k) s += c.jac[i][k] * restInv[k][j];
            f[i][j] = s;
        }
    double rot[3][3];
    polarRotation(f, rot);
    // Column-major out, matching fem.cpp's M3.
    for (int col = 0; col < 3; ++col)
        for (int row = 0; row < 3; ++row) out[col * 3 + row] = rot[row][col];
}
}  // namespace

void elementRotation(const RestForms& forms, const double current[kDof], double out[9]) {
    if (!forms.ok) {
        for (int i = 0; i < 9; ++i) out[i] = (i % 4 == 0) ? 1.0 : 0.0;
        return;
    }
    double restInv[3][3];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) restInv[i][j] = forms.restJacobianInverse[i * 3 + j];
    rotationFrom(restInv, current, out);
}

void elementRotation(const double rest[kDof], const double current[kDof], double out[9]) {
    Point r;
    evaluate(rest, 0.0, 0.0, 0.0, r);
    double restInv[3][3];
    if (invert3(r.jac, restInv) == 0.0) {
        for (int i = 0; i < 9; ++i) out[i] = (i % 4 == 0) ? 1.0 : 0.0;
        return;
    }
    rotationFrom(restInv, current, out);
}

namespace {
// The co-rotational force off a rotation that has already been formed. Both
// `internalForce` overloads route through it, so there is one expression of
// `f = -R K (R^T x - X)` rather than two that could drift apart.
void internalForceFrom(const double r[9], const double stiffness[kDof * kDof],
                       const double rest[kDof], const double current[kDof], double out[kDof]) {
    // u = R^T x - X, node by node.
    double u[kDof];
    for (int a = 0; a < kNodes; ++a)
        for (int i = 0; i < 3; ++i) {
            double s = 0.0;
            for (int k = 0; k < 3; ++k) s += r[i * 3 + k] * current[a * 3 + k];  // R^T x
            u[a * 3 + i] = s - rest[a * 3 + i];
        }

    double ku[kDof];
    for (int i = 0; i < kDof; ++i) {
        double s = 0.0;
        for (int j = 0; j < kDof; ++j) s += stiffness[i * kDof + j] * u[j];
        ku[i] = s;
    }
    for (int a = 0; a < kNodes; ++a)
        for (int i = 0; i < 3; ++i) {
            double s = 0.0;
            for (int k = 0; k < 3; ++k) s += r[k * 3 + i] * ku[a * 3 + k];  // R (K u)
            out[a * 3 + i] = -s;
        }
}
}  // namespace

void internalForce(const double stiffness[kDof * kDof], const double rest[kDof],
                   const double current[kDof], double out[kDof]) {
    double r[9];
    elementRotation(rest, current, r);
    internalForceFrom(r, stiffness, rest, current, out);
}

void internalForce(const RestForms& forms, const double stiffness[kDof * kDof],
                   const double rest[kDof], const double current[kDof], double out[kDof]) {
    double r[9];
    elementRotation(forms, current, r);
    internalForceFrom(r, stiffness, rest, current, out);
}

// --- Plasticity and tearing ---------------------------------------------------

void elementSize(const double nodes[kDof], double* inPlane, double* thickness) {
    if (inPlane != nullptr) *inPlane = 0.0;
    if (thickness != nullptr) *thickness = 0.0;

    // Mid-surface quad, as two triangles, so a warped element is still exact --
    // the same construction `PlatePanel::area()` uses.
    double mid[4][3];
    for (int a = 0; a < 4; ++a)
        for (int i = 0; i < 3; ++i)
            mid[a][i] = 0.5 * (nodes[a * 3 + i] + nodes[(a + 4) * 3 + i]);

    double area = 0.0;
    for (int t = 0; t < 2; ++t) {
        const int b = t + 1, c = t + 2;
        double e1[3], e2[3];
        for (int i = 0; i < 3; ++i) {
            e1[i] = mid[b][i] - mid[0][i];
            e2[i] = mid[c][i] - mid[0][i];
        }
        const double cross[3] = {e1[1] * e2[2] - e1[2] * e2[1], e1[2] * e2[0] - e1[0] * e2[2],
                                 e1[0] * e2[1] - e1[1] * e2[0]};
        area += 0.5 * std::sqrt(cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2]);
    }
    if (!(area > 0.0)) return;

    RestForms forms;
    computeForms(nodes, Formulation::Displacement, forms);
    if (!forms.ok) return;
    double volume = 0.0;
    for (int gp = 0; gp < kGauss; ++gp) volume += forms.weight[gp];

    if (inPlane != nullptr) *inPlane = std::sqrt(area);
    // Volume over mid-surface area, not the length of an edge: a sheared element
    // has the thickness of the perpendicular between its faces, and its slanted
    // edge is longer than that by the shear.
    if (thickness != nullptr) *thickness = volume / area;
}

void initialisePlasticState(const double nodes[kDof], const plasticity::Material& material,
                            ElementPlasticState& state) {
    state = ElementPlasticState{};
    double inPlane = 0.0, thickness = 0.0;
    elementSize(nodes, &inPlane, &thickness);
    state.failureStrain =
        plasticity::regularisedFailureStrain(material.failure, inPlane, thickness);
}

PlasticUpdate elementPlasticUpdate(const double rest[kDof], const double current[kDof],
                                   const plasticity::Material& material, Formulation form,
                                   ElementPlasticState& state, double force[kDof],
                                   double stress[kGauss * 6], const double eigenstrain[6]) {
    RestForms forms;
    computeForms(rest, form, forms);
    return elementPlasticUpdate(forms, rest, current, material, state, force, stress, eigenstrain);
}

PlasticUpdate elementPlasticUpdate(const RestForms& forms, const double rest[kDof],
                                   const double current[kDof],
                                   const plasticity::Material& material,
                                   ElementPlasticState& state, double force[kDof],
                                   double stress[kGauss * 6], const double eigenstrain[6]) {
    PlasticUpdate result;
    std::fill(force, force + kDof, 0.0);
    if (stress != nullptr) std::fill(stress, stress + kGauss * 6, 0.0);

    if (!forms.ok) return result;

    // Co-rotated displacement: the plastic history lives in the material frame, so
    // a finite rotation must not touch it. u = R^T x - X, as `internalForce`.
    double r[9];
    elementRotation(forms, current, r);
    double u[kDof];
    for (int a = 0; a < kNodes; ++a)
        for (int i = 0; i < 3; ++i) {
            double s = 0.0;
            for (int k = 0; k < 3; ++k) s += r[i * 3 + k] * current[a * 3 + k];
            u[a * 3 + i] = s - rest[a * 3 + i];
        }

    // **Once any integration point has torn, the enhanced strains are dropped and
    // the element finishes its life as the ANS hex.**
    //
    // r(alpha) = int G^T sigma dV = 0 says the enhanced modes carry no stress over
    // the element. That is a statement about a continuum, and an element with a
    // dead integration point in it is not one: Kaa loses rank as the tangent at the
    // dead points goes to zero, the Newton stops converging, and alpha wanders.
    // Measured on a plate torn under an in-plane strain gradient, the consequence
    // was not a wobble but a **stall** -- with four of eight points gone, the
    // surviving four were driven to a triaxiality below the damage cutoff, their
    // damage froze at 0.78 while their plastic strain went on from 0.49 to 0.89,
    // and the element never finished tearing at all. Dropping the enhanced modes
    // costs a small jump in stress at the step a point dies, on an element that is
    // in the act of failing, and it is the difference between a tear that completes
    // and one that does not.
    int easCount = forms.easCount;
    bool degraded = false;
    for (int gp = 0; gp < kGauss; ++gp) degraded = degraded || state.point[gp].failed;
    if (degraded) easCount = 0;

    double volume = 0.0;
    for (int gp = 0; gp < kGauss; ++gp) volume += forms.weight[gp];

    // Every Newton iterate restarts from the history at the *start* of the step;
    // updating in place would make the answer depend on how many iterations it
    // took, which is the classic way an element-level solve stops being a solve.
    plasticity::State start[kGauss];
    for (int gp = 0; gp < kGauss; ++gp) start[gp] = state.point[gp];

    double alpha[kEas] = {};
    for (int k = 0; k < easCount; ++k) alpha[k] = state.enhanced[k];

    plasticity::State trial[kGauss];
    plasticity::Increment increments[kGauss];
    double gaussStress[kGauss][6];
    double tangent[kGauss][36];

    // The residual has units of stress times volume. Yield strength times element
    // volume is the natural scale, and it is a property of the problem rather than
    // of the iterate, so the tolerance does not move as the element unloads.
    const double scale = material.flow.yieldStrength * volume;
    constexpr int kMaxIterations = 40;

    // Two attempts at most. The second is entered only when an integration point
    // died *inside* this step, which is precisely the step in which Kaa was losing
    // rank while the Newton was solving on it -- so whatever it converged to is not
    // to be committed. Redoing the element without the enhanced modes costs one
    // extra pass, once in the element's life.
    for (int attempt = 0; attempt < 2; ++attempt) {
        result.converged = false;
        result.enhancedResidual = 0.0;
        result.enhancedWork = 0.0;

        for (int iteration = 0;; ++iteration) {
            for (int gp = 0; gp < kGauss; ++gp) {
                double strain[6];
                for (int i = 0; i < 6; ++i) {
                    double s = 0.0;
                    for (int j = 0; j < kDof; ++j) s += forms.b[gp][i][j] * u[j];
                    for (int k = 0; k < easCount; ++k) s += forms.g[gp][i][k] * alpha[k];
                    // The one subtraction. `plasticity::update` integrates the law
                    // from the stored history to the strain it is given, so what it
                    // must be given is the *mechanical* strain -- and the enhanced
                    // Newton below then converges on the alpha that makes the real
                    // stress orthogonal to G, with no separate eigenstrain term of
                    // its own.
                    strain[i] = eigenstrain != nullptr ? s - eigenstrain[i] : s;
                }
                trial[gp] = start[gp];
                increments[gp] = plasticity::update(material, state.failureStrain, strain, trial[gp],
                                                    gaussStress[gp], tangent[gp]);
            }
            result.iterations = iteration + 1;

            if (easCount == 0) {
                result.converged = true;
                break;
            }

            double residual[kEas] = {};
            for (int gp = 0; gp < kGauss; ++gp)
                for (int k = 0; k < easCount; ++k) {
                    double s = 0.0;
                    for (int i = 0; i < 6; ++i) s += forms.g[gp][i][k] * gaussStress[gp][i];
                    residual[k] += forms.weight[gp] * s;
                }
            double norm = 0.0;
            for (int k = 0; k < easCount; ++k) norm += residual[k] * residual[k];
            norm = std::sqrt(norm);
            result.enhancedResidual = norm;
            if (norm <= 1e-12 * scale) {
                result.converged = true;
                break;
            }
            if (iteration + 1 >= kMaxIterations) break;

            double kaa[kEas * kEas] = {};
            for (int gp = 0; gp < kGauss; ++gp) {
                double cg[6][kEas];
                for (int i = 0; i < 6; ++i)
                    for (int k = 0; k < easCount; ++k) {
                        double s = 0.0;
                        for (int m = 0; m < 6; ++m) s += tangent[gp][i * 6 + m] * forms.g[gp][m][k];
                        cg[i][k] = s;
                    }
                for (int p = 0; p < easCount; ++p)
                    for (int q = 0; q < easCount; ++q) {
                        double s = 0.0;
                        for (int i = 0; i < 6; ++i) s += forms.g[gp][i][p] * cg[i][q];
                        kaa[p * kEas + q] += forms.weight[gp] * s;
                    }
            }
            double delta[kEas];
            for (int k = 0; k < easCount; ++k) delta[k] = -residual[k];
            // A non-positive pivot here means every integration point has torn and the
            // element has no stiffness left, which is a state to report rather than to
            // iterate on.
            if (!solveSmall(kaa, easCount, delta, 1)) break;

            // Converge on the correction's **energy**, not on the residual's magnitude.
            //
            // ||r|| is not a scale-free measure here and it is worth being explicit
            // about why. The enhanced thickness modes carry E_zeta,zeta, so their
            // columns of G are scaled by the Voigt transform's 1/t^2 -- of order 3e7
            // for 20 mm plate -- and Kaa inherits the square of that. Measured on a
            // bent element: a residual of 1e-3 corresponds to an error in alpha of
            // 6e-18 against an alpha of 2.7e-6, twelve significant digits. Chasing
            // ||r|| below its own floor there costs forty iterations and moves nothing;
            // the 40-iteration answer and the 4-iteration answer agree to every digit
            // printed.
            //
            // delta . r is the work the correction would do, in joules, and sigma_y * V
            // is the element's yield energy, so the ratio is dimensionless and
            // independent of how the enhanced modes happen to be normalised. Stopping
            // *before* applying the correction keeps the committed history matching the
            // alpha it was computed at.
            double work = 0.0;
            for (int k = 0; k < easCount; ++k) work += delta[k] * residual[k];
            result.enhancedWork = std::abs(work);
            if (result.enhancedWork <= 1e-16 * scale) {
                result.converged = true;
                break;
            }
            for (int k = 0; k < easCount; ++k) alpha[k] += delta[k];
        }

        bool diedThisStep = false;
        for (int gp = 0; gp < kGauss; ++gp) diedThisStep = diedThisStep || trial[gp].failed;
        if (!(diedThisStep && easCount > 0)) break;
        easCount = 0;
        for (int k = 0; k < kEas; ++k) alpha[k] = 0.0;
    }

    for (int gp = 0; gp < kGauss; ++gp) {
        state.point[gp] = trial[gp];
        result.dissipation += forms.weight[gp] * increments[gp].dissipation;
        if (increments[gp].yielded) ++result.yieldedPoints;
        if (trial[gp].failed) ++result.failedPoints;
        if (stress != nullptr)
            for (int i = 0; i < 6; ++i) stress[gp * 6 + i] = gaussStress[gp][i];
    }
    // All seven, not just easCount of them: a degraded element must report zeros
    // rather than the values it had before it started tearing.
    for (int k = 0; k < kEas; ++k) state.enhanced[k] = alpha[k];
    state.torn = result.failedPoints == kGauss;

    double internal[kDof] = {};
    for (int gp = 0; gp < kGauss; ++gp)
        for (int j = 0; j < kDof; ++j) {
            double s = 0.0;
            for (int i = 0; i < 6; ++i) s += forms.b[gp][i][j] * gaussStress[gp][i];
            internal[j] += forms.weight[gp] * s;
        }
    for (int a = 0; a < kNodes; ++a)
        for (int i = 0; i < 3; ++i) {
            double s = 0.0;
            for (int k = 0; k < 3; ++k) s += r[k * 3 + i] * internal[a * 3 + k];
            force[a * 3 + i] = -s;
        }
    return result;
}

double smallestJacobian(const double nodes[kDof]) {
    double worst = std::numeric_limits<double>::infinity();
    for (int a = 0; a < kNodes; ++a) {
        Point p;
        evaluate(nodes, kXi[a], kEta[a], kZta[a], p);
        worst = std::min(worst, determinant3(p.jac));
    }
    return worst;
}

ElementShape elementShape(const double nodes[kDof]) {
    ElementShape out;

    // Exact coincidence, not a tolerance. A collapsed element's two nodes are the
    // *same mesh node* scattered twice, so their coordinates are bit-identical; a
    // sliver whose corners are merely close still has a small positive determinant
    // there and is still judged on it. A tolerance here would turn "this element
    // is a wedge" into "this element is nearly a wedge", which is the loosening
    // this whole classification exists to avoid.
    for (int a = 0; a < kNodes; ++a)
        for (int b = a + 1; b < kNodes; ++b) {
            const double dx = nodes[a * 3 + 0] - nodes[b * 3 + 0];
            const double dy = nodes[a * 3 + 1] - nodes[b * 3 + 1];
            const double dz = nodes[a * 3 + 2] - nodes[b * 3 + 2];
            if (dx == 0.0 && dy == 0.0 && dz == 0.0) out.collapsed[a] = out.collapsed[b] = true;
        }
    for (int a = 0; a < kNodes; ++a)
        if (out.collapsed[a]) ++out.collapsedNodes;

    out.smallestNodal = std::numeric_limits<double>::infinity();
    bool soundCorners = true;
    for (int a = 0; a < kNodes; ++a) {
        Point p;
        evaluate(nodes, kXi[a], kEta[a], kZta[a], p);
        const double det = determinant3(p.jac);
        out.smallestNodal = std::min(out.smallestNodal, det);
        if (!(det > 0.0) && !out.collapsed[a]) soundCorners = false;
    }

    // The centre and the 2x2x2 rule: exactly where `computeForms` evaluates, and
    // exactly what it requires to be positive. Asking anywhere else would be
    // asking a different question from the one the solve asks.
    const double q = 1.0 / std::sqrt(3.0);
    out.smallestGauss = std::numeric_limits<double>::infinity();
    for (int gp = -1; gp < kGauss; ++gp) {
        const double xi = gp < 0 ? 0.0 : ((gp & 1) ? q : -q);
        const double eta = gp < 0 ? 0.0 : ((gp & 2) ? q : -q);
        const double zta = gp < 0 ? 0.0 : ((gp & 4) ? q : -q);
        Point p;
        evaluate(nodes, xi, eta, zta, p);
        out.smallestGauss = std::min(out.smallestGauss, determinant3(p.jac));
    }

    out.integrable = soundCorners && out.smallestGauss > 0.0;
    return out;
}

double criticalTimestep(const double nodes[kDof], const StructuralMaterial& material,
                        Formulation form, double safety) {
    double k[kDof * kDof];
    elementStiffness(nodes, material, form, k);
    double m[kNodes];
    elementMass(nodes, material.density, m);
    for (int a = 0; a < kNodes; ++a)
        if (!(m[a] > 0.0)) return 0.0;

    // Power iteration on the symmetric M^-1/2 K M^-1/2, whose largest eigenvalue
    // is omega_max^2. Symmetrising first matters: the unsymmetric M^-1 K has the
    // same spectrum but power iteration on it converges on the wrong vector norm.
    double d[kDof];
    for (int a = 0; a < kNodes; ++a)
        for (int i = 0; i < 3; ++i) d[a * 3 + i] = 1.0 / std::sqrt(m[a]);

    double v[kDof];
    for (int i = 0; i < kDof; ++i) v[i] = 1.0 + 0.37 * i - 0.11 * i * i;  // no symmetry to be trapped by
    double lambda = 0.0;
    for (int iteration = 0; iteration < 400; ++iteration) {
        double w[kDof];
        for (int i = 0; i < kDof; ++i) {
            double s = 0.0;
            for (int j = 0; j < kDof; ++j) s += k[i * kDof + j] * (d[j] * v[j]);
            w[i] = d[i] * s;
        }
        double numerator = 0.0, denominator = 0.0;
        for (int i = 0; i < kDof; ++i) {
            numerator += v[i] * w[i];
            denominator += v[i] * v[i];
        }
        const double next = denominator > 0.0 ? numerator / denominator : 0.0;
        double norm = 0.0;
        for (int i = 0; i < kDof; ++i) norm += w[i] * w[i];
        norm = std::sqrt(norm);
        if (!(norm > 0.0)) return 0.0;
        for (int i = 0; i < kDof; ++i) v[i] = w[i] / norm;
        if (iteration > 8 && std::abs(next - lambda) <= 1e-10 * std::abs(next)) {
            lambda = next;
            break;
        }
        lambda = next;
    }
    if (!(lambda > 0.0)) return 0.0;
    return safety * 2.0 / std::sqrt(lambda);
}

// --- Meshes -------------------------------------------------------------------

void HexMesh::pin(std::size_t node, int axis, double value) {
    if (fixed.size() != nodeCount() * 3) fixed.assign(nodeCount() * 3, 0u);
    if (prescribed.size() != nodeCount() * 3) prescribed.assign(nodeCount() * 3, 0.0);
    fixed[node * 3 + static_cast<std::size_t>(axis)] = 1u;
    prescribed[node * 3 + static_cast<std::size_t>(axis)] = value;
}

void HexMesh::gather(std::size_t element, const std::vector<double>& nodal,
                     double out[kDof]) const {
    for (int a = 0; a < kNodes; ++a) {
        const std::size_t n = index[element * kNodes + static_cast<std::size_t>(a)];
        for (int i = 0; i < 3; ++i) out[a * 3 + i] = nodal[n * 3 + static_cast<std::size_t>(i)];
    }
}

HexMesh makePlateMesh(double lx, double ly, double thickness, int nx, int ny, int nz) {
    HexMesh mesh;
    const int sx = nx + 1, sy = ny + 1, sz = nz + 1;
    const auto at = [&](int i, int j, int k) {
        return static_cast<std::uint32_t>((i * sy + j) * sz + k);
    };

    mesh.position.resize(static_cast<std::size_t>(sx) * sy * sz * 3);
    for (int i = 0; i < sx; ++i)
        for (int j = 0; j < sy; ++j)
            for (int k = 0; k < sz; ++k) {
                double* p = &mesh.position[at(i, j, k) * 3];
                p[0] = lx * static_cast<double>(i) / nx;
                p[1] = ly * static_cast<double>(j) / ny;
                // Mid-surface at z = 0, the convention scantlings.hpp fixes for
                // plate panels, so bending is symmetric about the element centre.
                p[2] = -0.5 * thickness + thickness * static_cast<double>(k) / nz;
            }

    for (int i = 0; i < nx; ++i)
        for (int j = 0; j < ny; ++j)
            for (int k = 0; k < nz; ++k) {
                const std::uint32_t v[kNodes] = {
                    at(i, j, k),             at(i + 1, j, k),
                    at(i + 1, j + 1, k),     at(i, j + 1, k),
                    at(i, j, k + 1),         at(i + 1, j, k + 1),
                    at(i + 1, j + 1, k + 1), at(i, j + 1, k + 1)};
                mesh.index.insert(mesh.index.end(), v, v + kNodes);
            }

    mesh.fixed.assign(mesh.nodeCount() * 3, 0u);
    mesh.prescribed.assign(mesh.nodeCount() * 3, 0.0);
    return mesh;
}

double criticalTimestep(const HexMesh& mesh, const StructuralMaterial& material,
                        Formulation form, double safety) {
    double smallest = std::numeric_limits<double>::infinity();
    std::vector<double> nodal(mesh.position);
    for (std::size_t e = 0; e < mesh.elementCount(); ++e) {
        double nodes[kDof];
        mesh.gather(e, nodal, nodes);
        const double dt = criticalTimestep(nodes, material, form, safety);
        if (dt > 0.0) smallest = std::min(smallest, dt);
    }
    return std::isfinite(smallest) ? smallest : 0.0;
}

// --- Banded solver ------------------------------------------------------------

BandedSpd::BandedSpd(std::size_t dofCount, std::size_t halfBandwidth)
    : n_(dofCount), b_(halfBandwidth), a_(dofCount * (halfBandwidth + 1), 0.0) {}

void BandedSpd::add(std::size_t row, std::size_t column, double value) {
    if (row < column) return;  // lower triangle only; see the header
    const std::size_t offset = row - column;
    if (offset > b_) return;
    a_[row * (b_ + 1) + offset] += value;
}

bool BandedSpd::factor() {
    for (std::size_t i = 0; i < n_; ++i) {
        const std::size_t lo = i >= b_ ? i - b_ : 0;
        for (std::size_t j = lo; j <= i; ++j) {
            double s = a_[i * (b_ + 1) + (i - j)];
            for (std::size_t k = lo; k < j; ++k)
                s -= a_[i * (b_ + 1) + (i - k)] * a_[j * (b_ + 1) + (j - k)];
            if (j < i) {
                a_[i * (b_ + 1) + (i - j)] = s / a_[j * (b_ + 1)];
            } else {
                if (!(s > 0.0)) return false;
                a_[i * (b_ + 1)] = std::sqrt(s);
            }
        }
    }
    return true;
}

void BandedSpd::solve(std::vector<double>& rightHandSide) const {
    for (std::size_t i = 0; i < n_; ++i) {
        const std::size_t lo = i >= b_ ? i - b_ : 0;
        double s = rightHandSide[i];
        for (std::size_t k = lo; k < i; ++k) s -= a_[i * (b_ + 1) + (i - k)] * rightHandSide[k];
        rightHandSide[i] = s / a_[i * (b_ + 1)];
    }
    for (std::size_t ii = n_; ii > 0; --ii) {
        const std::size_t i = ii - 1;
        double s = rightHandSide[i];
        const std::size_t hi = std::min(n_, i + b_ + 1);
        for (std::size_t k = i + 1; k < hi; ++k) s -= a_[k * (b_ + 1) + (k - i)] * rightHandSide[k];
        rightHandSide[i] = s / a_[i * (b_ + 1)];
    }
}

DofExpansion::DofExpansion(std::size_t dofCount, const std::vector<Mpc>& constrained) {
    isSlave_.assign(dofCount, 0u);
    std::vector<std::size_t> width(dofCount, 1);
    for (std::size_t c = 0; c < constrained.size(); ++c) {
        const Mpc& mpc = constrained[c];
        const std::string where = "constraint " + std::to_string(c);
        if (mpc.slave >= dofCount) {
            ok_ = false;
            problem_ = where + " eliminates degree of freedom " + std::to_string(mpc.slave) +
                       ", which is not in a system of " + std::to_string(dofCount);
            return;
        }
        if (mpc.master.size() != mpc.weight.size() || mpc.master.empty()) {
            ok_ = false;
            problem_ = where + " has " + std::to_string(mpc.master.size()) + " masters and " +
                       std::to_string(mpc.weight.size()) + " weights";
            return;
        }
        if (isSlave_[mpc.slave]) {
            ok_ = false;
            problem_ = where + " eliminates degree of freedom " + std::to_string(mpc.slave) +
                       " a second time; a degree of freedom has one expansion or none";
            return;
        }
        for (std::uint32_t m : mpc.master) {
            if (m >= dofCount) {
                ok_ = false;
                problem_ = where + " names master degree of freedom " + std::to_string(m) +
                           ", which is not in a system of " + std::to_string(dofCount);
                return;
            }
            if (m == mpc.slave) {
                ok_ = false;
                problem_ = where + " names its own slave " + std::to_string(m) + " as a master";
                return;
            }
        }
        isSlave_[mpc.slave] = 1u;
        width[mpc.slave] = mpc.master.size();
        ++eliminated_;
    }
    // A second pass, because "is this master a slave" is only answerable once every
    // slave is known -- checking it as the constraints arrive would accept a chain
    // written in one order and refuse the same chain written in the other.
    for (const Mpc& mpc : constrained)
        for (std::uint32_t m : mpc.master)
            if (isSlave_[m]) {
                ok_ = false;
                problem_ = "degree of freedom " + std::to_string(m) +
                           " is both a master and a slave: chained constraints are refused rather"
                           " than composed, because a chain resolved silently is a modelling error"
                           " that assembles";
                return;
            }

    start_.assign(dofCount + 1, 0);
    for (std::size_t d = 0; d < dofCount; ++d) start_[d + 1] = start_[d] + width[d];
    term_.assign(start_[dofCount], Term{0, 0.0});
    for (std::size_t d = 0; d < dofCount; ++d)
        if (!isSlave_[d]) term_[start_[d]] = Term{static_cast<std::uint32_t>(d), 1.0};
    for (const Mpc& mpc : constrained) {
        std::size_t at = start_[mpc.slave];
        for (std::size_t a = 0; a < mpc.master.size(); ++a)
            term_[at++] = Term{mpc.master[a], mpc.weight[a]};
    }
}

void DofExpansion::recover(std::vector<double>& values) const {
    for (std::size_t d = 0; d < isSlave_.size(); ++d) {
        if (!isSlave_[d]) continue;
        double sum = 0;
        for (const Term* t = begin(static_cast<std::uint32_t>(d));
             t != end(static_cast<std::uint32_t>(d)); ++t)
            sum += t->weight * values[t->dof];
        values[d] = sum;
    }
}

bool solveStatic(const HexMesh& mesh, const StructuralMaterial& material, Formulation form,
                 const std::vector<double>& load, std::vector<double>& displacement,
                 std::string* problem) {
    static const std::vector<DofBlock> none;
    return solveStatic(mesh, material, form, none, load, displacement, problem);
}

bool solveStatic(const HexMesh& mesh, const StructuralMaterial& material, Formulation form,
                 const std::vector<DofBlock>& extra, const std::vector<double>& load,
                 std::vector<double>& displacement, std::string* problem) {
    static const std::vector<Mpc> none;
    return solveStatic(mesh, material, form, extra, none, load, displacement, problem);
}

bool solveStatic(const HexMesh& mesh, const StructuralMaterial& material, Formulation form,
                 const std::vector<DofBlock>& extra, const std::vector<Mpc>& constrained,
                 const std::vector<double>& load, std::vector<double>& displacement,
                 std::string* problem) {
    const std::size_t nodes = mesh.nodeCount();
    const std::size_t elements = mesh.elementCount();

    // Geometry first, before anything that could return early. A mesh with every
    // degree of freedom prescribed has nothing to solve, but if one of its
    // elements is inside out the caller has a broken mesh and needs to hear so --
    // reporting success because there happened to be no work is the failing-open
    // pattern this repo has been bitten by before.
    for (std::size_t e = 0; e < elements; ++e) {
        double nodePos[kDof];
        mesh.gather(e, mesh.position, nodePos);
        // `integrable`, not `smallestJacobian > 0`: a collapsed hexahedron is a
        // wedge and integrates, an inverted one does not. See `ElementShape`.
        if (!elementShape(nodePos).integrable) {
            if (problem) *problem = "element " + std::to_string(e) + " is inverted or degenerate";
            return false;
        }
    }

    const DofExpansion expansion(nodes * 3, constrained);
    if (!expansion.ok()) {
        if (problem) *problem = expansion.problem();
        return false;
    }
    for (std::size_t d = 0; d < nodes * 3; ++d)
        if (expansion.eliminated(static_cast<std::uint32_t>(d)) && d < mesh.fixed.size() &&
            mesh.fixed[d]) {
            if (problem)
                *problem = "degree of freedom " + std::to_string(d) +
                           " is both prescribed by the mesh and eliminated by a constraint";
            return false;
        }

    displacement.assign(nodes * 3, 0.0);
    for (std::size_t d = 0; d < nodes * 3; ++d)
        if (d < mesh.fixed.size() && mesh.fixed[d]) displacement[d] = mesh.prescribed[d];

    // Free DOF numbering, and the bandwidth that numbering actually delivers. An
    // eliminated degree of freedom is not free and not prescribed: it has no
    // unknown at all, and its row is carried by its masters.
    std::vector<std::ptrdiff_t> map(nodes * 3, -1);
    std::size_t free = 0;
    for (std::size_t d = 0; d < nodes * 3; ++d)
        if (!(d < mesh.fixed.size() && mesh.fixed[d]) &&
            !expansion.eliminated(static_cast<std::uint32_t>(d)))
            map[d] = static_cast<std::ptrdiff_t>(free++);
    if (free == 0) return true;

    // The reach of one global degree of freedom, after expansion: the free slots
    // its row and column actually touch. A constrained one reaches its masters and
    // not itself, which is what makes the band below cover the transformed system
    // rather than the untransformed one.
    const auto reach = [&](std::size_t global, std::size_t& lo, std::size_t& hi, bool& any) {
        for (const DofExpansion::Term* t = expansion.begin(static_cast<std::uint32_t>(global));
             t != expansion.end(static_cast<std::uint32_t>(global)); ++t) {
            const std::ptrdiff_t d = map[t->dof];
            if (d < 0) continue;
            const auto u = static_cast<std::size_t>(d);
            lo = std::min(lo, u);
            hi = std::max(hi, u);
            any = true;
        }
    };

    std::size_t band = 0;
    for (std::size_t e = 0; e < elements; ++e) {
        std::size_t lo = free, hi = 0;
        bool any = false;
        for (int a = 0; a < kNodes; ++a) {
            const std::size_t n = mesh.index[e * kNodes + static_cast<std::size_t>(a)];
            for (int i = 0; i < 3; ++i) reach(n * 3 + static_cast<std::size_t>(i), lo, hi, any);
        }
        if (any) band = std::max(band, hi - lo);
    }
    // The extra blocks widen the band on the same rule. Leaving them out would
    // silently drop every entry that fell outside it, which reads as a slightly
    // soft answer rather than as a missing term.
    for (const DofBlock& block : extra) {
        // Skipped on the same rule the assembly loop uses, or a block that will
        // never be scattered still widens the band and is paid for in `n*(b+1)`
        // doubles and an `n*b^2` factorisation. Harmless to the answer -- assembly
        // coverage is a strict subset of this -- and pure waste.
        if (block.stiffness.size() < block.dof.size() * block.dof.size()) continue;
        std::size_t lo = free, hi = 0;
        bool any = false;
        for (std::uint32_t global : block.dof) {
            if (global >= nodes * 3) continue;
            reach(global, lo, hi, any);
        }
        if (any) band = std::max(band, hi - lo);
    }

    BandedSpd system(free, band);
    std::vector<double> rhs(free, 0.0);
    // A load applied at an eliminated degree of freedom is distributed to its
    // masters by the transpose of the constraint, which is what keeps the virtual
    // work identical -- the same statement `constraint::applyTiedForce` makes.
    for (std::size_t d = 0; d < nodes * 3; ++d) {
        const double f = d < load.size() ? load[d] : 0.0;
        if (f == 0.0) continue;
        for (const DofExpansion::Term* t = expansion.begin(static_cast<std::uint32_t>(d));
             t != expansion.end(static_cast<std::uint32_t>(d)); ++t)
            if (map[t->dof] >= 0) rhs[static_cast<std::size_t>(map[t->dof])] += t->weight * f;
    }

    // One scatter for both sources, so an element and a block cannot disagree about
    // what a constrained degree of freedom means.
    const auto scatter = [&](std::size_t row, std::size_t column, double value) {
        for (const DofExpansion::Term* p = expansion.begin(static_cast<std::uint32_t>(row));
             p != expansion.end(static_cast<std::uint32_t>(row)); ++p) {
            const std::ptrdiff_t rowDof = map[p->dof];
            if (rowDof < 0) continue;
            for (const DofExpansion::Term* q = expansion.begin(static_cast<std::uint32_t>(column));
                 q != expansion.end(static_cast<std::uint32_t>(column)); ++q) {
                const std::ptrdiff_t colDof = map[q->dof];
                const double term = p->weight * q->weight * value;
                if (colDof >= 0) {
                    system.add(static_cast<std::size_t>(rowDof), static_cast<std::size_t>(colDof),
                               term);
                } else {
                    // Prescribed DOF move to the right-hand side exactly, rather
                    // than through a penalty stiffness whose size would then be a
                    // tolerance the patch test had to live with.
                    rhs[static_cast<std::size_t>(rowDof)] -= term * displacement[q->dof];
                }
            }
        }
    };

    for (std::size_t e = 0; e < elements; ++e) {
        double nodePos[kDof];
        mesh.gather(e, mesh.position, nodePos);
        double ke[kDof * kDof];
        elementStiffness(nodePos, material, form, ke);

        std::size_t dof[kDof];
        for (int a = 0; a < kNodes; ++a) {
            const std::size_t n = mesh.index[e * kNodes + static_cast<std::size_t>(a)];
            for (int i = 0; i < 3; ++i) dof[a * 3 + i] = n * 3 + static_cast<std::size_t>(i);
        }
        for (int p = 0; p < kDof; ++p)
            for (int qq = 0; qq < kDof; ++qq) scatter(dof[p], dof[qq], ke[p * kDof + qq]);
    }

    for (const DofBlock& block : extra) {
        const std::size_t n = block.dof.size();
        if (block.stiffness.size() < n * n) continue;
        for (std::size_t p = 0; p < n; ++p) {
            if (block.dof[p] >= nodes * 3) continue;
            for (std::size_t qq = 0; qq < n; ++qq) {
                if (block.dof[qq] >= nodes * 3) continue;
                scatter(block.dof[p], block.dof[qq], block.stiffness[p * n + qq]);
            }
        }
    }

    if (!system.factor()) {
        if (problem)
            *problem = "stiffness is not positive definite: the constraints leave a rigid body mode";
        return false;
    }
    system.solve(rhs);
    for (std::size_t d = 0; d < nodes * 3; ++d)
        if (map[d] >= 0) displacement[d] = rhs[static_cast<std::size_t>(map[d])];
    expansion.recover(displacement);
    return true;
}

std::vector<double> uniformPressureLoad(const HexMesh& mesh, double pressure) {
    std::vector<double> load(mesh.nodeCount() * 3, 0.0);
    static constexpr int kFaces[6][4] = {{0, 1, 2, 3}, {4, 5, 6, 7}, {0, 1, 5, 4},
                                         {3, 2, 6, 7}, {0, 3, 7, 4}, {1, 2, 6, 5}};

    // A face carried by two elements is interior. Keyed on the sorted node set, so
    // the two elements' opposite windings still match.
    std::map<std::array<std::uint32_t, 4>, int> shared;
    for (std::size_t e = 0; e < mesh.elementCount(); ++e)
        for (const auto& face : kFaces) {
            std::array<std::uint32_t, 4> key{};
            for (int i = 0; i < 4; ++i)
                key[static_cast<std::size_t>(i)] =
                    mesh.index[e * kNodes + static_cast<std::size_t>(face[i])];
            std::sort(key.begin(), key.end());
            ++shared[key];
        }

    const double q = 1.0 / std::sqrt(3.0);
    static constexpr double kCorner[4][2] = {{-1, -1}, {+1, -1}, {+1, +1}, {-1, +1}};
    for (std::size_t e = 0; e < mesh.elementCount(); ++e) {
        std::uint32_t n[4];
        std::array<std::uint32_t, 4> key{};
        for (int i = 0; i < 4; ++i) {
            n[i] = mesh.index[e * kNodes + static_cast<std::size_t>(4 + i)];
            key[static_cast<std::size_t>(i)] = n[i];
        }
        std::sort(key.begin(), key.end());
        if (shared[key] != 1) continue;

        for (int gs = 0; gs < 2; ++gs)
            for (int gt = 0; gt < 2; ++gt) {
                const double s = gs ? q : -q, t = gt ? q : -q;
                double shape[4], ds[4], dt[4];
                for (int c = 0; c < 4; ++c) {
                    shape[c] = 0.25 * (1.0 + s * kCorner[c][0]) * (1.0 + t * kCorner[c][1]);
                    ds[c] = 0.25 * kCorner[c][0] * (1.0 + t * kCorner[c][1]);
                    dt[c] = 0.25 * kCorner[c][1] * (1.0 + s * kCorner[c][0]);
                }
                double xs[3] = {0, 0, 0}, xt[3] = {0, 0, 0};
                for (int c = 0; c < 4; ++c)
                    for (int i = 0; i < 3; ++i) {
                        xs[i] += ds[c] * mesh.position[n[c] * 3 + static_cast<std::size_t>(i)];
                        xt[i] += dt[c] * mesh.position[n[c] * 3 + static_cast<std::size_t>(i)];
                    }
                // (xs x xt) ds dt is the outward area vector of the +zeta face;
                // pressure pushes the other way.
                const double area[3] = {xs[1] * xt[2] - xs[2] * xt[1],
                                        xs[2] * xt[0] - xs[0] * xt[2],
                                        xs[0] * xt[1] - xs[1] * xt[0]};
                for (int c = 0; c < 4; ++c)
                    for (int i = 0; i < 3; ++i)
                        load[n[c] * 3 + static_cast<std::size_t>(i)] -=
                            pressure * shape[c] * area[i];
            }
    }
    return load;
}

std::vector<double> thermalLoad(const HexMesh& mesh, const StructuralMaterial& material,
                                Formulation form, const std::vector<double>& eigenstrain) {
    std::vector<double> load;
    if (eigenstrain.size() != mesh.elementCount() * 6) return load;
    load.assign(mesh.nodeCount() * 3, 0.0);
    for (std::size_t e = 0; e < mesh.elementCount(); ++e) {
        double nodes[kDof];
        for (int a = 0; a < kNodes; ++a) {
            const std::size_t n = mesh.index[e * kNodes + static_cast<std::size_t>(a)];
            for (int i = 0; i < 3; ++i)
                nodes[a * 3 + i] = mesh.position[n * 3 + static_cast<std::size_t>(i)];
        }
        double f[kDof];
        elementThermalLoad(nodes, material, form, &eigenstrain[e * 6], f);
        for (int a = 0; a < kNodes; ++a) {
            const std::size_t n = mesh.index[e * kNodes + static_cast<std::size_t>(a)];
            for (int i = 0; i < 3; ++i)
                load[n * 3 + static_cast<std::size_t>(i)] += f[a * 3 + i];
        }
    }
    return load;
}

}  // namespace sim::solidshell
