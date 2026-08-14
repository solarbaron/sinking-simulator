// SPDX-License-Identifier: MIT
#include "radiation.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <string>

namespace sim {
namespace {

using Complex = std::complex<double>;

double sqr(double x) { return x * x; }

// --- exp(z) E1(z) ------------------------------------------------------------
//
// The free-surface Green's function below is an exponential-integral evaluation
// in disguise, and it is always wanted multiplied by exp(z). Forming the product
// directly is what keeps it finite: E1 grows like e^{-z}/z far into the left
// half-plane, so E1 alone overflows where the product is of order 1/z.
//
// Only Re(z) <= 0 is ever asked for -- z = nu (y_field + y_source) + i nu dx and
// both depths are negative. The branch cut on the negative real axis is
// approached from above, which is the convention the +i pi sgn(dx) term in
// waveIntegral() is written to cancel against.
Complex expTimesE1(Complex z) {
    const double magnitude = std::abs(z);
    if (magnitude < 1e-300) return {0.0, 0.0};

    // The power series is well conditioned either where |z| is small or where
    // Re(z) is strongly negative (there the answer is exponentially large and
    // the terms cannot cancel it away). The gap -- large |Im z|, small Re z --
    // is where the continued fraction is at its best, so the two cover the plane
    // between them.
    if (magnitude <= 16.0 || z.real() < -8.0) {
        constexpr double kEulerGamma = 0.57721566490153286061;
        Complex term{1.0, 0.0};
        Complex sum{0.0, 0.0};
        for (int n = 1; n <= 400; ++n) {
            term *= -z / static_cast<double>(n);          // (-z)^n / n!
            const Complex add = -term / static_cast<double>(n);
            sum += add;
            if (n > 4 && std::abs(add) < 1e-18 * (std::abs(sum) + 1e-300)) break;
        }
        return std::exp(z) * (-kEulerGamma - std::log(z) + sum);
    }

    // Modified Lentz evaluation of exp(z) E1(z) = 1/(z+1- 1^2/(z+3- 2^2/(z+5- ...
    constexpr double kTiny = 1e-300;
    Complex f{kTiny, 0.0};
    Complex c = f;
    Complex d{0.0, 0.0};
    for (int i = 1; i <= 800; ++i) {
        const double an = (i == 1) ? 1.0 : -sqr(static_cast<double>(i - 1));
        const Complex b = z + Complex(2.0 * i - 1.0, 0.0);
        d = b + an * d;
        if (std::abs(d) < kTiny) d = kTiny;
        c = b + an / c;
        if (std::abs(c) < kTiny) c = kTiny;
        d = 1.0 / d;
        const Complex delta = c * d;
        f *= delta;
        if (std::abs(delta - 1.0) < 1e-15) break;
    }
    return f;
}

// PV integral_0^inf e^{k sy} e^{i k dx} / (k - nu) dk, closed form. `sy` is the
// *sum* of the field and source depths and is negative inside the fluid, so the
// integral is singular only when both points sit on the free surface.
Complex waveIntegral(double dx, double sy, double nu) {
    const Complex z{nu * sy, nu * dx};
    const double sign = (dx >= 0.0) ? 1.0 : -1.0;
    return expTimesE1(z) + Complex(0.0, kPi) * std::exp(z) * sign;
}

// --- Flat constant-strength source panel -------------------------------------

struct Panel {
    double ax = 0, ay = 0, bx = 0, by = 0;   // endpoints
    double mx = 0, my = 0;                   // midpoint (collocation point)
    double tx = 1, ty = 0;                   // unit tangent
    double nx = 0, ny = -1;                  // unit normal, out of the body
    double length = 0;
};

Panel makePanel(double ax, double ay, double bx, double by) {
    Panel p;
    p.ax = ax; p.ay = ay; p.bx = bx; p.by = by;
    p.length = std::hypot(bx - ax, by - ay);
    if (p.length > 0) {
        p.tx = (bx - ax) / p.length;
        p.ty = (by - ay) / p.length;
    }
    // The contour is wound keel -> waterline on the starboard (x > 0) side, for
    // which (t_y, -t_x) points out of the body.
    p.nx = p.ty;
    p.ny = -p.tx;
    p.mx = 0.5 * (ax + bx);
    p.my = 0.5 * (ay + by);
    return p;
}

Panel mirrorAcrossCentreline(const Panel& p) { return makePanel(-p.ax, p.ay, -p.bx, p.by); }
Panel mirrorAcrossSurface(const Panel& p) { return makePanel(p.ax, -p.ay, p.bx, -p.by); }

// phi = integral_panel ln|Q - P| ds, and its gradient. `onPanel` forces the
// offset to exactly zero so that the self-influence picks up the +pi jump
// analytically instead of through a signed zero.
void logPanel(const Panel& p, double qx, double qy, bool onPanel, double& phi, double& vx,
              double& vy) {
    const double dx = qx - p.ax;
    const double dy = qy - p.ay;
    const double s = dx * p.tx + dy * p.ty;
    const double h = onPanel ? 0.0 : dx * p.nx + dy * p.ny;
    const double r1 = std::max(std::hypot(s, h), 1e-300);
    const double r2 = std::max(std::hypot(s - p.length, h), 1e-300);
    const double theta1 = std::atan2(h, s);
    const double theta2 = std::atan2(h, s - p.length);
    phi = (p.length - s) * std::log(r2) + s * std::log(r1) - p.length + h * (theta2 - theta1);
    const double alongPanel = std::log(r1 / r2);
    const double acrossPanel = theta2 - theta1;   // exactly pi at the panel's own midpoint
    vx = alongPanel * p.tx + acrossPanel * p.nx;
    vy = alongPanel * p.ty + acrossPanel * p.ny;
}

// Eight-point Gauss-Legendre on [-1, 1]; enough for the smooth wave part of the
// Green's function over one panel.
constexpr int kGaussCount = 8;
// Eight-point nodes first, then four-point: panelInfluence() picks a slice.
constexpr double kGaussNode[kGaussCount + kGaussCount / 2] = {
    -0.96028985649753623, -0.79666647741362674, -0.52553240991632899, -0.18343464249564980,
    0.18343464249564980,  0.52553240991632899,  0.79666647741362674,  0.96028985649753623,
    -0.86113631159405258, -0.33998104358485626, 0.33998104358485626,  0.86113631159405258};
constexpr double kGaussWeight[kGaussCount + kGaussCount / 2] = {
    0.10122853629037626, 0.22238103445337447, 0.31370664587788729, 0.36268378337836198,
    0.36268378337836198, 0.31370664587788729, 0.22238103445337447, 0.10122853629037626,
    0.34785484513745385, 0.65214515486254614, 0.65214515486254614,  0.34785484513745385};

struct Influence {
    Complex potential;
    Complex gradientX;
    Complex gradientY;
};

// The deep-water free-surface source, integrated over one panel:
//
//     G = ln r - ln r' - 2 Re{J} + 2 pi i e^{nu (y + eta)} cos(nu (x - xi))
//
// with r' the distance to the source's mirror image in the free surface and J
// the principal-value integral above. Im{J} is *not* a solution of the
// free-surface condition -- only its real part is -- so the outgoing behaviour
// has to come from the regular standing wave, and getting that wrong produces a
// purely real source distribution and identically zero damping, which is exactly
// as wrong as it is quiet.
//
// `rigidLid` is the omega -> infinity limit, where the condition degenerates to
// phi = 0 on y = 0 and the wave terms vanish identically.
Influence panelInfluence(const Panel& p, double qx, double qy, double nu, bool onPanel,
                         bool rigidLid) {
    // The wave part varies on the shorter of the radiated wavelength and the
    // distance to the panel, so a well-separated panel in long waves needs half
    // the quadrature of a near one. This is the whole inner loop of the method
    // and the halving is worth having.
    const bool nearField =
        onPanel || std::hypot(qx - p.mx, qy - p.my) < 3.0 * p.length || nu * p.length > 0.5;
    const int nodes = nearField ? kGaussCount : kGaussCount / 2;
    const int first = nearField ? 0 : kGaussCount;
    double phi = 0, vx = 0, vy = 0;
    logPanel(p, qx, qy, onPanel, phi, vx, vy);
    const Panel image = mirrorAcrossSurface(p);
    double phiImage = 0, vxImage = 0, vyImage = 0;
    logPanel(image, qx, qy, false, phiImage, vxImage, vyImage);

    Influence out;
    out.potential = Complex(phi - phiImage, 0.0);
    out.gradientX = Complex(vx - vxImage, 0.0);
    out.gradientY = Complex(vy - vyImage, 0.0);
    if (rigidLid) return out;

    Complex wave{0, 0}, waveX{0, 0}, waveY{0, 0};
    for (int g = first; g < first + nodes; ++g) {
        const double u = 0.5 * (kGaussNode[g] + 1.0) * p.length;
        const double xi = p.ax + u * p.tx;
        const double eta = p.ay + u * p.ty;
        const double dx = qx - xi;
        const double sy = qy + eta;
        const Complex j0 = waveIntegral(dx, sy, nu);
        // d/dx J = i J1, d/dy J = J1, with J1 = 1/w + nu J and w = -sy - i dx.
        const Complex j1 = 1.0 / Complex(-sy, -dx) + nu * j0;
        const double weight = kGaussWeight[g] * 0.5 * p.length;
        const double decay = std::exp(nu * sy);
        wave += weight * (Complex(-2.0 * j0.real(), 0.0) +
                          Complex(0.0, 2.0 * kPi) * decay * std::cos(nu * dx));
        waveX += weight * (Complex(-2.0 * (Complex(0, 1) * j1).real(), 0.0) +
                           Complex(0.0, 2.0 * kPi) * decay * (-nu * std::sin(nu * dx)));
        waveY += weight * (Complex(-2.0 * j1.real(), 0.0) +
                           Complex(0.0, 2.0 * kPi) * decay * nu * std::cos(nu * dx));
    }
    out.potential += wave;
    out.gradientX += waveX;
    out.gradientY += waveY;
    return out;
}

// --- Small dense linear algebra ----------------------------------------------

// In-place LU with partial pivoting; `pivot` records the row swaps.
bool luFactor(std::vector<Complex>& a, int n, std::vector<int>& pivot) {
    pivot.resize(static_cast<std::size_t>(n));
    for (int k = 0; k < n; ++k) {
        int best = k;
        for (int i = k + 1; i < n; ++i)
            if (std::abs(a[static_cast<std::size_t>(i) * n + k]) >
                std::abs(a[static_cast<std::size_t>(best) * n + k]))
                best = i;
        if (std::abs(a[static_cast<std::size_t>(best) * n + k]) < 1e-300) return false;
        pivot[static_cast<std::size_t>(k)] = best;
        if (best != k)
            for (int j = 0; j < n; ++j)
                std::swap(a[static_cast<std::size_t>(k) * n + j],
                          a[static_cast<std::size_t>(best) * n + j]);
        const Complex diagonal = a[static_cast<std::size_t>(k) * n + k];
        for (int i = k + 1; i < n; ++i) {
            const Complex factor = a[static_cast<std::size_t>(i) * n + k] / diagonal;
            a[static_cast<std::size_t>(i) * n + k] = factor;
            for (int j = k + 1; j < n; ++j)
                a[static_cast<std::size_t>(i) * n + j] -=
                    factor * a[static_cast<std::size_t>(k) * n + j];
        }
    }
    return true;
}

void luSolve(const std::vector<Complex>& lu, int n, const std::vector<int>& pivot,
             std::vector<Complex>& b) {
    for (int k = 0; k < n; ++k) {
        const int p = pivot[static_cast<std::size_t>(k)];
        if (p != k) std::swap(b[static_cast<std::size_t>(k)], b[static_cast<std::size_t>(p)]);
        for (int i = k + 1; i < n; ++i)
            b[static_cast<std::size_t>(i)] -=
                lu[static_cast<std::size_t>(i) * n + k] * b[static_cast<std::size_t>(k)];
    }
    for (int i = n - 1; i >= 0; --i) {
        Complex sum = b[static_cast<std::size_t>(i)];
        for (int j = i + 1; j < n; ++j)
            sum -= lu[static_cast<std::size_t>(i) * n + j] * b[static_cast<std::size_t>(j)];
        b[static_cast<std::size_t>(i)] = sum / lu[static_cast<std::size_t>(i) * n + i];
    }
}

// Least squares by Householder QR: min ||A x - b||, A stored row-major.
bool leastSquares(std::vector<double> a, std::vector<double> b, int rows, int cols,
                  std::vector<double>& x) {
    if (rows < cols) return false;
    std::vector<double> v(static_cast<std::size_t>(rows));
    for (int k = 0; k < cols; ++k) {
        double norm = 0;
        for (int i = k; i < rows; ++i) norm += sqr(a[static_cast<std::size_t>(i) * cols + k]);
        norm = std::sqrt(norm);
        if (norm < 1e-300) return false;
        if (a[static_cast<std::size_t>(k) * cols + k] > 0) norm = -norm;
        double vv = 0;
        for (int i = k; i < rows; ++i) {
            v[static_cast<std::size_t>(i)] = a[static_cast<std::size_t>(i) * cols + k];
            if (i == k) v[static_cast<std::size_t>(i)] -= norm;
            vv += sqr(v[static_cast<std::size_t>(i)]);
        }
        if (vv < 1e-300) continue;
        for (int j = k; j < cols; ++j) {
            double s = 0;
            for (int i = k; i < rows; ++i)
                s += v[static_cast<std::size_t>(i)] * a[static_cast<std::size_t>(i) * cols + j];
            s *= 2.0 / vv;
            for (int i = k; i < rows; ++i)
                a[static_cast<std::size_t>(i) * cols + j] -= s * v[static_cast<std::size_t>(i)];
        }
        double s = 0;
        for (int i = k; i < rows; ++i)
            s += v[static_cast<std::size_t>(i)] * b[static_cast<std::size_t>(i)];
        s *= 2.0 / vv;
        for (int i = k; i < rows; ++i) b[static_cast<std::size_t>(i)] -= s * v[static_cast<std::size_t>(i)];
    }
    x.assign(static_cast<std::size_t>(cols), 0.0);
    for (int i = cols - 1; i >= 0; --i) {
        double sum = b[static_cast<std::size_t>(i)];
        for (int j = i + 1; j < cols; ++j)
            sum -= a[static_cast<std::size_t>(i) * cols + j] * x[static_cast<std::size_t>(j)];
        const double diagonal = a[static_cast<std::size_t>(i) * cols + i];
        if (std::abs(diagonal) < 1e-300) return false;
        x[static_cast<std::size_t>(i)] = sum / diagonal;
    }
    return true;
}

// Durand-Kerner: all roots of z^n + c[0] z^{n-1} + ... + c[n-1] at once.
std::vector<Complex> polynomialRoots(const std::vector<double>& c) {
    const int n = static_cast<int>(c.size());
    std::vector<Complex> root(static_cast<std::size_t>(n));
    Complex seed{0.4, 0.9};
    Complex power{1.0, 0.0};
    for (int i = 0; i < n; ++i) {
        root[static_cast<std::size_t>(i)] = power;
        power *= seed;
    }
    for (int sweep = 0; sweep < 600; ++sweep) {
        double motion = 0;
        for (int i = 0; i < n; ++i) {
            Complex value{1.0, 0.0};
            for (int j = 0; j < n; ++j)
                value = value * root[static_cast<std::size_t>(i)] + c[static_cast<std::size_t>(j)];
            Complex denominator{1.0, 0.0};
            for (int j = 0; j < n; ++j)
                if (j != i)
                    denominator *= root[static_cast<std::size_t>(i)] - root[static_cast<std::size_t>(j)];
            if (std::abs(denominator) < 1e-300) continue;
            const Complex step = value / denominator;
            root[static_cast<std::size_t>(i)] -= step;
            motion = std::max(motion, std::abs(step));
        }
        if (motion < 1e-14) break;
    }
    return root;
}

// --- Section contour ---------------------------------------------------------

// Panels on the starboard half only, keel (theta = 0) to waterline (pi/2).
// Cosine clustering puts them where the curvature and the free surface are.
std::vector<Panel> halfContour(const LewisSection& section, int panels) {
    std::vector<Panel> out;
    if (panels < 2) panels = 2;
    std::vector<double> px(static_cast<std::size_t>(panels) + 1);
    std::vector<double> py(static_cast<std::size_t>(panels) + 1);
    for (int i = 0; i <= panels; ++i) {
        const double u = static_cast<double>(i) / panels;
        const double theta = 0.5 * kPi * u;
        const Vec3 p = section.point(theta);
        px[static_cast<std::size_t>(i)] = p.x;
        py[static_cast<std::size_t>(i)] = p.y;
    }
    for (int i = 0; i < panels; ++i)
        out.push_back(makePanel(px[static_cast<std::size_t>(i)], py[static_cast<std::size_t>(i)],
                                px[static_cast<std::size_t>(i) + 1],
                                py[static_cast<std::size_t>(i) + 1]));
    return out;
}

// Generalised normal of a section mode at a panel midpoint, about the waterline
// centreline. Roll is the in-plane rotation (-y, x) . n.
double generalisedNormal(const Panel& p, int mode) {
    switch (mode) {
        case kSectionSway: return p.nx;
        case kSectionHeave: return p.ny;
        default: return -p.my * p.nx + p.mx * p.ny;
    }
}

// Sway and roll are odd in the transverse coordinate, heave is even. The mirror
// source therefore carries -1 for the first two and +1 for the third, which
// makes the symmetry exact rather than merely well converged.
double mirrorSign(int mode) { return mode == kSectionHeave ? 1.0 : -1.0; }

SectionCoefficients solveSection(const LewisSection& section, double omega, double density,
                                 int panels, bool rigidLid) {
    SectionCoefficients out;
    if (section.beam <= 0 || section.draft <= 0 || density <= 0) return out;
    if (!rigidLid && omega <= 0) return out;

    const double nu = rigidLid ? 0.0 : omega * omega / kGravity;
    const std::vector<Panel> panel = halfContour(section, panels);
    const int n = static_cast<int>(panel.size());

    // **Irregular frequencies.** A source distribution over a closed contour
    // reproduces the exterior problem correctly except at the eigenfrequencies of
    // the *interior* Dirichlet problem, where the system is near singular and the
    // answer is nonsense. They are a property of the geometry, not of the
    // discretisation, and refining the panels narrows them without removing them.
    // See the note in stripTheoryTable(), which detects and repairs them.
    const std::size_t stride = static_cast<std::size_t>(n);
    std::vector<Complex> systemEven(stride * stride), systemOdd(stride * stride);
    std::vector<Complex> potentialEven(stride * stride), potentialOdd(stride * stride);
    for (int i = 0; i < n; ++i) {
        const Panel& field = panel[static_cast<std::size_t>(i)];
        for (int j = 0; j < n; ++j) {
            const Influence self = panelInfluence(panel[static_cast<std::size_t>(j)], field.mx,
                                                  field.my, nu, i == j, rigidLid);
            const Influence other =
                panelInfluence(mirrorAcrossCentreline(panel[static_cast<std::size_t>(j)]), field.mx,
                               field.my, nu, false, rigidLid);
            const std::size_t at = static_cast<std::size_t>(i) * stride + j;
            systemEven[at] = (self.gradientX + other.gradientX) * field.nx +
                             (self.gradientY + other.gradientY) * field.ny;
            systemOdd[at] = (self.gradientX - other.gradientX) * field.nx +
                            (self.gradientY - other.gradientY) * field.ny;
            potentialEven[at] = self.potential + other.potential;
            potentialOdd[at] = self.potential - other.potential;
        }
    }

    std::vector<int> pivotEven, pivotOdd;
    std::vector<Complex> luEven = systemEven, luOdd = systemOdd;
    if (!luFactor(luEven, n, pivotEven) || !luFactor(luOdd, n, pivotOdd)) return out;

    // Solve all three modes, then pair every mode against every other through the
    // pressure integral. The whole contour is twice the half for any product of
    // two modes of the same symmetry class, and exactly zero across classes.
    std::vector<std::vector<Complex>> potentialOnHull(3);
    std::vector<Complex> farField(3, Complex(0, 0));
    for (int mode = 0; mode < 3; ++mode) {
        const bool even = (mode == kSectionHeave);
        std::vector<Complex> rhs(stride, Complex(0, 0));
        for (int i = 0; i < n; ++i)
            rhs[static_cast<std::size_t>(i)] =
                Complex(generalisedNormal(panel[static_cast<std::size_t>(i)], mode), 0.0);
        luSolve(even ? luEven : luOdd, n, even ? pivotEven : pivotOdd, rhs);

        std::vector<Complex>& phi = potentialOnHull[static_cast<std::size_t>(mode)];
        phi.assign(stride, Complex(0, 0));
        const std::vector<Complex>& potential = even ? potentialEven : potentialOdd;
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                phi[static_cast<std::size_t>(i)] +=
                    rhs[static_cast<std::size_t>(j)] *
                    potential[static_cast<std::size_t>(i) * stride + j];

        if (rigidLid) continue;
        // Radiated wave amplitude at x -> +infinity, from the source strengths.
        Complex amplitude{0, 0};
        for (int j = 0; j < n; ++j) {
            const Panel& q = panel[static_cast<std::size_t>(j)];
            Complex segment{0, 0};
            for (int g = 0; g < kGaussCount; ++g) {
                const double u = 0.5 * (kGaussNode[g] + 1.0) * q.length;
                const double xi = q.ax + u * q.tx;
                const double eta = q.ay + u * q.ty;
                const double w = kGaussWeight[g] * 0.5 * q.length * std::exp(nu * eta);
                segment += w * (std::exp(Complex(0, nu * xi)) +
                                mirrorSign(mode) * std::exp(Complex(0, -nu * xi)));
            }
            amplitude += rhs[static_cast<std::size_t>(j)] * segment;
        }
        farField[static_cast<std::size_t>(mode)] = Complex(0.0, 2.0 * kPi) * amplitude;
    }

    for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
            if ((a == kSectionHeave) != (b == kSectionHeave)) continue;   // odd, integrates out
            Complex integral{0, 0};
            for (int i = 0; i < n; ++i)
                integral += potentialOnHull[static_cast<std::size_t>(a)][static_cast<std::size_t>(i)] *
                            generalisedNormal(panel[static_cast<std::size_t>(i)], b) *
                            panel[static_cast<std::size_t>(i)].length;
            integral *= 2.0;
            out.addedMass[a][b] = -density * integral.real();
            out.damping[a][b] = rigidLid ? 0.0 : omega * density * integral.imag();
        }
    }

    if (!rigidLid) {
        // Energy check on heave: damping from the far field must equal damping
        // from the pressure integral. Two entirely different routes to the same
        // number, so their disagreement measures the discretisation directly.
        const double waveAmplitude = omega * std::abs(farField[kSectionHeave]) / kGravity;
        const double perDisplacement = waveAmplitude * omega;
        const double farDamping =
            density * kGravity * kGravity * sqr(perDisplacement) / (omega * omega * omega);
        const double near = out.damping[kSectionHeave][kSectionHeave];
        if (farDamping > 0) out.energyResidual = std::abs(near / farDamping - 1.0);
    }
    return out;
}

// --- Ogilvie transforms ------------------------------------------------------
//
// Both legs treat their input as piecewise linear between samples and integrate
// each interval against cos or sin *exactly*. That matters: the integrands
// oscillate arbitrarily fast at large t (or large omega), and a quadrature rule
// that samples them would need a mesh refined without limit, while the exact
// per-interval formula is as accurate at t = 100 s as at t = 0.
struct LinearSegment { double lo = 0, hi = 0, valueLo = 0, valueHi = 0; };

double integrateLinearTimesCos(const LinearSegment& s, double argument) {
    const double width = s.hi - s.lo;
    if (width <= 0) return 0.0;
    const double slope = (s.valueHi - s.valueLo) / width;
    if (std::abs(argument) < 1e-12)
        return s.valueLo * width + 0.5 * slope * width * width;
    const double t = argument;
    const double sinHi = std::sin(t * s.hi), sinLo = std::sin(t * s.lo);
    const double cosHi = std::cos(t * s.hi), cosLo = std::cos(t * s.lo);
    return s.valueLo * (sinHi - sinLo) / t +
           slope * (width * sinHi / t + (cosHi - cosLo) / (t * t));
}

double integrateLinearTimesSin(const LinearSegment& s, double argument) {
    const double width = s.hi - s.lo;
    if (width <= 0) return 0.0;
    const double slope = (s.valueHi - s.valueLo) / width;
    if (std::abs(argument) < 1e-12) return 0.0;
    const double t = argument;
    const double sinHi = std::sin(t * s.hi), sinLo = std::sin(t * s.lo);
    const double cosHi = std::cos(t * s.hi), cosLo = std::cos(t * s.lo);
    return s.valueLo * (cosLo - cosHi) / t +
           slope * (-width * cosHi / t + (sinHi - sinLo) / (t * t));
}

}  // namespace

// --- Lewis forms -------------------------------------------------------------

double LewisSection::sectionArea() const {
    return 0.5 * kPi * sqr(scale) * (1.0 - sqr(a1) - 3.0 * sqr(a3));
}

Vec3 LewisSection::point(double theta) const {
    const double halfBreadth = scale * ((1.0 + a1) * std::sin(theta) - a3 * std::sin(3.0 * theta));
    const double depth = scale * ((1.0 - a1) * std::cos(theta) + a3 * std::cos(3.0 * theta));
    return {halfBreadth, -depth, 0.0};
}

Vec3 LewisSection::tangent(double theta) const {
    const double dHalfBreadth =
        scale * ((1.0 + a1) * std::cos(theta) - 3.0 * a3 * std::cos(3.0 * theta));
    const double dDepth = scale * (-(1.0 - a1) * std::sin(theta) - 3.0 * a3 * std::sin(3.0 * theta));
    return {dHalfBreadth, -dDepth, 0.0};
}

// The three-parameter Lewis family in closed form. Writing H = (B/2)/T, the
// beam/draft ratio fixes a1 = c1 (1 + a3) with c1 = (H - 1)/(H + 1), and the
// area coefficient then reduces to a quadratic in u = 1 + a3:
//
//     (k + c1^2 + 3) u^2 - 6 u + 2 = 0,   k = (4 sigma / pi)(1 - c1^2)
//
// of whose two roots only the larger reproduces the circle (a3 = 0) at
// sigma = pi/4 with H = 1; the smaller lands on a3 = -1/2 there and is the
// classic way to get a plausible-looking section that is not the one asked for.
LewisSection lewisSection(double beam, double draft, double areaCoefficient) {
    LewisSection s;
    s.beam = beam;
    s.draft = draft;
    s.areaCoefficient = areaCoefficient;
    if (beam <= 0 || draft <= 0 || areaCoefficient <= 0) return s;

    const double h = (0.5 * beam) / draft;
    const double c1 = (h - 1.0) / (h + 1.0);
    const double k = (4.0 * areaCoefficient / kPi) * (1.0 - sqr(c1));
    const double d = k + sqr(c1) + 3.0;
    const double discriminant = 9.0 - 2.0 * d;
    // Outside the Lewis family the quadratic has no real root; clamp to the
    // fullest attainable section and let validateLewisSection() say so.
    const double u = (3.0 + std::sqrt(std::max(0.0, discriminant))) / d;
    s.a3 = u - 1.0;
    s.a1 = c1 * u;
    const double denominator = 1.0 + s.a1 + s.a3;
    s.scale = std::abs(denominator) > 1e-12 ? (0.5 * beam) / denominator : 0.0;
    return s;
}

std::vector<std::string> validateLewisSection(const LewisSection& section) {
    std::vector<std::string> problems;
    if (section.beam <= 0) problems.push_back("beam is not positive");
    if (section.draft <= 0) problems.push_back("draft is not positive");
    if (section.areaCoefficient <= 0) problems.push_back("sectional area coefficient is not positive");
    if (!problems.empty()) return problems;

    if (section.scale <= 0) {
        problems.push_back("the Lewis mapping degenerated");
        return problems;
    }
    const double wanted = section.areaCoefficient * section.beam * section.draft;
    const double got = section.sectionArea();
    if (std::abs(got - wanted) > 1e-9 * wanted)
        problems.push_back("sectional area coefficient " + std::to_string(section.areaCoefficient) +
                           " is outside the Lewis family for this beam/draft ratio; the closest "
                           "attainable is " + std::to_string(got / (section.beam * section.draft)));

    // The mapping folds over -- the "section" crosses itself -- when the
    // half-breadth stops increasing with theta before the waterline. A folded
    // contour still panels and still solves, and the answer is meaningless.
    double previous = -1.0;
    for (int i = 0; i <= 200; ++i) {
        const double x = section.point(0.5 * kPi * i / 200.0).x;
        if (x < previous - 1e-12) {
            problems.push_back("the Lewis form is not univalent: the section folds over itself");
            break;
        }
        previous = x;
    }
    return problems;
}

// --- Section coefficients ----------------------------------------------------

SectionCoefficients sectionCoefficients(const LewisSection& section, double omega, double density,
                                        int panels) {
    return solveSection(section, omega, density, panels, false);
}

SectionCoefficients sectionCoefficientsInfiniteFrequency(const LewisSection& section,
                                                         double density, int panels) {
    return solveSection(section, 0.0, density, panels, true);
}

// --- Strip theory ------------------------------------------------------------

std::vector<double> radiationFrequencyGrid(double omegaMin, double omegaMax, int count) {
    std::vector<double> grid;
    if (count < 2 || omegaMin <= 0 || omegaMax <= omegaMin) return grid;
    const double ratio = std::pow(omegaMax / omegaMin, 1.0 / (count - 1));
    double omega = omegaMin;
    for (int i = 0; i < count; ++i) {
        grid.push_back(omega);
        omega *= ratio;
    }
    grid.back() = omegaMax;
    return grid;
}

Matrix6 transferAddedMass(const Matrix6& matrix, const Vec3& offset) {
    // T = [[I, d~], [0, I]]; see the derivation in the header. Written out rather
    // than looped over as a general 6x6 product only in the sense that the two
    // multiplications are kept separate -- the fully expanded quadruple loop is
    // 1296 products for a matrix that is mostly identity, and this runs once a
    // tick inside the rigid-body integrator.
    double t[6][6] = {};
    for (int i = 0; i < 6; ++i) t[i][i] = 1.0;
    t[0][4] = -offset.z; t[0][5] =  offset.y;
    t[1][3] =  offset.z; t[1][5] = -offset.x;
    t[2][3] = -offset.y; t[2][4] =  offset.x;

    double half[6][6] = {};   // A * T
    for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j) {
            double sum = 0;
            for (int k = 0; k < 6; ++k) sum += matrix[static_cast<std::size_t>(i)]
                                                     [static_cast<std::size_t>(k)] * t[k][j];
            half[i][j] = sum;
        }

    Matrix6 out{};
    for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j) {
            double sum = 0;
            for (int k = 0; k < 6; ++k) sum += t[k][i] * half[k][j];
            out[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] = sum;
        }
    return out;
}

double RadiationTable::dampingAt(int i, int j, double omegaQuery) const {
    const int n = size();
    if (n == 0 || omegaQuery <= omega.front() || omegaQuery >= omega.back()) return 0.0;
    const std::size_t upper = static_cast<std::size_t>(
        std::lower_bound(omega.begin(), omega.end(), omegaQuery) - omega.begin());
    const std::size_t lower = upper - 1;
    const double t = (omegaQuery - omega[lower]) / (omega[upper] - omega[lower]);
    const double a = damping[lower][static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
    const double b = damping[upper][static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
    return a + t * (b - a);
}

namespace {

// Accumulate one station's 3x3 sectional matrix (sway, heave, roll about the
// *baseline* origin) into the 6x6 hull matrix, weighted by its share of length.
void accumulateStation(Matrix6& out, const double sectional[3][3], double x, double weight) {
    // Mode -> DOF, with the lever that strip theory gives each coupling.
    //   sway  -> sway (1) and yaw (5) with lever +x
    //   heave -> heave (2) and pitch (4) with lever -x
    //   roll  -> roll (3)
    const int dof[3][2] = {{1, 5}, {2, 4}, {3, -1}};
    const double lever[3][2] = {{1.0, x}, {1.0, -x}, {1.0, 0.0}};
    for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b) {
            if (sectional[a][b] == 0.0) continue;
            for (int p = 0; p < 2; ++p) {
                if (dof[a][p] < 0) continue;
                for (int q = 0; q < 2; ++q) {
                    if (dof[b][q] < 0) continue;
                    out[static_cast<std::size_t>(dof[a][p])][static_cast<std::size_t>(dof[b][q])] +=
                        weight * lever[a][p] * lever[b][q] * sectional[a][b];
                }
            }
        }
}

// Move a sectional matrix from the waterline centreline to the baseline
// centreline: the roll mode picks up -draft times the sway mode.
void transferToBaseline(double m[3][3], double draft) {
    const double s[3][3] = {{1, 0, -draft}, {0, 1, 0}, {0, 0, 1}};
    double tmp[3][3] = {};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int a = 0; a < 3; ++a)
                for (int b = 0; b < 3; ++b) tmp[i][j] += s[a][i] * m[a][b] * s[b][j];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) m[i][j] = tmp[i][j];
}

// A section solve whose near-field and far-field damping disagree by more than
// this is not a discretisation error, it is an irregular frequency: measured
// backgrounds reach 0.14 on the ferry's table -- 4.5e-3 on a smooth cylinder,
// which is where the "under 0.01" this used to claim came from -- and the spikes
// run to 60. So the margin between the background and this threshold is about
// 1.4x, not the 20x the old figure implied.
constexpr double kIrregularResidual = 0.2;

// Trapezium weight for station i: half the distance to each neighbour.
double stationWeight(const std::vector<RadiationStation>& stations, std::size_t i) {
    const std::size_t n = stations.size();
    if (n < 2) return 0.0;
    if (i == 0) return 0.5 * (stations[1].x - stations[0].x);
    if (i == n - 1) return 0.5 * (stations[n - 1].x - stations[n - 2].x);
    return 0.5 * (stations[i + 1].x - stations[i - 1].x);
}

}  // namespace

RadiationTable stripTheoryTable(const RadiationHull& hull, const std::vector<double>& omega) {
    RadiationTable table;
    if (hull.stations.size() < 2 || omega.empty()) return table;
    table.omega = omega;
    table.addedMass.assign(omega.size(), Matrix6{});
    table.damping.assign(omega.size(), Matrix6{});

    std::vector<LewisSection> section;
    section.reserve(hull.stations.size());
    for (const RadiationStation& s : hull.stations)
        section.push_back(lewisSection(s.beam, s.draft, s.areaCoefficient));

    for (std::size_t i = 0; i < hull.stations.size(); ++i) {
        const double weight = stationWeight(hull.stations, i);
        if (weight <= 0 || section[i].scale <= 0) continue;

        SectionCoefficients infinite = sectionCoefficientsInfiniteFrequency(
            section[i], hull.density, hull.panelsPerHalfSection);
        transferToBaseline(infinite.addedMass, hull.draft);
        accumulateStation(table.addedMassInfinite, infinite.addedMass, hull.stations[i].x, weight);

        // **Irregular frequencies, detected and repaired.** The source
        // formulation used per section is exact except at the eigenfrequencies of
        // the interior Dirichlet problem, where it returns nonsense -- measured
        // on a 25 x 6.5 m section as spikes at omega = 1.33, 1.94, 2.48, 2.94
        // rad/s, matching the closed form nu_m = (m pi / B) coth(m pi T / B) for
        // odd m, and coming out with *negative damping*: a ship extracting energy
        // from still water. Refining the panels narrows the spikes and does not
        // remove them, because they are a property of the geometry.
        //
        // The proper cure is an extended integral equation with unknowns on an
        // interior lid; that is not implemented here. What is implemented is
        // detection -- the near-field/far-field energy balance fails by orders of
        // magnitude exactly there, and by under a percent everywhere else, so the
        // signal is unambiguous -- followed by interpolation from the neighbours
        // that passed. The count is reported so that a hull whose whole grid is
        // being repaired cannot pass for a clean one.
        std::vector<SectionCoefficients> byFrequency(omega.size());
        std::vector<bool> rejected(omega.size(), false);
        for (std::size_t f = 0; f < omega.size(); ++f) {
            byFrequency[f] = sectionCoefficients(section[i], omega[f], hull.density,
                                                 hull.panelsPerHalfSection);
            const SectionCoefficients& c = byFrequency[f];
            rejected[f] = c.energyResidual > kIrregularResidual || c.damping[0][0] < 0 ||
                          c.damping[1][1] < 0 || c.damping[2][2] < 0;
            table.totalSolves += 1;
        }
        for (std::size_t f = 0; f < omega.size(); ++f) {
            if (!rejected[f]) {
                table.worstEnergyResidual =
                    std::max(table.worstEnergyResidual, byFrequency[f].energyResidual);
                continue;
            }
            table.repairedSolves += 1;
            std::size_t lo = f, hi = f;
            bool haveLo = false, haveHi = false;
            while (lo > 0) {
                --lo;
                if (!rejected[lo]) { haveLo = true; break; }
            }
            while (hi + 1 < omega.size()) {
                ++hi;
                if (!rejected[hi]) { haveHi = true; break; }
            }
            if (!haveLo && !haveHi) continue;   // nothing clean to lean on
            if (!haveLo) byFrequency[f] = byFrequency[hi];
            else if (!haveHi) byFrequency[f] = byFrequency[lo];
            else {
                const double t = (omega[f] - omega[lo]) / (omega[hi] - omega[lo]);
                for (int a = 0; a < 3; ++a)
                    for (int b = 0; b < 3; ++b) {
                        byFrequency[f].addedMass[a][b] =
                            byFrequency[lo].addedMass[a][b] +
                            t * (byFrequency[hi].addedMass[a][b] - byFrequency[lo].addedMass[a][b]);
                        byFrequency[f].damping[a][b] =
                            byFrequency[lo].damping[a][b] +
                            t * (byFrequency[hi].damping[a][b] - byFrequency[lo].damping[a][b]);
                    }
            }
        }
        for (std::size_t f = 0; f < omega.size(); ++f) {
            SectionCoefficients c = byFrequency[f];
            transferToBaseline(c.addedMass, hull.draft);
            transferToBaseline(c.damping, hull.draft);
            accumulateStation(table.addedMass[f], c.addedMass, hull.stations[i].x, weight);
            accumulateStation(table.damping[f], c.damping, hull.stations[i].x, weight);
        }
    }
    return table;
}

// --- Cummins / Ogilvie -------------------------------------------------------

std::vector<double> retardationFunction(const RadiationTable& table, int i, int j, double dt,
                                        int count) {
    std::vector<double> k;
    if (count <= 0 || dt <= 0 || table.size() < 2) return k;
    k.assign(static_cast<std::size_t>(count), 0.0);

    // B is taken as piecewise linear on the grid, running to zero at omega = 0
    // (which it does, exactly) and cut off above the last grid point. The cut-off
    // is the one approximation in this transform and the tests measure what it
    // leaves behind.
    std::vector<LinearSegment> segment;
    segment.push_back({0.0, table.omega.front(), 0.0,
                       table.damping.front()[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)]});
    for (std::size_t s = 0; s + 1 < table.omega.size(); ++s)
        segment.push_back({table.omega[s], table.omega[s + 1],
                           table.damping[s][static_cast<std::size_t>(i)][static_cast<std::size_t>(j)],
                           table.damping[s + 1][static_cast<std::size_t>(i)][static_cast<std::size_t>(j)]});

    for (int n = 0; n < count; ++n) {
        const double t = n * dt;
        double sum = 0;
        for (const LinearSegment& s : segment) sum += integrateLinearTimesCos(s, t);
        k[static_cast<std::size_t>(n)] = 2.0 / kPi * sum;
    }
    return k;
}

double dampingFromRetardation(const std::vector<double>& k, double dt, double omega) {
    if (k.size() < 2 || dt <= 0) return 0.0;
    double sum = 0;
    for (std::size_t n = 0; n + 1 < k.size(); ++n)
        sum += integrateLinearTimesCos(
            {n * dt, (n + 1) * dt, k[n], k[n + 1]}, omega);
    return sum;
}

double addedMassFromRetardation(const std::vector<double>& k, double dt, double omega,
                                double addedMassInfinite) {
    if (k.size() < 2 || dt <= 0 || omega <= 0) return addedMassInfinite;
    double sum = 0;
    for (std::size_t n = 0; n + 1 < k.size(); ++n)
        sum += integrateLinearTimesSin({n * dt, (n + 1) * dt, k[n], k[n + 1]}, omega);
    return addedMassInfinite - sum / omega;
}

double infiniteAddedMassFromRetardation(const std::vector<double>& k, double dt, double omega,
                                        double addedMassAtOmega) {
    if (k.size() < 2 || dt <= 0 || omega <= 0) return addedMassAtOmega;
    double sum = 0;
    for (std::size_t n = 0; n + 1 < k.size(); ++n)
        sum += integrateLinearTimesSin({n * dt, (n + 1) * dt, k[n], k[n + 1]}, omega);
    return addedMassAtOmega + sum / omega;
}

double memoryDecayTime(const std::vector<double>& k, double dt, double fraction) {
    if (k.empty() || dt <= 0) return 0.0;
    double peak = 0;
    for (double value : k) peak = std::max(peak, std::abs(value));
    if (peak <= 0) return 0.0;
    const double threshold = fraction * peak;
    std::size_t last = 0;
    for (std::size_t n = 0; n < k.size(); ++n)
        if (std::abs(k[n]) > threshold) last = n;
    return static_cast<double>(last) * dt;
}

// --- State space -------------------------------------------------------------

int RadiationStateSpace::stateCount() const {
    int n = 0;
    for (const Mode& m : modes) n += (m.frequency != 0.0) ? 2 : 1;
    return n;
}

double RadiationStateSpace::impulseResponse(double t) const {
    double sum = 0;
    for (const Mode& m : modes) {
        const double envelope = std::exp(-m.decay * t);
        sum += (m.frequency != 0.0)
                   ? envelope * (m.c0 * std::cos(m.frequency * t) - m.c1 * std::sin(m.frequency * t))
                   : envelope * m.c0;
    }
    return sum;
}

bool RadiationStateSpace::stable() const {
    for (const Mode& m : modes)
        if (!(m.decay > 0.0)) return false;
    return !modes.empty();
}

// Least-squares Prony. Fit a linear recurrence to the samples, root the
// characteristic polynomial to get the poles, then solve linearly for the
// residues. Exact when K really is a sum of damped sinusoids, which is what
// tests/test_radiation.cpp checks before this is turned loose on a hull.
StateSpaceFit fitStateSpace(const std::vector<double>& k, double dt, int order) {
    StateSpaceFit fit;
    const int samples = static_cast<int>(k.size());
    if (order < 1 || dt <= 0 || samples < 4 * order) return fit;

    double peak = 0;
    for (double value : k) peak = std::max(peak, std::abs(value));
    if (peak <= 0) return fit;

    // Linear prediction: k[m] = -sum a_i k[m-i], overdetermined, least squares.
    const int rows = samples - order;
    std::vector<double> a(static_cast<std::size_t>(rows) * order);
    std::vector<double> rhs(static_cast<std::size_t>(rows));
    for (int m = 0; m < rows; ++m) {
        for (int i = 0; i < order; ++i)
            a[static_cast<std::size_t>(m) * order + i] =
                k[static_cast<std::size_t>(m + order - 1 - i)] / peak;
        rhs[static_cast<std::size_t>(m)] = -k[static_cast<std::size_t>(m + order)] / peak;
    }
    std::vector<double> coefficient;
    if (!leastSquares(a, rhs, rows, order, coefficient)) return fit;

    const std::vector<Complex> root = polynomialRoots(coefficient);

    // Poles: z = exp(lambda dt). Anything on or outside the unit circle is
    // reflected back inside -- an unstable radiation model is not merely
    // inaccurate, it diverges -- and complex roots are kept once per conjugate
    // pair.
    struct Pole { double decay, frequency; };
    std::vector<Pole> pole;
    for (const Complex& z : root) {
        double magnitude = std::abs(z);
        if (magnitude < 1e-12) continue;
        if (magnitude >= 1.0) magnitude = 1.0 / magnitude;      // reflect
        const double decay = -std::log(magnitude) / dt;
        const double frequency = std::abs(std::arg(z)) / dt;
        if (!std::isfinite(decay) || !std::isfinite(frequency)) continue;
        if (frequency > 1e-9 && std::arg(z) < 0) continue;      // keep one of each pair
        pole.push_back({decay, frequency > 1e-9 ? frequency : 0.0});
    }
    if (pole.empty()) return fit;

    // Residues by linear least squares against the sampled K.
    int columns = 0;
    for (const Pole& p : pole) columns += (p.frequency != 0.0) ? 2 : 1;
    if (samples < columns) return fit;
    std::vector<double> basis(static_cast<std::size_t>(samples) * columns);
    std::vector<double> target(static_cast<std::size_t>(samples));
    for (int n = 0; n < samples; ++n) {
        const double t = n * dt;
        int column = 0;
        for (const Pole& p : pole) {
            const double envelope = std::exp(-p.decay * t);
            basis[static_cast<std::size_t>(n) * columns + column++] =
                envelope * ((p.frequency != 0.0) ? std::cos(p.frequency * t) : 1.0);
            if (p.frequency != 0.0)
                basis[static_cast<std::size_t>(n) * columns + column++] =
                    envelope * std::sin(p.frequency * t);
        }
        target[static_cast<std::size_t>(n)] = k[static_cast<std::size_t>(n)];
    }
    std::vector<double> residue;
    if (!leastSquares(basis, target, samples, columns, residue)) return fit;

    int column = 0;
    for (const Pole& p : pole) {
        RadiationStateSpace::Mode mode;
        mode.decay = p.decay;
        mode.frequency = p.frequency;
        mode.c0 = residue[static_cast<std::size_t>(column++)];
        // C exp(A t) B = e^{-sigma t}(c0 cos - c1 sin), so the sine coefficient
        // enters negated.
        mode.c1 = (p.frequency != 0.0) ? -residue[static_cast<std::size_t>(column++)] : 0.0;
        fit.model.modes.push_back(mode);
    }

    double sumSquare = 0, sumSquareK = 0;
    for (int n = 0; n < samples; ++n) {
        const double error =
            fit.model.impulseResponse(n * dt) - k[static_cast<std::size_t>(n)];
        sumSquare += sqr(error);
        sumSquareK += sqr(k[static_cast<std::size_t>(n)]);
        fit.peakError = std::max(fit.peakError, std::abs(error));
    }
    fit.rmsError = std::sqrt(sumSquare / samples);
    fit.relativeRms = sumSquareK > 0 ? fit.rmsError / std::sqrt(sumSquareK / samples) : 0.0;
    fit.converged = fit.model.stable() && std::isfinite(fit.rmsError);
    return fit;
}

// --- Runtime evaluator -------------------------------------------------------

RadiationForce::RadiationForce(const RadiationTable& table, double dt, int samples, int order) {
    addedMassInfinite_ = table.addedMassInfinite;
    if (table.size() < 2 || dt <= 0 || samples <= 0) return;

    // Scale is set by the largest damping anywhere in the table, so that a
    // coupling term three orders down does not get its own state-space model for
    // what is numerical dust.
    double largest = 0;
    for (const Matrix6& m : table.damping)
        for (int i = 0; i < 6; ++i)
            for (int j = 0; j < 6; ++j)
                largest = std::max(largest, std::abs(m[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)]));
    if (largest <= 0) return;

    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            double entryPeak = 0;
            for (const Matrix6& m : table.damping)
                entryPeak = std::max(entryPeak,
                                     std::abs(m[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)]));
            if (entryPeak < 1e-6 * largest) continue;

            const std::vector<double> k = retardationFunction(table, i, j, dt, samples);
            const StateSpaceFit fit = fitStateSpace(k, dt, order);
            if (!fit.converged) continue;
            Entry entry;
            entry.i = i;
            entry.j = j;
            entry.model = fit.model;
            entry.state.assign(static_cast<std::size_t>(fit.model.stateCount()), 0.0);
            entries_.push_back(entry);
            if (entryPeak > 0.1 * largest)
                worstRelativeRms_ = std::max(worstRelativeRms_, fit.relativeRms);
        }
    }
}

void RadiationForce::reset() {
    for (Entry& e : entries_) std::fill(e.state.begin(), e.state.end(), 0.0);
}

// Exact zero-order-hold update of each block. Closed form, so the step is
// unconditionally stable and independent of how dt compares with the fitted
// sampling interval -- an explicit integrator here would go unstable on exactly
// the stiff, fast-decaying modes the fit produces.
void RadiationForce::step(const std::array<double, 6>& velocity, double dt) {
    if (dt <= 0) return;
    for (Entry& e : entries_) {
        const double v = velocity[static_cast<std::size_t>(e.j)];
        std::size_t at = 0;
        for (const RadiationStateSpace::Mode& m : e.model.modes) {
            const double decayFactor = std::exp(-m.decay * dt);
            if (m.frequency != 0.0) {
                const double c = std::cos(m.frequency * dt);
                const double s = std::sin(m.frequency * dt);
                const double x0 = e.state[at], x1 = e.state[at + 1];
                const double u = decayFactor * c - 1.0;
                const double w = -decayFactor * s;
                const double denominator = sqr(m.decay) + sqr(m.frequency);
                const double g0 = (-m.decay * u - m.frequency * w) / denominator;
                const double g1 = (m.frequency * u - m.decay * w) / denominator;
                e.state[at] = decayFactor * (c * x0 + s * x1) + g0 * v;
                e.state[at + 1] = decayFactor * (-s * x0 + c * x1) + g1 * v;
                at += 2;
            } else {
                const double gain =
                    (m.decay > 1e-12) ? (1.0 - decayFactor) / m.decay : dt;
                e.state[at] = decayFactor * e.state[at] + gain * v;
                at += 1;
            }
        }
    }
}

std::array<double, 6> RadiationForce::memoryForce() const {
    std::array<double, 6> force{};
    for (const Entry& e : entries_) {
        double sum = 0;
        std::size_t at = 0;
        for (const RadiationStateSpace::Mode& m : e.model.modes) {
            if (m.frequency != 0.0) {
                sum += m.c0 * e.state[at] + m.c1 * e.state[at + 1];
                at += 2;
            } else {
                sum += m.c0 * e.state[at];
                at += 1;
            }
        }
        force[static_cast<std::size_t>(e.i)] += sum;
    }
    return force;
}

// --- Stations from a hull mesh -----------------------------------------------

RadiationHull radiationHullFromMesh(const TriMesh& hull, double waterlineZ, int stationCount,
                                    double density) {
    RadiationHull out;
    out.density = density;
    if (hull.verts.empty() || stationCount < 2) return out;

    Vec3 lo = hull.verts[0], hi = hull.verts[0];
    for (const Vec3& v : hull.verts) {
        lo = {std::min(lo.x, v.x), std::min(lo.y, v.y), std::min(lo.z, v.z)};
        hi = {std::max(hi.x, v.x), std::max(hi.y, v.y), std::max(hi.z, v.z)};
    }
    if (!(waterlineZ > lo.z)) return out;
    out.draft = waterlineZ - lo.z;

    // Inset from the extreme ends. A slab at the stem or transom collects almost
    // no volume, so beam and area both go to zero and their ratio -- the area
    // coefficient -- is whatever the rounding says. Strip theory has nothing
    // useful to contribute there in any case.
    const double length = hi.x - lo.x;
    const double inset = 0.5 * length / stationCount;
    const double first = lo.x + inset;
    const double span = (length - 2.0 * inset);
    const double thickness = span / (stationCount - 1);

    const double wide = (hi.y - lo.y) + 1.0;
    for (int i = 0; i < stationCount; ++i) {
        const double x = first + span * i / (stationCount - 1);
        const double half = 0.5 * thickness;
        // The wetted slab: everything between two transverse planes and below the
        // waterline. clipToBox is the same routine the compartment carving uses.
        const TriMesh slab = clipToBox(hull, {x - half, lo.y - wide, lo.z - wide},
                                       {x + half, hi.y + wide, waterlineZ});
        if (slab.verts.empty()) continue;

        double beamLo = slab.verts[0].y, beamHi = slab.verts[0].y, keel = slab.verts[0].z;
        for (const Vec3& v : slab.verts) {
            beamLo = std::min(beamLo, v.y);
            beamHi = std::max(beamHi, v.y);
            keel = std::min(keel, v.z);
        }

        RadiationStation station;
        station.x = x;
        station.beam = beamHi - beamLo;
        station.draft = waterlineZ - keel;
        // Sectional area is the slab's volume over its thickness. Taking it from
        // the same integrator the hydrostatics use means a hull that displaces
        // what it should also has stations that add up to what they should.
        const double area = integrate(slab).volume / thickness;
        if (station.beam <= 0 || station.draft <= 0 || area <= 0) continue;
        station.areaCoefficient = area / (station.beam * station.draft);
        out.stations.push_back(station);
    }
    return out;
}

// --- Validity ----------------------------------------------------------------

std::vector<std::string> validateRadiationHull(const RadiationHull& hull,
                                               const std::vector<double>& omega) {
    std::vector<std::string> problems;
    if (hull.stations.size() < 2) problems.push_back("fewer than two stations");
    if (hull.draft <= 0) problems.push_back("draft is not positive");
    if (hull.density <= 0) problems.push_back("density is not positive");
    if (hull.panelsPerHalfSection < 8)
        problems.push_back("fewer than 8 panels per half section is too coarse to trust");
    if (omega.size() < 2) problems.push_back("fewer than two frequencies");
    if (!problems.empty()) return problems;

    for (std::size_t i = 0; i + 1 < hull.stations.size(); ++i)
        if (!(hull.stations[i + 1].x > hull.stations[i].x)) {
            problems.push_back("stations are not strictly increasing in x");
            break;
        }
    for (std::size_t i = 0; i + 1 < omega.size(); ++i)
        if (!(omega[i + 1] > omega[i])) {
            problems.push_back("frequencies are not strictly increasing");
            break;
        }

    for (const RadiationStation& s : hull.stations) {
        const LewisSection section = lewisSection(s.beam, s.draft, s.areaCoefficient);
        for (const std::string& why : validateLewisSection(section))
            problems.push_back("station at x = " + std::to_string(s.x) + ": " + why);
    }

    // Strip theory is a slenderness assumption: the flow at a station is taken
    // to be two-dimensional, which needs the section to be small against the
    // length and against the wavelength.
    const double length = hull.stations.back().x - hull.stations.front().x;
    double widest = 0;
    for (const RadiationStation& s : hull.stations) widest = std::max(widest, s.beam);
    if (length > 0 && widest / length > 0.25)
        problems.push_back("beam/length is " + std::to_string(widest / length) +
                           "; strip theory assumes a slender hull, conventionally B/L below 0.25");
    // Strip theory is squeezed from both ends: the wavelength has to be short
    // against the length for the sections to be independent, and long against
    // the section for the flow at a station to be two-dimensional. A beamy ship
    // cannot satisfy both comfortably, so the flag is set where the shorter
    // requirement actually breaks -- a wavelength under a beam -- rather than at
    // a comfortable margin that would condemn every useful band.
    const double shortest = 2.0 * kPi * kGravity / sqr(omega.back());
    if (shortest < widest)
        problems.push_back("the highest frequency has a wavelength of " + std::to_string(shortest) +
                           " m against a beam of " + std::to_string(widest) +
                           " m; strip theory needs the wavelength long against the section");
    if (hull.stations.front().x * hull.stations.back().x >= 0)
        problems.push_back("stations do not straddle midship, so the pitch and yaw levers are "
                           "measured about the wrong point");
    return problems;
}

}  // namespace sim
