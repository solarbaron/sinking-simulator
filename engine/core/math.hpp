// SPDX-License-Identifier: MIT
// Minimal SI-unit math core. Right-handed, z-up.
// Ship body frame: +x forward (bow), +y to port, +z up. Origin at midship baseline.
#pragma once

#include <algorithm>
#include <array>
#include <cmath>

namespace sim {

inline constexpr double kPi           = 3.14159265358979323846;
inline constexpr double kRadToDeg     = 180.0 / kPi;
inline constexpr double kDegToRad     = kPi / 180.0;
inline constexpr double kGravity      = 9.80665;   // m/s^2
inline constexpr double kRhoSeawater  = 1025.0;    // kg/m^3
inline constexpr double kRhoFresh     = 998.2;     // kg/m^3
inline constexpr double kPatm         = 101325.0;  // Pa
inline constexpr double kRAir         = 287.052;   // J/(kg K)
inline constexpr double kGammaAir     = 1.4;       // adiabatic index
inline constexpr double kTAmbient     = 288.15;    // K

struct Vec3 {
    double x = 0, y = 0, z = 0;

    constexpr Vec3() = default;
    constexpr Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

    constexpr Vec3 operator+(const Vec3& o) const { return {x + o.x, y + o.y, z + o.z}; }
    constexpr Vec3 operator-(const Vec3& o) const { return {x - o.x, y - o.y, z - o.z}; }
    constexpr Vec3 operator-()              const { return {-x, -y, -z}; }
    constexpr Vec3 operator*(double s)      const { return {x * s, y * s, z * s}; }
    constexpr Vec3 operator/(double s)      const { return {x / s, y / s, z / s}; }
    constexpr Vec3& operator+=(const Vec3& o) { x += o.x; y += o.y; z += o.z; return *this; }
    constexpr Vec3& operator-=(const Vec3& o) { x -= o.x; y -= o.y; z -= o.z; return *this; }
    constexpr Vec3& operator*=(double s)      { x *= s; y *= s; z *= s; return *this; }

    constexpr double& operator[](int i)       { return i == 0 ? x : (i == 1 ? y : z); }
    constexpr double  operator[](int i) const { return i == 0 ? x : (i == 1 ? y : z); }
};

constexpr Vec3 operator*(double s, const Vec3& v) { return v * s; }
constexpr double dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
constexpr Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
inline double length(const Vec3& v)  { return std::sqrt(dot(v, v)); }
inline double length2(const Vec3& v) { return dot(v, v); }
inline Vec3 normalize(const Vec3& v) {
    const double l = length(v);
    return l > 1e-300 ? v / l : Vec3{0, 0, 0};
}

// Row-major 3x3.
struct Mat3 {
    std::array<double, 9> m{1, 0, 0, 0, 1, 0, 0, 0, 1};

    constexpr double& operator()(int r, int c)       { return m[r * 3 + c]; }
    constexpr double  operator()(int r, int c) const { return m[r * 3 + c]; }

    constexpr Vec3 operator*(const Vec3& v) const {
        return {m[0] * v.x + m[1] * v.y + m[2] * v.z,
                m[3] * v.x + m[4] * v.y + m[5] * v.z,
                m[6] * v.x + m[7] * v.y + m[8] * v.z};
    }
    constexpr Mat3 operator*(const Mat3& o) const {
        Mat3 r;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) {
                double s = 0;
                for (int k = 0; k < 3; ++k) s += (*this)(i, k) * o(k, j);
                r(i, j) = s;
            }
        return r;
    }
    constexpr Mat3 operator+(const Mat3& o) const {
        Mat3 r;
        for (int i = 0; i < 9; ++i) r.m[i] = m[i] + o.m[i];
        return r;
    }
    constexpr Mat3 operator*(double s) const {
        Mat3 r;
        for (int i = 0; i < 9; ++i) r.m[i] = m[i] * s;
        return r;
    }
    constexpr Mat3 transposed() const {
        Mat3 r;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) r(i, j) = (*this)(j, i);
        return r;
    }
    static constexpr Mat3 identity() { return {}; }
    static constexpr Mat3 zero() { Mat3 r; r.m.fill(0); return r; }
};

inline Mat3 inverse(const Mat3& a) {
    const double c00 = a(1, 1) * a(2, 2) - a(1, 2) * a(2, 1);
    const double c01 = a(1, 2) * a(2, 0) - a(1, 0) * a(2, 2);
    const double c02 = a(1, 0) * a(2, 1) - a(1, 1) * a(2, 0);
    const double det = a(0, 0) * c00 + a(0, 1) * c01 + a(0, 2) * c02;
    Mat3 r = Mat3::zero();
    if (std::abs(det) < 1e-300) return r;
    const double id = 1.0 / det;
    r(0, 0) = c00 * id;
    r(1, 0) = c01 * id;
    r(2, 0) = c02 * id;
    r(0, 1) = (a(0, 2) * a(2, 1) - a(0, 1) * a(2, 2)) * id;
    r(1, 1) = (a(0, 0) * a(2, 2) - a(0, 2) * a(2, 0)) * id;
    r(2, 1) = (a(0, 1) * a(2, 0) - a(0, 0) * a(2, 1)) * id;
    r(0, 2) = (a(0, 1) * a(1, 2) - a(0, 2) * a(1, 1)) * id;
    r(1, 2) = (a(0, 2) * a(1, 0) - a(0, 0) * a(1, 2)) * id;
    r(2, 2) = (a(0, 0) * a(1, 1) - a(0, 1) * a(1, 0)) * id;
    return r;
}

// Parallel-axis / point-mass inertia contribution of mass m at offset r about the origin.
inline Mat3 pointInertia(double m, const Vec3& r) {
    Mat3 I = Mat3::zero();
    I(0, 0) = m * (r.y * r.y + r.z * r.z);
    I(1, 1) = m * (r.x * r.x + r.z * r.z);
    I(2, 2) = m * (r.x * r.x + r.y * r.y);
    I(0, 1) = I(1, 0) = -m * r.x * r.y;
    I(0, 2) = I(2, 0) = -m * r.x * r.z;
    I(1, 2) = I(2, 1) = -m * r.y * r.z;
    return I;
}

// Column-major 4x4, matching GLSL's mat4 memory layout so a matrix can be pushed
// straight to a shader with no transpose. m[column * 4 + row].
struct Mat4 {
    std::array<double, 16> m{1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};

    constexpr double& operator()(int row, int column) { return m[column * 4 + row]; }
    constexpr double  operator()(int row, int column) const { return m[column * 4 + row]; }

    static constexpr Mat4 identity() { return {}; }

    constexpr Mat4 operator*(const Mat4& other) const {
        Mat4 result;
        for (int c = 0; c < 4; ++c)
            for (int r = 0; r < 4; ++r) {
                double sum = 0;
                for (int k = 0; k < 4; ++k) sum += (*this)(r, k) * other(k, c);
                result(r, c) = sum;
            }
        return result;
    }

    // Full homogeneous transform, returning w alongside so callers can do the
    // perspective divide themselves and detect points behind the eye.
    constexpr void transform(const Vec3& point, double out[4]) const {
        for (int r = 0; r < 4; ++r)
            out[r] = (*this)(r, 0) * point.x + (*this)(r, 1) * point.y +
                     (*this)(r, 2) * point.z + (*this)(r, 3);
    }

    void toFloats(float out[16]) const {
        for (int i = 0; i < 16; ++i) out[i] = static_cast<float>(m[i]);
    }
};

// Right-handed view matrix: the camera sits at `eye` looking toward `target`.
inline Mat4 lookAt(const Vec3& eye, const Vec3& target, const Vec3& up) {
    const Vec3 forward = normalize(target - eye);
    const Vec3 right = normalize(cross(forward, up));
    const Vec3 trueUp = cross(right, forward);

    Mat4 view;
    view(0, 0) = right.x;   view(0, 1) = right.y;   view(0, 2) = right.z;
    view(1, 0) = trueUp.x;  view(1, 1) = trueUp.y;  view(1, 2) = trueUp.z;
    view(2, 0) = -forward.x; view(2, 1) = -forward.y; view(2, 2) = -forward.z;
    view(0, 3) = -dot(right, eye);
    view(1, 3) = -dot(trueUp, eye);
    view(2, 3) = dot(forward, eye);
    view(3, 0) = 0; view(3, 1) = 0; view(3, 2) = 0; view(3, 3) = 1;
    return view;
}

// Perspective projection into **Vulkan** clip space: x and y in [-1, 1] with y
// pointing *down* the screen, and z in [0, 1] rather than OpenGL's [-1, 1].
// Getting either convention wrong renders a vertically mirrored image or a depth
// buffer that only uses half its range -- both of which look almost right.
inline Mat4 perspective(double fovYRadians, double aspect, double nearPlane, double farPlane) {
    const double f = 1.0 / std::tan(fovYRadians * 0.5);
    Mat4 projection;
    projection.m.fill(0.0);
    projection(0, 0) = f / aspect;
    projection(1, 1) = -f;  // Vulkan's y axis points down
    projection(2, 2) = farPlane / (nearPlane - farPlane);
    projection(2, 3) = nearPlane * farPlane / (nearPlane - farPlane);
    projection(3, 2) = -1.0;
    return projection;
}

// Clip-space position to pixel coordinates, doing the perspective divide.
// Returns false for points at or behind the eye, where the divide is meaningless.
inline bool clipToPixel(const double clip[4], double width, double height, double& outX,
                        double& outY) {
    if (clip[3] <= 1e-12) return false;
    outX = (clip[0] / clip[3] * 0.5 + 0.5) * width;
    outY = (clip[1] / clip[3] * 0.5 + 0.5) * height;
    return true;
}

struct Quat {
    double w = 1, x = 0, y = 0, z = 0;

    constexpr Quat() = default;
    constexpr Quat(double w_, double x_, double y_, double z_) : w(w_), x(x_), y(y_), z(z_) {}

    static Quat fromAxisAngle(const Vec3& axis, double angle) {
        const Vec3 a = normalize(axis);
        const double s = std::sin(0.5 * angle);
        return {std::cos(0.5 * angle), a.x * s, a.y * s, a.z * s};
    }

    constexpr Quat operator*(const Quat& o) const {
        return {w * o.w - x * o.x - y * o.y - z * o.z,
                w * o.x + x * o.w + y * o.z - z * o.y,
                w * o.y - x * o.z + y * o.w + z * o.x,
                w * o.z + x * o.y - y * o.x + z * o.w};
    }

    void normalizeInPlace() {
        const double l = std::sqrt(w * w + x * x + y * y + z * z);
        if (l > 1e-300) { w /= l; x /= l; y /= l; z /= l; }
    }

    Mat3 toMat3() const {
        Mat3 r;
        const double xx = x * x, yy = y * y, zz = z * z;
        const double xy = x * y, xz = x * z, yz = y * z;
        const double wx = w * x, wy = w * y, wz = w * z;
        r(0, 0) = 1 - 2 * (yy + zz); r(0, 1) = 2 * (xy - wz);     r(0, 2) = 2 * (xz + wy);
        r(1, 0) = 2 * (xy + wz);     r(1, 1) = 1 - 2 * (xx + zz); r(1, 2) = 2 * (yz - wx);
        r(2, 0) = 2 * (xz - wy);     r(2, 1) = 2 * (yz + wx);     r(2, 2) = 1 - 2 * (xx + yy);
        return r;
    }

    // Integrate this orientation by angular velocity omega (world frame) over dt.
    Quat integrated(const Vec3& omega, double dt) const {
        const Quat dq = Quat{0, omega.x * 0.5 * dt, omega.y * 0.5 * dt, omega.z * 0.5 * dt} * (*this);
        Quat q{w + dq.w, x + dq.x, y + dq.y, z + dq.z};
        q.normalizeInPlace();
        return q;
    }
};

// Heel (roll about x) and trim (pitch about y) in radians, extracted from a body->world rotation.
// Heel is positive to starboard; trim is positive bow-down.
inline void heelTrimFromRotation(const Mat3& R, double& heelRad, double& trimRad) {
    const Vec3 up = R * Vec3{0, 0, 1};       // body up, expressed in world
    const Vec3 fwd = R * Vec3{1, 0, 0};      // body forward, in world
    heelRad = std::atan2(up.y, up.z);        // +y is port, so port-down => positive; negate below
    heelRad = -heelRad;                      // now positive = starboard down
    trimRad = -std::asin(std::clamp(fwd.z, -1.0, 1.0));  // positive = bow down
}

}  // namespace sim
