#version 450
// Two homogeneous slabs of an emitting, absorbing, non-scattering grey gas.
//
// Every line of this is written out again in double in engine/gpu/smoke.cpp and
// asserted against from tests/test_smoke_render.cpp, which predicts pixels from
// that file rather than from this one. Keep the two in step: the value of the
// arrangement is that they disagree loudly when one drifts. The same trade
// material.hpp and hull.frag already have.
//
// **There is no ray march.** A two-zone model is exactly uniform inside each
// layer, so the optical depth along a ray is `k` times a path length and the path
// length is the analytic intersection of a segment with a half space. Marching
// would be sampling a function that is piecewise constant with one known
// breakpoint.

layout(location = 0) in vec3 fragWorld;
layout(location = 1) flat in int fragVolume;

layout(location = 0) out vec4 outColour;

struct Volume {
    vec4 loInterface;  // xyz body-frame low corner, w the layer interface height
    vec4 hi;           // xyz body-frame high corner, w unused
    vec4 row0;         // xyz row 0 of the body -> world rotation, w translation.x
    vec4 row1;
    vec4 row2;
    vec4 upper;        // rgb emitted radiance, a extinction coefficient
    vec4 lower;
};

layout(std430, set = 0, binding = 0) readonly buffer Volumes {
    Volume items[];
} volumes;

// The lit pass's depth attachment, sampled rather than tested against. A depth
// *test* would reject the whole fragment where the hull is nearer than the far
// face of the box; what is wanted is for the ray to stop at the hull, which is a
// clamp on the integration limit and not a rejection.
layout(set = 0, binding = 1) uniform sampler2D opaqueDepth;

layout(push_constant) uniform Push {
    mat4 modelViewProjection;
    vec4 eyeA;      // xyz world eye position, w the depth constant a
    vec4 forwardB;  // xyz unit view direction, w the depth constant b
} push;

// A direction component a reciprocal can be taken of. engine/gpu/smoke.cpp
// substitutes the same constant, so the two agree about a ray exactly parallel to
// a face instead of each being separately plausible. An infinity here would meet
// a zero difference and produce a NaN.
const float kParallel = 1e-30;

float safeInverse(float d) {
    return 1.0 / (abs(d) < kParallel ? (d < 0.0 ? -kParallel : kParallel) : d);
}

void main() {
    Volume v = volumes.items[fragVolume];

    const vec3 eye = push.eyeA.xyz;
    const vec3 forward = push.forwardB.xyz;
    const vec3 direction = normalize(fragWorld - eye);

    // Where the ray stops: the lit pass's own depth, turned back into a distance
    // along the view axis and then into a distance along *this* ray. An empty
    // pixel holds the clear depth of 1.0, which inverts to the far plane -- to
    // about eight parts in a million, because `depth + a` is a cancellation there
    // and the constants arrive as float32. It is the worst-conditioned point of
    // the inversion and it is also the one place nothing is drawn.
    const float depth = texelFetch(opaqueDepth, ivec2(gl_FragCoord.xy), 0).x;
    const float axial = push.forwardB.w / (depth + push.eyeA.w);
    const float along = max(dot(direction, forward), 1e-6);
    const float maxDistance = axial / along;

    // Into the body frame. R is stored by rows, so its transpose is the matrix
    // whose *columns* are those rows -- and because R is a rotation and the
    // direction is a unit vector, the ray parameter is metres in both frames.
    const mat3 bodyFromWorld = mat3(v.row0.xyz, v.row1.xyz, v.row2.xyz);
    const vec3 translation = vec3(v.row0.w, v.row1.w, v.row2.w);
    const vec3 origin = bodyFromWorld * (eye - translation);
    const vec3 ray = bodyFromWorld * direction;

    const vec3 inverse = vec3(safeInverse(ray.x), safeInverse(ray.y), safeInverse(ray.z));
    const vec3 a = (v.loInterface.xyz - origin) * inverse;
    const vec3 b = (v.hi.xyz - origin) * inverse;
    const vec3 entry = min(a, b);
    const vec3 exit = max(a, b);
    const float tEnter = max(max(entry.x, entry.y), max(entry.z, 0.0));
    const float tExit = min(min(exit.x, exit.y), min(exit.z, maxDistance));
    if (tExit <= tEnter) discard;

    const float span = tExit - tEnter;

    // The layer split. z is monotone along the ray, so each layer is one
    // contiguous interval and the order they are met in is the sign of ray.z.
    float upperPath;
    bool upperFirst = true;
    if (abs(ray.z) < kParallel) {
        upperPath = origin.z >= v.loInterface.w ? span : 0.0;
    } else {
        const float crossing = (v.loInterface.w - origin.z) / ray.z;
        if (ray.z < 0.0) {
            upperFirst = true;
            upperPath = clamp(min(tExit, crossing), tEnter, tExit) - tEnter;
        } else {
            upperFirst = false;
            upperPath = tExit - clamp(max(tEnter, crossing), tEnter, tExit);
        }
    }
    upperPath = clamp(upperPath, 0.0, span);
    const float lowerPath = span - upperPath;

    // Beer-Lambert, and Kirchhoff on the same exponential: a slab that transmits
    // everything emits nothing, however hot it is.
    const float tUpper = exp(-v.upper.a * upperPath);
    const float tLower = exp(-v.lower.a * lowerPath);
    const vec3 eUpper = v.upper.rgb * (1.0 - tUpper);
    const vec3 eLower = v.lower.rgb * (1.0 - tLower);
    const vec3 source = upperFirst ? eUpper + tUpper * eLower : eLower + tLower * eUpper;

    // Premultiplied alpha, blended ONE / ONE_MINUS_SRC_ALPHA, so the destination
    // is multiplied by exactly `exp(-k_u d_u) exp(-k_l d_l)` and the transmittance
    // is the hardware's arithmetic rather than a second copy of it here. At zero
    // extinction that factor is exactly 1.0 and the source is exactly 0.0, which
    // is what makes an unsmoked frame bit-identical rather than merely close.
    outColour = vec4(source, 1.0 - tUpper * tLower);
}
