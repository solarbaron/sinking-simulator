#version 450
// Metallic-roughness shading, plus two readback channels.
//
// The BRDF is written out in full in engine/gpu/material.hpp and is repeated
// independently in tests/test_hull_render.cpp, which predicts pixels from it
// rather than asking this shader what it thinks. Keep the three in step: the
// value of the arrangement is that two of them are written from the formula and
// disagree loudly when one drifts.

layout(location = 0) in vec3 fragNormal;
layout(location = 1) in vec3 fragWorld;
layout(location = 2) flat in uint fragMaterial;

layout(location = 0) out vec4 outColour;

// std430 packs two vec4s tightly, so this row is 32 bytes and matches
// gpu::GpuMaterial byte for byte.
struct Material {
    vec4 baseColour;  // rgb linear reflectance, w opacity (not yet rendered)
    vec4 params;      // x roughness, y metalness, zw reserved
};

layout(std430, set = 0, binding = 0) readonly buffer Materials {
    Material items[];
} materials;

layout(push_constant) uniform Push {
    mat4 modelViewProjection;
    vec4 sun;        // xyz unit vector toward the sun, w mode
    vec4 sunColour;  // rgb radiance, w exposure
    vec4 sky;        // rgb hemispheric ambient, w unused
    vec4 eye;        // xyz world eye position, w unused
} push;

const float kPi = 3.141592653589793;

// 16-bit code, high byte in red and low byte in green, blue as a tag so a reader
// can tell drawn geometry from the background without knowing the clear colour.
// Both bytes survive the UNORM8 store exactly: the store is round(v * 255) and
// n / 255.0 round-trips for every integer n in [0, 255].
vec4 encode16(float code) {
    float high = floor(code * (1.0 / 256.0));
    return vec4(high / 255.0, (code - high * 256.0) / 255.0, 1.0, 1.0);
}

void main() {
    int mode = int(push.sun.w + 0.5);

    // Material channel. Not a debug colour: it turns "does the ship occlude the
    // sea behind it" into an exact integer question at a pixel sim::clipToPixel
    // picked in advance, rather than a comparison of two shaded colours that
    // might legitimately be close.
    if (mode == 1) {
        outColour = encode16(float(fragMaterial));
        return;
    }
    // Depth channel: the value the depth test actually used, so a test can hold
    // it against the clip-space z the camera matrix predicts.
    if (mode == 2) {
        outColour = encode16(floor(clamp(gl_FragCoord.z, 0.0, 1.0) * 65535.0 + 0.5));
        return;
    }

    Material m = materials.items[fragMaterial];

    vec3 n = normalize(fragNormal);
    // Nothing is culled, for the same reason the sea is not: a hull is seen from
    // inside once it is cut away or flooded, and a back face there is a real
    // surface rather than an error. The far side is lit as the face it presents.
    if (!gl_FrontFacing) n = -n;

    vec3 v = normalize(push.eye.xyz - fragWorld);
    vec3 l = push.sun.xyz;
    vec3 h = normalize(l + v);

    float ndl = max(dot(n, l), 0.0);
    float ndv = max(dot(n, v), 1e-4);
    float ndh = max(dot(n, h), 0.0);
    float vdh = max(dot(v, h), 0.0);

    float rough = clamp(m.params.x, 0.03, 1.0);
    float alpha = rough * rough;
    float a2 = alpha * alpha;

    // GGX / Trowbridge-Reitz.
    float denominator = ndh * ndh * (a2 - 1.0) + 1.0;
    float D = a2 / (kPi * denominator * denominator);

    // Smith height-correlated visibility, with the 1 / (4 n.l n.v) of the
    // microfacet denominator folded in -- which is also what keeps this finite as
    // n.l goes to zero, where the separate form is 0 / 0.
    float gv = ndl * sqrt(ndv * ndv * (1.0 - a2) + a2);
    float gl = ndv * sqrt(ndl * ndl * (1.0 - a2) + a2);
    float Vis = 0.5 / max(gv + gl, 1e-9);

    vec3 f0 = mix(vec3(0.04), m.baseColour.rgb, m.params.y);
    vec3 F = f0 + (1.0 - f0) * pow(1.0 - vdh, 5.0);

    // Exactly Lambert: no (1 - F) factor. That is what makes the diffuse term a
    // pure cosine law, which tests/test_hull_render.cpp sweeps rather than spot
    // checks. material.hpp records the trade.
    vec3 diffuse = m.baseColour.rgb * (1.0 - m.params.y) / kPi;
    vec3 specular = vec3(D * Vis) * F;

    // The same hemispheric sky term the sea uses, so a hull and the water beside
    // it agree about where the sky is.
    vec3 sky = m.baseColour.rgb * push.sky.rgb * (0.5 + 0.5 * n.z);

    vec3 radiance = sky + (diffuse + specular) * ndl * push.sunColour.rgb;
    outColour = vec4(clamp(radiance * push.sunColour.w, 0.0, 1.0), 1.0);
}
