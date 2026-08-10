#version 450
// The prism a two-zone gas space occupies, as seen from wherever the camera is.
//
// **No vertex buffer.** The box is generated from gl_VertexIndex against the
// volume's own body-frame corners, so the only geometry the volumetric pass
// uploads is the volume table itself -- seven vec4s per burning compartment,
// against 36 vertices of nothing.
//
// The rasterizer culls the faces whose outward normal points at the eye, so each
// covered pixel receives exactly **one** fragment, from the far side of the box.
// That is what makes the composite in smoke.frag a single application of
// Beer-Lambert rather than one per face -- and it keeps working when the camera is
// inside the box, where the far faces are the only ones in front of the near
// plane at all.
//
// Pascal target: a plain vertex stage. No mesh shaders.

layout(location = 0) out vec3 fragWorld;
layout(location = 1) flat out int fragVolume;

// Mirrors gpu::GpuSmokeVolume in engine/gpu/smoke_gpu.cpp: seven vec4s, 112
// bytes, which std430 packs tightly.
struct Volume {
    vec4 loInterface;  // xyz body-frame low corner, w the layer interface height
    vec4 hi;           // xyz body-frame high corner, w unused
    vec4 row0;         // xyz row 0 of the body -> world rotation, w translation.x
    vec4 row1;         // ...                  1,                 w translation.y
    vec4 row2;         // ...                  2,                 w translation.z
    vec4 upper;        // rgb emitted radiance, a extinction coefficient
    vec4 lower;
};

layout(std430, set = 0, binding = 0) readonly buffer Volumes {
    Volume items[];
} volumes;

layout(push_constant) uniform Push {
    mat4 modelViewProjection;
    vec4 eyeA;      // xyz world eye position, w the depth constant a
    vec4 forwardB;  // xyz unit view direction, w the depth constant b
} push;

// Six quads as twelve triangles, each wound counter-clockwise seen from outside
// the box. Corner c takes hi on axis i where bit i of c is set.
const int kCorner[36] = int[36](
    0, 4, 6,  0, 6, 2,      // -x
    1, 3, 7,  1, 7, 5,      // +x
    0, 1, 5,  0, 5, 4,      // -y
    2, 6, 7,  2, 7, 3,      // +y
    0, 2, 3,  0, 3, 1,      // -z
    4, 5, 7,  4, 7, 6);     // +z

void main() {
    Volume v = volumes.items[gl_InstanceIndex];
    int c = kCorner[gl_VertexIndex];
    vec3 body = vec3((c & 1) != 0 ? v.hi.x : v.loInterface.x,
                     (c & 2) != 0 ? v.hi.y : v.loInterface.y,
                     (c & 4) != 0 ? v.hi.z : v.loInterface.z);
    // R is stored by rows, so the body -> world product is three dot products.
    vec3 world = vec3(dot(v.row0.xyz, body), dot(v.row1.xyz, body), dot(v.row2.xyz, body)) +
                 vec3(v.row0.w, v.row1.w, v.row2.w);

    fragWorld = world;
    fragVolume = gl_InstanceIndex;
    gl_Position = push.modelViewProjection * vec4(world, 1.0);
}
