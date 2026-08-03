#version 450
// The lit solid: hull, and anything else drawn into the same depth buffer.
//
// Positions arrive in **world space** -- engine/gpu/hull.cpp bakes the ship's
// rigid-body transform into the vertices, because the sea is world space and one
// depth buffer means one space. There is no model matrix and no per-object
// descriptor.
//
// The material is a flat integer index into a storage buffer, not an
// interpolated colour: a mod adds a material by adding a row to that buffer, so
// nothing here or in the fragment stage enumerates them.
//
// Pascal target: a plain vertex stage over indexed triangles. No mesh shaders.

layout(location = 0) in vec3 inPosition;
layout(location = 1) in vec3 inNormal;
layout(location = 2) in uint inMaterial;

layout(location = 0) out vec3 fragNormal;
layout(location = 1) out vec3 fragWorld;
layout(location = 2) flat out uint fragMaterial;

layout(push_constant) uniform Push {
    mat4 modelViewProjection;
    vec4 sun;        // xyz unit vector toward the sun, w mode
    vec4 sunColour;  // rgb radiance, w exposure
    vec4 sky;        // rgb hemispheric ambient, w unused
    vec4 eye;        // xyz world eye position, w unused
} push;

void main() {
    gl_Position = push.modelViewProjection * vec4(inPosition, 1.0);
    fragNormal = inNormal;
    // The eye vector is per fragment, not per vertex: the specular lobe is narrow
    // enough that interpolating it across a 60 m hull plate is visible.
    fragWorld = inPosition;
    fragMaterial = inMaterial;
}
