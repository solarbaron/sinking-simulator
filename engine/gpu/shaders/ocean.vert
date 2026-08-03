#version 450
// The sea surface grid.
//
// The vertices arrive already displaced: engine/gpu/ocean.cpp evaluates the same
// sim::WaveField the physics reads, so there is no second wave model here that
// could drift from it. Positions are world space and there is no model matrix --
// the sea does not move, the ship does.
//
// Pascal target: a plain vertex stage over an indexed grid. No mesh shaders.

layout(location = 0) in vec3 inPosition;
layout(location = 1) in vec3 inNormal;

layout(location = 0) out vec3 fragNormal;
layout(location = 1) out float fragHeight;

layout(push_constant) uniform Push {
    mat4 modelViewProjection;
    vec4 sun;      // xyz unit vector toward the sun, w strength
    vec4 water;    // rgb water colour, w ambient
    vec4 encode;   // x mode, y zMin, z 1 / span, w unused
} push;

void main() {
    gl_Position = push.modelViewProjection * vec4(inPosition, 1.0);
    fragNormal = inNormal;
    // World elevation, so the fragment stage can report the geometry it is
    // actually shading. Perspective-correct interpolation makes this exactly the
    // piecewise-linear mesh height at the point the fragment covers.
    fragHeight = inPosition.z;
}
