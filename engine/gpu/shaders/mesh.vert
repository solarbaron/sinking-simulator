#version 450
// Flat-shaded indexed triangles. The whole transform arrives as one push-constant
// matrix, so a draw needs no descriptor set at all -- which keeps the debug
// renderer independent of the bindless work that comes later.

layout(location = 0) in vec3 inPosition;
layout(location = 1) in vec3 inColour;

layout(location = 0) out vec3 fragColour;

layout(push_constant) uniform Push {
    mat4 modelViewProjection;
} push;

void main() {
    gl_Position = push.modelViewProjection * vec4(inPosition, 1.0);
    fragColour = inColour;
}
