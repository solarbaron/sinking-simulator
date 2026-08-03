#version 450

layout(location = 0) in vec3 fragNormal;
layout(location = 1) in float fragHeight;

layout(location = 0) out vec4 outColour;

layout(push_constant) uniform Push {
    mat4 modelViewProjection;
    vec4 sun;      // xyz unit vector toward the sun, w strength
    vec4 water;    // rgb water colour, w ambient
    vec4 encode;   // x mode (0 shaded, 1 elevation), y zMin, z 1 / span, w unused
} push;

void main() {
    if (push.encode.x > 0.5) {
        // Elevation channel: world z as a 16-bit code, high byte in red and low
        // byte in green, with blue as a surface tag so a reader can tell sea from
        // background without knowing the clear colour.
        //
        // Both bytes survive the UNORM8 store exactly. The store is
        // round(v * 255) and n / 255.0 round-trips for every integer n in
        // [0, 255], so what comes back is the code that went in -- which is what
        // lets tests/test_ocean.cpp hold rendered geometry against
        // WaveField::elevation() to a quarter of a millimetre rather than to a
        // colour ramp's eyeballed accuracy.
        float u = clamp((fragHeight - push.encode.y) * push.encode.z, 0.0, 1.0);
        float code = floor(u * 65535.0 + 0.5);
        float high = floor(code * (1.0 / 256.0));
        outColour = vec4(high / 255.0, (code - high * 256.0) / 255.0, 1.0, 1.0);
        return;
    }

    // Deliberately view-independent: Lambert against a directional sun plus a
    // hemispheric sky term, and no specular or Fresnel. A flat sea then has
    // exactly one colour over its whole surface, which is a closed form a test
    // can state and check to the last bit -- and the lighting still answers to
    // the waves, because the normal is the spectrum's own analytic slope.
    // Specular is what makes water look like water and it costs that assertion,
    // so it arrives with the reflection work rather than ahead of it.
    vec3 n = normalize(fragNormal);
    // The surface is seen from underneath too, from a flooded compartment or from
    // below the waterline, so nothing is culled and the far side is lit as the
    // face it actually presents.
    if (!gl_FrontFacing) n = -n;

    float ndotl = max(dot(n, push.sun.xyz), 0.0);
    float sky = 0.5 + 0.5 * n.z;
    vec3 lit = push.water.rgb * (push.water.w * sky + push.sun.w * ndotl);
    outColour = vec4(clamp(lit, 0.0, 1.0), 1.0);
}
