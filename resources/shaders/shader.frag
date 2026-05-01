#version 330 core
out vec4 fragColor;

// Additional information for lighting
in vec3 normal_worldSpace;
in vec3 position_worldSpace;

uniform int wire = 0;
uniform float red = 1.0;
uniform float green = 1.0;
uniform float blue = 1.0;
uniform float alpha = 1.0;
uniform mat4 view;

void main() {
    if (wire == 1) {
        fragColor = vec4(0.0, 0.0, 0.0, 1);
        return;
    }

    vec3 N = normalize(normal_worldSpace);
    vec3 camPos = -mat3(view) * view[3].xyz;
    vec3 V = normalize(camPos - position_worldSpace);

    // Material color
    vec3 baseColor = vec3(red, green, blue);

    // Warm light
    vec3 warmLightDir = normalize(vec3(2.0, 4.0, -3.0));
    vec3 warmLightColor = vec3(1.0, 0.95, 0.8);
    float warmDiffuse = max(dot(N, warmLightDir), 0.0);
    float warmSpec = pow(max(dot(N, normalize(warmLightDir + V)), 0.0), 64.0);

    // Cool light
    vec3 coolLightDir = normalize(vec3(-1.0, 0.5, 2.0));
    vec3 coolLightColor = vec3(0.3, 0.55, 0.75);
    float coolDiffuse = max(dot(N, coolLightDir), 0.0);

    // Arbitrary blue-ish ambient color
    vec3 ambient = vec3(0.1, 0.1, 0.16);

    // Add it all together
    vec3 color = baseColor * (ambient
                 + warmDiffuse * warmLightColor * 0.7
                 + coolDiffuse * coolLightColor * 0.3)
                 + warmSpec * warmLightColor * 0.3;

    // Dampen the color a bit
    color = color / (color + vec3(0.5));

    // Distance fog
    float dist = length(camPos - position_worldSpace);
    vec3 fogColor = vec3(0.05, 0.05, 0.1);
    float fogFactor = clamp(1.0 - exp(-0.005 * dist * dist), 0.0, 1.0);
    color = mix(color, fogColor, fogFactor);

    fragColor = vec4(color, 1.0);
}
