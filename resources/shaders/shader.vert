#version 330 core

layout(location = 0) in vec3 position;
layout(location = 1) in vec3 normal;
layout(location = 2) in vec2 texCoord;

uniform mat4 proj;
uniform mat4 view;
uniform mat4 model;

uniform mat3 inverseTransposeModel;

out vec3 normal_worldSpace;
out vec3 position_worldSpace;
out vec2 fragTexCoord;

void main() {
    normal_worldSpace   = normalize(inverseTransposeModel * normal);
    position_worldSpace = position;
    fragTexCoord        = texCoord;

    gl_Position = proj * view * model * vec4(position, 1.0);
}
