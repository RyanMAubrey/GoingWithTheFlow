#version 330 core

layout(location = 0) in vec3 position; // Position of the vertex
layout(location = 1) in vec3 normal;   // Normal of the vertex

uniform mat4 proj;
uniform mat4 view;
uniform mat4 model;

uniform mat3 inverseTransposeModel;

out vec3 normal_worldSpace;
out vec3 position_worldSpace;

void main() {
    normal_worldSpace   = normalize(inverseTransposeModel * normal);
    position_worldSpace = position;

    gl_Position = proj * view * model * vec4(position, 1.0);
}
