#version 330 core

layout(location = 0) in vec3 inPosition; 
layout(location = 1) in vec3 inNormal;   

out VS_OUT {
    vec3 FragPos;
    vec3 Normal;
} vs_out;

uniform mat4 ModelView;
uniform mat4 Perspective;
uniform mat3 NormalMatrix;

void main() {
    vec4 pos4 = vec4(inPosition, 1.0);
    vs_out.FragPos = (ModelView * pos4).xyz;
    vs_out.Normal = normalize(NormalMatrix * inNormal);
    gl_Position = Perspective * ModelView * vec4(inPosition, 1.0);
}
