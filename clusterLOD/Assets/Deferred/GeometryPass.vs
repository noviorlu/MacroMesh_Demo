#version 330 core

layout(location = 0) in vec3 inPosition; 
layout(location = 1) in vec3 inNormal;   
layout(location = 2) in vec2 inTexCoords;

out VS_OUT {
    vec3 FragPos;
    vec3 Normal;
    vec2 TexCoords;
} vs_out;

uniform mat4 ModelView;
uniform mat4 Perspective;
uniform mat3 NormalMatrix;

void main() {
    vec4 pos4 = vec4(inPosition, 1.0);
    vs_out.FragPos = (ModelView * pos4).xyz;
    vs_out.Normal = normalize(NormalMatrix * inNormal);
    vs_out.TexCoords = inTexCoords;
    gl_Position = Perspective * ModelView * vec4(inPosition, 1.0);
}
