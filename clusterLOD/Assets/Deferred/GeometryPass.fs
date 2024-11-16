#version 330 core

layout (location = 0) out vec4 gPosition;  
layout (location = 1) out vec4 geoNormal;
layout (location = 2) out vec4 gAlbedoID;
in VS_OUT {
    vec3 FragPos;
    vec3 Normal;
    vec2 TexCoords;
} fs_in;

struct Material {
    vec3 kd;
};
uniform Material material;
uniform int nodeId;

void main() {
    gPosition = vec4(fs_in.FragPos, gl_FragCoord.z);
    geoNormal = vec4(normalize(fs_in.Normal), nodeId);
    gAlbedoID = vec4(material.kd, 1.0);

    // float clr = float(fs_in.TexCoords.x / 255.0);
    // gAlbedoID = vec4(clr, clr, clr, 1.0);
}
