#version 330 core

out vec4 FragColor;
in vec2 TexCoords;
in vec3 TransformedLightPositions[20];

uniform sampler2D gPosition;
uniform sampler2D gNormal;
uniform sampler2D gAlbedoID;

struct Light {
    vec3 Color;
};

uniform int numLights;
uniform Light lights[256];

void main() {
    vec4 FragPos = texture(gPosition, TexCoords).rgba;
    vec4 Normal = texture(gNormal, TexCoords).rgba;
    vec4 Albedo = texture(gAlbedoID, TexCoords).rgba;

    // depth testing
    if(FragPos.a == 1.0){
        FragColor = vec4(0.0);
        return;
    }
    else{
        FragColor = vec4(Albedo.rgb, 1.0);
        return;
    }

    vec3 lighting = vec3(0.0);

    for (int i = 0; i < numLights; ++i) {
        vec3 lightDir = normalize(TransformedLightPositions[i] - FragPos.rgb);
        float diff = max(dot(Normal.rgb, lightDir), 0.0);
        vec3 diffuse = diff * lights[i].Color;

        lighting += diffuse * Albedo.rgb;
    }

    FragColor = vec4(lighting, 1.0);
}
