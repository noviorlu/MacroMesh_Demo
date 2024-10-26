#version 330 core

layout (location = 0) in vec2 inPosition;
layout (location = 1) in vec2 inTexCoords;

out vec2 TexCoords;
out vec3 TransformedLightPositions[20];

uniform mat4 View;

uniform int numLights;
uniform vec3 lightPositions[20];

void main() {
    TexCoords = inTexCoords;
    
    for (int i = 0; i < numLights; ++i) {
        vec4 transformedPos = View * vec4(lightPositions[i], 1.0);
        TransformedLightPositions[i] = transformedPos.xyz / transformedPos.w;
    }

    gl_Position = vec4(inPosition, 0.0, 1.0);
}
