
#include "Mesh.hpp"
glm::vec3 HSVtoRGB(float h, float s, float v) {
    h = fmod(h, 1.0f) * 6.0f;
    int i = static_cast<int>(floor(h));
    float f = h - i;
    float p = v * (1.0f - s);
    float q = v * (1.0f - f * s);
    float t = v * (1.0f - (1.0f - f) * s);

    switch (i) {
    case 0: return glm::vec3(v, t, p);
    case 1: return glm::vec3(q, v, p);
    case 2: return glm::vec3(p, v, t);
    case 3: return glm::vec3(p, q, v);
    case 4: return glm::vec3(t, p, v);
    case 5: return glm::vec3(v, p, q);
    default: return glm::vec3(0.0f, 0.0f, 0.0f);
    }
}

glm::vec3 genRdColor(){
    // generate a random hue (h) in [0, 1)
    float h = static_cast<float>(rand()) / RAND_MAX;
    float s = 0.8f;
    float v = 0.8f;

    // convert HSV to RGB
    return HSVtoRGB(h, s, v);
}

Cluster::Cluster(float Error) : Error(Error) {
    rdColor = genRdColor();
}

void Cluster::draw(const ShaderProgram& shader) const {
    shader.SetUniform3fv("material.kd", rdColor);
    Mesh::draw(shader);
}

ClusterGroup::ClusterGroup(float Error) : Error(Error) {}

void ClusterGroup::setGroupColor() {
    float h = static_cast<float>(rand()) / RAND_MAX;
    float s = 0.8f;
    for (auto cluster : clusters) {
        float v = 0.5f + static_cast<float>(rand()) / (2.0f * RAND_MAX); // v in range [0.5, 1.0)
        cluster->rdColor = HSVtoRGB(h, s, v);
    }
}