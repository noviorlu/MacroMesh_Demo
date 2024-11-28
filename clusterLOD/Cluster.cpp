
#include "Mesh.hpp"

glm::vec3 HSVtoRGB(float h, float s, float v) {
    h = fmod(h, 1.0f) * 6.0f;  // 将 h 限制在 [0, 6)
    int i = static_cast<int>(floor(h));
    float f = h - i; // 小数部分
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
        default: return glm::vec3(0.0f, 0.0f, 0.0f); // 容错
    }
}

glm::vec3 genRdColor(){
    // generate a random hue (h) in [0, 1)
    float h = static_cast<float>(rand()) / RAND_MAX; // 随机色相
    float s = 0.8f; // 固定饱和度，接近鲜艳的颜色
    float v = 0.8f; // 固定亮度

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