
#include "Mesh.hpp"

Cluster::Cluster(float Error) : Error(Error) {
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
    }
}