#pragma once

#include <glm/glm.hpp>
using namespace glm;

struct Vertex {
    glm::vec3 position;
    glm::vec3 normal;
    glm::vec2 uv;

    bool operator==(const Vertex& other) const {
        return position == other.position && normal == other.normal && uv == other.uv;
    }
    
    Vertex(const glm::vec3& position, const glm::vec3& normal, const glm::vec2& uv) 
        : position(position), normal(normal), uv(uv) {}
};

struct VertexHash {
    std::size_t operator()(const Vertex& vertex) const {
        std::size_t h1 = std::hash<float>()(vertex.position.x) ^ (std::hash<float>()(vertex.position.y) << 1) ^ (std::hash<float>()(vertex.position.z) << 2);
        std::size_t h2 = std::hash<float>()(vertex.normal.x) ^ (std::hash<float>()(vertex.normal.y) << 1) ^ (std::hash<float>()(vertex.normal.z) << 2);
        std::size_t h3 = std::hash<float>()(vertex.uv.x) ^ (std::hash<float>()(vertex.uv.y) << 1);
        return h1 ^ (h2 << 1) ^ (h3 << 2); // Combine the hashes
    }
};