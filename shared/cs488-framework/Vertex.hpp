#pragma once

#include <glm/glm.hpp>

#include <functional>

using namespace glm;

struct Vertex {
    glm::vec3 position;
    glm::vec3 normal;
    glm::vec2 uv;

    Vertex() : position(0.0f), normal(0.0f), uv(0.0f) {}

    Vertex(const glm::vec3& position, const glm::vec3& normal, const glm::vec2& uv)
        : position(position), normal(normal), uv(uv) {}

    bool operator==(const Vertex& other) const {
        return position == other.position;
    }

    struct Hash {
        std::size_t operator()(const Vertex& vertex) const {
            std::size_t h1 = std::hash<float>()(vertex.position.x) ^ 
                             (std::hash<float>()(vertex.position.y) << 1) ^ 
                             (std::hash<float>()(vertex.position.z) << 2);
            return h1;
        }
    };
};


// Specialization of std::hash for Vertex
namespace std {
    template <>
    struct hash<Vertex> {
        std::size_t operator()(const Vertex& vertex) const {
            return Vertex::Hash()(vertex);
        }
    };

    template <>
    struct hash<glm::vec3> {
        std::size_t operator()(const glm::vec3& v) const noexcept {
            return std::hash<float>()(v.x) ^ (std::hash<float>()(v.y) << 1) ^ (std::hash<float>()(v.z) << 2);
        }
    };
}

struct pair_hash {
    template <typename T1, typename T2>
    std::size_t operator()(const std::pair<T1, T2>& pair) const {
        T1 smaller = std::min(pair.first, pair.second);
        T2 larger = std::max(pair.first, pair.second);
        return std::hash<T1>()(smaller) ^ (std::hash<T2>()(larger) << 1);
    }
};
