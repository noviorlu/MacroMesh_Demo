#pragma once

#include <glm/glm.hpp>
#include <functional>

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
        std::size_t operator()(const glm::vec3& vec) const {
            auto hashX = std::hash<float>{}(vec.x);
            auto hashY = std::hash<float>{}(vec.y);
            auto hashZ = std::hash<float>{}(vec.z);
            return hashX ^ (hashY << 1) ^ (hashZ << 2);
        }
    };

    template <>
    struct equal_to<glm::vec3> {
        bool operator()(const glm::vec3& lhs, const glm::vec3& rhs) const {
            return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z;
        }
    };

    template <typename T1, typename T2>
    struct hash<std::pair<T1, T2>> {
        std::size_t operator()(const std::pair<T1, T2>& pair) const {
            return std::hash<T1>()(pair.first) ^ (std::hash<T2>()(pair.second) << 1);
        }
    };
}

