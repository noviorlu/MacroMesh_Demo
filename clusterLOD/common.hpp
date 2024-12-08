
#pragma once
#include <glm/glm.hpp>
#include <functional>

// Specialization of std::hash for Vertex
namespace std {
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

