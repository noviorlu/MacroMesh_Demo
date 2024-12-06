#pragma once
#include <glm/glm.hpp>

#include <queue>
#include <unordered_set>
#include <vector>
#include <array>
#include <algorithm>
#include <string>

#include <chrono>

class SimpleVertex{
public:
    glm::vec3 position;
    glm::vec3 normal;
    glm::vec2 uv;

    SimpleVertex(glm::vec3 position, glm::vec3 normal, glm::vec2 uv) : position(position), normal(normal), uv(uv) {}
};

class SimpleFace{
public: 
    SimpleVertex* v1;
    SimpleVertex* v2;
    SimpleVertex* v3;
    std::vector<SimpleFace*> adjacentFace;
    SimpleFace(SimpleVertex* v1, SimpleVertex* v2, SimpleVertex* v3) : v1(v1), v2(v2), v3(v3) {}

    int clusterId;

	static void swap(SimpleFace& a, SimpleFace& b) {
		std::swap(a.clusterId, b.clusterId);
        
        std::swap(a.v1, b.v1);
		std::swap(a.v2, b.v2);
		std::swap(a.v3, b.v3);

        for (auto* AadjFace : a.adjacentFace) {
            for (int i = 0; i < AadjFace->adjacentFace.size(); i++) {
                if (AadjFace->adjacentFace[i] == &a) {
                    AadjFace->adjacentFace[i] = &b;
                }
            }
        }
        
        for (auto* BadjFace : b.adjacentFace) {
            for (int i = 0; i < BadjFace->adjacentFace.size(); i++) {
                if (BadjFace->adjacentFace[i] == &b) {
                    BadjFace->adjacentFace[i] = &a;
                }
            }
        }

		std::swap(a.adjacentFace, b.adjacentFace);
	}
};


class SimpleMesh {
public: 
    std::string m_name;
    std::vector<SimpleVertex> m_vertices;
    std::vector<SimpleFace> m_faces;

    void importMesh(const std::string& objFilePath);
    void exportMesh(int startIdx, int endIdx, const std::string& objFilePath);
    void exportMesh(const std::string& objFilePath, int clusterId);
    void exportClusterGroup(const std::string& lodFolderPath);

    void splitterRecur(int start, int end, int depth);
    void splitter();

    void partition_loop(const std::string& objFilePath, const std::string& lodFolderPath);

    std::vector<int> m_clusterGroupOffsets;
};

