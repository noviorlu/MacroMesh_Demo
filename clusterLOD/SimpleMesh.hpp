#pragma once
#include <glm/glm.hpp>

#include <unordered_set>
#include <vector>
#include <algorithm>
#include <string>

#include <omp.h>

class MeshSimplifier;

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
    std::unordered_set<SimpleFace*> adjacentFace;
    SimpleFace(SimpleVertex* v1, SimpleVertex* v2, SimpleVertex* v3) : v1(v1), v2(v2), v3(v3) {
        omp_init_lock(&lock);
		adjacentFace = std::unordered_set<SimpleFace*>();
    }

    int clusterId;
    omp_lock_t lock;

	static void swap(SimpleFace& a, SimpleFace& b) {
		std::swap(a.clusterId, b.clusterId);
        
        std::swap(a.v1, b.v1);
		std::swap(a.v2, b.v2);
		std::swap(a.v3, b.v3);

        for (SimpleFace* AadjFace : a.adjacentFace) {
            omp_set_lock(&AadjFace->lock);
            {
			std::unordered_set<SimpleFace*>& adjAadj = AadjFace->adjacentFace;
            adjAadj.erase(&a);
            adjAadj.insert(&b);
            }
            omp_unset_lock(&AadjFace->lock);
        }
        
        for (auto* BadjFace : b.adjacentFace) {
            omp_set_lock(&BadjFace->lock);
            {
            std::unordered_set<SimpleFace*>& adjBadj = BadjFace->adjacentFace;
            adjBadj.erase(&b);
            adjBadj.insert(&a);
            }
            omp_unset_lock(&BadjFace->lock);
        }

		std::swap(a.adjacentFace, b.adjacentFace);
	}
};


class SimpleMesh {
public: 
    std::string m_name;
    std::vector<SimpleVertex> m_vertices;
    std::vector<SimpleFace> m_faces;

    std::vector<std::pair<SimpleVertex*, SimpleVertex*>> m_boundary;

    void importMesh(const std::string& objFilePath);
    void exportMesh(int startIdx, int endIdx, const std::string& objFilePath);
    void exportMesh(const std::string& objFilePath, int clusterId);
    void exportClusterGroup(const std::string& lodFolderPath);

    void splitterRecur(int start, int end, int depth);
    void splitter();
    void grouper();

    void exportMeshSimplifier(MeshSimplifier& simplifier, int startIdx, int endIdx);
    float QEM(int start, int end, const std::string& lodFolderPath, float ratio);
    void importMeshSimplifier(const MeshSimplifier& simplifier);

    void partition_loop(const std::string& objFilePath, const std::string& lodFolderPath);

    std::vector<int> m_clusterOffsets;
    std::vector<int> m_clusterGroupOffsets;
};

