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
    unsigned int clusterId = 0;
    unsigned int sequenceId = 0;

    SimpleVertex* v1;
    SimpleVertex* v2;
    SimpleVertex* v3;
    std::unordered_set<SimpleFace*> adjacentFace;
    SimpleFace(SimpleVertex* v1, SimpleVertex* v2, SimpleVertex* v3, unsigned int sequenceId) : v1(v1), v2(v2), v3(v3), sequenceId(sequenceId) {
		adjacentFace = std::unordered_set<SimpleFace*>();
    }
};


class SimpleMesh {
public: 
    std::string m_name;
    std::vector<SimpleVertex> m_vertices;
    std::vector<SimpleFace*> m_faces;

    std::vector<std::pair<SimpleVertex*, SimpleVertex*>> m_boundary;

    void importMesh(const std::string& objFilePath);
    void exportMesh(int startIdx, int endIdx, const std::string& objFilePath);
    void exportMesh(const std::vector<std::pair<int, int>>& indexRanges, const std::string& objFilePath);
    void exportMesh(const std::string& objFilePath, int clusterId);
    void exportCluster(const std::string& lodFolderPath);
    void exportClusterGroup(const std::string& lodFolderPath);

    void splitterRecur(unsigned int start, unsigned int end, int depth);
    void splitter();
    void grouperRecur(unsigned int start, unsigned int end, int depth);
    void grouper();

    void exportMeshSimplifier(MeshSimplifier& simplifier, int startIdx, int endIdx);
    float QEM(int start, int end, const std::string& lodFolderPath, float ratio);
    void importMeshSimplifier(const MeshSimplifier& simplifier);

    void partition_loop(const std::string& objFilePath, const std::string& lodFolderPath);

    std::vector<unsigned int> m_clusterOffsets;
    class Cluster {
    public:
        unsigned int startIdx;
        unsigned int endIdx;

        unsigned int sequenceId;

        std::unordered_set<Cluster*> adjClusters;

        Cluster(unsigned int startIdx, unsigned int endIdx, unsigned int sequenceId) : startIdx(startIdx), endIdx(endIdx), sequenceId(sequenceId) {}
    };
    std::vector<Cluster*> m_clusters;

    std::vector<unsigned int> m_clusterGroupOffsets;
};

