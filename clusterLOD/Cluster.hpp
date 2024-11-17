#pragma once

#include <iostream>
#include <vector>
#include <memory>
#include "cs488-framework/GlErrorCheck.hpp"
#include "Mesh.hpp"

#define MAX_TRI_IN_CLUSTER 256
#define MAXN_CLUSTER_IN_CLUSTERGROUP 32

class Cluster : public Mesh {
public:
    float Error;
    
    glm::vec3 rdColor;

    Cluster(
        float Error, 
        const Mesh& ref, 
        const std::vector<unsigned int>& triIndices
    );

    ~Cluster() {}

	void draw(const ShaderProgram& shader) const override;
};

// class ClusterGroup {
// public:
//     float Error;
//     std::vector<std::shared_ptr<Cluster>> clusters;

//     void addCluster(std::shared_ptr<Cluster> cluster) {
//         clusters.push_back(cluster);
//     }
// };

void MeshSplitter(Mesh& mesh);

std::vector<std::vector<size_t>> 
BuildAdjacencyList(const std::vector<unsigned int>& m_indexData);