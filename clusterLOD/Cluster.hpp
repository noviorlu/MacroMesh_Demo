#pragma once

#include <iostream>
#include <vector>
#include <memory>
#include "cs488-framework/GlErrorCheck.hpp"
#include "Mesh.hpp"
#include <metis.h>

#define MAX_TRI_IN_CLUSTER 256
#define MAXN_CLUSTER_IN_CLUSTERGROUP 32


class Cluster : public Mesh {
public:
    float Error;
    
    glm::vec3 rdColor;
    std::vector<Cluster*> adjacent_clusters;

    Cluster(
        float Error, 
        const Mesh& ref, 
        const std::vector<unsigned int>& triIndices
    );

    ~Cluster() {}

	void draw(const ShaderProgram& shader) const override;
};

class ClusterGroup : public Mesh {
public:
    float Error;
    
    std::vector<Cluster*> clusters;
};

void MeshSplitter(Mesh& mesh);

void clusterGrouping(const Mesh& mesh_ref, std::vector<ClusterGroup*>& cluster_groups);

std::vector<std::vector<size_t>> 
BuildAdjacencyList(const std::vector<unsigned int>& m_indexData);