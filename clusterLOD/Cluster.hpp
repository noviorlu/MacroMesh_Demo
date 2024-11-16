#pragma once

#include <iostream>
#include <vector>
#include <memory>

#include "Mesh.hpp"

#define N 256

class Cluster : public std::enable_shared_from_this<Cluster> {
public:
    float Error;
    
    Mesh* m_mesh;

    Cluster(float Error, Mesh* mesh) : Error(Error), m_mesh(mesh) {}

    ~Cluster() {if(m_mesh != nullptr) delete m_mesh;}
};

class ClusterGroup {
public:
    float Error;
    std::vector<std::shared_ptr<Cluster>> clusters;

    void addCluster(std::shared_ptr<Cluster> cluster) {
        clusters.push_back(cluster);
    }
};

void MeshSplitter(Mesh& mesh, std::vector<Cluster>& clusters, int num_parts);

std::vector<std::vector<size_t>> 
BuildAdjacencyList(const std::vector<unsigned int>& m_indexData);