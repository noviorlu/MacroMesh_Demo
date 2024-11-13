#pragma once

#include <iostream>
#include <vector>
#include <memory>

#include "Mesh.hpp"

class Cluster : public std::enable_shared_from_this<Cluster> {
public:
    float Error;
    
    Mesh m_mesh;

    Cluster(float Error) : Error(Error) {}

    ~Cluster() {}
};

class ClusterGroup {
public:
    float Error;
    std::vector<std::shared_ptr<Cluster>> clusters;

    void addCluster(std::shared_ptr<Cluster> cluster) {
        clusters.push_back(cluster);
    }
};
