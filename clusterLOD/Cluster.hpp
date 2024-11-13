#ifndef CLUSTER_HPP
#define CLUSTER_HPP

#include <iostream>
#include <vector>
#include <memory>

#include "Mesh.hpp"

class Cluster : public std::enable_shared_from_this<Cluster> {
public:
    float Error;
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

#endif // CLUSTER_HPP
