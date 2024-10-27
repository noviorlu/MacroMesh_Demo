#ifndef CLUSTER_HPP
#define CLUSTER_HPP

#include <iostream>
#include <vector>
#include <memory>

#include "cs488-framework/MeshConsolidator.hpp"

class Cluster : public std::enable_shared_from_this<Cluster> {
public:
    int id;
    std::vector<std::weak_ptr<Cluster>> parents;
    std::vector<std::shared_ptr<Cluster>> children;

    Cluster(int id) : id(id) {
        std::cout << "Cluster " << id << " created." << std::endl;
    }

    ~Cluster() {
        std::cout << "Cluster " << id << " destroyed." << std::endl;
    }

    void addChild(std::shared_ptr<Cluster> child) {
        children.push_back(child);
        child->parents.push_back(shared_from_this());
    }
};

class ClusterGroup {
public:
    std::vector<std::shared_ptr<Cluster>> clusters;

    void addCluster(std::shared_ptr<Cluster> cluster) {
        clusters.push_back(cluster);
    }

    void printClusters() {
        for (const auto& cluster : clusters) {
            std::cout << "Cluster " << cluster->id << ": ";
            std::cout << "Children [";
            for (const auto& child : cluster->children) {
                std::cout << child->id << " ";
            }
            std::cout << "] Parents [";
            for (const auto& parent : cluster->parents) {
                if (auto p = parent.lock()) {
                    std::cout << p->id << " ";
                }
            }
            std::cout << "]" << std::endl;
        }
    }
};

#endif // CLUSTER_HPP