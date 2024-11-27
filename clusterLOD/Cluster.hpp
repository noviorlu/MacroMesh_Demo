// #pragma once

// #include <iostream>
// #include <vector>
// #include <memory>
// #include "cs488-framework/GlErrorCheck.hpp"
// #include "Mesh.hpp"
// #include <metis.h>

// #include <unordered_set>
// #include <queue>
// #include <string>
// #include <glm/glm.hpp>
// #include <glm/gtx/matrix_operation.hpp>

// #define MAX_TRI_IN_CLUSTER 256
// #define MAXN_CLUSTER_IN_CLUSTERGROUP 32


// struct Edge {
//     unsigned int v1, v2;

//     Edge(unsigned int a, unsigned int b) {
//         v1 = std::min(a, b);
//         v2 = std::max(a, b);
//     }

//     bool operator==(const Edge& other) const {
//         return v1 == other.v1 && v2 == other.v2;
//     }
// };

// struct EdgeHash {
//     size_t operator()(const Edge& edge) const {
//         return std::hash<unsigned int>()(edge.v1) ^ std::hash<unsigned int>()(edge.v2);
//     }
// };

// class Cluster : public Mesh {
// public:
//     float Error;
    
//     glm::vec3 rdColor;

//     std::vector<Cluster*> adjacent_clusters;

//     Cluster(
//         float Error, 
//         const Mesh& ref, 
//         const std::vector<unsigned int>& triIndices
//     );

//     ~Cluster() {}

// 	void draw(const ShaderProgram& shader) const override;
// };

// class ClusterGroup {
// public:
//     float Error;
    
//     std::vector<Cluster*> clusters;

//     void GenMesh(Mesh& mesh);
// };

// void MeshSplitter(Mesh& mesh);

// void clusterGrouping(const Mesh& mesh_ref, std::vector<ClusterGroup*>& cluster_groups);

// std::vector<std::vector<size_t>> 
// BuildAdjacencyList(const std::vector<unsigned int>& m_indexData);