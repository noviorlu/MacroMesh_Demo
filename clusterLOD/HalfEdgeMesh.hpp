#pragma once
#include <vector>
#include <glm/glm.hpp>

#include <vector>
#include <unordered_map>
#include <unordered_set>

#include "Mesh.hpp"
#include "cs488-framework/Vertex.hpp"

struct HalfVertex;
struct HalfEdge;
struct Face;

struct HalfVertex : public Vertex {
    HalfEdge* edge;
    glm::mat4 quadric;

    HalfVertex(const glm::vec3& position, const glm::vec3& normal, const glm::vec2& uv) 
        : Vertex(position, normal, uv) 
    {
        edge = nullptr;
        quadric = glm::mat4(0.0f);
    }
};

struct HalfEdge {
    HalfVertex* origin;
    HalfEdge* twin;
    HalfEdge* next;
    HalfEdge* prev;
    Face* face;
};

struct Face {
    int index;
    int clusterIndex;
    HalfEdge* edge;

    Face(int index, HalfEdge* edge) : index(index), edge(edge) {
        clusterIndex = -1;
    }
};

class HalfEdgeMesh {
public:
    HalfEdgeMesh(const Mesh& mesh);
    ~HalfEdgeMesh();
    void exportMesh(Mesh& mesh);
private:
    std::vector<HalfVertex*> vertices;
    std::vector<HalfEdge*> edges;
    std::vector<Face*> faces;


    // void computeInitialQuadrics();
    // void computeEdgeCosts();
    // void contractEdge(const HalfEdge* edge);
    // glm::vec3 computeOptimalPosition(const HalfEdge* edge, bool& valid);

    // void simplifyMesh(size_t targetVertexCount);



    
    // std::priority_queue<std::pair<float, HalfEdge>> edgeQueue;



public:
    void partition_loop();
private:
    void BuildAdjacencyListFromHalfEdgeMesh(std::vector<std::vector<size_t>>& adjacency_list);
    void HalfEdgeMeshSplitter();
    void BuildClusterAdjacency(std::unordered_map<int, std::unordered_set<int>>& cluster_adjacency);
    void GroupClustersWithMETIS();
    
    std::unordered_map<int, int> clusterToGroupMap;
};
