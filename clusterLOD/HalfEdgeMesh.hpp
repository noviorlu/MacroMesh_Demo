#pragma once
#include "Mesh.hpp"
#include "cs488-framework/Vertex.hpp"

#include <glm/glm.hpp>

#include <queue>

struct HalfVertex;
struct HalfEdge;
struct Edge;
struct Face;

struct HalfVertex : public Vertex {
    HalfEdge* edge;
    glm::mat4 quadric;

    HalfVertex(const glm::vec3& position, const glm::vec3& normal, const glm::vec2& uv) 
        : Vertex(position, normal, uv), edge(nullptr), quadric(0.0f) {}

    HalfVertex(const Vertex& vertex)
        : Vertex(vertex), edge(nullptr), quadric(0.0f) {}
};

struct Face {
    HalfEdge* edge;
    Face(HalfEdge* edge) : edge(edge) {}
};

struct HalfEdge {
    HalfVertex* origin;
    HalfEdge* twin;
    HalfEdge* next;
    HalfEdge* prev;
    Face* face;
};

struct EdgeCost {
    float cost;
    glm::vec3 optimalPosition;
    HalfEdge* edge;

    bool operator<(const EdgeCost& other) const {
        return cost < other.cost;
    }
};
typedef std::priority_queue<EdgeCost, std::vector<EdgeCost>, std::greater<EdgeCost>> QEMqueue;

class HalfEdgeMesh {
public:
    HalfEdgeMesh(const Mesh& mesh);
    ~HalfEdgeMesh();
    void exportMesh(Mesh& mesh);
    void exportMesh(std::vector<Cluster>& clusters, std::vector<ClusterGroup>& clusterGroups);

private:
    std::vector<HalfVertex> m_vertices;
    std::vector<HalfEdge*> m_edges;
    std::vector<Face> m_faces;

public:
    void partition_loop();

// private:
//     void computeInitialQuadrics();
//     void computeEdgeCost(HalfEdge* edge, EdgeCost& edgeCost);
//     void simplifyMesh();
//     void mergeVertices(HalfVertex* v1, HalfVertex* v2, const glm::vec3& newPosition);

private:
    void BuildAdjacencyListForRange(
        std::vector<std::vector<size_t>>& adjacency_list,
        size_t startIdx,
        size_t endIdx
    );
    void HalfEdgeMeshSplitter();
    void HalfEdgeMeshSplitterRecursive(
        size_t startIdx,
        size_t endIdx,
        bool isParentClusterGroup
    );

    std::vector<size_t> m_clusterOffsets;
    std::vector<size_t> m_clusterGroupOffsets;
};
