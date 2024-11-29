#pragma once
#include <vector>
#include <glm/glm.hpp>

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <set>

#include "Mesh.hpp"
#include "cs488-framework/Vertex.hpp"

struct HalfVertex;
struct HalfEdge;
struct Face;

struct HalfVertex : public Vertex {
    HalfEdge* edge;
    glm::mat4 quadric;

    HalfVertex(const glm::vec3& position, const glm::vec3& normal, const glm::vec2& uv) 
        : Vertex(position, normal, uv), edge(nullptr), quadric(0.0f) {}

    HalfVertex(const Vertex& vertex)
        : Vertex(vertex), edge(nullptr), quadric(0.0f) {}
};

struct HalfEdge {
    HalfVertex* origin;
    HalfEdge* twin;
    HalfEdge* next;
    HalfEdge* prev;
    Face* face;
};

struct Face {
    HalfEdge* edge;
    Face(HalfEdge* edge) : edge(edge) {}
};

class HalfEdgeMesh {
public:
    HalfEdgeMesh(const Mesh& mesh);
    ~HalfEdgeMesh();
    void exportMesh(Mesh& mesh);

private:
    std::vector<HalfVertex> m_vertices;
    std::vector<HalfEdge*> m_edges;
    std::vector<Face> m_faces;

    std::vector<size_t> m_clusterOffsets;
    std::vector<size_t> m_clusterGroupOffsets;


    // void computeInitialQuadrics();
    // void computeEdgeCosts();
    // void contractEdge(const HalfEdge* edge);
    // glm::vec3 computeOptimalPosition(const HalfEdge* edge, bool& valid);

    // void simplifyMesh(size_t targetVertexCount);

    // std::priority_queue<std::pair<float, HalfEdge>> edgeQueue;

public:
    void partition_loop();
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
    
    std::unordered_map<int, int> clusterToGroupMap;
};
