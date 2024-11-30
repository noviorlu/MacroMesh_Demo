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

    HalfVertex() : Vertex(), edge(nullptr), quadric(0.0f) {}

    HalfVertex(const glm::vec3& position, const glm::vec3& normal, const glm::vec2& uv) 
        : Vertex(position, normal, uv), edge(nullptr), quadric(0.0f) {}

    HalfVertex(const Vertex& vertex)
        : Vertex(vertex), edge(nullptr), quadric(0.0f) {}

    void cleanAdjacentVerticesQuadric();

    std::vector<HalfVertex*> getAdjacentVertices();
};

struct Face {
    HalfEdge* edge;
    Face(HalfEdge* edge) : edge(edge) {}

    void updateFaceQuadric();
};

struct HalfEdge {
    HalfVertex* origin;
    HalfEdge* twin;
    HalfEdge* next;
    HalfEdge* prev;
    Face* face;

    struct Cost{
        float val;
        float lerpValue; 
        float distanceToLine;
        
        bool operator<(const Cost& other) const {
            return val < other.val;
        }
    };
    Cost* cost;

    HalfEdge() {
        origin = nullptr;
        twin = nullptr;
        next = nullptr;
        prev = nullptr;
        face = nullptr;
        cost = nullptr;
    }
    ~HalfEdge() {
        if(cost != nullptr) {
            delete cost;
            cost = nullptr;
            twin->cost = nullptr;
        }
        
    }
    void LerpVertex(Vertex& vert);
    void computeEdgeCost(bool force);
};

struct CompareHalfEdgeCost {
    bool operator()(HalfEdge* a, HalfEdge* b) const {
        if (a->cost == nullptr || b->cost == nullptr) {
            return a->cost == nullptr;
        }
        return a->cost->val < b->cost->val;
    }
};

using QEMQueue = std::priority_queue<HalfEdge*, std::vector<HalfEdge*>, CompareHalfEdgeCost>;


class HalfEdgeMesh {
public:
    HalfEdgeMesh(){}
    HalfEdgeMesh(const Mesh& mesh);
    HalfEdgeMesh(const std::string& objFilePath);
    ~HalfEdgeMesh();
    void exportMesh(Mesh& mesh);
    void exportMesh(std::vector<Cluster>& clusters, std::vector<ClusterGroup>& clusterGroups);

private:
    std::vector<HalfVertex*> m_vertices;
    std::vector<HalfEdge*> m_edges;
    std::vector<Face*> m_faces;

public:
    void partition_loop();
    void QEM();
private:
    void initCostComputation();
    HalfVertex* mergeEdge(HalfEdge* edge);
    void recomputeCost(HalfVertex* vertex);
    void clusterGroupQEM(QEMQueue& queue, int targetMerge);
    

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



/// Debugging functions
    void HalfEdgeMeshValidation();
    void HalfEdgeMeshPrint();
    friend void QEMTestcase();
};
