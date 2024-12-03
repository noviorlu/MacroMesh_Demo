#pragma once
#include "Mesh.hpp"
#include "cs488-framework/Vertex.hpp"

#include <glm/glm.hpp>

#include <queue>
#include <unordered_set>
#include <vector>
#include <array>
#include <algorithm>

#include <chrono>


struct HalfVertex;
struct HalfEdge;
struct Face;
struct Cost {
    float val;
    glm::vec3 optimalPosition;

    bool operator<(const Cost& other) const {
        return val < other.val;
    }
};

struct HalfVertex : public Vertex {
    bool isValid = true;
    
    HalfEdge* edge = nullptr;
    glm::mat4 quadric;

    HalfVertex() : Vertex(), edge(nullptr), quadric(0.0f) {}

    HalfVertex(const glm::vec3& position, const glm::vec3& normal, const glm::vec2& uv) 
        : Vertex(position, normal, uv), edge(nullptr), quadric(0.0f) {}

    HalfVertex(const Vertex& vertex)
        : Vertex(vertex), edge(nullptr), quadric(0.0f) {}

    void cleanAdjacentVerticesQuadric();

    std::vector<HalfVertex*> getAdjacentVertices();
};
struct HalfEdge {
	bool isValid = true;

    HalfVertex* origin = nullptr;
    HalfEdge* twin = nullptr;
    HalfEdge* next = nullptr;
    HalfEdge* prev = nullptr;
    Face* face = nullptr;
    Cost cost;

    HalfEdge() {
        origin = nullptr;
        twin = nullptr;
        next = nullptr;
        prev = nullptr;
        face = nullptr;
    }
    ~HalfEdge() {}
    void LerpVertex(Vertex& vert);
    void computeEdgeCost(bool force);
};

struct HashHalfEdge {
    size_t operator()(const HalfEdge* edge) const {
        if (edge == nullptr || edge->origin == nullptr || edge->twin == nullptr || edge->twin->origin == nullptr) {
            return 0;
        }
        HalfVertex* v1 = edge->origin;
        HalfVertex* v2 = edge->twin->origin;

        size_t h1 = std::hash<HalfVertex*>()(v1);
        size_t h2 = std::hash<HalfVertex*>()(v2);
        return h1 ^ (h2 * 31);
    }
};
struct EqualHalfEdge {
    bool operator()(const HalfEdge* a, const HalfEdge* b) const {
        if (a == nullptr || b == nullptr) return false;
        HalfVertex* a1 = a->origin;
        HalfVertex* a2 = a->twin->origin;
        HalfVertex* b1 = b->origin;
        HalfVertex* b2 = b->twin->origin;

        // Check if the sets {a1, a2} and {b1, b2} are equal (ignoring order)
        return (a1 == b1 && a2 == b2);
    }
};
struct CompareHalfEdgeCost {
    bool operator()(HalfEdge* a, HalfEdge* b) const {
        return a->cost.val < b->cost.val;
    }
};

using QEMQueue = std::priority_queue<HalfEdge*, std::vector<HalfEdge*>, CompareHalfEdgeCost>;

struct Face {
    bool isValid = true;

    HalfEdge* edge = nullptr;
    Face(HalfEdge* edge) : edge(edge) {}

    void updateFaceQuadric();

    friend std::ostream& operator<<(std::ostream& os, const Face& face) {
        os << "Face Info:\n";
        os << "Face pointer: " << &face << "\n";
        os << " - isValid: " << (face.isValid ? "true" : "false") << "\n";
        os << " - edge: " << face.edge;
        os << " - next: " << face.edge->next;
        os << " - prev: " << face.edge->prev << "\n";
        return os;
    }
};

class HalfEdgeMesh {
public:
    std::string m_name;

    HalfEdgeMesh(){}
    HalfEdgeMesh(const Mesh& mesh);
    ~HalfEdgeMesh();
    void importMesh(const std::string& objFilePath);
    void exportMesh(Mesh& mesh);
    void exportMeshToObjFiles(const std::string& folderPath);
    void exportMesh(std::vector<Cluster>& clusters, std::vector<ClusterGroup>& clusterGroups);
    void exportMesh(const std::string& objFilePath);

    std::vector<HalfVertex*> m_vertices;
    std::vector<HalfEdge*> m_edges;
    std::vector<Face*> m_faces;

public:
    void QEM();

    void initCostComputation();
    void edgeCollapse(HalfEdge* edge);
    void recomputeCost(HalfVertex* vertex);
    void clusterGroupQEM(QEMQueue& queue, int targetMerge);
    std::unordered_set<HalfEdge*> m_garbageEdgeCollector;

public:
    void BuildAdjacencyListForCluster(
        std::vector<std::vector<size_t>>& adjacency_list
    );
    void BuildAdjacencyListForRange(
        std::vector<std::vector<size_t>>& adjacency_list,
        size_t startIdx,
        size_t endIdx
    );
    void HalfEdgeMeshSplitter();
    void HalfEdgeMeshSplitterRecursive(
        size_t startIdx,
        size_t endIdx
    );

    std::vector<size_t> m_clusterOffsets;
    std::vector<size_t> m_clusterGroupResult; // same size as m_clusterOffsets
    int m_clusterGroupCount = 0;
/// Debugging functions
    void HalfEdgeMeshValidation();
    void HalfEdgeMeshPrint();
    friend void QEMTestcase();
    friend void QEMTestcase1();
};

using TripleHalfVertex = std::array<HalfVertex*, 3>;
struct TripleHalfVertexHash {
    size_t operator()(const TripleHalfVertex& vertices) const {
        TripleHalfVertex sortedVertices = vertices;
        std::sort(sortedVertices.begin(), sortedVertices.end());

        size_t seed = 0;
        for (HalfVertex* v : sortedVertices) {
            seed ^= std::hash<HalfVertex*>{}(v)+0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

struct TripleHalfVertexEqual {
    bool operator()(const TripleHalfVertex& lhs, const TripleHalfVertex& rhs) const {
        TripleHalfVertex sortedLhs = lhs;
        TripleHalfVertex sortedRhs = rhs;
        std::sort(sortedLhs.begin(), sortedLhs.end());
        std::sort(sortedRhs.begin(), sortedRhs.end());
        return sortedLhs == sortedRhs;
    }
};



struct EVertex;
struct EEdge;
struct EFace;
struct VertexPair {
    EVertex* v1;
    EVertex* v2;
    VertexPair(const VertexPair& other) : v1(other.v1), v2(other.v2) {}
    VertexPair(EVertex* vertex1, EVertex* vertex2) {
		assert(vertex1 != nullptr && vertex2 != nullptr);
		assert(vertex1 != vertex2);
        if (vertex1 < vertex2) {
            v1 = vertex1;
            v2 = vertex2;
        }
        else {
            v1 = vertex2;
            v2 = vertex1;
        }
    }
    bool operator==(const VertexPair& other) const {
        return v1 == other.v1 && v2 == other.v2;
    }

    // assign operator
	VertexPair& operator=(const VertexPair& other) {
		v1 = other.v1;
		v2 = other.v2;
		return *this;
	}
};

struct EVertex : public Vertex {
    bool isFakeBoundary = false;
    glm::mat4 quadric;
    std::unordered_set<EEdge*> edges;
    void addEdge(EEdge* edge) {edges.insert(edge);}
    void removeEdge(EEdge* edge) {edges.erase(edge);}
    EEdge* popEdge() {
        if (edges.empty()) return nullptr;
        EEdge* edge = *edges.begin();
        edges.erase(edges.begin());
        return edge;
    }
	EVertex(glm::vec3 position, glm::vec3 normal, glm::vec2 uv) : Vertex(position, normal, uv) {
        edges.clear();
    }
};

struct EEdge {
    bool isDirty = false; // for piority Queue Lazy deletion & reinsert
	bool isValid = true;
    bool isBoundary = false;
    bool isFakeBoundary = false;
	VertexPair vertices;
	std::unordered_set<EFace*> faces; // should be size of 2, might have cases where Face back to back thus Face size of 4
    
    Cost cost;
    void addFace(EFace* face) { faces.insert(face); }
    bool containFace(EFace* face) { return faces.find(face) != faces.end(); }
    void removeFace(EFace* face) { faces.erase(face); }
    EFace* popFace() {
        if (faces.empty()) return nullptr;
        EFace* face = *faces.begin();
        faces.erase(faces.begin());
        return face;
    }
    EEdge(const VertexPair& pair, EFace* face) : vertices(pair) {
        faces.clear();
        addFace(face);
		cost.val = 0.0f;
		cost.optimalPosition = glm::vec3(0.0f);
    }

    void computeEdgeCost();
};

struct EFace {
    bool isValid = true;
    std::array<EVertex*, 3> vertices; // should be size of 3
	EFace(EVertex* v0, EVertex* v1, EVertex* v2) {
		vertices[0] = v0;
		vertices[1] = v1;
		vertices[2] = v2;

        adjacentFaces.reserve(3);
	}
    void updateFaceQuadric();

    std::vector<EFace*> adjacentFaces; // at most 3, for partition propose
};

class EEdgePriorityQueue {
    struct EEdgeCompare {
        bool operator()(const EEdge* lhs, const EEdge* rhs) const {
            return lhs->cost.val > rhs->cost.val;
        }
    };
    std::priority_queue<EEdge*, std::vector<EEdge*>, EEdgeCompare> pq;
public:
    void push(EEdge* edge) { pq.push(edge); }
    EEdge* pop() {
        while (!pq.empty()) {
            EEdge* topEdge = pq.top();
            pq.pop();
            if (!topEdge->isValid) continue;
            if (topEdge->isDirty) {
                topEdge->isDirty = false;
                push(topEdge);
                continue;
            }
            return topEdge;
        }
        return nullptr;
    }
    bool empty() const { return pq.empty(); }
};


class EMesh {
public:
	std::vector<EVertex*> m_vertices;
	std::vector<EFace*> m_faces;
    struct VertexPairHash {
        std::size_t operator()(const VertexPair& pair) const {
            std::size_t h1 = std::hash<const EVertex*>()(pair.v1);
            std::size_t h2 = std::hash<const EVertex*>()(pair.v2);
            return h1 ^ (h2 << 1);
        }
    };
	std::unordered_map<VertexPair, EEdge*, VertexPairHash> m_edgeMap;

public:
    std::string m_name;
	void importEMesh(const std::string& objFilePath);
	void exportEMesh(const std::string& objFilePath);
    void exportEMeshToObjFiles(const std::string& folderPath);
    void exportEMesh(std::vector<Cluster>& clusters, std::vector<ClusterGroup>& clusterGroups);
public:
    void buildAdjacencyFaces();
    void buildAdjacencyListForRange(
        std::vector<std::vector<size_t>>& adjacency_list,
        size_t startIdx,
        size_t endIdx
    );
    void eMeshSplitterRecursive(
        size_t startIdx,
        size_t endIdx
    );
    void buildAdjacencyListForCluster(
        std::vector<std::vector<size_t>>& adjacency_list
    );

    void eMeshSplitter();
    std::vector<size_t> m_clusterOffsets;
    std::vector<size_t> m_clusterGroupResult; // same size as m_clusterOffsets
    int m_clusterGroupCount = 0;

public:
    void QEM(float ratio);
    void edgeCollapse(EEdge* edge);
    int m_validFaces = 0;

public:
    void validate();
    void print();
};


void ClusterGrouper(
    const std::vector<std::vector<size_t>>& adjacencyList,
    std::vector<size_t>& clusterGroupResult,
    int& clusterGroupCount
);
