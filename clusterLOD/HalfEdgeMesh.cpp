#include "HalfEdgeMesh.hpp"

#include <unordered_map>
#include <utility>

struct pair_hash {
    template <class T1, class T2>
    std::size_t operator()(const std::pair<T1, T2>& pair) const {
        auto hash1 = std::hash<T1>{}(pair.first);
        auto hash2 = std::hash<T2>{}(pair.second);
        return hash1 ^ (hash2 << 1);
    }
};

HalfEdgeMesh::HalfEdgeMesh(const Mesh& mesh) {
    for (size_t i = 0; i < mesh.m_vertexData.size(); ++i) {
        const Vertex& vertexData = mesh.m_vertexData[i];
        HalfVertex* vertex = new HalfVertex(vertexData.position, vertexData.normal, vertexData.uv);
        vertices.push_back(vertex);
    }

    std::unordered_map<std::pair<int, int>, HalfEdge*, pair_hash> edgeMap;

    for (size_t i = 0; i < mesh.m_indexData.size(); i += 3) {
        int v0 = mesh.m_indexData[i];
        int v1 = mesh.m_indexData[i + 1];
        int v2 = mesh.m_indexData[i + 2];

        HalfEdge* edge0 = new HalfEdge();
        HalfEdge* edge1 = new HalfEdge();
        HalfEdge* edge2 = new HalfEdge();

        edge0->origin = vertices[v0];
        edge1->origin = vertices[v1];
        edge2->origin = vertices[v2];

        edge0->next = edge1;
        edge1->next = edge2;
        edge2->next = edge0;

        edge0->prev = edge2;
        edge1->prev = edge0;
        edge2->prev = edge1;

        Face *face = new Face(faces.size(), edge0);
        faces.push_back(face);

        edge0->face = faces.back();
        edge1->face = faces.back();
        edge2->face = faces.back();

        edges.push_back(edge0);
        edges.push_back(edge1);
        edges.push_back(edge2);

        auto addEdge = [&](int start, int end, HalfEdge* edge) {
            auto twinIt = edgeMap.find({end, start});
            if (twinIt != edgeMap.end()) {
                edge->twin = twinIt->second;
                twinIt->second->twin = edge;
            } else {
                edge->twin = nullptr;
                edgeMap[{start, end}] = edge;
            }
        };

        addEdge(v0, v1, edge0);
        addEdge(v1, v2, edge1);
        addEdge(v2, v0, edge2);
    }

    for (auto& edge : edges) {
        if (!edge->origin->edge) {
            edge->origin->edge = edge;
        }
    }
}

HalfEdgeMesh::~HalfEdgeMesh() {
    for (auto edge : edges) {
        delete edge;
    }
    for(auto vertex : vertices) {
        delete vertex;
    }
    for(auto face : faces) {
        delete face;
    }
}

void HalfEdgeMesh::exportMesh(Mesh& mesh) {
    mesh.m_vertexData.clear();
    mesh.m_indexData.clear();

    std::unordered_map<const HalfVertex*, size_t> vertex_to_index;
    for (size_t i = 0; i < vertices.size(); ++i) {
        vertex_to_index[vertices[i]] = i;
        Vertex v(vertices[i]->position, vertices[i]->normal, vertices[i]->uv);
        mesh.m_vertexData.push_back(v);
    }

    for (auto& face : faces) {
        if(clusterToGroupMap[face->clusterIndex] != 0) continue;
        HalfEdge* edge = face->edge;
        do {
            size_t index = vertex_to_index[edge->origin];
            mesh.m_indexData.push_back(static_cast<unsigned int>(index));
            edge = edge->next;
        } while (edge != face->edge);
    }
}

// /// @brief Compute the initial quadrics for each vertex in the mesh.
// void HalfEdgeMesh::computeInitialQuadrics() {
//     // eterate all faces
//     for (auto& face : faces) {
//         HalfEdge* edge = face.edge;

//         // compute the normal of the face
//         glm::vec3 v0 = edge->origin->position;
//         glm::vec3 v1 = edge->next->origin->position;
//         glm::vec3 v2 = edge->next->next->origin->position;

//         glm::vec3 normal = glm::cross(v1 - v0, v2 - v0);
//         float area = glm::length(normal);

//         // ignore faces with zero area
//         if (area < 1e-6f) {
//             continue;
//         }

//         normal = glm::normalize(normal);

//         // compute the plane equation(ax + by + cz + d = 0)
//         float d = -glm::dot(normal, v0);
//         glm::vec4 plane(normal, d);

//         // compute the quadric matrix from the plane equation
//         // [a^2, ab, ac, ad]
//         // [ab, b^2, bc, bd]
//         // [ac, bc, c^2, cd]
//         // [ad, bd, cd, d^2]
//         glm::mat4 quadric = glm::outerProduct(plane, plane);

//         // add the quadric matrix to all vertices of the face
//         HalfEdge* currentEdge = edge;
//         do {
//             currentEdge->origin->quadric += quadric;
//             currentEdge = currentEdge->next;
//         } while (currentEdge != edge);
//     }
// }
