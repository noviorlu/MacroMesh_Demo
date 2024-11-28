#include "HalfEdgeMesh.hpp"

#include <unordered_map>
#include <utility>
#include <iostream>
#include <algorithm>

// #define HALF_EDGE_MESH_SANITY_CHECK

HalfEdgeMesh::HalfEdgeMesh(const Mesh& mesh) {
    faces.reserve(mesh.m_indexData.size() / 3);
    vertices.reserve(mesh.m_vertexData.size());
    edges.reserve(mesh.m_indexData.size());
    
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

        faces.push_back(Face(edge0));

        edge0->face = &faces.back();
        edge1->face = &faces.back();
        edge2->face = &faces.back();

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

#ifdef HALF_EDGE_MESH_SANITY_CHECK
    for (size_t i = 0; i < edges.size(); ++i) {
        const HalfEdge* edge = edges[i];

        // 检查 edge 的 face 是否为空
        if (!edge->face) {
            std::cerr << "Edge " << i << " has no face!" << std::endl;
        }

        // 检查 twin 的一致性
        if (!edge->twin) {
            std::cerr << "Edge " << i << " has no twin!" << std::endl;
        } else {
            // 检查双向一致性
            if (edge->twin->twin != edge) {
                std::cerr << "Edge " << i 
                        << " twin relationship is inconsistent! "
                        << "twin->twin != edge." << std::endl;
            }
            // 检查 twin 的 face 是否为空
            if (!edge->twin->face) {
                std::cerr << "Edge " << i 
                        << " has a twin with no face!" << std::endl;
            }
            // 检查方向是否匹配
            if (edge->origin != edge->twin->next->origin) {
                std::cerr << "Edge " << i 
                        << " twin direction mismatch! "
                        << "edge->origin != edge->twin->next->origin." << std::endl;
            }
        }
    }

    // 检查每个顶点是否与至少一条边连接
    for (size_t i = 0; i < vertices.size(); ++i) {
        const HalfVertex* vertex = vertices[i];
        if (!vertex->edge) {
            std::cerr << "Vertex " << i << " has no edge!" << std::endl;
        }
    }
#endif

    edges.shrink_to_fit();
}

HalfEdgeMesh::~HalfEdgeMesh() {
    for (auto edge : edges) {
        delete edge;
    }
    for(auto vertex : vertices) {
        delete vertex;
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
        HalfEdge* edge = face.edge;
        do {
            size_t index = vertex_to_index[edge->origin];
            mesh.m_indexData.push_back(static_cast<unsigned int>(index));
            edge = edge->next;
        } while (edge != face.edge);
    }
}

// void HalfEdgeMesh::exportMesh(std::vector<Cluster>& clusterList, std::vector<ClusterGroup>& clusterGroupList) {

//     std::unordered_map<const HalfVertex*, size_t> vertex_to_index;
//     for (size_t i = 0; i < vertices.size(); ++i) {
//         vertex_to_index[vertices[i]] = i;
//         Vertex v(vertices[i]->position, vertices[i]->normal, vertices[i]->uv);
//         mesh.m_vertexData.push_back(v);
//     }

//     for (auto& face : faces) {
//         HalfEdge* edge = face.edge;
//         do {
//             size_t index = vertex_to_index[edge->origin];
//             mesh.m_indexData.push_back(static_cast<unsigned int>(index));
//             edge = edge->next;
//         } while (edge != face.edge);
//     }
// }