#include "HalfEdgeMesh.hpp"

#include <unordered_map>
#include <iostream>

HalfEdgeMesh::HalfEdgeMesh(const Mesh& mesh) {
    m_faces.reserve(mesh.m_indexData.size() / 3);
    m_vertices.reserve(mesh.m_vertexData.size());
    m_edges.reserve(mesh.m_indexData.size());
    
    for (size_t i = 0; i < mesh.m_vertexData.size(); ++i) {
        m_vertices.push_back(new HalfVertex(mesh.m_vertexData[i]));
    }

    std::unordered_map<std::pair<int, int>, HalfEdge*, pair_hash> edgeMap;

    for (size_t i = 0; i < mesh.m_indexData.size(); i += 3) {
        int v0 = mesh.m_indexData[i];
        int v1 = mesh.m_indexData[i + 1];
        int v2 = mesh.m_indexData[i + 2];

        HalfEdge* edge0 = new HalfEdge();
        HalfEdge* edge1 = new HalfEdge();
        HalfEdge* edge2 = new HalfEdge();

        edge0->origin = m_vertices[v0];
        edge1->origin = m_vertices[v1];
        edge2->origin = m_vertices[v2];

        edge0->next = edge1;
        edge1->next = edge2;
        edge2->next = edge0;

        edge0->prev = edge2;
        edge1->prev = edge0;
        edge2->prev = edge1;

        Face* face = new Face(edge0);
        m_faces.push_back(face);

        edge0->face = face;
        edge1->face = face;
        edge2->face = face;

        m_edges.push_back(edge0);
        m_edges.push_back(edge1);
        m_edges.push_back(edge2);

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

    for (auto& edge : m_edges) {
        if (!edge->origin->edge) {
            edge->origin->edge = edge;
        }
    }

    m_edges.shrink_to_fit();
}

HalfEdgeMesh::~HalfEdgeMesh() {
    for (auto edge : m_edges) {
        delete edge;
    }
}

void HalfEdgeMesh::exportMesh(std::vector<Cluster>& clusters, std::vector<ClusterGroup>& clusterGroups) {
    clusters.clear();
    clusterGroups.clear();

    for (size_t i = 0; i < m_clusterOffsets.size(); ++i) {
        size_t startFace = m_clusterOffsets[i];
        size_t endFace = (i + 1 < m_clusterOffsets.size()) ? m_clusterOffsets[i + 1] : m_faces.size();

        Cluster cluster(0.0f);
        std::unordered_map<HalfVertex*, unsigned int> vertexIndexMap; // Key changed to Vertex*.
        unsigned int clusterVertexIndex = 0;

        for (size_t faceIdx = startFace; faceIdx < endFace; ++faceIdx) {
            HalfEdge* edge = m_faces[faceIdx]->edge;
            do {
                HalfVertex* vertex = edge->origin;
                if (vertexIndexMap.find(vertex) == vertexIndexMap.end()) {
                    vertexIndexMap[vertex] = clusterVertexIndex++;
                    cluster.m_vertexData.emplace_back(
                        vertex->position,
                        vertex->normal,
                        vertex->uv
                    );
                }
                cluster.m_indexData.push_back(vertexIndexMap[vertex]);
                edge = edge->next;
            } while (edge != m_faces[faceIdx]->edge);
        }
        clusters.push_back(std::move(cluster));
    }

    for (size_t i = 0; i < m_clusterGroupOffsets.size(); ++i) {
        size_t startCluster = m_clusterGroupOffsets[i];
        size_t endCluster = (i + 1 < m_clusterGroupOffsets.size()) ? m_clusterGroupOffsets[i + 1] : clusters.size();

        ClusterGroup group;
        for (size_t clusterIdx = startCluster; clusterIdx < endCluster; ++clusterIdx) {
            group.clusters.push_back(&clusters[clusterIdx]);
        }

        clusterGroups.push_back(std::move(group));
    }
}

void HalfEdgeMesh::HalfEdgeMeshValidation() {
    for (size_t i = 0; i < m_vertices.size(); ++i) {
        HalfVertex* vertex = m_vertices[i];
        if (!vertex) {
            std::cerr << "[ERROR] Vertex at index " << i << " is null." << std::endl;
            continue;
        }
        if (!vertex->edge) {
            std::cerr << "[ERROR] Vertex at index " << i << " has a null edge pointer." << std::endl;
        } else if (vertex->edge->origin != vertex) {
            std::cerr << "[ERROR] Vertex at index " << i 
                      << " is not the origin of its associated edge." << std::endl;
        }
    }

    for (size_t i = 0; i < m_edges.size(); ++i) {
        HalfEdge* edge = m_edges[i];
        if (!edge) {
            std::cerr << "[ERROR] Edge at index " << i << " is null." << std::endl;
            continue;
        }
        if (!edge->origin) {
            std::cerr << "[ERROR] Edge at index " << i << " has a null origin vertex." << std::endl;
        }
        if (!edge->twin) {
            std::cerr << "[ERROR] Edge at index " << i << " has no twin edge." << std::endl;
        } else if (edge->twin->twin != edge) {
            std::cerr << "[ERROR] Edge at index " << i 
                      << " has a twin whose twin does not point back to this edge." << std::endl;
        }
        if (!edge->next) {
            std::cerr << "[ERROR] Edge at index " << i << " has no next edge." << std::endl;
        }
        if (!edge->prev) {
            std::cerr << "[ERROR] Edge at index " << i << " has no previous edge." << std::endl;
        }
    }

    for (size_t i = 0; i < m_faces.size(); ++i) {
        Face* face = m_faces[i];
        if (!face) {
            std::cerr << "[ERROR] Face at index " << i << " is null." << std::endl;
            continue;
        }
        if (!face->edge) {
            std::cerr << "[ERROR] Face at index " << i << " has no associated edge." << std::endl;
        } else {
            HalfEdge* startEdge = face->edge;
            HalfEdge* currentEdge = startEdge;
            size_t edgeCount = 0;

            do {
                if (!currentEdge) {
                    std::cerr << "[ERROR] Null edge encountered in face at index " << i << "." << std::endl;
                    break;
                }
                if (currentEdge->face != face) {
                    std::cerr << "[ERROR] Edge in face at index " << i 
                              << " does not point back to the correct face." << std::endl;
                }
                currentEdge = currentEdge->next;
                edgeCount++;
                if (edgeCount > m_edges.size()) {
                    std::cerr << "[ERROR] Infinite loop detected in face at index " << i << "." << std::endl;
                    break;
                }
            } while (currentEdge != startEdge);
        }
    }
}

void HalfEdgeMesh::HalfEdgeMeshPrint() {
    std::cout << "Vertices:" << std::endl;
    for (size_t i = 0; i < m_vertices.size(); ++i) {
        HalfVertex* vertex = m_vertices[i];
        std::cout << "Vertex " << i << ": position = (" 
                  << vertex->position.x << ", " << vertex->position.y << ", " << vertex->position.z 
                  << "), pointer = " << vertex << std::endl;
    }

    std::cout << "Edges:" << std::endl;
    for (size_t i = 0; i < m_edges.size(); ++i) {
        HalfEdge* edge = m_edges[i];
        std::cout << "Edge " << i << ": origin position = (" 
                  << edge->origin->position.x << ", " << edge->origin->position.y << ", " << edge->origin->position.z 
                  << "), pointer = " << edge << ", twin = " << edge->twin 
                  << ", next = " << edge->next << ", prev = " << edge->prev << std::endl;
    }

    std::cout << "Faces:" << std::endl;
    for (size_t i = 0; i < m_faces.size(); ++i) {
        Face* face = m_faces[i];
        std::cout << "Face " << i << ": edge pointer = " << face->edge << std::endl;
    }
}
