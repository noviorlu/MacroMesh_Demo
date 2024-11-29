#include "HalfEdgeMesh.hpp"

#include <unordered_map>

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

        m_faces.push_back(Face(edge0));

        edge0->face = &m_faces.back();
        edge1->face = &m_faces.back();
        edge2->face = &m_faces.back();

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
            HalfEdge* edge = m_faces[faceIdx].edge;
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
            } while (edge != m_faces[faceIdx].edge);
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
