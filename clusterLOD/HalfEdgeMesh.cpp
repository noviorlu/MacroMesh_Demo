#include "HalfEdgeMesh.hpp"
#include "tiny_obj_loader.h"

#include <unordered_map>
#include <unordered_set>

#include <iostream>

#define DEBUG_HALFEDGEMESH


HalfEdgeMesh::HalfEdgeMesh(const std::string& objFilePath) {
    tinyobj::ObjReader reader;

    // 读取 .obj 文件
    if (!reader.ParseFromFile(objFilePath)) {
        if (!reader.Error().empty()) {
            std::cerr << "TinyObjReader: " << reader.Error() << std::endl;
        }
        throw std::runtime_error("Failed to load .obj file");
    }

    if (!reader.Warning().empty()) {
        std::cout << "TinyObjReader: " << reader.Warning() << std::endl;
    }

    // 获取 .obj 文件的属性和形状
    const tinyobj::attrib_t& attrib = reader.GetAttrib();
    const std::vector<tinyobj::shape_t>& shapes = reader.GetShapes();
    const std::vector<tinyobj::material_t>& materials = reader.GetMaterials();

    // 清空之前的顶点、边和面
    m_vertices.clear();
    m_edges.clear();
    m_faces.clear();

    // 创建顶点
    for (size_t i = 0; i < attrib.vertices.size() / 3; ++i) {
        glm::vec3 position = glm::vec3(
            attrib.vertices[3 * i + 0],
            attrib.vertices[3 * i + 1],
            attrib.vertices[3 * i + 2]);

        glm::vec3 normal = glm::vec3(0.0f); // 如果有法线数据，可以读取
        if (!attrib.normals.empty()) {
            normal = glm::vec3(
                attrib.normals[3 * i + 0],
                attrib.normals[3 * i + 1],
                attrib.normals[3 * i + 2]);
        }

        glm::vec2 uv = glm::vec2(0.0f); // 如果有纹理坐标数据，可以读取
        if (!attrib.texcoords.empty()) {
            uv = glm::vec2(
                attrib.texcoords[2 * i + 0],
                attrib.texcoords[2 * i + 1]);
        }

        m_vertices.push_back(new HalfVertex(position, normal, uv));
    }

    // 创建边和面
    for (const auto& shape : shapes) {
        size_t index_offset = 0;

        for (size_t f = 0; f < shape.mesh.num_face_vertices.size(); ++f) {
            size_t fv = shape.mesh.num_face_vertices[f];
            std::vector<HalfEdge*> faceEdges;
            Face* face = new Face(nullptr);

            for (size_t v = 0; v < fv; ++v) {
                tinyobj::index_t idx = shape.mesh.indices[index_offset + v];
                HalfVertex* vertex = m_vertices[idx.vertex_index];

                HalfEdge* edge = new HalfEdge();
                edge->origin = vertex;
                edge->face = face;

                // 设置顶点的边
                if (vertex->edge == nullptr) {
                    vertex->edge = edge;
                }

                faceEdges.push_back(edge);
            }
            index_offset += fv;

            // 链接环形边表
            for (size_t v = 0; v < fv; ++v) {
                size_t next = (v + 1) % fv;
                faceEdges[v]->next = faceEdges[next];
                faceEdges[next]->prev = faceEdges[v];
            }

            face->edge = faceEdges[0];
            m_faces.push_back(face);
            m_edges.insert(m_edges.end(), faceEdges.begin(), faceEdges.end());
        }
    }
}


HalfEdgeMesh::HalfEdgeMesh(const Mesh& mesh) {
    m_faces.reserve(mesh.m_indexData.size() / 3);
    m_vertices.reserve(mesh.m_vertexData.size());
    m_edges.reserve(mesh.m_indexData.size());

#ifdef DEBUG_HALFEDGEMESH
    // 验证顶点数据的唯一性
    std::unordered_set<glm::vec3> uniqueVertices;
    for (const auto& vertex : mesh.m_vertexData) {
        if (!uniqueVertices.insert(vertex.position).second) {
            std::cerr << "[ERROR] Duplicate vertex found at position: (" 
                      << vertex.position.x << ", " << vertex.position.y << ", " << vertex.position.z << ")" << std::endl;
            assert(false && "Duplicate vertex found in mesh vertex data.");
        }
    }

    // 验证索引数据是否有效
    bool hasInvalidIndex = false;
    for (const auto& index : mesh.m_indexData) {
        if (index >= mesh.m_vertexData.size() || index < 0) {
            std::cerr << "[ERROR] Invalid index found in index data: " << index << std::endl;
            hasInvalidIndex = true;
            
        }
    }
    assert(!hasInvalidIndex && "Invalid vertex index in mesh index data.");

    // 验证每条边是否有双向性
    bool hasInvalidEdge = false;
    std::unordered_map<std::pair<int, int>, int, pair_hash> edgeCount;
    for (size_t i = 0; i < mesh.m_indexData.size(); i += 3) {
        int v0 = mesh.m_indexData[i];
        int v1 = mesh.m_indexData[i + 1];
        int v2 = mesh.m_indexData[i + 2];

        auto addEdge = [&](int start, int end) {
            auto edge = std::make_pair(std::min(start, end), std::max(start, end));
            edgeCount[edge]++;
        };

        addEdge(v0, v1);
        addEdge(v1, v2);
        addEdge(v2, v0);
    }

    for (const auto& [edge, count] : edgeCount) {
        if (count != 2) {
            std::cerr << "[ERROR] Edge (" << edge.first << ", " << edge.second 
                      << ") is shared by " << count << " faces, expected 2." << std::endl;
            hasInvalidEdge = true;
        }
    }
    assert(!hasInvalidEdge && "Invalid edge sharing in mesh index data.");
#endif

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

#ifdef DEBUG_HALFEDGEMESH
    std::cout << "[DEBUG_HALFEDGEMESH] HalfEdgeMesh created with " << m_vertices.size() << " vertices, "
              << m_edges.size() << " edges, and " << m_faces.size() << " faces." << std::endl;

    size_t twinErrors = 0;

    for (size_t i = 0; i < m_edges.size(); ++i) {
        if (!m_edges[i]->twin) {
            std::cerr << "[ERROR] Edge at index " << i << " has no twin edge." << std::endl;
            twinErrors++;
        }
    }

    if (twinErrors > 0) {
        std::cerr << "[ERROR] HalfEdgeMesh validation failed: " << twinErrors << " edges are missing twins." << std::endl;
        assert(false && "HalfEdgeMesh validation failed: Missing twin edges.");
    } else {
        std::cout << "[DEBUG_HALFEDGEMESH] HalfEdgeMesh validation passed: All edges have twins." << std::endl;
    }
#endif
}

HalfEdgeMesh::~HalfEdgeMesh() {
    for (auto edge : m_edges) {
        delete edge;
    }
}

void HalfEdgeMesh::exportMesh(Mesh& mesh) {
    mesh.m_vertexData.clear();
    mesh.m_indexData.clear();
    std::unordered_map<HalfVertex*, unsigned int> vertexIndexMap;

    for (auto vertex : m_vertices) {
        mesh.m_vertexData.emplace_back(
            vertex->position,
            vertex->normal,
            vertex->uv
        );
    }

    for (auto face : m_faces) {
        HalfEdge* edge = face->edge;
        do {
            if (vertexIndexMap.find(edge->origin) == vertexIndexMap.end()) {
                vertexIndexMap[edge->origin] = vertexIndexMap.size();
            }
            mesh.m_indexData.push_back(vertexIndexMap[edge->origin]);
            edge = edge->next;
        } while (edge != face->edge);
    }
}


void HalfEdgeMesh::exportMesh(std::vector<Cluster>& clusters, std::vector<ClusterGroup>& clusterGroups) {
    clusters.clear();
    clusterGroups.clear();

    for (size_t i = 0; i < m_clusterOffsets.size(); ++i) {
        size_t startFace = m_clusterOffsets[i];
        size_t endFace = (i + 1 < m_clusterOffsets.size()) ? m_clusterOffsets[i + 1] : m_faces.size();

        Cluster cluster(0.0f);
        std::unordered_map<HalfVertex*, unsigned int> vertexIndexMap;
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
    std::vector<std::string> errorMessages;

    for (size_t i = 0; i < m_vertices.size(); ++i) {
        HalfVertex* vertex = m_vertices[i];
        if (!vertex) {
            errorMessages.emplace_back("[ERROR] Vertex at index " + std::to_string(i) + " is null.");
            continue;
        }
        if (!vertex->edge) {
            errorMessages.emplace_back("[ERROR] Vertex at index " + std::to_string(i) + " has a null edge pointer.");
        } else if (vertex->edge->origin != vertex) {
            errorMessages.emplace_back("[ERROR] Vertex at index " + std::to_string(i) 
                + " is not the origin of its associated edge.");
        }
    }

    for (size_t i = 0; i < m_edges.size(); ++i) {
        HalfEdge* edge = m_edges[i];
        if (!edge) {
            errorMessages.emplace_back("[ERROR] Edge at index " + std::to_string(i) + " is null.");
            continue;
        }
        if (!edge->origin) {
            errorMessages.emplace_back("[ERROR] Edge at index " + std::to_string(i) + " has a null origin vertex.");
        }
        if (!edge->twin) {
            errorMessages.emplace_back("[ERROR] Edge at index " + std::to_string(i) + " has no twin edge.");
        } else if (edge->twin->twin != edge) {
            errorMessages.emplace_back("[ERROR] Edge at index " + std::to_string(i) 
                + " has a twin whose twin does not point back to this edge.");
        }
        if (!edge->next) {
            errorMessages.emplace_back("[ERROR] Edge at index " + std::to_string(i) + " has no next edge.");
        }
        if (!edge->prev) {
            errorMessages.emplace_back("[ERROR] Edge at index " + std::to_string(i) + " has no previous edge.");
        }
    }

    for (size_t i = 0; i < m_faces.size(); ++i) {
        Face* face = m_faces[i];
        if (!face) {
            errorMessages.emplace_back("[ERROR] Face at index " + std::to_string(i) + " is null.");
            continue;
        }
        if (!face->edge) {
            errorMessages.emplace_back("[ERROR] Face at index " + std::to_string(i) + " has no associated edge.");
        } else {
            HalfEdge* startEdge = face->edge;
            HalfEdge* currentEdge = startEdge;
            size_t edgeCount = 0;

            do {
                if (!currentEdge) {
                    errorMessages.emplace_back("[ERROR] Null edge encountered in face at index " + std::to_string(i) + ".");
                    break;
                }
                if (currentEdge->face != face) {
                    errorMessages.emplace_back("[ERROR] Edge in face at index " + std::to_string(i) 
                        + " does not point back to the correct face.");
                }
                currentEdge = currentEdge->next;
                edgeCount++;
                if (edgeCount > m_edges.size()) {
                    errorMessages.emplace_back("[ERROR] Infinite loop detected in face at index " + std::to_string(i) + ".");
                    break;
                }
            } while (currentEdge != startEdge);
        }
    }

    // If there are any errors, print them and assert
    if (!errorMessages.empty()) {
        for (const auto& msg : errorMessages) {
            std::cerr << msg << std::endl;
        }
        assert(false && "HalfEdgeMesh validation failed. Check error messages.");
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
