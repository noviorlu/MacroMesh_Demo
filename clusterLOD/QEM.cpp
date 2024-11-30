
#include "HalfEdgeMesh.hpp"

#include <unordered_set>
#include <iostream>
#include <algorithm>

void HalfVertex::cleanAdjacentVerticesQuadric() {
    HalfEdge* edge = this->edge;
    HalfEdge* startEdge = edge;
    do {
        edge->origin->quadric = glm::mat4(0.0f);
        edge = edge->twin->next;
    } while (edge != startEdge);
}

/// @brief before using this function, make sure the vertices quadric is initialized to 0
void Face::updateFaceQuadric(){
    HalfEdge* edge = this->edge;
    HalfEdge* startEdge = edge;
    glm::vec3 v0 = startEdge->origin->position;
    glm::vec3 v1 = startEdge->next->origin->position;
    glm::vec3 v2 = startEdge->next->next->origin->position;

    glm::vec3 normal = glm::cross(v1 - v0, v2 - v0);
    float area = glm::length(normal);

    // Ignore faces with zero area
    if (area < 1e-6f) {
        return;
    }

    normal = glm::normalize(normal);

    // Compute the plane equation(ax + by + cz + d = 0)
    float d = -glm::dot(normal, v0);
    glm::vec4 plane(normal, d);

    // Compute the quadric matrix from the plane equation
    glm::mat4 quadric = glm::outerProduct(plane, plane);

    do {
        edge->origin->quadric += quadric;
        edge = edge->next;
    } while (edge != startEdge);
}

void HalfEdge::computeEdgeCost(bool force = false) {
    if (!force && cost != nullptr) return;
    if(force && this < twin) return;
    if(cost != nullptr) { delete cost; cost = nullptr; this->twin->cost = nullptr; }

    HalfVertex* v1 = origin;
    HalfVertex* v2 = twin->origin;

    glm::mat4 Q = v1->quadric + v2->quadric;

    glm::mat4 Q_sub = Q;
    Q_sub[3] = glm::vec4(0, 0, 0, 1);

    float lerpValue;
    float distanceToLine;

    glm::vec3 optimalPosition;
    if (glm::determinant(Q_sub) != 0.0f) {
        glm::vec4 v_optimal = glm::inverse(Q_sub) * glm::vec4(0, 0, 0, 1);
        optimalPosition = glm::vec3(v_optimal);

        glm::vec3 v1ToV2 = v2->position - v1->position;
        float lengthSquared = glm::dot(v1ToV2, v1ToV2);
        if (lengthSquared > 0.0f) {
            lerpValue = glm::dot(optimalPosition - v1->position, v1ToV2) / lengthSquared;
        } else {
            lerpValue = 0.5f;
        }
        lerpValue = glm::clamp(lerpValue, 0.0f, 1.0f);

        glm::vec3 v1ToOptimal = optimalPosition - v1->position;
        distanceToLine = glm::length(glm::cross(v1ToV2, v1ToOptimal)) / glm::length(v1ToV2);
    } else {
        optimalPosition = (v1->position + v2->position) / 2.0f;
        lerpValue = 0.5f;
        distanceToLine = 0.0f;
    }

    glm::vec4 v_opt(optimalPosition, 1.0f);
    float costValue = glm::dot(v_opt, Q * v_opt);
    costValue = glm::max(costValue, 0.0f);

    this->cost = new Cost{costValue, lerpValue, distanceToLine};
    this->twin->cost = this->cost;
}

void HalfEdge::LerpVertex(Vertex& vert) {
    if (cost == nullptr) {
        std::cerr << "[QEM:ERROR] Cost is not computed for this edge" << std::endl;
        exit(1);
    }
    HalfVertex& v1 = *this->origin;
    HalfVertex& v2 = *this->twin->origin;

    glm::vec3 v1ToV2 = v2.position - v1.position;
    glm::vec3 v1ToOptimal = glm::normalize(glm::cross(v1ToV2, glm::vec3(0, 1, 0)));
    if (glm::length(v1ToOptimal) < 1e-6f) {
        v1ToOptimal = glm::normalize(glm::cross(v1ToV2, glm::vec3(1, 0, 0)));
    }

    glm::vec3 lerpedPosition = glm::mix(v1.position, v2.position, cost->lerpValue);
    vert.position = lerpedPosition + cost->distanceToLine * v1ToOptimal;

    glm::vec3 lerpedNormal = glm::normalize(glm::mix(v1.normal, v2.normal, cost->lerpValue));
    vert.normal = glm::normalize(lerpedNormal + cost->distanceToLine * glm::normalize(v1ToOptimal));

    vert.uv = glm::mix(v1.uv, v2.uv, cost->lerpValue);
}

/// @brief Compute the initial quadrics for each vertex in the mesh.
void HalfEdgeMesh::initCostComputation() {
    for (auto& vertex : m_vertices) {
        vertex->quadric = glm::mat4(0.0f);
    }
    
    for (auto& face : m_faces) {
        face->updateFaceQuadric();
    }

    for (auto& edge : m_edges) {
        edge->computeEdgeCost();
    }
}

void HalfEdgeMesh::mergeEdge(HalfEdge* edge) {
    std::cout << "Starting mergeEdge..." << std::endl;

    // 1. remove v2 and set v1 to the new vertex
    HalfVertex* v1 = edge->origin;
    edge->LerpVertex(*v1);
    assert(edge->twin != nullptr);
    HalfVertex* v2 = edge->twin->origin;
    std::cout << "Merging vertices: " << v1 << " and " << v2 << std::endl;

    m_vertices.erase(std::remove(m_vertices.begin(), m_vertices.end(), v2), m_vertices.end());

    auto processFace = [&](Face* face) {
        HalfEdge* innerEdge = nullptr;

        HalfEdge* startEdge = face->edge;
        HalfEdge* currentEdge = startEdge;

        std::vector<HalfEdge*> innerEdges;

        do {
            if (currentEdge->origin != v1 && currentEdge->origin != v2) {
                innerEdge = currentEdge;
            }
            innerEdges.push_back(currentEdge);
            m_edges.erase(std::remove(m_edges.begin(), m_edges.end(), currentEdge), m_edges.end());
            currentEdge = currentEdge->next;
        } while (currentEdge != startEdge);

        // 找到 innerEdge 的 twin 并重新连接
        HalfEdge* innerEdges_twin = innerEdge->twin;
        HalfEdge* prevEdge_twin = innerEdge->prev->twin;

        innerEdges_twin->twin = prevEdge_twin;
        prevEdge_twin->twin = innerEdges_twin;

        // 删除 Face 和关联的边
        m_faces.erase(std::remove(m_faces.begin(), m_faces.end(), face), m_faces.end());
        delete face;

        for (auto edge : innerEdges) {
            delete edge;
        }
    };
    // Process faces and retrieve updated twin edges
    if(edge->face != nullptr) processFace(edge->face);
    if(edge->twin->face != nullptr) processFace(edge->twin->face);

    delete v2;
}





