
#include "HalfEdgeMesh.hpp"

#include <unordered_set>
#include <iostream>
#include <algorithm>

#include <glm/gtx/string_cast.hpp>

void HalfVertex::cleanAdjacentVerticesQuadric() {
    HalfEdge* edge = this->edge;
    HalfEdge* startEdge = edge;
    do {
        edge->origin->quadric = glm::mat4(0.0f);
        edge = edge->twin->next;
    } while (edge != startEdge);
}

std::vector<HalfVertex*> HalfVertex::getAdjacentVertices(){
    std::vector<HalfVertex*> vertices;
    HalfEdge* edge = this->edge;
    HalfEdge* startEdge = edge;
    do {
        vertices.push_back(edge->twin->origin);
        edge = edge->twin->next;
    } while (edge != startEdge);
    return vertices;
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
    // Assert origin and twin are not null
    assert(origin && "Origin vertex is null.");
    assert(twin && "Twin edge is null.");
    assert(twin->origin && "Twin's origin vertex is null.");

    if (!force && cost != nullptr) {
        return;
    }

    if (force && this > twin) {
        return;
    }

    if (cost != nullptr) {
        delete cost;
        cost = nullptr;
        twin->cost = nullptr;
    }

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

        // Debugging optimal position
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

    // Ensure costValue is non-negative
    costValue = glm::max(costValue, 0.0f);

    // Debugging final cost value
    this->cost = new Cost{costValue, lerpValue, distanceToLine};
    this->twin->cost = this->cost;

    // Assert cost is successfully set
    assert(this->cost != nullptr && "[QEM:ERROR] Cost allocation failed.");
    assert(this->twin->cost == this->cost && "[QEM:ERROR] Twin's cost should match current edge's cost.");
}


void HalfEdge::LerpVertex(Vertex& vert) {
    if (cost == nullptr) {
        assert(cost != nullptr && "[QEM:ERROR] Cost is not computed for this edge");
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

HalfVertex* HalfEdgeMesh::mergeEdge(HalfEdge* edge) {
    // 1. remove v2 and set v1 to the new vertex
    HalfVertex* v1 = edge->origin;
    edge->LerpVertex(*v1);
    assert(edge->twin != nullptr && "[QEM:ERROR] Twin edge is null.");
    HalfVertex* v2 = edge->twin->origin;

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

        HalfEdge* innerEdges_twin = innerEdge->twin;
        HalfEdge* prevEdge_twin = innerEdge->prev->twin;

        innerEdges_twin->twin = prevEdge_twin;
        prevEdge_twin->twin = innerEdges_twin;

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

    return v1;
}

void HalfEdgeMesh::recomputeCost(HalfVertex* vertex) {
    vertex->cleanAdjacentVerticesQuadric();

    HalfEdge* startEdge = vertex->edge;
    HalfEdge* currentEdge = startEdge;

    do {
        Face* face = currentEdge->face;
        if (face) {
            face->updateFaceQuadric();
        }
        currentEdge = currentEdge->twin->next;
    } while (currentEdge != startEdge);

    currentEdge = startEdge;
    do {
        currentEdge->computeEdgeCost(true);
        currentEdge = currentEdge->twin->next;
    } while (currentEdge != startEdge);
}


void HalfEdgeMesh::clusterGroupQEM(QEMQueue& queue, int targetAmount) {
    float totalError = 0.0f;
    int targetEdgeReduct = targetAmount / 2;
    while (m_vertices.size() > targetEdgeReduct) {
        if (queue.empty()) {
            std::cerr << "Priority queue is empty before reaching targetMerge." << std::endl;
            break;
        }

        HalfEdge* edge = queue.top();
        queue.pop();

        if (edge->cost) {
            totalError += edge->cost->val;
        }

        HalfVertex* v1 = mergeEdge(edge);
        recomputeCost(v1);
    }

    std::cout << "Total simplification error: " << totalError << std::endl;
}

void HalfEdgeMesh::QEM(){
    initCostComputation();
    
    // generate the QEMQueue
    QEMQueue queue;
    for (HalfEdge* edge : m_edges) {
        if (edge < edge->twin) {
            queue.push(edge);
        }
    }

    int targetAmount = m_vertices.size() / 2;
    clusterGroupQEM(queue, targetAmount);
}

