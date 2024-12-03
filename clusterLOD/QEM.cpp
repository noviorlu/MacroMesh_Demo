﻿
#include "HalfEdgeMesh.hpp"

#include <unordered_set>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <thread>

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
void Face::updateFaceQuadric() {
    if (!edge) return;

    glm::vec3 v0 = edge->origin->position;
    glm::vec3 v1 = edge->next->origin->position;
    glm::vec3 v2 = edge->next->next->origin->position;

    glm::vec3 normal = glm::cross(v1 - v0, v2 - v0);
    float area = glm::length(normal);

    if (area < 1e-6f) return;

    normal = glm::normalize(normal);

    float d = -glm::dot(normal, v0);
    glm::vec4 plane(normal, d);

    glm::mat4 quadric = glm::outerProduct(plane, plane);

    HalfEdge* currentEdge = edge;
    do {
        if (currentEdge && currentEdge->origin) {
            currentEdge->origin->quadric += quadric;
        }
        currentEdge = currentEdge->next;
    } while (currentEdge != edge);
}

void HalfEdge::computeEdgeCost(bool force = false) {
    HalfVertex* v1 = origin;
    HalfVertex* v2 = twin->origin;

    glm::mat4 Q = v1->quadric + v2->quadric;

    glm::mat4 Q_sub = Q;
    Q_sub[3] = glm::vec4(0, 0, 0, 1);

    glm::vec3 optimalPosition = (v1->position + v2->position) / 2.0f;

    glm::vec4 v_opt(optimalPosition, 1.0f);
    float costValue = glm::dot(v_opt, Q * v_opt);

    costValue = glm::max(costValue, 0.0f);

	this->cost = Cost{ costValue, optimalPosition };
}

void HalfEdge::LerpVertex(Vertex& vert) {
    HalfVertex& v1 = *this->origin;
    assert(this->twin != nullptr && "[QEM:ERROR] Twin edge is null.");
    assert(this->twin->origin != nullptr && "[QEM:ERROR] Twin's origin vertex is null.");
    HalfVertex& v2 = *this->twin->origin;

    vert.position = cost.optimalPosition;

    glm::vec3 lerpedNormal = glm::normalize(glm::mix(v1.normal, v2.normal, 0.5f));
    vert.normal = glm::normalize(lerpedNormal);

    vert.uv = glm::mix(v1.uv, v2.uv, 0.5f);
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

void HalfEdgeMesh::edgeCollapse(HalfEdge* edge) {
    assert(edge->isValid);
    HalfEdge*start, *cur;
	int counter = 0;

    HalfVertex* v1 = edge->origin;
	HalfVertex* v2 = edge->twin->origin;

    if (!v1->edge->isValid || !v2->edge->isValid) {
        edge->isValid = false;
        edge->twin->isValid = false;
        return;
    }

    if (edge->twin->next == edge->prev->twin) {
        // special case: Face back-to-back
        std::cout << "[SPECIAL CASE]: Face back-to-back:(" << edge->face << ")(" << edge->twin->face << ")" << std::endl;

        v2->position = (v1->position + v2->position) / 2.0f;

        HalfVertex* v3 = edge->prev->origin;
        if (v3->edge == edge->prev || v3->edge == edge->twin->prev) {
            std::cout << "TRIGGER SET V3 EDGE" << std::endl;
            v3->edge = edge->next->twin;
        }

        edge->next->twin->twin = edge->twin->prev->twin;
        edge->twin->prev->twin->twin = edge->next->twin;

		edge->isValid = false;
		edge->next->isValid = false;
		edge->prev->isValid = false;

		edge->twin->isValid = false;
		edge->twin->next->isValid = false;
		edge->twin->prev->isValid = false;

        edge->face->isValid = false;
        edge->twin->face->isValid = false;

        v1->isValid = false;
        if (v2->edge == edge->twin) {
            std::cout << "TRIGGER SET V2 EDGE" << std::endl;
            v2->edge = edge->twin->prev->twin;
        }
        //assert(v2->edge->isValid);
        //assert(v3->edge->isValid);

        //HalfEdgeMeshValidation();
        return;
    }

    cur = start = v2->edge; counter = 0;
    do {
        cur = cur->twin->next;
        if (counter++ > 1000) {
            // huge vertex loop no merge
            edge->isValid = false;
            edge->twin->isValid = false;
            return;
        }
    } while (cur != start);

    cur = start = v2->edge; counter = 0;
    do {
        if (cur->isValid) cur->origin = v1;
        cur = cur->twin->next;
    } while (cur != start);

    v1->position = (v1->position + v2->position) / 2.0f;

    if (v1->edge == edge || v1->edge == edge->twin->next) {
        std::cout << "TRIGGER SET V1 EDGE" << std::endl;
        v1->edge = edge->prev->twin;
        if (edge->prev->twin == edge->twin) {
            std::cout << "Detect Special Case:" << std::endl;
			std::cout << *(edge->face) << *(edge->twin->face) << std::endl;
			assert(false);
        }
    }

	HalfVertex* v1adj = edge->prev->origin;
	HalfVertex* v2adj = edge->twin->prev->origin;

    if (v1adj->edge == edge->prev) {
        std::cout << "TRIGGER SET V1ADJ EDGE" << std::endl;
        v1adj->edge = edge->next->twin;
    }
    if (v2adj->edge == edge->twin->prev) {
        std::cout << "TRIGGER SET V2ADJ EDGE" << std::endl;
        v2adj->edge = edge->twin->next->twin;
    }

    edge->next->twin->twin = edge->prev->twin;
    edge->prev->twin->twin = edge->next->twin;

    edge->twin->next->twin->twin = edge->twin->prev->twin;
    edge->twin->prev->twin->twin = edge->twin->next->twin;

    Face* face1 = edge->face;
    Face* face2 = edge->twin->face;

    cur = start = face1->edge; counter = 0;

    edge->isValid = false;
	edge->next->isValid = false;
	edge->prev->isValid = false;
	edge->twin->isValid = false;
    assert(edge->twin != v1->edge);
	edge->twin->next->isValid = false;
	edge->twin->prev->isValid = false;

    face1->isValid = false;
    face2->isValid = false;
    v2->isValid = false;

 //   assert(v1->isValid);
 //   assert(v1->edge->isValid);
 //   assert(v1adj->edge->isValid);
 //   assert(v2adj->edge->isValid);

 //   cur = start = v1->edge; counter = 0;
 //   std::unordered_set<HalfEdge*> setting;
 //   do {
 //       cur = cur->twin->next;
 //       assert(cur->isValid);
	//	if (setting.find(cur) != setting.end()) assert(false);
 //       setting.insert(cur);
 //       assert(counter++ < 1000);
 //   } while (cur != start);

	//setting.clear();
	//cur = start = v1adj->edge; counter = 0;
	//do {
	//	cur = cur->twin->next;
	//	assert(cur->isValid);
	//	if (setting.find(cur) != setting.end()) assert(false);
	//	setting.insert(cur);
	//	assert(counter++ < 1000);
	//} while (cur != start);

	//setting.clear();
	//cur = start = v2adj->edge; counter = 0;
	//do {
	//	cur = cur->twin->next;
	//	assert(cur->isValid);
	//	if (setting.find(cur) != setting.end()) assert(false);
	//	setting.insert(cur);
	//	assert(counter++ < 1000);
	//} while (cur != start);

	//HalfEdgeMeshValidation();
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

void HalfEdgeMesh::QEM(){
    int validFaces = m_faces.size();
	int targetFaces = validFaces / 2;
    while (validFaces > targetFaces) {
        initCostComputation();
        // generate the QEMQueue
        QEMQueue queue;
        for (HalfEdge* edge : m_edges) {
            if ( edge->isValid && edge < edge->twin) {
                queue.push(edge);
            }
        }
		std::cout << "Queue size: " << queue.size() << std::endl;
        edgeCollapse(queue.top());

		validFaces = 0;
        for (auto& face : m_faces) {
            if (face->isValid) validFaces++;
        }
    }
}















void processClusterGroup(
    EEdgePriorityQueue& pq,
    int& faceCount,
    int targetFaceCount,
    EMesh* mesh
);
void parallelProcessClusterGroups(
    std::vector<EEdgePriorityQueue>& clusterGroupQueues,
    std::vector<int>& clusterGroupFaceCount,
    const std::vector<int>& clusterGroupTargetFaceCount,
    EMesh* mesh
);

void EMesh::edgeCollapse(EEdge* edge){
    if(edge->isBoundary || edge->isFakeBoundary){
        edge->isDirty = true;
        return;
    }

    // merging v2 into v1 of edge
    EVertex* v1 = edge->vertices.v1;
    EVertex* v2 = edge->vertices.v2;
    if (v2->isFakeBoundary) {
		std::swap(v1, v2);
    }

    v1->position = edge->cost.optimalPosition;

    // remove faces from the adjacent vertices
    while(edge->faces.size() != 0){
        EFace* face = *edge->faces.begin();

        VertexPair pair1(face->vertices[0], face->vertices[1]);
        VertexPair pair2(face->vertices[1], face->vertices[2]);
        VertexPair pair3(face->vertices[2], face->vertices[0]);

        // find the edges from m_edgeMap
        if (m_edgeMap.find(pair1) != m_edgeMap.end()) {
            EEdge* edge1 = m_edgeMap[pair1];
            if (edge1->containFace(face)) edge1->removeFace(face);
        }
		if (m_edgeMap.find(pair2) != m_edgeMap.end()) {
			EEdge* edge2 = m_edgeMap[pair2];
			if (edge2->containFace(face)) edge2->removeFace(face);
		}
		if (m_edgeMap.find(pair3) != m_edgeMap.end()) {
			EEdge* edge3 = m_edgeMap[pair3];
			if (edge3->containFace(face)) edge3->removeFace(face);
		}
        
        face->isValid = false;
    }
    edge->isValid = false;

    // set all edges of v2 to v1
    while(v2->edges.size() != 0){
        EEdge* adjEdge = v2->popEdge();

        EVertex* v = adjEdge->vertices.v1 == v2 ? adjEdge->vertices.v2 : adjEdge->vertices.v1;
        if(v == v1) continue;

        m_edgeMap.erase(adjEdge->vertices);
        VertexPair pair(v1, v);

        // update the faces of the adjacent edge
        for(auto& face : adjEdge->faces){
            for(int i = 0; i < 3; ++i){
                if(face->vertices[i] == v2){
                    face->vertices[i] = v1;
                }
            }
        }
        if (m_edgeMap.find(pair) != m_edgeMap.end()) { // if exist, merge the duplicate edges by copying faces to the existing edge
			EEdge* existEdge = m_edgeMap[pair];
            while(EFace* popface = adjEdge->popFace()){
                existEdge->addFace(popface);
            }
            adjEdge->isValid = false;
        }
        else{
            m_edgeMap[pair] = adjEdge;
            adjEdge->vertices = pair;
            v1->addEdge(adjEdge);
        }
    }

    std::unordered_set<EFace*> faces;
    v1->quadric = glm::mat4(0.0f);
    for(auto egde : v1->edges){
        if (egde->isBoundary || !egde->isValid || edge->isDirty || edge->isFakeBoundary) continue;
        if(edge->vertices.v1 == v1) egde->vertices.v2->quadric = glm::mat4(0.0f);
        else egde->vertices.v1->quadric = glm::mat4(0.0f);
		if (egde->faces.size() == 0) continue;
        for(auto& face : egde->faces){
            faces.insert(face);
        }
    }

    for(auto& face : faces){
        face->updateFaceQuadric();
    }

	for (auto& edge : v1->edges) {
		edge->computeEdgeCost();
	}

}

void EFace::updateFaceQuadric(){
    glm::mat4 Q = glm::mat4(0.0f);
    glm::vec3 normal = glm::cross(
        vertices[1]->position - vertices[0]->position, 
        vertices[2]->position - vertices[0]->position
    );
    normal = glm::normalize(normal);
    float d = -glm::dot(normal, vertices[0]->position);
    glm::vec4 plane(normal, d);
    Q = glm::outerProduct(plane, plane);
    for (int i = 0; i < 3; ++i) {
        vertices[i]->quadric += Q;
    }
}

void EEdge::computeEdgeCost(){
    
    if(isBoundary || isFakeBoundary){
        cost = Cost{ std::numeric_limits<float>::infinity(), glm::vec3(0.0f) };
        return;
    }

    EVertex* v1 = vertices.v1;
    EVertex* v2 = vertices.v2;

    

    glm::mat4 Q = v1->quadric + v2->quadric;

    glm::mat4 Q_sub = Q;
    Q_sub[3] = glm::vec4(0, 0, 0, 1);

    glm::vec3 optimalPosition = (v1->position + v2->position) / 2.0f;
    if(v1->isFakeBoundary){
        optimalPosition = v1->position;
    }
    if(v2->isFakeBoundary){
        optimalPosition = v2->position;
    }


    glm::vec4 v_opt(optimalPosition, 1.0f);
    float costValue = glm::dot(v_opt, Q * v_opt);

    costValue = glm::max(costValue, 0.0f);

    cost = Cost{ costValue, optimalPosition };
}

void EMesh::QEM(float ratio){
    auto start = std::chrono::high_resolution_clock::now();

    for(EVertex* vertex : m_vertices){
        vertex->quadric = glm::mat4(0.0f);
    }

    // calculate quadric for each vertex by loop through faces
    std::unordered_map<const EFace*, size_t> faceToIndexMap;
    for (size_t i = 0; i < m_faces.size(); ++i) {
        faceToIndexMap[m_faces[i]] = i;
        m_faces[i]->updateFaceQuadric();
    }

    // calculate edge cost
    for (auto& [pair, edge] : m_edgeMap) {
        if (!edge->isValid) continue;
        edge->computeEdgeCost();
    }

    std::vector<EEdgePriorityQueue> clusterGroupQueues(m_clusterGroupCount);
    auto findClusterIndex = [this](int faceIdx) -> size_t {
        auto it = std::upper_bound(m_clusterOffsets.begin(), m_clusterOffsets.end(), faceIdx);
        return (it == m_clusterOffsets.begin()) ? 0 : (it - m_clusterOffsets.begin() - 1);
    };
    for (auto& edgePair : m_edgeMap) {
        EEdge* edge = edgePair.second;
        if (edge && edge->faces.size() == 2) {
            auto it = edge->faces.begin();
            EFace* face1 = *it++;
            EFace* face2 = *it;
            
            if (face1 && face2) {
                size_t clusterIdx1 = findClusterIndex(faceToIndexMap[face1]);
                size_t clusterIdx2 = findClusterIndex(faceToIndexMap[face2]);

                size_t groupIdx1 = m_clusterGroupResult[clusterIdx1];
                size_t groupIdx2 = m_clusterGroupResult[clusterIdx2];

                if (groupIdx1 != groupIdx2) {
                    edge->isFakeBoundary = true;
                    edge->vertices.v1->isFakeBoundary = true;
                    edge->vertices.v2->isFakeBoundary = true;
                } else {
                    clusterGroupQueues[groupIdx1].push(edge);
                }
            }
        }
    }

    std::vector<int> clusterGroupFaceCount(m_clusterGroupCount, 0);
    std::vector<int> clusterGroupTargetFaceCount(m_clusterGroupCount, 0);
    for (int i = 0; i < m_clusterOffsets.size(); i++) {
        if (i + 1 < m_clusterOffsets.size()) {
            clusterGroupFaceCount[m_clusterGroupResult[i]] += m_clusterOffsets[i + 1] - m_clusterOffsets[i];
        } else {
            clusterGroupFaceCount[m_clusterGroupResult[i]] += m_faces.size() - m_clusterOffsets[i];
        }
    }
    for (int i = 0; i < m_clusterGroupCount; i++) {
        clusterGroupTargetFaceCount[i] = clusterGroupFaceCount[i] * ratio;
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "QEM Initialization Time: " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << "s" << std::endl;

    start = std::chrono::high_resolution_clock::now();

    parallelProcessClusterGroups(clusterGroupQueues, clusterGroupFaceCount, clusterGroupTargetFaceCount, this);

    end = std::chrono::high_resolution_clock::now();
    std::cout << "QEM Collapse Time: " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << "s" << std::endl;
}

void processClusterGroup(
    EEdgePriorityQueue& pq,
    int& faceCount,
    int targetFaceCount,
    EMesh* mesh
) {
    while (faceCount > targetFaceCount) {
        if (faceCount % 1000 == 0) {
            std::cout << "Target Faces: " << targetFaceCount
                      << ", Current Faces: " << faceCount << std::endl;
        }

        EEdge* edge = pq.pop();
        if (edge) {
            mesh->edgeCollapse(edge);
            faceCount--;
        }
    }
}

void parallelProcessClusterGroups(
    std::vector<EEdgePriorityQueue>& clusterGroupQueues,
    std::vector<int>& clusterGroupFaceCount,
    const std::vector<int>& clusterGroupTargetFaceCount,
    EMesh* mesh
) {
    std::vector<std::thread> threads;

    for (size_t i = 0; i < clusterGroupQueues.size(); ++i) {
        threads.emplace_back(processClusterGroup,
                             std::ref(clusterGroupQueues[i]),
                             std::ref(clusterGroupFaceCount[i]),
                             clusterGroupTargetFaceCount[i],
                             mesh);
    }

    // 等待所有线程完成
    for (auto& t : threads) {
        t.join();
    }
}
