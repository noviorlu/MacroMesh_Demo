
// #include "HalfEdgeMesh.hpp"

// #include <unordered_set>

// /// @brief Compute the quadric matrix for a given face.
// /// @param v0 Position of the first vertex of the face.
// /// @param v1 Position of the second vertex of the face.
// /// @param v2 Position of the third vertex of the face.
// /// @return The computed quadric matrix.
// glm::mat4 computeFaceQuadric(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2) {
//     // Compute the normal of the face
//     glm::vec3 normal = glm::cross(v1 - v0, v2 - v0);
//     float area = glm::length(normal);

//     // Ignore faces with zero area
//     if (area < 1e-6f) {
//         return glm::mat4(0.0f);
//     }

//     normal = glm::normalize(normal);

//     // Compute the plane equation(ax + by + cz + d = 0)
//     float d = -glm::dot(normal, v0);
//     glm::vec4 plane(normal, d);

//     // Compute the quadric matrix from the plane equation
//     return glm::outerProduct(plane, plane);
// }

// /// @brief Compute the initial quadrics for each vertex in the mesh.
// void HalfEdgeMesh::computeInitialQuadrics() {
//     for (auto& face : m_faces) {
//         HalfEdge* edge = face.edge;

//         // Extract vertex positions
//         glm::vec3 v0 = edge->origin->position;
//         glm::vec3 v1 = edge->next->origin->position;
//         glm::vec3 v2 = edge->next->next->origin->position;

//         // Compute the quadric for the face
//         glm::mat4 quadric = computeFaceQuadric(v0, v1, v2);

//         // Add the quadric matrix to all vertices of the face
//         HalfEdge* currentEdge = edge;
//         do {
//             currentEdge->origin->quadric += quadric;
//             currentEdge = currentEdge->next;
//         } while (currentEdge != edge);
//     }
// }

// void HalfEdgeMesh::mergeVertices(HalfVertex* v1, HalfVertex* v2, const glm::vec3& newPosition) {
//     // Update v1's position and quadric
//     v1->position = newPosition;
//     v1->quadric += v2->quadric;

//     HalfEdge* startEdge = v2->edge;
//     HalfEdge* currentEdge = startEdge;

//     do {
//         currentEdge->origin = v1;
//         currentEdge = currentEdge->twin->next;
//     } while (currentEdge != startEdge);

//     // Remove redundant edges
//     for (auto it = m_edges.begin(); it != m_edges.end(); ) {
//         HalfEdge* edge = *it;
//         if (edge->origin == edge->twin->origin) {
//             delete edge->twin;
//             delete edge;
//             it = m_edges.erase(it);
//         } else {
//             ++it;
//         }
//     }

//     // Recompute quadrics for affected vertices
//     for (HalfEdge* edge : m_edges) {
//         if (edge->origin == v1) {
//             HalfEdge* startEdge = edge->face->edge;
//             glm::vec3 v0 = startEdge->origin->position;
//             glm::vec3 v1 = startEdge->next->origin->position;
//             glm::vec3 v2 = startEdge->next->next->origin->position;

//             // Compute the quadric for the face
//             glm::mat4 faceQuadric = computeFaceQuadric(v0, v1, v2);
//             edge->origin->quadric += faceQuadric;
//         }
//     }
// }

// void HalfEdgeMesh::computeEdgeCost(HalfEdge* edge, EdgeCost& edgeCost) {
//     HalfVertex* v1 = edge->origin;
//     HalfVertex* v2 = edge->twin->origin;

//     glm::mat4 Q = v1->quadric + v2->quadric;

//     glm::mat4 Q_sub = Q;
//     Q_sub[3] = glm::vec4(0, 0, 0, 1);
//     glm::vec3 optimalPosition;
//     bool valid = glm::inverse(Q_sub) != glm::mat4(0.0f);
//     if (valid) {
//         glm::vec4 v_optimal = glm::inverse(Q_sub) * glm::vec4(0, 0, 0, 1);
//         optimalPosition = glm::vec3(v_optimal);
//     } else {
//         optimalPosition = (v1->position + v2->position) / 2.0f;
//     }

//     glm::vec4 v_opt(optimalPosition, 1.0f);
//     float cost = glm::dot(v_opt, Q * v_opt);

//     edgeCost.cost = cost;
//     edgeCost.optimalPosition = optimalPosition;
//     edgeCost.edge = edge;
// }

// void HalfEdgeMesh::simplifyMesh(){
//     // 1. split the edges based on the face clustergroup
//     std::vector<std::unordered_set<HalfEdge*>> clusterGroupEdges(m_clusterGroupOffsets.size());
//     for (auto& edge : m_edges) {
//         HalfEdge* twin = edge->twin;
//         if (twin != nullptr && twin < edge) {
//             continue;
//         }

//         // Get face indices
//         size_t faceIdx = edge->face - &m_faces[0];
//         size_t twinFaceIdx = twin->face - &m_faces[0];

//         // Lambda to find the cluster group index
//         auto findClusterGroupIndex = [&](size_t faceIdx) -> size_t {
//             auto it = std::upper_bound(
//                 m_clusterGroupOffsets.begin(), 
//                 m_clusterGroupOffsets.end(), 
//                 faceIdx
//             );
//             return it == m_clusterGroupOffsets.begin() ? 0 : std::distance(m_clusterGroupOffsets.begin(), it) - 1;
//         };

//         // Find cluster group indices for both faces
//         size_t clusterGroupIdx1 = findClusterGroupIndex(faceIdx);
//         size_t clusterGroupIdx2 = findClusterGroupIndex(twinFaceIdx);

//         // Only add the edge if both faces belong to the same cluster group
//         if (clusterGroupIdx1 == clusterGroupIdx2) {
//             clusterGroupEdges[clusterGroupIdx1].insert(edge);
//         }
//     }

//     // 2. Compute the cost for each edge and create a priority queue
//     std::vector<QEMqueue> queues(clusterGroupEdges.size());
//     for (size_t i = 0; i < clusterGroupEdges.size(); ++i) {
//         for (auto edge : clusterGroupEdges[i]) {
//             EdgeCost edgeCost;
//             computeEdgeCost(edge, edgeCost);
//             queues[i].push(edgeCost);
//         }
//     }

//     // free clusterGroupEdges
//     clusterGroupEdges.clear();

//     // 3. Simplify the mesh based on the priority queue

// }



