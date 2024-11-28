
// #include "HalfEdgeMesh.hpp"

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
