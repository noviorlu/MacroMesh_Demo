#include "../HalfEdgeMesh.hpp"
#include <iostream>
#include <cassert>

#include <glm/gtx/string_cast.hpp>


void QEMTestcase();

void QEMTestcase() {
    HalfEdgeMesh mesh;
    std::cout << "[DEBUG] Created HalfEdgeMesh." << std::endl;

    // 初始化顶点
    HalfVertex* v0 = new HalfVertex(glm::vec3(0, 0, 0), glm::vec3(0, 0, 1), glm::vec2(0, 0));
    HalfVertex* v1 = new HalfVertex(glm::vec3(1, 0, 0), glm::vec3(0, 0, 1), glm::vec2(1, 0));
    HalfVertex* v2 = new HalfVertex(glm::vec3(1, 1, 0), glm::vec3(0, 0, 1), glm::vec2(1, 1));
    HalfVertex* v3 = new HalfVertex(glm::vec3(0, 1, 0), glm::vec3(0, 0, 1), glm::vec2(0, 1));
    std::cout << "[DEBUG] Initialized 4 vertices." << std::endl;

    // 初始化边
    HalfEdge* e0 = new HalfEdge();
    HalfEdge* e1 = new HalfEdge();
    HalfEdge* e2 = new HalfEdge();
    HalfEdge* e3 = new HalfEdge();

    HalfEdge* e0_twin = new HalfEdge();
    HalfEdge* e1_twin = new HalfEdge();
    HalfEdge* e2_twin = new HalfEdge();
    HalfEdge* e3_twin = new HalfEdge();

    e0->origin = v0; e0->twin = e0_twin; e0_twin->twin = e0; e0_twin->origin = v1;
    e1->origin = v1; e1->twin = e1_twin; e1_twin->twin = e1; e1_twin->origin = v2;
    e2->origin = v2; e2->twin = e2_twin; e2_twin->twin = e2; e2_twin->origin = v3;
    e3->origin = v3; e3->twin = e3_twin; e3_twin->twin = e3; e3_twin->origin = v0;

    e0->next = e1; e0->prev = e3;
    e1->next = e2; e1->prev = e0;
    e2->next = e3; e2->prev = e1;
    e3->next = e0; e3->prev = e2;

    e0_twin->next = e3_twin; e0_twin->prev = e1_twin;
    e1_twin->next = e0_twin; e1_twin->prev = e2_twin;
    e2_twin->next = e1_twin; e2_twin->prev = e3_twin;
    e3_twin->next = e2_twin; e3_twin->prev = e0_twin;

    Face* f0 = new Face(e0);

    e0->face = f0; e1->face = f0; e2->face = f0; e3->face = f0;

    e0_twin->face = nullptr;
    e1_twin->face = nullptr;
    e2_twin->face = nullptr;
    e3_twin->face = nullptr;

    v0->edge = e0;
    v1->edge = e1;
    v2->edge = e2;
    v3->edge = e3;

    mesh.m_vertices = {v0, v1, v2, v3};
    std::cout << "[DEBUG] Added vertices to the mesh. Vertex count: " << mesh.m_vertices.size() << std::endl;

    mesh.m_edges = {e0, e1, e2, e3, e0_twin, e1_twin, e2_twin, e3_twin};
    std::cout << "[DEBUG] Added edges to the mesh. Edge count: " << mesh.m_edges.size() << std::endl;

    mesh.m_faces = {f0};
    std::cout << "[DEBUG] Added faces to the mesh. Face count: " << mesh.m_faces.size() << std::endl;

    mesh.HalfEdgeMeshValidation();
    mesh.HalfEdgeMeshPrint();
    std::cout << "[DEBUG] HalfEdgeMesh validation passed after initialization." << std::endl;

    mesh.initCostComputation();
    std::cout << "[DEBUG] Initialized cost computation." << std::endl;

    // 合并 e0 和 e0_twin
    std::cout << "[DEBUG] Merging edge e0 and its twin." << std::endl;
    mesh.mergeEdge(e0);

    // 检查网格结构是否仍然有效
    try {
        mesh.HalfEdgeMeshValidation();
        mesh.HalfEdgeMeshPrint();
        std::cout << "[DEBUG] HalfEdgeMesh validation passed after merging edge." << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "[ERROR] HalfEdgeMesh validation failed: " << e.what() << std::endl;
    }
}

int main() {
    std::cout << "Running QEMTest..." << std::endl;
    try {
        QEMTestcase();
    } catch (const std::exception& e) {
        std::cerr << "[ERROR] Exception caught: " << e.what() << std::endl;
        return 1;
    } catch (...) {
        std::cerr << "[ERROR] Unknown exception caught." << std::endl;
        return 1;
    }
    return 0;
}
