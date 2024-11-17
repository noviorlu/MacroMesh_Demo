#pragma once

#include <iostream>
#include <vector>
#include <memory>
#include "cs488-framework/GlErrorCheck.hpp"
#include "Mesh.hpp"

#define N 256

class Cluster : public Mesh {
public:
    float Error;
    
    glm::vec3 rdColor;

    Cluster(
        float Error, 
        const Mesh& ref, 
        const std::vector<unsigned int>& triIndices
    );

    ~Cluster() {}

	void draw(const ShaderProgram& shader) const override{
        CHECK_GL_ERRORS;
        shader.SetUniform3fv("material.kd", rdColor);
        CHECK_GL_ERRORS;
        Mesh::draw(shader);
    }
};

// class ClusterGroup {
// public:
//     float Error;
//     std::vector<std::shared_ptr<Cluster>> clusters;

//     void addCluster(std::shared_ptr<Cluster> cluster) {
//         clusters.push_back(cluster);
//     }
// };

void MeshSplitter(Mesh& mesh, int num_parts);

std::vector<std::vector<size_t>> 
BuildAdjacencyList(const std::vector<unsigned int>& m_indexData);