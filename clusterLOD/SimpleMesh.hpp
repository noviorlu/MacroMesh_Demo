#pragma once
#include <glm/glm.hpp>

#include "common.hpp"

#include "cs488-framework/OpenGLImport.hpp"
#include "cs488-framework/ShaderProgram.hpp"

#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <string>

#include <omp.h>


glm::vec3 HSVtoRGB(float h, float s, float v);
glm::vec3 genRdColor();

class MeshSimplifier;

class SimpleVertex{
public:
    glm::vec3 position;
    glm::vec3 normal;
    glm::vec2 uv;

    SimpleVertex(glm::vec3 position, glm::vec3 normal, glm::vec2 uv) : position(position), normal(normal), uv(uv) {}
};

class SimpleFace{
public: 
    unsigned int sequenceId = 0;

    unsigned int v1;
    unsigned int v2;
    unsigned int v3;

    std::unordered_set<SimpleFace*> adjacentFace;
    SimpleFace(unsigned int v1, unsigned int v2, unsigned int v3, unsigned int sequenceId) : v1(v1), v2(v2), v3(v3), sequenceId(sequenceId) {
		adjacentFace = std::unordered_set<SimpleFace*>();
    }
};

class LodMesh;
class SimpleMesh {
public: 
    std::string m_name;
    std::unordered_map<glm::vec3, unsigned int> m_vertexMap;
    std::unordered_map<std::pair<unsigned int, unsigned int>, std::vector<SimpleFace*>> m_edgeMap;

    std::vector<SimpleVertex> m_vertices;
    std::vector<SimpleFace*> m_faces;

    unsigned int createVertex(glm::vec3 position, glm::vec3 normal, glm::vec2 uv); 
    void createEdge(unsigned int v1, unsigned int v2, SimpleFace* face);

    void importMesh(const std::string& objFilePath);
    void exportMesh(const std::vector<std::pair<int, int>>& indexRanges, const std::string& objFilePath);
    void exportMesh(const std::string& objFilePath);
    void exportClusterGroup(const std::string& lodFolderPath);

    void splitterRecur(unsigned int start, unsigned int end, int depth);
    void splitter();
    void grouperRecur(unsigned int start, unsigned int end, int depth);
    void grouper();

    void QEM(SimpleMesh* targetMesh, const std::string& lodFolderPath, float ratio);
    void exportMeshSimplifier(MeshSimplifier& simplifier, const std::vector<std::pair<int, int>>& indexRanges);
     void importMeshSimplifier(const MeshSimplifier& simplifier);

    static void partition_loop(LodMesh& lodMesh, const std::string& objFilePath, const std::string& lodFolderPath);

    std::vector<unsigned int> m_clusterOffsets;
    std::vector<unsigned int> m_clusterGroupOffsets;
    std::vector<float> m_clusterGroupErrors;

    class ClusterGroup;
    class Cluster {
    public:
        unsigned int startIdx;
        unsigned int endIdx;

        unsigned int sequenceId;

        float error = 0.0f;

        std::unordered_set<Cluster*> adjClusters;

        Cluster(
            unsigned int startIdx, unsigned int endIdx, 
            unsigned int sequenceId, float error)
             : startIdx(startIdx), endIdx(endIdx), 
        sequenceId(sequenceId), error(error) {}
    
        void uploadToGPU();
        void removeFromGPU();
        void draw(const ShaderProgram& shader) const;

        GLuint m_vbo;
        GLuint m_vao;
        GLuint m_ibo;
    };
    class ClusterGroup {
        public:
        float error = 0.0f;
        std::vector<Cluster*> m_clusterlist;

        std::vector<std::pair<int, int>> getClusterRanges() {
			std::vector<std::pair<int, int>> indexRanges;
			for (Cluster* cluster : m_clusterlist) {
				indexRanges.push_back(std::make_pair(cluster->startIdx, cluster->endIdx));
			}
			return indexRanges;
        }

        ClusterGroup(float error, const std::vector<Cluster*>& clusterlist) : error(error) {
            m_clusterlist = clusterlist;
        }
    };
    
    std::vector<Cluster*> m_clusters;
    std::vector<ClusterGroup*> m_clusterGroups;
};


class LodMesh{
public:
    std::vector<SimpleMesh*> lodMesh;
    void printLODInformation();
};
