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

#include "Mesh.hpp"


class MeshSimplifier;

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
    class ClusterGroup;
    struct IntermData{
        unsigned int startIdx;
        unsigned int endIdx;
        float error;
        ClusterGroup* childGroup;
        IntermData(unsigned int s, unsigned int e, float err, ClusterGroup* cg)
        : startIdx(s), endIdx(e), error(err), childGroup(cg) {}
    };
    typedef std::vector<IntermData> IntermDataList;
    class Cluster {
    public:
        unsigned int startIdx;
        unsigned int endIdx;

        unsigned int sequenceId;

        float error = 0.0f;

        ClusterGroup* childGroup = nullptr;


        std::unordered_set<Cluster*> adjClusters;

        Cluster(
            unsigned int startIdx, unsigned int endIdx, 
            unsigned int sequenceId, float error,
            ClusterGroup* childGroup
            
            )
        : startIdx(startIdx), endIdx(endIdx), 
        sequenceId(sequenceId), error(error), childGroup(childGroup) {}
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
public: 
    std::string m_name;
    std::unordered_map<glm::vec3, unsigned int> m_vertexMap;
    std::unordered_map<std::pair<unsigned int, unsigned int>, std::vector<SimpleFace*>> m_edgeMap;

    std::vector<Vertex> m_vertices;
    std::vector<SimpleFace*> m_faces;

    unsigned int createVertex(glm::vec3 position, glm::vec3 normal, glm::vec2 uv); 
    void createEdge(unsigned int v1, unsigned int v2, SimpleFace* face);

    void importMesh(const std::string& objFilePath);
    void exportMesh(const std::vector<std::pair<int, int>>& indexRanges, const std::string& objFilePath);
    void exportMesh(const std::string& objFilePath);
    void exportMesh(unsigned int startIdx, unsigned int endIdx, Mesh* mesh);

    void splitterRecur(std::vector<unsigned int>& clusterOffsets, unsigned int start, unsigned int end, int depth);
    void splitter(const SimpleMesh::IntermDataList& intermDataList);
    void splitter();
    void grouperRecur(std::vector<unsigned int>& clusterGroupOffsets, unsigned int start, unsigned int end, int depth);
    void grouper();

    SimpleMesh::IntermDataList QEM(SimpleMesh* targetMesh, const std::string& lodFolderPath, float ratio);
    void exportMeshSimplifier(MeshSimplifier& simplifier, const std::vector<std::pair<int, int>>& indexRanges);
     void importMeshSimplifier(const MeshSimplifier& simplifier);

    static void partition_loop(LodMesh& lodMesh, const std::string& objFilePath, const std::string& lodFolderPath);
};


class LodMesh{
public:
    std::vector<SimpleMesh*> lodMesh;
    void printLODInformation();

    void exportLodRuntimeMesh(LodRuntimeMesh& mesh);
};
