#include "SimpleMesh.hpp"

#include "cs488-framework/ObjFileDecoder.hpp"

#include <metis.h>
#include <filesystem>
#include <fstream>
#include <iostream>

#define MAX_TRI_IN_CLUSTER 256
#define MAX_CLUSTER_IN_CLUSTERGROUP 32
#define MAX_TRI_IN_CLUSTERGROUP (MAX_TRI_IN_CLUSTER * MAX_CLUSTER_IN_CLUSTERGROUP)

void SimpleMesh::importMesh(const std::string& objFilePath) {
	
    std::vector<glm::vec3> positions;
    std::vector<glm::vec3> normals;
    std::vector<glm::vec2> uvCoords;

    ObjFileDecoder::decode(objFilePath.c_str(), m_name, positions, normals, uvCoords);

    std::unordered_map<std::pair<SimpleVertex*, SimpleVertex*>, std::vector<SimpleFace*>> edgeMap;
    std::unordered_map<glm::vec3, SimpleVertex*> vertexMap;
    std::function<SimpleVertex*(int)> createVertex;
    if (!uvCoords.empty()) {
        createVertex = [&](int i) -> SimpleVertex* {
            if (vertexMap.find(positions[i]) == vertexMap.end()) {
                m_vertices.emplace_back(positions[i], normals[i], uvCoords[i]);
                vertexMap[positions[i]] = &m_vertices.back();
            }
            return vertexMap[positions[i]];
        };
    } else {
        createVertex = [&](int i) -> SimpleVertex* {
            if (vertexMap.find(positions[i]) == vertexMap.end()) {
                m_vertices.emplace_back(positions[i], normals[i], glm::vec2(0));
                vertexMap[positions[i]] = &m_vertices.back();
            }
            return vertexMap[positions[i]];
        };
    }
    auto createEdge = [&](SimpleVertex* v1, SimpleVertex* v2, SimpleFace* face) {
        if(v1 > v2) std::swap(v1, v2);
        auto edge = std::make_pair(v1, v2);
        if(edgeMap.find(edge) == edgeMap.end()){
            edgeMap[edge] = std::vector<SimpleFace*>();
        }
        edgeMap[edge].push_back(face);
    };

    m_faces.clear();
    m_vertices.clear();
    m_faces.reserve(positions.size() / 3);
    m_vertices.reserve(positions.size());

    for(int i = 0; i < positions.size(); i+=3){
        SimpleVertex *v1, *v2, *v3;
        v1 = createVertex(i);
        v2 = createVertex(i+1);
        v3 = createVertex(i+2);

        m_faces.push_back(SimpleFace(v1, v2, v3));

        createEdge(v1, v2, &m_faces.back());
        createEdge(v2, v3, &m_faces.back());
        createEdge(v3, v1, &m_faces.back());
    }

    // construct AdjacentFaces in all SimpleFace
    for(auto& edge : edgeMap){
        auto& faces = edge.second;
        for(int i = 0; i < faces.size(); i++){
            for(int j = 0; j < faces.size(); j++){
				if (i == j) continue;
                faces[i]->adjacentFace.push_back(faces[j]);
                faces[j]->adjacentFace.push_back(faces[i]);
            }
        }
    }
}

void SimpleMesh::splitterRecur(int startIdx, int endIdx, int depth){
	if (depth == 1) {
		m_clusterGroupOffsets.push_back(startIdx);
		return;
	}
    
    size_t triCount = endIdx - startIdx + 1;

    //if(triCount <= MAX_TRI_IN_CLUSTERGROUP){
    //    //#pragma omp critical
    //    m_clusterGroupOffsets.push_back(startIdx);
    //    return;
    //}
    
    std::vector<std::vector<size_t>> adj_list;
    adj_list.resize(triCount);

    for (int i = 0; i < triCount; i++) {
        for (auto adjFace : m_faces[i + startIdx].adjacentFace) {
            int faceIdx = adjFace - &m_faces[0];
			if (faceIdx >= startIdx && faceIdx <= endIdx)
				adj_list[i].push_back(faceIdx - startIdx);
        }
    }

    // METIS partitioning

    size_t totalEdges = 0;
    for (size_t i = 0; i < triCount; ++i) {
        totalEdges += adj_list[i].size();
    }

    idx_t numConstraints = 1;
    std::vector<idx_t> xadj(triCount + 1, 0);
    std::vector<idx_t> adjncy(totalEdges);

    for (size_t i = 0; i < triCount; ++i) {
        xadj[i + 1] = xadj[i] + adj_list[i].size();
        for (size_t neighbor : adj_list[i]) {
            adjncy.push_back(static_cast<idx_t>(neighbor));
        }
    }

    idx_t numParts = 2;
    std::vector<idx_t> partitionResult(triCount);
    idx_t edgeCut;

    int result = METIS_PartGraphRecursive(
        reinterpret_cast<idx_t*>(&triCount),           
        &numConstraints,        
        xadj.data(),            
        adjncy.data(), 
        NULL, NULL, NULL,
        &numParts,              
        NULL, NULL, NULL,      
        &edgeCut,               
        partitionResult.data()  
    );
    if (result != METIS_OK) { assert(false); }

    int ptr0 = 0;
    int ptr1 = partitionResult.size()-1;
    int lastZero = ptr0;
    while (ptr0 <= ptr1) {
        if (partitionResult[ptr0] == 1 && partitionResult[ptr1] == 0) {
            std::swap(partitionResult[ptr0], partitionResult[ptr1]);
            SimpleFace::swap(m_faces[ptr0 + startIdx], m_faces[ptr1 + startIdx]);
            lastZero = ptr0;
        }
        if (partitionResult[ptr0] == 0) {
            lastZero = ptr0;
            ptr0++;
        }
		if (partitionResult[ptr1] == 1) ptr1--;
    }

	adj_list.clear();
	xadj.clear();
	adjncy.clear();
	partitionResult.clear();


    size_t midIdx = lastZero;
    //#pragma omp parallel sections
    {
        //#pragma omp section
        {splitterRecur(startIdx, midIdx, depth+1);}

        //#pragma omp section
        {splitterRecur(midIdx+1, endIdx, depth+1);}
    }
}

void SimpleMesh::splitter(){
    m_clusterGroupOffsets.clear();
    splitterRecur(0, m_faces.size() - 1, 0);
    std::sort(m_clusterGroupOffsets.begin(), m_clusterGroupOffsets.end());
    m_clusterGroupOffsets.push_back(m_faces.size());
}

void SimpleMesh::exportClusterGroup(const std::string& lodFolderPath){
    if (!std::filesystem::exists(lodFolderPath)) {
        std::filesystem::create_directories(lodFolderPath);
    }

    // print clusterGroupResults
	std::cout << "ClusterGroupOffsets: ";
	for (size_t i = 0; i < m_clusterGroupOffsets.size(); ++i) {
		std::cout << m_clusterGroupOffsets[i] << " ";
	}
	std::cout << std::endl;

    //#pragma omp parallel for
    for (size_t clusterIndex = 0; clusterIndex < m_clusterGroupOffsets.size() - 1; ++clusterIndex) {
        size_t startIdx = m_clusterGroupOffsets[clusterIndex];
        size_t endIdx = m_clusterGroupOffsets[clusterIndex + 1] - 1;

        // construct file name
        std::string fileName = lodFolderPath + "/clusterGroup_" + std::to_string(clusterIndex) + ".obj";
		exportMesh(startIdx, endIdx, fileName);
    }
}

void SimpleMesh::exportMesh(int startIdx, int endIdx, const std::string& objFilePath) {
    std::ofstream outFile(objFilePath);
    if (!outFile.is_open()) {
        //#pragma omp critical
        std::cerr << "Failed to open file: " << objFilePath << std::endl;
        return;
    }

    int vertexCount = 1;
    std::unordered_map<SimpleVertex*, int> vertexMapping;
	std::vector<SimpleVertex*> vertices;
	vertexMapping.reserve(endIdx - startIdx + 1);
	vertices.reserve(endIdx - startIdx + 1);
    auto createVertexMapping = [&](SimpleVertex* v) {
        if (vertexMapping.find(v) == vertexMapping.end()) {
			vertices.push_back(v);
            vertexMapping[v] = vertexCount;
            vertexCount++;
        }
    };
    for (int i = startIdx; i <= endIdx; i++) {
        const SimpleFace& face = m_faces[i];
        createVertexMapping(face.v1);
        createVertexMapping(face.v2);
        createVertexMapping(face.v3);
    }

    for (SimpleVertex* vertex : vertices) {
        outFile << "v " << vertex->position.x << " " << vertex->position.y << " " << vertex->position.z << std::endl;
        outFile << "vn " << vertex->normal.x << " " << vertex->normal.y << " " << vertex->normal.z << std::endl;
        outFile << "vt " << vertex->uv.x << " " << vertex->uv.y << std::endl;
    }

    for (int i = startIdx; i <= endIdx; i++) {
        const SimpleFace& face = m_faces[i];
        int v1 = vertexMapping[face.v1];
        int v2 = vertexMapping[face.v2];
        int v3 = vertexMapping[face.v3];
        outFile << "f " << v1 << "/" << v1 << "/" << v1 << " " << v2 << "/" << v2 << "/" << v2 << " " << v3 << "/" << v3 << "/" << v3 << std::endl;
    }
}

void SimpleMesh::exportMesh(const std::string& objFilePath, int clusterId) {
    std::ofstream outFile(objFilePath);
    if (!outFile.is_open()) {
        std::cerr << "Failed to open file: " << objFilePath << std::endl;
        return;
    }

	for (const auto& vertex : m_vertices) {
		outFile << "v " << vertex.position.x << " " << vertex.position.y << " " << vertex.position.z << std::endl;
		outFile << "vn " << vertex.normal.x << " " << vertex.normal.y << " " << vertex.normal.z << std::endl;
		outFile << "vt " << vertex.uv.x << " " << vertex.uv.y << std::endl;
	}

	for (const auto& face : m_faces) {
		if (face.clusterId != clusterId) {
			continue;
		}
		int v1 = face.v1 - &m_vertices[0] + 1;
		int v2 = face.v2 - &m_vertices[0] + 1;
		int v3 = face.v3 - &m_vertices[0] + 1;
		outFile << "f " << v1 << "/" << v1 << "/" << v1 << " " << v2 << "/" << v2 << "/" << v2 << " " << v3 << "/" << v3 << "/" << v3 << std::endl;
	}
}

void SimpleMesh::partition_loop(const std::string& objFilePath, const std::string& lodFolderPath){
    importMesh(objFilePath);
    splitter();
    exportClusterGroup(lodFolderPath);
	//exportMesh(0, m_faces.size() - 1, lodFolderPath + "/clusterGroup_0.obj");
}
