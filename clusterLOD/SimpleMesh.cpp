#include "SimpleMesh.hpp"

#include "cs488-framework/ObjFileDecoder.hpp"

#include "QEM/Simplify.hpp"

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
            for(int j = i+1; j < faces.size(); j++){
                faces[i]->adjacentFace.insert(faces[j]);
                faces[j]->adjacentFace.insert(faces[i]);
            }
        }
    }
}

void SimpleMesh::exportClusterGroup(const std::string& lodFolderPath) {
    if (!std::filesystem::exists(lodFolderPath)) {
        std::filesystem::create_directories(lodFolderPath);
    }

#pragma omp parallel for
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
#pragma omp critical
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



void SimpleMesh::splitterRecur(int startIdx, int endIdx, int depth){
    size_t triCount = endIdx - startIdx + 1;

    if(triCount <= MAX_TRI_IN_CLUSTER){
        #pragma omp critical
        m_clusterOffsets.push_back(startIdx);
        return;
    }
    
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
    idx_t numConstraints = 1;
    std::vector<idx_t> xadj(triCount + 1, 0);
    std::vector<idx_t> adjncy;

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

    for(int i = 0; i < triCount; i++){
        m_faces[i + startIdx].clusterId = partitionResult[i];
    }

   int ptr0 = 0;
   int ptr1 = triCount-1;
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

	xadj.clear();
	adjncy.clear();
	partitionResult.clear();


   size_t midIdx = lastZero + startIdx;
   #pragma omp parallel sections
   {
       #pragma omp section
       {splitterRecur(startIdx, midIdx, depth+1);}

       #pragma omp section
       {splitterRecur(midIdx+1, endIdx, depth+1);}
   }
}

void SimpleMesh::splitter(){
    m_clusterGroupOffsets.clear();
    splitterRecur(0, m_faces.size() - 1, 0);
    std::sort(m_clusterGroupOffsets.begin(), m_clusterGroupOffsets.end());
    m_clusterGroupOffsets.push_back(m_faces.size());
}

void SimpleMesh::grouper(){

}


void SimpleMesh::partition_loop(const std::string& objFilePath, const std::string& lodFolderPath){
    auto starttime = std::chrono::high_resolution_clock::now();
    importMesh(objFilePath);
    auto endtime = std::chrono::high_resolution_clock::now();
    std::cout << "Import Mesh Time: " 
              << std::fixed << std::setprecision(5)
              << std::chrono::duration<double>(endtime - starttime).count()
              << "s" << std::endl;

	m_clusterGroupOffsets.clear();
	m_clusterOffsets.clear();

    starttime = std::chrono::high_resolution_clock::now();
    splitter();
    endtime = std::chrono::high_resolution_clock::now();
    std::cout << "Splitter Time: " 
              << std::fixed << std::setprecision(5)
              << std::chrono::duration<double>(endtime - starttime).count()
              << "s" << std::endl;

    // starttime = std::chrono::high_resolution_clock::now();
    // exportClusterGroup(lodFolderPath);
    // endtime = std::chrono::high_resolution_clock::now();
    //     std::cout << "Export Cluster Group Time: " 
    //           << std::fixed << std::setprecision(5)
    //           << std::chrono::duration<double>(endtime - starttime).count()
    //           << "s" << std::endl;

    // QEM(0, m_faces.size() - 1, lodFolderPath, 0.5);
    starttime = std::chrono::high_resolution_clock::now();
    if (!std::filesystem::exists(lodFolderPath)) {
        std::filesystem::create_directories(lodFolderPath);
    }

    #pragma omp parallel for
    for (size_t clusterIndex = 0; clusterIndex < m_clusterGroupOffsets.size() - 1; ++clusterIndex) {
        size_t startIdx = m_clusterGroupOffsets[clusterIndex];
        size_t endIdx = m_clusterGroupOffsets[clusterIndex + 1] - 1;

        // construct file name
        std::string fileName = lodFolderPath + "/clusterGroup_" + std::to_string(clusterIndex) + ".obj";
        QEM(startIdx, endIdx, fileName, 0.5);
    }
    endtime = std::chrono::high_resolution_clock::now();
    std::cout << "QEM Cluster Group + Export Time: " 
            << std::fixed << std::setprecision(5)
            << std::chrono::duration<double>(endtime - starttime).count()
            << "s" << std::endl;
}


void SimpleMesh::exportMeshSimplifier(MeshSimplifier& simplifier, int startIdx, int endIdx) {
    simplifier.vertices.clear();
    simplifier.triangles.clear();
    simplifier.refs.clear();

    std::unordered_map<int, int> vertexIndexMap;
    int simplifierVertexIndex = 0;

    for (int i = startIdx; i <= endIdx; ++i) {
        const auto& face = m_faces[i];

        auto processVertex = [&](SimpleVertex* vertex) {
            int originalIndex = std::distance(&m_vertices[0], vertex);
            if (vertexIndexMap.find(originalIndex) == vertexIndexMap.end()) {
                MeshSimplifier::Vertex simplifiedVertex;
                simplifiedVertex.p = vec3f(vertex->position.x, vertex->position.y, vertex->position.z);
                simplifiedVertex.tstart = 0;
                simplifiedVertex.tcount = 0;
                simplifiedVertex.border = 0;

                vertexIndexMap[originalIndex] = simplifierVertexIndex;
                simplifier.vertices.push_back(simplifiedVertex);
                ++simplifierVertexIndex;
            }
        };

        processVertex(face.v1);
        processVertex(face.v2);
        processVertex(face.v3);
    }

    for (int i = startIdx; i <= endIdx && i < m_faces.size(); ++i) {
        const auto& face = m_faces[i];
        MeshSimplifier::Triangle simplifiedTriangle;

        simplifiedTriangle.v[0] = vertexIndexMap[std::distance(&m_vertices[0], face.v1)];
        simplifiedTriangle.v[1] = vertexIndexMap[std::distance(&m_vertices[0], face.v2)];
        simplifiedTriangle.v[2] = vertexIndexMap[std::distance(&m_vertices[0], face.v3)];

        glm::vec3 faceNormal = glm::normalize(glm::cross(
            face.v2->position - face.v1->position,
            face.v3->position - face.v1->position
        ));
        simplifiedTriangle.n = vec3f(faceNormal.x, faceNormal.y, faceNormal.z);

        simplifiedTriangle.uvs[0] = vec3f(face.v1->uv.x, face.v1->uv.y, 0.0);
        simplifiedTriangle.uvs[1] = vec3f(face.v2->uv.x, face.v2->uv.y, 0.0);
        simplifiedTriangle.uvs[2] = vec3f(face.v3->uv.x, face.v3->uv.y, 0.0);

        simplifiedTriangle.deleted = 0;
        simplifiedTriangle.dirty = 0;
        simplifiedTriangle.attr = MeshSimplifier::Attributes::TEXCOORD | MeshSimplifier::Attributes::NORMAL;
        simplifiedTriangle.material = -1;

        simplifier.triangles.push_back(simplifiedTriangle);
    }
}

float SimpleMesh::QEM(int start, int end, const std::string& objFilePath, float ratio = 0.5) {
    MeshSimplifier simplifier;
    exportMeshSimplifier(simplifier, start, end);

    int target_count = simplifier.triangles.size() >> 1;

    if (ratio > 1.0) ratio = 1.0;
    if (ratio <= 0.0) {
        std::cout << "Ratio must be BETWEEN zero and one." << std::endl;
        return EXIT_FAILURE;
    }
    target_count = round((float)simplifier.triangles.size() * ratio);

    if (target_count < 4) {
        std::cout << "Object will not survive such extreme decimation" << std::endl;
        return EXIT_FAILURE;
    }

    int startSize = simplifier.triangles.size();
    simplifier.simplify_mesh(target_count, 7.0, true);

    if (simplifier.triangles.size() >= startSize) {
        std::cout << "Unable to reduce mesh." << std::endl;
        return EXIT_FAILURE;
    }

    simplifier.write_obj((objFilePath).c_str());

    #pragma omp critical
    std::cout << "Output: " << simplifier.vertices.size() << " vertices, " << simplifier.triangles.size() 
              << " triangles (" << (float)simplifier.triangles.size() / (float)startSize 
              << " reduction; " << ((float)(clock() - start)) / CLOCKS_PER_SEC << " sec)"
              << "Total error: " << simplifier.total_error << std::endl;
    return simplifier.total_error;
}

void SimpleMesh::importMeshSimplifier(const MeshSimplifier& simplifier) {
    m_vertices.clear();
    m_faces.clear();

    for (const auto& simplifiedVertex : simplifier.vertices) {
        glm::vec3 position(simplifiedVertex.p.x, simplifiedVertex.p.y, simplifiedVertex.p.z);
        glm::vec3 normal(0.0f, 0.0f, 0.0f);
        glm::vec2 uv(0.0f, 0.0f);

        m_vertices.emplace_back(position, normal, uv);
    }

    for (const auto& simplifiedTriangle : simplifier.triangles) {
        if (simplifiedTriangle.deleted) continue;

        SimpleFace face(
            &m_vertices[simplifiedTriangle.v[0]],
            &m_vertices[simplifiedTriangle.v[1]],
            &m_vertices[simplifiedTriangle.v[2]]
        );

        face.v1->uv = glm::vec2(simplifiedTriangle.uvs[0].x, simplifiedTriangle.uvs[0].y);
        face.v2->uv = glm::vec2(simplifiedTriangle.uvs[1].x, simplifiedTriangle.uvs[1].y);
        face.v3->uv = glm::vec2(simplifiedTriangle.uvs[2].x, simplifiedTriangle.uvs[2].y);

        glm::vec3 faceNormal(simplifiedTriangle.n.x, simplifiedTriangle.n.y, simplifiedTriangle.n.z);
        face.v1->normal += faceNormal;
        face.v2->normal += faceNormal;
        face.v3->normal += faceNormal;

        m_faces.push_back(face);
    }

    for (auto& vertex : m_vertices) {
        if (glm::length(vertex.normal) > 0.0f) {
            vertex.normal = glm::normalize(vertex.normal);
        } else {
            vertex.normal = glm::vec3(0.0f, 1.0f, 0.0f);
        }
    }
}
