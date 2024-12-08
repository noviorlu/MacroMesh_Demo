#include "SimpleMesh.hpp"

#include "cs488-framework/ObjFileDecoder.hpp"

#pragma push_macro("min")
#pragma push_macro("max")
#undef min
#undef max
#include "QEM/Simplify.hpp"
#pragma pop_macro("min")
#pragma pop_macro("max")

#include <metis.h>
#include <filesystem>
#include <functional>
#include <fstream>
#include <iostream>

#define MAX_TRI_IN_CLUSTER 256
#define MAX_CLUSTER_IN_CLUSTERGROUP 32
#define MAX_TRI_IN_CLUSTERGROUP (MAX_TRI_IN_CLUSTER * MAX_CLUSTER_IN_CLUSTERGROUP)

glm::vec3 HSVtoRGB(float h, float s, float v) {
    h = fmod(h, 1.0f) * 6.0f;
    int i = static_cast<int>(floor(h));
    float f = h - i;
    float p = v * (1.0f - s);
    float q = v * (1.0f - f * s);
    float t = v * (1.0f - (1.0f - f) * s);

    switch (i) {
    case 0: return glm::vec3(v, t, p);
    case 1: return glm::vec3(q, v, p);
    case 2: return glm::vec3(p, v, t);
    case 3: return glm::vec3(p, q, v);
    case 4: return glm::vec3(t, p, v);
    case 5: return glm::vec3(v, p, q);
    default: return glm::vec3(0.0f, 0.0f, 0.0f);
    }
}

glm::vec3 genRdColor(){
    // generate a random hue (h) in [0, 1)
    float h = static_cast<float>(rand()) / RAND_MAX;
    float s = 0.8f;
    float v = 0.8f;

    // convert HSV to RGB
    return HSVtoRGB(h, s, v);
}


unsigned int SimpleMesh::createVertex(glm::vec3 position, glm::vec3 normal, glm::vec2 uv) {
    if (m_vertexMap.find(position) == m_vertexMap.end()) {
        m_vertices.emplace_back(position, normal, uv);
        m_vertexMap[position] = m_vertices.size() - 1;
    }
    return m_vertexMap[position];
}

void SimpleMesh::createEdge(unsigned int v1, unsigned int v2, SimpleFace* face) {
    if (v1 > v2) std::swap(v1, v2);
    auto edge = std::make_pair(v1, v2);
    if (m_edgeMap.find(edge) == m_edgeMap.end()) {
        m_edgeMap[edge] = std::vector<SimpleFace*>();
    }
    m_edgeMap[edge].push_back(face);
}

void SimpleMesh::importMesh(const std::string& objFilePath) {
    std::vector<glm::vec3> positions;
    std::vector<glm::vec3> normals;
    std::vector<glm::vec2> uvCoords;
    ObjFileDecoder::decode(objFilePath.c_str(), m_name, positions, normals, uvCoords);
    
    std::function<unsigned int(int)> createVertex;
    if (!uvCoords.empty()) {
        createVertex = [&](int i) -> unsigned int {
            return SimpleMesh::createVertex(positions[i], normals[i], uvCoords[i]);
        };
    } else {
        createVertex = [&](int i) -> unsigned int {
            return SimpleMesh::createVertex(positions[i], normals[i], glm::vec2(0.0f));
        };
    }

    m_faces.clear();
    m_vertices.clear();
    m_faces.reserve(positions.size() / 3);
    m_vertices.reserve(positions.size());

    for(int i = 0; i < positions.size(); i+=3){
        unsigned int v1, v2, v3;
        v1 = createVertex(i);
        v2 = createVertex(i+1);
        v3 = createVertex(i+2);

        SimpleFace* face = new SimpleFace(v1, v2, v3, i / 3);
        m_faces.push_back(face);

        createEdge(v1, v2, face);
        createEdge(v2, v3, face);
        createEdge(v3, v1, face);
    }

    // construct AdjacentFaces in all SimpleFace
    for(auto& edge : m_edgeMap){
        auto& faces = edge.second;
        for(int i = 0; i < faces.size(); i++){
            for(int j = i+1; j < faces.size(); j++){
                faces[i]->adjacentFace.insert(faces[j]);
                faces[j]->adjacentFace.insert(faces[i]);
            }
        }
    }

    m_vertexMap.clear();
    m_edgeMap.clear();

    m_vertices.shrink_to_fit();
    m_faces.shrink_to_fit();

    m_clusterGroupOffsets.push_back(0);
    m_clusterGroupOffsets.push_back(m_faces.size());

    m_clusterGroupErrors.push_back(0);
    m_clusterGroups.push_back(new ClusterGroup(0, std::vector<Cluster*>()));
}

void SimpleMesh::exportMesh(const std::vector<std::pair<int, int>>& indexRanges, const std::string& objFilePath) {
    std::ofstream outFile(objFilePath);
    if (!outFile.is_open()) {
        #pragma omp critical
        std::cerr << "Failed to open file: " << objFilePath << std::endl;
        return;
    }
    int triSize = 0;
    for(auto& range : indexRanges) triSize += range.second - range.first + 1;
    

    int vertexCount = 1;
    std::unordered_map<unsigned int, int> vertexMapping;
    std::vector<unsigned int> vertices;
    vertexMapping.reserve(triSize * 3);
    vertices.reserve(triSize * 3);
    auto createVertexMapping = [&](unsigned int v) {
        if (vertexMapping.find(v) == vertexMapping.end()) {
            vertices.push_back(v);
            vertexMapping[v] = vertexCount;
            vertexCount++;
        }
    };
    
    for(auto&range : indexRanges)
    {
        int startIdx = range.first;
        int endIdx = range.second;
        for (int i = startIdx; i <= endIdx; i++) {
            SimpleFace* face = m_faces[i];
            createVertexMapping(face->v1);
            createVertexMapping(face->v2);
            createVertexMapping(face->v3);
        }
    }

    for (unsigned int vertexId : vertices) {
        SimpleVertex* vertex = &m_vertices[vertexId];
        outFile << "v " << vertex->position.x << " " << vertex->position.y << " " << vertex->position.z << std::endl;
        outFile << "vn " << vertex->normal.x << " " << vertex->normal.y << " " << vertex->normal.z << std::endl;
        outFile << "vt " << vertex->uv.x << " " << vertex->uv.y << std::endl;
    }

    for(auto&range : indexRanges)
    {
        int startIdx = range.first;
        int endIdx = range.second;
        for (int i = startIdx; i <= endIdx; i++) {
            SimpleFace* face = m_faces[i];
            int v1 = vertexMapping[face->v1];
            int v2 = vertexMapping[face->v2];
            int v3 = vertexMapping[face->v3];
            outFile << "f " << v1 << "/" << v1 << "/" << v1 << " " << v2 << "/" << v2 << "/" << v2 << " " << v3 << "/" << v3 << "/" << v3 << std::endl;
        }
    }
}

void SimpleMesh::exportMesh(const std::string& objFilePath) {
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

    for (SimpleFace* face : m_faces) {
        int v1 = face->v1 + 1;
        int v2 = face->v2 + 1;
        int v3 = face->v3 + 1;
        outFile << "f " << v1 << "/" << v1 << "/" << v1 << " " << v2 << "/" << v2 << "/" << v2 << " " << v3 << "/" << v3 << "/" << v3 << std::endl;
    }
}

void SimpleMesh::exportClusterGroup(const std::string& lodFolderPath) {
    if (!std::filesystem::exists(lodFolderPath)) {
        std::filesystem::create_directories(lodFolderPath);
    }
    #pragma omp parallel for
    for(int i = 0; i < m_clusterGroupOffsets.size() - 1; i++){
        size_t startIdx = m_clusterGroupOffsets[i];
        size_t endIdx = m_clusterGroupOffsets[i+1] - 1;

        // construct file name
        std::string fileName = lodFolderPath + "/clusterGroup_" + std::to_string(i) + ".obj";
        
        std::vector<std::pair<int, int>> indexRanges;
        for(int j = startIdx; j <= endIdx; j++){
            Cluster* cluster = m_clusters[j];
            indexRanges.push_back(std::make_pair(cluster->startIdx, cluster->endIdx));
        }
        exportMesh(indexRanges, fileName);
    }
}

void SimpleMesh::splitterRecur(unsigned int startIdx, unsigned int endIdx, int depth){
    size_t triCount = endIdx - startIdx + 1;
     if(triCount <= MAX_TRI_IN_CLUSTER){
    //if(depth == 3){
        #pragma omp critical
        m_clusterOffsets.push_back(startIdx);
        return;
    }
    
    std::vector<std::vector<size_t>> adj_list;
    adj_list.resize(triCount);
    for(int i = 0; i < triCount; ++i){
        SimpleFace* fc = m_faces[i+startIdx];
        for(SimpleFace* adjfc : fc->adjacentFace){
            if(adjfc->sequenceId > endIdx || adjfc->sequenceId < startIdx) continue;
            adj_list[i].push_back(adjfc->sequenceId - startIdx);
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

   int ptr0 = 0;
   int ptr1 = triCount-1;
   int lastZero = ptr0;
    while (ptr0 <= ptr1) {
        if (partitionResult[ptr0] == 1 && partitionResult[ptr1] == 0) {
            std::swap(partitionResult[ptr0], partitionResult[ptr1]);
            std::swap(m_faces[ptr0 + startIdx]->sequenceId, m_faces[ptr1 + startIdx]->sequenceId);
            std::swap(m_faces[ptr0 + startIdx], m_faces[ptr1 + startIdx]);
        }
        if (partitionResult[ptr0] == 0) {
            lastZero = ptr0;
            ptr0++;
        }
        if (partitionResult[ptr1] == 1){
            ptr1--;
        }
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

int findLowerBoundIndex(const std::vector<unsigned int>& mylist, int idx) {
    int left = 0;
    int right = mylist.size();
    int mid;

    while (left < right) {
        mid = left + (right - left) / 2;
        if (mylist[mid] < idx) {
            left = mid + 1;
        }
        else {
            right = mid;
        }
    }

    if (left > 0 && mylist[left - 1] <= idx) {
        return left - 1;
    }

    return left;
}

void SimpleMesh::splitter(){
    m_clusterOffsets.clear();
    #pragma omp parallel for
    for (int i = 0; i < m_clusterGroupOffsets.size() - 1; i++) {
		int startIdx = m_clusterGroupOffsets[i];
		int endIdx = m_clusterGroupOffsets[i + 1] - 1;

        splitterRecur(startIdx, endIdx, 0);
    }
    std::sort(m_clusterOffsets.begin(), m_clusterOffsets.end());
    m_clusterOffsets.push_back(m_faces.size());

    m_clusters.resize(m_clusterOffsets.size() - 1);

    #pragma omp parallel for
    for(int i = 0; i < m_clusterOffsets.size() - 1; i++){
        unsigned int startIdx = m_clusterOffsets[i];
        unsigned int endIdx = m_clusterOffsets[i+1] - 1;

        // find clusterGroup error and assign to cluster
        int idx = findLowerBoundIndex(m_clusterGroupOffsets, startIdx);
        Cluster* cluster = new Cluster(startIdx, endIdx, i, m_clusterGroupErrors[idx]);
        m_clusters[i] = cluster;

        for(int j = startIdx; j <= endIdx; j++){
            m_faces[j]->sequenceId = i;
        }
    }

    #pragma omp parallel for
    for(Cluster* cluster : m_clusters){
        for(int i = cluster->startIdx; i <= cluster->endIdx; i++){
            SimpleFace* face = m_faces[i];
            for(SimpleFace* adjFace : face->adjacentFace){
                if(adjFace->sequenceId != face->sequenceId){
                    cluster->adjClusters.insert(m_clusters[adjFace->sequenceId]);
                }
            }
        }
    }
}

void SimpleMesh::grouperRecur(unsigned int startIdx, unsigned int endIdx, int depth){
    size_t clusterCount = endIdx - startIdx + 1;
    
    //if(clusterCount <= 1){
     if(clusterCount <= MAX_CLUSTER_IN_CLUSTERGROUP){
        #pragma omp critical
        m_clusterGroupOffsets.push_back(startIdx);
        return;
    }

    std::vector<std::vector<size_t>> adj_list;
    adj_list.resize(clusterCount);
    for(int i = 0; i < clusterCount; ++i){
        Cluster* cluster = m_clusters[i+startIdx];
        for(Cluster* adjCluster : cluster->adjClusters){
            if(adjCluster->sequenceId > endIdx || adjCluster->sequenceId < startIdx) continue;
            adj_list[i].push_back(adjCluster->sequenceId - startIdx);
        }
    }

     // METIS partitioning
    idx_t numConstraints = 1;
    std::vector<idx_t> xadj(clusterCount + 1, 0);
    std::vector<idx_t> adjncy;

    for (size_t i = 0; i < clusterCount; ++i) {
        xadj[i + 1] = xadj[i] + adj_list[i].size();
        for (size_t neighbor : adj_list[i]) {
            adjncy.push_back(static_cast<idx_t>(neighbor));
        }
    }

    idx_t numParts = 2;
    std::vector<idx_t> partitionResult(clusterCount);
    idx_t edgeCut;

    int result = METIS_PartGraphRecursive(
        reinterpret_cast<idx_t*>(&clusterCount),           
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
   int ptr1 = clusterCount-1;
   int lastZero = ptr0;
    while (ptr0 <= ptr1) {
        if (partitionResult[ptr0] == 1 && partitionResult[ptr1] == 0) {
            std::swap(partitionResult[ptr0], partitionResult[ptr1]);
            std::swap(m_clusters[ptr0 + startIdx]->sequenceId, m_clusters[ptr1 + startIdx]->sequenceId);
            std::swap(m_clusters[ptr0 + startIdx], m_clusters[ptr1 + startIdx]);
        }
        if (partitionResult[ptr0] == 0) {
            lastZero = ptr0;
            ptr0++;
        }
        if (partitionResult[ptr1] == 1){
            ptr1--;
        }
    }

    xadj.clear();
    adjncy.clear();
    partitionResult.clear();
    size_t midIdx = lastZero + startIdx;

    #pragma omp parallel sections
    {
        #pragma omp section
        {grouperRecur(startIdx, midIdx, depth+1);}

        #pragma omp section
        {grouperRecur(midIdx+1, endIdx, depth+1);}
    }
}

void SimpleMesh::grouper(){
    m_clusterGroupOffsets.clear();
    m_clusterGroupErrors.clear();

    grouperRecur(0, m_clusters.size() - 1, 0);
    
    std::sort(m_clusterGroupOffsets.begin(), m_clusterGroupOffsets.end());
    m_clusterGroupOffsets.push_back(m_clusters.size());

    m_clusterGroups.clear();
    m_clusterGroups.resize(m_clusterGroupOffsets.size() - 1);

    #pragma omp parallel for
    for(int i = 0; i < m_clusterGroupOffsets.size() - 1; i++){
        unsigned int startIdx = m_clusterGroupOffsets[i];
        unsigned int endIdx = m_clusterGroupOffsets[i+1] - 1;
        
        float maxError = 0.0f;
        std::vector<Cluster*> clusters;
        for(int j = startIdx; j <= endIdx; j++){
            clusters.push_back(m_clusters[j]);
            maxError = m_clusters[j]->error > maxError ? m_clusters[j]->error : maxError;
        }
        ClusterGroup* clusterGroup = new ClusterGroup(maxError, clusters);
        m_clusterGroups[i] = clusterGroup;
    }
}

void SimpleMesh::partition_loop(LodMesh& lod, const std::string& objFilePath, const std::string& lodFolderPath) {
    std::vector<SimpleMesh*>& lodMesh = lod.lodMesh;
    
    lodMesh.resize(9);
	lodMesh[0] = new SimpleMesh();

    auto starttime = std::chrono::high_resolution_clock::now();
    lodMesh[0]->importMesh(objFilePath);
    auto endtime = std::chrono::high_resolution_clock::now();
    std::cout << "Import Mesh Time: "
        << std::fixed << std::setprecision(5)
        << std::chrono::duration<double>(endtime - starttime).count()
        << "s" << std::endl;

    for (int i = 0;; i++) {
        SimpleMesh* srcMesh = lodMesh[i];
        auto folderPath = lodFolderPath + "_" + std::to_string(i);
        if (!std::filesystem::exists(folderPath)) {
            std::filesystem::create_directories(folderPath);
        }
        
        starttime = std::chrono::high_resolution_clock::now();
        srcMesh->splitter();
        endtime = std::chrono::high_resolution_clock::now();
        std::cout << "Splitter Time: "
            << std::fixed << std::setprecision(5)
            << std::chrono::duration<double>(endtime - starttime).count()
            << "s" << std::endl;

        starttime = std::chrono::high_resolution_clock::now();
        srcMesh->grouper();
        endtime = std::chrono::high_resolution_clock::now();
        std::cout << "Grouper Time: "
            << std::fixed << std::setprecision(5)
            << std::chrono::duration<double>(endtime - starttime).count()
            << "s" << std::endl;

        if (i >= 8 || srcMesh->m_faces.size() < MAX_TRI_IN_CLUSTER) break;
        lodMesh[i + 1] = new SimpleMesh();
        SimpleMesh* targetMesh = lodMesh[i + 1];
        starttime = std::chrono::high_resolution_clock::now();
        srcMesh->QEM(targetMesh, folderPath, 0.5);
        endtime = std::chrono::high_resolution_clock::now();
        std::cout << "QEM Cluster Group + Export Time: "
            << std::fixed << std::setprecision(5)
            << std::chrono::duration<double>(endtime - starttime).count()
            << "s" << std::endl;
    }
}


void SimpleMesh::exportMeshSimplifier(MeshSimplifier& simplifier, const std::vector<std::pair<int, int>>& indexRanges) {
    simplifier.vertices.clear();
    simplifier.triangles.clear();
    simplifier.refs.clear();

    std::unordered_map<int, int> vertexIndexMap;
    int simplifierVertexIndex = 0;
    auto createVertex = [&](int vertIdx) {
        if (vertexIndexMap.find(vertIdx) == vertexIndexMap.end()) {
			SimpleVertex& vertex = m_vertices[vertIdx];
            
            MeshSimplifier::Vertex simplifiedVertex;
            simplifiedVertex.p = vec3f(vertex.position.x, vertex.position.y, vertex.position.z);
            simplifiedVertex.tstart = 0;
            simplifiedVertex.tcount = 0;
            simplifiedVertex.border = 0;

            vertexIndexMap[vertIdx] = simplifierVertexIndex;
            simplifier.vertices.push_back(simplifiedVertex);
            ++simplifierVertexIndex;
        }
    };

    for (auto& range : indexRanges)
    {
        int startIdx = range.first;
        int endIdx = range.second;
        for (int i = startIdx; i <= endIdx; ++i) {
            SimpleFace* face = m_faces[i];
            createVertex(face->v1);
            createVertex(face->v2);
            createVertex(face->v3);
        }
    }

    for(auto& range : indexRanges) 
    {
        int startIdx = range.first;
        int endIdx = range.second;
        for (int i = startIdx; i <= endIdx && i < m_faces.size(); ++i) {
            const auto& face = m_faces[i];
            MeshSimplifier::Triangle simplifiedTriangle;

            simplifiedTriangle.v[0] = vertexIndexMap[face->v1];
            simplifiedTriangle.v[1] = vertexIndexMap[face->v2];
            simplifiedTriangle.v[2] = vertexIndexMap[face->v3];

            SimpleVertex& v1 = m_vertices[face->v1];
            SimpleVertex& v2 = m_vertices[face->v2];
            SimpleVertex& v3 = m_vertices[face->v3];
            glm::vec3 faceNormal = glm::normalize(glm::cross(
                v2.position - v1.position,
                v3.position - v1.position
            ));
            simplifiedTriangle.n = vec3f(faceNormal.x, faceNormal.y, faceNormal.z);

            simplifiedTriangle.deleted = 0;
            simplifiedTriangle.dirty = 0;
            simplifiedTriangle.attr = MeshSimplifier::Attributes::TEXCOORD | MeshSimplifier::Attributes::NORMAL;
            simplifiedTriangle.material = -1;

            simplifier.triangles.push_back(simplifiedTriangle);
        }
    }
}

void SimpleMesh::QEM(SimpleMesh* targetMesh, const std::string& lodFolderPath, float ratio = 0.5) {
    targetMesh->m_clusterGroupErrors.reserve(m_clusterGroups.size());
    targetMesh->m_clusterGroupOffsets.reserve(m_clusterGroups.size()+1);

    targetMesh->m_vertices.reserve(m_vertices.size());
    targetMesh->m_vertexMap.reserve(m_vertices.size());
    targetMesh->m_faces.reserve(m_faces.size() * ratio);
    targetMesh->m_edgeMap.reserve(m_vertices.size() * 3);

    #pragma omp parallel for
	for (int i = 0; i < m_clusterGroups.size(); i++) {
		ClusterGroup* clusterGroup = m_clusterGroups[i];
        MeshSimplifier simplifier;
        
        std::vector<std::pair<int, int>> indexRanges = clusterGroup->getClusterRanges();
        exportMeshSimplifier(simplifier, indexRanges);
    
        int target_count = simplifier.triangles.size() >> 1;

        ratio = std::clamp(ratio, 0.0f, 1.0f);
        target_count = round((float)simplifier.triangles.size() * ratio);

        if (target_count >= 4){
            int startSize = simplifier.triangles.size();
            simplifier.simplify_mesh(target_count, 7.0, true);
        }
        
        std::string objFilePath = lodFolderPath + "/clusterGroup_QEM_" + std::to_string(i) + ".obj";
        simplifier.write_obj((objFilePath).c_str());

        clusterGroup->error += simplifier.total_error;
        #pragma omp critical
        {
            targetMesh->m_clusterGroupErrors.push_back(clusterGroup->error);
            targetMesh->m_clusterGroupOffsets.push_back(targetMesh->m_faces.size());
            targetMesh->importMeshSimplifier(simplifier);
        }
    }

    // construct AdjacentFaces in all SimpleFace
    for(auto& edge : targetMesh->m_edgeMap){
        auto& faces = edge.second;
        for(int i = 0; i < faces.size(); i++){
            for(int j = i+1; j < faces.size(); j++){
                faces[i]->adjacentFace.insert(faces[j]);
                faces[j]->adjacentFace.insert(faces[i]);
            }
        }
    }

    targetMesh->m_clusterGroupOffsets.push_back(targetMesh->m_faces.size());

    targetMesh->m_vertexMap.clear();
    targetMesh->m_edgeMap.clear();

    targetMesh->m_vertices.shrink_to_fit();
    targetMesh->m_faces.shrink_to_fit();

    std::cout << "Reuction from: " << this->m_faces.size() << " to " << targetMesh->m_faces.size() << std::endl;
}

void SimpleMesh::importMeshSimplifier(const MeshSimplifier& simplifier) {
    int vertSize = simplifier.vertices.size();
    int triSize = simplifier.triangles.size();

    std::vector<glm::vec3> vertices(vertSize);
    std::vector<glm::vec3> vertex_normals(vertSize);
    std::vector<glm::vec2> uvs(vertSize);
    std::vector<bool> usedVertices(vertSize, false);

    for(int i = 0; i < vertSize; i++){
        vertices[i] = glm::vec3(simplifier.vertices[i].p.x, simplifier.vertices[i].p.y, simplifier.vertices[i].p.z);
    }

    for(int i = 0; i < vertSize; i++) {
        vertex_normals[i] = glm::vec3(0.0, 0.0, 0.0);
    }

    for(int i = 0; i < triSize; i++){
        if(simplifier.triangles[i].deleted) continue;
        glm::vec3 face_normal = glm::vec3(simplifier.triangles[i].n.x, simplifier.triangles[i].n.y, simplifier.triangles[i].n.z);

        int v1 = simplifier.triangles[i].v[0];
        int v2 = simplifier.triangles[i].v[1];
        int v3 = simplifier.triangles[i].v[2];

        vertex_normals[v1] += face_normal;
        vertex_normals[v2] += face_normal;
        vertex_normals[v3] += face_normal;

        usedVertices[v1] = true;
        usedVertices[v2] = true;
        usedVertices[v3] = true;

        uvs[v1] = glm::vec2(simplifier.triangles[i].uvs[0].x, simplifier.triangles[i].uvs[0].y);
        uvs[v2] = glm::vec2(simplifier.triangles[i].uvs[1].x, simplifier.triangles[i].uvs[1].y);
        uvs[v3] = glm::vec2(simplifier.triangles[i].uvs[2].x, simplifier.triangles[i].uvs[2].y);
    }

    for(int i = 0; i < vertSize; i++){
        float len = glm::length(vertex_normals[i]);
        if(len > 0.0){
            vertex_normals[i] /= len;
        } else {
            vertex_normals[i] = glm::vec3(0.0, 1.0, 0.0);
        }
    }

    std::vector<int> vertexIndexMap;
    for(int i = 0; i < vertSize; i++){
        if(!usedVertices[i]) vertexIndexMap.push_back(-1);
        else {
            unsigned int vertex = createVertex(vertices[i], vertex_normals[i], glm::vec2(0));
			vertexIndexMap.push_back(vertex);
        }
    }

    vertices.clear();
    vertex_normals.clear();
    uvs.clear();
    usedVertices.clear();

	int prevTriSize = m_faces.size();

    for(int i = 0; i < triSize; i++){
        if(simplifier.triangles[i].deleted) continue;
        int v1 = simplifier.triangles[i].v[0];
        int v2 = simplifier.triangles[i].v[1];
        int v3 = simplifier.triangles[i].v[2];

        int vertex1 = vertexIndexMap[v1];
        int vertex2 = vertexIndexMap[v2];
        int vertex3 = vertexIndexMap[v3];

        assert(vertex1 != -1 && vertex2 != -1 && vertex3 != -1);

        SimpleFace* face = new SimpleFace(vertex1, vertex2, vertex3, i + prevTriSize);
        m_faces.push_back(face);

		createEdge(vertex1, vertex2, face);
		createEdge(vertex2, vertex3, face);
		createEdge(vertex3, vertex1, face);
    }
}

void LodMesh::printLODInformation() {
    for (int i = 0; i < lodMesh.size(); i++) {
        SimpleMesh* mesh = lodMesh[i];
        if(mesh == nullptr) break;
        std::cout << "LOD " << i << ":\n";
        std::cout << "Vertices: " << mesh->m_vertices.size() << std::endl;
        std::cout << "Faces: " << mesh->m_faces.size() << std::endl;
        std::cout << "Clusters: " << mesh->m_clusters.size() << std::endl;
        std::cout << "ClusterGroups: " << mesh->m_clusterGroups.size() << std::endl;
        std::cout << "ClusterGroupSize & Error: ";
        for(int j = 0; j < mesh->m_clusterGroups.size(); j++){
            int triCount = 0;
            for(auto clusters : mesh->m_clusterGroups[j]->m_clusterlist){
                triCount = clusters->endIdx - clusters->startIdx + 1;
            }
            std::cout << "[" << j << "]: (" << triCount << " : "<< mesh->m_clusterGroups[j]->error << ") ";
        }
        std::cout << std::endl;
    }
}
