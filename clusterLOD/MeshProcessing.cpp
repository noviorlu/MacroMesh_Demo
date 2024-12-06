#include <metis.h>
#include <filesystem>

#include <queue>
#include <iostream>
#include <algorithm>
#include <map>
#include <fstream>
#include <string>
#include <sstream>
#include <map>
#include <cstdio>
#include <cstdlib>
#ifdef _WIN32
#define popen _popen
#define pclose _pclose
#endif

#include "HalfEdgeMesh.hpp"

#define MAX_TRI_IN_CLUSTER 256
#define MAX_CLUSTER_IN_CLUSTERGROUP 32
#define MAX_TRI_IN_CLUSTERGROUP (MAX_TRI_IN_CLUSTER * MAX_CLUSTER_IN_CLUSTERGROUP)

void HalfEdgeMesh::BuildAdjacencyListForCluster(
    std::vector<std::vector<size_t>>& adjacency_list
) {
    adjacency_list.clear();
	adjacency_list.resize(m_clusterOffsets.size());
    std::vector<std::unordered_set<size_t>> adjacency_set(m_clusterOffsets.size());

	std::unordered_map<const Face*, size_t> faceToIndexMap;
    for (size_t i = 0; i < m_faces.size(); ++i) {
        faceToIndexMap[m_faces[i]] = i;
    }

    auto findClusterIndex = [this, &faceToIndexMap](const Face* face) -> size_t {
        size_t faceIdx = faceToIndexMap[face];
        auto clusterIt = std::lower_bound(m_clusterOffsets.begin(), m_clusterOffsets.end(), faceIdx);
        return (clusterIt == m_clusterOffsets.end() || *clusterIt > faceIdx) ? (clusterIt - m_clusterOffsets.begin() - 1) : (clusterIt - m_clusterOffsets.begin());
    };

    for(int i = 0; i < m_edges.size(); ++i){
        HalfEdge* edge = m_edges[i];
        if (!edge->isValid) continue;
		if (edge > edge->twin) continue;

        Face* face1 = edge->face;
        Face* face2 = edge->twin->face;

        if (face1 && face2) {
            size_t clusterIdx1 = findClusterIndex(face1);
            size_t clusterIdx2 = findClusterIndex(face2);

            if (clusterIdx1 != clusterIdx2) {
                adjacency_set[clusterIdx1].insert(clusterIdx2);
                adjacency_set[clusterIdx2].insert(clusterIdx1);
            }
        }
    }

    for (size_t i = 0; i < m_clusterOffsets.size(); ++i) {
        adjacency_list[i].reserve(adjacency_set[i].size());
        adjacency_list[i] = std::vector<size_t>(adjacency_set[i].begin(), adjacency_set[i].end());
    }
}


void HalfEdgeMesh::BuildAdjacencyListForRange(
    std::vector<std::vector<size_t>>& adjacency_list,
    size_t startIdx,
    size_t endIdx
) {
    size_t rangeSize = endIdx - startIdx + 1;

    std::unordered_map<const Face*, size_t> faceToIndexMap;
    faceToIndexMap.reserve(rangeSize);
    for (size_t i = 0; i < rangeSize; ++i) {
        faceToIndexMap[m_faces[i]] = i;
    }

    adjacency_list.resize(rangeSize);
    for (size_t i = startIdx; i <= endIdx; ++i) {
        const Face* face = m_faces[i];
        const HalfEdge* edge = face->edge;
        do {
            const HalfEdge* twin = edge->twin;
            if (twin && twin->face) {
                auto it = faceToIndexMap.find(twin->face);
                if (it != faceToIndexMap.end()) {
                    adjacency_list[i - startIdx].push_back(it->second);
                }
            }
            edge = edge->next;
        } while (edge != face->edge);
    }
}

    void HalfEdgeMesh::HalfEdgeMeshSplitterRecursive(
        size_t startIdx,
        size_t endIdx,
        bool isParentClusterGroup = false,
        int depth = 0
    ) {
        if(depth > 3){
            #pragma omp critical
            m_clusterGroupOffset.push_back(startIdx);
            return;
        }
        
        // if (startIdx >= endIdx) { return; }
        idx_t numElements = static_cast<idx_t>(endIdx - startIdx + 1);
        // if (numElements <= MAX_TRI_IN_CLUSTER) {
        //     m_clusterOffsets.push_back(startIdx);
        //     return;
        // }

        // if(numElements <= MAX_TRI_IN_CLUSTERGROUP && !isParentClusterGroup){
        //     m_clusterGroupOffset.push_back(startIdx);
        //     isParentClusterGroup = true;
        //     return;
        // }

        std::vector<std::vector<size_t>> adjacencyList;
        BuildAdjacencyListForRange(adjacencyList, startIdx, endIdx);

        idx_t numConstraints = 1;
        std::vector<idx_t> xadj(numElements + 1, 0);
        std::vector<idx_t> adjncy;

        for (size_t i = 0; i < numElements; ++i) {
            xadj[i + 1] = xadj[i] + adjacencyList[i].size();
            for (size_t neighbor : adjacencyList[i]) {
                adjncy.push_back(static_cast<idx_t>(neighbor));
            }
        }

        idx_t numParts = 2;
        std::vector<idx_t> partitionResult(numElements);
        idx_t edgeCut;

        int result = METIS_PartGraphRecursive(
            &numElements,           
            &numConstraints,        
            xadj.data(),            
            adjncy.data(),          
            NULL, NULL, NULL,       
            &numParts,              
            NULL, NULL, NULL,      
            &edgeCut,               
            partitionResult.data()  
        );

        if (result != METIS_OK) {
            std::cerr << "Partitioning failed with error code: " << result << std::endl;
            return;
        }

        size_t part0ptr = startIdx;
        size_t part1ptr = endIdx;

        while (part0ptr <= part1ptr) {
            if (partitionResult[part0ptr - startIdx] == 0) {
                ++part0ptr;
            } else {
                while (part1ptr > startIdx && partitionResult[part1ptr - startIdx] == 1) {
                    --part1ptr;
                }
                if (part0ptr < part1ptr) {
                    std::swap(m_faces[part0ptr], m_faces[part1ptr]);
                    ++part0ptr;
                    --part1ptr;
                }
            }
        }

        size_t midIdx = part0ptr;
        #pragma omp parallel sections
        {
            #pragma omp section
            {
                HalfEdgeMeshSplitterRecursive(startIdx, midIdx - 1, isParentClusterGroup, depth+1);
            }
            #pragma omp section
            {
                HalfEdgeMeshSplitterRecursive(midIdx, endIdx, isParentClusterGroup, depth+1);
            }
        }
    }

    void HalfEdgeMesh::HalfEdgeMeshSplitter() {
        m_clusterOffsets.clear();
        auto start = std::chrono::high_resolution_clock::now();
        m_clusterGroupOffset.clear();
        HalfEdgeMeshSplitterRecursive(0, m_faces.size() - 1, false, 0);
        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "Mesh Partitioning Time: " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << "s" << std::endl;
        std::cout << "CLusterGroup Offset: ";
        
        std::sort(m_clusterGroupOffset.begin(), m_clusterGroupOffset.end());

        m_clusterGroupOffset.push_back(m_faces.size());
        for(auto val : m_clusterGroupOffset) std::cout << val << " ";
        std::cout << std::endl;
        // start = std::chrono::high_resolution_clock::now();
        // std::vector<std::vector<size_t>> adjacencyList;
        // BuildAdjacencyListForCluster(adjacencyList);
        // ClusterGrouper(adjacencyList, m_clusterGroupResult, m_clusterGroupCount);
        // end = std::chrono::high_resolution_clock::now();
        // std::cout << "Cluster Grouping Time: " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << "s" << std::endl;
    }







void ClusterGrouper(
    const std::vector<std::vector<size_t>>& adjacencyList,
    std::vector<size_t>&  clusterGroupResult,
    int& clusterGroupCount
) 
{
    idx_t numClusters = adjacencyList.size();

    std::vector<idx_t> xadj(numClusters + 1, 0);
    std::vector<idx_t> adjncy;

    for (size_t i = 0; i < numClusters; ++i) {
        xadj[i + 1] = xadj[i] + adjacencyList[i].size();
        for (size_t neighbor : adjacencyList[i]) {
            adjncy.push_back(static_cast<idx_t>(neighbor));
        }
    }

    idx_t numParts = static_cast<idx_t>(
        (numClusters + MAX_CLUSTER_IN_CLUSTERGROUP - 1) / MAX_CLUSTER_IN_CLUSTERGROUP
    );
    if (numParts > numClusters) {
        numParts = numClusters;
    }

    idx_t ncon = 1;
    std::vector<real_t> tpwgts(numParts * ncon, 1.0 / numParts);
    std::vector<real_t> ubvec(ncon, 1.05);

    std::vector<idx_t> partitionResult(numClusters, 0);
    idx_t edgeCut;

    int result = METIS_PartGraphKway(
        &numClusters,
        &ncon,
        xadj.data(),
        adjncy.data(),
        NULL,
        NULL,
        NULL,
        &numParts,
        tpwgts.data(),
        ubvec.data(),
        NULL,
        &edgeCut,
        partitionResult.data()
    );

    if (result != METIS_OK) {
        std::cerr << "METIS partitioning failed with error code: " << result << std::endl;
        return;
    }

    clusterGroupResult.clear();
    clusterGroupResult.reserve(numClusters);
    for (size_t i = 0; i < numClusters; ++i) {
        clusterGroupResult.push_back(partitionResult[i]);
    }
    clusterGroupCount = numParts;
}

void FastQEM(const std::string& lodFolder) {
    std::vector<std::filesystem::path> objFiles;
    for (const auto& entry : std::filesystem::directory_iterator(lodFolder)) {
        if (entry.path().extension() == ".obj") {
            objFiles.push_back(entry.path());
        }
    }

    std::string simplifyExe = "C:/projects/MacroMesh_Demo/clusterLOD/QEM/simplify.exe";
    std::string tempOutputFolder = lodFolder + "/Temp";
    std::map<std::string, double> fileErrorMap;

    // 创建输出文件夹和临时文件夹
    std::filesystem::create_directories(tempOutputFolder);

    #pragma omp parallel for
    for (int i = 0; i < objFiles.size(); ++i) {
        const auto& objFile = objFiles[i];
        std::string inputPath = objFile.generic_string();

        // 将简化后的文件输出到新的 outputFolder
        std::string tempOutputPath = (std::filesystem::path(tempOutputFolder) / objFile.filename()).generic_string();
        std::string finalOutputPath = (std::filesystem::path(lodFolder) / objFile.filename()).generic_string();

        std::string command = simplifyExe + " " + inputPath + " " + tempOutputPath + " 0.5";

        std::string output;
        char buffer[128];
        FILE* pipe = popen(command.c_str(), "r");
        if (!pipe) {
            std::cerr << "Error: failed to run command." << std::endl;
            continue;
        }

        while (fgets(buffer, sizeof(buffer), pipe) != nullptr) {
            output += buffer;
        }
        pclose(pipe);

        std::istringstream stream(output);
        std::string line;
        double totalError = -1;
        while (std::getline(stream, line)) {
            if (line.find("Total error:") != std::string::npos) {
                std::istringstream lineStream(line);
                std::string label;
                lineStream >> label >> label >> totalError;
                break;
            }
        }

        #pragma omp critical
        {
            if (totalError >= 0) {
                fileErrorMap[objFile.filename().string()] = totalError;
            } else {
                std::cerr << "Warning: Total error not found for " << inputPath << std::endl;
            }
        }

        // try {
        //     // 将临时文件移动到最终输出文件夹中
        //     std::filesystem::rename(tempOutputPath, finalOutputPath);
        // } catch (const std::exception& e) {
        //     std::cerr << "Error moving file to output folder " << finalOutputPath << ": " << e.what() << std::endl;
        // }
    }

    // 删除临时文件夹
    // std::filesystem::remove_all(tempOutputFolder);

    // 输出生成的文件信息
    std::vector<HalfEdgeMesh> meshes;
    for (const auto& [filename, error] : fileErrorMap) {
        std::cout << "File: " << filename << ", Total Error: " << error << std::endl;
        // HalfEdgeMesh mesh;
        // mesh.importMesh((std::filesystem::path(outputFolder) / filename).generic_string(), error);
        // meshes.push_back(std::move(mesh));
    }
}















void EMesh::buildAdjacencyFaces(){
    int i = 0;
    for (auto& face : m_faces) {
        VertexPair pair1(face->vertices[0], face->vertices[1]);
        VertexPair pair2(face->vertices[1], face->vertices[2]);
        VertexPair pair3(face->vertices[2], face->vertices[0]);

        for(auto pair : {pair1, pair2, pair3}) {
            for(auto adjFace : m_edgeMap[pair]->faces){
                if(face == adjFace) continue;
                face->adjacentFaces.push_back(adjFace);
            }
        }
    }
}

void EMesh::buildAdjacencyListForRange(
    std::vector<std::vector<size_t>>& adjacency_list,
    size_t startIdx,
    size_t endIdx
) {
    std::unordered_map<const EFace*, size_t> faceToIndexMap;
    for (size_t i = startIdx; i <= endIdx; ++i) {
        faceToIndexMap[m_faces[i]] = i - startIdx;
    }

    size_t rangeSize = endIdx - startIdx + 1;
    adjacency_list.resize(rangeSize);
    for (size_t i = startIdx; i <= endIdx; ++i) {
        const EFace* face = m_faces[i];
        // use the precomputed adjacentFaces in the face
        for (const EFace* adjFace : face->adjacentFaces) {
            auto it = faceToIndexMap.find(adjFace);
            if (it != faceToIndexMap.end()) {
                adjacency_list[i - startIdx].push_back(it->second);
            }
        }
    }
}

void EMesh::eMeshSplitterRecursive(
    size_t startIdx,
    size_t endIdx
) {
    if (startIdx >= endIdx) { return; }

    idx_t numElements = static_cast<idx_t>(endIdx - startIdx + 1);
    if (numElements <= MAX_TRI_IN_CLUSTER) {
        m_clusterOffsets.push_back(startIdx);
        return;
    }

    std::vector<std::vector<size_t>> adjacencyList;
    buildAdjacencyListForRange(adjacencyList, startIdx, endIdx);

    idx_t numConstraints = 1;
    std::vector<idx_t> xadj(numElements + 1, 0);
    std::vector<idx_t> adjncy;

    for (size_t i = 0; i < numElements; ++i) {
        xadj[i + 1] = xadj[i] + adjacencyList[i].size();
        for (size_t neighbor : adjacencyList[i]) {
            adjncy.push_back(static_cast<idx_t>(neighbor));
        }
    }

    idx_t numParts = 2;
    std::vector<idx_t> partitionResult(numElements);
    idx_t edgeCut;

    int result = METIS_PartGraphRecursive(
        &numElements,
        &numConstraints,
        xadj.data(),
        adjncy.data(),
        NULL, NULL, NULL,
        &numParts,
        NULL, NULL, NULL,
        &edgeCut,
        partitionResult.data()
    );

    if (result != METIS_OK) {
        std::cerr << "Partitioning failed with error code: " << result << std::endl;
        return;
    }

    size_t part0ptr = startIdx;
    size_t part1ptr = endIdx;

    while (part0ptr <= part1ptr) {
        if (partitionResult[part0ptr - startIdx] == 0) {
            ++part0ptr;
        } else {
            while (part1ptr > startIdx && partitionResult[part1ptr - startIdx] == 1) {
                --part1ptr;
            }
            if (part0ptr < part1ptr) {
                std::swap(m_faces[part0ptr], m_faces[part1ptr]);
                ++part0ptr;
                --part1ptr;
            }
        }
    }

    size_t midIdx = part0ptr;

    eMeshSplitterRecursive(startIdx, midIdx - 1);
    eMeshSplitterRecursive(midIdx, endIdx);
}

void EMesh::eMeshSplitter() {
    m_clusterOffsets.clear();
    auto start = std::chrono::high_resolution_clock::now();
    
    buildAdjacencyFaces();
    eMeshSplitterRecursive(0, m_faces.size() - 1);
    
    
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Mesh Partitioning Time: " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << "s" << std::endl;

    start = std::chrono::high_resolution_clock::now();
    std::vector<std::vector<size_t>> adjacencyList;
    buildAdjacencyListForCluster(adjacencyList);
    ClusterGrouper(adjacencyList, m_clusterGroupResult, m_clusterGroupCount);
    end = std::chrono::high_resolution_clock::now();
    std::cout << "Cluster Grouping Time: " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << "s" << std::endl;

}

void EMesh::buildAdjacencyListForCluster(
    std::vector<std::vector<size_t>>& adjacency_list
) {
    std::unordered_map<const EFace*, size_t> faceToIndexMap;
    for (size_t i = 0; i < m_faces.size(); ++i) {
        faceToIndexMap[m_faces[i]] = i;
    }
    
    std::vector<std::unordered_set<size_t>> adjacencySetList(m_clusterOffsets.size());

    auto findClusterIndex = [this](int faceIdx) -> size_t {
        auto it = std::lower_bound(m_clusterOffsets.begin(), m_clusterOffsets.end(), faceIdx);
        return (it == m_clusterOffsets.end() || *it > faceIdx) ? (it - m_clusterOffsets.begin() - 1) : (it - m_clusterOffsets.begin());
    };

    for (size_t i = 0; i < m_clusterOffsets.size(); ++i) {
        size_t startIdx = m_clusterOffsets[i];
        size_t endIdx = (i + 1 < m_clusterOffsets.size()) ? m_clusterOffsets[i + 1] : m_faces.size();

        for (size_t j = startIdx; j < endIdx; ++j) {
            for (auto adjFace : m_faces[j]->adjacentFaces) {
                size_t adjClusterIdx = findClusterIndex(faceToIndexMap[adjFace]);
                if (i != adjClusterIdx) {
                    adjacencySetList[i].insert(adjClusterIdx);
                }
            }
        }
    }

    adjacency_list.resize(m_clusterOffsets.size());
    for (size_t i = 0; i < adjacencySetList.size(); ++i) {
        adjacency_list[i] = std::vector<size_t>(adjacencySetList[i].begin(), adjacencySetList[i].end());
    }
}

