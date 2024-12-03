#include <metis.h>

#include <queue>
#include <iostream>
#include <algorithm>
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
    std::unordered_map<const Face*, size_t> faceToIndexMap;
    for (size_t i = startIdx; i <= endIdx; ++i) {
        faceToIndexMap[m_faces[i]] = i - startIdx;
    }

    size_t rangeSize = endIdx - startIdx + 1;
    adjacency_list.resize(rangeSize);

    for (size_t i = startIdx; i <= endIdx; ++i) {
        const HalfEdge* edge = m_faces[i]->edge;

        do {
            const HalfEdge* twin = edge->twin;
            if (twin && twin->face) {
                auto it = faceToIndexMap.find(twin->face);
                if (it != faceToIndexMap.end()) {
                    adjacency_list[i - startIdx].push_back(it->second);
                }
            }
            edge = edge->next;
        } while (edge != m_faces[i]->edge);
    }
}

void HalfEdgeMesh::HalfEdgeMeshSplitterRecursive(
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
            HalfEdgeMeshSplitterRecursive(startIdx, midIdx - 1);
        }
        #pragma omp section
        {
            HalfEdgeMeshSplitterRecursive(midIdx, endIdx);
        }
    }
}

void HalfEdgeMesh::HalfEdgeMeshSplitter() {
    m_clusterOffsets.clear();
    auto start = std::chrono::high_resolution_clock::now();
    HalfEdgeMeshSplitterRecursive(0, m_faces.size() - 1);
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Mesh Partitioning Time: " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << "s" << std::endl;

    start = std::chrono::high_resolution_clock::now();
    std::vector<std::vector<size_t>> adjacencyList;
    BuildAdjacencyListForCluster(adjacencyList);
    ClusterGrouper(adjacencyList, m_clusterGroupResult, m_clusterGroupCount);
    end = std::chrono::high_resolution_clock::now();
    std::cout << "Cluster Grouping Time: " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << "s" << std::endl;
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
