#include <metis.h>

#include <queue>
#include <iostream>
#include <algorithm>
#include "HalfEdgeMesh.hpp"

// #define DEBUG_MESHPROCESSING

#define MAX_TRI_IN_CLUSTER 256
#define MAX_CLUSTER_IN_CLUSTERGROUP 32
#define MAX_TRI_IN_CLUSTERGROUP (MAX_TRI_IN_CLUSTER * MAX_CLUSTER_IN_CLUSTERGROUP)

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

#ifdef DEBUG_MESHPROCESSING
int depth = 0;
int aimdepth = 3;
#endif

void HalfEdgeMesh::HalfEdgeMeshSplitterRecursive(
    size_t startIdx,
    size_t endIdx,
    bool isParentClusterGroup = false
) {

#ifdef DEBUG_MESHPROCESSING
    ++depth;
    if(depth == aimdepth){
        std::cout << "startIdx: " << startIdx << " endIdx: " << endIdx << std::endl;
        m_clusterOffsets.push_back(startIdx);
        return;
    }
#endif

    if (startIdx >= endIdx) { return; }
    
    idx_t numElements = static_cast<idx_t>(endIdx - startIdx + 1);
    if (numElements <= MAX_TRI_IN_CLUSTER) {
        m_clusterOffsets.push_back(startIdx);
        return;
    }

    if (!isParentClusterGroup && numElements <= MAX_TRI_IN_CLUSTERGROUP) {
        m_clusterGroupOffsets.push_back(startIdx);
        isParentClusterGroup = true;
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

    HalfEdgeMeshSplitterRecursive(startIdx, midIdx - 1, isParentClusterGroup);
#ifdef DEBUG_MESHPROCESSING
    depth--;
#endif
    HalfEdgeMeshSplitterRecursive(midIdx, endIdx, isParentClusterGroup);
#ifdef DEBUG_MESHPROCESSING
    depth--;
#endif
}

void HalfEdgeMesh::HalfEdgeMeshSplitter() {
    m_clusterOffsets.clear();
    m_clusterGroupOffsets.clear();

    HalfEdgeMeshSplitterRecursive(0, m_faces.size() - 1);

#ifdef DEBUG_MESHPROCESSING
    //print the cluster and cluster group offsets
    std::cout << "Cluster Offsets: ";
    for (size_t i = 0; i < m_clusterOffsets.size(); ++i) {
        std::cout << m_clusterOffsets[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "Cluster Group Offsets: ";
    for (size_t i = 0; i < m_clusterGroupOffsets.size(); ++i) {
        std::cout << m_clusterGroupOffsets[i] << " ";
    }
    std::cout << std::endl;
#endif
}


void HalfEdgeMesh::partition_loop(){

    HalfEdgeMeshSplitter();
}