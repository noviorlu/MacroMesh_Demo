#include <metis.h>

#include <queue>
#include <iostream>
#include <algorithm>
#include "HalfEdgeMesh.hpp"


#define MAX_TRI_IN_CLUSTER 256
#define MAX_CLUSTER_IN_CLUSTERGROUP 32
#define MAX_TRI_IN_CLUSTERGROUP (MAX_TRI_IN_CLUSTER * MAX_CLUSTER_IN_CLUSTERGROUP)

void HalfEdgeMesh::BuildAdjacencyListForRange(
    std::vector<std::vector<size_t>>& adjacency_list,
    size_t startIdx,
    size_t endIdx
) {
    size_t rangeSize = endIdx - startIdx + 1;
    adjacency_list.resize(rangeSize);

    for (size_t i = startIdx; i <= endIdx; ++i) {
        const HalfEdge* edge = m_faces[i].edge;

        do {
            const HalfEdge* twin = edge->twin;
            if (twin && twin->face) {
                size_t twinFaceIndex = static_cast<size_t>(twin->face - &m_faces[0]);

                if (twinFaceIndex >= startIdx && twinFaceIndex <= endIdx) {
                    adjacency_list[i - startIdx].push_back(twinFaceIndex - startIdx);
                }
            }
            edge = edge->next;
        } while (edge != m_faces[i].edge);
    }
}

int depth = 0;
int aimdepth = 2;

void HalfEdgeMesh::HalfEdgeMeshSplitterRecursive(
    size_t startIdx,
    size_t endIdx,
    bool isParentClusterGroup = false
) {
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
                // before swapping update the halfedge face pointers

                HalfEdge* edge0 = m_faces[part0ptr].edge;
                do {
                    edge0->face = &m_faces[part1ptr];
                    edge0 = edge0->next;
                } while (edge0 != m_faces[part0ptr].edge);

                HalfEdge* edge1 = m_faces[part1ptr].edge;
                do {
                    edge1->face = &m_faces[part0ptr];
                    edge1 = edge1->next;
                } while (edge1 != m_faces[part1ptr].edge);
                
                std::swap(m_faces[part0ptr], m_faces[part1ptr]);
                
                ++part0ptr;
                --part1ptr;
            }
        }
    }

    size_t midIdx = part0ptr;

    HalfEdgeMeshSplitterRecursive(startIdx, midIdx - 1, isParentClusterGroup);
    HalfEdgeMeshSplitterRecursive(midIdx, endIdx, isParentClusterGroup);
}

void HalfEdgeMesh::HalfEdgeMeshSplitter() {
    m_clusterOffsets.clear();
    m_clusterGroupOffsets.clear();

    HalfEdgeMeshSplitterRecursive(0, m_faces.size() - 1);

    // print the cluster and cluster group offsets
    std::cout << "Cluster Offsets: ";
    for (size_t offset : m_clusterOffsets) {
        std::cout << offset << " ";
    }

    std::cout << std::endl;

    std::cout << "Cluster Group Offsets: ";
    for (size_t offset : m_clusterGroupOffsets) {
        std::cout << offset << " ";
    }

    std::cout << std::endl;
}


void HalfEdgeMesh::partition_loop(){

    HalfEdgeMeshSplitter();
}