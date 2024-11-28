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
    size_t endIdx) 
{
    // Adjust the size of the adjacency list to match the range
    size_t rangeSize = endIdx - startIdx;
    adjacency_list.resize(rangeSize);

    // Iterate over the specified range of faces
    for (size_t i = startIdx; i < endIdx; ++i) {
        const Face& face = faces[i];
        const HalfEdge* edge = face.edge;

        do {
            if (edge->twin && edge->twin->face) {
                // Compute the index of the twin's face
                ptrdiff_t twinFaceIndex = edge->twin->face - faces.data();
                
                // Only consider neighbors within the specified range
                if (twinFaceIndex >= static_cast<ptrdiff_t>(startIdx) &&
                    twinFaceIndex < static_cast<ptrdiff_t>(endIdx)) 
                {
                    size_t localIndex = twinFaceIndex - startIdx;
                    adjacency_list[i - startIdx].push_back(localIndex);
                }
            }
            edge = edge->next;
        } while (edge != face.edge);
    }
}

void HalfEdgeMesh::HalfEdgeMeshSplitterRecursive(
    size_t startIdx,
    size_t endIdx,
    bool isParentClusterGroup = false
) {
    idx_t numElements = static_cast<idx_t>(endIdx - startIdx);
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

    // Use std::partition to group faces by partitionResult
    auto partitionPoint = std::partition(
        faces.begin() + startIdx, faces.begin() + endIdx,
        [&](const Face& face) {
            size_t idx = &face - &faces[0];
            return partitionResult[idx - startIdx] == 0;
        }
    );

    // Calculate split indices
    size_t midIdx = std::distance(faces.begin(), partitionPoint);

    HalfEdgeMeshSplitterRecursive(startIdx, midIdx, isParentClusterGroup);
    HalfEdgeMeshSplitterRecursive(midIdx, endIdx, isParentClusterGroup);
}


void HalfEdgeMesh::HalfEdgeMeshSplitter() {
    m_clusterOffsets.clear();
    m_clusterGroupOffsets.clear();

    HalfEdgeMeshSplitterRecursive(0, faces.size());

    // print the splitted cluster and cluster group
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