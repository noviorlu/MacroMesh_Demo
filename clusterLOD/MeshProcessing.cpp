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


void HalfEdgeMesh::BuildClusterAdjacency(std::unordered_map<int, std::unordered_set<int>>& cluster_adjacency) {
    for (const HalfEdge* edge : this->edges) {
        const Face* face1 = edge->face;
        const Face* face2 = edge->twin ? edge->twin->face : nullptr;

        if (face1 && face2) {
            int cluster1 = face1->clusterIndex;
            int cluster2 = face2->clusterIndex;

            if (cluster1 != cluster2) {
                cluster_adjacency[cluster1].insert(cluster2);
                cluster_adjacency[cluster2].insert(cluster1);
            }
        }
    }
}


void HalfEdgeMesh::GroupClustersWithMETIS() {
    std::unordered_map<int, std::unordered_set<int>> cluster_adjacency;
    BuildClusterAdjacency(cluster_adjacency);
    
    idx_t num_clusters = cluster_adjacency.size(); 
    if (num_clusters == 0) {
        std::cerr << "No clusters to group!" << std::endl;
        return;
    }

    std::vector<idx_t> xadj(num_clusters + 1, 0);
    std::vector<idx_t> adjncy;
    std::vector<idx_t> adjwgt;

    std::unordered_map<int, idx_t> cluster_to_idx;
    std::vector<int> idx_to_cluster(num_clusters);
    idx_t num_constraints = 1;
    idx_t num_groups = static_cast<idx_t>(std::ceil(static_cast<float>(num_clusters) / MAX_CLUSTER_IN_CLUSTERGROUP));
    std::vector<idx_t> partition_result(num_clusters, 0);
    idx_t objval;
    
    idx_t idx = 0;
    for (const auto& [cluster_id, neighbors] : cluster_adjacency) {
        cluster_to_idx[cluster_id] = idx;
        idx_to_cluster[idx] = cluster_id;
        idx++;
    }

    idx = 0;
    for (const auto& [cluster_id, neighbors] : cluster_adjacency) {
        xadj[idx + 1] = xadj[idx] + neighbors.size(); // 更新 xadj

        for (int neighbor : neighbors) {
            auto it = cluster_to_idx.find(neighbor);
            if (it != cluster_to_idx.end()) {
                adjncy.push_back(it->second);
                adjwgt.push_back(1);
            } else {
                throw std::runtime_error("Neighbor not found in cluster adjacency map!");
            }
        }

        idx++;
    }


    // 调用 METIS_PartGraphKway
    int result = METIS_PartGraphKway(
        &num_clusters,         // 节点数量
        &num_constraints,      // 约束数量
        xadj.data(),           // 邻接起始索引
        adjncy.data(),         // 邻接节点
        nullptr,               // 顶点权重（可为 NULL）
        nullptr,               // 节点大小（可为 NULL）
        adjwgt.data(),         // 边权重
        &num_groups,           // 分区数量
        nullptr,               // 每个分区的目标权重（可为 NULL）
        nullptr,               // 不平衡向量（可为 NULL）
        nullptr,               // 选项（可为 NULL）
        &objval,               // 输出目标值
        partition_result.data() // 分区结果
    );

    // 检查 METIS 调用结果
    if (result == METIS_OK) {
        std::cout << "Partitioning succeeded. Objective value: " << objval << std::endl;
    } else {
        std::cerr << "METIS failed with error code: " << result << std::endl;
    }

    for (idx_t i = 0; i < num_clusters; ++i) {
        clusterToGroupMap[idx_to_cluster[i]] = partition_result[i];
    }
}


void HalfEdgeMesh::partition_loop(){

    HalfEdgeMeshSplitter();

    GroupClustersWithMETIS();

}