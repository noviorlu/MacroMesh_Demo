#include <metis.h>

#include <queue>
#include <iostream>
#include "HalfEdgeMesh.hpp"


#define MAX_TRI_IN_CLUSTER 256
#define MAXN_CLUSTER_IN_CLUSTERGROUP 32


void HalfEdgeMesh::BuildAdjacencyListFromHalfEdgeMesh(std::vector<std::vector<size_t>>& adjacency_list) {
    adjacency_list.resize(this->faces.size());

    for (size_t i = 0; i < this->faces.size(); ++i) {
        const Face& face = *(this->faces[i]);
        
        const HalfEdge* edge = face.edge;
        do {
            if (edge->twin) {
                size_t adjacent_face_index = edge->twin->face->index;
                adjacency_list[i].push_back(adjacent_face_index);
            }
            edge = edge->next;
        } while (edge != face.edge);
    }
}

void HalfEdgeMesh::HalfEdgeMeshSplitter() {
    std::vector<std::vector<size_t>> adjacency_list;
    BuildAdjacencyListFromHalfEdgeMesh(adjacency_list);
    idx_t num_elements = this->faces.size();
    idx_t num_constraints = 1;

    std::vector<idx_t> xadj(num_elements + 1, 0);
    std::vector<idx_t> adjncy;
    std::vector<idx_t> partition_result(num_elements);
    idx_t objval;

    for (size_t i = 0; i < num_elements; ++i) {
        xadj[i + 1] = xadj[i] + adjacency_list[i].size();
        for (size_t neighbor : adjacency_list[i]) {
            adjncy.push_back(neighbor);
        }
    }

    std::unordered_map<idx_t, std::vector<size_t>> partition_map;

    // METIS Splitter
    {
        int num_parts = static_cast<int>(std::ceil(static_cast<float>(num_elements) / MAX_TRI_IN_CLUSTER));
        int result = METIS_PartGraphKway(
            &num_elements,
            &num_constraints,
            xadj.data(),
            adjncy.data(),
            NULL, NULL, NULL,
            &num_parts,
            NULL, NULL, NULL,
            &objval,
            partition_result.data()
        );

        if (result != METIS_OK) {
            std::cerr << "Partitioning failed with error code: " << result << std::endl;
            return;
        }

        for (size_t i = 0; i < num_elements; ++i) {
            this->faces[i]->clusterIndex = partition_result[i];
        }
    }
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
    idx_t num_groups = static_cast<idx_t>(std::ceil(static_cast<float>(num_clusters) / MAXN_CLUSTER_IN_CLUSTERGROUP));
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