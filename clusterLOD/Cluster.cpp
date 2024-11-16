#include "Cluster.hpp"
#include <iostream>
#include <metis.h>
#include <unordered_map>
#include <set>
#include <algorithm>

void MeshSplitter(Mesh& mesh, std::vector<Cluster>& clusters, int num_parts) {
    // 创建邻接列表
    auto adjacency_list = BuildAdjacencyList(mesh.m_indexData);
    idx_t num_elements = mesh.m_indexData.size() / 3; // 节点数量（三角形数量）
    idx_t num_constraints = 1;                        // 通常为 1

    std::vector<idx_t> xadj(num_elements + 1, 0);     // 邻接起始索引数组
    std::vector<idx_t> adjncy;                        // 邻接节点数组

    std::vector<idx_t> partition_result(num_elements); // 存储分区结果
    idx_t objval;                                      // 分区目标值

    // 构建 xadj 和 adjncy 数组
    for (size_t i = 0; i < num_elements; ++i) {
        xadj[i + 1] = xadj[i] + adjacency_list[i].size();
        for (size_t neighbor : adjacency_list[i]) {
            adjncy.push_back(neighbor);
        }
    }

    // 初始化分区结果数组
    partition_result.resize(num_elements);

    // 调用 METIS_PartGraphKway
    int result = METIS_PartGraphKway(
        &num_elements,             // 节点数量
        &num_constraints,          // 约束数量
        xadj.data(),               // 邻接关系的起始索引
        adjncy.data(),             // 邻接节点
        NULL,                      // 顶点权重（可为 NULL）
        NULL,                      // 节点大小（可为 NULL）
        NULL,                      // 边权重（可为 NULL）
        &num_parts,                // 目标簇数量
        NULL,                      // 每个分区的目标权重（可为 NULL）
        NULL,                      // 不平衡向量（可为 NULL）
        NULL,                      // 选项（可为 NULL）
        &objval,                   // 输出目标值
        partition_result.data()    // 输出的分区结果
    );

    if (result == METIS_OK) {
        std::cout << "Partitioning successful!" << std::endl;
    } else {
        std::cerr << "Partitioning failed with error code: " << result << std::endl;
    }

    // group partitioning
    std::unordered_map<idx_t, std::vector<unsigned int>> partition_map;
    for (size_t i = 0; i < num_elements; ++i) {
        partition_map[partition_result[i]].push_back(i);
    }

    // const std::vector<unsigned int>& triList = partition_map[0];
    // // Remove all other triangles in m_indexData except the triList
    
    // std::vector<unsigned int> new_indexData;
    // for (unsigned int tri : triList) {
    //     new_indexData.push_back(mesh.m_indexData[3 * tri]);
    //     new_indexData.push_back(mesh.m_indexData[3 * tri + 1]);
    //     new_indexData.push_back(mesh.m_indexData[3 * tri + 2]);
    // }
    // mesh.m_indexData = std::move(new_indexData);

    // // print the triangle indices data
    // for (size_t i = 0; i < mesh.m_indexData.size(); i += 3) {
    //     std::cout << mesh.m_indexData[i] << " " 
    //               << mesh.m_indexData[i + 1] << " " 
    //               << mesh.m_indexData[i + 2] << std::endl;
    // }


    // for(auto& entry : partition_map) {
    //     Mesh* new_mesh = new Mesh(mesh, entry.second);
    //     clusters.push_back(Cluster(0, new_mesh));
    // }
}

struct pair_hash {
    template <class T1, class T2>
    std::size_t operator()(const std::pair<T1, T2>& p) const {
        auto hash1 = std::hash<T1>{}(p.first);
        auto hash2 = std::hash<T2>{}(p.second);
        return hash1 ^ (hash2 << 1); // Combine the two hashes
    }
};

std::vector<std::vector<size_t>> 
BuildAdjacencyList(const std::vector<unsigned int>& m_indexData) 
{
    std::vector<std::vector<size_t>> adjacency_list;
    std::unordered_map<std::pair<size_t, size_t>, std::vector<size_t>, pair_hash> edge_to_triangles;

    // 每三个索引表示一个三角形
    size_t num_triangles = m_indexData.size() / 3;
    
    // 辅助函数，用于生成有序的边
    auto make_edge = [](size_t a, size_t b) {
        return std::make_pair(std::min(a, b), std::max(a, b));
    };

    // 构建边到三角形的映射
    for (size_t i = 0; i < num_triangles; ++i) {
        size_t v1 = m_indexData[3 * i];
        size_t v2 = m_indexData[3 * i + 1];
        size_t v3 = m_indexData[3 * i + 2];

        // 创建三条边并记录它们所属的三角形
        auto edge1 = make_edge(v1, v2);
        auto edge2 = make_edge(v2, v3);
        auto edge3 = make_edge(v3, v1);

        edge_to_triangles[edge1].push_back(i);
        edge_to_triangles[edge2].push_back(i);
        edge_to_triangles[edge3].push_back(i);
    }

    // 初始化邻接列表
    adjacency_list.resize(num_triangles);

    // 使用边信息构建邻接列表
    for (const auto& entry : edge_to_triangles) {
        const std::vector<size_t>& tris = entry.second;

        if (tris.size() > 1) {
            for (size_t j = 0; j < tris.size(); ++j) {
                for (size_t k = j + 1; k < tris.size(); ++k) {
                    size_t tri1 = tris[j];
                    size_t tri2 = tris[k];
                    
                    // 双向添加邻接关系
                    adjacency_list[tri1].push_back(tri2);
                    adjacency_list[tri2].push_back(tri1);
                }
            }
        }
    }

    // 去重，移除重复的邻接项
    for (auto& neighbors : adjacency_list) {
        std::set<size_t> unique_neighbors(neighbors.begin(), neighbors.end());
        neighbors.assign(unique_neighbors.begin(), unique_neighbors.end());
    }

    return adjacency_list;
}