#include "Cluster.hpp"
#include <iostream>
#include <metis.h>
#include <tuple>
#include <unordered_map>
#include <set>


void MeshSplitter(Mesh& mesh, std::vector<Cluster>& clusters, int num_parts) {
    // 创建邻接列表
    auto adjacency_list = BuildAdjacencyList(mesh.m_indexData);
    std::vector<int> partition_result; // 存储分区结果

    // 计算三角形数量
    idx_t num_elements = mesh.m_indexData.size() / 3;
    idx_t num_constraints = 1;
    idx_t objval;

    // 初始化 METIS 的 xadj 和 adjncy 数组
    std::vector<idx_t> xadj(num_elements + 1, 0);
    std::vector<idx_t> adjncy;

    // 构建 xadj 和 adjncy 数组
    for (size_t i = 0; i < num_elements; ++i) {
        xadj[i + 1] = xadj[i] + adjacency_list[i].size();
        for (size_t neighbor : adjacency_list[i]) {
            adjncy.push_back(neighbor);
        }
    }

    // 初始化分区结果数组
    partition_result.resize(num_elements);

    // 调用 METIS_PartGraphKway 进行分区
    int result = METIS_PartGraphKway(&num_elements,       // 节点数量（三角形数量）
                                     &num_constraints,    // 约束数量
                                     xadj.data(),         // 邻接关系的起始索引
                                     adjncy.data(),       // 邻接节点
                                     NULL, NULL, NULL,    // 权重参数（可设为 NULL）
                                     &num_parts,          // 目标簇数量
                                     NULL, NULL,          // 选项（可设为 NULL）
                                     &objval,             // 输出目标值
                                     partition_result.data()); // 输出的分区结果

    if (result != METIS_OK) {
        std::cerr << "Partitioning failed with error code: " << result << std::endl;
        return;
    }
    
    std::cout << "Partitioning successful!" << std::endl;

    // 将分区结果组织到 clusters 中
    clusters.resize(num_parts);
    for (size_t i = 0; i < partition_result.size(); ++i) {
        int cluster_id = partition_result[i];
        clusters[cluster_id].triangles.push_back(mesh.m_indexData[i * 3]);
        clusters[cluster_id].triangles.push_back(mesh.m_indexData[i * 3 + 1]);
        clusters[cluster_id].triangles.push_back(mesh.m_indexData[i * 3 + 2]);
    }
}
std::vector<std::vector<size_t>> 
BuildAdjacencyList(const std::vector<unsigned int>& m_indexData) 
{
    std::vector<std::vector<size_t>> adjacency_list;
    std::unordered_map<std::tuple<size_t, size_t>, std::vector<size_t>> edge_to_triangles;

    // 每三个索引表示一个三角形
    size_t num_triangles = m_indexData.size() / 3;
    
    // 辅助函数，用于生成有序的边
    auto make_edge = [](size_t a, size_t b) {
        return std::make_tuple(std::min(a, b), std::max(a, b));
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