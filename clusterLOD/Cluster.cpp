#include "Cluster.hpp"
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <algorithm>

glm::vec3 HSVtoRGB(float h, float s, float v) {
    h = fmod(h, 1.0f) * 6.0f;  // 将 h 限制在 [0, 6)
    int i = static_cast<int>(floor(h));
    float f = h - i; // 小数部分
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
        default: return glm::vec3(0.0f, 0.0f, 0.0f); // 容错
    }
}

glm::vec3 genRdColor(){
    // generate a random hue (h) in [0, 1)
    float h = static_cast<float>(rand()) / RAND_MAX; // 随机色相
    float s = 0.8f; // 固定饱和度，接近鲜艳的颜色
    float v = 0.8f; // 固定亮度

    // convert HSV to RGB
    return HSVtoRGB(h, s, v);
}

Cluster::Cluster(float Error, const Mesh& ref, const std::vector<unsigned int>& triIndices) : Error(Error), Mesh(ref, triIndices) {
    rdColor = genRdColor();
}

void Cluster::draw(const ShaderProgram& shader) const {
    CHECK_GL_ERRORS;
    shader.SetUniform3fv("material.kd", rdColor);
    CHECK_GL_ERRORS;
    Mesh::draw(shader);
}
  
void MeshSplitter(Mesh& mesh) {
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

    float counting_parts = float(mesh.m_indexData.size() / 3.0f / MAX_TRI_IN_CLUSTER);
    int num_parts = static_cast<int>(std::ceil(static_cast<float>(mesh.m_indexData.size()) / 3.0f / MAX_TRI_IN_CLUSTER));

    // output information about howmany cluster it's going to generate
    std::cout << "Generating " << num_parts << " clusters" << std::endl;

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

    // group partitioning REGION
    std::unordered_map<idx_t, std::vector<unsigned int>> partition_map;
    for (size_t i = 0; i < num_elements; ++i) {
        partition_map[partition_result[i]].push_back(i);
    }

    { /// Cluster Creation REGION
        // Step 1: 构建 Cluster Adjacency 和 ID 映射
        std::unordered_map<idx_t, std::unordered_set<idx_t>> cluster_adjacency;
        for (size_t tri = 0; tri < partition_result.size(); ++tri) {
            idx_t current_cluster = partition_result[tri];

            for (size_t neighbor : adjacency_list[tri]) {
                idx_t neighbor_cluster = partition_result[neighbor];
                if (current_cluster != neighbor_cluster) {
                    cluster_adjacency[current_cluster].insert(neighbor_cluster);
                }
            }
        }

        // Step 2: 创建 Clusters 并存入 mesh.m_clusterList
        std::unordered_map<idx_t, Cluster*> cluster_map;
        for (const auto& [cluster_id, triangles] : partition_map) {
            auto* cluster = new Cluster(0, mesh, triangles);
            cluster->adjacent_clusters.reserve(cluster_adjacency[cluster_id].size());
            mesh.m_clusterList.push_back(cluster);
            cluster_map[cluster_id] = cluster;
        }

        // Step 3: 构建邻接关系
        for (const auto& [cluster_id, neighbors] : cluster_adjacency) {
            Cluster* current_cluster = cluster_map[cluster_id];

            for (idx_t neighbor_id : neighbors) {
                current_cluster->adjacent_clusters.push_back(cluster_map[neighbor_id]);
            }
        }
    } /// Cluster Creation REGION End

    // 清空 mesh 原始数据
    mesh.m_vertexNormalData.clear();
    mesh.m_vertexPositionData.clear();
    mesh.m_vertexUVData.clear();
    mesh.m_indexData.clear();

    mesh.m_vertexNormalData.shrink_to_fit();
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

void clusterGrouping(const Mesh& mesh_ref, std::vector<ClusterGroup*>& cluster_groups) {
    const auto& clusters = mesh_ref.m_clusterList; // 获取 clusters 引用
    idx_t num_clusters = clusters.size();         // cluster 的数量

    if (num_clusters == 0) {
        std::cerr << "No clusters to group!" << std::endl;
        return;
    }

    // 构建 METIS 的邻接矩阵
    std::vector<idx_t> xadj(num_clusters + 1, 0);  // 邻接起始索引
    std::vector<idx_t> adjncy;                     // 邻接节点
    std::vector<idx_t> adjwgt;                     // 边权重（这里用 1）
    // 定义其他 METIS 参数
    idx_t num_constraints = 1;
    idx_t num_groups = static_cast<idx_t>(std::ceil(static_cast<float>(num_clusters) / MAXN_CLUSTER_IN_CLUSTERGROUP));
    std::vector<idx_t> partition_result(num_clusters, 0);
    idx_t objval;

    // 构建邻接矩阵
    for (size_t i = 0; i < clusters.size(); ++i) {
        Cluster* current_cluster = static_cast<Cluster*>(clusters[i]);
        xadj[i + 1] = xadj[i] + current_cluster->adjacent_clusters.size();

        for (Cluster* neighbor_cluster : current_cluster->adjacent_clusters) {
            // 获取邻接 cluster 的索引
            auto it = std::find(clusters.begin(), clusters.end(), neighbor_cluster);

            idx_t neighbor_idx = std::distance(clusters.begin(), it);
            adjncy.push_back(neighbor_idx);
            adjwgt.push_back(1); // 假设所有边权重为 1
        }
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

    // create ClusterGroup, set the error to the largest error in the clusters
    cluster_groups.clear();
    cluster_groups.reserve(num_groups);
    for (idx_t i = 0; i < num_groups; ++i) {
        ClusterGroup* cluster_group = new ClusterGroup();
        cluster_group->Error = 0.0f;
        cluster_groups.push_back(cluster_group);
    }

    // 将 cluster 分配到 cluster group
    for (size_t i = 0; i < num_clusters; ++i) {
        Cluster* current_cluster = static_cast<Cluster*>(clusters[i]);
        idx_t group_id = partition_result[i];
        ClusterGroup* current_group = cluster_groups[group_id];
        // 将 cluster 加入到 group 中
        current_group->clusters.push_back(current_cluster);
        // 更新 group 的 Error
        current_group->Error = std::max(current_group->Error, current_cluster->Error);
    }

    // // [DEBUG]: generate a random color for each clustergroup and assign to each cluster
    // for (ClusterGroup* group : cluster_groups) {
    //     glm::vec3 rdColor = genRdColor();
    //     for (Cluster* cluster : group->clusters) {
    //         cluster->rdColor = rdColor;
    //     }
    // }
}


