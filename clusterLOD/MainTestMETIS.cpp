// #include <metis.h>
// #include <iostream>

// int main() {
//     // 定义图的基本属性
//     idx_t nVertices = 6;  // 顶点数
//     idx_t nEdges = 7;     // 边数
//     idx_t nParts = 2;     // 要划分的部分数
//     idx_t nCon = 1;       // 图的约束数，通常为 1

//     // 定义图的结构
//     idx_t xadj[7] = {0, 2, 5, 7, 9, 12, 14};  // 邻接表索引
//     idx_t adjncy[14] = {1, 3, 0, 2, 4, 1, 5, 2, 4, 0, 3, 5, 1, 4};  // 邻接表

//     // 定义划分输出数组
//     idx_t part[6];

//     // 定义 tpwgts（每个部分的目标权重）和 ubvec（不平衡向量）
//     real_t tpwgts[2] = {0.5, 0.5};  // 两部分各占50%
//     real_t ubvec[1] = {1.05};  // 允许 5% 的不平衡

//     // 定义一些可选参数（可以根据需要调整）
//     idx_t options[METIS_NOPTIONS];
//     METIS_SetDefaultOptions(options);
//     options[METIS_OPTION_NUMBERING] = 0;  // 使用C风格的从0开始的编号

//     // 定义 edgecut 用于存储切割的边数
//     idx_t edgecut;

//     // 调用 METIS_PartGraphKway 进行划分
//     int result = METIS_PartGraphKway(&nVertices, &nCon, xadj, adjncy, 
//                                      nullptr, nullptr, nullptr, &nParts, 
//                                      tpwgts, ubvec, options, &edgecut, part);

//     // 检查结果
//     if (result == METIS_OK) {
//         std::cout << "Graph partitioning successful. Edgecut: " << edgecut << std::endl;
//         std::cout << "Partition results:\n";
//         for (int i = 0; i < nVertices; i++) {
//             std::cout << "Vertex " << i << " is in partition " << part[i] << std::endl;
//         }
//     } else {
//         std::cerr << "METIS failed to partition the graph.\n";
//     }

//     return 0;
// }
