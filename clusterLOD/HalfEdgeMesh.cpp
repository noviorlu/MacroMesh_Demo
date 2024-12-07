// #include "HalfEdgeMesh.hpp"

//  //#define DEBUG_HALFEDGEMESH
// #include <fstream>
// #include <sstream>
// #include <filesystem>
// #include <unordered_map>
// #include <unordered_set>
// #include <iostream>
// #include <glm/gtx/string_cast.hpp>
// #include "cs488-framework/ObjFileDecoder.hpp"

// #include <omp.h>

// HalfEdgeMesh::HalfEdgeMesh(const Mesh& mesh) {
//     m_faces.reserve(mesh.m_indexData.size() / 3);
//     m_vertices.reserve(mesh.m_vertexData.size());
//     m_edges.reserve(mesh.m_indexData.size());

// 	std::unordered_set<glm::vec3> uniqueVertices;
// 	for (auto& vertex : mesh.m_vertexData) {
// 		if (uniqueVertices.find(vertex.position) == uniqueVertices.end()) {
// 			uniqueVertices.insert(vertex.position);
// 		}
// 		else {
// 			std::cout << "Duplicate vertex found: " << glm::to_string(vertex.position) << std::endl;
//             assert(false);
// 		}
// 	}

//     for (size_t i = 0; i < mesh.m_vertexData.size(); ++i) {
//         m_vertices.push_back(new HalfVertex(mesh.m_vertexData[i]));
//     }

//     std::unordered_map<std::pair<int, int>, HalfEdge*> edgeMap;

//     for (size_t i = 0; i < mesh.m_indexData.size(); i += 3) {
//         int v0 = mesh.m_indexData[i];
//         int v1 = mesh.m_indexData[i + 1];
//         int v2 = mesh.m_indexData[i + 2];

//         HalfEdge* edge0 = new HalfEdge();
//         HalfEdge* edge1 = new HalfEdge();
//         HalfEdge* edge2 = new HalfEdge();

//         edge0->origin = m_vertices[v0];
//         edge1->origin = m_vertices[v1];
//         edge2->origin = m_vertices[v2];

//         edge0->next = edge1;
//         edge1->next = edge2;
//         edge2->next = edge0;

//         edge0->prev = edge2;
//         edge1->prev = edge0;
//         edge2->prev = edge1;

//         Face* face = new Face(edge0);
//         m_faces.push_back(face);

//         edge0->face = face;
//         edge1->face = face;
//         edge2->face = face;
//         assert(edge0 != nullptr && "Edge pointer 0 is null!");
//         assert(edge1 != nullptr && "Edge pointer 1 is null!");
//         assert(edge2 != nullptr && "Edge pointer 2 is null!");
//         m_edges.push_back(edge0);
//         m_edges.push_back(edge1);
//         m_edges.push_back(edge2);

//         auto addEdge = [&](int start, int end, HalfEdge* edge) {
//             auto twinIt = edgeMap.find({end, start});
//             if (twinIt != edgeMap.end()) {
//                 edge->twin = twinIt->second;
//                 twinIt->second->twin = edge;
//             } else {
//                 edge->twin = nullptr;
//                 edgeMap[{start, end}] = edge;
//             }
//         };

//         addEdge(v0, v1, edge0);
//         addEdge(v1, v2, edge1);
//         addEdge(v2, v0, edge2);
//     }
    


//     std::unordered_map<HalfVertex*, HalfEdge*> startVertexToBoundaryEdge;
// 	int size = m_edges.size();
//     for (int i = 0; i < size; i++) {
//         HalfEdge* edge = m_edges[i];

// 		if (edge->origin->edge != edge) edge->origin->edge = edge;

// 		if (edge->twin != nullptr) continue;
// 		HalfEdge* twinEdge = new HalfEdge();
// 		twinEdge->origin = edge->next->origin;
// 		twinEdge->twin = edge;
// 		edge->twin = twinEdge;
// 		twinEdge->face = nullptr;
// 		m_edges.push_back(twinEdge);
// 		startVertexToBoundaryEdge[twinEdge->origin] = twinEdge;
//     }

//     for (auto& [vertex, edge] : startVertexToBoundaryEdge) {
//         HalfEdge* twinEdge = edge;
//         HalfEdge* nextEdge = startVertexToBoundaryEdge[edge->twin->origin];
//         if (nextEdge) {
//             twinEdge->next = nextEdge;
//             nextEdge->prev = twinEdge;
//         }
//     }
// }

// HalfEdgeMesh::~HalfEdgeMesh() {
//     for (auto edge : m_edges) {
//         delete edge;
//     }

// 	for (auto vert : m_vertices) {
// 		delete vert;
// 	}
	
//     for (auto face : m_faces) {
// 		delete face;
// 	}
// }

// void HalfEdgeMesh::importMesh(const std::string& objFilePath, float error) {
//     m_error = error;

//     std::vector<glm::vec3> positions;
//     std::vector<glm::vec3> normals;
//     std::vector<glm::vec2> uvCoords;

//     auto starttime = std::chrono::high_resolution_clock::now();
//     ObjFileDecoder::decode(objFilePath.c_str(), m_name, positions, normals, uvCoords);

//     auto endtime = std::chrono::high_resolution_clock::now();
//     std::chrono::duration<double> elapsed = endtime - starttime;
//     std::cout << "Time taken to decode: " << elapsed.count() << "s" << std::endl;

//     starttime = std::chrono::high_resolution_clock::now();
// 	std::unordered_map<glm::vec3, HalfVertex*> vertexMap;
// 	std::unordered_map<std::pair<HalfVertex*, HalfVertex*>, HalfEdge*> edgeMap;
// 	for (size_t i = 0; i < positions.size(); i+=3) {
//         HalfVertex* vertex0;
//         HalfVertex* vertex1;
//         HalfVertex* vertex2;
//         if (vertexMap.find(positions[i]) != vertexMap.end()) vertex0 = vertexMap[positions[i]];
//         else {
// 			if (uvCoords.size() > 0) vertex0 = new HalfVertex(positions[i], normals[i], uvCoords[i]);
// 			else vertex0 = new HalfVertex(positions[i], normals[i], glm::vec2(0.0f));
//             vertexMap[positions[i]] = vertex0;
//         }
// 		if (vertexMap.find(positions[i + 1]) != vertexMap.end()) vertex1 = vertexMap[positions[i + 1]];
// 		else {
//             if (uvCoords.size() > 0) vertex1 = new HalfVertex(positions[i + 1], normals[i + 1], uvCoords[i + 1]);
// 			else vertex1 = new HalfVertex(positions[i + 1], normals[i + 1], glm::vec2(0.0f));
// 			vertexMap[positions[i + 1]] = vertex1;
// 		}
// 		if (vertexMap.find(positions[i + 2]) != vertexMap.end()) vertex2 = vertexMap[positions[i + 2]];
// 		else {
//             if (uvCoords.size() > 0) vertex2 = new HalfVertex(positions[i + 2], normals[i + 2], uvCoords[i + 2]);
// 			else vertex2 = new HalfVertex(positions[i + 2], normals[i + 2], glm::vec2(0.0f));
// 			vertexMap[positions[i + 2]] = vertex2;
// 		}

// 		assert(edgeMap.find({ vertex0, vertex1 }) == edgeMap.end());
// 		assert(edgeMap.find({ vertex1, vertex2 }) == edgeMap.end());
// 		assert(edgeMap.find({ vertex2, vertex0 }) == edgeMap.end());

//         HalfEdge* edge0 = new HalfEdge();
//         HalfEdge* edge1 = new HalfEdge();
//         HalfEdge* edge2 = new HalfEdge();
// 		edgeMap[{vertex0, vertex1}] = edge0;
// 		edgeMap[{vertex1, vertex2}] = edge1;
// 		edgeMap[{vertex2, vertex0}] = edge2;

// 		Face* face = new Face(edge0);
//         edge0->face = face;
//         edge1->face = face;
//         edge2->face = face;
// 		m_faces.push_back(face);

// 		edge0->origin = vertex0;
// 		edge1->origin = vertex1;
// 		edge2->origin = vertex2;

//         vertex0->edge = edge0;
//         vertex1->edge = edge1;
//         vertex2->edge = edge2;

//         if (edgeMap.find({ vertex1, vertex0 }) != edgeMap.end()) {
// 			edge0->twin = edgeMap[{vertex1, vertex0}];
// 			edgeMap[{vertex1, vertex0}]->twin = edge0;
//         }
// 		if (edgeMap.find({ vertex2, vertex1 }) != edgeMap.end()) {
// 			edge1->twin = edgeMap[{vertex2, vertex1}];
// 			edgeMap[{vertex2, vertex1}]->twin = edge1;
// 		}
// 		if (edgeMap.find({ vertex0, vertex2 }) != edgeMap.end()) {
// 			edge2->twin = edgeMap[{vertex0, vertex2}];
// 			edgeMap[{vertex0, vertex2}]->twin = edge2;
// 		}

// 		edge0->next = edge1;
// 		edge1->next = edge2;
// 		edge2->next = edge0;

// 		edge0->prev = edge2;
// 		edge1->prev = edge0;
// 		edge2->prev = edge1;
// 	}

// 	for (auto& [vertexPair, edge] : edgeMap) {
// 		m_edges.push_back(edge);
// 	}
// 	for (auto& [pos, vert] : vertexMap) {
// 		m_vertices.push_back(vert);
// 	}

//     // find all edges that doesnot have twin
// 	std::unordered_map<HalfVertex*, HalfEdge*> startVertexToBoundaryEdge;
//     for (int i = m_edges.size()-1; i >= 0; --i) {
// 		auto edge = m_edges[i];
//         if (edge->twin != nullptr) continue;
//         HalfEdge* twinEdge = new HalfEdge();

//         twinEdge->twin = edge;
//         edge->twin = twinEdge;

//         twinEdge->origin = edge->next->origin;

//         startVertexToBoundaryEdge[twinEdge->origin] = twinEdge;
//         m_edges.push_back(twinEdge);
//     }

// 	for (auto& [vertex, edge] : startVertexToBoundaryEdge) {
// 		HalfEdge* twinEdge = edge;
// 		HalfEdge* nextEdge = startVertexToBoundaryEdge[edge->twin->origin];
		
// 		assert(nextEdge != nullptr);
//         twinEdge->next = nextEdge;
// 		nextEdge->prev = twinEdge;
// 	}

//     std::cout << "Imported EMesh from " << objFilePath << std::endl;
//     endtime = std::chrono::high_resolution_clock::now();
//     elapsed = endtime - starttime;
//     std::cout << "Time taken to import: " << elapsed.count() << "s" << std::endl;
// }

// void HalfEdgeMesh::exportMesh(Mesh& mesh) {
//     mesh.m_vertexData.clear();
//     mesh.m_indexData.clear();
//     std::unordered_map<HalfVertex*, unsigned int> vertexIndexMap;
//     for (auto vertex : m_vertices) {
// 		if (vertex->isValid == false) continue;
//         mesh.m_vertexData.emplace_back(
//             vertex->position,
//             vertex->normal,
//             vertex->uv
//         );
//         vertexIndexMap[vertex] = mesh.m_vertexData.size();
//     }

//     for (auto face : m_faces) {
// 		if (face->isValid == false) continue;
//         HalfEdge* edge = face->edge;
//         do {
//             assert(edge->isValid);
//             assert(vertexIndexMap.find(edge->origin) == vertexIndexMap.end());
//             mesh.m_indexData.push_back(vertexIndexMap[edge->origin]);
//             edge = edge->next;
//         } while (edge != face->edge);
//     }
// }

// void HalfEdgeMesh::exportMesh(std::vector<Cluster>& clusters, std::vector<ClusterGroup>& clusterGroups) {
//     clusters.clear();
//     clusterGroups.clear();

//     std::vector<Cluster> localClusters(m_clusterOffsets.size());

//     #pragma omp parallel for
//     for (size_t i = 0; i < m_clusterOffsets.size(); ++i) {
//         size_t startFace = m_clusterOffsets[i];
//         size_t endFace = (i + 1 < m_clusterOffsets.size()) ? m_clusterOffsets[i + 1] : m_faces.size();

//         Cluster cluster(0.0f);
//         std::unordered_map<HalfVertex*, unsigned int> vertexIndexMap; // 每个线程自己的map
//         unsigned int clusterVertexIndex = 0;

//         for (size_t faceIdx = startFace; faceIdx < endFace; ++faceIdx) {
//             HalfEdge* edge = m_faces[faceIdx]->edge;
//             do {
//                 HalfVertex* vertex = edge->origin;
//                 if (vertexIndexMap.find(vertex) == vertexIndexMap.end()) {
//                     vertexIndexMap[vertex] = clusterVertexIndex++;
//                     cluster.m_vertexData.emplace_back(
//                         vertex->position,
//                         vertex->normal,
//                         vertex->uv
//                     );
//                 }
//                 cluster.m_indexData.push_back(vertexIndexMap[vertex]);
//                 edge = edge->next;
//             } while (edge != m_faces[faceIdx]->edge);
//         }
//         localClusters[i] = std::move(cluster);
//     }

//     for (size_t i = 0; i < localClusters.size(); ++i) {
//         clusters.push_back(std::move(localClusters[i]));
//     }

//     std::unordered_map<int, float> clusterGroupClr;
//     for (size_t i = 0; i < m_clusterGroupResult.size(); ++i) {
//         int group = m_clusterGroupResult[i];
//         if (clusterGroupClr.find(group) == clusterGroupClr.end()) {
//             clusterGroupClr[group] = static_cast<float>(rand()) / RAND_MAX;
//         }
//     }

//     #pragma omp parallel for
//     for (int i = 0; i < clusters.size(); i++) {
//         float s = 0.8f;
//         float v = 0.5f + static_cast<float>(rand()) / (2.0f * RAND_MAX);
//         clusters[i].rdColor = HSVtoRGB(clusterGroupClr[m_clusterGroupResult[i]], s, v);
//     }
// }

// void HalfEdgeMesh::exportMeshToObjFiles(const std::string& folderPath) {
//     if (!std::filesystem::exists(folderPath)) {
//         std::filesystem::create_directories(folderPath);
//     }

//     #pragma omp parallel for
//     for (size_t clusterIndex = 0; clusterIndex < m_clusterGroupOffset.size() - 1; ++clusterIndex) {
//         size_t startIdx = m_clusterGroupOffset[clusterIndex];
//         size_t endIdx = m_clusterGroupOffset[clusterIndex + 1];

//         // 构造文件名
//         std::string fileName = folderPath + "/clusterGroup_" + std::to_string(clusterIndex) + ".obj";
//         std::ofstream outFile(fileName);

//         if (!outFile.is_open()) {
//             #pragma omp critical
//             std::cerr << "Failed to open file: " << fileName << std::endl;
//             continue;
//         }

//         // 存储唯一顶点
//         std::unordered_map<HalfVertex*, size_t> vertexMap;
//         size_t vertexIndex = 1;

//         // 写入顶点信息
//         for (size_t i = startIdx; i < endIdx; ++i) {
//             const auto& face = m_faces[i];
//             const HalfEdge* edge = face->edge;
//             do {
//                 HalfVertex* vertex = edge->origin;
//                 if (vertexMap.find(vertex) == vertexMap.end()) {
//                     vertexMap[vertex] = vertexIndex++;
//                     outFile << "v " << vertex->position.x << " " << vertex->position.y << " " << vertex->position.z << "\n";
//                 }
//                 edge = edge->next;
//             } while (edge != face->edge);
//         }

//         // 写入面信息
//         for (size_t i = startIdx; i < endIdx; ++i) {
//             const auto& face = m_faces[i];
//             outFile << "f";
//             const HalfEdge* edge = face->edge;
//             do {
//                 outFile << " " << vertexMap[edge->origin];
//                 edge = edge->next;
//             } while (edge != face->edge);
//             outFile << "\n";
//         }

//         outFile.close();
//     }
// }

// void HalfEdgeMesh::exportMesh(const std::string& objFilePath) {
//     std::ofstream objFile(objFilePath);
//     if (!objFile.is_open()) {
//         throw std::runtime_error("Failed to open file: " + objFilePath);
//     }

//     std::unordered_map<HalfVertex*, unsigned int> vertexIndexMap;
//     unsigned int currentVertexIndex = 1;

//     for (auto vertex : m_vertices) {
//         objFile << "v "
//             << vertex->position.x << " "
//             << vertex->position.y << " "
//             << vertex->position.z << "\n";
//         vertexIndexMap[vertex] = currentVertexIndex++;
//     }

//     for (auto vertex : m_vertices) {
//         objFile << "vt "
//             << vertex->uv.x << " "
//             << vertex->uv.y << "\n";
//     }

//     for (auto vertex : m_vertices) {
//         objFile << "vn "
//             << vertex->normal.x << " "
//             << vertex->normal.y << " "
//             << vertex->normal.z << "\n";
//     }

//     for (auto face : m_faces) {
//         objFile << "f";
//         HalfEdge* edge = face->edge;
//         do {
//             unsigned int vertexIndex = vertexIndexMap[edge->origin];
// 			if (vertexIndex == 0) {
// 				std::cerr << "[ERROR] Vertex index not found in map." << std::endl;
//                 assert(false);
// 			}
//             objFile << " " << vertexIndex << "/" << vertexIndex << "/" << vertexIndex;
//             edge = edge->next;
//         } while (edge != face->edge);
//         objFile << "\n";
//     }

//     objFile.close();
// }

// void HalfEdgeMesh::HalfEdgeMeshValidation() {
//     std::vector<std::string> errorMessages;

// 	std::unordered_set<HalfVertex*> vertexSet;
// 	std::unordered_set<HalfEdge*> edgeSet;
// 	std::unordered_set<Face*> faceSet;

// 	for (auto vertex : m_vertices) {
//         if (vertex->isValid == false) continue;
// 		if (vertexSet.find(vertex) != vertexSet.end()) {
// 			errorMessages.emplace_back("[ERROR] Duplicate vertex found in the mesh.");
// 		}
// 		vertexSet.insert(vertex);
// 	}

// 	for (auto edge : m_edges) {
// 		if (edge->isValid == false) continue;
// 		if (edgeSet.find(edge) != edgeSet.end()) {
// 			errorMessages.emplace_back("[ERROR] Duplicate edge found in the mesh.");
// 		}
// 		edgeSet.insert(edge);
// 	}

// 	for (auto face : m_faces) {
// 		if (face->isValid == false) continue;
// 		if (faceSet.find(face) != faceSet.end()) {
// 			errorMessages.emplace_back("[ERROR] Duplicate face found in the mesh.");
// 		}
// 		faceSet.insert(face);
// 	}

//     for (size_t i = 0; i < m_vertices.size(); ++i) {
//         HalfVertex* vertex = m_vertices[i];
//         if (vertex->isValid == false) continue;
//         if (!vertex) {
//             errorMessages.emplace_back("[ERROR] Vertex at index " + std::to_string(i) + " is null.");
//             continue;
//         }
//         if (!vertex->edge) {
//             errorMessages.emplace_back("[ERROR] Vertex at index " + std::to_string(i) + " has a null edge pointer.");
//         } else if (vertex->edge->origin != vertex) {
//             errorMessages.emplace_back("[ERROR] Vertex at index " + std::to_string(i) 
//                 + " is not the origin of its associated edge.");
// 		}
// 		else if (edgeSet.find(vertex->edge) == edgeSet.end()) {
// 			//errorMessages.emplace_back("[ERROR] Vertex at index " + std::to_string(i)
// 			//	+ " has an edge that is not in the edge list.");
//             errorMessages.emplace_back("[ERROR] Vertex at index " + std::to_string(i)
//                 + " has an edge that is not in the edge list.");
// 		}
//         // check infinity Loop
//         std::unordered_set<HalfEdge*> setting;
//         HalfEdge* start, * cur;
//         start = cur = vertex->edge;
//         int counter = 0;
//         do {
//             assert(cur->origin == vertex);
//             assert(cur->isValid);
//             cur = cur->twin->next;
//             if (setting.find(cur) != setting.end()) assert(false);
//             setting.insert(cur);
//             if (counter++ > 1000) {
//                 errorMessages.emplace_back("[ERROR] Vertex " + std::to_string(i) + " has an infinity loop.");
//                 break;
//             }
//         } while (cur != start);
//     }

//     for (size_t i = 0; i < m_edges.size(); ++i) {
//         HalfEdge* edge = m_edges[i];

//         if (edge->isValid == false) continue;
//         if (!edge) {
//             errorMessages.emplace_back("[ERROR] Edge at index " + std::to_string(i) + " is null.");
//             continue;
//         }
//         if (!edge->origin) {
//             errorMessages.emplace_back("[ERROR] Edge at index " + std::to_string(i) + " has a null origin vertex.");
//         }
// 		if (vertexSet.find(edge->origin) == vertexSet.end()) {
// 			errorMessages.emplace_back("[ERROR] Edge at index " + std::to_string(i) + " has an origin vertex that is not in the vertex list.");
// 		}
// 		//if (edge->face != nullptr && faceSet.find(edge->face) == faceSet.end()) {
// 		//	errorMessages.emplace_back("[ERROR] Edge at index " + std::to_string(i) + " has a face that is not in the face list.");
// 		//}
//         if (!edge->twin) {
//             errorMessages.emplace_back("[ERROR] Edge at index " + std::to_string(i) + " has no twin edge.");
//         } else if (edge->twin->twin != edge) {
//             errorMessages.emplace_back("[ERROR] Edge at index " + std::to_string(i) 
//                 + " has a twin whose twin does not point back to this edge.");
//         }
//         if (!edge->next) {
//             errorMessages.emplace_back("[ERROR] Edge at index " + std::to_string(i) + " has no next edge.");
//         }
//         if (!edge->prev) {
//             errorMessages.emplace_back("[ERROR] Edge at index " + std::to_string(i) + " has no previous edge.");
//         } 
//     }


//     std::unordered_map<TripleHalfVertex, std::vector<Face*>, TripleHalfVertexHash, TripleHalfVertexEqual> faceMap;

//     for (size_t i = 0; i < m_faces.size(); ++i) {
//         Face* face = m_faces[i];

//         if (face->isValid == false) continue;
//         if (!face) {
//             errorMessages.emplace_back("[ERROR] Face at index " + std::to_string(i) + " is null.");
//             continue;
//         }
//         if (!face->edge) {
//             errorMessages.emplace_back("[ERROR] Face at index " + std::to_string(i) + " has no associated edge.");
//         } else {
//             HalfEdge* startEdge = face->edge;
//             HalfEdge* currentEdge = startEdge;
//             size_t edgeCount = 0;
//             do {
//                 if (!currentEdge) {
//                     errorMessages.emplace_back("[ERROR] Null edge encountered in face at index " + std::to_string(i) + ".");
//                     break;
//                 }
//                 if (currentEdge->face != face) {
//                     errorMessages.emplace_back("[ERROR] Edge in face at index " + std::to_string(i) 
//                         + " does not point back to the correct face.");
//                 }
// 				if (edgeSet.find(currentEdge) == edgeSet.end()) {
// 					errorMessages.emplace_back("[ERROR] Edge in face at index " + std::to_string(i)
// 						+ " is not in the edge list.");
// 				}
//                 currentEdge = currentEdge->next;
//                 edgeCount++;
//                 if (edgeCount > 3) {
//                     errorMessages.emplace_back("[ERROR] not triangle " + std::to_string(i) + ".");
//                     break;
//                 }
//             } while (currentEdge != startEdge);
//         }
//     }

//     // If there are any errors, print them and assert
//     if (!errorMessages.empty()) {
//         for (const auto& msg : errorMessages) {
//             std::cerr << msg << std::endl;
//         }
//         assert(false && "HalfEdgeMesh validation failed. Check error messages.");
//     }
// }

// void HalfEdgeMesh::HalfEdgeMeshPrint() {
//     std::unordered_map<HalfVertex*, size_t> vertexIndexMap;
//     for (size_t i = 0; i < m_vertices.size(); ++i) {
//         vertexIndexMap[m_vertices[i]] = i;
//     }

//     std::unordered_map<HalfEdge*, size_t> edgeIndexMap;
//     for (size_t i = 0; i < m_edges.size(); ++i) {
//         edgeIndexMap[m_edges[i]] = i;
//     }

//     std::unordered_map<Face*, size_t> faceIndexMap;
//     for (size_t i = 0; i < m_faces.size(); ++i) {
//         faceIndexMap[m_faces[i]] = i;
//     }

//     std::cout << "Vertices:" << std::endl;
//     for (size_t i = 0; i < m_vertices.size(); ++i) {
//         HalfVertex* vertex = m_vertices[i];
// 		if (vertex->isValid == false) continue;
//         std::cout << "Vertex " << i << "("<< vertex << "): position = ("
//             << vertex->position.x << ", " << vertex->position.y << ", " << vertex->position.z
//             << "), edge = " << (vertex->edge ? edgeIndexMap[vertex->edge] : -1) << std::endl;
//     }

//     std::cout << "Edges:" << std::endl;
//     for (size_t i = 0; i < m_edges.size(); ++i) {
//         HalfEdge* edge = m_edges[i];
//         if (edge->isValid == false) continue;
//         std::cout << "Edge " << i << "(" << edge << "): origin = "
//             << (edge->origin ? vertexIndexMap[edge->origin] : -1)
//             << ", twin = " << (edge->twin ? edgeIndexMap[edge->twin] : -1)
//             << ", next = " << (edge->next ? edgeIndexMap[edge->next] : -1)
//             << ", prev = " << (edge->prev ? edgeIndexMap[edge->prev] : -1) << std::endl;
//     }

//     std::cout << "Faces:" << std::endl;
//     for (size_t i = 0; i < m_faces.size(); ++i) {
//         Face* face = m_faces[i];
//         if (face->isValid == false) continue;

//         std::cout << "Face " << i << ": edge = "
//             << (face->edge ? edgeIndexMap[face->edge] : -1);
    
//         // print inner edge of this face
// 		HalfEdge* startEdge = face->edge;
// 		HalfEdge* currentEdge = startEdge;
//         std::cout << " Edge Loop: ";
//         do {
//             std::cout << edgeIndexMap[currentEdge]  << " ";
// 			currentEdge = currentEdge->next;
// 		} while (currentEdge != startEdge);
//         std::cout << std::endl;
//     }
// }






// void EMesh::importEMesh(const std::string& objFilePath) {
//     std::vector<glm::vec3> positions;
//     std::vector<glm::vec3> normals;
//     std::vector<glm::vec2> uvCoords;

//     auto starttime = std::chrono::high_resolution_clock::now();
//     ObjFileDecoder::decode(objFilePath.c_str(), m_name, positions, normals, uvCoords);

//     auto endtime = std::chrono::high_resolution_clock::now();
//     std::chrono::duration<double> elapsed = endtime - starttime;
//     std::cout << "Time taken to decode: " << elapsed.count() << "s" << std::endl;

//     starttime = std::chrono::high_resolution_clock::now();

//     std::unordered_map<glm::vec3, EVertex*> vertexMap;

//     m_edgeMap.clear();
//     m_faces.clear();
//     m_vertices.clear();

//     m_edgeMap.reserve(positions.size());
//     m_faces.reserve(positions.size() / 3);
//     m_vertices.reserve(positions.size());

//     auto getOrCreateVertex = [&](const glm::vec3& pos, const glm::vec3& norm, const glm::vec2& uv) -> EVertex* {
//         if (vertexMap.find(pos) == vertexMap.end()) {
// 			EVertex* vertex = new EVertex(pos, norm, uv);
//             m_vertices.emplace_back(vertex);
//             vertexMap[pos] = vertex;
//         }
//         return vertexMap[pos];
//     };
//     auto createEdge = [&](EVertex* v1, EVertex* v2, EFace* face) {
//         VertexPair pair(v1, v2);
// 		EEdge* edge = nullptr;
//         if (m_edgeMap.find(pair) == m_edgeMap.end()) {
// 			edge = new EEdge(pair, face);
//             m_edgeMap.emplace(pair, edge);
//             v1->addEdge(edge);
//             v2->addEdge(edge);
// 		}
// 		else {
// 			edge = m_edgeMap[pair];

// 			edge->addFace(face);
// 		}
//     };

//     if (uvCoords.size() > 0) {
//         for (size_t i = 0; i < positions.size(); i += 3) {
//             EVertex* v0 = getOrCreateVertex(positions[i], normals[i], uvCoords[i]);
//             EVertex* v1 = getOrCreateVertex(positions[i + 1], normals[i + 1], uvCoords[i + 1]);
//             EVertex* v2 = getOrCreateVertex(positions[i + 2], normals[i + 2], uvCoords[i + 2]);

//             EFace* face = new EFace(v0, v1, v2);
//             m_faces.push_back(face);

//             createEdge(v0, v1, face);
//             createEdge(v1, v2, face);
//             createEdge(v2, v0, face);
//         }
//     }
//     else {
//         for (size_t i = 0; i < positions.size(); i += 3) {
//             EVertex* v0 = getOrCreateVertex(positions[i], normals[i], glm::vec2(0));
//             EVertex* v1 = getOrCreateVertex(positions[i + 1], normals[i + 1], glm::vec2(0));
//             EVertex* v2 = getOrCreateVertex(positions[i + 2], normals[i + 2], glm::vec2(0));
            
//             EFace* face = new EFace(v0, v1, v2);
//             m_faces.push_back(face);

//             createEdge(v0, v1, face);
//             createEdge(v1, v2, face);
//             createEdge(v2, v0, face);
//         }
//     }

//     for(auto& [pair, edge] : m_edgeMap) {
//         if(edge->faces.size() == 1) {
//             edge->isBoundary = true;
//         }
//     }

//     std::cout << "Imported EMesh from " << objFilePath << std::endl;
//     endtime = std::chrono::high_resolution_clock::now();
//     elapsed = endtime - starttime;
//     std::cout << "Time taken to import: " << elapsed.count() << "s" << std::endl;
// }

// void EMesh::exportEMesh(const std::string& objFilePath) {
//     std::ofstream outFile(objFilePath);
//     if (!outFile.is_open()) {
//         std::cerr << "Failed to open file for writing: " << objFilePath << std::endl;
//         return;
//     }

// 	std::unordered_map<EVertex*, int> vertexIndexMap;
//     int counter = 1;
//     for (EVertex* vertex : m_vertices) {
//         if (vertex->edges.size() == 0) continue;
//         vertexIndexMap[vertex] = counter;
//         counter++;
//         outFile << "v " << vertex->position.x << " " << vertex->position.y << " " << vertex->position.z << "\n";
//     }

//     for (EVertex* vertex : m_vertices) {
//         if (vertex->edges.size() == 0) continue;
//         outFile << "vn " << vertex->normal.x << " " << vertex->normal.y << " " << vertex->normal.z << "\n";
//     }

//     for (EVertex* vertex : m_vertices) {
//         if (vertex->edges.size() == 0) continue;
//         outFile << "vt " << vertex->uv.x << " " << vertex->uv.y << "\n";
//     }

//     for (const auto& face : m_faces) {
//         if (!face->isValid) continue;

//         outFile << "f ";
//         for (int i = 0; i < 3; ++i) {
//             auto vertex = face->vertices[i];
//             int vertexIndex = vertexIndexMap[vertex];
//             outFile << vertexIndex << "/" << vertexIndex << "/" << vertexIndex << " ";
//         }
//         outFile << "\n";
//     }

//     outFile.close();
// }

// void EMesh::exportEMeshToObjFiles(const std::string& folderPath) {
//     if (!std::filesystem::exists(folderPath)) {
//         std::filesystem::create_directories(folderPath);
//     }

//     std::unordered_map<int, std::vector<EFace*>> clusterGroupMap;

//     for (size_t i = 0; i < m_clusterOffsets.size() - 1; ++i) {
//         size_t startIdx = m_clusterOffsets[i];
//         size_t endIdx = m_clusterOffsets[i + 1];
//         int group = m_clusterGroupResult[i];

//         for (size_t j = startIdx; j < endIdx; ++j) {
//             clusterGroupMap[group].push_back(m_faces[j]);
//         }
//     }

//     for (const auto& [group, faces] : clusterGroupMap) {
//         std::string fileName = folderPath + "/clusterGroup_" + std::to_string(group) + ".obj";
//         std::ofstream outFile(fileName);

//         if (!outFile.is_open()) {
//             std::cerr << "Failed to open file: " << fileName << std::endl;
//             continue;
//         }

//         std::unordered_map<EVertex*, size_t> vertexMap;
//         size_t vertexIndex = 1;
//         for (const auto& face : faces) {
//             for (const auto& vertex : face->vertices) {
//                 if (vertexMap.find(vertex) == vertexMap.end()) {
//                     vertexMap[vertex] = vertexIndex++;
//                     outFile << "v " << vertex->position.x << " " << vertex->position.y << " " << vertex->position.z << "\n";
//                 }
//             }
//         }

//         for (const auto& face : faces) {
//             outFile << "f";
//             for (const auto& vertex : face->vertices) {
//                 outFile << " " << vertexMap[vertex];
//             }
//             outFile << "\n";
//         }

//         outFile.close();
//         std::cout << "Exported clusterGroup " << group << " to " << fileName << std::endl;
//     }
// }

// void EMesh::exportEMesh(std::vector<Cluster>& clusters, std::vector<ClusterGroup>& clusterGroups) {
//     clusters.clear();
//     clusterGroups.clear();

//     for (size_t i = 0; i < m_clusterOffsets.size(); ++i) {
//         size_t startFace = m_clusterOffsets[i];
//         size_t endFace = (i + 1 < m_clusterOffsets.size()) ? m_clusterOffsets[i + 1] : m_faces.size();

//         Cluster cluster(0.0f);
//         std::unordered_map<EVertex*, unsigned int> vertexIndexMap;
//         unsigned int clusterVertexIndex = 0;

//         for (size_t faceIdx = startFace; faceIdx < endFace; ++faceIdx) {
//             EFace* face = m_faces[faceIdx];

//             for (EVertex* vertex : face->vertices) {
//                 if (vertexIndexMap.find(vertex) == vertexIndexMap.end()) {
//                     vertexIndexMap[vertex] = clusterVertexIndex++;
//                     cluster.m_vertexData.emplace_back(
//                         vertex->position,
//                         vertex->normal,
//                         vertex->uv
//                     );
//                 }
//                 cluster.m_indexData.push_back(vertexIndexMap[vertex]);
//             }
//         }
//         clusters.push_back(std::move(cluster));
//     }

//     std::unordered_map<int, float> clusterGroupClr;
// 	for (size_t i = 0; i < m_clusterGroupResult.size(); ++i) {
// 		int group = m_clusterGroupResult[i];
// 		if (clusterGroupClr.find(group) == clusterGroupClr.end()) {
//             clusterGroupClr[group] = static_cast<float>(rand()) / RAND_MAX;
// 		}
// 	}

//     for (int i = 0; i < clusters.size(); i++) {
//         float s = 0.8f;
//         float v = 0.5f + static_cast<float>(rand()) / (2.0f * RAND_MAX);
// 		clusters[i].rdColor = HSVtoRGB(clusterGroupClr[m_clusterGroupResult[i]], s, v);
//     }

//     //clusterGroups.reserve(m_clusterGroupCount);
//     //for(int i = 0; i < m_clusterGroupCount; ++i) {
//     //    ClusterGroup group(0.0f);
//     //    clusterGroups.push_back(std::move(group));
//     //}
//     //for(int i = 0; i < clusters.size(); ++i) {
//     //    clusterGroups[m_clusterGroupResult[i]].clusters.push_back(&clusters[i]);
//     //}

//     //for (auto& group : clusterGroups) {
//     //    group.setGroupColor();
//     //}
// }

// void EMesh::validate(){
//     // face should not have empty vertices
//     for (auto& face : m_faces) {
//         if (face->isValid == false) continue;
//         if (face->vertices[0] == nullptr || face->vertices[1] == nullptr || face->vertices[2] == nullptr) {
//             std::cerr << "[ERROR] Face has empty vertices." << std::endl;
//             assert(false);
//         }

//         // validate if face's edge contains face
//         for (int i = 0; i < 3; ++i) {
//             EVertex* v1 = face->vertices[i];
//             EVertex* v2 = face->vertices[(i + 1) % 3];
//             VertexPair pair(v1, v2);
// 			if (m_edgeMap.find(pair) == m_edgeMap.end()) {
// 				std::cerr << "[ERROR] Face has an edge that is not in the edge list." << std::endl;
// 				assert(false);
// 			}

//             if(!m_edgeMap[pair]->containFace(face)) {
//                 std::cerr << "[<ERROR>] Face has an edge that does not contain the face. " << face << std::endl;
//                 assert(false);
//             }
//         }
//     }

//     // edge should not have empty vertices
//     for (auto& [pair, edge] : m_edgeMap) {
//         if (edge->isValid == false) continue;
//         if (edge->vertices.v1 == nullptr || edge->vertices.v2 == nullptr) {
//             std::cerr << "[ERROR] Edge has empty vertices." << std::endl;
//             assert(false);
//         }

//         for (auto& face : edge->faces) {
//             if (face == nullptr) {
//                 std::cerr << "[ERROR] Edge has empty face." << std::endl;
//                 assert(false);
//             }

//             // check if face in m_faces
//             if (std::find(m_faces.begin(), m_faces.end(), face) == m_faces.end()) {
//                 std::cerr << "[ERROR] Edge has a face that is not in the face list." << std::endl;
//                 assert(false);
//             }
            

//             if(face->vertices[0] == nullptr || face->vertices[1] == nullptr || face->vertices[2] == nullptr) {
//                 std::cerr << "[ERROR] Face has empty vertices." << std::endl;
//                 assert(false);
//             }
//         }
//     }
// }

// void EMesh::print() {
//     std::cout << "Vertices:" << std::endl;
//     for (const auto& vertex : m_vertices) {
//         std::cout << "Vertex " << vertex << ": ";
//         for (const auto& edge : vertex->edges) {
//             std::cout << edge << " ";
//         }
//         std::cout << std::endl;
//     }

//     std::cout << "Edges:" << std::endl;
//     for (const auto& [pair, edge] : m_edgeMap) {
//         if (edge->isValid == false) continue;
//         std::cout << "Edge " << edge << ": (" << pair.v1 << ", " << pair.v2 << ") Faces: ";
//         for (const auto& face : edge->faces) {
//             std::cout << face << " ";
//         }
//         std::cout << std::endl;
//     }

//     std::cout << "Faces:" << std::endl;
//     for (const auto& face : m_faces) {
//         if (face->isValid == false) continue;
//         std::cout << "Face " << face << ": ";
//         std::cout << m_edgeMap[VertexPair(face->vertices[0], face->vertices[1])] << " ";
//         std::cout << m_edgeMap[VertexPair(face->vertices[1], face->vertices[2])] << " ";
//         std::cout << m_edgeMap[VertexPair(face->vertices[2], face->vertices[0])] << " ";
//         std::cout << std::endl;
//     }
// }
