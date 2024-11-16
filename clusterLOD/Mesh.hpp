#pragma once

#include "cs488-framework/BatchInfo.hpp"
#include "cs488-framework/OpenGLImport.hpp"

#include "cs488-framework/Vertex.hpp"
#include <glm/glm.hpp>

#include <initializer_list>
#include <vector>
#include <unordered_map>
#include <string>

#define ENABLE_IBO true

// String identifier for a mesh.
typedef std::string  MeshId;

// File path to a .obj file.
typedef std::string ObjFilePath;


// MeshInfoMap is an associative container that maps a unique MeshId to a BatchInfo
// object. Each BatchInfo object contains an index offset and the number of indices
// required to render the mesh with identifier MeshId.
class Mesh;
typedef std::unordered_map<MeshId, Mesh*>  MeshInfoMap;

class Cluster;

/*
* Class for consolidating all vertex data within a list of .obj files.
*/
class Mesh {
public:
	Mesh();

	Mesh(std::initializer_list<ObjFilePath>  objFileList);

	Mesh(const Mesh& other, const std::vector<unsigned int>& triangleIndices);

	~Mesh();

	void uploadToGPU();

	void removeFromGPU();

	void draw() const;

	static MeshInfoMap s_meshInfoMap;

	std::vector<glm::vec3> m_vertexPositionData;
	std::vector<glm::vec3> m_vertexNormalData;
	std::vector<glm::vec2> m_vertexUVData;

	std::vector<unsigned int> m_indexData;

	std::vector<Cluster> m_clusterList;

	GLuint m_vbo;
	GLuint m_vao;
	GLuint m_ibo;
};


