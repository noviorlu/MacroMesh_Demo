#pragma once

#include "cs488-framework/BatchInfo.hpp"
#include "cs488-framework/OpenGLImport.hpp"
#include "cs488-framework/ShaderProgram.hpp"

#include "cs488-framework/Vertex.hpp"
#include <glm/glm.hpp>

#include <initializer_list>
#include <vector>
#include <unordered_map>
#include <string>
#include <iomanip> // For formatted output
#include <sstream> // For easier testing/debugging

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
class ClusterGroup;

/*
* Class for consolidating all vertex data within a list of .obj files.
*/
class Mesh {
public:
	Mesh();

	Mesh(std::initializer_list<ObjFilePath>  objFileList);

	Mesh(const Mesh& other, const std::vector<unsigned int>& triangleIndices);

	Mesh(const std::vector<Mesh>& mergeMeshes);

	~Mesh();

	void uploadToGPU();

	void removeFromGPU();

	virtual void draw(const ShaderProgram& shader) const;

	static MeshInfoMap s_meshInfoMap;

	std::vector<Vertex> m_vertexData;

	std::vector<unsigned int> m_indexData;

	std::vector<Cluster> m_clusterList;
	std::vector<ClusterGroup> m_clusterGroupList;

	GLuint m_vbo;
	GLuint m_vao;
	GLuint m_ibo;
};


class Cluster : public Mesh {
public:
    float Error;
    glm::vec3 rdColor;

    Cluster(float Error);
    ~Cluster() {}

	void draw(const ShaderProgram& shader) const override;
};

class ClusterGroup {
public:
    float Error;
    std::vector<Cluster*> clusters;
};
