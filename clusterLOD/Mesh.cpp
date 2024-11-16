#include "Mesh.hpp"
#include "Cluster.hpp"
using namespace glm;
using namespace std;

#include "cs488-framework/Exception.hpp"
#include "cs488-framework/ObjFileDecoder.hpp"

#include "cs488-framework/OpenGLImport.hpp"
#include "cs488-framework/GlErrorCheck.hpp"

#include <unordered_map>
#include <iostream>

MeshInfoMap Mesh::s_meshInfoMap;

//----------------------------------------------------------------------------------------
// Default constructor
Mesh::Mesh()
{

}

//----------------------------------------------------------------------------------------
// Destructor
Mesh::~Mesh()
{

}

//----------------------------------------------------------------------------------------
template <typename T>
static void appendVector (
		std::vector<T> & dest,
		const std::vector<T> & source
) {
	// Increase capacity to hold source.size() more elements
	dest.reserve(dest.size() + source.size());

	dest.insert(dest.end(), source.begin(), source.end());
}


//----------------------------------------------------------------------------------------
Mesh::Mesh(
		std::initializer_list<ObjFilePath> objFileList
) {
	if(objFileList.size() != 1) {
		throw Exception("Error within Mesh: objFileList.size() != 1\n");
	}
	MeshId meshId;
	BatchInfo batchInfo;
	
	const ObjFilePath & objFile = *objFileList.begin();

#if ENABLE_IBO == true
	ObjFileDecoder::decode(objFile.c_str(), meshId, m_vertexPositionData, m_vertexNormalData, m_vertexUVData, m_indexData);
	batchInfo.numIndices = m_indexData.size();
#else
	ObjFileDecoder::decode(objFile.c_str(), meshId, m_vertexPositionData, m_vertexNormalData, m_vertexUVData);
	batchInfo.numIndices = m_vertexPositionData.size();
#endif

	batchInfo.startIndex = 0;
	s_meshInfoMap[meshId] = this;
}

namespace std {
    template<>
    struct hash<glm::vec3> {
        std::size_t operator()(const glm::vec3& v) const noexcept {
            // 简单组合哈希函数
            std::hash<float> floatHasher;
            size_t h1 = floatHasher(v.x);
            size_t h2 = floatHasher(v.y);
            size_t h3 = floatHasher(v.z);
            return h1 ^ (h2 << 1) ^ (h3 << 2); // 合并哈希值
        }
    };
}

Mesh::Mesh(const Mesh& other, const std::vector<unsigned int>& triangleIndices) {
    std::unordered_map<glm::vec3, unsigned int> vertexMap;
    
    for (unsigned int triIdx : triangleIndices) {
        for (size_t i = 0; i < 3; ++i) {
            unsigned int idx = other.m_indexData[3 * triIdx + i];
            
            glm::vec3 vertex = other.m_vertexPositionData[idx];
            if (vertexMap.find(vertex) == vertexMap.end()) {
                unsigned int newIdx = m_vertexPositionData.size();
                vertexMap[vertex] = newIdx;
                m_vertexPositionData.push_back(vertex);
                m_vertexNormalData.push_back(other.m_vertexNormalData[idx]);
                m_vertexUVData.push_back(other.m_vertexUVData[idx]);
            }
            m_indexData.push_back(vertexMap[vertex]);
        }
    }
}

void Mesh::uploadToGPU() {
    std::vector<Vertex> vertexData;
    vertexData.reserve(m_vertexPositionData.size());

    for (unsigned int i = 0; i < m_vertexPositionData.size(); i++) {
        Vertex vertex = {
            m_vertexPositionData[i],
            m_vertexNormalData[i],
            m_vertexUVData[i]
        };
        vertexData.push_back(vertex);
    }

    glGenVertexArrays(1, &m_vao);
    glBindVertexArray(m_vao);

    glGenBuffers(1, &m_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, m_vbo);
    glBufferData(GL_ARRAY_BUFFER, vertexData.size() * sizeof(Vertex), vertexData.data(), GL_STATIC_DRAW);

    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)0);

    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, normal));

    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, uv));

#if ENABLE_IBO == true
    glGenBuffers(1, &m_ibo);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_ibo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, m_indexData.size() * sizeof(unsigned int), m_indexData.data(), GL_STATIC_DRAW);
#endif

    glBindBuffer(GL_ARRAY_BUFFER, 0);

    glBindVertexArray(0);
    CHECK_GL_ERRORS;
}

void Mesh::removeFromGPU() {
    glDeleteBuffers(1, &m_vbo);
    glDeleteVertexArrays(1, &m_vao);
#if ENABLE_IBO == true
    glDeleteBuffers(1, &m_ibo);
#endif
}

void Mesh::draw() const {
    if(m_clusterList.size() > 0) {
        for (const Cluster& cluster : m_clusterList) {
            cluster.m_mesh->draw();
        }
        return;
    }

    glBindVertexArray(m_vao);
    
    #if ENABLE_IBO == true
        glDrawElements(GL_TRIANGLES, m_indexData.size(), GL_UNSIGNED_INT, 0);
    #else
        glDrawArrays(GL_TRIANGLES, 0, m_vertexPositionData.size());
    #endif
    
    glBindVertexArray(0);
}