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
	std::vector<glm::vec3> m_vertexPositionData;
	std::vector<glm::vec3> m_vertexNormalData;
	std::vector<glm::vec2> m_vertexUVData;
    
    ObjFileDecoder::decode(objFile.c_str(), meshId, m_vertexPositionData, m_vertexNormalData, m_vertexUVData, m_indexData);
	
    m_vertexData.reserve(m_vertexPositionData.size());

    for (unsigned int i = 0; i < m_vertexPositionData.size(); i++) {
        Vertex vertex = {
            m_vertexPositionData[i],
            m_vertexNormalData[i],
            m_vertexUVData[i]
        };
        m_vertexData.push_back(vertex);
    }
    
    batchInfo.numIndices = m_indexData.size();
#else
	ObjFileDecoder::decode(objFile.c_str(), meshId, m_vertexPositionData, m_vertexNormalData, m_vertexUVData);
	batchInfo.numIndices = m_vertexPositionData.size();
#endif

	batchInfo.startIndex = 0;
	s_meshInfoMap[meshId] = this;
}

Mesh::Mesh(const Mesh& mesh, const std::vector<unsigned int>& triList) {
    std::unordered_map<Vertex, unsigned int, VertexHash> vertexMap;
    for (unsigned int index : triList) {
        if (index >= mesh.m_vertexData.size()) {
            throw std::out_of_range("Index in mesh.m_indexData is out of range of m_vertexData.");
        }
        const Vertex& vertex = mesh.m_vertexData[index];

        if (vertexMap.find(vertex) == vertexMap.end()) {
            vertexMap[vertex] = m_vertexData.size();
            m_vertexData.push_back(vertex);
        }

        m_indexData.push_back(vertexMap[vertex]);
    }

    m_vertexData.shrink_to_fit();
    m_indexData.shrink_to_fit();
}

Mesh::Mesh(const std::vector<Mesh>& mergeMeshes) {
    std::unordered_map<Vertex, unsigned int, VertexHash> vertexMap;

    for (const Mesh& mesh : mergeMeshes) {
        for (const Vertex& vertex : mesh.m_vertexData) {
            if (vertexMap.find(vertex) == vertexMap.end()) {
                vertexMap[vertex] = m_vertexData.size();
                m_vertexData.push_back(vertex);
            }
        }

        for (unsigned int index : mesh.m_indexData) {
            if (index >= mesh.m_vertexData.size()) {
                throw std::out_of_range("Index in mesh.m_indexData is out of range of m_vertexData.");
            }
            m_indexData.push_back(vertexMap[mesh.m_vertexData[index]]);
        }
    }

    m_vertexData.shrink_to_fit();
    m_indexData.shrink_to_fit();
}

void Mesh::uploadToGPU() {
    if(m_clusterList.size() > 0) {
        for(Mesh* cluster : m_clusterList) {
            cluster->uploadToGPU();
        }
        return;
    }

    

    glGenVertexArrays(1, &m_vao);
    glBindVertexArray(m_vao);

    glGenBuffers(1, &m_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, m_vbo);
    glBufferData(GL_ARRAY_BUFFER, m_vertexData.size() * sizeof(Vertex), m_vertexData.data(), GL_STATIC_DRAW);

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

void Mesh::draw(const ShaderProgram& shader) const {
    if(m_clusterList.size() > 0) {
        for (Mesh* cluster : m_clusterList) {
            cluster->draw(shader);
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

