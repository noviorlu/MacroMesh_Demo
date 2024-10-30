#include "Mesh.hpp"
using namespace glm;
using namespace std;

#include "cs488-framework/Exception.hpp"
#include "cs488-framework/ObjFileDecoder.hpp"

#include "cs488-framework/OpenGLImport.hpp"
#include "cs488-framework/GlErrorCheck.hpp"

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

//----------------------------------------------------------------------------------------
// Returns the starting memory location for vertex position data.
const float * Mesh::getVertexPositionDataPtr() const {
	return &(m_vertexPositionData[0].x);
}

//----------------------------------------------------------------------------------------
// Returns the starting memory location for vertex normal data.
const float * Mesh::getVertexNormalDataPtr() const {
    return &(m_vertexNormalData[0].x);
}

//----------------------------------------------------------------------------------------
// Returns the total number of bytes of all vertex position data.
size_t Mesh::getNumVertexPositionBytes() const {
	return m_vertexPositionData.size() * sizeof(vec3);
}

//----------------------------------------------------------------------------------------
// Returns the total number of bytes of all vertex normal data.
size_t Mesh::getNumVertexNormalBytes() const {
	return m_vertexNormalData.size() * sizeof(vec3);
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


void Mesh::draw() const {
    glBindVertexArray(m_vao);
    
    #if ENABLE_IBO == true
        glDrawElements(GL_TRIANGLES, m_indexData.size(), GL_UNSIGNED_INT, 0);
    #else
        glDrawArrays(GL_TRIANGLES, 0, m_vertexPositionData.size());
    #endif
    
    glBindVertexArray(0);
}