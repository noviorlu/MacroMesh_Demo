// Termm-Fall 2024

#include "GeometryNode.hpp"

#include "cs488-framework/GlErrorCheck.hpp"
#include "cs488-framework/MathUtils.hpp"

#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/io.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <iostream>

using namespace glm;

//---------------------------------------------------------------------------------------
GeometryNode::GeometryNode(
		const std::string & meshId,
		const std::string & name
)
	: SceneNode(name),
	  meshId(meshId)
{
	m_nodeType = NodeType::GeometryNode;
}

void GeometryNode::draw(
	const glm::mat4 & modelMatrix, const glm::mat4 &viewMatrix, 
	const ShaderProgram &shader) const
{
	shader.enable();
	
	//-- Set ModelView matrix:
	mat4 modelView = viewMatrix * modelMatrix * trans;
	auto id = shader.getUniformLocation("ModelView");
	glUniformMatrix4fv(id, 1, GL_FALSE, &modelView[0][0]);

	//-- Set NormMatrix:
	mat3 normalMatrix = glm::transpose(glm::inverse(mat3(modelView)));
	auto id2 = shader.getUniformLocation("NormalMatrix");
	glUniformMatrix3fv(id2, 1, GL_FALSE, &normalMatrix[0][0]);

	//-- Set Material values:
	if(selected) {
		auto id3 = shader.getUniformLocation("material.kd");
		glUniform3fv(id3, 1, &selectedColor[0]);
	}
	else{
		auto id3 = shader.getUniformLocation("material.kd");
		glUniform3fv(id3, 1, &material.kd[0]);
	}
 
	//-- Set SceneNode ID:
	auto id4 = shader.getUniformLocation("nodeId");
	glUniform1i(id4, m_nodeId);

	Mesh::s_meshInfoMap[meshId]->draw();

	shader.disable();
}
