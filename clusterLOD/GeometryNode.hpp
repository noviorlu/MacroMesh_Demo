// Termm-Fall 2024

#pragma once

#include "SceneNode.hpp"

class GeometryNode : public SceneNode {
public:
	GeometryNode(
		const std::string & meshId,
		const std::string & name
	);

	void draw(
		const glm::mat4 & modelMatrix,
		const glm::mat4 & viewMatrix,
		const ShaderProgram & shader
	) const override;

	Material material;

	// Mesh Identifier. This must correspond to an object name of
	// a loaded .obj file.
	std::string meshId;
	glm::vec3 selectedColor = glm::vec3(0.0, 0.0, 0.5);
};
