// Termm-Fall 2024

#pragma once

#include "SceneNode.hpp"

class JointNode : public SceneNode {
public:
	JointNode(const std::string & name);
	virtual ~JointNode();

	void set_joint_x(double min, double init, double max);
	void set_joint_y(double min, double init, double max);

    SceneNode* findJointNodes(const int& nodeId, SceneNode* closestJointNode) override;

    void rotate(char axis, float angle) override;

	void draw(
		const glm::mat4& modelMatrix, const glm::mat4& viewMatrix, 
		const ShaderProgram& shader) const override;


	struct JointRange {
		double min, init, max;
	};

	void restoreInitialTrans() override { 
        rotate('x', m_joint_x.init - current_x);
		rotate('y', m_joint_y.init - current_y);
		current_x = m_joint_x.init;
		current_y = m_joint_y.init;
		JointNode::rotate('x', current_x);
		JointNode::rotate('y', current_y);
		SceneNode::restoreInitialTrans();
    }
    void storeInitialTrans() override {
		current_x = m_joint_x.init;
		current_y = m_joint_y.init;
		rotate('x', current_x);
		rotate('y', current_y);
        SceneNode::storeInitialTrans();
    }

	JointRange m_joint_x, m_joint_y;
	float current_x = 0, current_y = 0;

	glm::mat4 rot;
};
