// Termm-Fall 2024

#pragma once

#include "cs488-framework/CS488Window.hpp"
#include "cs488-framework/OpenGLImport.hpp"
#include "cs488-framework/ShaderProgram.hpp"
#include "Mesh.hpp"

#include "SceneNode.hpp"
#include "JointNode.hpp"
#include "GBuffer.hpp"

#include <glm/glm.hpp>
#include <memory>
#include <vector>
#include <unordered_set>
#include <cmath>

using namespace std;

struct LightSource {
	glm::vec3 position;
	glm::vec3 rgbIntensity;
};


class clusterLOD : public CS488Window {
public:
	clusterLOD(const std::string & luaSceneFile);
	virtual ~clusterLOD();

protected:
	virtual void init() override;
	virtual void appLogic() override;
	virtual void guiLogic() override;
	virtual void draw() override;
	virtual void cleanup() override;

	//-- Virtual callback methods
	virtual bool cursorEnterWindowEvent(int entered) override;
	virtual bool mouseMoveEvent(double xPos, double yPos) override;
	virtual bool mouseButtonInputEvent(int button, int actions, int mods) override;
	virtual bool mouseScrollEvent(double xOffSet, double yOffSet) override;
	virtual bool windowResizeEvent(int width, int height) override;
	virtual bool keyInputEvent(int key, int action, int mods) override;

	//-- One time initialization methods:
	void processLuaSceneFile(const std::string & filename);
	void createShaderProgram();
	void enableVertexShaderInputSlots();
	void uploadVertexDataToVbos();
	void mapVboDataToVertexShaderInputLocations();
	void initViewMatrix();
	void initLightSources();
	void OneLightSource();
	void addRandomLightSource();
	void addTenRandomLightSource();
	void removeLightSource(int i);
	void removeTenLightSource();
	void dynamicLightSource();
	bool isDynamicLightSource = false;

	void initPerspectiveMatrix();
	void uploadCommonSceneUniforms();
	void renderSceneGraph(const SceneNode &node);
	void renderArcCircle();

	glm::mat4 m_perpsective;
	glm::mat4 m_view;

	std::vector<LightSource> m_lights;


	//-- GL resources for mesh geometry data:
	GLuint m_vao_meshData;
	GLuint m_vbo_vertexPositions;
	GLuint m_vbo_vertexNormals;
	GLint m_positionAttribLocation;
	GLint m_normalAttribLocation;

	//-- GL resources for trackball circle geometry:
	GLuint m_vbo_arcCircle;
	GLuint m_vao_arcCircle;
	GLint m_arc_positionAttribLocation;
	ShaderProgram m_shader_arcCircle;

	GBuffer m_gBuffer;
	ShaderProgram m_geometryPass;
	ShaderProgram m_lightingPass;

	LodRuntimeMesh* m_meshConsolidator;
	float m_errorThreshold = 0.0f;

	std::string m_luaSceneFile;

	std::shared_ptr<SceneNode> m_rootNode;

	// Control Stuffs----------------------------------------------------------------------------------------
	bool m_enableCircle;
	bool m_enableZBuffer;
	bool m_enableBackFaceCull;
	bool m_enableFrontFaceCull;

	void resetControls();

	enum ControlMode {
		ORIENTATION,
		JOINTS
	};

	ControlMode m_controlMode;
	double prevXPos, prevYPos;
	bool m_mouseLeftPressed;
	bool m_mouseMiddlePressed;
	bool m_mouseRightPressed;
	bool m_insideCircle;
	void reset();


	std::unordered_set<JointNode*> m_jointNodes;
	SceneNode* selectedNode = nullptr;
	int pickedID = 0;
	bool selectedPick = false;
};
