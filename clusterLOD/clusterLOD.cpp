// Termm-Fall 2024

#include "clusterLOD.hpp"
#include "scene_lua.hpp"
using namespace std;

#include "cs488-framework/GlErrorCheck.hpp"
#include "cs488-framework/MathUtils.hpp"
#include "GeometryNode.hpp"
#include "JointNode.hpp"

#include "HalfEdgeMesh.hpp"

#include <imgui/imgui.h>

#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/io.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <math.h>

std::string AssetFilePath = "C:/projects/MacroMesh_Demo/clusterLOD/Assets/";
std::string ModelFilePath = "C:/projects/MacroMesh_Demo/models/";

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

using namespace glm;

static bool show_gui = true;

const size_t CIRCLE_PTS = 48;

//----------------------------------------------------------------------------------------
// Constructor
clusterLOD::clusterLOD(const std::string & luaSceneFile)
	: m_luaSceneFile(luaSceneFile),
	  m_positionAttribLocation(0),
	  m_normalAttribLocation(0),
	  m_vao_meshData(0),
	  m_vbo_vertexPositions(0),
	  m_vbo_vertexNormals(0),
	  m_vao_arcCircle(0),
	  m_vbo_arcCircle(0)
{

}

//----------------------------------------------------------------------------------------
// Destructor
clusterLOD::~clusterLOD()
{

}

//----------------------------------------------------------------------------------------
/*
 * Called once, at program start.
 */
void clusterLOD::init()
{
	// Set the background colour.
	glClearColor(0.0, 0.0, 0.0, 1.0);

	createShaderProgram();

	glGenVertexArrays(1, &m_vao_arcCircle);
	glGenVertexArrays(1, &m_vao_meshData);
	enableVertexShaderInputSlots();

	processLuaSceneFile(m_luaSceneFile);

	m_meshConsolidator = new Mesh();
	HalfEdgeMesh heMesh;
	heMesh.importMesh(ModelFilePath + "bunny/bunny.obj");
	heMesh.HalfEdgeMeshSplitter();
	Mesh::s_meshInfoMap[heMesh.m_name] = m_meshConsolidator;
	heMesh.exportMesh(m_meshConsolidator->m_clusterList, m_meshConsolidator->m_clusterGroupList);
	heMesh.exportMeshToObjFiles(ModelFilePath + "bunny/LOD1");
	FastQEM(ModelFilePath + "bunny/LOD1");

	// Acquire the MeshInfoMap from the Mesh.
	m_meshConsolidator->uploadToGPU();

	uploadVertexDataToVbos();
	mapVboDataToVertexShaderInputLocations();

	initPerspectiveMatrix();

	initViewMatrix();

	m_rootNode->storeInitialTrans();
	reset();

	OneLightSource();

	m_gBuffer.initialize(m_framebufferWidth, m_framebufferHeight);
	// Exiting the current scope calls delete automatically on meshConsolidator freeing
	// all vertex data resources.  This is fine since we already copied this data to
	// VBOs on the GPU.  We have no use for storing vertex data on the CPU side beyond
	// this point.
}

//----------------------------------------------------------------------------------------
void clusterLOD::processLuaSceneFile(const std::string & filename) {
	// This version of the code treats the Lua file as an Asset,
	// so that you'd launch the program with just the filename
	// of a puppet in the Assets/ directory.
	// std::string assetFilePath = (filename.c_str());
	// m_rootNode = std::shared_ptr<SceneNode>(import_lua(assetFilePath));

	// This version of the code treats the main program argument
	// as a straightforward pathname.
	m_rootNode = std::shared_ptr<SceneNode>(import_lua(filename));
	if (!m_rootNode) {
		std::cerr << "Could Not Open " << filename << std::endl;
	}
}

//----------------------------------------------------------------------------------------
void clusterLOD::createShaderProgram()
{
	m_geometryPass.generateProgramObject();
	m_geometryPass.attachVertexShader( (AssetFilePath + "Deferred/GeometryPass.vs").c_str() );
	m_geometryPass.attachFragmentShader( (AssetFilePath + "Deferred/GeometryPass.fs").c_str() );
	m_geometryPass.link();

	m_lightingPass.generateProgramObject();
	m_lightingPass.attachVertexShader( (AssetFilePath + "Deferred/LightPass.vs").c_str() );
	m_lightingPass.attachFragmentShader( (AssetFilePath + "Deferred/LightPass.fs").c_str() );
	m_lightingPass.link();

	m_shader_arcCircle.generateProgramObject();
	m_shader_arcCircle.attachVertexShader( (AssetFilePath + "arc_VertexShader.vs").c_str() );
	m_shader_arcCircle.attachFragmentShader( (AssetFilePath + "arc_FragmentShader.fs").c_str() );
	m_shader_arcCircle.link();
}

//----------------------------------------------------------------------------------------
void clusterLOD::enableVertexShaderInputSlots()
{
	//-- Enable input slots for m_vao_arcCircle:
	{
		glBindVertexArray(m_vao_arcCircle);

		// Enable the vertex shader attribute location for "position" when rendering.
		m_arc_positionAttribLocation = m_shader_arcCircle.getAttribLocation("position");
		glEnableVertexAttribArray(m_arc_positionAttribLocation);
	}

	// Restore defaults
	glBindVertexArray(0);
}

//----------------------------------------------------------------------------------------
void clusterLOD::uploadVertexDataToVbos () {
	// Generate VBO to store the trackball circle.
	{
		glGenBuffers( 1, &m_vbo_arcCircle );
		glBindBuffer( GL_ARRAY_BUFFER, m_vbo_arcCircle );

		float *pts = new float[ 2 * CIRCLE_PTS ];
		for( size_t idx = 0; idx < CIRCLE_PTS; ++idx ) {
			float ang = 2.0 * M_PI * float(idx) / CIRCLE_PTS;
			pts[2*idx] = cos( ang );
			pts[2*idx+1] = sin( ang );
		}

		glBufferData(GL_ARRAY_BUFFER, 2*CIRCLE_PTS*sizeof(float), pts, GL_STATIC_DRAW);

		glBindBuffer(GL_ARRAY_BUFFER, 0);
		CHECK_GL_ERRORS;
	}
}

//----------------------------------------------------------------------------------------
void clusterLOD::mapVboDataToVertexShaderInputLocations()
{
	// Bind VAO in order to record the data mapping.
	glBindVertexArray(m_vao_arcCircle);

	// Tell GL how to map data from the vertex buffer "m_vbo_arcCircle" into the
	// "position" vertex attribute location for any bound vertex shader program.
	glBindBuffer(GL_ARRAY_BUFFER, m_vbo_arcCircle);
	glVertexAttribPointer(m_arc_positionAttribLocation, 2, GL_FLOAT, GL_FALSE, 0, nullptr);

	//-- Unbind target, and restore default values:
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindVertexArray(0);

	CHECK_GL_ERRORS;
}

//----------------------------------------------------------------------------------------
void clusterLOD::initPerspectiveMatrix()
{
	float aspect = ((float)m_windowWidth) / m_windowHeight;
	m_perpsective = glm::perspective(degreesToRadians(60.0f), aspect, 0.1f, 100.0f);
}


//----------------------------------------------------------------------------------------
void clusterLOD::initViewMatrix() {
	m_view = glm::lookAt(vec3(0.0f, 0.0f, 0.0f), vec3(0.0f, 0.0f, -1.0f),
			vec3(0.0f, 1.0f, 0.0f));
}

//----------------------------------------------------------------------------------------
void clusterLOD::initLightSources() {
	m_lights.clear();
    //World-space position
    for (int x = -1; x <= 1; x += 2) {
        for (int y = -1; y <= 1; y += 2) {
            for (int z = -1; z <= 1; z += 2) {
                LightSource light;
                light.position = vec3(5.0f * x, 5.0f * y, 5.0f * z);
                light.position = glm::vec3(m_rootNode->get_transform() * glm::vec4(light.position, 1.0f));

                if (x < 0) {
                    light.rgbIntensity = vec3(0.17, 0.64, 0.85);
                } else {
                    light.rgbIntensity = vec3(0.81, 0.57, 0.05);
                }

                m_lights.push_back(light);
            }
        }
    }
}

void clusterLOD::OneLightSource(){
	// setup a simple point light
	m_lights.clear();
	LightSource light;
	light.position = vec3(3.0f, 3.0f, 3.0f);
	light.rgbIntensity = vec3(0.8f); // light
	m_lights.push_back(light);
}

void clusterLOD::addRandomLightSource(){
	// Randomly generate a light source inside -10 to 10 cube and with random color
	if(m_lights.size() >= 256){
		std::cout << "Cannot add more than 256 light sources" << std::endl;
		return;
	}
	
	LightSource light;
	light.position = vec3(rand() % 20 - 10, rand() % 20 - 10, rand() % 20 - 10);
	light.rgbIntensity = vec3((rand() % 100) / 100.0f, (rand() % 100) / 100.0f, (rand() % 100) / 100.0f) * 0.5f;


	m_lights.push_back(light);
}

void clusterLOD::addTenRandomLightSource(){
	for(int i = 0; i < 10; i++){
		addRandomLightSource();
	}
}

void clusterLOD::removeLightSource(int i){
	if(i >= m_lights.size() || i < 0){
		std::cout << "Cannot remove light source at index " << i << std::endl;
		return;
	}
	// Remove the ith light source
	m_lights.erase(m_lights.begin() + i);
}

void clusterLOD::removeTenLightSource(){
	for(int i = 0; i < 10; i++){
		if(m_lights.size() == 0){
			std::cout << "No more light sources to remove" << std::endl;
			return;
		}
		m_lights.pop_back();
	}
}

void clusterLOD::dynamicLightSource(){
	// Move the light source in a circle
	for(int i = 0; i < m_lights.size(); i++){
		// dynamically move the light source in a spherical path with different angles
		float time = glfwGetTime();
		float radius = 10.0f;
		float speed = 1.0f;
		float angleOffset = i * (M_PI / 4); // different angle for each light source
		m_lights[i].position.x = radius * cos(speed * time + angleOffset) * sin(speed * time);
		m_lights[i].position.z = radius * sin(speed * time + angleOffset) * sin(speed * time);
		m_lights[i].position.y = radius * cos(speed * time); // spherical coordinate for y
	}
}

//----------------------------------------------------------------------------------------
void clusterLOD::uploadCommonSceneUniforms() {
	int pick = 0;
	if(m_controlMode == ControlMode::JOINTS) pick = m_rootNode->totalSceneNodes();
	m_geometryPass.enable();
	{
		//-- Set Perpsective matrix uniform for the scene:
		m_geometryPass.SetUniformMat4f("Perspective", m_perpsective);
		CHECK_GL_ERRORS;
	}
	m_geometryPass.disable();

	m_lightingPass.enable();
	{
		m_lightingPass.SetUniformMat4f("View", m_view);

		m_lightingPass.SetUniform1i("gPosition", 0);

		m_lightingPass.SetUniform1i("gNormal", 1);

		m_lightingPass.SetUniform1i("gAlbedoID", 2);

		m_lightingPass.SetUniform1i("numLights", m_lights.size());


		for (int i = 0; i < m_lights.size(); ++i) {
			std::string lightPosStr = "lightPositions[" + std::to_string(i) + "]";
			auto &light = m_lights[i];

			m_lightingPass.SetUniform3fv(lightPosStr.c_str(), light.position);
		}

		for (int i = 0; i < m_lights.size(); ++i) {
			std::string lightColorStr = "lights[" + std::to_string(i) + "].Color";
			auto &light = m_lights[i];

			m_lightingPass.SetUniform3fv(lightColorStr.c_str(), light.rgbIntensity);
		}
	}
	m_lightingPass.disable();
}

//----------------------------------------------------------------------------------------
/*
 * Called once per frame, before guiLogic().
 */
void clusterLOD::appLogic()
{
	// Place per frame, application logic here ...
	uploadCommonSceneUniforms();

	if(pickedID != 0){
		SceneNode* selectedNode = m_rootNode->findJointNodes(pickedID, nullptr);
		if(selectedNode != nullptr){
			// find if selectedNode exists in m_jointNodes, if yes delete, if no push_back
			auto jNode = static_cast<JointNode*>(selectedNode);

			auto it = m_jointNodes.find(jNode);
			if (it != m_jointNodes.end()) {
				m_jointNodes.erase(it);
			} else {
				
				m_jointNodes.insert(jNode);
			}
		}

		pickedID = 0;
	}

	if(isDynamicLightSource){
		dynamicLightSource();	
	}
}

//----------------------------------------------------------------------------------------
/*
 * Called once per frame, after appLogic(), but before the draw() method.
 */
void clusterLOD::guiLogic()
{
	if( !show_gui ) {
		return;
	}

	static bool firstRun(true);
	if (firstRun) {
		ImGui::SetNextWindowPos(ImVec2(50, 50));
		firstRun = false;
	}

	static bool showDebugWindow(true);
	ImGuiWindowFlags windowFlags(ImGuiWindowFlags_AlwaysAutoResize);
	float opacity(0.5f);

	ImGui::Begin("Properties", &showDebugWindow, ImVec2(100,100), opacity,
			windowFlags);

		// Add more gui elements here here ...
		ImGui::Checkbox("Enable Circle", &m_enableCircle);
		ImGui::Checkbox("Enable Z-Buffer", &m_enableZBuffer);
		ImGui::Checkbox("Enable Backface Culling", &m_enableBackFaceCull);
		ImGui::Checkbox("Enable Frontface Culling", &m_enableFrontFaceCull);

		ImGui::RadioButton("Orientation", (int*)&m_controlMode, ControlMode::ORIENTATION);
		ImGui::RadioButton("Joints", (int*)&m_controlMode, ControlMode::JOINTS);

		if(m_controlMode == ControlMode::JOINTS){
			ImGui::Text("Selected Joints: ");
			for(auto joint : m_jointNodes){
				ImGui::Text(joint->m_name.c_str());
			}
		}

		if (ImGui::Button("Reset")) {
			reset();
		}

		// Create Button, and check if it was clicked:
		if( ImGui::Button( "Quit Application" ) ) {
			glfwSetWindowShouldClose(m_window, GL_TRUE);
		}

		if (ImGui::IsAnyItemHovered()) {
			m_mouseLeftPressed = false;
			m_mouseMiddlePressed = false;
			m_mouseRightPressed = false;
		}

		ImGui::Text( "Framerate: %.1f FPS", ImGui::GetIO().Framerate );

	ImGui::End();

	ImGui::Begin("Light Source", &showDebugWindow, ImVec2(100,100), opacity,
			windowFlags);

		if (ImGui::Button("One Light Source")) {
			OneLightSource();
		}

		if (ImGui::Button("Multi Light Source")) {
			initLightSources();
		}


		if (ImGui::Button("Add Light Source")) {
			addRandomLightSource();
		}

		if (ImGui::Button("Add 10 Light Sources")) {
			addTenRandomLightSource();
		}

		if (ImGui::Button("Remove 10 Light Sources")) {
			removeTenLightSource();
		}
		ImGui::Checkbox("Dynamic Light Source", &isDynamicLightSource);
		
		for (int i = 0; i < m_lights.size(); ++i) {
			std::string lightPosStr = "Light " + std::to_string(i);
			ImGui::SliderFloat3(lightPosStr.c_str(), &m_lights[i].position.x, -10.0f, 10.0f);
			ImGui::ColorEdit3(lightPosStr.c_str(), &m_lights[i].rgbIntensity.x);
			if (ImGui::Button("Remove Light Source")) {
				removeLightSource(i);
			}
		}

	ImGui::End();
}

//----------------------------------------------------------------------------------------
/*
 * Called once per frame, after guiLogic().
 */
void clusterLOD::draw() {
	renderSceneGraph(*m_rootNode);

	if(m_enableCircle) {
		renderArcCircle();
	}
}

//----------------------------------------------------------------------------------------
void clusterLOD::renderSceneGraph(const SceneNode & root) {
	// Geometry pass
	m_gBuffer.bind();
	glClearColor(0.0, 0.0, 0.0, 1.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
	// REALLY IMPORTANT!!!!! if blend on, the GBuffer texture automatically unify
	glDisable(GL_BLEND);

	if(m_enableZBuffer) 
		glEnable(GL_DEPTH_TEST);
	else 
		glDisable(GL_DEPTH_TEST);
	
	if(m_enableBackFaceCull || m_enableFrontFaceCull) {
		glEnable(GL_CULL_FACE);
		if(m_enableBackFaceCull && m_enableFrontFaceCull) {
			glCullFace(GL_FRONT_AND_BACK);
		}
		else if(m_enableBackFaceCull) {
			glCullFace(GL_BACK);
		}
		else {
			glCullFace(GL_FRONT);
		}
	}
	else {
		glDisable(GL_CULL_FACE);
	}
	CHECK_GL_ERRORS;

	glBindVertexArray(m_vao_meshData);
	m_geometryPass.enable();
	root.draw(glm::mat4(1.0f), m_view, m_geometryPass);
	m_geometryPass.disable();
	glBindVertexArray(0);
	CHECK_GL_ERRORS;

	// for color picking propose
	if(m_controlMode == ControlMode::JOINTS && selectedPick){
		// get cursor position
		double xpos, ypos;
		glfwGetCursorPos(m_window, &xpos, &ypos);
		int viewport[4];
		glGetIntegerv(GL_VIEWPORT, viewport); 
		int adjustedY = viewport[3] - (int)ypos - 1;

		float pixelColor[4];
		glReadBuffer(GL_COLOR_ATTACHMENT1);
		glReadPixels((int)xpos, adjustedY, 1, 1, GL_RGBA, GL_FLOAT, pixelColor);

		// print pixelColor
		std::cout << "pixelColor: " << pixelColor[0] << " " << pixelColor[1] << " " << pixelColor[2] << " " << pixelColor[3] << std::endl;
		pickedID = (int)pixelColor[3];

		selectedPick = false;
	}

	m_gBuffer.unbind();
	// Lighting pass
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	m_lightingPass.enable();

	glDisable(GL_DEPTH_TEST);
	glDisable(GL_CULL_FACE);

	m_gBuffer.draw(m_lightingPass);
	CHECK_GL_ERRORS;
	m_lightingPass.disable();
}

//----------------------------------------------------------------------------------------
// Draw the trackball circle.
void clusterLOD::renderArcCircle() {
	glBindVertexArray(m_vao_arcCircle);

	m_shader_arcCircle.enable();
		GLint m_location = m_shader_arcCircle.getUniformLocation( "M" );
		float aspect = float(m_framebufferWidth)/float(m_framebufferHeight);
		glm::mat4 M;
		if( aspect > 1.0 ) {
			M = glm::scale( glm::mat4(), glm::vec3( 0.5/aspect, 0.5, 1.0 ) );
		} else {
			M = glm::scale( glm::mat4(), glm::vec3( 0.5, 0.5*aspect, 1.0 ) );
		}
		glUniformMatrix4fv( m_location, 1, GL_FALSE, value_ptr( M ) );
		glDrawArrays( GL_LINE_LOOP, 0, CIRCLE_PTS );
	m_shader_arcCircle.disable();

	glBindVertexArray(0);
	CHECK_GL_ERRORS;
}

//----------------------------------------------------------------------------------------
/*
 * Called once, after program is signaled to terminate.
 */
void clusterLOD::cleanup()
{

}

//----------------------------------------------------------------------------------------
/*
 * Event handler.  Handles cursor entering the window area events.
 */
bool clusterLOD::cursorEnterWindowEvent (
		int entered
) {
	bool eventHandled(false);

	// Fill in with event handling code...

	return eventHandled;
}

//----------------------------------------------------------------------------------------
/*
 * Event handler.  Handles mouse cursor movement events.
 */
bool clusterLOD::mouseMoveEvent (
		double xPos,
		double yPos
) {
	bool eventHandled(false);

	// if mouse on imgui do nothing
	if(ImGui::GetIO().WantCaptureMouse) return eventHandled;

	// Fill in with event handling code...
	double deltaX = xPos - prevXPos;
	double deltaY = yPos - prevYPos;
	
	if( m_controlMode == ControlMode::ORIENTATION ) {
		if(m_mouseLeftPressed){
			vec3 translation(deltaX * 0.01f, -deltaY * 0.01f, 0.0f);
			m_rootNode->translate(translation);
		}
		else if(m_mouseMiddlePressed){
			vec3 translation(0.0f, 0.0f, -deltaY * 0.01f);
			m_rootNode->translate(translation);
		}
		else if(m_mouseRightPressed){
			float angleX = deltaX * 0.5f;
			float angleY = deltaY * 0.5f;

			if(m_insideCircle){
				m_rootNode->rotate('y', angleX);
				m_rootNode->rotate('x', angleY);
			}
			else{
				m_rootNode->rotate('z', angleY);
			}
		}
	}

	if( m_controlMode == ControlMode::JOINTS) {
		if (m_mouseMiddlePressed) {
			for (auto joint : m_jointNodes) {
				float deltaAngle = -deltaY * 0.5f;
				joint->rotate('x', deltaAngle);
			}
		} else if (m_mouseRightPressed) {
			for (auto joint : m_jointNodes) {
				float deltaAngle = deltaX * 0.5f;
				joint->rotate('y', deltaAngle);
			}
		}
	}


	prevXPos = xPos;
	prevYPos = yPos;
	eventHandled = true;
	return eventHandled;
}

//----------------------------------------------------------------------------------------
/*
 * Event handler.  Handles mouse button events.
 */
bool clusterLOD::mouseButtonInputEvent (
		int button,
		int actions,
		int mods
) {
	bool eventHandled(false);

	// Fill in with event handling code...
	if( actions == GLFW_PRESS ) {
		double xpos, ypos;
		glfwGetCursorPos(m_window, &xpos, &ypos);

		prevXPos = xpos;
		prevYPos = ypos;

		if( button == GLFW_MOUSE_BUTTON_LEFT ) {
			m_mouseLeftPressed = true;
			selectedPick = true;
		}
		else if( button == GLFW_MOUSE_BUTTON_RIGHT ) {
			m_mouseRightPressed = true;

			// check if mouse inside circle
			float x = (2.0f * prevXPos) / m_framebufferWidth - 1.0f;
			float y = 1.0f - (2.0f * prevYPos) / m_framebufferHeight;

			float aspect = float(m_framebufferWidth) / float(m_framebufferHeight);
			if (aspect > 1.0f) {
				x *= aspect;
			} else {
				y /= aspect;
			}

			float distance = sqrt(x * x + y * y);
			if (distance <= 0.5f) {
				m_insideCircle = true;
			} else {
				m_insideCircle = false;
			}
		}
		else if( button == GLFW_MOUSE_BUTTON_MIDDLE ) {
			m_mouseMiddlePressed = true;
		}
	}
	else if( actions == GLFW_RELEASE ) {
		if( button == GLFW_MOUSE_BUTTON_LEFT ) {
			m_mouseLeftPressed = false;
		}
		else if( button == GLFW_MOUSE_BUTTON_RIGHT ) {
			m_mouseRightPressed = false;
		}
		else if( button == GLFW_MOUSE_BUTTON_MIDDLE ) {
			m_mouseMiddlePressed = false;
		}
	}

	return eventHandled;
}

//----------------------------------------------------------------------------------------
/*
 * Event handler.  Handles mouse scroll wheel events.
 */
bool clusterLOD::mouseScrollEvent (
		double xOffSet,
		double yOffSet
) {
	bool eventHandled(false);

	// Fill in with event handling code...

	return eventHandled;
}

//----------------------------------------------------------------------------------------
/*
 * Event handler.  Handles window resize events.
 */
bool clusterLOD::windowResizeEvent (
		int width,
		int height
) {
	bool eventHandled(false);
	initPerspectiveMatrix();
	return eventHandled;
}

//----------------------------------------------------------------------------------------
/*
 * Event handler.  Handles key input events.
 */
bool clusterLOD::keyInputEvent (
		int key,
		int action,
		int mods
) {
	bool eventHandled(false);

	if( action == GLFW_PRESS ) {
		if( key == GLFW_KEY_M ) {
			show_gui = !show_gui;
			eventHandled = true;
		}
	}
	// Fill in with event handling code...
	if (action == GLFW_PRESS) {
		switch (key) {
			case GLFW_KEY_A:
				reset();
				eventHandled = true;
				break;
			case GLFW_KEY_Q:
				glfwSetWindowShouldClose(m_window, GL_TRUE);
				eventHandled = true;
				break;
			case GLFW_KEY_C:
				m_enableCircle = !m_enableCircle;
				eventHandled = true;
				break;
			case GLFW_KEY_Z:
				m_enableZBuffer = !m_enableZBuffer;
				eventHandled = true;
				break;
			case GLFW_KEY_B:
				m_enableBackFaceCull = !m_enableBackFaceCull;
				eventHandled = true;
				break;
			case GLFW_KEY_F:
				m_enableFrontFaceCull = !m_enableFrontFaceCull;
				eventHandled = true;
				break;
			case GLFW_KEY_P:
				m_controlMode = ControlMode::ORIENTATION;
				eventHandled = true;
				break;
			case GLFW_KEY_J:
				m_controlMode = ControlMode::JOINTS;
				eventHandled = true;
				break;
		}
	}

	return eventHandled;
}

//----------------------------------------------------------------------------------------
void clusterLOD::resetControls(){
	m_enableCircle = true;
	m_enableZBuffer = true;
	m_enableBackFaceCull = false;
	m_enableFrontFaceCull = false;

	m_controlMode = ControlMode::ORIENTATION;

	m_mouseLeftPressed = false;
	m_mouseMiddlePressed = false;
	m_mouseRightPressed = false;

	for(auto joint : m_jointNodes){
		for(auto child : joint->children){
			child->selected = false;
		}
	}
	m_jointNodes.clear();
}

void clusterLOD::reset(){
	resetControls();
	m_rootNode->restoreInitialTrans();
}