/*
 * ShaderProgram
 */

#pragma once

#include "OpenGLImport.hpp"

#include <glm/glm.hpp>
#include <string>
#include <unordered_map>

class ShaderProgram {
public:
    ShaderProgram();

    ~ShaderProgram();

    void generateProgramObject();

    void attachVertexShader(const char * filePath);
    
    void attachFragmentShader(const char * filePath);
    
    void attachGeometryShader(const char * filePath);

    void link();

    void enable() const;

    void disable() const;

    void recompileShaders();

    GLuint getProgramObject() const;

    // external / internal (Error/Warning)
    GLint getUniformLocation(const std::string& uniformName) const;
    GLint GetUniformLocation(const std::string& uniformName) const;

    void SetUniform1i(const std::string& name, int value) const;
    void SetUniform1iv(const std::string& name, int count, int* value) const;
    void SetUniformMat4f(const std::string& name, const glm::mat4& matrix) const;
    void SetUniformMat3f(const std::string& name, const glm::mat3& matrix) const;
    void SetUniform4fv(const std::string& name, const glm::vec4& value) const;
    void SetUniform4f(const std::string& name, float v0, float v1, float v2, float v3) const;
    void SetUniform3fv(const std::string& name, const glm::vec3& value) const;
    void SetUniform3f(const std::string& name, float v0, float v1, float v2) const;
    void SetUniform1f(const std::string& name, float value) const;

    GLint getAttribLocation(const char * attributeName) const;
    struct Shader {
        GLuint shaderObject;
        std::string filePath;

        Shader();
    };

    Shader vertexShader;
    Shader fragmentShader;
    Shader geometryShader;

    GLuint programObject;
    GLuint prevProgramObject;
    GLuint activeProgram;

    void extractSourceCode(std::string & shaderSource, const std::string & filePath);
    
    void extractSourceCodeAndCompile(const Shader &shader);

    GLuint createShader(GLenum shaderType);

    void compileShader(GLuint shaderObject, const std::string & shader);

    void checkCompilationStatus(GLuint shaderObject);

    void checkLinkStatus();

    void deleteShaders();

    mutable std::unordered_map<std::string, int> m_UniformLocationCache;
};

