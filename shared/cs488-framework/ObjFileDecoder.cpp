#include "ObjFileDecoder.hpp"
using namespace glm;

#include <fstream>
#include <sstream>
#include <iostream>
#include <cstring>
#include <unordered_map>
#include <tuple>
using namespace std;

#include "cs488-framework/Exception.hpp"
#include "cs488-framework/Vertex.hpp"

//---------------------------------------------------------------------------------------
void ObjFileDecoder::decode(
		const char * objFilePath,
		std::string & objectName,
        std::vector<vec3> & positions,
        std::vector<vec3> & normals,
        std::vector<vec2> & uvCoords
) {

	// Empty containers, and start fresh before inserting data from .obj file
	positions.clear();
	normals.clear();
	uvCoords.clear();

    ifstream in(objFilePath, std::ios::in);
    in.exceptions(std::ifstream::badbit);

    if (!in) {
        stringstream errorMessage;
        errorMessage << "Unable to open .obj file " << objFilePath
            << " within method ObjFileDecoder::decode" << endl;

        throw Exception(errorMessage.str().c_str());
    }

    string currentLine;
    int positionIndexA, positionIndexB, positionIndexC;
    int normalIndexA, normalIndexB, normalIndexC;
    int uvCoordIndexA, uvCoordIndexB, uvCoordIndexC;
    vector<vec3> temp_positions;
    vector<vec3> temp_normals;
    vector<vec2> temp_uvCoords;

	objectName = "";

    while (!in.eof()) {
        try {
            getline(in, currentLine);
        } catch (const ifstream::failure &e) {
            in.close();
            stringstream errorMessage;
            errorMessage << "Error calling getline() -- " << e.what() << endl;
            throw Exception(errorMessage.str());
        }
	    if (currentLine.substr(0, 2) == "o ") {
		    // Get entire line excluding first 2 chars.
		    istringstream s(currentLine.substr(2));
		    s >> objectName;


	    } else if (currentLine.substr(0, 2) == "v ") {
            // Vertex data on this line.
            // Get entire line excluding first 2 chars.
            istringstream s(currentLine.substr(2));
            glm::vec3 vertex;
            s >> vertex.x;
            s >> vertex.y;
            s >> vertex.z;
            temp_positions.push_back(vertex);

        } else if (currentLine.substr(0, 3) == "vn ") {
            // Normal data on this line.
            // Get entire line excluding first 2 chars.
            istringstream s(currentLine.substr(2));
            vec3 normal;
            s >> normal.x;
            s >> normal.y;
            s >> normal.z;
            temp_normals.push_back(normal);

        } else if (currentLine.substr(0, 3) == "vt ") {
            // Texture coordinate data on this line.
            // Get entire line excluding first 2 chars.
            istringstream s(currentLine.substr(2));
            vec2 textureCoord;
            s >> textureCoord.s;
            s >> textureCoord.t;
            temp_uvCoords.push_back(textureCoord);

        } else if (currentLine.substr(0, 2) == "f ") {
            // Face index data on this line.

            int index;

            // sscanf will return the number of matched index values it found
            // from the pattern.
            int numberOfIndexMatches = sscanf(currentLine.c_str(), "f %d/%d/%d",
                                              &index, &index, &index);

            if (numberOfIndexMatches == 3) {
                // Line contains indices of the pattern vertex/uv-cord/normal.
                sscanf(currentLine.c_str(), "f %d/%d/%d %d/%d/%d %d/%d/%d",
                       &positionIndexA, &uvCoordIndexA, &normalIndexA,
                       &positionIndexB, &uvCoordIndexB, &normalIndexB,
                       &positionIndexC, &uvCoordIndexC, &normalIndexC);

                // .obj file uses indices that start at 1, so subtract 1 so they start at 0.
                uvCoordIndexA--;
                uvCoordIndexB--;
                uvCoordIndexC--;

                uvCoords.push_back(temp_uvCoords[uvCoordIndexA]);
                uvCoords.push_back(temp_uvCoords[uvCoordIndexB]);
                uvCoords.push_back(temp_uvCoords[uvCoordIndexC]);

            } else {
                // Line contains indices of the pattern vertex//normal.
                sscanf(currentLine.c_str(), "f %d//%d %d//%d %d//%d",
		               &positionIndexA, &normalIndexA,
                       &positionIndexB, &normalIndexB,
                       &positionIndexC, &normalIndexC);
            }

            positionIndexA--;
            positionIndexB--;
            positionIndexC--;
            normalIndexA--;
            normalIndexB--;
            normalIndexC--;

            positions.push_back(temp_positions[positionIndexA]);
            positions.push_back(temp_positions[positionIndexB]);
            positions.push_back(temp_positions[positionIndexC]);

            normals.push_back(temp_normals[normalIndexA]);
            normals.push_back(temp_normals[normalIndexB]);
            normals.push_back(temp_normals[normalIndexC]);
        }
    }

    in.close();

	if (objectName.compare("") == 0) {
		// No 'o' object name tag defined in .obj file, so use the file name
		// minus the '.obj' ending as the objectName.
		const char * ptr = strrchr(objFilePath, '/');
		objectName.assign(ptr+1);
		size_t pos = objectName.find('.');
		objectName.resize(pos);
	}
}

//---------------------------------------------------------------------------------------
void ObjFileDecoder::decode(
		const char * objFilePath,
		std::string & objectName,
        std::vector<vec3> & positions,
        std::vector<vec3> & normals
) {
    std::vector<vec2> uvCoords;
    decode(objFilePath, objectName, positions, normals, uvCoords);
}





void ObjFileDecoder::decode(
    const char* objFilePath,
    std::string& objectName,
    std::vector<Vertex>& vertices,
    std::vector<unsigned int>& indices
) {
    vertices.clear();
    indices.clear();

    std::ifstream in(objFilePath, std::ios::in);
    in.exceptions(std::ifstream::badbit);

    if (!in) {
        throw std::runtime_error("Unable to open .obj file: " + std::string(objFilePath));
    }

    std::string currentLine;
    std::vector<glm::vec3> temp_positions;
    std::vector<glm::vec3> temp_normals;
    std::vector<glm::vec2> temp_uvCoords;

    objectName.clear();

    // Store unique vertices with accumulative data
    std::unordered_map<Vertex, std::pair<unsigned int, float>> uniqueVertexMap;

    auto getOrAddVertexIndex = [&](int posIdx, int uvIdx, int normIdx, float weight) -> unsigned int {
        glm::vec3 position = temp_positions[posIdx];
        glm::vec3 normal = temp_normals[normIdx];
        glm::vec2 uv = temp_uvCoords[uvIdx];

        Vertex vertex(position, normal, uv);
        auto it = uniqueVertexMap.find(vertex);
        if (it != uniqueVertexMap.end()) {
            // Update the existing vertex with weighted attributes
            unsigned int existingIndex = it->second.first;
            float& accumulatedWeight = it->second.second;
            float newWeight = accumulatedWeight + weight;

            // Update the existing vertex attributes with area-weighted averages
            vertices[existingIndex].normal = (vertices[existingIndex].normal * accumulatedWeight + normal * weight) / newWeight;
            vertices[existingIndex].uv = (vertices[existingIndex].uv * accumulatedWeight + uv * weight) / newWeight;

            accumulatedWeight = newWeight;
            return existingIndex;
        }

        // If vertex is unique, add it to the vertex list and map
        unsigned int newIndex = vertices.size();
        vertices.push_back(vertex);
        uniqueVertexMap[vertex] = {newIndex, weight};
        return newIndex;
    };


    while (std::getline(in, currentLine)) {
        if (currentLine.substr(0, 2) == "o ") {
            objectName = currentLine.substr(2);
        } else if (currentLine.substr(0, 2) == "v ") {
            glm::vec3 position;
            sscanf(currentLine.c_str(), "v %f %f %f", &position.x, &position.y, &position.z);
            temp_positions.push_back(position);
        } else if (currentLine.substr(0, 3) == "vn ") {
            glm::vec3 normal;
            sscanf(currentLine.c_str(), "vn %f %f %f", &normal.x, &normal.y, &normal.z);
            temp_normals.push_back(normal);
        } else if (currentLine.substr(0, 3) == "vt ") {
            glm::vec2 uv;
            sscanf(currentLine.c_str(), "vt %f %f", &uv.x, &uv.y);
            temp_uvCoords.push_back(uv);
        } else if (currentLine.substr(0, 2) == "f ") {
            int posIdx[3], uvIdx[3], normIdx[3];
            sscanf(currentLine.c_str(),
                   "f %d/%d/%d %d/%d/%d %d/%d/%d",
                   &posIdx[0], &uvIdx[0], &normIdx[0],
                   &posIdx[1], &uvIdx[1], &normIdx[1],
                   &posIdx[2], &uvIdx[2], &normIdx[2]);

            for (int i = 0; i < 3; i++) {
                posIdx[i]--;
                uvIdx[i]--;
                normIdx[i]--;
            }

            // Calculate the triangle size (area)
            glm::vec3 edge1 = temp_positions[posIdx[1]] - temp_positions[posIdx[0]];
            glm::vec3 edge2 = temp_positions[posIdx[2]] - temp_positions[posIdx[0]];
            float triangleArea = glm::length(glm::cross(edge1, edge2)) * 0.5f;

            // Pass the triangle area as weight
            indices.push_back(getOrAddVertexIndex(posIdx[0], uvIdx[0], normIdx[0], triangleArea));
            indices.push_back(getOrAddVertexIndex(posIdx[1], uvIdx[1], normIdx[1], triangleArea));
            indices.push_back(getOrAddVertexIndex(posIdx[2], uvIdx[2], normIdx[2], triangleArea));
        }
    }

    in.close();

    if (objectName.empty()) {
        const char* baseName = strrchr(objFilePath, '/');
        objectName = baseName ? baseName + 1 : objFilePath;
        objectName = objectName.substr(0, objectName.find('.'));
    }
}

