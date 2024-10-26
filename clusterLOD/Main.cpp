// Term-Fall 2024

#include "clusterLOD.hpp"

#include <iostream>
using namespace std;

int main( int argc, char **argv ) 
{
	if (argc > 1) {
		std::string luaSceneFile(argv[1]);
		std::string title("F24 Assignment 3 * (");
		title += luaSceneFile;
		title += ")";

		CS488Window::launch(argc, argv, new clusterLOD(luaSceneFile), 1024, 768, title);

	} else {
		cout << "Must supply Lua file as First argument to program.\n";
        cout << "For example:\n";
        cout << "./clusterLOD Assets/simpleScene.lua\n";
	}

	return 0;
}
