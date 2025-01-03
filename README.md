# Cluster-based Dynamic LOD generation & selection

https://github.com/user-attachments/assets/996cf7c3-2f4e-46a7-ac60-fadd8ceb9999

---

## Dependencies
* OpenGL 3.2+
* GLFW
    * http://www.glfw.org/
* Lua
    * http://www.lua.org/
* Premake4
    * https://github.com/premake/premake-4.x/wiki
    * http://premake.github.io/download.html
* GLM
    * http://glm.g-truc.net/0.9.7/index.html
* ImGui
    * https://github.com/ocornut/imgui
* METIS

---


## Building Projects
We use **premake4** as our cross-platform build system. First you will need to build all
the static libraries that the projects depend on. To build the libraries, open up a
terminal, and **cd** to the top level of the CS488 project directory and then run the
following:

    $ premake4 gmake
    $ make

This will build the following static libraries, and place them in the top level **lib**
folder of your cs488 project directory.
* libcs488-framework.a
* libglfw3.a
* libimgui.a

Next we can build a specific project.  To do this, **cd** into one of the project folders,
say **A0** for example, and run the following terminal commands in order to compile the A0 executable using all .cpp files in the A0 directory:

    $ cd A0/
    $ premake4 gmake
    $ make


----

## Windows
Sorry for all of the hardcore Microsoft fans out there.  We have not had time to test the build system on Windows yet. Currently our build steps work for OSX and Linux, but all the steps should be the same on Windows, except you will need different libraries to link against when building your project executables.  Some good news is that premake4 can output a Visual Studio .sln file by running:

    $ premake4 gmake

 This should point you in the general direction.

 if imgui has some fault like,
 
    $ process_begin: CreateProcess(NULL, cc -MD -MP -DDEBUG -I../shared -I../shared/gl3w -I../shared/imgui -I../shared/include -I../shared/imgui/examples/opengl3_example -I../shared/imgui/examples/libs/gl3w -g -o Debug/imgui/gl3w.o -MF Debug/imgui/gl3w.d -c ../shared/gl3w/GL/gl3w.c, ...) failed.
    $ make (e=2): The system cannot find the file specified.
    $ make[1]: *** [imgui.make:133: Debug/imgui/gl3w.o] Error 2
    $ make: *** [Makefile:37: imgui] Error 2

goto file **build/imgui.make**, change **$(CC)** into **gcc**

I strongly Recemmond using vcpkg to install all necessary packages.
vcpkg install egl-registry gklib glew glfw3 glm lua metis opengl-registry opengl stb tinyobjloader vcpkg-cmake-config vcpkg-cmake
