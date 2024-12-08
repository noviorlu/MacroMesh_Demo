-- Include directories
includeDirList = {
    "../shared",
    "../shared/include",
    "../shared/gl3w",
    "../shared/imgui",
    "C:/projects/vcpkg/installed/x64-windows/include" -- Add vcpkg include path
}

-- Library directories
libDirectories = {
    "C:/projects/vcpkg/installed/x64-windows-static/lib",
    "C:/projects/MacroMesh_Demo/lib"
}

-- Libraries to link
linkLibs = {
    "cs488-framework", -- Project framework library
    "imgui", -- ImGui library
    "glfw3", -- GLFW library
    "opengl32", -- OpenGL library
    "gdi32", -- Windows Graphics Device Interface
    "user32", -- Windows User Interface
    "kernel32", -- Windows Kernel API
    "shell32", -- Windows Shell API
    "Imm32", -- Input Method Manager
    "lua", -- Lua library
    "tinyobjloader", -- TinyObjLoader library
    "metis", -- METIS library
    "gklib" -- GKlib for METIS
}

-- Solution definition
workspace "CS488-Projects"
    configurations { "Debug", "Release" }
    platforms { "x64" } -- Target 64-bit architecture
    location "vsproject" -- Visual Studio project files location

    filter "platforms:x64"
        architecture "x64"

-- Cluster project
project "cluster"
    kind "ConsoleApp"
    language "C++"
    cppdialect "C++17"
    location "build/cluster"
    objdir "build/cluster/obj"
    targetdir "bin/cluster"

    includedirs (includeDirList)
    libdirs (libDirectories)
    links (linkLibs)

    files { "*.cpp" }

    filter "configurations:Debug"
        defines { "DEBUG" }
        symbols "On"

    filter "configurations:Release"
        defines { "NDEBUG" }
        optimize "On"

