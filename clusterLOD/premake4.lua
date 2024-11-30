-- Shared configurations
includeDirList = {
    "../shared",
    "../shared/include",
    "../shared/gl3w",
    "../shared/imgui",
    "../shared/glfw-3.3.8/include",
}

libDirectories = { 
    "../lib",
    "../shared/glfw-3.3.8/build/src",
}

if os.get() == "macosx" then
    linkLibs = {
        "cs488-framework",
        "imgui",
        "glfw3",
        "lua"
    }
    linkOptionList = { "-framework IOKit", "-framework Cocoa", "-framework CoreVideo", "-framework OpenGL" }
end

if os.get() == "linux" then
    linkLibs = {
        "cs488-framework",
        "imgui",
        "glfw3",
        "lua",
        "GL",
        "Xinerama",
        "Xcursor",
        "Xxf86vm",
        "Xi",
        "Xrandr",
        "X11",
        "stdc++",
        "dl",
        "pthread"
    }
    linkOptionList = {}
end

if os.get() == "windows" then
    linkLibs = {
        "cs488-framework",
        "imgui",
        "glfw3",
        "opengl32",
        "gdi32",
        "user32",
        "kernel32",
        "shell32",
        "Imm32",
        "lua",
        "metis",
        "gklib",
        "mingwex",
        "msvcrt"
    }
    linkoptions { "-lmingwex", "-lmsvcrt" }
    table.insert(includeDirList, "C:/projects/vcpkg/installed/x64-windows/include")
    table.insert(includeDirList, "C:/projects/vcpkg/installed/x86-mingw-static/include")
    table.insert(libDirectories, "C:/projects/vcpkg/installed/x64-windows/lib")
    table.insert(libDirectories, "C:/projects/vcpkg/installed/x86-mingw-static/lib")
    buildoptions { "-std=c++17" }
    toolset "gcc"
    linkOptionList = {}
end

buildOptions = {"-std=c++17"}

-- Handle custom build flags
newoption {
    trigger = "build-cluster",
    description = "Build only the cluster project"
}

newoption {
    trigger = "build-qemtest",
    description = "Build only the QEMTest project"
}

-- Solution definition
solution "CS488-Projects"
    configurations { "Debug", "Release" }

    -- Cluster project
    if not _OPTIONS["build-qemtest"] then
        project "cluster"
            kind "ConsoleApp"
            language "C++"
            cppdialect "C++17"
            location "build/cluster"
            objdir "build/cluster/obj"
            targetdir "bin/cluster"
            buildoptions (buildOptions)
            libdirs (libDirectories)
            links (linkLibs)
            linkoptions (linkOptionList)
            includedirs (includeDirList)
            files { "*.cpp" }

            configuration "Debug"
                defines { "DEBUG" }
                flags { "Symbols" }

            configuration "Release"
                defines { "NDEBUG" }
                flags { "Optimize" }
    end

    -- QEMTest project
    if not _OPTIONS["build-cluster"] then
        project "QEMTest"
            kind "ConsoleApp"
            language "C++"
            cppdialect "C++17"
            location "build/QEMTest"
            objdir "build/QEMTest/obj"
            targetdir "bin/QEMTest"
            buildoptions (buildOptions)
            libdirs (libDirectories)
            links (linkLibs)
            linkoptions (linkOptionList)
            includedirs (includeDirList)
            files { "./QEMTest/QEMTest.cpp", "./HalfEdgeMesh.cpp", "./MeshProcessing.cpp", "./QEM.cpp", "./Mesh.cpp", "./Cluster.cpp" }

            configuration "Debug"
                defines { "DEBUG" }
                flags { "Symbols" }

            configuration "Release"
                defines { "NDEBUG" }
                flags { "Optimize" }
    end
