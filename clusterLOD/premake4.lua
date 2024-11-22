-- Additional Linux libs: "X11", "Xxf86vm", "Xi", "Xrandr", "stdc++"
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
        "mingwex", -- 必须在末尾
        "msvcrt"   -- 必须在末尾
    }
    linkoptions { "-lmingwex", "-lmsvcrt" }
    table.insert(includeDirList, "C:/projects/vcpkg/installed/x64-windows/include")
    table.insert(includeDirList, "C:/projects/vcpkg/installed/x86-mingw-static/include")
    table.insert(libDirectories, "C:/projects/vcpkg/installed/x64-windows/lib")
    table.insert(libDirectories, "C:/projects/vcpkg/installed/x86-mingw-static/lib")
    -- 指定 GCC 编译器
    buildoptions { "-std=c++17" }
    toolset "gcc" -- 明确指定使用 GCC 工具链
end


-- Build Options:
if os.get() == "macosx" then
    linkOptionList = { "-framework IOKit", "-framework Cocoa", "-framework CoreVideo", "-framework OpenGL" }
end

buildOptions = {"-std=c++17"}

solution "CS488-Projects"
    configurations { "Debug", "Release" }

    project "cluster"
        kind "ConsoleApp"
        language "C++"
        cppdialect "C++17"
        location "build"
        objdir "build"
        targetdir "."
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
