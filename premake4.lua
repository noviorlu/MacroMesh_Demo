workspace "CS488_Framework"
    configurations { "Debug", "Release" }
    platforms { "x86", "x64" }
    location "vsproject"

    -- 设置架构
    filter "platforms:x86"
        architecture "x86"
    filter "platforms:x64"
        architecture "x64"

-- 包含目录和库路径
includeDirList = { 
    "shared",
    "shared/gl3w",
    "shared/imgui",
    "shared/include"
}

libDirectories = {
    "lib"
}

buildOptions = {"-std=c++17"}

-- 获取当前系统平台
PLATFORM = os.target()

-- 添加 vcpkg 路径（根据需要调整路径）
if PLATFORM == "windows" then
    local vcpkgBasePath = "C:/projects/vcpkg/installed/x64-windows"
    table.insert(includeDirList, vcpkgBasePath .. "/include")
    table.insert(libDirectories, vcpkgBasePath .. "/lib")
end

-- GLFW 静态库生成
if PLATFORM == "windows" and not os.isfile("lib/glfw3.lib") then
    os.execute("cmake -S shared/glfw-3.3.8 -B shared/glfw-3.3.8/build -G \"Visual Studio 16 2019\" -A x64")
    os.execute("cmake --build shared/glfw-3.3.8/build --config Release")
    os.execute("copy shared\\glfw-3.3.8\\build\\src\\Release\\glfw3.lib lib\\glfw3.lib")
end

-- Lua 静态库生成
if PLATFORM == "windows" and not os.isfile("lib/lua.lib") then
    local vcpkgLuaPath = "C:/projects/vcpkg/installed/x64-windows"
    if os.isdir(vcpkgLuaPath) then
        print("Using Lua from vcpkg.")
        os.execute("mkdir lib")
        os.execute("copy " .. vcpkgLuaPath .. "\\lib\\lua.lib lib\\lua.lib")
        table.insert(includeDirList, vcpkgLuaPath .. "\\include")
    else
        os.execute("cmake -S shared/lua-5.4.6 -B shared/lua-5.4.6/build -G \"Visual Studio 16 2019\" -A x64")
        os.execute("cmake --build shared/lua-5.4.6/build --config Release")
        os.execute("copy shared\\lua-5.4.6\\src\\Release\\lua.lib lib\\lua.lib")
    end
end

-- 静态库解决方案
solution "BuildStaticLibs"
    configurations { "Debug", "Release" }
    platforms { "x86", "x64" }

    filter "configurations:Debug"
        runtime "Debug"
        symbols "On"

    filter "configurations:Release"
        runtime "Release"
        optimize "On"

-- cs488-framework 项目
project "cs488-framework"
    kind "StaticLib"
    language "C++"
    location "build/cs488-framework"
    objdir "build/cs488-framework/obj"
    targetdir "lib"
    buildoptions (buildOptions)
    includedirs (includeDirList)
    files { "shared/cs488-framework/*.cpp", "shared/cs488-framework/*.h" }

-- imgui 项目
project "imgui"
    kind "StaticLib"
    language "C++"
    location "build/imgui"
    objdir "build/imgui/obj"
    targetdir "lib"
    includedirs (includeDirList)
    includedirs {
        "shared/imgui/examples/opengl3_example",
        "shared/imgui/examples/libs/gl3w/"
    }
    files { 
        "shared/imgui/*.cpp",
        "shared/gl3w/GL/gl3w.c"
    }

-- lodepng 项目
project "lodepng"
    kind "StaticLib"
    language "C++"
    location "build/lodepng"
    objdir "build/lodepng/obj"
    targetdir "lib"
    includedirs (includeDirList)
    files { "shared/lodepng/lodepng.cpp" }
