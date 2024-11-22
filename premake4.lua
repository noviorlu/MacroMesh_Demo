includeDirList = includeDirList or { 
    "shared",
    "shared/gl3w",
    "shared/imgui",
    "shared/include"
}

libDirectories = libDirectories or {}

buildOptions = {"-std=c++17"}

-- Get the current OS platform
PLATFORM = os.target()  -- 使用 os.target() 替代 os.get()

-- Add vcpkg paths for Lua
if PLATFORM == "windows" then
    local vcpkgBasePath = "C:/projects/vcpkg/installed/x86-mingw-static"
    table.insert(includeDirList, vcpkgBasePath .. "/include")
    table.insert(libDirectories, vcpkgBasePath .. "/lib")
end

-- Build glfw3 static library and copy it into <cs488_root>/lib if it is not
-- already present.
if not os.isfile("lib/libglfw3.a") then
    os.chdir("shared/glfw-3.3.8")
    os.mkdir("build")
    os.chdir("build")
    os.execute("cmake ../")
    os.execute("make")
    os.chdir("../../../")
    os.mkdir("lib")
    os.execute("cp shared/glfw-3.3.8/build/src/libglfw3.a lib/")
end

if not os.isfile("lib/liblua.a") then
    local vcpkgLuaPath = "C:/projects/vcpkg/installed/x86-mingw-static"
    if PLATFORM == "windows" and os.isdir(vcpkgLuaPath) then
        print("Using Lua from vcpkg.")
        os.execute("mkdir lib")
        os.execute("copy " .. vcpkgLuaPath .. "\\lib\\liblua.a lib\\liblua.a")
        table.insert(includeDirList, vcpkgLuaPath .. "\\include")
    else
        os.chdir("shared/lua-5.4.6")
        local result
        if PLATFORM == "macosx" then
            result = os.execute("make macosx")
        elseif PLATFORM == "linux" then
            result = os.execute("make linux")
        elseif PLATFORM == "windows" then
            result = os.execute("make mingw")
        end

        print("Make result: ", result)

        os.chdir("../../")
        if PLATFORM == "windows" then
            os.execute("copy shared\\lua-5.4.6\\src\\liblua.a lib\\liblua.a")
        else
            os.execute("cp shared/lua-5.4.6/src/liblua.a lib/")
        end
    end
end

solution "BuildStaticLibs"
    configurations { "Debug", "Release" }

    configuration "Debug"
        defines { "DEBUG" }
        flags { "Symbols" }

    configuration "Release"
        defines { "NDEBUG" }
        flags { "Optimize" }

    -- Builds cs488-framework static library
    project "cs488-framework"
        kind "StaticLib"
        language "C++"
        location "build"
        objdir "build"
        targetdir "lib"
        buildoptions (buildOptions)
        includedirs (includeDirList)
        files { "shared/cs488-framework/*.cpp" }

    -- Build imgui static library
    project "imgui"
        kind "StaticLib"
        language "C++"
        location "build"
        objdir "build"
        targetdir "lib"
        includedirs (includeDirList)
        includedirs {
            "shared/imgui/examples/opengl3_example",
            "shared/imgui/examples/libs/gl3w/",
        }
        files { 
            "shared/imgui/*.cpp",
            "shared/gl3w/GL/gl3w.c"
        }

    -- Build lodepng static library
    project "lodepng"
        kind "StaticLib"
        language "C++"
        location "build"
        objdir "build"
        targetdir "lib"
        includedirs (includeDirList)
        includedirs {
            "shared/lodepng"
        }
        files { 
            "shared/lodepng/lodepng.cpp"
        }
