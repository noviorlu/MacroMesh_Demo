premake5 vs2022

msbuild vsproject\CS488_Framework.sln /p:Configuration=Debug /p:Platform=x64
msbuild vsproject\CS488_Framework.sln /p:Configuration=Release /p:Platform=x64

msbuild BuildStaticLibs.sln /property:Configuration=Debug /property:Platform=x64 /t:Rebuild
msbuild BuildStaticLibs.sln /property:Configuration=Release /property:Platform=x64 /t:Rebuild

msbuild BuildStaticLibs.sln /property:Configuration=Debug /property:Platform=x64 /t:Build

cd clusterLOD
premake5 vs2022
msbuild vsproject\CS488-Projects.sln /property:Configuration=Debug /property:Platform=x64 /t:Rebuild
msbuild vsproject\CS488-Projects.sln /property:Configuration=Release /property:Platform=x64 /t:Rebuild



msbuild vsproject\CS488-Projects.sln /property:Configuration=Release /property:Platform=x64 /t:Build

msbuild vsproject\CS488-Projects.sln /property:Configuration=Debug /property:Platform=x64 /t:Build
.\bin\cluster\cluster.exe .\Assets\simpleScene.lua

cd ..
msbuild vsproject\CS488_Framework.sln /p:Configuration=Release /p:Platform=x64
msbuild BuildStaticLibs.sln /property:Configuration=Release /property:Platform=x64 /t:Rebuild
cd clusterLOD
msbuild vsproject\CS488-Projects.sln /property:Configuration=Release /property:Platform=x64 /t:Rebuild
.\bin\cluster\cluster.exe .\Assets\simpleScene.lua


cd ..
msbuild vsproject\CS488_Framework.sln /p:Configuration=Debug /p:Platform=x64
msbuild BuildStaticLibs.sln /property:Configuration=Debug /property:Platform=x64 /t:Rebuild
cd clusterLOD
msbuild vsproject\CS488-Projects.sln /property:Configuration=Debug /property:Platform=x64 /t:Rebuild
.\bin\cluster\cluster.exe .\Assets\simpleScene.lua

msbuild vsproject\CS488-Projects.sln /property:Configuration=Debug /property:Platform=x64 /t:Build
