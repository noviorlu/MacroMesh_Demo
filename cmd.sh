premake5 vs2022
msbuild BuildStaticLibs.sln /property:Configuration=Debug /property:Platform=x64 /t:Rebuild
msbuild BuildStaticLibs.sln /property:Configuration=Debug /property:Platform=x64 /t:Build

cd clusterLOD
premake5 vs2022
msbuild vsproject\CS488-Projects.sln /property:Configuration=Debug /property:Platform=x64 /t:Rebuild

msbuild vsproject\CS488-Projects.sln /property:Configuration=Debug /property:Platform=x64 /t:Build
.\bin\cluster\cluster.exe .\Assets\simpleScene.lua