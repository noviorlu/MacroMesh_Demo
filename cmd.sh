premake5 vs2022
msbuild BuildStaticLibs.sln /property:Configuration=Debug /property:Platform=x64 /t:Rebuild

cd clusterLOD
premake5 vs2022
msbuild vsproject\CS488-Projects.sln /property:Configuration=Debug /property:Platform=x64 /t:Rebuild