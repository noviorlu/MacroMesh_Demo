@echo off
set CONFIG=Debug
set PLATFORM=x64
set OPTIMIZATION=/O2
set BUILD=Build

msbuild vsproject\CS488_Framework.sln /p:Configuration=%CONFIG% /p:Platform=%PLATFORM% /p:CL=%OPTIMIZATION% /t:%BUILD%
msbuild BuildStaticLibs.sln /p:Configuration=%CONFIG% /p:Platform=%PLATFORM% /p:CL=%OPTIMIZATION% /t:%BUILD%

cd clusterLOD
msbuild vsproject\CS488-Projects.sln /p:Configuration=%CONFIG% /p:Platform=%PLATFORM% /p:CL=%OPTIMIZATION% /t:%BUILD%

@REM .\bin\cluster\cluster.exe .\Assets\simpleScene.lua

