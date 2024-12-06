@echo off
set CONFIG=Debug
set PLATFORM=x64
set OPTIMIZATION=/O2
set BUILD=Build

msbuild vsproject\CS488_Framework.sln /p:Configuration=%CONFIG% /p:Platform=%PLATFORM% /p:CL=%OPTIMIZATION%
msbuild BuildStaticLibs.sln /p:Configuration=%CONFIG% /p:Platform=%PLATFORM% /p:CL=%OPTIMIZATION% /t:%BUILD%

cd clusterLOD
msbuild vsproject\CS488-Projects.sln /p:Configuration=%CONFIG% /p:Platform=%PLATFORM% /p:CL=%OPTIMIZATION% /t:%BUILD%

.\bin\cluster\cluster.exe .\Assets\simpleScene.lua

@REM @echo off
@REM set CONFIG=Release
@REM set PLATFORM=x64
@REM set OPTIMIZATION=/O3

@REM msbuild vsproject\CS488_Framework.sln /p:Configuration=%CONFIG% /p:Platform=%PLATFORM% /p:CL=%OPTIMIZATION%
@REM msbuild BuildStaticLibs.sln /p:Configuration=%CONFIG% /p:Platform=%PLATFORM% /p:CL=%OPTIMIZATION% /t:%BUILD%

@REM cd clusterLOD
@REM msbuild vsproject\CS488-Projects.sln /p:Configuration=%CONFIG% /p:Platform=%PLATFORM% /p:CL=%OPTIMIZATION% /t:%BUILD%

@REM .\bin\cluster\cluster.exe .\Assets\simpleScene.lua
