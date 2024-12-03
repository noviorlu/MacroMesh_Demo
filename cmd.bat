@REM @echo off
@REM set CONFIG=Release
@REM set PLATFORM=x64
@REM set OPTIMIZATION=/O2

@REM msbuild vsproject\CS488_Framework.sln /p:Configuration=%CONFIG% /p:Platform=%PLATFORM% /p:CL=%OPTIMIZATION%
@REM msbuild BuildStaticLibs.sln /p:Configuration=%CONFIG% /p:Platform=%PLATFORM% /p:CL=%OPTIMIZATION% /t:Rebuild

@REM cd clusterLOD
@REM msbuild vsproject\CS488-Projects.sln /p:Configuration=%CONFIG% /p:Platform=%PLATFORM% /p:CL=%OPTIMIZATION% /t:Rebuild

@REM .\bin\cluster\cluster.exe .\Assets\simpleScene.lua

@echo off
set CONFIG=Release
set PLATFORM=x64
set OPTIMIZATION=/O3

msbuild vsproject\CS488_Framework.sln /p:Configuration=%CONFIG% /p:Platform=%PLATFORM% /p:CL=%OPTIMIZATION%
msbuild BuildStaticLibs.sln /p:Configuration=%CONFIG% /p:Platform=%PLATFORM% /p:CL=%OPTIMIZATION% /t:Build

cd clusterLOD
msbuild vsproject\CS488-Projects.sln /p:Configuration=%CONFIG% /p:Platform=%PLATFORM% /p:CL=%OPTIMIZATION% /t:Build

.\bin\cluster\cluster.exe .\Assets\simpleScene.lua
