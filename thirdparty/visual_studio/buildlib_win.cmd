@echo off

rem Build script for Visual Studio
rem ------------------------------
rem This file is a batch script that will compile one or all library/-ies
rem using the Visual Studio development environment.
rem
rem Script Parameters:
rem ------------------
rem  %1: build configuration; either "dbg" or "opt"
rem  %2: target platform; either "Win32" or "x64"
rem  %3: optional: project name; if not given all projects are compiled

rem ===========================================================================

set DEVENV=""

rem Test for Visual Studio 2010
if "%VS100COMNTOOLS%" neq "" set DEVENV="%VS100COMNTOOLS%..\IDE\devenv.exe"

rem Ensure that we have a path to devenv.exe
if %DEVENV% == "" goto novs

rem ===========================================================================

rem Check if the third parameter is given
if "%3" == "" goto libsall

rem Compile a single project from the solution

echo Compiling '%3' in '%1' mode for '%2'...

call %DEVENV% .\visual_studio\thirdparty.sln /build "%1|%2" /project %3
if %errorlevel% == 0 goto end

echo.
echo ERROR: Failed to compile '%3' in '%1' Mode for '%2' !
echo.
pause
goto end

rem ===========================================================================
:libsall

rem Compile all projects from the solution

echo Compiling all libraries in '%1' mode for '%2'...

call %DEVENV% .\visual_studio\thirdparty.sln /build "%1|%2"
if %errorlevel% == 0 goto end

echo.
echo ERROR: Failed to compile libraries in '%1' Mode for '%2' !
echo.
pause
goto end

rem ===========================================================================
:novs

echo.
echo ERROR: No compatible Visual Studio installation found
echo.
pause
goto end

rem ===========================================================================
:end