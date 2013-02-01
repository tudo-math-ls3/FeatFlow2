@echo off

rem Build script for Visual Studio
rem ------------------------------
rem This file is a batch script that will compile one or all library/-ies
rem using the Visual Studio development environment.
rem
rem Script Parameters:
rem ------------------
rem  %1: build configuration; either "Debug" or "Release"
rem  %2: target platform; either "Win32" or "Win64"
rem  %3: optional: project name; if not given all projects are compiled

rem ===========================================================================

set DEVENV=""

rem Test for Visual Studio 2010
if "%VS100COMNTOOLS%" neq "" set DEVENV="%VS100COMNTOOLS%..\IDE\devenv.exe"
if "%VS110COMNTOOLS%" neq "" set DEVENV="%VS110COMNTOOLS%..\IDE\devenv.exe"
if "%VS120COMNTOOLS%" neq "" set DEVENV="%VS120COMNTOOLS%..\IDE\devenv.exe"
if "%VS130COMNTOOLS%" neq "" set DEVENV="%VS130COMNTOOLS%..\IDE\devenv.exe"
if "%VS140COMNTOOLS%" neq "" set DEVENV="%VS140COMNTOOLS%..\IDE\devenv.exe"
if "%VS150COMNTOOLS%" neq "" set DEVENV="%VS150COMNTOOLS%..\IDE\devenv.exe"

if "%VS200COMNTOOLS%" neq "" set DEVENV="%VS200COMNTOOLS%..\IDE\devenv.exe"
if "%VS210COMNTOOLS%" neq "" set DEVENV="%VS210COMNTOOLS%..\IDE\devenv.exe"

rem Ensure that we have a path to devenv.exe
if %DEVENV% == "" goto novs

rem ===========================================================================

rem Check if the third parameter is given
if "%3" == "" goto libsall


rem Compile a single project from the solution

echo Compiling '%3' in '%1' mode for '%2'...

%DEVENV% .\visual_studio\libraries.sln /build "%1|%2" /project %3
if %errorlevel% == 0 goto end

echo.
echo ERROR: Failed to compile '%3' in '%1' Mode for '%2' !
goto end

rem ===========================================================================
:libsall

rem Compile all projects from the solution

echo Compiling all libraries in '%1' mode for '%2'...

%DEVENV% .\visual_studio\libraries.sln /build "%1|%2"
if %errorlevel% == 0 goto end

echo.
echo ERROR: Failed to compile libraries in '%1' Mode for '%2' !
goto end

rem ===========================================================================
:novs

echo.
echo ERROR: No compatible Visual Studio installation found
echo.
goto end

rem ===========================================================================
:end