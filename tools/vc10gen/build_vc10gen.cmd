@echo off

set DEVENV=""

rem Test for Visual Studio 2010
if "%VS100COMNTOOLS%" neq "" set DEVENV="%VS100COMNTOOLS%..\IDE\devenv.exe"

rem Ensure that we have a path to devenv.exe
if %DEVENV% == "" goto novs

rem ===========================================================================

echo Compiling 'vc10gen'...

%DEVENV% vc10gen.sln /build "Debug|Win32"
if %errorlevel% == 0 goto end

echo.
echo ERROR: Failed to compile 'vc10gen' !
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