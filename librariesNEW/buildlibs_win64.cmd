@echo off

rem ***************************************************************************
rem Compile Win64 libraries
rem ***************************************************************************

call .\bin\buildlib_win.cmd Debug Win64
call .\bin\buildlib_win.cmd Release Win64


rem ***************************************************************************
rem Clean up temporary files
rem ***************************************************************************

call .\bin\cleanobj_win.cmd
