@echo off

rem ***************************************************************************
rem Compile Win64 libraries
rem ***************************************************************************

call .\visual_studio\buildlib_win.cmd dbg x64
call .\visual_studio\buildlib_win.cmd opt x64


rem ***************************************************************************
rem Clean up temporary files
rem ***************************************************************************

call .\visual_studio\cleanobj_win.cmd
