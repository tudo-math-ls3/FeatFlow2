@echo off

rem ***************************************************************************
rem Compile Win32 libraries
rem ***************************************************************************

call .\visual_studio\buildlib_win.cmd dbg Win32
call .\visual_studio\buildlib_win.cmd opt Win32


rem ***************************************************************************
rem Clean up temporary files
rem ***************************************************************************

call .\visual_studio\cleanobj_win.cmd
