@echo off

rem ***************************************************************************
rem Compile Win32 libraries
rem ***************************************************************************

call .\bin\buildlib_win.cmd dbg Win32
call .\bin\buildlib_win.cmd opt Win32


rem ***************************************************************************
rem Clean up temporary files
rem ***************************************************************************

call .\bin\cleanobj_win.cmd
