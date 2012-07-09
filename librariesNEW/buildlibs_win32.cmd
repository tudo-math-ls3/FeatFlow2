@echo off

rem ***************************************************************************
rem Compile Win32 libraries
rem ***************************************************************************

call .\bin\buildlib_win.cmd Debug Win32
call .\bin\buildlib_win.cmd Release Win32


rem ***************************************************************************
rem Clean up temporary files
rem ***************************************************************************

call .\bin\cleanobj_win.cmd
