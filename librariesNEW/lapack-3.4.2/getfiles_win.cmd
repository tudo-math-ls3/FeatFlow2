@echo off

rem Download and extract the file

set URL=http://www.netlib.org/lapack/lapack-3.4.2.tgz

call .\bin\getfile_win.cmd %URL%
