@echo off

rem Download and extract the file

set URL=http://www.cise.ufl.edu/research/sparse/umfpack/UMFPACK-5.6.0.tar.gz
set URL2=http://www.cise.ufl.edu/research/sparse/SuiteSparse_config/SuiteSparse_config-4.0.0.tar.gz

@echo off
call .\bin\getfile_win.cmd %URL%
call .\bin\getfile_win.cmd %URL2%
