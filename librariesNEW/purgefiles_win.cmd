@echo off

echo Purging AMD...
call .\AMD\purgefiles_win.cmd
echo Purging BLAS...
call .\BLAS\purgefiles_win.cmd
echo Purging lapack-3.4.2...
call .\lapack-3.4.2\purgefiles_win.cmd
echo Purging SuiteSparse_config
call .\SuiteSparse_config\purgefiles_win.cmd
echo Purging UMFPACK
call .\UMFPACK\purgefiles_win.cmd
