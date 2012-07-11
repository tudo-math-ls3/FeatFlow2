@echo off

rem Automatic Library downloader
rem
rem This script downloads and extracts all necessary library files
rem from the internet. It calls the corresponding "getfiles_win.bat" scripts
rem in the directories of the libraries. THese subscripts then download
rem the library/-ies and extract them.

call .\AMD\getfiles_win.cmd
call .\BLAS\getfiles_win.cmd
call .\lapack-3.4.1\getfiles_win.cmd
call .\SuiteSparse_config\getfiles_win.cmd
call .\UMFPACK\getfiles_win.cmd
