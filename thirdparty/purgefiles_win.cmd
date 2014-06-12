@echo off

rem ============================================================================
echo Purging AMD...
rmdir /s /q .\AMD\Demo
rmdir /s /q .\AMD\Doc
rmdir /s /q .\AMD\Include
rmdir /s /q .\AMD\Lib
rmdir /s /q .\AMD\MATLAB
rmdir /s /q .\AMD\Source
del /q .\AMD\Makefile
del /q .\AMD\README.txt

rem ============================================================================
echo Purging BLAS...
del /q .\BLAS\*.f
del /q .\BLAS\make.inc
del /q .\BLAS\Makefile

rem ============================================================================
echo Purging lapack-3.5.0...
rmdir /s /q .\lapack-3.5.0\BLAS
rmdir /s /q .\lapack-3.5.0\CMAKE
rmdir /s /q .\lapack-3.5.0\DOCS
rmdir /s /q .\lapack-3.5.0\INSTALL
rmdir /s /q .\lapack-3.5.0\lapacke
rmdir /s /q .\lapack-3.5.0\SRC
rmdir /s /q .\lapack-3.5.0\TESTING
del /q .\lapack-3.5.0\CMakeLists.txt
del /q .\lapack-3.5.0\CTestConfig.cmake
del /q .\lapack-3.5.0\CTestCustom.cmake.in
del /q .\lapack-3.5.0\lapack.pc.in
del /q .\lapack-3.5.0\lapack_build.cmake
del /q .\lapack-3.5.0\lapack_testing.py
del /q .\lapack-3.5.0\LICENSE
del /q .\lapack-3.5.0\make.inc.example
del /q .\lapack-3.5.0\Makefile
del /q .\lapack-3.5.0\README

rem ============================================================================
echo Purging SuiteSparse_config...
rmdir /s /q .\SuiteSparse_config\xerbla
del /q .\SuiteSparse_config\Makefile
del /q .\SuiteSparse_config\README.txt
del /q .\SuiteSparse_config\SuiteSparse_config.c
del /q .\SuiteSparse_config\SuiteSparse_config.h
del /q .\SuiteSparse_config\SuiteSparse_config.mk
del /q .\SuiteSparse_config\SuiteSparse_config_GPU.mk
del /q .\SuiteSparse_config\SuiteSparse_config_Mac.mk

rem ============================================================================
echo Purging UMFPACK...
rmdir /s /q .\UMFPACK\Demo
rmdir /s /q .\UMFPACK\Doc
rmdir /s /q .\UMFPACK\Include
rmdir /s /q .\UMFPACK\Lib
rmdir /s /q .\UMFPACK\MATLAB
rmdir /s /q .\UMFPACK\Source
rmdir /s /q .\UMFPACK\Tcov
del /q .\UMFPACK\Makefile
del /q .\UMFPACK\README.txt
