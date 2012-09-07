@echo off

rem Delete all LAPACK files

rmdir /s /q .\lapack-3.4.1\BLAS
rmdir /s /q .\lapack-3.4.1\CMAKE
rmdir /s /q .\lapack-3.4.1\DOCS
rmdir /s /q .\lapack-3.4.1\INSTALL
rmdir /s /q .\lapack-3.4.1\lapacke
rmdir /s /q .\lapack-3.4.1\SRC
rmdir /s /q .\lapack-3.4.1\TESTING

del /q .\lapack-3.4.1\CMakeLists.txt
del /q .\lapack-3.4.1\CTestConfig.cmake
del /q .\lapack-3.4.1\CTestCustom.cmake.in
del /q .\lapack-3.4.1\lapack.pc.in
del /q .\lapack-3.4.1\lapack_build.cmake
del /q .\lapack-3.4.1\lapack_testing.py
del /q .\lapack-3.4.1\LICENSE
del /q .\lapack-3.4.1\make.inc.example
del /q .\lapack-3.4.1\Makefile
del /q .\lapack-3.4.1\README

