@echo off

rem Delete all LAPACK files

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

