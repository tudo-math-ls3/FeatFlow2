@echo off

rem Delete all LAPACK files

rmdir /s /q .\lapack-3.4.2\BLAS
rmdir /s /q .\lapack-3.4.2\CMAKE
rmdir /s /q .\lapack-3.4.2\DOCS
rmdir /s /q .\lapack-3.4.2\INSTALL
rmdir /s /q .\lapack-3.4.2\lapacke
rmdir /s /q .\lapack-3.4.2\SRC
rmdir /s /q .\lapack-3.4.2\TESTING

del /q .\lapack-3.4.2\CMakeLists.txt
del /q .\lapack-3.4.2\CTestConfig.cmake
del /q .\lapack-3.4.2\CTestCustom.cmake.in
del /q .\lapack-3.4.2\lapack.pc.in
del /q .\lapack-3.4.2\lapack_build.cmake
del /q .\lapack-3.4.2\lapack_testing.py
del /q .\lapack-3.4.2\LICENSE
del /q .\lapack-3.4.2\make.inc.example
del /q .\lapack-3.4.2\Makefile
del /q .\lapack-3.4.2\README

