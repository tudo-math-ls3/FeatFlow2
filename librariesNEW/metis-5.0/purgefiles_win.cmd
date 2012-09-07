@echo off

rem Delete all Metis files

rmdir /s /q .\metis-5.0\build
rmdir /s /q .\metis-5.0\GKlib
rmdir /s /q .\metis-5.0\graphs
rmdir /s /q .\metis-5.0\include
rmdir /s /q .\metis-5.0\libmetis
rmdir /s /q .\metis-5.0\manual
rmdir /s /q .\metis-5.0\programs

del /q .\metis-5.0\BUILD.txt
del /q .\metis-5.0\BUILD-Windows.txt
del /q .\metis-5.0\Changelog
del /q .\metis-5.0\CMakeLists.txt
del /q .\metis-5.0\Install.txt
del /q .\metis-5.0\LICENSE.txt
del /q .\metis-5.0\Makefile
del /q .\metis-5.0\vsgen.bat
