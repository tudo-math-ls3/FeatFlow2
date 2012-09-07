@echo off

rem Delete all UMFPACK files

rmdir /s /q .\UMFPACK\Demo
rmdir /s /q .\UMFPACK\Doc
rmdir /s /q .\UMFPACK\Include
rmdir /s /q .\UMFPACK\Lib
rmdir /s /q .\UMFPACK\MATLAB
rmdir /s /q .\UMFPACK\Source
rmdir /s /q .\UMFPACK\Tcov

del /q .\UMFPACK\Makefile
del /q .\UMFPACK\README.txt
