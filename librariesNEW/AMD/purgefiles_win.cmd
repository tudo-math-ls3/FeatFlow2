@echo off

rem Delete all AMD files

rmdir /s /q .\AMD\Demo
rmdir /s /q .\AMD\Doc
rmdir /s /q .\AMD\Include
rmdir /s /q .\AMD\Lib
rmdir /s /q .\AMD\MATLAB
rmdir /s /q .\AMD\Source
del /q .\AMD\Makefile
del /q .\AMD\README.txt
