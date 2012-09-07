@echo off

rem Delete all SPLIB files

rmdir /s /q .\splib\blas
rmdir /s /q .\splib\include
rmdir /s /q .\splib\precond
rmdir /s /q .\splib\solvers
rmdir /s /q .\splib\tools

del /q .\splib\driver.f
del /q .\splib\LICENSE
del /q .\splib\make.driver
del /q .\splib\make.inc
del /q .\splib\make.splib
del /q .\splib\makefile.lnk
del /q .\splib\matrix
del /q .\splib\methods
del /q .\splib\README
del /q .\splib\rmfiles
del /q .\splib\samples
del /q .\splib\splib.F
del /q .\splib\splib.f.lnk
