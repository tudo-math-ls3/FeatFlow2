@echo off

rem Delete all SuiteSparse files

rmdir /s /q .\SuiteSparse_config\xerbla

del /q .\SuiteSparse_config\Makefile
del /q .\SuiteSparse_config\README.txt
del /q .\SuiteSparse_config\SuiteSparse_config.c
del /q .\SuiteSparse_config\SuiteSparse_config.h
del /q .\SuiteSparse_config\SuiteSparse_config.mk
del /q .\SuiteSparse_config\SuiteSparse_config_GPU.mk
del /q .\SuiteSparse_config\SuiteSparse_config_Mac.mk
