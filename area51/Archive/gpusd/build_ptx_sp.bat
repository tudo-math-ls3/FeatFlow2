@echo off
rem pre-process file
cl /EP /C ./ptx/gpusd_ptx.h >./ptx/gpusd_sp.ptx
rem compile file
nvcc -cubin -o ./ptx/gpusd_sp.cubin ./ptx/gpusd_sp.ptx

