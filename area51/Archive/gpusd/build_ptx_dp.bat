@echo off
rem pre-process file
cl /EP /C /D__DOUBLE_PRECISION__ ./ptx/gpusd_ptx.h >./ptx/gpusd_dp.ptx
rem compile file
nvcc -cubin -o ./ptx/gpusd_dp.cubin ./ptx/gpusd_dp.ptx

