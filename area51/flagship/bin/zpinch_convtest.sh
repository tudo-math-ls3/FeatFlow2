#!/bin/bash

# Use 4 OpenMP threads
export OMP_NUM_THREADS=4

# Use Rusanov-type dissipation
export idissipationtype=3

# Use low-order method
export istabilisation=1

export res=6x120_rz
export dstep=1e-4
(cd .. && ./flagship-pc64-nehalem-linux-intel-mkl data/benchmark/zpinch/2d/zpinch-rz.dat)

export res=6x240_rz
export dstep=5e-5
(cd .. && ./flagship-pc64-nehalem-linux-intel-mkl data/benchmark/zpinch/2d/zpinch-rz.dat)

export res=6x480_rz
export dstep=2.5e-5
(cd .. && ./flagship-pc64-nehalem-linux-intel-mkl data/benchmark/zpinch/2d/zpinch-rz.dat)

export res=6x960_rz
export dstep=1.25e-5
(cd .. && ./flagship-pc64-nehalem-linux-intel-mkl data/benchmark/zpinch/2d/zpinch-rz.dat)

export res=6x1920_rz
export dstep=6.25e-6
(cd .. && ./flagship-pc64-nehalem-linux-intel-mkl data/benchmark/zpinch/2d/zpinch-rz.dat)

export res=6x3840_rz
export dstep=3.125e-6
(cd .. && ./flagship-pc64-nehalem-linux-intel-mkl data/benchmark/zpinch/2d/zpinch-rz.dat)

export res=6x7680_rz
export dstep=1.5625e-6
(cd .. && ./flagship-pc64-nehalem-linux-intel-mkl data/benchmark/zpinch/2d/zpinch-rz.dat)