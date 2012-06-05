#!/bin/bash

# Use 4 OpenMP threads
export OMP_NUM_THREADS=4

# Use Rusanov-type dissipation
export idissipationtype=3

# Use Rusanov-type dissipation in preconditioner
export isystemprecond=3

# Use low-order method
export istabilisation=1

# Level 7
# -------
export ilev=7
export dstep=1e-4

(cd .. && ./flagship-pc64-nehalem-linux-intel-mkl data/benchmark/hydro/2d/doublemach.dat)

# Level 8
# -------
export ilev=8
export dstep=5e-5

(cd .. && ./flagship-pc64-nehalem-linux-intel-mkl data/benchmark/hydro/2d/doublemach.dat)

# Level 9
# -------
export ilev=9
export dstep=2.5e-5

(cd .. && ./flagship-pc64-nehalem-linux-intel-mkl data/benchmark/hydro/2d/doublemach.dat)

# Level 10
# -------
export ilev=10
export dstep=1.25e-5

(cd .. && ./flagship-pc64-nehalem-linux-intel-mkl data/benchmark/hydro/2d/doublemach.dat)