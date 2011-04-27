#!/bin/bash

# Use 4 OpenMP threads
export OMP_NUM_THREADS=4

# Use Rusanov-type dissipation
export idissipationtype=3

# Use Rusanov-type dissipation in preconditioner
export isystemprecond=1

# Use low-order method
export istabilisation=12



# Grid 1
# -------
export ilev=1
export dstep=1e-3
export gridfile=zpinch_6x120_rz

(cd .. && ./flagship-pc64-nehalem-linux-intel-goto2 data/benchmark/hydro/2d/shocktube.dat)


# Grid 2
# -------
export ilev=1
export dstep=5e-4
export gridfile=zpinch_6x240_rz

(cd .. && ./flagship-pc64-nehalem-linux-intel-goto2 data/benchmark/hydro/2d/shocktube.dat)


# Grid 3
# -------
export ilev=1
export dstep=2.5e-4
export gridfile=zpinch_6x480_rz

(cd .. && ./flagship-pc64-nehalem-linux-intel-goto2 data/benchmark/hydro/2d/shocktube.dat)


# Grid 4
# -------
export ilev=1
export dstep=1.25e-4
export gridfile=zpinch_6x960_rz

(cd .. && ./flagship-pc64-nehalem-linux-intel-goto2 data/benchmark/hydro/2d/shocktube.dat)


# Grid 5
# -------
export ilev=1
export dstep=6.25e-5
export gridfile=zpinch_6x1920_rz

(cd .. && ./flagship-pc64-nehalem-linux-intel-goto2 data/benchmark/hydro/2d/shocktube.dat)

# Grid 6
# -------
export ilev=1
export dstep=3.125e-5
export gridfile=zpinch_6x3840_rz

(cd .. && ./flagship-pc64-nehalem-linux-intel-goto2 data/benchmark/hydro/2d/shocktube.dat)


# Grid 7
# -------
export ilev=1
export dstep=1.5625e-5
export gridfile=zpinch_6x7680_rz

(cd .. && ./flagship-pc64-nehalem-linux-intel-goto2 data/benchmark/hydro/2d/shocktube.dat)
