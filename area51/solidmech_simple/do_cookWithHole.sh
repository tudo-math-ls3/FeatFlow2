#!/bin/bash

# set grid base name
grid="cookWithHole"

# choose shear modulus
mus="80.194"

# choose Poisson ratio
nus="0.3"

# choose MG levels
mgs="03 04 05 06 07 08 09"
#mgs="02 03 04 05 06 07 08"

#-----------------------------

# set temporary base name of dat-file and result file
datFile="solidmech2d"
resultFile="result"


# loop over all mu
for mu in ${mus}; do

# loop over all nu
for nu in ${nus}; do

# loop over all MG levels
for mg in ${mgs}; do

#nlmin=2
nlmin=$(expr ${mg} - 1)

# create the temporary dat file
cat > dat/${datFile}.dat <<END_OF_DATA
# PRM-file of the domain
sgridFilePRM = './pre/${grid}.prm'

# TRI-file of the mesh
sgridFileTRI = './pre/${grid}.tri'

# Element type to use for the discretisation (Q1, Q2)
selementType = Q1

# boundaries
nboundaries= 2

# boundary segments
NboundarySegments(2) =
4
4

# set boundary conditions for first boundary (Neumann = N, Dirichlet = D)
Sbc1(4) =
N
N
N
D

# set boundary conditions for second boundary (Neumann = N, Dirichlet = D)
Sbc2(4) =
N
N
N
N

# Minimum level of the discretisation
NLMIN = ${nlmin}

# Maximum level of the discretisation
NLMAX = ${mg}

# type of configuration (possible values: REAL, ANALYTICAL)
ctypeOfSimulation = REAL

# set function IDs (only needed in case of ctypeOfSimulation .eq. SIMUL_ANALYTICAL)
cfuncID_u1 = 53
cfuncID_u2 = 52

# deformation(possible values: Y (YES), N (NO))
Deformation = Y

# calculate sol on a point(possible values: Y (YES), N (NO))
inquirePoint = N

# Points where sol will be calculated (only needed in case of inquirePoint .eq. Y)
inquirePointX = 48.0
inquirePointY = 60.0

# Reference sol to calculate error (only needed in case of inquirePoint .eq. Y)
refSolU1 = -20.648992
refSolU2 = 27.642747

# max number of iterations
niterations = 1000

# type of solver (possible values: DIRECT_SOLVER,BICGSTAB_SOLVER,MG_SOLVER, CG_SOLVER)
ctypeOfSolver = MG_SOLVER

# type of smoother (possible values: JACOBI, ILU)
ctypeOfSmoother = JACOBI

# Cycle identifier (0=F-cycle, 1=V-cycle, 2=W-cycle)
ccycle = 1

# number of smoothing steps
nsmoothingSteps = 128

# damping parameter
ddamp = 0.45

# tolerance
dtolerance = 1E-10

# material parameters (Poisson ratio nu and shear modulus mu and damping damp)
dnu = ${nu}
dmu = ${mu}

# set constant RHS values (only needed in case of ctypeOfSimulation .eq. SIMUL_REAL)
drhsVol1   = 0
drhsVol2   = 0

DrhsBoundx1(4) =
0.0
0.0
0.0
0.0

DrhsBoundx2(4) = 
0.0
0.0
0.0
0.0

DrhsBoundy1(4) = 
0.0
15.625
0.0
0.0

DrhsBoundy2(4) = 
0.0
0.0
0.0
0.0

END_OF_DATA
# dat file has been created now


# start the program, write screen output to result file
echo "*** Start computation on level ${mg}..."

./solidmech | tee ${resultFile}.txt

echo "*** ... computation finished!"

# Make sure log directory exists
test -d temp/${grid} || mkdir -p temp/${grid}

# move/rename dat file and result file
datFile2="temp/${grid}/${datFile}_mu${mu}_nu${nu}_mg${mg}"
resultFile2="temp/${grid}/${resultFile}_mu${mu}_nu${nu}_mg${mg}"

\mv dat/${datFile}.dat ${datFile2}.dat
\mv ${resultFile}.txt ${resultFile2}.txt

echo "*** Moved dat file to ${datFile2}.dat."
echo "*** Moved result file to ${resultFile2}.txt."

echo ""
echo "*** -----------------------------"
echo ""

done # mgs
done # nus
done # mus

