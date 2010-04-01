#!/bin/bash

# set grid base name
grid="cook"

# choose shear modulus
mus="80.194"

# choose Poisson ratio
nus="0.3"

# choose MG levels
mgs="03"
#mgs="02 03 04 05 06 07 08 09 10"

#-----------------------------

# set temporary base name of dat-file and result file
datFile="elasticity_2d_disp_smallDeform_static"
resultFile="result"


# loop over all mu
for mu in ${mus}; do

# loop over all nu
for nu in ${nus}; do

# loop over all MG levels
for mg in ${mgs}; do

# create the temporary dat file
cat > dat/${datFile}.dat <<END_OF_DATA
# PRM-file of the domain
gridFilePRM = './pre/${grid}.prm'

# TRI-file of the mesh
gridFileTRI = './pre/${grid}.tri'

# boundary segments
numBoundarySegments(1) =
4

# type of equatio to solve ('Poisson' or 'elasticity')
equation = elasticity

# material parameters (Poisson ratio nu and shear modulus mu)
nu = ${nu}
mu = ${mu}

# boundary conditions
# bc[i](#segments x #blocks)
# 'D' or 'N', for each boundary i, each segment and each component
bc1(8) =
N
N

N
N

N
N

D
D

# type of simulation ('real' or 'analytic')
simulation = real

# surface forces in x- and y-direction
# forceSurface[i](#segments x dim), in case of real simulation, for each boundary i
# and each segment
forceSurface1(8) =
0.0
0.0

0.0
15.625

0.0
0.0

0.0
0.0

# given constant volume force in x- and y-direction in case of real simulation
forceVolumeX   = 0.0
forceVolumeY   = 0.0

# ID of analytical function for u1 and u2 in case of analytical simulation
funcID_u1 = 0
funcID_u2 = 0

# finite element discretisation ('Q1' or 'Q2')
element = Q1

# minimum and maximum grid level
levelMin = 1
levelMax = ${mg}

# solver ('DIRECT', 'CG', 'BICGSTAB', 'MG', 'CG_MG', 'MG_CG' or 'MG_BICGSTAB')
solver = MG

# maximum number of iterations
numIter = 1000

# relative stopping criterion
tolerance = 1.0e-8

# elementary preconditioner/smoother 'Jacobi' or 'ILU' (only for MG solver)
elementaryPrec = ILU

# MG cycle ('V', 'F' or 'W') (only for MG solver)
mgCycle = V

# number of smoothing steps (only for MG solver)
numSmoothingSteps = 2

# damping parameter (only for MG solver)
damping = 0.7

# show deformation in visual output ('YES' or 'NO')
showDeformation = YES

# x- and y-coordinate of points where the FE solution is to be evaluated
evalPoints(4) =
48.0
60.0

48.0
44.0

# reference solution values for u1 and u2 in evaluation points
refSols(2) =
-20.648992
27.642747

END_OF_DATA
# dat file has been created now


# start the program, write screen output to result file
echo "*** Start computation on level ${mg}..."

./elasticity | tee ${resultFile}.log

echo "*** ... computation finished!"

# move/rename dat file and result file
datFile2="log/${grid}_mu${mu}_nu${nu}_mg${mg}"
resultFile2="log/${grid}_mu${mu}_nu${nu}_mg${mg}"

\mv dat/${datFile}.dat ${datFile2}.dat
\mv ${resultFile}.log ${resultFile2}.log

echo "*** Moved dat file to ${datFile2}.dat."
echo "*** Moved result file to ${resultFile2}.log."

echo ""
echo "*** -----------------------------"
echo ""

done # mgs
done # nus
done # mus

