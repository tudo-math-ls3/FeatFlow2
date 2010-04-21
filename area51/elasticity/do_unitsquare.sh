#!/bin/bash

# set grid base name (expected: ./pre/<basename>.tri, ./pre/<basename>.prm)
grid="unitsquare"

# choose shear modulus
mus="0.5"

# choose Poisson ratio
nus="0.3"

# set solver file base name (expected: ./dat/<basename>.dat)
solver="BICGSTAB"
#solver="MG"
#solver="BICGSTAB_MG"
#solver="MG_BICGSTAB"
#solver="BICGSTAB_MG_BICGSTAB"

# choose MG levels
levelMaxs="04"
#levelMaxs="02 03 04 05 06 07 08"
levelMin="01"

#-----------------------------

# set temporary base name of dat-file and result file
datFile="elasticity_2d_unitsquare"
resultFile="result"


# loop over all mu
for mu in ${mus}; do

# loop over all nu
for nu in ${nus}; do

# loop over all MG levels
for levelMax in ${levelMaxs}; do

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
D
D

D
D

N
N

N
N

# type of simulation ('real' or 'analytic')
simulation = analytic

# surface forces in x- and y-direction
# forceSurface[i](#segments x dim), in case of real simulation, for each boundary i
# and each segment
forceSurface1(8) =
0.0
0.0

0.0
0.0

0.0
0.0

0.0
0.0

# given constant volume force in x- and y-direction in case of real simulation
forceVolumeX   = 0.0
forceVolumeY   = 0.0

# ID of analytical function for u1 and u2 in case of analytical simulation
funcID_u1 = 4
funcID_u2 = 52

# finite element discretisation ('Q1' or 'Q2')
element = Q1

# minimum and maximum grid level
levelMin = ${levelMin}
levelMax = ${levelMax}

# solver
solverFile = ./dat/${solver}.dat

# show deformation in visual output ('YES' or 'NO')
showDeformation = YES

## x- and y-coordinate of points where the FE solution is to be evaluated
#evalPoints(2) =
#0.5
#0.5

# reference solution values for u1 and u2 in evaluation points
#refSols(2) =
#-20.648992
#27.642747

END_OF_DATA
# dat file has been created now


# start the program, write screen output to result file
echo "*** ************************************"
echo "*** Start computation on level ${levelMax}..."
echo "*** ************************************"

./elasticity dat/${datFile}.dat | tee ${resultFile}.log

echo "*** ************************************"
echo "*** ... computation finished!"
echo "*** ************************************"

# move/rename dat file and result file
datFile2="log/${grid}_mu${mu}_nu${nu}_mg${levelMax}"
resultFile2="log/${grid}_mu${mu}_nu${nu}_mg${levelMax}"

\mv dat/${datFile}.dat ${datFile2}.dat
\mv ${resultFile}.log ${resultFile2}.log

echo "*** Moved dat file to ${datFile2}.dat."
echo "*** Moved result file to ${resultFile2}.log."

echo ""
echo "*** ----------------------------------------------------------"
echo ""

done # levelMaxs
done # nus
done # mus

