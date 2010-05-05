#!/bin/bash

# set grid base name (expected: ./pre/<basename>.tri, ./pre/<basename>.prm)
grid="unitsquare"

# choose shear modulus
mus="0.5"

# choose Poisson ratio
nus="0.5"

# set solver file base name (expected: ./dat/<basename>.dat)
solver="UMFPACK"
#solver="BICGSTAB"
#solver="MG"
#solver="BICGSTAB_MG"
#solver="MG_BICGSTAB"
#solver="BICGSTAB_MG_BICGSTAB"

# choose MG levels
#levelMaxs="08"
levelMaxs="02 03 04 05"
levelMin="0"

#-----------------------------

# set default log directory
logdir="log"
# make sure default log directory exists
test -d "log" || mkdir -p "log"

# To be able to run several instances of the same binary concurrently, one can call
# this script via 'do.sh <foo>'. Then, dat files and log files are created in the
# subdirectory './log/<foo>'.
if test $# -gt 0; then
  logdir="$logdir/$1"
  # make sure log subdirectory exists
  test -d $logdir || mkdir -p $logdir
fi

# loop over all mu
for mu in ${mus}; do

# loop over all nu
for nu in ${nus}; do

# set base name of dat-file and result file
datFile="$logdir/${grid}_mu${mu}_nu${nu}_mixed"
resultFile="$logdir/${grid}_mu${mu}_nu${nu}_mixed"

# remove existing files
test -f $datFile.dat && \rm $datFile.dat
test -f $resultFile.log && \rm $resultFile.log
# create empty result file
touch $resultFile.log

# loop over all MG levels
for levelMax in ${levelMaxs}; do

# create the temporary dat file
cat > ${datFile}.dat <<END_OF_DATA
# PRM-file of the domain
gridFilePRM = './pre/${grid}.prm'

# TRI-file of the mesh
gridFileTRI = './pre/${grid}.tri'

# type of equatio to solve ('Poisson' or 'elasticity')
equation = elasticity

# FE formulation ('displ' (= pure displacement), 'mixed' (=mixed u/p), 'Stokes')
formulation = mixed

# material parameters (Poisson ratio nu and shear modulus mu)
nu = ${nu}
mu = ${mu}

# type of simulation ('real' or 'analytic')
simulation = analytic

# boundary conditions
# bc[i](#segments x #blocks)
# type of BC for the three components (values can be omitted when simulation = analytical) 
bc1(4) =
'D' 'D' 'N'
'N' 'N' 'N'
'D' 'D' 'N'
'D' 'D' 'N'

# given constant volume force in x- and y-direction in case of real simulation
forceVolumeX   = 0.0
forceVolumeY   = 0.0

# ID of analytical function for u1 and u2 in case of analytical simulation
funcID_u1 = 12
funcID_u2 = 14
funcID_p  = 28
#funcID_u1 = 4
#funcID_u2 = 52
#funcID_p  = 0

# finite element discretisation ('Q1' or 'Q2')
element = Q2

# FE discretisation of the pressure space ('Q1' or 'Q2')
# (only necessary in case of the mixed formulation or Stokes) 
elementPress = P1

# minimum and maximum grid level
levelMin = ${levelMin}
levelMax = ${levelMax}

# solver
solverFile = ./dat/${solver}.dat

# show deformation in visual output ('YES' or 'NO')
showDeformation = NO

# x- and y-coordinate of points where the FE solution is to be evaluated
#evalPoints(1) =
#0.6   0.6

# reference solution values for u1 and u2 in evaluation points
#refSols(1) =
#3.5999625438E-02  3.2721335614E-02  1.1203535942E-01

END_OF_DATA
# dat file has been created now


# start the program, write screen output to result file
echo "*** ************************************"
echo "*** Start computation on level ${levelMax}..."
echo "*** ************************************"
./elasticity ${datFile}.dat | tee ${resultFile}.log
echo "*** ************************************"
echo "*** ... computation finished!"
echo "*** ************************************"
echo "*** stored dat file ${datFile}.dat"
echo ""
echo "--------------------------------------------------------------------------------"
echo ""

# collect the logs of all MG levels in one file
echo "#######################################################" >> ${resultFile}.allMG.log
echo "###            MG level $levelMax                          ###" >> ${resultFile}.allMG.log
echo "#######################################################" >> ${resultFile}.allMG.log
cat ${resultFile}.log >> ${resultFile}.allMG.log
echo "" >> ${resultFile}.allMG.log

done # levelMaxs

# move/rename dat file and result file
\mv ${resultFile}.allMG.log ${resultFile}.log
echo "*** stored the results of all levels in ${resultFile}.log"

done # nus
done # mus

