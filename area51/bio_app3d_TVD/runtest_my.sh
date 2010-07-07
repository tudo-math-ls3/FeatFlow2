#!/bin/bash

for NLMAX in 2 3 4 5; do
#for PHI in 0.1 0.2 0.3; do
#for REL in 0.1 0.3 0.4; do

################### simulation ohne PSEP ################### 
cat >data/chemotaxis.dat <<END_OF_DATA
# This is the data file for the scalar convection-diffusion-reaction
# equation of the problem CODIRE.
#
#------------------------------------------------------------------------------
# General section

[GENERAL]

# Level where to solve. >= NLMIN !
NLMAX = ${NLMAX}
# Output level. 
# =2 export the last solution vectors only
# =1 export data to gmv files
# =0 no gmv files
OUTPUT = 3
ISTEP_GMV = 2

# gmv output folder. this is just needed by chemotaxiscoupld2
# 0 and nonzero integers alternate the output folder
# =0 outputfolder ~/gmvcpld/
# =1 outputfolder ~/gmvcpld2/
GMVFOLDER = 0

# Setting the initial solution vectors
# =1 referred to as u_0 = u_reference
# =0 referred to u_0 to be constant. The const value is set in COEFFICIENTS
INITSOL = 1

# Whether or not using l2_projection to obtain the ICs or interpolation
L2PROJ = 1

# Whether or not checking neg. solution values 
# If =1 as soon as we compute neg. values the sim stops
checkneg = 0

# Threshold for neg. values
# Should be set properly to 
negthres = -0.1

# Tolerance for the defect correction abort criterium
defectTol = 0.00001
#------------------------------------------------------------------------------
# The 'equation' section configures the different constants in the equation.
#
# The equation reads as follows:
#
#       d/dt u = D_1 * Laplace u - CHI*grad *(u * grad c) + f 
#       d/dt c = D_2 * Laplace c - ALPHA*u*c + BETA*u + g
#
#       the params a_,b_ refer to the initial condition to u and c
#       this is only used for the nontest Chemotaxis Pb
#       e.g.  u_0 = a*dexp(-b*((x-0.5_DP)*(x-0.5_DP)+(y-0.5_DP)*(y-0.5_DP)))
#               c_0 = a*dexp(-b*((x)*(x)+(y)*(y)))

[COEFFICIENTS]

n=1
gamma=1
PHI=1.0
r=1.0
w=1.0
sigma = 1.0
chi = 1.0
d_1 = 1.0
alpha=1.0
beta = 1.0
d_2 = 1.0

# The following coefficients are only used for the nontest Chemotaxis Pb
a_cells = 1000.0
a_chemo = 500.0
b_cells = 100.0
b_chemo = 50.0

# The following parameter is used to initialize the cell-concentration to a given mass
mass = 26

# The following stands for the initial sol vector c_0. One expect this vector to approach
# the stat. Limes-sol c_num ~ c_analy = 1. The sol vector u_0 is automaticly set to 
# u_0 = 2.0 = u_analy . Else please specify.
c_0  =  0.0
 u_0  = 0.0

# The following belongs to the chemotaxiscoupldRHS.f90:
# If we want to start with the initial reference fct, we're free to 
# set the scaling params, in order to set
# e.g. u_0 = SCALE_U * x*(x-1)*y*(y-1)
SCALE_C = 1.0
SCALE_U = 2.0
#------------------------------------------------------------------------------
# This section configures the time stepping of the heat conduction equation

[TIMESTEPPING]

# Maximum number of time steps.
# Time iteration stops if #iterations >= niterations or time >= dtimemax.
ntimesteps = 4

# If != 0 the sim will run at this amount of steps neglecting loop ctrl , e.g. error ctrl
steps = 800

# Time step size.
dtstep = 0.1

# Start time
starttime = 0.0

# maximum iterations for defect correction
# used for chemotaxis_cherkur_bilf_nonlin
maxiterationdef = 20

[ERROR]

# the error threshold for steady state
tol = 0.01

[NORM]

# mentioning which norm should be used for the errorcontrol
controlnorm = 0.4

[TEST]

# mentioning which norm should be used for the errorcontrol
convecRelaxation = 1.0
END_OF_DATA

#PHI=`echo $PHI | sed 's/\./_/g'`
#REL=`echo $REL | sed 's/\./_/g'`
# start simulation
#./bio_app
./bio_app > test_lev${NLMAX}_FCT

done
#done
#done