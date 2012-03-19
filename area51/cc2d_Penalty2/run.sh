#!/bin/bash

for NLMAX in 7;do
for DLAMBDA in 1e3 1e4 1e5 1e6 1e7 1e8;do

cat >./data/discretisation.dat <<END_OF_DATA

# ------------------------------------------------------------------------
# Spatial discretization
#
# This file contains parameters about the spatial discretization
# of the problem.
# ------------------------------------------------------------------------

###################
[CC-DISCRETISATION]
###################

# Type of the initial solution.
# =0: Use zero solution.
# =1: Read from file. isolutionStart specifies the level of the initial 
#     solution, sinitialSolutionFilename specifies the filename to be read.
#     The file is expected as formatted text file.
# =2: Read from file. isolutionStart specifies the level of the initial 
#     solution, sinitialSolutionFilename specifies the filename to be read.
#     The file is expected as unformatted binary file.
#     (smaller, but processor dependent and not human readable!)
# =3: Analytical initial solution. Create the solution by calling
#     the callback functions coeff_AnalyticSolution_X, 
#     coeff_AnalyticSolution_Y and coeff_AnalyticSolution_P
#     in cccallback.f90.

ctypeInitialSolution = 0

# If ctypeInitialSolution=1: Level of the initial solution in the file
# to be read.
# =1,2,3,...: File sinitialSolutionFilename specifies a solution on 
#             level iinitialSolutionLevel.
#             iinitialSolutionlevel must be in the range
#             nlmin <= iinitialSolutionlevel <= nlmax
# =0,-1,-2,...: File sinitialSolutionFilename specifies a solution on 
#               level (NLMAX+iinitialSolutionLevel).

iinitialSolutionlevel = 0

# Name/path of file with the initial solution

sinitialSolutionFilename = ''

# Element type to use for the initial solution.
# =-1: The initial solution is given in the same FEM space as specified
#      by ielementType below
# >=0: The initial solution is given in the FEM space 
#      ielementTypeInitialSolution. All FEM spaces that are possible
#      for ielementType are available here.

ielementTypeInitialSolution = -1

# Write final solution vector to a file. After a stationary simulation,
# the stationary solution is written out. After a nonstationary simulation,
# the solution after the last timestep or in every timestep is written out.
#
# =0: don't write
# =1: write on level specified by iwriteSolutionLevel, use formatted output.
#     For a nonstationary simulation: Write out the solution after the last 
#     timestep.
# =2: write on level specified by iwriteSolutionLevel, use unformatted output.
#     For a nonstationary simulation: Write out the solution after the last 
#     timestep.
# =3: write on level specified by iwriteSolutionLevel, use formatted output.
#     For a nonstationary simulation: Write out the solution every iwriteSolDeltaSteps 
#     timesteps or every dwriteSolDeltaTime interval of simulation time.
# =4: write on level specified by iwriteSolutionLevel, use unformatted output.
#     For a nonstationary simulation: Write out the solution every iwriteSolDeltaSteps 
#     timesteps or every dwriteSolDeltaTime interval of simulation time.
#
# To write out raw solution data during a nonstationary simulation,
# it's also possible to use the 'film output' feature in the postprocessing.
# While this feature here writes out the solution to one single file
# (always overwriting the previous result), 'film output' will write the
# solution to a sequence of files. The format of both files is the same,
# both can be used as initial solution to be read with ctypeInitialSolution<>0.

cwriteFinalSolution = 0

# Specifies the level where to write the final solution.
# 0=don't write
# =1,2,3,...: Write to file swriteSolutionFilename on
#             level iwriteSolutionLevel.
#             iwriteSolutionLevel must be in the range
#             nlmin <= iwriteSolutionLevel <= nlmax
# =0,-1,-2,...: Write to file swriteSolutionFilename on
#               level (NLMAX+iwriteSolutionLevel).

iwriteSolutionLevel = 0

# Name/path of file for writing the final solution

swriteSolutionFilename = ''

# Timestep counter for writing out solutions if cwriteFinalSolution=3/4.
# The solution will be written out every iwriteSolDeltaSteps timesteps.
# =1: write out every timestep.
# This switch works cumulative to iwriteSolDeltaTime, so the solution is
# written out if either the condition given by iwriteSolDeltaSteps or 
# dwriteSolDeltaTime is reached.

iwriteSolDeltaSteps = 1

# Timestep interval for writing out solutions if cwriteFinalSolution=3/4.
# The solution will be written out after each time interval of length
# dwriteSolDeltaTime. 
# =0.0: write out every timestep.
# This switch works cumulative to iwriteSolDeltaSteps, so the solution is
# written out if either the condition given by iwriteSolDeltaSteps or 
# dwriteSolDeltaTime is reached.

dwriteSolDeltaTime = 0.0

# ------------------------------------------------------------------------

# minimum mg-level; <=0: Use level (NLMAX - |NLMIN|) as coarse grid

NLMIN              = 3

# maximum mg-level = level of computation of nonlinear iteration

NLMAX              = ${NLMAX}

# ------------------------------------------------------------------------

# Activate fbm

ifbm = 0

# Type of equation (ISTOKES); 0=Navier Stokes, 1=Stokes calculation

iequation    = 0

# Subtype of the above equation; 0=Gradient tensor, 1=Deformation tensor
# A value > 0 requires to set iupwind >= 4 to activate the extended assembly!

isubEquation = 0

# Model to use for the viscosity.
# A value > 0 requires to set iupwind >= 4 to activate the extended assembly!
# =0: Constant viscosity NU = 1/RE (standard)
# =1: Power law: nu = nu_0 * z^(dviscoexponent/2 - 1), nu_0 = 1/RE, z=||D(u)||^2+dviscoEps
# =2: Bingham fluid: nu = nu_0 + sqrt(2)/2 * dviscoyield / sqrt(|D(u)||^2+dviscoEps^2), nu_0 = 1/RE
# =3: General viscoplastic fluid: nu = nu_0 +  + sqrt(2)/2 * dviscoyield * z^(dviscoexponent/2 - 1), 
#     nu_0 = 1/RE, z = ||D(u)||^2 + dviscoEps^2
# >3: User defined viscosity defined in the callback function 
#     cccallback.f90: getNonconstantViscosity

cviscoModel = 0

# Exponent parameter for the viscosity model

dviscoexponent = 2.0

# Epsilon regularisation for the viscosity model

dviscoEps = 0.01

# Yield stress constant for Bingham fluid

dviscoyield = 1.0

# Viscosity parmeter 1/NU if the viscosity is constant

RE = 1000

# Type of bdry. conditions.
# 0=DIRICHLET stationary, 
# 1=Pressure drop nonstationary, 
# 2=NEUMANN part. nonstationary (standard)
# Must be >= 1 for moving fictitious boundary components and/or 
# nonstationary boundary conditions.

iBoundary = 2

# Activates the moving-frames formulation.
# 0=deactivate
# 1=activate
# Can be used to track e.g. falling particles.
# When activated, the callback routine getMovingFrameVelocity must return
# the velocity and acceleration of the moving frame in every timestep.
# The acceleration is added to the RHS of the Navier-Stokes equation
# while the velocity influences the Dirichlet boundary condition
# and the posiprocessing.

imovingFrame = 0

# ------------------------------------------------------------------------

# Type of right hand side. 
# 0=zero
# 1=steady inhomog., analytically defined by coeff_RHS_x/coeff_RHS_y 
#   before the simulation
# 2=nonsteady force, analytically defined by coeff_RHS_x/coeff_RHS_y 
#   in every timestep of the simulation
# 3=Stationary RHS prescribed by a file.
# 4=Nonstationary RHS prescribed by a sequence of files.
#   The files are linearly interpolated in time if necessary.
#   The parameters sfilenameRHS, drhsXXXX and irhsXXXX below
#   define the file sequence.

iRHS = 2

# For irhs=4: Basic filename of the file sequence defining the right hand
# side. An index ".00000", ".00001",... is appended to this filename.
# The files are uniformly distributed in the time interval 
# [drhsTimeInit..drhsTimeMax] and interpolated in time if necessary.

sfilenameRHS = ""

# Index of the first RHS-file. Standard=0 which starts the file sequence
# with extension ".00000".

irhsFirstIndex = 0

# Number of files in the file sequence.

irhsFileCount = 0

# Defines if the files in the file sequence are formatted.
# =0: unformatted. =1: Formatted (standard)

irhsFormatted = 1

# Defines the starting point in time of the right hand side defined by
# a file sequence. Standard=default initial time.

drhsTimeInit = %{TIME-DISCRETISATION.dtimeInit}

# Defines the maximum simulation time for the file sequence defining
# the RHS. Standard=Default max. time.

drhsTimeMax = %{TIME-DISCRETISATION.dtimeMax}

# Multiplier for the RHS, X-component. If a RHS was read from a file sequence, 
# it is multiplied by this value before being used.

drhsMultiplyX = 1.0

# Multiplier for the RHS, Y-component. If a RHS was read from a file sequence, 
# it is multiplied by this value before being used.

drhsMultiplyY = 1.0

# ------------------------------------------------------------------------

# Type of element pair to use for the discretisation.
#  0 = Q1~(E031) / Q1~(E031) / Q0
#  1 = Q1~(E030) / Q1~(E030) / Q0
#  2 = Q1~(EM31) / Q1~(EM31) / Q0
#  3 = Q1~(EM30) / Q1~(EM30) / Q0 = standard
#  4 = Q2 (E013) / Q2 (E013) / QP1
#  5 = Q1~(EM30) / Q1~(EM30) / Q0 unpivoted (much faster than 3 but less stable)
#  6 = Q1~(EM30) / Q1~(EM30) / Q0 unscaled (slightly faster than 3 but less stable)
#  7 = Q1~(EM30) / Q1~(EM30) / Q0 new interface implementation
#  8 = Q1~(EB30) / Q1~(EB30) / Q0 (experimental, works only with NLMIN=NLMAX and general VANKA/UMFPACK)
# 10 = Q2~(EB50) / Q2~(EB50) / QP1
# 11 = Q2~(EM50) / Q2~(EM50) / QP1
# 12 = Q2 (E013) / Q2 (E013) / QP1 (nonparametric)
# 13 = Q2 (E013) / Q2 (E013) / QP1 (nonparametric, solves local system; slower than 12 but more stable)
# 14 = Q1~(EB30) / Q1~(EB30) / QP1 (parametric, experimental)
# 15 = Q1~(EB30) / Q1~(EB30) / QP1 (nonparametric, experimental)
# 20 = Q1        / Q1        / Q1 (unstabilised!)
# 30 = P1~(E020) / P1~(E020) / P0 

# (EM30 = nonparametric, nonconformal Rannacher-Turek element)
# (EB30 = parametric, nonconformal Rannacher-Turek element with bubble)
# (QP1  = Quadrilateral discontinuous P1 element)
# (E020 = non-conforming Crouzeix-Raviart element)

ielementType = 2

# ------------------------------------------------------------------------

# Type of stabilization of the convective term. 
# 0=Streamline Diffusion (Q1~, Q2)
# 1=upwind
# 2=unified edge oriented jump stabilisation
# 3=fast unified edge oriented jump stabilisation (precomputed matrix)
# 4=Streamline Diffusion, element independent implementation
# 5=unified edge oriented jump stabilisation with the nonlinear
#   matrix being assembled by the element independent implementation
#   Streamline Diffusion method (Suporting mixed TRI/QUAD-meshes)
# REMARK: To use Q2 with EOJ stabilisation and VANKA smoother, damp
# the VANKA with domega=0.5.

iUpwind = 3

# Relaxation parameter for upwind/Streamline Diffusion/Jump stabilisation.
# Standard values: Streamline diffusion=1.0, Upwind=0.1, Jump stabil=0.01

dUpsam = 0.01

# 2nd Relaxation parameter in the jump stabilisation. Standard = 0.

dUpsamStar = 0.0

# Exponent for edge length weight in the jump stabilisation. Standard = 2
# (corresponds to a weight dupsam*h^2).

deojEdgeExp = 2.0

# Calculation of local H for Streamline Diffusion
# 0= set local H to sqrt(area(T))
# 1= set local H to the length of the way a particle travels through
#    a cell

iLocalH = 0

# Element type to use for the unified edge oriented jump stabilisation.
# Allows to use a different finite element for setting up the stabilisation
# than for discretising the equation. Only used if iUpwind=2.
# =-1: Use ielementType.
# >=0: Use a different FE type for the stabilisation than for the equation.
#      See the documentation of ielementType which numbers are allowed here.

ielementTypeStabil = -1

# ------------------------------------------------------------------------

# Generation of mass matrix. 0=lumped mass,1=real mass matrix
iMass              = 1

# Type of lumping of the mass matrix. 0=usual lumping,1=diagonal lumping
iMassLumpType      = 0

# Construction of the divergence/gradient (B-) matrices.
# =0: Create gradient matrix by weak derivative integral (standard, div=B^T)
# =1: Create gradient matrix by strong derivative integral (div!=B^T).
#     Only possible for non-constant pressure functions (e.g. Q2/QP1 discretisation)

istrongDerivativeBmatrix = 0

# cubature formula for Mass matrix. E.g. "AUTO_G3" for 3-point Gauss.
# ="": use default cubature formula.

scubMass           = 

# cubature formula for Stokes/Laplacian matrix. E.g. "AUTO_G3" for 3-point Gauss.
# ="": use default cubature formula.

scubStokes         = 

# cubature formula for Pressure matrices B. E.g. "AUTO_G3" for 3-point Gauss.
# ="": use default cubature formula.

scubB              = 

# cubature formula for RHS F. E.g. "AUTO_G3" for 3-point Gauss.
# ="": use default cubature formula.

scubF              = 

# cubature formula for the EOJ stabilisation (1D along edges).
# ="": use default cubature formula.

scubEOJ            = 

# ------------------------------------------------------------------------
# Parameters to configure matrix generation on multiple levels.
# This stabilises the preconditioner of the nonlinear iteration
# in case of high aspect ratios in the cells of the mesh.
# Works only with Q1~ discretisation.

# Configure adaptive matrix generation.
# 0=no adaptive matrix gen.=standard; 
# 1=switch to constant mat. gen., depending on dAdMatThreshold; 
# 2=switch mat. gen. depenting on size of ncurrent and neighbour element

iAdaptiveMatrix    = 0

# Treshold parameter: 
# switch construction of coarse grid matrices from standard finite element 
# approach to locally constant interpolation for all rows belonging to elements 
# with aspect ratio >= dAdMatThreshold; 
# 0D0=all elements; 20D0=standard; -1D0=infinity=no matrix modifications

dAdMatThreshold    = 20.0


# ========================================================================
# The following section describes the interlevel projection (prolongation/
# restriction) that should be used when a solution must be projected
# to a higher or lower level in the discretisation
# (e.g. when solving with a multigrid solver).

#############
[CC-PROLREST]
#############

# Projection type for velocity and pressure.
# 0 = hard-coded projection (default)
# 1 = matrix-based projection
iprojTypeVelocity = 0
iprojTypePressure = 0

# Order of the prolongation/restriction to use for velocity components.
# -1=default, 0=constant, 1=linear, 2=quadratic,...

iinterpolationOrderVel   = -1

# Order of the prolongation/restriction to use for pressure components.
# -1=default, 0=constant, 1=linear, 2=quadratic,...

iinterpolationOrderPress = -1

# Prolongation/restriction variant for nonconforming elements 
# E030/EM30/E031/EM31 for the velocity.
# = 0: Use default prol/rest,
# = 1: Use standard prol/rest, equally weighted (1/2 from left, 1/2 from right),
# = 2: Use extended prol/rest, equally weighted (1/2 from left, 1/2 from right),
# = 3: Use extended prol/rest, weighted by element size (L2 projection),
# = 4: Use extended prol/rest, weighted by element size of neighbour element
# To activate extended prol/rest, set this to >= 2 after initialising the
# interlevel projection structure!

iinterpolationVariantVel = 0

# Configuration parameter for extended prolongation of E030/EM30/E031/EM31
# element for the velocity. Only valid if iinterpolationVariant >= 2.
# Aspect-ratio indicator; controls switching to constant prolongation.
# <=1: switch depending on aspect ratio of current element (standard),
#  =2: switch depending on aspect ratio of current element and
#      neighbour element

iintARIndicatorEX3YVel   = 1

# Configuration parameter for extended prolongation of E030/EM30/E031/EM31
# element for the velocity. Only valid if iprolEX3Yvariant >= 2.
# Upper bound aspect ratio; for all elements with higher AR
# the prolongation is switched to constant prol/rest.
# standard = 20.0

dintARboundEX3YVel       = 20.0

END_OF_DATA

cat >./data/postprocessing.dat <<END_OF_DATA

# ------------------------------------------------------------------------
# Postprocessing
#
# The parameters in this file specify the postprocessing of the solutions
# that come from the stationary or nonstationary solver.
# ------------------------------------------------------------------------

###################
[CC-POSTPROCESSING]
###################

# Whether to apply postprocessing to a time interpolated solution or not.
# This is used for error analysis as well as drag/lift calculation but
# not for UCD/film output.
# =0: Apply postprocessing to the raw solution vector.
# =1 (standard): Apply postprocessing to a time-interpolated solution
#   with velocity and pressure matching the same point in time.
#   Necessary for higher order time stepping schemes like Crank-Nicolson.

ipostprocTimeInterpSolution = 1

# ------------------------------------------------------------------------
# UCD output (GMV, AVS, Paraview...)

# Type of output file to generate from solutions.
# 0=disabled
# 1=GMV
# 2=AVS
# 3=Paraview (VTK)
# 4=Matlab
# 5=binary GMV

ioutputUCD = 0

# Level where to write UCD output. Standard=0 = maximum level.
# =1,2,3,...  : Write on level ilevelUCD.
# =0,-1,-2,...: Write on level NLMAX - |ilevelUCD|
# Standard is =0 which writes the solution at the maximum level.
# For a value <> 0, the solution is projected down to the desired
# level and written out there (saves disc space).

ilevelUCD = 1

# Filename for UCD output.
# In a nonstationary simulation, a number '.0001','.0002',... is appended
# to this.
sfilenameUCD = '%{spostdirectory}/u.${Dlambda}.${NLMAX}.vtk'

# ------------------------------------------------------------------------
# Nonstationary UCD output parameters

# Start time for UCD output.

dminTimeUCD = -1.0E100

# End time for UCD output.

dmaxTimeUCD = 1.0E100

# Time difference for UCD output.
# =0.0: generate UCD output for every timestep in the range 
#       dminTimeUCD..dmaxTimeUCD

dtimeDifferenceUCD = 0.0

# For a nonstationary simulation if dtimeDifferenceUCD>0: Interpolate the 
# solution of subsequent timesteps to match the time scale defined by
# dtimeDifferenceUCD.
# =0: Postprocess the solution at the current time. (May lead to
#     non-equidistant timesteps in the UCD export)
# =1: Standard. Interpolate the solution linearly to get an equidistant 
#     timestepping for the UCD export. (May lead to slight time interpolation
#     errors due to linear time interpolation.)
# If dtimeDifferenceUCD=0, the solutions are written out as they come.

iinterpolateSolutionUCD = 1

# Start number for file suffix when writing UCD files (1='.0001',...).
# Whether this is used or not depends on the type of output
# (e.g. GMV uses this)

istartSuffixUCD = 0

# ------------------------------------------------------------------------
# Film output.
# In a nonstationary simulation, this writes the raw solution vectors 
# to a sequence of files. For stationary simulations, this has no effect.

# Type of output file to generate from solutions.
# 0=disabled (standard)
# 1=formatted film output
# 2=unformatted film output

ioutputFilm = 0

# Level where to write Film output. Standard = 0 = maximum level.
# =1,2,3,...  : Write on level ilevelFilm.
# =0,-1,-2,...: Write on level NLMAX - |ilevelFilm|
# Standard is =0 which writes the solution at the maximum level.
# For a value <> 0, the solution is projected down to the desired
# level and written out there (saves disc space).

ilevelFilm = 0

# Basic filename for Film output. A number '.0001','.0002',... is appended
# to this.
sfilenameFilm = '%{ssolutiondirectory}/solution'

# Start time for Film output.

dminTimeFilm = -1.0E100

# End time for Film output.

dmaxTimeFilm = 1.0E100

# Time difference for Film output.
# =0.0: generate Film output for every timestep in the range 
#       dminTimeFilm..dmaxTimeFilm

dtimeDifferenceFilm = 0.0

# For a nonstationary simulation if dtimeDifferenceFilm>0: Interpolate the 
# solution of subsequent timesteps to match the time scale defined by
# dtimeDifferenceFilm.
# =0: Standard. Postprocess the solution at the current time. (May lead to
#     non-equidistant timesteps in the UCD export)
# =1: Interpolate the solution linearly to get an equidistant timestepping
#     for the film export. (May lead to slight time interpolation errors
#     due to linear time interpolation.)
# If dtimeDifferenceFilm=0, the solutions are written out as they come.

iinterpolateSolutionFilm = 0

# Start number for file suffix when writing Film files (1='.0001',...).

istartSuffixFilm = 0

# ------------------------------------------------------------------------
# Error analysis

# Calculate L2-error to reference solution.
# =0: Don't calculate
# =1: Calculate the error

ierrorAnalysisL2 = 0

# Calculate H1-error to reference solution.
# =0: Don't calculate
# =1: Calculate the error

ierrorAnalysisH1 = 0

# Calculate the kinetic energy.
# =0: Don't calculate
# =1: Calculate the energy

icalcKineticEnergy = 0

# Cubature formula for error analysis. E.g. "AUTO_G3" for 3-point Gauss.
# ="": use default cubature formula.

scubError          = 

# Whether to write the L2/H1-error and/or kinetic energy to a file

iwriteErrorAnalysisL2 = 0
iwriteErrorAnalysisH1 = 0
iwriteKineticEnergy = 0

# Filename for the L2/H1-error/kinetic energy if iwriteerrorAnalysisXXXX = 1

sfilenameErrorAnalysisL2 = '%{ssolutiondirectory}/error_L2'
sfilenameErrorAnalysisH1 = '%{ssolutiondirectory}/error_H12'
sfilenameKineticEnergy = '%{ssolutiondirectory}/energy'

# ------------------------------------------------------------------------
# Body forces

# Calculate the body forces
# =0: Don't calculate
# =1: Calculate with standard FEAT-1 method. Constant viscosity.
# =2: Calculate with extended method, supporting mixed meshes,
#     general elements and nonlinear viscosity (standard)
# =3: Calculate using volume forces, standard FEAT-1 method. 
#     Constant viscosity.
# =4: Calculate using volume forces with extended method, supporting mixed
#     meshes, general elements and nonlinear viscosity.

icalcBodyForces = 3

# Specifies the tensor structure to use for the computation.
# =-1: automatic
# = 0: simple gradient tensor
# = 1: full gradient tensor
# = 2: deformation tensor

ibodyForcesFormulation = -1

# Number of the boundary component where to calculate the body forces.
# For flow around cylinder, this is =2 e.g.

ibodyForcesBdComponent = 2

# 1st coefficient in the boundary integral of the drag coefficient.
# If this is commented out, 1/RE is assumed.

# dbdForcesCoeff1 = 0.001

# 2nd coefficient in the boundary integral of the drag coefficient.
# If this is commented out, 0.004 is assumed (corresonds to flow 
# around cylinder with RE=1000: Umean=0.2, len=0.1
# -> coeff = Umean^2*len = 0.04*0.1 = 0.004 )

# dbdForcesCoeff2 = 0.004

# Whether to write the body forces to a file

iwriteBodyForces = 1

# Filename for the body forces

sfilenameBodyForces = '%{ssolutiondirectory}/bdforces_${DLAMBDA}_${NLMAX}'

# ------------------------------------------------------------------------
# Point values

# Calculate the value of different quantities in a number of points
# y evaluating the solution. The number in braces "(n)" defines the number
# of points to evaluate. Afterwards follows a list of points in the
# format
#   "x y type der"
# with x=x-coordinate, y=y-coordinate, type=type of quantity to evaluate,
# der=derivative quantifier. Here:
#   type=1: x-velocity, =2: y-velocity, =3: pressure
#   der =0: function value, =1: x-derivative, =2: y-derivative
#
# Example: Don't evaluate anything:
#   cevaluatePointValues(0) =
#
# Example: Evaluate the pressure in (0.15,0.2) and (0.25,0.2)
#   cevaluatePointValues(2) =
#     0.15 0.2 3 0
#     0.25 0.2 3 0

 cevaluatePointValues(8) =
      0.15 0.2 3 0
      0.25 0.2 3 0
      0.1375 0.2 3 0
      0.2625 0.2 3 0
      0.4 0.2 1 0
      0.4 0.2 2 0
      0.65 0.2 1 0
      0.65 0.2 2 0

# Whether to write the point values to a file

iwritePointValues = 1

# Filename for the point values

sfilenamePointValues = '%{ssolutiondirectory}/pointvalues_${DLAMBDA}_${NLMAX}'

# ------------------------------------------------------------------------
# Flux values

# Calculate the flux through one or multilple lines.
# The number in braces "(n)" defines the number of lines where to
# evaluate the flux. Afterwards follows a list of points in the
# format
#   "x-start y-start x-end y-end"
# that defines the starting and ending point of the line.
#
# Example: Don't evaluate anything:
#   cevaluateFluxValues(0) =
#
# Example: Evaluate the flux through the line (0.5,0.5)->(0.5,0.0)
#   cevaluateFluxValues(1) =
#     0.5 0.5 0.5 0.0

cevaluateFluxValues(0) =

# Whether to write the point values to a file

iwriteFluxValues = 0

# Filename for the point values

sfilenameFluxValues = '%{ssolutiondirectory}/flux'


END_OF_DATA

cat >./data/output.dat <<END_OF_DATA

# ------------------------------------------------------------------------
# Output parameters
# 
# This file contains parameters configuring the general output
# of the program.
# ------------------------------------------------------------------------

[GENERALOUTPUT]

# Level of output when reading parameters from the INI files.
# =0=1: no output, =2: print parameters of .INI files to terminal
MSHOW_Initialisation = 0

# level of output to terminal
# =0: no output, =1: basic output; =2: standard output; 
# >=3: extended output with more and more details
MT_OutputLevel = 1

# Log file for messages. ='': No log file output.
smsgLog = '%{slogdirectory}/output_${DLAMBDA}_${NLMAX}.log'

# Log file for error messages; usually coincides with smsgLog to print
# errors into the same log file as standard messages.
# ='': Use the same log file as smsgLog
serrorLog = ''

# Log file for messages used in regression test benchmarks.
# ='': No benchmark output (standard)
sbenchLog = ''

END_OF_DATA

cat >./data/penalty.dat <<END_OF_DATA
############
[CC-PENALTY]
############

# Activate penalty method and choose the way of calculating penalty parameter
#   - 0.0 - not active
#   - 1.0 - active

dPenalty = 1.0

# Penalty method
#   - 1 - full-Lambda (standard)
#   - 2 - fractional Lambda

ipenalty = 1

# Penalty parameter (standard 1000)
dlambda = 1e6

# itypePenaltyAssem indicates how to assemble Penalty matrix
#   - 1 means only one bilinear form with a cub. formula and non constant coeff cc_Lambda
#   - 2 means two bilinear forms combining simple cub. formula with adaptive cub. formula

itypePenaltyAssem = 1

# Generation of penalty matrix. 0=real mass,1 = HRZ mass
iPenalty_lump = 0

# Element type 
#  0 = Q1~(E031) / Q1~(E031) / Q0
#  1 = Q1~(E030) / Q1~(E030) / Q0
#  2 = Q1~(EM31) / Q1~(EM31) / Q0
#  3 = Q1~(EM30) / Q1~(EM30) / Q0 = standard
#  4 = Q2 (E013) / Q2 (E013) / QP1
#  5 = Q1~(EM30) / Q1~(EM30) / Q0 unpivoted (much faster than 3 but less stable)
# ... (see discretisation.dat)

ielementType_Penalty = 2

# cubature formula for Penalty matrix. 
# ="": use default cubature formula.

scubPenalty = G3X3

# adaptive cubature formula for Penalty matrix. Only for itypePenaltyAssem = 2 
# ="": use default cubature formula.

scubPenalty_sum = TRZ

# ilocalrefinement - indicates the local refinement we want (level n -> pow(2,2*n))

ilocalrefinement     = 5

END_OF_DATA

./cc2d

done
done