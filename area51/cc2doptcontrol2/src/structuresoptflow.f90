!##############################################################################
!# ****************************************************************************
!# <name> structuresoptc </name>
!# ****************************************************************************
!#
!# <purpose>
!# Underlying structures of the space-time optimal control solver.
!# </purpose>
!##############################################################################

module structuresoptflow

  use fsystem
  use storage
  use boundary
  use triangulation
  use cubature
  use paramlist
  use discretebc
  use discretefbc
  use fparser
  use linearsystemscalar
  use multilevelprojection
  
  use collection

  use spatialdiscretisation
  use timediscretisation
  use timescalehierarchy

  use meshhierarchy
  use fespacehierarchybase
  use fespacehierarchy
  use spacetimehierarchy

  use spacetimevectors
  use analyticsolution
  use spacetimeinterlevelprojection
  
  use assemblytemplates
  use assemblytemplatesoptc

  use constantsoptc
  use structuresoptc
  
  implicit none
  
  private
  
  public :: t_optcPostprocessing
  public :: t_settings_optflow
  public :: t_settings_discr
  
!<types>

!<typeblock>

  ! Postprocessing structure. Whenever the postprocessing calculation routines are
  ! called, they fills this structure using the current solution vector (and other
  ! information if necessary). The information in this structure can then be used
  ! for GMV output e.g.
  type t_optcPostprocessing

    ! <!-- Input parameters -->
    
    ! Physics of the problem
    type(t_settings_physics), pointer :: p_rphysics
  
    ! Type of output file to generate from solutions.
    ! 0=disabled
    ! 1=GMV
    ! 2=AVS
    ! 3=Paraview (VTK)
    ! 4=Matlab
    integer :: ioutputUCD = 0
    
    ! Filename for UCD output.
    ! A timestep index '.0001','.0002',... is appended to this.
    character(len=SYS_STRLEN) :: sfilenameUCD = "./gmv/u"
    
    ! Output format for the final solution.
    ! =0: don't write
    ! =1: write out, use formatted output (default).
    ! =2: write out, use unformatted output.
    integer :: cwriteFinalSolution = 1

    ! Filename of a file sequence where the final solution is saved to.
    ! ="": Disable.
    character(len=SYS_STRLEN) :: sfinalSolutionFileName = ""

    ! Output format for the final control.
    ! =0: don't write
    ! =1: write out, use formatted output (default).
    ! =2: write out, use unformatted output.
    integer :: cwriteFinalControl = 1

    ! Filename of a file sequence where the final control is saved to.
    ! ="": Disable.
    character(len=SYS_STRLEN) :: sfinalControlFileName = ""
    
    ! Whether to calculate the values of the optimisation functional
    ! J(.) as well as ||y-y0|| etc. during the postprocessing of space-time vectors.
    integer :: icalcFunctionalValues = 0

    ! Whether to calculate the error to the analytic reference function
    ! ranalyticRefFunction during the postprocessing of space-time vectors.
    integer :: icalcError = 0

    ! Analytic reference function to be used during error calculation
    ! if icalcError > 0.
    type(t_anSolution) :: ranalyticRefFunction
    
    ! Whether to calculate drag/lift forces.
    integer :: icalcForces = 0
    
    ! Boundary component where to calculate body forces
    integer :: ibodyForcesBdComponent = 0
    
    ! 1st coefficient in the boundary integral of the drag coefficient.
    ! If this is commented out, 1/RE is assumed.
    real(DP) :: dbdForcesCoeff1 = 0.0_DP
    
    ! 2nd coefficient in the boundary integral of the drag coefficient.
    ! If this is commented out, 0.004 is assumed (corresonds to flow
    ! around cylinder with RE=1000: Umean=0.2, len=0.1
    ! -> coeff = Umean^2*len = 0.04*0.1 = 0.004 )
    real(DP) :: dbdForcesCoeff2 = 0.0_DP

    ! Whether to write drag/lift forces to hard disc.
    integer :: iwriteBodyForces = 0

    ! Filename for the body forces
    character(len=SYS_STRLEN) :: sfilenameBodyForces = ""

    ! Whether to calculate flux
    integer :: icalcFlux = 0

    ! Start/End coordinates of the lines along which to calculate the flux.
    real(DP), dimension(4) :: Dfluxline = (/0.0,0.0,0.0,0.0/)

    ! Whether to write the flux to a file
    integer :: iwriteFlux = 0

    ! Filename for the flux
    character(len=SYS_STRLEN) :: sfilenameFlux = ""

    ! Whether to calculate flux
    integer :: icalcKineticEnergy = 0

    ! Whether to write the flux to a file
    integer :: iwriteKineticEnergy = 0

    ! Filename for the flux
    character(len=SYS_STRLEN) :: sfilenameKineticEnergy = ""
    ! <!-- the following parameters are automatically maintained during a simulation -->
    
    ! Space that is available in rsolution. One of the CCSPACE_xxxx constants.
    integer :: cspace = CCSPACE_PRIMALDUAL
    
    ! Underlying space discretisation of the primal space
    type(t_blockDiscretisation), pointer :: p_rspaceDiscrPrimal

    ! Underlying space discretisation
    type(t_blockDiscretisation), pointer :: p_rspaceDiscr

    ! Underlying time discretisation
    type(t_timeDiscretisation), pointer :: p_rtimeDiscr
    
    ! Space discretisation based on P1/Q1/P0/Q0 for visualisation output.
    type(t_blockDiscretisation) :: rspaceDiscrLinear

    ! Boundary conditions to use.
    type(t_optcBDC), pointer :: p_rboundaryConditions => null()
    
    ! Points coordinates where to evaluate point values. =NULL: Do not evaluate
    ! point values.
    real(DP), dimension(:,:), pointer :: p_DcoordsPointEval => null()
    
    ! Type of the point value to evaluate. Every entry corresponds to one
    ! point coordinate in p_DcoordsPointEval. The tuples are formed by
    ! (type,der) with
    !   type=1: primal x-velocity, =2: primal y-velocity, =3: primal pressure,
    !       =4: dual x-velocity, =5: dual y-velocity, =6: dual pressure,
    !   der =0: function value, =1: x-derivative, =2: y-derivative
    integer, dimension(:,:), pointer :: p_ItypePointEval => null()
    
    ! Whether or not to write the point values to a file.
    integer :: iwritePointValues = 0
    
    ! Filename for the point values if iwritePointValues <> 0.
    character(len=SYS_STRLEN) :: sfilenamePointValues = ""
        
!    ! A discretisation structure that describes a piecewise constant discretisation
!    ! (usually P0 or Q0).
!    type(t_spatialDiscretisation) :: rdiscrConstant
!
!    ! A discretisation structure that describes a piecewise linear discretisation
!    ! (usually P1 or Q1).
!    type(t_spatialDiscretisation) :: rdiscrLinear
!
!    ! A discretisation structure that describes a piecewise quadratic discretisation
!    ! (usually P2 or Q2).
!    type(t_spatialDiscretisation) :: rdiscrQuadratic
!
!    ! A vector that describes the X-velocity field in the vertices
!    type(t_vectorScalar) :: rvectorVelX
!
!    ! A vector that describes the Y-velocity field in the vertices
!    type(t_vectorScalar) :: rvectorVelY
!
!    ! A vector that describes the pressure field in the vertices
!    type(t_vectorScalar) :: rvectorPressure
!
!    ! A vector that describes the pressure field in the cells
!    type(t_vectorScalar) :: rvectorPressureCells
!
!    ! A vector that describes the streamfunction
!    type(t_vectorScalar) :: rvectorStreamfunction
!
!    ! A vector that describes the H1-error of the velocity field in the vertices
!    type(t_vectorScalar) :: rvectorH1err
!
!    ! A vector that describes the H1-error of the pressure in the cells
!    type(t_vectorScalar) :: rvectorH1errCells
!
  end type

!</typeblock>

!<typeblock>
  
  ! Structure encapsuling the main discretisation
  type t_settings_discr
  
    ! Type of element pair to use for the discretisation.
    ! 0 = Q1~(E031) / Q1~(E031) / Q0
    ! 1 = Q1~(E030) / Q1~(E030) / Q0
    ! 2 = Q1~(EM31) / Q1~(EM31) / Q0
    ! 3 = Q1~(EM30) / Q1~(EM30) / Q0 = standard
    ! 4 = Q2 (E013) / Q2 (E013) / QP1
    ! 5 = Q1~(EM30) / Q1~(EM30) / Q0 unpivoted (much faster than 3 but less stable)
    ! 6 = Q1~(EM30) / Q1~(EM30) / Q0 unscaled (slightly faster than 3 but less stable)
    ! (EM30 = nonparametric, nonconformal Rannacher-Turek element)
    ! (QP1  = Quadrilateral discontinuous P1 element)
    integer :: ielementType = 3
    
    ! cubature formula for Mass matrix
    integer(i32) :: icubMass = CUB_GEN_AUTO

    ! cubature formula for Stokes/Laplacian matrix
    integer(i32) :: icubStokes = CUB_GEN_AUTO

    ! cubature formula for Pressure matrices B
    integer(i32) :: icubB = CUB_GEN_AUTO

    ! cubature formula for RHS F
    integer(i32) :: icubF = CUB_GEN_AUTO

  end type

!</typeblock>

!<typeblock>

  ! Structure collecting all settings of the space-time optflow solver.
  type t_settings_optflow

    !<!--
    ! ################
    ! INPUT PARAMETERS
    ! ################
    !
    ! The following parameters must be set by the main program in order to
    ! run the solver.
    ! -->

    ! An object for saving the domain
    type(t_boundary) :: rboundary

    ! Physics of the problem
    type(t_settings_physics) :: rphysicsPrimal

    ! Parameters of the optimal control problem
    type(t_settings_optcontrol) :: rsettingsOptControl

    ! Settings controlling the spatial discretisation (element, cubature)
    type(t_settings_discr) :: rsettingsSpaceDiscr

    ! Stabilisation parameters for the primal and dual system.
    type(t_settings_stabil) :: rstabilPrimal
    type(t_settings_stabil) :: rstabilDual

    ! Space coarse mesh without refinement
    type(t_triangulation) :: rtriaCoarse
    
    ! Time coarse mesh without refinement
    type(t_timeDiscretisation) :: rtimeCoarse

    ! Settings that define the refinement in space
    type(t_settings_refinement) :: rrefinementSpace

    ! Settings that define the refinement in time
    type(t_settings_refinement) :: rrefinementTime
    
    ! A mesh hierarchy with all available space meshes.
    type(t_meshHierarchy) :: rmeshHierarchy
    
    ! A hierarchy of time levels
    type(t_timescaleHierarchy) :: rtimeHierarchy
    
    ! A level info hierarchy for the assembly of stuff on all levels.
    type(t_staticSpaceAsmHierarchy) :: rspaceAsmHierarchy

    ! A level info hierarchy for the assembly of optimal control stuff on all levels.
    type(t_staticSpaceAsmHierarchyOptC) :: rspaceAsmHierarchyOptC

    ! A hierarchy of space levels for velocity+pressure (primal/dual space)
    type(t_feHierarchy) :: rfeHierPrimal
    
    ! A hierarchy of space levels for velocity+pressure (primal+dual space)
    type(t_feHierarchy) :: rfeHierPrimalDual

    ! Projection hierarchy for the interlevel projection in space (primal+dual space).
    type(t_interlevelProjectionHier) :: rprjHierSpacePrimal
    
    ! Projection hierarchy for the interlevel projection in space (primal/dual space).
    type(t_interlevelProjectionHier) :: rprjHierSpacePrimalDual

    ! A space-time hierarchy based on the primal/dual space
    type(t_spaceTimeHierarchy) :: rspaceTimeHierPrimal
    
    ! A space-time hierarchy based on the primal+dual space
    type(t_spaceTimeHierarchy) :: rspaceTimeHierPrimalDual

    ! Projection hierarchy for the interlevel projection in space/time (primal+dual space).
    type(t_sptiProjHierarchy) :: rprjHierSpaceTimePrimal
    
    ! Projection hierarchy for the interlevel projection in space/time (primal/dual space).
    type(t_sptiProjHierarchy) :: rprjHierSpaceTimePrimalDual

    ! All debug flags used by the application
    type(t_optcDebugFlags) :: rdebugFlags

    !<!-- ------------------------------------- -->
    !<!-- INITIAL CONDITION, BOUNDARY CONDITION -->
    !<!-- ------------------------------------- -->

    ! Analytic solution defining the initial condition.
    type(t_anSolution) :: rinitialCondition

    ! The boundary conditions
    type(t_optcBDC) :: roptcBDC

    !<!-- ------------------------------- -->
    !<!-- APPLICATION SPECIFIC PARAMETERS -->
    !<!-- ------------------------------- -->

    ! An application specific parameter list that allows to pass parameters
    ! to callback routines.
    type(t_parlist), pointer :: p_rparlist => null()
    
    ! User defined global data.
    type(t_globalData) :: rglobalData
    
  end type
  
!</typeblock>

!!<typeblock>
!
!  ! Application-specific type block for the Nav.St. problem
!  type t_problem
!
!    ! Output level during the initialisation phase.
!    integer :: MSHOW_Initialisation
!
!    ! Output level of the application.
!    integer :: MT_OutputLevel
!
!    ! Physics of the problem
!    type(t_settings_physics) :: rphysicsPrimal
!
!    ! A parameter block for everything that controls the optimal control
!    ! problem.
!    type(t_settings_optcontrol) :: roptcontrol
!
!    ! An object for saving the domain:
!    type(t_boundary) :: rboundary
!
!    ! A collection object that saves structural data and some
!    ! problem-dependent information which is e.g. passed to
!    ! callback routines.
!    type(t_collection)                    :: rcollection
!
!    ! A param list that saves all parameters from the DAT/INI file(s).
!    type(t_parlist)                       :: rparamList
!
!    ! A t_problem_oneshot structure that is only valid during the
!    ! execution of the command parser.
!    type(t_problem_oneshot)               :: rdataOneshot
!
!  end type
!
!!</typeblock>

!</types>

end module
