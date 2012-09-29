!##############################################################################
!# ****************************************************************************
!# <name> newtoniterationlinear </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module realises the linear subsolver of the Newton iteration in the 
!# control space of the optimal control problem.
!# </purpose>
!##############################################################################

module newtoniterationlinear

  use fsystem
  use genoutput
  use paramlist
  use statistics
  
  use spatialdiscretisation
  use timediscretisation
  use linearsystemscalar
  use linearsystemblock
  
  use scalarpde
  use linearformevaluation
  use bilinearformevaluation
  use feevaluation2
  use blockmatassemblybase
  use blockmatassembly
  use collection
  use iterationcontrol
  use matrixmodification
  
  use spacetimevectors
  use analyticsolution
  
  use structuresdiscretisation
  use structuresoptcontrol
  use structuresgeneral
  use structuresoptflow
  use structuresoperatorasm
  use assemblytemplates
  
  use spacematvecassembly
  use spacelinearsolver
  use spacesolver
  
  use kktsystemspaces
  use kktsystem
  use kktsystemhierarchy
  
  use spacetimehierarchy
  use spacetimeinterlevelprj
  
  use structuresnewton
  use postprocessing
  
  implicit none
  
  private

!<constants>

!<constantblock description = "Solver identifier">

  ! Undefined solver
  integer, parameter, public :: NLIN_SOLVER_UNDEFINED = 0

  ! Richardson iteration
  integer, parameter, public :: NLIN_SOLVER_RICHARDSON = 1
  
  ! Multigrid iteration
  integer, parameter, public :: NLIN_SOLVER_MULTIGRID = 2

  ! CG iteration
  integer, parameter, public :: NLIN_SOLVER_CG = 3

  ! BiCGStab iteration
  integer, parameter, public :: NLIN_SOLVER_BICGSTAB = 4

!</constantblock>

!<constantblock description = "Solver types; used for statistics.">

  ! Just the calculation of the residual
  integer, parameter, public :: NLIN_STYPE_RESCALC = 0

  ! A general linear solver
  integer, parameter, public :: NLIN_STYPE_LINSOL = 1

  ! A subsolver in a multigrid iteration
  integer, parameter, public :: NLIN_STYPE_MGSUBSOLVER = 2

!</constantblock>

!</constants>

  public :: t_linsolMultigrid
  public :: t_linsolParameters

!<type>

  !<typeblock>

  ! Linear solver statistics
  type t_newtonlinSolverStat
  
    !<-- ALL GRIDS -->
    
    ! Number of iterations necessary for the linear solver
    integer :: niterations = 0
    
    ! Total time necessary for the linear solver
    type(t_timer) :: rtotalTime
    
    ! Time necessary for solving the linearised forward problem on all grids
    type(t_timer) :: rtimeForwardLin

    ! Time necessary for solving the linearised backward problem on all grids
    type(t_timer) :: rtimeBackwardLin

    ! Time necessary for the creation of defects on all grids
    type(t_timer) :: rtimeDefect
    
    ! Time necessary for smoothing on all grids
    type(t_timer) :: rtimeSmoothing
    
    ! Time necessary for prolongation/restriction
    type(t_timer) :: rtimeProlRest

    ! Statistics of the subsolver in space.
    type(t_spaceslSolverStat) :: rspaceslSolverStat
  
    !<-- FINEST GRID -->

    ! Time necessary for solving the linearised forward problem, finest grid
    type(t_timer) :: rtimeForwardLinFine

    ! Time necessary for solving the linearised backward problem, finest grid
    type(t_timer) :: rtimeBackwardLinFine

    ! Time necessary for the creation of defects, finest grid
    type(t_timer) :: rtimeDefectFine
    
    ! Time necessary for smoothing, finest grid
    type(t_timer) :: rtimeSmoothingFine
    
    ! Statistics of the subsolver in space on the highest level
    type(t_spaceslSolverStat) :: rspaceslSolverStatFine

    !<-- COARSE GRID -->
    
    ! Number of iterations necessary for the coarse grid solver
    integer :: niterationsCoarse = 0

    ! Time necessary for solving the linearised forward problem, coarsest grid
    type(t_timer) :: rtimeForwardLinCoarse

    ! Time necessary for solving the linearised backward problem, coarsest grid
    type(t_timer) :: rtimeBackwardLinCoarse

    ! Time necessary for the creation of defects, coarsest grid
    type(t_timer) :: rtimeDefectCoarse

    ! Time necessary for coarse grid solving
    type(t_timer) :: rtimeSolverCoarse

    ! Statistics of the subsolver in space on the lowest level
    type(t_spaceslSolverStat) :: rspaceslSolverStatCoarse

  end type

  !</typeblock>

  public :: t_newtonlinSolverStat

  !<typeblock>
  
  ! Multigrid parameters and structures
  type t_linsolMultigrid

    ! Identifies the cycle.
    ! =0: F-cycle,
    ! =1: V-cycle,
    ! =2: W-cycle
    integer :: ccycle = 0

    ! Specifies the smoother.
    ! =1: Damped Richardson iteration
    ! =2: CG
    ! =3: BiCGStab
    integer :: csmoother = 0

    ! Specifies the coarse grid solver.
    ! =1: Damped Richardson iteration
    ! =2: CG
    ! =3: BiCGStab
    integer :: ccoarseGridSolver = 0

    ! Number of presmoothing steps
    integer :: nsmPre = 0

    ! Number of postsmoothing steps
    integer :: nsmPost = 4
    
    ! Coarse grid correction damping parameter
    real(DP) :: dcoarseGridCorrectionWeight = 1.0_DP
    
    ! Section in the parameter list with parameters for the smoother
    character(LEN=PARLST_MLSECTION) :: ssectionSmoother
    
    ! Section in the parameter list with parameters for the coarse grid solver
    character(LEN=PARLST_MLSECTION) :: ssectionCoarseGridSolver
    
    ! Link to the parameter list
    type(t_parlist), pointer :: p_rparList => null()
    
    ! For every level in the hierarchy, a linear solver structure
    ! that configures the smoother. The solver structure at level 1
    ! defines the coarse grid solver
    type(t_linsolParameters), dimension(:), pointer :: p_rsubSolvers => null()
  
    ! KKT system hierarchy for solutions on the different levels
    type(t_kktsystemDirDerivHierarchy) :: rsolutionHier

    ! KKT system hierarchy for defect vectors on the different levels
    type(t_kktsystemDirDerivHierarchy) :: rdefectHier

    ! KKT system hierarchy for RHS vectors on the different levels
    type(t_kktsystemDirDerivHierarchy) :: rrhsHier
  
    ! KKT system hierarchy for temporary vectors on the different levels
    type(t_kktsystemDirDerivHierarchy) :: rtempHier

    ! Cycle counter
    integer, dimension(:), pointer :: p_IcycleCount => null()

  end type
  
  !</typeblock>
  
  !<typeblock>
  
  ! CG parameters and structures
  type t_linsolCG
    ! Temporary memory 
    type(t_controlSpace) :: rr

    ! Temporary memory 
    type(t_controlSpace) :: rAp
    
    ! More temporary memory
    type(t_kktsystemDirDeriv) :: rp
    
    ! Use real residual in the iteration.
    logical :: brealres = .true.
  end type

  !</typeblock>
  
  public :: t_linsolCG

! *****************************************************************************

!<typeblock>
  
  ! BiCGSTab parameters and substructures
  type t_linsolBiCGStab
  
    ! Temporary vectors to use during the solution process
    type(t_controlSpace), dimension(3) :: RtempVectors

    ! More temporary memory
    type(t_kktsystemDirDeriv) :: rp

    ! More temporary memory
    type(t_kktsystemDirDeriv) :: rr
    
  end type

!</typeblock>
  
  public :: t_linsolBiCGStab
  
! *****************************************************************************
  
  !<typeblock>
  
  ! Richardson parameters and structures
  type t_linsolRichardson
  
    ! Temporary memory on all levels for calculations
    type(t_controlSpace) :: rtempVector
  
  end type

  !</typeblock>
  
  public :: t_linsolRichardson

! *****************************************************************************
  
  !<typeblock>
  
  ! Contains the parameters of the Newton iteration.
  type t_linsolParameters
  
    ! <!-- --------------------------------------- -->
    ! <!-- GENRERAL PARAMETERS AND SOLVER SETTINGS -->
    ! <!-- --------------------------------------- -->

    ! Specifies the type of the solver.
    ! One of the NLIN_SOLVER_xxxx constants.
    integer :: csolverType = NLIN_SOLVER_UNDEFINED

    ! General newton parameters
    type(t_newtonPrecParameters) :: rprecParameters

    ! Iteration control parameters
    type(t_iterationControl) :: riter

    ! Parameters of the OptFlow solver
    type(t_settings_optflow), pointer :: p_rsettingsSolver => null()
  
    ! <!-- ----------------------------- -->
    ! <!-- SUBSOLVERS AND OTHER SETTINGS -->
    ! <!-- ----------------------------- -->

    ! Hierarchy of solvers in space for all levels.
    ! Linearised primal equation.
    type(t_spaceSolverHierarchy), pointer :: p_rsolverHierPrimalLin => null()

    ! Hierarchy of solvers in space for all levels.
    ! Linearised dual equation.
    type(t_spaceSolverHierarchy), pointer :: p_rsolverHierDualLin => null()

    ! Defines a policy how to generate the initial condition of a timestep.
    ! =0: Always take zero
    ! =1: Propagate the solution of the previous/next timestep to the
    !     current one. (Default)
    ! =2: Take the solution of the last space-time iteration
    ! Warning: Avoid the setting cspatialInitCondPolicy = 0/1 if the linear
    ! solver in space does not solve up to a high accuracy.
    ! Otherwise, the solution of the linearised primal equation in the first
    ! timestep is not better than the stopping criterion of the space solver
    ! and thus, the complete space-time solution will not be better!
    ! Remember, there is no nonlinear loop around each timestep!
    ! That means, the stopping criterion of the linear space-time solver is
    ! overwritten by the stopping criterion of the space-solver, which
    ! is usually an undesired behaviour: Although the linear space-time solver
    ! solves up to, e.g., 1E-15, the total solution of the space-time solver
    ! is not better than depsrel(linear space-solver)!!!
    integer :: cspatialInitCondPolicy = SPINITCOND_PREVITERATE
    
    ! Equation flags that specify modifications to the equation to solve.
    ! One of the SPACESLH_EQNF_xxxx constants.
    integer :: ceqnflags = SPACESLH_EQNF_DEFAULT

    ! Parameters for the multigrid solver
    type(t_linsolMultigrid), pointer :: p_rsubnodeMultigrid => null()

    ! Parameters for the Richardson solver
    type(t_linsolRichardson), pointer :: p_rsubnodeRichardson => null()

    ! Parameters for the CG solver
    type(t_linsolCG), pointer :: p_rsubnodeCG => null()
    
    ! Parameters for the BiCGStab solver
    type(t_linsolBiCGStab), pointer :: p_rsubnodeBiCGStab => null()
    
    ! <!-- -------------- -->
    ! <!-- TEMPORARY DATA -->
    ! <!-- -------------- -->
    
    ! Underlying KKT system hierarchy that defines the shape of the solutions.
    type(t_kktsystemHierarchy), pointer :: p_rkktsystemHierarchy => null()

    ! Projection hierarchy for the interlevel projection in space/time, primal space.
    type(t_sptiProjHierarchyBlock), pointer :: p_rprjHierSpaceTimePrimal => null()

    ! Projection hierarchy for the interlevel projection in space/time, dual space.
    type(t_sptiProjHierarchyBlock), pointer :: p_rprjHierSpaceTimeDual => null()

    ! Projection hierarchy for the interlevel projection in space/time, control space.
    type(t_sptiProjHierarchyBlock), pointer :: p_rprjHierSpaceTimeControl => null()
    
  end type
  
  !</typeblock>
  
!</types>

  ! Apply preconditioning of a defect with the Newton preconditioner
  ! in the control space.
  public :: newtonlin_precond
  
  ! Sets the stopping criterion for the use application of the adaptive
  ! Newton algorithm
  public :: newtonlin_adNewton_setEps
  
  ! Basic initialisation of the Newton solver
  public :: newtonlin_init
  
  ! Structural initialisation
  public :: newtonlin_initStructure
  
  ! Propagates the fine-grid solution to the coarse grid solutions
  public :: newtonlin_initNonlinearData
  
  ! Final initialisation
  public :: newtonlin_initData
  
  ! Cleanup of data
  public :: newtonlin_doneData

  ! Cleanup of structures
  public :: newtonlin_doneStructure

  ! Final cleanup
  public :: newtonlin_done
  
  ! Clear a statistics block
  public :: newtonlin_clearStatistics
  
  ! SUm up two statistics blocks
  public :: newtonlin_sumStatistics

contains

  ! ***************************************************************************

!<subroutine>

  subroutine newtonlin_getResidual (&
      rlinsolParam,rkktsystemDirDeriv,rrhs,rresidual,rstatistics,dres,iresnorm)
  
!<description>
  ! Calculates the residual in the linearised control equation.
!</description>
  
!<inputoutput>
  ! Parameters for the iteration.
  type(t_linsolParameters), intent(inout) :: rlinsolParam
  
  ! Structure defining the linearised KKT system.
  type(t_kktsystemDirDeriv), intent(inout), target :: rkktsystemDirDeriv
  
  ! Right-hand side of the linearised control equation
  type(t_controlSpace), intent(inout), target :: rrhs
  
  ! On output, this structure receives a representation of the search
  ! direction / residual in the Newton iteration
  type(t_controlSpace), intent(inout) :: rresidual

  ! type of norm. A LINALG_NORMxxxx constant.
  integer, intent(in), optional :: iresnorm
!</inputoutput>

!<output>
  ! Solver statistics. Statistics of the residual calculation are added to this.
  type(t_newtonlinSolverStat), intent(out) :: rstatistics

  ! L2-Norm of the residual
  real(DP), intent(out), optional :: dres
!</output>

!</subroutine>

    ! local variables
    type(t_spaceslSolverStat) :: rlocalStat

    call stat_startTimer (rstatistics%rtotalTime)
    
    ! -------------------------------------------------------------
    ! Step 1: Solve the primal and dual system.
    ! -------------------------------------------------------------

    ! --------------------
    ! Forward equation
    ! --------------------

    if (rlinsolParam%rprecParameters%ioutputLevel .ge. 3) then
      call output_line ("Linear space-time Residual: Solving the primal equation")
    end if

    ! Solve the primal equation, update the primal solution.
    output_iautoOutputIndent = output_iautoOutputIndent + 2

    call kkt_solvePrimalDirDeriv (rkktsystemDirDeriv,&
        rlinsolParam%p_rsolverHierPrimalLin,rlinsolParam%cspatialInitCondPolicy,&
        rlinsolParam%ceqnflags,rlocalStat)

    output_iautoOutputIndent = output_iautoOutputIndent - 2
    
    if (rlinsolParam%rprecParameters%ioutputLevel .ge. 4) then
      call output_line ("Linear space-time Residual: Time for solving      : "//&
          trim(sys_sdL(rlocalStat%rtotalTime%delapsedReal,10)))
      call output_line ("Linear space-time Residual: Time for space-defects: "//&
          trim(sys_sdL(rlocalStat%rtimeDefect%delapsedReal,10)))
      call output_line ("Linear space-time Residual:   Time for space-RHS  : "//&
          trim(sys_sdL(rlocalStat%rtimeRHS%delapsedReal,10)))
      call output_line ("Linear space-time Residual: Time for mat. assembly: "//&
          trim(sys_sdL(rlocalStat%rtimeMatrixAssembly%delapsedReal,10)))
      call output_line ("Linear space-time Residual: Time for factorisation: "//&
          trim(sys_sdL(rlocalStat%rlssSolverStat%rtimeSymbolicFactorisation%delapsedReal+&
                       rlocalStat%rlssSolverStat%rtimeNumericFactorisation%delapsedReal,10)))
      call output_line ("Linear space-time Residual: Time for space-solver : "//&
          trim(sys_sdL(rlocalStat%rlssSolverStat%rtotalTime%delapsedReal,10)))
    end if

    ! Add time to the time of the forward equation
    call stat_addTimers (rlocalStat%rtotalTime,rstatistics%rtimeForwardLin)
    call spacesl_sumStatistics(rlocalStat,rstatistics%rspaceslSolverStat)
    
    ! --------------------
    ! Backward equation
    ! --------------------

    if (rlinsolParam%rprecParameters%ioutputLevel .ge. 3) then
      call output_line ("Linear space-time Residual: Solving the dual equation")
    end if

    ! Solve the dual equation, update the dual solution.
    output_iautoOutputIndent = output_iautoOutputIndent + 2

    call kkt_solveDualDirDeriv (rkktsystemDirDeriv,&
        rlinsolParam%p_rsolverHierDualLin,rlinsolParam%cspatialInitCondPolicy,&
        rlinsolParam%ceqnflags,rlocalStat)

    output_iautoOutputIndent = output_iautoOutputIndent - 2
    
    if (rlinsolParam%rprecParameters%ioutputLevel .ge. 4) then
      call output_line ("Linear space-time Residual: Time for solving      : "//&
          trim(sys_sdL(rlocalStat%rtotalTime%delapsedReal,10)))
      call output_line ("Linear space-time Residual: Time for space-defects: "//&
          trim(sys_sdL(rlocalStat%rtimeDefect%delapsedReal,10)))
      call output_line ("Linear space-time Residual:   Time for space-RHS  : "//&
          trim(sys_sdL(rlocalStat%rtimeRHS%delapsedReal,10)))
      call output_line ("Linear space-time Residual: Time for mat. assembly: "//&
          trim(sys_sdL(rlocalStat%rtimeMatrixAssembly%delapsedReal,10)))
      call output_line ("Linear space-time Residual: Time for factorisation: "//&
          trim(sys_sdL(rlocalStat%rlssSolverStat%rtimeSymbolicFactorisation%delapsedReal+&
                       rlocalStat%rlssSolverStat%rtimeNumericFactorisation%delapsedReal,10)))
      call output_line ("Linear space-time Residual: Time for space-solver : "//&
          trim(sys_sdL(rlocalStat%rlssSolverStat%rtotalTime%delapsedReal,10)))
    end if

    ! Add time to the time of the backward equation
    call stat_addTimers (rlocalStat%rtotalTime,rstatistics%rtimeBackwardLin)
    call spacesl_sumStatistics(rlocalStat,rstatistics%rspaceslSolverStat)

    ! -------------------------------------------------------------
    ! Step 2: Calculate the residual
    ! -------------------------------------------------------------

    ! Take the solution of the linearised primal/dual system and
    ! calculate the residual in the control space.
    call kkt_calcControlResDirDeriv (rkktsystemDirDeriv,rrhs,rresidual,dres,iresnorm)

    call stat_stopTimer (rstatistics%rtotalTime)
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonlin_applyOperator (rlinsolParam,rkktsystemDirDeriv,rrhs,rstatistics)
  
!<description>
  ! Applies the operator of the linearised control equation.
!</description>
  
!<inputoutput>
  ! Parameters for the iteration.
  type(t_linsolParameters), intent(inout) :: rlinsolParam
  
  ! Structure defining the linearised KKT system.
  type(t_kktsystemDirDeriv), intent(inout), target :: rkktsystemDirDeriv
!</inputoutput>

!<output>
  ! Solver statistics
  type(t_newtonlinSolverStat), intent(out) :: rstatistics

  ! Receives the result: d = J''(u)g
  type(t_controlSpace), intent(inout), target :: rrhs
!</output>

!</subroutine>

    ! local variables
    type(t_spaceslSolverStat) :: rlocalStat

    call stat_startTimer (rstatistics%rtotalTime)

    ! -------------------------------------------------------------
    ! Step 1: Solve the primal and dual system.
    ! -------------------------------------------------------------

    ! --------------------
    ! Forward equation
    ! --------------------

    if (rlinsolParam%rprecParameters%ioutputLevel .ge. 3) then
      call output_line ("Linear space-time Residual: Solving the primal equation")
    end if

    ! Solve the primal equation, update the primal solution.
    output_iautoOutputIndent = output_iautoOutputIndent + 2

    call kkt_solvePrimalDirDeriv (rkktsystemDirDeriv,&
        rlinsolParam%p_rsolverHierPrimalLin,rlinsolParam%cspatialInitCondPolicy,&
        rlinsolParam%ceqnflags,rlocalStat)
    
    output_iautoOutputIndent = output_iautoOutputIndent - 2

    if (rlinsolParam%rprecParameters%ioutputLevel .ge. 4) then
      call output_line ("Linear space-time Residual: Time for solving      : "//&
          trim(sys_sdL(rlocalStat%rtotalTime%delapsedReal,10)))
      call output_line ("Linear space-time Residual: Time for space-defects: "//&
          trim(sys_sdL(rlocalStat%rtimeDefect%delapsedReal,10)))
      call output_line ("Linear space-time Residual:   Time for space-RHS  : "//&
          trim(sys_sdL(rlocalStat%rtimeRHS%delapsedReal,10)))
      call output_line ("Linear space-time Residual: Time for mat. assembly: "//&
          trim(sys_sdL(rlocalStat%rtimeMatrixAssembly%delapsedReal,10)))
      call output_line ("Linear space-time Residual: Time for factorisation: "//&
          trim(sys_sdL(rlocalStat%rlssSolverStat%rtimeSymbolicFactorisation%delapsedReal+&
                       rlocalStat%rlssSolverStat%rtimeNumericFactorisation%delapsedReal,10)))
      call output_line ("Linear space-time Residual: Time for space-solver : "//&
          trim(sys_sdL(rlocalStat%rlssSolverStat%rtotalTime%delapsedReal,10)))
    end if

    ! Add time to the time of the forward equation
    call stat_addTimers (rlocalStat%rtotalTime,rstatistics%rtimeForwardLin)
    call spacesl_sumStatistics(rlocalStat,rstatistics%rspaceslSolverStat)

    ! --------------------
    ! Backward equation
    ! --------------------

    if (rlinsolParam%rprecParameters%ioutputLevel .ge. 3) then
      call output_line ("Linear space-time Residual: Solving the dual equation")
    end if

    ! Solve the dual equation, update the dual solution.
    output_iautoOutputIndent = output_iautoOutputIndent + 2

    call kkt_solveDualDirDeriv (rkktsystemDirDeriv,&
        rlinsolParam%p_rsolverHierDualLin,rlinsolParam%cspatialInitCondPolicy,&
        rlinsolParam%ceqnflags,rlocalStat)

    output_iautoOutputIndent = output_iautoOutputIndent - 2

    if (rlinsolParam%rprecParameters%ioutputLevel .ge. 4) then
      call output_line ("Linear space-time Residual: Time for solving      : "//&
          trim(sys_sdL(rlocalStat%rtotalTime%delapsedReal,10)))
      call output_line ("Linear space-time Residual: Time for space-defects: "//&
          trim(sys_sdL(rlocalStat%rtimeDefect%delapsedReal,10)))
      call output_line ("Linear space-time Residual:   Time for space-RHS  : "//&
          trim(sys_sdL(rlocalStat%rtimeRHS%delapsedReal,10)))
      call output_line ("Linear space-time Residual: Time for mat. assembly: "//&
          trim(sys_sdL(rlocalStat%rtimeMatrixAssembly%delapsedReal,10)))
      call output_line ("Linear space-time Residual: Time for factorisation: "//&
          trim(sys_sdL(rlocalStat%rlssSolverStat%rtimeSymbolicFactorisation%delapsedReal+&
                       rlocalStat%rlssSolverStat%rtimeNumericFactorisation%delapsedReal,10)))
      call output_line ("Linear space-time Residual: Time for space-solver : "//&
          trim(sys_sdL(rlocalStat%rlssSolverStat%rtotalTime%delapsedReal,10)))
    end if

    ! Add time to the time of the backward equation
    call stat_addTimers (rlocalStat%rtotalTime,rstatistics%rtimeBackwardLin)
    call spacesl_sumStatistics(rlocalStat,rstatistics%rspaceslSolverStat)

    ! -------------------------------------------------------------
    ! Step 2: Calculate the residual
    ! -------------------------------------------------------------

    ! Take the solution of the linearised primal/dual system and
    ! calculate the residual in the control space.
    call kkt_applyControlDirDeriv (rkktsystemDirDeriv,rrhs)

    call stat_stopTimer (rstatistics%rtotalTime)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonlin_smoothCorrection (rsmootherParams,nsm,&
      rkktsystemDirDeriv,rrhs,rstatistics)
  
!<description>
  ! Applies the space-time smoother.
!</description>

!<input>
  ! Number of smoothing steps
  integer, intent(in) :: nsm
!</input>

!<inputoutput>
  ! Parameters of the smoother
  type(t_linsolParameters), intent(inout) :: rsmootherParams
  
  ! Structure defining the linearised KKT system.
  ! The linearised control in this structure is used as
  ! initial and target vector.
  ! The control on the maximum level receives the result.
  type(t_kktsystemDirDeriv), intent(inout), target :: rkktsystemDirDeriv
  
  ! Right-hand side of the linearised control equation
  type(t_controlSpace), intent(inout), target :: rrhs
!</inputoutput>

!<ouptut>
  ! Solver statistics. 
  type(t_newtonlinSolverStat), intent(out) :: rstatistics
!</ouptut>

!</subroutine>

    ! Configure the solver to apply exactly nsm steps
    if (nsm .le. 0) return
    rsmootherParams%riter%nminiterations = nsm
    rsmootherParams%riter%nmaxiterations = nsm
    
    ! Call the smoother.
    select case (rsmootherParams%csolverType)
    
    ! --------------------------------------
    ! Richardson iteration
    ! --------------------------------------
    case (NLIN_SOLVER_RICHARDSON)
      ! Call the Richardson iteration to calculate an update
      ! for the control.
      call newtonlin_richardson (rsmootherParams,rkktsystemDirDeriv,rrhs,rstatistics)

    ! --------------------------------------
    ! CG iteration
    ! --------------------------------------
    case (NLIN_SOLVER_CG)
      ! Call the Richardson iteration to calculate an update
      ! for the control.
      call newtonlin_cg (rsmootherParams,rkktsystemDirDeriv,rrhs,rstatistics)

    ! --------------------------------------
    ! BiCGSTab iteration
    ! --------------------------------------
    case (NLIN_SOLVER_BICGSTAB)
      ! Call the Richardson iteration to calculate an update
      ! for the control.
      call newtonlin_bicgstab (rsmootherParams,rkktsystemDirDeriv,rrhs,rstatistics)

    ! --------------------------------------
    ! unknown iteration
    ! --------------------------------------
    case default
              
      call output_line ("Invalid smoother.", &
          OU_CLASS_ERROR,OU_MODE_STD,"newtonlin_smoothCorrection")
      call sys_halt()

    end select
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonlin_richardson (rlinsolParam,rkktsystemDirDeriv,rrhs,rstatistics)
  
!<description>
  ! Applies a Richardson iteration in the control space
  ! for the linearised KKT system.
!</description>

!<inputoutput>
  ! Parameters for the iteration.
  ! The output parameters are changed according to the iteration.
  type(t_linsolParameters), intent(inout) :: rlinsolParam

  ! Structure defining the linearised KKT system.
  ! The linearised control in this structure is used as
  ! initial and target vector.
  ! The control on the maximum level receives the result.
  type(t_kktsystemDirDeriv), intent(inout), target :: rkktsystemDirDeriv
  
  ! Right-hand side of the linearised control equation
  type(t_controlSpace), intent(inout), target :: rrhs
!</inputoutput>

!<ouptut>
  ! Solver statistics. 
  type(t_newtonlinSolverStat), intent(out) :: rstatistics
!</ouptut>

!</subroutine>

    ! Temporary control vector
    type(t_controlSpace), pointer :: p_rtempVector
    type(t_newtonlinSolverStat) :: rlocalStat
    real(DP) :: dres

    ! Measure iteration time
    call stat_startTimer (rstatistics%rtotalTime)

    ! Get the temp vector
    p_rtempVector => rlinsolParam%p_rsubnodeRichardson%rtempVector

    ! Richardson is a do loop...
    call itc_initIteration(rlinsolParam%riter)

    do while (.true.)
    
      ! -------------------------------------------------------------
      ! Get the current residual / search direction
      ! -------------------------------------------------------------

      ! Compute the residual and its norm.
      output_iautoOutputIndent = output_iautoOutputIndent + 2
      call newtonlin_getResidual (rlinsolParam,rkktsystemDirDeriv,rrhs,&
          p_rtempVector,rlocalStat,dres,&
          rlinsolParam%rprecParameters%iresnorm)
      output_iautoOutputIndent = output_iautoOutputIndent - 2
      
      call newtonlin_sumStatistics(rlocalStat,rstatistics,NLIN_STYPE_RESCALC)
      
      if (rlinsolParam%riter%cstatus .eq. ITC_STATUS_UNDEFINED) then
        ! Remember the initial residual
        call itc_initResidual(rlinsolParam%riter,dres)
      else
        ! Push the residual, increase the iteration counter
        call itc_pushResidual(rlinsolParam%riter,dres)
      end if
      
      if (rlinsolParam%rprecParameters%ioutputLevel .ge. 2) then
        call output_line ("Space-time Richardson: Iteration "// &
            trim(sys_siL(rlinsolParam%riter%niterations,10))// &
            ", ||res(u)|| = "// &
            trim(sys_sdEL(dres,10)))
      end if

      ! -------------------------------------------------------------
      ! Check for convergence / divergence / ...
      ! -------------------------------------------------------------
      if (rlinsolParam%riter%cstatus .ne. ITC_STATUS_CONTINUE) exit
      
      ! -------------------------------------------------------------
      ! Update of the solution
      ! -------------------------------------------------------------

      ! Add the residual (damped) to the current control.
      ! This updates the current control.
      call kktsp_controlLinearComb (p_rtempVector,&
          rlinsolParam%rprecParameters%domega,&
          rkktsystemDirDeriv%p_rcontrolLin,1.0_DP)
    
      ! -------------------------------------------------------------
      ! Proceed with the iteration
      ! -------------------------------------------------------------
    
    end do

    call stat_stopTimer (rstatistics%rtotalTime)

    ! Statistics
    rstatistics%niterations = rstatistics%niterations + rlinsolParam%riter%niterations

    if (rlinsolParam%rprecParameters%ioutputLevel .ge. 2) then
      call output_lbrk()
      call output_line ("Space-time Richardson: Statistics.")
      call output_lbrk()
      call itc_printStatistics(rlinsolParam%riter)
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonlin_cg (rlinsolParam,rkktsystemDirDeriv,rrhs,rstatistics)
  
!<description>
  ! Applies a CG iteration in the control space
  ! for the linearised KKT system.
!</description>

!<inputoutput>
  ! Parameters for the iteration.
  ! The output parameters are changed according to the iteration.
  type(t_linsolParameters), intent(inout), target :: rlinsolParam

  ! Structure defining the linearised KKT system.
  ! The linearised control in this structure is used as
  ! initial and target vector.
  ! The control on the maximum level receives the result.
  type(t_kktsystemDirDeriv), intent(inout), target :: rkktsystemDirDeriv
  
  ! Right-hand side of the linearised control equation
  type(t_controlSpace), intent(inout), target :: rrhs
!</inputoutput>

!<ouptut>
  ! Solver statistics. 
  type(t_newtonlinSolverStat), intent(out) :: rstatistics
!</ouptut>

!</subroutine>

    real(DP) :: dalpha,dbeta,dgamma,dgammaOld,dtemp,dres
    type(t_controlSpace), pointer :: p_rr,p_rAp
    type(t_kktsystemDirDeriv), pointer :: p_rp
    type(t_newtonlinSolverStat) :: rlocalStat
    type(t_iterationControl) :: rlocaliter
    logical :: brealres

    ! Measure the total time
    call stat_startTimer (rstatistics%rtotalTime)

    ! Get some temp memoty
    p_rAp => rlinsolParam%p_rsubnodeCG%rAp
    p_rr => rlinsolParam%p_rsubnodeCG%rr
    p_rp => rlinsolParam%p_rsubnodeCG%rp

    ! Initialise used vectors with zero
    call kktsp_clearControl(p_rr)
    call kktsp_clearControl(p_rp%p_rcontrolLin)
    call kktsp_clearControl(p_rAp)
    
    ! Initialization
    dalpha = 1.0_DP
    dbeta  = 1.0_DP
    dgamma = 1.0_DP
    dgammaOld = 1.0_DP
    
    brealres = rlinsolParam%p_rsubnodeCG%brealRes

    ! Create the initial defect in rd
    output_iautoOutputIndent = output_iautoOutputIndent + 2
    call newtonlin_getResidual (rlinsolParam,rkktsystemDirDeriv,rrhs,p_rr,rlocalStat)
    output_iautoOutputIndent = output_iautoOutputIndent - 2
        
    call newtonlin_sumStatistics(rlocalStat,rstatistics,NLIN_STYPE_RESCALC)
        
    call kktsp_controlCopy (p_rr,p_rp%p_rcontrolLin)

    ! Scalar product of rp.
    dgamma = kktsp_scalarProductControl(p_rr,p_rr)

    ! Apply the loop is a do loop...
    call itc_initIteration(rlinsolParam%riter)

    do while (.true.)
    
!      if ((mod(rlinsolParam%riter%niterations,5) .eq. 0) .and.&
!          (rlinsolParam%riter%niterations .ne. 0)) then
!        ! Restart
!
!        ! Initialization
!        dalpha = 1.0_DP
!        dbeta  = 1.0_DP
!        dgamma = 1.0_DP
!        dgammaOld = 1.0_DP
!        
!        ! DEBUG!!!
!        !call kktsp_controlCopy (p_rr,p_rp%p_rcontrolLin)
!        
!        ! Create the initial defect in rd
!        call output_line ("Space-time CG: Restart.")
!        !output_iautoOutputIndent = output_iautoOutputIndent + 2
!        !call newtonlin_getResidual (rlinsolParam,rkktsystemDirDeriv,rrhs,p_rr,rlocalStat)
!        !output_iautoOutputIndent = output_iautoOutputIndent - 2
!            
!        !call newtonlin_sumStatistics(rlocalStat,rstatistics,NLIN_STYPE_RESCALC)
!
!        ! DEBUG!!!
!        !call kktsp_controlCompare (p_rr,p_rp%p_rcontrolLin)
!            
!        call kktsp_controlCopy (p_rr,p_rp%p_rcontrolLin)
!
!        ! Scalar product of rp.
!        dgamma = kktsp_scalarProductControl(p_rr,p_rr)
!
!        ! Initialise used vectors with zero
!        !call kktsp_clearControl(p_rr)
!        call kktsp_clearControl(p_rAp)
!
!      end if
    
      ! -------------------------------------------------------------
      ! Norm of the residual
      ! -------------------------------------------------------------

      call kkt_controlResidualNorm (&
          rkktsystemDirDeriv%p_rkktsystem%p_roperatorAsmHier%ranalyticData,&
          p_rr,dres,rlinsolParam%rprecParameters%iresnorm)

      if (rlinsolParam%riter%cstatus .eq. ITC_STATUS_UNDEFINED) then
        ! Remember the initial residual
        call itc_initResidual(rlinsolParam%riter,dres)
      else
        ! Push the residual, increase the iteration counter
        call itc_pushResidual(rlinsolParam%riter,dres)
      end if

      if (rlinsolParam%rprecParameters%ioutputLevel .ge. 2) then
        call output_line ("Space-time CG: Iteration "// &
            trim(sys_siL(rlinsolParam%riter%niterations,10))// &
            ", ||res(u)|| = "// &
            trim(sys_sdEL(dres,10)))
      end if
      
      ! -------------------------------------------------------------
      ! Check for convergence / divergence / ...
      ! -------------------------------------------------------------
      if ((.not. brealres) .and. (rlinsolParam%riter%cstatus .eq. ITC_STATUS_CONVERGED)) then
        
        ! Apply a restart to check if we really reached the residual

        ! Create the defect in rd
        call output_line ("Space-time CG: Restart for final residual check.")
            
        output_iautoOutputIndent = output_iautoOutputIndent + 2
        call newtonlin_getResidual (rlinsolParam,rkktsystemDirDeriv,rrhs,p_rr,rlocalStat)
        call kkt_controlResidualNorm (&
            rkktsystemDirDeriv%p_rkktsystem%p_roperatorAsmHier%ranalyticData,&
            p_rr,dres,rlinsolParam%rprecParameters%iresnorm)
        output_iautoOutputIndent = output_iautoOutputIndent - 2
            
        call newtonlin_sumStatistics(rlocalStat,rstatistics,NLIN_STYPE_RESCALC)
        
        rlocaliter = rlinsolParam%riter
        call itc_resetResidualQueue (rlocaliter)
        call itc_repushResidual(rlocaliter,dres)

        if (rlinsolParam%rprecParameters%ioutputLevel .ge. 2) then
          call output_line ("Space-time CG: Iteration "// &
              trim(sys_siL(rlinsolParam%riter%niterations,10))// &
              ", ||res(u)|| = "// &
              trim(sys_sdEL(dres,10)))
        end if

        if (rlocaliter%cstatus .eq. ITC_STATUS_CONTINUE) then
          ! Restart, continue the iteration
          call output_line ("Space-time CG: Not converged. Continuing...")

          ! Restart
          dalpha = 1.0_DP
          dbeta  = 1.0_DP
          dgamma = 1.0_DP
          dgammaOld = 1.0_DP
          
          !brealres = .true.
          
          call kktsp_controlCopy (p_rr,p_rp%p_rcontrolLin)

          ! Scalar product of rp.
          dgamma = kktsp_scalarProductControl(p_rr,p_rr)

          ! Initialise used vectors with zero
          call kktsp_clearControl(p_rAp)

          ! Reset the statistics
          rlinsolParam%riter = rlocaliter

        end if

      end if
      
      if (rlinsolParam%riter%cstatus .ne. ITC_STATUS_CONTINUE) exit
      
      ! -------------------------------------------------------------
      ! Update of the solution
      ! -------------------------------------------------------------
      ! From:
      !    http://en.wikipedia.org/wiki/Conjugate_gradient_method
      ! and the Num-I script of Prof, Turek.

      ! Calculate Ap.
      output_iautoOutputIndent = output_iautoOutputIndent + 2
      call newtonlin_applyOperator (rlinsolParam,p_rp,p_rAp,rlocalStat)
      output_iautoOutputIndent = output_iautoOutputIndent - 2
      
      call newtonlin_sumStatistics(rlocalStat,rstatistics,NLIN_STYPE_RESCALC)

      ! Calculate the parameter ALPHA
      dtemp = kktsp_scalarProductControl(p_rp%p_rcontrolLin,p_rAp)
      if (abs(dtemp) .gt. SYS_MINREAL_DP) then
        dalpha = dgamma / dtemp
        !dalpha = kktsp_scalarProductControl(p_rp%p_rcontrolLin,p_rr)
        !dalpha = dalpha / dtemp
      else
        dalpha = 0.0_DP
      end if
      
      if (rlinsolParam%rprecParameters%ioutputLevel .ge. 3) then
        call output_line ("Space-time CG: <p,Ap>    = "// &
            trim(sys_sdEL(dtemp,10)))
        call output_line ("Space-time CG: GAMMA     = "// &
            trim(sys_sdEL(dgamma,10)))
        call output_line ("Space-time CG: ALPHA     = "// &
            trim(sys_sdEL(dalpha,10)))
      end if

      ! Add the residual (damped) to the current control.
      ! This updates the current control.
      call kktsp_controlLinearComb (p_rp%p_rcontrolLin,dalpha,&
          rkktsystemDirDeriv%p_rcontrolLin,1.0_DP)

      ! Update the residual
      if (.not. brealres) then
        ! Pseudo-residual
        call kktsp_controlLinearComb (p_rAp,-dalpha,p_rr,1.0_DP)
      else
        ! Real residual
        output_iautoOutputIndent = output_iautoOutputIndent + 2
        call newtonlin_getResidual (rlinsolParam,rkktsystemDirDeriv,rrhs,p_rr,rlocalStat)
        output_iautoOutputIndent = output_iautoOutputIndent - 2
      end if
          
      ! Calculate beta
      dgammaOld = dgamma
      dgamma = kktsp_scalarProductControl(p_rr,p_rr)
      !dgamma = kktsp_scalarProductControl(p_rr,p_rp%p_rcontrolLin)
      dbeta = dgamma / dgammaOld
      
      if (rlinsolParam%rprecParameters%ioutputLevel .ge. 3) then
        call output_line ("Space-time CG: GAMMA_NEW = "// &
            trim(sys_sdEL(dgamma,10)))
        call output_line ("Space-time CG: BETA      = "// &
            trim(sys_sdEL(dbeta,10)))
      end if
      
      ! Update rp.
      call kktsp_controlLinearComb (p_rr,1.0_DP,p_rp%p_rcontrolLin,dbeta)
    
      ! -------------------------------------------------------------
      ! Proceed with the iteration
      ! -------------------------------------------------------------
    
    end do
    
    call stat_stopTimer (rstatistics%rtotalTime)

    ! Statistics
    rstatistics%niterations = rstatistics%niterations + rlinsolParam%riter%niterations
    
    if (rlinsolParam%rprecParameters%ioutputLevel .ge. 2) then
      call output_lbrk()
      call output_line ("Space-time CG: Statistics.")
      call output_lbrk()
      call itc_printStatistics(rlinsolParam%riter)
    end if
    
!    ! DEBUG!!!
!    call newtonlin_getResidual (rlinsolParam,rkktsystemDirDeriv,rrhs,p_rr,rlocalStat)
!    call kkt_controlResidualNorm (&
!        rkktsystemDirDeriv%p_rkktsystem%p_roperatorAsmHier%ranalyticData,&
!        p_rr,dres,rlinsolParam%rprecParameters%iresnorm)
!    if (rlinsolParam%rprecParameters%ioutputLevel .ge. 2) then
!      call output_line ("Space-time CG: Finished "// &
!          trim(sys_siL(rlinsolParam%riter%niterations,10))// &
!          ", ||res(u)|| = "// &
!          trim(sys_sdEL(dres,10)))
!    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonlin_bicgstab (rlinsolParam,rkktsystemDirDeriv,rrhs,rstatistics)
  
!<description>
  ! Applies a BiCGStab iteration in the control space
  ! for the linearised KKT system.
!</description>

!<inputoutput>
  ! Parameters for the iteration.
  ! The output parameters are changed according to the iteration.
  type(t_linsolParameters), intent(inout), target :: rlinsolParam

  ! Structure defining the linearised KKT system.
  ! The linearised control in this structure is used as
  ! initial and target vector.
  ! The control on the maximum level receives the result.
  type(t_kktsystemDirDeriv), intent(inout), target :: rkktsystemDirDeriv
  
  ! Right-hand side of the linearised control equation
  type(t_controlSpace), intent(inout), target :: rrhs
!</inputoutput>

!<ouptut>
  ! Solver statistics. 
  type(t_newtonlinSolverStat), intent(out) :: rstatistics
!</ouptut>

!</subroutine>

  ! local variables
  real(DP) :: dalpha,dbeta,domega0,domega1,domega2,dres
  real(DP) :: drho1,drho0
  type(t_newtonlinSolverStat) :: rlocalStat
  type(t_iterationControl) :: rlocaliter

  ! Our structure
  type(t_linsolBiCGStab), pointer :: p_rsubnode
  
  ! Pointers to temporary vectors - named for easier access
  type(t_controlSpace), pointer :: p_rr0,p_rpA,p_rsA
  type(t_kktsystemDirDeriv), pointer :: p_rp,p_rr
  
    ! Getch some information
    p_rsubnode => rlinsolParam%p_rsubnodeBiCGStab

    ! Set pointers to the temporary vectors
    p_rr   => p_rsubnode%rp
    p_rp   => p_rsubnode%rr

    p_rr0  => p_rsubnode%RtempVectors(1)
    p_rpA  => p_rsubnode%RtempVectors(2)
    p_rsA  => p_rsubnode%RtempVectors(3)
    
    ! Measure the total time
    call stat_startTimer (rstatistics%rtotalTime)

    ! rrhs is our RHS. rkktsystemDirDeriv%p_rcontrolLin points to a new vector which will be our
    ! iteration vector. 
      
    ! Initialise used vectors with zero
    call kktsp_clearControl(p_rp%p_rcontrolLin)
    call kktsp_clearControl(p_rpA)
    
    ! Initialise the iteration vector with zero.

    ! Initialization
    drho0  = 1.0_DP
    dalpha = 1.0_DP
    domega0 = 1.0_DP

    ! Start the iteration
    call itc_initIteration(rlinsolParam%riter)

    do

      ! -------------------------------------------------------------
      ! Get/Save the norm of the current residual
      ! -------------------------------------------------------------
      if (rlinsolParam%riter%cstatus .eq. ITC_STATUS_UNDEFINED) then

        ! Create the initial defect in rd
        output_iautoOutputIndent = output_iautoOutputIndent + 2
        call newtonlin_getResidual (rlinsolParam,&
            rkktsystemDirDeriv,rrhs,p_rr%p_rcontrolLin,rlocalStat)
        output_iautoOutputIndent = output_iautoOutputIndent - 2
            
        call newtonlin_sumStatistics(rlocalStat,rstatistics,NLIN_STYPE_RESCALC)

        ! Get the norm of the residuum
        dres = kktsp_getNormControl ( &
            p_rr%p_rcontrolLin,rlinsolParam%rprecParameters%iresnorm)
            
        if (.not.((dres .ge. 1E-99_DP) .and. &
                  (dres .le. 1E99_DP))) dres = 0.0_DP

        ! Remember the initial residual
        call itc_initResidual(rlinsolParam%riter,dres)

        if (rlinsolParam%rprecparameters%ioutputLevel .ge. 2) then
          call output_line ("Space-time BiCGStab: Iteration "// &
              trim(sys_siL(rlinsolParam%riter%niterations,10))//&
              ",  ||res(u)|| = "//trim(sys_sdEL(dres,10)) )
        end if

        ! If this already stops the iteration, cancel here.
        if (rlinsolParam%riter%cstatus .ne. ITC_STATUS_CONTINUE) exit

        ! -------------------------------------------------------------
        ! Initialisation of the actual iteration
        ! -------------------------------------------------------------

        call kktsp_controlCopy(p_rr%p_rcontrolLin,p_rr0)

      else

        ! Push the residual, increase the iteration counter
        call itc_pushResidual(rlinsolParam%riter,dres)

        if (rlinsolParam%rprecparameters%ioutputLevel .ge. 2) then
          call output_line ("Space-time BiCGStab: Iteration "// &
              trim(sys_siL(rlinsolParam%riter%niterations,10))//&
              ",  ||res(u)|| = "//trim(sys_sdEL(dres,10)) )
        end if

        ! -------------------------------------------------------------
        ! Check for convergence / divergence / ...
        ! -------------------------------------------------------------
        if (rlinsolParam%riter%cstatus .eq. ITC_STATUS_CONVERGED) then
          
          ! Apply a restart to check if we really reached the residual

          ! Create the defect in rd
          call output_line ("Space-time BiCGStab: Restart for final residual check.")
              
          output_iautoOutputIndent = output_iautoOutputIndent + 2
          call newtonlin_getResidual (rlinsolParam,&
              rkktsystemDirDeriv,rrhs,p_rr%p_rcontrolLin,rlocalStat)
          output_iautoOutputIndent = output_iautoOutputIndent - 2
              
          call newtonlin_sumStatistics(rlocalStat,rstatistics,NLIN_STYPE_RESCALC)
          
          dres = kktsp_getNormControl ( &
              p_rr%p_rcontrolLin,rlinsolParam%rprecParameters%iresnorm)

          rlocaliter = rlinsolParam%riter
          call itc_resetResidualQueue (rlocaliter)
          call itc_repushResidual(rlocaliter,dres)

          if (rlinsolParam%rprecParameters%ioutputLevel .ge. 2) then
            call output_line ("Space-time BiCGStab: Iteration "// &
                trim(sys_siL(rlinsolParam%riter%niterations,10))// &
                ", ||res(u)|| = "// &
                trim(sys_sdEL(dres,10)))
          end if

          if (rlocaliter%cstatus .eq. ITC_STATUS_CONTINUE) then
            ! Restart, continue the iteration
            call output_line ("Space-time BiCGStab: Not converged. Continuing...")

            ! Initialise the iteration vector with zero.
            call kktsp_clearControl(p_rp%p_rcontrolLin)
            call kktsp_clearControl(p_rpA)
            
            ! Reinitialization
            drho0  = 1.0_DP
            dalpha = 1.0_DP
            domega0 = 1.0_DP

            call kktsp_controlCopy(p_rr%p_rcontrolLin,p_rr0)
            
            ! Reset the statistics
            rlinsolParam%riter = rlocaliter
            
          end if

        end if
        
        if (rlinsolParam%riter%cstatus .ne. ITC_STATUS_CONTINUE) exit
              
      end if

      ! -------------------------------------------------------------
      ! BiCGStab iteration
      ! -------------------------------------------------------------

      drho1 = kktsp_scalarProductControl (p_rr0,p_rr%p_rcontrolLin)

      if (drho0*domega0 .eq. 0.0_DP) then
        ! Should not happen
        if (rlinsolParam%rprecparameters%ioutputLevel .ge. 2) then
          call output_line ("Space-time BiCGStab: Iteration prematurely stopped! "//&
                "Correction vector is zero!")
        end if

        ! Some tuning for the output, then cancel.
        rlinsolParam%riter%cstatus = ITC_STATUS_STAGNATED
        exit
        
      end if

      dbeta=(drho1*dalpha)/(drho0*domega0)
      drho0 = drho1

      call kktsp_controlLinearComb (p_rr%p_rcontrolLin  ,1.0_DP,p_rp%p_rcontrolLin,dbeta)
      call kktsp_controlLinearComb (p_rpA ,-dbeta*domega0,p_rp%p_rcontrolLin,1.0_DP)

      output_iautoOutputIndent = output_iautoOutputIndent + 2
      call newtonlin_applyOperator (rlinsolParam,p_rp,p_rpA,rlocalStat)
      output_iautoOutputIndent = output_iautoOutputIndent - 2
      
      call newtonlin_sumStatistics(rlocalStat,rstatistics,NLIN_STYPE_RESCALC)

      dalpha = kktsp_scalarProductControl (p_rr0,p_rpA)
      
      if (abs(dalpha) .eq. 0.0_DP) then
        ! We are below machine exactness - we can not do anything more...
        ! May happen with very small problems with very few unknowns!
        if (rlinsolParam%rprecparameters%ioutputLevel .ge. 2) then
          call output_line ("Space-time BiCGStab: Convergence failed, ALPHA=0!")
        end if
        rlinsolParam%riter%cstatus = ITC_STATUS_STAGNATED
        exit
      end if
        
      dalpha = drho1/dalpha

      call kktsp_controlLinearComb (p_rpA,-dalpha,p_rr%p_rcontrolLin,1.0_DP)

      output_iautoOutputIndent = output_iautoOutputIndent + 2
      call newtonlin_applyOperator (rlinsolParam,p_rr,p_rsA,rlocalStat)
      output_iautoOutputIndent = output_iautoOutputIndent - 2
      
      call newtonlin_sumStatistics(rlocalStat,rstatistics,NLIN_STYPE_RESCALC)

      domega1 = kktsp_scalarProductControl (p_rsA,p_rr%p_rcontrolLin)
      domega2 = kktsp_scalarProductControl (p_rsA,p_rsA)
      
      if (domega1 .eq. 0.0_DP) then
        domega0 = 0.0_DP
      else
        if (domega2 .eq. 0.0_DP) then
          if (rlinsolParam%rprecparameters%ioutputLevel .ge. 2) then
            call output_line ("Space-time BiCGStab: Convergence failed: omega=0!")
          end if
          rlinsolParam%riter%cstatus = ITC_STATUS_STAGNATED
          exit
        end if
        domega0 = domega1/domega2
      end if

      call kktsp_controlLinearComb (&
          p_rp%p_rcontrolLin ,dalpha, rkktsystemDirDeriv%p_rcontrolLin,1.0_DP)
      call kktsp_controlLinearComb (&
          p_rr%p_rcontrolLin ,domega0,rkktsystemDirDeriv%p_rcontrolLin,1.0_DP)

      call kktsp_controlLinearComb (p_rsA,-domega0,p_rr%p_rcontrolLin,1.0_DP)

      ! Get the norm of the new residuum
      dres = kktsp_getNormControl (p_rr%p_rcontrolLin,rlinsolParam%rprecParameters%iresnorm)
     
      ! That is it - next iteration!
    end do

    call stat_stopTimer (rstatistics%rtotalTime)

    ! Statistics
    rstatistics%niterations = rstatistics%niterations + rlinsolParam%riter%niterations

    if (rlinsolParam%rprecparameters%ioutputLevel .ge. 2) then
      call output_lbrk()
      call output_line ("Space-time BiCGStab statistics:")
      call output_lbrk()
      call itc_printStatistics(rlinsolParam%riter)
    end if
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonlin_initMG (rlinsolMultigrid,rparlist,ssection)
  
!<description>
  ! Basic initialisation of the multigrid preconditioner.
!</description>

!<input>
  ! Parameter list that contains the data for the preconditioning.
  type(t_parlist), intent(in), target :: rparlist
  
  ! Entry Section in the parameter list containing the data of the 
  ! preconditioner.
  character(len=*), intent(in) :: ssection
!</input>

!<inputoutput>
  ! Structure to be initialised.
  type(t_linsolMultigrid), intent(out) :: rlinsolMultigrid
!</inputoutput>

!</subroutine>

    ! Read parameters from the DAT file.
    call parlst_getvalue_int (rparlist, ssection, "ccycle", &
        rlinsolMultigrid%ccycle, rlinsolMultigrid%ccycle)
    
    call parlst_getvalue_int (rparlist, ssection, "csmoother", &
        rlinsolMultigrid%csmoother, rlinsolMultigrid%csmoother)

    call parlst_getvalue_int (rparlist, ssection, "ccoarseGridSolver", &
        rlinsolMultigrid%ccoarseGridSolver, rlinsolMultigrid%ccoarseGridSolver)

    call parlst_getvalue_int (rparlist, ssection, "nsmPre", &
        rlinsolMultigrid%nsmPre, rlinsolMultigrid%nsmPre)
        
    call parlst_getvalue_int (rparlist, ssection, "nsmPost", &
        rlinsolMultigrid%nsmPost, rlinsolMultigrid%nsmPost)

    call parlst_getvalue_double (rparlist, ssection, "dcoarseGridCorrectionWeight", &
        rlinsolMultigrid%dcoarseGridCorrectionWeight, &
        rlinsolMultigrid%dcoarseGridCorrectionWeight)

    call parlst_getvalue_string (rparlist, ssection, "ssectionSmoother", &
        rlinsolMultigrid%ssectionSmoother,"SPACETIME-SMOOTHER",&
        bdequote=.true.)

    call parlst_getvalue_string (rparlist, ssection, "ssectionCoarseGridSolver", &
        rlinsolMultigrid%ssectionCoarseGridSolver,"SPACETIME-COARSEGRIDSOLVER",&
        bdequote=.true.)
        
    ! Remember the parameter list for later initialisations
    ! of the subsolvers.
    rlinsolMultigrid%p_rparList => rparList

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonlin_initExtParams (rsolver,ssection,rparamList,rparentSolver)
  
!<description>
  ! Initialises extended solver parameters according to a parameter list.
!</description>
  
!<input>
  ! Parameter list with the parameters configuring the nonlinear solver
  type(t_parlist), intent(in) :: rparamList

  ! Name of the section in the parameter list containing the parameters
  ! of the nonlinear solver.
  character(LEN=*), intent(in) :: ssection

  ! OPTIONAL: Solver structure of the main solver.
  ! If present, default parameters are taken from here.
  type(t_linsolParameters), intent(in), optionaL :: rparentSolver
!</input>

!<output>
  ! Solver structure of the subsolver to be set up
  type(t_linsolParameters), intent(inout) :: rsolver
!</output>

!</subroutine>

    ! local variables
    type(t_parlstSection), pointer :: p_rsection
    integer :: irealres

    call parlst_querysection(rparamList, ssection, p_rsection)

    if (.not. associated(p_rsection)) then
      call output_line ("Cannot create linear solver; no section """//&
          trim(ssection)//"""!", &
          OU_CLASS_ERROR,OU_MODE_STD,"newtonlin_initExtParams")
      call sys_halt()
    end if

    ! Generation of the initial condition    
    if (present(rparentSolver)) then
      rsolver%cspatialInitCondPolicy = rparentSolver%cspatialInitCondPolicy
    end if
    
    ! May be overwritten by a parameter in the section.
    call parlst_getvalue_int (p_rsection, "cspatialInitCondPolicy", &  
        rsolver%cspatialInitCondPolicy, rsolver%cspatialInitCondPolicy)

    if (rsolver%csolverType .eq. NLIN_SOLVER_CG) then
      ! CG parameters
      call parlst_getvalue_int (p_rsection, "brealres", &  
          irealres, 0)
      rsolver%p_rsubnodeCG%brealres = irealres .ne. 0
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonlin_initStrucMultigrid (rsolver,rlinsolMultigrid,rkktsystemHierarchy)
  
!<description>
  ! Initialises a multigrid substructure.
!</description>

!<input>
  ! Solver structure of the main solver.
  type(t_linsolParameters), intent(in) :: rsolver

  ! Defines the basic hierarchy of the solutions of the KKT system.
  ! This can be a "template" structure, i.e., memory for the solutions
  ! in rkktsystemHierarchy does not have to be allocated.
  type(t_kktsystemHierarchy), intent(in), target :: rkktsystemHierarchy
!</input>

!<inputoutput>
  ! Structure to be initialised.
  type(t_linsolMultigrid), intent(inout) :: rlinsolMultigrid
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: ilevel, nlevels

    nlevels = rkktsystemHierarchy%nlevels

    ! Allocate the smoothers / coarse grid solvers
    allocate (rlinsolMultigrid%p_rsubSolvers(nlevels))
    
    ! Propagate structures
    do ilevel = 1,nlevels
      rlinsolMultigrid%p_rsubSolvers(ilevel)%p_rsolverHierPrimalLin => rsolver%p_rsolverHierPrimalLin
      rlinsolMultigrid%p_rsubSolvers(ilevel)%p_rsolverHierDualLin => rsolver%p_rsolverHierDualLin
      rlinsolMultigrid%p_rsubSolvers(ilevel)%p_rkktsystemHierarchy => rsolver%p_rkktsystemHierarchy
    end do
    
    ! Level 1 is the coarse grid solver
    select case (rlinsolMultigrid%ccoarseGridSolver)
    
    ! -------------------------------
    ! Richardson solver
    ! -------------------------------
    case (1)
    
      ! Partial initialisation of the corresponding linear solver structure
      allocate(rlinsolMultigrid%p_rsubSolvers(1)%p_rsubnodeRichardson)

      rlinsolMultigrid%p_rsubSolvers(1)%csolverType = NLIN_SOLVER_RICHARDSON
      call newtonlin_initBasicParams (&
          rlinsolMultigrid%p_rsubSolvers(1)%rprecParameters,&
          rlinsolMultigrid%p_rsubSolvers(1)%riter,&
          rlinsolMultigrid%ssectionCoarseGridSolver,rlinsolMultigrid%p_rparList)
      call newtonlin_initExtParams (&
          rlinsolMultigrid%p_rsubSolvers(1),&
          rlinsolMultigrid%ssectionCoarseGridSolver,rlinsolMultigrid%p_rparList,rsolver)

    ! -------------------------------
    ! CG solver
    ! -------------------------------
    case (2)
    
      ! Partial initialisation of the corresponding linear solver structure
      allocate(rlinsolMultigrid%p_rsubSolvers(1)%p_rsubnodeCG)

      rlinsolMultigrid%p_rsubSolvers(1)%csolverType = NLIN_SOLVER_CG
      call newtonlin_initBasicParams (&
          rlinsolMultigrid%p_rsubSolvers(1)%rprecParameters,&
          rlinsolMultigrid%p_rsubSolvers(1)%riter,&
          rlinsolMultigrid%ssectionCoarseGridSolver,rlinsolMultigrid%p_rparList)
      call newtonlin_initExtParams (&
          rlinsolMultigrid%p_rsubSolvers(1),&
          rlinsolMultigrid%ssectionCoarseGridSolver,rlinsolMultigrid%p_rparList,rsolver)

    ! -------------------------------
    ! BiCGStab solver
    ! -------------------------------
    case (3)
    
      ! Partial initialisation of the corresponding linear solver structure
      allocate(rlinsolMultigrid%p_rsubSolvers(1)%p_rsubnodeBiCGStab)

      rlinsolMultigrid%p_rsubSolvers(1)%csolverType = NLIN_SOLVER_BICGSTAB
      call newtonlin_initBasicParams (&
          rlinsolMultigrid%p_rsubSolvers(1)%rprecParameters,&
          rlinsolMultigrid%p_rsubSolvers(1)%riter,&
          rlinsolMultigrid%ssectionCoarseGridSolver,rlinsolMultigrid%p_rparList)
      call newtonlin_initExtParams (&
          rlinsolMultigrid%p_rsubSolvers(1),&
          rlinsolMultigrid%ssectionCoarseGridSolver,rlinsolMultigrid%p_rparList,rsolver)

    case default
      call output_line ("Invalid coarse grid solver.", &
          OU_CLASS_ERROR,OU_MODE_STD,"newtonlin_initStrucMultigrid")
      call sys_halt()

    end select
    
    ! Level 2..nlevels are smoothers
    do ilevel = 2,nlevels

      ! Level 1 is the coarse grid solver
      select case (rlinsolMultigrid%csmoother)
      
      ! -------------------------------
      ! Richardson smoother
      ! -------------------------------
      case (1)

        ! Partial initialisation of the corresponding linear solver structure
        allocate(rlinsolMultigrid%p_rsubSolvers(ilevel)%p_rsubnodeRichardson)

        rlinsolMultigrid%p_rsubSolvers(ilevel)%csolverType = NLIN_SOLVER_RICHARDSON
        call newtonlin_initBasicParams (&
            rlinsolMultigrid%p_rsubSolvers(ilevel)%rprecParameters,&
            rlinsolMultigrid%p_rsubSolvers(ilevel)%riter,&
            rlinsolMultigrid%ssectionSmoother,rlinsolMultigrid%p_rparList)
        call newtonlin_initExtParams (&
            rlinsolMultigrid%p_rsubSolvers(ilevel),&
            rlinsolMultigrid%ssectionSmoother,rlinsolMultigrid%p_rparList,rsolver)
        
      ! -------------------------------
      ! CG smoother
      ! -------------------------------
      case (2)

        ! Partial initialisation of the corresponding linear solver structure
        allocate(rlinsolMultigrid%p_rsubSolvers(ilevel)%p_rsubnodeCG)

        rlinsolMultigrid%p_rsubSolvers(ilevel)%csolverType = NLIN_SOLVER_CG
        call newtonlin_initBasicParams (&
            rlinsolMultigrid%p_rsubSolvers(ilevel)%rprecParameters,&
            rlinsolMultigrid%p_rsubSolvers(ilevel)%riter,&
            rlinsolMultigrid%ssectionSmoother,rlinsolMultigrid%p_rparList)
        call newtonlin_initExtParams (&
            rlinsolMultigrid%p_rsubSolvers(ilevel),&
            rlinsolMultigrid%ssectionSmoother,rlinsolMultigrid%p_rparList,rsolver)

      ! -------------------------------
      ! BiCGStab smoother
      ! -------------------------------
      case (3)

        ! Partial initialisation of the corresponding linear solver structure
        allocate(rlinsolMultigrid%p_rsubSolvers(ilevel)%p_rsubnodeBiCGStab)

        rlinsolMultigrid%p_rsubSolvers(ilevel)%csolverType = NLIN_SOLVER_BICGSTAB
        call newtonlin_initBasicParams (&
            rlinsolMultigrid%p_rsubSolvers(ilevel)%rprecParameters,&
            rlinsolMultigrid%p_rsubSolvers(ilevel)%riter,&
            rlinsolMultigrid%ssectionSmoother,rlinsolMultigrid%p_rparList)
        call newtonlin_initExtParams (&
            rlinsolMultigrid%p_rsubSolvers(ilevel),&
            rlinsolMultigrid%ssectionSmoother,rlinsolMultigrid%p_rparList,rsolver)
        
      case default
        call output_line ("Invalid smoother.", &
            OU_CLASS_ERROR,OU_MODE_STD,"newtonlin_initStrucMultigrid")
        call sys_halt()

      end select
    
    end do
    
    ! Initialise the structures of the subsolvers
    do ilevel = 1,nlevels
      call newtonlin_initStructure (rlinsolMultigrid%p_rsubSolvers(ilevel),&
          rkktsystemHierarchy,rsolver%p_rprjHierSpaceTimePrimal,&
          rsolver%p_rprjHierSpaceTimeDual,rsolver%p_rprjHierSpaceTimeControl,ilevel)
    end do
    
    ! Allocate the cycle counter.
    allocate(rlinsolMultigrid%p_IcycleCount(nlevels))
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonlin_initDataMultigrid (rsolver,rlinsolMultigrid,rsolutionHierarchy)
  
!<description>
  ! Initialises a multigrid substructure.
!</description>

!<input>
  ! Solver structure of the main solver.
  type(t_linsolParameters), intent(in) :: rsolver

  ! Defines a hierarchy of the solutions of the KKT system.
  type(t_kktsystemHierarchy), intent(inout), target :: rsolutionHierarchy
!</input>

!<inputoutput>
  ! Structure to be initialised.
  type(t_linsolMultigrid), intent(inout) :: rlinsolMultigrid
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: ilevel

    ! Take a look to all levels
    do ilevel = 1,size(rlinsolMultigrid%p_rsubSolvers)

      ! Initialise the data of the subsolvers
      call newtonlin_initData (rlinsolMultigrid%p_rsubSolvers(ilevel),rsolutionHierarchy,ilevel)
      
      ! Modification of the equation
      rlinsolMultigrid%p_rsubSolvers(ilevel)%ceqnflags = rsolver%ceqnflags

    end do

    ! Create temp vectors for the iteration based on the given KKT system hierarchy
    call kkth_initHierarchyDirDeriv (rlinsolMultigrid%rsolutionHier,rsolutionHierarchy)
    call kkth_initHierarchyDirDeriv (rlinsolMultigrid%rdefectHier,rsolutionHierarchy)
    call kkth_initHierarchyDirDeriv (rlinsolMultigrid%rrhsHier,rsolutionHierarchy)
    call kkth_initHierarchyDirDeriv (rlinsolMultigrid%rtempHier,rsolutionHierarchy)
   
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonlin_doneDataMultigrid (rlinsolMultigrid)
  
!<description>
  ! Cleanup of the data initalised in newtonlin_initStructure.
!</description>

!<inputoutput>
  ! Structure to be cleaned up.
  type(t_linsolMultigrid), intent(inout) :: rlinsolMultigrid
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: ilevel

    ! Take a look to all levels
    do ilevel = 1,size(rlinsolMultigrid%p_rsubSolvers)

      call newtonlin_doneData (rlinsolMultigrid%p_rsubSolvers(ilevel))
    
    end do

    ! Release temp data
    call kkth_doneHierarchyDirDeriv (rlinsolMultigrid%rsolutionHier)
    call kkth_doneHierarchyDirDeriv (rlinsolMultigrid%rdefectHier)
    call kkth_doneHierarchyDirDeriv (rlinsolMultigrid%rrhsHier)
    call kkth_doneHierarchyDirDeriv (rlinsolMultigrid%rtempHier)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonlin_doneStrucMultigrid (rlinsolMultigrid)
  
!<description>
  ! Cleanup of the data initalised in newtonlin_initStructure.
!</description>

!<inputoutput>
  ! Structure to be cleaned up.
  type(t_linsolMultigrid), intent(inout) :: rlinsolMultigrid
!</inputoutput>

!</subroutine>

    integer :: ilevel

    ! Release subsolvers
    do ilevel = 1,size(rlinsolMultigrid%p_rsubSolvers)
    
      ! Release structural data
      call newtonlin_doneStructure (rlinsolMultigrid%p_rsubSolvers(ilevel))

      ! Final cleanup
      select case (rlinsolMultigrid%p_rsubSolvers(ilevel)%csolverType)
      
      ! -------------------------------
      ! Richardson smoother
      ! -------------------------------
      case (NLIN_SOLVER_RICHARDSON)

        deallocate(rlinsolMultigrid%p_rsubSolvers(ilevel)%p_rsubnodeRichardson)

      ! -------------------------------
      ! CG smoother
      ! -------------------------------
      case (NLIN_SOLVER_CG)

        deallocate(rlinsolMultigrid%p_rsubSolvers(ilevel)%p_rsubnodeCG)
            
      ! -------------------------------
      ! BICGSTAB smoother
      ! -------------------------------
      case (NLIN_SOLVER_BICGSTAB)

        deallocate(rlinsolMultigrid%p_rsubSolvers(ilevel)%p_rsubnodeBiCGStab)

      end select
    
    end do

    ! Release memory
    deallocate(rlinsolMultigrid%p_IcycleCount)
    deallocate(rlinsolMultigrid%p_rsubSolvers)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonlin_multigrid (rlinsolParam,rkktsysDirDerivHierarchy,rrhs,rstatistics,ilevel)
  
!<description>
  ! Applies a Multigrid iteration in the control space
  ! for the linearised KKT system.
!</description>

!<input>
  ! Level in the hierarchy.
  integer, intent(in) :: ilevel
!</input>

!<inputoutput>
  ! Parameters for the iteration.
  ! The output parameters are changed according to the iteration.
  type(t_linsolParameters), intent(inout) :: rlinsolParam

  ! Hierarchy of directional derivatives of the KKT system.
  ! The control on the maximum level receives the result.
  type(t_kktsystemDirDerivHierarchy), intent(inout) :: rkktsysDirDerivHierarchy

  ! Right-hand side of the linearised control equation
  type(t_controlSpace), intent(inout), target :: rrhs
!</inputoutput>

!<ouptut>
  ! Solver statistics. 
  type(t_newtonlinSolverStat), intent(out) :: rstatistics
!</ouptut>

!</subroutine>

    ! Call the MG driver on the maximum level.
    call newtonlin_multigridDriver (&
        rlinsolParam,ilevel,rlinsolParam%p_rkktsystemHierarchy%nlevels,&
        rkktsysDirDerivHierarchy,rrhs,rstatistics)
        
    ! Statistics
    if (rlinsolParam%rprecParameters%ioutputLevel .ge. 2) then
      call output_lbrk()
      call output_line ("Space-time Multigrid: Statistics.")
      call output_lbrk()
      call itc_printStatistics(rlinsolParam%riter)
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  recursive subroutine MG_initCycles (ccycle,ilevel,nlevels,IcycleCount)
  
!<description>
  ! Initialises the cycle counter on level 1..ilevel.
!</description>

!<input>
  ! Cycle type.
  ! =0: F-cycle, =1: V-cycle, =2: W-cycle
  integer, intent(in) :: ccycle

  ! Current level
  integer, intent(in) :: ilevel
  
  ! Maximum level
  integer, intent(in) :: nlevels
!</input>

!<inputoutput>
  ! Cycle counter.
  integer, dimension(:), intent(inout) :: IcycleCount
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i

    if (ccycle .eq. 0) then
      ! F-cycle. Cycle counter is =2 on all levels
      ! if ilevel=nlevels and =1 otherwise.
      if (ilevel .eq. nlevels) then
        do i=1,ilevel
          IcycleCount(i) = 2
        end do
      else
        do i=1,ilevel
          IcycleCount(i) = 1
        end do
      end if
      
    else
      ! V-cyclke, W-cycle,...
      do i=1,ilevel
        IcycleCount(i) = ccycle
      end do
    end if    

  end subroutine

  ! ***************************************************************************

!<subroutine>

  recursive subroutine newtonlin_multigridDriver (rlinsolParam,ilevel,nlevels,&
      rkktsysDirDerivHierarchy,rrhs,rstatistics)
  
!<description>
  ! This is the actual multigrid driver which applies the
  ! multilevel iteration.
!</description>

!<input>
  ! Current level
  integer, intent(in) :: ilevel
  
  ! Maximum level in the hierarchy
  integer, intent(in) :: nlevels
!</input>

!<inputoutput>
  ! Parameters for the iteration.
  ! The output parameters are changed according to the iteration.
  type(t_linsolParameters), intent(inout), target :: rlinsolParam

  ! Hierarchy of directional derivatives of the KKT system.
  ! The control on level ilevel receives the result.
  type(t_kktsystemDirDerivHierarchy), intent(inout) :: rkktsysDirDerivHierarchy
  
  ! Right-hand side of the linearised control equation
  type(t_controlSpace), intent(inout), target :: rrhs
!</inputoutput>

!<ouptut>
  ! Solver statistics. 
  type(t_newtonlinSolverStat), intent(out) :: rstatistics
!</ouptut>

!</subroutine>

    ! local variables
    type(t_linsolMultigrid), pointer :: p_rsubnodeMultigrid
    real(DP) :: dres,dalpha
    type(t_kktsystemDirDeriv), pointer :: p_rkktSysRhsCoarse, p_rkktSysSolCoarse, p_rkktSysSolFine
    type(t_kktsystemDirDeriv), pointer :: p_rkktSysSolution, p_rkktSysDefect, p_rkktSysTemp
    type(t_newtonlinSolverStat) :: rlocalStat
    type(t_iterationControl), pointer :: p_riter
    type(t_iterationControl), target :: riterLocal
    
!    ! DEBUG!!!
!    type(t_spacetimeOperatorAsm) :: roperatorAsm
    
    ! Get multigrid parameters.
    p_rsubnodeMultigrid => rlinsolParam%p_rsubnodeMultigrid

!    ! DEBUG!!!
!    call optcpp_spaceTimeVisualisation (rlinsolParam%p_rsettingsSolver%rpostproc,rrhs%p_rvectorAccess,&
!        CCSPACE_INTERMEDCONTROL,rlinsolParam%p_rsettingsSolver%rphysics,&
!        rlinsolParam%p_rsettingsSolver%rsettingsOptControl,10+ilevel)        
!    call sptivec_saveToFileSequence (&
!        rrhs%p_rvector,"(""ns/controllin."//trim(sys_siL(10+ilevel,10))//"."",I5.5)",.true.)
!    call kkth_getKKTsystemDirDeriv (rkktsysDirDerivHierarchy,ilevel,p_rkktSysSolution)
!    
!    
!    call optcpp_spaceTimeVisualisation (rlinsolParam%p_rsettingsSolver%rpostproc,&
!        p_rkktSysSolution%p_rkktSystem%p_rprimalSol%p_rvectorAccess,CCSPACE_PRIMAL,&
!        rlinsolParam%p_rsettingsSolver%rphysics,rlinsolParam%p_rsettingsSolver%rsettingsOptControl,10+ilevel)        
!    
!    call optcpp_spaceTimeVisualisation (rlinsolParam%p_rsettingsSolver%rpostproc,&
!        p_rkktSysSolution%p_rkktSystem%p_rdualSol%p_rvectorAccess,CCSPACE_DUAL,&
!        rlinsolParam%p_rsettingsSolver%rphysics,rlinsolParam%p_rsettingsSolver%rsettingsOptControl,10+ilevel)        
!
!    call sptivec_saveToFileSequence (&
!        p_rkktSysSolution%p_rkktSystem%p_rprimalSol%p_rvector,"(""ns/primal."//trim(sys_siL(10+ilevel,10))//"."",I5.5)",.true.)
!    call sptivec_saveToFileSequence (&
!        p_rkktSysSolution%p_rkktSystem%p_rdualSol%p_rvector,"(""ns/dual."//trim(sys_siL(10+ilevel,10))//"."",I5.5)",.true.)
!    call sptivec_saveToFileSequence (&
!        p_rkktSysSolution%p_rkktSystem%p_rintermedControl%p_rvector,"(""ns/intermedc."//trim(sys_siL(10+ilevel,10))//"."",I5.5)",.true.)
!    call sptivec_saveToFileSequence (&
!        p_rkktSysSolution%p_rkktSystem%p_rcontrol%p_rvector,"(""ns/control."//trim(sys_siL(10+ilevel,10))//"."",I5.5)",.true.)

    ! On level 1 apply the coarse grid solver
    if (ilevel .eq. 1) then

!      call sptivec_loadFromFileSequence (&
!          rrhs%p_rvector,"(""ns/controllin."//trim(sys_siL(10+ilevel,10))//"."",I5.5)",0,20,1,.true.)
!      call sptivec_scaleVector (rrhs%p_rvector,0.9_DP)
!      call sptivec_loadfromfilesequence (&
!          p_rkktsyssolution%p_rkktsystem%p_rprimalsol%p_rvector,"(""ns/primal."//trim(sys_sil(10+ilevel,10))//"."",i5.5)",0,20,1,.true.)
!      call sptivec_loadfromfilesequence (&
!          p_rkktsyssolution%p_rkktsystem%p_rdualsol%p_rvector,"(""ns/dual."//trim(sys_sil(10+ilevel,10))//"."",i5.5)",0,20,1,.true.)
!      call sptivec_loadfromfilesequence (&
!          p_rkktsyssolution%p_rkktsystem%p_rintermedcontrol%p_rvector,"(""ns/intermedc."//trim(sys_sil(10+ilevel,10))//"."",i5.5)",0,20,1,.true.)
!      call sptivec_loadfromfilesequence (&
!          p_rkktsyssolution%p_rkktsystem%p_rcontrol%p_rvector,"(""ns/control."//trim(sys_sil(10+ilevel,10))//"."",i5.5)",0,20,1,.true.)
!      
!      call optcpp_spaceTimeVisualisation (rlinsolParam%p_rsettingsSolver%rpostproc,&
!          rrhs%p_rvectorAccess,CCSPACE_INTERMEDCONTROL,rlinsolParam%p_rsettingsSolver%rphysics,&
!          rlinsolParam%p_rsettingsSolver%rsettingsOptControl,30+ilevel)
      
      if (nlevels .eq. 1) then
        if (rlinsolParam%rprecParameters%ioutputLevel .ge. 1) then
          call output_line (&
              "Space-time Multigrid: Only one level. Switching back to 1-level solver.")
        end if
      end if

      ! Solve. Do not work in-situ. We take the solution
      ! directly from rkktsysDirDerivHierarchy after the solution process
      ! is completed.
      output_iautoOutputIndent = output_iautoOutputIndent + 2
      call newtonlin_precond (p_rsubnodeMultigrid%p_rsubSolvers(ilevel),&
          rkktsysDirDerivHierarchy,rrhs,rstatistics,ilevel,binsitu=.false.)
      output_iautoOutputIndent = output_iautoOutputIndent - 2
    
      if (nlevels .eq. 1) then
        ! Return statistics data of the subsolver.
        call itc_copyStatistics (p_rsubnodeMultigrid%p_rsubSolvers(ilevel)%riter,rlinsolParam%riter)
      end if
    
!      ! DEBUG!!!
!      call kkth_getKKTsystemDirDeriv (rkktsysDirDerivHierarchy,ilevel,p_rkktSysSolution)
!      call optcpp_spaceTimeVisualisation (rlinsolParam%p_rsettingsSolver%rpostproc,&
!          p_rkktSysSolution%p_rcontrolLin%p_rvectorAccess,CCSPACE_CONTROL,&
!          rlinsolParam%p_rsettingsSolver%rphysics,rlinsolParam%p_rsettingsSolver%rsettingsOptControl,20+ilevel)
!
!      call optcpp_spaceTimeVisualisation (rlinsolParam%p_rsettingsSolver%rpostproc,&
!          p_rkktSysSolution%p_rprimalSolLin%p_rvectorAccess,CCSPACE_PRIMAL,&
!          rlinsolParam%p_rsettingsSolver%rphysics,rlinsolParam%p_rsettingsSolver%rsettingsOptControl,20+ilevel)
!      
!      call optcpp_spaceTimeVisualisation (rlinsolParam%p_rsettingsSolver%rpostproc,&
!          p_rkktSysSolution%p_rdualSolLin%p_rvectorAccess,CCSPACE_DUAL,&
!          rlinsolParam%p_rsettingsSolver%rphysics,rlinsolParam%p_rsettingsSolver%rsettingsOptControl,20+ilevel)


      ! Finish
      return
    end if
    
    ! Total time
    call stat_startTimer(rstatistics%rtotalTime)

    ! On level 2..nlevels, apply the smoothing and coarse grid correction process.

    if (ilevel .eq. nlevels) then
      ! Initialise the cycle counter for all levels.
      call MG_initCycles (&
          p_rsubnodeMultigrid%ccycle,nlevels,nlevels,p_rsubnodeMultigrid%p_IcycleCount)

      ! On level nlevels, apply as many iterations as configured.
      p_riter => rlinsolParam%riter
    else
      ! On all the other levels, apply as many iterations as configured in the cycle.
      !
      ! Get a copy of the iteration control structure and fix the number of
      ! iterations.
      riterLocal = rlinsolParam%riter
      riterLocal%nminIterations = p_rsubnodeMultigrid%p_IcycleCount(ilevel)
      riterLocal%nmaxIterations = p_rsubnodeMultigrid%p_IcycleCount(ilevel)
      p_riter => riterLocal
    end if
    
    ! Apply the multigrid iteration

    call kkth_getKKTsystemDirDeriv (rkktsysDirDerivHierarchy,ilevel,p_rkktSysSolution)
    call kkth_getKKTsystemDirDeriv (p_rsubnodeMultigrid%rrhsHier,ilevel-1,p_rkktSysRhsCoarse)
    call kkth_getKKTsystemDirDeriv (p_rsubnodeMultigrid%rsolutionHier,ilevel-1,p_rkktSysSolCoarse)
    call kkth_getKKTsystemDirDeriv (p_rsubnodeMultigrid%rsolutionHier,ilevel,p_rkktSysSolFine)
    call kkth_getKKTsystemDirDeriv (p_rsubnodeMultigrid%rdefectHier,ilevel,p_rkktSysDefect)
    call kkth_getKKTsystemDirDeriv (p_rsubnodeMultigrid%rtempHier,ilevel,p_rkktSysTemp)
    
    ! Zero initial solution
    call kkt_clearDirDeriv (p_rkktSysSolution)

    ! Start the iteration
    call itc_initIteration(p_riter)

    do
    
      if (ilevel .lt. nlevels) then
        if (p_riter%cstatus .ne. ITC_STATUS_UNDEFINED) then
          ! On level < nlevels, do not calculate the next residual but
          ! exit immediately if the maximum number of iterations
          ! (= cycle count) will be reached.
          if (p_riter%niterations+1 .ge. p_riter%nmaxIterations) exit
        end if
      end if
      
      ! Get the residual
      output_iautoOutputIndent = output_iautoOutputIndent + 2
      
      call newtonlin_getResidual (rlinsolParam,&
          p_rkktSysSolution,rrhs,&
          p_rkktSysDefect%p_rcontrolLin,rlocalStat,&
          dres,rlinsolParam%rprecParameters%iresnorm)
      
      output_iautoOutputIndent = output_iautoOutputIndent - 2
      
      call newtonlin_sumStatistics(rlocalStat,rstatistics,NLIN_STYPE_RESCALC)

      ! Some checks and output on the maximum level
      if (ilevel .eq. nlevels) then
      
        if (p_riter%cstatus .eq. ITC_STATUS_UNDEFINED) then
          ! Remember the initial residual
          call itc_initResidual(p_riter,dres)
        else
          ! Push the residual, increase the iteration counter
          call itc_pushResidual(p_riter,dres)
        end if

        if (rlinsolParam%rprecParameters%ioutputLevel .ge. 2) then
          call output_line ("Space-time MG: Iteration "// &
              trim(sys_siL(p_riter%niterations,10))// &
              ", ||res(u)|| = "// &
              trim(sys_sdEL(dres,10)))
        end if
        
        ! -------------------------------------------------------------
        ! Check for convergence / divergence / ...
        ! -------------------------------------------------------------
        if (p_riter%cstatus .ne. ITC_STATUS_CONTINUE) exit

      else        
        if (p_riter%cstatus .eq. ITC_STATUS_UNDEFINED) then
          ! Remember the initial residual
          call itc_initResidual(p_riter,dres)
        else
          ! Push the residual, increase the iteration counter
          call itc_pushResidual(p_riter,dres)
        end if

        if (rlinsolParam%rprecParameters%ioutputLevel .ge. 2) then
          call output_line ("Space-time MG: Level "//trim(sys_siL(ilevel,10))// &
              ", Step "//trim(sys_siL(p_riter%niterations,10))// &
              ", ||res(u)|| = "// &
              trim(sys_sdEL(dres,10)))
        end if

      end if
      
      ! Apply presmoothing
      if (p_rsubnodeMultigrid%nsmpre .ne. 0) then
        if (rlinsolParam%rprecParameters%ioutputLevel .ge. 3) then
          call output_line ("Space-time Multigrid: Invoking pre-smoother.")
        end if
        
        output_iautoOutputIndent = output_iautoOutputIndent + 2
        call newtonlin_smoothCorrection (p_rsubnodeMultigrid%p_rsubSolvers(ilevel),&
              p_rsubnodeMultigrid%nsmpre,&
              p_rkktSysSolution,rrhs,rlocalStat)
        output_iautoOutputIndent = output_iautoOutputIndent - 2

        ! Statistics
        call newtonlin_sumStatistics (rlocalStat,rstatistics,NLIN_STYPE_MGSUBSOLVER,ilevel,nlevels)
        
        ! Get the new residual -- it has changed due to the smoothing.
        output_iautoOutputIndent = output_iautoOutputIndent + 2
        call newtonlin_getResidual (rlinsolParam,&
            p_rkktSysSolution,rrhs,&
            p_rkktSysDefect%p_rcontrolLin,rlocalStat,&
            dres,rlinsolParam%rprecParameters%iresnorm)
        output_iautoOutputIndent = output_iautoOutputIndent - 2

        call newtonlin_sumStatistics(rlocalStat,rstatistics,NLIN_STYPE_RESCALC)
        
      end if
      
      if (rlinsolParam%rprecParameters%ioutputLevel .ge. 3) then
        call output_line (&
            "Space-time Multigrid: Switching to level "//trim(sys_siL(ilevel-1,10))//".")
      end if

      ! Restriction of the defect
      call stat_startTimer(rstatistics%rtimeProlRest)
      
      call sptipr_performRestriction (&
          rlinsolParam%p_rprjHierSpaceTimeControl,ilevel,CCSPACE_CONTROL,&
          p_rkktSysRhsCoarse%p_rcontrolLin%p_rvectorAccess, &
          p_rkktSysDefect%p_rcontrolLin%p_rvectorAccess)
      
!      ! DEBUG!!!
!      call stoh_getOpAsm_slvtlv (roperatorAsm,&
!          p_rkktSysRhsCoarse%p_rkktsystem%p_roperatorAsmHier,&
!          p_rkktSysRhsCoarse%p_rkktsystem%ispacelevel,&
!          p_rkktSysRhsCoarse%p_rkktsystem%itimelevel)
!      call smva_projectToDiv0 (roperatorAsm%p_rasmTemplates,&
!          p_rkktSysRhsCoarse%p_rdualSolLin%p_rvectorAccess,&
!          p_rkktSysRhsCoarse%p_rcontrolLin%p_rvectorAccess)
      
      call stat_stopTimer(rstatistics%rtimeProlRest)
      
      if (ilevel-1 .eq. 1) then
        if (rlinsolParam%rprecParameters%ioutputLevel .ge. 3) then
          call output_line (&
              "Space-time Multigrid: Invoking coarse-grid solver.")
        end if
      end if

      ! Coarse grid correction
      output_iautoOutputIndent = output_iautoOutputIndent + 2
      call newtonlin_multigridDriver (rlinsolParam,ilevel-1,nlevels,&
          p_rsubnodeMultigrid%rsolutionHier,&
          p_rkktSysRhsCoarse%p_rcontrolLin,rlocalStat)
      output_iautoOutputIndent = output_iautoOutputIndent - 2
      
      ! Statistics
      call newtonlin_sumStatistics (rlocalStat,rstatistics,NLIN_STYPE_MGSUBSOLVER,ilevel-1,nlevels)
      
      if (rlinsolParam%rprecParameters%ioutputLevel .ge. 3) then
        call output_line (&
            "Space-time Multigrid: Switching to level "//trim(sys_siL(ilevel,10))//".")
      end if

      ! Prolongation of the correction.
      ! Prolongate to the "defect" vector on the current level and use
      ! this for the coarse grid correction. 
      ! Do not use the "solution" vector as target of the prolongation since
      ! that may be in use due to the multigrid iteration.
      call stat_startTimer(rstatistics%rtimeProlRest)
      
      call sptipr_performProlongation (rlinsolParam%p_rprjHierSpaceTimeControl,ilevel,&
          p_rkktSysSolCoarse%p_rcontrolLin%p_rvectorAccess,&
          p_rkktSysDefect%p_rcontrolLin%p_rvectorAccess)
      
      call stat_stopTimer(rstatistics%rtimeProlRest)

      dalpha = 1.0_DP
!      call newtonlin_calcCorrEnergyMin (rlinsolParam,&
!          p_rkktSysSolution,p_rkktSysDefect,rrhs,&
!          p_rkktSysTemp%p_rcontrolLin,dalpha)
!
!      if (rlinsolParam%rprecParameters%ioutputLevel .ge. 2) then
!        call output_line (&
!            "Space-time Multigrid: Coarse grid correction weight: "//trim(sys_sdEL(dalpha,10)))
!      end if

      ! And the actual correction...
      call kktsp_controlLinearComb (&
          p_rkktSysDefect%p_rcontrolLin,&
          p_rsubnodeMultigrid%dcoarseGridCorrectionWeight,&
          p_rkktSysSolution%p_rcontrolLin,&
          dalpha)

      ! Apply postsmoothing
      if (p_rsubnodeMultigrid%nsmpost .ne. 0) then
        if (rlinsolParam%rprecParameters%ioutputLevel .ge. 3) then
          call output_line (&
              "Space-time Multigrid: Invoking post-smoother.")
        end if

        output_iautoOutputIndent = output_iautoOutputIndent + 2
        call newtonlin_smoothCorrection (p_rsubnodeMultigrid%p_rsubSolvers(ilevel),&
              p_rsubnodeMultigrid%nsmpost,&
              p_rkktSysSolution,rrhs,rlocalStat)
        output_iautoOutputIndent = output_iautoOutputIndent - 2
              
        ! Statistics
        call newtonlin_sumStatistics (rlocalStat,rstatistics,NLIN_STYPE_MGSUBSOLVER,ilevel,nlevels)

        ! Calculation of the new residual is done in the beginning of the loop.
      end if

      ! Next iteration
      if (ilevel .eq. nlevels) then
        ! Re-initialise the cycles
        call MG_initCycles (&
            p_rsubnodeMultigrid%ccycle,ilevel,nlevels,p_rsubnodeMultigrid%p_IcycleCount)
      end if

    end do
    
    rstatistics%niterations = rstatistics%niterations + p_riter%niterations

    ! Initialise the cycle counter for the next run.
    call MG_initCycles (&
        p_rsubnodeMultigrid%ccycle,ilevel,nlevels,p_rsubnodeMultigrid%p_IcycleCount)
        
    call stat_stopTimer(rstatistics%rtotalTime)

!    ! DEBUG!!!
!    call kkth_getKKTsystemDirDeriv (rkktsysDirDerivHierarchy,ilevel,p_rkktSysSolution)
!    call optcpp_spaceTimeVisualisation (rlinsolParam%p_rsettingsSolver%rpostproc,&
!        p_rkktSysSolution%p_rcontrolLin%p_rvectorAccess,CCSPACE_CONTROL,&
!        rlinsolParam%p_rsettingsSolver%rphysics,rlinsolParam%p_rsettingsSolver%rsettingsOptControl,10+ilevel)  
        
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonlin_calcCorrEnergyMin (rlinsolParam,&
      rvector,rcorrVector,rrhs,rtempVector,dalpha,rstatistics)

!<description>
  ! This routine calculates the optimal coarse grid correction parameter
  ! alpha adaptively, using the energy minimisation formula.
!</description>

!<inputoutput>
  ! Parameters for the iteration.
  ! The output parameters are changed according to the iteration.
  type(t_linsolParameters), intent(inout), target :: rlinsolParam

  ! The current (uncorrected) solution vector
  type(t_kktsystemDirDeriv), intent(inout) :: rvector

  ! The correction vector which war calculated with the coarse grid
  ! and is to be added to the solution vector:
  !   x = x + alpha * correction
  ! (with correction=<tex>$P^{-1}(b-Ax)$</tex> and <tex>$P^{-1}=$</tex> multigrid on the
  ! coarse level.
  type(t_kktsystemDirDeriv), intent(inout) :: rcorrVector

  ! The current RHS vector
  type(t_controlSpace), intent(inout) :: rrhs

!</input>

!<inputoutput>
  ! A temporary vector
  type(t_controlSpace), intent(inout) :: rtempVector
!</inputoutput>

!<output>
  ! The optimal correction parameter <tex>$\alpha$</tex> for the coarse grid correction.
  real(DP), intent(out) :: dalpha

  ! Solver statistics. Statistics of the residual calculation are added to this.
  type(t_newtonlinSolverStat), intent(out), optional :: rstatistics
!</output>

!</subroutine>

    ! local variables
    type(t_newtonlinSolverStat) :: rlocatStat
    real(DP) :: a,b

    ! We calculate the optimal alpha by energy minimisation, i.e.
    ! (c.f. p. 206 in Turek`s book):
    !
    !             ( f_k - A_k x_k  ,  corr_k )
    ! alpha_k := -------------------------------------
    !             ( A_k corr_k     ,  corr_k )
    !
    ! For this purpose, we need the temporary vector.
    !
    ! Calculate nominator of the fraction

    call newtonlin_getResidual (rlinsolParam,rvector,rrhs,rtempVector,rlocatStat)
    a = kktsp_scalarProductControl(rtempVector,rcorrVector%p_rcontrolLin)

    call newtonlin_applyOperator (rlinsolParam,rcorrVector,rtempVector,rlocatStat)
    b = kktsp_scalarProductControl(rtempVector,rcorrVector%p_rcontrolLin)

    ! Return the alpha.
    if (b .ne. 0.0_DP) then
      dalpha = a/b
    else
      dalpha = 1.0_DP
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonlin_calcCorrDefMin (rlinsolParam,&
      rvector,rcorrVector,rrhs,rtempVector,rtempVector2,dalpha,rstatistics)

!<description>
  ! This routine calculates the optimal coarse grid correction parameter
  ! alpha adaptively, using the energy minimisation formula.
!</description>

!<inputoutput>
  ! Parameters for the iteration.
  ! The output parameters are changed according to the iteration.
  type(t_linsolParameters), intent(inout), target :: rlinsolParam

  ! The current (uncorrected) solution vector
  type(t_kktsystemDirDeriv), intent(inout) :: rvector

  ! The correction vector which war calculated with the coarse grid
  ! and is to be added to the solution vector:
  !   x = x + alpha * correction
  ! (with correction=<tex>$P^{-1}(b-Ax)$</tex> and <tex>$P^{-1}=$</tex> multigrid on the
  ! coarse level.
  type(t_kktsystemDirDeriv), intent(inout) :: rcorrVector

  ! The current RHS vector
  type(t_controlSpace), intent(inout) :: rrhs

!</input>

!<inputoutput>
  ! A temporary vector
  type(t_controlSpace), intent(inout) :: rtempVector

  ! A second temporary vector
  type(t_controlSpace), intent(inout) :: rtempVector2
!</inputoutput>

!<output>
  ! The optimal correction parameter <tex>$\alpha$</tex> for the coarse grid correction.
  real(DP), intent(out) :: dalpha

  ! Solver statistics. Statistics of the residual calculation are added to this.
  type(t_newtonlinSolverStat), intent(out) :: rstatistics
!</output>

!</subroutine>

    ! local variables
    real(DP) :: a,b
    type(t_newtonlinSolverStat) :: rlocatStat

    ! We calculate the optimal alpha by energy minimisation, i.e.
    ! (c.f. p. 206 in Turek`s book):
    !
    !             ( f_k - A_k x_k  ,  A_k corr_k )
    ! alpha_k := -------------------------------------
    !             ( A_k corr_k     ,  A_k corr_k )
    !
    ! For this purpose, we need the temporary vector.
    !
    ! Calculate nominator of the fraction

    call newtonlin_getResidual (rlinsolParam,rvector,rrhs,rtempVector,rlocatStat)
    call newtonlin_applyOperator (rlinsolParam,rcorrVector,rtempVector2,rlocatStat)

    ! Calculate the nominator / demoninator of the fraction
    a = kktsp_scalarProductControl(rtempVector,rtempVector2)
    b = kktsp_scalarProductControl(rtempVector2,rtempVector2)

    ! Return the alpha.
    if (b .ne. 0.0_DP) then
      dalpha = a/b
    else
      dalpha = 1.0_DP
    end if

  end subroutine
  
! ***************************************************************************

!<subroutine>

  subroutine smva_projectSpaceToDiv0 (rasmTemplates,rdiscrDual,rvector)
  
!<description>
  ! Projects a vector into the divergence free subspace.
!</description>

!<input>
  ! Assembly templates
  type(t_staticSpaceAsmTemplates), intent(in) :: rasmTemplates
  
  ! Block discretisation of the dual space
  type(t_blockDiscretisation), intent(in) :: rdiscrDual
!</input>

!<inputoutput>
  ! Vector to be projected
  type(t_vectorBlock), intent(inout) :: rvector
!</inputoutput>

!</subroutine>
  
    ! local variables
    type(t_vectorBlock) :: rrhs,rtemp,rsol
    type(t_matrixBlock) :: rmatrix
    type(t_linsolNode), pointer :: p_rsolverNode
    integer :: ierror
    real(DP), dimension(:), pointer :: p_Dx,p_Db
    integer, dimension(1) :: Irows
    
    ! Set up the RHS
    call lsysbl_createVectorBlock (rdiscrDual,rrhs)
    call lsysbl_createVectorBlock (rdiscrDual,rtemp)
    call lsysbl_createVectorBlock (rdiscrDual,rsol)
    
    ! DEBUG!!!
    call lsysbl_getbase_double (rrhs,p_Db)
    call lsysbl_getbase_double (rsol,p_Dx)
    
    ! Matrix-vector multiplication with the mass matix
    call lsyssc_scalarMatVec (rasmTemplates%rmatrixMass,&
        rvector%RvectorBlock(1),rrhs%RvectorBlock(1),1.0_DP,0.0_DP)
    call lsyssc_scalarMatVec (rasmTemplates%rmatrixMass,&
        rvector%RvectorBlock(2),rrhs%RvectorBlock(2),1.0_DP,0.0_DP)
    call lsyssc_copyVector (rvector%RvectorBlock(1),rsol%RvectorBlock(1))
    call lsyssc_copyVector (rvector%RvectorBlock(2),rsol%RvectorBlock(2))
    call lsyssc_clearVector (rsol%RvectorBlock(3))
    call lsyssc_clearVector (rrhs%RvectorBlock(3))

    ! Set up a block matrix
    call lsysbl_createEmptyMatrix (rmatrix,3,3)
    call lsyssc_duplicateMatrix (rasmTemplates%rmatrixMass,rmatrix%RmatrixBlock(1,1),&
        LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
    call lsyssc_duplicateMatrix (rasmTemplates%rmatrixMass,rmatrix%RmatrixBlock(2,2),&
        LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
    call lsyssc_duplicateMatrix (rasmTemplates%rmatrixB1,rmatrix%RmatrixBlock(1,3),&
        LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
    call lsyssc_duplicateMatrix (rasmTemplates%rmatrixB2,rmatrix%RmatrixBlock(2,3),&
        LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
    call lsyssc_duplicateMatrix (rasmTemplates%rmatrixD1,rmatrix%RmatrixBlock(3,1),&
        LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
    call lsyssc_duplicateMatrix (rasmTemplates%rmatrixD2,rmatrix%RmatrixBlock(3,2),&
        LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
    call lsyssc_duplicateMatrix (rasmTemplates%rmatrixMassPressureExtStruc,rmatrix%RmatrixBlock(3,3),&
        LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
    call lsyssc_clearMatrix (rmatrix%RmatrixBlock(3,3))
    Irows = (/1/)
    call mmod_replaceLinesByZero (rmatrix%RmatrixBlock(3,1),Irows)
    call mmod_replaceLinesByZero (rmatrix%RmatrixBlock(3,2),Irows)
    call mmod_replaceLineByLumpedMass (rmatrix%RmatrixBlock(3,3),1,&
        rasmTemplates%rmatrixMassPressureLumpInt)
    call lsysbl_updateMatStrucInfo (rmatrix)
    
    ! Solve
    call linsol_initUMFPACK4 (p_rsolverNode)
    call linsol_setMatrix (p_rsolverNode,rmatrix)
    call linsol_initStructure (p_rsolverNode, ierror)
    call linsol_initData (p_rsolverNode, ierror)
    call linsol_solveAdaptively (p_rsolverNode, rsol,rrhs,rtemp)
    call linsol_doneData (p_rsolverNode)
    call linsol_doneStructure (p_rsolverNode)
    call linsol_releaseSolver (p_rsolverNode)

    ! Get the solution    
    call lsyssc_copyVector (rsol%RvectorBlock(1),rvector%RvectorBlock(1))
    call lsyssc_copyVector (rsol%RvectorBlock(2),rvector%RvectorBlock(2))

    ! Clean up
    call lsysbl_releaseMatrix (rmatrix)
    call lsysbl_releaseVector (rrhs)
    call lsysbl_releaseVector (rtemp)
    call lsysbl_releaseVector (rsol)

  end subroutine

! ***************************************************************************

!<subroutine>

  subroutine smva_projectToDiv0 (rasmTemplates,rdualSol,rcontrol)
  
!<description>
  ! Projects a control vector into the divergence free subspace.
!</description>

!<input>
  ! Assembly templates
  type(t_staticSpaceAsmTemplates), intent(in) :: rasmTemplates

  ! Dual solution, corresponding to the control
  type(t_spaceTimeVectorAccess), intent(in) :: rdualSol
!</input>

!<inputoutput>
  ! Control vector to be projected
  type(t_spaceTimeVectorAccess), intent(inout) :: rcontrol
!</inputoutput>

!</subroutine>
  
    ! local variables
    integer :: istep
    type(t_vectorBlock), pointer :: p_rvector
    
    do istep=1,rcontrol%p_rspaceTimeVector%NEQtime
      ! Get the subvector
      call sptivec_getVectorFromPool (rcontrol, istep, p_rvector)
      
      ! Project into the divergence free space
      call smva_projectSpaceToDiv0 (rasmTemplates,rdualSol%p_rspaceDiscr,p_rvector)
      
      ! Save back
      call sptivec_commitVecInPool (rcontrol, istep)
      
    end do
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  recursive subroutine newtonlin_precond (&
      rlinsolParam,rkktsysDirDerivHierarchy,rnewtonDir,rstatistics,ilevel,binsitu)
  
!<description>
  ! Calculates the Newton search direction by applying the Newton
  ! preconditioner to the given residual rnewtonDir.
  ! This calculates a correction g with
  !      J''(u) g = d
  ! where d=rnewtonDir. On return of this routine, there is rnewtonDir=g.
!</description>

!<input>
  ! Level in the hierarchy. If not present, the maximum level is assumed.
  integer, intent(in), optional :: ilevel
!</input>

!<inputoutput>
  ! Parameters for the iteration.
  ! The output parameters are changed according to the iteration.
  type(t_linsolParameters), intent(inout) :: rlinsolParam

  ! Hierarchy of directional derivatives of the KKT system.
  ! Used for temporary calculations.
  ! The solution at level ilevel receives the result.
  type(t_kktsystemDirDerivHierarchy), intent(inout) :: rkktsysDirDerivHierarchy

  ! This structure contains the search direction to be
  ! preconditioned. Will be replaced by the precoditioned search direction
  ! if binSitu=true.
  type(t_controlSpace), intent(inout) :: rnewtonDir
  
  ! OPTIONAL: Whether to work in-situ. If set to TRUE(default),
  ! rnewtonDir is replaced by the result.
  logical, intent(in), optional :: binsitu
!</inputoutput>

!<output>
  ! Solver statistics.
  type(t_newtonlinSolverStat), intent(out) :: rstatistics
!</output>

!</subroutine>

    integer :: ilv
    type(t_kktsystemDirDeriv), pointer :: p_rkktsysDirDeriv
    type(t_timer) :: rtimer
    
    call stat_clearTimer (rtimer)
    call stat_startTimer (rtimer)
    
    ! The calculation is done on the topmost level
    ilv = rkktsysDirDerivHierarchy%nlevels
    if (present(ilevel)) ilv = ilevel
    
    p_rkktsysDirDeriv => rkktsysDirDerivHierarchy%p_RkktSysDirDeriv(ilv)

    ! For the calculation of the Newon search direction, we have to solve
    ! the linear system
    !
    !    J''(u) g = d
    !
    ! On the entry of this routine, there is d=rnewtonDir. When the routine
    ! is left, rnewtonDir receives g.
    !
    ! For this purpose, we need a linear solver based on J''(u) in the
    ! control space. For distributed control, this is multigrid-based,
    ! for boundary control (or general control with a limited number of
    ! control variables) a single-grid solver.
    !
    ! The system to solve is rather similar to the nonlinear system,
    ! but involves a couple of additional terms. It is still an iteration
    ! in u and thus, involves a forward and a backward solve for every 
    ! update -- probably on different refinement levels.
    !
    ! In the following, rnewtonDir is used as right-hand side of the
    ! equation. The solution is calculated in 
    !     rkktsystemDirDeriv%p_rcontrolLin
    ! which starts with zero.
    
    call kkt_clearDirDeriv (p_rkktsysDirDeriv)

    select case (rlinsolParam%csolverType)
    
    ! --------------------------------------
    ! Richardson iteration
    ! --------------------------------------
    case (NLIN_SOLVER_RICHARDSON)
      ! Call the Richardson iteration to calculate an update
      ! for the control.
      call newtonlin_richardson (rlinsolParam,p_rkktsysDirDeriv,rnewtonDir,rstatistics)
          
    ! --------------------------------------
    ! CG iteration
    ! --------------------------------------
    case (NLIN_SOLVER_CG)
      ! Call the Richardson iteration to calculate an update
      ! for the control.
      call newtonlin_cg (rlinsolParam,p_rkktsysDirDeriv,rnewtonDir,rstatistics)
          
    ! --------------------------------------
    ! BiCGStab iteration
    ! --------------------------------------
    case (NLIN_SOLVER_BICGSTAB)
      ! Call the Richardson iteration to calculate an update
      ! for the control.
      call newtonlin_bicgstab (rlinsolParam,p_rkktsysDirDeriv,rnewtonDir,rstatistics)
          
    ! --------------------------------------
    ! Multigrid iteration
    ! --------------------------------------
    case (NLIN_SOLVER_MULTIGRID)
      ! Invoke multigrid.
      call newtonlin_multigrid (rlinsolParam,rkktsysDirDerivHierarchy,rnewtonDir,rstatistics,ilv)
      
    end select
        
    ! Overwrite the rnewtonDir with the update.
    if (present(binsitu)) then
      if (binsitu) then
        call kktsp_controlLinearComb (&
            p_rkktsysDirDeriv%p_rcontrolLin,1.0_DP,rnewtonDir,0.0_DP)
      end if
    else
      call kktsp_controlLinearComb (&
          p_rkktsysDirDeriv%p_rcontrolLin,1.0_DP,rnewtonDir,0.0_DP)
    end if
    
    ! Add the time which this routine needed.    
    call stat_stopTimer (rtimer)
    call stat_subTimers (rstatistics%rtotalTime,rtimer)
    
    call stat_addTimers (rtimer,rstatistics%rtotalTime)
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonlin_init (rlinsolParam,rsettingsSolver,&
      rsolverHierPrimalLin,rsolverHierDualLin,rparlist,ssection)
  
!<description>
  ! Basic initialisation of the preconditioner.
!</description>

!<input>
  ! Parameters of the OptFlow solver
  type(t_settings_optflow), intent(in), target :: rsettingsSolver

  ! Hierarchy of solvers in space for all levels.
  ! Linearised primal equation.
  type(t_spaceSolverHierarchy), target :: rsolverHierPrimalLin
  
  ! Hierarchy of solvers in space for all levels.
  ! Linearised dual equation.
  type(t_spaceSolverHierarchy), target :: rsolverHierDualLin

  ! Parameter list that contains the data for the preconditioning.
  type(t_parlist), intent(in) :: rparlist
  
  ! Entry Section in the parameter list containing the data of the 
  ! preconditioner.
  character(len=*), intent(in) :: ssection
!</input>

!<inputoutput>
  ! Structure to be initialised.
  type(t_linsolParameters), intent(out) :: rlinsolParam
!</inputoutput>

!</subroutine>

    character(len=SYS_STRLEN) :: ssectionSpaceTimeMG
   
    ! Remember the settings for later use
    rlinsolParam%p_rsettingsSolver => rsettingsSolver
    rlinsolParam%p_rsolverHierPrimalLin => rsolverHierPrimalLin
    rlinsolParam%p_rsolverHierDualLin => rsolverHierDualLin
   
    ! Get the corresponding parameters from the parameter list.
    call newtonlin_initBasicParams (&
        rlinsolParam%rprecParameters,rlinsolParam%riter,ssection,rparlist)
    call newtonlin_initExtParams (rlinsolParam,ssection,rparlist)

    ! Solver type
    call parlst_getvalue_int (rparlist, ssection, "csolverType", &
        rlinsolParam%csolverType, rlinsolParam%csolverType)

    select case (rlinsolParam%csolverType)
    
    ! -----------------------------------------------------
    ! Multigrid inítialisation
    ! -----------------------------------------------------
    case (NLIN_SOLVER_MULTIGRID)
    
      ! Get the section configuring the space-time MG
      call parlst_getvalue_string (rparlist, ssection, &
          "ssectionSpaceTimeMG", ssectionSpaceTimeMG, "SPACETIME-LINSOLVER",bdequote=.true.)
          
      ! Initialise the MG parameters
      allocate(rlinsolParam%p_rsubnodeMultigrid)
      call newtonlin_initMG (rlinsolParam%p_rsubnodeMultigrid,rparlist,ssectionSpaceTimeMG)

    ! -----------------------------------------------------
    ! Richardson inítialisation
    ! -----------------------------------------------------
    case (NLIN_SOLVER_RICHARDSON)
    
      ! Initialise the MG parameters
      allocate(rlinsolParam%p_rsubnodeRichardson)

    ! -----------------------------------------------------
    ! CG inítialisation
    ! -----------------------------------------------------
    case (NLIN_SOLVER_CG)
    
      ! Initialise the MG parameters
      allocate(rlinsolParam%p_rsubnodeCG)
    
    ! -----------------------------------------------------
    ! BiCGStab inítialisation
    ! -----------------------------------------------------
    case (NLIN_SOLVER_BICGSTAB)
    
      ! Initialise the MG parameters
      allocate(rlinsolParam%p_rsubnodeBiCGStab)
    
    end select

  end subroutine

  ! ***************************************************************************

!<subroutine>

  recursive subroutine newtonlin_initStructure (rlinsolParam,rkktsystemHierarchy,&
      rprjHierSpaceTimePrimal,rprjHierSpaceTimeDual,rprjHierSpaceTimeControl,ilevel)
  
!<description>
  ! Structural initialisation of the preconditioner
!</description>

!<input>
  ! Defines the basic hierarchy of the solutions of the KKT system.
  ! This can be a "template" structure, i.e., memory for the solutions
  ! in rkktsystemHierarchy does not have to be allocated.
  type(t_kktsystemHierarchy), intent(in), target :: rkktsystemHierarchy

  ! Projection hierarchy for the interlevel projection in space/time, primal space.
  type(t_sptiProjHierarchyBlock), intent(in), target :: rprjHierSpaceTimePrimal

  ! Projection hierarchy for the interlevel projection in space/time, dual space.
  type(t_sptiProjHierarchyBlock), intent(in), target :: rprjHierSpaceTimeDual

  ! Projection hierarchy for the interlevel projection in space/time, control space.
  type(t_sptiProjHierarchyBlock), intent(in), target :: rprjHierSpaceTimeControl
  
  ! OPTIONAL: Level in the hierarchy. If not specified, the maximum level is assumed.
  integer, intent(in), optional :: ilevel
!</input>

!<inputoutput>
  ! Structure to be initialised.
  type(t_linsolParameters), intent(inout) :: rlinsolParam
!</inputoutput>

!</subroutine>

    integer :: ilev,i
   
    ilev = rkktsystemHierarchy%nlevels
    if (present(ilevel)) ilev = ilevel

    ! Remember the structure of the solutions.
    rlinsolParam%p_rkktsystemHierarchy => rkktsystemHierarchy
    
    ! Save the projection hierarchies which allows us to apply
    ! prolongation/restriction.
    rlinsolParam%p_rprjHierSpaceTimePrimal => rprjHierSpaceTimePrimal
    rlinsolParam%p_rprjHierSpaceTimeDual => rprjHierSpaceTimeDual
    rlinsolParam%p_rprjHierSpaceTimeControl => rprjHierSpaceTimeControl

    ! Initialise subsolvers
    select case (rlinsolParam%csolverType)
    
    ! -----------------------------------------------------
    ! Multigrid inítialisation
    ! -----------------------------------------------------
    case (NLIN_SOLVER_MULTIGRID)

      call newtonlin_initStrucMultigrid (&
          rlinsolParam,rlinsolParam%p_rsubnodeMultigrid,rkktsystemHierarchy)
      
    ! -----------------------------------------------------
    ! Richardson inítialisation
    ! -----------------------------------------------------
    case (NLIN_SOLVER_RICHARDSON)
    
      call kkth_initControl (&
          rlinsolParam%p_rsubnodeRichardson%rtempVector,rkktsystemHierarchy,ilev)

    ! -----------------------------------------------------
    ! CG inítialisation
    ! -----------------------------------------------------
    case (NLIN_SOLVER_CG)
    
      call kkth_initControl (&
          rlinsolParam%p_rsubnodeCG%rr,rkktsystemHierarchy,ilev)
      call kkth_initControl (&
          rlinsolParam%p_rsubnodeCG%rAp,rkktsystemHierarchy,ilev)

    ! -----------------------------------------------------
    ! BiCGStab inítialisation
    ! -----------------------------------------------------
    case (NLIN_SOLVER_BICGSTAB)
    
      do i=1,size(rlinsolParam%p_rsubnodeBiCGStab%RtempVectors)
        call kkth_initControl (&
            rlinsolParam%p_rsubnodeBiCGStab%RtempVectors(i),rkktsystemHierarchy,ilev)
      end do

    end select
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonlin_initNonlinearData (rlinsolParam,rsolutionHierarchy)
  
!<description>
  ! Propagates the nonlinearity from the highest to all levels in the hierarchy.
!</description>

!<inputoutput>
  ! Defines a hierarchy of the solutions of the KKT system.
  ! On the highest level, the solution is expected.
  ! The solution is propagated to all lower levels.
  type(t_kktsystemHierarchy), intent(inout), target :: rsolutionHierarchy

  ! Parameters of the linera solver.
  type(t_linsolParameters), intent(inout) :: rlinsolParam
!</inputoutput>

!</subroutine>

    integer :: ilevel
    type(t_kktsystem), pointer :: p_rkktSystem,p_rkktSystemCoarse

    ! Propagate the solution to the lower levels.
    do ilevel = rsolutionHierarchy%nlevels,2,-1

      call kkth_getKKTsystem (rsolutionHierarchy,ilevel-1,p_rkktSystemCoarse)
      call kkth_getKKTsystem (rsolutionHierarchy,ilevel,p_rkktSystem)
      
      call sptipr_performInterpolation (rlinsolParam%p_RprjHierSpaceTimePrimal,ilevel,&
          p_rkktSystemCoarse%p_rprimalSol%p_rvectorAccess, &
          p_rkktSystem%p_rprimalSol%p_rvectorAccess)

      call sptipr_performInterpolation (rlinsolParam%p_RprjHierSpaceTimeDual,ilevel,&
          p_rkktSystemCoarse%p_rdualSol%p_rvectorAccess, &
          p_rkktSystem%p_rdualSol%p_rvectorAccess)

      call sptipr_performInterpolation (rlinsolParam%p_RprjHierSpaceTimeControl,ilevel,&
          p_rkktSystemCoarse%p_rcontrol%p_rvectorAccess, &
          p_rkktSystem%p_rcontrol%p_rvectorAccess)

      call sptipr_performInterpolation (rlinsolParam%p_RprjHierSpaceTimeControl,ilevel,&
          p_rkktSystemCoarse%p_rintermedControl%p_rvectorAccess, &
          p_rkktSystem%p_rintermedControl%p_rvectorAccess)
    
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  recursive subroutine newtonlin_initData (rlinsolParam,rsolutionHierarchy,ilevel)
  
!<description>
  ! Final preparation of the preconditioner for preconditioning.
!</description>

!<inputoutput>
  ! Defines a hierarchy of the solutions of the KKT system.
  ! On the highest level, the solution is expected.
  ! The solution is propagated to all lower levels.
  type(t_kktsystemHierarchy), intent(inout), target :: rsolutionHierarchy

  ! Structure to be initialised.
  type(t_linsolParameters), intent(inout) :: rlinsolParam
  
  ! OPTIONAL: Level in the hierarchy. If not specified, the maximum level is assumed.
  integer, intent(in), optional :: ilevel
!</inputoutput>

!</subroutine>

    type(t_kktsystem), pointer :: p_rkktSystem
    integer :: ilev
   
    ilev = rsolutionHierarchy%nlevels
    if (present(ilevel)) ilev = ilevel
   
    ! Initialise subsolvers
    select case (rlinsolParam%csolverType)
    
    ! -----------------------------------------------------
    ! Multigrid inítialisation
    ! -----------------------------------------------------
    case (NLIN_SOLVER_MULTIGRID)

      call newtonlin_initDataMultigrid (&
          rlinsolParam,rlinsolParam%p_rsubnodeMultigrid,rsolutionHierarchy)
      
    ! -------------------------------
    ! CG smoother/solver
    ! -------------------------------
    case (NLIN_SOLVER_CG)

      ! Get the KKT system solution on that level
      call kkth_getKKTsystem (&
          rsolutionHierarchy,ilev,p_rkktsystem)

      ! Initialise the derivative in that point
      call kkt_initKKTsystemDirDeriv (&
          rlinsolParam%p_rsubnodeCG%rp,p_rkktsystem)

    ! -------------------------------
    ! BiCGStab smoother/solver
    ! -------------------------------
    case (NLIN_SOLVER_BICGSTAB)

      ! Get the KKT system solution on that level
      call kkth_getKKTsystem (&
          rsolutionHierarchy,ilev,p_rkktsystem)

      ! Initialise the derivative in that point
      call kkt_initKKTsystemDirDeriv (&
          rlinsolParam%p_rsubnodeBiCGStab%rp,p_rkktsystem)
      call kkt_initKKTsystemDirDeriv (&
          rlinsolParam%p_rsubnodeBiCGStab%rr,p_rkktsystem)
      
    end select
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  recursive subroutine newtonlin_doneData (rlinsolParam)
  
!<description>
  ! Cleanup of the data initalised in newtonlin_initData.
!</description>

!<inputoutput>
  ! Structure to be cleaned up.
  type(t_linsolParameters), intent(inout) :: rlinsolParam
!</inputoutput>

!</subroutine>

    ! Clean up subsolvers
    select case (rlinsolParam%csolverType)
    
    ! -----------------------------------------------------
    ! Multigrid
    ! -----------------------------------------------------
    case (NLIN_SOLVER_MULTIGRID)

      call newtonlin_doneDataMultigrid (rlinsolParam%p_rsubnodeMultigrid)
      
    ! -------------------------------
    ! CG smoother/solver
    ! -------------------------------
    case (NLIN_SOLVER_CG)

      ! Initialise the derivative in that point
      call kkt_doneKKTsystemDirDeriv (rlinsolParam%p_rsubnodeCG%rp)

    ! -------------------------------
    ! BiCGStab smoother/solver
    ! -------------------------------
    case (NLIN_SOLVER_BICGSTAB)

      ! Initialise the derivative in that point
      call kkt_doneKKTsystemDirDeriv (rlinsolParam%p_rsubnodeBiCGStab%rr)
      call kkt_doneKKTsystemDirDeriv (rlinsolParam%p_rsubnodeBiCGStab%rp)

    end select

  end subroutine

  ! ***************************************************************************

!<subroutine>

  recursive subroutine newtonlin_doneStructure (rlinsolParam)
  
!<description>
  ! Cleanup of the data initalised in newtonlin_initStructure.
!</description>

!<inputoutput>
  ! Structure to be cleaned up.
  type(t_linsolParameters), intent(inout) :: rlinsolParam
!</inputoutput>

!</subroutine>

    integer :: i

    ! Clean up subsolvers
    select case (rlinsolParam%csolverType)
    
    ! -----------------------------------------------------
    ! Multigrid
    ! -----------------------------------------------------
    case (NLIN_SOLVER_MULTIGRID)

      call newtonlin_doneStrucMultigrid (rlinsolParam%p_rsubnodeMultigrid)
      
    ! -----------------------------------------------------
    ! Richardson inítialisation
    ! -----------------------------------------------------
    case (NLIN_SOLVER_RICHARDSON)
    
      call kktsp_doneControlVector (rlinsolParam%p_rsubnodeRichardson%rtempVector)

    ! -----------------------------------------------------
    ! CG inítialisation
    ! -----------------------------------------------------
    case (NLIN_SOLVER_CG)
    
      call kktsp_doneControlVector (rlinsolParam%p_rsubnodeCG%rAp)
      call kktsp_doneControlVector (rlinsolParam%p_rsubnodeCG%rr)

    ! -----------------------------------------------------
    ! BiCGStab inítialisation
    ! -----------------------------------------------------
    case (NLIN_SOLVER_BICGSTAB)
    
      do i=1,size(rlinsolParam%p_rsubnodeBiCGStab%RtempVectors)
        call kktsp_doneControlVector (rlinsolParam%p_rsubnodeBiCGStab%RtempVectors(i))
      end do

    end select

  end subroutine

  ! ***************************************************************************

!<subroutine>

  recursive subroutine newtonlin_done (rlinsolParam)
  
!<description>
  ! Clean up a preconditioner.
!</description>

!<inputoutput>
  ! Structure to be cleaned up.
  type(t_linsolParameters), intent(inout) :: rlinsolParam
!</inputoutput>

!</subroutine>

    ! Release memory
    select case (rlinsolParam%csolverType)
    
    ! -----------------------------------------------------
    ! Multigrid inítialisation
    ! -----------------------------------------------------
    case (NLIN_SOLVER_MULTIGRID)

      deallocate(rlinsolParam%p_rsubnodeMultigrid)
      
    ! -----------------------------------------------------
    ! Richardson inítialisation
    ! -----------------------------------------------------
    case (NLIN_SOLVER_RICHARDSON)
    
      ! Initialise the MG parameters
      deallocate(rlinsolParam%p_rsubnodeRichardson)

    ! -----------------------------------------------------
    ! CG inítialisation
    ! -----------------------------------------------------
    case (NLIN_SOLVER_CG)
    
      ! Initialise the MG parameters
      deallocate(rlinsolParam%p_rsubnodeCG)

    ! -----------------------------------------------------
    ! BiCGSTab inítialisation
    ! -----------------------------------------------------
    case (NLIN_SOLVER_BICGSTAB)
    
      ! Initialise the MG parameters
      deallocate(rlinsolParam%p_rsubnodeBiCGStab)

    end select
   
  end subroutine


  ! ***************************************************************************

!<subroutine>

  subroutine newtonlin_adNewton_setEps (rlinsolParam,depsAbs,depsRel)
  
!<description>
  ! Sets the stopping criterion for the use application of the adaptive
  ! Newton algorithm.
  !
  ! WARNING!!! THIS ROUTINE HAS A SIDE EFFECT!
  ! IT SETS THE STOPPING CRITERION OF ALL SOLVERS IN SPACE TO AN APPROPRIATE
  ! VALUE! IF THE SPACE SOVLERS AE USED SOMEWHERE ELSE, THE STOPPING CRITERION
  ! IS LOST AND THUS, THEY MAY NOT BEHAVE AS EXPECTED!!!
!</description>

!<inputoutput>
  ! Solver structure.
  type(t_linsolParameters), intent(inout) :: rlinsolParam
  
  ! New absolute stopping criterion. Absolute residual.
  ! =0.0: Switch off the check.
  real(DP), intent(in) :: depsAbs

  ! New absolute stopping criterion. Relative residual.
  ! =0.0: Switch off the check.
  real(DP), intent(in) :: depsRel
!</inputoutput>

!</subroutine>

    ! Set the stopping criteria
    rlinsolParam%riter%dtolRel = depsRel
    rlinsolParam%riter%dtolAbs = depsAbs
    
    if (depsAbs .eq. 0.0_DP) then
      rlinsolParam%riter%ctolMode = ITC_STOP_MODE_REL
    else
      rlinsolParam%riter%ctolMode = ITC_STOP_MODE_ABS
    end if
    
    ! Special settings
    select case (rlinsolParam%csolverType)
    
    ! -------------------------------------------
    ! Multigrid
    ! -------------------------------------------
    case (NLIN_SOLVER_MULTIGRID)
    
      ! Is this a one-level solver?
      if (rlinsolParam%p_rkktsystemHierarchy%nlevels .eq. 1) then
        ! The coarse grid solver is used instead of the solver
        ! identified by rlinsolParam. Therefore, impose the stopping
        ! criteria in the coarse grid solver.
        rlinsolParam%p_rsubnodeMultigrid%p_RsubSolvers(1)%riter%dtolRel = depsRel
        rlinsolParam%p_rsubnodeMultigrid%p_RsubSolvers(1)%riter%dtolAbs = depsAbs

        if (depsAbs .eq. 0.0_DP) then
          rlinsolParam%p_rsubnodeMultigrid%p_RsubSolvers(1)%riter%ctolMode = ITC_STOP_MODE_REL
        else
          rlinsolParam%p_rsubnodeMultigrid%p_RsubSolvers(1)%riter%ctolMode = ITC_STOP_MODE_ABS
        end if
        
      else
        
        ! Tweak the coarse grid solver. Gain two digits, but not more than
        ! Specified in the absolute stopping criterion (plus one digit).
        rlinsolParam%p_rsubnodeMultigrid%p_RsubSolvers(1)%riter%ctolMode = ITC_STOP_MODE_ONE_OF
        rlinsolParam%p_rsubnodeMultigrid%p_RsubSolvers(1)%riter%dtolRel = 0.01_DP
        rlinsolParam%p_rsubnodeMultigrid%p_RsubSolvers(1)%riter%dtolAbs = depsAbs * 0.1_DP
        
      end if
    
    end select

!    Do not pass the stopping criterion to the subsolvers in space.
!    This does not work in the current implementation.
!    The residual of the subsolvers in space is measured with a different 
!    measure than the space-time residual: The space-time residual conrtols
!    the control u whereas the space-residual controls y and lambda. The size of
!    the latter ones usually have a completely different value than the norm
!    of the u residual. Thus, taking the norm of the u-residual to control
!    the residual of y and lambda is simply completely wrong!
!    ! Configure the subsolver in space
!    call spaceslh_setEps (&
!        rlinsolParam%p_rsolverHierPrimalLin,depsAbs,depsRel)
!    call spaceslh_setEps (&
!        rlinsolParam%p_rsolverHierDualLin,depsAbs,depsRel)
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonlin_clearStatistics(rstatistics)
  
!<description>
  ! Resets a statistic structure.
!</description>

!<inputoutput>
  ! Structure to be reset.
  type(t_newtonlinSolverStat), intent(inout) :: rstatistics
!</inputoutput>

!</subroutine>

    rstatistics%niterations = 0
    rstatistics%niterationsCoarse = 0

    call stat_clearTimer(rstatistics%rtotalTime)
    call stat_clearTimer(rstatistics%rtimeForwardLin)
    call stat_clearTimer(rstatistics%rtimeBackwardLin)
    call stat_clearTimer(rstatistics%rtimeDefect)
    call stat_clearTimer(rstatistics%rtimeSmoothing)
    call stat_clearTimer(rstatistics%rtimeProlRest)
    
    call spacesl_clearStatistics (rstatistics%rspaceslSolverStat)
    
    call stat_clearTimer(rstatistics%rtimeForwardLinFine)
    call stat_clearTimer(rstatistics%rtimeBackwardLinFine)
    call stat_clearTimer(rstatistics%rtimeDefectFine)
    call stat_clearTimer(rstatistics%rtimeSmoothingFine)
    
    call spacesl_clearStatistics (rstatistics%rspaceslSolverStatFine)
    
    call stat_clearTimer(rstatistics%rtimeForwardLinCoarse)
    call stat_clearTimer(rstatistics%rtimeBackwardLinCoarse)
    call stat_clearTimer(rstatistics%rtimeDefectCoarse)
    call stat_clearTimer(rstatistics%rtimeSolverCoarse)
    
    call spacesl_clearStatistics (rstatistics%rspaceslSolverStatCoarse)
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonlin_sumStatistics(rstatistics1,rstatistics2,csolverType,&
      ilevel,nlevels)
  
!<description>
  ! Sums up the data of rstatistics1 to the data in rstatistics2.
!</description>

!<input>
  ! Source structure
  type(t_newtonlinSolverStat), intent(in) :: rstatistics1
  
  ! Solver type. One of the NLIN_STYPE_xxxx constants.
  integer, intent(in) :: csolverType
  
  ! OPTIONAL: Current level. Can be omitted for csolverType=NLIN_STYPE_LINSOL.
  integer, intent(in), optional :: ilevel
  
  ! Total number of levels in the hierarchy. 
  ! Can be omitted for csolverType=NLIN_STYPE_LINSOL.
  integer, intent(in), optional :: nlevels
!</input>

!<inputoutput>
  ! Destination structure.
  type(t_newtonlinSolverStat), intent(inout) :: rstatistics2
!</inputoutput>

!</subroutine>

    ! What to add where depends on the solver type.
    
    select case (csolverType)
    
    ! ----------------------------------------------------------
    ! Full Linear solver. Sum up all parameters as they are.
    ! ----------------------------------------------------------
    case (NLIN_STYPE_LINSOL,NLIN_STYPE_RESCALC)
    
      rstatistics2%niterations = &
          rstatistics2%niterations + rstatistics1%niterations
      rstatistics2%niterationsCoarse = &
          rstatistics2%niterationsCoarse + rstatistics1%niterationsCoarse

      if (csolverType .eq. NLIN_STYPE_RESCALC) then
        ! Total time goes to the time for the residual calculation
        call stat_addTimers(rstatistics1%rtotalTime,rstatistics2%rtimeDefect)
      else 
        ! Sum up the total time
        call stat_addTimers(rstatistics1%rtotalTime,rstatistics2%rtotalTime)
      end if
      
      call stat_addTimers(rstatistics1%rtimeForwardLin,rstatistics2%rtimeForwardLin)
      call stat_addTimers(rstatistics1%rtimeBackwardLin,rstatistics2%rtimeBackwardLin)
      call stat_addTimers(rstatistics1%rtimeSmoothing,rstatistics2%rtimeSmoothing)
      call stat_addTimers(rstatistics1%rtimeDefect,rstatistics2%rtimeDefect)
      call stat_addTimers(rstatistics1%rtimeProlRest,rstatistics2%rtimeProlRest)
      
      call spacesl_sumStatistics (&
          rstatistics1%rspaceslSolverStat,rstatistics2%rspaceslSolverStat)
      
      call stat_addTimers(rstatistics1%rtimeForwardLinFine,rstatistics2%rtimeForwardLinFine)
      call stat_addTimers(rstatistics1%rtimeBackwardLinFine,rstatistics2%rtimeBackwardLinFine)
      call stat_addTimers(rstatistics1%rtimeDefectFine,rstatistics2%rtimeDefectFine)
      call stat_addTimers(rstatistics1%rtimeSmoothingFine,rstatistics2%rtimeSmoothingFine)
      
      call spacesl_sumStatistics (&
          rstatistics1%rspaceslSolverStatFine,rstatistics2%rspaceslSolverStatFine)
      
      call stat_addTimers(rstatistics1%rtimeForwardLinCoarse,rstatistics2%rtimeForwardLinCoarse)
      call stat_addTimers(rstatistics1%rtimeBackwardLinCoarse,rstatistics2%rtimeBackwardLinCoarse)
      call stat_addTimers(rstatistics1%rtimeDefectCoarse,rstatistics2%rtimeDefectCoarse)
      call stat_addTimers(rstatistics1%rtimeSolverCoarse,rstatistics2%rtimeSolverCoarse)
      
      call spacesl_sumStatistics (&
          rstatistics1%rspaceslSolverStatCoarse,rstatistics2%rspaceslSolverStatCoarse)

    ! ----------------------------------------------------------
    ! Subsolver in a multigrid solver.
    ! On level 1, this is the coarse grid solver.
    ! On level 2..nlevels, this is a smoother.
    ! ----------------------------------------------------------
    case (NLIN_STYPE_MGSUBSOLVER)
    
      if (ilevel .eq. 1) then
        ! ---------------------------------------
        ! General solver statistics
        ! ---------------------------------------
        
        call stat_addTimers(rstatistics1%rtimeForwardLin,rstatistics2%rtimeForwardLin)
        call stat_addTimers(rstatistics1%rtimeBackwardLin,rstatistics2%rtimeBackwardLin)
        call stat_addTimers(rstatistics1%rtimeDefect,rstatistics2%rtimeDefect)
        call stat_addTimers(rstatistics1%rtimeSmoothing,rstatistics2%rtimeSmoothing)
        
        call spacesl_sumStatistics (&
            rstatistics1%rspaceslSolverStat,rstatistics2%rspaceslSolverStat)      

        ! ---------------------------------------
        ! Coarse grid solver statistics
        ! ---------------------------------------
        call stat_addTimers(rstatistics1%rtotalTime,rstatistics2%rtimeSolverCoarse)
        call stat_addTimers(rstatistics1%rtimeForwardLin,rstatistics2%rtimeForwardLinCoarse)
        call stat_addTimers(rstatistics1%rtimeBackwardLin,rstatistics2%rtimeBackwardLinCoarse)
        call stat_addTimers(rstatistics1%rtimeDefect,rstatistics2%rtimeDefectCoarse)
        
        call spacesl_sumStatistics (&
            rstatistics1%rspaceslSolverStat,rstatistics2%rspaceslSolverStatCoarse)      
        
      else 
        ! ---------------------------------------
        ! General smoother statistics
        ! ---------------------------------------
        ! Note that the total time is not summed up here.
        ! It is always calculated outside of the iteration.
        
        call stat_addTimers(rstatistics1%rtotalTime,rstatistics2%rtimeSmoothing)
        call stat_addTimers(rstatistics1%rtimeForwardLin,rstatistics2%rtimeForwardLin)
        call stat_addTimers(rstatistics1%rtimeBackwardLin,rstatistics2%rtimeBackwardLin)
        call stat_addTimers(rstatistics1%rtimeDefect,rstatistics2%rtimeDefect)
        
        call spacesl_sumStatistics (&
            rstatistics1%rspaceslSolverStat,rstatistics2%rspaceslSolverStat)      
      
        if ((ilevel .eq. nlevels)) then
          ! ---------------------------------------
          ! Fine grid smoother statistics
          ! ---------------------------------------
          call stat_addTimers(rstatistics1%rtotalTime,rstatistics2%rtimeSmoothingFine)
          call stat_addTimers(rstatistics1%rtimeForwardLin,rstatistics2%rtimeForwardLinFine)
          call stat_addTimers(rstatistics1%rtimeBackwardLin,rstatistics2%rtimeBackwardLinFine)
          call stat_addTimers(rstatistics1%rtimeDefect,rstatistics2%rtimeDefectFine)
          
          call spacesl_sumStatistics (&
              rstatistics1%rspaceslSolverStat,rstatistics2%rspaceslSolverStatFine)      
        end if
        
      end if
      
    end select

  end subroutine

end module
