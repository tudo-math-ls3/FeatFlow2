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
  use spacetimeinterlevelprojection
  
  use structuresnewton
  
  implicit none
  
  private

!<constants>

!<constantblock description = "Solver identifier">

  ! Richardson iteration
  integer, parameter, public :: NLIN_SOLVER_RICHARDSON = 0
  
  ! Multigrid iteration
  integer, parameter, public :: NLIN_SOLVER_MULTIGRID = 1

  ! CG iteration
  integer, parameter, public :: NLIN_SOLVER_CG = 2

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
    ! =0: Damped Richardson iteration
    integer :: csmoother = 0

    ! Specifies the coarse grid solver.
    ! =0: Damped Richardson iteration
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
  end type

  !</typeblock>
  
  public :: t_linsolCG
  
  !<typeblock>
  
  ! Richardson parameters and structures
  type t_linsolRichardson
  
    ! Temporary memory on all levels for calculations
    type(t_controlSpace) :: rtempVector
  
  end type

  !</typeblock>
  
  public :: t_linsolRichardson
  
  !<typeblock>
  
  ! Contains the parameters of the Newton iteration.
  type t_linsolParameters
  
    ! <!-- --------------------------------------- -->
    ! <!-- GENRERAL PARAMETERS AND SOLVER SETTINGS -->
    ! <!-- --------------------------------------- -->

    ! Specifies the type of the solver.
    ! =0: Damped Richardson iteration with damping parameter omega.
    ! =1: Space-time multigrid method applied to the underlying 
    !     space-time hierarchy
    integer :: csolverType = 0

    ! General newton parameters
    type(t_newtonPrecParameters) :: rprecParameters

    ! Parameters of the OptFlow solver
    type(t_settings_optflow), pointer :: p_rsettingsSolver => null()
  
!    ! <!-- ----------------------------- -->
!    ! <!-- STATISTICS                    -->
!    ! <!-- ----------------------------- -->
!    
!    ! Total computation time
!    type(t_timer) :: rtotalTime
!    
!    ! Total time for smoothing
!    type(t_timer) :: rtimeSmooth
!
!    ! Total time for smoothing on the finest mesh
!    type(t_timer) :: rtimeSmoothFinest
!
!    ! Total time for coarse grid solving
!    type(t_timer) :: rtimeCoarseGrid
!
!    ! Total time for prolongation/restriction
!    type(t_timer) :: rtimeProlRest

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
    
    ! Parameters for the multigrid solver
    type(t_linsolMultigrid), pointer :: p_rsubnodeMultigrid => null()

    ! Parameters for the Richardson solver
    type(t_linsolRichardson), pointer :: p_rsubnodeRichardson => null()

    ! Parameters for the CG solver
    type(t_linsolCG), pointer :: p_rsubnodeCG => null()
    
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

    if (rlinsolParam%rprecParameters%ioutputLevel .ge. 2) then
      call output_line ("Linear space-time Residual: Solving the primal equation")
    end if

    ! Solve the primal equation, update the primal solution.
    output_iautoOutputIndent = output_iautoOutputIndent + 2

    call kkt_solvePrimalDirDeriv (rkktsystemDirDeriv,&
        rlinsolParam%p_rsolverHierPrimalLin,rlinsolParam%cspatialInitCondPolicy,&
        rlocalStat)

    output_iautoOutputIndent = output_iautoOutputIndent - 2
    
    call spacesl_sumStatistics(rlocalStat,rstatistics%rspaceslSolverStat)
    
    if (rlinsolParam%rprecParameters%ioutputLevel .ge. 3) then
      call output_line ("Linear space-time Residual: Time for solving: "//&
          trim(sys_sdL(rlocalStat%rtotalTime%delapsedReal,10)))
    end if

    ! Add time to the time of the forward equation
    call stat_addTimers (rlocalStat%rtotalTime,rstatistics%rtimeForwardLin)

    ! --------------------
    ! Backward equation
    ! --------------------

    if (rlinsolParam%rprecParameters%ioutputLevel .ge. 2) then
      call output_line ("Linear space-time Residual: Solving the dual equation")
    end if

    ! Solve the dual equation, update the dual solution.
    output_iautoOutputIndent = output_iautoOutputIndent + 2

    call kkt_solveDualDirDeriv (rkktsystemDirDeriv,&
        rlinsolParam%p_rsolverHierDualLin,rlinsolParam%cspatialInitCondPolicy,&
        rlocalStat)

    output_iautoOutputIndent = output_iautoOutputIndent - 2
    
    call spacesl_sumStatistics(rlocalStat,rstatistics%rspaceslSolverStat)
    
    if (rlinsolParam%rprecParameters%ioutputLevel .ge. 3) then
      call output_line ("Linear space-time Residual: Time for solving: "//&
          trim(sys_sdL(rlocalStat%rtotalTime%delapsedReal,10)))
    end if

    ! Add time to the time of the backward equation
    call stat_addTimers (rlocalStat%rtotalTime,rstatistics%rtimeBackwardLin)

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

    ! -------------------------------------------------------------
    ! Step 1: Solve the primal and dual system.
    ! -------------------------------------------------------------

    ! --------------------
    ! Forward equation
    ! --------------------

    if (rlinsolParam%rprecParameters%ioutputLevel .ge. 2) then
      call output_line ("Linear space-time Residual: Solving the primal equation")
    end if

    ! Solve the primal equation, update the primal solution.
    output_iautoOutputIndent = output_iautoOutputIndent + 2

    call kkt_solvePrimalDirDeriv (rkktsystemDirDeriv,&
        rlinsolParam%p_rsolverHierPrimalLin,rlinsolParam%cspatialInitCondPolicy,&
        rstatistics%rspaceslSolverStat)
    
    output_iautoOutputIndent = output_iautoOutputIndent - 2

    if (rlinsolParam%rprecParameters%ioutputLevel .ge. 3) then
      call output_line ("Linear space-time Residual: Time for solving: "//&
          trim(sys_sdL(rlocalStat%rtotalTime%delapsedReal,10)))
    end if

    ! Add time to the time of the forward equation
    call stat_addTimers (rlocalStat%rtotalTime,rstatistics%rtimeForwardLin)

    ! --------------------
    ! Backward equation
    ! --------------------

    if (rlinsolParam%rprecParameters%ioutputLevel .ge. 2) then
      call output_line ("Linear space-time Residual: Solving the dual equation")
    end if

    ! Solve the dual equation, update the dual solution.
    output_iautoOutputIndent = output_iautoOutputIndent + 2

    call kkt_solveDualDirDeriv (rkktsystemDirDeriv,&
        rlinsolParam%p_rsolverHierDualLin,rlinsolParam%cspatialInitCondPolicy,&
        rstatistics%rspaceslSolverStat)

    output_iautoOutputIndent = output_iautoOutputIndent - 2

    if (rlinsolParam%rprecParameters%ioutputLevel .ge. 3) then
      call output_line ("Linear space-time Residual: Time for solving: "//&
          trim(sys_sdL(rlocalStat%rtotalTime%delapsedReal,10)))
    end if

    ! Add time to the time of the backward equation
    call stat_addTimers (rlocalStat%rtotalTime,rstatistics%rtimeBackwardLin)

    ! -------------------------------------------------------------
    ! Step 2: Calculate the residual
    ! -------------------------------------------------------------

    ! Take the solution of the linearised primal/dual system and
    ! calculate the residual in the control space.
    call kkt_applyControlDirDeriv (rkktsystemDirDeriv,rrhs)

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
    rsmootherParams%rprecParameters%nminiterations = nsm
    rsmootherParams%rprecParameters%nmaxiterations = nsm
    
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

    ! Measure iteration time
    call stat_startTimer (rstatistics%rtotalTime)

    ! Get the temp vector
    p_rtempVector => rlinsolParam%p_rsubnodeRichardson%rtempVector

    ! Richardson is a do loop...
    rlinsolParam%rprecParameters%niterations = 0

    do while (.true.)
    
      ! -------------------------------------------------------------
      ! Get the current residual / search direction
      ! -------------------------------------------------------------

      ! Compute the residual and its norm.
      output_iautoOutputIndent = output_iautoOutputIndent + 2
      call newtonlin_getResidual (rlinsolParam,rkktsystemDirDeriv,rrhs,&
          p_rtempVector,rlocalStat,rlinsolParam%rprecParameters%dresFinal,&
          rlinsolParam%rprecParameters%iresnorm)
      output_iautoOutputIndent = output_iautoOutputIndent - 2
      
      call newtonlin_sumStatistics(rlocalStat,rstatistics,NLIN_STYPE_RESCALC)
      
      if (rlinsolParam%rprecParameters%niterations .eq. 0) then
        ! Remember the initial residual
        rlinsolParam%rprecParameters%dresInit = rlinsolParam%rprecParameters%dresFinal
      end if

      if (rlinsolParam%rprecParameters%ioutputLevel .ge. 2) then
        call output_line ("Space-time Richardson: Iteration "// &
            trim(sys_siL(rlinsolParam%rprecParameters%niterations,10))// &
            ", ||res(u)|| = "// &
            trim(sys_sdEL(rlinsolParam%rprecParameters%dresFinal,10)))
      end if

      ! -------------------------------------------------------------
      ! Check for convergence
      ! -------------------------------------------------------------
      if (newtonlin_checkConvergence(rlinsolParam%rprecParameters)) exit
      
      ! -------------------------------------------------------------
      ! Check for divergence
      ! -------------------------------------------------------------
      if (newtonlin_checkDivergence(rlinsolParam%rprecParameters)) exit
      
      ! -------------------------------------------------------------
      ! Check other stopping criteria
      ! -------------------------------------------------------------
      if (newtonlin_checkIterationStop(rlinsolParam%rprecParameters)) exit
      
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
      ! Next iteration
      rlinsolParam%rprecParameters%niterations = &
          rlinsolParam%rprecParameters%niterations + 1
    
    end do

    call stat_stopTimer (rstatistics%rtotalTime)

    ! Statistics
    if (rlinsolParam%rprecParameters%dresInit .gt. &
        rlinsolParam%rprecParameters%drhsZero) then
      rlinsolParam%rprecParameters%dconvergenceRate = &
                  (rlinsolParam%rprecParameters%dresFinal / &
                   rlinsolParam%rprecParameters%dresInit) ** &
                  (1.0_DP/real(rlinsolParam%rprecParameters%niterations,DP))
    else
      rlinsolParam%rprecParameters%dconvergenceRate = 0.0_DP
    end if

    rstatistics%niterations = rstatistics%niterations + rlinsolParam%rprecParameters%niterations

    if (rlinsolParam%rprecParameters%ioutputLevel .ge. 2) then
      call output_lbrk()
      call output_line ("Space-time Richardson: Statistics.")
      call output_lbrk()
      call output_line ("#Iterations             : "//&
            trim(sys_siL(rlinsolParam%rprecParameters%niterations,10)) )
      call output_line ("!!INITIAL RES!!         : "//&
            trim(sys_sdEL(rlinsolParam%rprecParameters%dresInit,15)) )
      call output_line ("!!RES!!                 : "//&
            trim(sys_sdEL(rlinsolParam%rprecParameters%dresFinal,15)) )
      if (rlinsolParam%rprecParameters%dresInit .gt. &
          rlinsolParam%rprecParameters%drhsZero) then
        call output_line ("!!RES!!/!!INITIAL RES!! : "//&
          trim(sys_sdEL(rlinsolParam%rprecParameters%dresFinal / &
                        rlinsolParam%rprecParameters%dresInit,15)) )
      else
        call output_line ("!!RES!!/!!INITIAL RES!! : "//&
              trim(sys_sdEL(0.0_DP,15)) )
      end if
      call output_lbrk ()
      call output_line ('Rate of convergence     : '//&
            trim(sys_sdEL(rlinsolParam%rprecParameters%dconvergenceRate,15)) )

    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonlin_cg (rlinsolParam,rkktsystemDirDeriv,rrhs,rstatistics)
  
!<description>
  ! Applies a Richardson iteration in the control space
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

    real(DP) :: dalpha,dbeta,dgamma,dgammaOld
    type(t_controlSpace), pointer :: p_rr,p_rAp
    type(t_kktsystemDirDeriv), pointer :: p_rp
    type(t_newtonlinSolverStat) :: rlocalStat

    ! Measure total time
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

    ! Create the initial defect in rd
    output_iautoOutputIndent = output_iautoOutputIndent + 2
    call newtonlin_getResidual (rlinsolParam,rkktsystemDirDeriv,rrhs,p_rr,rlocalStat)
    output_iautoOutputIndent = output_iautoOutputIndent - 2
        
    call newtonlin_sumStatistics(rlocalStat,rstatistics,NLIN_STYPE_RESCALC)
        
    call kktsp_controlCopy (p_rr,p_rp%p_rcontrolLin)

    ! Scalar product of rp.
    dgamma = kktsp_scalarProductControl(p_rr,p_rr)

    ! Richardson is a do loop...
    rlinsolParam%rprecParameters%niterations = 0

    do while (.true.)
    
      ! -------------------------------------------------------------
      ! Norm of the residual
      ! -------------------------------------------------------------

      rlinsolParam%rprecParameters%dresFinal = &
          kktsp_getNormControl (p_rr,rlinsolParam%rprecParameters%iresnorm)

      if (rlinsolParam%rprecParameters%niterations .eq. 0) then
        ! Remember the initial residual
        rlinsolParam%rprecParameters%dresInit = rlinsolParam%rprecParameters%dresFinal
      end if

      if (rlinsolParam%rprecParameters%ioutputLevel .ge. 2) then
        call output_line ("Space-time CG: Iteration "// &
            trim(sys_siL(rlinsolParam%rprecParameters%niterations,10))// &
            ", ||res(u)|| = "// &
            trim(sys_sdEL(rlinsolParam%rprecParameters%dresFinal,10)))
      end if

      ! -------------------------------------------------------------
      ! Check for convergence
      ! -------------------------------------------------------------
      if (newtonlin_checkConvergence(rlinsolParam%rprecParameters)) exit
      
      ! -------------------------------------------------------------
      ! Check for divergence
      ! -------------------------------------------------------------
      if (newtonlin_checkDivergence(rlinsolParam%rprecParameters)) exit
      
      ! -------------------------------------------------------------
      ! Check other stopping criteria
      ! -------------------------------------------------------------
      if (newtonlin_checkIterationStop(rlinsolParam%rprecParameters)) exit
      
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
      
      call newtonlin_sumStatistics(rlocalStat,rstatistics,NLIN_STYPE_LINSOL)

      ! Calculate the parameter ALPHA
      dalpha = dgamma / kktsp_scalarProductControl(p_rp%p_rcontrolLin,p_rAp)

      ! Add the residual (damped) to the current control.
      ! This updates the current control.
      call kktsp_controlLinearComb (p_rp%p_rcontrolLin,dalpha,&
          rkktsystemDirDeriv%p_rcontrolLin,1.0_DP)

      ! Update the residual
      call kktsp_controlLinearComb (p_rAp,-dalpha,p_rr,1.0_DP)
          
      ! Calculate beta
      dgammaOld = dgamma
      dgamma = kktsp_scalarProductControl(p_rr,p_rr)
      dbeta = dgamma / dgammaOld
      
      ! Update rp.
      call kktsp_controlLinearComb (p_rr,1.0_DP,p_rp%p_rcontrolLin,dbeta)
    
      ! -------------------------------------------------------------
      ! Proceed with the iteration
      ! -------------------------------------------------------------
      ! Next iteration
      rlinsolParam%rprecParameters%niterations = &
          rlinsolParam%rprecParameters%niterations + 1
    
    end do
    
    call stat_stopTimer (rstatistics%rtotalTime)

    ! Statistics
    if (rlinsolParam%rprecParameters%dresInit .gt. &
        rlinsolParam%rprecParameters%drhsZero) then
      rlinsolParam%rprecParameters%dconvergenceRate = &
                  (rlinsolParam%rprecParameters%dresFinal / &
                   rlinsolParam%rprecParameters%dresInit) ** &
                  (1.0_DP/real(rlinsolParam%rprecParameters%niterations,DP))
    else
      rlinsolParam%rprecParameters%dconvergenceRate = 0.0_DP
    end if

    rstatistics%niterations = rstatistics%niterations + rlinsolParam%rprecParameters%niterations
    
    if (rlinsolParam%rprecParameters%ioutputLevel .ge. 2) then
      call output_lbrk()
      call output_line ("Space-time CG: Statistics.")
      call output_lbrk()
      call output_line ("#Iterations             : "//&
            trim(sys_siL(rlinsolParam%rprecParameters%niterations,10)) )
      call output_line ("!!INITIAL RES!!         : "//&
            trim(sys_sdEL(rlinsolParam%rprecParameters%dresInit,15)) )
      call output_line ("!!RES!!                 : "//&
            trim(sys_sdEL(rlinsolParam%rprecParameters%dresFinal,15)) )
      if (rlinsolParam%rprecParameters%dresInit .gt. &
          rlinsolParam%rprecParameters%drhsZero) then
        call output_line ("!!RES!!/!!INITIAL RES!! : "//&
          trim(sys_sdEL(rlinsolParam%rprecParameters%dresFinal / &
                        rlinsolParam%rprecParameters%dresInit,15)) )
      else
        call output_line ("!!RES!!/!!INITIAL RES!! : "//&
              trim(sys_sdEL(0.0_DP,15)) )
      end if
      call output_lbrk ()
      call output_line ('Rate of convergence     : '//&
            trim(sys_sdEL(rlinsolParam%rprecParameters%dconvergenceRate,15)) )

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

  subroutine newtonlin_initStrucMultigrid (rsolver,rlinsolMultigrid,rkktsystemHierarchy)
  
!<description>
  ! Initialises a multigrid substructure.
!</description>

!<input>
  ! Solver structure of the main solver.
  type(t_linsolParameters), intent(in) :: rsolver

  ! Defines the basic hierarchy of the solutions of the KKT system.
  ! This can be a 'template' structure, i.e., memory for the solutions
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
    case (NLIN_SOLVER_RICHARDSON)
    
      ! Partial initialisation of the corresponding linear solver structure
      rlinsolMultigrid%p_rsubSolvers(1)%csolverType = NLIN_SOLVER_RICHARDSON
      call newtonlin_initBasicParams (&
          rlinsolMultigrid%p_rsubSolvers(1)%rprecParameters,&
          rlinsolMultigrid%ssectionCoarseGridSolver,rlinsolMultigrid%p_rparList)

      allocate(rlinsolMultigrid%p_rsubSolvers(1)%p_rsubnodeRichardson)

    ! -------------------------------
    ! CG solver
    ! -------------------------------
    case (NLIN_SOLVER_CG)
    
      ! Partial initialisation of the corresponding linear solver structure
      rlinsolMultigrid%p_rsubSolvers(1)%csolverType = NLIN_SOLVER_CG
      call newtonlin_initBasicParams (&
          rlinsolMultigrid%p_rsubSolvers(1)%rprecParameters,&
          rlinsolMultigrid%ssectionCoarseGridSolver,rlinsolMultigrid%p_rparList)
      allocate(rlinsolMultigrid%p_rsubSolvers(1)%p_rsubnodeCG)

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
      case (NLIN_SOLVER_RICHARDSON)

        ! Partial initialisation of the corresponding linear solver structure
        rlinsolMultigrid%p_rsubSolvers(ilevel)%csolverType = NLIN_SOLVER_RICHARDSON
        call newtonlin_initBasicParams (&
            rlinsolMultigrid%p_rsubSolvers(ilevel)%rprecParameters,&
            rlinsolMultigrid%ssectionSmoother,rlinsolMultigrid%p_rparList)
        
        allocate(rlinsolMultigrid%p_rsubSolvers(ilevel)%p_rsubnodeRichardson)

      ! -------------------------------
      ! CG smoother
      ! -------------------------------
      case (NLIN_SOLVER_CG)

        ! Partial initialisation of the corresponding linear solver structure
        rlinsolMultigrid%p_rsubSolvers(ilevel)%csolverType = NLIN_SOLVER_CG
        call newtonlin_initBasicParams (&
            rlinsolMultigrid%p_rsubSolvers(ilevel)%rprecParameters,&
            rlinsolMultigrid%ssectionSmoother,rlinsolMultigrid%p_rparList)
        
        allocate(rlinsolMultigrid%p_rsubSolvers(ilevel)%p_rsubnodeCG)

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

    end do

    ! Create temp vectors for the iteration based on the given KKT system hierarchy
    call kkth_initHierarchyDirDeriv (rlinsolMultigrid%rsolutionHier,rsolutionHierarchy)
    call kkth_initHierarchyDirDeriv (rlinsolMultigrid%rdefectHier,rsolutionHierarchy)
    call kkth_initHierarchyDirDeriv (rlinsolMultigrid%rrhsHier,rsolutionHierarchy)
   
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
    if (rlinsolParam%rprecParameters%dresInit .gt. &
        rlinsolParam%rprecParameters%drhsZero) then
      rlinsolParam%rprecParameters%dconvergenceRate = &
                  (rlinsolParam%rprecParameters%dresFinal / &
                   rlinsolParam%rprecParameters%dresInit) ** &
                  (1.0_DP/real(rlinsolParam%rprecParameters%niterations,DP))
    else
      rlinsolParam%rprecParameters%dconvergenceRate = 0.0_DP
    end if

    if (rlinsolParam%rprecParameters%ioutputLevel .ge. 2) then
      call output_lbrk()
      call output_line ("Space-time Multigrid: Statistics.")
      call output_lbrk()
      call output_line ("#Iterations             : "//&
            trim(sys_siL(rlinsolParam%rprecParameters%niterations,10)) )
      call output_line ("!!INITIAL RES!!         : "//&
            trim(sys_sdEL(rlinsolParam%rprecParameters%dresInit,15)) )
      call output_line ("!!RES!!                 : "//&
            trim(sys_sdEL(rlinsolParam%rprecParameters%dresFinal,15)) )
      if (rlinsolParam%rprecParameters%dresInit .gt. &
          rlinsolParam%rprecParameters%drhsZero) then
        call output_line ("!!RES!!/!!INITIAL RES!! : "//&
          trim(sys_sdEL(rlinsolParam%rprecParameters%dresFinal / &
                        rlinsolParam%rprecParameters%dresInit,15)) )
      else
        call output_line ("!!RES!!/!!INITIAL RES!! : "//&
              trim(sys_sdEL(0.0_DP,15)) )
      end if
      call output_lbrk ()
      call output_line ('Rate of convergence     : '//&
            trim(sys_sdEL(rlinsolParam%rprecParameters%dconvergenceRate,15)) )

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
    integer :: nminIterations,nmaxIterations,ite
    real(DP) :: dresInit,dresFinal
    type(t_kktsystemDirDeriv), pointer :: p_rkktSysRhsCoarse, p_rkktSysSolCoarse, p_rkktSysSolFine
    type(t_kktsystemDirDeriv), pointer :: p_rkktSysSolution, p_rkktSysDefect
    type(t_newtonlinSolverStat) :: rlocalStat
    
    ! Get multigrid parameters.
    p_rsubnodeMultigrid => rlinsolParam%p_rsubnodeMultigrid

    ! On level 1 apply the coarse grid solver
    if (ilevel .eq. 1) then
      
      if (nlevels .eq. 1) then
        if (rlinsolParam%rprecParameters%ioutputLevel .ge. 1) then
          call output_line (&
              "Space-time Multigrid: Only one level. Switching back to 1-level solver.")
        end if
      end if

      output_iautoOutputIndent = output_iautoOutputIndent + 2
      call newtonlin_precond (p_rsubnodeMultigrid%p_rsubSolvers(ilevel),&
          rkktsysDirDerivHierarchy,rrhs,rstatistics,ilevel)
      output_iautoOutputIndent = output_iautoOutputIndent - 2
    
      rlinsolParam%rprecParameters%dresInit = &
          p_rsubnodeMultigrid%p_rsubSolvers(ilevel)%rprecParameters%dresInit
      rlinsolParam%rprecParameters%dresFinal = &
          p_rsubnodeMultigrid%p_rsubSolvers(ilevel)%rprecParameters%dresFinal
      rlinsolParam%rprecParameters%niterations = &
          p_rsubnodeMultigrid%p_rsubSolvers(ilevel)%rprecParameters%niterations
    
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
      nminIterations = rlinsolParam%rprecParameters%nminIterations
      nmaxIterations = rlinsolParam%rprecParameters%nmaxIterations
    else
      ! On all the other levels, apply as many iterations as configured in the cycle.
      nminIterations = p_rsubnodeMultigrid%p_IcycleCount(ilevel)
      nmaxIterations = p_rsubnodeMultigrid%p_IcycleCount(ilevel)
    end if
    
    ! Apply the multigrid iteration

    call kkth_getKKTsystemDirDeriv (rkktsysDirDerivHierarchy,ilevel,p_rkktSysSolution)
    call kkth_getKKTsystemDirDeriv (p_rsubnodeMultigrid%rrhsHier,ilevel-1,p_rkktSysRhsCoarse)
    call kkth_getKKTsystemDirDeriv (p_rsubnodeMultigrid%rsolutionHier,ilevel-1,p_rkktSysSolCoarse)
    call kkth_getKKTsystemDirDeriv (p_rsubnodeMultigrid%rsolutionHier,ilevel,p_rkktSysSolFine)
    call kkth_getKKTsystemDirDeriv (p_rsubnodeMultigrid%rdefectHier,ilevel,p_rkktSysDefect)
    
    ite = 0
    call kkt_clearDirDeriv (p_rkktSysSolution)

    do
    
      if (ilevel .lt. nlevels) then
        ! On level < nlevels, do not calculate the next residual but
        ! exit immediately if the maximum number of iterations
        ! (= cycle count) is reached.
        if (ite .ge. nmaxIterations) exit
      end if
      
      ! Get the residual
      output_iautoOutputIndent = output_iautoOutputIndent + 2
      
      call newtonlin_getResidual (rlinsolParam,&
          p_rkktSysSolution,rrhs,&
          p_rkktSysDefect%p_rcontrolLin,rlocalStat,&
          dresFinal,rlinsolParam%rprecParameters%iresnorm)
      
      output_iautoOutputIndent = output_iautoOutputIndent - 2
      
      call newtonlin_sumStatistics(rlocalStat,rstatistics,NLIN_STYPE_RESCALC)

      ! Remember the initial residual
      if (ite .eq. 0) dresInit = dresFinal

      ! Some checks and output on the maximum level
      if (ilevel .eq. nlevels) then
      
        if (rlinsolParam%rprecParameters%ioutputLevel .ge. 2) then
          call output_line ("Space-time MG: Iteration "// &
              trim(sys_siL(ite,10))// &
              ", ||res(u)|| = "// &
              trim(sys_sdEL(dresFinal,10)))
        end if
        
        ! Save statistics on the maximum level
        rlinsolParam%rprecParameters%dresInit = dresInit
        rlinsolParam%rprecParameters%dresFinal = dresFinal
        rlinsolParam%rprecParameters%niterations = ite

        ! -------------------------------------------------------------
        ! Check for convergence
        ! -------------------------------------------------------------
        if (newtonlin_checkConvergence(rlinsolParam%rprecParameters)) exit
        
        ! -------------------------------------------------------------
        ! Check for divergence
        ! -------------------------------------------------------------
        if (newtonlin_checkDivergence(rlinsolParam%rprecParameters)) exit
        
        ! -------------------------------------------------------------
        ! Check other stopping criteria
        ! -------------------------------------------------------------
        if (newtonlin_checkIterationStop(rlinsolParam%rprecParameters)) exit
      
      end if
      
      ! Apply presmoothing
      if (p_rsubnodeMultigrid%nsmpre .ne. 0) then
        if (rlinsolParam%rprecParameters%ioutputLevel .ge. 3) then
          call output_line (&
              "Space-time Multigrid: Invoking pre-smoother.")
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
            dresFinal,rlinsolParam%rprecParameters%iresnorm)
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

      ! And the actual correction...
      call kktsp_controlLinearComb (&
          p_rkktSysDefect%p_rcontrolLin,&
          p_rsubnodeMultigrid%dcoarseGridCorrectionWeight,&
          p_rkktSysSolution%p_rcontrolLin,&
          1.0_DP)

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
      ite = ite + 1

      if (ilevel .eq. nlevels) then
        ! Re-initialise the cycles
        call MG_initCycles (&
            p_rsubnodeMultigrid%ccycle,ilevel,nlevels,p_rsubnodeMultigrid%p_IcycleCount)
      end if

    end do
    
    rstatistics%niterations = rstatistics%niterations + ite

    ! Initialise the cycle counter for the next run.
    call MG_initCycles (&
        p_rsubnodeMultigrid%ccycle,ilevel,nlevels,p_rsubnodeMultigrid%p_IcycleCount)
        
    call stat_stopTimer(rstatistics%rtotalTime)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  recursive subroutine newtonlin_precond (&
      rlinsolParam,rkktsysDirDerivHierarchy,rnewtonDir,rstatistics,ilevel)
  
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
  type(t_kktsystemDirDerivHierarchy), intent(inout) :: rkktsysDirDerivHierarchy

  ! This structure contains the search direction to be
  ! preconditioned. Will be replaced by the precoditioned search direction.
  type(t_controlSpace), intent(inout) :: rnewtonDir
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
    ! Multigrid iteration
    ! --------------------------------------
    case (NLIN_SOLVER_MULTIGRID)
      ! Invoke multigrid.
      call newtonlin_multigrid (rlinsolParam,rkktsysDirDerivHierarchy,rnewtonDir,rstatistics,ilv)
      
    end select
        
    ! Overwrite the rnewtonDir with the update.
    call kktsp_controlLinearComb (&
        p_rkktsysDirDeriv%p_rcontrolLin,1.0_DP,rnewtonDir,0.0_DP)
    
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
    call newtonlin_initBasicParams (rlinsolParam%rprecParameters,ssection,rparlist)

    ! Solver type
    call parlst_getvalue_int (rparlist, ssection, "csolverType", &
        rlinsolParam%csolverType, rlinsolParam%csolverType)

    ! Generation of the initial condition
    call parlst_getvalue_int (rparlist, ssection, "cspatialInitCondPolicy", &
        rlinsolParam%cspatialInitCondPolicy, rlinsolParam%cspatialInitCondPolicy)
        
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
  ! This can be a 'template' structure, i.e., memory for the solutions
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

    integer :: ilev
   
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
  ! IS LOST AND THUS, THEY MAY BEHAVE NOT AS EXPECTED!!!
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
    rlinsolParam%rprecParameters%depsRel = depsRel
    rlinsolParam%rprecParameters%depsAbs = depsAbs

    ! Configure the subsolver in space
    call spaceslh_setEps (&
        rlinsolParam%p_rsolverHierPrimalLin,depsAbs,depsRel)
    call spaceslh_setEps (&
        rlinsolParam%p_rsolverHierDualLin,depsAbs,depsRel)
    
  contains
  
    subroutine setEpsSpaceSolver (rlinsol,depsAbs)
    
    type(t_linsolSpace), intent(inout) :: rlinsol
    real(DP) :: depsAbs
    
      ! Solver type?
      select case (rlinsol%isolverType)
      
      ! ------------------------------
      ! UMFPACK. Nothing to do.
      ! ------------------------------
      case (LSS_LINSOL_UMFPACK)
      
      ! ------------------------------
      ! Multigrid
      ! ------------------------------
      case (LSS_LINSOL_MG)
      
        ! Stopping criterion of the MG solver in space
        rlinsol%p_rsolverNode%depsRel = depsRel
        rlinsol%p_rsolverNode%depsAbs = depsAbs
        
        ! Coarse grid solver does not have to be changed.
        ! If it has the wrong stopping criterion, the iteration 
        ! just takes longer...
      
      case default
        call output_line ("Unknown solver in space.", &
            OU_CLASS_ERROR,OU_MODE_STD,"newtonlin_adNewton_setEps")
        call sys_halt()
      
      end select
    
    end subroutine

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
