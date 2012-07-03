!##############################################################################
!# ****************************************************************************
!# <name> newtoniteration </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module realises a (semismooth) Newton iteration in the control 
!# space of the optimal control problem.
!# </purpose>
!##############################################################################

module newtoniteration

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
  use structuresnewton
  use assemblytemplates
  
  use spacematvecassembly
  use spacelinearsolver
  use spacesolver
  
  use spacetimehierarchy
  use spacetimeinterlevelprojection
  use kktsystemspaces
  use kktsystem
  use kktsystemhierarchy
  
  use newtoniterationlinear

  implicit none
  
  private

!<types>

!<typeblock>

  ! Linear solver statistics
  type t_newtonitSolverStat
  
    ! Number of iterations necessary for the solver
    integer :: niterations = 0
    
    ! Total time necessary for the solver
    type(t_timer) :: rtotalTime

    ! Time necessary for solving the nonlinear forward problem
    type(t_timer) :: rtimeForward

    ! Time necessary for solving the backward problem
    type(t_timer) :: rtimeBackward
    
    ! Time necessary for the creation of nonlinear defects
    type(t_timer) :: rtimeDefect

    ! Time necessary for the prolongation of the solution to all levels
    type(t_timer) :: rtimeProlRest
    
    ! Time necessary for postprocessing
    type(t_timer) :: rtimePostprocessing

    ! Statistics of the subsolver in space used for calculating the residual.
    type(t_spaceslSolverStat) :: rspaceslSolverStat

    ! Statistics of the linear subsolver
    type(t_newtonlinSolverStat) :: rnewtonlinSolverStat
    
  end type

!</typeblock>

  public :: t_newtonitSolverStat

!<typeblock>
  
  ! Parameters for the space-time Newton algorithm.
  type t_spacetimeNewton

    ! <!-- --------------------------------------- -->
    ! <!-- GENRERAL PARAMETERS AND SOLVER SETTINGS -->
    ! <!-- --------------------------------------- -->

    ! General newton parameters
    type(t_newtonParameters) :: rnewtonParams

    ! Parameters of the OptFlow solver
    type(t_settings_optflow), pointer :: p_rsettingsSolver => null()

!    ! <!-- ----------------------------- -->
!    ! <!-- STATISTICS                    -->
!    ! <!-- ----------------------------- -->
!    
!    ! Total computation time
!    type(t_timer) :: rtotalTime
!
!    ! Total time for nonlinear defects
!    type(t_timer) :: rtimeNonlinearDefects
!    
!    ! Total time for solviong linear subproblems
!    type(t_timer) :: rtimePreconditioning
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

    ! Defines a policy how to generate the initial condition of a timestep.
    ! =0: Always take zero
    ! =1: Propagate the solution of the previous/next timestep to the
    !     current one. (Default)
    ! =2: Take the solution of the last space-time iteration
    integer :: cspatialInitCondPolicy = SPINITCOND_PREVTIMESTEP

    ! Parameters for the linear space-time subsolver.
    type(t_linsolParameters) :: rlinsolParam
    
    ! Hierarchy of linear solvers in space for all levels.
    ! Nonlinear primal equation.
    type(t_linsolHierarchySpace), pointer :: p_rlinsolHierPrimal => null()

    ! Hierarchy of linear solvers in space for all levels.
    ! Dual equation.
    type(t_linsolHierarchySpace), pointer :: p_rlinsolHierDual => null()

    ! Hierarchy of linear solvers in space for all levels.
    ! Linearised primal equation.
    type(t_linsolHierarchySpace), pointer :: p_rlinsolHierPrimalLin => null()

    ! Hierarchy of linear solvers in space for all levels.
    ! Linearised dual equation.
    type(t_linsolHierarchySpace), pointer :: p_rlinsolHierDualLin => null()

    ! Hierarchy of solvers in space for all levels.
    ! Nonlinear primal equation.
    type(t_spaceSolverHierarchy), pointer :: p_rsolverHierPrimal => null()

    ! Hierarchy of solvers in space for all levels.
    ! Dual equation.
    type(t_spaceSolverHierarchy), pointer :: p_rsolverHierDual => null()

    ! Hierarchy of solvers in space for all levels.
    ! Linearised primal equation.
    type(t_spaceSolverHierarchy), pointer :: p_rsolverHierPrimalLin => null()

    ! Hierarchy of solvers in space for all levels.
    ! Linearised dual equation.
    type(t_spaceSolverHierarchy), pointer :: p_rsolverHierDualLin => null()
    
    ! <!-- -------------- -->
    ! <!-- TEMPORARY DATA -->
    ! <!-- -------------- -->
    
    ! Underlying KKT system hierarchy that defines the shape of the solutions.
    type(t_kktsystemHierarchy), pointer :: p_rkktsystemHierarchy => null()
    
    ! Hierarchy of solutions (nonlinearities) for the nonlinear iteration
    type(t_kktsystemHierarchy), pointer :: p_rsolutionHierarchy => null()
        
    ! Hierarchy of directional derivatives. Calculated during the Newton iteration.
    type(t_kktsystemDirDerivHierarchy), pointer :: p_rdirDerivHierarchy => null()
    
    ! Descent direction
    type(t_controlSpace), pointer :: p_rdescentDir => null()
    
  end type

!</typeblock>

   public :: t_spacetimeNewton

!</types>
  
  ! Basic initialisation of the Newton solver
  public :: newtonit_init
  
  ! Structural initialisation
  public :: newtonit_initStructure
  
  ! Apply a Newton iteration in the control space
  public :: newtonit_solve
  
  ! Cleanup of structures
  public :: newtonit_doneStructure

  ! Final cleanup
  public :: newtonit_done
  
contains

  ! ***************************************************************************

!<subroutine>

  subroutine newtonit_init (rsolver,rsettingsSolver,ssection,rparamList)
  
!<description>
  ! Initialises the solver parameters according to a parameter list.
!</description>
  
!<input>
  ! Parameters of the OptFlow solver
  type(t_settings_optflow), intent(in), target :: rsettingsSolver
  
  ! Parameter list with the parameters configuring the nonlinear solver
  type(t_parlist), intent(in) :: rparamList

  ! Name of the section in the parameter list containing the parameters
  ! of the nonlinear solver.
  character(LEN=*), intent(in) :: ssection
!</input>

!<output>
  ! Solver structure receiving the parameters
  type(t_spacetimeNewton), intent(inout) :: rsolver
!</output>

!</subroutine>

    character(LEN=SYS_STRLEN) :: ssolverNonlin,ssolverLin
    character(LEN=SYS_STRLEN) :: ssolverSpaceForward,ssolverSpaceBackward
    character(LEN=SYS_STRLEN) :: ssolverSpaceForwardLin,ssolverSpaceBackwardLin

    ! Remember the solver settings for later use
    rsolver%p_rsettingsSolver => rsettingsSolver

    ! Initialise basic parameters
    call newtonit_initBasicParams (rsolver%rnewtonParams,ssection,rparamList)

    ! Get the sections with the parameters for the nonlinear / linear
    ! solver in space
    call parlst_getvalue_int (rparamList, ssection, "cspatialInitCondPolicy", &
        rsolver%cspatialInitCondPolicy, rsolver%cspatialInitCondPolicy)

    call parlst_getvalue_string (rparamList, ssection, &
        "ssectionNonlinSolverSpace", ssolverNonlin, "CC-NONLINEARSOLVER",bdequote=.true.)

    call parlst_getvalue_string (rparamList, ssection, &
        "ssectionLinSolverSpace", ssolverLin, "SPACETIME-LINSOLVER",bdequote=.true.)
        
    call parlst_getvalue_string (rparamList, ssection, &
        "ssectionLinSolverSpaceForward", ssolverSpaceForward, "CC-LINEARSOLVER",bdequote=.true.)

    call parlst_getvalue_string (rparamList, ssection, &
        "ssectionLinSolverSpaceBackward", ssolverSpaceBackward, "CC-LINEARSOLVER",bdequote=.true.)

    call parlst_getvalue_string (rparamList, ssection, &
        "ssectionLinSolverSpaceForwardLin", ssolverSpaceForwardlin, "CC-LINEARSOLVER",bdequote=.true.)

    call parlst_getvalue_string (rparamList, ssection, &
        "ssectionLinSolverSpaceBackwardLin", ssolverSpaceBackwardLin, "CC-LINEARSOLVER",bdequote=.true.)

    ! Create linear solvers in space for linera subproblems in the
    ! primal/dual space.

    allocate(rsolver%p_rlinsolHierPrimal)
    allocate(rsolver%p_rlinsolHierDual)
    allocate(rsolver%p_rlinsolHierPrimalLin)
    allocate(rsolver%p_rlinsolHierDualLin)

    ! Create solver structures for the same levels.

    ! Forward equation on all levels
    call lssh_createLinsolHierarchy (rsolver%p_rlinsolHierPrimal,&
        rsettingsSolver%rfeHierarchyPrimal,rsettingsSolver%rprjHierSpacePrimal,&
        1,0,rsettingsSolver%rphysics%cequation,rsettingsSolver%p_rparlist,&
        ssolverSpaceForward,rsettingsSolver%rdebugFlags)

    ! Backward equation on all levels
    call lssh_createLinsolHierarchy (rsolver%p_rlinsolHierDual,&
        rsettingsSolver%rfeHierarchyDual,rsettingsSolver%rprjHierSpaceDual,&
        1,0,rsettingsSolver%rphysics%cequation,rsettingsSolver%p_rparlist,&
        ssolverSpaceBackward,rsettingsSolver%rdebugFlags)

    ! Linearised forward equation on all levels
    call lssh_createLinsolHierarchy (rsolver%p_rlinsolHierPrimalLin,&
        rsettingsSolver%rfeHierarchyPrimal,rsettingsSolver%rprjHierSpacePrimal,&
        1,0,rsettingsSolver%rphysics%cequation,rsettingsSolver%p_rparlist,&
        ssolverSpaceForwardLin,rsettingsSolver%rdebugFlags)

    ! Linearised backward equation on all levels
    call lssh_createLinsolHierarchy (rsolver%p_rlinsolHierDualLin,&
        rsettingsSolver%rfeHierarchyDual,rsettingsSolver%rprjHierSpaceDual,&
        1,0,rsettingsSolver%rphysics%cequation,rsettingsSolver%p_rparlist,&
        ssolverSpaceBackwardLin,rsettingsSolver%rdebugFlags)
    
    ! Create the corresponding solver hierarchies.
    allocate(rsolver%p_rsolverHierPrimal)
    allocate(rsolver%p_rsolverHierDual)
    allocate(rsolver%p_rsolverHierPrimalLin)
    allocate(rsolver%p_rsolverHierDualLin)

    ! Forward equation. Created on all levels but only used on the highest one.
    caLL spaceslh_init (rsolver%p_rsolverHierPrimal,&
        OPTP_PRIMAL,rsolver%p_rlinsolHierPrimal,ssolverNonlin,rparamList)

    ! Backward equation, only linear. Created on all levels but only used on the highest one.
    ! Also fetch parameters of the nonlinear solver from the data file.
    ! The parameter are not used except for the output level, which determins
    ! the amoount of output of the solver.
    caLL spaceslh_init (rsolver%p_rsolverHierDual,&
        OPTP_DUAL,rsolver%p_rlinsolHierDual,ssolverNonlin,rparamList)

    ! The definition of the lineraised forward/backward equation depends upon
    ! whether we use the full Newton approach or not.
    select case (rsolver%rnewtonParams%ctypeIteration)

    ! --------------
    ! Partial Newton
    ! --------------
    case (1)
      ! Linearised forward equation, only linear. Used on all levels.
      ! Uses the same linear solver as the forward solver.
      caLL spaceslh_init (rsolver%p_rsolverHierPrimalLin,&
          OPTP_PRIMALLIN_SIMPLE,rsolver%p_rlinsolHierPrimalLin,ssolverNonlin,rparamList)

      ! Linearised forward equation, only linear. Used on all levels.
      ! Uses the same linear solver as the backward solver.#
      caLL spaceslh_init (rsolver%p_rsolverHierDualLin,&
          OPTP_DUALLIN_SIMPLE,rsolver%p_rlinsolHierDualLin,ssolverNonlin,rparamList)
    
    ! ----------------------------
    ! Full Newton, adaptive Newton
    ! ----------------------------
    case (2,3)
      ! Linearised forward equation, only linear. Used on all levels.
      ! Uses the same linear solver as the forward solver.
      caLL spaceslh_init (rsolver%p_rsolverHierPrimalLin,&
          OPTP_PRIMALLIN,rsolver%p_rlinsolHierPrimalLin,ssolverNonlin,rparamList)

      ! Linearised forward equation, only linear. Used on all levels.
      ! Uses the same linear solver as the backward solver.
      caLL spaceslh_init (rsolver%p_rsolverHierDualLin,&
          OPTP_DUALLIN,rsolver%p_rlinsolHierDualLin,ssolverNonlin,rparamList)
          
    case default
      call output_line ("Invalid nonlinear iteration",&
          OU_CLASS_ERROR,OU_MODE_STD,"newtonit_init")
      call sys_halt()
      
    end select

    ! Initialise the linear subsolver
    call newtonlin_init (rsolver%rlinsolParam,rsettingsSolver,&
        rsolver%p_rsolverHierPrimalLin,rsolver%p_rsolverHierDualLin,&
        rparamList,ssolverLin)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonit_getResidual (rsolver,rkktsystem,rresidual,dres,iresnorm,rstatistics)
  
!<description>
  ! Calculates the basic (unprecondiotioned) search direction of the 
  ! Newton algorithm.
!</description>
  
!<input>
  ! Type of norm. A LINALG_NORMxxxx constant.
  integer, intent(in) :: iresnorm
!</input>

!<inputoutput>
  ! Parameters for the Newton iteration.
  type(t_spacetimeNewton), intent(inout) :: rsolver

  ! Structure defining the KKT system. The control in this structure
  ! defines the current 'state' of the Newton algorithm.
  ! On output, the primal and dual solution are updated.
  type(t_kktsystem), intent(inout) :: rkktsystem
  
  ! On output, this structure receives a representation of the search
  ! direction / residual in the Newton iteration.
  type(t_controlSpace), intent(inout) :: rresidual

  ! Statistic structure which receives the statistics of the iteration.
  type(t_newtonitSolverStat), intent(inout) :: rstatistics
!</inputoutput>

!<output>
  ! L2-Norm of the residual
  real(DP), intent(out) :: dres
!</output>

!</subroutine>

    ! local variables
    type(t_vectorBlock), pointer :: p_rvector
    real(DP), dimension(:), pointer :: p_Ddata
    type(t_timer) :: rtimer
    
    call sptivec_getVectorFromPool(rkktsystem%p_rdualSol%p_rvectorAccess,2,p_rvector)
    call lsysbl_getbase_double (p_rvector,p_Ddata)

    ! -------------------------------------------------------------
    ! Step 1: Solve the primal and dual system.
    ! -------------------------------------------------------------

    if (rsolver%rnewtonParams%ioutputLevel .ge. 2) then
      call output_line ("Nonlin. space-time Residual: Solving the primal equation")
    end if

    ! Solve the primal equation, update the primal solution.
    call stat_clearTimer (rtimer)
    call stat_startTimer (rtimer)
    output_iautoOutputIndent = output_iautoOutputIndent + 2

    call kkt_solvePrimal (rkktsystem,rsolver%p_rsolverHierPrimal,&
        rsolver%cspatialInitCondPolicy,rstatistics%rspaceslSolverStat)

    output_iautoOutputIndent = output_iautoOutputIndent - 2
    call stat_stopTimer (rtimer)

    if (rsolver%rnewtonParams%ioutputLevel .ge. 3) then
      call output_line ("Nonlin. space-time Residual: Time for solving: "//&
          trim(sys_sdL(rtimer%delapsedReal,10)))
    end if
    
    ! Add time to the time of the forward equation
    call stat_addTimers (rtimer,rstatistics%rtimeForward)

    if (rsolver%rnewtonParams%ioutputLevel .ge. 2) then
      call output_line ("Nonlin. space-time Residual: Solving the dual equation")
    end if

    ! Solve the dual equation, update the dual solution.
    call stat_clearTimer (rtimer)
    call stat_startTimer (rtimer)
    output_iautoOutputIndent = output_iautoOutputIndent + 2

    call kkt_solveDual (rkktsystem,rsolver%p_rsolverHierDual,&
        rsolver%cspatialInitCondPolicy,rstatistics%rspaceslSolverStat)

    output_iautoOutputIndent = output_iautoOutputIndent - 2
    call stat_stopTimer (rtimer)
    
    ! Add time to the time of the backward equation
    call stat_addTimers (rtimer,rstatistics%rtimeBackward)

    if (rsolver%rnewtonParams%ioutputLevel .ge. 3) then
      call output_line ("Nonlin. space-time Residual: Time for solving: "//&
          trim(sys_sdL(rtimer%delapsedReal,10)))
    end if

    ! -------------------------------------------------------------
    ! Step 2: Calculate the search direction   
    ! -------------------------------------------------------------

    ! The search direction is just the residual in the control equation.
    call kkt_calcControlRes (rkktsystem,rresidual,dres,iresnorm)

  end subroutine

    ! ***************************************************************************

!<subroutine>

  subroutine newtonit_updateControl (rsolver,rkktsystem,rcorrection)
  
!<description>
  ! Calculates the basic (unprecondiotioned) search direction of the 
  ! Newton algorithm.
!</description>  

!<input>
  ! Parameters for the Newton iteration.
  ! The output parameters are changed according to the iteration.
  type(t_spacetimeNewton), intent(in) :: rsolver

  ! The preconditioned search direction
  type(t_controlSpace), intent(inout) :: rcorrection
!</input>
  
!<inputoutput>
  ! Structure defining the KKT system. The control in this structure
  ! is updated according to the search direction.
  type(t_kktsystem), intent(inout) :: rkktsystem
!</inputoutput>

!</subroutine>

    ! Currectly, this is just a linear combination of the control variables.
    !
    !    u_n+1  =  u_n + g_n
    !
    ! Later, a step length control can be added here.
    
    call kktsp_controlLinearComb (&
        rcorrection,rsolver%rnewtonParams%domega,rkktsystem%p_rcontrol,1.0_DP)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonit_adNewton_setEps (rsolver,radNewtonParams,dresInit,dresLastIte)
  
!<description>
  ! Realises the adaptive Newton algorithm. Sets the stopping criterions of
  ! all solvers to appropriate values.
  !
  ! WARNING!!! THIS ROUTINE HAS A SIDE EFFECT!
  ! IT SETS THE STOPPING CRITERION OF ALL SOLVERS IN SPACE TO AN APPROPRIATE
  ! VALUE! IF THE SPACE SOVLERS AE USED SOMEWHERE ELSE, THE STOPPING CRITERION
  ! IS LOST AND THUS, THEY MAY BEHAVE NOT AS EXPECTED!!!
!</description>

!<input>
  ! Parameters of the adaptive Newton algotithm.
  type(t_ccDynamicNewtonControl), intent(in) :: radNewtonParams
  
  ! Initial residual.
  ! May be set to 0.0 if there is no initial residual.
  real(DP), intent(in) :: dresInit
  
  ! Residual obtained in the last nonlinear iteration.
  ! In the first call, this should be set to dresInit.
  real(DP), intent(in) :: dresLastIte
!</input>

!<inputoutput>
  ! Parameters for the Newton iteration.
  ! The output parameters are changed according to the iteration.
  type(t_spacetimeNewton), intent(inout) :: rsolver
!</inputoutput>

!</subroutine>

    real(DP) :: depsAbs,depsRel,ddigitsGained, ddigitsToGain
    
    if ((dresInit .eq. 0.0_DP) .and. (dresLastIte .eq. 0.0_DP)) then
      
      ! Gain two digits for the initial residual, that is enough
      depsAbs = 0.0_DP
      depsRel = 1.0E-2_DP
      
      if (rsolver%rnewtonParams%ioutputLevel .ge. 2) then
        call output_line ("Adaptive Newton: New stopping criterion. ||res_rel|| < "//&
            trim(sys_sdEL(depsRel,10)))
      end if

    else
    
      ! We have to determine a new dresAbs for all solver components:
      ! - At least, gain as many digits as configured in the adaptive-Newton 
      !   structure has to be gained
      ! - The number of digits to gain in the next iteration has to be an
      !   appropriate multiple of the digits already gained.
      ! So...
      
      ddigitsGained = dresLastIte/dresInit
      
      ddigitsToGain = min(radNewtonParams%dinexactNewtonEpsRel*ddigitsGained,&
          ddigitsGained ** radNewtonParams%dinexactNewtonExponent)
          
      depsRel = 0.0_DP
      depsAbs = max(dresInit * ddigitsToGain,radNewtonParams%dinexactNewtonEpsAbs)
      
      if (rsolver%rnewtonParams%ioutputLevel .ge. 2) then
        call output_line ("Adaptive Newton: New stopping criterion. ||res|| < "//&
            trim(sys_sdEL(depsAbs,10)))
      end if
    end if
    
    ! Initialise the linear subsolver(s).
    call newtonlin_adNewton_setEps (rsolver%rlinsolParam,depsAbs,depsRel)
    
    ! Initialise the nonlinear and linear solver in space which are used
    ! for the calvulation of the current residual.
    !call spaceslh_setEps (rsolver%p_rsolverHierPrimal,depsAbs,depsRel)
    !call spaceslh_setEps (rsolver%p_rsolverHierDual,depsAbs,depsRel)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonit_solve (rsolver,rsolution,rstatistics)
  
!<description>
  ! Applies a Newton iteration to solve the space-time system.
!</description>
  
!<inputoutput>
  ! Parameters for the Newton iteration.
  ! The output parameters are changed according to the iteration.
  type(t_spacetimeNewton), intent(inout) :: rsolver

  ! Structure defining the solutions of the KKT system.
  type(t_kktsystem), intent(inout), target :: rsolution
!</inputoutput>

!<inputoutput>
  ! Statistic structure which receives the statistics of the iteration.
  type(t_newtonitSolverStat), intent(inout) :: rstatistics
!</inputoutput>

!</subroutine>
   
    ! local variables
    type(t_controlSpace), pointer :: p_rdescentDir
    type(t_kktsystem), pointer :: p_rsolution
    type(t_kktsystemDirDeriv), pointer :: p_rsolutionDirDeriv
    type(t_timer) :: rtotalTime
    type(t_newtonlinSolverStat) :: rstatisticsLinSol
    real(DP) :: delapsedReal
    
    ! Clear output statistics
    call newtonit_clearStatistics(rstatistics)
    
    ! Measure the total computational time
    call stat_startTimer(rtotalTime)
    
    ! Initialise data for the nonlinear iteration
    call newtonit_initData (rsolver,rsolution)

    ! Get a pointer to the solution and directional derivative on the maximum level
    call kkth_getKKTsystem (rsolver%p_rsolutionHierarchy,0,p_rsolution)
    call kkth_getKKTsystemDirDeriv (rsolver%p_rdirDerivHierarchy,0,p_rsolutionDirDeriv)
    
    ! Temporary vector(s)
    p_rdescentDir => rsolver%p_rdescentDir
    
    ! Prepare a structure that encapsules the directional derivative.
    
    ! Apply the Newton iteration
    rsolver%rnewtonParams%niterations = 0
    rsolver%rnewtonParams%dresInit = 0.0_DP
    rsolver%rnewtonParams%dresFinal = 0.0_DP
    
    do while (.true.)
    
      ! The Newton iteration reads
      !
      !    u_n+1  =  u_n  -  [J''(u_n)]^-1  J'(u_n)
      !
      ! or in other words,
      !
      !    u_n+1  =  u_n  +  [J''(u_n)]^-1  d_n
      !
      ! with the residual   d_n = -J'(u_n)   specifying a 'search direction'.
      
      ! -------------------------------------------------------------
      ! Get the current residual / search direction
      ! -------------------------------------------------------------

      if (rsolver%rnewtonParams%ioutputLevel .ge. 2) then
        call output_line ("Space-time Newton: Calculating the residual")
      end if

      ! Measure the time for the computation of the residual
      call stat_startTimer(rstatistics%rtimeDefect)

      ! Compute the basic (unpreconditioned) search direction d_n.
      output_iautoOutputIndent = output_iautoOutputIndent + 2
      
      call newtonit_getResidual (rsolver,p_rsolution,p_rdescentDir,&
          rsolver%rnewtonParams%dresFinal,rsolver%rnewtonParams%iresnorm,rstatistics)
      
      output_iautoOutputIndent = output_iautoOutputIndent - 2
      call stat_stopTimer(rstatistics%rtimeDefect)

      if (rsolver%rnewtonParams%niterations .eq. 0) then
        ! Remember the initial residual
        rsolver%rnewtonParams%dresInit = rsolver%rnewtonParams%dresFinal
      end if

      if (rsolver%rnewtonParams%ioutputLevel .ge. 2) then
        call output_line ("Space-time Newton: Iteration "// &
            trim(sys_siL(rsolver%rnewtonParams%niterations,10))// &
            ", ||res(u)|| = "// &
            trim(sys_sdEL(rsolver%rnewtonParams%dresFinal,10)))
      end if

      ! -------------------------------------------------------------
      ! Check for convergence
      ! -------------------------------------------------------------
      if (newtonit_checkConvergence (rsolver%rnewtonParams)) exit
      
      ! -------------------------------------------------------------
      ! Check for divergence
      ! -------------------------------------------------------------
      if (newtonit_checkDivergence (rsolver%rnewtonParams)) exit
      
      ! -------------------------------------------------------------
      ! Check other stopping criteria
      ! -------------------------------------------------------------
      if (newtonit_checkIterationStop (&
          rsolver%rnewtonParams,rsolver%rnewtonParams%niterations)) exit
      
      ! -------------------------------------------------------------
      ! Adaptive Newton for the next iteration?
      ! -------------------------------------------------------------
      if (rsolver%rnewtonParams%ctypeIteration .eq. 3) then
        call newtonit_adNewton_setEps (&
            rsolver,rsolver%rnewtonParams%radaptiveNewton,&
            rsolver%rnewtonParams%dresInit,rsolver%rnewtonParams%dresFinal)
      end if
    
      ! -------------------------------------------------------------
      ! Preconditioning with the Newton matrix
      ! -------------------------------------------------------------
      
      if (rsolver%rnewtonParams%ioutputLevel .ge. 2) then
        call output_line ("Space-time Newton: Preconditioning")
      end if

      ! Initialise data arrays in the linear subsolver.
      call stat_startTimer(rstatistics%rtimeProlRest)
      call newtonlin_initNonlinearData (rsolver%rlinsolParam,rsolver%p_rsolutionHierarchy)
      call stat_stopTimer(rstatistics%rtimeProlRest)
      
      call newtonlin_initData (rsolver%rlinsolParam,rsolver%p_rsolutionHierarchy)

      ! Actual Newton iteration. Apply the Newton preconditioner
      ! to get the Newton search direction:
      !
      !    J''(u_n) g_n  =  d_n
      !
      ! The control on the maximum level of p_rdirDerivHierarchy
      ! (identified by p_rsolutionDirDeriv%p_rcontrolLin) receives the result.
      call newtonlin_clearStatistics (rstatisticsLinSol)

      output_iautoOutputIndent = output_iautoOutputIndent + 2
      
      call newtonlin_precond (rsolver%rlinsolParam,&
          rsolver%p_rdirDerivHierarchy,p_rdescentDir,rstatisticsLinSol)
      
      output_iautoOutputIndent = output_iautoOutputIndent - 2
      
      ! Statistics
      call newtonlin_sumStatistics (&
          rstatisticsLinSol,rstatistics%rnewtonlinSolverStat,NLIN_STYPE_LINSOL)
      
      ! Clean up data in the linear subsolver
      call newtonlin_doneData (rsolver%rlinsolParam)

      ! -------------------------------------------------------------
      ! Update of the solution
      ! -------------------------------------------------------------

      ! Update the control according to the search direction:
      !
      !    u_n+1  =  u_n  +  g_n
      !
      ! or to any configured step-length control rule.
      call newtonit_updateControl (&
          rsolver,p_rsolution,p_rsolutionDirDeriv%p_rcontrolLin)

!    ! DEBUG!!!      
!    call kktsp_dualLinearComb (&
!        p_rsolutionDirDeriv%p_rdualSolLin,p_rsolution%p_rdualSol,1.0_DP,1.0_DP)
!    call kkt_calcControlRes (p_rsolution,p_rdescentDir,&
!        rsolver%rnewtonParams%dresFinal,rsolver%rnewtonParams%iresnorm)
!    !-> dres = u+du + 1/alpha (lambda+dlambda) =0 -> is ok!

      ! -------------------------------------------------------------
      ! Proceed with the next iteration
      ! -------------------------------------------------------------
      ! Next iteration
      rsolver%rnewtonParams%niterations = rsolver%rnewtonParams%niterations + 1
    
      if (rsolver%rnewtonParams%ioutputLevel .ge. 2) then
        call output_line ("Space-time Newton: Time for the linear solver = "//&
            sys_sdL(rstatisticsLinSol%rtotalTime%delapsedReal,10))
        call stat_sampleTimer(rtotalTime,delapsedReal)
        call output_line ("Space-time Newton: Computation time           = "//&
            sys_sdL(delapsedReal,10))
      end if

    end do
    
    ! Release data
    call newtonit_doneData (rsolver)
    
    call stat_stopTimer(rtotalTime)
    call stat_addTimers(rtotalTime,rstatistics%rtotalTime)
    
    rstatistics%niterations = rstatistics%niterations + rsolver%rnewtonParams%niterations
    
    ! Statistics
    if (rsolver%rnewtonParams%ioutputLevel .ge. 2) then
      call output_lbrk()
      call output_line ("Space-time Newton: Statistics.")
      call output_lbrk()
      call output_line ("#Nonlinear Iterations   : "//&
            trim(sys_siL(rsolver%rnewtonParams%niterations,10)) )
      call output_line ("#Iterations Precond.    : "//&
            trim(sys_siL(rsolver%rnewtonParams%nlinearIterations,10)) )
      call output_line ("!!INITIAL RES!!         : "//&
            trim(sys_sdEL(rsolver%rnewtonParams%dresInit,15)) )
      call output_line ("!!RES!!                 : "//&
            trim(sys_sdEL(rsolver%rnewtonParams%dresFinal,15)) )
      if (rsolver%rnewtonParams%dresInit .gt. &
          rsolver%rnewtonParams%drhsZero) then
        call output_line ("!!RES!!/!!INITIAL RES!! : "//&
          trim(sys_sdEL(rsolver%rnewtonParams%dresFinal / &
                        rsolver%rnewtonParams%dresInit,15)) )
      else
        call output_line ("!!RES!!/!!INITIAL RES!! : "//&
              trim(sys_sdEL(0.0_DP,15)) )
      end if

    end if
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonit_initStructure (rsolver,rkktsystemHierarchy,&
      rprjHierSpaceTimePrimal,rprjHierSpaceTimeDual,rprjHierSpaceTimeControl)
  
!<description>
  ! Structural initialisation of the Newton solver.
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
!</input>

!<inputoutput>
  ! Structure to be initialised.
  type(t_spacetimeNewton), intent(inout) :: rsolver
!</inputoutput>

!</subroutine>

    ! Remember the structure of the solutions.
    rsolver%p_rkktsystemHierarchy => rkktsystemHierarchy
    
    ! Initialise the structures of the linear subsolver.
    call newtonlin_initStructure (rsolver%rlinsolParam,rkktsystemHierarchy,&
        rprjHierSpaceTimePrimal,rprjHierSpaceTimeDual,rprjHierSpaceTimeControl)
   
  end subroutine


  ! ***************************************************************************

!<subroutine>

  subroutine newtonit_initData (rsolver,rsolution)
  
!<description>
  ! Final preparation of the Newton solver.
!</description>

!<input>
  ! Structure defining the solutions of the KKT system.
  type(t_kktsystem), intent(in), target :: rsolution
!</input>

!<inputoutput>
  ! Structure to be initialised.
  type(t_spacetimeNewton), intent(inout) :: rsolver
!</inputoutput>

!</subroutine>

    ! Allocate memory for the descent direction vector
    allocate(rsolver%p_rdescentDir)

    ! Create a vector holding the descent direction (highest level)
    call kkth_initControl (rsolver%p_rdescentDir,rsolver%p_rkktsystemHierarchy,0)
    
    ! Create temporary memory for the nonlinearity and the search direction.
    allocate(rsolver%p_rsolutionHierarchy)

    call kkth_initHierarchy (rsolver%p_rsolutionHierarchy,&
        rsolver%p_rkktsystemHierarchy%p_roperatorAsmHier,&
        rsolver%p_rkktsystemHierarchy%p_rspaceTimeHierPrimal,&
        rsolver%p_rkktsystemHierarchy%p_rspaceTimeHierDual,&
        rsolver%p_rkktsystemHierarchy%p_rspaceTimeHierControl,.true.,rsolution)

    ! create tempoprary memory for the search direction connected with
    ! the above solution hierarchy.
    allocate(rsolver%p_rdirDerivHierarchy)
    
    call kkth_initHierarchyDirDeriv (&
        rsolver%p_rdirDerivHierarchy,rsolver%p_rsolutionHierarchy)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonit_doneData (rsolver)
  
!<description>
  ! Cleanup of the data initalised in newtonit_initData.
!</description>

!<inputoutput>
  ! Structure to be cleaned up.
  type(t_spacetimeNewton), intent(inout) :: rsolver
!</inputoutput>

!</subroutine>

    ! Release the descent direction, memory for the directional derivative,...
    call kktsp_doneControlVector (rsolver%p_rdescentDir)
    deallocate(rsolver%p_rdescentDir)
    
    call kkth_doneHierarchyDirDeriv (rsolver%p_rdirDerivHierarchy)
    deallocate(rsolver%p_rdirDerivHierarchy)

    call kkth_doneHierarchy (rsolver%p_rsolutionHierarchy)
    deallocate(rsolver%p_rsolutionHierarchy)
   
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonit_doneStructure (rsolver)
  
!<description>
  ! Cleanup of the data initalised in newtonit_initStructure.
!</description>

!<inputoutput>
  ! Structure to be cleaned up.
  type(t_spacetimeNewton), intent(inout) :: rsolver
!</inputoutput>

!</subroutine>
   
    ! Release structures in the subsolver
    call newtonlin_doneStructure (rsolver%rlinsolParam)
   
    ! Detach the solution structure
    nullify(rsolver%p_rkktsystemHierarchy)
   
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonit_done (rsolver)
  
!<description>
  ! Clean up the Newton iteration.
!</description>

!<inputoutput>
  ! Structure to be cleaned up.
  type(t_spacetimeNewton), intent(inout) :: rsolver
!</inputoutput>

!</subroutine>

    ! Release the linear subsolver.
    call newtonlin_done (rsolver%rlinsolParam)

    ! Release the linear solvers in space.
    call spaceslh_done (rsolver%p_rsolverHierDualLin)
    call spaceslh_done (rsolver%p_rsolverHierPrimalLin)
    call spaceslh_done (rsolver%p_rsolverHierDual)
    call spaceslh_done (rsolver%p_rsolverHierPrimal)

    call lssh_releaseLinsolHierarchy (rsolver%p_rlinsolHierDualLin)
    call lssh_releaseLinsolHierarchy (rsolver%p_rlinsolHierPrimalLin)
    call lssh_releaseLinsolHierarchy (rsolver%p_rlinsolHierDual)
    call lssh_releaseLinsolHierarchy (rsolver%p_rlinsolHierPrimal)

    deallocate (rsolver%p_rsolverHierDualLin)
    deallocate (rsolver%p_rsolverHierPrimalLin)
    deallocate (rsolver%p_rsolverHierDual)
    deallocate (rsolver%p_rsolverHierPrimal)

    deallocate (rsolver%p_rlinsolHierDualLin)
    deallocate (rsolver%p_rlinsolHierPrimalLin)
    deallocate (rsolver%p_rlinsolHierDual)
    deallocate (rsolver%p_rlinsolHierPrimal)
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonit_clearStatistics(rstatistics)
  
!<description>
  ! Resets a statistic structure.
!</description>

!<inputoutput>
  ! Structure to be reset.
  type(t_newtonitSolverStat), intent(inout) :: rstatistics
!</inputoutput>

!</subroutine>

    rstatistics%niterations = 0
    
    call stat_clearTimer(rstatistics%rtotalTime)
    call stat_clearTimer(rstatistics%rtimeForward)
    call stat_clearTimer(rstatistics%rtimeBackward)
    call stat_clearTimer(rstatistics%rtimeDefect)
    call stat_clearTimer(rstatistics%rtimePostprocessing)
    
    call spacesl_clearStatistics (rstatistics%rspaceslSolverStat)
    call newtonlin_clearStatistics (rstatistics%rnewtonlinSolverStat)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonit_sumStatistics(rstatistics1,rstatistics2)
  
!<description>
  ! Sums up the data of rstatistics1 to the data in rstatistics2.
!</description>

!<input>
  ! Source structure
  type(t_newtonitSolverStat), intent(in) :: rstatistics1
!</input>

!<inputoutput>
  ! Destination structure.
  type(t_newtonitSolverStat), intent(inout) :: rstatistics2
!</inputoutput>

!</subroutine>

    rstatistics2%niterations = rstatistics2%niterations + rstatistics1%niterations
    
    call stat_addTimers(rstatistics1%rtotalTime,         rstatistics2%rtotalTime)
    call stat_addTimers(rstatistics1%rtimeForward,       rstatistics2%rtimeForward)
    call stat_addTimers(rstatistics1%rtimeBackward,      rstatistics2%rtimeBackward)
    call stat_addTimers(rstatistics1%rtimeDefect,        rstatistics2%rtimeDefect)
    call stat_addTimers(rstatistics1%rtimeProlRest,      rstatistics2%rtimeProlRest)
    call stat_addTimers(rstatistics1%rtimePostprocessing,rstatistics2%rtimePostprocessing)
    
    call spacesl_sumStatistics (&
        rstatistics1%rspaceslSolverStat,rstatistics2%rspaceslSolverStat)

    call newtonlin_sumStatistics (&
        rstatistics1%rnewtonlinSolverStat,rstatistics2%rnewtonlinSolverStat,NLIN_STYPE_LINSOL)

  end subroutine

end module
