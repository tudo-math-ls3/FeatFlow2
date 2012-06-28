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
  
  ! Parameters for the space-time Newton algorithm.
  type t_spacetimeNewton

    ! <!-- --------------------------------------- -->
    ! <!-- GENRERAL PARAMETERS AND SOLVER SETTINGS -->
    ! <!-- --------------------------------------- -->

    ! General newton parameters
    type(t_newtonParameters) :: rnewtonParams

    ! Parameters of the OptFlow solver
    type(t_settings_optflow), pointer :: p_rsettingsSolver => null()

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
        ssolverSpaceForward,rsettingsSolver%rdebugFlags)

    ! Linearised backward equation on all levels
    call lssh_createLinsolHierarchy (rsolver%p_rlinsolHierDualLin,&
        rsettingsSolver%rfeHierarchyDual,rsettingsSolver%rprjHierSpaceDual,&
        1,0,rsettingsSolver%rphysics%cequation,rsettingsSolver%p_rparlist,&
        ssolverSpaceBackward,rsettingsSolver%rdebugFlags)
    
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

  subroutine newtonit_getResidual (rsolver,rkktsystem,rresidual,dres,iresnorm)
  
!<description>
  ! Calculates the basic (unprecondiotioned) search direction of the 
  ! Newton algorithm.
!</description>
  
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

  ! Type of norm. A LINALG_NORMxxxx constant.
  integer, intent(in) :: iresnorm
!</inputoutput>

!<output>
  ! L2-Norm of the residual
  real(DP), intent(out) :: dres
!</output>

!</subroutine>
    type(t_vectorBlock), pointer :: p_rvector
    real(DP), dimension(:), pointer :: p_Ddata
    
    call sptivec_getVectorFromPool(rkktsystem%p_rdualSol%p_rvectorAccess,2,p_rvector)
    call lsysbl_getbase_double (p_rvector,p_Ddata)

    ! -------------------------------------------------------------
    ! Step 1: Solve the primal and dual system.
    ! -------------------------------------------------------------

    if (rsolver%rnewtonParams%ioutputLevel .ge. 2) then
      call output_line ("  Nonlin. space-time Residual: Solving the primal equation")
    end if

    ! Solve the primal equation, update the primal solution.
    call kkt_solvePrimal (rkktsystem,rsolver%p_rsolverHierPrimal,&
        rsolver%cspatialInitCondPolicy)

    if (rsolver%rnewtonParams%ioutputLevel .ge. 2) then
      call output_line ("  Nonlin. space-time Residual: Solving the dual equation")
    end if

    ! Solve the dual equation, update the dual solution.
    call kkt_solveDual (rkktsystem,rsolver%p_rsolverHierDual,&
        rsolver%cspatialInitCondPolicy)

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

  subroutine newtonit_solve (rsolver,rsolution)
  
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

!</subroutine>
   
    ! local variables
    type(t_controlSpace), pointer :: p_rdescentDir
    type(t_kktsystem), pointer :: p_rsolution
    type(t_kktsystemDirDeriv), pointer :: p_rsolutionDirDeriv
    
    ! Initialise data for the nonlinear iteration
    call newtonit_initData (rsolver,rsolution)

    ! Get a pointer to the solution and directional derivative on the maximum level
    call kkth_getKKTsystem (rsolver%p_rsolutionHierarchy,0,p_rsolution)
    call kkth_getKKTsystemDirDeriv (rsolver%p_rdirDerivHierarchy,0,p_rsolutionDirDeriv)
    
    ! Temporary vector(s)
    p_rdescentDir => rsolver%p_rdescentDir
    
    ! Prepare a structure that encapsules the directional derivative.
    
    ! Apply the Newton iteration
    rsolver%rnewtonParams%nnonlinearIterations = 0
    
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

      ! Compute the basic (unpreconditioned) search direction d_n.
      call newtonit_getResidual (rsolver,p_rsolution,p_rdescentDir,&
          rsolver%rnewtonParams%dresFinal,rsolver%rnewtonParams%iresnorm)

      if (rsolver%rnewtonParams%nnonlinearIterations .eq. 0) then
        ! Remember the initial residual
        rsolver%rnewtonParams%dresInit = rsolver%rnewtonParams%dresFinal
      end if

      if (rsolver%rnewtonParams%ioutputLevel .ge. 2) then
        call output_line ("Space-time Newton: Iteration "// &
            trim(sys_siL(rsolver%rnewtonParams%nnonlinearIterations,10))// &
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
      if (newtonit_checkIterationStop (rsolver%rnewtonParams)) exit
      
      ! -------------------------------------------------------------
      ! Preconditioning with the Newton matrix
      ! -------------------------------------------------------------
      
      if (rsolver%rnewtonParams%ioutputLevel .ge. 2) then
        call output_line ("Space-time Newton: Preconditioning")
      end if

      ! Initialise data arrays in the linear subsolver.
      call newtonlin_initNonlinearData (rsolver%rlinsolParam,rsolver%p_rsolutionHierarchy)
      call newtonlin_initData (rsolver%rlinsolParam,rsolver%p_rsolutionHierarchy)

      ! Actual Newton iteration. Apply the Newton preconditioner
      ! to get the Newton search direction:
      !
      !    J''(u_n) g_n  =  d_n
      !
      ! The control on the maximum level of p_rdirDerivHierarchy
      ! (identified by p_rsolutionDirDeriv%p_rcontrolLin) receives the result.
      output_iautoOutputIndent = output_iautoOutputIndent + 2
      call newtonlin_precond (rsolver%rlinsolParam,&
          rsolver%p_rdirDerivHierarchy,p_rdescentDir)
      output_iautoOutputIndent = output_iautoOutputIndent - 2
      
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
      rsolver%rnewtonParams%nnonlinearIterations = &
          rsolver%rnewtonParams%nnonlinearIterations + 1
    
    end do
    
    ! Release data
    call newtonit_doneData (rsolver)
    
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

end module
