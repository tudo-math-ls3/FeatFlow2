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

    ! Parameters for the linear space-time subsolver.
    type(t_linsolParameters) :: rlinsolParam
    
    ! Hierarchy of linear solvers in space for all levels.
    ! Nonlinear primal equation.
    type(t_linsolHierarchySpace), pointer :: p_rlinsolHierPrimal => null()

    ! Hierarchy of linear solvers in space for all levels.
    ! Dual equation.
    type(t_linsolHierarchySpace), pointer :: p_rlinsolHierDual => null()

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

    integer :: cequation
    character(LEN=SYS_STRLEN) :: ssolverNonlin, ssolverLin

    ! Remember the solver settings for later use
    rsolver%p_rsettingsSolver => rsettingsSolver

    ! Initialise basic parameters
    call newtonit_initBasicParams (rsolver%rnewtonParams,ssection,rparamList)

    ! Get the sections with the parameters for the nonlinear / linear
    ! solver in space
    call parlst_getvalue_string (rparamList, ssection, &
        "ssectionNonlinSolverSpace", ssolverNonlin, "CC-NONLINEARSOLVER",bdequote=.true.)
        
    call parlst_getvalue_string (rparamList, ssection, &
        "ssectionLinSolverSpace", ssolverLin, "CC-LINEARSOLVER",bdequote=.true.)

    ! Type of equation?
    cequation = LSS_EQN_GENERAL
    select case (rsettingsSolver%rphysics%cequation)
    case (0,1)
      ! Stokes / Navier-Stokes
      cequation = LSS_EQN_STNAVST2D
    end select
    
    ! Create linear solvers in space for linera subproblems in the
    ! primal/dual space.
    !
    ! Create solver structures for the same levels.
    !
    ! Forward equation on all levels
    allocate(rsolver%p_rlinsolHierPrimal)
    call lssh_createLinsolHierarchy (rsolver%p_rlinsolHierPrimal,&
        rsettingsSolver%rfeHierarchyPrimal,rsettingsSolver%rprjHierSpacePrimal,&
        1,0,cequation,rsettingsSolver%p_rparlist,ssolverLin,rsettingsSolver%rdebugFlags)

    ! Backward equation on all levels
    allocate(rsolver%p_rlinsolHierDual)
    call lssh_createLinsolHierarchy (rsolver%p_rlinsolHierDual,&
        rsettingsSolver%rfeHierarchyDual,rsettingsSolver%rprjHierSpaceDual,&
        1,0,cequation,rsettingsSolver%p_rparlist,ssolverLin,rsettingsSolver%rdebugFlags)
    
    ! Create the corresponding solver hierarchies.
    !
    ! Forward equation. Created on all levels but only used on the highest one.
    allocate(rsolver%p_rsolverHierPrimal)
    caLL spaceslh_init (rsolver%p_rsolverHierPrimal,rsettingsSolver,&
        OPTP_PRIMAL,rsolver%p_rlinsolHierPrimal,rsettingsSolver%roptcBDCSpaceHierarchy,&
        ssolverNonlin,rparamList)

    ! Backward equation, only linear. Created on all levels but only used on the highest one.
    allocate(rsolver%p_rsolverHierDual)
    caLL spaceslh_init (rsolver%p_rsolverHierDual,rsettingsSolver,&
        OPTP_DUAL,rsolver%p_rlinsolHierDual,rsettingsSolver%roptcBDCSpaceHierarchy)

    ! Linearised forward equation, only linear. Used on all levels.
    ! Uses the same linear solver as the forward solver.
    allocate(rsolver%p_rsolverHierPrimalLin)
    caLL spaceslh_init (rsolver%p_rsolverHierPrimalLin,rsettingsSolver,&
        OPTP_PRIMALLIN,rsolver%p_rlinsolHierPrimal,rsettingsSolver%roptcBDCSpaceHierarchy)

    ! Linearised forward equation, only linear. Used on all levels.
    ! Uses the same linear solver as the backward solver.#
    allocate(rsolver%p_rsolverHierDualLin)
    caLL spaceslh_init (rsolver%p_rsolverHierDualLin,rsettingsSolver,&
        OPTP_DUALLIN,rsolver%p_rlinsolHierDual,rsettingsSolver%roptcBDCSpaceHierarchy)

    ! Initialise the linear subsolver
    call newtonlin_init (rsolver%rlinsolParam,rsettingsSolver,&
        rsolver%p_rsolverHierPrimalLin,rsolver%p_rsolverHierDualLin,rparamList,ssolverLin)
        
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonit_getResidual (rsolver,rkktsystem,rresidual,dres)
  
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
!</inputoutput>

!<output>
  ! L2-Norm of the residual
  real(DP), intent(out) :: dres
!</output>

!</subroutine>

    ! -------------------------------------------------------------
    ! Step 1: Solve the primal and dual system.
    ! -------------------------------------------------------------

    ! Solve the primal equation, update the primal solution.
    call kkt_solvePrimal (rkktsystem,rsolver%p_rsolverHierPrimal)

    ! Solve the dual equation, update the dual solution.
    call kkt_solveDual (rkktsystem,rsolver%p_rsolverHierDual)

    ! -------------------------------------------------------------
    ! Step 2: Calculate the search direction   
    ! -------------------------------------------------------------

    ! The search direction is just the residual in the control equation.
    call kkt_calcControlRes (rkktsystem,rresidual,dres)

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

      ! Compute the basic (unpreconditioned) search direction d_n.
      call newtonit_getResidual (rsolver,p_rsolution,p_rdescentDir,&
          rsolver%rnewtonParams%dresFinal)

      if (rsolver%rnewtonParams%nnonlinearIterations .eq. 1) then
        ! Remember the initial residual
        rsolver%rnewtonParams%dresInit = rsolver%rnewtonParams%dresFinal
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
      
      ! Actual Newton iteration. Apply the Newton preconditioner
      ! to get the Newton search direction:
      !
      !    J''(u_n) g_n  =  d_n
      !
      ! The control on the maximum level of p_rdirDerivHierarchy
      ! (identified by p_rsolutionDirDeriv%p_rcontrolLin) receives the result.
      call newtonlin_precond (rsolver%rlinsolParam,&
          rsolver%p_rdirDerivHierarchy,p_rdescentDir)
      
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

  subroutine newtonit_initStructure (rsolver,rkktsystemHierarchy)
  
!<description>
  ! Structural initialisation of the Newton solver.
!</description>

!<input>
  ! Defines the basic hierarchy of the solutions of the KKT system.
  ! This can be a 'template' structure, i.e., memory for the solutions
  ! in rkktsystemHierarchy does not have to be allocated.
  type(t_kktsystemHierarchy), intent(in), target :: rkktsystemHierarchy
!</input>

!<inputoutput>
  ! Structure to be initialised.
  type(t_spacetimeNewton), intent(inout) :: rsolver
!</inputoutput>

!</subroutine>

    ! Remember the structure of the solutions.
    rsolver%p_rkktsystemHierarchy => rkktsystemHierarchy
    
    ! Initialise the structures of the linear subsolver.
    call newtonlin_initStructure (rsolver%rlinsolParam,rkktsystemHierarchy)
   
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

    call lssh_releaseLinsolHierarchy (rsolver%p_rlinsolHierDual)
    call lssh_releaseLinsolHierarchy (rsolver%p_rlinsolHierPrimal)

    deallocate (rsolver%p_rsolverHierDualLin)
    deallocate (rsolver%p_rsolverHierPrimalLin)
    deallocate (rsolver%p_rsolverHierDual)
    deallocate (rsolver%p_rsolverHierPrimal)

    deallocate (rsolver%p_rlinsolHierDual)
    deallocate (rsolver%p_rlinsolHierPrimal)
    
  end subroutine

end module
