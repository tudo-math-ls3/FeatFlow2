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
  use assemblytemplates
  
  use spacematvecassembly
  use spacelinearsolver
  use spacesolver
  
  use kktsystemspaces
  use kktsystem
  
  use spacetimehierarchy
  
  use structuresnewton
  
  implicit none
  
  private

!<type>

  !<typeblock>
  
  ! Contains the parameters of the Newton iteration.
  type t_linsolParameters
  
    ! <!-- --------------------------------------- -->
    ! <!-- GENRERAL PARAMETERS AND SOLVER SETTINGS -->
    ! <!-- --------------------------------------- -->

    ! General newton parameters
    type(t_newtonPrecParameters) :: rprecParameters

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
    
    ! <!-- -------------- -->
    ! <!-- TEMPORARY DATA -->
    ! <!-- -------------- -->
    
    ! Temporary memory for calculations
    type(t_controlSpace), pointer :: p_rtempVector => null()

  end type
  
  !</typeblock>
  
  public :: t_linsolParameters

!</types>

  ! Apply preconditioning of a defect with the Newton preconditioner
  ! in the control space.
  public :: newtonlin_precond
  
  ! Basic initialisation of the Newton solver
  public :: newtonlin_init
  
  ! Structural initialisation
  public :: newtonlin_initStructure
  
  ! Final initialisation
  public :: newtonlin_initData
  
  ! Cleanup of data
  public :: newtonlin_doneData

  ! Cleanup of structures
  public :: newtonlin_doneStructure

  ! Final cleanup
  public :: newtonlin_done

contains

  ! ***************************************************************************

!<subroutine>

  subroutine newtonlin_getResidual (rlinsolParam,rkktsystemDirDeriv,rrhs,rresidual,dres)
  
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
    call kkt_solvePrimalDirDeriv (rkktsystemDirDeriv,&
        rlinsolParam%p_rsolverHierPrimalLin)
    
    ! Solve the dual equation, update the dual solution.
    call kkt_solveDualDirDeriv (rkktsystemDirDeriv,&
        rlinsolParam%p_rsolverHierDualLin)

    ! -------------------------------------------------------------
    ! Step 2: Calculate the residual
    ! -------------------------------------------------------------

    ! Take the solution of the linearised primal/dual system and
    ! calculate the residual in the control space.
    call kkt_calcControlResDirDeriv (rkktsystemDirDeriv,rrhs,rresidual,dres)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonlin_richardson (rlinsolParam,rkktsystemDirDeriv,rrhs)
  
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
  type(t_kktsystemDirDeriv), intent(inout), target :: rkktsystemDirDeriv
  
  ! Right-hand side of the linearised control equation
  type(t_controlSpace), intent(inout), target :: rrhs
!</inputoutput>

!</subroutine>

    ! Richardson is a do loop...
    rlinsolParam%rprecParameters%niterations = 0

    do while (.true.)
    
      ! -------------------------------------------------------------
      ! Get the current residual / search direction
      ! -------------------------------------------------------------

      ! Compute the residual and its norm.
      call newtonlin_getResidual (rlinsolParam,rkktsystemDirDeriv,rrhs,&
          rlinsolParam%p_rtempVector,rlinsolParam%rprecParameters%dresFinal)

      if (rlinsolParam%rprecParameters%niterations .eq. 1) then
        ! Remember the initial residual
        rlinsolParam%rprecParameters%dresInit = rlinsolParam%rprecParameters%dresFinal
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
      call kktsp_controlLinearComb (rlinsolParam%p_rtempVector,&
          rlinsolParam%rprecParameters%domega,&
          rkktsystemDirDeriv%p_rcontrolLin,1.0_DP)
    
      ! -------------------------------------------------------------
      ! Proceed with the iteration
      ! -------------------------------------------------------------
      ! Next iteration
      rlinsolParam%rprecParameters%niterations = &
          rlinsolParam%rprecParameters%niterations + 1
    
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonlin_precond (rlinsolParam,rkktsystemDirDeriv,rnewtonDir)
  
!<description>
  ! Calculates the Newton search direction by applying the Newton
  ! preconditioner to the given residual rnewtonDir.
  ! This calculates a correction g with
  !      J''(u) g = d
  ! where d=rnewtonDir. On return of this routine, there is rnewtonDir=g.
!</description>

!<inputoutput>
  ! Parameters for the iteration.
  ! The output parameters are changed according to the iteration.
  type(t_linsolParameters), intent(inout) :: rlinsolParam

  ! Structure defining directional derivatives of the KKT system.
  ! Used for temporary calculations.
  type(t_kktsystemDirDeriv), intent(inout) :: rkktsystemDirDeriv 

  ! This structure contains the search direction to be
  ! preconditioned. Will be replaced by the precoditioned search direction.
  type(t_controlSpace), intent(inout) :: rnewtonDir
!</inputoutput>

!</subroutine>

    integer :: ispacelevel

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
    
    call kktsp_clearPrimal (rkktsystemDirDeriv%p_rprimalSolLin)
    call kktsp_clearDual (rkktsystemDirDeriv%p_rdualLinSol)
    call kktsp_clearControl (rkktsystemDirDeriv%p_rcontrolLin)
    
    ! Call the Richardson iteration to calculate an update
    ! for the control.
    call newtonlin_richardson (rlinsolParam,rkktsystemDirDeriv,rnewtonDir)
    
    ! Overwrite the rnewtonDir with the update.
    call kktsp_controlLinearComb (&
        rkktsystemDirDeriv%p_rcontrolLin,1.0_DP,rnewtonDir,0.0_DP)

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
   
    ! Remember the settings for later use
    rlinsolParam%p_rsettingsSolver => rsettingsSolver
    rlinsolParam%p_rsolverHierPrimalLin => rsolverHierPrimalLin
    rlinsolParam%p_rsolverHierDualLin => rsolverHierDualLin
   
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonlin_initStructure (rlinsolParam)
  
!<description>
  ! Structural initialisation of the preconditioner
!</description>

!<inputoutput>
  ! Structure to be initialised.
  type(t_linsolParameters), intent(inout) :: rlinsolParam
!</inputoutput>

!</subroutine>
   
  end subroutine


  ! ***************************************************************************

!<subroutine>

  subroutine newtonlin_initData (rlinsolParam)
  
!<description>
  ! Final preparation of the preconditioner for preconditioning.
!</description>

!<inputoutput>
  ! Structure to be initialised.
  type(t_linsolParameters), intent(inout) :: rlinsolParam
!</inputoutput>

!</subroutine>
   
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonlin_doneData (rlinsolParam)
  
!<description>
  ! Cleanup of the data initalised in newtonlin_initData.
!</description>

!<inputoutput>
  ! Structure to be cleaned up.
  type(t_linsolParameters), intent(inout) :: rlinsolParam
!</inputoutput>

!</subroutine>
   
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonlin_doneStructure (rlinsolParam)
  
!<description>
  ! Cleanup of the data initalised in newtonlin_initStructure.
!</description>

!<inputoutput>
  ! Structure to be cleaned up.
  type(t_linsolParameters), intent(inout) :: rlinsolParam
!</inputoutput>

!</subroutine>
   
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonlin_done (rlinsolParam)
  
!<description>
  ! Clean up a preconditioner.
!</description>

!<inputoutput>
  ! Structure to be cleaned up.
  type(t_linsolParameters), intent(inout) :: rlinsolParam
!</inputoutput>

!</subroutine>
   
  end subroutine

end module
