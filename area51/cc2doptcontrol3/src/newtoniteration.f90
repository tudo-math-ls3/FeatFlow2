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
  use structuresnewton
  use assemblytemplates
  
  use spacematvecassembly
  
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
    
    ! <!-- -------------- -->
    ! <!-- TEMPORARY DATA -->
    ! <!-- -------------- -->
    
    ! Underlying KKT system hierarchy that defines the shape of the solutions.
    type(t_kktsystemHierarchy), pointer :: p_rkktsystemHierarchy => null()
    
    ! Descent direction
    type(t_controlSpace), pointer :: p_rdescentDir => null()
    
  end type

!</typeblock>

   public :: t_spacetimeNewton

!</types>
  
  ! Basic initialisation of the Newton solver
  public :: newtonit_initParams
  
  ! Structural initialisation
  public :: newtonit_initStructure
  
  ! Final initialisation
  public :: newtonit_initData
  
  ! Apply a Newton iteration in the control space
  public :: newtonit_solve
  
  ! Cleanup of data
  public :: newtonit_doneData

  ! Cleanup of structures
  public :: newtonit_doneStructure

  ! Final cleanup
  public :: newtonit_done
  
contains

  ! ***************************************************************************

!<subroutine>

  subroutine newtonit_initParams (rsolver,rsettingsSolver,ssection,rparamList)
  
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

    ! Remember the solver settings for later use
    rsolver%p_rsettingsSolver => rsettingsSolver

    ! Initialise basic parameters
    call newtonit_initBasicParams (rsolver%rnewtonParams,ssection,rparamList)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonit_getResidual (rkktsystem,rresidual,dres)
  
!<description>
  ! Calculates the basic (unprecondiotioned) search direction of the 
  ! Newton algorithm.
!</description>
  
!<inputoutput>
  ! Structure defining the KKT system. The control in this structure
  ! defines the current 'state' of the Newton algorithm.
  ! On output, the 
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
    call kkt_solvePrimal (rkktsystem)
    
    ! Solve the dual equation, update the dual solution.
    call kkt_solveDual (rkktsystem)

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

  subroutine newtonit_solve (rsolver,rkktsystemHierarchy)
  
!<inputoutput>
  ! Parameters for the Newton iteration.
  ! The output parameters are changed according to the iteration.
  type(t_spacetimeNewton), intent(inout) :: rsolver

  ! Structure defining a hierarchy of solutions of the KKT system.
  ! The solution on the maximum level is taken as initial
  ! values. On exit, the structure contains improved solutions
  ! on the maximum level. The content of the lower levels is
  ! undefined.
  type(t_kktsystemHierarchy), intent(inout), target :: rkktsystemHierarchy
!</inputoutput>

!</subroutine>
   
    ! local variables
    type(t_kktsystemDirDeriv) :: rkktsystemDirDeriv
    type(t_controlSpace), pointer :: p_rdescentDir
    type(t_kktsystem), pointer :: p_rkktsystem
    
    ! Get a pointer to the solution of the maximum level
    p_rkktsystem => rkktsystemHierarchy%p_RkktSystems(rkktsystemHierarchy%nlevels)
    
    ! Temporary vector
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
      call newtonit_getResidual (p_rkktsystem,p_rdescentDir,&
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
      call newtonlin_precondNewton (rsolver%rlinsolParam,&
          rkktsystemDirDeriv,p_rdescentDir)
      
      ! -------------------------------------------------------------
      ! Update of the solution
      ! -------------------------------------------------------------

      ! Update the control according to the search direction:
      !
      !    u_n+1  =  u_n  +  g_n
      !
      ! or to any configured step-length control rule.
      call newtonit_updateControl (rsolver,p_rkktsystem,p_rdescentDir)
      
      ! -------------------------------------------------------------
      ! Proceed with the next iteration
      ! -------------------------------------------------------------
      ! Next iteration
      rsolver%rnewtonParams%nnonlinearIterations = &
          rsolver%rnewtonParams%nnonlinearIterations + 1
    
    end do
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonit_init (rsolver,rparlist,ssection)
  
!<description>
  ! Basic initialisation of the Newton solver.
!</description>

!<input>
  ! Parameter list that contains the data for the preconditioning.
  type(t_parlist), intent(in) :: rparlist
  
  ! Entry Section in the parameter list containing the data of the 
  ! preconditioner.
  character(len=*), intent(in) :: ssection
!</input>

!<inputoutput>
  ! Structure to be initialised.
  type(t_spacetimeNewton), intent(out) :: rsolver
!</inputoutput>

!</subroutine>
   
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

    ! Remember the structure of the solutions
    rsolver%p_rkktsystemHierarchy => rkktsystemHierarchy
   
  end subroutine


  ! ***************************************************************************

!<subroutine>

  subroutine newtonit_initData (rsolver)
  
!<description>
  ! Final preparation of the Newton solver.
!</description>

!<inputoutput>
  ! Structure to be initialised.
  type(t_spacetimeNewton), intent(inout) :: rsolver
!</inputoutput>

!</subroutine>

    ! Allocate memory for the descent direction vector
    allocate(rsolver%p_rdescentDir)

    ! Create a vector holding the descent direction (highest level)
    call kkth_initControl (rsolver%p_rdescentDir,rsolver%p_rkktsystemHierarchy,0)
   
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

    ! Release the descent direction
    call kktsp_doneControlVector (rsolver%p_rdescentDir)
  
    deallocate(rsolver%p_rdescentDir)
   
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
   
  end subroutine

end module
