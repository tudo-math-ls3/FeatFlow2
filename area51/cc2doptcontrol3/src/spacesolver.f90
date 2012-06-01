!##############################################################################
!# ****************************************************************************
!# <name> spacesolver </name>
!# ****************************************************************************
!#
!# <purpose>
!# Contains routines for solving nonlinear/linear systems in space.
!#
!# For a forward equation, the space solver applies a nonlinear loop.
!# For a backward equation or a linearised forward/backward equation,
!# the space solver directly delegates the solution task to a linear
!# subsolver.
!# </purpose>
!##############################################################################

module spacesolver

  use fsystem
  use genoutput
  
  use spatialdiscretisation
  use timediscretisation
  use linearalgebra
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
  use structuresboundaryconditions
  use assemblytemplates
  
  use spatialbc
  
  use structuresoperatorasm
  use spacematvecassembly
  
  use kktsystemspaces
  
  use spacelinearsolver
  use structuresnewton

  implicit none
  
  private

!<types>

!<typeblock>

  ! This structure encapsules a hierarchy of linear solvers for all
  ! levels in a space discrtisation hierarchy.
  type t_spaceSolverHierarchy
  
    ! <!-- --------------------------- -->
    ! <!-- NONLINEAR SOLVER PARAMETERS -->
    ! <!-- --------------------------- -->

    ! General Newton parameters for nonlinear iterations.
    ! This is level independent.
    type(t_newtonParameters) :: rnewtonParams

    ! Temporary vectors for nonlinear iterations
    type(t_vectorBlock), pointer :: p_rd => null()
    type(t_vectorBlock), pointer :: p_rx => null()

    ! <!-- ------------------------ -->
    ! <!-- LINEAR SOLVER STRUCTURES -->
    ! <!-- ------------------------ -->

    ! Hierarchy of linear solvers in space used for solving
    ! auxiliary linear subproblems. The solver on level ilevel
    ! will be used to solve linear subproblems.
    type(t_linsolHierarchySpace), pointer :: p_rlssHierarchy => null()

    ! <!-- ------------------------------ -->
    ! <!-- BOUNDARY CONDITIONS STRUCTURES -->
    ! <!-- ------------------------------ -->

    ! Boundary condition hierarchy for all space levels, primal and dual space
    type(t_optcBDCSpaceHierarchy), pointer :: p_roptcBDCSpaceHierarchy => null()

    ! <!-- ---------------- -->
    ! <!-- GENERAL SETTINGS -->
    ! <!-- ---------------- -->

    ! Type of equation to be solved here. This can be 
    ! OPTP_PRIMAL, OPTP_DUAL, OPTP_PRIMALLIN or OPTP_DUALLIN,
    ! depending on which equation to solve.
    integer :: copType = 0

    ! Parameters of the OptFlow solver
    type(t_settings_optflow), pointer :: p_rsettingsSolver => null()

  end type

!</typeblock>

  public :: t_spaceSolverHierarchy

!</types>

  ! Basic initialisation of a space solver
  public :: spaceslh_init

  ! Structural initialisation
  public :: spaceslh_initStructure
  
  ! Final initialisation
  public :: spaceslh_initData
  
  ! Solves the spatial system
  public :: spaceslh_solve
  
  ! Cleanup of data
  public :: spaceslh_doneData

  ! Cleanup of structures
  public :: spaceslh_doneStructure

  ! Final cleanup
  public :: spaceslh_done

contains

  ! ***************************************************************************

!<subroutine>

  subroutine spaceslh_init (rsolver,rsettingsSolver,copType,&
      rlssHierarchy,roptcBDCSpaceHierarchy,ssection,rparamList)
  
!<description>
  ! Initialises the solver parameters according to a parameter list.
  ! For each level in rlssHierarchy, a space solver will be created.
!</description>
  
!<input>
  ! Parameters of the OptFlow solver
  type(t_settings_optflow), intent(in), target :: rsettingsSolver
  
  ! Type of equation to be solved here. This can be 
  ! OPTP_PRIMAL, OPTP_DUAL, OPTP_PRIMALLIN or OPTP_DUALLIN,
  ! depending on which equation to solve.
  integer, intent(in) :: copType

  ! Hierarchy of linear solvers in space used for solving
  ! auxiliary linear subproblems. The solver on level ilevel
  ! will be used to solve linear subproblems.
  type(t_linsolHierarchySpace), target :: rlssHierarchy

  ! Boundary condition hierarchy for all space levels, primal and dual space
  type(t_optcBDCSpaceHierarchy), target :: roptcBDCSpaceHierarchy

  ! OPTIONAL: Name of the section in the parameter list containing the parameters
  ! of the nonlinear solver.
  ! Can be omitted for linear solvers in space.
  character(LEN=*), intent(in), optional :: ssection

  ! OPTIONAL: Parameter list with the parameters configuring the nonlinear solver
  ! Can be omitted for linear solvers in space.
  type(t_parlist), intent(in), optional :: rparamList
!</input>

!<output>
  ! Solver structure receiving the parameters
  type(t_spaceSolverHierarchy), intent(inout) :: rsolver
!</output>

!</subroutine>

    ! Remember the solver settings for later use
    rsolver%p_rsettingsSolver => rsettingsSolver
    rsolver%p_rlssHierarchy => rlssHierarchy
    rsolver%p_roptcBDCSpaceHierarchy => roptcBDCSpaceHierarchy
    rsolver%copType = copType

    if (present(ssection) .and. present(rparamList)) then
    
      ! Initialise basic parameters of the nonlinear solver
      call newtonit_initBasicParams (rsolver%rnewtonParams,ssection,rparamList)
    
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spaceslh_initStructure (rsolver, ilevel, ierror)
  
!<description>
  ! Structural initialisation of the Newton solver.
!</description>

!<input>
  ! Level in the global space hierarchy, the solver should be applied to.
  integer, intent(in) :: ilevel
!</input>

!<inputoutput>
  ! Structure to be initialised.
  type(t_spaceSolverHierarchy), intent(inout) :: rsolver
!</inputoutput>

!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  integer, intent(out) :: ierror
!</output>

!</subroutine>

    ! Initialise the structures of the associated linear subsolver
    call lssh_initStructure (&
        rsolver%p_rlsshierarchy,ilevel,ierror)

  end subroutine


  ! ***************************************************************************

!<subroutine>

  subroutine spaceslh_initData (rsolver, ilevel, cflags, ierror)
  
!<description>
  ! Final preparation of the Newton solver.
!</description>

!<input>
  ! Level in the global space hierarchy, the solver should be applied to.
  integer, intent(in) :: ilevel

  ! A combination of LSS_SLFLAGS_xxxx constants defining the way,
  ! the solver performs filtering.
  integer(I32), intent(in) :: cflags
!</input>

!<inputoutput>
  ! Structure to be initialised.
  type(t_spaceSolverHierarchy), intent(inout) :: rsolver
!</inputoutput>

!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  integer, intent(out) :: ierror
!</output>

!</subroutine>

    ! Initialise the structures of the associated linear subsolver
    call lssh_initData (&
        rsolver%p_rlsshierarchy,ilevel,cflags,ierror)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spaceslh_doneData (rsolver, ilevel)
  
!<description>
  ! Cleanup of the data initalised in spaceslh_initData.
!</description>

!<input>
  ! Level in the global space hierarchy to be cleaned up.
  integer, intent(in) :: ilevel
!</input>

!<inputoutput>
  ! Structure to be cleaned up.
  type(t_spaceSolverHierarchy), intent(inout) :: rsolver
!</inputoutput>

!</subroutine>

    ! Clean up the structures of the associated linear subsolver
    call lssh_doneData (rsolver%p_rlsshierarchy,ilevel)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spaceslh_doneStructure (rsolver,ilevel)
  
!<description>
  ! Cleanup of the data initalised in spaceslh_initStructure.
!</description>

!<input>
  ! Level in the global space hierarchy to be cleaned up.
  integer, intent(in) :: ilevel
!</input>

!<inputoutput>
  ! Structure to be cleaned up.
  type(t_spaceSolverHierarchy), intent(inout) :: rsolver
!</inputoutput>

!</subroutine>
   
    ! Clean up the structures of the associated linear subsolver
    call lssh_doneStructure (rsolver%p_rlsshierarchy,ilevel)
   
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spaceslh_done (rsolver)
  
!<description>
  ! Clean up the Newton iteration.
!</description>

!<inputoutput>
  ! Structure to be cleaned up.
  type(t_spaceSolverHierarchy), intent(inout) :: rsolver
!</inputoutput>

!</subroutine>

    type (t_spaceSolverHierarchy) :: rtemplate

    ! Initialise with standard parameters    
    rsolver = rtemplate

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spaceslh_solve (rsolver,ilevel,rd)
  
!<description>
  ! Solves the spatial linear/nonlinear system.
!</description>

!<input>
  ! Level in the global space hierarchy corresponding to rd.
  integer, intent(in) :: ilevel
!</input>

!<inputoutput>
  ! Parameters for the Newton iteration.
  ! The output parameters are changed according to the iteration.
  type(t_spaceSolverHierarchy), intent(inout) :: rsolver

  ! Vector to be preconditioned. Used as right-hand side.
  ! Will be replaced by the solution.
  type(t_vectorBlock), intent(inout), target :: rd
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_vectorBlock), pointer :: p_rx,p_rd
   
    ! Which type of operator is to be solved?
    select case (rsolver%coptype)

    ! ---------------------------------------------------------------
    ! Primal equation. Nonlinear loop.
    ! ---------------------------------------------------------------
    case (OPTP_PRIMAL)

      ! This is the most complicated part. On level ilevel, we have
      ! to apply a nonlinear loop.
      !
      ! At first, get a temp vector we can use for creating the
      ! nonlinear defect.
      p_rd => rsolver%p_rd
      p_rx => rsolver%p_rx
      
      ! ...
      
      
      
      

    ! ---------------------------------------------------------------
    ! Dual equation. Linear loop.
    ! ---------------------------------------------------------------
    case (OPTP_DUAL)
      ! Call the linear solver, it does the job for us.
      call lssh_precondDefect (rsolver%p_rlsshierarchy,ilevel,rd)
    
    ! ---------------------------------------------------------------
    ! Linearised primal equation. Linear loop.
    ! ---------------------------------------------------------------
    case (OPTP_PRIMALLIN,OPTP_PRIMALLIN_SIMPLE)
      ! Call the linear solver, it does the job for us.
      call lssh_precondDefect (rsolver%p_rlsshierarchy,ilevel,rd)

    ! ---------------------------------------------------------------
    ! Linearised dual equation. Linear loop.
    ! ---------------------------------------------------------------
    case (OPTP_DUALLIN)
      ! Call the linear solver, it does the job for us.
      call lssh_precondDefect (rsolver%p_rlsshierarchy,ilevel,rd)

    end select

!
!
!
!    ! local variables
!    type(t_kktsystemDirDeriv) :: rkktsystemDirDeriv
!    type(t_controlSpace), pointer :: p_rdescentDir
!    type(t_kktsystem), pointer :: p_rkktsystem
!    
!    ! Get a pointer to the solution of the maximum level
!    p_rkktsystem => rkktsystemHierarchy%p_RkktSystems(rkktsystemHierarchy%nlevels)
!    
!    ! Temporary vector
!    p_rdescentDir => rsolver%p_rdescentDir
!    
!    ! Prepare a structure that encapsules the directional derivative.
!    
!    ! Apply the Newton iteration
!    rsolver%rnewtonParams%nnonlinearIterations = 0
!    
!    do while (.true.)
!    
!      ! The Newton iteration reads
!      !
!      !    u_n+1  =  u_n  -  [J''(u_n)]^-1  J'(u_n)
!      !
!      ! or in other words,
!      !
!      !    u_n+1  =  u_n  +  [J''(u_n)]^-1  d_n
!      !
!      ! with the residual   d_n = -J'(u_n)   specifying a 'search direction'.
!      
!      ! -------------------------------------------------------------
!      ! Get the current residual / search direction
!      ! -------------------------------------------------------------
!
!      ! Compute the basic (unpreconditioned) search direction d_n.
!      call newtonit_getResidual (p_rkktsystem,p_rdescentDir,&
!          rsolver%rnewtonParams%dresFinal)
!
!      if (rsolver%rnewtonParams%nnonlinearIterations .eq. 1) then
!        ! Remember the initial residual
!        rsolver%rnewtonParams%dresInit = rsolver%rnewtonParams%dresFinal
!      end if
!
!      ! -------------------------------------------------------------
!      ! Check for convergence
!      ! -------------------------------------------------------------
!      if (newtonit_checkConvergence (rsolver%rnewtonParams)) exit
!      
!      ! -------------------------------------------------------------
!      ! Check for divergence
!      ! -------------------------------------------------------------
!      if (newtonit_checkDivergence (rsolver%rnewtonParams)) exit
!      
!      ! -------------------------------------------------------------
!      ! Check other stopping criteria
!      ! -------------------------------------------------------------
!      if (newtonit_checkIterationStop (rsolver%rnewtonParams)) exit
!      
!      ! -------------------------------------------------------------
!      ! Preconditioning with the Newton matrix
!      ! -------------------------------------------------------------
!      
!      ! Actual Newton iteration. Apply the Newton preconditioner
!      ! to get the Newton search direction:
!      !
!      !    J''(u_n) g_n  =  d_n
!      !
!      call newtonlin_precondNewton (rsolver%rlinsolParam,&
!          rkktsystemDirDeriv,p_rdescentDir)
!      
!      ! -------------------------------------------------------------
!      ! Update of the solution
!      ! -------------------------------------------------------------
!
!      ! Update the control according to the search direction:
!      !
!      !    u_n+1  =  u_n  +  g_n
!      !
!      ! or to any configured step-length control rule.
!      call newtonit_updateControl (rsolver,p_rkktsystem,p_rdescentDir)
!      
!      ! -------------------------------------------------------------
!      ! Proceed with the next iteration
!      ! -------------------------------------------------------------
!      ! Next iteration
!      rsolver%rnewtonParams%nnonlinearIterations = &
!          rsolver%rnewtonParams%nnonlinearIterations + 1
!    
!    end do
    
  end subroutine


end module
