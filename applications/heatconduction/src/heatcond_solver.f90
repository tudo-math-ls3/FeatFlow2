!##############################################################################
!# ****************************************************************************
!# <name> heatcond_solver </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains routines to maintain the linear solver used in the
!# nonstationary heat conduction problem. The following routines can be found
!# here:
!#
!# 1.) hc5_initSolver
!#     -> Initialise the linear solver, allocate memory.
!#
!# 2.) hc5_doneSolver
!#     -> Clean up the linear solver, release memory.
!# </purpose>
!##############################################################################

module heatcond_solver

  use fsystem
  use storage
  use linearsolver
  use boundary
  use bilinearformevaluation
  use linearformevaluation
  use cubature
  use matrixfilters
  use vectorfilters
  use bcassembly
  use triangulation
  use spatialdiscretisation
  use coarsegridcorrection
  use ucd
  use timestepping
  use genoutput
  
  use collection
  use paramlist
    
  use heatcond_callback
  
  use heatcond_basic
  
  implicit none
  
contains

  ! ***************************************************************************

!<subroutine>

  subroutine hc5_initSolver (rproblem)
  
!<description>
  ! Initialises the linear solver according to the problem rproblem.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
    integer :: ilvmin,ilvmax
    integer :: i

    ! Error indicator during initialisation of the solver
    integer :: ierror
  
    ! A solver node that accepts parameters for the linear solver
    type(t_linsolNode), pointer :: p_rsolverNode,p_rsmoother
    type(t_linsolNode), pointer :: p_rcoarseGridSolver,p_rpreconditioner

    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    type(t_matrixBlock), dimension(1:rproblem%ilvmax) :: Rmatrices
    
    ! One level of multigrid
    type(t_linsolMG2LevelInfo), pointer :: p_rlevelInfo
    
    ilvmin = rproblem%ilvmin
    ilvmax = rproblem%ilvmax
    
    ! Now we have to build up the level information for multigrid.
    !
    ! At first, initialise a standard interlevel projection structure. We
    ! can use the same structure for all levels.
    call mlprj_initProjectionDiscr (rproblem%rprojection,&
         rproblem%RlevelInfo(ilvmax)%p_rdiscretisation)
    
    ! Create a Multigrid-solver. Attach the above filter chain
    ! to the solver, so that the solver automatically filters
    ! the vector during the solution process.
    call linsol_initMultigrid2 (p_rsolverNode,ilvmax-ilvmin+1)
    
    ! Then set up smoothers / coarse grid solver:
    do i=ilvmin,ilvmax
      
      ! And add this multi-grid level. We will use the same smoother
      ! for pre- and post-smoothing.
      call linsol_getMultigrid2Level (p_rsolverNode,i-ilvmin+1,p_rlevelInfo)

      ! Attach the filter chain which imposes boundary conditions on that level.
      p_rlevelInfo%p_RfilterChain => rproblem%RlevelInfo(i)%RfilterChain

      ! During the linear solver, the boundary conditions must
      ! frequently be imposed to the vectors. This is done using
      ! a filter chain. As the linear solver does not work with
      ! the actual solution vectors but with defect vectors instead,
      ! a filter for implementing the real boundary conditions
      ! would be wrong.
      ! Therefore, create a filter chain with one filter only,
      ! which implements Dirichlet-conditions into a defect vector.
      call filter_clearFilterChain (rproblem%RlevelInfo(i)%RfilterChain,&
          rproblem%RlevelInfo(i)%nfilters)
      call filter_newFilterDiscBCDef (rproblem%RlevelInfo(i)%RfilterChain,&
          rproblem%RlevelInfo(i)%nfilters,rproblem%RlevelInfo(i)%p_rdiscreteBC)

      ! On the coarsest grid, set up a coarse grid solver and no smoother
      ! On finer grids, set up a smoother but no coarse grid solver.
      nullify(p_rpreconditioner)
      nullify(p_rsmoother)
      nullify(p_rcoarseGridSolver)
      if (i .eq. ilvmin) then
        ! Set up a BiCGStab solver with ILU preconditioning as coarse grid solver
        ! would be:
        ! call linsol_initMILUs1x1 (p_rpreconditioner,0,0.0_DP)
        ! call linsol_initBiCGStab (p_rcoarseGridSolver,p_rpreconditioner,&
        !     p_rlevelInfo%p_RfilterChain)
        
        ! Set up UMFPACK coarse grid solver.
        call linsol_initUMFPACK4 (p_rcoarseGridSolver)
        
        ! Initialise the level with the coarse grid solver.
        p_rlevelInfo%p_rcoarseGridSolver => p_rcoarseGridSolver

      else
        ! Setting up Jacobi smoother for multigrid would be:
        ! call linsol_initJacobi (p_rsmoother)

        ! Set up an ILU smoother for multigrid with damping parameter 0.7,
        ! 4 smoothing steps:
        call linsol_initMILUs1x1 (p_rsmoother,0,0.0_DP)
        call linsol_convertToSmoother (p_rsmoother,4,0.7_DP)
        
        ! Use the smoother for pre- and postsmoothing
        p_rlevelInfo%p_rpresmoother => p_rsmoother
        p_rlevelInfo%p_rpostsmoother => p_rsmoother

      end if
    
    end do
    
    ! Set the output level of the solver to 2 for some output
    p_rsolverNode%ioutputLevel = 2

    ! Attach the system matrices to the solver.
    !
    ! We copy our matrices to a big matrix array and transfer that
    ! to the setMatrices routines. This intitialises then the matrices
    ! on all levels according to that array.
    Rmatrices(ilvmin:ilvmax) = rproblem%RlevelInfo(ilvmin:ilvmax)%rmatrix
    call linsol_setMatrices(p_RsolverNode,Rmatrices(ilvmin:ilvmax))
    
    ! Save the solver node in the problem structure, finish
    rproblem%p_rsolverNode => p_rsolverNode

    ! Allocate memory, initialise solver structures according to the
    ! linear system we just attached.
    call linsol_initStructure (p_rsolverNode,ierror)

    if (ierror .ne. LINSOL_ERR_NOERROR) then
      call output_line("Matrix structure invalid!",OU_CLASS_ERROR)
      call sys_halt()
    end if
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine hc5_doneSolver (rproblem)
  
!<description>
  ! Releases the solver from the problem structure.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!</subroutine>
 
    ! Release solver data and structure
    call linsol_doneData (rproblem%p_rsolverNode)
    call linsol_doneStructure (rproblem%p_rsolverNode)
    
    ! Release the solver node and all subnodes attached to it (if at all):
    call linsol_releaseSolver (rproblem%p_rsolverNode)
    
    ! Release the multilevel projection structure.
    call mlprj_doneProjection (rproblem%rprojection)

  end subroutine

end module
