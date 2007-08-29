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

MODULE heatcond_solver

  USE fsystem
  USE storage
  USE linearsolver
  USE boundary
  USE bilinearformevaluation
  USE linearformevaluation
  USE cubature
  USE matrixfilters
  USE vectorfilters
  USE bcassembly
  USE triangulation
  USE spatialdiscretisation
  USE sortstrategy
  USE coarsegridcorrection
  USE ucd
  USE timestepping
  USE genoutput
  
  USE collection
  USE paramlist
    
  USE heatcond_callback
  
  USE heatcond_basic
  
  IMPLICIT NONE
  
CONTAINS

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE hc5_initSolver (rproblem)
  
!<description>
  ! Initialises the linear solver according to the problem rproblem.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
    INTEGER :: ilvmin,ilvmax
    INTEGER :: i

    ! Error indicator during initialisation of the solver
    INTEGER :: ierror    
  
    ! A filter chain to filter the vectors and the matrix during the
    ! solution process.
    TYPE(t_filterChain), DIMENSION(:), POINTER :: p_RfilterChain

    ! A solver node that accepts parameters for the linear solver    
    TYPE(t_linsolNode), POINTER :: p_rsolverNode,p_rsmoother
    TYPE(t_linsolNode), POINTER :: p_rcoarseGridSolver,p_rpreconditioner

    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    TYPE(t_matrixBlock), DIMENSION(NNLEV) :: Rmatrices
    
    ! One level of multigrid
    TYPE(t_linsolMGLevelInfo), POINTER :: p_rlevelInfo
    
    ilvmin = rproblem%ilvmin
    ilvmax = rproblem%ilvmax
    
    ! During the linear solver, the boundary conditions must
    ! frequently be imposed to the vectors. This is done using
    ! a filter chain. As the linear solver does not work with 
    ! the actual solution vectors but with defect vectors instead,
    ! a filter for implementing the real boundary conditions 
    ! would be wrong.
    ! Therefore, create a filter chain with one filter only,
    ! which implements Dirichlet-conditions into a defect vector.
    rproblem%RfilterChain(1)%ifilterType = FILTER_DISCBCDEFREAL

    ! Now we have to build up the level information for multigrid.
    !
    ! At first, initialise a standard interlevel projection structure. We
    ! can use the same structure for all levels.
    CALL mlprj_initProjectionDiscr (rproblem%rprojection,&
         rproblem%RlevelInfo(ilvmax)%p_rdiscretisation)
    
    ! Create a Multigrid-solver. Attach the above filter chain
    ! to the solver, so that the solver automatically filters
    ! the vector during the solution process.
    p_RfilterChain => rproblem%RfilterChain
    CALL linsol_initMultigrid (p_rsolverNode,p_RfilterChain)
    
    ! Then set up smoothers / coarse grid solver:
    DO i=ilvmin,ilvmax
      
      ! On the coarsest grid, set up a coarse grid solver and no smoother
      ! On finer grids, set up a smoother but no coarse grid solver.
      NULLIFY(p_rpreconditioner)
      NULLIFY(p_rsmoother)
      NULLIFY(p_rcoarseGridSolver)
      IF (i .EQ. ilvmin) THEN
        ! Set up a BiCGStab solver with ILU preconditioning as coarse grid solver
        ! would be:
        ! CALL linsol_initMILUs1x1 (p_rpreconditioner,0,0.0_DP)
        ! CALL linsol_initBiCGStab (p_rcoarseGridSolver,p_rpreconditioner,p_RfilterChain)
        
        ! Set up UMFPACK coarse grid solver.
        CALL linsol_initUMFPACK4 (p_rcoarseGridSolver)

      ELSE
        ! Setting up Jacobi smoother for multigrid would be:
        ! CALL linsol_initJacobi (p_rsmoother)

        ! Setting up Jin-Wei-Tam smoother for multigrid would be:
        ! CALL linsol_initJinWeiTam (p_rsmoother)

        ! Set up an ILU smoother for multigrid with damping parameter 0.7,
        ! 4 smoothing steps:
         CALL linsol_initMILUs1x1 (p_rsmoother,0,0.0_DP)
         CALL linsol_convertToSmoother (p_rsmoother,4,0.7_DP)
        
      END IF
    
      ! Add the level.
      CALL linsol_addMultigridLevel (p_rlevelInfo,p_rsolverNode, rproblem%rprojection,&
                                     p_rsmoother,p_rsmoother,p_rcoarseGridSolver)
    END DO
    
    ! Set the output level of the solver to 2 for some output
    p_rsolverNode%ioutputLevel = 2

    ! Attach the system matrices to the solver.
    !
    ! We copy our matrices to a big matrix array and transfer that
    ! to the setMatrices routines. This intitialises then the matrices
    ! on all levels according to that array.
    Rmatrices(ilvmin:ilvmax) = rproblem%RlevelInfo(ilvmin:ilvmax)%rmatrix
    CALL linsol_setMatrices(p_RsolverNode,Rmatrices(ilvmin:ilvmax))
    
    ! Save the solver node in the problem structure, finish
    rproblem%p_rsolverNode => p_rsolverNode

    ! Allocate memory, initialise solver structures according to the
    ! linear system we just attached.
    CALL linsol_initStructure (p_rsolverNode,ierror)
    IF (ierror .NE. LINSOL_ERR_NOERROR) STOP
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE hc5_doneSolver (rproblem)
  
!<description>
  ! Releases the solver from the problem structure.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

!</subroutine>
 
    ! Release solver data and structure
    CALL linsol_doneData (rproblem%p_rsolverNode)
    CALL linsol_doneStructure (rproblem%p_rsolverNode)
    
    ! Release the solver node and all subnodes attached to it (if at all):
    CALL linsol_releaseSolver (rproblem%p_rsolverNode)
    
    ! Release the multilevel projection structure.
    CALL mlprj_doneProjection (rproblem%rprojection)

  END SUBROUTINE

END MODULE