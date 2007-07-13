!##############################################################################
!# ****************************************************************************
!# <name> cc2dmediumm2timesupersystem </name>
!# ****************************************************************************
!#
!# <purpose>
!# </purpose>
!##############################################################################

MODULE cc2dmediumm2spacetimesolver

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
  USE coarsegridcorrection
  USE spdiscprojection
  USE nonlinearsolver
  USE paramlist
  USE linearsolverautoinitialise
  USE matrixrestriction
  USE paramlist
  USE timestepping
  USE l2projection
  
  USE collection
  USE convection
    
  USE cc2dmediumm2basic
  USE cc2dmedium_callback

  USE cc2dmediumm2nonlinearcore
  USE cc2dmediumm2nonlinearcoreinit
  USE cc2dmediumm2stationary
  USE adaptivetimestep
  USE cc2dmediumm2timeanalysis
  USE cc2dmediumm2boundary
  USE cc2dmediumm2discretisation
  USE cc2dmediumm2postprocessing
  USE cc2dmediumm2matvecassembly
  USE cc2dmediumm2spacetimediscret
  USE cc2dmediumm2nonlinearcore
  
  USE spacetimevectors
  USE spacetimeinterlevelprojection
  USE dofmapping
  
  USE matrixio
    
  IMPLICIT NONE

!<constantblock description="Algorithm identifiers">

  ! Undefined algorithm
  INTEGER, PARAMETER :: SPTILS_ALG_UNDEFINED     = 0
  
  ! Preconditioned defect correction (Richardson iteration);
  ! $x_{n+1} = x_n + \omega P^{-1} (b-Ax)$
  INTEGER, PARAMETER :: SPTILS_ALG_DEFCORR       = 1
  
  ! Block Jacobi preconditioner which applies a spatial preconditioner
  ! to each diagonal block of the coupled space-time system.
  INTEGER, PARAMETER :: SPTILS_ALG_BLOCKJACOBI   = 2
  
  ! Multigrid iteration
  INTEGER, PARAMETER :: SPTILS_ALG_MULTIGRID     = 9

!</constantblock>

!<types>

!<typeblock>
  
  ! The following structure defines a 'solver node' for the linear
  ! solver / preconditioner of the global system.
  ! As solver, currently a single-grid preconditioner,
  ! defect correction solver and multigrid solver/preconditioner is 
  ! supported. The multigrid preconditioner solves the system
  ! simultaneously in space/time!
  !
  ! Like in the standard linear solver library, the linear solver
  ! nodes for the global space/time system can be 'hung' into each
  ! other to configure the solver.
  
  TYPE t_sptilsNode
    
    ! OUTPUT: Result
    ! The result of the solution process.
    ! =0: success. =1: iteration broke down, diverging, =2: error in the parameters,
    ! <0: algorithm-specific error
    INTEGER                    :: iresult
    
    ! OUTPUT: Number of performed iterations, if the solver
    ! is of iterative nature.
    ! Is to 1 by the solver if not used (indicating at least 1 performed 
    ! iteration, which is always the case).
    INTEGER                    :: iiterations
    
    ! OUTPUT PARAMETER FOR SOLVERS WITH RESIDUAL CHECK: 
    ! Norm of initial residuum
    REAL(DP)                        :: dinitialDefect

    ! OUTPUT PARAMETER FOR SOLVERS WITH RESIDUAL CHECK: 
    ! Norm of final residuum
    REAL(DP)                        :: dfinalDefect

    ! OUTPUT PARAMETER FOR ITERATIVE SOLVERS WITH RESIDUAL CHECK: 
    ! Convergence rate
    REAL(DP)                        :: dconvergenceRate

    ! OUTPUT PARAMETER FOR ITERATIVE SOLVERS WITH RESIDUAL CHECK: 
    ! Asymptotic convergence rate
    REAL(DP)                        :: dasymptoticConvergenceRate

    ! OUTPUT PARAMETER:
    ! Total time for solver
    REAL(DP)                        :: dtimeTotal

    ! OUTPUT PARAMETER FOR SOLVERS THAT SUPPORT FILTERING:
    ! Total time for filtering
    REAL(DP)                        :: dtimeFiltering
    
    ! INPUT PARAMETER:
    ! General solver parameter; solver specific use.
    ! Standard value = 1.0 (corresponds to 'no damping' e.g. with the defect
    ! correction iteration)
    REAL(DP)                        :: domega  = 1.0_DP

    ! INPUT PARAMETER FOR ITERATIVE SOLVERS: 
    ! Relative stopping criterion. Stop iteration if
    ! !!defect!! < EPSREL * !!initial defect!!.
    ! =0: ignore, use absolute stopping criterion; standard = 1E-5
    ! Remark: don't set depsAbs=depsRel=0!
    REAL(DP)                        :: depsRel = 1E-5_DP

    ! INPUT PARAMETER FOR ITERATIVE SOLVERS: 
    ! Absolute stopping criterion. Stop iteration if
    ! !!defect!! < EPSREL.
    ! =0: ignore, use relative stopping criterion; standard = 1E-5
    ! Remark: don't set depsAbs=depsRel=0!
    REAL(DP)                        :: depsAbs = 1E-5_DP

    ! INPUT PARAMETER FOR ITERATIVE SOLVERS: 
    ! Relative divergence criterion.  Treat iteration as
    ! diverged if
    !   !!defect!! >= DIVREL * !!initial defect!!
    ! A value of SYS_INFINITY disables the relative divergence check.
    ! standard = 1E3
    REAL(DP)                        :: ddivRel = 1E3_DP

    ! INPUT PARAMETER FOR ITERATIVE SOLVERS: 
    ! Absolute divergence criterion.  Treat iteration as
    ! diverged if
    !   !!defect!! >= DIVREL
    ! A value of SYS_INFINITY disables the absolute divergence check.
    ! standard = SYS_INFINITY
    REAL(DP)                        :: ddivAbs = SYS_INFINITY

    ! INPUT PARAMETER FOR ITERATIVE SOLVERS: 
    ! RHS-vector is treated as zero if max(defect) < drhsZero
    REAL(DP)                        :: drhsZero = 1E-90_DP

    ! INPUT PARAMETER FOR ITERATIVE SOLVERS: 
    ! Type of stopping criterion to use. One of the
    ! LINSOL_STOP_xxxx constants.
    INTEGER                    :: istoppingCriterion = LINSOL_STOP_STANDARD

    ! INPUT PARAMETER FOR ITERATIVE SOLVERS: 
    ! Minimum number of iterations top perform
    INTEGER                    :: nminIterations = 1

    ! INPUT PARAMETER FOR ITERATIVE SOLVERS: 
    ! Maximum number of iterations top perform
    INTEGER                    :: nmaxIterations = 50
    
    ! INPUT PARAMETER FOR SOLVERS WITH RESIDUAL CHECK: 
    ! Perform residual checks
    ! YES: check residuals (default), 
    ! NO: don't check residuals, simply perform as many iterations as 
    ! configured by nminIterations. Is typically set to NO for smoothing
    ! with a solver.
    INTEGER                    :: iresCheck = YES

    ! INPUT PARAMETER FOR SOLVERS WITH RESIDUAL CHECK: 
    ! Type of norm to use in the residual checking (cf. linearalgebra.f90).
    ! =0: euclidian norm, =1: l1-norm, =2: l2-norm, =3: MAX-norm
    INTEGER                    :: iresNorm = 2
    
    ! INPUT PARAMETER: Output level
    ! This determines the output level of the solver.
    ! =0: no output, =1: basic output, =2, extended output
    INTEGER                    :: ioutputLevel = 2

    ! INPUT PARAMETER FOR ITERATIVE SOLVERS WITH RESIDUAL CHECK:
    ! Number of iterations to perform before printing out the
    ! norm of the residual to screen.
    ! =1: Print residual in every iteration
    INTEGER                    :: niteResOutput = 1

    ! READ ONLY: Algorithm identifier. 
    ! One of the SPTILS_ALG_xxxx constants.
    ! Depending on the value, a solver-specific structure may be 
    ! assigned to this structure below.
    INTEGER                    :: calgorithm = SPTILS_ALG_UNDEFINED
    
    ! A pointer to the problem structure that defines the problem.
    TYPE(t_problem), POINTER   :: p_rproblem
    
    ! STATUS FOR ITERATIVE SOLVERS: Current iteration
    INTEGER                    :: icurrentIteration

    ! Space time discretisation structure that allows apply the system matrix.
    TYPE(t_ccoptSpaceTimeDiscretisation) :: rspaceTimeDiscr
    
    ! Pointer to a structure for the Defect correction solver; NULL() if not set
    TYPE (t_sptilsSubnodeDefCorr), POINTER        :: p_rsubnodeDefCorr     => NULL()

    ! Pointer to a structure for the block Jacobi preconditioner
    TYPE(t_sptilsSubnodeBlockJacobi), POINTER     :: p_rsubnodeBlockJacobi => NULL()

    ! Pointer to a structure for the multigrid preconditioner
    TYPE(t_sptilsSubnodeMultigrid), POINTER       :: p_rsubnodeMultigrid => NULL()

  END TYPE
  
!</typeblock>

! *****************************************************************************

!<typeblock>
  
  ! This structure realises the subnode for the BiCGStab solver.
  ! The entry p_rpreconditioner points either to NULL() or to another
  ! t_sptilsNode structure for the solver that realises the 
  ! preconditioning.
  
  TYPE t_sptilsSubnodeDefCorr
  
    ! A pointer to the solver node for the preconditioner or NULL(),
    ! if no preconditioner is used.
    TYPE(t_sptilsNode), POINTER       :: p_rpreconditioner            => NULL()
    
    ! A temporary space-time vector.
    TYPE(t_spacetimeVector)           :: rtempVector
    
    ! A temporary space-time vector.
    TYPE(t_spacetimeVector)           :: rtempVector2
    
  END TYPE
  
!</typeblock>

! *****************************************************************************

!<typeblock>
  
  ! This structure realises the subnode for the Block Jacobi solver.
  
  TYPE t_sptilsSubnodeBlockJacobi

    ! Pointer to a spatial preconditioner structure that defines the
    ! preconditioning in each substep of the global system.
    TYPE(t_ccspatialPreconditioner), POINTER :: p_rspatialPreconditioner
    
  END TYPE
  
!</typeblock>

! *****************************************************************************

!<typeblock>
  
  ! Level data for multigrid. This structure forms an entry in a linked list
  ! of levels which multigrid uses for solving the given problem.
  ! The solver subnode t_sptilsSubnodeMultigrid holds pointers to the head
  ! and the tail of the list.
  
  TYPE t_sptilsMGLevelInfo
  
    ! Space time discretisation structure that allows apply the system matrix.
    TYPE(t_ccoptSpaceTimeDiscretisation) :: rspaceTimeDiscr
    
    ! A RHS vector for that level. The memory for this is allocated
    ! in initStructure and released in doneStructure.
    TYPE(t_spacetimeVector)             :: rrhsVector

    ! A solution vector for that level. The memory for this is allocated
    ! in initStructure and released in doneStructure.
    TYPE(t_spacetimeVector)             :: rsolutionVector

    ! A temporary vector for that level. The memory for this is allocated
    ! in initStructure and released in doneStructure.
    TYPE(t_spacetimeVector)             :: rtempVector

    ! A pointer to the solver node for the presmoothing or NULL(),
    ! if no presmoother is used.
    TYPE(t_sptilsNode), POINTER         :: p_rpresmoother         => NULL()

    ! A pointer to the solver node for the postsmoothing or NULL(),
    ! if no presmoother is used.
    TYPE(t_sptilsNode), POINTER         :: p_rpostsmoother        => NULL()

    ! A pointer to the solver node for the coarse grid solver if
    ! the level corresponding to this structure represents the coarse
    ! grid. NULL() otherwise.
    TYPE(t_sptilsNode), POINTER         :: p_rcoarseGridSolver    => NULL()
    
    ! An interlevel projection structure that configures the transport
    ! of defect/solution vectors from one level to another.
    ! For the coarse grid, this structure is ignored.
    ! For a finer grid (e.g. level 4), this defines the grid transfer
    ! between the current and the lower level (so here between level 3 and
    ! level 4).
    TYPE(t_sptiProjection)   :: rinterlevelProjection
    
    ! STATUS/INTERNAL: A temporary vector used for prolongation/restriction
    TYPE(t_vectorBlock)      :: rprjVector
    
    ! STATUS/INTERNAL: MG cycle information.
    ! Number of cycles to perform on this level.
    INTEGER                        :: ncycles

    ! STATUS/INTERNAL: MG cycle information.
    ! Number of remaining cycles to perform on this level.
    INTEGER                        :: ncyclesRemaining
    
    ! STATUS/INTERNAL: Number of current cycle on that level. 
    ! Only used if adaptive cycles are activated by setting 
    ! depsRelCycle < 1E99_DP in t_sptilsSubnodeMultigrid.
    INTEGER                        :: icycleCount   = 0
    
  END TYPE
  
!</typeblock>

  ! ***************************************************************************

!<typeblock>
  
  ! This structure realises the subnode for the Multigrid solver.
  ! There are pointers to the head and the tail of a linked list of
  ! t_mgLevelInfo structures which hold the data of all levels.
  ! The tail of the list corresponds to the maximum level where
  ! multigrid solves the problem.
  
  TYPE t_sptilsSubnodeMultigrid
  
    ! INPUT PARAMETER: Cycle identifier. 
    !  0=F-cycle, 
    !  1=V-cycle, 
    !  2=W-cycle.
    INTEGER                       :: icycle                   = 0
    
    ! Array of t_sptilsMGLevelInfo structures for all the levels the MG solver
    ! should handle.
    TYPE(t_sptilsMGLevelInfo), DIMENSION(:), POINTER :: p_Rlevels => NULL()
    
    ! INPUT PARAMETER: Minimum level in p_Rlevels.
    INTEGER                       :: NLMIN                    = 1
    
    ! INPUT PARAMETER: Maximum level in p_Rlevels.
    INTEGER                       :: NLMAX                    = 1
    
  END TYPE
  
!</typeblock>

!</types>

CONTAINS

  ! ***************************************************************************
  
!<subroutine>
  
  RECURSIVE SUBROUTINE sptils_setMatrices (rsolverNode,Rmatrices)
  
!<description>
  
  ! Initialises the system matrices in the solver node rsolverNode and
  ! in all nodes attached to it (preconditioners,...). For this purpose,
  ! the corresponding initialisation routine of each solver is called.
  !
  ! The last element in Rmatrices defines the system matrix of the
  ! linear system that is to solve. If there's a solver involved
  ! that needs matrices on multiple grids (e.g. multigrid), Rmatrices
  ! must contain the matrices of all levels. In this case, Rmatrices(1)
  ! defines the matrix on the coarsest level to solve, Rmatrices(2)
  ! the matrix of the first refinement etc. up to the level where
  ! to solve the problem.
  !
  ! The matrix structures passed in Rmatrices are copied to the system
  ! matrix structures in the solver node. If the system matrix changes
  ! (or more precisely, if handles or structural data of the system matrix
  ! changes), sptils_setMatrices has to be called again to update all
  ! the system matrix/matrices in the solver!
  
!</description>
  
!<input>
  ! Array of space-time discretisation structures that specify the system matrices 
  ! on all levels of the discretisation.
  ! This is passed through all initialisation routines, but actually used 
  ! only by the multigrid initialisation routine.
  TYPE(t_ccoptSpaceTimeDiscretisation), DIMENSION(:), INTENT(IN) :: Rmatrices
!</input>
  
!<inputoutput>
  ! The solver node which should be initialised
  TYPE(t_sptilsNode), INTENT(INOUT)             :: rsolverNode
!</inputoutput>
  
!</subroutine>

  ! Copy the matrix structure on the finest level to the rsolverNode
  ! structure. This corresponds to the system we want to solve.
  rsolverNode%rspaceTimeDiscr = Rmatrices(UBOUND(Rmatrices,1))
  
  ! Depending on the solver type, call the corresponding initialisation
  ! routine. For all single-grid solvers, pass the matrix on the highest
  ! level in RlinsysInfo. For multigrid we pass all matrices.
  ! That way the solvers can decide to do anything else with the matrix
  ! or to initialise their subsolvers.
  !
  ! For single grid solvers, the whole Rmatrices array is usually also
  ! passed to the initialisation routine. This allows to initialise
  ! the sub-nodes of each solver.
  
  SELECT CASE(rsolverNode%calgorithm)
  CASE (SPTILS_ALG_DEFCORR)
    CALL sptils_setMatrixDefCorr (rsolverNode,Rmatrices)
  CASE (SPTILS_ALG_BLOCKJACOBI)
    CALL sptils_setMatrixBlockJacobi (rsolverNode,Rmatrices)
  CASE (SPTILS_ALG_MULTIGRID)
    CALL sptils_setMatrixMultigrid (rsolverNode,Rmatrices)
  CASE DEFAULT
  END SELECT

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>
  
  RECURSIVE SUBROUTINE sptils_initStructure (rsolverNode, ierror)
  
!<description>
  
  ! Initialises the problem structure in the solver node rsolverNode by calling
  ! the initialisation routine of the appropriate solver. The solver
  ! initialisation routine itself can call this procedure to initialise
  ! its sub-solver nodes.
  
!</description>
  
!<inputoutput>
  ! The solver node which should be initialised
  TYPE(t_sptilsNode), INTENT(INOUT)                     :: rsolverNode
!</inputoutput>

!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to 
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  INTEGER, INTENT(OUT) :: ierror
!</output>
  
!</subroutine>

    ! A-priori we have no error...
    ierror = LINSOL_ERR_NOERROR
    
    ! Call the structure-init routine of the specific solver
    
    SELECT CASE(rsolverNode%calgorithm)
    CASE (SPTILS_ALG_DEFCORR)
      CALL sptils_initStructureDefCorr (rsolverNode,ierror)
    CASE (SPTILS_ALG_BLOCKJACOBI)
      CALL sptils_initStructureBlockJacobi (rsolverNode,ierror)
    CASE (SPTILS_ALG_MULTIGRID)
      CALL sptils_initStructureMultigrid (rsolverNode,ierror)
    END SELECT
  
  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>
  
  RECURSIVE SUBROUTINE sptils_initData (rsolverNode, ierror)
  
!<description>
  ! Initialises the problem data in the solver node rsolverNode by calling
  ! the initialisation routine of the appropriate solver. The solver
  ! initialisation routine itself can call this procedure to initialise
  ! its sub-solver nodes.
  ! The initialisation of the problem structure allowes the solver component
  ! to perform some 'precalculation', e.g. the UMFPACK4 or ILU solver can 
  ! perform a numerical factorisation. The problem structure usually does
  ! not change during a simulation, except when the grid moves e.g.
!</description>
  
!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to 
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  INTEGER, INTENT(OUT) :: ierror
!</output>

!<inputoutput>
  ! The solver node which should be initialised
  TYPE(t_sptilsNode), INTENT(INOUT)                     :: rsolverNode
!</inputoutput>

!</subroutine>

    ! A-priori we have no error...
    ierror = LINSOL_ERR_NOERROR

    ! Call the data-init routine of the specific solver
    
    SELECT CASE(rsolverNode%calgorithm)
    CASE (LINSOL_ALG_DEFCORR)
      CALL sptils_initDataDefCorr (rsolverNode,ierror)
    CASE (SPTILS_ALG_BLOCKJACOBI)
      CALL sptils_initDataBlockJacobi (rsolverNode,ierror)
    CASE (SPTILS_ALG_MULTIGRID)
      CALL sptils_initDataMultigrid (rsolverNode,ierror)
    END SELECT
  
  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>
  
  RECURSIVE SUBROUTINE sptils_updateStructure (rsolverNode, ierror)
  
!<description>
  
  ! Reinitialises the problem structure in the solver node rsolverNode
  
!</description>
  
!<inputoutput>
  ! The solver node which should be reinitialised
  TYPE(t_sptilsNode), INTENT(INOUT)                     :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to 
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  INTEGER, INTENT(OUT) :: ierror
!</output>
  
!</subroutine>

    ! A-priori we have no error...
    ierror = LINSOL_ERR_NOERROR

    ! For now, call the structure-done- and structure-init routines.
    ! Maybe rewritten later for higher performance or special cases...
    
    CALL sptils_doneStructure (rsolverNode)
    CALL sptils_initStructure (rsolverNode,ierror)
  
  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>
  
  RECURSIVE SUBROUTINE sptils_updateData (rsolverNode,ierror)
  
!<description>
  ! Reinitialises the problem data in the solver node rsolverNode.
!</description>
  
!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to 
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  INTEGER, INTENT(OUT) :: ierror
!</output>

!<inputoutput>
  ! The solver node containing the solver confuguration
  TYPE(t_sptilsNode), INTENT(INOUT)                     :: rsolverNode
!</inputoutput>
  
!</subroutine>

    ! A-priori we have no error...
    ierror = LINSOL_ERR_NOERROR

    ! For now, call the data-done- and data-init routines.
    ! Maybe rewritten later for higher performance or special cases...
    
    CALL sptils_doneData (rsolverNode)
    CALL sptils_initData (rsolverNode,ierror)

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>
  
  RECURSIVE SUBROUTINE sptils_doneData (rsolverNode)
  
!<description>
  ! Releases the problem data in the solver node rsolverNode.
!</description>
  
!<inputoutput>
  
  ! The solver node containing the solver confuguration
  TYPE(t_sptilsNode), INTENT(INOUT)                     :: rsolverNode
  
!</inputoutput>
  
!</subroutine>

    ! Call the data-done routine of the specific solver
    
    SELECT CASE(rsolverNode%calgorithm)
    CASE (LINSOL_ALG_DEFCORR)
      CALL sptils_doneDataDefCorr (rsolverNode)
    CASE (SPTILS_ALG_BLOCKJACOBI)
      CALL sptils_doneDataBlockJacobi (rsolverNode)
    CASE (SPTILS_ALG_MULTIGRID)
      CALL sptils_doneDataMultigrid (rsolverNode)
    END SELECT

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>
  
  RECURSIVE SUBROUTINE sptils_doneStructure (rsolverNode)
  
!<description>
  ! Releases the problem structure in the solver node rsolverNode.
  ! This is done by calling the appropriate routine of the
  ! actual solver.
!</description>
  
!<inputoutput>
  ! The solver node which should be reinitialised
  TYPE(t_sptilsNode), INTENT(INOUT)                     :: rsolverNode
!</inputoutput>
  
!</subroutine>

    ! Call the data-done routine of the specific solver
    
    SELECT CASE(rsolverNode%calgorithm)
    CASE (LINSOL_ALG_DEFCORR)
      CALL sptils_doneStructureDefCorr (rsolverNode)
    CASE (SPTILS_ALG_BLOCKJACOBI)
      CALL sptils_doneStructureBlockJacobi (rsolverNode)
    CASE (SPTILS_ALG_MULTIGRID)
      CALL sptils_doneStructureMultigrid (rsolverNode)
    END SELECT

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>
  
  RECURSIVE SUBROUTINE sptils_releaseSolver (p_rsolverNode,bkeepSolverNode)
  
!<description>
  ! This routine releases the solver node rsolverNode and all attached
  ! sub-solver nodes (for preconditioning, smoothing etc. if specified)
  ! from the heap.
!</description>
  
!<input>
  ! OPTIONAL: If set to TRUE, the structure p_rsolverNode is not released
  ! from memory. If set to FALSE or not existent (the usual setting), the 
  ! structure p_rsolverNode will also be removed from the heap after 
  ! cleaning up.
  ! Remark: The subnodes of the solver (if there are any) are always
  ! removed from the heap!
  LOGICAL, INTENT(IN), OPTIONAL :: bkeepSolverNode
!</input>

!<inputoutput>
  ! The solver node which is to be released. If the node contains attached
  ! subsolvers (preconditioners, smoothers,...) they are also released.
  ! On return, rsolverNode is NULL().
  TYPE(t_sptilsNode), POINTER     :: p_rsolverNode
!</inputoutput>
  
!</subroutine>

    IF (.NOT. ASSOCIATED(p_rsolverNode)) THEN
      ! Warning message, return
      PRINT *,'sptils_releaseSolver warning: Solver note not assigned!'
      RETURN
    END IF

    ! Depending on the solver type, call the corresponding done-routine
    SELECT CASE(p_rsolverNode%calgorithm)
    CASE (LINSOL_ALG_DEFCORR)
      CALL sptils_doneDefCorr (p_rsolverNode)
    CASE (SPTILS_ALG_BLOCKJACOBI)
      CALL sptils_doneBlockJacobi (p_rsolverNode)
    CASE (SPTILS_ALG_MULTIGRID)
      CALL sptils_doneMultigrid (p_rsolverNode)
    CASE DEFAULT
    END SELECT
    
    ! Clean up the associated matrix structure.
    ! Of course, the memory of the matrix is not released from memory, because
    ! if there's a matrix attached, it belongs to the application, not to the
    ! solver!
    !CALL lsysbl_releaseMatrix(p_rsolverNode%rsystemMatrix)
    
    ! Finally release the node itself (if we are allowed to).
    ! Deallocate the structure (if we are allowed to), finish.
    IF (.NOT. PRESENT(bkeepSolverNode)) THEN
      DEALLOCATE(p_rsolverNode)
    ELSE
      IF (.NOT. bkeepSolverNode) THEN
        DEALLOCATE(p_rsolverNode)
      END IF
    END IF
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE sptils_precondDefect (rsolverNode,rd)
  
!<description>
  ! This routine applies the solver configured in rsolverNode
  ! as preconditioner to the defect vector rd.
  ! rd will be overwritten by the preconditioned defect.
  !
  ! The matrix must have been attached to the system before calling
  ! this routine, and the initStructure/initData routines
  ! must have been called to prepare the solver for solving
  ! the problem.
  !
  ! The vector rd must be compatible to the system matrix (same number
  ! of equations, same sorting strategy,...).
!</description>
  
!<inputoutput>
  ! The solver node containing the solver configuration
  TYPE(t_sptilsNode), INTENT(INOUT)                :: rsolverNode
  
  ! On call to this routine: The defect vector to be preconditioned.
  ! Will be overwritten by the preconditioned defect.
  TYPE(t_spacetimeVector), INTENT(INOUT)               :: rd
!</inputoutput>
  
!</subroutine>

    ! Select the solver as configured in rsolverNode and let it perform
    ! the actual preconditioning task.
    
    SELECT CASE(rsolverNode%calgorithm)
    CASE (SPTILS_ALG_DEFCORR)
      CALL sptils_precDefCorr (rsolverNode,rd)
    CASE (SPTILS_ALG_BLOCKJACOBI)
      CALL sptils_precBlockJacobi (rsolverNode,rd)
    CASE (SPTILS_ALG_MULTIGRID)
      CALL sptils_precMultigrid (rsolverNode,rd)
    CASE DEFAULT
      PRINT *,'Unknown space-time preconditioner!'
      CALL sys_halt()
    END SELECT

  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE sptils_initSolverGeneral (rproblem,p_rsolverNode)
  
!<description>
  ! Creates a new solver node p_rsolverNode on the heap with default
  ! settings and returns a pointer to it.
!</description>
  
!<input>
  ! The problem structure that defines the problem.
  ! A pointer to this is saved in the solver node, so the problem structure
  ! must not be released before the solver node is released.
  TYPE(t_problem), TARGET :: rproblem
!</input>
  
!<output>
  ! A pointer to a new solver node on the heap.
  TYPE(t_sptilsNode), POINTER         :: p_rsolverNode
!</output>
  
!</subroutine>

    ! Allocate the node - the default initialisation will set most of the
    ! parameters in the structure.
    ALLOCATE(p_rsolverNode)
    
    ! Save a pointer to the problem structure in the solver node
    p_rsolverNode%p_rproblem => rproblem

  END SUBROUTINE
  
! *****************************************************************************
! Routines for the Defect Correction iteration
! *****************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE sptils_initDefCorr (rproblem,p_rsolverNode,p_rpreconditioner)
  
!<description>
  ! Creates a t_sptilsNode solver structure for the defect correction iteration.
  ! The node can be used to directly solve a problem or to be attached 
  ! as solver or preconditioner to another solver structure. The node can be 
  ! deleted by sptils_releaseSolver.
  !
  ! The defect correction performs nmaxIterations iterations of the type
  !    $$ x_{n+1}  =  x_n  +  (b-Ax_n) $$
  ! with $x_0:=0$. 
  ! It's possible to include a damping parameter to this operation by 
  ! changing rsolverNode%domega to a value $\not =1$. In this case, the
  ! defect correction iteration changes to the Richardson iteration
  !
  !    $$ x_{n+1}  =  x_n  +  \omega(b-Ax_n) $$
  !
  ! By specifying an additional preconditioner, it's possible to perform
  ! the preconditioned defect correction iteration
  !
  !    $$ x_{n+1}  =  x_n  +  \omega P^{-1} (b-Ax_n) $$
!</description>
  
!<input>
  ! The problem structure that defines the problem.
  ! A pointer to this is saved in the solver node, so the problem structure
  ! must not be released before the solver node is released.
  TYPE(t_problem), TARGET :: rproblem

  ! OPTIONAL: A pointer to the solver structure of a solver that should be 
  ! used for preconditioning. If not given or set to NULL(), no preconditioning 
  ! will be used.
  TYPE(t_sptilsNode), POINTER, OPTIONAL   :: p_rpreconditioner
!</input>
  
!<output>
  ! A pointer to a t_sptilsNode structure. Is set by the routine, any previous
  ! value of the pointer is destroyed.
  TYPE(t_sptilsNode), POINTER         :: p_rsolverNode
!</output>
  
!</subroutine>
  
    ! Create a default solver structure
    CALL sptils_initSolverGeneral(rproblem,p_rsolverNode)
    
    ! Initialise the type of the solver
    p_rsolverNode%calgorithm = SPTILS_ALG_DEFCORR
    
    ! Allocate a subnode for our solver.
    ! This initialises most of the variables with default values appropriate
    ! to this solver.
    ALLOCATE(p_rsolverNode%p_rsubnodeDefCorr)
    
    ! Attach the preconditioner if given. 
    
    IF (PRESENT(p_rpreconditioner)) THEN 
      p_rsolverNode%p_rsubnodeDefCorr%p_rpreconditioner => p_rpreconditioner
    END IF

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE sptils_doneDefCorr (rsolverNode)
  
!<description>
  ! This routine releases all temporary memory for the efect correction solver 
  ! from the heap.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure of UMFPACK4 which is to be cleaned up.
  TYPE(t_sptilsNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!</subroutine>
  
    ! Check if there's a preconditioner attached. If yes, release it.
    IF (ASSOCIATED(rsolverNode%p_rsubnodeDefCorr%p_rpreconditioner)) THEN
      CALL sptils_releaseSolver(rsolverNode%p_rsubnodeDefCorr%p_rpreconditioner)
    END IF

    ! Release memory if still associated
    CALL sptils_doneDataDefCorr (rsolverNode)
    CALL sptils_doneStructureDefCorr (rsolverNode)
    
    ! Release the subnode structure
    DEALLOCATE(rsolverNode%p_rsubnodeDefCorr)
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE sptils_setMatrixDefCorr (rsolverNode,Rmatrices)
  
!<description>
  ! This routine is called if the system matrix changes.
  ! The routine calls sptils_setMatrices for the preconditioner
  ! to inform also that one about the change of the matrix pointer.
!</description>
  
!<input>
  ! An array of system matrices which is simply passed to the initialisation 
  ! routine of the preconditioner.
  TYPE(t_ccoptSpaceTimeDiscretisation), DIMENSION(:), INTENT(IN) :: Rmatrices
!</input>
  
!<inputoutput>
  ! The t_sptilsNode structure of the solver
  TYPE(t_sptilsNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!</subroutine>

    ! Do we have a preconditioner given? If yes, call the general initialisation
    ! routine to initialise it.
    
    IF (ASSOCIATED(rsolverNode%p_rsubnodeDefCorr%p_rpreconditioner)) THEN
      CALL sptils_setMatrices (rsolverNode%p_rsubnodeDefCorr%p_rpreconditioner, &
                               Rmatrices)
    END IF

  END SUBROUTINE

! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE sptils_initStructureDefCorr (rsolverNode,ierror)
  
!<description>
  ! Solver preparation. Perform symbolic factorisation (not of the defect
  ! correcion solver, but of subsolvers). Allocate temporary memory.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure of the UMFPACK4 solver
  TYPE(t_sptilsNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to 
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  INTEGER, INTENT(OUT) :: ierror
!</output>
  
!</subroutine>
    
    ! A-priori we have no error...
    ierror = LINSOL_ERR_NOERROR
    
    ! Allocate memory for the two temp vectors.
    CALL sptivec_initVector (rsolverNode%p_rsubnodeDefCorr%rtempVector,&
        dof_igetNDofGlobBlock(rsolverNode%rspaceTimeDiscr%p_rlevelInfo%p_rdiscretisation),&
        rsolverNode%rspaceTimeDiscr%niterations)

    CALL sptivec_initVector (rsolverNode%p_rsubnodeDefCorr%rtempVector2,&
        dof_igetNDofGlobBlock(rsolverNode%rspaceTimeDiscr%p_rlevelInfo%p_rdiscretisation),&
        rsolverNode%rspaceTimeDiscr%niterations)

    ! We simply pass isubgroup to the subsolvers when calling them.
    ! Inside of this routine, there's not much to do with isubgroup,
    ! as there is no special information (like a factorisation)
    ! associated with this type of solver, which has to be allocated,
    ! released or updated...

    ! Call the init routine of the preconditioner.
    IF (ASSOCIATED(rsolverNode%p_rsubnodeDefCorr%p_rpreconditioner)) THEN
      CALL sptils_initStructure (rsolverNode%p_rsubnodeDefCorr%p_rpreconditioner,ierror)
    END IF
    
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE sptils_initDataDefCorr (rsolverNode, ierror)
  
!<description>
  ! Calls the initData subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with sptils_initData.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure of the UMFPACK4 solver
  TYPE(t_sptilsNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to 
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  INTEGER, INTENT(OUT) :: ierror
!</output>

!</subroutine>

    ! Call the init routine of the preconditioner.
    IF (ASSOCIATED(rsolverNode%p_rsubnodeDefCorr%p_rpreconditioner)) THEN
      CALL sptils_initData (rsolverNode%p_rsubnodeDefCorr%p_rpreconditioner,ierror)
    END IF
    
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE sptils_doneDataDefCorr (rsolverNode)
  
!<description>
  ! Calls the doneData subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with sptils_doneData.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure of the UMFPACK4 solver
  TYPE(t_sptilsNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!</subroutine>

    ! Call the done routine of the preconditioner.
    IF (ASSOCIATED(rsolverNode%p_rsubnodeDefCorr%p_rpreconditioner)) THEN
      CALL sptils_doneData (rsolverNode%p_rsubnodeDefCorr%p_rpreconditioner)
    END IF
    
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE sptils_doneStructureDefCorr (rsolverNode)
  
!<description>
  ! Calls the doneStructure subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with sptils_doneStructure.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure of the UMFPACK4 solver
  TYPE(t_sptilsNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!</subroutine>

    ! Release the temp vectors
    CALL sptivec_releaseVector(rsolverNode%p_rsubnodeDefCorr%rtempVector2)
    CALL sptivec_releaseVector(rsolverNode%p_rsubnodeDefCorr%rtempVector)

    ! Call the init routine of the preconditioner.
    IF (ASSOCIATED(rsolverNode%p_rsubnodeDefCorr%p_rpreconditioner)) THEN
      CALL sptils_doneStructure (rsolverNode%p_rsubnodeDefCorr%p_rpreconditioner)
    END IF
    
  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>
  
  RECURSIVE SUBROUTINE sptils_precDefCorr (rsolverNode,rd)
  
!<description>
  ! Applies the Defect Correction preconditioner $P \approx A$ to the defect 
  ! vector rd and solves $Pd_{new} = d$.
  ! rd will be overwritten by the preconditioned defect.
  !
  ! The matrix must have been attached to the system before calling
  ! this routine, and the initStructure/initData routines
  ! must have been called to prepare the solver for solving
  ! the problem.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure of the solver
  TYPE(t_sptilsNode), INTENT(INOUT), TARGET :: rsolverNode

  ! On call to this routine: The defect vector to be preconditioned.
  ! Will be overwritten by the preconditioned defect.
  TYPE(t_spacetimeVector), INTENT(INOUT)    :: rd
!</inputoutput>
  
!</subroutine>

  ! local variables
  INTEGER :: ite
  REAL(DP) :: dres,dfr
  TYPE(t_spacetimeVector), POINTER :: p_rx,p_rdef
  TYPE(t_sptilsNode), POINTER :: p_rprecSubnode
  
  ! Damping parameter
  REAL(DP) :: domega

  ! The system matrix
  TYPE(t_ccoptSpaceTimeDiscretisation), POINTER :: p_rmatrix
  
  ! Minimum number of iterations, print-sequence for residuals
  INTEGER :: nminIterations, niteResOutput
  
  ! Whether to filter/prcondition
  LOGICAL bprec
  
  ! The local subnode
  TYPE(t_sptilsSubnodeDefCorr), POINTER :: p_rsubnode

    ! Status reset
    rsolverNode%iresult = 0
    
    ! Getch some information
    p_rsubnode => rsolverNode%p_rsubnodeDefCorr
    p_rmatrix => rsolverNode%rspaceTimeDiscr

    ! Check the parameters
    IF (rd%ntimesteps .EQ. 0) THEN
    
      ! Parameters wrong
      rsolverNode%iresult = 2
      RETURN
    END IF

    ! Minimum number of iterations
 
    nminIterations = MAX(rsolverNode%nminIterations,0)
    
    ! Damping parameter
    domega = rsolverNode%domega
      
    ! Use preconditioning? Filtering?

    bprec = ASSOCIATED(rsolverNode%p_rsubnodeDefCorr%p_rpreconditioner)
    
    IF (bprec) THEN
      p_rprecSubnode => p_rsubnode%p_rpreconditioner
    END IF

    ! Iteration when the residuum is printed:

    niteResOutput = MAX(1,rsolverNode%niteResOutput)

    ! We want to perform <= nmaxIterations of the form
    !
    !    $$ x_{n+1}  =  x_n  +  \omega P^{-1} (b-Ax) $$
    !
    ! At first, set up our temporary vectors, which holds the current
    ! 'solution' and the current 'defect'. We already allocated them 
    ! during the initialisation phase - now we want to use them!
    
    p_rx   => rsolverNode%p_rsubnodeDefCorr%rtempVector
    p_rdef => rsolverNode%p_rsubnodeDefCorr%rtempVector2
    
    ! rd is our RHS. p_rx points to a new vector which will be our
    ! iteration vector. At the end of this routine, we replace
    ! rd by p_rx.
    ! Clear our iteration vector p_rx.
    CALL sptivec_clearVector (p_rx)
  
    ! Copy our RHS rd to p_rdef. As the iteration vector is 0, this
    ! is also our initial defect.

    CALL sptivec_copyVector(rd,p_rdef)
    IF (bprec) THEN
      ! Perform preconditioning with the assigned preconditioning
      ! solver structure.
      CALL sptils_precondDefect (p_rprecSubnode,p_rdef)
    END IF
    
    ! Get the norm of the residuum
    dres = sptivec_vectorNorm (p_rdef,rsolverNode%iresNorm)
    IF (.NOT.((dres .GE. 1E-99_DP) .AND. &
              (dres .LE. 1E99_DP))) dres = 0.0_DP

    ! Initialize starting residuum
      
    rsolverNode%dinitialDefect = dres

    ! Check if out initial defect is zero. This may happen if the filtering
    ! routine filters "everything out"!
    ! In that case we can directly stop our computation.

    IF ( rsolverNode%dinitialDefect .LT. rsolverNode%drhsZero ) THEN
     
      ! final defect is 0, as initialised in the output variable above

      CALL sptivec_clearVector(p_rx)
      ite = 0
      rsolverNode%dfinalDefect = dres
          
    ELSE

      IF (rsolverNode%ioutputLevel .GE. 2) THEN
        CALL output_line ('DefCorr: Iteration '// &
             TRIM(sys_siL(0,10))//',  !!RES!! = '//&
             TRIM(sys_sdEL(rsolverNode%dinitialDefect,15)) )
      END IF

      ! Perform at most nmaxIterations loops to get a new vector

      DO ite = 1,rsolverNode%nmaxIterations
      
        rsolverNode%icurrentIteration = ite
        
        ! In p_rdef, we now have the current residuum $P^{-1} (b-Ax)$.
        ! Add it (damped) to the current iterate p_x to get
        !   $$ x  :=  x  +  \omega P^{-1} (b-Ax) $$

        CALL sptivec_vectorLinearComb (p_rdef ,p_rx,domega,1.0_DP)

        ! Calculate the residuum for the next step : (b-Ax)
        CALL sptivec_copyVector (rd,p_rdef)
        CALL c2d2_spaceTimeMatVec (rsolverNode%p_rproblem, p_rmatrix, &
            p_rx,p_rdef, -1.0_DP,1.0_DP)
        IF (bprec) THEN
          ! Perform preconditioning with the assigned preconditioning
          ! solver structure.
          CALL sptils_precondDefect (p_rprecSubnode,p_rdef)
        END IF
        
        ! Get the norm of the new (final?) residuum
        dfr = sptivec_vectorNorm (p_rdef,rsolverNode%iresNorm)
     
        rsolverNode%dfinalDefect = dfr

        ! Test if the iteration is diverged
        IF (rsolverNode%depsRel .NE. SYS_INFINITY) THEN
          IF ( .NOT. (dfr .LE. rsolverNode%dinitialDefect*rsolverNode%ddivRel) ) THEN
            CALL output_line ('DefCorr: Solution diverging!')
            rsolverNode%iresult = 1
            EXIT
          END IF
        END IF
     
        ! At least perform nminIterations iterations
        IF (ite .GE. nminIterations) THEN
        
          ! Check if the iteration converged
          !
          ! Absolute convergence criterion? Check the norm directly.
          IF (rsolverNode%depsAbs .NE. 0.0_DP) THEN
            IF (.NOT. (dfr .GT. rsolverNode%depsAbs)) THEN
              EXIT
            END IF
          END IF
          
          ! Relative convergence criterion? Multiply with initial residuum
          ! and check the norm. 
          IF (rsolverNode%depsRel .NE. 0.0_DP) THEN
            IF (.NOT. (dfr .GT. rsolverNode%depsRel * rsolverNode%dinitialDefect)) THEN
              EXIT
            END IF
          END IF
          
        END IF

        ! print out the current residuum

        IF ((rsolverNode%ioutputLevel .GE. 2) .AND. &
            (MOD(ite,niteResOutput).EQ.0)) THEN
          CALL output_line ('DefCorr: Iteration '// &
              TRIM(sys_siL(ITE,10))//',  !!RES!! = '//&
              TRIM(sys_sdEL(rsolverNode%dfinalDefect,15)) )
        END IF

      END DO

      ! Set ITE to NIT to prevent printing of "NIT+1" of the loop was
      ! completed

      IF (ite .GT. rsolverNode%nmaxIterations) &
        ite = rsolverNode%nmaxIterations

      ! Finish - either with an error or if converged.
      ! Print the last residuum.

      IF ((rsolverNode%ioutputLevel .GE. 2) .AND. &
          (ite .GE. 1) .AND. (ITE .LT. rsolverNode%nmaxIterations) .AND. &
          (rsolverNode%iresult .GE. 0)) THEN
        CALL output_line ('DefCorr: Iteration '// &
            TRIM(sys_siL(ITE,10))//',  !!RES!! = '//&
            TRIM(sys_sdEL(rsolverNode%dfinalDefect,15)) )
      END IF

    END IF

    rsolverNode%iiterations = ite
    
    ! Overwrite our previous RHS by the new correction vector p_rx.
    ! This completes the preconditioning.
    CALL sptivec_copyVector (p_rx,rd)
      
    ! Don't calculate anything if the final residuum is out of bounds -
    ! would result in NaN's,...
      
    IF (rsolverNode%dfinalDefect .LT. 1E99_DP) THEN
    
      ! If the initial defect was zero, the solver immediately
      ! exits - and so the final residuum is zero and we performed
      ! no steps; so the resulting convergence rate stays zero.
      ! In the other case the convergence rate computes as
      ! (final defect/initial defect) ** 1/nit :

      IF (rsolverNode%dfinalDefect .GT. rsolverNode%drhsZero) THEN
        rsolverNode%dconvergenceRate = &
                    (rsolverNode%dfinalDefect / rsolverNode%dinitialDefect) ** &
                    (1.0_DP/REAL(rsolverNode%iiterations,DP))
      END IF
      
      IF (rsolverNode%ioutputLevel .GE. 2) THEN
        CALL output_lbrk()
        CALL output_line ('DefCorr statistics:')
        CALL output_lbrk()
        CALL output_line ('Iterations              : '//&
             TRIM(sys_siL(rsolverNode%iiterations,10)) )
        CALL output_line ('!!INITIAL RES!!         : '//&
             TRIM(sys_sdEL(rsolverNode%dinitialDefect,15)) )
        CALL output_line ('!!RES!!                 : '//&
             TRIM(sys_sdEL(rsolverNode%dfinalDefect,15)) )
        IF (rsolverNode%dinitialDefect .GT. rsolverNode%drhsZero) THEN     
          CALL output_line ('!!RES!!/!!INITIAL RES!! : '//&
            TRIM(sys_sdEL(rsolverNode%dfinalDefect / rsolverNode%dinitialDefect,15)) )
        ELSE
          CALL output_line ('!!RES!!/!!INITIAL RES!! : '//&
               TRIM(sys_sdEL(0.0_DP,15)) )
        END IF
        CALL output_lbrk ()
        CALL output_line ('Rate of convergence     : '//&
             TRIM(sys_sdEL(rsolverNode%dconvergenceRate,15)) )

      END IF

      IF (rsolverNode%ioutputLevel .EQ. 1) THEN
        CALL output_line (&
              'DefCorr: Iterations/Rate of convergence: '//&
              TRIM(sys_siL(rsolverNode%iiterations,10))//' /'//&
              TRIM(sys_sdEL(rsolverNode%dconvergenceRate,15)) )
      END IF
      
    ELSE
      ! DEF=Infinity; RHO=Infinity, set to 1
      rsolverNode%dconvergenceRate = 1.0_DP
      rsolverNode%dasymptoticConvergenceRate = 1.0_DP
    END IF  
  
  END SUBROUTINE
  
! *****************************************************************************
! Block Jacobi preconditioner
! *****************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE sptils_initBlockJacobi (rproblem,p_rsolverNode,rspatialPrecond)
  
!<description>
  ! Creates a t_sptilsNode solver structure for the block Jacobi preconditioner.
!</description>
  
!<input>
  ! The problem structure that defines the problem.
  ! A pointer to this is saved in the solver node, so the problem structure
  ! must not be released before the solver node is released.
  TYPE(t_problem), TARGET :: rproblem
  
  ! A spatial preconditioner structure that defines how to perform the preconditioning
  ! in each substep. A pointer to this is noted in p_rsolverNode, so rspatialPrecond
  ! should not be released before the solver is destroyed.
  TYPE(t_ccspatialPreconditioner), INTENT(IN), TARGET :: rspatialPrecond
!</input>
  
!<output>
  ! A pointer to a t_sptilsNode structure. Is set by the routine, any previous
  ! value of the pointer is destroyed.
  TYPE(t_sptilsNode), POINTER         :: p_rsolverNode
!</output>
  
!</subroutine>
  
    ! Create a default solver structure
    CALL sptils_initSolverGeneral(rproblem,p_rsolverNode)
    
    ! Initialise the type of the solver
    p_rsolverNode%calgorithm = SPTILS_ALG_BLOCKJACOBI
    
    ! Allocate a subnode for our solver.
    ! This initialises most of the variables with default values appropriate
    ! to this solver.
    ALLOCATE(p_rsolverNode%p_rsubnodeBlockJacobi)
    
    ! Attach the preconditioner if given. 
    p_rsolverNode%p_rsubnodeBlockJacobi%p_rspatialPreconditioner => rspatialPrecond
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE sptils_doneBlockJacobi (rsolverNode)
  
!<description>
  ! This routine releases all temporary memory for the efect correction solver 
  ! from the heap.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure of UMFPACK4 which is to be cleaned up.
  TYPE(t_sptilsNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!</subroutine>
  
    ! Release memory if still associated
    CALL sptils_doneDataBlockJacobi (rsolverNode)
    CALL sptils_doneStructureBlockJacobi (rsolverNode)
    
    ! Release the subnode structure
    DEALLOCATE(rsolverNode%p_rsubnodeBlockJacobi)
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE sptils_setMatrixBlockJacobi (rsolverNode,Rmatrices)
  
!<description>
  ! This routine is called if the system matrix changes.
  ! The routine calls sptils_setMatrices for the preconditioner
  ! to inform also that one about the change of the matrix pointer.
!</description>
  
!<input>
  ! An array of system matrices which is simply passed to the initialisation 
  ! routine of the preconditioner.
  TYPE(t_ccoptSpaceTimeDiscretisation), DIMENSION(:), INTENT(IN) :: Rmatrices
!</input>
  
!<inputoutput>
  ! The t_sptilsNode structure of the solver
  TYPE(t_sptilsNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!</subroutine>

    ! Do we have a preconditioner given? If yes, call the general initialisation
    ! routine to initialise it.
    
  END SUBROUTINE

! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE sptils_initStructureBlockJacobi (rsolverNode,ierror)
  
!<description>
  ! Solver preparation. Perform symbolic factorisation (not of the defect
  ! correcion solver, but of subsolvers). Allocate temporary memory.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure of the UMFPACK4 solver
  TYPE(t_sptilsNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to 
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  INTEGER, INTENT(OUT) :: ierror
!</output>
  
!</subroutine>
    
    ! A-priori we have no error...
    ierror = LINSOL_ERR_NOERROR
    
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE sptils_initDataBlockJacobi (rsolverNode, ierror)
  
!<description>
  ! Calls the initData subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with sptils_initData.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure of the UMFPACK4 solver
  TYPE(t_sptilsNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to 
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  INTEGER, INTENT(OUT) :: ierror
!</output>

!</subroutine>

    ! A-priori we have no error...
    ierror = LINSOL_ERR_NOERROR

  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE sptils_doneDataBlockJacobi (rsolverNode)
  
!<description>
  ! Calls the doneData subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with sptils_doneData.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure of the UMFPACK4 solver
  TYPE(t_sptilsNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!</subroutine>

  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE sptils_doneStructureBlockJacobi (rsolverNode)
  
!<description>
  ! Calls the doneStructure subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with sptils_doneStructure.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure of the UMFPACK4 solver
  TYPE(t_sptilsNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!</subroutine>

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>
  
  RECURSIVE SUBROUTINE sptils_precBlockJacobi (rsolverNode,rd)
  
!<description>
  ! Applies the Block Jacobi preconditioner $P \approx A$ to the defect 
  ! vector rd and solves $Pd_{new} = d$.
  ! rd will be overwritten by the preconditioned defect.
  !
  ! The term 'Block Jacobi' means that a spatial preconditioner is applied
  ! to the diagonal block of each time step in the global space-time
  ! system. Usually, this spatial preconditioner is a linear solver
  ! (multigrid or whatever). More precisely, any spatial preconditioner that
  ! can be used with c2d2_precondDefect is suitable.
  !
  ! The matrix must have been attached to the system before calling
  ! this routine, and the initStructure/initData routines
  ! must have been called to prepare the solver for solving
  ! the problem.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure of the solver
  TYPE(t_sptilsNode), INTENT(INOUT), TARGET :: rsolverNode

  ! On call to this routine: The defect vector to be preconditioned.
  ! Will be overwritten by the preconditioned defect.
  TYPE(t_spacetimeVector), INTENT(INOUT)    :: rd
!</inputoutput>
  
!</subroutine>

    ! local variables
    INTEGER :: isubstep
    REAL(DP) :: dtheta,dtstep
    LOGICAL :: bsuccess
    TYPE(t_ccmatrixComponents) :: rmatrixComponents
    TYPE(t_ccoptSpaceTimeDiscretisation), POINTER :: p_rspaceTimeDiscr
    TYPE(t_vectorBlock) :: rtempVectorD,rtempVectorX
    
    ! DEBUG!!!
    !REAL(DP), DIMENSION(:), POINTER :: p_Dx,p_Dd
    
    ! Get a pointer to the space-time discretisation structure that defines
    ! how to apply the global system matrix.
    p_rspaceTimeDiscr => rsolverNode%rspaceTimeDiscr
    
    dtheta = rsolverNode%p_rproblem%rtimedependence%dtimeStepTheta
    dtstep = p_rspaceTimeDiscr%dtstep

    ! Create temp vectors for X, B and D.
    CALL lsysbl_createVecBlockByDiscr (p_rspaceTimeDiscr%p_rlevelInfo%p_rdiscretisation,&
        rtempVectorD,.TRUE.)
    CALL lsysbl_createVecBlockByDiscr (p_rspaceTimeDiscr%p_rlevelInfo%p_rdiscretisation,&
        rtempVectorX,.TRUE.)
        
    ! Attach the boundary conditions to the temp vectors.
    rtempVectorD%p_rdiscreteBC => p_rspaceTimeDiscr%p_rlevelInfo%p_rdiscreteBC
    rtempVectorD%p_rdiscreteBCfict => p_rspaceTimeDiscr%p_rlevelInfo%p_rdiscreteFBC

    rtempVectorX%p_rdiscreteBC => p_rspaceTimeDiscr%p_rlevelInfo%p_rdiscreteBC
    rtempVectorX%p_rdiscreteBCfict => p_rspaceTimeDiscr%p_rlevelInfo%p_rdiscreteFBC

    ! The weights in the rmatrixComponents structure are later initialised
    ! according to the actual situation when the matrix is to be used.
    rmatrixComponents%p_rdiscretisation         => &
        p_rspaceTimeDiscr%p_rlevelInfo%p_rdiscretisation
    rmatrixComponents%p_rmatrixStokes           => &
        p_rspaceTimeDiscr%p_rlevelInfo%rmatrixStokes          
    rmatrixComponents%p_rmatrixB1             => &
        p_rspaceTimeDiscr%p_rlevelInfo%rmatrixB1              
    rmatrixComponents%p_rmatrixB2             => &
        p_rspaceTimeDiscr%p_rlevelInfo%rmatrixB2              
    rmatrixComponents%p_rmatrixMass           => &
        p_rspaceTimeDiscr%p_rlevelInfo%rmatrixMass            
    rmatrixComponents%p_rmatrixIdentityPressure => &
        p_rspaceTimeDiscr%p_rlevelInfo%rmatrixIdentityPressure
        
    rmatrixComponents%iupwind = &
        collct_getvalue_int (rsolverNode%p_rproblem%rcollection,'IUPWIND')
    rmatrixComponents%dnu = &
        collct_getvalue_real (rsolverNode%p_rproblem%rcollection,'NU')
    rmatrixComponents%dupsam = &
        collct_getvalue_real (rsolverNode%p_rproblem%rcollection,'UPSAM')
    
    ! ----------------------------------------------------------------------
    ! We use a block-Jacobi scheme for preconditioning...
    !
    ! For this purpose, loop through the substeps.
    
    DO isubstep = 0,p_rspaceTimeDiscr%niterations
    
      ! Current time step?
      rsolverNode%p_rproblem%rtimedependence%dtime = &
          rsolverNode%p_rproblem%rtimedependence%dtimeInit + isubstep * dtstep

      CALL output_separator (OU_SEP_MINUS)

      CALL output_line ('Block-Jacobi preconditioning of timestep: '//&
          TRIM(sys_siL(isubstep,10))//&
          ' Time: '//TRIM(sys_sdEL(rsolverNode%p_rproblem%rtimedependence%dtime,10)))
    
      ! -----
      ! Discretise the boundary conditions at the new point in time -- 
      ! if the boundary conditions are nonconstant in time!
      IF (collct_getvalue_int (rsolverNode%p_rproblem%rcollection,'IBOUNDARY') &
          .NE. 0) THEN
        CALL c2d2_updateDiscreteBC (rsolverNode%p_rproblem, .FALSE.)
      END IF

      ! DEBUG!!!      
      !CALL lsysbl_getbase_double (rtempVectorX,p_Dx)
      !CALL lsysbl_getbase_double (rtempVectorD,p_Dd)

      ! Read in the RHS/solution/defect vector of the current timestep.
      ! If no solution is specified, we have a linear problem and thus
      ! the content of rtempVector is not relevant; actually it's even
      ! zero by initialisation.
      IF (ASSOCIATED(p_rspaceTimeDiscr%p_rsolution)) THEN
        CALL sptivec_getTimestepData (p_rspaceTimeDiscr%p_rsolution, &
            isubstep, rtempVectorX)
      END IF
      CALL sptivec_getTimestepData (rd, isubstep, rtempVectorD)

      ! Set up the matrix weights for the diagonal matrix
      CALL c2d2_setupMatrixWeights (rsolverNode%p_rproblem,p_rspaceTimeDiscr,dtheta,&
        isubstep,0,rmatrixComponents)
        
      ! Perform preconditioning of the spatial defect with the method provided by the
      ! core equation module.
      CALL c2d2_precondDefect (&
          rsolverNode%p_rsubnodeBlockJacobi%p_rspatialPreconditioner,&
          rmatrixComponents,&
          rtempVectorD,rtempVectorX,bsuccess,rsolverNode%p_rproblem%rcollection)      
    
      ! Save back the preconditioned defect.
      CALL sptivec_setTimestepData (rd, isubstep, rtempVectorD)
      
    END DO
    
    CALL lsysbl_releaseVector (rtempVectorX)
    CALL lsysbl_releaseVector (rtempVectorD)
    
  END SUBROUTINE

! *****************************************************************************
! Multigrid preconditioner
! *****************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE sptils_initMultigrid (rproblem,NLMIN,NLMAX,p_rsolverNode)
  
!<description>
  ! Creates a t_sptilsNode solver structure for the multigrid preconditioner.
!</description>
  
!<input>
  ! The problem structure that defines the problem.
  ! A pointer to this is saved in the solver node, so the problem structure
  ! must not be released before the solver node is released.
  TYPE(t_problem), TARGET :: rproblem
  
  ! Minimum level; level where to apply the coarse grid solver.
  ! Must be >= rproblem%NLMIN!
  INTEGER, INTENT(IN) :: NLMIN
  
  ! Maximum level; level where the preconditioner is applied.
  ! Must be <= rproblem%NLMAX!
  INTEGER, INTENT(IN) :: NLMAX
!</input>
  
!<output>
  ! A pointer to a t_sptilsNode structure. Is set by the routine, any previous
  ! value of the pointer is destroyed.
  TYPE(t_sptilsNode), POINTER         :: p_rsolverNode
!</output>
  
!</subroutine>
  
    ! Create a default solver structure
    CALL sptils_initSolverGeneral(rproblem,p_rsolverNode)
    
    ! Initialise the type of the solver
    p_rsolverNode%calgorithm = SPTILS_ALG_MULTIGRID
    
    ! Allocate a subnode for our solver.
    ! This initialises most of the variables with default values appropriate
    ! to this solver.
    ALLOCATE(p_rsolverNode%p_rsubnodeMultigrid)
    
    ! Set parameters and allocate memory for the level information.
    p_rsolverNode%p_rsubnodeMultigrid%NLMIN = NLMIN
    p_rsolverNode%p_rsubnodeMultigrid%NLMAX = NLMAX
    ALLOCATE(p_rsolverNode%p_rsubnodeMultigrid%p_Rlevels(1:NLMAX-NLMIN+1))
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE sptils_doneMultigrid (rsolverNode)
  
!<description>
  ! This routine releases all temporary memory for the efect correction solver 
  ! from the heap.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure of UMFPACK4 which is to be cleaned up.
  TYPE(t_sptilsNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!</subroutine>
  
    ! Release memory if still associated
    CALL sptils_doneDataMultigrid (rsolverNode)
    CALL sptils_doneStructureMultigrid (rsolverNode)
    
    ! Release the subnode structure
    DEALLOCATE(rsolverNode%p_rsubnodeMultigrid%p_Rlevels)
    DEALLOCATE(rsolverNode%p_rsubnodeMultigrid)
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE sptils_setMatrixMultigrid (rsolverNode,Rmatrices)
  
!<description>
  ! This routine is called if the system matrix changes.
  ! The routine calls sptils_setMatrices for the preconditioner
  ! to inform also that one about the change of the matrix pointer.
!</description>
  
!<input>
  ! An array of system matrices which is simply passed to the initialisation 
  ! routine of the preconditioner.
  TYPE(t_ccoptSpaceTimeDiscretisation), DIMENSION(:), INTENT(IN) :: Rmatrices
!</input>
  
!<inputoutput>
  ! The t_sptilsNode structure of the solver
  TYPE(t_sptilsNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!</subroutine>

    INTEGER :: ilev

    ! Loop through the level. 
    DO ilev=1,SIZE(Rmatrices)
      ! On each level, set the matrix.
      rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%rspaceTimeDiscr = &
          Rmatrices(ilev)
      
      ! Call the setmatrices-routine for the presmoother/postsmoother/
      ! coarse grid solver
      IF (ASSOCIATED(rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%p_rpresmoother)) THEN
        CALL sptils_setMatrices(&
            rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%p_rpresmoother,&
            Rmatrices(1:ilev))
      END IF
      
      IF (ASSOCIATED(rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%p_rpostsmoother)) THEN
        CALL sptils_setMatrices(&
            rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%p_rpostsmoother,&
            Rmatrices(1:ilev))
      END IF
      
      IF (ASSOCIATED(rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%p_rcoarseGridSolver)) THEN
        CALL sptils_setMatrices(&
            rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%p_rcoarseGridSolver,&
            Rmatrices(1:ilev))
      END IF
    END DO
    
  END SUBROUTINE

! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE sptils_initStructureMultigrid (rsolverNode,ierror)
  
!<description>
  ! Solver preparation. Perform symbolic factorisation (not of the defect
  ! correcion solver, but of subsolvers). Allocate temporary memory.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure of the UMFPACK4 solver
  TYPE(t_sptilsNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to 
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  INTEGER, INTENT(OUT) :: ierror
!</output>
  
!</subroutine>

    ! local variables
    INTEGER :: ilev,NLMAX,ntimesteps
    INTEGER(PREC_VECIDX) :: NEQ
    
    ! A-priori we have no error...
    ierror = LINSOL_ERR_NOERROR
    
    NLMAX = SIZE(rsolverNode%p_rsubnodeMultigrid%p_Rlevels)

    ! On each level, call the initStructure routine of the
    ! presmoother/postsmoother/coarse grid solver


    ! Loop through the level. 
    DO ilev=1,NLMAX
      ! Call the setmatrices-routine for the presmoother/postsmoother/
      ! coarse grid solver
      IF (ASSOCIATED(rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%p_rpresmoother)) THEN
        CALL sptils_initStructure(&
            rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%p_rpresmoother,ierror)
        IF (ierror .NE. ierror) RETURN
      END IF
      
      IF (ASSOCIATED(rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%p_rpostsmoother)) THEN
        CALL sptils_initStructure(&
            rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%p_rpostsmoother,ierror)
        IF (ierror .NE. ierror) RETURN
      END IF
      
      IF (ASSOCIATED(rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%p_rcoarseGridSolver)) THEN
        CALL sptils_initStructure(&
            rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%p_rcoarseGridSolver,ierror)
        IF (ierror .NE. ierror) RETURN
      END IF
      
      ! Generate an interlevel projection structure for that level
      CALL sptipr_initProjection (&
          rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%rinterlevelProjection,&
          rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%&
              rspaceTimeDiscr%p_rlevelInfo%p_rdiscretisation)
              
      ntimesteps = rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%&
              rspaceTimeDiscr%niterations
      NEQ = dof_igetNDofGlobBlock(rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%&
              rspaceTimeDiscr%p_rlevelInfo%p_rdiscretisation)
              
      ! On all levels except for the maximum one, create a solution vector
      IF (ilev .LT. NLMAX) THEN
        CALL sptivec_initVector (&
            rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%rsolutionVector,&
            NEQ,ntimesteps)
      END IF
      
      ! On all levels except for the first one, create a RHS and a temp vector
      IF (ilev .GT. 1) THEN
        CALL sptivec_initVector (&
            rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%rrhsVector,&
            NEQ,ntimesteps)
        CALL sptivec_initVector (&
            rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%rtempVector,&
            NEQ,ntimesteps)
      END IF
      
      ! Create a block temp vector for the interlevel projection
      CALL lsysbl_createVecBlockByDiscr(&
          rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%rspaceTimeDiscr%&
          p_rlevelInfo%p_rdiscretisation,&
          rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%rprjVector,.FALSE.)

    END DO
    
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE sptils_initDataMultigrid (rsolverNode, ierror)
  
!<description>
  ! Calls the initData subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with sptils_initData.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure of the UMFPACK4 solver
  TYPE(t_sptilsNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to 
  ! LINSOL_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  INTEGER, INTENT(OUT) :: ierror
!</output>

!</subroutine>

    ! local variables
    INTEGER :: ilev
    
    ! A-priori we have no error...
    ierror = LINSOL_ERR_NOERROR
    
    ! On each level, call the initStructure routine of the
    ! presmoother/postsmoother/coarse grid solver

    ! Loop through the level. 
    DO ilev=1,SIZE(rsolverNode%p_rsubnodeMultigrid%p_Rlevels)
      ! Call the setmatrices-routine for the presmoother/postsmoother/
      ! coarse grid solver
      IF (ASSOCIATED(rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%p_rpresmoother)) THEN
        CALL sptils_initData(&
            rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%p_rpresmoother,ierror)
        IF (ierror .NE. ierror) RETURN
      END IF
      
      IF (ASSOCIATED(rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%p_rpostsmoother)) THEN
        CALL sptils_initData(&
            rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%p_rpostsmoother,ierror)
        IF (ierror .NE. ierror) RETURN
      END IF
      
      IF (ASSOCIATED(rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%p_rcoarseGridSolver)) THEN
        CALL sptils_initData(&
            rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%p_rcoarseGridSolver,ierror)
        IF (ierror .NE. ierror) RETURN
      END IF
    END DO

  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE sptils_doneDataMultigrid (rsolverNode)
  
!<description>
  ! Calls the doneData subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with sptils_doneData.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure of the UMFPACK4 solver
  TYPE(t_sptilsNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!</subroutine>

    ! local variables
    INTEGER :: ilev
    
    ! On each level, call the initStructure routine of the
    ! presmoother/postsmoother/coarse grid solver

    ! Loop through the level. 
    DO ilev=1,SIZE(rsolverNode%p_rsubnodeMultigrid%p_Rlevels)
      ! Call the setmatrices-routine for the presmoother/postsmoother/
      ! coarse grid solver
      IF (ASSOCIATED(rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%p_rpresmoother)) THEN
        CALL sptils_doneData(&
            rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%p_rpresmoother)
      END IF
      
      IF (ASSOCIATED(rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%p_rpostsmoother)) THEN
        CALL sptils_doneData(&
            rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%p_rpostsmoother)
      END IF
      
      IF (ASSOCIATED(rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%p_rcoarseGridSolver)) THEN
        CALL sptils_doneData(&
            rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%p_rcoarseGridSolver)
      END IF
    END DO

  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE sptils_doneStructureMultigrid (rsolverNode)
  
!<description>
  ! Calls the doneStructure subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with sptils_doneStructure.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure of the UMFPACK4 solver
  TYPE(t_sptilsNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!</subroutine>

    ! local variables
    INTEGER :: ilev,NLMAX
    
    ! On each level, call the initStructure routine of the
    ! presmoother/postsmoother/coarse grid solver

    ! Loop through the level. 
    NLMAX = SIZE(rsolverNode%p_rsubnodeMultigrid%p_Rlevels)
    DO ilev=1,NLMAX
      ! Call the setmatrices-routine for the presmoother/postsmoother/
      ! coarse grid solver
      IF (ASSOCIATED(rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%p_rpresmoother)) THEN
        CALL sptils_doneStructure(&
            rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%p_rpresmoother)
      END IF
      
      IF (ASSOCIATED(rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%p_rpostsmoother)) THEN
        CALL sptils_doneStructure(&
            rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%p_rpostsmoother)
      END IF
      
      IF (ASSOCIATED(rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%p_rcoarseGridSolver)) THEN
        CALL sptils_doneStructure(&
            rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%p_rcoarseGridSolver)
      END IF


      ! Release the projection structure
      CALL sptipr_doneProjection (&
          rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%rinterlevelProjection)
              
      ! Release vectors
              
      IF (ilev .LT. NLMAX) THEN
        CALL sptivec_releaseVector (&
            rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%rsolutionVector)
      END IF
      
      IF (ilev .GT. 1) THEN
        CALL sptivec_releaseVector (&
            rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%rrhsVector)
        CALL sptivec_releaseVector (&
            rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%rtempVector)
      END IF

      ! Release the temp vector for prolongation/restriction
      CALL lsysbl_releaseVector (&
        rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%rprjVector)

    END DO

  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE linsol_setMultigridLevel (rsolverNode,ilevel,&
                    rinterlevelProjection,&
                    p_rpresmoother,p_rpostsmoother,p_rcoarseGridSolver)
                    
!<description>
  ! This routine initialises the data of a multigrid level in the multigrid 
  ! solver rsolverNode. 
  ! The given coarse-grid solver and smoother-structures are attached
  ! to the level. 
  !
  ! It's allowed to use p_rpresmoother=p_rpostsmoother. A value of
  ! NULL() is also permissable, this deactivates the corresponding smoother.
  !
  ! It's allowed to use the same rinterlevelProjection on all levels,
  ! since the structure is level independent (as long as the
  ! spatial discretisation structures on different levels are 'compatible'
  ! what they have to be anyway).
!</description>
  
!<inputoutput>
  ! The solver structure of the multigrid solver, where the level
  ! should be added to. Must already be initialised for the multigrid
  ! solver.
  TYPE(t_sptilsNode)                             :: rsolverNode
!</inputoutput>
  
!<input>
  ! Level in the MG solver to be initialised.
  INTEGER, INTENT(IN) :: ilevel

  ! An interlevel projection structure that configures the projection
  ! between the solutions on a finer and a coarser grid. The structure
  ! must have been initialised with mlprj_initProjection.
  !
  ! Note that this structure is level-independent (as long as the
  ! spatial discretisation structures on different levels are 'compatible'
  ! what they have to be anyway), so the same structure can be used
  ! to initialise all levels!
  TYPE(t_sptiProjection) :: rinterlevelProjection
  
  ! Optional: A pointer to the solver structure of a solver that should be 
  ! used for presmoothing. This structure is used as a template to create an
  ! appropriate solver structure for all the levels. The structure itself is
  ! used at the finest level.
  ! If not given or set to NULL(), no presmoother will be used.
  TYPE(t_sptilsNode), POINTER, OPTIONAL   :: p_rpresmoother
  
  ! Optional: A pointer to the solver structure of a solver that should be 
  ! used for postsmoothing. This structure is used as a template to create an
  ! appropriate solver structure for all the levels. The structure itself is
  ! used at the finest level.
  ! If not given or set to NULL(), no presmoother will be used.
  TYPE(t_sptilsNode), POINTER, OPTIONAL   :: p_rpostsmoother

  ! Optional: A pointer to the solver structure of a solver that should be 
  ! used for coarse grid solving. 
  ! Should only be given for the very first level.
  TYPE(t_sptilsNode), POINTER, OPTIONAL   :: p_rcoarseGridSolver
!</input>
  
!</subroutine>

    ! local variables
    TYPE(t_sptilsMGLevelInfo), POINTER :: p_rlevelInfo
    
    ! Make sure the solver node is configured for multigrid
    IF ((rsolverNode%calgorithm .NE. LINSOL_ALG_MULTIGRID) .OR. &
        (.NOT. ASSOCIATED(rsolverNode%p_rsubnodeMultigrid))) THEN
      PRINT *,'Error: Multigrid structure not initialised'
      CALL sys_halt()
    END IF
    
    p_rlevelInfo => rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilevel)
    
    ! Attach the sub-solvers
    IF (ASSOCIATED(p_rcoarseGridSolver)) THEN
      p_rlevelInfo%p_rcoarseGridSolver => p_rcoarseGridSolver
    END IF

    IF (ASSOCIATED(p_rpresmoother)) THEN
      p_rlevelInfo%p_rpresmoother => p_rpresmoother
    END IF

    IF (ASSOCIATED(p_rpostsmoother)) THEN
      p_rlevelInfo%p_rpostsmoother => p_rpostsmoother
    END IF
    
    ! Initialise the interlevel projection structure,
    ! copy the data of our given template.
    p_rlevelInfo%rinterlevelProjection = rinterlevelProjection

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>
  
  SUBROUTINE sptils_convertToSmoother (rsolverNode,nsmoothingSteps,domega)
  
!<description>
  ! Converts a t_linsolNode to a smoother structure. A smoother is a solver
  ! that performs a fixed number of iterations without respecting any
  ! residuum. nsmoothingSteps is the number of steps the smoother should
  ! perform.
!</description>
  
!<input>
  ! Number of steps the smoother should perform
  INTEGER, INTENT(IN)          :: nsmoothingSteps
  
  ! OPTIONAL: Damping parameter.
  ! The parameter rsolverNode%domega is set to this value in order
  ! to set up the damping.
  REAL(DP), INTENT(IN), OPTIONAL :: domega
!</input>
  
!<inputoutput>
  ! The t_sptilsNode structure of the solver which should be configured as smoother.
  TYPE(t_sptilsNode), INTENT(INOUT), TARGET :: rsolverNode
!</inputoutput>
  
!</subroutine>

    rsolverNode%depsRel = 0.0_DP
    rsolverNode%depsAbs = 0.0_DP
    rsolverNode%nminIterations = nsmoothingSteps
    rsolverNode%nmaxIterations = nsmoothingSteps
    
    IF (PRESENT(domega)) rsolverNode%domega = domega

  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE sptils_smoothCorrection (rsolverNode,rx,rb,rtemp)
  
!<description>
  ! This routine performs a smoothing process on the vector rx
  ! belonging to a linear system $Ax=b$.
  ! rsolverNode identifies a solver structure that is converted to a
  ! smoother using linsol_convertToSmoother: The number of smoothing
  ! steps is written to rsolverNode%nmaxIterations and 
  ! rsolverNode%nmaxIterations so that the solver performs a definite
  ! number of iterations regardless of the residual.
  !
  ! rx is assumed to be of 'defect' type. There are two types of smoothing
  ! processes, depending on which type of preconditioner rsolverNode is:
  ! 1.) If rsolverNode is an iterative solver, linsol_smoothCorrection
  !     calls the associated solver P to compute $x=P^{-1}b$ using
  !     rsolverNode%nmaxIterations iterations.
  ! 2.) If rsolverNode is a 1-step solver, linsol_smoothCorrection
  !     performs rsolverNode%nmaxIterations defect correction steps
  !      $$ x := x + P^{-1} (b-Ax) $$
  !     with $x_0=0$.
  !
  ! The matrix must have been attached to the system before calling
  ! this routine, and the initStructure/initData routines
  ! must have been called to prepare the solver for solving
  ! the problem.
!</description>
  
!<input>
  ! The RHS vector of the system
  TYPE(t_spacetimeVector), INTENT(IN), TARGET        :: rb
!</input>
  
!<inputoutput>
  ! The t_sptilsNode structure of the solver
  TYPE(t_sptilsNode), INTENT(INOUT), TARGET :: rsolverNode
  
  ! The initial solution vector; receives the solution of the system
  TYPE(t_spacetimeVector), INTENT(INOUT)                 :: rx
  
  ! A temporary vector of the same size and structure as rx.
  TYPE(t_spacetimeVector), INTENT(INOUT)                 :: rtemp
!</inputoutput>
  
!</subroutine>

    INTEGER :: i
    INTEGER :: iiterations
    REAL(DP) :: dres
    TYPE(t_ccoptSpaceTimeDiscretisation), POINTER :: p_rmatrix
    !DEBUG: REAL(DP), DIMENSION(:), POINTER :: p_Ddata,p_Ddata2
    
    ! Cancel if nmaxIterations = number of smoothing steps is =0.
    IF (rsolverNode%nmaxIterations .LE. 0) RETURN
    
    ! This is a 1-step solver, we have to emulate the smoothing
    ! iterations. Perform rsolverNode%nmaxIterations steps of the
    ! form
    !     $$ x_{n+1} = x_n + P^{-1}(b-Ax_n) $$
    ! with $x_0 = 0$.
    
    !DEBUG: CALL lsysbl_getbase_double (rx,p_Ddata)
    !DEBUG: CALL lsysbl_getbase_double (rtemp,p_Ddata2)
    
    ! Apply nmaxIterations times defect correction to the given solution rx.
    p_rmatrix => rsolverNode%rspaceTimeDiscr
    
    iiterations = rsolverNode%nmaxIterations
    DO i=1,iiterations
    
      CALL sptivec_copyVector(rb,rtemp)
      CALL c2d2_spaceTimeMatVec (rsolverNode%p_rproblem, p_rmatrix, &
          rx,rtemp, -1.0_DP,1.0_DP,dres)
      
      IF (rsolverNode%ioutputLevel .GE. 2) THEN
        IF (.NOT.((dres .GE. 1E-99_DP) .AND. (dres .LE. 1E99_DP))) dres = 0.0_DP
                  
        CALL output_line ('Smoother: Step '//TRIM(sys_siL(i-1,10))//&
            ' !!RES!! = '//TRIM(sys_sdEL(dres,15)) )
      END IF
      
      CALL sptils_precondDefect(rsolverNode,rtemp)
      CALL sptivec_vectorLinearComb (rtemp,rx,1.0_DP,1.0_DP)
      
    END DO

    ! Probably print the final residuum
    IF (rsolverNode%ioutputLevel .GE. 2) THEN
      CALL sptivec_copyVector(rb,rtemp)
      CALL c2d2_spaceTimeMatVec (rsolverNode%p_rproblem, p_rmatrix, &
          rx,rtemp, -1.0_DP,1.0_DP,dres)
      
      IF (.NOT.((dres .GE. 1E-99_DP) .AND. (dres .LE. 1E99_DP))) dres = 0.0_DP
                
      CALL output_line ('Smoother: Step '//TRIM(sys_siL(i-1,10))//&
          ' !!RES!! = '//TRIM(sys_sdEL(dres,15)) )
    END IF
    
  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>
  
  RECURSIVE SUBROUTINE sptils_precMultigrid (rsolverNode,rd)
  
!<description>
  ! Applies the Multigrid preconditioner $P \approx A$ to the defect 
  ! vector rd and solves $Pd_{new} = d$.
  ! rd will be overwritten by the preconditioned defect.
  !
  ! The matrices must have been attached to the system before calling
  ! this routine, and the initStructure/initData routines
  ! must have been called to prepare the solver for solving
  ! the problem.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure of the solver
  TYPE(t_sptilsNode), INTENT(INOUT), TARGET :: rsolverNode

  ! On call to this routine: The defect vector to be preconditioned.
  ! Will be overwritten by the preconditioned defect.
  TYPE(t_spacetimeVector), INTENT(INOUT)    :: rd
!</inputoutput>
  
!</subroutine>

    !CALL sptils_precondDefect (rsolverNode%p_rsubnodeMultigrid%p_Rlevels(1)%p_rcoarseGridSolver,&
    !  rd)

  ! local variables
  INTEGER :: nminIterations,nmaxIterations,niteResOutput
  INTEGER :: ite,ilev,nlmax,i,nlmin
  REAL(DP) :: dres,dstep
  
  ! The system matrix on the current level
  TYPE(t_ccoptSpaceTimeDiscretisation), POINTER :: p_rmatrix
  
  ! Our MG structure
  TYPE(t_sptilsSubnodeMultigrid), POINTER :: p_rsubnode
  
    ! Solve the system!
    
    ! Getch some information
    p_rsubnode => rsolverNode%p_rsubnodeMultigrid
    
    ! Get the system matrix on the finest level:
    p_rmatrix => rsolverNode%rspaceTimeDiscr

    IF (p_rsubnode%icycle .LT. 0) THEN
      ! Wrong cycle
      PRINT *,'Multigrid: Wrong cycle!'
      rsolverNode%iresult = 2
      RETURN
    END IF
    
    ! Minimum/maximum number of iterations
    nminIterations = MAX(rsolverNode%nminIterations,0)
    nmaxIterations = MAX(rsolverNode%nmaxIterations,0)
      
    ! Iteration when the residuum is printed:
    niteResOutput = MAX(1,rsolverNode%niteResOutput)

    ! Status reset
    rsolverNode%iresult = 0
    rsolverNode%icurrentIteration = 0
    rsolverNode%dinitialDefect = 0.0_DP
    rsolverNode%dfinalDefect = 0.0_DP
    rsolverNode%dconvergenceRate = 0.0_DP
    rsolverNode%dasymptoticConvergenceRate = 0.0_DP
    
    ! We start on the maximum level. 
    ilev = p_rsubnode%NLMAX
    nlmax = p_rsubnode%NLMAX
    nlmin = p_rsubnode%NLMIN
    
    ! Is there only one level? Can be seen if the current level
    ! already contains a coarse grid solver.
    IF (nlmin .EQ. nlmax) THEN
    
      IF (rsolverNode%ioutputLevel .GT. 1) THEN
        CALL output_line ('Multigrid: Only one level. '//&
             'Switching back to standard solver.')
      END IF
      CALL sptils_precondDefect(p_rsubnode%p_Rlevels(ilev)%p_rcoarseGridSolver,rd)
      
      ! Take the statistics from the coarse grid solver.
      rsolverNode%dinitialDefect = &
        p_rsubnode%p_Rlevels(ilev)%p_rcoarseGridSolver%dinitialDefect
      rsolverNode%dfinalDefect = &
        p_rsubnode%p_Rlevels(ilev)%p_rcoarseGridSolver%dfinalDefect
      rsolverNode%dconvergenceRate = &
        p_rsubnode%p_Rlevels(ilev)%p_rcoarseGridSolver%dconvergenceRate
      rsolverNode%iiterations = &
        p_rsubnode%p_Rlevels(ilev)%p_rcoarseGridSolver%iiterations
      
    ELSE
      ! Get the norm of the initial residuum.
      ! As the initial iteration vector is zero, this is the norm
      ! of the RHS:
      dres = sptivec_vectorNorm (rd,rsolverNode%iresNorm)
      IF (.NOT.((dres .GE. 1E-99_DP) .AND. &
                (dres .LE. 1E99_DP))) dres = 0.0_DP

      ! Initialize starting residuum
      rsolverNode%dinitialDefect = dres

      ! Check if out initial defect is zero. This may happen if the filtering
      ! routine filters "everything out"!
      ! In that case we can directly stop our computation.

      IF ( rsolverNode%dinitialDefect .LT. rsolverNode%drhsZero ) THEN
      
        ! final defect is 0, as initialised in the output variable above
        CALL sptivec_clearVector(rd)
        rsolverNode%dfinalDefect = dres
        rsolverNode%dfinalDefect = dres
        rsolverNode%dconvergenceRate = 0.0_DP
        rsolverNode%iiterations = 0
        
      ELSE
        
        ! At first, reset the cycle counters of all levels
        DO i=nlmin,nlmax
          IF (p_rsubnode%icycle .EQ. 0) THEN
            p_rsubnode%p_Rlevels(i)%ncycles = 2
          ELSE
            p_rsubnode%p_Rlevels(i)%ncycles = p_rsubnode%icycle
          END IF  
        END DO
        
        ! Afterwards, p_rcurrentLevel points to the maximum level again.
        ! There, we set ncycles to 1.
        p_rsubnode%p_Rlevels(ilev)%ncycles = 1
        
        ! Print out the initial residuum

        IF (rsolverNode%ioutputLevel .GE. 2) THEN
          CALL output_line ('Multigrid: Iteration '// &
              TRIM(sys_siL(0,10))//',  !!RES!! = '//&
              TRIM(sys_sdEL(rsolverNode%dinitialDefect,15)) )
        END IF
        
        ! Copy the initial RHS to the RHS vector on the maximum level.
        CALL sptivec_copyVector(rd,p_rsubnode%p_Rlevels(ilev)%rrhsVector)
        
        ! Replace the solution vector on the finest level by rd.
        ! Afterwards, rd and the solution vector on the finest level
        ! share the same memory location, so we use rd as iteration
        ! vector directly.
        p_rsubnode%p_Rlevels(ilev)%rsolutionVector = rd
        
        ! Clear the initial solution vector.
        CALL sptivec_clearVector (p_rsubnode%p_Rlevels(ilev)%rsolutionVector)
        
        ! Start multigrid iteration; perform at most nmaxiterations iterations.
        DO ite = 1, nmaxiterations
        
          rsolverNode%icurrentIteration = ite
          
          ! Initialize cycle counters for all levels.
          DO i=nlmin,nlmax
            p_rsubnode%p_Rlevels(i)%ncyclesRemaining = &
                p_rsubnode%p_Rlevels(i)%ncycles
          END DO
        
          ! p_rcurrentLevel now points to the maximum level.
          ilev = nlmax

          ! Get the system matrix on the finest level:
          p_rmatrix => p_rsubnode%p_Rlevels(ilev)%rspaceTimeDiscr
          
          IF (rsolverNode%ioutputLevel .GE. 3) THEN
            CALL output_line ('Multigrid: Current mesh level: '//TRIM(sys_siL(ilev,5)))
          END IF
          
          ! Build the defect...
          CALL sptivec_copyVector (&
              p_rsubnode%p_Rlevels(ilev)%rrhsVector,&
              p_rsubnode%p_Rlevels(ilev)%rtempVector)
          IF (ite .NE. 1) THEN   ! initial solution vector is zero!
            CALL c2d2_spaceTimeMatVec (rsolverNode%p_rproblem, p_rmatrix, &
                p_rsubnode%p_Rlevels(ilev)%rsolutionVector,&
                p_rsubnode%p_Rlevels(ilev)%rtempVector, -1.0_DP,1.0_DP)
          END IF
          
          cycleloop: DO  ! Loop for the cycles
          
            ! On the maximum level we already built out defect vector. If we are    
            ! on a lower level than NLMAX, perform smoothing+restriction down to the
            ! coarse level. We identify the coarse level by checking if
            ! the current level has a coarse grid solver.
            
            DO WHILE (ilev .GT. nlmin)
            
              ! Perform the pre-smoothing with the current solution vector
              IF (ASSOCIATED(p_rsubnode%p_Rlevels(ilev)%p_rpreSmoother)) THEN
                CALL sptils_smoothCorrection (&
                          p_rsubnode%p_Rlevels(ilev)%p_rpreSmoother,&
                          p_rsubnode%p_Rlevels(ilev)%rsolutionVector,&
                          p_rsubnode%p_Rlevels(ilev)%rrhsVector,&
                          p_rsubnode%p_Rlevels(ilev)%rtempVector)
              END IF
            
              ! Build the defect vector
              CALL sptivec_copyVector (p_rsubnode%p_Rlevels(ilev)%rrhsVector,&
                  p_rsubnode%p_Rlevels(ilev)%rtempVector)
              CALL c2d2_spaceTimeMatVec (rsolverNode%p_rproblem, p_rmatrix, &
                  p_rsubnode%p_Rlevels(ilev)%rsolutionVector,&
                  p_rsubnode%p_Rlevels(ilev)%rtempVector, -1.0_DP,1.0_DP,dres)
              
              ! Extended output
              IF (ASSOCIATED(p_rsubnode%p_Rlevels(ilev)%p_rpreSmoother) .AND. &
                  (rsolverNode%ioutputLevel .GE. 3) .AND. &
                  (MOD(ite,niteResOutput) .EQ. 0)) THEN
                  
                IF (.NOT.((dres .GE. 1E-99_DP) .AND. &
                          (dres .LE. 1E99_DP))) dres = 0.0_DP
                          
                CALL output_line ('Multigrid: Level '//TRIM(sys_siL(ilev,5))//&
                    ' after presm.:     !!RES!! = '//TRIM(sys_sdEL(dres,15)) )
              END IF

              ! Restriction of the defect. The restricted defect is placed
              ! in the right hand side vector of the lower level.
              ! The projection parameters from level ilev to level
              ! ilev-1 is configured in the rprojection structure of the
              ! current level.
              !
              ! When restricting to the coarse grid, we directy restrict into
              ! the solution vector. It's used there as RHS and replaced in-situ
              ! by the solution by the coarse grid solver. So don't need the RHS vector
              ! on the coarse grid and save one vector-copy.
              ! 
              ! Otherwise, we restrict to the RHS on the lower level and continue
              ! the smoothing process there.
              IF (ilev .GT. nlmin+1) THEN
                ! We don't project to the coarse grid
                CALL sptipr_performRestriction (&
                      p_rsubnode%p_Rlevels(ilev)%rinterlevelProjection,&
                      p_rsubnode%p_Rlevels(ilev-1)%rrhsVector, &
                      p_rsubnode%p_Rlevels(ilev)%rtempVector, &
                      p_rsubnode%p_Rlevels(ilev-1)%rprjVector, &
                      p_rsubnode%p_Rlevels(ilev)%rprjVector)

                ! Choose zero as initial vector on lower level. 
                CALL sptivec_clearVector (p_rsubnode%p_Rlevels(ilev-1)%rsolutionVector)
                
                ! Extended output and/or adaptive cycles
                IF (rsolverNode%ioutputLevel .GE. 3) THEN
                
                  dres = sptivec_vectorNorm (p_rsubnode%p_Rlevels(ilev-1)%rrhsVector,&
                      rsolverNode%iresNorm)
                  IF (.NOT.((dres .GE. 1E-99_DP) .AND. &
                            (dres .LE. 1E99_DP))) dres = 0.0_DP
                            
                  ! If the output level is high enough, print that residuum norm.   
                  IF (MOD(ite,niteResOutput).EQ. 0) THEN
                    CALL output_line ('Multigrid: Level '//TRIM(sys_siL(ilev-1,5))//&
                        ' after restrict.:  !!RES!! = '//TRIM(sys_sdEL(dres,15)) )
                  END IF
                  
                END IF

              ELSE
              
                ! The vector is to be restricted to the coarse grid.
                CALL sptipr_performRestriction (&
                      p_rsubnode%p_Rlevels(ilev)%rinterlevelProjection,&
                      p_rsubnode%p_Rlevels(ilev-1)%rsolutionVector, &
                      p_rsubnode%p_Rlevels(ilev)%rtempVector, &
                      p_rsubnode%p_Rlevels(ilev-1)%rprjVector, &
                      p_rsubnode%p_Rlevels(ilev)%rprjVector)

                ! Extended output
                IF ((rsolverNode%ioutputLevel .GE. 3) .AND. &
                    (MOD(ite,niteResOutput).EQ.0)) THEN
                    
                  dres = sptivec_vectorNorm (p_rsubnode%p_Rlevels(ilev-1)%rsolutionVector,&
                      rsolverNode%iresNorm)
                  IF (.NOT.((dres .GE. 1E-99_DP) .AND. &
                            (dres .LE. 1E99_DP))) dres = 0.0_DP
                            
                  CALL output_line ('Multigrid: Level '//TRIM(sys_siL(ilev-1,5))//&
                      ' after restrict.:  !!RES!! = '//TRIM(sys_sdEL(dres,15)) )
                END IF

              END IF              
            
              ! Go down one level
              ilev = ilev - 1
              p_rmatrix => p_rsubnode%p_Rlevels(ilev)%rspaceTimeDiscr

              IF (rsolverNode%ioutputLevel .GE. 3) THEN
                CALL output_line ('Multigrid: Current mesh level: '//TRIM(sys_siL(ilev,5)))
              END IF
              
              ! If we are not on the lowest level, repeat the smoothing of 
              ! the solution/restriction of the new defect in the next loop 
              ! pass...
            END DO   ! ilev > minimum level
            
            ! Now we reached the coarse grid.
            
            IF (rsolverNode%ioutputLevel .GE. 3) THEN
              CALL output_line ('Multigrid: Invoking coarse grid solver.')
            END IF
            
            ! Solve the system on lowest level by preconditioning
            ! of the RHS=defect vector.
            CALL sptils_precondDefect(p_rsubnode%p_Rlevels(ilev)%p_rcoarseGridSolver,&
                                      p_rsubnode%p_Rlevels(ilev)%rsolutionVector)
            
            ! Now we have the solution vector on the lowest level - we have to go
            ! upwards now... but probably not to NLMAX! That depends on the cycle.
            !
            DO WHILE (ilev .LT. nlmax)
              ilev = ilev + 1
              p_rmatrix => p_rsubnode%p_Rlevels(ilev)%rspaceTimeDiscr
              
              IF (rsolverNode%ioutputLevel .GE. 3) THEN
                CALL output_line ('Multigrid: Current mesh level: '//TRIM(sys_siL(ilev,5)))
              END IF

              ! Prolongate the solution vector from the coarser level
              ! to the temp vector on the finer level.
              CALL sptipr_performProlongation (&
                    p_rsubnode%p_Rlevels(ilev)%rinterlevelProjection,&
                    p_rsubnode%p_Rlevels(ilev-1)%rsolutionVector, &
                    p_rsubnode%p_Rlevels(ilev)%rtempVector, &
                    p_rsubnode%p_Rlevels(ilev-1)%rprjVector, &
                    p_rsubnode%p_Rlevels(ilev)%rprjVector)

              ! Step length control. By default, choose step length 1.0.
              dstep = 1.0_DP
              
              ! Perform the coarse grid correction by adding the coarse grid
              ! solution (with the calculated step-length parameter) to
              ! the current solution.
              
              CALL sptivec_vectorLinearComb (p_rsubnode%p_Rlevels(ilev)%rtempVector,&
                 p_rsubnode%p_Rlevels(ilev)%rsolutionVector,&
                 dstep,1.0_DP)
                            
              ! Extended output
              IF ((rsolverNode%ioutputLevel .GE. 3) .AND. &
                  (MOD(ite,niteResOutput).EQ.0)) THEN
                  
                CALL sptivec_copyVector (p_rsubnode%p_Rlevels(ilev)%rrhsVector,&
                                         p_rsubnode%p_Rlevels(ilev)%rtempVector)
                CALL c2d2_spaceTimeMatVec (rsolverNode%p_rproblem, p_rmatrix, &
                    p_rsubnode%p_Rlevels(ilev)%rsolutionVector,&
                    p_rsubnode%p_Rlevels(ilev)%rtempVector, -1.0_DP,1.0_DP,dres)

                IF (.NOT.((dres .GE. 1E-99_DP) .AND. &
                          (dres .LE. 1E99_DP))) dres = 0.0_DP
                          
                CALL output_line ('Multigrid: Level '//TRIM(sys_siL(ilev,5))//&
                    ' after c.g.corr.:  !!RES!! = '//TRIM(sys_sdEL(dres,15)) )
              END IF
                                            
              ! Perform the post-smoothing with the current solution vector
              IF (ASSOCIATED(p_rsubnode%p_Rlevels(ilev)%p_rpostSmoother)) THEN
                CALL sptils_smoothCorrection (&
                          p_rsubnode%p_Rlevels(ilev)%p_rpostSmoother,&
                          p_rsubnode%p_Rlevels(ilev)%rsolutionVector,&
                          p_rsubnode%p_Rlevels(ilev)%rrhsVector,&
                          p_rsubnode%p_Rlevels(ilev)%rtempVector)

                ! Extended output
                IF ((rsolverNode%ioutputLevel .GE. 3) .AND. &
                    (MOD(ite,niteResOutput).EQ.0)) THEN
                    
                CALL sptivec_copyVector (p_rsubnode%p_Rlevels(ilev)%rrhsVector,&
                                         p_rsubnode%p_Rlevels(ilev)%rtempVector)
                CALL c2d2_spaceTimeMatVec (rsolverNode%p_rproblem, p_rmatrix, &
                    p_rsubnode%p_Rlevels(ilev)%rsolutionVector,&
                    p_rsubnode%p_Rlevels(ilev)%rtempVector, -1.0_DP,1.0_DP,dres)

                  IF (.NOT.((dres .GE. 1E-99_DP) .AND. &
                            (dres .LE. 1E99_DP))) dres = 0.0_DP
                            
                  CALL output_line ('Multigrid: Level '//TRIM(sys_siL(ilev,5))//&
                      ' after postsm.:    !!RES!! = '//TRIM(sys_sdEL(dres,15)) )
                END IF
                                              
              END IF

              ! Update the iteration counter(s) for realising the MG-cycle(s).
              ! Then either repeat this loop to perform the next prolongation or   
              ! repeat the cycleloop to do perform a next MG sweep on the current      
              ! level.                                                        
              !                                                               
              ! Here icycle defines how the cycle-counters are updated.           
              ! For a W-cycle the cycle counter is resetted to 2 if the sweep is 
              ! fulfilled on the current level, for F-cycle it's set to 1 to not
              ! perform more that 1 cycle on the current level anymore.       
              
              p_rsubnode%p_Rlevels(ilev)%ncyclesRemaining = &
                  p_rsubnode%p_Rlevels(ilev)%ncyclesRemaining-1
              IF (p_rsubnode%p_Rlevels(ilev)%ncyclesRemaining .LE. 0) THEN
                IF (p_rsubnode%icycle .EQ. 0) THEN
                  p_rsubnode%p_Rlevels(ilev)%ncyclesRemaining = 1
                ELSE
                  ! Cycle finished. Reset counter for next cycle.
                  p_rsubnode%p_Rlevels(ilev)%ncyclesRemaining = &
                      p_rsubnode%p_Rlevels(ilev)%ncycles
                END IF
              ELSE

                IF (rsolverNode%ioutputLevel .GE. 3) THEN
                  CALL output_line ('Multigrid: Cycle on level '&
                      //TRIM(sys_siL(ilev,5))//' finished.')
                END IF

                ! Next cycle; go down starting from the current level
                CYCLE cycleloop

              END IF
              
            END DO ! ilev < nlmax
            
            ! We finally reached the maximum level. Quit the cycle loop
            ! to go on with the next iteration.
            EXIT cycleloop
          
          END DO cycleloop
          
          ! We have (hopefully) successfully performed one MG-sweep, starting
          ! and ending on the finest level. As we are now on the finest level
          ! again, we can update our defect vector to test the current
          ! residuum...
          !
          ! Calculate the residuum and its norm.
          CALL sptivec_copyVector (p_rsubnode%p_Rlevels(ilev)%rrhsVector,&
                                    p_rsubnode%p_Rlevels(ilev)%rtempVector)
          CALL c2d2_spaceTimeMatVec (rsolverNode%p_rproblem, p_rmatrix, &
              p_rsubnode%p_Rlevels(ilev)%rsolutionVector,&
              p_rsubnode%p_Rlevels(ilev)%rtempVector, -1.0_DP,1.0_DP,dres)

          IF (.NOT.((dres .GE. 1E-99_DP) .AND. &
                    (dres .LE. 1E99_DP))) dres = 0.0_DP
          
          rsolverNode%dfinalDefect = dres
          
          ! Test if the iteration is diverged
          IF (rsolverNode%depsRel .NE. SYS_INFINITY) THEN
            IF ( .NOT. (dres .LE. rsolverNode%dinitialDefect*rsolverNode%ddivRel) ) THEN
              CALL output_line ('Multigrid: Solution diverging!')
              rsolverNode%iresult = 1
              EXIT
            END IF
          END IF

          ! At least perform nminIterations iterations
          IF (ite .GE. nminIterations) THEN
          
            ! Check if the iteration converged
            !
            ! Absolute convergence criterion? Check the norm directly.
            IF (rsolverNode%depsAbs .NE. 0.0_DP) THEN
              IF (.NOT. (dres .GT. rsolverNode%depsAbs)) THEN
                EXIT
              END IF
            END IF
            
            ! Relative convergence criterion? Multiply with initial residuum
            ! and check the norm. 
            IF (rsolverNode%depsRel .NE. 0.0_DP) THEN
              IF (.NOT. (dres .GT. rsolverNode%depsRel * rsolverNode%dinitialDefect)) THEN
                EXIT
              END IF
            END IF
            
          END IF

          ! print out the current residuum

          IF ((rsolverNode%ioutputLevel .GE. 2) .AND. &
              (MOD(ite,niteResOutput).EQ.0)) THEN
            CALL output_line ('Multigrid: Iteration '// &
                TRIM(sys_siL(ITE,10))//',  !!RES!! = '//&
                TRIM(sys_sdEL(rsolverNode%dfinalDefect,15)) )
          END IF
          
        END DO  ! ite
        
        ! Set ITE to NIT to prevent printing of "NIT+1" of the loop was
        ! completed

        IF (ite .GT. rsolverNode%nmaxIterations) &
          ite = rsolverNode%nmaxIterations

        ! Finish - either with an error or if converged.
        ! Print the last residuum.

        IF ((rsolverNode%ioutputLevel .GE. 2) .AND. &
            (ite .GE. 1) .AND. (ITE .LT. rsolverNode%nmaxIterations) .AND. &
            (rsolverNode%iresult .GE. 0)) THEN
          CALL output_line ('Multigrid: Iteration '// &
              TRIM(sys_siL(ITE,10))//',  !!RES!! = '//&
              TRIM(sys_sdEL(rsolverNode%dfinalDefect,15)) )
        END IF
        
        ! Final number of iterations
        rsolverNode%iiterations = ite
      
      END IF

      ! Finally, we look at some statistics:
      !
      ! Don't calculate anything if the final residuum is out of bounds -
      ! would result in NaN's,...
        
      IF (rsolverNode%dfinalDefect .LT. 1E99_DP) THEN
      
        ! If the initial defect was zero, the solver immediately
        ! exits - and so the final residuum is zero and we performed
        ! no steps; so the resulting multigrid convergence rate stays zero.
        ! In the other case the multigrid convergence rate computes as
        ! (final defect/initial defect) ** 1/nit :

        IF (rsolverNode%dinitialDefect .GT. rsolverNode%drhsZero) THEN
          rsolverNode%dconvergenceRate = &
                      (rsolverNode%dfinalDefect / rsolverNode%dinitialDefect) ** &
                      (1.0_DP/REAL(rsolverNode%iiterations,DP))
        END IF
      END IF
    
    END IF

    ! As the solution vector on the finest level shared its memory with rd,
    ! we just calculated the new correction vector!
      
    IF (rsolverNode%dfinalDefect .LT. 1E99_DP) THEN
      
      IF (rsolverNode%ioutputLevel .GE. 2) THEN
        CALL output_lbrk()
        CALL output_line ('Multigrid statistics:')
        CALL output_lbrk()
        CALL output_line ('Iterations              : '//&
             TRIM(sys_siL(rsolverNode%iiterations,10)) )
        CALL output_line ('!!INITIAL RES!!         : '//&
             TRIM(sys_sdEL(rsolverNode%dinitialDefect,15)) )
        CALL output_line ('!!RES!!                 : '//&
             TRIM(sys_sdEL(rsolverNode%dfinalDefect,15)) )
        IF (rsolverNode%dinitialDefect .GT. rsolverNode%drhsZero) THEN     
          CALL output_line ('!!RES!!/!!INITIAL RES!! : '//&
            TRIM(sys_sdEL(rsolverNode%dfinalDefect / rsolverNode%dinitialDefect,15)) )
        ELSE
          CALL output_line ('!!RES!!/!!INITIAL RES!! : '//&
               TRIM(sys_sdEL(0.0_DP,15)) )
        END IF
        CALL output_lbrk ()
        CALL output_line ('Rate of convergence     : '//&
             TRIM(sys_sdEL(rsolverNode%dconvergenceRate,15)) )

      END IF

      IF (rsolverNode%ioutputLevel .EQ. 1) THEN
        CALL output_line (&
              'MG: Iterations/Rate of convergence: '//&
              TRIM(sys_siL(rsolverNode%iiterations,10))//' /'//&
              TRIM(sys_sdEL(rsolverNode%dconvergenceRate,15)) )
      END IF

    ELSE
      ! DEF=Infinity; RHO=Infinity, set to 1
      rsolverNode%dconvergenceRate = 1.0_DP
      rsolverNode%dasymptoticConvergenceRate = 1.0_DP
    END IF  
  
  END SUBROUTINE

END MODULE
