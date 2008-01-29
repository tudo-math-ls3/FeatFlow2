!##############################################################################
!# ****************************************************************************
!# <name> cc2dmediumm2spacetimesolver </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module realises a linear preconditioner for space-time coupled
!# systems and is somehow an extension to the usual preconditioner for
!# linear systems.
!#
!# The preconditioner module has a similar structure and can similarly be used
!# as the usual linearsolver module. The t_sptilsNode structure is the
!# equvalient to the t_solverNode structure in linearsolver.f90.
!# To initialise / use the space-time coupled preconditioner, proceed as
!# follows:
!#
!#   ! Create a nullified pointer to a t_sptilsNode structure
!#
!#   TYPE(t_sptilsNode), POINTER :: p_rsolver => NULL()
!#
!#   ! Declare an array for the space-time coupled system matrices
!#
!#   TYPE(t_ccoptSpaceTimeDiscretisation), DIMENSION(NLMIN:NLMAX) :: Rmatrices
!#
!#   ! Initialise the preconditioner. This produces you a pointer to a
!#   ! t_sptilsNode structure that identifies the solver.
!#
!#   CALL sptils_initXXXX (p_rsolver,...)
!#
!#   ! Attach the system matrices (a set of t_ccoptSpaceTimeDiscretisation
!#   ! structures) to the solver
!#
!#   ... -> Rmatrices(NLMIN:NLMAX)
!#   CALL sptils_setMatrices (...,Rmatrices,p_rsolver)
!#
!#   ! Call the initialisation routines
!#
!#   CALL sptils_initStructure (p_rsolver)
!#   CALL sptils_initData (p_rsolver)
!#
!#   ! Preform the preconditioning for a given space-time defect vector rd
!#
!#   CALL sptils_precondDefect (p_rsolver,rd)
!#
!#   ! Clean up the structure of the solver
!#
!#   CALL sptils_doneData (p_rsolver)
!#   CALL sptils_doneStructure (p_rsolver)
!#   CALL sptils_releaseSolver (p_rsolver)
!#   
!# The usage of the Multigrid preconditioner is slightly more complicated.
!# You have to initialise a smoother and a coarse-grid solver and attach
!# it to MG. Here an example of how to create a MG preconditioner with
!# Block-Jacobi smoother. The coarse-grid solver is Defect correction
!# with Block-Jacobi preconditioning:
!#
!#   ! Create a nullified pointer to a t_sptilsNode structure.
!#   ! Preconditioner, smoother, coarse grid solver, main solver.
!#
!#   TYPE(t_sptilsNode), POINTER :: p_rsolver => NULL()
!#   TYPE(t_sptilsNode), POINTER :: p_rcgrSolver => NULL()
!#   TYPE(t_sptilsNode), POINTER :: p_rprecond => NULL()
!#   TYPE(t_sptilsNode), POINTER :: p_rSmoother => NULL()
!#
!#   ! Declare an array for the space-time coupled system matrices
!#
!#   TYPE(t_ccoptSpaceTimeDiscretisation), DIMENSION(NLMIN:NLMAX) :: Rmatrices
!#
!#   ! Initialise the MG preconditioner. This produces you a pointer to a
!#   ! t_sptilsNode structure that identifies the solver.
!#   ! TIMENLMIN and TIMENLMAX are the time refinement levels. The solution
!#   ! is calculated on time level TIMENLMAX.
!#
!#   CALL sptils_initMultigrid (...,TIMENLMIN,TIMENLMAX,p_rsolver)
!#
!#   ! Attach level information to MG. On every level, create a separate
!#   ! smoother. On the coarse grid, create a coarse grid solver
!#
!#   DO ilev=TIMENLMIN,TIMENLMAX
!#
!#     NULLIFY(p_rsmoother)
!#     NULLIFY(p_rcgrSolver)
!#
!#     IF (ilev .EQ. 1) THEN
!#       ! On the minimum level, create a coarse grid solver.
!#       !
!#       ! Block Jacobi preconditioner
!#       CALL sptils_initBlockJacobi (rproblem,p_rprecond,RspatialPrecond(ilev))
!#       
!#       ! Defect correction solver
!#       CALL sptils_initDefCorr (rproblem,p_rcgrSolver,p_rprecond)
!#     ELSE
!#       ! On higher levels, create a Block Jacobi (pre- and post-) smoother
!#       CALL sptils_initBlockJacobi (rproblem,p_rsmoother,...)
!#       CALL sptils_convertToSmoother (p_rsmoother...)
!#     END IF
!#      
!#     ! Iinally initialise the level ilev with that       
!#     CALL sptils_setMultigridLevel (p_rsolver,ilev,...,&
!#         p_rsmoother,p_rsmoother,p_rcgrSolver)
!#   END DO
!#
!#   ! Attach the system matrices (a set of t_ccoptSpaceTimeDiscretisation
!#   ! structures) to the solver
!#
!#   ... -> Rmatrices(NLMIN:NLMAX)
!#   CALL sptils_setMatrices (p_rsolver,Rmatrices)
!#
!#   ! Call the initialisation routines
!#
!#   CALL sptils_initStructure (p_rsolver)
!#   CALL sptils_initData (p_rsolver)
!#
!#   ! Preform the preconditioning for a given space-time defect vector rd
!#
!#   CALL sptils_precondDefect (p_rsolver,rd)
!#
!#   ! Clean up the structure of the solver
!#
!#   CALL sptils_doneData (p_rsolver)
!#   CALL sptils_doneStructure (p_rsolver)
!#   CALL sptils_releaseSolver (p_rsolver)
!# 
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
  USE cc2dmediumm2timeanalysis
  USE cc2dmediumm2boundary
  USE cc2dmediumm2discretisation
  USE cc2dmediumm2postprocessing
  USE cc2dmediumm2matvecassembly
  USE cc2dmediumm2nonlinearcore
  USE spacetimediscretisation
  USE spacetimelinearsystem
  
  USE spacetimevectors
  USE spacetimeinterlevelprojection
  USE dofmapping
  USE timeboundaryconditions
  
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
  
  ! Block Gauss-Seidel preconditioner
  INTEGER, PARAMETER :: SPTILS_ALG_BLOCKSOR      = 4
  
  ! CG iteration (preconditioned) 
  INTEGER, PARAMETER :: SPTILS_ALG_CG            = 6

  ! BiCGStab iteration (preconditioned) 
  INTEGER, PARAMETER :: SPTILS_ALG_BICGSTAB      = 7

  ! UMFPACK solver
  INTEGER, PARAMETER :: SPTILS_ALG_UMFPACK4      = 11

  ! Block orward-Backward Gauss-Seidel solver
  INTEGER, PARAMETER :: SPTILS_ALG_BLOCKFBGS    = 12

  ! Multigrid iteration
  INTEGER, PARAMETER :: SPTILS_ALG_MULTIGRID     = 9

!</constantblock>

! *****************************************************************************

!<constantblock description="Error constants returned by initialisation routines">

  ! Initialisation routine went fine
  INTEGER, PARAMETER :: SPTILS_ERR_NOERROR       = 0
  
!</constantblock>

! *****************************************************************************

!<constantblock description="Bitfield identifiers for the ability of a solver">

  ! Solver can handle multiple levels
  INTEGER(I32), PARAMETER :: SPTILS_ABIL_MULTILEVEL   = 2**2
  
  ! Solver allows checking the defect during the iteration.
  ! Solvers not capable of this perform only a fixed number of solution
  ! steps (e.g. UMFPACK performs always one step).
  INTEGER(I32), PARAMETER :: SPTILS_ABIL_CHECKDEF     = 2**3
  
  ! Solver is a direct solver (e.g. UMFPACK, ILU).
  ! Otherwise the solver is of iterative nature and might perform
  ! multiple steps to solve the problem.
  INTEGER(I32), PARAMETER :: SPTILS_ABIL_DIRECT       = 2**4
  
  ! Solver might use subsolvers (preconditioners, smoothers,...)
  INTEGER(I32), PARAMETER :: SPTILS_ABIL_USESUBSOLVER = 2**5
  
!</constantblock>

! *****************************************************************************

!<constantblock description="Identifiers for stopping criterium istoppingCriterion.">

  ! Use standard stopping criterion.
  ! If depsRel>0: use relative stopping criterion.
  ! If depsAbs>0: use abs stopping criterion.
  ! If both are > 0: use both, i.e. the iteration stops when both,
  !    the relative AND the absolute stopping criterium holds
  INTEGER, PARAMETER :: SPTILS_STOP_STANDARD     = 0

  ! Use 'minimum' stopping criterion.
  ! If depsRel>0: use relative stopping criterion.
  ! If depsAbs>0: use abs stopping criterion.
  ! If both are > 0: use one of them, i.e. the iteration stops when the
  !    either the relative OR the absolute stopping criterium holds
  INTEGER, PARAMETER :: SPTILS_STOP_ONEOF        = 1
  
!</constantblock>

!</constants>

! *****************************************************************************

!<types>

!<typeblock>
  
  ! The following structure defines a 'solver node' for the 
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
    ! SPTILS_STOP_xxxx constants.
    INTEGER                    :: istoppingCriterion = SPTILS_STOP_STANDARD

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
    
    ! READ ONLY: Solver ability tag. 
    ! Bitfield. A combination of SPTILS_ABIL_xxxx
    ! flags that specify the ability of the actual solver (e.g.
    ! whether it can handle block-matrices or only scalar matrices,...).
    ! Calling a solver that is not able to handle a specific problem
    ! will cause an error by the framework.
    INTEGER(I32)                    :: ccapability = 0

    ! A pointer to the problem structure that defines the problem.
    TYPE(t_problem), POINTER   :: p_rproblem
    
    ! STATUS FOR ITERATIVE SOLVERS: Current iteration
    INTEGER                    :: icurrentIteration

    ! Structure of the space time matrix that the solver should use.
    TYPE(t_ccoptSpaceTimeMatrix) :: rmatrix
    
    ! Pointer to a structure for the Defect correction solver; NULL() if not set
    TYPE (t_sptilsSubnodeDefCorr), POINTER        :: p_rsubnodeDefCorr     => NULL()

    ! Pointer to a structure for the CG solver; NULL() if not set
    TYPE (t_sptilsSubnodeCG), POINTER             :: p_rsubnodeCG     => NULL()

    ! Pointer to a structure for the BiCGStab solver; NULL() if not set
    TYPE (t_sptilsSubnodeBiCGStab), POINTER       :: p_rsubnodeBiCGStab => NULL()

    ! Pointer to a structure for the block Jacobi preconditioner
    TYPE(t_sptilsSubnodeBlockJacobi), POINTER     :: p_rsubnodeBlockJacobi => NULL()

    ! Pointer to a structure for the block Gauss-Seidel preconditioner
    TYPE(t_sptilsSubnodeBlockSOR), POINTER        :: p_rsubnodeBlockSOR => NULL()

    ! Pointer to a structure for the multigrid preconditioner
    TYPE(t_sptilsSubnodeMultigrid), POINTER       :: p_rsubnodeMultigrid => NULL()

    ! Pointer to a structure for the UMFPACK4 preconditioner
    TYPE(t_sptilsSubnodeUMFPACK4), POINTER        :: p_rsubnodeUMFPACK4 => NULL()

    ! Pointer to a structure for the UMFPACK4 preconditioner
    TYPE(t_sptilsSubnodeBlockFBGS), POINTER       :: p_rsubnodeBlockFBGS => NULL()

  END TYPE
  
!</typeblock>

! *****************************************************************************

!<typeblock>
  
  ! This structure realises the subnode for the Defect correction solver.
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
    
    ! A temp vector in space.
    TYPE(t_vectorBlock)               :: rtempVectorSpace
  END TYPE
  
!</typeblock>

! *****************************************************************************

!<typeblock>
  
  ! This structure realises the subnode for the CG solver.
  ! The entry p_rpreconditioner points either to NULL() or to another
  ! t_linsolNode structure for the solver that realises the 
  ! preconditioning.
  
  TYPE t_sptilsSubnodeCG
  
    ! Temporary vectors to use during the solution process
    TYPE(t_spacetimeVector), DIMENSION(4) :: RtempVectors

    ! A pointer to the solver node for the preconditioner or NULL(),
    ! if no preconditioner is used.
    TYPE(t_sptilsNode), POINTER       :: p_rpreconditioner            => NULL()
  
    ! A temporary space-time vector.
    TYPE(t_vectorBlock)               :: rtempVectorSpace
  
  END TYPE
  
!</typeblock>

! *****************************************************************************

!<typeblock>
  
  ! This structure realises the subnode for the BiCGStab solver.
  ! The entry p_rpreconditioner points either to NULL() or to another
  ! t_linsolNode structure for the solver that realises the 
  ! preconditioning.
  
  TYPE t_sptilsSubnodeBiCGStab
  
    ! Temporary vectors to use during the solution process
    TYPE(t_spaceTimeVector), DIMENSION(6) :: RtempVectors
    
    ! Another temporary vector for right- and symmetrical preconditioning
    TYPE(t_vectorBlock) :: rprecondTemp
    
    ! A pointer to the solver node for the preconditioner or NULL(),
    ! if no preconditioner is used.
    TYPE(t_sptilsNode), POINTER       :: p_rpreconditioner            => NULL()

    ! A temporary space-time vector.
    TYPE(t_vectorBlock)               :: rtempVectorSpace
    
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
  
  ! This structure realises the subnode for the Block Gauss-Seidel solver.
  
  TYPE t_sptilsSubnodeBlockSOR

    ! Pointer to a spatial preconditioner structure that defines the
    ! preconditioning in each substep.
    TYPE(t_ccspatialPreconditioner), POINTER :: p_rspatialPreconditioner

  END TYPE
  
!</typeblock>

! *****************************************************************************

!<typeblock>
  
  ! This structure realises the subnode for the Block Gauss-Seidel solver.
  
  TYPE t_sptilsSubnodeBlockFBGS

    ! Pointer to a spatial preconditioner structure that defines the
    ! preconditioning in each substep.
    TYPE(t_ccspatialPreconditioner), POINTER :: p_rspatialPreconditioner

    ! A temporary space-time vector.
    TYPE(t_spacetimeVector)           :: rtempVector
    
    ! SOR-Damping parameter for damping the time dependence
    REAL(DP) :: domegaSOR

  END TYPE
  
!</typeblock>

! *****************************************************************************

!<typeblock>
  
  ! Level data for multigrid. This structure forms an entry in a linked list
  ! of levels which multigrid uses for solving the given problem.
  ! The solver subnode t_sptilsSubnodeMultigrid holds pointers to the head
  ! and the tail of the list.
  
  TYPE t_sptilsMGLevelInfo
  
    ! Structure that defines the space time matrix on that level.
    TYPE(t_ccoptSpaceTimeMatrix)        :: rmatrix
    
    ! A RHS vector for that level. The memory for this is allocated
    ! in initStructure and released in doneStructure.
    TYPE(t_spacetimeVector)             :: rrhsVector

    ! A solution vector for that level. The memory for this is allocated
    ! in initStructure and released in doneStructure.
    TYPE(t_spacetimeVector)             :: rsolutionVector

    ! A temporary vector for the adaptive coarse grid correction
    TYPE(t_spacetimeVector)             :: rtempCGCvector

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
    
    ! Pointer to a structure for the CG solver; NULL() if not set
    TYPE (t_sptilsSubnodeCG), POINTER   :: p_rsubnodeCG          => NULL()

    ! An interlevel projection structure that configures the transport
    ! of defect/solution vectors from one level to another.
    ! For the coarse grid, this structure is ignored.
    ! For a finer grid (e.g. level 4), this defines the grid transfer
    ! between the current and the lower level (so here between level 3 and
    ! level 4).
    TYPE(t_sptiProjection)   :: rinterlevelProjection
    
    ! Relative daptive cycle convergence criterion for coarse levels.
    ! This value is usually =1E99_DP which deactivates adaptive cycles.
    ! The user can set this variable on initialisation of Multigrid to
    ! a value < 1E99_DP. In this case, the complete multigrid cycle on all 
    ! levels except for the fine grid is repeated until
    !  |res. after postsmoothing| < depsRelCycle * |initial res on that level|.
    ! This allows 'adaptive cycles' which e.g. gain one digit on a coarse
    ! level before prolongating the solution to the fine grid.
    ! This is an extension to the usual F/V/W-cycle scheme.
    REAL(DP)                      :: depsRelCycle             = 1E99_DP

    ! Absolute adaptive cycle convergence criterion for coarse levels.
    ! This value is usually =1E99_DP which deactivates adaptive cycles.
    ! The user can set this variable on initialisation of Multigrid to
    ! a value < 1E99_DP. In this case, the complete multigrid cycle on all 
    ! levels except for the fine grid is repeated until
    !  |res. after postsmoothing| < depsAbsCycle.
    ! This allows 'adaptive cycles' which e.g. gain one digit on a coarse
    ! level before prolongating the solution to the fine grid.
    ! This is an extension to the usual F/V/W-cycle scheme.
    REAL(DP)                      :: depsAbsCycle             = 1E99_DP
    
    ! Maximum number of adaptive cycles performed on that level.
    ! Only used if adaptive cycles are activated by setting 
    ! depsRelCycle < 1E99_DP. -1=infinity=standard
    INTEGER                       :: nmaxAdaptiveCycles       = -1
    
    ! STATUS/INTERNAL: A temporary vector used for prolongation/restriction
    TYPE(t_vectorBlock)      :: rprjVector
    
    ! STATUS/INTERNAL: MG cycle information.
    ! Number of cycles to perform on this level.
    INTEGER                        :: ncycles

    ! STATUS/INTERNAL: MG cycle information.
    ! Number of remaining cycles to perform on this level.
    INTEGER                        :: ncyclesRemaining
    
    ! STATUS/INTERNAL: initial residuum when a solution is restricted to this
    ! level. Only used if adaptive cycles are activated by setting 
    ! depsRelCycle < 1E99_DP in t_linsolSubnodeMultigrid.
    REAL(DP)                       :: dinitResCycle = 0.0_DP
    
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
    
    ! INPUT PARAMETER: Minimum level in p_Rlevels.
    INTEGER                       :: NLMIN                    = 1
    
    ! INPUT PARAMETER: Maximum level in p_Rlevels.
    INTEGER                       :: NLMAX                    = 1

    ! INPUT PARAMETER: Minimum value for adaptive coarse grid correction.
    ! A value of dalphamin=dalphamax deactivates the adaptive coarse
    ! grid correction. Standard = 1.0.
    REAL(DP)                      :: dalphamin                = 1.0_DP
    
    ! INPUT PARAMETER: Maximum value for adaptive coarse grid correction.
    ! A value of dalphamin=dalphamax deactivates the adaptive coarse
    ! grid correction. Standard = 1.0.
    REAL(DP)                      :: dalphamax                = 1.0_DP
    
    ! Array of t_sptilsMGLevelInfo structures for all the levels the MG solver
    ! should handle.
    TYPE(t_sptilsMGLevelInfo), DIMENSION(:), POINTER :: p_Rlevels => NULL()
    
  END TYPE
  
!</typeblock>

! *****************************************************************************

!<typeblock>
  
  ! This structure realises the subnode for the UMFPACK4 solver.
  
  TYPE t_sptilsSubnodeUMFPACK4
  
    ! Control structure for UMFPACK4; contains parameter for the solver
    REAL(DP), DIMENSION(20) :: Dcontrol

    ! Handle for symbolic factorisation.
    ! This is not a FEAT-Handle!
    INTEGER(I32) :: isymbolic = 0

    ! Handle for numeric factorisation
    ! This is not a FEAT-Handle!
    INTEGER(I32) :: inumeric = 0

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
  TYPE(t_ccoptSpaceTimeMatrix), DIMENSION(:), INTENT(IN) :: Rmatrices
!</input>
  
!<inputoutput>
  ! The solver node which should be initialised
  TYPE(t_sptilsNode), INTENT(INOUT)             :: rsolverNode
!</inputoutput>
  
!</subroutine>

  ! Copy the matrix structure on the finest level to the rsolverNode
  ! structure. This corresponds to the system we want to solve.
  rsolverNode%rmatrix = Rmatrices(UBOUND(Rmatrices,1))
  
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
  CASE (SPTILS_ALG_CG)
    CALL sptils_setMatrixCG (rsolverNode,Rmatrices)
  CASE (SPTILS_ALG_BiCGStab)
    CALL sptils_setMatrixBiCGStab (rsolverNode,Rmatrices)
  CASE (SPTILS_ALG_BLOCKJACOBI)
    CALL sptils_setMatrixBlockJacobi (rsolverNode,Rmatrices)
  CASE (SPTILS_ALG_BLOCKSOR)
    CALL sptils_setMatrixBlockSOR (rsolverNode,Rmatrices)
  CASE (SPTILS_ALG_BLOCKFBGS)
    CALL sptils_setMatrixBlockFBGS (rsolverNode,Rmatrices)
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
  ! One of the SPTILS_ERR_XXXX constants. A value different to 
  ! SPTILS_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  INTEGER, INTENT(OUT) :: ierror
!</output>
  
!</subroutine>

    ! A-priori we have no error...
    ierror = SPTILS_ERR_NOERROR
    
    ! Call the structure-init routine of the specific solver
    
    SELECT CASE(rsolverNode%calgorithm)
    CASE (SPTILS_ALG_DEFCORR)
      CALL sptils_initStructureDefCorr (rsolverNode,ierror)
    CASE (SPTILS_ALG_CG)
      CALL sptils_initStructureCG (rsolverNode,ierror)
    CASE (SPTILS_ALG_BICGSTAB)
      CALL sptils_initStructureBiCGStab (rsolverNode,ierror)
    CASE (SPTILS_ALG_BLOCKJACOBI)
      CALL sptils_initStructureBlockJacobi (rsolverNode,ierror)
    CASE (SPTILS_ALG_BLOCKSOR)
      CALL sptils_initStructureBlockSOR (rsolverNode,ierror)
    CASE (SPTILS_ALG_BLOCKFBGS)
      CALL sptils_initStructureBlockFBGS (rsolverNode,ierror)
    CASE (SPTILS_ALG_MULTIGRID)
      CALL sptils_initStructureMultigrid (rsolverNode,ierror)
    CASE (SPTILS_ALG_UMFPACK4)
      CALL sptils_initStructureUMFPACK4 (rsolverNode,ierror)
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
  ! One of the SPTILS_ERR_XXXX constants. A value different to 
  ! SPTILS_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  INTEGER, INTENT(OUT) :: ierror
!</output>

!<inputoutput>
  ! The solver node which should be initialised
  TYPE(t_sptilsNode), INTENT(INOUT)                     :: rsolverNode
!</inputoutput>

!</subroutine>

    ! A-priori we have no error...
    ierror = SPTILS_ERR_NOERROR

    ! Call the data-init routine of the specific solver
    
    SELECT CASE(rsolverNode%calgorithm)
    CASE (SPTILS_ALG_DEFCORR)
      CALL sptils_initDataDefCorr (rsolverNode,ierror)
    CASE (SPTILS_ALG_CG)
      CALL sptils_initDataCG (rsolverNode,ierror)
    CASE (SPTILS_ALG_BICGSTAB)
      CALL sptils_initDataBiCGStab (rsolverNode,ierror)
    CASE (SPTILS_ALG_BLOCKJACOBI)
      CALL sptils_initDataBlockJacobi (rsolverNode,ierror)
    CASE (SPTILS_ALG_BLOCKSOR)
      CALL sptils_initDataBlockSOR (rsolverNode,ierror)
    CASE (SPTILS_ALG_BLOCKFBGS)
      CALL sptils_initDataBlockFBGS (rsolverNode,ierror)
    CASE (SPTILS_ALG_MULTIGRID)
      CALL sptils_initDataMultigrid (rsolverNode,ierror)
    CASE (SPTILS_ALG_UMFPACK4)
      CALL sptils_initDataUMFPACK4 (rsolverNode,ierror)
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
  ! One of the SPTILS_ERR_XXXX constants. A value different to 
  ! SPTILS_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  INTEGER, INTENT(OUT) :: ierror
!</output>
  
!</subroutine>

    ! A-priori we have no error...
    ierror = SPTILS_ERR_NOERROR

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
  ! One of the SPTILS_ERR_XXXX constants. A value different to 
  ! SPTILS_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  INTEGER, INTENT(OUT) :: ierror
!</output>

!<inputoutput>
  ! The solver node containing the solver confuguration
  TYPE(t_sptilsNode), INTENT(INOUT)                     :: rsolverNode
!</inputoutput>
  
!</subroutine>

    ! A-priori we have no error...
    ierror = SPTILS_ERR_NOERROR

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
    CASE (SPTILS_ALG_DEFCORR)
      CALL sptils_doneDataDefCorr (rsolverNode)
    CASE (SPTILS_ALG_CG)
      CALL sptils_doneDataCG (rsolverNode)
    CASE (SPTILS_ALG_BICGSTAB)
      CALL sptils_doneDataBiCGStab (rsolverNode)
    CASE (SPTILS_ALG_BLOCKJACOBI)
      CALL sptils_doneDataBlockJacobi (rsolverNode)
    CASE (SPTILS_ALG_BLOCKSOR)
      CALL sptils_doneDataBlockSOR (rsolverNode)
    CASE (SPTILS_ALG_BLOCKFBGS)
      CALL sptils_doneDataBlockFBGS (rsolverNode)
    CASE (SPTILS_ALG_MULTIGRID)
      CALL sptils_doneDataMultigrid (rsolverNode)
    CASE (SPTILS_ALG_UMFPACK4)
      CALL sptils_doneDataUMFPACK4 (rsolverNode)
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
    CASE (SPTILS_ALG_DEFCORR)
      CALL sptils_doneStructureDefCorr (rsolverNode)
    CASE (SPTILS_ALG_CG)
      CALL sptils_doneStructureCG (rsolverNode)
    CASE (SPTILS_ALG_BICGSTAB)
      CALL sptils_doneStructureBiCGStab (rsolverNode)
    CASE (SPTILS_ALG_BLOCKJACOBI)
      CALL sptils_doneStructureBlockJacobi (rsolverNode)
    CASE (SPTILS_ALG_BLOCKSOR)
      CALL sptils_doneStructureBlockSOR (rsolverNode)
    CASE (SPTILS_ALG_BLOCKFBGS)
      CALL sptils_doneStructureBlockFBGS (rsolverNode)
    CASE (SPTILS_ALG_MULTIGRID)
      CALL sptils_doneStructureMultigrid (rsolverNode)
    CASE (SPTILS_ALG_UMFPACK4)
      CALL sptils_doneStructureUMFPACK4 (rsolverNode)
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
    CASE (SPTILS_ALG_DEFCORR)
      CALL sptils_doneDefCorr (p_rsolverNode)
    CASE (SPTILS_ALG_CG)
      CALL sptils_doneCG (p_rsolverNode)
    CASE (SPTILS_ALG_BICGSTAB)
      CALL sptils_doneBiCGStab (p_rsolverNode)
    CASE (SPTILS_ALG_BLOCKJACOBI)
      CALL sptils_doneBlockJacobi (p_rsolverNode)
    CASE (SPTILS_ALG_BLOCKSOR)
      CALL sptils_doneBlockSOR (p_rsolverNode)
    CASE (SPTILS_ALG_BLOCKFBGS)
      CALL sptils_doneBlockFBGS (p_rsolverNode)
    CASE (SPTILS_ALG_MULTIGRID)
      CALL sptils_doneMultigrid (p_rsolverNode)
    CASE (SPTILS_ALG_UMFPACK4)
      CALL sptils_doneUMFPACK4 (p_rsolverNode)
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

!<function>
  
  LOGICAL FUNCTION sptils_testConvergence (rsolverNode, dvecNorm, rdef) RESULT(loutput)
  
!<description>
  ! Tests a defect vector rdef whether it is in a defined tolerance configured 
  ! in the solver node, so the iteration of an iterative solver can be
  ! can be treated as 'converged'.
  ! The iteration is treated as 'converged' if both, the relative and the
  ! absolute convergence criterion are fulfilled (or only one of them,
  ! respectively, if the other is switched off).
  !
  ! The solver must have initialised the 'dinitialDefect' variable
  ! of the solver structure for this routine to work properly!
!</description>
  
!<result>
  ! Boolean value. =TRUE if the convergence criterion is reached; 
  ! =FALSE otherwise.
!</result>
  
!<input>
  ! The solver node that contains the convergence criterion
  TYPE(t_sptilsNode), INTENT(IN) :: rsolverNode
  
  ! OPTIONAL: The defect vector which norm should be tested.
  ! If existent, the norm of the vector is returned in dvecNorm.
  ! If not existent, the routine assumes that dvecNrm is the norm
  ! of the vector and checks convergence depending on dvecNorm.
  TYPE(t_spaceTimeVector), INTENT(IN), OPTIONAL :: rdef
!</input>

!<inputoutput>
  ! Norm of the defect vector. 
  ! If rdef if present, the routine will calculate the norm of rdef and return
  ! it in dvecNorm.
  ! If rdef is not present, dvecNorm is assumed to be a valid norm of a
  ! vector and convergence is tested using dvecNorm.
  REAL(DP), INTENT(INOUT) :: dvecNorm
!</inputoutput>
  
!</function>

    ! Calculate the norm of the vector or take the one given
    ! as parameter
    IF (PRESENT(rdef)) THEN
      dvecNorm = sptivec_vectorNorm (rdef,rsolverNode%iresNorm)
    END IF
    
    SELECT CASE (rsolverNode%istoppingCriterion)
    
    CASE (SPTILS_STOP_ONEOF)
      ! Iteration stops if either the absolute or the relative criterium holds.
      loutput = .FALSE.
      
      ! Absolute convergence criterion? Check the norm directly.
      IF (rsolverNode%depsAbs .NE. 0.0_DP) THEN
        IF (.NOT. (dvecNorm .GT. rsolverNode%depsAbs)) THEN
          loutput = .TRUE.
          RETURN
        END IF
      END IF
      
      ! Relative convergence criterion? Multiply with initial residuum
      ! and check the norm. 
      IF (rsolverNode%depsRel .NE. 0.0_DP) THEN
        IF (.NOT. &
            (dvecNorm .GT. rsolverNode%depsRel * rsolverNode%dinitialDefect)) THEN
          loutput = .TRUE.
          RETURN
        END IF
      END IF
    
    CASE DEFAULT
      ! Standard stopping criterion.
      ! Iteration stops if both the absolute and the relative criterium holds.
      loutput = .TRUE.
      
      ! Absolute convergence criterion? Check the norm directly.
      IF (rsolverNode%depsAbs .NE. 0.0_DP) THEN
        IF (dvecNorm .GT. rsolverNode%depsAbs) THEN
          loutput = .FALSE.
          RETURN
        END IF
      END IF
      
      ! Relative convergence criterion? Multiply with initial residuum
      ! and check the norm. 
      IF (rsolverNode%depsRel .NE. 0.0_DP) THEN
        IF (dvecNorm .GT. rsolverNode%depsRel * rsolverNode%dinitialDefect) THEN
          loutput = .FALSE.
          RETURN
        END IF
      END IF
    END SELECT
    
  END FUNCTION
  
  ! ***************************************************************************

!<function>
  
  LOGICAL FUNCTION sptils_testDivergence (rsolverNode, dvecNorm, rdef) RESULT(loutput)
  
!<description>
  ! Tests a defect vector rx whether it is out of a defined tolerance configured 
  ! in the solver node, so the iteration of an iterative solver can be
  ! can be treated as 'diverged'.
  ! The iteration is treated as 'diverged' if one criterion, the relative or the
  ! absolute divergence criterion is fulfilled.
  !
  ! The solver must have initialised the 'dinitialDefect' variable
  ! of the solver structure for this routine to work properly!
!</description>
  
!<result>
  ! Boolean value. =TRUE if the divergence criterion is reached; 
  ! =FALSE otherwise.
!</result>
  
!<input>
  ! The solver node that contains the divergence criterion
  TYPE(t_sptilsNode), INTENT(IN) :: rsolverNode
  
  ! OPTIONAL: The defect vector which norm should be tested.
  ! If existent, the norm of the vector is returned in dvecNorm.
  ! If not existent, the routine assumes that dvecNrm is the norm
  ! of the vector and checks divergence depending on dvecNorm.
  TYPE(t_spaceTimeVector), INTENT(IN), OPTIONAL :: rdef
!</input>

!<inputoutput>
  ! Norm of the defect vector. 
  ! If rdef if present, the routine will calculate the norm of rdef and return
  ! it in dvecNorm.
  ! If rdef is not present, dvecNorm is assumed to be a valid norm of a
  ! vector and divergence is tested using dvecNorm.
  REAL(DP), INTENT(INOUT) :: dvecNorm
!</inputoutput>

!</function>

    ! Calculate the norm of the vector if not given
    ! as parameter
    IF (PRESENT(rdef)) THEN
      dvecNorm = sptivec_vectorNorm (rdef,rsolverNode%iresNorm)
    END IF
    
    loutput = .FALSE.
    
    ! Absolute divergence criterion? Check the norm directly.
    IF (rsolverNode%ddivAbs .NE. SYS_INFINITY) THEN
     
      ! use NOT here - gives a better handling of special cases like NaN!
      IF ( .NOT. (dvecNorm .LE. rsolverNode%ddivAbs)) THEN
        loutput = .TRUE.
      END IF
      
    END IF
    
    ! Relative divergence criterion? Multiply with initial residuum
    ! and check the norm. 
    IF (rsolverNode%depsRel .NE. SYS_INFINITY) THEN
      IF ( .NOT. (dvecNorm .LE. rsolverNode%dinitialDefect*rsolverNode%ddivRel) ) THEN
        loutput = .TRUE.
      END IF
    END IF
  
  END FUNCTION
  
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
    CASE (SPTILS_ALG_CG)
      CALL sptils_precCG (rsolverNode,rd)
    CASE (SPTILS_ALG_BICGSTAB)
      CALL sptils_precBiCGStab (rsolverNode,rd)
    CASE (SPTILS_ALG_BLOCKJACOBI)
      CALL sptils_precBlockJacobi (rsolverNode,rd)
    CASE (SPTILS_ALG_BLOCKSOR)
      CALL sptils_precBlockSOR (rsolverNode,rd)
    CASE (SPTILS_ALG_BLOCKFBGS)
      CALL sptils_precBlockFBGS (rsolverNode,rd)
    CASE (SPTILS_ALG_MULTIGRID)
      CALL sptils_precMultigrid (rsolverNode,rd)
    CASE (SPTILS_ALG_UMFPACK4)
      CALL sptils_precUMFPACK4 (rsolverNode,rd)
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
    
    ! Initialise the ability bitfield with the ability of this solver:
    p_rsolverNode%ccapability = SPTILS_ABIL_CHECKDEF + &
                                SPTILS_ABIL_USESUBSOLVER

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
  ! The t_sptilsNode structure which is to be cleaned up.
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
  TYPE(t_ccoptSpaceTimeMatrix), DIMENSION(:), INTENT(IN) :: Rmatrices
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
  ! The t_sptilsNode structure 
  TYPE(t_sptilsNode), INTENT(INOUT), TARGET :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the SPTILS_ERR_XXXX constants. A value different to 
  ! SPTILS_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  INTEGER, INTENT(OUT) :: ierror
!</output>
  
!</subroutine>

    ! Local variables
    TYPE(t_ccoptSpaceTimeDiscretisation), POINTER :: p_rspaceTimeDiscr
    
    p_rspaceTimeDiscr => rsolverNode%rmatrix%p_rspaceTimeDiscretisation
    
    ! A-priori we have no error...
    ierror = SPTILS_ERR_NOERROR
    
    ! Allocate memory for the two temp vectors.
    CALL sptivec_initVector (rsolverNode%p_rsubnodeDefCorr%rtempVector,&
        p_rspaceTimeDiscr%NEQtime,&
        p_rspaceTimeDiscr%p_rlevelInfo%p_rdiscretisation)

    CALL sptivec_initVector (rsolverNode%p_rsubnodeDefCorr%rtempVector2,&
        p_rspaceTimeDiscr%NEQtime,&
        p_rspaceTimeDiscr%p_rlevelInfo%p_rdiscretisation)
        
    ! and memory for a spatial temp vector.
    CALL lsysbl_createVecBlockByDiscr (&
        p_rspaceTimeDiscr%p_rlevelInfo%p_rdiscretisation,&
        rsolverNode%p_rsubnodeDefCorr%rtempVectorSpace)
    rsolverNode%p_rsubnodeDefCorr%rtempVectorSpace%p_rdiscreteBC => &
        p_rspaceTimeDiscr%p_rlevelInfo%p_rdiscreteBC
    rsolverNode%p_rsubnodeDefCorr%rtempVectorSpace%p_rdiscreteBCfict => &
        p_rspaceTimeDiscr%p_rlevelInfo%p_rdiscreteFBC

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
  ! The t_sptilsNode structure 
  TYPE(t_sptilsNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the SPTILS_ERR_XXXX constants. A value different to 
  ! SPTILS_ERR_NOERROR indicates that an error happened during the
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
  ! The t_sptilsNode structure 
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
  ! The t_sptilsNode structure 
  TYPE(t_sptilsNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!</subroutine>

    ! Release the temp vectors.
    ! Note that the vectors may already be released from a previous call.
    IF (rsolverNode%p_rsubnodeDefCorr%rtempVector2%NEQtime .NE. 0) &
      CALL sptivec_releaseVector(rsolverNode%p_rsubnodeDefCorr%rtempVector2)
    IF (rsolverNode%p_rsubnodeDefCorr%rtempVector%NEQtime .NE. 0) &
      CALL sptivec_releaseVector(rsolverNode%p_rsubnodeDefCorr%rtempVector)
    IF (rsolverNode%p_rsubnodeDefCorr%rtempVectorSpace%NEQ .NE. 0) &
      CALL lsysbl_releaseVector(rsolverNode%p_rsubnodeDefCorr%rtempVectorSpace)

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
  REAL(DP) :: dres,dfr,dresnorm
  TYPE(t_spacetimeVector), POINTER :: p_rx,p_rdef
  TYPE(t_sptilsNode), POINTER :: p_rprecSubnode
  
  ! Damping parameter
  REAL(DP) :: domega

  ! The system matrix
  TYPE(t_ccoptSpaceTimeMatrix), POINTER :: p_rmatrix
  
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
    p_rmatrix => rsolverNode%rmatrix

    ! Check the parameters
    IF (rd%NEQtime .EQ. 0) THEN
    
      ! Parameters wrong
      rsolverNode%iresult = 2
      RETURN
    END IF

    ! Minimum number of iterations
 
    nminIterations = MAX(rsolverNode%nminIterations,0)
    
    ! Damping parameter
    domega = rsolverNode%domega
      
    ! Use preconditioning? 
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
        CALL output_line ('Space-Time-DefCorr: Iteration '// &
             TRIM(sys_siL(0,10))//',  !!RES!! = '//&
             TRIM(sys_sdEL(rsolverNode%dinitialDefect,15)) )
      END IF

      ! Perform at most nmaxIterations loops to get a new vector

      DO ite = 1,rsolverNode%nmaxIterations
      
        rsolverNode%icurrentIteration = ite
        
        IF (bprec) THEN
          ! Perform preconditioning with the assigned preconditioning
          ! solver structure.
          CALL sptils_precondDefect (p_rprecSubnode,p_rdef)
        END IF
        
        ! In p_rdef, we now have the current residuum $P^{-1} (b-Ax)$.
        ! Add it (damped) to the current iterate p_x to get
        !   $$ x  :=  x  +  \omega P^{-1} (b-Ax) $$

        CALL sptivec_vectorLinearComb (p_rdef ,p_rx,domega,1.0_DP)

        ! Calculate the residuum for the next step : (b-Ax)
        CALL sptivec_copyVector (rd,p_rdef)
        CALL c2d2_spaceTimeMatVec (rsolverNode%p_rproblem, p_rmatrix, &
            p_rx,p_rdef, -1.0_DP,1.0_DP,SPTID_FILTER_DEFECT,dresnorm,.FALSE.)
        
        ! Filter the defect for boundary conditions in space and time.
        !CALL tbc_implementInitCondDefect (&
        !    p_rmatrix%p_rspaceTimeDiscretisation,p_rdef,p_rsubnode%rtempVectorSpace)
        !CALL tbc_implementBCdefect (rsolverNode%p_rproblem,&
        !   p_rmatrix%p_rspaceTimeDiscretisation,p_rdef,p_rsubnode%rtempVectorSpace)
        
        ! Get the norm of the new (final?) residuum
        dfr = sptivec_vectorNorm (p_rdef,rsolverNode%iresNorm)
     
        rsolverNode%dfinalDefect = dfr

        ! Test if the iteration is diverged
        IF (sptils_testDivergence(rsolverNode,dfr)) THEN
          CALL output_line ('Space-Time-DefCorr: Solution diverging!')
          rsolverNode%iresult = 1
          EXIT
        END IF
     
        ! At least perform nminIterations iterations
        IF (ite .GE. nminIterations) THEN
        
          ! Check if the iteration converged
          IF (sptils_testConvergence(rsolverNode,dfr)) EXIT
          
        END IF

        ! print out the current residuum

        IF ((rsolverNode%ioutputLevel .GE. 2) .AND. &
            (MOD(ite,niteResOutput).EQ.0)) THEN
          CALL output_line ('Space-Time-DefCorr: Iteration '// &
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
        CALL output_line ('Space-Time-DefCorr: Iteration '// &
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
        CALL output_line ('Space-Time-DefCorr statistics:')
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
              'Space-Time-DefCorr: Iterations/Rate of convergence: '//&
              TRIM(sys_siL(rsolverNode%iiterations,10))//' /'//&
              TRIM(sys_sdEL(rsolverNode%dconvergenceRate,15)) )
      END IF
      
    ELSE
      ! DEF=Infinity; RHO=Infinity, set to 1
      rsolverNode%dconvergenceRate = 1.0_DP
    END IF  
  
  END SUBROUTINE
  
! *****************************************************************************
! Block Jacobi preconditioner
! *****************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE sptils_initBlockJacobi (rproblem,p_rsolverNode,domega,rspatialPrecond)
  
!<description>
  ! Creates a t_sptilsNode solver structure for the block Jacobi preconditioner.
!</description>
  
!<input>
  ! The problem structure that defines the problem.
  ! A pointer to this is saved in the solver node, so the problem structure
  ! must not be released before the solver node is released.
  TYPE(t_problem), TARGET :: rproblem
  
  ! Damping parameter
  REAL(DP), INTENT(IN) :: domega
  
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
    
    ! Initialise the ability bitfield with the ability of this solver:
    p_rsolverNode%ccapability = SPTILS_ABIL_DIRECT

    p_rsolverNode%domega = domega
    
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
  ! The t_sptilsNode structure which is to be cleaned up.
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
  TYPE(t_ccoptSpaceTimeMatrix), DIMENSION(:), INTENT(IN) :: Rmatrices
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
  ! The t_sptilsNode structure 
  TYPE(t_sptilsNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the SPTILS_ERR_XXXX constants. A value different to 
  ! SPTILS_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  INTEGER, INTENT(OUT) :: ierror
!</output>
  
!</subroutine>
    
    ! A-priori we have no error...
    ierror = SPTILS_ERR_NOERROR
    
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
  ! The t_sptilsNode structure 
  TYPE(t_sptilsNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the SPTILS_ERR_XXXX constants. A value different to 
  ! SPTILS_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  INTEGER, INTENT(OUT) :: ierror
!</output>

!</subroutine>

    ! A-priori we have no error...
    ierror = SPTILS_ERR_NOERROR

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
  ! The t_sptilsNode structure 
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
  ! The t_sptilsNode structure 
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
    INTEGER :: isubstep,iidxNonlin
    REAL(DP) :: dtheta,dtstep
    LOGICAL :: bsuccess
    TYPE(t_ccmatrixComponents) :: rmatrixComponents
    TYPE(t_ccoptSpaceTimeDiscretisation), POINTER :: p_rspaceTimeDiscr
    TYPE(t_ccoptSpaceTimeMatrix), POINTER :: p_rspaceTimeMatrix
    TYPE(t_vectorBlock) :: rtempVectorD,rtempVectorX
    
    ! DEBUG!!!
    REAL(DP), DIMENSION(:), POINTER :: p_Dx,p_Dd,p_Dsol
    
    ! Get a pointer to the space-time discretisation structure that defines
    ! how to apply the global system matrix.
    p_rspaceTimeDiscr => rsolverNode%rmatrix%p_rspaceTimeDiscretisation
    p_rspaceTimeMatrix => rsolverNode%rmatrix
    
    dtheta = rsolverNode%p_rproblem%rtimedependence%dtimeStepTheta
    dtstep = p_rspaceTimeDiscr%rtimeDiscr%dtstep

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
        
    rmatrixComponents%dnu = &
        collct_getvalue_real (rsolverNode%p_rproblem%rcollection,'NU')
    rmatrixComponents%iupwind1 = &
        collct_getvalue_int (rsolverNode%p_rproblem%rcollection,'IUPWIND1')
    rmatrixComponents%dupsam1 = &
        collct_getvalue_real (rsolverNode%p_rproblem%rcollection,'UPSAM1')
    rmatrixComponents%iupwind2 = &
        collct_getvalue_int (rsolverNode%p_rproblem%rcollection,'IUPWIND2')
    rmatrixComponents%dupsam2 = &
        collct_getvalue_real (rsolverNode%p_rproblem%rcollection,'UPSAM2')
    
    ! ----------------------------------------------------------------------
    ! We use a block-Jacobi scheme for preconditioning...
    !
    ! For this purpose, loop through the substeps.
    
    DO isubstep = 0,p_rspaceTimeDiscr%NEQtime-1
    
      ! Current time step?
      rsolverNode%p_rproblem%rtimedependence%dtime = &
          rsolverNode%p_rproblem%rtimedependence%dtimeInit + isubstep * dtstep

      IF (rsolverNode%ioutputLevel .GE. 1) THEN
        CALL output_line ('Space-Time-Block-Jacobi preconditioning of timestep: '//&
            TRIM(sys_siL(isubstep,10))//&
            ', Time: '//TRIM(sys_sdL(rsolverNode%p_rproblem%rtimedependence%dtime,10)))
      END IF
    
      ! -----
      ! Discretise the boundary conditions at the new point in time -- 
      ! if the boundary conditions are nonconstant in time!
      IF (collct_getvalue_int (rsolverNode%p_rproblem%rcollection,'IBOUNDARY') &
          .NE. 0) THEN
        CALL c2d2_updateDiscreteBC (rsolverNode%p_rproblem, .FALSE.)
      END IF

      ! DEBUG!!!      
      CALL lsysbl_getbase_double (rtempVectorX,p_Dx)
      CALL lsysbl_getbase_double (rtempVectorD,p_Dd)

      ! Read in the RHS/solution/defect vector of the current timestep.
      ! If no solution is specified, we have a linear problem and thus
      ! the content of rtempVector is not relevant; actually it's even
      ! zero by initialisation.
      IF (ASSOCIATED(p_rspaceTimeMatrix%p_rsolution)) THEN
        CALL sptivec_getTimestepData (p_rspaceTimeMatrix%p_rsolution, &
            1+isubstep, rtempVectorX)
        
        ! DEBUG!!!
        CALL lsysbl_getbase_double (rtempVectorX,p_Dsol)
      END IF
      CALL sptivec_getTimestepData (rd, 1+isubstep, rtempVectorD)

      ! Set up the matrix weights for the diagonal matrix
      CALL c2d2_setupMatrixWeights (rsolverNode%p_rproblem,p_rspaceTimeMatrix,dtheta,&
        isubstep,0,rmatrixComponents,iidxNonlin)
        
      ! Perform preconditioning of the spatial defect with the method provided by the
      ! core equation module.
      CALL c2d2_precondDefect (&
          rsolverNode%p_rsubnodeBlockJacobi%p_rspatialPreconditioner,&
          rmatrixComponents,&
          rtempVectorD,rtempVectorX,bsuccess,rsolverNode%p_rproblem%rcollection)      
    
      ! Scale by omega
      CALL lsysbl_scaleVector (rtempVectorD,rsolverNode%domega)
    
      ! Save back the preconditioned defect.
      CALL sptivec_setTimestepData (rd, 1+isubstep, rtempVectorD)
      
    END DO
    
    CALL lsysbl_releaseVector (rtempVectorX)
    CALL lsysbl_releaseVector (rtempVectorD)
    
  END SUBROUTINE

! *****************************************************************************
! Block Gauss-Seidel preconditioner
! *****************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE sptils_initBlockSOR (rproblem,p_rsolverNode,&
      domega,rspatialPrecond)
  
!<description>
  ! Creates a t_sptilsNode solver structure for the block 
  ! Gauss Seidel preconditioner.
!</description>
  
!<input>
  ! The problem structure that defines the problem.
  ! A pointer to this is saved in the solver node, so the problem structure
  ! must not be released before the solver node is released.
  TYPE(t_problem), TARGET :: rproblem
  
  ! Relaxation parameter. =1.0: Block Gauss-Seidel.
  REAL(DP), INTENT(IN) :: domega
  
  ! A spatial preconditioner structure that defines how to perform the preconditioning
  ! in each substep.
  ! A pointer to this is noted in p_rsolverNode, so rspatialPrecond
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
    p_rsolverNode%calgorithm = SPTILS_ALG_BLOCKSOR
    
    ! Initialise the ability bitfield with the ability of this solver:
    p_rsolverNode%ccapability = SPTILS_ABIL_DIRECT

    ! Save the relaxation parameter
    p_rsolverNode%domega = domega
    
    ! Allocate a subnode for our solver.
    ! This initialises most of the variables with default values appropriate
    ! to this solver.
    ALLOCATE(p_rsolverNode%p_rsubnodeBlockSOR)
    
    ! Attach the preconditioner if given. 
    p_rsolverNode%p_rsubnodeBlockSOR%p_rspatialPreconditioner => &
        rspatialPrecond
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE sptils_doneBlockSOR (rsolverNode)
  
!<description>
  ! This routine releases all temporary memory for the solver 
  ! from the heap.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure which is to be cleaned up.
  TYPE(t_sptilsNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!</subroutine>
  
    ! Release memory if still associated
    CALL sptils_doneDataBlockSOR (rsolverNode)
    CALL sptils_doneStructureBlockSOR (rsolverNode)
    
    ! Release the subnode structure
    DEALLOCATE(rsolverNode%p_rsubnodeBlockSOR)
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE sptils_setMatrixBlockSOR (rsolverNode,Rmatrices)
  
!<description>
  ! This routine is called if the system matrix changes.
  ! The routine calls sptils_setMatrices for the preconditioner
  ! to inform also that one about the change of the matrix pointer.
!</description>
  
!<input>
  ! An array of system matrices which is simply passed to the initialisation 
  ! routine of the preconditioner.
  TYPE(t_ccoptSpaceTimeMatrix), DIMENSION(:), INTENT(IN) :: Rmatrices
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
  
  RECURSIVE SUBROUTINE sptils_initStructureBlockSOR (rsolverNode,ierror)
  
!<description>
  ! Solver preparation. Perform symbolic factorisation (not of the defect
  ! correcion solver, but of subsolvers). Allocate temporary memory.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure 
  TYPE(t_sptilsNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the SPTILS_ERR_XXXX constants. A value different to 
  ! SPTILS_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  INTEGER, INTENT(OUT) :: ierror
!</output>
  
!</subroutine>
    
    ! A-priori we have no error...
    ierror = SPTILS_ERR_NOERROR
    
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE sptils_initDataBlockSOR (rsolverNode, ierror)
  
!<description>
  ! Calls the initData subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with sptils_initData.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure 
  TYPE(t_sptilsNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the SPTILS_ERR_XXXX constants. A value different to 
  ! SPTILS_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  INTEGER, INTENT(OUT) :: ierror
!</output>

!</subroutine>

    ! A-priori we have no error...
    ierror = SPTILS_ERR_NOERROR

  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE sptils_doneDataBlockSOR (rsolverNode)
  
!<description>
  ! Calls the doneData subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with sptils_doneData.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure 
  TYPE(t_sptilsNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!</subroutine>

  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE sptils_doneStructureBlockSOR (rsolverNode)
  
!<description>
  ! Calls the doneStructure subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with sptils_doneStructure.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure 
  TYPE(t_sptilsNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!</subroutine>

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>
  
  RECURSIVE SUBROUTINE sptils_precBlockSOR (rsolverNode,rd)
  
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
    INTEGER :: isubstep,iidxNonlin
    REAL(DP) :: dtheta,dtstep
    LOGICAL :: bsuccess
    TYPE(t_ccmatrixComponents) :: rmatrixComponents
    TYPE(t_ccoptSpaceTimeDiscretisation), POINTER :: p_rspaceTimeDiscr
    TYPE(t_ccoptSpaceTimeMatrix), POINTER :: p_rspaceTimeMatrix
    TYPE(t_vectorBlock) :: rtempVectorD1,rtempVectorD2,rtempVectorD3,rtempVectorX
    
    ! DEBUG!!!
    REAL(DP), DIMENSION(:), POINTER :: p_Dx,p_Dd
        
    ! Get a pointer to the space-time discretisation structure that defines
    ! how to apply the global system matrix.
    p_rspaceTimeDiscr => rsolverNode%rmatrix%p_rspaceTimeDiscretisation
    p_rspaceTimeMatrix => rsolverNode%rmatrix
    
    dtheta = rsolverNode%p_rproblem%rtimedependence%dtimeStepTheta
    dtstep = p_rspaceTimeDiscr%rtimeDiscr%dtstep

    ! Create temp vectors
    CALL lsysbl_createVecBlockByDiscr (p_rspaceTimeDiscr%p_rlevelInfo%p_rdiscretisation,&
        rtempVectorD1,.TRUE.)
    CALL lsysbl_createVecBlockByDiscr (p_rspaceTimeDiscr%p_rlevelInfo%p_rdiscretisation,&
        rtempVectorD2,.TRUE.)
    CALL lsysbl_createVecBlockByDiscr (p_rspaceTimeDiscr%p_rlevelInfo%p_rdiscretisation,&
        rtempVectorD3,.TRUE.)

    CALL lsysbl_createVecBlockByDiscr (p_rspaceTimeDiscr%p_rlevelInfo%p_rdiscretisation,&
        rtempVectorX,.TRUE.)
        
    ! Attach the boundary conditions to the temp vectors.
    rtempVectorD2%p_rdiscreteBC => p_rspaceTimeDiscr%p_rlevelInfo%p_rdiscreteBC
    rtempVectorD2%p_rdiscreteBCfict => p_rspaceTimeDiscr%p_rlevelInfo%p_rdiscreteFBC

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
        
    rmatrixComponents%dnu = &
        collct_getvalue_real (rsolverNode%p_rproblem%rcollection,'NU')
    rmatrixComponents%iupwind1 = &
        collct_getvalue_int (rsolverNode%p_rproblem%rcollection,'IUPWIND1')
    rmatrixComponents%dupsam1 = &
        collct_getvalue_real (rsolverNode%p_rproblem%rcollection,'UPSAM1')
    rmatrixComponents%iupwind2 = &
        collct_getvalue_int (rsolverNode%p_rproblem%rcollection,'IUPWIND2')
    rmatrixComponents%dupsam2 = &
        collct_getvalue_real (rsolverNode%p_rproblem%rcollection,'UPSAM2')

    ! ----------------------------------------------------------------------
    ! We use a block-Gauss-Seidel scheme for preconditioning...
    !
    ! The system has the form
    !
    ! AA                x1 = dx1
    ! AA   M            l1   dl1
    ! M   AA            x2   dx2
    !     AA   M        l2   dl2
    !     M   AA        x3   dx3
    !         AA   M    l3   dl3
    !         M   AA    x4   dx4
    !             AA    l4   dl4
    !
    ! We have to note that the primal system advances forward in time
    ! while the dual system marches backward in time!
    ! For a Gauss-Seidel approach, we ignore the dual velocity, which
    ! gives us a lower-triangular system
    !
    ! AA                x1 = dx1
    ! AA                l1   dl1
    ! M   AA            x2   dx2
    !     AA            l2   dl2
    !     M   AA        x3   dx3
    !         AA        l3   dl3
    !         M   AA    x4   dx4
    !             AA    l4   dl4
    !
    ! So we throw the primal velocity of the previous timestep to the RHS 
    ! and then use Block Jacobi preconditioning to get an update.
    ! Here, the primal velocity of the previous timestep is weighted
    ! by omega.
    ! Example for x2/l2:
    !
    !                   x1 =    
    !                   l1      
    ! M   AA            x2   dx2
    !     AA            l2   dl2
    !                   x3      
    !                   l3      
    !
    ! =>
    !
    !     AA            x2 = dx2 - domega Mx1 
    !     AA            l2   dl2             
    !
    ! =>
    !
    !  x2^new = AA^-1 ( dx2 - domega Mx1 )
    !  l2^new   AA    ( dl2              )
    !
    ! domega=1 gives the usual Gauss-Seidel approach.
    !
    ! Ok, that's not the reall SOR algorithm (a constant multiplier is still
    ! missing), but a similar thing...
    
    ! Load the solution of the 0th timestep.
    CALL sptivec_getTimestepData (rd, 1+0, rtempVectorD2)
    
    ! Loop through the substeps we have to update
    DO isubstep = 0,p_rspaceTimeDiscr%NEQtime-1
    
      ! Current time step?
      rsolverNode%p_rproblem%rtimedependence%dtime = &
          rsolverNode%p_rproblem%rtimedependence%dtimeInit + isubstep * dtstep

      IF (rsolverNode%ioutputLevel .GE. 1) THEN
        CALL output_line ('Space-Time-Block-SOR preconditioning of timestep: '//&
            TRIM(sys_siL(isubstep,10))//&
            ', Time: '//TRIM(sys_sdL(rsolverNode%p_rproblem%rtimedependence%dtime,10)))
      END IF
    
      ! -----
      ! Discretise the boundary conditions at the new point in time -- 
      ! if the boundary conditions are nonconstant in time!
      IF (collct_getvalue_int (rsolverNode%p_rproblem%rcollection,'IBOUNDARY') &
          .NE. 0) THEN
        CALL c2d2_updateDiscreteBC (rsolverNode%p_rproblem, .FALSE.)
      END IF

      ! DEBUG!!!      
      CALL lsysbl_getbase_double (rtempVectorX,p_Dx)
      CALL lsysbl_getbase_double (rtempVectorD2,p_Dd)

      ! Read in the RHS/solution/defect vector of the current timestep.
      ! If no solution is specified, we have a linear problem and thus
      ! the content of rtempVector is not relevant; actually it's even
      ! zero by initialisation.
      IF (ASSOCIATED(p_rspaceTimeMatrix%p_rsolution)) THEN
        CALL sptivec_getTimestepData (p_rspaceTimeMatrix%p_rsolution, &
            1+isubstep, rtempVectorX)
      END IF
      

      ! Is this the first timestep or not?
      IF (isubstep .GT. 0) THEN
      
        ! Create d2 = RHS - Mx1 
        CALL c2d2_setupMatrixWeights (rsolverNode%p_rproblem,p_rspaceTimeMatrix,dtheta,&
          isubstep,-1,rmatrixComponents,iidxNonlin)
          
        CALL c2d2_assembleDefect (rmatrixComponents,rtempVectorD1,rtempVectorD2,rsolverNode%domega)
        
      END IF

      ! Is this the last timestep or not?
      IF (isubstep .LT. p_rspaceTimeDiscr%NEQtime-1) THEN

        ! Read the RHS of the next timestep
        CALL sptivec_getTimestepData (rd, 1+isubstep+1, rtempVectorD3)
      
      END IF

      ! Set up the matrix weights for the diagonal matrix
      CALL c2d2_setupMatrixWeights (rsolverNode%p_rproblem,p_rspaceTimeMatrix,dtheta,&
        isubstep,0,rmatrixComponents,iidxNonlin)

      ! Perform preconditioning of the spatial defect with the method 
      ! provided by the core equation module.
      CALL c2d2_precondDefect (&
          rsolverNode%p_rsubnodeBlockSOR%p_rspatialPreconditioner,&
          rmatrixComponents,&
          rtempVectorD2,rtempVectorX,bsuccess,rsolverNode%p_rproblem%rcollection)      
    
      ! Save back the preconditioned defect.
      CALL sptivec_setTimestepData (rd, 1+isubstep, rtempVectorD2)
      
      ! Shift the RHS vectors: 1 <- 2 <- 3
      CALL lsysbl_copyVector (rtempVectorD2,rtempVectorD1)
      CALL lsysbl_copyVector (rtempVectorD3,rtempVectorD2)
      
    END DO

    CALL lsysbl_releaseVector (rtempVectorX)
    CALL lsysbl_releaseVector (rtempVectorD3)
    CALL lsysbl_releaseVector (rtempVectorD2)
    CALL lsysbl_releaseVector (rtempVectorD1)
    
  END SUBROUTINE

! *****************************************************************************
! Block FBGS preconditioner
! *****************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE sptils_initBlockFBGS (rproblem,p_rsolverNode,&
      domega,domegaSOR,rspatialPrecond)
  
!<description>
  ! Creates a t_sptilsNode solver structure for the block 
  ! FBGS preconditioner.
!</description>
  
!<input>
  ! The problem structure that defines the problem.
  ! A pointer to this is saved in the solver node, so the problem structure
  ! must not be released before the solver node is released.
  TYPE(t_problem), TARGET :: rproblem
  
  ! Relaxation parameter. =1.0: Full block FBGS.
  REAL(DP), INTENT(IN) :: domega

  ! Relaxation parameter for the time dependence. =1.0=Standard.
  REAL(DP), INTENT(IN) :: domegaSOR
  
  ! A spatial preconditioner structure that defines how to perform the preconditioning
  ! in each substep.
  ! A pointer to this is noted in p_rsolverNode, so rspatialPrecond
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
    p_rsolverNode%calgorithm = SPTILS_ALG_BLOCKFBGS
    
    ! Initialise the ability bitfield with the ability of this solver:
    p_rsolverNode%ccapability = 0

    ! Save the relaxation parameter
    p_rsolverNode%domega = domega
    
    ! Allocate a subnode for our solver.
    ! This initialises most of the variables with default values appropriate
    ! to this solver.
    ALLOCATE(p_rsolverNode%p_rsubnodeBlockFBGS)
    
    ! Save the relaxation parameter for the SOR-aproach
    p_rsolverNode%p_rsubnodeBlockFBGS%domegaSOR = domegaSOR

    ! Attach the preconditioner if given. 
    p_rsolverNode%p_rsubnodeBlockFBGS%p_rspatialPreconditioner => &
        rspatialPrecond
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE sptils_doneBlockFBGS (rsolverNode)
  
!<description>
  ! This routine releases all temporary memory for the solver 
  ! from the heap.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure which is to be cleaned up.
  TYPE(t_sptilsNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!</subroutine>
  
    ! Release memory if still associated
    CALL sptils_doneDataBlockSOR (rsolverNode)
    CALL sptils_doneStructureBlockSOR (rsolverNode)
    
    ! Release the subnode structure
    DEALLOCATE(rsolverNode%p_rsubnodeBlockFBGS)
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE sptils_setMatrixBlockFBGS (rsolverNode,Rmatrices)
  
!<description>
  ! This routine is called if the system matrix changes.
  ! The routine calls sptils_setMatrices for the preconditioner
  ! to inform also that one about the change of the matrix pointer.
!</description>
  
!<input>
  ! An array of system matrices which is simply passed to the initialisation 
  ! routine of the preconditioner.
  TYPE(t_ccoptSpaceTimeMatrix), DIMENSION(:), INTENT(IN) :: Rmatrices
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
  
  RECURSIVE SUBROUTINE sptils_initStructureBlockFBGS (rsolverNode,ierror)
  
!<description>
  ! Solver preparation. Perform symbolic factorisation (not of the defect
  ! correcion solver, but of subsolvers). Allocate temporary memory.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure 
  TYPE(t_sptilsNode), INTENT(INOUT), TARGET :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the SPTILS_ERR_XXXX constants. A value different to 
  ! SPTILS_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  INTEGER, INTENT(OUT) :: ierror
!</output>
  
!</subroutine>

    TYPE(t_ccoptspaceTimeDiscretisation), POINTER :: p_rspaceTimeDiscr
    
    p_rspaceTimeDiscr => rsolverNode%rmatrix%p_rspaceTimeDiscretisation

    ! Allocate memory for a temp vector.
    CALL sptivec_initVector (rsolverNode%p_rsubnodeBlockFBGS%rtempVector,&
        p_rspaceTimeDiscr%NEQtime,&
        p_rspaceTimeDiscr%p_rlevelInfo%p_rdiscretisation)

    ! A-priori we have no error...
    ierror = SPTILS_ERR_NOERROR
    
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE sptils_initDataBlockFBGS (rsolverNode, ierror)
  
!<description>
  ! Calls the initData subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with sptils_initData.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure 
  TYPE(t_sptilsNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the SPTILS_ERR_XXXX constants. A value different to 
  ! SPTILS_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  INTEGER, INTENT(OUT) :: ierror
!</output>

!</subroutine>

    ! A-priori we have no error...
    ierror = SPTILS_ERR_NOERROR

  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE sptils_doneDataBlockFBGS (rsolverNode)
  
!<description>
  ! Calls the doneData subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with sptils_doneData.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure 
  TYPE(t_sptilsNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!</subroutine>

  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE sptils_doneStructureBlockFBGS (rsolverNode)
  
!<description>
  ! Calls the doneStructure subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with sptils_doneStructure.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure 
  TYPE(t_sptilsNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!</subroutine>

    IF (rsolverNode%p_rsubnodeBlockFBGS%rtempVector%NEQtime .NE. 0) &
      CALL sptivec_releaseVector(rsolverNode%p_rsubnodeBlockFBGS%rtempVector)

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>
  
  RECURSIVE SUBROUTINE sptils_precBlockFBGS (rsolverNode,rd)
  
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
    INTEGER :: isubstep,iiteration,iidxNonlin
    REAL(DP) :: dtheta,dtstep,dtime
    LOGICAL :: bsuccess
    TYPE(t_ccmatrixComponents) :: rmatrixComponents
    TYPE(t_ccoptSpaceTimeDiscretisation), POINTER :: p_rspaceTimeDiscr
    TYPE(t_ccoptSpaceTimeMatrix), POINTER :: p_rspaceTimeMatrix
    TYPE(t_vectorBlock) :: rtempVectorD1,rtempVectorD2,rtempVectorD3
    TYPE(t_vectorBlock) :: rtempVectorX1,rtempVectorX2,rtempVectorX3
    TYPE(t_vectorBlock) :: rtempVectorSol
    TYPE(t_vectorBlock) :: rtempVectorRHS
    TYPE(t_spacetimeVector), POINTER :: p_rx
    LOGICAL :: bcalcNorm
    REAL(DP) :: domegaSOR,dres,dresInit
    
    ! DEBUG!!!
    REAL(DP), DIMENSION(:), POINTER :: p_Dx1,p_Dx2,p_Dx3,p_Dd1,p_Dd2,p_Dd3,p_Dd,p_Dsol
        
    ! Get a pointer to the space-time discretisation structure that defines
    ! how to apply the global system matrix.
    p_rspaceTimeDiscr => rsolverNode%rmatrix%p_rspaceTimeDiscretisation
    p_rspaceTimeMatrix => rsolverNode%rmatrix
    
    dtheta = rsolverNode%p_rproblem%rtimedependence%dtimeStepTheta
    dtstep = p_rspaceTimeDiscr%rtimeDiscr%dtstep
    
    domegaSOR = rsolverNode%p_rsubnodeBlockFBGS%domegaSOR

    ! Create temp vectors
    CALL lsysbl_createVecBlockByDiscr (p_rspaceTimeDiscr%p_rlevelInfo%p_rdiscretisation,&
        rtempVectorD1,.TRUE.)
    CALL lsysbl_createVecBlockByDiscr (p_rspaceTimeDiscr%p_rlevelInfo%p_rdiscretisation,&
        rtempVectorD2,.TRUE.)
    CALL lsysbl_createVecBlockByDiscr (p_rspaceTimeDiscr%p_rlevelInfo%p_rdiscretisation,&
        rtempVectorD3,.TRUE.)
    CALL lsysbl_createVecBlockByDiscr (p_rspaceTimeDiscr%p_rlevelInfo%p_rdiscretisation,&
        rtempVectorRHS,.TRUE.)

    CALL lsysbl_createVecBlockByDiscr (p_rspaceTimeDiscr%p_rlevelInfo%p_rdiscretisation,&
        rtempVectorX1,.TRUE.)
    CALL lsysbl_createVecBlockByDiscr (p_rspaceTimeDiscr%p_rlevelInfo%p_rdiscretisation,&
        rtempVectorX2,.TRUE.)
    CALL lsysbl_createVecBlockByDiscr (p_rspaceTimeDiscr%p_rlevelInfo%p_rdiscretisation,&
        rtempVectorX3,.TRUE.)

    ! Solution vector -- for setting up the defect in nonlinear problems.
    CALL lsysbl_createVecBlockByDiscr (p_rspaceTimeDiscr%p_rlevelInfo%p_rdiscretisation,&
        rtempVectorSol,.TRUE.)
        
    ! DEBUG!!!      
    CALL lsysbl_getbase_double (rtempVectorX1,p_Dx1)
    CALL lsysbl_getbase_double (rtempVectorX2,p_Dx2)
    CALL lsysbl_getbase_double (rtempVectorX3,p_Dx3)
    CALL lsysbl_getbase_double (rtempVectorD1,p_Dd1)
    CALL lsysbl_getbase_double (rtempVectorD2,p_Dd2)
    CALL lsysbl_getbase_double (rtempVectorD3,p_Dd3)
    CALL lsysbl_getbase_double (rtempVectorRHS,p_Dd)
        
    ! Attach the boundary conditions to the temp vectors.
    rtempVectorD1%p_rdiscreteBC => p_rspaceTimeDiscr%p_rlevelInfo%p_rdiscreteBC
    rtempVectorD1%p_rdiscreteBCfict => p_rspaceTimeDiscr%p_rlevelInfo%p_rdiscreteFBC

    rtempVectorD2%p_rdiscreteBC => p_rspaceTimeDiscr%p_rlevelInfo%p_rdiscreteBC
    rtempVectorD2%p_rdiscreteBCfict => p_rspaceTimeDiscr%p_rlevelInfo%p_rdiscreteFBC

    rtempVectorD3%p_rdiscreteBC => p_rspaceTimeDiscr%p_rlevelInfo%p_rdiscreteBC
    rtempVectorD3%p_rdiscreteBCfict => p_rspaceTimeDiscr%p_rlevelInfo%p_rdiscreteFBC

    rtempVectorSol%p_rdiscreteBC => p_rspaceTimeDiscr%p_rlevelInfo%p_rdiscreteBC
    rtempVectorSol%p_rdiscreteBCfict => p_rspaceTimeDiscr%p_rlevelInfo%p_rdiscreteFBC

    rtempVectorRHS%p_rdiscreteBC => p_rspaceTimeDiscr%p_rlevelInfo%p_rdiscreteBC
    rtempVectorRHS%p_rdiscreteBCfict => p_rspaceTimeDiscr%p_rlevelInfo%p_rdiscreteFBC

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
        
    rmatrixComponents%dnu = &
        collct_getvalue_real (rsolverNode%p_rproblem%rcollection,'NU')
    rmatrixComponents%iupwind1 = &
        collct_getvalue_int (rsolverNode%p_rproblem%rcollection,'IUPWIND1')
    rmatrixComponents%dupsam1 = &
        collct_getvalue_real (rsolverNode%p_rproblem%rcollection,'UPSAM1')
    rmatrixComponents%iupwind2 = &
        collct_getvalue_int (rsolverNode%p_rproblem%rcollection,'IUPWIND2')
    rmatrixComponents%dupsam2 = &
        collct_getvalue_real (rsolverNode%p_rproblem%rcollection,'UPSAM2')

    ! Probably we have to calclate the norm of the residual while calculating...
    bcalcNorm = (rsolverNode%nminIterations .NE. rsolverNode%nmaxIterations) .OR.&
                (rsolverNode%ioutputLevel .GE. 2)

    ! ----------------------------------------------------------------------
    ! We use a Block-FBGS scheme for preconditioning.
    !
    ! This is an iterative approach that performs nminIterations global
    ! defect correction steps.
    !
    ! The system has the form
    !
    ! AA                x1 = dx1
    ! AA   M            l1   dl1
    ! M   AA            x2   dx2
    !     AA   M        l2   dl2
    !     M   AA        x3   dx3
    !         AA   M    l3   dl3
    !         M   AA    x4   dx4
    !             AA    l4   dl4
    !
    ! We have to note that the primal system advances forward in time
    ! while the dual system marches backward in time!
    ! In a first step, we reduce the system to the current time step.
    ! Example for x2/l2:
    !
    !                   x1 =    
    !                   l1      
    ! M   AA            x2   dx2
    !     AA   M        l2   dl2
    !                   x3      
    !                   l3      
    !
    ! We throw everything to the RHS to get a 'local' defect -- the current
    ! velocity (here x2/l2) as well as the primal velocity terms of the previous
    ! time step as well as the dual velocity of the next timestep.
    ! =>
    !
    !                   bx2 := dx2 - Mx1       - AA x2
    !                   bl2    dl2       - Ml3   AA l2
    !
    ! We perform preconditioning with the current system matrix and add the
    ! preconditioned defect to the current solution vector:
    ! =>
    !
    !  x2^new = x2 + AA^-1 ( dx2 - Mx1       - AA x2 ) = x2 + AA^-1 ( bx2 )
    !  l2^new   l2   AA    ( dl2       - Ml3   AA l2 )   l2   AA    ( bl2 )
    !
    ! The whole thing is an iterative process starting with x=l=0.
    ! We repeat this process nminIteration times to get a new preconditioned
    ! defect ( x,l ) which then replaces rd.
    
    p_rx   => rsolverNode%p_rsubnodeBlockFBGS%rtempVector
    
    ! Ok, let's start. Initial solution is zero.
    CALL sptivec_clearVector (p_rx)
    
    ! Norm of the initial residuum
    IF (bcalcNorm) THEN
      dresInit = sptivec_vectorNorm (rd,LINALG_NORML2)
      IF (dresInit .EQ. 0.0_DP) dresInit = 1.0_DP
      rsolverNode%dinitialDefect = dresInit
      dres = dresInit
    END IF
    
    DO iiteration = 1,rsolverNode%nmaxIterations

      ! Probably print the current residuum (of the previous step)
      IF (rsolverNode%nminIterations .NE. rsolverNode%nmaxIterations) THEN
        IF (rsolverNode%ioutputLevel .GE. 2) THEN
          dres = SQRT(dres / REAL(p_rspaceTimeDiscr%NEQtime,DP))
          CALL output_line ('Space-Time-Block-FBGS: Iteration '// &
              TRIM(sys_siL(iiteration-1,10))//',  !!RES!! = '//&
              TRIM(sys_sdEL(dres,15)) )
        END IF
        
        ! Check for convergence
        IF ((iiteration .GT. rsolverNode%nminIterations) .AND. &
            (iiteration .LT. rsolverNode%nmaxIterations)) THEN
          IF (sptils_testConvergence (rsolverNode, dres)) EXIT
          IF (sptils_testDivergence (rsolverNode, dres)) EXIT
        END IF
      END IF

      ! Filter the current solution for boundary conditions in space and time.
      ! rtempVectorRHS is used as temp vector here.
      CALL tbc_implementInitCondDefect (&
          p_rspaceTimeDiscr,p_rx,rtempVectorRHS)
      CALL tbc_implementBCdefect (rsolverNode%p_rproblem,&
         p_rspaceTimeDiscr,p_rx,rtempVectorRHS)

      ! -----
      ! Backward in time
      ! -----

      ! Load the RHS and solution of the 0th timestep.
      CALL sptivec_getTimestepData (rd, &
          p_rspaceTimeDiscr%NEQtime, rtempVectorD2)
      
      ! Current iterate
      CALL sptivec_getTimestepData (p_rx, &
          p_rspaceTimeDiscr%NEQtime, rtempVectorX2)

      ! Loop through the substeps we have to update
      DO isubstep = p_rspaceTimeDiscr%NEQtime-1,0,-1
      
        ! Current point in time
        dtime = p_rspaceTimeDiscr%rtimeDiscr%dtimeInit + isubstep * dtstep

        rsolverNode%p_rproblem%rtimedependence%dtime = dtime

        IF (rsolverNode%ioutputLevel .GE. 1) THEN
          CALL output_line ('Space-Time-Block-FBGS preconditioning of timestep: '//&
              TRIM(sys_siL(isubstep,10))//&
              ', Time: '//TRIM(sys_sdL(rsolverNode%p_rproblem%rtimedependence%dtime,10)))
        END IF
      
        ! -----
        ! Discretise the boundary conditions at the new point in time -- 
        ! if the boundary conditions are nonconstant in time!
        IF (collct_getvalue_int (rsolverNode%p_rproblem%rcollection,'IBOUNDARY') &
            .NE. 0) THEN
          CALL c2d2_updateDiscreteBC (rsolverNode%p_rproblem, .FALSE.)
        END IF

        ! The RHS which is put into the preconditioner is set up in 
        ! rtempVectorRHS to prevent rtempVectorD2 from getting destroyed.
        CALL lsysbl_copyVector (rtempVectorD2,rtempVectorRHS)
        
        ! Is this the first timestep or not?
        IF (isubstep .GT. 0) THEN
        
          ! Read the RHS and solution of the next timestep
          CALL sptivec_getTimestepData (rd, 1+isubstep-1, rtempVectorD1)
          CALL sptivec_getTimestepData (p_rx, 1+isubstep-1, rtempVectorX1)

          ! Create d2 = RHS - Mx1 
          CALL c2d2_setupMatrixWeights (rsolverNode%p_rproblem,p_rspaceTimeMatrix,dtheta,&
            isubstep,-1,rmatrixComponents,iidxNonlin)
            
          CALL c2d2_assembleDefect (rmatrixComponents,rtempVectorX1,&
            rtempVectorRHS,domegaSOR)
          
        END IF

        ! Is this the last timestep or not?
        IF (isubstep .LT. p_rspaceTimeDiscr%NEQtime-1) THEN
          
          ! Create d2 = RHS - Ml3 
          CALL c2d2_setupMatrixWeights (rsolverNode%p_rproblem,p_rspaceTimeMatrix,dtheta,&
            isubstep,1,rmatrixComponents,iidxNonlin)
            
          CALL c2d2_assembleDefect (rmatrixComponents,rtempVectorX3,&
            rtempVectorRHS,domegaSOR)
          
        END IF

        ! Read in the solution vector of the current timestep (for nonlinear problems).
        IF (ASSOCIATED(p_rspaceTimeMatrix%p_rsolution)) THEN
          CALL sptivec_getTimestepData (p_rspaceTimeMatrix%p_rsolution, &
              1+isubstep, rtempVectorSol)

          ! DEBUG!!!
          CALL lsysbl_getbase_double (rtempVectorSol,p_Dsol)
        END IF  

        ! Set up the matrix weights for the diagonal matrix
        CALL c2d2_setupMatrixWeights (rsolverNode%p_rproblem,p_rspaceTimeMatrix,dtheta,&
          isubstep,0,rmatrixComponents,iidxNonlin)
          
        ! Create d2 = RHS - A(solution) X2
        CALL c2d2_assembleDefect (rmatrixComponents,rtempVectorX2,rtempVectorRHS,&
            1.0_DP,rtempVectorSol)
            
        ! Filter the defect for BC's and initial conditions if necessary
        IF (isubstep .EQ. 0) THEN
          CALL tbc_implementInitCondDefSingle (p_rspaceTimeDiscr, rtempVectorRHS)
        ELSE IF (isubstep .EQ. p_rspaceTimeDiscr%NEQtime-1) THEN
          CALL tbc_implementTermCondDefSingle (p_rspaceTimeDiscr, rtempVectorRHS)
        END IF

        CALL tbc_implementSpatialBCdefect (&
            rsolverNode%p_rproblem,isubstep,dtime,p_rspaceTimeDiscr,rtempVectorRHS)
        
        ! Ok, we have the local defect.
        !
        ! Perform preconditioning of the spatial defect with the method 
        ! provided by the core equation module.
        CALL c2d2_precondDefect (&
            rsolverNode%p_rsubnodeBlockFBGS%p_rspatialPreconditioner,&
            rmatrixComponents,&
            rtempVectorRHS,rtempVectorSol,bsuccess,rsolverNode%p_rproblem%rcollection)      
      
        ! Add that defect to the current solution -- damped by domega.
        CALL lsysbl_vectorLinearComb (rtempVectorRHS,rtempVectorX2,&
            rsolverNode%domega,1.0_DP)
      
        ! Save the new solution.
        CALL sptivec_setTimestepData (p_rx, 1+isubstep, rtempVectorX2)
        
        ! Shift the RHS/solution vectors: 1 -> 2 -> 3
        CALL lsysbl_copyVector (rtempVectorD2,rtempVectorD3)
        CALL lsysbl_copyVector (rtempVectorD1,rtempVectorD2)

        CALL lsysbl_copyVector (rtempVectorX2,rtempVectorX3)
        CALL lsysbl_copyVector (rtempVectorX1,rtempVectorX2)

      END DO
      
      ! -----
      ! Forward in time
      ! -----

      ! Norm of the residuum
      dres = 0.0_DP

      ! Load the RHS and solution of the 0th timestep.
      CALL sptivec_getTimestepData (rd, 1+0, rtempVectorD2)
      
      ! Current iterate
      CALL sptivec_getTimestepData (p_rx, 1+0, rtempVectorX2)

      ! Loop through the substeps we have to update
      DO isubstep = 0,p_rspaceTimeDiscr%NEQtime-1
      
        ! Current point in time
        dtime = p_rspaceTimeDiscr%rtimeDiscr%dtimeInit + isubstep * dtstep

        rsolverNode%p_rproblem%rtimedependence%dtime = dtime

        IF (rsolverNode%ioutputLevel .GE. 1) THEN
          CALL output_line ('Space-Time-Block-FBGS preconditioning of timestep: '//&
              TRIM(sys_siL(isubstep,10))//&
              ', Time: '//TRIM(sys_sdL(rsolverNode%p_rproblem%rtimedependence%dtime,10)))
        END IF
      
        ! -----
        ! Discretise the boundary conditions at the new point in time -- 
        ! if the boundary conditions are nonconstant in time!
        IF (collct_getvalue_int (rsolverNode%p_rproblem%rcollection,'IBOUNDARY') &
            .NE. 0) THEN
          CALL c2d2_updateDiscreteBC (rsolverNode%p_rproblem, .FALSE.)
        END IF

        ! The RHS which is put into the preconditioner is set up in 
        ! rtempVectorRHS to prevent rtempVectorD2 from getting destroyed.
        CALL lsysbl_copyVector (rtempVectorD2,rtempVectorRHS)
        
        ! Is this the first timestep or not?
        IF (isubstep .GT. 0) THEN
        
          ! Create d2 = RHS - Mx1 
          CALL c2d2_setupMatrixWeights (rsolverNode%p_rproblem,p_rspaceTimeMatrix,dtheta,&
            isubstep,-1,rmatrixComponents,iidxNonlin)
            
          CALL c2d2_assembleDefect (rmatrixComponents,rtempVectorX1,&
            rtempVectorRHS,domegaSOR)
          
        END IF

        ! Is this the last timestep or not?
        IF (isubstep .LT. p_rspaceTimeDiscr%NEQtime-1) THEN
          
          ! Read the RHS and solution of the next timestep
          CALL sptivec_getTimestepData (rd, 1+isubstep+1, rtempVectorD3)
          CALL sptivec_getTimestepData (p_rx, 1+isubstep+1, rtempVectorX3)
        
          ! Create d2 = RHS - Ml3 
          CALL c2d2_setupMatrixWeights (rsolverNode%p_rproblem,p_rspaceTimeMatrix,dtheta,&
            isubstep,1,rmatrixComponents,iidxNonlin)
            
          CALL c2d2_assembleDefect (rmatrixComponents,rtempVectorX3,&
            rtempVectorRHS,domegaSOR)
          
        END IF

        ! Read in the solution vector of the current timestep (for nonlinear problems).
        IF (ASSOCIATED(p_rspaceTimeMatrix%p_rsolution)) THEN
          CALL sptivec_getTimestepData (p_rspaceTimeMatrix%p_rsolution, &
              1+isubstep, rtempVectorSol)
        END IF

        ! Set up the matrix weights for the diagonal matrix
        CALL c2d2_setupMatrixWeights (rsolverNode%p_rproblem,p_rspaceTimeMatrix,dtheta,&
          isubstep,0,rmatrixComponents,iidxNonlin)
          
        ! Create d2 = RHS - A(solution) X2
        CALL c2d2_assembleDefect (rmatrixComponents,rtempVectorX2,rtempVectorRHS,&
            1.0_DP,rtempVectorSol)
            
        ! Filter the defect for BC's and initial conditions if necessary
        IF (isubstep .EQ. 0) THEN
          CALL tbc_implementInitCondDefSingle (p_rspaceTimeDiscr, rtempVectorRHS)
        ELSE IF (isubstep .EQ. p_rspaceTimeDiscr%NEQtime-1) THEN
          CALL tbc_implementTermCondDefSingle (p_rspaceTimeDiscr, rtempVectorRHS)
        END IF

        CALL tbc_implementSpatialBCdefect (&
            rsolverNode%p_rproblem,isubstep,dtime,p_rspaceTimeDiscr,rtempVectorRHS)
        
        ! Ok, we have the local defect.
        ! Sum up the norm to the norm of the global vector.
        IF (bcalcNorm) &
          dres = dres + lsysbl_vectorNorm (rtempVectorRHS,LINALG_NORML2)**2
        
        ! Perform preconditioning of the spatial defect with the method 
        ! provided by the core equation module.
        CALL c2d2_precondDefect (&
            rsolverNode%p_rsubnodeBlockFBGS%p_rspatialPreconditioner,&
            rmatrixComponents,&
            rtempVectorRHS,rtempVectorSol,bsuccess,rsolverNode%p_rproblem%rcollection)      
      
        ! Add that defect to the current solution -- damped by domega.
        CALL lsysbl_vectorLinearComb (rtempVectorRHS,rtempVectorX2,&
            rsolverNode%domega,1.0_DP)
      
        ! Save the new solution.
        CALL sptivec_setTimestepData (p_rx, 1+isubstep, rtempVectorX2)
        
        ! Shift the RHS/solution vectors: 1 <- 2 <- 3
        CALL lsysbl_copyVector (rtempVectorD2,rtempVectorD1)
        CALL lsysbl_copyVector (rtempVectorD3,rtempVectorD2)

        CALL lsysbl_copyVector (rtempVectorX2,rtempVectorX1)
        CALL lsysbl_copyVector (rtempVectorX3,rtempVectorX2)

      END DO
      
    END DO ! iiteration
    
    ! Overwrite the rd by our solution.
    CALL sptivec_copyVector (p_rx,rd)
    
    ! Release memory, finish.
    CALL lsysbl_releaseVector (rtempVectorRHS)
    CALL lsysbl_releaseVector (rtempVectorSol)
    CALL lsysbl_releaseVector (rtempVectorD3)
    CALL lsysbl_releaseVector (rtempVectorD2)
    CALL lsysbl_releaseVector (rtempVectorD1)
    CALL lsysbl_releaseVector (rtempVectorX3)
    CALL lsysbl_releaseVector (rtempVectorX2)
    CALL lsysbl_releaseVector (rtempVectorX1)
    
  END SUBROUTINE

! *****************************************************************************
! Routines for the CG solver
! *****************************************************************************

!<subroutine>
  
  SUBROUTINE sptils_initCG (rproblem,p_rsolverNode,p_rpreconditioner)
  
!<description>
  ! Creates a t_linsolNode solver structure for the CG solver. The node
  ! can be used to directly solve a problem or to be attached as solver
  ! or preconditioner to another solver structure. The node can be deleted
  ! by linsol_releaseSolver.
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
  ! A pointer to a t_linsolNode structure. Is set by the routine, any previous
  ! value of the pointer is destroyed.
  TYPE(t_sptilsNode), POINTER         :: p_rsolverNode
!</output>
  
!</subroutine>

    ! Create a default solver structure
    CALL sptils_initSolverGeneral(rproblem,p_rsolverNode)
    
    ! Initialise the type of the solver
    p_rsolverNode%calgorithm = SPTILS_ALG_CG
    
    ! Initialise the ability bitfield with the ability of this solver:
    p_rsolverNode%ccapability = LINSOL_ABIL_CHECKDEF + &
                                LINSOL_ABIL_USESUBSOLVER
    
    ! Allocate the subnode for CG.
    ! This initialises most of the variables with default values appropriate
    ! to this solver.
    ALLOCATE(p_rsolverNode%p_rsubnodeCG)
    
    ! Attach the preconditioner if given. 
    
    IF (PRESENT(p_rpreconditioner)) THEN 
      p_rsolverNode%p_rsubnodeCG%p_rpreconditioner => p_rpreconditioner
    END IF

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>
  
  RECURSIVE SUBROUTINE sptils_setMatrixCG (rsolverNode,Rmatrices)
  
!<description>
  
  ! This routine is called if the system matrix changes.
  ! The routine calls linsol_setMatrices for the preconditioner of CG
  ! to inform also that one about the change of the matrix pointer.
  
!</description>
  
!<input>
  ! An array of system matrices which is simply passed to the initialisation 
  ! routine of the preconditioner.
  TYPE(t_ccoptSpaceTimeMatrix), DIMENSION(:), INTENT(IN)   :: Rmatrices
!</input>
  
!<inputoutput>
  ! The t_linsolNode structure of the CG solver
  TYPE(t_sptilsNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!</subroutine>

    ! Do we have a preconditioner given? If yes, call the general initialisation
    ! routine to initialise it.
    
    IF (ASSOCIATED(rsolverNode%p_rsubnodeCG%p_rpreconditioner)) THEN
      CALL sptils_setMatrices (rsolverNode%p_rsubnodeCG%p_rpreconditioner, &
                               Rmatrices)
    END IF

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>
  
  RECURSIVE SUBROUTINE sptils_initStructureCG (rsolverNode, ierror)
  
!<description>
  ! Calls the initStructure subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_initStructure.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the CG solver
  TYPE(t_sptilsNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the SPTILS_ERR_XXXX constants. A value different to 
  ! SPTILS_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  INTEGER, INTENT(OUT) :: ierror
!</output>

!</subroutine>

    ! local variables
    INTEGER :: i,NEQtime
    TYPE(t_sptilsSubnodeCG), POINTER :: p_rsubnode
    
    ! A-priori we have no error...
    ierror = SPTILS_ERR_NOERROR

    ! We simply pass isubgroup to the subsolvers when calling them.
    ! Inside of this routine, there's not much to do with isubgroup,
    ! as there is no special information (like a factorisation)
    ! associated with this type of solver, which has to be allocated,
    ! released or updated...

    ! Call the init routine of the preconditioner.
    IF (ASSOCIATED(rsolverNode%p_rsubnodeCG%p_rpreconditioner)) THEN
      CALL sptils_initStructure (rsolverNode%p_rsubnodeCG%p_rpreconditioner,ierror)
    END IF
    
    ! CG needs 3 temporary vectors + 1 for preconditioning. 
    ! Allocate that here! Use the default data type prescribed in the solver 
    ! structure for allocating the temp vectors.
    p_rsubnode => rsolverNode%p_rsubnodeCG
    NEQtime = rsolverNode%rmatrix%p_rspaceTimeDiscretisation%NEQtime
    DO i=1,4
      CALL sptivec_initVector (p_rsubnode%RtempVectors(i),NEQtime,&
        rsolverNode%rmatrix%p_rspaceTimeDiscretisation%p_rlevelInfo%p_rdiscretisation)
    END DO
  
    ! Allocate memory for a spatial temp vector.
    CALL lsysbl_createVecBlockByDiscr (&
        rsolverNode%rmatrix%p_rspaceTimeDiscretisation%p_rlevelInfo%p_rdiscretisation,&
        rsolverNode%p_rsubnodeCG%rtempVectorSpace)
    rsolverNode%p_rsubnodeCG%rtempVectorSpace%p_rdiscreteBC => &
        rsolverNode%rmatrix%p_rspaceTimeDiscretisation%p_rlevelInfo%p_rdiscreteBC
    rsolverNode%p_rsubnodeCG%rtempVectorSpace%p_rdiscreteBCfict => &
        rsolverNode%rmatrix%p_rspaceTimeDiscretisation%p_rlevelInfo%p_rdiscreteFBC

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>
  
  RECURSIVE SUBROUTINE sptils_initDataCG (rsolverNode, ierror)
  
!<description>
  ! Calls the initData subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_initData.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the CG solver
  TYPE(t_sptilsNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the SPTILS_ERR_XXXX constants. A value different to 
  ! SPTILS_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  INTEGER, INTENT(OUT) :: ierror
!</output>

!</subroutine>

    ! A-priori we have no error...
    ierror = SPTILS_ERR_NOERROR

    ! We simply pass isubgroup to the subsolvers when calling them.
    ! Inside of this routine, there's not much to do with isubgroup,
    ! as there is no special information (like a factorisation)
    ! associated with this type of solver, which has to be allocated,
    ! released or updated...

    ! Call the init routine of the preconditioner.
    IF (ASSOCIATED(rsolverNode%p_rsubnodeCG%p_rpreconditioner)) THEN
      CALL sptils_initData (rsolverNode%p_rsubnodeCG%p_rpreconditioner, &
                            ierror)
    END IF
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE sptils_doneDataCG (rsolverNode)
  
!<description>
  ! Calls the doneData subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_doneData.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the CG solver
  TYPE(t_sptilsNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!</subroutine>

    ! We simply pass isubgroup to the subsolvers when calling them.
    ! Inside of this routine, there's not much to do with isubgroup,
    ! as there is no special information (like a factorisation)
    ! associated with this type of solver, which has to be allocated,
    ! released or updated...

    ! Call the done routine of the preconditioner.
    IF (ASSOCIATED(rsolverNode%p_rsubnodeCG%p_rpreconditioner)) THEN
      CALL sptils_doneData (rsolverNode%p_rsubnodeCG%p_rpreconditioner)
    END IF
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE sptils_doneStructureCG (rsolverNode)
  
!<description>
  ! Calls the doneStructure subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_doneStructure.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the CG solver
  TYPE(t_sptilsNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!</subroutine>

    ! local variables
    INTEGER :: i
    TYPE(t_sptilsSubnodeCG), POINTER :: p_rsubnode
    
    ! We simply pass isubgroup to the subsolvers when calling them.
    ! Inside of this routine, there's not much to do with isubgroup,
    ! as there is no special information (like a factorisation)
    ! associated with this type of solver, which has to be allocated,
    ! released or updated...

    ! Call the init routine of the preconditioner.
    IF (ASSOCIATED(rsolverNode%p_rsubnodeCG%p_rpreconditioner)) THEN
      CALL sptils_doneStructure (rsolverNode%p_rsubnodeCG%p_rpreconditioner)
    END IF
    
    ! Release temporary data if associated
    p_rsubnode => rsolverNode%p_rsubnodeCG
    IF (p_rsubnode%RtempVectors(1)%NEQtime .NE. 0) THEN
      DO i=4,1,-1
        CALL sptivec_releaseVector (p_rsubnode%RtempVectors(i))
      END DO
    END IF
    
    ! Release the spatial temp vector
    IF (rsolverNode%p_rsubnodeCG%rtempVectorSpace%NEQ .NE. 0) THEN
      CALL lsysbl_releaseVector (rsolverNode%p_rsubnodeCG%rtempVectorSpace)
    END IF
      
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE sptils_doneCG (rsolverNode)
  
!<description>
  ! This routine releases all temporary memory for the CG solver from
  ! the heap. In particular, if a preconditioner is attached to the solver
  ! structure, it's also released from the heap by calling 
  ! linsol_releaseSolver for it.
  ! This DONE routine is declared as RECURSIVE to permit a clean
  ! interaction with linsol_releaseSolver.
!</description>
  
!<input>
  ! A pointer to a t_linsolNode structure of the CG solver.
  TYPE(t_sptilsNode), POINTER         :: rsolverNode
!</input>
  
!</subroutine>

    ! Check if there's a preconditioner attached. If yes, release it.
    IF (ASSOCIATED(rsolverNode%p_rsubnodeCG%p_rpreconditioner)) THEN
      CALL sptils_releaseSolver(rsolverNode%p_rsubnodeCG%p_rpreconditioner)
    END IF
    
    ! Release memory if still associated
    CALL sptils_doneDataCG (rsolverNode)
    CALL sptils_doneStructureCG (rsolverNode)
    
    ! Release the CG subnode
    DEALLOCATE(rsolverNode%p_rsubnodeCG)

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>
  
  RECURSIVE SUBROUTINE sptils_precCG (rsolverNode,rd)
  
!<description>
  ! Applies CG preconditioner $P \approx A$ to the defect 
  ! vector rd and solves $Pd_{new} = d$.
  ! rd will be overwritten by the preconditioned defect.
  !
  ! The matrix must have been attached to the system before calling
  ! this routine, and the initStructure/initData routines
  ! must have been called to prepare the solver for solving
  ! the problem.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the CG solver
  TYPE(t_sptilsNode), INTENT(INOUT), TARGET :: rsolverNode
   
  ! On call to this routine: The defect vector to be preconditioned.
  ! Will be overwritten by the preconditioned defect.
  TYPE(t_spacetimeVector), INTENT(INOUT)        :: rd
!</inputoutput>
  
!</subroutine>

  ! local variables
  REAL(DP) :: dalpha,dbeta,dgamma,dgammaOld,dres,dfr
  INTEGER :: ite

  ! The system matrix
  TYPE(t_ccoptSpaceTimeDiscretisation), POINTER :: p_rspaceTimeDiscr
  TYPE(t_ccoptSpaceTimeMatrix), POINTER :: p_rmatrix
  
  ! Minimum number of iterations, print-sequence for residuals
  INTEGER :: nminIterations, niteResOutput
  
  ! Whether to precondition
  LOGICAL bprec
  
  ! Our structure
  TYPE(t_sptilsSubnodeCG), POINTER :: p_rsubnode
  
  ! Pointers to temporary vectors - named for easier access
  TYPE(t_spacetimeVector), POINTER :: p_DR,p_DP,p_DD,p_rx
  TYPE(t_sptilsNode), POINTER :: p_rprecSubnode
  
  ! DEBUG!!!
  REAL(DP), DIMENSION(:), POINTER :: p_DRdata,p_DPdata,p_DDdata,p_rxdata,p_rddata
  
    ! Solve the system!
  
    ! Status reset
    rsolverNode%iresult = 0
    
    ! Get some information
    p_rsubnode => rsolverNode%p_rsubnodeCG
    p_rspaceTimeDiscr => rsolverNode%rmatrix%p_rspaceTimeDiscretisation
    p_rmatrix => rsolverNode%rmatrix

    ! Check the parameters
    IF ((rd%NEQtime .EQ. 0) .OR. &
        (p_rspaceTimeDiscr%rtimeDiscr%nintervals .EQ. 0) ) THEN
    
      ! Parameters wrong
      rsolverNode%iresult = 2
      RETURN
    END IF

    ! Minimum number of iterations
 
    nminIterations = MAX(rsolverNode%nminIterations,0)
      
    ! Use preconditioning? Filtering?

    bprec = ASSOCIATED(rsolverNode%p_rsubnodeCG%p_rpreconditioner)
    
    ! Iteration when the residuum is printed:

    niteResOutput = MAX(1,rsolverNode%niteResOutput)

    ! Set pointers to the temporary vectors
    ! residuum vector
    p_DR   => p_rsubnode%RtempVectors(1)
    ! preconditioned residuum vector
    p_DP   => p_rsubnode%RtempVectors(2)
    ! direction vector
    p_DD   => p_rsubnode%RtempVectors(3)
    ! solution vector
    p_rx   => p_rsubnode%RtempVectors(4)
    
    IF (bprec) THEN
      p_rprecSubnode => p_rsubnode%p_rpreconditioner
    END IF
    
    ! rd is our RHS. p_rx points to a new vector which will be our
    ! iteration vector. At the end of this routine, we replace
    ! rd by p_rx.
    ! Clear our iteration vector p_rx.
    CALL sptivec_clearVector (p_rx)
      
    ! Initialize used vectors with zero
      
    CALL sptivec_clearVector(p_DR)
    CALL sptivec_clearVector(p_DP)
    CALL sptivec_clearVector(p_DD)
    
    ! Initialization

    dalpha = 1.0_DP
    dbeta  = 1.0_DP
    dgamma = 1.0_DP
    dgammaOld = 1.0_DP

    ! Copy our RHS rd to p_DR. As the iteration vector is 0, this
    ! is also our initial defect.

    CALL sptivec_copyVector(rd,p_DR)
    
    ! Filter the defect for boundary conditions in space and time.
    CALL tbc_implementInitCondDefect (&
        p_rspaceTimeDiscr,p_DR,p_rsubnode%rtempVectorSpace)
    CALL tbc_implementBCdefect (rsolverNode%p_rproblem,&
        p_rspaceTimeDiscr,p_DR,p_rsubnode%rtempVectorSpace)
    
    ! Get the norm of the residuum
    dres = sptivec_vectorNorm (p_DR,rsolverNode%iresNorm)
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
        CALL output_line ('Space-Time-CG: Iteration '// &
             TRIM(sys_siL(0,10))//',  !!RES!! = '//&
             TRIM(sys_sdEL(rsolverNode%dinitialDefect,15)) )
      END IF

      ! Copy the residuum vector p_DR to the preconditioned one.
      CALL sptivec_copyVector(p_DR,p_DP)

      IF (bprec) THEN
        ! Perform preconditioning with the assigned preconditioning
        ! solver structure.
        CALL sptils_precondDefect (p_rprecSubnode,p_DP)
      END IF

      ! Calculate gamma
      dgamma = sptivec_scalarProduct (p_DR, p_DP)
      
      ! Perform at most nmaxIterations loops to get a new vector
      
      ! Copy the preconditioned residual vector to the direction vector
      CALL sptivec_copyVector (p_DP, p_DD)

      DO ite = 1,rsolverNode%nmaxIterations
      
        rsolverNode%icurrentIteration = ite
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! STEP 1: Calculate alpha
        !
        !            < g_k , h_k >           gamma
        ! alpha := ----------------- = -----------------
        !          < d_k , A * d_k >   < d_k , A * d_k >
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! We will now abuse the p_DP vector to save the result of
        ! the matrix-vector-multiplication A * d_k.
        ! Using this trick, we can save one temporary vector.
        CALL c2d2_spaceTimeMatVec (rsolverNode%p_rproblem, p_rmatrix, &
            p_DD,p_DP, 1.0_DP,0.0_DP,SPTID_FILTER_NONE)
        
        ! Calculate alpha for the CG iteration
        dalpha = dgamma / sptivec_scalarProduct (p_DD, p_DP)

        IF (dalpha .EQ. 0.0_DP) THEN
          ! We are below machine exactness - we can't do anything more...
          ! May happen with very small problems with very few unknowns!
          IF (rsolverNode%ioutputLevel .GE. 2) THEN
            CALL output_line('Space-Time-CG: Convergence failed, ALPHA=0!')
            rsolverNode%iresult = -2
            EXIT
          END IF
        END IF
        
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! STEP 2: Calculate x_{k+1}
        !
        ! x_{k+1} := x_k + alpha * d_k
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        CALL sptivec_vectorLinearComb (p_DD, p_rx, dalpha, 1.0_DP)

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! STEP 3: Calculate g_{k+1}
        !
        ! g_{k+1} := g_k - alpha * A * d_k
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Since we have abused p_DP for saving A * d_k, we can use it now
        ! for the linear combination of g_{k+1}
        CALL sptivec_vectorLinearComb (p_DP, p_DR, -dalpha, 1.0_DP)

        ! Filter the defect for boundary conditions in space and time.
        CALL tbc_implementInitCondDefect (&
            p_rspaceTimeDiscr,p_DR,p_rsubnode%rtempVectorSpace)
        CALL tbc_implementBCdefect (rsolverNode%p_rproblem,&
           p_rspaceTimeDiscr,p_DR,p_rsubnode%rtempVectorSpace)

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! STEP 4: Calculate ||g_{k+1}|| and write some output
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Get the norm of the new (final?) residuum
        dfr = sptivec_vectorNorm (p_DR,rsolverNode%iresNorm)
     
        rsolverNode%dfinalDefect = dfr

        ! Test if the iteration is diverged
        IF (sptils_testDivergence(rsolverNode,dfr)) THEN
          CALL output_line('Space-Time-CG: Solution diverging!')
          rsolverNode%iresult = 1
          EXIT
        END IF
     
        ! At least perform nminIterations iterations
        IF (ite .GE. nminIterations) THEN
        
          ! Check if the iteration converged
          IF (sptils_testConvergence(rsolverNode,dfr)) EXIT
          
        END IF

        ! print out the current residuum

        IF ((rsolverNode%ioutputLevel .GE. 2) .AND. &
            (MOD(ite,niteResOutput).EQ.0)) THEN
          CALL output_line ('Space-Time-CG: Iteration '// &
              TRIM(sys_siL(ITE,10))//',  !!RES!! = '//&
              TRIM(sys_sdEL(rsolverNode%dfinalDefect,15)) )
        END IF
        
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! STEP 5: Calculate h_{k+1}
        !
        ! In other words: Copy the residual vector and apply the
        ! preconditioner, if given.
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Copy the residuum vector p_DR to the preconditioned one.
        CALL sptivec_copyVector(p_DR,p_DP)
    
        IF (bprec) THEN
          ! Perform preconditioning with the assigned preconditioning
          ! solver structure.
          CALL sptils_precondDefect (p_rprecSubnode,p_DP)
        END IF
        
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! STEP 6: Calculate gamma
        !
        ! gamma' := gamma
        ! gamma  := < g_{k+1}, h_{k+1} >
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        dgammaOld = dgamma
        dgamma = sptivec_scalarProduct (p_DR, p_DP)
        
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! STEP 7: Calculate beta
        !
        !          gamma
        ! beta := --------
        !          gamma'
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        dbeta = dgamma / dgammaOld

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! STEP 8: Calculate d_{k+1}
        !
        ! d_{k+1} := beta * d_k + h_{k+1}
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        CALL sptivec_vectorLinearComb (p_DP, p_DD, 1.0_DP, dbeta)
        
        ! That's it - next iteration!
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
        CALL output_line ('Space-Time-CG: Iteration '// &
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
        CALL output_line ('Space-Time-CG statistics:')
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
              'Space-Time-CG: Iterations/Rate of convergence: '//&
              TRIM(sys_siL(rsolverNode%iiterations,10))//' /'//&
              TRIM(sys_sdEL(rsolverNode%dconvergenceRate,15)) )
      END IF

    ELSE
      ! DEF=Infinity; RHO=Infinity, set to 1
      rsolverNode%dconvergenceRate = 1.0_DP
    END IF  
  
  END SUBROUTINE

! *****************************************************************************
! Routines for the UMFPACK4 preconditioner
! *****************************************************************************

!<subroutine>
  
  SUBROUTINE sptils_initUMFPACK4 (rproblem,rsolverNode)
  
!<description>
  
  ! Creates a t_sptilsNode solver structure for the UMFPACK4 solver. The node
  ! can be used to directly solve a problem or to be attached as solver
  ! or preconditioner to another solver structure. The node can be deleted
  ! by sptils_releaseSolver.
  
!</description>
  
!<input>
  ! The problem structure that defines the problem.
  ! A pointer to this is saved in the solver node, so the problem structure
  ! must not be released before the solver node is released.
  TYPE(t_problem), TARGET :: rproblem
!</input>
  
!<output>
  
  ! A pointer to a t_sptilsNode structure. Is set by the routine, any previous
  ! value of the pointer is destroyed.
  TYPE(t_sptilsNode), POINTER         :: rsolverNode
   
!</output>
  
!</subroutine>

  ! Create a default solver structure
  
  CALL sptils_initSolverGeneral(rproblem,rsolverNode)
  
  ! Initialise the type of the solver
  rsolverNode%calgorithm = SPTILS_ALG_UMFPACK4
  
  ! Initialise the ability bitfield with the ability of this solver:
  rsolverNode%ccapability = SPTILS_ABIL_DIRECT
  
  ! Allocate the subnode for UMFPACK4.
  ! This initialises most of the variables with default values appropriate
  ! to this solver.
  ALLOCATE(rsolverNode%p_rsubnodeUMFPACK4)
  
  ! Initialize the UMFPACK4 control structure:
  CALL UMF4DEF(rsolverNode%p_rsubnodeUMFPACK4%Dcontrol)
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE sptils_initStructureUMFPACK4 (rsolverNode,ierror)
  
!<description>
  ! Performs a symbolic factorisation on the assigned matrix.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure of the UMFPACK4 solver
  TYPE(t_sptilsNode), INTENT(INOUT),TARGET  :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the SPTILS_ERR_XXXX constants. A value different to 
  ! SPTILS_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  INTEGER, INTENT(OUT) :: ierror
!</output>
  
!</subroutine>

  ! A-priori we have no error...
  ierror = SPTILS_ERR_NOERROR

  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE sptils_initDataUMFPACK4 (rsolverNode, ierror)
  
!<description>
  ! Performs a numeric factorisation on the assigned matrix.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure of the UMFPACK4 solver
  TYPE(t_sptilsNode), INTENT(INOUT), TARGET :: rsolverNode
!</inputoutput>

!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to 
  ! SPTILS_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  INTEGER, INTENT(OUT) :: ierror
!</output>
  
!</subroutine>

    ! local variables
    TYPE(t_matrixScalar), POINTER :: p_rmatrix
    TYPE(t_matrixBlock), TARGET :: rmatrixGlobal
    TYPE(t_matrixBlock) :: rtempMatrix
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_Kld
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Kcol
    REAL(DP), DIMENSION(:), POINTER :: p_DA

    ! Status variables of UMFPACK4; receives the UMFPACK-specific return code
    ! of a call to the solver routines.
    REAL(DP), DIMENSION(90) :: Dinfo
    
    ! A-priori we have no error...
    ierror = SPTILS_ERR_NOERROR

    ! We have to create a global matrix first!
    ! Call the assembly routine to assemble the global block matrix.
    ! Afterwards, reshape the data to form a scalar matrix which
    ! can be feed to UMFPACK.
    CALL assembleGlobalSpaceTimeMatrix (rsolverNode%p_rproblem,&
        rsolverNode%rmatrix,rmatrixGlobal)
    CALL glsys_assembleGlobal (rmatrixGlobal, rtempMatrix, .TRUE., .TRUE.)
    
    ! The global block matrix is not needed anymore.
    CALL lsysbl_releaseMatrix (rmatrixGlobal)
    
    ! Implement Dirichlet boundary conditions into the pressure matrix.
    ! This is necessary if we have pure Dirichlet BC's in the velocity,
    ! so the pressure is not uniquely defined. In this case, we fix
    ! the pressure in the first node by writing a unit line into
    ! the global matrix for the first pressure DOF in every timestep.
    CALL pressureDirichlet (rsolverNode%p_rproblem,&
        rsolverNode%rmatrix,rtempMatrix)

    !CALL matio_writeBlockMatrixHR (rtempMatrix, 'matrix',&
    !                               .TRUE., 0, 'matrix.txt', '(E15.5)')
    !CALL matio_writeBlockMatrixMaple (rtempMatrix, 'A',&
    !                             0, 'matrix.txt', '(E20.10)')

    ! Now start to modify the temp matrix for UMFPACK's needs.  
    p_rmatrix => rtempMatrix%RmatrixBlock(1,1)
    
    !CALL matio_spyMatrix('matrix','matrix',p_rmatrix,.TRUE.)

    !!! DEBUG
    !CALL lsyssc_getbase_double (rtempMatrix,p_DA)
    !WHERE (abs(p_Da) .LT. 1.0E-12_DP) p_Da = 0.0_DP
    !CALL matio_writeMatrixHR (p_rmatrix, 'matrix',&
    !                          .TRUE., 0, 'matrix.txt', '(D20.10)')
    !CALL sys_halt()

    ! Modify Kcol/Kld of the matrix. Subtract 1 to get them 0-based.
    CALL lsyssc_addIndex (p_rmatrix%h_Kcol,-1_I32)
    CALL lsyssc_addIndex (p_rmatrix%h_Kld,-1_I32)
    
    ! Get the data arrays.
    CALL lsyssc_getbase_Kcol (p_rmatrix,p_Kcol)
    CALL lsyssc_getbase_Kld (p_rmatrix,p_Kld)
    CALL lsyssc_getbase_double (p_rmatrix,p_DA)
    
    ! Perform a symbolic factorization...
    CALL UMF4SYM(rtempMatrix%NEQ,rtempMatrix%NCOLS,p_Kld,p_Kcol,p_Da, &
                rsolverNode%p_rsubnodeUMFPACK4%isymbolic,&
                rsolverNode%p_rsubnodeUMFPACK4%Dcontrol,&
                Dinfo)
                 
    ! Check Dinfo(1) if there is an error
    SELECT CASE (INT(Dinfo(1)))
    CASE (0)
    
      ! Perform a numeric factorization...
      CALL UMF4NUM(p_Kld,p_Kcol,p_Da, &
              rsolverNode%p_rsubnodeUMFPACK4%isymbolic,&
              rsolverNode%p_rsubnodeUMFPACK4%inumeric,&
              rsolverNode%p_rsubnodeUMFPACK4%Dcontrol,&
              Dinfo)
              
      ! Check Dinfo(1) if there is an error
      SELECT CASE (INT(Dinfo(1)))
      CASE (0)
        ! ok.
      CASE (1)
        ! Singular matrix
        ierror = LINSOL_ERR_SINGULAR
      CASE (-1)
        ! no memory
        ierror = LINSOL_ERR_NOMEMORY
      CASE (-11)
        ! no memory
        ierror = LINSOL_ERR_MATRIXHASCHANGED
      CASE DEFAULT
        ! don't know what went wrong
        ierror = LINSOL_ERR_INITERROR
      END SELECT
    
    CASE (1)
      ! Singular matrix
      ierror = LINSOL_ERR_SINGULAR
      RETURN
    CASE (-1)
      ! no memory
      ierror = LINSOL_ERR_NOMEMORY
    CASE DEFAULT
      ! don't know what went wrong
      ierror = LINSOL_ERR_INITERROR
    END SELECT

    ! Throw away the temporary matrices
    CALL lsysbl_releaseMatrix (rtempMatrix)

  CONTAINS

    ! ---------------------------------------------------------------

    SUBROUTINE pressureDirichlet (rproblem,rsupermatrix,rmatrix)
    
    ! If we have a Pure-Dirichlet problem, this replaces the row of the first
    ! pressure DOF in every timestep by zero.
    
    ! The problem structure that defines the problem.
    TYPE(t_problem), INTENT(INOUT) :: rproblem

    ! The source space-time matrix where rmatrix was generated from
    TYPE(t_ccoptSpaceTimeMatrix), INTENT(IN) :: rsupermatrix

    ! The global space time matrix to be modified.
    ! Must be a 1x1 block matrix with the submatrix (1,1) representing
    ! the global space time matrix.
    TYPE(t_matrixBlock), INTENT(INOUT) :: rmatrix
    
      INTEGER(PREC_VECIDX), DIMENSION(:), ALLOCATABLE :: Iidx
      INTEGER(PREC_VECIDX) :: ivelSize,ipSize,neq
      INTEGER :: i
      TYPE(t_ccoptSpaceTimeDiscretisation), POINTER :: p_rspaceTimeDiscr
    
      p_rspaceTimeDiscr => rsupermatrix%p_rspaceTimeDiscretisation
    
      ! Nothing to do if there is Nuemann boundary here.
      IF (p_rspaceTimeDiscr%p_rlevelInfo%bhasNeumannBoundary) RETURN
    
      neq = p_rspaceTimeDiscr%p_rlevelInfo%rpreallocatedSystemMatrix%neq
      ivelSize = 2 * p_rspaceTimeDiscr%p_rlevelInfo%rmatrixStokes%NEQ
      ipSize = p_rspaceTimeDiscr%p_rlevelInfo%rmatrixB1%NCOLS
    
      ! Create an array containing all the rows where a unit vector is to be
      ! imposed. Size = 2 * number of timesteps in the supermatrix - 1
      ! (primal + dual pressure, not the initial condition)
      ALLOCATE(Iidx(0:2*(p_rspaceTimeDiscr%NEQtime-1)))
      
      ! Add the row of the first primal/dual pressure DOF to that array.
      Iidx(0) = neq-ipSize+1  ! 0th time step -> only dual pressure because of init. cond.
      DO i=1,p_rspaceTimeDiscr%NEQtime-1
        Iidx(2*i-1) = i*neq+ivelSize+1
        Iidx(2*i  ) = i*neq+neq-ipSize+1
      END DO
      
      ! Filter the matrix, implement the unit vectors.
      CALL mmod_replaceLinesByUnit (rmatrix%RmatrixBlock(1,1),Iidx)
      
      ! Release memory, finish.
      DEALLOCATE(Iidx)
    
    END SUBROUTINE

  END SUBROUTINE

  ! ***************************************************************************

  SUBROUTINE assembleGlobalSpaceTimeMatrix (rproblem,rsupermatrix,rmatrix)
  
  ! Assembles a block matrix rmatrix from a space-time matrix rsupermatrix
  ! by plugging together all blocks.
  ! Implement boundarys conditions into rmatrix during the assembly.
  ! WARNING: SPACE CONSUMING!!!
  
  ! The problem structure that defines the problem.
  TYPE(t_problem), INTENT(INOUT) :: rproblem

  ! The source space-time matrix
  TYPE(t_ccoptSpaceTimeMatrix), INTENT(IN), TARGET :: rsupermatrix
  
  ! The destination block matrix
  TYPE(t_matrixBlock), INTENT(OUT) :: rmatrix

    ! local variables  
    TYPE(t_matrixBlock) :: rblockTemp
    TYPE(t_ccmatrixComponents) :: rmatrixComponents
    TYPE(t_vectorBlock) :: rvector
    INTEGER :: isubstep,ileft,iright,ix,iidxNonlin
    INTEGER(PREC_VECIDX), DIMENSION(1) :: Irows
    INTEGER(PREC_VECIDX) :: idiag
    TYPE(t_ccoptSpaceTimeDiscretisation), POINTER :: p_rspaceTimeDiscr
  
    p_rspaceTimeDiscr => rsupermatrix%p_rspaceTimeDiscretisation
   
    ! Create a global matrix:
    CALL lsysbl_createEmptyMatrix (rmatrix,6*(p_rspaceTimeDiscr%NEQtime))
  
    ! Basic initialisation of rmatrixComponents with the pointers to the
    ! matrices / discretisation structures on the current level.
    !
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
    rmatrixComponents%dnu = collct_getvalue_real (rproblem%rcollection,'NU')
    rmatrixComponents%iupwind1 = collct_getvalue_int (rproblem%rcollection,'IUPWIND1')
    rmatrixComponents%dupsam1 = collct_getvalue_real (rproblem%rcollection,'UPSAM1')
    rmatrixComponents%iupwind2 = collct_getvalue_int (rproblem%rcollection,'IUPWIND2')
    rmatrixComponents%dupsam2 = collct_getvalue_real (rproblem%rcollection,'UPSAM2')

    ! Copy references to the preallocated system matrix. Use that as space
    ! for storing matrix data.
    CALL lsysbl_duplicateMatrix (p_rspaceTimeDiscr%p_rlevelInfo%rpreallocatedSystemMatrix,&
        rblockTemp,LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
        
    ! Create a vector for evaluating the nonlinearity.
    ! (The vector is always created but only used when there is a nonlinearity!)
    CALL lsysbl_createVecBlockIndMat(rblockTemp,rvector)
    
    ! Loop through the substeps
    DO isubstep = 0,p_rspaceTimeDiscr%NEQtime-1
    
      ! Current point in time
      rproblem%rtimedependence%dtime = &
          rproblem%rtimedependence%dtimeInit + isubstep*p_rspaceTimeDiscr%rtimeDiscr%dtstep
      rproblem%rtimedependence%itimestep = isubstep

      ! -----
      ! Discretise the boundary conditions at the new point in time -- 
      ! if the boundary conditions are nonconstant in time!
      IF (collct_getvalue_int (rproblem%rcollection,'IBOUNDARY') .NE. 0) THEN
        CALL c2d2_updateDiscreteBC (rproblem, .FALSE.)
      END IF
      
      ! Assemble diagonal blocks as well as the band above and below the diagonal.
      ileft = -1
      iright = 1
      IF (isubstep .EQ. 0) ileft = 0
      IF (isubstep .EQ. p_rspaceTimeDiscr%NEQtime-1) iright = 0
      
      ! Loop over the matrix bands in the current row isubstep
      DO ix = ileft,iright
      
        ! Set up the matrix weights of that submatrix.
        CALL c2d2_setupMatrixWeights (rproblem,rsupermatrix,&
          p_rspaceTimeDiscr%rtimeDiscr%dtheta,isubstep,ix,rmatrixComponents,iidxNonlin)
          
        ! If there is a nonlinearity involved, get the evaluation point.
        IF (iidxNonlin .GT. 0)  THEN
          CALL sptivec_getTimestepData (rsupermatrix%p_rsolution, 1+iidxNonlin-1, rvector)
        END IF
      
        ! Assemble the matrix in rblockTemp.
        ! Note that the weights of the matrices must be set before, otherwise
        ! the assembly routines would complain about missing matrices :-)
        rblockTemp%RmatrixBlock(1,1)%dscaleFactor = 1.0_DP
        rblockTemp%RmatrixBlock(2,2)%dscaleFactor = 1.0_DP
        rblockTemp%RmatrixBlock(4,4)%dscaleFactor = 1.0_DP
        rblockTemp%RmatrixBlock(5,5)%dscaleFactor = 1.0_DP
        CALL c2d2_assembleMatrix (CCMASM_COMPUTE,CCMASM_MTP_AUTOMATIC,&
            rblockTemp,rmatrixComponents,rvector,ctypePrimalDual=0) 

        ! Switch of matrices that aren't needed.
        SELECT CASE (ix)
        CASE (-1)
          ! Specify the matrix as 'off-diagonal' matrix because it's not on the
          ! main diagonal of the supermatrix.
          rblockTemp%imatrixSpec = LSYSBS_MSPEC_OFFDIAGSUBMATRIX
          
          ! Mass matrices in the dual equation not needed since the band
          ! below the diagonal couples only the timesteps of the primal equation.
          ! This saves some memory as the matrices are =0 anyway.
          rblockTemp%RmatrixBlock(1,1)%dscaleFactor = 1.0_DP
          rblockTemp%RmatrixBlock(2,2)%dscaleFactor = 1.0_DP
          rblockTemp%RmatrixBlock(4,4)%dscaleFactor = 0.0_DP
          rblockTemp%RmatrixBlock(5,5)%dscaleFactor = 0.0_DP
        CASE (1)
          ! Specify the matrix as 'off-diagonal' matrix because it's not on the
          ! main diagonal of the supermatrix.
          rblockTemp%imatrixSpec = LSYSBS_MSPEC_OFFDIAGSUBMATRIX
          
          ! Mass matrices in the primal equation not needed since the band
          ! above the diagonal couples only the timesteps of the dual equation.
          ! This saves some memory as the matrices are =0 anyway.
          rblockTemp%RmatrixBlock(1,1)%dscaleFactor = 0.0_DP
          rblockTemp%RmatrixBlock(2,2)%dscaleFactor = 0.0_DP
          rblockTemp%RmatrixBlock(4,4)%dscaleFactor = 1.0_DP
          rblockTemp%RmatrixBlock(5,5)%dscaleFactor = 1.0_DP
        CASE DEFAULT
          rblockTemp%imatrixSpec = LSYSBS_MSPEC_GENERAL
          
          rblockTemp%RmatrixBlock(1,1)%dscaleFactor = 1.0_DP
          rblockTemp%RmatrixBlock(2,2)%dscaleFactor = 1.0_DP
          rblockTemp%RmatrixBlock(4,4)%dscaleFactor = 1.0_DP
          rblockTemp%RmatrixBlock(5,5)%dscaleFactor = 1.0_DP
        END SELECT

        ! Include the boundary conditions into that matrix.
        CALL matfil_discreteBC (rblockTemp,p_rspaceTimeDiscr%p_rlevelInfo%p_rdiscreteBC)
        CALL matfil_discreteFBC (rblockTemp,p_rspaceTimeDiscr%p_rlevelInfo%p_rdiscreteFBC)

        ! Include the current matrix into the global matrix 
        CALL insertMatrix (rblockTemp,rmatrix,(isubstep+ix)*6+1,isubstep*6+1,.FALSE.)
        
        ! DEBUG!!!
        !IF (isubstep .EQ. p_rspaceTimeDiscr%niterations) THEN
        !  CALL matio_writeBlockMatrixHR (rblockTemp, 'matrix',&
        !    .TRUE., 0, 'matrix'//TRIM(sys_siL(ix,10))//'.txt', '(E20.10)')
        !END IF
        
      END DO  

      IF (.NOT. p_rspaceTimeDiscr%p_rlevelInfo%bhasNeumannBoundary) THEN
        ! Insert a 'point matrix' containing a zero in the pressure block.
        ! This allows the boundary condition implementation routine to
        ! insert a unit vector for the pressure if we have a pure-Dirichlet
        ! problem.
        ! Don't insert the matrix if there is already one! An existing matrix
        ! is always an identity matrix, as it happens to appear in the first
        ! and last time strep.
        ! IF (isubstep .NE. 0) THEN
        IF (.NOT. lsysbl_isSubmatrixPresent (rmatrix,isubstep*6+3,isubstep*6+3)) THEN
          CALL createPointMatrix (rmatrix%RmatrixBlock(isubstep*6+3,isubstep*6+3),&
              p_rspaceTimeDiscr%p_rlevelInfo%rmatrixIdentityPressure%NEQ,1)
        END IF
        IF (.NOT. lsysbl_isSubmatrixPresent (rmatrix,isubstep*6+6,isubstep*6+6)) THEN
          CALL createPointMatrix (rmatrix%RmatrixBlock(isubstep*6+6,isubstep*6+6),&
              p_rspaceTimeDiscr%p_rlevelInfo%rmatrixIdentityPressure%NEQ,1)
        END IF
        
      END IF
      
    END DO
    
    ! Release the temp vector
    CALL lsysbl_releaseVector (rvector)
    
    ! Release the temp matrix
    CALL lsysbl_releaseMatrix (rblockTemp)

    ! Update structural information of the global matrix.
    CALL lsysbl_updateMatStrucInfo (rmatrix)
  
  CONTAINS

    ! ---------------------------------------------------------------
    
    SUBROUTINE insertMatrix (rsource,rdest,ileft,itop,bcopyStructure)
    
    ! Includes rsource into rdest at position ileft,itop
    TYPE(t_matrixBlock), INTENT(IN) :: rsource
    TYPE(t_matrixBlock), INTENT(INOUT) :: rdest
    INTEGER, INTENT(IN) :: ileft
    INTEGER, INTENT(IN) :: itop
    ! Duplicate the structure. FALSE=share the structure.
    LOGICAL, INTENT(IN) :: bcopyStructure 
    
    INTEGER :: i,j,ccopy
    
    ccopy = LSYSSC_DUP_SHARE
    IF (bcopyStructure) ccopy = LSYSSC_DUP_COPY
    
    DO j=1,rsource%ndiagBlocks
      DO i=1,rsource%ndiagBlocks
        IF (lsysbl_isSubmatrixPresent (rsource,i,j)) THEN
          IF (lsysbl_isSubmatrixPresent (rdest,i+itop-1,j+ileft-1)) THEN
            CALL lsyssc_releaseMatrix (rdest%RmatrixBlock(j+ileft-1,i+itop-1))
          END IF
          CALL lsyssc_duplicateMatrix (rsource%RmatrixBlock(i,j),&
              rdest%RmatrixBlock(i+itop-1,j+ileft-1),&
              ccopy,LSYSSC_DUP_COPY)
        END IF
      END DO
    END DO
        
    END SUBROUTINE

    ! ---------------------------------------------------------------
    
    SUBROUTINE createSumMatrix (rmatrix,NCOLS,NEQ,irow)
    
    ! Creates a structure-9 matrix that contains one row, filled by 'ones'.
    ! This corresponds to an additional equation that sums up all
    ! vector entries, which corresponds to setting up the integral
    ! of a function.
    ! If the matrix does not exist, it's created.
    ! If the matrix exists, it's overwritten.
    
    ! Matrix to be set up
    TYPE(t_matrixScalar), INTENT(INOUT) :: rmatrix
    
    ! Number of rows/columns in the matrix
    INTEGER(PREC_VECIDX), INTENT(IN) :: NCOLS
    INTEGER(PREC_VECIDX), INTENT(IN) :: NEQ
    
    ! Number of the row that should contain the line '1...1'.
    INTEGER(PREC_VECIDX), INTENT(IN) :: irow
    
      ! local variables
      INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_Kld
      INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_Kdiagonal
      REAL(DP), DIMENSION(:), POINTER :: p_Ddata
      INTEGER(PREC_VECIDX) :: i
      
      ! The matrix should get the shape:
      !   1 1 1 1 ... 1 1 1
      !   0 0 0 0 ... 0 0 0      
      !   ...
      !   0 0 0 0 ... 0 0 0      
      !
      ! Create the matrix by hand. If necessary, allocate memory.
      IF (rmatrix%h_Kld .EQ. ST_NOHANDLE) THEN
        CALL storage_new ('createSumMatrix', 'Kld', NEQ+1, &
                          ST_INT, rmatrix%h_Kld,ST_NEWBLOCK_NOINIT)
        CALL storage_new ('createSumMatrix', 'Kdiagonal', NEQ+1, &
                          ST_INT, rmatrix%h_Kdiagonal,ST_NEWBLOCK_NOINIT)
        CALL storage_new ('createSumMatrix', 'Kcol', NCOLS, &
                          ST_INT, rmatrix%h_Kcol,ST_NEWBLOCK_NOINIT)
                          
        ! DA has probably a zero in front and at the end as dummy
        ! for diagonal entries.
        CALL storage_new ('createSumMatrix', 'Da', NCOLS, &
                          ST_DOUBLE, rmatrix%h_Da,ST_NEWBLOCK_NOINIT)

        ! Initialise matrix format and other parameters, that's it.
        rmatrix%NEQ = NEQ
        rmatrix%NCOLS = NCOLS
        rmatrix%cmatrixFormat = LSYSSC_MATRIX9
        rmatrix%NA = NCOLS
      
      END IF
      
      CALL lsyssc_getbase_double (rmatrix,p_Ddata)
      CALL lsyssc_getbase_Kcol (rmatrix,p_Kcol)
      CALL lsyssc_getbase_Kld (rmatrix,p_Kld)
      CALL lsyssc_getbase_Kdiagonal (rmatrix,p_Kdiagonal)
      
      ! irow'th row
      p_Ddata(:) = 1.0_DP
      
      p_Kld(1:irow) = 1
      p_Kld(irow+1:) = NCOLS+1
      
      p_Kdiagonal(1:irow-1) = 1
      p_Kdiagonal(irow) = irow
      
      DO i=1,NCOLS
        p_Kcol(i) = i
      END DO

    END SUBROUTINE

    ! ---------------------------------------------------------------
    
    SUBROUTINE createPointMatrix (rmatrix,NEQ,irow)
    
    ! Creates a structure-9 matrix that contains exactly one 'zero'
    ! on the diagonal of row irow. 
    
    ! Matrix to be set up
    TYPE(t_matrixScalar), INTENT(INOUT) :: rmatrix
    
    ! Number of rows/columns in the matrix
    INTEGER(PREC_VECIDX), INTENT(IN) :: NEQ
    
    ! Number of the row that should contain the 'one'
    INTEGER(PREC_VECIDX), INTENT(IN) :: irow
    
      ! local variables
      INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Kcol
      INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_Kld
      INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_Kdiagonal
      REAL(DP), DIMENSION(:), POINTER :: p_Ddata
      INTEGER(PREC_VECIDX) :: i
      
      ! Create the matrix by hand. If necessary, allocate memory.
      IF (rmatrix%h_Kld .EQ. ST_NOHANDLE) THEN
        CALL storage_new ('createPointMatrix', 'Kld', NEQ+1, &
                          ST_INT, rmatrix%h_Kld,ST_NEWBLOCK_NOINIT)
        CALL storage_new ('createPointMatrix', 'Kdiagonal', NEQ, &
                          ST_INT, rmatrix%h_Kdiagonal,ST_NEWBLOCK_NOINIT)
        CALL storage_new ('createPointMatrix', 'Kcol', 1, &
                          ST_INT, rmatrix%h_Kcol,ST_NEWBLOCK_NOINIT)
                          
        ! DA has probably a zero in front and at the end as dummy
        ! for diagonal entries.
        CALL storage_new ('createPointMatrix', 'Da', 1, &
                          ST_DOUBLE, rmatrix%h_Da,ST_NEWBLOCK_NOINIT)

        ! Initialise matrix format and other parameters, that's it.
        rmatrix%NEQ = NEQ
        rmatrix%NCOLS = NEQ
        rmatrix%cmatrixFormat = LSYSSC_MATRIX9
        rmatrix%NA = 1
        rmatrix%dscaleFactor = 1.0_DP
      
      END IF
      
      CALL lsyssc_getbase_double (rmatrix,p_Ddata)
      CALL lsyssc_getbase_Kcol (rmatrix,p_Kcol)
      CALL lsyssc_getbase_Kld (rmatrix,p_Kld)
      CALL lsyssc_getbase_Kdiagonal (rmatrix,p_Kdiagonal)
      
      ! irow'th row
      p_Ddata(1) = 0.0_DP
      p_Kcol(1) = irow

      p_Kld(1:irow) = 1
      p_Kld(irow+1:) = 2
      
      p_Kdiagonal(:) = 1
      
    END SUBROUTINE
    
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE sptils_doneDataUMFPACK4 (rsolverNode)
  
!<description>
  ! Releases the memory of the numeric factorisation of the given matrix.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure of the UMFPACK4 solver
  TYPE(t_sptilsNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!</subroutine>

  ! Release the numerical factorisation if associated
  IF (rsolverNode%p_rsubnodeUMFPACK4%inumeric .NE. 0) THEN
    CALL UMF4FNUM(rsolverNode%p_rsubnodeUMFPACK4%inumeric)
    rsolverNode%p_rsubnodeUMFPACK4%inumeric = 0
  END IF  
  
  ! Release the symbolical factorisation if associated
  IF (rsolverNode%p_rsubnodeUMFPACK4%isymbolic .NE. 0) THEN
    CALL UMF4FSYM(rsolverNode%p_rsubnodeUMFPACK4%isymbolic)
    rsolverNode%p_rsubnodeUMFPACK4%isymbolic = 0
  END IF
  
END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE sptils_doneStructureUMFPACK4 (rsolverNode)
  
!<description>
  ! Releases the memory of the numeric factorisation of the given matrix.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the UMFPACK4 solver
  TYPE(t_sptilsNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!</subroutine>

  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE sptils_doneUMFPACK4 (rsolverNode)
  
!<description>
  ! This routine releases all temporary memory for the UMFPACK4 solver from
  ! the heap.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure of UMFPACK4 which is to be cleaned up.
  TYPE(t_sptilsNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!</subroutine>
  
  ! Release symbolical and numerical factorisation if still associated...
  CALL sptils_doneDataUMFPACK4 (rsolverNode)
  CALL sptils_doneStructureUMFPACK4 (rsolverNode)
  
  DEALLOCATE(rsolverNode%p_rsubnodeUMFPACK4)
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE sptils_precUMFPACK4 (rsolverNode,rd)
  
!<description>
  ! Applies UMFPACK preconditioner $P \approx A$ to the defect 
  ! vector rd and solves $Pd_{new} = d$.
  ! rd will be overwritten by the preconditioned defect.
  !
  ! The matrix must have been attached to the system before calling
  ! this routine, and the initStructure/initData routines
  ! must have been called to prepare the solver for solving
  ! the problem.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure of the UMFPACK4 solver
  TYPE(t_sptilsNode), INTENT(INOUT)         :: rsolverNode

  ! On call to this routine: The defect vector to be preconditioned.
  ! Will be overwritten by the preconditioned defect.
  TYPE(t_spacetimeVector), INTENT(INOUT)        :: rd
!</inputoutput>
  
!</subroutine>
 
  ! local variables
  INTEGER :: KSYS
  REAL(DP), DIMENSION(:), POINTER :: p_Dx,p_Db
  TYPE(t_vectorBlock) :: rb,rx
  TYPE(t_sptilsSubnodeUMFPACK4), POINTER :: p_rsubnode

  ! Status variables of UMFPACK4; receives the UMFPACK-specific return code
  ! of a call to the solver routines.
  REAL(DP), DIMENSION(90) :: Dinfo
  
    ! Status reset
    rsolverNode%iresult = 0
    
    ! Getch some information
    p_rsubnode => rsolverNode%p_rsubnodeUMFPACK4
    
    ! Copy the RHS rd to the temp vector; it will be overwritten
    ! by the solution vector
    CALL sptivec_convertSupervecToVector (rd,rb)

    ! Get the array    
    CALL lsysbl_getbase_double(rb,p_Db)

    ! Allocate memory for the destination
    CALL sptivec_convertSupervecToVector (rd,rx)

    ! Get the RHS and solution vector data
    CALL lsysbl_getbase_double(rx,p_Dx)

    ! Solve the system
    ! Solve the system. Note that UMFPACK expects the matrix in
    ! CSR format, which is transposed to our matrix format 9 --
    ! So we solve the transposed system:
    KSYS = 1
    CALL UMF4SOL(KSYS,p_Dx,p_Db,rsolverNode%p_rsubnodeUMFPACK4%inumeric,&
                 rsolverNode%p_rsubnodeUMFPACK4%Dcontrol,Dinfo)

    ! Save the result
    CALL sptivec_convertVectorToSupervec (rx,rd)

    ! Scale by omega
    CALL sptivec_scaleVector (rd,rsolverNode%domega)
                
    ! Release temp vectors
    CALL lsysbl_releaseVector (rx)
    CALL lsysbl_releaseVector (rb)
  
    ! Check the solver status
    SELECT CASE (INT(Dinfo(1)))
    CASE (0) 
      ! All ok.
      rsolverNode%iiterations = 1
      rsolverNode%dconvergenceRate = 0.0_DP
    CASE DEFAULT
      ! We had an error. Don't know which one.
      rsolverNode%iresult = -1
    END SELECT
    
  END SUBROUTINE

! *****************************************************************************
! Routines for the BiCGStab solver
! *****************************************************************************

!<subroutine>
  
  SUBROUTINE sptils_initBiCGStab (rproblem,p_rsolverNode,p_rpreconditioner)
  
!<description>
  ! Creates a t_sptilsNode solver structure for the BiCGStab solver. The node
  ! can be used to directly solve a problem or to be attached as solver
  ! or preconditioner to another solver structure. The node can be deleted
  ! by linsol_releaseSolver.
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
    p_rsolverNode%calgorithm = SPTILS_ALG_BICGSTAB
    
    ! Initialise the ability bitfield with the ability of this solver:
    p_rsolverNode%ccapability = SPTILS_ABIL_CHECKDEF + &
                                SPTILS_ABIL_USESUBSOLVER
    
    ! Allocate the subnode for BiCGStab.
    ! This initialises most of the variables with default values appropriate
    ! to this solver.
    ALLOCATE(p_rsolverNode%p_rsubnodeBiCGStab)
    
    ! Attach the preconditioner if given. 
    
    IF (PRESENT(p_rpreconditioner)) THEN 
      p_rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner => p_rpreconditioner
    END IF

  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE sptils_setMatrixBiCGStab (rsolverNode,Rmatrices)
  
!<description>
  
  ! This routine is called if the system matrix changes.
  ! The routine calls sptils_setMatrices for the preconditioner of BiCGStab
  ! to inform also that one about the change of the matrix pointer.
  
!</description>
  
!<input>
  ! An array of system matrices which is simply passed to the initialisation 
  ! routine of the preconditioner.
  TYPE(t_ccoptSpaceTimeMatrix), DIMENSION(:), INTENT(IN)   :: Rmatrices
!</input>
  
!<inputoutput>
  ! The t_sptilsNode structure of the BiCGStab solver
  TYPE(t_sptilsNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!</subroutine>

    ! Do we have a preconditioner given? If yes, call the general initialisation
    ! routine to initialise it.
    
    IF (ASSOCIATED(rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner)) THEN
      CALL sptils_setMatrices (rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner, &
                              Rmatrices)
    END IF

  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE sptils_initStructureBiCGStab (rsolverNode, ierror)
  
!<description>
  ! Calls the initStructure subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_initStructure.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure of the UMFPACK4 solver
  TYPE(t_sptilsNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to 
  ! SPTILS_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  INTEGER, INTENT(OUT) :: ierror
!</output>

!</subroutine>

    ! local variables
    INTEGER :: isubgroup,i,NEQtime
    TYPE(t_sptilsSubnodeBiCGStab), POINTER :: p_rsubnode
    
    ! A-priori we have no error...
    ierror = SPTILS_ERR_NOERROR

    ! We simply pass isubgroup to the subsolvers when calling them.
    ! Inside of this routine, there's not much to do with isubgroup,
    ! as there is no special information (like a factorisation)
    ! associated with this type of solver, which has to be allocated,
    ! released or updated...

    ! Call the init routine of the preconditioner.
    IF (ASSOCIATED(rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner)) THEN
      CALL sptils_initStructure (rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner, &
                                isubgroup)
    END IF
    
    ! BiCGStab needs 5 temporary vectors + 1 for preconditioning. 
    ! Allocate that here! Use the default data type prescribed in the solver 
    ! structure for allocating the temp vectors.
    p_rsubnode => rsolverNode%p_rsubnodeBiCGStab
    NEQtime = rsolverNode%rmatrix%p_rspaceTimeDiscretisation%NEQtime
    DO i=1,6
      CALL sptivec_initVector (p_rsubnode%RtempVectors(i),NEQtime,&
        rsolverNode%rmatrix%p_rspaceTimeDiscretisation%p_rlevelInfo%p_rdiscretisation)
    END DO
    
    ! Allocate memory for a spatial temp vector.
    CALL lsysbl_createVecBlockByDiscr (&
        rsolverNode%rmatrix%p_rspaceTimeDiscretisation%p_rlevelInfo%p_rdiscretisation,&
        rsolverNode%p_rsubnodeBiCGStab%rtempVectorSpace)
    rsolverNode%p_rsubnodeBiCGStab%rtempVectorSpace%p_rdiscreteBC => &
        rsolverNode%rmatrix%p_rspaceTimeDiscretisation%p_rlevelInfo%p_rdiscreteBC
    rsolverNode%p_rsubnodeBiCGStab%rtempVectorSpace%p_rdiscreteBCfict => &
        rsolverNode%rmatrix%p_rspaceTimeDiscretisation%p_rlevelInfo%p_rdiscreteFBC

  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE sptils_initDataBiCGStab (rsolverNode, ierror)
  
!<description>
  ! Calls the initData subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_initData.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure of the UMFPACK4 solver
  TYPE(t_sptilsNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to 
  ! SPTILS_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  INTEGER, INTENT(OUT) :: ierror
!</output>

!</subroutine>

    ! A-priori we have no error...
    ierror = SPTILS_ERR_NOERROR

    ! We simply pass isubgroup to the subsolvers when calling them.
    ! Inside of this routine, there's not much to do with isubgroup,
    ! as there is no special information (like a factorisation)
    ! associated with this type of solver, which has to be allocated,
    ! released or updated...

    ! Call the init routine of the preconditioner.
    IF (ASSOCIATED(rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner)) THEN
      CALL sptils_initData (rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner, ierror)
    END IF
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE sptils_doneDataBiCGStab (rsolverNode, isolverSubgroup)
  
!<description>
  ! Calls the doneData subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_doneData.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure of the UMFPACK4 solver
  TYPE(t_sptilsNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific 
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  INTEGER, OPTIONAL, INTENT(IN)                    :: isolverSubgroup
!</input>

!</subroutine>

    ! We simply pass isubgroup to the subsolvers when calling them.
    ! Inside of this routine, there's not much to do with isubgroup,
    ! as there is no special information (like a factorisation)
    ! associated with this type of solver, which has to be allocated,
    ! released or updated...

    ! Call the done routine of the preconditioner.
    IF (ASSOCIATED(rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner)) THEN
      CALL sptils_doneData (rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner)
    END IF
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE sptils_doneStructureBiCGStab (rsolverNode, isolverSubgroup)
  
!<description>
  ! Calls the doneStructure subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_doneStructure.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure of the UMFPACK4 solver
  TYPE(t_sptilsNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific 
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  INTEGER, OPTIONAL, INTENT(IN)                    :: isolverSubgroup
!</input>

!</subroutine>

    ! local variables
    INTEGER :: i
    TYPE(t_sptilsSubnodeBiCGStab), POINTER :: p_rsubnode
    
    ! We simply pass isubgroup to the subsolvers when calling them.
    ! Inside of this routine, there's not much to do with isubgroup,
    ! as there is no special information (like a factorisation)
    ! associated with this type of solver, which has to be allocated,
    ! released or updated...

    ! Call the init routine of the preconditioner.
    IF (ASSOCIATED(rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner)) THEN
      CALL sptils_doneStructure (rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner)
    END IF
    
    ! Release temporary data if associated
    p_rsubnode => rsolverNode%p_rsubnodeBiCGStab
    IF (p_rsubnode%RtempVectors(1)%NEQ .NE. 0) THEN
      DO i=6,1,-1
        CALL sptivec_releaseVector(p_rsubnode%RtempVectors(i))
      END DO
    END IF
    
    CALL lsysbl_releaseVector (rsolverNode%p_rsubnodeBiCGStab%rtempVectorSpace)
    
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE sptils_doneBiCGStab (rsolverNode)
  
!<description>
  ! This routine releases all temporary memory for the BiCGStab solver from
  ! the heap. In particular, if a preconditioner is attached to the solver
  ! structure, it's also released from the heap by calling 
  ! linsol_releaseSolver for it.
  ! This DONE routine is declared as RECURSIVE to permit a clean
  ! interaction with linsol_releaseSolver.
!</description>
  
!<input>
  ! A pointer to a t_sptilsNode structure of the BiCGStab solver.
  TYPE(t_sptilsNode), POINTER         :: rsolverNode
!</input>
  
!</subroutine>

    ! Check if there's a preconditioner attached. If yes, release it.
    IF (ASSOCIATED(rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner)) THEN
      CALL sptils_releaseSolver(rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner)
    END IF
    
    ! Release memory if still associated
    CALL sptils_doneDataBiCGStab (rsolverNode)
    CALL sptils_doneStructureBiCGStab (rsolverNode)
    
    ! Release the BiCGStab subnode
    DEALLOCATE(rsolverNode%p_rsubnodeBiCGStab)

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  RECURSIVE SUBROUTINE sptils_precBiCGStab (rsolverNode,rd)
  
!<description>
  ! Applies BiCGStab preconditioner $P \approx A$ to the defect 
  ! vector rd and solves $Pd_{new} = d$.
  ! rd will be overwritten by the preconditioned defect.
  !
  ! The matrix must have been attached to the system before calling
  ! this routine, and the initStructure/initData routines
  ! must have been called to prepare the solver for solving
  ! the problem.
  !
  ! The implementation follows the original paper introducing BiCGStab:
  !   van der Vorst, H.A.; BiCGStab: A Fast and Smoothly Converging
  !   Variant of Bi-CG for the Solution of Nonsymmetric Linear Systems;
  !   SIAM J. Sci. Stat. Comput. 1992, Vol. 13, No. 2, pp. 631-644
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure of the BiCGStab solver
  TYPE(t_sptilsNode), INTENT(INOUT), TARGET :: rsolverNode
   
  ! On call to this routine: The defect vector to be preconditioned.
  ! Will be overwritten by the preconditioned defect.
  TYPE(t_spaceTimeVector), INTENT(INOUT)        :: rd
!</inputoutput>
  
!</subroutine>

  ! local variables
  REAL(DP) :: dalpha,dbeta,domega0,domega1,domega2,dres
  REAL(DP) :: drho1,drho0,dfr
  INTEGER :: ite,i

  ! The system matrix
  TYPE(t_ccoptSpaceTimeDiscretisation), POINTER :: p_rspaceTimeDiscr
  TYPE(t_ccoptSpaceTimeMatrix), POINTER :: p_rmatrix
  
  ! Minimum number of iterations, print-sequence for residuals
  INTEGER :: nminIterations, niteResOutput
  
  ! Whether to filter/prcondition
  LOGICAL bprec,bfilter
  
  ! Our structure
  TYPE(t_sptilsSubnodeBiCGStab), POINTER :: p_rsubnode
  
  ! Pointers to temporary vectors - named for easier access
  TYPE(t_spaceTimeVector), POINTER :: p_DR,p_DR0,p_DP,p_DPA,p_DSA,p_rx
  TYPE(t_sptilsNode), POINTER :: p_rprecSubnode
  
    ! Solve the system!
  
    ! Status reset
    rsolverNode%iresult = 0
    
    ! Getch some information
    p_rsubnode => rsolverNode%p_rsubnodeBiCGStab
    p_rspaceTimeDiscr => rsolverNode%rmatrix%p_rspaceTimeDiscretisation
    p_rmatrix => rsolverNode%rmatrix

    ! Check the parameters
    IF ((rd%NEQtime .EQ. 0) .OR. &
        (p_rspaceTimeDiscr%rtimeDiscr%nintervals .EQ. 0)) THEN
    
      ! Parameters wrong
      rsolverNode%iresult = 2
      RETURN
    END IF

    ! Minimum number of iterations
 
    nminIterations = MAX(rsolverNode%nminIterations,0)
      
    ! Use preconditioning? Filtering?

    bprec = ASSOCIATED(rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner)
    
    ! Iteration when the residuum is printed:

    niteResOutput = MAX(1,rsolverNode%niteResOutput)

    ! Set pointers to the temporary vectors
    p_DR   => p_rsubnode%RtempVectors(1)
    p_DR0  => p_rsubnode%RtempVectors(2)
    p_DP   => p_rsubnode%RtempVectors(3)
    p_DPA  => p_rsubnode%RtempVectors(4)
    p_DSA  => p_rsubnode%RtempVectors(5)
    p_rx   => p_rsubnode%RtempVectors(6)
    
    IF (bprec) THEN
      p_rprecSubnode => p_rsubnode%p_rpreconditioner
    END IF
    
    ! rd is our RHS. p_rx points to a new vector which will be our
    ! iteration vector. At the end of this routine, we replace
    ! rd by p_rx.
    ! Clear our iteration vector p_rx.
    CALL sptivec_clearVector (p_rx)
      
    ! Initialize used vectors with zero
      
    CALL sptivec_clearVector(p_DP)
    CALL sptivec_clearVector(p_DPA)
    
    ! Initialise the iteration vector with zero.

    ! Initialization

    drho0  = 1.0_DP
    dalpha = 1.0_DP
    domega0 = 1.0_DP

    ! Copy our RHS rd to p_DR. As the iteration vector is 0, this
    ! is also our initial defect.

    CALL sptivec_copyVector(rd,p_DR)
    
    ! Filter the defect for boundary conditions in space and time.
    CALL tbc_implementInitCondDefect (&
        p_rspaceTimeDiscr,p_DR,p_rsubnode%rtempVectorSpace)
    CALL tbc_implementBCdefect (rsolverNode%p_rproblem,&
        p_rspaceTimeDiscr,p_DR,p_rsubnode%rtempVectorSpace)

    IF (bprec) THEN
      ! Perform preconditioning with the assigned preconditioning
      ! solver structure.
      CALL sptils_precondDefect (p_rprecSubnode,p_DR)
    END IF
    
    ! Get the norm of the residuum
    dres = sptivec_vectorNorm (p_DR,rsolverNode%iresNorm)
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
        CALL output_line ('Space-Time-BiCGStab: Iteration '// &
             TRIM(sys_siL(0,10))//',  !!RES!! = '//&
             TRIM(sys_sdEL(rsolverNode%dinitialDefect,15)) )
      END IF

      CALL sptivec_copyVector(p_DR,p_DR0)

      ! Perform at most nmaxIterations loops to get a new vector

      DO ite = 1,rsolverNode%nmaxIterations
      
        rsolverNode%icurrentIteration = ite

        drho1 = sptivec_scalarProduct (p_DR0,p_DR) 

        IF (drho0*domega0 .EQ. 0.0_DP) THEN
          ! Should not happen
          IF (rsolverNode%ioutputLevel .GE. 2) THEN
            CALL output_line ('Space-Time-BiCGStab: Iteration prematurely stopped! '//&
                 'Correction vector is zero!')
          END IF

          ! Some tuning for the output, then cancel.

          rsolverNode%iresult = -1
          rsolverNode%iiterations = ITE-1
          EXIT
          
        END IF

        dbeta=(drho1*dalpha)/(drho0*domega0)
        drho0 = drho1

        CALL sptivec_vectorLinearComb (p_DR ,p_DP,1.0_DP,dbeta)
        CALL sptivec_vectorLinearComb (p_DPA ,p_DP,-dbeta*domega0,1.0_DP)

        CALL c2d2_spaceTimeMatVec (rsolverNode%p_rproblem, p_rmatrix, &
            p_DP,p_DPA, 1.0_DP,0.0_DP,SPTID_FILTER_NONE)
        
        ! Filter the defect for boundary conditions in space and time.
        CALL tbc_implementInitCondDefect (&
            p_rspaceTimeDiscr,p_DPA,p_rsubnode%rtempVectorSpace)
        CALL tbc_implementBCdefect (rsolverNode%p_rproblem,&
            p_rspaceTimeDiscr,p_DPA,p_rsubnode%rtempVectorSpace)

        IF (bprec) THEN
          ! Perform preconditioning with the assigned preconditioning
          ! solver structure.
          CALL sptils_precondDefect (p_rprecSubnode,p_DPA)
        END IF

        dalpha = sptivec_scalarProduct (p_DR0,p_DPA)
        
        IF (dalpha .EQ. 0.0_DP) THEN
          ! We are below machine exactness - we can't do anything more...
          ! May happen with very small problems with very few unknowns!
          IF (rsolverNode%ioutputLevel .GE. 2) THEN
            CALL output_line ('Space-Time-BiCGStab: Convergence failed, ALPHA=0!')
            rsolverNode%iresult = -2
            EXIT
          END IF
        END IF
        
        dalpha = drho1/dalpha

        CALL sptivec_vectorLinearComb (p_DPA,p_DR,-dalpha,1.0_DP)

        CALL c2d2_spaceTimeMatVec (rsolverNode%p_rproblem, p_rmatrix, &
            p_DR,p_DSA, 1.0_DP,0.0_DP,SPTID_FILTER_NONE)
        
        ! Filter the defect for boundary conditions in space and time.
        CALL tbc_implementInitCondDefect (&
            p_rspaceTimeDiscr,p_DSA,p_rsubnode%rtempVectorSpace)
        CALL tbc_implementBCdefect (rsolverNode%p_rproblem,&
            p_rspaceTimeDiscr,p_DSA,p_rsubnode%rtempVectorSpace)

        IF (bprec) THEN
          ! Perform preconditioning with the assigned preconditioning
          ! solver structure.
          CALL sptils_precondDefect (p_rprecSubnode,p_DSA)
        END IF
        
        domega1 = sptivec_scalarProduct (p_DSA,p_DR)
        domega2 = sptivec_scalarProduct (p_DSA,p_DSA)
        
        IF (domega1 .EQ. 0.0_DP) THEN
          domega0 = 0.0_DP
        ELSE
          IF (domega2 .EQ. 0.0_DP) THEN
            IF (rsolverNode%ioutputLevel .GE. 2) THEN
              CALL output_line ('Space-Time-BiCGStab: Convergence failed: omega=0!')
              rsolverNode%iresult = -2
              EXIT
            END IF
          END IF
          domega0 = domega1/domega2
        END IF

        CALL sptivec_vectorLinearComb (p_DP ,p_rx,dalpha,1.0_DP)
        CALL sptivec_vectorLinearComb (p_DR ,p_rx,domega0,1.0_DP)

        CALL sptivec_vectorLinearComb (p_DSA,p_DR,-domega0,1.0_DP)

        ! Get the norm of the new (final?) residuum
        dfr = sptivec_vectorNorm (p_DR,rsolverNode%iresNorm)
     
        rsolverNode%dfinalDefect = dfr

        ! Test if the iteration is diverged
        IF (sptils_testDivergence(rsolverNode,dfr)) THEN
          CALL output_line ('Space-Time-BiCGStab: Solution diverging!')
          rsolverNode%iresult = 1
          EXIT
        END IF
     
        ! At least perform nminIterations iterations
        IF (ite .GE. nminIterations) THEN
        
          ! Check if the iteration converged
          IF (sptils_testConvergence(rsolverNode,dfr)) EXIT
          
        END IF

        ! print out the current residuum

        IF ((rsolverNode%ioutputLevel .GE. 2) .AND. &
            (MOD(ite,niteResOutput).EQ.0)) THEN
          CALL output_line ('Space-Time-BiCGStab: Iteration '// &
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
        CALL output_line ('Space-Time-BiCGStab: Iteration '// &
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
        CALL output_line ('Space-Time-BiCGStab statistics:')
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
              'Space-Time-BiCGStab: Iterations/Rate of convergence: '//&
              TRIM(sys_siL(rsolverNode%iiterations,10))//' /'//&
              TRIM(sys_sdEL(rsolverNode%dconvergenceRate,15)) )
      END IF
      
    ELSE
      ! DEF=Infinity; RHO=Infinity, set to 1
      rsolverNode%dconvergenceRate = 1.0_DP
    END IF  
  
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
    
    ! Initialise the ability bitfield with the ability of this solver:
    p_rsolverNode%ccapability = SPTILS_ABIL_MULTILEVEL + SPTILS_ABIL_CHECKDEF     + &
                                SPTILS_ABIL_USESUBSOLVER

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
  ! The t_sptilsNode structure which is to be cleaned up.
  TYPE(t_sptilsNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!</subroutine>
    INTEGER :: ilev
  
    ! Release memory if still associated
    CALL sptils_doneDataMultigrid (rsolverNode)
    CALL sptils_doneStructureMultigrid (rsolverNode)
    
    ! Release level information
    DO ilev=rsolverNode%p_rsubnodeMultigrid%NLMAX,rsolverNode%p_rsubnodeMultigrid%NLMIN,-1
    
      ! Pre- and postsmoother may be identical; release them only once!
      IF (ASSOCIATED(rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%p_rpresmoother)) THEN
        IF (.NOT. ASSOCIATED(&
            rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%p_rpresmoother,&
            rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%p_rpostsmoother)) THEN
          CALL sptils_releaseSolver(&
              rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%p_rpresmoother)
        ELSE
          NULLIFY(rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%p_rpresmoother)
        END IF
      END IF

      IF (ASSOCIATED(&
          rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%p_rpostsmoother)) THEN
        CALL sptils_releaseSolver(&
            rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%p_rpostsmoother)
      END IF
      
      IF (ASSOCIATED(&
          rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%p_rcoarseGridSolver)) THEN
        CALL sptils_releaseSolver(&
            rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%p_rcoarseGridSolver)
      END IF
      
    END DO
    
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
  TYPE(t_ccoptSpaceTimeMatrix), DIMENSION(:), INTENT(IN) :: Rmatrices
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
      rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%rmatrix = &
          Rmatrices(ilev)
      
      ! Call the setmatrices-routine for the presmoother/postsmoother/
      ! coarse grid solver
      IF (.NOT. ASSOCIATED(&
          rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%p_rpresmoother,&
          rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%p_rpostsmoother)) THEN
        ! May be the presmoother and postsmoother are identical; release them
        ! only once!
        IF (ASSOCIATED(rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%p_rpresmoother)) THEN
          CALL sptils_setMatrices(&
              rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%p_rpresmoother,&
              Rmatrices(1:ilev))
        END IF
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
  ! The t_sptilsNode structure 
  TYPE(t_sptilsNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the SPTILS_ERR_XXXX constants. A value different to 
  ! SPTILS_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  INTEGER, INTENT(OUT) :: ierror
!</output>
  
!</subroutine>

    ! local variables
    INTEGER :: ilev,NLMAX,NEQtime
    INTEGER(PREC_VECIDX) :: NEQ
    TYPE(t_sptilsMGLevelInfo), POINTER :: p_rmgLevel
    
    ! A-priori we have no error...
    ierror = SPTILS_ERR_NOERROR
    
    ! Maximum time level.
    NLMAX = SIZE(rsolverNode%p_rsubnodeMultigrid%p_Rlevels)

    ! On each level, call the initStructure routine of the
    ! presmoother/postsmoother/coarse grid solver

    ! Loop through the level. 
    DO ilev=1,NLMAX
    
      p_rmgLevel => rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)
    
      ! Call the setmatrices-routine for the presmoother/postsmoother/
      ! coarse grid solver
      IF (.NOT. ASSOCIATED(p_rmgLevel%p_rpresmoother,p_rmgLevel%p_rpostsmoother)) THEN
        ! May be the presmoother and postsmoother are identical; initialise them
        ! only once!
        IF (ASSOCIATED(p_rmgLevel%p_rpresmoother)) THEN
          CALL sptils_initStructure(p_rmgLevel%p_rpresmoother,ierror)
          IF (ierror .NE. ierror) RETURN
        END IF
      END IF
      
      IF (ASSOCIATED(p_rmgLevel%p_rpostsmoother)) THEN
        CALL sptils_initStructure(p_rmgLevel%p_rpostsmoother,ierror)
        IF (ierror .NE. ierror) RETURN
      END IF
      
      IF (ASSOCIATED(p_rmgLevel%p_rcoarseGridSolver)) THEN
        CALL sptils_initStructure(p_rmgLevel%p_rcoarseGridSolver,ierror)
        IF (ierror .NE. ierror) RETURN
      END IF
      
      ! Generate an interlevel projection structure for that level
      !CALL sptipr_initProjection (&
      !    rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%rinterlevelProjection,&
      !    rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%&
      !        rspaceTimeDiscr%p_rlevelInfo%p_rdiscretisation)
              
      NEQtime = p_rmgLevel%rmatrix%p_rspaceTimeDiscretisation%NEQtime
              
      ! On all levels except for the maximum one, create a solution vector
      IF (ilev .LT. NLMAX) THEN
        CALL sptivec_initVector (&
            p_rmgLevel%rsolutionVector,NEQtime,p_rmgLevel%&
            rmatrix%p_rspaceTimeDiscretisation%p_rlevelInfo%p_rdiscretisation)
      END IF
      
      ! On all levels except for the first one, create a RHS and a temp vector
      IF (ilev .GT. 1) THEN
        CALL sptivec_initVector (&
            p_rmgLevel%rrhsVector,NEQtime,p_rmgLevel%&
            rmatrix%p_rspaceTimeDiscretisation%p_rlevelInfo%p_rdiscretisation)
        CALL sptivec_initVector (&
            p_rmgLevel%rtempVector,NEQtime,p_rmgLevel%&
            rmatrix%p_rspaceTimeDiscretisation%p_rlevelInfo%p_rdiscretisation)
      END IF
      
      ! If adaptive coarse grid correction is activated, we need a temporary
      ! vector for the adaptive coarse grid correction on every level.
      IF (rsolverNode%p_rsubnodeMultigrid%dalphamin .NE. &
          rsolverNode%p_rsubnodeMultigrid%dalphamax) THEN
        CALL sptivec_initVector (&
            p_rmgLevel%rtempCGCvector,NEQtime,p_rmgLevel%&
            rmatrix%p_rspaceTimeDiscretisation%p_rlevelInfo%p_rdiscretisation)
      END IF
      
      ! Create a block temp vector for the interlevel projection
      CALL lsysbl_createVecBlockByDiscr(&
          p_rmgLevel%rmatrix%p_rspaceTimeDiscretisation%p_rlevelInfo%p_rdiscretisation,&
          p_rmgLevel%rprjVector,.FALSE.)

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
  ! The t_sptilsNode structure 
  TYPE(t_sptilsNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the SPTILS_ERR_XXXX constants. A value different to 
  ! SPTILS_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  INTEGER, INTENT(OUT) :: ierror
!</output>

!</subroutine>

    ! local variables
    INTEGER :: ilev
    TYPE(t_sptilsMGLevelInfo), POINTER :: p_rmgLevel
    
    ! A-priori we have no error...
    ierror = SPTILS_ERR_NOERROR
    
    ! On each level, call the initStructure routine of the
    ! presmoother/postsmoother/coarse grid solver

    ! Loop through the level. 
    DO ilev=1,SIZE(rsolverNode%p_rsubnodeMultigrid%p_Rlevels)
     
      p_rmgLevel => rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)
    
      ! Call the setmatrices-routine for the presmoother/postsmoother/
      ! coarse grid solver
      IF (.NOT. ASSOCIATED(p_rmgLevel%p_rpresmoother,p_rmgLevel%p_rpostsmoother)) THEN
        ! May be the presmoother and postsmoother are identical; initialise them
        ! only once!
        IF (ASSOCIATED(p_rmgLevel%p_rpresmoother)) THEN
          CALL sptils_initData(p_rmgLevel%p_rpresmoother,ierror)
          IF (ierror .NE. ierror) RETURN
        END IF
      END IF
      
      IF (ASSOCIATED(p_rmgLevel%p_rpostsmoother)) THEN
        CALL sptils_initData(p_rmgLevel%p_rpostsmoother,ierror)
        IF (ierror .NE. ierror) RETURN
      END IF
      
      IF (ASSOCIATED(p_rmgLevel%p_rcoarseGridSolver)) THEN
        CALL sptils_initData(p_rmgLevel%p_rcoarseGridSolver,ierror)
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
  ! The t_sptilsNode structure 
  TYPE(t_sptilsNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!</subroutine>

    ! local variables
    INTEGER :: ilev
    TYPE(t_sptilsMGLevelInfo), POINTER :: p_rmgLevel
    
    ! On each level, call the initStructure routine of the
    ! presmoother/postsmoother/coarse grid solver

    ! Loop through the level. 
    DO ilev=1,SIZE(rsolverNode%p_rsubnodeMultigrid%p_Rlevels)
    
      p_rmgLevel => rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)
    
      ! Call the setmatrices-routine for the presmoother/postsmoother/
      ! coarse grid solver
      IF (.NOT. ASSOCIATED(p_rmgLevel%p_rpresmoother,p_rmgLevel%p_rpostsmoother)) THEN
        ! May be the presmoother and postsmoother are identical; release them
        ! only once!
        IF (ASSOCIATED(p_rmgLevel%p_rpresmoother)) THEN
          CALL sptils_doneData(p_rmgLevel%p_rpresmoother)
        END IF
      END IF
      
      IF (ASSOCIATED(p_rmgLevel%p_rpostsmoother)) THEN
        CALL sptils_doneData(p_rmgLevel%p_rpostsmoother)
      END IF
      
      IF (ASSOCIATED(p_rmgLevel%p_rcoarseGridSolver)) THEN
        CALL sptils_doneData(p_rmgLevel%p_rcoarseGridSolver)
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
  ! The t_sptilsNode structure 
  TYPE(t_sptilsNode), INTENT(INOUT)         :: rsolverNode
!</inputoutput>
  
!</subroutine>

    ! local variables
    INTEGER :: ilev,NLMAX
    TYPE(t_sptilsMGLevelInfo), POINTER :: p_rmgLevel
    
    ! On each level, call the initStructure routine of the
    ! presmoother/postsmoother/coarse grid solver

    ! Loop through the level. 
    NLMAX = SIZE(rsolverNode%p_rsubnodeMultigrid%p_Rlevels)
    DO ilev=1,NLMAX
    
      p_rmgLevel => rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)
    
      ! Call the setmatrices-routine for the presmoother/postsmoother/
      ! coarse grid solver
      IF (.NOT. ASSOCIATED(p_rmgLevel%p_rpresmoother,p_rmgLevel%p_rpostsmoother)) THEN
        ! May be the presmoother and postsmoother are identical; release them
        ! only once!
        IF (ASSOCIATED(p_rmgLevel%p_rpresmoother)) THEN
          CALL sptils_doneStructure(p_rmgLevel%p_rpresmoother)
        END IF
      END IF
      
      IF (ASSOCIATED(p_rmgLevel%p_rpostsmoother)) THEN
        CALL sptils_doneStructure(p_rmgLevel%p_rpostsmoother)
      END IF
      
      IF (ASSOCIATED(p_rmgLevel%p_rcoarseGridSolver)) THEN
        CALL sptils_doneStructure(p_rmgLevel%p_rcoarseGridSolver)
      END IF


      ! Release the projection structure
      CALL sptipr_doneProjection (p_rmgLevel%rinterlevelProjection)
              
      ! Release vectors
              
      IF (ilev .LT. NLMAX) THEN
        CALL sptivec_releaseVector (p_rmgLevel%rsolutionVector)
      END IF
      
      IF (ilev .GT. 1) THEN
        CALL sptivec_releaseVector (p_rmgLevel%rrhsVector)
        CALL sptivec_releaseVector (p_rmgLevel%rtempVector)
      END IF
      
      ! If adaptive coarse grid correction is activated, we need a temporary
      ! vector for the adaptive coarse grid correction on every level.
      IF (rsolverNode%p_rsubnodeMultigrid%dalphamin .NE. &
          rsolverNode%p_rsubnodeMultigrid%dalphamax) THEN
        CALL sptivec_releaseVector (p_rmgLevel%rtempCGCvector)
      END IF

      ! Release the temp vector for prolongation/restriction
      CALL lsysbl_releaseVector (p_rmgLevel%rprjVector)

    END DO

  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE sptils_setMultigridLevel (rsolverNode,ilevel,&
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
    IF ((rsolverNode%calgorithm .NE. SPTILS_ALG_MULTIGRID) .OR. &
        (.NOT. ASSOCIATED(rsolverNode%p_rsubnodeMultigrid))) THEN
      PRINT *,'Error: Multigrid structure not initialised'
      CALL sys_halt()
    END IF
    
    p_rlevelInfo => rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilevel)
    
    ! Attach the sub-solvers
    IF (PRESENT(p_rcoarseGridSolver)) THEN
      p_rlevelInfo%p_rcoarseGridSolver => p_rcoarseGridSolver
    ELSE
      NULLIFY(p_rlevelInfo%p_rcoarseGridSolver)
    END IF

    IF (PRESENT(p_rpresmoother)) THEN
      p_rlevelInfo%p_rpresmoother => p_rpresmoother
    ELSE
      NULLIFY(p_rlevelInfo%p_rpresmoother)
    END IF

    IF (PRESENT(p_rpostsmoother)) THEN
      p_rlevelInfo%p_rpostsmoother => p_rpostsmoother
    ELSE
      NULLIFY(p_rlevelInfo%p_rpostsmoother)
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
  
  RECURSIVE SUBROUTINE sptils_smoothCorrection (rsolverNode,rx,rb,rtemp,rspatialTemp)
  
!<description>
  ! This routine performs a smoothing process on the vector rx
  ! belonging to a linear system $Ax=b$.
  ! rsolverNode identifies a solver structure that is converted to a
  ! smoother using sptils_convertToSmoother: The number of smoothing
  ! steps is written to rsolverNode%nmaxIterations and 
  ! rsolverNode%nmaxIterations so that the solver performs a definite
  ! number of iterations regardless of the residual.
  !
  ! rx is assumed to be of 'defect' type. There are two types of smoothing
  ! processes, depending on which type of preconditioner rsolverNode is:
  ! 1.) If rsolverNode is an iterative solver, sptils_smoothCorrection
  !     calls the associated solver P to compute $x=P^{-1}b$ using
  !     rsolverNode%nmaxIterations iterations.
  ! 2.) If rsolverNode is a 1-step solver, sptils_smoothCorrection
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

  ! A spatial temporary vector of the same size and structure as every timestep in rx.
  TYPE(t_vectorBlock), INTENT(INOUT)                     :: rspatialTemp
!</inputoutput>
  
!</subroutine>

    INTEGER :: i
    INTEGER :: iiterations
    REAL(DP) :: dres,dresInit
    TYPE(t_ccoptSpaceTimeMatrix), POINTER :: p_rmatrix
    TYPE(t_ccoptSpaceTimeDiscretisation), POINTER :: p_rspaceTimeDiscr
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
    p_rmatrix => rsolverNode%rmatrix
    p_rspaceTimeDiscr => rsolverNode%rmatrix%p_rspaceTimeDiscretisation
    
    ! Do we have an iterative or one-step solver given?
    ! A 1-step solver performs the following loop nmaxIterations times, while an iterative
    ! solver is called only once and performs nmaxIterations steps internally.
    ! (This is a convention. Calling an iterative solver i times with j internal steps
    ! would also be possible, but we don't implement that here.)
    IF (IAND(rsolverNode%ccapability,SPTILS_ABIL_DIRECT) .NE. 0) THEN
      iiterations = rsolverNode%nmaxIterations
    ELSE
      iiterations = 1
    END IF
    
    DO i=1,iiterations
    
      CALL sptivec_copyVector(rb,rtemp)
      CALL c2d2_spaceTimeMatVec (rsolverNode%p_rproblem, p_rmatrix, &
          rx,rtemp, -1.0_DP,1.0_DP,SPTID_FILTER_DEFECT,dres)
      
      IF (iiterations .EQ. 1) dresInit = dres
      IF (dresInit .EQ. 0) dresInit = 1.0_DP
      
      ! Implement boundary conditions into the defect
      CALL tbc_implementInitCondDefect (&
          p_rspaceTimeDiscr,rtemp,rspatialTemp)
      CALL tbc_implementBCdefect (rsolverNode%p_rproblem,&
          p_rspaceTimeDiscr,rtemp,rspatialTemp)
      
      IF (rsolverNode%ioutputLevel .GE. 2) THEN
        IF (.NOT.((dres .GE. 1E-99_DP) .AND. (dres .LE. 1E99_DP))) dres = 0.0_DP
                  
        CALL output_line ('Space-Time-Smoother: Step '//TRIM(sys_siL(i-1,10))//&
            ' !!RES!! = '//TRIM(sys_sdEL(dres,15)) )
      END IF
      
      ! Check for convergence
      IF (rsolverNode%istoppingCriterion .EQ. 0) THEN
        IF ((dres .LT. rsolverNode%depsAbs) .AND. &
            (dres .LT. rsolverNode%depsRel*dresInit)) EXIT
      ELSE
        IF ((dres .LT. rsolverNode%depsAbs) .OR. &
            (dres .LT. rsolverNode%depsRel*dresInit)) EXIT
      END IF
      
      CALL sptils_precondDefect(rsolverNode,rtemp)
      CALL sptivec_vectorLinearComb (rtemp,rx,1.0_DP,1.0_DP)
      
    END DO

    ! Probably print the final residuum
    IF (rsolverNode%ioutputLevel .GE. 2) THEN
      IF (i .GE. iiterations) THEN
        ! We only have to recalculate the residuum if we haven't stopped
        ! prematurely.
        CALL sptivec_copyVector(rb,rtemp)
        CALL c2d2_spaceTimeMatVec (rsolverNode%p_rproblem, p_rmatrix, &
            rx,rtemp, -1.0_DP,1.0_DP,SPTID_FILTER_DEFECT,dres)
        IF (.NOT.((dres .GE. 1E-99_DP) .AND. (dres .LE. 1E99_DP))) dres = 0.0_DP
                  
        CALL output_line ('Space-Time-Smoother: Step '//TRIM(sys_siL(i-1,10))//&
            ' !!RES!! = '//TRIM(sys_sdEL(dres,15)) )
      END IF
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
  TYPE(t_ccoptSpaceTimeMatrix), POINTER :: p_rmatrix
  TYPE(t_ccoptSpaceTimeDiscretisation), POINTER :: p_rspaceTimeDiscr
  
  ! Our MG structure
  TYPE(t_sptilsSubnodeMultigrid), POINTER :: p_rsubnode
  
    ! Solve the system!
    
    ! Getch some information
    p_rsubnode => rsolverNode%p_rsubnodeMultigrid
    
    ! Get the system matrix on the finest level:
    p_rmatrix => rsolverNode%rmatrix
    p_rspaceTimeDiscr => rsolverNode%rmatrix%p_rspaceTimeDiscretisation

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
    
    ! We start on the maximum level. 
    ilev = p_rsubnode%NLMAX
    nlmax = p_rsubnode%NLMAX
    nlmin = p_rsubnode%NLMIN
    
    ! Is there only one level? Can be seen if the current level
    ! already contains a coarse grid solver.
    IF (nlmin .EQ. nlmax) THEN
    
      IF (rsolverNode%ioutputLevel .GT. 1) THEN
        CALL output_line ('Space-Time-Multigrid: Only one level. '//&
             'Switching back to standard solver.')
      END IF

      ! DEBUG!!!
      !CALL sptivec_saveToFileSequence (rd,&
      !    '(''./debugdata/vector2.txt.'',I5.5)',.TRUE.)

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
          CALL output_line ('Space-Time-Multigrid: Iteration '// &
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
          p_rmatrix => p_rsubnode%p_Rlevels(ilev)%rmatrix
          
          IF (rsolverNode%ioutputLevel .GE. 3) THEN
            CALL output_line (&
              'Space-Time-Multigrid: Current mesh level: '//TRIM(sys_siL(ilev,5)))
          END IF
          
          ! Build the defect...
          CALL sptivec_copyVector (&
              p_rsubnode%p_Rlevels(ilev)%rrhsVector,&
              p_rsubnode%p_Rlevels(ilev)%rtempVector)
          IF (ite .NE. 1) THEN   ! initial solution vector is zero!
            CALL c2d2_spaceTimeMatVec (rsolverNode%p_rproblem, p_rmatrix, &
                p_rsubnode%p_Rlevels(ilev)%rsolutionVector,&
                p_rsubnode%p_Rlevels(ilev)%rtempVector, -1.0_DP,1.0_DP,&
                SPTID_FILTER_DEFECT)
          
            ! Implement boundary conditions into the defect
            CALL tbc_implementInitCondDefect (&
                p_rsubnode%p_Rlevels(ilev)%rmatrix%p_rspaceTimeDiscretisation,&
                p_rsubnode%p_Rlevels(ilev)%rtempVector,&
                p_rsubnode%p_Rlevels(ilev)%rprjVector)
            CALL tbc_implementBCdefect (rsolverNode%p_rproblem,&
                p_rsubnode%p_Rlevels(ilev)%rmatrix%p_rspaceTimeDiscretisation,&
                p_rsubnode%p_Rlevels(ilev)%rtempVector,&
                p_rsubnode%p_Rlevels(ilev)%rprjVector)
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
                          p_rsubnode%p_Rlevels(ilev)%rtempVector,&
                          p_rsubnode%p_Rlevels(ilev)%rprjVector)
              END IF
            
              ! Build the defect vector
              CALL sptivec_copyVector (p_rsubnode%p_Rlevels(ilev)%rrhsVector,&
                  p_rsubnode%p_Rlevels(ilev)%rtempVector)
              CALL c2d2_spaceTimeMatVec (rsolverNode%p_rproblem, p_rmatrix, &
                  p_rsubnode%p_Rlevels(ilev)%rsolutionVector,&
                  p_rsubnode%p_Rlevels(ilev)%rtempVector, -1.0_DP,1.0_DP,&
                  SPTID_FILTER_DEFECT,dres)
              
              ! Extended output
              IF (ASSOCIATED(p_rsubnode%p_Rlevels(ilev)%p_rpreSmoother) .AND. &
                  (rsolverNode%ioutputLevel .GE. 3) .AND. &
                  (MOD(ite,niteResOutput) .EQ. 0)) THEN
                  
                IF (.NOT.((dres .GE. 1E-99_DP) .AND. &
                          (dres .LE. 1E99_DP))) dres = 0.0_DP
                          
                CALL output_line ('Space-Time-Multigrid: Level '//TRIM(sys_siL(ilev,5))//&
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

                ! Implement boundary conditions into the defect
                CALL tbc_implementInitCondDefect (&
                    p_rsubnode%p_Rlevels(ilev-1)%rmatrix%p_rspaceTimeDiscretisation,&
                    p_rsubnode%p_Rlevels(ilev-1)%rrhsVector,&
                    p_rsubnode%p_Rlevels(ilev-1)%rprjVector)
                CALL tbc_implementBCdefect (rsolverNode%p_rproblem,&
                    p_rsubnode%p_Rlevels(ilev-1)%rmatrix%p_rspaceTimeDiscretisation,&
                    p_rsubnode%p_Rlevels(ilev-1)%rrhsVector,&
                    p_rsubnode%p_Rlevels(ilev-1)%rprjVector)

                ! Choose zero as initial vector on lower level. 
                CALL sptivec_clearVector (p_rsubnode%p_Rlevels(ilev-1)%rsolutionVector)
                
                ! Extended output and/or adaptive cycles
                IF (rsolverNode%ioutputLevel .GE. 3) THEN
                
                  dres = sptivec_vectorNorm (p_rsubnode%p_Rlevels(ilev-1)%rrhsVector,&
                      rsolverNode%iresNorm)
                  IF (.NOT.((dres .GE. 1E-99_DP) .AND. &
                            (dres .LE. 1E99_DP))) dres = 0.0_DP
                            
                  ! In case adaptive cycles are activated, save the 'initial' residual
                  ! of that level into the level structure. Then we can later check
                  ! if we have to repeat the cycle on the coarse mesh.
                  IF ((p_rsubnode%p_Rlevels(ilev-1)%depsRelCycle .NE. 1E99_DP) .OR.&
                      (p_rsubnode%p_Rlevels(ilev-1)%depsAbsCycle .NE. 1E99_DP)) THEN
                    p_rsubnode%p_Rlevels(ilev-1)%dinitResCycle = dres
                    p_rsubnode%p_Rlevels(ilev-1)%icycleCount = 1
                  END IF

                  ! If the output level is high enough, print that residuum norm.   
                  IF (MOD(ite,niteResOutput).EQ. 0) THEN
                    CALL output_line ('Space-Time-Multigrid: Level '//TRIM(sys_siL(ilev-1,5))//&
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

                ! DEBUG!!!
                !CALL sptivec_saveToFileSequence (p_rsubnode%p_Rlevels(ilev-1)%rsolutionVector,&
                !    '(''./debugdata/vector1.txt.'',I5.5)',.TRUE.)

                ! Implement boundary conditions into the defect
                CALL tbc_implementInitCondDefect (&
                    p_rsubnode%p_Rlevels(ilev-1)%rmatrix%p_rspaceTimeDiscretisation,&
                    p_rsubnode%p_Rlevels(ilev-1)%rsolutionVector,&
                    p_rsubnode%p_Rlevels(ilev-1)%rprjVector)
                CALL tbc_implementBCdefect (rsolverNode%p_rproblem,&
                    p_rsubnode%p_Rlevels(ilev-1)%rmatrix%p_rspaceTimeDiscretisation,&
                    p_rsubnode%p_Rlevels(ilev-1)%rsolutionVector,&
                    p_rsubnode%p_Rlevels(ilev-1)%rprjVector)

                ! Extended output
                IF ((rsolverNode%ioutputLevel .GE. 3) .AND. &
                    (MOD(ite,niteResOutput).EQ.0)) THEN
                    
                  dres = sptivec_vectorNorm (p_rsubnode%p_Rlevels(ilev-1)%rsolutionVector,&
                      rsolverNode%iresNorm)
                  IF (.NOT.((dres .GE. 1E-99_DP) .AND. &
                            (dres .LE. 1E99_DP))) dres = 0.0_DP
                            
                  CALL output_line ('Space-Time-Multigrid: Level '//TRIM(sys_siL(ilev-1,5))//&
                      ' after restrict.:  !!RES!! = '//TRIM(sys_sdEL(dres,15)) )
                END IF

              END IF              
            
              ! Go down one level
              ilev = ilev - 1
              p_rmatrix => p_rsubnode%p_Rlevels(ilev)%rmatrix

              IF (rsolverNode%ioutputLevel .GE. 3) THEN
                CALL output_line ('Space-Time-Multigrid: Current mesh level: '//TRIM(sys_siL(ilev,5)))
              END IF
              
              ! If we are not on the lowest level, repeat the smoothing of 
              ! the solution/restriction of the new defect in the next loop 
              ! pass...
            END DO   ! ilev > minimum level
            
            ! Now we reached the coarse grid.
            
            IF (rsolverNode%ioutputLevel .GE. 3) THEN
              CALL output_line ('Space-Time-Multigrid: Invoking coarse grid solver.')
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
              p_rmatrix => p_rsubnode%p_Rlevels(ilev)%rmatrix
              
              IF (rsolverNode%ioutputLevel .GE. 3) THEN
                CALL output_line ('Space-Time-Multigrid: Current mesh level: '&
                    //TRIM(sys_siL(ilev,5)))
              END IF

              ! Prolongate the solution vector from the coarser level
              ! to the temp vector on the finer level.
              CALL sptipr_performProlongation (&
                    p_rsubnode%p_Rlevels(ilev)%rinterlevelProjection,&
                    p_rsubnode%p_Rlevels(ilev-1)%rsolutionVector, &
                    p_rsubnode%p_Rlevels(ilev)%rtempVector, &
                    p_rsubnode%p_Rlevels(ilev-1)%rprjVector, &
                    p_rsubnode%p_Rlevels(ilev)%rprjVector,&
                    p_rmatrix%p_rspaceTimeDiscretisation,&
                    p_rsubnode%p_Rlevels(ilev-1)%rmatrix%p_rspaceTimeDiscretisation,&
                    rsolverNode%p_rproblem)
                    
              ! Implement boundary conditions into the vector.
              ! It's still a defect, although a preconditioned one.
              CALL tbc_implementInitCondDefect (&
                  p_rsubnode%p_Rlevels(ilev)%rmatrix%p_rspaceTimeDiscretisation,&
                  p_rsubnode%p_Rlevels(ilev)%rtempVector, &
                  p_rsubnode%p_Rlevels(ilev)%rprjVector)
              CALL tbc_implementBCdefect (rsolverNode%p_rproblem,&
                  p_rsubnode%p_Rlevels(ilev)%rmatrix%p_rspaceTimeDiscretisation,&
                  p_rsubnode%p_Rlevels(ilev)%rtempVector,&
                  p_rsubnode%p_Rlevels(ilev)%rprjVector)

              ! Step length control. By default, choose step length dalphamin
              ! which is usually = 1.0.
              dstep = p_rsubnode%dalphamin ! 1.0_DP
              
              ! If adaptive coarse grid correction is activated, get the
              ! step length parameter by energy minimisation.
              IF (p_rsubnode%dalphamin .NE. p_rsubnode%dalphamax) THEN
              
                ! Calculate the optimal alpha.
                CALL calcCGCenergyMin (rsolverNode%p_rproblem,&
                    p_rsubnode%p_Rlevels(ilev)%rsolutionVector,&
                    p_rsubnode%p_Rlevels(ilev)%rtempVector,&
                    p_rsubnode%p_Rlevels(ilev)%rrhsVector,&
                    p_rmatrix,&
                    p_rsubnode%p_Rlevels(ilev)%rtempCGCvector,&
                    p_rsubnode%p_Rlevels(ilev)%rprjVector,&
                    p_rsubnode%dalphamin, p_rsubnode%dalphamax, dstep)
                
                ! Some output
                IF (rsolverNode%ioutputLevel .GE. 2) THEN
                  CALL output_line ('Coarse grid correction factor: '//&
                    'dstep = '//TRIM(sys_sdEL(dstep,15)) )
                END IF
                
              END IF
              
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
                    p_rsubnode%p_Rlevels(ilev)%rtempVector, -1.0_DP,1.0_DP,&
                    SPTID_FILTER_DEFECT,dres)

                IF (.NOT.((dres .GE. 1E-99_DP) .AND. &
                          (dres .LE. 1E99_DP))) dres = 0.0_DP
                          
                CALL output_line ('Space-Time-Multigrid: Level '//TRIM(sys_siL(ilev,5))//&
                    ' after c.g.corr.:  !!RES!! = '//TRIM(sys_sdEL(dres,15)) )
              END IF
                                            
              ! Perform the post-smoothing with the current solution vector
              IF (ASSOCIATED(p_rsubnode%p_Rlevels(ilev)%p_rpostSmoother)) THEN
                CALL sptils_smoothCorrection (&
                          p_rsubnode%p_Rlevels(ilev)%p_rpostSmoother,&
                          p_rsubnode%p_Rlevels(ilev)%rsolutionVector,&
                          p_rsubnode%p_Rlevels(ilev)%rrhsVector,&
                          p_rsubnode%p_Rlevels(ilev)%rtempVector,&
                          p_rsubnode%p_Rlevels(ilev)%rprjVector)

                ! Extended output
                IF ((rsolverNode%ioutputLevel .GE. 3) .AND. &
                    (MOD(ite,niteResOutput).EQ.0)) THEN
                    
                  CALL sptivec_copyVector (p_rsubnode%p_Rlevels(ilev)%rrhsVector,&
                                          p_rsubnode%p_Rlevels(ilev)%rtempVector)
                  CALL c2d2_spaceTimeMatVec (rsolverNode%p_rproblem, p_rmatrix, &
                      p_rsubnode%p_Rlevels(ilev)%rsolutionVector,&
                      p_rsubnode%p_Rlevels(ilev)%rtempVector, -1.0_DP,1.0_DP,&
                      SPTID_FILTER_DEFECT,dres)

                  IF (.NOT.((dres .GE. 1E-99_DP) .AND. &
                            (dres .LE. 1E99_DP))) dres = 0.0_DP
                            
                  CALL output_line ('Space-Time-Multigrid: Level '//TRIM(sys_siL(ilev,5))//&
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

                  IF (((p_rsubnode%p_Rlevels(ilev)%depsRelCycle .NE. 1E99_DP) .OR. &
                       (p_rsubnode%p_Rlevels(ilev)%depsRelCycle .NE. 1E99_DP)) .AND. &
                      (ilev .LT. NLMAX)) THEN
                      
                    ! Adaptive cycles activated. 
                    !
                    ! We are on a level < nlmax.
                    ! At first, calculate the residuum on that level.
                    CALL sptivec_copyVector (p_rsubnode%p_Rlevels(ilev)%rrhsVector,&
                                            p_rsubnode%p_Rlevels(ilev)%rtempVector)
                    CALL c2d2_spaceTimeMatVec (rsolverNode%p_rproblem, p_rmatrix, &
                        p_rsubnode%p_Rlevels(ilev)%rsolutionVector,&
                        p_rsubnode%p_Rlevels(ilev)%rtempVector, -1.0_DP,1.0_DP,&
                        SPTID_FILTER_DEFECT,dres)

                    IF (.NOT.((dres .GE. 1E-99_DP) .AND. &
                              (dres .LE. 1E99_DP))) dres = 0.0_DP
                              
                    ! Compare it with the initial residuum. If it's not small enough
                    ! and if we haven't reached the maximum number of cycles,
                    ! repeat the complete cycle.
                    IF ( ((p_rsubnode%p_Rlevels(ilev)%nmaxAdaptiveCycles .LE. -1) &
                          .OR. &
                          (p_rsubnode%p_Rlevels(ilev)%icycleCount .LT. &
                           p_rsubnode%p_Rlevels(ilev)%nmaxAdaptiveCycles)) &
                        .AND. &
                        ((dres .GT. p_rsubnode%p_Rlevels(ilev)%depsRelCycle * &
                                   p_rsubnode%p_Rlevels(ilev)%dinitResCycle) .AND. &
                         (dres .GT. p_rsubnode%p_Rlevels(ilev)%depsAbsCycle) ) ) THEN

                      IF (rsolverNode%ioutputLevel .GE. 3) THEN
                        CALL output_line ( &
                          TRIM( &
                          sys_siL(p_rsubnode%p_Rlevels(ilev)%icycleCount,10)) &
                          //'''th repetition of cycle on level '// &
                          TRIM(sys_siL(ilev,5))//'.')
                      END IF

                      p_rsubnode%p_Rlevels(ilev)%icycleCount = &
                          p_rsubnode%p_Rlevels(ilev)%icycleCount+1
                      CYCLE cycleloop
                    END IF
                    
                    ! Otherwise: The cycle(s) is/are finished; 
                    ! the END DO goes up one level.
                    
                  END IF

                END IF
              ELSE

                IF (rsolverNode%ioutputLevel .GE. 3) THEN
                  CALL output_line ('Space-Time-Multigrid: Cycle on level '&
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
              p_rsubnode%p_Rlevels(ilev)%rtempVector, -1.0_DP,1.0_DP,&
              SPTID_FILTER_DEFECT,dres)

          IF (.NOT.((dres .GE. 1E-99_DP) .AND. &
                    (dres .LE. 1E99_DP))) dres = 0.0_DP
          
          rsolverNode%dfinalDefect = dres
          
          ! Test if the iteration is diverged
          IF (sptils_testDivergence(rsolverNode,dres)) THEN
            CALL output_line ('Space-Time-Multigrid: Solution diverging!')
            rsolverNode%iresult = 1
            EXIT
          END IF

          ! At least perform nminIterations iterations
          IF (ite .GE. nminIterations) THEN
          
            ! Check if the iteration converged
            IF (sptils_testConvergence(rsolverNode,dres)) EXIT
            
          END IF

          ! print out the current residuum

          IF ((rsolverNode%ioutputLevel .GE. 2) .AND. &
              (MOD(ite,niteResOutput).EQ.0)) THEN
            CALL output_line ('Space-Time-Multigrid: Iteration '// &
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
          CALL output_line ('Space-Time-Multigrid: Iteration '// &
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
        CALL output_line ('Space-Time-Multigrid statistics:')
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
    END IF  
  
  CONTAINS
  
    ! ---------------------------------------------------------------
  
    SUBROUTINE calcCGCenergyMin (rproblem,&
                    rsolutionVector,rcorrectionVector,rrhsVector,rmatrix,&
                    rtempVector,rtempSpaceVector,dalphamin, dalphamax, dstep)
                  
    ! Calculates the optimal coarse grid correction factor using
    ! energy minimisation.
    
    ! Problem structure
    TYPE(t_problem), INTENT(INOUT) :: rproblem
    
    ! Uncorrected solution vector
    TYPE(t_spacetimeVector), INTENT(IN) :: rsolutionVector
    
    ! Correction vector
    TYPE(t_spacetimeVector), INTENT(IN) :: rcorrectionVector

    ! RHS vector
    TYPE(t_spacetimeVector), INTENT(IN) :: rrhsVector
    
    ! Corresponding space-time matrix
    TYPE(t_ccoptSpaceTimeMatrix), INTENT(IN) :: rmatrix
    
    ! Temp vector
    TYPE(t_spacetimeVector), INTENT(INOUT) :: rtempVector
    
    ! A temp vector, only in space
    TYPE(t_vectorBlock), INTENT(INOUT) :: rtempSpaceVector
    
    ! Minimum/maximum value for the correction factor
    REAL(DP), INTENT(IN) :: dalphamin,dalphamax
    
    ! OUTPUT: Calculated step length parameter
    REAL(DP), INTENT(OUT) :: dstep
    
      ! local variables
      REAL(DP) :: dnom, ddenom
    
      ! The formula is just as for standard multigrid.
      ! With u=uncorrected solution, c=correction vector, 
      ! f=rhs vector, calculate
      !
      !           ( f-Au, c )
      !  dstep = -------------
      !            ( Ac, c )
      !
      ! Calculate the nominator
      CALL sptivec_copyVector (rrhsVector,rtempVector)
      CALL c2d2_spaceTimeMatVec (rproblem, rmatrix, &
          rsolutionVector,rtempVector, -1.0_DP,1.0_DP,SPTID_FILTER_DEFECT)

      ! Implement boundary conditions into the vector.
      !CALL tbc_implementInitCondDefect (&
      !    rmatrix%p_rspaceTimeDiscretisation,&
      !    rtempVector, rtempSpaceVector)
      !CALL tbc_implementBCdefect (rproblem,&
      !    rmatrix%p_rspaceTimeDiscretisation,&
      !    rtempVector, rtempSpaceVector)

      ! Get the nominator      
      dnom = sptivec_scalarProduct (rtempVector,rcorrectionVector)

      ! Calculate the denominator:
      
      CALL sptivec_clearVector (rtempVector)
      CALL c2d2_spaceTimeMatVec (rproblem, rmatrix, &
          rcorrectionVector,rtempVector, 1.0_DP,0.0_DP,SPTID_FILTER_DEFECT)
      
      ! Implement boundary conditions into the vector.
      !CALL tbc_implementInitCondDefect (&
      !    rmatrix%p_rspaceTimeDiscretisation,&
      !    rtempVector, rtempSpaceVector)
      !CALL tbc_implementBCdefect (rproblem,&
      !    rmatrix%p_rspaceTimeDiscretisation,&
      !    rtempVector, rtempSpaceVector)

      ! Get the denominator      
      ddenom = sptivec_scalarProduct (rtempVector,rcorrectionVector)
      
      ! Trick to avoid div/0
      IF (ddenom .EQ. 0.0_DP) ddenom = 1.0_DP
      
      ! Result...
      dstep = dnom / ddenom
                   
      ! Force it to be in the allowed range.
      dstep = MIN(dalphamax,MAX(dalphamin,dstep))
                    
    END SUBROUTINE
  
  END SUBROUTINE

END MODULE
