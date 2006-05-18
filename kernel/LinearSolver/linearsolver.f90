!##############################################################################
!# ****************************************************************************
!# <name> linearsolver </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains a collection of linear solvers of various types.
!# The solvers are highly integrated and allow to be nested into each
!# other in a very general order. On the other hand, the solvers are
!# kept "black box like", so they can be configured from outside for
!# a large varienty of problems.
!#
!# The solvers basically deal with block systems of the following form:
!#
!# 1.) General linear systems of the form
!#         $$Ax = b$$
!#     with A being an arbitrary block matrix
!# 2.) Saddle-point block matrices of the form
!#           [A  B1] (x) = (f1)
!#           [B2  0] (y)   (f2)
!#     with a regular block matrix A and general block matrices B1 and B2
!# 3.) Saddle-point-like block matrices of the form
!#           [A  B1] (x) = (f1)
!#           [B2  C] (y)   (f2)
!#     with a regular block matrix A, general block matrices B1 and B2
!#     and a 'stabilisation' matrix $C \sim 0$.
!#
!# Cases 2.) + 3.) are special cases of 1.) in that sense, that specialised
!# solvers are available for these kind of systems which allow faster solving.
!#
!# To solve a problem, one has basically to call the following routines
!# 1.) linsol_initXXXX                - Initialises a solver, returns a solver
!#                                      structure identifying the solver
!# 2.) linsol_setMatrices             - Attach the system matrix/matrices to the
!#                                      solver
!# 3.) linsol_initStructure           - Allow the solvers to perform problem-
!#                                      structure specific initialisation
!#                                      (e.g. symbolical factorisation).
!#                                      During this phase, each solver allocates
!#                                      temporary memory it needs for the solution
!#                                      process
!# 3.) linsol_initData                - Allow the solvers to perform problem-
!#                                      data specific initialisation
!#                                      (e.g. numerical factorisation)
!# 4.) linsol_performSolve            - Solve the problem for a defect vector
!# or  linsol_performSolveAdaptively  - Solve the problem with an initial 
!#                                      solution vector
!# 5.) linsol_doneData                - Release problem-data specific information
!# 6.) linsol_doneStructure           - Release problem-structure specific 
!#                                      information. Release temporary memory.
!# 7.) linsol_releaseSolver           - Clean up solver structures, remove the
!#                                      solver structure from the heap.
!# </purpose>
!##############################################################################

MODULE linearsolver

  USE fsystem
  USE linearsystemblock
  USE filtersupport
  
  IMPLICIT NONE

! *****************************************************************************
! *****************************************************************************
! *****************************************************************************

!<constants>

!<constantblock description="Algorithm identifiers">

  ! Undefined algorithm
  INTEGER, PARAMETER :: LINSOL_ALG_UNDEFINED     = 0
  
  ! Richardson iteration $x_{n+1} = x_n + \omega(b-Ax)$
  INTEGER, PARAMETER :: LINSOL_ALG_RICHARDSON    = 1
  
  ! Jacobi iteration $x_{n+1} = x_n + \omega D^{-1}(b-Ax)$
  INTEGER, PARAMETER :: LINSOL_ALG_JACOBI        = 2

  ! Preconditioned defect correction $x_{n+1} = x_n + \omega P^{-1}(b-Ax)$
  INTEGER, PARAMETER :: LINSOL_ALG_DEFCOR        = 3
  
  ! SOR/GS iteration $x_{n+1} = x_n + (L+\omega D)^{-1}(b-Ax)$
  INTEGER, PARAMETER :: LINSOL_ALG_SOR           = 4
  
  ! SSOR iteration
  INTEGER, PARAMETER :: LINSOL_ALG_SSOR          = 5
  
  ! CG iteration (preconditioned) 
  INTEGER, PARAMETER :: LINSOL_ALG_CG            = 6
  
  ! BiCGStab iteration (preconditioned) 
  INTEGER, PARAMETER :: LINSOL_ALG_BICGSTAB      = 7

  ! GMRES iteration (preconditioned) 
  INTEGER, PARAMETER :: LINSOL_ALG_GMRES         = 8
  
  ! Multigrid iteration
  INTEGER, PARAMETER :: LINSOL_ALG_MULTIGRID     = 9

  ! UMFPACK2
  INTEGER, PARAMETER :: LINSOL_ALG_UMFPACK2      = 10

  ! UMFPACK4
  INTEGER, PARAMETER :: LINSOL_ALG_UMFPACK4      = 11
  
  ! ILU(0) iteration (scalar system)
  INTEGER, PARAMETER :: LINSOL_ALG_ILU0          = 50
  
  ! (M)ILU(s) iteration (scalar system)
  INTEGER, PARAMETER :: LINSOL_ALG_MILUS         = 51
  
  ! SPAI(0) iteration (scalar system)
  INTEGER, PARAMETER :: LINSOL_ALG_SPAI0         = 52

  ! SPAI(k) iteration (scalar system)
  INTEGER, PARAMETER :: LINSOL_ALG_SPAIK         = 53
  
  ! VANCA iteration (general)
  INTEGER, PARAMETER :: LINSOL_ALG_VANCA         = 54
  
  ! VANCA iteration (2D Nav.St., pure $\tilde Q_1/P_0$ discretisation)
  INTEGER, PARAMETER :: LINSOL_ALG_VANCAQ1TP02DNS = 55
  
!</constantblock>

!<constantblock description="Identifiers for stopping criteria">

  ! Use standard stopping criterion.
  ! If depsRel>0: use relative stopping criterion.
  ! If depsAbs>0: use abs stopping criterion.
  ! If both are > 0: use both, i.e. stop if both criteria hold
  INTEGER, PARAMETER :: LINSOL_STOP_STANDARD     = 0
  
!</constantblock>

! *****************************************************************************

!<constantblock description="Bitfield identifiers for the ability of a linear solver">

  ! Solver can handle scalar systems
  INTEGER, PARAMETER :: LINSOL_ABIL_SCALAR       = 2**0

  ! Solver can handle block systems
  INTEGER, PARAMETER :: LINSOL_ABIL_BLOCK        = 2**1

  ! Solver can handle multiple levels
  INTEGER, PARAMETER :: LINSOL_ABIL_MULTILEVEL   = 2**2
  
  ! Solver allows checking the defect during the iteration.
  ! Solvers not capable of this perform only a fixed number of solution
  ! steps (e.g. UMFPACK performs always one step).
  INTEGER, PARAMETER :: LINSOL_ABIL_CHECKDEF     = 2**3
  
  ! Solver is a direct solver (e.g. UMFPACK, ILU).
  ! Otherwise the solver is of iterative nature and might perform
  ! multiple steps to solve the problem.
  INTEGER, PARAMETER :: LINSOL_ABIL_DIRECT       = 2**4
  
  ! Solver might use subsolvers (preconditioners, smoothers,...)
  INTEGER, PARAMETER :: LINSOL_ABIL_USESUBSOLVER = 2**5
  
  ! Solver supports filtering
  INTEGER, PARAMETER :: LINSOL_ABIL_USEFILTER    = 2**6
  
!</constantblock>

!</constants>

! *****************************************************************************
! *****************************************************************************
! *****************************************************************************

!<types>
  
!<typeblock>
  
  ! Object that gives the connection to the linear system and the 
  ! discretisation. Contains stuff like pointers to system matrices,
  ! and other level-dependent information, which does not change
  ! during the solution process.
  
  TYPE t_LinearSystemInfo
    
    ! Pointer to the system matrix. This is always of block-type
    TYPE(t_matrixBlock), POINTER            :: p_rsystemMatrix    => NULL()
    
  END TYPE
  
!</typeblock>

! *****************************************************************************

!<typeblock>
  
  ! This is the central structure which repesents a solver.
  ! It collects all information for the solution process.
  !
  ! The caller who wants to solve a problem has to initialise
  ! such a structure with all information that is necessary
  ! to solve the linear system. For more complicated solvers
  ! which contain sub-solvers like preconditioners or smoothers
  ! (like in MG), such a structure must be created also for the
  ! sub-solver and attached to the main solver. This way, the
  ! structures form a ''solver tree'' with one solver-structure
  ! forming the root and the other structures, which are attached
  ! to this, forming the branches/nodes.
  !
  ! The structure contains general solver parameters which appear is
  ! many solvers. Whether or not a special parameter is respected
  ! by a solver depends on the solver (e.g. UMFPACK won't respect
  ! the ''maximum iterations'' parameter). Further solver-specific
  ! parameters which fit not into the general scheme are kept in
  ! solver specific structures. 
  !
  ! The structure contains INPUT, OUTPUT and STATUS parameters. INPUT
  ! parameters can be changed by the caller prior to solving.
  ! OUTPUT parameters indicate the final status of the solution process.
  ! STATUS parameters are only valid during the solution
  ! process. They should not be accessed from outside, as they are
  ! maintained internally by the solver.
  ! Parameters marked as READ ONLY are set during the initialisation phase
  ! of the solver and must not be changed by the caller.
  
  TYPE t_linsolNode
    
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
    ! General solver parameter; solver specific use
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
    REAL(DP)                        :: ddivRel = SYS_INFINITY

    ! INPUT PARAMETER FOR ITERATIVE SOLVERS: 
    ! Absolute divergence criterion.  Treat iteration as
    ! diverged if
    !   !!defect!! >= DIVREL
    ! A value of SYS_INFINITY disables the absolute divergence check.
    ! standard = SYS_INFINITY
    REAL(DP)                        :: ddivAbs = SYS_INFINITY

    ! INPUT PARAMETER FOR ITERATIVE SOLVERS: 
    ! RHS-vector is treated as zero if max(defect) < drhsZero
    REAL(DP)                        :: drhsZero = 1E-12_DP

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
    
    ! INPUT PARAMETER FOR ITERATIVE SOLVERS: 
    ! Number of iterations that should be used to calculate
    ! the asymptotic convergtence rate, if the solver supports
    ! to generate an asymptotic convergence rate of the last
    ! couple of steps. =0: don't determine asymptotic convergence rate 
    INTEGER                    :: niteAsymptoticCVR = 0
    
    ! INPUT PARAMETER: Output level
    ! This determines the output level of the solver.
    ! =0: no output, =1: basic output, =2, extended output
    INTEGER                    :: ioutputLevel = 2

    ! INPUT PARAMETER FOR ITERATIVE SOLVERS WITH RESIDUAL CHECK:
    ! Number of iterations to perform before printing out the
    ! norm of the residual to screen.
    ! =1: Print residual in every iteration
    INTEGER                    :: niteResOutput = 1

    ! INPUT PARAMETER: Solver subgroup
    ! By default, every solver in a soler tree belongs to solver
    ! subgroup 0. This means, the solver is initialised  the default
    ! call to linsol_initProblemStructure and linsol_initProblemData.
    ! By assigning a solver to a different subgroup , the solver will
    ! not initialise itself by a call to these routines. Instead,
    ! it will initialise itself if these routines are called for
    ! the specified subgroup number.
    ! Example: Multigrid solver, UMFPACK coarse grid solver,
    !  ILU(3) smoother.
    !  MG-Solver, UMFPACK solver: subgroup number := 0 (default)
    !  ILU(3) smoother:           subgroup number := 1
    ! linsol_initProblemStructure and linsol_initProblemData will
    ! initialise only MG and UMFPACK, not ILU(3). Instead,
    ! ILU(3) is initialised if linsol_initProblemStructure and 
    ! linsol_initProblemData are called with 'isolverSubgroup=1'.
    ! So, this mechanism provides a possibility to prevent some
    ! solvers from being initialised or allows to initialise
    ! some solvers manually.
    INTEGER                     :: isolverSubgroup = 0
    
    ! INPUT PARAMETER: t_LinearSystemInfo structure that holds
    ! information about the linear system (system matrix,...)
    ! level where to solve (i.e. usually the finest level). 
    ! For multilevel-capable algorithms, information about the other 
    ! levels are saved in local solver-specific structures.
    TYPE(t_LinearSystemInfo)          :: rlinearSystemInfo           

    ! READ ONLY: Algorithm identifier. 
    ! One of the LINSOL_ALG_xxxx constants.
    ! Depending on the value, a solver-specific structure may be 
    ! assigned to this structure below.
    INTEGER                    :: calgorithm = LINSOL_ALG_UNDEFINED
    
    ! READ ONLY: Solver ability tag. 
    ! Bitfield. A combination of LINSOL_ABIL_xxxx
    ! flags that specify the ability of the actual solver (e.g.
    ! whether it can handle block-matrices or only scalar matrices,...).
    ! Calling a solver that is not able to handle a specific problem
    ! will cause an error by the framework.
    INTEGER(I32)                    :: ccapability = 0
    
    ! STATUS FOR ITERATIVE SOLVERS: Current iteration
    INTEGER                    :: icurrentIteration
    
    ! Pointer to a structure for the BiCGStab solver; NULL() if not set
    TYPE (t_linsolSubnodeBiCGStab), POINTER       :: p_rsubnodeBiCGStab    => NULL()

    ! Pointer to a structure for the UMFPACK4 solver; NULL() if not set
    TYPE (t_linsolSubnodeUMFPACK4), POINTER       :: p_rsubnodeUMFPACK4    => NULL()

    ! Pointer to a structure for the VANCA-Q1TP0NS2D solver; NULL() if not set
    TYPE (t_linsolSubnodeVANCAQ1TP0NS2D), POINTER :: p_rsubnodeVANCAQ1TP0NS2D => NULL()
    
    ! Pointer to a structure for the Multigrid solver; NULL() if not set
    TYPE (t_linsolSubnodeMultigrid), POINTER      :: p_rsubnodeMultigrid   => NULL()

    ! Pointer to a structure for the ILU0 1x1 solver; NULL() if not set
    TYPE (t_linsolSubnodeILU01x1), POINTER        :: p_rsubnodeILU01x1     => NULL()

  END TYPE
  
!</typeblock>
  
! *****************************************************************************

!<typeblock>
  
  ! This structure realises the subnode for the scalar ILU(0) solver in
  ! FEAT style.
  
  TYPE t_linsolSubnodeILU01x1
  
    ! A pointer to the ILU0-decomposition of the main matrix.
    ! Shares the structure with the system matrix.
    TYPE(t_matrixScalar)              :: p_DiluDecomposition      
    
  END TYPE
  
!</typeblock>
  
! *****************************************************************************

!<typeblock>
  
  ! This structure realises the subnode for the BiCGStab solver.
  ! The entry p_rpreconditioner points either to NULL() or to another
  ! t_linsolNode structure for the solver that realises the 
  ! preconditioning.
  
  TYPE t_linsolSubnodeBiCGStab
  
    ! Temporary vectors to use during the solution process
    TYPE(t_vectorBlock), DIMENSION(6) :: RtempVectors

    ! A pointer to the solver node for the preconditioner or NULL(),
    ! if no preconditioner is used.
    TYPE(t_linsolNode), POINTER       :: p_rpreconditioner            => NULL()
    
    ! A pointer to a filter chain, as this solver supports filtering.
    ! The filter chain must be configured for being applied to defect vectors.
    TYPE(t_filterChain), DIMENSION(:), POINTER      :: p_RfilterChain => NULL()
  
  END TYPE
  
!</typeblock>

! *****************************************************************************

!<typeblock>
  
  ! This structure realises the subnode for the UMFPACK4 solver.
  
  TYPE t_linsolSubnodeUMFPACK4
  
    ! Control structure for UMFPACK4; contains parameter for the solver
    REAL(DP), DIMENSION(20) :: Control

    ! Status variables of UMFPACK4; receives the UMFPACK-specific return code
    ! of a call to the solver routines.
    REAL(DP), DIMENSION(90) :: Info

    ! Handle for symbolic factorisation
    INTEGER(I32) :: symbolic = ST_NOHANDLE

    ! Handle for numeric factorisation
    INTEGER(I32) :: numeric = ST_NOHANDLE

  END TYPE
  
!</typeblock>

! *****************************************************************************

!<typeblock>
  
  ! This structure realises the subnode for the VANCA solver.
  ! Specialised VANCA, $\tilde Q_1/P_0$, 2D Navier Stokes.
  
  TYPE t_linsolSubnodeVANCAQ1TP0NS2D

    INTEGER DUMMY  

  END TYPE
  
!</typeblock>

! *****************************************************************************

!<typeblock>

  ! This structure collects timing information for the multigrid solver.

  TYPE t_linsolMGTiming
    
    ! Total time for prolongation
    REAL(DP)                           :: d_tmProlongation

    ! Total time for restriction
    REAL(DP)                           :: d_tmRestriction
  
    ! Total time for defect calculation
    REAL(DP)                           :: d_tmDefCalc
    
    ! Total time for smoothing
    REAL(DP)                           :: d_tmSmooth
    
    ! Total time for coarse grid solving
    REAL(DP)                           :: d_tmCoarseGridSolve

    ! Total time for coarse grid corrections
    REAL(DP)                           :: d_tmCorrection
  
  END TYPE

!</typeblock>

! *****************************************************************************

!<typeblock>
  
  ! Level data for multigrid. This structure forms an entry in a linked list
  ! of levels which multigrid uses for solving the given problem.
  ! The solver subnode t_linsolSubnodeMultigrid holds pointers to the head
  ! and the tail of the list.
  
  TYPE t_linsolMGLevelInfo
  
    ! A level identifier. Can be set to identify the global number of the
    ! level, this structure refers to. Only for output purposes.
    INTEGER                        :: ilevel                 = 0
    
    ! Pointer to previous (lower) level or NULL() if the node is the 
    ! head of the list
    TYPE(t_linsolMGLevelInfo), POINTER  :: p_rprevLevel           => NULL()

    ! Pointer to next (higher) level or NULL() if the node is the 
    ! tail of the list
    TYPE(t_linsolMGLevelInfo), POINTER  :: p_rnextLevel           => NULL()
  
    ! t_LinearSystemInfo structure that holds information
    ! about the linear system (system matrices,...)
    TYPE(t_LinearSystemInfo)            :: rlinearSystemInfo
    
    ! A temporary vector of the same size and structure as the solution vector
    ! corresponding to this level. The memorya for this is allocated
    ! in initStructure and released in doneStructure.
    TYPE(t_vectorBlock)                 :: rtemp

    ! A pointer to the solver node for the presmoothing or NULL(),
    ! if no presmoother is used.
    TYPE(t_linsolNode), POINTER         :: p_rpresmoother         => NULL()

    ! A pointer to the solver node for the postsmoothing or NULL(),
    ! if no presmoother is used.
    TYPE(t_linsolNode), POINTER         :: p_rpostsmoother        => NULL()

    ! A pointer to the solver node for the coarse grid solver if
    ! the level corresponding to this structure represents the coarse
    ! grid. NULL() otherwise.
    TYPE(t_linsolNode), POINTER         :: p_rcoarseGridSolver    => NULL()
    
    ! A pointer to a filter chain, as this solver supports filtering.
    ! The filter chain must be configured for being applied to defect vectors.
    TYPE(t_filterChain), DIMENSION(:),POINTER      :: p_RfilterChain => NULL()
    
    ! STATUS/INTERNAL: MG cycle information.
    ! Number of cycles to perform on this level.
    INTEGER                        :: ncycles

    ! STATUS/INTERNAL: MG cycle information.
    ! Number of remaining cycles to perform on this level.
    INTEGER                        :: ncyclesRemaining
    
  END TYPE
  
!</typeblock>

  ! ***************************************************************************

!<typeblock>
  
  ! This structure realises the subnode for the Multigrid solver.
  ! There are pointers to the head and the tail of a linked list of
  ! t_mgLevelInfo structures which hold the data of all levels.
  ! The tail of the list corresponds to the maximum level where
  ! multigrid solves the problem.
  
  TYPE t_linsolSubnodeMultigrid
  
    ! Cycle identifier. 0=F-cycle, 1=V-cycle, 2=W-cycle
    INTEGER                       :: icycle                   = 0
    
    ! Number of levels in the linked list of multigrid levels
    INTEGER                       :: nlevels                  = 0
    
    ! Pointer to the head of the linked list of levels; corresponds
    ! to the lowest level.
    TYPE(t_linsolMGLevelInfo), POINTER :: p_rlevelInfoHead         => NULL()

    ! Pointer to the tail of the linked list of levels; corresponds
    ! to the highest level, i.e. the level where the system should be
    ! solved.
    TYPE(t_linsolMGLevelInfo), POINTER :: p_rlevelInfoTail         => NULL()
    
    ! Minimum step length for optimal coarse grid correction
    REAL(DP)                           :: dminStep                 = 1.0_DP

    ! Maximum step length for optimal coarse grid correction
    REAL(DP)                           :: dmaxStep                 = 1.0_DP
    
    ! A pointer to a filter chain, as this solver supports filtering.
    ! The filter chain must be configured for being applied to defect vectors.
    TYPE(t_filterChain), DIMENSION(:),POINTER :: p_RfilterChain     => NULL()
  
  END TYPE
  
!</typeblock>

!</types>

  INTERFACE linsol_setMatrices
    MODULE PROCEDURE linsol_setMatricesIndirect
    MODULE PROCEDURE linsol_setMatricesDirect
    MODULE PROCEDURE linsol_setOnelevelMatrixDirect
  END INTERFACE

! *****************************************************************************
! *****************************************************************************
! *****************************************************************************

CONTAINS

  ! ***************************************************************************

  ! ***************************************************************************
  
!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_setMatricesIndirect (rsolverNode,RlinsysInfo)
  
!<description>
  
  ! Initialises the system matrices in the solver node rsolverNode and
  ! in all nodes attached to it (preconditioners,...). For this purpose,
  ! the corresponding initialisation routine of each solver is called.
  
!</description>
  
!<input>
  
  ! Array of pointers to system matrices on all levels of the discretisation.
  ! This is passed through all initialisation routines, but actually used 
  ! only by the multigrid initialisation routine.
  TYPE(t_LinearSystemInfo), DIMENSION(:), INTENT(IN) :: RlinsysInfo
  
!</input>
  
!<inputoutput>
  
  ! The solver node which should be initialised
  TYPE(t_linsolNode), INTENT(INOUT)                     :: rsolverNode
  
!</inputoutput>
  
!</subroutine>

  ! Copy the matrix structure on the finest level to the rsolverNode
  ! structure. This corresponds to the system we want to solve.
  ! (Remark: This copies the pointers in the structures, not the
  !  data the pointer points to!)
  
  rsolverNode%rlinearSystemInfo = RlinsysInfo(UBOUND(RlinsysInfo,1))
  
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
  CASE (LINSOL_ALG_UMFPACK4)
    ! UMFPACK needs no matrix initialisation routine, as it does not
    ! contain subsolvers. An attached matrix is processed in the 
    ! symbolical/numerical factorisation.
  CASE (LINSOL_ALG_VANCAQ1TP02DNS)
    ! VANCA needs no matrix initialisation routine, as it does not
    ! contain subsolvers. An attached matrix is processed in the 
    ! symbolical/numerical factorisation.
  CASE (LINSOL_ALG_BICGSTAB)
    CALL linsol_setMatrixBiCGStab (rsolverNode, &
                                    RlinsysInfo(UBOUND(RlinsysInfo,1)),RlinsysInfo)
  CASE (LINSOL_ALG_MULTIGRID)
    CALL linsol_setMatrixMultigrid (rsolverNode, RlinsysInfo)
  CASE DEFAULT
  END SELECT

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>
  
  SUBROUTINE linsol_setMatricesDirect (rsolverNode,Rmatrices)
  
!<description>
  
  ! Initialises the system matrices in the solver node rsolverNode and
  ! in all nodes attached to it (preconditioners,...). For this purpose,
  ! the corresponding initialisation routine of each solver is called.
  !
  ! Warning: This routine should be used with care! The solver will save
  !   pointers to the matrices in the matrix array Rmatrices. This array
  !   of matrices must be available until the system is solved!!!
  !   In contrast, using linsol_setMatricesIndirect allows to attach
  !   matrices anywhere in memory, without the restriction of being
  !   restricted to an array of t_matrixBlock structures!
  
!</description>
  
!<input>
  
  ! Array of pointers to system matrices on all levels of the discretisation.
  ! This is passed through all initialisation routines, but actually used 
  ! only by the multigrid initialisation routine.
  TYPE(t_matrixBlock), DIMENSION(:), INTENT(IN),TARGET :: Rmatrices
  
!</input>
  
!<inputoutput>
  
  ! The solver node which should be initialised
  TYPE(t_linsolNode), INTENT(INOUT)                     :: rsolverNode
  
!</inputoutput>

!</subroutine>

  ! local variables
  
  ! t_linearSolverInfo structure array that receives the matrix pointers
  TYPE(t_linearSystemInfo), DIMENSION(SIZE(Rmatrices)) :: RlinsysInfo
  INTEGER :: i
  
  ! Copy the matrix pointers.
  DO i=1,SIZE(Rmatrices)
    RlinsysInfo(i)%p_rsystemMatrix => Rmatrices(i)
  END DO
  
  ! The rest of RlinsysInfo(.) will stay uninitialised - ok, there's no rest
  ! up to now :-)
  
  ! Call the standard initialisation routine
  CALL linsol_setMatricesIndirect (rsolverNode,RlinsysInfo)

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>
  
  SUBROUTINE linsol_setOnelevelMatrixDirect (rsolverNode,rmatrix)
  
!<description>
  
  ! This is a special-case subroutine for the case that there is only
  ! one level to solve (e.g. when solving directly with UMFPACK).
  ! rmatrix is exactly one matrix of the level where to solve.
  
!</description>
  
!<input>
  
  ! Array of pointers to system matrices on all levels of the discretisation.
  ! This is passed through all initialisation routines, but actually used 
  ! only by the multigrid initialisation routine.
  TYPE(t_matrixBlock), INTENT(IN),TARGET :: rmatrix
  
!</input>
  
!<inputoutput>
  
  ! The solver node which should be initialised
  TYPE(t_linsolNode), INTENT(INOUT)                     :: rsolverNode
  
!</inputoutput>

!</subroutine>

  ! local variables
  
  ! t_linearSolverInfo structure array that receives the matrix pointer;
  ! it has length 1 here in this special situation.
  TYPE(t_linearSystemInfo), DIMENSION(1) :: RlinsysInfo
  
  ! Copy the matrix pointer.
  RlinsysInfo(1)%p_rsystemMatrix => rmatrix
  
  ! The rest of RlinsysInfo(.) will stay uninitialised - ok, there's no rest
  ! up to now :-)
  
  ! Call the standard initialisation routine
  CALL linsol_setMatricesIndirect (rsolverNode,RlinsysInfo)

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_initStructure (rsolverNode, isolverSubgroup)
  
!<description>
  
  ! Initialises the problem structure in the solver node rsolverNode by calling
  ! the initialisation routine of the appropriate solver. The solver
  ! initialisation routine itself can call this procedure to initialise
  ! its sub-solver nodes.
  ! The initialisation of the problem structure allowes the solver component
  ! to perform some 'precalculation', e.g. the UMFPACK4 or ILU solver can 
  ! perform a symbolical factorisation. The problem structure usually does
  ! not change during a simulation, except when the grid moves e.g..
  
!</description>
  
!<inputoutput>
  
  ! The solver node which should be initialised
  TYPE(t_linsolNode), INTENT(INOUT)                     :: rsolverNode
  
!</inputoutput>

!<input>
    
  ! Optional parameter. isolverSubgroup allows to specify a specific 
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are initialised.
  INTEGER, OPTIONAL, INTENT(IN)                    :: isolverSubgroup
  
!</input>
  
!</subroutine>

  ! local variables
  INTEGER :: isubgroup
  
  ! by default, initialise solver subroup 0
  isubgroup = 0
  IF (PRESENT(isolversubgroup)) isubgroup = isolverSubgroup

  ! Call the structure-init routine of the specific solver
  
  SELECT CASE(rsolverNode%calgorithm)
  CASE (LINSOL_ALG_VANCA)
    ! VANCA: no init-routine.
  CASE (LINSOL_ALG_VANCAQ1TP02DNS)
    ! VANCA: no init-routine.
  CASE (LINSOL_ALG_UMFPACK4)
    CALL linsol_initStructureUMFPACK4 (rsolverNode,isubgroup)
  CASE (LINSOL_ALG_BICGSTAB)
    CALL linsol_initStructureBiCGStab (rsolverNode,isubgroup)
  CASE (LINSOL_ALG_MULTIGRID)
    CALL linsol_initStructureMultigrid (rsolverNode,isubgroup)
  END SELECT
  
  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_initData (rsolverNode, isolverSubgroup)
  
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
  
!<inputoutput>

!<input>
    
  ! Optional parameter. isolverSubgroup allows to specify a specific 
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are initialised.
  INTEGER, OPTIONAL, INTENT(IN)                    :: isolverSubgroup
  
!</input>

  !<inputpoutput>
    
  ! The solver node which should be initialised
  TYPE(t_linsolNode), INTENT(INOUT)                     :: rsolverNode
  
!</inputoutput>

!</subroutine>

  ! local variables
  INTEGER :: isubgroup
  
  ! by default, initialise solver subroup 0
  isubgroup = 0
  IF (PRESENT(isolversubgroup)) isubgroup = isolverSubgroup

  ! Call the data-init routine of the specific solver
  
  SELECT CASE(rsolverNode%calgorithm)
  CASE (LINSOL_ALG_VANCA)
    ! VANCA: no init-routine.
  CASE (LINSOL_ALG_VANCAQ1TP02DNS)
    ! VANCA: no init-routine.
  CASE (LINSOL_ALG_UMFPACK4)
    CALL linsol_initDataUMFPACK4 (rsolverNode,isubgroup)
  CASE (LINSOL_ALG_BICGSTAB)
    CALL linsol_initDataBiCGStab (rsolverNode,isubgroup)
  CASE (LINSOL_ALG_MULTIGRID)
    CALL linsol_initDataMultigrid (rsolverNode,isubgroup)
  END SELECT
  
  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_updateStructure (rsolverNode, isolverSubgroup)
  
!<description>
  
  ! Reinitialises the problem structure in the solver node rsolverNode
  
!</description>
  
!<inputoutput>
  
  ! The solver node which should be reinitialised
  TYPE(t_linsolNode), INTENT(INOUT)                     :: rsolverNode
  
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
  INTEGER :: isubgroup
  
  ! by default, initialise solver subroup 0
  isubgroup = 0
  IF (PRESENT(isolversubgroup)) isubgroup = isolverSubgroup

  ! For now, call the structure-done- and structure-init routines.
  ! Maybe rewritten later for higher performance or special cases...
  
  CALL linsol_doneStructure (rsolverNode,isubgroup)
  CALL linsol_initStructure (rsolverNode,isubgroup)
  
  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_updateData (rsolverNode,isolverSubgroup)
  
!<description>
  
  ! Reinitialises the problem data in the solver node rsolverNode.
  
!</description>
  
!<inputoutput>
  
  ! The solver node containing the solver confuguration
  TYPE(t_linsolNode), INTENT(INOUT)                     :: rsolverNode
  
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
  INTEGER :: isubgroup
  
  ! by default, initialise solver subroup 0
  isubgroup = 0
  IF (PRESENT(isolversubgroup)) isubgroup = isolverSubgroup

  ! For now, call the data-done- and data-init routines.
  ! Maybe rewritten later for higher performance or special cases...
  
  CALL linsol_doneData (rsolverNode, isubgroup)
  CALL linsol_initData (rsolverNode, isubgroup)

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_doneStructure (rsolverNode, isolverSubgroup)
  
!<description>
  
  ! Releases the problem structure in the solver node rsolverNode.
  ! This is done by calling the appropriate routine of the
  ! actual solver.
  
!</description>
  
!<inputoutput>
  
  ! The solver node which should be reinitialised
  TYPE(t_linsolNode), INTENT(INOUT)                     :: rsolverNode
  
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
  INTEGER :: isubgroup
  
  ! by default, handle solver subroup 0
  isubgroup = 0
  IF (PRESENT(isolversubgroup)) isubgroup = isolverSubgroup

  ! Call the data-done routine of the specific solver
  
  SELECT CASE(rsolverNode%calgorithm)
  CASE (LINSOL_ALG_VANCA)
    ! VANCA: no done-routine.
  CASE (LINSOL_ALG_VANCAQ1TP02DNS)
    ! VANCA: no done-routine.
  CASE (LINSOL_ALG_UMFPACK4)
    CALL linsol_doneStructureUMFPACK4 (rsolverNode)
  CASE (LINSOL_ALG_BICGSTAB)
    CALL linsol_doneStructureBiCGStab (rsolverNode)
  CASE (LINSOL_ALG_MULTIGRID)
    CALL linsol_doneStructureMultigrid (rsolverNode)
  END SELECT

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_doneData (rsolverNode,isolverSubgroup)
  
!<description>
  
  ! Releases the problem data in the solver node rsolverNode.
  
!</description>
  
!<inputoutput>
  
  ! The solver node containing the solver confuguration
  TYPE(t_linsolNode), INTENT(INOUT)                     :: rsolverNode
  
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
  INTEGER :: isubgroup
  
  ! by default, handle solver subroup 0
  isubgroup = 0
  IF (PRESENT(isolversubgroup)) isubgroup = isolverSubgroup

  ! Call the data-done routine of the specific solver
  
  SELECT CASE(rsolverNode%calgorithm)
  CASE (LINSOL_ALG_VANCA)
    ! VANCA: no done-routine.
  CASE (LINSOL_ALG_VANCAQ1TP02DNS)
    ! VANCA: no done-routine.
  CASE (LINSOL_ALG_UMFPACK4)
    CALL linsol_doneDataUMFPACK4 (rsolverNode)
  CASE (LINSOL_ALG_BICGSTAB)
    CALL linsol_doneDataBiCGStab (rsolverNode)
  CASE (LINSOL_ALG_MULTIGRID)
    CALL linsol_doneDataMultigrid (rsolverNode)
  END SELECT
  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_releaseSolver (p_rsolverNode,bkeepSolverNode)
  
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
  TYPE(t_linsolNode), POINTER     :: p_rsolverNode
!</inputoutput>
  
!</subroutine>

  IF (.NOT. ASSOCIATED(p_rsolverNode)) THEN
    ! Warning message, return
    PRINT *,'linsol_releaseSolver warning: Solver note not assigned!'
    RETURN
  END IF

  ! Depending on the solver type, call the corresponding done-routine
  SELECT CASE(p_rsolverNode%calgorithm)
  CASE (LINSOL_ALG_VANCA)
    ! VANCA: no done-routine.
  CASE (LINSOL_ALG_VANCAQ1TP02DNS)
    ! VANCA: no done-routine.
  CASE (LINSOL_ALG_UMFPACK4)
    CALL linsol_doneUMFPACK4 (p_rsolverNode)
  CASE (LINSOL_ALG_BICGSTAB)
    CALL linsol_doneBiCGStab (p_rsolverNode)
  CASE (LINSOL_ALG_MULTIGRID)
    CALL linsol_doneMultigrid (p_rsolverNode)
  CASE DEFAULT
  END SELECT
  
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
  
  LOGICAL FUNCTION linsol_testConvergence (rsolverNode, dvecNorm, rdef) RESULT(loutput)
  
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
  TYPE(t_linsolNode), INTENT(IN) :: rsolverNode
  
  ! OPTIONAL: The defect vector which norm should be tested.
  ! If existent, the norm of the vector is returned in dvecNorm.
  ! If not existent, the routine assumes that ddvecNrm is the norm
  ! of the vector and checks convergence depending on dvecNorm.
  TYPE(t_vectorBlock), INTENT(IN), OPTIONAL :: rdef
  
!</input>

!<inputoutput>
  
  ! Norm of the defect vector. 
  ! If rdef if present, the routine will calculate the norm of dx and return
  ! it in dvecNorm.
  ! If rdef is not present, dvecNorm is assumed to be a valid norm of a
  ! vector and convergence is tested using dvecNorm.
  REAL(DP), INTENT(INOUT) :: dvecNorm

!</inputoutput>
  
!</function>

  ! Calculate the norm of the vector or take the one given
  ! as parameter
  IF (PRESENT(rdef)) THEN
    dvecNorm = lsysbl_vectorNorm (rdef,rsolverNode%iresNorm)
  END IF
  
  loutput = .TRUE.
  
  ! Absolute convergence criterion? Check the norm directly.
  IF (rsolverNode%depsAbs .NE. 0.0_DP) THEN
    IF (dvecNorm .GT. rsolverNode%depsAbs) THEN
      loutput = .FALSE.
    END IF
  END IF
  
  ! Relative convergence criterion? Multiply with initial residuum
  ! and check the norm. 
  IF (rsolverNode%depsRel .NE. 0.0_DP) THEN
    IF (dvecNorm .GT. rsolverNode%depsRel * rsolverNode%dinitialDefect) THEN
      loutput = .FALSE.
    END IF
  END IF
  
  END FUNCTION
  
  ! ***************************************************************************

!<function>
  
  LOGICAL FUNCTION linsol_testDivergence (rsolverNode, dvecNorm, rdef) RESULT(loutput)
  
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
  TYPE(t_linsolNode), INTENT(IN) :: rsolverNode
  
  ! OPTIONAL: The defect vector which norm should be tested.
  ! If existent, the norm of the vector is returned in dvecNorm.
  ! If not existent, the routine assumes that ddvecNrm is the norm
  ! of the vector and checks divergence depending on dvecNorm.
  TYPE(t_vectorBlock), INTENT(IN), OPTIONAL :: rdef
  
!</input>

!<inputoutput>
  
  ! Norm of the defect vector. 
  ! If rdef if present, the routine will calculate the norm of dx and return
  ! it in dvecNorm.
  ! If rdef is not present, dvecNorm is assumed to be a valid norm of a
  ! vector and divergence is tested using dvecNorm.
  REAL(DP), INTENT(INOUT) :: dvecNorm

!</inputoutput>
!</function>

  ! Calculate the norm of the vector if not given
  ! as parameter
  IF (PRESENT(rdef)) THEN
    dvecNorm = lsysbl_vectorNorm (rdef,rsolverNode%iresNorm)
  END IF
  
  loutput = .FALSE.
  
  ! Absolute convergence criterion? Check the norm directly.
  IF (rsolverNode%ddivAbs .NE. SYS_INFINITY) THEN
   
    ! use NOT here - gives a better handling of special cases like NaN!
    IF ( .NOT. (dvecNorm .LE. rsolverNode%ddivAbs)) THEN
      loutput = .TRUE.
    END IF
    
  END IF
  
  ! Relative convergence criterion? Multiply with initial residuum
  ! and check the norm. 
  IF (rsolverNode%depsRel .NE. SYS_INFINITY) THEN
    IF ( .NOT. (dvecNorm*rsolverNode%dinitialDefect .LE. &
                rsolverNode%ddivRel) ) THEN
      loutput = .TRUE.
    END IF
  END IF
  
  END FUNCTION
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_performSolve (rsolverNode,rx,rb,rtemp)
  
!<description>
  
  ! This routine starts the solution process to solve a linear system of the
  ! form $Ax=b$. The solver configuration must be given by rsolverNode.
  ! rb is the right hand side vector of the system.
  ! rx is assumed to be =0 and receives the solution vector.
  ! The matrix must have been attached to the system before calling
  ! this routine.
  !
  ! rtemp is a temporary vector of the same size as rx which is used
  ! for intermediate computation. Whether is's used or not dependens
  ! on the actual solver.
  
!</description>
  
!<input>
  
  ! The RHS vector of the system
  TYPE(t_vectorBlock), INTENT(IN), TARGET            :: rb

!</input>
  
!<inputoutput>

  ! The solver node containing the solver configuration
  TYPE(t_linsolNode), INTENT(INOUT)                :: rsolverNode
  
  ! The initial solution vector; receives the solution of the system
  TYPE(t_vectorBlock), INTENT(INOUT)               :: rx

  ! A temporary vector of the same size and structure as rx.
  TYPE(t_vectorBlock), INTENT(INOUT)               :: rtemp
  
!</inputoutput>
  
!</subroutine>

  ! Call the solver routine of the specific solver to calculate
  ! the correction vector.
  ! Use p_rdef as RHS and p_rcorr as solution vector - which
  ! either point to the real RHS/solution or to the defect/correction
  ! vector we just allocated.
  
  SELECT CASE(rsolverNode%calgorithm)
  CASE (LINSOL_ALG_VANCA)
    ! not yet implemented
  CASE (LINSOL_ALG_VANCAQ1TP02DNS)
    CALL linsol_solveVANCAQ1TP0NS2D (rsolverNode,rx,rb,rtemp)
  CASE (LINSOL_ALG_UMFPACK4)
    CALL linsol_solveUMFPACK4 (rsolverNode,rx,rb)
  CASE (LINSOL_ALG_BICGSTAB)
    CALL linsol_solveBiCGStab (rsolverNode,rx,rb,rtemp)
  CASE (LINSOL_ALG_MULTIGRID)
    CALL linsol_solveMultigrid (rsolverNode,rx,rb,rtemp)
  END SELECT

  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE linsol_performSolveAdaptively (rsolverNode,rMatrix,rx,rb,rtemp)
  
!<description>
  
  ! This routine starts the solution process to solve a linear system of the
  ! form $Ax=b$. The solver configuration must be given by rsolverNode.
  ! rb is the right hand side vector of the system. rx is an
  ! initial solution vector, which will be overwritten by the solution.
  !
  ! The solver does not implement any boundary conditions of the real
  ! system into the solution vector, this has to be done by the caller
  ! if necessary!
  !
  ! The matrix must have been attached to the system before calling
  ! this routine.
  !
  ! rtemp is a temporary vector of the same size as rx which is used
  ! for intermediate computation. Whether is's used or not dependens
  ! on the actual solver.
  
!</description>
  
!<input>
  
  ! The RHS vector of the system
  TYPE(t_vectorBlock), INTENT(IN), TARGET            :: rb

  ! The system matrix
  TYPE(t_matrixBlock), INTENT(IN)                    :: rMatrix
  
!</input>
  
!<inputoutput>
  
  ! The solver node containing the solver configuration
  TYPE(t_linsolNode), INTENT(INOUT)                  :: rsolverNode
  
  ! The initial solution vector; receives the solution of the system
  TYPE(t_vectorBlock), INTENT(INOUT)                 :: rx
  
  ! A temporary vector of the same size and structure as rx.
  TYPE(t_vectorBlock), INTENT(INOUT)                 :: rtemp

!</inputoutput>
  
!</subroutine>

  ! Method-specific remarks:
  ! The linear system $Ax=b$ is reformulated into a one-step defect-correction 
  ! approach
  !     $$ x  =  x_0  +  A^{-1}  ( b - A x_0 ) $$
  ! The standard solver above is then used to solve
  !     $$ Ay = b-Ax_0 $$
  ! and the solution $x$ is then calculated by $x=x+y$. 
  
  ! local variables
  
  ! Defect vector, to be build from rb and rx
  TYPE(t_vectorBlock)               :: rdef,rcorr
  
  ! Calculate the defect:
  
  ! Create rdef as temporary vector based on rx.
  CALL lsysbl_createVectorBlock (rx, rdef, .FALSE.)
  
  ! To build (b-Ax), copy the RHS to the temporary vector
  CALL lsysbl_vectorCopy (rb,rdef)
  
  CALL lsysbl_blockMatVec (rMatrix, rx, rdef, 1.0_DP, 1.0_DP)
  
  ! Allocate memory for the correction vector,
  ! Fill it with zero.
  CALL lsysbl_createVectorBlock (rx, rcorr, .TRUE.)
  
  ! Call linsol_performSolve to solve the subproblem $Ay = b-Ax$.
  CALL linsol_performSolve (rsolverNode,rx,rb,rtemp)
  
  ! Add the correction vector to the solution vector and release the memory.
  ! In case we have one... If the initial vector was assumed as zero, we don't
  ! have a correction vector, the result is directly in rx.
  
  ! Correct the solution vector: x=x+y
  CALL lsysbl_vectorLinearComb (rdef,rx,1.0_DP,1.0_DP)

  ! Release memory
  CALL lsysbl_releaseVector (rcorr)
  CALL lsysbl_releaseVector (rdef)
  
  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>
  
  SUBROUTINE linsol_convertToSmoother (rsolverNode,nsmoothingSteps)
  
!<description>
  
  ! Converts a t_linsolNode to a smoother structure. A smoother is a solver
  ! that performs a fixed number of iterations without respecting any
  ! residuum. nsmoothingSteps is the number of steps the smoother should
  ! perform.
  
!</description>
  
!<input>
  
  ! Number of steps the smoother should perform
  INTEGER, INTENT(IN)          :: nsmoothingSteps
  
!</input>
  
!<inputoutput>
  
  ! Solver node which should be configured as smoother.
  TYPE(t_linsolNode), INTENT(INOUT) :: rsolverNode
  
!</inputoutput>
  
!</subroutine>

  rsolverNode%depsRel = 0.0_DP
  rsolverNode%depsAbs = 0.0_DP
  rsolverNode%nminIterations = nsmoothingSteps
  rsolverNode%nmaxIterations = nsmoothingSteps
  rsolverNode%iresCheck      = NO

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>
  
  SUBROUTINE linsol_initSolverGeneral (p_rsolverNode)
  
!<description>
  
  ! Creates a new solver node p_rsolverNode on the heap with default
  ! settings and returns a pointer to it.
  
!</description>
  
!<output>
  
  ! A pointer to a new solver node on the heap.
  TYPE(t_linsolNode), POINTER         :: p_rsolverNode
  
!</output>
  
!</subroutine>

  ! Allocate the node - the default initialisation will set most of the
  ! parameters in the structure.
  ALLOCATE(p_rsolverNode)

  END SUBROUTINE
  
! *****************************************************************************
! Routines for the VANCA CC2D solver
! *****************************************************************************

!<subroutine>
  
  SUBROUTINE linsol_initVANCAQ1TP0NS2D (p_rsolverNode)
  
!<description>
  
  ! Creates a t_linsolNode solver structure for the VANCA solver. The node
  ! can be used to directly solve a problem or to be attached as solver
  ! or preconditioner to another solver structure. The node can be deleted
  ! by linsol_releaseSolver.
  !
  ! This VANCA solver has no done-routine as there is no dynamic information
  ! allocated.
  
!</description>
  
!<output>
  
  ! A pointer to a t_linsolNode structure. Is set by the routine, any previous
  ! value of the pointer is destroyed.
  TYPE(t_linsolNode), POINTER         :: p_rsolverNode
   
!</output>
  
!</subroutine>
  
  ! Create a default solver structure
  
  CALL linsol_initSolverGeneral(p_rsolverNode)
  
  ! Initialise the type of the solver
  p_rsolverNode%calgorithm = LINSOL_ALG_VANCAQ1TP02DNS 
  
  ! Initialise the ability bitfield with the ability of this solver:
  p_rsolverNode%ccapability = LINSOL_ABIL_SCALAR + LINSOL_ABIL_BLOCK + LINSOL_ABIL_CHECKDEF
  
  ! Allocate the subnode for VANCA.
  ! This initialises most of the variables with default values appropriate
  ! to this solver.
  ALLOCATE(p_rsolverNode%p_rsubnodeVANCAQ1TP0NS2D)
  
  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>
  
  SUBROUTINE linsol_solveVANCAQ1TP0NS2D (rsolverNode,rx,rb,rtemp)
  
!<description>
  
  ! Solves the linear system $Ax=b$ with VANCA. The matrix $A$ must be
  ! attached to the solver previously by linsol_setMatrices.
  
!</description>
  
!<input>
  
  ! Right hand side of the system
  TYPE(t_vectorBlock), INTENT(IN)           :: rb
  
!</input>
  
!<inputoutput>
  
  ! The t_linsolNode structure of the UMFPACK4 solver
  TYPE(t_linsolNode), INTENT(INOUT)         :: rsolverNode

  ! Receives the solution vector
  TYPE(t_vectorBlock), INTENT(INOUT)        :: rx
  
  ! A temporary vector of the same size and structure as rx.
  TYPE(t_vectorBlock), INTENT(INOUT)        :: rtemp

!</inputoutput>
  
!</subroutine>

  ! solve the system
  PRINT *,'to be implemented...'
  
  END SUBROUTINE
  
! *****************************************************************************
! Routines for the UMFPACK4 solver
! *****************************************************************************

!<subroutine>
  
  SUBROUTINE linsol_initUMFPACK4 (rsolverNode)
  
!<description>
  
  ! Creates a t_linsolNode solver structure for the UMFPACK4 solver. The node
  ! can be used to directly solve a problem or to be attached as solver
  ! or preconditioner to another solver structure. The node can be deleted
  ! by linsol_releaseSolver.
  
!</description>
  
!<output>
  
  ! A pointer to a t_linsolNode structure. Is set by the routine, any previous
  ! value of the pointer is destroyed.
  TYPE(t_linsolNode), POINTER         :: rsolverNode
   
!</output>
  
!</subroutine>

  ! Create a default solver structure
  
  CALL linsol_initSolverGeneral(rsolverNode)
  
  ! Initialise the type of the solver
  rsolverNode%calgorithm = LINSOL_ALG_UMFPACK4
  
  ! Initialise the ability bitfield with the ability of this solver:
  rsolverNode%ccapability = LINSOL_ABIL_SCALAR + LINSOL_ABIL_DIRECT
  
  ! Allocate the subnode for UMFPACK4.
  ! This initialises most of the variables with default values appropriate
  ! to this solver.
  ALLOCATE(rsolverNode%p_rsubnodeUMFPACK4)
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE linsol_initStructureUMFPACK4 (rsolverNode,isolverSubgroup)
  
!<description>
  
  ! Performs a symbolic factorisation on the assigned matrix.
  
!</description>
  
!<inputoutput>
  
  ! The t_linsolNode structure of the UMFPACK4 solver
  TYPE(t_linsolNode), INTENT(INOUT)         :: rsolverNode
   
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
  INTEGER :: isubgroup
  
  ! by default, initialise solver subroup 0
  isubgroup = 0
  IF (PRESENT(isolversubgroup)) isubgroup = isolverSubgroup

  ! Stop if there's no matrix assigned
  
  IF (.NOT. ASSOCIATED(rsolverNode%rlinearSystemInfo%p_rsystemMatrix)) THEN
    PRINT *,'Error: No matrix associated!'
    STOP
  END IF
  
  ! If isubgroup coincides with isolverSubgroup from the solver
  ! structure, we must initialise - otherwise we simply
  ! skip the initialisation.
  IF (isubgroup .EQ. rsolverNode%isolverSubgroup) THEN
  
    ! perform a symbolic factorization...
    PRINT *,'to be implemented'
    
  END IF
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE linsol_initDataUMFPACK4 (rsolverNode, isolverSubgroup)
  
!<description>
  
  ! Performs a numeric factorisation on the assigned matrix.
  
!</description>
  
!<inputoutput>
  
  ! The t_linsolNode structure of the UMFPACK4 solver
  TYPE(t_linsolNode), INTENT(INOUT)         :: rsolverNode
   
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
  INTEGER :: isubgroup
  
  ! by default, initialise solver subroup 0
  isubgroup = 0
  IF (PRESENT(isolversubgroup)) isubgroup = isolverSubgroup

  ! Stop if there's no matrix assigned
  
  IF (.NOT. ASSOCIATED(rsolverNode%rlinearSystemInfo%p_rsystemMatrix)) THEN
    PRINT *,'Error: No matrix associated!'
    STOP
  END IF
  
  ! If isubgroup coincides with isolverSubgroup from the solver
  ! structure, we must initialise - otherwise we simply
  ! skip the initialisation.
  IF (isubgroup .EQ. rsolverNode%isolverSubgroup) THEN
  
    ! perform a numerical factorization...
    PRINT *,'to be implemented'
    
  END IF
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE linsol_doneDataUMFPACK4 (rsolverNode,isolverSubgroup)
  
!<description>
  
  ! Releases the memory of teh numeric factorisation of the given matrix.
  
!</description>
  
!<inputoutput>
  
  ! The t_linsolNode structure of the UMFPACK4 solver
  TYPE(t_linsolNode), INTENT(INOUT)         :: rsolverNode
   
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
  INTEGER :: isubgroup
  
  ! by default, initialise solver subroup 0
  isubgroup = 0
  IF (PRESENT(isolversubgroup)) isubgroup = isolverSubgroup

  ! If isubgroup coincides with isolverSubgroup from the solver
  ! structure, we must proceed - otherwise we simply
  ! skip the handling.
  IF (isubgroup .EQ. rsolverNode%isolverSubgroup) THEN
  
    ! Release the numerical factorisation
    PRINT *,'to be implemented'
  
  END IF
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE linsol_doneStructureUMFPACK4 (rsolverNode, isolverSUbgroup)
  
!<description>
  
  ! Releases the memory of teh numeric factorisation of the given matrix.
  
!</description>
  
!<inputoutput>
  
  ! The t_linsolNode structure of the UMFPACK4 solver
  TYPE(t_linsolNode), INTENT(INOUT)         :: rsolverNode
   
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
  INTEGER :: isubgroup
  
  ! by default, initialise solver subroup 0
  isubgroup = 0
  IF (PRESENT(isolversubgroup)) isubgroup = isolverSubgroup

  ! If isubgroup coincides with isolverSubgroup from the solver
  ! structure, we must proceed - otherwise we simply
  ! skip the handling.
  IF (isubgroup .EQ. rsolverNode%isolverSubgroup) THEN
    
    ! Release the symbolical factorisation
    PRINT *,'to be implemented'
    
  END IF
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE linsol_doneUMFPACK4 (rsolverNode, isolverSubgroup)
  
!<description>
  
  ! This routine releases all temporary memory for the UMFPACK4 solver from
  ! the heap.
  
!</description>
  
!<inputoutput>
  
  ! The t_linsolNode structure of UMFPACK4 which is to be cleaned up.
  TYPE(t_linsolNode), INTENT(INOUT)         :: rsolverNode
   
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
  INTEGER :: isubgroup
  
  ! by default, initialise solver subroup 0
  isubgroup = 0
  IF (PRESENT(isolversubgroup)) isubgroup = isolverSubgroup

  ! Release symbolical and numerical factorisation if still associated...
  CALL linsol_doneDataUMFPACK4 (rsolverNode, isubgroup)
  CALL linsol_doneStructureUMFPACK4 (rsolverNode, isubgroup)
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE linsol_solveUMFPACK4 (rsolverNode,rx,rb)
  
!<description>
  
  ! Solves the linear system $Ax=b$ with UMFPACK4. The matrix $A$ must be
  ! attached to the solver previously by linsol_setMatrices.
  ! Symbolical and numerical factorisation must already be performed
  ! on the matrix.
  
!</description>
  
!<input>
  
  ! Right hand side of the system
  TYPE(t_vectorBlock), INTENT(IN)           :: rb
  
!</input>
  
!<inputoutput>
  
  ! The t_linsolNode structure of the UMFPACK4 solver
  TYPE(t_linsolNode), INTENT(INOUT)         :: rsolverNode

  ! Receives the solution vector
  TYPE(t_vectorBlock), INTENT(INOUT)        :: rx
  
!</inputoutput>
  
!</subroutine>

  ! solve the system
  PRINT *,'to be implemented...'
  
  END SUBROUTINE
  
! *****************************************************************************
! Routines for the BiCGStab solver
! *****************************************************************************

!<subroutine>
  
  SUBROUTINE linsol_initBiCGStab (p_rsolverNode,p_rpreconditioner,p_Rfilter)
  
!<description>
  
  ! Creates a t_linsolNode solver structure for the BiCGStab solver. The node
  ! can be used to directly solve a problem or to be attached as solver
  ! or preconditioner to another solver structure. The node can be deleted
  ! by linsol_releaseSolver.
  
!</description>
  
!<input>
  
  ! Optional: A pointer to the solver structure of a solver that should be 
  ! used for preconditioning. If not given or set to NULL(), no preconditioning 
  ! will be used.
  TYPE(t_linsolNode), POINTER, OPTIONAL   :: p_rpreconditioner
  
  ! Optional: A pointer to a filter chain (i.e. an array of t_filterChain
  ! structures) if filtering should be applied to the vector during the 
  ! iteration. If not given or set to NULL(), no filtering will be used.
  ! The filter chain (i.e. the array) must exist until the system is solved!
  ! The filter chain must be configured for being applied to defect vectors.
  TYPE(t_filterChain), DIMENSION(:), POINTER, OPTIONAL   :: p_Rfilter
  
!</input>
  
!<output>
  
  ! A pointer to a t_linsolNode structure. Is set by the routine, any previous
  ! value of the pointer is destroyed.
  TYPE(t_linsolNode), POINTER         :: p_rsolverNode
   
!</output>
  
!</subroutine>

  ! Create a default solver structure
  CALL linsol_initSolverGeneral(p_rsolverNode)
  
  ! Initialise the type of the solver
  p_rsolverNode%calgorithm = LINSOL_ALG_BICGSTAB
  
  ! Initialise the ability bitfield with the ability of this solver:
  p_rsolverNode%ccapability = LINSOL_ABIL_SCALAR + LINSOL_ABIL_BLOCK    + &
                              LINSOL_ABIL_CHECKDEF + &
                              LINSOL_ABIL_USESUBSOLVER + &
                              LINSOL_ABIL_USEFILTER
  
  ! Allocate the subnode for BiCGStab.
  ! This initialises most of the variables with default values appropriate
  ! to this solver.
  ALLOCATE(p_rsolverNode%p_rsubnodeBiCGStab)
  
  ! Attach the preconditioner if given. 
  
  IF (PRESENT(p_rpreconditioner)) THEN 
    p_rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner => p_rpreconditioner
  END IF

  ! Attach the filter if given. 
  
  IF (PRESENT(p_Rfilter)) THEN
    p_rsolverNode%p_rsubnodeBiCGStab%p_RfilterChain => p_Rfilter
  END IF
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_setMatrixBiCGStab (rsolverNode,rsystemMatrix, &
                                                  Rmatrices)
  
!<description>
  
  ! This routine is called if the pointer to the system matrix changes.
  ! The routine calls linsol_setMatrices for the preconditioner of BiCGStab
  ! to inform also that one about the change of the matrix pointer.
  
!</description>
  
!<input>
  
  ! The system matrix that is to be assigned to BiCGStab
  TYPE(t_LinearSystemInfo), INTENT(IN), TARGET   :: rsystemMatrix

  ! An array of system matrices which is simply passed to the initialisation 
  ! routine of the preconditioner.
  TYPE(t_LinearSystemInfo), INTENT(IN), DIMENSION(:)   :: Rmatrices

!</input>
  
!<inputoutput>
  
  ! The t_linsolNode structure of the UMFPACK4 solver
  TYPE(t_linsolNode), INTENT(INOUT)         :: rsolverNode
   
!</inputoutput>
  
!</subroutine>

  ! Do we have a preconditioner given? If yes, call the general initialisation
  ! routine to initialise it.
  
  IF (ASSOCIATED(rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner)) THEN
    CALL linsol_setMatrices (rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner, &
                             Rmatrices)
  END IF

  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_initStructureBiCGStab (rsolverNode, isolverSubgroup)
  
!<description>
  
  ! Calls the initStructure subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_initStructure.
  
!</description>
  
!<inputoutput>
  
  ! The t_linsolNode structure of the UMFPACK4 solver
  TYPE(t_linsolNode), INTENT(INOUT)         :: rsolverNode
   
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
  INTEGER :: isubgroup,i
  TYPE(t_linsolSubnodeBiCGStab), POINTER :: p_rsubnode
  
  ! BiCGStab needs 5 temporary vectors + 1 for preconditioning. Allocate that here!
  p_rsubnode => rsolverNode%p_rsubnodeBiCGStab
  DO i=1,6
    CALL lsysbl_createVecBlockIndMat (rsolverNode%rlinearSystemInfo%p_rsystemMatrix, &
         p_rsubnode%RtempVectors(i),.FALSE.)
  END DO
  
  ! by default, initialise solver subroup 0
  isubgroup = 0
  IF (PRESENT(isolversubgroup)) isubgroup = isolverSubgroup

  ! We simply pass isubgroup to the subsolvers when calling them.
  ! Inside of this routine, there's not much to do with isubgroup,
  ! as there is no special information (like a factorisation)
  ! associated with this type of solver, which has to be allocated,
  ! released or updated...

  ! Call the init routine of the preconditioner.
  IF (ASSOCIATED(rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner)) THEN
    CALL linsol_initStructure (rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner, &
                               isubgroup)
  END IF
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE linsol_initDataBiCGStab (rsolverNode, isolverSubgroup)
  
!<description>
  
  ! Calls the initData subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_initData.
  
!</description>
  
!<inputoutput>
  
  ! The t_linsolNode structure of the UMFPACK4 solver
  TYPE(t_linsolNode), INTENT(INOUT)         :: rsolverNode
   
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
  INTEGER :: isubgroup
  
  ! by default, initialise solver subroup 0
  isubgroup = 0
  IF (PRESENT(isolversubgroup)) isubgroup = isolverSubgroup

  ! We simply pass isubgroup to the subsolvers when calling them.
  ! Inside of this routine, there's not much to do with isubgroup,
  ! as there is no special information (like a factorisation)
  ! associated with this type of solver, which has to be allocated,
  ! released or updated...

  ! Call the init routine of the preconditioner.
  IF (ASSOCIATED(rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner)) THEN
    CALL linsol_initData (rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner, &
                          isubgroup)
  END IF
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE linsol_doneDataBiCGStab (rsolverNode, isolverSUbgroup)
  
!<description>
  
  ! Calls the doneData subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_doneData.
  
!</description>
  
!<inputoutput>
  
  ! The t_linsolNode structure of the UMFPACK4 solver
  TYPE(t_linsolNode), INTENT(INOUT)         :: rsolverNode
   
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
  INTEGER :: isubgroup
  
  ! by default, initialise solver subroup 0
  isubgroup = 0
  IF (PRESENT(isolversubgroup)) isubgroup = isolverSubgroup

  ! We simply pass isubgroup to the subsolvers when calling them.
  ! Inside of this routine, there's not much to do with isubgroup,
  ! as there is no special information (like a factorisation)
  ! associated with this type of solver, which has to be allocated,
  ! released or updated...

  ! Call the done routine of the preconditioner.
  IF (ASSOCIATED(rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner)) THEN
    CALL linsol_doneData (rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner, &
                          isubgroup)
  END IF
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE linsol_doneStructureBiCGStab (rsolverNode, isolverSubgroup)
  
!<description>
  
  ! Calls the doneStructure subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_doneStructure.
  
!</description>
  
!<inputoutput>
  
  ! The t_linsolNode structure of the UMFPACK4 solver
  TYPE(t_linsolNode), INTENT(INOUT)         :: rsolverNode
   
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
  INTEGER :: isubgroup,i
  TYPE(t_linsolSubnodeBiCGStab), POINTER :: p_rsubnode
  
  ! by default, initialise solver subroup 0
  isubgroup = 0
  IF (PRESENT(isolversubgroup)) isubgroup = isolverSubgroup

  ! We simply pass isubgroup to the subsolvers when calling them.
  ! Inside of this routine, there's not much to do with isubgroup,
  ! as there is no special information (like a factorisation)
  ! associated with this type of solver, which has to be allocated,
  ! released or updated...

  ! Call the init routine of the preconditioner.
  IF (ASSOCIATED(rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner)) THEN
    CALL linsol_doneStructure (rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner, &
                               isubgroup)
  END IF
  
  ! Release temporary data
  p_rsubnode => rsolverNode%p_rsubnodeBiCGStab
  DO i=6,1,-1
    CALL lsysbl_releaseVector (p_rsubnode%RtempVectors(i))
  END DO
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_doneBiCGStab (rsolverNode)
  
!<description>
  
  ! This routine releases all temporary memory for the BiCGStab solver from
  ! the heap. In particular, if a preconditioner is attached to the solver
  ! structure, it's also released from the heap by calling 
  ! linsol_releaseSolver for it.
  ! This DONE routine is declared as RECURSIVE to prevent a clean
  ! interaction with linsol_releaseSolver.
  
!</description>
  
!<input>
  
  ! A pointer to a t_linsolNode structure of the BiCGStab solver.
  TYPE(t_linsolNode), POINTER         :: rsolverNode
   
!</input>
  
!</subroutine>

  ! Check if there's a preconditioner attached. If yes, release it.
  IF (ASSOCIATED(rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner)) THEN
    CALL linsol_releaseSolver(rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner)
  END IF
  
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE linsol_solveBiCGStab (rsolverNode,rx,rb,rtemp)
  
!<description>
  
  ! Solves the linear system $Ax=b$ with BiCGStab. The matrix $A$ must be
  ! attached to the solver previously by linsol_setMatrices.
  ! linsol_initProblemStructure and linsol_initProblemData must have
  ! been called so that BiCGStab is prepared to solve the system. 
  !
  ! The implementation follows the original paper introducing BiCGStab:
  !   van der Vorst, H.A.; BiCGStab: A Fast and Smoothly Converging
  !   Variant of Bi-CG for the Solution of Nonsymmetric Linear Systems;
  !   SIAM J. Sci. Stat. Comput. 1992, Vol. 13, No. 2, pp. 631-644
  
!</description>
  
!<input>
  
  ! Right hand side of the system
  TYPE(t_vectorBlock), INTENT(IN)           :: rb
  
!</input>
  
!<inputoutput>
  
  ! The t_linsolNode structure of the UMFPACK4 solver
  TYPE(t_linsolNode), INTENT(INOUT)         :: rsolverNode
   
  ! Initial solution vector; receives the final solution vector when the
  ! iteration finishes.
  TYPE(t_vectorBlock), INTENT(INOUT)        :: rx
   
  ! A temporary vector of the same size and structure as rx.
  TYPE(t_vectorBlock), INTENT(INOUT)        :: rtemp

!</inputoutput>
  
!</subroutine>

  ! local variables
  REAL(DP) :: dalpha,dbeta,domega0,domega1,domega2,dres
  REAL(DP) :: drho1,drho0,dfr
  INTEGER :: ireslength,ite,i

  ! The queue saves the current residual and the two previous residuals.
  REAL(DP), DIMENSION(32) :: Dresqueue
  
  ! The system matrix
  TYPE(t_matrixBlock), POINTER :: p_rmatrix
  
  ! Minimum number of iterations, print-sequence for residuals
  INTEGER :: nminIterations, niteResOutput
  
  ! Whether to filter/prcondition
  LOGICAL bprec,bfilter
  
  ! Our structure
  TYPE(t_linsolSubnodeBiCGStab), POINTER :: p_rsubnode
  
  ! Pointers to temporary vectors - named for easier access
  TYPE(t_vectorBlock), POINTER :: p_DR,p_DR0,p_DP,p_DPA,p_DSA,p_rprec
  TYPE(t_linsolNode), POINTER :: p_rprecSubnode
  TYPE(t_filterChain), DIMENSION(:), POINTER :: p_RfilterChain
  
  ! Solve the system!
  
    ! Status reset
    rsolverNode%iresult = 0
    
    ! Getch some information
    p_rsubnode => rsolverNode%p_rsubnodeBiCGStab
    p_rmatrix => rsolverNode%rlinearSystemInfo%p_rsystemMatrix

    ! Check the parameters
    IF ((rx%NEQ .EQ. 0) .OR. (rb%NEQ .EQ. 0) .OR. (rtemp%NEQ .EQ. 0) .OR. &
        (rb%NEQ .NE. rx%NEQ) .OR. (rtemp%NEQ .NE. rx%NEQ) .OR. &
        (p_rmatrix%NEQ .EQ. 0) .OR. (p_rmatrix%NEQ .NE. rx%NEQ) ) THEN
    
      ! Parameters wrong
      rsolverNode%iresult = 2
      RETURN
    END IF

    ! Length of the queue of last residuals for the computation of
    ! the asymptotic convergence rate

    ireslength = MAX(0,MIN(32,rsolverNode%niteAsymptoticCVR))

    ! Minimum number of iterations
 
    nminIterations = MAX(rsolverNode%nminIterations,0)
      
    ! Use preconditioning? Filtering?

    bprec = ASSOCIATED(rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner)
    bfilter = ASSOCIATED(rsolverNode%p_rsubnodeBiCGStab%p_RfilterChain)
    
    ! Iteration when the residuum is printed:

    niteResOutput = MAX(1,rsolverNode%niteResOutput)

    ! Set pointers to the temporary vectors
    p_DR   => p_rsubnode%RtempVectors(1)
    p_DR0  => p_rsubnode%RtempVectors(2)
    p_DP   => p_rsubnode%RtempVectors(3)
    p_DPA  => p_rsubnode%RtempVectors(4)
    p_DSA  => p_rsubnode%RtempVectors(5)
    
    ! All vectors share the same boundary conditions as rx!
    ! So assign now all discretisation-related information (boundary
    ! conditions,...) to the temporary vectors.
    CALL lsysbl_assignDiscretIndirect (rx,p_DR )
    CALL lsysbl_assignDiscretIndirect (rx,p_DR0)
    CALL lsysbl_assignDiscretIndirect (rx,p_DP )
    CALL lsysbl_assignDiscretIndirect (rx,p_DPA)
    CALL lsysbl_assignDiscretIndirect (rx,p_DSA)
    
    IF (bprec) THEN
      p_rprec => p_rsubnode%RtempVectors(6)
      p_rprecSubnode => p_rsubnode%p_rpreconditioner
    END IF
    IF (bfilter) THEN
      p_RfilterChain => p_rsubnode%p_RfilterChain
    END IF
      
    ! Initialize used vectors with zero
      
    CALL lsysbl_vectorClear(p_DP)
    CALL lsysbl_vectorClear(p_DPA)

    ! Initialization

    drho0  = 1.0_DP
    dalpha = 1.0_DP
    domega0 = 1.0_DP

    ! Build the initial defect precond*(b-Ax) using
    ! Matrix-vector multiplication with preconditioning.

    CALL lsysbl_vectorCopy(rb,p_DR)
    CALL lsysbl_blockMatVec (p_rmatrix, rx,p_DR, -1.0_DP, 1.0_DP)
    IF (bfilter) THEN
      ! Apply the filter chain to the vector
      CALL filter_applyFilterChainVec (p_DR, p_RfilterChain)
    END IF
    IF (bprec) THEN
      ! Perform preconditioning with the assigned preconditioning
      ! solver structure.
      CALL lsysbl_vectorCopy(p_DR,p_rprec)
      CALL linsol_performSolve (p_rprecSubnode,p_DR,p_Rprec,rtemp)
    END IF
    
    ! Get the norm of the residuum
    dres = lsysbl_vectorNorm (p_DR,rsolverNode%iresNorm)
    IF (.NOT.((dres .GE. 1D-99) .AND. &
              (dres .LE. 1D99))) dres = 0.0_DP

    ! Initialize starting residuum
      
    rsolverNode%dinitialDefect = dres

    ! initialize the queue of the last residuals with RES

    Dresqueue = dres

    ! Check if out initial defect is zero. This may happen if the filtering
    ! routine filters "everything out"!
    ! In that case we can directly stop our computation.

    IF ( rsolverNode%dinitialDefect .LT. rsolverNode%drhsZero ) THEN
     
      ! final defect is 0, as initialised in the output variable above

      CALL lsysbl_vectorClear(rx)
      ite = 0
      rsolverNode%dfinalDefect = dres
          
    ELSE

      IF (rsolverNode%ioutputLevel .GE. 2) THEN
        PRINT *,&
          'BiCGStab: Iteration ',0,',  !!RES!! = ',rsolverNode%dinitialDefect
      END IF

      CALL lsysbl_vectorCopy(p_DR,p_DR0)

      ! Perform at most nmaxIterations loops to get a new vector

      DO ite = 1,rsolverNode%nmaxIterations
      
        rsolverNode%icurrentIteration = ite

        drho1 = lsysbl_scalarProduct (p_DR0,p_DR) 

        IF (drho0*domega0 .EQ. 0.0_DP) THEN
          ! Should not happen
          IF (rsolverNode%ioutputLevel .GE. 2) THEN
            PRINT *,&
      'BiCGStab: Iteration prematurely stopped! Correction vector is zero!'
          END IF

          ! Some tuning for the output, then cancel.

          rsolverNode%iresult = -1
          rsolverNode%iiterations = ITE-1
          EXIT
          
        END IF

        dbeta=(drho1*dalpha)/(drho0*domega0)
        drho0 = drho1

        CALL lsysbl_vectorLinearComb (p_DR ,p_DP,1.0_DP,dbeta)
        CALL lsysbl_vectorLinearComb (p_DPA ,p_DP,-dbeta*domega0,1.0_DP)

        CALL lsysbl_blockMatVec (p_rmatrix, p_DP,p_DPA, 1.0_DP,0.0_DP)
        IF (bfilter) THEN
          ! Apply the filter chain to the vector
          CALL filter_applyFilterChainVec (p_DPA, p_RfilterChain)
        END IF
        IF (bprec) THEN
          ! Perform preconditioning with the assigned preconditioning
          ! solver structure.
          CALL lsysbl_vectorCopy(p_DPA,p_rprec)
          CALL linsol_performSolve (p_rprecSubnode,p_DPA,p_Rprec,rtemp)
        END IF

        dalpha = lsysbl_scalarProduct (p_DR0,p_DPA)
        
        IF (dalpha .EQ. 0.0_DP) THEN
          ! We are below machine exactness - we can't do anything more...
          ! May happen with very small problems with very few unknowns!
          IF (rsolverNode%ioutputLevel .GE. 2) THEN
            PRINT *,'BiCGStab: Convergence failed!'
            rsolverNode%iresult = 1
            EXIT
          END IF
        END IF
        
        dalpha = drho1/dalpha

        CALL lsysbl_vectorLinearComb (p_DPA,p_DR,-dalpha,1.0_DP)

        CALL lsysbl_blockMatVec (p_rmatrix, p_DR,p_DSA, 1.0_DP,0.0_DP)
        
        IF (bfilter) THEN
          ! Apply the filter chain to the vector
          CALL filter_applyFilterChainVec (p_DSA, p_RfilterChain)
        END IF
        IF (bprec) THEN
          ! Perform preconditioning with the assigned preconditioning
          ! solver structure.
          CALL lsysbl_vectorCopy(p_DSA,p_rprec)
          CALL linsol_performSolve (p_rprecSubnode,p_DSA,p_Rprec,rtemp)
        END IF
        
        domega1 = lsysbl_scalarProduct (p_DSA,p_DR)
        domega2 = lsysbl_scalarProduct (p_DSA,p_DSA)
        
        IF (domega1 .EQ. 0.0_DP) THEN
          domega0 = 0.0_DP
        ELSE
          IF (domega2 .EQ. 0.0_DP) THEN
            IF (rsolverNode%ioutputLevel .GE. 2) THEN
              PRINT *,'BiCGStab: Convergence failed!'
              rsolverNode%iresult = 1
              EXIT
            END IF
          END IF
          domega0 = domega1/domega2
        END IF

        CALL lsysbl_vectorLinearComb (p_DP ,rx,dalpha,1.0_DP)
        CALL lsysbl_vectorLinearComb (p_DR ,rx,domega0,1.0_DP)

        CALL lsysbl_vectorLinearComb (p_DSA,p_DR,-domega0,1.0_DP)

        ! Get the norm of the new (final?) residuum
        dfr = lsysbl_vectorNorm (p_DR,rsolverNode%iresNorm)
     
        ! Shift the queue with the last residuals and add the new
        ! residual to it
        dresqueue(1:ireslength) = EOSHIFT(dresqueue(1:ireslength),1,dfr)

        rsolverNode%dfinalDefect = dfr
     
        ! At least perform nminIterations iterations

        IF (ite .GE. nminIterations) THEN
        
          ! Check if the iteration converged
          IF (linsol_testConvergence(rsolverNode,dfr)) EXIT
          
        END IF

        ! print out the current residuum

        IF ((rsolverNode%ioutputLevel .GE. 2) .AND. &
            (MOD(ite,niteResOutput).EQ.0)) THEN
          !WRITE (MTERM,'(A,I7,A,D25.16)') 
          PRINT *,'BiCGStab: Iteration ',ITE,',  !!RES!! = ',rsolverNode%dfinalDefect
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
        !WRITE (MTERM,'(A,I7,A,D25.16)') 
        PRINT *,'BiCGStab: Iteration ',ITE,',  !!RES!! = ',rsolverNode%dfinalDefect
      END IF

    END IF

    rsolverNode%iiterations = ite
      
    ! Don't calculate anything if the final residuum is out of bounds -
    ! would result in NaN's,...
      
      IF (rsolverNode%dfinalDefect .LT. 1E99_DP) THEN
      
        ! Calculate asymptotic convergence rate
      
        IF (dresqueue(1) .GE. 1E-70_DP) THEN
          I = MAX(1,MIN(rsolverNode%iiterations,ireslength-1))
          rsolverNode%dasymptoticConvergenceRate = &
            (rsolverNode%dfinalDefect / dresqueue(1))**(1.0_DP/REAL(I,DP))
        END IF

        ! If the initial defect was zero, the solver immediately
        ! exits - and so the final residuum is zero and we performed
        ! no steps; so the resulting multigrid convergence rate stays zero.
        ! In the other case the multigrid convergence rate computes as
        ! (final defect/initial defect) ** 1/nit :

        IF (rsolverNode%dfinalDefect .GT. rsolverNode%drhsZero) THEN
          rsolverNode%dconvergenceRate = &
                      (rsolverNode%dfinalDefect / rsolverNode%dinitialDefect) ** &
                      (1.0_DP/REAL(rsolverNode%iiterations,DP))
        END IF
        
        IF (rsolverNode%ioutputLevel .GE. 2) THEN
          !WRITE (MTERM,'(A)') ''
          PRINT *
          !WRITE (MTERM,'(A)') 
          PRINT *,'BiCGStab statistics:'
          !WRITE (MTERM,'(A)') ''
          PRINT *
          !WRITE (MTERM,'(A,I5)')     'Iterations              : ',
     !*            IPARAM(OITE)
          PRINT *,'Iterations              : ',rsolverNode%iiterations
          !WRITE (MTERM,'(A,D24.12)') '!!INITIAL RES!!         : ',
     !*            DPARAM(ODEFINI)
          PRINT *,'!!INITIAL RES!!         : ',rsolverNode%dinitialDefect
          !WRITE (MTERM,'(A,D24.12)') '!!RES!!                 : ',
     !*            DPARAM(ODEFFIN)
          PRINT *,'!!RES!!                 : ',rsolverNode%dfinalDefect
          IF (rsolverNode%dinitialDefect .GT. rsolverNode%drhsZero) THEN     
            !WRITE (MTERM,'(A,D24.12)') '!!RES!!/!!INITIAL RES!! : ',
     !*              rparam%dfinalDefect / rparam%dinitialDefect
            PRINT *,'!!RES!!/!!INITIAL RES!! : ',&
                    rsolverNode%dfinalDefect / rsolverNode%dinitialDefect
          ELSE
            !WRITE (MTERM,'(A,D24.12)') '!!RES!!/!!INITIAL RES!! : ',
     !*              0D0
            PRINT*,'!!RES!!/!!INITIAL RES!! : ',0.0_DP
          END IF
          !WRITE (MTERM,'(A)') ''
          PRINT *
          !WRITE (MTERM,'(A,D24.12)') 'Rate of convergence     : ',
     !*            DPARAM(ORHO)
          PRINT *,'Rate of convergence     : ',rsolverNode%dconvergenceRate

        END IF

        IF (rsolverNode%ioutputLevel .EQ. 1) THEN
!          WRITE (MTERM,'(A,I5,A,D24.12)') 
!     *          'BiCGStab: Iterations/Rate of convergence: ',
!     *          IPARAM(OITE),' /',DPARAM(ORHO)
!          WRITE (MTERM,'(A,I5,A,D24.12)') 
!     *          'BiCGStab: Iterations/Rate of convergence: ',
!     *          IPARAM(OITE),' /',DPARAM(ORHO)
          PRINT *,&
                'BiCGStab: Iterations/Rate of convergence: ',&
                rsolverNode%iiterations,' /',rsolverNode%dconvergenceRate
        END IF
        
      ELSE
        ! DEF=Infinity; RHO=Infinity, set to 1
        rsolverNode%dconvergenceRate = 1.0_DP
        rsolverNode%dasymptoticConvergenceRate = 1.0_DP
      END IF  
  
  END SUBROUTINE
  
! *****************************************************************************
! Routines for the Multigrid solver
! *****************************************************************************

!<subroutine>
  
  SUBROUTINE linsol_addMultigridLevel (p_rlevelInfo,rsolverNode, &
                    p_rpresmoother,p_rpostsmoother,p_rcoarseGridSolver,iappend)
                    
!<description>
  
  ! This routine adds a new level to the linked list of levels in the multigrid
  ! solver. The given coarse-grid solver and smoother-structures are attached
  ! to the level. The routine returns a pointer to the new level-info structure
  ! in p_rlevelInfo to allow the caller to modify the standard settings of
  ! that level if necessary.
  
!</description>
  
!<input>
  
  ! Optional: A pointer to the solver structure of a solver that should be 
  ! used for presmoothing. This structure is used as a template to create an
  ! appropriate solver structure for all the levels. The structure itself is
  ! used at the finest level.
  ! If not given or set to NULL(), no presmoother will be used.
  TYPE(t_linsolNode), POINTER, OPTIONAL   :: p_rpresmoother
  
  ! Optional: A pointer to the solver structure of a solver that should be 
  ! used for postsmoothing. This structure is used as a template to create an
  ! appropriate solver structure for all the levels. The structure itself is
  ! used at the finest level.
  ! If not given or set to NULL(), no presmoother will be used.
  TYPE(t_linsolNode), POINTER, OPTIONAL   :: p_rpostsmoother

  ! Optional: A pointer to the solver structure of a solver that should be 
  ! used for coarse grid solving. 
  ! Should only be given for the very first level.
  TYPE(t_linsolNode), POINTER, OPTIONAL   :: p_rcoarseGridSolver
  
  ! Optional: Position where to put the new structure.
  ! YES or not specified: append the structure to the end of the list as new
  ! higher level.
  ! NO: Insert the structure as new lowest level. The caller has to make sure
  ! that the coarse grid solver on the previous lowest level is removed!
  INTEGER, INTENT(IN), OPTIONAL                  :: iappend
  
!</input>
  
!<inputoutput>
  
  ! The solver structure of the multigrid solver, where the level
  ! should be added to. Must already be initialised for the multigrid
  ! solver.
  TYPE(t_linsolNode), POINTER, OPTIONAL :: rsolverNode
  
!</inputoutput>
  
!<output>
  
!</subroutine>

  ! local variables
  INTEGER :: i
  
  ! The t_levelInfo structure for the new level that was added to the
  ! multigrid solver.
  TYPE(t_linsolMGLevelInfo), POINTER     :: p_rlevelInfo
  
!</output>
  
  ! Make sure the solver node is configured for multigrid
  IF ((rsolverNode%calgorithm .NE. LINSOL_ALG_MULTIGRID) .OR. &
      (.NOT. ASSOCIATED(rsolverNode%p_rsubnodeMultigrid))) THEN
    PRINT *,'Error: Multigrid structure not initialised'
    STOP
  END IF
  
  ! Create a new level-info structure
  NULLIFY(p_rlevelInfo)
  ALLOCATE(p_rlevelInfo)
  
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

  ! Attach the level-info structure to the linked list.
  ! Does the list exist?
  IF (.NOT. ASSOCIATED(rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoHead)) THEN
  
    ! List does not exist. Create a new one.
    rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoHead => p_rlevelInfo
    rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoTail => p_rlevelInfo
    rsolverNode%p_rsubnodeMultigrid%nlevels = 1
  
  ELSE
  
    ! Do we have to insert it as new head?
    i = YES
    IF (PRESENT(iappend)) i=iappend
    
    IF (i .EQ. NO) THEN
    
      ! New lowest level
      p_rlevelInfo%p_rnextLevel => rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoHead
      rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoHead => p_rlevelInfo
    
    ELSE
    
      ! New highest level
      p_rlevelInfo%p_rprevLevel => rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoTail
      rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoTail => p_rlevelInfo
      
    END IF
    
    ! Increase the number of existing levels
    rsolverNode%p_rsubnodeMultigrid%nlevels = &
      rsolverNode%p_rsubnodeMultigrid%nlevels + 1
  
  END IF

  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE linsol_initMultigrid (p_rsolverNode,p_Rfilter)
  
!<description>
  
  ! Creates a t_linsolNode solver structure for the Multigrid solver. The node
  ! can be used to directly solve a problem or to be attached as solver
  ! or preconditioner to another solver structure. The node can be deleted
  ! by linsol_releaseSolver.
  ! Before the solver can used for solving the problem, the caller must add
  ! information about all levels (matrices,...) to the solver. 
  ! This can be done by linsol_addMultigridLevel.
  
!</description>
  
!<input>
  
  ! Optional: A pointer to a filter chain (i.e. an array of t_filterChain
  ! structures) if filtering should be applied to the vector during the 
  ! iteration. If not given or set to NULL(), no filtering will be used.
  ! The filter chain (i.e. the array) must exist until the system is solved!
  ! The filter chain must be configured for being applied to defect vectors.
  TYPE(t_filterChain), DIMENSION(:), POINTER, OPTIONAL   :: p_Rfilter
  
!</input>
  
!<output>
  
  ! A pointer to a t_linsolNode structure. Is set by the routine, any previous
  ! value of the pointer is destroyed.
  TYPE(t_linsolNode), POINTER             :: p_rsolverNode
   
!</output>
  
!</subroutine>

  ! Create a default solver structure
  CALL linsol_initSolverGeneral(p_rsolverNode)
  
  ! Initialise the type of the solver
  p_rsolverNode%calgorithm = LINSOL_ALG_MULTIGRID
  
  ! Initialise the ability bitfield with the ability of this solver:
  p_rsolverNode%ccapability = LINSOL_ABIL_SCALAR     + LINSOL_ABIL_BLOCK        + &
                              LINSOL_ABIL_MULTILEVEL + LINSOL_ABIL_CHECKDEF     + &
                              LINSOL_ABIL_USESUBSOLVER + &
                              LINSOL_ABIL_USEFILTER
  
  ! Allocate the subnode for Multigrid.
  ! This initialises most of the variables with default values appropriate
  ! to this solver.
  ALLOCATE(p_rsolverNode%p_rsubnodeMultigrid)
  
  ! Attach the filter if given. 
  
  IF (PRESENT(p_Rfilter)) THEN
    p_rsolverNode%p_rsubnodeMultigrid%p_RfilterChain => p_Rfilter
  END IF
  
  END SUBROUTINE
  
! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_setMatrixMultigrid (rsolverNode,Rmatrices)
  
!<description>
  
  ! Assigns the system matrix rsystemMatrix to the Multigrid solver on 
  ! all levels.
  
!</description>
  
!<input>
  
  ! An array of system matrices on all levels.
  ! Each level in multigrid is initialised separately.
  TYPE(t_LinearSystemInfo), INTENT(IN), DIMENSION(:)   :: Rmatrices

!</input>
  
!<inputoutput>
  
  ! The t_linsolNode structure of the UMFPACK4 solver
  TYPE(t_linsolNode), INTENT(INOUT)         :: rsolverNode
   
!</inputoutput>
  
!</subroutine>

  ! local variables
  TYPE(t_linsolMGLevelInfo), POINTER :: p_rcurrentLevel
  INTEGER                       :: ilevel,nlmin
  
  ! Make sure the solver node is configured for multigrid
  IF ((rsolverNode%calgorithm .NE. LINSOL_ALG_MULTIGRID) .OR. &
      (.NOT. ASSOCIATED(rsolverNode%p_rsubnodeMultigrid))) THEN
    PRINT *,'Error: Multigrid structure not initialised'
    STOP
  END IF
  
  ! Make sure we have the right amount of matrices
  IF (SIZE(Rmatrices) .NE. rsolverNode%p_rsubnodeMultigrid%nlevels) THEN
    PRINT *,'Error: Wrong number of matrices'
    STOP
  END IF

  ! Check for every level, if there's a presmoother, postsmoother or
  ! coarse grid solver structure attached. If yes, attach the matrix
  ! to it.
  ! For this purpose, call the general initialisation routine, but
  ! pass only that part of the Rmatrices array that belongs
  ! to the range of levels between the coarse grid and the current grid.
  
  p_rcurrentLevel => rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoHead
  nlmin  = LBOUND(Rmatrices,1)
  ilevel = LBOUND(Rmatrices,1)
  DO WHILE(ASSOCIATED(p_rcurrentLevel))
    IF (ASSOCIATED(p_rcurrentLevel%p_rpreSmoother)) THEN
      CALL linsol_setMatrices (p_rcurrentLevel%p_rpreSmoother, &
                                Rmatrices(nlmin:ilevel))
    END IF
    IF (ASSOCIATED(p_rcurrentLevel%p_rpostSmoother)) THEN
      CALL linsol_setMatrices (p_rcurrentLevel%p_rpostSmoother, &
                                Rmatrices(nlmin:ilevel))
    END IF
    IF (ASSOCIATED(p_rcurrentLevel%p_rcoarseGridSolver)) THEN
      CALL linsol_setMatrices (p_rcurrentLevel%p_rcoarseGridSolver, &
                                Rmatrices(nlmin:ilevel))
    END IF
    p_rcurrentLevel => p_rcurrentLevel%p_rnextLevel
    ilevel = ilevel + 1
  END DO
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_initStructureMultigrid (rsolverNode,isolverSubgroup)
  
!<description>
  
  ! Calls the initStructure subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_initStructure.
  
!</description>
  
!<inputoutput>
  
  ! The t_linsolNode structure of the UMFPACK4 solver
  TYPE(t_linsolNode), INTENT(INOUT)         :: rsolverNode
   
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
  TYPE(t_linsolMGLevelInfo), POINTER :: p_rcurrentLevel
  INTEGER :: isubgroup
  
  ! by default, initialise solver subroup 0
  isubgroup = 0
  IF (PRESENT(isolversubgroup)) isubgroup = isolverSubgroup

  ! We simply pass isubgroup to the subsolvers when calling them.
  ! Inside of this routine, there's not much to do with isubgroup,
  ! as there is no special information (like a factorisation)
  ! associated with this type of solver, which has to be allocated,
  ! released or updated...

  ! Call the init routine of the preconditioner.

  ! Make sure the solver node is configured for multigrid
  IF ((rsolverNode%calgorithm .NE. LINSOL_ALG_MULTIGRID) .OR. &
      (.NOT. ASSOCIATED(rsolverNode%p_rsubnodeMultigrid))) THEN
    PRINT *,'Error: Multigrid structure not initialised'
    STOP
  END IF

  ! Check for every level, if there's a presmoother, postsmoother or
  ! coarse grid solver structure attached. 
  
  p_rcurrentLevel => rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoHead
  DO WHILE(ASSOCIATED(p_rcurrentLevel))
    IF (ASSOCIATED(p_rcurrentLevel%p_rpreSmoother)) THEN
      CALL linsol_initStructure(p_rcurrentLevel%p_rpreSmoother,isubgroup)
    END IF
    IF (ASSOCIATED(p_rcurrentLevel%p_rpostSmoother)) THEN
      CALL linsol_initStructure(p_rcurrentLevel%p_rpostSmoother,isubgroup)
    END IF
    IF (ASSOCIATED(p_rcurrentLevel%p_rcoarseGridSolver)) THEN
      CALL linsol_initStructure(p_rcurrentLevel%p_rcoarseGridSolver,isubgroup)
    END IF
    p_rcurrentLevel => p_rcurrentLevel%p_rnextLevel
  END DO
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE linsol_initDataMultigrid (rsolverNode,isolverSUbgroup)
  
!<description>
  
  ! Calls the initData subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_initData.
  
!</description>
  
!<inputoutput>
  
  ! The t_linsolNode structure of the UMFPACK4 solver
  TYPE(t_linsolNode), INTENT(INOUT)         :: rsolverNode
   
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
  TYPE(t_linsolMGLevelInfo), POINTER :: p_rcurrentLevel
  INTEGER :: isubgroup
  
  ! by default, initialise solver subroup 0
  isubgroup = 0
  IF (PRESENT(isolversubgroup)) isubgroup = isolverSubgroup
  
  ! We simply pass isubgroup to the subsolvers when calling them.
  ! Inside of this routine, there's not much to do with isubgroup,
  ! as there is no special information (like a factorisation)
  ! associated with this type of solver, which has to be allocated,
  ! released or updated...

  ! Call the init routine of the preconditioner.
  
  ! Make sure the solver node is configured for multigrid
  IF ((rsolverNode%calgorithm .NE. LINSOL_ALG_MULTIGRID) .OR. &
      (.NOT. ASSOCIATED(rsolverNode%p_rsubnodeMultigrid))) THEN
    PRINT *,'Error: Multigrid structure not initialised'
    STOP
  END IF

  ! Check for every level, if there's a presmoother, postsmoother or
  ! coarse grid solver structure attached. 
  
  p_rcurrentLevel => rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoHead
  DO WHILE(ASSOCIATED(p_rcurrentLevel))
    IF (ASSOCIATED(p_rcurrentLevel%p_rpreSmoother)) THEN
      CALL linsol_initData(p_rcurrentLevel%p_rpreSmoother,isubgroup)
    END IF
    IF (ASSOCIATED(p_rcurrentLevel%p_rpostSmoother)) THEN
      CALL linsol_initData(p_rcurrentLevel%p_rpostSmoother,isubgroup)
    END IF
    IF (ASSOCIATED(p_rcurrentLevel%p_rcoarseGridSolver)) THEN
      CALL linsol_initData(p_rcurrentLevel%p_rcoarseGridSolver,isubgroup)
    END IF
    p_rcurrentLevel => p_rcurrentLevel%p_rnextLevel
  END DO
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE linsol_doneDataMultigrid (rsolverNode,isolverSubgroup)
  
!<description>
  
  ! Calls the doneData subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_doneData.
  
!</description>
  
!<inputoutput>
  
  ! The t_linsolNode structure of the UMFPACK4 solver
  TYPE(t_linsolNode), INTENT(INOUT)         :: rsolverNode
   
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
  TYPE(t_linsolMGLevelInfo), POINTER :: p_rcurrentLevel
  INTEGER :: isubgroup
  
  ! by default, initialise solver subroup 0
  isubgroup = 0
  IF (PRESENT(isolversubgroup)) isubgroup = isolverSubgroup

  ! We simply pass isubgroup to the subsolvers when calling them.
  ! Inside of this routine, there's not much to do with isubgroup,
  ! as there is no special information (like a factorisation)
  ! associated with this type of solver, which has to be allocated,
  ! released or updated...

  ! Call the init routine of the preconditioner.
  
  ! Make sure the solver node is configured for multigrid
  IF ((rsolverNode%calgorithm .NE. LINSOL_ALG_MULTIGRID) .OR. &
      (.NOT. ASSOCIATED(rsolverNode%p_rsubnodeMultigrid))) THEN
    PRINT *,'Error: Multigrid structure not initialised'
    STOP
  END IF

  ! Check for every level, if there's a presmoother, postsmoother or
  ! coarse grid solver structure attached. 
  
  p_rcurrentLevel => rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoHead
  DO WHILE(ASSOCIATED(p_rcurrentLevel))
    IF (ASSOCIATED(p_rcurrentLevel%p_rpreSmoother)) THEN
      CALL linsol_doneData(p_rcurrentLevel%p_rpreSmoother,isubgroup)
    END IF
    IF (ASSOCIATED(p_rcurrentLevel%p_rpostSmoother)) THEN
      CALL linsol_doneData(p_rcurrentLevel%p_rpostSmoother,isubgroup)
    END IF
    IF (ASSOCIATED(p_rcurrentLevel%p_rcoarseGridSolver)) THEN
      CALL linsol_doneData(p_rcurrentLevel%p_rcoarseGridSolver,isubgroup)
    END IF
    p_rcurrentLevel => p_rcurrentLevel%p_rnextLevel
  END DO
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE linsol_doneStructureMultigrid (rsolverNode,isolverSubgroup)
  
!<description>
  
  ! Calls the doneStructure subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_doneStructure.
  
!</description>
  
!<inputoutput>
  
  ! The t_linsolNode structure of the UMFPACK4 solver
  TYPE(t_linsolNode), INTENT(INOUT)         :: rsolverNode
   
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
  TYPE(t_linsolMGLevelInfo), POINTER :: p_rcurrentLevel
  INTEGER :: isubgroup
  
  ! by default, initialise solver subroup 0
  isubgroup = 0
  IF (PRESENT(isolversubgroup)) isubgroup = isolverSubgroup

  ! We simply pass isubgroup to the subsolvers when calling them.
  ! Inside of this routine, there's not much to do with isubgroup,
  ! as there is no special information (like a factorisation)
  ! associated with this type of solver, which has to be allocated,
  ! released or updated...

  ! Call the init routine of the preconditioner.
  
  ! Make sure the solver node is configured for multigrid
  IF ((rsolverNode%calgorithm .NE. LINSOL_ALG_MULTIGRID) .OR. &
      (.NOT. ASSOCIATED(rsolverNode%p_rsubnodeMultigrid))) THEN
    PRINT *,'Error: Multigrid structure not initialised'
    STOP
  END IF

  ! Check for every level, if there's a presmoother, postsmoother or
  ! coarse grid solver structure attached. 
  
  p_rcurrentLevel => rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoHead
  DO WHILE(ASSOCIATED(p_rcurrentLevel))
    IF (ASSOCIATED(p_rcurrentLevel%p_rpreSmoother)) THEN
      CALL linsol_doneStructure(p_rcurrentLevel%p_rpreSmoother,isubgroup)
    END IF
    IF (ASSOCIATED(p_rcurrentLevel%p_rpostSmoother)) THEN
      CALL linsol_doneStructure(p_rcurrentLevel%p_rpostSmoother,isubgroup)
    END IF
    IF (ASSOCIATED(p_rcurrentLevel%p_rcoarseGridSolver)) THEN
      CALL linsol_doneStructure(p_rcurrentLevel%p_rcoarseGridSolver,isubgroup)
    END IF
    p_rcurrentLevel => p_rcurrentLevel%p_rnextLevel
  END DO
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE linsol_doneMultigrid (rsolverNode)
  
!<description>
  
  ! This routine releases all temporary memory for the multigrid solver from
  ! the heap. In particular, this releases the solver structures of all
  ! subsolvers (smoother, coarse grid solver).
  ! This DONE routine is declared as RECURSIVE to prevent a clean
  ! interaction with linsol_releaseSolver.
  
!</description>
  
!<input>
  
  ! A pointer to a t_linsolNode structure of the Multigrid solver.
  TYPE(t_linsolNode), INTENT(OUT)                  :: rsolverNode
   
!</input>
  
!</subroutine>

  ! local variables
  TYPE(t_linsolMGLevelInfo), POINTER :: p_rcurrentLevel
  
  ! Make sure the solver node is configured for multigrid
  IF ((rsolverNode%calgorithm .NE. LINSOL_ALG_MULTIGRID) .OR. &
      (.NOT. ASSOCIATED(rsolverNode%p_rsubnodeMultigrid))) THEN
    PRINT *,'Error: Multigrid structure not initialised'
    STOP
  END IF

  ! Check for every level, if there's a presmoother, postsmoother or
  ! coarse grid solver structure attached. If yes, release it.
  
  p_rcurrentLevel => rsolverNode%p_rsubnodeMultigrid%p_rlevelInfoHead
  DO WHILE(ASSOCIATED(p_rcurrentLevel))
    IF (ASSOCIATED(p_rcurrentLevel%p_rpreSmoother)) THEN
      CALL linsol_releaseSolver(p_rcurrentLevel%p_rpreSmoother)
    END IF
    IF (ASSOCIATED(p_rcurrentLevel%p_rpostSmoother)) THEN
      CALL linsol_releaseSolver(p_rcurrentLevel%p_rpostSmoother)
    END IF
    IF (ASSOCIATED(p_rcurrentLevel%p_rcoarseGridSolver)) THEN
      CALL linsol_releaseSolver(p_rcurrentLevel%p_rcoarseGridSolver)
    END IF
    p_rcurrentLevel => p_rcurrentLevel%p_rnextLevel
  END DO
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE linsol_solveMultigrid (rsolverNode,rx,rb,rtemp)
  
!<description>
  
  ! Solves the linear system $Ax=b$ with Multigrid. The matrices $A$ on all
  ! levels must be attached to the solver previously by linsol_setMatrices.
  ! rb is the RHS and rx is the initial iteration vector.
  
!</description>
  
!<input>
  
  ! Right hand side of the system
  TYPE(t_vectorBlock), INTENT(IN),TARGET    :: rb
  
!</input>
  
!<inputoutput>
  
  ! The t_linsolNode structure of the UMFPACK4 solver
  TYPE(t_linsolNode), INTENT(INOUT)         :: rsolverNode
   
  ! Initial solution vector; receives the final solution vector when the
  ! iteration finishes.
  TYPE(t_vectorBlock), INTENT(INOUT)        :: rx
  
  ! A temporary vector of the same size and structure as rx.
  TYPE(t_vectorBlock), INTENT(INOUT)        :: rtemp

!</inputoutput>
  
!</subroutine>

  ! solve the system
  PRINT *,'to be implemented...'
  
  END SUBROUTINE
  
END MODULE
