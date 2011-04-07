!##############################################################################
!# ****************************************************************************
!# <name> linearspacetimesolver </name>
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

module linearspacetimesolver

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
  use spdiscprojection
  use nonlinearsolver
  use paramlist
  use linearsolverautoinitialise
  use matrixrestriction
  use paramlist
  use timestepping
  use globalsystem
  
  use collection
  use convection
    
  use basicstructures
  use user_callback

  use spacepreconditioner
  use spacepreconditionerinit
  use stationaryoptcsolver
  use timeanalysis
  use spatialbc
  use spacediscretisation
  use postprocessing
  use spacematvecassembly
  use spacepreconditioner
  use spacetimediscretisation
  use spacetimelinearsystem
  use statistics
  
  use spacetimevectors
  use spacetimeinterlevelprojection
  use dofmapping
  use timeboundaryconditions
  use spacetimevanca
  
  use matrixio
    
  implicit none

!<constantblock description="Algorithm identifiers">

  ! Undefined algorithm
  integer, parameter :: SPTILS_ALG_UNDEFINED     = 0
  
  ! Preconditioned defect correction (Richardson iteration);
  ! $x_{n+1} = x_n + \omega P^{-1} (b-Ax)$
  integer, parameter :: SPTILS_ALG_DEFCORR       = 1
  
  ! Block Jacobi preconditioner which applies a spatial preconditioner
  ! to each diagonal block of the coupled space-time system.
  integer, parameter :: SPTILS_ALG_BLOCKJACOBI   = 2
  
  ! Block Gauss-Seidel preconditioner
  integer, parameter :: SPTILS_ALG_BlockFBSOR      = 4
  
  ! CG iteration (preconditioned) 
  integer, parameter :: SPTILS_ALG_CG            = 6

  ! BiCGStab iteration (preconditioned) 
  integer, parameter :: SPTILS_ALG_BICGSTAB      = 7

  ! UMFPACK solver
  integer, parameter :: SPTILS_ALG_UMFPACK4      = 11

  ! Block orward-Backward Gauss-Seidel solver
  integer, parameter :: SPTILS_ALG_BLOCKFBGS    = 12

  ! Multigrid iteration
  integer, parameter :: SPTILS_ALG_MULTIGRID     = 9

!</constantblock>

! *****************************************************************************

!<constantblock description="Error constants returned by initialisation routines">

  ! Initialisation routine went fine
  integer, parameter :: SPTILS_ERR_NOERROR       = 0
  
!</constantblock>

! *****************************************************************************

!<constantblock description="Bitfield identifiers for the ability of a solver">

  ! Solver can handle multiple levels
  integer(I32), parameter :: SPTILS_ABIL_MULTILEVEL   = 2**2
  
  ! Solver allows checking the defect during the iteration.
  ! Solvers not capable of this perform only a fixed number of solution
  ! steps (e.g. UMFPACK performs always one step).
  integer(I32), parameter :: SPTILS_ABIL_CHECKDEF     = 2**3
  
  ! Solver is a direct solver (e.g. UMFPACK, ILU).
  ! Otherwise the solver is of iterative nature and might perform
  ! multiple steps to solve the problem.
  integer(I32), parameter :: SPTILS_ABIL_DIRECT       = 2**4
  
  ! Solver might use subsolvers (preconditioners, smoothers,...)
  integer(I32), parameter :: SPTILS_ABIL_USESUBSOLVER = 2**5
  
!</constantblock>

! *****************************************************************************

!<constantblock description="Identifiers for stopping criterium istoppingCriterion.">

  ! Use standard stopping criterion.
  ! If depsRel>0: use relative stopping criterion.
  ! If depsAbs>0: use abs stopping criterion.
  ! If both are > 0: use both, i.e. the iteration stops when both,
  !    the relative AND the absolute stopping criterium holds
  integer, parameter :: SPTILS_STOP_STANDARD     = 0

  ! Use 'minimum' stopping criterion.
  ! If depsRel>0: use relative stopping criterion.
  ! If depsAbs>0: use abs stopping criterion.
  ! If both are > 0: use one of them, i.e. the iteration stops when the
  !    either the relative OR the absolute stopping criterium holds
  integer, parameter :: SPTILS_STOP_ONEOF        = 1
  
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
  
  type t_sptilsNode
    
    ! OUTPUT: Result
    ! The result of the solution process.
    ! =0: success. =1: iteration broke down, diverging, =2: error in the parameters,
    ! <0: algorithm-specific error
    integer                    :: iresult
    
    ! OUTPUT: Number of performed iterations, if the solver
    ! is of iterative nature.
    ! Is to 1 by the solver if not used (indicating at least 1 performed 
    ! iteration, which is always the case).
    integer                    :: iiterations
    
    ! OUTPUT PARAMETER FOR SOLVERS WITH RESIDUAL CHECK: 
    ! Norm of initial residuum
    real(DP)                        :: dinitialDefect

    ! OUTPUT PARAMETER FOR SOLVERS WITH RESIDUAL CHECK: 
    ! Norm of final residuum
    real(DP)                        :: dfinalDefect

    ! OUTPUT PARAMETER FOR ITERATIVE SOLVERS WITH RESIDUAL CHECK: 
    ! Convergence rate
    real(DP)                        :: dconvergenceRate

    ! OUTPUT PARAMETER:
    ! Total time for solver
    real(DP)                        :: dtimeTotal

    ! OUTPUT PARAMETER FOR SOLVERS THAT SUPPORT FILTERING:
    ! Total time for filtering
    real(DP)                        :: dtimeFiltering
    
    ! INPUT PARAMETER:
    ! General solver parameter; solver specific use.
    ! Standard value = 1.0 (corresponds to 'no damping' e.g. with the defect
    ! correction iteration)
    real(DP)                        :: domega  = 1.0_DP

    ! INPUT PARAMETER FOR ITERATIVE SOLVERS: 
    ! Relative stopping criterion. Stop iteration if
    ! !!defect!! < EPSREL * !!initial defect!!.
    ! =0: ignore, use absolute stopping criterion; standard = 1E-5
    ! Remark: don't set depsAbs=depsRel=0!
    real(DP)                        :: depsRel = 1E-5_DP

    ! INPUT PARAMETER FOR ITERATIVE SOLVERS: 
    ! Absolute stopping criterion. Stop iteration if
    ! !!defect!! < EPSREL.
    ! =0: ignore, use relative stopping criterion; standard = 1E-5
    ! Remark: don't set depsAbs=depsRel=0!
    real(DP)                        :: depsAbs = 1E-5_DP
    
    ! INPUT PARAMETER FOR ITERATIVE SOLVERS:
    ! Difference in the residual stopping criterion. Stop iteration if
    ! !!defect_new!! - !!defect_old!! < depsDiff * !!defect_old!!.
    ! This stopping criterion holds additionally to depsAbs/depsRel
    ! in an OR sense -- i.e. the iteration is stopped if
    ! !!defect_new!! - !!defect_old!! < depsDiff * !!defect_old!!
    ! holds OR if the stopping criterion given by depsAbs/depsRel
    ! holds!
    real(DP)                        :: depsDiff = 0.0_DP

    ! INPUT PARAMETER FOR ITERATIVE SOLVERS: 
    ! Relative divergence criterion.  Treat iteration as
    ! diverged if
    !   !!defect!! >= DIVREL * !!initial defect!!
    ! A value of SYS_INFINITY_DP disables the relative divergence check.
    ! standard = 1E3
    real(DP)                        :: ddivRel = 1E3_DP

    ! INPUT PARAMETER FOR ITERATIVE SOLVERS: 
    ! Absolute divergence criterion.  Treat iteration as
    ! diverged if
    !   !!defect!! >= DIVREL
    ! A value of SYS_INFINITY_DP disables the absolute divergence check.
    ! standard = SYS_INFINITY_DP
    real(DP)                        :: ddivAbs = SYS_INFINITY_DP

    ! INPUT PARAMETER FOR ITERATIVE SOLVERS: 
    ! RHS-vector is treated as zero if max(defect) < drhsZero
    real(DP)                        :: drhsZero = 1E-90_DP

    ! INPUT PARAMETER FOR ITERATIVE SOLVERS: 
    ! Type of stopping criterion to use. One of the
    ! SPTILS_STOP_xxxx constants.
    integer                    :: istoppingCriterion = SPTILS_STOP_STANDARD

    ! INPUT PARAMETER FOR ITERATIVE SOLVERS: 
    ! Minimum number of iterations top perform
    integer                    :: nminIterations = 1

    ! INPUT PARAMETER FOR ITERATIVE SOLVERS: 
    ! Maximum number of iterations top perform
    integer                    :: nmaxIterations = 50
    
    ! INPUT PARAMETER FOR SOLVERS WITH RESIDUAL CHECK: 
    ! Perform residual checks
    ! YES: check residuals (default), 
    ! NO: don't check residuals, simply perform as many iterations as 
    ! configured by nminIterations. Is typically set to NO for smoothing
    ! with a solver.
    integer                    :: iresCheck = YES

    ! INPUT PARAMETER FOR SOLVERS WITH RESIDUAL CHECK: 
    ! Type of norm to use in the residual checking (cf. linearalgebra.f90).
    ! =0: euclidian norm, =1: l1-norm, =2: l2-norm, =3: MAX-norm
    integer                    :: iresNorm = 2
    
    ! INPUT PARAMETER: Output level
    ! This determines the output level of the solver.
    ! =0: no output, =1: basic output, =2, extended output
    integer                    :: ioutputLevel = 2

    ! INPUT PARAMETER FOR ITERATIVE SOLVERS WITH RESIDUAL CHECK:
    ! Number of iterations to perform before printing out the
    ! norm of the residual to screen.
    ! =1: Print residual in every iteration
    integer                    :: niteResOutput = 1

    ! READ ONLY: Algorithm identifier. 
    ! One of the SPTILS_ALG_xxxx constants.
    ! Depending on the value, a solver-specific structure may be 
    ! assigned to this structure below.
    integer                    :: calgorithm = SPTILS_ALG_UNDEFINED
    
    ! READ ONLY: Solver ability tag. 
    ! Bitfield. A combination of SPTILS_ABIL_xxxx
    ! flags that specify the ability of the actual solver (e.g.
    ! whether it can handle block-matrices or only scalar matrices,...).
    ! Calling a solver that is not able to handle a specific problem
    ! will cause an error by the framework.
    integer(I32)                    :: ccapability = 0

    ! A pointer to the problem structure that defines the problem.
    type(t_problem), pointer   :: p_rproblem
    
    ! STATUS FOR ITERATIVE SOLVERS: Current iteration
    integer                    :: icurrentIteration
    
    ! STATUS FOR ITERATIVE SOLVERS: Last defect
    real(DP)                   :: dlastDefect
    
    ! STATISTICS OUTPUT: Time for solving problems in space.
    type(t_timer)              :: rtimeSpacePrecond
    
    ! Structure of the space time matrix that the solver should use.
    type(t_ccoptSpaceTimeMatrix) :: rmatrix
    
    ! Pointer to a structure for the Defect correction solver; NULL() if not set
    type (t_sptilsSubnodeDefCorr), pointer        :: p_rsubnodeDefCorr     => null()

    ! Pointer to a structure for the CG solver; NULL() if not set
    type (t_sptilsSubnodeCG), pointer             :: p_rsubnodeCG     => null()

    ! Pointer to a structure for the BiCGStab solver; NULL() if not set
    type (t_sptilsSubnodeBiCGStab), pointer       :: p_rsubnodeBiCGStab => null()

    ! Pointer to a structure for the block Jacobi preconditioner
    type(t_sptilsSubnodeBlockJacobi), pointer     :: p_rsubnodeBlockJacobi => null()

    ! Pointer to a structure for the block Gauss-Seidel preconditioner
    type(t_sptilsSubnodeBlockFBSOR), pointer        :: p_rsubnodeBlockFBSOR => null()

    ! Pointer to a structure for the multigrid preconditioner
    type(t_sptilsSubnodeMultigrid), pointer       :: p_rsubnodeMultigrid => null()

    ! Pointer to a structure for the UMFPACK4 preconditioner
    type(t_sptilsSubnodeUMFPACK4), pointer        :: p_rsubnodeUMFPACK4 => null()

    ! Pointer to a structure for the UMFPACK4 preconditioner
    type(t_sptilsSubnodeBlockFBGS), pointer       :: p_rsubnodeBlockFBGS => null()

  end type
  
!</typeblock>

! *****************************************************************************

!<typeblock>
  
  ! This structure realises the subnode for the Defect correction solver.
  ! The entry p_rpreconditioner points either to NULL() or to another
  ! t_sptilsNode structure for the solver that realises the 
  ! preconditioning.
  
  type t_sptilsSubnodeDefCorr
  
    ! A pointer to the solver node for the preconditioner or NULL(),
    ! if no preconditioner is used.
    type(t_sptilsNode), pointer       :: p_rpreconditioner            => null()
    
    ! A temporary space-time vector.
    type(t_spacetimeVector)           :: rtempVector
    
    ! A temporary space-time vector.
    type(t_spacetimeVector)           :: rtempVector2
    
  end type
  
!</typeblock>

! *****************************************************************************

!<typeblock>
  
  ! This structure realises the subnode for the CG solver.
  ! The entry p_rpreconditioner points either to NULL() or to another
  ! t_linsolNode structure for the solver that realises the 
  ! preconditioning.
  
  type t_sptilsSubnodeCG
  
    ! Temporary vectors to use during the solution process
    type(t_spacetimeVector), dimension(4) :: RtempVectors

    ! A pointer to the solver node for the preconditioner or NULL(),
    ! if no preconditioner is used.
    type(t_sptilsNode), pointer       :: p_rpreconditioner            => null()
  
    ! A temporary space-time vector.
    type(t_vectorBlock)               :: rtempVectorSpace
  
  end type
  
!</typeblock>

! *****************************************************************************

!<typeblock>
  
  ! This structure realises the subnode for the BiCGStab solver.
  ! The entry p_rpreconditioner points either to NULL() or to another
  ! t_linsolNode structure for the solver that realises the 
  ! preconditioning.
  
  type t_sptilsSubnodeBiCGStab
  
    ! Temporary vectors to use during the solution process
    type(t_spaceTimeVector), dimension(6) :: RtempVectors
    
    ! Another temporary vector for right- and symmetrical preconditioning
    type(t_vectorBlock) :: rprecondTemp
    
    ! A pointer to the solver node for the preconditioner or NULL(),
    ! if no preconditioner is used.
    type(t_sptilsNode), pointer       :: p_rpreconditioner            => null()

    ! A temporary space-time vector.
    type(t_vectorBlock)               :: rtempVectorSpace
    
  end type
  
!</typeblock>

! *****************************************************************************

!<typeblock>
  
  ! This structure realises the subnode for the Block Jacobi solver.
  
  type t_sptilsSubnodeBlockJacobi

    ! Pointer to a spatial preconditioner structure that defines the
    ! preconditioning in each substep of the global system.
    type(t_ccspatialPreconditioner), pointer :: p_rspatialPreconditioner

  end type
  
!</typeblock>

! *****************************************************************************

!<typeblock>
  
  ! This structure realises the subnode for the Block Gauss-Seidel solver.
  
  type t_sptilsSubnodeBlockFBSOR

    ! Pointer to a spatial preconditioner structure that defines the
    ! preconditioning in each substep.
    type(t_ccspatialPreconditioner), pointer :: p_rspatialPreconditioner

    ! A temporary space-time vector.
    type(t_spacetimeVector)           :: rtempVector

    ! SOR-Damping parameter for damping the time dependence
    real(DP) :: domegaSOR
    
    ! Only partial update of the primal/dual solution
    logical  :: bpartialUpdate

  end type
  
!</typeblock>

! *****************************************************************************

!<typeblock>
  
  ! This structure realises the subnode for the Block Gauss-Seidel solver.
  
  type t_sptilsSubnodeBlockFBGS

    ! Pointer to a spatial preconditioner structure that defines the
    ! preconditioning in each substep.
    type(t_ccspatialPreconditioner), pointer :: p_rspatialPreconditioner

    ! A temporary space-time vector.
    type(t_spacetimeVector)           :: rtempVector
    
    ! SOR-Damping parameter for damping the time dependence
    real(DP) :: domegaSOR

  end type
  
!</typeblock>

! *****************************************************************************

!<typeblock>
  
  ! Level data for multigrid. This structure forms an entry in a linked list
  ! of levels which multigrid uses for solving the given problem.
  ! The solver subnode t_sptilsSubnodeMultigrid holds pointers to the head
  ! and the tail of the list.
  
  type t_sptilsMGLevelInfo
  
    ! Structure that defines the space time matrix on that level.
    type(t_ccoptSpaceTimeMatrix)        :: rmatrix
    
    ! A RHS vector for that level. The memory for this is allocated
    ! in initStructure and released in doneStructure.
    type(t_spacetimeVector)             :: rrhsVector

    ! A solution vector for that level. The memory for this is allocated
    ! in initStructure and released in doneStructure.
    type(t_spacetimeVector)             :: rsolutionVector

    ! A temporary vector for the adaptive coarse grid correction
    type(t_spacetimeVector)             :: rtempCGCvector

    ! A temporary vector for that level. The memory for this is allocated
    ! in initStructure and released in doneStructure.
    type(t_spacetimeVector)             :: rtempVector

    ! A pointer to the solver node for the presmoothing or NULL(),
    ! if no presmoother is used.
    type(t_sptilsNode), pointer         :: p_rpresmoother         => null()

    ! A pointer to the solver node for the postsmoothing or NULL(),
    ! if no presmoother is used.
    type(t_sptilsNode), pointer         :: p_rpostsmoother        => null()

    ! A pointer to the solver node for the coarse grid solver if
    ! the level corresponding to this structure represents the coarse
    ! grid. NULL() otherwise.
    type(t_sptilsNode), pointer         :: p_rcoarseGridSolver    => null()
    
    ! Pointer to a structure for the CG solver; NULL() if not set
    type (t_sptilsSubnodeCG), pointer   :: p_rsubnodeCG          => null()

    ! An interlevel projection structure that configures the transport
    ! of defect/solution vectors from one level to another.
    ! For the coarse grid, this structure is ignored.
    ! For a finer grid (e.g. level 4), this defines the grid transfer
    ! between the current and the lower level (so here between level 3 and
    ! level 4).
    type(t_sptiProjection)   :: rinterlevelProjection
    
    ! Relative daptive cycle convergence criterion for coarse levels.
    ! This value is usually =1E99_DP which deactivates adaptive cycles.
    ! The user can set this variable on initialisation of Multigrid to
    ! a value < 1E99_DP. In this case, the complete multigrid cycle on all 
    ! levels except for the fine grid is repeated until
    !  |res. after postsmoothing| < depsRelCycle * |initial res on that level|.
    ! This allows 'adaptive cycles' which e.g. gain one digit on a coarse
    ! level before prolongating the solution to the fine grid.
    ! This is an extension to the usual F/V/W-cycle scheme.
    real(DP)                      :: depsRelCycle             = 1E99_DP

    ! Absolute adaptive cycle convergence criterion for coarse levels.
    ! This value is usually =1E99_DP which deactivates adaptive cycles.
    ! The user can set this variable on initialisation of Multigrid to
    ! a value < 1E99_DP. In this case, the complete multigrid cycle on all 
    ! levels except for the fine grid is repeated until
    !  |res. after postsmoothing| < depsAbsCycle.
    ! This allows 'adaptive cycles' which e.g. gain one digit on a coarse
    ! level before prolongating the solution to the fine grid.
    ! This is an extension to the usual F/V/W-cycle scheme.
    real(DP)                      :: depsAbsCycle             = 1E99_DP
    
    ! Maximum number of adaptive cycles performed on that level.
    ! Only used if adaptive cycles are activated by setting 
    ! depsRelCycle < 1E99_DP. -1=infinity=standard
    integer                       :: nmaxAdaptiveCycles       = -1
    
    ! STATUS/INTERNAL: A temporary vector used for prolongation/restriction
    type(t_vectorBlock)      :: rprjVector
    
    ! STATUS/INTERNAL: MG cycle information.
    ! Number of cycles to perform on this level.
    integer                        :: ncycles

    ! STATUS/INTERNAL: MG cycle information.
    ! Number of remaining cycles to perform on this level.
    integer                        :: ncyclesRemaining
    
    ! STATUS/INTERNAL: initial residuum when a solution is restricted to this
    ! level. Only used if adaptive cycles are activated by setting 
    ! depsRelCycle < 1E99_DP in t_linsolSubnodeMultigrid.
    real(DP)                       :: dinitResCycle = 0.0_DP
    
    ! STATUS/INTERNAL: Number of current cycle on that level. 
    ! Only used if adaptive cycles are activated by setting 
    ! depsRelCycle < 1E99_DP in t_sptilsSubnodeMultigrid.
    integer                        :: icycleCount   = 0
    
  end type
  
!</typeblock>

  ! ***************************************************************************

!<typeblock>
  
  ! This structure realises the subnode for the Multigrid solver.
  ! There are pointers to the head and the tail of a linked list of
  ! t_mgLevelInfo structures which hold the data of all levels.
  ! The tail of the list corresponds to the maximum level where
  ! multigrid solves the problem.
  
  type t_sptilsSubnodeMultigrid
  
    ! INPUT PARAMETER: Cycle identifier. 
    !  0=F-cycle, 
    !  1=V-cycle, 
    !  2=W-cycle.
    integer                       :: icycle                   = 0
    
    ! INPUT PARAMETER: Minimum level in p_Rlevels.
    integer                       :: NLMIN                    = 1
    
    ! INPUT PARAMETER: Maximum level in p_Rlevels.
    integer                       :: NLMAX                    = 1

    ! INPUT PARAMETER: Minimum value for adaptive coarse grid correction.
    ! A value of dalphamin=dalphamax deactivates the adaptive coarse
    ! grid correction. Standard = 1.0.
    real(DP)                      :: dalphamin                = 1.0_DP
    
    ! INPUT PARAMETER: Maximum value for adaptive coarse grid correction.
    ! A value of dalphamin=dalphamax deactivates the adaptive coarse
    ! grid correction. Standard = 1.0.
    real(DP)                      :: dalphamax                = 1.0_DP
    
    ! A list of weights for the different equations in each timestep.
    ! If not associated, a standard value of 1.0 is assumed for every equation.
    ! If associated, the calculation routine for the optimal
    ! coarse grid correction will multiply residuals of the corresponding
    ! equation with this factor.
    ! Example: 2D Navier-Stokes: (1,1,-1)
    ! -> The pressure equation is multiplied by -1 before taking the norm
    !    of the residual.
    ! Can be used to symmetrise equations.
    real(DP), dimension(:), pointer :: p_DequationWeights => null()

    ! Array of t_sptilsMGLevelInfo structures for all the levels the MG solver
    ! should handle.
    type(t_sptilsMGLevelInfo), dimension(:), pointer :: p_Rlevels => null()
    
    ! STATISTICS OUTPUT: Total time needed for smoothing operations
    type(t_timer) :: rtimeSmoothing
        
    ! STATISTICS OUTPUT: Total time needed for the coarse grid solver
    type(t_timer) :: rtimeCoarseGridSolver
    
    ! STATISTICS OUTPUT: Time needed for linear algebra stuff (matrix-vector, 
    ! vector-copy, prolongation/restriction,...)
    type(t_timer) :: rtimeLinearAlgebra
    
    ! STATISTICS OUTPUT: Time needed for prolongation/restriction
    type(t_timer) :: rtimeProlRest
    
  end type
  
!</typeblock>

! *****************************************************************************

!<typeblock>
  
  ! This structure realises the subnode for the UMFPACK4 solver.
  
  type t_sptilsSubnodeUMFPACK4

    ! DEBUG flag: With this flag, UMFPACK can be told to write out the system
    ! matrix in initData after assembling it.
    ! =0: Don't write anything.
    ! =1: Write out the matrix to a text file.
    ! =2: Write out the matrix to a text file and stop the program afterwards.
    integer :: cwriteMatrix = 0
    
    ! File name of the file receiving the system matrix in human readable
    ! form if bwriteMatrix=true.
    character(len=SYS_STRLEN) :: sfilename = "matrix.txt"

    ! When writing a matrix: X/Y-position of the first block  or (0,0) for auto.
    integer :: ixfirst = 0
    integer :: iyfirst = 0

    ! When writing a matrix: X/Y-position of the last block  or (0,0) for auto.
    integer :: ixlast = 0
    integer :: iylast = 0

    ! <!--
    ! The following variables are internal variables and not to be used
    ! by any application
    ! -->
  
    ! Control structure for UMFPACK4; contains parameter for the solver
    real(DP), dimension(20) :: Dcontrol

    ! Handle for symbolic factorisation.
    ! This is not a FEAT-Handle!
    integer(I32) :: isymbolic = 0

    ! Handle for numeric factorisation
    ! This is not a FEAT-Handle!
    integer(I32) :: inumeric = 0

  end type
  
!</typeblock>

!</types>

contains

  ! ***************************************************************************
  
!<subroutine>
  
  recursive subroutine sptils_setMatrices (rsolverNode,Rmatrices)
  
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
  type(t_ccoptSpaceTimeMatrix), dimension(:), intent(IN) :: Rmatrices
!</input>
  
!<inputoutput>
  ! The solver node which should be initialised
  type(t_sptilsNode), intent(INOUT)             :: rsolverNode
!</inputoutput>
  
!</subroutine>

  ! Copy the matrix structure on the finest level to the rsolverNode
  ! structure. This corresponds to the system we want to solve.
  rsolverNode%rmatrix = Rmatrices(ubound(Rmatrices,1))
  
  ! Depending on the solver type, call the corresponding initialisation
  ! routine. For all single-grid solvers, pass the matrix on the highest
  ! level in RlinsysInfo. For multigrid we pass all matrices.
  ! That way the solvers can decide to do anything else with the matrix
  ! or to initialise their subsolvers.
  !
  ! For single grid solvers, the whole Rmatrices array is usually also
  ! passed to the initialisation routine. This allows to initialise
  ! the sub-nodes of each solver.
  
  select case(rsolverNode%calgorithm)
  case (SPTILS_ALG_DEFCORR)
    call sptils_setMatrixDefCorr (rsolverNode,Rmatrices)
  case (SPTILS_ALG_CG)
    call sptils_setMatrixCG (rsolverNode,Rmatrices)
  case (SPTILS_ALG_BiCGStab)
    call sptils_setMatrixBiCGStab (rsolverNode,Rmatrices)
  case (SPTILS_ALG_BLOCKJACOBI)
    call sptils_setMatrixBlockJacobi (rsolverNode,Rmatrices)
  case (SPTILS_ALG_BlockFBSOR)
    call sptils_setMatrixBlockFBSOR (rsolverNode,Rmatrices)
  case (SPTILS_ALG_BLOCKFBGS)
    call sptils_setMatrixBlockFBGS (rsolverNode,Rmatrices)
  case (SPTILS_ALG_MULTIGRID)
    call sptils_setMatrixMultigrid (rsolverNode,Rmatrices)
  case DEFAULT
  end select

  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>
  
  recursive subroutine sptils_initStructure (rsolverNode, ierror)
  
!<description>
  
  ! Initialises the problem structure in the solver node rsolverNode by calling
  ! the initialisation routine of the appropriate solver. The solver
  ! initialisation routine itself can call this procedure to initialise
  ! its sub-solver nodes.
  
!</description>
  
!<inputoutput>
  ! The solver node which should be initialised
  type(t_sptilsNode), intent(INOUT)                     :: rsolverNode
!</inputoutput>

!<output>
  ! One of the SPTILS_ERR_XXXX constants. A value different to 
  ! SPTILS_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  integer, intent(OUT) :: ierror
!</output>
  
!</subroutine>

    ! A-priori we have no error...
    ierror = SPTILS_ERR_NOERROR
    
    ! Call the structure-init routine of the specific solver
    
    select case(rsolverNode%calgorithm)
    case (SPTILS_ALG_DEFCORR)
      call sptils_initStructureDefCorr (rsolverNode,ierror)
    case (SPTILS_ALG_CG)
      call sptils_initStructureCG (rsolverNode,ierror)
    case (SPTILS_ALG_BICGSTAB)
      call sptils_initStructureBiCGStab (rsolverNode,ierror)
    case (SPTILS_ALG_BLOCKJACOBI)
      call sptils_initStructureBlockJacobi (rsolverNode,ierror)
    case (SPTILS_ALG_BlockFBSOR)
      call sptils_initStructureBlockFBSOR (rsolverNode,ierror)
    case (SPTILS_ALG_BLOCKFBGS)
      call sptils_initStructureBlockFBGS (rsolverNode,ierror)
    case (SPTILS_ALG_MULTIGRID)
      call sptils_initStructureMultigrid (rsolverNode,ierror)
    case (SPTILS_ALG_UMFPACK4)
      call sptils_initStructureUMFPACK4 (rsolverNode,ierror)
    end select
  
  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>
  
  recursive subroutine sptils_initData (rsolverNode, ierror)
  
!<description>
  ! Initialises the problem data in the solver node rsolverNode by calling
  ! the initialisation routine of the appropriate solver. The solver
  ! initialisation routine itself can call this procedure to initialise
  ! its sub-solver nodes.
  ! The initialisation of the problem structure allows the solver component
  ! to perform some 'precalculation', e.g. the UMFPACK4 or ILU solver can 
  ! perform a numerical factorisation. The problem structure usually does
  ! not change during a simulation, except when the grid moves e.g.
!</description>
  
!<output>
  ! One of the SPTILS_ERR_XXXX constants. A value different to 
  ! SPTILS_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  integer, intent(OUT) :: ierror
!</output>

!<inputoutput>
  ! The solver node which should be initialised
  type(t_sptilsNode), intent(INOUT)                     :: rsolverNode
!</inputoutput>

!</subroutine>

    ! A-priori we have no error...
    ierror = SPTILS_ERR_NOERROR

    ! Call the data-init routine of the specific solver
    
    select case(rsolverNode%calgorithm)
    case (SPTILS_ALG_DEFCORR)
      call sptils_initDataDefCorr (rsolverNode,ierror)
    case (SPTILS_ALG_CG)
      call sptils_initDataCG (rsolverNode,ierror)
    case (SPTILS_ALG_BICGSTAB)
      call sptils_initDataBiCGStab (rsolverNode,ierror)
    case (SPTILS_ALG_BLOCKJACOBI)
      call sptils_initDataBlockJacobi (rsolverNode,ierror)
    case (SPTILS_ALG_BlockFBSOR)
      call sptils_initDataBlockFBSOR (rsolverNode,ierror)
    case (SPTILS_ALG_BLOCKFBGS)
      call sptils_initDataBlockFBGS (rsolverNode,ierror)
    case (SPTILS_ALG_MULTIGRID)
      call sptils_initDataMultigrid (rsolverNode,ierror)
    case (SPTILS_ALG_UMFPACK4)
      call sptils_initDataUMFPACK4 (rsolverNode,ierror)
    end select
  
  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>
  
  recursive subroutine sptils_updateStructure (rsolverNode, ierror)
  
!<description>
  
  ! Reinitialises the problem structure in the solver node rsolverNode
  
!</description>
  
!<inputoutput>
  ! The solver node which should be reinitialised
  type(t_sptilsNode), intent(INOUT)                     :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the SPTILS_ERR_XXXX constants. A value different to 
  ! SPTILS_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  integer, intent(OUT) :: ierror
!</output>
  
!</subroutine>

    ! A-priori we have no error...
    ierror = SPTILS_ERR_NOERROR

    ! For now, call the structure-done- and structure-init routines.
    ! Maybe rewritten later for higher performance or special cases...
    
    call sptils_doneStructure (rsolverNode)
    call sptils_initStructure (rsolverNode,ierror)
  
  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>
  
  recursive subroutine sptils_updateData (rsolverNode,ierror)
  
!<description>
  ! Reinitialises the problem data in the solver node rsolverNode.
!</description>
  
!<output>
  ! One of the SPTILS_ERR_XXXX constants. A value different to 
  ! SPTILS_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  integer, intent(OUT) :: ierror
!</output>

!<inputoutput>
  ! The solver node containing the solver confuguration
  type(t_sptilsNode), intent(INOUT)                     :: rsolverNode
!</inputoutput>
  
!</subroutine>

    ! A-priori we have no error...
    ierror = SPTILS_ERR_NOERROR

    ! For now, call the data-done- and data-init routines.
    ! Maybe rewritten later for higher performance or special cases...
    
    call sptils_doneData (rsolverNode)
    call sptils_initData (rsolverNode,ierror)

  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>
  
  recursive subroutine sptils_doneData (rsolverNode)
  
!<description>
  ! Releases the problem data in the solver node rsolverNode.
!</description>
  
!<inputoutput>
  
  ! The solver node containing the solver confuguration
  type(t_sptilsNode), intent(INOUT)                     :: rsolverNode
  
!</inputoutput>
  
!</subroutine>

    ! Call the data-done routine of the specific solver
    
    select case(rsolverNode%calgorithm)
    case (SPTILS_ALG_DEFCORR)
      call sptils_doneDataDefCorr (rsolverNode)
    case (SPTILS_ALG_CG)
      call sptils_doneDataCG (rsolverNode)
    case (SPTILS_ALG_BICGSTAB)
      call sptils_doneDataBiCGStab (rsolverNode)
    case (SPTILS_ALG_BLOCKJACOBI)
      call sptils_doneDataBlockJacobi (rsolverNode)
    case (SPTILS_ALG_BlockFBSOR)
      call sptils_doneDataBlockFBSOR (rsolverNode)
    case (SPTILS_ALG_BLOCKFBGS)
      call sptils_doneDataBlockFBGS (rsolverNode)
    case (SPTILS_ALG_MULTIGRID)
      call sptils_doneDataMultigrid (rsolverNode)
    case (SPTILS_ALG_UMFPACK4)
      call sptils_doneDataUMFPACK4 (rsolverNode)
    end select

  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>
  
  recursive subroutine sptils_doneStructure (rsolverNode)
  
!<description>
  ! Releases the problem structure in the solver node rsolverNode.
  ! This is done by calling the appropriate routine of the
  ! actual solver.
!</description>
  
!<inputoutput>
  ! The solver node which should be reinitialised
  type(t_sptilsNode), intent(INOUT)                     :: rsolverNode
!</inputoutput>
  
!</subroutine>

    ! Call the data-done routine of the specific solver
    
    select case(rsolverNode%calgorithm)
    case (SPTILS_ALG_DEFCORR)
      call sptils_doneStructureDefCorr (rsolverNode)
    case (SPTILS_ALG_CG)
      call sptils_doneStructureCG (rsolverNode)
    case (SPTILS_ALG_BICGSTAB)
      call sptils_doneStructureBiCGStab (rsolverNode)
    case (SPTILS_ALG_BLOCKJACOBI)
      call sptils_doneStructureBlockJacobi (rsolverNode)
    case (SPTILS_ALG_BlockFBSOR)
      call sptils_doneStructureBlockFBSOR (rsolverNode)
    case (SPTILS_ALG_BLOCKFBGS)
      call sptils_doneStructureBlockFBGS (rsolverNode)
    case (SPTILS_ALG_MULTIGRID)
      call sptils_doneStructureMultigrid (rsolverNode)
    case (SPTILS_ALG_UMFPACK4)
      call sptils_doneStructureUMFPACK4 (rsolverNode)
    end select

  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>
  
  recursive subroutine sptils_releaseSolver (p_rsolverNode,bkeepSolverNode)
  
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
  logical, intent(IN), optional :: bkeepSolverNode
!</input>

!<inputoutput>
  ! The solver node which is to be released. If the node contains attached
  ! subsolvers (preconditioners, smoothers,...) they are also released.
  ! On return, rsolverNode is NULL().
  type(t_sptilsNode), pointer     :: p_rsolverNode
!</inputoutput>
  
!</subroutine>

    if (.not. associated(p_rsolverNode)) then
      ! Warning message, return
      print *,'sptils_releaseSolver warning: Solver note not assigned!'
      return
    end if

    ! Depending on the solver type, call the corresponding done-routine
    select case(p_rsolverNode%calgorithm)
    case (SPTILS_ALG_DEFCORR)
      call sptils_doneDefCorr (p_rsolverNode)
    case (SPTILS_ALG_CG)
      call sptils_doneCG (p_rsolverNode)
    case (SPTILS_ALG_BICGSTAB)
      call sptils_doneBiCGStab (p_rsolverNode)
    case (SPTILS_ALG_BLOCKJACOBI)
      call sptils_doneBlockJacobi (p_rsolverNode)
    case (SPTILS_ALG_BlockFBSOR)
      call sptils_doneBlockFBSOR (p_rsolverNode)
    case (SPTILS_ALG_BLOCKFBGS)
      call sptils_doneBlockFBGS (p_rsolverNode)
    case (SPTILS_ALG_MULTIGRID)
      call sptils_doneMultigrid (p_rsolverNode)
    case (SPTILS_ALG_UMFPACK4)
      call sptils_doneUMFPACK4 (p_rsolverNode)
    case DEFAULT
    end select
    
    ! Clean up the associated matrix structure.
    ! Of course, the memory of the matrix is not released from memory, because
    ! if there's a matrix attached, it belongs to the application, not to the
    ! solver!
    !CALL lsysbl_releaseMatrix(p_rsolverNode%rsystemMatrix)
    
    ! Finally release the node itself (if we are allowed to).
    ! Deallocate the structure (if we are allowed to), finish.
    if (.not. present(bkeepSolverNode)) then
      deallocate(p_rsolverNode)
    else
      if (.not. bkeepSolverNode) then
        deallocate(p_rsolverNode)
      end if
    end if
  
  end subroutine
  
  ! ***************************************************************************

!<function>
  
  logical function sptils_testConvergence (rsolverNode, dvecNorm, rdef) result(loutput)
  
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
  type(t_sptilsNode), intent(IN) :: rsolverNode
  
  ! OPTIONAL: The defect vector which norm should be tested.
  ! If existent, the norm of the vector is returned in dvecNorm.
  ! If not existent, the routine assumes that dvecNrm is the norm
  ! of the vector and checks convergence depending on dvecNorm.
  type(t_spaceTimeVector), intent(IN), optional :: rdef
!</input>

!<inputoutput>
  ! Norm of the defect vector. 
  ! If rdef if present, the routine will calculate the norm of rdef and return
  ! it in dvecNorm.
  ! If rdef is not present, dvecNorm is assumed to be a valid norm of a
  ! vector and convergence is tested using dvecNorm.
  real(DP), intent(INOUT) :: dvecNorm
!</inputoutput>
  
!</function>

    ! Calculate the norm of the vector or take the one given
    ! as parameter
    if (present(rdef)) then
      dvecNorm = sptivec_vectorNorm (rdef,rsolverNode%iresNorm)
    end if
    
    ! Relative difference in the residuals small enough?
    if (rsolverNode%depsDiff .ne. 0.0_DP) then
      if (abs(dvecNorm-rsolverNode%dlastDefect) .lt. &
              rsolverNode%depsDiff*rsolverNode%dlastDefect) then
        loutput = .true.
        return
      end if
    end if
    
    select case (rsolverNode%istoppingCriterion)
    
    case (SPTILS_STOP_ONEOF)
      ! Iteration stops if either the absolute or the relative criterium holds.
      loutput = .false.
      
      ! Absolute convergence criterion? Check the norm directly.
      if (rsolverNode%depsAbs .ne. 0.0_DP) then
        if (.not. (dvecNorm .gt. rsolverNode%depsAbs)) then
          loutput = .true.
          return
        end if
      end if
      
      ! Relative convergence criterion? Multiply with initial residuum
      ! and check the norm. 
      if (rsolverNode%depsRel .ne. 0.0_DP) then
        if (.not. &
            (dvecNorm .gt. rsolverNode%depsRel * rsolverNode%dinitialDefect)) then
          loutput = .true.
          return
        end if
      end if
      
    case DEFAULT
      ! Standard stopping criterion.
      ! Iteration stops if both the absolute and the relative criterium holds.
      loutput = .true.
      
      ! Absolute convergence criterion? Check the norm directly.
      if (rsolverNode%depsAbs .ne. 0.0_DP) then
        if (dvecNorm .gt. rsolverNode%depsAbs) then
          loutput = .false.
          return
        end if
      end if
      
      ! Relative convergence criterion? Multiply with initial residuum
      ! and check the norm. 
      if (rsolverNode%depsRel .ne. 0.0_DP) then
        if (dvecNorm .gt. rsolverNode%depsRel * rsolverNode%dinitialDefect) then
          loutput = .false.
          return
        end if
      end if
    end select
    
  end function
  
  ! ***************************************************************************

!<function>
  
  logical function sptils_testDivergence (rsolverNode, dvecNorm, rdef) result(loutput)
  
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
  type(t_sptilsNode), intent(IN) :: rsolverNode
  
  ! OPTIONAL: The defect vector which norm should be tested.
  ! If existent, the norm of the vector is returned in dvecNorm.
  ! If not existent, the routine assumes that dvecNrm is the norm
  ! of the vector and checks divergence depending on dvecNorm.
  type(t_spaceTimeVector), intent(IN), optional :: rdef
!</input>

!<inputoutput>
  ! Norm of the defect vector. 
  ! If rdef if present, the routine will calculate the norm of rdef and return
  ! it in dvecNorm.
  ! If rdef is not present, dvecNorm is assumed to be a valid norm of a
  ! vector and divergence is tested using dvecNorm.
  real(DP), intent(INOUT) :: dvecNorm
!</inputoutput>

!</function>

    ! Calculate the norm of the vector if not given
    ! as parameter
    if (present(rdef)) then
      dvecNorm = sptivec_vectorNorm (rdef,rsolverNode%iresNorm)
    end if
    
    loutput = .false.
    
    ! Absolute divergence criterion? Check the norm directly.
    if (rsolverNode%ddivAbs .ne. SYS_INFINITY_DP) then
     
      ! use NOT here - gives a better handling of special cases like NaN!
      if ( .not. (dvecNorm .le. rsolverNode%ddivAbs)) then
        loutput = .true.
      end if
      
    end if
    
    ! Relative divergence criterion? Multiply with initial residuum
    ! and check the norm. 
    if (rsolverNode%depsRel .ne. SYS_INFINITY_DP) then
      if ( .not. (dvecNorm .le. rsolverNode%dinitialDefect*rsolverNode%ddivRel) ) then
        loutput = .true.
      end if
    end if
  
  end function
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine sptils_precondDefect (rsolverNode,rd)
  
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
  type(t_sptilsNode), intent(INOUT)                :: rsolverNode
  
  ! On call to this routine: The defect vector to be preconditioned.
  ! Will be overwritten by the preconditioned defect.
  type(t_spacetimeVector), intent(INOUT)               :: rd
!</inputoutput>
  
!</subroutine>

    ! Select the solver as configured in rsolverNode and let it perform
    ! the actual preconditioning task.
    
    select case(rsolverNode%calgorithm)
    case (SPTILS_ALG_DEFCORR)
      call sptils_precDefCorr (rsolverNode,rd)
    case (SPTILS_ALG_CG)
      call sptils_precCG (rsolverNode,rd)
    case (SPTILS_ALG_BICGSTAB)
      call sptils_precBiCGStab (rsolverNode,rd)
    case (SPTILS_ALG_BLOCKJACOBI)
      call sptils_precBlockJacobi (rsolverNode,rd)
    case (SPTILS_ALG_BlockFBSOR)
      call sptils_precBlockFBSOR (rsolverNode,rd)
    case (SPTILS_ALG_BLOCKFBGS)
      call sptils_precBlockFBGS (rsolverNode,rd)
    case (SPTILS_ALG_MULTIGRID)
      call sptils_precMultigrid (rsolverNode,rd)
    case (SPTILS_ALG_UMFPACK4)
      call sptils_precUMFPACK4 (rsolverNode,rd)
    case DEFAULT
      print *,'Unknown space-time preconditioner!'
      call sys_halt()
    end select

  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  subroutine sptils_initSolverGeneral (rproblem,p_rsolverNode)
  
!<description>
  ! Creates a new solver node p_rsolverNode on the heap with default
  ! settings and returns a pointer to it.
!</description>
  
!<input>
  ! The problem structure that defines the problem.
  ! A pointer to this is saved in the solver node, so the problem structure
  ! must not be released before the solver node is released.
  type(t_problem), target :: rproblem
!</input>
  
!<output>
  ! A pointer to a new solver node on the heap.
  type(t_sptilsNode), pointer         :: p_rsolverNode
!</output>
  
!</subroutine>

    ! Allocate the node - the default initialisation will set most of the
    ! parameters in the structure.
    allocate(p_rsolverNode)
    
    ! Save a pointer to the problem structure in the solver node
    p_rsolverNode%p_rproblem => rproblem

  end subroutine
  
! *****************************************************************************
! Routines for the Defect Correction iteration
! *****************************************************************************

!<subroutine>
  
  recursive subroutine sptils_initDefCorr (rproblem,p_rsolverNode,p_rpreconditioner)
  
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
  type(t_problem), target :: rproblem

  ! OPTIONAL: A pointer to the solver structure of a solver that should be 
  ! used for preconditioning. If not given or set to NULL(), no preconditioning 
  ! will be used.
  type(t_sptilsNode), pointer, optional   :: p_rpreconditioner
!</input>
  
!<output>
  ! A pointer to a t_sptilsNode structure. Is set by the routine, any previous
  ! value of the pointer is destroyed.
  type(t_sptilsNode), pointer         :: p_rsolverNode
!</output>
  
!</subroutine>
  
    ! Create a default solver structure
    call sptils_initSolverGeneral(rproblem,p_rsolverNode)
    
    ! Initialise the type of the solver
    p_rsolverNode%calgorithm = SPTILS_ALG_DEFCORR
    
    ! Initialise the ability bitfield with the ability of this solver:
    p_rsolverNode%ccapability = SPTILS_ABIL_CHECKDEF + &
                                SPTILS_ABIL_USESUBSOLVER

    ! Allocate a subnode for our solver.
    ! This initialises most of the variables with default values appropriate
    ! to this solver.
    allocate(p_rsolverNode%p_rsubnodeDefCorr)
    
    ! Attach the preconditioner if given. 
    
    if (present(p_rpreconditioner)) then 
      p_rsolverNode%p_rsubnodeDefCorr%p_rpreconditioner => p_rpreconditioner
    end if
    
  end subroutine

  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine sptils_doneDefCorr (rsolverNode)
  
!<description>
  ! This routine releases all temporary memory for the efect correction solver 
  ! from the heap.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure which is to be cleaned up.
  type(t_sptilsNode), intent(INOUT)         :: rsolverNode
!</inputoutput>
  
!</subroutine>
  
    ! Check if there's a preconditioner attached. If yes, release it.
    if (associated(rsolverNode%p_rsubnodeDefCorr%p_rpreconditioner)) then
      call sptils_releaseSolver(rsolverNode%p_rsubnodeDefCorr%p_rpreconditioner)
    end if

    ! Release memory if still associated
    call sptils_doneDataDefCorr (rsolverNode)
    call sptils_doneStructureDefCorr (rsolverNode)
    
    ! Release the subnode structure
    deallocate(rsolverNode%p_rsubnodeDefCorr)
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine sptils_setMatrixDefCorr (rsolverNode,Rmatrices)
  
!<description>
  ! This routine is called if the system matrix changes.
  ! The routine calls sptils_setMatrices for the preconditioner
  ! to inform also that one about the change of the matrix pointer.
!</description>
  
!<input>
  ! An array of system matrices which is simply passed to the initialisation 
  ! routine of the preconditioner.
  type(t_ccoptSpaceTimeMatrix), dimension(:), intent(IN) :: Rmatrices
!</input>
  
!<inputoutput>
  ! The t_sptilsNode structure of the solver
  type(t_sptilsNode), intent(INOUT)         :: rsolverNode
!</inputoutput>
  
!</subroutine>

    ! Do we have a preconditioner given? If yes, call the general initialisation
    ! routine to initialise it.
    
    if (associated(rsolverNode%p_rsubnodeDefCorr%p_rpreconditioner)) then
      call sptils_setMatrices (rsolverNode%p_rsubnodeDefCorr%p_rpreconditioner, &
                               Rmatrices)
    end if

  end subroutine

! ***************************************************************************

!<subroutine>
  
  recursive subroutine sptils_initStructureDefCorr (rsolverNode,ierror)
  
!<description>
  ! Solver preparation. Perform symbolic factorisation (not of the defect
  ! correcion solver, but of subsolvers). Allocate temporary memory.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure 
  type(t_sptilsNode), intent(INOUT), target :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the SPTILS_ERR_XXXX constants. A value different to 
  ! SPTILS_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  integer, intent(OUT) :: ierror
!</output>
  
!</subroutine>

    ! Local variables
    type(t_ccoptSpaceTimeDiscretisation), pointer :: p_rspaceTimeDiscr
    
    p_rspaceTimeDiscr => rsolverNode%rmatrix%p_rspaceTimeDiscr
    
    ! A-priori we have no error...
    ierror = SPTILS_ERR_NOERROR
    
    ! Allocate memory for the two temp vectors.
    call sptivec_initVector (rsolverNode%p_rsubnodeDefCorr%rtempVector,&
        p_rspaceTimeDiscr%NEQtime,&
        p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation)

    call sptivec_initVector (rsolverNode%p_rsubnodeDefCorr%rtempVector2,&
        p_rspaceTimeDiscr%NEQtime,&
        p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation)
        
    ! We simply pass isubgroup to the subsolvers when calling them.
    ! Inside of this routine, there's not much to do with isubgroup,
    ! as there is no special information (like a factorisation)
    ! associated with this type of solver, which has to be allocated,
    ! released or updated...

    ! Call the init routine of the preconditioner.
    if (associated(rsolverNode%p_rsubnodeDefCorr%p_rpreconditioner)) then
      call sptils_initStructure (rsolverNode%p_rsubnodeDefCorr%p_rpreconditioner,ierror)
    end if
    
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine sptils_initDataDefCorr (rsolverNode, ierror)
  
!<description>
  ! Calls the initData subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with sptils_initData.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure 
  type(t_sptilsNode), intent(INOUT)         :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the SPTILS_ERR_XXXX constants. A value different to 
  ! SPTILS_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  integer, intent(OUT) :: ierror
!</output>

!</subroutine>

    ! Call the init routine of the preconditioner.
    if (associated(rsolverNode%p_rsubnodeDefCorr%p_rpreconditioner)) then
      call sptils_initData (rsolverNode%p_rsubnodeDefCorr%p_rpreconditioner,ierror)
    end if
    
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine sptils_doneDataDefCorr (rsolverNode)
  
!<description>
  ! Calls the doneData subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with sptils_doneData.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure 
  type(t_sptilsNode), intent(INOUT)         :: rsolverNode
!</inputoutput>
  
!</subroutine>

    ! Call the done routine of the preconditioner.
    if (associated(rsolverNode%p_rsubnodeDefCorr%p_rpreconditioner)) then
      call sptils_doneData (rsolverNode%p_rsubnodeDefCorr%p_rpreconditioner)
    end if
    
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine sptils_doneStructureDefCorr (rsolverNode)
  
!<description>
  ! Calls the doneStructure subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with sptils_doneStructure.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure 
  type(t_sptilsNode), intent(INOUT)         :: rsolverNode
!</inputoutput>
  
!</subroutine>

    ! Release the temp vectors.
    ! Note that the vectors may already be released from a previous call.
    if (rsolverNode%p_rsubnodeDefCorr%rtempVector2%NEQtime .ne. 0) &
      call sptivec_releaseVector(rsolverNode%p_rsubnodeDefCorr%rtempVector2)
    if (rsolverNode%p_rsubnodeDefCorr%rtempVector%NEQtime .ne. 0) &
      call sptivec_releaseVector(rsolverNode%p_rsubnodeDefCorr%rtempVector)

    ! Call the init routine of the preconditioner.
    if (associated(rsolverNode%p_rsubnodeDefCorr%p_rpreconditioner)) then
      call sptils_doneStructure (rsolverNode%p_rsubnodeDefCorr%p_rpreconditioner)
    end if
    
  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>
  
  recursive subroutine sptils_precDefCorr (rsolverNode,rd)
  
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
  type(t_sptilsNode), intent(INOUT), target :: rsolverNode

  ! On call to this routine: The defect vector to be preconditioned.
  ! Will be overwritten by the preconditioned defect.
  type(t_spacetimeVector), intent(INOUT)    :: rd
!</inputoutput>
  
!</subroutine>

  ! local variables
  integer :: ite
  real(DP) :: dres,dfr,dresnorm
  type(t_spacetimeVector), pointer :: p_rx,p_rdef
  type(t_sptilsNode), pointer :: p_rprecSubnode
  
  ! Damping parameter
  real(DP) :: domega

  ! The system matrix
  type(t_ccoptSpaceTimeMatrix), pointer :: p_rmatrix
  
  ! Minimum number of iterations, print-sequence for residuals
  integer :: nminIterations, niteResOutput
  
  ! Whether to filter/prcondition
  logical bprec
  
  ! The local subnode
  type(t_sptilsSubnodeDefCorr), pointer :: p_rsubnode

    ! Status reset
    rsolverNode%iresult = 0
    
    call stat_clearTimer (rsolverNode%rtimeSpacePrecond)
    
    ! Getch some information
    p_rsubnode => rsolverNode%p_rsubnodeDefCorr
    p_rmatrix => rsolverNode%rmatrix

    ! Check the parameters
    if (rd%NEQtime .eq. 0) then
    
      ! Parameters wrong
      rsolverNode%iresult = 2
      return
    end if

    ! Minimum number of iterations
 
    nminIterations = max(rsolverNode%nminIterations,0)
    
    ! Damping parameter
    domega = rsolverNode%domega
      
    ! Use preconditioning? 
    bprec = associated(rsolverNode%p_rsubnodeDefCorr%p_rpreconditioner)
    
    if (bprec) then
      p_rprecSubnode => p_rsubnode%p_rpreconditioner
    end if

    ! Iteration when the residuum is printed:

    niteResOutput = max(1,rsolverNode%niteResOutput)

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
    call sptivec_clearVector (p_rx)
  
    ! Copy our RHS rd to p_rdef. As the iteration vector is 0, this
    ! is also our initial defect.

    call sptivec_copyVector(rd,p_rdef)

    ! Get the norm of the residuum
    dres = sptivec_vectorNorm (p_rdef,rsolverNode%iresNorm)
    if (.not.((dres .ge. 1E-99_DP) .and. &
              (dres .le. 1E99_DP))) dres = 0.0_DP

    ! Initialize starting residuum
      
    rsolverNode%dinitialDefect = dres
    rsolverNode%dlastDefect = 0.0_DP

    ! Check if out initial defect is zero. This may happen if the filtering
    ! routine filters "everything out"!
    ! In that case we can directly stop our computation.

    if ( rsolverNode%dinitialDefect .lt. rsolverNode%drhsZero ) then
     
      ! final defect is 0, as initialised in the output variable above

      call sptivec_clearVector(p_rx)
      ite = 0
      rsolverNode%dfinalDefect = dres
          
    else

      if (rsolverNode%ioutputLevel .ge. 2) then
        call output_line ('Space-Time-DefCorr: Iteration '// &
             trim(sys_siL(0,10))//',  !!RES!! = '//&
             trim(sys_sdEL(rsolverNode%dinitialDefect,15)) )
      end if

      ! Perform at most nmaxIterations loops to get a new vector

      do ite = 1,rsolverNode%nmaxIterations
      
        rsolverNode%icurrentIteration = ite
        
        if (bprec) then
          ! Perform preconditioning with the assigned preconditioning
          ! solver structure.
          call sptils_precondDefect (p_rprecSubnode,p_rdef)
          call stat_addTimers (p_rprecSubnode%rtimeSpacePrecond,rsolverNode%rtimeSpacePrecond)
        end if
        
        ! In p_rdef, we now have the current residuum $P^{-1} (b-Ax)$.
        ! Add it (damped) to the current iterate p_x to get
        !   $$ x  :=  x  +  \omega P^{-1} (b-Ax) $$

        call sptivec_vectorLinearComb (p_rdef ,p_rx,domega,1.0_DP)

        ! Calculate the residuum for the next step : (b-Ax).
        ! Dimultaneously filter the defect for bondary conditions.
        call sptivec_copyVector (rd,p_rdef)
        call cc_spaceTimeMatVec (rsolverNode%p_rproblem, p_rmatrix, &
            p_rx,p_rdef, -1.0_DP,1.0_DP,SPTID_FILTER_DEFECT,dresnorm,.false.)
        
        ! Filter the defect for boundary conditions in space and time.
        !
        ! Here not necessary, since this is done in the matrix vector multiplication!
        !CALL tbc_implementInitCondDefect (&
        !    p_rmatrix%p_rspaceTimeDiscr,p_rdef,p_rsubnode%rtempVectorSpace)
        !CALL tbc_implementBCdefect (rsolverNode%p_rproblem,&
        !   p_rmatrix%p_rspaceTimeDiscr,p_rdef,p_rsubnode%rtempVectorSpace)
        
        ! Get the norm of the new (final?) residuum
        dfr = sptivec_vectorNorm (p_rdef,rsolverNode%iresNorm)
     
        rsolverNode%dlastDefect = rsolverNode%dfinalDefect
        rsolverNode%dfinalDefect = dfr

        ! Test if the iteration is diverged
        if (sptils_testDivergence(rsolverNode,dfr)) then
          call output_line ('Space-Time-DefCorr: Solution diverging!')
          rsolverNode%iresult = 1
          exit
        end if
     
        ! At least perform nminIterations iterations
        if (ite .ge. nminIterations) then
        
          ! Check if the iteration converged
          if (sptils_testConvergence(rsolverNode,dfr)) exit
          
        end if

        ! print out the current residuum

        if ((rsolverNode%ioutputLevel .ge. 2) .and. &
            (mod(ite,niteResOutput).eq.0)) then
          call output_line ('Space-Time-DefCorr: Iteration '// &
              trim(sys_siL(ITE,10))//',  !!RES!! = '//&
              trim(sys_sdEL(rsolverNode%dfinalDefect,15)) )
        end if

      end do

      ! Set ITE to NIT to prevent printing of "NIT+1" of the loop was
      ! completed

      if (ite .gt. rsolverNode%nmaxIterations) &
        ite = rsolverNode%nmaxIterations

      ! Finish - either with an error or if converged.
      ! Print the last residuum.

      if ((rsolverNode%ioutputLevel .ge. 2) .and. &
          (ite .ge. 1) .and. (ITE .lt. rsolverNode%nmaxIterations) .and. &
          (rsolverNode%iresult .ge. 0)) then
        call output_line ('Space-Time-DefCorr: Iteration '// &
            trim(sys_siL(ITE,10))//',  !!RES!! = '//&
            trim(sys_sdEL(rsolverNode%dfinalDefect,15)) )
      end if

    end if

    rsolverNode%iiterations = ite
    
    ! Overwrite our previous RHS by the new correction vector p_rx.
    ! This completes the preconditioning.
    call sptivec_copyVector (p_rx,rd)
      
    ! Don't calculate anything if the final residuum is out of bounds -
    ! would result in NaN's,...
      
    if (rsolverNode%dfinalDefect .lt. 1E99_DP) then
    
      ! If the initial defect was zero, the solver immediately
      ! exits - and so the final residuum is zero and we performed
      ! no steps; so the resulting convergence rate stays zero.
      ! In the other case the convergence rate computes as
      ! (final defect/initial defect) ** 1/nit :

      if (rsolverNode%dfinalDefect .gt. rsolverNode%drhsZero) then
        rsolverNode%dconvergenceRate = &
                    (rsolverNode%dfinalDefect / rsolverNode%dinitialDefect) ** &
                    (1.0_DP/real(rsolverNode%iiterations,DP))
      end if
      
      if (rsolverNode%ioutputLevel .ge. 2) then
        call output_lbrk()
        call output_line ('Space-Time-DefCorr statistics:')
        call output_lbrk()
        call output_line ('Iterations              : '//&
             trim(sys_siL(rsolverNode%iiterations,10)) )
        call output_line ('!!INITIAL RES!!         : '//&
             trim(sys_sdEL(rsolverNode%dinitialDefect,15)) )
        call output_line ('!!RES!!                 : '//&
             trim(sys_sdEL(rsolverNode%dfinalDefect,15)) )
        if (rsolverNode%dinitialDefect .gt. rsolverNode%drhsZero) then     
          call output_line ('!!RES!!/!!INITIAL RES!! : '//&
            trim(sys_sdEL(rsolverNode%dfinalDefect / rsolverNode%dinitialDefect,15)) )
        else
          call output_line ('!!RES!!/!!INITIAL RES!! : '//&
               trim(sys_sdEL(0.0_DP,15)) )
        end if
        call output_lbrk ()
        call output_line ('Rate of convergence     : '//&
             trim(sys_sdEL(rsolverNode%dconvergenceRate,15)) )

      end if

      if (rsolverNode%ioutputLevel .eq. 1) then
        call output_line (&
              'Space-Time-DefCorr: Iterations/Rate of convergence: '//&
              trim(sys_siL(rsolverNode%iiterations,10))//' /'//&
              trim(sys_sdEL(rsolverNode%dconvergenceRate,15)) )
      end if
      
    else
      ! DEF=Infinity; RHO=Infinity, set to 1
      rsolverNode%dconvergenceRate = 1.0_DP
    end if  
  
  end subroutine
  
! *****************************************************************************
! Block Jacobi preconditioner
! *****************************************************************************

!<subroutine>
  
  recursive subroutine sptils_initBlockJacobi (rproblem,p_rsolverNode,domega,rspatialPrecond)
  
!<description>
  ! Creates a t_sptilsNode solver structure for the block Jacobi preconditioner.
!</description>
  
!<input>
  ! The problem structure that defines the problem.
  ! A pointer to this is saved in the solver node, so the problem structure
  ! must not be released before the solver node is released.
  type(t_problem), target :: rproblem
  
  ! Damping parameter
  real(DP), intent(IN) :: domega
  
  ! A spatial preconditioner structure that defines how to perform the preconditioning
  ! in each substep. A pointer to this is noted in p_rsolverNode, so rspatialPrecond
  ! should not be released before the solver is destroyed.
  type(t_ccspatialPreconditioner), intent(IN), target :: rspatialPrecond
!</input>
  
!<output>
  ! A pointer to a t_sptilsNode structure. Is set by the routine, any previous
  ! value of the pointer is destroyed.
  type(t_sptilsNode), pointer         :: p_rsolverNode
!</output>
  
!</subroutine>
  
    ! Create a default solver structure
    call sptils_initSolverGeneral(rproblem,p_rsolverNode)
    
    ! Initialise the type of the solver
    p_rsolverNode%calgorithm = SPTILS_ALG_BLOCKJACOBI
    
    ! Initialise the ability bitfield with the ability of this solver:
    p_rsolverNode%ccapability = SPTILS_ABIL_DIRECT

    p_rsolverNode%domega = domega
    
    ! Allocate a subnode for our solver.
    ! This initialises most of the variables with default values appropriate
    ! to this solver.
    allocate(p_rsolverNode%p_rsubnodeBlockJacobi)
    
    ! Attach the preconditioner if given. 
    p_rsolverNode%p_rsubnodeBlockJacobi%p_rspatialPreconditioner => rspatialPrecond
    
  end subroutine

  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine sptils_doneBlockJacobi (rsolverNode)
  
!<description>
  ! This routine releases all temporary memory for the efect correction solver 
  ! from the heap.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure which is to be cleaned up.
  type(t_sptilsNode), intent(INOUT)         :: rsolverNode
!</inputoutput>
  
!</subroutine>
  
    ! Release memory if still associated
    call sptils_doneDataBlockJacobi (rsolverNode)
    call sptils_doneStructureBlockJacobi (rsolverNode)
    
    ! Release the subnode structure
    deallocate(rsolverNode%p_rsubnodeBlockJacobi)
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine sptils_setMatrixBlockJacobi (rsolverNode,Rmatrices)
  
!<description>
  ! This routine is called if the system matrix changes.
  ! The routine calls sptils_setMatrices for the preconditioner
  ! to inform also that one about the change of the matrix pointer.
!</description>
  
!<input>
  ! An array of system matrices which is simply passed to the initialisation 
  ! routine of the preconditioner.
  type(t_ccoptSpaceTimeMatrix), dimension(:), intent(IN) :: Rmatrices
!</input>
  
!<inputoutput>
  ! The t_sptilsNode structure of the solver
  type(t_sptilsNode), intent(INOUT)         :: rsolverNode
!</inputoutput>
  
!</subroutine>

    ! Do we have a preconditioner given? If yes, call the general initialisation
    ! routine to initialise it.
    
  end subroutine

! ***************************************************************************

!<subroutine>
  
  recursive subroutine sptils_initStructureBlockJacobi (rsolverNode,ierror)
  
!<description>
  ! Solver preparation. Perform symbolic factorisation (not of the defect
  ! correcion solver, but of subsolvers). Allocate temporary memory.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure 
  type(t_sptilsNode), intent(INOUT)         :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the SPTILS_ERR_XXXX constants. A value different to 
  ! SPTILS_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  integer, intent(OUT) :: ierror
!</output>
  
!</subroutine>
    
    ! A-priori we have no error...
    ierror = SPTILS_ERR_NOERROR
    
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine sptils_initDataBlockJacobi (rsolverNode, ierror)
  
!<description>
  ! Calls the initData subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with sptils_initData.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure 
  type(t_sptilsNode), intent(INOUT)         :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the SPTILS_ERR_XXXX constants. A value different to 
  ! SPTILS_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  integer, intent(OUT) :: ierror
!</output>

!</subroutine>

    ! A-priori we have no error...
    ierror = SPTILS_ERR_NOERROR

  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine sptils_doneDataBlockJacobi (rsolverNode)
  
!<description>
  ! Calls the doneData subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with sptils_doneData.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure 
  type(t_sptilsNode), intent(INOUT)         :: rsolverNode
!</inputoutput>
  
!</subroutine>

  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine sptils_doneStructureBlockJacobi (rsolverNode)
  
!<description>
  ! Calls the doneStructure subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with sptils_doneStructure.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure 
  type(t_sptilsNode), intent(INOUT)         :: rsolverNode
!</inputoutput>
  
!</subroutine>

  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>
  
  recursive subroutine sptils_precBlockJacobi (rsolverNode,rd)
  
!<description>
  ! Applies the Block Jacobi preconditioner $P \approx A$ to the defect 
  ! vector rd and solves $Pd_{new} = d$.
  ! rd will be overwritten by the preconditioned defect.
  !
  ! The term 'Block Jacobi' means that a spatial preconditioner is applied
  ! to the diagonal block of each time step in the global space-time
  ! system. Usually, this spatial preconditioner is a linear solver
  ! (multigrid or whatever). More precisely, any spatial preconditioner that
  ! can be used with cc_precondDefect is suitable.
  !
  ! The matrix must have been attached to the system before calling
  ! this routine, and the initStructure/initData routines
  ! must have been called to prepare the solver for solving
  ! the problem.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure of the solver
  type(t_sptilsNode), intent(INOUT), target :: rsolverNode

  ! On call to this routine: The defect vector to be preconditioned.
  ! Will be overwritten by the preconditioned defect.
  type(t_spacetimeVector), intent(INOUT)    :: rd
!</inputoutput>
  
!</subroutine>

    ! local variables
    integer :: isubstep, nlmax
    real(DP) :: dtheta,dtstep,dtime
    logical :: bsuccess
    type(t_nonlinearSpatialMatrix) :: rnonlinearSpatialMatrix
    type(t_ccoptSpaceTimeDiscretisation), pointer :: p_rspaceTimeDiscr
    type(t_ccoptSpaceTimeMatrix), pointer :: p_rspaceTimeMatrix
    type(t_vectorBlock) :: rtempVectorD,rtempVectorX1,rtempVectorX2,rtempVectorX3
    type(t_ccspatialPreconditioner), pointer :: p_rspatialPreconditioner
    type(t_cccoreEquationOneLevel), pointer :: p_rlevelInfo
    
    ! DEBUG!!!
    real(DP), dimension(:), pointer :: p_Dx,p_Dd,p_Dsol
    
    call stat_clearTimer (rsolverNode%rtimeSpacePrecond)
    
    ! Get a pointer to the space-time discretisation structure that defines
    ! how to apply the global system matrix.
    p_rspaceTimeDiscr => rsolverNode%rmatrix%p_rspaceTimeDiscr
    p_rspaceTimeMatrix => rsolverNode%rmatrix
    
    ! Get a pointer to our spatial preconditioner
    p_rspatialPreconditioner => &
        rsolverNode%p_rsubnodeBlockJacobi%p_rspatialPreconditioner
        
    ! Get the maximum level. This level corresponds to the above
    ! space-time discretisation
    nlmax = p_rspatialPreconditioner%NLMAX
    
    ! Get preconditioner information about the maximum level
    p_rlevelInfo => p_rspatialPreconditioner%RcoreEquation(nlmax)
    
    dtheta = rsolverNode%p_rproblem%rtimedependence%dtimeStepTheta
    dtstep = p_rspaceTimeDiscr%rtimeDiscr%dtstep

    ! Create temp vectors for X, B and D.
    call lsysbl_createVecBlockByDiscr (p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
        rtempVectorD,.true.)
    call lsysbl_createVecBlockByDiscr (p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
        rtempVectorX1,.true.)
    call lsysbl_createVecBlockByDiscr (p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
        rtempVectorX2,.true.)
    call lsysbl_createVecBlockByDiscr (p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
        rtempVectorX3,.true.)
        
    ! Basic initialisation of our nonlinear matrix in space - for defect correction
    ! and matrix assembly in order to be used in a spatial preconditioner.
    call cc_initNonlinMatrix (rnonlinearSpatialMatrix,rsolverNode%p_rproblem,&
        p_rspaceTimeMatrix%p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
        p_rspaceTimeMatrix%p_rspaceTimeDiscr%p_rlevelInfo%rstaticInfo)

    call cc_preparePrecondMatrixAssembly (rnonlinearSpatialMatrix,&
      nlmax,p_rspatialPreconditioner%nlmin,p_rspatialPreconditioner%nlmax,&
      p_rspatialPreconditioner%rprecSpecials)
      
    ! ----------------------------------------------------------------------
    ! We use a block-Jacobi scheme for preconditioning...
    !
    ! For this purpose, loop through the substeps.
    
    do isubstep = 0,p_rspaceTimeDiscr%NEQtime-1
    
      ! Current time step?
      dtime = &
          rsolverNode%p_rproblem%rtimedependence%dtimeInit + isubstep * dtstep

      if (rsolverNode%ioutputLevel .ge. 1) then
        call output_line ('Space-Time-Block-Jacobi preconditioning of timestep: '//&
            trim(sys_siL(isubstep,10))//&
            ', Time: '//trim(sys_sdL(dtime,10)))
      end if
    
      ! -----
      ! Update the boundary conditions in the preconditioner
      call cc_updatePreconditionerBC (rsolverNode%p_rproblem,&
          rsolverNode%p_rsubnodeBlockJacobi%p_rspatialPreconditioner,dtime)

      ! DEBUG!!!      
      call lsysbl_getbase_double (rtempVectorX2,p_Dx)
      call lsysbl_getbase_double (rtempVectorD,p_Dd)

      ! Read in the RHS/solution/defect vector of the current timestep.
      ! If no solution is specified, we have a linear problem and thus
      ! the content of rtempVector is not relevant; actually it's even
      ! zero by initialisation.
      if (associated(p_rspaceTimeMatrix%p_rsolution)) then
        if (isubstep .gt. 0) then
          call sptivec_getTimestepData (p_rspaceTimeMatrix%p_rsolution, &
              1+isubstep-1, rtempVectorX1)
        end if
        call sptivec_getTimestepData (p_rspaceTimeMatrix%p_rsolution, &
            1+isubstep, rtempVectorX2)
        if (isubstep .lt. p_rspaceTimeDiscr%NEQtime-1) then
          call sptivec_getTimestepData (p_rspaceTimeMatrix%p_rsolution, &
              1+isubstep+1, rtempVectorX3)
        end if
        
        ! DEBUG!!!
        call lsysbl_getbase_double (rtempVectorX2,p_Dsol)
      end if
      call sptivec_getTimestepData (rd, 1+isubstep, rtempVectorD)

      ! Set up the matrix weights for the diagonal matrix
      call cc_initNonlinMatrix (rnonlinearSpatialMatrix,rsolverNode%p_rproblem,&
          p_rspaceTimeMatrix%p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
          p_rspaceTimeMatrix%p_rspaceTimeDiscr%p_rlevelInfo%rstaticInfo)
      call cc_setupMatrixWeights (rsolverNode%p_rproblem,p_rspaceTimeMatrix,dtheta,&
        isubstep,0,rnonlinearSpatialMatrix)
        
      ! Perform preconditioning of the spatial defect with the method provided by the
      ! core equation module.
      call stat_startTimer (rsolverNode%rtimeSpacePrecond)
      call cc_precondSpaceDefect (&
          rsolverNode%p_rsubnodeBlockJacobi%p_rspatialPreconditioner,&
          rnonlinearSpatialMatrix,&
          rtempVectorD,rtempVectorX1,rtempVectorX2,rtempVectorX3,&
          bsuccess,rsolverNode%p_rproblem%rcollection)      
      call stat_stopTimer (rsolverNode%rtimeSpacePrecond)
    
      ! Scale by omega
      call lsysbl_scaleVector (rtempVectorD,rsolverNode%domega)
    
      ! Save back the preconditioned defect.
      call sptivec_setTimestepData (rd, 1+isubstep, rtempVectorD)
      
    end do
    
    call lsysbl_releaseVector (rtempVectorX3)
    call lsysbl_releaseVector (rtempVectorX2)
    call lsysbl_releaseVector (rtempVectorX1)
    call lsysbl_releaseVector (rtempVectorD)
    
  end subroutine

! *****************************************************************************
! Block Gauss-Seidel preconditioner
! *****************************************************************************

!<subroutine>
  
  recursive subroutine sptils_initBlockFBSOR (rproblem,p_rsolverNode,&
      domega,domegaSOR,rspatialPrecond,bpartialUpdate)
  
!<description>
  ! Creates a t_sptilsNode solver structure for the block 
  ! Gauss Seidel preconditioner.
!</description>
  
!<input>
  ! The problem structure that defines the problem.
  ! A pointer to this is saved in the solver node, so the problem structure
  ! must not be released before the solver node is released.
  type(t_problem), target :: rproblem
  
  ! Damping parameter. =1.0: No damping
  real(DP), intent(IN) :: domega

  ! Relaxation parameter. =1.0: Block Gauss-Seidel, =0.0: Block-Jacobi
  real(DP), intent(IN) :: domegaSOR
  
  ! A spatial preconditioner structure that defines how to perform the preconditioning
  ! in each substep.
  ! A pointer to this is noted in p_rsolverNode, so rspatialPrecond
  ! should not be released before the solver is destroyed.
  type(t_ccspatialPreconditioner), intent(IN), target :: rspatialPrecond
  
  ! When updating the solution, update only the primal solution on the forward
  ! and the dual solution on the backward sweep.
  logical, intent(in) :: bpartialUpdate
!</input>
  
!<output>
  ! A pointer to a t_sptilsNode structure. Is set by the routine, any previous
  ! value of the pointer is destroyed.
  type(t_sptilsNode), pointer         :: p_rsolverNode
!</output>
  
!</subroutine>
  
    ! Create a default solver structure
    call sptils_initSolverGeneral(rproblem,p_rsolverNode)
    
    ! Initialise the type of the solver
    p_rsolverNode%calgorithm = SPTILS_ALG_BlockFBSOR
    
    ! Initialise the ability bitfield with the ability of this solver:
    p_rsolverNode%ccapability = 0

    ! Save the relaxation parameter
    p_rsolverNode%domega = domega
    
    ! Allocate a subnode for our solver.
    ! This initialises most of the variables with default values appropriate
    ! to this solver.
    allocate(p_rsolverNode%p_rsubnodeBlockFBSOR)
    
    ! Save the relaxation parameter for the SOR-aproach
    p_rsolverNode%p_rsubnodeBlockFBSOR%domegaSOR = domegaSOR

    ! Attach the preconditioner if given. 
    p_rsolverNode%p_rsubnodeBlockFBSOR%p_rspatialPreconditioner => &
        rspatialPrecond
    
    ! Remember if we only have to do a partial update.    
    p_rsolverNode%p_rsubnodeBlockFBSOR%bpartialUpdate = bpartialUpdate
    
  end subroutine

  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine sptils_doneBlockFBSOR (rsolverNode)
  
!<description>
  ! This routine releases all temporary memory for the solver 
  ! from the heap.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure which is to be cleaned up.
  type(t_sptilsNode), intent(INOUT)         :: rsolverNode
!</inputoutput>
  
!</subroutine>
  
    ! Release memory if still associated
    call sptils_doneDataBlockFBSOR (rsolverNode)
    call sptils_doneStructureBlockFBSOR (rsolverNode)
    
    ! Release the subnode structure
    deallocate(rsolverNode%p_rsubnodeBlockFBSOR)
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine sptils_setMatrixBlockFBSOR (rsolverNode,Rmatrices)
  
!<description>
  ! This routine is called if the system matrix changes.
  ! The routine calls sptils_setMatrices for the preconditioner
  ! to inform also that one about the change of the matrix pointer.
!</description>
  
!<input>
  ! An array of system matrices which is simply passed to the initialisation 
  ! routine of the preconditioner.
  type(t_ccoptSpaceTimeMatrix), dimension(:), intent(IN) :: Rmatrices
!</input>
  
!<inputoutput>
  ! The t_sptilsNode structure of the solver
  type(t_sptilsNode), intent(INOUT)         :: rsolverNode
!</inputoutput>
  
!</subroutine>

    ! Do we have a preconditioner given? If yes, call the general initialisation
    ! routine to initialise it.
    
  end subroutine

! ***************************************************************************

!<subroutine>
  
  recursive subroutine sptils_initStructureBlockFBSOR (rsolverNode,ierror)
  
!<description>
  ! Solver preparation. Perform symbolic factorisation (not of the defect
  ! correcion solver, but of subsolvers). Allocate temporary memory.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure 
  type(t_sptilsNode), intent(INOUT)         :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the SPTILS_ERR_XXXX constants. A value different to 
  ! SPTILS_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  integer, intent(OUT) :: ierror
!</output>
  
!</subroutine>
    
    type(t_ccoptspaceTimeDiscretisation), pointer :: p_rspaceTimeDiscr
    
    p_rspaceTimeDiscr => rsolverNode%rmatrix%p_rspaceTimeDiscr

    ! Allocate memory for temp vectors.
    call sptivec_initVector (rsolverNode%p_rsubnodeBlockFBSOR%rtempVector,&
        p_rspaceTimeDiscr%NEQtime,&
        p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation)

    ! A-priori we have no error...
    ierror = SPTILS_ERR_NOERROR
    
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine sptils_initDataBlockFBSOR (rsolverNode, ierror)
  
!<description>
  ! Calls the initData subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with sptils_initData.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure 
  type(t_sptilsNode), intent(INOUT)         :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the SPTILS_ERR_XXXX constants. A value different to 
  ! SPTILS_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  integer, intent(OUT) :: ierror
!</output>

!</subroutine>

    ! A-priori we have no error...
    ierror = SPTILS_ERR_NOERROR

  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine sptils_doneDataBlockFBSOR (rsolverNode)
  
!<description>
  ! Calls the doneData subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with sptils_doneData.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure 
  type(t_sptilsNode), intent(INOUT)         :: rsolverNode
!</inputoutput>
  
!</subroutine>

  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine sptils_doneStructureBlockFBSOR (rsolverNode)
  
!<description>
  ! Calls the doneStructure subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with sptils_doneStructure.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure 
  type(t_sptilsNode), intent(INOUT)         :: rsolverNode
!</inputoutput>
  
!</subroutine>

    if (rsolverNode%p_rsubnodeBlockFBSOR%rtempVector%NEQ .ne. 0) then
      call sptivec_releaseVector (rsolverNode%p_rsubnodeBlockFBSOR%rtempVector)
    end if

  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>
  
  recursive subroutine sptils_precBlockFBSOR (rsolverNode,rd)
  
!<description>
  ! Applies the Block Jacobi preconditioner $P \approx A$ to the defect 
  ! vector rd and solves $Pd_{new} = d$.
  ! rd will be overwritten by the preconditioned defect.
  !
  ! The term 'Block Jacobi' means that a spatial preconditioner is applied
  ! to the diagonal block of each time step in the global space-time
  ! system. Usually, this spatial preconditioner is a linear solver
  ! (multigrid or whatever). More precisely, any spatial preconditioner that
  ! can be used with cc_precondSpaceDefect is suitable.
  !
  ! The matrix must have been attached to the system before calling
  ! this routine, and the initStructure/initData routines
  ! must have been called to prepare the solver for solving
  ! the problem.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure of the solver
  type(t_sptilsNode), intent(INOUT), target :: rsolverNode

  ! On call to this routine: The defect vector to be preconditioned.
  ! Will be overwritten by the preconditioned defect.
  type(t_spacetimeVector), intent(INOUT)    :: rd
!</inputoutput>
  
!</subroutine>

    ! local variables
    integer :: isubstep,iiteration,nlmax
    real(DP) :: dtheta,dtstep,dtime
    logical :: bsuccess
    type(t_nonlinearSpatialMatrix) :: rnonlinearSpatialMatrix
    type(t_ccoptSpaceTimeDiscretisation), pointer :: p_rspaceTimeDiscr
    type(t_ccoptSpaceTimeMatrix), pointer :: p_rspaceTimeMatrix
    type(t_vectorBlock) :: rtempVectorD1,rtempVectorD2,rtempVectorD3
    type(t_vectorBlock) :: rtempVectorX1,rtempVectorX2,rtempVectorX3
    type(t_vectorBlock) :: rtempVector2X1,rtempVector2X2,rtempVector2X3
    type(t_vectorBlock), dimension(3) :: rtempVectorSol
    type(t_vectorBlock) :: rtempVectorRHS
    type(t_spacetimeVector), pointer :: p_rx
    logical :: bcalcNorm
    real(DP) :: domegaSOR,dres,dresInit
    type(t_ccspatialPreconditioner), pointer :: p_rspatialPreconditioner
    type(t_cccoreEquationOneLevel), pointer :: p_rlevelInfo
    
    ! DEBUG!!!
    real(DP), dimension(:), pointer :: p_Dx1,p_Dx2,p_Dx3,p_Dd1,p_Dd2,p_Dd3,p_Dd,p_Dsol
    
    ! Timer init
    call stat_clearTimer (rsolverNode%rtimeSpacePrecond)
    
    ! Get a pointer to the space-time discretisation structure that defines
    ! how to apply the global system matrix.
    p_rspaceTimeDiscr => rsolverNode%rmatrix%p_rspaceTimeDiscr
    p_rspaceTimeMatrix => rsolverNode%rmatrix
    
    ! Get a pointer to our spatial preconditioner
    p_rspatialPreconditioner => &
        rsolverNode%p_rsubnodeBlockFBSOR%p_rspatialPreconditioner
        
    ! Get the maximum level. This level corresponds to the above
    ! space-time discretisation
    nlmax = p_rspatialPreconditioner%NLMAX
    
    ! Get preconditioner information about the maximum level
    p_rlevelInfo => p_rspatialPreconditioner%RcoreEquation(nlmax)
    
    dtheta = rsolverNode%p_rproblem%rtimedependence%dtimeStepTheta
    dtstep = p_rspaceTimeDiscr%rtimeDiscr%dtstep
    
    domegaSOR = rsolverNode%p_rsubnodeBlockFBSOR%domegaSOR

    ! Create temp vectors
    call lsysbl_createVecBlockByDiscr (p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
        rtempVectorD1,.true.)
    call lsysbl_createVecBlockByDiscr (p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
        rtempVectorD2,.true.)
    call lsysbl_createVecBlockByDiscr (p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
        rtempVectorD3,.true.)
    call lsysbl_createVecBlockByDiscr (p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
        rtempVectorRHS,.true.)

    call lsysbl_createVecBlockByDiscr (p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
        rtempVectorX1,.true.)
    call lsysbl_createVecBlockByDiscr (p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
        rtempVectorX2,.true.)
    call lsysbl_createVecBlockByDiscr (p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
        rtempVectorX3,.true.)

    call lsysbl_createVecBlockByDiscr (p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
        rtempVector2X1,.true.)
    call lsysbl_createVecBlockByDiscr (p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
        rtempVector2X2,.true.)
    call lsysbl_createVecBlockByDiscr (p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
        rtempVector2X3,.true.)

    ! Solution vectors -- for setting up the defect in nonlinear problems.
    ! For the previous (1), current (2) and next (3) time step.
    call lsysbl_createVecBlockByDiscr (p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
        rtempVectorSol(1),.true.)
    call lsysbl_createVecBlockByDiscr (p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
        rtempVectorSol(2),.true.)
    call lsysbl_createVecBlockByDiscr (p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
        rtempVectorSol(3),.true.)
        
    ! DEBUG!!!      
    call lsysbl_getbase_double (rtempVectorX1,p_Dx1)
    call lsysbl_getbase_double (rtempVectorX2,p_Dx2)
    call lsysbl_getbase_double (rtempVectorX3,p_Dx3)
    call lsysbl_getbase_double (rtempVectorD1,p_Dd1)
    call lsysbl_getbase_double (rtempVectorD2,p_Dd2)
    call lsysbl_getbase_double (rtempVectorD3,p_Dd3)
    call lsysbl_getbase_double (rtempVectorRHS,p_Dd)
        
    ! Probably we have to calclate the norm of the residual while calculating...
    bcalcNorm = (rsolverNode%nminIterations .ne. rsolverNode%nmaxIterations) .or.&
                (rsolverNode%ioutputLevel .ge. 2)

    ! Basic initialisation of our nonlinear matrix in space - for defect correction
    ! and matrix assembly in order to be used in a spatial preconditioner.
    call cc_initNonlinMatrix (rnonlinearSpatialMatrix,rsolverNode%p_rproblem,&
        p_rspaceTimeMatrix%p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
        p_rspaceTimeMatrix%p_rspaceTimeDiscr%p_rlevelInfo%rstaticInfo)

    call cc_preparePrecondMatrixAssembly (rnonlinearSpatialMatrix,&
      nlmax,p_rspatialPreconditioner%nlmin,p_rspatialPreconditioner%nlmax,&
      p_rspatialPreconditioner%rprecSpecials)

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
    
    p_rx   => rsolverNode%p_rsubnodeBlockFBSOR%rtempVector
    
    ! Ok, let's start. Initial solution is zero.
    call sptivec_clearVector (p_rx)
    
    ! Norm of the initial residuum
    if (bcalcNorm) then
      dresInit = sptivec_vectorNorm (rd,LINALG_NORML2)
      if (dresInit .eq. 0.0_DP) dresInit = 1.0_DP
      rsolverNode%dinitialDefect = dresInit
      rsolverNode%dlastDefect = 0.0_DP
      dres = dresInit
    end if
    
    do iiteration = 1,rsolverNode%nmaxIterations
    
      ! Probably print the current residuum (of the previous step)
      if (rsolverNode%nminIterations .ne. rsolverNode%nmaxIterations) then
        if (rsolverNode%ioutputLevel .ge. 2) then
          dres = sqrt(dres / real(p_rspaceTimeDiscr%NEQtime,DP))
          call output_line ('Space-Time-Block-FBGS: Iteration '// &
              trim(sys_siL(iiteration-1,10))//',  !!RES!! = '//&
              trim(sys_sdEL(dres,15)) )
        end if
        
        ! Check for convergence
        if ((iiteration .gt. rsolverNode%nminIterations) .and. &
            (iiteration .lt. rsolverNode%nmaxIterations)) then
          if (sptils_testConvergence (rsolverNode, dres)) exit
          if (sptils_testDivergence (rsolverNode, dres)) exit
        end if
      end if

      ! Filter the current solution for boundary conditions in space and time.
      ! rtempVectorRHS is used as temp vector here.
      call tbc_implementInitCondDefect (&
          p_rspaceTimeDiscr,p_rx,rtempVectorRHS)
      call tbc_implementBCdefect (rsolverNode%p_rproblem,&
         p_rspaceTimeDiscr,p_rx,rtempVectorRHS)

      ! -----
      ! Backward in time
      ! -----

      ! Load the last solution vectors for handling the nonlinearity.
      ! rtempVectorSol(2) holds the data for the current timestep.
      ! rtempVectorSol(1) holds the data from the previous and
      ! rtempVectorSol(3) that of the next timestep.
      if (associated(p_rspaceTimeMatrix%p_rsolution)) then
        call sptivec_getTimestepData (p_rspaceTimeMatrix%p_rsolution, &
            p_rspaceTimeDiscr%NEQtime, rtempVectorSol(2))
      end if

      ! Load the RHS and solution of the last timestep.
      call sptivec_getTimestepData (rd, &
          p_rspaceTimeDiscr%NEQtime, rtempVectorD2)
      
      ! Current iterate
      call sptivec_getTimestepData (p_rx, &
          p_rspaceTimeDiscr%NEQtime, rtempVectorX2)

      ! The 2nd temp vector holds the data as well, so to keep a backup of the
      ! three timesteps which change during the iteration.
      call lsysbl_copyVector (rtempVectorX2,rtempVector2X2)

      ! Loop through the substeps we have to update
      do isubstep = p_rspaceTimeDiscr%NEQtime-1,0,-1
      
        ! Current point in time
        dtime = p_rspaceTimeDiscr%rtimeDiscr%dtimeInit + isubstep * dtstep

        if (rsolverNode%ioutputLevel .ge. 1) then
          call output_line ('Space-Time-Block-FBGS preconditioning of timestep: '//&
              trim(sys_siL(isubstep,10))//&
              ', Time: '//trim(sys_sdL(dtime,10)))
        end if
      
        ! Update the boundary conditions in the preconditioner
        call cc_updatePreconditionerBC (rsolverNode%p_rproblem,&
            rsolverNode%p_rsubnodeBlockFBSOR%p_rspatialPreconditioner,dtime)

        ! The RHS which is put into the preconditioner is set up in 
        ! rtempVectorRHS to prevent rtempVectorD2 from getting destroyed.
        call lsysbl_copyVector (rtempVectorD2,rtempVectorRHS)
        
        ! Is this the first timestep or not?
        if (isubstep .gt. 0) then
        
          ! Get the evaluation point for the nonlinearity in the previous (=next) timestep
          ! into rtempVectorSol(1)
          if (associated(p_rspaceTimeMatrix%p_rsolution)) then
            call sptivec_getTimestepData (p_rspaceTimeMatrix%p_rsolution, &
                1+isubstep-1, rtempVectorSol(1))
          end if
        
          ! Read the RHS and solution of the next timestep
          call sptivec_getTimestepData (rd, 1+isubstep-1, rtempVectorD1)
          call sptivec_getTimestepData (p_rx, 1+isubstep-1, rtempVectorX1)

          ! Create d2 = RHS - Mx1 
          call cc_initNonlinMatrix (rnonlinearSpatialMatrix,rsolverNode%p_rproblem,&
              p_rspaceTimeMatrix%p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
              p_rspaceTimeMatrix%p_rspaceTimeDiscr%p_rlevelInfo%rstaticInfo)
          call cc_setupMatrixWeights (rsolverNode%p_rproblem,p_rspaceTimeMatrix,dtheta,&
            isubstep,-1,rnonlinearSpatialMatrix)
          
          call cc_assembleDefect (rnonlinearSpatialMatrix,rtempVectorX1,&
            rtempVectorRHS,1.0_DP,rtempVectorSol(1),rtempVectorSol(2),rtempVectorSol(3))
          
        end if

        ! Is this the last timestep or not?
        if (isubstep .lt. p_rspaceTimeDiscr%NEQtime-1) then
          
          ! Create d2 = RHS - omega*Ml3 - (1-omega)*Ml3_old
          call cc_initNonlinMatrix (rnonlinearSpatialMatrix,rsolverNode%p_rproblem,&
              p_rspaceTimeMatrix%p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
              p_rspaceTimeMatrix%p_rspaceTimeDiscr%p_rlevelInfo%rstaticInfo)
          call cc_setupMatrixWeights (rsolverNode%p_rproblem,p_rspaceTimeMatrix,dtheta,&
            isubstep,1,rnonlinearSpatialMatrix)
            
          call cc_assembleDefect (rnonlinearSpatialMatrix,rtempVectorX3,&
            rtempVectorRHS,domegaSOR,rtempVectorSol(1),rtempVectorSol(2),rtempVectorSol(3))

          call cc_assembleDefect (rnonlinearSpatialMatrix,rtempVector2X3,&
            rtempVectorRHS,1.0_DP-domegaSOR,rtempVectorSol(1),rtempVectorSol(2),rtempVectorSol(3))
          
        end if

        ! Set up the matrix weights for the diagonal matrix
        call cc_initNonlinMatrix (rnonlinearSpatialMatrix,rsolverNode%p_rproblem,&
            p_rspaceTimeMatrix%p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
            p_rspaceTimeMatrix%p_rspaceTimeDiscr%p_rlevelInfo%rstaticInfo)
        call cc_setupMatrixWeights (rsolverNode%p_rproblem,p_rspaceTimeMatrix,dtheta,&
          isubstep,0,rnonlinearSpatialMatrix)
          
        ! Create d2 = RHS - A(solution) X2
        call cc_assembleDefect (rnonlinearSpatialMatrix,rtempVectorX2,rtempVectorRHS,&
            1.0_DP,rtempVectorSol(1),rtempVectorSol(2),rtempVectorSol(3))
            
        ! Filter the defect for BC's and initial conditions if necessary
        if (isubstep .eq. 0) then
          call tbc_implementInitCondDefSingle (p_rspaceTimeDiscr, rtempVectorRHS)
        else if (isubstep .eq. p_rspaceTimeDiscr%NEQtime-1) then
          call tbc_implementTermCondDefSingle (p_rspaceTimeDiscr, rtempVectorRHS)
        end if

        call tbc_implementSpatialBCdefect (&
            rsolverNode%p_rproblem,dtime,p_rspaceTimeDiscr,rtempVectorRHS)
        
        ! Ok, we have the local defect.
        !
        ! Perform preconditioning of the spatial defect with the method 
        ! provided by the core equation module.
        call stat_startTimer (rsolverNode%rtimeSpacePrecond)
        call cc_precondSpaceDefect (&
            rsolverNode%p_rsubnodeBlockFBSOR%p_rspatialPreconditioner,&
            rnonlinearSpatialMatrix,rtempVectorRHS,&
            rtempVectorSol(1),rtempVectorSol(2),rtempVectorSol(3),&
            bsuccess,rsolverNode%p_rproblem%rcollection)      
        call stat_stopTimer (rsolverNode%rtimeSpacePrecond)
        
        ! Add that defect to the current solution -- damped by domega.
        if (.not. rsolverNode%p_rsubnodeBlockFBSOR%bpartialUpdate) then
          call lsysbl_vectorLinearComb (rtempVectorRHS,rtempVectorX2,&
              rsolverNode%domega,1.0_DP)
        else
          call lsyssc_vectorLinearComb (rtempVectorRHS%RvectorBlock(4),&
              rtempVectorX2%RvectorBlock(4),rsolverNode%domega,1.0_DP)
          call lsyssc_vectorLinearComb (rtempVectorRHS%RvectorBlock(5),&
              rtempVectorX2%RvectorBlock(5),rsolverNode%domega,1.0_DP)
          call lsyssc_vectorLinearComb (rtempVectorRHS%RvectorBlock(6),&
              rtempVectorX2%RvectorBlock(6),rsolverNode%domega,1.0_DP)
        end if
      
        ! Save the new solution.
        call sptivec_setTimestepData (p_rx, 1+isubstep, rtempVectorX2)
        
        ! Shift the RHS/solution vectors: 1 -> 2 -> 3
        call lsysbl_copyVector (rtempVectorD2,rtempVectorD3)
        call lsysbl_copyVector (rtempVectorD1,rtempVectorD2)

        call lsysbl_copyVector (rtempVectorX2,rtempVectorX3)
        call lsysbl_copyVector (rtempVectorX1,rtempVectorX2)

        ! The 2nd temp vector holds the data as well, so to keep a backup of the
        ! three timesteps which change during the iteration.
        call lsysbl_copyVector (rtempVector2X2,rtempVector2X3)
        call lsysbl_copyVector (rtempVectorX2,rtempVector2X2)

        if (associated(p_rspaceTimeMatrix%p_rsolution)) then
          call lsysbl_copyVector (rtempVectorSol(2),rtempVectorSol(3))
          call lsysbl_copyVector (rtempVectorSol(1),rtempVectorSol(2))
        end if

      end do
      
      ! -----
      ! Forward in time
      ! -----
      
      ! rtempVectorSol(1) is undefined, but we don't need it.
      ! rtempVectorSol(2) holds the data of the 0th timestep, 
      ! rtempVectorSol(3) of the 1st one.
      
      if (associated(p_rspaceTimeMatrix%p_rsolution)) then
        call sptivec_getTimestepData (p_rspaceTimeMatrix%p_rsolution, &
            1, rtempVectorSol(2))
        ! rtempVectorSol(3) is set later.
      end if
      
      ! Norm of the residuum
      dres = 0.0_DP

      ! Load the RHS and solution of the 0th timestep.
      call sptivec_getTimestepData (rd, 1+0, rtempVectorD2)
      
      ! Current iterate
      call sptivec_getTimestepData (p_rx, 1+0, rtempVectorX2)

      ! The 2nd temp vector holds the data as well, so to keep a backup of the
      ! three timesteps which change during the iteration.
      call lsysbl_copyVector (rtempVectorX2,rtempVector2X2)

      ! Loop through the substeps we have to update
      do isubstep = 0,p_rspaceTimeDiscr%NEQtime-1
      
        ! Current point in time
        dtime = p_rspaceTimeDiscr%rtimeDiscr%dtimeInit + isubstep * dtstep

        if (rsolverNode%ioutputLevel .ge. 1) then
          call output_line ('Space-Time-Block-FBGS preconditioning of timestep: '//&
              trim(sys_siL(isubstep,10))//&
              ', Time: '//trim(sys_sdL(dtime,10)))
        end if
      
        ! Update the boundary conditions in the preconditioner
        call cc_updatePreconditionerBC (rsolverNode%p_rproblem,&
            rsolverNode%p_rsubnodeBlockFBSOR%p_rspatialPreconditioner,dtime)

        ! The RHS which is put into the preconditioner is set up in 
        ! rtempVectorRHS to prevent rtempVectorD2 from getting destroyed.
        call lsysbl_copyVector (rtempVectorD2,rtempVectorRHS)
        
        ! Is this the first timestep or not?
        if (isubstep .gt. 0) then
        
          ! Create d2 = RHS - omega*Mx1 - (1-omega)*Mx1_old
          call cc_initNonlinMatrix (rnonlinearSpatialMatrix,rsolverNode%p_rproblem,&
              p_rspaceTimeMatrix%p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
              p_rspaceTimeMatrix%p_rspaceTimeDiscr%p_rlevelInfo%rstaticInfo)
          call cc_setupMatrixWeights (rsolverNode%p_rproblem,p_rspaceTimeMatrix,dtheta,&
            isubstep,-1,rnonlinearSpatialMatrix)
          
          call cc_assembleDefect (rnonlinearSpatialMatrix,rtempVectorX1,&
            rtempVectorRHS,domegaSOR,rtempVectorSol(1),rtempVectorSol(2),rtempVectorSol(3))

          call cc_assembleDefect (rnonlinearSpatialMatrix,rtempVector2X1,&
            rtempVectorRHS,1.0_DP-domegaSOR,rtempVectorSol(1),rtempVectorSol(2),rtempVectorSol(3))
          
        end if

        ! Is this the last timestep or not?
        if (isubstep .lt. p_rspaceTimeDiscr%NEQtime-1) then
          
          ! Get the evaluation point for the nonlinearity in the next timestep
          ! into rtempVectorSol(3).
          if (associated(p_rspaceTimeMatrix%p_rsolution)) then
            call sptivec_getTimestepData (p_rspaceTimeMatrix%p_rsolution, &
                1+isubstep+1, rtempVectorSol(3))
          end if

          ! Read the RHS and solution of the next timestep
          call sptivec_getTimestepData (rd, 1+isubstep+1, rtempVectorD3)
          call sptivec_getTimestepData (p_rx, 1+isubstep+1, rtempVectorX3)
        
          ! Create d2 = RHS - Ml3 
          call cc_initNonlinMatrix (rnonlinearSpatialMatrix,rsolverNode%p_rproblem,&
              p_rspaceTimeMatrix%p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
              p_rspaceTimeMatrix%p_rspaceTimeDiscr%p_rlevelInfo%rstaticInfo)
          call cc_setupMatrixWeights (rsolverNode%p_rproblem,p_rspaceTimeMatrix,dtheta,&
            isubstep,1,rnonlinearSpatialMatrix)
          
          call cc_assembleDefect (rnonlinearSpatialMatrix,rtempVectorX3,&
            rtempVectorRHS,1.0_DP,rtempVectorSol(1),rtempVectorSol(2),rtempVectorSol(3))
          
        end if

        ! Set up the matrix weights for the diagonal matrix
        call cc_initNonlinMatrix (rnonlinearSpatialMatrix,rsolverNode%p_rproblem,&
            p_rspaceTimeMatrix%p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
            p_rspaceTimeMatrix%p_rspaceTimeDiscr%p_rlevelInfo%rstaticInfo)
        call cc_setupMatrixWeights (rsolverNode%p_rproblem,p_rspaceTimeMatrix,dtheta,&
          isubstep,0,rnonlinearSpatialMatrix)
          
        ! Create d2 = RHS - A(solution) X2
        call cc_assembleDefect (rnonlinearSpatialMatrix,rtempVectorX2,rtempVectorRHS,&
            1.0_DP,rtempVectorSol(1),rtempVectorSol(2),rtempVectorSol(3))
            
        ! Filter the defect for BC's and initial conditions if necessary
        if (isubstep .eq. 0) then
          call tbc_implementInitCondDefSingle (p_rspaceTimeDiscr, rtempVectorRHS)
        else if (isubstep .eq. p_rspaceTimeDiscr%NEQtime-1) then
          call tbc_implementTermCondDefSingle (p_rspaceTimeDiscr, rtempVectorRHS)
        end if

        call tbc_implementSpatialBCdefect (&
            rsolverNode%p_rproblem,dtime,p_rspaceTimeDiscr,rtempVectorRHS)
        
        ! Ok, we have the local defect.
        ! Sum up the norm to the norm of the global vector.
        if (bcalcNorm) &
          dres = dres + lsysbl_vectorNorm (rtempVectorRHS,LINALG_NORML2)**2
        
        ! Perform preconditioning of the spatial defect with the method 
        ! provided by the core equation module.
        call stat_startTimer (rsolverNode%rtimeSpacePrecond)
        call cc_precondSpaceDefect (&
            rsolverNode%p_rsubnodeBlockFBSOR%p_rspatialPreconditioner,&
            rnonlinearSpatialMatrix,rtempVectorRHS,&
            rtempVectorSol(1),rtempVectorSol(2),rtempVectorSol(3),&
            bsuccess,rsolverNode%p_rproblem%rcollection)      
        call stat_stopTimer (rsolverNode%rtimeSpacePrecond)
      
        ! Add that defect to the current solution -- damped by domega.
        if (.not. rsolverNode%p_rsubnodeBlockFBSOR%bpartialUpdate) then
          call lsysbl_vectorLinearComb (rtempVectorRHS,rtempVectorX2,&
              rsolverNode%domega,1.0_DP)
        else
          call lsyssc_vectorLinearComb (rtempVectorRHS%RvectorBlock(1),&
              rtempVectorX2%RvectorBlock(1),rsolverNode%domega,1.0_DP)
          call lsyssc_vectorLinearComb (rtempVectorRHS%RvectorBlock(2),&
              rtempVectorX2%RvectorBlock(2),rsolverNode%domega,1.0_DP)
          call lsyssc_vectorLinearComb (rtempVectorRHS%RvectorBlock(3),&
              rtempVectorX2%RvectorBlock(3),rsolverNode%domega,1.0_DP)
        end if
      
        ! Save the new solution.
        call sptivec_setTimestepData (p_rx, 1+isubstep, rtempVectorX2)
        
        ! Shift the RHS/solution vectors: 1 <- 2 <- 3
        call lsysbl_copyVector (rtempVectorD2,rtempVectorD1)
        call lsysbl_copyVector (rtempVectorD3,rtempVectorD2)

        call lsysbl_copyVector (rtempVectorX2,rtempVectorX1)
        call lsysbl_copyVector (rtempVectorX3,rtempVectorX2)

        ! The 2nd temp vector holds the data as well, so to keep a backup of the
        ! three timesteps which change during the iteration.
        call lsysbl_copyVector (rtempVector2X2,rtempVector2X1)
        call lsysbl_copyVector (rtempVectorX2,rtempVector2X2)

        if (associated(p_rspaceTimeMatrix%p_rsolution)) then
          call lsysbl_copyVector (rtempVectorSol(2),rtempVectorSol(1))
          call lsysbl_copyVector (rtempVectorSol(3),rtempVectorSol(2))
        end if

      end do
      
    end do ! iiteration
    
    ! Overwrite the rd by our solution.
    call sptivec_copyVector (p_rx,rd)
    
    ! Release memory, finish.
    call lsysbl_releaseVector (rtempVectorRHS)
    call lsysbl_releaseVector (rtempVectorSol(3))
    call lsysbl_releaseVector (rtempVectorSol(2))
    call lsysbl_releaseVector (rtempVectorSol(1))
    call lsysbl_releaseVector (rtempVectorD3)
    call lsysbl_releaseVector (rtempVectorD2)
    call lsysbl_releaseVector (rtempVectorD1)
    call lsysbl_releaseVector (rtempVector2X3)
    call lsysbl_releaseVector (rtempVector2X2)
    call lsysbl_releaseVector (rtempVector2X1)
    call lsysbl_releaseVector (rtempVectorX3)
    call lsysbl_releaseVector (rtempVectorX2)
    call lsysbl_releaseVector (rtempVectorX1)
    
  end subroutine

! *****************************************************************************
! Block FBGS preconditioner
! *****************************************************************************

!<subroutine>
  
  recursive subroutine sptils_initBlockFBGS (rproblem,p_rsolverNode,&
      domega,domegaSOR,rspatialPrecond)
  
!<description>
  ! Creates a t_sptilsNode solver structure for the block 
  ! FBGS preconditioner.
!</description>
  
!<input>
  ! The problem structure that defines the problem.
  ! A pointer to this is saved in the solver node, so the problem structure
  ! must not be released before the solver node is released.
  type(t_problem), target :: rproblem
  
  ! Relaxation parameter. =1.0: Full block FBGS.
  real(DP), intent(IN) :: domega

  ! Relaxation parameter for the time dependence. =1.0=Standard.
  real(DP), intent(IN) :: domegaSOR
  
  ! A spatial preconditioner structure that defines how to perform the preconditioning
  ! in each substep.
  ! A pointer to this is noted in p_rsolverNode, so rspatialPrecond
  ! should not be released before the solver is destroyed.
  type(t_ccspatialPreconditioner), intent(IN), target :: rspatialPrecond
!</input>
  
!<output>
  ! A pointer to a t_sptilsNode structure. Is set by the routine, any previous
  ! value of the pointer is destroyed.
  type(t_sptilsNode), pointer         :: p_rsolverNode
!</output>
  
!</subroutine>
  
    ! Create a default solver structure
    call sptils_initSolverGeneral(rproblem,p_rsolverNode)
    
    ! Initialise the type of the solver
    p_rsolverNode%calgorithm = SPTILS_ALG_BLOCKFBGS
    
    ! Initialise the ability bitfield with the ability of this solver:
    p_rsolverNode%ccapability = 0

    ! Save the relaxation parameter
    p_rsolverNode%domega = domega
    
    ! Allocate a subnode for our solver.
    ! This initialises most of the variables with default values appropriate
    ! to this solver.
    allocate(p_rsolverNode%p_rsubnodeBlockFBGS)
    
    ! Save the relaxation parameter for the SOR-aproach
    p_rsolverNode%p_rsubnodeBlockFBGS%domegaSOR = domegaSOR

    ! Attach the preconditioner if given. 
    p_rsolverNode%p_rsubnodeBlockFBGS%p_rspatialPreconditioner => &
        rspatialPrecond
    
  end subroutine

  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine sptils_doneBlockFBGS (rsolverNode)
  
!<description>
  ! This routine releases all temporary memory for the solver 
  ! from the heap.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure which is to be cleaned up.
  type(t_sptilsNode), intent(INOUT)         :: rsolverNode
!</inputoutput>
  
!</subroutine>
  
    ! Release memory if still associated
    call sptils_doneDataBlockFBGS (rsolverNode)
    call sptils_doneStructureBlockFBGS (rsolverNode)
    
    ! Release the subnode structure
    deallocate(rsolverNode%p_rsubnodeBlockFBGS)
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine sptils_setMatrixBlockFBGS (rsolverNode,Rmatrices)
  
!<description>
  ! This routine is called if the system matrix changes.
  ! The routine calls sptils_setMatrices for the preconditioner
  ! to inform also that one about the change of the matrix pointer.
!</description>
  
!<input>
  ! An array of system matrices which is simply passed to the initialisation 
  ! routine of the preconditioner.
  type(t_ccoptSpaceTimeMatrix), dimension(:), intent(IN) :: Rmatrices
!</input>
  
!<inputoutput>
  ! The t_sptilsNode structure of the solver
  type(t_sptilsNode), intent(INOUT)         :: rsolverNode
!</inputoutput>
  
!</subroutine>

    ! Do we have a preconditioner given? If yes, call the general initialisation
    ! routine to initialise it.
    
  end subroutine

! ***************************************************************************

!<subroutine>
  
  recursive subroutine sptils_initStructureBlockFBGS (rsolverNode,ierror)
  
!<description>
  ! Solver preparation. Perform symbolic factorisation (not of the defect
  ! correcion solver, but of subsolvers). Allocate temporary memory.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure 
  type(t_sptilsNode), intent(INOUT), target :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the SPTILS_ERR_XXXX constants. A value different to 
  ! SPTILS_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  integer, intent(OUT) :: ierror
!</output>
  
!</subroutine>

    type(t_ccoptspaceTimeDiscretisation), pointer :: p_rspaceTimeDiscr
    
    p_rspaceTimeDiscr => rsolverNode%rmatrix%p_rspaceTimeDiscr

    ! Allocate memory for a temp vector.
    call sptivec_initVector (rsolverNode%p_rsubnodeBlockFBGS%rtempVector,&
        p_rspaceTimeDiscr%NEQtime,&
        p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation)

    ! A-priori we have no error...
    ierror = SPTILS_ERR_NOERROR
    
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine sptils_initDataBlockFBGS (rsolverNode, ierror)
  
!<description>
  ! Calls the initData subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with sptils_initData.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure 
  type(t_sptilsNode), intent(INOUT)         :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the SPTILS_ERR_XXXX constants. A value different to 
  ! SPTILS_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  integer, intent(OUT) :: ierror
!</output>

!</subroutine>

    ! A-priori we have no error...
    ierror = SPTILS_ERR_NOERROR

  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine sptils_doneDataBlockFBGS (rsolverNode)
  
!<description>
  ! Calls the doneData subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with sptils_doneData.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure 
  type(t_sptilsNode), intent(INOUT)         :: rsolverNode
!</inputoutput>
  
!</subroutine>

  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine sptils_doneStructureBlockFBGS (rsolverNode)
  
!<description>
  ! Calls the doneStructure subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with sptils_doneStructure.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure 
  type(t_sptilsNode), intent(INOUT)         :: rsolverNode
!</inputoutput>
  
!</subroutine>

    if (rsolverNode%p_rsubnodeBlockFBGS%rtempVector%NEQtime .ne. 0) &
      call sptivec_releaseVector(rsolverNode%p_rsubnodeBlockFBGS%rtempVector)

  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>
  
  recursive subroutine sptils_precBlockFBGS (rsolverNode,rd)
  
!<description>
  ! Applies the Block Jacobi preconditioner $P \approx A$ to the defect 
  ! vector rd and solves $Pd_{new} = d$.
  ! rd will be overwritten by the preconditioned defect.
  !
  ! The term 'Block Jacobi' means that a spatial preconditioner is applied
  ! to the diagonal block of each time step in the global space-time
  ! system. Usually, this spatial preconditioner is a linear solver
  ! (multigrid or whatever). More precisely, any spatial preconditioner that
  ! can be used with cc_precondSpaceDefect is suitable.
  !
  ! The matrix must have been attached to the system before calling
  ! this routine, and the initStructure/initData routines
  ! must have been called to prepare the solver for solving
  ! the problem.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure of the solver
  type(t_sptilsNode), intent(INOUT), target :: rsolverNode

  ! On call to this routine: The defect vector to be preconditioned.
  ! Will be overwritten by the preconditioned defect.
  type(t_spacetimeVector), intent(INOUT)    :: rd
!</inputoutput>
  
!</subroutine>

    ! local variables
    integer :: isubstep,iiteration,nlmax
    real(DP) :: dtheta,dtstep,dtime
    logical :: bsuccess
    type(t_nonlinearSpatialMatrix) :: rnonlinearSpatialMatrix
    type(t_ccoptSpaceTimeDiscretisation), pointer :: p_rspaceTimeDiscr
    type(t_ccoptSpaceTimeMatrix), pointer :: p_rspaceTimeMatrix
    type(t_vectorBlock) :: rtempVectorD1,rtempVectorD2,rtempVectorD3
    type(t_vectorBlock) :: rtempVectorX1,rtempVectorX2,rtempVectorX3
    type(t_vectorBlock), dimension(3) :: rtempVectorSol
    type(t_vectorBlock) :: rtempVectorRHS
    type(t_spacetimeVector), pointer :: p_rx
    logical :: bcalcNorm
    real(DP) :: domegaSOR,dres,dresInit
    type(t_ccspatialPreconditioner), pointer :: p_rspatialPreconditioner
    type(t_cccoreEquationOneLevel), pointer :: p_rlevelInfo
    
    ! DEBUG!!!
    real(DP), dimension(:), pointer :: p_Dx1,p_Dx2,p_Dx3,p_Dd1,p_Dd2,p_Dd3,p_Dd,p_Dsol
    
    ! Timer init
    call stat_clearTimer (rsolverNode%rtimeSpacePrecond)
    
    ! Get a pointer to the space-time discretisation structure that defines
    ! how to apply the global system matrix.
    p_rspaceTimeDiscr => rsolverNode%rmatrix%p_rspaceTimeDiscr
    p_rspaceTimeMatrix => rsolverNode%rmatrix
    
    ! Get a pointer to our spatial preconditioner
    p_rspatialPreconditioner => &
        rsolverNode%p_rsubnodeBlockFBGS%p_rspatialPreconditioner
        
    ! Get the maximum level. This level corresponds to the above
    ! space-time discretisation
    nlmax = p_rspatialPreconditioner%NLMAX
    
    ! Get preconditioner information about the maximum level
    p_rlevelInfo => p_rspatialPreconditioner%RcoreEquation(nlmax)
    
    dtheta = rsolverNode%p_rproblem%rtimedependence%dtimeStepTheta
    dtstep = p_rspaceTimeDiscr%rtimeDiscr%dtstep
    
    domegaSOR = rsolverNode%p_rsubnodeBlockFBGS%domegaSOR

    ! Create temp vectors
    call lsysbl_createVecBlockByDiscr (p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
        rtempVectorD1,.true.)
    call lsysbl_createVecBlockByDiscr (p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
        rtempVectorD2,.true.)
    call lsysbl_createVecBlockByDiscr (p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
        rtempVectorD3,.true.)
    call lsysbl_createVecBlockByDiscr (p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
        rtempVectorRHS,.true.)

    call lsysbl_createVecBlockByDiscr (p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
        rtempVectorX1,.true.)
    call lsysbl_createVecBlockByDiscr (p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
        rtempVectorX2,.true.)
    call lsysbl_createVecBlockByDiscr (p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
        rtempVectorX3,.true.)

    ! Solution vectors -- for setting up the defect in nonlinear problems.
    ! For the previous (1), current (2) and next (3) time step.
    call lsysbl_createVecBlockByDiscr (p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
        rtempVectorSol(1),.true.)
    call lsysbl_createVecBlockByDiscr (p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
        rtempVectorSol(2),.true.)
    call lsysbl_createVecBlockByDiscr (p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
        rtempVectorSol(3),.true.)
        
    ! DEBUG!!!      
    call lsysbl_getbase_double (rtempVectorX1,p_Dx1)
    call lsysbl_getbase_double (rtempVectorX2,p_Dx2)
    call lsysbl_getbase_double (rtempVectorX3,p_Dx3)
    call lsysbl_getbase_double (rtempVectorD1,p_Dd1)
    call lsysbl_getbase_double (rtempVectorD2,p_Dd2)
    call lsysbl_getbase_double (rtempVectorD3,p_Dd3)
    call lsysbl_getbase_double (rtempVectorRHS,p_Dd)
        
    ! Probably we have to calclate the norm of the residual while calculating...
    bcalcNorm = (rsolverNode%nminIterations .ne. rsolverNode%nmaxIterations) .or.&
                (rsolverNode%ioutputLevel .ge. 2)

    ! Basic initialisation of our nonlinear matrix in space - for defect correction
    ! and matrix assembly in order to be used in a spatial preconditioner.
    call cc_initNonlinMatrix (rnonlinearSpatialMatrix,rsolverNode%p_rproblem,&
        p_rspaceTimeMatrix%p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
        p_rspaceTimeMatrix%p_rspaceTimeDiscr%p_rlevelInfo%rstaticInfo)

    call cc_preparePrecondMatrixAssembly (rnonlinearSpatialMatrix,&
      nlmax,p_rspatialPreconditioner%nlmin,p_rspatialPreconditioner%nlmax,&
      p_rspatialPreconditioner%rprecSpecials)

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
    call sptivec_clearVector (p_rx)
    
    ! Norm of the initial residuum
    if (bcalcNorm) then
      dresInit = sptivec_vectorNorm (rd,LINALG_NORML2)
      if (dresInit .eq. 0.0_DP) dresInit = 1.0_DP
      rsolverNode%dinitialDefect = dresInit
      rsolverNode%dlastDefect = 0.0_DP
      dres = dresInit
    end if
    


!      ! DEBUG!!!
!      !...
!      CALL ccopt_precSpaceTimeVanca (rsolverNode%p_rproblem,p_rspaceTimeMatrix,p_rx,rd,&
!          rsolverNode%domega,0,1)
!      CALL sptivec_copyVector (p_rx,rd)
!      CALL sptivec_scaleVector (rd,0.7_DP)

!      CALL lsysbl_releaseVector (rtempVectorRHS)
!      CALL lsysbl_releaseVector (rtempVectorSol(3))
!      CALL lsysbl_releaseVector (rtempVectorSol(2))
!      CALL lsysbl_releaseVector (rtempVectorSol(1))
!      CALL lsysbl_releaseVector (rtempVectorD3)
!      CALL lsysbl_releaseVector (rtempVectorD2)
!      CALL lsysbl_releaseVector (rtempVectorD1)
!      CALL lsysbl_releaseVector (rtempVectorX3)
!      CALL lsysbl_releaseVector (rtempVectorX2)
!      CALL lsysbl_releaseVector (rtempVectorX1)
!      RETURN



    do iiteration = 1,rsolverNode%nmaxIterations

      ! Probably print the current residuum (of the previous step)
      if (rsolverNode%nminIterations .ne. rsolverNode%nmaxIterations) then
        if (rsolverNode%ioutputLevel .ge. 2) then
          dres = sqrt(dres / real(p_rspaceTimeDiscr%NEQtime,DP))
          call output_line ('Space-Time-Block-FBGS: Iteration '// &
              trim(sys_siL(iiteration-1,10))//',  !!RES!! = '//&
              trim(sys_sdEL(dres,15)) )
        end if
        
        ! Check for convergence
        if ((iiteration .gt. rsolverNode%nminIterations) .and. &
            (iiteration .lt. rsolverNode%nmaxIterations)) then
          if (sptils_testConvergence (rsolverNode, dres)) exit
          if (sptils_testDivergence (rsolverNode, dres)) exit
        end if
      end if

      ! Filter the current solution for boundary conditions in space and time.
      ! rtempVectorRHS is used as temp vector here.
      call tbc_implementInitCondDefect (&
          p_rspaceTimeDiscr,p_rx,rtempVectorRHS)
      call tbc_implementBCdefect (rsolverNode%p_rproblem,&
         p_rspaceTimeDiscr,p_rx,rtempVectorRHS)

      ! -----
      ! Backward in time
      ! -----

      ! Load the last solution vectors for handling the nonlinearity.
      ! rtempVectorSol(2) holds the data for the current timestep.
      ! rtempVectorSol(1) holds the data from the previous and
      ! rtempVectorSol(3) that of the next timestep.
      if (associated(p_rspaceTimeMatrix%p_rsolution)) then
        call sptivec_getTimestepData (p_rspaceTimeMatrix%p_rsolution, &
            p_rspaceTimeDiscr%NEQtime, rtempVectorSol(2))
      end if

      ! Load the RHS and solution of the last timestep.
      call sptivec_getTimestepData (rd, &
          p_rspaceTimeDiscr%NEQtime, rtempVectorD2)
      
      ! Current iterate
      call sptivec_getTimestepData (p_rx, &
          p_rspaceTimeDiscr%NEQtime, rtempVectorX2)

      ! Loop through the substeps we have to update
      do isubstep = p_rspaceTimeDiscr%NEQtime-1,0,-1
      
        ! Current point in time
        dtime = p_rspaceTimeDiscr%rtimeDiscr%dtimeInit + isubstep * dtstep

        if (rsolverNode%ioutputLevel .ge. 1) then
          call output_line ('Space-Time-Block-FBGS preconditioning of timestep: '//&
              trim(sys_siL(isubstep,10))//&
              ', Time: '//trim(sys_sdL(dtime,10)))
        end if
      
        ! Update the boundary conditions in the preconditioner
        call cc_updatePreconditionerBC (rsolverNode%p_rproblem,&
            rsolverNode%p_rsubnodeBlockFBGS%p_rspatialPreconditioner,dtime)

        ! The RHS which is put into the preconditioner is set up in 
        ! rtempVectorRHS to prevent rtempVectorD2 from getting destroyed.
        call lsysbl_copyVector (rtempVectorD2,rtempVectorRHS)
        
        ! Is this the first timestep or not?
        if (isubstep .gt. 0) then
        
          ! Get the evaluation point for the nonlinearity in the previous (=next) timestep
          ! into rtempVectorSol(1)
          if (associated(p_rspaceTimeMatrix%p_rsolution)) then
            call sptivec_getTimestepData (p_rspaceTimeMatrix%p_rsolution, &
                1+isubstep-1, rtempVectorSol(1))
          end if
        
          ! Read the RHS and solution of the next timestep
          call sptivec_getTimestepData (rd, 1+isubstep-1, rtempVectorD1)
          call sptivec_getTimestepData (p_rx, 1+isubstep-1, rtempVectorX1)

          ! Create d2 = RHS - Mx1 
          call cc_initNonlinMatrix (rnonlinearSpatialMatrix,rsolverNode%p_rproblem,&
              p_rspaceTimeMatrix%p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
              p_rspaceTimeMatrix%p_rspaceTimeDiscr%p_rlevelInfo%rstaticInfo)
          call cc_setupMatrixWeights (rsolverNode%p_rproblem,p_rspaceTimeMatrix,dtheta,&
            isubstep,-1,rnonlinearSpatialMatrix)
          
          call cc_assembleDefect (rnonlinearSpatialMatrix,rtempVectorX1,&
            rtempVectorRHS,domegaSOR,rtempVectorSol(1),rtempVectorSol(2),rtempVectorSol(3))
          
        end if

        ! Is this the last timestep or not?
        if (isubstep .lt. p_rspaceTimeDiscr%NEQtime-1) then
          
          ! Create d2 = RHS - Ml3 
          call cc_initNonlinMatrix (rnonlinearSpatialMatrix,rsolverNode%p_rproblem,&
              p_rspaceTimeMatrix%p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
              p_rspaceTimeMatrix%p_rspaceTimeDiscr%p_rlevelInfo%rstaticInfo)
          call cc_setupMatrixWeights (rsolverNode%p_rproblem,p_rspaceTimeMatrix,dtheta,&
            isubstep,1,rnonlinearSpatialMatrix)
            
          call cc_assembleDefect (rnonlinearSpatialMatrix,rtempVectorX3,&
            rtempVectorRHS,domegaSOR,rtempVectorSol(1),rtempVectorSol(2),rtempVectorSol(3))
          
        end if

        ! Set up the matrix weights for the diagonal matrix
        call cc_initNonlinMatrix (rnonlinearSpatialMatrix,rsolverNode%p_rproblem,&
            p_rspaceTimeMatrix%p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
            p_rspaceTimeMatrix%p_rspaceTimeDiscr%p_rlevelInfo%rstaticInfo)
        call cc_setupMatrixWeights (rsolverNode%p_rproblem,p_rspaceTimeMatrix,dtheta,&
          isubstep,0,rnonlinearSpatialMatrix)
          
        ! Create d2 = RHS - A(solution) X2
        call cc_assembleDefect (rnonlinearSpatialMatrix,rtempVectorX2,rtempVectorRHS,&
            1.0_DP,rtempVectorSol(1),rtempVectorSol(2),rtempVectorSol(2))
            
        ! Filter the defect for BC's and initial conditions if necessary
        if (isubstep .eq. 0) then
          call tbc_implementInitCondDefSingle (p_rspaceTimeDiscr, rtempVectorRHS)
        else if (isubstep .eq. p_rspaceTimeDiscr%NEQtime-1) then
          call tbc_implementTermCondDefSingle (p_rspaceTimeDiscr, rtempVectorRHS)
        end if

        call tbc_implementSpatialBCdefect (&
            rsolverNode%p_rproblem,dtime,p_rspaceTimeDiscr,rtempVectorRHS)
        
        ! Ok, we have the local defect.
        !
        ! Perform preconditioning of the spatial defect with the method 
        ! provided by the core equation module.
        call stat_startTimer (rsolverNode%rtimeSpacePrecond)
        call cc_precondSpaceDefect (&
            rsolverNode%p_rsubnodeBlockFBGS%p_rspatialPreconditioner,&
            rnonlinearSpatialMatrix,rtempVectorRHS,&
            rtempVectorSol(1),rtempVectorSol(2),rtempVectorSol(3),&
            bsuccess,rsolverNode%p_rproblem%rcollection)      
        call stat_stopTimer (rsolverNode%rtimeSpacePrecond)
      
        ! Add that defect to the current solution -- damped by domega.
        call lsysbl_vectorLinearComb (rtempVectorRHS,rtempVectorX2,&
            rsolverNode%domega,1.0_DP)
      
        ! Save the new solution.
        call sptivec_setTimestepData (p_rx, 1+isubstep, rtempVectorX2)
        
        ! Shift the RHS/solution vectors: 1 -> 2 -> 3
        call lsysbl_copyVector (rtempVectorD2,rtempVectorD3)
        call lsysbl_copyVector (rtempVectorD1,rtempVectorD2)

        call lsysbl_copyVector (rtempVectorX2,rtempVectorX3)
        call lsysbl_copyVector (rtempVectorX1,rtempVectorX2)

        if (associated(p_rspaceTimeMatrix%p_rsolution)) then
          call lsysbl_copyVector (rtempVectorSol(2),rtempVectorSol(3))
          call lsysbl_copyVector (rtempVectorSol(1),rtempVectorSol(2))
        end if

      end do
      
      ! -----
      ! Forward in time
      ! -----
      
      ! rtempVectorSol(1) is undefined, but we don't need it.
      ! rtempVectorSol(2) holds the data of the 0th timestep, 
      ! rtempVectorSol(3) of the 1st one.
      
      ! Norm of the residuum
      dres = 0.0_DP

      ! Load the RHS and solution of the 0th timestep.
      call sptivec_getTimestepData (rd, 1+0, rtempVectorD2)
      
      ! Current iterate
      call sptivec_getTimestepData (p_rx, 1+0, rtempVectorX2)

      ! Loop through the substeps we have to update
      do isubstep = 0,p_rspaceTimeDiscr%NEQtime-1
      
        ! Current point in time
        dtime = p_rspaceTimeDiscr%rtimeDiscr%dtimeInit + isubstep * dtstep

        if (rsolverNode%ioutputLevel .ge. 1) then
          call output_line ('Space-Time-Block-FBGS preconditioning of timestep: '//&
              trim(sys_siL(isubstep,10))//&
              ', Time: '//trim(sys_sdL(dtime,10)))
        end if
      
        ! Update the boundary conditions in the preconditioner
        call cc_updatePreconditionerBC (rsolverNode%p_rproblem,&
            rsolverNode%p_rsubnodeBlockFBGS%p_rspatialPreconditioner,dtime)

        ! The RHS which is put into the preconditioner is set up in 
        ! rtempVectorRHS to prevent rtempVectorD2 from getting destroyed.
        call lsysbl_copyVector (rtempVectorD2,rtempVectorRHS)
        
        ! Is this the first timestep or not?
        if (isubstep .gt. 0) then
        
          ! Create d2 = RHS - Mx1 
          call cc_initNonlinMatrix (rnonlinearSpatialMatrix,rsolverNode%p_rproblem,&
              p_rspaceTimeMatrix%p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
              p_rspaceTimeMatrix%p_rspaceTimeDiscr%p_rlevelInfo%rstaticInfo)
          call cc_setupMatrixWeights (rsolverNode%p_rproblem,p_rspaceTimeMatrix,dtheta,&
            isubstep,-1,rnonlinearSpatialMatrix)
          
          call cc_assembleDefect (rnonlinearSpatialMatrix,rtempVectorX1,&
            rtempVectorRHS,domegaSOR,rtempVectorSol(1),rtempVectorSol(2),rtempVectorSol(3))
          
        end if

        ! Is this the last timestep or not?
        if (isubstep .lt. p_rspaceTimeDiscr%NEQtime-1) then
          
          if (isubstep .gt. 0) then
          
            ! Get the evaluation point for the nonlinearity in the next timestep
            ! into rtempVectorSol(3).
            ! In the 0th timestep, we don't need that; there, rtempVectorSol(3)
            ! is set from the previous backward loop.
            if (associated(p_rspaceTimeMatrix%p_rsolution)) then
              call sptivec_getTimestepData (p_rspaceTimeMatrix%p_rsolution, &
                  1+isubstep+1, rtempVectorSol(3))
            end if
          end if
          
          ! Read the RHS and solution of the next timestep
          call sptivec_getTimestepData (rd, 1+isubstep+1, rtempVectorD3)
          call sptivec_getTimestepData (p_rx, 1+isubstep+1, rtempVectorX3)
        
          ! Create d2 = RHS - Ml3 
          call cc_initNonlinMatrix (rnonlinearSpatialMatrix,rsolverNode%p_rproblem,&
              p_rspaceTimeMatrix%p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
              p_rspaceTimeMatrix%p_rspaceTimeDiscr%p_rlevelInfo%rstaticInfo)
          call cc_setupMatrixWeights (rsolverNode%p_rproblem,p_rspaceTimeMatrix,dtheta,&
            isubstep,1,rnonlinearSpatialMatrix)
          
          call cc_assembleDefect (rnonlinearSpatialMatrix,rtempVectorX3,&
            rtempVectorRHS,domegaSOR,&
            rtempVectorSol(1),rtempVectorSol(2),rtempVectorSol(3))
          
        end if

        ! Set up the matrix weights for the diagonal matrix
        call cc_initNonlinMatrix (rnonlinearSpatialMatrix,rsolverNode%p_rproblem,&
            p_rspaceTimeMatrix%p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
            p_rspaceTimeMatrix%p_rspaceTimeDiscr%p_rlevelInfo%rstaticInfo)
        call cc_setupMatrixWeights (rsolverNode%p_rproblem,p_rspaceTimeMatrix,dtheta,&
          isubstep,0,rnonlinearSpatialMatrix)
          
        ! Create d2 = RHS - A(solution) X2
        call cc_assembleDefect (rnonlinearSpatialMatrix,rtempVectorX2,rtempVectorRHS,&
            1.0_DP,rtempVectorSol(1),rtempVectorSol(2),rtempVectorSol(3))
            
        ! Filter the defect for BC's and initial conditions if necessary
        if (isubstep .eq. 0) then
          call tbc_implementInitCondDefSingle (p_rspaceTimeDiscr, rtempVectorRHS)
        else if (isubstep .eq. p_rspaceTimeDiscr%NEQtime-1) then
          call tbc_implementTermCondDefSingle (p_rspaceTimeDiscr, rtempVectorRHS)
        end if

        call tbc_implementSpatialBCdefect (&
            rsolverNode%p_rproblem,dtime,p_rspaceTimeDiscr,rtempVectorRHS)
        
        ! Ok, we have the local defect.
        ! Sum up the norm to the norm of the global vector.
        if (bcalcNorm) &
          dres = dres + lsysbl_vectorNorm (rtempVectorRHS,LINALG_NORML2)**2
        
        ! Perform preconditioning of the spatial defect with the method 
        ! provided by the core equation module.
        call stat_startTimer (rsolverNode%rtimeSpacePrecond)
        call cc_precondSpaceDefect (&
            rsolverNode%p_rsubnodeBlockFBGS%p_rspatialPreconditioner,&
            rnonlinearSpatialMatrix,rtempVectorRHS,&
            rtempVectorSol(1),rtempVectorSol(2),rtempVectorSol(3),&
            bsuccess,rsolverNode%p_rproblem%rcollection)      
        call stat_stopTimer (rsolverNode%rtimeSpacePrecond)
      
        ! Add that defect to the current solution -- damped by domega.
        call lsysbl_vectorLinearComb (rtempVectorRHS,rtempVectorX2,&
            rsolverNode%domega,1.0_DP)
      
        ! Save the new solution.
        call sptivec_setTimestepData (p_rx, 1+isubstep, rtempVectorX2)
        
        ! Shift the RHS/solution vectors: 1 <- 2 <- 3
        call lsysbl_copyVector (rtempVectorD2,rtempVectorD1)
        call lsysbl_copyVector (rtempVectorD3,rtempVectorD2)

        call lsysbl_copyVector (rtempVectorX2,rtempVectorX1)
        call lsysbl_copyVector (rtempVectorX3,rtempVectorX2)

        if (associated(p_rspaceTimeMatrix%p_rsolution)) then
          call lsysbl_copyVector (rtempVectorSol(2),rtempVectorSol(1))
          call lsysbl_copyVector (rtempVectorSol(3),rtempVectorSol(2))
        end if

      end do
      
    end do ! iiteration
    
    ! Overwrite the rd by our solution.
    call sptivec_copyVector (p_rx,rd)
    
    ! Release memory, finish.
    call lsysbl_releaseVector (rtempVectorRHS)
    call lsysbl_releaseVector (rtempVectorSol(3))
    call lsysbl_releaseVector (rtempVectorSol(2))
    call lsysbl_releaseVector (rtempVectorSol(1))
    call lsysbl_releaseVector (rtempVectorD3)
    call lsysbl_releaseVector (rtempVectorD2)
    call lsysbl_releaseVector (rtempVectorD1)
    call lsysbl_releaseVector (rtempVectorX3)
    call lsysbl_releaseVector (rtempVectorX2)
    call lsysbl_releaseVector (rtempVectorX1)
    
  end subroutine

! *****************************************************************************
! Routines for the CG solver
! *****************************************************************************

!<subroutine>
  
  subroutine sptils_initCG (rproblem,p_rsolverNode,p_rpreconditioner)
  
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
  type(t_problem), target :: rproblem

  ! OPTIONAL: A pointer to the solver structure of a solver that should be 
  ! used for preconditioning. If not given or set to NULL(), no preconditioning 
  ! will be used.
  type(t_sptilsNode), pointer, optional   :: p_rpreconditioner
!</input>
  
!<output>
  ! A pointer to a t_linsolNode structure. Is set by the routine, any previous
  ! value of the pointer is destroyed.
  type(t_sptilsNode), pointer         :: p_rsolverNode
!</output>
  
!</subroutine>

    ! Create a default solver structure
    call sptils_initSolverGeneral(rproblem,p_rsolverNode)
    
    ! Initialise the type of the solver
    p_rsolverNode%calgorithm = SPTILS_ALG_CG
    
    ! Initialise the ability bitfield with the ability of this solver:
    p_rsolverNode%ccapability = LINSOL_ABIL_CHECKDEF + &
                                LINSOL_ABIL_USESUBSOLVER
    
    ! Allocate the subnode for CG.
    ! This initialises most of the variables with default values appropriate
    ! to this solver.
    allocate(p_rsolverNode%p_rsubnodeCG)
    
    ! Attach the preconditioner if given. 
    
    if (present(p_rpreconditioner)) then 
      p_rsolverNode%p_rsubnodeCG%p_rpreconditioner => p_rpreconditioner
    end if

  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>
  
  recursive subroutine sptils_setMatrixCG (rsolverNode,Rmatrices)
  
!<description>
  
  ! This routine is called if the system matrix changes.
  ! The routine calls linsol_setMatrices for the preconditioner of CG
  ! to inform also that one about the change of the matrix pointer.
  
!</description>
  
!<input>
  ! An array of system matrices which is simply passed to the initialisation 
  ! routine of the preconditioner.
  type(t_ccoptSpaceTimeMatrix), dimension(:), intent(IN)   :: Rmatrices
!</input>
  
!<inputoutput>
  ! The t_linsolNode structure of the CG solver
  type(t_sptilsNode), intent(INOUT)         :: rsolverNode
!</inputoutput>
  
!</subroutine>

    ! Do we have a preconditioner given? If yes, call the general initialisation
    ! routine to initialise it.
    
    if (associated(rsolverNode%p_rsubnodeCG%p_rpreconditioner)) then
      call sptils_setMatrices (rsolverNode%p_rsubnodeCG%p_rpreconditioner, &
                               Rmatrices)
    end if

  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>
  
  recursive subroutine sptils_initStructureCG (rsolverNode, ierror)
  
!<description>
  ! Calls the initStructure subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_initStructure.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the CG solver
  type(t_sptilsNode), intent(INOUT)         :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the SPTILS_ERR_XXXX constants. A value different to 
  ! SPTILS_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  integer, intent(OUT) :: ierror
!</output>

!</subroutine>

    ! local variables
    integer :: i,NEQtime
    type(t_sptilsSubnodeCG), pointer :: p_rsubnode
    
    ! A-priori we have no error...
    ierror = SPTILS_ERR_NOERROR

    ! We simply pass isubgroup to the subsolvers when calling them.
    ! Inside of this routine, there's not much to do with isubgroup,
    ! as there is no special information (like a factorisation)
    ! associated with this type of solver, which has to be allocated,
    ! released or updated...

    ! Call the init routine of the preconditioner.
    if (associated(rsolverNode%p_rsubnodeCG%p_rpreconditioner)) then
      call sptils_initStructure (rsolverNode%p_rsubnodeCG%p_rpreconditioner,ierror)
    end if
    
    ! CG needs 3 temporary vectors + 1 for preconditioning. 
    ! Allocate that here! Use the default data type prescribed in the solver 
    ! structure for allocating the temp vectors.
    p_rsubnode => rsolverNode%p_rsubnodeCG
    NEQtime = rsolverNode%rmatrix%p_rspaceTimeDiscr%NEQtime
    do i=1,4
      call sptivec_initVector (p_rsubnode%RtempVectors(i),NEQtime,&
        rsolverNode%rmatrix%p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation)
    end do
  
    ! Allocate memory for a spatial temp vector.
    call lsysbl_createVecBlockByDiscr (&
        rsolverNode%rmatrix%p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
        rsolverNode%p_rsubnodeCG%rtempVectorSpace)

  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>
  
  recursive subroutine sptils_initDataCG (rsolverNode, ierror)
  
!<description>
  ! Calls the initData subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_initData.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the CG solver
  type(t_sptilsNode), intent(INOUT)         :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the SPTILS_ERR_XXXX constants. A value different to 
  ! SPTILS_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  integer, intent(OUT) :: ierror
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
    if (associated(rsolverNode%p_rsubnodeCG%p_rpreconditioner)) then
      call sptils_initData (rsolverNode%p_rsubnodeCG%p_rpreconditioner, &
                            ierror)
    end if
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine sptils_doneDataCG (rsolverNode)
  
!<description>
  ! Calls the doneData subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_doneData.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the CG solver
  type(t_sptilsNode), intent(INOUT)         :: rsolverNode
!</inputoutput>
  
!</subroutine>

    ! We simply pass isubgroup to the subsolvers when calling them.
    ! Inside of this routine, there's not much to do with isubgroup,
    ! as there is no special information (like a factorisation)
    ! associated with this type of solver, which has to be allocated,
    ! released or updated...

    ! Call the done routine of the preconditioner.
    if (associated(rsolverNode%p_rsubnodeCG%p_rpreconditioner)) then
      call sptils_doneData (rsolverNode%p_rsubnodeCG%p_rpreconditioner)
    end if
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine sptils_doneStructureCG (rsolverNode)
  
!<description>
  ! Calls the doneStructure subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_doneStructure.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the CG solver
  type(t_sptilsNode), intent(INOUT)         :: rsolverNode
!</inputoutput>
  
!</subroutine>

    ! local variables
    integer :: i
    type(t_sptilsSubnodeCG), pointer :: p_rsubnode
    
    ! We simply pass isubgroup to the subsolvers when calling them.
    ! Inside of this routine, there's not much to do with isubgroup,
    ! as there is no special information (like a factorisation)
    ! associated with this type of solver, which has to be allocated,
    ! released or updated...

    ! Call the init routine of the preconditioner.
    if (associated(rsolverNode%p_rsubnodeCG%p_rpreconditioner)) then
      call sptils_doneStructure (rsolverNode%p_rsubnodeCG%p_rpreconditioner)
    end if
    
    ! Release temporary data if associated
    p_rsubnode => rsolverNode%p_rsubnodeCG
    if (p_rsubnode%RtempVectors(1)%NEQtime .ne. 0) then
      do i=4,1,-1
        call sptivec_releaseVector (p_rsubnode%RtempVectors(i))
      end do
    end if
    
    ! Release the spatial temp vector
    if (rsolverNode%p_rsubnodeCG%rtempVectorSpace%NEQ .ne. 0) then
      call lsysbl_releaseVector (rsolverNode%p_rsubnodeCG%rtempVectorSpace)
    end if
    
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine sptils_doneCG (rsolverNode)
  
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
  type(t_sptilsNode), pointer         :: rsolverNode
!</input>
  
!</subroutine>

    ! Check if there's a preconditioner attached. If yes, release it.
    if (associated(rsolverNode%p_rsubnodeCG%p_rpreconditioner)) then
      call sptils_releaseSolver(rsolverNode%p_rsubnodeCG%p_rpreconditioner)
    end if
    
    ! Release memory if still associated
    call sptils_doneDataCG (rsolverNode)
    call sptils_doneStructureCG (rsolverNode)
    
    ! Release the CG subnode
    deallocate(rsolverNode%p_rsubnodeCG)

  end subroutine

  ! ***************************************************************************
  
!<subroutine>
  
  recursive subroutine sptils_precCG (rsolverNode,rd)
  
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
  type(t_sptilsNode), intent(INOUT), target :: rsolverNode
   
  ! On call to this routine: The defect vector to be preconditioned.
  ! Will be overwritten by the preconditioned defect.
  type(t_spacetimeVector), intent(INOUT)        :: rd
!</inputoutput>
  
!</subroutine>

  ! local variables
  real(DP) :: dalpha,dbeta,dgamma,dgammaOld,dres,dfr
  integer :: ite

  ! The system matrix
  type(t_ccoptSpaceTimeDiscretisation), pointer :: p_rspaceTimeDiscr
  type(t_ccoptSpaceTimeMatrix), pointer :: p_rmatrix
  
  ! Minimum number of iterations, print-sequence for residuals
  integer :: nminIterations, niteResOutput
  
  ! Whether to precondition
  logical bprec
  
  ! Our structure
  type(t_sptilsSubnodeCG), pointer :: p_rsubnode
  
  ! Pointers to temporary vectors - named for easier access
  type(t_spacetimeVector), pointer :: p_DR,p_DP,p_DD,p_rx
  type(t_sptilsNode), pointer :: p_rprecSubnode
  
  ! DEBUG!!!
  real(DP), dimension(:), pointer :: p_DRdata,p_DPdata,p_DDdata,p_rxdata,p_rddata
  
    ! Solve the system!
  
    ! Status reset
    rsolverNode%iresult = 0
    call stat_clearTimer (rsolverNode%rtimeSpacePrecond)
    
    ! Get some information
    p_rsubnode => rsolverNode%p_rsubnodeCG
    p_rspaceTimeDiscr => rsolverNode%rmatrix%p_rspaceTimeDiscr
    p_rmatrix => rsolverNode%rmatrix

    ! Check the parameters
    if ((rd%NEQtime .eq. 0) .or. &
        (p_rspaceTimeDiscr%rtimeDiscr%nintervals .eq. 0) ) then
    
      ! Parameters wrong
      rsolverNode%iresult = 2
      return
    end if

    ! Minimum number of iterations
 
    nminIterations = max(rsolverNode%nminIterations,0)
      
    ! Use preconditioning? Filtering?

    bprec = associated(rsolverNode%p_rsubnodeCG%p_rpreconditioner)
    
    ! Iteration when the residuum is printed:

    niteResOutput = max(1,rsolverNode%niteResOutput)

    ! Set pointers to the temporary vectors
    ! residuum vector
    p_DR   => p_rsubnode%RtempVectors(1)
    ! preconditioned residuum vector
    p_DP   => p_rsubnode%RtempVectors(2)
    ! direction vector
    p_DD   => p_rsubnode%RtempVectors(3)
    ! solution vector
    p_rx   => p_rsubnode%RtempVectors(4)
    
    if (bprec) then
      p_rprecSubnode => p_rsubnode%p_rpreconditioner
    end if
    
    ! rd is our RHS. p_rx points to a new vector which will be our
    ! iteration vector. At the end of this routine, we replace
    ! rd by p_rx.
    ! Clear our iteration vector p_rx.
    call sptivec_clearVector (p_rx)
      
    ! Initialize used vectors with zero
      
    call sptivec_clearVector(p_DR)
    call sptivec_clearVector(p_DP)
    call sptivec_clearVector(p_DD)
    
    ! Initialization

    dalpha = 1.0_DP
    dbeta  = 1.0_DP
    dgamma = 1.0_DP
    dgammaOld = 1.0_DP

    ! Copy our RHS rd to p_DR. As the iteration vector is 0, this
    ! is also our initial defect.

    call sptivec_copyVector(rd,p_DR)
    
    ! Filter the defect for boundary conditions in space and time.
    call tbc_implementInitCondDefect (&
        p_rspaceTimeDiscr,p_DR,p_rsubnode%rtempVectorSpace)
    call tbc_implementBCdefect (rsolverNode%p_rproblem,&
        p_rspaceTimeDiscr,p_DR,p_rsubnode%rtempVectorSpace)
    
    ! Get the norm of the residuum
    dres = sptivec_vectorNorm (p_DR,rsolverNode%iresNorm)
    if (.not.((dres .ge. 1E-99_DP) .and. &
              (dres .le. 1E99_DP))) dres = 0.0_DP

    ! Initialize starting residuum
      
    rsolverNode%dinitialDefect = dres
    rsolverNode%dlastDefect = 0.0_DP

    ! Check if out initial defect is zero. This may happen if the filtering
    ! routine filters "everything out"!
    ! In that case we can directly stop our computation.

    if ( rsolverNode%dinitialDefect .lt. rsolverNode%drhsZero ) then
     
      ! final defect is 0, as initialised in the output variable above

      call sptivec_clearVector(p_rx)
      ite = 0
      rsolverNode%dfinalDefect = dres
          
    else

      if (rsolverNode%ioutputLevel .ge. 2) then
        call output_line ('Space-Time-CG: Iteration '// &
             trim(sys_siL(0,10))//',  !!RES!! = '//&
             trim(sys_sdEL(rsolverNode%dinitialDefect,15)) )
      end if

      ! Copy the residuum vector p_DR to the preconditioned one.
      call sptivec_copyVector(p_DR,p_DP)

      if (bprec) then
        ! Perform preconditioning with the assigned preconditioning
        ! solver structure.
        call sptils_precondDefect (p_rprecSubnode,p_DP)
        call stat_addTimers (p_rprecSubnode%rtimeSpacePrecond,rsolverNode%rtimeSpacePrecond)
      end if

      ! Calculate gamma
      dgamma = sptivec_scalarProduct (p_DR, p_DP)
      
      ! Perform at most nmaxIterations loops to get a new vector
      
      ! Copy the preconditioned residual vector to the direction vector
      call sptivec_copyVector (p_DP, p_DD)

      do ite = 1,rsolverNode%nmaxIterations
      
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
        call cc_spaceTimeMatVec (rsolverNode%p_rproblem, p_rmatrix, &
            p_DD,p_DP, 1.0_DP,0.0_DP,SPTID_FILTER_NONE)
        
        ! Calculate alpha for the CG iteration
        dalpha = dgamma / sptivec_scalarProduct (p_DD, p_DP)

        if (dalpha .eq. 0.0_DP) then
          ! We are below machine exactness - we can't do anything more...
          ! May happen with very small problems with very few unknowns!
          if (rsolverNode%ioutputLevel .ge. 2) then
            call output_line('Space-Time-CG: Convergence failed, ALPHA=0!')
            rsolverNode%iresult = -2
            exit
          end if
        end if
        
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! STEP 2: Calculate x_{k+1}
        !
        ! x_{k+1} := x_k + alpha * d_k
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        call sptivec_vectorLinearComb (p_DD, p_rx, dalpha, 1.0_DP)

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! STEP 3: Calculate g_{k+1}
        !
        ! g_{k+1} := g_k - alpha * A * d_k
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Since we have abused p_DP for saving A * d_k, we can use it now
        ! for the linear combination of g_{k+1}
        call sptivec_vectorLinearComb (p_DP, p_DR, -dalpha, 1.0_DP)

        ! Filter the defect for boundary conditions in space and time.
        call tbc_implementInitCondDefect (&
            p_rspaceTimeDiscr,p_DR,p_rsubnode%rtempVectorSpace)
        call tbc_implementBCdefect (rsolverNode%p_rproblem,&
           p_rspaceTimeDiscr,p_DR,p_rsubnode%rtempVectorSpace)

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! STEP 4: Calculate ||g_{k+1}|| and write some output
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Get the norm of the new (final?) residuum
        dfr = sptivec_vectorNorm (p_DR,rsolverNode%iresNorm)
     
        rsolverNode%dlastDefect = rsolverNode%dfinalDefect
        rsolverNode%dfinalDefect = dfr

        ! Test if the iteration is diverged
        if (sptils_testDivergence(rsolverNode,dfr)) then
          call output_line('Space-Time-CG: Solution diverging!')
          rsolverNode%iresult = 1
          exit
        end if
     
        ! At least perform nminIterations iterations
        if (ite .ge. nminIterations) then
        
          ! Check if the iteration converged
          if (sptils_testConvergence(rsolverNode,dfr)) exit
          
        end if

        ! print out the current residuum

        if ((rsolverNode%ioutputLevel .ge. 2) .and. &
            (mod(ite,niteResOutput).eq.0)) then
          call output_line ('Space-Time-CG: Iteration '// &
              trim(sys_siL(ITE,10))//',  !!RES!! = '//&
              trim(sys_sdEL(rsolverNode%dfinalDefect,15)) )
        end if
        
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! STEP 5: Calculate h_{k+1}
        !
        ! In other words: Copy the residual vector and apply the
        ! preconditioner, if given.
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        ! Copy the residuum vector p_DR to the preconditioned one.
        call sptivec_copyVector(p_DR,p_DP)
    
        if (bprec) then
          ! Perform preconditioning with the assigned preconditioning
          ! solver structure.
          call sptils_precondDefect (p_rprecSubnode,p_DP)
          call stat_addTimers (p_rprecSubnode%rtimeSpacePrecond,rsolverNode%rtimeSpacePrecond)
        end if
        
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
        call sptivec_vectorLinearComb (p_DP, p_DD, 1.0_DP, dbeta)
        
        ! That's it - next iteration!
      end do

      ! Set ITE to NIT to prevent printing of "NIT+1" of the loop was
      ! completed

      if (ite .gt. rsolverNode%nmaxIterations) &
        ite = rsolverNode%nmaxIterations

      ! Finish - either with an error or if converged.
      ! Print the last residuum.


      if ((rsolverNode%ioutputLevel .ge. 2) .and. &
          (ite .ge. 1) .and. (ITE .lt. rsolverNode%nmaxIterations) .and. &
          (rsolverNode%iresult .ge. 0)) then
        call output_line ('Space-Time-CG: Iteration '// &
            trim(sys_siL(ITE,10))//',  !!RES!! = '//&
            trim(sys_sdEL(rsolverNode%dfinalDefect,15)) )
      end if

    end if

    rsolverNode%iiterations = ite
    
    ! Overwrite our previous RHS by the new correction vector p_rx.
    ! This completes the preconditioning.
    call sptivec_copyVector (p_rx,rd)
      
    ! Don't calculate anything if the final residuum is out of bounds -
    ! would result in NaN's,...
      
    if (rsolverNode%dfinalDefect .lt. 1E99_DP) then
    
      ! If the initial defect was zero, the solver immediately
      ! exits - and so the final residuum is zero and we performed
      ! no steps; so the resulting convergence rate stays zero.
      ! In the other case the convergence rate computes as
      ! (final defect/initial defect) ** 1/nit :

      if (rsolverNode%dfinalDefect .gt. rsolverNode%drhsZero) then
        rsolverNode%dconvergenceRate = &
                    (rsolverNode%dfinalDefect / rsolverNode%dinitialDefect) ** &
                    (1.0_DP/real(rsolverNode%iiterations,DP))
      end if

      if (rsolverNode%ioutputLevel .ge. 2) then
        call output_lbrk()
        call output_line ('Space-Time-CG statistics:')
        call output_lbrk()
        call output_line ('Iterations              : '//&
             trim(sys_siL(rsolverNode%iiterations,10)) )
        call output_line ('!!INITIAL RES!!         : '//&
             trim(sys_sdEL(rsolverNode%dinitialDefect,15)) )
        call output_line ('!!RES!!                 : '//&
             trim(sys_sdEL(rsolverNode%dfinalDefect,15)) )
        if (rsolverNode%dinitialDefect .gt. rsolverNode%drhsZero) then     
          call output_line ('!!RES!!/!!INITIAL RES!! : '//&
            trim(sys_sdEL(rsolverNode%dfinalDefect / rsolverNode%dinitialDefect,15)) )
        else
          call output_line ('!!RES!!/!!INITIAL RES!! : '//&
               trim(sys_sdEL(0.0_DP,15)) )
        end if
        call output_lbrk ()
        call output_line ('Rate of convergence     : '//&
             trim(sys_sdEL(rsolverNode%dconvergenceRate,15)) )

      end if

      if (rsolverNode%ioutputLevel .eq. 1) then
        call output_line (&
              'Space-Time-CG: Iterations/Rate of convergence: '//&
              trim(sys_siL(rsolverNode%iiterations,10))//' /'//&
              trim(sys_sdEL(rsolverNode%dconvergenceRate,15)) )
      end if

    else
      ! DEF=Infinity; RHO=Infinity, set to 1
      rsolverNode%dconvergenceRate = 1.0_DP
    end if  
  
  end subroutine

! *****************************************************************************
! Routines for the UMFPACK4 preconditioner
! *****************************************************************************

!<subroutine>
  
  subroutine sptils_initUMFPACK4 (rproblem,rsolverNode)
  
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
  type(t_problem), target :: rproblem
!</input>
  
!<output>
  
  ! A pointer to a t_sptilsNode structure. Is set by the routine, any previous
  ! value of the pointer is destroyed.
  type(t_sptilsNode), pointer         :: rsolverNode
   
!</output>
  
!</subroutine>

  ! Create a default solver structure
  
  call sptils_initSolverGeneral(rproblem,rsolverNode)
  
  ! Initialise the type of the solver
  rsolverNode%calgorithm = SPTILS_ALG_UMFPACK4
  
  ! Initialise the ability bitfield with the ability of this solver:
  rsolverNode%ccapability = SPTILS_ABIL_DIRECT
  
  ! Allocate the subnode for UMFPACK4.
  ! This initialises most of the variables with default values appropriate
  ! to this solver.
  allocate(rsolverNode%p_rsubnodeUMFPACK4)
  
  ! Initialize the UMFPACK4 control structure:
  call UMF4DEF(rsolverNode%p_rsubnodeUMFPACK4%Dcontrol)
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine sptils_initStructureUMFPACK4 (rsolverNode,ierror)
  
!<description>
  ! Performs a symbolic factorisation on the assigned matrix.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure of the UMFPACK4 solver
  type(t_sptilsNode), intent(INOUT),target  :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the SPTILS_ERR_XXXX constants. A value different to 
  ! SPTILS_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  integer, intent(OUT) :: ierror
!</output>
  
!</subroutine>

  ! A-priori we have no error...
  ierror = SPTILS_ERR_NOERROR

  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine sptils_initDataUMFPACK4 (rsolverNode, ierror)
  
!<description>
  ! Performs a numeric factorisation on the assigned matrix.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure of the UMFPACK4 solver
  type(t_sptilsNode), intent(INOUT), target :: rsolverNode
!</inputoutput>

!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to 
  ! SPTILS_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  integer, intent(OUT) :: ierror
!</output>
  
!</subroutine>

    ! local variables
    type(t_matrixScalar), pointer :: p_rmatrix
    type(t_matrixBlock), target :: rmatrixGlobal
    type(t_matrixBlock) :: rtempMatrix
    integer, dimension(:), pointer :: p_Kld
    integer, dimension(:), pointer :: p_Kcol
    real(DP), dimension(:), pointer :: p_DA
    integer :: ixfirst,ixlast,iyfirst,iylast
    type(t_matrixBlock) :: rsubmatrix,rtempSubmatrix
    integer, dimension(:), allocatable :: IpureDirichletTimesteps
    integer :: npureDirichletTimesteps

    ! Status variables of UMFPACK4; receives the UMFPACK-specific return code
    ! of a call to the solver routines.
    real(DP), dimension(90) :: Dinfo
    
    ! A-priori we have no error...
    ierror = SPTILS_ERR_NOERROR
    
    allocate(IpureDirichletTimesteps(rsolverNode%rmatrix%p_rspaceTimeDiscr%NEQtime))

    ! We have to create a global matrix first!
    ! Call the assembly routine to assemble the global block matrix.
    ! Afterwards, reshape the data to form a scalar matrix which
    ! can be feed to UMFPACK.
    call assembleGlobalSpaceTimeMatrix (rsolverNode%p_rproblem,&
        rsolverNode%rmatrix,rmatrixGlobal,&
        IpureDirichletTimesteps,npureDirichletTimesteps)

    ! Extract a submatrix for file output?
    if (rsolverNode%p_rsubnodeUMFPACK4%cwriteMatrix .gt. 0) then
      
      ixfirst = rsolverNode%p_rsubnodeUMFPACK4%ixfirst
      ixlast  = rsolverNode%p_rsubnodeUMFPACK4%ixlast
      iyfirst = rsolverNode%p_rsubnodeUMFPACK4%iyfirst
      iylast  = rsolverNode%p_rsubnodeUMFPACK4%iylast
      if ((ixfirst .ne. 0) .or. (ixlast .ne. 0) .or. &
          (iyfirst .ne. 0) .or. (iylast .ne. 0)) then
        if (ixfirst .eq. 0) ixfirst = 1
        if (iyfirst .eq. 0) iyfirst = 1
        if (ixlast .eq. 0) ixlast = rmatrixGlobal%nblocksPerRow
        if (iylast .eq. 0) iylast = rmatrixGlobal%nblocksPerCol
        call lsysbl_deriveSubmatrix (rmatrixGlobal,rsubmatrix,&
            LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE,&
            iyfirst,iylast,ixfirst,ixlast)
            
        ! Assemble a separate submatrix
        call glsys_assembleGlobal (rsubmatrix, rtempsubMatrix, .true., .true.)

        ! The global block matrix is not needed anymore.
        call lsysbl_releaseMatrix (rsubmatrix)
      end if
    end if

    call glsys_assembleGlobal (rmatrixGlobal, rtempMatrix, .true., .true.)
    
    ! The global block matrix is not needed anymore.
    call lsysbl_releaseMatrix (rmatrixGlobal)
    
    ! Implement Dirichlet boundary conditions into the pressure matrix.
    ! This is necessary if we have pure Dirichlet BC's in the velocity,
    ! so the pressure is not uniquely defined. In this case, we fix
    ! the pressure in the first node by writing a unit line into
    ! the global matrix for the first pressure DOF in every timestep.
    call pressureDirichlet (rsolverNode%p_rproblem,&
        rsolverNode%rmatrix,rtempMatrix,&
        IpureDirichletTimesteps,npureDirichletTimesteps)

    !CALL matio_writeBlockMatrixHR (rtempMatrix, 'matrix',&
    !                               .TRUE., 0, 'matrix.txt', '(E15.5)')
    !CALL matio_writeBlockMatrixMaple (rtempMatrix, 'A',&
    !                             0, 'matrix.txt', '(E20.10)')

    if (rsolverNode%p_rsubnodeUMFPACK4%cwriteMatrix .gt. 0) then
      
      ! Extract a submatrix?
      if ((ixfirst .ne. 0) .or. (ixlast .ne. 0) .or. &
          (iyfirst .ne. 0) .or. (iylast .ne. 0)) then
        p_rmatrix => rtempsubMatrix%RmatrixBlock(1,1)
      else
        p_rmatrix => rtempMatrix%RmatrixBlock(1,1)
      end if
    
      call lsyssc_getbase_double (p_rmatrix,p_DA)
      where (abs(p_Da) .lt. 1.0E-12_DP) p_Da = 0.0_DP
      call matio_writeMatrixHR (p_rmatrix, 'matrix',&
                                .true., 0, rsolverNode%p_rsubnodeUMFPACK4%sfilename,&
                                '(E15.5)')
      if (rsolverNode%p_rsubnodeUMFPACK4%cwriteMatrix .eq. 2) then
        call sys_halt()
      end if
  
      ! Cleanup    
      if ((ixfirst .ne. 0) .or. (ixlast .ne. 0) .or. &
          (iyfirst .ne. 0) .or. (iylast .ne. 0)) then
        ! Change the pointer back to the original matrix
        call lsysbl_releaseMatrix (rtempsubMatrix)
      end if
      
    end if

    ! Now start to modify the temp matrix for UMFPACK's needs.  
    p_rmatrix => rtempMatrix%RmatrixBlock(1,1)
    
    !CALL matio_spyMatrix('matrix','matrix',p_rmatrix,.TRUE.)

    ! Modify Kcol/Kld of the matrix. Subtract 1 to get them 0-based.
    call lsyssc_addIndex (p_rmatrix%h_Kcol,-1_I32)
    call lsyssc_addIndex (p_rmatrix%h_Kld,-1_I32)
    
    ! Get the data arrays.
    call lsyssc_getbase_Kcol (p_rmatrix,p_Kcol)
    call lsyssc_getbase_Kld (p_rmatrix,p_Kld)
    call lsyssc_getbase_double (p_rmatrix,p_DA)
    
    ! Perform a symbolic factorization...
    call UMF4SYM(rtempMatrix%NEQ,rtempMatrix%NCOLS,p_Kld,p_Kcol,p_Da, &
                rsolverNode%p_rsubnodeUMFPACK4%isymbolic,&
                rsolverNode%p_rsubnodeUMFPACK4%Dcontrol,&
                Dinfo)
                 
    ! Check Dinfo(1) if there is an error
    select case (int(Dinfo(1)))
    case (0)
    
      ! Perform a numeric factorization...
      call UMF4NUM(p_Kld,p_Kcol,p_Da, &
              rsolverNode%p_rsubnodeUMFPACK4%isymbolic,&
              rsolverNode%p_rsubnodeUMFPACK4%inumeric,&
              rsolverNode%p_rsubnodeUMFPACK4%Dcontrol,&
              Dinfo)
              
      ! Check Dinfo(1) if there is an error
      select case (int(Dinfo(1)))
      case (0)
        ! ok.
      case (1)
        ! Singular matrix
        ierror = LINSOL_ERR_SINGULAR
      case (-1)
        ! no memory
        ierror = LINSOL_ERR_NOMEMORY
      case (-11)
        ! no memory
        ierror = LINSOL_ERR_MATRIXHASCHANGED
      case DEFAULT
        ! don't know what went wrong
        ierror = LINSOL_ERR_INITERROR
      end select
    
    case (1)
      ! Singular matrix
      ierror = LINSOL_ERR_SINGULAR
      return
    case (-1)
      ! no memory
      ierror = LINSOL_ERR_NOMEMORY
    case DEFAULT
      ! don't know what went wrong
      ierror = LINSOL_ERR_INITERROR
    end select

    ! Throw away the temporary matrices
    call lsysbl_releaseMatrix (rtempMatrix)
    
    deallocate(IpureDirichletTimesteps)

  contains

    ! ---------------------------------------------------------------

    subroutine pressureDirichlet (rproblem,rsupermatrix,rmatrix,Itimesteps,ntimesteps)
    
    ! If we have a Pure-Dirichlet problem, this replaces the row of the first
    ! pressure DOF in every timestep by zero.
    
    ! The problem structure that defines the problem.
    type(t_problem), intent(INOUT) :: rproblem

    ! The source space-time matrix where rmatrix was generated from
    type(t_ccoptSpaceTimeMatrix), intent(IN) :: rsupermatrix

    ! The global space time matrix to be modified.
    ! Must be a 1x1 block matrix with the submatrix (1,1) representing
    ! the global space time matrix.
    type(t_matrixBlock), intent(INOUT) :: rmatrix
    
    ! List of timesteps which have no Neumann boundary in the BC`s.
    integer, dimension(:), intent(in) :: Itimesteps
    
    ! Number of timesteps in Itimesteps
    integer, intent(in) :: ntimesteps
    
      integer, dimension(:), allocatable :: Iidx
      integer :: ivelSize,ipSize,neq
      integer :: i,idof
      type(t_ccoptSpaceTimeDiscretisation), pointer :: p_rspaceTimeDiscr
    
      ! Nothing to do if there is Neumann boundary here.
      if (ntimesteps .eq. 0) return
    
      p_rspaceTimeDiscr => rsupermatrix%p_rspaceTimeDiscr
    
      neq = dof_igetNDofGlobBlock(p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation)
      ivelSize = 2 * p_rspaceTimeDiscr%p_rlevelInfo%rstaticInfo%rmatrixStokes%NEQ
      ipSize = p_rspaceTimeDiscr%p_rlevelInfo%rstaticInfo%rmatrixB1%NCOLS
    
      ! Create an array containing all the rows where a unit vector is to be
      ! imposed. Size = 2 * number of timesteps
      ! (primal + dual pressure, not the initial condition)
      allocate(Iidx(2*ntimesteps))
      !allocate(Iidx(0:2*(p_rspaceTimeDiscr%NEQtime-1)))
      
      do i=1,ntimesteps
        idof = Itimesteps(i)-1
        if (idof .eq. 0) then
          ! 0th time step -> only dual pressure because of init. cond.;
          ! implement the unit vector twice, so there's no special handling ti
          ! be done to the array indices.
          Iidx(2*i-1) = idof*neq+ivelSize+1 ! idof*neq+neq-ipSize+1
          Iidx(2*i  ) = idof*neq+neq-ipSize+1
        else
          Iidx(2*i-1) = idof*neq+ivelSize+1
          Iidx(2*i  ) = idof*neq+neq-ipSize+1
        end if
      end do
      
!      ! Add the row of the first primal/dual pressure DOF to that array.
!      Iidx(0) = neq-ipSize+1  ! 0th time step -> only dual pressure because of init. cond.
!      do i=1,p_rspaceTimeDiscr%NEQtime-1
!        Iidx(2*i-1) = i*neq+ivelSize+1
!        Iidx(2*i  ) = i*neq+neq-ipSize+1
!      end do
      
      ! Filter the matrix, implement the unit vectors.
      call mmod_replaceLinesByUnit (rmatrix%RmatrixBlock(1,1),Iidx)
      
      ! Release memory, finish.
      deallocate(Iidx)
    
    end subroutine

  end subroutine

  ! ***************************************************************************

  subroutine assembleGlobalSpaceTimeMatrix (rproblem,rsupermatrix,rmatrix,&
      IpureDirichletTimesteps,npureDirichletTimesteps)
  
  ! Assembles a block matrix rmatrix from a space-time matrix rsupermatrix
  ! by plugging together all blocks.
  ! Implement boundarys conditions into rmatrix during the assembly.
  ! WARNING: SPACE CONSUMING!!!
  
  ! The problem structure that defines the problem.
  type(t_problem), intent(INOUT) :: rproblem

  ! The source space-time matrix
  type(t_ccoptSpaceTimeMatrix), intent(IN), target :: rsupermatrix
  
  ! The destination block matrix
  type(t_matrixBlock), intent(OUT) :: rmatrix
  
  ! Returns a list of all timesteps which have pure Dirichlet BC`s.
  integer, dimension(:), intent(inout) :: IpureDirichletTimesteps
  
  ! Number of timesteps with pure Dirichlet BC`s.
  integer, intent(out) :: npureDirichletTimesteps

    ! local variables  
    type(t_matrixBlock) :: rblockTemp
    type(t_nonlinearSpatialMatrix) :: rnonlinearSpatialMatrix
    type(t_vectorBlock) :: rvector1,rvector2,rvector3
    integer :: isubstep,ileft,iright,ix
    integer, dimension(1) :: Irows
    integer :: idiag
    type(t_ccoptSpaceTimeDiscretisation), pointer :: p_rspaceTimeDiscr
    type(t_ccPreconditionerSpecials) :: rprecSpecials
    real(dp) :: dtime
    type(t_discreteBC) :: rdiscreteBC
    type(t_discreteFBC) :: rdiscreteFBC
    logical :: bhasNeumann
    integer :: ivelSize,ipSize,neq
    integer, dimension(2) :: Iidx

    ! DEBUG!!!
    real(dp), dimension(:), pointer :: p_Ddata1,p_Ddata2,p_Ddata3
      
    ! Initialise the collection for the assembly.
    call cc_initCollectForAssembly (rproblem,dtime,rproblem%rcollection)

    p_rspaceTimeDiscr => rsupermatrix%p_rspaceTimeDiscr
   
    neq = dof_igetNDofGlobBlock(p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation)
    ivelSize = 2 * p_rspaceTimeDiscr%p_rlevelInfo%rstaticInfo%rmatrixStokes%NEQ
    ipSize = p_rspaceTimeDiscr%p_rlevelInfo%rstaticInfo%rmatrixB1%NCOLS
   
    ! Create a global matrix:
    call lsysbl_createEmptyMatrix (rmatrix,6*(p_rspaceTimeDiscr%NEQtime))

    ! Initialise the boundary conditions
    call bcasm_initDiscreteBC(rdiscreteBC)
    call bcasm_initDiscreteFBC(rdiscreteFBC)
  
    ! Get a temporary system matrix. Use the default preconditioner specials, that's
    ! enough for us.
    call cc_allocPrecSystemMatrix (rproblem,rprecSpecials,&
        p_rspaceTimeDiscr%ilevel,rproblem%nlmin,rproblem%nlmax,&
        p_rspaceTimeDiscr%p_rlevelInfo,CCMASM_MTP_AUTOMATIC,&
        rblockTemp)

    ! Attach the boundary conditions
    call lsysbl_assignDiscreteBC(rblockTemp,rdiscreteBC)
    call lsysbl_assignDiscreteFBC(rblockTemp,rdiscreteFBC)

    ! Basic initialisation of our nonlinear matrix in space - for defect correction
    ! and matrix assembly in order to be used in a spatial preconditioner.
    call cc_initNonlinMatrix (rnonlinearSpatialMatrix,rproblem,&
        p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
        p_rspaceTimeDiscr%p_rlevelInfo%rstaticInfo)

    ! Create a vector for evaluating the nonlinearity in the previous, current and
    ! next timestep.
    ! (The vector is always created but only used when there is a nonlinearity!)
    call lsysbl_createVecBlockIndMat(rblockTemp,rvector1)
    call lsysbl_createVecBlockIndMat(rblockTemp,rvector2)
    call lsysbl_createVecBlockIndMat(rblockTemp,rvector3)
    
    call lsysbl_getbase_double (rvector1,p_Ddata1)
    call lsysbl_getbase_double (rvector2,p_Ddata2)
    call lsysbl_getbase_double (rvector3,p_Ddata3)
    
    ! Put the 'current' solution to the 'middle' vector.
    if (associated(rsupermatrix%p_rsolution))  then
      call sptivec_getTimestepData (rsupermatrix%p_rsolution, 1, rvector2)
    end if
    
    npureDirichletTimesteps = 0
    
    ! Loop through the substeps
    do isubstep = 0,p_rspaceTimeDiscr%NEQtime-1

      ! Load the data of the 'next' timestep to rvector3.    
      if (associated(rsupermatrix%p_rsolution) .and. &
          (isubstep .ne. p_rspaceTimeDiscr%NEQtime-1))  then
        call sptivec_getTimestepData (rsupermatrix%p_rsolution, 1+isubstep+1, rvector3)
      end if
    
      ! Current point in time
      dtime = &
          rproblem%rtimedependence%dtimeInit + isubstep*p_rspaceTimeDiscr%rtimeDiscr%dtstep

      ! -----
      ! Discretise the boundary conditions at the new point in time.
      call bcasm_clearDiscreteBC(rdiscreteBC)
      call bcasm_clearDiscreteFBC(rdiscreteFBC)
      call cc_assembleBDconditions (rproblem,dtime,&
          p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
          CCDISCBC_PRIMALDUAL,rdiscreteBC,rproblem%rcollection,bhasNeumann)
      call cc_assembleFBDconditions (rproblem,dtime,&
          p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
          CCDISCBC_PRIMALDUAL,rdiscreteFBC,rproblem%rcollection)
      
      ! Assemble diagonal blocks as well as the band above and below the diagonal.
      ileft = -1
      iright = 1
      if (isubstep .eq. 0) ileft = 0
      if (isubstep .eq. p_rspaceTimeDiscr%NEQtime-1) iright = 0
      
      ! Loop over the matrix bands in the current row isubstep
      do ix = ileft,iright
      
        ! Set up the matrix weights of that submatrix.
        call cc_initNonlinMatrix (rnonlinearSpatialMatrix,rproblem,&
            p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
            p_rspaceTimeDiscr%p_rlevelInfo%rstaticInfo)
        call cc_setupMatrixWeights (rproblem,rsupermatrix,&
          p_rspaceTimeDiscr%rtimeDiscr%dtheta,isubstep,ix,rnonlinearSpatialMatrix)
          
        ! Assemble the matrix in rblockTemp.
        ! Note that the weights of the matrices must be set before, otherwise
        ! the assembly routines would complain about missing matrices :-)
        call smva_assembleMatrix (CCMASM_COMPUTE,CCMASM_MTP_AUTOMATIC,&
            rblockTemp,rnonlinearSpatialMatrix,rvector1,rvector2,rvector3) 

        ! Switch of matrices that aren't needed.
        select case (ix)
        case (-1)
          ! Specify the matrix as 'off-diagonal' matrix because it's not on the
          ! main diagonal of the supermatrix.
          rblockTemp%imatrixSpec = LSYSBS_MSPEC_OFFDIAGSUBMATRIX
          
          ! Mass matrices in the dual equation not needed since the band
          ! below the diagonal couples only the timesteps of the primal equation.
          ! This saves some memory as the matrices are =0 anyway.
!          rblockTemp%RmatrixBlock(1,1)%dscaleFactor = 1.0_DP
!          rblockTemp%RmatrixBlock(2,2)%dscaleFactor = 1.0_DP
!          rblockTemp%RmatrixBlock(4,4)%dscaleFactor = 0.0_DP
!          rblockTemp%RmatrixBlock(5,5)%dscaleFactor = 0.0_DP
        case (1)
          ! Specify the matrix as 'off-diagonal' matrix because it's not on the
          ! main diagonal of the supermatrix.
          rblockTemp%imatrixSpec = LSYSBS_MSPEC_OFFDIAGSUBMATRIX
          
          ! Mass matrices in the primal equation not needed since the band
          ! above the diagonal couples only the timesteps of the dual equation.
          ! This saves some memory as the matrices are =0 anyway.
!          rblockTemp%RmatrixBlock(1,1)%dscaleFactor = 0.0_DP
!          rblockTemp%RmatrixBlock(2,2)%dscaleFactor = 0.0_DP
!          rblockTemp%RmatrixBlock(4,4)%dscaleFactor = 1.0_DP
!          rblockTemp%RmatrixBlock(5,5)%dscaleFactor = 1.0_DP
        case DEFAULT
          rblockTemp%imatrixSpec = LSYSBS_MSPEC_GENERAL
          
!          rblockTemp%RmatrixBlock(1,1)%dscaleFactor = 1.0_DP
!          rblockTemp%RmatrixBlock(2,2)%dscaleFactor = 1.0_DP
!          rblockTemp%RmatrixBlock(4,4)%dscaleFactor = 1.0_DP
!          rblockTemp%RmatrixBlock(5,5)%dscaleFactor = 1.0_DP
        end select

        ! Include the boundary conditions into that matrix.
        call matfil_discreteBC (rblockTemp,rdiscreteBC)
        call matfil_discreteFBC (rblockTemp,rdiscreteFBC)

        ! Include the current matrix into the global matrix 
        call insertMatrix (rblockTemp,rmatrix,(isubstep+ix)*6+1,isubstep*6+1,.false.)
        
        ! DEBUG!!!
        !IF (isubstep .EQ. p_rspaceTimeDiscr%niterations) THEN
        !  CALL matio_writeBlockMatrixHR (rblockTemp, 'matrix',&
        !    .TRUE., 0, 'matrix'//TRIM(sys_siL(ix,10))//'.txt', '(E20.10)')
        !END IF
        
      end do  

      if (.not. bhasNeumann) then
        ! Insert a 'point matrix' containing a zero in the pressure block.
        ! This allows the boundary condition implementation routine to
        ! insert a unit vector for the pressure if we have a pure-Dirichlet
        ! problem.
        ! Don't insert the matrix if there is already one! An existing matrix
        ! is always an identity matrix, as it happens to appear in the first
        ! and last time strep.
        ! IF (isubstep .NE. 0) THEN
        if (.not. lsysbl_isSubmatrixPresent (rmatrix,isubstep*6+3,isubstep*6+3)) then
          call createPointMatrix (rmatrix%RmatrixBlock(isubstep*6+3,isubstep*6+3),&
              p_rspaceTimeDiscr%p_rlevelInfo%rstaticInfo%rmatrixTemplateFEMPressure%NEQ,1)
        end if
        if (.not. lsysbl_isSubmatrixPresent (rmatrix,isubstep*6+6,isubstep*6+6)) then
          call createPointMatrix (rmatrix%RmatrixBlock(isubstep*6+6,isubstep*6+6),&
              p_rspaceTimeDiscr%p_rlevelInfo%rstaticInfo%rmatrixTemplateFEMPressure%NEQ,1)
        end if
        
        ! Remember this timestep as "pure Dirichlet".
        npureDirichletTimesteps = npureDirichletTimesteps + 1
        IpureDirichletTimesteps(npureDirichletTimesteps) = isubstep+1
        
      end if
      
      ! Cycle the solution vectors: 1 <- 2 <- 3
      if (associated(rsupermatrix%p_rsolution)) then
        call lsysbl_copyVector (rvector2,rvector1)
        call lsysbl_copyVector (rvector3,rvector2)
        ! Now rvector3 is free and can take the data of the next timestep.
      end if
      
    end do
    
    ! Release the temp vectors
    call lsysbl_releaseVector (rvector3)
    call lsysbl_releaseVector (rvector2)
    call lsysbl_releaseVector (rvector1)
    
    ! Release the temp matrix
    call lsysbl_releaseMatrix (rblockTemp)

    ! Release the boundary conditions
    call bcasm_releaseDiscreteBC(rdiscreteBC)
    call bcasm_releaseDiscreteFBC(rdiscreteFBC)

    ! Update structural information of the global matrix.
    call lsysbl_updateMatStrucInfo (rmatrix)
  
    ! Clean up the collection (as we are done with the assembly, that's it.
    call cc_doneCollectForAssembly (rproblem,rproblem%rcollection)

  contains

    ! ---------------------------------------------------------------
    
    subroutine insertMatrix (rsource,rdest,ileft,itop,bcopyStructure)
    
    ! Includes rsource into rdest at position ileft,itop
    type(t_matrixBlock), intent(IN) :: rsource
    type(t_matrixBlock), intent(INOUT) :: rdest
    integer, intent(IN) :: ileft
    integer, intent(IN) :: itop
    ! Duplicate the structure. FALSE=share the structure.
    logical, intent(IN) :: bcopyStructure 
    
    integer :: i,j,ccopy
    
    ccopy = LSYSSC_DUP_SHARE
    if (bcopyStructure) ccopy = LSYSSC_DUP_COPY
    
    do j=1,rsource%nblocksPerRow
      do i=1,rsource%nblocksPerCol
        if (lsysbl_isSubmatrixPresent (rsource,i,j)) then
          if (lsysbl_isSubmatrixPresent (rdest,i+itop-1,j+ileft-1)) then
            call lsyssc_releaseMatrix (rdest%RmatrixBlock(j+ileft-1,i+itop-1))
          end if
          call lsyssc_duplicateMatrix (rsource%RmatrixBlock(i,j),&
              rdest%RmatrixBlock(i+itop-1,j+ileft-1),&
              ccopy,LSYSSC_DUP_COPY)
        end if
      end do
    end do
        
    end subroutine

    ! ---------------------------------------------------------------
    
    subroutine createSumMatrix (rmatrix,NCOLS,NEQ,irow)
    
    ! Creates a structure-9 matrix that contains one row, filled by 'ones'.
    ! This corresponds to an additional equation that sums up all
    ! vector entries, which corresponds to setting up the integral
    ! of a function.
    ! If the matrix does not exist, it's created.
    ! If the matrix exists, it's overwritten.
    
    ! Matrix to be set up
    type(t_matrixScalar), intent(INOUT) :: rmatrix
    
    ! Number of rows/columns in the matrix
    integer, intent(IN) :: NCOLS
    integer, intent(IN) :: NEQ
    
    ! Number of the row that should contain the line '1...1'.
    integer, intent(IN) :: irow
    
      ! local variables
      integer, dimension(:), pointer :: p_Kcol
      integer, dimension(:), pointer :: p_Kld
      integer, dimension(:), pointer :: p_Kdiagonal
      real(DP), dimension(:), pointer :: p_Ddata
      integer :: i
      
      ! The matrix should get the shape:
      !   1 1 1 1 ... 1 1 1
      !   0 0 0 0 ... 0 0 0      
      !   ...
      !   0 0 0 0 ... 0 0 0      
      !
      ! Create the matrix by hand. If necessary, allocate memory.
      if (rmatrix%h_Kld .eq. ST_NOHANDLE) then
        call storage_new ('createSumMatrix', 'Kld', NEQ+1, &
                          ST_INT, rmatrix%h_Kld,ST_NEWBLOCK_NOINIT)
        call storage_new ('createSumMatrix', 'Kdiagonal', NEQ+1, &
                          ST_INT, rmatrix%h_Kdiagonal,ST_NEWBLOCK_NOINIT)
        call storage_new ('createSumMatrix', 'Kcol', NCOLS, &
                          ST_INT, rmatrix%h_Kcol,ST_NEWBLOCK_NOINIT)
                          
        ! DA has probably a zero in front and at the end as dummy
        ! for diagonal entries.
        call storage_new ('createSumMatrix', 'Da', NCOLS, &
                          ST_DOUBLE, rmatrix%h_Da,ST_NEWBLOCK_NOINIT)

        ! Initialise matrix format and other parameters, that's it.
        rmatrix%NEQ = NEQ
        rmatrix%NCOLS = NCOLS
        rmatrix%cmatrixFormat = LSYSSC_MATRIX9
        rmatrix%NA = NCOLS
      
      end if
      
      call lsyssc_getbase_double (rmatrix,p_Ddata)
      call lsyssc_getbase_Kcol (rmatrix,p_Kcol)
      call lsyssc_getbase_Kld (rmatrix,p_Kld)
      call lsyssc_getbase_Kdiagonal (rmatrix,p_Kdiagonal)
      
      ! irow'th row
      p_Ddata(:) = 1.0_DP
      
      p_Kld(1:irow) = 1
      p_Kld(irow+1:) = NCOLS+1
      
      p_Kdiagonal(1:irow-1) = 1
      p_Kdiagonal(irow) = irow
      
      do i=1,NCOLS
        p_Kcol(i) = i
      end do

    end subroutine

    ! ---------------------------------------------------------------
    
    subroutine createPointMatrix (rmatrix,NEQ,irow)
    
    ! Creates a structure-9 matrix that contains exactly one 'zero'
    ! on the diagonal of row irow. 
    
    ! Matrix to be set up
    type(t_matrixScalar), intent(INOUT) :: rmatrix
    
    ! Number of rows/columns in the matrix
    integer, intent(IN) :: NEQ
    
    ! Number of the row that should contain the 'one'
    integer, intent(IN) :: irow
    
      ! local variables
      integer, dimension(:), pointer :: p_Kcol
      integer, dimension(:), pointer :: p_Kld
      integer, dimension(:), pointer :: p_Kdiagonal
      real(DP), dimension(:), pointer :: p_Ddata
      integer :: i
      
      ! Create the matrix by hand. If necessary, allocate memory.
      if (rmatrix%h_Kld .eq. ST_NOHANDLE) then
        call storage_new ('createPointMatrix', 'Kld', NEQ+1, &
                          ST_INT, rmatrix%h_Kld,ST_NEWBLOCK_NOINIT)
        call storage_new ('createPointMatrix', 'Kdiagonal', NEQ, &
                          ST_INT, rmatrix%h_Kdiagonal,ST_NEWBLOCK_NOINIT)
        call storage_new ('createPointMatrix', 'Kcol', 1, &
                          ST_INT, rmatrix%h_Kcol,ST_NEWBLOCK_NOINIT)
                          
        ! DA has probably a zero in front and at the end as dummy
        ! for diagonal entries.
        call storage_new ('createPointMatrix', 'Da', 1, &
                          ST_DOUBLE, rmatrix%h_Da,ST_NEWBLOCK_NOINIT)

        ! Initialise matrix format and other parameters, that's it.
        rmatrix%NEQ = NEQ
        rmatrix%NCOLS = NEQ
        rmatrix%cmatrixFormat = LSYSSC_MATRIX9
        rmatrix%NA = 1
        rmatrix%dscaleFactor = 1.0_DP
      
      end if
      
      call lsyssc_getbase_double (rmatrix,p_Ddata)
      call lsyssc_getbase_Kcol (rmatrix,p_Kcol)
      call lsyssc_getbase_Kld (rmatrix,p_Kld)
      call lsyssc_getbase_Kdiagonal (rmatrix,p_Kdiagonal)
      
      ! irow'th row
      p_Ddata(1) = 0.0_DP
      p_Kcol(1) = irow

      p_Kld(1:irow) = 1
      p_Kld(irow+1:) = 2
      
      p_Kdiagonal(:) = 1
      
    end subroutine
    
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine sptils_doneDataUMFPACK4 (rsolverNode)
  
!<description>
  ! Releases the memory of the numeric factorisation of the given matrix.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure of the UMFPACK4 solver
  type(t_sptilsNode), intent(INOUT)         :: rsolverNode
!</inputoutput>
  
!</subroutine>

  ! Release the numerical factorisation if associated
  if (rsolverNode%p_rsubnodeUMFPACK4%inumeric .ne. 0) then
    call UMF4FNUM(rsolverNode%p_rsubnodeUMFPACK4%inumeric)
    rsolverNode%p_rsubnodeUMFPACK4%inumeric = 0
  end if  
  
  ! Release the symbolical factorisation if associated
  if (rsolverNode%p_rsubnodeUMFPACK4%isymbolic .ne. 0) then
    call UMF4FSYM(rsolverNode%p_rsubnodeUMFPACK4%isymbolic)
    rsolverNode%p_rsubnodeUMFPACK4%isymbolic = 0
  end if
  
end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine sptils_doneStructureUMFPACK4 (rsolverNode)
  
!<description>
  ! Releases the memory of the numeric factorisation of the given matrix.
!</description>
  
!<inputoutput>
  ! The t_linsolNode structure of the UMFPACK4 solver
  type(t_sptilsNode), intent(INOUT)         :: rsolverNode
!</inputoutput>
  
!</subroutine>

  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine sptils_doneUMFPACK4 (rsolverNode)
  
!<description>
  ! This routine releases all temporary memory for the UMFPACK4 solver from
  ! the heap.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure of UMFPACK4 which is to be cleaned up.
  type(t_sptilsNode), intent(INOUT)         :: rsolverNode
!</inputoutput>
  
!</subroutine>
  
  ! Release symbolical and numerical factorisation if still associated...
  call sptils_doneDataUMFPACK4 (rsolverNode)
  call sptils_doneStructureUMFPACK4 (rsolverNode)
  
  deallocate(rsolverNode%p_rsubnodeUMFPACK4)
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine sptils_precUMFPACK4 (rsolverNode,rd)
  
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
  type(t_sptilsNode), intent(INOUT)         :: rsolverNode

  ! On call to this routine: The defect vector to be preconditioned.
  ! Will be overwritten by the preconditioned defect.
  type(t_spacetimeVector), intent(INOUT)        :: rd
!</inputoutput>
  
!</subroutine>
 
  ! local variables
  integer :: KSYS
  real(DP), dimension(:), pointer :: p_Dx,p_Db
  type(t_vectorBlock) :: rb,rx
  type(t_sptilsSubnodeUMFPACK4), pointer :: p_rsubnode

  ! Status variables of UMFPACK4; receives the UMFPACK-specific return code
  ! of a call to the solver routines.
  real(DP), dimension(90) :: Dinfo
  
!  type(t_matrixBlock) :: rmatrixGlobal,rmatrixGlobalSc
!  type(t_vectorBlock) :: rrhs,rsol
!  TYPE(t_spacetimeVector) :: rsol2
!  real(dp) :: dnorm
  
    ! Status reset
    rsolverNode%iresult = 0
    
    ! Getch some information
    p_rsubnode => rsolverNode%p_rsubnodeUMFPACK4
    
    ! Copy the RHS rd to the temp vector; it will be overwritten
    ! by the solution vector
    call sptivec_convertSupervecToVector (rd,rb)

    ! Get the array    
    call lsysbl_getbase_double(rb,p_Db)

    ! Allocate memory for the destination
    call sptivec_convertSupervecToVector (rd,rx)

    ! Get the RHS and solution vector data
    call lsysbl_getbase_double(rx,p_Dx)

    ! Solve the system
    ! Solve the system. Note that UMFPACK expects the matrix in
    ! CSR format, which is transposed to our matrix format 9 --
    ! So we solve the transposed system:
    KSYS = 1
    call UMF4SOL(KSYS,p_Dx,p_Db,rsolverNode%p_rsubnodeUMFPACK4%inumeric,&
                 rsolverNode%p_rsubnodeUMFPACK4%Dcontrol,Dinfo)

!    !!! DEBUG
!    CALL assembleGlobalSpaceTimeMatrix (rsolverNode%p_rproblem,&
!        rsolverNode%rmatrix,rmatrixGlobal)
!    call lsysbl_createVecBlockDirect (rrhs, (/rb%NEQ/),.false.)
!    call lsysbl_createVecBlockDirect (rsol, (/rb%NEQ/),.false.)
!    CALL lsysbl_getbase_double(rrhs,p_Db)
!    call storage_copy (rb%h_Ddata,rrhs%h_Ddata)
!    call storage_copy (rx%h_Ddata,rsol%h_Ddata)
!    call glsys_assembleGlobal (rmatrixGlobal,rmatrixGlobalSc, .true., .true.)
!    call lsysbl_blockMatVec (rmatrixGlobalSc,rsol,rrhs,-1.0_DP,1.0_DP)
!    call lsysbl_releaseMatrix (rmatrixGlobal)
!    call lsysbl_releaseMatrix (rmatrixGlobalSc)
!    call lsysbl_releaseVector (rrhs)
!    
!    call sptivec_copyVector (rd,rsol2)
!    CALL sptivec_convertVectorToSupervec (rx,rsol2)
!    call cc_spaceTimeMatVec (rsolverNode%p_rproblem, rsolverNode%rmatrix, &
!      rsol2, rd, -1.0_DP, 1.0_DP, SPTID_FILTER_DEFECT, dnorm,.true.)
!    call sptivec_releaseVector (rsol2)


    ! Save the result
    call sptivec_convertVectorToSupervec (rx,rd)

    ! Scale by omega
    call sptivec_scaleVector (rd,rsolverNode%domega)
                
    ! Release temp vectors
    call lsysbl_releaseVector (rx)
    call lsysbl_releaseVector (rb)
  
    ! Check the solver status
    select case (int(Dinfo(1)))
    case (0) 
      ! All ok.
      rsolverNode%iiterations = 1
      rsolverNode%dconvergenceRate = 0.0_DP
    case DEFAULT
      ! We had an error. Don't know which one.
      rsolverNode%iresult = -1
    end select
    
  end subroutine

! *****************************************************************************
! Routines for the BiCGStab solver
! *****************************************************************************

!<subroutine>
  
  subroutine sptils_initBiCGStab (rproblem,p_rsolverNode,p_rpreconditioner)
  
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
  type(t_problem), target :: rproblem

  ! OPTIONAL: A pointer to the solver structure of a solver that should be 
  ! used for preconditioning. If not given or set to NULL(), no preconditioning 
  ! will be used.
  type(t_sptilsNode), pointer, optional   :: p_rpreconditioner
!</input>
  
!<output>
  ! A pointer to a t_sptilsNode structure. Is set by the routine, any previous
  ! value of the pointer is destroyed.
  type(t_sptilsNode), pointer         :: p_rsolverNode
!</output>
  
!</subroutine>

    ! Create a default solver structure
    call sptils_initSolverGeneral(rproblem,p_rsolverNode)
    
    ! Initialise the type of the solver
    p_rsolverNode%calgorithm = SPTILS_ALG_BICGSTAB
    
    ! Initialise the ability bitfield with the ability of this solver:
    p_rsolverNode%ccapability = SPTILS_ABIL_CHECKDEF + &
                                SPTILS_ABIL_USESUBSOLVER
    
    ! Allocate the subnode for BiCGStab.
    ! This initialises most of the variables with default values appropriate
    ! to this solver.
    allocate(p_rsolverNode%p_rsubnodeBiCGStab)
    
    ! Attach the preconditioner if given. 
    
    if (present(p_rpreconditioner)) then 
      p_rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner => p_rpreconditioner
    end if

  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine sptils_setMatrixBiCGStab (rsolverNode,Rmatrices)
  
!<description>
  
  ! This routine is called if the system matrix changes.
  ! The routine calls sptils_setMatrices for the preconditioner of BiCGStab
  ! to inform also that one about the change of the matrix pointer.
  
!</description>
  
!<input>
  ! An array of system matrices which is simply passed to the initialisation 
  ! routine of the preconditioner.
  type(t_ccoptSpaceTimeMatrix), dimension(:), intent(IN)   :: Rmatrices
!</input>
  
!<inputoutput>
  ! The t_sptilsNode structure of the BiCGStab solver
  type(t_sptilsNode), intent(INOUT)         :: rsolverNode
!</inputoutput>
  
!</subroutine>

    ! Do we have a preconditioner given? If yes, call the general initialisation
    ! routine to initialise it.
    
    if (associated(rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner)) then
      call sptils_setMatrices (rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner, &
                              Rmatrices)
    end if

  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine sptils_initStructureBiCGStab (rsolverNode, ierror)
  
!<description>
  ! Calls the initStructure subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_initStructure.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure of the UMFPACK4 solver
  type(t_sptilsNode), intent(INOUT)         :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to 
  ! SPTILS_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  integer, intent(OUT) :: ierror
!</output>

!</subroutine>

    ! local variables
    integer :: isubgroup,i,NEQtime
    type(t_sptilsSubnodeBiCGStab), pointer :: p_rsubnode
    
    ! A-priori we have no error...
    ierror = SPTILS_ERR_NOERROR

    ! We simply pass isubgroup to the subsolvers when calling them.
    ! Inside of this routine, there's not much to do with isubgroup,
    ! as there is no special information (like a factorisation)
    ! associated with this type of solver, which has to be allocated,
    ! released or updated...

    ! Call the init routine of the preconditioner.
    if (associated(rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner)) then
      call sptils_initStructure (rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner, &
                                isubgroup)
    end if
    
    ! BiCGStab needs 5 temporary vectors + 1 for preconditioning. 
    ! Allocate that here! Use the default data type prescribed in the solver 
    ! structure for allocating the temp vectors.
    p_rsubnode => rsolverNode%p_rsubnodeBiCGStab
    NEQtime = rsolverNode%rmatrix%p_rspaceTimeDiscr%NEQtime
    do i=1,6
      call sptivec_initVector (p_rsubnode%RtempVectors(i),NEQtime,&
        rsolverNode%rmatrix%p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation)
    end do
    
    ! Allocate memory for a spatial temp vector.
    call lsysbl_createVecBlockByDiscr (&
        rsolverNode%rmatrix%p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
        rsolverNode%p_rsubnodeBiCGStab%rtempVectorSpace)
    
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine sptils_initDataBiCGStab (rsolverNode, ierror)
  
!<description>
  ! Calls the initData subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_initData.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure of the UMFPACK4 solver
  type(t_sptilsNode), intent(INOUT)         :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the LINSOL_ERR_XXXX constants. A value different to 
  ! SPTILS_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  integer, intent(OUT) :: ierror
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
    if (associated(rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner)) then
      call sptils_initData (rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner, ierror)
    end if
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine sptils_doneDataBiCGStab (rsolverNode, isolverSubgroup)
  
!<description>
  ! Calls the doneData subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_doneData.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure of the UMFPACK4 solver
  type(t_sptilsNode), intent(INOUT)         :: rsolverNode
!</inputoutput>
  
!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific 
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  integer, optional, intent(IN)                    :: isolverSubgroup
!</input>

!</subroutine>

    ! We simply pass isubgroup to the subsolvers when calling them.
    ! Inside of this routine, there's not much to do with isubgroup,
    ! as there is no special information (like a factorisation)
    ! associated with this type of solver, which has to be allocated,
    ! released or updated...

    ! Call the done routine of the preconditioner.
    if (associated(rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner)) then
      call sptils_doneData (rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner)
    end if
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine sptils_doneStructureBiCGStab (rsolverNode, isolverSubgroup)
  
!<description>
  ! Calls the doneStructure subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with linsol_doneStructure.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure of the UMFPACK4 solver
  type(t_sptilsNode), intent(INOUT)         :: rsolverNode
!</inputoutput>
  
!<input>
  ! Optional parameter. isolverSubgroup allows to specify a specific 
  ! subgroup of solvers in the solver tree to be processed. By default,
  ! all solvers in subgroup 0 (the default solver group) are processed,
  ! solvers in other solver subgroups are ignored.
  ! If isolverSubgroup != 0, only the solvers belonging to subgroup
  ! isolverSubgroup are processed.
  integer, optional, intent(IN)                    :: isolverSubgroup
!</input>

!</subroutine>

    ! local variables
    integer :: i
    type(t_sptilsSubnodeBiCGStab), pointer :: p_rsubnode
    
    ! We simply pass isubgroup to the subsolvers when calling them.
    ! Inside of this routine, there's not much to do with isubgroup,
    ! as there is no special information (like a factorisation)
    ! associated with this type of solver, which has to be allocated,
    ! released or updated...

    ! Call the init routine of the preconditioner.
    if (associated(rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner)) then
      call sptils_doneStructure (rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner)
    end if
    
    ! Release temporary data if associated
    p_rsubnode => rsolverNode%p_rsubnodeBiCGStab
    if (p_rsubnode%RtempVectors(1)%NEQ .ne. 0) then
      do i=6,1,-1
        call sptivec_releaseVector(p_rsubnode%RtempVectors(i))
      end do
    end if
    
    if (rsolverNode%p_rsubnodeBiCGStab%rtempVectorSpace%NEQ .ne. 0) then
      call lsysbl_releaseVector (rsolverNode%p_rsubnodeBiCGStab%rtempVectorSpace)
    end if
    
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine sptils_doneBiCGStab (rsolverNode)
  
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
  type(t_sptilsNode), pointer         :: rsolverNode
!</input>
  
!</subroutine>

    ! Check if there's a preconditioner attached. If yes, release it.
    if (associated(rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner)) then
      call sptils_releaseSolver(rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner)
    end if
    
    ! Release memory if still associated
    call sptils_doneDataBiCGStab (rsolverNode)
    call sptils_doneStructureBiCGStab (rsolverNode)
    
    ! Release the BiCGStab subnode
    deallocate(rsolverNode%p_rsubnodeBiCGStab)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  recursive subroutine sptils_precBiCGStab (rsolverNode,rd)
  
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
  type(t_sptilsNode), intent(INOUT), target :: rsolverNode
   
  ! On call to this routine: The defect vector to be preconditioned.
  ! Will be overwritten by the preconditioned defect.
  type(t_spaceTimeVector), intent(INOUT)        :: rd
!</inputoutput>
  
!</subroutine>

  ! local variables
  real(DP) :: dalpha,dbeta,domega0,domega1,domega2,dres
  real(DP) :: drho1,drho0,dfr,dresscale,dresunprec
  integer :: ite

  ! The system matrix
  type(t_ccoptSpaceTimeDiscretisation), pointer :: p_rspaceTimeDiscr
  type(t_ccoptSpaceTimeMatrix), pointer :: p_rmatrix
  
  ! Minimum number of iterations, print-sequence for residuals
  integer :: nminIterations, niteResOutput
  
  ! Whether to filter/prcondition
  logical bprec,bfilter
  
  ! Our structure
  type(t_sptilsSubnodeBiCGStab), pointer :: p_rsubnode
  
  ! Pointers to temporary vectors - named for easier access
  type(t_spaceTimeVector), pointer :: p_DR,p_DR0,p_DP,p_DPA,p_DSA,p_rx
  type(t_sptilsNode), pointer :: p_rprecSubnode
  
    ! Solve the system!
  
    ! Status reset
    rsolverNode%iresult = 0
    call stat_clearTimer (rsolverNode%rtimeSpacePrecond)
    
    ! Getch some information
    p_rsubnode => rsolverNode%p_rsubnodeBiCGStab
    p_rspaceTimeDiscr => rsolverNode%rmatrix%p_rspaceTimeDiscr
    p_rmatrix => rsolverNode%rmatrix

    ! Check the parameters
    if ((rd%NEQtime .eq. 0) .or. &
        (p_rspaceTimeDiscr%rtimeDiscr%nintervals .eq. 0)) then
    
      ! Parameters wrong
      rsolverNode%iresult = 2
      return
    end if

    ! Minimum number of iterations
 
    nminIterations = max(rsolverNode%nminIterations,0)
      
    ! Use preconditioning? Filtering?

    bprec = associated(rsolverNode%p_rsubnodeBiCGStab%p_rpreconditioner)
    
    ! Iteration when the residuum is printed:

    niteResOutput = max(1,rsolverNode%niteResOutput)

    ! Set pointers to the temporary vectors
    p_DR   => p_rsubnode%RtempVectors(1)
    p_DR0  => p_rsubnode%RtempVectors(2)
    p_DP   => p_rsubnode%RtempVectors(3)
    p_DPA  => p_rsubnode%RtempVectors(4)
    p_DSA  => p_rsubnode%RtempVectors(5)
    p_rx   => p_rsubnode%RtempVectors(6)
    
    if (bprec) then
      p_rprecSubnode => p_rsubnode%p_rpreconditioner
    end if
    
    ! rd is our RHS. p_rx points to a new vector which will be our
    ! iteration vector. At the end of this routine, we replace
    ! rd by p_rx.
    ! Clear our iteration vector p_rx.
    call sptivec_clearVector (p_rx)
      
    ! Initialize used vectors with zero
      
    call sptivec_clearVector(p_DP)
    call sptivec_clearVector(p_DPA)
    
    ! Initialise the iteration vector with zero.

    ! Initialization

    drho0  = 1.0_DP
    dalpha = 1.0_DP
    domega0 = 1.0_DP

    ! Copy our RHS rd to p_DR. As the iteration vector is 0, this
    ! is also our initial defect.

    call sptivec_copyVector(rd,p_DR)
    
    ! Filter the defect for boundary conditions in space and time.
    call tbc_implementInitCondDefect (&
        p_rspaceTimeDiscr,p_DR,p_rsubnode%rtempVectorSpace)
    call tbc_implementBCdefect (rsolverNode%p_rproblem,&
        p_rspaceTimeDiscr,p_DR,p_rsubnode%rtempVectorSpace)

    ! Get the norm of the residuum
    dres = sptivec_vectorNorm (p_DR,rsolverNode%iresNorm)
    if (.not.((dres .ge. 1E-99_DP) .and. &
              (dres .le. 1E99_DP))) dres = 0.0_DP

    if (bprec) then
      ! Perform preconditioning with the assigned preconditioning
      ! solver structure.
      call sptils_precondDefect (p_rprecSubnode,p_DR)
      call stat_addTimers (p_rprecSubnode%rtimeSpacePrecond,rsolverNode%rtimeSpacePrecond)
      
      ! We scane the absolute stopping criterion by the difference
      ! between the preconditioned and unpreconditioned defect --
      ! to encounter the difference in the residuals.
      ! This is of course an approximation to 
      dresunprec = dres
      dres = sptivec_vectorNorm (p_DR,rsolverNode%iresNorm)
      
      call output_line ('Space-Time-BiCGStab: Iteration '// &
          trim(sys_siL(0,10))//',  !!RES(unscaled)!! = '//&
          trim(sys_sdEL(dres,15)) )
      
      if (.not.((dres .ge. 1E-99_DP) .and. &
                (dres .le. 1E99_DP))) dres = 1.0_DP
      dresscale = dresunprec / dres
    else
      dresscale = 1.0_DP
    end if
    
    ! Initialize starting residuum
      
    rsolverNode%dinitialDefect = dres
    rsolverNode%dlastDefect = 0.0_DP

    ! Check if out initial defect is zero. This may happen if the filtering
    ! routine filters "everything out"!
    ! In that case we can directly stop our computation.

    if ( rsolverNode%dinitialDefect .lt. rsolverNode%drhsZero ) then
     
      ! final defect is 0, as initialised in the output variable above

      call sptivec_clearVector(p_rx)
      ite = 0
      rsolverNode%dfinalDefect = dres
          
    else

      if (rsolverNode%ioutputLevel .ge. 2) then
        if (bprec) then
          call output_line ('Space-Time-BiCGStab: Iteration '// &
              trim(sys_siL(0,10))//',  !!RES(scaled)!! = '//&
              trim(sys_sdEL(rsolverNode%dinitialDefect*dresscale,15)) )
        else
          call output_line ('Space-Time-BiCGStab: Iteration '// &
              trim(sys_siL(0,10))//',  !!RES!! = '//&
              trim(sys_sdEL(rsolverNode%dinitialDefect,15)) )
        end if
      end if

      call sptivec_copyVector(p_DR,p_DR0)

      ! Perform at most nmaxIterations loops to get a new vector

      do ite = 1,rsolverNode%nmaxIterations
      
        rsolverNode%icurrentIteration = ite

        drho1 = sptivec_scalarProduct (p_DR0,p_DR) 

        if (drho0*domega0 .eq. 0.0_DP) then
          ! Should not happen
          if (rsolverNode%ioutputLevel .ge. 2) then
            call output_line ('Space-Time-BiCGStab: Iteration prematurely stopped! '//&
                 'Correction vector is zero!')
          end if

          ! Some tuning for the output, then cancel.

          rsolverNode%iresult = -1
          rsolverNode%iiterations = ITE-1
          exit
          
        end if

        dbeta=(drho1*dalpha)/(drho0*domega0)
        drho0 = drho1

        call sptivec_vectorLinearComb (p_DR ,p_DP,1.0_DP,dbeta)
        call sptivec_vectorLinearComb (p_DPA ,p_DP,-dbeta*domega0,1.0_DP)

        call cc_spaceTimeMatVec (rsolverNode%p_rproblem, p_rmatrix, &
            p_DP,p_DPA, 1.0_DP,0.0_DP,SPTID_FILTER_NONE)
        
        ! Filter the defect for boundary conditions in space and time.
        call tbc_implementInitCondDefect (&
            p_rspaceTimeDiscr,p_DPA,p_rsubnode%rtempVectorSpace)
        call tbc_implementBCdefect (rsolverNode%p_rproblem,&
            p_rspaceTimeDiscr,p_DPA,p_rsubnode%rtempVectorSpace)

        if (bprec) then
          ! Perform preconditioning with the assigned preconditioning
          ! solver structure.
          call sptils_precondDefect (p_rprecSubnode,p_DPA)
          call stat_addTimers (p_rprecSubnode%rtimeSpacePrecond,rsolverNode%rtimeSpacePrecond)
        end if

        dalpha = sptivec_scalarProduct (p_DR0,p_DPA)
        
        if (dalpha .eq. 0.0_DP) then
          ! We are below machine exactness - we can't do anything more...
          ! May happen with very small problems with very few unknowns!
          if (rsolverNode%ioutputLevel .ge. 2) then
            call output_line ('Space-Time-BiCGStab: Convergence failed, ALPHA=0!')
            rsolverNode%iresult = -2
            exit
          end if
        end if
        
        dalpha = drho1/dalpha

        call sptivec_vectorLinearComb (p_DPA,p_DR,-dalpha,1.0_DP)

        call cc_spaceTimeMatVec (rsolverNode%p_rproblem, p_rmatrix, &
            p_DR,p_DSA, 1.0_DP,0.0_DP,SPTID_FILTER_NONE)
        
        ! Filter the defect for boundary conditions in space and time.
        call tbc_implementInitCondDefect (&
            p_rspaceTimeDiscr,p_DSA,p_rsubnode%rtempVectorSpace)
        call tbc_implementBCdefect (rsolverNode%p_rproblem,&
            p_rspaceTimeDiscr,p_DSA,p_rsubnode%rtempVectorSpace)

        if (bprec) then
          ! Perform preconditioning with the assigned preconditioning
          ! solver structure.
          call sptils_precondDefect (p_rprecSubnode,p_DSA)
          call stat_addTimers (p_rprecSubnode%rtimeSpacePrecond,rsolverNode%rtimeSpacePrecond)
        end if
        
        domega1 = sptivec_scalarProduct (p_DSA,p_DR)
        domega2 = sptivec_scalarProduct (p_DSA,p_DSA)
        
        if (domega1 .eq. 0.0_DP) then
          domega0 = 0.0_DP
        else
          if (domega2 .eq. 0.0_DP) then
            if (rsolverNode%ioutputLevel .ge. 2) then
              call output_line ('Space-Time-BiCGStab: Convergence failed: omega=0!')
              rsolverNode%iresult = -2
              exit
            end if
          end if
          domega0 = domega1/domega2
        end if

        call sptivec_vectorLinearComb (p_DP ,p_rx,dalpha,1.0_DP)
        call sptivec_vectorLinearComb (p_DR ,p_rx,domega0,1.0_DP)

        call sptivec_vectorLinearComb (p_DSA,p_DR,-domega0,1.0_DP)

        ! Get the norm of the new (final?) residuum
        dfr = sptivec_vectorNorm (p_DR,rsolverNode%iresNorm)
     
        rsolverNode%dlastDefect = rsolverNode%dfinalDefect
        rsolverNode%dfinalDefect = dfr

        ! Test if the iteration is diverged
        dfr = dfr*dresscale
        if (sptils_testDivergence(rsolverNode,dfr)) then
          call output_line ('Space-Time-BiCGStab: Solution diverging!')
          rsolverNode%iresult = 1
          exit
        end if
     
        ! At least perform nminIterations iterations
        if (ite .ge. nminIterations) then
        
          ! Check if the iteration converged
          if (sptils_testConvergence(rsolverNode,dfr)) exit
          
        end if

        ! print out the current residuum

        if ((rsolverNode%ioutputLevel .ge. 2) .and. &
            (mod(ite,niteResOutput).eq.0)) then
          if (bprec) then
            call output_line ('Space-Time-BiCGStab: Iteration '// &
                trim(sys_siL(ITE,10))//',  !!RES(scaled)!! = '//&
                trim(sys_sdEL(rsolverNode%dfinalDefect*dresscale,15)) )
          else
            call output_line ('Space-Time-BiCGStab: Iteration '// &
                trim(sys_siL(ITE,10))//',  !!RES!! = '//&
                trim(sys_sdEL(rsolverNode%dfinalDefect,15)) )
          end if
        end if

      end do

      ! Set ITE to NIT to prevent printing of "NIT+1" of the loop was
      ! completed

      if (ite .gt. rsolverNode%nmaxIterations) &
        ite = rsolverNode%nmaxIterations

      ! Finish - either with an error or if converged.
      ! Print the last residuum.

      if ((rsolverNode%ioutputLevel .ge. 2) .and. &
          (ite .ge. 1) .and. (ITE .lt. rsolverNode%nmaxIterations) .and. &
          (rsolverNode%iresult .ge. 0)) then
        if (bprec) then
          call output_line ('Space-Time-BiCGStab: Iteration '// &
              trim(sys_siL(ITE,10))//',  !!RES(scaled)!! = '//&
              trim(sys_sdEL(rsolverNode%dfinalDefect*dresscale,15)) )
        else
          call output_line ('Space-Time-BiCGStab: Iteration '// &
              trim(sys_siL(ITE,10))//',  !!RES!! = '//&
              trim(sys_sdEL(rsolverNode%dfinalDefect,15)) )
        end if
      end if

    end if

    rsolverNode%iiterations = ite
    
    ! Overwrite our previous RHS by the new correction vector p_rx.
    ! This completes the preconditioning.
    call sptivec_copyVector (p_rx,rd)
      
    ! Don't calculate anything if the final residuum is out of bounds -
    ! would result in NaN's,...
      
    if (rsolverNode%dfinalDefect .lt. 1E99_DP) then
    
      ! If the initial defect was zero, the solver immediately
      ! exits - and so the final residuum is zero and we performed
      ! no steps; so the resulting convergence rate stays zero.
      ! In the other case the convergence rate computes as
      ! (final defect/initial defect) ** 1/nit :

      if (rsolverNode%dfinalDefect .gt. rsolverNode%drhsZero) then
        rsolverNode%dconvergenceRate = &
                    (rsolverNode%dfinalDefect / rsolverNode%dinitialDefect) ** &
                    (1.0_DP/real(rsolverNode%iiterations,DP))
      end if

      if (rsolverNode%ioutputLevel .ge. 2) then
        call output_lbrk()
        call output_line ('Space-Time-BiCGStab statistics:')
        call output_lbrk()
        call output_line ('Iterations              : '//&
             trim(sys_siL(rsolverNode%iiterations,10)) )
        call output_line ('!!INITIAL RES!!         : '//&
             trim(sys_sdEL(rsolverNode%dinitialDefect,15)) )
        call output_line ('!!RES!!                 : '//&
             trim(sys_sdEL(rsolverNode%dfinalDefect,15)) )
        if (rsolverNode%dinitialDefect .gt. rsolverNode%drhsZero) then     
          call output_line ('!!RES!!/!!INITIAL RES!! : '//&
            trim(sys_sdEL(rsolverNode%dfinalDefect / rsolverNode%dinitialDefect,15)) )
        else
          call output_line ('!!RES!!/!!INITIAL RES!! : '//&
               trim(sys_sdEL(0.0_DP,15)) )
        end if
        call output_lbrk ()
        call output_line ('Rate of convergence     : '//&
             trim(sys_sdEL(rsolverNode%dconvergenceRate,15)) )

      end if

      if (rsolverNode%ioutputLevel .eq. 1) then
        call output_line (&
              'Space-Time-BiCGStab: Iterations/Rate of convergence: '//&
              trim(sys_siL(rsolverNode%iiterations,10))//' /'//&
              trim(sys_sdEL(rsolverNode%dconvergenceRate,15)) )
      end if
      
    else
      ! DEF=Infinity; RHO=Infinity, set to 1
      rsolverNode%dconvergenceRate = 1.0_DP
    end if  
  
  end subroutine
  
! *****************************************************************************
! Multigrid preconditioner
! *****************************************************************************

!<subroutine>
  
  recursive subroutine sptils_initMultigrid (rproblem,NLMIN,NLMAX,p_rsolverNode)
  
!<description>
  ! Creates a t_sptilsNode solver structure for the multigrid preconditioner.
!</description>
  
!<input>
  ! The problem structure that defines the problem.
  ! A pointer to this is saved in the solver node, so the problem structure
  ! must not be released before the solver node is released.
  type(t_problem), target :: rproblem
  
  ! Minimum level; level where to apply the coarse grid solver.
  ! Must be >= rproblem%NLMIN!
  integer, intent(IN) :: NLMIN
  
  ! Maximum level; level where the preconditioner is applied.
  ! Must be <= rproblem%NLMAX!
  integer, intent(IN) :: NLMAX
!</input>
  
!<output>
  ! A pointer to a t_sptilsNode structure. Is set by the routine, any previous
  ! value of the pointer is destroyed.
  type(t_sptilsNode), pointer         :: p_rsolverNode
!</output>
  
!</subroutine>
  
    ! Create a default solver structure
    call sptils_initSolverGeneral(rproblem,p_rsolverNode)
    
    ! Initialise the type of the solver
    p_rsolverNode%calgorithm = SPTILS_ALG_MULTIGRID
    
    ! Initialise the ability bitfield with the ability of this solver:
    p_rsolverNode%ccapability = SPTILS_ABIL_MULTILEVEL + SPTILS_ABIL_CHECKDEF     + &
                                SPTILS_ABIL_USESUBSOLVER

    ! Allocate a subnode for our solver.
    ! This initialises most of the variables with default values appropriate
    ! to this solver.
    allocate(p_rsolverNode%p_rsubnodeMultigrid)
    
    ! Set parameters and allocate memory for the level information.
    p_rsolverNode%p_rsubnodeMultigrid%NLMIN = NLMIN
    p_rsolverNode%p_rsubnodeMultigrid%NLMAX = NLMAX
    allocate(p_rsolverNode%p_rsubnodeMultigrid%p_Rlevels(1:NLMAX-NLMIN+1))
    
  end subroutine

  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine sptils_doneMultigrid (rsolverNode)
  
!<description>
  ! This routine releases all temporary memory for the efect correction solver 
  ! from the heap.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure which is to be cleaned up.
  type(t_sptilsNode), intent(INOUT)         :: rsolverNode
!</inputoutput>
  
!</subroutine>
    integer :: ilev
  
    ! Release memory if still associated
    call sptils_doneDataMultigrid (rsolverNode)
    call sptils_doneStructureMultigrid (rsolverNode)
    
    ! Release level information
    do ilev=rsolverNode%p_rsubnodeMultigrid%NLMAX,rsolverNode%p_rsubnodeMultigrid%NLMIN,-1
    
      ! Pre- and postsmoother may be identical; release them only once!
      if (associated(rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%p_rpresmoother)) then
        if (.not. associated(&
            rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%p_rpresmoother,&
            rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%p_rpostsmoother)) then
          call sptils_releaseSolver(&
              rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%p_rpresmoother)
        else
          nullify(rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%p_rpresmoother)
        end if
      end if

      if (associated(&
          rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%p_rpostsmoother)) then
        call sptils_releaseSolver(&
            rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%p_rpostsmoother)
      end if
      
      if (associated(&
          rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%p_rcoarseGridSolver)) then
        call sptils_releaseSolver(&
            rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%p_rcoarseGridSolver)
      end if
      
    end do
    
    ! Release the subnode structure
    deallocate(rsolverNode%p_rsubnodeMultigrid%p_Rlevels)
    deallocate(rsolverNode%p_rsubnodeMultigrid)
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine sptils_setMatrixMultigrid (rsolverNode,Rmatrices)
  
!<description>
  ! This routine is called if the system matrix changes.
  ! The routine calls sptils_setMatrices for the preconditioner
  ! to inform also that one about the change of the matrix pointer.
!</description>
  
!<input>
  ! An array of system matrices which is simply passed to the initialisation 
  ! routine of the preconditioner.
  type(t_ccoptSpaceTimeMatrix), dimension(:), intent(IN) :: Rmatrices
!</input>
  
!<inputoutput>
  ! The t_sptilsNode structure of the solver
  type(t_sptilsNode), intent(INOUT)         :: rsolverNode
!</inputoutput>
  
!</subroutine>

    integer :: ilev

    ! Loop through the level. 
    do ilev=1,size(Rmatrices)
      ! On each level, set the matrix.
      rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%rmatrix = &
          Rmatrices(ilev)
      
      ! Call the setmatrices-routine for the presmoother/postsmoother/
      ! coarse grid solver
      if (.not. associated(&
          rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%p_rpresmoother,&
          rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%p_rpostsmoother)) then
        ! May be the presmoother and postsmoother are identical; release them
        ! only once!
        if (associated(rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%p_rpresmoother)) then
          call sptils_setMatrices(&
              rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%p_rpresmoother,&
              Rmatrices(1:ilev))
        end if
      end if
      
      if (associated(rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%p_rpostsmoother)) then
        call sptils_setMatrices(&
            rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%p_rpostsmoother,&
            Rmatrices(1:ilev))
      end if
      
      if (associated(rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%p_rcoarseGridSolver)) then
        call sptils_setMatrices(&
            rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%p_rcoarseGridSolver,&
            Rmatrices(1:ilev))
      end if
    end do
    
  end subroutine

! ***************************************************************************

!<subroutine>
  
  recursive subroutine sptils_initStructureMultigrid (rsolverNode,ierror)
  
!<description>
  ! Solver preparation. Perform symbolic factorisation (not of the defect
  ! correcion solver, but of subsolvers). Allocate temporary memory.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure 
  type(t_sptilsNode), intent(INOUT)         :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the SPTILS_ERR_XXXX constants. A value different to 
  ! SPTILS_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  integer, intent(OUT) :: ierror
!</output>
  
!</subroutine>

    ! local variables
    integer :: ilev,NLMAX,NEQtime
    type(t_sptilsMGLevelInfo), pointer :: p_rmgLevel
    
    ! A-priori we have no error...
    ierror = SPTILS_ERR_NOERROR
    
    ! Maximum time level.
    NLMAX = size(rsolverNode%p_rsubnodeMultigrid%p_Rlevels)

    ! On each level, call the initStructure routine of the
    ! presmoother/postsmoother/coarse grid solver

    ! Loop through the level. 
    do ilev=1,NLMAX
    
      p_rmgLevel => rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)
    
      ! Call the setmatrices-routine for the presmoother/postsmoother/
      ! coarse grid solver
      if (.not. associated(p_rmgLevel%p_rpresmoother,p_rmgLevel%p_rpostsmoother)) then
        ! May be the presmoother and postsmoother are identical; initialise them
        ! only once!
        if (associated(p_rmgLevel%p_rpresmoother)) then
          call sptils_initStructure(p_rmgLevel%p_rpresmoother,ierror)
          if (ierror .ne. ierror) return
        end if
      end if
      
      if (associated(p_rmgLevel%p_rpostsmoother)) then
        call sptils_initStructure(p_rmgLevel%p_rpostsmoother,ierror)
        if (ierror .ne. ierror) return
      end if
      
      if (associated(p_rmgLevel%p_rcoarseGridSolver)) then
        call sptils_initStructure(p_rmgLevel%p_rcoarseGridSolver,ierror)
        if (ierror .ne. ierror) return
      end if
      
      ! Generate an interlevel projection structure for that level
      !CALL sptipr_initProjection (&
      !    rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%rinterlevelProjection,&
      !    rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)%&
      !        rspaceTimeDiscr%p_rlevelInfo%p_rdiscretisation)
              
      NEQtime = p_rmgLevel%rmatrix%p_rspaceTimeDiscr%NEQtime
              
      ! On all levels except for the maximum one, create a solution vector
      if (ilev .lt. NLMAX) then
        call sptivec_initVector (&
            p_rmgLevel%rsolutionVector,NEQtime,p_rmgLevel%&
            rmatrix%p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation)
      end if
      
      ! On all levels except for the first one, create a RHS and a temp vector
      if (ilev .gt. 1) then
        call sptivec_initVector (&
            p_rmgLevel%rrhsVector,NEQtime,p_rmgLevel%&
            rmatrix%p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation)
        call sptivec_initVector (&
            p_rmgLevel%rtempVector,NEQtime,p_rmgLevel%&
            rmatrix%p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation)
      end if
      
      ! If adaptive coarse grid correction is activated, we need a temporary
      ! vector for the adaptive coarse grid correction on every level.
      if (rsolverNode%p_rsubnodeMultigrid%dalphamin .ne. &
          rsolverNode%p_rsubnodeMultigrid%dalphamax) then
        call sptivec_initVector (&
            p_rmgLevel%rtempCGCvector,NEQtime,p_rmgLevel%&
            rmatrix%p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation)
      end if
      
      ! Create a block temp vector for the interlevel projection
      call lsysbl_createVecBlockByDiscr(&
          p_rmgLevel%rmatrix%p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
          p_rmgLevel%rprjVector,.false.)

    end do
    
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine sptils_initDataMultigrid (rsolverNode, ierror)
  
!<description>
  ! Calls the initData subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with sptils_initData.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure 
  type(t_sptilsNode), intent(INOUT)         :: rsolverNode
!</inputoutput>
  
!<output>
  ! One of the SPTILS_ERR_XXXX constants. A value different to 
  ! SPTILS_ERR_NOERROR indicates that an error happened during the
  ! initialisation phase.
  integer, intent(OUT) :: ierror
!</output>

!</subroutine>

    ! local variables
    integer :: ilev
    type(t_sptilsMGLevelInfo), pointer :: p_rmgLevel
    
    ! A-priori we have no error...
    ierror = SPTILS_ERR_NOERROR
    
    ! On each level, call the initStructure routine of the
    ! presmoother/postsmoother/coarse grid solver

    ! Loop through the level. 
    do ilev=1,size(rsolverNode%p_rsubnodeMultigrid%p_Rlevels)
     
      p_rmgLevel => rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)
    
      ! Call the setmatrices-routine for the presmoother/postsmoother/
      ! coarse grid solver
      if (.not. associated(p_rmgLevel%p_rpresmoother,p_rmgLevel%p_rpostsmoother)) then
        ! May be the presmoother and postsmoother are identical; initialise them
        ! only once!
        if (associated(p_rmgLevel%p_rpresmoother)) then
          call sptils_initData(p_rmgLevel%p_rpresmoother,ierror)
          if (ierror .ne. ierror) return
        end if
      end if
      
      if (associated(p_rmgLevel%p_rpostsmoother)) then
        call sptils_initData(p_rmgLevel%p_rpostsmoother,ierror)
        if (ierror .ne. ierror) return
      end if
      
      if (associated(p_rmgLevel%p_rcoarseGridSolver)) then
        call sptils_initData(p_rmgLevel%p_rcoarseGridSolver,ierror)
        if (ierror .ne. ierror) return
      end if
    end do

  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine sptils_doneDataMultigrid (rsolverNode)
  
!<description>
  ! Calls the doneData subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with sptils_doneData.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure 
  type(t_sptilsNode), intent(INOUT)         :: rsolverNode
!</inputoutput>
  
!</subroutine>

    ! local variables
    integer :: ilev
    type(t_sptilsMGLevelInfo), pointer :: p_rmgLevel
    
    ! On each level, call the initStructure routine of the
    ! presmoother/postsmoother/coarse grid solver

    ! Loop through the level. 
    do ilev=1,size(rsolverNode%p_rsubnodeMultigrid%p_Rlevels)
    
      p_rmgLevel => rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)
    
      ! Call the setmatrices-routine for the presmoother/postsmoother/
      ! coarse grid solver
      if (.not. associated(p_rmgLevel%p_rpresmoother,p_rmgLevel%p_rpostsmoother)) then
        ! May be the presmoother and postsmoother are identical; release them
        ! only once!
        if (associated(p_rmgLevel%p_rpresmoother)) then
          call sptils_doneData(p_rmgLevel%p_rpresmoother)
        end if
      end if
      
      if (associated(p_rmgLevel%p_rpostsmoother)) then
        call sptils_doneData(p_rmgLevel%p_rpostsmoother)
      end if
      
      if (associated(p_rmgLevel%p_rcoarseGridSolver)) then
        call sptils_doneData(p_rmgLevel%p_rcoarseGridSolver)
      end if
    end do

  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine sptils_doneStructureMultigrid (rsolverNode)
  
!<description>
  ! Calls the doneStructure subroutine of the subsolver.
  ! Maybe the subsolver needs that...
  ! The routine is declared RECURSIVE to get a clean interaction
  ! with sptils_doneStructure.
!</description>
  
!<inputoutput>
  ! The t_sptilsNode structure 
  type(t_sptilsNode), intent(INOUT)         :: rsolverNode
!</inputoutput>
  
!</subroutine>

    ! local variables
    integer :: ilev,NLMAX
    type(t_sptilsMGLevelInfo), pointer :: p_rmgLevel
    
    ! On each level, call the initStructure routine of the
    ! presmoother/postsmoother/coarse grid solver

    ! Loop through the level. 
    NLMAX = size(rsolverNode%p_rsubnodeMultigrid%p_Rlevels)
    do ilev=1,NLMAX
    
      p_rmgLevel => rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilev)
    
      ! Call the setmatrices-routine for the presmoother/postsmoother/
      ! coarse grid solver
      if (.not. associated(p_rmgLevel%p_rpresmoother,p_rmgLevel%p_rpostsmoother)) then
        ! May be the presmoother and postsmoother are identical; release them
        ! only once!
        if (associated(p_rmgLevel%p_rpresmoother)) then
          call sptils_doneStructure(p_rmgLevel%p_rpresmoother)
        end if
      end if
      
      if (associated(p_rmgLevel%p_rpostsmoother)) then
        call sptils_doneStructure(p_rmgLevel%p_rpostsmoother)
      end if
      
      if (associated(p_rmgLevel%p_rcoarseGridSolver)) then
        call sptils_doneStructure(p_rmgLevel%p_rcoarseGridSolver)
      end if


      ! Release the projection structure
      call sptipr_doneProjection (p_rmgLevel%rinterlevelProjection)
              
      ! Release vectors
              
      if (ilev .lt. NLMAX) then
        call sptivec_releaseVector (p_rmgLevel%rsolutionVector)
      end if
      
      if (ilev .gt. 1) then
        call sptivec_releaseVector (p_rmgLevel%rrhsVector)
        call sptivec_releaseVector (p_rmgLevel%rtempVector)
      end if
      
      ! If adaptive coarse grid correction is activated, we need a temporary
      ! vector for the adaptive coarse grid correction on every level.
      if (rsolverNode%p_rsubnodeMultigrid%dalphamin .ne. &
          rsolverNode%p_rsubnodeMultigrid%dalphamax) then
        call sptivec_releaseVector (p_rmgLevel%rtempCGCvector)
      end if

      ! Release the temp vector for prolongation/restriction
      call lsysbl_releaseVector (p_rmgLevel%rprjVector)

    end do

  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  subroutine sptils_setMultigridLevel (rsolverNode,ilevel,&
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
  type(t_sptilsNode)                             :: rsolverNode
!</inputoutput>
  
!<input>
  ! Level in the MG solver to be initialised.
  integer, intent(IN) :: ilevel

  ! An interlevel projection structure that configures the projection
  ! between the solutions on a finer and a coarser grid. The structure
  ! must have been initialised with mlprj_initProjection.
  !
  ! Note that this structure is level-independent (as long as the
  ! spatial discretisation structures on different levels are 'compatible'
  ! what they have to be anyway), so the same structure can be used
  ! to initialise all levels!
  type(t_sptiProjection) :: rinterlevelProjection
  
  ! Optional: A pointer to the solver structure of a solver that should be 
  ! used for presmoothing. This structure is used as a template to create an
  ! appropriate solver structure for all the levels. The structure itself is
  ! used at the finest level.
  ! If not given or set to NULL(), no presmoother will be used.
  type(t_sptilsNode), pointer, optional   :: p_rpresmoother
  
  ! Optional: A pointer to the solver structure of a solver that should be 
  ! used for postsmoothing. This structure is used as a template to create an
  ! appropriate solver structure for all the levels. The structure itself is
  ! used at the finest level.
  ! If not given or set to NULL(), no presmoother will be used.
  type(t_sptilsNode), pointer, optional   :: p_rpostsmoother

  ! Optional: A pointer to the solver structure of a solver that should be 
  ! used for coarse grid solving. 
  ! Should only be given for the very first level.
  type(t_sptilsNode), pointer, optional   :: p_rcoarseGridSolver
!</input>
  
!</subroutine>

    ! local variables
    type(t_sptilsMGLevelInfo), pointer :: p_rlevelInfo
    
    ! Make sure the solver node is configured for multigrid
    if ((rsolverNode%calgorithm .ne. SPTILS_ALG_MULTIGRID) .or. &
        (.not. associated(rsolverNode%p_rsubnodeMultigrid))) then
      print *,'Error: Multigrid structure not initialised'
      call sys_halt()
    end if
    
    p_rlevelInfo => rsolverNode%p_rsubnodeMultigrid%p_Rlevels(ilevel)
    
    ! Attach the sub-solvers
    if (present(p_rcoarseGridSolver)) then
      p_rlevelInfo%p_rcoarseGridSolver => p_rcoarseGridSolver
    else
      nullify(p_rlevelInfo%p_rcoarseGridSolver)
    end if

    if (present(p_rpresmoother)) then
      p_rlevelInfo%p_rpresmoother => p_rpresmoother
    else
      nullify(p_rlevelInfo%p_rpresmoother)
    end if

    if (present(p_rpostsmoother)) then
      p_rlevelInfo%p_rpostsmoother => p_rpostsmoother
    else
      nullify(p_rlevelInfo%p_rpostsmoother)
    end if
    
    ! Initialise the interlevel projection structure,
    ! copy the data of our given template.
    p_rlevelInfo%rinterlevelProjection = rinterlevelProjection

  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>
  
  subroutine sptils_convertToSmoother (rsolverNode,nsmoothingSteps,domega)
  
!<description>
  ! Converts a t_linsolNode to a smoother structure. A smoother is a solver
  ! that performs a fixed number of iterations without respecting any
  ! residuum. nsmoothingSteps is the number of steps the smoother should
  ! perform.
!</description>
  
!<input>
  ! Number of steps the smoother should perform
  integer, intent(IN)          :: nsmoothingSteps
  
  ! OPTIONAL: Damping parameter.
  ! The parameter rsolverNode%domega is set to this value in order
  ! to set up the damping.
  real(DP), intent(IN), optional :: domega
!</input>
  
!<inputoutput>
  ! The t_sptilsNode structure of the solver which should be configured as smoother.
  type(t_sptilsNode), intent(INOUT), target :: rsolverNode
!</inputoutput>
  
!</subroutine>

    rsolverNode%depsRel = 0.0_DP
    rsolverNode%depsAbs = 0.0_DP
    rsolverNode%nminIterations = nsmoothingSteps
    rsolverNode%nmaxIterations = nsmoothingSteps
    
    if (present(domega)) rsolverNode%domega = domega

  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine sptils_smoothCorrection (rsolverNode,rx,rb,rtemp,rspatialTemp)
  
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
  type(t_spacetimeVector), intent(IN), target        :: rb
!</input>
  
!<inputoutput>
  ! The t_sptilsNode structure of the solver
  type(t_sptilsNode), intent(INOUT), target :: rsolverNode
  
  ! The initial solution vector; receives the solution of the system
  type(t_spacetimeVector), intent(INOUT)                 :: rx
  
  ! A temporary vector of the same size and structure as rx.
  type(t_spacetimeVector), intent(INOUT)                 :: rtemp

  ! A spatial temporary vector of the same size and structure as every timestep in rx.
  type(t_vectorBlock), intent(INOUT)                     :: rspatialTemp
!</inputoutput>
  
!</subroutine>

    integer :: i
    integer :: iiterations
    real(DP) :: dres,dresInit
    type(t_ccoptSpaceTimeMatrix), pointer :: p_rmatrix
    type(t_ccoptSpaceTimeDiscretisation), pointer :: p_rspaceTimeDiscr
    type(t_timer) :: rtimer
    !DEBUG: REAL(DP), DIMENSION(:), POINTER :: p_Ddata,p_Ddata2
    
    call stat_clearTimer (rtimer)
    
    ! Cancel if nmaxIterations = number of smoothing steps is =0.
    if (rsolverNode%nmaxIterations .le. 0) return
    
    ! This is a 1-step solver, we have to emulate the smoothing
    ! iterations. Perform rsolverNode%nmaxIterations steps of the
    ! form
    !     $$ x_{n+1} = x_n + P^{-1}(b-Ax_n) $$
    ! with $x_0 = 0$.
    
    !DEBUG: CALL lsysbl_getbase_double (rx,p_Ddata)
    !DEBUG: CALL lsysbl_getbase_double (rtemp,p_Ddata2)
    
    ! Apply nmaxIterations times defect correction to the given solution rx.
    p_rmatrix => rsolverNode%rmatrix
    p_rspaceTimeDiscr => rsolverNode%rmatrix%p_rspaceTimeDiscr
    
    ! Do we have an iterative or one-step solver given?
    ! A 1-step solver performs the following loop nmaxIterations times, while an iterative
    ! solver is called only once and performs nmaxIterations steps internally.
    ! (This is a convention. Calling an iterative solver i times with j internal steps
    ! would also be possible, but we don't implement that here.)
    if (iand(rsolverNode%ccapability,SPTILS_ABIL_DIRECT) .ne. 0) then
      iiterations = rsolverNode%nmaxIterations
    else
      iiterations = 1
    end if
    
    ! DEBUG!!!
    !call sptivec_saveToFileSequence (rx,'(''smoothin.txt.'',I5.5)',.true.)
    
    do i=1,iiterations
    
      call sptivec_copyVector(rb,rtemp)
      call cc_spaceTimeMatVec (rsolverNode%p_rproblem, p_rmatrix, &
          rx,rtemp, -1.0_DP,1.0_DP,SPTID_FILTER_DEFECT,dres)
      
      if (iiterations .eq. 1) dresInit = dres
      if (dresInit .eq. 0) dresInit = 1.0_DP
      
      ! Implement boundary conditions into the defect
      call tbc_implementInitCondDefect (&
          p_rspaceTimeDiscr,rtemp,rspatialTemp)
      call tbc_implementBCdefect (rsolverNode%p_rproblem,&
          p_rspaceTimeDiscr,rtemp,rspatialTemp)
      
      if (rsolverNode%ioutputLevel .eq. 1) then
        call output_line ('Space-Time-Smoother: Step '//trim(sys_siL(i-1,10)))
      end if

      if (rsolverNode%ioutputLevel .ge. 2) then
        if (.not.((dres .ge. 1E-99_DP) .and. (dres .le. 1E99_DP))) dres = 0.0_DP
        call output_line ('Space-Time-Smoother: Step '//trim(sys_siL(i-1,10))//&
            ' !!RES!! = '//trim(sys_sdEL(dres,15)) )
      end if
      
      ! Check for convergence
      if (rsolverNode%istoppingCriterion .eq. 0) then
        if ((dres .lt. rsolverNode%depsAbs) .and. &
            (dres .lt. rsolverNode%depsRel*dresInit)) exit
      else
        if ((dres .lt. rsolverNode%depsAbs) .or. &
            (dres .lt. rsolverNode%depsRel*dresInit)) exit
      end if
      
      ! Stop the time for space-preconditioning and sum it up. Return the sum as result.
      call sptils_precondDefect(rsolverNode,rtemp)
      call stat_addTimers (rsolverNode%rtimeSpacePrecond,rtimer)
      call sptivec_vectorLinearComb (rtemp,rx,1.0_DP,1.0_DP)
      
    end do
    
    ! DEBUG!!!
    !call sptivec_saveToFileSequence (rx,'(''smoothout.txt.'',I5.5)',.true.)
    
    rsolverNode%rtimeSpacePrecond = stat_rcloneTimer (rtimer)

    ! Probably print the final residuum
    if (rsolverNode%ioutputLevel .ge. 2) then
      if (i .ge. iiterations) then
        ! We only have to recalculate the residuum if we haven't stopped
        ! prematurely.
        call sptivec_copyVector(rb,rtemp)
        call cc_spaceTimeMatVec (rsolverNode%p_rproblem, p_rmatrix, &
            rx,rtemp, -1.0_DP,1.0_DP,SPTID_FILTER_DEFECT,dres)
        if (.not.((dres .ge. 1E-99_DP) .and. (dres .le. 1E99_DP))) dres = 0.0_DP
                  
        if (iand(rsolverNode%ccapability,SPTILS_ABIL_DIRECT) .ne. 0) then
          call output_line ('Space-Time-Smoother: Step '//trim(sys_siL(i-1,10))//&
            ' !!RES!! = '//trim(sys_sdEL(dres,15)) )
        else
          call output_line ('Space-Time-Smoother: Step '//&
            trim(sys_siL(rsolverNode%nmaxIterations,10))//&
            ' !!RES!! = '//trim(sys_sdEL(dres,15)) )
        end if
      end if
    end if
    
  end subroutine

  ! ***************************************************************************
  
!<subroutine>
  
  recursive subroutine sptils_precMultigrid (rsolverNode,rd)
  
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
  type(t_sptilsNode), intent(INOUT), target :: rsolverNode

  ! On call to this routine: The defect vector to be preconditioned.
  ! Will be overwritten by the preconditioned defect.
  type(t_spacetimeVector), intent(INOUT)    :: rd
!</inputoutput>
  
!</subroutine>

    !CALL sptils_precondDefect (rsolverNode%p_rsubnodeMultigrid%p_Rlevels(1)%p_rcoarseGridSolver,&
    !  rd)

  ! local variables
  integer :: nminIterations,nmaxIterations,niteResOutput
  integer :: ite,ilev,nlmax,i,nlmin
  real(DP) :: dres,dstep
  
  ! The system matrix on the current level
  type(t_ccoptSpaceTimeMatrix), pointer :: p_rmatrix
  type(t_ccoptSpaceTimeDiscretisation), pointer :: p_rspaceTimeDiscr
  
  ! Our MG structure
  type(t_sptilsSubnodeMultigrid), pointer :: p_rsubnode
  
  ! For statistics:
  type(t_timer) :: rtimer
  
    ! Solve the system!
    
    ! Getch some information
    p_rsubnode => rsolverNode%p_rsubnodeMultigrid
    
    ! Initialise timers
    call stat_clearTimer (p_rsubnode%rtimeSmoothing)
    call stat_clearTimer (p_rsubnode%rtimeCoarseGridSolver)
    call stat_clearTimer (p_rsubnode%rtimeLinearAlgebra)
    call stat_clearTimer (p_rsubnode%rtimeProlRest)
    
    call stat_clearTimer (rsolverNode%rtimeSpacePrecond)

    ! Get the system matrix on the finest level:
    p_rmatrix => rsolverNode%rmatrix
    p_rspaceTimeDiscr => rsolverNode%rmatrix%p_rspaceTimeDiscr

    if (p_rsubnode%icycle .lt. 0) then
      ! Wrong cycle
      print *,'Multigrid: Wrong cycle!'
      rsolverNode%iresult = 2
      return
    end if
    
    ! Minimum/maximum number of iterations
    nminIterations = max(rsolverNode%nminIterations,0)
    nmaxIterations = max(rsolverNode%nmaxIterations,0)
      
    ! Iteration when the residuum is printed:
    niteResOutput = max(1,rsolverNode%niteResOutput)

    ! Status reset
    rsolverNode%iresult = 0
    rsolverNode%icurrentIteration = 0
    rsolverNode%dinitialDefect = 0.0_DP
    rsolverNode%dlastDefect = 0.0_DP
    rsolverNode%dfinalDefect = 0.0_DP
    rsolverNode%dconvergenceRate = 0.0_DP
    
    ! We start on the maximum level. 
    ilev = p_rsubnode%NLMAX
    nlmax = p_rsubnode%NLMAX
    nlmin = p_rsubnode%NLMIN
    
    ! Is there only one level? Can be seen if the current level
    ! already contains a coarse grid solver.
    if (nlmin .eq. nlmax) then
    
      if (rsolverNode%ioutputLevel .gt. 1) then
        call output_line ('Space-Time-Multigrid: Only one level. '//&
             'Switching back to standard solver.')
      end if

      ! DEBUG!!!
      !CALL sptivec_saveToFileSequence (rd,&
      !    '(''./debugdata/vector2.txt.'',I5.5)',.TRUE.)

      call stat_startTimer (p_rsubnode%rtimeCoarseGridSolver)
      call sptils_precondDefect(p_rsubnode%p_Rlevels(ilev)%p_rcoarseGridSolver,rd)
      call stat_stopTimer (p_rsubnode%rtimeCoarseGridSolver)
      
      ! Sum up the time for space preconditioning
      call stat_addTimers (&
          p_rsubnode%p_Rlevels(ilev)%p_rcoarseGridSolver%rtimeSpacePrecond,&
          rsolverNode%rtimeSpacePrecond)
      
      ! Take the statistics from the coarse grid solver.
      rsolverNode%dinitialDefect = &
        p_rsubnode%p_Rlevels(ilev)%p_rcoarseGridSolver%dinitialDefect
      rsolverNode%dlastDefect = &
        p_rsubnode%p_Rlevels(ilev)%p_rcoarseGridSolver%dlastDefect
      rsolverNode%dfinalDefect = &
        p_rsubnode%p_Rlevels(ilev)%p_rcoarseGridSolver%dfinalDefect
      rsolverNode%dconvergenceRate = &
        p_rsubnode%p_Rlevels(ilev)%p_rcoarseGridSolver%dconvergenceRate
      rsolverNode%iiterations = &
        p_rsubnode%p_Rlevels(ilev)%p_rcoarseGridSolver%iiterations
      
    else
      ! rtimeLinearAlgebra calculates the total time. By subtracting all
      ! other time information, we'll get at the end the time for
      ! linear algebra.
      call stat_startTimer (p_rsubnode%rtimeLinearAlgebra)
    
      ! Get the norm of the initial residuum.
      ! As the initial iteration vector is zero, this is the norm
      ! of the RHS:
      dres = sptivec_vectorNorm (rd,rsolverNode%iresNorm)
      if (.not.((dres .ge. 1E-99_DP) .and. &
                (dres .le. 1E99_DP))) dres = 0.0_DP

      ! Initialize starting residuum
      rsolverNode%dinitialDefect = dres
      rsolverNode%dlastDefect = 0.0_DP

      ! Check if out initial defect is zero. This may happen if the filtering
      ! routine filters "everything out"!
      ! In that case we can directly stop our computation.

      if ( rsolverNode%dinitialDefect .lt. rsolverNode%drhsZero ) then
      
        ! final defect is 0, as initialised in the output variable above
        call sptivec_clearVector(rd)
        rsolverNode%dlastDefect = dres
        rsolverNode%dfinalDefect = dres
        rsolverNode%dconvergenceRate = 0.0_DP
        rsolverNode%iiterations = 0
        
      else
        
        ! At first, reset the cycle counters of all levels
        do i=nlmin,nlmax
          if (p_rsubnode%icycle .eq. 0) then
            p_rsubnode%p_Rlevels(i)%ncycles = 2
          else
            p_rsubnode%p_Rlevels(i)%ncycles = p_rsubnode%icycle
          end if  
        end do
        
        ! Afterwards, p_rcurrentLevel points to the maximum level again.
        ! There, we set ncycles to 1.
        p_rsubnode%p_Rlevels(ilev)%ncycles = 1
        
        ! Print out the initial residuum

        if (rsolverNode%ioutputLevel .ge. 2) then
          call output_line ('Space-Time-Multigrid: Iteration '// &
              trim(sys_siL(0,10))//',  !!RES!! = '//&
              trim(sys_sdEL(rsolverNode%dinitialDefect,15)) )
        end if
        
        ! Copy the initial RHS to the RHS vector on the maximum level.
        call sptivec_copyVector(rd,p_rsubnode%p_Rlevels(ilev)%rrhsVector)
        
        ! Replace the solution vector on the finest level by rd.
        ! Afterwards, rd and the solution vector on the finest level
        ! share the same memory location, so we use rd as iteration
        ! vector directly.
        p_rsubnode%p_Rlevels(ilev)%rsolutionVector = rd
        
        ! Clear the initial solution vector.
        call sptivec_clearVector (p_rsubnode%p_Rlevels(ilev)%rsolutionVector)
        
        ! Start multigrid iteration; perform at most nmaxiterations iterations.
        do ite = 1, nmaxiterations
        
          rsolverNode%icurrentIteration = ite
          
          ! Initialize cycle counters for all levels.
          do i=nlmin,nlmax
            p_rsubnode%p_Rlevels(i)%ncyclesRemaining = &
                p_rsubnode%p_Rlevels(i)%ncycles
          end do
        
          ! p_rcurrentLevel now points to the maximum level.
          ilev = nlmax

          ! Get the system matrix on the finest level:
          p_rmatrix => p_rsubnode%p_Rlevels(ilev)%rmatrix
          
          if (rsolverNode%ioutputLevel .ge. 3) then
            call output_line (&
              'Space-Time-Multigrid: Current mesh level: '//trim(sys_siL(ilev,5)))
          end if
          
          ! Build the defect...
          call sptivec_copyVector (&
              p_rsubnode%p_Rlevels(ilev)%rrhsVector,&
              p_rsubnode%p_Rlevels(ilev)%rtempVector)
          if (ite .ne. 1) then   ! initial solution vector is zero!
            call cc_spaceTimeMatVec (rsolverNode%p_rproblem, p_rmatrix, &
                p_rsubnode%p_Rlevels(ilev)%rsolutionVector,&
                p_rsubnode%p_Rlevels(ilev)%rtempVector, -1.0_DP,1.0_DP,&
                SPTID_FILTER_DEFECT)
          
            ! Implement boundary conditions into the defect
            call tbc_implementInitCondDefect (&
                p_rsubnode%p_Rlevels(ilev)%rmatrix%p_rspaceTimeDiscr,&
                p_rsubnode%p_Rlevels(ilev)%rtempVector,&
                p_rsubnode%p_Rlevels(ilev)%rprjVector)
            call tbc_implementBCdefect (rsolverNode%p_rproblem,&
                p_rsubnode%p_Rlevels(ilev)%rmatrix%p_rspaceTimeDiscr,&
                p_rsubnode%p_Rlevels(ilev)%rtempVector,&
                p_rsubnode%p_Rlevels(ilev)%rprjVector)
          end if
          
          cycleloop: do  ! Loop for the cycles
          
            ! On the maximum level we already built out defect vector. If we are    
            ! on a lower level than NLMAX, perform smoothing+restriction down to the
            ! coarse level. We identify the coarse level by checking if
            ! the current level has a coarse grid solver.
            
            do while (ilev .gt. nlmin)
            
              ! Perform the pre-smoothing with the current solution vector
              if (associated(p_rsubnode%p_Rlevels(ilev)%p_rpreSmoother)) then
                call stat_startTimer (p_rsubnode%rtimeSmoothing)
                call sptils_smoothCorrection (&
                          p_rsubnode%p_Rlevels(ilev)%p_rpreSmoother,&
                          p_rsubnode%p_Rlevels(ilev)%rsolutionVector,&
                          p_rsubnode%p_Rlevels(ilev)%rrhsVector,&
                          p_rsubnode%p_Rlevels(ilev)%rtempVector,&
                          p_rsubnode%p_Rlevels(ilev)%rprjVector)
                call stat_stopTimer (p_rsubnode%rtimeSmoothing)

                call stat_addTimers (&
                    p_rsubnode%p_Rlevels(ilev)%p_rpreSmoother%rtimeSpacePrecond,&
                    rsolverNode%rtimeSpacePrecond)

              end if
            
              ! Build the defect vector
              call sptivec_copyVector (p_rsubnode%p_Rlevels(ilev)%rrhsVector,&
                  p_rsubnode%p_Rlevels(ilev)%rtempVector)
              call cc_spaceTimeMatVec (rsolverNode%p_rproblem, p_rmatrix, &
                  p_rsubnode%p_Rlevels(ilev)%rsolutionVector,&
                  p_rsubnode%p_Rlevels(ilev)%rtempVector, -1.0_DP,1.0_DP,&
                  SPTID_FILTER_DEFECT,dres)
              
              ! Extended output
              if (associated(p_rsubnode%p_Rlevels(ilev)%p_rpreSmoother) .and. &
                  (rsolverNode%ioutputLevel .ge. 3) .and. &
                  (mod(ite,niteResOutput) .eq. 0)) then
                  
                if (.not.((dres .ge. 1E-99_DP) .and. &
                          (dres .le. 1E99_DP))) dres = 0.0_DP
                          
                call output_line ('Space-Time-Multigrid: Level '//trim(sys_siL(ilev,5))//&
                    ' after presm.:     !!RES!! = '//trim(sys_sdEL(dres,15)) )
              end if

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
              if (ilev .gt. nlmin+1) then
                call stat_startTimer (p_rsubnode%rtimeProlRest)
              
                ! We don't project to the coarse grid
                call sptipr_performRestriction (&
                      p_rsubnode%p_Rlevels(ilev)%rinterlevelProjection,&
                      p_rsubnode%p_Rlevels(ilev-1)%rrhsVector, &
                      p_rsubnode%p_Rlevels(ilev)%rtempVector, &
                      p_rsubnode%p_Rlevels(ilev-1)%rprjVector, &
                      p_rsubnode%p_Rlevels(ilev)%rprjVector)

                ! Implement boundary conditions into the defect
                call tbc_implementInitCondDefect (&
                    p_rsubnode%p_Rlevels(ilev-1)%rmatrix%p_rspaceTimeDiscr,&
                    p_rsubnode%p_Rlevels(ilev-1)%rrhsVector,&
                    p_rsubnode%p_Rlevels(ilev-1)%rprjVector)
                call tbc_implementBCdefect (rsolverNode%p_rproblem,&
                    p_rsubnode%p_Rlevels(ilev-1)%rmatrix%p_rspaceTimeDiscr,&
                    p_rsubnode%p_Rlevels(ilev-1)%rrhsVector,&
                    p_rsubnode%p_Rlevels(ilev-1)%rprjVector)
                    
                call stat_stopTimer (p_rsubnode%rtimeProlRest)

                ! Choose zero as initial vector on lower level. 
                call sptivec_clearVector (p_rsubnode%p_Rlevels(ilev-1)%rsolutionVector)
                
                ! Extended output and/or adaptive cycles
                if (rsolverNode%ioutputLevel .ge. 3) then
                
                  dres = sptivec_vectorNorm (p_rsubnode%p_Rlevels(ilev-1)%rrhsVector,&
                      rsolverNode%iresNorm)
                  if (.not.((dres .ge. 1E-99_DP) .and. &
                            (dres .le. 1E99_DP))) dres = 0.0_DP
                            
                  ! In case adaptive cycles are activated, save the 'initial' residual
                  ! of that level into the level structure. Then we can later check
                  ! if we have to repeat the cycle on the coarse mesh.
                  if ((p_rsubnode%p_Rlevels(ilev-1)%depsRelCycle .ne. 1E99_DP) .or.&
                      (p_rsubnode%p_Rlevels(ilev-1)%depsAbsCycle .ne. 1E99_DP)) then
                    p_rsubnode%p_Rlevels(ilev-1)%dinitResCycle = dres
                    p_rsubnode%p_Rlevels(ilev-1)%icycleCount = 1
                  end if

                  ! If the output level is high enough, print that residuum norm.   
                  if (mod(ite,niteResOutput).eq. 0) then
                    call output_line ('Space-Time-Multigrid: Level '//trim(sys_siL(ilev-1,5))//&
                        ' after restrict.:  !!RES!! = '//trim(sys_sdEL(dres,15)) )
                  end if
                  
                end if

              else
              
                call stat_startTimer (p_rsubnode%rtimeProlRest)
              
                ! The vector is to be restricted to the coarse grid.
                call sptipr_performRestriction (&
                      p_rsubnode%p_Rlevels(ilev)%rinterlevelProjection,&
                      p_rsubnode%p_Rlevels(ilev-1)%rsolutionVector, &
                      p_rsubnode%p_Rlevels(ilev)%rtempVector, &
                      p_rsubnode%p_Rlevels(ilev-1)%rprjVector, &
                      p_rsubnode%p_Rlevels(ilev)%rprjVector)

                ! DEBUG!!!
                !CALL sptivec_saveToFileSequence (p_rsubnode%p_Rlevels(ilev-1)%rsolutionVector,&
                !    '(''./debugdata/vector1.txt.'',I5.5)',.TRUE.)

                ! Implement boundary conditions into the defect
                call tbc_implementInitCondDefect (&
                    p_rsubnode%p_Rlevels(ilev-1)%rmatrix%p_rspaceTimeDiscr,&
                    p_rsubnode%p_Rlevels(ilev-1)%rsolutionVector,&
                    p_rsubnode%p_Rlevels(ilev-1)%rprjVector)
                call tbc_implementBCdefect (rsolverNode%p_rproblem,&
                    p_rsubnode%p_Rlevels(ilev-1)%rmatrix%p_rspaceTimeDiscr,&
                    p_rsubnode%p_Rlevels(ilev-1)%rsolutionVector,&
                    p_rsubnode%p_Rlevels(ilev-1)%rprjVector)

                call stat_stopTimer (p_rsubnode%rtimeProlRest)

                ! Extended output
                if ((rsolverNode%ioutputLevel .ge. 3) .and. &
                    (mod(ite,niteResOutput).eq.0)) then
                    
                  dres = sptivec_vectorNorm (p_rsubnode%p_Rlevels(ilev-1)%rsolutionVector,&
                      rsolverNode%iresNorm)
                  if (.not.((dres .ge. 1E-99_DP) .and. &
                            (dres .le. 1E99_DP))) dres = 0.0_DP
                            
                  call output_line ('Space-Time-Multigrid: Level '//trim(sys_siL(ilev-1,5))//&
                      ' after restrict.:  !!RES!! = '//trim(sys_sdEL(dres,15)) )
                end if

              end if              
            
              ! Go down one level
              ilev = ilev - 1
              p_rmatrix => p_rsubnode%p_Rlevels(ilev)%rmatrix

              if (rsolverNode%ioutputLevel .ge. 3) then
                call output_line ('Space-Time-Multigrid: Current mesh level: '//trim(sys_siL(ilev,5)))
              end if
              
              ! If we are not on the lowest level, repeat the smoothing of 
              ! the solution/restriction of the new defect in the next loop 
              ! pass...
            end do   ! ilev > minimum level
            
            ! Now we reached the coarse grid.
            
            if (rsolverNode%ioutputLevel .ge. 3) then
              call output_line ('Space-Time-Multigrid: Invoking coarse grid solver.')
            end if
            
            call stat_startTimer (p_rsubnode%rtimeCoarseGridSolver)
            
            ! Solve the system on lowest level by preconditioning
            ! of the RHS=defect vector.
            call sptils_precondDefect(p_rsubnode%p_Rlevels(ilev)%p_rcoarseGridSolver,&
                                      p_rsubnode%p_Rlevels(ilev)%rsolutionVector)
                                      
            call stat_stopTimer (p_rsubnode%rtimeCoarseGridSolver)

            ! Sum up the time for space preconditioning
            call stat_addTimers (&
                p_rsubnode%p_Rlevels(ilev)%p_rcoarseGridSolver%rtimeSpacePrecond,&
                rsolverNode%rtimeSpacePrecond)
            
            ! Now we have the solution vector on the lowest level - we have to go
            ! upwards now... but probably not to NLMAX! That depends on the cycle.
            !
            do while (ilev .lt. nlmax)
              ilev = ilev + 1
              p_rmatrix => p_rsubnode%p_Rlevels(ilev)%rmatrix
              
              if (rsolverNode%ioutputLevel .ge. 3) then
                call output_line ('Space-Time-Multigrid: Current mesh level: '&
                    //trim(sys_siL(ilev,5)))
              end if

              call stat_startTimer (p_rsubnode%rtimeProlRest)
              
              ! Prolongate the solution vector from the coarser level
              ! to the temp vector on the finer level.
              call sptipr_performProlongation (&
                    p_rsubnode%p_Rlevels(ilev)%rinterlevelProjection,&
                    p_rsubnode%p_Rlevels(ilev-1)%rsolutionVector, &
                    p_rsubnode%p_Rlevels(ilev)%rtempVector, &
                    p_rsubnode%p_Rlevels(ilev-1)%rprjVector, &
                    p_rsubnode%p_Rlevels(ilev)%rprjVector,&
                    p_rmatrix%p_rspaceTimeDiscr,&
                    p_rsubnode%p_Rlevels(ilev-1)%rmatrix%p_rspaceTimeDiscr,&
                    rsolverNode%p_rproblem)
                    
              ! Implement boundary conditions into the vector.
              ! It's still a defect, although a preconditioned one.
              call tbc_implementInitCondDefect (&
                  p_rsubnode%p_Rlevels(ilev)%rmatrix%p_rspaceTimeDiscr,&
                  p_rsubnode%p_Rlevels(ilev)%rtempVector, &
                  p_rsubnode%p_Rlevels(ilev)%rprjVector)
              call tbc_implementBCdefect (rsolverNode%p_rproblem,&
                  p_rsubnode%p_Rlevels(ilev)%rmatrix%p_rspaceTimeDiscr,&
                  p_rsubnode%p_Rlevels(ilev)%rtempVector,&
                  p_rsubnode%p_Rlevels(ilev)%rprjVector)

              call stat_stopTimer (p_rsubnode%rtimeProlRest)

              ! Step length control. By default, choose step length dalphamin
              ! which is usually = 1.0.
              dstep = p_rsubnode%dalphamin ! 1.0_DP
              
              ! If adaptive coarse grid correction is activated, get the
              ! step length parameter by energy minimisation.
              if (p_rsubnode%dalphamin .ne. p_rsubnode%dalphamax) then
              
                ! Calculate the optimal alpha.
                if (.not. associated(p_rsubnode%p_DequationWeights)) then
                  call calcCGCenergyMin (rsolverNode%p_rproblem,&
                      p_rsubnode%p_Rlevels(ilev)%rsolutionVector,&
                      p_rsubnode%p_Rlevels(ilev)%rtempVector,&
                      p_rsubnode%p_Rlevels(ilev)%rrhsVector,&
                      p_rmatrix,&
                      p_rsubnode%p_Rlevels(ilev)%rtempCGCvector,&
                      p_rsubnode%p_Rlevels(ilev)%rprjVector,&
                      p_rsubnode%dalphamin, p_rsubnode%dalphamax, dstep)
                else
                  call calcCGCenergyMinWeighted (rsolverNode%p_rproblem,&
                      p_rsubnode%p_Rlevels(ilev)%rsolutionVector,&
                      p_rsubnode%p_Rlevels(ilev)%rtempVector,&
                      p_rsubnode%p_Rlevels(ilev)%rrhsVector,&
                      p_rmatrix,&
                      p_rsubnode%p_Rlevels(ilev)%rtempCGCvector,&
                      p_rsubnode%p_Rlevels(ilev)%rprjVector,&
                      p_rsubnode%dalphamin, p_rsubnode%dalphamax, &
                      p_rsubnode%p_DequationWeights,dstep)
                end if
                
                ! Some output
                if (rsolverNode%ioutputLevel .ge. 2) then
                  call output_line ('Coarse grid correction factor: '//&
                    'dstep = '//trim(sys_sdEL(dstep,15)) )
                end if
                
              end if
              
              ! Perform the coarse grid correction by adding the coarse grid
              ! solution (with the calculated step-length parameter) to
              ! the current solution.
              
              call sptivec_vectorLinearComb (p_rsubnode%p_Rlevels(ilev)%rtempVector,&
                 p_rsubnode%p_Rlevels(ilev)%rsolutionVector,&
                 dstep,1.0_DP)
                            
              ! Extended output
              if ((rsolverNode%ioutputLevel .ge. 3) .and. &
                  (mod(ite,niteResOutput).eq.0)) then
                  
                call sptivec_copyVector (p_rsubnode%p_Rlevels(ilev)%rrhsVector,&
                                         p_rsubnode%p_Rlevels(ilev)%rtempVector)
                call cc_spaceTimeMatVec (rsolverNode%p_rproblem, p_rmatrix, &
                    p_rsubnode%p_Rlevels(ilev)%rsolutionVector,&
                    p_rsubnode%p_Rlevels(ilev)%rtempVector, -1.0_DP,1.0_DP,&
                    SPTID_FILTER_DEFECT,dres)

                if (.not.((dres .ge. 1E-99_DP) .and. &
                          (dres .le. 1E99_DP))) dres = 0.0_DP
                          
                call output_line ('Space-Time-Multigrid: Level '//trim(sys_siL(ilev,5))//&
                    ' after c.g.corr.:  !!RES!! = '//trim(sys_sdEL(dres,15)) )
              end if
                                            
              ! Perform the post-smoothing with the current solution vector
              if (associated(p_rsubnode%p_Rlevels(ilev)%p_rpostSmoother)) then
                call stat_startTimer (p_rsubnode%rtimeSmoothing)
                call sptils_smoothCorrection (&
                          p_rsubnode%p_Rlevels(ilev)%p_rpostSmoother,&
                          p_rsubnode%p_Rlevels(ilev)%rsolutionVector,&
                          p_rsubnode%p_Rlevels(ilev)%rrhsVector,&
                          p_rsubnode%p_Rlevels(ilev)%rtempVector,&
                          p_rsubnode%p_Rlevels(ilev)%rprjVector)
                call stat_stopTimer (p_rsubnode%rtimeSmoothing)

                call stat_addTimers (&
                    p_rsubnode%p_Rlevels(ilev)%p_rpostSmoother%rtimeSpacePrecond,&
                    rsolverNode%rtimeSpacePrecond)

                ! Extended output
                if ((rsolverNode%ioutputLevel .ge. 3) .and. &
                    (mod(ite,niteResOutput).eq.0)) then
                    
                  call sptivec_copyVector (p_rsubnode%p_Rlevels(ilev)%rrhsVector,&
                                          p_rsubnode%p_Rlevels(ilev)%rtempVector)
                  call cc_spaceTimeMatVec (rsolverNode%p_rproblem, p_rmatrix, &
                      p_rsubnode%p_Rlevels(ilev)%rsolutionVector,&
                      p_rsubnode%p_Rlevels(ilev)%rtempVector, -1.0_DP,1.0_DP,&
                      SPTID_FILTER_DEFECT,dres)

                  if (.not.((dres .ge. 1E-99_DP) .and. &
                            (dres .le. 1E99_DP))) dres = 0.0_DP
                            
                  call output_line ('Space-Time-Multigrid: Level '//trim(sys_siL(ilev,5))//&
                      ' after postsm.:    !!RES!! = '//trim(sys_sdEL(dres,15)) )
                end if
                                              
              end if

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
              if (p_rsubnode%p_Rlevels(ilev)%ncyclesRemaining .le. 0) then
                if (p_rsubnode%icycle .eq. 0) then
                  p_rsubnode%p_Rlevels(ilev)%ncyclesRemaining = 1
                else
                  ! Cycle finished. Reset counter for next cycle.
                  p_rsubnode%p_Rlevels(ilev)%ncyclesRemaining = &
                      p_rsubnode%p_Rlevels(ilev)%ncycles

                  if (((p_rsubnode%p_Rlevels(ilev)%depsRelCycle .ne. 1E99_DP) .or. &
                       (p_rsubnode%p_Rlevels(ilev)%depsRelCycle .ne. 1E99_DP)) .and. &
                      (ilev .lt. NLMAX)) then
                      
                    ! Adaptive cycles activated. 
                    !
                    ! We are on a level < nlmax.
                    ! At first, calculate the residuum on that level.
                    call sptivec_copyVector (p_rsubnode%p_Rlevels(ilev)%rrhsVector,&
                                            p_rsubnode%p_Rlevels(ilev)%rtempVector)
                    call cc_spaceTimeMatVec (rsolverNode%p_rproblem, p_rmatrix, &
                        p_rsubnode%p_Rlevels(ilev)%rsolutionVector,&
                        p_rsubnode%p_Rlevels(ilev)%rtempVector, -1.0_DP,1.0_DP,&
                        SPTID_FILTER_DEFECT,dres)

                    if (.not.((dres .ge. 1E-99_DP) .and. &
                              (dres .le. 1E99_DP))) dres = 0.0_DP
                              
                    ! Compare it with the initial residuum. If it's not small enough
                    ! and if we haven't reached the maximum number of cycles,
                    ! repeat the complete cycle.
                    if ( ((p_rsubnode%p_Rlevels(ilev)%nmaxAdaptiveCycles .le. -1) &
                          .or. &
                          (p_rsubnode%p_Rlevels(ilev)%icycleCount .lt. &
                           p_rsubnode%p_Rlevels(ilev)%nmaxAdaptiveCycles)) &
                        .and. &
                        ((dres .gt. p_rsubnode%p_Rlevels(ilev)%depsRelCycle * &
                                   p_rsubnode%p_Rlevels(ilev)%dinitResCycle) .and. &
                         (dres .gt. p_rsubnode%p_Rlevels(ilev)%depsAbsCycle) ) ) then

                      if (rsolverNode%ioutputLevel .ge. 3) then
                        call output_line ( &
                          trim( &
                          sys_siL(p_rsubnode%p_Rlevels(ilev)%icycleCount,10)) &
                          //'''th repetition of cycle on level '// &
                          trim(sys_siL(ilev,5))//'.')
                      end if

                      p_rsubnode%p_Rlevels(ilev)%icycleCount = &
                          p_rsubnode%p_Rlevels(ilev)%icycleCount+1
                      cycle cycleloop
                    end if
                    
                    ! Otherwise: The cycle(s) is/are finished; 
                    ! the END DO goes up one level.
                    
                  end if

                end if
              else

                if (rsolverNode%ioutputLevel .ge. 3) then
                  call output_line ('Space-Time-Multigrid: Cycle on level '&
                      //trim(sys_siL(ilev,5))//' finished.')
                end if

                ! Next cycle; go down starting from the current level
                cycle cycleloop

              end if
              
            end do ! ilev < nlmax
            
            ! We finally reached the maximum level. Quit the cycle loop
            ! to go on with the next iteration.
            exit cycleloop
          
          end do cycleloop
          
          ! We have (hopefully) successfully performed one MG-sweep, starting
          ! and ending on the finest level. As we are now on the finest level
          ! again, we can update our defect vector to test the current
          ! residuum...
          !
          ! Calculate the residuum and its norm.
          call sptivec_copyVector (p_rsubnode%p_Rlevels(ilev)%rrhsVector,&
                                    p_rsubnode%p_Rlevels(ilev)%rtempVector)
          call cc_spaceTimeMatVec (rsolverNode%p_rproblem, p_rmatrix, &
              p_rsubnode%p_Rlevels(ilev)%rsolutionVector,&
              p_rsubnode%p_Rlevels(ilev)%rtempVector, -1.0_DP,1.0_DP,&
              SPTID_FILTER_DEFECT,dres)

          if (.not.((dres .ge. 1E-99_DP) .and. &
                    (dres .le. 1E99_DP))) dres = 0.0_DP
          
          rsolverNode%dlastDefect = rsolverNode%dfinalDefect
          rsolverNode%dfinalDefect = dres
          
          ! Test if the iteration is diverged
          if (sptils_testDivergence(rsolverNode,dres)) then
            call output_line ('Space-Time-Multigrid: Solution diverging!')
            rsolverNode%iresult = 1
            exit
          end if

          ! At least perform nminIterations iterations
          if (ite .ge. nminIterations) then
          
            ! Check if the iteration converged
            if (sptils_testConvergence(rsolverNode,dres)) exit
            
          end if

          ! print out the current residuum

          if ((rsolverNode%ioutputLevel .ge. 2) .and. &
              (mod(ite,niteResOutput).eq.0)) then
            call output_line ('Space-Time-Multigrid: Iteration '// &
                trim(sys_siL(ITE,10))//',  !!RES!! = '//&
                trim(sys_sdEL(rsolverNode%dfinalDefect,15)) )
          end if
          
        end do  ! ite
        
        ! Set ITE to NIT to prevent printing of "NIT+1" of the loop was
        ! completed

        if (ite .gt. rsolverNode%nmaxIterations) &
          ite = rsolverNode%nmaxIterations

        ! Finish - either with an error or if converged.
        ! Print the last residuum.

        if ((rsolverNode%ioutputLevel .ge. 2) .and. &
            (ite .ge. 1) .and. (ITE .lt. rsolverNode%nmaxIterations) .and. &
            (rsolverNode%iresult .ge. 0)) then
          call output_line ('Space-Time-Multigrid: Iteration '// &
              trim(sys_siL(ITE,10))//',  !!RES!! = '//&
              trim(sys_sdEL(rsolverNode%dfinalDefect,15)) )
        end if
        
        ! Final number of iterations
        rsolverNode%iiterations = ite
      
      end if

      ! Finally, we look at some statistics:
      !
      ! Don't calculate anything if the final residuum is out of bounds -
      ! would result in NaN's,...
        
      if (rsolverNode%dfinalDefect .lt. 1E99_DP) then
      
        ! If the initial defect was zero, the solver immediately
        ! exits - and so the final residuum is zero and we performed
        ! no steps; so the resulting multigrid convergence rate stays zero.
        ! In the other case the multigrid convergence rate computes as
        ! (final defect/initial defect) ** 1/nit :

        if (rsolverNode%dinitialDefect .gt. rsolverNode%drhsZero) then
          rsolverNode%dconvergenceRate = &
                      (rsolverNode%dfinalDefect / rsolverNode%dinitialDefect) ** &
                      (1.0_DP/real(rsolverNode%iiterations,DP))
        end if
      end if
    
      ! Calculate the time for linear algebra
      call stat_stopTimer (p_rsubnode%rtimeLinearAlgebra)
      call stat_subTimers (p_rsubnode%rtimeSmoothing,p_rsubnode%rtimeLinearAlgebra)
      call stat_subTimers (p_rsubnode%rtimeCoarseGridSolver,p_rsubnode%rtimeLinearAlgebra)
      call stat_subTimers (p_rsubnode%rtimeProlRest,p_rsubnode%rtimeLinearAlgebra)

    end if
    
    ! As the solution vector on the finest level shared its memory with rd,
    ! we just calculated the new correction vector!
      
    if (rsolverNode%dfinalDefect .lt. 1E99_DP) then
      
      if (rsolverNode%ioutputLevel .ge. 2) then
        call output_lbrk()
        call output_line ('Space-Time-Multigrid statistics:')
        call output_lbrk()
        call output_line ('Iterations              : '//&
             trim(sys_siL(rsolverNode%iiterations,10)) )
        call output_line ('!!INITIAL RES!!         : '//&
             trim(sys_sdEL(rsolverNode%dinitialDefect,15)) )
        call output_line ('!!RES!!                 : '//&
             trim(sys_sdEL(rsolverNode%dfinalDefect,15)) )
        if (rsolverNode%dinitialDefect .gt. rsolverNode%drhsZero) then     
          call output_line ('!!RES!!/!!INITIAL RES!! : '//&
            trim(sys_sdEL(rsolverNode%dfinalDefect / rsolverNode%dinitialDefect,15)) )
        else
          call output_line ('!!RES!!/!!INITIAL RES!! : '//&
               trim(sys_sdEL(0.0_DP,15)) )
        end if
        call output_lbrk ()
        call output_line ('Rate of convergence     : '//&
             trim(sys_sdEL(rsolverNode%dconvergenceRate,15)) )

      end if

      if (rsolverNode%ioutputLevel .eq. 1) then
        call output_line (&
              'MG: Iterations/Rate of convergence: '//&
              trim(sys_siL(rsolverNode%iiterations,10))//' /'//&
              trim(sys_sdEL(rsolverNode%dconvergenceRate,15)) )
      end if

    else
      ! DEF=Infinity; RHO=Infinity, set to 1
      rsolverNode%dconvergenceRate = 1.0_DP
    end if  
  
  contains
  
    ! ---------------------------------------------------------------
  
    subroutine calcCGCenergyMin (rproblem,&
                    rsolutionVector,rcorrectionVector,rrhsVector,rmatrix,&
                    rtempVector,rtempSpaceVector,dalphamin, dalphamax, dstep)
                  
    ! Calculates the optimal coarse grid correction factor using
    ! energy minimisation.
    
    ! Problem structure
    type(t_problem), intent(INOUT) :: rproblem
    
    ! Uncorrected solution vector
    type(t_spacetimeVector), intent(IN) :: rsolutionVector
    
    ! Correction vector
    type(t_spacetimeVector), intent(IN) :: rcorrectionVector

    ! RHS vector
    type(t_spacetimeVector), intent(IN) :: rrhsVector
    
    ! Corresponding space-time matrix
    type(t_ccoptSpaceTimeMatrix), intent(IN) :: rmatrix
    
    ! Temp vector
    type(t_spacetimeVector), intent(INOUT) :: rtempVector
    
    ! A temp vector, only in space
    type(t_vectorBlock), intent(INOUT) :: rtempSpaceVector
    
    ! Minimum/maximum value for the correction factor
    real(DP), intent(IN) :: dalphamin,dalphamax
    
    ! OUTPUT: Calculated step length parameter
    real(DP), intent(OUT) :: dstep
    
      ! local variables
      real(DP) :: dnom, ddenom
    
      ! The formula is just as for standard multigrid.
      ! With u=uncorrected solution, c=correction vector, 
      ! f=rhs vector, calculate
      !
      !           ( f-Au, c )
      !  dstep = -------------
      !            ( Ac, c )
      !
      ! Calculate the nominator
      call sptivec_copyVector (rrhsVector,rtempVector)
      call cc_spaceTimeMatVec (rproblem, rmatrix, &
          rsolutionVector,rtempVector, -1.0_DP,1.0_DP,SPTID_FILTER_DEFECT)

      ! Implement boundary conditions into the vector.
      !CALL tbc_implementInitCondDefect (&
      !    rmatrix%p_rspaceTimeDiscr,&
      !    rtempVector, rtempSpaceVector)
      !CALL tbc_implementBCdefect (rproblem,&
      !    rmatrix%p_rspaceTimeDiscr,&
      !    rtempVector, rtempSpaceVector)

      ! Get the nominator
      dnom = sptivec_scalarProduct (rtempVector,rcorrectionVector)

      ! Calculate the denominator:
      
      call sptivec_clearVector (rtempVector)
      call cc_spaceTimeMatVec (rproblem, rmatrix, &
          rcorrectionVector,rtempVector, 1.0_DP,0.0_DP,SPTID_FILTER_DEFECT)
      
      ! Implement boundary conditions into the vector.
      !CALL tbc_implementInitCondDefect (&
      !    rmatrix%p_rspaceTimeDiscr,&
      !    rtempVector, rtempSpaceVector)
      !CALL tbc_implementBCdefect (rproblem,&
      !    rmatrix%p_rspaceTimeDiscr,&
      !    rtempVector, rtempSpaceVector)

      ! Get the denominator      
      ddenom = sptivec_scalarProduct (rtempVector,rcorrectionVector)
      
      ! Trick to avoid div/0
      if (ddenom .eq. 0.0_DP) ddenom = 1.0_DP
      
      ! Result...
      dstep = dnom / ddenom
                   
      ! Force it to be in the allowed range.
      dstep = min(dalphamax,max(dalphamin,dstep))
                    
    end subroutine
  
    ! ---------------------------------------------------------------
  
    subroutine calcCGCenergyMinWeighted (rproblem,&
                    rsolutionVector,rcorrectionVector,rrhsVector,rmatrix,&
                    rtempVector,rtempSpaceVector,dalphamin, dalphamax, &
                    Dweights,dstep)
                  
    ! Calculates the optimal coarse grid correction factor using
    ! energy minimisation. Use weighted energy minimisation that
    ! multiplies the residualy by a given factor before setting up
    ! the correction factor.
    
    ! Problem structure
    type(t_problem), intent(INOUT) :: rproblem
    
    ! Uncorrected solution vector
    type(t_spacetimeVector), intent(IN) :: rsolutionVector
    
    ! Correction vector
    type(t_spacetimeVector), intent(IN) :: rcorrectionVector

    ! RHS vector
    type(t_spacetimeVector), intent(IN) :: rrhsVector
    
    ! Corresponding space-time matrix
    type(t_ccoptSpaceTimeMatrix), intent(IN) :: rmatrix
    
    ! Temp vector
    type(t_spacetimeVector), intent(INOUT) :: rtempVector
    
    ! A temp vector, only in space
    type(t_vectorBlock), intent(INOUT) :: rtempSpaceVector
    
    ! Minimum/maximum value for the correction factor
    real(DP), intent(IN) :: dalphamin,dalphamax
    
    ! A weight factor for each equation of a timestep. This is applied
    ! to all timesteps. The weighting factor is multiplied by the residual
    ! before taking scalar products. Standard value is an array (1,1,1,...).
    real(DP), dimension(:), intent(IN) :: Dweights
    
    ! OUTPUT: Calculated step length parameter
    real(DP), intent(OUT) :: dstep
    
      ! local variables
      real(DP) :: dnom, ddenom
    
      ! The formula is just as for standard multigrid.
      ! With u=uncorrected solution, c=correction vector, 
      ! f=rhs vector, calculate
      !
      !           ( f-Au, c )
      !  dstep = -------------
      !            ( Ac, c )
      !
      ! Calculate the nominator
      call sptivec_copyVector (rrhsVector,rtempVector)
      call cc_spaceTimeMatVec (rproblem, rmatrix, &
          rsolutionVector,rtempVector, -1.0_DP,1.0_DP,SPTID_FILTER_DEFECT)

      ! Implement boundary conditions into the vector.
      !CALL tbc_implementInitCondDefect (&
      !    rmatrix%p_rspaceTimeDiscr,&
      !    rtempVector, rtempSpaceVector)
      !CALL tbc_implementBCdefect (rproblem,&
      !    rmatrix%p_rspaceTimeDiscr,&
      !    rtempVector, rtempSpaceVector)

      ! Get the nominator
      dnom = sptivec_scalarProductWeighted (rtempVector,rcorrectionVector,Dweights)
      ! Calculate the denominator:
      
      call sptivec_clearVector (rtempVector)
      call cc_spaceTimeMatVec (rproblem, rmatrix, &
          rcorrectionVector,rtempVector, 1.0_DP,0.0_DP,SPTID_FILTER_DEFECT)
      
      ! Implement boundary conditions into the vector.
      !CALL tbc_implementInitCondDefect (&
      !    rmatrix%p_rspaceTimeDiscr,&
      !    rtempVector, rtempSpaceVector)
      !CALL tbc_implementBCdefect (rproblem,&
      !    rmatrix%p_rspaceTimeDiscr,&
      !    rtempVector, rtempSpaceVector)

      ! Get the denominator      
      ddenom = sptivec_scalarProductWeighted (rtempVector,rcorrectionVector,Dweights)
      
      ! Trick to avoid div/0
      if (ddenom .eq. 0.0_DP) ddenom = 1.0_DP
      
      ! Result...
      dstep = dnom / ddenom
                   
      ! Force it to be in the allowed range.
      dstep = min(dalphamax,max(dalphamin,dstep))
                    
    end subroutine

  end subroutine

end module
