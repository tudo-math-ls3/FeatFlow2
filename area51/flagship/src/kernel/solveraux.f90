!##############################################################################
!# ****************************************************************************
!# <name> solveraux </name>
!# ****************************************************************************
!#
!# <purpose>
!#
!# This module provides the basic data structures and subroutines for
!# handling both linear/nonlinear and one-level/multilevel solvers.
!# A calling program should create a new solver structure and fill the
!# required parameters either by reading a user-supplied parameters or
!# via direct adjustment.
!#
!# Furthermore, each structure of type t_solver provides information
!# about convergence and norms after the solver has terminated.
!#
!# The following routines are available:
!#
!# ----------------------------------------------------------------------------
!#
!# 1.) solver_createSolver = solver_createSolverDirect /
!#                           solver_createSolverIndirect
!#     -> Creates a new solver from parameter list
!#
!# 2.) solver_releaseSolver
!#     -> Releases an existing solver structure
!#
!# 3.) solver_infoSolver
!#     -> Outputs information about the solver structure
!#
!# 4.) solver_resetSolver
!#      -> Resets the solver structure to initial values
!#
!# 5.) solver_removeTempFromSolver
!#     -> Removes temporal storage from solver structure
!#
!# 6.) solver_adjustHierarchy
!#     -> Adjusts the internal hierarchy of the complete solver structure
!#
!# 7.) solver_showHierarchy
!#     -> Shows the solver hierarchy in human readable form
!#
!# 8) solver_statistics
!#    -> Computes statistics for a solver structure
!#
!# 9.) solver_getStatus
!#     -> Returns the status of solver structure in human readable form
!#
!# 10.) solver_getNextSolver
!#      -> Gets next solver structure in hierarchy
!#
!# 11.) solver_getNextSolverByType
!#      -> Returns next solver with specified type
!#
!# 12.) solver_getMinimumMultigridlevel
!#      -> Returns the minimum multigrid level required in the whole structure
!#
!# 13.) solver_getMaximumMultigridlevel
!#      -> Returns the maximum multigrid level required in the whole structure
!#
!# 14.) solver_setMinimumMultigridlevel
!#      -> Sets the minimum multigrid level required in the whole structure
!#
!# 15.) solver_setMaximumMultigridlevel
!#      -> Sets the maximum multigrid level required in the whole structure
!#
!# 16.) solver_setBoundaryCondition
!#      -> Associates boundary conditions structure to the solver
!#
!# 17.) solver_setSolverMatrix = solver_setSolverMatrixScalar /
!#                               solver_setSolverMatrixBlock
!#      -> Associates a matrix to the solver
!#
!# 18.) solver_setPrecondMatrix = solver_setPrecondMatrixScalar /
!#                                solver_setPrecondMatrixBlock
!#      -> Associates a matrix to the preconditioner
!#
!# 19.) solver_setSmootherMatrix = solver_setSmootherMatrixScalar /
!#                                 solver_setSmootherMatrixBlock
!#      -> Associates a matrix to the smoother
!#
!# 20.) solver_updateStructure
!#      -> Updates the structure of the solver recursively
!#
!# 21.) solver_updateContent
!#      -> Updates the content of the solver recursively
!#
!# 22.) solver_prolongationScalar
!#      -> Prolongates some scalar-vector from coarse to fine grid
!#
!# 23.) solver_prolongationBlock
!#      -> Prolongates some block-vector from coarse to fine grid
!#
!# 24.) solver_restrictionScalar
!#      -> Restricts some scalar-vector from fine to coarse grid
!#
!# 25.) solver_restrictionBlock
!#      -> Restricts some block-vector from fine to coarse grid
!#
!# 26.) solver_copySolver
!#      -> Copies input/output parameters from one solver to another solver
!#
!# 27.) solver_testConvergence
!#      -> Checks all convergence criteria for a particular solver
!#
!# 28.) solver_testDivergence
!#      -> Checks all divergence criteria for a particular solver
!#
!# 29.) solver_testStagnation
!#      -> Checks all stagnation criteria for a particular solver
!#
!# The following auxiliary routines are available:
!#
!# 1.) solver_initUMFPACK
!#     -> Initialize the UMFPACK preconditioner
!#
!# 2.) solver_initILU
!#     -> Initialize the ILU-preconditioner
!#
!# </purpose>
!##############################################################################

module solveraux

  use boundarycondaux
  use boundaryfilter
  use fsystem
  use genoutput
  use linearalgebra
  use linearsystemblock
  use linearsystemscalar
  use paramlist
  use storage
  use triangulation

  implicit none

  private
  public :: t_solver
  public :: t_solverMultigrid
  public :: t_solverUMFPACK
  public :: t_solverJacobi
  public :: t_solverSSOR
  public :: t_solverBiCGSTAB
  public :: t_solverGMRES
  public :: t_solverILU
  public :: t_solverDefcor
  public :: t_solverNewton

  public :: solver_createSolver
  public :: solver_releaseSolver
  public :: solver_resetSolver
  public :: solver_infoSolver
  public :: solver_copySolver
  public :: solver_removeTempFromSolver
  public :: solver_adjustHierarchy
  public :: solver_showHierarchy
  public :: solver_statistics
  public :: solver_getStatus
  public :: solver_getNextSolver
  public :: solver_getNextSolverByType
  public :: solver_getNextSolverByTypes
  public :: solver_getMinimumMultigridlevel
  public :: solver_getMaximumMultigridlevel
  public :: solver_setMinimumMultigridlevel
  public :: solver_setMaximumMultigridlevel
  public :: solver_setBoundaryCondition
  public :: solver_setSolverMatrix
  public :: solver_setPrecondMatrix
  public :: solver_setSmootherMatrix
  public :: solver_updateStructure
  public :: solver_updateContent
  public :: solver_prolongationScalar
  public :: solver_prolongationBlock
  public :: solver_restrictionScalar
  public :: solver_restrictionBlock
  public :: solver_testConvergence
  public :: solver_testDivergence
  public :: solver_testStagnation

  ! ****************************************************************************

  interface solver_createSolver
    module procedure solver_createSolverDirect
    module procedure solver_createSolverIndirect
  end interface

  interface solver_setSolverMatrix
    module procedure solver_setSolverMatrixScalar
    module procedure solver_setSolverMatrixBlock
  end interface

  interface solver_setPrecondMatrix
    module procedure solver_setPrecondMatrixScalar
    module procedure solver_setPrecondMatrixBlock
  end interface

  interface solver_setSmootherMatrix
    module procedure solver_setSmootherMatrixScalar
    module procedure solver_setSmootherMatrixBlock
  end interface

  ! ****************************************************************************

!<constants>

!<constantblock description="Global solver types">

  ! Note: Solver types are the most general descriptions of a solver.
  ! This module only supports linear/nonlinear single-/multigrid solvers

  ! Undefined solver type
  integer, parameter, public :: SV_UNDEFINED    = 0

  ! Nonlinear solver type
  integer, parameter, public :: SV_NONLINEAR    = 1

  ! Linear solver type
  integer, parameter, public :: SV_LINEAR       = 2

  ! Nonlinear multigrid solver type
  integer, parameter, public :: SV_NONLINEARMG  = 3

  ! Linear multigrid solver type
  integer, parameter, public :: SV_LINEARMG     = 4

  ! Full multigrid solver for steady state computations
  integer, parameter, public :: SV_FMG          = 5

  ! Block of coupled solvers
  integer, parameter, public :: SV_COUPLED      = 6
!</constantblock>

  ! ****************************************************************************

!<constantblock description="Global flags for information output">

  ! Silent run, no output at all
  integer, parameter, public :: SV_IOLEVEL_SILENT  = 0

  ! Output only errors
  integer, parameter, public :: SV_IOLEVEL_ERROR   = 1

  ! Output only errors and warnings
  integer, parameter, public :: SV_IOLEVEL_WARNING = 2

  ! Output errors, warnings and information
  integer, parameter, public :: SV_IOLEVEL_INFO    = 3

  ! Output errors, warnings and verbose information
  integer, parameter, public :: SV_IOLEVEL_VERBOSE = 4

!</constantblock>

  ! ****************************************************************************

!<constantblock description="Global status messages of the solver">

  ! Solver converged
  integer, parameter, public :: SV_CONVERGED    = 0

  ! Solver detected zero right-hand side
  integer, parameter, public :: SV_ZERO_RHS     = 1

  ! Solver detected zero defect
  integer, parameter, public :: SV_ZERO_DEF     = 2

  ! Solver detected too large defect
  integer, parameter, public :: SV_INF_DEF      = 3

  ! Solver detected too large increase of defect
  integer, parameter, public :: SV_INCR_DEF     = 4

  ! Solver diverged
  integer, parameter, public :: SV_DIVERGED     = 5

  ! Solver stagnets
  integer, parameter, public :: SV_STAGNATED    = 6
!</constantblock>

  ! ****************************************************************************

!<constantblock description="Identifiers for stopping criterium istoppingCriterion.">

  ! Use 'maximum' stopping criterion.
  ! If depsRel>0: use relative stopping criterion.
  ! If depsAbs>0: use abs stopping criterion.
  ! If both are > 0: use both, i.e. the iteration stops when both,
  !    the relative AND the absolute stopping criterium holds
  integer, parameter, public :: SV_STOP_ALL = 0

  ! Use 'minimum' stopping criterion.
  ! If depsRel>0: use relative stopping criterion.
  ! If depsAbs>0: use abs stopping criterion.
  ! If both are > 0: use one of them, i.e. the iteration stops when the
  !    either the relative OR the absolute stopping criterium holds
  integer, parameter, public :: SV_STOP_ANY = 1

!</constantblock>

  ! ****************************************************************************

!<constantblock description="Flags for the solver specification bitfield">

  ! Solver structure needs update
  integer(I32), parameter :: SV_SSPEC_STRUCTURENEEDSUPDATE = 2**0

  ! Solver content needs update
  integer(I32), parameter :: SV_SSPEC_CONTENTNEEDSUPDATE   = 2**1

  ! Solver needs update (both structure and content)
  integer(I32), parameter :: SV_SSPEC_NEEDSUPDATE          = SV_SSPEC_STRUCTURENEEDSUPDATE+&
                                                             SV_SSPEC_CONTENTNEEDSUPDATE

  ! Empty solver structure
  integer(I32), parameter :: SV_SSPEC_EMPTY                = SV_SSPEC_NEEDSUPDATE
!</constantblock>

  ! ****************************************************************************

!<constantblock description="Global nonlinear solution algorithms">

  ! Fixed-point iteration
  integer, parameter, public :: NLSOL_SOLVER_FIXEDPOINT    = 101
!</constantblock>

  ! ****************************************************************************

!<constantblock description="Global nonlinear preconditioners">

  ! Preconditioning by segregated/block-diagonal approach
  integer, parameter, public :: NLSOL_PRECOND_BLOCKD        = 1

  ! Preconditioning by defect-correction approach
  integer, parameter, public :: NLSOL_PRECOND_DEFCOR        = 2

  ! Preconditioning by (inexact) Newton approach
  integer, parameter, public :: NLSOL_PRECOND_NEWTON        = 3
  integer, parameter, public :: NLSOL_PRECOND_NEWTON_FAILED = -3
!</constantblock>

  ! ****************************************************************************

!<constantblock description="Global nonlinear smoothers">

  ! Fixed-point iteration
  integer, parameter, public :: NLSOL_SMOOTHER_FIXEDPOINT   = NLSOL_SOLVER_FIXEDPOINT
!</constantblock>

  ! ****************************************************************************

!<constantblock description="Flags for the nonlinear solver operation specification bitfield">

  ! Nonlinear solver requires calculation of preconditioner
  integer(I32), parameter, public :: NLSOL_OPSPEC_CALCPRECOND    = 2**0

  ! Nonlinear solver requires calculation of constant rhs vector
  integer(I32), parameter, public :: NLSOL_OPSPEC_CALCRHS        = 2**1

  ! Nonlinear solver requires calculation of residual vector
  integer(I32), parameter, public :: NLSOL_OPSPEC_CALCRESIDUAL   = 2**2

  ! Nonlinear solver requires calculation of Jacobian operator
  integer(I32), parameter, public :: NLSOL_OPSPEC_CALCJACOBIAN   = 2**3

  ! Nonlinear solver requires application of Jacobian operator
  integer(I32), parameter, public :: NLSOL_OPSPEC_APPLYJACOBIAN  = 2**4
!</constantblock>

  ! ****************************************************************************

!<constantblock description="Global linear solution algorithms">

  ! no solver
  integer, parameter, public :: LINSOL_SOLVER_NONE     = 0

  ! (block-)Jacobi iteration $x_1 = x_0 + \omega D^{-1} (b-Ax_0)$
  integer, parameter, public :: LINSOL_SOLVER_JACOBI   = 2

  ! (block-)SOR/GS iteration $x_1 = x_0 + (L+\omega D)^{-1}(b-Ax_0)$
  integer, parameter, public :: LINSOL_SOLVER_SOR      = 4

  ! (block-)SSOR iteration
  integer, parameter, public :: LINSOL_SOLVER_SSOR     = 5

  ! BiCGStab iteration (preconditioned)
  integer, parameter, public :: LINSOL_SOLVER_BICGSTAB = 7

  ! GMRES iteration (preconditioned)
  integer, parameter, public :: LINSOL_SOLVER_GMRES    = 8

  ! Identifier for Umfpack4
  integer, parameter, public :: LINSOL_SOLVER_UMFPACK4 = 11

  ! MILU(s) solver (scalar system)
  integer, parameter, public :: LINSOL_SOLVER_ILU      = 50
!</constantblock>

  ! ****************************************************************************

!<constantblock description="Global linear smoothers">

  ! no smoother
  integer, parameter, public :: LINSOL_SMOOTHER_NONE     = LINSOL_SOLVER_NONE

  ! (block-)Jacobi smoother
  integer, parameter, public :: LINSOL_SMOOTHER_JACOBI   = LINSOL_SOLVER_JACOBI

  ! (block-)SOR/GS smoother
  integer, parameter, public :: LINSOL_SMOOTHER_SOR      = LINSOL_SOLVER_SOR

  ! (bock-)BiCGSTAB smoother
  integer, parameter, public :: LINSOL_SMOOTHER_BiCGSTAB = LINSOL_SOLVER_BiCGSTAB

  ! (block-)GMRES smoother
  integer, parameter, public :: LINSOL_SMOOTHER_GMRES    = LINSOL_SOLVER_GMRES

  ! (block-)SSOR smoother
  integer, parameter, public :: LINSOL_SMOOTHER_SSOR     = LINSOL_SOLVER_SSOR

  ! ILU(0) smoother (scalar system)
  integer, parameter, public :: LINSOL_SMOOTHER_ILU      = LINSOL_SOLVER_ILU
!</constantblock>

!</constants>

  ! *****************************************************************************

!<types>

!<typeblock>

  ! This data structure contains all settings/parameters for linear/nonlinear
  ! single/multigrid-solvers. Hence, it is the most abstract solver structure.
  !
  ! For single grid solvers, it contains all specifig solver and preconditioner
  ! nodes which need to be allocated/deallocated depending on the solver/
  ! preconditioner type. In addition, it contains a subnode p_solverMultigrid
  ! which has to be allocated/deallocated in case the solver is applied
  ! to multiple grid levels.

  type t_solver

    ! INPUT PARAMETER FOR ITERATIVE SOLVERS:
    ! Identifier of the solver
    character(SYS_STRLEN) :: ssolverName = ''

    ! INPUT PARAMETER FOR ITERATIVE SOLVERS
    ! Identifies the type of the solver
    ! Valid values: SV_UNDEFINED, SV_FMG, SV_NONLINEARMG,
    !               SV_NONLINEAR, SV_LINEARMG, SV_LINEAR
    integer :: csolverType = SV_UNDEFINED

    ! INPUT PARAMETER FOR ITERATIVE SOLVERS
    ! Solver specification tag. This is a bitfield coming from an
    ! OR combination of difference SV_SSPEC_xxxx constants and it
    ! specifies various details of the solver structure.
    integer(I32) :: isolverSpec = SV_SSPEC_EMPTY

    ! INPUT PARAMETER FOR ITERATIVE SOLVERS
    ! Type of solver
    integer :: isolver = SV_UNDEFINED

    ! INPUT PARAMETER FOR ITERATIVE SOLVERS
    ! Type of preconditioner
    ! Valid values: SV_UNDEFINED, NLSOL_PRECOND_BLOCKD,
    !               NLSOL_PRECOND_DEFCOR, NLSOL_PRECOND_NEWTON
    integer :: iprecond = SV_UNDEFINED

    ! INPUT PARAMETER FOR ITERATIVE SOLVERS
    ! Information output level
    integer :: ioutputLevel = 0

    ! INPUT PARAMETER FOR ITERATIVE SOLVERS
    integer :: istoppingCriterion = SV_STOP_ALL

    ! OUTPUT PARAMETER FOR ITERATIVE SOLVERS:
    ! Status of the solver
    integer :: istatus = 0

    ! INPUT PARAMETER FOR ITERATIVE SOLVERS
    ! Type of norm to use in the residual checking
    integer :: iresNorm = 0

    ! INPUT PARAMETER FOR ITERATIVE SOLVERS:
    ! Minimum number of iterations top perform
    integer :: nminIterations = 0

    ! INPUT PARAMETER FOR ITERATIVE SOLVERS:
    ! Maximum number of iterations top perform
    integer :: nmaxIterations = 0

    ! OUTPUT PARAMETER FOR ITERATIVE SOLVERS:
    ! Number of performed iterations
    integer :: iiterations = 0

    ! OUTPUT PARAMETER FOR ITERATIVE SOLVERS:
    ! Total number of performed iterations
    integer :: niterations = 0

    ! INPUT PARAMETER FOR ITERATIVE SOLVERS:
    ! Relaxation factor; standard value = 1.0
    ! If DOMEGA < 0 than implicit under-relaxation is perform by
    ! ABS(DOMEGA). Otherwise, the value DOMEGA is used to perform
    ! explicit under-relaxation of the system.
    real(DP) :: domega = 1.0_DP

    ! INPUT PARAMETER FOR ITERATIVE SOLVERS:
    ! Relative stopping criterion. Stop iteration if
    ! !!defect!! < EPSREL * !!initial defect!!.
    ! =0: ignore, use absolute stopping criterion
    ! Remark: do not set depsAbs=depsRel=0!
    real(DP) :: depsRel = 0.0_DP

    ! INPUT PARAMETER FOR ITERATIVE SOLVERS:
    ! Absolute stopping criterion. Stop iteration if
    ! !!defect!! < EPSREL.
    ! =0: ignore, use relative stopping criterion
    ! Remark: do not set depsAbs=depsRel=0!
    real(DP) :: depsAbs = 0.0_DP

    ! INPUT PARAMETER FOR ITERATIVE SOLVERS:
    ! Relative divergence criterion. Treat iteration as
    ! diverged if
    !   !!defect!! >= DIVREL * !!initial defect!!
    ! A value of SYS_INFINITY disables the relative divergence check.
    real(DP) :: ddivRel = SYS_INFINITY

    ! INPUT PARAMETER FOR ITERATIVE SOLVERS:
    ! Absolute divergence criterion. Treat iteration as
    ! diverged if
    !   !!defect!! >= DIVREL
    ! A value of SYS_INFINITY disables the absolute divergence check.
    real(DP) :: ddivAbs = SYS_INFINITY

    ! INPUT PARAMETER FOR ITERATIVE SOLVERS:
    ! RHS-vector is treated as zero if max(rhs) < drhsZero
    real(DP) :: drhsZero = 0.0_DP

    ! INPUT PARAMETER FOR ITERATIVE SOLVERS:
    ! Defect-vector is treated as zero if !!defect!! < ddefZero
    real(DP) :: ddefZero = 0.0_DP

    ! INPUT PARAMETER FOR ITERATIVE SOLVERS:
    ! Relative stopping criterion. Stop iteration of
    ! !!current solution - last solution!! <= EPSSTAG
    ! =0: ignore, do not check for stagnation
    real(DP) :: depsStag = 0.0_DP

    ! OUTPUT PARAMETER FOR ITERATIVE SOLVERS:
    ! Norm of initial right-hand side
    real(DP) :: dinitialRHS = 0.0_DP

    ! OUTPUT PARAMETER FOR ITERATIVE SOLVERS:
    ! Norm of initial residuum
    real(DP) :: dinitialDefect = 0.0_DP

    ! OUTPUT PARAMETER FOR ITERATIVE SOLVERS:
    ! Norm of final residuum
    real(DP) :: dfinalDefect = 0.0_DP

    ! OUTPUT PARAMETER FOR ITERATIVE SOLVERS:
    ! Convergence rate
    real(DP) :: dconvergenceRate = 0.0_DP

    ! INPUT: system boundary, this boundary condition is used in
    ! all nonlinear/linear iterations to impose boundary conditions
    type(t_boundaryCondition), pointer :: rboundaryCondition => null()

    ! INTERNAL: substructure for generic solvers
    type(t_solver), dimension(:), pointer :: p_solverSubnode => null()

    ! INTERNAL: substructure for Multigrid solver
    type(t_solverMultigrid), pointer :: p_solverMultigrid => null()

    ! INTERNAL: substructure for UMFPACK solver
    type(t_solverUMFPACK), pointer :: p_solverUMFPACK => null()

    ! INTERNAL: substructure for Jacobi solver
    type(t_solverJacobi), pointer :: p_solverJacobi => null()

    ! INTERNAL: substructure for (S)SOR solver
    type(t_solverSSOR), pointer :: p_solverSSOR => null()

    ! INTERNAL: substructure for BiCGSTAB solver
    type(t_solverBiCGSTAB), pointer :: p_solverBiCGSTAB => null()

    ! INTERNAL: substructure for GMRES solver
    type(t_solverGMRES), pointer :: p_solverGMRES => null()

    ! INTERNAL: substructure for ILU solver
    type(t_solverILU), pointer :: p_solverILU => null()

    ! INTERNAL: substructure for Defect correction
    type(t_solverDefcor), pointer :: p_solverDefcor => null()

    ! INTERNAL: substructure for Newton`s method
    type(t_solverNewton), pointer :: p_solverNewton => null()
  end type t_solver

!</typeblock>

  ! *****************************************************************************

!<typeblock>

  ! This data structure contains all local data for a multigrid solver

  type t_solverMultigrid

    ! INPUT: minimum multigrid level
    integer :: nlmin = 0

    ! INTERNAL: initial setting of NLMIN
    integer :: initialNlmin = 0

    ! INPUT: maximum multigrid level
    integer :: nlmax = 0

    ! INTERNAL: initial setting if NLMAX
    integer :: initialNlmax = 0

    ! INPUT: minimum number of multigrid steps
    integer :: ilmin = 0

    ! INPUT: maximum number of multigrid steps
    integer :: ilmax = 0

    ! INPUT: type of multigrid cycle
    integer :: icycle = 0

    ! INPUT: Identifier of the smoother
    character(SYS_STRLEN) :: ssmootherName = ''

    ! INPUT: Identifies the type of smoother
    ! Valid values: SV_UNDEFINED, SV_NONLINEAR, SV_LINEAR
    integer :: csmootherType = SV_UNDEFINED

    ! INPUT: number of pre-smoothing steps
    integer :: npresmooth = 0

    ! INPUT: number of post-smoothing steps
    integer :: npostsmooth = 0

    ! INPUT: factor for smoothing steps
    integer :: nsmoothFactor = 0

    ! INTERNAL: structure for coarsegrid solver
    type(t_solver), pointer :: p_solverCoarsegrid => null()

    ! INTERNAL: structure for smoothing
    type(t_solver), dimension(:), pointer :: p_smoother => null()

    ! INTERNAL: system matrices
    type(t_matrixBlock), dimension(:), pointer :: rmatrix => null()

    ! INTERNAL: temporal vectors
    type(t_vectorBlock), dimension(:), pointer :: RtempVectors => null()
  end type t_solverMultigrid

!</typeblock>

  ! *****************************************************************************

!<typeblock>

  ! This data structure contains all local data for direct solving via UMFPACK

  type t_solverUMFPACK

    ! Cotrol structure for UMFPACK4; contains parameters for the solver
    real(DP), dimension(20) :: Dcontrol

    ! Handle for symbolic factorisation
    ! This is not a FEAT-Handle!
    integer(I32) :: isymbolic = 0

    ! Array of handles for symbolic factorisation
    ! These are no FEAT-Handles!
    integer(I32), dimension(:), pointer :: IsymbolicBlock => null()

    ! Handle for numeric factorisation
    ! This is not a FEAT-Handle!
    integer(I32) :: inumeric = 0

    ! Array of handles for numeric factorisation
    ! These are no FEAT-Handles!
    integer(I32), dimension(:), pointer :: InumericBlock => null()

    ! INTERNAL: system matrix
    type(t_matrixBlock) :: rmatrix

    ! INTERNAL: temporal vectors
    type(t_vectorBlock) :: rtempVector
  end type t_solverUMFPACK

!</typeblock>

  ! *****************************************************************************

!<typeblock>

  ! This data structure contains all local data for a Jacobi algorithm

  type t_solverJacobi

    ! System matrix
    type(t_matrixBlock) :: rmatrix

    ! Temporal vector
    type(t_vectorBlock) :: rtempVector
  end type t_solverJacobi

!</typeblock>

  ! *****************************************************************************

!<typeblock>

  ! This data structure contains all local data for a (S)SOR algorithm

  type t_solverSSOR

    ! System matrix
    type(t_matrixBlock) :: rmatrix

    ! Temporal vector
    type(t_vectorBlock) :: rtempVector
  end type t_solverSSOR

!</typeblock>

  ! *****************************************************************************

!<typeblock>

  ! This data structure contains all local data for a BiCGSTAB algorithm

  type t_solverBiCGSTAB

    ! System matrix
    type(t_matrixBlock) :: rmatrix

    ! Temporal vectors
    type(t_vectorBlock), dimension(8) :: RtempVectors

    ! Preconditioner
    type(t_solver), pointer :: p_precond => null()
  end type t_solverBiCGSTAB

!</typeblock>

  ! *****************************************************************************

!<typeblock>

  ! This data structure contains all local data for a GMRES algorithm

  type t_solverGMRES

    ! System matrix
    type(t_matrixBlock) :: rmatrix

    ! Temporal vectors
    type(t_vectorBlock), dimension(:), pointer :: rv
    type(t_vectorBlock), dimension(:), pointer :: rz

    ! Handles for temporal vectors
    integer :: h_Dh = ST_NOHANDLE
    integer :: h_Dc = ST_NOHANDLE
    integer :: h_Ds = ST_NOHANDLE
    integer :: h_Dq = ST_NOHANDLE

    ! Preconditioner
    type(t_solver), pointer :: p_precond => null()

  end type t_solverGMRES

!</typeblock>

  ! *****************************************************************************

!<typeblock>

  ! This data structure contains all local data for ILU algorithm

  type t_solverILU

    ! INPUT: For IFILL > 0, it defines the fill-in level for
    !        for the incomplete LU decomposition
    !        For IFILL = 0, the standard ILU decomposition is performed
    !        For IFILL =-1, the standard ILU decomposition is performed with shifting
    !        For IFILL =-2, the modified ILU decomposition is performed with shifting
    integer :: ifill = 0

    ! INPUT: Relaxation factor for (M)ILU(s)
    real(DP) :: domega = 1.0_DP

    ! Handle to array [1..*] of integer.
    ! Workspace containing the decomposed matrix.
    integer :: h_Idata = ST_NOHANDLE

    ! Handle to array [1..NA] of doubles.
    integer :: h_Ddata = ST_NOHANDLE

    ! Parameter 1 returned by SPLIB for factorised matrix
    integer :: lu = 0

    ! Parameter 2 returned by SPLIB for factorised matrix
    integer :: jlu = 0

    ! Parameter 3 returned by SPLIB for factorised matrix
    integer :: ilup = 0

    ! INPUT: tolerance for ILU-factorisation
    real(DP) :: depsILU = 0.0_DP

    ! INPUT: ILU decomposition of the system matrix
    type(t_matrixBlock) :: rmatrix

    ! Temporal vector
    type(t_vectorBlock) :: rtempVector

    ! Structure for block ILU
    type(t_solverILU), dimension(:), pointer :: p_rsolverBlockILU => null()
  end type t_solverILU

!</typeblock>

  ! *****************************************************************************

!<typeblock>

  ! This data structure contains all local data for a defect correction algorithm

  type t_solverDefcor

    ! Temporal vectors
    type(t_vectorBlock), dimension(3) :: RtempVectors
  end type t_solverDefcor

!</typeblock>

  ! *****************************************************************************

!<typeblock>

  ! This data structure contains all local data for Newton`s algorithm

  type t_solverNewton

    ! INPUT: check sufficient decrease condition?
    integer :: icheckSufficientDecrease = 0

    ! INPUT: maximum number of backtracking steps
    integer :: nmaxBacktrackingSteps = 0

    ! INPUT: number of steps before Jacobian matrix needs update
    integer :: iupdateFrequency = 1

    ! INPUT: strategy for chosing the forcing term
    real(DP) :: dforcingStrategy = 0.0_DP

    ! INPUT: strategy for computing the perturbation parameter
    real(DP) :: dperturbationStrategy = 0.0_DP

    ! Temporal vectors
    type(t_vectorBlock), dimension(5) :: RtempVectors
  end type t_solverNewton

!</typeblock>

!</types>

contains

  ! *****************************************************************************

!<subroutine>

  recursive subroutine solver_createSolverDirect(rparlist, ssectionName, rsolver)

!<description>
    ! This subroutine creates a new solver structure from a given
    ! parameter list. Note that this routine is called recursively
    ! in order  to create the entire solver structure top-down.
    !
    ! Remark, this subroutine should only be called for the top-level
    ! solver. All solver subnodes are allocated and created recursively.
!</description>

!<input>
    ! Parameter list containing all data
    type(t_parlist), intent(in) :: rparlist

    ! Section name of the parameter list containing solver data
    character(LEN=*), intent(in) :: ssectionName
!</input>

!<output>
    ! Solver structure
    type(t_solver), intent(out) :: rsolver
!</output>
!</subroutine>

    ! local variables
    character(LEN=SYS_STRLEN) :: ssolverName,sprecondName,ssmootherName
    integer, dimension(2) :: Isize
    integer :: csolverType,nKrylov,isolver,ismoother,nsolver

    ! Get solver type from parameter list
    call parlst_getvalue_int(rparlist, ssectionName, "csolvertype", csolverType)

    ! The INTENT(out) already initializes rsolver with the most
    ! important information. In the first part mandatory and optional
    ! parameters which apply to some type of solvers are read from
    ! the parameter list and applied to the solver.
    rsolver%ssolverName = trim(adjustl(ssectionName))
    rsolver%csolverType = csolverType

    ! What kind of solver are we
    select case(csolverType)

    case (SV_COUPLED)
      ! ----------------------------------------------------------------------------------
      ! Create a coupled solver:
      ! ~~~~~~~~~~~~~~~~~~~~~~~~
      ! This type of solver has an array of p_solverSubnode subnodes, where each
      ! solver subnode represents the complete solver structure for an individual
      ! component of the coupled problem.
      nsolver = parlst_querysubstrings(rparlist, ssectionName, "ssolvername")
      allocate(rsolver%p_solverSubnode(nsolver))
      do isolver = 1, nsolver
        call parlst_getvalue_string(rparlist, ssectionName, "ssolvername",&
                                    ssolverName, isubstring=isolver)
        call solver_createSolver(rparlist, ssolverName,&
                                 rsolver%p_solverSubnode(isolver))
      end do

    case (SV_FMG)
      ! ----------------------------------------------------------------------------------
      ! Create a full multigrid solver:
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! This type of solver has a p_solverMultigrid subnode and a p_solverCoarsegrid
      ! subnode which is located below the multigrid solver subnode. Both structures
      ! are allocated and filled with required values.
      allocate(rsolver%p_solverMultigrid)
      call parlst_getvalue_int(rparlist, ssectionName, "nlmin",&
                               rsolver%p_solverMultigrid%nlmin)
      call parlst_getvalue_int(rparlist, ssectionName, "nlmax",&
                               rsolver%p_solverMultigrid%nlmax)
      call parlst_getvalue_int(rparlist, ssectionName, "icycle", &
                               rsolver%p_solverMultigrid%icycle)
      call parlst_getvalue_int(rparlist, ssectionName, "ilmin", &
                               rsolver%p_solverMultigrid%ilmin)
      call parlst_getvalue_int(rparlist, ssectionName, "ilmax", &
                               rsolver%p_solverMultigrid%ilmax)

      ! Set initial levels
      rsolver%p_solverMultigrid%initialNlmin = rsolver%p_solverMultigrid%nlmin
      rsolver%p_solverMultigrid%initialNlmax = rsolver%p_solverMultigrid%nlmax

      ! Coarsegrid solver:
      ! ~~~~~~~~~~~~~~~~~~
      ! Ok, we need a so-called coarse grid solver which is used to solve the
      ! nonlinear algebraic system on all multigrid levels from NLMIN up to NLMAX.
      ! Note that it suffices to have one nonlinear coarse grid solver for all
      ! multigrid levels since the nonlinear solver has no matrices and all
      ! temporal vectors are resised in the nonlinear solver routines.
      allocate(rsolver%p_solverMultigrid%p_solverCoarsegrid)

      ! Get type of solver and create solver structure recursively:
      ! For full multigrid only nonlinear multi- or single-grid solvers are
      ! admissible as coarse grid solvers. All other solver type must be rejected.
      call parlst_getvalue_int(rparlist, ssectionName, "isolver", isolver)
      select case(isolver)

      case (SV_NONLINEARMG, SV_NONLINEAR)
        call parlst_getvalue_string(rparlist, ssectionName, "ssolvername", ssolverName)
        call solver_createSolver(rparlist, ssolverName,&
                                 rsolver%p_solverMultigrid%p_solverCoarsegrid)

      case DEFAULT
        call output_line('Invalid nonlinear solver for full multigrid solver!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'solver_createSolverDirect')
        call sys_halt()
      end select
      ! Now, we are done with the full multigrid solver !!!


    case (SV_NONLINEARMG)
      ! ----------------------------------------------------------------------------------
      ! Create a nonlinear multigrid solver:
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! This type of solver has a p_solverMultigrid subnode and a p_solverCoarsegrid
      ! subnode which is located below the multigrid solver subnode. Both structures
      ! are allocated and filled with required values.
      ! In addition, a nonlinear multigrid solver may have a nonlinear smoother which
      ! is used to improve the approximation of the nonlinear solution on each level.
      allocate(rsolver%p_solverMultigrid)
      call parlst_getvalue_int(rparlist, ssectionName, "nlmin", &
                               rsolver%p_solverMultigrid%nlmin)
      call parlst_getvalue_int(rparlist, ssectionName, "nlmax", &
                               rsolver%p_solverMultigrid%nlmax)
      call parlst_getvalue_int(rparlist, ssectionName, "icycle", &
                               rsolver%p_solverMultigrid%icycle)
      call parlst_getvalue_int(rparlist, ssectionName, "ilmin", &
                               rsolver%p_solverMultigrid%ilmin)
      call parlst_getvalue_int(rparlist, ssectionName, "ilmax", &
                               rsolver%p_solverMultigrid%ilmax)

      ! Set initial levels
      rsolver%p_solverMultigrid%initialNlmin = rsolver%p_solverMultigrid%nlmin
      rsolver%p_solverMultigrid%initialNlmax = rsolver%p_solverMultigrid%nlmax

      ! Coarsegrid solver:
      ! ~~~~~~~~~~~~~~~~~~
      ! Ok, we need a so-called coarse grid solver which is used to solve the
      ! nonlinear algebraic system on the coarsest grid.
      allocate(rsolver%p_solverMultigrid%p_solverCoarsegrid)

      ! Get type of solver and create solver structure recursively:
      ! For nonlinear multigrid only nonlinear multi- or single-grid solvers are
      ! admissible as coarse grid solvers. All other solver types must be rejected.
      call parlst_getvalue_int(rparlist, ssectionName, "isolver", isolver)
      select case(isolver)

      case (SV_FMG, SV_NONLINEARMG, SV_NONLINEAR)
        call parlst_getvalue_string(rparlist, ssectionName, "ssolvername", ssolverName)
        call solver_createSolver(rparlist, ssolverName,&
                                 rsolver%p_solverMultigrid%p_solverCoarsegrid)

      case DEFAULT
        call output_line('Invalid coarse grid solver for nonlinear multigrid solver!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'solver_createSolverDirect')
        call sys_halt()
      end select

      ! Smoother:
      ! ~~~~~~~~~
      ! Optionally, we need a pre-/postsmoothing algorithm to smooth
      ! the nonlinear algebraic equation. Check if either pre- or
      ! postsmoothing should be applied and allocate smoother
      ! structute accordingly.
      call parlst_getvalue_int(rparlist, ssectionName, "npresmooth", &
                               rsolver%p_solverMultigrid%npresmooth)
      call parlst_getvalue_int(rparlist, ssectionName, "npostsmooth", &
                               rsolver%p_solverMultigrid%npostsmooth)

      ! Create smoother recursively
      if ((rsolver%p_solverMultigrid%npresmooth  .ne. 0) .or. &
          (rsolver%p_solverMultigrid%npostsmooth .ne. 0)) then
        call parlst_getvalue_int(rparlist, ssectionName, "nsmoothFactor", &
                                 rsolver%p_solverMultigrid%nsmoothFactor)
        call parlst_getvalue_int(rparlist, ssectionName, "ismoother", &
                                 ismoother)
        call parlst_getvalue_string(rparlist, ssectionName, "ssmootherName", &
                                    ssmootherName)

        ! For nonlinear multigrid only nonlinear multigrid- or single
        ! -grid solvers are admissible as smoothers. All other solver
        ! types must be rejected.
        select case(ismoother)

        case (SV_NONLINEAR, SV_NONLINEARMG)
          rsolver%p_solverMultigrid%csmootherType = ismoother
          call create_smoother(rparlist, ssmootherName, rsolver%p_solverMultigrid)

        case DEFAULT
          call output_line('Invalid smoother for nonlinear multigrid solver!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'solver_createSolverDirect')
          call sys_halt()
        end select
      end if
      ! Now, we are done with the nonlinear multigrid solver !!!


    case (SV_LINEARMG)
      ! ----------------------------------------------------------------------------------
      ! Create a linear multigrid solver:
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! This type of solver has a p_solverMultigrid subnode and a p_solverCoarsegrid
      ! subnode which is located below the multigrid solver subnode. Both structures
      ! are allocated and filled with required values.
      ! In addition, a linear multigrid solver may have a linear smoother which is used
      ! to improve the approximation of the linear solution on each level.
      allocate(rsolver%p_solverMultigrid)
      call parlst_getvalue_int(rparlist, ssectionName, "nlmin", &
                               rsolver%p_solverMultigrid%nlmin)
      call parlst_getvalue_int(rparlist, ssectionName, "nlmax", &
                               rsolver%p_solverMultigrid%nlmax)
      call parlst_getvalue_int(rparlist, ssectionName, "icycle", &
                               rsolver%p_solverMultigrid%icycle)
      call parlst_getvalue_int(rparlist, ssectionName, "ilmin", &
                               rsolver%p_solverMultigrid%ilmin)
      call parlst_getvalue_int(rparlist, ssectionName, "ilmax", &
                               rsolver%p_solverMultigrid%ilmax)

      ! Set initial levels
      rsolver%p_solverMultigrid%initialNlmin = rsolver%p_solverMultigrid%nlmin
      rsolver%p_solverMultigrid%initialNlmax = rsolver%p_solverMultigrid%nlmax

      ! Coarsegrid solver:
      ! ~~~~~~~~~~~~~~~~~~
      ! Ok, we need a so-called coarse grid solver which is used to solve the
      ! linear algebraic equation on the coarsest grid.
      allocate(rsolver%p_solverMultigrid%p_solverCoarsegrid)

      ! Get type of solver and create solver structure recursively:
      ! For linear multigrid only linear multi- or single-grid solvers are
      ! admissible as coarse grid solvers. All other solver types must be rejected.
      call parlst_getvalue_int(rparlist, ssectionName, "isolver", isolver)
      select case(isolver)

      case (SV_LINEARMG, SV_LINEAR)
        call parlst_getvalue_string(rparlist, ssectionName, "ssolverName", ssolverName)
        call solver_createSolver(rparlist, ssolverName,&
                                 rsolver%p_solverMultigrid%p_solverCoarsegrid)

      case DEFAULT
        call output_line('Invalid coarse grid solver for linear multigrid solver!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'solver_createSolverDirect')
        call sys_halt()
      end select

      ! Smoother:
      ! ~~~~~~~~~
      ! Optionally, we need a pre-/postsmoothing algorithm to smooth the linear algebraic equation.
      ! Check if either pre- or postsmoothing should be applied and allocate smoother accordingly.
      call parlst_getvalue_int(rparlist, ssectionName, "npresmooth", &
                               rsolver%p_solverMultigrid%npresmooth)
      call parlst_getvalue_int(rparlist, ssectionName, "npostsmooth", &
                               rsolver%p_solverMultigrid%npostsmooth)

      ! Create smoother recursively
      if ((rsolver%p_solverMultigrid%npresmooth  .ne. 0) .or. &
          (rsolver%p_solverMultigrid%npostsmooth .ne. 0)) then
        call parlst_getvalue_int(rparlist, ssectionName, "nsmoothfactor", &
                                 rsolver%p_solverMultigrid%nsmoothfactor)
        call parlst_getvalue_int(rparlist, ssectionName, "ismoother", &
                                 ismoother)
        call parlst_getvalue_string(rparlist, ssectionName, "ssmootherName", &
                                    ssmootherName)

        ! For linear multigrid only linear multigrid- or single-grid solvers are
        ! admissible as smoothers. All other solver types must be rejected
        select case(ismoother)

        case (SV_LINEAR, SV_LINEARMG)
          call create_smoother(rparlist, ssmootherName, rsolver%p_solverMultigrid)

        case DEFAULT
          call output_line('Invalid smoother for nonlinear multigrid solver!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'solver_createSolverDirect')
          call sys_halt()
        end select
      end if
      ! Now, we are done with the linear multigrid solver !!!


    case (SV_NONLINEAR)
      ! ----------------------------------------------------------------------------------
      ! Create a nonlinear single-grid solver:
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! This type of solver has no explicit solver subnode but a nonlinear
      ! preconditioner structure. According to the type of nonlinear preconditioning
      ! the corresponding structure is allocated and filled with values
      call parlst_getvalue_int(rparlist, ssectionName, "isolver", &
                               rsolver%isolver)

      ! Ok, now we set up the individual solver structure
      select case (rsolver%isolver)
      case (NLSOL_SOLVER_FIXEDPOINT)
        ! Ok, set up the nonlinear preconditioner
        call create_preconditioner(rparlist, ssectionName, rsolver)

      case DEFAULT
        call output_line('Unsupported nonlinear solver!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'solver_createSolverDirect')
        call sys_halt()
      end select

      ! Finally, we need a subnode for the linear solver. Here, we do not know
      ! if a linear multigrid or single grid solver should be employed. Hence,
      ! we just get the name of the linear solver from the parameter list and
      ! create the linear solver recursively:
      allocate(rsolver%p_solverSubnode(1))
      call parlst_getvalue_string(rparlist, ssectionName, "ssolverName", &
                                  ssolverName)
      call solver_createSolver(rparlist, ssolverName, rsolver%p_solverSubnode(1))
      ! Now, we are done with the nonlinear singlegrid solver !!!

    case (SV_LINEAR)
      ! ----------------------------------------------------------------------------------
      ! Create a linear singlegrid solver:
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! This type of solver has several subnodes for the individual single grid solver.
      ! According to the desired solver, the corresponding subnode is allocated and
      ! filled with values.
      call parlst_getvalue_int(rparlist, ssectionName, "isolver", &
                               rsolver%isolver)

      ! Ok, now set up the individual solver structures
      select case (rsolver%isolver)
      case (LINSOL_SOLVER_UMFPACK4)
        allocate(rsolver%p_solverUMFPACK)
        ! Initialize control structure
        call UMF4DEF(rsolver%p_solverUMFPACK%Dcontrol)

      case (LINSOL_SOLVER_JACOBI)
        allocate(rsolver%p_solverJacobi)

      case (LINSOL_SOLVER_SOR,LINSOL_SOLVER_SSOR)
        allocate(rsolver%p_solverSSOR)

      case (LINSOL_SOLVER_BICGSTAB)
        allocate(rsolver%p_solverBiCGSTAB)

        ! Check for preconditioner
        call parlst_getvalue_string(rparlist, ssectionName, "sprecondName", &
                                    sprecondName)

        if (trim(adjustl(sprecondName)) .ne. "") then
          allocate(rsolver%p_solverBiCGSTAB%p_precond)
          call create_preconditioner(rparlist, sprecondName,&
                                     rsolver%p_solverBiCGSTAB%p_precond)
        end if

      case (LINSOL_SOLVER_GMRES)
        allocate(rsolver%p_solverGMRES)

        ! Get dimension of Krylov subspace
        call parlst_getvalue_int(rparlist, ssectionName, "nkrylov", &
                                 nKrylov)
        allocate(rsolver%p_solverGMRES%rv(nKrylov+1))
        allocate(rsolver%p_solverGMRES%rz(nKrylov))

        Isize=nKrylov+1
        call storage_new('solver_createDirect', 'Dh', Isize, ST_DOUBLE,&
                         rsolver%p_solverGMRES%h_Dh, ST_NEWBLOCK_NOINIT)
        call storage_new('solver_createDirect', 'Dc', nKrylov, ST_DOUBLE,&
                         rsolver%p_solverGMRES%h_Dc, ST_NEWBLOCK_NOINIT)
        call storage_new('solver_createDirect', 'Ds', nKrylov, ST_DOUBLE,&
                         rsolver%p_solverGMRES%h_Ds, ST_NEWBLOCK_NOINIT)
        call storage_new('solver_createDirect', 'Dq', nKrylov+1, ST_DOUBLE,&
                         rsolver%p_solverGMRES%h_Dq, ST_NEWBLOCK_NOINIT)

        ! Check for preconditioner
        call parlst_getvalue_string(rparlist, ssectionName, "sprecondName", &
                                    sprecondName)

        if (trim(adjustl(sprecondName)) .ne. "") then
          allocate(rsolver%p_solverGMRES%p_precond)
          call create_preconditioner(rparlist, sprecondName,&
                                     rsolver%p_solverGMRES%p_precond)
        end if

      case (LINSOL_SOLVER_ILU)
        allocate(rsolver%p_solverILU)
        call parlst_getvalue_double(rparlist, ssectionName, "depsILU", &
                                    rsolver%p_solverILU%depsILU)
        call parlst_getvalue_double(rparlist, ssectionName, "domega", &
                                    rsolver%p_solverILU%domega)
        call parlst_getvalue_int(rparlist, ssectionName, "ifill", &
                                 rsolver%p_solverILU%ifill)

      case DEFAULT
        call output_line('Unsupported linear solver!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'solver_createSolverDirect')
        call sys_halt()
      end select

    case DEFAULT
      call output_line('Unsupported solver type!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'solver_createSolverDirect')
      call sys_halt()
    end select

    ! In the second part, parameters which apply to all solvers are read
    call parlst_getvalue_int(rparlist, ssectionName, "ioutputlevel", &
                             rsolver%ioutputlevel)
    call parlst_getvalue_int(rparlist, ssectionName, "iresNorm", &
                             rsolver%iresNorm)
    call parlst_getvalue_int(rparlist, ssectionName, "istoppingCriterion", &
                             rsolver%istoppingCriterion)
    call parlst_getvalue_int(rparlist, ssectionName, "nminIterations", &
                             rsolver%nminIterations)
    call parlst_getvalue_int(rparlist, ssectionName, "nmaxIterations", &
                             rsolver%nmaxIterations)
    call parlst_getvalue_double(rparlist, ssectionName, "domega", &
                                rsolver%domega)
    call parlst_getvalue_double(rparlist, ssectionName, "depsRel", &
                                rsolver%depsRel)
    call parlst_getvalue_double(rparlist, ssectionName, "depsAbs", &
                                rsolver%depsAbs)
    call parlst_getvalue_double(rparlist, ssectionName, "depsStag", &
                                rsolver%depsStag)
    call parlst_getvalue_double(rparlist, ssectionName, "ddivRel", &
                                rsolver%ddivRel)
    call parlst_getvalue_double(rparlist, ssectionName, "ddivAbs", &
                                rsolver%ddivAbs)
    call parlst_getvalue_double(rparlist, ssectionName, "drhsZero", &
                                rsolver%drhsZero)
    call parlst_getvalue_double(rparlist, ssectionName, "ddefZero", &
                                rsolver%ddefZero)

  contains

    ! Here, the real working routine follows

    !*************************************************************
    ! This subroutine turns an structure of type t_solver into a
    ! preconditioner. Depending on the type of solver, various
    ! parameters are modified so as to meet the requirements of
    ! a preconditioner of read in from parameter list

    subroutine create_preconditioner(rparlist, ssectionName, rsolver)
      ! Parameter list containing all data
      type(t_parlist), intent(in) :: rparlist

      ! Section name of the parameter list containing solver data
      character(LEN=*), intent(in) :: ssectionName

      ! Solver structure
      type(t_solver), intent(inout) :: rsolver

      ! local variables
      integer :: csolverType


      ! Get preconditioner type from parameter list
      call parlst_getvalue_int(rparlist, ssectionName, "csolverType", csolverType)

      ! Check if given solver matches the current solver type
      if ((rsolver%csolverType .ne. SV_UNDEFINED) .and.&
          (rsolver%csolverType .ne. csolverType)) then
        call output_line('The given solver does not match specified type!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'create_preconditioner')
        call sys_halt()
      end if

      ! In the first part mandatory and optional parameters which apply to some type
      ! of solvers are read from the parameter list and applied to the solver.
      rsolver%ssolverName = trim(adjustl(ssectionName))
      rsolver%csolverType = csolverType

      ! What kind of solver are we?
      select case(csolverType)

      case (SV_NONLINEARMG)
        ! ----------------------------------------------------------------------------------
        ! Create a nonlinear multigrid preconditioner:
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! We have to create a preconditioner as a real solver subnode
        call solver_createSolver(rparlist, ssectionName, rsolver)

      case (SV_NONLINEAR)
        ! ----------------------------------------------------------------------------------
        ! Create a nonlinear single-grid preconditioner:
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        call parlst_getvalue_int(rparlist, ssectionName, "iprecond", &
                                 rsolver%iprecond)

        ! What kind of preconditioner should be created?
        select case(rsolver%iprecond)
        case (NLSOL_PRECOND_BLOCKD,&
              NLSOL_PRECOND_DEFCOR)
          allocate(rsolver%p_solverDefcor)

        case (NLSOL_PRECOND_NEWTON)
          allocate(rsolver%p_solverNewton)
          call parlst_getvalue_int(rparlist, ssectionName, "icheckSufficientDecrease",&
                                   rsolver%p_solverNewton%icheckSufficientDecrease)
          call parlst_getvalue_int(rparlist, ssectionName, "nmaxBacktrackingSteps",&
                                   rsolver%p_solverNewton%nmaxBacktrackingSteps)
          call parlst_getvalue_int(rparlist, ssectionName, "iupdateFrequency",&
                                   rsolver%p_solverNewton%iupdateFrequency)
          call parlst_getvalue_double(rparlist, ssectionName, "dforcingStrategy",&
                                      rsolver%p_solverNewton%dforcingStrategy)
          call parlst_getvalue_double(rparlist, ssectionName, "dperturbationStrategy",&
                                      rsolver%p_solverNewton%dperturbationStrategy)

        case DEFAULT
          call output_line('Unsupported nonlinear preconditioner!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'create_preconditioner')
          call sys_halt()
        end select

      case (SV_LINEARMG)
        ! ----------------------------------------------------------------------------------
        ! Create a linear multigrid preconditioner:
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! We have to create a preconditioner as a real solver subnode
        call solver_createSolver(rparlist, ssectionName, rsolver)

      case (SV_LINEAR)
        ! ----------------------------------------------------------------------------------
        ! Create a linear preconditioner:
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        call parlst_getvalue_int(rparlist, ssectionName, "isolver", &
                                 rsolver%isolver)

        ! What kind of preconditioner should be created?
        select case(rsolver%isolver)
        case (LINSOL_SOLVER_NONE)
          ! do nothing

        case (LINSOL_SOLVER_JACOBI)
          allocate(rsolver%p_solverJacobi)

        case (LINSOL_SOLVER_SOR,LINSOL_SOLVER_SSOR)
          allocate(rsolver%p_solverSSOR)

        case (LINSOL_SOLVER_ILU)
          allocate(rsolver%p_solverILU)
          call parlst_getvalue_double(rparlist, ssectionName, "depsILU", &
                                      rsolver%p_solverILU%depsILU)
          call parlst_getvalue_double(rparlist, ssectionName, "domega", &
                                      rsolver%p_solverILU%domega)
          call parlst_getvalue_int(rparlist, ssectionName, "ifill", &
                                   rsolver%p_solverILU%ifill)

        case DEFAULT
          ! Ok, we could not identify one of the simple preconditioners.
          ! Let us try to create the preconditioner as real solver subnode
          call solver_createSolver(rparlist, ssectionName, rsolver)
        end select

        ! Let us (re-)set some data for the linear preconditioner
        rsolver%nminIterations   = 1
        rsolver%nmaxIterations   = 1
        rsolver%domega           = 1.0_DP
        rsolver%drhsZero         = 0.0_DP
        rsolver%ddefZero         = 0.0_DP
        rsolver%depsAbs          = 0.0_DP
        rsolver%depsRel          = 0.0_DP
        rsolver%depsStag         = 0.0_DP
        rsolver%ddivRel          = SYS_INFINITY
        rsolver%ddivAbs          = SYS_INFINITY

      case DEFAULT
        call output_line('Unsupported preconditioner type!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'create_preconditioner')
        call sys_halt()
      end select

    end subroutine create_preconditioner

    !*************************************************************
    ! This subroutine creates the smoother subnode for a given
    ! structure of type t_solverMultigrid. Each smoother subnode
    ! is created individually from the parameter list.

    subroutine create_smoother(rparlist, ssectionName, rsolver)
      ! Parameter list containing all data
      type(t_parlist), intent(in) :: rparlist

      ! Section name of the parameter list containing smoother data
      character(LEN=*), intent(in) :: ssectionName

      ! Multigrid solver structure
      type(t_solverMultigrid), intent(inout) :: rsolver

      ! local variables
      integer :: i,csolverType

      ! Get smoother type from parameter list
      call parlst_getvalue_int(rparlist, ssectionName, "csolverType", &
                               csolverType)

      ! In the first part mandatory and optional parameters which apply to some type
      ! of solvers are read from the parameter list and applied to the solverrr.
      rsolver%ssmootherName = trim(adjustl(ssectionName))
      rsolver%csmootherType = csolverType

      ! What kind of smoother are we?
      select case(csolverType)
      case (SV_NONLINEAR)
        ! ----------------------------------------------------------------------------------
        ! Create a nonlinear smoother:
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! Create nonlinear smoother recursively. Note that we create only one (!) template
        ! smoother which will be used (later) do generate smoother subnodes on all levels.
        if (rsolver%nlmin .lt. rsolver%nlmax) then
          allocate(rsolver%p_smoother(rsolver%nlmin+1:rsolver%nlmax))
          do i = rsolver%nlmin+1, rsolver%nlmax
            call solver_createSolver(rparlist, ssectionName, rsolver%p_smoother(i))
          end do
        end if

      case (SV_LINEAR)
        ! ----------------------------------------------------------------------------------
        ! Create a linear smoother:
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~

        ! Create linear smoother recursively. Note that we create an array of solvers
        ! for each level. Each level is created and initialised individually.
        if (rsolver%nlmin .lt. rsolver%nlmax) then
          allocate(rsolver%p_smoother(rsolver%nlmin+1:rsolver%nlmax))
          do i = rsolver%nlmin+1, rsolver%nlmax
            call solver_createSolver(rparlist, ssectionName, rsolver%p_smoother(i))
          end do
        end if

      case DEFAULT
        call output_line('Unsupported type of smoother!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'create_smoother')
        call sys_halt()
      end select

    end subroutine create_smoother

  end subroutine solver_createSolverDirect

  ! ***************************************************************************

!<subroutine>

  recursive subroutine solver_createSolverIndirect(rsolver, rsolverTemplate)

!<description>
    ! This subroutine creates a new solver structure by cloning an existing solver structure
!</description>

!<input>
    ! Template solver structure
    type(t_solver), intent(in) :: rsolverTemplate
!</input>

!<output>
    ! Solver structure
    type(t_solver), intent(out) :: rsolver
!</output>
!</subroutine>

    ! local variables
    integer :: i

    ! The INTENT(out) already initializes rsolver with the most
    ! important information. The rest comes now
    rsolver = rsolverTemplate

    ! Mark the solver as "empty" so that it structure and
    ! its content is updated automatically
    rsolver%isolverSpec = ior(rsolver%isolverSpec, SV_SSPEC_EMPTY)

    ! Set boundary pointer
    if (associated(rsolverTemplate%rboundaryCondition))&
        rsolver%rboundaryCondition => rsolverTemplate%rboundaryCondition

    ! Next, we have to proceed to the subnodes and create them recursively

    ! Generic solver subnode
    if (associated(rsolverTemplate%p_solverSubnode)) then
      allocate(rsolver%p_solverSubnode(size(rsolverTemplate%p_solverSubnode)))
      do i = 1, size(rsolverTemplate%p_solverSubnode)
        call solver_createSolver(rsolver%p_solverSubnode(i),&
                                 rsolverTemplate%p_solverSubnode(i))
      end do
    else
      nullify(rsolver%p_solverSubnode)
    end if

    ! Multigrid solver
    if (associated(rsolverTemplate%p_solverMultigrid)) then
      allocate(rsolver%p_solverMultigrid)
      call create_solverMultigrid(rsolver%p_solverMultigrid,&
                                  rsolverTemplate%p_solverMultigrid)
    else
      nullify(rsolver%p_solverMultigrid)
    end if

    ! UMFPACK solver
    if (associated(rsolverTemplate%p_solverUMFPACK)) then
      allocate(rsolver%p_solverUMFPACK)
      call create_solverUMFPACK(rsolver%p_solverUMFPACK,&
                                rsolverTemplate%p_solverUMFPACK)
    else
      nullify(rsolver%p_solverUMFPACK)
    end if

    ! Jacobi solver
    if (associated(rsolverTemplate%p_solverJacobi)) then
      allocate(rsolver%p_solverJacobi)
      call create_solverJacobi(rsolver%p_solverJacobi,&
                               rsolverTemplate%p_solverJacobi)
    else
      nullify(rsolver%p_solverJacobi)
    end if

    ! (S)SOR solver
    if (associated(rsolverTemplate%p_solverSSOR)) then
      allocate(rsolver%p_solverSSOR)
      call create_solverSSOR(rsolver%p_solverSSOR,&
                             rsolverTemplate%p_solverSSOR)
    else
      nullify(rsolver%p_solverSSOR)
    end if

    ! BiCGSTAB solver
    if (associated(rsolverTemplate%p_solverBiCGSTAB)) then
      allocate(rsolver%p_solverBiCGSTAB)
      call create_solverBiCGSTAB(rsolver%p_solverBiCGSTAB,&
                                 rsolverTemplate%p_solverBiCGSTAB)
    else
      nullify(rsolver%p_solverBiCGSTAB)
    end if

    ! GMRES solver
    if (associated(rsolverTemplate%p_solverGMRES)) then
      allocate(rsolver%p_solverGMRES)
      call create_solverGMRES(rsolver%p_solverGMRES,&
                              rsolverTemplate%p_solverGMRES)
    else
      nullify(rsolver%p_solverGMRES)
    end if

    ! ILU solver
    if (associated(rsolverTemplate%p_solverILU)) then
      allocate(rsolver%p_solverILU)
      call create_solverILU(rsolver%p_solverILU,&
                            rsolverTemplate%p_solverILU)
    else
      nullify(rsolver%p_solverILU)
    end if

    ! Defect correction algorithm
    if (associated(rsolverTemplate%p_solverDefcor)) then
      allocate(rsolver%p_solverDefcor)
      call create_solverDefcor(rsolver%p_solverDefcor,&
                               rsolverTemplate%p_solverDefcor)
    else
      nullify(rsolver%p_solverDefcor)
    end if

    ! Newton`s algorithm
    if (associated(rsolverTemplate%p_solverNewton)) then
      allocate(rsolver%p_solverNewton)
      call create_solverNewton(rsolver%p_solverNewton,&
                               rsolverTemplate%p_solverNewton)
    else
      nullify(rsolver%p_solverNewton)
    end if

  contains

    ! Here, the real working routine follows

    !*************************************************************
    ! Create a multigrid solver

    subroutine create_solverMultigrid(rsolver, rsolverTemplate)
      type(t_solverMultigrid), intent(in) :: rsolverTemplate
      type(t_solverMultigrid), intent(out) :: rsolver

      ! local variables
      integer :: i,ilbound,iubound

      ! The INTENT(out) already initializes rsolver with the most important
      ! information. The rest comes now
      rsolver = rsolverTemplate

      ! Create coarsegrid solver if required
      if (associated(rsolverTemplate%p_solverCoarsegrid)) then
        allocate(rsolver%p_solverCoarsegrid)
        call solver_createSolver(rsolver%p_solverCoarsegrid,&
                                 rsolverTemplate%p_solverCoarsegrid)
      else
        nullify(rsolver%p_solverCoarsegrid)
      end if

      ! Create smoother if required
      if (associated(rsolverTemplate%p_smoother)) then
        ilbound = lbound(rsolverTemplate%p_smoother,1)
        iubound = ubound(rsolverTemplate%p_smoother,1)
        allocate(rsolver%p_smoother(ilbound:iubound))
        do i = ilbound, iubound
          call solver_createSolver(rsolver%p_smoother(i),&
                                   rsolverTemplate%p_smoother(i))
        end do
      else
        nullify(rsolver%p_smoother)
      end if

      ! Nullify pointers
      nullify(rsolver%rmatrix)
      nullify(rsolver%RtempVectors)
    end subroutine create_solverMultigrid

    !*************************************************************
    ! Create an UMFPACK solver

    subroutine create_solverUMFPACK(rsolver, rsolverTemplate)
      type(t_solverUMFPACK), intent(in) :: rsolverTemplate
      type(t_solverUMFPACK), intent(out) :: rsolver

      ! Do nothing, the INTENT(out) initializes the solver

    end subroutine create_solverUMFPACK

    !*************************************************************
    ! Create a Jacobi solver

    subroutine create_solverJacobi(rsolver, rsolverTemplate)
      type(t_solverJacobi), intent(in) :: rsolverTemplate
      type(t_solverJacobi), intent(out) :: rsolver

      ! Do nothing, the INTENT(out) initializes the solver

    end subroutine create_solverJacobi

    !*************************************************************
    ! Create an (S)SOR solver

    subroutine create_solverSSOR(rsolver, rsolverTemplate)
      type(t_solverSSOR), intent(in) :: rsolverTemplate
      type(t_solverSSOR), intent(out) :: rsolver

      ! Do nothing, the INTENT(out) initializes the solver

    end subroutine create_solverSSOR

    !*************************************************************
    ! Create a BiCGSTAB solver

    subroutine create_solverBiCGSTAB(rsolver, rsolverTemplate)
      type(t_solverBiCGSTAB), intent(in) :: rsolverTemplate
      type(t_solverBiCGSTAB), intent(out) :: rsolver

      ! Create preconditioner if required
      if (associated(rsolverTemplate%p_precond)) then
        allocate(rsolver%p_precond)
        call solver_createSolver(rsolver%p_precond,&
                                 rsolverTemplate%p_precond)
      else
        nullify(rsolver%p_precond)
      end if
    end subroutine create_solverBiCGSTAB

    !*************************************************************
    ! Create a GMRES solver

    subroutine create_solverGMRES(rsolver, rsolverTemplate)
      type(t_solverGMRES), intent(in) :: rsolverTemplate
      type(t_solverGMRES), intent(out) :: rsolver

      ! local variables
      integer :: ilbound,iubound

      ! Create subarrays for temporal vectors
      if (associated(rsolverTemplate%rv)) then
        ilbound = lbound(rsolverTemplate%rv,1)
        iubound = ubound(rsolverTemplate%rv,1)
        allocate(rsolver%rv(ilbound:iubound))
      else
        nullify(rsolver%rv)
      end if

      if (associated(rsolverTemplate%rz)) then
        ilbound = lbound(rsolverTemplate%rz,1)
        iubound = ubound(rsolverTemplate%rz,1)
        allocate(rsolver%rz(ilbound:iubound))
      else
        nullify(rsolver%rz)
      end if

      ! Create handles. Well, we actually copy the handles.
      if (rsolverTemplate%h_Dh .ne. ST_NOHANDLE)&
          call storage_copy(rsolverTemplate%h_Dh, rsolver%h_Dh)
      if (rsolverTemplate%h_Dc .ne. ST_NOHANDLE)&
          call storage_copy(rsolverTemplate%h_Dc, rsolver%h_Dc)
      if (rsolverTemplate%h_Ds .ne. ST_NOHANDLE)&
          call storage_copy(rsolverTemplate%h_Ds, rsolver%h_Ds)
      if (rsolverTemplate%h_Dq .ne. ST_NOHANDLE)&
          call storage_copy(rsolverTemplate%h_Dq, rsolver%h_Dq)

      ! Create preconditioner if required
      if (associated(rsolverTemplate%p_precond)) then
        allocate(rsolver%p_precond)
        call solver_createSolver(rsolver%p_precond,&
                                 rsolverTemplate%p_precond)
      else
        nullify(rsolver%p_precond)
      end if
    end subroutine create_solverGMRES

    !*************************************************************
    ! Create an ILU solver

    subroutine create_solverILU(rsolver, rsolverTemplate)
      type(t_solverILU), intent(in) :: rsolverTemplate
      type(t_solverILU), intent(out) :: rsolver

      ! The INTENT(out) initializes the solver already.
      ! The rest is done now.
      rsolver%ifill   = rsolverTemplate%ifill
      rsolver%domega  = rsolverTemplate%domega
      rsolver%depsILU = rsolverTemplate%depsILU
    end subroutine create_solverILU

    !*************************************************************
    ! Create a defect correction algorithm

    subroutine create_solverDefcor(rsolver, rsolverTemplate)
      type(t_solverDefcor), intent(in) :: rsolverTemplate
      type(t_solverDefcor), intent(out) :: rsolver

      ! Do nothing, the INTENT(out) initializes the solver

    end subroutine create_solverDefcor

    !*************************************************************
    ! Create a Newton algorithm

    subroutine create_solverNewton(rsolver, rsolverTemplate)
      type(t_solverNewton), intent(in) :: rsolverTemplate
      type(t_solverNewton), intent(out) :: rsolver

      ! The INTENT(out) initializes the solver already. The rest is done now.
      rsolver%icheckSufficientDecrease = rsolverTemplate%icheckSufficientDecrease
      rsolver%nmaxBacktrackingSteps    = rsolverTemplate%nmaxBacktrackingSteps
      rsolver%iupdateFrequency         = rsolverTemplate%iupdateFrequency
      rsolver%dforcingStrategy         = rsolverTemplate%dforcingStrategy
      rsolver%dperturbationStrategy    = rsolver%dperturbationStrategy
    end subroutine create_solverNewton

  end subroutine solver_createSolverIndirect

  ! ***************************************************************************

!<subroutine>

  recursive subroutine solver_releaseSolver(rsolver)

!<description>
    ! This subroutine releases a given solver structure,
    ! whereby all of its subnodes are released recursively.
!</description>

!<inputoutput>
    ! Solver structure
    type(t_solver), intent(inout) :: rsolver
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: i

    ! Nullify boundary condition pointer
    nullify(rsolver%rboundaryCondition)

    ! Proceed to the subnodes. Note that the procedure is the same for all
    ! subnodes. First, the subnode itself is release by calling the corresponding
    ! subroutine. If the solver hierarchy does not stop at the given subnode, then
    ! the release process is continued recursively.
    ! Afterwards the subnode is  physically deallocated and its pointer is nullifies.

    ! Generic solver subnode
    if (associated(rsolver%p_solverSubnode)) then
      do i = 1, size(rsolver%p_solverSubnode)
        call solver_releaseSolver(rsolver%p_solverSubnode(i))
      end do
      deallocate(rsolver%p_solverSubnode)
      nullify(rsolver%p_solverSubnode)
    end if

    ! Multigrid solver
    if (associated(rsolver%p_solverMultigrid)) then
      call release_solverMultigrid(rsolver%p_solverMultigrid)
      deallocate(rsolver%p_solverMultigrid)
      nullify(rsolver%p_solverMultigrid)
    end if

    ! UMFPACK solver
    if (associated(rsolver%p_solverUMFPACK)) then
      call release_solverUMFPACK(rsolver%p_solverUMFPACK)
      deallocate(rsolver%p_solverUMFPACK)
      nullify(rsolver%p_solverUMFPACK)
    end if

    ! Jacobi solver
    if (associated(rsolver%p_solverJacobi)) then
      call release_solverJacobi(rsolver%p_solverJacobi)
      deallocate(rsolver%p_solverJacobi)
      nullify(rsolver%p_solverJacobi)
    end if

    ! (S)SOR solver
    if (associated(rsolver%p_solverSSOR)) then
      call release_solverSSOR(rsolver%p_solverSSOR)
      deallocate(rsolver%p_solverSSOR)
      nullify(rsolver%p_solverSSOR)
    end if

    ! BiCGSTAB solver
    if (associated(rsolver%p_solverBiCGSTAB)) then
      call release_solverBiCGSTAB(rsolver%p_solverBiCGSTAB)
      deallocate(rsolver%p_solverBiCGSTAB)
      nullify(rsolver%p_solverBiCGSTAB)
    end if

    ! GMRES solver
    if (associated(rsolver%p_solverGMRES)) then
      call release_solverGMRES(rsolver%p_solverGMRES)
      deallocate(rsolver%p_solverGMRES)
      nullify(rsolver%p_solverGMRES)
    end if

    ! ILU solver
    if (associated(rsolver%p_solverILU)) then
      call release_solverILU(rsolver%p_solverILU)
      deallocate(rsolver%p_solverILU)
      nullify(rsolver%p_solverILU)
    end if

    ! Defect correction algorithm
    if (associated(rsolver%p_solverDefcor)) then
      call release_solverDefcor(rsolver%p_solverDefcor)
      deallocate(rsolver%p_solverDefcor)
      nullify(rsolver%p_solverDefcor)
    end if

    ! Newton`s algorithm
    if (associated(rsolver%p_solverNewton)) then
      call release_solverNewton(rsolver%p_solverNewton)
      deallocate(rsolver%p_solverNewton)
      nullify(rsolver%p_solverNewton)
    end if

  contains

    ! Here, the real working routine follows

    !*************************************************************
    ! Release multigrid solver

    subroutine release_solverMultigrid(rsolver)
      type(t_solverMultigrid), intent(inout) :: rsolver

      ! local variables
      integer :: i

      ! Release matrices
      if (associated(rsolver%rmatrix)) then
        do i = lbound(rsolver%rmatrix,1),&
               ubound(rsolver%rmatrix,1)
          call lsysbl_releaseMatrix(rsolver%rmatrix(i))
        end do
        deallocate(rsolver%rmatrix)
        nullify(rsolver%rmatrix)
      end if

      ! Release temporal vectors
      if (associated(rsolver%RtempVectors)) then
        do i = lbound(rsolver%RtempVectors,1),&
               ubound(rsolver%RtempVectors,1)
          call lsysbl_releaseVector(rsolver%RtempVectors(i))
        end do
        deallocate(rsolver%RtempVectors)
        nullify(rsolver%RtempVectors)
      end if

      ! Release coarse grid solver if required
      if (associated(rsolver%p_solverCoarsegrid)) then
        call solver_releaseSolver(rsolver%p_solverCoarsegrid)
        deallocate(rsolver%p_solverCoarsegrid)
        nullify(rsolver%p_solverCoarsegrid)
      end if

      ! Release all smoother subnodes
      if (associated(rsolver%p_smoother)) then
        do i = lbound(rsolver%p_smoother,1),&
               ubound(rsolver%p_smoother,1)
          call solver_releaseSolver(rsolver%p_smoother(i))
        end do
        deallocate(rsolver%p_smoother)
        nullify(rsolver%p_smoother)
      end if
    end subroutine release_solverMultigrid

    !*************************************************************
    ! Release UMFPACK solver

    subroutine release_solverUMFPACK(rsolver)
      type(t_solverUMFPACK), intent(inout) :: rsolver

      integer :: i

      ! Release UMFPACK handles
      if (rsolver%isymbolic .ne. 0) call UMF4FSYM(rsolver%isymbolic)
      if (rsolver%inumeric  .ne. 0) call UMF4FNUM(rsolver%inumeric)

      if (associated(rsolver%IsymbolicBlock)) then
        do i = lbound(rsolver%IsymbolicBlock,1),&
               ubound(Rsolver%IsymbolicBlock,1)
          if (rsolver%IsymbolicBlock(i) .ne. 0) call UMF4FSYM(rsolver%IsymbolicBlock(i))
        end do
        deallocate(rsolver%IsymbolicBlock)
      end if

      if (associated(rsolver%InumericBlock)) then
        do i = lbound(rsolver%InumericBlock,1),&
               ubound(Rsolver%InumericBlock,1)
          if (rsolver%InumericBlock(i) .ne. 0) call UMF4FNUM(rsolver%InumericBlock(i))
        end do
        deallocate(rsolver%InumericBlock)
      end if

      ! Release matrix
      call lsysbl_releaseMatrix(rsolver%rmatrix)

      ! Release temporal vector
      call lsysbl_releaseVector(rsolver%rtempVector)
    end subroutine release_solverUMFPACK

    !*************************************************************
    ! Release Jacobi solver

    subroutine release_solverJacobi(rsolver)
      type(t_solverJacobi), intent(inout) :: rsolver

      ! Release matrix
      call lsysbl_releaseMatrix(rsolver%rmatrix)

      ! Release temporal vector
      call lsysbl_releaseVector(rsolver%rtempVector)
    end subroutine release_solverJacobi

    !*************************************************************
    ! Release (S)SOR solver

    subroutine release_solverSSOR(rsolver)
      type(t_solverSSOR), intent(inout) :: rsolver

      ! Release matrix
      call lsysbl_releaseMatrix(rsolver%rmatrix)

      ! Release temporal vector
      call lsysbl_releaseVector(rsolver%rtempVector)
    end subroutine release_solverSSOR

    !*************************************************************
    ! Release BiCGSTAB solver

    subroutine release_solverBiCGSTAB(rsolver)
      type(t_solverBiCGSTAB), intent(inout) :: rsolver

      ! local variables
      integer :: i

      ! Release matrix
      call lsysbl_releaseMatrix(rsolver%rmatrix)

      ! Release temporal vectors
      do i = lbound(rsolver%RtempVectors,1),&
             ubound(rsolver%RtempVectors,1)
        call lsysbl_releaseVector(rsolver%RtempVectors(i))
      end do

      ! Release preconditioner if required
      if (associated(rsolver%p_precond)) then
        call solver_releaseSolver(rsolver%p_precond)
        deallocate(rsolver%p_precond)
        nullify(rsolver%p_precond)
      end if
    end subroutine release_solverBiCGSTAB

    !*************************************************************
    ! Release GMRES solver

    subroutine release_solverGMRES(rsolver)
      type(t_solverGMRES), intent(inout) :: rsolver

      ! local variables
      integer :: i

      ! Release matrix
      call lsysbl_releaseMatrix(rsolver%rmatrix)

      ! Release handles
      call storage_free(rsolver%h_Dh)
      call storage_free(rsolver%h_Dc)
      call storage_free(rsolver%h_Ds)
      call storage_free(rsolver%h_Dq)

      ! Release temporal vectors
      do i = lbound(rsolver%rv,1),&
             ubound(rsolver%rv,1)
        call lsysbl_releaseVector(rsolver%rv(i))
      end do

      do i = lbound(rsolver%rz,1),&
             ubound(rsolver%rz,1)
        call lsysbl_releaseVector(rsolver%rz(i))
      end do

      ! Deallocate subarray of temporal vectors
      deallocate(rsolver%rv)
      deallocate(rsolver%rz)
      nullify(rsolver%rv)
      nullify(rsolver%rz)

      ! Release preconditioner if required
      if (associated(rsolver%p_precond)) then
        call solver_releaseSolver(rsolver%p_precond)
        deallocate(rsolver%p_precond)
        nullify(rsolver%p_precond)
      end if
    end subroutine release_solverGMRES

    !*************************************************************
    ! Release ILU solver

    recursive subroutine release_solverILU(rsolver)
      type(t_solverILU), intent(inout) :: rsolver

      integer :: i

      ! Release structure for block ILU solver
      if (associated(rsolver%p_rsolverBlockILU)) then
        do i = lbound(rsolver%p_rsolverBlockILU,1),&
               ubound(rsolver%p_rsolverBlockILU,1)
          call release_solverILU(rsolver%p_rsolverBlockILU(i))
        end do
        deallocate(rsolver%p_rsolverBlockILU)
      end if

      ! Release handles
      if (rsolver%h_Idata .ne. ST_NOHANDLE) call storage_free(rsolver%h_Idata)
      if (rsolver%h_Ddata .ne. ST_NOHANDLE) call storage_free(rsolver%h_Ddata)

      ! Release matrix
      call lsysbl_releaseMatrix(rsolver%rmatrix)

      ! Release temporal vector
      call lsysbl_releaseVector(rsolver%rtempVector)
    end subroutine release_solverILU

    !*************************************************************
    ! Release defect correction algorithm

    subroutine release_solverDefcor(rprecond)
      type(t_solverDefcor), intent(inout) :: rprecond

      ! local variables
      integer :: i

      ! Release temporal vectors
      do i = lbound(rprecond%RtempVectors,1),&
             ubound(rprecond%RtempVectors,1)
        call lsysbl_releaseVector(rprecond%RtempVectors(i))
      end do
    end subroutine release_solverDefcor

    !*************************************************************
    ! Release Newton`s algorithm

    subroutine release_solverNewton(rprecond)
      type(t_solverNewton), intent(inout) :: rprecond

      ! local variables
      integer :: i

      ! Release temporal vectors
      do i = lbound(rprecond%RtempVectors,1),&
             ubound(rprecond%RtempVectors,1)
        call lsysbl_releaseVector(rprecond%RtempVectors(i))
      end do
    end subroutine release_solverNewton

  end subroutine solver_releaseSolver

  ! ***************************************************************************

!<subroutine>

  recursive subroutine solver_infoSolver(rsolver, bprintInternal)

!<description>
    ! This subroutine prints information about the solver. If bprintInternal is
    ! specified than a detailed output of the internal solver structure us printed.
!</description>

!<input>
    ! Solver structure
    type(t_solver), intent(in) :: rsolver

    ! OPTIONAL: If true, a detailed output of the internal solver structure is printed
    logical, intent(in), optional :: bprintInternal
!</input>
!</subroutine>

    ! local variable
    integer :: i
    logical :: bprint


    ! What should be printed?
    bprint = .false.
    if (present(bprintInternal)) bprint=bprintInternal

    ! Output general information
    call output_line('Solver:')
    call output_line('-------')
    call output_line('Name: '//trim(rsolver%ssolverName))

    select case(rsolver%csolverType)
    case (SV_COUPLED)
      call output_line('Number of coupled components:                 '//&
          trim(sys_siL(size(rsolver%p_solverSubnode),3)))

    case (SV_FMG)
      call output_line('Number of full multigrid steps:               '//&
          trim(sys_siL(rsolver%niterations,15)))
      if (associated(rsolver%p_solverMultigrid%p_solverCoarsegrid)) then
        if (rsolver%niterations .ne. 0)&
            call output_line('Number of nonlinear steps per multigrid step: '&
            //trim(sys_siL(int(rsolver%p_solverMultigrid%p_solverCoarsegrid%niterations/&
            real(rsolver%niterations, DP)),15)))
      end if

    case (SV_NONLINEARMG)
      call output_line('Number of nonlinear multigrid steps:          '//&
          trim(sys_siL(rsolver%niterations, 15)))
      if (associated(rsolver%p_solverMultigrid%p_solverCoarsegrid)) then
        if (rsolver%niterations .ne. 0)&
            call output_line('Number of nonlinear steps per multigrid step: '//&
            trim(sys_siL(int(rsolver%p_solverMultigrid%p_solverCoarsegrid%niterations/&
            real(rsolver%niterations,DP)),15)))
      end if

    case (SV_LINEARMG)
      call output_line('Number of linear steps:                       '//&
          trim(sys_siL(rsolver%niterations,15)))
      if (associated(rsolver%p_solverMultigrid%p_solverCoarsegrid)) then
        if (rsolver%niterations .ne. 0)&
            call output_line('Number of linear steps per multigrid step:    '//&
            trim(sys_siL(int(rsolver%p_solverMultigrid%p_solverCoarsegrid%niterations/&
            real(rsolver%niterations,DP)),15)))
      end if

    case(SV_NONLINEAR)
      call output_line('Number of nonlinear steps:                    '//&
          trim(sys_siL(rsolver%niterations,15)))
      if (associated(rsolver%p_solverSubnode)) then
        if (rsolver%niterations .ne. 0)&
            call output_line('Number of linear steps per nonlinear step:    '//&
            trim(sys_siL(int(rsolver%p_solverSubnode(1)%niterations/&
            real(rsolver%niterations,DP)),15)))
      end if

    case(SV_LINEAR)
      call output_line('Number of linear steps:                       '//&
          trim(sys_siL(rsolver%niterations,15)))

    case default
      call output_line('Unsupported solver type!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'solver_infoSolver')
      call sys_halt()
    end select

    ! Output internal information
    if (bprint) then
      call output_line('csolverType:                                  '//trim(sys_siL(rsolver%csolverType,3)))
      call output_line('isolver:                                      '//trim(sys_siL(rsolver%isolver,3)))
      call output_line('iprecond:                                     '//trim(sys_siL(rsolver%iprecond,3)))
      call output_line('ioutputlevel:                                 '//trim(sys_siL(rsolver%ioutputLevel,3)))
      call output_line('iresNorm:                                     '//trim(sys_siL(rsolver%iresNorm,3)))
      call output_line('istoppingCriterion:                           '//trim(sys_siL(rsolver%istoppingCriterion,3)))
      call output_line('nminIterations:                               '//trim(sys_siL(rsolver%nminIterations,5)))
      call output_line('nmaxIterations:                               '//trim(sys_siL(rsolver%nmaxIterations,5)))
      call output_line('domega:                                       '//trim(sys_sdL(abs(rsolver%domega),2))//&
                                                                         trim(merge(' [implicit]', ' [explicit]',&
                                                                         rsolver%domega .lt. 0.0_DP)))
      call output_line('depsRel:                                      '//trim(sys_sdL(rsolver%depsRel,5)))
      call output_line('depsAbs:                                      '//trim(sys_sdL(rsolver%depsAbs,5)))
      call output_line('ddivRel:                                      '//trim(sys_sdL(rsolver%ddivRel,5)))
      call output_line('ddivAbs:                                      '//trim(sys_sdL(rsolver%ddivAbs,5)))
      call output_line('drhsZero:                                     '//trim(sys_sdL(rsolver%drhsZero,5)))
      call output_line('ddefZero:                                     '//trim(sys_sdL(rsolver%ddefZero,5)))
      call output_line('depsStag:                                     '//trim(sys_sdL(rsolver%depsStag,5)))
      call output_line('boundary condition:                           '//merge('set    ', 'missing',&
                                                                         associated(rsolver%rboundaryCondition)))

      call output_line('istatus:                                      '//trim(solver_getStatus(rsolver)))
      call output_line('iiterations:                                  '//trim(sys_siL(rsolver%iiterations,5)))
      call output_line('niterations:                                  '//trim(sys_siL(rsolver%niterations,5)))
      call output_line('dinitialRHS:                                  '//trim(sys_sdL(rsolver%dinitialRHS,5)))
      call output_line('dinitialDefect:                               '//trim(sys_sdL(rsolver%dinitialDefect,5)))
      call output_line('dfinalDefect:                                 '//trim(sys_sdL(rsolver%dfinalDefect,5)))
      call output_line('dconvergenceRate:                             '//trim(sys_sdL(rsolver%dconvergenceRate,5)))

      ! UMFPACK4 solver
      if (associated(rsolver%p_solverUMFPACK)) then
        call output_lbrk()
        call output_line('>>> UMFPACK subnode:')
        call output_line('--------------------')
        call output_line('isymbolic: '//trim(sys_siL(int(rsolver%p_solverUMFPACK%isymbolic),5)))
        call output_line('inumeric:  '//trim(sys_siL(int(rsolver%p_solverUMFPACK%inumeric),5)))
        call output_line('>>> system matrix:')
        call output_line('------------------')
        call lsysbl_infoMatrix(rsolver%p_solverUMFPACK%rmatrix)
        call output_line('>>> temporal vector:')
        call output_line('--------------------')
        call lsysbl_infoVector(rsolver%p_solverUMFPACK%rtempVector)
      end if

      ! Jacobi solver
      if (associated(rsolver%p_solverJacobi)) then
        call output_lbrk()
        call output_line('>>> Jacobi solver:')
        call output_line('------------------')
        call output_line('>>> System matrix:')
        call output_line('------------------')
        call lsysbl_infoMatrix(rsolver%p_solverJacobi%rmatrix)
      end if

      ! (S)SOR solver
      if (associated(rsolver%p_solverSSOR)) then
        call output_lbrk()
        call output_line('>>> (S)SOR solver:')
        call output_line('------------------')
        call output_line('>>> System matrix:')
        call output_line('------------------')
        call lsysbl_infoMatrix(rsolver%p_solverSSOR%rmatrix)
      end if

      ! BiCGSTAB solver
      if (associated(rsolver%p_solverBiCGSTAB)) then
        call output_lbrk()
        call output_line('>>> BiCGSTAB solver:')
        call output_line('--------------------')
        call output_line('>>> Temporal vectors:')
        call output_line('---------------------')
        do i = lbound(rsolver%p_solverBiCGSTAB%RtempVectors,1),&
               ubound(rsolver%p_solverBiCGSTAB%RtempVectors,1)
          call lsysbl_infoVector(rsolver%p_solverBiCGSTAB%RtempVectors(i))
        end do
        call output_line('>>> System matrix:')
        call output_line('------------------')
        call lsysbl_infoMatrix(rsolver%p_solverBiCGSTAB%rmatrix)
        if (associated(rsolver%p_solverBiCGSTAB%p_precond)) then
          call output_line('>>> Preconditioner:')
          call output_line('-------------------')
          call solver_infoSolver(rsolver%p_solverBiCGSTAB%p_precond, bprint)
        end if
      end if

      ! GMRES solver
      if (associated(rsolver%p_solverGMRES)) then
        call output_lbrk()
        call output_line('>>> GMRES solver:')
        call output_line('-----------------')
        call output_line('h_Dh: '//trim(sys_siL(rsolver%p_solverGMRES%h_Dh,5)))
        call output_line('h_Dc: '//trim(sys_siL(rsolver%p_solverGMRES%h_Dc,5)))
        call output_line('h_Ds: '//trim(sys_siL(rsolver%p_solverGMRES%h_Ds,5)))
        call output_line('h_Dq: '//trim(sys_siL(rsolver%p_solverGMRES%h_Dq,5)))
        if (associated(rsolver%p_solverGMRES%rv)) then
          call output_line('>>> Temporal vector RV (first entry):')
          call output_line('-------------------------------------')
          call lsysbl_infoVector(rsolver%p_solverGMRES%rv(1))
        end if
        if (associated(rsolver%p_solverGMRES%rz)) then
          call output_line('>>> Temporal vector RZ (first entry):')
          call output_line('-------------------------------------')
          call lsysbl_infoVector(rsolver%p_solverGMRES%rz(1))
        end if
        call output_line('>>> System matrix:')
        call output_linE('------------------')
        call lsysbl_infoMatrix(rsolver%p_solverGMRES%rmatrix)
        if (associated(rsolver%p_solverGMRES%p_precond)) then
          call output_line('>>> Preconditioner:')
          call output_line('-------------------')
          call solver_infoSolver(rsolver%p_solverGMRES%p_precond, bprint)
        end if
      end if

      ! ILU solver
      if (associated(rsolver%p_solverILU)) then
        call output_lbrk()
        call output_line('>>> ILU solver:')
        call output_line('---------------')
        call output_line('ifill:   '//trim(sys_siL(rsolver%p_solverILU%ifill,3)))
        call output_line('domega:  '//trim(sys_sdL(rsolver%p_solverILU%domega,5)))
        call output_line('depsILU: '//trim(sys_sdL(rsolver%p_solverILU%depsILU,5)))
        call output_line('>>> system matrix:')
        call output_line('------------------')
        call lsysbl_infoMatrix(rsolver%p_solverILU%rmatrix)
      end if

      ! Defcor preconditioner
      if (associated(rsolver%p_solverDefcor)) then
        call output_lbrk()
        call output_line('>>> Defect correction preconditioner:')
        call output_line('-------------------------------------')
        call output_line('>>> Temporal vectors:')
        call output_line('---------------------')
        do i = lbound(rsolver%p_solverDefcor%RtempVectors,1),&
               ubound(rsolver%p_solverDefcor%RtempVectors,1)
          call lsysbl_infoVector(rsolver%p_solverDefcor%RtempVectors(i))
        end do
      end if

      ! Newton preconditioner
      if (associated(rsolver%p_solverNewton)) then
        call output_lbrk()
        call output_line('>>> Newton preconditioner:')
        call output_line('--------------------------')
        call output_line('icheckSufficientDecrease: '//trim(sys_siL(rsolver%p_solverNewton%icheckSufficientDecrease,3)))
        call output_line('nmaxBacktrackingSteps:    '//trim(sys_siL(rsolver%p_solverNewton%nmaxBacktrackingSteps,5)))
        call output_line('iupdateFrequency:         '//trim(sys_siL(rsolver%p_solverNewton%iupdateFrequency,5)))
        call output_line('dforcingStrategy:         '//trim(sys_sdL(rsolver%p_solverNewton%dforcingStrategy,5)))
        call output_line('dperturbationStrategy:    '//trim(sys_sdL(rsolver%p_solverNewton%dperturbationStrategy,5)))
        call output_line('>>> Temporal vectors:')
        call output_line('---------------------')
        do i = lbound(rsolver%p_solverNewton%RtempVectors,1),&
               ubound(rsolver%p_solverNewton%RtempVectors,1)
          call lsysbl_infoVector(rsolver%p_solverNewton%RtempVectors(i))
        end do
      end if

      ! Solver subnode
      if (associated(rsolver%p_solverSubnode)) then
        do i = 1, size(rsolver%p_solverSubnode)
          call output_lbrk()
          call output_line('>>> Generic solver subnode: '//trim(sys_siL(i,3)))
          call output_line('-------------------------------')
          call solver_infoSolver(rsolver%p_solverSubnode(i), bprint)
        end do
      end if

      ! Multigrid solver
      if (associated(rsolver%p_solverMultigrid)) then
        call output_lbrk()
        call output_line('>>> Multigrid solver:')
        select case(rsolver%p_solverMultigrid%icycle)
        case (0)
          call output_line('cycle: F-cycle')
        case (-1,1)
          call output_line('cycle: V-cycle')
        case (-2,2)
          call output_line('cycle: W-cycle')
        case DEFAULT
          call output_line('cycle: unknown !!!')
        end select
        call output_line('nlmin:         '//trim(sys_siL(rsolver%p_solverMultigrid%initialNlmin,5))//', '//&
                                            trim(sys_siL(rsolver%p_solverMultigrid%nlmin,5)))
        call output_line('nlmax          '//trim(sys_siL(rsolver%p_solverMultigrid%initialNlmax,5))//', '//&
                                            trim(sys_siL(rsolver%p_solverMultigrid%nlmax,5)))
        call output_line('ilmin:         '//trim(sys_siL(rsolver%p_solverMultigrid%ilmin,5)))
        call output_line('ilmax:         '//trim(sys_siL(rsolver%p_solverMultigrid%ilmax,5)))
        call output_line('csmootherType: '//trim(sys_siL(rsolver%p_solverMultigrid%csmootherType,5)))
        call output_line('npresmooth:    '//trim(sys_siL(rsolver%p_solverMultigrid%npresmooth,5)))
        call output_line('npostmooth:    '//trim(sys_siL(rsolver%p_solverMultigrid%npostsmooth,5)))
        if (associated(rsolver%p_solverMultigrid%rmatrix)) then
          call output_line('>>> System matrices:')
          call output_line('--------------------')
          do i = lbound(rsolver%p_solverMultigrid%rmatrix,1),&
                 ubound(rsolver%p_solverMultigrid%rmatrix,1)
            call lsysbl_infoMatrix(rsolver%p_solverMultigrid%rmatrix(i))
          end do
        end if
        call output_line('>>> Temporal vectors:')
        call output_line('---------------------')
        if (associated(rsolver%p_solverMultigrid%RtempVectors)) then
          do i = lbound(rsolver%p_solverMultigrid%RtempVectors,1),&
                 ubound(rsolver%p_solverMultigrid%RtempVectors,1)
            call lsysbl_infoVector(rsolver%p_solverMultigrid%RtempVectors(i))
          end do
        end if

        ! Smoothers
        if (associated(rsolver%p_solverMultigrid%p_smoother)) then
          do i = lbound(rsolver%p_solverMultigrid%p_smoother,1),&
                 ubound(rsolver%p_solverMultigrid%p_smoother,1)
            call output_lbrk()
            call output_line('>>> Smoother on Level '//trim(sys_siL(i,3)))
            call solver_infoSolver(rsolver%p_solverMultigrid%p_smoother(i), bprint)
          end do
        end if

        ! Coarse-grid solver
        if (associated(rsolver%p_solverMultigrid%p_solverCoarsegrid)) then
          call output_lbrk()
          call output_line('>>> Coarsegrid solver:')
          call output_line('----------------------')
          call solver_infoSolver(rsolver%p_solverMultigrid%p_solverCoarsegrid, bprint)
        end if
      end if

    else

      ! UMFPACK4 solver
      if (associated(rsolver%p_solverUMFPACK)) then
        call output_lbrk()
        call output_line('>>> UMFPACK subnode:')
        call output_line('--------------------')
      end if

      ! Jacobi solver
      if (associated(rsolver%p_solverJacobi)) then
        call output_lbrk()
        call output_line('>>> Jacobi solver:')
        call output_line('------------------')
      end if

      ! (S)SOR solver
      if (associated(rsolver%p_solverSSOR)) then
        call output_lbrk()
        call output_line('>>> (S)SOR solver:')
        call output_line('------------------')
      end if

      ! BiCGSTAB solver
      if (associated(rsolver%p_solverBiCGSTAB)) then
        call output_lbrk()
        call output_line('>>> BiCGSTAB solver:')
        call output_line('--------------------')
        if (associated(rsolver%p_solverBiCGSTAB%p_precond)) then
        call output_lbrk()
          call output_line('>>> Preconditioner:')
          call output_line('-------------------')
          call solver_infoSolver(rsolver%p_solverBiCGSTAB%p_precond)
        end if
      end if

      ! GMRES solver
      if (associated(rsolver%p_solverGMRES)) then
        call output_lbrk()
        call output_line('>>> GMRES solver:')
        call output_line('-----------------')
        if (associated(rsolver%p_solverGMRES%p_precond)) then
          call output_lbrk()
          call output_line('>>> Preconditioner:')
          call output_line('-------------------')
          call solver_infoSolver(rsolver%p_solverGMRES%p_precond)
        end if
      end if

      ! ILU solver
      if (associated(rsolver%p_solverILU)) then
        call output_lbrk()
        call output_line('>>> ILU solver:')
        call output_line('---------------')
      end if

      ! Defcor preconditioner
      if (associated(rsolver%p_solverDefcor)) then
        call output_lbrk()
        call output_line('>>> Defect correction preconditioner:')
        call output_line('-------------------------------------')
      end if

      ! Newton preconditioner
      if (associated(rsolver%p_solverNewton)) then
        call output_lbrk()
        call output_line('>>> Newton preconditioner:')
        call output_line('--------------------------')
      end if

      ! Solver subnode
      if (associated(rsolver%p_solverSubnode)) then
        do i = 1, size(rsolver%p_solverSubnode)
          call output_lbrk()
          call output_line('>>> Generic solver subnode: '//trim(sys_siL(i,3)))
          call output_line('-------------------------------')
          call solver_infoSolver(rsolver%p_solverSubnode(i), bprint)
        end do
      end if

      ! Multigrid solver
      if (associated(rsolver%p_solverMultigrid)) then
        call output_lbrk()
        call output_line('>>> Multigrid solver:')
        call output_line('---------------------')

        ! Smoothers
        if (associated(rsolver%p_solverMultigrid%p_smoother)) then
          do i = lbound(rsolver%p_solverMultigrid%p_smoother,1),&
                 ubound(rsolver%p_solverMultigrid%p_smoother,1)
            call output_lbrk()
            call output_line('>>> Smoother on Level '//trim(sys_siL(i,3)))
            call output_line('-------------------------')
            call solver_infoSolver(rsolver%p_solverMultigrid%p_smoother(i), bprint)
          end do
        end if

        ! Coarse-grid solver
        if (associated(rsolver%p_solverMultigrid%p_solverCoarsegrid)) then
          call output_lbrk()
          call output_line('>>> Coarsegrid solver:')
          call output_line('----------------------')
          call solver_infoSolver(rsolver%p_solverMultigrid%p_solverCoarsegrid, bprint)
        end if
      end if

    end if

    call output_separator(OU_SEP_PLUS)
    call output_lbrk()
  end subroutine solver_infoSolver

  ! ***************************************************************************

!<subroutine>

  recursive subroutine solver_resetSolver(rsolver, bresetStatistics)

!<description>
    ! This subroutine resets the solver to its initial values
!</description>

!<input>
    ! If true, the statistical output parameters are reset
    logical, intent(in) :: bresetStatistics
!</input>

!<inputoutput>
    ! solver structure
    type(t_solver), intent(inout) :: rsolver
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: i

    ! Reset statistical data (if required)
    if (bresetStatistics) then
      rsolver%dinitialRHS      = 0.0_DP
      rsolver%dinitialDefect   = 0.0_DP
      rsolver%dfinalDefect     = 0.0_DP
      rsolver%dconvergenceRate = 0.0_DP
      rsolver%istatus          = 0
      rsolver%iiterations      = 0
      rsolver%niterations      = 0
    end if

    ! Reset generic solver subnode
    if (associated(rsolver%p_solverSubnode)) then
      do i = 1, size(rsolver%p_solverSubnode)
        call solver_resetSolver(rsolver%p_solverSubnode(i), bresetStatistics)
      end do
    end if

    ! Reset multigrid solver
    if (associated(rsolver%p_solverMultigrid)) then

      ! Reset multigrid levels
      rsolver%p_solverMultigrid%nlmin = rsolver%p_solverMultigrid%initialNlmin
      rsolver%p_solverMultigrid%nlmax = rsolver%p_solverMultigrid%initialNlmax

      ! Reset coarsegrid solver
      if (associated(rsolver%p_solverMultigrid%p_solverCoarsegrid)) then
        call solver_resetSolver(rsolver%p_solverMultigrid%p_solverCoarsegrid,&
                                bresetStatistics)
      end if

      ! Reset smoother
      if (associated(rsolver%p_solverMultigrid%p_smoother)) then
        do i = lbound(rsolver%p_solverMultigrid%p_smoother,1),&
               ubound(rsolver%p_solverMultigrid%p_smoother,1)
          call solver_resetSolver(rsolver%p_solverMultigrid%p_smoother(i),&
                                  bresetStatistics)
        end do
      end if
    end if

    ! Reset preconditioners
    if (associated(rsolver%p_solverBiCGSTAB)) then
      if (associated(rsolver%p_solverBiCGSTAB%p_precond)) then
        call solver_resetSolver(rsolver%p_solverBiCGSTAB%p_precond,&
                                bresetStatistics)
      end if
    end if

    if (associated(rsolver%p_solverGMRES)) then
      if (associated(rsolver%p_solverGMRES%p_precond)) then
        call solver_resetSolver(rsolver%p_solverGMRES%p_precond,&
                                bresetStatistics)
      end if
    end if
  end subroutine solver_resetSolver

  ! ***************************************************************************

!<subroutine>

  recursive subroutine solver_removeTempFromSolver(rsolver)

!<description>
    ! This subroutine recursively removes all temporal memory from the
    ! given solver structure and all of its associated subnodes
!</description>

!<inputoutput>
    ! solver structure
    type(t_solver), intent(inout) :: rsolver
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: i

    ! This type of solver does not have and temporals.
    ! Hence, we can directly proceed to the subnodes.

    ! Generic solver subnode
    if (associated(rsolver%p_solverSubnode)) then
      do i = 1, size(rsolver%p_solverSubnode)
        call solver_removeTempFromSolver(rsolver%p_solverSubnode(i))
      end do
    end if

    ! Multigrid solver
    if (associated(rsolver%p_solverMultigrid)) then
      call remove_solverMultigrid(rsolver%p_solverMultigrid)
    end if

    ! UMFPACK solver
    if (associated(rsolver%p_solverUMFPACK)) then
      call remove_solverUMFPACK(rsolver%p_solverUMFPACK)
    end if

    ! Jacobi solver
    if (associated(rsolver%p_solverJacobi)) then
      call remove_solverJacobi(rsolver%p_solverJacobi)
    end if

    ! (S)SOR solver
    if (associated(rsolver%p_solverSSOR)) then
      call remove_solverSSOR(rsolver%p_solverSSOR)
    end if

    ! BiCGSTAB solver
    if (associated(rsolver%p_solverBiCGSTAB)) then
      call remove_solverBiCGSTAB(rsolver%p_solverBiCGSTAB)
    end if

    ! GMRES solver
    if (associated(rsolver%p_solverGMRES)) then
      call remove_solverGMRES(rsolver%p_solverGMRES)
    end if

    ! ILU solver
    if (associated(rsolver%p_solverILU)) then
      call remove_solverILU(rsolver%p_solverILU)
    end if

    ! Defcor preconditioner
    if (associated(rsolver%p_solverDefcor)) then
      call remove_solverDefcor(rsolver%p_solverDefcor)
    end if

    ! Newton preconditioner
    if (associated(rsolver%p_solverNewton)) then
      call remove_solverNewton(rsolver%p_solverNewton)
    end if

  contains

    ! Here, the real working routine follows

    !*************************************************************
    ! Remove temporal data from multigrid solver

    subroutine remove_solverMultigrid(rsolver)
      type(t_solverMultigrid), intent(inout) :: rsolver

      ! local variables
      integer :: i

      ! Release matrix
      if (associated(rsolver%rmatrix)) then
        do i = lbound(rsolver%rmatrix,1),&
               ubound(rsolver%rmatrix,1)
          call lsysbl_releaseMatrix(rsolver%rmatrix(i))
        end do
      end if

      ! Release temporal vectors
      if (associated(rsolver%RtempVectors)) then
        do i = lbound(rsolver%RtempVectors,1),&
               ubound(rsolver%RtempVectors,1)
          call lsysbl_releaseVector(rsolver%RtempVectors(i))
        end do
      end if

      ! Remove temporal from coarsegrid solver
      if (associated(rsolver%p_solverCoarsegrid)) then
        call solver_removeTempFromSolver(rsolver%p_solverCoarsegrid)
      end if

      ! Remove temporals from smoother
      if (associated(rsolver%p_smoother)) then
        do i = lbound(rsolver%p_smoother,1),&
               ubound(rsolver%p_smoother,1)
          call solver_removeTempFromSolver(rsolver%p_smoother(i))
        end do
      end if
    end subroutine remove_solverMultigrid

    !*************************************************************
    ! Remove temporal data from UMFPACK solver

    subroutine remove_solverUMFPACK(rsolver)
      type(t_solverUMFPACK), intent(inout) :: rsolver

      integer :: i

      ! Release UMFPACK handles
      if (rsolver%isymbolic .ne. 0) call UMF4FSYM(rsolver%isymbolic)
      if (rsolver%inumeric  .ne. 0) call UMF4FNUM(rsolver%inumeric)

      if (associated(rsolver%IsymbolicBlock)) then
        do i = lbound(rsolver%IsymbolicBlock,1),&
               ubound(Rsolver%IsymbolicBlock,1)
          if (rsolver%IsymbolicBlock(i) .ne. 0) call UMF4FSYM(rsolver%IsymbolicBlock(i))
        end do
      end if

      if (associated(rsolver%InumericBlock)) then
        do i = lbound(rsolver%InumericBlock,1),&
               ubound(Rsolver%InumericBlock,1)
          if (rsolver%InumericBlock(i) .ne. 0) call UMF4FNUM(rsolver%InumericBlock(i))
        end do
      end if

      ! Release matrix
      call lsysbl_releaseMatrix(rsolver%rmatrix)

      ! Release temporal vector
      call lsysbl_releaseVector(rsolver%rtempVector)
    end subroutine remove_solverUMFPACK

    !*************************************************************
    ! Remove temporal data from Jacobi solver

    subroutine remove_solverJacobi(rsolver)
      type(t_solverJacobi), intent(inout) :: rsolver

      ! Release matrix
      call lsysbl_releaseMatrix(rsolver%rmatrix)

      ! Release temporal vector
      call lsysbl_releaseVector(rsolver%rtempVector)
    end subroutine remove_solverJacobi

    !*************************************************************
    ! Remove temporal data from (S)SOR solver

    subroutine remove_solverSSOR(rsolver)
      type(t_solverSSOR), intent(inout) :: rsolver

      ! Release matrix
      call lsysbl_releaseMatrix(rsolver%rmatrix)

      ! Release temporal vector
      call lsysbl_releaseVector(rsolver%rtempVector)
    end subroutine remove_solverSSOR

    !*************************************************************
    ! Remove temporal data from BiCGSTAB solver

    subroutine remove_solverBiCGSTAB(rsolver)
      type(t_solverBiCGSTAB), intent(inout) :: rsolver

      ! local variables
      integer :: i

      ! Release matrix
      call lsysbl_releaseMatrix(rsolver%rmatrix)

      ! Release temporal vectors
      do i = lbound(rsolver%RtempVectors,1),&
             ubound(rsolver%RtempVectors,1)
        call lsysbl_releaseVector(rsolver%RtempVectors(i))
      end do

      ! Remove temporals from preconditioner
      if (associated(rsolver%p_precond)) then
        call solver_removeTempFromSolver(rsolver%p_precond)
      end if
    end subroutine remove_solverBiCGSTAB

    !*************************************************************
    ! Remove temporal data from GMRES solver

    subroutine remove_solverGMRES(rsolver)
      type(t_solverGMRES), intent(inout) :: rsolver

      ! local variables
      integer :: i

      ! Release matrix
      call lsysbl_releaseMatrix(rsolver%rmatrix)

      ! Release temporal vectors
      do i = lbound(rsolver%rv,1),&
             ubound(rsolver%rv,1)
        call lsysbl_releaseVector(rsolver%rv(i))
      end do

      do i = lbound(rsolver%rz,1),&
             ubound(rsolver%rz,1)
        call lsysbl_releaseVector(rsolver%rz(i))
      end do

      ! Remove temporals from preconditioner
      if (associated(rsolver%p_precond)) then
        call solver_removeTempFromSolver(rsolver%p_precond)
      end if
    end subroutine remove_solverGMRES

    !*************************************************************
    ! Remove temporal data from ILU solver

    subroutine remove_solverILU(rsolver)
      type(t_solverILU), intent(inout) :: rsolver

      ! Release ILU handles
      if (rsolver%h_Idata .ne. ST_NOHANDLE) call storage_free(rsolver%h_Idata)
      if (rsolver%h_Ddata .ne. ST_NOHANDLE) call storage_free(rsolver%h_Ddata)

      ! Release matrix
      call lsysbl_releaseMatrix(rsolver%rmatrix)

      ! Release temporal vector
      call lsysbl_releaseVector(rsolver%rtempVector)
    end subroutine remove_solverILU

    !*************************************************************
    ! Remove temporal data from defect correction algorithm

    subroutine remove_solverDefcor(rsolver)
      type(t_solverDefcor), intent(inout) :: rsolver

      ! local variables
      integer :: i

      ! Release temporal vectors
      do i = lbound(rsolver%RtempVectors,1),&
             ubound(rsolver%RtempVectors,1)
        call lsysbl_releaseVector(rsolver%RtempVectors(i))
      end do
    end subroutine remove_solverDefcor

    !*************************************************************
    ! Remove temporal data from Newton`s algorithm

    subroutine remove_solverNewton(rsolver)
      type(t_solverNewton), intent(inout) :: rsolver

      ! local variables
      integer :: i

      ! Release temporal vectors
      do i = lbound(rsolver%RtempVectors,1),&
             ubound(rsolver%RtempVectors,1)
        call lsysbl_releaseVector(rsolver%RtempVectors(i))
      end do
    end subroutine remove_solverNewton

  end subroutine solver_removeTempFromSolver

  ! ***************************************************************************

!<subroutine>

  recursive subroutine solver_adjustHierarchy(rsolver, nlminOpt, nlmaxOpt)

!<description>
    ! This subroutine adjusts the multigrid hierarchy of the complete
    ! solver structure. This may be required, if the user provides
    ! mixed up settings.
!</description>

!<input>
    ! OPTIONAL: minimum multigrid level to use
    integer, intent(in), optional :: nlminOpt

    ! OPTIONAL: maximum multigrid level to use
    integer, intent(in), optional :: nlmaxOpt
!</input>

!<inputoutput>
    ! Solver structure
    type(t_solver), intent(inout) :: rsolver
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_solverMultigrid), pointer :: p_solverMultigrid
    type(t_solver) :: rsolverTemplateSmoother
    integer :: i, nlmin, nlmax, ilbound, iubound

    ! Determine NLMIN and NLMAX for structure
    if (present(nlminOpt)) then
      nlmin = nlminOpt
    else
      nlmin = solver_getMinimumMultigridlevel(rsolver)
    end if
    if (present(nlmaxOpt)) then
      nlmax = nlmaxOpt
    else
      nlmax = solver_getMaximumMultigridlevel(rsolver)
    end if

    ! Do we have a solver subnode?
    if (associated(rsolver%p_solverSubnode)) then
      do i = 1, size(rsolver%p_solverSubnode)
        call solver_adjustHierarchy(rsolver%p_solverSubnode(i), nlmin, nlmax)
      end do
    end if

    ! Do we have a multigrid structure?
    if (associated(rsolver%p_solverMultigrid)) then
      p_solverMultigrid => rsolver%p_solverMultigrid

      ! Should multigrid be used or is NLMIN = NLMAX
      if (p_solverMultigrid%initialNlmin .eq. &
          p_solverMultigrid%initialNlmax) then

        ! Check maximum multigrid level
        if (p_solverMultigrid%nlmax .ne. nlmax) then
          p_solverMultigrid%nlmax = nlmax
          rsolver%isolverSpec = ior(rsolver%isolverSpec,&
                                    SV_SSPEC_STRUCTURENEEDSUPDATE)
        end if

        ! Check minimum multigrid level
        if (p_solverMultigrid%nlmin .ne. nlmax) then
          p_solverMultigrid%nlmin = nlmax
          rsolver%isolverSpec = ior(rsolver%isolverSpec,&
                                    SV_SSPEC_STRUCTURENEEDSUPDATE)
        end if
      else

        ! Check maximum multigrid level
        if (p_solverMultigrid%nlmax .ne. nlmax) then
          p_solverMultigrid%nlmax = nlmax
          rsolver%isolverSpec = ior(rsolver%isolverSpec,&
                                    SV_SSPEC_STRUCTURENEEDSUPDATE)
        end if

        ! Check minimum multigrid level
        if (p_solverMultigrid%nlmin .ne. min(nlmin, p_solverMultigrid%nlmin, nlmax)) then
          p_solverMultigrid%nlmin = min(nlmin, p_solverMultigrid%nlmin, nlmax)
          rsolver%isolverSpec = ior(rsolver%isolverSpec,&
                                    SV_SSPEC_STRUCTURENEEDSUPDATE)
        end if
      end if

      ! Do we have a coarse grid solver?
      if (associated(p_solverMultigrid%p_solverCoarsegrid)) then
        call solver_adjustHierarchy(p_solverMultigrid%p_solverCoarsegrid, nlmin, nlmax)
      end if

      ! Do we have a smoother?
      if (associated(p_solverMultigrid%p_smoother)) then
        ilbound = lbound(p_solverMultigrid%p_smoother, 1)
        iubound = ubound(p_solverMultigrid%p_smoother, 1)

        ! Do we have to adjust smoother?
        if ((p_solverMultigrid%nlmin+1 .ne. ilbound) .or.&
            (p_solverMultigrid%nlmax   .ne. iubound)) then

          ! Create template solver from smoother
          call solver_createSolver(rsolverTemplateSmoother,&
                                   p_solverMultigrid%p_smoother(iubound))

          ! Release all existing smoothers
          do i = ilbound, iubound
            call solver_releaseSolver(p_solverMultigrid%p_smoother(i))
          end do
          deallocate(p_solverMultigrid%p_smoother)

          ! Create new array of smoothers from template solver
          if (p_solverMultigrid%nlmin .lt. p_solverMultigrid%nlmax) then
            allocate(p_solverMultigrid%p_smoother(&
                p_solverMultigrid%nlmin:p_solverMultigrid%nlmax))
            do i = p_solverMultigrid%nlmin+1, p_solverMultigrid%nlmax
              call solver_createSolver(p_solverMultigrid%p_smoother(i),&
                                       rsolverTemplateSmoother)
            end do
          end if

          ! Release template smoother
          call solver_releaseSolver(rsolverTemplateSmoother)
          rsolver%isolverSpec = ior(rsolver%isolverSpec,&
                                    SV_SSPEC_STRUCTURENEEDSUPDATE)
        end if

        ! Finally, adjust hierarchy in smoothers, whereby the maximum
        ! multigrid level cannot exceed the level of the smoother
        do i = p_solverMultigrid%nlmin+1, p_solverMultigrid%nlmax
          call solver_adjustHierarchy(p_solverMultigrid%p_smoother(i),&
                                      p_solverMultigrid%nlmin, min(i, p_solverMultigrid%nlmax))
        end do
      end if
    end if

    ! Do we have a single-grid solver which has a preconditioner?

    ! BiCGSTAB ?
    if (associated(rsolver%p_solverBiCGSTAB)) then
      if (associated(rsolver%p_solverBiCGSTAB%p_precond)) then
        call solver_adjustHierarchy(rsolver%p_solverBiCGSTAB%p_precond)
      end if
    end if

    ! GMRES ?
    if (associated(rsolver%p_solverGMRES)) then
      if (associated(rsolver%p_solverGMRES%p_precond)) then
        call solver_adjustHierarchy(rsolver%p_solverGMRES%p_precond)
      end if
    end if

  end subroutine solver_adjustHierarchy

  ! ***************************************************************************

!<subroutine>

  subroutine solver_showHierarchy(rsolver)

!<description>
    ! This subroutine shows the solver hierarchy in human-readable format
!</description>

!<input>
    ! solver of the top-level
    type(t_solver), intent(in) :: rsolver
!</input>
!</subroutine>

    call output_separator(OU_SEP_PLUS)
    call output_line('Solver hierarchy for top-level solver: '//trim(rsolver%ssolverName))
    call output_separator(OU_SEP_PLUS)
    call showLevel(rsolver, '  ')

  contains

    ! Here, the real working routine follows

    !*************************************************************
    ! This subroutine shows the current level of th solver hierarchy
    ! and proceeds to the next sublevel

    recursive subroutine showLevel(rsolver,cindent)
      type(t_solver),   intent(in) :: rsolver
      character(LEN=*), intent(in) :: cindent

      integer :: i

      select case(rsolver%csolverType)
      case (SV_FMG)
        call startIndent(cindent, "Full Multigrid-Solver")
        if (associated(rsolver%p_solverMultigrid)) then
          call continueIndent(cindent, "NLMIN: "//&
              trim(adjustl(sys_si(rsolver%p_solverMultigrid%nlmin,3)))//" ("//&
              trim(adjustl(sys_si(rsolver%p_solverMultigrid%initialNlmin,3)))//")")
          call continueIndent(cindent, "NLMAX: "//&
              trim(adjustl(sys_si(rsolver%p_solverMultigrid%nlmax,3)))//" ("//&
              trim(adjustl(sys_si(rsolver%p_solverMultigrid%initialNlmax,3)))//")")

          if (associated(rsolver%p_solverMultigrid%p_smoother)) then
            do i = lbound(rsolver%p_solverMultigrid%p_smoother,1),&
                   ubound(rsolver%p_solverMultigrid%p_smoother,1)

              call continueIndent(cindent, "Smoother on level "//&
                  trim(adjustl(sys_si(i,3))))
              call showLevel(rsolver%p_solverMultigrid%p_smoother(i), cindent//"  |")
            end do
          end if

          if (associated(rsolver%p_solverMultigrid%p_solverCoarsegrid)) then
            call continueIndent(cindent, "Coarsegrid-Solver")
            call showLevel(rsolver%p_solverMultigrid%p_solverCoarsegrid,cindent//"  |")
          end if
        end if
        call stopIndent(cindent)

      case (SV_NONLINEARMG)
        call startIndent(cindent, "Nonlinear Multigrid-Solver")
        if (associated(rsolver%p_solverMultigrid)) then
          call continueIndent(cindent, "NLMIN: "//&
              trim(adjustl(sys_si(rsolver%p_solverMultigrid%nlmin,3)))//" ("//&
              trim(adjustl(sys_si(rsolver%p_solverMultigrid%initialNlmin,3)))//")")
          call continueIndent(cindent, "NLMAX: "//&
              trim(adjustl(sys_si(rsolver%p_solverMultigrid%nlmax,3)))//" ("//&
              trim(adjustl(sys_si(rsolver%p_solverMultigrid%initialNlmax,3)))//")")

          if (associated(rsolver%p_solverMultigrid%p_smoother)) then
            do i = lbound(rsolver%p_solverMultigrid%p_smoother,1),&
                   ubound(rsolver%p_solverMultigrid%p_smoother,1)

              call continueIndent(cindent, "Smoother on level "//&
                  trim(adjustl(sys_si(i,3))))
              call showLevel(rsolver%p_solverMultigrid%p_smoother(i), cindent//"  |")
            end do
          end if

          if (associated(rsolver%p_solverMultigrid%p_solverCoarsegrid)) then
            call continueIndent(cindent, "Coarsegrid-Solver")
            call showLevel(rsolver%p_solverMultigrid%p_solverCoarsegrid,cindent//"  |")
          end if
        end if
        call stopIndent(cindent)

      case (SV_LINEARMG)
        call startIndent(cindent, "Linear Multigrid-Solver")
        if (associated(rsolver%p_solverMultigrid)) then
          call continueIndent(cindent, "NLMIN: "//&
              trim(adjustl(sys_si(rsolver%p_solverMultigrid%nlmin,3)))//" ("//&
              trim(adjustl(sys_si(rsolver%p_solverMultigrid%initialNlmin,3)))//")")
          call continueIndent(cindent, "NLMAX: "//&
              trim(adjustl(sys_si(rsolver%p_solverMultigrid%nlmax,3)))//" ("//&
              trim(adjustl(sys_si(rsolver%p_solverMultigrid%initialNlmax,3)))//")")

          if (associated(rsolver%p_solverMultigrid%p_smoother)) then
            do i = lbound(rsolver%p_solverMultigrid%p_smoother,1),&
                   ubound(rsolver%p_solverMultigrid%p_smoother,1)

              call continueIndent(cindent, "Smoother on level "//&
                  trim(adjustl(sys_si(i,3))))
              call showLevel(rsolver%p_solverMultigrid%p_smoother(i), cindent//"  |")
            end do
          end if

          if (associated(rsolver%p_solverMultigrid%p_solverCoarsegrid)) then
            call continueIndent(cindent, "Coarsegrid-Solver")
            call showLevel(rsolver%p_solverMultigrid%p_solverCoarsegrid,cindent//"  |")
          end if
        end if
        call stopIndent(cindent)

      case(SV_NONLINEAR)
        call startIndent(cindent, "Nonlinear Solver")
        if (associated(rsolver%p_solverSubnode)) then
          call continueIndent(cindent, "Subsolver")
          call showLevel(rsolver%p_solverSubnode(1),cindent//"  |")
        end if

        if (associated(rsolver%p_solverDefcor)) then
          call continueIndent(cindent, "Defect-correction Solver")
        end if

        if (associated(rsolver%p_solverNewton)) then
          call continueIndent(cindent, "Newton Solver")
        end if
        call stopIndent(cindent)

      case(SV_LINEAR)
        call startIndent(cindent, "Linear Solver")
        if (associated(rsolver%p_solverSubnode)) then
          call continueIndent(cindent, "Subsolver")
          call showLevel(rsolver%p_solverSubnode(1),cindent//"  |")
        end if

        if (associated(rsolver%p_solverUMFPACK)) then
          call continueIndent(cindent, "UMFPACK Solver")
        end if

        if (associated(rsolver%p_solverJacobi)) then
          call continueIndent(cindent, "Jacobi Solver")
        end if

        if (associated(rsolver%p_solverSSOR)) then
          call continueIndent(cindent, "SSOR Solver")
        end if

        if (associated(rsolver%p_solverBiCGSTAB)) then
          call continueIndent(cindent, "BiCGSTAB Solver")
          if (associated(rsolver%p_solverBiCGSTAB%p_precond)) then
            call continueIndent(cindent, "BiCGSTAB Preconditioner")
            call showLevel(rsolver%p_solverBiCGSTAB%p_precond,cindent//"  |")
          end if
        end if

        if (associated(rsolver%p_solverGMRES)) then
          call continueIndent(cindent, "GMRES Solver")
          if (associated(rsolver%p_solverGMRES%p_precond)) then
            call continueIndent(cindent, "GMRES Preconditioner")
            call showLevel(rsolver%p_solverGMRES%p_precond,cindent//"  |")
          end if
        end if

        if (associated(rsolver%p_solverILU)) then
          call continueIndent(cindent, "ILU Solver")
        end if
        call stopIndent(cindent)

      end select

    end subroutine showLevel

    !*************************************************************
    ! This subroutine starts new indentation, i.e. "+--"

    subroutine startIndent(cindent, cmsg)
      character(LEN=*), intent(in) :: cindent
      character(LEN=*), intent(in) :: cmsg

      call output_line(cindent//'+--'//cmsg)

    end subroutine startIndent

    !*************************************************************
    ! This subroutine continues indentation, i.e. "  |"

    subroutine continueIndent(cindent, cmsg)
      character(LEN=*), intent(in) :: cindent
      character(LEN=*), intent(in) :: cmsg

      call output_line(cindent//'  |'//cmsg)

    end subroutine continueIndent

    !*************************************************************
    ! This subroutine stopt indentation, i.e. "  +"

    subroutine stopIndent(cindent)
      character(LEN=*), intent(in) :: cindent

      call output_line(cindent//'  +')

    end subroutine stopIndent

  end subroutine solver_showHierarchy

  ! ***************************************************************************

!<subroutine>

  subroutine solver_statistics(rsolver, iiterations)

!<description>
    ! For a given solver structure, this subroutine computes solver
    ! statistics such as total number of iterations, the rate of
    ! convergence and the status. If the optional parameter
    ! rsolverSource is given, then its data is copied to rsolver
    ! and used to compute the solver statistics for rsolver.
!</description>

!<input>
    ! Number of iterations
    integer, intent(in) :: iiterations
!</input>

!<inputoutput>
    ! Solver structure
    type(t_solver), intent(inout) :: rsolver
!</inputoutput>
!</subroutine>

    ! Set number of iterations
    rsolver%iiterations = min(iiterations, rsolver%nmaxIterations)

    ! Update total number of iterations
    rsolver%niterations = rsolver%niterations + rsolver%iiterations


    ! Compute rate of convergence
    if (rsolver%dinitialDefect .gt. rsolver%dfinalDefect) then
      rsolver%dconvergenceRate = (rsolver%dfinalDefect/&
                                  rsolver%dinitialDefect)**(1._DP/iiterations)
    else
      rsolver%dconvergenceRate = 1._DP
    end if


    ! Determine status of solver
    if (rsolver%istatus .eq. SV_STAGNATED) return

    if ((rsolver%dconvergenceRate .ge. 1._DP) .or.&
        (iiterations .ge. rsolver%nmaxIterations)) then
      rsolver%istatus = SV_DIVERGED
    else
      rsolver%istatus = SV_CONVERGED
    end if

  end subroutine solver_statistics

  ! ***************************************************************************

!<function>

  pure function solver_getStatus(rsolver) result(cstat)

!<description>
    ! This function returns a character string which represents the status
!</description>

!<input>
    ! Solver structure
    type(t_solver), intent(in) :: rsolver
!</input>

!<result>
    ! Status of the solver in human readable format
    character(LEN=SYS_STRLEN) :: cstat
!</result>
!</function>

    select case(rsolver%istatus)
    case(SV_CONVERGED)
      cstat = "converged"

    case(SV_ZERO_RHS)
      cstat = "failed [ zero right-hand side ]"

    case(SV_ZERO_DEF)
      cstat = "failed [ zero initial defect ]"

    case(SV_INF_DEF)
      cstat = "failed [ initial defect too large ]"

    case(SV_INCR_DEF)
      cstat = "failed [ defect increased ] "

    case(SV_DIVERGED)
      cstat = "diverged"

    case(SV_STAGNATED)
      cstat = "failed [ stagnation ]"

    case DEFAULT
      cstat = "undefined status"
    end select

  end function solver_getStatus

  ! ***************************************************************************

!<function>

  function solver_getNextSolver(rsolver, isubsolver) result(p_rsubsolver)

!<description>
    ! This functions sets the pointer to the next solver structure in the solver
    ! hierarchy. If the type of the given solver denotes a single grid solver,
    ! then the corresponding solver subnode is returned if available.
    ! For multigrid-type solvers the corresponding coarsegrid solver is
    ! returned. If there is no subsolver, then NULL ist returned.
!</description>

!<input>
    ! Solver that is used as base solver
    type(t_solver), intent(in) :: rsolver

    ! OPTIONAL: Number of the subsolver if there are more than one
    integer, intent(in), optional :: isubsolver
!</input>

!<result>
    ! Pointer to the subsolver
    ! If no subseolver exists, then NULL is returned
    type(t_solver), pointer :: p_rsubsolver
!</result>
!</function>

    ! local variables
    integer :: isub

    ! Initialize subsolver number
    isub = 1
    if (present(isubsolver)) isub = isubsolver

    ! What type of solver are we?
    select case(rsolver%csolverType)

    case (SV_COUPLED)
      if (associated(rsolver%p_solverSubnode)) then
        if (size(rsolver%p_solverSubnode) .ge. isub) then
          p_rsubsolver => rsolver%p_solverSubnode(isub)
        else
          nullify(p_rsubsolver)
        end if
      else
        nullify(p_rsubsolver)
      end if

    case (SV_FMG, SV_NONLINEARMG, SV_LINEARMG)
      if (associated(rsolver%p_solverMultigrid)) then
        if (associated(rsolver%p_solverMultigrid%p_solverCoarsegrid)) then
          p_rsubsolver => rsolver%p_solverMultigrid%p_solverCoarsegrid
        else
          nullify(p_rsubsolver)
        end if
      else
        nullify(p_rsubsolver)
      end if

    case (SV_NONLINEAR, SV_LINEAR)
      if (associated(rsolver%p_solverSubnode)) then
        p_rsubsolver => rsolver%p_solverSubnode(1)
      else
        nullify(p_rsubsolver)
      end if

    case DEFAULT
      nullify(p_rsubsolver)
    end select

  end function solver_getNextSolver

  ! ***************************************************************************

!<function>

  function solver_getNextSolverByTypes(rsolver, CsolverType, isubsolver)&
      result(p_rsubsolver)

!<description>
    ! This functions sets the pointer to the next solver structure in the
    ! solver hierarchy which satisfies the specified type isolver. If
    ! there is no subsolver which satisfies the prescribed type, then
    ! NULL ist returned.
!</description>

!<input>
    ! Solver that is used as base solver
    type(t_solver), intent(in), target :: rsolver

    ! Array of types of solver to return
    integer, dimension(:), intent(in) :: CsolverType

    ! OPTIONAL: Number of the subsolver if there are more than one
    integer, intent(in), optional :: isubsolver
!</input>

!<result>
    ! Pointer to the subsolver
    ! If no subseolver exists, then NULL is returned
    type(t_solver), pointer :: p_rsubsolver
!</result>
!</function>

    p_rsubsolver => rsolver
    do while (associated(p_rsubsolver))
      if (any(p_rsubsolver%csolverType .eq. CsolverType)) exit
      p_rsubsolver => solver_getNextSolver(p_rsubsolver, isubsolver)
    end do

  end function solver_getNextSolverByTypes

  ! ***************************************************************************

!<function>

  function solver_getNextSolverByType(rsolver, csolverType, isubsolver)&
      result(p_rsubsolver)

!<description>
    ! This functions sets the pointer to the next solver structure in the
    ! solver hierarchy which satisfies the specified type isolver. If
    ! there is no subsolver which satisfies the prescribed type, then
    ! NULL ist returned.
!</description>

!<input>
    ! Solver that is used as base solver
    type(t_solver), intent(in), target :: rsolver

    ! Type of solver to return
    integer, intent(in) :: csolverType

    ! OPTIONAL: Number of the subsolver if there are more than one
    integer, intent(in), optional :: isubsolver
!</input>

!<result>
    ! Pointer to the subsolver
    ! If no subseolver exists, then NULL is returned
    type(t_solver), pointer :: p_rsubsolver
!</result>
!</function>

    p_rsubsolver => rsolver
    do while (associated(p_rsubsolver))
      if (p_rsubsolver%csolverType .eq. csolverType) exit
      p_rsubsolver => solver_getNextSolver(p_rsubsolver, isubsolver)
    end do

  end function solver_getNextSolverByType

  ! ***************************************************************************

!<function>

  function solver_getMinimumMultigridlevel(rsolver, nlminAlt) result(nlmin)

!<description>
    ! This function returns the minimal multigrid level required throughout the
    ! entire solver structure. Note that preconditioner and smoother structures
    ! are not taken into account as they are only auxiliary solvers.
    !
    ! If the optional parameter nlminAlt is given, then this value is adopted
    ! if no valid nlmin could be found in the solver structure
!</description>

!<input>
    ! solver
    type(t_solver), intent(in) :: rsolver

    ! OPTIONAL: alternative nlmin which is used, if no nlmin could be found
    integer, intent(in), optional :: nlminAlt
!</input>

!<result>
    ! minimum multigrid level
    integer :: nlmin
!</result>
!</function>

    ! Determine NLMIN
    nlmin = get_nlmin(rsolver)

    ! Adjust NLMIN if required
    if (nlmin .eq. huge(1)) then
      if (present(nlminAlt)) then
        nlmin = nlminAlt
      else
        nlmin = 1
      end if
    end if

  contains

    ! Here, the real working routine follows

    !*************************************************************
    ! This function returns the minimum mulitigrid level

    recursive function get_nlmin(rsolver) result(nlmin)
      type(t_solver), intent(in) :: rsolver
      integer :: nlmin

      ! local varialbles
      integer :: i

      ! Initialisation
      nlmin = huge(1)

      ! Do we have a generic solver subnode?
      if (associated(rsolver%p_solverSubnode)) then
        do i = 1, size(rsolver%p_solverSubnode)
          nlmin = min(nlmin, get_nlmin(rsolver%p_solverSubnode(i)))
        end do
      end if

      ! Do we have multigrid subnode?
      if (associated(rsolver%p_solverMultigrid)) then
        nlmin = min(nlmin, rsolver%p_solverMultigrid%nlmin)

        ! Do we have coarsegrid subnode?
        if (associated(rsolver%p_solverMultigrid%p_solverCoarsegrid)) then
          nlmin = min(nlmin, get_nlmin(rsolver%p_solverMultigrid%p_solverCoarsegrid))
        end if
      end if
    end function get_nlmin

  end function solver_getMinimumMultigridlevel

  ! ***************************************************************************

!<function>

  function solver_getMaximumMultigridlevel(rsolver, nlmaxAlt) result(nlmax)

!<description>
    ! This function returns the maximal multigrid level required throughout the
    ! entire solver structure. Note that preconditioner and smoother structures
    ! are not taken into account as they are only auxiliary solvers.
    !
    ! If the optional parameter nlmaxAlt is given, then this value is adopted
    ! if no valid nlmax could be found in the solver structure
!</description>

!<input>
    ! single-grid structure
    type(t_solver), intent(in) :: rsolver

    ! OPTIONAL: alternative nlmax which is used, if no nlmax could be found
    integer, intent(in), optional :: nlmaxAlt
!</input>

!<result>
    ! maximal multigrid level
    integer :: nlmax
!</result>
!</function>

    ! Determine NLMAX
    nlmax = get_nlmax(rsolver)

    ! Adjust NLMAX if required
    if (nlmax .eq. 0) then
      if (present(nlmaxAlt)) then
        nlmax = nlmaxAlt
      else
        nlmax = 1
      end if
    end if

  contains

    ! Here, the real working routine follows

    !*************************************************************
    ! This function returns the maximum mulitigrid level

    recursive function get_nlmax(rsolver) result(nlmax)
      type(t_solver), intent(in) :: rsolver
      integer :: nlmax

      ! local variables
      integer :: i

      ! Initialisation
      nlmax = 0

      ! Do we have a generic solver subnode?
      if (associated(rsolver%p_solverSubnode)) then
        do i = 1, size(rsolver%p_solverSubnode)
          nlmax = max(nlmax, get_nlmax(rsolver%p_solverSubnode(i)))
        end do
      end if

      ! Do we have multigrid subnode?
      if (associated(rsolver%p_solverMultigrid)) then
        nlmax = max(nlmax, rsolver%p_solverMultigrid%nlmax)

        ! Do we have coarsegrid subnode?
        if (associated(rsolver%p_solverMultigrid%p_solverCoarsegrid)) then
          nlmax = max(nlmax, get_nlmax(rsolver%p_solverMultigrid%p_solverCoarsegrid))
        end if
      end if
    end function get_nlmax

  end function solver_getMaximumMultigridlevel

  ! ***************************************************************************

!<subroutine>

  subroutine solver_setMinimumMultigridlevel(rsolver, nlmin)

!<description>
    ! This subroutine sets the minimum multigrid level
    ! required throughout the entire solver structure
!</description>

!<input>
    ! minimal multigrid level
    integer, intent(in) :: nlmin
!</input>

!<inputoutput>
    ! solver
    type(t_solver), intent(inout) :: rsolver
!</inputoutput>
!</subroutine>

    ! Set minimum multigrid level
    call set_nlmin(rsolver, nlmin)

  contains

    ! Here, the real working routine follows

    !*************************************************************
    ! This subroutine sets the minimum multigrid leve recursively

    recursive subroutine set_nlmin(rsolver, nlmin)
      type(t_solver), intent(inout) :: rsolver
      integer, intent(in) :: nlmin

      ! local variables
      integer :: i

      ! Do we have a generic solver subnode?
      if (associated(rsolver%p_solverSubnode)) then
        do i = 1, size(rsolver%p_solverSubnode)
          call set_nlmin(rsolver%p_solverSubnode(i), nlmin)
        end do
      end if

      ! Do we have multigrid subnode?
      if (associated(rsolver%p_solverMultigrid)) then
        rsolver%p_solverMultigrid%nlmin = nlmin

        ! Do we have coarsegrid subnode?
        if (associated(rsolver%p_solverMultigrid%p_solverCoarsegrid)) then
          call set_nlmin(rsolver%p_solverMultigrid%p_solverCoarsegrid, nlmin)
        end if

        ! Do we have a smoother subnode?
        if (associated(rsolver%p_solverMultigrid%p_smoother)) then
          do i = lbound(rsolver%p_solverMultigrid%p_smoother,1),&
                 ubound(rsolver%p_solverMultigrid%p_smoother,1)
            call set_nlmin(rsolver%p_solverMultigrid%p_smoother(i), nlmin)
          end do
        end if
      end if

      ! Do we have a single-grid solver which has a preconditioner?

      ! BiCGSTAB ?
      if (associated(rsolver%p_solverBiCGSTAB)) then
        if (associated(rsolver%p_solverBiCGSTAB%p_precond)) then
          call set_nlmin(rsolver%p_solverBiCGSTAB%p_precond, nlmin)
        end if
      end if

      ! GMRES ?
      if (associated(rsolver%p_solverGMRES)) then
        if (associated(rsolver%p_solverGMRES%p_precond)) then
          call set_nlmin(rsolver%p_solverGMRES%p_precond, nlmin)
        end if
      end if
    end subroutine set_nlmin

  end subroutine solver_setMinimumMultigridlevel

  ! ***************************************************************************

!<subroutine>

  subroutine solver_setMaximumMultigridlevel(rsolver, nlmax)

!<description>
    ! This subroutine sets the maximum multigrid level
    ! required throughout the entire solver structure
!</description>

!<input>
    ! maximal multigrid level
    integer, intent(in) :: nlmax
!</input>

!<inputoutput>
    ! solver
    type(t_solver), intent(inout) :: rsolver
!</inputoutput>
!</subroutine>

    ! Set maximum multigrid level
    call set_nlmax(rsolver, nlmax)

  contains

    ! Here, the real working routine follows

    !*************************************************************
    ! This subroutine sets the maximum multigrid leve recursively

    recursive subroutine set_nlmax(rsolver, nlmax)
      type(t_solver), intent(inout) :: rsolver
      integer, intent(in) :: nlmax

      ! local variables
      integer :: i

      ! Do we have a generic solver subnode?
      if (associated(rsolver%p_solverSubnode)) then
        do i = 1, size(rsolver%p_solverSubnode)
          call set_nlmax(rsolver%p_solverSubnode(i), nlmax)
        end do
      end if

      ! Do we have multigrid subnode?
      if (associated(rsolver%p_solverMultigrid)) then
        rsolver%p_solverMultigrid%nlmax = nlmax

        ! Do we have coarsegrid subnode?
        if (associated(rsolver%p_solverMultigrid%p_solverCoarsegrid)) then
          call set_nlmax(rsolver%p_solverMultigrid%p_solverCoarsegrid, nlmax)
        end if

        ! Do we have a smoother subnode?
        if (associated(rsolver%p_solverMultigrid%p_smoother)) then
          do i = lbound(rsolver%p_solverMultigrid%p_smoother,1),&
                 ubound(rsolver%p_solverMultigrid%p_smoother,1)
            call set_nlmax(rsolver%p_solverMultigrid%p_smoother(i), nlmax)
          end do
        end if
      end if

      ! Do we have a single-grid solver which has a preconditioner?

      ! BiCGSTAB ?
      if (associated(rsolver%p_solverBiCGSTAB)) then
        if (associated(rsolver%p_solverBiCGSTAB%p_precond)) then
          call set_nlmax(rsolver%p_solverBiCGSTAB%p_precond, nlmax)
        end if
      end if

      ! GMRES ?
      if (associated(rsolver%p_solverGMRES)) then
        if (associated(rsolver%p_solverGMRES%p_precond)) then
          call set_nlmax(rsolver%p_solverGMRES%p_precond, nlmax)
        end if
      end if
    end subroutine set_nlmax

  end subroutine solver_setMaximumMultigridlevel

  ! ***************************************************************************

!<subroutine>

  recursive subroutine solver_setBoundaryCondition(rsolver, rboundaryCondition, brecursive)

!<description>
    ! This subroutine associates a given boundary condition to the solver.
    ! If brecursive=.TRUE. then the same boundary condition is used for
    ! all solvers in the underlying solver hierarchy. Otherwise, the
    ! boundary conditions is only associated to the top-level solver.
!</description>

!<input>
    ! boundary condition
    type(t_boundaryCondition), intent(in), target :: rboundaryCondition

    ! Flag: If TRUE, then boundary conditions are set recursively for all
    !       subsolvers in the solver hierarchy.
    logical, intent(in) :: brecursive
!</input>

!<inputoutput>
    ! Solver structure
    type(t_solver), intent(inout) :: rsolver
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: i

    ! Set pointer to boundary condition
    rsolver%rboundaryCondition => rboundaryCondition

    ! Set boundary conditions recursively?
    if (.not.brecursive) return

    ! Do we have a generic solver subnode?
    if (associated(rsolver%p_solverSubnode)) then
      do i = 1, size(rsolver%p_solverSubnode)
        call solver_setBoundaryCondition(rsolver%p_solverSubnode(i),&
            rboundaryCondition, brecursive)
      end do
    end if

    ! Do we have a multigrid subnode?
    if (associated(rsolver%p_solverMultigrid)) then

      ! Do we have a coarsegrid solver?
      if (associated(rsolver%p_solverMultigrid%p_solverCoarsegrid)) then
        call solver_setBoundaryCondition(&
            rsolver%p_solverMultigrid%p_solverCoarsegrid, rboundaryCondition, brecursive)
      end if

      ! Do we have a smoother?
      if (associated(rsolver%p_solverMultigrid%p_smoother)) then
        do i = lbound(rsolver%p_solverMultigrid%p_smoother,1),&
               ubound(rsolver%p_solverMultigrid%p_smoother,1)
          call solver_setBoundaryCondition(&
              rsolver%p_solverMultigrid%p_smoother(i), rboundaryCondition, brecursive)
        end do
      end if
    end if

    ! Do we have a single-grid solver which has a preconditioner

    ! BiCGSTAB ?
    if (associated(rsolver%p_solverBiCGSTAB)) then
      if (associated(rsolver%p_solverBiCGSTAB%p_precond)) then
        call solver_setBoundaryCondition(&
            rsolver%p_solverBiCGSTAB%p_precond, rboundaryCondition, brecursive)
      end if
    end if

    ! GMRES ?
    if (associated(rsolver%p_solverGMRES)) then
      if (associated(rsolver%p_solverGMRES%p_precond)) then
        call solver_setBoundaryCondition(&
            rsolver%p_solverGMRES%p_precond, rboundaryCondition, brecursive)
      end if
    end if

  end subroutine solver_setBoundaryCondition

  ! ***************************************************************************

!<subroutine>

  subroutine solver_setSolverMatrixScalar(rsolver, rmatrix, ilev)

!<description>
    ! This subroutine associates a scalar matrix to the solver.
!</description>

!<input>
    ! scalar matrix
    type(t_matrixScalar), intent(in) :: rmatrix

    ! OPTIONAL: level of the multigrid solver
    integer, intent(in), optional :: ilev
!</input>

!<inputoutput>
    ! solver structure
    type(t_solver), intent(inout) :: rsolver
!</inputoutput>
!</subroutine>

    ! local variable
    type(t_solverMultigrid), pointer :: p_solverMultigrid


    ! What kind of solver are we?
    select case(rsolver%csolverType)
    case (SV_FMG, SV_NONLINEARMG, SV_NONLINEAR)
      ! ------------------------------------------------------------------------
      ! There are no matrices for nonlinear solvers


    case (SV_LINEARMG)
      ! ------------------------------------------------------------------------
      ! Set system matrices for linear multigrid solver
      p_solverMultigrid => rsolver%p_solverMultigrid

      ! Check if multigrid level is specified
      if (.not.present(ilev)) then
        call output_line('Multigrid level is missing!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'solver_setSolverMatrixScalar')
        call sys_halt()
      end if

      ! Check if subarray of matrices exists
      if (associated(p_solverMultigrid%Rmatrix)) then

        ! Check if subarray of matrices has correct size
        if ((lbound(p_solverMultigrid%Rmatrix,1) .le. ilev) .and.&
            (ubound(p_solverMultigrid%Rmatrix,1) .ge. ilev)) then

          ! Release existing matrix and set new matrix
          call lsysbl_releaseMatrix(p_solverMultigrid%Rmatrix(ilev))
          call lsysbl_createMatFromScalar(rmatrix,&
                                          p_solverMultigrid%Rmatrix(ilev))

          ! Mark solver for update: structure + content
          rsolver%isolverSpec = ior(rsolver%isolverSpec,&
                                    SV_SSPEC_NEEDSUPDATE)
        else
          call output_line('Specified multigrid level does not exist!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'solver_setSolverMatrixScalar')
          call sys_halt()
        end if
      else
        call output_line('Multigrid solver subnode does not exist!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'solver_setSolverMatrixScalar')
        call sys_halt()
      end if


    case (SV_LINEAR)
      ! ------------------------------------------------------------------------

      ! UMFPACK subnode
      if (associated(rsolver%p_solverUMFPACK)) then

        ! Release existing matrix and set new matrix
        call lsysbl_releaseMatrix(rsolver%p_solverUMFPACK%rmatrix)
        call lsysbl_createMatFromScalar(rmatrix,&
                                        rsolver%p_solverUMFPACK%rmatrix)

        ! Mark solver for update: structure + content
        rsolver%isolverSpec = ior(rsolver%isolverSpec,&
                                  SV_SSPEC_NEEDSUPDATE)
      end if

      ! Jacobi subnode
      if (associated(rsolver%p_solverJacobi)) then

        ! Release existing matrix and set new matrix
        call lsysbl_releaseMatrix(rsolver%p_solverJacobi%rmatrix)
        call lsysbl_createMatFromScalar(rmatrix,&
                                        rsolver%p_solverJacobi%rmatrix)

        ! Mark solver for update: structure + content
        rsolver%isolverSpec = ior(rsolver%isolverSpec,&
                                  SV_SSPEC_NEEDSUPDATE)
      end if

      ! (S)SOR subnode
      if (associated(rsolver%p_solverSSOR)) then

        ! Release existing matrix and set new matrix
        call lsysbl_releaseMatrix(rsolver%p_solverSSOR%rmatrix)
        call lsysbl_createMatFromScalar(rmatrix,&
                                        rsolver%p_solverSSOR%rmatrix)

        ! Mark solver for update: structure + content
        rsolver%isolverSpec = ior(rsolver%isolverSpec,&
                                  SV_SSPEC_NEEDSUPDATE)
      end if

      ! BiCGSTAB subnode
      if (associated(rsolver%p_solverBiCGSTAB)) then

        ! Release existing matrix and set new matrix
        call lsysbl_releaseMatrix(rsolver%p_solverBiCGSTAB%rmatrix)
        call lsysbl_createMatFromScalar(rmatrix,&
                                        rsolver%p_solverBiCGSTAB%rmatrix)

        ! Mark solver for update: structure + content
        rsolver%isolverSpec = ior(rsolver%isolverSpec,&
                                  SV_SSPEC_NEEDSUPDATE)
      end if

      ! GMRES subnode
      if (associated(rsolver%p_solverGMRES)) then

        ! Release existing matrix and set new matrix
        call lsysbl_releaseMatrix(rsolver%p_solverGMRES%rmatrix)
        call lsysbl_createMatFromScalar(rmatrix,&
                                        rsolver%p_solverGMRES%rmatrix)

        ! Mark solver for update: structure + content
        rsolver%isolverSpec = ior(rsolver%isolverSpec,&
                                  SV_SSPEC_NEEDSUPDATE)
      end if

      ! ILU subnode
      if (associated(rsolver%p_solverILU)) then

        ! Release existing matrix and set new matrix
        call lsysbl_releaseMatrix(rsolver%p_solverILU%rmatrix)
        call lsysbl_createMatFromScalar(rmatrix,&
                                        rsolver%p_solverILU%rmatrix)

        ! Mark solver for update: structure + content
        rsolver%isolverSpec = ior(rsolver%isolverSpec,&
                                  SV_SSPEC_NEEDSUPDATE)
      end if

    case DEFAULT
      call output_line('Unsupported solver type!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'solver_setSolverMatrixScalar')
      call sys_halt()
    end select

  end subroutine solver_setSolverMatrixScalar

  ! ***************************************************************************

!<subroutine>

  subroutine solver_setSolverMatrixBlock(rsolver, rmatrix, ilev)

!<description>
    ! This subroutine associates a block matrix to the solver.
!</description>

!<input>
    ! block matrix
    type(t_matrixBlock), intent(in) :: rmatrix

    ! OPTIONAL: level of the multigrid solver
    integer, intent(in), optional :: ilev
!</input>

!<inputoutput>
    ! solver structure
    type(t_solver), intent(inout) :: rsolver
!</inputoutput>
!</subroutine>

    ! local variable
    type(t_solverMultigrid), pointer :: p_solverMultigrid


    ! What kind of solver are we?
    select case(rsolver%csolverType)
    case (SV_FMG, SV_NONLINEARMG, SV_NONLINEAR)
      ! ------------------------------------------------------------------------
      ! There are no matrices for nonlinear solvers


    case (SV_LINEARMG)
      ! ------------------------------------------------------------------------
      ! Set system matrices for linear multigrid solver
      p_solverMultigrid => rsolver%p_solverMultigrid

      ! Check if multigrid level is specified
      if (.not.present(ilev)) then
        call output_line('Multigrid level is missing!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'solver_setSolverMatrixBlock')
        call sys_halt()
      end if

      ! Check if subarray of matrices exists
      if (associated(p_solverMultigrid%Rmatrix)) then

        ! Check if subarray of matrices has correct size
        if ((lbound(p_solverMultigrid%Rmatrix,1) .le. ilev) .and.&
            (ubound(p_solverMultigrid%Rmatrix,1) .ge. ilev)) then

          ! Release existing matrix and set new matrix
          call lsysbl_releaseMatrix(p_solverMultigrid%Rmatrix(ilev))
          call lsysbl_duplicateMatrix(rmatrix,&
                                      p_solverMultigrid%rmatrix(ilev),&
                                      LSYSSC_DUP_SHARE, LSYSSC_DUP_SHARE)

          ! Mark solver for update: structure + content
          rsolver%isolverSpec = ior(rsolver%isolverSpec,&
                                    SV_SSPEC_NEEDSUPDATE)
        else
          call output_line('Specified multigrid level does not exist!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'solver_setSolverMatrixBlock')
          call sys_halt()
        end if
      else
        call output_line('Multigrid solver subnode does not exist!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'solver_setSolverMatrixBlock')
        call sys_halt()
      end if


    case (SV_LINEAR)
      ! ------------------------------------------------------------------------

      ! UMFPACK subnode
      if (associated(rsolver%p_solverUMFPACK)) then

        ! Release existing matrix and set new matrix
        call lsysbl_releaseMatrix(rsolver%p_solverUMFPACK%rmatrix)
        call lsysbl_duplicateMatrix(rmatrix,&
                                    rsolver%p_solverUMFPACK%rmatrix,&
                                    LSYSSC_DUP_SHARE, LSYSSC_DUP_SHARE)

        ! Mark solver for update: structure + content
        rsolver%isolverSpec = ior(rsolver%isolverSpec,&
                                  SV_SSPEC_NEEDSUPDATE)
      end if

      ! Jacobi subnode
      if (associated(rsolver%p_solverJacobi)) then

        ! Release existing matrix and set new matrix
        call lsysbl_releaseMatrix(rsolver%p_solverJacobi%rmatrix)
        call lsysbl_duplicateMatrix(rmatrix,&
                                    rsolver%p_solverJacobi%rmatrix,&
                                    LSYSSC_DUP_SHARE, LSYSSC_DUP_SHARE)

        ! Mark solver for update: structure + content
        rsolver%isolverSpec = ior(rsolver%isolverSpec,&
                                  SV_SSPEC_NEEDSUPDATE)
      end if

      ! (S)SOR subnode
      if (associated(rsolver%p_solverSSOR)) then

        ! Release existing matrix and set new matrix
        call lsysbl_releaseMatrix(rsolver%p_solverSSOR%rmatrix)
        call lsysbl_duplicateMatrix(rmatrix,&
                                    rsolver%p_solverSSOR%rmatrix,&
                                    LSYSSC_DUP_SHARE, LSYSSC_DUP_SHARE)

        ! Mark solver for update: structure + content
        rsolver%isolverSpec = ior(rsolver%isolverSpec,&
                                  SV_SSPEC_NEEDSUPDATE)
      end if

      ! BiCGSTAB subnode
      if (associated(rsolver%p_solverBiCGSTAB)) then

        ! Release existing matrix and set new matrix
        call lsysbl_releaseMatrix(rsolver%p_solverBiCGSTAB%rmatrix)
        call lsysbl_duplicateMatrix(rmatrix,&
                                    rsolver%p_solverBiCGSTAB%rmatrix,&
                                    LSYSSC_DUP_SHARE, LSYSSC_DUP_SHARE)

        ! Mark solver for update: structure + content
        rsolver%isolverSpec = ior(rsolver%isolverSpec,&
                                  SV_SSPEC_NEEDSUPDATE)
      end if

      ! GMRES subnode
      if (associated(rsolver%p_solverGMRES)) then

        ! Release existing matrix and set new matrix
        call lsysbl_releaseMatrix(rsolver%p_solverGMRES%rmatrix)
        call lsysbl_duplicateMatrix(rmatrix,&
                                    rsolver%p_solverGMRES%rmatrix,&
                                    LSYSSC_DUP_SHARE, LSYSSC_DUP_SHARE)

        ! Mark solver for update: structure + content
        rsolver%isolverSpec = ior(rsolver%isolverSpec,&
                                  SV_SSPEC_NEEDSUPDATE)
      end if

      ! ILU subnode
      if (associated(rsolver%p_solverILU)) then

        ! Release existing matrix and set new matrix
        call lsysbl_releaseMatrix(rsolver%p_solverILU%rmatrix)
        call lsysbl_duplicateMatrix(rmatrix,&
                                    rsolver%p_solverILU%rmatrix,&
                                    LSYSSC_DUP_SHARE, LSYSSC_DUP_SHARE)

        ! Mark solver for update: structure + content
        rsolver%isolverSpec = ior(rsolver%isolverSpec,&
                                  SV_SSPEC_NEEDSUPDATE)
      end if

    case DEFAULT
      call output_line('Unsupported solver type!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'solver_setSolverMatrixBlock')
      call sys_halt()
    end select

  end subroutine solver_setSolverMatrixBlock

  ! ***************************************************************************

!<subroutine>

  subroutine solver_setPrecondMatrixScalar(rsolver, rmatrix, ilev)

!<description>
    ! This subroutine associates a scalar matrix to the preconditioner of the solver.
!</description>

!<input>
    ! scalar matrix
    type(t_matrixScalar), intent(in) :: rmatrix

    ! OPTIONAL: level of the multigrid solver
    integer, intent(in), optional :: ilev
!</input>

!<inputoutput>
    ! solver structure
    type(t_solver), intent(inout) :: rsolver
!</inputoutput>
!</subroutine>

    ! What kind of solver are we?
    select case(rsolver%csolverType)
    case (SV_FMG, SV_NONLINEARMG, SV_LINEARMG, SV_NONLINEAR)
      ! ------------------------------------------------------------------------
      ! No preconditioner required


    case (SV_LINEAR)
      ! ------------------------------------------------------------------------
      ! BiCGSTAB solver
      if (associated(rsolver%p_solverBiCGSTAB)) then
        if (associated(rsolver%p_solverBiCGSTAB%p_precond)) then
          call solver_setSolverMatrix(rsolver%p_solverBiCGSTAB%p_precond,&
                                      rmatrix, ilev)
        end if
      end if

      ! GMRES solver
      if (associated(rsolver%p_solverGMRES)) then
        if (associated(rsolver%p_solverGMRES%p_precond)) then
          call solver_setSolverMatrix(rsolver%p_solverGMRES%p_precond,&
                                      rmatrix, ilev)
        end if
      end if


    case DEFAULT
      call output_line('Unsupported solver type!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'solver_setPrecondMatrixScalar')
      call sys_halt()
    end select
  end subroutine solver_setPrecondMatrixScalar


  ! ***************************************************************************

!<subroutine>

  subroutine solver_setPrecondMatrixBlock(rsolver, rmatrix, ilev)

!<description>
    ! This subroutine associates a block matrix to the preconditioner of the solver.
!</description>

!<input>
    ! block matrix
    type(t_matrixBlock), intent(in) :: rmatrix

    ! OPTIONAL: level of the multigrid solver
    integer, intent(in), optional :: ilev
!</input>

!<inputoutput>
    ! solver structure
    type(t_solver), intent(inout) :: rsolver
!</inputoutput>
!</subroutine>

    ! What kind of solver are we?
    select case(rsolver%csolverType)
    case (SV_FMG, SV_NONLINEARMG, SV_LINEARMG, SV_NONLINEAR)
      ! ------------------------------------------------------------------------
      ! No preconditioner required


    case (SV_LINEAR)
      ! ------------------------------------------------------------------------
      ! BiCGSTAB solver
      if (associated(rsolver%p_solverBiCGSTAB)) then
        if (associated(rsolver%p_solverBiCGSTAB%p_precond)) then
          call solver_setSolverMatrix(rsolver%p_solverBiCGSTAB%p_precond,&
                                      rmatrix, ilev)
        end if
      end if

      ! GMRES solver
      if (associated(rsolver%p_solverGMRES)) then
        if (associated(rsolver%p_solverGMRES%p_precond)) then
          call solver_setSolverMatrix(rsolver%p_solverGMRES%p_precond,&
                                      rmatrix, ilev)
        end if
      end if


    case DEFAULT
      call output_line('Unsupported solver type!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'solver_setPrecondMatrixBlock')
      call sys_halt()
    end select
  end subroutine solver_setPrecondMatrixBlock

  ! ***************************************************************************

!<subroutine>

  subroutine solver_setSmootherMatrixScalar(rsolver, rmatrix, ilev)

!<description>
    ! This subroutine associates a scalar matrix to the smoother of the solver.
!</description>

!<input>
    ! scalar matrix
    type(t_matrixScalar), intent(in) :: rmatrix

    ! OPTIONAL: level of the multigrid solver
    integer, intent(in), optional :: ilev
!</input>

!<inputoutput>
    ! solver structure
    type(t_solver), intent(inout) :: rsolver
!</inputoutput>
!</subroutine>

    ! local variable
    type(t_solverMultigrid), pointer :: p_solverMultigrid
    type(t_solver), dimension(:), pointer :: p_smoother

    ! What kind of solver are we?
    select case(rsolver%csolverType)
    case (SV_FMG, SV_NONLINEARMG)
      ! ------------------------------------------------------------------------
      ! There are no matrices for nonlinear solvers


    case (SV_NONLINEAR, SV_LINEAR)
      ! ------------------------------------------------------------------------
      ! There are no smoothers for single-grid solvers


    case (SV_LINEARMG)
      ! ------------------------------------------------------------------------
      ! Linear multigrid solver
      if (associated(rsolver%p_solverMultigrid)) then

        ! Check if multigrid level is specified
        if (.not.present(ilev)) then
          call output_line('Multigrid level is missing!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'solver_setSmootherMatrixScalar')
          call sys_halt()
        end if

        ! Set pointer
        p_solverMultigrid => rsolver%p_solverMultigrid

        ! Check if subarray of smoothers exists
        if (associated(p_solverMultigrid%p_smoother)) then

          ! Set pointer
          p_smoother => p_solverMultigrid%p_smoother

          ! Check if subarray of smoothers has correct size
          if ((lbound(p_smoother,1) .le. ilev) .and.&
              (ubound(p_smoother,1) .ge. ilev)) then
            call solver_setSolverMatrix(p_smoother(ilev), rmatrix)
          else
            call output_line('Specified multigrid level does not exist!',&
                             OU_CLASS_ERROR,OU_MODE_STD,'solver_setSmootherMatrixScalar')
            call sys_halt()
          end if
        end if
      end if


    case DEFAULT
      call output_line('Unsupported solver type!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'solver_setSmootherMatrixScalar')
      call sys_halt()
    end select

  end subroutine solver_setSmootherMatrixScalar

  ! ***************************************************************************

!<subroutine>

  subroutine solver_setSmootherMatrixBlock(rsolver, rmatrix, ilev)

!<description>
    ! This subroutine associates a block matrix to the smoother of the solver.
!</description>

!<input>
    ! block matrix
    type(t_matrixBlock), intent(in) :: rmatrix

    ! OPTIONAL: level of the multigrid solver
    integer, intent(in), optional :: ilev
!</input>

!<inputoutput>
    ! solver structure
    type(t_solver), intent(inout) :: rsolver
!</inputoutput>
!</subroutine>

    ! local variable
    type(t_solverMultigrid), pointer :: p_solverMultigrid
    type(t_solver), dimension(:), pointer :: p_smoother

    ! What kind of solver are we?
    select case(rsolver%csolverType)
    case (SV_FMG, SV_NONLINEARMG)
      ! ------------------------------------------------------------------------
      ! There are no matrices for nonlinear solvers


    case (SV_NONLINEAR, SV_LINEAR)
      ! ------------------------------------------------------------------------
      ! There are no smoothers for single-grid solvers


    case (SV_LINEARMG)
      ! ------------------------------------------------------------------------
      ! Linear multigrid solver
      if (associated(rsolver%p_solverMultigrid)) then

        ! Check if multigrid level is specified
        if (.not.present(ilev)) then
          call output_line('Multigrid level is missing!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'solver_setSmootherMatrixBlock')
          call sys_halt()
        end if

        ! Set pointer
        p_solverMultigrid => rsolver%p_solverMultigrid

        ! Check if subarray of smoothers exists
        if (associated(p_solverMultigrid%p_smoother)) then

          ! Set pointer
          p_smoother => p_solverMultigrid%p_smoother

          ! Check if subarray of smoothers has correct size
          if ((lbound(p_smoother,1) .le. ilev) .and.&
              (ubound(p_smoother,1) .ge. ilev)) then
            call solver_setSolverMatrix(p_smoother(ilev), rmatrix)
          else
            call output_line('Specified multigrid level does not exist!',&
                             OU_CLASS_ERROR,OU_MODE_STD,'solver_setSmootherMatrixBlock')
            call sys_halt()
          end if
        end if
      end if


    case DEFAULT
      call output_line('Unsupported solver type!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'solver_setSmootherMatrixBlock')
      call sys_halt()
    end select

  end subroutine solver_setSmootherMatrixBlock

  ! ***************************************************************************

!<subroutine>

  recursive subroutine solver_updateStructure(rsolver)

!<description>
    ! This subroutine updates the internal structure of the solver.
    ! Remark: The task of this subroutine is twofold.
    !
    ! 1.) If the internal structures like temporal vectors, smoothers,
    !     etc. are not allocated, or if their shape does not match
    !     that of the encompassing structure, then this routine
    !     allocates the internal structures. Consequently, this
    !     subroutine should be called immediately, once the top-most
    !     solver structure has been created from parameter file in
    !     order to build the internal structure.
    !
    ! 2.) If the internal structures exist and if, e.g., matrices are
    !     attached, then the temporal vectors, etc. are created
    !     /resised to the correct dimensions. Consequently, the
    !     subroutione should be called whenever the structure of the
    !     problem, e.g., triangulation, discretisation, changes and
    !     requires some modifications of the solver structure.
!</description>

!<inputoutput>
    ! Solver structure
    type(t_solver), intent(inout) :: rsolver
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_solverMultigrid), pointer :: p_solverMultigrid
    integer :: i,ilev,ilbound,iubound,nlmin,nlmax

    ! What kind of solver are we?
    select case(rsolver%csolverType)

    case (SV_COUPLED)
      ! ------------------------------------------------------------------------
      ! Coupled solver - do nothing for the solver itself

      ! Remove update marker from specification bitfields
      rsolver%isolverSpec = iand(rsolver%isolverSpec,&
                                 not(SV_SSPEC_STRUCTURENEEDSUPDATE))

      ! Proceed to solver subnodes
      if (associated(rsolver%p_solverSubnode)) then
        do i = 1, size(rsolver%p_solverSubnode)
          call solver_updateStructure(rsolver%p_solverSubnode(i))
        end do
      end if


    case (SV_FMG)
      ! ------------------------------------------------------------------------
      ! Full multigrid solver - do nothing for the solver itself

      ! Remove update marker from specification bitfields
      rsolver%isolverSpec = iand(rsolver%isolverSpec,&
                                 not(SV_SSPEC_STRUCTURENEEDSUPDATE))

      ! Proceed to coarsegrid solver
      if (associated(rsolver%p_solverMultigrid)) then
        if (associated(rsolver%p_solverMultigrid%p_solverCoarsegrid)) then
          call solver_updateStructure(rsolver%p_solverMultigrid%p_solverCoarsegrid)
        end if
      end if


    case (SV_NONLINEARMG)
      ! ------------------------------------------------------------------------
      ! Nonlinear multigrid solver
      if (associated(rsolver%p_solverMultigrid)) then

        ! Set pointer
        p_solverMultigrid => rsolver%p_solverMultigrid

        ! Check if structure needs update
        if (iand(rsolver%isolverSpec, SV_SSPEC_STRUCTURENEEDSUPDATE) .ne. 0) then

          ! Step 1: Prepare structure
          nlmin = p_solverMultigrid%nlmin
          nlmax = p_solverMultigrid%nlmax

          ! Check if subarray for temporal vectors exists.
          ! Remark: On each level NLMIN+1,...,NLMAX, we need three
          ! temporal vectors for residual on finegrid, residual on
          ! coarsegrid and auxiliary solution vector plus one auxiliary
          ! solution vector on the coarsest level, i.e., NLMIN.
          ! Hence, the array must have dimensions (NLMIN : 3*NLMAX-2*NLMIN)
          if (associated(p_solverMultigrid%RtempVectors)) then

            ! Do we have to reallocate the subarray of matrices?
            ilbound = lbound(p_solverMultigrid%RtempVectors,1)
            iubound = ubound(p_solverMultigrid%RtempVectors,1)

            if ((ilbound > nlmin) .or. (iubound < 3*nlmax-2*nlmin)) then
              do i = ilbound, iubound
                call lsysbl_releaseVector(p_solverMultigrid%RtempVectors(i))
              end do
              deallocate(p_solverMultigrid%RtempVectors)

              ! Allocate subarray(NLMIN:3*NLMAX-2*NLMIN) for temporal vectors
              allocate(p_solverMultigrid%RtempVectors(nlmin:3*nlmax-2*nlmin))
            end if

          else
            ! Allocate subarray(NLMIN:3*NLMAX-2*NLMIN) for temporal vectors
            allocate(p_solverMultigrid%RtempVectors(nlmin:3*nlmax-2*nlmin))
          end if
          ! The initialisation of temporal arrays will be done in the
          ! nonlinear multigrid solver itself since no matrices are
          ! available for nonlinear solvers.

          ! Remove update marker from specification bitfields
          rsolver%isolverSpec = iand(rsolver%isolverSpec,&
                                     not(SV_SSPEC_STRUCTURENEEDSUPDATE))
        end if
      end if


    case (SV_LINEARMG)
      ! ------------------------------------------------------------------------
      ! Linear multigrid solver
      if (associated(rsolver%p_solverMultigrid)) then

        ! Set pointer
        p_solverMultigrid => rsolver%p_solverMultigrid

        ! Check if structure needs update
        if (iand(rsolver%isolverSpec, SV_SSPEC_STRUCTURENEEDSUPDATE) .ne. 0) then

          ! Step 1: prepare structure
          nlmin = p_solverMultigrid%nlmin
          nlmax = p_solverMultigrid%nlmax

          ! Check if subarray for system matrices exist
          if (associated(p_solverMultigrid%rmatrix)) then

            ! Do we have to reallocate the subarray of matrices?
            ilbound = lbound(p_solverMultigrid%rmatrix,1)
            iubound = ubound(p_solverMultigrid%rmatrix,1)

            if ((ilbound > nlmin) .or. (iubound < nlmax)) then
              do i = ilbound, iubound
                call lsysbl_releaseMatrix(p_solverMultigrid%rmatrix(i))
              end do
              deallocate(p_solverMultigrid%rmatrix)

              ! Allocate subarray(NLMIN:NLMAX) for system matrices
              allocate(p_solverMultigrid%rmatrix(nlmin:nlmax))
            end if

          else
            ! Allocate subarray(NLMIN:NLMAX) for system matrices
            allocate(p_solverMultigrid%rmatrix(nlmin:nlmax))
          end if

          ! Check if subarray for temporal vectors exists.
          ! Remark: On each level NLMIN+1,...,NLMAX, we need three
          ! temporal vectors for residual on finegrid, residual on
          ! coarsegrid and auxiliary solution vector plus one auxiliary
          ! solution vector on the coarsest level, i.e., NLMIN.
          ! Hence, the array must have dimensions (NLMIN : 3*NLMAX-2*NLMIN
          if (associated(p_solverMultigrid%RtempVectors)) then

            ! Do we have to reallocate the subarray of matrices?
            ilbound = lbound(p_solverMultigrid%RtempVectors,1)
            iubound = ubound(p_solverMultigrid%RtempVectors,1)

            if ((ilbound > nlmin) .or. (iubound < 3*nlmax-2*nlmin)) then
              do i = ilbound, iubound
                call lsysbl_releaseVector(p_solverMultigrid%RtempVectors(i))
              end do
              deallocate(p_solverMultigrid%RtempVectors)

              ! Allocate subarray(NLMIN:3*NLMAX-2*NLMIN) for temporal vectors
              allocate(p_solverMultigrid%RtempVectors(nlmin:3*nlmax-2*nlmin))
            end if

          else
            ! Allocate subarray(NLMIN:3*NLMAX-2*NLMIN) for temporal vectors
            allocate(p_solverMultigrid%RtempVectors(nlmin:3*nlmax-2*nlmin))
          end if


          ! Ok, now the subarrays definitively exist. Let us try to initialize the
          ! temporal vectors from the attached matrices (if any).
          ! We need the following layout:
          !   raux(NLMIN), raux(NLMIN),   rresc(NLMIN),   rresf(NLMIN+1),
          !                raux(NLMIN+1), rresc(NLMIN+1), rresf(NLMIN+2),
          !                ...            ...             ...
          !                raux(NLMAX-1), rresc(NLMAX-1), rresf(NLMAX)

          ! Loop over levels NLMIN,...,NLMAX-1
          do ilev = nlmin, nlmax-1

            ! Check if we have a usable matrix for current level
            if (p_solverMultigrid%rmatrix(ilev)%NEQ > 0) then

              ! We need three vectors on each level except for the highest level
              do i = 0, 2
                if (p_solverMultigrid%RtempVectors(nlmin+3*(ilev-nlmin)+i)%NEQ > 0) then
                  ! Resize existing vector
                  call lsysbl_resizeVecBlockIndMat(p_solverMultigrid%rmatrix(ilev),&
                      p_solverMultigrid%RtempVectors(nlmin+3*(ilev-nlmin)+i), .false.)
                else
                  ! Create vector from matrix
                  call lsysbl_createVecBlockIndMat(p_solverMultigrid%rmatrix(ilev),&
                      p_solverMultigrid%RtempVectors(nlmin+3*(ilev-nlmin)+i), .false.)
                end if
              end do
            end if
          end do

          ! Finally, create one temporal vectors on the finest grid
          if (p_solverMultigrid%rmatrix(nlmax)%NEQ > 0) then
            if (p_solverMultigrid%RtempVectors(3*nlmax-2*nlmin)%NEQ > 0) then
              ! Resize existing vector
              call lsysbl_resizeVecBlockIndMat(p_solverMultigrid%rmatrix(nlmax),&
                  p_solverMultigrid%RtempVectors(3*nlmax-2*nlmin), .false.)
            else
              ! Create vector from matrix
              call lsysbl_createVecBlockIndMat(p_solverMultigrid%rmatrix(nlmax),&
                  p_solverMultigrid%RtempVectors(3*nlmax-2*nlmin), .false.)
            end if

            ! Remove update marker from specification bitfields
            rsolver%isolverSpec = iand(rsolver%isolverSpec,&
                                       not(SV_SSPEC_STRUCTURENEEDSUPDATE))
          end if
        end if


        ! Update coarsegrid solver
        if (associated(p_solverMultigrid%p_solverCoarsegrid)) then
          call solver_updateStructure(p_solverMultigrid%p_solverCoarsegrid)
        end if


        ! Update smoother subnode
        if (associated(p_solverMultigrid%p_smoother)) then
          call smoother_updateStructure(p_solverMultigrid)
        end if
      end if


    case (SV_NONLINEAR)
      ! ------------------------------------------------------------------------
      ! Nonlinear single-grid solver - do nothing for the solver itself

      ! Remove update marker from specification bitfields
      rsolver%isolverSpec = iand(rsolver%isolverSpec,&
                                 not(SV_SSPEC_STRUCTURENEEDSUPDATE))

      ! Proceed to the generic solver subnode
      if (associated(rsolver%p_solverSubnode)) then
        call solver_updateStructure(rsolver%p_solverSubnode(1))
      end if


    case (SV_LINEAR)
      ! ------------------------------------------------------------------------
      ! Linear single-grid solver
      !
      ! 1.-2. Step: Prepare individual solver subnodes and, eventually,
      !             solver structure for preconditioners

      ! Check if structure needs update
      if (iand(rsolver%isolverSpec, SV_SSPEC_STRUCTURENEEDSUPDATE) .ne. 0) then

        ! UMFPACK solver
        if (associated(rsolver%p_solverUMFPACK)) then

          ! Check if matrix is attached
          if (rsolver%p_solverUMFPACK%rmatrix%NEQ > 0) then

            ! We need one temporal vector
            if (rsolver%p_solverUMFPACK%rtempVector%NEQ > 0) then
              ! Resize existing vector
              call lsysbl_resizeVecBlockIndMat(rsolver%p_solverUMFPACK%rmatrix,&
                                               rsolver%p_solverUMFPACK%rtempVector, .false.)
            else
              ! Create vector from matrix
              call lsysbl_createVecBlockIndMat(rsolver%p_solverUMFPACK%rmatrix,&
                                               rsolver%p_solverUMFPACK%rtempVector, .false.)
            end if

          end if
        end if   ! UMFPACK solver


        ! Jacobi solver
        if (associated(rsolver%p_solverJacobi)) then

          ! Check if matrix is attached
          if (rsolver%p_solverJacobi%rmatrix%NEQ > 0) then

            ! We need one temporal vector
            if (rsolver%p_solverJacobi%rtempVector%NEQ > 0) then
              ! Resize existing vector
              call lsysbl_resizeVecBlockIndMat(rsolver%p_solverJacobi%rmatrix,&
                                               rsolver%p_solverJacobi%rtempVector, .false.)
            else
              ! Create vector from matrix
              call lsysbl_createVecBlockIndMat(rsolver%p_solverJacobi%rmatrix,&
                                               rsolver%p_solverJacobi%rtempVector, .false.)
            end if

          end if
        end if   ! Jacobi solver


        ! (S)SOR solver - nothing to update since no temporal vectors
        !                 are required and no preconditioning is available
        if (associated(rsolver%p_solverSSOR)) then

          ! Check if matrix is attached
          if (rsolver%p_solverSSOR%rmatrix%NEQ > 0) then

            ! We need one temporal vector
            if (rsolver%p_solverSSOR%rtempVector%NEQ > 0) then
              ! Resize existing vector
              call lsysbl_resizeVecBlockIndMat(rsolver%p_solverSSOR%rmatrix,&
                                               rsolver%p_solverSSOR%rtempVector, .false.)
            else
              ! Create vector from matrix
              call lsysbl_createVecBlockIndMat(rsolver%p_solverSSOR%rmatrix,&
                                               rsolver%p_solverSSOR%rtempVector, .false.)
            end if

          end if
        end if   ! SSOR solver


        ! BiCGSTAB solver
        if (associated(rsolver%p_solverBiCGSTAB)) then

          ! Check if matrix is attached
          if (rsolver%p_solverBiCGSTAB%rmatrix%NEQ > 0) then

            ! We need some temporal vectors
            if (rsolver%p_solverBiCGSTAB%RtempVectors(1)%NEQ > 0) then

              ! Resize existing vectors
              do i = lbound(rsolver%p_solverBiCGSTAB%RtempVectors,1),&
                     ubound(rsolver%p_solverBiCGSTAB%RtempVectors,1)
                call lsysbl_resizeVecBlockIndMat(rsolver%p_solverBiCGSTAB%rmatrix,&
                                                 rsolver%p_solverBiCGSTAB%RtempVectors(i), .false.)
              end do

            else

              ! Create vectors from matrix
              do i = lbound(rsolver%p_solverBiCGSTAB%RtempVectors,1),&
                     ubound(rsolver%p_solverBiCGSTAB%RtempVectors,1)
                call lsysbl_createVecBlockIndMat(rsolver%p_solverBiCGSTAB%rmatrix,&
                                                 rsolver%p_solverBiCGSTAB%RtempVectors(i), .false.)
              end do
            end if
          end if

          ! Check if preconditioner is attached
          if (associated(rsolver%p_solverBiCGSTAB%p_precond)) then
            call solver_updateStructure(rsolver%p_solverBiCGSTAB%p_precond)
          end if
        end if   ! BiCGSTAB solver


        ! GMRES solver
        if (associated(rsolver%p_solverGMRES)) then

          ! Check if matrix is attached
          if (rsolver%p_solverGMRES%rmatrix%NEQ > 0) then

            ! We need some temporal vectors
            if (rsolver%p_solverGMRES%rv(1)%NEQ > 0) then
              ! Resize existing vectors rv
              do i = lbound(rsolver%p_solverGMRES%rv,1),&
                     ubound(rsolver%p_solverGMRES%rv,1)
                call lsysbl_resizeVecBlockIndMat(rsolver%p_solverGMRES%rmatrix,&
                                                 rsolver%p_solverGMRES%rv(i), .false.)
              end do
            else
              ! Create vectors rz from matrix
              do i = lbound(rsolver%p_solverGMRES%rv,1),&
                     ubound(rsolver%p_solverGMRES%rv,1)
                call lsysbl_createVecBlockIndMat(rsolver%p_solverGMRES%rmatrix,&
                                                 rsolver%p_solverGMRES%rv(i), .false.)
              end do
            end if

            if (rsolver%p_solverGMRES%rz(1)%NEQ > 0) then
              ! Resize existing vectors rz
              do i = lbound(rsolver%p_solverGMRES%rz,1),&
                     ubound(rsolver%p_solverGMRES%rz,1)
                call lsysbl_resizeVecBlockIndMat(rsolver%p_solverGMRES%rmatrix,&
                                                 rsolver%p_solverGMRES%rz(i), .false.)
              end do
            else
              ! Create vectors rz from matrix
              do i = lbound(rsolver%p_solverGMRES%rz,1),&
                     ubound(rsolver%p_solverGMRES%rz,1)
                call lsysbl_createVecBlockIndMat(rsolver%p_solverGMRES%rmatrix,&
                                                 rsolver%p_solverGMRES%rz(i), .false.)
              end do
            end if
          end if

          ! Check if preconditioner is attached
          if (associated(rsolver%p_solverGMRES%p_precond)) then
            call solver_updateStructure(rsolver%p_solverGMRES%p_precond)
          end if
        end if   ! GMRES solver


        ! ILU solver
        if (associated(rsolver%p_solverILU)) then

          ! Check if matrix is attached
          if (rsolver%p_solverILU%rmatrix%NEQ > 0) then

            ! We need one temporal vector
            if (rsolver%p_solverILU%rtempVector%NEQ > 0) then
              ! Resize existing vector
              call lsysbl_resizeVecBlockIndMat(rsolver%p_solverILU%rmatrix,&
                                               rsolver%p_solverILU%rtempVector, .false.)
            else
              ! Create vector from matrix
              call lsysbl_createVecBlockIndMat(rsolver%p_solverILU%rmatrix,&
                                               rsolver%p_solverILU%rtempVector, .false.)
            end if

          end if
        end if   ! ILU solver

        ! Remove update marker from specification bitfields
        rsolver%isolverSpec = iand(rsolver%isolverSpec,&
                                   not(SV_SSPEC_STRUCTURENEEDSUPDATE))

      end if


    case DEFAULT
      call output_line('Unsupported solver type!',&
          OU_CLASS_ERROR,OU_MODE_STD,'solver_updateStructure')
      call sys_halt()
    end select

  contains

    ! Here, the real working routine follows

    !*************************************************************
    ! Update the structure of the smoother

    subroutine smoother_updateStructure(rsolver)
      type(t_solverMultigrid), intent(inout) :: rsolver

      ! local variables
      type(t_solver) :: rsolverTemplateSmoother
      integer :: ilbound,iubound

      ! What kind of smoother are we?
      select case(rsolver%csmootherType)

      case (SV_NONLINEAR, SV_LINEAR)
        ! ----------------------------------------------------------------------
        ! Nonlinear/Linear smoother subnode

        ! Check if subarray of smoothers exists
        if (associated(rsolver%p_smoother)) then

          ! Do we have to reallocate the subarray of smoothers?
          ilbound = lbound(rsolver%p_smoother,1)
          iubound = ubound(rsolver%p_smoother,1)

          if ((ilbound > rsolver%nlmin+1) .or. (iubound < rsolver%nlmax)) then

            ! Create template solver from smoother
            call solver_createSolver(rsolverTemplateSmoother,&
                                     rsolver%p_smoother(iubound))

            ! Release all existing smoother
            do i = ilbound, iubound
              call solver_releaseSolver(rsolver%p_smoother(i))
            end do
            deallocate(rsolver%p_smoother)

            ! Create new array of smoothers from template solver
            if (rsolver%nlmin+1 .lt. rsolver%nlmax) then
              allocate(rsolver%p_smoother(rsolver%nlmin+1:rsolver%nlmax))
              do i = rsolver%nlmin+1, rsolver%nlmax
                call solver_createSolver(rsolver%p_smoother(i),&
                                         rsolverTemplateSmoother)
              end do
            end if

            ! Release template smoother
            call solver_releaseSolver(rsolverTemplateSmoother)

            ! That is it, the subarray of smoothers exists and has the correct dimensions
          end if

        else
          call output_line('Smoother structure does not exist!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'smoother_updateStructure')
          call sys_halt()
        end if

        ! Ok, now the subarray of smoothers definitively exists and has the correct dimension.
        ! Initialize the smoothers on all levels, whereby each individual smoother is
        ! treated as a standard solver
        do i = lbound(rsolver%p_smoother,1), &
               ubound(rsolver%p_smoother,1)
          call solver_updateStructure(rsolver%p_smoother(i))
        end do


      case DEFAULT
        call output_line('Unsupported solver type!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'smoother_updateStructure')
        call sys_halt()
      end select

    end subroutine smoother_updateStructure

  end subroutine solver_updateStructure

  ! ***************************************************************************

!<subroutine>

  recursive subroutine solver_updateContent(rsolver)

!<description>
    ! This subroutine updates the internal content of the solver.
    ! As a prerequisite the matrices have to be attached to the
    ! corresponding solver nodes.
!</description>

!<inputoutput>
    ! Solver structure
    type(t_solver), intent(inout) :: rsolver
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_solverMultigrid), pointer :: p_solverMultigrid
    integer :: i

    ! What kind of solver are we?
    select case(rsolver%csolverType)
    case (SV_FMG)
      ! ------------------------------------------------------------------------
      ! Full multigrid subnode - do nothing for the solver itself

      ! Remove update marker from specification bitfields
      rsolver%isolverSpec = iand(rsolver%isolverSpec,&
                                 not(SV_SSPEC_CONTENTNEEDSUPDATE))

      ! Proceed to coarsegrid solver
      if (associated(rsolver%p_solverMultigrid)) then
        if (associated(rsolver%p_solverMultigrid%p_solverCoarsegrid)) then
          call solver_updateContent(rsolver%p_solverMultigrid%p_solverCoarsegrid)
        end if
      end if


    case (SV_NONLINEARMG)
      ! ------------------------------------------------------------------------
      ! Nonlinear multigrid solver - do nothing for the solver itself

      ! Remove update marker from specification bitfields
      rsolver%isolverSpec = iand(rsolver%isolverSpec,&
                                 not(SV_SSPEC_CONTENTNEEDSUPDATE))

      if (associated(rsolver%p_solverMultigrid)) then

        ! Set pointer
        p_solverMultigrid => rsolver%p_solverMultigrid

        ! Proceed to coarse grid solver
        if (associated(p_solverMultigrid%p_solverCoarsegrid)) then
          call solver_updateContent(p_solverMultigrid%p_solverCoarsegrid)
        end if

        ! ... and smoother
        if (associated(p_solverMultigrid%p_smoother)) then
          do i = lbound(p_solverMultigrid%p_smoother,1),&
                 ubound(p_solverMultigrid%p_smoother,1)
            call solver_updateContent(p_solverMultigrid%p_smoother(i))
          end do
        end if
      end if


    case (SV_NONLINEAR)
      ! ------------------------------------------------------------------------
      ! Nonlinear single-grid solver - do nothing for the solver itself

      ! Remove update marker from specification bitfields
      rsolver%isolverSpec = iand(rsolver%isolverSpec,&
                                 not(SV_SSPEC_CONTENTNEEDSUPDATE))

      ! Porceed to the generic solver subnode
      if (associated(rsolver%p_solverSubnode)) then
        call solver_updateContent(rsolver%p_solverSubnode(1))
      end if


    case (SV_LINEARMG)
      ! ------------------------------------------------------------------------
      ! Linear multigrid solver - do nothing for the solver itself

      ! Remove update marker from specification bitfields
      rsolver%isolverSpec = iand(rsolver%isolverSpec,&
                                 not(SV_SSPEC_CONTENTNEEDSUPDATE))

      if (associated(rsolver%p_solverMultigrid)) then

        ! Set pointer
        p_solverMultigrid => rsolver%p_solverMultigrid

        ! Proceed to coarse grid solver
        if (associated(p_solverMultigrid%p_solverCoarsegrid)) then
          call solver_updateContent(p_solverMultigrid%p_solverCoarsegrid)
        end if

        ! ... and smoother
        if (associated(p_solverMultigrid%p_smoother)) then
          do i = lbound(p_solverMultigrid%p_smoother,1),&
                 ubound(p_solverMultigrid%p_smoother,1)
            call solver_updateContent(p_solverMultigrid%p_smoother(i))
          end do
        end if
      end if


    case (SV_LINEAR)
      ! ------------------------------------------------------------------------
      ! Linear single-grid solver

      ! Check if content needs update
      if (iand(rsolver%isolverSpec, SV_SSPEC_CONTENTNEEDSUPDATE) .ne. 0) then

        ! Initialize UMFPACK solver
        if (associated(rsolver%p_solverUMFPACK)) then
          call solver_initUMFPACK(rsolver%p_solverUMFPACK)
        end if

        ! Initialize ILU factorisation
        if (associated(rsolver%p_solverILU)) then
          call solver_initILU(rsolver%p_solverILU)
        end if

        ! Initialize preconditioner of BiCGSTAB solver
        if (associated(rsolver%p_solverBiCGSTAB)) then
          if (associated(rsolver%p_solverBiCGSTAB%p_precond)) then
            call solver_updateContent(rsolver%p_solverBiCGSTAB%p_precond)
          end if
        end if

        ! Initialize preconditioner of GMRES solver
        if (associated(rsolver%p_solverGMRES)) then
          if (associated(rsolver%p_solverGMRES%p_precond)) then
            call solver_updateContent(rsolver%p_solverGMRES%p_precond)
          end if
        end if

        ! Remove update marker from specification bitfields
        rsolver%isolverSpec = iand(rsolver%isolverSpec,&
                                   not(SV_SSPEC_CONTENTNEEDSUPDATE))
      end if


    case DEFAULT
      call output_line('Unsupported solver type!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'solver_updateContent')
      call sys_halt()
    end select

  end subroutine solver_updateContent

  ! ***************************************************************************

!<subroutine>

  subroutine solver_prolongationScalar(rtriangulationCoarse, rtriangulationFine,&
                                       rxCoarse, rxFine)

!<description>
    ! This subroutine performs prolongation from coarse grid to fine grid, i.e.
    ! rx(Fine) := prolongation(rx(Coarse)).
    !
    ! Note that prolongation can be performed in situ if both vectors are the same.
!</description>

!<input>
    ! Coarse triangulation structure
    type(t_triangulation), intent(in) :: rtriangulationCoarse

    ! Fine triangulation structure
    type(t_triangulation), intent(in) :: rtriangulationFine

    ! Coarse mesh vector
    type(t_vectorScalar), intent(in) :: rxCoarse
!</input>

!<inputoutput>
    ! Fine mesh vector
    type(t_vectorScalar), intent(inout) :: rxFine
!</inputoutput>
!</subroutine>

    ! Local variables
    integer, dimension(:,:), pointer :: p_IneighboursAtElementFine
    integer, dimension(:,:), pointer :: p_IneighboursAtElementCoarse
    integer, dimension(:,:), pointer :: p_IverticesAtElementFine
    integer, dimension(:,:), pointer :: p_IverticesAtElementCoarse
    real(DP), dimension(:), pointer :: p_DxFine,p_DxCoarse
    integer :: nvar

    ! Set global pointers
    call storage_getbase_int2d (rtriangulationFine%h_IneighboursAtElement,&
                                p_IneighboursAtElementFine)
    call storage_getbase_int2d (rtriangulationCoarse%h_IneighboursAtElement,&
                                p_IneighboursAtElementCoarse)
    call storage_getbase_int2d (rtriangulationFine%h_IverticesAtElement,&
                                p_IverticesAtElementFine)
    call storage_getbase_int2d (rtriangulationCoarse%h_IverticesAtElement,&
                                p_IverticesAtElementCoarse)

    nvar = rxCoarse%NVAR

    ! Check if scalar vectors are compatible
    if (rxFine%NVAR .ne. rxCoarse%NVAR) then
      call output_line('Scalar vectors must be compatible!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'solver_prolongationScalar')
      call sys_halt()
    end if

    ! Set pointers
    call lsyssc_getbase_double(rxFine, p_DxFine)
    call lsyssc_getbase_double(rxCoarse, p_DxCoarse)

    ! Explicitely adopt values from coarser level
    call lalg_copyVectorDble(p_DxCoarse, p_DxFine)

    ! Do the prolongation
    call do_prolongationP1Q1(nvar, rtriangulationCoarse%NEL,&
                             p_IneighboursAtElementFine, p_IneighboursAtElementCoarse,&
                             p_IverticesAtElementFine, p_IverticesAtElementCoarse, &
                             p_DxFine, p_DxCoarse)

  contains

    ! Here, the real working routine follows

    !*************************************************************
    ! Perform prolongation from coarse to fine grid
    ! for P1 and Q1 finite elements

    subroutine do_prolongationP1Q1(nvar, nel, IneighboursAtElementFine, &
                                   IneighboursAtElementCoarse, IverticesAtElementFine, &
                                   IverticesAtElementCoarse, DxFine, DxCoarse)
      integer, intent(in) :: nvar
      integer, intent(in) :: nel
      integer, dimension(:,:), intent(in) :: IneighboursAtElementFine
      integer, dimension(:,:), intent(in) :: IneighboursAtElementCoarse
      integer, dimension(:,:), intent(in) :: IverticesAtElementFine
      integer, dimension(:,:), intent(in) :: IverticesAtElementCoarse
      real(DP), dimension(nvar,*), intent(in) :: DxCoarse
      real(DP), dimension(nvar,*), intent(inout) :: DxFine

      real(DP), dimension(nvar,TRIA_MAXNVE2D) :: Dxloc
      integer :: iel,iel1,iel2,iel3,iel4,nve

      ! Loop over all elements of the coarse grid
      do iel = 1, nel

        ! Determine number and values of local vertices
        nve = tria_getNVE(IverticesAtElementCoarse, iel)
        Dxloc(:,1:nve) = DxCoarse(:,IverticesAtElementCoarse(1:nve,iel))

        ! What kind of element are we?
        select case(nve)
        case(TRIA_NVETRI2D)
          ! linear triangle
          if (IneighboursAtElementCoarse(1,iel) < iel)&
              DxFine(:,IverticesAtElementFine(1,iel)) = 0.5_DP*(Dxloc(:,1)+Dxloc(:,2))

          if (IneighboursAtElementCoarse(2,iel) < iel)&
              DxFine(:,IverticesAtElementFine(2,iel)) = 0.5_DP*(Dxloc(:,2)+Dxloc(:,3))

          if (IneighboursAtElementCoarse(3,iel) < iel)&
              DxFine(:,IverticesAtElementFine(3,iel)) = 0.5_DP*(Dxloc(:,3)+Dxloc(:,1))


        case(TRIA_NVEQUAD2D)
          ! bilinear quadrilateral
          iel1 = iel
          iel2 = IneighboursAtElementFine(2,iel1)
          iel3 = IneighboursAtElementFine(2,iel2)
          iel4 = IneighboursAtElementFine(2,iel3)

          if (IneighboursAtElementCoarse(1,iel ) < iel)&
              DxFine(:,IverticesAtElementFine(2,iel1)) = 0.5_DP*(Dxloc(:,1)+Dxloc(:,2))

          if (IneighboursAtElementCoarse(2,iel ) < iel)&
              DxFine(:,IverticesAtElementFine(2,iel2)) = 0.5_DP*(Dxloc(:,2)+Dxloc(:,3))

          if (IneighboursAtElementCoarse(3,iel ) < iel)&
              DxFine(:,IverticesAtElementFine(2,iel3)) = 0.5_DP*(Dxloc(:,3)+Dxloc(:,4))

          if (IneighboursAtElementCoarse(4,iel ) < iel)&
              DxFine(:,IverticesAtElementFine(2,iel4)) = 0.5_DP*(Dxloc(:,4)+Dxloc(:,1))

          DxFine(:,IverticesAtElementFine(3,iel)) = 0.25_DP*sum(Dxloc,2)

        case DEFAULT
          call output_line('Invalid number of vertices per element!',&
              OU_CLASS_ERROR,OU_MODE_STD,'do_prolongationP1Q1')
          call sys_halt()
        end select
      end do
    end subroutine do_prolongationP1Q1

  end subroutine solver_prolongationScalar

  ! ***************************************************************************

!<subroutine>

  subroutine solver_prolongationBlock(rtriangulationCoarse, rtriangulationFine,&
                                      rxCoarse, rxFine)

!<description>
    ! This subroutine performs prolongation from coarse grid to fine grid, i.e.
    ! rx(Fine) := prolongation(rx(Coarse)).
    !
    ! Note that prolongation can be performed in situ if both vectors are the same.
!</description>

!<input>
    ! Coarse triangulation structure
    type(t_triangulation), intent(in) :: rtriangulationCoarse

    ! Fine triangulation structure
    type(t_triangulation), intent(in) :: rtriangulationFine

    ! Coarse mesh vector
    type(t_vectorBlock), intent(in) :: rxCoarse
!</input>

!<inputoutput>
    ! Fine mesh vector
    type(t_vectorBlock), intent(inout) :: rxFine
!</inputoutput>
!</subroutine>

    ! Local variables
    integer :: iblock

    if (rxFine%nblocks .ne. rxCoarse%nblocks) then
      call output_line('Block vectors must be compatible!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'solver_prolongationBlock')
      call sys_halt()
    end if

    ! Treat each block separately
    do iblock = 1, rxCoarse%nblocks
      call solver_prolongationScalar(rtriangulationCoarse, rtriangulationFine,&
                                     rxCoarse%RvectorBlock(iblock),&
                                     rxFine%RvectorBlock(iblock))
    end do

  end subroutine solver_prolongationBlock

  ! ***************************************************************************

!<subroutine>

  subroutine solver_restrictionScalar(rtriangulationFine, rtriangulationCoarse,&
                                      rxFine, rxCoarse)

!<description>
    ! This subroutine performs restriction from fine grid to coarse grid, i.e.
    ! rx(Coarse) := restriction(rx(Fine)).
    !
    ! Note that restriction can be performed in situ if both vectors are the same.
!</description>

!<input>
    ! Coarse triangulation structure
    type(t_triangulation), intent(in) :: rtriangulationFine

    ! Fine triangulation structure
    type(t_triangulation), intent(in) :: rtriangulationCoarse

    ! Coarse mesh vector
    type(t_vectorScalar), intent(in) :: rxFine
!</input>

!<inputoutput>
    ! Fine mesh vector
    type(t_vectorScalar), intent(inout) :: rxCoarse
!</inputoutput>
!</subroutine>

    ! local variables
    integer, dimension(:,:), pointer :: p_IneighboursAtElementFine
    integer, dimension(:,:), pointer :: p_IverticesAtElementFine
    real(DP), dimension(:), pointer :: p_DxFine,p_DxCoarse
    integer :: nvar

    ! Set global pointers
    call storage_getbase_int2D (rtriangulationFine%h_IneighboursAtElement,&
                                p_IneighboursAtElementFine)
    call storage_getbase_int2D (rtriangulationFine%h_IverticesAtElement,&
                                p_IverticesAtElementFine)

    nvar = rxFine%NVAR

    ! Check if scalar vectors are compatible
    if (rxFine%NVAR .ne. rxCoarse%NVAR) then
      call output_line('Scalar vectors must be compatible!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'solver_restrictionScalar')
      call sys_halt()
    end if

    ! Set pointers
    call lsyssc_getbase_double(rxFine, p_DxFine)
    call lsyssc_getbase_double(rxCoarse, p_DxCoarse)

    ! Explicitely adopt values 1:NVT_Coarse from coarser level.
    call lalg_copyVectorDble(p_DxFine(1:nvar*rtriangulationCoarse%NVT), p_DxCoarse)

    call do_restrictionP1Q1(nvar, rtriangulationFine%NEL, rtriangulationCoarse%NEL, &
                            p_IneighboursAtElementFine, p_IverticesAtElementFine,&
                            p_DxFine, p_DxCoarse)

    ! If restriction is performed in situ, then you must nor forget so
    ! adjust the dimension of the coarse grid vector outside of this routine

  contains

    ! Here, the real working routine follows

    !*************************************************************
    ! Perform restriction from fine to coarse grid
    ! for P1 and Q1 finite elements

    subroutine do_restrictionP1Q1(nvar, nelFine, nelCoarse, &
                                  IneighboursAtElementFine, IverticesAtElementFine,&
                                  DxFine, DxCoarse)
      integer, intent(in) :: nvar
      integer, intent(in) :: nelFine,nelCoarse
      integer, dimension(:,:), intent(in) :: IneighboursAtElementFine
      integer, dimension(:,:), intent(in) :: IverticesAtElementFine
      real(DP), dimension(nvar,*), intent(in) :: DxFine
      real(DP), dimension(nvar,*), intent(inout) :: DxCoarse

      integer, dimension(TRIA_MAXNVE2D) :: iloc
      integer :: iel,nve

      ! Loop over all elements of the fine grid
      do iel = 1, nelFine

        ! Determine number local vertices
        nve = tria_getNVE(IverticesAtElementFine, iel)

        ! Check if current element is a triangle
        ! Then, ignore element if IEL < NEL of the coarse grid
        if (nve .eq. TRIA_NVETRI2D .and. iel < nelCoarse+1) cycle

        ! Determine number of local vertices
        iloc(1:nve) = IverticesAtElementFine(1:nve,iel)

        ! What kind of element are we?
        select case(nve)
        case(TRIA_NVETRI2D)
          ! linear triangle
          if (IneighboursAtElementFine(1,iel) < iel)&
              DxCoarse(:,iloc(1)) = DxCoarse(:,iloc(1)) + 0.5_DP*DxFine(:,iloc(2))

          if (IneighboursAtElementFine(3,iel) < iel)&
              DxCoarse(:,iloc(1)) = DxCoarse(:,iloc(1)) + 0.5_DP*DxFine(:,iloc(3))


        case(TRIA_NVEQUAD2D)
          ! bilinear quadrilateral
          DxCoarse(:,iloc(1)) = DxCoarse(:,iloc(1)) + 0.25_DP*sum(DxFine(:,iloc(2:4)),2)

          if (IneighboursAtElementFine(1,iel) .eq. 0)&
              DxCoarse(:,iloc(1)) = DxCoarse(:,iloc(1))+0.25_DP*DxFine(:,iloc(2))

          if (IneighboursAtElementFine(4,iel) .eq. 0)&
              DxCoarse(:,iloc(1)) = DxCoarse(:,iloc(1)) + 0.25_DP*DxFine(:,iloc(4))


        case DEFAULT
          call output_line('Invalid number of vertices per element!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'do_restrictionP1Q1')
          call sys_halt()
        end select
      end do
    end subroutine do_restrictionP1Q1

  end subroutine solver_restrictionScalar

  ! ***************************************************************************

!<subroutine>

  subroutine solver_restrictionBlock(rtriangulationFine, rtriangulationCoarse,&
                                     rxFine, rxCoarse)

!<description>
    ! This subroutine performs restriction from fine grid to coarse grid, i.e.
    ! rx(Coarse) := restriction(rx(Fine)).
    !
    ! Note that restriction can be performed in situ if both vectors are the same.
!</description>

!<input>
    ! Coarse triangulation structure
    type(t_triangulation), intent(in) :: rtriangulationFine

    ! Fine triangulation structure
    type(t_triangulation), intent(in) :: rtriangulationCoarse

    ! Coarse mesh vector
    type(t_vectorBlock), intent(in) :: rxFine
!</input>

!<inputoutput>
    ! Fine mesh vector
    type(t_vectorBlock), intent(inout) :: rxCoarse
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: iblock

    if (rxFine%nblocks .ne. rxCoarse%nblocks) then
      call output_line('Block vectors must be compatible!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'solver_restrictionBlock')
      call sys_halt()
    end if

    ! Treat each block separately
    do iblock = 1, rxFine%nblocks

      call solver_restrictionScalar(rtriangulationFine, rtriangulationCoarse,&
                                    rxFine%RvectorBlock(iblock),&
                                    rxCoarse%RvectorBlock(iblock))
    end do

  end subroutine solver_restrictionBlock

  ! ***************************************************************************

!<subroutine>

  subroutine solver_copySolver(rsolverSource, rsolverDest,&
                               bcopyInput, bcopyOutput)

!<description>
    ! This subroutine copies some data from rsolverSource to rsolverDest.
    !
    ! We follow the following convention:
    !  Input   = all input parameters of the solver, e.g. tolerances
    !  Output  = all output parameters of the solver, e.g. residuals
    !
    ! Note that this routine does not transfer the structure of the
    ! source solver to the destination solver, that is, the type of solver,
    ! preconditioner, etc. will not be copied. If this should be done
    ! the solver_createSolverIndirect routine should be used instead.
!</description>

!<input>
    ! The source solver structure
    type(t_solver), intent(in) :: rsolverSource

    ! Copy flag for input parameters
    logical, intent(in) :: bcopyInput

    ! Copy flag for output parameters
    logical, intent(in) :: bcopyOutput
!</input>

!<inputoutput>
    ! The destination solver structure
    type(t_solver), intent(inout) :: rsolverDest
!</inputoutput>

!</subroutine>

    ! Copy input parameters?
    if (bcopyInput) then
      rsolverDest%iresNorm           = rsolverSource%iresNorm
      rsolverDest%istoppingCriterion = rsolverSource%istoppingCriterion
      rsolverDest%nminIterations     = rsolverSource%nminIterations
      rsolverDest%nmaxIterations     = rsolverSource%nmaxIterations
      rsolverDest%domega             = rsolverSource%domega
      rsolverDest%depsRel            = rsolverSource%depsRel
      rsolverDest%depsAbs            = rsolverSource%depsAbs
      rsolverDest%depsStag           = rsolverSource%depsStag
      rsolverDest%ddivRel            = rsolverSource%ddivRel
      rsolverDest%ddivAbs            = rsolverSource%ddivAbs
      rsolverDest%drhsZero           = rsolverSource%drhsZero
      rsolverDest%ddefZero           = rsolverSource%ddefZero
    end if

    ! Copy output parameters?
    if (bcopyOutput) then
      rsolverDest%istatus            = rsolverSource%istatus
      rsolverDest%iiterations        = rsolverSource%iiterations
      rsolverDest%dinitialRHS        = rsolverSource%dinitialRHS
      rsolverDest%dinitialDefect     = rsolverSource%dinitialDefect
      rsolverDest%dfinalDefect       = rsolverSource%dfinalDefect
      rsolverDest%dconvergenceRate   = rsolverSource%dconvergenceRate
    end if

  end subroutine solver_copySolver

  ! ***************************************************************************

!<function>

  function solver_testConvergence(rsolver) result(bconverged)

!<description>
    ! This function check all convergence criteria for the solver
!</description>

!<input>
    ! solver structure
    type(t_solver), intent(in) :: rsolver
!</input>

!<result>
    ! Boolean value.
    ! =TRUE if the convergence criterion is reached;
    ! =FALSE otherwise.
    logical :: bconverged
!</result>
!</function>

    select case (rsolver%istoppingCriterion)

    case (SV_STOP_ANY)
      ! Iteration stops if either the absolute or the relative
      ! criterium holds.
      bconverged = .false.

      ! Absolute convergence criterion? Check the norm directly.
      if (rsolver%depsAbs .ne. 0.0_DP) then
        if (.not. (rsolver%dfinalDefect .gt. rsolver%depsAbs)) then
          bconverged = .true.
          return
        end if
      end if

      ! Relative convergence criterion? Multiply with initial residuum
      ! and check the norm.
      if (rsolver%depsRel .ne. 0.0_DP) then
        if (.not. (rsolver%dfinalDefect .gt. rsolver%depsRel *&
                   rsolver%dinitialDefect)) then
          bconverged = .true.
          return
        end if
      end if

    case (SV_STOP_ALL)
      ! Iteration stops if both the absolute and the relative criterium holds.
      bconverged = .true.

      ! Absolute convergence criterion? Check the norm directly.
      if (rsolver%depsAbs .ne. 0.0_DP) then
        if (rsolver%dfinalDefect .gt. rsolver%depsAbs) then
          bconverged = .false.
          return
        end if
      end if

      ! Relative convergence criterion? Multiply with initial residuum
      ! and check the norm.
      if (rsolver%depsRel .ne. 0.0_DP) then
        if (rsolver%dfinalDefect .gt. rsolver%depsRel *&
            rsolver%dinitialDefect) then
          bconverged = .false.
          return
        end if
      end if

    case default
      call output_line('Invalid stopping criterion',&
          OU_CLASS_ERROR,OU_MODE_STD,'solver_testConvergence')
      call sys_halt()
    end select

  end function solver_testConvergence

  ! ***************************************************************************

!<function>

  function solver_testDivergence(rsolver) result(bdiverged)

!<description>
    ! This function check all divergence criteria for the solver
!</description>

!<input>
    ! solver structure
    type(t_solver), intent(in) :: rsolver
!</input>

!<result>
    ! Boolean value.
    ! =TRUE if the divergence criterion holds;
    ! =FALSE otherwise.
    logical :: bdiverged
!</result>
!</function>

    bdiverged = .false.

    ! Absolute divergence criterion? Check the norm directly.
    if (rsolver%ddivAbs .ne. SYS_INFINITY) then

      ! use NOT here - gives a better handling of special cases like NaN!
      if ( .not. (rsolver%dfinalDefect .le. rsolver%ddivAbs)) then
        bdiverged = .true.
        return
      end if

    end if

    ! Relative divergence criterion? Multiply with initial residuum
    ! and check the norm.
    if (rsolver%depsRel .ne. SYS_INFINITY) then

      if ( .not. (rsolver%dfinalDefect .le. rsolver%dinitialDefect*rsolver%ddivRel) ) then
        bdiverged = .true.
        return
      end if

    end if

  end function solver_testDivergence

  ! ***************************************************************************

!<function>

  function solver_testStagnation(rsolver, doldDefect) result(bstagnated)

!<description>
    ! This function check all stagnation criteria for the solver
!</description>

!<input>
    ! solver structure
    type(t_solver), intent(in) :: rsolver

    ! norm of previous defect
    real(DP), intent(in) :: doldDefect
!</input>

!<result>
    ! Boolean value.
    ! =TRUE if the stagnation criterion holds;
    ! =FALSE otherwise.
    logical :: bstagnated
!</result>
!</function>

    bstagnated = .false.

    ! Stagnation criterion?
    if (rsolver%depsStag .ne. 0.0_DP) then
      if (rsolver%dfinalDefect .gt. rsolver%depsStag*doldDefect) then
        bstagnated = .true.
        return
      end if
    end if

  end function solver_testStagnation

  ! ***************************************************************************
  ! Auxiliary routines
  ! ***************************************************************************

!<subroutine>

  recursive subroutine solver_initILU(rsolver)

!<description>
    ! This subroutine initializes the ILU solver.
    ! The subroutine requires the matrix to be set before.
    ! Otherwise, it will terminate with an error.
!</description>

!<inputoutput>
    ! ILU-solver structure
    type(t_solverILU), intent(inout) :: rsolver
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_matrixScalar), pointer :: p_rmatrix
    integer :: i,isize,nblock


    ! Do we have a block-matrix?
    if ((rsolver%rmatrix%nblocksPerCol .ne. 1) .or.&
        (rsolver%rmatrix%nblocksPerRow .ne. 1)) then

      nblock = min(rsolver%rmatrix%nblocksPerCol,&
                   rsolver%rmatrix%nblocksPerRow)

      ! Check if block structure exists
      if (associated(rsolver%p_rsolverBlockILU)) then

        ! Check if block structure has correct dimensions
        isize = size(rsolver%p_rsolverBlockILU)
        if (nblock .ne. isize) then
          do i = 1, isize
            call release_solverILU(rsolver%p_rsolverBlockILU(i))
          end do
          deallocate(rsolver%p_rsolverBlockILU)
          allocate(rsolver%p_rsolverBlockILU(nblock))
        end if

      else

        allocate(rsolver%p_rsolverBlockILU(nblock))

      end if

      ! Loop over all diagonal blocks
      do i = 1, nblock

        ! Fill solver
        rsolver%p_rsolverBlockILU(i)%ifill   = rsolver%ifill
        rsolver%p_rsolverBlockILU(i)%domega  = rsolver%domega
        rsolver%p_rsolverBlockILU(i)%depsILU = rsolver%depsILU

        ! Attach submatrix from diagonal block
        call lsysbl_createMatFromScalar(rsolver%rmatrix%RmatrixBlock(i,i),&
                                        rsolver%p_rsolverBlockILU(i)%rmatrix)

        ! Initialize ILU solver for submatrix
        call solver_initILU(rsolver%p_rsolverBlockILU(i))
      end do

    else

      ! Set pointer to scalar matrix
      p_rmatrix => rsolver%rmatrix%RmatrixBlock(1,1)
      call do_scalarILU(rsolver, p_rmatrix)

    end if

  contains

    ! Here, the working routine follows

    !**************************************************************
    ! Release ILU solver
    recursive subroutine release_solverILU(rsolver)

      type(t_solverILU), intent(inout) :: rsolver

      integer :: i

      ! Release structure for block ILU solver
      if (associated(rsolver%p_rsolverBlockILU)) then
        do i = lbound(rsolver%p_rsolverBlockILU,1),&
               ubound(rsolver%p_rsolverBlockILU,1)
          call release_solverILU(rsolver%p_rsolverBlockILU(i))
        end do
        deallocate(rsolver%p_rsolverBlockILU)
      end if

      ! Release handles
      if (rsolver%h_Idata .ne. ST_NOHANDLE) call storage_free(rsolver%h_Idata)
      if (rsolver%h_Ddata .ne. ST_NOHANDLE) call storage_free(rsolver%h_Ddata)

      ! Release matrix
      call lsysbl_releaseMatrix(rsolver%rmatrix)

      ! Release temporal vector
      call lsysbl_releaseVector(rsolver%rtempVector)
    end subroutine release_solverILU


    !**************************************************************
    ! Initialize (M)ILU(s) solver for scalar matrices.
    ! This subroutine initializes the (shifted, modified) ILU
    ! decomposition or the ILU(s) decomposition with fill-in.
    subroutine do_scalarILU(rsolver, rmatrix)

      type(t_solverILU), intent(inout) :: rsolver
      type(t_matrixScalar), intent(INout) :: rmatrix

      real(DP), dimension(:), pointer :: p_DA,p_Ddata
      integer, dimension(:), pointer :: p_Kld
      integer, dimension(:), pointer :: p_Kcol
      integer, dimension(:), pointer :: p_Kdiagonal
      integer :: isize


      ! Set pointer
      call lsyssc_getbase_double(rmatrix, p_DA)

      ! What kind of ILU algorithm should be used?
      if (rsolver%ifill .le. 0) then

        ! What kind of matrix are we?
        select case(rmatrix%cmatrixFormat)

        case (LSYSSC_MATRIX7, LSYSSC_MATRIX9)

          ! (Re-)allocate required memory for matrix data
          if (rsolver%h_Ddata .eq. ST_NOHANDLE) then
            call storage_new('do_scalarILU', 'h_Ddata', rmatrix%NA,&
                             ST_DOUBLE, rsolver%h_Ddata, ST_NEWBLOCK_NOINIT)
          else
            call storage_getsize(rsolver%h_Ddata, isize)
            if (isize < rmatrix%NA) then
              call storage_realloc('do_scalarILU', rmatrix%NA,&
                                   rsolver%h_Ddata, ST_NEWBLOCK_NOINIT, .false.)
            end if
          end if


          ! Copy matrix data
          call storage_getbase_double(rsolver%h_Ddata, p_Ddata, rmatrix%NA)
          call lalg_copyVectorDble(p_DA, p_Ddata)

          ! Set pointers
          call lsyssc_getbase_Kld (rmatrix, p_Kld)
          call lsyssc_getbase_Kcol(rmatrix, p_Kcol)


          if (rmatrix%cmatrixFormat .eq. LSYSSC_MATRIX7) then

            ! Perform (M)ILU(0) factorisation for matrix format 7
            call do_mat7MILU0(p_Ddata, p_Kcol, p_Kld,&
                              rmatrix%NEQ, rmatrix%NVAR, abs(rsolver%ifill),&
                              abs(rsolver%domega), rsolver%depsILU)

          else

            ! Perform (M)ILU(0) factorisation for matrix format 9
            call lsyssc_getbase_Kdiagonal(rmatrix, p_Kdiagonal)

            call do_mat9MILU0(p_Ddata, p_Kcol, p_Kld, p_Kdiagonal,&
                              rmatrix%NEQ, rmatrix%NVAR, abs(rsolver%ifill),&
                              abs(rsolver%domega), rsolver%depsILU)

          end if


        case (LSYSSC_MATRIX7INTL, LSYSSC_MATRIX9INTL)

          ! What kind of interleave format are we?
          select case(rmatrix%cinterleavematrixFormat)

          case (LSYSSC_MATRIXD)

            ! (Re-)allocate required memory for matrix data
            if (rsolver%h_Ddata .eq. ST_NOHANDLE) then
              call storage_new('do_scalarILU', 'h_Ddata', rmatrix%NA*rmatrix%NVAR,&
                               ST_DOUBLE, rsolver%h_Ddata, ST_NEWBLOCK_NOINIT)
            else
              call storage_getsize(rsolver%h_Ddata, isize)
              if (isize < rmatrix%NA*rmatrix%NVAR) then
                call storage_realloc('do_scalarILU', rmatrix%NA*rmatrix%NVAR,&
                                     rsolver%h_Ddata, ST_NEWBLOCK_NOINIT, .false.)
              end if
            end if


            ! Copy matrix data
            call storage_getbase_double(rsolver%h_Ddata, p_Ddata,&
                                        rmatrix%NA*rmatrix%NVAR)
            call lalg_copyVectorDble(p_DA, p_Ddata)

            ! Set pointers
            call lsyssc_getbase_Kld (rmatrix, p_Kld)
            call lsyssc_getbase_Kcol(rmatrix, p_Kcol)

            if (rmatrix%cmatrixFormat .eq. LSYSSC_MATRIX7INTL) then

              ! Perform (M)ILU(0) factorisation for matrix format 7
              call do_mat7MILU0(p_Ddata, p_Kcol, p_Kld,&
                                rmatrix%NEQ, rmatrix%NVAR, abs(rsolver%ifill),&
                                abs(rsolver%domega), rsolver%depsILU)

            else

              ! Perform (M)ILU(0) factorisation for matrix format 9
              call lsyssc_getbase_Kdiagonal(rmatrix, p_Kdiagonal)
              call do_mat9MILU0(p_Ddata, p_Kcol, p_Kld, p_Kdiagonal,&
                                rmatrix%NEQ, rmatrix%NVAR, abs(rsolver%ifill),&
                                abs(rsolver%domega), rsolver%depsILU)

            end if


          case (LSYSSC_MATRIX1)

            ! (Re-)allocate required memory for matrix data
            if (rsolver%h_Ddata .eq. ST_NOHANDLE) then
              call storage_new('do_scalarILU', 'h_Ddata',&
                               rmatrix%NEQ*rmatrix%NVAR*rmatrix%NVAR,&
                               ST_DOUBLE, rsolver%h_Ddata, ST_NEWBLOCK_NOINIT)
            else
              call storage_getsize(rsolver%h_Ddata, isize)
              if (isize < rmatrix%NEQ*rmatrix%NVAR*rmatrix%NVAR) then
                call storage_realloc('do_scalarILU',&
                                     rmatrix%NEQ*rmatrix%NVAR*rmatrix%NVAR,&
                                     rsolver%h_Ddata, ST_NEWBLOCK_NOINIT, .false.)
              end if
            end if

            ! Copy matrix data (explicitly)
            call storage_getbase_double(rsolver%h_Ddata, p_Ddata,&
                                        rmatrix%NEQ*rmatrix%NVAR*rmatrix%NVAR)

            if (rmatrix%cmatrixFormat .eq. LSYSSC_MATRIX7INTL) then

              ! Perform BILU(0) factorisation for matrix format 7
              call lsyssc_getbase_Kld (rmatrix, p_Kld)
              call do_mat79Intl1BILU0(p_DA, p_Ddata, p_Kld,&
                                      rmatrix%NEQ, rmatrix%NVAR)

            else

              ! Perform BILU(0) factorisation for matrix format 9
              call lsyssc_getbase_Kdiagonal(rmatrix, p_Kdiagonal)
              call do_mat79Intl1BILU0(p_DA, p_Ddata, p_Kdiagonal,&
                                      rmatrix%NEQ, rmatrix%NVAR)

            end if

          case DEFAULT
            call output_line('Unsupported interleave matrix format!',&
                             OU_CLASS_ERROR,OU_MODE_STD,'do_scalarILU')
            call sys_halt()
          end select

        case DEFAULT
          call output_line('Unsupported matrix format!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'do_scalarILU')
          call sys_halt()
        end select

      else

        ! What kind of matrix are we?
        select case(rmatrix%cmatrixFormat)

        case (LSYSSC_MATRIX7,LSYSSC_MATRIX9)

          ! Set pointers
          call lsyssc_getbase_Kld(rmatrix, p_Kld)
          call lsyssc_getbase_Kcol(rmatrix, p_Kcol)

          ! Perform (M)ILU(s), s>0 factorisation for matrix format 7/9
          call do_mat79MILUs(p_DA, p_Kcol, p_Kld,&
                             rmatrix%NEQ, rmatrix%NA, rsolver%ifill,&
                             abs(rsolver%domega), rsolver%lu, rsolver%jlu,&
                             rsolver%ilup, rsolver%h_Idata)

        case DEFAULT
          call output_line('Unsupported matrix format!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'do_scalarILU')
          call sys_halt()
        end select

      end if
    end subroutine do_scalarILU


    !**************************************************************
    ! Calculate (M)ILU(s) decomposition of a matrix in format CSR7/CSR9

    subroutine do_mat79MILUs(Da, Kcol, Kld, neq, na, ifill,&
                             domega, lu, jlu, ilup, h_Idata)
      real(DP), dimension(:), intent(inout) :: Da
      integer, dimension(:), intent(in) :: Kcol
      integer, dimension(:), intent(in) :: Kld
      integer, intent(in) :: neq
      integer, intent(in) :: na
      integer, intent(in) :: ifill
      real(DP), intent(in) :: domega
      integer(I32), intent(out) :: lu,jlu,ilup
      integer, intent(inout) :: h_Idata

      ! Declare our ILUS-routine from SPLIB as interface to be sure,
      ! parameter interfaces are checked by compiler
      interface
      subroutine ilus(n,a,colind,rwptr,&
                      s,relax,&
                      lu,jlu,ilup,&
                      iwork,maxstr,&
                      ierr,mneed)
        use fsystem

        ! Integer precision for ILU solver
        integer, parameter :: LINSOL_PRECOND_ILUINT = I32

        ! Double precision for ILU solver
        integer, parameter :: LINSOL_PRECOND_ILUDP  = DP

        integer(LINSOL_PRECOND_ILUINT) n, iwork(*),s,ierr,rwptr(*),colind(*)
        real(LINSOL_PRECOND_ILUDP)     a(*),relax
        integer(LINSOL_PRECOND_ILUINT) mneed,maxstr,nzlu,remain
        logical milu
        integer(LINSOL_PRECOND_ILUINT) lu,jlu,ilup
      end subroutine ilus
    end interface

      ! local variables
      integer(I32), dimension(:), pointer :: p_Idata
      integer(I32) :: maxstr,ierr,mneed

      ! Calculate a memory guess for how much memory the matrix needs.
      ! If the handle h_Idata is not associated to some memory block, then
      ! allocate a dummy block of size 1. Otherwise, try to reuse the existing
      ! memory block first and perform reallocation only if mandatory.
      if (h_Idata .eq. ST_NOHANDLE) then
        call storage_new('solver_initILU', 'p_Idata', int(1,I32),&
                         ST_INT, h_Idata, ST_NEWBLOCK_NOINIT)
      end if
      call storage_getbase_int(h_Idata, p_Idata)
      maxstr = 0
      Da     = 1

      ! Calculate the (M)ILU(s) matrix using SPLIB.
      call ILUS(neq, Da, Kcol, Kld, ifill, domega,&
                lu, jlu, ilup, p_Idata, maxstr, ierr, mneed)
      maxstr = max(mneed, 3*NA+(3+4*ifill)*NEQ)+10000

      if (ierr .ne. 0) then
        try: do

          ! Reallocate the memory
          call storage_realloc('solver_initILU', int(maxstr, I32),&
                                h_Idata, ST_NEWBLOCK_NOINIT, .false.)
          call storage_getbase_int(h_Idata, p_Idata)

          ! Calculate the ILU(s) matrix using SPLIB
          call ILUS(neq, Da, Kcol, Kld, ifill, domega,&
                    lu, jlu, ilup, p_Idata, maxstr, ierr, mneed)

          ! Error?
          select case(ierr)
          case (:-1)
            call output_line('Singular matrix!',&
                             OU_CLASS_ERROR,OU_MODE_STD,'solver_initILU')
            call sys_halt()

          case (0)
            ! ok
            exit try

          case (1)
            ! Insufficient memory
            maxstr = maxstr+max(mneed, NEQ)

          case DEFAULT
            ! Unknown error
            call output_line('Internal error!',&
                             OU_CLASS_ERROR,OU_MODE_STD,'solver_initILU')
            call sys_halt()
          end select
        end do try
      end if

      ! Check that not too much memory is wasted
      if ((h_Idata .ne. ST_NOHANDLE) .and. (mneed < maxstr/2)) then
        call storage_realloc('solver_initILU', int(mneed,I32),&
                             h_Idata, ST_NEWBLOCK_NOINIT, .true.)
      end if
    end subroutine do_mat79MILUs


    !**************************************************************
    ! Calculate (shifted, modified) ILU(0) decomposition,
    ! whereby the matrix is stored in format CSR7

    subroutine do_mat7MILU0(Da, Kcol, Kld, neq, nvar, ilu, alpha, tol)

      real(DP), intent(in) :: alpha,tol
      integer, dimension(:), intent(in) :: Kcol
      integer, dimension(:), intent(in) :: Kld
      integer, intent(in) :: neq
      integer, intent(in) :: nvar
      integer, intent(in) :: ilu

      real(DP), dimension(nvar,*), intent(inout) :: Da

      ! local variables
      real(DP), dimension(nvar) :: A,A_ij
      real(DP) :: alpha1
      integer :: ip,jc,ild,jld,icol,jcol,jcol0,ieq,ivar


      if (ilu .eq. 1) then

        ! Constant shifting for ILU
        if (alpha .ne. 0._DP .and. alpha .ne. -1._DP) then
          alpha1 = 1._DP/(1._DP+alpha)

!$omp  parallel do default(shared) private(ieq,ild)
          do ieq = 1, neq
            do ild = Kld(ieq)+1, Kld(ieq+1)-1
              Da(:,ild) = Da(:,ild)*alpha1
            end do
          end do
!$omp  end parallel do
        end if

      else

        ! Constant shift for MILUE
        if (alpha .ne. 0._DP) then
          alpha1 = 1._DP+alpha
!$omp  parallel do default(shared) private(ieq,ild)
          do ieq = 1, neq
            ild = Kld(ieq)
            Da(:,ild) = Da(:,ild)*alpha1
          end do
!$omp  end parallel do
        end if

      end if


      if (ilu .eq. 2) then
        !-----------------------------------------------------------------------
        ! Modified incomplete LU decomposition
        !-----------------------------------------------------------------------

        ! Incomplete Gauss elimination
        ild = 1
        outer1: do ieq = 2, neq
          A = Da(:,ild)

          do ivar = 1, nvar
            if (abs(A(ivar)) .lt. tol) Da(ivar,ild) = sign(tol, A(ivar))
          end do

          ild = Kld(ieq)
          if (Kcol(ild+1) .ge. ieq) cycle outer1

          inner1: do jcol = ild+1, Kld(ieq+1)-1
            icol = Kcol(jcol)
            if (icol .ge. ieq) cycle outer1
            jld = Kld(icol)

            ! Compute A_ij, that is, the entry of the eliminated matrix at position (IEQ,ICOL)
            A_ij = Da(:,jcol)/Da(:,jld)
            Da(:,jcol) = A_ij

            ! Store pointer
            ip = Kld(ieq+1)-1

            ! Loop over subdiagonal entries in line ICOL
            sub1: do jc = Kld(icol+1)-1, jld, -1
              jcol0 = Kcol(jc)

              ! First check for diagonal entry - diagonals are stored separately
              if (jcol0 .eq. ieq) then

                ! Insert at diagonal position
                Da(:,ild) = Da(:,ild)-A_ij*Da(:,jc)

              else

                if ((jcol0 .lt. ieq)  .and.&
                    (jcol0 .le. icol)) cycle inner1

                ! Look for an entry at position (ICOL,JCOL0)
                search1: do
                  if (Kcol(ip) .lt. jcol0) then

                    ! Insert at diagonal position
                    Da(:,ild) = Da(:,ild)-A_ij*Da(:,jc)

                    ! That is it exit searching
                    exit search1

                  elseif (Kcol(ip) .eq. jcol0) then

                    ! Off-diagonal entry (ICOL,JCOL0) found
                    Da(:,ip) = Da(:,ip)-A_ij*Da(:,jc)
                    ip = ip-1

                    ! That is it exit searching
                    exit search1

                  else

                    ! Entry not yet found
                    ip = ip-1
                    if (ip .gt. ild) cycle search1

                    ! Insert at diagonal position
                    Da(:,ild) = Da(:,ild)-A_ij*Da(:,jc)

                    ! That is it exit searching
                    exit search1
                  end if
                end do search1
              end if
            end do sub1
          end do inner1
        end do outer1

      else

        !-----------------------------------------------------------------------
        ! Shifted incomplete LU decomposition
        !-----------------------------------------------------------------------

        ! Incomplete Gauss elimination
        ild = 1
        outer2: do ieq = 2, neq
          A = Da(:,ild)

          do ivar = 1, nvar
            if (abs(A(ivar)) .lt. tol) Da(ivar,ild) = sign(tol, A(ivar))
          end do

          ild = Kld(ieq)
          if (Kcol(ild+1) .ge. ieq) cycle outer2

          inner2: do jcol = ild+1, Kld(ieq+1)-1
            icol = Kcol(jcol)
            if (icol .ge. ieq) cycle outer2
            jld = Kld(icol)

            ! Compute A_ij, that is, the entry of the eliminated matrix at position (IEQ,ICOL)
            A_ij = Da(:,jcol)/Da(:,jld)
            Da(:,jcol) = A_ij

            ! Store pointer
            ip = Kld(ieq+1)-1

            ! Loop over subdiagonal entries in line ICOL
            sub2: do jc = Kld(icol+1)-1, jld, -1
              jcol0 = Kcol(jc)

              ! First check for diagonal entry - diagonals are stored separately
              if (jcol0 .eq. ieq) then

                ! Insert at diagonal position
                Da(:,ild) = Da(:,ild)-A_ij*Da(:,jc)

              else

                if ((jcol0 .lt. ieq)  .and.&
                    (jcol0 .le. icol)) cycle inner2

                ! Look for an entry at position (ICOL,JCOL0)
                search2: do
                  if (Kcol(ip) .eq. jcol0) then

                    ! Off-diagonal entry (ICOL,JCOL0) found
                    Da(:,ip) = Da(:,ip)-A_ij*Da(:,jc)
                    ip = ip-1

                    ! That is it, exit searching
                    exit search2

                  elseif (Kcol(ip) .gt. jcol0) then
                    ip = ip-1

                  else

                    ! That is it, exit searching
                    exit search2
                  end if
                end do search2
              end if
            end do sub2
          end do inner2
        end do outer2
      end if
    end subroutine do_mat7MILU0


    !**************************************************************
    ! Calculate (shifted, modified) ILU(0) decomposition,
    ! whereby the matrix is stored in format CSR9

    subroutine do_mat9MILU0(Da, Kcol, Kld, Kdiagonal, neq, nvar, ilu, alpha, tol)

      real(DP), intent(in) :: alpha,tol
      integer, dimension(:), intent(in) :: Kcol
      integer, dimension(:), intent(in) :: Kld
      integer, dimension(:), intent(in) :: Kdiagonal
      integer, intent(in) :: neq
      integer, intent(in) :: nvar
      integer, intent(in) :: ilu

      real(DP), dimension(nvar,*), intent(inout) :: Da

      ! local variables
      real(DP), dimension(nvar) :: A,A_ij
      real(DP) :: alpha1
      integer :: ip,jc,ild,jld,icol,jcol,jcol0,ieq,ivar


      if (ilu .eq. 1) then

        ! Constant shifting for ILU
        if (alpha .ne. 0._DP .and. alpha .ne. -1._DP) then
          alpha1 = 1._DP/(1._DP+alpha)

!$omp  parallel do default(shared) private(ieq,ild)
          do ieq = 1,neq
            do ild = Kld(ieq), Kdiagonal(ieq)-1
              Da(:,ild) = Da(:,ild)*alpha1
            end do
            do ild = Kdiagonal(ieq)+1, Kld(ieq+1)-1
              Da(:,ild) = Da(:,ild)*alpha1
            end do
          end do
!$omp  end parallel do
        end if

      else

        ! Constant shift for MILU
        if (alpha .ne. 0._DP) then
          alpha1 = 1._DP+alpha
!$omp  parallel do default(shared) private(ieq,ild)
          do ieq = 1, neq
            ild = Kdiagonal(ieq)
            Da(:,ild) = Da(:,ild)*alpha1
          end do
!$omp  end parallel do
        end if

      end if


      if (ilu .eq. 2) then
        !-----------------------------------------------------------------------
        ! Modified incomplete LU decomposition
        !-----------------------------------------------------------------------

        ! Incomplete Gauss elimination
        ild = Kdiagonal(1)
        outer1: do ieq = 2, neq
          A = Da(:,ild)

          do ivar = 1, nvar
            if (abs(A(ivar)) .lt. tol) Da(ivar,ild) = sign(tol, A(ivar))
          end do

          ild = Kdiagonal(ieq)
          if (Kcol(Kld(ieq)) .ge. ieq) cycle outer1

          inner1: do jcol = Kld(ieq), Kld(ieq+1)-1
            icol = Kcol(jcol)
            if (icol .ge. ieq) cycle outer1
            jld = Kdiagonal(icol)

            ! Compute A_ij, that is, the entry of the eliminated matrix at position (IEQ,ICOL)
            A_ij = Da(:,jcol)/Da(:,jld)
            Da(:,jcol) = A_ij

            ! Store pointer
            ip = Kld(ieq+1)-1

            ! Loop over subdiagonal entries in line ICOL
            sub1: do jc = Kld(icol+1)-1, jld, -1
              jcol0 = Kcol(jc)

              ! First check for diagonal entry - diagonals are stored separately
              if (jcol0 .eq. ieq) then

                ! Insert at diagonal position
                Da(:,ild) = Da(:,ild)-A_ij*Da(:,jc)

              else

                if ((jcol0 .lt. ieq)  .and.&
                    (jcol0 .le. icol)) cycle inner1

                ! Look for an entry at position (ICOL,JCOL0)
                search1: do
                  if (Kcol(ip) .lt. jcol0) then

                    ! Insert at diagonal position
                    Da(:,ild) = Da(:,ild)-A_ij*Da(:,jc)

                    ! That is it exit searching
                    exit search1

                  elseif (Kcol(ip) .eq. jcol0) then

                    ! Off-diagonal entry (ICOL,JCOL0) found
                    Da(:,ip) = Da(:,ip)-A_ij*Da(:,jc)
                    ip = ip-1

                    ! That is it exit searching
                    exit search1

                  else

                    ! Entry not yet found
                    ip = ip-1
                    if (ip .gt. ild) cycle search1

                    ! Insert at diagonal position
                    Da(:,ild) = Da(:,ild)-A_ij*Da(:,jc)

                    ! That is it exit searching
                    exit search1
                  end if
                end do search1
              end if
            end do sub1
          end do inner1
        end do outer1

      else

        !-----------------------------------------------------------------------
        ! Shifted incomplete LU decomposition
        !-----------------------------------------------------------------------

        ! Incomplete Gauss elimination
        ild = Kdiagonal(1)
        outer2: do ieq = 2, neq
          A = Da(:,ild)

          do ivar = 1, nvar
            if (abs(A(ivar)) .lt. tol) Da(ivar,ild) = sign(tol, A(ivar))
          end do

          ild = Kdiagonal(ieq)
          if (Kcol(Kld(ieq)) .ge. ieq) cycle outer2

          inner2: do jcol = Kld(ieq), Kld(ieq+1)-1
            icol = Kcol(jcol)
            if (icol .ge. ieq) cycle outer2
            jld = Kdiagonal(icol)

            ! Compute A_ij, that is, the entry of the eliminated matrix at position (IEQ,ICOL)
            A_ij = Da(:,jcol)/Da(:,jld)
            Da(:,jcol) = A_ij

            ! Store pointer
            ip = Kld(ieq+1)-1

            ! Loop over subdiagonal entries in line ICOL
            sub2: do jc = Kld(icol+1)-1, jld, -1
              jcol0 = Kcol(jc)

              ! First check for diagonal entry - diagonals are stored separately
              if (jcol0 .eq. ieq) then

                ! Insert at diagonal position
                Da(:,ild) = Da(:,ild)-A_ij*Da(:,jc)

              else

                if ((jcol0 .lt. ieq)  .and.&
                    (jcol0 .le. icol)) cycle inner2

                ! Look for an entry at position (ICOL,JCOL0)
                search2: do
                  if (Kcol(ip) .eq. jcol0) then

                    ! Off-diagonal entry (ICOL,JCOL0) found
                    Da(:,ip) = Da(:,ip)-A_ij*Da(:,jc)
                    ip = ip-1

                    ! That is it, exit searching
                    exit search2

                  elseif (Kcol(ip) .gt. jcol0) then
                    ip = ip-1

                  else

                    ! That is it, exit searching
                    exit search2
                  end if
                end do search2
              end if
            end do sub2
          end do inner2
        end do outer2
      end if
    end subroutine do_mat9MILU0


    !**************************************************************
    ! Calculate BILU(0) decomposition of a matrix in format CSR7/CSR9,
    ! whereby the local blocks are full matrices

    subroutine do_mat79Intl1BILU0(DDa, Da, Kdiagonal, neq, nvar)

      real(DP), dimension(nvar,nvar,*), intent(in) :: DDa
      real(DP), dimension(nvar,nvar,*), intent(inout) :: Da
      integer, dimension(:), intent(in) :: Kdiagonal
      integer, intent(in) :: neq
      integer, intent(in) :: nvar

      ! local variables
      real(DP) :: a_ij
      integer :: ild,ieq,ivar,jvar,kvar


      ! Loop over all diagonal blocks
      do ieq = 1, neq

        ! Get position of diagonal block
        ild = Kdiagonal(ieq)

        ! Copy data from original matrix
        Da(:,:,ieq) = DDa(:,:,ild)

        ! Perform complete LU decomposition (without pivoting)

        ! Loop over columns of Crout`s method
        do jvar = 1, nvar
          do ivar = 1, jvar-1
            a_ij = Da(ivar,jvar,ieq)
            do kvar = 1, ivar-1
              a_ij = a_ij-Da(ivar,kvar,ieq)*Da(kvar,jvar,ieq)
            end do
            Da(ivar,jvar,ieq) = a_ij
          end do


          do ivar = jvar, nvar
            a_ij = Da(ivar,jvar,ieq)
            do kvar = 1, jvar-1
              a_ij = a_ij-Da(ivar,kvar,ieq)*Da(kvar,jvar,ieq)
            end do
            Da(ivar,jvar,ieq) = a_ij
          end do

          ! Divide by the diagonal element
          do ivar = jvar+1, nvar
            Da(ivar,jvar,ieq) = Da(ivar,jvar,ieq)/Da(jvar,jvar,ieq)
          end do
        end do
      end do
    end subroutine do_mat79Intl1BILU0
  end subroutine solver_initILU

  ! ***************************************************************************

!<subroutine>

  subroutine solver_initUMFPACK(rsolver)

!<description>
    ! This subroutine initializes the UMFPACK subsystem
!</description>

!<inputoutput>
    ! Umfpack-solver structure
    type(t_solverUMFPACK), intent(inout) :: rsolver
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_matrixScalar) :: rmatrixScalar
    integer, dimension(:), pointer :: p_Kcol
    integer, dimension(:), pointer :: p_Kld
    real(DP), dimension(:), pointer :: p_Da


    ! Status variables of UMFPACK4; receives the UMFPACK-specific return code
    ! of a call to the solver routines.
    real(DP), dimension(90) :: Dinfo

    ! Check if matrix is associated
    if ((rsolver%rmatrix%nblocksPerCol .eq. 0) .or.&
        (rsolver%rmatrix%nblocksPerRow .eq. 0)) then
      call output_line('Matrix must be associated to solver node!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'solver_initUMFPACK')
      call sys_halt()
    end if

    ! Prepare the system matrix for UMFPACK and convert it into a scalar matrix
    call initPrepareMatrix(rsolver%rmatrix, rmatrixScalar)

    ! Set pointers to row- and column indices and data array
    call lsyssc_getbase_Kld   (rmatrixScalar, p_Kld)
    call lsyssc_getbase_Kcol  (rmatrixScalar, p_Kcol)
    call lsyssc_getbase_double(rmatrixScalar, p_Da)


    ! Perform symbolical factorisation
    if (rsolver%isymbolic .ne. 0) call UMF4FSYM(rsolver%isymbolic)
    call UMF4SYM(rmatrixScalar%NEQ, rmatrixScalar%NCOLS,&
                 p_Kld, p_Kcol, p_Da, rsolver%isymbolic,&
                 rsolver%Dcontrol, Dinfo)

    ! Check Dinfo(1) if there is an error
    select case(int(Dinfo(1)))
    case (0)
      ! ok.

    case (1)
      ! Singular matrix
      call output_line('Matrix is singular!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'solver_initUMFPACK')
      call sys_halt()

    case (-1)
      ! Insufficient memory
      call output_line('Insufficient memory!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'solver_initUMFPACK')
      call sys_halt()

    case DEFAULT
      ! Internal error
      call output_line('Internal UMFPACK error!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'solver_initUMFPACK')
      call sys_halt()
    end select


    ! Perform numerical factorisation
    if (rsolver%inumeric .ne. 0) call UMF4FNUM(rsolver%inumeric)
    call UMF4NUM(p_Kld, p_Kcol, p_Da, rsolver%isymbolic,&
                 rsolver%inumeric, rsolver%Dcontrol, Dinfo)

    ! Check Dinfo(1) if there is an error
    select case(int(Dinfo(1)))
    case (0)
      ! ok.

    case (1)
      ! Singular matrix
      call output_line('Matrix is singular!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'solver_initUMFPACK')
      call sys_halt()

    case (-1)
      ! Insufficient memory
      call output_line('Insufficient memory!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'solver_initUMFPACK')
      call sys_halt()

    case DEFAULT
      ! Unknown error
      call output_line('Internal UMFPACK error!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'solver_initUMFPACK')
      call sys_halt()
    end select

    ! Release temporary matrix
    call lsyssc_releaseMatrix(rmatrixScalar)

  contains

    ! Here, the real working routine follows

    !*************************************************************
    ! Prepare the matrix for UMFACK

    subroutine initPrepareMatrix(rmatrixSrc, rmatrixDest)
      type(t_matrixBlock), intent(in) :: rmatrixSrc
      type(t_matrixScalar), intent(inout) :: rmatrixDest


      ! Check if matrix is a 1x1 block matrix
      if ((rmatrixSrc%nblocksPerCol .eq. 1) .and.&
          (rmatrixSrc%nblocksPerRow .eq. 1)) then

        ! Initialize working matrix
        select case(rmatrixSrc%RmatrixBlock(1,1)%cmatrixFormat)

        case (LSYSSC_MATRIX9)
          ! Format 9 is exactly the UMFPACK matrix.
          ! Make a copy of the matrix structure, but use the same matrix entries.
          call lsyssc_duplicateMatrix(rmatrixSrc%RmatrixBlock(1,1), rmatrixDest,&
                                      LSYSSC_DUP_COPY, LSYSSC_DUP_SHARE)

        case (LSYSSC_MATRIX7)
          ! For format 7, we have to modify the matrix slightly.
          ! Make a copy of the complete matrix
          call lsyssc_duplicateMatrix(rmatrixSrc%RmatrixBlock(1,1), rmatrixDest,&
                                      LSYSSC_DUP_COPY, LSYSSC_DUP_COPY)

          ! Convert to format 9
          call lsyssc_convertMatrix(rmatrixDest, LSYSSC_MATRIX9, .true.)

        case DEFAULT
          call output_line('Unsupported matrix format!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'initPrepareMatrix')
          call sys_halt()
        end select

        ! Since UMFPACK4 is written in C, all arrays start at position 0 instead of 1
        call lsyssc_addIndex(rmatrixDest%h_Kld,  -1_I32, ilength=rmatrixDest%NEQ+1)
        call lsyssc_addIndex(rmatrixDest%h_Kcol, -1_I32, ilength=rmatrixDest%NA)

      else

        call output_line('Block matrices are not supported!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'initPrepareMatrix')
        call sys_halt()

      end if

    end subroutine initPrepareMatrix

  end subroutine solver_initUMFPACK

end module solveraux
