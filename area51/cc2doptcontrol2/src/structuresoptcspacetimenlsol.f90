!##############################################################################
!# ****************************************************************************
!# <name> structuresoptcspacetimenlsol </name>
!# ****************************************************************************
!#
!# <purpose>
!# Underlying parameter structures of the nonlinear space-time solver.
!#
!# Routines in this module:
!#
!# 1.) stnlsinit_initSolver
!#     -> Initialises a nonlinear solver based on a nonlinear solver parameter
!#        structure
!#
!# 2.) stnlsinit_doneSolver
!#     -> Releases a nonlinear solver.
!# </purpose>
!##############################################################################

module structuresoptcspacetimenlsol

  use fsystem
  use genoutput
  use storage
  use boundary
  use triangulation
  use paramlist
  use discretebc
  use discretefbc
  use fparser
  use linearsystemscalar
  use statistics
  
  use collection
  
  use timediscretisation
  use fespacehierarchybase
  use fespacehierarchy
  use spacetimevectors
  use spacetimehierarchy

  use structuresspacetimelinsol
  use structuresoptflow
  use spacetimelinearsystem
  use spacetimelinearsolver
  
  use spacetimeneumannbc
  
  implicit none
  
  private
  
!<constants>

!<constantblock description="Identifiers for the different types of nonlinear solvers.">

  ! No solver, dummy
  integer, parameter, public :: CCNLS_NONE         = 0

  ! Preconditioning by linear solver, solving the linearised system
  integer, parameter, public :: CCNLS_LINEARSOLVER  = 1

  ! Preconditioning by Newton-Iteration
  integer, parameter, public :: CCNLS_NEWTON        = 2
  
  ! Preconditioning by inexact/adaptive Newton iteration
  integer, parameter, public :: CCNLS_INEXACTNEWTON = 3

!</constantblock>

!</constants>

!<types>

!<typeblock>

  ! Collects all parameters of the nonlinear space-time solver
  type t_nlstsolver
  
    ! Type of preconditioner to use in every timestep, = type of linear
    ! subproblem to solve in the nonlinear space time defect correction loop.
    ! One of the CCNLS_xxxx constants.
    ! =1: Preconditioning by solving the standard system.
    ! =2: Preconditioning by Newton-Iteration.
    !     This is a combination of a linear solver (itypePreconditioning=1)
    !     in conjunction with an extended matrix to improve convergence speed
    !     (at the cost of assembly time and stability of the linear solver/
    !     preconditioner).
    ! =3: Inexact Newton. Stopping criterion of the linear solver is chosen
    !     adaptively. The stopping criteria of the linear solvers are not used.

    integer :: ctypeNonlinearIteration = CCNLS_NEWTON

    ! Whether to postprocess intermediate solutions.
    ! =1: After each nonlinear step, apply postprocessing to the solution.

    integer :: cpostprocessIterates = 0

    ! Minimum number of steps

    integer :: nminIterations     = 1

    ! Maximum number of steps
          
    integer :: nmaxIterations     = 10

    ! Damping of residuals, i.e. reduction of relative error
    ! on finest grid; smaller -> more iterations

    real(dp) :: depsRel            = 1E-8

    ! Limit for residuals, i.e. absolute error on finest grid;
    ! The linear solver stops if both, absolute error < depsAbs and
    ! rel. error < depsRel

    real(dp) :: depsAbs            = 1E-0

    ! General damping parameter when adding the preconditioned defect
    ! to the solution in the nonlinear space-time iteration

    real(dp) :: domega             = 1.0

    ! Limit for differences in the residuals
    real(dp) :: depsDiff           = 1E-5

    ! Output level of the solver

    integer :: ioutputLevel       = 2

    ! Maximum number of defect correction iterations before starting Newton.
    ! Standard=0=use Newton / inexact Newton immediately.

    integer :: nmaxFixedPointIterations = 0

    ! Relative convergence criterion to use for defect correction
    ! before Newton is started.

    real(DP) :: depsRelFixedPoint = 1E-1

    ! Stopping criterion for the linear solver in the inexact Newton iteration.
    ! Controls the minimum number of digits to gain in the linear solver
    ! per Newton iteration. Only used if ctypePreconditioner = 3.

    real(dp) :: dinexactNewtonEpsRel = 1.0E-2

    ! Exponent to control the stopping criterion for the linear solver in
    ! an inexact Newton iteration. =2 result in quadratic convergence,
    ! =1.5 in superlinear convergence. Only used if ctypePreconditioner = 3.

    real(dp) :: dinexactNewtonExponent = 2.0
  
    ! <!-- Parameters automatically maintained by the solver -->
    
    ! Reference to the problem structure.
    type(t_settings_optflow), pointer :: p_rsettings => null()
    
    ! Space-time solver node for the linear space-time preconditioner.
    ! Points either to p_rmgSolver or p_rsgSolver.
    type(t_sptilsNode), pointer :: p_rspaceTimePrec => null()
    
    ! Space-time solver node of an MG solver -- of present.
    type(t_sptilsNode), pointer :: p_rmgSolver => null()

    ! Space-time solver node of a coarse grid solver if a multigrid solver
    ! is involved.
    type(t_sptilsNode), pointer :: p_rcgrSolver => null()
    
    ! Nonlinear space-time matrix.
    !type(t_ccoptSpaceTimeMatrix), pointer :: p_rmatrix => null()
    
    ! Preconditioner matrices on all levels
    type(t_ccoptSpaceTimeMatrix), dimension(:), pointer :: p_RprecMatrices => null()
  
    ! Evaluation points of the nonlinearities on all levels
    type(t_spacetimeVector), dimension(:), pointer :: p_Rsolutions => null()

    ! Neumann boundary conditions on all levels; for nonlinear boundarty conditions
    type(t_sptiNeumannBoundary), dimension(:), pointer :: p_rsptiNeumannBC => null()
  
    ! <!-- Output -->

    ! Number of linear iterations of all spatial preconditioners.
    integer :: nlinearIterationsSpace = 0

    ! Number of linear iterations
    integer :: nlinearIterations = 0

    ! Number of nonlinear iterations
    integer :: nnonlinearIterations = 0
    
    ! Number of iterations on the coarse grid
    integer :: ncoarsegridIterations = 0
    
    ! STATISTICS: Total time for nonlinear iteration
    type(t_timer) :: rtimerNonlinear
    
    ! STATISTICS: Total time for linear preconditioning
    type(t_timer) :: rtimerPreconditioner
  
    ! STATISTICS: Total time needed for smoothing operations
    type(t_timer) :: rtimeSmoothing
        
    ! STATISTICS: Total time needed for the coarse grid solver
    type(t_timer) :: rtimeCoarseGridSolver
    
    ! STATISTICS: Time needed for linear algebra stuff (matrix-vector,
    ! vector-copy, prolongation/restriction,...)
    type(t_timer) :: rtimeLinearAlgebra
    
    ! STATISTICS: Time needed for prolongation/restriction
    type(t_timer) :: rtimeProlRest
    
    ! STATISTICS: Time for solving problems in space.
    type(t_timer) :: rtimeSpacePrecond

    ! STATISTICS: Time for initialisation / factorisation of the space time system;
    ! global and in one step
    type(t_timer) :: rtimeFactorisation

    ! STATISTICS: Time for postprocessing
    type(t_timer) :: rtimePostprocessing
    
  end type

!</typeblock>

  public :: t_nlstsolver
  public :: stnlsinit_initSolver
  public :: stnlsinit_doneSolver
  public :: stnlsinit_printSolverStatistics

contains

  ! ***************************************************************************
  
!<subroutine>

  subroutine stnlsinit_initSolver (rsettings,ispaceTimeLevel,rsettingsPrecond,rsolver)

!<description>
  ! Creates a nonlinear solver node rsolversettings based on the settings in
  ! rsettings and rprecsettings
!</description>

!<input>
  ! Global settings structure.
  type(t_settings_optflow), intent(inout), target :: rsettings
  
  ! Absolute level of the solver (relative to the global hierarchies in
  ! rsettings).
  integer, intent(in) :: ispaceTimeLevel

  ! Parameters configuring the linear space-time preconditioner
  type(t_settings_nlstprec), intent(in) :: rsettingsPrecond
!</input>

!<output>
  ! A solver node representing a nonlinear space-time solver.
  type(t_nlstsolver), intent(out) :: rsolver
!</output>

!</subroutine>

    ! local variables
    integer :: ilev,ispaceLevel
    type(t_feSpaceLevel), pointer :: p_rfeSpaceLevel
    type(t_timeDiscretisation), pointer :: p_rtimeDiscr

    ! Allocate memory for the preconditioner matrices/vectors on all levels
    allocate(rsolver%p_RprecMatrices(ispaceTimeLevel))
    allocate(rsolver%p_Rsolutions(ispaceTimeLevel))
    allocate(rsolver%p_rsptiNeumannBC(ispaceTimeLevel))

    ! Allocate temp vectors on all levels for the nonlinearity
    do ilev = 1,ispaceTimeLevel
    
      ! Get the level
      call sth_getLevel (rsettings%rspaceTimeHierPrimalDual,ilev,&
          p_rfeSpaceLevel,p_rtimeDiscr,ispaceLevel)
      
      ! Create the vector.
      call sptivec_initVector (rsolver%p_Rsolutions(ilev),&
          p_rtimeDiscr,p_rfeSpaceLevel%p_rdiscretisation)
          
      ! Prepare arrays for the Neumann boundary conditions.
      ! Used for nonlinear boundary conditions.
      call stnm_createNeumannBoundary (&
          p_rfeSpaceLevel%p_rdiscretisation,p_rtimeDiscr,&
          rsettings%rspaceAsmHierarchy%p_RasmTemplList(ispaceLevel),&
          rsolver%p_rsptiNeumannBC(ilev))
    
    end do
        
    ! Get the linear preconditioner.
    call stlsinit_getSolver (rsettings,ispaceTimeLevel,rsettingsPrecond,&
        rsolver%p_rspaceTimePrec,rsolver%p_rmgSolver,rsolver%p_rcgrSolver)
    
    ! Initialise projection hierarchies for multilevel algorithms.
    if (associated (rsolver%p_rmgSolver)) then
      call sptils_setHierarchyMultigrid (rsolver%p_rmgSolver,&
          rsettings%rspaceTimeHierPrimalDual,rsettings%rprjHierSpaceTimePrimalDual)
    end if
    
  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine stnlsinit_doneSolver (rsolver)

!<description>
  ! Cleans up the nonlinear space-time solver.
!</description>

!<inputoutput>
  ! Nonlinear solver structure to clean up.
  type(t_nlstsolver), intent(inout) :: rsolver
!</inputoutput>

!</subroutine>
    ! local variables
    integer :: ilev

    ! Release the linear preconditioner.
    call sptils_releaseSolver(rsolver%p_rspaceTimePrec)
    
    ! Other pointers can be set to null() as they are part of p_rspaceTimePrec
    ! and automatically released.
    rsolver%p_rmgSolver => null()
    rsolver%p_rcgrSolver => null()

    ! Release temp vectors and Neumann BC's.
    do ilev = 1,size(rsolver%p_Rsolutions)
      call stnm_releaseNeumannBoundary (rsolver%p_rsptiNeumannBC(ilev))
      call sptivec_releaseVector(rsolver%p_Rsolutions(ilev))
    end do

    ! Release preconditioner matrices
    deallocate(rsolver%p_Rsolutions)
    deallocate(rsolver%p_rsptiNeumannBC)
    deallocate(rsolver%p_RprecMatrices)

  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine stnlsinit_printSolverStatistics (rsolver)

!<description>
  ! Prints statistics about the solver to the terminal
!</description>

!<inputoutput>
  ! Nonlinear solver structure
  type(t_nlstsolver), intent(in) :: rsolver
!</inputoutput>

!</subroutine>

    ! Print some statistical data
    call output_line ("Nonlinear solver statistics:")
    call output_lbrk()
    call output_line ('Total time for nonlinear iteration = '// &
        trim(sys_sdL(rsolver%rtimerNonlinear%delapsedReal,10)))
    call output_line ('Total time for factorisation       = '// &
        trim(sys_sdL(rsolver%rtimeFactorisation%delapsedReal,10)))
    call output_line ('Total time for preconditioning     = '// &
        trim(sys_sdL(rsolver%rtimerPreconditioner%delapsedReal,10)))
    call output_line ('#nonlinear iterations              = '//&
        trim(sys_siL(rsolver%nnonlinearIterations,10)))
    call output_line ('#iterations preconditioner         = '//&
        trim(sys_siL(rsolver%nlinearIterations,10)))
    call output_line ('#iterations coarse-grid solver     = '//&
        trim(sys_siL(rsolver%ncoarsegridIterations,10)))
    call output_line ('#iterations space-preconditioners  = '//&
        trim(sys_siL(rsolver%nlinearIterationsSpace,10)))
    call output_line ('Total time for postprocessing      = '// &
        trim(sys_sdL(rsolver%rtimePostprocessing%delapsedReal,10)))

    call output_separator (OU_SEP_MINUS)

    call output_line ("Preconditioner statistics:")
    call output_lbrk()
    call output_line ("Total time for non-c.grid solving  = "//&
        sys_sdL(rsolver%rtimerNonlinear%delapsedReal - &
                rsolver%rtimeCoarseGridSolver%delapsedReal - &
                rsolver%rtimePostprocessing%delapsedReal,10))
    call output_line ("Total time for smoothing           = "//&
        sys_sdL(rsolver%rtimeSmoothing%delapsedReal,10))
    call output_line ("Total time for coarse grid solving = "//&
        sys_sdL(rsolver%rtimeCoarseGridSolver%delapsedReal,10))
    call output_line ("Total time for linear algebra      = "//&
        sys_sdL(rsolver%rtimeLinearAlgebra%delapsedReal,10))
    call output_line ("Total time for prol/rest           = "//&
        sys_sdL(rsolver%rtimeProlRest%delapsedReal,10))
    call output_lbrk()
    call output_line ("Total time for prec. in space      = "//&
        sys_sdL(rsolver%rtimeSpacePrecond%delapsedReal,10))

  end subroutine

end module
