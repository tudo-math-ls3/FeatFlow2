!##############################################################################
!# ****************************************************************************
!# <name> CHnonlinearcore </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains the routines to solve the stationary core equation
!# of the problem with a nonlinear solver.
!# (This is comparable to the NSDEF routine in older FEAT versions)
!#
!# The discretised core equation reads at the moment:
!#
!#  $$        A \phi  +   B w   = f_1 $$
!#  $$        D \phi  +   C w   = f_2 $$
!#
!# with
!#
!#   $$ A = \alpha M  +  \gamma Conv $$
!#   $$ B = \eta Lap W               $$
!#   $$ C = \beta M                  $$
!#   $$ D = \tau N(c) + \theta Lap c $$
!#
!# and
!#   $\alpha$ = 0/1     - switches the mass matrix on/off;
!#   $\gamma$ = 0/1     - Switches the convective term on/off;
!#                          =0 for stationary problem,
!#   $\eta$   = 0/1     - Switches the 'B'-term on/off,
!#   $\theta$           - weight for the Laplace matrix,
!#   $\tau$   = 0/1     - Switches the 'B^T'-term on/off,
!#   $\beta$  =0/1      - In front of C, switches the 'C' -term on/off
!#
!# (y,p) is the phase variable/chemPotential solution pair.
!#
!# The core equation is abstractly written a nonlinear system of the form
!#
!#  $$ A(x)x = b $$
!#
!# and solved with the nonlinear solver from the kernel, using the defect
!# correction approach
!#
!#  $$  x_{n+1}  =  x_n  +  \omega_n C^{-1} ( b - A(x_n) x_n )  $$
!#
!# where $C^{-1}$ means to apply a suitable preconditioner (inverse mass
!# matrix, apply the linearised $A(x_n)^-1$ with multigrid, apply Newton or 
!# do something similar). 
!#
!# The following routines can be found here:
!#
!#  1.) CH_getNonlinearSolver
!#      -> Initialise nonlinear solver configuration with information 
!#         from INI/DAT files
!#
!#  2) CH_solveCoreEquation
!#      -> Starts the nonlinear iteration to solve the core equation.
!#
!# callback routines for the nonlinear solver:
!#
!#  6.) CH_getDefect
!#      -> callback routine. Calculate nonlinear defect
!#
!#  7.) CH_getOptimalDamping
!#     -> Auxiliary routine. Calculate optimal damping parameter
!#
!#  8.) CH_precondDefect
!#      -> callback routine. Preconditioning of nonlinear defect
!#
!#  9.) CH_resNormCheck
!#      -> callback routine: Check the residuals for convergence
!#
!# 10.) CH_nonlinearSolver
!#      -> The actual nonlinear solver; similar to NSDEF in old CH versions.
!#
!# 11.) CH_initPreconditioner  <---> CH_releasePreconditioner
!#     -> Initialize preconditioner.
!#
!# 12.) 
!#
!# To solve a system with the core equation, one has to deal with two
!# main structures. On one hand, one has a nonlinear iteration structure 
!# t_nlsolNode from the kernel; this is initialised by CH_getNonlinearSolver.
!# On the other hand, one has to maintain a 'nonlinear iteration structure'
!# of type t_CHNonlinearIteration, which configures the core equation and
!# specifies parameters for the solver how to work.
!#
!# The basic procedure for dealing with the core equation is as follows:
!#
!#  a) CH_getNonlinearSolver  -> Basic initialisation on the nonlinear solver
!#
!#  b) Basic initialisation of the core equation structures
!#
!#  c) Initialise further parameters in the core equation structure manually
!#     (e.g. preconditioner, pointer to matrices, ...).
!#     It's important, that the 'outer' application initialises pointers to
!#     matrices, otherwise nothing will work!
!#     This all has to be done with the nonlinear-iteration-structure directly.
!#
!#  d) Initialises the constants in the core equation structure according to
!#     the concrete situation.
!#
!#  e) CH_solveCoreEquation -> Solve the core equation with the nonlinear
!#                               solver structure from CH_getNonlinearSolver.
!#
!#  f) Release core equation structure
!#
!# Stebs b), c), d) and f) are done in the routines in the file 
!# CHnonlinearcoreinit.f90. This file here contains only the very core
!# routines of the nonlinear iteration, leaving all initialisation to
!# CHnonlinearcoreinit.f90.
!#
!# </purpose>
!##############################################################################

module CahnHilliard_nonlinearcore

  use fsystem
  use storage
  use linearsolver
  use boundary
  use linearalgebra
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
  use trilinearformevaluation
  use matrixio
  use statistics
  use collection
  use convection
  use linearsystemscalar
  use linearsystemblock
  
  use CahnHilliard_matvec
  use CahnHilliard_callback
  
!~~~~~~~~~~~~use NS problem~~~~~~~~~~~~~~~~~~
! use ccbasic
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  implicit none
  
!<constants>

!<constantblock description="Preconditioner identifiers for the defect in the nonlinear iteration">

  ! No preconditioning
  integer, parameter :: CHPREC_NONE         = -1

  ! Preconditioning with inverse mass matrix (not yet implemented)
  integer, parameter :: CHPREC_INVERSEMASS   = 0

  ! Preconditioning by linear solver, solving the linearised system
  integer, parameter :: CHPREC_LINEARSOLVER  = 1

  ! Preconditioning by Newton-Iteration
  integer, parameter :: CHPREC_NEWTON        = 2

  ! Preconditioning by dynamic Newton-Iteration (uses defect correction
  ! and switches automatically to Newton if the error is small enough)
  integer, parameter :: CHPREC_NEWTONDYNAMIC = 3

!</constantblock>
!</constants>
  
!<types>

!<typeblock>

  ! This structure controls the Newton iteration -- i.e. the preconditioning
  ! with the Frechet derivative of the CH equation, which
  ! can lead to quadratic covergence of the nonlinear solver.
  ! As Newton works only in the basin of attraction of the solution,
  ! the parameters in this structure allow to define a switching criterion
  ! when to use Newton. In the first couple of iterations, defect correction
  ! is used, while the iteration switches to Newton if the residuum is small
  ! enough.
  ! This block is used if CHPREC_NEWTONDYNAMIC is used as preconditioner.
  type t_CHDynamicNewtonControl
  
    ! Minimum number of fix point iterations before switching to
    ! preconditioning with the Newton matrix. (IFIXMIN)

    integer :: nminFixPointIterations = 0

    ! Maximum number of fix point iterations before switching to
    ! preconditioning with the Newton matrix. (IFIXMAX)

    integer :: nmaxFixPointIterations = 999

    ! Norm of absolute residuum before applying Newton. 
    ! Newton is only applied
    ! if   ||absolute residuum|| < depsAbsNewton
    ! and  ||relative residuum|| < depsRelNewton.
    ! Otherwise, the usual fix point iteration is used.
    ! Stamndard value = 1E-5.

    real(DP) :: depsAbsNewton = 1.0E-5_DP

    ! Norm of relative residuum before applying Newton. 
    ! Newton is only applied
    ! if   ||absolute residuum|| < depsAbsNewton
    ! and  ||relative residuum|| < depsRelNewton.
    ! Otherwise, the usual fix point iteration is used.
    ! Standard value = 1E99 -> The absolute residuum counts.

    real(DP) :: depsRelNewton = 1.0E99_DP
  
    ! Whether to use the inexact Newton iteration or not.
    ! The inexact Newton controls the stopping criterion of the linear
    ! solver according to the nonlinear residual.
    
    integer :: cinexactNewton = 1

    ! Stopping criterion for the linear solver in the inexact Newton iteration.
    ! Controls the minimum number of digits to gain in the linear solver
    ! per Newton iteration. Only used if cinexactNewton = 1.
    
    real(dp) :: dinexactNewtonEpsRel = 1.0E-2_DP

    ! Exponent to control the stopping criterion for the linear solver in
    ! an inexact Newton iteration. =2 result in quadratic convergence,
    ! =1.5 in superlinear convergence. Only used if cinexactNewton = 1.
    
    real(dp) :: dinexactNewtonExponent = 2.0_DP
  
  end type
  
!</typeblock>


!<typeblock>

  ! Preconditioner structure for CHxD. This structure saves the configuration of the
  ! preconditioner that is used during the nonlinear iteration. 
  
  type t_CHPreconditioner
  
    ! Type of preconditioner.
    ! This is one of the CHPREC_xxxx flags as defined above (CHPREC_INVERSEMASS for
    ! preconditioning with inverse mass matrix, CHPREC_LINEARSOLVER for solving a linear
    ! system, CHPREC_NEWTON for a Newton iteration,...)
    integer :: ctypePreconditioning = CHPREC_NONE
    
    ! Pointer to linear solver node if a linear solver is the preconditioner.
    ! (Thus, this applies for the defect correction and the Newton preconditioner).
    type(t_linsolNode), pointer :: p_rsolverNode => NULL()

    ! If the linerar solver contains a multigrid preconditioner, this is a pointer
    ! to the coarse grid solver. Otherwise, the pointer is not associated.
    type(t_linsolNode), pointer :: p_rcgrSolver => NULL()
    
    ! An interlevel projection structure for changing levels
    type(t_interlevelProjectionBlock), pointer :: p_rprojection => NULL()
    
    ! Configuration block for the adaptive Newton preconditioner.
    ! Is only valid if ctypePreconditioning=CHPREC_NEWTONDYNAMIC!
    type(t_CHDynamicNewtonControl) :: radaptiveNewton

    ! Temporary scalar vector; used for calculating the nonlinear matrix
    ! on lower levels / projecting the solution from higher to lower levels.
    type(t_vectorScalar), pointer :: p_rtempVectorSc => NULL()

    ! Temporary scalar vector; used for calculating the optimal damping
    ! parameter.
    type(t_vectorScalar), pointer :: p_rtempVectorSc2 => NULL()

  end type

!</typeblock>

!<typeblock>

  ! This type is used to save some situation specific assembly information
  ! during the setup phase of the nonlinear solver. Here it's noted, if
  ! and whose matrices exist and/or must be assmebled transposed to be
  ! compatible with the preconditioner and more. It's more or less
  ! a collection if different flags.
  type t_CHPreconditionerSpecials
  
    ! Whether to use 'adaptive matrices', i.e. set up coarse grid matrices
    ! with the help of fine grid matrices. This is used for very special
    ! discretisations only (e.g. Q1~/Q0). =0: deactivate
    integer :: iadaptiveMatrices    = 0
    
    ! A configuration parameter for adaptive matrices.
    real(DP) :: dadMatThreshold     = 0.0_DP

    ! If the preconditioner is a linear solver:
    ! Type of solver.
    ! =0: Gauss elimination (UMFPACK)
    ! =1: Multigrid solver
    integer :: isolverType = 0
    
    ! If the preconditioner is the linear multigrid solver:
    ! Type of smoother.
    ! =0: general VANCA (slow, but independent of the discretisation and of the problem)
    ! =1: general VANCA; 'direct' method, bypassing the defect correction approach.
    !     (-> specialised variant of 0, but slightly faster)
    ! =2: Simple Jacobi-like VANCA, 2D Cahn Hilliard problem, general discretisation
    !     (i.e. automatically chooses the best suitable VANCA variant).
    ! =3: Simple Jacobi-like VANCA, 2D Cahn Hilliard problem, general discretisation
    !     (i.e. automatically chooses the best suitable VANCA variant).
    !     'direct' method, bypassing the defect correction approach.
    !     (-> specialised variant of 8, but faster)
    ! =4: Full VANCA, 2D Cahn Hilliard problem, general discretisation
    !     (i.e. automatically chooses the best suitable VANCA variant).
    ! =5: Full VANCA, 2D Cahn Hilliard problem, general discretisation
    !     (i.e. automatically chooses the best suitable VANCA variant).
    !     'direct' method, bypassing the defect correction approach.
    !     (-> specialised variant of 10, but faster)
    integer :: ismootherType = 3
    
    ! If the preconditioner is the linear multigrid solver:
    ! Type of coarse grid solver.    
    ! =0: Gauss elimination (UMFPACK)
    ! =1: Defect correction with diagonal VANCA preconditioning.
    ! =2: BiCGStab with diagonal VANCA preconditioning
    integer :: icoarseGridSolverType = 1
        
  end type

!</typeblock>

!<typeblock>

  ! Represents the core equation on one level of the discretisation.
  ! Collects all information that are necessary to assemble the 
  ! (linearised) system matrix and RHS vector.
  type t_CHcoreEquationOneLevel
  
    ! Pointer to the discretisation structure of that level
    type(t_blockDiscretisation), pointer :: p_rdiscretisation => null()

    ! A matrix for that specific level (=nu*Laplace)
    type(t_matrixScalar), pointer :: p_rmatrixA => null()

    ! B-matrix for that specific level. 
    type(t_matrixScalar), pointer :: p_rmatrixB => null()

    ! C-matrix for that specific level. 
    type(t_matrixScalar), pointer :: p_rmatrixC => null()

    ! D-matrix for that specific level. 
    type(t_matrixScalar), pointer :: p_rmatrixD => null()

    ! Mass matrix
    type(t_matrixScalar), pointer :: p_rmatrixMass => null()

    ! Laplace matrix
    type(t_matrixScalar), pointer :: p_rmatrixLaplace => null()

    ! Convective matrix
    type(t_matrixScalar), pointer :: p_rmatrixConv => null()
    
    ! Temporary vector for the interpolation of a solution to a lower level.
    ! Exists only on levels NLMIN..NLMAX-1 !
    type(t_vectorBlock), pointer :: p_rtempVector => null()

    ! Block matrix, which is used in the defect correction / Newton
    ! algorithm as preconditioner matrix of the correspnding underlying
    ! linear sytem. Is usually the (linearise) system matrix or
    ! a Newton matrix. This matrix is changed during the
    ! nonlinear iteration and used e.g. if a linear solver (Multigrid) is
    ! used for preconditioning.
    type(t_matrixBlock), pointer :: p_rmatrixPreconditioner => null()

  end type

!</typeblock>

!<typeblock>

  ! This configuration block configures all parameters that are needed
  ! by the callback routines to perform the nonlinear iteration.
  ! On start of the program, a structure of this type is initialised.
  ! The entries of this structure are saved to the collection structure.
  ! When a callback routine is called, the structure is rebuild from
  ! the collection. When he nonlinear iteration is finished, the
  ! parameters are removed from the colletion again.
  type t_CHNonlinearIteration
  
    ! ALPHA-parameter that controls the weight of the mass matrix in the
    ! core equation. =0.0 for stationary simulations.
    real(DP) :: dalpha = 0.0_DP
    
	! GAMMA-parameter that controls the weight in front of the
    ! Convetive term...
    real(DP) :: dgamma = 0.0_DP

    ! ETA-parameter that switch the B-term on/off.
    real(DP) :: deta = 0.0_DP
    
    ! THETA-parameter that controls the weight of the A matrix
    ! in the core equation. =1.0 for stationary simulations.
    real(DP) :: dtheta = 0.0_DP
    
    ! TAU-parameter that switch the D-term on/off
    real(DP) :: dtau = 0.0_DP

    ! BETA-parameter that switch the B-term on/off.
    real(DP) :: dbeta = 0.0_DP
    
    ! Minimum allowed damping parameter; OMGMIN
    real(DP) :: domegaMin = 0.0_DP
    
    ! Maximum allowed damping parameter; OMGMAX
    real(DP) :: domegaMax = 2.0_DP
    
    ! Output level of the nonlinear iteration
    integer :: MT_OutputLevel = 2
    
    ! Minimum discretisation level
    integer :: NLMIN = 1
    
    ! Maximum discretisation level
    integer :: NLMAX = 1
    
    ! Pointer to the solution vector that is changed during the iteration
    type(t_vectorBlock), pointer :: p_rsolution => null()
    
    ! Pointer to the RHS vector
    type(t_vectorBlock), pointer :: p_rrhs => null()
    
    ! A filter chain that is used for implementing boundary conditions into
    ! (linear and nonlinear) defect vectors.
    type(t_filterChain), dimension(:), pointer :: p_RfilterChain
    
    ! A t_CHPreconditioner saving information about the preconditioner.
    type(t_CHPreconditioner) :: rpreconditioner
    
    ! A t_CHPreconditionerSpecials structure that saves information about
    ! special 'tweaks' in matrices such that everything works.
    type(t_CHPreconditionerSpecials) :: rprecSpecials
    
    ! An array of t_CHcoreEquationOneLevel structures for all levels
    ! of the discretisation.
    type(t_CHcoreEquationOneLevel), dimension(:), pointer :: RcoreEquation => null()
    
    ! Auxiliary variable: Saves the initial defect in the nonlinear iteration
    real(DP), dimension(2) :: DresidualInit = 0.0_DP
    
    ! Auxiliary variable: Saves the last defect in the nonlinear iteration
    real(DP), dimension(2) :: DresidualOld = 0.0_DP

    ! Auxiliary variable: Norm of the relative change = norm of the 
    ! preconditioned residual in the nonlinear iteration
    real(DP), dimension(2) :: DresidualCorr = 0.0_DP
    
    ! Auxiliary variable: Convergence criteria of the nonlinear solver
    real(DP), dimension(5) :: DepsNL = 0.0_DP
    
    ! Auxiliary variable: Last calculated damping parameter
    real(DP) :: domegaNL = 0.0_DP
    
    ! Auxiliary variable: Convergence rate of linear solver (if this is
    ! applied as preconditioner).
    real(DP) :: drhoLinearSolver = 0.0_DP
    
  end type

!</typeblock>

!</types>

contains

  ! ***************************************************************************
  ! callback routines for the nonlinear solver
  ! ***************************************************************************

!<subroutine>

  subroutine CH_solveCoreEquation (rCHproblem,rnonlinearIteration,rnlSolver,&
      rvector,rrhs,rtempBlock, rNSproblem, rNSvector, rNSrhs)
  
!<description>
  ! This routine invokes the nonlinear solver to solve the core equation
  ! as configured in the core equation structure.
!</description>

!<input>
  ! The right-hand-side vector to use in the equation.
  type(t_vectorBlock), intent(IN) :: rrhs
!</input>
  
!<inputoutput>
  ! The problem structure characterising the whole problem.
  type(t_CHproblem), intent(INOUT)                :: rCHproblem

  ! A nonlinear-iteration structure that configures the core equation to solve.
  ! Can be initialised e.g. by CH_createNonlinearLoop + manual setup of
  ! the coefficients of the terms.
  type(t_CHNonlinearIteration), intent(INOUT) :: rnonlinearIteration

  ! A nonlinear solver configuration.
  ! Can be initialised e.g. by using CH_getNonlinearSolver.
  type(t_nlsolNode), intent(INOUT) :: rnlSolver
  
  ! Initial solution vector. Is replaced by the new solution vector.
  type(t_vectorBlock), intent(INOUT) :: rvector

  ! A temporary block vector for the nonlinear iteration.
  ! OPTIONAL: If not specified, a temporary vector is automatically allocated.
  type(t_vectorBlock), intent(INOUT), target, optional :: rtempBlock
!~~~~~~~~~~~NS problem~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !optional
  type(t_problem), intent(INOUT), optional :: rNSproblem
  type(t_vectorBlock), intent(INOUT), target, optional :: rNSvector
  type(t_vectorBlock), intent(INOUT), target, optional :: rNSrhs
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_vectorBlock), pointer :: p_rtempBlock
    type(t_timer) :: rtimer

    if (.not. present(rtempBlock)) then
      ! Create a temporary vector we need for the nonlinear iteration.
      allocate (p_rtempBlock)
      call lsysbl_createVecBlockIndirect (rrhs, p_rtempBlock, .false.)
    else 
      p_rtempBlock => rtempBlock
    end if

    call stat_clearTimer(rtimer)
    call stat_startTimer(rtimer)

    ! call the nonlinear solver. For preconditioning and defect calculation,
    ! the solver calls our callback routines.
    call CH_nonlinearSolver(rCHproblem,rnonlinearIteration,rnlSolver,&
        rvector,rrhs,p_rtempBlock, rNSproblem, rNSvector, rNSrhs)

    ! Gather statistics
    call stat_stopTimer(rtimer)
!     rCHproblem%rstatistics%dtimeNonlinearSolver = &
!      rCHproblem%rstatistics%dtimeNonlinearSolver + rtimer%delapsedreal
      
!     rCHproblem%rstatistics%nnonlinearIterations = &
!      rCHproblem%rstatistics%nnonlinearIterations + rnlSolver%iiterations

    if (.not. present(rtempBlock)) then
      ! Release the temporary vector
      call lsysbl_releaseVector (p_rtempBlock)
      deallocate (p_rtempBlock)
    end if

  end subroutine

! ***************************************************************************
!<subroutine>

  subroutine CH_nonlinearSolver (rCHproblem,rnonlinearIteration,&
      rsolverNode,rx,rb,rd, rNSproblem, rNSvector, rNSrhs)
             
!<description>
  ! This routine invokes the nonlinear defect correction iteration
  !     $$  x_{n+1}  =  x_n  +  J^{-1} ( b - A(x_n) x_n )  $$
  !
  ! It's a modification of the routine nlsol_performSolve in the kernel
  ! to allow passing the problem structure to the different callback routines.
  !
  ! The defect correction loop is split into three tasks:
  !
  ! 1.) Calculation of nonlinear defect
  !         $$  d_n  :=  b - A(x_n)x_n  $$
  !     For this purpose, the routine CH_getDefect is called.
  !
  ! 2.) Preconditioning of the nonlinear defect
  !         $$  u_n  :=  J^{-1} d_n  $$
  !     with a preconditioner $J^{-1}$. The preconditioner is
  !     realised by the routine CH_precondDefect.
  !
  ! 3.) Correction
  !         $$  x_{n+1}  :=  x_n  +  \omega u_n  $$
  !     $\omega$ is usually = 1. The defined routine CH_precondDefect
  !     can modify $\omega$ in every iteration.
  !
  ! The callback routine CH_resNormCheck is called in every iteration.
  ! This routine can check the residuum for convergence/divergence and can print
  ! information about the norm of the residuum to screen.
  !
  ! At the beginning of the nonlinear iteration, the routines 
  ! CH_getResiduum and CH_resNormCheck are called once with ite=0 to calculate
  ! the initial defect, its norm and to check if already the initial vector
  ! is the solution.
  !
  ! If a linear solver is chosen for preconditioning, the nonlinear solver
  ! assumes that this is completely initialised by the application and ready 
  ! for action; it will directly call linsol_precondDefect without any further
  ! initialisation!
!</description>

!<inputoutput>
  ! The problem structure characterising the whole problem.
  type(t_CHproblem), intent(INOUT)                :: rCHproblem

  ! Reference to the nonlinear iteration structure that configures the
  ! main nonlinear equation. Intermediate data is changed during the iteration.
  type(t_CHnonlinearIteration), intent(INOUT)   :: rnonlinearIteration

  ! The nonlinear solver node that configures the solution process.
  type(t_nlsolNode), intent(INOUT)              :: rsolverNode
  
  ! INPUT: Initial solution vector.
  ! OUTPUT: Final iteration vector.
  type(t_vectorBlock), intent(INOUT)            :: rx
             
  ! Temporary vector. Must be of the same size/type as rx/rb.
  type(t_vectorBlock), intent(INOUT)            :: rd
!</inputoutput>
  
!<input>
  ! Right hand side vector of the equation.
  type(t_vectorBlock), intent(IN)               :: rb
!~~~~~~~~~~~NS problem~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !optional
  type(t_problem), intent(INOUT), optional :: rNSproblem
  type(t_vectorBlock), intent(INOUT), target, optional :: rNSvector
  type(t_vectorBlock), intent(IN), target, optional :: rNSrhs
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!</input>

!</subroutine>

    ! local variables
    integer :: ite
    !real(DP), dimension(NLSOL_MAXEQUATIONSERROR) :: DvecNorm
    type(t_vectorBlock) :: rtemp
    real(DP) :: domega
    integer :: nblocks
    logical :: bconvergence,bdivergence,bsuccess
  
    ! In case our preconditioner is a matrix-vector multiplication,
    ! allocate memory for another temporary vector used 
    ! during the MV multiplication.
    if (rsolverNode%cpreconditioner .eq. NLSOL_PREC_MATRIX) then
      call lsysbl_createVecBlockIndirect (rx,rtemp,.false.)
    end if
    
    ! Status reset
    rsolverNode%iresult = 0

 !     print *, rd%RvectorBlock(1)%h_Ddata
 !   print *, rb%RvectorBlock(1)%h_Ddata
 !   print *, 'before precon'

! MCai, debug 
!    print *, 'before get Defect'
    ! Calculate the initial nonlinear defect to rd:   d = b-A(x)x
    call CH_getDefect(rCHproblem,rnonlinearIteration,0,rx,rb,rd,rNSproblem, rNSvector)
!    print *, rd%RvectorBlock(1)%h_Ddata
!    print *, rb%RvectorBlock(1)%h_Ddata
!    print *, 'after get defec'
	
    ite = 0
    rsolverNode%icurrentIteration = ite

    ! The standard convergence/divergence test supports only up to 
    ! NLSOL_MAXEQUATIONSERROR equations.
    nblocks = min(rb%nblocks,NLSOL_MAXEQUATIONSERROR)

    ! Initial test for convergence/divergence.
    call CH_resNormCheck (rCHproblem,rnonlinearIteration,&
        ite,rx,rb,rd,bconvergence,bdivergence)

    ! Get the initial residuum; CH_resNormCheck saved that to DresidualOld.
    rsolverNode%DinitialDefect(1) = sqrt(rnonlinearIteration%DresidualInit(1)**2 + &
                                         rnonlinearIteration%DresidualInit(2)**2)
    rsolverNode%DfinalDefect(1) = rsolverNode%DinitialDefect(1)

!MCai, to make sure precondDefect  is called, we set bconvergence=.false.
!    print *, 'we make sure that precondDefect is called'
!     bconvergence =.false.

    ! Perform at least nminIterations iterations
    if (ite .lt. rsolverNode%nminIterations) bconvergence = .false.
  
    ! Check for divergence
    if (bdivergence) rsolverNode%iresult = 1
      
    if ((.not. bconvergence) .and. (.not. bdivergence)) then
    
      ! Let's do the nonlinear loop...
      !
      ! Initialise the domega-value for the damping of the correction
      ! as prescribed by the parameters of the solver.
      ! The callback routines might change it...
      domega = rsolverNode%domega
      
      ! Perform at most nmaxIterations iterations
      do ite = 1,rsolverNode%nmaxIterations
        rsolverNode%icurrentIteration = ite
      
        ! Perform preconditioning on the defect vector:  u = J^{-1} d
        bsuccess = .true.

        ! call CH_precondDefect to do the preconditioning.
        ! The routine is allowed to change domega during the
        ! iteration if necessary. The nonlinear solver here does not touch
        ! domega anymore, so the callback routine is the only one changing it.
        call CH_precondDefect (rCHproblem,rnonlinearIteration,&
            ite,rd,rx,rb,domega,bsuccess,rNSproblem, rNSvector)

! MCai, debug
!         print *, 'end of precond'
        ! If bsuccess=false, the preconditioner had an error.
        if (.not. bsuccess) then
          call output_line ('NLSOL: Iteration '//&
              trim(sys_siL(ite,10))//' canceled as the preconditioner went down!')
          rsolverNode%iresult = 3
          exit
        end if
        
        ! If domega=0.0, the solution vector would stay unchanged. In this
        ! case, the nonlinear solver would not proceed at all, and the next
        ! iteration would behave exactly as before!
        ! So in this case, there's nothing to do, we can stop the iteration.
        if (domega .eq. 0.0_DP) then
          call output_line ('NLSOL: Iteration '//&
              trim(sys_siL(ite,10))//' canceled as there is no progress anymore!')
          exit
        else
          ! Add the correction vector in rd to rx;
          ! damped by the damping parameter:           x := x + domega u
          call lsysbl_vectorLinearComb (rd,rx,domega,1.0_DP)

          ! Calculate the new nonlinear defect to rd:  d = b-A(x)x
          call CH_getDefect (rCHproblem,rnonlinearIteration,ite,rx,rb,rd,&
		                      rNSproblem, rNSvector)
          ! Check the defect for convergence.
          call CH_resNormCheck (rCHproblem,rnonlinearIteration,&
              ite,rx,rb,rd,bconvergence,bdivergence)
          ! Get the new residual; CH_resNormCheck saved that to DresidualOld.
          rsolverNode%DfinalDefect(1) = sqrt(rnonlinearIteration%DresidualOld(1)**2 + &
                                             rnonlinearIteration%DresidualOld(2)**2)

          ! Perform at least nminIterations iterations
          if (ite .lt. rsolverNode%nminIterations) bconvergence = .false.

          ! Check for convergence
          if (bconvergence) exit

          ! Check for divergence
          if (bdivergence) then
            rsolverNode%iresult = 1
            exit
          end if
        
        end if
        
      end do ! ite

!      print *,'we check convergence and divergence'
!      print *, bconvergence
!      print *, bdivergence

    end if ! not converged and not diverged

    ! Set ITE to NIT to prevent printing of "NIT+1" of the loop was
    ! completed

    if (ite .gt. rsolverNode%nmaxIterations) &
      ite = rsolverNode%nmaxIterations
      
    if (.not. bconvergence) then 
      ! Convergence criterion not reached, but solution did not diverge.
      rsolverNode%iresult = -1
    end if

    rsolverNode%iiterations = ite
    
    ! Release temporary memory
    if (rsolverNode%cpreconditioner .eq. NLSOL_PREC_MATRIX) then
      call lsysbl_releaseVector (rtemp)
    end if

    ! Nonlinear loop finished.

  end subroutine

  ! ***************************************************************************

  !<subroutine>
  
  subroutine CH_getDefect (rCHproblem,rnonlinearIteration,ite,rx,rb,rd, &
                            rNSproblem, rNSvector)
  
  use linearsystemblock
  use collection    
    
  !<description>
    ! FOR NONLINEAR ITERATION:
    ! Defect vector calculation callback routine. Based on the current iteration 
    ! vector rx and the right hand side vector rb, this routine has to compute the 
    ! defect vector rd.
  !</description>

  !<input>
    ! Number of current iteration. 0=build initial defect
  integer, intent(IN)                           :: ite

  ! Current iteration vector
  type(t_vectorBlock), intent(IN),target        :: rx

  ! Right hand side vector of the equation.
  type(t_vectorBlock), intent(IN)               :: rb
  !</input>
               
  !<inputoutput>
  ! The problem structure characterising the whole problem.
  type(t_CHproblem), intent(INOUT)                :: rCHproblem

  ! Reference to the nonlinear iteration structure that configures the
  ! main nonlinear equation. Intermediate data is changed during the iteration.
  type(t_CHnonlinearIteration), intent(INOUT)   :: rnonlinearIteration

  ! Defect vector b-A(x)x. This must be filled with data by the callback routine.
  type(t_vectorBlock), intent(INOUT)            :: rd


!~~~~~~~~~~~NS problem~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! OPTIONAL: 
  type(t_problem), intent(INOUT), optional :: rNSproblem
  type(t_vectorBlock), intent(INOUT), target, optional :: rNSvector
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !</inputoutput>
  
  !</subroutine>
    
	! local variables  
    ! The nonlinear iteration structure
    type(t_nonlinearCHMatrix) :: rnonlinearCHMatrix
    integer :: ilvmax
    type(t_filterChain), dimension(:), pointer :: p_RfilterChain
  
      ! a timer structure
    type(t_timer) :: rtimer
                  
    ! DEBUG!!!
    !real(DP), DIMENSION(:), POINTER :: p_Ddata,p_Ddata2

    call stat_clearTimer(rtimer)
    call stat_startTimer(rtimer)

    ! Build the nonlinear defect
    call lsysbl_copyVector (rb,rd)

    ! DEBUG!!!
    !call lsysbl_getbase_double (rx,p_Ddata)
    !call lsysbl_getbase_double (rd,p_Ddata2)

    ilvmax = rnonlinearIteration%NLMAX

    ! Initialise the matrix assembly structure rnonlinearCHMatrix to describe the
    ! matrix we want to have.
    rnonlinearCHMatrix%dalpha = rnonlinearIteration%dalpha
    rnonlinearCHMatrix%dgamma = rnonlinearIteration%dgamma
    rnonlinearCHMatrix%deta = rnonlinearIteration%deta
    rnonlinearCHMatrix%dtau = rnonlinearIteration%dtau
   rnonlinearCHMatrix%dtheta = rnonlinearIteration%dtheta
    rnonlinearCHMatrix%dbeta = rnonlinearIteration%dbeta

    rnonlinearCHMatrix%p_rdiscretisation => &
        rnonlinearIteration%RcoreEquation(ilvmax)%p_rdiscretisation
!    rnonlinearCHMatrix%p_rmatrixA => &
!         rnonlinearIteration%RcoreEquation(ilvmax)%p_rmatrixA
!     rnonlinearCHMatrix%p_rmatrixB => &
!         rnonlinearIteration%RcoreEquation(ilvmax)%p_rmatrixB
!     rnonlinearCHMatrix%p_rmatrixC => &
!         rnonlinearIteration%RcoreEquation(ilvmax)%p_rmatrixC
!     rnonlinearCHMatrix%p_rmatrixD => &
!         rnonlinearIteration%RcoreEquation(ilvmax)%p_rmatrixD

    rnonlinearCHMatrix%p_rmatrixMass => &
        rnonlinearIteration%RcoreEquation(ilvmax)%p_rmatrixMass
    rnonlinearCHMatrix%p_rmatrixLaplace => &
        rnonlinearIteration%RcoreEquation(ilvmax)%p_rmatrixLaplace
!    rnonlinearCHMatrix%p_rmatrixConv => &
!        rnonlinearIteration%RcoreEquation(ilvmax)%p_rmatrixConv

    ! returns  rd := dcx A(ry) rx + dcd rd
    !call CH_nonlinearMatMul (rnonlinearCHMatrix,rx,rd,-1.0_DP,1.0_DP)

    call CH_nonlinearMatMul (rnonlinearCHMatrix,rx,rd,-1.0_DP,1.0_DP,rx, &
                                 rNSproblem, rNSvector)

      p_RfilterChain => rnonlinearIteration%p_RfilterChain
      if (associated(p_RfilterChain)) then    
        call filter_applyFilterChainVec (rd, p_RfilterChain)
      end if
      
      ! Gather statistics
    call stat_stopTimer(rtimer)
!      rCHproblem%rstatistics%dtimeDefectCalculation = &
!        rCHproblem%rstatistics%dtimeDefectCalculation + rtimer%delapsedreal
      
  end subroutine
  ! ***************************************************************************

  !<subroutine>

  subroutine CH_precondDefect (rCHproblem,rnonlinearIteration,&
                          ite,rd,rx,rb,domega,bsuccess, &
                           rNSproblem, rNSvector)
  
  use linearsystemblock
  use collection
    
  !<description>
    ! FOR NONLINEAR ITERATION:
    ! Defect vector calculation callback routine. Based on the current iteration 
    ! vector rx and the right hand side vector rb, this routine has to compute the 
    ! defect vector rd. 
  !</description>

  !<inputoutput>
    ! The problem structure characterising the whole problem.
    type(t_CHproblem), intent(INOUT)                :: rCHproblem

    ! Reference to the nonlinear iteration structure that configures the
    ! main nonlinear equation. Intermediate data is changed during the iteration.
    type(t_CHnonlinearIteration), intent(INOUT), target   :: rnonlinearIteration

    ! Number of current iteration. 
    integer, intent(IN)                           :: ite

    ! Defect vector b-A(x)x. This must be replaced by J^{-1} rd by a preconditioner.
    type(t_vectorBlock), intent(INOUT)            :: rd

    ! Damping parameter. Is set to rsolverNode%domega (usually = 1.0_DP)
    ! on the first call to the callback routine.
    ! The callback routine can modify this parameter according to any suitable
    ! algorithm to calculate an 'optimal damping' parameter. The nonlinear loop
    ! will then use this for adding rd to the solution vector:
    ! $$ x_{n+1} = x_n + domega*rd $$
    ! domega will stay at this value until it's changed again.
    real(DP), intent(INOUT)                       :: domega

    ! If the preconditioning was a success. Is normally automatically set to
    ! TRUE. If there is an error in the preconditioner, this flag can be
    ! set to FALSE. In this case, the nonlinear solver breaks down with
    ! the error flag set to 'preconditioner broke down'.
    logical, intent(INOUT)                        :: bsuccess
  !</inputoutput>
  
  !<input>
    ! Current iteration vector
    type(t_vectorBlock), intent(IN), target       :: rx

    ! Current right hand side of the nonlinear system
    type(t_vectorBlock), intent(IN), target       :: rb
!~~~~~~~~~~~NS problem~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !optional
  type(t_problem), intent(INOUT), optional :: rNSproblem
  type(t_vectorBlock), intent(INOUT), target, optional :: rNSvector
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !</input>
  
  !</subroutine>
    
	! local variables 
    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    type(t_vectorScalar), pointer :: p_rvectorTemp,p_rvectorTemp2
    type(t_vectorBlock) :: rtemp1,rtemp2
    integer :: ierror
    type(t_linsolNode), pointer :: p_rsolverNode,p_rcgrSolver
    integer, dimension(2) :: Cnorms
    
    integer :: i
    real(DP) :: dresInit,dres,dtempdef
    logical :: bassembleNewton
    type(t_matrixBlock), dimension(:), pointer :: Rmatrices
    type(t_CHDynamicNewtonControl), pointer :: p_rnewton
    type(t_filterChain), dimension(:), pointer :: p_RfilterChain

    ! a local timer object
    type(t_timer) :: rtimer

    ! DEBUG!!!
    !real(dp), dimension(:), pointer :: p_vec,p_def,p_da
    !call lsysbl_getbase_double (rd,p_def)
    !call lsysbl_getbase_double (rx,p_vec)

!~~~~~~~~MCai, debug
!      print *, 'we make sure that precondiint is used'
!      print *, rnonlinearIteration%rpreconditioner%ctypePreconditioning
!      print *, CHPREC_LINEARSOLVER

!      rnonlinearIteration%rpreconditioner%ctypePreconditioning=CHPREC_LINEARSOLVER

      select case (rnonlinearIteration%rpreconditioner%ctypePreconditioning)
      case (CHPREC_NONE)
        ! No preconditioning
        domega = max(rnonlinearIteration%domegamin, &
                    min(rnonlinearIteration%domegamax,domega))

      case (CHPREC_LINEARSOLVER,CHPREC_NEWTON,CHPREC_NEWTONDYNAMIC)
        ! Preconditioning with a linear solver.
        !
        ! At first, assemble the preconditioner matrices on all levels
        ! and incorporate all boundary conditions.
        !
        ! Should we assemble the Newton part?

        bassembleNewton = .false.

        if (rnonlinearIteration%rpreconditioner%ctypePreconditioning .eq. &
            CHPREC_NEWTON) then

          ! Use Newton in any case.
          bassembleNewton = .true.

        else if (rnonlinearIteration%rpreconditioner%ctypePreconditioning .eq. &
            CHPREC_NEWTONDYNAMIC) then

          ! Adaptive Newton. Check the iteration and the residual whether to use
          ! Newton or not. But at first, get a shortcut to the parameter structure
          ! of the adaptive Newton...

          p_rnewton => rnonlinearIteration%rpreconditioner%radaptiveNewton

          if (ite .gt. p_rnewton%nmaxFixPointIterations) then
            ! Force Newton to be used.
            bassembleNewton = .true.
          else
            if (ite .gt. p_rnewton%nminFixPointIterations) then
              ! In this case, the residuum of the last iterate decides on 
              ! whether to use Newton or not.
              dresInit = sqrt(rnonlinearIteration%DresidualInit(1)**2 + &
                            rnonlinearIteration%DresidualInit(2)**2)
              dres    = sqrt(rnonlinearIteration%DresidualOld(1)**2 + &
                            rnonlinearIteration%DresidualOld(2)**2)
              if ((dres .lt. p_rnewton%depsAbsNewton) .and. &
                  (dres .lt. p_rnewton%depsRelNewton * dresInit)) then
                bassembleNewton = .true.
              end if
            end if

            ! Otherwise: Use fixpoint iteration...
          end if
        end if
        call stat_clearTimer(rtimer)
        call stat_startTimer(rtimer)

        ! Assemble the preconditioner matrices in rnonlinearIteration
        ! on all levels that the solver uses.
!MCai, debug
!        print *, 'assembleLinsoMatrix is called'
        call assembleLinsolMatrices (rCHproblem, rnonlinearIteration,rCHproblem%rcollection,&
            bassembleNewton,rx,rnonlinearIteration%NLMIN,rnonlinearIteration%NLMAX, &
             rNSproblem, rNSvector)

        ! Gather statistics
        call stat_stopTimer(rtimer)
 !       rCHproblem%rstatistics%dtimeMatrixAssembly = &
 !         rCHproblem%rstatistics%dtimeMatrixAssembly + rtimer%delapsedreal

        call stat_clearTimer(rtimer)
        call stat_startTimer(rtimer)
        
        ! Our 'parent' (the caller of the nonlinear solver) has prepared
        ! a preconditioner node for us (a linear solver with symbolically
        ! factorised matrices). Get this from the collection.
      
!       print *, 'attach solver node'
!       print *, rnonlinearIteration%rpreconditioner%ctypePreconditioning
        p_rsolverNode => rnonlinearIteration%rpreconditioner%p_rsolverNode
! MCai, the following code is for Dynamic Newton. ..
!         if (rnonlinearIteration%rpreconditioner%ctypePreconditioning .eq. &
!             CHPREC_NEWTONDYNAMIC) then
!           ! Do we have to apply the inexact Newton?
!           if (p_rnewton%cinexactNewton .ne. 0) then
!           
!             ! Adaptive stopping criterion / inexact Newton active.
!             !
!             ! Determine the stopping criterion for the linear solver.
!             ! This is an adaptive stopping criterion depending on the current
!             ! defect. In detail, for the inexact Newton we have
!             !
!             !   |b-Ax_{i+1}|         ( |b-Ax_i| ) exp             ( |b-Ax_i| )
!             !   ------------ = min { ( -------- )     , depsrel * ( -------- ) }
!             !     |b-Ax_0|           ( |b-Ax_0| )                 ( |b-Ax_0| )
!             !
!             ! see e.g. [Michael Hinze, Habilitation, p. 51]
!             !
!             ! If Newton is not active, we taje the formula
!             !
!             !   |b-Ax_{i+1}|             ( |b-Ax_i| )
!             !   ------------ = depsrel * ( -------- ) 
!             !     |b-Ax_0|               ( |b-Ax_0| )
!             !
!             ! to always gain depsrel.
!             ! Switch off the relative stopping criterion in the linear solver:
!             
!             p_rsolverNode%istoppingCriterion = 0
!             
!             ! Just for safetyness, gain at least one digit.
!             p_rsolverNode%depsRel = 1.0E-1_DP
!             
!             ! Calculate the new absolute stopping criterion:
!             dresInit = sqrt(rnonlinearIteration%DresidualInit(1)**2 + &
!                           rnonlinearIteration%DresidualInit(2)**2)
!             dres    = sqrt(rnonlinearIteration%DresidualOld(1)**2 + &
!                           rnonlinearIteration%DresidualOld(2)**2)
!             
!             dtempdef = dres / dresInit
!             
!             if (bassembleNewton) then
!               p_rsolverNode%depsAbs = &
!                   MIN(dtempDef**p_rnewton%dinexactNewtonExponent, &
!                       p_rnewton%dinexactNewtonEpsRel*dtempdef) * dresInit
!             else      
!               p_rsolverNode%depsAbs = p_rnewton%dinexactNewtonEpsRel*dtempdef*dresInit
!             end if
!             
!             ! If we have a multigrid solver, we also have to take care for
!             ! the coarse grid solver!
!             i = rnonlinearIteration%NLMIN
!             p_rcgrSolver => rnonlinearIteration%rpreconditioner%p_rcgrSolver
!             if (associated(p_rcgrSolver)) then
!               ! For the coarse grid solver, we choose the same stopping criterion.
!               ! But just for safetyness, the coarse grid solver should gain at least
!               ! one digit!
!               p_rcgrSolver%istoppingCriterion = 0
!               p_rcgrSolver%depsRel = 1.0E-1_DP
!               p_rcgrSolver%depsAbs = p_rsolverNode%depsAbs
!             end if
!             
!           end if
!           
!         end if

        ! Re-attach the system matrices to the solver.
        ! Note that no pointers and no handles are changed, so we can savely do
        ! that without calling linsol_doneStructure/linsol_doneStructure.
        ! This simply informs the solver about possible new scaling factors
        ! in the matrices in case they have changed...
!MCai
!        allocate(Rmatrices(rnonlinearIteration%NLMIN:rnonlinearIteration%NLMAX))
        allocate(Rmatrices(rnonlinearIteration%NLMAX:rnonlinearIteration%NLMAX))
       
        do i=rnonlinearIteration%NLMAX,rnonlinearIteration%NLMAX
!        do i=rnonlinearIteration%NLMIN,rnonlinearIteration%NLMAX
          call lsysbl_duplicateMatrix ( &
            rnonlinearIteration%RcoreEquation(i)%p_rmatrixPreconditioner, &
            Rmatrices(i), LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
        end do
        
        call linsol_setMatrices(&
            rnonlinearIteration%rpreconditioner%p_rsolverNode,Rmatrices(:))

        ! The solver got the matrices; clean up Rmatrices, it was only of temporary
        ! nature...
        do i=rnonlinearIteration%NLMAX,rnonlinearIteration%NLMAX
!        do i=rnonlinearIteration%NLMIN,rnonlinearIteration%NLMAX
          call lsysbl_releaseMatrix (Rmatrices(i))
        end do
        deallocate(Rmatrices)

        ! Initialise data of the solver. This in fact performs a numeric
        ! factorisation of the matrices in UMFPACK-like solvers.
        call linsol_initData (p_rsolverNode, ierror)
        if (ierror .ne. LINSOL_ERR_NOERROR) stop
        
        ! Gather statistics
        call stat_stopTimer(rtimer)
!        rCHproblem%rstatistics%dtimeLinearSolverFactorisation = &
!          rCHproblem%rstatistics%dtimeLinearSolverFactorisation + rtimer%delapsedreal
        
        call stat_clearTimer(rtimer)
        call stat_startTimer(rtimer)
        
        ! Finally solve the system. As we want to solve Ax=b with
        ! b being the real RHS and x being the real solution vector,
        ! we use linsol_solveAdaptively. If b is a defect
        ! RHS and x a defect update to be added to a solution vector,
        ! we would have to use linsol_precondDefect instead.
        call linsol_precondDefect (p_rsolverNode,rd)

        ! Gather statistics
        call stat_stopTimer(rtimer)
!        rCHproblem%rstatistics%dtimeLinearSolver = &
!          rCHproblem%rstatistics%dtimeLinearSolver + rtimer%delapsedreal
!        rCHproblem%rstatistics%nlinearIterations = &
!          rCHproblem%rstatistics%nlinearIterations + p_rsolverNode%iiterations
        
        ! Remember convergence rate for output
        rnonlinearIteration%drhoLinearSolver = p_rsolverNode%dconvergenceRate

        ! Release the numeric factorisation of the matrix.
        ! We don't release the symbolic factorisation, as we can use them
        ! for the next iteration.
        call linsol_doneData (p_rsolverNode)
        
        ! Did the preconditioner work?
        bsuccess = p_rsolverNode%iresult .eq. 0
        
      end select

! Debug
!       print *, 'what is preconditionging?'
!       print *, rnonlinearIteration%rpreconditioner%ctypePreconditioning

      ! Finally calculate a new damping parameter domega.
      if (rnonlinearIteration%rpreconditioner%ctypePreconditioning .ne. CHPREC_NONE) then
        
        ! For this purpose, we need two temporary vectors.
        ! On one hand, we have p_rvectorTemp.
        ! Get the second temporary vector from the collection as it was
        ! prepared by our 'parent' that invoked the nonlinear solver.
        p_rvectorTemp => rnonlinearIteration%rpreconditioner%p_rtempVectorSc
        p_rvectorTemp2 => rnonlinearIteration%rpreconditioner%p_rtempVectorSc2
        
        ! Both temp vectors are scalar, but we need block-vectors in the
        ! structure of rx/rb. Derive block vectors in that structure that
        ! share their memory with the scalar temp vectors. Note that the
        ! temp vectors are created large enough by our parent!
        call lsysbl_createVecFromScalar (p_rvectorTemp,rtemp1)
        call lsysbl_enforceStructure (rb,rtemp1)

        call lsysbl_createVecFromScalar (p_rvectorTemp2,rtemp2)
        call lsysbl_enforceStructure (rb,rtemp2)

        ! Calculate the omega
        call CH_getOptimalDamping (rCHproblem,rnonlinearIteration,&
            rd,rx,rb,rtemp1,rtemp2,domega)

        ! Remember damping parameter for output
        rnonlinearIteration%domegaNL = domega

        ! Release the temp block vectors. This only cleans up the structure.
        ! The data is not released from heap as it belongs to the
        ! scalar temp vectors.
        call lsysbl_releaseVector (rtemp2)
        call lsysbl_releaseVector (rtemp1)
      end if
      
      ! Calculate the max-norm of the correction vector.
      ! This is used for the stopping criterium in CH_resNormCheck!
      Cnorms(:) = LINALG_NORMMAX
      call lsysbl_vectorNormBlock(rd,Cnorms,rnonlinearIteration%DresidualCorr)
      
      if ((.not. bsuccess) .and. (domega .ge. 0.001_DP)) then
        ! The preconditioner did actually not work, but the solution is not
        ! 'too bad'. So we accept it still.
        bsuccess = .true.
      end if
      
      if (bsuccess) then
        ! Filter the final defect
        p_RfilterChain => rnonlinearIteration%p_RfilterChain
        call filter_applyFilterChainVec (rd, p_RfilterChain)
        call vecfil_discreteNLSlipBCdef (rd)
      end if

    contains
    
      subroutine assembleLinsolMatrices (rCHproblem, rnonlinearIteration,rcollection,&
          bassembleNewton,rx,NLMIN,NLMAX, rNSproblem, rNSvector)

      use linearsystemblock
      use collection

      ! Assembles on every level a matrix for the linear-solver/Newton preconditioner.
      ! bnewton allows to specify whether the Newton matrix or only the standard
      ! system matrix is evaluated. The output is written to the p_rpreconditioner 
      ! matrices specified in the rnonlinearIteration structure.
 
! MCai,
      type(t_CHproblem), intent(INOUT) :: rCHproblem

      ! Reference to a collection structure that contains all parameters of the
      ! discretisation (for nonlinearity, etc.).
      type(t_collection), intent(INOUT)                :: rcollection

      ! Nonlinear iteration structure where to write the linearised system matrix to.
      type(t_CHnonlinearIteration), intent(INOUT)      :: rnonlinearIteration

      ! TRUE  = Assemble the Newton preconditioner.
      ! FALSE = Assemble the standard defect correction preconditioner
      !         (i.e. the linearised system matrix).
      logical, intent(IN) :: bassembleNewton
      
      ! Minimum level of the preconditioner that is to be initialised.
      integer, intent(IN)                              :: NLMIN
      
      ! Maximum level of the preconditioner that is to be initialised.
      ! This must corresponds to the last matrix in Rmatrices.
      integer, intent(IN)                              :: NLMAX
      
      ! Current iteration vector. 
      type(t_vectorBlock), intent(IN), target          :: rx

!~~~~~~~~~~~~~~~parameter from NS problem~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      type(t_problem), intent(INOUT), optional :: rNSproblem
      type(t_vectorBlock), intent(INOUT), optional :: rNSvector
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! local variables
      real(DP) :: dnewton
      integer :: ilev,icol
      integer(PREC_MATIDX), dimension(1) :: Irows = (/1/)
      type(t_matrixBlock), pointer :: p_rmatrix,p_rmatrixFine
      type(t_vectorScalar), pointer :: p_rvectorTemp
      type(t_vectorBlock), pointer :: p_rvectorFine,p_rvectorCoarse
      type(t_nonlinearCHMatrix) :: rnonlinearCHMatrix
      type(t_matrixScalar), pointer :: p_rmatrixTemplateFEM
      type(t_interlevelProjectionBlock), pointer :: p_rprojection
      type(t_timer) :: rtimer
      ! A filter chain for the linear solver
      type(t_filterChain), dimension(:), pointer :: p_RfilterChain

      ! DEBUG!!!
    !    real(dp), dimension(:), pointer :: p_vec,p_def,p_da
    !    call lsysbl_getbase_double (rd,p_def)
    !    call lsysbl_getbase_double (rx,p_vec)

        ! Get the interlevel projection structure and the temporary vector
        ! from the collection.
        ! Our 'parent' prepared there how to interpolate the solution on the
        ! fine grid to coarser grids.
        p_rprojection => rnonlinearIteration%rpreconditioner%p_rprojection
        p_rvectorTemp => rnonlinearIteration%rpreconditioner%p_rtempVectorSc

        ! Get the filter chain. We need tghat later to filter the matrices.        
        p_RfilterChain => rnonlinearIteration%p_RfilterChain

        ! On all levels, we have to set up the nonlinear system matrix,
        ! so that the linear solver can be applied to it.
        
        nullify(p_rmatrix)
        ilev=NLMAX

!        do ilev=NLMAX,NLMIN,-1
        
          ! Get the matrix on the current level.
          ! Shift the previous matrix to the pointer of the fine grid matrix.
          p_rmatrixFine => p_rmatrix

!MCai, before use the following code We need to make sure that RcoreEquation(ilev)%p_rmatrixPreconditioner
! is initialized or assembled. 
          allocate(rnonlinearIteration%RcoreEquation(ilev)%p_rmatrixPreconditioner)
          call lsysbl_createMatBlockByDiscr (rnonlinearIteration%RcoreEquation(ilev)%p_rdiscretisation,&
                            rnonlinearIteration%RcoreEquation(ilev)%p_rmatrixPreconditioner)

      p_rmatrixTemplateFEM => rCHproblem%RlevelInfo(NLMAX)%rmatrixTemplateFEM

      call lsyssc_duplicateMatrix (p_rmatrixTemplateFEM,&
        rnonlinearIteration%RcoreEquation(ilev)%p_rmatrixPreconditioner%RmatrixBlock(1,1), &
        LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
      call lsyssc_duplicateMatrix (p_rmatrixTemplateFEM,&
        rnonlinearIteration%RcoreEquation(ilev)%p_rmatrixPreconditioner%RmatrixBlock(1,2), &
        LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
      call lsyssc_duplicateMatrix (p_rmatrixTemplateFEM,&
        rnonlinearIteration%RcoreEquation(ilev)%p_rmatrixPreconditioner%RmatrixBlock(2,1), &
        LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
      call lsyssc_duplicateMatrix (p_rmatrixTemplateFEM,&
        rnonlinearIteration%RcoreEquation(ilev)%p_rmatrixPreconditioner%RmatrixBlock(2,2), &
        LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

          p_rmatrix => rnonlinearIteration%RcoreEquation(ilev)%p_rmatrixPreconditioner
! Debug
!          print *, rnonlinearIteration%RcoreEquation(ilev)%p_rmatrixPreconditioner%rmatrixBlock(1,1)%h_Da
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          ! On the highest level, we use rx as solution to build the nonlinear
          ! matrix. On lower levels, we have to create a solution
          ! on that level from a fine-grid solution before we can use
          ! it to build the matrix!
          if (ilev .eq. NLMAX) then
            p_rvectorCoarse => rx
          end if

          ! Initialise the matrix assembly structure rnonlinearCHMatrix to describe the
          ! matrix we want to have.
          rnonlinearCHMatrix%dalpha = rnonlinearIteration%dalpha
          rnonlinearCHMatrix%dgamma = rnonlinearIteration%dgamma
          rnonlinearCHMatrix%deta = rnonlinearIteration%deta
          rnonlinearCHMatrix%dtau = rnonlinearIteration%dtau
		  rnonlinearCHMatrix%dtheta = rnonlinearIteration%dtheta
          rnonlinearCHMatrix%dbeta = rnonlinearIteration%dbeta

          rnonlinearCHMatrix%iadaptiveMatrices = &
              rnonlinearIteration%rprecSpecials%iadaptiveMatrices
          rnonlinearCHMatrix%dadmatthreshold = &
              rnonlinearIteration%rprecSpecials%dadmatthreshold

          rnonlinearCHMatrix%p_rdiscretisation => &
              rnonlinearIteration%RcoreEquation(ilev)%p_rdiscretisation
          rnonlinearCHMatrix%p_rmatrixMass => &
              rnonlinearIteration%RcoreEquation(ilev)%p_rmatrixMass
          rnonlinearCHMatrix%p_rmatrixLaplace => &
              rnonlinearIteration%RcoreEquation(ilev)%p_rmatrixLaplace
!          rnonlinearCHMatrix%p_rmatrixConv => &
!              rnonlinearIteration%RcoreEquation(ilev)%p_rmatrixConv

!           rnonlinearCHMatrix%p_rmatrixA => &
!               rnonlinearIteration%RcoreEquation(ilev)%p_rmatrixA
!           rnonlinearCHMatrix%p_rmatrixB => &
!               rnonlinearIteration%RcoreEquation(ilev)%p_rmatrixB
!           rnonlinearCHMatrix%p_rmatrixC => &
!               rnonlinearIteration%RcoreEquation(ilev)%p_rmatrixC
!           rnonlinearCHMatrix%p_rmatrixD => &
!               rnonlinearIteration%RcoreEquation(ilev)%p_rmatrixD
      
 
!          print *, 'we should have no problem in assemb, we need preconditioning'
!          print *, ilev
  
          ! Assemble the matrix.
          ! If we are on a lower level, we can specify a 'fine-grid' matrix.
          if (ilev .eq. NLMAX) then
             ! we remark the following, it comes from cc2d
            !call CH_assembleMatrix (CHMASM_COMPUTE,CHMASM_MTP_AUTOMATIC,&
            !    p_rmatrix,rnonlinearCHMatrix,p_rvectorCoarse)
   !           stop

            call parlst_readfromfile(rCHproblem%rparamList, 'data/CH_solver.dat')
          
            call CH_assembleMatrix (rCHproblem, p_rmatrix,&
                         rnonlinearCHMatrix,rx, rCHproblem%rparamList,&
                                      rNSproblem, rNSvector)
          end if

          ! Boundary conditions
          ! ---------------------------------------------------

          if (associated(p_RfilterChain)) then
            ! Apply the filter chain to the matrix.
            ! As the filter consists only of an implementation filter for
            ! boundary conditions, this implements the boundary conditions
            ! into the system matrix.
            call filter_applyFilterChainMat (p_rmatrix, p_RfilterChain)
          else
            ! call the matrix filter for the boundary conditions to include the BC's
            ! into the matrix.
            call matfil_discreteBC (p_rmatrix)
            call matfil_discreteFBC (p_rmatrix)
          end if
            
!        end do
        
      end subroutine
      
    end subroutine
  ! ***************************************************************************

!<subroutine>

    subroutine CH_initNonlinearLoop (rCHproblem,nlmin,nlmax,rvector,rrhs,&
        rnonlinearIteration,sname)
  
  !<description>
    ! Initialises the given nonlinear iteration structure rnonlinearIteration.
    ! Creates the structure with CH_createNonlinearLoop and saves all
    ! problem dependent parameters and matrices in it.
    ! The routine automatically calls CH_createNonlinearLoop to initialise
    ! the structure with internal parameters.
    ! Note: This does not initialise the preconditioner in that structure!
  !</description>
  
  !<input>
    ! A problem structure saving problem-dependent information.
    type(t_CHproblem), intent(INOUT), target :: rCHproblem

    ! Minimum refinement level in the rCHproblem structure that is allowed to be used
    ! by the preconditioners.
    integer, intent(IN) :: nlmin
  
    ! Maximum refinement level in the rCHproblem structure that is allowed to be used
    ! by the preconditioners. This level must correspond to rvector and rrhs.
    integer, intent(IN) :: nlmax

    ! The solution vector which is modified later during the nonlinear iteration.
    type(t_vectorBlock), intent(IN), target :: rvector

    ! The right-hand-side vector to use in the equation
    type(t_vectorBlock), intent(IN), target :: rrhs

    ! Name of the section in the parameter list containing the parameters
    ! of the nonlinear solver.
    character(LEN=*), intent(IN) :: sname
  !</input>

  !<output>
    ! Nonlinar iteration structure saving data for the callback routines.
    ! Is filled with data.
    type(t_CHnonlinearIteration), intent(OUT) :: rnonlinearIteration
  !</output>

  !</subroutine>

      ! local variables
      integer :: ilevel
      logical :: bneumann
      type(t_parlstSection), pointer :: p_rsection

      rnonlinearIteration%NLMIN = NLMIN
      rnonlinearIteration%NLMAX = NLMAX

      ! Initialise the matrix pointers on all levels that we have to maintain.
      allocate(rnonlinearIteration%RcoreEquation(NLMIN:NLMAX))

!      rnonlinearIteration%MT_OutputLevel = rCHproblem%MT_OutputLevel
      rnonlinearIteration%NLMAX = rCHproblem%NLMAX
    
      ! Get the minimum/maximum damping parameter from the parameter list, save
      ! them to the nonlinear iteration structure (which is now initialised).
      call parlst_getvalue_double (rCHproblem%rparamList, 'CC2D-NONLINEAR', &
                                 'domegaMin', rnonlinearIteration%domegaMin, 0.0_DP)
                              
      call parlst_getvalue_double (rCHproblem%rparamList, 'CC2D-NONLINEAR', &
                                 'domegaMax', rnonlinearIteration%domegaMax, 2.0_DP)
    
      ! Save pointers to the RHS and solution vector
      rnonlinearIteration%p_rsolution => rvector
      rnonlinearIteration%p_rrhs => rrhs
    
      ! Set the preconditioner to 'nothing'
      rnonlinearIteration%rpreconditioner%ctypePreconditioning = -1
    
      ! Deactivate any 'tweak' flags in the final-assembly structure
      rnonlinearIteration%rprecSpecials%iadaptiveMatrices = 0
    
      ! Assign the matrix pointers in the nonlinear iteration structure to
      ! all our matrices that we want to use.
      do ilevel = nlmin,nlmax
        rnonlinearIteration%RcoreEquation(ilevel)%p_rdiscretisation => &
          rCHproblem%RlevelInfo(ilevel)%rdiscretisation
        rnonlinearIteration%RcoreEquation(ilevel)%p_rmatrixMass => &
          rCHproblem%RlevelInfo(ilevel)%rmatrixMass
        rnonlinearIteration%RcoreEquation(ilevel)%p_rmatrixLaplace => &
          rCHproblem%RlevelInfo(ilevel)%rmatrixLaplace
!        rnonlinearIteration%RcoreEquation(ilevel)%p_rmatrixConv => &
!          rCHproblem%RlevelInfo(ilevel)%rmatrixConv

!         rnonlinearIteration%RcoreEquation(ilevel)%p_rmatrixA => &
!           rCHproblem%RlevelInfo(ilevel)%rmatrixA
!         rnonlinearIteration%RcoreEquation(ilevel)%p_rmatrixB => &
!           rCHproblem%RlevelInfo(ilevel)%rmatrixB
!         rnonlinearIteration%RcoreEquation(ilevel)%p_rmatrixC => &
!           rCHproblem%RlevelInfo(ilevel)%rmatrixC
!         rnonlinearIteration%RcoreEquation(ilevel)%p_rmatrixD => &
!           rCHproblem%RlevelInfo(ilevel)%rmatrixD


!        rnonlinearIteration%RcoreEquation(ilevel)%p_rtempVector => &
!          rCHproblem%RlevelInfo(ilevel)%rtempVector

      end do
      
      ! Clear auxiliary variables for the nonlinear iteration
      rnonlinearIteration%DresidualInit = 0.0_DP
      rnonlinearIteration%DresidualOld  = 0.0_DP
    
      call parlst_querysection(rCHproblem%rparamList, sname, p_rsection) 

      if (.not. associated(p_rsection)) then
        call output_line ('Cannot create nonlinear solver; no section '''//&
          trim(sname)//'''!', &
          OU_CLASS_ERROR,OU_MODE_STD,'CH_initNonlinearLoop')
        call sys_halt()
      end if

      ! Get stopping criteria of the nonlinear iteration
      call parlst_getvalue_double (p_rsection, 'depsD', &
                                 rnonlinearIteration%DepsNL(1), 0.1_DP)

      call parlst_getvalue_double (p_rsection, 'depsUR', &
                                 rnonlinearIteration%DepsNL(3), 0.1_DP)

      call parlst_getvalue_double (p_rsection, 'depsPR', &
                                 rnonlinearIteration%DepsNL(4), 0.1_DP)

      call parlst_getvalue_double (p_rsection, 'ddmpD', &
                                 rnonlinearIteration%DepsNL(5), 0.1_DP)
      
      ! Set up a filter that modifies the block vectors/matrix
      ! according to boundary conditions. This filter chain is applied to each
      ! defect vector during the linear and nonlinear iteration.
      allocate(rnonlinearIteration%p_RfilterChain(3))
    
      ! Initialise the first filter of the filter chain as boundary
      ! implementation filter for defect vectors:
      rnonlinearIteration%p_RfilterChain(1)%ifilterType = FILTER_DISCBCDEFreal

      ! Do we have Neumann boundary?
      !
      ! The bhasNeumannBoundary flag of the higher level decides about that...
!      bneumann = rCHproblem%RlevelInfo(rCHproblem%NLMAX)%bhasNeumannBoundary
!      rnonlinearIteration%p_RfilterChain(3)%ifilterType = FILTER_DONOTHING

    end subroutine
  ! ***************************************************************************

!<subroutine>

    subroutine CH_doneNonlinearLoop (rnonlinearIteration)
  
!<description>
  ! Releases memory allocated in cc_initNonlinearLoop.
  ! The routine automatically calls cc_releaseNonlinearLoop to release
  ! internal parameters.
!</description>

!<inputoutput>
  ! The nonlinear iteration structure that should be cleaned up.
    type(t_CHNonlinearIteration), intent(INOUT) :: rnonlinearIteration
!</inputoutput>

!</subroutine>
    
    ! Release the filter chain for the defect vectors.
      if (associated(rnonlinearIteration%p_RfilterChain)) &
        deallocate(rnonlinearIteration%p_RfilterChain)
      
      if (associated(rnonlinearIteration%RcoreEquation)) &
      deallocate(rnonlinearIteration%RcoreEquation)

      rnonlinearIteration%NLMIN = 0
      rnonlinearIteration%NLMAX = 0

    end subroutine
   

  ! ***************************************************************************

  !<subroutine>

    subroutine CH_getOptimalDamping (rCHproblem,rnonlinearIteration,rd,rx,rb,&
        rtemp1,rtemp2,domega, rNSproblem, rNSvector)
  
  !<description>
    ! This subroutine is called inside of the nonlinear loop, to be precise,
    ! inside of CH_precondDefect. It calculates an optiman damping parameter
    ! for the nonlinear defect correction.
    !
    ! The nonlinear loop reads:
    !
    !     $$ u_(n+1) = u_n + OMEGA * C^{-1}d_n $$
    !
    ! with $d_n$ the nonlinear defect and $C^{-1}$ a preconditioner (usually
    ! the linearised system).
    ! Based on the current solution $u_n$, the defect vector $d_n$, the RHS 
    ! vector $f_n$ and the previous parameter OMEGA, a new 
    ! OMEGA=domega value is calculated.
    !
    ! The nonlinear system matrix on the finest level in the collection is 
    ! overwritten by $A(u_n+domega_{old}*C^{-1}d_n)$.
  !</description>

  !<input>
    ! Current iteration vector
    type(t_vectorBlock), intent(IN)               :: rx

    ! Current RHS vector of the nonlinear equation
    type(t_vectorBlock), intent(IN)               :: rb

    ! Defect vector b-A(x)x. 
    type(t_vectorBlock), intent(IN)               :: rd
  !</input>

  !<inputoutput>
    ! The problem structure characterising the whole problem.
    type(t_CHproblem), intent(INOUT)                :: rCHproblem

    ! Reference to the nonlinear iteration structure that configures the
    ! main nonlinear equation. Intermediate data is changed during the iteration.
    type(t_CHnonlinearIteration), intent(INOUT)   :: rnonlinearIteration

    ! A temporary vector in the structure of rx
    type(t_vectorBlock), intent(INOUT)            :: rtemp1

    ! A 2nd temporary vector in the structure of rx
    type(t_vectorBlock), intent(INOUT)            :: rtemp2

    ! Damping parameter. On entry: an initial value given e.g. by the
    ! previous step.
    ! On return: The new damping parameter.
    real(DP), intent(INOUT)                       :: domega
!~~~~~~~~~~~~~~~~~~~NS problem~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  type(t_problem), intent(INOUT), optional :: rNSproblem
  type(t_vectorBlock), intent(INOUT), optional :: rNSvector
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !</inputoutput>
  
  !</subroutine>

    ! local variables
    integer :: ilvmax
    real(DP) :: dskv1,dskv2
    type(t_matrixBlock) :: rmatrix
    type(t_timer) :: rtimer
    ! A filter chain for the linear solver
    type(t_filterChain), dimension(:), pointer :: p_RfilterChain
    
    type(t_nonlinearCHMatrix) :: rnonlinearCHMatrix

!    ! DEBUG!!!:
!    real(dp), dimension(:), pointer :: p_vec,p_def,p_temp1,p_temp2,p_da
!    call lsysbl_getbase_double (rd,p_def)
!    call lsysbl_getbase_double (rx,p_vec)
!    call lsysbl_getbase_double (rtemp1,p_temp1)
!    call lsysbl_getbase_double (rtemp2,p_temp2)
!    ilvmax = collct_getvalue_int (p_rcollection,'NLMAX')

      ! Is there anything to do?
      if (rnonlinearIteration%domegaMin .ge. rnonlinearIteration%domegaMax) then
        ! No - cancel.
        domega = rnonlinearIteration%domegaMin
        return
      end if
      
      call stat_clearTimer(rtimer)
      call stat_startTimer(rtimer)
      
      ! Get minimum/maximum level from the collection
      ilvmax = rnonlinearIteration%NLMAX
      
      p_RfilterChain => rnonlinearIteration%p_RfilterChain
      
      ! We now want to calculate a new OMEGA parameter
      ! with OMGMIN < OMEGA < OMGMAX.
      !
      ! The defect correction for a problem like T(u)u=f has the form
      !
      !       u_(n+1)  =  u_n  +  OMEGA * C * ( f - T(u_n)u_n )
      !                =  u_n  +  OMEGA * d_n
      !
      ! with an appropriate preconditioner C, which we don't care here.
      ! In our case, this iteration system can be written as:
      !
      ! (u1)     (u1)                     ( (f1)   [ A   B] (u1) )
      ! (u2)  := (u2)  + OMEGA * p^{-1} * ( (f2) - [ D   C] (u2) )
    !
      !                                   |------------------------------|
      !                                              = d_n
      !                            |-------------------------------------|
      !                                        = Y = (y1,y2,yp) = rd
      !
      ! with KST1=KST1(u1,u2,p) and Y=rd being the solution from
      ! the Oseen equation with 
      !
      !                  [ A  B ]
      !    C = T(u_n) = 
      !                  [ D  C ]
      !
      ! The parameter OMEGA is calculated as the result of the 1D
      ! minimisation problem:
      !
      !   OMEGA = min_omega || T(u^l+omega*Y)*(u^l+omega*Y) - f ||_E
      !
      !           < T(u^l+omegaold*Y)Y , f - T(u^l+omegaold*Y)u^l >
      !        ~= -------------------------------------------------
      !              < T(u^l+omegaold*Y)Y , T(u^l+omegaold*Y)Y >
      !
      ! when choosing omegaold=previous omega, which is a good choice
      ! as one can see by linearisation (see p. 170, Turek's book).
      !
      ! Here, ||.||_E denotes the the Euclidian norm to the Euclidian 
      ! scalar product <.,.>.
      
      ! ==================================================================
      ! First term of scalar product in the nominator
      !
      ! Calculate the new nonlinear block A at the
      ! point rtemp1 = u_n + omegaold*Y
      ! ==================================================================
      !
      ! At first, calculate the point rtemp1 = u_n+omegaold*Y where
      ! to evaluate the matrix.

      call lsysbl_copyVector(rd,rtemp1)
      call lsysbl_vectorLinearComb (rx,rtemp1,1.0_DP,domega)

      ! Should we discretise the Cahn Hilliard nonlinearity?
      ! Would mean to rebuild the system matrix. Otherwise we can reuse
      ! our previous diffusion matrix without rebuilding it.
      
      ! Re-assemble the nonlinear system matrix on the maximum level.
      ! Initialise the matrix assembly structure rnonlinearCHMatrix to describe the
      ! matrix we want to have.

      rnonlinearCHMatrix%dalpha = rnonlinearIteration%dalpha
      rnonlinearCHMatrix%dgamma = rnonlinearIteration%dgamma
      rnonlinearCHMatrix%deta = rnonlinearIteration%deta
      rnonlinearCHMatrix%dtau = rnonlinearIteration%dtau
      rnonlinearCHMatrix%dtheta = rnonlinearIteration%dtheta
      rnonlinearCHMatrix%dbeta = rnonlinearIteration%dbeta     
     
	  rnonlinearCHMatrix%iadaptiveMatrices = &
      rnonlinearIteration%rprecSpecials%iadaptiveMatrices
      rnonlinearCHMatrix%dadmatthreshold = &
          rnonlinearIteration%rprecSpecials%dadmatthreshold
      rnonlinearCHMatrix%p_rdiscretisation => &
          rnonlinearIteration%RcoreEquation(ilvmax)%p_rdiscretisation
      rnonlinearCHMatrix%p_rmatrixMass => &
          rnonlinearIteration%RcoreEquation(ilvmax)%p_rmatrixMass
      rnonlinearCHMatrix%p_rmatrixLaplace => &
          rnonlinearIteration%RcoreEquation(ilvmax)%p_rmatrixLaplace
!      rnonlinearCHMatrix%p_rmatrixConv => &
!          rnonlinearIteration%RcoreEquation(ilvmax)%p_rmatrixConv

!      rnonlinearCHMatrix%p_rmatrixA => &
!           rnonlinearIteration%RcoreEquation(ilvmax)%p_rmatrixA
!       rnonlinearCHMatrix%p_rmatrixB => &
!           rnonlinearIteration%RcoreEquation(ilvmax)%p_rmatrixB
!       rnonlinearCHMatrix%p_rmatrixC=> &
!           rnonlinearIteration%RcoreEquation(ilvmax)%p_rmatrixC
!       rnonlinearCHMatrix%p_rmatrixD => &
!           rnonlinearIteration%RcoreEquation(ilvmax)%p_rmatrixD

      ! Assemble the matrix.        
      !call CH_assembleMatrix (CHMASM_COMPUTE,CHMASM_MTP_AUTOMATIC,&
      !    rmatrix,rnonlinearCHMatrix,rtemp1)


         call CH_assembleMatrix (rCHproblem, rmatrix,&
                          rnonlinearCHMatrix,rx,rCHproblem%rparamList,&
                                    rNSproblem, rNSvector)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      ! We don't have to implement any boundary conditions into the matrix
      ! as we apply an appropriate filter to the defect vector after
      ! each matrix-vector-multiplication below!
      
      ! ==================================================================
      ! Second term of the scalar product in the nominator
      ! Calculate the defect rtemp2 = F-T*u_n.
      ! ==================================================================

      call lsysbl_copyVector (rb,rtemp2)
      call lsysbl_blockMatVec (rmatrix, rx, rtemp2, -1.0_DP, 1.0_DP)
      
      ! This is a defect vector - filter it! This e.g. implements boundary
      ! conditions.
      if (associated(p_RfilterChain)) then
        call filter_applyFilterChainVec (rtemp2, p_RfilterChain)
      end if
      
      ! Filter the resulting defect vector through the slip-boundary-
      ! condition vector filter for implementing nonlinear slip boundary
      ! conditions into a defect vector.
      call vecfil_discreteNLSlipBCdef (rtemp2)
      
      ! ==================================================================
      ! For all terms in the fraction:
      ! Calculate the value  rtemp1 = T*Y
      ! ==================================================================

      call lsysbl_blockMatVec (rmatrix, rd, rtemp1, 1.0_DP, 0.0_DP)
      
      ! This is a defect vector against 0 - filter it! This e.g. 
      ! implements boundary conditions.
      if (associated(p_RfilterChain)) then
        call filter_applyFilterChainVec (rtemp1, p_RfilterChain)
      end if
      
      ! Filter the resulting defect vector through the slip-boundary-
      ! condition vector filter for implementing nonlinear slip boundary
      ! conditions into a defect vector.
      call vecfil_discreteNLSlipBCdef (rtemp1)

      ! Release the matrix again
      call lsysbl_releaseMatrix (rmatrix)
      
      ! ==================================================================
      ! Calculation of the fraction terms.
      ! Calculate nominator:    dskv1:= (T*Y,D)   = (rtemp1,rtemp2)
      ! Calculate denominator:  dskv2:= (T*Y,T*Y) = (rtemp1,rtemp1)
      ! ==================================================================
      
      dskv1 = lsysbl_scalarProduct (rtemp1, rtemp2)
      dskv2 = lsysbl_scalarProduct (rtemp1, rtemp1)
      
      if (dskv2 .lt. 1.0E-40_DP) then
        call output_line ('dskv2 nearly zero. Optimal damping parameter singular.', &
            OU_CLASS_ERROR,OU_MODE_STD,'CH_getOptimalDamping')
        call output_line ('Is the triangulation ok??? .tri-file destroyed?', &
            OU_CLASS_ERROR,OU_MODE_STD,'CH_getOptimalDamping')
        call output_line ('Boundary conditions set up properly?', &
            OU_CLASS_ERROR,OU_MODE_STD,'CH_getOptimalDamping')
        call sys_halt()
        stop
      end if
      
      ! Ok, we have the nominator and the denominator. Divide them
      ! by each other to calculate the new OMEGA.
      
      domega = dskv1 / dskv2
      
      ! And make sure it's in the allowed range:
      
      domega = max(rnonlinearIteration%domegamin, &
                   min(rnonlinearIteration%domegamax,domega))
      
      ! That's it, we have our new Omega.

      ! Gather statistics
      call stat_stopTimer(rtimer)
!      rCHproblem%rstatistics%dtimeOptimalCorrection = &
!        rCHproblem%rstatistics%dtimeDefectCalculation + rtimer%delapsedreal
      
    end subroutine


  ! ***************************************************************************

    subroutine CH_resNormCheck (rCHproblem,rnonlinearIteration,&
        ite,rx,rb,rd,bconvergence,bdivergence)
  
    use linearsystemblock
    use collection

  !<description>
    ! Residual norm calculation & printing routine.
    ! This routine is called each time the norm of the residuum was calculated.
    ! It has to check the current residuum for convergence and/or divergence
    ! and can print the residuum to screen.
  !</description>

  !<inputoutput>
    ! The problem structure characterising the whole problem.
    type(t_CHproblem), intent(INOUT)                :: rCHproblem

    ! Reference to the nonlinear iteration structure that configures the
    ! main nonlinear equation. Intermediate data is changed during the iteration.
    type(t_CHnonlinearIteration), intent(INOUT)   :: rnonlinearIteration

    ! Number of current iteration. Is set to 0 when the callback routine
    ! is called the first time. In this situation, rd describes the initial
    ! defect vector.
    integer, intent(IN)                           :: ite
  
    ! Current iteration vector
    type(t_vectorBlock), intent(IN), target       :: rx

    ! Right hand side vector of the equation.
    type(t_vectorBlock), intent(IN), target       :: rb

    ! Defect vector b-A(x)x.
    type(t_vectorBlock), intent(IN), target       :: rd
  !</inputoutput>
  
  !<output>
    ! Must be set to TRUE by the callback routine if the residuum rd
    ! is within a desired tolerance, so that the solver should treat
    ! the iteration as 'converged'.
    logical, intent(OUT)                        :: bconvergence
  
    ! Must be set to TRUE by the callback routine if the residuum rd
    ! is out of a desired tolerance, so that the solver should treat
    ! the iteration as 'diverged'.
    logical, intent(OUT)                        :: bdivergence
  !</output>

      ! local variables
      real(DP), dimension(3) :: Dresiduals
      real(DP) :: dresOld,drhoNL,ddelP,ddelU,dtmp,dresU,dresP,dres, dresINIT
      real(DP) :: depsD,depsUR,depsPR,depsRES
      integer, dimension(2) :: Cnorms

      ! Calculate norms of the solution/defect vector
      call CH_getDefectNorm (rx,rb,rd,Dresiduals)
      Dresiduals(3) = sqrt(Dresiduals(1)**2 + Dresiduals(2)**2)

      ! In the first iteration (initial defect), print the norm of the defect
      ! and save the norm of the initial residuum to the structure
      if (ite .eq. 0) then
      
        call output_separator (OU_SEP_MINUS)     
        call output_line (' IT  RELPhi   RELChP   DEFPhi   DEFChP ' // &
                          ' DEF-TOT   RHONL   OMEGNL   RHOMG')
        call output_separator (OU_SEP_MINUS)     
        call output_line ('  0                  '// &
            trim(sys_sdEP(Dresiduals(1),9,2))//&
            trim(sys_sdEP(Dresiduals(2),9,2))//&
	    trim(sys_sdEP(Dresiduals(3),9,2)))

        call output_separator (OU_SEP_MINUS)

        rnonlinearIteration%DresidualInit (1:2) = Dresiduals(1:2)
        rnonlinearIteration%DresidualOld (1:2) = Dresiduals(1:2)

        bconvergence = .false.
        bdivergence = .false.
      
      else
        ! In the other iterations, calculate the relative change and
        ! test the convergence criteria.
      
        ! Old defect:
        dresOld = sqrt(rnonlinearIteration%DresidualOld(1)**2 + &
                       rnonlinearIteration%DresidualOld(2)**2)
        
        ! Replace the 'old' residual by the current one
        rnonlinearIteration%DresidualOld(1:2) = Dresiduals(1:2)
        
        ! Calculate norms of the solution/defect vector, calculated above
	dresU = Dresiduals(1)
        dresP = Dresiduals(2)
        dres  = sqrt(dresU**2 + dresP**2)
        
        ! Calculate relative maximum changes 
        ! This simply calculates some postprocessing values of the relative
        ! change in the solution.
        !
        ! Maximum norm of solution vector:
        !
        !              || YP ||_max       || Pnew - Pold ||_max
        !   DELP := ------------------- = ---------------------
        !              || P ||_max           || Pnew ||_max
        !
        !
        ! Relative change of solution vector:
        !
        !            || (Y1,Y2) ||_max    || Unew - Uold ||_max
        !   DELU := ------------------- = --------------------- 
        !           || (KU1,KU2) ||_max       || Unew ||_max
        !
        ! The norms || YP ||_max, || Yi ||_max are saved in the nonlinear
        ! iteration structure from the last preconditioning!
        ! The other MAX-norms have to be calculated from U

! MCai, 

        Cnorms(:) = LINALG_NORMMAX
        call lsysbl_vectorNormBlock (rx,Cnorms,Dresiduals)

        dtmp = max(Dresiduals(1),Dresiduals(2))
        if (dtmp .lt. 1.0E-8_DP) dtmp = 1.0_DP
		ddelU = rnonlinearIteration%DresidualCorr(1)/dtmp
		ddelP = rnonlinearIteration%DresidualCorr(2)/dtmp
        
        ! Check if the nonlinear iteration can prematurely terminate.
        !        
        ! Get the stopping criteria from the parameters.
        ! Use the DepsNL data according to the initialisation above.
        dresINIT = sqrt(rnonlinearIteration%DresidualInit(1)**2 + &
                        rnonlinearIteration%DresidualInit(2)**2)

        ! dresInit=0 may hardly occur -- except when we expect 'no flow'.
        ! But to prevent a check against "something<=0" in this case below,
        ! set dresInit to something <> 0.
        if (dresINIT .eq. 0.0_DP) dresINIT = 1.0_DP

        depsD   = rnonlinearIteration%DepsNL(1)
        depsUR  = rnonlinearIteration%DepsNL(3)
        depsPR  = rnonlinearIteration%DepsNL(4)
        depsRES = rnonlinearIteration%DepsNL(5)*dresINIT

       
        ! All residual information calculated.
        ! Check for divergence; use a 'NOT' for better NaN handling.
        bdivergence = .not. (dres/dresOld .lt. 1E3)
! MCai, debug
!        print *, 'we need to rewrite this part to check stopping criterion'

!        print *, dres/dresOld 
!        print *, 'problem is here'
!        print *, bdivergence
! 
!         print *, (ddelU .le. depsUR)
!         print *, (ddelP .le. depsPR)
!         print *, (dresU .le. depsD)
!         print *, (dresU .le. depsD)
!         print *, (dres .le. depsRES)
!         print *, dresINIT
!         print *, depsRES
        ! Check for convergence
        if((ddelU .le. depsUR).and.(ddelP .le. depsPR) .and.  &
		   (dresU .le. depsD) .and.(dresP .le. depsD).and. &
                   (dres .le. depsRES)) then
		  bconvergence = .true.
        else
          bconvergence = .false.
        end if

        if ((ddelU .lt. SYS_EPSreal*1E2_DP) .and. &
            (ddelP .lt. SYS_EPSreal*1E2_DP)) then
          ! We are hard on machine exactness, so stop the iteraton
          bconvergence =.true.
        end if
        
        ! Print residual information
        call output_line ( &
            trim(sys_si(ite,3))//' '//&
            trim(sys_sdEP(ddelU,9,2))// &
            trim(sys_sdEP(ddelP,9,2))// &
            trim(sys_sdEP(dresU,9,2))// &
            trim(sys_sdEP(dresP,9,2))// &
            trim(sys_sdEP(dres,9,2))// &
            trim(sys_sdEP(drhoNL,9,2))// &
            trim(sys_sdEP(rnonlinearIteration%domegaNL,9,2))// &
            trim(sys_sdEP(rnonlinearIteration%drhoLinearSolver,9,2)) &
            )
        
      end if
      
    end subroutine

! ***************************************************************************

!<subroutine>

  subroutine CH_getDefectNorm (rvector,rrhs,rdefect,Dresiduals)
  
!<description>
  ! Calculates a couple of norms from a given solution, defect and RHS vector
  ! of the nonlinear system. This can be used to check convergence criteria
  ! etc.
!</description>

!<input>
  ! The solution vector which is modified later during the nonlinear iteration.
  type(t_vectorBlock), intent(IN) :: rvector

  ! The right-hand-side vector to use in the equation
  type(t_vectorBlock), intent(IN) :: rrhs
  
  ! A defect vector calculated with rvector and rrhs
  type(t_vectorBlock), intent(IN) :: rdefect
!</input>

!<output>
  ! An array receiving different defect norms calculated by the above vectors.
  ! Dresiduals(1) = RESPhi   = ||defect_phi|| / ||rhs|| = phi residual
  ! Dresiduals(2) = RESChP = ||defect_ChP|| / ||rhs||   = chemical potential residual
  real(DP), dimension(:), intent(OUT) :: Dresiduals
!</output>

!</subroutine>

    ! local variables
    real(DP) :: dresF,DresTmp(2)
    integer, dimension(2) :: Cnorms

    !-----------------------------------------------------------------------
    !     Compute the relative l2-norms  RESU
    !-----------------------------------------------------------------------

    Cnorms(:) = LINALG_NORMEUCLID

    ! RESF := sqrt ( ||F1||_E**2 + ||F2||_E**2 )

    call lsysbl_vectorNormBlock (rrhs,Cnorms,DresTmp)
    dresF = max(DresTmp(1), DresTmp(2))
    
    if (dresF .lt. 1.0E-8_DP) dresF = 1.0_DP

    !               || (D1,D2) ||_E
    ! RESU = -----------------------------
    !        sqrt( ||F1||_E**2+||F2||_E**2 )

    call lsysbl_vectorNormBlock (rdefect,Cnorms,DresTmp)

    Dresiduals(1) = DresTmp(1)/dresF
    Dresiduals(2) = DresTmp(2)/dresF
    ! DNORMU = || (U1,U2) ||_l2 
!     print *, 'check why divergence'
!     print *, DresTmp(1)
!     print *, DresTmp(2)

  end subroutine

! ***************************************************************************
! Routines to solve the core equation
! ***************************************************************************

!<subroutine>

  subroutine CH_getNonlinearSolver (rnlSolver, rparamList, sname)
  
!<description>
  ! Creates a nonlinear solver node rnlSolver and initialises it with parameters
  ! from the INI/DAT files given in the parameter list rparamList.
  ! sname is the name of a section in rparamList that configures the
  ! parameter of the nonlinear solver.
!</description>

!<input>
  ! Parameter list that contains the parameters from the INI/DAT file(s).
  type(t_parlist), intent(IN) :: rparamList
  
  ! Name of the section in the parameter list containing the parameters
  ! of the nonlinear solver.
  character(LEN=*), intent(IN) :: sname
!</input>

!<output>
  ! A t_nlsolNode structure that contains the configuration of the nonlinear
  ! solver. The parameters are initialised according to the information
  ! in the section sname of the parameter list rparamList
  type(t_nlsolNode) :: rnlSolver
!</output>

!</subroutine>

    ! local variables
    type(t_parlstSection), pointer :: p_rsection

    ! Check that there is a section called sname - otherwise we
    ! cannot create anything!
    
    call parlst_querysection(rparamList, sname, p_rsection) 

    if (.not. associated(p_rsection)) then
      call output_line ('Cannot create nonlinear solver; no section '''//&
          trim(sname)//'''!', &
          OU_CLASS_ERROR,OU_MODE_STD,'CH_getNonlinearSolver')
      call sys_halt()
    end if
    
    ! Parse the given parameters now to initialise the solver node.
    ! This is now CHxD-specific!
    
    call parlst_getvalue_int (p_rsection, 'nminIterations', &
                              rnlSolver%nminIterations, rnlSolver%nminIterations)

    call parlst_getvalue_int (p_rsection, 'nmaxIterations', &
                              rnlSolver%nmaxIterations, rnlSolver%nmaxIterations)

    call parlst_getvalue_int (p_rsection, 'ioutputLevel', &
                              rnlSolver%ioutputLevel, rnlSolver%ioutputLevel)

    call parlst_getvalue_double (p_rsection, 'depsUR', &
                                 rnlSolver%DepsRel(1), rnlSolver%DepsRel(1))
    rnlSolver%DepsRel(2) = rnlSolver%DepsRel(1)

    call parlst_getvalue_double (p_rsection, 'depsD', &
                                 rnlSolver%DepsAbs(1), rnlSolver%DepsAbs(1))
    rnlSolver%DepsAbs(2) = rnlSolver%DepsAbs(1)

    ! Initial damping parameter.
    call parlst_getvalue_double (p_rsection, 'domegaIni', &
                                 rnlSolver%domega, rnlSolver%domega)

  end subroutine

 
  ! ***************************************************************************

!<subroutine>

  subroutine CH_doneLinearSolver (rnonlinearIteration)
  
!<description>
  ! Releases information of the linear solver preconditioner from the
  ! structure rnonlinearIteration. Cleans up what was configured
  ! in cc_initLinearSolver.
!</description>

!<inputoutput>
  ! The nonlinear iteration structure to be cleaned up.
  type(t_CHnonlinearIteration), intent(INOUT) :: rnonlinearIteration
!</inputoutput>

!</subroutine>

    ! Release the solver and all subsolvers.
    call linsol_releaseSolver(rnonlinearIteration%rpreconditioner%p_rsolverNode)

  end subroutine

 ! ***************************************************************************

!<subroutine>

  subroutine CH_initLinearSolver (rproblem,rnonlinearIteration,ssection)
  
!<description>
  ! This routine initialises a linear solver structure to be used
  ! for preconditioning of the nonlinear defect.
!</description>

!<input>
  ! Name of the section in the parameter list that contains the parameters
  ! of the linear solver.
  character(LEN=*), intent(IN) :: ssection
!</input>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_CHproblem), intent(INOUT), target :: rproblem

  ! The nonlinear iteration structure to which a preconditioner should
  ! be initialised. The parameters for the linear solver are written to
  ! the rpreconditioner substructure.
  ! The rprecSpecials structure is configured according to the linear
  ! solver.
  type(t_CHnonlinearIteration), intent(INOUT) :: rnonlinearIteration
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_parlstSection), pointer :: p_rsection
    type(t_parlist), pointer :: p_rparamList
    integer :: nlevels, ilev, nsm
    
    integer :: isolverType,ismootherType,icoarseGridSolverType
    character(LEN=SYS_STRLEN) :: sstring,ssolverSection,ssmootherSection
    character(LEN=SYS_STRLEN) :: scoarseGridSolverSection,spreconditionerSection
    type(t_linsolMG2LevelInfo), pointer :: p_rlevelInfo
    type(t_linsolNode), pointer :: p_rpreconditioner, p_rsmoother
    type(t_linsolNode), pointer :: p_rsolverNode

    ! Check that there is a section called ssolverName - otherwise we
    ! cannot create anything!

!MCai, the following part can be read from data file    
!    p_rparamList => rproblem%rparamList
        
!    call parlst_querysection(p_rparamList, ssection, p_rsection) 
    
!    if (.not. associated(p_rsection)) then
!      call output_line ('Cannot create linear solver; no section '''//trim(ssection)//&
!                        '''!', OU_CLASS_ERROR,OU_MODE_STD,'cc_initLinearSolver')
!      call sys_halt()
!    end if
    
    ! Get the parameters that configure the solver type
    
!    call parlst_getvalue_int (p_rsection, 'isolverType', isolverType, 1)
!    call parlst_getvalue_int (p_rsection, 'ismootherType', ismootherType, 3)
!    call parlst_getvalue_int (p_rsection, 'icoarseGridSolverType', &
!        icoarseGridSolverType, 1)
    isolverType=0
	ismootherType=1
	icoarseGridSolverType=1
	    
    rnonlinearIteration%rprecSpecials%isolverType = isolverType
    rnonlinearIteration%rprecSpecials%ismootherType = ismootherType
    rnonlinearIteration%rprecSpecials%icoarseGridSolverType = icoarseGridSolverType

!    call parlst_getvalue_string (p_rsection, 'ssolverSection', sstring,'')
!    read (sstring,*) ssolverSection
!    call parlst_getvalue_string (p_rsection, 'ssmootherSection', sstring,'')
!    read (sstring,*) ssmootherSection
!    call parlst_getvalue_string (p_rsection, 'scoarseGridSolverSection', sstring,'')
!    read (sstring,*) scoarseGridSolverSection
    
    ! Which type of solver do we have?
    
    select case (isolverType)
    
    case (0)
      ! This is the UMFPACK solver. Very easy to initialise. No parameters at all.
      call linsol_initUMFPACK4 (p_rsolverNode)

    case (1)
      print *, 'we only have Gauss elimination for preconditioning part, MG is not implemented'
    end select    

    ! Put the final solver node to the preconditioner structure.
    rnonlinearIteration%rpreconditioner%p_rsolverNode => p_rsolverNode

  end subroutine


  ! ***************************************************************************

!<subroutine>

  subroutine CH_initPreconditioner (rCHproblem,rnonlinearIteration,rCHvector,rrhs,&
		  rNSproblem, rNSvector)
  
!<description>
  ! This routine prepares the preconditioner that us used during the
  ! nonlinear iteration. The structure rpreconditioner will be initialised
  ! based on the information in rCHproblem.
  ! Necessary variables will be added to the collection structure in
  ! rCHproblem%rcollection to be available in the callback routines.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_CHproblem), intent(INOUT), target :: rCHproblem
!</inputoutput>

!<input>
  ! The current solution vector.
  type(t_vectorBlock), intent(IN) :: rCHvector

  ! The right-hand-side vector to use in the equation
  type(t_vectorBlock), intent(IN) :: rrhs
!</input>

!<inputoutput>
  ! Nonlinar iteration structure saving data for the callback routines.
  ! This is configured according to the preconditioner as specified in
  ! the DAT files.
  type(t_CHnonlinearIteration), intent(INOUT) :: rnonlinearIteration

!~~~~~~~~~~~~~~~~~~~NS problem~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  type(t_problem), intent(INOUT), optional :: rNSproblem
  type(t_vectorBlock), intent(INOUT), optional :: rNSvector
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: NLMIN,NLMAX
    integer :: i
    integer(PREC_VECIDX) :: imaxmem
    character(LEN=PARLST_MLDATA) :: ssolverName,sstring

!MCai
    ! First initialisation.
    ! Has to be set to TRUE on the first call. Initialises the preconditioner
    ! with the structure of the matrices.
    logical :: binit=.TRUE.

    ! Whether the structure of the system matrices is new.
    ! This variable has to be set to TRUE whenever there was a structure in 
    ! the system matrices. This reinitialises the linear solver.
    logical :: bstructuralUpdate=.TRUE.
  
    ! A level-info structure specifying the matrices of the problem.
    type(t_CHproblem_lvl), target :: rlevelInfo
!MCai, we do not need cmatrixType at current stage. 
!    integer :: cmatrixType=CHMASM_MTP_FULLTENSOR

    type(t_nonlinearCHmatrix) :: rnonlinearCHMatrix

    ! A pointer to the discretisation structure with the data.
    type(t_blockDiscretisation), POINTER :: p_rdiscretisation

    type(t_matrixScalar), pointer :: p_rmatrixTemplateFEM

    type(t_matrixBlock), dimension(:), pointer :: Rmatrices
    type(t_matrixBlock) :: rmatrix
    integer :: ierror
    type(t_matrixBlock), pointer ::  p_rmatrixPreconditioner

!</input>

    ! At first, ask the parameters in the INI/DAT file which type of 
    ! preconditioner is to be used. The data in the preconditioner structure
    ! is to be initialised appropriately!
!MCai
    call parlst_readfromfile(rCHproblem%rparamList, 'data/nonlinsol_cc2d.dat')
!~MCai
    call parlst_getvalue_int (rCHproblem%rparamList, 'CC2D-NONLINEAR', &
        'itypePreconditioning', &
        rnonlinearIteration%rpreconditioner%ctypePreconditioning, 1)
 
! MCai, debug....
!    print *, 'whether we have precond'
!    print *, rnonlinearIteration%rpreconditioner%ctypePreconditioning

    select case (rnonlinearIteration%rpreconditioner%ctypePreconditioning)
    case (CHPREC_NONE)
      ! No preconditioner
    case (CHPREC_LINEARSOLVER,CHPREC_NEWTON,CHPREC_NEWTONDYNAMIC)
      ! Ok, we have to initialise a linear solver for solving the linearised
      ! problem.
      !
      ! Which levels have we to take care of during the solution process?
      NLMIN = rnonlinearIteration%NLMIN
      NLMAX = rnonlinearIteration%NLMAX
      
      ! Figure out the name of the section that contains the information
      ! about the linear subsolver. Ask the parameter list from the INI/DAT file
      ! for the 'slinearSolver' value
      call parlst_getvalue_string (rCHproblem%rparamList, 'CC2D-NONLINEAR', &
                                  'slinearSolver', sstring, '')
      ssolverName = ''
      if (sstring .ne. '') read (sstring,*) ssolverName
      if (ssolverName .eq. '') then
        call output_line ('No linear subsolver!', &
            OU_CLASS_ERROR,OU_MODE_STD,'CH_initPreconditioner')
        call sys_halt()
      end if

      ! Initialise a standard interlevel projection structure. We
      ! can use the same structure for all levels. Therefore it's enough
      ! to initialise one structure using the RHS vector on the finest
      ! level to specify the shape of the PDE-discretisation.
      allocate(rnonlinearIteration%rpreconditioner%p_rprojection)
      call mlprj_initProjectionVec (&
          rnonlinearIteration%rpreconditioner%p_rprojection,rrhs)
      
      ! Initialise the projection structure with data from the INI/DAT
      ! files. This allows to configure prolongation/restriction.
!      call CH_getProlRest (rnonlinearIteration%rpreconditioner%p_rprojection, &
!          rCHproblem%rparamList,  'CC-PROLREST')
 
!~~~~~~~~~~~Think about, what kind of solver need to be used?~~~~~~~~~~~~~~~~~
      ! Initialise the linear solver as configured in the DAT file.
      call CH_initLinearSolver (rCHproblem,rnonlinearIteration,ssolverName)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~      
      ! How much memory is necessary for performing the level change?
      ! We ourself must build nonlinear matrices on multiple levels and have
      ! to interpolate the solution vector from finer level to coarser ones.
      ! We need temporary memory for this purpose...

      imaxmem = 0
!MCai, originally, it is 
!      do i=NLMIN+1,NLMAX
      do i=NLMAX,NLMAX
 
        ! Pass the system metrices on the coarse/fine grid to 
        ! mlprj_getTempMemoryMat to specify the discretisation structures
        ! of all equations in the PDE there.
        imaxmem = max(imaxmem,mlprj_getTempMemoryDirect (&
            rnonlinearIteration%rpreconditioner%p_rprojection,&
            rCHproblem%RlevelInfo(i-1)% &
              rdiscretisation%RspatialDiscr(1:rrhs%nblocks),&
            rCHproblem%RlevelInfo(i)% &
              rdiscretisation%RspatialDiscr(1:rrhs%nblocks)))
      end do

      ! Set up a scalar temporary vector that we need for building up nonlinear
      ! matrices. It must be at least as large as MAXMEM and NEQ(finest level),
      ! as we use it for resorting vectors, too.
      allocate(rnonlinearIteration%rpreconditioner%p_rtempVectorSc)
      call lsyssc_createVector (rnonlinearIteration%rpreconditioner%p_rtempVectorSc,&
                                max(imaxmem,rrhs%NEQ),.false.)
      
      ! Set up a second temporary vector that we need for calculating
      ! the optimal defect correction.
      allocate(rnonlinearIteration%rpreconditioner%p_rtempVectorSc2)
      call lsyssc_createVector (rnonlinearIteration%rpreconditioner%p_rtempVectorSc2,&
                                rrhs%NEQ,.false.,ST_DOUBLE)

	 ! Initialise the matrices.
      call CH_updatePreconditioner (rCHproblem,rnonlinearIteration,&
          rCHvector,rrhs,.true.,.true., &
		  rNSproblem, rNSvector)
     ! or 
     ! call CH_updatePreconditioner (rCHproblem,rnonlinearIteration,&
     !     rCHvector,rrhs,binit,bstructuralUpdate)    
    case DEFAULT
      
      ! Unknown preconditioner
      call output_line ('Unknown preconditioner for nonlinear iteration!', &
          OU_CLASS_ERROR,OU_MODE_STD,'cc_initPreconditioner')
      call sys_halt()
      
    end select

  end subroutine


! ***************************************************************************
  subroutine CH_updatePreconditioner (rCHproblem,rnonlinearIteration,&
                              rCHvector,rrhs,binit,bstructuralUpdate,&
	                                  rNSproblem, rNSvector)
!<description>
  ! This routine has to be called whenever the system matrices change.
  ! It initialises (depending on the system matrices) the matrices of the
  ! preconditioner or performs an update of them.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_CHproblem), intent(INOUT), target :: rCHproblem
!</inputoutput>

!<input>
  ! The current solution vector.
  type(t_vectorBlock), intent(IN) :: rCHvector

  ! The right-hand-side vector to use in the equation
  type(t_vectorBlock), intent(IN) :: rrhs

  ! First initialisation.
  ! Has to be set to TRUE on the first call. Initialises the preconditioner
  ! with the structure of the matrices.
  logical, intent(IN) :: binit

  ! Whether the structure of the system matrices is new.
  ! This variable has to be set to TRUE whenever there was a structure in 
  ! the system matrices. This reinitialises the linear solver.
  logical, intent(IN) :: bstructuralUpdate
!</input>

!<inputoutput>
  ! Nonlinar iteration structure saving data for the callback routines.
  ! Preconditioner data is saved here.
  type(t_CHnonlinearIteration), intent(INOUT) :: rnonlinearIteration

!~~~~~~~~~~~~~~~~~~~NS problem~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  type(t_problem), intent(INOUT), optional :: rNSproblem
  type(t_vectorBlock), intent(INOUT), optional :: rNSvector
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: NLMIN,NLMAX
    integer :: i
!MCai, we do not need cmatrixType so far, but if we use different preconditioners,
!We may need it. 
!    integer :: cmatrixType=CHMASM_MTP_FULLTENSOR
!    character(LEN=PARLST_MLDATA) :: snewton
!    type(t_timer) :: rtimer

    ! Error indicator during initialisation of the solver
    integer :: ierror    
  
    ! A pointer to the matrix of the preconditioner
    type(t_matrixBlock), pointer :: p_rmatrixPreconditioner

    type(t_nonlinearCHmatrix) :: rnonlinearCHMatrix
   ! A pointer to the discretisation structure with the data.
    type(t_blockDiscretisation), POINTER :: p_rdiscretisation

    type(t_matrixScalar), pointer :: p_rmatrixTemplateFEM
    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    type(t_matrixBlock), dimension(:), pointer :: Rmatrices


  !~~~~~~~~~~~~~~~~Begin of update preconditioner~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Initialise the matrices.
!MCai, we updatePreconditioner, in which CH_assembleMatrix is called....
      !call CH_updatePreconditioner (rCHproblem,rnonlinearIteration,&
      !    rCHvector,rrhs,.true.,.true.)
!      print *, rnonlinearIteration%rpreconditioner%ctypePreconditioning
      ! Which levels have we to take care of during the solution process?
      NLMIN = rnonlinearIteration%NLMIN
      NLMAX = rnonlinearIteration%NLMAX

      select case (rnonlinearIteration%rpreconditioner%ctypePreconditioning)
      case (CHPREC_NONE)
      ! No preconditioner
      case (CHPREC_LINEARSOLVER,CHPREC_NEWTON,CHPREC_NEWTONDYNAMIC)
      
        ! Which levels have we to take care of during the solution process?
        NLMIN = rnonlinearIteration%NLMIN
        NLMAX = rnonlinearIteration%NLMAX
      
        ! Initialise the preconditioner matrices on all levels.
        !Originally, it is i=NLMIN,...
		!do i=NLMIN,NLMAX
        do i=NLMAX,NLMAX
      
          ! Prepare the preconditioner matrices level i. 
          if (binit .or. bstructuralUpdate) then

!MCai, we do not need cmatrixType, if we use different preconditioners, we may need it.
!            ! What type of matrix do we have? Is it a 'primitive' matrix
!            ! (with A11 and A22 being the same) or a 'full-tensor matrix'
!            ! with Aij all independent from each other
!            if ((rnonlinearIteration%rpreconditioner%ctypePreconditioning .eq. &
!                CHPREC_NEWTON) .or. &
!               (rnonlinearIteration%rpreconditioner%ctypePreconditioning .eq. &
!               CHPREC_NEWTONDYNAMIC)) then
!              cmatrixType = CHMASM_MTP_FULLTENSOR
!            else
!              cmatrixType = CHMASM_MTP_AUTOMATIC
!            end if

            ! Allocate memory for the matrix or release the existing matrix,
            ! if there's a structural update.
            if (binit) then
              allocate(rnonlinearIteration%RcoreEquation(i)%p_rmatrixPreconditioner)
            else
              call lsysbl_releaseMatrix(&
                  rnonlinearIteration%RcoreEquation(i)%p_rmatrixPreconditioner)
            end if 
!~~~~~~~~~~~~~~~~~~Here we need to assemble matrix for preconditioner~~~~~~~~~~~~~~~
! In fact, we do not use cmatrixType here, we replace cc_allocSystemMatrix by the 
! code between '!~~~~~'
        !    call CH_allocSystemMatrix (rCHproblem,rCHproblem%RlevelInfo(i),cmatrixType,&
        !        rnonlinearIteration%RcoreEquation(i)%p_rmatrixPreconditioner)
        !         cc_allocSystemMatrix (rCHproblem,rlevelInfo,cmatrixType,rmatrix)
     
            ! Ask the problem structure to give us the discretisation structure
            p_rdiscretisation => rCHproblem%rlevelInfo(NLMAX)%rdiscretisation
    
            ! Get a pointer to the template FEM matrix.
            p_rmatrixTemplateFEM => rCHproblem%rlevelInfo(NLMAX)%rmatrixTemplateFEM

            ! Initialise the block matrix with default values based on
            ! the discretisation.
            call lsysbl_createMatBlockByDiscr (p_rdiscretisation,& 
                       rnonlinearIteration%RcoreEquation(i)%p_rmatrixPreconditioner) 

            rnonlinearCHMatrix%dalpha = rnonlinearIteration%dalpha   ! A = \alpha M
			! A = \alpha M + \gamma Conv
            rnonlinearCHMatrix%dgamma = rnonlinearIteration%dgamma
			! B = \eta M(c) \Lap
	        rnonlinearCHMatrix%deta = rnonlinearIteration%deta  
            rnonlinearCHMatrix%dbeta =rnonlinearIteration%dbeta      ! (2,2) block, C
			rnonlinearCHMatrix%dtau = rnonlinearIteration%dtau       ! D = \tau N(.) + \theta Lap
            rnonlinearCHMatrix%dtheta = rnonlinearIteration%dtheta     ! D = \tau N(.) + \theta Lap

            rnonlinearCHMatrix%p_rdiscretisation => rCHproblem%rlevelInfo(NLMAX)%rdiscretisation
            rnonlinearCHMatrix%p_rmatrixTemplateFEM => rCHproblem%rlevelInfo(NLMAX)%rmatrixTemplateFEM

!     	    rnonlinearCHMatrix%p_rmatrixA => rCHproblem%rlevelInfo(NLMAX)%rmatrixA
!             rnonlinearCHMatrix%p_rmatrixB => rCHproblem%rlevelInfo(NLMAX)%rmatrixB
! 	        rnonlinearCHMatrix%p_rmatrixC => rCHproblem%rlevelInfo(NLMAX)%rmatrixC
!             rnonlinearCHMatrix%p_rmatrixD => rCHproblem%rlevelInfo(NLMAX)%rmatrixD
            rnonlinearCHMatrix%p_rmatrixMass => rCHproblem%rlevelInfo(NLMAX)%rmatrixMass
            rnonlinearCHMatrix%p_rmatrixLaplace => rCHproblem%rlevelInfo(NLMAX)%rmatrixLaplace
!            rnonlinearCHMatrix%p_rmatrixConv => rCHproblem%rlevelInfo(NLMAX)%rmatrixConv
        
! How to assemble matrix??
!            call cc_assembleMatrix (CCMASM_ALLOCMEM,cmatrixType,&
!             rmatrix,rnonlinearCCMatrix)
!            print *, 'before assemble Matrixt'

            call CH_assembleMatrix (rCHproblem, &
			   rnonlinearIteration%RcoreEquation(i)%p_rmatrixPreconditioner,&
                         rnonlinearCHMatrix,rCHvector, rCHproblem%rparamList, &
						         rNSproblem, rNSvector)

!~~~~~~~~~~~~~~~~~~The above, similar to cc_allocSystemMatrix~~~~~~~~~~~~~~~~~~~~~~~~~
            ! Attach boundary conditions
            rnonlinearIteration%RcoreEquation(i)%p_rmatrixPreconditioner%p_rdiscreteBC &
                => rCHproblem%RlevelInfo(i)%p_rdiscreteBC
         
          end if
        
          ! On the current level, set up a global preconditioner matrix.
          p_rmatrixPreconditioner => &
            rnonlinearIteration%RcoreEquation(i)%p_rmatrixPreconditioner        
        end do
      
!      call stat_clearTimer(rtimer)
!      call stat_startTimer(rtimer)
      
        ! Attach the system matrices to the solver.
        !
        ! For this purpose, copy the matrix structures from the preconditioner
        ! matrices to Rmatrix.
        allocate(Rmatrices(1:NLMAX))
!MCai, originally, it is i=NLMIN, NLMAX
!       do i=NLMIN,NLMAX
        do i=NLMAX,NLMAX
          call lsysbl_duplicateMatrix ( &
            rnonlinearIteration%RcoreEquation(i)%p_rmatrixPreconditioner, &
            Rmatrices(i), LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
        end do
!MCai,.....      
!       call linsol_setMatrices(&
!           rnonlinearIteration%rpreconditioner%p_rsolverNode,Rmatrices(NLMIN:NLMAX))
      
        call linsol_setMatrices(&
           rnonlinearIteration%rpreconditioner%p_rsolverNode,Rmatrices(NLMAX:NLMAX))
		      
       ! The solver got the matrices; clean up Rmatrices, it was only of temporary
       ! nature...
!MCai, ....
!      do i=NLMIN,NLMAX
       do i=NLMAX,NLMAX
         call lsysbl_releaseMatrix (Rmatrices(i))
       end do
       deallocate(Rmatrices)
      
       ! Initialise structure/data of the solver. This allows the
       ! solver to allocate memory / perform some precalculation
       ! to the problem.
       if (binit) then
         call linsol_initStructure (rnonlinearIteration%rpreconditioner%p_rsolverNode,&
             ierror)
         if (ierror .ne. LINSOL_ERR_NOERROR) stop

       else if (bstructuralUpdate) then       
         call linsol_updateStructure (rnonlinearIteration%rpreconditioner%p_rsolverNode,&
             ierror)
         if (ierror .ne. LINSOL_ERR_NOERROR) stop
       end if
      
      ! Gather statistics
!      call stat_stopTimer(rtimer)
!      rCHproblem%rstatistics%dtimeLinearSolverFactorisation = &
!        rCHproblem%rstatistics%dtimeLinearSolverFactorisation + rtimer%delapsedreal
      case DEFAULT
      
      ! Unknown preconditioner
        call output_line ('Unknown preconditioner for nonlinear iteration!', &
          OU_CLASS_ERROR,OU_MODE_STD,'CH_updatePreconditioner')
        call sys_halt()
      
    end select

  end subroutine

! ***************************************************************************

!<subroutine>

  subroutine CH_releasePreconditioner (rnonlinearIteration)
  
!<description>
  ! This routine releases the preconditioner for the nonlinear iteration
  ! which was prepared in CH_initPreconditioner. Memory is released
  ! from heap.
!</description>

!<inputoutput>
  ! Nonlinar iteration structure saving data for the callback routines.
  ! The preconditioner data is removed from that,
  type(t_CHnonlinearIteration), intent(INOUT) :: rnonlinearIteration
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i

    ! Which preconditioner do we have?    
    select case (rnonlinearIteration%rpreconditioner%ctypePreconditioning)
    case (CHPREC_NONE)
      ! No preconditioning
    case (CHPREC_LINEARSOLVER,CHPREC_NEWTON,CHPREC_NEWTONDYNAMIC)
      ! Preconditioner was a linear solver structure.
      !
      ! Release the preconditioner matrix on every level
      !MCai, we only init for NLMAX...
!	  do i=rnonlinearIteration%NLMIN,rnonlinearIteration%NLMAX
	  do i=rnonlinearIteration%NLMAX,rnonlinearIteration%NLMAX
      
        call lsysbl_releaseMatrix ( &
          rnonlinearIteration%RcoreEquation(i)%p_rmatrixPreconditioner)
        deallocate(rnonlinearIteration%RcoreEquation(i)%p_rmatrixPreconditioner)
        
        ! Also release A, B, C, D everywhere where they exist.
!      	  call lsyssc_releaseMatrix (rnonlinearIteration%RcoreEquation(i)%p_rmatrixA)
!         call lsyssc_releaseMatrix (rnonlinearIteration%RcoreEquation(i)%p_rmatrixB)
!         call lsyssc_releaseMatrix (rnonlinearIteration%RcoreEquation(i)%p_rmatrixC)
!         call lsyssc_releaseMatrix (rnonlinearIteration%RcoreEquation(i)%p_rmatrixD)
        
      end do
      
      ! Release the temporary vector(s)
      call lsyssc_releaseVector (rnonlinearIteration%rpreconditioner%p_rtempVectorSc)
      deallocate(rnonlinearIteration%rpreconditioner%p_rtempVectorSc)
      
      call lsyssc_releaseVector (rnonlinearIteration%rpreconditioner%p_rtempVectorSc2)
      deallocate(rnonlinearIteration%rpreconditioner%p_rtempVectorSc2)
      
      ! Clean up data about the projection etc.
      call mlprj_doneProjection(rnonlinearIteration%rpreconditioner%p_rprojection)
      deallocate(rnonlinearIteration%rpreconditioner%p_rprojection)

      ! Clean up the linear solver, release all memory, remove the solver node
      ! from memory.
      call linsol_doneStructure (rnonlinearIteration%rpreconditioner%p_rsolverNode)
      call CH_doneLinearSolver (rnonlinearIteration)
      
    case DEFAULT
      
      ! Unknown preconditioner
      call output_line ('Unknown preconditioner for nonlinear iteration!', &
          OU_CLASS_ERROR,OU_MODE_STD,'CH_releasePreconditioner')
      call sys_halt()
      
    end select

  end subroutine


end module
