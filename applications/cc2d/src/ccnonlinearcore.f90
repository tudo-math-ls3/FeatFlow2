!##############################################################################
!# ****************************************************************************
!# <name> ccnonlinearcore </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains the routines to solve the stationary core equation
!# of the problem with a nonlinear solver.
!# (This is comparable to the NSDEF routine in older FEAT versions)
!#
!# The discretised core equation reads at the moment:
!#
!#  $$        A_1 y   +  \eta B p   = f_1 $$
!#  $$   \tau B^T y                 = f_2 $$
!#
!# with
!#
!#   $$ A_1 = \alpha M  +  \theta L  +  \gamma N(y) $$
!#
!# and
!#
!#   $M$     = mass matrix,
!#   $L$     = Stokes matrix ($\nu$*Laplace),
!#   $N(y)$  = Nonlinearity includung stabilisation,
!#
!#   $\alpha$ = 0/1     - switches the mass matrix on/off;
!#                          =0 for stationary problem,
!#   $\theta$           - weight for the Laplace matrix,
!#   $\gamma$ = 0/1     - Switches the nonlinearity on/off;
!#                          =0 for Stokes system,
!#   $\eta$   = 0/1     - Switches the 'B'-term on/off,
!#   $\tau$   = 0/1     - Switches the 'B^T'-term on/off,
!#
!# (y,p) is the velocity/pressure solution pair.
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
!#  1.) cc_getNonlinearSolver
!#      -> Initialise nonlinear solver configuration with information
!#         from INI/DAT files
!#
!#  2) cc_solveCoreEquation
!#      -> Starts the nonlinear iteration to solve the core equation.
!#
!# Callback routines for the nonlinear solver:
!#
!#  6.) cc_getDefect
!#      -> Callback routine. Calculate nonlinear defect
!#
!#  7.) cc_getOptimalDamping
!#     -> Auxiliary routine. Calculate optimal damping parameter
!#
!#  8.) cc_precondDefect
!#      -> Callback routine. Preconditioning of nonlinear defect
!#
!#  9.) cc_resNormCheck
!#      -> Callback routine: Check the residuals for convergence
!#
!# 10.) cc_nonlinearSolver
!#      -> The actual nonlinear solver; similar to NSDEF in old CC versions.
!#
!# Auxiliary routines:
!#
!# 1.) cc_initNonlinMatrix
!#     -> Initialise the basic structure of a nonlinear matrix.
!#
!# 2.) cc_prepareNonlinMatrixAssembly
!#     -> Prepare a nonlinear matrix for the assembly in order to be
!#        used in a preconditioner.
!#
!# To solve a system with the core equation, one has to deal with two
!# main structures. On one hand, one has a nonlinear iteration structure
!# t_nlsolNode from the kernel; this is initialised by cc_getNonlinearSolver.
!# On the other hand, one has to maintain a 'nonlinear iteration structure'
!# of type t_ccNonlinearIteration, which configures the core equation and
!# specifies parameters for the solver how to work.
!#
!# The basic procedure for dealing with the core equation is as follows:
!#
!#  a) cc_getNonlinearSolver  -> Basic initialisation on the nonlinear solver
!#
!#  b) Basic initialisation of the core equation structures
!#
!#  c) Initialise further parameters in the core equation structure manually
!#     (e.g. preconditioner, pointer to matrices, ...).
!#     It is important, that the 'outer' application initialises pointers to
!#     matrices, otherwise nothing will work!
!#     This all has to be done with the nonlinear-iteration-structure directly.
!#
!#  d) Initialises the constants in the core equation structure according to
!#     the concrete situation.
!#
!#  e) cc_solveCoreEquation -> Solve the core equation with the nonlinear
!#                               solver structure from cc_getNonlinearSolver.
!#
!#  f) Release core equation structure
!#
!# Stebs b), c), d) and f) are done in the routines in the file
!# ccnonlinearcoreinit.f90. This file here contains only the very core
!# routines of the nonlinear iteration, leaving all initialisation to
!# ccnonlinearcoreinit.f90.
!#
!# </purpose>
!##############################################################################

module ccnonlinearcore

  use fsystem
  use storage
  use linearsolver
  use boundary
  use linearalgebra
  use cubature
  use matrixfilters
  use vectorfilters
  use discretebc
  use bcassembly
  use triangulation
  use spatialdiscretisation
  use coarsegridcorrection
  use spdiscprojection
  use nonlinearsolver
  use paramlist
  use multilevelprojection
  use filtersupport
  use linearsolverautoinitialise
  use bilinearformevaluation
  use linearformevaluation
  use matrixrestriction
  use matrixmodification
  use trilinearformevaluation
  use matrixio
  use statistics
  use collection
  use convection
  
  use ccmatvecassembly
    
  use cccallback
  
  implicit none
  
!<constants>

!<constantblock description="Preconditioner identifiers for the defect in the nonlinear iteration">

  ! No preconditioning
  integer, parameter :: CCPREC_NONE         = -1

  ! Preconditioning with inverse mass matrix (not yet implemented)
  integer, parameter :: CCPREC_INVERSEMASS   = 0

  ! Preconditioning by linear solver, solving the linearised system
  integer, parameter :: CCPREC_LINEARSOLVER  = 1

  ! Preconditioning by Newton-Iteration
  integer, parameter :: CCPREC_NEWTON        = 2

  ! Preconditioning by dynamic Newton-Iteration (uses defect correction
  ! and switches automatically to Newton if the error is small enough)
  integer, parameter :: CCPREC_NEWTONDYNAMIC = 3

!</constantblock>

!</constants>

  
!<types>

!<typeblock>

  ! This structure controls the Newton iteration -- i.e. the preconditioning
  ! with the Frechet derivative of the Navier--Stokes equation, which
  ! can lead to quadratic covergence of the nonlinear solver.
  ! As Newton works only in the basin of attraction of the solution,
  ! the parameters in this structure allow to define a switching criterion
  ! when to use Newton. In the first couple of iterations, defect correction
  ! is used, while the iteration switches to Newton if the residuum is small
  ! enough.
  ! This block is used if CCPREC_NEWTONDYNAMIC is used as preconditioner.
  type t_ccDynamicNewtonControl
  
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

  ! Preconditioner structure for CCxD. This structure saves the configuration of the
  ! preconditioner that is used during the nonlinear iteration.
  
  type t_ccPreconditioner
  
    ! Type of preconditioner.
    ! This is one of the CCPREC_xxxx flags as defined above (CCPREC_INVERSEMASS for
    ! preconditioning with inverse mass matrix, CCPREC_LINEARSOLVER for solving a linear
    ! system, CCPREC_NEWTON for a Newton iteration,...)
    integer :: ctypePreconditioning = CCPREC_NONE
    
    ! Pointer to linear solver node if a linear solver is the preconditioner.
    ! (Thus, this applies for the defect correction and the Newton preconditioner).
    type(t_linsolNode), pointer :: p_rsolverNode => null()

    ! If the linerar solver contains a multigrid preconditioner, this is a pointer
    ! to the coarse grid solver. Otherwise, the pointer is not associated.
    type(t_linsolNode), pointer :: p_rcgrSolver => null()
    
    ! Configuration block for the adaptive Newton preconditioner.
    ! Is only valid if ctypePreconditioning=CCPREC_NEWTONDYNAMIC!
    type(t_ccDynamicNewtonControl) :: radaptiveNewton

    ! Temporary scalar vector; used for calculating the nonlinear matrix
    ! on lower levels / projecting the solution from higher to lower levels.
    type(t_vectorScalar), pointer :: p_rtempVectorSc => null()

    ! Temporary scalar vector; used for calculating the optimal damping
    ! parameter.
    type(t_vectorScalar), pointer :: p_rtempVectorSc2 => null()

  end type

!</typeblock>

!<typeblock>

  ! This type is used to save some situation specific assembly information
  ! during the setup phase of the nonlinear solver. Here it is noted, if
  ! and whose matrices exist and/or must be assmebled transposed to be
  ! compatible with the preconditioner and more. It is more or less
  ! a collection if different flags.
  type t_ccPreconditionerSpecials
  
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
    ! =2: Simple Jacobi-like VANCA, 2D Navier Stokes problem, general discretisation
    !     (i.e. automatically chooses the best suitable VANCA variant).
    ! =3: Simple Jacobi-like VANCA, 2D Navier Stokes problem, general discretisation
    !     (i.e. automatically chooses the best suitable VANCA variant).
    !     'direct' method, bypassing the defect correction approach.
    !     (-> specialised variant of 8, but faster)
    ! =4: Full VANCA, 2D Navier Stokes problem, general discretisation
    !     (i.e. automatically chooses the best suitable VANCA variant).
    ! =5: Full VANCA, 2D Navier Stokes problem, general discretisation
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
        
    ! This flag is set to .TRUE. the formulation needs an additional
    ! C-block instead of the 0-block.
    logical :: bneedPressureDiagonalBlock = .false.

    ! This flag is set to .TRUE. if there are no Neumann boundary
    ! components. In that case, the pressure matrices of direct
    ! solvers must be changed.
    logical :: bpressureIndefinite = .false.
    
    ! Set to TRUE if the preconditioner needs virtually transposed B matrices
    ! as D matrices on all levels except for the coarse mesh.
    logical :: bneedVirtTransposedD = .false.
    
    ! Set to TRUE if the preconditioner needs virtually transposed B matrices
    ! as D matrices on the coarse mesh.
    logical :: bneedVirtTransposedDonCoarse = .false.
    
  end type

!</typeblock>

!<typeblock>

  ! Represents the core equation on one level of the discretisation.
  ! Collects all information that are necessary to assemble the
  ! (linearised) system matrix and RHS vector.
  type t_cccoreEquationOneLevel
  
    ! Pointer to the discretisation structure of that level
    type(t_blockDiscretisation), pointer :: p_rdiscretisation => null()

    ! Pointer to template information on this level.
    type(t_asmTemplates), pointer :: p_rasmTempl => null()

    ! Pointer to dynamic information on this level.
    type(t_dynamicLevelInfo), pointer :: p_rdynamicInfo => null()

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
  
    ! An interlevel projection structure for changing levels
    type(t_interlevelProjectionBlock), pointer :: p_rprojection => null()

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
  type t_ccNonlinearIteration
  
    ! ALPHA-parameter that controls the weight of the mass matrix in the
    ! core equation. =0.0 for stationary simulations.
    real(DP) :: dalpha = 0.0_DP
    
    ! THETA-parameter that controls the weight of the Stokes matrix
    ! in the core equation. =1.0 for stationary simulations.
    real(DP) :: dtheta = 0.0_DP
    
    ! GAMMA-parameter that controls the weight in front of the
    ! nonlinearity. =1.0 for Navier-Stokes, =0.0 for Stokes equation.
    real(DP) :: dgamma = 0.0_DP

    ! ETA-parameter that switch the B-term on/off.
    real(DP) :: deta = 0.0_DP
    
    ! TAU-parameter that switch the B^T-term on/off
    real(DP) :: dtau = 0.0_DP
    
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
    
    ! A t_ccPreconditioner saving information about the preconditioner.
    type(t_ccPreconditioner) :: rpreconditioner
    
    ! A t_ccPreconditionerSpecials structure that saves information about
    ! special 'tweaks' in matrices such that everything works.
    type(t_ccPreconditionerSpecials) :: rprecSpecials
    
    ! An array of t_cccoreEquationOneLevel structures for all levels
    ! of the discretisation.
    type(t_cccoreEquationOneLevel), dimension(:), pointer :: RcoreEquation => null()
    
    ! Auxiliary variable: Saves the initial defect in the nonlinear iteration
    real(DP), dimension(2) :: DresidualInit = 0.0_DP
    
    ! Auxiliary variable: Saves the last defect in the nonlinear iteration
    real(DP), dimension(2) :: DresidualOld = 0.0_DP

    ! Auxiliary variable: Norm of the relative change = norm of the
    ! preconditioned residual in the nonlinear iteration
    real(DP), dimension(3) :: DresidualCorr = 0.0_DP
    
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

  !<subroutine>
  
    subroutine cc_initNonlinMatrix (rnonlinearCCMatrix,rproblem,&
        rdiscretisation,rasmTempl,rdynamicLevelInfo)
  
  !<description>
    ! Initialises the rnonlinearCCMatrix structure with parameters and pointers
    ! from the main problem and precalculated information.
  !</description>

  !<input>
    ! Global problem structure.
    type(t_problem), intent(in), target :: rproblem
    
    ! Discretisation of the level where the matrix is to be assembled.
    type(t_blockDiscretisation), intent(in), target :: rdiscretisation
    
    ! Template structures on this level.
    type(t_asmTemplates), intent(in), target :: rasmTempl

    ! Dynamic structures on this level. Usually change in every timestep.
    type(t_dynamicLevelInfo), intent(in), target :: rdynamicLevelInfo
  !</input>
  
  !<inputoutput>
    ! Nonlinear matrix structure.
    ! Basic parameters in this structure are filled with data.
    type(t_nonlinearCCMatrix), intent(inout) :: rnonlinearCCMatrix
  !</inputoutput>
               
  !</subroutine>
      
      ! Initialise the matrix assembly structure rnonlinearCCMatrix
      ! with basic global information.
      !
      ! 1.) Model, stabilisation
      rnonlinearCCMatrix%p_rstabilisation => rproblem%rstabilisation
      rnonlinearCCMatrix%p_rphysics => rproblem%rphysics
      
      ! 2.) Pointers to global precalculated matrices.
      rnonlinearCCMatrix%p_rdiscretisation => rdiscretisation
      rnonlinearCCMatrix%p_rasmTempl => rasmTempl
      rnonlinearCCMatrix%p_rdynamicInfo => rdynamicLevelInfo
      
    end subroutine

  ! ***************************************************************************

  !<subroutine>
  
    subroutine cc_prepareNonlinMatrixAssembly (rnonlinearCCMatrix,&
        ilev,nlmin,nlmax,rprecSpecials)
  
  !<description>
    ! Prepares a rnonlinearCCMatrix structure for the assembly according
    ! to a preconditioner. rprecSpecials specifies a couple of preconditioner
    ! flags that configure the shape of the system matrix that the preconditioner
    ! needs.
    !
    ! cc_initNonlinMatrix must have been called prior to this routine to
    ! initialise the basic matrix. cc_prepareNonlinMatrixAssembly will then
    ! add assembly-specific parameters of the preconditioner.
  !</description>

  !<input>
    ! Current assembly level.
    integer, intent(in) :: ilev
    
    ! Minimum assembly level.
    integer, intent(in) :: nlmin
    
    ! Maximum assembly level.
    integer, intent(in) :: nlmax
  
    ! Structure with assembly-specific parameters of the preconditioner.
    type(t_ccPreconditionerSpecials), intent(in) :: rprecSpecials
  !</input>
  
  !<inputoutput>
    ! Nonlinear matrix structure.
    ! Basic parameters in this structure are filled with data.
    type(t_nonlinearCCMatrix), intent(inout) :: rnonlinearCCMatrix
  !</inputoutput>
               
  !</subroutine>
      
      ! Parameters for adaptive matrices for Q1~ with anisotropic elements
      rnonlinearCCMatrix%iadaptiveMatrices = rprecSpecials%iadaptiveMatrices
      rnonlinearCCMatrix%dadmatthreshold = rprecSpecials%dadmatthreshold
      
      ! Depending on the level, we have to set up information about
      ! transposing B-matrices.
      if (ilev .eq. nlmin) then
        rnonlinearCCMatrix%bvirtualTransposedD = rprecSpecials%bneedVirtTransposedDonCoarse
      else
        rnonlinearCCMatrix%bvirtualTransposedD = rprecSpecials%bneedVirtTransposedD
      end if
      
    end subroutine

  ! ***************************************************************************
  ! Callback routines for the nonlinear solver
  ! ***************************************************************************

  !<subroutine>
  
    subroutine cc_getDefect (rproblem,rnonlinearIteration,rsolverNode,ite,rx,rb,rd)
  
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
    integer, intent(in)                           :: ite

    ! Current iteration vector
    type(t_vectorBlock), intent(in),target        :: rx

    ! Right hand side vector of the equation.
    type(t_vectorBlock), intent(in)               :: rb
  !</input>
               
  !<inputoutput>
    ! The problem structure characterising the whole problem.
    type(t_problem), intent(inout)                :: rproblem

    ! Reference to the nonlinear iteration structure that configures the
    ! main nonlinear equation. Intermediate data is changed during the iteration.
    type(t_ccnonlinearIteration), intent(inout)   :: rnonlinearIteration

    ! The nonlinear solver node that configures the solution process.
    type(t_nlsolNode), intent(inout)              :: rsolverNode

    ! Defect vector b-A(x)x. This must be filled with data by the callback routine.
    type(t_vectorBlock), intent(inout)            :: rd
  !</inputoutput>
  
  !</subroutine>
      
      ! The nonlinear iteration structure
      type(t_nonlinearCCMatrix) :: rnonlinearCCMatrix
      integer :: ilvmax
      type(t_filterChain), dimension(:), pointer :: p_RfilterChain
  
      ! a timer structure
      type(t_timer) :: rtimer
                  
      ! DEBUG!!!
      !REAL(DP), DIMENSION(:), POINTER :: p_Ddata,p_Ddata2

      call stat_clearTimer(rtimer)
      call stat_startTimer(rtimer)

      ! Build the nonlinear defect
      call lsysbl_copyVector (rb,rd)
      
      ! DEBUG!!!
      !CALL lsysbl_getbase_double (rx,p_Ddata)
      !CALL lsysbl_getbase_double (rd,p_Ddata2)
      
      ilvmax = rnonlinearIteration%NLMAX
          
      ! Initialise the matrix assembly structure rnonlinearCCMatrix to describe the
      ! matrix we want to have.
      call cc_initNonlinMatrix (rnonlinearCCMatrix,rproblem,&
        rnonlinearIteration%RcoreEquation(ilvmax)%p_rdiscretisation,&
        rnonlinearIteration%RcoreEquation(ilvmax)%p_rasmTempl,&
        rnonlinearIteration%RcoreEquation(ilvmax)%p_rdynamicInfo)

      rnonlinearCCMatrix%dalpha = rnonlinearIteration%dalpha
      rnonlinearCCMatrix%dtheta = rnonlinearIteration%dtheta
      rnonlinearCCMatrix%dgamma = rnonlinearIteration%dgamma
      rnonlinearCCMatrix%deta = rnonlinearIteration%deta
      rnonlinearCCMatrix%dtau = rnonlinearIteration%dtau
      
      call cc_nonlinearMatMul (rnonlinearCCMatrix,rx,rd,-1.0_DP,1.0_DP,rproblem)
      
      p_RfilterChain => rnonlinearIteration%p_RfilterChain
      if (associated(p_RfilterChain)) then
        call filter_applyFilterChainVec (rd, p_RfilterChain)
      end if
    
      call vecfil_discreteNLSlipBCdef (rd)
      
      ! Gather statistics
      call stat_stopTimer(rtimer)
      rproblem%rstatistics%dtimeDefectCalculation = &
        rproblem%rstatistics%dtimeDefectCalculation + rtimer%delapsedReal
      
    end subroutine
    
  ! ***************************************************************************

  !<subroutine>

    subroutine cc_getOptimalDamping (rproblem,rnonlinearIteration,rsolverNode,&
        rd,rx,rb,rtemp1,rtemp2,domega)
  
  !<description>
    ! This subroutine is called inside of the nonlinear loop, to be precise,
    ! inside of cc_precondDefect. It calculates an optiman damping parameter
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
    type(t_vectorBlock), intent(in)               :: rx

    ! Current RHS vector of the nonlinear equation
    type(t_vectorBlock), intent(in)               :: rb

    ! Defect vector b-A(x)x.
    type(t_vectorBlock), intent(in)               :: rd
  !</input>

  !<inputoutput>
    ! The problem structure characterising the whole problem.
    type(t_problem), intent(inout)                :: rproblem

    ! Reference to the nonlinear iteration structure that configures the
    ! main nonlinear equation. Intermediate data is changed during the iteration.
    type(t_ccnonlinearIteration), intent(inout)   :: rnonlinearIteration

    ! The nonlinear solver node that configures the solution process.
    type(t_nlsolNode), intent(inout)              :: rsolverNode

    ! A temporary vector in the structure of rx
    type(t_vectorBlock), intent(inout)            :: rtemp1

    ! A 2nd temporary vector in the structure of rx
    type(t_vectorBlock), intent(inout)            :: rtemp2

    ! Damping parameter. On entry: an initial value given e.g. by the
    ! previous step.
    ! On return: The new damping parameter.
    real(DP), intent(inout)                       :: domega
  !</inputoutput>
  
  !</subroutine>

    ! local variables
    integer :: ilvmax
    real(DP) :: dskv1,dskv2
    type(t_matrixBlock) :: rmatrix
    type(t_timer) :: rtimer
    ! A filter chain for the linear solver
    type(t_filterChain), dimension(:), pointer :: p_RfilterChain
    
    type(t_nonlinearCCMatrix) :: rnonlinearCCMatrix

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
      ! with an appropriate preconditioner C, which we do not care here.
      ! In our case, this iteration system can be written as:
      !
      ! (u1)     (u1)                     ( (f1)   [ A         B1] (u1) )
      ! (u2)  := (u2)  + OMEGA * C^{-1} * ( (f2) - [      A    B2] (u2) )
      ! (p )     (p )                     ( (fp)   [ D1   D2   0 ] (p ) )
      !
      !                                   |------------------------------|
      !                                              = d_n
      !                            |-------------------------------------|
      !                                        = Y = (y1,y2,yp) = rd
      !
      ! with KST1=KST1(u1,u2,p) and Y=rd being the solution from
      ! the Oseen equation with
      !
      !                  [ A         B1 ]
      !    C = T(u_n) =  [      A    B2 ]
      !                  [ D1   D2   0  ]
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
      ! as one can see by linearisation (see p. 170, Turek`s book).
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

      ! Should we discretise the Navier-Stokes nonlinearity?
      ! Would mean to rebuild the system matrix. Otherwise we can reuse
      ! our previous diffusion matrix without rebuilding it.
      
      ! Re-assemble the nonlinear system matrix on the maximum level.
      ! Initialise the matrix assembly structure rnonlinearCCMatrix to describe the
      ! matrix we want to have.
 
      call cc_initNonlinMatrix (rnonlinearCCMatrix,rproblem,&
        rnonlinearIteration%RcoreEquation(ilvmax)%p_rdiscretisation,&
        rnonlinearIteration%RcoreEquation(ilvmax)%p_rasmTempl,&
        rnonlinearIteration%RcoreEquation(ilvmax)%p_rdynamicInfo)
        
      call cc_prepareNonlinMatrixAssembly (rnonlinearCCMatrix,&
        rnonlinearIteration%NLMAX,&
        rnonlinearIteration%NLMIN,rnonlinearIteration%NLMAX,&
        rnonlinearIteration%rprecSpecials)
        
      rnonlinearCCMatrix%dalpha = rnonlinearIteration%dalpha
      rnonlinearCCMatrix%dtheta = rnonlinearIteration%dtheta
      rnonlinearCCMatrix%dgamma = rnonlinearIteration%dgamma
      rnonlinearCCMatrix%deta = rnonlinearIteration%deta
      rnonlinearCCMatrix%dtau = rnonlinearIteration%dtau

      ! Assemble the matrix.
      call cc_assembleMatrix (CCMASM_COMPUTE,CCMASM_MTP_AUTOMATIC,&
          rmatrix,rnonlinearCCMatrix,rproblem,rtemp1)
      
      ! We do not have to implement any boundary conditions into the matrix
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
            OU_CLASS_ERROR,rsolverNode%coutputMode,'cc_getOptimalDamping')
        call output_line ('Is the triangulation ok??? .tri-file destroyed?', &
            OU_CLASS_ERROR,rsolverNode%coutputMode,'cc_getOptimalDamping')
        call output_line ('Boundary conditions set up properly?', &
            OU_CLASS_ERROR,rsolverNode%coutputMode,'cc_getOptimalDamping')
        call sys_halt()
      end if
      
      ! Ok, we have the nominator and the denominator. Divide them
      ! by each other to calculate the new OMEGA.
      
      domega = dskv1 / dskv2
      
      ! And make sure it is in the allowed range:
      
      domega = max(rnonlinearIteration%domegamin, &
                   min(rnonlinearIteration%domegamax,domega))
      
      ! That is it, we have our new Omega.

      ! Gather statistics
      call stat_stopTimer(rtimer)
      rproblem%rstatistics%dtimeOptimalCorrection = &
        rproblem%rstatistics%dtimeDefectCalculation + rtimer%delapsedReal
      
    end subroutine

  ! ***************************************************************************

  !<subroutine>

    subroutine cc_precondDefect (rproblem,rnonlinearIteration,rsolverNode,&
        ite,rd,rx,rb,domega,bsuccess)
  
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
    type(t_problem), intent(inout)                :: rproblem

    ! Reference to the nonlinear iteration structure that configures the
    ! main nonlinear equation. Intermediate data is changed during the iteration.
    type(t_ccnonlinearIteration), intent(inout), target   :: rnonlinearIteration

    ! The nonlinear solver node that configures the solution process.
    type(t_nlsolNode), intent(inout)              :: rsolverNode

    ! Number of current iteration.
    integer, intent(in)                           :: ite

    ! Defect vector b-A(x)x. This must be replaced by J^{-1} rd by a preconditioner.
    type(t_vectorBlock), intent(inout)            :: rd

    ! Damping parameter. Is set to rsolverNode%domega (usually = 1.0_DP)
    ! on the first call to the callback routine.
    ! The callback routine can modify this parameter according to any suitable
    ! algorithm to calculate an 'optimal damping' parameter. The nonlinear loop
    ! will then use this for adding rd to the solution vector:
    ! $$ x_{n+1} = x_n + domega*rd $$
    ! domega will stay at this value until it is changed again.
    real(DP), intent(inout)                       :: domega

    ! If the preconditioning was a success. Is normally automatically set to
    ! TRUE. If there is an error in the preconditioner, this flag can be
    ! set to FALSE. In this case, the nonlinear solver breaks down with
    ! the error flag set to 'preconditioner broke down'.
    logical, intent(inout)                        :: bsuccess
  !</inputoutput>
  
  !<input>
    ! Current iteration vector
    type(t_vectorBlock), intent(in), target       :: rx

    ! Current right hand side of the nonlinear system
    type(t_vectorBlock), intent(in), target       :: rb
  !</input>
  
  !</subroutine>
  
    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    type(t_vectorScalar), pointer :: p_rvectorTemp,p_rvectorTemp2
    type(t_vectorBlock) :: rtemp1,rtemp2
    integer :: ierror
    type(t_linsolNode), pointer :: p_rsolverNode,p_rcgrSolver
    integer, dimension(3) :: Cnorms
    
    integer :: i
    real(DP) :: dresInit,dres,dtempdef
    logical :: bassembleNewton
    type(t_matrixBlock), dimension(:), pointer :: Rmatrices
    type(t_ccDynamicNewtonControl), pointer :: p_rnewton
    type(t_filterChain), dimension(:), pointer :: p_RfilterChain

    ! a local timer object
    type(t_timer) :: rtimer

    ! DEBUG!!!
    !real(dp), dimension(:), pointer :: p_vec,p_def,p_da
    !call lsysbl_getbase_double (rd,p_def)
    !call lsysbl_getbase_double (rx,p_vec)

      select case (rnonlinearIteration%rpreconditioner%ctypePreconditioning)
      case (CCPREC_NONE)
        ! No preconditioning
        domega = max(rnonlinearIteration%domegamin, &
                    min(rnonlinearIteration%domegamax,domega))
      case (CCPREC_LINEARSOLVER,CCPREC_NEWTON,CCPREC_NEWTONDYNAMIC)
        ! Preconditioning with a linear solver.
        !
        ! At first, assemble the preconditioner matrices on all levels
        ! and incorporate all boundary conditions.
        !
        ! Should we assemble the Newton part?
        
        bassembleNewton = .false.
        
        if (rnonlinearIteration%rpreconditioner%ctypePreconditioning .eq. &
            CCPREC_NEWTON) then
            
          ! Use Newton in any case.
          bassembleNewton = .true.
          
        else if (rnonlinearIteration%rpreconditioner%ctypePreconditioning .eq. &
            CCPREC_NEWTONDYNAMIC) then
            
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
        call assembleLinsolMatrices (rnonlinearIteration,rproblem%rcollection,&
            bassembleNewton,rx,rnonlinearIteration%NLMIN,rnonlinearIteration%NLMAX)
        
        ! Gather statistics
        call stat_stopTimer(rtimer)
        rproblem%rstatistics%dtimeMatrixAssembly = &
          rproblem%rstatistics%dtimeMatrixAssembly + rtimer%delapsedReal

        call stat_clearTimer(rtimer)
        call stat_startTimer(rtimer)
        
        ! Our 'parent' (the caller of the nonlinear solver) has prepared
        ! a preconditioner node for us (a linear solver with symbolically
        ! factorised matrices). Get this from the collection.
      
        p_rsolverNode => rnonlinearIteration%rpreconditioner%p_rsolverNode

        if (rnonlinearIteration%rpreconditioner%ctypePreconditioning .eq. &
            CCPREC_NEWTONDYNAMIC) then
          ! Do we have to apply the inexact Newton?
          if (p_rnewton%cinexactNewton .ne. 0) then
          
            ! Adaptive stopping criterion / inexact Newton active.
            !
            ! Determine the stopping criterion for the linear solver.
            ! This is an adaptive stopping criterion depending on the current
            ! defect. In detail, for the inexact Newton we have
            !
            !   |b-Ax_{i+1}|         ( |b-Ax_i| ) exp             ( |b-Ax_i| )
            !   ------------ = min { ( -------- )     , depsrel * ( -------- ) }
            !     |b-Ax_0|           ( |b-Ax_0| )                 ( |b-Ax_0| )
            !
            ! see e.g. [Michael Hinze, Habilitation, p. 51]
            !
            ! If Newton is not active, we taje the formula
            !
            !   |b-Ax_{i+1}|             ( |b-Ax_i| )
            !   ------------ = depsrel * ( -------- )
            !     |b-Ax_0|               ( |b-Ax_0| )
            !
            ! to always gain depsrel.
            ! Switch off the relative stopping criterion in the linear solver:
            
            p_rsolverNode%istoppingCriterion = 0
            
            ! Just for safetyness, gain at least one digit.
            p_rsolverNode%depsRel = 1.0E-1_DP
            
            ! Calculate the new absolute stopping criterion:
            dresInit = sqrt(rnonlinearIteration%DresidualInit(1)**2 + &
                          rnonlinearIteration%DresidualInit(2)**2)
            dres    = sqrt(rnonlinearIteration%DresidualOld(1)**2 + &
                          rnonlinearIteration%DresidualOld(2)**2)
            
            dtempdef = dres / dresInit
            
            if (bassembleNewton) then
              p_rsolverNode%depsAbs = &
                  min(dtempDef**p_rnewton%dinexactNewtonExponent, &
                      p_rnewton%dinexactNewtonEpsRel*dtempdef) * dresInit
            else
              p_rsolverNode%depsAbs = p_rnewton%dinexactNewtonEpsRel*dtempdef*dresInit
            end if
            
            ! If we have a multigrid solver, we also have to take care for
            ! the coarse grid solver!
            i = rnonlinearIteration%NLMIN
            p_rcgrSolver => rnonlinearIteration%rpreconditioner%p_rcgrSolver
            if (associated(p_rcgrSolver)) then
              ! For the coarse grid solver, we choose the same stopping criterion.
              ! But just for safetyness, the coarse grid solver should gain at least
              ! one digit!
              p_rcgrSolver%istoppingCriterion = 0
              p_rcgrSolver%depsRel = 1.0E-1_DP
              p_rcgrSolver%depsAbs = p_rsolverNode%depsAbs
            end if
            
          end if
          
        end if
          
        ! Re-attach the system matrices to the solver.
        ! Note that no pointers and no handles are changed, so we can savely do
        ! that without calling linsol_doneStructure/linsol_doneStructure.
        ! This simply informs the solver about possible new scaling factors
        ! in the matrices in case they have changed...
        allocate(Rmatrices(rnonlinearIteration%NLMIN:rnonlinearIteration%NLMAX))
        do i=rnonlinearIteration%NLMIN,rnonlinearIteration%NLMAX
          call lsysbl_duplicateMatrix ( &
            rnonlinearIteration%RcoreEquation(i)%p_rmatrixPreconditioner, &
            Rmatrices(i), LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
        end do
        
        call linsol_setMatrices(&
            rnonlinearIteration%rpreconditioner%p_rsolverNode,Rmatrices(:))
            
        ! The solver got the matrices; clean up Rmatrices, it was only of temporary
        ! nature...
        do i=rnonlinearIteration%NLMIN,rnonlinearIteration%NLMAX
          call lsysbl_releaseMatrix (Rmatrices(i))
        end do
        deallocate(Rmatrices)

        ! Initialise data of the solver. This in fact performs a numeric
        ! factorisation of the matrices in UMFPACK-like solvers.
        call linsol_initData (p_rsolverNode, ierror)
        if (ierror .ne. LINSOL_ERR_NOERROR) then
          call output_line ('linsol_initData failed! Matrix singular!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'cc_precondDefect')
          call sys_halt()
        end if
        
        ! Gather statistics
        call stat_stopTimer(rtimer)
        rproblem%rstatistics%dtimeLinearSolverFactorisation = &
          rproblem%rstatistics%dtimeLinearSolverFactorisation + rtimer%delapsedReal
        
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
        rproblem%rstatistics%dtimeLinearSolver = &
          rproblem%rstatistics%dtimeLinearSolver + rtimer%delapsedReal
        rproblem%rstatistics%nlinearIterations = &
          rproblem%rstatistics%nlinearIterations + p_rsolverNode%iiterations
        
        ! Remember convergence rate for output
        rnonlinearIteration%drhoLinearSolver = p_rsolverNode%dconvergenceRate

        ! Release the numeric factorisation of the matrix.
        ! We do not release the symbolic factorisation, as we can use them
        ! for the next iteration.
        call linsol_doneData (p_rsolverNode)
        
        ! Did the preconditioner work?
        bsuccess = p_rsolverNode%iresult .eq. 0
        
      end select
      
      ! Finally calculate a new damping parameter domega.
      if (rnonlinearIteration%rpreconditioner%ctypePreconditioning .ne. CCPREC_NONE) then
        
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
        call cc_getOptimalDamping (rproblem,rnonlinearIteration,rsolverNode,&
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
      ! This is used for the stopping criterium in cc_resNormCheck!
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
    
      subroutine assembleLinsolMatrices (rnonlinearIteration,rcollection,&
          bassembleNewton,rx,NLMIN,NLMAX)

      use linearsystemblock
      use collection

      ! Assembles on every level a matrix for the linear-solver/Newton preconditioner.
      ! bnewton allows to specify whether the Newton matrix or only the standard
      ! system matrix is evaluated. The output is written to the p_rpreconditioner
      ! matrices specified in the rnonlinearIteration structure.

      ! Reference to a collection structure that contains all parameters of the
      ! discretisation (for nonlinearity, etc.).
      type(t_collection), intent(inout)                :: rcollection

      ! Nonlinear iteration structure where to write the linearised system matrix to.
      type(t_ccnonlinearIteration), intent(inout)      :: rnonlinearIteration

      ! TRUE  = Assemble the Newton preconditioner.
      ! FALSE = Assemble the standard defect correction preconditioner
      !         (i.e. the linearised system matrix).
      logical, intent(in) :: bassembleNewton
      
      ! Minimum level of the preconditioner that is to be initialised.
      integer, intent(in)                              :: NLMIN
      
      ! Maximum level of the preconditioner that is to be initialised.
      ! This must corresponds to the last matrix in Rmatrices.
      integer, intent(in)                              :: NLMAX
      
      ! Current iteration vector.
      type(t_vectorBlock), intent(in), target          :: rx

      ! local variables
      real(DP) :: dnewton
      integer :: ilev,icol
      integer, dimension(1) :: Irows = (/1/)
      type(t_matrixBlock), pointer :: p_rmatrix,p_rmatrixFine
      type(t_vectorScalar), pointer :: p_rvectorTemp
      type(t_vectorBlock), pointer :: p_rvectorFine,p_rvectorCoarse
      type(t_nonlinearCCMatrix) :: rnonlinearCCMatrix
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
        p_rvectorTemp => rnonlinearIteration%rpreconditioner%p_rtempVectorSc

        ! Get the filter chain. We need tghat later to filter the matrices.
        p_RfilterChain => rnonlinearIteration%p_RfilterChain

        ! On all levels, we have to set up the nonlinear system matrix,
        ! so that the linear solver can be applied to it.
        
        nullify(p_rmatrix)

        do ilev=NLMAX,NLMIN,-1
        
          ! Get the matrix on the current level.
          ! Shift the previous matrix to the pointer of the fine grid matrix.
          p_rmatrixFine => p_rmatrix
          p_rmatrix => rnonlinearIteration%RcoreEquation(ilev)%p_rmatrixPreconditioner
        
          ! On the highest level, we use rx as solution to build the nonlinear
          ! matrix. On lower levels, we have to create a solution
          ! on that level from a fine-grid solution before we can use
          ! it to build the matrix!
          if (ilev .eq. NLMAX) then
          
            p_rvectorCoarse => rx
            
          else
            ! We have to discretise a level hierarchy and are on a level < NLMAX.
            
            ! Get the projection structure for this level.
            p_rprojection => rnonlinearIteration%RcoreEquation(ilev+1)%p_rprojection

            ! Get the temporary vector on level i. Will receive the solution
            ! vector on that level.
            p_rvectorCoarse => rnonlinearIteration%RcoreEquation(ilev)%p_rtempVector
            
            ! Get the solution vector on level i+1. This is either the temporary
            ! vector on that level, or the solution vector on the maximum level.
            if (ilev .lt. NLMAX-1) then
              p_rvectorFine => rnonlinearIteration%RcoreEquation(ilev+1)%p_rtempVector
            else
              p_rvectorFine => rx
            end if

            ! Interpolate the solution from the finer grid to the coarser grid.
            ! The interpolation is configured in the interlevel projection
            ! structure we got from the collection.
            call mlprj_performInterpolation (p_rprojection,p_rvectorCoarse, &
                                             p_rvectorFine,p_rvectorTemp)

            ! Apply the filter chain to the temp vector.
            ! This implements the boundary conditions that are attached to it.
            ! NOTE: Deactivated for standard CC2D compatibility -- and because
            ! it has to be checked whether the correct boundary conditions
            ! are attached to that vector!
            ! CALL filter_applyFilterChainVec (p_rvectorCoarse, p_RfilterChain)

          end if

          ! Initialise the matrix assembly structure rnonlinearCCMatrix to describe the
          ! matrix we want to have.
          
          call cc_initNonlinMatrix (rnonlinearCCMatrix,rproblem,&
            rnonlinearIteration%RcoreEquation(ilev)%p_rdiscretisation,&
            rnonlinearIteration%RcoreEquation(ilev)%p_rasmTempl,&
            rnonlinearIteration%RcoreEquation(ilev)%p_rdynamicInfo)
            
          call cc_prepareNonlinMatrixAssembly (rnonlinearCCMatrix,&
            ilev,nlmin,nlmax,rnonlinearIteration%rprecSpecials)
          
          rnonlinearCCMatrix%dalpha = rnonlinearIteration%dalpha
          rnonlinearCCMatrix%dtheta = rnonlinearIteration%dtheta
          rnonlinearCCMatrix%dgamma = rnonlinearIteration%dgamma
          if (bassembleNewton) rnonlinearCCMatrix%dnewton = rnonlinearIteration%dgamma
          rnonlinearCCMatrix%deta = rnonlinearIteration%deta
          rnonlinearCCMatrix%dtau = rnonlinearIteration%dtau

          ! Assemble the matrix.
          ! If we are on a lower level, we can specify a 'fine-grid' matrix.
          if (ilev .eq. NLMAX) then
            call cc_assembleMatrix (CCMASM_COMPUTE,CCMASM_MTP_AUTOMATIC,&
                p_rmatrix,rnonlinearCCMatrix,rproblem,p_rvectorCoarse)
          else
            call cc_assembleMatrix (CCMASM_COMPUTE,CCMASM_MTP_AUTOMATIC,&
                p_rmatrix,rnonlinearCCMatrix,rproblem,p_rvectorCoarse,p_rmatrixFine)
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
            ! Call the matrix filter for the boundary conditions to include the BC`s
            ! into the matrix.
            call matfil_discreteBC (p_rmatrix)
            call matfil_discreteFBC (p_rmatrix)
          end if
            
          ! 'Nonlinear' boundary conditions like slip boundary conditions
          ! are not implemented with a filter chain into a matrix.
          ! Call the appropriate matrix filter of 'nonlinear' boundary
          ! conditions manually:
          call matfil_discreteNLSlipBC (p_rmatrix,.true.)
            
        end do
        
        if (rnonlinearIteration%rprecSpecials%bpressureIndefinite) then
          
          ! The 3,3-matrix must exist! This is ensured by the initialisation routine.
          !
          ! We have a pure Dirichlet problem. This may give us some difficulties
          ! in the case, the preconditioner uses a direct solver (UMFPACK).
          ! In this case, we have to include a unit vector to the pressure
          ! matrix to make the problem definite!
          if (rnonlinearIteration%rprecSpecials%isolverType .eq. 0) then
            ! Include a unit vector
            call mmod_replaceLinesByUnitBlk (&
                rnonlinearIteration%RcoreEquation(NLMAX)%p_rmatrixPreconditioner,3,Irows)
          end if
          
          if ((rnonlinearIteration%rprecSpecials%isolverType .eq. 1) .or. &
              (rnonlinearIteration%rprecSpecials%isolverType .eq. 2)) then
          
            ! If we have a MG solver, We also check the coarse grid solver for
            ! the same thing!
            ! What we do not check is the smoother, thus we assume that smoothers
            ! are always solvers that allow the applicance of a filter chain.
            if (rnonlinearIteration%rprecSpecials%icoarseGridSolverType .eq. 0) then
              ! Include a unit vector
              call mmod_replaceLinesByUnitBlk (&
                  rnonlinearIteration%RcoreEquation(NLMIN)%p_rmatrixPreconditioner,3,Irows)
            end if
            
          end if
            
        else
          
          ! The 3,3-block must be a zero-matrix. So if it is present, clear it.
          if (lsysbl_isSubmatrixPresent(p_rmatrix,3,3)) &
            call lsyssc_clearMatrix (p_rmatrix%RmatrixBlock(3,3))
          
        end if

      end subroutine
      
    end subroutine

  ! ***************************************************************************

    subroutine cc_resNormCheck (rproblem,rnonlinearIteration,rsolverNode,&
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
    type(t_problem), intent(inout)                :: rproblem

    ! Reference to the nonlinear iteration structure that configures the
    ! main nonlinear equation. Intermediate data is changed during the iteration.
    type(t_ccnonlinearIteration), intent(inout)   :: rnonlinearIteration

    ! The nonlinear solver node that configures the solution process.
    type(t_nlsolNode), intent(inout)              :: rsolverNode

    ! Number of current iteration. Is set to 0 when the callback routine
    ! is called the first time. In this situation, rd describes the initial
    ! defect vector.
    integer, intent(in)                           :: ite
  
    ! Current iteration vector
    type(t_vectorBlock), intent(in), target       :: rx

    ! Right hand side vector of the equation.
    type(t_vectorBlock), intent(in), target       :: rb

    ! Defect vector b-A(x)x.
    type(t_vectorBlock), intent(in), target       :: rd
  !</inputoutput>
  
  !<output>
    ! Must be set to TRUE by the callback routine if the residuum rd
    ! is within a desired tolerance, so that the solver should treat
    ! the iteration as 'converged'.
    logical, intent(out)                        :: bconvergence
  
    ! Must be set to TRUE by the callback routine if the residuum rd
    ! is out of a desired tolerance, so that the solver should treat
    ! the iteration as 'diverged'.
    logical, intent(out)                        :: bdivergence
  !</output>

      ! local variables
      real(DP), dimension(3) :: Dresiduals
      real(DP) :: dresOld,drhoNL,ddelP,ddelU,dtmp,dresU,dresDIV,dres,dresINIT
      real(DP) :: depsD,depsDiv,depsUR,depsPR,depsRES
      integer, dimension(3) :: Cnorms

      ! Calculate norms of the solution/defect vector
      call cc_getDefectNorm (rx,rb,rd,Dresiduals)
      
      Dresiduals(3) = sqrt(Dresiduals(1)**2 + Dresiduals(2)**2)

      ! In the first iteration (initial defect), print the norm of the defect
      ! and save the norm of the initial residuum to the structure
      if (ite .eq. 0) then
      
        call output_separator (OU_SEP_MINUS,coutputMode=rsolverNode%coutputMode)
        call output_line (' IT  RELU     RELP     DEF-U    DEF-DIV'// &
                          '  DEF-TOT  RHONL    OMEGNL   RHOMG',&
                          coutputMode=rsolverNode%coutputMode)
        call output_separator (OU_SEP_MINUS,coutputMode=rsolverNode%coutputMode)
        call output_line ('  0                   '// &
            trim(sys_sdEP(Dresiduals(1),9,2))//&
            trim(sys_sdEP(Dresiduals(2),9,2))//&
            trim(sys_sdEP(Dresiduals(3),9,2)),coutputMode=rsolverNode%coutputMode)
        call output_separator (OU_SEP_MINUS,coutputMode=rsolverNode%coutputMode)

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
        
        ! Initial defect
        dresINIT = sqrt(rnonlinearIteration%DresidualInit(1)**2 + &
                        rnonlinearIteration%DresidualInit(2)**2)
                        
        ! dresInit=0 may hardly occur -- except when we expect 'no flow'.
        ! But to prevent a check against "something<=0" in this case below,
        ! set dresInit to something <> 0.
        if (dresINIT .eq. 0.0_DP) dresINIT = 1.0_DP

        ! Replace the 'old' residual by the current one
        rnonlinearIteration%DresidualOld(1:2) = Dresiduals(1:2)

        ! Nonlinear convergence rate
        drhoNL = (Dresiduals(3)/dresINIT) ** (1.0_DP/real(ite,DP))
        
        ! Calculate norms of the solution/defect vector, calculated above
        dresU   = Dresiduals(1)
        dresDIV = Dresiduals(2)
        dres    = sqrt(dresU**2 + dresDIV**2)
        
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
        ! The other MAX-norms have to be calculated from U...
        
        Cnorms(:) = LINALG_NORMMAX
        call lsysbl_vectorNormBlock (rx,Cnorms,Dresiduals)

        dtmp = max(Dresiduals(1),Dresiduals(2))
        if (dtmp .lt. 1.0E-8_DP) dtmp = 1.0_DP
        ddelU = max(rnonlinearIteration%DresidualCorr(1),&
                    rnonlinearIteration%DresidualCorr(2))/dtmp
        
        dtmp = Dresiduals(3)
        if (dtmp .lt. 1.0E-8_DP) dtmp = 1.0_DP
        ddelP = rnonlinearIteration%DresidualCorr(3)/dtmp
        
        ! Check if the nonlinear iteration can prematurely terminate.
        !
        ! Get the stopping criteria from the parameters.
        ! Use the DepsNL data according to the initialisation above.
        depsD   = rnonlinearIteration%DepsNL(1)
        depsDiv = rnonlinearIteration%DepsNL(2)
        depsUR  = rnonlinearIteration%DepsNL(3)
        depsPR  = rnonlinearIteration%DepsNL(4)
        depsRES = rnonlinearIteration%DepsNL(5)*dresINIT ! -> ddampingD
        
        ! All residual information calculated.
        ! Check for divergence; use a 'NOT' for better NaN handling.
        bdivergence = .not. (dres/dresINIT .lt. 1E5)
        
        ! Check for convergence
        if((ddelU .le. depsUR).and.(ddelP .le. depsPR)   .and. &
           (dresU .le. depsD) .and.(dresDiv .le. depsDiv).and. &
           (dres .le. depsRES)) then
          bconvergence = .true.
        else
          bconvergence = .false.
        end if

        if ((ddelU .lt. SYS_EPSREAL_DP*1E2_DP) .and. &
            (ddelP .lt. SYS_EPSREAL_DP*1E2_DP)) then
          ! We are hard on machine exactness, so stop the iteraton
          bconvergence =.true.
        end if
        
        ! Print residual information
        call output_line ( &
            trim(sys_si(ite,3))//' '//&
            trim(sys_sdEP(ddelU,9,2))// &
            trim(sys_sdEP(ddelP,9,2))// &
            trim(sys_sdEP(dresU,9,2))// &
            trim(sys_sdEP(dresDIV,9,2))// &
            trim(sys_sdEP(dres,9,2))// &
            trim(sys_sdEP(drhoNL,9,2))// &
            trim(sys_sdEP(rnonlinearIteration%domegaNL,9,2))// &
            trim(sys_sdEP(rnonlinearIteration%drhoLinearSolver,9,2)), &
            coutputMode=rsolverNode%coutputMode)
        
      end if
      
    end subroutine

! ***************************************************************************

!<subroutine>

  subroutine cc_getDefectNorm (rvector,rrhs,rdefect,Dresiduals)
  
!<description>
  ! Calculates a couple of norms from a given solution, defect and RHS vector
  ! of the nonlinear system. This can be used to check convergence criteria
  ! etc.
!</description>

!<input>
  ! The solution vector which is modified later during the nonlinear iteration.
  type(t_vectorBlock), intent(in) :: rvector

  ! The right-hand-side vector to use in the equation
  type(t_vectorBlock), intent(in) :: rrhs
  
  ! A defect vector calculated with rvector and rrhs
  type(t_vectorBlock), intent(in) :: rdefect
!</input>

!<output>
  ! An array receiving different defect norms calculated by the above vectors.
  ! Dresiduals(1) = RESU   = ||defect_u|| / ||rhs|| = velocity residual
  ! Dresiduals(2) = RESDIV = ||p|| / ||u||          = divergence residual
  real(DP), dimension(:), intent(out) :: Dresiduals
!</output>

!</subroutine>

    ! local variables
    real(DP) :: dresF,DresTmp(2),dnormU
    integer, dimension(2) :: Cnorms

    !-----------------------------------------------------------------------
    !     Compute the relative l2-norms  RESU,RESDIV
    !-----------------------------------------------------------------------

    Cnorms(:) = LINALG_NORMEUCLID

    ! RESF := max ( ||F1||_E , ||F2||_E )

    call lsysbl_vectorNormBlock (rrhs,Cnorms,DresTmp)
    dresF = max(DresTmp(1),DresTmp(2))
    if (dresF .lt. 1.0E-8_DP) dresF = 1.0_DP

    !               || (D1,D2) ||_E
    ! RESU = -----------------------------
    !        max ( ||F1||_E , ||F2||_E )

    call lsysbl_vectorNormBlock (rdefect,Cnorms,DresTmp)
    Dresiduals(1) = sqrt(DresTmp(1)**2+DresTmp(2)**2)/dresF

    ! DNORMU = || (U1,U2) ||_l2

    call lsysbl_vectorNormBlock (rvector,Cnorms,DresTmp)
    dnormU = sqrt(DresTmp(1)**2+DresTmp(2)**2)
    if (dnormU .lt. 1.0E-8_DP) dnormU = 1.0_DP

    !             || DP ||_E
    ! RESDIV = ----------------
    !          || (U1,U2) ||_E

    Dresiduals(2) = &
        lsyssc_vectorNorm (rdefect%RvectorBlock(3),LINALG_NORMEUCLID) / dnormU
        
  end subroutine

! ***************************************************************************
! Routines to solve the core equation
! ***************************************************************************

!<subroutine>

  subroutine cc_getNonlinearSolver (rnlSolver, rparamList, sname)
  
!<description>
  ! Creates a nonlinear solver node rnlSolver and initialises it with parameters
  ! from the INI/DAT files given in the parameter list rparamList.
  ! sname is the name of a section in rparamList that configures the
  ! parameter of the nonlinear solver.
!</description>

!<input>
  ! Parameter list that contains the parameters from the INI/DAT file(s).
  type(t_parlist), intent(in) :: rparamList
  
  ! Name of the section in the parameter list containing the parameters
  ! of the nonlinear solver.
  character(LEN=*), intent(in) :: sname
!</input>

!<output>
  ! A t_nlsolNode structure that contains the configuration of the nonlinear
  ! solver. The parameters are initialised according to the information
  ! in the section sname of the parameter list rparamList
  type(t_nlsolNode), intent(inout) :: rnlSolver
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
          OU_CLASS_ERROR,OU_MODE_STD,'cc_getNonlinearSolver')
      call sys_halt()
    end if
    
    ! Parse the given parameters now to initialise the solver node.
    ! This is now CCxD-specific!
    
    call parlst_getvalue_int (p_rsection, 'nminIterations', &
                              rnlSolver%nminIterations, rnlSolver%nminIterations)

    call parlst_getvalue_int (p_rsection, 'nmaxIterations', &
                              rnlSolver%nmaxIterations, rnlSolver%nmaxIterations)

    call parlst_getvalue_int (p_rsection, 'ioutputLevel', &
                              rnlSolver%ioutputLevel, rnlSolver%ioutputLevel)

    call parlst_getvalue_double (p_rsection, 'depsUR', &
                                 rnlSolver%DepsRel(1), rnlSolver%DepsRel(1))
    rnlSolver%DepsRel(2) = rnlSolver%DepsRel(1)

    call parlst_getvalue_double (p_rsection, 'depsPR', &
                                 rnlSolver%DepsRel(3), rnlSolver%DepsRel(3))

    call parlst_getvalue_double (p_rsection, 'depsD', &
                                 rnlSolver%DepsAbs(1), rnlSolver%DepsAbs(1))
    rnlSolver%DepsAbs(2) = rnlSolver%DepsAbs(1)

    call parlst_getvalue_double (p_rsection, 'depsDiv', &
                                 rnlSolver%DepsAbs(3), rnlSolver%DepsAbs(3))

    ! Initial damping parameter.
    call parlst_getvalue_double (p_rsection, 'domegaIni', &
                                 rnlSolver%domega, rnlSolver%domega)

    ! We write out the data of the nonlinear solver to the benchmark
    ! log file as well.
    rnlSolver%coutputMode = OU_MODE_STD+OU_MODE_BENCHLOG

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_nonlinearSolver (rproblem,rnonlinearIteration,&
      rsolverNode,rx,rb,rd)
             
!<description>
  ! This routine invokes the nonlinear defect correction iteration
  !     $$  x_{n+1}  =  x_n  +  J^{-1} ( b - A(x_n) x_n )  $$
  !
  ! It is a modification of the routine nlsol_performSolve in the kernel
  ! to allow passing the problem structure to the different callback routines.
  !
  ! The defect correction loop is split into three tasks:
  !
  ! 1.) Calculation of nonlinear defect
  !         $$  d_n  :=  b - A(x_n)x_n  $$
  !     For this purpose, the routine cc_getDefect is called.
  !
  ! 2.) Preconditioning of the nonlinear defect
  !         $$  u_n  :=  J^{-1} d_n  $$
  !     with a preconditioner $J^{-1}$. The preconditioner is
  !     realised by the routine cc_precondDefect.
  !
  ! 3.) Correction
  !         $$  x_{n+1}  :=  x_n  +  \omega u_n  $$
  !     $\omega$ is usually = 1. The defined routine cc_precondDefect
  !     can modify $\omega$ in every iteration.
  !
  ! The callback routine cc_resNormCheck is called in every iteration.
  ! This routine can check the residuum for convergence/divergence and can print
  ! information about the norm of the residuum to screen.
  !
  ! At the beginning of the nonlinear iteration, the routines
  ! cc_getResiduum and cc_resNormCheck are called once with ite=0 to calculate
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
  type(t_problem), intent(inout)                :: rproblem

  ! Reference to the nonlinear iteration structure that configures the
  ! main nonlinear equation. Intermediate data is changed during the iteration.
  type(t_ccnonlinearIteration), intent(inout)   :: rnonlinearIteration

  ! The nonlinear solver node that configures the solution process.
  type(t_nlsolNode), intent(inout)              :: rsolverNode
  
  ! INPUT: Initial solution vector.
  ! OUTPUT: Final iteration vector.
  type(t_vectorBlock), intent(inout)            :: rx
             
  ! Temporary vector. Must be of the same size/type as rx/rb.
  type(t_vectorBlock), intent(inout)            :: rd
!</inputoutput>
  
!<input>
  ! Right hand side vector of the equation.
  type(t_vectorBlock), intent(in)               :: rb
!</input>

!</subroutine>

  ! local variables
  integer :: ite
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
    
    ! Calculate the initial nonlinear defect to rd:   d = b-A(x)x
    call cc_getDefect (rproblem,rnonlinearIteration,rsolverNode,0,rx,rb,rd)
    
    ite = 0
    rsolverNode%icurrentIteration = ite

    ! The standard convergence/divergence test supports only up to
    ! NLSOL_MAXEQUATIONSERROR equations.
    nblocks = min(rb%nblocks,NLSOL_MAXEQUATIONSERROR)
    
    ! Initial test for convergence/divergence.
    call cc_resNormCheck (rproblem,rnonlinearIteration,rsolverNode,&
        ite,rx,rb,rd,bconvergence,bdivergence)
        
    ! Get the initial residuum; cc_resNormCheck saved that to DresidualOld.
    rsolverNode%DinitialDefect(1) = sqrt(rnonlinearIteration%DresidualInit(1)**2 + &
                                         rnonlinearIteration%DresidualInit(2)**2)
    rsolverNode%DfinalDefect(1) = rsolverNode%DinitialDefect(1)

    ! Perform at least nminIterations iterations
    if (ite .lt. rsolverNode%nminIterations) bconvergence = .false.
  
    ! Check for divergence
    if (bdivergence) rsolverNode%iresult = 1
      
    if ((.not. bconvergence) .and. (.not. bdivergence)) then
    
      ! Let us do the nonlinear loop...
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

        ! Call cc_precondDefect to do the preconditioning.
        ! The routine is allowed to change domega during the
        ! iteration if necessary. The nonlinear solver here does not touch
        ! domega anymore, so the callback routine is the only one changing it.
        call cc_precondDefect (rproblem,rnonlinearIteration,rsolverNode,&
            ite,rd,rx,rb,domega,bsuccess)
        
        ! If bsuccess=false, the preconditioner had an error.
        if (.not. bsuccess) then
          if (rsolverNode%ioutputLevel .ge. 0) then
            call output_line ('NLSOL: Iteration '//&
                trim(sys_siL(ite,10))//' canceled as the preconditioner went down!',&
                coutputMode=rsolverNode%coutputMode)
          end if
          rsolverNode%iresult = 3
          exit
        end if
        
        ! If domega=0.0, the solution vector would stay unchanged. In this
        ! case, the nonlinear solver would not proceed at all, and the next
        ! iteration would behave exactly as before!
        ! So in this case, there is nothing to do, we can stop the iteration.
        if (domega .eq. 0.0_DP) then
          if (rsolverNode%ioutputLevel .ge. 1) then
            call output_line ('NLSOL: Iteration '//&
                trim(sys_siL(ite,10))//' canceled as there is no progress anymore!',&
                coutputMode=rsolverNode%coutputMode)
          end if
          exit
        else
          ! Add the correction vector in rd to rx;
          ! damped by the damping parameter:           x := x + domega u
          call lsysbl_vectorLinearComb (rd,rx,domega,1.0_DP)
          
          ! Calculate the new nonlinear defect to rd:  d = b-A(x)x
          call cc_getDefect (rproblem,rnonlinearIteration,rsolverNode,ite,rx,rb,rd)

          ! Check the defect for convergence.
          call cc_resNormCheck (rproblem,rnonlinearIteration,rsolverNode,&
              ite,rx,rb,rd,bconvergence,bdivergence)
              
          ! Get the new residual; cc_resNormCheck saved that to DresidualOld.
          rsolverNode%DfinalDefect(1) = sqrt(rnonlinearIteration%DresidualOld(1)**2 + &
                                             rnonlinearIteration%DresidualOld(2)**2)

          ! Perform at least nminIterations iterations
          if (ite .lt. rsolverNode%nminIterations) bconvergence = .false.
        
          ! Check for convergence
          if (bconvergence) exit
        
          ! Check for divergence
          if (bdivergence) then
            if (rsolverNode%ioutputLevel .ge. 0) then
              call output_line ('NLSOL: Iteration '//&
                  trim(sys_siL(ite,10))//' canceled, divergence detected!',&
                  coutputMode=rsolverNode%coutputMode)
            end if
            rsolverNode%iresult = 1
            exit
          end if
        
        end if
        
      end do ! ite
      
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

  subroutine cc_solveCoreEquation (rproblem,rnonlinearIteration,rnlSolver,&
      rvector,rrhs,rtempBlock)
  
!<description>
  ! This routine invokes the nonlinear solver to solve the core equation
  ! as configured in the core equation structure.
!</description>

!<input>
  ! The right-hand-side vector to use in the equation.
  type(t_vectorBlock), intent(in) :: rrhs
!</input>
  
!<inputoutput>
  ! The problem structure characterising the whole problem.
  type(t_problem), intent(inout)                :: rproblem

  ! A nonlinear-iteration structure that configures the core equation to solve.
  ! Can be initialised e.g. by cc_createNonlinearLoop + manual setup of
  ! the coefficients of the terms.
  type(t_ccNonlinearIteration), intent(inout) :: rnonlinearIteration

  ! A nonlinear solver configuration.
  ! Can be initialised e.g. by using cc_getNonlinearSolver.
  type(t_nlsolNode), intent(inout) :: rnlSolver
  
  ! Initial solution vector. Is replaced by the new solution vector.
  type(t_vectorBlock), intent(inout) :: rvector

  ! A temporary block vector for the nonlinear iteration.
  ! OPTIONAL: If not specified, a temporary vector is automatically allocated.
  type(t_vectorBlock), intent(inout), target, optional :: rtempBlock
  
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

    ! Call the nonlinear solver. For preconditioning and defect calculation,
    ! the solver calls our callback routines.
    call cc_nonlinearSolver(rproblem,rnonlinearIteration,rnlSolver,&
        rvector,rrhs,p_rtempBlock)

    ! Gather statistics
    call stat_stopTimer(rtimer)
    rproblem%rstatistics%dtimeNonlinearSolver = &
      rproblem%rstatistics%dtimeNonlinearSolver + rtimer%delapsedReal
      
    rproblem%rstatistics%nnonlinearIterations = &
      rproblem%rstatistics%nnonlinearIterations + rnlSolver%iiterations

    if (.not. present(rtempBlock)) then
      ! Release the temporary vector
      call lsysbl_releaseVector (p_rtempBlock)
      deallocate (p_rtempBlock)
    end if
        
  end subroutine

end module
