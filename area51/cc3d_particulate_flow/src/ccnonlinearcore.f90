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
!#     It's important, that the 'outer' application initialises pointers to
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
    type(t_linsolNode), pointer :: p_rsolverNode
    
    ! An interlevel projection structure for changing levels
    type(t_interlevelProjectionBlock), pointer :: p_rprojection
    
    ! Configuration block for the adaptive Newton preconditioner.
    ! Is only valid if ctypePreconditioning=CCPREC_NEWTONDYNAMIC!
    type(t_ccDynamicNewtonControl) :: radaptiveNewton

    ! Temporary scalar vector; used for calculating the nonlinear matrix
    ! on lower levels / projecting the solution from higher to lower levels.
    type(t_vectorScalar), pointer :: p_rtempVectorSc

    ! Temporary scalar vector; used for calculating the optimal damping
    ! parameter.
    type(t_vectorScalar), pointer :: p_rtempVectorSc2

  end type

!</typeblock>

!<typeblock>

  ! This type is used to save some situation specific assembly information
  ! during the setup phase of the nonlinear solver. Here it's noted, if
  ! and whose matrices exist and/or must be assmebled transposed to be
  ! compatible with the preconditioner and more. It's more or less
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
    ! =2: Simple Jacobi-like VANCA, 3D Navier Stokes problem, general discretisation
    !     (i.e. automatically chooses the best suitable VANCA variant).
    ! =3: Simple Jacobi-like VANCA, 3D Navier Stokes problem, general discretisation
    !     (i.e. automatically chooses the best suitable VANCA variant).
    !     'direct' method, bypassing the defect correction approach.
    !     (-> specialised variant of 8, but faster)
    ! =4: Full VANCA, 3D Navier Stokes problem, general discretisation
    !     (i.e. automatically chooses the best suitable VANCA variant).
    ! =5: Full VANCA, 3D Navier Stokes problem, general discretisation
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
        
    ! This flag is set to .TRUE. if there are no Neumann boundary
    ! components. In that case, the pressure matrices of direct
    ! solvers must be changed.
    logical :: bpressureGloballyIndefinite = .false.
    
  end type

!</typeblock>

!<typeblock>

  ! Represents the core equation on one level of the discretisation.
  ! Collects all information that are necessary to assemble the 
  ! (linearised) system matrix and RHS vector.
  type t_cccoreEquationOneLevel
  
    ! Pointer to the discretisation structure of that level
    type(t_blockDiscretisation), pointer :: p_rdiscretisation => null()

    ! Stokes matrix for that specific level (=nu*Laplace)
    type(t_matrixScalar), pointer :: p_rmatrixStokes => null()

    ! B1-matrix for that specific level. 
    type(t_matrixScalar), pointer :: p_rmatrixB1 => null()

    ! B2-matrix for that specific level. 
    type(t_matrixScalar), pointer :: p_rmatrixB2 => null()

    ! B3-matrix for that specific level. 
    type(t_matrixScalar), pointer :: p_rmatrixB3 => null()

    ! Mass matrix
    type(t_matrixScalar), pointer :: p_rmatrixMass => null()
    
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
  
    ! Velocity coupling matrix $A_{12}$.
    ! Exists only of deformation tensor or Newton iteration is used in the
    ! nonlinear iteration.
    type(t_matrixScalar), pointer :: p_rmatrixVelocityCoupling12 => null()

    ! Velocity coupling matrix $A_{21}$.
    ! Exists only of deformation tensor or Newton iteration is used in the
    ! nonlinear iteration.
    type(t_matrixScalar), pointer :: p_rmatrixVelocityCoupling21 => null()

    ! Velocity coupling matrix $A_{13}$.
    ! Exists only of deformation tensor or Newton iteration is used in the
    ! nonlinear iteration.
    type(t_matrixScalar), pointer :: p_rmatrixVelocityCoupling13 => null()

    ! Velocity coupling matrix $A_{31}$.
    ! Exists only of deformation tensor or Newton iteration is used in the
    ! nonlinear iteration.
    type(t_matrixScalar), pointer :: p_rmatrixVelocityCoupling31 => null()

    ! Velocity coupling matrix $A_{23}$.
    ! Exists only of deformation tensor or Newton iteration is used in the
    ! nonlinear iteration.
    type(t_matrixScalar), pointer :: p_rmatrixVelocityCoupling23 => null()

    ! Velocity coupling matrix $A_{32}$.
    ! Exists only of deformation tensor or Newton iteration is used in the
    ! nonlinear iteration.
    type(t_matrixScalar), pointer :: p_rmatrixVelocityCoupling32 => null()
    
    ! Pointer to a B1^T-matrix.
    ! This pointer may point to NULL(). In this case, B1^T is created
    ! by 'virtually transposing' the B1 matrix.
    !
    ! Note: This information is automatically created when the preconditioner
    ! is initialised! The main application does not have to initialise it!
    type(t_matrixScalar), pointer :: p_rmatrixB1T => null()

    ! Pointer to a B2^T-matrix.
    ! This pointer may point to NULL(). In this case, B2^T is created
    ! by 'virtually transposing' the B2 matrix.
    !
    ! Note: This information is automatically created when the preconditioner
    ! is initialised! The main application does not have to initialise it!
    type(t_matrixScalar), pointer :: p_rmatrixB2T => null()

    ! Pointer to a B3^T-matrix.
    ! This pointer may point to NULL(). In this case, B3^T is created
    ! by 'virtually transposing' the B3 matrix.
    !
    ! Note: This information is automatically created when the preconditioner
    ! is initialised! The main application does not have to initialise it!
    type(t_matrixScalar), pointer :: p_rmatrixB3T => null()

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
    real(DP), dimension(4) :: DresidualCorr = 0.0_DP
    
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
  ! Callback routines for the nonlinear solver
  ! ***************************************************************************

  !<subroutine>
  
    subroutine cc_getDefect (rproblem,rnonlinearIteration,ite,rx,rb,rd)
  
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
    type(t_problem), intent(INOUT)                :: rproblem

    ! Reference to the nonlinear iteration structure that configures the
    ! main nonlinear equation. Intermediate data is changed during the iteration.
    type(t_ccnonlinearIteration), intent(INOUT)   :: rnonlinearIteration

    ! Defect vector b-A(x)x. This must be filled with data by the callback routine.
    type(t_vectorBlock), intent(INOUT)            :: rd
  !</inputoutput>
  
  !</subroutine>
      
      ! The nonlinear iteration structure
      type(t_nonlinearCCMatrix) :: rnonlinearCCMatrix
      integer :: ilvmax
      type(t_filterChain), dimension(:), pointer :: p_RfilterChain
  
      ! a timer structure
      type(t_timer) :: rtimer
                  
      ! Build the nonlinear defect
      call lsysbl_copyVector (rb,rd)
      
      ilvmax = rnonlinearIteration%NLMAX
          
      ! Initialise the matrix assembly structure rnonlinearCCMatrix to describe the
      ! matrix we want to have.
      rnonlinearCCMatrix%dalpha = rnonlinearIteration%dalpha
      rnonlinearCCMatrix%dtheta = rnonlinearIteration%dtheta
      rnonlinearCCMatrix%dgamma = rnonlinearIteration%dgamma
      rnonlinearCCMatrix%deta = 1.0_DP
      rnonlinearCCMatrix%dtau = 1.0_DP
      rnonlinearCCMatrix%iupwind = rproblem%rstabilisation%iupwind
      rnonlinearCCMatrix%dnu = rproblem%dnu
      rnonlinearCCMatrix%dupsam = rproblem%rstabilisation%dupsam
      rnonlinearCCMatrix%clocalH = rproblem%rstabilisation%clocalH
      rnonlinearCCMatrix%p_rdiscretisation => &
          rnonlinearIteration%RcoreEquation(ilvmax)%p_rdiscretisation
      rnonlinearCCMatrix%p_rmatrixStokes => &
          rnonlinearIteration%RcoreEquation(ilvmax)%p_rmatrixStokes
      rnonlinearCCMatrix%p_rmatrixB1 => &
          rnonlinearIteration%RcoreEquation(ilvmax)%p_rmatrixB1
      rnonlinearCCMatrix%p_rmatrixB2 => &
          rnonlinearIteration%RcoreEquation(ilvmax)%p_rmatrixB2
      rnonlinearCCMatrix%p_rmatrixB3 => &
          rnonlinearIteration%RcoreEquation(ilvmax)%p_rmatrixB3
      rnonlinearCCMatrix%p_rmatrixB1T => &
          rnonlinearIteration%RcoreEquation(ilvmax)%p_rmatrixB1T
      rnonlinearCCMatrix%p_rmatrixB2T => &
          rnonlinearIteration%RcoreEquation(ilvmax)%p_rmatrixB2T
      rnonlinearCCMatrix%p_rmatrixB3T => &
          rnonlinearIteration%RcoreEquation(ilvmax)%p_rmatrixB3T
      rnonlinearCCMatrix%p_rmatrixMass => &
          rnonlinearIteration%RcoreEquation(ilvmax)%p_rmatrixMass

      call cc_nonlinearMatMul (rnonlinearCCMatrix,rx,rd,-1.0_DP,1.0_DP)        
      
      p_RfilterChain => rnonlinearIteration%p_RfilterChain
      if (associated(p_RfilterChain)) then    
        call filter_applyFilterChainVec (rd, p_RfilterChain)
      end if
    
      call vecfil_discreteNLSlipBCdef (rd)
      
    end subroutine
    
  ! ***************************************************************************

  !<subroutine>

    subroutine cc_getOptimalDamping (rproblem,rnonlinearIteration,rd,rx,rb,&
        rtemp1,rtemp2,domega)
  
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
    type(t_vectorBlock), intent(IN)               :: rx

    ! Current RHS vector of the nonlinear equation
    type(t_vectorBlock), intent(IN)               :: rb

    ! Defect vector b-A(x)x. 
    type(t_vectorBlock), intent(IN)               :: rd
  !</input>

  !<inputoutput>
    ! The problem structure characterising the whole problem.
    type(t_problem), intent(INOUT)                :: rproblem

    ! Reference to the nonlinear iteration structure that configures the
    ! main nonlinear equation. Intermediate data is changed during the iteration.
    type(t_ccnonlinearIteration), intent(INOUT)   :: rnonlinearIteration

    ! A temporary vector in the structure of rx
    type(t_vectorBlock), intent(INOUT)            :: rtemp1

    ! A 2nd temporary vector in the structure of rx
    type(t_vectorBlock), intent(INOUT)            :: rtemp2

    ! Damping parameter. On entry: an initial value given e.g. by the
    ! previous step.
    ! On return: The new damping parameter.
    real(DP), intent(INOUT)                       :: domega
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

      ! Is there anything to do?
      if (rnonlinearIteration%domegaMin .ge. rnonlinearIteration%domegaMax) then
        ! No - cancel.
        domega = rnonlinearIteration%domegaMin
        return
      end if
      
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
      ! (u1)     (u1)                     ( (f1)   [ A         B1] (u1) )
      ! (u2)  := (u2)  + OMEGA * C^{-1} * ( (f2) - [      A    B2] (u2) )
      ! (p )     (p )                     ( (fp)   [ B1^T B2^T 0 ] (p ) )
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
      !                  [ B1^T B2^T 0  ]
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

      ! Should we discretise the Navier-Stokes nonlinearity?
      ! Would mean to rebuild the system matrix. Otherwise we can reuse
      ! our previous diffusion matrix without rebuilding it.
      
      ! Re-assemble the nonlinear system matrix on the maximum level.
      ! Initialise the matrix assembly structure rnonlinearCCMatrix to describe the
      ! matrix we want to have.
      rnonlinearCCMatrix%dalpha = rnonlinearIteration%dalpha
      rnonlinearCCMatrix%dtheta = rnonlinearIteration%dtheta
      rnonlinearCCMatrix%dgamma = rnonlinearIteration%dgamma
      rnonlinearCCMatrix%deta = 1.0_DP
      rnonlinearCCMatrix%dtau = 1.0_DP
      rnonlinearCCMatrix%iupwind = rproblem%rstabilisation%iupwind
      rnonlinearCCMatrix%dnu = rproblem%dnu
      rnonlinearCCMatrix%dupsam = rproblem%rstabilisation%dupsam
      rnonlinearCCMatrix%clocalH = rproblem%rstabilisation%clocalH
      rnonlinearCCMatrix%iadaptiveMatrices = &
          rnonlinearIteration%rprecSpecials%iadaptiveMatrices
      rnonlinearCCMatrix%dadmatthreshold = &
          rnonlinearIteration%rprecSpecials%dadmatthreshold
      rnonlinearCCMatrix%p_rdiscretisation => &
          rnonlinearIteration%RcoreEquation(ilvmax)%p_rdiscretisation
      rnonlinearCCMatrix%p_rmatrixStokes => &
          rnonlinearIteration%RcoreEquation(ilvmax)%p_rmatrixStokes
      rnonlinearCCMatrix%p_rmatrixB1 => &
          rnonlinearIteration%RcoreEquation(ilvmax)%p_rmatrixB1
      rnonlinearCCMatrix%p_rmatrixB2 => &
          rnonlinearIteration%RcoreEquation(ilvmax)%p_rmatrixB2
      rnonlinearCCMatrix%p_rmatrixB3 => &
          rnonlinearIteration%RcoreEquation(ilvmax)%p_rmatrixB3
      rnonlinearCCMatrix%p_rmatrixMass => &
          rnonlinearIteration%RcoreEquation(ilvmax)%p_rmatrixMass

      ! Assemble the matrix.        
      call cc_assembleMatrix (CCMASM_COMPUTE,CCMASM_MTP_AUTOMATIC,&
          rmatrix,rnonlinearCCMatrix,rtemp1)
      
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
      !
      ! Note: implementation commented out!
      ! Seems to work better without!
      ! CALL vecfil_discreteNLSlipBCdef (rtemp2)
      
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
      
      ! Release the matrix again
      call lsysbl_releaseMatrix (rmatrix)
      
      ! Filter the resulting defect vector through the slip-boundary-
      ! condition vector filter for implementing nonlinear slip boundary
      ! conditions into a defect vector.
      !
      ! Note: implementation commented out!
      ! Seems to work better without!
      ! CALL vecfil_discreteNLSlipBCdef (rtemp1)

      ! ==================================================================
      ! Calculation of the fraction terms.
      ! Calculate nominator:    dskv1:= (T*Y,D)   = (rtemp1,rtemp2)
      ! Calculate denominator:  dskv2:= (T*Y,T*Y) = (rtemp1,rtemp1)
      ! ==================================================================
      
      dskv1 = lsysbl_scalarProduct (rtemp1, rtemp2)
      dskv2 = lsysbl_scalarProduct (rtemp1, rtemp1)
      
      if (dskv2 .lt. 1.0E-40_DP) then
        call output_line ('dskv2 nearly zero. Optimal damping parameter singular.', &
            OU_CLASS_ERROR,OU_MODE_STD,'cc_getOptimalDamping')
        call output_line ('Is the triangulation ok??? .tri-file destroyed?', &
            OU_CLASS_ERROR,OU_MODE_STD,'cc_getOptimalDamping')
        call output_line ('Boundary conditions set up properly?', &
            OU_CLASS_ERROR,OU_MODE_STD,'cc_getOptimalDamping')
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
      
    end subroutine

  ! ***************************************************************************

  !<subroutine>

    subroutine cc_precondDefect (rproblem,rnonlinearIteration,&
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
    type(t_problem), intent(INOUT)                :: rproblem

    ! Reference to the nonlinear iteration structure that configures the
    ! main nonlinear equation. Intermediate data is changed during the iteration.
    type(t_ccnonlinearIteration), intent(INOUT), target   :: rnonlinearIteration

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
  !</input>
  
  !</subroutine>
  
    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    type(t_vectorScalar), pointer :: p_rvectorTemp,p_rvectorTemp2
    type(t_vectorBlock) :: rtemp1,rtemp2
    integer :: ierror
    type(t_linsolNode), pointer :: p_rsolverNode
    integer, dimension(4) :: Cnorms
    
    integer :: i
    real(DP) :: dresInit,dres
    logical :: bassembleNewton
    type(t_matrixBlock), dimension(:), pointer :: Rmatrices
    type(t_ccDynamicNewtonControl), pointer :: p_rnewton
    type(t_filterChain), dimension(:), pointer :: p_RfilterChain

    ! a local timer object
    type(t_timer) :: rtimer

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
        
        ! Assemble the preconditioner matrices in rnonlinearIteration
        ! on all levels that the solver uses.
        call assembleLinsolMatrices (rnonlinearIteration,rproblem%rcollection,&
            bassembleNewton,rx,rnonlinearIteration%NLMIN,rnonlinearIteration%NLMAX)
        
        ! Our 'parent' (the caller of the nonlinear solver) has prepared
        ! a preconditioner node for us (a linear solver with symbolically
        ! factorised matrices). Get this from the collection.
      
        p_rsolverNode => rnonlinearIteration%rpreconditioner%p_rsolverNode

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
        if (ierror .ne. LINSOL_ERR_NOERROR) stop
        
        ! Finally solve the system. As we want to solve Ax=b with
        ! b being the real RHS and x being the real solution vector,
        ! we use linsol_solveAdaptively. If b is a defect
        ! RHS and x a defect update to be added to a solution vector,
        ! we would have to use linsol_precondDefect instead.
        call linsol_precondDefect (p_rsolverNode,rd)
        
        ! Remember convergence rate for output
        rnonlinearIteration%drhoLinearSolver = p_rsolverNode%dconvergenceRate

        ! Release the numeric factorisation of the matrix.
        ! We don't release the symbolic factorisation, as we can use them
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
        call cc_getOptimalDamping (rproblem,rnonlinearIteration,&
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
      type(t_collection), intent(INOUT)                :: rcollection

      ! Nonlinear iteration structure where to write the linearised system matrix to.
      type(t_ccnonlinearIteration), intent(INOUT)      :: rnonlinearIteration

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

      ! local variables
      real(DP) :: dnewton
      integer :: ilev,icol
      integer(PREC_MATIDX), dimension(1) :: Irows = (/1/)
      type(t_matrixBlock), pointer :: p_rmatrix,p_rmatrixFine
      type(t_vectorScalar), pointer :: p_rvectorTemp
      type(t_vectorBlock), pointer :: p_rvectorFine,p_rvectorCoarse
      type(t_nonlinearCCMatrix) :: rnonlinearCCMatrix
      type(t_interlevelProjectionBlock), pointer :: p_rprojection
      type(t_timer) :: rtimer
      ! A filter chain for the linear solver
      type(t_filterChain), dimension(:), pointer :: p_RfilterChain
      
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
          rnonlinearCCMatrix%dalpha = rnonlinearIteration%dalpha
          rnonlinearCCMatrix%dtheta = rnonlinearIteration%dtheta
          rnonlinearCCMatrix%dgamma = rnonlinearIteration%dgamma
          if (bassembleNewton) rnonlinearCCMatrix%dnewton = rnonlinearIteration%dgamma
          rnonlinearCCMatrix%deta = 1.0_DP
          rnonlinearCCMatrix%dtau = 1.0_DP
          rnonlinearCCMatrix%iupwind = rproblem%rstabilisation%iupwind
          rnonlinearCCMatrix%dnu = rproblem%dnu
          rnonlinearCCMatrix%dupsam = rproblem%rstabilisation%dupsam
          rnonlinearCCMatrix%clocalH = rproblem%rstabilisation%clocalH
          rnonlinearCCMatrix%iadaptiveMatrices = &
              rnonlinearIteration%rprecSpecials%iadaptiveMatrices
          rnonlinearCCMatrix%dadmatthreshold = &
              rnonlinearIteration%rprecSpecials%dadmatthreshold
          rnonlinearCCMatrix%p_rdiscretisation => &
              rnonlinearIteration%RcoreEquation(ilev)%p_rdiscretisation
          rnonlinearCCMatrix%p_rmatrixStokes => &
              rnonlinearIteration%RcoreEquation(ilev)%p_rmatrixStokes
          rnonlinearCCMatrix%p_rmatrixB1 => &
              rnonlinearIteration%RcoreEquation(ilev)%p_rmatrixB1
          rnonlinearCCMatrix%p_rmatrixB2 => &
              rnonlinearIteration%RcoreEquation(ilev)%p_rmatrixB2
          rnonlinearCCMatrix%p_rmatrixB3 => &
              rnonlinearIteration%RcoreEquation(ilev)%p_rmatrixB3
          rnonlinearCCMatrix%p_rmatrixB1T => &
              rnonlinearIteration%RcoreEquation(ilev)%p_rmatrixB1T
          rnonlinearCCMatrix%p_rmatrixB2T => &
              rnonlinearIteration%RcoreEquation(ilev)%p_rmatrixB2T
          rnonlinearCCMatrix%p_rmatrixB3T => &
              rnonlinearIteration%RcoreEquation(ilev)%p_rmatrixB3T
          rnonlinearCCMatrix%p_rmatrixMass => &
              rnonlinearIteration%RcoreEquation(ilev)%p_rmatrixMass

          ! Assemble the matrix.
          ! If we are on a lower level, we can specify a 'fine-grid' matrix.
          if (ilev .eq. NLMAX) then
            call cc_assembleMatrix (CCMASM_COMPUTE,CCMASM_MTP_AUTOMATIC,&
                p_rmatrix,rnonlinearCCMatrix,p_rvectorCoarse)
          else
            call cc_assembleMatrix (CCMASM_COMPUTE,CCMASM_MTP_AUTOMATIC,&
                p_rmatrix,rnonlinearCCMatrix,p_rvectorCoarse,p_rmatrixFine)
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
            ! Call the matrix filter for the boundary conditions to include the BC's
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
        
        if (rnonlinearIteration%rprecSpecials%bpressureGloballyIndefinite) then
          
          ! The 4,4-matrix must exist! This is ensured by the initialisation routine.
          !
          ! We have a pure Dirichlet problem. This may give us some difficulties
          ! in the case, the preconditioner uses a direct solver (UMFPACK).
          ! In this case, we have to include a unit vector to the pressure
          ! matrix to make the problem definite!
          if (rnonlinearIteration%rprecSpecials%isolverType .eq. 0) then
            p_rmatrix => rnonlinearIteration%RcoreEquation(NLMAX)%p_rmatrixPreconditioner
            
            ! Include a unit vector
            call mmod_replaceLinesByZero(p_rmatrix%RmatrixBlock(4,1),Irows)
            call mmod_replaceLinesByZero(p_rmatrix%RmatrixBlock(4,2),Irows)
            call mmod_replaceLinesByZero(p_rmatrix%RmatrixBlock(4,3),Irows)
            call mmod_replaceLinesByUnit(p_rmatrix%RmatrixBlock(4,4),Irows)
            
          end if
          
          if (rnonlinearIteration%rprecSpecials%isolverType .eq. 1) then
          
            ! If we have a MG solver, We also check the coarse grid solver for 
            ! the same thing!
            ! What we don't check is the smoother, thus we assume that smoothers
            ! are always solvers that allow the applicance of a filter chain.
            if (rnonlinearIteration%rprecSpecials%icoarseGridSolverType .eq. 0) then
              p_rmatrix => rnonlinearIteration%RcoreEquation(NLMIN)%p_rmatrixPreconditioner
              
              ! Include a unit vector
              call mmod_replaceLinesByZero(p_rmatrix%RmatrixBlock(4,1),Irows)
              call mmod_replaceLinesByZero(p_rmatrix%RmatrixBlock(4,2),Irows)
              call mmod_replaceLinesByZero(p_rmatrix%RmatrixBlock(4,3),Irows)
              call mmod_replaceLinesByUnit(p_rmatrix%RmatrixBlock(4,4),Irows)
              
            end if
            
          end if
            
        else
          
          ! The 4,4-block must be a zero-matrix. So if it's present, clear it.
          if (lsysbl_isSubmatrixPresent(p_rmatrix,4,4)) &
            call lsyssc_clearMatrix (p_rmatrix%RmatrixBlock(4,4))
          
        end if        

      end subroutine
      
    end subroutine

  ! ***************************************************************************

    subroutine cc_resNormCheck (rproblem,rnonlinearIteration,&
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
    type(t_problem), intent(INOUT)                :: rproblem

    ! Reference to the nonlinear iteration structure that configures the
    ! main nonlinear equation. Intermediate data is changed during the iteration.
    type(t_ccnonlinearIteration), intent(INOUT)   :: rnonlinearIteration

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
      real(DP), dimension(4) :: Dresiduals
      real(DP) :: dresOld,drhoNL,ddelP,ddelU,dtmp,dresU,dresDIV,dres,dresINIT
      real(DP) :: depsD,depsDiv,depsUR,depsPR,depsRES
      integer, dimension(4) :: Cnorms

      ! Calculate norms of the solution/defect vector
      call cc_getDefectNorm (rx,rb,rd,Dresiduals)
      Dresiduals(3) = sqrt(Dresiduals(1)**2 + Dresiduals(2)**2 )

      ! In the first iteration (initial defect), print the norm of the defect
      ! and save the norm of the initial residuum to the structure
      if (ite .eq. 0) then
      
        call output_separator (OU_SEP_MINUS)     
        call output_line (' IT  RELU     RELP     DEF-U    DEF-DIV'// &
                          '  DEF-TOT  RHONL    OMEGNL   RHOMG')
        call output_separator (OU_SEP_MINUS)     
        call output_line ('  0                   '// &
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

        ! Nonlinear convergence rate
        drhoNL = (Dresiduals(3)/dresOld) ** (1.0_DP/real(ite,DP))
        
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

        dtmp = max(Dresiduals(1),Dresiduals(2),Dresiduals(3))
        if (dtmp .lt. 1.0E-8_DP) dtmp = 1.0_DP
        ddelU = max(rnonlinearIteration%DresidualCorr(1),&
                    rnonlinearIteration%DresidualCorr(2),&
                    rnonlinearIteration%DresidualCorr(3))/dtmp
        
        dtmp = Dresiduals(4)
        if (dtmp .lt. 1.0E-8_DP) dtmp = 1.0_DP
        ddelP = rnonlinearIteration%DresidualCorr(4)/dtmp
        
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
        depsDiv = rnonlinearIteration%DepsNL(2)
        depsUR  = rnonlinearIteration%DepsNL(3)
        depsPR  = rnonlinearIteration%DepsNL(4)
        depsRES = rnonlinearIteration%DepsNL(5)*dresINIT
        
        ! All residual information calculated.
        ! Check for divergence; use a 'NOT' for better NaN handling.
        bdivergence = .not. (dres/dresOld .lt. 1E2)
        
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
            trim(sys_sdEP(rnonlinearIteration%drhoLinearSolver,9,2)) &
            )
        
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
  type(t_vectorBlock), intent(IN) :: rvector

  ! The right-hand-side vector to use in the equation
  type(t_vectorBlock), intent(IN) :: rrhs
  
  ! A defect vector calculated with rvector and rrhs
  type(t_vectorBlock), intent(IN) :: rdefect
!</input>

!<output>
  ! An array receiving different defect norms calculated by the above vectors.
  ! Dresiduals(1) = RESU   = ||defect_u|| / ||rhs|| = velocity residual
  ! Dresiduals(2) = RESDIV = ||p|| / ||u||          = divergence residual
  real(DP), dimension(:), intent(OUT) :: Dresiduals
!</output>

!</subroutine>

    ! local variables
    real(DP) :: dresF,DresTmp(3),dnormU
    integer, dimension(3) :: Cnorms

    !-----------------------------------------------------------------------
    !     Compute the relative l2-norms  RESU,RESDIV
    !-----------------------------------------------------------------------

    Cnorms(:) = LINALG_NORMEUCLID

    ! RESF := max ( ||F1||_E , ||F2||_E, ||F3||_E )

    call lsysbl_vectorNormBlock (rrhs,Cnorms,DresTmp)
    dresF = max(DresTmp(1),DresTmp(2),DresTmp(3))
    if (dresF .lt. 1.0E-8_DP) dresF = 1.0_DP

    !               || (D1,D2,D3) ||_E
    ! RESU = -----------------------------
    !        max ( ||F1||_E , ||F2||_E, ||F3||_E )

    call lsysbl_vectorNormBlock (rdefect,Cnorms,DresTmp)
    Dresiduals(1) = sqrt(DresTmp(1)**2+DresTmp(2)**2+DresTmp(3)**2)/dresF

    ! DNORMU = || (U1,U2,U3) ||_l2 

    call lsysbl_vectorNormBlock (rvector,Cnorms,DresTmp)
    dnormU = sqrt(DresTmp(1)**2+DresTmp(2)**2+DresTmp(3)**2)
    if (dnormU .lt. 1.0E-8_DP) dnormU = 1.0_DP

    !             || DP ||_E
    ! RESDIV = ----------------
    !          || (U1,U2,U3) ||_E

    Dresiduals(2) = &
        lsyssc_vectorNorm (rdefect%RvectorBlock(4),LINALG_NORMEUCLID) / dnormU
        
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

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_nonlinearSolver (rproblem,rnonlinearIteration,&
      rsolverNode,rx,rb,rd)
             
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
  type(t_problem), intent(INOUT)                :: rproblem

  ! Reference to the nonlinear iteration structure that configures the
  ! main nonlinear equation. Intermediate data is changed during the iteration.
  type(t_ccnonlinearIteration), intent(INOUT)   :: rnonlinearIteration

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
!</input>

!</subroutine>

  ! local variables
  type(t_collection), pointer :: p_rcollection
  integer :: ite
  real(DP), dimension(NLSOL_MAXEQUATIONSERROR) :: DvecNorm
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
    call cc_getDefect (rproblem,rnonlinearIteration,0,rx,rb,rd)
    
    ite = 0
    rsolverNode%icurrentIteration = ite

    ! The standard convergence/divergence test supports only up to 
    ! NLSOL_MAXEQUATIONSERROR equations.
    nblocks = min(rb%nblocks,NLSOL_MAXEQUATIONSERROR)
    
    ! Initial test for convergence/divergence.
    call cc_resNormCheck (rproblem,rnonlinearIteration,&
        ite,rx,rb,rd,bconvergence,bdivergence)

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

        ! Call cc_precondDefect to do the preconditioning.
        ! The routine is allowed to change domega during the
        ! iteration if necessary. The nonlinear solver here does not touch
        ! domega anymore, so the callback routine is the only one changing it.
        call cc_precondDefect (rproblem,rnonlinearIteration,&
            ite,rd,rx,rb,domega,bsuccess)
        
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
          call cc_getDefect (rproblem,rnonlinearIteration,ite,rx,rb,rd)

          ! Check the defect for convergence.
          call cc_resNormCheck (rproblem,rnonlinearIteration,&
              ite,rx,rb,rd,bconvergence,bdivergence)

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
  type(t_vectorBlock), intent(IN) :: rrhs
!</input>
  
!<inputoutput>
  ! The problem structure characterising the whole problem.
  type(t_problem), intent(INOUT)                :: rproblem

  ! A nonlinear-iteration structure that configures the core equation to solve.
  ! Can be initialised e.g. by cc_createNonlinearLoop + manual setup of
  ! the coefficients of the terms.
  type(t_ccNonlinearIteration), intent(INOUT) :: rnonlinearIteration

  ! A nonlinear solver configuration.
  ! Can be initialised e.g. by using cc_getNonlinearSolver.
  type(t_nlsolNode), intent(INOUT) :: rnlSolver
  
  ! Initial solution vector. Is replaced by the new solution vector.
  type(t_vectorBlock), intent(INOUT) :: rvector

  ! A temporary block vector for the nonlinear iteration.
  ! OPTIONAL: If not specified, a temporary vector is automatically allocated.
  type(t_vectorBlock), intent(INOUT), target, optional :: rtempBlock
  
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_vectorBlock), pointer :: p_rtempBlock

    if (.not. present(rtempBlock)) then
      ! Create a temporary vector we need for the nonlinear iteration.
      allocate (p_rtempBlock)
      call lsysbl_createVecBlockIndirect (rrhs, p_rtempBlock, .false.)
    else 
      p_rtempBlock => rtempBlock
    end if

    ! Call the nonlinear solver. For preconditioning and defect calculation,
    ! the solver calls our callback routines.
    call cc_nonlinearSolver(rproblem,rnonlinearIteration,rnlSolver,&
        rvector,rrhs,p_rtempBlock)

    if (.not. present(rtempBlock)) then
      ! Release the temporary vector
      call lsysbl_releaseVector (p_rtempBlock)
      deallocate (p_rtempBlock)
    end if
        
  end subroutine

end module
