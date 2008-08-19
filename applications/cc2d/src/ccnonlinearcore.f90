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

MODULE ccnonlinearcore

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
  USE trilinearformevaluation
  USE matrixio
  USE statistics
  USE collection
  USE convection
  
  USE ccmatvecassembly
    
  USE cccallback
  
  IMPLICIT NONE
  
!<constants>

!<constantblock description="Preconditioner identifiers for the defect in the nonlinear iteration">

  ! No preconditioning
  INTEGER, PARAMETER :: CCPREC_NONE         = -1

  ! Preconditioning with inverse mass matrix (not yet implemented)
  INTEGER, PARAMETER :: CCPREC_INVERSEMASS   = 0

  ! Preconditioning by linear solver, solving the linearised system
  INTEGER, PARAMETER :: CCPREC_LINEARSOLVER  = 1

  ! Preconditioning by Newton-Iteration
  INTEGER, PARAMETER :: CCPREC_NEWTON        = 2

  ! Preconditioning by dynamic Newton-Iteration (uses defect correction
  ! and switches automatically to Newton if the error is small enough)
  INTEGER, PARAMETER :: CCPREC_NEWTONDYNAMIC = 3

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
  TYPE t_ccDynamicNewtonControl
  
    ! Minimum number of usul fix point iteration before to switch to
    ! preconfitioning with the Newton matrix. (IFIXMIN)

    INTEGER :: nminFixPointIterations = 0

    ! Maximum number of usul fix point iteration before to switch to
    ! preconfitioning with the Newton matrix. (IFIXMAX)

    INTEGER :: nmaxFixPointIterations = 999

    ! Norm of absolute residuum before applying Newton. 
    ! Newton is only applied
    ! if   ||absolute residuum|| < depsAbsNewton
    ! and  ||relative residuum|| < depsRelNewton.
    ! Otherwise, the usual fix point iteration is used.
    ! Stamndard value = 1E-5.

    REAL(DP) :: depsAbsNewton = 1.0E-5_DP

    ! Norm of relative residuum before applying Newton. 
    ! Newton is only applied
    ! if   ||absolute residuum|| < depsAbsNewton
    ! and  ||relative residuum|| < depsRelNewton.
    ! Otherwise, the usual fix point iteration is used.
    ! Standard value = 1E99 -> The absolute residuum counts.

    REAL(DP) :: depsRelNewton = 1.0E99_DP
  
  END TYPE
  
!</typeblock>


!<typeblock>

  ! Preconditioner structure for CCxD. This structure saves the configuration of the
  ! preconditioner that is used during the nonlinear iteration. 
  
  TYPE t_ccPreconditioner
  
    ! Type of preconditioner.
    ! This is one of the CCPREC_xxxx flags as defined above (CCPREC_INVERSEMASS for
    ! preconditioning with inverse mass matrix, CCPREC_LINEARSOLVER for solving a linear
    ! system, CCPREC_NEWTON for a Newton iteration,...)
    INTEGER :: ctypePreconditioning = CCPREC_NONE
    
    ! Pointer to linear solver node if a linear solver is the preconditioner.
    ! (Thus, this applies for the defect correction and the Newton preconditioner).
    TYPE(t_linsolNode), POINTER :: p_rsolverNode
    
    ! An interlevel projection structure for changing levels
    TYPE(t_interlevelProjectionBlock), POINTER :: p_rprojection
    
    ! Configuration block for the adaptive Newton preconditioner.
    ! Is only valid if ctypePreconditioning=CCPREC_NEWTONDYNAMIC!
    TYPE(t_ccDynamicNewtonControl) :: radaptiveNewton

    ! Temporary scalar vector; used for calculating the nonlinear matrix
    ! on lower levels / projecting the solution from higher to lower levels.
    TYPE(t_vectorScalar), POINTER :: p_rtempVectorSc

    ! Temporary scalar vector; used for calculating the optimal damping
    ! parameter.
    TYPE(t_vectorScalar), POINTER :: p_rtempVectorSc2

  END TYPE

!</typeblock>

!<typeblock>

  ! This type is used to save some situation specific assembly information
  ! during the setup phase of the nonlinear solver. Here it's noted, if
  ! and whose matrices exist and/or must be assmebled transposed to be
  ! compatible with the preconditioner and more. It's more or less
  ! a collection if different flags.
  TYPE t_ccPreconditionerSpecials
  
    ! Whether to use 'adaptive matrices', i.e. set up coarse grid matrices
    ! with the help of fine grid matrices. This is used for very special
    ! discretisations only (e.g. Q1~/Q0). =0: deactivate
    INTEGER :: iadaptiveMatrices    = 0
    
    ! A configuration parameter for adaptive matrices.
    REAL(DP) :: dadMatThreshold     = 0.0_DP

    ! If the preconditioner is a linear solver:
    ! Type of solver.
    ! =0: Gauss elimination (UMFPACK)
    ! =1: Multigrid solver
    INTEGER :: isolverType = 0
    
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
    INTEGER :: ismootherType = 3
    
    ! If the preconditioner is the linear multigrid solver:
    ! Type of coarse grid solver.    
    ! =0: Gauss elimination (UMFPACK)
    ! =1: Defect correction with diagonal VANCA preconditioning.
    ! =2: BiCGStab with diagonal VANCA preconditioning
    INTEGER :: icoarseGridSolverType = 1
        
    ! This flag is set to .TRUE. if there are no Neumann boundary
    ! components. In that case, the pressure matrices of direct
    ! solvers must be changed.
    LOGICAL :: bpressureGloballyIndefinite = .FALSE.
    
  END TYPE

!</typeblock>

!<typeblock>

  ! Represents the core equation on one level of the discretisation.
  ! Collects all information that are necessary to assemble the 
  ! (linearised) system matrix and RHS vector.
  TYPE t_cccoreEquationOneLevel
  
    ! Pointer to the discretisation structure of that level
    TYPE(t_blockDiscretisation), POINTER :: p_rdiscretisation => NULL()

    ! Stokes matrix for that specific level (=nu*Laplace)
    TYPE(t_matrixScalar), POINTER :: p_rmatrixStokes => NULL()

    ! B1-matrix for that specific level. 
    TYPE(t_matrixScalar), POINTER :: p_rmatrixB1 => NULL()

    ! B2-matrix for that specific level. 
    TYPE(t_matrixScalar), POINTER :: p_rmatrixB2 => NULL()

    ! Mass matrix
    TYPE(t_matrixScalar), POINTER :: p_rmatrixMass => NULL()
    
    ! Temporary vector for the interpolation of a solution to a lower level.
    ! Exists only on levels NLMIN..NLMAX-1 !
    TYPE(t_vectorBlock), POINTER :: p_rtempVector => NULL()

    ! Block matrix, which is used in the defect correction / Newton
    ! algorithm as preconditioner matrix of the correspnding underlying
    ! linear sytem. Is usually the (linearise) system matrix or
    ! a Newton matrix. This matrix is changed during the
    ! nonlinear iteration and used e.g. if a linear solver (Multigrid) is
    ! used for preconditioning.
    TYPE(t_matrixBlock), POINTER :: p_rmatrixPreconditioner => NULL()
  
    ! Velocity coupling matrix $A_{12}$.
    ! Exists only of deformation tensor or Newton iteration is used in the
    ! nonlinear iteration.
    TYPE(t_matrixScalar), POINTER :: p_rmatrixVelocityCoupling12 => NULL()

    ! Velocity coupling matrix $A_{21}$.
    ! Exists only of deformation tensor or Newton iteration is used in the
    ! nonlinear iteration.
    TYPE(t_matrixScalar), POINTER :: p_rmatrixVelocityCoupling21 => NULL()
    
    ! Pointer to a B1^T-matrix.
    ! This pointer may point to NULL(). In this case, B1^T is created
    ! by 'virtually transposing' the B1 matrix.
    !
    ! Note: This information is automatically created when the preconditioner
    ! is initialised! The main application does not have to initialise it!
    TYPE(t_matrixScalar), POINTER :: p_rmatrixB1T => NULL()

    ! Pointer to a B2-matrix.
    ! This pointer may point to NULL(). In this case, B2^T is created
    ! by 'virtually transposing' the B2 matrix.
    !
    ! Note: This information is automatically created when the preconditioner
    ! is initialised! The main application does not have to initialise it!
    TYPE(t_matrixScalar), POINTER :: p_rmatrixB2T => NULL()

  END TYPE

!</typeblock>

!<typeblock>

  ! This configuration block configures all parameters that are needed
  ! by the callback routines to perform the nonlinear iteration.
  ! On start of the program, a structure of this type is initialised.
  ! The entries of this structure are saved to the collection structure.
  ! When a callback routine is called, the structure is rebuild from
  ! the collection. When he nonlinear iteration is finished, the
  ! parameters are removed from the colletion again.
  TYPE t_ccNonlinearIteration
  
    ! ALPHA-parameter that controls the weight of the mass matrix in the
    ! core equation. =0.0 for stationary simulations.
    REAL(DP) :: dalpha = 0.0_DP
    
    ! THETA-parameter that controls the weight of the Stokes matrix
    ! in the core equation. =1.0 for stationary simulations.
    REAL(DP) :: dtheta = 0.0_DP
    
    ! GAMMA-parameter that controls the weight in front of the
    ! nonlinearity. =1.0 for Navier-Stokes, =0.0 for Stokes equation.
    REAL(DP) :: dgamma = 0.0_DP

    ! ETA-parameter that switch the B-term on/off.
    REAL(DP) :: deta = 0.0_DP
    
    ! TAU-parameter that switch the B^T-term on/off
    REAL(DP) :: dtau = 0.0_DP
    
    
    ! Minimum allowed damping parameter; OMGMIN
    REAL(DP) :: domegaMin = 0.0_DP
    
    ! Maximum allowed damping parameter; OMGMAX
    REAL(DP) :: domegaMax = 2.0_DP
    
    ! Output level of the nonlinear iteration
    INTEGER :: MT_OutputLevel = 2
    
    ! Minimum discretisation level
    INTEGER :: NLMIN = 1
    
    ! Maximum discretisation level
    INTEGER :: NLMAX = 1
    
    ! Pointer to the solution vector that is changed during the iteration
    TYPE(t_vectorBlock), POINTER :: p_rsolution => NULL()
    
    ! Pointer to the RHS vector
    TYPE(t_vectorBlock), POINTER :: p_rrhs => NULL()
    
    ! A filter chain that is used for implementing boundary conditions into
    ! (linear and nonlinear) defect vectors.
    TYPE(t_filterChain), DIMENSION(:), POINTER :: p_RfilterChain
    
    ! A t_ccPreconditioner saving information about the preconditioner.
    TYPE(t_ccPreconditioner) :: rpreconditioner
    
    ! A t_ccPreconditionerSpecials structure that saves information about
    ! special 'tweaks' in matrices such that everything works.
    TYPE(t_ccPreconditionerSpecials) :: rprecSpecials
    
    ! An array of t_cccoreEquationOneLevel structures for all levels
    ! of the discretisation.
    TYPE(t_cccoreEquationOneLevel), DIMENSION(:), POINTER :: RcoreEquation => NULL()
    
    ! Auxiliary variable: Saves the initial defect in the nonlinear iteration
    REAL(DP), DIMENSION(2) :: DresidualInit = 0.0_DP
    
    ! Auxiliary variable: Saves the last defect in the nonlinear iteration
    REAL(DP), DIMENSION(2) :: DresidualOld = 0.0_DP

    ! Auxiliary variable: Norm of the relative change = norm of the 
    ! preconditioned residual in the nonlinear iteration
    REAL(DP), DIMENSION(3) :: DresidualCorr = 0.0_DP
    
    ! Auxiliary variable: Convergence criteria of the nonlinear solver
    REAL(DP), DIMENSION(5) :: DepsNL = 0.0_DP
    
    ! Auxiliary variable: Last calculated damping parameter
    REAL(DP) :: domegaNL = 0.0_DP
    
    ! Auxiliary variable: Convergence rate of linear solver (if this is
    ! applied as preconditioner).
    REAL(DP) :: drhoLinearSolver = 0.0_DP
    
  END TYPE

!</typeblock>

!</types>

CONTAINS
  
  ! ***************************************************************************
  ! Callback routines for the nonlinear solver
  ! ***************************************************************************

  !<subroutine>
  
    SUBROUTINE cc_getDefect (rproblem,rnonlinearIteration,ite,rx,rb,rd)
  
    USE linearsystemblock
    USE collection
    
    
  !<description>
    ! FOR NONLINEAR ITERATION:
    ! Defect vector calculation callback routine. Based on the current iteration 
    ! vector rx and the right hand side vector rb, this routine has to compute the 
    ! defect vector rd.
  !</description>

  !<input>
    ! Number of current iteration. 0=build initial defect
    INTEGER, INTENT(IN)                           :: ite

    ! Current iteration vector
    TYPE(t_vectorBlock), INTENT(IN),TARGET        :: rx

    ! Right hand side vector of the equation.
    TYPE(t_vectorBlock), INTENT(IN)               :: rb
  !</input>
               
  !<inputoutput>
    ! The problem structure characterising the whole problem.
    TYPE(t_problem), INTENT(INOUT)                :: rproblem

    ! Reference to the nonlinear iteration structure that configures the
    ! main nonlinear equation. Intermediate data is changed during the iteration.
    TYPE(t_ccnonlinearIteration), INTENT(INOUT)   :: rnonlinearIteration

    ! Defect vector b-A(x)x. This must be filled with data by the callback routine.
    TYPE(t_vectorBlock), INTENT(INOUT)            :: rd
  !</inputoutput>
  
  !</subroutine>
      
      ! The nonlinear iteration structure
      TYPE(t_nonlinearCCMatrix) :: rnonlinearCCMatrix
      INTEGER :: ilvmax
      TYPE(t_filterChain), DIMENSION(:), POINTER :: p_RfilterChain
  
      ! a timer structure
      TYPE(t_timer) :: rtimer
                  
      ! DEBUG!!!
      !REAL(DP), DIMENSION(:), POINTER :: p_Ddata,p_Ddata2

      ! Build the nonlinear defect
      CALL lsysbl_copyVector (rb,rd)
      
      ! DEBUG!!!
      !CALL lsysbl_getbase_double (rx,p_Ddata)
      !CALL lsysbl_getbase_double (rd,p_Ddata2)
      
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
      rnonlinearCCMatrix%p_rdiscretisation => &
          rnonlinearIteration%RcoreEquation(ilvmax)%p_rdiscretisation
      rnonlinearCCMatrix%p_rmatrixStokes => &
          rnonlinearIteration%RcoreEquation(ilvmax)%p_rmatrixStokes
      rnonlinearCCMatrix%p_rmatrixB1 => &
          rnonlinearIteration%RcoreEquation(ilvmax)%p_rmatrixB1
      rnonlinearCCMatrix%p_rmatrixB2 => &
          rnonlinearIteration%RcoreEquation(ilvmax)%p_rmatrixB2
      rnonlinearCCMatrix%p_rmatrixB1T => &
          rnonlinearIteration%RcoreEquation(ilvmax)%p_rmatrixB1T
      rnonlinearCCMatrix%p_rmatrixB2T => &
          rnonlinearIteration%RcoreEquation(ilvmax)%p_rmatrixB2T
      rnonlinearCCMatrix%p_rmatrixMass => &
          rnonlinearIteration%RcoreEquation(ilvmax)%p_rmatrixMass

      CALL cc_nonlinearMatMul (rnonlinearCCMatrix,rx,rd,-1.0_DP,1.0_DP)        
      
      p_RfilterChain => rnonlinearIteration%p_RfilterChain
      IF (ASSOCIATED(p_RfilterChain)) THEN    
        CALL filter_applyFilterChainVec (rd, p_RfilterChain)
      END IF
    
      CALL vecfil_discreteNLSlipBCdef (rd)
      
    END SUBROUTINE
    
  ! ***************************************************************************

  !<subroutine>

    SUBROUTINE cc_getOptimalDamping (rproblem,rnonlinearIteration,rd,rx,rb,&
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
    TYPE(t_vectorBlock), INTENT(IN)               :: rx

    ! Current RHS vector of the nonlinear equation
    TYPE(t_vectorBlock), INTENT(IN)               :: rb

    ! Defect vector b-A(x)x. 
    TYPE(t_vectorBlock), INTENT(IN)               :: rd
  !</input>

  !<inputoutput>
    ! The problem structure characterising the whole problem.
    TYPE(t_problem), INTENT(INOUT)                :: rproblem

    ! Reference to the nonlinear iteration structure that configures the
    ! main nonlinear equation. Intermediate data is changed during the iteration.
    TYPE(t_ccnonlinearIteration), INTENT(INOUT)   :: rnonlinearIteration

    ! A temporary vector in the structure of rx
    TYPE(t_vectorBlock), INTENT(INOUT)            :: rtemp1

    ! A 2nd temporary vector in the structure of rx
    TYPE(t_vectorBlock), INTENT(INOUT)            :: rtemp2

    ! Damping parameter. On entry: an initial value given e.g. by the
    ! previous step.
    ! On return: The new damping parameter.
    REAL(DP), INTENT(INOUT)                       :: domega
  !</inputoutput>
  
  !</subroutine>

    ! local variables
    INTEGER :: ilvmax
    REAL(DP) :: dskv1,dskv2
    TYPE(t_matrixBlock) :: rmatrix
    TYPE(t_timer) :: rtimer
    ! A filter chain for the linear solver
    TYPE(t_filterChain), DIMENSION(:), POINTER :: p_RfilterChain
    
    TYPE(t_nonlinearCCMatrix) :: rnonlinearCCMatrix

!    ! DEBUG!!!:
!    real(dp), dimension(:), pointer :: p_vec,p_def,p_temp1,p_temp2,p_da
!    call lsysbl_getbase_double (rd,p_def)
!    call lsysbl_getbase_double (rx,p_vec)
!    call lsysbl_getbase_double (rtemp1,p_temp1)
!    call lsysbl_getbase_double (rtemp2,p_temp2)
!    ilvmax = collct_getvalue_int (p_rcollection,'NLMAX')

      ! Is there anything to do?
      IF (rnonlinearIteration%domegaMin .GE. rnonlinearIteration%domegaMax) THEN
        ! No - cancel.
        domega = rnonlinearIteration%domegaMin
        RETURN
      END IF
      
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

      CALL lsysbl_copyVector(rd,rtemp1)
      CALL lsysbl_vectorLinearComb (rx,rtemp1,1.0_DP,domega)

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
      rnonlinearCCMatrix%p_rmatrixMass => &
          rnonlinearIteration%RcoreEquation(ilvmax)%p_rmatrixMass

      ! Assemble the matrix.        
      CALL cc_assembleMatrix (CCMASM_COMPUTE,CCMASM_MTP_AUTOMATIC,&
          rmatrix,rnonlinearCCMatrix,rtemp1)
      
      ! We don't have to implement any boundary conditions into the matrix
      ! as we apply an appropriate filter to the defect vector after
      ! each matrix-vector-multiplication below!
      
      ! ==================================================================
      ! Second term of the scalar product in the nominator
      ! Calculate the defect rtemp2 = F-T*u_n.
      ! ==================================================================

      CALL lsysbl_copyVector (rb,rtemp2)
      CALL lsysbl_blockMatVec (rmatrix, rx, rtemp2, -1.0_DP, 1.0_DP)
      
      ! This is a defect vector - filter it! This e.g. implements boundary
      ! conditions.
      IF (ASSOCIATED(p_RfilterChain)) THEN
        CALL filter_applyFilterChainVec (rtemp2, p_RfilterChain)
      END IF
      
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

      CALL lsysbl_blockMatVec (rmatrix, rd, rtemp1, 1.0_DP, 0.0_DP)
      
      ! This is a defect vector against 0 - filter it! This e.g. 
      ! implements boundary conditions.
      IF (ASSOCIATED(p_RfilterChain)) THEN
        CALL filter_applyFilterChainVec (rtemp1, p_RfilterChain)
      END IF
      
      ! Release the matrix again
      CALL lsysbl_releaseMatrix (rmatrix)
      
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
      
      IF (dskv2 .LT. 1.0E-40_DP) THEN
        CALL output_line ('dskv2 nearly zero. Optimal damping parameter singular.', &
            OU_CLASS_ERROR,OU_MODE_STD,'cc_getOptimalDamping')
        CALL output_line ('Is the triangulation ok??? .tri-file destroyed?', &
            OU_CLASS_ERROR,OU_MODE_STD,'cc_getOptimalDamping')
        CALL output_line ('Boundary conditions set up properly?', &
            OU_CLASS_ERROR,OU_MODE_STD,'cc_getOptimalDamping')
        CALL sys_halt()
        STOP
      END IF
      
      ! Ok, we have the nominator and the denominator. Divide them
      ! by each other to calculate the new OMEGA.
      
      domega = dskv1 / dskv2
      
      ! And make sure it's in the allowed range:
      
      domega = max(rnonlinearIteration%domegamin, &
                   min(rnonlinearIteration%domegamax,domega))
      
      ! That's it, we have our new Omega.
      
    END SUBROUTINE

  ! ***************************************************************************

  !<subroutine>

    SUBROUTINE cc_precondDefect (rproblem,rnonlinearIteration,&
        ite,rd,rx,rb,domega,bsuccess)
  
    USE linearsystemblock
    USE collection
    
  !<description>
    ! FOR NONLINEAR ITERATION:
    ! Defect vector calculation callback routine. Based on the current iteration 
    ! vector rx and the right hand side vector rb, this routine has to compute the 
    ! defect vector rd. 
  !</description>

  !<inputoutput>
    ! The problem structure characterising the whole problem.
    TYPE(t_problem), INTENT(INOUT)                :: rproblem

    ! Reference to the nonlinear iteration structure that configures the
    ! main nonlinear equation. Intermediate data is changed during the iteration.
    TYPE(t_ccnonlinearIteration), INTENT(INOUT), TARGET   :: rnonlinearIteration

    ! Number of current iteration. 
    INTEGER, INTENT(IN)                           :: ite

    ! Defect vector b-A(x)x. This must be replaced by J^{-1} rd by a preconditioner.
    TYPE(t_vectorBlock), INTENT(INOUT)            :: rd

    ! Damping parameter. Is set to rsolverNode%domega (usually = 1.0_DP)
    ! on the first call to the callback routine.
    ! The callback routine can modify this parameter according to any suitable
    ! algorithm to calculate an 'optimal damping' parameter. The nonlinear loop
    ! will then use this for adding rd to the solution vector:
    ! $$ x_{n+1} = x_n + domega*rd $$
    ! domega will stay at this value until it's changed again.
    REAL(DP), INTENT(INOUT)                       :: domega

    ! If the preconditioning was a success. Is normally automatically set to
    ! TRUE. If there is an error in the preconditioner, this flag can be
    ! set to FALSE. In this case, the nonlinear solver breaks down with
    ! the error flag set to 'preconditioner broke down'.
    LOGICAL, INTENT(INOUT)                        :: bsuccess
  !</inputoutput>
  
  !<input>
    ! Current iteration vector
    TYPE(t_vectorBlock), INTENT(IN), TARGET       :: rx

    ! Current right hand side of the nonlinear system
    TYPE(t_vectorBlock), INTENT(IN), TARGET       :: rb
  !</input>
  
  !</subroutine>
  
    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    TYPE(t_vectorScalar), POINTER :: p_rvectorTemp,p_rvectorTemp2
    TYPE(t_vectorBlock) :: rtemp1,rtemp2
    INTEGER :: ierror
    TYPE(t_linsolNode), POINTER :: p_rsolverNode
    INTEGER, DIMENSION(3) :: Cnorms
    
    INTEGER :: i
    REAL(DP) :: dresInit,dres
    LOGICAL :: bassembleNewton
    TYPE(t_matrixBlock), DIMENSION(:), POINTER :: Rmatrices
    TYPE(t_ccDynamicNewtonControl), POINTER :: p_rnewton
    TYPE(t_filterChain), DIMENSION(:), POINTER :: p_RfilterChain

    ! a local timer object
    TYPE(t_timer) :: rtimer

    ! DEBUG!!!
    !real(dp), dimension(:), pointer :: p_vec,p_def,p_da
    !call lsysbl_getbase_double (rd,p_def)
    !call lsysbl_getbase_double (rx,p_vec)

      SELECT CASE (rnonlinearIteration%rpreconditioner%ctypePreconditioning)
      CASE (CCPREC_NONE)
        ! No preconditioning
        domega = max(rnonlinearIteration%domegamin, &
                    min(rnonlinearIteration%domegamax,domega))
      CASE (CCPREC_LINEARSOLVER,CCPREC_NEWTON,CCPREC_NEWTONDYNAMIC)
        ! Preconditioning with a linear solver.
        !
        ! At first, assemble the preconditioner matrices on all levels
        ! and incorporate all boundary conditions.
        !
        ! Should we assemble the Newton part?
        
        bassembleNewton = .FALSE.
        
        IF (rnonlinearIteration%rpreconditioner%ctypePreconditioning .EQ. &
            CCPREC_NEWTON) THEN
            
          ! Use Newton in any case.
          bassembleNewton = .TRUE.
          
        ELSE IF (rnonlinearIteration%rpreconditioner%ctypePreconditioning .EQ. &
            CCPREC_NEWTONDYNAMIC) THEN
            
          ! Adaptive Newton. Check the iteration and the residual whether to use
          ! Newton or not. But at first, get a shortcut to the parameter structure
          ! of the adaptive Newton...
            
          p_rnewton => rnonlinearIteration%rpreconditioner%radaptiveNewton
          
          IF (ite .GT. p_rnewton%nmaxFixPointIterations) THEN
            ! Force Newton to be used.
            bassembleNewton = .TRUE.
          ELSE
            IF (ite .GT. p_rnewton%nminFixPointIterations) THEN
              ! In this case, the residuum of the last iterate decides on 
              ! whether to use Newton or not.
              dresInit = SQRT(rnonlinearIteration%DresidualInit(1)**2 + &
                            rnonlinearIteration%DresidualInit(2)**2)
              dres    = SQRT(rnonlinearIteration%DresidualOld(1)**2 + &
                            rnonlinearIteration%DresidualOld(2)**2)
              IF ((dres .LT. p_rnewton%depsAbsNewton) .AND. &
                  (dres .LT. p_rnewton%depsRelNewton * dresInit)) THEN
                bassembleNewton = .TRUE.
              END IF
            END IF
            
            ! Otherwise: Use fixpoint iteration...
            
          END IF
          
        END IF
        
        ! Assemble the preconditioner matrices in rnonlinearIteration
        ! on all levels that the solver uses.
        CALL assembleLinsolMatrices (rnonlinearIteration,rproblem%rcollection,&
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
        ALLOCATE(Rmatrices(rnonlinearIteration%NLMIN:rnonlinearIteration%NLMAX))
        DO i=rnonlinearIteration%NLMIN,rnonlinearIteration%NLMAX
          CALL lsysbl_duplicateMatrix ( &
            rnonlinearIteration%RcoreEquation(i)%p_rmatrixPreconditioner, &
            Rmatrices(i), LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
        END DO
        
        CALL linsol_setMatrices(&
            rnonlinearIteration%rpreconditioner%p_rsolverNode,Rmatrices(:))
            
        ! The solver got the matrices; clean up Rmatrices, it was only of temporary
        ! nature...
        DO i=rnonlinearIteration%NLMIN,rnonlinearIteration%NLMAX
          CALL lsysbl_releaseMatrix (Rmatrices(i))
        END DO
        DEALLOCATE(Rmatrices)

        ! Initialise data of the solver. This in fact performs a numeric
        ! factorisation of the matrices in UMFPACK-like solvers.
        CALL linsol_initData (p_rsolverNode, ierror)
        IF (ierror .NE. LINSOL_ERR_NOERROR) STOP
        
        ! Finally solve the system. As we want to solve Ax=b with
        ! b being the real RHS and x being the real solution vector,
        ! we use linsol_solveAdaptively. If b is a defect
        ! RHS and x a defect update to be added to a solution vector,
        ! we would have to use linsol_precondDefect instead.
        CALL linsol_precondDefect (p_rsolverNode,rd)
        
        ! Remember convergence rate for output
        rnonlinearIteration%drhoLinearSolver = p_rsolverNode%dconvergenceRate

        ! Release the numeric factorisation of the matrix.
        ! We don't release the symbolic factorisation, as we can use them
        ! for the next iteration.
        CALL linsol_doneData (p_rsolverNode)
        
        ! Did the preconditioner work?
        bsuccess = p_rsolverNode%iresult .EQ. 0
        
      END SELECT
      
      ! Finally calculate a new damping parameter domega.
      IF (rnonlinearIteration%rpreconditioner%ctypePreconditioning .NE. CCPREC_NONE) THEN
        
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
        CALL lsysbl_createVecFromScalar (p_rvectorTemp,rtemp1)
        CALL lsysbl_enforceStructure (rb,rtemp1)

        CALL lsysbl_createVecFromScalar (p_rvectorTemp2,rtemp2)
        CALL lsysbl_enforceStructure (rb,rtemp2)

        ! Calculate the omega
        CALL cc_getOptimalDamping (rproblem,rnonlinearIteration,&
            rd,rx,rb,rtemp1,rtemp2,domega)

        ! Remember damping parameter for output
        rnonlinearIteration%domegaNL = domega

        ! Release the temp block vectors. This only cleans up the structure.
        ! The data is not released from heap as it belongs to the
        ! scalar temp vectors.
        CALL lsysbl_releaseVector (rtemp2)
        CALL lsysbl_releaseVector (rtemp1)
      END IF
      
      ! Calculate the max-norm of the correction vector.
      ! This is used for the stopping criterium in cc_resNormCheck!
      Cnorms(:) = LINALG_NORMMAX
      rnonlinearIteration%DresidualCorr(:) = lsysbl_vectorNormBlock(rd,Cnorms)
      
      IF ((.NOT. bsuccess) .AND. (domega .GE. 0.001_DP)) THEN
        ! The preconditioner did actually not work, but the solution is not
        ! 'too bad'. So we accept it still.
        bsuccess = .TRUE.
      END IF
      
      IF (bsuccess) THEN
        ! Filter the final defect
        p_RfilterChain => rnonlinearIteration%p_RfilterChain
        CALL filter_applyFilterChainVec (rd, p_RfilterChain)
      END IF

    CONTAINS
    
      SUBROUTINE assembleLinsolMatrices (rnonlinearIteration,rcollection,&
          bassembleNewton,rx,NLMIN,NLMAX)

      USE linearsystemblock
      USE collection

      ! Assembles on every level a matrix for the linear-solver/Newton preconditioner.
      ! bnewton allows to specify whether the Newton matrix or only the standard
      ! system matrix is evaluated. The output is written to the p_rpreconditioner 
      ! matrices specified in the rnonlinearIteration structure.

      ! Reference to a collection structure that contains all parameters of the
      ! discretisation (for nonlinearity, etc.).
      TYPE(t_collection), INTENT(INOUT)                :: rcollection

      ! Nonlinear iteration structure where to write the linearised system matrix to.
      TYPE(t_ccnonlinearIteration), INTENT(INOUT)      :: rnonlinearIteration

      ! TRUE  = Assemble the Newton preconditioner.
      ! FALSE = Assemble the standard defect correction preconditioner
      !         (i.e. the linearised system matrix).
      LOGICAL, INTENT(IN) :: bassembleNewton
      
      ! Minimum level of the preconditioner that is to be initialised.
      INTEGER, INTENT(IN)                              :: NLMIN
      
      ! Maximum level of the preconditioner that is to be initialised.
      ! This must corresponds to the last matrix in Rmatrices.
      INTEGER, INTENT(IN)                              :: NLMAX
      
      ! Current iteration vector. 
      TYPE(t_vectorBlock), INTENT(IN), TARGET          :: rx

      ! local variables
      REAL(DP) :: dnewton
      INTEGER :: ilev,icol
      INTEGER(PREC_MATIDX), DIMENSION(1) :: Irows = (/1/)
      TYPE(t_matrixBlock), POINTER :: p_rmatrix,p_rmatrixFine
      TYPE(t_vectorScalar), POINTER :: p_rvectorTemp
      TYPE(t_vectorBlock), POINTER :: p_rvectorFine,p_rvectorCoarse
      TYPE(t_nonlinearCCMatrix) :: rnonlinearCCMatrix
      TYPE(t_interlevelProjectionBlock), POINTER :: p_rprojection
      TYPE(t_timer) :: rtimer
      ! A filter chain for the linear solver
      TYPE(t_filterChain), DIMENSION(:), POINTER :: p_RfilterChain
      
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
        
        NULLIFY(p_rmatrix)

        DO ilev=NLMAX,NLMIN,-1
        
          ! Get the matrix on the current level.
          ! Shift the previous matrix to the pointer of the fine grid matrix.
          p_rmatrixFine => p_rmatrix
          p_rmatrix => rnonlinearIteration%RcoreEquation(ilev)%p_rmatrixPreconditioner
        
          ! On the highest level, we use rx as solution to build the nonlinear
          ! matrix. On lower levels, we have to create a solution
          ! on that level from a fine-grid solution before we can use
          ! it to build the matrix!
          IF (ilev .EQ. NLMAX) THEN
          
            p_rvectorCoarse => rx
            
          ELSE
            ! We have to discretise a level hierarchy and are on a level < NLMAX.
            
            ! Get the temporary vector on level i. Will receive the solution
            ! vector on that level. 
            p_rvectorCoarse => rnonlinearIteration%RcoreEquation(ilev)%p_rtempVector
            
            ! Get the solution vector on level i+1. This is either the temporary
            ! vector on that level, or the solution vector on the maximum level.
            IF (ilev .LT. NLMAX-1) THEN
              p_rvectorFine => rnonlinearIteration%RcoreEquation(ilev+1)%p_rtempVector
            ELSE
              p_rvectorFine => rx
            END IF

            ! Interpolate the solution from the finer grid to the coarser grid.
            ! The interpolation is configured in the interlevel projection
            ! structure we got from the collection.
            CALL mlprj_performInterpolation (p_rprojection,p_rvectorCoarse, &
                                             p_rvectorFine,p_rvectorTemp)

            ! Apply the filter chain to the temp vector.
            ! This implements the boundary conditions that are attached to it.
            ! NOTE: Deactivated for standard CC2D compatibility -- and because
            ! it has to be checked whether the correct boundary conditions
            ! are attached to that vector!
            ! CALL filter_applyFilterChainVec (p_rvectorCoarse, p_RfilterChain)

          END IF

          ! Initialise the matrix assembly structure rnonlinearCCMatrix to describe the
          ! matrix we want to have.
          rnonlinearCCMatrix%dalpha = rnonlinearIteration%dalpha
          rnonlinearCCMatrix%dtheta = rnonlinearIteration%dtheta
          rnonlinearCCMatrix%dgamma = rnonlinearIteration%dgamma
          IF (bassembleNewton) rnonlinearCCMatrix%dnewton = rnonlinearIteration%dgamma
          rnonlinearCCMatrix%deta = 1.0_DP
          rnonlinearCCMatrix%dtau = 1.0_DP
          rnonlinearCCMatrix%iupwind = rproblem%rstabilisation%iupwind
          rnonlinearCCMatrix%dnu = rproblem%dnu
          rnonlinearCCMatrix%dupsam = rproblem%rstabilisation%dupsam
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
          rnonlinearCCMatrix%p_rmatrixB1T => &
              rnonlinearIteration%RcoreEquation(ilev)%p_rmatrixB1T
          rnonlinearCCMatrix%p_rmatrixB2T => &
              rnonlinearIteration%RcoreEquation(ilev)%p_rmatrixB2T
          rnonlinearCCMatrix%p_rmatrixMass => &
              rnonlinearIteration%RcoreEquation(ilev)%p_rmatrixMass

          ! Assemble the matrix.
          ! If we are on a lower level, we can specify a 'fine-grid' matrix.
          IF (ilev .EQ. NLMAX) THEN
            CALL cc_assembleMatrix (CCMASM_COMPUTE,CCMASM_MTP_AUTOMATIC,&
                p_rmatrix,rnonlinearCCMatrix,p_rvectorCoarse)
          ELSE
            CALL cc_assembleMatrix (CCMASM_COMPUTE,CCMASM_MTP_AUTOMATIC,&
                p_rmatrix,rnonlinearCCMatrix,p_rvectorCoarse,p_rmatrixFine)
          END IF

          ! Boundary conditions
          ! ---------------------------------------------------

          IF (ASSOCIATED(p_RfilterChain)) THEN
            ! Apply the filter chain to the matrix.
            ! As the filter consists only of an implementation filter for
            ! boundary conditions, this implements the boundary conditions
            ! into the system matrix.
            CALL filter_applyFilterChainMat (p_rmatrix, p_RfilterChain)
          ELSE
            ! Call the matrix filter for the boundary conditions to include the BC's
            ! into the matrix.
            CALL matfil_discreteBC (p_rmatrix)
            CALL matfil_discreteFBC (p_rmatrix)
          END IF
            
          ! 'Nonlinear' boundary conditions like slip boundary conditions
          ! are not implemented with a filter chain into a matrix.
          ! Call the appropriate matrix filter of 'nonlinear' boundary
          ! conditions manually:
          CALL matfil_discreteNLSlipBC (p_rmatrix,.TRUE.)
            
        END DO
        
        IF (rnonlinearIteration%rprecSpecials%bpressureGloballyIndefinite) THEN
          
          ! The 3,3-matrix must exist! This is ensured by the initialisation routine.
          !
          ! We have a pure Dirichlet problem. This may give us some difficulties
          ! in the case, the preconditioner uses a direct solver (UMFPACK).
          ! In this case, we have to include a unit vector to the pressure
          ! matrix to make the problem definite!
          IF (rnonlinearIteration%rprecSpecials%isolverType .EQ. 0) THEN
            p_rmatrix => rnonlinearIteration%RcoreEquation(NLMAX)%p_rmatrixPreconditioner
            
            ! Include a unit vector
            CALL mmod_replaceLinesByZero(p_rmatrix%RmatrixBlock(3,1),Irows)
            CALL mmod_replaceLinesByZero(p_rmatrix%RmatrixBlock(3,2),Irows)
            CALL mmod_replaceLinesByUnit(p_rmatrix%RmatrixBlock(3,3),Irows)
            
          END IF
          
          IF (rnonlinearIteration%rprecSpecials%isolverType .EQ. 1) THEN
          
            ! If we have a MG solver, We also check the coarse grid solver for 
            ! the same thing!
            ! What we don't check is the smoother, thus we assume that smoothers
            ! are always solvers that allow the applicance of a filter chain.
            IF (rnonlinearIteration%rprecSpecials%icoarseGridSolverType .EQ. 0) THEN
              p_rmatrix => rnonlinearIteration%RcoreEquation(NLMIN)%p_rmatrixPreconditioner
              
              ! Include a unit vector
              CALL mmod_replaceLinesByZero(p_rmatrix%RmatrixBlock(3,1),Irows)
              CALL mmod_replaceLinesByZero(p_rmatrix%RmatrixBlock(3,2),Irows)
              CALL mmod_replaceLinesByUnit(p_rmatrix%RmatrixBlock(3,3),Irows)
              
            END IF
            
          END IF
            
        ELSE
          
          ! The 3,3-block must be a zero-matrix. So if it's present, clear it.
          IF (lsysbl_isSubmatrixPresent(p_rmatrix,3,3)) &
            CALL lsyssc_clearMatrix (p_rmatrix%RmatrixBlock(3,3))
          
        END IF        

      END SUBROUTINE
      
    END SUBROUTINE

  ! ***************************************************************************

    SUBROUTINE cc_resNormCheck (rproblem,rnonlinearIteration,&
        ite,rx,rb,rd,bconvergence,bdivergence)
  
    USE linearsystemblock
    USE collection
    
  !<description>
    ! Residual norm calculation & printing routine.
    ! This routine is called each time the norm of the residuum was calculated.
    ! It has to check the current residuum for convergence and/or divergence
    ! and can print the residuum to screen.
  !</description>

  !<inputoutput>
    ! The problem structure characterising the whole problem.
    TYPE(t_problem), INTENT(INOUT)                :: rproblem

    ! Reference to the nonlinear iteration structure that configures the
    ! main nonlinear equation. Intermediate data is changed during the iteration.
    TYPE(t_ccnonlinearIteration), INTENT(INOUT)   :: rnonlinearIteration

    ! Number of current iteration. Is set to 0 when the callback routine
    ! is called the first time. In this situation, rd describes the initial
    ! defect vector.
    INTEGER, INTENT(IN)                           :: ite
  
    ! Current iteration vector
    TYPE(t_vectorBlock), INTENT(IN), TARGET       :: rx

    ! Right hand side vector of the equation.
    TYPE(t_vectorBlock), INTENT(IN), TARGET       :: rb

    ! Defect vector b-A(x)x.
    TYPE(t_vectorBlock), INTENT(IN), TARGET       :: rd
  !</inputoutput>
  
  !<output>
    ! Must be set to TRUE by the callback routine if the residuum rd
    ! is within a desired tolerance, so that the solver should treat
    ! the iteration as 'converged'.
    LOGICAL, INTENT(OUT)                        :: bconvergence
  
    ! Must be set to TRUE by the callback routine if the residuum rd
    ! is out of a desired tolerance, so that the solver should treat
    ! the iteration as 'diverged'.
    LOGICAL, INTENT(OUT)                        :: bdivergence
  !</output>

      ! local variables
      REAL(DP), DIMENSION(3) :: Dresiduals
      REAL(DP) :: dresOld,drhoNL,ddelP,ddelU,dtmp,dresU,dresDIV,dres,dresINIT
      REAL(DP) :: depsD,depsDiv,depsUR,depsPR,depsRES
      INTEGER, DIMENSION(3) :: Cnorms

      ! Calculate norms of the solution/defect vector
      CALL cc_getDefectNorm (rx,rb,rd,Dresiduals)
      Dresiduals(3) = SQRT(Dresiduals(1)**2 + Dresiduals(2)**2)

      ! In the first iteration (initial defect), print the norm of the defect
      ! and save the norm of the initial residuum to the structure
      IF (ite .EQ. 0) THEN
      
        CALL output_separator (OU_SEP_MINUS)     
        CALL output_line (' IT  RELU     RELP     DEF-U    DEF-DIV'// &
                          '  DEF-TOT  RHONL    OMEGNL   RHOMG')
        CALL output_separator (OU_SEP_MINUS)     
        CALL output_line ('  0                   '// &
            TRIM(sys_sdEP(Dresiduals(1),9,2))//&
            TRIM(sys_sdEP(Dresiduals(2),9,2))//&
            TRIM(sys_sdEP(Dresiduals(3),9,2)))
        CALL output_separator (OU_SEP_MINUS)     

        rnonlinearIteration%DresidualInit (1:2) = Dresiduals(1:2)
        rnonlinearIteration%DresidualOld (1:2) = Dresiduals(1:2)

        bconvergence = .FALSE.
        bdivergence = .FALSE.
      
      ELSE
        ! In the other iterations, calculate the relative change and
        ! test the convergence criteria.
      
        ! Old defect:
        dresOld = SQRT(rnonlinearIteration%DresidualOld(1)**2 + &
                       rnonlinearIteration%DresidualOld(2)**2)
        
        ! Replace the 'old' residual by the current one
        rnonlinearIteration%DresidualOld(1:2) = Dresiduals(1:2)

        ! Nonlinear convergence rate
        drhoNL = (Dresiduals(3)/dresOld) ** (1.0_DP/REAL(ite,DP))
        
        ! Calculate norms of the solution/defect vector, calculated above
        dresU   = Dresiduals(1)
        dresDIV = Dresiduals(2)
        dres    = SQRT(dresU**2 + dresDIV**2)
        
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
        Dresiduals(:) = lsysbl_vectorNormBlock (rx,Cnorms)

        dtmp = MAX(Dresiduals(1),Dresiduals(2))
        IF (dtmp .LT. 1.0E-8_DP) dtmp = 1.0_DP
        ddelU = MAX(rnonlinearIteration%DresidualCorr(1),&
                    rnonlinearIteration%DresidualCorr(2))/dtmp
        
        dtmp = Dresiduals(3)
        IF (dtmp .LT. 1.0E-8_DP) dtmp = 1.0_DP
        ddelP = rnonlinearIteration%DresidualCorr(3)/dtmp
        
        ! Check if the nonlinear iteration can prematurely terminate.
        !        
        ! Get the stopping criteria from the parameters.
        ! Use the DepsNL data according to the initialisation above.
        dresINIT = SQRT(rnonlinearIteration%DresidualInit(1)**2 + &
                        rnonlinearIteration%DresidualInit(2)**2)
                        
        ! dresInit=0 may hardly occur -- except when we expect 'no flow'.
        ! But to prevent a check against "something<=0" in this case below,
        ! set dresInit to something <> 0.
        IF (dresINIT .EQ. 0.0_DP) dresINIT = 1.0_DP

        depsD   = rnonlinearIteration%DepsNL(1)
        depsDiv = rnonlinearIteration%DepsNL(2)
        depsUR  = rnonlinearIteration%DepsNL(3)
        depsPR  = rnonlinearIteration%DepsNL(4)
        depsRES = rnonlinearIteration%DepsNL(5)*dresINIT
        
        ! All residual information calculated.
        ! Check for divergence; use a 'NOT' for better NaN handling.
        bdivergence = .NOT. (dres/dresOld .LT. 1E2)
        
        ! Check for convergence
        IF((ddelU .LE. depsUR).AND.(ddelP .LE. depsPR)   .AND. &
           (dresU .LE. depsD) .AND.(dresDiv .LE. depsDiv).AND. &
           (dres .LE. depsRES)) THEN
          bconvergence = .TRUE.
        ELSE
          bconvergence = .FALSE.
        END IF

        IF ((ddelU .LT. SYS_EPSREAL*1E2_DP) .AND. &
            (ddelP .LT. SYS_EPSREAL*1E2_DP)) THEN
          ! We are hard on machine exactness, so stop the iteraton
          bconvergence =.TRUE.
        END IF
        
        ! Print residual information
        CALL output_line ( &
            TRIM(sys_si(ite,3))//' '//&
            TRIM(sys_sdEP(ddelU,9,2))// &
            TRIM(sys_sdEP(ddelP,9,2))// &
            TRIM(sys_sdEP(dresU,9,2))// &
            TRIM(sys_sdEP(dresDIV,9,2))// &
            TRIM(sys_sdEP(dres,9,2))// &
            TRIM(sys_sdEP(drhoNL,9,2))// &
            TRIM(sys_sdEP(rnonlinearIteration%domegaNL,9,2))// &
            TRIM(sys_sdEP(rnonlinearIteration%drhoLinearSolver,9,2)) &
            )
        
      END IF
      
    END SUBROUTINE

! ***************************************************************************

!<subroutine>

  SUBROUTINE cc_getDefectNorm (rvector,rrhs,rdefect,Dresiduals)
  
!<description>
  ! Calculates a couple of norms from a given solution, defect and RHS vector
  ! of the nonlinear system. This can be used to check convergence criteria
  ! etc.
!</description>

!<input>
  ! The solution vector which is modified later during the nonlinear iteration.
  TYPE(t_vectorBlock), INTENT(IN) :: rvector

  ! The right-hand-side vector to use in the equation
  TYPE(t_vectorBlock), INTENT(IN) :: rrhs
  
  ! A defect vector calculated with rvector and rrhs
  TYPE(t_vectorBlock), INTENT(IN) :: rdefect
!</input>

!<output>
  ! An array receiving different defect norms calculated by the above vectors.
  ! Dresiduals(1) = RESU   = ||defect_u|| / ||rhs|| = velocity residual
  ! Dresiduals(2) = RESDIV = ||p|| / ||u||          = divergence residual
  REAL(DP), DIMENSION(:), INTENT(OUT) :: Dresiduals
!</output>

!</subroutine>

    ! local variables
    REAL(DP) :: dresF,DresTmp(2),dnormU
    INTEGER, DIMENSION(2) :: Cnorms

    !-----------------------------------------------------------------------
    !     Compute the relative l2-norms  RESU,RESDIV
    !-----------------------------------------------------------------------

    Cnorms(:) = LINALG_NORMEUCLID

    ! RESF := max ( ||F1||_E , ||F2||_E )

    DresTmp = lsysbl_vectorNormBlock (rrhs,Cnorms)
    dresF = MAX(DresTmp(1),DresTmp(2))
    IF (dresF .LT. 1.0E-8_DP) dresF = 1.0_DP

    !               || (D1,D2) ||_E
    ! RESU = -----------------------------
    !        max ( ||F1||_E , ||F2||_E )

    DresTmp = lsysbl_vectorNormBlock (rdefect,Cnorms)
    Dresiduals(1) = SQRT(DresTmp(1)**2+DresTmp(2)**2)/dresF

    ! DNORMU = || (U1,U2) ||_l2 

    DresTmp = lsysbl_vectorNormBlock (rvector,Cnorms)
    dnormU = SQRT(DresTmp(1)**2+DresTmp(2)**2)
    IF (dnormU .LT. 1.0E-8_DP) dnormU = 1.0_DP

    !             || DP ||_E
    ! RESDIV = ----------------
    !          || (U1,U2) ||_E

    Dresiduals(2) = &
        lsyssc_vectorNorm (rdefect%RvectorBlock(3),LINALG_NORMEUCLID) / dnormU
        
  END SUBROUTINE

! ***************************************************************************
! Routines to solve the core equation
! ***************************************************************************

!<subroutine>

  SUBROUTINE cc_getNonlinearSolver (rnlSolver, rparamList, sname)
  
!<description>
  ! Creates a nonlinear solver node rnlSolver and initialises it with parameters
  ! from the INI/DAT files given in the parameter list rparamList.
  ! sname is the name of a section in rparamList that configures the
  ! parameter of the nonlinear solver.
!</description>

!<input>
  ! Parameter list that contains the parameters from the INI/DAT file(s).
  TYPE(t_parlist), INTENT(IN) :: rparamList
  
  ! Name of the section in the parameter list containing the parameters
  ! of the nonlinear solver.
  CHARACTER(LEN=*), INTENT(IN) :: sname
!</input>

!<output>
  ! A t_nlsolNode structure that contains the configuration of the nonlinear
  ! solver. The parameters are initialised according to the information
  ! in the section sname of the parameter list rparamList
  TYPE(t_nlsolNode) :: rnlSolver
!</output>

!</subroutine>

    ! local variables
    TYPE(t_parlstSection), POINTER :: p_rsection

    ! Check that there is a section called sname - otherwise we
    ! cannot create anything!
    
    CALL parlst_querysection(rparamList, sname, p_rsection) 

    IF (.NOT. ASSOCIATED(p_rsection)) THEN
      CALL output_line ('Cannot create nonlinear solver; no section '''//&
          TRIM(sname)//'''!', &
          OU_CLASS_ERROR,OU_MODE_STD,'cc_getNonlinearSolver')
      CALL sys_halt()
    END IF
    
    ! Parse the given parameters now to initialise the solver node.
    ! This is now CCxD-specific!
    
    CALL parlst_getvalue_int (p_rsection, 'nminIterations', &
                              rnlSolver%nminIterations, rnlSolver%nminIterations)

    CALL parlst_getvalue_int (p_rsection, 'nmaxIterations', &
                              rnlSolver%nmaxIterations, rnlSolver%nmaxIterations)

    CALL parlst_getvalue_int (p_rsection, 'ioutputLevel', &
                              rnlSolver%ioutputLevel, rnlSolver%ioutputLevel)

    CALL parlst_getvalue_double (p_rsection, 'depsUR', &
                                 rnlSolver%DepsRel(1), rnlSolver%DepsRel(1))
    rnlSolver%DepsRel(2) = rnlSolver%DepsRel(1)

    CALL parlst_getvalue_double (p_rsection, 'depsPR', &
                                 rnlSolver%DepsRel(3), rnlSolver%DepsRel(3))

    CALL parlst_getvalue_double (p_rsection, 'depsD', &
                                 rnlSolver%DepsAbs(1), rnlSolver%DepsAbs(1))
    rnlSolver%DepsAbs(2) = rnlSolver%DepsAbs(1)

    CALL parlst_getvalue_double (p_rsection, 'depsDiv', &
                                 rnlSolver%DepsAbs(3), rnlSolver%DepsAbs(3))

    ! Initial damping parameter.
    CALL parlst_getvalue_double (p_rsection, 'domegaIni', &
                                 rnlSolver%domega, rnlSolver%domega)

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE cc_nonlinearSolver (rproblem,rnonlinearIteration,&
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
  TYPE(t_problem), INTENT(INOUT)                :: rproblem

  ! Reference to the nonlinear iteration structure that configures the
  ! main nonlinear equation. Intermediate data is changed during the iteration.
  TYPE(t_ccnonlinearIteration), INTENT(INOUT)   :: rnonlinearIteration

  ! The nonlinear solver node that configures the solution process.
  TYPE(t_nlsolNode), INTENT(INOUT)              :: rsolverNode
  
  ! INPUT: Initial solution vector.
  ! OUTPUT: Final iteration vector.
  TYPE(t_vectorBlock), INTENT(INOUT)            :: rx
             
  ! Temporary vector. Must be of the same size/type as rx/rb.
  TYPE(t_vectorBlock), INTENT(INOUT)            :: rd
!</inputoutput>
  
!<input>
  ! Right hand side vector of the equation.
  TYPE(t_vectorBlock), INTENT(IN)               :: rb
!</input>

!</subroutine>

  ! local variables
  INTEGER :: ite
  REAL(DP), DIMENSION(NLSOL_MAXEQUATIONSERROR) :: DvecNorm
  TYPE(t_vectorBlock) :: rtemp
  REAL(DP) :: domega
  INTEGER :: nblocks
  LOGICAL :: bconvergence,bdivergence,bsuccess
  
    ! In case our preconditioner is a matrix-vector multiplication,
    ! allocate memory for another temporary vector used 
    ! during the MV multiplication.
    IF (rsolverNode%cpreconditioner .EQ. NLSOL_PREC_MATRIX) THEN
      CALL lsysbl_createVecBlockIndirect (rx,rtemp,.FALSE.)
    END IF
    
    ! Status reset
    rsolverNode%iresult = 0
    
    ! Calculate the initial nonlinear defect to rd:   d = b-A(x)x
    CALL cc_getDefect (rproblem,rnonlinearIteration,0,rx,rb,rd)
    
    ite = 0
    rsolverNode%icurrentIteration = ite

    ! The standard convergence/divergence test supports only up to 
    ! NLSOL_MAXEQUATIONSERROR equations.
    nblocks = MIN(rb%nblocks,NLSOL_MAXEQUATIONSERROR)
    
    ! Initial test for convergence/divergence.
    CALL cc_resNormCheck (rproblem,rnonlinearIteration,&
        ite,rx,rb,rd,bconvergence,bdivergence)

    ! Perform at least nminIterations iterations
    IF (ite .LT. rsolverNode%nminIterations) bconvergence = .FALSE.
  
    ! Check for divergence
    IF (bdivergence) rsolverNode%iresult = 1
      
    IF ((.NOT. bconvergence) .AND. (.NOT. bdivergence)) THEN
    
      ! Let's do the nonlinear loop...
      !
      ! Initialise the domega-value for the damping of the correction
      ! as prescribed by the parameters of the solver.
      ! The callback routines might change it...
      domega = rsolverNode%domega
      
      ! Perform at most nmaxIterations iterations
      DO ite = 1,rsolverNode%nmaxIterations
        rsolverNode%icurrentIteration = ite
      
        ! Perform preconditioning on the defect vector:  u = J^{-1} d
        bsuccess = .TRUE.

        ! Call cc_precondDefect to do the preconditioning.
        ! The routine is allowed to change domega during the
        ! iteration if necessary. The nonlinear solver here does not touch
        ! domega anymore, so the callback routine is the only one changing it.
        CALL cc_precondDefect (rproblem,rnonlinearIteration,&
            ite,rd,rx,rb,domega,bsuccess)
        
        ! If bsuccess=false, the preconditioner had an error.
        IF (.NOT. bsuccess) THEN
          CALL output_line ('NLSOL: Iteration '//&
              TRIM(sys_siL(ite,10))//' canceled as the preconditioner went down!')
          rsolverNode%iresult = 3
          EXIT
        END IF
        
        ! If domega=0.0, the solution vector would stay unchanged. In this
        ! case, the nonlinear solver would not proceed at all, and the next
        ! iteration would behave exactly as before!
        ! So in this case, there's nothing to do, we can stop the iteration.
        IF (domega .EQ. 0.0_DP) THEN
          CALL output_line ('NLSOL: Iteration '//&
              TRIM(sys_siL(ite,10))//' canceled as there is no progress anymore!')
          EXIT
        ELSE
          ! Add the correction vector in rd to rx;
          ! damped by the damping parameter:           x := x + domega u
          CALL lsysbl_vectorLinearComb (rd,rx,domega,1.0_DP)
          
          ! Calculate the new nonlinear defect to rd:  d = b-A(x)x
          CALL cc_getDefect (rproblem,rnonlinearIteration,ite,rx,rb,rd)

          ! Check the defect for convergence.
          CALL cc_resNormCheck (rproblem,rnonlinearIteration,&
              ite,rx,rb,rd,bconvergence,bdivergence)

          ! Perform at least nminIterations iterations
          IF (ite .LT. rsolverNode%nminIterations) bconvergence = .FALSE.
        
          ! Check for convergence
          IF (bconvergence) EXIT
        
          ! Check for divergence
          IF (bdivergence) THEN
            rsolverNode%iresult = 1
            EXIT
          END IF
        
        END IF
        
      END DO ! ite
      
    END IF ! not converged and not diverged

    ! Set ITE to NIT to prevent printing of "NIT+1" of the loop was
    ! completed

    IF (ite .GT. rsolverNode%nmaxIterations) &
      ite = rsolverNode%nmaxIterations
      
    IF (.NOT. bconvergence) THEN 
      ! Convergence criterion not reached, but solution did not diverge.
      rsolverNode%iresult = -1
    END IF

    rsolverNode%iiterations = ite
    
    ! Release temporary memory
    IF (rsolverNode%cpreconditioner .EQ. NLSOL_PREC_MATRIX) THEN
      CALL lsysbl_releaseVector (rtemp)
    END IF
  
    ! Nonlinear loop finished.

  END SUBROUTINE

! ***************************************************************************

!<subroutine>

  SUBROUTINE cc_solveCoreEquation (rproblem,rnonlinearIteration,rnlSolver,&
      rvector,rrhs,rtempBlock)
  
!<description>
  ! This routine invokes the nonlinear solver to solve the core equation
  ! as configured in the core equation structure.
!</description>

!<input>
  ! The right-hand-side vector to use in the equation.
  TYPE(t_vectorBlock), INTENT(IN) :: rrhs
!</input>
  
!<inputoutput>
  ! The problem structure characterising the whole problem.
  TYPE(t_problem), INTENT(INOUT)                :: rproblem

  ! A nonlinear-iteration structure that configures the core equation to solve.
  ! Can be initialised e.g. by cc_createNonlinearLoop + manual setup of
  ! the coefficients of the terms.
  TYPE(t_ccNonlinearIteration), INTENT(INOUT) :: rnonlinearIteration

  ! A nonlinear solver configuration.
  ! Can be initialised e.g. by using cc_getNonlinearSolver.
  TYPE(t_nlsolNode), INTENT(INOUT) :: rnlSolver
  
  ! Initial solution vector. Is replaced by the new solution vector.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rvector

  ! A temporary block vector for the nonlinear iteration.
  ! OPTIONAL: If not specified, a temporary vector is automatically allocated.
  TYPE(t_vectorBlock), INTENT(INOUT), TARGET, OPTIONAL :: rtempBlock
  
!</inputoutput>

!</subroutine>

    ! local variables
    TYPE(t_vectorBlock), POINTER :: p_rtempBlock

    IF (.NOT. PRESENT(rtempBlock)) THEN
      ! Create a temporary vector we need for the nonlinear iteration.
      ALLOCATE (p_rtempBlock)
      CALL lsysbl_createVecBlockIndirect (rrhs, p_rtempBlock, .FALSE.)
    ELSE 
      p_rtempBlock => rtempBlock
    END IF

    ! Call the nonlinear solver. For preconditioning and defect calculation,
    ! the solver calls our callback routines.
    CALL cc_nonlinearSolver(rproblem,rnonlinearIteration,rnlSolver,&
        rvector,rrhs,p_rtempBlock)

    IF (.NOT. PRESENT(rtempBlock)) THEN
      ! Release the temporary vector
      CALL lsysbl_releaseVector (p_rtempBlock)
      DEALLOCATE (p_rtempBlock)
    END IF
        
  END SUBROUTINE

END MODULE
