!##############################################################################
!# ****************************************************************************
!# <name> ACnonlinearcore </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains the routines to solve the stationary core equation
!# of the problem with a nonlinear solver.
!# (This is comparable to the NSDEF routine in older FEAT versions)
!#
!# The discretised core equation reads at the moment:
!#
!#  $$        A_1 y    = f_1 $$
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
!#  1.) AC_getNonlinearSolver
!#      -> Initialise nonlinear solver configuration with information
!#         from INI/DAT files
!#
!#  2) AC_solveCoreEquation
!#      -> Starts the nonlinear iteration to solve the core equation.
!#
!# Callback routines for the nonlinear solver:
!#
!#  6.) AC_getDefect
!#      -> Callback routine. Calculate nonlinear defect
!#
!#  7.) AC_getOptimalDamping
!#     -> Auxiliary routine. Calculate optimal damping parameter
!#
!#  8.) AC_precondDefect
!#      -> Callback routine. Preconditioning of nonlinear defect
!#
!#  9.) AC_resNormCheck
!#      -> Callback routine: Check the residuals for convergence
!#
!# 10.) AC_nonlinearSolver
!#      -> The actual nonlinear solver; similar to NSDEF in old AC versions.
!#
!# To solve a system with the core equation, one has to deal with two
!# main structures. On one hand, one has a nonlinear iteration structure
!# t_nlsolNode from the kernel; this is initialised by AC_getNonlinearSolver.
!# On the other hand, one has to maintain a 'nonlinear iteration structure'
!# of type t_ACNonlinearIteration, which configures the core equation and
!# specifies parameters for the solver how to work.
!#
!# The basic procedure for dealing with the core equation is as follows:
!#
!#  a) AC_getNonlinearSolver  -> Basic initialisation on the nonlinear solver
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
!#  e) AC_solveCoreEquation -> Solve the core equation with the nonlinear
!#                               solver structure from AC_getNonlinearSolver.
!#
!#  f) Release core equation structure
!#
!# Stebs b), c), d) and f) are done in the routines in the file
!# ACnonlinearcoreinit.f90. This file here contains only the very core
!# routines of the nonlinear iteration, leaving all initialisation to
!# ACnonlinearcoreinit.f90.
!#
!# </purpose>
!##############################################################################

module AllenCahn_nonlinearcore

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
  use trilinearformevaluation
  use matrixio
  use statistics
  use collection
  use convection
  use linearsystemblock
  
  use AllenCahn_matvec
  use AllenCahn_callback
  
  implicit none
  
!<constants>

!<constantblock description="Preconditioner identifiers for the defect in the nonlinear iteration">

  ! No preconditioning
  integer, parameter :: ACPREC_NONE         = -1

  ! Preconditioning with inverse mass matrix (not yet implemented)
  integer, parameter :: ACPREC_INVERSEMASS   = 0

  ! Preconditioning by linear solver, solving the linearised system
  integer, parameter :: ACPREC_LINEARSOLVER  = 1

  ! Preconditioning by Newton-Iteration
  integer, parameter :: ACPREC_NEWTON        = 2

!</constantblock>

!</constants>
! We do not need CoreEquationOnelevel and t_ACpreconditioner

!<typeblock>
  ! This configuration block configures all parameters that are needed
  ! by the callback routines to perform the nonlinear iteration.
  ! On start of the program, a structure of this type is initialised.
  ! The entries of this structure are saved to the collection structure.
  ! When a callback routine is called, the structure is rebuild from
  ! the collection. When he nonlinear iteration is finished, the
  ! parameters are removed from the colletion again.
  type t_ACNonlinearIteration
  
    ! ALPHA-parameter that controls the weight of the mass matrix in the
    ! core equation. =0.0 for stationary simulations.
    real(DP) :: dalpha = 0.0_DP
    
    ! THETA-parameter that controls the weight of the Laplace term
    real(DP) :: dtheta = 0.0_DP

	! CONV-parameter that controls the wright of the Convective term
    real(DP) :: dconv = 0.0_DP

    ! GAMMA-parameter that controls the weight in front of the
    ! nonlinearity.
    real(DP) :: dgamma = 0.0_DP
  
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
    
	! Pointer to linear solver node if a linear solver is the preconditioner.
    ! (Thus, this applies for the defect correction and the Newton preconditioner).
    type(t_linsolNode), pointer :: p_rsolverNode => NULL()

    ! Pointer to the solution vector that is changed during the iteration
    type(t_vectorBlock), pointer :: p_rsolution => null()
    
    ! Pointer to the RHS vector
    type(t_vectorBlock), pointer :: p_rrhs => null()
    
    ! A filter chain that is used for implementing boundary conditions into
    ! (linear and nonlinear) defect vectors.
    type(t_filterChain), dimension(:), pointer :: p_RfilterChain
 
    ! Auxiliary variable: Saves the initial defect in the nonlinear iteration
    real(DP), dimension(1) :: DresidualInit = 0.0_DP
    
    ! Auxiliary variable: Saves the last defect in the nonlinear iteration
    real(DP), dimension(1) :: DresidualOld = 0.0_DP

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

    ! Temporary scalar vector; used for calculating the nonlinear matrix
    ! on lower levels / projecting the solution from higher to lower levels.
    type(t_vectorScalar), pointer :: p_rtempVectorSc => NULL()

    ! Temporary scalar vector; used for calculating the optimal damping
    ! parameter.
    type(t_vectorScalar), pointer :: p_rtempVectorSc2 => NULL()
    
  end type

!</typeblock>
!</types>

contains
  ! ***************************************************************************
!<subroutine>

  subroutine AC_nonlinearMatMul (rACproblem,rx,rACnonlinearIteration, rd,dcx,dcd)

!<description>
  ! This routine performs a matrix-vector multiplication with a nonlinear
  ! matrix:
  !return      rd := dcx A(ry) rx + dcd rd     dcd=1.0_DP
  ! with the system matrix A(.) defined by the configuration in rACproblem.
  ! The caller must initialise the rACproblem according to how the
  ! matrix should look like.
  !
  ! The parameter ry is optional. If specified, this parameter defines where to
  ! evaluate the nonlinearity (if the system matrix $A$ contains a nonlinearity).
  ! If not specified, ry=rx is assumed.
  !
  ! The routine will not include any boundary conditions in the defect.
!</description>

  ! A problem structure providing all necessary 'source' information
  ! about how to set up the matrix.
  type(t_ACproblem), intent(INOUT) :: rACproblem

  ! This vector specifies the 'x' that is multiplied to the matrix.
  type(t_vectorBlock), intent(IN), target :: rx

  type(t_ACnonlinearIteration), intent(INOUT) :: rACnonlinearIteration

  ! Multiplication factor in front of the term 'A(ry) rx'.
  real(DP), intent(IN) :: dcx

  ! Multiplication factor in front of the term 'rd'.
  real(DP), intent(IN) :: dcd

!</input>

!<inputoutput>
  ! Destination vector. cx*A(ry)*rx is subtracted from this vector.
  type(t_vectorBlock), intent(INOUT) :: rd

!</inputoutput>
  
!</subroutine>

    ! local variables
    type(t_matrixBlock) :: rmatrix
    
!     ! DEBUG!!!
!     real(dp), dimension(:), pointer :: p_DdataX,p_DdataD
!
!     call lsysbl_getbase_double (rx,p_DdataX)
!     call lsysbl_getbase_double (rd,p_DdataD)

    ! Probably weight the input vector.
    if (dcd .ne. 1.0_DP) then
      call lsysbl_scaleVector (rd,dcd)
    end if

    ! The system matrix looks like:
    ! A = A = \alpha M + \theta Laplace + \conv Conv + \gamma Nonlinear
	
    !MCai~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!------------------Start to cal dcx * rmatrix * rx---------------------
	! (1,1) block
    ! Subtract the mass matrix stuff?
    if (rACnonlinearIteration%dalpha .ne. 0.0_DP) then
      call lsyssc_scalarMatVec (&
	  rACproblem%RlevelInfo(rACproblem%NLMAX)%rmatrixMass%rmatrixBlock(1,1), &
          rx%RvectorBlock(1), rd%RvectorBlock(1), &
          dcx*rACnonlinearIteration%dalpha, 1.0_DP)
    end if

    if (rACnonlinearIteration%dtheta .ne. 0.0_DP) then
   	  call lsyssc_scalarMatVec (&
	  rACproblem%RlevelInfo(rACproblem%NLMAX)%rmatrixStatic%rmatrixBlock(1,1),&
          rx%RvectorBlock(1), rd%RvectorBlock(1), &
          dcx*rACnonlinearIteration%dtheta, 1.0_DP)
    end if

    if (rACnonlinearIteration%dconv .ne. 0.0_DP) then
      call lsyssc_scalarMatVec (&
	  rACproblem%RlevelInfo(rACproblem%NLMAX)%rmatrixConv%rmatrixBlock(1,1), &
          rx%RvectorBlock(1), rd%RvectorBlock(1), &
          dcx*rACnonlinearIteration%dconv, 1.0_DP)
    end if

    ! Nonlinear term
    if (rACnonlinearIteration%dgamma .ne. 0.0_DP) then
      !MCai, we first remark this to make sure nonlinear solver is correct.
!     ! nonlinear term: rd=rd + (f(c), \ki)/eps**2
      call nonlinear_calcDefect (rACproblem, rx%RvectorBlock(1), &
	       rd%RvectorBlock(1), dcx*rACnonlinearIteration%dgamma, 1.0_DP)
    end if

  end subroutine

  ! ***************************************************************************
!<subroutine>
!MCai
  subroutine nonlinear_calcDefect (rACproblem, rvectorScalar, &
		         rdefect, dcx, dcd)

!<description>
  ! Calculates:
  !      rd=dcd * rd + dcx * N(rx)
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_ACproblem), intent(INOUT) :: rACproblem

  ! A pointer to the RHS vector.
  type(t_vectorScalar), intent(IN)  :: rvectorScalar
  type(t_vectorScalar), intent(INOUT) :: rdefect
   
  real(DP), intent(in) :: dcx
  real(DP), intent(in) :: dcd
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!</inputoutput>

!</subroutine>

  ! local variables
    INTEGER :: i
  
    ! A linear form describing the analytic problem to solve
    type(t_linearForm) :: rlinform

    ! A pointer to the discretisation structure with the data.
    type(t_blockDiscretisation), pointer :: p_rdiscretisation
	! A tmp vector created from rvectorScalar
    type(t_vectorBlock), target :: rvectorBlocktmp
    type(t_vectorScalar) :: rdefecttmp
    ! A tmp collection for passing data
    type(t_collection) :: rcollectionTmp

!~~~Mcai,
!   Notice that there is a nonlinear term f(\phi)
    !
    ! At first set up the corresponding linear form (f,\Phi_j):

    ! Create two temp vectors.
     call lsysbl_createVecFromScalar (rvectorScalar,rvectorBlocktmp)
     call lsyssc_duplicateVector (rdefect,rdefecttmp,&
	                       LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)

    rlinform%itermCount = 1
    rlinform%Idescriptors(1) = DER_FUNC

    ! get the discretisation
    p_rdiscretisation =>rACproblem%RlevelInfo(rACproblem%NLMAX)%p_rdiscretisation

    call collct_init(rcollectionTmp)

    ! We pass our collection structure as well to this routine,
    ! so the callback routine has access to everything what is

    rcollectionTmp%p_rvectorQuickAccess1 => rvectorBlocktmp

    ! We use rcollection to pass rvectorBlocktmp into nonlinear coeff.
     call linf_buildVectorScalar (p_rdiscretisation%RspatialDiscr(1),rlinform,.false.,&
                rdefecttmp,coeff_nonlinear,&
                rcollectionTmp)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! rd=dcd * rd + dcx * rdefecttmp
    call lsyssc_vectorLinearComb (rdefecttmp, &
	          rdefect,dcx,dcd)
    
	! Release tem vectors.
    call lsysbl_releaseVector(rvectorBlocktmp)
    call lsyssc_releaseVector(rdefecttmp)

    call collct_done(rcollectionTmp)
  end subroutine

! ***************************************************************************
  subroutine AC_updatePreconditioner (rACproblem,rACnonlinearIteration,&
                              rACvector)
                   
!<description>
  ! This routine has to be called whenever the system matrices change.
  ! It initialises (depending on the system matrices) the matrices of the
  ! preconditioner or performs an update of them.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_ACproblem), intent(INOUT) :: rACproblem
!</inputoutput>

!<input>
  ! The current solution vector.
  type(t_vectorBlock), intent(INOUT) :: rACvector
!   ! The right-hand-side vector to use in the equation
!   type(t_vectorBlock), intent(IN) :: rrhs


!<inputoutput>
  ! Nonlinar iteration structure saving data for the callback routines.
  ! Preconditioner data is saved here.
  type(t_ACnonlinearIteration), intent(INOUT) :: rACnonlinearIteration

!</inputoutput>
!</subroutine>

    ! local variables
    integer :: NLMIN,NLMAX
    integer :: i
!MCai, we do not need cmatrixType so far, but if we use different preconditioners,
!We may need it.

    ! Error indicator during initialisation of the solver
    integer :: ierror
  
    ! A pointer to the matrix of the preconditioner
    type(t_matrixBlock), pointer :: p_rmatrixPreconditioner

   ! A pointer to the discretisation structure with the data.
    type(t_blockDiscretisation), POINTER :: p_rdiscretisation

    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    type(t_matrixBlock), dimension(:), pointer :: Rmatrices

  !~~~~~~~~~~~~~~~~Begin of updating preconditioner~~~~~~~~~~~~~~~~~~~~~~~~~~~
      NLMIN = rACproblem%NLMIN
      NLMAX = rACproblem%NLMAX

      allocate(rACnonlinearIteration%p_rtempVectorSc)
      call lsyssc_createVector (rACnonlinearIteration%p_rtempVectorSc,&
                                rACvector%NEQ,.false.)
      
      ! Set up a second temporary vector that we need for calculating
      ! the optimal defect correction.
      allocate(rACnonlinearIteration%p_rtempVectorSc2)
      call lsyssc_createVector (rACnonlinearIteration%p_rtempVectorSc2,&
                                rACvector%NEQ,.false.,ST_DOUBLE)

      call AC_assembleMatJacobian(rACproblem, rACvector)
      ! In this subroutine, we directly linear combine various matrices.
      ! we use rACproblem%rmatrix, we do not need p_rmatrixPreconditioner
!~~~Here we need to assemble matrix for preconditioner, and save them in rACproblem~~~~~~~~~~
      do i=NLMIN,NLMAX
         ! First clear system matrix
        call lsysbl_clearMatrix(rACproblem%RlevelInfo(i)%rmatrix)
        if (rACnonlinearIteration%dalpha .ne. 0.0_DP) then
          !Mass matrix
          !call lsyssc_duplicateMatrix (rACproblem%RlevelInfo(i)%rmatrixMass%RmatrixBlock(1,1),&
          !               rACproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(1,1),&
          !               LSYSSC_DUP_IGNORE,LSYSSC_DUP_COPY)
          call lsyssc_matrixLinearComb (rACproblem%RlevelInfo(i)%rmatrixMass%RmatrixBlock(1,1),&
                         rACnonlinearIteration%dalpha,&
                         rACproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(1,1),&
                         1.0_DP,&
                         rACproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(1,1),&
                         .FALSE.,.FALSE.,.TRUE.,.TRUE.)
        end if

        if (rACnonlinearIteration%dtheta .ne. 0.0_DP) then
		  ! Laplace Matrix
          call lsyssc_matrixLinearComb (rACproblem%RlevelInfo(i)%rmatrixStatic%RmatrixBlock(1,1),&
                         rACnonlinearIteration%dtheta,&
                         rACproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(1,1),&
                         1.0_DP,&
                         rACproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(1,1),&
                         .FALSE.,.FALSE.,.TRUE.,.TRUE.)
        end if
 
        if (rACnonlinearIteration%dconv .ne. 0.0_DP) then
          call lsyssc_matrixLinearComb (rACproblem%RlevelInfo(i)%rmatrixConv%RmatrixBlock(1,1),&
                         rACnonlinearIteration%dconv,&
                         rACproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(1,1),&
                         1.0_DP,&
                         rACproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(1,1),&
                         .FALSE.,.FALSE.,.TRUE.,.TRUE.)
        end if

        ! Nonlinear term
        if (rACnonlinearIteration%dgamma .ne. 0.0_DP) then
          !MCai, we first remark this to make sure nonlinear solver is correct.
          !     ! nonlinear term: rd=rd - (f(c), \ki)/eps**2
          ! Jacobian  Matrix, take care, we may have many choices of preconditioners: Jacobian or...
          ! Here, the system matrix is time independent.
          call lsyssc_matrixLinearComb (rACproblem%RlevelInfo(i)%rmatrixJacobian%RmatrixBlock(1,1),&
                         rACnonlinearIteration%dgamma,&
                         rACproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(1,1),&
                         1.0_DP,&
                         rACproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(1,1),&
                         .FALSE.,.FALSE.,.TRUE.,.TRUE.)
        end if
          ! Attach boundary conditions ??? How to implement it in this ?
      end do

      ! Attach the system matrices to the solver.
      !
      ! For this purpose, copy the matrix structures from the preconditioner
      ! matrices to Rmatrix.
      allocate(Rmatrices(1:NLMAX))
!MCai, May 2, we use the following
      do i=NLMIN,NLMAX
        call lsysbl_duplicateMatrix ( &
               rACproblem%rlevelInfo(i)%rmatrix, &
               Rmatrices(i), LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
      end do

      rACnonlinearIteration%p_rsolverNode => rACproblem%p_rsolverNode
      call linsol_setMatrices(&
            rACnonlinearIteration%p_rsolverNode,Rmatrices(NLMIN:NLMAX))

      ! The solver got the matrices; clean up Rmatrices, it was only of temporary
      ! nature...
      do i=NLMIN,NLMAX
        call lsysbl_releaseMatrix (Rmatrices(i))
      end do
      deallocate(Rmatrices)

      ! Initialise structure/data of the solver. This allows the
      ! solver to allocate memory / perform some precalculation
      ! to the problem.
      call linsol_initStructure (rACnonlinearIteration%p_rsolverNode,&
              ierror)

      if (ierror .ne. LINSOL_ERR_NOERROR) stop

  end subroutine
! ***************************************************************************
! Routines to solve the core equation
! ***************************************************************************
! ***************************************************************************

!<subroutine>

  subroutine AC_solveCoreEquation (rACproblem,rnonlinearIteration,rnlSolver,&
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
  type(t_ACproblem), intent(INOUT)                :: rACproblem

  ! A nonlinear-iteration structure that configures the core equation to solve.
  ! Can be initialised e.g. by AC_createNonlinearLoop + manual setup of
  ! the coefficients of the terms.
  type(t_ACNonlinearIteration), intent(INOUT) :: rnonlinearIteration

  ! A nonlinear solver configuration.
  ! Can be initialised e.g. by using AC_getNonlinearSolver.
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
    call AC_nonlinearSolver(rACproblem,rnonlinearIteration,rnlSolver,&
        rvector,rrhs,p_rtempBlock)

    ! Gather statistics
    call stat_stopTimer(rtimer)

    if (.not. present(rtempBlock)) then
      ! Release the temporary vector
      call lsysbl_releaseVector (p_rtempBlock)
      deallocate (p_rtempBlock)
    end if
        
  end subroutine
!<subroutine>

  subroutine AC_nonlinearSolver (rACproblem,rnonlinearIteration,&
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
  !     For this purpose, the routine AC_getDefect is called.
  !
  ! 2.) Preconditioning of the nonlinear defect
  !         $$  u_n  :=  J^{-1} d_n  $$
  !     with a preconditioner $J^{-1}$. The preconditioner is
  !     realised by the routine AC_precondDefect.
  !
  ! 3.) Correction
  !         $$  x_{n+1}  :=  x_n  +  \omega u_n  $$
  !     $\omega$ is usually = 1. The defined routine AC_precondDefect
  !     can modify $\omega$ in every iteration.
  !
  ! The callback routine AC_resNormCheck is called in every iteration.
  ! This routine can check the residuum for convergence/divergence and can print
  ! information about the norm of the residuum to screen.
  !
  ! At the beginning of the nonlinear iteration, the routines
  ! AC_getResiduum and AC_resNormCheck are called once with ite=0 to calculate
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
  type(t_ACproblem), intent(INOUT)                :: rACproblem

  ! Reference to the nonlinear iteration structure that configures the
  ! main nonlinear equation. Intermediate data is changed during the iteration.
  type(t_ACnonlinearIteration), intent(INOUT)   :: rnonlinearIteration

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
    
 !  print *, 'before getDefect'
    ! Calculate the initial nonlinear defect to rd:   d = b-A(x)x
    call AC_getDefect (rACproblem,rnonlinearIteration,0,rx,rb,rd)
 !   print *, 'after getdefect'
 !   stop
    ite = 0
    rsolverNode%icurrentIteration = ite

    ! The standard convergence/divergence test supports only up to
    ! NLSOL_MAXEQUATIONSERROR equations.
    nblocks = min(rb%nblocks,NLSOL_MAXEQUATIONSERROR)
    
    ! Initial test for convergence/divergence.
    call AC_resNormCheck (rACproblem,rnonlinearIteration,&
        ite,rx,rb,rd,bconvergence,bdivergence)


    ! Get the initial residuum; AC_resNormCheck saved that to DresidualOld.
    rsolverNode%DinitialDefect(1) = abs(rnonlinearIteration%DresidualInit(1))
    rsolverNode%DfinalDefect(1) = rsolverNode%DinitialDefect(1)
  
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

        ! Call AC_precondDefect to do the preconditioning.
        ! The routine is allowed to change domega during the
        ! iteration if necessary. The nonlinear solver here does not touch
        ! domega anymore, so the callback routine is the only one changing it.
        call AC_precondDefect (rACproblem,rnonlinearIteration,&
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
          call AC_getDefect (rACproblem,rnonlinearIteration,ite,rx,rb,rd)

          ! Check the defect for convergence.
          call AC_resNormCheck (rACproblem,rnonlinearIteration,&
              ite,rx,rb,rd,bconvergence,bdivergence)
              
          ! Get the new residual; AC_resNormCheck saved that to DresidualOld.
          rsolverNode%DfinalDefect(1) = abs(rnonlinearIteration%DresidualOld(1))

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
!******************************************************************************
  !<subroutine>

    subroutine AC_getOptimalDamping (rACproblem,rnonlinearIteration,rd,rx,rb,&
        rtemp1,rtemp2,domega)
  !<description>
    ! This subroutine is called inside of the nonlinear loop, to be precise,
    ! inside of AC_precondDefect. It calculates an optiman damping parameter
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
    type(t_vectorBlock), intent(INOUT)               :: rx

    ! Current RHS vector of the nonlinear equation
    type(t_vectorBlock), intent(IN)               :: rb

    ! Defect vector b-A(x)x.
    type(t_vectorBlock), intent(IN)               :: rd
  !</input>

  !<inputoutput>
    ! The problem structure characterising the whole problem.
    type(t_ACproblem), intent(INOUT)                :: rACproblem

    ! Reference to the nonlinear iteration structure that configures the
    ! main nonlinear equation. Intermediate data is changed during the iteration.
    type(t_ACnonlinearIteration), intent(INOUT)   :: rnonlinearIteration

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
      ! (u1)  := (u1)         + \omega C^{-1}( (f1)   [ A        ] (u1) )
      !                                      |--------------------------|
      !                                              = d_n
      !                               |---------------------------------|
      !                                        = Y = (y1,y2,yp) = rd
      !
      !    C = T(u_n)
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
      ! Initialise the matrix assembly structure rACnonlinearIteration to describe the
      ! matrix we want to have.

! MCai, this is different from cc2d or CH2d.
      ! Assemble the matrix: to get rmatrix
      call AC_updatePreconditioner(rACproblem,rnonlinearIteration, rx)
      call lsysbl_copyMatrix(rACproblem%RlevelInfo(rACproblem%NLMAX)%rmatrix, rmatrix)

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
            OU_CLASS_ERROR,OU_MODE_STD,'AC_getOptimalDamping')
        call output_line ('Is the triangulation ok??? .tri-file destroyed?', &
            OU_CLASS_ERROR,OU_MODE_STD,'AC_getOptimalDamping')
        call output_line ('Boundary conditions set up properly?', &
            OU_CLASS_ERROR,OU_MODE_STD,'AC_getOptimalDamping')
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
    end subroutine
  ! ***************************************************************************


    subroutine AC_resNormCheck (rACproblem,rnonlinearIteration,&
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
    type(t_ACproblem), intent(INOUT)                :: rACproblem

    ! Reference to the nonlinear iteration structure that configures the
    ! main nonlinear equation. Intermediate data is changed during the iteration.
    type(t_ACnonlinearIteration), intent(INOUT)   :: rnonlinearIteration

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
      real(DP), dimension(1) :: Dresiduals
      real(DP) :: dresOld,drhoNL,ddelU,dtmp,dresU,dresINIT
      real(DP) :: depsD,depsDiv,depsUR, depsRES
      integer, dimension(1) :: Cnorms

      ! Calculate norms of the solution/defect vector
      call AC_getDefectNorm (rx,rb,rd,Dresiduals)
      
      ! In the first iteration (initial defect), print the norm of the defect
      ! and save the norm of the initial residuum to the structure
      if (ite .eq. 0) then
      
        call output_separator (OU_SEP_MINUS)
        call output_line (' IT  RELU      DEF-U   '// &
                          '   RHONL    OMEGNL   RHOMG')
        call output_separator (OU_SEP_MINUS)
        call output_line ('  0         '// &
            trim(sys_sdEP(Dresiduals(1),9,2)))
        call output_separator (OU_SEP_MINUS)

        rnonlinearIteration%DresidualInit (1) = Dresiduals(1)
        rnonlinearIteration%DresidualOld  (1) = Dresiduals(1)

        bconvergence = .false.
        bdivergence = .false.
      
      else
        ! In the other iterations, calculate the relative change and
        ! test the convergence criteria.
      
        ! Old defect:
        dresOld = abs(rnonlinearIteration%DresidualOld(1))
        
        ! Replace the 'old' residual by the current one
        rnonlinearIteration%DresidualOld(1) = Dresiduals(1)

        ! Nonlinear convergence rate
        drhoNL = (Dresiduals(1)/dresOld) ** (1.0_DP/real(ite,DP))
        
        ! Calculate norms of the solution/defect vector, calculated above
        dresU   = Dresiduals(1)
     
        ! Calculate relative maximum changes
        ! This simply calculates some postprocessing values of the relative
        ! change in the solution.
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

        dtmp = abs(Dresiduals(1))
        if (dtmp .lt. 1.0E-8_DP) dtmp = 1.0_DP
        ddelU = abs(rnonlinearIteration%DresidualCorr(1))/dtmp
        
        dtmp = Dresiduals(1)
       
        ! Check if the nonlinear iteration can prematurely terminate.
        !
        ! Get the stopping criteria from the parameters.
        ! Use the DepsNL data aACording to the initialisation above.
        dresINIT = abs(rnonlinearIteration%DresidualInit(1))
                        
        ! dresInit=0 may hardly occur -- except when we expect 'no flow'.
        ! But to prevent a check against "something<=0" in this case below,
        ! set dresInit to something <> 0.
        if (dresINIT .eq. 0.0_DP) dresINIT = 1.0_DP

        depsD   = rnonlinearIteration%DepsNL(1)
        depsUR  = rnonlinearIteration%DepsNL(3)
        depsRES = rnonlinearIteration%DepsNL(5)*dresINIT
        
        ! All residual information calculated.
        ! Check for divergence; use a 'NOT' for better NaN handling.
        bdivergence = .not. (dresU/dresOld .lt. 1E2)
        
        ! Check for convergence
        if((ddelU .le. depsUR).and. &
           (dresU .le. depsD)) then
          bconvergence = .true.
        else
          bconvergence = .false.
        end if

        if ((ddelU .lt. SYS_EPSREAL_DP*1E2_DP)) then
          ! We are hard on machine exactness, so stop the iteraton
          bconvergence =.true.
        end if
        
        ! Print residual information
        call output_line ( &
            trim(sys_si(ite,3))//' '//&
            trim(sys_sdEP(ddelU,9,2))// &
            trim(sys_sdEP(dresU,9,2))// &
            trim(sys_sdEP(drhoNL,9,2))// &
            trim(sys_sdEP(rnonlinearIteration%domegaNL,9,2))// &
            trim(sys_sdEP(rnonlinearIteration%drhoLinearSolver,9,2)) &
            )
        
      end if
      
    end subroutine

! ***************************************************************************

!<subroutine>

  subroutine AC_getDefectNorm (rvector,rrhs,rdefect,Dresiduals)
  
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
  real(DP), dimension(:), intent(OUT) :: Dresiduals
!</output>

!</subroutine>

    ! local variables
    real(DP) :: dresF,DresTmp(1),dnormU
    integer, dimension(1) :: Cnorms

    !-----------------------------------------------------------------------
    !     Compute the relative l2-norms  RESU,RESDIV
    !-----------------------------------------------------------------------

    Cnorms(:) = LINALG_NORMEUCLID

    ! RESF := abs ( ||F1||_E)

    call lsysbl_vectorNormBlock (rrhs,Cnorms,DresTmp)
    dresF = abs(DresTmp(1))
    
    if (dresF .lt. 1.0E-8_DP) dresF = 1.0_DP

    ! RESU = || D1 ||

    call lsysbl_vectorNormBlock (rdefect,Cnorms,DresTmp)
    Dresiduals(1) = abs(DresTmp(1))/dresF

    ! DNORMU = || U1 ||_l2

    call lsysbl_vectorNormBlock (rvector,Cnorms,DresTmp)
    dnormU = abs(DresTmp(1))
    if (dnormU .lt. 1.0E-8_DP) dnormU = 1.0_DP

  end subroutine

  ! ***************************************************************************
  ! Callback routines for the nonlinear solver
  ! ***************************************************************************

  !<subroutine>
  
    subroutine AC_getDefect (rACproblem,rnonlinearIteration,ite,rx,rb,rd)
    ! return rd= rb - A(rx)rx
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
    type(t_ACproblem), intent(INOUT)                :: rACproblem

    ! Reference to the nonlinear iteration structure that configures the
    ! main nonlinear equation. Intermediate data is changed during the iteration.
    type(t_ACnonlinearIteration), intent(INOUT)   :: rnonlinearIteration

    ! Defect vector b-A(x)x. This must be filled with data by the callback routine.
    type(t_vectorBlock), intent(INOUT)            :: rd
  !</inputoutput>
  
  !</subroutine>
      
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
          
      ! returns  rd := dcx A(ry) rx + dcd rd
      call AC_nonlinearMatMul (rACproblem,rx, &
	                    rnonlinearIteration,rd,-1.0_DP,1.0_DP)
       	                     
      
      p_RfilterChain => rnonlinearIteration%p_RfilterChain
      if (associated(p_RfilterChain)) then
        call filter_applyFilterChainVec (rd, p_RfilterChain)
      end if
      
      ! Gather statistics
      call stat_stopTimer(rtimer)
      
    end subroutine

  !*****************************************************************************
  ! Compute the nonlinear matrix times a vector, it is used to form defect
  !*****************************************************************************
  !MCai, we first remark this part.
  !MCai, May 2, we deleted the following code (similar to cc2d assembleVelDefect)
  !shall we implement stabilization?
  !    subroutine assembleACDefect (rACproblem,&
  !        rmatrix,rvector,rdefect,rvelocityVector,dvectorWeight)
        
  ! With a matrix 'A' of the theoretical form
    !
    !       A := dalpha*M + dtheta*Laplace + dgamma*N(p_rvector) +
    !            dnewton*N*(p_rvector)
    !
  ! ***************************************************************************

  !<subroutine>

    subroutine AC_precondDefect (rACproblem,rnonlinearIteration,&
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
    type(t_ACproblem), intent(INOUT)                :: rACproblem

    ! Reference to the nonlinear iteration structure that configures the
    ! main nonlinear equation. Intermediate data is changed during the iteration.
    type(t_ACnonlinearIteration), intent(INOUT), target   :: rnonlinearIteration

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
    type(t_vectorBlock), intent(INOUT), target       :: rx

    ! Current right hand side of the nonlinear system
    type(t_vectorBlock), intent(IN), target       :: rb
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
    type(t_filterChain), dimension(:), pointer :: p_RfilterChain

    ! a local timer object
    type(t_timer) :: rtimer

    ! DEBUG!!!
    !real(dp), dimension(:), pointer :: p_vec,p_def,p_da
    !call lsysbl_getbase_double (rd,p_def)
    !call lsysbl_getbase_double (rx,p_vec)

      ! Preconditioning with a linear solver.
      !
      ! At first, assemble the preconditioner matrices on all levels
      ! and incorporate all boundary conditions.
      !

      bassembleNewton = .false.

      dresInit = abs(rnonlinearIteration%DresidualInit(1))
      dres    = abs(rnonlinearIteration%DresidualOld(1))

      call stat_clearTimer(rtimer)
      call stat_startTimer(rtimer)

      ! Assemble the preconditioner matrices in rnonlinearIteration
      ! on all levels that the solver uses.
      call assembleLinsolMatrices (rACproblem, rnonlinearIteration,rACproblem%rcollection,&
            bassembleNewton,rx,rnonlinearIteration%NLMIN,rnonlinearIteration%NLMAX)

      ! Gather statistics
      call stat_stopTimer(rtimer)
      call stat_clearTimer(rtimer)
      call stat_startTimer(rtimer)

      ! Our 'parent' (the caller of the nonlinear solver) has prepared
      ! a preconditioner node for us (a linear solver with symbolically
      ! factorised matrices). Get this from the collection.

      p_rsolverNode => rnonlinearIteration%p_rsolverNode

      ! Re-attach the system matrices to the solver.
      ! Note that no pointers and no handles are changed, so we can savely do
      ! that without calling linsol_doneStructure/linsol_doneStructure.
      ! This simply informs the solver about possible new scaling factors
      ! in the matrices in case they have changed...\

      allocate(Rmatrices(rnonlinearIteration%NLMIN:rnonlinearIteration%NLMAX))
      do i=rnonlinearIteration%NLMIN,rnonlinearIteration%NLMAX
         call lsysbl_duplicateMatrix ( &
           rACproblem%RlevelInfo(i)%rmatrix, &
           Rmatrices(i), LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
      end do

      call linsol_setMatrices(&
           rnonlinearIteration%p_rsolverNode,Rmatrices(:))

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

      ! Gather statistics
      call stat_stopTimer(rtimer)
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
      ! Remember convergence rate for output
      rnonlinearIteration%drhoLinearSolver = p_rsolverNode%dconvergenceRate
       
      ! Release the numeric factorisation of the matrix.
      ! We don't release the symbolic factorisation, as we can use them
      ! for the next iteration.
      call linsol_doneData (p_rsolverNode)
      ! Did the preconditioner work?
      bsuccess = p_rsolverNode%iresult .eq. 0
           
!MCai, we need assembleLinsolMatrices because, we assemble Jacobian in terms of
! previous iteration, rather than previous time step solution.

      ! If there is preconditioner for nonlinear system, we need to cal optimal omega
        p_rvectorTemp => rnonlinearIteration%p_rtempVectorSc
        p_rvectorTemp2 => rnonlinearIteration%p_rtempVectorSc2
      ! Finally calculate a new damping parameter domega.
!      if (rnonlinearIteration%rpreconditioner%ctypePreconditioning .ne. ACPREC_NONE) then
      ! Both temp vectors are scalar, but we need block-vectors in the
      ! structure of rx/rb. Derive block vectors in that structure that
      ! share their memory with the scalar temp vectors. Note that the
      ! temp vectors are created large enough by our parent!
      call lsysbl_createVecFromScalar (p_rvectorTemp,rtemp1)
      call lsysbl_enforceStructure (rb,rtemp1)
      call lsysbl_createVecFromScalar (p_rvectorTemp2,rtemp2)
      call lsysbl_enforceStructure (rb,rtemp2)

      ! Calculate the omega
! MCai, ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      call AC_getOptimalDamping (rACproblem,rnonlinearIteration,&
         rd,rx,rb,rtemp1,rtemp2,domega)
!	  domega =1.0_DP
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! Remember damping parameter for output
      rnonlinearIteration%domegaNL = domega

      ! Release the temp block vectors. This only cleans up the structure.
      ! The data is not released from heap as it belongs to the
      ! scalar temp vectors.
      call lsysbl_releaseVector (rtemp2)
      call lsysbl_releaseVector (rtemp1)
!      end if
      
      ! Calculate the max-norm of the correction vector.
      ! This is used for the stopping criterium in AC_resNormCheck!
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

      subroutine assembleLinsolMatrices (rACproblem, rnonlinearIteration,rcollection,&
          bassembleNewton,rx,NLMIN,NLMAX)

      use linearsystemblock
      use collection

      ! Assembles on every level a matrix for the linear-solver/Newton preconditioner.
      ! bnewton allows to specify whether the Newton matrix or only the standard
      ! system matrix is evaluated. The output is written to the p_rpreconditioner
      ! matrices specified in the rnonlinearIteration structure.

      ! The problem structure characterising the whole problem.
      type(t_ACproblem), intent(INOUT)                :: rACproblem

      ! Reference to a collection structure that contains all parameters of the
      ! discretisation (for nonlinearity, etc.).
      type(t_collection), intent(INOUT)                :: rcollection

      ! Nonlinear iteration structure where to write the linearised system matrix to.
      type(t_ACnonlinearIteration), intent(INOUT)      :: rnonlinearIteration

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
      type(t_vectorBlock), intent(INOUT), target          :: rx

      ! local variables
      real(DP) :: dnewton
      integer :: ilev,icol
      integer(PREC_MATIDX), dimension(1) :: Irows = (/1/)
      type(t_matrixBlock), pointer :: p_rmatrix,p_rmatrixFine
      type(t_vectorScalar), pointer :: p_rvectorTemp
      type(t_vectorBlock), pointer :: p_rvectorFine,p_rvectorCoarse
      type(t_interlevelProjectionBlock), pointer :: p_rprojection
      type(t_timer) :: rtimer
      ! A filter chain for the linear solver
      type(t_filterChain), dimension(:), pointer :: p_RfilterChain
      
      ! DEBUG!!!
    !    real(dp), dimension(:), pointer :: p_vec,p_def,p_da
    !    call lsysbl_getbase_double (rd,p_def)
    !    call lsysbl_getbase_double (rx,p_vec)

        ! Get the filter chain. We need tghat later to filter the matrices.
        p_RfilterChain => rnonlinearIteration%p_RfilterChain

        ! On all levels, we have to set up the nonlinear system matrix,
        ! so that the linear solver can be applied to it.
! no rrhs
       call AC_updatePreconditioner (rACproblem,rnonlinearIteration, rx)

      end subroutine

    end subroutine

  !<subroutine>

  subroutine AC_getNonlinearSolver (rnlSolver)
  
!<description>
  ! MCai modified the code so that it does not depending on input DAT file
  ! Creates a nonlinear solver node rnlSolver and initialises it with parameters
  ! from the INI/DAT files given in the parameter list rparamList.
  ! sname is the name of a section in rparamList that configures the
  ! parameter of the nonlinear solver.
!</description>

  ! A t_nlsolNode structure that contains the configuration of the nonlinear
  ! solver. The parameters are initialised aACording to the information
  ! in the section sname of the parameter list rparamList
  type(t_nlsolNode) :: rnlSolver
!</output>

!</subroutine>
    
    ! Parse the given parameters now to initialise the solver node.
    ! This is now ACxD-specific!
    
    rnlSolver%nminIterations=1
    rnlSolver%nmaxIterations=20
    rnlSolver%ioutputLevel=2
    rnlSolver%DepsRel(1)=1.0E-5
    rnlSolver%DepsAbs(1)=1.0E-5
    rnlSolver%domega=1.0

  end subroutine

  ! ***************************************************************************

!<subroutine>

    subroutine AC_initNonlinearLoop (rACproblem,nlmin,nlmax,rvector,rrhs,&
        rnonlinearIteration,sname)
  
  !<description>
    ! Initialises the given nonlinear iteration structure rnonlinearIteration.
    ! Creates the structure with AC_createNonlinearLoop and saves all
    ! problem dependent parameters and matrices in it.
    ! The routine automatically calls AC_createNonlinearLoop to initialise
    ! the structure with internal parameters.
    ! Note: This does not initialise the preconditioner in that structure!
  !</description>
  
  !<input>
    ! A problem structure saving problem-dependent information.
    type(t_ACproblem), intent(INOUT), target :: rACproblem

    ! Minimum refinement level in the rACproblem structure that is allowed to be used
    ! by the preconditioners.
    integer, intent(IN) :: nlmin
  
    ! Maximum refinement level in the rACproblem structure that is allowed to be used
    ! by the preconditioners. This level must correspond to rvector and rrhs.
    integer, intent(IN) :: nlmax

    ! The solution vector which is modified later during the nonlinear iteration.
    type(t_vectorBlock), intent(IN), target :: rvector

    ! The right-hand-side vector to use in the equation
    type(t_vectorBlock), intent(IN), target :: rrhs

    ! Name of the section in the parameter list containing the parameters
    ! of the nonlinear solver.
    character(LEN=*), intent(IN), optional :: sname
  !</input>

  !<output>
    ! Nonlinar iteration structure saving data for the callback routines.
    ! Is filled with data.
    type(t_ACnonlinearIteration), intent(OUT) :: rnonlinearIteration
  !</output>

  !</subroutine>

      ! local variables
      integer :: ilevel
      logical :: bneumann
      type(t_parlstSection), pointer :: p_rsection

      rnonlinearIteration%NLMIN = NLMIN
      rnonlinearIteration%NLMAX = NLMAX

!      rnonlinearIteration%MT_OutputLevel = rACproblem%MT_OutputLevel
      rnonlinearIteration%NLMAX = rACproblem%NLMAX

      rnonlinearIteration%domegaMin=1.0
      rnonlinearIteration%domegaMax=1.0
      ! Save pointers to the RHS and solution vector
      rnonlinearIteration%p_rsolution => rvector
      rnonlinearIteration%p_rrhs => rrhs
    
      ! Set the preconditioner to 'nothing'

      ! MCai: we do not need to assign the matrix pointers in the nonlinear iteration

      ! Clear auxiliary variables for the nonlinear iteration
      rnonlinearIteration%DresidualInit = 0.0_DP
      rnonlinearIteration%DresidualOld  = 0.0_DP

      ! Get stopping criteria of the nonlinear iteration
      rnonlinearIteration%DepsNL(1)=1.0E-5
      rnonlinearIteration%DepsNL(3)=1.0E-5
      rnonlinearIteration%DepsNL(5)=1.0E-5

      ! Set up a filter that modifies the block vectors/matrix
      ! according to boundary conditions. This filter chain is applied to each
      ! defect vector during the linear and nonlinear iteration.
      allocate(rnonlinearIteration%p_RfilterChain(3))
    
      ! Initialise the first filter of the filter chain as boundary
      ! implementation filter for defect vectors:
      rnonlinearIteration%p_RfilterChain(1)%ifilterType = FILTER_DISCBCDEFreal

!MCai, May 2, take care of the following code. How implement Neumann BC?
      ! Do we have Neumann boundary?
      !
      ! The bhasNeumannBoundary flag of the higher level decides about that...
!      bneumann = rACproblem%RlevelInfo(rACproblem%NLMAX)%bhasNeumannBoundary
!      rnonlinearIteration%p_RfilterChain(3)%ifilterType = FILTER_DONOTHING

    end subroutine
  ! ***************************************************************************

!<subroutine>

    subroutine AC_doneNonlinearLoop (rnonlinearIteration)
  
!<description>
  ! Releases memory allocated in cc_initNonlinearLoop.
  ! The routine automatically calls cc_releaseNonlinearLoop to release
  ! internal parameters.
!</description>

!<inputoutput>
  ! The nonlinear iteration structure that should be cleaned up.
    type(t_ACNonlinearIteration), intent(INOUT) :: rnonlinearIteration
!</inputoutput>

!</subroutine>
    
    ! Release the filter chain for the defect vectors.
      if (associated(rnonlinearIteration%p_RfilterChain)) &
        deallocate(rnonlinearIteration%p_RfilterChain)

      rnonlinearIteration%NLMIN = 0
      rnonlinearIteration%NLMAX = 0

    end subroutine
   


end module
