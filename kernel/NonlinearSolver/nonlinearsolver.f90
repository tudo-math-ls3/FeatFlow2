!##############################################################################
!# ****************************************************************************
!# <name> nonlinearsolver </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains a black box nonlinear solver using a defect correction
!# loop. The solver basically performs the following iteration:
!# <tex>
!#     $$  x_{n+1}  =  x_n  +  J^{-1} ( b - A(x_n) x_n )  $$
!# </tex>
!#
!# The defect correction loop is split into three tasks:
!#
!# 1.) Calculation of nonlinear defect
!#
!#        <tex> $$  d_n  :=  b - A(x_n)x_n  $$ </tex>
!#
!#     The calculation of the nonlinear defect is performed user defined.
!#     The nonlinear loop itself does not know anything about the matrix.
!#     In fact, the matrix does not have to exist - this depends on the
!#     problem.
!#
!# 2.) Preconditioning of the nonlinear defect
!#
!#        <tex> $$  u_n  :=  J^{-1} d_n  $$ </tex>
!#
!#     with a preconditioner <tex>$J^{-1}$</tex>. This can be a Jacobian,
!#     a linear solver, a matrix, the inverse of a lumped mass matrix or
!#     similar.
!#
!# 3.) Correction
!#
!#        <tex> $$  x_{n+1}  :=  x_n  +  \omega u_n  $$ </tex>
!#
!#     <tex>$\omega$</tex> is usually = 1.
!#
!#
!# The following routines can be found here:
!#
!# 1.) nlsol_testConvergence
!#     -> Auxiliary routine. Test during the iteration for convergence.
!#
!# 2.) nlsol_testDivergence
!#     -> Auxiliary routine. Test during the iteration for divergence.
!#
!# 3.) nlsol_setPrecMatrix
!#     -> Install a (block) matrix as preconditioner for the defect
!#
!# 4.) nlsol_setPrecMatrixSc
!#     -> Install a scalar matrix as preconditioner for the defect
!#
!# 5.) nlsol_setPrecMass
!#     -> Install a (block) lumped mass matrix as preconditioner for
!#        the defect
!#
!# 6.) nlsol_setPrecMatrixSc
!#     -> Install a scalar lumped mass matrix as preconditioner for
!#        the defect
!#
!# 7.) nlsol_setPrecLinsol
!#     -> Installs a linear solver node as preconditioner for the defect
!#
!# 8.) nlsol_performSolve
!#     -> Invokes the nonlinear solver
!#
!# 9.) nlsol_performSolveSc
!#     -> Invokes the nonlinear solver for a scalar solution/rhs vector
!# </purpose>
!##############################################################################

module nonlinearsolver

!$use omp_lib
  use fsystem
  use genoutput
  use storage
  use linearsolver
  use spatialdiscretisation
  use linearsystemscalar
  use linearsystemblock
  use collection
  
  implicit none
  
  private

  public :: nlsol_testConvergence
  public :: nlsol_testDivergence
  public :: nlsol_setPrecMatrix
  public :: nlsol_setPrecMatrixSc
  public :: nlsol_setPrecMass
  public :: nlsol_setPrecLinsol
  public :: nlsol_performSolve
  public :: nlsol_performSolveSc
  
! *****************************************************************************
! *****************************************************************************
! *****************************************************************************

!<constants>

!<constantblock description="Identifiers for preconditioner type">
  
  ! User defined or no preconditioner (if callback routine is not given).
  integer, parameter, public :: NLSOL_PREC_USERDEF  = 0
  
  ! Preconditioner = Multiplication with a matrix
  integer, parameter, public :: NLSOL_PREC_MATRIX   = 1

  ! Preconditioner = Multiplication with inverse lumped mass matrix
  integer, parameter, public :: NLSOL_PREC_LMASS    = 2
  
  ! Preconditioner = Call to linear solver
  integer, parameter, public :: NLSOL_PREC_LINSOL   = 3
!</constantblock>

!<constantblock description="Identifiers for stopping criterium istoppingCrterium.">

  ! Use standard stopping criterion.
  ! If depsRel>0: use relative stopping criterion.
  ! If depsAbs>0: use abs stopping criterion.
  ! If both are > 0: use both, i.e. the iteration stops when both,
  !    the relative AND the absolute stopping criterium holds
  integer, parameter, public :: NLSOL_STOP_STANDARD     = 0

  ! Use 'minimum' stopping criterion.
  ! If depsRel>0: use relative stopping criterion.
  ! If depsAbs>0: use abs stopping criterion.
  ! If both are > 0: use one of them, i.e. the iteration stops when the
  !    either the relative OR the absolute stopping criterium holds
  integer, parameter, public :: NLSOL_STOP_ONEOF        = 1
  
!</constantblock>

!<constantblock description="Other constants.">
  ! Maximum number of equations supported by the standard error control routine
  ! when checking convergence/divergence criteria.
  integer, parameter, public :: NLSOL_MAXEQUATIONSERROR = 16
!</constantblock>

!</constants>

! *****************************************************************************
! *****************************************************************************
! *****************************************************************************

!<types>

!<typeblock>
  
  ! This is the central structure which repesents the nonlinear solver.
  ! It collects all information for the solution process.
  !
  ! The caller who wants to solve a problem has to initialise
  ! such a structure with all information that is necessary
  ! to solve the nonlinear system.
  !
  ! The structure contains INPUT, OUTPUT and STATUS parameters. INPUT
  ! parameters can be changed by the caller prior to solving.
  ! OUTPUT parameters indicate the final status of the solution process.
  ! STATUS parameters are only valid during the solution
  ! process. They should not be accessed from outside, as they are
  ! maintained internally by the solver.
  ! Parameters marked as READ ONLY are set during the initialisation phase
  ! of the solver and must not be changed by the caller.
  
  type t_nlsolNode
    
    ! OUTPUT: Result
    ! The result of the solution process.
    ! =-1: convergence criterion not reached.
    ! =0: success.
    ! =1: iteration broke down, diverging.
    ! =2: iteration broke down, preconditioner did not work.
    ! =3: error in the parameters.
    integer                    :: iresult
    
    ! OUTPUT: Number of performed iterations, if the solver
    ! is of iterative nature.
    ! Is to 1 by the solver if not used (indicating at least 1 performed
    ! iteration, which is always the case).
    integer                    :: iiterations
    
    ! OUTPUT PARAMETER:
    ! Norm of initial residuum for each subvector of the given block vector.
    ! Only valid, if the fcb_resNormCheck callback procedure is not specified
    ! in the call to the solver, otherwise undefined!
    real(DP), dimension(NLSOL_MAXEQUATIONSERROR)  :: DinitialDefect

    ! OUTPUT PARAMETER:
    ! Norm of final residuum for each subvector of the given block vector.
    ! Only valid, if the fcb_resNormCheck callback procedure is not specified
    ! in the call to the solver, otherwise undefined!
    real(DP), dimension(NLSOL_MAXEQUATIONSERROR)  :: DfinalDefect

    ! OUTPUT PARAMETER:
    ! Norm of initial residuum for the complete solution vector.
    ! Only valid, if the fcb_resNormCheck callback procedure is not specified
    ! in the call to the solver, otherwise undefined!
    real(DP) :: dinitialDefectTotal = 0.0_DP

    ! OUTPUT PARAMETER:
    ! Norm of final residuum for the complete solution vector.
    ! Only valid, if the fcb_resNormCheck callback procedure is not specified
    ! in the call to the solver, otherwise undefined!
    real(DP) :: dfinalDefectTotal = 0.0_DP

    ! OUTPUT PARAMETER:
    ! Total time for nonlinear solver
    real(DP)                        :: dtimeTotal

    ! OUTPUT PARAMETER:
    ! Total time for calculating nonlinear defects
    real(DP)                        :: dtimeNLdefect

    ! OUTPUT PARAMETER:
    ! Total time for nonlinear preconditioning
    real(DP)                        :: dtimeNLpreconditioning
    
    ! INPUT PARAMETER:
    ! Damping parameter for defect correction.
    ! Standard value = 1.0 (corresponds to 'no damping')
    real(DP)                        :: domega  = 1.0_DP

    ! INPUT PARAMETER:
    ! Standard relative stopping criterion for each subvector of a block vector.
    ! Only valid, if there is no callback routine for the calculation of
    ! residuals passed to the nonlinear solver.
    ! Stop iteration if everywhere
    !   !!defect!! < EPSREL * !!initial defect!!.
    ! =0: ignore, use absolute stopping criterion; standard = 1E-5
    ! Remark: do not set depsAbs=depsRel=0!
    real(DP), dimension(NLSOL_MAXEQUATIONSERROR) :: DepsRel = 1E-5_DP

    ! INPUT PARAMETER:
    ! Standard absolute stopping criterion for each subvector of a block vector.
    ! Only valid, if there is no callback routine for the calculation of
    ! residuals passed to the nonlinear solver.
    ! Stop iteration if everywhere
    !   !!defect!! < EPSREL.
    ! =0: ignore, use relative stopping criterion; standard = 1E-5
    ! Remark: do not set depsAbs=depsRel=0!
    real(DP), dimension(NLSOL_MAXEQUATIONSERROR) :: DepsAbs = 1E-5_DP

    ! INPUT PARAMETER:
    ! Standard relative stopping criterion for the complete solution vector.
    ! Only valid, if there is no callback routine for the calculation of
    ! residuals passed to the nonlinear solver.
    ! Stop iteration if everywhere
    !   !!defect!! < EPSREL * !!initial defect!!.
    ! =0: ignore, use absolute stopping criterion; standard = 0
    ! Remark: do not set depsAbs=depsRel=0!
    real(DP) :: depsRelTotal = 0.0_DP

    ! INPUT PARAMETER:
    ! Standard absolute stopping criterion for the complete solution vector.
    ! Only valid, if there is no callback routine for the calculation of
    ! residuals passed to the nonlinear solver.
    ! Stop iteration if everywhere
    !   !!defect!! < EPSREL.
    ! =0: ignore, use relative stopping criterion; standard = 0
    ! Remark: do not set depsAbsTotal=depsRelTotal=0!
    real(DP) :: depsAbsTotal = 0.0_DP

    ! INPUT PARAMETER:
    ! Standard relative divergence criterion for each subvector of a block vector.
    ! Only valid, if there is no callback routine for the calculation of
    ! residuals passed to the nonlinear solver.
    ! Treat iteration as diverged if anywhere
    !   !!defect!! >= DIVREL * !!initial defect!!
    ! A value of SYS_INFINITY_DP disables the relative divergence check.
    ! standard = 1E3
    real(DP), dimension(NLSOL_MAXEQUATIONSERROR) :: DdivRel = 1E3_DP

    ! INPUT PARAMETER:
    ! Standard absolute divergence criterion for each subvector of a block vector.
    ! Only valid, if there is no callback routine for the calculation of
    ! residuals passed to the nonlinear solver.
    ! Treat iteration as diverged if anywhere
    !   !!defect!! >= DIVABS
    ! A value of SYS_INFINITY_DP disables the absolute divergence check.
    ! standard = SYS_INFINITY_DP
    real(DP), dimension(NLSOL_MAXEQUATIONSERROR) :: DdivAbs = SYS_INFINITY_DP

    ! INPUT PARAMETER:
    ! Type of stopping criterion to use for standard convergence test. One of the
    ! NLSOL_STOP_xxxx constants.
    ! Note: This parameter is only evaluated in the stanard convergence test.
    ! If the caller of the nonlinear solver specifies a callback routine fcb_resNormCheck
    ! for checking the convergence, that callback routine must implement its own
    ! logic to handle relative and absolute convrgence criteria!
    integer                    :: istoppingCriterion = NLSOL_STOP_STANDARD

    ! INPUT PARAMETER:
    ! Minimum number of iterations top perform
    integer                    :: nminIterations = 1

    ! INPUT PARAMETER:
    ! Maximum number of iterations top perform
    integer                    :: nmaxIterations = 50
    
    ! INPUT PARAMETER FOR SOLVERS WITH RESIDUAL CHECK:
    ! For every subvector: Type of norm to use in the residual checking
    ! (cf. linearalgebra.f90).
    ! =0: euclidian norm, =1: l1-norm, =2: l2-norm, =3: MAX-norm
    integer, dimension(NLSOL_MAXEQUATIONSERROR) :: IresNorm = 2

    ! INPUT PARAMETER FOR SOLVERS WITH RESIDUAL CHECK:
    ! Type of norm to use in the residual checking of the total vector
    ! (cf. linearalgebra.f90).
    ! =0: euclidian norm, =1: l1-norm, =2: l2-norm, =3: MAX-norm
    integer :: iresNormTotal = 2
    
    ! INPUT PARAMETER: Output level
    ! This determines the output level of the solver.
    ! =0: no output, =1: basic output, =2, extended output
    integer                    :: ioutputLevel = 0
    
    ! INPUT PARAMETER: Output mode. Used for printing messages.
    ! =OU_MODE_STD: Print messages to the terminal and probably to a log
    ! file (if a log file is opened).
    integer(I32)               :: coutputmode = OU_MODE_STD

    ! INPUT PARAMETER:
    ! Type of preconditioner to use. One of the NLSOL_PREC_xxxx constants.
    integer                    :: cpreconditioner = NLSOL_PREC_USERDEF
    
    ! INPUT PARAMETER:
    ! If cpreconditioner=NLSOL_PREC_MATRIX block matrix that is multiplied
    ! to the defect vector.
    ! If cpreconditioner=NLSOL_PREC_LMASS: Lumped mass matrix. Only matrix
    ! format 7 and 9 are supported.
    type(t_matrixBlock)        :: rmatrixPreconditioner
    
    ! INPUT PARAMETER:
    ! If cpreconditioner=NLSOL_PREC_LINSOL: Pointer to linear solver node
    ! defining the preconditioner
    type(t_linsolNode), pointer :: p_rlinsolNode => null()
    
    ! STATUS:
    ! Current iteration
    integer                     :: icurrentIteration

  end type
  
  public :: t_nlsolNode
  
!</typeblock>
  
!</types>

! *****************************************************************************
! *****************************************************************************
! *****************************************************************************

contains

  ! ***************************************************************************

!<function>
  
  logical function nlsol_testConvergence (rsolverNode, DvecNorm, dvecNormTotal,&
      nvectorblocks, rdef) result(loutput)
  
!<description>
  
  ! Tests a defect vector rdef whether it is in a defined tolerance configured
  ! in the solver node, so the iteration of an iterative solver can be
  ! can be treated as 'converged'.
  ! The iteration is treated as 'converged' if both, the relative and the
  ! absolute convergence criterion are fulfilled (or only one of them,
  ! respectively, if the other is switched off).
  !
  ! The solver must have initialised the 'DinitialDefect' variable
  ! of the solver structure for this routine to work properly!
  
!</description>
  
!<result>
  ! Boolean value. =TRUE if the convergence criterion is reached;
  ! =FALSE otherwise.
!</result>
  
!<input>
  
  ! The solver node of the nonlinear solver that contains the
  ! convergence criterion
  type(t_nlsolNode), intent(in) :: rsolverNode
  
  ! Number of blocks in the solution vector/equation
  integer, intent(in) :: nvectorblocks
  
  ! OPTIONAL: The defect vector which norm should be tested.
  ! If existent, the norm of the subvectors is returned in DvecNorm.
  ! If not existent, the routine assumes that DvecNorm is the norm
  ! of the vector and checks convergence depending on DvecNorm.
  type(t_vectorBlock), intent(in), optional :: rdef
  
!</input>

!<inputoutput>
  
  ! Norm of the subvectors on the defect vector.
  ! If rdef if present, the routine will calculate the norm of the subvectors
  ! and return them in dvecNorm.
  ! If rdef is not present, DvecNorm is assumed to be a valid norm of a
  ! vector and convergence is tested using DvecNorm.
  real(DP), dimension(:), intent(inout) :: DvecNorm

  ! Total Norm of the defect vector.
  ! If rdef if present, the routine will calculate the norm of the subvectors
  ! and return them in dvecNormTotal.
  ! If rdef is not present, DvecNorm is assumed to be a valid norm of a
  ! vector and convergence is tested using DvecNorm.
  real(DP), intent(inout) :: dvecNormTotal

!</inputoutput>
  
!</function>

    ! local variables
    integer :: i,nblocks
    logical :: bok

    nblocks = min(NLSOL_MAXEQUATIONSERROR,nvectorblocks)

    ! Calculate the norm of the vector or take the one given
    ! as parameter
    if (present(rdef)) then
      call lsysbl_vectorNormBlock (rdef,rsolverNode%IresNorm(1:nblocks),&
          DvecNorm(1:nblocks))
      where (.not.((DvecNorm .ge. 1E-99_DP) .and. (DvecNorm .le. 1E99_DP)))
        DvecNorm = 0.0_DP
      end where
      dvecNormTotal = lsysbl_vectorNorm (rdef,rsolverNode%iresNormTotal)
      if (.not.((dvecNormTotal .ge. 1E-99_DP) .and. (dvecNormTotal .le. 1E99_DP))) &
        dvecNormTotal = 0.0_DP
    end if
  
    select case (rsolverNode%istoppingCriterion)
    
    case (LINSOL_STOP_ONEOF)
      ! Iteration stops if either the absolute or the relative criterium holds.
      loutput = .false.
      
      ! Absolute convergence criterion? Check the norm directly.
      if (rsolverNode%DepsAbs(i) .ne. 0.0_DP) then
        bok = .false.
        do i=1,nblocks
          bok = bok .or. (.not. (DvecNorm(i) .gt. rsolverNode%DepsAbs(i)))
        end do
        if (bok) then
          loutput = .true.
          return
        end if
      end if
        
      ! Relative convergence criterion? Multiply with initial residuum
      ! and check the norm.
      ! If the initial defect is 0, this rule is not to be applied!
      ! (May happen in the first step of 'flow around cylinder'
      !  where we have only convection in the X-direction, not in the Y-direction
      !  and a still fluid.)
      if ((rsolverNode%DepsRel(i) .ne. 0.0_DP) .and. &
          (rsolverNode%dinitialDefect(i) .gt. SYS_EPSREAL_DP)) then
        bok = .false.
        do i=1,nblocks
          bok = bok .or. (.not. (DvecNorm(i) .gt. &
              rsolverNode%depsRel(i) * rsolverNode%dinitialDefect(i)))
        end do
        if (bok) then
          loutput = .true.
          return
        end if
        
      end if

    case DEFAULT
      ! Standard stopping criterion.
      ! Iteration stops if both the absolute and the relative criterium holds.
      loutput = .true.
      
      do i=1,nblocks
        ! Absolute convergence criterion? Check the norm directly.
        if (rsolverNode%DepsAbs(i) .ne. 0.0_DP) then
          if (DvecNorm(i) .gt. rsolverNode%DepsAbs(i)) then
            loutput = .false.
            return
          end if
        end if
        
        ! Relative convergence criterion? Multiply with initial residuum
        ! and check the norm.
        ! If the initial defect is 0, this rule is not to be applied!
        ! (May happen in the first step of 'flow around cylinder'
        !  where we have only convection in the X-direction, not in the Y-direction
        !  and a still fluid.)
        if ((rsolverNode%DepsRel(i) .ne. 0.0_DP) .and. &
            (rsolverNode%dinitialDefect(i) .gt. SYS_EPSREAL_DP)) then
          if (DvecNorm(i) .gt. &
              rsolverNode%depsRel(i) * rsolverNode%dinitialDefect(i)) then
            loutput = .false.
            return
          end if
        end if
        
      end do ! i
    end select

    ! Absolute convergence criterion?
    if (rsolverNode%depsAbsTotal .ne. 0.0_DP) then
      if (.not. (DvecNormTotal .gt. rsolverNode%depsAbsTotal)) then
        loutput = .true.
      end if
    end if

    ! Relative convergence of the total vector.
    if ((rsolverNode%depsRelTotal .ne. 0.0_DP) .and. &
        (rsolverNode%dinitialDefectTotal .gt. SYS_EPSREAL_DP)) then
      if (.not. (dvecNormTotal .gt. &
            rsolverNode%depsRelTotal * rsolverNode%dinitialDefectTotal)) then
        loutput = .true.
      end if
      
    end if
  
  end function
  
  ! ***************************************************************************

!<function>
  
  logical function nlsol_testDivergence (rsolverNode, DvecNorm, nvectorblocks, rdef) &
          result(loutput)
  
!<description>
  
  ! Tests a defect vector rx whether it is out of a defined tolerance configured
  ! in the solver node, so the iteration of an iterative solver can be
  ! can be treated as 'diverged'.
  ! The iteration is treated as 'diverged' if one criterion, the relative or the
  ! absolute divergence criterion is fulfilled.
  !
  ! The solver must have initialised the 'DinitialDefect' variable
  ! of the solver structure for this routine to work properly!

!</description>
  
  !<result>
  ! Boolean value. =TRUE if the divergence criterion is reached;
  ! =FALSE otherwise.
  !</result>
  
!<input>
  
  ! The solver node of the nonlinear solver that contains the
  ! convergence criterion
  type(t_nlsolNode), intent(in) :: rsolverNode
  
  ! Number of blocks in the solution vector/equation
  integer, intent(in) :: nvectorblocks
  
  ! OPTIONAL: The defect vector which norm should be tested.
  ! If existent, the norm of the subvectors is returned in DvecNorm.
  ! If not existent, the routine assumes that DvecNrm is the norm
  ! of the vector and checks convergence depending on DvecNorm.
  type(t_vectorBlock), intent(in), optional :: rdef
  
!</input>

!<inputoutput>
  
  ! Norm of the subvectors on the defect vector.
  ! If rdef if present, the routine will calculate the norm of the subvectors
  ! and return them in dvecNorm.
  ! If rdef is not present, DvecNorm is assumed to be a valid norm of a
  ! vector and convergence is tested using DvecNorm.
  real(DP), dimension(:), intent(inout) :: DvecNorm

!</inputoutput>
  
!</function>

  ! local variables
  integer :: i,nblocks

  nblocks = min(NLSOL_MAXEQUATIONSERROR,nvectorblocks)

  ! Calculate the norm of the vector if not given
  ! as parameter
  if (present(rdef)) then
    call lsysbl_vectorNormBlock (rdef,rsolverNode%IresNorm(1:nblocks),&
        DvecNorm(1:nblocks))
    where (.not.((DvecNorm .ge. 1D-99) .and. (DvecNorm .le. 1D99)))
      DvecNorm = 0.0_DP
    end where
  end if
  
  loutput = .false.
  
  do i=1,nblocks
  
    ! Absolute divergence criterion? Check the norm directly.
    if (rsolverNode%DdivAbs(i) .ne. SYS_INFINITY_DP) then
     
      ! use NOT here - gives a better handling of special cases like NaN!
      if ( .not. (DvecNorm(i) .le. rsolverNode%DdivAbs(i))) then
        loutput = .true.
        return
      end if
      
    end if
    
    ! Relative divergence criterion? Multiply with initial residuum
    ! and check the norm.
    ! If the initial defect is 0, this rule is not to be applied!
    ! (May happen in the first step of 'flow around cylinder'
    !  where we have only convection in the X-direction, not in the Y-direction
    !  and a still fluid.)
    if ((rsolverNode%DinitialDefect(i) .ne. 0.0_DP) .and. &
        (rsolverNode%DepsRel(i) .ne. SYS_INFINITY_DP)) then
      if ( .not. (DvecNorm(i) .le. &
          rsolverNode%DinitialDefect(i) * rsolverNode%DdivRel(i)) ) then
        loutput = .true.
        return
      end if
    end if
  
  end do ! i
  
  end function

  ! ***************************************************************************

!<subroutine>

  subroutine nlsol_setPrecMatrix (rsolverNode,rmatrix)
  
!<description>
  ! This routine installs a standard matrix preconditioner into the nonlinear
  ! subroutine. rmatrix =: P is a block matrix which is multiplied to the defect
  ! vector in every iteration:
  !   <tex> $$  x_{n+1}  =  x_n  +  P (b-A(x)x)  $$ </tex>
!</description>

!<input>
  ! The block matrix.
  type(t_matrixBlock), intent(in) :: rmatrix
!</input>

!<inputoutput>
  ! The solver node of the nonlinear iteration loop where to set the
  ! preconditioner.
  type(t_nlsolNode), intent(inout) :: rsolverNode
!</inputoutput>

!</subroutine>

    ! Install the matrix. Copy the matrix structure into the solver node.
    rsolverNode%rmatrixPreconditioner = rmatrix
    rsolverNode%rmatrixPreconditioner%RmatrixBlock%imatrixSpec = &
      ior(rsolverNode%rmatrixPreconditioner%RmatrixBlock%imatrixSpec,LSYSSC_MSPEC_ISCOPY)
      
    ! Set the preconditioner type flag appropriately.
    rsolverNode%cpreconditioner = NLSOL_PREC_MATRIX

  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine nlsol_setPrecMatrixSc (rsolverNode,rmatrix)
  
!<description>
  ! This routine installs a standard matrix preconditioner into the nonlinear
  ! subroutine. rmatrix is a scalar matrix. The matrix is added as 1x1 block
  ! matrix preconditioner P := [rmatrix] to the solver node and is multiplied
  ! to the defect vector in every iteration:
  !    <tex> $$  x_{n+1}  =  x_n  +  P (b-A(x)x)  $$ </tex>
!</description>

!<input>
  ! The scalar matrix.
  type(t_matrixScalar), intent(in) :: rmatrix
!</input>

!<inputoutput>
  ! The solver node of the nonlinear iteration loop where to set the
  ! preconditioner.
  type(t_nlsolNode), intent(inout) :: rsolverNode
!</inputoutput>

!</subroutine>

    ! Create a 1x1 block matrix from rmatrix and save it to the solver node.
    call lsysbl_createMatFromScalar (rmatrix,rsolverNode%rmatrixPreconditioner)
    rsolverNode%rmatrixPreconditioner%RmatrixBlock%imatrixSpec = &
      ior(rsolverNode%rmatrixPreconditioner%RmatrixBlock%imatrixSpec,LSYSSC_MSPEC_ISCOPY)
      
    ! Set the preconditioner type flag appropriately.
    rsolverNode%cpreconditioner = NLSOL_PREC_MATRIX

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine nlsol_setPrecMass (rsolverNode,rmatrix)
  
!<description>
  ! This routine installs a lumped mass matrix preconditioner into the nonlinear
  ! subroutine. rmatrix =: M is a block matrix which contains only diagonal
  ! blocks with lumped mass matrices. this matrix is is multiplied to the defect
  ! vector in every iteration:
  !    <tex> $$  x_{n+1}  =  x_n  +  M^{-1} (b-A(x)x)  $$ </tex>
!</description>

!<input>
  ! The block matrix. Must be a matrix only containing diagonal blocks.
  ! the diagonal blocks must be lumped mass matrices.
  type(t_matrixBlock), intent(in) :: rmatrix
!</input>

!<inputoutput>
  ! The solver node of the nonlinear iteration loop where to set the
  ! preconditioner.
  type(t_nlsolNode), intent(inout) :: rsolverNode
!</inputoutput>

!</subroutine>

    ! Install the matrix.
    call nlsol_setPrecMatrix (rsolverNode,rmatrix)
      
    ! Set the preconditioner type flag appropriately.
    rsolverNode%cpreconditioner = NLSOL_PREC_LMASS

  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine nlsol_setPrecMassSc (rsolverNode,rmatrix)
  
!<description>
  ! This routine installs a standard matrix preconditioner into the nonlinear
  ! subroutine. rmatrix is a scalar matrix. The matrix is added as 1x1 block
  ! matrix preconditioner P := [rmatrix] to the solver node and its inverse
  ! is multiplied to the defect vector in every iteration:
  !    <tex> $$  x_{n+1}  =  x_n  +  M^{-1} (b-A(x)x)  $$ </tex>
!</description>

!<input>
  ! The lumped mass matrix.
  type(t_matrixScalar), intent(in) :: rmatrix
!</input>

!<inputoutput>
  ! The solver node of the nonlinear iteration loop where to set the
  ! preconditioner.
  type(t_nlsolNode), intent(inout) :: rsolverNode
!</inputoutput>

!</subroutine>

    ! Install the matrix
    call nlsol_setPrecMatrixSc (rsolverNode,rmatrix)
      
    ! Set the preconditioner type flag appropriately.
    rsolverNode%cpreconditioner = NLSOL_PREC_LMASS

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine nlsol_setPrecLinsol (rsolverNode,rlinsolNode)
  
!<description>
  ! This routine installs a linear solver node as preconditioner to the
  ! nonlinear solver. The linear solver <tex>$J^{-1}$</tex> is called after
  ! each defect calculation to precondition the defect:
  !    <tex> $$  x_{n+1}  =  x_n  +  J^{-1} (b-A(x)x)  $$ </tex>
!</description>

!<input>
  ! The linear solver node identifying the preconditioner
  type(t_linsolNode), intent(in), target :: rlinsolNode
!</input>

!<inputoutput>
  ! The solver node of the nonlinear iteration loop where to set the
  ! preconditioner.
  type(t_nlsolNode), intent(inout) :: rsolverNode
!</inputoutput>

!</subroutine>

    ! Install the preconditioner
    rsolverNode%p_rlinsolNode => rlinsolNode
      
    ! Set the preconditioner type flag appropriately.
    rsolverNode%cpreconditioner = NLSOL_PREC_LINSOL

  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine nlsol_performSolve (rsolverNode,rx,rb,rd,&
                                 fcb_getDefect,fcb_precondDefect,fcb_resNormCheck,&
                                 rcollection)
             
!<description>
  ! This routine invokes the nonlinear defect correction iteration
  !    <tex> $$  x_{n+1}  =  x_n  +  J^{-1} ( b - A(x_n) x_n )  $$ </tex>
  !
  ! The defect correction loop is split into three tasks:
  !
  ! 1.) Calculation of nonlinear defect
  !        <tex> $$  d_n  :=  b - A(x_n)x_n  $$ </tex>
  !     For this purpose, the callback routine fcb_getDefect is called.
  !
  ! 2.) Preconditioning of the nonlinear defect
  !        <tex> $$  u_n  :=  J^{-1} d_n  $$ </tex>
  !     with a preconditioner <tex>$J^{-1}$</tex>. The preconditioner can be
  !     a) a matrix (if rsolverNode%cpreconditioner=NLSOL_PREC_MATRIX) or
  !     b) a lumped mass matrix (if rsolverNode%cpreconditioner=NLSOL_PREC_LMASS) or
  !     c) a linear solver (if rsolverNode%cpreconditioner=NLSOL_PREC_LINSOL) or
  !     d) no preconditioner (if rsolverNode%cpreconditioner=NLSOL_PREC_USERDEF and
  !        fcb_precondDefect is not present) or
  !     e) a user defined preconditioner, realised by the callback routine
  !        fcb_precondDefect (if rsolverNode%cpreconditioner=NLSOL_PREC_USERDEF and
  !        fcb_precondDefect is present). In this case, fcb_precondDefect is
  !        called for the preconditioning.
  !
  ! 3.) Correction
  !         <tex> $$  x_{n+1}  :=  x_n  +  \omega u_n  $$ </tex>
  !     <tex>$\omega$</tex> is usually = 1. A user defined callback routine fcb_precondDefect
  !     can modify <tex>$\omega$</tex> in every iteration.
  !
  ! If present, the callback routine fcb_resNormCheck is called in every iteration.
  ! This routine can check the residuum for convergence/divergence and can print
  ! information about the norm of the residuum to screen.
  !
  ! At the beginning of the nonlinear iteration, the callback routines
  ! fcb_getResiduum and fcb_resNormCheck are called once with ite=0 to calculate
  ! the initial defect, its norm and to check if already the initial vector
  ! is the solution.
  !
  ! If a linear solver is chosen for preconditioning, the nonlinear solver
  ! assumes that this is completely initialised by the application and ready
  ! for action; it will directly call linsol_precondDefect without any further
  ! initialisation!
!</description>

!<inputoutput>
  ! The nonlinear solver node that configures the solution process.
  type(t_nlsolNode), intent(inout)              :: rsolverNode
  
  ! INPUT: Initial solution vector.
  ! OUTPUT: Final iteration vector.
  type(t_vectorBlock), intent(inout)            :: rx
             
  ! Temporary vector. Must be of the same size/type as rx/rb.
  type(t_vectorBlock), intent(inout)            :: rd
  
  ! OPTIONAL: Collection structure that saves problem-dependent information.
  ! This is passed without being changed to the callback routines of this
  ! algorithm.
  type(t_collection), intent(inout), target, optional :: rcollection
  
!</inputoutput>
  
!<input>
  ! Right hand side vector of the equation.
  type(t_vectorBlock), intent(in)               :: rb
  
  ! Defect vector calculation routine. Based on the current iteration vector
  ! rx and the right hand side vector rb, this routine has to compute the
  ! defect vector rd.
  include 'intf_nlsolcallback.inc'
  
  ! OPTIONAL: Preconditioning routine. This routine accepts a defect vector rd
  ! and replaces it by a preconditioned defect vector <tex>$J^{-1} rd$</tex>.
  ! If this parameter is not present, the preconditioner is either a matrix
  ! or there is no preconditioner (depending on the variable
  ! rsolverNode\%cpreconditioner).
  optional :: fcb_precondDefect
  
  ! OPTIONAL: Residual norm calculation and printing routine.
  ! If not present, the standard absolute/relative stopping criteria of the
  !  solver node with the configuration of the nonlinear solver are used.
  ! If present, this routine is called each time the norm of the residuum was
  !  calculated. It has to check the current residuum for convergence and/or
  !  divergence and can print the residuum to screen.
  optional :: fcb_resNormCheck
  
!</input>

!</subroutine>

  ! local variables
  type(t_collection), pointer :: p_rcollection
  integer :: ite,i
  real(DP), dimension(NLSOL_MAXEQUATIONSERROR) :: DvecNorm
  real(DP) :: dvecNormTotal
  type(t_vectorBlock) :: rtemp
  real(DP) :: domega
  integer :: nblocks
  logical :: bconvergence,bdivergence,bsuccess
  
    ! Do we have a collection?
    nullify(p_rcollection)
    if (present(rcollection)) p_rcollection => rcollection
    
    ! In case our preconditioner is a matrix-vector multiplication,
    ! allocate memory for another temporary vector used
    ! during the MV multiplication.
    if (rsolverNode%cpreconditioner .eq. NLSOL_PREC_MATRIX) then
      call lsysbl_createVecBlockIndirect (rx,rtemp,.false.)
    end if
    
    ! Status reset
    rsolverNode%iresult = 0
    
    ! Calculate the initial nonlinear defect to rd:   d = b-A(x)x
    call fcb_getDefect (0,rx,rb,rd,p_rcollection)
    
    ite = 0
    rsolverNode%icurrentIteration = ite

    ! The standard convergence/divergence test supports only up to
    ! NLSOL_MAXEQUATIONSERROR equations.
    nblocks = min(rb%nblocks,NLSOL_MAXEQUATIONSERROR)
    
    ! Initial test for convergence/divergence.
    if (present(fcb_resNormCheck)) then
      ! Calculate the defect with the callback routine
      call fcb_resNormCheck (ite,rx,rb,rd,bconvergence,bdivergence,p_rcollection)
    else
      ! Calculate the norm of the defect:
      DvecNorm = 0.0_DP
      call lsysbl_vectorNormBlock (rd,rsolverNode%IresNorm(1:nblocks),&
          DvecNorm(1:nblocks))
      where (.not.((DvecNorm .ge. 1D-99) .and. (DvecNorm .le. 1D99)))
        DvecNorm = 0.0_DP
      end where
      rsolverNode%DinitialDefect = DvecNorm
      rsolverNode%DfinalDefect = DvecNorm

      dvecNormTotal = lsysbl_vectorNorm (rd,rsolverNode%iresNormTotal)
      if (.not.((dvecNormTotal .ge. 1D-99) .and. (dvecNormTotal .le. 1D99))) &
        dvecNormTotal = 0.0_DP
      rsolverNode%dinitialDefectTotal = dvecNormTotal
      rsolverNode%dfinalDefectTotal = dvecNormTotal
      
      bconvergence = nlsol_testConvergence (rsolverNode, DvecNorm, dvecNormTotal, nblocks)
      bdivergence  = nlsol_testDivergence (rsolverNode, DvecNorm, nblocks)
      
      if (rsolverNode%ioutputLevel .ge. 2) then
        call output_line ('NLSOL: Iteration '//&
             trim(sys_siL(ite,10))//', !!RES!! =',bnolinebreak=.true.,&
             coutputMode=rsolverNode%coutputMode)
        do i=1,nblocks
          call output_line (' '//trim(sys_sdEL(DvecNorm(i),15)),bnolinebreak=.true.,&
              coutputMode=rsolverNode%coutputMode,cdateTimeLogPolicy=OU_DTP_NONE)
        end do
        call output_lbrk()
      end if
    end if

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
        select case (rsolverNode%cpreconditioner)
        case (NLSOL_PREC_MATRIX)
          ! Multiplication with a matrix.
          call lsysbl_copyVector (rd,rtemp)
          call lsysbl_blockMatVec (rsolverNode%rmatrixPreconditioner, &
                                   rtemp, rd, 1.0_DP, 0.0_DP)
                                   
        case (NLSOL_PREC_LMASS)
          ! Multiply with the inverse of the diagonal in the block matrix.
          call lsysbl_invertedDiagMatVec (rsolverNode%rmatrixPreconditioner,&
                                          rd,1.0_DP,rd)
                                          
        case (NLSOL_PREC_LINSOL)
          ! Preconditioner is a linear solver with fixed matrix.
          call linsol_precondDefect(rsolverNode%p_rlinsolNode,rd)
          bsuccess = rsolverNode%p_rlinsolNode%iresult .eq. 0
          
        case DEFAULT
          ! User defined or no preconditioning.
          ! Is a callback function available?
          if (present(fcb_precondDefect)) then
            ! The callback routine is allowed to change domega during the
            ! iteration if necessary. The nonlinear solver here does not touch
            ! domega anymore, so the callback routine is the only one changing it.
            call fcb_precondDefect (ite,rd,rx,rb,domega,bsuccess,p_rcollection)
          end if
          
        end select
        
        ! If bsuccess=false, the preconditioner had an error.
        if (.not. bsuccess) then
          call output_line ('NLSOL: Iteration '//&
              trim(sys_siL(ite,10))//' canceled as the preconditioner went down!',&
              coutputMode=rsolverNode%coutputMode)
          rsolverNode%iresult = 3
          exit
        end if
        
        ! If domega=0.0, the solution vector would stay unchanged. In this
        ! case, the nonlinear solver would not proceed at all, and the next
        ! iteration would behave exactly as before!
        ! So in this case, there is nothing to do, we can stop the iteration.
        if (domega .eq. 0.0_DP) then
          call output_line ('NLSOL: Iteration '//&
              trim(sys_siL(ite,10))//' canceled as there is no progress anymore!',&
              coutputMode=rsolverNode%coutputMode)
          exit
        else
          ! Add the correction vector in rd to rx;
          ! damped by the damping parameter:           x := x + domega u
          call lsysbl_vectorLinearComb (rd,rx,domega,1.0_DP)
          
          ! Calculate the new nonlinear defect to rd:  d = b-A(x)x
          call fcb_getDefect (ite,rx,rb,rd,p_rcollection)

          ! Check the defect for convergence
          if (present(fcb_resNormCheck)) then
            ! Calculate the defect with the callback routine
            call fcb_resNormCheck (ite,rx,rb,rd,bconvergence,bdivergence,p_rcollection)
          else
            ! Calculate the norm of the defect:
            call lsysbl_vectorNormBlock (rd,rsolverNode%IresNorm(1:nblocks),&
                DvecNorm(1:nblocks))
            where (.not.((DvecNorm .ge. 1E-99_DP) .and. (DvecNorm .le. 1E99_DP)))
              DvecNorm = 0.0_DP
            end where
            rsolverNode%DfinalDefect = DvecNorm

            dvecNormTotal = lsysbl_vectorNorm (rd,rsolverNode%iresNormTotal)
            if (.not.((dvecNormTotal .ge. 1E-99_DP) .and. (dvecNormTotal .le. 1E99_DP))) &
              dvecNormTotal = 0.0_DP
            rsolverNode%dfinalDefectTotal = dvecNormTotal
            
            bconvergence = nlsol_testConvergence (rsolverNode, DvecNorm, dvecNormTotal, nblocks)
            bdivergence  = nlsol_testDivergence (rsolverNode, DvecNorm, nblocks)

            if (rsolverNode%ioutputLevel .ge. 2) then
              call output_line ('NLSOL: Iteration '//&
                  trim(sys_siL(ite,10))//', !!RES!! =',bnolinebreak=.true.,&
                  coutputMode=rsolverNode%coutputMode)
              do i=1,nblocks
                call output_line (' '//trim(sys_sdEL(DvecNorm(i),15)),bnolinebreak=.true.,&
                coutputMode=rsolverNode%coutputMode,cdateTimeLogPolicy=OU_DTP_NONE)
              end do
              call output_lbrk(cdateTimeLogPolicy=OU_DTP_NONE)
            end if
          end if

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

  subroutine nlsol_performSolveSc (rsolverNode,rx,rb,rd,&
                          fcb_getDefect,fcb_precondDefect,fcb_resNormCheck,&
                          rcollection)
             
!<description>
  ! This routine performs the same task as nlsol_performSolve, but for
  ! scalar systems instead of block systems. For this purpose,
  ! the scalar vectors rx,rb and rd are converted temporarily into block
  ! vectors; then nlsol_performSolve is called to perform the nonlinear
  ! iteration.
!</description>

!<inputoutput>
  ! The nonlinear solver node that configures the solution process.
  type(t_nlsolNode), intent(inout)              :: rsolverNode
  
  ! INPUT: Initial solution vector.
  ! OUTPUT: Final iteration vector.
  type(t_vectorScalar), intent(inout)            :: rx
             
  ! Temporary vector. Must be of the same size/type as rx/rb.
  type(t_vectorScalar), intent(inout)            :: rd
  
  ! OPTIONAL: Collection structure that saves problem-dependent information.
  ! This is passed without being changed to the callback routines of this
  ! algorithm.
  type(t_collection), intent(inout), target, optional :: rcollection
  
!</inputoutput>
  
!<input>
  ! Right hand side vector of the equation.
  type(t_vectorScalar), intent(in)               :: rb
  
  ! Defect vector calculation routine. Based on the current iteration vector
  ! rx and the right hand side vector rb, this routine has to compute the
  ! defect vector rd.
  include 'intf_nlsolcallback.inc'
  
  ! OPTIONAL: Preconditioning routine. This routine accepts a defect vector rd
  ! and replaces it by a preconditioned defect vector <tex>$J^{-1} rd$</tex>.
  ! If this parameter is not present, the preconditioner is either a matrix
  ! or there is no preconditioner (depending on the variable
  ! rsolverNode\%cpreconditioner).
  optional :: fcb_precondDefect
  
  ! OPTIONAL: Residual norm calculation and printing routine.
  ! If not present, the standard absolute/relative stopping criteria of the
  !  solver node with the configuration of the nonlinear solver are used.
  ! If present, this routine is called each time the norm of the residuum was
  !  calculated. It has to check the current residuum for convergence and/or
  !  divergence and can print the residuum to screen.
  optional :: fcb_resNormCheck
  
!</input>

!</subroutine>

    ! local variables
    type(t_vectorBlock) :: rxBlock,rbBlock,rdBlock
    
    ! Convert the vectors on-the-fly to block vectors.
    ! The new vectors share the same memory as the old, so the solver will use
    ! and overwrite the old input vectors.
    call lsysbl_createVecFromScalar (rx,rxBlock)
    call lsysbl_createVecFromScalar (rb,rbBlock)
    call lsysbl_createVecFromScalar (rd,rdBlock)
    
    ! Invoke the solver - that is all.
    call nlsol_performSolve (rsolverNode,rxBlock,rbBlock,rdBlock,&
                            fcb_getDefect,fcb_precondDefect,fcb_resNormCheck,&
                            rcollection)
                          
    call lsysbl_releaseVector (rdBlock)
    call lsysbl_releaseVector (rbBlock)
    call lsysbl_releaseVector (rxBlock)
                        
  end subroutine
  
end module
