!##############################################################################
!# ****************************************************************************
!# <name> nonlinearsolver </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains a black box nonlinear solver using a defect correction
!# loop. The solver basically performs the following iteration:
!#
!#     $$  x_{n+1}  =  x_n  +  J^{-1} ( b - A(x_n) x_n )  $$
!#
!# The defect correction loop is split into three tasks:
!#
!# 1.) Calculation of nonlinear defect
!#
!#         $$  d_n  :=  b - A(x_n)x_n  $$
!#
!#     The calculation of the nonlinear defect is performed user defined.
!#     The nonlinear loop itself does not know anything about the matrix.
!#     In fact, the matrix does not have to exist - this depends on the 
!#     problem.
!#
!# 2.) Preconditioning of the nonlinear defect
!#
!#         $$  u_n  :=  J^{-1} d_n  $$
!#
!#     with a preconditioner $J^{-1}$. This can be a Jacobian, a linear solver,
!#     a matrix, the inverse of a lumped mass matrix or similar.
!#
!# 3.) Correction
!#
!#         $$  x_{n+1}  :=  x_n  +  \omega u_n  $$
!#
!#     $\omega$ is usually = 1.
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

MODULE nonlinearsolver

  USE fsystem
  USE storage
  USE linearsolver
  USE spatialdiscretisation
  USE linearsystemscalar
  USE linearsystemblock
  USE collection
  
  IMPLICIT NONE

! *****************************************************************************
! *****************************************************************************
! *****************************************************************************

!<constants>

!<constantblock description="Identifiers for preconditioner type">
  
  ! User defined or no preconditioner (if callback routine is not given).
  INTEGER, PARAMETER :: NLSOL_PREC_USERDEF  = 0
  
  ! Preconditioner = Multiplication with a matrix
  INTEGER, PARAMETER :: NLSOL_PREC_MATRIX   = 1

  ! Preconditioner = Multiplication with inverse lumped mass matrix
  INTEGER, PARAMETER :: NLSOL_PREC_LMASS    = 2
  
  ! Preconditioner = Call to linear solver
  INTEGER, PARAMETER :: NLSOL_PREC_LINSOL   = 3
!</constantblock>

!<constantblock description="Identifiers for stopping criteria">

  ! Use standard stopping criterion.
  ! If depsRel>0: use relative stopping criterion.
  ! If depsAbs>0: use abs stopping criterion.
  ! If both are > 0: use both, i.e. stop if both criteria hold
  INTEGER, PARAMETER :: NLSOL_STOP_STANDARD = 0
  
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
  
  TYPE t_nlsolNode
    
    ! OUTPUT: Result
    ! The result of the solution process.
    ! =-1: convergence criterion not reached.
    ! =0: success. 
    ! =1: iteration broke down, diverging.
    ! =2: iteration broke down, preconditioner did not work.
    ! =3: error in the parameters.
    INTEGER                    :: iresult
    
    ! OUTPUT: Number of performed iterations, if the solver
    ! is of iterative nature.
    ! Is to 1 by the solver if not used (indicating at least 1 performed 
    ! iteration, which is always the case).
    INTEGER                    :: iiterations
    
    ! OUTPUT PARAMETER:
    ! Norm of initial residuum for each subvector of the given block vector.
    ! Only valid, if the fcb_resNormCheck callback procedure is not specified
    ! in the call to the solver, otherwise undefined!
    REAL(DP), DIMENSION(SPDISC_MAXEQUATIONS)  :: DinitialDefect

    ! OUTPUT PARAMETER:
    ! Norm of final residuum for each subvector of the given block vector.
    ! Only valid, if the fcb_resNormCheck callback procedure is not specified
    ! in the call to the solver, otherwise undefined!
    REAL(DP), DIMENSION(SPDISC_MAXEQUATIONS)  :: DfinalDefect

    ! OUTPUT PARAMETER:
    ! Total time for nonlinear solver
    REAL(DP)                        :: dtimeTotal

    ! OUTPUT PARAMETER:
    ! Total time for calculating nonlinear defects
    REAL(DP)                        :: dtimeNLdefect

    ! OUTPUT PARAMETER:
    ! Total time for nonlinear preconditioning
    REAL(DP)                        :: dtimeNLpreconditioning
    
    ! INPUT PARAMETER:
    ! Damping parameter for defect correction.
    ! Standard value = 1.0 (corresponds to 'no damping')
    REAL(DP)                        :: domega  = 1.0_DP

    ! INPUT PARAMETER:
    ! Standard relative stopping criterion for each subvector of a block vector.
    ! Only valid, if there is no callback routine for the calculation of
    ! residuals passed to the nonlinear solver.
    ! Stop iteration if everywhere
    !   !!defect!! < EPSREL * !!initial defect!!.
    ! =0: ignore, use absolute stopping criterion; standard = 1E-5
    ! Remark: don't set depsAbs=depsRel=0!
    REAL(DP), DIMENSION(SPDISC_MAXEQUATIONS) :: DepsRel = 1E-5_DP

    ! INPUT PARAMETER:
    ! Standard absolute stopping criterion for each subvector of a block vector. 
    ! Only valid, if there is no callback routine for the calculation of
    ! residuals passed to the nonlinear solver.
    ! Stop iteration if everywhere
    !   !!defect!! < EPSREL.
    ! =0: ignore, use relative stopping criterion; standard = 1E-5
    ! Remark: don't set depsAbs=depsRel=0!
    REAL(DP), DIMENSION(SPDISC_MAXEQUATIONS) :: DepsAbs = 1E-5_DP

    ! INPUT PARAMETER:
    ! Standard relative divergence criterion for each subvector of a block vector. 
    ! Only valid, if there is no callback routine for the calculation of
    ! residuals passed to the nonlinear solver.
    ! Treat iteration as diverged if anywhere
    !   !!defect!! >= DIVREL * !!initial defect!!
    ! A value of SYS_INFINITY disables the relative divergence check.
    ! standard = 1E3
    REAL(DP), DIMENSION(SPDISC_MAXEQUATIONS) :: DdivRel = 1E3_DP

    ! INPUT PARAMETER:
    ! Standard absolute divergence criterion for each subvector of a block vector. 
    ! Only valid, if there is no callback routine for the calculation of
    ! residuals passed to the nonlinear solver.
    ! Treat iteration as diverged if anywhere
    !   !!defect!! >= DIVREL
    ! A value of SYS_INFINITY disables the absolute divergence check.
    ! standard = SYS_INFINITY
    REAL(DP), DIMENSION(SPDISC_MAXEQUATIONS) :: DdivAbs = SYS_INFINITY

    ! INPUT PARAMETER: 
    ! Type of stopping criterion to use. One of the
    ! NLSOL_STOP_xxxx constants.
    INTEGER                    :: istoppingCriterion = NLSOL_STOP_STANDARD

    ! INPUT PARAMETER: 
    ! Minimum number of iterations top perform
    INTEGER                    :: nminIterations = 1

    ! INPUT PARAMETER: 
    ! Maximum number of iterations top perform
    INTEGER                    :: nmaxIterations = 50
    
    ! INPUT PARAMETER FOR SOLVERS WITH RESIDUAL CHECK: 
    ! For every subvector: Type of norm to use in the residual checking 
    ! (cf. linearalgebra.f90).
    ! =0: euclidian norm, =1: l1-norm, =2: l2-norm, =3: MAX-norm
    INTEGER, DIMENSION(SPDISC_MAXEQUATIONS) :: IresNorm = 2
    
    ! INPUT PARAMETER: Output level
    ! This determines the output level of the solver.
    ! =0: no output, =1: basic output, =2, extended output
    INTEGER                    :: ioutputLevel = 0

    ! INPUT PARAMETER:
    ! Type of preconditioner to use. One of the NLSOL_PREC_xxxx constants.
    INTEGER                    :: cpreconditioner = NLSOL_PREC_USERDEF
    
    ! INPUT PARAMETER: Output save.
    ! If = 0, autosave is deactivated.
    ! If > 0, the current iteration vector is saved every iautosave iterations
    !  to disc. The basic filename is given by sautosaveName.
    INTEGER                    :: iautosave    = 0
    
    ! INPUT PARAMETER
    ! Basic filename if autosave is activated (iautosave > 0). The number of the
    ! current iteration is appended to this filename, so the filename
    ! for writing a solution in iteration ite will be '[sautosaveName].[ite]'.
    CHARACTER(LEN=SYS_STRLEN)  :: sautosaveName = 'solution'
    
    ! INPUT PARAMETER: 
    ! If cpreconditioner=NLSOL_PREC_MATRIX block matrix that is multiplied 
    ! to the defect vector.
    ! If cpreconditioner=NLSOL_PREC_LMASS: Lumped mass matrix. Only matrix
    ! format 7 and 9 are supported.
    TYPE(t_matrixBlock)        :: rmatrixPreconditioner
    
    ! INPUT PARAMETER:
    ! If cpreconditioner=NLSOL_PREC_LINSOL: Pointer to linear solver node
    ! defining the preconditioner
    TYPE(t_linsolNode), POINTER :: p_rlinsolNode => NULL()
    
    ! STATUS: 
    ! Current iteration
    INTEGER                     :: icurrentIteration

  END TYPE
  
!</typeblock>
  
!</types>

! *****************************************************************************
! *****************************************************************************
! *****************************************************************************

CONTAINS

  ! ***************************************************************************

!<function>
  
  LOGICAL FUNCTION nlsol_testConvergence (rsolverNode, DvecNorm, nblocks, rdef) &
          RESULT(loutput)
  
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
  TYPE(t_nlsolNode), INTENT(IN) :: rsolverNode
  
  ! Number of blocks in the solution vector/equation
  INTEGER, INTENT(IN) :: nblocks
  
  ! OPTIONAL: The defect vector which norm should be tested.
  ! If existent, the norm of the subvectors is returned in DvecNorm.
  ! If not existent, the routine assumes that DvecNrm is the norm
  ! of the vector and checks convergence depending on DvecNorm.
  TYPE(t_vectorBlock), INTENT(IN), OPTIONAL :: rdef
  
!</input>

!<inputoutput>
  
  ! Norm of the subvectors on the defect vector.
  ! If rdef if present, the routine will calculate the norm of the subvectors
  ! and return them in dvecNorm.
  ! If rdef is not present, DvecNorm is assumed to be a valid norm of a
  ! vector and convergence is tested using DvecNorm.
  REAL(DP), DIMENSION(:), INTENT(INOUT) :: DvecNorm

!</inputoutput>
  
!</function>

  ! local variables
  INTEGER :: i

  ! Calculate the norm of the vector or take the one given
  ! as parameter
  IF (PRESENT(rdef)) THEN
    DvecNorm = lsysbl_vectorNormBlock (rdef,rsolverNode%IresNorm)
    WHERE (.NOT.((DvecNorm .GE. 1D-99) .AND. (DvecNorm .LE. 1D99))) 
      DvecNorm = 0.0_DP
    END WHERE
  END IF
  
  loutput = .TRUE.
  
  DO i=1,nblocks
    !  Absolute convergence criterion? Check the norm directly.
    IF (rsolverNode%DepsAbs(i) .NE. 0.0_DP) THEN
      IF (DvecNorm(i) .GT. rsolverNode%DepsAbs(i)) THEN
        loutput = .FALSE.
        RETURN
      END IF
    END IF
    
    ! Relative convergence criterion? Multiply with initial residuum
    ! and check the norm. 
    ! If the initial defect is 0, this rule is not to be applied!
    ! (May happen in the first step of 'flow around cylinder'
    !  where we have only convection in the X-direction, not in the Y-direction
    !  and a still fluid.)
    IF ((rsolverNode%DepsRel(i) .NE. 0.0_DP) .AND. &
        (rsolverNode%dinitialDefect(i) .GT. SYS_EPSREAL)) THEN
      IF (DvecNorm(i) .GT. &
          rsolverNode%depsRel(i) * rsolverNode%dinitialDefect(i)) THEN
        loutput = .FALSE.
        RETURN
      END IF
    END IF
    
  END DO ! i
  
  END FUNCTION
  
  ! ***************************************************************************

!<function>
  
  LOGICAL FUNCTION nlsol_testDivergence (rsolverNode, DvecNorm, nblocks, rdef) &
          RESULT(loutput)
  
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
  TYPE(t_nlsolNode), INTENT(IN) :: rsolverNode
  
  ! Number of blocks in the solution vector/equation
  INTEGER, INTENT(IN) :: nblocks
  
  ! OPTIONAL: The defect vector which norm should be tested.
  ! If existent, the norm of the subvectors is returned in DvecNorm.
  ! If not existent, the routine assumes that DvecNrm is the norm
  ! of the vector and checks convergence depending on DvecNorm.
  TYPE(t_vectorBlock), INTENT(IN), OPTIONAL :: rdef
  
!</input>

!<inputoutput>
  
  ! Norm of the subvectors on the defect vector.
  ! If rdef if present, the routine will calculate the norm of the subvectors
  ! and return them in dvecNorm.
  ! If rdef is not present, DvecNorm is assumed to be a valid norm of a
  ! vector and convergence is tested using DvecNorm.
  REAL(DP), DIMENSION(:), INTENT(INOUT) :: DvecNorm

!</inputoutput>
  
!</function>

  ! local variables
  INTEGER :: i

  ! Calculate the norm of the vector if not given
  ! as parameter
  IF (PRESENT(rdef)) THEN
    DvecNorm = lsysbl_vectorNormBlock (rdef,rsolverNode%IresNorm)
    WHERE (.NOT.((DvecNorm .GE. 1D-99) .AND. (DvecNorm .LE. 1D99))) 
      DvecNorm = 0.0_DP
    END WHERE
  END IF
  
  loutput = .FALSE.
  
  DO i=1,nblocks
  
    ! Absolute divergence criterion? Check the norm directly.
    IF (rsolverNode%DdivAbs(i) .NE. SYS_INFINITY) THEN
     
      ! use NOT here - gives a better handling of special cases like NaN!
      IF ( .NOT. (DvecNorm(i) .LE. rsolverNode%DdivAbs(i))) THEN
        loutput = .TRUE.
        RETURN
      END IF
      
    END IF
    
    ! Relative divergence criterion? Multiply with initial residuum
    ! and check the norm. 
    ! If the initial defect is 0, this rule is not to be applied!
    ! (May happen in the first step of 'flow around cylinder'
    !  where we have only convection in the X-direction, not in the Y-direction
    !  and a still fluid.)
    IF ((rsolverNode%DinitialDefect(i) .NE. 0.0_DP) .AND. &
        (rsolverNode%DepsRel(i) .NE. SYS_INFINITY)) THEN
      IF ( .NOT. (DvecNorm(i) .LE. &
          rsolverNode%DinitialDefect(i) * rsolverNode%DdivRel(i)) ) THEN
        loutput = .TRUE.
        RETURN
      END IF
    END IF
  
  END DO ! i
  
  END FUNCTION

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE nlsol_setPrecMatrix (rsolverNode,rmatrix)
  
!<description>
  ! This routine installs a standard matrix preconditioner into the nonlinear
  ! subroutine. rmatrix =: P is a block matrix which is multiplied to the defect
  ! vector in every iteration:
  !    $$  x_{n+1}  =  x_n  +  P (b-A(x)x)  $$
!</description>

!<input>
  ! The block matrix.
  TYPE(t_matrixBlock), INTENT(IN) :: rmatrix
!</input>

!<inputoutput>
  ! The solver node of the nonlinear iteration loop where to set the 
  ! preconditioner.
  TYPE(t_nlsolNode), INTENT(INOUT) :: rsolverNode
!</inputoutput>

!</subroutine>

    ! Install the matrix. Copy the matrix structure into the solver node.
    rsolverNode%rmatrixPreconditioner = rmatrix
    rsolverNode%rmatrixPreconditioner%RmatrixBlock%imatrixSpec = &
      IOR(rsolverNode%rmatrixPreconditioner%RmatrixBlock%imatrixSpec,LSYSSC_MSPEC_ISCOPY)
      
    ! Set the preconditioner type flag appropriately.
    rsolverNode%cpreconditioner = NLSOL_PREC_MATRIX

  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE nlsol_setPrecMatrixSc (rsolverNode,rmatrix)
  
!<description>
  ! This routine installs a standard matrix preconditioner into the nonlinear
  ! subroutine. rmatrix is a scalar matrix. The matrix is added as 1x1 block
  ! matrix preconditioner P := [rmatrix] to the solver node and is multiplied 
  ! to the defect vector in every iteration:
  !    $$  x_{n+1}  =  x_n  +  P (b-A(x)x)  $$
!</description>

!<input>
  ! The scalar matrix.
  TYPE(t_matrixScalar), INTENT(IN) :: rmatrix
!</input>

!<inputoutput>
  ! The solver node of the nonlinear iteration loop where to set the 
  ! preconditioner.
  TYPE(t_nlsolNode), INTENT(INOUT) :: rsolverNode
!</inputoutput>

!</subroutine>

    ! Create a 1x1 block matrix from rmatrix and save it to the solver node.
    CALL lsysbl_createMatFromScalar (rmatrix,rsolverNode%rmatrixPreconditioner)
    rsolverNode%rmatrixPreconditioner%RmatrixBlock%imatrixSpec = &
      IOR(rsolverNode%rmatrixPreconditioner%RmatrixBlock%imatrixSpec,LSYSSC_MSPEC_ISCOPY)
      
    ! Set the preconditioner type flag appropriately.
    rsolverNode%cpreconditioner = NLSOL_PREC_MATRIX

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE nlsol_setPrecMass (rsolverNode,rmatrix)
  
!<description>
  ! This routine installs a lumped mass matrix preconditioner into the nonlinear
  ! subroutine. rmatrix =: M is a block matrix which contains only diagonal 
  ! blocks with lumped mass matrices. this matrix is is multiplied to the defect
  ! vector in every iteration:
  !    $$  x_{n+1}  =  x_n  +  M^{-1} (b-A(x)x)  $$
!</description>

!<input>
  ! The block matrix. Must be a matrix only containing diagonal blocks.
  ! the diagonal blocks must be lumped mass matrices.
  TYPE(t_matrixBlock), INTENT(IN) :: rmatrix
!</input>

!<inputoutput>
  ! The solver node of the nonlinear iteration loop where to set the 
  ! preconditioner.
  TYPE(t_nlsolNode), INTENT(INOUT) :: rsolverNode
!</inputoutput>

!</subroutine>

    ! Install the matrix. 
    CALL nlsol_setPrecMatrix (rsolverNode,rmatrix)
      
    ! Set the preconditioner type flag appropriately.
    rsolverNode%cpreconditioner = NLSOL_PREC_LMASS

  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE nlsol_setPrecMassSc (rsolverNode,rmatrix)
  
!<description>
  ! This routine installs a standard matrix preconditioner into the nonlinear
  ! subroutine. rmatrix is a scalar matrix. The matrix is added as 1x1 block
  ! matrix preconditioner P := [rmatrix] to the solver node and its inverse
  ! is multiplied to the defect vector in every iteration:
  !    $$  x_{n+1}  =  x_n  +  M^{-1} (b-A(x)x)  $$
!</description>

!<input>
  ! The lumped mass matrix.
  TYPE(t_matrixScalar), INTENT(IN) :: rmatrix
!</input>

!<inputoutput>
  ! The solver node of the nonlinear iteration loop where to set the 
  ! preconditioner.
  TYPE(t_nlsolNode), INTENT(INOUT) :: rsolverNode
!</inputoutput>

!</subroutine>

    ! Install the matrix
    CALL nlsol_setPrecMatrixSc (rsolverNode,rmatrix)
      
    ! Set the preconditioner type flag appropriately.
    rsolverNode%cpreconditioner = NLSOL_PREC_LMASS

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE nlsol_setPrecLinsol (rsolverNode,rlinsolNode)
  
!<description>
  ! This routine installs a linear solver node as preconditioner to the
  ! nonlinear solver. The linear solver $J^{-1}$ is called after each
  ! defect calculation to precondition the defect:
  !    $$  x_{n+1}  =  x_n  +  J^{-1} (b-A(x)x)  $$
!</description>

!<input>
  ! The linear solver node identifying the preconditioner
  TYPE(t_linsolNode), INTENT(IN), TARGET :: rlinsolNode
!</input>

!<inputoutput>
  ! The solver node of the nonlinear iteration loop where to set the 
  ! preconditioner.
  TYPE(t_nlsolNode), INTENT(INOUT) :: rsolverNode
!</inputoutput>

!</subroutine>

    ! Install the preconditioner
    rsolverNode%p_rlinsolNode => rlinsolNode
      
    ! Set the preconditioner type flag appropriately.
    rsolverNode%cpreconditioner = NLSOL_PREC_LINSOL

  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE nlsol_performSolve (rsolverNode,rx,rb,rd,&
                                 fcb_getDefect,fcb_precondDefect,fcb_resNormCheck,&
                                 rcollection)
             
!<description>
  ! This routine invokes the nonlinear defect correction iteration
  !     $$  x_{n+1}  =  x_n  +  J^{-1} ( b - A(x_n) x_n )  $$
  !
  ! The defect correction loop is split into three tasks:
  !
  ! 1.) Calculation of nonlinear defect
  !         $$  d_n  :=  b - A(x_n)x_n  $$
  !     For this purpose, the callback routine fcb_getDefect is called.
  !
  ! 2.) Preconditioning of the nonlinear defect
  !         $$  u_n  :=  J^{-1} d_n  $$
  !     with a preconditioner $J^{-1}$. The preconditioner can be
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
  !         $$  x_{n+1}  :=  x_n  +  \omega u_n  $$
  !     $\omega$ is usually = 1. A user defined callback routine fcb_precondDefect
  !     can modify $\omega$ in every iteration.
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
  TYPE(t_nlsolNode), INTENT(INOUT)              :: rsolverNode
  
  ! INPUT: Initial solution vector.
  ! OUTPUT: Final iteration vector.
  TYPE(t_vectorBlock), INTENT(INOUT)            :: rx
             
  ! Temporary vector. Must be of the same size/type as rx/rb.
  TYPE(t_vectorBlock), INTENT(INOUT)            :: rd
  
  ! OPTIONAL: Collection structure that saves problem-dependent information.
  ! This is passed without being changed to the callback routines of this
  ! algorithm.
  TYPE(t_collection), INTENT(INOUT), TARGET, OPTIONAL :: rcollection
  
!</inputoutput>
  
!<input>
  ! Right hand side vector of the equation.
  TYPE(t_vectorBlock), INTENT(IN)               :: rb
  
  ! Defect vector calculation routine. Based on the current iteration vector
  ! rx and the right hand side vector rb, this routine has to compute the 
  ! defect vector rd.
  INCLUDE 'intf_nlsolcallback.inc'
  
  ! OPTIONAL: Preconditioning routine. This routine accepts a defect vector rd
  ! and replaces it by a preconditioned defect vector $J^{-1} rd$. 
  ! If this parameter is not present, the preconditioner is either a matrix
  ! or there is no preconditioner (depending on the variable
  ! rsolverNode\%cpreconditioner).
  OPTIONAL :: fcb_precondDefect
  
  ! OPTIONAL: Residual norm calculation and printing routine.
  ! If not present, the standard absolute/relative stopping criteria of the
  !  solver node with the configuration of the nonlinear solver are used.
  ! If present, this routine is called each time the norm of the residuum was 
  !  calculated. It has to check the current residuum for convergence and/or 
  !  divergence and can print the residuum to screen.
  OPTIONAL :: fcb_resNormCheck
  
!</input>

!</subroutine>

  ! local variables
  TYPE(t_collection), POINTER :: p_rcollection
  INTEGER :: ite,i
  REAL(DP), DIMENSION(SPDISC_MAXEQUATIONS) :: DvecNorm
  TYPE(t_vectorBlock) :: rtemp
  REAL(DP) :: domega
  LOGICAL :: bconvergence,bdivergence,bsuccess
  
    ! Do we have a collection?
    NULLIFY(p_rcollection)
    IF (PRESENT(rcollection)) p_rcollection => rcollection
    
    ! In case our preconditioner is a matrix-vector multiplication,
    ! allocate memory for another temporary vector used 
    ! during the MV multiplication.
    IF (rsolverNode%cpreconditioner .EQ. NLSOL_PREC_MATRIX) THEN
      CALL lsysbl_createVecBlockIndirect (rx,rtemp,.FALSE.)
    END IF
    
    ! Status reset
    rsolverNode%iresult = 0
    
    ! Calculate the initial nonlinear defect to rd:   d = b-A(x)x
    CALL fcb_getDefect (0,rx,rb,rd,p_rcollection)
    
    ite = 0
    rsolverNode%icurrentIteration = ite
    
    ! Initial test for convergence/divergence.
    IF (PRESENT(fcb_resNormCheck)) THEN
      ! Calculate the defect with the callback routine
      CALL fcb_resNormCheck (ite,rx,rb,rd,bconvergence,bdivergence,p_rcollection)
    ELSE
      ! Calculate the norm of the defect:
      DvecNorm = 0.0_DP
      DvecNorm = lsysbl_vectorNormBlock (rd,rsolverNode%IresNorm(1:rd%nblocks))
      WHERE (.NOT.((DvecNorm .GE. 1D-99) .AND. (DvecNorm .LE. 1D99))) 
        DvecNorm = 0.0_DP
      END WHERE
      rsolverNode%DinitialDefect = DvecNorm
      rsolverNode%DfinalDefect = DvecNorm
      
      bconvergence = nlsol_testConvergence (rsolverNode, DvecNorm, rd%nblocks)
      bdivergence  = nlsol_testDivergence (rsolverNode, DvecNorm, rd%nblocks)
      
      IF (rsolverNode%ioutputLevel .GE. 2) THEN
        CALL output_line ('NLSOL: Iteration '//&
             TRIM(sys_siL(ite,10))//', !!RES!! =',bnolinebreak=.TRUE.)
        DO i=1,rb%nblocks
          CALL output_line (' '//TRIM(sys_sdEL(DvecNorm(i),15)),bnolinebreak=.TRUE.)
        END DO
        CALL output_lbrk()
      END IF
    END IF

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
        SELECT CASE (rsolverNode%cpreconditioner)
        CASE (NLSOL_PREC_MATRIX)
          ! Multiplication with a matrix.
          CALL lsysbl_copyVector (rd,rtemp)
          CALL lsysbl_blockMatVec (rsolverNode%rmatrixPreconditioner, &
                                   rtemp, rd, 1.0_DP, 0.0_DP)
                                   
        CASE (NLSOL_PREC_LMASS)
          ! Multiply with the inverse of the diagonal in the block matrix.
          CALL lsysbl_invertedDiagMatVec (rsolverNode%rmatrixPreconditioner,&
                                          rd,1.0_DP,rd)
                                          
        CASE (NLSOL_PREC_LINSOL)
          ! Preconditioner is a linear solver with fixed matrix.
          CALL linsol_precondDefect(rsolverNode%p_rlinsolNode,rd)
          bsuccess = rsolverNode%p_rlinsolNode%iresult .EQ. 0
          
        CASE DEFAULT
          ! User defined or no preconditioning.
          ! Is a callback function available?
          IF (PRESENT(fcb_precondDefect)) THEN
            ! The callback routine is allowed to change domega during the
            ! iteration if necessary. The nonlinear solver here does not touch
            ! domega anymore, so the callback routine is the only one changing it.
            CALL fcb_precondDefect (ite,rd,rx,rb,domega,bsuccess,p_rcollection)
          END IF
          
        END SELECT
        
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
          
          ! Autosave: Save the current iteration vector to disc if desired.
          IF (rsolverNode%iautosave .NE. 0) THEN
            ! Don't change the order of these if's - produces DIVISION BY ZERO on
            ! some compilers with optimisation!
            IF (MOD(ite,rsolverNode%iautosave) .EQ. 0) THEN
              CALL output_line ('Iteration: '//&
                  TRIM(sys_siL(ITE,10))//': Autosave not yet implemented.')
            END IF
          END IF

          ! Calculate the new nonlinear defect to rd:  d = b-A(x)x
          CALL fcb_getDefect (ite,rx,rb,rd,p_rcollection)

          ! Check the defect for convergence        
          IF (PRESENT(fcb_resNormCheck)) THEN
            ! Calculate the defect with the callback routine
            CALL fcb_resNormCheck (ite,rx,rb,rd,bconvergence,bdivergence,p_rcollection)
          ELSE
            ! Calculate the norm of the defect:
            DvecNorm(1:rd%nblocks) = lsysbl_vectorNormBlock (rd,rsolverNode%IresNorm)
            WHERE (.NOT.((DvecNorm .GE. 1E-99_DP) .AND. (DvecNorm .LE. 1E99_DP))) 
              DvecNorm = 0.0_DP
            END WHERE
            rsolverNode%DfinalDefect = DvecNorm
            
            bconvergence = nlsol_testConvergence (rsolverNode, DvecNorm, rd%nblocks)
            bdivergence  = nlsol_testDivergence (rsolverNode, DvecNorm, rd%nblocks)

            IF (rsolverNode%ioutputLevel .GE. 2) THEN
              CALL output_line ('NLSOL: Iteration '//&
                  TRIM(sys_siL(ite,10))//', !!RES!! =',bnolinebreak=.TRUE.)
              DO i=1,rb%nblocks
                CALL output_line (' '//TRIM(sys_sdEL(DvecNorm(i),15)),bnolinebreak=.TRUE.)
              END DO
              CALL output_lbrk()
            END IF
          END IF

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

  SUBROUTINE nlsol_performSolveSc (rsolverNode,rx,rb,rd,&
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
  TYPE(t_nlsolNode), INTENT(INOUT)              :: rsolverNode
  
  ! INPUT: Initial solution vector.
  ! OUTPUT: Final iteration vector.
  TYPE(t_vectorScalar), INTENT(INOUT)            :: rx
             
  ! Temporary vector. Must be of the same size/type as rx/rb.
  TYPE(t_vectorScalar), INTENT(INOUT)            :: rd
  
  ! OPTIONAL: Collection structure that saves problem-dependent information.
  ! This is passed without being changed to the callback routines of this
  ! algorithm.
  TYPE(t_collection), INTENT(INOUT), TARGET, OPTIONAL :: rcollection
  
!</inputoutput>
  
!<input>
  ! Right hand side vector of the equation.
  TYPE(t_vectorScalar), INTENT(IN)               :: rb
  
  ! Defect vector calculation routine. Based on the current iteration vector
  ! rx and the right hand side vector rb, this routine has to compute the 
  ! defect vector rd.
  INCLUDE 'intf_nlsolcallback.inc'
  
  ! OPTIONAL: Preconditioning routine. This routine accepts a defect vector rd
  ! and replaces it by a preconditioned defect vector $J^{-1} rd$. 
  ! If this parameter is not present, the preconditioner is either a matrix
  ! or there is no preconditioner (depending on the variable
  ! rsolverNode\%cpreconditioner).
  OPTIONAL :: fcb_precondDefect
  
  ! OPTIONAL: Residual norm calculation and printing routine.
  ! If not present, the standard absolute/relative stopping criteria of the
  !  solver node with the configuration of the nonlinear solver are used.
  ! If present, this routine is called each time the norm of the residuum was 
  !  calculated. It has to check the current residuum for convergence and/or 
  !  divergence and can print the residuum to screen.
  OPTIONAL :: fcb_resNormCheck
  
!</input>

!</subroutine>

  ! local variables
  TYPE(t_vectorBlock) :: rxBlock,rbBlock,rdBlock
  
  ! Convert the vectors on-the-fly to block vectors.
  ! The new vectors share the same memory as the old, so the solver will use
  ! and overwrite the old input vectors.
  CALL lsysbl_createVecFromScalar (rx,rxBlock)
  CALL lsysbl_createVecFromScalar (rb,rbBlock)
  CALL lsysbl_createVecFromScalar (rd,rdBlock)
  
  ! Invoke the solver - that's all. 
  CALL nlsol_performSolve (rsolverNode,rxBlock,rbBlock,rdBlock,&
                          fcb_getDefect,fcb_precondDefect,fcb_resNormCheck,&
                          rcollection)
                        
  END SUBROUTINE
  
END MODULE
