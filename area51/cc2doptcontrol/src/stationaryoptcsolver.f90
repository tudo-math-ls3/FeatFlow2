!##############################################################################
!# ****************************************************************************
!# <name> stationaryoptcsolver </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module invokes the nonlinear solver to solve the basic CC2D
!# optimal control problem.
!#
!# 1.) cc_solve
!#     -> Invoke the nonlinear solver to solve the stationary coupled system.
!#
!# </purpose>
!##############################################################################

module stationaryoptcsolver

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
  
  use collection
  use convection
    
  use basicstructures
  use user_callback
  
  use spacematvecassembly
  use spacepreconditioner
  use spacepreconditionerinit
  use postprocessing
  
  implicit none

contains
  
  ! ***************************************************************************

!<subroutine>

  subroutine cc_solve (rproblem,rvector,rrhs)
  
!<description>
  ! Solves the given problem by applying a nonlinear solver with a preconditioner
  ! configured in the INI/DAT files.
  !
  ! On call to this routine, rproblem%rrhs configures the right hand side of
  ! the nonlinear iteration, while rproblem%rvector specifies the initial iteration
  ! vector.
  ! With $u$:=rvector, $f$:=rrhs and the (nonlinear) matrix $A(u)$ configured
  ! in rproblem, the routine calls the nonlinear solver to solve $A(u)u = f$.
  ! During the nonlinear iteration, the nonlinear matrix (matrices on all levels),
  ! and the vector rvector in rproblem will be changed.
  ! rproblem%rvector receives the final solution vector of the nonlinear iteration.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(INOUT), target :: rproblem

  ! The solution vector which is to be used as initial vector for the nonlinear
  ! iteration. After the iteration, this is replaced by the new solution vector.
  type(t_vectorBlock), intent(INOUT) :: rvector
!</inputoutput>

!<input>
  ! The right-hand-side vector to use in the equation
  type(t_vectorBlock), intent(IN) :: rrhs
!</input>

!</subroutine>

!    ! local variables
!    REAL(DP) :: dalphaC
!    TYPE(t_nonlinearSpatialMatrix) :: rmatrix
!    TYPE(t_ccspatialPreconditioner) :: rpreconditioner
!    TYPE(t_vectorBlock) :: rd
!    LOGICAL :: bsuccess
!
!    ! The nonlinear solver configuration
!    TYPE(t_nlsolNode) :: rnlSol
!
!    ! Some parameters for the nonlinear loop
!    INTEGER :: nminIterations,nmaxIterations,iglobIter
!    REAL(DP) :: depsRel,depsAbs,domega,ddefNorm,dinitDefNorm
!
!    !!!!!!!!!!!!!!
!    ! NOT TESTED !
!    !!!!!!!!!!!!!!
!
!    ! Set up all the weights in the core equation according to the current timestep.
!    CALL parlst_getvalue_double (rproblem%rparamList,'OPTIMALCONTROL',&
!                                'dalphaC',dalphaC,0.1_DP)
!    rmatrix%dnu = rproblem%dnu
!
!    rmatrix%diota1  = 0.0_DP
!    rmatrix%dkappa1 = 0.0_DP
!    rmatrix%dalpha1 = 0.0_DP
!    rmatrix%dtheta1 = 1.0_DP
!    rmatrix%dgamma1 = REAL(1-rproblem%iequation,DP)
!    rmatrix%deta1   = 1.0_DP
!    rmatrix%dtau1   = 1.0_DP
!    rmatrix%dmu1    = 1.0_DP/dalphaC
!
!    rmatrix%diota2  = 0.0_DP
!    rmatrix%dkappa2 = 0.0_DP
!    rmatrix%dalpha2 = 0.0_DP
!    rmatrix%dtheta2 = 1.0_DP
!    rmatrix%dgamma2 = REAL(1-rproblem%iequation,DP)
!    rmatrix%deta2   = 1.0_DP
!    rmatrix%dtau2   = 1.0_DP
!    rmatrix%dmu2    = -1.0_DP
!
!    ! Initialise pointers to the matrix building blocks
!    rmatrix%p_rpreallocatedMatrix => &
!      rproblem%RlevelInfo(rproblem%NLMAX)%rpreallocatedSystemMatrix
!
!    rmatrix%p_rdiscretisation => &
!      rproblem%RlevelInfo(rproblem%NLMAX)%p_rdiscretisation
!
!    rmatrix%p_rmatrixTemplateFEM => &
!      rproblem%RlevelInfo(rproblem%NLMAX)%rmatrixTemplateFEM
!
!    rmatrix%p_rmatrixTemplateGradient => &
!      rproblem%RlevelInfo(rproblem%NLMAX)%rmatrixTemplateGradient
!
!    rmatrix%p_rmatrixStokes => &
!      rproblem%RlevelInfo(rproblem%NLMAX)%rmatrixStokes
!
!    rmatrix%p_rmatrixB1 => &
!      rproblem%RlevelInfo(rproblem%NLMAX)%rmatrixB1
!
!    rmatrix%p_rmatrixB2 => &
!      rproblem%RlevelInfo(rproblem%NLMAX)%rmatrixB2
!
!    rmatrix%p_rmatrixMass => &
!      rproblem%RlevelInfo(rproblem%NLMAX)%rmatrixMass
!
!    rmatrix%p_rmatrixIdentityPressure => &
!      rproblem%RlevelInfo(rproblem%NLMAX)%rmatrixIdentityPressure
!
!
!    ! Initialise the preconditioner for the nonlinear iteration
!    CALL cc_initPreconditioner (rproblem,rproblem%NLMIN,rproblem%NLMAX,&
!        rpreconditioner,0)
!    CALL cc_configPreconditioner (rproblem,rpreconditioner)
!
!    ! Print out the value of the optimal control functional for the
!    ! initial solution vector directly prior to the solution process.
!    CALL cc_printControlFunctionalStat (rproblem,rvector)
!
!    ! Get configuration parameters from the DAT file
!    CALL parlst_getvalue_int (rproblem%rparamList, 'CC2D-NONLINEAR', &
!                              'nminIterations', nminIterations, 1)
!    CALL parlst_getvalue_int (rproblem%rparamList, 'CC2D-NONLINEAR', &
!                             'nmaxIterations', nmaxIterations, 10)
!    CALL parlst_getvalue_double (rproblem%rparamList, 'CC2D-NONLINEAR', &
!                                 'depsRel', depsRel, 1E-5_DP)
!    CALL parlst_getvalue_double (rproblem%rparamList, 'CC2D-NONLINEAR', &
!                                 'depsAbs', depsAbs, 1E-5_DP)
!    CALL parlst_getvalue_double (rproblem%rparamList, 'CC2D-NONLINEAR', &
!                                 'domega', domega, 1.0_DP)
!
!    ! allocate a temp vector
!    CALL lsysbl_createVecBlockIndirect (rvector,rd)
!
!    ! Get the initial nonlinear defect
!    CALL lsysbl_copyVector (rrhs,rd)
!    CALL cc_assembleDefect (rmatrix,rvector,rd)
!
!    ! ... and its norm
!    ddefNorm = lsysbl_vectorNorm (rd,LINALG_NORML2)
!    dinitDefNorm = ddefNorm
!
!    ! Nonlinear loop
!    iglobIter = 0
!
!    CALL output_separator (OU_SEP_EQUAL)
!    CALL output_line ('Iteration: '//sys_si(iglobIter,10)//&
!        ' ||Res|| = '//sys_sdEP(ddefNorm,20,10))
!
!    DO WHILE ((iglobIter .LT. nminIterations) .OR. &
!              ((ddefNorm .GT. depsRel*dinitDefNorm) .AND. (ddefNorm .GT. depsAbs) &
!              .AND. (iglobIter .LT. nmaxIterations)))
!
!      iglobIter = iglobIter+1
!
!      ! Preconditioning of the defect: d=C^{-1}d
!      CALL cc_precondDefect (rpreconditioner,rmatrix,rd,rvector,bsuccess,&
!          rproblem%rcollection)
!
!      ! Add the defect: x = x + omega*d
!      CALL lsysbl_vectorLinearComb (rd,rvector,domega,1.0_DP)
!
!      ! Assemble the new defect: d=b-Ax
!      CALL lsysbl_copyVector (rrhs,rd)
!      CALL cc_assembleDefect (rmatrix,rvector,rd)
!
!      ! ... and its norm
!      ddefNorm = lsysbl_vectorNorm (rd,LINALG_NORML2)
!
!      CALL output_line ('Iteration: '//sys_si(iglobIter,10)//&
!          ' ||Res|| = '//sys_sdEP(ddefNorm,20,10))
!
!    END DO
!
!    CALL output_separator (OU_SEP_EQUAL)
!
!    ! Release the preconditioner and the temp vector
!    CALL cc_releasePreconditioner (rpreconditioner)
!
!    CALL lsysbl_releaseVector (rd)
!
!    CALL output_lbrk()
!    CALL output_line ('Nonlinear solver statistics')
!    CALL output_line ('---------------------------')
!    CALL output_line ('Initial defect: '//TRIM(sys_sdEL(rnlSol%DinitialDefect(1),15)))
!    CALL output_line ('Final defect:  '//TRIM(sys_sdEL(rnlSol%DfinalDefect(1),15)))
!    CALL output_line ('#Iterations:   '//TRIM(sys_siL(rnlSol%iiterations,10)))
!
  end subroutine

end module
