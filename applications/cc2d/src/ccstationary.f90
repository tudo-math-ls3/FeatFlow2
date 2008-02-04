!##############################################################################
!# ****************************************************************************
!# <name> ccstationary </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module invokes the nonlinear solver to solve the basic CC2D problem.
!#
!# 1.) cc_solve
!#     -> Invoke the nonlinear solver to solve the stationary coupled system.
!#
!# </purpose>
!##############################################################################

MODULE ccstationary

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
  
  USE collection
  USE convection
    
  USE ccbasic
  USE cccallback
  
  USE ccnonlinearcore
  USE ccnonlinearcoreinit
  
  IMPLICIT NONE

CONTAINS
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE cc_solve (rproblem,rvector,rrhs)
  
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
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem

  ! The solution vector which is to be used as initial vector for the nonlinear
  ! iteration. After the iteration, this is replaced by the new solution vector.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rvector
!</inputoutput>

!<input>
  ! The right-hand-side vector to use in the equation
  TYPE(t_vectorBlock), INTENT(IN) :: rrhs
!</input>

!</subroutine>

    ! local variables
    TYPE(t_ccNonlinearIteration) :: rnonlinearIteration

    ! The nonlinear solver configuration
    TYPE(t_nlsolNode) :: rnlSol
    
    ! Initialise the nonlinear solver node rnlSol with parameters from
    ! the INI/DAT files.
    CALL cc_getNonlinearSolver (rnlSol, rproblem%rparamList, 'CC2D-NONLINEAR')

    ! Initialise the nonlinear loop. This is to prepare everything for
    ! or callback routines that are called from the nonlinear solver.
    ! The preconditioner in that structure is initialised later.
    CALL cc_initNonlinearLoop (&
        rproblem,rproblem%NLMIN,rproblem%NLMAX,rvector,rrhs,&
        rnonlinearIteration,'CC2D-NONLINEAR')
        
    ! Set up all the weights in the core equation according to the current timestep.
    rnonlinearIteration%dalpha = 0.0_DP
    rnonlinearIteration%dtheta = 1.0_DP
    rnonlinearIteration%dgamma = REAL(1-rproblem%iequation,DP)
    rnonlinearIteration%deta   = 1.0_DP
    rnonlinearIteration%dtau   = 1.0_DP

    ! Initialise the preconditioner for the nonlinear iteration
    CALL cc_initPreconditioner (rproblem,&
        rnonlinearIteration,rvector,rrhs)

    ! Call the nonlinear solver to solve the core equation.
    CALL cc_solveCoreEquation (rproblem,rnonlinearIteration,rnlSol,&
        rvector,rrhs)             
             
    ! Release the preconditioner
    CALL cc_releasePreconditioner (rnonlinearIteration)
    
    ! Release parameters of the nonlinear loop, final clean up
    CALL cc_doneNonlinearLoop (rnonlinearIteration)
             
    CALL output_lbrk()
    CALL output_line ('Nonlinear solver statistics')
    CALL output_line ('---------------------------')
    CALL output_line ('Intial defect: '//TRIM(sys_sdEL(rnlSol%DinitialDefect(1),15)))
    CALL output_line ('Final defect:  '//TRIM(sys_sdEL(rnlSol%DfinalDefect(1),15)))
    CALL output_line ('#Iterations:   '//TRIM(sys_siL(rnlSol%iiterations,10)))
    
  END SUBROUTINE

END MODULE
