!##############################################################################
!# ****************************************************************************
!# <name> cc2dminim2stationary </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module invokes the nonlinear solver to solve the basic CC2D problem.
!#
!# 1.) c2d2_solve
!#     -> Invoke the nonlinear solver to solve the stationary coupled system.
!#
!# </purpose>
!##############################################################################

MODULE cc2dmediumm2stationary

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
    
  USE cc2dmediumm2basic
  USE cc2dmedium_callback
  
  USE cc2dmediumm2nonlinearcore
  
  IMPLICIT NONE

CONTAINS
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE c2d2_solve (rproblem,rvector,rrhs)
  
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
  ! A problem astructure saving problem-dependent information.
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
    TYPE(t_vectorBlock), TARGET :: rtempBlock
    TYPE(t_ccNonlinearIteration) :: rnonlinearIteration

    ! The nonlinear solver configuration
    TYPE(t_nlsolNode) :: rnlSol
    
    ! Initialise the nonlinear solver node rnlSol with parameters from
    ! the INI/DAT files.
    CALL c2d2_getNonlinearSolver (rnlSol, rproblem%rparamList, 'CC2D-NONLINEAR')

    ! Initialise the nonlinear loop. This is to prepare everything for
    ! or callback routines that are called from the nonlinear solver.
    ! The preconditioner in that structure is initialised later.
    CALL c2d2_initNonlinearLoop (&
        rproblem,rvector,rrhs,rnonlinearIteration,'CC2D-NONLINEAR')
        
    ! Initialise the core equation to the stationary (Navier)-Stokes.
    CALL c2d2_setupCoreEquation (rnonlinearIteration,&
        0.0_DP,1.0_DP,REAL(1-collct_getvalue_int (rproblem%rcollection,'ISTOKES'),DP))
    
    ! Check the matrices if they are compatible to our
    ! preconditioner. If not, we have to modify the matrices a little
    ! bit to make it compatible.
    CALL c2d2_checkAssembly (rproblem,rrhs,rnonlinearIteration)
    
    ! Make the matrices compatible if they are not
    CALL c2d2_finaliseMatrices (rnonlinearIteration)
    
    ! Initialise the preconditioner for the nonlinear iteration
    CALL c2d2_preparePreconditioner (rproblem,&
        rnonlinearIteration%rpreconditioner,rvector,rrhs)

    ! Save the nonlinear-iteration structure to the collection.
    ! It's rebuild from the collection in the callback routines.
    CALL c2d2_saveNonlinearLoop (rnonlinearIteration,rproblem%rcollection)
    
    ! Create a temporary vector we need for the nonlinear iteration.
    CALL lsysbl_createVecBlockIndirect (rrhs, rtempBlock, .FALSE.)

    ! Call the nonlinear solver. For preconditioning and defect calculation,
    ! the solver calls our callback routines.
    CALL nlsol_performSolve(rnlSol,rvector,rrhs,rtempBlock,&
                            c2d2_getDefect,c2d2_precondDefect,c2d2_resNormCheck,&
                            rcollection=rproblem%rcollection)

    ! Release the temporary vector
    CALL lsysbl_releaseVector (rtempBlock)
    
    ! Release parameters of the nonlinear loop
    CALL c2d2_doneNonlinearLoop (rnonlinearIteration,rproblem%rcollection)
             
    ! Release the preconditioner
    CALL c2d2_releasePreconditioner (rproblem,rnonlinearIteration%rpreconditioner)
    
    CALL output_lbrk()
    CALL output_line ('Nonlinear solver statistics')
    CALL output_line ('---------------------------')
    CALL output_line ('Intial defect: '//TRIM(sys_sdEL(rnlSol%DinitialDefect(1),15)))
    CALL output_line ('Final defect:  '//TRIM(sys_sdEL(rnlSol%DfinalDefect(1),15)))
    CALL output_line ('#Iterations:   '//TRIM(sys_siL(rnlSol%iiterations,10)))
    
  END SUBROUTINE

END MODULE
