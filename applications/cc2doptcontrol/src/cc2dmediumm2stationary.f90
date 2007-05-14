!##############################################################################
!# ****************************************************************************
!# <name> cc2dminim2stationary </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module invokes the nonlinear solver to solve the basic CC2D 
!# optimal control problem.
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
  USE cc2dmediumm2nonlinearcoreinit
  USE cc2dmediumm2postprocessing
  
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
    REAL(DP) :: dalphaC
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
        rproblem,rproblem%NLMIN,rproblem%NLMAX,rvector,rrhs,&
        rnonlinearIteration,'CC2D-NONLINEAR')
        
    ! Set up all the weights in the core equation according to the current timestep.
    CALL parlst_getvalue_double (rproblem%rparamList,'OPTIMALCONTROL',&
                                'dalphaC',dalphaC,0.1_DP)
    rnonlinearIteration%diota1  = 0.0_DP
    rnonlinearIteration%dkappa1 = 0.0_DP
    rnonlinearIteration%dalpha1 = 0.0_DP
    rnonlinearIteration%dtheta1 = 1.0_DP
    rnonlinearIteration%dgamma1 = REAL(1-rproblem%iequation,DP)
    rnonlinearIteration%deta1   = 1.0_DP
    rnonlinearIteration%dtau1   = 1.0_DP
    rnonlinearIteration%dmu1    = 1.0_DP/dalphaC

    rnonlinearIteration%diota2  = 0.0_DP
    rnonlinearIteration%dkappa2 = 0.0_DP
    rnonlinearIteration%dalpha2 = 0.0_DP
    rnonlinearIteration%dtheta2 = 1.0_DP
    rnonlinearIteration%dgamma2 = REAL(1-rproblem%iequation,DP)
    rnonlinearIteration%deta2   = 1.0_DP
    rnonlinearIteration%dtau2   = 1.0_DP
    rnonlinearIteration%dmu2    = -1.0_DP
    
    ! Check the matrices if they are compatible to our
    ! preconditioner. If not, we later have to modify the matrices a little
    ! bit to make it compatible. 
    ! The result of this matrix analysis is saved to the rfinalAssembly structure 
    ! in rnonlinearIteration and allows us later to switch between these two
    ! matrix representations: Compatibility to the discretisation routines
    ! and compatibity to the preconditioner.
    ! The c2d2_checkAssembly routine below uses this information to perform
    ! the actual modification in the matrices.
    CALL c2d2_checkAssembly (rproblem,rnonlinearIteration,rrhs,&
        rnonlinearIteration%rfinalAssembly)
    
    ! Using rfinalAssembly as computed above, make the matrices compatible 
    ! to our preconditioner if they are not.
    CALL c2d2_finaliseMatrices (rnonlinearIteration)
    
    ! Initialise the preconditioner for the nonlinear iteration
    CALL c2d2_preparePreconditioner (rproblem,&
        rnonlinearIteration,rvector,rrhs)
        
    ! Print out the value of the optimal control functional for the
    ! initial solution vector directly prior to the solution process.
    CALL c2d2_printControlFunctionalStat (rproblem,rvector)

    ! Call the nonlinear solver to solve the core equation.
    CALL c2d2_solveCoreEquation (rnlSol,rnonlinearIteration,&
        rvector,rrhs,rproblem%rcollection)             
             
    ! Release the preconditioner
    CALL c2d2_releasePreconditioner (rnonlinearIteration)
    
    ! Release parameters of the nonlinear loop, final clean up
    CALL c2d2_doneNonlinearLoop (rnonlinearIteration)
             
    CALL output_lbrk()
    CALL output_line ('Nonlinear solver statistics')
    CALL output_line ('---------------------------')
    CALL output_line ('Intial defect: '//TRIM(sys_sdEL(rnlSol%DinitialDefect(1),15)))
    CALL output_line ('Final defect:  '//TRIM(sys_sdEL(rnlSol%DfinalDefect(1),15)))
    CALL output_line ('#Iterations:   '//TRIM(sys_siL(rnlSol%iiterations,10)))
    
  END SUBROUTINE

END MODULE
