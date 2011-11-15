!##############################################################################
!# ****************************************************************************
!# <name> ccstationary </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module invokes the nonlinear solver to solve the basic CC2D problem.
!#
!# 1.) cc_solveStationary
!#     -> Invoke the nonlinear solver to solve the stationary coupled system.
!#
!# </purpose>
!##############################################################################

module ccstationary

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
    
  use ccbasic
  use cccallback
  
  use ccnonlinearcore
  use ccnonlinearcoreinit
  
  implicit none

contains
  
  ! ***************************************************************************

!<subroutine>

  subroutine cc_solveStationary (rproblem,rvector,rrhs)
  
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
  type(t_problem), intent(inout), target :: rproblem

  ! The solution vector which is to be used as initial vector for the nonlinear
  ! iteration. After the iteration, this is replaced by the new solution vector.
  type(t_vectorBlock), intent(inout) :: rvector
!</inputoutput>

!<input>
  ! The right-hand-side vector to use in the equation
  type(t_vectorBlock), intent(in) :: rrhs
!</input>

!</subroutine>

    ! local variables
    type(t_ccNonlinearIteration) :: rnonlinearIteration

    ! The nonlinear solver configuration
    type(t_nlsolNode) :: rnlSol
    
    ! Initialise the nonlinear solver node rnlSol with parameters from
    ! the INI/DAT files.
    call cc_getNonlinearSolver (rnlSol, rproblem%rparamList, 'CC2D-NONLINEAR')

    ! Initialise the nonlinear loop. This is to prepare everything for
    ! or callback routines that are called from the nonlinear solver.
    ! The preconditioner in that structure is initialised later.
    call cc_initNonlinearLoop (&
        rproblem,rproblem%NLMIN,rproblem%NLMAX,rvector,rrhs,&
        rnonlinearIteration,'CC2D-NONLINEAR')
        
    ! Set up all the weights in the core equation according to the current timestep.
    rnonlinearIteration%dalpha1 = 0.0_DP
    rnonlinearIteration%dalpha2 = 0.0_DP

    rnonlinearIteration%dtheta1 = 1.0_DP
    rnonlinearIteration%dtheta2 = 1.0_DP

    rnonlinearIteration%dgamma1 = real(1-rproblem%rphysics%iequation1,DP)
    rnonlinearIteration%dgamma2 = real(1-rproblem%rphysics%iequation2,DP)

    rnonlinearIteration%deta1   = 1.0_DP
    rnonlinearIteration%deta2   = 1.0_DP

    rnonlinearIteration%dtau1   = rproblem%rphysics%dmu1
    rnonlinearIteration%dtau2   = rproblem%rphysics%dmu2

    ! Initialise the preconditioner for the nonlinear iteration
    call cc_initPreconditioner (rproblem,&
        rnonlinearIteration,rvector,rrhs)

    ! Call the nonlinear solver to solve the core equation.
    call cc_solveCoreEquation (rproblem,rnonlinearIteration,rnlSol,&
        rvector,rrhs)             
             
    ! Release the preconditioner
    call cc_releasePreconditioner (rnonlinearIteration)
    
    ! Release parameters of the nonlinear loop, final clean up
    call cc_doneNonlinearLoop (rnonlinearIteration)
             
    call output_lbrk()
    call output_line ('Nonlinear solver statistics',coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG)
    call output_line ('---------------------------',coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG)
    call output_line ('Initial defect: '//trim(sys_sdEL(rnlSol%DinitialDefect(1),15)),&
        coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG)
    call output_line ('Final defect:  '//trim(sys_sdEL(rnlSol%DfinalDefect(1),15)),&
        coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG)
    call output_line ('#Iterations:   '//trim(sys_siL(rnlSol%iiterations,10)),&
        coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG)
    
  end subroutine

end module
