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

module cc2dmediumm2stationary

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
    
  use cc2dmediumm2basic
  use cc2dmedium_callback
  
  use cc2dmediumm2nonlinearcore
  use cc2dmediumm2nonlinearcoreinit
  
  implicit none

contains
  
  ! ***************************************************************************

!<subroutine>

  subroutine c2d2_solve (rproblem,rvector,rrhs)
  
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

    ! local variables
    type(t_ccNonlinearIteration) :: rnonlinearIteration

    ! The nonlinear solver configuration
    type(t_nlsolNode) :: rnlSol
    
    ! Initialise the nonlinear solver node rnlSol with parameters from
    ! the INI/DAT files.
    call c2d2_getNonlinearSolver (rnlSol, rproblem%rparamList, 'CC2D-NONLINEAR')

    ! Initialise the nonlinear loop. This is to prepare everything for
    ! or callback routines that are called from the nonlinear solver.
    ! The preconditioner in that structure is initialised later.
    call c2d2_initNonlinearLoop (&
        rproblem,rproblem%NLMIN,rproblem%NLMAX,rvector,rrhs,&
        rnonlinearIteration,'CC2D-NONLINEAR')
        
    ! Set up all the weights in the core equation according to the current timestep.
    rnonlinearIteration%dalpha = 0.0_DP
    rnonlinearIteration%dtheta = 1.0_DP
    rnonlinearIteration%dgamma = real(1-rproblem%iequation,DP)
    rnonlinearIteration%deta   = 1.0_DP
    rnonlinearIteration%dtau   = 1.0_DP

    ! Initialise the preconditioner for the nonlinear iteration
    call c2d2_initPreconditioner (rproblem,&
        rnonlinearIteration,rvector,rrhs)

    ! Call the nonlinear solver to solve the core equation.
    call c2d2_solveCoreEquation (rnlSol,rnonlinearIteration,&
        rvector,rrhs,rproblem%rcollection)
             
    ! Release the preconditioner
    call c2d2_releasePreconditioner (rnonlinearIteration)
    
    ! Release parameters of the nonlinear loop, final clean up
    call c2d2_doneNonlinearLoop (rnonlinearIteration)
             
    call output_lbrk()
    call output_line ('Nonlinear solver statistics')
    call output_line ('---------------------------')
    call output_line ('Initial defect: '//trim(sys_sdEL(rnlSol%DinitialDefect(1),15)))
    call output_line ('Final defect:  '//trim(sys_sdEL(rnlSol%DfinalDefect(1),15)))
    call output_line ('#Iterations:   '//trim(sys_siL(rnlSol%iiterations,10)))
    
  end subroutine

end module
