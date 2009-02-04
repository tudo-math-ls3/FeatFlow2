!##############################################################################
!# ****************************************************************************
!# <name> euler_callback </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains all callback functions which are required to solve the 
!# compressible Euler/Navier staokes equations in arbitrary spatial dimensions.
!#
!# The following callback functions are available:
!#
!# 1.) euler_calcPreconditioner
!#     -> Calculates the nonlinear preconditioner
!# </purpose>
!##############################################################################

module euler_callback

  use afcstabilisation
  use boundaryfilter
  use collection
  use euler_basic
  use euler_callback2d
  use fsystem
  use genoutput
  use linearsystemblock
  use linearsystemscalar
  use problem
  use solver
  use statistics
  use storage

  implicit none

  private
  public :: euler_calcPreconditioner

contains

  !*****************************************************************************

!<subroutine>

  subroutine euler_calcPreconditioner(rproblemLevel, rtimestep, rsolver,&
                                       rsol, rcollection)

!<description>
    ! This subroutine calculates the nonlinear preconditioner and
    ! configures the linear solver structure accordingly. Depending on
    ! the nonlinear solver, the low-order evolution operator or the
    ! Jacobian operator is adopted as nonlinear preconditioner.
!</description>

!<input>
    ! time-stepping structure
    type(t_timestep), intent(IN) :: rtimestep

    ! solution vector
    type(t_vectorBlock), intent(IN) :: rsol
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(INOUT) :: rproblemLevel

    ! solver structure
    type(t_solver), intent(INOUT) :: rsolver

    ! collection
    type(t_collection), intent(INOUT) :: rcollection
!</inputoutput>
!</subroutine>
    
    ! local variables
    type(t_timer), pointer :: rtimer

    integer :: icoupled

    ! Start time measurement for matrix evaluation
    rtimer => collct_getvalue_timer(rcollection, 'timerAssemblyMatrix')
    call stat_startTimer(rtimer, STAT_TIMERSHORT)

    !---------------------------------------------------------------------------
    ! Assemble divergence operator
    !---------------------------------------------------------------------------

    icoupled = collct_getvalue_int(rcollection, 'icoupled')

    ! What kind of coupling is applied?
    select case(icoupled)
      
    ! Assemble block-diagonal divergence operator
    case (FLOW_SEGREGATED)

    ! Assemble full block transport operator
    case (FLOW_ALLCOUPLED)

    case DEFAULT
      call output_line('Invalid flow coupling!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'euler_calcPreconditioner')
      call sys_halt()
    end select

  end subroutine euler_calcPreconditioner

end module euler_callback
