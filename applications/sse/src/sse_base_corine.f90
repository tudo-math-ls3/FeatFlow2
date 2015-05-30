!##############################################################################
!# ****************************************************************************
!# <name> sse_base_corine </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module provides the basic routines and constants for Corine`s
!# problem.
!#
!# 1.) sse_initParamCorine
!#     -> Initialises the parameters of Corine`s problem
!#
!# 2.) sse_infoCorine
!#     -> Output information about Corine`s problem
!#
!# </purpose>
!##############################################################################

module sse_base_corine

  use fsystem
  use genoutput
  use paramlist

  use sse_base

  implicit none

  private
  public :: sse_initParamCorine
  public :: sse_infoCorine

!<constants>

!<constantblock description="Constants for problem subtypes">

  ! Corine`s problem (add more if needed)
  integer, parameter, public :: CORINE_STD = 0
!</constantblock>

!</constants>

!<publicvars>

  ! Bed friction in the momentum equation (=r) 
  real(DP), public :: dr                  = 0.0_DP

  ! Gravitational constant in the momentum equation (=delta)
  real(DP), public :: dDelta              = 0.0_DP

  ! Contribution of bedload transport in the bed evolution equation (=mu)
  real(DP), public :: dmu                 = 0.0_DP

  ! Constant to prevent the bed shear stress from blowing up (=h0) 
  real(DP), public :: dh0                 = 0.0_DP

  ! Ratio of diffusive timescale and the tidal period (=kappa) 
  real(DP), public :: dkappa              = 0.0_DP

  ! Ratio of the deposition timescale and the tidal period (=a) 
  real(DP), public :: da                  = 0.0_DP

  ! Width of the basin (=W) 
  real(DP), public :: dW                  = 0.0_DP

!</publicvars>

contains

  ! ***************************************************************************

!<subroutine>

  subroutine sse_initParamCorine(cproblemType,rparlist)

!<description>
    ! This subroutine initialises the global parameters of Corine`s problem
!</description>

!<input>
    ! Problem type
    integer, intent(in) :: cproblemType

    ! Parameter list
    type(t_parlist), intent(in) :: rparlist
!</input>
!</subroutine>

    ! local variables
    character(len=SYS_STRLEN) :: sparam,ssection

    ! Read config section
    ssection = sse_getSection(cproblemType)
    call parlst_getvalue_string(rparlist, ssection, 'problemparam', sparam)

    ! Read parameters from parameter list (non-existing parameters are replaced by maximum value)
    call parlst_getvalue_double(rparlist, sparam, 'dDelta', dDelta, SYS_MAXREAL_DP)
    call parlst_getvalue_double(rparlist, sparam, 'dW',     dW,     SYS_MAXREAL_DP)
    call parlst_getvalue_double(rparlist, sparam, 'da',     da,     SYS_MAXREAL_DP)
    call parlst_getvalue_double(rparlist, sparam, 'dh0',    dh0,    SYS_MAXREAL_DP)
    call parlst_getvalue_double(rparlist, sparam, 'dkappa', dkappa, SYS_MAXREAL_DP)
    call parlst_getvalue_double(rparlist, sparam, 'dmu',    dmu,    SYS_MAXREAL_DP)
    call parlst_getvalue_double(rparlist, sparam, 'dr',     dr,     SYS_MAXREAL_DP)

  end subroutine sse_initParamCorine

  ! ***************************************************************************

!<subroutine>

  subroutine sse_infoCorine(cproblemType)

!<description>   
    ! This subroutine outputs information about the global parameters
    ! of Corine`s problem
!</description>

!<input>
    ! Problem type
    integer, intent(in) :: cproblemType
!</input>
!</subroutine>

    select case (cproblemType)
    case (CORINE_1D)
      call output_line('PROBLEMTYPE........: Corine`s problem in 1D')
    case (CORINE_2D)
      call output_line('PROBLEMTYPE........: Corine`s problem in 2D')
    case default
      call output_line("Invalid problem type", &
          OU_CLASS_ERROR,OU_MODE_STD,"sse_infoCorine")
      call sys_halt()
    end select

    call output_line('Delta..............: '//trim(adjustl(sys_sdE(dDelta,5))))
    call output_line('W..................: '//trim(adjustl(sys_sdE(dW,5))))
    call output_line('a..................: '//trim(adjustl(sys_sdE(da,5))))
    call output_line('h0.................: '//trim(adjustl(sys_sdE(dh0,5))))
    call output_line('kappa..............: '//trim(adjustl(sys_sdE(dkappa,5))))
    call output_line('mu.................: '//trim(adjustl(sys_sdE(dmu,5))))
    call output_line('r..................: '//trim(adjustl(sys_sdE(dr,5))))

  end subroutine sse_infoCorine

end module sse_base_corine
