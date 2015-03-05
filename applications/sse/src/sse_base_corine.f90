!##############################################################################
!# ****************************************************************************
!# <name> sse_base_corine </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module provides the basic routines and constants for Corine`s
!# problem.
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
  real(DP), public :: dr        = 0.0_DP

  ! Gravitational constant in the momentum equation (=delta)
  real(DP), public :: dDelta    = 0.0_DP

  ! Contribution of bedload transport in the bed evolution equation (=mu)
  real(DP), public :: dmu       = 0.0_DP

  ! Constant to prevent the bed shear stress from blowing up (=h0) 
  real(DP), public :: dh0       = 0.0_DP

  ! Ratio of diffusive timescale and the tidal period (=kappa) 
  real(DP), public :: dkappa    = 0.0_DP

  ! Ratio of the deposition timescale and the tidal period (=a) 
  real(DP), public :: da        = 0.0_DP

  ! Width of the basin (=W) 
  real(DP), public :: dW        = 0.0_DP

!</publicvars>

contains

  ! ***************************************************************************

!<subroutine>

  subroutine sse_initParamCorine(cproblemtype,rparlist)

!<description>
    ! This subroutine initialises the global parameters of Corine`s problem
!</description>

!<input>
    ! Problem type
    integer, intent(in) :: cproblemtype
    
    ! Parameter list
    type(t_parlist), intent(in) :: rparlist
!</input>
!</subroutine>

    ! local variables
    character(len=SYS_STRLEN) :: sconfig,ssection

    ! Read config section
    ssection = sse_getSection(cproblemtype)
    call parlst_getvalue_string(rparlist, ssection, 'problemconfig', sconfig)
    
    ! Read parameters from parameter list (non-existing parameters are replaced by maximum value)
    call parlst_getvalue_double(rparlist, sconfig, 'dr',     dr,     SYS_MAXREAL_DP)
    call parlst_getvalue_double(rparlist, sconfig, 'dDelta', dDelta, SYS_MAXREAL_DP)
    call parlst_getvalue_double(rparlist, sconfig, 'dmu',    dmu,    SYS_MAXREAL_DP)
    call parlst_getvalue_double(rparlist, sconfig, 'dh0',    dh0,    SYS_MAXREAL_DP)
    call parlst_getvalue_double(rparlist, sconfig, 'dkappa', dkappa, SYS_MAXREAL_DP)
    call parlst_getvalue_double(rparlist, sconfig, 'da',     da,     SYS_MAXREAL_DP)
    call parlst_getvalue_double(rparlist, sconfig, 'dW',     dW,     SYS_MAXREAL_DP)
    
  end subroutine sse_initParamCorine

  ! ***************************************************************************

!<subroutine>

  subroutine sse_infoCorine(cproblemtype,cproblemsubtype)

!<description>   
    ! This subroutine outputs information about the global parameters
    ! of Corine`s problem
!</description>

!<input>
    ! Problem type
    integer, intent(in) :: cproblemtype
    
    ! Problem subtype
    integer, intent(in) :: cproblemsubtype
!</input>
!</subroutine>

    select case (cproblemtype)
    case (CORINE_1D)
      call output_line('PROBLEMTYPE........: Corine`s problem in 1D')
    case (CORINE_2D)
      call output_line('PROBLEMTYPE........: Corine`s problem in 2D')
    case default
      call output_line("Invalid problem type", &
        OU_CLASS_ERROR,OU_MODE_STD,"sse_infoCorine")
      call sys_halt()
    end select

    select case (cproblemsubtype)
    case (CORINE_STD)
      call output_line('PROBLEMSUBTYPE.....: Standard benchmark')
    case default
      call output_line("Invalid problem subtype", &
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
