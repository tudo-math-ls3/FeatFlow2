!##############################################################################
!# ****************************************************************************
!# <name> sse_base_poisson </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module provides the basic routines and constants for the
!# Poisson problem. 
!# </purpose>
!##############################################################################

module sse_base_poisson

  use fsystem
  use genoutput
  use paramlist

  use sse_base
  
  implicit none

  private
  public :: sse_initParamPoisson
  public :: sse_infoPoisson

!<constants>

!<constantblock description="Constants for problem subtypes">

  ! Pure Dirichlet Poisson problem
  integer, parameter, public :: POISSON_DIRICHLET         = 0

  ! Mixed Dirichlet-Neumann problem
  integer, parameter, public :: POISSON_DIRICHLET_NEUMANN = 1

  ! Pure Neumann Poisson problem
  integer, parameter, public :: POISSON_NEUMANN           = 2
!</constantblock>

!</constants>
  
!<publicvars>

  ! Scaling parameter
  real(DP), public :: dpoisson   = 0.0_DP
  
  ! Dirichlet boundary value (=u_D)
  real(DP), public :: ddirichlet = 0.0_DP

  ! Neumann boundary value (=g)
  real(DP), public :: dneumann   = 0.0_DP

!</publicvars>

contains

  ! ***************************************************************************

!<subroutine>

  subroutine sse_initParamPoisson(cproblemtype,rparlist)

!<description>
    ! This subroutine initialises the global parameters of the Poison problem
!</description>

!<input>
    ! Problem type
    integer, intent(in) :: cproblemtype
    
    ! Parameter list
    type(t_parlist), intent(in) :: rparlist
!</input>
!</subroutine>

    ! local variable
    character(len=SYS_STRLEN) :: sconfig,ssection

    ! Read config section
    ssection = sse_getSection(cproblemtype)
    call parlst_getvalue_string(rparlist, ssection, 'problemconfig', sconfig)
    
    ! Read parameters from parameter list (non-existing parameters are replaced by maximum value)
    call parlst_getvalue_double(rparlist, sconfig, 'dpoisson',   dpoisson,   SYS_MAXREAL_DP)
    call parlst_getvalue_double(rparlist, sconfig, 'ddirichlet', ddirichlet, SYS_MAXREAL_DP)
    call parlst_getvalue_double(rparlist, sconfig, 'dneumann',   dneumann,   SYS_MAXREAL_DP)
    
  end subroutine sse_initParamPoisson

  ! ***************************************************************************

!<subroutine>

  subroutine sse_infoPoisson(cproblemtype,cproblemsubtype)

!<description>   
    ! This subroutine outputs information about the global parameters
    ! of the Poisson problem
!</description>

!<input>
    ! Problem type
    integer, intent(in) :: cproblemtype
    
    ! Problem subtype
    integer, intent(in) :: cproblemsubtype
!</input>
!</subroutine>

    select case (cproblemtype)
    case (POISSON_SCALAR)
      call output_line('PROBLEMTYPE........: Poisson problem')
      call output_line('FORMULATION........: second-order equation')
    case (POISSON_SYSTEM)
      call output_line('PROBLEMTYPE........: Poisson problem')
      call output_line('FORMULATION........: first-order formulation')
    case default
      call output_line("Invalid problem type", &
        OU_CLASS_ERROR,OU_MODE_STD,"sse_infoPoisson")
      call sys_halt()
    end select

    select case (cproblemsubtype)
    case (POISSON_DIRICHLET)
      call output_line('PROBLEMSUBTYPE.....: Pure Dirichlet problem')
      call output_line('Poisson............: '//trim(adjustl(sys_sdE(dpoisson,5))))
      call output_line('Dirichlet..........: '//trim(adjustl(sys_sdE(ddirichlet,5))))
      
    case (POISSON_DIRICHLET_NEUMANN)
      call output_line('PROBLEMSUBTYPE.....: Mixed Dirichlet-Neumann problem')
      call output_line('Poisson............: '//trim(adjustl(sys_sdE(dpoisson,5))))
      call output_line('Dirichlet..........: '//trim(adjustl(sys_sdE(ddirichlet,5))))
      call output_line('Neumann............: '//trim(adjustl(sys_sdE(dneumann,5))))
      
    case (POISSON_NEUMANN)
      call output_line('PROBLEMSUBTYPE.....: Pure Neumann problem')
      call output_line('Poisson............: '//trim(adjustl(sys_sdE(dpoisson,5))))
      call output_line('Neumann............: '//trim(adjustl(sys_sdE(dneumann,5))))
      
    case default
      call output_line("Invalid problem subtype", &
          OU_CLASS_ERROR,OU_MODE_STD,"sse_infoPoisson")
      call sys_halt()
    end select

  end subroutine sse_infoPoisson
    
end module sse_base_poisson
