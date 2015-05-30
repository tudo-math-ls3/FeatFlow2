!##############################################################################
!# ****************************************************************************
!# <name> sse_base_poisson </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module provides the basic routines and constants for the
!# Poisson problem.
!#
!# 1.) sse_initParamPoisson
!#     -> Initialises the parameters of the Poisson problem
!#
!# 2.) sse_infoPoisson
!#     -> Output information about the Poisson problem
!#
!# </purpose>
!##############################################################################

module sse_base_poisson

  use fparser
  use fsystem
  use genoutput
  use paramlist

  use sse_base

  implicit none

  private
  public :: sse_initParamPoisson
  public :: sse_infoPoisson

!<publicvars>

  ! Scaling parameter
  real(DP), public :: dpoisson  = 1.0_DP

  ! Expression for scaling parameter
  integer, public ::  cpoisson  = 0

  ! Expression for right-hand side
  integer, public ::  crhs      = 0

  ! Expression for solution and its derivatives
  integer, public :: csol       = 0
  integer, public :: csol_x     = 0
  integer, public :: csol_y     = 0
  integer, public :: csol_xx    = 0
  integer, public :: csol_xy    = 0
  integer, public :: csol_yy    = 0
  integer, public :: csol_xxx   = 0
  integer, public :: csol_xxy   = 0
  integer, public :: csol_xyy   = 0
  integer, public :: csol_yyy   = 0

!</publicvars>

contains

  ! ***************************************************************************

!<subroutine>

  subroutine sse_initParamPoisson(cproblemType,rparlist)

!<description>
    ! This subroutine initialises the global parameters of the Poison problem
!</description>

!<input>
    ! Problem type
    integer, intent(in) :: cproblemType

    ! Parameter list
    type(t_parlist), intent(in) :: rparlist
!</input>
!</subroutine>

    ! local variable
    character(len=SYS_STRLEN) :: sparam,ssection,sexpression

    ! Read config section
    ssection = sse_getSection(cproblemType)
    call parlst_getvalue_string(rparlist, ssection, 'problemparam', sparam)

    ! Read parameters from parameter list (non-existing parameters are replaced by maximum value)
    call parlst_getvalue_double(rparlist, sparam, 'dpoisson', dpoisson,    SYS_MAXREAL_DP)

    ! Parse expression for scaling parameter: cpoisson
    call parlst_getvalue_string(rparlist, sparam, 'cpoisson', sexpression, '')
    if (trim(adjustl(sexpression)) .ne. '') then
      call fparser_parseFunction(rfparser, 'cpoisson', sexpression, (/'x','y'/), icomp=cpoisson)
    end if

    ! Parse expression for right-hand side: crhs
    call parlst_getvalue_string(rparlist, sparam, 'crhs', sexpression, '')
    if (trim(adjustl(sexpression)) .ne. '') then
      call fparser_parseFunction(rfparser, 'crhs', sexpression, (/'x','y'/), icomp=crhs)
    end if

    ! Parse expression for solution: csol
    call parlst_getvalue_string(rparlist, sparam, 'csol', sexpression, '')
    if (trim(adjustl(sexpression)) .ne. '') then
      call fparser_parseFunction(rfparser, 'csol', sexpression, (/'x','y'/), icomp=csol)
    end if

    ! Parse expression for x-derivative of solution: csol_x
    call parlst_getvalue_string(rparlist, sparam, 'csol_x', sexpression, '')
    if (trim(adjustl(sexpression)) .ne. '') then
      call fparser_parseFunction(rfparser, 'csol_x', sexpression, (/'x','y'/), icomp=csol_x)
    end if

    ! Parse expression for y-derivative of solution: csol_y
    call parlst_getvalue_string(rparlist, sparam, 'csol_y', sexpression, '')
    if (trim(adjustl(sexpression)) .ne. '') then
      call fparser_parseFunction(rfparser, 'csol_y', sexpression, (/'x','y'/), icomp=csol_y)
    end if

    ! Parse expression for xx-derivative of solution: csol_xx
    call parlst_getvalue_string(rparlist, sparam, 'csol_xx', sexpression, '')
    if (trim(adjustl(sexpression)) .ne. '') then
      call fparser_parseFunction(rfparser, 'csol_xx', sexpression, (/'x','y'/), icomp=csol_xx)
    end if

    ! Parse expression for xy-derivative of solution: csol_xy
    call parlst_getvalue_string(rparlist, sparam, 'csol_xy', sexpression, '')
    if (trim(adjustl(sexpression)) .ne. '') then
      call fparser_parseFunction(rfparser, 'csol_xy', sexpression, (/'x','y'/), icomp=csol_xy)
    end if

    ! Parse expression for yy-derivative of solution: csol_yy
    call parlst_getvalue_string(rparlist, sparam, 'csol_yy', sexpression, '')
    if (trim(adjustl(sexpression)) .ne. '') then
      call fparser_parseFunction(rfparser, 'csol_yy', sexpression, (/'x','y'/), icomp=csol_yy)
    end if

    ! Parse expression for xxx-derivative of solution: csol_xxx
    call parlst_getvalue_string(rparlist, sparam, 'csol_xxx', sexpression, '')
    if (trim(adjustl(sexpression)) .ne. '') then
      call fparser_parseFunction(rfparser, 'csol_xxx', sexpression, (/'x','y'/), icomp=csol_xxx)
    end if

    ! Parse expression for xxy-derivative of solution: csol_xxy
    call parlst_getvalue_string(rparlist, sparam, 'csol_xxy', sexpression, '')
    if (trim(adjustl(sexpression)) .ne. '') then
      call fparser_parseFunction(rfparser, 'csol_xxy', sexpression, (/'x','y'/), icomp=csol_xxy)
    end if

    ! Parse expression for xyy-derivative of solution: csol_xyy
    call parlst_getvalue_string(rparlist, sparam, 'csol_xyy', sexpression, '')
    if (trim(adjustl(sexpression)) .ne. '') then
      call fparser_parseFunction(rfparser, 'csol_xyy', sexpression, (/'x','y'/), icomp=csol_xyy)
    end if

    ! Parse expression for yyy-derivative of solution: csol_yyy
    call parlst_getvalue_string(rparlist, sparam, 'csol_yyy', sexpression, '')
    if (trim(adjustl(sexpression)) .ne. '') then
      call fparser_parseFunction(rfparser, 'csol_yyy', sexpression, (/'x','y'/), icomp=csol_yyy)
    end if

  end subroutine sse_initParamPoisson

  ! ***************************************************************************

!<subroutine>

  subroutine sse_infoPoisson(cproblemType)

!<description>   
    ! This subroutine outputs information about the global parameters
    ! of the Poisson problem
!</description>

!<input>
    ! Problem type
    integer, intent(in) :: cproblemType
!</input>
!</subroutine>

    select case (cproblemType)
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

    call output_line('DPOISSON...........: '//trim(adjustl(sys_sdE(dpoisson,5))))
    call output_line('CPOISSON...........: '//trim(adjustl(sys_si(cpoisson,3))))
    call output_line('CRHS...............: '//trim(adjustl(sys_si(crhs,3))))
    call output_line('CSOL...............: '//trim(adjustl(sys_si(csol,3))))
    call output_line('CSOL_X.............: '//trim(adjustl(sys_si(csol_x,3))))
    call output_line('CSOL_Y.............: '//trim(adjustl(sys_si(csol_y,3))))
    call output_line('CSOL_XX............: '//trim(adjustl(sys_si(csol_xx,3))))
    call output_line('CSOL_XY............: '//trim(adjustl(sys_si(csol_xy,3))))
    call output_line('CSOL_YY............: '//trim(adjustl(sys_si(csol_yy,3))))
    call output_line('CSOL_XXX...........: '//trim(adjustl(sys_si(csol_xxx,3))))
    call output_line('CSOL_XXY...........: '//trim(adjustl(sys_si(csol_xxy,3))))
    call output_line('CSOL_XYY...........: '//trim(adjustl(sys_si(csol_xyy,3))))
    call output_line('CSOL_YYY...........: '//trim(adjustl(sys_si(csol_yyy,3))))

  end subroutine sse_infoPoisson

end module sse_base_poisson
