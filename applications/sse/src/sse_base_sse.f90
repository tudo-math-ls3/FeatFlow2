!##############################################################################
!# ****************************************************************************
!# <name> sse_base_sse </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module provides the basic routines and constants for the
!# SSE problem.
!#
!# 1.) sse_initParamSSE
!#     -> Initialises the parameters of the SSE problem
!#
!# 2.) sse_infoSSE
!#     -> Output information about the SSE  problem
!#
!# </purpose>
!##############################################################################

module sse_base_sse

  use fparser
  use fsystem
  use genoutput
  use paramlist

  use sse_base

  implicit none

  private
  public :: sse_initParamSSE
  public :: sse_infoSSE

!<constants>

!<constantblock description="Constants for bottom profiles">

  ! SSE linear bottom profile
  integer, parameter, public :: SSE_BOTTOMPROFILE_LINEAR    = 0

  ! SSE constant bottom profile
  integer, parameter, public :: SSE_BOTTOMPROFILE_CONSTANT  = 1

  ! SSE parabolic bottom profile
  integer, parameter, public :: SSE_BOTTOMPROFILE_PARABOLIC = 2
!</constantblock>

!<constantblock description="Constants for stress">

  ! SSE variable stress
  integer, parameter, public :: SSE_STRESS_VARIABLE     = 0

  ! SSE constant stress
  integer, parameter, public :: SSE_STRESS_CONSTANT     = 1

  ! SSE stress proportional to bathymetry
  integer, parameter, public :: SSE_STRESS_PROPORTIONAL = 2
!</constantblock>

!<constantblock description="Constants for eddy viscosity">

  ! SSE variable eddy viscosity
  integer, parameter, public :: SSE_VISCOSITY_VARIABLE     = 0

  ! SSE constant eddy viscosity
  integer, parameter, public :: SSE_VISCOSITY_CONSTANT     = 1

  ! SSE eddy viscosity proportional to bathymetry
  integer, parameter, public :: SSE_VISCOSITY_PROPORTIONAL = 2
!</constantblock>

!</constants>

!<publicvars>

  ! Length of the channel (=L)
  real(DP), public :: dlength             = 0.0_DP

  ! Convergent length of the channel (=Lb)
  real(DP), public :: dlengthB            = 0.0_DP

  ! Width of the entrance of the channel (=B)
  real(DP), public :: dwidth              = 0.0_DP

  ! Type of the bed profile
  integer, public :: ibathymetryType      = 0

  ! Mean depth of the channel (H and H0)
  real(DP), public :: dheight             = 0.0_DP
  real(DP), public :: dheight0            = 0.0_DP

  ! Depth at the end in the scaled domain (=a = dheight0/dheight)
  real(DP), public :: dheightRatio        = 0.0_DP

  ! Constant forcing at open boundary (=M2)
  real(DP), public :: dforcing            = 0.0_DP

  ! Coriolis acceleration (=f)
  real(DP), public :: dcoraccel           = 0.0_DP

  ! Frequency of the tidal constituent (=omega)
  real(DP), public :: dtidalfreq          = 0.0_DP

  ! Type of the bottom stress
  integer, public :: istressType          = 0

  ! Bottom stress (=s0)
  real(DP), public :: dstress             = 0.0_DP

  ! Type of the vertical eddy viscosity
  integer, public :: iviscosityType       = 0

  ! Vertical eddy viscosity (=Av0)
  real(DP), public :: dviscosity          = 0.0_DP

  ! Gravitational acceleration (=g)
  real(DP), public :: dgravaccel          = 0.0_DP

  ! Expression for solution and its derivatives
  integer, public :: csolReal             = 0
  integer, public :: csolImag             = 0
  integer, public :: csolReal_x           = 0
  integer, public :: csolImag_x           = 0
  integer, public :: csolReal_y           = 0
  integer, public :: csolImag_y           = 0
  integer, public :: csolReal_xx          = 0
  integer, public :: csolImag_xx          = 0
  integer, public :: csolReal_xy          = 0
  integer, public :: csolImag_xy          = 0
  integer, public :: csolReal_yy          = 0
  integer, public :: csolImag_yy          = 0
  integer, public :: csolReal_xxx         = 0
  integer, public :: csolImag_xxx         = 0
  integer, public :: csolReal_xxy         = 0
  integer, public :: csolImag_xxy         = 0
  integer, public :: csolReal_xyy         = 0
  integer, public :: csolImag_xyy         = 0
  integer, public :: csolReal_yyy         = 0
  integer, public :: csolImag_yyy         = 0

  ! Expression for bottom profile
  integer, public :: cbathymetry          = 0

  ! Expression for bottom stress
  integer, public :: cstress              = 0

  ! Expression for vertical eddy viscosity
  integer, public :: cviscosity           = 0
  
!</publicvars>

contains

  ! ***************************************************************************

!<subroutine>

  subroutine sse_initParamSSE(cproblemType,rparlist)

!<description>
    ! This subroutine initialises the global parameters of the SSE problem
!</description>

!<input>
    ! Problem type
    integer, intent(in) :: cproblemType

    ! Parameter list
    type(t_parlist), intent(in) :: rparlist
!</input>
!</subroutine>

    ! local variables
    character(len=SYS_STRLEN) :: sparam,ssection,sexpression

    ! Read config section
    ssection = sse_getSection(cproblemType)
    call parlst_getvalue_string(rparlist, ssection, 'problemparam', sparam)

    ! Read parameters from parameter list (non-existing parameters are replaced by maximum value)
    call parlst_getvalue_double(rparlist,  sparam, 'dcoraccel',      dcoraccel,   SYS_MAXREAL_DP)
    call parlst_getvalue_double(rparlist,  sparam, 'dforcing',       dforcing,    SYS_MAXREAL_DP)
    call parlst_getvalue_double(rparlist,  sparam, 'dgravaccel',     dgravaccel,  SYS_MAXREAL_DP)
    call parlst_getvalue_double(rparlist,  sparam, 'dheight',        dheight,     SYS_MAXREAL_DP)
    call parlst_getvalue_double(rparlist,  sparam, 'dheight0',       dheight0,    SYS_MAXREAL_DP)
    call parlst_getvalue_double(rparlist,  sparam, 'dheightRatio',   dheightRatio,SYS_MAXREAL_DP)
    call parlst_getvalue_double(rparlist,  sparam, 'dlength',        dlength,     SYS_MAXREAL_DP)
    call parlst_getvalue_double(rparlist,  sparam, 'dlengthB',       dlengthB,    SYS_MAXREAL_DP)
    call parlst_getvalue_double(rparlist,  sparam, 'dstress',        dstress,     SYS_MAXREAL_DP)
    call parlst_getvalue_double(rparlist,  sparam, 'dtidalfreq',     dtidalfreq,  SYS_MAXREAL_DP)
    call parlst_getvalue_double(rparlist,  sparam, 'dviscosity',     dviscosity,  SYS_MAXREAL_DP)
    call parlst_getvalue_double(rparlist,  sparam, 'dwidth',         dwidth,      SYS_MAXREAL_DP)
    call parlst_getvalue_int(rparlist,     sparam, 'ibathymetryType',ibathymetryType, SYS_MAXINT)
    call parlst_getvalue_int(rparlist,     sparam, 'istressType',    istressType,     SYS_MAXINT)
    call parlst_getvalue_int(rparlist,     sparam, 'iviscosityType', iviscosityType,  SYS_MAXINT)

    ! Parse expression for real part of solution: csolReal
    call parlst_getvalue_string(rparlist, sparam, 'csolReal', sexpression, '')
    if (trim(adjustl(sexpression)) .ne. '') then
      call fparser_parseFunction(rfparser, 'csolReal', sexpression, (/'x','y'/), icomp=csolReal)
    end if

    ! Parse expression for imaginary part of solution: csolImag
    call parlst_getvalue_string(rparlist, sparam, 'csolImag', sexpression, '')
    if (trim(adjustl(sexpression)) .ne. '') then
      call fparser_parseFunction(rfparser, 'csolImag', sexpression, (/'x','y'/), icomp=csolImag)
    end if

    ! Parse expression for real part of x-derivative of solution: csolReal_x
    call parlst_getvalue_string(rparlist, sparam, 'csolReal_x', sexpression, '')
    if (trim(adjustl(sexpression)) .ne. '') then
      call fparser_parseFunction(rfparser, 'csolReal_x', sexpression, (/'x','y'/), icomp=csolReal_x)
    end if

    ! Parse expression for imaginary part of x-derivative of solution: csolImag_x
    call parlst_getvalue_string(rparlist, sparam, 'csolImag_x', sexpression, '')
    if (trim(adjustl(sexpression)) .ne. '') then
      call fparser_parseFunction(rfparser, 'csolImag_x', sexpression, (/'x','y'/), icomp=csolImag_x)
    end if

    ! Parse expression for real part of y-derivative of solution: csolReal_y
    call parlst_getvalue_string(rparlist, sparam, 'csolReal_y', sexpression, '')
    if (trim(adjustl(sexpression)) .ne. '') then
      call fparser_parseFunction(rfparser, 'csolReal_y', sexpression, (/'x','y'/), icomp=csolReal_y)
    end if

    ! Parse expression for imaginary part of y-derivative of solution: csolImag_y
    call parlst_getvalue_string(rparlist, sparam, 'csolImag_y', sexpression, '')
    if (trim(adjustl(sexpression)) .ne. '') then
      call fparser_parseFunction(rfparser, 'csolImag_y', sexpression, (/'x','y'/), icomp=csolImag_y)
    end if

    ! Parse expression for real part of xx-derivative of solution: csolReal_xx
    call parlst_getvalue_string(rparlist, sparam, 'csolReal_xx', sexpression, '')
    if (trim(adjustl(sexpression)) .ne. '') then
      call fparser_parseFunction(rfparser, 'csolReal_xx', sexpression, (/'x','y'/), icomp=csolReal_xx)
    end if
        
    ! Parse expression for imaginary part of xx-derivative of solution: csolImag_xx
    call parlst_getvalue_string(rparlist, sparam, 'csolImag_xx', sexpression, '')
    if (trim(adjustl(sexpression)) .ne. '') then
      call fparser_parseFunction(rfparser, 'csolImag_xx', sexpression, (/'x','y'/), icomp=csolImag_xx)
    end if

    ! Parse expression for real part of xy-derivative of solution: csolReal_xy
    call parlst_getvalue_string(rparlist, sparam, 'csolReal_xy', sexpression, '')
    if (trim(adjustl(sexpression)) .ne. '') then
      call fparser_parseFunction(rfparser, 'csolReal_xy', sexpression, (/'x','y'/), icomp=csolReal_xy)
    end if

     ! Parse expression for imaginary part of xy-derivative of solution: csolImag_xy
    call parlst_getvalue_string(rparlist, sparam, 'csolImag_xy', sexpression, '')
    if (trim(adjustl(sexpression)) .ne. '') then
      call fparser_parseFunction(rfparser, 'csolImag_xy', sexpression, (/'x','y'/), icomp=csolImag_xy)
    end if

    ! Parse expression for real part of yy-derivative of solution: csolReal_yy
    call parlst_getvalue_string(rparlist, sparam, 'csolReal_yy', sexpression, '')
    if (trim(adjustl(sexpression)) .ne. '') then
      call fparser_parseFunction(rfparser, 'csolReal_yy', sexpression, (/'x','y'/), icomp=csolReal_yy)
    end if

    ! Parse expression for imaginary part of yy-derivative of solution: csolImag_yy
    call parlst_getvalue_string(rparlist, sparam, 'csolImag_yy', sexpression, '')
    if (trim(adjustl(sexpression)) .ne. '') then
      call fparser_parseFunction(rfparser, 'csolImag_yy', sexpression, (/'x','y'/), icomp=csolImag_yy)
    end if

    ! Parse expression for real part of xxx-derivative of solution: csolReal_xxx
    call parlst_getvalue_string(rparlist, sparam, 'csolReal_xxx', sexpression, '')
    if (trim(adjustl(sexpression)) .ne. '') then
      call fparser_parseFunction(rfparser, 'csolReal_xxx', sexpression, (/'x','y'/), icomp=csolReal_xxx)
    end if

    ! Parse expression for imaginary part of xxx-derivative of solution: csolImag_xxx
    call parlst_getvalue_string(rparlist, sparam, 'csolImag_xxx', sexpression, '')
    if (trim(adjustl(sexpression)) .ne. '') then
      call fparser_parseFunction(rfparser, 'csolImag_xxx', sexpression, (/'x','y'/), icomp=csolImag_xxx)
    end if

    ! Parse expression for real part of xxy-derivative of solution: csolReal_xxy
    call parlst_getvalue_string(rparlist, sparam, 'csolReal_xxy', sexpression, '')
    if (trim(adjustl(sexpression)) .ne. '') then
      call fparser_parseFunction(rfparser, 'csolReal_xxy', sexpression, (/'x','y'/), icomp=csolReal_xxy)
    end if

    ! Parse expression for imaginary part of xxy-derivative of solution: csolImag_xxy
    call parlst_getvalue_string(rparlist, sparam, 'csolImag_xxy', sexpression, '')
    if (trim(adjustl(sexpression)) .ne. '') then
      call fparser_parseFunction(rfparser, 'csolImag_xxy', sexpression, (/'x','y'/), icomp=csolImag_xxy)
    end if

    ! Parse expression for real part of xyy-derivative of solution: csolReal_xyy
    call parlst_getvalue_string(rparlist, sparam, 'csolReal_xyy', sexpression, '')
    if (trim(adjustl(sexpression)) .ne. '') then
      call fparser_parseFunction(rfparser, 'csolReal_xyy', sexpression, (/'x','y'/), icomp=csolReal_xyy)
    end if

    ! Parse expression for imaginary part of xyy-derivative of solution: csolImag_xyy
    call parlst_getvalue_string(rparlist, sparam, 'csolImag_xyy', sexpression, '')
    if (trim(adjustl(sexpression)) .ne. '') then
      call fparser_parseFunction(rfparser, 'csolImag_xyy', sexpression, (/'x','y'/), icomp=csolImag_xyy)
    end if

    ! Parse expression for real part of yyy-derivative of solution: csolReal_yyy
    call parlst_getvalue_string(rparlist, sparam, 'csolReal_yyy', sexpression, '')
    if (trim(adjustl(sexpression)) .ne. '') then
      call fparser_parseFunction(rfparser, 'csolReal_yyy', sexpression, (/'x','y'/), icomp=csolReal_yyy)
    end if

    ! Parse expression for imaginary part of yyy-derivative of solution: csolImag_yyy
    call parlst_getvalue_string(rparlist, sparam, 'csolImag_yyy', sexpression, '')
    if (trim(adjustl(sexpression)) .ne. '') then
      call fparser_parseFunction(rfparser, 'csolImag_yyy', sexpression, (/'x','y'/), icomp=csolImag_yyy)
    end if

    ! Parse expression for bathymetry profile
    call parlst_getvalue_string(rparlist, sparam, 'cbathymetry', sexpression, '')
    if (trim(adjustl(sexpression)) .ne. '') then
      call fparser_parseFunction(rfparser, 'cbathymetry', sexpression, (/'x','y'/), icomp=cbathymetry)
    end if

    ! Parse expression for bottom stress
    call parlst_getvalue_string(rparlist, sparam, 'cstress', sexpression, '')
    if (trim(adjustl(sexpression)) .ne. '') then
      call fparser_parseFunction(rfparser, 'cstress', sexpression, (/'x','y'/), icomp=cstress)
    end if

    ! Parse expression for vertical eddy viscosity
    call parlst_getvalue_string(rparlist, sparam, 'cviscosity', sexpression, '')
    if (trim(adjustl(sexpression)) .ne. '') then
      call fparser_parseFunction(rfparser, 'cviscosity', sexpression, (/'x','y'/), icomp=cviscosity)
    end if
    
  end subroutine sse_initParamSSE

  ! ***************************************************************************

!<subroutine>

  subroutine sse_infoSSE(cproblemType)

!<description>   
    ! This subroutine outputs information about the global parameters
    ! of the SSE problem
!</description>

!<input>
    ! Problem type
    integer, intent(in) :: cproblemType
!</input>
!</subroutine>

    select case (cproblemType)
    case (SSE_SCALAR)
      call output_line('PROBLEMTYPE........: SSE problem')
      call output_line('FORMULATION........: second-order equation')
    case (SSE_SYSTEM1)
      call output_line('PROBLEMTYPE........: SSE problem')
      call output_line('FORMULATION........: first-order formulation No.1')
    case (SSE_SYSTEM2)
      call output_line('PROBLEMTYPE........: SSE problem')
      call output_line('FORMULATION........: first-order formulation No.2')
    case default
      call output_line("Invalid problem type", &
          OU_CLASS_ERROR,OU_MODE_STD,"sse_infoSSE")
      call sys_halt()
    end select

    call output_line('B..................: '//trim(adjustl(sys_sdE(dwidth,5))))
    call output_line('H0/H...............: '//trim(adjustl(sys_sdE(dheightRatio,5))))
    call output_line('H0.................: '//trim(adjustl(sys_sdE(dheight0,5))))
    call output_line('H..................: '//trim(adjustl(sys_sdE(dheight,5))))
    call output_line('L..................: '//trim(adjustl(sys_sdE(dlength,5))))
    call output_line('LB.................: '//trim(adjustl(sys_sdE(dlengthB,5))))
    call output_line('M2.................: '//trim(adjustl(sys_sdE(dforcing,5))))
    call output_line('F..................: '//trim(adjustl(sys_sdE(dcoraccel,5))))
    call output_line('G..................: '//trim(adjustl(sys_sdE(dgravaccel,5))))
    call output_line('OMEGA..............: '//trim(adjustl(sys_sdE(dtidalfreq,5))))
    call output_line('CBATHYMETRY........: '//trim(adjustl(sys_si(cbathymetry,3))))
    call output_line('CSTRESS............: '//trim(adjustl(sys_si(cstress,3))))
    call output_line('CVISCOSITY.........: '//trim(adjustl(sys_si(cviscosity,3))))
    call output_line('CSOLREAL...........: '//trim(adjustl(sys_si(csolReal,3))))
    call output_line('CSOLIMAG...........: '//trim(adjustl(sys_si(csolImag,3))))
    call output_line('CSOLREAL_X.........: '//trim(adjustl(sys_si(csolReal_x,3))))
    call output_line('CSOLIMAG_X.........: '//trim(adjustl(sys_si(csolImag_x,3))))
    call output_line('CSOLREAL_Y.........: '//trim(adjustl(sys_si(csolReal_y,3))))
    call output_line('CSOLIMAG_Y.........: '//trim(adjustl(sys_si(csolImag_y,3))))
    call output_line('CSOLREAL_XX........: '//trim(adjustl(sys_si(csolReal_xx,3))))
    call output_line('CSOLIMAG_XX........: '//trim(adjustl(sys_si(csolImag_xx,3))))
    call output_line('CSOLREAL_XY........: '//trim(adjustl(sys_si(csolReal_xy,3))))
    call output_line('CSOLIMAG_XY........: '//trim(adjustl(sys_si(csolImag_xy,3))))
    call output_line('CSOLREAL_YY........: '//trim(adjustl(sys_si(csolReal_yy,3))))
    call output_line('CSOLIMAG_YY........: '//trim(adjustl(sys_si(csolImag_yy,3))))
    call output_line('CSOLREAL_XXX.......: '//trim(adjustl(sys_si(csolReal_xxx,3))))
    call output_line('CSOLIMAG_XXX.......: '//trim(adjustl(sys_si(csolImag_xxx,3))))
    call output_line('CSOLREAL_XXY.......: '//trim(adjustl(sys_si(csolReal_xxy,3))))
    call output_line('CSOLIMAG_XXY.......: '//trim(adjustl(sys_si(csolImag_xxy,3))))
    call output_line('CSOLREAL_XYY.......: '//trim(adjustl(sys_si(csolReal_xyy,3))))
    call output_line('CSOLIMAG_XYY.......: '//trim(adjustl(sys_si(csolImag_xyy,3))))
    call output_line('CSOLREAL_YYY.......: '//trim(adjustl(sys_si(csolReal_yyy,3))))
    call output_line('CSOLIMAG_YYY.......: '//trim(adjustl(sys_si(csolImag_yyy,3))))

  end subroutine sse_infoSSE

end module sse_base_sse
