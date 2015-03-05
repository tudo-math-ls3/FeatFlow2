!##############################################################################
!# ****************************************************************************
!# <name> sse_base </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module provides the basic routines and constants.
!# </purpose>
!##############################################################################

module sse_base

  use fsystem
  use genoutput

  implicit none

  private

  public sse_getNDIM
  public sse_getNVAR
  public sse_getSection
  
#ifdef USE_COMPILER_INTEL
  public :: sinh
  public :: cosh
#endif

!<constants>

!<constantblock description="Constants for problem types">

  ! Compute standard Poisson problem in 2D
  integer, parameter, public :: POISSON_SCALAR = 0

  ! Compute standard Poisson problem in 2D as first-order system
  integer, parameter, public :: POISSON_SYSTEM = 1
  
  ! Compute SSE solution from scalar problem in 2D
  integer, parameter, public :: SSE_SCALAR     = 2

  ! Compute SSE solution from first-order system ($\sigma=A\nabla N$) in 2D
  integer, parameter, public :: SSE_SYSTEM1    = 3

  ! Compute SSE solution from first-order system ($\sigma=\nabla N$) in 2D
  integer, parameter, public :: SSE_SYSTEM2    = 4

  ! Compute Corine`s problem in 1D
  integer, parameter, public :: CORINE_1D      = 5

  ! Compute Corine`s problem in 2D
  integer, parameter, public :: CORINE_2D      = 6
!</constantblock>

  
!<constantblock description="Constants for complex numbers">

  ! Real part of complex number
  complex(DP), parameter, public :: creal = cmplx(1.0_DP,0.0_DP)

  ! Imaginary part of complex number
  complex(DP), parameter, public :: cimg = cmplx(0.0_DP,1.0_DP)

!</constantblock>

!</constants>
  
contains

  ! ***************************************************************************
  
#ifdef USE_COMPILER_INTEL
!<function>

  elemental function sinh(cx)

!<description>
    ! Complex valued hyperbolic sine functions (available in Fortran 2008)
!</description>

!<input>
    complex(DP), intent(in) :: cx
!</input>

!<result>
    complex(DP) :: sinh
!</result>
!</function>

    sinh = -cmplx(0.0_DP,1.0_DP) * sin(cmplx(0.0_DP,1.0_DP)*cx) 

  end function
#endif
  
  ! ***************************************************************************
  
#ifdef USE_COMPILER_INTEL
!<function>

  elemental function cosh(cx)

!<description>
    ! Complex valued hyperbolic sine functions (available in Fortran 2008)
!</description>

!<input>
    complex(DP), intent(in) :: cx
!</input>

!<result>
    complex(DP) :: cosh
!</result>
!</function>

    cosh = cos(cmplx(0.0_DP,1.0_DP)*cx) 

  end function
#endif

  ! ***************************************************************************

!<function>

  function sse_getNDIM(cproblemtype) result(ndim)

!<description>
    ! This function returns the number of spatial dimensions
!</description>

!<input>
    ! Problem type
    integer, intent(in) :: cproblemtype
!</input>

!<result>
    ! Number of spatial dimensions
    integer :: ndim
!</result>
!</function>

    select case(cproblemtype)
    case(POISSON_SCALAR, POISSON_SYSTEM)
      ndim = 2

    case(SSE_SCALAR, SSE_SYSTEM1,SSE_SYSTEM2)
      ndim = 2

    case (CORINE_1D)
      ndim = 1

    case (CORINE_2D)
      ndim = 2

    case default
      ndim = 0
      call output_line("Invalid problem type", &
          OU_CLASS_ERROR,OU_MODE_STD,"sse_getNDIM")
      call sys_halt()
    end select
  end function sse_getNDIM
  
  ! ***************************************************************************

!<function>

  function sse_getNVAR(cproblemtype) result(nvar)

!<description>
    ! This function returns the number of variables
!</description>

!<input>
    ! Problem type
    integer, intent(in) :: cproblemtype
!</input>

!<result>
    ! Number of variables
    integer :: nvar
!</result>
!</function>

    select case(cproblemtype)
    case(POISSON_SCALAR)
      nvar = 1

    case(POISSON_SYSTEM)
      nvar = 1 + sse_getNDIM(cproblemtype)

    case(SSE_SCALAR)
      nvar = 1 * 2

    case(SSE_SYSTEM1,SSE_SYSTEM2)
      nvar = (1 + sse_getNDIM(cproblemtype)) * 2

    case (CORINE_1D)
      nvar = 6

    case (CORINE_2D)
      nvar = 0

    case default
      nvar = 0

      call output_line("Invalid problem type", &
          OU_CLASS_ERROR,OU_MODE_STD,"sse_getNVAR")
      call sys_halt()
    end select
  end function sse_getNVAR

  ! ***************************************************************************

!<function>

  function sse_getSection(cproblemtype) result(sstring)

!<description>
    ! This function returns the section name
!</description>

!<input>
    ! Problem type
    integer, intent(in) :: cproblemtype
!</input>

!<result>
    ! Section name
    character(SYS_STRLEN) :: sstring
!</result>
!</function>

    select case(cproblemtype)
    case(POISSON_SCALAR, POISSON_SYSTEM)
      sstring = 'POISSON'

    case(SSE_SCALAR, SSE_SYSTEM1,SSE_SYSTEM2)
      sstring = 'SSE'

    case (CORINE_1D,CORINE_2D)
      sstring = 'CORINE'

    case default
      sstring = ''

      call output_line("Invalid problem type", &
          OU_CLASS_ERROR,OU_MODE_STD,"sse_getSection")
      call sys_halt()
    end select
  end function sse_getSection
  
end module sse_base
