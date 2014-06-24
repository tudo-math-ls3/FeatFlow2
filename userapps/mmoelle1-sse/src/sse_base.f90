!##############################################################################
!# ****************************************************************************
!# <name> sse_base </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module provides the basic routines and constants for solving
!# the elliptic equation for sea surface elevation.
!# </purpose>
!##############################################################################

module sse_base

  use fsystem

  implicit none

  public :: sse_bottomProfile
  private

!<constants>

!<constantblock description="Constants for problem type">

  ! Comput standard Poisson problem
  integer, parameter, public :: POISSON_SCALAR = 0

  ! Comput standard Poisson problem as first-order system
  integer, parameter, public :: POISSON_SYSTEM = 1

  ! Compute SSE solution from scalar problem
  integer, parameter, public :: SSE_SCALAR     = 2

  ! Compute SSE solution from first-order system ($\sigma=A\nabla N$)
  integer, parameter, public :: SSE_SYSTEM1    = 3

  ! Compute SSE solution from first-order system ($\sigma=\nabla N$)
  integer, parameter, public :: SSE_SYSTEM2    = 4

!</constantblock>

!<constantblock description="Constants for complex numbers">

  ! Real part of complex number
  complex(DP), parameter, public :: creal = cmplx(1.0_DP,0.0_DP)

  ! Imaginary part of complex number
  complex(DP), parameter, public :: cimg = cmplx(0.0_DP,1.0_DP)

!</constantblock>

!<constantblock description="Problem parameters">

  ! Tidal frequency
  real(DP), parameter, public :: dtidalfreq = 1.454441043328608E-4_DP

  ! Gravitational acceleration
  real(DP), parameter, public :: dgravaccel = 10.0_DP

  ! Bottom stress
  real(DP), parameter, public :: dstress    = 0.049_DP

  ! Eddy viscosity
  real(DP), parameter, public :: dviscosity = 0.012_DP

  ! Coriolis acceleration
  real(DP), parameter, public :: dcoraccel  = 0.0_DP

  ! Constants r1 and r2
  complex(DP), parameter, public :: cr1 = sqrt( cimg*(dcoraccel-dtidalfreq)/dviscosity)
  complex(DP), parameter, public :: cr2 = sqrt(-cimg*(dcoraccel+dtidalfreq)/dviscosity)


!</constantblock>

!</constants>

contains

!<function>

  elemental function sse_bottomProfile(dx,dy) result(dh)

!<description>
    ! This function computes the bottom profile as a function of the
    ! Cartisian coordinates (x,y)
!</description>

!<input>
    ! Cartesian coordinates
    real(DP), intent(in) :: dx,dy
!</input>

!<result>
    ! Height of the bottom profile
    real(DP) :: dh
!</result>

!</function>

    dh = 12.2_DP

  end function sse_bottomProfile

end module sse_base
