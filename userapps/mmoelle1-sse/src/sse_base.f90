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

  private
  public :: sse_bottomProfile
  public :: sse_bottomStress
  public :: sse_EddyViscosity

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

#ifdef CASE_ALEX
!<constantblock description="Problem parameters for Alex benchmark">

  ! Length of the channel
  real(DP), parameter, public :: dlength = 6.0E4_DP
  
  ! Convergent length of the channel
  real(DP), parameter, public :: dlengthB = 1.0e20_DP

  ! Width of the entrance of the channel
  real(DP), parameter, public :: dwidth = 15000.0_DP

  ! Type of the width of the channel
  integer, parameter, public :: iwidthType = 2

  ! Type of the bed profile
  integer, parameter, public :: ibathymetryType = 0

  ! Mean depth of the channel
  real(DP), parameter, public :: dheight  = 7.7_DP
  real(DP), parameter, public :: dheight0 = 7.7_DP

  ! Depth at the end in the scaled domain i.e. a = dheight0/dheight
  real(DP), parameter, public :: dheightRatio  = 0.0_DP

  ! Constant forcing at ioen boundary
  real(DP), parameter, public :: dforcing = 1.0_DP

  ! Coriolis acceleration
  real(DP), parameter, public :: dcoraccel  = 0.0_DP
  
  ! Frequency of the tidal constituent
  real(DP), parameter, public :: dtidalfreq = 1.454441043328608E-4_DP

  ! Type of the bottom stress
  integer, parameter, public :: istressType = 1

  ! Bottom stress
  real(DP), parameter, private :: dstress = 0.00000098_DP

  ! Type of the vertical eddy viscosity
  integer, parameter, public :: iviscosityType = 1

  ! Vertical eddy viscosity
  real(DP), parameter, public :: dviscosity = 0.019_DP

!</constantblock>
#else
#ifdef CASE_MARCHI
!<constantblock description="Problem parameters for Marchi benchmark">

  ! Length of the channel
  real(DP), parameter, public :: dlength = 1.0_DP
  
  ! Convergent length of the channel
  real(DP), parameter, public :: dlengthB = 0.0_DP

  ! Width of the entrance of the channel
  real(DP), parameter, public :: dwidth = 1.0_DP

  ! Type of the width of the channel
  integer, parameter, public :: iwidthType = 1

  ! Type of the bed profile
  integer, parameter, public :: ibathymetryType = 0

  ! Mean depth of the channel
  real(DP), parameter, public :: dheight  = 0.0_DP
  real(DP), parameter, public :: dheight0 = 0.0_DP

  ! Depth at the end in the scaled domain i.e. a = dheight0/dheight
  real(DP), parameter, public :: dheightRatio  = 0.0_DP

  ! Constant forcing at ioen boundary
  real(DP), parameter, public :: dforcing = 0.0_DP

  ! Coriolis acceleration
  real(DP), parameter, public :: dcoraccel  = 0.0_DP
  
  ! Frequency of the tidal constituent
  real(DP), parameter, public :: dtidalfreq = 0.0_DP

  ! Type of the bottom stress
  integer, parameter, public :: istressType = 1

  ! Bottom stress
  real(DP), parameter, private :: dstress = 0.0_DP

  ! Type of the vertical eddy viscosity
  integer, parameter, public :: iviscosityType = 1

  ! Vertical eddy viscosity
  real(DP), parameter, public :: dviscosity = 0.0_DP

!</constantblock>
#else
#ifdef CASE_WALTERS
!<constantblock description="Problem parameters for Walters benchmark">

  ! Length of the channel
  real(DP), parameter, public :: dlength = 4.0E3_DP
  
  ! Convergent length of the channel
  real(DP), parameter, public :: dlengthB = 0.0_DP

  ! Width of the entrance of the channel
  real(DP), parameter, public :: dwidth = 2.0E3_DP

  ! Type of the width of the channel
  integer, parameter, public :: iwidthType = 1

  ! Type of the bed profile
  integer, parameter, public :: ibathymetryType = 1

  ! Mean depth of the channel
  real(DP), parameter, public :: dheight  = 10.0_DP
  real(DP), parameter, public :: dheight0 = 0.0_DP

  ! Depth at the end in the scaled domain i.e. a = dheight0/dheight
  real(DP), parameter, public :: dheightRatio  = 0.0_DP

  ! Constant forcing at ioen boundary
  real(DP), parameter, public :: dforcing = 0.1_DP

  ! Coriolis acceleration
  real(DP), parameter, public :: dcoraccel  = 0.0_DP
  
  ! Frequency of the tidal constituent
  real(DP), parameter, public :: dtidalfreq = 0.00174532925199433_DP

  ! Type of the bottom stress
  integer, parameter, public :: istressType = 1

  ! Bottom stress
  real(DP), parameter, private :: dstress = 0.0_DP

  ! Type of the vertical eddy viscosity
  integer, parameter, public :: iviscosityType = 1

  ! Vertical eddy viscosity
  real(DP), parameter, public :: dviscosity = 0.012_DP

!</constantblock>
#else
#ifdef CASE_WINANT
!<constantblock description="Problem parameters for Winant benchmark">

  ! Length of the channel
  real(DP), parameter, public :: dlength = 6.0E4_DP
  
  ! Convergent length of the channel
  real(DP), parameter, public :: dlengthB = 0.0_DP

  ! Width of the entrance of the channel
  real(DP), parameter, public :: dwidth = 1000.0_DP

  ! Type of the width of the channel
  integer, parameter, public :: iwidthType = 1

  ! Type of the bed profile
  integer, parameter, public :: ibathymetryType = 2

  ! Mean depth of the channel
  real(DP), parameter, public :: dheight  = 10.0_DP
  real(DP), parameter, public :: dheight0 = 0.0_DP

  ! Depth at the end in the scaled domain i.e. a = dheight0/dheight
  real(DP), parameter, public :: dheightRatio  = 0.01_DP

  ! Constant forcing at ioen boundary
  real(DP), parameter, public :: dforcing = 1.0_DP

  ! Coriolis acceleration
  real(DP), parameter, public :: dcoraccel  = 0.0_DP
  
  ! Frequency of the tidal constituent
  real(DP), parameter, public :: dtidalfreq = 1.454441043328608E-4_DP

  ! Type of the bottom stress
  integer, parameter, public :: istressType = 1

  ! Bottom stress
  real(DP), parameter, private :: dstress = 1.0E20_DP

  ! Type of the vertical eddy viscosity
  integer, parameter, public :: iviscosityType = 1

  ! Vertical eddy viscosity
  real(DP), parameter, public :: dviscosity = 1.0E-3_DP

!</constantblock>
#else

#error 'Test case is undefined.'

#endif
#endif
#endif
#endif

!<constantblock description="Problem parameters for all benchmarks">

  ! Gravitational acceleration
  real(DP), parameter, public :: dgravaccel = 10.0_DP

  ! Constants r1 and r2
  real(DP), parameter, public :: dr1 = sqrt( cimg*(dcoraccel-dtidalfreq)/dviscosity)
  real(DP), parameter, public :: dr2 = sqrt(-cimg*(dcoraccel+dtidalfreq)/dviscosity)

!</constantblock>

!</constants>

contains

  ! ***************************************************************************

!<function>

  elemental function sse_bottomProfile(dx,dy) result(dh)

!<description>
    ! This function computes the bathymetry as a function of the
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

    if     (ibathymetryType .eq. 0) then

      ! linear bottom profile
      dh = dheight0 + (dheight-dheight0)*(1-dx/dlength)

    elseif (ibathymetryType .eq. 1) then

      ! constant bottom profile
      dh = dheight

    elseif (ibathymetryType .eq. 2) then

      ! parabolic bottom profile
      dh = dheight*(dheightRatio+(1.0_DP-dheightRatio)*1-(dy/dlengthB)**2)
    end if

  end function sse_bottomProfile

  ! ***************************************************************************

!<function>

  elemental function sse_bottomStress(dx,dy) result(ds)

!<description>
    ! This function computes the bottom stress as a function of the
    ! Cartesian coordinates (x,y)
!</description>

!<input>
    ! Cartesian coordinates
    real(DP), intent(in) :: dx,dy
!</input>

!<result>
    ! Bottom stress
    real(DP) :: ds
!</result>

!</function>

    ! local parameters
    real(DP), parameter :: dx1 = -15000.0_DP
    real(DP), parameter :: dx2 = -10000.0_DP
    real(DP), parameter :: ds1 =  0.1_DP
    real(DP), parameter :: ds2 =  0.00049_DP

    real(DP), parameter :: da = (ds1-ds2)/(dx1-dx2)
    real(DP), parameter :: db = (ds2*dx1-ds1*dx2)/(dx1-dx2)

    if     (istressType .eq. 0) then

      ! variable stress
      if (dx .le. dx1) then
        ds = ds1
      elseif (dx .ge. dx2) then
        ds = ds2
      else
        ds = da*dx+db
      end if
          
    elseif (istressType .eq. 1) then
      
      ! constant stress
      ds = dstress

    elseif (istressType .eq. 2) then

      ! stress proportional to bathymetry
      ds = dstress * sse_bottomProfile(dx,dy) / dheight

    end if

  end function sse_bottomStress

  ! ***************************************************************************

!<function>

  elemental function sse_EddyViscosity(dx,dy) result(dAv)

!<description>
    ! This function computes the vertical eddy viscosity as a function
    ! of the Cartesian coordinates (x,y)
!</description>

!<input>
    ! Cartesian coordinates
    real(DP), intent(in) :: dx,dy
!</input>

!<result>
    ! Vertical eddy viscosity
    real(DP) :: dAv
!</result>

!</function>

    ! local parameters
    real(DP), parameter :: dx1 = -15000.0_DP
    real(DP), parameter :: dx2 = -10000.0_DP
    real(DP), parameter :: dAv1 = 1.0_DP
    real(DP), parameter :: dAv2 = 0.0012_DP

    real(DP), parameter :: da = (dAv1-dAv2)/(dx1-dx2)
    real(DP), parameter :: db = (dAv2*dx1-dAv1*dx2)/(dx1-dx2)

    if     (iviscosityType .eq. 0) then

      ! variable eddy viscosity
      if (dx .le. dx1) then
        dAv = dAv1
      elseif (dx .ge. dx2) then
        dAv = dAv2
      else
        dAv = da*dx+db
      end if
          
    elseif (iviscosityType .eq. 1) then
      
      ! constant eddy viscosity
      dAv = dviscosity

    elseif (iviscosityType .eq. 2) then

      ! eddy viscosity proportional to bathymetry
      dAv = dviscosity * sse_bottomProfile(dx,dy) / dheight

    end if

  end function sse_EddyViscosity

end module sse_base
