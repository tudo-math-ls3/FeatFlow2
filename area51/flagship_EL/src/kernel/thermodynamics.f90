!##############################################################################
!# ****************************************************************************
!# <name> thermodynamics </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains some special conversion routines which are used
!# to handle thermodynamic applications in one, two and three dimensions.
!#
!# The following routines are available
!#
!#  1.) thdyn_speedofSound
!#      -> Calculates the speed of sound
!#
!#  2.) thdyn_Machnumber
!#      -> Calculates the Mach number
!#
!#  3.) thdyn_totalEnergy
!#      -> Calculates the total energy
!#
!#  4.) thdyn_pressure
!#      -> Calculates the pressure
!#
!# </purpose>
!##############################################################################
module thermodynamics
  
  use fsystem

  implicit none

  private

  public :: thdyn_speedOfSound
  public :: thdyn_Machnumber
  public :: thdyn_totalEnergy
  public :: thdyn_pressure

  ! *****************************************************************************

  interface thdyn_Machnumber
    module procedure thdyn_Machnumber1D
    module procedure thdyn_Machnumber2D
    module procedure thdyn_Machnumber3D
  end interface

  interface thdyn_totalEnergy
    module procedure thdyn_totalEnergy1D
    module procedure thdyn_totalEnergy2D
    module procedure thdyn_totalEnergy3D
  end interface

  interface thdyn_pressure
    module procedure thdyn_pressure1D
    module procedure thdyn_pressure2D
    module procedure thdyn_pressure3D
  end interface

  ! *****************************************************************************

!<constants>

!<constantblock description="Ratio of specific heats">
  
  real(DP), parameter, public :: GAMMA_AIR        = 1.4_DP

!</constantblock>


!<constantblock description="Specific heats at constant pressure and volume">
  
  real(DP), parameter, public :: CV_AIR           = 717_DP
  real(DP), parameter, public :: CP_AIR           = 1004_DP

!</constantblock>


!<constantblock description="Universal gas constants">

  real(DP), parameter, public :: RGAS_AIR         = 287_DP

!</constantblock>
!</constants>

contains

  !*****************************************************************************

!<function>
  
  elemental function thdyn_speedOfSound(gamma, p, rho) result(c)

!<description>
    ! This function computes the speed of sound
    ! $$c=\sqrt{\gamma p}{\rho}$$
!</description>

!<input>
    ! ratio of specific heats
    real(DP), intent(in) :: gamma
    
    ! pressure
    real(DP), intent(in) :: p

    ! density
    real(DP), intent(in) :: rho
!</input>

!<result>
    ! speed of sound
    real(DP) :: c
!</result>
!</function>

    c=sqrt(gamma*p/rho)
  end function thdyn_speedOfSound

  !*****************************************************************************

!<function>

  elemental function thdyn_Machnumber1D(gamma, p, rho, u) result(M)

!<description>
    ! This function computes the Mach number in 1D
    ! $$M=\sqrt{\frac{\rho}{\gamma p}}\|\{bf v}\|$$
!</description>

!<input>
    ! ratio of specific heats
    real(DP), intent(in) :: gamma
    
    ! pressure
    real(DP), intent(in) :: p

    ! density
    real(DP), intent(in) :: rho

    ! x-velocity component
    real(DP), intent(in) :: u
!</input>

!<result>
    ! Mach number
    real(DP) :: M
!</result>
!</function>

    M = sqrt(max((rho*u*u)/(GAMMA*p), SYS_EPSREAL))

  end function thdyn_Machnumber1D

  !*****************************************************************************

!<function>

  elemental function thdyn_Machnumber2D(gamma, p, rho, u, v) result(M)

!<description>
    ! This function computes the Mach number in 2D
    ! $$M=\sqrt{\frac{\rho}{\gamma p}}\|\{bf v}\|$$
!</description>

!<input>
    ! ratio of specific heats
    real(DP), intent(in) :: gamma
    
    ! pressure
    real(DP), intent(in) :: p

    ! density
    real(DP), intent(in) :: rho

    ! x-velocity component
    real(DP), intent(in) :: u
    
    ! y-velocity component
    real(DP), intent(in) :: v
!</input>

!<result>
    ! Mach number
    real(DP) :: M
!</result>
!</function>

    M = sqrt(max((rho*u*u+rho*v*v)/(GAMMA*p), SYS_EPSREAL))

  end function thdyn_Machnumber2D

  !*****************************************************************************

!<function>

  elemental function thdyn_Machnumber3D(gamma, p, rho, u, v, w) result(M)

!<description>
    ! This function computes the Mach number in 3D
    ! $$M=\sqrt{\frac{\rho}{\gamma p}}\|\{bf v}\|$$
!</description>

!<input>
    ! ratio of specific heats
    real(DP), intent(in) :: gamma
    
    ! pressure
    real(DP), intent(in) :: p

    ! density
    real(DP), intent(in) :: rho

    ! x-velocity component
    real(DP), intent(in) :: u
    
    ! y-velocity component
    real(DP), intent(in) :: v

    ! z-velocity component
    real(DP), intent(in) :: w
!</input>

!<result>
    ! Mach number
    real(DP) :: M
!</result>
!</function>

    M = sqrt(max((rho*u*u+rho*v*v+rho*w*w)/(GAMMA*p), SYS_EPSREAL))

  end function thdyn_Machnumber3D

  !*****************************************************************************

!<function>

  elemental function thdyn_totalEnergy1D(gamma, p, rho, u) result(E)

!<description>
    ! This function computes the total energy in 1D
    ! $$E=\frac{p}{\gamma(\gamma-1)}+\frac12\|{\bf v\|^2$$
!</description>

!<input>
    ! ratio of specific heats
    real(DP), intent(in) :: gamma

    ! pressure
    real(DP), intent(in) :: p

    ! density
    real(DP), intent(in) :: rho

    ! x-velocity component
    real(DP), intent(in) :: u
!</input>

!<result>
    ! total energy
    real(DP) :: E
!</result>
!</function>

    E = p/(rho*(gamma-1._DP))+0.5_DP*(u*u)
  end function thdyn_totalEnergy1D

  !*****************************************************************************

!<function>

  elemental function thdyn_totalEnergy2D(gamma, p, rho, u, v) result(E)

!<description>
    ! This function computes the total energy in 2D
    ! $$E=\frac{p}{\gamma(\gamma-1)}+\frac12\|{\bf v\|^2$$
!</description>

!<input>
    ! ratio of specific heats
    real(DP), intent(in) :: gamma

    ! pressure
    real(DP), intent(in) :: p

    ! density
    real(DP), intent(in) :: rho

    ! x-velocity component
    real(DP), intent(in) :: u

    ! y-velocity component
    real(DP), intent(in) :: v
!</input>

!<result>
    ! total energy
    real(DP) :: E
!</result>
!</function>

    E = p/(rho*(gamma-1._DP))+0.5_DP*(u*u+v*v)
  end function thdyn_totalEnergy2D

!*****************************************************************************

!<function>

  elemental function thdyn_totalEnergy3D(gamma, p, rho, u, v, w) result(E)

!<description>
    ! This function computes the total energy in 3D
    ! $$E=\frac{p}{\gamma(\gamma-1)}+\frac12\|{\bf v\|^2$$
!</description>

!<input>
    ! ratio of specific heats
    real(DP), intent(in) :: gamma

    ! pressure
    real(DP), intent(in) :: p

    ! density
    real(DP), intent(in) :: rho

    ! x-velocity component
    real(DP), intent(in) :: u

    ! y-velocity component
    real(DP), intent(in) :: v

    ! z-velocity component
    real(DP), intent(in) :: w
!</input>

!<result>
    ! total energy
    real(DP) :: E
!</result>
!</function>

    E = p/(rho*(gamma-1._DP))+0.5_DP*(u*u+v*v+w*w)
  end function thdyn_totalEnergy3D

  !*****************************************************************************

!<function>

  elemental function thdyn_pressure1D(gamma, E, rho, u) result(p)

!<description>
    ! This function computes the pressure in 1D
    ! $$p=(\gamma-1)\rho(E-\frac12\|{\bf v}\}^2)$$
!</description>

!<input>
    ! ratio of specific heats
    real(DP), intent(in) :: gamma

    ! total energy
    real(DP), intent(in) :: E

    ! density
    real(DP), intent(in) :: rho

    ! x-velocity component
    real(DP), intent(in) :: u
!</input>

!<result>
    ! pressure
    real(DP) :: p
!</result>
!</function>

    p = (gamma-1)*rho*(E-0.5_DP*(u*u))
  end function thdyn_pressure1D

  !*****************************************************************************

!<function>

  elemental function thdyn_pressure2D(gamma, E, rho, u, v) result(p)

!<description>
    ! This function computes the pressure in 2D
    ! $$p=(\gamma-1)\rho(E-\frac12\|{\bf v}\}^2)$$
!</description>

!<input>
    ! ratio of specific heats
    real(DP), intent(in) :: gamma

    ! total energy
    real(DP), intent(in) :: E

    ! density
    real(DP), intent(in) :: rho

    ! x-velocity component
    real(DP), intent(in) :: u

    ! y-velocity component
    real(DP), intent(in) :: v
!</input>

!<result>
    ! pressure
    real(DP) :: p
!</result>
!</function>

    p = (gamma-1)*rho*(E-0.5_DP*(u*u+v*v))
  end function thdyn_pressure2D

  !*****************************************************************************

!<function>

  elemental function thdyn_pressure3D(gamma, E, rho, u, v, w) result(p)

!<description>
    ! This function computes the pressure in 3D
    ! $$p=(\gamma-1)\rho(E-\frac12\|{\bf v}\}^2)$$
!</description>

!<input>
    ! ratio of specific heats
    real(DP), intent(in) :: gamma

    ! total energy
    real(DP), intent(in) :: E

    ! density
    real(DP), intent(in) :: rho

    ! x-velocity component
    real(DP), intent(in) :: u

    ! y-velocity component
    real(DP), intent(in) :: v

    ! z-velocity component
    real(DP), intent(in) :: w
!</input>

!<result>
    ! pressure
    real(DP) :: p
!</result>
!</function>

    p = (gamma-1)*rho*(E-0.5_DP*(u*u+v*v+w*w))
  end function thdyn_pressure3D
end module thermodynamics
