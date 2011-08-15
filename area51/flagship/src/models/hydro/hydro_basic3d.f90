!##############################################################################
!# ****************************************************************************
!# <name> hydro_basic3d </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains the basic routines required to solve the
!# compressible Euler/Navier-Stokes equations in 3D
!#
!# The following routines are available:
!#
!# 1.) hydro_getVarInterleaveFormat3d
!#     -> Extracts a single variable from the scalar vector of
!#        conservative variables stored in interleave format
!#
!# 2.) hydro_getVarBlockFormat3d
!#     -> Extracts a single variable from the vector of conservative
!#        variables stored in block format
!#
!# </purpose>
!##############################################################################

module hydro_basic3d

#define HYDRO_NDIM 3
#include "hydro.h"

  use fsystem
  use genoutput
  use linearsystemblock
  use linearsystemscalar

  implicit none

  private

  public :: hydro_getVarInterleaveFormat3d
  public :: hydro_getVarBlockFormat3d

contains

  !*****************************************************************************

!<subroutine>

  subroutine hydro_getVarInterleaveFormat3d(neq, nvar, cvariable, Ddata, Dvalue)

!<description>
    ! This subroutine extracs a single variable from the vector of
    ! conservative variabels which is stored in interleave format in 3D
!</description>

!<input>
    ! Number of equations
    integer, intent(in) :: neq

    ! Number of variables
    integer, intent(in) :: nvar

    ! Identifier for the variable
    character(LEN=*), intent(in) :: cvariable

    ! Vector of conservative variables
    real(DP), dimension(nvar,neq), intent(in) :: Ddata
!</input>

!<output>
    ! Extracted single variable
    real(DP), dimension(:), intent(out) :: Dvalue
!</output>
!</subroutine>

    ! local variables
    integer :: ieq

    
    if (trim(cvariable) .eq. 'density') then
      do ieq = 1, neq
        Dvalue(ieq) = DENSITY2(Ddata,IDX2_FORWARD,ieq,0,0)
      end do
      
    elseif (trim(cvariable) .eq. 'velocity_magnitude') then
      do ieq = 1, neq
        Dvalue(ieq) = VELMAGNITUDE2(Ddata,IDX2_FORWARD,ieq,0,0)
      end do

    elseif (trim(cvariable) .eq. 'velocity_x') then
      do ieq = 1, neq
        Dvalue(ieq) = XVELOCITY2(Ddata,IDX2_FORWARD,ieq,0,0)
      end do

    elseif (trim(cvariable) .eq. 'velocity_y') then
      do ieq = 1, neq
        Dvalue(ieq) = YVELOCITY2(Ddata,IDX2_FORWARD,ieq,0,0)
      end do

    elseif (trim(cvariable) .eq. 'velocity_z') then
      do ieq = 1, neq
        Dvalue(ieq) = ZVELOCITY2(Ddata,IDX2_FORWARD,ieq,0,0)
      end do

    elseif (trim(cvariable) .eq. 'momentum_x') then
      do ieq = 1, neq
        Dvalue(ieq) = XMOMENTUM2(Ddata,IDX2_FORWARD,ieq,0,0)
      end do

    elseif (trim(cvariable) .eq. 'momentum_y') then
      do ieq = 1, neq
        Dvalue(ieq) = YMOMENTUM2(Ddata,IDX2_FORWARD,ieq,0,0)
      end do

    elseif (trim(cvariable) .eq. 'momentum_z') then
      do ieq = 1, neq
        Dvalue(ieq) = ZMOMENTUM2(Ddata,IDX2_FORWARD,ieq,0,0)
      end do

    elseif (trim(cvariable) .eq. 'energy') then
      do ieq = 1, neq
        Dvalue(ieq) = SPECIFICTOTALENERGY2(Ddata,IDX2_FORWARD,ieq,0,0)
      end do
      
    elseif (trim(cvariable) .eq. 'total_energy') then
      do ieq = 1, neq
        Dvalue(ieq) = TOTALENERGY2(Ddata,IDX2_FORWARD,ieq,0,0)
      end do

    elseif (trim(cvariable) .eq. 'internal_energy') then
      do ieq = 1, neq
        Dvalue(ieq) = INTERNALENERGY2(Ddata,IDX2_FORWARD,ieq,0,0)
      end do
      
    elseif (trim(cvariable) .eq. 'kinetic_energy') then
      do ieq = 1, neq
        Dvalue(ieq) = KINETICENERGY2(Ddata,IDX2_FORWARD,ieq,0,0)
      end do

    elseif (trim(cvariable) .eq. 'pressure') then
      do ieq = 1, neq
        Dvalue(ieq) = PRESSURE2(Ddata,IDX2_FORWARD,ieq,0,0)
      end do

    elseif (trim(cvariable) .eq. 'machnumber') then
      do ieq = 1, neq
        Dvalue(ieq) = MACHNUMBER2(Ddata,IDX2_FORWARD,ieq,0,0)
      end do

    elseif (trim(cvariable) .eq. 'speedofsound') then
      do ieq = 1, neq
        Dvalue(ieq) = SOUNDSPEED2(Ddata,IDX2_FORWARD,ieq,0,0)
      end do
      
    else
      
      call output_line('Invalid variable name!',&
          OU_CLASS_ERROR,OU_MODE_STD,'hydro_getVarInterleaveFormat3d')
      call sys_halt()
      
    end if
    
  end subroutine hydro_getVarInterleaveFormat3d

  !*****************************************************************************

!<subroutine>

  subroutine hydro_getVarBlockformat3d(neq, nvar, cvariable, Ddata, Dvalue)

!<description>
    ! This subroutine extracs a single variable from the vector of
    ! conservative variabels which is stored in block format in 3D
!</description>

!<input>
    ! Number of equations
    integer, intent(in) :: neq

    ! Number of variables
    integer, intent(in) :: nvar

    ! Identifier for the variable
    character(LEN=*), intent(in) :: cvariable

    ! Vector of conservative variables
    real(DP), dimension(neq,nvar), intent(in) :: Ddata
!</input>

!<output>
    ! Extracted single variable
    real(DP), dimension(:), intent(out) :: Dvalue
!</output>
!</subroutine>

    ! local variables
    integer :: ieq


    if (trim(cvariable) .eq. 'density') then
      do ieq = 1, neq
        Dvalue(ieq) = DENSITY2(Ddata,IDX2_REVERSE,ieq,0,0)
      end do
      
    elseif (trim(cvariable) .eq. 'velocity_magnitude') then
      do ieq = 1, neq
        Dvalue(ieq) = VELMAGNITUDE2(Ddata,IDX2_REVERSE,ieq,0,0)
      end do

    elseif (trim(cvariable) .eq. 'velocity_x') then
      do ieq = 1, neq
        Dvalue(ieq) = XVELOCITY2(Ddata,IDX2_REVERSE,ieq,0,0)
      end do

    elseif (trim(cvariable) .eq. 'velocity_y') then
      do ieq = 1, neq
        Dvalue(ieq) = YVELOCITY2(Ddata,IDX2_REVERSE,ieq,0,0)
      end do

    elseif (trim(cvariable) .eq. 'velocity_z') then
      do ieq = 1, neq
        Dvalue(ieq) = ZVELOCITY2(Ddata,IDX2_REVERSE,ieq,0,0)
      end do

    elseif (trim(cvariable) .eq. 'momentum_x') then
      do ieq = 1, neq
        Dvalue(ieq) = XMOMENTUM2(Ddata,IDX2_REVERSE,ieq,0,0)
      end do

    elseif (trim(cvariable) .eq. 'momentum_y') then
      do ieq = 1, neq
        Dvalue(ieq) = YMOMENTUM2(Ddata,IDX2_REVERSE,ieq,0,0)
      end do

    elseif (trim(cvariable) .eq. 'momentum_z') then
      do ieq = 1, neq
        Dvalue(ieq) = ZMOMENTUM2(Ddata,IDX2_REVERSE,ieq,0,0)
      end do

    elseif (trim(cvariable) .eq. 'energy') then
      do ieq = 1, neq
        Dvalue(ieq) = SPECIFICTOTALENERGY2(Ddata,IDX2_REVERSE,ieq,0,0)
      end do
      
    elseif (trim(cvariable) .eq. 'total_energy') then
      do ieq = 1, neq
        Dvalue(ieq) = TOTALENERGY2(Ddata,IDX2_REVERSE,ieq,0,0)
      end do

    elseif (trim(cvariable) .eq. 'internal_energy') then
      do ieq = 1, neq
        Dvalue(ieq) = INTERNALENERGY2(Ddata,IDX2_REVERSE,ieq,0,0)
      end do
      
    elseif (trim(cvariable) .eq. 'kinetic_energy') then
      do ieq = 1, neq
        Dvalue(ieq) = KINETICENERGY2(Ddata,IDX2_REVERSE,ieq,0,0)
      end do

    elseif (trim(cvariable) .eq. 'pressure') then
      do ieq = 1, neq
        Dvalue(ieq) = PRESSURE2(Ddata,IDX2_REVERSE,ieq,0,0)
      end do

    elseif (trim(cvariable) .eq. 'machnumber') then
      do ieq = 1, neq
        Dvalue(ieq) = MACHNUMBER2(Ddata,IDX2_REVERSE,ieq,0,0)
      end do

    elseif (trim(cvariable) .eq. 'speedofsound') then
      do ieq = 1, neq
        Dvalue(ieq) = SOUNDSPEED2(Ddata,IDX2_REVERSE,ieq,0,0)
      end do
      
    else
      
      call output_line('Invalid variable name!',&
          OU_CLASS_ERROR,OU_MODE_STD,'hydro_getVarBlockformat3d')
      call sys_halt()
      
    end if    

  end subroutine hydro_getVarBlockformat3d

end module hydro_basic3d
