!##############################################################################
!# ****************************************************************************
!# <name> mhd_basic </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains the basic routines and provides all global variables
!# which are required to solve the magneto hydrodynamic equations
!#
!# </purpose>
!##############################################################################

module mhd_basic

  use codire_basic
  use euler_basic

  implicit none

  private
  public :: t_mhdSimple


!<types>
  
!<typeblock>

  ! This structure contains all required data to describe an instance
  ! of the simplified magneto hydrodynamic flow benchmark application
  type t_mhdSimple

    ! Application descriptor for the compressible Euler model
    type(t_euler) :: rappDescrEuler

    ! Application descriptor for the scalar transport model
    type(t_codire) :: rappDescrTransport
    
  end type t_mhdSimple

!</typeblock>

!</types>

end module mhd_basic
