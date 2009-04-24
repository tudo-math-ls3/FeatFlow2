!##############################################################################
!# ****************************************************************************
!# <name> zpinch_basic </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains the basic routines and provides all global variables
!# which are required to solve the magneto hydrodynamic equations
!#
!# </purpose>
!##############################################################################

module zpinch_basic

  use transport_basic
  use euler_basic

  implicit none

  private
  public :: t_Zpinch


!<types>
  
!<typeblock>

  ! This structure contains all required data to describe an instance
  ! of the simplified magneto hydrodynamic flow benchmark application
  type t_Zpinch

    ! Application descriptor for the compressible Euler model
    type(t_euler) :: rappDescrEuler

    ! Application descriptor for the scalar transport model
    type(t_transport) :: rappDescrTransport
    
  end type t_Zpinch

!</typeblock>

!</types>

end module zpinch_basic
