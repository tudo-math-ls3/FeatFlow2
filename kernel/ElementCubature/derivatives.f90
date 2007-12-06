!##############################################################################
!# ****************************************************************************
!# <name> derivatives </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains some definition to define derivative types.
!# Basically there are three different types of derivative quantifiers:
!# for 1D, 2D and 3D discretisations. The definitions are used e.g. to
!# specify which types of derivatives a Finite Element function should
!# calculate; see element.f90 for a deeper look how to work with these
!# constants!
!# </purpose>
!##############################################################################

MODULE derivatives

  USE fsystem

  IMPLICIT NONE
  
!<constants>

!<constantblock description="Descriptors to identify derivative types in 1D">

  ! function value in term
  INTEGER, PARAMETER :: DER_FUNC1D     = 1
  
  ! 1st derivative in term
  INTEGER, PARAMETER :: DER_DERIV1D_X  = 2

  ! 2nd derivative in term
  INTEGER, PARAMETER :: DER_DERIV1D_XX = 3

!</constantblock>

!<constantblock description="Descriptors to identify derivative types in 2D">

  ! function value in term
  INTEGER, PARAMETER :: DER_FUNC       = 1
  INTEGER, PARAMETER :: DER_FUNC2D     = 1
  
  ! x derivative in term
  INTEGER, PARAMETER :: DER_DERIV_X    = 2
  INTEGER, PARAMETER :: DER_DERIV2D_X  = 2

  ! y derivative in term
  INTEGER, PARAMETER :: DER_DERIV_Y    = 3
  INTEGER, PARAMETER :: DER_DERIV2D_Y  = 3

  ! 2nd x derivative in term
  INTEGER, PARAMETER :: DER_DERIV_XX   = 4
  INTEGER, PARAMETER :: DER_DERIV2D_XX = 4

  ! xy derivative in term
  INTEGER, PARAMETER :: DER_DERIV_XY   = 5
  INTEGER, PARAMETER :: DER_DERIV2D_XY = 5

  ! 2nd y derivative in term
  INTEGER, PARAMETER :: DER_DERIV_YY   = 6
  INTEGER, PARAMETER :: DER_DERIV2D_YY = 6

!</constantblock>

!<constantblock description="Descriptors to identify derivative types in 3D">

  ! function value in term
  INTEGER, PARAMETER :: DER_FUNC3D     = 1
  
  ! x derivative in term
  INTEGER, PARAMETER :: DER_DERIV3D_X  = 2

  ! y derivative in term
  INTEGER, PARAMETER :: DER_DERIV3D_Y  = 3

  ! z derivative in term
  INTEGER, PARAMETER :: DER_DERIV3D_Z  = 4

  ! 2nd x derivative in term
  INTEGER, PARAMETER :: DER_DERIV3D_XX = 5

  ! xy derivative in term
  INTEGER, PARAMETER :: DER_DERIV3D_XY = 6

  ! 2nd y derivative in term
  INTEGER, PARAMETER :: DER_DERIV3D_YY = 7

  ! xz derivative in term
  INTEGER, PARAMETER :: DER_DERIV3D_XZ = 8

  ! yz derivative in term
  INTEGER, PARAMETER :: DER_DERIV3D_YZ = 9

  ! zz derivative in term
  INTEGER, PARAMETER :: DER_DERIV3D_ZZ = 10

!</constantblock>

!<constantblock description="General constants.">

  ! Number of derivative identifiers.
  INTEGER, PARAMETER :: DER_MAXNDER  = 10

!</constantblock>

!</constants>

END MODULE