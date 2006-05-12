!##############################################################################
!# ****************************************************************************
!# <name> derivatives </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains some definition to define derivative types.
!# </purpose>
!##############################################################################

!# Current version: $Id: element.f90,v 1.35 2006/03/09 14:55:43 grabo Exp $ 

MODULE derivatives

  USE fsystem

  IMPLICIT NONE
  
!<constants>
!<constantblock description="Descriptors to identify derivative types">

  ! function value in term
  INTEGER, PARAMETER :: DER_FUNC     = 1
  
  ! x derivative in term
  INTEGER, PARAMETER :: DER_DERIV_X  = 2

  ! y derivative in term
  INTEGER, PARAMETER :: DER_DERIV_Y  = 3

  ! 2nd x derivative in term
  INTEGER, PARAMETER :: DER_DERIV_XX = 4

  ! xy derivative in term
  INTEGER, PARAMETER :: DER_DERIV_XY = 5

  ! 2nd y derivative in term
  INTEGER, PARAMETER :: DER_DERIV_YY = 6

  ! Number of derivative identifiers
  INTEGER, PARAMETER :: DER_MAXNDER  = 6


!</constantblock>

!</constants>

END MODULE