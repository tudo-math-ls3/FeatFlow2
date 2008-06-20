!##############################################################################
!# ****************************************************************************
!# <name> hadaptaux3d </name>
!# ****************************************************************************
!#
!# <purpose>
!#
!# WARNING: Do not USE this module in your applications unless you really
!#          know what you are doing. This module does no error checking!!!
!#
!# This module contains all auxiliary routines which are required for
!# performing h-adaptivity in 3D. Unlike other modules, all subroutines
!# are declared PUBLIC since they are used by module HADAPTIVITY.
!#
!# The following routines are available:
!#
!# </purpose>
!##############################################################################

MODULE hadaptaux3d

  USE collection
  USE fsystem
  USE hadaptaux
  USE hadaptaux2d

  IMPLICIT NONE

  PUBLIC

END MODULE hadaptaux3d
