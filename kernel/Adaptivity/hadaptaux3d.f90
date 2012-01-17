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

module hadaptaux3d

!$use omp_lib
  use collection
  use fsystem
  use hadaptaux

  implicit none

  private

end module hadaptaux3d
