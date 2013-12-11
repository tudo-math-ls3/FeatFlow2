!##############################################################################
!# ****************************************************************************
!# <name> mapbase </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module file contains basic constants used for maps.
!#
!# </purpose>
!##############################################################################

module mapbase

  use fsystem

  implicit none

  private

!<constants>

!<constantblock description="Flags for the iterator specification bitfield">

  ! Reverse iterator
  integer(I32), parameter, public :: MAP_MSPEC_REVERSE = 2_I32**0

  ! Iterator refers to element which already exists
  integer(I32), parameter, public :: MAP_MSPEC_EXISTS  = 2_I32**1

!</constantblock>

!<constantblock description="Internal tags for map status">

  ! Tag for empty tree
  integer, parameter, public :: MNULL  = 0

  ! Tag for root of tree
  integer, parameter, public :: MROOT  = 0

  ! Tag for next free item
  integer, parameter, public :: MFREE  = 0

  ! Tag for left child node
  integer, parameter, public :: MLEFT  = 0

  ! Tag for right child node
  integer, parameter, public :: MRIGHT = 1

!</constantblock>
!</constants>

end module mapbase
