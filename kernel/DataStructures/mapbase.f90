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

!<constants>

!<constantblock description="Flags for the iterator specification bitfield">

  ! Reverse iterator
  integer(I32), parameter :: MAP_MSPEC_REVERSE = 2**0

  ! Iterator refers to element which already exists
  integer(I32), parameter :: MAP_MSPEC_EXISTS  = 2**1

!</constantblock>

!<constantblock description="Internal tags for map status">

  ! Tag for empty tree
  integer, parameter :: MNULL   = 0

  ! Tag for root of tree
  integer, parameter :: MROOT   = 0
  
  ! Tag for next free item
  integer, parameter :: MFREE   = 0

  ! Tag for left child node
  integer, parameter :: MLEFT   = 0

  ! Tag for right child node
  integer, parameter :: MRIGHT  = 1

!</constantblock>
!</constants>

end module mapbase
