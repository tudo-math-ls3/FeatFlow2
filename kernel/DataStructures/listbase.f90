!##############################################################################
!# ****************************************************************************
!# <name> listbase </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module file contains basic constants used for linked lists.
!#
!# </purpose>
!##############################################################################

module listbase

  use fsystem

  implicit none

!<constants>
!<constantblock description="Flags for the iterator specification bitfield">

  ! Reverse iterator
  integer(I32), parameter :: LIST_LSPEC_REVERSE = 2**0

  ! Iterator refers to virtual position. This flag is set if the key
  ! value is not found in an ordered list and the returned iterator
  ! refers to the element which follows the desired element in
  ! lexicographical order
  integer(I32), parameter :: LIST_LSPEC_VIRTUAL = 2**1

!</constantblock>

!<constantblock description="Internal tags for list status">
  
  ! Tag for empty list
  integer, parameter :: LNULL =  0
  
  ! Tag for head of list
  integer, parameter :: LHEAD = -2

  ! Tag for tail of list
  integer, parameter :: LTAIL = -1

  ! Tag for next free item
  integer, parameter :: LFREE =  0

!</constantblock>
!</constants>

end module listbase
