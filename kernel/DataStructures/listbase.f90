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

  implicit none

!<constants>
!<constantblock description="Global flags for list ordering">

  ! Identifier for unordered list
  integer, parameter :: LIST_UNORDERED  = 0
  
  ! Identifier for increasingly ordered list
  integer, parameter :: LIST_INCREASING = 1

  ! Identifier for decreasingly ordered list
  integer, parameter :: LIST_DECREASING = 2

  ! Identifier for ordered list
  integer, parameter :: LIST_CSR7       = 3

!</constantblock>

!<constantblock description="Global flags for list likn-type">

  ! Identifier fir single-linked list
  integer, parameter :: LIST_SINGLELINKED = 1
  
  ! Identifier for double-linked list
  integer, parameter :: LIST_DOUBLELINKED = 2

!</constantblock>

!<constantblock description="Global flags for list operations">

  ! Identifier for "not found in list"
  integer, parameter :: LIST_NOT_FOUND = -1

  ! Identifier for "found in list"
  integer, parameter :: LIST_FOUND     =  0
  
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
