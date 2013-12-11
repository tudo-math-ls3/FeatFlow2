!##############################################################################
!# ****************************************************************************
!# <name> arraylistbase </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module file contains basic constants used for arraylists.
!#
!# </purpose>
!##############################################################################

module arraylistbase

  use fsystem

  implicit none

  private

!<constants>
!<constantblock description="Flags for the iterator specification bitfield">

  ! Reverse iterator
  integer(I32), parameter, public :: ALST_LSPEC_REVERSE = 2_I32**0

  ! Iterator refers to virtual position. This flag is set if the key
  ! value is not found in an ordered list and the returned iterator
  ! refers to the element which follows the desired element in
  ! lexicographical order
  integer(I32), parameter, public :: ALST_LSPEC_VIRTUAL = 2_I32**1

!</constantblock>

!<constantblock description="Internal tags for arraylist status">

  ! Tag for empty arraylist
  integer, parameter, public :: ALST_NULL = 0

  ! Tag for next free position in storage of arraylist
  integer, parameter, public :: ALST_FREE = 0

  ! Tag for head of each list
  integer, parameter, public :: ALST_HEAD = 1

  ! Tag for tail of each list
  integer, parameter, public :: ALST_TAIL = 2

  ! Tag for number of entries stored in each list
  integer, parameter, public :: ALST_NA   = 3

!</constantblock>
!</constants>

end module arraylistbase
