!##############################################################################
!# ****************************************************************************
!# <name> uuid </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains routines to handle Universally Unique Identifier.
!#
!# The following routines can be found here
!#
!#  1.) uuid_createUUID = uuid_createUUID_directly /
!#                        uuid_createUUID_indirectly
!#      -> Create new UUID 
!#
!#  2.) uuid_isEqual
!#      -> Check if two UUIDs are equal
!#
!#  3.) uuid_isNil
!#      -> Check if UUID equals the nil UUID
!#
!#  4.) uuid_isGreater
!#      -> Check if one UUID is greater than another
!#
!#  5.) uuid_isGreaterEqual
!#      -> Check if one UUID is greater or equal than another
!#
!#  6.) uuid_isSmaller
!#      -> Check if one UUID is smaller than another
!#
!#  7.) uuid_isSmallerEqual
!#      -> Check if one UUID is smaller or equal than another
!#
!#  8.) uuid_isEqual
!#      -> Check if two UUIDs are equal
!#
!#  9.) uuid_conv2String
!#      -> Convert UUID to string representation
!#
!# </purpose>
!##############################################################################
module uuid

  use fsystem
  use genoutput
  
  implicit none
  

!<types>

!<typeblock>
  
  ! Type for universally unique identifier
  type t_uuid
    
    ! version of the UUID
    integer :: cversion

    ! the UUID data
    integer(1), dimension(16) :: data
  end type t_uuid

!</typeblock>

!</types>

!************************************************************************
!************************************************************************
!************************************************************************

  interface uuid_createUUID
    module procedure uuid_createUUID_directly
    module procedure uuid_createUUID_indirectly
  end interface
  
!************************************************************************
!************************************************************************
!************************************************************************

  interface operator(.eq.)
    module procedure uuid_isEqual
  end interface

  interface operator(.gt.)
    module procedure uuid_isGreater
  end interface

  interface operator(.ge.)
    module procedure uuid_isGreaterEqual
  end interface

  interface operator(.lt.)
    module procedure uuid_isSmaller
  end interface
  
  interface operator(.le.)
    module procedure uuid_isSmallerEqual
  end interface

  interface operator(.ne.)
    module procedure uuid_isNotEqual
  end interface

contains

!************************************************************************

!<subroutine>

  subroutine uuid_createUUID_directly(cversion, ruuid)

!<description>
    ! This subroutine creates a universally unique identifier directly
!</description>

!<input>
    ! version for the UUID
    integer, intent(in) :: cversion
!</input>

!<output>
    ! the UUID
    type(t_uuid), intent(out) :: ruuid
!</output>

!</subroutine>

    ! local variables
    integer :: i,byte0,byte
    
    select case(cversion)
    case (0)
      ! nil UUID
      ruuid%cversion = 0
      ruuid%data     = 0
      
    case (4)
      ! Pseudo-random UUID following the algorithm outlined on pages 10 and 11 of the
      ! "UUIDs and GUIDs" Internet Draft found at
      ! http://www-old.ics.uci.edu/pub/ietf/webdav/uuid-guid/draft-leach-uuids-guids-01.txt
      ruuid%cversion = 4
      
      ! Compute seed from system clock
      call system_clock(byte0)
      
      do i = 1, 16
        byte = MOD((57*byte0+1), 128)
        byte0 = MOD((57*byte+1), 128)
        ruuid%data(i) = byte0
      end do

      ! Include version: Take 7th byte and perform "and" operation with 0x0f
      !                  followed by an "or" operation with 0x40
      byte = ruuid%data(7)
      byte = iand(byte, 15)
      byte = ior(byte, 64)
      ruuid%data(7) = byte
      
      ! Include variant: Take 9th byte and perform "and" operation with 0x3f
      !                  followed by an "or" operation with 0x80
      byte = ruuid%data(9)
      byte = iand(byte, 63)
      byte = ior(byte, 128)
      ruuid%data(9) = byte
      
    case default
      call output_line('Unsupported UUID version!',&
          OU_CLASS_ERROR,OU_MODE_STD,'uuid_createUUID_directly')
      call sys_halt()
    end select
    
  end subroutine uuid_createUUID_directly

!************************************************************************

!<subroutine>

  subroutine uuid_createUUID_indirectly(svalue, ruuid)

!<description>
    ! This subroutine creates a universally unique identifier
    ! from a given string
!</description>

!<input>
    ! the string representation of the UUID
    character(len=36), intent(in) :: svalue
!</input>

!<output>
    ! the UUID
    type(t_uuid), intent(out)     :: ruuid
!</outptu>

!</subroutine>
    
    ! Set UUID data
    read(svalue( 1: 2), '(Z2)') ruuid%data(1)
    read(svalue( 3: 4), '(Z2)') ruuid%data(2)
    read(svalue( 5: 6), '(Z2)') ruuid%data(3)
    read(svalue( 7: 8), '(Z2)') ruuid%data(4)

    if (svalue(9:9) .ne. '-') then
      call output_line('String does not represent UUID!',&
          OU_CLASS_ERROR,OU_MODE_STD,'uuid_createUUID_indirectly')
      call sys_halt()
    end if

    read(svalue(10:11),'(Z2)') ruuid%data(5)
    read(svalue(12:13),'(Z2)') ruuid%data(6)

    if (svalue(14:14) .ne. '-') then
      call output_line('String does not represent UUID!',&
          OU_CLASS_ERROR,OU_MODE_STD,'uuid_createUUID_indirectly')
      call sys_halt()
    end if

    read(svalue(15:16),'(Z2)') ruuid%data(7)
    read(svalue(17:18),'(Z2)') ruuid%data(8)

    if (svalue(19:19) .ne. '-') then
      call output_line('String does not represent UUID!',&
          OU_CLASS_ERROR,OU_MODE_STD,'uuid_createUUID_indirectly')
      call sys_halt()
    end if

    read(svalue(20:21),'(Z2)') ruuid%data(9)
    read(svalue(22:23),'(Z2)') ruuid%data(10)

    if (svalue(24:24) .ne. '-') then
      call output_line('String does not represent UUID!',&
          OU_CLASS_ERROR,OU_MODE_STD,'uuid_createUUID_indirectly')
      call sys_halt()
    end if

    read(svalue(25:26),'(Z2)') ruuid%data(11)
    read(svalue(27:28),'(Z2)') ruuid%data(12)
    read(svalue(29:30),'(Z2)') ruuid%data(13)
    read(svalue(31:32),'(Z2)') ruuid%data(14)
    read(svalue(33:34),'(Z2)') ruuid%data(15)
    read(svalue(35:36),'(Z2)') ruuid%data(16)

    ! Set version
    ruuid%cversion = 4

  end subroutine uuid_createUUID_indirectly
  
!************************************************************************

!<function>

  function uuid_isEqual(ruuid1, ruuid2) result (bisEqual)

!<description>
    ! This function compares two UUIDs and returns .TRUE.
    ! if they are equal, otherwise .FALSE. is returned
!</description>

!<input>
    ! the first UUID
    type(t_uuid), intent(in) :: ruuid1

    ! the second UUID
    type(t_uuid), intent(in) :: ruuid2
!</input>

!<result>
    logical :: bisEqual
!</result>

!</function>

    ! local variable
    integer :: i
    
    bisEqual = ruuid1%cversion .eq. ruuid2%cversion
    
    do i = 1, 16
      bisEqual = bisEqual .and. (ruuid1%data(i) .eq. ruuid2%data(i))
    end do
  end function uuid_isEqual

!************************************************************************

!<function>

  function uuid_isNil(ruuid) result (bisNil)

!<description>
    ! This function checks if the UUID equals the nil UUID.
!</description>

!<input>
    ! the UUID
    type(t_uuid), intent(in) :: ruuid
!</input>

!<result>
    logical :: bisNil
!</result>

!</function>

    ! local variable
    type(t_uuid) :: ruuidTmp
    
    ruuidTmp%cversion = ruuid%cversion
    ruuidTmp%data = 0

    bisNil = uuid_isEqual(ruuid, ruuidTmp)

  end function uuid_isNil

!************************************************************************

!<function>

  function uuid_isGreater(ruuid1, ruuid2) result (bisGreater)

!<description>
    ! This function compares two UUIDs and returns .TRUE.
    ! if the first one is greater than the second one.
!</description>

!<input>
    ! the first UUID
    type(t_uuid), intent(in) :: ruuid1

    ! the second UUID
    type(t_uuid), intent(in) :: ruuid2
!</input>

!<result>
    logical :: bisGreater
!</result>

!</function>

    ! local variable
    integer :: i
    
    bisGreater = .true.

    do i = 1, 16
      bisGreater = bisGreater .and. (ruuid1%data(i) .gt. ruuid2%data(i))
    end do
  end function uuid_isGreater

!************************************************************************

!<function>

  function uuid_isGreaterEqual(ruuid1, ruuid2) result (bisGreaterEqual)

!<description>
    ! This function compares two UUIDs and returns .TRUE.
    ! if the first one is greater or equal than the second one.
!</description>

!<input>
    ! the first UUID
    type(t_uuid), intent(in) :: ruuid1

    ! the second UUID
    type(t_uuid), intent(in) :: ruuid2
!</input>

!<result>
    logical :: bisGreaterEqual
!</result>

!</function>

    ! local variable
    integer :: i
    
    bisGreaterEqual = .true.

    do i = 1, 16
      bisGreaterEqual = bisGreaterEqual .and. (ruuid1%data(i) .ge. ruuid2%data(i))
    end do
  end function uuid_isGreaterEqual

!************************************************************************

!<function>

  function uuid_isSmaller(ruuid1, ruuid2) result (bisSmaller)

!<description>
    ! This function compares two UUIDs and returns .TRUE.
    ! if the first one is smaller than the second one.
!</description>

!<input>
    ! the first UUID
    type(t_uuid), intent(in) :: ruuid1

    ! the second UUID
    type(t_uuid), intent(in) :: ruuid2
!</input>

!<result>
    logical :: bisSmaller
!</result>

!</function>

    ! local variable
    integer :: i
    
    bisSmaller = .true.

    do i = 1, 16
      bisSmaller = bisSmaller .and. (ruuid1%data(i) .lt. ruuid2%data(i))
    end do
  end function uuid_isSmaller

!************************************************************************

!<function>

  function uuid_isSmallerEqual(ruuid1, ruuid2) result (bisSmallerEqual)

!<description>
    ! This function compares two UUIDs and returns .TRUE.
    ! if the first one is smaller or equal than the second one.
!</description>

!<input>
    ! the first UUID
    type(t_uuid), intent(in) :: ruuid1

    ! the second UUID
    type(t_uuid), intent(in) :: ruuid2
!</input>

!<result>
    logical :: bisSmallerEqual
!</result>

!</function>

    ! local variable
    integer :: i
    
    bisSmallerEqual = .true.

    do i = 1, 16
      bisSmallerEqual = bisSmallerEqual .and. (ruuid1%data(i) .le. ruuid2%data(i))
    end do
  end function uuid_isSmallerEqual

!************************************************************************

!<function>

  function uuid_isNotEqual(ruuid1, ruuid2) result (bisNotEqual)

!<description>
    ! This function compares two UUIDs and returns .TRUE.
    ! if they are note equal, otherwise .FALSE. is returned
!</description>

!<input>
    ! the first UUID
    type(t_uuid), intent(in) :: ruuid1

    ! the second UUID
    type(t_uuid), intent(in) :: ruuid2
!</input>

!<result>
    logical :: bisNotEqual
!</result>

!</function>

    bisNotEqual = .not.uuid_isEqual(ruuid1, ruuid2)
  end function uuid_isNotEqual

!************************************************************************

!<function>

  function uuid_conv2String(ruuid) result (svalue)

!<description>
    ! This function converts a UUID to string representation
!</description>

!<input>
    ! the UUID
    type(t_uuid), intent(in) :: ruuid
!</input>

!<result>
    character(len=36) :: svalue
!</result>

!</function>

    write(svalue( 1: 2),'(Z2.2)') ruuid%data(1)
    write(svalue( 3: 4),'(Z2.2)') ruuid%data(2)
    write(svalue( 5: 6),'(Z2.2)') ruuid%data(3)
    write(svalue( 7: 8),'(Z2.2)') ruuid%data(4)
    
    svalue(9:9) = '-'

    write(svalue(10:11),'(Z2.2)') ruuid%data(5)
    write(svalue(12:13),'(Z2.2)') ruuid%data(6)

    svalue(14:14) = '-'

    write(svalue(15:16),'(Z2.2)') ruuid%data(7)
    write(svalue(17:18),'(Z2.2)') ruuid%data(8)

    svalue(19:19) = '-'

    write(svalue(20:21),'(Z2.2)') ruuid%data(9)
    write(svalue(22:23),'(Z2.2)') ruuid%data(10)

    svalue(24:24) = '-'

    write(svalue(25:26),'(Z2.2)') ruuid%data(11)
    write(svalue(27:28),'(Z2.2)') ruuid%data(12)
    write(svalue(29:30),'(Z2.2)') ruuid%data(13)
    write(svalue(31:32),'(Z2.2)') ruuid%data(14)
    write(svalue(33:34),'(Z2.2)') ruuid%data(15)
    write(svalue(35:36),'(Z2.2)') ruuid%data(16)

  end function uuid_conv2String
end module uuid
