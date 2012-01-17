!##############################################################################
!# ****************************************************************************
!# <name> uuid </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains routines to handle Universally Unique Identifier.
!#
!# Citation from "http://en.wikipedia.org/wiki/UUID":
!#
!#    A Universally Unique Identifier (UUID) is an identifier standard used
!#    in software construction, standardised by the Open Software Foundation
!#    (OSF) as part of the Distributed Computing Environment (DCE). The intent
!#    of UUIDs is to enable distributed systems to uniquely identify information
!#    without significant central coordination. Thus, anyone can create a UUID
!#    and use it to identify something with reasonable confidence that the
!#    identifier will never be unintentionally used by anyone for anything else.
!#    Information labeled with UUIDs can therefore be later combined into a single
!#    database without needing to resolve name conflicts. The most widespread use
!#    of this standard is in Microsoft`s Globally Unique Identifiers (GUIDs).
!#    Other significant users include Linux`s ext2/ext3 filesystem, LUKS encrypted
!#    partitions, GNOME, KDE, and Mac OS X, all of which use implementations
!#    derived from the uuid library found in the e2fsprogs package.
!#
!#    A UUID is a 16-byte (128-bit) number. The number of theoretically
!#    possible UUIDs is therefore $2^{16*8} = 2^128 = 256^16$ or about $3.4*10^38$.
!#    This means that 1 trillion UUIDs would have to be created
!#    every nanosecond for 10 billion years to exhaust the number of
!#    UUIDs.
!#
!#    In its canonical form, a UUID consists of 32 hexadecimal digits,
!#    displayed in 5 groups separated by hyphens, in the form 8-4-4-4-
!#    12 for a total of 36 characters. For example:
!#
!# <verb>
!#        550e8400-e29b-41d4-a716-446655440000
!# </verb>
!#
!#    A UUID may also be used with a specific identifier intentionally
!#    used repeatedly to identify the same thing in different contexts.
!#    For example, in Microsoft`s Component Object Model, every
!#    component must implement the IUnknown interface, which is done by
!#    creating a UUID representing IUnknown. In all cases wherever
!#    IUnknown is used, whether it is being used by a process trying to
!#    access the IUnknown interface in a component, or by a component
!#    implementing the IUnknown interface, it is always referenced by
!#    the same identifier: 00000000-0000-0000-C000-000000000046.
!#
!#    <bf> Version 1 (MAC address) </bf>
!#
!#    Conceptually, the original (version 1) generation scheme for
!#    UUIDs was to concatenate the UUID version with the MAC address of
!#    the computer that is generating the UUID, and with the number of
!#    100-nanosecond intervals since the adoption of the Gregorian
!#    calendar. In practice, the actual algorithm is more complicated.
!#    This scheme has been criticised in that it is not sufficiently
!#    'opaque'; it reveals both the identity of the computer that
!#    generated the UUID and the time at which it did so.
!#
!#    <bf> Version 2 (DCE Security) </bf>
!#
!#    Version 2 UUIDs are similar to Version 1 UUIDs, with the upper
!#    byte of the clock sequence replaced by the identifier for a
!#    "local domain" (typically either the "POSIX UID domain" or the
!#    "POSIX GID domain") and the first 4 bytes of the timestamp
!#    replaced by the user`s POSIX UID or GID (with the "local domain"
!#    identifier indicating which it is).[1][2]
!#
!#    <bf> Version 3 (MD5 hash) </bf>
!#
!#    Version 3 UUIDs use a scheme deriving a UUID via MD5 from a URL,
!#    a fully qualified domain name, an Object identifier, a
!#    distinguished name (DN as used in Lightweight Directory Access
!#    Protocol), or on names in unspecified namespaces. Version 3 UUIDs
!#    have the form xxxxxxxx-xxxx-3xxx-xxxx-xxxxxxxxxxxx with
!#    hexadecimal digits x.
!#
!#    To determine the version 3 UUID of a given name the UUID of the
!#    namespace, e.g. 6ba7b810-9dad-11d1-80b4-00c04fd430c8 for a
!#    domain, is transformed to a string of bytes corresponding to its
!#    hexadecimal digits, concatenated with the input name, hashed with
!#    MD5 yielding 128 bits. Six bits are replaced by fixed values,
!#    four of these bits indicate the version, 0011 for version 3.
!#    Finally the fixed hash is transformed back into the hexadecimal
!#    form with hyphens separating the parts relevant in other UUID
!#    versions.
!#
!#    <bf> Version 4 (random) </bf>
!#
!#    Version 4 UUIDs use a scheme relying only on random numbers. This
!#    algorithm sets the version number, as well as two reserved bits.
!#    All other bits are set using a random or pseudorandom data
!#    source.
!#
!#    <bf> Version 5 (SHA-1 hash) </bf>
!#
!#    Version 5 UUIDs use a scheme with SHA-1 hashing, otherwise it is
!#    the same idea as in version 3. RFC 4122 states that version 5 is
!#    preferred over version 3 name based UUIDs.
!#
!# The following routines can be found here:
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

!$use omp_lib
  use fsystem
  use genoutput
  
  implicit none
  
  private

!<types>

!<typeblock>
  
  ! Type for universally unique identifier
  type t_uuid
    ! the UUID data
    integer, dimension(16) :: data = 0
  end type t_uuid
  
  public :: t_uuid

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
  
  public :: operator(.eq.)

  interface operator(.gt.)
    module procedure uuid_isGreater
  end interface
  
  public :: operator(.gt.)

  interface operator(.ge.)
    module procedure uuid_isGreaterEqual
  end interface
  
  public :: operator(.ge.)

  interface operator(.lt.)
    module procedure uuid_isSmaller
  end interface
  
  public :: operator(.lt.)
  
  interface operator(.le.)
    module procedure uuid_isSmallerEqual
  end interface
  
  public :: operator(.le.)

  interface operator(.ne.)
    module procedure uuid_isNotEqual
  end interface
  
  public :: operator(.ne.)

  public :: uuid_createUUID, uuid_createUUID_directly,uuid_createUUID_indirectly
  public :: uuid_isEqual
  public :: uuid_isNil
  public :: uuid_isGreater
  public :: uuid_isGreaterEqual
  public :: uuid_isSmaller
  public :: uuid_isSmallerEqual
  public :: uuid_conv2String
  
contains

!************************************************************************

!<subroutine>

  subroutine uuid_createUUID_directly(cversion, ruuid)

!<description>
    ! This subroutine creates a universally unique identifier directly
!</description>

!<input>
    ! version for the UUID.
    ! =0: Create zero-UUID.
    integer, intent(in) :: cversion
!</input>

!<output>
    ! the UUID
    type(t_uuid), intent(out) :: ruuid
!</output>

!</subroutine>

    ! local variables
    real(DP), dimension(16) :: Drandom
    integer :: i
    
    select case(cversion)
    case (0)
      ! nil UUID
      ruuid%data     = 0
      
    case (4)
      ! Pseudo-random UUID following the algorithm outlined on pages 10 and 11 of the
      ! "UUIDs and GUIDs" Internet Draft found at
      ! http://www.ics.uci.edu/pub/ietf/webdav/uuid-guid/draft-leach-uuids-guids-01.txt
      
      ! Generate 16 random double values
      call random_number(Drandom)
      
      do i = 1, 16
        ruuid%data(i) = int(128*Drandom(i))
      end do
      
      ! Include version: Take 7th byte and perform "and" operation with 0x0f
      !                  followed by an "or" operation with 0x40
      ruuid%data(7) = iand(ruuid%data(7), 15)
      ruuid%data(7) = ior(ruuid%data(7), 64)
      
      ! Include variant: Take 9th byte and perform "and" operation with 0x3f
      !                  followed by an "or" operation with 0x80
      ruuid%data(9) = iand(ruuid%data(9), 63)
      ruuid%data(9) = ior(ruuid%data(9), 128)
            
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
    
    bisEqual = .true.
    
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
