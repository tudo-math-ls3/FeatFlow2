!##############################################################################
!# ****************************************************************************
!# <name> numbersets </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains routines to maintain sets of numbers.
!#
!# The basic type here is t_directAccessIntSet which capsules a set of
!# integer numbers that can directly be accessed. Elements can be added and
!# removed from the set up to its maximum size, and it can be checked
!# if an element is in the set. 
!# Note: This structure is realised as a bitfield for memory efficient
!# storage of the set.
!#
!# The module contains the following subroutines:
!#
!# 1.) nsets_initDASet
!#     -> Initialises a direct-access set for integer numbers.
!#
!# 2.) nsets_doneDASet
!#     -> Releases a direct-access set for integer numbers.
!#
!# 3.) nsets_addToDASet
!#     -> Adds a number to a direct-access set for integer numbers.
!#
!# 4.) nsets_removeFromDASet
!#     -> Removes a number from a direct-access set for integer numbers.
!#
!# 5.) nsets_DASetContains
!#     -> Checks if a number is in a direct-access set for integer numbers.
!# 
!# 6.) nsets_clearDASet
!#     -> Removes all elements from a direct-access set for integer numbers.
!#
!# 7.) nsets_invertDASet
!#     -> Inverts a direct-access set for integer numbers.
!#
!# 8.) nsets_addAllDASet
!#     -> Adds all possible elements to a direct-access set for integer 
!#        numbers.
!# </purpose>
!##############################################################################

module numbersets

  use fsystem
  use storage
  
  implicit none
  
  private

!<types>

!<typeblock>

  ! A direct access set that allows to directly access a set of integer numbers.
  type t_directAccessIntSet
  
    ! A handle identifying the memory of the set.
    integer :: ihandle
    
    ! Maximum number of elements in the set.
    integer :: nmaxEntries
    
    ! Number of bits in an integer
    integer :: ibitsPerInteger
    
    ! Shift divisor to calculate the integer number from a position.
    integer :: ishiftDivisor
    
    ! Bitmask to calculate the bit position from a position
    integer :: ibitmask
    
    ! A pointer to the set data
    integer, dimension(:), pointer :: p_Idata
    
  end type
  
!</typeblock>

  public :: t_directAccessIntSet

!</types>

  public :: nsets_initDASet
  public :: nsets_doneDASet
  public :: nsets_addToDASet
  public :: nsets_removeFromDASet
  public :: nsets_DASetContains
  public :: nsets_clearDASet
  public :: nsets_invertDASet
  public :: nsets_addAllDASet

contains

  !****************************************************************************

!<subroutine>
  subroutine nsets_initDASet (rset,nelements)
!<description>
  ! Initialises a direct access set for at most icount elements.
!</description>

!<input>
  ! Maximum number of elements in the set.
  integer, intent(in) :: nelements
!</input>

!<output>
  ! The initialised set.
  type(t_directAccessIntSet), intent(out) :: rset
!</output>

!</subroutine>

    integer, parameter :: i = 0
    integer :: j
    
    ! Calculate the divisors for a quick access to the array.
    rset%ishiftDivisor = 1
    j = bit_size(i)
    do while (.not. btest(j,1))
      rset%ishiftDivisor = rset%ishiftDivisor + 1
      j = ishft(j,-1)
    end do
    rset%ibitmask = ishft(1,rset%ishiftDivisor)-1
    rset%ishiftDivisor = -rset%ishiftDivisor
    
    ! Allocate as much memory as necessary. Fill it with 0.
    call storage_new ('sets_initDASet', 'set', (nelements+bit_size(i)-1)/bit_size(i),&
        ST_INT, rset%ihandle, ST_NEWBLOCK_ZERO)
    call storage_getbase_int(rset%ihandle,rset%p_Idata)
    
    rset%nmaxEntries = nelements
  
  end subroutine
  
  !****************************************************************************

!<subroutine>
  subroutine nsets_doneDASet (rset)
!<description>
  ! Releases a direct access set..
!</description>

!<inputoutput>
  ! The set to release
  type(t_directAccessIntSet), intent(inout) :: rset
!</inputoutput>

!</subroutine>

    nullify(rset%p_Idata)
    call storage_free (rset%ihandle)

  end subroutine
  
  !****************************************************************************
  
!<subroutine>
  subroutine nsets_addToDASet (rset,ielement)
!<description>
  ! Adds an element to the set rset if it's not already part of it.
!</description>

!<inputoutput>
  ! The set to modify.
  type(t_directAccessIntSet), intent(inout) :: rset
!</inputoutput>

!<input>
  ! Number of the element to add to the set.
  integer, intent(in) :: ielement
!</input>

!</subroutine>

    ! Set the bit. 
    rset%p_Idata(1+ishft(ielement,rset%ishiftDivisor)) = &
        ibset(rset%p_Idata(1+ishft(ielement,rset%ishiftDivisor)),&
            iand(ielement,rset%ibitmask)-1)
  
  end subroutine
  
  !****************************************************************************

!<subroutine>
  subroutine nsets_removeFromDASet (rset,ielement)
!<description>
  ! Removes an element from the set it it's part of it.
!</description>

!<inputoutput>
  ! The set to modify.
  type(t_directAccessIntSet), intent(inout) :: rset
!</inputoutput>

!<input>
  ! Number of the element to remove from the set.
  integer, intent(in) :: ielement
!</input>
!</subroutine>
  
    rset%p_Idata(1+ishft(ielement,rset%ishiftDivisor)) = &
        ibclr(rset%p_Idata(1+ishft(ielement,rset%ishiftDivisor)),&
            iand(ielement,rset%ibitmask)-1)
  
  end subroutine
  
  !****************************************************************************

!<function>
  integer function nsets_DASetContains (rset,ielement)
!<description>
  ! Determines if a set contains the element ielement
!</description>

!<inputoutput>
  ! The set to check.
  type(t_directAccessIntSet), intent(in) :: rset
!</inputoutput>

!<input>
  ! Number of the element to check.
  integer, intent(in) :: ielement
!</input>

!<returns>
  ! 1, if the element is in the set, 0 otherwise.
!</returns>

!</function>

    nsets_DASetContains = ibits(&
        rset%p_Idata(1+ishft(ielement,rset%ishiftDivisor)),&
        iand(ielement,rset%ibitmask)-1,1)
  
  end function

  !****************************************************************************

!<subroutine>
  subroutine nsets_clearDASet (rset)
!<description>
  ! Removes all elements from a direct-access set.
!</description>

!<inputoutput>
  ! The set to modify.
  type(t_directAccessIntSet), intent(in) :: rset
!</inputoutput>

!</subroutine>

    call storage_clear(rset%ihandle)
  
  end subroutine

  !****************************************************************************

!<subroutine>
  subroutine nsets_invertDASet (rset)
!<description>
  ! Inverts a direct-access set. All elements in the set are removed,
  ! all elements not in the set are added.
!</description>

!<inputoutput>
  ! The set to modify.
  type(t_directAccessIntSet), intent(in) :: rset
!</inputoutput>

!</subroutine>

    integer :: i
    
    do i=1,size(rset%p_Idata)
      rset%p_Idata(i) = NOT(rset%p_Idata(i))
    end do
  
  end subroutine
  
  !****************************************************************************

!<subroutine>
  subroutine nsets_addAllDASet (rset)
!<description>
  ! Adds all possible elements to a a direct-access set.
!</description>

!<inputoutput>
  ! The set to modify.
  type(t_directAccessIntSet), intent(in) :: rset
!</inputoutput>

!</subroutine>

    integer :: i
    
    do i=1,size(rset%p_Idata)
      rset%p_Idata(i) = NOT(0)
    end do
  
  end subroutine

end module
