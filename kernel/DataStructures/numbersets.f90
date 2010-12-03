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
!#
!# 9.) nsets_nelementsInSet
!#     -> Returns the number of elements in the set.
!#
!# 10.) nsets_getElements
!#      -> Creates a list of all elements in the set.
!#
!# 11.) nsets_putElements
!#      -> Puts elements from a list into the set.
!# </purpose>
!##############################################################################

module numbersets

  use fsystem
  use genoutput
  use storage
  
  implicit none
  
  private

!<types>

!<typeblock>

  ! A direct access set that allows to directly access a set of integer numbers.
  type t_directAccessIntSet
  
    ! A handle identifying the memory of the set.
    integer :: ihandle
    
    ! Specifies whether the memory should automatically be deallocated.
    logical :: bautodealloc
    
    ! Maximum number of elements in the set.
    integer :: nmaxEntries
    
    ! Number of integers in the data array.
    integer :: nintegers
    
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
  public :: nsets_nelementsInSet
  public :: nsets_getElements
  public :: nsets_putElements

contains

  !****************************************************************************

!<subroutine>
  subroutine nsets_initDASet (rset,nelements,ihandle)
!<description>
  ! Initialises a direct access set for at most icount elements.
!</description>

!<input>
  ! Maximum number of elements in the set.
  integer, intent(in) :: nelements
  
  ! OPTIONAL: A handle to an integer array, large enough to hold the set.
  ! If not present, a new memory block is automatically allocated.
  ! If present, the memory block is assumed to belong to the caller and
  ! is not automatically released upon the destruction of the set structure.
  integer, intent(in), optional :: ihandle
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
    rset%ibitsPerInteger = bit_size(i)
    j = rset%ibitsPerInteger
    do while (.not. btest(j,1))
      rset%ishiftDivisor = rset%ishiftDivisor + 1
      j = ishft(j,-1)
    end do
    rset%ibitmask = ishft(1,rset%ishiftDivisor)-1
    rset%ishiftDivisor = -rset%ishiftDivisor
    rset%nintegers = (nelements+rset%ibitsPerInteger-1)/rset%ibitsPerInteger
    
    if (.not. present(ihandle)) then
      ! Allocate as much memory as necessary. Fill it with 0.
      call storage_new ('sets_initDASet', 'set', rset%nintegers,&
          ST_INT, rset%ihandle, ST_NEWBLOCK_ZERO)
      rset%bautodealloc = .true.
    else
      ! Check if the memory block is large enough.
      call storage_getsize(ihandle,j)
      if (j*rset%ibitsPerInteger .lt. nelements) then
        call output_line ('Specified memory block too small.', &
                          OU_CLASS_ERROR,OU_MODE_STD,'nsets_initDASet')
        call sys_halt()
      end if
      
      ! We use that memory block
      rset%bautodealloc = .true.
      rset%ihandle = ihandle
      call storage_clear (ihandle)
    end if
    
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
    
    ! Only deallocate if this array belongs to the structure.
    if (rset%bautodealloc) then
      call storage_free (rset%ihandle)
    end if

  end subroutine
  
  !****************************************************************************
  
!<subroutine>
  subroutine nsets_addToDASet (rset,ielement)
!<description>
  ! Adds an element to the set rset if it is not already part of it.
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
    rset%p_Idata(1+ishft(ielement-1,rset%ishiftDivisor)) = &
        ibset(rset%p_Idata(1+ishft(ielement-1,rset%ishiftDivisor)),&
            iand(ielement-1,rset%ibitmask))
  
  end subroutine
  
  !****************************************************************************

!<subroutine>
  subroutine nsets_removeFromDASet (rset,ielement)
!<description>
  ! Removes an element from the set it it is part of it.
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
  
    rset%p_Idata(1+ishft(ielement-1,rset%ishiftDivisor)) = &
        ibclr(rset%p_Idata(1+ishft(ielement-1,rset%ishiftDivisor)),&
            iand(ielement-1,rset%ibitmask))
  
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
        rset%p_Idata(1+ishft(ielement-1,rset%ishiftDivisor)),&
        iand(ielement-1,rset%ibitmask),1)
  
  end function

  !****************************************************************************

!<subroutine>
  subroutine nsets_clearDASet (rset)
!<description>
  ! Removes all elements from a direct-access set.
!</description>

!<inputoutput>
  ! The set to modify.
  type(t_directAccessIntSet), intent(inout) :: rset
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
  type(t_directAccessIntSet), intent(inout) :: rset
!</inputoutput>

!</subroutine>

    integer :: i
    
    do i=1,size(rset%p_Idata)
      rset%p_Idata(i) = not(rset%p_Idata(i))
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
  type(t_directAccessIntSet), intent(inout) :: rset
!</inputoutput>

!</subroutine>

    integer :: i
    
    do i=1,size(rset%p_Idata)
      rset%p_Idata(i) = not(0)
    end do
  
  end subroutine

  !****************************************************************************

!<function>
  integer function nsets_nelementsInSet (rset)
!<description>
  ! Determines the number of elements in the set.
!</description>

!<input>
  ! The set to check.
  type(t_directAccessIntSet), intent(in) :: rset
!</input>

!<result>
  ! Number of elements in the set.
!</result>

!</function>
    integer :: i,j,k,ncount

    if (rset%nmaxEntries .eq. rset%ibitsPerInteger*rset%nintegers) then
    
      ! Loop through all integers
      ncount = 0
      do i = 1,rset%nintegers
        ! Get the integer
        j = rset%p_Idata(i)
        
        ! Shift it and count every one.
        do k=1,rset%ibitsPerInteger
          ncount = ncount + iand(j,1)
          j = ishft(j,-1)
        end do
      end do
      
    else
    
      ! Loop through all integers except the last one.
      ncount = 0
      do i = 1,rset%nintegers-1
        ! Get the integer
        j = rset%p_Idata(i)
        
        ! Shift it and count every one.
        do k=1,rset%ibitsPerInteger
          ncount = ncount + iand(j,1)
          j = ishft(j,-1)
        end do
      end do
      
      ! Get the last integer
      j = rset%p_Idata(i)
      
      ! Shift it and count every one.
      do k=1,iand(rset%nmaxEntries,rset%ibitmask)
        ncount = ncount + iand(j,1)
        j = ishft(j,-1)
      end do
    end if
    
    nsets_nelementsInSet = ncount
  
  end function

  !****************************************************************************

!<subroutine>
  subroutine nsets_getElements (rset,Ilist,nelementsInList)
!<description>
  ! Creates a list of all elements in the set.
!</description>

!<input>
  ! The set to check.
  type(t_directAccessIntSet), intent(in) :: rset
!</input>

!<output>
  ! An array receiving the elements in the set.
  ! The array must be large enough to hold all elements.
  integer, dimension(:), intent(out) :: Ilist
  
  ! OPTIONAL: Returns the number of elements in the list.
  integer, intent(out), optional :: nelementsInList
!</output>

!<result>
  ! Number of elements in the set.
!</result>

!</subroutine>

    integer :: i,j,k,ncount,iidx

    if (rset%nmaxEntries .eq. rset%ibitsPerInteger*rset%nintegers) then
    
      ! Loop through all integers
      iidx = 0
      ncount = 0
      do i = 1,rset%nintegers
        ! Get the integer
        j = rset%p_Idata(i)
        
        ! Get every element from the bitfield.
        do k=0,rset%ibitsPerInteger-1
          iidx = iidx + 1
          if (btest(j,k)) then
            ncount = ncount+1
            Ilist(ncount) = iidx
          end if
        end do
      end do
      
    else
    
      ! Loop through all integers except the last one.
      ncount = 0
      iidx = 0
      do i = 1,rset%nintegers-1
        ! Get the integer
        j = rset%p_Idata(i)
        
        ! Get every element from the bitfield.
        do k=0,rset%ibitsPerInteger-1
          iidx = iidx + 1
          if (btest(j,k)) then
            ncount = ncount+1
            Ilist(ncount) = iidx
          end if
        end do
      end do
      
      ! Get the last integer
      j = rset%p_Idata(i)
      
      ! Shift it and get the elements.
      do k=0,iand(rset%nmaxEntries,rset%ibitmask)-1
        iidx = iidx + 1
        if (btest(j,k)) then
          ncount = ncount+1
          Ilist(ncount) = iidx
        end if
      end do
    end if
    
    if (present(nelementsInList)) &
      nelementsInList = ncount
    
  end subroutine

  !****************************************************************************

!<subroutine>
  subroutine nsets_putElements (rset,Ilist,nelementsInList)
!<description>
  ! Puts all elements in the list Ilist into the set rset.
!</description>

!<input>
  ! An array with elements for the set.
  integer, dimension(:), intent(in) :: Ilist
  
  ! OPTIONAL: Length of the list Ilist.
  integer, intent(in), optional :: nelementsInList
!</input>

!<inputoutput>
  ! The set to modify.
  type(t_directAccessIntSet), intent(inout) :: rset
!</inputoutput>

!<result>
  ! Number of elements in the set.
!</result>

!</subroutine>

    integer :: i,j,ielement
    
    j = size(Ilist)
    if (present(nelementsInList)) j = nelementsInList
    
    ! Loop through the elements and put them into the set.
    do i = 1,j
      ielement = Ilist(i)
      rset%p_Idata(1+ishft(ielement-1,rset%ishiftDivisor)) = &
          ibset(rset%p_Idata(1+ishft(ielement-1,rset%ishiftDivisor)),&
              iand(ielement-1,rset%ibitmask))
    end do

  end subroutine

end module
