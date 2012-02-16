!-*- mode: f90; -*-

#ifndef _LIST_H_
#define _LIST_H_

!##############################################################################
!# ****************************************************************************
!# <name> template_TD(list,T,D) </name>
!# ****************************************************************************
!#
!# <purpose>
!#
!# This module header file implements a List.
!#
!# Lists are a kind of sequence container. As such, their elements are
!# ordered following a linear sequence.
!#
!# List containers are implemented as doubly linked lists; linked
!# lists can store each of the elements they contain in different and
!# unrelated storage locations. The ordering is kept by the
!# association to each element of a link to the element preceding it
!# and a link to the element following it.
!#
!# This provides the following advantages to list containers:
!# - Efficient insertion and removal of elements anywhere in the
!#   container (constant time).
!# - Efficient moving elements and block of elements within the
!#   container or even between different containers (constant time).
!# - Iterating over the elements in forward order (linear time) or
!#   reverse order (doubly linked lists only)
!#
!# Compared to other base standard sequence containers, lists perform
!# generally better in inserting, extracting and moving elements in
!# any position within the container, and therefore also in algorithms
!# that make intensive use of these, like sorting algorithms.
!#
!# The main drawback of lists compared to other sequence containers is
!# that they lack direct access to the elements by their position; For
!# example, to access the sixth element in a list one has to iterate
!# from a known position (like the beginning or the end) to that
!# position, which takes linear time in the distance between
!# these. They also consume some extra memory to keep the linking
!# information associated to each element (which may be an important
!# factor for large lists of small-sized elements).
!#
!# Storage is handled automatically by the list object, allowing it to
!# be expanded and contracted automatically as needed.
!#
!# -----------------------------------------------------------------------------
!# This module is designed with greatest compatibility to the C++
!# Standard Template Library (STL) concerning naming conventions.
!# -----------------------------------------------------------------------------
!#
!# The following routines are available:
!#
!#  1.) list_create (constructor in STL)
!#      -> Creates an empty list
!#
!#  2.) list_release (destructor in STL)
!#      -> Releases a list
!#
!#  3.) list_resize (resize in STL)
!#      -> Reallocates memory for a list
!#
!#  4.) list_clear (clear in STL)
!#      -> Removes all content from a list
!#
!#  5.) list_copy (no equivalence in STL)
!#      -> Copies data to/from a list
!#
!#  6.) list_swap (swap in STL)
!#      -> Swaps two lists
!#
!#  7.) list_begin / list_rbegin (begin / rbegin in STL)
!#      -> Returns a (reverse) iterator referring to the first/last element 
!#         in a list
!#
!#  8.) list_end / list_rend (end / rend in STL)
!#      -> Returns a (reverse) iterator referring to the element right after
!#         the last/right before the first element in a list
!#
!#  9.) list_next (++iterator in STL)
!#      -> Increases the (reverse) iterator
!#
!# 10.) list_prior (--iterator in STL)
!#      -> Decreases the (reverse) iterator
!#
!# 11.) list_get (no equivalence in STL)
!#      -> Gets key and value at given position in a list
!#
!# 12.) list_assign (assign in STL)
!#      -> Assigns key and value data to a list dropping existing content
!#
!# 13.) list_push_front (push_front in STL)
!#      -> Pushes data to the front of a list
!#
!# 14.) list_push_back (push_back in STL)
!#      -> Pushes data to the back of a list
!#
!# 15.) list_pop_front (pop_front in STL)
!#      -> Removes the first element from a list
!#
!# 16.) list_pop_back (pop_back in STL)
!#      -> Removes the last element from a list
!#
!# 17.) list_front (front in STL)
!#      -> Gets key and value of the first element in a list
!#
!# 18.) list_back (front in STL)
!#      -> Gets key and value of the last element in a list
!#
!# 19.) list_insert (insert in STL)
!#      -> Inserts data into a list
!#
!# 20.) list_erase (erase in STL)
!#      -> Deletes data from a list
!#
!# 21.) list_size (size in STL)
!#      -> Returns the size of a list
!#
!# 22.) list_max_size (max_size in STL)
!#      -> Returns the maximum size of a list
!#
!# 23.) list_empty (empty in STL)
!#      -> Tests if the list is empty
!#
!# 24.) list_find (no equivalence in STL)
!#      -> Searches for data in a list
!#
!# 25.) list_print (no equivalence in STL)
!#      -> Prints content of a list
!#
!# 26.) list_info (no equivalence in STL)
!#      -> Prints information about a list
!#
!# 27.) list_duplicate (no equivalence in STL)
!#      -> Creates a backup of a list
!#
!# 28.) list_restore (no equivalence in STL)
!#      -> Restores a list from a previous backup
!#
!# 29.) list_reverse (reverse in STL)
!#      -> Reverses the order of elements in a list
!#
!# 30.) list_sort (sort in STL)
!#      -> Sorts the entries in a list
!#
!# 31.) list_isNull
!#      -> Tests if the iterator is NULL
!#
!# 32.) list_hasSpec
!#      -> Tests if the iterator has specification flags
!#
!#
!# The following operations in STL are not supported yet
!#
!# splice, remove, remove_if, unique, merge
!#
!#
!# The following operators are available:
!#
!# 1.) "="  assignment operator
!#
!#
!# The following operators in STL are not supported yet
!#
!# lexicographical comparisons; ==, /=, <, <=, >, >=
!#
!#
!# Example:
!#
!# Standard usage of a list is as follows:
!# 
!# ! Initialise an unorderd list for 10 items
!# call list_create(rlist, 10)
!#
!# ! Insert new items at the end of the list
!# call list_push_back(rlist, 1)
!# call list_push_back(rlist, 2)
!# call list_push_back(rlist, 3)
!# call list_push_back(rlist, 4)
!#
!# ! Perform forward iteration through list items
!# rit = list_begin(rlist)
!# do while (rit /= list_end(rlist))
!#   call list_get(rlist, rit, p_key)
!#   print *, "key=",p_key
!#   call list_next(rit)
!# end do
!# 
!# ! Perform reverse iteration through list items
!# rit = list_rbegin(rlist)
!# do while (rit /= list_rend(rlist))
!#   call list_get(rlist, rit, p_key)
!#   print *, "key=",p_key
!#   call list_next(rit)
!# end do
!#
!# ! Release list
!# call list_release(rlist)
!#
!# </purpose>
!##############################################################################

#include "template.h"

  implicit none

  private
  public :: template_TD(t_list,T,D)
  public :: template_TD(it_list,T,D)
  public :: list_create
  public :: list_release
  public :: list_resize
  public :: list_clear
  public :: list_copy
  public :: list_swap
  public :: list_begin
  public :: list_rbegin
  public :: list_end
  public :: list_rend
  public :: list_next
  public :: list_prior
  public :: list_get
  public :: list_assign
  public :: list_push_front
  public :: list_push_back
  public :: list_pop_front
  public :: list_pop_back
  public :: list_front
  public :: list_back
  public :: list_insert
  public :: list_erase
  public :: list_size
  public :: list_max_size
  public :: list_empty
  public :: list_find
  public :: list_print
  public :: list_info
  public :: list_duplicate
  public :: list_restore
  public :: list_reverse
  public :: list_sort
  public :: list_isNull
  public :: list_hasSpec
  
  public assignment(=)
  public operator(==)
  public operator(/=)
  public operator(<)
  public operator(<=)
  public operator(>)
  public operator(>=)

  interface list_create
    module procedure template_TD(list_create,T,D)
  end interface

  interface list_release
    module procedure template_TD(list_release,T,D)
  end interface

  interface list_resize
    module procedure template_TD(list_resize,T,D)
  end interface

  interface list_clear
    module procedure template_TD(list_clear,T,D)
  end interface

  interface list_copy
    module procedure template_TD(list_cpy1,T,D)
    module procedure template_TD(list_cpy2,T,D)
    module procedure template_TD(list_cpy3,T,D)
    module procedure template_TD(list_cpy4,T,D)
  end interface

  interface list_swap
    module procedure template_TD(list_swap,T,D)
  end interface

  interface list_begin
    module procedure template_TD(list_begin,T,D)
  end interface

  interface list_rbegin
    module procedure template_TD(list_rbegin,T,D)
  end interface
  
  interface list_end
    module procedure template_TD(list_end,T,D)
  end interface

  interface list_rend
    module procedure template_TD(list_rend,T,D)
  end interface

  interface list_next
    module procedure template_TD(list_next,T,D)
  end interface

  interface list_prior
    module procedure template_TD(list_prior,T,D)
  end interface

  interface list_get
    module procedure template_TD(list_get1,T,D)
    module procedure template_TD(list_get2,T,D)
  end interface

  interface list_assign
    module procedure template_TD(list_assign1,T,D)
    module procedure template_TD(list_assign2,T,D)
  end interface

  interface list_push_front
    module procedure template_TD(list_push_front,T,D)
  end interface

  interface list_push_back
    module procedure template_TD(list_push_back,T,D)
  end interface

  interface list_pop_front
    module procedure template_TD(list_pop_front,T,D)
  end interface
  
  interface list_pop_back
    module procedure template_TD(list_pop_back,T,D)
  end interface

  interface list_front
    module procedure template_TD(list_front1,T,D)
    module procedure template_TD(list_front2,T,D)
  end interface
  
  interface list_back
    module procedure template_TD(list_back1,T,D)
    module procedure template_TD(list_back2,T,D)
  end interface

  interface list_insert
    module procedure template_TD(list_insert1,T,D)
    module procedure template_TD(list_insert2,T,D)
    module procedure template_TD(list_insert3,T,D)
  end interface

  interface list_erase
    module procedure template_TD(list_erase1,T,D)
    module procedure template_TD(list_erase2,T,D)
  end interface

  interface list_size
    module procedure template_TD(list_size,T,D)
  end interface
  
  interface list_max_size
    module procedure template_TD(list_max_size,T,D)
  end interface

  interface list_empty
    module procedure template_TD(list_empty,T,D)
  end interface
  
  interface list_find
    module procedure template_TD(list_find1,T,D)
    module procedure template_TD(list_find2,T,D)
  end interface

  interface list_print
    module procedure template_TD(list_print,T,D)
  end interface

  interface list_info
    module procedure template_TD(list_info,T,D)
  end interface
  
  interface list_duplicate
    module procedure template_TD(list_duplicate,T,D)
  end interface
  
  interface list_restore
    module procedure template_TD(list_restore,T,D)
  end interface
  
  interface list_reverse
    module procedure template_TD(list_reverse,T,D)
  end interface

  interface list_sort
    module procedure template_TD(list_sort,T,D)
  end interface

  interface list_isNull
    module procedure template_TD(list_isNull,T,D)
  end interface

  interface list_hasSpec
    module procedure template_TD(list_hasSpec,T,D)
  end interface


  interface assignment(=)
    module procedure template_TD(list_fassign,T,D)
  end interface

  
  interface operator(==)
    module procedure template_TD(it_list_eq,T,D)
  end interface
  
  interface operator(/=)
    module procedure template_TD(it_list_ne,T,D)
  end interface

  interface operator(<)
    module procedure template_TD(it_list_lt,T,D)
  end interface

  interface operator(<=)
    module procedure template_TD(it_list_le,T,D)
  end interface

  interface operator(>)
    module procedure template_TD(it_list_gt,T,D)
  end interface

  interface operator(>=)
    module procedure template_TD(it_list_ge,T,D)
  end interface

  !************************************************************************

!<types>

!<typeblock>
  
  ! A doubly linked list

  type template_TD(t_list,T,D)
    private

    ! Number of elements that are currently stored in the list
    integer :: NA = 0

    ! Total number of elements that can be stored in the list
    integer :: NNA = 0

    ! Total number of elements that can initially be stored in all
    ! lists. This information is needed to compute the growth of the
    ! arraylist after several resize operations
    integer :: NNA0 = 0

    ! Total number of resize operations
    integer :: nresize = 0

    ! Factor by which the list is enlarged if new storage is allocate
    real(DP) :: dfactor = 1.5_DP

#ifdef T_STORAGE
    ! Handle to the list key
    integer :: h_Key = ST_NOHANDLE
#endif

    ! Pointer to list key
    TTYPE(T_TYPE), dimension(:), pointer :: p_Key => null()

    ! Handle and pointer to the list next-structure
    integer :: h_Knext = ST_NOHANDLE
    integer, dimension(:), pointer :: p_Knext => null()

    ! Handle and pointer to the list previous-structure
    integer :: h_Kprev = ST_NOHANDLE
    integer, dimension(:), pointer :: p_Kprev => null()

#ifdef D
    ! Dimension of the auxiliary data values to be stored
    integer :: isizeData = 0
#ifdef D_STORAGE
    ! Handle to the list auxiliary data
    integer :: h_Data = ST_NOHANDLE
#endif

    ! Pointer to list auxiliary data
    DTYPE(D_TYPE), dimension(:,:), pointer :: p_Data => null()
#endif
    
  end type
  
!</typeblock>

!<typeblock>

  ! An iterator for a doubly linked list

  type template_TD(it_list,T,D)
    private

    ! Absolute position of the current element
    integer :: ipos = LNULL

    ! Specification flag. This is a bitfield coming from an OR
    ! combination of different LIST_LSPEC_xxxx constants and specifies
    ! various details of the list iterator.
    integer(I32) :: iSpec = 0_I32
    
    ! Pointer to the underlying doubly linked list
    type(template_TD(t_list,T,D)), pointer :: p_rlist => null()

  end type

!</typeblock>

!</types>

contains

  !************************************************************************
  
!<subroutine>

#ifdef D
  subroutine template_TD(list_create,T,D)(rlist, NNA, isizeData, dfactor)
#else
  subroutine template_TD(list_create,T,D)(rlist, NNA, dfactor)
#endif

!<description>
    ! This subroutine creates a new list
!</description>

!<input>
    ! Total number of elements that can be stored in list
    integer, intent(in) :: NNA

#ifdef D
    ! Dimension of the auxiliary data values to be stored
    integer, intent(in) :: isizeData
#endif

    ! OPTIONAL: Factor by which the list should be enlarged if memory
    ! needs to be reallocated
    real(DP), intent(in), optional :: dfactor
!</input>

!<output>
    ! The linked list
    type(template_TD(t_list,T,D)), intent(out) :: rlist
!</output>
!</subroutine>

#ifdef D
    ! local variables
    integer, dimension(2) :: Isize
#endif

    ! Set factor
    if (present(dfactor)) then
      if (dfactor > 1.0_DP) rlist%dfactor = dfactor
    end if
    
    ! Initialise list
    rlist%NNA  = max(0,NNA)
    rlist%NNA0 = max(0,NNA)
    
    ! Allocate memory and associate pointers
    call storage_new('list_create', 'Knext', LHEAD, rlist%NNA, ST_INT,&
        rlist%h_Knext, ST_NEWBLOCK_NOINIT)
    call storage_getbase_int(rlist%h_Knext, rlist%p_Knext)
    
    call storage_new('list_create','Kprev', rlist%NNA, ST_INT,&
        rlist%h_Kprev, ST_NEWBLOCK_NOINIT)
    call storage_getbase_int(rlist%h_Kprev, rlist%p_Kprev)

    ! Allocate memory for Key
#ifdef T_STORAGE
    call storage_new('list_create', 'Key', rlist%NNA, T_STORAGE,&
        rlist%h_Key, ST_NEWBLOCK_NOINIT)
    call storage_getbase(rlist%h_Key, rlist%p_Key)
#else
    allocate(rlist%p_Key(rlist%NNA))
#endif
    
#ifdef D
    ! Set size of auxiliary data
    rlist%isizeData = max(0,isizeData)

    ! Allocate memory for auxiliary data
    if (rlist%isizeData > 0) then
      Isize = (/rlist%isizeData, rlist%NNA/)

#ifdef D_STORAGE
      call storage_new('list_create', 'Data', Isize, D_STORAGE,&
          rlist%h_Data, ST_NEWBLOCK_NOINIT)
      call storage_getbase(rlist%h_Data, rlist%p_Data)
#else
      allocate(rlist%p_Data(Isize(1),Isize(2)))
#endif
    end if
#endif

    ! Clear list
    call list_clear(rlist)

  end subroutine

  !************************************************************************

!<subroutine>
  
  subroutine template_TD(list_release,T,D)(rlist)

!<description>
    ! This subroutine releases an existing list
!</description>

!<inputoutput>
    ! The linked list
    type(template_TD(t_list,T,D)), intent(inout) :: rlist
!</inputoutput>
!</subroutine>

    ! Release memory
    if (rlist%h_Knext .ne. ST_NOHANDLE) call storage_free(rlist%h_Knext)
    if (rlist%h_Kprev .ne. ST_NOHANDLE) call storage_free(rlist%h_Kprev)

    nullify(rlist%p_Knext)
    nullify(rlist%p_Kprev)

#ifdef T_STORAGE
    if (rlist%h_Key   .ne. ST_NOHANDLE) call storage_free(rlist%h_Key)
    nullify(rlist%p_Key)
#else
    if (associated(rlist%p_Key))&
        deallocate(rlist%p_Key)
#endif

#ifdef D
    ! Reset size of auxiliary data
    rlist%isizeData = 0
#ifdef D_STORAGE
    if (rlist%h_Data  .ne. ST_NOHANDLE) call storage_free(rlist%h_Data)
    nullify(rlist%p_Data)
#else
    if (associated(rlist%p_Data))&
        deallocate(rlist%p_Data)
#endif
#endif
    
    ! Reset linked list
    rlist%NA        = 0
    rlist%NNA       = 0
    rlist%dfactor   = 1.5_DP

  end subroutine

  !************************************************************************

!<subroutine>
  
  subroutine template_TD(list_resize,T,D)(rlist, NNA)

!<description>
    ! This subroutine reallocates memory for an existing list
!</description>

!<input>
    ! New number of total elements that can be stored in the list
    integer, intent(in) :: NNA
!</input>

!<inputoutput>
    ! The linked list
    type(template_TD(t_list,T,D)), intent(inout) :: rlist
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: NNAOld
    
#ifndef T_STORAGE
    TTYPE(T_TYPE), dimension(:), pointer :: p_Key
#endif

#ifdef D
#ifndef D_STORAGE
    DTYPE(D_TYPE), dimension(:,:), pointer :: p_Data
#endif
#endif

    ! Save old size
    NNAOld = rlist%NNA

    ! Set new size
    rlist%NNA = max(0,NNA)
    rlist%nresize = rlist%nresize+1

    ! Reallocate structures
    call storage_realloc('list_resize', rlist%NNA+3, rlist%h_Knext,&
        ST_NEWBLOCK_NOINIT, .true.)
    call storage_getbase_int(rlist%h_Knext, rlist%p_Knext)

    call storage_realloc('list_resize', rlist%NNA, rlist%h_Kprev,&
        ST_NEWBLOCK_NOINIT, .true.)
    call storage_getbase_int(rlist%h_Kprev, rlist%p_Kprev)

    ! Reallocate Key
#ifdef T_STORAGE
    call storage_realloc('list_resize', rlist%NNA, rlist%h_Key,&
        ST_NEWBLOCK_NOINIT, .true.)
    call storage_getbase(rlist%h_Key, rlist%p_Key)
#else
    allocate(p_Key(NNAOld))
    p_Key = rlist%p_Key
    deallocate(rlist%p_Key)
    allocate(rlist%p_Key(rlist%NNA))
    rlist%p_Key(1:NNAOld) = p_Key
    deallocate(p_Key)
#endif
    
#ifdef D
    ! Reallocate auxiliary data
    if (rlist%isizeData > 0) then
#ifdef D_STORAGE
      call storage_realloc('list_resize', rlist%NNA, rlist%h_Data,&
          ST_NEWBLOCK_NOINIT, .true.)
      call storage_getbase(rlist%h_Data, rlist%p_Data)
#else
      allocate(p_Data(size(rlist%p_Data,1),NNAOld))
      p_Data = rlist%p_Data
      deallocate(rlist%p_Data)
      allocate(rlist%p_Data(size(p_Data,1),rlist%NNA))
      rlist%p_Data(:,1:NNAOld) = p_Data
      deallocate(p_Data)
#endif
    end if
#endif

  end subroutine

  !************************************************************************

!<subroutine>
  
  pure subroutine template_TD(list_clear,T,D)(rlist)

!<description>
    ! This subroutine clears the content of the list
!</description>

!<inputoutput>
    ! The linked list
    type(template_TD(t_list,T,D)), intent(inout) :: rlist
!</inputoutput>
!</subroutine>

    ! Initialise list structure
    rlist%NA      = 0
    rlist%nresize = 0
    rlist%p_Knext(LFREE) = 1
    rlist%p_Knext(LHEAD) = LNULL
    rlist%p_Knext(LTAIL) = LNULL

  end subroutine

  !************************************************************************

!<subroutine>

#ifdef D
  subroutine template_TD(list_cpy1,T,D)(h_KeySrc, rlist, h_DataSrc)
#else
  subroutine template_TD(list_cpy1,T,D)(h_KeySrc, rlist)
#endif

!<description>
    ! This subroutine copies the content of the given storage
    ! handle(s) to the list
!</description>

!<input>
    ! Handle to the key values
    integer, intent(in) :: h_KeySrc

#ifdef D
    ! OPTIONAL: Handle to the data values
    integer, intent(in), optional :: h_DataSrc
#endif
!</input>

!<inputoutput>
    ! The linked list
    type(template_TD(t_list,T,D)), intent(inout) :: rlist
!</inputoutput>
!</subroutine>

#ifdef T_STORAGE
    ! local variables
    T_TYPE, dimension(:), pointer :: p_KeySrc

#ifdef D
#ifdef D_STORAGE
    D_TYPE, dimension(:,:), pointer :: p_DataSrc

    ! Get pointers and copy the key and data values to the list
    call storage_getbase(h_KeySrc, p_KeySrc)
    call storage_getbase(h_DataSrc, p_DataSrc)
    call list_copy(p_KeySrc, rlist, p_DataSrc)
#else
    call output_line('List does not support storage handles!',&
        OU_CLASS_ERROR,OU_MODE_STD,'list_cpy1')
    call sys_halt()
#endif

#else
    ! Get pointer and copy the key values to the list
    call storage_getbase(h_KeySrc, p_KeySrc)
    call list_copy(p_KeySrc, rlist)
#endif

#else
    call output_line('List does not support storage handles!',&
        OU_CLASS_ERROR,OU_MODE_STD,'list_cpy1')
    call sys_halt()
#endif

  end subroutine

  !************************************************************************

!<subroutine>

#ifdef D
  subroutine template_TD(list_cpy2,T,D)(KeySrc, rlist, DataSrc)
#else
  subroutine template_TD(list_cpy2,T,D)(KeySrc, rlist)
#endif

!<description>
    ! This subroutine copies the content of the given array(s) to the list
!</description>

!<input>
    ! Array with key values
    TTYPE(T_TYPE), dimension(:), intent(in) :: KeySrc

#ifdef D
    ! OPTIONAL: Array with data values
    DTYPE(D_TYPE), dimension(:,:), intent(in), optional :: DataSrc
#endif
!</input>

!<inputoutput>
    ! The linked list
    type(template_TD(t_list,T,D)), intent(inout) :: rlist
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: ipos

#ifdef D
    if (present(DataSrc)) then
      do ipos = 1, size(KeySrc)
        call list_push_back(rlist, KeySrc(ipos), dataSrc(:,ipos))      
      end do
    else
      do ipos = 1, size(KeySrc)
        call list_push_back(rlist, KeySrc(ipos))
      end do
    end if
#else
    do ipos = 1, size(KeySrc)
      call list_push_back(rlist, KeySrc(ipos))
    end do
#endif

  end subroutine

  !************************************************************************

!<subroutine>

#ifdef D
  subroutine template_TD(list_cpy3,T,D)(rlist, h_KeyDest, h_DataDest)
#else
  subroutine template_TD(list_cpy3,T,D)(rlist, h_KeyDest)
#endif

!<description>
    ! This subroutine copies the content of the given list to the handle(s)
!</description>

!<input>
    ! The linked list
    type(template_TD(t_list,T,D)), intent(in) :: rlist
!</input>

!<inputoutput>
    ! Handle to the key values
    integer, intent(inout) :: h_KeyDest

#ifdef D
    ! OPTIONAL: Handle to the data values
    integer, intent(inout), optional :: h_DataDest
#endif
!</inputoutput>
!</subroutine>

    ! local variable
    integer :: isize

#ifdef T_STORAGE
    T_TYPE, dimension(:), pointer :: p_KeyDest

#ifdef D
    integer, dimension(2) :: Isize2
#ifdef D_STORAGE
    D_TYPE, dimension(:,:), pointer :: p_DataDest

    if (present(h_DataDest)) then
    
      ! Check key handle
      if (h_KeyDest .eq. ST_NOHANDLE) then
        call storage_new('list_cpy3', 'Key',&
            rlist%NA, T_STORAGE, h_KeyDest, ST_NEWBLOCK_NOINIT)
      else
        call storage_getsize(h_KeyDest, isize)
        if (isize < rlist%NA) then
          call storage_realloc('list_cpy3',&
              rlist%NA, h_KeyDest, ST_NEWBLOCK_NOINIT, .false.)
        end if
      end if
      call storage_getbase(h_KeyDest, p_KeyDest)
      
      if (rlist%isizeData <= 0) then
        if (h_DataDest .ne. ST_NOHANDLE) call storage_free(h_DataDest)
        call output_line('List does not provide auxiliary data!',&
            OU_CLASS_WARNING,OU_MODE_STD,'list_cpy3')

        ! Copy key values from list
        call list_copy(rlist, p_KeyDest)
        
        ! That`s it
        return
      end if

      ! Check data handle
      if (h_DataDest .eq. ST_NOHANDLE) then
        call storage_new('list_cpy3', 'Data',&
            (/rlist%isizeData,rlist%NA/),&
            D_STORAGE, h_DataDest, ST_NEWBLOCK_NOINIT)
      else
        call storage_getsize(h_DataDest, Isize2)

        if (Isize2(1) .ne. rlist%isizeData) then
          call output_line('Size of data array is not compatible!',&
              OU_CLASS_ERROR,OU_MODE_STD,'list_cpy3')
        end if
        
        if (Isize2(2) < rlist%NA) then
          call storage_realloc('list_cpy3',&
              rlist%NA, h_DataDest, ST_NEWBLOCK_NOINIT, .false.)
        end if
      end if
      call storage_getbase(h_DataDest, p_DataDest)
              
      ! Copy key and data values from list
      call list_copy(rlist, p_KeyDest, p_DataDest)
      
    else

      ! Check key handle
      if (h_KeyDest .eq. ST_NOHANDLE) then
        call storage_new('list_cpy3', 'Key',&
            rlist%NA, T_STORAGE, h_KeyDest, ST_NEWBLOCK_NOINIT)
      else
        call storage_getsize(h_KeyDest, isize)
        if (isize < rlist%NA) then
          call storage_realloc('list_cpy3',&
              rlist%NA, h_KeyDest, ST_NEWBLOCK_NOINIT, .false.)
        end if
      end if
      call storage_getbase(h_KeyDest, p_KeyDest)
      
      ! Copy key values from list
      call list_copy(rlist, p_KeyDest)
      
    end if
#else
    call output_line('List does not support storage handles!',&
        OU_CLASS_ERROR,OU_MODE_STD,'list_cpy3')
    call sys_halt()
#endif

#else

    ! Check key handle
    if (h_KeyDest .eq. ST_NOHANDLE) then
      call storage_new('list_cpy3', 'Key',&
          rlist%NA, T_STORAGE, h_KeyDest, ST_NEWBLOCK_NOINIT)
    else
      call storage_getsize(h_KeyDest, isize)
      if (isize < rlist%NA) then
        call storage_realloc('list_cpy3',&
            rlist%NA, h_KeyDest, ST_NEWBLOCK_NOINIT, .false.)
      end if
    end if
    call storage_getbase(h_KeyDest, p_KeyDest)
    
    ! Copy key values from list
    call list_copy(rlist, p_KeyDest)

#endif

#else
    call output_line('List does not support storage handles!',&
        OU_CLASS_ERROR,OU_MODE_STD,'list_cpy3')
    call sys_halt()
#endif

  end subroutine
  
  !************************************************************************

!<subroutine>

#ifdef D
  subroutine template_TD(list_cpy4,T,D)(rlist, KeyDest, DataDest)
#else
  subroutine template_TD(list_cpy4,T,D)(rlist, KeyDest)
#endif

!<description>
    ! This subroutine copies the content of the given list to the array(s)
!</description>

!<input>
    ! The linked list
    type(template_TD(t_list,T,D)), intent(in) :: rlist
!</input>

!<inputoutput>
    ! Array with key values
    TTYPE(T_TYPE), dimension(:), intent(inout) :: KeyDest

#ifdef D
    ! OPTIONAL: Array with data values
    DTYPE(D_TYPE), dimension(:,:), intent(inout), optional :: DataDest
#endif
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: ipos,jpos

    ipos = rlist%p_Knext(LHEAD)
    jpos = 0

#ifdef D

    if (present(DataDest)) then

      kd_copy: do
        jpos = jpos+1
        KeyDest(jpos)    = rlist%p_Key(ipos)
        DataDest(:,jpos) = rlist%p_Data(:,ipos)
        if (ipos .eq. rlist%p_Knext(LTAIL)) exit kd_copy
        ipos = rlist%p_Knext(ipos)
      end do kd_copy

    else

      k_copy: do
        jpos = jpos+1
        KeyDest(jpos) = rlist%p_Key(ipos)
        if (ipos .eq. rlist%p_Knext(LTAIL)) exit k_copy
        ipos = rlist%p_Knext(ipos)
      end do k_copy

    end if

#else

    copy: do
      jpos = jpos+1
      KeyDest(jpos) = rlist%p_Key(ipos)
      if (ipos .eq. rlist%p_Knext(LTAIL)) exit copy
      ipos = rlist%p_Knext(ipos)
    end do copy

#endif

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine template_TD(list_swap,T,D)(rlist1, rlist2)

!<description>
    ! This subroutine swaps content of two lists
!</description>
    
!<inputoutput>
    ! First linked list
    type(template_TD(t_list,T,D)), intent(inout) :: rlist1

    ! Second linked list
    type(template_TD(t_list,T,D)), intent(inout) :: rlist2
!</inputoutput>
!</subroutine>

    ! local variables
    type(template_TD(t_list,T,D)) :: rlist
    
    ! Swap lists
    rlist  = rlist1
    rlist1 = rlist2
    rlist2 = rlist

    ! Reassociated pointers of first list
    call storage_getbase_int(rlist1%h_Knext, rlist1%p_Knext)
    call storage_getbase_int(rlist1%h_Kprev, rlist1%p_Kprev)

    ! Reassociated pointers of second list
    call storage_getbase_int(rlist2%h_Knext, rlist2%p_Knext)
    call storage_getbase_int(rlist2%h_Kprev, rlist2%p_Kprev)

#ifdef T_STORAGE
    ! Reassociated pointers to key values
    call storage_getbase(rlist1%h_Key, rlist1%p_Key)
    call storage_getbase(rlist2%h_Key, rlist2%p_Key)
#endif

#if defined D && defined D_STORAGE
    ! Reassociated pointers to data values
    if (rlist1%h_Data .ne. ST_NOHANDLE) &
        call storage_getbase(rlist1%h_Data, rlist1%p_Data)
    if (rlist2%h_Data .ne. ST_NOHANDLE) &
        call storage_getbase(rlist2%h_Data, rlist2%p_Data)
#endif
    
  end subroutine

  !************************************************************************

!<function>
  
  function template_TD(list_begin,T,D)(rlist) result(riterator)

!<description>
    ! This function returns an iterator referring to the first
    ! element of the list
!</description>

!<input>
    ! The linked list
    type(template_TD(t_list,T,D)), intent(in), target :: rlist
!</input>

!<result>
    ! The iterator
    type(template_TD(it_list,T,D)) :: riterator
!</result>

!</function>
    
    ! Attach list to iterator
    riterator%p_rlist => rlist

    ! Initialise iterator
    riterator%ipos  = rlist%p_Knext(LHEAD)
    riterator%iSpec = 0_I32

  end function

  !************************************************************************

!<function>
  
  function template_TD(list_rbegin,T,D)(rlist) result(riterator)

!<description>
    ! This function returns a reverse iterator positioned to the last
    ! element of the list
!</description>

!<input>
    ! The linked list
    type(template_TD(t_list,T,D)), intent(in), target :: rlist
!</input>

!<result>
    ! The iterator
    type(template_TD(it_list,T,D)) :: riterator
!</result>

!</function>
    
    ! Attach list to iterator
    riterator%p_rlist => rlist

    ! Initialise iterator
    riterator%ipos  = rlist%p_Knext(LTAIL)
    riterator%iSpec = LIST_LSPEC_REVERSE

  end function

  !************************************************************************

!<function>
  
  function template_TD(list_end,T,D)(rlist) result(riterator)

!<description>
    ! This function returns an iterator referring to the past-the-end
    ! element of the list
!</description>

!<input>
    ! The linked list
    type(template_TD(t_list,T,D)), intent(in), target :: rlist
!</input>

!<result>
    ! The iterator
    type(template_TD(it_list,T,D)) :: riterator
!</result>
!</function>
    
    ! Attach list to iterator
    riterator%p_rlist => rlist

    ! Initialise iterator
    riterator%ipos  = LNULL
    riterator%iSpec = 0_I32

  end function

  !************************************************************************

!<function>
  
  function template_TD(list_rend,T,D)(rlist) result(riterator)

!<description>
    ! This function returns a reverse iterator referring to the
    ! element right before the first element of the list
!</description>

!<input>
    ! The linked list
    type(template_TD(t_list,T,D)), intent(in), target :: rlist
!</input>

!<result>
    ! The iterator
    type(template_TD(it_list,T,D)) :: riterator
!</result>
!</function>
    
    ! Attach list to iterator
    riterator%p_rlist => rlist

    ! Initialise iterator
    riterator%ipos  = LNULL
    riterator%iSpec = LIST_LSPEC_REVERSE

  end function

  !************************************************************************

!<subroutine>
  
  pure subroutine template_TD(list_next,T,D)(riterator)

!<description>
    ! This subroutine increments the list iterator by one.
!</description>

!<inputoutput>
    ! The iterator
    type(template_TD(it_list,T,D)), intent(inout) :: riterator
!</inputoutput>
!</subroutine>
    
    if (.not.list_isNull(riterator)) then
      if (iand(riterator%iSpec, LIST_LSPEC_REVERSE).eq.0) then
        riterator%ipos = riterator%p_rlist%p_Knext(riterator%ipos)
      else
        riterator%ipos = riterator%p_rlist%p_Kprev(riterator%ipos)
      end if
    else
      if (iand(riterator%iSpec, LIST_LSPEC_REVERSE).eq.0) then
        riterator%ipos = riterator%p_rlist%p_Knext(LHEAD)
      else
        riterator%ipos = riterator%p_rlist%p_Knext(LTAIL)
      end if
    end if
    
  end subroutine

  !************************************************************************

!<subroutine>
  
  pure subroutine template_TD(list_prior,T,D)(riterator)

!<description>
    ! This subroutine decrements the list iterator by one.
!</description>

!<inputoutput>
    ! The iterator
    type(template_TD(it_list,T,D)), intent(inout) :: riterator
!</inputoutput>
!</subroutine>
    
    if (.not.list_isNull(riterator)) then
      if (iand(riterator%iSpec, LIST_LSPEC_REVERSE) .eq. 0) then
        riterator%ipos = riterator%p_rlist%p_Kprev(riterator%ipos)
      else
        riterator%ipos = riterator%p_rlist%p_Knext(riterator%ipos)
      end if
    else
      if (iand(riterator%iSpec, LIST_LSPEC_REVERSE) .eq. 0) then
        riterator%ipos = riterator%p_rlist%p_Knext(LTAIL)
      else
        riterator%ipos = riterator%p_rlist%p_Knext(LHEAD)
      end if
    end if
    
  end subroutine

  !************************************************************************

!<subroutine>

#ifdef D
  subroutine template_TD(list_get1,T,D)(rlist, rposition, p_key, p_data)
#else
  subroutine template_TD(list_get1,T,D)(rlist, rposition, p_key)
#endif

!<description>
    ! This subroutine returns pointers to the key and data stored at
    ! the position addressed by the iterator.
!</description>

!<input>
    ! The linked list
    type(template_TD(t_list,T,D)), intent(in) :: rlist

    ! The iterator
    type(template_TD(it_list,T,D)), intent(in) :: rposition
!</input>

!<output>
    ! Pointer to the key value
    TTYPE(T_TYPE), intent(out), pointer :: p_key

#ifdef D
    ! OPTIONAL: Pointer to the data
    DTYPE(D_TYPE), dimension(:), intent(out), pointer, optional :: p_data
#endif
!</output>
!</subroutine>
    
    ! Get key
    p_key => rlist%p_key(rposition%ipos)

#ifdef D
    ! Get data
    if (present(p_data)) then
      if (rlist%isizeData > 0) then
        p_data => rlist%p_Data(:,rposition%ipos)
      else
        nullify(p_data)
      end if
    end if
#endif

  end subroutine

  !************************************************************************

!<function>

  function template_TD(list_get2,T,D)(rlist, rposition) result(key)

!<description>
    ! This functions return the key stored at the position addressed
    ! by the iterator.
!</description>

!<input>
    ! The linked list
    type(template_TD(t_list,T,D)), intent(in) :: rlist

    ! The iterator
    type(template_TD(it_list,T,D)), intent(in) :: rposition
!</input>

!<result>
    ! key value
    TTYPE(T_TYPE) :: key
!</result>
!</function>
    
    ! Get key
    key = rlist%p_key(rposition%ipos)

  end function

  !************************************************************************

!<subroutine>

#ifdef D
  subroutine template_TD(list_assign1,T,D)(rlist, n, key, data)
#else
  subroutine template_TD(list_assign1,T,D)(rlist, n, key)
#endif

!<description>
    ! This subroutine removes all existing data from the list and
    ! assigns n copies of the given key and data values
!</description>

!<input>
    ! Number of copies
    integer, intent(in) :: n

    ! key value
    TTYPE(T_TYPE), intent(in) :: key

#ifdef D
    ! OPTIONAL: Data
    DTYPE(D_TYPE), dimension(:), intent(in), optional :: data
#endif
!</input>

!<inputoutput>
    ! The linked list
    type(template_TD(t_list,T,D)), intent(inout) :: rlist
!</inputoutput>
!</subroutine>

    ! local variable
    integer :: i

    ! Clear list
    call list_clear(rlist)

    ! Push copies to the list
    do i = 1, n
#ifdef D
      call list_push_back(rlist, key, data)
#else
      call list_push_back(rlist, key)
#endif
    end do

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine template_TD(list_assign2,T,D)(rlist, rfirst, rlast)

!<description>
    ! This subroutine removes all existing data from the list and
    ! assigns the content in the range [rfirst,rlast) to the list.
!</description>

!<input>
    ! Iterator referring to the first element
    type(template_TD(it_list,T,D)), intent(in) :: rfirst

    ! Iterator referring to the past-the-end element
    type(template_TD(it_list,T,D)), intent(in) :: rlast
!</input>

!<inputoutput>
    ! The linked list
    type(template_TD(t_list,T,D)), intent(inout) :: rlist
!</inputoutput>
!</subroutine>

    ! local variable
    type(template_TD(it_list,T,D)) :: riterator
    TTYPE(T_TYPE), pointer :: p_key

#ifdef D
    DTYPE(D_TYPE), dimension(:), pointer :: p_data
#endif

    ! Clear list
    call list_clear(rlist)

    ! Push content to the list
    riterator = rfirst
    do while(riterator /= rlast)

#ifdef D
      call list_get(riterator%p_rlist, riterator, p_key, p_data)
      call list_push_back(rlist, p_key, p_data)
#else
      call list_get(riterator%p_rlist, riterator, p_key)
      call list_push_back(rlist, p_key)
#endif
      call list_next(riterator)

    end do
    
  end subroutine

  !************************************************************************

!<subroutine>

#ifdef D
  subroutine template_TD(list_push_front,T,D)(rlist, key, data)
#else
  subroutine template_TD(list_push_front,T,D)(rlist, key)
#endif

!<description>
    ! This subroutine prepends a new element to the list
!</description>

!<input>
    ! key value
    TTYPE(T_TYPE), intent(in) :: key

#ifdef D
    ! OPTIONAL: Data
    DTYPE(D_TYPE), dimension(:), intent(in), optional :: data
#endif
!</input>

!<inputoutput>
    ! The linked list
    type(template_TD(t_list,T,D)), intent(inout) :: rlist
!</inputoutput>
!</subroutine>

    ! local variable
    integer :: ipos

    ! Check if list needs to be resized
    rlist%NA = rlist%NA+1
    ipos     = rlist%p_Knext(LFREE)
    if (abs(ipos) > rlist%NNA) then
      call list_resize(rlist, ceiling(rlist%dfactor*rlist%NNA))
    end if

    ! Set next free position
    if (ipos > 0) then
      rlist%p_Knext(LFREE) = ipos+1
    else
      ipos = abs(ipos)
      rlist%p_Knext(LFREE) = rlist%p_Knext(ipos)
    end if

    ! Set head and tail
    if (rlist%p_Knext(LHEAD) .eq. LNULL) then
      ! Push to empty list
      rlist%p_Knext(LHEAD) = ipos
      rlist%p_Knext(LTAIL) = ipos
      rlist%p_Knext(ipos)  = LNULL
      rlist%p_Kprev(ipos)  = LNULL
    else
      ! Push to non-empty list
      rlist%p_Kprev(rlist%p_Knext(LHEAD)) = ipos
      rlist%p_Knext(ipos)  = rlist%p_Knext(LHEAD)
      rlist%p_Knext(LHEAD) = ipos
      rlist%p_Kprev(ipos)  = LNULL
    end if

    ! Store key value
    rlist%p_Key(ipos) = key

#ifdef D
    ! Store data value
    if (present(data)) then
      rlist%p_Data(:,ipos) = data
    end if
#endif
    
  end subroutine

  !************************************************************************

!<subroutine>

#ifdef D
  subroutine template_TD(list_push_back,T,D)(rlist, key, data)
#else
  subroutine template_TD(list_push_back,T,D)(rlist, key)
#endif

!<description>
    ! This subroutine appends a new element to the list
!</description>

!<input>
    ! key value
    TTYPE(T_TYPE), intent(in) :: key

#ifdef D
    ! OPTIONAL: Data
    DTYPE(D_TYPE), dimension(:), intent(in), optional :: data
#endif
!</input>

!<inputoutput>
    ! The linked list
    type(template_TD(t_list,T,D)), intent(inout) :: rlist
!</inputoutput>
!</subroutine>

    ! local variable
    integer :: ipos

    ! Check if list needs to be resized
    rlist%NA = rlist%NA+1
    ipos     = rlist%p_Knext(LFREE)
    if (abs(ipos) > rlist%NNA) then
      call list_resize(rlist, ceiling(rlist%dfactor*rlist%NNA))
    end if

    ! Set next free position
    if (ipos > 0) then
      rlist%p_Knext(LFREE) = ipos+1
    else
      ipos = abs(ipos)
      rlist%p_Knext(LFREE) = rlist%p_Knext(ipos)
    end if
    
    ! Set head and tail
    if (rlist%p_Knext(LHEAD) .eq. LNULL) then
      ! Push to empty list
      rlist%p_Knext(LHEAD) = ipos
      rlist%p_Knext(LTAIL) = ipos
      rlist%p_Knext(ipos)  = LNULL
      rlist%p_Kprev(ipos)  = LNULL
    else
      ! Push to non-empty list
      rlist%p_Kprev(ipos)  = rlist%p_Knext(LTAIL)
      rlist%p_Knext(rlist%p_Knext(LTAIL)) = ipos
      rlist%p_Knext(LTAIL) = ipos
      rlist%p_Knext(ipos)  = LNULL
    end if

    ! Store key value
    rlist%p_Key(ipos) = key

#ifdef D
    ! Store data value
    if (present(data)) then
      rlist%p_Data(:,ipos) = data
    end if
#endif

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine template_TD(list_pop_front,T,D)(rlist)

!<description>
    ! This subroutine removes the first element from the list
!</description>

!<inputoutput>
    ! The linked list
    type(template_TD(t_list,T,D)), intent(inout) :: rlist
!</inputoutput>
!</subroutine>
    
    ! local variable
    type(template_TD(it_list,T,D)) :: riterator
    
    ! Check if list is empty
    if (list_empty(rlist)) return

    riterator = list_begin(rlist)
    riterator = list_erase(rlist, riterator)

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine template_TD(list_pop_back,T,D)(rlist)

!<description>
    ! This subroutine removes the last element from the list
!</description>

!<inputoutput>
    ! The linked list
    type(template_TD(t_list,T,D)), intent(inout) :: rlist
!</inputoutput>
!</subroutine>

    ! local variable
    type(template_TD(it_list,T,D)) :: riterator
    
    ! Check if list is empty
    if (list_empty(rlist)) return

    riterator = list_rbegin(rlist)
    riterator = list_erase(rlist, riterator)

  end subroutine

  !************************************************************************

!<subroutine>

#ifdef D
  subroutine template_TD(list_front1,T,D)(rlist, p_key, p_data)
#else
  subroutine template_TD(list_front1,T,D)(rlist, p_key)
#endif

!<description>
    ! This subroutine returns pointers to the key and data stored in
    ! the first element
!</description>

!<input>
    ! The linked list
    type(template_TD(t_list,T,D)), intent(in) :: rlist
!</input>

!<output>
    ! key value
    TTYPE(T_TYPE), intent(out), pointer :: p_key

#ifdef D
    ! OPTIONAL: Data
    DTYPE(D_TYPE), dimension(:), intent(out), pointer, optional :: p_data
#endif
!</output>
!</subroutine>

    ! local variable
    type(template_TD(it_list,T,D)) :: riterator

    ! Check if list is empty
    if (list_empty(rlist)) then
      
      nullify(p_key)
#ifdef D
      nullify(p_data)
#endif

    else
      
      riterator = list_begin(rlist)

#ifdef D
      call list_get(rlist, riterator, p_key, p_data)
#else
      call list_get(rlist, riterator, p_key)
#endif

    end if

  end subroutine

  !************************************************************************

!<function>

  function template_TD(list_front2,T,D)(rlist) result(key)

!<description>
    ! This function returns the key stored in the first element
!</description>

!<input>
    ! The linked list
    type(template_TD(t_list,T,D)), intent(in) :: rlist
!</input>

!<result>
    ! key value
    TTYPE(T_TYPE) :: key
!</result>
!</function>

    ! local variable
    type(template_TD(it_list,T,D)) :: riterator

    ! Check if list is empty
    if (list_empty(rlist)) return
    
    riterator = list_begin(rlist)
    key = list_get(rlist, riterator)

  end function

  !************************************************************************

!<subroutine>

#ifdef D
  subroutine template_TD(list_back1,T,D)(rlist, p_key, p_data)
#else
  subroutine template_TD(list_back1,T,D)(rlist, p_key)
#endif

!<description>
    ! This subroutine returns pointers to the key and data stored in
    ! the first element
!</description>

!<input>
    ! The linked list
    type(template_TD(t_list,T,D)), intent(in) :: rlist
!</input>

!<output>
    ! key value
    TTYPE(T_TYPE), intent(out), pointer :: p_key

#ifdef D
    ! OPTIONAL: Data
    DTYPE(D_TYPE), dimension(:), intent(out), pointer, optional :: p_data
#endif
!</output>
!</subroutine>

    ! local variable
    type(template_TD(it_list,T,D)) :: riterator

    ! Check if list is empty
    if (list_empty(rlist)) then
      
      nullify(p_key)
#ifdef D
      nullify(p_data)
#endif

    else
      
      riterator = list_rbegin(rlist)

#ifdef D
      call list_get(rlist, riterator, p_key, p_data)
#else
      call list_get(rlist, riterator, p_key)
#endif

    end if

  end subroutine

  !************************************************************************
  
!<function>

  function template_TD(list_back2,T,D)(rlist) result(key)

!<description>
    ! This function returns the key stored in the last element
!</description>

!<input>
    ! The linked list
    type(template_TD(t_list,T,D)), intent(in) :: rlist
!</input>

!<result>
    ! key value
    TTYPE(T_TYPE) :: key
!</result>
!</function>

    ! local variable
    type(template_TD(it_list,T,D)) :: riterator

    ! Check if list is empty
    if (list_empty(rlist)) return

    riterator = list_rbegin(rlist)
    key = list_get(rlist, riterator)

  end function
  
  !************************************************************************

!<function>

#ifdef D
  function template_TD(list_insert1,T,D)(rlist, rposition, key, data)&
                                         result(riterator)
#else
  function template_TD(list_insert1,T,D)(rlist, rposition, key)&
                                         result(riterator)
#endif

!<description>
    ! This function inserts a new element into the list before the element
    ! at rposition and returns an iterator pointing to the new element
!</description>

!<input>
    ! The iterator
    type(template_TD(it_list,T,D)), intent(in) :: rposition

    ! key value
    TTYPE(T_TYPE), intent(in) :: key

#ifdef D
    ! OPTIONAL: Data
    DTYPE(D_TYPE), dimension(:), intent(in), optional :: data
#endif
!</input>

!<inputoutput>
    ! The linked list
    type(template_TD(t_list,T,D)), intent(inout), target :: rlist
!</inputoutput>

!<result>
    ! The new iterator
    type(template_TD(it_list,T,D)) :: riterator
!</result>
!</function>

    ! local variable
    integer :: ipred,ipos
    
    ! Insert new element at the end of the list?
    if (rposition%ipos .eq. LNULL) then
#ifdef D
      call list_push_back(rlist, key, data)
#else
      call list_push_back(rlist, key)
#endif
      ! Create new iterator
      riterator       = rposition
      riterator%ipos  = rlist%p_Knext(LHEAD)
      riterator%iSpec = iand(riterator%iSpec,&
                             not(LIST_LSPEC_REVERSE))
      ! That`s it
      return
    end if

    ! Get position of predecessor
    ipred = rlist%p_Kprev(rposition%ipos)
    
    ! Insert new element at the beginning of the list?
    if (ipred .eq. LNULL) then
#ifdef D
      call list_push_front(rlist, key, data)
#else
      call list_push_front(rlist, key)
#endif
      ! Create new iterator
      riterator       = rposition
      riterator%ipos  = rlist%p_Knext(LTAIL)
      riterator%iSpec = iand(riterator%iSpec,&
                             not(LIST_LSPEC_REVERSE))
      ! That`s it
      return
    end if

    ! Check if list needs to be resized
    rlist%NA = rlist%NA+1
    ipos     = rlist%p_Knext(LFREE)
    if (abs(ipos) > rlist%NNA) then
      call list_resize(rlist, ceiling(rlist%dfactor*rlist%NNA))
    end if

    ! Set next free position
    if (ipos > 0) then
      rlist%p_Knext(LFREE) = ipos+1
    else
      ipos = abs(ipos)
      rlist%p_Knext(LFREE) = rlist%p_Knext(ipos)
    end if
    
    ! Insert element between its predecessor IPRED and 
    ! the element at position rposition
    rlist%p_Kprev(rposition%ipos) = ipos
    rlist%p_Knext(ipos)  = rlist%p_Knext(ipred)
    rlist%p_Knext(ipred) = ipos
    rlist%p_Kprev(ipos)  = ipred

    ! Store key value
    rlist%p_Key(ipos) = key

#ifdef D
    ! Store data value
    if (present(data)) then
      rlist%p_Data(:,ipos) = data
    end if
#endif

    ! Create new iterator
    riterator       = rposition
    riterator%ipos  = ipos
    riterator%iSpec = iand(riterator%iSpec,&
                           not(LIST_LSPEC_REVERSE))

  end function

  !************************************************************************

!<subroutine>

#ifdef D
  subroutine template_TD(list_insert2,T,D)(rlist, rposition, n, key, data)
#else
  subroutine template_TD(list_insert2,T,D)(rlist, rposition, n, key)
#endif

!<description>
    ! This subroutine inserts n copies of the given key and data
    ! before the element at rposition.
!</description>

!<input>
    ! The iterator
    type(template_TD(it_list,T,D)), intent(in) :: rposition

    ! Number of copies
    integer, intent(in) :: n

    ! key value
    TTYPE(T_TYPE), intent(in) :: key

#ifdef D
    ! OPTIONAL: Data
    DTYPE(D_TYPE), dimension(:), intent(in), optional :: data
#endif
!</input>

!<inputoutput>
    ! The linked list
    type(template_TD(t_list,T,D)), intent(inout) :: rlist
!</inputoutput>
!</subroutine>

    ! local variable
    type(template_TD(it_list,T,D)) :: riterator
    integer :: i

    do i = 1, n
#ifdef D
      riterator = list_insert(rlist, rposition, key, data)
#else
      riterator = list_insert(rlist, rposition, key)
#endif
    end do

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine template_TD(list_insert3,T,D)(rlist, rposition, rfirst, rlast)

!<description>
    ! This subroutine inserts content in the range [rfirst,rlast)
    ! before the element at rposition.
!</description>

!<input>
    ! The iterator
    type(template_TD(it_list,T,D)), intent(in) :: rposition

    ! Iterator referring to the first element
    type(template_TD(it_list,T,D)), intent(in) :: rfirst

    ! Iterator referring to the past-the-last element
    type(template_TD(it_list,T,D)), intent(in) :: rlast
!</input>

!<inputoutput>
    ! The linked list
    type(template_TD(t_list,T,D)), intent(inout) :: rlist
!</inputoutput>
!</subroutine>

    ! local variable
    type(template_TD(it_list,T,D)) :: riterator,riter1
    TTYPE(T_TYPE), pointer :: p_key

#ifdef D
    DTYPE(D_TYPE), dimension(:), pointer :: p_data
#endif

    ! Push content to the list
    riterator = rfirst
    do while(riterator /= rlast)

#ifdef D
      call list_get(riterator%p_rlist, riterator, p_key, p_data)
      riter1 = list_insert(rlist, rposition, p_key, p_data)
#else
      call list_get(riterator%p_rlist, riterator, p_key)
      riter1 = list_insert(rlist, rposition, p_key)
#endif
      call list_next(riterator)
    end do

  end subroutine

  !************************************************************************

!<function>

  function template_TD(list_erase1,T,D)(rlist, rposition) result(riterator)

!<description>
    ! This function removes the element from the list at rposition and
    ! returns an iterator pointing to the element following the
    ! removed element.
!</description>

!<input>
    ! The iterator
    type(template_TD(it_list,T,D)) :: rposition
!</input>

!<inputoutput>
    ! The linked list
    type(template_TD(t_list,T,D)), intent(inout) :: rlist
!</inputoutput>

!<result>
    ! The new iterator
    type(template_TD(it_list,T,D)) :: riterator
!</result>
!</function>

    ! local variables
    integer :: ipred,ipos,inext
    
    ! Delete element
    rlist%NA = rlist%NA-1

    ipos  = rposition%ipos
    ipred = rlist%p_Kprev(ipos)

    ! Are we at the head of the list?
    if (ipred .eq. LNULL) then
      ! Erase first list element
      inext = rlist%p_Knext(ipos)
      rlist%p_Knext(LHEAD) = inext
    else
      ! Erase element inside the list
      inext = rlist%p_Knext(ipos)
      rlist%p_Knext(ipred) = inext
    end if
    
    ! Are we at the tail of the list?
    if (inext .eq. LNULL) then
      rlist%p_Knext(LTAIL) = ipred
    else
      rlist%p_Kprev(inext) = ipred
    end if
    
    ! Update free position
    rlist%p_Knext(ipos)  = rlist%p_Knext(LFREE)
    rlist%p_Knext(LFREE) = -ipos

    ! Create new iterator
    riterator       = rposition
    riterator%ipos  = inext
    riterator%iSpec = iand(riterator%iSpec,&
                           not(LIST_LSPEC_REVERSE))

  end function

  !************************************************************************

!<function>

  function template_TD(list_erase2,T,D)(rlist, rfirst, rlast) result(riterator)

!<description>
    ! This function removes the elements [rfirst,rlast) from the list
    ! and returns an iterator referring to the element following the
    ! last removed element.
!</description>

!<input>
    ! Iterator referring to the first element
    type(template_TD(it_list,T,D)), intent(in) :: rfirst

    ! Iterator referring to the past-the-last element
    type(template_TD(it_list,T,D)), intent(in) :: rlast
!</input>

!<inputoutput>
    ! The linked list
    type(template_TD(t_list,T,D)), intent(inout) :: rlist
!</inputoutput>

!<result>
    ! The new iterator
    type(template_TD(it_list,T,D)) :: riterator
!</result>
!</function>

    ! Remove elements
    riterator = rfirst
    do while(riterator /= rlast)
      riterator = list_erase(rlist, riterator)
    end do

  end function

  !************************************************************************

!<function>

  recursive function template_TD(list_find1,T,D)(rlist, key, rpositionGuess)&
                                                 result(riterator)

!<description>
    ! This function searches for a given key value in the list and
    ! returns its position in the list. If the key value is not
    ! present in the list then NULL is returned.
!</description>

!<input>
    ! key value
    TTYPE(T_TYPE), intent(in) :: key

    ! The linked list
    type(template_TD(t_list,T,D)), intent(in) :: rlist

    ! OPTIONAL: initial guess of the position
    type(template_TD(it_list,T,D)), intent(in), optional :: rpositionGuess
!</input>

!<result>
    ! The iterator
    type(template_TD(it_list,T,D)) :: riterator
!</result>
!</function>

    ! local variable
    logical :: bsuccess

    ! Do we have an initial guess
    if (present(rpositionGuess)) then
      riterator = rpositionGuess
    else
      riterator = list_begin(rlist)
    end if
    
    ! Initialisation
    bsuccess = .false.

    do while(.not.list_isNull(riterator))
      if (rlist%p_Key(riterator%ipos) .eq. key) then
        bsuccess = .true.; exit
      end if
      riterator%ipos = rlist%p_Knext(riterator%ipos)
    end do
    
    if (.not.bsuccess) then
      riterator%ipos = LNULL
      
      ! If we started searching with an initial guess but we were
      ! unable to find the element then re-search from scratch
      if (present(rpositionGuess))&
          riterator = list_find(rlist, key)
    end if

  end function

  !************************************************************************

!<function>

  recursive function template_TD(list_find2,T,D)(rlist, key, bisReverse,&
                                                 rpositionGuess) result(riterator)

!<description>
    ! This function searches for a given key value in a sorted list
    ! and returns its position in the list. If the key value is not
    ! present in the list then riterator refers to the element which
    ! lexicographically follows the key which was not found.
    !
    ! Note that this function requires the list to be sorted without
    ! checking this. Setting the parameter bisReverse to .true. 
    ! indicates that the list is sorted in reverse order.
!</description>

!<input>
    ! key value
    TTYPE(T_TYPE), intent(in) :: key

    ! The linked list
    type(template_TD(t_list,T,D)), intent(in) :: rlist

    ! Flag: if TRUE then list is assumed in reverse order
    logical, intent(in) :: bisReverse

    ! OPTIONAL: initial guess of the position
    type(template_TD(it_list,T,D)), intent(in), optional :: rpositionGuess
!</input>

!<result>
    ! The iterator
    type(template_TD(it_list,T,D)) :: riterator
!</result>
!</function>

    ! local variable
    logical :: bsuccess

    ! Do we have an initial guess
    if (present(rpositionGuess)) then
      riterator = rpositionGuess
    else
      riterator = list_begin(rlist)
    end if
    
    ! Initialisation
    bsuccess = .false.


    if (bisReverse) then
      ! Search for key value in reverse order
      do while(.not.list_isNull(riterator))
        if (rlist%p_Key(riterator%ipos) .eq. key) then
          bsuccess = .true.; exit
        elseif (rlist%p_Key(riterator%ipos) .lt. key) then
          bsuccess = .false.; exit
        end if
        riterator%ipos = rlist%p_Knext(riterator%ipos)
      end do
    else
      ! Search for key value in default order
      do while(.not.list_isNull(riterator))
        if (rlist%p_Key(riterator%ipos) .eq. key) then
          bsuccess = .true.; exit
        elseif (rlist%p_Key(riterator%ipos) .gt. key) then
          bsuccess = .false.; exit
        end if
        riterator%ipos = rlist%p_Knext(riterator%ipos)
      end do     
    end if

    if (.not.bsuccess) then
      ! Special treatment of position
      riterator%iSpec = riterator%iSpec + LIST_LSPEC_VIRTUAL
      
      ! If we started searching with an initial guess but we were
      ! unable to find the element then re-search from scratch
      if (present(rpositionGuess))&
          riterator = list_find(rlist, key, bisReverse)
    end if
    
  end function

  !************************************************************************

  subroutine template_TD(list_print,T,D)(rlist)

!<description>
    ! This subroutine prints the content of the list
!</description>

!<input>
    ! The linked list
    type(template_TD(t_list,T,D)), intent(in) :: rlist
!</input>
!</subroutine>

#ifdef T_STORAGE

    ! local variable
    type(template_TD(it_list,T,D)) :: riterator

    riterator = list_begin(rlist)
    do while (riterator /= list_end(rlist))
      write(*,*) rlist%p_Key(riterator%ipos)
      call list_next(riterator)
    end do
    
#else
    
    call output_line('Unable to print list with derived data type!',&
        OU_CLASS_WARNING,OU_MODE_STD,'list_print')

#endif

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine template_TD(list_info,T,D)(rlist)

!<description>
    ! This subroutine prints information about the list
!</description>

!<input>
    ! list
    type(template_TD(t_list,T,D)), intent(in) :: rlist
!</input>
!</subroutine>

    call output_line('List:')
    call output_line('----------')
    call output_line('NA:       '//trim(sys_siL(rlist%NA,15)))
    call output_line('NNA:      '//trim(sys_siL(rlist%NNA,15)))
    call output_line('NNA0:     '//trim(sys_siL(rlist%NNA0,15)))
    call output_line('nresize:  '//trim(sys_siL(rlist%nresize,15)))
    call output_line('dfactor:  '//trim(sys_sdL(rlist%dfactor,2)))
    call output_line('h_Knext:  '//trim(sys_siL(rlist%h_Knext,15)))
    call output_line('h_Kprev:  '//trim(sys_siL(rlist%h_Kprev,15)))

#ifdef T_STORAGE
    call output_line('h_Key:    '//trim(sys_siL(rlist%h_Key,15)))
#endif
#if defined D && defined D_STORAGE
    call output_line('h_Data:   '//trim(sys_siL(rlist%h_Data,15)))
#endif

    call output_lbrk()

    call output_line('Current data  memory usage: '//&
        trim(sys_sdL(100*rlist%NA/&
        real(max(1,rlist%NNA),DP),2))//'%')
    call output_line('Total data memory usage:    '//&
        trim(sys_sdL(100*rlist%NNA/&
        real(max(1,rlist%NNA0),DP),2))//'%')

    call output_lbrk()

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine template_TD(list_duplicate,T,D)(rlist, rlistBackup)

!<description>
    ! This subroutine makes a copy of an list in memory.
    ! It does not make sense to share some information between lists,
    ! so each  is physically copied from the source rlist
    ! to the destination rlistBackup.
!</description>

!<input>
    ! The rraylist for which a backup copy should be generated
    type(template_TD(t_list,T,D)), intent(in) :: rlist
!</input>

!<inputoutput>
    ! Backup copy list
    type(template_TD(t_list,T,D)), intent(inout) :: rlistBackup
!</inputoutput>
!</subroutine>

    ! Release backup list
    call list_release(rlistBackup)

    ! Copy all data
    rlistBackup = rlist

    ! Reset handles
    rlistBackup%h_Knext = ST_NOHANDLE
    rlistBackup%h_Kprev = ST_NOHANDLE

    ! Copy storage blocks
    if (rlist%h_Kprev .ne. ST_NOHANDLE) then
      call storage_copy(rlist%h_Kprev, rlistBackup%h_Kprev)
      call storage_getbase_int(rlistBackup%h_Kprev,&
          rlistBackup%p_Kprev)
    end if

    if (rlist%h_Knext .ne. ST_NOHANDLE) then
      call storage_copy(rlist%h_Knext, rlistBackup%h_Knext)
      call storage_getbase_int(rlistBackup%h_Knext,&
                               rlistBackup%p_Knext)
    end if

#ifdef T_STORAGE
    rlistBackup%h_Key = ST_NOHANDLE
    if (rlist%h_Key .ne. ST_NOHANDLE) then
      call storage_copy(rlist%h_Key, rlistBackup%h_Key)
      call storage_getbase(rlistBackup%h_Key,&
                           rlistBackup%p_Key)
    end if
#else
    if (associated(rlist%p_Key)) then
      allocate(rlistBackup%p_Key(size(rlist%p_Key)))
      rlistBackup%p_Key = rlist%p_Key
    else
      nullify(rlistBackup%p_Key)
    end if
#endif

#ifdef D
#ifdef D_STORAGE
    rlistBackup%h_Data = ST_NOHANDLE
    if (rlist%h_Data .ne. ST_NOHANDLE) then
      call storage_copy(rlist%h_Data, rlistBackup%h_Data)
      call storage_getbase(rlistBackup%h_Data,&
                           rlistBackup%p_Data)
    end if
#else
    if (associated(rlist%p_Data)) then
      allocate(rlistBackup%p_Data(size(rlist%p_Data,1),&
                                  size(rlist%p_Data,2)))
      rlistBackup%p_Data = rlist%p_Data
    else
      nullify(rlistBackup%p_Data)
    end if
#endif
#endif

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine template_TD(list_restore,T,D)(rlistBackup, rlist)

!<description>
    ! This subroutine restores an list from a previous backup.
    ! The format and ordering of both lists must be the same.
!</description>

!<input>
    ! Backup copy of the list
    type(template_TD(t_list,T,D)), intent(in) :: rlistBackup
!</input>

!<inputoutput>
    ! Destination list
    type(template_TD(t_list,T,D)), intent(inout) :: rlist
!</inputoutput>
!</subroutine>

    ! Release list
    call list_release(rlist)
    
    ! Duplicate the backup
    call list_duplicate(rlistBackup, rlist)

  end subroutine

  !************************************************************************

!<function>

  pure function template_TD(list_size,T,D)(rlist) result(isize)

!<description>
    ! Returns the size of the list
!</description>

!<input>
    ! The linked list
    type(template_TD(t_list,T,D)), intent(in) :: rlist
!</input>

!<result>
    ! Size of the list
    integer :: isize
!</result>
!</function>

    isize = rlist%NA

  end function

  !************************************************************************

!<function>

  pure function template_TD(list_max_size,T,D)(rlist) result(imaxsize)

!<description>
    ! Returns the maximum size of the list
!</description>

!<input>
    ! The linked list
    type(template_TD(t_list,T,D)), intent(in) :: rlist
!</input>

!<result>
    ! Maximum size of the list
    integer :: imaxsize
!</result>
!</function>

    imaxsize = rlist%NNA

  end function

  !************************************************************************

!<function>

  pure function template_TD(list_empty,T,D)(rlist) result(bempty)

!<description>
    ! Checks if the list is empty
!</description>

!<input>
    ! The linked list
    type(template_TD(t_list,T,D)), intent(in) :: rlist
!</input>

!<result>
    ! Flag: is true if the list is empty
    logical :: bempty
!</result>
!</function>

    bempty = (rlist%NA .eq. 0)

  end function

  !************************************************************************

!<subroutine>

  subroutine template_TD(list_fassign,T,D)(rlistDest, rlistSrc)

!<description>
    ! Assigns the content of rlistSrc to rlistDest
!</description>

!<input>
    ! Source list
    type(template_TD(t_list,T,D)), intent(in) :: rlistSrc
!</input>

!<output>
    ! Destination list
    type(template_TD(t_list,T,D)), intent(out) :: rlistDest
!</output>
!</subroutine>

    ! local variable
    integer :: i

    ! Create empty list
#ifdef D
    call list_create(rlistDest, rlistSrc%NNA, rlistSrc%isizeData,&
        rlistSrc%dfactor)
#else
    call list_create(rlistDest, rlistSrc%NNA, rlistSrc%dfactor)
#endif

    ! Set size
    rlistDest%NA = rlistSrc%NA
    
    ! Set structure
    rlistDest%p_Knext(LHEAD) = rlistSrc%p_Knext(LHEAD)
    rlistDest%p_Knext(LTAIL) = rlistSrc%p_Knext(LTAIL)
    rlistDest%p_Knext(LFREE) = rlistSrc%p_Knext(LFREE)

    do i = 1, rlistSrc%NA
      rlistDest%p_Knext(i) = rlistSrc%p_Knext(i)
      rlistDest%p_Kprev(i) = rlistSrc%p_Kprev(i)
      rlistDest%p_Key(i)   = rlistSrc%p_Key(i)
    end do
    
#ifdef D
    do i = 1, rlistSrc%NA
      rlistDest%p_Data(:,i) = rlistSrc%p_Data(:,i)
    end do
#endif
    
  end subroutine

  !************************************************************************

!<subroutine>

  subroutine template_TD(list_reverse,T,D)(rlist)

!<description>
    ! This subroutine reverses the ordering of the list
!</description>

!<inputoutput>
    ! The linked list to be reverted
    type(template_TD(t_list,T,D)), intent(inout) :: rlist
!</inputoutput>
!</subroutine>

    ! local variable
    type(template_TD(it_list,T,D)) :: riterator
    TTYPE(T_TYPE), pointer :: p_key

#ifdef D
    DTYPE(D_TYPE), dimension(:), pointer :: p_data
#endif

    riterator = list_begin(rlist)
    call list_next(riterator)
    do while (riterator /= list_end(rlist))

#ifdef D
      call list_get(rlist, riterator, p_key, p_data)
      call list_push_front(rlist, p_key, p_data)
#else
      call list_get(rlist, riterator, p_key)
      call list_push_front(rlist, p_key)
#endif
      riterator = list_erase(rlist, riterator)
    end do

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine template_TD(list_sort,T,D)(rlist)

!<description>
    ! This subroutine sorts the elements in the list
!</description>

!<inputoutput>
    ! The linked list to be sorted
    type(template_TD(t_list,T,D)), intent(inout) :: rlist
!</inputoutput>
!</subroutine>

    ! local variable
    logical :: bswapped
    type(template_TD(it_list,T,D)) :: rposition,riterator
    TTYPE(T_TYPE), pointer :: p_key1,p_key2
    TTYPE(T_TYPE) :: key

#ifdef D
    DTYPE(D_TYPE), dimension(:), pointer :: p_data1,p_data2
    DTYPE(D_TYPE), dimension(:), allocatable :: data

    ! Allocate temporal memory
    if (rlist%isizeData > 0) allocate(data(rlist%isizeData))
#endif

    rposition = list_rbegin(rlist)

    ! Bubble sort algorithm
    outer: do
      bswapped  = .false.
      riterator = list_begin(rlist)
      inner: do while (riterator /= rposition)
        
#ifdef D
        if (rlist%isizeData > 0) then

          call list_get(rlist,riterator,p_key1,p_data1)
          call list_next(riterator)
          call list_get(rlist,riterator,p_key2,p_data2)
          if (p_key1 > p_key2) then
            key    = p_key2
            p_key2 = p_key1
            p_key1 = key
            data    = p_data2
            p_data2 = p_data1
            p_data1 = data
            bswapped = .true.
          end if

        else

          call list_get(rlist,riterator,p_key1)
          call list_next(riterator)
          call list_get(rlist,riterator,p_key2)
          if (p_key1 > p_key2) then
            key    = p_key2
            p_key2 = p_key1
            p_key1 = key
            bswapped = .true.
          end if

        end if
#else
        call list_get(rlist,riterator,p_key1)
        call list_next(riterator)
        call list_get(rlist,riterator,p_key2)
        if (p_key1 > p_key2) then
          key    = p_key2
          p_key2 = p_key1
          p_key1 = key
          bswapped = .true.
        end if
#endif

      end do inner
      call list_next(rposition)
      if (.not.bswapped) exit
    end do outer

#ifdef D
    if (rlist%isizeData > 0) deallocate(data)
#endif

  end subroutine

  !************************************************************************

!<function>

  pure function template_TD(list_isNull,T,D)(riterator) result(bisNull)

!<description>
    ! Checks if the iterator is NULL
!</description>

!<input>
    ! Iterator
    type(template_TD(it_list,T,D)), intent(in) :: riterator
!</input>

!<result>
    ! True if the iterator is NULL.
    logical :: bisNull
!</result>
!</function>

    bisNull = ((riterator%ipos .eq. LNULL) .or.&
           iand(riterator%iSpec, LIST_LSPEC_VIRTUAL) .ne. 0)

  end function

  !************************************************************************

!<function>

  pure function template_TD(list_hasSpec,T,D)(riterator,iSpec) result(bhasSpec)

!<description>
    ! Checks if the iterator has the given specification flag
!</description>

!<input>
    ! Iterator
    type(template_TD(it_list,T,D)), intent(in) :: riterator

    ! Specification flag
    integer(I32), intent(in) :: iSpec
!</input>

!<result>
    ! True if the iterator has given specification flag
    logical :: bhasSpec
!</result>
!</function>

    bhasSpec = (iand(riterator%iSpec, iSpec) .ne. 0)

  end function

  !************************************************************************

!<function>

  pure function template_TD(it_list_eq,T,D)(riterator1,riterator2) result(beq)

!<description>
    ! Compare two iterators for equality
!</description>

!<input>
    ! Iterators
    type(template_TD(it_list,T,D)), intent(in) :: riterator1,riterator2
!</input>

!<result>
    ! True if both iterators are equal
    logical :: beq
!</result>
!</function>

    beq = (riterator1%ipos == riterator2%ipos)

  end function

  !************************************************************************

!<function>

  pure function template_TD(it_list_ne,T,D)(riterator1,riterator2) result(bne)

!<description>
    ! Compare two iterators for inequality
!</description>

!<input>
    ! Iterators
    type(template_TD(it_list,T,D)), intent(in) :: riterator1,riterator2
!</input>

!<result>
    ! True if both iterators are not equal
    logical :: bne
!</result>
!</function>

    bne = (riterator1%ipos /= riterator2%ipos)

  end function

  !************************************************************************

!<function>

  pure function template_TD(it_list_lt,T,D)(riterator1,riterator2) result(blt)

!<description>
    ! Checks lexicographical ordering of two iterators
!</description>

!<input>
    ! Iterators
    type(template_TD(it_list,T,D)), intent(in) :: riterator1,riterator2
!</input>

!<result>
    ! True if the lexicographical ordering of iterator1 is smaller
    ! than that of iterator2
    logical :: blt
!</result>
!</function>

    blt = (riterator1%ipos < riterator2%ipos)

  end function

  !************************************************************************

!<function>

  pure function template_TD(it_list_le,T,D)(riterator1,riterator2) result(ble)

!<description>
    ! Checks lexicographical ordering of two iterators
!</description>

!<input>
    ! Iterators
    type(template_TD(it_list,T,D)), intent(in) :: riterator1,riterator2
!</input>

!<result>
    ! True if the lexicographical ordering of iterator1 is smaller
    ! than or equal to that of iterator2
    logical :: ble
!</result>
!</function>

    ble = (riterator1%ipos <= riterator2%ipos)

  end function

  !************************************************************************

!<function>

  pure function template_TD(it_list_gt,T,D)(riterator1,riterator2) result(bgt)

!<description>
    ! Checks lexicographical ordering of two iterators
!</description>

!<input>
    ! Iterators
    type(template_TD(it_list,T,D)), intent(in) :: riterator1,riterator2
!</input>

!<result>
    ! True if the lexicographical ordering of iterator1 is greater
    ! than that of iterator2
    logical :: bgt
!</result>
!</function>

    bgt = (riterator1%ipos > riterator2%ipos)

  end function

  !************************************************************************

!<function>

  pure function template_TD(it_list_ge,T,D)(riterator1,riterator2) result(bge)

!<description>
    ! Checks lexicographical ordering of two iterators
!</description>

!<input>
    ! Iterators
    type(template_TD(it_list,T,D)), intent(in) :: riterator1,riterator2
!</input>

!<result>
    ! True if the lexicographical ordering of iterator1 is greater
    ! than or equal to that of iterator2
    logical :: bge
!</result>
!</function>

    bge = (riterator1%ipos >= riterator2%ipos)

  end function

#endif
