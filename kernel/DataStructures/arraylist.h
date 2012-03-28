#ifndef _ARRAYLIST_H_
#define _ARRAYLIST_H_

!##############################################################################
!# ****************************************************************************
!# <name> FEAT2_PP_TEMPLATE_TD(arraylist,T,D) </name>
!# ****************************************************************************
!#
!# <purpose>
!#
!# This module header file implements an array of (linear) linked lists, also
!# known as arraylist, for T_TYPE key values and D_TYPE auxiliary data.
!#
!# -----------------------------------------------------------------------------
!# This module is designed with greatest compatibility to the C++
!# Standard Template Library (STL) concerning naming conventions.
!# -----------------------------------------------------------------------------
!#
!# The following routines are available:
!#
!# 1.) alst_create / alst_createTbl (constructor in STL)
!#     -> Creates an empty arraylist /
!#        Creates an empty table in an arraylist
!#
!# 2.) alst_release / alst_releaseTbl (destructor in STL)
!#     -> Releases an arraylist /
!#        Releases a table in an arraylist
!#
!# 3.) alst_resize / alst_resizeTbl (resize in STL)
!#     -> Reallocates memory for an arraylist /
!#        Reallocates memory for a table in an arraylist
!#
!# 4.) alst_copy / alst_copyTbl (no equivalence in STL)
!#     -> Copies data to/from an arraylist /
!#        Copies data to/from a table of an arraylist
!#
!# 5.) alst_clear / alst_clearTbl (clear in STL)
!#     -> Removes all content from an arraylist /
!#        Remvoes all content from a table in an arraylist
!#
!# 6.) alst_swap / alst_swapTbl (swap in STL)
!#     -> Swaps content of two different arraylists /
!#        Swaps content of two different tables in an arraylist
!#
!# 7.) alst_begin / alst_rbegin (begin / rbegin in STL)
!#      -> Returns a (reverse) iterator referring to the first/last element 
!#         in an arraylist
!#
!# 8.) alst_end / alst_rend (end / rendin STL)
!#      -> Returns a (reverse) iterator referring to the element right after
!#         the last/right before the first element in an arraylist
!#
!# 9.) alst_next (++iterator in STL)
!#      -> Increases the (reverse) iterator
!#
!# 10.) alst_prior (--iterator in STL)
!#      -> Decreases the (reverse) iterator
!#
!# 11.) alst_get (no equivalence in STL)
!#      -> Gets key value at given position in arraylist
!#
!# 12.) alst_getbase_key / alst_getbase_data (no equivalence in STL)
!#      -> Gets pointer to key value / auxiliary data value
!#         at given position in arraylist
!#
!# 13.) alst_assign (assign in STL)
!#      -> Assigns key and value data to an arraylist dropping existing content
!#
!# 14.) alst_push_front (push_front in STL)
!#     -> Pushes data to the front of an arraylist
!#
!# 15.) alst_push_back (push_back in STL)
!#      -> Pushes data to the back of an arraylist
!#
!# 16.) alst_pop_front (pop_front in STL)
!#      -> Removes the first element from an arraylist
!#
!# 17.) alst_pop_back (pop_back in STL)
!#      -> Removes the last element from an arraylist
!#
!# 18.) alst_front / alst_getbase_front(front in STL)
!#      -> Gets key and value of the first element in an arraylist
!#
!# 19.) alst_back / alst_getbase_back (back in STL)
!#      -> Gets key and value of the last element in an arraylist
!#
!# 20.) alst_insert (insert in STL)
!#      -> Inserts data into an arraylist
!#
!# 21.) alst_erase (erase in STL)
!#      -> Deletes data from an arraylist
!#
!# 22.) alst_size (size in STL)
!#      -> Returns the size of an arraylist
!#
!# 23.) alst_max_size (max_size in STL)
!#      -> Returns the maximum size of an arraylist
!#
!# 24.) alst_ntable (no equivalence in STL)
!#      -> Returns the number of tables of an arraylist
!#
!# 25.) alst_empty / alst_emptyTbl (empty in STL)
!#      -> Tests if the arraylist is empty
!#
!# 26.) alst_find (no equivalence in STL)
!#      -> Searches for data in arraylist
!#
!# 27.) alst_print (no equivalence in STL)
!#      -> Prints content of arraylist
!#
!# 28.) alst_info (no equivalence in STL)
!#      -> Prints information about arraylist
!#
!# 29.) alst_duplicate (no equivalence in STL)
!#      -> Creates a backup of an arraylist
!#
!# 30.) alst_restore (no equivalence in STL)
!#      -> Restores an arraylist from a previous backup
!#
!# 31.) alst_reverse alst_reverseTbl (reverse in STL)
!#      -> Reverses the order of elements in a list
!#
!# 32.) alst_sort / alst_sortTbl (sort in STL)
!#      -> Sorts the entries in an arraylist
!#
!# 33.) alst_isNull
!#      -> Tests if the iterator is NULL
!#
!# 34.) alst_hasSpec
!#      -> Tests if the iterator has specification flag
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
!# </purpose>
!##############################################################################

#include "../template.h"

  implicit none

  private
  public :: FEAT2_PP_TEMPLATE_TD(t_arraylist,T,D)
  public :: FEAT2_PP_TEMPLATE_TD(it_arraylist,T,D)
  public :: alst_create,alst_createTbl
  public :: alst_release,alst_releaseTbl
  public :: alst_resize,alst_resizeTbl
  public :: alst_clear,alst_clearTbl
  public :: alst_copy,alst_copyTbl
  public :: alst_swap,alst_swapTbl
  public :: alst_begin
  public :: alst_rbegin
  public :: alst_end
  public :: alst_rend
  public :: alst_next
  public :: alst_prior
  public :: alst_get
  public :: alst_getbase_key
#ifdef D
  public :: alst_getbase_data
#endif
  public :: alst_assign
  public :: alst_push_front
  public :: alst_push_back
  public :: alst_pop_front
  public :: alst_pop_back
  public :: alst_front
  public :: alst_getbase_front
  public :: alst_back
  public :: alst_getbase_back
  public :: alst_insert
  public :: alst_erase
  public :: alst_size
  public :: alst_max_size
  public :: alst_ntable
  public :: alst_empty,alst_emptyTbl
  public :: alst_find
  public :: alst_print,alst_printTbl
  public :: alst_info
  public :: alst_duplicate
  public :: alst_restore
  public :: alst_reverse,alst_reverseTbl
  public :: alst_sort,alst_sortTbl
  public :: alst_isNull
  public :: alst_hasSpec

  public assignment(=)
  public operator(==)
  public operator(/=)
  public operator(<)
  public operator(<=)
  public operator(>)
  public operator(>=)

  interface alst_create
    module procedure FEAT2_PP_TEMPLATE_TD(alst_create,T,D)
    module procedure FEAT2_PP_TEMPLATE_TD(alst_createTbl,T,D)
  end interface

  interface alst_createTbl
    module procedure FEAT2_PP_TEMPLATE_TD(alst_createTbl,T,D)
  end interface

  interface alst_release
    module procedure FEAT2_PP_TEMPLATE_TD(alst_release,T,D)
    module procedure FEAT2_PP_TEMPLATE_TD(alst_releaseTbl,T,D)
  end interface

  interface alst_releaseTbl
    module procedure FEAT2_PP_TEMPLATE_TD(alst_releaseTbl,T,D)
  end interface

  interface alst_resize
    module procedure FEAT2_PP_TEMPLATE_TD(alst_resize,T,D)
  end interface

  interface alst_resizeTbl
    module procedure FEAT2_PP_TEMPLATE_TD(alst_resizeTbl,T,D)
  end interface

  interface alst_clear
    module procedure FEAT2_PP_TEMPLATE_TD(alst_clear,T,D)
    module procedure FEAT2_PP_TEMPLATE_TD(alst_clearTbl,T,D)
  end interface

  interface alst_clearTbl
    module procedure FEAT2_PP_TEMPLATE_TD(alst_clearTbl,T,D)
  end interface

  interface alst_copy
    module procedure FEAT2_PP_TEMPLATE_TD(alst_cpy1,T,D)
    module procedure FEAT2_PP_TEMPLATE_TD(alst_cpy2,T,D)
    module procedure FEAT2_PP_TEMPLATE_TD(alst_cpy3,T,D)
    module procedure FEAT2_PP_TEMPLATE_TD(alst_cpy4,T,D)
  end interface

  interface alst_copyTbl
    module procedure FEAT2_PP_TEMPLATE_TD(alst_cpy1Tbl,T,D)
    module procedure FEAT2_PP_TEMPLATE_TD(alst_cpy2Tbl,T,D)
    module procedure FEAT2_PP_TEMPLATE_TD(alst_cpy3Tbl,T,D)
    module procedure FEAT2_PP_TEMPLATE_TD(alst_cpy4Tbl,T,D)
  end interface

  interface alst_swap
    module procedure FEAT2_PP_TEMPLATE_TD(alst_swap,T,D)
    module procedure FEAT2_PP_TEMPLATE_TD(alst_swapTbl,T,D)
  end interface

  interface alst_swapTbl
    module procedure FEAT2_PP_TEMPLATE_TD(alst_swapTbl,T,D)
  end interface

  interface alst_begin
    module procedure FEAT2_PP_TEMPLATE_TD(alst_begin,T,D)
  end interface

  interface alst_rbegin
    module procedure FEAT2_PP_TEMPLATE_TD(alst_rbegin,T,D)
  end interface

  interface alst_end
    module procedure FEAT2_PP_TEMPLATE_TD(alst_end,T,D)
  end interface

  interface alst_rend
    module procedure FEAT2_PP_TEMPLATE_TD(alst_rend,T,D)
  end interface

  interface alst_next
    module procedure FEAT2_PP_TEMPLATE_TD(alst_next,T,D)
  end interface

  interface alst_prior
    module procedure FEAT2_PP_TEMPLATE_TD(alst_prior,T,D)
  end interface

  interface alst_get
    module procedure FEAT2_PP_TEMPLATE_TD(alst_get,T,D)
  end interface

  interface alst_getbase_key
    module procedure FEAT2_PP_TEMPLATE_TD(alst_getbase_key,T,D)
  end interface

#ifdef D
  interface alst_getbase_data
    module procedure FEAT2_PP_TEMPLATE_TD(alst_getbase_data,T,D)
  end interface
#endif

  interface alst_assign
    module procedure FEAT2_PP_TEMPLATE_TD(alst_assign1,T,D)
    module procedure FEAT2_PP_TEMPLATE_TD(alst_assign2,T,D)
  end interface

  interface alst_push_front
    module procedure FEAT2_PP_TEMPLATE_TD(alst_push_front,T,D)
  end interface

  interface alst_push_back
    module procedure FEAT2_PP_TEMPLATE_TD(alst_push_back,T,D)
  end interface

  interface alst_pop_front
    module procedure FEAT2_PP_TEMPLATE_TD(alst_pop_front,T,D)
  end interface

  interface alst_pop_back
    module procedure FEAT2_PP_TEMPLATE_TD(alst_pop_back,T,D)
  end interface

  interface alst_front
    module procedure FEAT2_PP_TEMPLATE_TD(alst_front,T,D)
  end interface

  interface alst_getbase_front
    module procedure FEAT2_PP_TEMPLATE_TD(alst_getbase_front,T,D)
  end interface

  interface alst_back
    module procedure FEAT2_PP_TEMPLATE_TD(alst_back,T,D)
  end interface

  interface alst_getbase_back
    module procedure FEAT2_PP_TEMPLATE_TD(alst_getbase_back,T,D)
  end interface

  interface alst_insert
    module procedure FEAT2_PP_TEMPLATE_TD(alst_insert1,T,D)
    module procedure FEAT2_PP_TEMPLATE_TD(alst_insert2,T,D)
    module procedure FEAT2_PP_TEMPLATE_TD(alst_insert3,T,D)
  end interface

  interface alst_erase
    module procedure FEAT2_PP_TEMPLATE_TD(alst_erase1,T,D)
    module procedure FEAT2_PP_TEMPLATE_TD(alst_erase2,T,D)
  end interface

  interface alst_size
    module procedure FEAT2_PP_TEMPLATE_TD(alst_size,T,D)
  end interface

  interface alst_max_size
    module procedure FEAT2_PP_TEMPLATE_TD(alst_max_size,T,D)
  end interface

  interface alst_ntable
    module procedure FEAT2_PP_TEMPLATE_TD(alst_ntable,T,D)
  end interface

  interface alst_empty
    module procedure FEAT2_PP_TEMPLATE_TD(alst_empty,T,D)
    module procedure FEAT2_PP_TEMPLATE_TD(alst_emptyTbl,T,D)
  end interface

  interface alst_emptyTbl
    module procedure FEAT2_PP_TEMPLATE_TD(alst_emptyTbl,T,D)
  end interface

  interface alst_find
    module procedure FEAT2_PP_TEMPLATE_TD(alst_find1,T,D)
    module procedure FEAT2_PP_TEMPLATE_TD(alst_find2,T,D)
  end interface

  interface alst_print
    module procedure FEAT2_PP_TEMPLATE_TD(alst_print,T,D)
    module procedure FEAT2_PP_TEMPLATE_TD(alst_printTbl,T,D)
  end interface

  interface alst_printTbl
    module procedure FEAT2_PP_TEMPLATE_TD(alst_printTbl,T,D)
  end interface

  interface alst_info
    module procedure FEAT2_PP_TEMPLATE_TD(alst_info,T,D)
  end interface

  interface alst_duplicate
    module procedure FEAT2_PP_TEMPLATE_TD(alst_duplicate,T,D)
  end interface

  interface alst_restore
    module procedure FEAT2_PP_TEMPLATE_TD(alst_restore,T,D)
  end interface

  interface alst_reverse
    module procedure FEAT2_PP_TEMPLATE_TD(alst_reverse,T,D)
    module procedure FEAT2_PP_TEMPLATE_TD(alst_reverseTbl,T,D)
  end interface

  interface alst_reverseTbl
    module procedure FEAT2_PP_TEMPLATE_TD(alst_reverseTbl,T,D)
  end interface

  interface alst_sort
    module procedure FEAT2_PP_TEMPLATE_TD(alst_sort,T,D)
    module procedure FEAT2_PP_TEMPLATE_TD(alst_sortTbl,T,D)
  end interface

  interface alst_sortTbl
    module procedure FEAT2_PP_TEMPLATE_TD(alst_sortTbl,T,D)
  end interface

  interface alst_isNull
    module procedure FEAT2_PP_TEMPLATE_TD(alst_isNull,T,D)
  end interface

  interface alst_hasSpec
    module procedure FEAT2_PP_TEMPLATE_TD(alst_hasSpec,T,D)
  end interface


  interface assignment(=)
    module procedure FEAT2_PP_TEMPLATE_TD(alst_fassign,T,D)
  end interface


  interface operator(==)
    module procedure FEAT2_PP_TEMPLATE_TD(it_alst_eq,T,D)
  end interface

  interface operator(/=)
    module procedure FEAT2_PP_TEMPLATE_TD(it_alst_ne,T,D)
  end interface

  interface operator(<)
    module procedure FEAT2_PP_TEMPLATE_TD(it_alst_lt,T,D)
  end interface

  interface operator(<=)
    module procedure FEAT2_PP_TEMPLATE_TD(it_alst_le,T,D)
  end interface

  interface operator(>)
    module procedure FEAT2_PP_TEMPLATE_TD(it_alst_gt,T,D)
  end interface

  interface operator(>=)
    module procedure FEAT2_PP_TEMPLATE_TD(it_alst_ge,T,D)
  end interface

  !************************************************************************

!<types>

!<typeblock>

  ! An arraylist

  type FEAT2_PP_TEMPLATE_TD(t_arraylist,T,D)
    private

    ! Number of tables that are currently stored in the arraylist
    integer :: ntable = 0

    ! Total number of tables that can be stored in the arraylist
    integer :: nntable = 0

    ! Total number of tables that can initially be stored in the
    ! arraylist. This information is needed to compute the growth of
    ! the arraylist after several resize operations
    integer :: nntable0 = 0

    ! Number of items that are currently stored in all lists
    integer :: NA = 0

    ! Total number of items that can be stored in all lists
    integer :: NNA = 0

    ! Total number of items that can initially be stored in all
    ! lists. This information is needed to compute the growth of the
    ! arraylist after several resize operations
    integer :: NNA0 = 0

    ! Number of resize operations performed with the arraylist
    integer :: nresize = 0

    ! Factor by which the arraylist is enlarged if new storage is allocate
    real(DP) :: dfactor = 1.5_DP

    ! Handle and pointer to the lookup table
    integer :: h_Ktable = ST_NOHANDLE
    integer, dimension(:,:), pointer :: p_Ktable => null()

    ! Handle and pointer to the arraylist structure
    integer :: h_Knext = ST_NOHANDLE
    integer :: h_Kprev = ST_NOHANDLE
    integer, dimension(:), pointer ::   p_Knext => null()
    integer, dimension(:), pointer ::   p_Kprev => null()

#ifdef T_STORAGE
    ! Handle to the arraylist key
    integer :: h_Key = ST_NOHANDLE
#endif

    ! Pointer to the key values of the arraylist 
    FEAT2_PP_TTYPE(T_TYPE), dimension(:), pointer :: p_Key => null()

#ifdef D
    ! Dimension of the auxiliary data values to be stored
    integer :: isizeData = 0
#ifdef D_STORAGE
    ! Handle to the list auxiliary data
    integer :: h_Data = ST_NOHANDLE
#endif

    ! Pointer to arraylist auxiliary data
    FEAT2_PP_DTYPE(D_TYPE), dimension(:,:), pointer :: p_Data => null()
#endif

  end type

!</typeblock>

!<typeblock>

  ! An iterator for an arraylist

  type FEAT2_PP_TEMPLATE_TD(it_arraylist,T,D)
    private

    ! Number of the table
    integer :: itable = 0

    ! Absolute position of the current element
    integer :: ipos = ALST_NULL

    ! Specification flag. This is a bitfield coming from an OR
    ! combination of different ALST_LSPEC_xxxx constants and specifies
    ! various details of the list iterator.
    integer(I32) :: iSpec = 0_I32

    ! Pointer to the underlying arraylist
    type(FEAT2_PP_TEMPLATE_TD(t_arraylist,T,D)), pointer :: p_rarraylist => null()

  end type

!</typeblock>

!</types>

contains

  !************************************************************************

!<subroutine>

#ifdef D
  subroutine FEAT2_PP_TEMPLATE_TD(alst_create,T,D)(rarraylist, nntable, NNA,&
      isizeData, dfactor)
#else
  subroutine FEAT2_PP_TEMPLATE_TD(alst_create,T,D)(rarraylist, nntable, NNA, dfactor)
#endif

!<description>
    ! This subroutine creates a new arraylist
!</description>

!<input>
    ! Total number of tables
    integer, intent(in) :: nntable

    ! Total number of items that can be stored in all lists
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
    ! The arraylist
    type(FEAT2_PP_TEMPLATE_TD(t_arraylist,T,D)), intent(out) :: rarraylist
!</output>
!</subroutine>

    ! local variables
    integer, dimension(2) :: Isize

    ! Set factor
    if (present(dfactor)) then
      if (dfactor > 1.0_DP) rarraylist%dfactor=dfactor
    end if

    ! Initialise list
    rarraylist%nntable  = max(0,nntable)
    rarraylist%nntable0 = max(0,nntable)
    rarraylist%NNA      = max(0,NNA)
    rarraylist%NNA0     = max(0,NNA)

    ! Allocate memory and associate pointers
    Isize=(/ALST_NA,rarraylist%nntable/)
    call storage_new('alst_create', 'Ktable', Isize,&
        ST_INT, rarraylist%h_Ktable, ST_NEWBLOCK_NOINIT)
    call storage_getbase_int2D(rarraylist%h_Ktable, rarraylist%p_Ktable)

    call storage_new('alst_create', 'Knext', ALST_FREE,&
        rarraylist%NNA, ST_INT, rarraylist%h_Knext, ST_NEWBLOCK_NOINIT)
    call storage_getbase_int(rarraylist%h_Knext, rarraylist%p_Knext)

    call storage_new('alst_create','Kprev', rarraylist%NNA, ST_INT,&
        rarraylist%h_Kprev, ST_NEWBLOCK_NOINIT)
    call storage_getbase_int(rarraylist%h_Kprev, rarraylist%p_Kprev)

    ! Allocate memory for Key
#ifdef T_STORAGE
    call storage_new('alst_create', 'Key', rarraylist%NNA, T_STORAGE,&
        rarraylist%h_Key, ST_NEWBLOCK_NOINIT)
    call storage_getbase(rarraylist%h_Key, rarraylist%p_Key)
#else
    allocate(rarraylist%p_Key(rarraylist%NNA))
#endif

#ifdef D
    ! Set size of auxiliary data
    rarraylist%isizeData = max(0,isizeData)

    ! Allocate memory for auxiliary data
    if (rarraylist%isizeData > 0) then
      Isize = (/rarraylist%isizeData, rarraylist%NNA/)

#ifdef D_STORAGE
      call storage_new('arraylist_create', 'Data', Isize, D_STORAGE,&
          rarraylist%h_Data, ST_NEWBLOCK_NOINIT)
      call storage_getbase(rarraylist%h_Data, rarraylist%p_Data)
#else
      allocate(rarraylist%p_Data(Isize(1),Isize(2)))
#endif
    end if
#endif

    ! Clear arraylist
    call alst_clear(rarraylist)

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine FEAT2_PP_TEMPLATE_TD(alst_createTbl,T,D)(rarraylist, itable)

!<description>
    ! This subroutine creates new table entries up to position itable.
!</description>

!<input>
    ! Number of last table
    integer, intent(in) :: itable
!</input>

!<inputoutput>
    ! The arraylist
    type(FEAT2_PP_TEMPLATE_TD(t_arraylist,T,D)), intent(inout) :: rarraylist
!</inputoutput>
!</subroutine>

    ! local variable
    integer :: i

    ! Check if table number is valid
    if (itable < 1) then
      call output_line('Invalid table number',&
          OU_CLASS_ERROR,OU_MODE_STD,'alst_createTable')
      call sys_halt()
    end if

    ! Resize tables if required
    if (itable > rarraylist%nntable) call alst_resizeTbl(&
        rarraylist, ceiling(itable*rarraylist%dfactor))

    ! Initialise structures
    do i = rarraylist%ntable+1, itable
      rarraylist%p_Ktable(ALST_HEAD,i) = ALST_NULL
      rarraylist%p_Ktable(ALST_TAIL,i) = ALST_NULL
      rarraylist%p_Ktable(ALST_NA,i)   = 0
    end do

    ! Set new table size
    rarraylist%ntable = max(rarraylist%ntable, itable)

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine FEAT2_PP_TEMPLATE_TD(alst_release,T,D)(rarraylist)

!<description>
    ! This subroutine releases an existing arraylist
!</description>

!<inputoutput>
    type(FEAT2_PP_TEMPLATE_TD(t_arraylist,T,D)), intent(inout) :: rarraylist
!</inputoutput>
!</subroutine>

    ! Release memory
    if (rarraylist%h_Ktable .ne. ST_NOHANDLE) call storage_free(rarraylist%h_Ktable)
    if (rarraylist%h_Knext  .ne. ST_NOHANDLE) call storage_free(rarraylist%h_Knext)
    if (rarraylist%h_Kprev  .ne. ST_NOHANDLE) call storage_free(rarraylist%h_Kprev)

    nullify(rarraylist%p_Ktable)
    nullify(rarraylist%p_Knext)
    nullify(rarraylist%p_Kprev)

#ifdef T_STORAGE
    if (rarraylist%h_Key    .ne. ST_NOHANDLE) call storage_free(rarraylist%h_Key)
    nullify(rarraylist%p_Key)
#else
    if (associated(rarraylist%p_Key))&
        deallocate(rarraylist%p_Key)
#endif

#ifdef D
    ! Reset size of auxiliary data
    rarraylist%isizeData = 0
#ifdef D_STORAGE
    if (rarraylist%h_Data   .ne. ST_NOHANDLE) call storage_free(rarraylist%h_Data)
    nullify(rarraylist%p_Data)
#else
    if (associated(rarraylist%p_Data))&
        deallocate(rarraylist%p_Data)
#endif
#endif

    ! Reset list
    rarraylist%ntable      = 0
    rarraylist%nntable     = 0
    rarraylist%nntable0    = 0
    rarraylist%NA          = 0
    rarraylist%NNA         = 0
    rarraylist%NNA0        = 0
    rarraylist%dfactor     = 1.5_DP
    rarraylist%nresize     = 0

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine FEAT2_PP_TEMPLATE_TD(alst_releaseTbl,T,D)(rarraylist, itable)

!<description>
    ! This subroutine releases a table from the arraylist
!</description>

!<input>
    ! Number of the table to be released
    integer, intent(in) :: itable
!</input>

!<inputoutput>
    ! The arraylist
    type(FEAT2_PP_TEMPLATE_TD(t_arraylist,T,D)), intent(inout) :: rarraylist
!</inputoutput>
!</subroutine>

    ! Check if table exists
    if (itable < 1 .or. itable > rarraylist%ntable) then
      call output_line('Invalid table number!',&
          OU_CLASS_ERROR,OU_MODE_STD,'alst_releaseTbl')
      call sys_halt()
    end if

    ! Decrease number of entries by the number of entries present
    ! in the table which is released
    rarraylist%NA = rarraylist%NA - rarraylist%p_Ktable(ALST_NA,itable)

    ! Reset table
    rarraylist%p_Ktable(ALST_HEAD,itable) = ALST_NULL
    rarraylist%p_Ktable(ALST_HEAD,itable) = ALST_NULL
    rarraylist%p_Ktable(ALST_NA,itable)   = 0

    ! Decrease number of tables if the last table has been deleted
    if (itable .eq. rarraylist%ntable)&
        rarraylist%ntable = rarraylist%ntable-1

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine FEAT2_PP_TEMPLATE_TD(alst_resize,T,D)(rarraylist, NNA)

!<description>
    ! This subroutine reallocates memory for an existing list
!</description>

!<input>
    ! New number of total items that can be stored in the list
    integer, intent(in) :: NNA
!</input>

!<inputoutput>
    ! The arraylist
    type(FEAT2_PP_TEMPLATE_TD(t_arraylist,T,D)), intent(inout) :: rarraylist
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: NNAOld

#ifndef T_STORAGE
    FEAT2_PP_TTYPE(T_TYPE), dimension(:), pointer :: p_Key
#endif

#ifdef D
#ifndef D_STORAGE
    FEAT2_PP_DTYPE(D_TYPE), dimension(:,:), pointer :: p_Data
#endif
#endif

    ! Save old size
    NNAOld = rarraylist%NNA

    ! Set new size and increase counter
    rarraylist%NNA = max(0,NNA)
    rarraylist%nresize = rarraylist%nresize+1

    call storage_realloc('alst_resize', rarraylist%NNA+1,&
        rarraylist%h_Knext, ST_NEWBLOCK_NOINIT, .true.)
    call storage_getbase_int(rarraylist%h_Knext, rarraylist%p_Knext)

    call storage_realloc('alst_resize', rarraylist%NNA,&
        rarraylist%h_Kprev, ST_NEWBLOCK_NOINIT, .true.)
    call storage_getbase_int(rarraylist%h_Kprev, rarraylist%p_Kprev)

    ! Reallocate Key
#ifdef T_STORAGE
    call storage_realloc('arraylist_resize', rarraylist%NNA,&
        rarraylist%h_Key, ST_NEWBLOCK_NOINIT, .true.)
    call storage_getbase(rarraylist%h_Key, rarraylist%p_Key)
#else
    allocate(p_Key(NNAOld))
    p_Key = rarraylist%p_Key
    deallocate(rarraylist%p_Key)
    allocate(rarraylist%p_Key(rarraylist%NNA))
    rarraylist%p_Key(1:NNAOld) = p_Key
    deallocate(p_Key)
#endif

#ifdef D
    ! Reallocate auxiliary data
    if (rarraylist%isizeData > 0) then
#ifdef D_STORAGE
      call storage_realloc('arraylist_resize', rarraylist%NNA,&
          rarraylist%h_Data, ST_NEWBLOCK_NOINIT, .true.)
      call storage_getbase(rarraylist%h_Data, rarraylist%p_Data)
#else
      allocate(p_Data(size(rarraylist%p_Data,1),NNAOld))
      p_Data = rarraylist%p_Data
      deallocate(rarraylist%p_Data)
      allocate(rarraylist%p_Data(size(p_Data,1),rarraylist%NNA))
      rarraylist%p_Data(:,1:NNAOld) = p_Data
      deallocate(p_Data)
#endif
    end if
#endif

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine FEAT2_PP_TEMPLATE_TD(alst_resizeTbl,T,D)(rarraylist, nntable)

!<description>
    ! This subroutine reallocates memory for the lookup table
!</description>

!<input>
    ! New number of tables
    integer, intent(in) :: nntable
!</input>

!<inputoutput>
    ! The arraylist
    type(FEAT2_PP_TEMPLATE_TD(t_arraylist,T,D)), intent(inout) :: rarraylist
!</inputoutput>
!</subroutine>

    ! Set new size
    rarraylist%nntable = nntable
    rarraylist%nresize = rarraylist%nresize+1

    call storage_realloc('alst_resizeTbl', rarraylist%nntable,&
        rarraylist%h_Ktable, ST_NEWBLOCK_NOINIT, .true.)
    call storage_getbase_int2D(rarraylist%h_Ktable, rarraylist%p_Ktable)

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine FEAT2_PP_TEMPLATE_TD(alst_clear,T,D)(rarraylist)

!<description>
    ! This subroutine clears the content of the arraylist
!</description>

!<inputoutput>
    ! The arraylist
    type(FEAT2_PP_TEMPLATE_TD(t_arraylist,T,D)), intent(inout) :: rarraylist
!</inputoutput>
!</subroutine>

    ! local variable
    integer :: itable

    do itable = 1, rarraylist%nntable
      call alst_clearTbl(rarraylist, itable)
    end do

    ! Initialise arraylist structure
    rarraylist%ntable   = 0
    rarraylist%NA       = 0
    rarraylist%nresize  = 0

    rarraylist%p_Knext(ALST_FREE) = 1

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine FEAT2_PP_TEMPLATE_TD(alst_clearTbl,T,D)(rarraylist, itable)

!<description>
    ! This subroutine clears the content of the table of the arraylist
!</description>

!<input>
    ! Number of the table to be cleared
    integer, intent(in) :: itable
!</input>

!<inputoutput>
    ! The arraylist
    type(FEAT2_PP_TEMPLATE_TD(t_arraylist,T,D)), intent(inout) :: rarraylist
!</inputoutput>
!</subroutine>

    ! Decrease number of entries by the number of entries present
    ! in the table which is released
    rarraylist%NA = rarraylist%NA - rarraylist%p_Ktable(ALST_NA,itable)

    ! Reset table
    rarraylist%p_Ktable(ALST_HEAD,itable) = ALST_NULL
    rarraylist%p_Ktable(ALST_HEAD,itable) = ALST_NULL
    rarraylist%p_Ktable(ALST_NA,itable)   = 0

  end subroutine

  !************************************************************************

!<subroutine>

#ifdef D
  subroutine FEAT2_PP_TEMPLATE_TD(alst_cpy1,T,D)(h_TableSrc, h_KeySrc, rarraylist, h_DataSrc)
#else
  subroutine FEAT2_PP_TEMPLATE_TD(alst_cpy1,T,D)(h_TableSrc, h_KeySrc, rarraylist)
#endif

!<description>
    ! This subroutine copies the content of the given handle(s) to the arraylist
!</description>

!<input>
    ! Handle to the table
    integer, intent(in) :: h_TableSrc

    ! Handle to the key value
    integer, intent(in) :: h_KeySrc

#ifdef D
    ! OPTIONAL: Handle to the data values
    integer, intent(in), optional :: h_DataSrc
#endif
!</input>

!<inputoutput>
    ! The arraylist
    type(FEAT2_PP_TEMPLATE_TD(t_arraylist,T,D)), intent(inout) :: rarraylist
!</inputoutput>
!</subroutine>

    ! local variables
    integer, dimension(:), pointer :: p_TableSrc

#ifdef T_STORAGE
    T_TYPE, dimension(:), pointer :: p_KeySrc

#ifdef D
#ifdef D_STORAGE
    D_TYPE, dimension(:,:), pointer :: p_DataSrc

    ! Get pointers and copy the key,value pairs to the arraylist
    call storage_getbase_int(h_TableSrc, p_TableSrc)
    call storage_getbase(h_KeySrc, p_KeySrc)
    call storage_getbase(h_DataSrc, p_DataSrc)
    call alst_copy(p_TableSrc, p_KeySrc, rarraylist, p_DataSrc)
#else
    call output_line('Arraylist does not support storage handles!',&
        OU_CLASS_ERROR,OU_MODE_STD,'alst_cpy1')
    call sys_halt()
#endif

#else
    ! Get pointer and copy the keys to the arraylist
    call storage_getbase_int(h_TableSrc, p_TableSrc)
    call storage_getbase(h_KeySrc, p_KeySrc)
    call alst_copy(p_TableSrc, p_KeySrc, rarraylist)
#endif

#else
    call output_line('Arraylist does not support storage handles!',&
        OU_CLASS_ERROR,OU_MODE_STD,'list_cpy1')
    call sys_halt()
#endif

  end subroutine

  !************************************************************************

!<subroutine>

#ifdef D
  subroutine FEAT2_PP_TEMPLATE_TD(alst_cpy2,T,D)(TableSrc, KeySrc, rarraylist, DataSrc)
#else
  subroutine FEAT2_PP_TEMPLATE_TD(alst_cpy2,T,D)(TableSrc, KeySrc, rarraylist)
#endif

!<description>
    ! This subroutine copies the content of the given array(s) to the arraylist
!</description>

!<input>
    ! Array with table structure
    integer, dimension(:), intent(in) :: TableSrc

    ! Array with key values
    FEAT2_PP_TTYPE(T_TYPE), dimension(:), intent(in) :: KeySrc

#ifdef D
    ! OPTIONAL: Array with data values
    FEAT2_PP_DTYPE(D_TYPE), dimension(:,:), intent(in), optional :: DataSrc
#endif
!</input>

!<inputoutput>
    ! The linked list
    type(FEAT2_PP_TEMPLATE_TD(t_arraylist,T,D)), intent(inout) :: rarraylist
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: ipos,itable   

#ifdef D
    if (present(DataSrc)) then
      do itable = 1, size(TableSrc)-1
        do ipos = TableSrc(itable), TableSrc(itable+1)-1
          call alst_push_back(rarraylist, itable, KeySrc(ipos), DataSrc(:,ipos))
        end do
      end do
    else
      do itable = 1, size(TableSrc)-1
        do ipos = TableSrc(itable), TableSrc(itable+1)-1
          call alst_push_back(rarraylist, itable, KeySrc(ipos))
        end do
      end do
    end if
#else
    do itable = 1, size(TableSrc)-1
      do ipos = TableSrc(itable), TableSrc(itable+1)-1
        call alst_push_back(rarraylist, itable, KeySrc(ipos))
      end do
    end do
#endif

  end subroutine

  !************************************************************************

!<subroutine>

#ifdef D
  subroutine FEAT2_PP_TEMPLATE_TD(alst_cpy3,T,D)(rarraylist, h_TableDest, h_KeyDest, h_DataDest)
#else
  subroutine FEAT2_PP_TEMPLATE_TD(alst_cpy3,T,D)(rarraylist, h_TableDest, h_KeyDest)
#endif

!<description>
    ! This subroutine copies the content of the given arraylist to the handles
!</description>

!<input>
    ! The arraylist
    type(FEAT2_PP_TEMPLATE_TD(t_arraylist,T,D)), intent(in) :: rarraylist
!</input>

!<inputoutput>
    ! Handle to the table structure
    integer, intent(inout) :: h_TableDest

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
    integer, dimension(:), pointer :: p_TableDest

#ifdef T_STORAGE
    T_TYPE, dimension(:), pointer :: p_KeyDest

#ifdef D
    integer, dimension(2) :: Isize2
#ifdef D_STORAGE
    D_TYPE, dimension(:,:), pointer :: p_DataDest

    ! Check table handle
    if (h_TableDest .eq. ST_NOHANDLE) then
      call storage_new('alst_cpy3', 'Table',&
          rarraylist%ntable+1, ST_INT, h_TableDest, ST_NEWBLOCK_NOINIT)
    else
      call storage_getsize(h_TableDest, isize)
      if (isize < rarraylist%ntable+1) then
        call storage_realloc('alst_cpy3',&
            rarraylist%ntable+1, h_TableDest, ST_NEWBLOCK_NOINIT, .false.)
      end if
    end if
    call storage_getbase_int(h_TableDest, p_TableDest)

    ! Check key handle
    if (h_KeyDest .eq. ST_NOHANDLE) then
      call storage_new('alst_cpy3', 'Key',&
          rarraylist%NA, T_STORAGE, h_KeyDest, ST_NEWBLOCK_NOINIT)
    else
      call storage_getsize(h_KeyDest, isize)
      if (isize < rarraylist%NA) then
        call storage_realloc('alst_cpy3',&
            rarraylist%NA, h_KeyDest, ST_NEWBLOCK_NOINIT, .false.)
      end if
    end if
    call storage_getbase(h_KeyDest, p_KeyDest)

    ! Do we have to copy auxiliary data?
    if (present(h_DataDest)) then

      if (rarraylist%isizeData <= 0) then
        if (h_DataDest .ne. ST_NOHANDLE) call storage_free(h_DataDest)
        call output_line('Arraylist does not provide auxiliary data!',&
            OU_CLASS_WARNING,OU_MODE_STD,'alst_cpy3')

        ! Copy key values from list
        call alst_copy(rarraylist, p_TableDest, p_KeyDest)

        ! That`s it
        return
      end if

      ! Check data handle
      if (h_DataDest .eq. ST_NOHANDLE) then
        call storage_new('alst_cpy3', 'Data',&
            (/rarraylist%isizeData,rarraylist%NA/),&
            D_STORAGE, h_DataDest, ST_NEWBLOCK_NOINIT)
      else
        call storage_getsize(h_DataDest, Isize2)

        if (Isize2(1) .ne. rarraylist%isizeData) then
          call output_line('Size of data array is not compatible!',&
              OU_CLASS_ERROR,OU_MODE_STD,'alst_cpy3')
        end if

        if (Isize2(2) < rarraylist%NA) then
          call storage_realloc('alst_cpy3',&
              rarraylist%NA, h_DataDest, ST_NEWBLOCK_NOINIT, .false.)
        end if
      end if
      call storage_getbase(h_DataDest, p_DataDest)

      ! Copy table, key and data values from arraylist
      call alst_copy(rarraylist, p_TableDest, p_KeyDest, p_DataDest)

    else

      ! Copy table, and key values from arraylist
      call alst_copy(rarraylist, p_TableDest, p_KeyDest)

    end if
#else
    call output_line('Arraylist does not support storage handles!',&
        OU_CLASS_ERROR,OU_MODE_STD,'alst_cpy3')
    call sys_halt()
#endif

#else

    ! Check table handle
    if (h_TableDest .eq. ST_NOHANDLE) then
      call storage_new('alst_cpy3', 'Table',&
          rarraylist%ntable+1, ST_INT, h_TableDest, ST_NEWBLOCK_NOINIT)
    else
      call storage_getsize(h_TableDest, isize)
      if (isize < rarraylist%ntable+1) then
        call storage_realloc('alst_cpy3',&
            rarraylist%ntable+1, h_TableDest, ST_NEWBLOCK_NOINIT, .false.)
      end if
    end if
    call storage_getbase_int(h_TableDest, p_TableDest)

    ! Check key handle
    if (h_KeyDest .eq. ST_NOHANDLE) then
      call storage_new('alst_cpy3', 'Key',&
          rarraylist%NA, T_STORAGE, h_KeyDest, ST_NEWBLOCK_NOINIT)
    else
      call storage_getsize(h_KeyDest, isize)
      if (isize < rarraylist%NA) then
        call storage_realloc('alst_cpy3',&
            rarraylist%NA, h_KeyDest, ST_NEWBLOCK_NOINIT, .false.)
      end if
    end if
    call storage_getbase(h_KeyDest, p_KeyDest)

    ! Copy table, and key values from arraylist
    call alst_copy(rarraylist, p_TableDest, p_KeyDest)

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
  subroutine FEAT2_PP_TEMPLATE_TD(alst_cpy4,T,D)(rarraylist, TableDest, KeyDest, DataDest)
#else
  subroutine FEAT2_PP_TEMPLATE_TD(alst_cpy4,T,D)(rarraylist, TableDest, KeyDest)
#endif

!<description>
    ! This subroutine copies the content of the given arraylist to the arrays
!</description>

!<input>
    ! The arraylist
    type(FEAT2_PP_TEMPLATE_TD(t_arraylist,T,D)), intent(in) :: rarraylist
!</input>

!<inputoutput>
    ! Array with table structure
    integer, dimension(:), intent(inout) :: TableDest

    ! Array with key values
    FEAT2_PP_TTYPE(T_TYPE), dimension(:), intent(inout) :: KeyDest

#ifdef D
    ! OPTIONAL: Array with data values
    FEAT2_PP_DTYPE(D_TYPE), dimension(:,:), intent(inout), optional :: DataDest
#endif
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: icount,itable,ntable,ipos

    ! Check if table array is valid
    ntable = size(TableDest)-1
    if (ntable+1 < rarraylist%ntable) then
      call output_line('Invalid dimension of table array!',&
          OU_CLASS_ERROR,OU_MODE_STD,'alst_cpy4')
      call sys_halt()
    end if

#ifdef D

    if (present(DataDest)) then

      icount = 1
      do itable = 1, ntable
        TableDest(itable) = icount

        ipos = rarraylist%p_Ktable(ALST_HEAD,itable)
        kd_copy: do while (ipos .ne. ALST_NULL)
          KeyDest(icount)    = rarraylist%p_Key(ipos)
          DataDest(:,icount) = rarraylist%p_Data(:,ipos)
          icount = icount+1
          if (ipos .eq. rarraylist%p_Ktable(ALST_TAIL,itable)) exit kd_copy
          ipos = rarraylist%p_Knext(ipos)
        end do kd_copy
      end do
      TableDest(ntable+1) = icount

    else

      icount = 1
      do itable = 1, ntable
        TableDest(itable) = icount

        ipos = rarraylist%p_Ktable(ALST_HEAD,itable)
        k_copy: do while (ipos .ne. ALST_NULL)
          KeyDest(icount)    = rarraylist%p_Key(ipos)
          icount = icount+1
          if (ipos .eq. rarraylist%p_Ktable(ALST_TAIL,itable)) exit k_copy
          ipos = rarraylist%p_Knext(ipos)
        end do k_copy
      end do
      TableDest(ntable+1) = icount

    end if

#else

    icount = 1
      do itable = 1, ntable
        TableDest(itable) = icount

        ipos = rarraylist%p_Ktable(ALST_HEAD,itable)
        copy: do while (ipos .ne. ALST_NULL)
          KeyDest(icount) = rarraylist%p_Key(ipos)
          icount = icount+1
          if (ipos .eq. rarraylist%p_Ktable(ALST_TAIL,itable)) exit copy
          ipos = rarraylist%p_Knext(ipos)
        end do copy
      end do
      TableDest(ntable+1) = icount

#endif

  end subroutine

  !************************************************************************

!<subroutine>

#ifdef D
  subroutine FEAT2_PP_TEMPLATE_TD(alst_cpy1Tbl,T,D)(h_KeySrc, itable, rarraylist, h_DataSrc)
#else
  subroutine FEAT2_PP_TEMPLATE_TD(alst_cpy1Tbl,T,D)(h_KeySrc, itable, rarraylist)
#endif

!<description>
    ! This subroutine copies the content of the given handle(s) to the
    ! list associated with the given table number
!</description>

!<input>
    ! Handle to the key value
    integer, intent(in) :: h_KeySrc

    ! Number of the table
    integer, intent(in) :: itable

#ifdef D
    ! OPTIONAL: Handle to the data values
    integer, intent(in), optional :: h_DataSrc
#endif
!</input>

!<inputoutput>
    ! The arraylist
    type(FEAT2_PP_TEMPLATE_TD(t_arraylist,T,D)), intent(inout) :: rarraylist
!</inputoutput>
!</subroutine>

#ifdef T_STORAGE
    ! local variables
    T_TYPE, dimension(:), pointer :: p_KeySrc

#ifdef D
#ifdef D_STORAGE
    D_TYPE, dimension(:,:), pointer :: p_DataSrc

    ! Get pointers and copy the key and data values to the arraylist
    call storage_getbase(h_KeySrc, p_KeySrc)

    if (present(h_DataSrc)) then
      call storage_getbase(h_DataSrc, p_DataSrc)
      call alst_copyTbl(p_KeySrc, itable, rarraylist, p_DataSrc)
    else
      call alst_copyTbl(p_KeySrc, itable, rarraylist)
    end if
#else
    call output_line('List does not support storage handles!',&
        OU_CLASS_ERROR,OU_MODE_STD,'alst_cpy1Tbl')
    call sys_halt()
#endif

#else
    ! Get pointer and copy the keys to the arraylist
    call storage_getbase(h_KeySrc, p_KeySrc)
    call alst_copyTbl(p_KeySrc, itable, rarraylist)
#endif

#else
    call output_line('List does not support storage handles!',&
        OU_CLASS_ERROR,OU_MODE_STD,'alst_cpy1Tbl')
    call sys_halt()
#endif

  end subroutine

  !************************************************************************

!<subroutine>

#ifdef D
  subroutine FEAT2_PP_TEMPLATE_TD(alst_cpy2Tbl,T,D)(KeySrc, itable, rarraylist, DataSrc)
#else
  subroutine FEAT2_PP_TEMPLATE_TD(alst_cpy2Tbl,T,D)(KeySrc, itable, rarraylist)
#endif

!<description>
    ! This subroutine copies the content of the given array(s) to the
    ! list associated with the given table number
!</description>

!<input>
    ! Array with key values
    FEAT2_PP_TTYPE(T_TYPE), dimension(:), intent(in) :: KeySrc

    ! Number of the table
    integer, intent(in) :: itable

#ifdef D
    ! OPTIONAL: Array with data values
    FEAT2_PP_DTYPE(D_TYPE), dimension(:,:), intent(in), optional :: DataSrc
#endif
!</input>

!<inputoutput>
    ! The linked list
    type(FEAT2_PP_TEMPLATE_TD(t_arraylist,T,D)), intent(inout) :: rarraylist
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: ipos

#ifdef D
    if (present(DataSrc)) then
      do ipos = 1, size(KeySrc)
        call alst_push_back(rarraylist, itable, KeySrc(ipos), DataSrc(:,ipos))      
      end do
    else
      do ipos = 1, size(KeySrc)
        call alst_push_back(rarraylist, itable, KeySrc(ipos))
      end do
    end if
#else
    do ipos = 1, size(KeySrc)
      call alst_push_back(rarraylist, itable, KeySrc(ipos))
    end do
#endif

  end subroutine

  !************************************************************************

!<subroutine>

#ifdef D
  subroutine FEAT2_PP_TEMPLATE_TD(alst_cpy3Tbl,T,D)(rarraylist, itable, h_KeyDest, h_DataDest, ncount)
#else
  subroutine FEAT2_PP_TEMPLATE_TD(alst_cpy3Tbl,T,D)(rarraylist, itable, h_KeyDest, ncount)
#endif

!<description>
    ! This subroutine copies the content of the list associated with
    ! the given table number to the given handle(s)
!</description>

!<input>
    ! The arraylist
    type(FEAT2_PP_TEMPLATE_TD(t_arraylist,T,D)), intent(in) :: rarraylist

    ! Number of the table
    integer, intent(in) :: itable
!</input>

!<inputoutput>
    ! Handle to the key values
    integer, intent(inout) :: h_KeyDest

#ifdef D
    ! OPTIONAL: Handle to the data values
    integer, intent(inout), optional :: h_DataDest
#endif
!</inputoutput>

!<output>
    ! OPTIONAL: number of entries
    integer, intent(out), optional :: ncount
!</output>
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
        call storage_new('alst_cpy3', 'Key',&
            rarraylist%NA, T_STORAGE, h_KeyDest, ST_NEWBLOCK_NOINIT)
      else
        call storage_getsize(h_KeyDest, isize)
        if (isize < rarraylist%NA) then
          call storage_realloc('alst_cpy3',&
              rarraylist%NA, h_KeyDest, ST_NEWBLOCK_NOINIT, .false.)
        end if
      end if
      call storage_getbase(h_KeyDest, p_KeyDest)

      if (rarraylist%isizeData <= 0) then
        call output_line('Arraylist does not provide auxiliary data!',&
            OU_CLASS_ERROR,OU_MODE_STD,'alst_cpy3Tbl')

        ! Copy key values from arraylist
        call alst_copyTbl(rarraylist, itable, p_KeyDest, ncount=ncount)

        ! That`s it
        return
      end if

      ! Check data handle
      if (h_DataDest .eq. ST_NOHANDLE) then
        call storage_new('alst_cpy3', 'Data',&
            (/rarraylist%isizeData,rarraylist%NA/),&
            D_STORAGE, h_DataDest, ST_NEWBLOCK_NOINIT)
      else
        call storage_getsize(h_DataDest, Isize2)

        if (Isize2(1) .ne. rarraylist%isizeData) then
          call output_line('Size of data array is not compatible!',&
              OU_CLASS_ERROR,OU_MODE_STD,'alst_cpy3')
        end if

        if (Isize2(2) < rarraylist%NA) then
          call storage_realloc('alst_cpy3',&
              rarraylist%NA, h_DataDest, ST_NEWBLOCK_NOINIT, .false.)
        end if
      end if
      call storage_getbase(h_DataDest, p_DataDest)

      ! Copy key and data values from arraylist
      call alst_copyTbl(rarraylist, itable, p_KeyDest, p_DataDest, ncount)

    else

      ! Check key handle
      if (h_KeyDest .eq. ST_NOHANDLE) then
        call storage_new('alst_cpy3', 'Key',&
            rarraylist%NA, T_STORAGE, h_KeyDest, ST_NEWBLOCK_NOINIT)
      else
        call storage_getsize(h_KeyDest, isize)
        if (isize < rarraylist%NA) then
          call storage_realloc('alst_cpy3',&
              rarraylist%NA, h_KeyDest, ST_NEWBLOCK_NOINIT, .false.)
        end if
      end if
      call storage_getbase(h_KeyDest, p_KeyDest)

      ! Copy key values from arraylist
      call alst_copyTbl(rarraylist, itable, p_KeyDest, ncount=ncount)

    end if
#else
    call output_line('Arraylist does not support storage handles!',&
        OU_CLASS_ERROR,OU_MODE_STD,'alst_cpy3Tbl')
    call sys_halt()
#endif

#else

     ! Check key handle
      if (h_KeyDest .eq. ST_NOHANDLE) then
        call storage_new('alst_cpy3', 'Key',&
            rarraylist%NA, T_STORAGE, h_KeyDest, ST_NEWBLOCK_NOINIT)
      else
        call storage_getsize(h_KeyDest, isize)
        if (isize < rarraylist%NA) then
          call storage_realloc('alst_cpy3',&
              rarraylist%NA, h_KeyDest, ST_NEWBLOCK_NOINIT, .false.)
        end if
      end if
      call storage_getbase(h_KeyDest, p_KeyDest)

    ! Copy key values from arraylist
    call alst_copyTbl(rarraylist, itable, p_KeyDest, ncount)

#endif

#else
    call output_line('Arraylist does not support storage handles!',&
        OU_CLASS_ERROR,OU_MODE_STD,'alst_cpy3Tbl')
    call sys_halt()
#endif

  end subroutine

  !************************************************************************

!<subroutine>

#ifdef D
  subroutine FEAT2_PP_TEMPLATE_TD(alst_cpy4Tbl,T,D)(rarraylist, itable, Key, data, ncount)
#else
  subroutine FEAT2_PP_TEMPLATE_TD(alst_cpy4Tbl,T,D)(rarraylist, itable, Key, ncount)
#endif

!<description>
    ! This subroutine copies the content of the list associated with
    ! the given table number to the given array(s)
!</description>

!<input>
    ! The arraylist
    type(FEAT2_PP_TEMPLATE_TD(t_arraylist,T,D)), intent(in) :: rarraylist

    ! Number of the table
    integer, intent(in) :: itable
!</input>

!<inputoutput>
    ! Array with key values
    FEAT2_PP_TTYPE(T_TYPE), dimension(:), intent(inout) :: Key

#ifdef D
    ! OPTIONAL: Array with data values
    FEAT2_PP_DTYPE(D_TYPE), dimension(:,:), intent(inout), optional :: data
#endif
!</inputoutput>

!<output>
    ! OPTIONAL: number of entries
    integer, intent(out), optional :: ncount
!</output>
!</subroutine>

    ! local variables
    integer :: ipos,icount

    ! Check if table is valid
    if (itable < 1 .or. itable > rarraylist%ntable) then
      call output_line('Invalid table number',&
          OU_CLASS_ERROR,OU_MODE_STD,'alst_cpy4Tbl')
      call sys_halt()
    end if

    icount = 0
    ipos   = rarraylist%p_Ktable(ALST_HEAD,itable)

#ifdef D

    if (present(data)) then

      kd_copy: do
        icount = icount+1
        Key(icount)    = rarraylist%p_Key(ipos)
        data(:,icount) = rarraylist%p_Data(:,ipos)
        if (ipos .eq. rarraylist%p_Ktable(ALST_TAIL,itable)) exit kd_copy
        ipos = rarraylist%p_Knext(ipos)
      end do kd_copy

    else

      k_copy: do
        icount = icount+1
        Key(icount) = rarraylist%p_Key(ipos)
        if (ipos .eq. rarraylist%p_Ktable(ALST_TAIL,itable)) exit k_copy
        ipos = rarraylist%p_Knext(ipos)
      end do k_copy

    end if

#else

    copy: do
      icount = icount+1
      Key(icount) = rarraylist%p_Key(ipos)
      if (ipos .eq. rarraylist%p_Ktable(ALST_TAIL,itable)) exit copy
      ipos = rarraylist%p_Knext(ipos)
    end do copy

#endif

    if (present(ncount)) ncount = icount

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine FEAT2_PP_TEMPLATE_TD(alst_swap,T,D)(rarraylist1, rarraylist2)

!<description>
    ! This subroutine swaps content of two arraylists
!</description>

!<inputoutput>
    ! First arraylist
    type(FEAT2_PP_TEMPLATE_TD(t_arraylist,T,D)), intent(inout) :: rarraylist1

    ! Second arraylist
    type(FEAT2_PP_TEMPLATE_TD(t_arraylist,T,D)), intent(inout) :: rarraylist2
!</inputoutput>
!</subroutine>

    ! local variables
    type(FEAT2_PP_TEMPLATE_TD(t_arraylist,T,D)) :: rarraylist

    ! Swap lists
    rarraylist  = rarraylist1
    rarraylist1 = rarraylist2
    rarraylist2 = rarraylist

    ! Reassociated pointers of first list
    call storage_getbase_int(rarraylist1%h_Knext, rarraylist1%p_Knext)
    call storage_getbase_int(rarraylist1%h_Kprev, rarraylist1%p_Kprev)

    ! Reassociated pointers of second list
    call storage_getbase_int(rarraylist2%h_Knext, rarraylist2%p_Knext)
    call storage_getbase_int(rarraylist2%h_Kprev, rarraylist2%p_Kprev)

#ifdef T_STORAGE
    ! Reassociated pointers to key values
    call storage_getbase(rarraylist1%h_Key, rarraylist1%p_Key)
    call storage_getbase(rarraylist2%h_Key, rarraylist2%p_Key)
#endif

#if defined D && defined D_STORAGE
    ! Reassociated pointers to data values
    if (rarraylist1%h_Data .ne. ST_NOHANDLE) &
        call storage_getbase(rarraylist1%h_Data, rarraylist1%p_Data)
    if (rarraylist2%h_Data .ne. ST_NOHANDLE) &
        call storage_getbase(rarraylist2%h_Data, rarraylist2%p_Data)
#endif

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine FEAT2_PP_TEMPLATE_TD(alst_swapTbl,T,D)(rarraylist, itable, jtable)

!<description>
    ! This subroutine swaps two tables and the associated list in the
    ! given arraylist
!</description>

!<input>
    ! Numbers of the tables to be swapped
    integer, intent(in) :: itable,jtable
!</input>

!<inputoutput>
    ! The arraylist
    type(FEAT2_PP_TEMPLATE_TD(t_arraylist,T,D)), intent(inout) :: rarraylist
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: ihead,itail,ina

    ! Check if table exists
    if (itable < 1 .or. itable > rarraylist%ntable) then
      call output_line('Invalid table number!',&
          OU_CLASS_ERROR,OU_MODE_STD,'alst_swapTbl')
      call sys_halt()
    end if

    ! Swap
    ihead = rarraylist%p_Ktable(ALST_HEAD,itable)
    itail = rarraylist%p_Ktable(ALST_TAIL,itable)
    ina   = rarraylist%p_Ktable(ALST_NA,  itable)

    rarraylist%p_Ktable(ALST_HEAD,itable) = rarraylist%p_Ktable(ALST_HEAD,jtable)
    rarraylist%p_Ktable(ALST_TAIL,itable) = rarraylist%p_Ktable(ALST_TAIL,jtable)
    rarraylist%p_Ktable(ALST_NA,  itable) = rarraylist%p_Ktable(ALST_NA,  jtable)

    rarraylist%p_Ktable(ALST_HEAD,jtable) = ihead
    rarraylist%p_Ktable(ALST_TAIL,jtable) = itail
    rarraylist%p_Ktable(ALST_NA,  jtable) = ina

  end subroutine

  !************************************************************************

!<function>

  function FEAT2_PP_TEMPLATE_TD(alst_begin,T,D)(rarraylist, itable) result(riterator)

!<description>
    ! This function returns an iterator referring to the first list
    ! item in the given table
!</description>

!<input>
    ! The arraylist
    type(FEAT2_PP_TEMPLATE_TD(t_arraylist,T,D)), intent(in), target :: rarraylist

    ! Number of the table
    integer, intent(in) :: itable
!</input>

!<result>
    ! The iterator
    type(FEAT2_PP_TEMPLATE_TD(it_arraylist,T,D)) :: riterator
!</result>
!</function>

    if (itable < 1 .or. itable > rarraylist%ntable) then
      call output_line('Invalid table number!',&
          OU_CLASS_WARNING,OU_MODE_STD,'alst_begin')
      riterator%itable = 0
      riterator%ipos   = 0
      riterator%iSpec  = 0_I32
      nullify(riterator%p_rarraylist)
      return
    end if

    ! Attach arraylist to iterator
    riterator%p_rarraylist => rarraylist

    ! Initialise iterator
    riterator%itable = itable
    riterator%ipos   = rarraylist%p_Ktable(ALST_HEAD,itable)
    riterator%iSpec  = 0_I32

  end function

  !************************************************************************

!<function>

  function FEAT2_PP_TEMPLATE_TD(alst_rbegin,T,D)(rarraylist, itable) result(riterator)

!<description>
    ! This function returns a reverse iterator positioned to the last
    ! element of the list in the given table
!</description>

!<input>
    ! The arraylist
    type(FEAT2_PP_TEMPLATE_TD(t_arraylist,T,D)), intent(in), target :: rarraylist

    ! Number of the table
    integer, intent(in) :: itable
!</input>

!<result>
    ! The iterator
    type(FEAT2_PP_TEMPLATE_TD(it_arraylist,T,D)) :: riterator
!</result>
!</function>

    if (itable < 1 .or. itable > rarraylist%ntable) then
      call output_line('Invalid table number!',&
          OU_CLASS_WARNING,OU_MODE_STD,'alst_rbegin')
      riterator%itable = 0
      riterator%ipos   = 0
      riterator%iSpec  = 0_I32
      nullify(riterator%p_rarraylist)
      return
    end if

    ! Attach arraylist to iterator
    riterator%p_rarraylist => rarraylist

    ! Initialise iterator
    riterator%itable = itable
    riterator%ipos   = rarraylist%p_Ktable(ALST_TAIL,itable)
    riterator%iSpec  = ALST_LSPEC_REVERSE

  end function

  !************************************************************************

!<function>

  function FEAT2_PP_TEMPLATE_TD(alst_end,T,D)(rarraylist, itable) result(riterator)

!<description>
    ! This function returns an iterator referring to the past-the-end
    ! element of the list in the given table
!</description>

!<input>
    ! The arraylist
    type(FEAT2_PP_TEMPLATE_TD(t_arraylist,T,D)), intent(in), target :: rarraylist

    ! Number of the table
    integer, intent(in) :: itable
!</input>

!<result>
    ! The iterator
    type(FEAT2_PP_TEMPLATE_TD(it_arraylist,T,D)) :: riterator
!</result>
!</function>

    if (itable < 1 .or. itable > rarraylist%ntable) then
      call output_line('Invalid table number!',&
          OU_CLASS_WARNING,OU_MODE_STD,'alst_end')
      riterator%itable = 0
      riterator%ipos   = 0
      riterator%iSpec  = 0_I32
      nullify(riterator%p_rarraylist)
      return
    end if

    ! Attach arraylist to iterator
    riterator%p_rarraylist => rarraylist

    ! Initialise iterator
    riterator%itable = itable
    riterator%ipos   = ALST_NULL
    riterator%iSpec  = 0_I32

  end function

  !************************************************************************

!<function>

  function FEAT2_PP_TEMPLATE_TD(alst_rend,T,D)(rarraylist, itable) result(riterator)

!<description>
    ! This function returns an iterator referring to the element right
    ! before the first element of the list in the given table
!</description>

!<input>
    ! The arraylist
    type(FEAT2_PP_TEMPLATE_TD(t_arraylist,T,D)), intent(in), target :: rarraylist

    ! Number of the table
    integer, intent(in) :: itable
!</input>

!<result>
    ! The iterator
    type(FEAT2_PP_TEMPLATE_TD(it_arraylist,T,D)) :: riterator
    !</result>
!</function>

    if (itable < 1 .or. itable > rarraylist%ntable) then
      call output_line('Invalid table number!',&
          OU_CLASS_WARNING,OU_MODE_STD,'alst_rend')
      riterator%itable = 0
      riterator%ipos   = 0
      riterator%iSpec  = 0_I32
      nullify(riterator%p_rarraylist)
      return
    end if

    ! Attach arraylist to iterator
    riterator%p_rarraylist => rarraylist

    ! Initialise iterator
    riterator%itable = itable
    riterator%ipos   = ALST_NULL
    riterator%iSpec  = ALST_LSPEC_REVERSE

  end function

  !************************************************************************

!<subroutine>

  pure subroutine FEAT2_PP_TEMPLATE_TD(alst_next,T,D)(riterator)

!<description>
    ! This subroutine increments the list iterator by one.
!</description>

!<inputoutput>
    ! The iterator
    type(FEAT2_PP_TEMPLATE_TD(it_arraylist,T,D)), intent(inout) :: riterator
!</inputoutput>
!</subroutine>

    if (riterator%ipos .ne. ALST_NULL) then
      if (iand(riterator%iSpec, ALST_LSPEC_REVERSE).eq.0) then
        riterator%ipos = riterator%p_rarraylist%p_Knext(riterator%ipos)
      else
        riterator%ipos = riterator%p_rarraylist%p_Kprev(riterator%ipos)
      end if
    else
      if (iand(riterator%iSpec, ALST_LSPEC_REVERSE).eq.0) then
        riterator%ipos = riterator%p_rarraylist%p_Ktable(ALST_HEAD,riterator%itable)
      else
        riterator%ipos = riterator%p_rarraylist%p_Ktable(ALST_TAIL,riterator%itable)
      end if
    end if

  end subroutine

  !************************************************************************

!<subroutine>

  pure subroutine FEAT2_PP_TEMPLATE_TD(alst_prior,T,D)(riterator)

!<description>
    ! This subroutine decrements the list iterator by one.
!</description>

!<inputoutput>
    ! The iterator
    type(FEAT2_PP_TEMPLATE_TD(it_arraylist,T,D)), intent(inout) :: riterator
!</inputoutput>
!</subroutine>

    if (riterator%ipos .ne. ALST_NULL) then
      if (iand(riterator%iSpec, ALST_LSPEC_REVERSE).eq.0) then
        riterator%ipos = riterator%p_rarraylist%p_Kprev(riterator%ipos)
      else
        riterator%ipos = riterator%p_rarraylist%p_Knext(riterator%ipos)
      end if
    else
      if (iand(riterator%iSpec, ALST_LSPEC_REVERSE).eq.0) then
        riterator%ipos = riterator%p_rarraylist%p_Ktable(ALST_TAIL,riterator%itable)
      else
        riterator%ipos = riterator%p_rarraylist%p_Ktable(ALST_HEAD,riterator%itable)       
      end if
    end if

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine FEAT2_PP_TEMPLATE_TD(alst_getbase_key,T,D)(rarraylist, rposition, p_key)

!<description>
    ! This subroutine returns pointers to the key value stored at
    ! the position addressed by the iterator.
!</description>

!<input>
    ! The arraylist
    type(FEAT2_PP_TEMPLATE_TD(t_arraylist,T,D)), intent(in) :: rarraylist

    ! The iterator
    type(FEAT2_PP_TEMPLATE_TD(it_arraylist,T,D)), intent(in) :: rposition
!</input>

!<output>
    ! Pointer to the key value
    FEAT2_PP_TTYPE(T_TYPE), pointer :: p_key
!</output>
!</subroutine>

    ! Get key
    p_key => rarraylist%p_key(rposition%ipos)

  end subroutine

  !************************************************************************

#ifdef D
!<subroutine>

  subroutine FEAT2_PP_TEMPLATE_TD(alst_getbase_data,T,D)(rarraylist, rposition, p_data)

!<description>
    ! This subroutine returns pointers to the key value stored at
    ! the position addressed by the iterator.
!</description>

!<input>
    ! The arraylist
    type(FEAT2_PP_TEMPLATE_TD(t_arraylist,T,D)), intent(in) :: rarraylist

    ! The iterator
    type(FEAT2_PP_TEMPLATE_TD(it_arraylist,T,D)), intent(in) :: rposition
!</input>

!<output>
    ! Pointer to the auxiliary data value
    FEAT2_PP_DTYPE(D_TYPE), dimension(:), pointer :: p_data
!</output>
!</subroutine>

    if(rarraylist%isizeData > 0) then
      p_data => rarraylist%p_Data(:,rposition%ipos)
    else
      nullify(p_data)
    end if

  end subroutine
#endif

  !************************************************************************

!<function>

  function FEAT2_PP_TEMPLATE_TD(alst_get,T,D)(rarraylist, rposition) result(key)

!<description>
    ! This functions return the key stored at the position addressed
    ! by the iterator.
!</description>

!<input>
    ! The arraylist
    type(FEAT2_PP_TEMPLATE_TD(t_arraylist,T,D)), intent(in) :: rarraylist

    ! The iterator
    type(FEAT2_PP_TEMPLATE_TD(it_arraylist,T,D)), intent(in) :: rposition
!</input>

!<result>
    ! key value
    FEAT2_PP_TTYPE(T_TYPE) :: key
!</result>
!</function>

    ! Get key
    key = rarraylist%p_key(rposition%ipos)

  end function

  !************************************************************************

!<subroutine>

#ifdef D
  subroutine FEAT2_PP_TEMPLATE_TD(alst_assign1,T,D)(rarraylist, itable, n, key, data)
#else
  subroutine FEAT2_PP_TEMPLATE_TD(alst_assign1,T,D)(rarraylist, itable, n, key)
#endif

!<description>
    ! This subroutine removes all existing data from the list in the
    ! given table and assigns n copies of the given key and data values
!</description>

!<input>
    ! Number of the table
    integer, intent(in) :: itable

    ! Number of copies
    integer, intent(in) :: n

    ! key value
    FEAT2_PP_TTYPE(T_TYPE), intent(in) :: key

#ifdef D
    ! OPTIONAL: Data
    FEAT2_PP_DTYPE(D_TYPE), dimension(:), intent(in), optional :: data
#endif
!</input>

!<inputoutput>
    ! The arraylist
    type(FEAT2_PP_TEMPLATE_TD(t_arraylist,T,D)), intent(inout) :: rarraylist
!</inputoutput>
!</subroutine>

    ! local variable
    integer :: i

    ! Clear table in arraylist
    call alst_clearTbl(rarraylist, itable)

    ! Push copies to the list
    do i = 1, n
#ifdef D
      call alst_push_back(rarraylist, itable, key, data)
#else
      call alst_push_back(rarraylist, itable, key)
#endif
    end do

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine FEAT2_PP_TEMPLATE_TD(alst_assign2,T,D)(rarraylist, itable, rfirst, rlast)

!<description>
    ! This subroutine removes all existing data from the list in the
    ! given table and assigns the content in the range [rfirst,rlast)
    ! to the list.
!</description>

!<input>
    ! Number of the table
    integer, intent(in) :: itable

    ! Iterator referring to the first element
    type(FEAT2_PP_TEMPLATE_TD(it_arraylist,T,D)), intent(in) :: rfirst

    ! Iterator referring to the past-the-end element
    type(FEAT2_PP_TEMPLATE_TD(it_arraylist,T,D)), intent(in) :: rlast
!</input>

!<inputoutput>
    ! The arraylist
    type(FEAT2_PP_TEMPLATE_TD(t_arraylist,T,D)), intent(inout) :: rarraylist
!</inputoutput>
!</subroutine>

    ! local variable
    type(FEAT2_PP_TEMPLATE_TD(it_arraylist,T,D)) :: riterator
    FEAT2_PP_TTYPE(T_TYPE), pointer :: p_key

#ifdef D
    FEAT2_PP_DTYPE(D_TYPE), dimension(:), pointer :: p_data
#endif

    ! Clear list
    call alst_clearTbl(rarraylist, itable)

    ! Push content to the list
    riterator = rfirst
    do while(riterator /= rlast)

#ifdef D
      call alst_getbase_key(riterator%p_rarraylist, riterator, p_key)
      call alst_getbase_data(riterator%p_rarraylist, riterator, p_data)
      call alst_push_back(rarraylist, itable, p_key, p_data)
#else
      call alst_getbase_key(riterator%p_rarraylist, riterator, p_key)
      call alst_push_back(rarraylist, itable, p_key)
#endif
      call alst_next(riterator)

    end do

  end subroutine

  !************************************************************************

!<subroutine>

#ifdef D
  subroutine FEAT2_PP_TEMPLATE_TD(alst_push_front,T,D)(rarraylist, itable, key, data)
#else
  subroutine FEAT2_PP_TEMPLATE_TD(alst_push_front,T,D)(rarraylist, itable, key)
#endif

!<description>
    ! This subroutine prepends a new element to the list of the given table
!</description>

!<input>
    ! Number of the table
    integer, intent(in) :: itable

    ! key value
    FEAT2_PP_TTYPE(T_TYPE), intent(in) :: key

#ifdef D
    ! OPTIONAL: Data
    FEAT2_PP_DTYPE(D_TYPE), dimension(:), intent(in), optional :: data
#endif
!</input>

!<inputoutput>
    ! The arraylist
    type(FEAT2_PP_TEMPLATE_TD(t_arraylist,T,D)), intent(inout) :: rarraylist
!</inputoutput>
!</subroutine>

    ! local variable
    integer :: ipos,ihead

    ! Check if tables need to be created
    if (itable < 1) then
      call output_line('Invalid table number!',&
          OU_CLASS_ERROR,OU_MODE_STD,'alst_push_front')
      call sys_halt()
    end if
    if (rarraylist%ntable < itable) call alst_createTbl(rarraylist, itable)

    ! Check if list needs to be enlarged
    rarraylist%NA = rarraylist%NA+1
    rarraylist%p_Ktable(ALST_NA,itable) = rarraylist%p_Ktable(ALST_NA,itable)+1
    ipos          = rarraylist%p_Knext(ALST_FREE)
    if (abs(ipos) > rarraylist%NNA) then
      call alst_resize(rarraylist, ceiling(rarraylist%dfactor*rarraylist%NNA))
    end if

    ! Set next free position
    if (ipos > 0) then
      rarraylist%p_Knext(ALST_FREE) = ipos+1
    else
      ipos = abs(ipos)
      rarraylist%p_Knext(ALST_FREE) = rarraylist%p_Knext(ipos)
    end if

    ! Set head and tail
    if (rarraylist%p_Ktable(ALST_HEAD,itable) .eq. ALST_NULL) then
      ! Push to empty list
      rarraylist%p_Ktable(ALST_HEAD,itable) = ipos
      rarraylist%p_Ktable(ALST_TAIL,itable) = ipos
      rarraylist%p_Knext(ipos)              = ALST_NULL
      rarraylist%p_Kprev(ipos)              = ALST_NULL
    else
      ! Push to non-empty list
      ihead = rarraylist%p_Ktable(ALST_HEAD,itable)
      rarraylist%p_Kprev(ihead)             = ipos
      rarraylist%p_Knext(ipos)              = ihead
      rarraylist%p_Ktable(ALST_HEAD,itable) = ipos
      rarraylist%p_Kprev(ipos)              = ALST_NULL
    end if

    ! Store key value
    rarraylist%p_Key(ipos) = key

#ifdef D
    ! Store data value
    if (present(data)) then
      rarraylist%p_Data(:,ipos) = data
    end if
#endif

  end subroutine

  !************************************************************************

!<subroutine>

#ifdef D
  subroutine FEAT2_PP_TEMPLATE_TD(alst_push_back,T,D)(rarraylist, itable, key, data)
#else
  subroutine FEAT2_PP_TEMPLATE_TD(alst_push_back,T,D)(rarraylist, itable, key)
#endif

!<description>
    ! This subroutine appends a new element to the list of the given table
!</description>

!<input>
    ! Number of the table
    integer, intent(in) :: itable

    ! key value
    FEAT2_PP_TTYPE(T_TYPE), intent(in) :: key

#ifdef D
    ! OPTIONAL: Data
    FEAT2_PP_DTYPE(D_TYPE), dimension(:), intent(in), optional :: data
#endif
!</input>

!<inputoutput>
    ! The arraylist
    type(FEAT2_PP_TEMPLATE_TD(t_arraylist,T,D)), intent(inout) :: rarraylist
!</inputoutput>
!</subroutine>

    ! local variable
    integer :: ipos,itail

    ! Check if tables need to be created
    if (itable < 1) then
      call output_line('Invalid table number!',&
          OU_CLASS_ERROR,OU_MODE_STD,'alst_push_back')
      call sys_halt()
    end if
    if (rarraylist%ntable < itable) call alst_createTbl(rarraylist, itable)

    ! Check if list needs to be enlarged
    rarraylist%NA = rarraylist%NA+1
    rarraylist%p_Ktable(ALST_NA,itable) = rarraylist%p_Ktable(ALST_NA,itable)+1
    ipos          = rarraylist%p_Knext(ALST_FREE)
    if (abs(ipos) > rarraylist%NNA) then
      call alst_resize(rarraylist, ceiling(rarraylist%dfactor*rarraylist%NNA))
    end if

    ! Set next free position
    if (ipos > 0) then
      rarraylist%p_Knext(ALST_FREE) = ipos+1
    else
      ipos = abs(ipos)
      rarraylist%p_Knext(ALST_FREE) = rarraylist%p_Knext(ipos)
    end if

    ! Set head and tail
    if (rarraylist%p_Ktable(ALST_HEAD,itable) .eq. ALST_NULL) then
      ! Push to empty list
      rarraylist%p_Ktable(ALST_HEAD,itable) = ipos
      rarraylist%p_Ktable(ALST_TAIL,itable) = ipos
      rarraylist%p_Knext(ipos)              = ALST_NULL
      rarraylist%p_Kprev(ipos)              = ALST_NULL
    else
      ! Push to non-empty list
      itail = rarraylist%p_Ktable(ALST_TAIL,itable)
      rarraylist%p_Kprev(ipos)              = itail
      rarraylist%p_Knext(itail)             = ipos
      rarraylist%p_Ktable(ALST_TAIL,itable) = ipos
      rarraylist%p_Knext(ipos)              = ALST_NULL
    end if

    ! Store key value
    rarraylist%p_Key(ipos) = key

#ifdef D
    ! Store data value
    if (present(data)) then
      rarraylist%p_Data(:,ipos) = data
    end if
#endif

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine FEAT2_PP_TEMPLATE_TD(alst_pop_front,T,D)(rarraylist, itable)

!<description>
    ! This subroutine removes the first element from the list
    ! of the given table
!</description>

!<input>
    ! Number of the table
    integer, intent(in) :: itable
!</input>

!<inputoutput>
    ! The arraylist
    type(FEAT2_PP_TEMPLATE_TD(t_arraylist,T,D)), intent(inout) :: rarraylist
!</inputoutput>
!</subroutine>

    ! local variable
    type(FEAT2_PP_TEMPLATE_TD(it_arraylist,T,D)) :: riterator

    ! Check if list is empty
    if (alst_emptyTbl(rarraylist, itable)) return

    riterator = alst_begin(rarraylist, itable)
    riterator = alst_erase(rarraylist, riterator)

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine FEAT2_PP_TEMPLATE_TD(alst_pop_back,T,D)(rarraylist, itable)

!<description>
    ! This subroutine removes the last element from the list
    ! of the given table
!</description>

!<input>
    ! Number of the table
    integer, intent(in) :: itable
!</input>

!<inputoutput>
    ! The arraylist
    type(FEAT2_PP_TEMPLATE_TD(t_arraylist,T,D)), intent(inout) :: rarraylist
!</inputoutput>
!</subroutine>

    ! local variable
    type(FEAT2_PP_TEMPLATE_TD(it_arraylist,T,D)) :: riterator

    ! Check if list is empty
    if (alst_emptyTbl(rarraylist, itable)) return

    riterator = alst_rbegin(rarraylist, itable)
    riterator = alst_erase(rarraylist, riterator)

  end subroutine

  !************************************************************************

!<subroutine>

#ifdef D
  subroutine FEAT2_PP_TEMPLATE_TD(alst_getbase_front,T,D)(rarraylist, itable, p_key, p_data)
#else
  subroutine FEAT2_PP_TEMPLATE_TD(alst_getbase_front,T,D)(rarraylist, itable, p_key)
#endif

!<description>
    ! This subroutine returns pointers to the key and data stored in
    ! the first element of the given table
!</description>

!<input>
    ! The arraylist
    type(FEAT2_PP_TEMPLATE_TD(t_arraylist,T,D)), intent(in) :: rarraylist

    ! Number of the table
    integer, intent(in) :: itable
!</input>

!<output>
    ! Pointer to key value
    FEAT2_PP_TTYPE(T_TYPE), pointer :: p_key

#ifdef D
    ! OPTIONAL: Pointer to data value
    FEAT2_PP_DTYPE(D_TYPE), dimension(:), pointer, optional :: p_data
#endif
!</output>
!</subroutine>

    ! local variable
    type(FEAT2_PP_TEMPLATE_TD(it_arraylist,T,D)) :: riterator

    ! Check if list is empty
    if (alst_emptyTbl(rarraylist, itable)) then

      nullify(p_key)
#ifdef D
      nullify(p_data)
#endif

    else

      riterator = alst_begin(rarraylist, itable)

#ifdef D
      call alst_getbase_key(rarraylist, riterator, p_key)
      call alst_getbase_data(rarraylist, riterator, p_data)
#else
      call alst_getbase_key(rarraylist, riterator, p_key)
#endif

    end if

  end subroutine

  !************************************************************************

!<function>

  function FEAT2_PP_TEMPLATE_TD(alst_front,T,D)(rarraylist, itable) result(key)

!<description>
    ! This function returns the key stored in the first element of the
    ! given table
!</description>

!<input>
    ! The arraylist
    type(FEAT2_PP_TEMPLATE_TD(t_arraylist,T,D)), intent(in) :: rarraylist

    ! Number of the table
    integer, intent(in) :: itable
!</input>

!<result>
    ! key value
    FEAT2_PP_TTYPE(T_TYPE) :: key
!</result>
!</function>

    ! local variable
    type(FEAT2_PP_TEMPLATE_TD(it_arraylist,T,D)) :: riterator

    ! Check if list is empty
    if (alst_emptyTbl(rarraylist, itable)) return

    riterator = alst_begin(rarraylist, itable)
    key = alst_get(rarraylist, riterator)

  end function

  !************************************************************************

!<subroutine>

#ifdef D
  subroutine FEAT2_PP_TEMPLATE_TD(alst_getbase_back,T,D)(rarraylist, itable, p_key, p_data)
#else
  subroutine FEAT2_PP_TEMPLATE_TD(alst_getbase_back,T,D)(rarraylist, itable, p_key)
#endif

!<description>
    ! This subroutine returns pointers to the key and data stored in
    ! the first element of the given table
!</description>

!<input>
    ! The arraylist
    type(FEAT2_PP_TEMPLATE_TD(t_arraylist,T,D)), intent(in) :: rarraylist

    ! Number of the table
    integer, intent(in) :: itable
!</input>

!<output>
    ! Pointer to key value
    FEAT2_PP_TTYPE(T_TYPE), pointer :: p_key

#ifdef D
    ! OPTIONAL: Pointer to data
    FEAT2_PP_DTYPE(D_TYPE), dimension(:), pointer, optional :: p_data
#endif
!</output>
!</subroutine>

    ! local variable
    type(FEAT2_PP_TEMPLATE_TD(it_arraylist,T,D)) :: riterator

    ! Check if list is empty
    if (alst_emptyTbl(rarraylist, itable)) then

      nullify(p_key)
#ifdef D
      nullify(p_data)
#endif

    else

      riterator = alst_rbegin(rarraylist, itable)

#ifdef D
      call alst_getbase_key(rarraylist, riterator, p_key)
      call alst_getbase_data(rarraylist, riterator, p_data)
#else
      call alst_getbase_key(rarraylist, riterator, p_key)
#endif

    end if

  end subroutine

  !************************************************************************

!<function>

  function FEAT2_PP_TEMPLATE_TD(alst_back,T,D)(rarraylist, itable) result(key)

!<description>
    ! This function returns the key stored in the last element of the
    ! given table
!</description>

!<input>
    ! The arraylist
    type(FEAT2_PP_TEMPLATE_TD(t_arraylist,T,D)), intent(in) :: rarraylist

    ! Number of the table
    integer, intent(in) :: itable
!</input>

!<result>
    ! key value
    FEAT2_PP_TTYPE(T_TYPE) :: key
!</result>
!</function>

    ! local variable
    type(FEAT2_PP_TEMPLATE_TD(it_arraylist,T,D)) :: riterator

    ! Check if list is empty
    if (alst_emptyTbl(rarraylist, itable)) return

    riterator = alst_rbegin(rarraylist, itable)
    key = alst_get(rarraylist, riterator)

  end function

  !************************************************************************

!<function>

#ifdef D
  function FEAT2_PP_TEMPLATE_TD(alst_insert1,T,D)(rarraylist, rposition, key, data)&
                                         result(riterator)
#else
  function FEAT2_PP_TEMPLATE_TD(alst_insert1,T,D)(rarraylist, rposition, key)&
                                         result(riterator)
#endif

!<description>
    ! This subroutine inserts a new element into the list of the given
    ! table before the element at rposition and returns an iterator
    ! pointing to the new element
!</description>

!<input>
    ! The iterator
    type(FEAT2_PP_TEMPLATE_TD(it_arraylist,T,D)), intent(in) :: rposition

    ! key value
    FEAT2_PP_TTYPE(T_TYPE), intent(in) :: key

#ifdef D
    ! OPTIONAL: Data
    FEAT2_PP_DTYPE(D_TYPE), dimension(:), intent(in), optional :: data
#endif
!</input>

!<inputoutput>
    ! The arraylist
    type(FEAT2_PP_TEMPLATE_TD(t_arraylist,T,D)), intent(inout) :: rarraylist
!</inputoutput>

!<result>
    ! The new iterator
    type(FEAT2_PP_TEMPLATE_TD(it_arraylist,T,D)) :: riterator
!</result>
!</function>

    ! local variable
    integer :: ipred,ipos,itable

    ! Get table number
    itable = rposition%itable

    ! Insert new element at the end of the list?
    if (rposition%ipos .eq. ALST_NULL) then
#ifdef D
      call alst_push_back(rarraylist, itable, key, data)
#else
      call alst_push_back(rarraylist, itable, key)
#endif
      ! Create new iterator
      riterator       = rposition
      riterator%ipos  = rarraylist%p_Ktable(ALST_HEAD,itable)
      riterator%iSpec = iand(riterator%iSpec,&
                             not(ALST_LSPEC_REVERSE))
      ! That`s it
      return
    end if

    ! Get position of predecessor
    ipred = rarraylist%p_Kprev(rposition%ipos)

    ! Insert new element at the beginning of the list?
    if (ipred .eq. ALST_NULL) then
#ifdef D
      call alst_push_front(rarraylist, itable, key, data)
#else
      call alst_push_front(rarraylist, itable, key)
#endif
      ! Create new iterator
      riterator       = rposition
      riterator%ipos  = rarraylist%p_Ktable(ALST_TAIL,itable)
      riterator%iSpec = iand(riterator%iSpec,&
                             not(ALST_LSPEC_REVERSE))
      ! That`s it
      return
    end if

    ! Check if tables need to be created
    if (rarraylist%ntable < itable)&
        call alst_createTbl(rarraylist, itable)

    ! Check if array list needs to be enlarged
    rarraylist%NA = rarraylist%NA+1
    ipos          = rarraylist%p_Knext(ALST_FREE)
    rarraylist%p_Ktable(ALST_NA,itable) = rarraylist%p_Ktable(ALST_NA,itable)+1
    if (abs(ipos) > rarraylist%NNA) then
      call alst_resize(rarraylist, ceiling(rarraylist%dfactor*rarraylist%NNA))
    end if

    ! Set next free position
    if (ipos > 0) then
      rarraylist%p_Knext(ALST_FREE) = ipos+1
    else
      ipos = abs(ipos)
      rarraylist%p_Knext(ALST_FREE) = rarraylist%p_Knext(ipos)
    end if

    ! Insert element between its predecessor IPRED and 
    ! the element at position rposition in table itable
    rarraylist%p_Kprev(rposition%ipos) = ipos
    rarraylist%p_Knext(ipos)  = rarraylist%p_Knext(ipred)
    rarraylist%p_Knext(ipred) = ipos
    rarraylist%p_Kprev(ipos)  = ipred

    ! Store key value
    rarraylist%p_Key(ipos) = key

#ifdef D
    ! Store data value
    if (present(data)) then
      rarraylist%p_Data(:,ipos) = data
    end if
#endif

    ! Create new iterator
    riterator       = rposition
    riterator%ipos  = ipos
    riterator%iSpec = iand(riterator%iSpec,&
                           not(ALST_LSPEC_REVERSE))

  end function

  !************************************************************************

!<function>

#ifdef D
  function FEAT2_PP_TEMPLATE_TD(alst_insert2,T,D)(rarraylist, rposition, n, key, data)&
                                         result(riterator)
#else
  function FEAT2_PP_TEMPLATE_TD(alst_insert2,T,D)(rarraylist, rposition, n, key)&
                                         result(riterator)
#endif

!<description>
    ! This subroutine inserts n copies of the given key and data
    ! before the element at rposition.
!</description>

!<input>
    ! The iterator
    type(FEAT2_PP_TEMPLATE_TD(it_arraylist,T,D)), intent(in) :: rposition

    ! Number of copies
    integer, intent(in) :: n

    ! key value
    FEAT2_PP_TTYPE(T_TYPE), intent(in) :: key

#ifdef D
    ! OPTIONAL: Data
    FEAT2_PP_DTYPE(D_TYPE), dimension(:), intent(in), optional :: data
#endif
!</input>

!<inputoutput>
    ! The arraylist
    type(FEAT2_PP_TEMPLATE_TD(t_arraylist,T,D)), intent(inout) :: rarraylist
!</inputoutput>

!<result>
    ! The new iterator
    type(FEAT2_PP_TEMPLATE_TD(it_arraylist,T,D)) :: riterator
!</result>
!</function>

    ! local variable
    integer :: i

    do i = 1, n
#ifdef D
      riterator = alst_insert(rarraylist, rposition, key, data)
#else
      riterator = alst_insert(rarraylist, rposition, key)
#endif
    end do

  end function

  !************************************************************************

!<function>

  function FEAT2_PP_TEMPLATE_TD(alst_insert3,T,D)(rarraylist, rposition, rfirst, rlast)&
                                         result(riterator)

!<description>
    ! This subroutine inserts content in the range [rfirst,rlast)
    ! before the element at rposition.
!</description>

!<input>
    ! The iterator
    type(FEAT2_PP_TEMPLATE_TD(it_arraylist,T,D)), intent(in) :: rposition

    ! Iterator referring to the first element
    type(FEAT2_PP_TEMPLATE_TD(it_arraylist,T,D)), intent(in) :: rfirst

    ! Iterator referring to the past-the-last element
    type(FEAT2_PP_TEMPLATE_TD(it_arraylist,T,D)), intent(in) :: rlast
!</input>

!<inputoutput>
    ! The arraylist
    type(FEAT2_PP_TEMPLATE_TD(t_arraylist,T,D)), intent(inout) :: rarraylist
!</inputoutput>

!<result>
    ! The new iterator
    type(FEAT2_PP_TEMPLATE_TD(it_arraylist,T,D)) :: riterator
!</result>
!</function>

    ! local variable
    type(FEAT2_PP_TEMPLATE_TD(it_arraylist,T,D)) :: riter
    FEAT2_PP_TTYPE(T_TYPE), pointer :: p_key

#ifdef D
    FEAT2_PP_DTYPE(D_TYPE), dimension(:), pointer :: p_data
#endif

    ! Push content to the list
    riter = rfirst
    do while(riter /= rlast)

#ifdef D
      call alst_getbase_key(riter%p_rarraylist, riter, p_key)
      call alst_getbase_data(riter%p_rarraylist, riter, p_data)
      riterator = alst_insert(rarraylist, rposition, p_key, p_data)
#else
      call alst_getbase_key(riter%p_rarraylist, riter, p_key)
      riterator = alst_insert(rarraylist, rposition, p_key)
#endif
      call alst_next(riter)
    end do

  end function

  !************************************************************************

!<function>

  function FEAT2_PP_TEMPLATE_TD(alst_erase1,T,D)(rarraylist, rposition) result(riterator)

!<description>
    ! This function removes the element from the arraylist at
    ! rposition and returns an iterator pointing to the element
    ! following the removed element.
!</description>

!<input>
    ! The iterator
    type(FEAT2_PP_TEMPLATE_TD(it_arraylist,T,D)) :: rposition
!</input>

!<inputoutput>
    ! The arraylist
    type(FEAT2_PP_TEMPLATE_TD(t_arraylist,T,D)), intent(inout) :: rarraylist
!</inputoutput>

!<result>
    ! The new iterator
    type(FEAT2_PP_TEMPLATE_TD(it_arraylist,T,D)) :: riterator
!</result>
!</function>

    ! local variables
    integer :: inext,ipos,ipred,itable

    ! Get position data
    ipos   = rposition%ipos
    ipred  = rarraylist%p_Kprev(ipos)
    itable = rposition%itable

    ! Delete entry from arraylist
    rarraylist%NA = rarraylist%NA-1
    rarraylist%p_Ktable(ALST_NA,itable) = rarraylist%p_Ktable(ALST_NA,itable)-1

    ! Are we at the head of the list?
    if (ipred .eq. ALST_NULL) then
      ! Erase first list item
      inext = rarraylist%p_Knext(ipos)
      rarraylist%p_Ktable(ALST_HEAD,itable) = inext
    else
      ! Erase item inside the list
      inext = rarraylist%p_Knext(ipos)
      rarraylist%p_Knext(ipred) = inext
    end if

    ! Are we at the tail of the list
    if (inext .eq. ALST_NULL) then
      rarraylist%p_Ktable(ALST_TAIL,itable) = ipred
    else
      rarraylist%p_Kprev(inext) = ipred
    end if

    ! Update free position
    rarraylist%p_Knext(ipos) = rarraylist%p_Knext(ALST_FREE)
    rarraylist%p_Knext(ALST_FREE) = -ipos

    ! Create new iterator
    riterator       = rposition
    riterator%ipos  = inext
    riterator%iSpec = iand(riterator%iSpec,&
                           not(ALST_LSPEC_REVERSE))

  end function

  !************************************************************************

!<function>

  function FEAT2_PP_TEMPLATE_TD(alst_erase2,T,D)(rarraylist, rfirst, rlast) result(riterator)

!<description>
    ! This function removes the elements [rfirst,rlast) from the list
    ! and returns an iterator referring to the element following the
    ! last removed element.
!</description>

!<input>
    ! Iterator referring to the first element
    type(FEAT2_PP_TEMPLATE_TD(it_arraylist,T,D)), intent(in) :: rfirst

    ! Iterator referring to the past-the-last element
    type(FEAT2_PP_TEMPLATE_TD(it_arraylist,T,D)), intent(in) :: rlast
!</input>

!<inputoutput>
    ! The arraylist
    type(FEAT2_PP_TEMPLATE_TD(t_arraylist,T,D)), intent(inout) :: rarraylist
!</inputoutput>

!<result>
    ! The new iterator
    type(FEAT2_PP_TEMPLATE_TD(it_arraylist,T,D)) :: riterator
!</result>
!</function>

    ! Remove elements
    riterator = rfirst
    do while(riterator /= rlast)
      riterator = alst_erase(rarraylist, riterator)
    end do

  end function

  !************************************************************************

!<function>

  recursive function FEAT2_PP_TEMPLATE_TD(alst_find1,T,D)(rarraylist, itable, key,&
                                                 rpositionGuess) result(riterator)

!<description>
    ! This function searches for a given key value in the list of the
    ! given table and returns its position in the list. If the key
    ! value is not present in the list then NULL is returned.
!</description>

!<input>
    ! key value
    FEAT2_PP_TTYPE(T_TYPE), intent(in) :: key

    ! Number of the table
    integer, intent(in) :: itable

    ! The arraylist
    type(FEAT2_PP_TEMPLATE_TD(t_arraylist,T,D)), intent(in) :: rarraylist

    ! OPTIONAL: initial guess of the position
    type(FEAT2_PP_TEMPLATE_TD(it_arraylist,T,D)), intent(in), optional :: rpositionGuess
!</input>

!<result>
    ! The iterator
    type(FEAT2_PP_TEMPLATE_TD(it_arraylist,T,D)) :: riterator
!</result>
!</function>

    ! local variables
    logical :: bsuccess

    ! Do we have an initial guess
    if (present(rpositionGuess)) then
      riterator = rpositionGuess
    else
      riterator = alst_begin(rarraylist, itable)
    end if

    ! Initialisation
    bsuccess = .false.

    do while(riterator%ipos .ne. ALST_NULL)
      if (rarraylist%p_Key(riterator%ipos) .eq. key) then
        bsuccess = .true.; exit
      end if
      riterator%ipos = rarraylist%p_Knext(riterator%ipos)
    end do

    if (.not.bsuccess) then
      riterator%ipos = ALST_NULL

      ! If we started searching with an initial guess but we were
      ! unable to find the element then re-search from scratch
      if (present(rpositionGuess))&
          riterator = alst_find(rarraylist, itable, key)
    end if

  end function

!************************************************************************

!<function>

  recursive function FEAT2_PP_TEMPLATE_TD(alst_find2,T,D)(rarraylist, itable, key,&
                                                 bisReverse, rpositionGuess)&
                                                 result(riterator)

!<description>

    ! This function searches for a given key value in a sorted list of
    ! the given table and returns its position in the list. If the key
    ! value is not present in the list then riterator refers to the
    ! element which lexicographically follows the key which was not
    ! found.
    !
    ! Note that this function requires the array list to be sorted
    ! without checking this. Setting the parameter bisReverse to
    ! .true.  indicates that the list is sorted in reverse order.
!</description>

!<input>
    ! key value
    FEAT2_PP_TTYPE(T_TYPE), intent(in) :: key

    ! Number of the table
    integer, intent(in) :: itable

    ! The arraylist
    type(FEAT2_PP_TEMPLATE_TD(t_arraylist,T,D)), intent(in) :: rarraylist

    ! Flag: if TRUE then list is assumed in reverse order
    logical, intent(in) :: bisReverse

    ! OPTIONAL: initial guess of the position
    type(FEAT2_PP_TEMPLATE_TD(it_arraylist,T,D)), intent(in), optional :: rpositionGuess
!</input>

!<result>
    ! The iterator
    type(FEAT2_PP_TEMPLATE_TD(it_arraylist,T,D)) :: riterator
!</result>
!</function>

    ! local variables
    logical :: bsuccess

    ! Do we have an initial guess
    if (present(rpositionGuess)) then
      riterator = rpositionGuess
    else
      riterator = alst_begin(rarraylist, itable)
    end if

    ! Initialisation
    bsuccess = .false.


    if (bisReverse) then
      ! Search for key value in reverse order
      do while(riterator%ipos .ne. ALST_NULL)
        if (rarraylist%p_Key(riterator%ipos) .eq. key) then
          bsuccess = .true.; exit
        elseif (rarraylist%p_Key(riterator%ipos) .lt. key) then
          bsuccess = .false.; exit
        end if
        riterator%ipos = rarraylist%p_Knext(riterator%ipos)
      end do
    else
      ! Search for key value in default order
      do while(riterator%ipos .ne. ALST_NULL)
        if (rarraylist%p_Key(riterator%ipos) .eq. key) then
          bsuccess = .true.; exit
        elseif (rarraylist%p_Key(riterator%ipos) .gt. key) then
          bsuccess = .false.; exit
        end if
        riterator%ipos = rarraylist%p_Knext(riterator%ipos)
      end do
    end if

    if (.not.bsuccess) then
      ! Special treatment of position
      riterator%iSpec = riterator%iSpec + ALST_LSPEC_VIRTUAL

      ! If we started searching with an initial guess but we were
      ! unable to find the element then re-search from scratch
      if (present(rpositionGuess))&
          riterator = alst_find(rarraylist, itable, key, bisReverse)
    end if

  end function

  !************************************************************************

!<subroutine>

  subroutine FEAT2_PP_TEMPLATE_TD(alst_print,T,D)(rarraylist)

!<description>
    ! This subroutine prints the content of the arraylist
!</description>

!<input>
    ! The arraylist
    type(FEAT2_PP_TEMPLATE_TD(t_arraylist,T,D)), intent(in) :: rarraylist
!</input>
!</subroutine>

#ifdef T_STORAGE

    ! local variable
    type(FEAT2_PP_TEMPLATE_TD(it_arraylist,T,D)) :: riterator
    integer :: itable

    ! Loop over all tables
    do itable = 1, rarraylist%ntable
      call output_line('Table number: ' // trim(sys_siL(itable,15)))
      call output_line('----------------------------------------')

      riterator = alst_begin(rarraylist, itable)
      do while (riterator /= alst_end(rarraylist, itable))
        write(*,*) rarraylist%p_Key(riterator%ipos)
        call alst_next(riterator)
      end do
    end do

#else

    call output_line('Unable to print arraylist with derived data type!',&
        OU_CLASS_WARNING,OU_MODE_STD,'alst_print')

#endif

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine FEAT2_PP_TEMPLATE_TD(alst_printTbl,T,D)(rarraylist, itable)

!<description>
    ! This subroutine prints the content of the given table of the arraylist
!</description>

!<input>
    ! The arraylist
    type(FEAT2_PP_TEMPLATE_TD(t_arraylist,T,D)), intent(in) :: rarraylist

    ! Number of the table
    integer, intent(in) :: itable
!</input>
!</subroutine>

#ifdef T_STORAGE

    ! local variable
    type(FEAT2_PP_TEMPLATE_TD(it_arraylist,T,D)) :: riterator

    riterator = alst_begin(rarraylist, itable)
    do while (riterator /= alst_end(rarraylist, itable))
      write(*,*) rarraylist%p_Key(riterator%ipos)
      call alst_next(riterator)
    end do

#else

    call output_line('Unable to print arraylist with derived data type!',&
        OU_CLASS_WARNING,OU_MODE_STD,'alst_printTbl')

#endif

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine FEAT2_PP_TEMPLATE_TD(alst_info,T,D)(rarraylist)

!<description>
    ! This subroutine prints information about the arraylist
!</description>

!<input>
    ! arraylist
    type(FEAT2_PP_TEMPLATE_TD(t_arraylist,T,D)), intent(in) :: rarraylist
!</input>
!</subroutine>

    call output_line('Arraylist:')
    call output_line('----------')
    call output_line('NA:       '//trim(sys_siL(rarraylist%NA,15)))
    call output_line('NNA:      '//trim(sys_siL(rarraylist%NNA,15)))
    call output_line('NNA0:     '//trim(sys_siL(rarraylist%NNA0,15)))
    call output_line('ntable:   '//trim(sys_siL(rarraylist%ntable,15)))
    call output_line('nntable:  '//trim(sys_siL(rarraylist%nntable,15)))
    call output_line('nntable0: '//trim(sys_siL(rarraylist%nntable0,15)))
    call output_line('nresize:  '//trim(sys_siL(rarraylist%nresize,15)))
    call output_line('dfactor:  '//trim(sys_sdL(rarraylist%dfactor,2)))
    call output_line('h_Ktable: '//trim(sys_siL(rarraylist%h_Ktable,15)))
    call output_line('h_Knext:  '//trim(sys_siL(rarraylist%h_Knext,15)))
    call output_line('h_Kprev:  '//trim(sys_siL(rarraylist%h_Kprev,15)))

#ifdef T_STORAGE
    call output_line('h_Key:    '//trim(sys_siL(rarraylist%h_Key,15)))
#endif
#if defined D && defined D_STORAGE
    call output_line('h_Data:   '//trim(sys_siL(rarraylist%h_Data,15)))
#endif

    call output_lbrk()

    call output_line('Current data  memory usage: '//&
        trim(sys_sdL(100*rarraylist%NA/&
        real(max(1,rarraylist%NNA),DP),2))//'%')
    call output_line('Current table memory usage: '//&
        trim(sys_sdL(100*rarraylist%ntable/&
        real(max(1,rarraylist%nntable),DP),2))//'%')
    call output_line('Total data memory usage:    '//&
        trim(sys_sdL(100*rarraylist%NNA/&
        real(max(1,rarraylist%NNA0),DP),2))//'%')
    call output_line('Total table memory usage:   '//&
        trim(sys_sdL(100*rarraylist%nntable/&
        real(max(1,rarraylist%nntable0),DP),2))//'%')

    call output_lbrk()

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine FEAT2_PP_TEMPLATE_TD(alst_duplicate,T,D)(rarraylist, rarraylistBackup)

!<description>
    ! This subroutine makes a copy of an arraylist in memory.
    ! It does not make sense to share some information between arraylists,
    ! so each array is physically copied from the source rarraylist
    ! to the destination rarraylistBackup.
!</description>

!<input>
    ! The rraylist for which a backup copy should be generated
    type(FEAT2_PP_TEMPLATE_TD(t_arraylist,T,D)), intent(in) :: rarraylist
!</input>

!<inputoutput>
    ! Backup copy arraylist
    type(FEAT2_PP_TEMPLATE_TD(t_arraylist,T,D)), intent(inout) :: rarraylistBackup
!</inputoutput>
!</subroutine>

    ! Release backup arraylist
    call alst_release(rarraylistBackup)

    ! Copy all data
    rarraylistBackup = rarraylist

    ! Reset handles
    rarraylistBackup%h_Ktable = ST_NOHANDLE
    rarraylistBackup%h_Knext  = ST_NOHANDLE
    rarraylistBackup%h_Kprev  = ST_NOHANDLE

    ! Copy storage blocks
    if (rarraylist%h_Ktable .ne. ST_NOHANDLE) then
      call storage_copy(rarraylist%h_Ktable, rarraylistBackup%h_Ktable)
      call storage_getbase_int2D(rarraylistBackup%h_Ktable,&
                                 rarraylistBackup%p_Ktable)
    end if

    if (rarraylist%h_Knext .ne. ST_NOHANDLE) then
      call storage_copy(rarraylist%h_Knext, rarraylistBackup%h_Knext)
      call storage_getbase_int(rarraylistBackup%h_Knext,&
                               rarraylistBackup%p_Knext)
    end if

    if (rarraylist%h_Kprev .ne. ST_NOHANDLE) then
      call storage_copy(rarraylist%h_Kprev, rarraylistBackup%h_Kprev)
      call storage_getbase_int(rarraylistBackup%h_Kprev,&
                               rarraylistBackup%p_Kprev)
    end if

#ifdef T_STORAGE
    rarraylistBackup%h_Key = ST_NOHANDLE
    if (rarraylist%h_Key .ne. ST_NOHANDLE) then
      call storage_copy(rarraylist%h_Key, rarraylistBackup%h_Key)
      call storage_getbase(rarraylistBackup%h_Key,&
                           rarraylistBackup%p_Key)
    end if
#else
    if (associated(rarraylist%p_Key)) then
      allocate(rarraylistBackup%p_Key(size(rarraylist%p_Key)))
      rarraylistBackup%p_Key = rarraylist%p_Key
    else
      nullify(rarraylistBackup%p_Key)
    end if
#endif

#ifdef D
#ifdef D_STORAGE
    rarraylistBackup%h_Data = ST_NOHANDLE
    if (rarraylist%h_Data .ne. ST_NOHANDLE) then
      call storage_copy(rarraylist%h_Data, rarraylistBackup%h_Data)
      call storage_getbase(rarraylistBackup%h_Data,&
                           rarraylistBackup%p_Data)
    end if
#else
    if (associated(rarraylist%p_Data)) then
      allocate(rarraylistBackup%p_Data(size(rarraylist%p_Data,1),&
                                       size(rarraylist%p_Data,2)))
      rarraylistBackup%p_Data = rarraylist%p_Data
    else
      nullify(rarraylistBackup%p_Data)
    end if
#endif
#endif

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine FEAT2_PP_TEMPLATE_TD(alst_restore,T,D)(rarraylistBackup, rarraylist)

!<description>
    ! This subroutine restores an arraylist from a previous backup.
!</description>

!<input>
    ! Backup copy of the arraylist
    type(FEAT2_PP_TEMPLATE_TD(t_arraylist,T,D)), intent(in) :: rarraylistBackup
!</input>

!<inputoutput>
    ! Destination arraylist
    type(FEAT2_PP_TEMPLATE_TD(t_arraylist,T,D)), intent(inout) :: rarraylist
!</inputoutput>
!</subroutine>

    ! Release arraylist
    call alst_release(rarraylist)

    ! Duplicate the backup
    call alst_duplicate(rarraylistBackup, rarraylist)

  end subroutine

  !************************************************************************

!<function>

  pure function FEAT2_PP_TEMPLATE_TD(alst_size,T,D)(rarraylist) result(isize)

!<description>
    ! Returns the size of the arraylist
!</description>

!<input>
    ! The arraylist
    type(FEAT2_PP_TEMPLATE_TD(t_arraylist,T,D)), intent(in) :: rarraylist
!</input>

!<result>
    ! Size of the arraylist
    integer :: isize
!</result>
!</function>

    isize = rarraylist%NA

  end function

  !************************************************************************

!<function>

  pure function FEAT2_PP_TEMPLATE_TD(alst_max_size,T,D)(rarraylist) result(imaxsize)

!<description>
    ! Returns the maximum size of the arraylist
!</description>

!<input>
    ! The arraylist
    type(FEAT2_PP_TEMPLATE_TD(t_arraylist,T,D)), intent(in) :: rarraylist
!</input>

!<result>
    ! Maximum size of the arraylist
    integer :: imaxsize
!</result>
!</function>

    imaxsize = rarraylist%NNA

  end function

  !************************************************************************

!<function>

  pure function FEAT2_PP_TEMPLATE_TD(alst_ntable,T,D)(rarraylist) result(ntable)

!<description>
    ! Returns the number of tables of the arraylist
!</description>

!<input>
    ! The arraylist
    type(FEAT2_PP_TEMPLATE_TD(t_arraylist,T,D)), intent(in) :: rarraylist
!</input>

!<result>
    ! Number of tables
    integer :: ntable
!</result>
!</function>

    ntable = rarraylist%ntable

  end function

  !************************************************************************

!<function>

  pure function FEAT2_PP_TEMPLATE_TD(alst_empty,T,D)(rarraylist) result(bempty)

!<description>
    ! Checks if the arraylist is empty
!</description>

!<input>
    ! The arraylist
    type(FEAT2_PP_TEMPLATE_TD(t_arraylist,T,D)), intent(in) :: rarraylist
!</input>

!<result>
    ! Flag: is true if the list is empty
    logical :: bempty
!</result>
!</function>

    bempty = (rarraylist%NA .eq. 0)

  end function

  !************************************************************************

!<function>

  pure function FEAT2_PP_TEMPLATE_TD(alst_emptyTbl,T,D)(rarraylist, itable) result(bempty)

!<description>
    ! Checks if the given table of the arraylist is empty
!</description>

!<input>
    ! The arraylist
    type(FEAT2_PP_TEMPLATE_TD(t_arraylist,T,D)), intent(in) :: rarraylist

    ! Number of the table
    integer, intent(in) :: itable
!</input>

!<result>
    ! Flag: is true if the list is empty
    logical :: bempty
!</result>
!</function>

    ! Check if table exists
    if (itable < 1 .or. itable > rarraylist%ntable) then
      bempty = .true.
    else
      bempty = (rarraylist%p_Ktable(ALST_NA,itable) .eq. 0)
    end if

  end function

  !************************************************************************

!<subroutine>

  subroutine FEAT2_PP_TEMPLATE_TD(alst_fassign,T,D)(rarraylistDest, rarraylistSrc)

!<description>
    ! Assigns the content of rarraylistSrc to rarraylistDest
!</description>

!<input>
    ! Source arraylist
    type(FEAT2_PP_TEMPLATE_TD(t_arraylist,T,D)), intent(in) :: rarraylistSrc
!</input>

!<output>
    ! Destination arraylist
    type(FEAT2_PP_TEMPLATE_TD(t_arraylist,T,D)), intent(out) :: rarraylistDest
!</output>
!</subroutine>


    ! Create empty arraylist
#ifdef D
    call alst_create(rarraylistDest, rarraylistSrc%nntable, rarraylistSrc%NNA,&
        rarraylistSrc%isizeData, rarraylistSrc%dfactor)
#else
    call alst_create(rarraylistDest, rarraylistSrc%nntable, rarraylistSrc%NNA,&
        rarraylistSrc%dfactor)
#endif

    ! Set size
    rarraylistDest%NA     = rarraylistsrc%NA
    rarraylistDest%ntable = rarraylistsrc%ntable

    ! Set structure
    rarraylistDest%p_Ktable = rarraylistSrc%p_Ktable
    rarraylistDest%p_Knext  = rarraylistSrc%p_Knext
    rarraylistDest%p_Kprev  = rarraylistSrc%p_Kprev
    rarraylistDest%p_Key    = rarraylistSrc%p_Key

#ifdef D
    rarraylistDest%p_Data   = rarraylistSrc%p_Data
#endif

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine FEAT2_PP_TEMPLATE_TD(alst_reverse,T,D)(rarraylist)

!<description>
    ! This subroutine reverses the ordering of the arraylist
!</description>

!<inputoutput>
    ! The arraylist to be reverted
    type(FEAT2_PP_TEMPLATE_TD(t_arraylist,T,D)), intent(inout) :: rarraylist
!</inputoutput>
!</subroutine>

    ! local variable
    integer :: itable

    ! Loop over all tables
    do itable = 1, rarraylist%ntable
      call alst_reverseTbl(rarraylist, itable)
    end do

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine FEAT2_PP_TEMPLATE_TD(alst_reverseTbl,T,D)(rarraylist, itable)

!<description>
    ! This subroutine reverses the ordering of the list of the given table
!</description>

!<input>
    ! Number of the table
    integer, intent(in) :: itable
!</input>

!<inputoutput>
    ! The arraylist to be reverted
    type(FEAT2_PP_TEMPLATE_TD(t_arraylist,T,D)), intent(inout) :: rarraylist
!</inputoutput>
!</subroutine>

    ! local variable
    type(FEAT2_PP_TEMPLATE_TD(it_arraylist,T,D)) :: riterator
    FEAT2_PP_TTYPE(T_TYPE), pointer :: p_key

#ifdef D
    FEAT2_PP_DTYPE(D_TYPE), dimension(:), pointer :: p_data
#endif

    riterator = alst_begin(rarraylist, itable)
    call alst_next(riterator)
    do while (riterator /= alst_end(rarraylist, itable))

#ifdef D
      call alst_getbase_key(rarraylist, riterator, p_key)
      call alst_getbase_data(rarraylist, riterator, p_data)
      call alst_push_front(rarraylist, itable, p_key, p_data)
#else
      call alst_getbase_key(rarraylist, riterator, p_key)
      call alst_push_front(rarraylist, itable, p_key)
#endif
      riterator = alst_erase(rarraylist, riterator)
    end do

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine FEAT2_PP_TEMPLATE_TD(alst_sort,T,D)(rarraylist)

!<description>
    ! This subroutine sorts the elements in the arraylist
!</description>

!<inputoutput>
    ! The arraylist to be sorted
    type(FEAT2_PP_TEMPLATE_TD(t_arraylist,T,D)), intent(inout) :: rarraylist
!</inputoutput>
!</subroutine>

    ! local variable
    integer :: itable

    ! Loop over all tables
    do itable = 1, rarraylist%ntable
      call alst_sortTbl(rarraylist, itable)
    end do

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine FEAT2_PP_TEMPLATE_TD(alst_sortTbl,T,D)(rarraylist, itable)

!<description>
    ! This subroutine sorts the elements in the list of the given table
!</description>

!<input>
    ! Number of the table
    integer, intent(in) :: itable
!</input>

!<inputoutput>
    ! The arraylist to be sorted
    type(FEAT2_PP_TEMPLATE_TD(t_arraylist,T,D)), intent(inout) :: rarraylist
!</inputoutput>
!</subroutine>

    ! local variable
    logical :: bswapped
    type(FEAT2_PP_TEMPLATE_TD(it_arraylist,T,D)) :: rposition,riterator
    FEAT2_PP_TTYPE(T_TYPE), pointer :: p_key1,p_key2
    FEAT2_PP_TTYPE(T_TYPE) :: key

#ifdef D
    FEAT2_PP_DTYPE(D_TYPE), dimension(:), pointer :: p_data1,p_data2
    FEAT2_PP_DTYPE(D_TYPE), dimension(:), allocatable :: data

    ! Allocate temporal memory
    if (rarraylist%isizeData > 0) allocate(data(rarraylist%isizeData))
#endif

    rposition = alst_rbegin(rarraylist, itable)

    ! Bubble sort algorithm
    outer: do
      bswapped  = .false.
      riterator = alst_begin(rarraylist, itable)
      inner: do while (riterator /= rposition)

#ifdef D
        if (rarraylist%isizeData > 0) then

          call alst_getbase_key(rarraylist,riterator,p_key1)
          call alst_getbase_data(rarraylist,riterator,p_data1)
          call alst_next(riterator)
          call alst_getbase_key(rarraylist,riterator,p_key2)
          call alst_getbase_data(rarraylist,riterator,p_data2)
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

          call alst_getbase_key(rarraylist,riterator,p_key1)
          call alst_next(riterator)
          call alst_getbase_key(rarraylist,riterator,p_key2)
          if (p_key1 > p_key2) then
            key    = p_key2
            p_key2 = p_key1
            p_key1 = key
            bswapped = .true.
          end if

        end if
#else
        call alst_getbase_key(rarraylist,riterator,p_key1)
        call alst_next(riterator)
        call alst_getbase_key(rarraylist,riterator,p_key2)
        if (p_key1 > p_key2) then
          key    = p_key2
          p_key2 = p_key1
          p_key1 = key
          bswapped = .true.
        end if
#endif

      end do inner
      call alst_next(rposition)
      if (.not.bswapped) exit
    end do outer

#ifdef D
    if (rarraylist%isizeData > 0) deallocate(data)
#endif

  end subroutine

  !************************************************************************

!<function>

  pure function FEAT2_PP_TEMPLATE_TD(alst_isNull,T,D)(riterator) result(bisNull)

!<description>
    ! Checks if the iterator is NULL
!</description>

!<input>
    ! Iterator
    type(FEAT2_PP_TEMPLATE_TD(it_arraylist,T,D)), intent(in) :: riterator
!</input>

!<result>
    ! True if the iterator is NULL.
    logical :: bisNull
!</result>
!</function>

    bisNull = ((riterator%ipos .eq. ALST_NULL) .or.&
           iand(riterator%iSpec, ALST_LSPEC_VIRTUAL) .ne. 0)

  end function

  !************************************************************************

!<function>

  pure function FEAT2_PP_TEMPLATE_TD(alst_hasSpec,T,D)(riterator,iSpec) result(bhasSpec)

!<description>
    ! Checks if the iterator has the given specification flag
!</description>

!<input>
    ! Iterator
    type(FEAT2_PP_TEMPLATE_TD(it_arraylist,T,D)), intent(in) :: riterator

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

  pure function FEAT2_PP_TEMPLATE_TD(it_alst_eq,T,D)(riterator1,riterator2) result(beq)

!<description>
    ! Compare two iterators for equality
!</description>

!<input>
    ! Iterators
    type(FEAT2_PP_TEMPLATE_TD(it_arraylist,T,D)), intent(in) :: riterator1,riterator2
!</input>

!<result>
    ! True if both iterators are equal
    logical :: beq
!</result>
!</function>

    ! Initialisation
    beq = (riterator1%ipos == riterator2%ipos)

  end function

  !************************************************************************

!<function>

  pure function FEAT2_PP_TEMPLATE_TD(it_alst_ne,T,D)(riterator1,riterator2) result(bne)

!<description>
    ! Compare two iterators for inequality
!</description>

!<input>
    ! Iterators
    type(FEAT2_PP_TEMPLATE_TD(it_arraylist,T,D)), intent(in) :: riterator1,riterator2
!</input>

!<result>
    ! True if both iterators are not equal
    logical :: bne
!</result>
!</function>

    ! Initialisation
    bne = (riterator1%ipos /= riterator2%ipos)

  end function

  !************************************************************************

!<function>

  pure function FEAT2_PP_TEMPLATE_TD(it_alst_lt,T,D)(riterator1,riterator2) result(blt)

!<description>
    ! Checks lexicographical ordering of two iterators
!</description>

!<input>
    ! Iterators
    type(FEAT2_PP_TEMPLATE_TD(it_arraylist,T,D)), intent(in) :: riterator1,riterator2
!</input>

!<result>
    ! True if the lexicographical ordering of iterator1 is smaller
    ! than that of iterator2
    logical :: blt
!</result>
!</function>

    ! Initialisation
    blt = (riterator1%ipos < riterator2%ipos)

  end function

  !************************************************************************

!<function>

  pure function FEAT2_PP_TEMPLATE_TD(it_alst_le,T,D)(riterator1,riterator2) result(ble)

!<description>
    ! Checks lexicographical ordering of two iterators
!</description>

!<input>
    ! Iterators
    type(FEAT2_PP_TEMPLATE_TD(it_arraylist,T,D)), intent(in) :: riterator1,riterator2
!</input>

!<result>
    ! True if the lexicographical ordering of iterator1 is smaller
    ! than or equal to that of iterator2
    logical :: ble
!</result>
!</function>

    ! Initialisation
    ble = (riterator1%ipos <= riterator2%ipos)

  end function

  !************************************************************************

!<function>

  pure function FEAT2_PP_TEMPLATE_TD(it_alst_gt,T,D)(riterator1,riterator2) result(bgt)

!<description>
    ! Checks lexicographical ordering of two iterators
!</description>

!<input>
    ! Iterators
    type(FEAT2_PP_TEMPLATE_TD(it_arraylist,T,D)), intent(in) :: riterator1,riterator2
!</input>

!<result>
    ! True if the lexicographical ordering of iterator1 is greater
    ! than that of iterator2
    logical :: bgt
!</result>
!</function>

    ! Initialisation
    bgt = (riterator1%ipos > riterator2%ipos)

  end function

  !************************************************************************

!<function>

  pure function FEAT2_PP_TEMPLATE_TD(it_alst_ge,T,D)(riterator1,riterator2) result(bge)

!<description>
    ! Checks lexicographical ordering of two iterators
!</description>

!<input>
    ! Iterators
    type(FEAT2_PP_TEMPLATE_TD(it_arraylist,T,D)), intent(in) :: riterator1,riterator2
!</input>

!<result>
    ! True if the lexicographical ordering of iterator1 is greater
    ! than or equal to that of iterator2
    logical :: bge
!</result>
!</function>

    ! Initialisation
    bge = (riterator1%ipos >= riterator2%ipos)

  end function

#endif
