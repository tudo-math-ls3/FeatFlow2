!-*- mode: f90; -*-

#ifndef _MAP_H_
#define _MAP_H_

!##############################################################################
!# ****************************************************************************
!# <name> FEAT2_PP_TEMPLATE_TD(map,T,D) </name>
!# ****************************************************************************
!#
!# <purpose>
!#
!# This module header file implements a map.
!#
!# Maps are a kind of container that contains a sorted list of unique
!# key-value pairs. Search, removal, and insertion operations have
!# algorithmic complexity. This implementatin is based on an AVL-tree.
!#
!# -----------------------------------------------------------------------------
!# This module is designed with greatest compatibility to the C++
!# Standard Template Library (STL) concerning naming conventions.
!# -----------------------------------------------------------------------------
!#
!# The following routines are available:
!#
!#  1.) map_create (constructor in STL)
!#      -> Creates an empty map
!#
!#  2.) map_release (destructor in STL)
!#      -> Releases a map
!#
!#  3.) map_resize (no equivalence in STL)
!#      -> Reallocates memory for a map
!#
!#  4.) map_clear (clear in STL)
!#      -> Removes all content from map
!#
!#  5.) map_copy (no equivalence in STL)
!#      -> Copies data to/from a map
!#
!# 6.) map_swap (swap in STL)
!#     -> Swaps two maps
!#
!# 7.) map_begin / map_rbegin (begin / rbegin in STL)
!#      -> Returns a (reverse) iterator referring to the first/last element 
!#         in a map
!#
!#  8.) map_end / map_rend (end / rend in STL)
!#      -> Returns a (reverse) iterator referring to the element right after
!#         the last/right before the first element in a map
!#
!#  9.) map_next (++iterator in STL)
!#      -> Increases the (reverse) iterator
!#
!# 10.) map_prior (--iterator in STL)
!#      -> Decreases the (reverse) iterator
!#
!# 11.) map_get (no equivalence in STL)
!#      -> Gets key value at given position in a map
!#
!# 12.) map_getbase_data (no equivalence in STL)
!#      -> Gets pointer to auxiliary data value
!#         at given position in a map
!#
!# 13.) map_insert (insert in STL)
!#      -> Inserts data into a map
!#
!# 14.) map_erase (erase in STL)
!#     -> Deletes data from a map
!#
!# 15.) map_size (size in STL)
!#      -> Returns the size of a map
!#
!# 16.) map_max_size (max_size in STL)
!#      -> Returns the maximum size of a map
!#
!# 17.) map_empty (empty in STL)
!#      -> Tests if the map is empty
!#
!# 18.) map_find (find in STL)
!#      -> Searches for data in a map
!#
!# 19.) map_print (no equivalence in STL)
!#      -> Prints content of a map
!#
!# 20.) map_info (no equivalence in STL)
!#      -> Outpus statistical information about a map
!#
!# 21.) map_duplicate (no equivalence in STL)
!#      -> Create a duplicate / backup of a map
!#
!# 22.) map_restore (no equivalence in STL)
!#      -> Restore a map from a previous backup
!#
!# 23.) map_isNull
!#      -> Tests if the iterator is NULL
!#
!# 24.) map_hasSpec
!#      -> Tests if the iterator has specification flags
!#
!#
!# The following operations in STL are not supported yet
!#
!# count, lower_bound, upper_bound, equal_range
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
  public :: FEAT2_PP_TEMPLATE_TD(t_map,T,D)
  public :: FEAT2_PP_TEMPLATE_TD(it_map,T,D)
  public :: map_create
  public :: map_release
  public :: map_resize
  public :: map_clear
  public :: map_copy
  public :: map_swap
  public :: map_begin
  public :: map_rbegin
  public :: map_end
  public :: map_rend
  public :: map_next
  public :: map_prior
  public :: map_get
#ifdef D
  public :: map_getbase_data
#endif
  public :: map_insert
  public :: map_erase
  public :: map_size
  public :: map_max_size
  public :: map_empty
  public :: map_find
  public :: map_print
  public :: map_info
  public :: map_duplicate
  public :: map_restore
  public :: map_isNull
  public :: map_hasSpec


  public assignment(=)
  public operator(==)
  public operator(/=)
  public operator(<)
  public operator(<=)
  public operator(>)
  public operator(>=)

  interface map_create
    module procedure FEAT2_PP_TEMPLATE_TD(map_create,T,D)
  end interface

  interface map_release
    module procedure FEAT2_PP_TEMPLATE_TD(map_release,T,D)
  end interface

  interface map_resize
    module procedure FEAT2_PP_TEMPLATE_TD(map_resize,T,D)
  end interface

  interface map_clear
    module procedure FEAT2_PP_TEMPLATE_TD(map_clear,T,D)
  end interface

  interface map_copy
    module procedure FEAT2_PP_TEMPLATE_TD(map_cpy1,T,D)
    module procedure FEAT2_PP_TEMPLATE_TD(map_cpy2,T,D)
    module procedure FEAT2_PP_TEMPLATE_TD(map_cpy3,T,D)
    module procedure FEAT2_PP_TEMPLATE_TD(map_cpy4,T,D)
  end interface

  interface map_swap
    module procedure FEAT2_PP_TEMPLATE_TD(map_swap,T,D)
  end interface
  
  interface map_begin
    module procedure FEAT2_PP_TEMPLATE_TD(map_begin,T,D)
  end interface

  interface map_rbegin
    module procedure FEAT2_PP_TEMPLATE_TD(map_rbegin,T,D)
  end interface
  
  interface map_end
    module procedure FEAT2_PP_TEMPLATE_TD(map_end,T,D)
  end interface

  interface map_rend
    module procedure FEAT2_PP_TEMPLATE_TD(map_rend,T,D)
  end interface

  interface map_next
    module procedure FEAT2_PP_TEMPLATE_TD(map_next,T,D)
  end interface

  interface map_prior
    module procedure FEAT2_PP_TEMPLATE_TD(map_prior,T,D)
  end interface

  interface map_get
    module procedure FEAT2_PP_TEMPLATE_TD(map_get,T,D)
  end interface

#ifdef D
  interface map_getbase_data
    module procedure FEAT2_PP_TEMPLATE_TD(map_getbase_data,T,D)
  end interface
#endif

  interface map_insert
    module procedure FEAT2_PP_TEMPLATE_TD(map_insert1,T,D)
    module procedure FEAT2_PP_TEMPLATE_TD(map_insert2,T,D)
  end interface
  
  interface map_erase
    module procedure FEAT2_PP_TEMPLATE_TD(map_erase1,T,D)
    module procedure FEAT2_PP_TEMPLATE_TD(map_erase2,T,D)
  end interface
  
  interface map_size
    module procedure FEAT2_PP_TEMPLATE_TD(map_size,T,D)
  end interface
  
  interface map_max_size
    module procedure FEAT2_PP_TEMPLATE_TD(map_max_size,T,D)
  end interface

  interface map_empty
    module procedure FEAT2_PP_TEMPLATE_TD(map_empty,T,D)
  end interface
  
  interface map_find
    module procedure FEAT2_PP_TEMPLATE_TD(map_find,T,D)
  end interface
 
  interface map_print
    module procedure FEAT2_PP_TEMPLATE_TD(map_print,T,D)
  end interface
  
  interface map_info
    module procedure FEAT2_PP_TEMPLATE_TD(map_info,T,D)
  end interface
  
  interface map_duplicate
    module procedure FEAT2_PP_TEMPLATE_TD(map_duplicate,T,D)
  end interface
  
  interface map_restore
    module procedure FEAT2_PP_TEMPLATE_TD(map_restore,T,D)
  end interface

  interface map_isNull
    module procedure FEAT2_PP_TEMPLATE_TD(map_isNull,T,D)
  end interface

  interface map_hasSpec
    module procedure FEAT2_PP_TEMPLATE_TD(map_hasSpec,T,D)
  end interface
  
  
  interface assignment(=)
    module procedure FEAT2_PP_TEMPLATE_TD(map_fassign,T,D)
  end interface


  interface operator(==)
    module procedure FEAT2_PP_TEMPLATE_TD(it_map_eq,T,D)
  end interface
  
  interface operator(/=)
    module procedure FEAT2_PP_TEMPLATE_TD(it_map_ne,T,D)
  end interface

  interface operator(<)
    module procedure FEAT2_PP_TEMPLATE_TD(it_map_lt,T,D)
  end interface

  interface operator(<=)
    module procedure FEAT2_PP_TEMPLATE_TD(it_map_le,T,D)
  end interface

  interface operator(>)
    module procedure FEAT2_PP_TEMPLATE_TD(it_map_gt,T,D)
  end interface

  interface operator(>=)
    module procedure FEAT2_PP_TEMPLATE_TD(it_map_ge,T,D)
  end interface

  !************************************************************************

!<types>

!<typeblock>

  type FEAT2_PP_TEMPLATE_TD(t_map,T,D)
    private
    
    ! The map is realised as AVL-tree. This special kind of binary
    ! search trees takes special care that the tree is well balanced,
    ! that is, the hieht of no two subtrees differs by more than two.

    ! Number of items that are currently stored in the map
    integer :: NA = 0

    ! Total number of items that can be stored in the map
    integer :: NNA = 0

    ! Total number of items that can initially be stored in the map.
    ! This information is needed to compute the growth of the map
    ! after several resize operations.
    integer :: NNA0 = 0

    ! Total number of resize operations
    integer :: nresize = 0
    
    ! Factor by which the map is enlarged if new storage is allocate
    real(DP) :: dfactor = 1.5_DP

#ifdef T_STORAGE
    ! Handle to the map key data
    integer :: h_Key = ST_NOHANDLE
#endif

    ! Pointer to the map key
    FEAT2_PP_TTYPE(T_TYPE), dimension(:), pointer :: p_Key => null()

    ! Handle and pointer to the balance data
    integer :: h_Kbal = ST_NOHANDLE
    integer, dimension(:), pointer :: p_Kbal => null()

    ! Handle and pointer to the parent data
    integer :: h_Kparent = ST_NOHANDLE
    integer, dimension(:), pointer :: p_Kparent => null()

    ! Handle and pointer to the children data
    integer :: h_Kchild = ST_NOHANDLE
    integer, dimension(:,:), pointer :: p_Kchild => null()
    
#ifdef D
    ! Dimension of the auxiliary data values to be stored
    integer :: isizeData = 0
#ifdef D_STORAGE
    ! Handle to the tree auxiliary data
    integer :: h_Data = ST_NOHANDLE
#endif

    ! Pointer to map auxiliary data
    FEAT2_PP_DTYPE(D_TYPE), dimension(:,:), pointer :: p_Data => null()
#endif

  end type
  
!</typeblock>

!<typeblock>

  ! An iterator for a map

  type FEAT2_PP_TEMPLATE_TD(it_map,T,D)
    private

    ! Absolute position of the current element
    integer :: ipos = MNULL

    ! Specification flag. This is a bitfield coming from an OR
    ! combination of different MAP_MSPEC_xxxx constants and specifies
    ! various details of the list iterator.
    integer(I32) :: iSpec = 0_I32

    ! Pointer to the underlying map
    type(FEAT2_PP_TEMPLATE_TD(t_map,T,D)), pointer :: p_rmap => null()

  end type

!</typeblock>

!</types>

contains

  !************************************************************************

!<subroutine>

#ifdef D
  subroutine FEAT2_PP_TEMPLATE_TD(map_create,T,D)(rmap, NNA, isizeData, dfactor)
#else
  subroutine FEAT2_PP_TEMPLATE_TD(map_create,T,D)(rmap, NNA, dfactor)
#endif

!<description>
    ! This subroutine creates a new map
!</description>

!<input>
    ! Total number of items that can be stored in the map
    integer, intent(in) :: NNA
    
#ifdef D
    ! Dimension of the auxiliary data values to be stored
    integer, intent(in) :: isizeData
#endif

    ! OPTIONAL: Factor by which the map should be enlarged if memory
    ! needs to be reallocated
    real(DP), intent(in), optional :: dfactor
!</input>

!<output>
    ! The map
    type(FEAT2_PP_TEMPLATE_TD(t_map,T,D)), intent(out) :: rmap
!</output>
!</subroutine>

    ! local variables
    integer, dimension(2) :: Isize1,Isize2
    
    ! Set factor
    if (present(dfactor)) then
      if (dfactor > 1.0_DP) rmap%dfactor=dfactor
    end if

    ! Initialise map
    rmap%NNA  = max(0,NNA)
    rmap%NNA0 = max(0,NNA)
    
    ! Allocate memory for internal tree structure
    Isize1 = (/MLEFT, MROOT/)
    Isize2 = (/MRIGHT,rmap%NNA/)
    call storage_new('map_create', 'Kchild', Isize1, Isize2,&
                     ST_INT, rmap%h_Kchild, ST_NEWBLOCK_ZERO)
    call storage_getbase_int2D(rmap%h_Kchild, rmap%p_Kchild)
    
    call storage_new('map_create', 'Kparent', MROOT, rmap%NNA,&
                     ST_INT, rmap%h_Kparent, ST_NEWBLOCK_ZERO)
    call storage_getbase_int(rmap%h_Kparent, rmap%p_Kparent)

    call storage_new('map_create', 'Kbal', rmap%NNA,&
                     ST_INT, rmap%h_Kbal, ST_NEWBLOCK_ZERO)
    call storage_getbase_int(rmap%h_Kbal, rmap%p_Kbal)

     ! Allocate memory for Key
#ifdef T_STORAGE
    call storage_new('map_create', 'Key', rmap%NNA, T_STORAGE,&
        rmap%h_Key, ST_NEWBLOCK_ZERO)
    call storage_getbase(rmap%h_Key, rmap%p_Key)
#else
    allocate(rmap%p_Key(rmap%NNA))
#endif

#ifdef D
    ! Set size of auxiliary data
    rmap%isizeData = max(0,isizeData)

    ! Allocate memory for auxiliary data
    if (rmap%isizeData > 0) then
      Isize1 = (/rmap%isizeData, rmap%NNA/)

#ifdef D_STORAGE
      call storage_new('map_create', 'Data', Isize1, D_STORAGE,&
          rmap%h_Data, ST_NEWBLOCK_NOINIT)
      call storage_getbase(rmap%h_Data, rmap%p_Data)
#else
      allocate(rmap%p_Data(Isize1(1),Isize1(2)))
#endif
    end if
#endif

    ! Clear map
    call map_clear(rmap)

  end subroutine

  !************************************************************************

!<subroutine>
  
  subroutine FEAT2_PP_TEMPLATE_TD(map_release,T,D)(rmap)

!<description>
    ! This subroutine releases an existing map
!</description>

!<inputoutput>
    ! The map
    type(FEAT2_PP_TEMPLATE_TD(t_map,T,D)), intent(inout) :: rmap
!</inputoutput>
!</subroutine>

    ! Release memory
    if (rmap%h_Kbal    .ne. ST_NOHANDLE) call storage_free(rmap%h_Kbal)
    if (rmap%h_Kchild  .ne. ST_NOHANDLE) call storage_free(rmap%h_Kchild)
    if (rmap%h_Kparent .ne. ST_NOHANDLE) call storage_free(rmap%h_Kparent)

    nullify(rmap%p_Kbal)
    nullify(rmap%p_Kchild)
    nullify(rmap%p_Kparent)

#ifdef T_STORAGE
    if (rmap%h_Key .ne. ST_NOHANDLE) call storage_free(rmap%h_Key)
    nullify(rmap%p_Key)
#else
    if (associated(rmap%p_Key)) deallocate(rmap%p_Key)
#endif

#ifdef D
    ! Reset size of auxiliary data
    rmap%isizeData = 0
#ifdef D_STORAGE
    if (rmap%h_Data   .ne. ST_NOHANDLE) call storage_free(rmap%h_Data)
    nullify(rmap%p_Data)
#else
    if (associated(rmap%p_Data))&
        deallocate(rmap%p_Data)
#endif
#endif

    ! Reset map
    rmap%NA      = 0
    rmap%NNA     = 0
    rmap%NNA0    = 0
    rmap%nresize = 0
    rmap%dfactor = 1.5_DP

  end subroutine

  !************************************************************************

!<subroutine>
  
  subroutine FEAT2_PP_TEMPLATE_TD(map_resize,T,D)(rmap, NNA)

!<description>
    ! This subroutine reallocates memory for an existing map
!</description>

!<input>
    ! New number of total items that can be stored in the map
    integer, intent(in) :: NNA
!</input>

!<inputoutput>
    ! The map
    type(FEAT2_PP_TEMPLATE_TD(t_map,T,D)), intent(inout) :: rmap
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
    NNAOld = rmap%NNA

    ! Set new size
    rmap%NNA = max(0,NNA)
    rmap%nresize = rmap%nresize+1

    ! Reallocate structures
    call storage_realloc('map_resize', rmap%NNA+1,&
                         rmap%h_Kchild, ST_NEWBLOCK_ZERO, .true.)
    call storage_getbase_int2D(rmap%h_Kchild, rmap%p_Kchild)
    
    call storage_realloc('map_resize', rmap%NNA+1,&
                         rmap%h_Kparent, ST_NEWBLOCK_ZERO, .true.)
    call storage_getbase_int(rmap%h_Kparent, rmap%p_Kparent)

    call storage_realloc('map_resize', rmap%NNA,&
                         rmap%h_Kbal, ST_NEWBLOCK_ZERO, .true.)
    call storage_getbase_int(rmap%h_Kbal, rmap%p_Kbal)

    ! Reallocate Key
#ifdef T_STORAGE
    call storage_realloc('map_resize', rmap%NNA, rmap%h_Key,&
        ST_NEWBLOCK_ZERO, .true.)
    call storage_getbase(rmap%h_Key, rmap%p_Key)
#else
    allocate(p_Key(NNAOld))
    p_Key = rmap%p_Key
    deallocate(rmap%p_Key)
    allocate(rmap%p_Key(rmap%NNA))
    rmap%p_Key(1:NNAOld) = p_Key
    deallocate(p_Key)
#endif
    
#ifdef D
    ! Reallocate auxiliary data
    if (rmap%isizeData > 0) then
#ifdef D_STORAGE
      call storage_realloc('map_resize', rmap%NNA, rmap%h_Data,&
          ST_NEWBLOCK_NOINIT, .true.)
      call storage_getbase(rmap%h_Data, rmap%p_Data)
#else
      allocate(p_Data(size(rmap%p_Data,1),NNAOld))
      p_Data = rmap%p_Data
      deallocate(rmap%p_Data)
      allocate(rmap%p_Data(size(p_Data,1),rmap%NNA))
      rmap%p_Data(:,1:NNAOld) = p_Data
      deallocate(p_Data)
#endif
    end if
#endif

  end subroutine

  !************************************************************************

!<subroutine>
  
  subroutine FEAT2_PP_TEMPLATE_TD(map_clear,T,D)(rmap)

!<description>
    ! This subroutine clears the content of the map
!</description>

!<inputoutput>
    ! The map
    type(FEAT2_PP_TEMPLATE_TD(t_map,T,D)), intent(inout) :: rmap
!</inputoutput>
!</subroutine>

    ! Initialise map
    rmap%NA      = 0
    rmap%nresize = 0
    rmap%p_Kchild(MFREE,MROOT) = 1
    rmap%p_Kparent(MROOT)      = MNULL

  end subroutine

  !************************************************************************

!<subroutine>

#ifdef D
  subroutine FEAT2_PP_TEMPLATE_TD(map_cpy1,T,D)(h_KeySrc, rmap, h_DataSrc)
#else
  subroutine FEAT2_PP_TEMPLATE_TD(map_cpy1,T,D)(h_KeySrc, rmap)
#endif

!<description>
    ! This subroutine copies the content of the given storage
    ! handle(s) to the map
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
    ! The map
    type(FEAT2_PP_TEMPLATE_TD(t_map,T,D)), intent(inout) :: rmap
!</inputoutput>
!</subroutine>

#ifdef T_STORAGE
    ! local variables
    T_TYPE, dimension(:), pointer :: p_KeySrc

#ifdef D
#ifdef D_STORAGE
    D_TYPE, dimension(:,:), pointer :: p_DataSrc

    ! Get pointers and copy the key and data values to the map
    call storage_getbase(h_KeySrc, p_KeySrc)
    call storage_getbase(h_DataSrc, p_DataSrc)
    call map_copy(p_KeySrc, rmap, p_DataSrc)
#else
    call output_line('Map does not support storage handles!',&
        OU_CLASS_ERROR,OU_MODE_STD,'map_cpy1')
    call sys_halt()
#endif

#else
    ! Get pointer and copy the key values to the map
    call storage_getbase(h_KeySrc, p_KeySrc)
    call map_copy(p_KeySrc, rmap)
#endif

#else
    call output_line('Map does not support storage handles!',&
        OU_CLASS_ERROR,OU_MODE_STD,'map_cpy1')
    call sys_halt()
#endif

  end subroutine

  !************************************************************************

!<subroutine>

#ifdef D
  subroutine FEAT2_PP_TEMPLATE_TD(map_cpy2,T,D)(KeySrc, rmap, DataSrc)
#else
  subroutine FEAT2_PP_TEMPLATE_TD(map_cpy2,T,D)(KeySrc, rmap)
#endif

!<description>
    ! This subroutine copies the content of the given array(s) to the map
!</description>

!<input>
    ! Array with key values
    FEAT2_PP_TTYPE(T_TYPE), dimension(:), intent(in) :: KeySrc

#ifdef D
    ! OPTIONAL: Array with data values
    FEAT2_PP_DTYPE(D_TYPE), dimension(:,:), intent(in), optional :: DataSrc
#endif
!</input>

!<inputoutput>
    ! The map
    type(FEAT2_PP_TEMPLATE_TD(t_map,T,D)), intent(inout) :: rmap
!</inputoutput>
!</subroutine>

    ! local variables
    type(FEAT2_PP_TEMPLATE_TD(it_map,T,D)) :: riterator
    integer :: i

#ifdef D
    if (present(DataSrc)) then
      do i = 1, size(KeySrc)
        riterator = map_insert(rmap, KeySrc(i), dataSrc(:,i))      
      end do
    else
      do i = 1, size(KeySrc)
        riterator = map_insert(rmap, KeySrc(i))
      end do
    end if
#else
    do i = 1, size(KeySrc)
      riterator = map_insert(rmap, KeySrc(i))
    end do
#endif

  end subroutine

  !************************************************************************

!<subroutine>

#ifdef D
  subroutine FEAT2_PP_TEMPLATE_TD(map_cpy3,T,D)(rmap, h_KeyDest, h_DataDest)
#else
  subroutine FEAT2_PP_TEMPLATE_TD(map_cpy3,T,D)(rmap, h_KeyDest)
#endif

!<description>
    ! This subroutine copies the content of the given map to the handle(s)
!</description>

!<input>
    ! The map
    type(FEAT2_PP_TEMPLATE_TD(t_map,T,D)), intent(in) :: rmap
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
        call storage_new('map_cpy3', 'Key',&
            rmap%NA, T_STORAGE, h_KeyDest, ST_NEWBLOCK_NOINIT)
      else
        call storage_getsize(h_KeyDest, isize)
        if (isize < rmap%NA) then
          call storage_realloc('map_cpy3',&
              rmap%NA, h_KeyDest, ST_NEWBLOCK_NOINIT, .false.)
        end if
      end if
      call storage_getbase(h_KeyDest, p_KeyDest)
      
      if (rmap%isizeData <= 0) then
        if (h_DataDest .ne. ST_NOHANDLE) call storage_free(h_DataDest)
        call output_line('Map does not provide auxiliary data!',&
            OU_CLASS_WARNING,OU_MODE_STD,'map_cpy3')

        ! Copy key values from the tree
        call map_copy(rmap, p_KeyDest)
        
        ! That`s it
        return
      end if

      ! Check data handle
      if (h_DataDest .eq. ST_NOHANDLE) then
        call storage_new('map_cpy3', 'Data',&
            (/rmap%isizeData,rmap%NA/),&
            D_STORAGE, h_DataDest, ST_NEWBLOCK_NOINIT)
      else
        call storage_getsize(h_DataDest, Isize2)

        if (Isize2(1) .ne. rmap%isizeData) then
          call output_line('Size of data array is not compatible!',&
              OU_CLASS_ERROR,OU_MODE_STD,'map_cpy3')
        end if
        
        if (Isize2(2) < rmap%NA) then
          call storage_realloc('map_cpy3',&
              rmap%NA, h_DataDest, ST_NEWBLOCK_NOINIT, .false.)
        end if
      end if
      call storage_getbase(h_DataDest, p_DataDest)
              
      ! Copy key and data values from the tree
      call map_copy(rmap, p_KeyDest, p_DataDest)
      
    else

      ! Check key handle
      if (h_KeyDest .eq. ST_NOHANDLE) then
        call storage_new('map_cpy3', 'Key',&
            rmap%NA, T_STORAGE, h_KeyDest, ST_NEWBLOCK_NOINIT)
      else
        call storage_getsize(h_KeyDest, isize)
        if (isize < rmap%NA) then
          call storage_realloc('map_cpy3',&
              rmap%NA, h_KeyDest, ST_NEWBLOCK_NOINIT, .false.)
        end if
      end if
      call storage_getbase(h_KeyDest, p_KeyDest)
      
      ! Copy key values from the tree
      call map_copy(rmap, p_KeyDest)
      
    end if
#else
    call output_line('Map does not support storage handles!',&
        OU_CLASS_ERROR,OU_MODE_STD,'map_cpy3')
    call sys_halt()
#endif

#else

    ! Check key handle
    if (h_KeyDest .eq. ST_NOHANDLE) then
      call storage_new('map_cpy3', 'Key',&
          rmap%NA, T_STORAGE, h_KeyDest, ST_NEWBLOCK_NOINIT)
    else
      call storage_getsize(h_KeyDest, isize)
      if (isize < rmap%NA) then
        call storage_realloc('map_cpy3',&
            rmap%NA, h_KeyDest, ST_NEWBLOCK_NOINIT, .false.)
      end if
    end if
    call storage_getbase(h_KeyDest, p_KeyDest)
    
    ! Copy key valus from the tree
    call map_copy(rmap, p_KeyDest)

#endif

#else
    call output_line('Map does not support storage handles!',&
        OU_CLASS_ERROR,OU_MODE_STD,'map_cpy3')
    call sys_halt()
#endif

  end subroutine

  !************************************************************************

!<subroutine>

#ifdef D
  subroutine FEAT2_PP_TEMPLATE_TD(map_cpy4,T,D)(rmap, KeyDest, DataDest)
#else
  subroutine FEAT2_PP_TEMPLATE_TD(map_cpy4,T,D)(rmap, KeyDest)
#endif

!<description>
    ! This subroutine copies the content of the given tree to the array(s)
!</description>

!<input>
    ! The map
    type(FEAT2_PP_TEMPLATE_TD(t_map,T,D)), intent(in) :: rmap
!</input>

!<inputoutput>
    ! Array with key values
    FEAT2_PP_TTYPE(T_TYPE), dimension(:), intent(inout) :: KeyDest

#ifdef D
    ! OPTIONAL: Array with data values
    FEAT2_PP_DTYPE(D_TYPE), dimension(:,:), intent(inout), optional :: DataDest
#endif
!</inputoutput>
!</subroutine>

    ! local variable
    integer j

    ! Check size of array
    if (size(KeyDest) < rmap%NA) then
      call output_line('Array for key values is too small!',&
          OU_CLASS_ERROR,OU_MODE_STD,'map_cpy4')
      call sys_halt()
    end if
  
#ifdef D

    if (present(DataDest)) then

      ! Check size of array
      if (size(DataDest) < rmap%NA) then
        call output_line('Array for data values is too small!',&
            OU_CLASS_ERROR,OU_MODE_STD,'map_cpy4')
        call sys_halt()
      end if

      j = 0 ! initialise global counter
      
      call copyKeyData(rmap%p_Kchild(MRIGHT, MROOT))

    else

      j = 0 ! initialise global counter
      
      call copyKey(rmap%p_Kchild(MRIGHT, MROOT))
      
    end if

#else

    j = 0 ! initialise global counter

    call copyKey(rmap%p_Kchild(MRIGHT, MROOT))

#endif

  contains
    
    ! Here, the real working routine follows
    
    !**************************************************************
    ! Copy the content of the key to an array (inorder traversal)
    
    recursive subroutine copyKey(i)
      integer, intent(in) :: i
      
      if (rmap%p_Kchild(MLEFT,i) .ne. MNULL)&
          call copyKey(rmap%p_Kchild(MLEFT,i))

      j = j+1
      KeyDest(j) = rmap%p_Key(i)

      if (rmap%p_Kchild(MRIGHT,i) .ne. MNULL)&
          call copyKey(rmap%p_Kchild(MRIGHT,i))

    end subroutine copyKey
    
#ifdef D
    
    !**************************************************************
    ! Copy the content of the key and data to arrays (inorder traversal)
    
    recursive subroutine copyKeyData(i)
      integer, intent(in) :: i
      
      if (rmap%p_Kchild(MLEFT,i) .ne. MNULL)&
          call copyKeyData(rmap%p_Kchild(MLEFT,i))

      j = j+1
      KeyDest(j)    = rmap%p_Key(i)
      DataDest(:,j) = rmap%p_Data(:,i)

      if (rmap%p_Kchild(MRIGHT,i) .ne. MNULL)&
          call copyKeyData(rmap%p_Kchild(MRIGHT,i))

    end subroutine copyKeyData

#endif

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine FEAT2_PP_TEMPLATE_TD(map_swap,T,D)(rmap1, rmap2)

!<description>
    ! This subroutine swaps content of two maps
!</description>
    
!<inputoutput>
    ! First map
    type(FEAT2_PP_TEMPLATE_TD(t_map,T,D)), intent(inout) :: rmap1

    ! Second map
    type(FEAT2_PP_TEMPLATE_TD(t_map,T,D)), intent(inout) :: rmap2
!</inputoutput>
!</subroutine>

    ! local variables
    type(FEAT2_PP_TEMPLATE_TD(t_map,T,D)) :: rmap
    
    ! Swap maps
    rmap  = rmap1
    rmap1 = rmap2
    rmap2 = rmap
    
    ! Reassociated pointers of first map
    call storage_getbase_int(rmap1%h_Kbal, rmap1%p_Kbal)
    call storage_getbase_int(rmap1%h_Kparent, rmap1%p_Kparent)
    call storage_getbase_int2d(rmap1%h_Kchild, rmap1%p_Kchild)

    ! Reassociated pointers of second map
    call storage_getbase_int(rmap2%h_Kbal, rmap2%p_Kbal)
    call storage_getbase_int(rmap2%h_Kparent, rmap2%p_Kparent)
    call storage_getbase_int2d(rmap2%h_Kchild, rmap2%p_Kchild)

#ifdef T_STORAGE
    ! Reassociated pointers to key values
    call storage_getbase(rmap1%h_Key, rmap1%p_Key)
    call storage_getbase(rmap2%h_Key, rmap2%p_Key)
#endif

#if defined D && defined D_STORAGE
    ! Reassociated pointers to data values
    if (rmap1%h_Data .ne. ST_NOHANDLE) &
        call storage_getbase(rmap1%h_Data, rmap1%p_Data)
    if (rmap2%h_Data .ne. ST_NOHANDLE) &
        call storage_getbase(rmap2%h_Data, rmap2%p_Data)
#endif
    
  end subroutine

  !************************************************************************

!<function>
  
  function FEAT2_PP_TEMPLATE_TD(map_begin,T,D)(rmap) result(riterator)

!<description>
    ! This function returns an iterator referring to the first
    ! element of the map
!</description>

!<input>
    ! The map
    type(FEAT2_PP_TEMPLATE_TD(t_map,T,D)), intent(in), target :: rmap
!</input>

!<result>
    ! The iterator
    type(FEAT2_PP_TEMPLATE_TD(it_map,T,D)) :: riterator
!</result>

!</function>
    
    ! Attach map to iterator
    riterator%p_rmap => rmap

    ! Initialise iterator
    riterator%ipos  = rmap%p_Kchild(MRIGHT, MROOT)
    riterator%iSpec = 0_I32

    ! Proceed to left-most child
    if (riterator%ipos .eq. MNULL) return
    do while (rmap%p_Kchild(MLEFT, riterator%ipos) .ne. MNULL)
      riterator%ipos = rmap%p_Kchild(MLEFT, riterator%ipos)
    end do

  end function

  !************************************************************************

!<function>
  
  function FEAT2_PP_TEMPLATE_TD(map_rbegin,T,D)(rmap) result(riterator)

!<description>
    ! This function returns a reverse iterator positioned to the last
    ! element of the map
!</description>

!<input>
    ! The map
    type(FEAT2_PP_TEMPLATE_TD(t_map,T,D)), intent(in), target :: rmap
!</input>

!<result>
    ! The iterator
    type(FEAT2_PP_TEMPLATE_TD(it_map,T,D)) :: riterator
!</result>

!</function>
    
    ! Attach map to iterator
    riterator%p_rmap => rmap

    ! Initialise iterator
    riterator%ipos  = rmap%p_Kchild(MRIGHT, MROOT)
    riterator%iSpec = MAP_MSPEC_REVERSE

    ! Proceed to right-most child
    if (riterator%ipos .eq. MNULL) return
    do while (rmap%p_Kchild(MRIGHT, riterator%ipos) .ne. MNULL)
      riterator%ipos = rmap%p_Kchild(MRIGHT, riterator%ipos)
    end do

  end function

  !************************************************************************

!<function>
  
  function FEAT2_PP_TEMPLATE_TD(map_end,T,D)(rmap) result(riterator)

!<description>
    ! This function returns an iterator referring to the past-the-end
    ! element of the map
!</description>

!<input>
    ! The map
    type(FEAT2_PP_TEMPLATE_TD(t_map,T,D)), intent(in), target :: rmap
!</input>

!<result>
    ! The iterator
    type(FEAT2_PP_TEMPLATE_TD(it_map,T,D)) :: riterator
!</result>
!</function>
    
    ! Attach map to iterator
    riterator%p_rmap => rmap

    ! Initialise iterator
    riterator%ipos  = MNULL
    riterator%iSpec = 0_I32

  end function

  !************************************************************************

!<function>
  
  function FEAT2_PP_TEMPLATE_TD(map_rend,T,D)(rmap) result(riterator)

!<description>
    ! This function returns a reverse iterator referring to the
    ! element right before the first element of the map
!</description>

!<input>
    ! The map
    type(FEAT2_PP_TEMPLATE_TD(t_map,T,D)), intent(in), target :: rmap
!</input>

!<result>
    ! The iterator
    type(FEAT2_PP_TEMPLATE_TD(it_map,T,D)) :: riterator
!</result>
!</function>
    
    ! Attach map to iterator
    riterator%p_rmap => rmap

    ! Initialise iterator
    riterator%ipos  = MNULL
    riterator%iSpec = MAP_MSPEC_REVERSE

  end function

  !************************************************************************

!<subroutine>
  
  subroutine FEAT2_PP_TEMPLATE_TD(map_next,T,D)(riterator)

!<description>
    ! This subroutine increments the map iterator by one.
!</description>

!<inputoutput>
    ! The iterator
    type(FEAT2_PP_TEMPLATE_TD(it_map,T,D)), intent(inout) :: riterator
!</inputoutput>
!</subroutine>
    
    ! local variable
    integer :: iparent

    if (iand(riterator%iSpec, MAP_MSPEC_REVERSE).eq.0) then

      if (.not.map_isNull(riterator)) then
        
        ! 1. If right subtree of node is not NULL, then successor lies
        ! in right subtree. Do following. Go to right subtree and return
        ! the node with minimum key value in right subtree.
        if (riterator%p_rmap%p_Kchild(MRIGHT,riterator%ipos) .ne. MNULL) then
          riterator%ipos = riterator%p_rmap%p_Kchild(MRIGHT,riterator%ipos)
          do while (riterator%p_rmap%p_Kchild(MLEFT,riterator%ipos) .ne. MNULL)
            riterator%ipos = riterator%p_rmap%p_Kchild(MLEFT,riterator%ipos)
          end do
          
          ! That`s it
          return
        end if
        
        ! If right subtree of node is NULL, then successor is one of the
        ! ancestors. Do following. Travel up using the parent pointer
        ! until you see a node which is left child of it’s parent. The
        ! parent of such a node is the successor.
        iparent = riterator%p_rmap%p_Kparent(riterator%ipos)
        do while ((iparent .ne. MNULL) .and.&
                  (riterator%ipos .eq. riterator%p_rmap%p_Kchild(MRIGHT,iparent)))
          riterator%ipos = iparent
          iparent = riterator%p_rmap%p_Kparent(iparent)
        end do
        riterator%ipos = iparent

      else

        ! Get smallest key in map
        riterator = map_begin(riterator%p_rmap)

      end if

    else

      if (.not.map_isNull(riterator)) then
        
        ! 1. If right subtree of node is not NULL, then successor lies
        ! in left subtree. Do following. Go to left subtree and return
        ! the node with maximum key value in left subtree.
        if (riterator%p_rmap%p_Kchild(MLEFT,riterator%ipos) .ne. MNULL) then
          riterator%ipos = riterator%p_rmap%p_Kchild(MLEFT,riterator%ipos)
          do while (riterator%p_rmap%p_Kchild(MRIGHT,riterator%ipos) .ne. MNULL)
            riterator%ipos = riterator%p_rmap%p_Kchild(MRIGHT,riterator%ipos)
          end do
          
          ! That`s it
          return
        end if
        
        ! If left subtree of node is NULL, then successor is one of the
        ! ancestors. Do following. Travel up using the parent pointer
        ! until you see a node which is right child of it’s parent. The
        ! parent of such a node is the successor.
        iparent = riterator%p_rmap%p_Kparent(riterator%ipos)
        do while ((iparent .ne. MNULL) .and.&
                  (riterator%ipos .eq. riterator%p_rmap%p_Kchild(MLEFT,iparent)))
          riterator%ipos = iparent
          iparent = riterator%p_rmap%p_Kparent(iparent)
        end do
        riterator%ipos = iparent

      else

        ! Get largest key in map
        riterator = map_rbegin(riterator%p_rmap)
        
      end if
     
    end if

  end subroutine

  !************************************************************************

!<subroutine>
  
  subroutine FEAT2_PP_TEMPLATE_TD(map_prior,T,D)(riterator)

!<description>
    ! This subroutine decrements the map iterator by one.
!</description>

!<inputoutput>
    ! The iterator
    type(FEAT2_PP_TEMPLATE_TD(it_map,T,D)), intent(inout) :: riterator
!</inputoutput>
!</subroutine>
    
    ! local variable
    integer :: iparent

    if (iand(riterator%iSpec, MAP_MSPEC_REVERSE).eq.0) then

      if (.not.map_isNull(riterator)) then

        ! 1. If right subtree of node is not NULL, then successor lies
        ! in left subtree. Do following. Go to left subtree and return
        ! the node with maximum key value in left subtree.
        if (riterator%p_rmap%p_Kchild(MLEFT,riterator%ipos) .ne. MNULL) then
          riterator%ipos = riterator%p_rmap%p_Kchild(MLEFT,riterator%ipos)
          do while (riterator%p_rmap%p_Kchild(MRIGHT,riterator%ipos) .ne. MNULL)
            riterator%ipos = riterator%p_rmap%p_Kchild(MRIGHT,riterator%ipos)
          end do
          
          ! That`s it
          return
        end if
        
        ! If left subtree of node is NULL, then successor is one of the
        ! ancestors. Do following. Travel up using the parent pointer
        ! until you see a node which is right child of it’s parent. The
        ! parent of such a node is the successor.
        iparent = riterator%p_rmap%p_Kparent(riterator%ipos)
        do while ((iparent .ne. MNULL) .and.&
                  (riterator%ipos .eq. riterator%p_rmap%p_Kchild(MLEFT,iparent)))
          riterator%ipos = iparent
          iparent = riterator%p_rmap%p_Kparent(iparent)
        end do
        riterator%ipos = iparent

      else

        ! Get largest key in map
        riterator = map_rbegin(riterator%p_rmap)

      end if

    else

      if (.not.map_isNull(riterator)) then
        
        ! 1. If right subtree of node is not NULL, then successor lies
        ! in right subtree. Do following. Go to right subtree and return
        ! the node with minimum key value in right subtree.
        if (riterator%p_rmap%p_Kchild(MRIGHT,riterator%ipos) .ne. MNULL) then
          riterator%ipos = riterator%p_rmap%p_Kchild(MRIGHT,riterator%ipos)
          do while (riterator%p_rmap%p_Kchild(MLEFT,riterator%ipos) .ne. MNULL)
            riterator%ipos = riterator%p_rmap%p_Kchild(MLEFT,riterator%ipos)
          end do
          
          ! That`s it
          return
        end if
        
        ! If right subtree of node is NULL, then successor is one of the
        ! ancestors. Do following. Travel up using the parent pointer
        ! until you see a node which is left child of it’s parent. The
        ! parent of such a node is the successor.
        iparent = riterator%p_rmap%p_Kparent(riterator%ipos)
        do while ((iparent .ne. MNULL) .and.&
                  (riterator%ipos .eq. riterator%p_rmap%p_Kchild(MRIGHT,iparent)))
          riterator%ipos = iparent
          iparent = riterator%p_rmap%p_Kparent(iparent)
        end do
        riterator%ipos = iparent
        
      else

        ! Get smallest key in map
        riterator = map_begin(riterator%p_rmap)

      end if

    end if
    
  end subroutine

  !************************************************************************

!<subroutine>

#ifdef D
  subroutine FEAT2_PP_TEMPLATE_TD(map_getbase_data,T,D)(rmap, rposition, p_data)

!<description>
    ! This subroutine returns pointers to the data stored at the
    ! position addressed by the iterator.
!</description>

!<input>
    ! The map
    type(FEAT2_PP_TEMPLATE_TD(t_map,T,D)), intent(in) :: rmap

    ! The iterator
    type(FEAT2_PP_TEMPLATE_TD(it_map,T,D)), intent(in) :: rposition
!</input>

!<output>
    ! Data pointer
    FEAT2_PP_DTYPE(D_TYPE), dimension(:), pointer :: p_data
!</output>
!</subroutine>

    ! Get data
    if (rmap%isizeData > 0)&
        p_data => rmap%p_Data(:,rposition%ipos)

  end subroutine
#endif

  !************************************************************************

!<function>

  function FEAT2_PP_TEMPLATE_TD(map_get,T,D)(rmap, rposition) result(key)

!<description>
    ! This functions return the key stored at the position addressed
    ! by the iterator.
!</description>

!<input>
    ! The map
    type(FEAT2_PP_TEMPLATE_TD(t_map,T,D)), intent(in) :: rmap

    ! The iterator
    type(FEAT2_PP_TEMPLATE_TD(it_map,T,D)), intent(in) :: rposition
!</input>

!<result>
    ! key value
    FEAT2_PP_TTYPE(T_TYPE) :: key
!</result>
!</function>
    
    ! Get key
    key = rmap%p_key(rposition%ipos)

  end function
  
  !************************************************************************

!<function>

#ifdef D
  function FEAT2_PP_TEMPLATE_TD(map_insert1,T,D)(rmap, key, data) result(riterator)
#else
  function FEAT2_PP_TEMPLATE_TD(map_insert1,T,D)(rmap, key) result(riterator)
#endif

!<description>
    ! This function inserts a new item into the map and returns an
    ! iterator referring to the position of the newly inserted
    ! element. If an element with the same key already exists, then
    ! the iterator refers to this position and only the auxiliary data
    ! is overwritten (if given).
!</description>

!<input>
    ! key value
    FEAT2_PP_TTYPE(T_TYPE), intent(in) :: key

#ifdef D
    ! OPTIONAL: Data
    FEAT2_PP_DTYPE(D_TYPE), dimension(:), intent(in), optional :: data
#endif
!</input>

!<inputoutput>
    ! The map
    type(FEAT2_PP_TEMPLATE_TD(t_map,T,D)), intent(inout), target :: rmap
!</inputoutput>

!<result>
    ! The new iterator
    type(FEAT2_PP_TEMPLATE_TD(it_map,T,D)) :: riterator
!</result>
!</function>

    ! local variables
    logical :: bexists
    integer :: ipos,idir,jpos

    ! Search for key in map and build-up the search path which may be
    ! required to rebalance the tree after insertion
    bexists = .false.
    jpos = MROOT
    idir = MRIGHT
    ipos = rmap%p_Kchild(idir,jpos)
    
    search: do
      if (ipos .eq. MNULL) exit search
      
      if (rmap%p_Key(ipos) .eq. key) then
        bexists = .true.
        exit search
      end if
      
      jpos = ipos
      idir = merge(MLEFT,MRIGHT,rmap%p_Key(jpos) > key)
      ipos = rmap%p_Kchild(idir,jpos)
    end do search

    ! Check if an element with the same key already exists
    if (bexists) then

      ! Get position of the corresponding leaf node
      ipos = rmap%p_Kchild(idir,jpos)

#ifdef D
      ! Update the auxiliary data (if present)
      if (present(data)) then
        if (rmap%isizeData .eq. size(data)) then
          rmap%p_Data(:,ipos) = data
        else
          call output_line('Size of auxiliary data mismatch!',&
              OU_CLASS_WARNING,OU_MODE_STD,'map_insert1')
        end if
      end if
#endif
      
    else

      ! Adjust size and position
      rmap%NA = rmap%NA+1
      ipos    = rmap%p_Kchild(MFREE,MROOT)

      ! Check if memory needs to be expanded
      if (abs(ipos) > rmap%NNA) &
          call map_resize(rmap,ceiling(rmap%dfactor*rmap%NNA))
      
      ! Compute next free position in memory
      if (ipos > 0) then
        rmap%p_Kchild(MFREE,MROOT) = ipos+1
      else
        ipos = abs(ipos)
        rmap%p_Kchild(MFREE,MROOT) = rmap%p_Kchild(MFREE,ipos)
      end if

      ! Create new leaf node and store key value
      rmap%p_Key(ipos)           = key
      rmap%p_Kbal(ipos)          = 0
      rmap%p_Kparent(ipos)       = jpos
      rmap%p_Kchild(MLEFT, ipos) = MNULL
      rmap%p_Kchild(MRIGHT,ipos) = MNULL
      
      ! Insert new node into the tree and into the search path
      rmap%p_Kchild(idir,jpos) = ipos

#ifdef D
      ! Store auxiliary data in memory (if present)
      if (present(data)) then
        if (rmap%isizeData .eq. size(data)) then
          rmap%p_Data(:,ipos) = data
        else
          call output_line('Size of auxiliary data mismatch!',&
              OU_CLASS_WARNING,OU_MODE_STD,'map_insert1')
        end if
      end if
#endif
      
      ! Invoke rebalance procedure
      if (jpos .ne. MROOT) call rebalanceAfterInsert(rmap, jpos, idir)
    end if
    
    ! Initialise the iterator referring to the current element
    riterator%p_rmap => rmap
    riterator%ipos   = ipos
    riterator%iSpec  = 0_I32

  end function

  !************************************************************************

!<function>

  function FEAT2_PP_TEMPLATE_TD(map_insert2,T,D)(rmap, rfirst, rlast) result(riterator)

!<description>
    ! This function inserts content in the range [rfirst,rlast)
    ! into the map.
!</description>

!<input>
    ! Iterator referring to the first element
    type(FEAT2_PP_TEMPLATE_TD(it_map,T,D)), intent(in) :: rfirst

    ! Iterator referring to the past-the-last element
    type(FEAT2_PP_TEMPLATE_TD(it_map,T,D)), intent(in) :: rlast
!</input>

!<inputoutput>
    ! The map
    type(FEAT2_PP_TEMPLATE_TD(t_map,T,D)), intent(inout) :: rmap
!</inputoutput>

!<result>
    ! The new iterator
    type(FEAT2_PP_TEMPLATE_TD(it_map,T,D)) :: riterator
!</result>
!</function>

    ! local variable
    type(FEAT2_PP_TEMPLATE_TD(it_map,T,D)) :: riter
    FEAT2_PP_TTYPE(T_TYPE) :: key

#ifdef D
    FEAT2_PP_DTYPE(D_TYPE), dimension(:), pointer :: p_data
#endif
    
    riter = rfirst
    do while(riter /= rlast)

      key = map_get(riter%p_rmap, riter)
      
#ifdef D
      call map_getbase_data(riter%p_rmap, riter, p_data)
      riterator = map_insert(rmap, key, p_data)
#else
      riterator = map_insert(rmap, key)
#endif
      call map_next(riter)
    end do

  end function

  !************************************************************************

!<subroutine>

  subroutine FEAT2_PP_TEMPLATE_TD(map_erase1,T,D)(rmap, rposition)

!<description>
    ! This subroutine deletes an element from the map
!</description>

!<input>
    ! The iterator
    type(FEAT2_PP_TEMPLATE_TD(it_map,T,D)), intent(in) :: rposition
!</input>

!<inputoutput>
    ! The map
    type(FEAT2_PP_TEMPLATE_TD(t_map,T,D)), intent(inout) :: rmap
!</inputoutput>
!</subroutine>
    
    ! local variables
    integer :: idir,ipos,ipred

    ! Get element to be deleted, its parent element and its direction
    ! (left/right) relative to the parent elment
    ipos  = rposition%ipos
    ipred = rmap%p_Kparent(ipos)
    idir  = merge(MLEFT, MRIGHT, rmap%p_Kchild(MLEFT,ipred) .eq. ipos)

    ! Compute new dimensions
    rmap%NA = rmap%NA-1
    
    ! Check if node to be deleted has two children.
    if ((rmap%p_Kchild(MLEFT, ipos) .ne. MNULL) .and. &
        (rmap%p_Kchild(MRIGHT,ipos) .ne. MNULL)) then

      ! Descent in the left(!) subtree with root rposition always to
      ! the right until an element IPOS without right child is found
      ! and interchange data between element IPOS and the element to
      ! be deleted which is positioned at rposition
      ipred = rposition%ipos
      idir  = MLEFT

      descend: do
        ipos  = rmap%p_Kchild(idir,ipred)

        if (rmap%p_Kchild(MRIGHT,ipos) .eq. MNULL) then
          ! Change key
          rmap%p_Key(rposition%ipos) = rmap%p_Key(ipos)
#ifdef D
          ! Change auxiliary data
          if (rmap%isizeData > 0)&
              rmap%p_Data(:,rposition%ipos) = rmap%p_Data(:,ipos)
#endif
          ! That`s it
          exit descend
        end if
        
        ! Proceed with right subtree
        ipred = ipos
        idir  = MRIGHT
      end do descend

      ! Now, element IPOS to be deleted has less than two childen
      ipos = rmap%p_Kchild(idir,ipred)
    end if
    
    if (rmap%p_Kchild(MLEFT,ipos) .eq. MNULL) then
      if (rmap%p_Kchild(MRIGHT,ipos) .eq. MNULL) then
        ! Element IPOS to be deleted is leaf
        rmap%p_Kchild(idir,ipred) = MNULL
      else
        ! Element IPOS to be deleted has right child
        rmap%p_Kchild(idir,ipred) = rmap%p_Kchild(MRIGHT,ipos)
        rmap%p_Kparent(rmap%p_Kchild(MRIGHT,ipos)) = ipred
      end if
    else
      ! Element IPOS to be deleted has right child
      rmap%p_Kchild(idir,ipred) = rmap%p_Kchild(MLEFT,ipos)
      rmap%p_Kparent(rmap%p_Kchild(MLEFT,ipos)) = ipred
    end if

    ! Mark element IPOS as deleted
    rmap%p_Kchild(MFREE,ipos)  = rmap%p_Kchild(MFREE,MROOT)
    rmap%p_Kchild(MFREE,MROOT) = -ipos
    
    ! Invoke rebalance procedure
    if (ipred .ne. MROOT) call rebalanceAfterDeletion(rmap, ipred, idir)

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine FEAT2_PP_TEMPLATE_TD(map_erase2,T,D)(rmap, rfirst, rlast)

!<description>
    ! This subroutine removes the elements [rfirst,rlast) from the map.
!</description>

!<input>
    ! Iterator referring to the first element
    type(FEAT2_PP_TEMPLATE_TD(it_map,T,D)), intent(in) :: rfirst

    ! Iterator referring to the past-the-last element
    type(FEAT2_PP_TEMPLATE_TD(it_map,T,D)), intent(in) :: rlast
!</input>

!<inputoutput>
    ! The map
    type(FEAT2_PP_TEMPLATE_TD(t_map,T,D)), intent(inout) :: rmap
!</inputoutput>
!</subroutine>

    ! local variable
    type(FEAT2_PP_TEMPLATE_TD(it_map,T,D)) :: riterator,ritTmp

    riterator = rfirst
    do while(riterator /= rlast)
      ritTmp = riterator
      call map_next(riterator)
      call map_erase(rmap, ritTmp)
    end do

  end subroutine

  !************************************************************************

!<function>

  function FEAT2_PP_TEMPLATE_TD(map_find,T,D)(rmap, key) result(riterator)

!<description>
    ! This function searches for a given key value in the map and
    ! returns an iterator referring to the element in the map if the
    ! key is present. Othewise NULL is returned.
!</description>

!<input>
    ! key value
    FEAT2_PP_TTYPE(T_TYPE), intent(in) :: key
!</input>

!<inputoutput>
    ! The map
    type(FEAT2_PP_TEMPLATE_TD(t_map,T,D)), intent(inout), target :: rmap
!</inputoutput>

!<result>
    ! The iterator
    type(FEAT2_PP_TEMPLATE_TD(it_map,T,D)) :: riterator
!</result>
!</function>

    ! Initialise iterator
    riterator%p_rmap => rmap
    riterator%ipos   = rmap%p_Kchild(MRIGHT, MROOT)
    riterator%iSpec  = 0_I32

    find: do while (.not.(map_isNull(riterator)))
      if (rmap%p_Key(riterator%ipos) .eq. key) then
        return
      elseif (rmap%p_Key(riterator%ipos) .gt. key) then
        riterator%ipos  = rmap%p_Kchild(MLEFT,riterator%ipos)
      else
        riterator%ipos  = rmap%p_Kchild(MRIGHT,riterator%ipos)
      end if
    end do find
    
    ! If we are here, then the element could not be found
    riterator%ipos  = MNULL

  end function

  !************************************************************************

!<subroutine>
  
  subroutine FEAT2_PP_TEMPLATE_TD(map_print,T,D)(rmap)

!<description>
    ! This subroutine prints the content of the map
!</description>

!<input>
    ! The map
    type(FEAT2_PP_TEMPLATE_TD(t_map,T,D)), intent(in) :: rmap
!</input>
!</subroutine>

    ! local variable
    type(FEAT2_PP_TEMPLATE_TD(it_map,T,D)) :: riterator

    riterator = map_begin(rmap)

    do while (.not.map_isNull(riterator))
      write(*,*) rmap%p_Key(riterator%ipos)
      call map_next(riterator)
    end do
    
  end subroutine

  !************************************************************************

!<subroutine>
  
  subroutine FEAT2_PP_TEMPLATE_TD(map_info,T,D)(rmap)

!<description>
    ! This subroutine outputs statistical info about the map
!</description>

!<input>
    ! The map
    type(FEAT2_PP_TEMPLATE_TD(t_map,T,D)), intent(in) :: rmap
!</input>
!</subroutine>

    call output_line('Tree statistics:')
    call output_line('----------------')
    call output_line('NA:       '//trim(sys_siL(rmap%NA,15)))
    call output_line('NNA:      '//trim(sys_siL(rmap%NNA,15)))
    call output_line('NNA0:     '//trim(sys_siL(rmap%NNA0,15)))
    call output_line('NRESIZE:  '//trim(sys_siL(rmap%NRESIZE,15)))
#ifdef T_STORAGE
    call output_line('h_Key:    '//trim(sys_siL(rmap%h_Key,15)))
#endif
    call output_line('h_Kbal:   '//trim(sys_siL(rmap%h_Kbal,15)))
    call output_line('h_Kchild: '//trim(sys_siL(rmap%h_Kchild,15)))
#if defined D && defined D_STORAGE
    call output_line('h_Data:   '//trim(sys_siL(rmap%h_Data,15)))
#endif
    call output_lbrk()
    call output_line('Current data  memory usage: '//&
        trim(sys_sdL(100*rmap%NA/real(rmap%NNA,DP),2))//'%')

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine FEAT2_PP_TEMPLATE_TD(map_duplicate,T,D)(rmap, rmapBackup)

!<description>
    ! This subroutine makes a copy of a map in memory. It
    ! does not make sense to share some information between binary
    ! trees, so each vectors is physically copied from the source tree
    ! to the destination tree.
!</description>

!<input>
    ! Source map
    type(FEAT2_PP_TEMPLATE_TD(t_map,T,D)), intent(in) :: rmap
!</input>

!<inputoutput>
    ! Destination map
    type(FEAT2_PP_TEMPLATE_TD(t_map,T,D)), intent(inout) :: rmapBackup
!</inputoutput>
!</subroutine>

    ! Release backup tree
    call map_release(rmapBackup)

    ! Copy all data
    rmapBackup = rmap

    ! Reset Handles
    rmapBackup%h_Kbal   = ST_NOHANDLE
    rmapBackup%h_Kchild = ST_NOHANDLE

    ! Copy storage blocks
    if (rmap%h_Kbal .ne. ST_NOHANDLE) then
      call storage_copy(rmap%h_Kbal, rmapBackup%h_Kbal)
      call storage_getbase_int(rmapBackup%h_Kbal, rmapBackup%p_Kbal)
    end if

    if (rmap%h_Kchild .ne. ST_NOHANDLE) then
      call storage_copy(rmap%h_Kchild, rmapBackup%h_Kchild)
      call storage_getbase_int2D(rmapBackup%h_Kchild, rmapBackup%p_Kchild)
    end if

#ifdef T_STORAGE
    rmapBackup%h_Key = ST_NOHANDLE
    if (rmap%h_Key .ne. ST_NOHANDLE) then
      call storage_copy(rmap%h_Key, rmapBackup%h_Key)
      call storage_getbase(rmapBackup%h_Key, rmapBackup%p_Key)
    endif
#else
    if (associated(rmap%p_Key)) then
      allocate(rmapBackup%p_Key(size(rmap%p_Key)))
      rmapBackup%p_Key = rmap%p_Key
    else
      nullify(rmapBackup%p_Key)
    end if
#endif

#ifdef D
#ifdef D_STORAGE
    rmapBackup%h_Data = ST_NOHANDLE
    if (rmap%h_Data .ne. ST_NOHANDLE) then
      call storage_copy(rmap%h_Data, rmapBackup%h_Data)
      call storage_getbase(rmapBackup%h_Data,&
                           rmapBackup%p_Data)
    end if
#else
    if (associated(rmap%p_Data)) then
      allocate(rmapBackup%p_Data(size(rmap%p_Data,1),&
                                   size(rmap%p_Data,2)))
      rmapBackup%p_Data = rmap%p_Data
    else
      nullify(rmapBackup%p_Data)
    end if
#endif
#endif
   
  end subroutine

  !************************************************************************

!<subroutine>

  subroutine FEAT2_PP_TEMPLATE_TD(map_restore,T,D)(rmapBackup, rmap)

!<description>
    ! This subroutine restores a map from a previous backup.
!</description>

!<input>
    ! Backup copy of the map
    type(FEAT2_PP_TEMPLATE_TD(t_map,T,D)), intent(in) :: rmapBackup
!</input>

!<inputoutput>
    ! Destination map
    type(FEAT2_PP_TEMPLATE_TD(t_map,T,D)), intent(inout) :: rmap
!</inputoutput>
!</subroutine>
    
    ! Release map
    call map_release(rmap)
    
    ! Duplicate the backup
    call map_duplicate(rmapBackup, rmap)

  end subroutine

  !************************************************************************

!<function>

  pure function FEAT2_PP_TEMPLATE_TD(map_size,T,D)(rmap) result(isize)

!<description>
    ! Returns the size of the map
!</description>

!<input>
    ! The map
    type(FEAT2_PP_TEMPLATE_TD(t_map,T,D)), intent(in) :: rmap
!</input>

!<result>
    ! Size of the map
    integer :: isize
!</result>
!</function>

    isize = rmap%NA

  end function

  !************************************************************************

!<function>

  pure function FEAT2_PP_TEMPLATE_TD(map_max_size,T,D)(rmap) result(imaxsize)

!<description>
    ! Returns the maximum size of the map
!</description>

!<input>
    ! The map
    type(FEAT2_PP_TEMPLATE_TD(t_map,T,D)), intent(in) :: rmap
!</input>

!<result>
    ! Maximum size of the map
    integer :: imaxsize
!</result>
!</function>

    imaxsize = rmap%NNA

  end function

  !************************************************************************

!<function>

  pure function FEAT2_PP_TEMPLATE_TD(map_empty,T,D)(rmap) result(bempty)

!<description>
    ! Checks if the map is empty
!</description>

!<input>
    ! The map
    type(FEAT2_PP_TEMPLATE_TD(t_map,T,D)), intent(in) :: rmap
!</input>

!<result>
    ! Flag: is true if the map is empty
    logical :: bempty
!</result>
!</function>

    bempty = (rmap%NA .eq. 0)

  end function

  !************************************************************************

!<subroutine>

  subroutine FEAT2_PP_TEMPLATE_TD(map_fassign,T,D)(rmapDest, rmapSrc)

!<description>
    ! Assigns the content of rmapSrc to rmapDest
!</description>

!<input>
    ! Source map
    type(FEAT2_PP_TEMPLATE_TD(t_map,T,D)), intent(in) :: rmapSrc
!</input>

!<output>
    ! Destination map
    type(FEAT2_PP_TEMPLATE_TD(t_map,T,D)), intent(out) :: rmapDest
!</output>
!</subroutine>

    ! Create empty map
#ifdef D
    call map_create(rmapDest, rmapSrc%NNA, rmapSrc%isizeData, rmapSrc%dfactor)
#else
    call map_create(rmapDest, rmapSrc%NNA, rmapSrc%dfactor)
#endif

    ! Set size
    rmapDest%NA = rmapSrc%NA
    
    ! Set structure
    rmapDest%p_Kbal   = rmapSrc%p_Kbal
    rmapDest%p_Kchild = rmapSrc%p_Kchild
    
#ifdef D
    rmapDest%p_Data = rmapSrc%p_Data
#endif
    
  end subroutine

  !************************************************************************

!<function>

  pure function FEAT2_PP_TEMPLATE_TD(map_isNull,T,D)(riterator) result(bisNull)

!<description>
    ! Checks if the iterator is NULL
!</description>

!<input>
    ! Iterator
    type(FEAT2_PP_TEMPLATE_TD(it_map,T,D)), intent(in) :: riterator
!</input>

!<result>
    ! True if the iterator is NULL.
    logical :: bisNull
!</result>
!</function>

    bisNull = (riterator%ipos == MNULL)

  end function

  !************************************************************************

!<function>

  pure function FEAT2_PP_TEMPLATE_TD(map_hasSpec,T,D)(riterator,iSpec) result(bhasSpec)

!<description>
    ! Checks if the iterator has the given specification flag
!</description>

!<input>
    ! Iterator
    type(FEAT2_PP_TEMPLATE_TD(it_map,T,D)), intent(in) :: riterator

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

  pure function FEAT2_PP_TEMPLATE_TD(it_map_eq,T,D)(riterator1,riterator2) result(beq)

!<description>
    ! Compare two iterators for equality
!</description>

!<input>
    ! Iterators
    type(FEAT2_PP_TEMPLATE_TD(it_map,T,D)), intent(in) :: riterator1,riterator2
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

  pure function FEAT2_PP_TEMPLATE_TD(it_map_ne,T,D)(riterator1,riterator2) result(bne)

!<description>
    ! Compare two iterators for inequality
!</description>

!<input>
    ! Iterators
    type(FEAT2_PP_TEMPLATE_TD(it_map,T,D)), intent(in) :: riterator1,riterator2
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

  pure function FEAT2_PP_TEMPLATE_TD(it_map_lt,T,D)(riterator1,riterator2) result(blt)

!<description>
    ! Checks lexicographical ordering of two iterators
!</description>

!<input>
    ! Iterators
    type(FEAT2_PP_TEMPLATE_TD(it_map,T,D)), intent(in) :: riterator1,riterator2
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

  pure function FEAT2_PP_TEMPLATE_TD(it_map_le,T,D)(riterator1,riterator2) result(ble)

!<description>
    ! Checks lexicographical ordering of two iterators
!</description>

!<input>
    ! Iterators
    type(FEAT2_PP_TEMPLATE_TD(it_map,T,D)), intent(in) :: riterator1,riterator2
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

  pure function FEAT2_PP_TEMPLATE_TD(it_map_gt,T,D)(riterator1,riterator2) result(bgt)

!<description>
    ! Checks lexicographical ordering of two iterators
!</description>

!<input>
    ! Iterators
    type(FEAT2_PP_TEMPLATE_TD(it_map,T,D)), intent(in) :: riterator1,riterator2
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

  pure function FEAT2_PP_TEMPLATE_TD(it_map_ge,T,D)(riterator1,riterator2) result(bge)

!<description>
    ! Checks lexicographical ordering of two iterators
!</description>

!<input>
    ! Iterators
    type(FEAT2_PP_TEMPLATE_TD(it_map,T,D)), intent(in) :: riterator1,riterator2
!</input>

!<result>
    ! True if the lexicographical ordering of iterator1 is greater
    ! than or equal to that of iterator2
    logical :: bge
!</result>
!</function>

    bge = (riterator1%ipos >= riterator2%ipos)

  end function
  
  !########################################################################
  !# Auxiliary routines private to this module
  !########################################################################

!<subroutine>
  
  recursive subroutine rebalanceAfterInsert(rmap, v, idir)

!<description>
    ! This subroutine rebalances the AVL tree after the insertion of
    ! node v as the idir's (left/right) child of its parent node.
    ! NOTE: This subroutine is private and must not be called by-hand.
!</description>

!<input>
    ! starting node
    integer, intent(in) :: v

    ! direction from which starting node is reached from its parent
    integer, intent(in) :: idir
!</input>

!<inputoutput>
    ! The binary search tree
    type(FEAT2_PP_TEMPLATE_TD(t_map,T,D)), intent(inout) :: rmap
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: d,r,w,x,y,z
    
    ! Which rotation needs to be performed?
    select case (idir)
      
    case (MLEFT)
      ! Node V has been reached from the left
      select case(rmap%p_Kbal(v))
        
      case (1) ! bal(v)=1
        rmap%p_Kbal(v) = 0
        
      case (0) ! bal(v)=0
        rmap%p_Kbal(v) = -1

        if (v .ne. rmap%p_Kchild(MRIGHT,MROOT)) then
          ! Get root of subtree to be rebalanced
          r = rmap%p_Kparent(v)
          d = merge(MLEFT, MRIGHT, rmap%p_Kchild(MLEFT,r) .eq. v)
          call rebalanceAfterInsert(rmap,r,d)
        end if
        
      case (-1) ! bal(v)=-1
        x = rmap%p_Kchild(MLEFT,v)
        w = rmap%p_Kchild(MRIGHT,x)
        
        select case(rmap%p_Kbal(x))
          
        case (-1,0) ! bal(x)=-1 or bal(x)=0
          ! Get root of subtree
          r = rmap%p_Kparent(v)
          d = merge(MLEFT, MRIGHT, rmap%p_Kchild(MLEFT,r) .eq. v)

          ! Single-Right-rotation
          rmap%p_Kchild(d,r)      = x
          rmap%p_Kchild(MLEFT,v)  = w
          rmap%p_Kchild(MRIGHT,x) = v

          rmap%p_Kparent(v)       = x
          rmap%p_Kparent(w)       = v
          rmap%p_Kparent(x)       = r

          rmap%p_Kbal(x)          = rmap%p_Kbal(x)+1
          rmap%p_Kbal(v)          = -rmap%p_Kbal(x)
          
        case (1) ! bal(x)=1
          ! Get root of subtree
          r = rmap%p_Kparent(v)
          d = merge(MLEFT, MRIGHT, rmap%p_Kchild(MLEFT,r) .eq. v)
          y = rmap%p_Kchild(MLEFT,w)
          z = rmap%p_Kchild(MRIGHT,w)

          ! Double-Left-Right-rotation
          rmap%p_Kchild(d,r)      = w
          rmap%p_Kchild(MLEFT,v)  = z
          rmap%p_Kchild(MRIGHT,x) = y
          rmap%p_Kchild(MLEFT,w)  = x
          rmap%p_Kchild(MRIGHT,w) = v

          rmap%p_Kparent(v)       = w
          rmap%p_Kparent(w)       = r
          rmap%p_Kparent(x)       = w
          rmap%p_Kparent(y)       = x
          rmap%p_Kparent(z)       = v

          rmap%p_Kbal(v)          = -min(0,rmap%p_Kbal(w))
          rmap%p_Kbal(x)          = -max(0,rmap%p_Kbal(w))
          rmap%p_Kbal(w)          = 0
        end select
        
      end select
      
    case (MRIGHT)
      ! Node V has been reached from the right
      select case(rmap%p_Kbal(v))
        
      case (-1) ! bal(v)=-1
        rmap%p_Kbal(v) = 0
        
      case (0) ! bal(v)=0
        rmap%p_Kbal(v) = 1
        if (v .ne. rmap%p_Kchild(MRIGHT,MROOT)) then
          ! Get root of subtree to be rebalanced
          r = rmap%p_Kparent(v)
          d = merge(MLEFT, MRIGHT, rmap%p_Kchild(MLEFT,r) .eq. v)
          call rebalanceAfterInsert(rmap,r,d)
        end if
        
      case (1) ! bal(v)=1
        x = rmap%p_Kchild(MRIGHT,v)
        w = rmap%p_Kchild(MLEFT,x)
        
        select case(rmap%p_Kbal(x))
          
        case (0,1) ! bal(x)=0 or bal(x)=1
          ! Get root of subtree
          r = rmap%p_Kparent(v)
          d = merge(MLEFT, MRIGHT, rmap%p_Kchild(MLEFT,r) .eq. v)

          ! Single-Left-rotation
          rmap%p_Kchild(d,r)      = x
          rmap%p_Kchild(MRIGHT,v) = w
          rmap%p_Kchild(MLEFT,x)  = v

          rmap%p_Kparent(v)       = x
          rmap%p_Kparent(w)       = v
          rmap%p_Kparent(x)       = r

          rmap%p_Kbal(x)          = rmap%p_Kbal(x)-1
          rmap%p_Kbal(v)          = -rmap%p_Kbal(x)
          
        case (-1) ! bal(x)=-1
          ! Get root of subtree
          r = rmap%p_Kparent(v)
          d = merge(MLEFT, MRIGHT, rmap%p_Kchild(MLEFT,r) .eq. v)
          y = rmap%p_Kchild(MLEFT,w)
          z = rmap%p_Kchild(MRIGHT,w)

          ! Double-Right-Left-rotation
          rmap%p_Kchild(d,r)      = w
          rmap%p_Kchild(MRIGHT,v) = y
          rmap%p_Kchild(MLEFT,x)  = z
          rmap%p_Kchild(MLEFT,w)  = v
          rmap%p_Kchild(MRIGHT,w) = x

          rmap%p_Kparent(v)       = w
          rmap%p_Kparent(w)       = r
          rmap%p_Kparent(x)       = w
          rmap%p_Kparent(y)       = v
          rmap%p_Kparent(z)       = x

          rmap%p_Kbal(v)          = -max(0,rmap%p_Kbal(w))
          rmap%p_Kbal(x)          = -min(0,rmap%p_Kbal(w))
          rmap%p_Kbal(w)          = 0
        end select
        
      end select
      
    end select

  end subroutine

  !************************************************************************

!<subroutine>
  
  recursive subroutine rebalanceAfterDeletion(rmap, v, idir)

!<description>
    ! This subroutine rebalances the AVL tree after the deletion of the node
    ! node v as the idir's (left/right) child of its parent node.

    ! NOTE: This subroutine is private and must not be called by-hand.
!</description>

!<input>
    ! starting node
    integer, intent(in) :: v

    ! direction from which starting node is reached from its parent
    integer, intent(in) :: idir
!</input>

!<inputoutput>
    ! The map
    type(FEAT2_PP_TEMPLATE_TD(t_map,T,D)), intent(inout) :: rmap
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: d,r,w,x,xbal,y,z
    
    ! Which rotation needs to be performed?
    select case (idir)
      
    case (MLEFT)
      ! Node V has been reached from the left
      select case (rmap%p_Kbal(v))
        
      case (0) ! bal(v)=0
        rmap%p_Kbal(v) = 1
        
      case (-1) ! bal(v)=-1
        rmap%p_Kbal(v) = 0

        if (v .ne. rmap%p_Kchild(MRIGHT,MROOT)) then
          ! Get root of subtree to be rebalanced
          r = rmap%p_Kparent(v)
          d = merge(MLEFT, MRIGHT, rmap%p_Kchild(MLEFT,r) .eq. v)
          call rebalanceAfterDeletion(rmap,r,d)
        end if
        
      case (1) ! bal(v)=1
        x = rmap%p_Kchild(MRIGHT,v)
        w = rmap%p_Kchild(MLEFT,x)

        ! Get root of subtree
        r = rmap%p_Kparent(v)
        d = merge(MLEFT, MRIGHT, rmap%p_Kchild(MLEFT,r) .eq. v)

        select case(rmap%p_Kbal(x))
          
        case (-1) ! bal(x)=-1
          y = rmap%p_Kchild(MLEFT,w)
          z = rmap%p_Kchild(MRIGHT,w)
          
          ! Double-Right-Left-rotation
          rmap%p_Kchild(MRIGHT,v) = y
          rmap%p_Kchild(MLEFT,x)  = z
          rmap%p_Kchild(MLEFT,w)  = v
          rmap%p_Kchild(MRIGHT,w) = x

          rmap%p_Kparent(v)       = w
          rmap%p_Kparent(w)       = r
          rmap%p_Kparent(x)       = w
          rmap%p_Kparent(y)       = v
          rmap%p_Kparent(z)       = x
          
          rmap%p_Kbal(v)          = -max(0,rmap%p_Kbal(w))
          rmap%p_Kbal(x)          = -min(0,rmap%p_Kbal(w))
          rmap%p_Kbal(w)          = 0
          
          if (v .eq. rmap%p_Kchild(MRIGHT,MROOT)) then
            rmap%p_Kchild(MRIGHT,MROOT) = w
          else
            rmap%p_Kchild(d,r) = w
            call rebalanceAfterDeletion(rmap,r,d)
          end if
          
        case DEFAULT ! bal(x)=0 or bal(x)=1

          ! Single-Left-rotation
          rmap%p_Kchild(MRIGHT,v) = w
          rmap%p_Kchild(MLEFT,x)  = v

          rmap%p_Kparent(v)       = x
          rmap%p_Kparent(w)       = v
          rmap%p_Kparent(x)       = r

          xbal                    = rmap%p_Kbal(x)
          rmap%p_Kbal(x)          = rmap%p_Kbal(x)-1
          rmap%p_Kbal(v)          = -rmap%p_Kbal(x)
          
          if (v .eq. rmap%p_Kchild(MRIGHT,MROOT)) then
            rmap%p_Kchild(MRIGHT,MROOT) = x
          else
            rmap%p_Kchild(d,r) = x
            if (xbal .eq. 1) call rebalanceAfterDeletion(rmap,r,d)
          end if
          
        end select
        
      end select
      
    case (MRIGHT)
      ! Node V has been reached from the right
      select case (rmap%p_Kbal(v))
        
      case (0) ! bal(v)=0
        rmap%p_Kbal(v) = -1
        
      case (1) ! bal(v)=1
        rmap%p_Kbal(v) = 0

        if (v .ne. rmap%p_Kchild(MRIGHT,MROOT)) then
          ! Get root of subtree to be rebalanced
          r = rmap%p_Kparent(v)
          d = merge(MLEFT, MRIGHT, rmap%p_Kchild(MLEFT,r) .eq. v)
          call rebalanceAfterDeletion(rmap,r,d)
        end if
        
      case (-1) ! bal(v)=-1
        x = rmap%p_Kchild(MLEFT,v)
        w = rmap%p_Kchild(MRIGHT,x)

        ! Get root of subtree
        r = rmap%p_Kparent(v)
        d = merge(MLEFT, MRIGHT, rmap%p_Kchild(MLEFT,r) .eq. v)
        
        select case(rmap%p_Kbal(x))
          
        case (1) ! bal(x)=1
          y = rmap%p_Kchild(MLEFT,w)
          z = rmap%p_Kchild(MRIGHT,w)

          ! Double-Left-Right-rotation
          rmap%p_Kchild(MLEFT,v)  = z
          rmap%p_Kchild(MRIGHT,x) = y
          rmap%p_Kchild(MLEFT,w)  = x
          rmap%p_Kchild(MRIGHT,w) = v

          rmap%p_Kparent(v)       = w
          rmap%p_Kparent(w)       = r
          rmap%p_Kparent(x)       = w
          rmap%p_Kparent(y)       = x
          rmap%p_Kparent(z)       = v

          rmap%p_Kbal(v)          = -min(0,rmap%p_Kbal(w))
          rmap%p_Kbal(x)          = -max(0,rmap%p_Kbal(w))
          rmap%p_Kbal(w)          = 0
          
          if (v .eq. rmap%p_Kchild(MRIGHT,MROOT)) then
            rmap%p_Kchild(MRIGHT,MROOT) = w
          else
            rmap%p_Kchild(d,r) = w
            call rebalanceAfterDeletion(rmap,r,d)
          end if
          
        case DEFAULT ! bal(x)=0 or bal(x)=-1
          ! Single-Right-rotation
          rmap%p_Kchild(MLEFT,v)  = w
          rmap%p_Kchild(MRIGHT,x) = v

          rmap%p_Kparent(v)       = x
          rmap%p_Kparent(w)       = v
          rmap%p_Kparent(x)       = r

          xbal                    = rmap%p_Kbal(x)
          rmap%p_Kbal(x)          = rmap%p_Kbal(x)+1
          rmap%p_Kbal(v)          = -rmap%p_Kbal(x)
          
          if (v .eq. rmap%p_Kchild(MRIGHT,MROOT)) then
            rmap%p_Kchild(MRIGHT,MROOT) = x
          else
            rmap%p_Kchild(d,r) = x
            if (xbal .eq. -1) call rebalanceAfterDeletion(rmap,r,d)
          end if
          
        end select
        
      end select
      
    end select

  end subroutine
  
  !************************************************************************

!<function>

  function checkConsistency(rmap) result(bisConsistent)

!<description>
    ! This function checks internal consistency such as correct parent
    ! child relations of the given map
!</description>

!<input>
    ! The map
    type(FEAT2_PP_TEMPLATE_TD(t_map,T,D)), intent(in) :: rmap
!</input>

!<result>
    ! True if the map is internally consistent
    logical :: bisConsistent
!</result>
!</function>

    ! Initialisation
    bisConsistent = .true.
    if (rmap%p_Kchild(MRIGHT,MROOT) .ne. MNULL)&
        bisConsistent = consistency(rmap%p_Kchild(MRIGHT,MROOT))
    
  contains

    recursive function consistency(ipos) result(bisConsistent)
      integer, intent(in) :: ipos
      logical :: bisConsistent

      ! Initialisation
      bisConsistent = .true.
      
      ! Left child?
      if (rmap%p_Kchild(MLEFT,ipos) .ne. MNULL) then
        
        ! Check consistency between left child and parent ndoe
        bisConsistent = (bisConsistent .and.&
            rmap%p_Kparent(rmap%p_Kchild(MLEFT,ipos)) .eq. ipos)

        ! Proceed with left subtree
        bisConsistent = (bisConsistent .and.&
            consistency(rmap%p_Kchild(MLEFT,ipos)))
      end if
      
      ! Right child?
      if (rmap%p_Kchild(MRIGHT,ipos) .ne. MNULL) then
        
        ! Check consistency between right child and parent ndoe
        bisConsistent = (bisConsistent .and.&
            rmap%p_Kparent(rmap%p_Kchild(MRIGHT,ipos)) .eq. ipos)

        ! Proceed with right subtree
        bisConsistent = (bisConsistent .and.&
            consistency(rmap%p_Kchild(MRIGHT,ipos)))
      end if
          
    end function

  end function

  !************************************************************************

!<function>
  
  elemental function log2(i)

!<description>
    ! Compute logarithm of inter i
!</description>

!<input>
    ! integer value
    integer, intent(in) :: i
!</input>

!<result>
    ! log_2(i)
    real(DP) :: log2
!</result>
!</function>
        
    log2 = log(real(i,DP))/log(2._DP)

  end function

  !************************************************************************

!<function>

  pure function getHeight(rmap, rposition) result(h)

!<description>
    ! This function computes the height of the subtree with root rposition
!</description>

!<input>
    ! The map
    type(FEAT2_PP_TEMPLATE_TD(t_map,T,D)), intent(in) :: rmap

    ! The iterator
    type(FEAT2_PP_TEMPLATE_TD(it_map,T,D)), intent(in) :: rposition
!</input>

!<result>
    integer :: h
!</result>
!</function>

    h = height(rposition%ipos)

  contains

    pure recursive function height(ipos) result(h)
      integer, intent(in) :: ipos
      integer :: h,hl,hr
      
      if (rmap%p_Kchild(MLEFT,ipos) .ne. MNULL) then
        hl = height(rmap%p_Kchild(MLEFT,ipos))
      else
        hl = 0
      end if
      
      if (rmap%p_Kchild(MRIGHT,ipos) .ne. MNULL) then
        hr = height(rmap%p_Kchild(MRIGHT,ipos))
      else
        hr = 0
      end if
      
      h = max(hl,hr)+1
    end function height

  end function

#endif
