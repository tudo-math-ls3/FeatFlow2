!-*- mode: f90; -*-

#ifndef _LIST_H_
#define _LIST_H_

!##############################################################################
!# ****************************************************************************
!# <name> TEMPLATE(list,T) </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module header file implements a linked list for T_TYPE data
!#
!#  1.) list_create
!#      -> Create an empty list
!#
!#  2.) list_release
!#      -> Release an existing list
!#
!#  3.) list_resize
!#      -> Reallocate memory for an existing list
!#
!#  4.) list_clear
!#      -> Remove all content from list
!#
!#  5.) list_copyTo
!#      -> Copy data to a list
!#
!#  6.) list_copyFrom
!#      -> Copy data from a list
!#
!#  7.) list_swap
!#      -> Swap two lists
!#
!#  8.) list_first
!#      -> Get the position of the first item in list
!#
!#  9.) list_last
!#      -> Get the position of the last item in list
!#
!# 10.) list_next
!#      -> Get the position of the next item in list
!#
!# 11.) list_prev
!#      -> Get the position of the previous item in list
!#
!# 12.) list_get
!#      -> Get key and value at given position in list
!#
!# 13.) list_prepend
!#      -> Prepend data to list
!#
!# 14.) list_append
!#      -> Append data to list
!#
!# 15.) list_insert
!#      -> Insert data into list
!#
!# 16.) list_delete
!#      -> Delete data from list
!#
!# 17.) list_search
!#      -> Search for data in list
!#
!# 18.) list_print
!#      -> Print content of list
!#
!# </purpose>
!##############################################################################

#include "template.h"

module TEMPLATE_D(list,T,D)

!$use omp_lib
  use fsystem
  use genoutput
  use listbase
  use storage

#ifndef T_STORAGE
  __external_use__(T_MODULE)
#endif

#ifdef D
#ifndef D_STORAGE
  __external_use__(D_MODULE)
#endif
#endif
  
  implicit none

  private
  public :: TEMPLATE_D(t_list,T,D)
  public :: list_create
  public :: list_release
  public :: list_resize
  public :: list_clear
  public :: list_copyTo
  public :: list_copyFrom
  public :: list_swap
  public :: list_first
  public :: list_last
  public :: list_next
  public :: list_prev
  public :: list_get
  public :: list_prepend
  public :: list_append
  public :: list_insert
  public :: list_delete
  public :: list_search
  public :: list_print

  public :: LIST_UNORDERED, LIST_INCREASING, LIST_DECREASING, LIST_CSR7
  public :: LIST_SINGLELINKED, LIST_DOUBLELINKED
  public :: LIST_NOT_FOUND, LIST_FOUND, LNULL

  interface list_create
    module procedure TEMPLATE_D(list_create,T,D)
  end interface

  interface list_release
    module procedure TEMPLATE_D(list_release,T,D)
  end interface

  interface list_resize
    module procedure TEMPLATE_D(list_resize,T,D)
  end interface

  interface list_clear
    module procedure TEMPLATE_D(list_clear,T,D)
  end interface

  interface list_copyTo
    module procedure TEMPLATE_D(list_cpyH2L,T,D)
    module procedure TEMPLATE_D(list_cpyA2L,T,D)
  end interface

  interface list_copyFrom
    module procedure TEMPLATE_D(list_cpyL2H,T,D)
    module procedure TEMPLATE_D(list_cpyL2A,T,D)
  end interface

  interface list_swap
    module procedure TEMPLATE_D(list_swap,T,D)
  end interface

  interface list_first
    module procedure TEMPLATE_D(list_first,T,D)
  end interface
  
  interface list_last
    module procedure TEMPLATE_D(list_last,T,D)
  end interface

  interface list_next
    module procedure TEMPLATE_D(list_next,T,D)
  end interface

  interface list_prev
    module procedure TEMPLATE_D(list_prev,T,D)
  end interface

  interface list_get
    module procedure TEMPLATE_D(list_get,T,D)
  end interface

  interface list_prepend
    module procedure TEMPLATE_D(list_prepend,T,D)
  end interface

  interface list_append
    module procedure TEMPLATE_D(list_append,T,D)
  end interface

  interface list_insert
    module procedure TEMPLATE_D(list_insert,T,D)
  end interface

  interface list_delete
    module procedure TEMPLATE_D(list_delete,T,D)
  end interface

  interface list_search
    module procedure TEMPLATE_D(list_search,T,D)
  end interface

  interface list_print
    module procedure TEMPLATE_D(list_print,T,D)
  end interface

  !************************************************************************

!<types>
!<typeblock>

  type TEMPLATE_D(t_list,T,D)
    private

    ! Type-Tag: Single-linked, Double-linked
    integer :: clinkType = LIST_UNORDERED

    ! Type of list ordering
    integer :: cordering = LIST_UNORDERED
    
    ! Position of the last item inserted into the list
    integer :: item

    ! Number of items that are currently stored in the list
    integer :: NA = 0

    ! Total number of items that can be stored in the list
    integer :: NNA = 0

    ! Total number of resize operations
    integer :: NRESIZE = 0

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

    ! Pointer to list key
    DTYPE(D_TYPE), dimension(:,:), pointer :: p_Data => null()
#endif
    
  end type
  
!</typeblock>
!</types>

contains

  !************************************************************************
  
!<subroutine>

#ifdef D
  subroutine TEMPLATE_D(list_create,T,D)(rlist, nna, isizeData,&
      cordering, clinkType, dfactor)
#else
  subroutine TEMPLATE_D(list_create,T,D)(rlist, nna,&
      cordering, clinkType, dfactor)
#endif

!<description>
    ! This subroutine creates a new list
!</description>

!<input>
    ! Total number of items that can be stored in list
    integer, intent(in) :: nna

#ifdef D
    ! Dimension of the auxiliary data values to be stored
    integer, intent(in) :: isizeData
#endif

    ! OPTIONAL: Format-tag. Type of list ordering
    integer, intent(in), optional :: cordering
    
    ! OPTIONAL: Type of linking (single/double). If not specified
    ! then a single-linked list is generated
    integer, intent(in), optional :: clinkType

    ! OPTIONAL: Factor by which the list should be enlarged if memory
    ! needs to be reallocated
    real(DP), intent(in), optional :: dfactor
!</input>

!<output>
    ! The linked list
    type(TEMPLATE_D(t_list,T,D)), intent(out) :: rlist
!</output>
!</subroutine>

    ! local variables
    integer, dimension(2) :: Isize

    ! Set factor
    if (present(dfactor)) then
      if (dfactor > 1.0_DP) rlist%dfactor = dfactor
    end if
    
    ! Set ordering
    if (present(cordering)) then
      rlist%cordering = cordering
    end if

    ! Initialise list
    rlist%item      = LHEAD
    rlist%na        = 0
    rlist%nna       = max(0,nna)
    rlist%clinkType = LIST_SINGLELINKED
    if (present(clinkType)) rlist%clinkType = clinkType
    
    ! Allocate memory and associate pointers
    call storage_new('list_create', 'Knext', LHEAD, rlist%nna, ST_INT,&
        rlist%h_Knext, ST_NEWBLOCK_NOINIT)
    call storage_getbase_int(rlist%h_Knext, rlist%p_Knext)
    
    ! Double-linked list?
    if (rlist%clinkType .eq. LIST_DOUBLELINKED) then
      call storage_new('list_create','Kprev', rlist%nna, ST_INT,&
          rlist%h_Kprev, ST_NEWBLOCK_NOINIT)
      call storage_getbase_int(rlist%h_Kprev, rlist%p_Kprev)
    end if

    ! Allocate memory for Key
#ifdef T_STORAGE
    call storage_new('list_create', 'Key', rlist%nna, T_STORAGE,&
        rlist%h_Key, ST_NEWBLOCK_NOINIT)
    call storage_getbase(rlist%h_Key, rlist%p_Key)
#else
    allocate(rlist%p_Key(rlist%nna))
#endif

    ! Initialise list structure
    rlist%p_Knext(LFREE) = 1
    rlist%p_Knext(LHEAD) = LNULL
    rlist%p_Knext(LTAIL) = LNULL

#ifdef D
    ! Set size of auxiliary data
    rlist%isizeData = max(0,isizeData)

    ! Allocate memory for auxiliary data
    if (rlist%isizeData > 0) then
      Isize = (/rlist%isizeData, rlist%nna/)

#ifdef D_STORAGE
      call storage_new('list_create', 'Data', Isize, D_STORAGE,&
          rlist%h_Data, ST_NEWBLOCK_NOINIT)
      call storage_getbase(rlist%h_Data, rlist%p_Data)
#else
      allocate(rlist%p_Data(Isize(1),Isize(2)))
#endif
    end if
#endif

  end subroutine

  !************************************************************************

!<subroutine>
  
  subroutine TEMPLATE_D(list_release,T,D)(rlist)

!<description>
    ! This subroutine releases an existing list
!</description>

!<inputoutput>
    ! The linked list
    type(TEMPLATE_D(t_list,T,D)), intent(inout) :: rlist
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
    deallocate(rlist%p_Key)
#endif

#ifdef D
    ! Reset size of auxiliary data
    rlist%isizeData = 0
#ifdef D_STORAGE
    if (rlist%h_Data  .ne. ST_NOHANDLE) call storage_free(rlist%h_Data)
    nullify(rlist%p_Data)
#else
    deallocate(rlist%p_Data)
#endif
#endif
    
    ! Reset linked list
    rlist%clinkType   = LIST_UNORDERED
    rlist%cordering   = LIST_UNORDERED
    rlist%item        = 0
    rlist%na          = 0
    rlist%nna         = 0
    rlist%dfactor     = 1.5_DP

  end subroutine

  !************************************************************************

!<subroutine>
  
  subroutine TEMPLATE_D(list_resize,T,D)(rlist, nna)

!<description>
    ! This subroutine reallocates memory for an existing list
!</description>

!<input>
    ! New number of total items that can be stored in the list
    integer, intent(in) :: nna
!</input>

!<inputoutput>
    ! The linked list
    type(TEMPLATE_D(t_list,T,D)), intent(inout) :: rlist
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: nnaOld
    
#ifndef T_STORAGE
    type(T_TYPE), dimension(:), pointer :: p_Key
#endif

#ifdef D
#ifndef D_STORAGE
    type(D_TYPE), dimension(:,:), pointer :: p_Data
#endif
#endif

    ! Save old size
    nnaOld = rlist%nna

    ! Set new size
    rlist%nna = nna

    ! Reallocate structures
    call storage_realloc('list_resize', rlist%nna+3, rlist%h_Knext,&
        ST_NEWBLOCK_NOINIT, .true.)
    call storage_getbase_int(rlist%h_Knext, rlist%p_Knext)

    if (rlist%clinkType .eq. LIST_DOUBLELINKED) then
      call storage_realloc('list_resize', rlist%nna, rlist%h_Kprev,&
          ST_NEWBLOCK_NOINIT, .true.)
      call storage_getbase_int(rlist%h_Kprev, rlist%p_Kprev)
    end if

    ! Reallocate Key
#ifdef T_STORAGE
    call storage_realloc('list_resize', rlist%nna, rlist%h_Key,&
        ST_NEWBLOCK_NOINIT, .true.)
    call storage_getbase(rlist%h_Key, rlist%p_Key)
#else
    allocate(p_Key(nnaOld))
    p_Key = rlist%p_Key
    deallocate(rlist%p_Key)
    allocate(rlist%p_Key(rlist%nna))
    rlist%p_Key(1:nnaOld) = p_Key
    deallocate(p_Key)
#endif
    
#ifdef D
    ! Reallocate auxiliary data
    if (rlist%isizeData > 0) then
#ifdef D_STORAGE
      call storage_realloc('list_resize', rlist%nna, rlist%h_Data,&
          ST_NEWBLOCK_NOINIT, .true.)
      call storage_getbase(rlist%h_Data, rlist%p_Data)
#else
      allocate(p_Data(size(rlist%p_Data,1),nnaOld))
      p_Data = rlist%p_Data
      deallocate(rlist%p_Data)
      allocate(rlist%p_Data(size(p_Data,1),rlist%nna))
      rlist%p_Data(:,1:nnaOld) = p_Data
      deallocate(p_Data)
#endif
    end if
#endif

  end subroutine

  !************************************************************************

!<subroutine>
  
  subroutine TEMPLATE_D(list_clear,T,D)(rlist)

!<description>
    ! This subroutine clears the content of the list
!</description>

!<inputoutput>
    ! The linked list
    type(TEMPLATE_D(t_list,T,D)), intent(inout) :: rlist
!</inputoutput>
!</subroutine>

    ! Initialise list structure
    rlist%item = LHEAD
    rlist%na   = 0
    rlist%p_Knext(LFREE) = 1
    rlist%p_Knext(LHEAD) = LNULL
    rlist%p_Knext(LTAIL) = LNULL

  end subroutine

  !************************************************************************

!<subroutine>

#ifdef D
  subroutine TEMPLATE_D(list_cpyH2L,T,D)(h_Key, rlist, h_Data)
#else
  subroutine TEMPLATE_D(list_cpyH2L,T,D)(h_Key, rlist)
#endif

!<description>
    ! This subroutine copies the content of the given storage
    ! handle(s) to the list
!</description>

!<input>
    ! Handle to the key values
    integer, intent(in) :: h_Key

#ifdef D
    ! OPTIONAL: Handle to the data values
    integer, intent(in), optional :: h_Data
#endif
!</input>

!<inputoutput>
    ! The linked list
    type(TEMPLATE_D(t_list,T,D)), intent(inout) :: rlist
!</inputoutput>
!</subroutine>

#ifdef T_STORAGE
    ! local variables
    T_TYPE, dimension(:), pointer :: p_Key

#ifdef D
#ifdef D_STORAGE
    D_TYPE, dimension(:,:), pointer :: p_Data

    ! Get pointers and copy the key,value pairs to the list
    call storage_getbase(h_Key, p_Key)
    call storage_getbase(h_Data, p_Data)
    call list_copyTo(p_Key, rlist, p_Data)
#else
    call output_line('List does not support storage handles!',&
        OU_CLASS_ERROR,OU_MODE_STD,'list_cpyH2L')
    call sys_halt()
#endif

#else
    ! Get pointer and copy the keys to the list
    call storage_getbase(h_Key, p_Key)
    call list_copyTo(p_Key, rlist)
#endif

#else
    call output_line('List does not support storage handles!',&
        OU_CLASS_ERROR,OU_MODE_STD,'list_cpyH2L')
    call sys_halt()
#endif

  end subroutine

  !************************************************************************

!<subroutine>

#ifdef D
  subroutine TEMPLATE_D(list_cpyA2L,T,D)(Key, rlist, Data)
#else
  subroutine TEMPLATE_D(list_cpyA2L,T,D)(Key, rlist)
#endif

!<description>
    ! This subroutine copies the content of the given array(s) to the list
!</description>

!<input>
    ! Array with key values
    TTYPE(T_TYPE), dimension(:), intent(in) :: Key

#ifdef D
    ! OPTIONAL: Array with data values
    DTYPE(D_TYPE), dimension(:,:), intent(in), optional :: Data
#endif
!</input>

!<inputoutput>
    ! The linked list
    type(TEMPLATE_D(t_list,T,D)), intent(inout) :: rlist
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: ipos,kpos    

#ifdef D
    if (present(Data)) then
      do ipos = 1, size(Key)
        call list_append(rlist, Key(ipos), kpos, Data(:,ipos))      
      end do
    else
      do ipos = 1, size(Key)
        call list_append(rlist, Key(ipos), kpos)
      end do
    end if
#else
    do ipos = 1, size(Key)
      call list_append(rlist, Key(ipos), kpos)
    end do
#endif

  end subroutine

  !************************************************************************

!<subroutine>

#ifdef D
  subroutine TEMPLATE_D(list_cpyL2H,T,D)(rlist, h_Key, h_Data)
#else
  subroutine TEMPLATE_D(list_cpyL2H,T,D)(rlist, h_Key)
#endif

!<description>
    ! This subroutine copies the content of the given list to the handle(s)
!</description>

!<input>
    ! The linked list
    type(TEMPLATE_D(t_list,T,D)), intent(in) :: rlist
!</input>

!<inputoutput>
    ! Handle to the key values
    integer, intent(inout) :: h_Key

#ifdef D
    ! OPTIONAL: Handle to the data values
    integer, intent(inout), optional :: h_Data
#endif
!</inputoutput>
!</subroutine>

#ifdef T_STORAGE
    ! local variable
    T_TYPE, dimension(:), pointer :: p_Key

#ifdef D
#ifdef D_STORAGE
    D_TYPE, dimension(:,:), pointer :: p_Data

    if (present(h_Data)) then
    
      ! Free memory?
      if (h_Key  .ne. ST_NOHANDLE) call storage_free(h_Key)
      if (h_Data .ne. ST_NOHANDLE) call storage_free(h_Data)

      ! Allocate new storage
      call storage_new('list_cpyL2H', 'Key', rlist%na,&
          T_STORAGE, h_Key, ST_NEWBLOCK_NOINIT)
      call storage_getbase(h_Key, p_Key)

      if (rlist%isizeData <= 0) then
        call output_line('List does not provide auxiliary data!',&
            OU_CLASS_ERROR,OU_MODE_STD,'list_cpyL2H')

        ! Copy key values from list
        call list_copyFrom(rlist, p_Key)
        
        ! That`s it
        return
      end if
        
      ! Allocate new storage
      call storage_new('list_cpyL2H', 'Data', (/rlist%isizeData,rlist%na/),&
          D_STORAGE, h_Key, ST_NEWBLOCK_NOINIT)
      call storage_getbase(h_Data, p_Data)
      
      ! Copy key and values pairs from list
      call list_copyFrom(rlist, p_Key, p_Data)
      
    else

      ! Free memory?
      if (h_Key .ne. ST_NOHANDLE) call storage_free(h_Key)
      
      ! Allocate new storage
      call storage_new('list_cpyL2H', 'Key', rlist%na,&
          T_STORAGE, h_Key, ST_NEWBLOCK_NOINIT)
      call storage_getbase(h_Key, p_Key)
      
      ! Copy list
      call list_copyFrom(rlist, p_Key)
      
    end if
#else
    call output_line('List does not support storage handles!',&
        OU_CLASS_ERROR,OU_MODE_STD,'list_cpyL2H')
    call sys_halt()
#endif

#else

    ! Free memory?
    if (h_Key .ne. ST_NOHANDLE) call storage_free(h_Key)
    
    ! Allocate new storage
    call storage_new('list_cpyL2H', 'Key', rlist%na,&
        T_STORAGE, h_Key, ST_NEWBLOCK_NOINIT)
    call storage_getbase(h_Key, p_Key)
    
    ! Copy list
    call list_copyFrom(rlist, p_Key)

#endif

#else
    call output_line('List does not support storage handles!',&
        OU_CLASS_ERROR,OU_MODE_STD,'list_cpyL2H')
    call sys_halt()
#endif

  end subroutine
  
  !************************************************************************

!<subroutine>

#ifdef D
  subroutine TEMPLATE_D(list_cpyL2A,T,D)(rlist, Key, Data)
#else
  subroutine TEMPLATE_D(list_cpyL2A,T,D)(rlist, Key)
#endif

!<description>
    ! This subroutine copies the content of the given list to the array(s)
!</description>

!<input>
    ! The linked list
    type(TEMPLATE_D(t_list,T,D)), intent(in) :: rlist
!</input>

!<inputoutput>
    ! Array with key values
    TTYPE(T_TYPE), dimension(:), intent(inout) :: Key

#ifdef D
    ! OPTIONAL: Array with data values
    DTYPE(D_TYPE), dimension(:,:), intent(inout), optional :: Data
#endif
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: ipos,jpos

    ipos = rlist%p_Knext(LHEAD)
    jpos = 0

#ifdef D

    if (present(Data)) then

      kd_copy: do
        jpos = jpos+1
        Key(jpos)    = rlist%p_Key(ipos)
        Data(:,jpos) = rlist%p_Data(:,ipos)
        if (ipos .eq. rlist%p_Knext(LTAIL)) exit kd_copy
        ipos = rlist%p_Knext(ipos)
      end do kd_copy

    else

      k_copy: do
        jpos = jpos+1
        Key(jpos) = rlist%p_Key(ipos)
        if (ipos .eq. rlist%p_Knext(LTAIL)) exit k_copy
        ipos = rlist%p_Knext(ipos)
      end do k_copy

    end if

#else

    copy: do
      jpos = jpos+1
      Key(jpos) = rlist%p_Key(ipos)
      if (ipos .eq. rlist%p_Knext(LTAIL)) exit copy
      ipos = rlist%p_Knext(ipos)
    end do copy

#endif

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine TEMPLATE_D(list_swap,T,D)(rlist1, rlist2)

!<description>
    ! This subroutine swaps content of two lists
!</description>
    
!<inputoutput>
    ! First linked list
    type(TEMPLATE_D(t_list,T,D)), intent(inout) :: rlist1

    ! Second linked list
    type(TEMPLATE_D(t_list,T,D)), intent(inout) :: rlist2
!</inputoutput>
!</subroutine>

    ! local variables
    type(TEMPLATE_D(t_list,T,D)) :: rlist
    
    ! Swap lists
    rlist  = rlist1
    rlist1 = rlist2
    rlist2 = rlist

    ! Reassociated pointers of first list
    call storage_getbase_int(rlist1%h_Knext, rlist1%p_Knext)
    if (rlist1%h_Kprev .ne. ST_NOHANDLE) &
        call storage_getbase_int(rlist1%h_Kprev, rlist1%p_Kprev)

    ! Reassociated pointers of second list
    call storage_getbase_int(rlist2%h_Knext, rlist2%p_Knext)
    if (rlist2%h_Kprev .ne. ST_NOHANDLE) &
        call storage_getbase_int(rlist2%h_Kprev, rlist2%p_Kprev)

#if defined D && defined D_STORAGE
    if (rlist1%h_Data .ne. ST_NOHANDLE) &
        call storage_getbase(rlist1%h_Data, rlist1%p_Data)
    if (rlist2%h_Data .ne. ST_NOHANDLE) &
        call storage_getbase(rlist2%h_Data, rlist2%p_Data)
#endif
    
  end subroutine

  !************************************************************************

!<function>
  
  pure function TEMPLATE_D(list_first,T,D)(rlist) result(ipos)

!<description>
    ! This function returns the position of the first list item
!</description>

!<input>
    ! The linked list
    type(TEMPLATE_D(t_list,T,D)), intent(in) :: rlist
!</input>

!<result>
    ! Position of first item
    integer :: ipos
!</result>
!</function>
    
    ipos = rlist%p_Knext(LHEAD)

  end function

  !************************************************************************

!<function>
  
  pure function TEMPLATE_D(list_last,T,D)(rlist) result(ipos)

!<description>
    ! This function returns the position of the last list item
!</description>

!<input>
    ! The linked list
    type(TEMPLATE_D(t_list,T,D)), intent(in) :: rlist
!</input>

!<result>
    ! Position of last item
    integer :: ipos
!</result>
!</function>
    
    ipos = rlist%p_Knext(LTAIL)

  end function

  !************************************************************************

!<function>
  
  function TEMPLATE_D(list_next,T,D)(rlist,breset) result(ipos)

!<description>
    ! This function returns the position of the next list item
    ! and resets the list if required. If the last list item is
    ! reached, then LNULL is returned and the list is reset
    ! implicitely.
!</description>

!<input>
    ! Reset list?
    logical, intent(in) :: breset
!</input>

!<inputoutput>
    ! The linked list
    type(TEMPLATE_D(t_list,T,D)), intent(inout) :: rlist
!</inputoutput>

!<result>
    ! Position of next item
    integer :: ipos
!</result>
!</function>
    
    ! Reset?
    if (breset) rlist%item=LHEAD

    ipos       = rlist%p_Knext(rlist%item)
    rlist%item = ipos

  end function

  !************************************************************************

!<function>
  
  function TEMPLATE_D(list_prev,T,D)(rlist,breset) result(ipos)

!<description>
    ! This function returns the position of the previous list item
    ! and resets the list if required. If the last list item is
    ! reached, then LNULL is returned and the list is reset
    ! implicitely.
!</description>

!<input>
    ! Reset list?
    logical, intent(in) :: breset
!</input>

!<inputoutput>
    ! The linked list
    type(TEMPLATE_D(t_list,T,D)), intent(inout) :: rlist
!</inputoutput>

!<result>
    ! Position of previous  item
    integer :: ipos
!</result>
!</function>
    
    if (rlist%clinkType .ne. LIST_DOUBLELINKED) then
      call output_line('This operation is only available for doubly-linked lists!',&
          OU_CLASS_ERROR,OU_MODE_STD,'list_prev')
      call sys_halt()
    end if

    ! Reset?
    if (breset) rlist%item=LTAIL
    
    if (rlist%item .eq. LNULL) then
      ipos = rlist%p_Knext(LTAIL)
    else
      ipos = rlist%p_Kprev(rlist%item)
    end if
    rlist%item = ipos

  end function

  !************************************************************************

!<subroutine>

#ifdef D
  subroutine TEMPLATE_D(list_get,T,D)(rlist, ipos, key, Data)
#else
  subroutine TEMPLATE_D(list_get,T,D)(rlist, ipos, key)
#endif

!<description>
    ! This subroutine returns the key and data stored at the given position
!</description>

!<input>
    ! The linked list
    type(TEMPLATE_D(t_list,T,D)), intent(inout) :: rlist

    ! Position of the data
    integer, intent(in) :: ipos
!</input>

!<output>
    ! key value
    TTYPE(T_TYPE), intent(out) :: key

#ifdef D
    ! OPTIONAL: Data
    DTYPE(D_TYPE), dimension(:), intent(out), optional :: Data
#endif
!</output>
!</subroutine>

    ! Check if position is valid
    if (ipos > rlist%na) then
      call output_line('Invalid position in list!',&
          OU_CLASS_ERROR,OU_MODE_STD,'list_get')
      call sys_halt()
    end if

    ! Get key
    key = rlist%p_key(ipos)

#ifdef D
    if (present(Data) .and. rlist%isizeData > 0)&
        Data = rlist%p_Data(:,ipos)
#endif

  end subroutine

  !************************************************************************

!<subroutine>

#ifdef D
  subroutine TEMPLATE_D(list_prepend,T,D)(rlist, key, ipos, Data)
#else
  subroutine TEMPLATE_D(list_prepend,T,D)(rlist, key, ipos)
#endif

!<description>
    ! This subroutine prepends a new item to the list
!</description>

!<input>
    ! key value
    TTYPE(T_TYPE), intent(in) :: key

#ifdef D
    ! OPTIONAL: Data
    DTYPE(D_TYPE), dimension(:), intent(in), optional :: Data
#endif
!</input>

!<inputoutput>
    ! The linked list
    type(TEMPLATE_D(t_list,T,D)), intent(inout) :: rlist
!</inputoutput>

!<output>
    ! Position of the prepended item
    integer, intent(out) :: ipos
!</output>
!</subroutine>

    ! Check if list needs to be resized
    rlist%na = rlist%na+1
    ipos     = rlist%p_Knext(LFREE)
    if (abs(ipos) > rlist%nna) then
      call list_resize(rlist, ceiling(rlist%dfactor*rlist%nna))
    end if

    ! Set next free position
    if (ipos > 0) then
      rlist%p_Knext(LFREE) = ipos+1
    else
      ipos = abs(ipos)
      rlist%p_Knext(LFREE) = rlist%p_Knext(ipos)
    end if

    ! Set head, tail and data
    select case (rlist%clinkType)
    case (LIST_SINGLELINKED)
      if (rlist%p_Knext(LHEAD) .eq. LNULL) then
        rlist%p_Knext(LHEAD) = ipos
        rlist%p_Knext(LTAIL) = ipos
        rlist%p_Knext(ipos)  = LNULL
      else
        rlist%p_Knext(ipos)  = rlist%p_Knext(LHEAD)
        rlist%p_Knext(LHEAD) = ipos
      end if
      
    case (LIST_DOUBLELINKED)
      if (rlist%p_Knext(LHEAD) .eq. LNULL) then
        rlist%p_Knext(LHEAD) = ipos
        rlist%p_Knext(LTAIL) = ipos
        rlist%p_Knext(ipos)  = LNULL
        rlist%p_Kprev(ipos)  = LNULL
      else
        rlist%p_Kprev(rlist%p_Knext(LHEAD)) = ipos
        rlist%p_Knext(ipos)  = rlist%p_Knext(LHEAD)
        rlist%p_Knext(LHEAD) = ipos
        rlist%p_Kprev(ipos)  = LNULL
      end if
      
    case DEFAULT
      call output_line('Invalid link type!',&
          OU_CLASS_ERROR,OU_MODE_STD,'list_prepend')
      call sys_halt()
    end select

    ! Store key value
    rlist%p_Key(ipos) = key

#ifdef D
    ! Store data value
    if (present(Data)) then
      rlist%p_Data(:,ipos) = Data
    end if
#endif

  end subroutine

  !************************************************************************

!<subroutine>

#ifdef D
  subroutine TEMPLATE_D(list_append,T,D)(rlist, key, ipos, Data)
#else
  subroutine TEMPLATE_D(list_append,T,D)(rlist, key, ipos)
#endif

!<description>
    ! This subroutine appends a new item to the list
!</description>

!<input>
    ! key value
    TTYPE(T_TYPE), intent(in) :: key

#ifdef D
    ! OPTIONAL: Data
    DTYPE(D_TYPE), dimension(:), intent(in), optional :: Data
#endif
!</input>

!<inputoutput>
    ! The linked list
    type(TEMPLATE_D(t_list,T,D)), intent(inout) :: rlist
!</inputoutput>

!<output>
    ! Position of the prepended item
    integer, intent(out) :: ipos
!</output>
!</subroutine>

    ! Check if list needs to be resized
    rlist%na = rlist%na+1
    ipos     = rlist%p_Knext(LFREE)
    if (abs(ipos) > rlist%nna) then
      call list_resize(rlist, ceiling(rlist%dfactor*rlist%nna))
    end if

    ! Set next free position
    if (ipos > 0) then
      rlist%p_Knext(LFREE) = ipos+1
    else
      ipos = abs(ipos)
      rlist%p_Knext(LFREE) = rlist%p_Knext(ipos)
    end if
    
    ! Set head, tail and data
    select case(rlist%clinkType)
    case (LIST_SINGLELINKED)
      if (rlist%p_Knext(LHEAD) .eq. LNULL) then
        rlist%p_Knext(LHEAD) = ipos
        rlist%p_Knext(LTAIL) = ipos
        rlist%p_Knext(ipos)  = LNULL
      else
        rlist%p_Knext(rlist%p_Knext(LTAIL)) = ipos
        rlist%p_Knext(LTAIL) = ipos
        rlist%p_Knext(ipos)  = LNULL
      end if
      
    case (LIST_DOUBLELINKED)
      if (rlist%p_Knext(LHEAD) .eq. LNULL) then
        rlist%p_Knext(LHEAD) = ipos
        rlist%p_Knext(LTAIL) = ipos
        rlist%p_Knext(ipos)  = LNULL
        rlist%p_Kprev(ipos)  = LNULL
      else
        rlist%p_Kprev(ipos)  = rlist%p_Knext(LTAIL)
        rlist%p_Knext(rlist%p_Knext(LTAIL)) = ipos
        rlist%p_Knext(LTAIL) = ipos
        rlist%p_Knext(ipos)  = LNULL
      end if
      
    case DEFAULT
      call output_line('Invalid link type!',&
          OU_CLASS_ERROR,OU_MODE_STD,'list_append')
      call sys_halt()
    end select

    ! Store key value
    rlist%p_Key(ipos) = key

#ifdef D
    ! Store data value
    if (present(Data)) then
      rlist%p_Data(:,ipos) = Data
    end if
#endif

  end subroutine

  !************************************************************************

!<subroutine>

#ifdef D
  subroutine TEMPLATE_D(list_insert,T,D)(rlist, key, ipred, ipos, Data)
#else
  subroutine TEMPLATE_D(list_insert,T,D)(rlist, key, ipred, ipos)
#endif

!<description>
    ! This subroutine inserts a new item into the list AFTER position ipred
!</description>

!<input>
    ! key value
    TTYPE(T_TYPE), intent(in) :: key

    ! position of predecessor
    integer, intent(in) :: ipred

#ifdef D
    ! OPTIONAL: Data
    DTYPE(D_TYPE), dimension(:), intent(in), optional :: Data
#endif
!</input>

!<inputoutput>
    ! The linked list
    type(TEMPLATE_D(t_list,T,D)), intent(inout) :: rlist
!</inputoutput>

!<output>
    ! Position of the prepended item
    integer, intent(out) :: ipos
!</output>
!</subroutine>

    ! Check if list needs to be resized
    rlist%na = rlist%na+1
    ipos     = rlist%p_Knext(LFREE)
    if (abs(ipos) > rlist%nna) then
      call list_resize(rlist, ceiling(rlist%dfactor*rlist%nna))
    end if

    ! Set next free position
    if (ipos > 0) then
      rlist%p_Knext(LFREE) = ipos+1
    else
      ipos = abs(ipos)
      rlist%p_Knext(LFREE) = rlist%p_Knext(ipos)
    end if
    
    ! Set head, tail and data
    select case(rlist%clinkType)
    case (LIST_SINGLELINKED)
      if (rlist%p_Knext(LHEAD) .eq. LNULL) then
        rlist%p_Knext(LHEAD) = ipos
        rlist%p_Knext(LTAIL) = ipos
        rlist%p_Knext(ipos)  = LNULL
      elseif (ipred .eq. rlist%p_Knext(LTAIL)) then
        rlist%p_Knext(ipred) = ipos
        rlist%p_Knext(LTAIL) = ipos
        rlist%p_Knext(ipos)  = LNULL
      else
        rlist%p_Knext(ipos)  = rlist%p_Knext(ipred)
        rlist%p_Knext(ipred) = ipos
      end if

    case (LIST_DOUBLELINKED)
      if (rlist%p_Knext(LHEAD) .eq. LNULL) then
        rlist%p_Knext(LHEAD) = ipos
        rlist%p_Knext(LTAIL) = ipos
        rlist%p_Knext(ipos)  = LNULL
        rlist%p_Kprev(ipos)  = LNULL
      elseif (ipred .eq. rlist%p_Knext(LTAIL)) then
        rlist%p_Kprev(ipos)  = rlist%p_Knext(LTAIL)
        rlist%p_Knext(ipred) = ipos
        rlist%p_Knext(LTAIL) = ipos
        rlist%p_Knext(ipos)  = LNULL
      else
        rlist%p_Kprev(rlist%p_Knext(ipred)) = ipos
        rlist%p_Knext(ipos)  = rlist%p_Knext(ipred)
        rlist%p_Knext(ipred) = ipos
        rlist%p_Kprev(ipos) = ipred
      end if

    case DEFAULT
      call output_line('Invalid link type!',&
          OU_CLASS_ERROR,OU_MODE_STD,'list_insert')
      call sys_halt()
    end select

    ! Store key value
    rlist%p_Key(ipos) = key

#ifdef D
    ! Store data value
    if (present(Data)) then
      rlist%p_Data(:,ipos) = Data
    end if
#endif

  end subroutine

  !************************************************************************

!<function>

  function TEMPLATE_D(list_delete,T,D)(rlist, key) result(f)

!<description>
    ! This function deletes an item from the list
!</description>

!<input>
    ! key value
    TTYPE(T_TYPE), intent(in) :: key
!</input>

!<inputoutput>
    ! The linked list
    type(TEMPLATE_D(t_list,T,D)), intent(inout) :: rlist
!</inputoutput>

!<result>
    ! Result of the deletion LIST_NOUT_FOUND / LIST_FOUND
    integer :: f
!</result>
!</function>

    ! local variables
    integer :: ipred,ipos

    ! Search for data
    f = list_search(rlist, key, ipred)
    if (f .eq. LIST_NOT_FOUND) return
    
    ! Delete item
    rlist%na = rlist%na-1
    ipos     = rlist%p_Knext(ipred)
    if (rlist%p_Knext(ipred) .eq. rlist%p_Knext(LTAIL))&
        rlist%p_Knext(LTAIL) = ipred

    ! Update free position
    rlist%p_Knext(ipred) = rlist%p_Knext(ipos)
    rlist%p_Knext(ipos)  = rlist%p_Knext(LFREE)
    rlist%p_Knext(LFREE) = -ipos

  end function

  !************************************************************************

!<function>

  function TEMPLATE_D(list_search,T,D)(rlist, key, ipred) result(f)

!<description>
    ! This function searches for a given key value in the list
!</description>

!<input>
    ! key value
    TTYPE(T_TYPE), intent(in) :: key

    ! The linked list
    type(TEMPLATE_D(t_list,T,D)), intent(in) :: rlist
!</input>

!<output>
    ! Position of the predecessor of the item found
    integer, intent(out) :: ipred
!</output>

!<result>
    ! Result of the search LIST_NOUT_FOUND / LIST_FOUND
    integer :: f
!</result>
!</function>

    ! local variables
    integer :: inext
    
    ! Initialization
    f     = LIST_NOT_FOUND
    ipred = LHEAD
    inext = rlist%p_Knext(ipred)

    ! Check if list is empty
    if (inext .eq. LNULL) return

    ! What kind of ordering are we
    select case(rlist%cordering)
    case (LIST_UNORDERED)
      do while(ipred.ne.rlist%p_Knext(LTAIL))
        inext = rlist%p_Knext(ipred)
        if (rlist%p_Key(inext) .eq. key) then
          f = LIST_FOUND; exit
        end if
        ipred = rlist%p_Knext(ipred)
      end do
      
    case (LIST_INCREASING)
      do while(ipred.ne.rlist%p_Knext(LTAIL))
        inext = rlist%p_Knext(ipred)
        if (rlist%p_Key(inext) .eq. key) then
          f = LIST_FOUND; exit
        end if
        
        if (rlist%p_Key(inext) > key) exit
        ipred = rlist%p_Knext(ipred)
      end do
      
    case (LIST_DECREASING)
      do while(ipred.ne.rlist%p_Knext(LTAIL))
        inext = rlist%p_Knext(ipred)
        if (rlist%p_Key(inext) .eq. key) then
          f = LIST_FOUND; exit
        end if
        
        if (rlist%p_Key(inext) < key) exit
        ipred = rlist%p_Knext(ipred)
      end do
      
    case (LIST_CSR7)
      do while(ipred.ne.rlist%p_Knext(LTAIL))
        inext = rlist%p_Knext(ipred)
        if (rlist%p_Key(inext) .eq. key) then
          f = LIST_FOUND; exit
        end if
        
        if ((ipred.ne.LHEAD) .and. rlist%p_Key(inext) > key) exit
        
        if (rlist%p_Knext(ipred) .eq. rlist%p_Knext(LTAIL)) then
          ipred = rlist%p_Knext(ipred); exit
        end if
        ipred = rlist%p_Knext(ipred)
      end do
    end select

  end function

  !************************************************************************

  subroutine TEMPLATE_D(list_print,T,D)(rlist)

!<description>
    ! This subroutine prints the content of the list
!</description>

!<input>
    ! The linked list
    type(TEMPLATE_D(t_list,T,D)), intent(in) :: rlist
!</input>
!</subroutine>

    ! local variable
    integer :: ipos
    
#ifdef T_STORAGE

    ipos = rlist%p_Knext(LHEAD)
    if (ipos .eq. LNULL) return
    
    do while (ipos .ne. rlist%p_Knext(LTAIL))
      write(*,*) ipos,":",rlist%p_Key(ipos)
      ipos = rlist%p_Knext(ipos)
    end do
    write(*,*) ipos,":",rlist%p_Key(ipos)

#else
    
    call output_line('Unable to print list with derived data type!',&
        OU_CLASS_WARNING,OU_MODE_STD,'list_print')

#endif

  end subroutine

end module

#endif
